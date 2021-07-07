#define _GNU_SOURCE /* asprintf() */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h> /* PATH_MAX */
#include <dirent.h>

#define DIR_PATH "./Resultados" // diretório dos arquivos gerados
#define N_HEADER_LN 6 // qtd de linhas reservadas ao cabeçalho

const char L3_FILE_PREFIX[] = "L3";
const char L2_FILE_PREFIX[] = "L2CACHE";
const char DP_FILE_PREFIX[] = "FLOPS_DP";

const char L3_FIELD[]  = "L3 bandwidth [MBytes/s]";
const char L2_FIELD[]  = "L2 miss ratio";
const char DP_FIELD1[] = "DP MFLOP/s";
const char DP_FIELD2[] = "AVX DP MFLOP/s";

// valores extraídos a partir dos campos
struct val {       
  double d[25]; // armazena valores
  size_t n;     // qtd de valores
};

// para qsort()
int cstrcmp(const void *a, const void *b) {
  return strcmp(*(const char**)a, *(const char**)b);
}

// imprime elementos de struct val lado a lado
void printTable(struct val *data) 
{
  for (int i=0; i < data->n-1; ++i)
    printf("%-10g ", data->d[i]);
  printf("%-10g", data->d[data->n-1]);
}

// recolhe o nome dos arquivos encontrados em folder, com extensão 'ext'
char **getFiles(char path[], char ext_filter[], size_t *n_files) 
{
  char **fnames=NULL;

  // obtem nome dos arquivos de extensão .txt do diretório
  struct dirent *dir;
  DIR *d = opendir(path);
  if (NULL == d) return NULL;

  char *ext;
  size_t i=0;
  while ((dir = readdir(d))) {
    ext = strchr(dir->d_name, '.');
    if (ext && 0 == strcmp(ext, ext_filter)) 
    {
      ++i;
      char **tmp = realloc(fnames, sizeof(char*) * i);
      if (NULL == tmp) return NULL;

      fnames = tmp;

      asprintf(&fnames[i-1], "%s/%s", path, dir->d_name);
      if (NULL == fnames[i-1]) return NULL;
    }
  }
  closedir(d);
    
  // ordena arquivos alfabeticamente
  qsort(fnames, i, sizeof(char*), &cstrcmp);

  *n_files = i;
  return fnames;
}

int main(void)
{
  size_t n_files=0; // qtd de arquivos
  char **fnames = getFiles(DIR_PATH, ".txt", &n_files); // extrair arquivos .txt
  assert(fnames != NULL && "Não foi possível receber 'fnames'");

  struct val L3_64={0}, L3_100={0}, L3_128={0}, L3_2000={0}, L3_2048={0};
  struct val L2_64={0}, L2_100={0}, L2_128={0}, L2_2000={0}, L2_2048={0};
  struct val DP_64={0}, DP_100={0}, DP_128={0}, DP_2000={0}, DP_2048={0};
  struct val AVX_64={0}, AVX_100={0}, AVX_128={0}, AVX_2000={0}, AVX_2048={0};

  // extrai informações dos arquivos
  FILE *f;
  char ln[4096], *aux;
  struct val *data;
  for (size_t i=0; i < n_files; ++i) 
  {
    f = fopen(fnames[i], "rb");
    assert(f != NULL && "Não foi possível abrir arquivo");

    // pula linhas de cabeçalho
    int ln_counter=0, c;
    while (ln_counter < N_HEADER_LN && EOF != (c = fgetc(f))) {
      if ('\n' == c) ++ln_counter;
    }

    const char *field;
    int offset=0;
    if (NULL != (aux = strstr(fnames[i], L3_FILE_PREFIX))) {
      switch (strtol(aux+sizeof(L3_FILE_PREFIX), NULL, 10)) {
      case 64:   data = &L3_64;   break;
      case 100:  data = &L3_100;  break;
      case 128:  data = &L3_128;  break;
      case 2000: data = &L3_2000; break;
      case 2048: data = &L3_2048; break;
      default:   exit(EXIT_FAILURE);
      }
      field = L3_FIELD;
      offset = sizeof(L3_FIELD);
    }
    else if (NULL != (aux = strstr(fnames[i], L2_FILE_PREFIX))) {
      switch (strtol(aux+sizeof(L2_FILE_PREFIX), NULL, 10)) {
      case 64:   data = &L2_64;   break;
      case 100:  data = &L2_100;  break;
      case 128:  data = &L2_128;  break;
      case 2000: data = &L2_2000; break;
      case 2048: data = &L2_2048; break;
      default:   exit(EXIT_FAILURE);
      }
      field = L2_FIELD;
      offset = sizeof(L2_FIELD);
    }
    else if (NULL != (aux = strstr(fnames[i], DP_FILE_PREFIX))) {
      int size = strtol(aux + sizeof(DP_FILE_PREFIX), NULL, 10);
      while (fgets(ln, sizeof(ln), f)) {
        if (NULL != (aux = strstr(ln, DP_FIELD2))) {
          switch (size) {
          case 64:   data = &AVX_64;   break;
          case 100:  data = &AVX_100;  break;
          case 128:  data = &AVX_128;  break;
          case 2000: data = &AVX_2000; break;
          case 2048: data = &AVX_2048; break;
          default:   exit(EXIT_FAILURE);
          }
          data->d[data->n] = strtod(aux + sizeof(DP_FIELD2), NULL);
          ++data->n;
        }
        else if (NULL != (aux = strstr(ln, DP_FIELD1))) {
          switch (size) {
          case 64:   data = &DP_64;   break;
          case 100:  data = &DP_100;  break;
          case 128:  data = &DP_128;  break;
          case 2000: data = &DP_2000; break;
          case 2048: data = &DP_2048; break;
          default:   exit(EXIT_FAILURE);
          }
          data->d[data->n] = strtod(aux + sizeof(DP_FIELD1), NULL);
          ++data->n;
        }
      }
      field = NULL;
      offset = 0;
    }
#if 0
    else {
      fprintf(stderr, "Arquivo inesperado: %s\n", fnames[i]);
      exit(EXIT_FAILURE);
    }
#endif
    if (!field) {
      fclose(f);
      continue;
    }

    while (fgets(ln, sizeof(ln), f)) {
      if (NULL != (aux = strstr(ln, field))) {
        data->d[data->n] = strtod(aux + offset, NULL);
        ++data->n;
      }
    }
  }

  puts("#N\t\t\tL3");
  printf("64\t");   printTable(&L3_64);   putchar('\n');
  printf("100\t");  printTable(&L3_100);  putchar('\n');
  printf("128\t");  printTable(&L3_128);  putchar('\n');
  printf("2000\t"); printTable(&L3_2000); putchar('\n');
  printf("2048\t"); printTable(&L3_2048); putchar('\n');
  putchar('\n');

  puts("#N\t\t\tL2CACHE");
  printf("64\t");   printTable(&L2_64);   putchar('\n');
  printf("100\t");  printTable(&L2_100);  putchar('\n');
  printf("128\t");  printTable(&L2_128);  putchar('\n');
  printf("2000\t"); printTable(&L2_2000); putchar('\n');
  printf("2048\t"); printTable(&L2_2048); putchar('\n');
  putchar('\n');

  puts("#N\t\t\tFLOPS_DP\t\t\t\t\tFLOPS_AVX");
  printf("64\t");   printTable(&L2_64);   putchar('\t'); printTable(&AVX_64);   putchar('\n');
  printf("100\t");  printTable(&L2_100);  putchar('\t'); printTable(&AVX_100);  putchar('\n');
  printf("128\t");  printTable(&L2_128);  putchar('\t'); printTable(&AVX_128);  putchar('\n');
  printf("2000\t"); printTable(&L2_2000); putchar('\t'); printTable(&AVX_2000); putchar('\n');
  printf("2048\t"); printTable(&L2_2048); putchar('\t'); printTable(&AVX_2048); putchar('\n');
  putchar('\n');

  return EXIT_SUCCESS;
}
