#define _GNU_SOURCE /* asprintf() */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h> /* PATH_MAX */
#include <dirent.h>

#define DIR_PATH "./Resultados" // diretório dos arquivos gerados
#define N_HEADER_LN 6 // qtd de linhas reservadas ao cabeçalho

#define RUNTIME_FIELD "Runtime (RDTSC) [s]"

const char L3_FILE_PREFIX[] = "L3";
const char L2_FILE_PREFIX[] = "L2CACHE";
const char DP_FILE_PREFIX[] = "FLOPS_DP";

const char L3_FIELD[]  = "L3 bandwidth [MBytes/s]";
const char L2_FIELD[]  = "L2 miss ratio";
const char DP_FIELD1[] = "DP MFLOP/s";
const char DP_FIELD2[] = "AVX DP MFLOP/s";

// valores extraídos a partir dos campos
struct val {       
  double d[16];          // armazena valores
  size_t n;              // qtd de valores
  double runtime_normal; // tempo de execução
  double runtime_otimiz; // tempo de execução
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

  struct val L3_64={0}, L3_100={0}, L3_128={0}, L3_1000={0}, L3_1024={0}, L3_2000={0}, L3_2048={0}, L3_3000={0}, L3_4096={0}, L3_5000={0};
  struct val L2_64={0}, L2_100={0}, L2_128={0}, L2_1000={0}, L2_1024={0},  L2_2000={0}, L2_2048={0}, L2_3000={0}, L2_4096={0}, L2_5000={0};
  struct val DP_64={0}, DP_100={0}, DP_128={0}, DP_1000={0}, DP_1024={0}, DP_2000={0}, DP_2048={0}, DP_3000={0}, DP_4096={0}, DP_5000={0};
  struct val AVX_64={0}, AVX_100={0}, AVX_128={0}, AVX_1000={0}, AVX_1024={0}, AVX_2000={0}, AVX_2048={0}, AVX_3000={0}, AVX_4096={0}, AVX_5000={0};

  enum { OPTIMIZED, NORMAL } tag;

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
    long size;
    if (NULL != (aux = strstr(fnames[i], L3_FILE_PREFIX))) {
      switch (size = strtol(aux+sizeof(L3_FILE_PREFIX), NULL, 10)) {
      case 64:   data = &L3_64;   break;
      case 100:  data = &L3_100;  break;
      case 128:  data = &L3_128;  break;
      case 1000: data = &L3_1000; break;
      case 1024: data = &L3_1024; break;
      case 2000: data = &L3_2000; break;
      case 2048: data = &L3_2048; break;
      case 3000: data = &L3_3000; break;
      case 4096: data = &L3_4096; break;
      case 5000: data = &L3_5000; break;
      default:   fprintf(stderr, "L3 - %ld", size); exit(EXIT_FAILURE);
      }
      field = L3_FIELD;
      offset = sizeof(L3_FIELD);
    }
    else if (NULL != (aux = strstr(fnames[i], L2_FILE_PREFIX))) {
      switch (size = strtol(aux+sizeof(L2_FILE_PREFIX), NULL, 10)) {
      case 64:   data = &L2_64;   break;
      case 100:  data = &L2_100;  break;
      case 128:  data = &L2_128;  break;
      case 1000: data = &L2_1000; break;
      case 1024: data = &L2_1024; break;
      case 2000: data = &L2_2000; break;
      case 2048: data = &L2_2048; break;
      case 3000: data = &L2_3000; break;
      case 4096: data = &L2_4096; break;
      case 5000: data = &L2_5000; break;
      default:   fprintf(stderr, "L2 - %ld", size); exit(EXIT_FAILURE);
      }
      field = L2_FIELD;
      offset = sizeof(L2_FIELD);
    }
    else if (NULL != (aux = strstr(fnames[i], DP_FILE_PREFIX))) {
      size = strtol(aux + sizeof(DP_FILE_PREFIX), NULL, 10);
      while (fgets(ln, sizeof(ln), f)) {
        if (NULL != (aux = strstr(ln, DP_FIELD2))) {
          switch (size) {
          case 64:   data = &AVX_64;   break;
          case 100:  data = &AVX_100;  break;
          case 128:  data = &AVX_128;  break;
          case 1000: data = &AVX_1000; break;
          case 1024: data = &AVX_1024; break;
          case 2000: data = &AVX_2000; break;
          case 2048: data = &AVX_2048; break;
          case 3000: data = &AVX_3000; break;
          case 4096: data = &AVX_4096; break;
          case 5000: data = &AVX_5000; break;
          default:   fprintf(stderr, "AVX - %ld", size); exit(EXIT_FAILURE);
          }
          data->d[data->n] = strtod(aux + sizeof(DP_FIELD2), NULL);
          ++data->n;
        }
        else if (NULL != (aux = strstr(ln, DP_FIELD1))) {
          switch (size) {
          case 64:   data = &DP_64;   break;
          case 100:  data = &DP_100;  break;
          case 128:  data = &DP_128;  break;
          case 1000: data = &DP_1000; break;
          case 1024: data = &DP_1024; break;
          case 2000: data = &DP_2000; break;
          case 2048: data = &DP_2048; break;
          case 3000: data = &DP_3000; break;
          case 4096: data = &DP_4096; break;
          case 5000: data = &DP_5000; break;
          default:   fprintf(stderr, "DP - %ld", size); exit(EXIT_FAILURE);
          }
          data->d[data->n] = strtod(aux + sizeof(DP_FIELD1), NULL);
          ++data->n;
        }
        else if (NULL != (aux = strstr(ln, RUNTIME_FIELD))) {
          switch (size) {
          case 64:   data = &DP_64;   break;
          case 100:  data = &DP_100;  break;
          case 128:  data = &DP_128;  break;
          case 1000: data = &DP_1000; break;
          case 1024: data = &DP_1024; break;
          case 2000: data = &DP_2000; break;
          case 2048: data = &DP_2048; break;
          case 3000: data = &DP_3000; break;
          case 4096: data = &DP_4096; break;
          case 5000: data = &DP_5000; break;
          default:   fprintf(stderr, "DP - %ld", size); exit(EXIT_FAILURE);
          }
          aux += sizeof(RUNTIME_FIELD);
          if (OPTIMIZED == tag)
            data->runtime_otimiz = strtod(aux, NULL);
          else
            data->runtime_normal = strtod(aux, NULL);
        }
        else if (NULL != (aux = strstr(ln, "MatRowVetOtimiz")) || NULL != (aux = strstr(ln , "MatMatRowOtimiz"))) {
          tag = OPTIMIZED;
        }
        else if (NULL != (aux = strstr(ln, "MatRowVet")) || NULL != (aux = strstr(ln , "MatMatRow"))) {
          tag = NORMAL;
        }
      }
      data->runtime_otimiz /= data->n;
      data->runtime_normal /= data->n;
      continue; /* EARLY CONTINUE */
    }

    while (fgets(ln, sizeof(ln), f)) {
      if (NULL != (aux = strstr(ln, RUNTIME_FIELD))) {
        aux += sizeof(RUNTIME_FIELD);
        if (OPTIMIZED == tag)
          data->runtime_otimiz = strtod(aux, NULL);
        else
          data->runtime_normal = strtod(aux, NULL);
      }
      else if (NULL != (aux = strstr(ln, field))) {
        data->d[data->n] = strtod(aux + offset, NULL);
        ++data->n;
      }
      else if (NULL != (aux = strstr(ln, "MatRowVetOtimiz")) || NULL != (aux = strstr(ln , "MatMatRowOtimiz"))) {
        tag = OPTIMIZED;
      }
      else if (NULL != (aux = strstr(ln, "MatRowVet")) || NULL != (aux = strstr(ln , "MatMatRow"))) {
        tag = NORMAL;
      }
    }
    data->runtime_otimiz /= data->n;
    data->runtime_normal /= data->n;
  }

  puts("#N\t\t\tRUNTIME OTIMIZADO\n\tL3\t   L2CACHE    FLOPS_DP   FLOPS_AVX");
  printf(
      "64\t%-10g %-10g %-10g %-10g\n"
      "100\t%-10g %-10g %-10g %-10g\n"
      "128\t%-10g %-10g %-10g %-10g\n"
      "1000\t%-10g %-10g %-10g %-10g\n"
      "1024\t%-10g %-10g %-10g %-10g\n"
      "2000\t%-10g %-10g %-10g %-10g\n"
      "2048\t%-10g %-10g %-10g %-10g\n"
      "3000\t%-10g %-10g %-10g %-10g\n"
      "4096\t%-10g %-10g %-10g %-10g\n"
      "5000\t%-10g %-10g %-10g %-10g\n\n",
      L3_64.runtime_otimiz, L2_64.runtime_otimiz, DP_64.runtime_otimiz, AVX_64.runtime_otimiz,
      L3_100.runtime_otimiz, L2_100.runtime_otimiz, DP_100.runtime_otimiz, AVX_100.runtime_otimiz,
      L3_128.runtime_otimiz, L2_128.runtime_otimiz, DP_128.runtime_otimiz, AVX_128.runtime_otimiz,
      L3_1000.runtime_otimiz, L2_1000.runtime_otimiz, DP_1000.runtime_otimiz, AVX_1000.runtime_otimiz,
      L3_1024.runtime_otimiz, L2_1024.runtime_otimiz, DP_1024.runtime_otimiz, AVX_1024.runtime_otimiz,
      L3_2000.runtime_otimiz, L2_2000.runtime_otimiz, DP_2000.runtime_otimiz, AVX_2000.runtime_otimiz,
      L3_2048.runtime_otimiz, L2_2048.runtime_otimiz, DP_2048.runtime_otimiz, AVX_2048.runtime_otimiz,
      L3_3000.runtime_otimiz, L2_3000.runtime_otimiz, DP_3000.runtime_otimiz, AVX_3000.runtime_otimiz,
      L3_4096.runtime_otimiz, L2_4096.runtime_otimiz, DP_4096.runtime_otimiz, AVX_4096.runtime_otimiz,
      L3_5000.runtime_otimiz, L2_5000.runtime_otimiz, DP_5000.runtime_otimiz, AVX_5000.runtime_otimiz);

  puts("#N\t\t\tRUNTIME NORMAL\n\tL3\t   L2CACHE    FLOPS_DP   FLOPS_AVX");
  printf(
      "64\t%-10g %-10g %-10g %-10g\n"
      "100\t%-10g %-10g %-10g %-10g\n"
      "128\t%-10g %-10g %-10g %-10g\n"
      "1000\t%-10g %-10g %-10g %-10g\n"
      "1024\t%-10g %-10g %-10g %-10g\n"
      "2000\t%-10g %-10g %-10g %-10g\n"
      "2048\t%-10g %-10g %-10g %-10g\n"
      "3000\t%-10g %-10g %-10g %-10g\n"
      "4096\t%-10g %-10g %-10g %-10g\n"
      "5000\t%-10g %-10g %-10g %-10g\n\n",
      L3_64.runtime_normal, L2_64.runtime_normal, DP_64.runtime_normal, AVX_64.runtime_normal,
      L3_100.runtime_normal, L2_100.runtime_normal, DP_100.runtime_normal, AVX_100.runtime_normal,
      L3_128.runtime_normal, L2_128.runtime_normal, DP_128.runtime_normal, AVX_128.runtime_normal,
      L3_1000.runtime_normal, L2_1000.runtime_normal, DP_1000.runtime_normal, AVX_1000.runtime_normal,
      L3_1024.runtime_normal, L2_1024.runtime_normal, DP_1024.runtime_normal, AVX_1024.runtime_normal,
      L3_2000.runtime_normal, L2_2000.runtime_normal, DP_2000.runtime_normal, AVX_2000.runtime_normal,
      L3_2048.runtime_normal, L2_2048.runtime_normal, DP_2048.runtime_normal, AVX_2048.runtime_normal,
      L3_3000.runtime_normal, L2_3000.runtime_normal, DP_3000.runtime_normal, AVX_3000.runtime_normal,
      L3_4096.runtime_normal, L2_4096.runtime_normal, DP_4096.runtime_normal, AVX_4096.runtime_normal,
      L3_5000.runtime_normal, L2_5000.runtime_normal, DP_5000.runtime_normal, AVX_5000.runtime_normal);

  puts("#N\t\t\tL2");
  puts("\t     MatRowVet             MatMatRow");
  puts("\t Com        Sem        Com        Sem");
  printf("64\t");   printTable(&L2_64);   putchar('\n');
  printf("100\t");  printTable(&L2_100);  putchar('\n');
  printf("128\t");  printTable(&L2_128);  putchar('\n');
  printf("1000\t"); printTable(&L2_1000); putchar('\n');
  printf("1024\t"); printTable(&L2_1024); putchar('\n');
  printf("2000\t"); printTable(&L2_2000); putchar('\n');
  printf("2048\t"); printTable(&L2_2048); putchar('\n');
  printf("3000\t"); printTable(&L2_3000); putchar('\n');
  printf("4096\t"); printTable(&L2_4096); putchar('\n');
  printf("5000\t"); printTable(&L2_5000); putchar('\n');
  putchar('\n');

  puts("#N\t\t\tL3");
  puts("\t     MatRowVet             MatMatRow");
  puts("\t Com        Sem        Com        Sem");
  printf("64\t");   printTable(&L3_64);   putchar('\n');
  printf("100\t");  printTable(&L3_100);  putchar('\n');
  printf("128\t");  printTable(&L3_128);  putchar('\n');
  printf("1000\t"); printTable(&L3_1000); putchar('\n');
  printf("1024\t"); printTable(&L3_1024); putchar('\n');
  printf("2000\t"); printTable(&L3_2000); putchar('\n');
  printf("2048\t"); printTable(&L3_2048); putchar('\n');
  printf("3000\t"); printTable(&L3_3000); putchar('\n');
  printf("4096\t"); printTable(&L3_4096); putchar('\n');
  printf("5000\t"); printTable(&L3_5000); putchar('\n');
  putchar('\n');

  printf("#N\t\t\tFLOPS_DP\t\t\t\t\tFLOPS_AVX\n");
  printf("\t     MatRowVet             MatMatRow");
  printf("\t\t   MatRowVet             MatMatRow\n");
  printf("\t Com        Sem        Com        Sem");
  printf("\t\tCom        Sem        Com        Sem\n");
  printf("64\t");   printTable(&DP_64);   putchar('\t'); printTable(&AVX_64);   putchar('\n');
  printf("100\t");  printTable(&DP_100);  putchar('\t'); printTable(&AVX_100);  putchar('\n');
  printf("128\t");  printTable(&DP_128);  putchar('\t'); printTable(&AVX_128);  putchar('\n');
  printf("1000\t");  printTable(&DP_1000);  putchar('\t'); printTable(&AVX_1000);  putchar('\n');
  printf("1024\t");  printTable(&DP_1024);  putchar('\t'); printTable(&AVX_1024);  putchar('\n');
  printf("2000\t"); printTable(&DP_2000); putchar('\t'); printTable(&AVX_2000); putchar('\n');
  printf("2048\t"); printTable(&DP_2048); putchar('\t'); printTable(&AVX_2048); putchar('\n');
  printf("3000\t"); printTable(&DP_3000); putchar('\t'); printTable(&AVX_3000); putchar('\n');
  printf("4096\t"); printTable(&DP_4096); putchar('\t'); printTable(&AVX_4096); putchar('\n');
  printf("5000\t"); printTable(&DP_5000); putchar('\t'); printTable(&AVX_5000); putchar('\n');
  putchar('\n');

  return EXIT_SUCCESS;
}
