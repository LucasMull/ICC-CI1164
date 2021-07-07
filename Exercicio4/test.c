#define _GNU_SOURCE /* asprintf() */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h> /* PATH_MAX */
#include <dirent.h>

#define DIR_PATH "./Resultados" // diretório dos arquivos gerados
#define N_HEADER_LN 6 // qtd de linhas reservadas ao cabeçalho

#define L3_FIELD   "L3 bandwidth [MBytes/s]"
#define L2_FIELD   "L2 miss ratio"
#define AOP_FIELD1 "DP MFLOP/s"
#define AOP_FIELD2 "AVX DP MFLOP/s"

// valores extraídos a partir dos campos
struct val {       
  char str[25][25]; // armazena valores como string
  size_t n;         // qtd de valores (max 25)
};

// para qsort()
static int cstrcmp(const void *a, const void *b) {
  return strcmp(*(const char**)a, *(const char**)b);
}

int main(void)
{
  char **fnames=NULL; // nome dos arquivos para extrair
  size_t n_files=0; // qtd de arquivos
  struct val L3={0}, L2={0}, AOP={0}; // armazena dados extraidos

  // obtem nome dos arquivos de extensão .txt do diretório
  struct dirent *dir;
  DIR *d = opendir(DIR_PATH);
  assert(d != NULL && "Arquivo inexistente");

  char *ext;
  while ((dir = readdir(d))) {
    ext = strchr(dir->d_name, '.');
    if (ext && 0 == strcmp(ext, ".txt")) 
    {
      ++n_files;
      char **tmp = realloc(fnames, sizeof(char*) * n_files);
      assert(tmp != NULL && "Sem memória");

      fnames = tmp;

      asprintf(&fnames[n_files-1], DIR_PATH"/%s", dir->d_name);
      assert(fnames[n_files-1] != NULL && "Sem memória");
    }
  }
  closedir(d);
    
  // ordena arquivos alfabeticamente
  qsort(fnames, n_files, sizeof(char*), &cstrcmp);

  // extrai informações dos arquivos
  FILE *f;
  char ln[4096], *aux;
  for (size_t i=0; i < n_files; ++i) {
    f = fopen(fnames[i], "rb");
    assert(f != NULL && "Não foi possível abrir arquivo");

    // pula linhas de cabeçalho
    int ln_counter=0, c;
    while (ln_counter < N_HEADER_LN && EOF != (c = fgetc(f))) {
      if ('\n' == c) 
        ++ln_counter;
    }

    if (strstr(fnames[i], "L3")) {
      while (fgets(ln, sizeof(ln), f)) {
        if (NULL != (aux = strstr(ln, L3_FIELD))) {
          aux += sizeof(L3_FIELD);
          strncpy(L3.str[L3.n], aux, sizeof(L3.str[0]));
          ++L3.n;
        }
      }
    }
    else if (strstr(fnames[i], "L2CACHE")) {
      while (fgets(ln, sizeof(ln), f)) {
        if (NULL != (aux = strstr(ln, L2_FIELD))) {
          aux += sizeof(L2_FIELD);
          strncpy(L2.str[L2.n], aux, sizeof(L2.str[0]));
          ++L2.n;
        }
      }
    }
    else if (strstr(fnames[i], "FLOPS_DP")) {
      while (fgets(ln, sizeof(ln), f)) {
        if (NULL != (aux = strstr(ln, AOP_FIELD1))) {
          aux += sizeof(AOP_FIELD1);
          strncpy(AOP.str[AOP.n], aux, sizeof(AOP.str[0]));
          ++AOP.n;
        }
        else if (NULL != (aux = strstr(ln, AOP_FIELD2))) {
          aux += sizeof(AOP_FIELD2);
          strncpy(AOP.str[AOP.n], aux, sizeof(AOP.str[0]));
          ++AOP.n;
        }
      }
    }
#if 0
    else {
      fprintf(stderr, "Arquivo inesperado: %s\n", fnames[i]);
      exit(EXIT_FAILURE);
    }
#endif
    fclose(f);
  }

  return 0;
}
