/**
 * Luan Machado Bernardt | GRR20190363
 * Lucas Müller          | GRR20197160
 */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

#include "matrixLib.h"


int main (int argc, char **argv) {

    t_sist *SL;
    double tempoTri, tempoUx, tempoLy;
    char opt;
    FILE *f_out = stdout;
    _Bool pivotP=0;

    while (-1 != (opt = getopt(argc,argv,"po:"))) {
        switch(opt) {
        case 'p':
            pivotP = 1;
            break;
        case 'o':
            f_out = fopen(optarg, "wb");
            break;
        default:
            fprintf(stderr,
              "Uso: %s [-p|-o <arquivo-saida>]\n"
              "\t-p com pivoteamento parcial\n"
              "\t-o imprimir resultados na saida especificada\n", 
              argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    float *vet;
    while (!feof(stdin))
    {
      SL = SL_leitura();
      if (!SL) return EXIT_FAILURE;

      vet = SL_interpolacao(SL, 0);
      if (!vet) return EXIT_FAILURE;

      SL_printMatrix(f_out, vet, SL->n, 1);
#if 0
      SL_printMatrix(f_out, SL->Int, SL->n, SL->n);
      SL_printMatrix(f_out, SL->L, SL->n, SL->n);
      SL_printMatrix(f_out, SL->U, SL->n, SL->n);
      SL_printMatrix(f_out, SL->A, SL->n, SL->m);
#endif
      SL_libera(SL);
      free(vet);
    }

    if (f_out != stdout) fclose(f_out);

    return EXIT_SUCCESS;
}
