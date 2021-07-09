/**
 * Luan Machado Bernardt | GRR20190363
 * Lucas Müller          | GRR20197160
 */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

#include "matrixLib.h"


int main (int argc, char **argv) {

    t_matrix *Mat;
    double tempoTri, tempoUx, tempoLy;
    char opt;
    FILE *f_out=stdout;
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

    while (!feof(stdin))
    {
      Mat = readMatrix();
      if (!Mat) exit(EXIT_FAILURE);
      
      fprintf(f_out,"%.20s\n", "Original ##################");
      printMatrix(f_out,Mat->A,Mat->n);
      
      Mat->Id = geraIdentidade(Mat->n);
      triangularizaMatrix(Mat,pivotP,&tempoTri);

      geraInversa(Mat,Mat->Id,&tempoLy,&tempoUx);
      fprintf(f_out,"%.20s\n", "Inversa ##################");
      printMatrix(f_out,Mat->Inv,Mat->n);

      fprintf(f_out,"###############\n");
      fprintf(f_out,"# Tempo Triangularização: %e ms\n",tempoTri);
      fprintf(f_out,"# Tempo cálculo de Y: %e ms\n",tempoLy);
      fprintf(f_out,"# Tempo cálculo de X: %e ms\n",tempoUx);
      for (unsigned int i=0; i<Mat->n; ++i) {
          fprintf(f_out,"# Norma L2 dos residuos (%d): ", i);
          fprintf(f_out,"%g\n",normaL2Residuo(Mat,Mat->Id[i],i));
      }
      fprintf(f_out,"###############\n");

      limpaStruct(Mat);
    }

    if (f_out != stdout) fclose(f_out);

    return EXIT_SUCCESS;
}
