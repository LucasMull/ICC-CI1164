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

    while (!feof(stdin))
    {
      SL = SL_leitura();
      if (SL) {
        fprintf(f_out,"\nN = %d",SL->n);
        fprintf(f_out,"\n");
        
        fprintf(f_out,"%.20s\n", "Original ##################");
        SL_printMatrix(f_out, SL->A, SL->n);
        
        SL->Id = SL_geraIdentidade(SL->n);
        if (SL_triangulariza(SL, pivotP, &tempoTri) != -1) {

            SL_geraInversa(SL,&tempoLy,&tempoUx);
            fprintf(f_out,"%.20s\n", "Inversa ##################");
            SL_printMatrix(f_out,SL->Inv,SL->n);

            fprintf(f_out,"###########\n");
            fprintf(f_out,"# Tempo Triangularizacao: %e ms\n",tempoTri);
            fprintf(f_out,"# Tempo calculo de Y: %e ms\n",tempoLy);
            fprintf(f_out,"# Tempo calculo de X: %e ms\n",tempoUx);
            for (unsigned int i=0; i<SL->n; ++i) {
                fprintf(f_out, "# Norma L2 dos residuos (%d): ", i);
                fprintf(f_out, "%g\n",SL_normaL2Residuo(SL, &SL->Id[SL->n*i], i));
            }
            fprintf(f_out, "###########\n");
        }
        SL_libera(SL);
      }
    }

    if (f_out != stdout) fclose(f_out);

    return EXIT_SUCCESS;
}
