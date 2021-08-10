/**
 * Luan Machado Bernardt | GRR20190363
 * Lucas Müller          | GRR20197160
 */

#include <stdlib.h>
#include <stdio.h>

#include "matrixLib.h"


int main (int argc, char **argv) {

    t_sist *SL;
    FILE *f_out = stdout;

    double *vet;
    while (!feof(stdin))
    {
        SL = SL_leitura();
        if (!SL) return EXIT_SUCCESS;
        
        double *B = SL_alocaMatrix(1, SL->n);
        if (!B) return EXIT_FAILURE;

        for (int i=0; i<SL->m; ++i) {
            vet = SL_interpolacao(SL, i, B);
            if (!vet) return EXIT_FAILURE;
            SL_printMatrix(f_out, vet, SL->n, 1);

            free(vet);

            vet = SL_ajusteDeCurvas(SL, i, B);
            if (!vet) return EXIT_FAILURE;
            SL_printMatrix(f_out, vet, SL->n, 1);

            free(vet);
        }
#if 0
        SL_printMatrix(f_out, SL->Int, SL->n, SL->n);
        SL_printMatrix(f_out, SL->L, SL->n, SL->n);
        SL_printMatrix(f_out, SL->U, SL->n, SL->n);
        SL_printMatrix(f_out, SL->A, SL->n, SL->m);
#endif
        free(B);
        SL_libera(SL);
    }

    if (f_out != stdout) fclose(f_out);

    return EXIT_SUCCESS;
}
