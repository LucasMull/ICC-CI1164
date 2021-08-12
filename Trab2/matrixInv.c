/**
 * Luan Machado Bernardt | GRR20190363
 * Lucas Müller          | GRR20197160
 */

#include <stdlib.h>
#include <stdio.h>

#ifndef _NO_LIKWID
#include <likwid.h>
#endif

#include "matrixLib.h"


int main (int argc, char **argv) {

    t_sist *SL;
    FILE *f_out = stdout;

#ifdef _NO_LIKWID
    while (!feof(stdin))
    {
        SL = SL_leitura();
        if (!SL) break;
        
        double *B = SL_alocaMatrix(1, SL->n);
        if (!B) return EXIT_FAILURE;

        double *pol = SL_alocaMatrix(1, SL->n);
        if (!pol) return EXIT_FAILURE;

        // guardar valores de x exp
        double *lookup = SL_alocaMatrix(SL->n, SL->n);
        if (!lookup) return EXIT_FAILURE;

        for (int i=0; i<SL->m; ++i) {
            if (SL_interpolacao(SL, i, B))
              return EXIT_FAILURE;
		
            // separa SL->Int em LU
            if (SL_triangulariza(SL, SL->Int, B)) 
              return EXIT_FAILURE;

            SL_substituicao(SL, B, pol);

            SL_printMatrix(f_out, pol, SL->n, 1);
	    
            if (SL_ajusteDeCurvas(SL, i, B, lookup))
              return EXIT_FAILURE;

            // separa SL->Ajc em LU
            if (SL_triangulariza_otimiz(SL, SL->Ajc, B))
              return EXIT_FAILURE;

            SL_substituicao(SL, B, pol);
		
            SL_printMatrix(f_out, pol, SL->n, 1);
        }

        free(lookup);
        free(pol);
        free(B);
        SL_libera(SL);
    }
#else
    LIKWID_MARKER_INIT;   
    while (!feof(stdin))
    {
        SL = SL_leitura();
        if (!SL) break;
        
        double *B = SL_alocaMatrix(1, SL->n);
        if (!B) return EXIT_FAILURE;

        double *pol = SL_alocaMatrix(1, SL->n);
        if (!pol) return EXIT_FAILURE;

        for (int i=0; i<SL->m; ++i) {
            LIKWID_MARKER_START("Interpolacao");
	    if (SL_interpolacao(SL, i, B)) return EXIT_FAILURE;
	    LIKWID_MARKER_STOP("Interpolacao");

            // separa SL->Int em LU
            LIKWID_MARKER_START("Triangulariza");
            if (SL_triangulariza(SL, SL->Int, B)) return EXIT_FAILURE;
            LIKWID_MARKER_STOP("Triangulariza");
            SL_substituicao(SL, B, pol);

            SL_printMatrix(f_out, pol, SL->n, 1);
	    
            LIKWID_MARKER_START("AjusteDeCurvas");
            if (SL_ajusteDeCurvas(SL, i, B)) return EXIT_FAILURE;
            LIKWID_MARKER_STOP("AjusteDeCurvas");

            // separa SL->Ajc em LU
            LIKWID_MARKER_START("TriangularizaOtimiz");
            if (SL_triangulariza_otimiz(SL, SL->Ajc, B)) return EXIT_FAILURE;
            LIKWID_MARKER_STOP("TriangularizaOtimiz");

            SL_substituicao(SL, B, pol);
		
            SL_printMatrix(f_out, pol, SL->n, 1);
        }
        free(B);
        SL_libera(SL);
    }
    LIKWID_MARKER_CLOSE;
#endif

    if (f_out != stdout) fclose(f_out);
    
    return EXIT_SUCCESS;
}
