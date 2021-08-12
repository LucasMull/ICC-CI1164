/**
 * Luan Machado Bernardt | GRR20190363
 * Lucas Müller          | GRR20197160
 */

#include <stdlib.h>
#include <stdio.h>

#ifndef _NO_LIKWID
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_START(a)
#define LIKWID_MARKER_STOP(a)
#endif

#include "matrixLib.h"


int main (int argc, char **argv) {

    t_sist *SL, *Int, *Ajc;
    double *pol, *lookup;

    LIKWID_MARKER_INIT;   
    while (!feof(stdin))
    {
        SL = SL_leitura();
        if (!SL) break;

        Int = SL_aloca(SL->n, SL->n);
        if (!Int) return EXIT_FAILURE;

        Ajc = SL_aloca(SL->n, SL->n);
        if (!Ajc) return EXIT_FAILURE;

        pol = SL_alocaMatrix(1, SL->n);
        if (!pol) return EXIT_FAILURE;

        // guardar valores de x exp
        lookup = SL_alocaMatrix(SL->n, SL->n);
        if (!lookup) return EXIT_FAILURE;


        for (int i=0; i<SL->m; ++i) {
            LIKWID_MARKER_START("Interpolacao");
            if (SL_interpolacao(SL, Int, i)) return EXIT_FAILURE;
            LIKWID_MARKER_STOP("Interpolacao");
		
            // separa SL->Int em LU
            if (SL_triangulariza_otimiz(Int)) return EXIT_FAILURE;

            SL_substituicao(Int, pol);
            SL_printMatrix(stderr, pol, SL->n, 1);
	    
            LIKWID_MARKER_START("AjusteDeCurvas");
            if (SL_ajusteDeCurvas(SL, Ajc, i, lookup)) return EXIT_FAILURE;
            LIKWID_MARKER_STOP("AjusteDeCurvas");

            LIKWID_MARKER_START("TriangularizaOtimiz");
            if (SL_triangulariza_otimiz(Ajc)) return EXIT_FAILURE;
            LIKWID_MARKER_STOP("TriangularizaOtimiz");

            SL_substituicao(Ajc, pol);

            SL_printMatrix(stderr, pol, SL->n, 1);

            LIKWID_MARKER_START("Triangulariza");
            if (SL_triangulariza(Ajc)) return EXIT_FAILURE;
            LIKWID_MARKER_STOP("Triangulariza");
        }

        free(lookup);
        free(pol);
        SL_libera(Ajc);
        SL_libera(Int);
        SL_libera(SL);
    }
    LIKWID_MARKER_CLOSE;

    return EXIT_SUCCESS;
}
