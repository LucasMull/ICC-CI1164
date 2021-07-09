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
    char *output = NULL;
    int pivotP = 0;

    while ((opt = getopt(argc,argv,"po:")) > 0) {
        switch(opt) {
            case 'p':
                pivotP = 1;
                break;
            case 'o':
                output = optarg;
                break;
        }
    }

    Mat = readMatrix();
    
    printf("%.20s\n", "Original ##################");
    printMatrix(Mat->A,Mat->n);
    
    Mat->Id = geraIdentidade(Mat->n);
    triangularizaMatrix(Mat,pivotP,&tempoTri);

    geraInversa(Mat,Mat->Id,&tempoLy,&tempoUx);
    printf("%.20s\n", "Inversa ##################");
    printMatrix(Mat->Inv,Mat->n);

    printf("###############\n");
    printf("# Tempo Triangularização: %g ms\n",tempoTri);
    printf("# Tempo cálculo de Y: %g ms\n",tempoLy);
    printf("# Tempo cálculo de X: %g ms\n",tempoUx);
    for (unsigned int i=0; i<Mat->n; ++i) {
        printf("# Norma L2 dos residuos (%d): ", i);
        printf("%g\n",normaL2Residuo(Mat,Mat->Id[i],i));
    }
     printf("###############\n");

    limpaStruct(Mat);

    return EXIT_SUCCESS;
}
