/**
 * Luan Machado Bernardt | GRR20190363
 * Lucas Müller          | GRR20197160
 */

#include <stdlib.h>
#include <stdio.h>

#include "matrixLib.h"

int main () {

    t_matrix *Mat;
    float **matId;
    double tempoTri, tempoUx, tempoLy;

    Mat = readMatrix();

    triangularizaMatrix(Mat,0,&tempoTri);
    printf("%.20s\n", "Original ##################");
    printMatrix(Mat->A,Mat->n);

    matId = geraIdentidade(Mat->n);

    geraInversa(Mat,matId,&tempoLy,&tempoUx);
    printf("%.20s\n", "Inversa ##################");
    printMatrix(Mat->Inv,Mat->n);
    printf("###############\n");
    printf("# Tempo Triangularização: %g ms\n",tempoTri);
    printf("# Tempo cálculo de Y: %g ms\n",tempoLy);
    printf("# Tempo cálculo de X: %g ms\n",tempoUx);
    for (unsigned int i=0; i<Mat->n; ++i) {
        printf("# Norma L2 dos residuos (%d): ", i);
        printf("%g\n",normaL2Residuo(Mat,matId[i],i));
    }
     printf("###############\n");
    
    free(matId);
    limpaStruct(Mat);

    return EXIT_SUCCESS;
}
