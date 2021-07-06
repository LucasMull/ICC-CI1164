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
    double tempo;

    Mat = readMatrix();

    triangularizaMatrix(Mat,0,&tempo);
    printf("%.20s\n", "Original ##################");
    printMatrix(Mat->A,Mat->n);

    printf("%.20s\n", "U ##################");
    printMatrix(Mat->U,Mat->n);

    printf("%.20s\n", "L ##################");
    printMatrix(Mat->L,Mat->n);

    matId = geraIdentidade(Mat->n);
    printf("%.20s\n", "Identidade ##################");
    printMatrix(matId,Mat->n);

    geraInversa(Mat,matId);
    printf("%.20s\n", "Inversa ##################");
    printMatrix(Mat->Inv,Mat->n);

    float **res = alocaMatrix(Mat->n);
    float norma = normaL2Residuo(Mat,matId,res);
    printf("%.20s\n", "Residuo ##################");
    printMatrix(res,Mat->n);
    printf("Norma L2: %g\n", norma);

    free(res);
    free(matId);
    limpaStruct(Mat);

    return EXIT_SUCCESS;
}
