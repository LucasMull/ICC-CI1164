/**
 * Luan Machado Bernardt | GRR20190363
 * Lucas Müller          | GRR20197160
 */

#include "matrixLib.h"
#include <stdlib.h>
#include <stdio.h>

int main () {

    t_matrix *newMatrix;
    double tempo;

    newMatrix = readMatrix();
    triangularizaMatrix(newMatrix,0,&tempo);
    printf("##################\n");
    printMatrix(newMatrix->A,newMatrix->n);
    printf("##################\n");
    printMatrix(newMatrix->U,newMatrix->n);
    printf("##################\n");
    printMatrix(newMatrix->L,newMatrix->n);
    return 0;
}
