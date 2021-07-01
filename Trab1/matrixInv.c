#include "matrixLib.h"
#include <stdlib.h>
#include <stdio.h>

int main () {

    t_matrix *newMatrix;

    newMatrix = readMatrix();
    printMatrix(newMatrix->A,newMatrix->n);
    return 0;
}