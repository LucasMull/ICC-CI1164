#include "matrixLib.h"
#include <stdio.h>
#include <stdlib.h>

// Aloca memoria para a struct t_matrix com tamanho n
// Retorna um ponteiro para t_matrix ou NULL se houver erro
t_matrix *alocaMatrix(unsigned int n) {

	t_matrix *newMatrix = malloc(sizeof(t_matrix));
	if (!newMatrix)
		return NULL;

	newMatrix->A = malloc(n * sizeof(float *));
	if (!newMatrix->A)
		return NULL;


	for (unsigned int i=0; i<n; ++i) {
		newMatrix->A[i] = malloc(n * sizeof(float));
		if (!newMatrix->A[i])
			return NULL;
	}
	
	return newMatrix; 
}

// Le valores de stdin e preenche newMatrix
// Retorna ponteiro para t_matrix ou NULL se der erro
t_matrix *readMatrix() {

	unsigned int n;
	float num;

	scanf("%u",&n);
	t_matrix *newMatrix = alocaMatrix(n);
	if (!newMatrix) {
		perror("Erro ao alocar matriz");
		return NULL;
	} 

	newMatrix->n = n;

	for (unsigned int i=0; i<n; ++i)
		for (unsigned int j=0; j<n; ++j) {
			scanf("%f",&num);
			newMatrix->A[i][j] = num;
		}

	return newMatrix;
}

void printMatrix(float **matrix, int n) {

	for (unsigned int i=0; i<n; ++i) {
		for (unsigned int j=0; j<n; ++j)
			printf("%g ",matrix[i][j]);
		printf("\n");
	}
}