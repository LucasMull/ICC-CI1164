/**
 * Luan Machado Bernardt | GRR20190363
 * Lucas Müller          | GRR20197160
 */

#include "matrixLib.h"
#include <stdio.h>
#include <stdlib.h>

/*!
  \brief Aloca matriz

  \param n tamanho da matriz

  \return ponteiro para matriz. NULL se houve erro de alocação
  */
static float** alocaMatrix(unsigned int n) {

  /* efetua alocação de matriz em 1D para facilitar limpeza */
  float **newMatrix = malloc(n*sizeof(float*) + n*n*sizeof(float));
  if (!newMatrix)
    return NULL;

  /* inicializa cada ponteiro para seu bloco de memória
   *        consecutivo alocado */
  float *addr = (float*)(newMatrix + n);
  for (unsigned int i=0; i < n; ++i) {
    newMatrix[i] = addr;
    addr += n;
  }
  return newMatrix;
}

// Aloca memoria para a struct t_matrix com tamanho n
// Retorna um ponteiro para t_matrix ou NULL se houver erro
t_matrix *alocaStruct(unsigned int n) {

	t_matrix *newMatrix = malloc(sizeof(t_matrix));
	if (!newMatrix)
		return NULL;

	newMatrix->A = alocaMatrix(n);
	if (!newMatrix->A)
		return NULL;
#if 0	
	newMatrix->L = malloc(n * sizeof(float));
	if (!newMatrix->L)
		return NULL;
	newMatrix->U = malloc(n * sizeof(float));
	if (!newMatrix->U)
		return NULL;
#endif
	newMatrix->L = alocaMatrix(n);
	if (!newMatrix->L)
		return NULL;
	newMatrix->U = alocaMatrix(n);
	if (!newMatrix->U)
		return NULL;

	return newMatrix; 
}

// Le valores de stdin e preenche newMatrix
// Retorna ponteiro para t_matrix ou NULL se der erro
t_matrix *readMatrix() {

	unsigned int n;
	float num;

	scanf("%u",&n);
	t_matrix *newMatrix = alocaStruct(n);
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
