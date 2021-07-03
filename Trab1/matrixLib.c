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
  float **newMatrix = calloc(1, n*sizeof(float*) + n*n*sizeof(float));
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

static float** copiaMatrix(t_matrix *Mat) {

	float **aux = alocaMatrix(Mat->n);
	if (!aux)
    	return NULL;

	for (unsigned int i=0; i<Mat->n; ++i)
		for (unsigned int j=0; j<Mat->n; ++j)
			aux[i][j] = Mat->A[i][j];
	
	return aux;
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
			printf("%-10g ",matrix[i][j]);
		printf("\n");
	}
}

// Triangulariza a matriz Mat->A n x n
// Separa Mat->A em L e U
// Retorna 0 se sucesso e -1 caso contrario
int triangularizaMatrix(t_matrix *Mat, int pivotP, double *tTotal) {
    
	float **copia = copiaMatrix(Mat);
    if (!copia) {
        perror("Erro Eliminacao Gauss: falha ao copiar sistema linear");
        return -1;
    }
    
    //*tTotal = timestamp();
    
    // Transforma a matriz em uma triangular com pivoteamento parcial
    for (int i=0; i<Mat->n; i++) {
        //if (pivotP) {
		//pivo = maxValue(copia,i);
        //if (pivo != i)
            //trocaLinha(copia,i,pivo);
		//}

		Mat->U[i][i] = 1;

        for (int j=i+1; j<Mat->n; j++) {
            double m = copia[j][i] / copia[i][i];
            copia[j][i] = 0.0f;
			Mat->U[j][i] = m;
            for (int k=i+1; k<Mat->n; k++)
                copia[j][k] -= copia[i][k] * m;
        }
    }

	Mat->L = copia;
    //*tTotal = timestamp() - *tTotal;
    return 0;
}