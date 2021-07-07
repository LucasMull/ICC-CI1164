/**
 * Luan Machado Bernardt | GRR20190363
 * Lucas Müller          | GRR20197160
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "matrixLib.h"

/*!
  \brief Aloca matriz

  \param n tamanho da matriz

  \return ponteiro para matriz. NULL se houve erro de alocação
*/
float** alocaMatrix(unsigned int n) {

  /* efetua alocação de matriz em 1D para facilitar limpeza */
  float **newMatrix = calloc(1, n*sizeof(float*) + n*n*sizeof(float));
  if (!newMatrix) return NULL;

  /* inicializa cada ponteiro para seu bloco de memória
   *        consecutivo alocado */
  float *addr = (float*)(newMatrix + n);
  for (unsigned int i=0; i < n; ++i) {
    newMatrix[i] = addr;
    addr += n;
  }
  return newMatrix;
}

/*!
  \brief Copia matriz

  \param Mat a matriz que será copiada

  \return cópia da matriz. NULL se houve erro de alocação
*/
static float** copiaMatrix(t_matrix *Mat) {

	float **aux = alocaMatrix(Mat->n);
	if (!aux) return NULL;

	for (unsigned int i=0; i<Mat->n; ++i)
		for (unsigned int j=0; j<Mat->n; ++j)
			aux[i][j] = Mat->A[i][j];
	return aux;
}

/*!
  \brief Aloca memória para a struct t_matrix
  \todo em caso de falha de alocação, é necessário limpar a memória previamente alocada (se houver)

  \param n tamanho da matriz

  \return ponteiro para t_matriz. NULL se houve erro de alocação
*/
t_matrix *alocaStruct(unsigned int n) {

	t_matrix *newMatrix = malloc(sizeof(t_matrix));
	if (!newMatrix) return NULL;

	newMatrix->A = alocaMatrix(n);
	if (!newMatrix->A) return NULL;

	newMatrix->Inv = alocaMatrix(n);
	if (!newMatrix->Inv) return NULL;

#if 0
	newMatrix->L = malloc(n * sizeof(float));
	if (!newMatrix->L) return NULL;
	newMatrix->U = malloc(n * sizeof(float));
	if (!newMatrix->U) return NULL;
#endif

	newMatrix->L = alocaMatrix(n);
	if (!newMatrix->L) return NULL;
	newMatrix->U = alocaMatrix(n);
	if (!newMatrix->U) return NULL;
	return newMatrix; 
}

/*!
  \brief Libera recursos alocados por alocaStruct()

  \param Mat matriz a ser liberada da memória
*/
void limpaStruct(t_matrix *Mat) {

  free(Mat->A);
  free(Mat->Inv);
  free(Mat->L);
  free(Mat->U);
  free(Mat);
}

/*!
  \brief Essa função calcula a norma L2 do resíduo de uma matriz 

  \param Mat ponteiro para a matriz
  \param matId ponteiro para a matriz identidade
  \param col Coluna de matInv a ser multiplicada

  \return Norma L2 do resíduo.
*/
float normaL2Residuo(t_matrix *Mat, float *matId, unsigned int col) {

  float sum = 0.0f;
  float res;
    
    for (unsigned int i=0; i<Mat->n; ++i) {
        res = 0.0f;
        for (int j=0; j<Mat->n; ++j)
            res += Mat->A[i][j] * Mat->Inv[j][col];
        res = matId[i] - res;
        sum += powf(res,2.0f);
    }
    return (sqrtf(sum));
}

/*!
  \brief Le valores de stdin para preencher t_matrix

  \return ponteiro para t_matriz. NULL se houve erro de alocação
*/
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

/*!
  \brief Imprime matriz

  \param matrix matriz a ser impressa
  \param n tamanho da matriz
*/
void printMatrix(float **matrix, int n) {

	printf("\n");
  for (unsigned int i=0; i<n; ++i) {
		for (unsigned int j=0; j<n; ++j)
			printf("%-10g ",matrix[i][j]);
		printf("\n");
	}
  printf("\n");
}

/*!
  \brief Triangulariza a matriz Mat->a de norma n
  \note separa Mat->a em L e U

  \param Mat matriz a ser triangularizada
  \param pivotP pivo parcial
  \param tTotal recebe tempo decorrido para cálculo
  \return 0 se sucesso e -1 em caso de falha
*/
int triangularizaMatrix(t_matrix *Mat, int pivotP, double *tTotal) {
    
    float **copia = copiaMatrix(Mat);
    if (!copia) {
        perror("Erro Eliminacao Gauss: falha ao copiar sistema linear");
        return -1;
    }
    
    *tTotal = timestamp();
    
    // Transforma a matriz em uma triangular com pivoteamento parcial
    for (int i=0; i<Mat->n; i++) 
    {
#if 0
        if (pivotP) {
          pivo = maxValue(copia,i);
          if (pivo != i)
              trocaLinha(copia,i,pivo);
        }
#endif
        Mat->L[i][i] = 1;
        for (int j=i+1; j<Mat->n; j++) {
            double m = copia[j][i] / copia[i][i];
            copia[j][i] = 0.0f;
            Mat->L[j][i] = m;
            for (int k=i+1; k<Mat->n; k++)
                copia[j][k] -= copia[i][k] * m;
        }
    }

    /// @todo Mat->U é alocado em alocaStruct, é preciso copiar conteúdo de 'copia' para Mat->U, ou então não alocar memória para Mat->U
    Mat->U = copia;

    *tTotal = timestamp() - *tTotal;

    return 0;
}

/*!
  \brief Gera uma matriz identidade a partir de uma norma n

  \param n tamanho da matriz
  \return matriz identidade n x n, NULL em caso de falha
*/
float **geraIdentidade(unsigned int n) {

	float **matId = alocaMatrix(n);
	if (!matId) {
		perror("Erro ao alocar matriz identidade");
		return NULL;
	}
	
	for (int i=0; i<n; ++i)
		matId[i][i] = 1.0f;

	return matId;
}

/*!
  \brief Gera a inversa da matriz

  Calcula o sistema Ly=I e Ux=y e armazena em Mat->Inv
  \param Mat matriz original a ser invertida
  \param matId matriz identidade
*/
void geraInversa(t_matrix *Mat, float **matId, double *timeLy, double *timeUx) {
	
  // Calcula Ly=I
	*timeLy = timestamp();
  for (int k=0; k<Mat->n; ++k) 
  {
		for (int i=0; i<Mat->n; ++i) {
			Mat->Inv[i][k] = matId[k][i];
			for (int j=i-1; j>=0; --j)
				Mat->Inv[i][k] -= Mat->L[i][j] * Mat->Inv[j][k];
			Mat->Inv[i][k] /= Mat->L[i][i];
		}
	}
  *timeLy = timestamp() - *timeLy;

  // Calcula Ux=y
  *timeUx = timestamp();
	for (int k=0; k<Mat->n; ++k) 
  {
		for (int i=Mat->n-1; i>=0; --i) {
			for (int j=i+1; j<Mat->n; ++j)
				Mat->Inv[i][k] -= Mat->U[i][j] * Mat->Inv[j][k];
			Mat->Inv[i][k] /= Mat->U[i][i];
		}
	}
  *timeUx = timestamp() - *timeUx;
}