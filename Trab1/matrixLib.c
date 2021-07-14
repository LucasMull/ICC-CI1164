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
  if (!newMatrix) {
    perror("Falha ao alocar matriz\n");
    return NULL;
  }
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

#if 0
	newMatrix->U = alocaMatrix(n);
	if (!newMatrix->U) return NULL;
#endif
	return newMatrix; 
}

/*!
  \brief Libera recursos alocados por alocaStruct()

  \param Mat matriz a ser liberada da memória
*/
void limpaStruct(t_matrix *Mat) {

  free(Mat->A);
  free(Mat->Inv);
  free(Mat->Id);
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

    char ln[1024];
    
    // extrai a linha contendo a ordem da matriz
    if (!fgets(ln, sizeof(ln), stdin)) {
        perror("Falha de leitura");
        return NULL;
    }
    unsigned int n = strtoul(ln, NULL, 10);
    if (!n) {
        fputs("Não foi possível obter ordem de matriz\n", stderr);
        return NULL;
    }

    // realiza alocação do Sistema Linear
    t_matrix *newMatrix = alocaStruct(n);
    if (!newMatrix) {
        fputs("Não foi possível alocar 'newMatrix'\n", stderr);
        return NULL;
    }
    newMatrix->n = n;

    // extrai as linhas contendo os elementos da matriz
    char *valorAtual, *valorProx;
    for (unsigned int i=0; i < newMatrix->n; ++i) {
        valorAtual = ln;
        if (!fgets(ln, sizeof(ln), stdin)) {
            perror("Falha de leitura");
            return NULL;
        }

        for (unsigned int j=0; j < newMatrix->n; ++j) {
            newMatrix->A[i][j] = strtof(valorAtual, &valorProx);
            if (!valorProx) {
                fputs("Não foi possível obter o próximo valor\n", stderr);
                return NULL;
            }
            valorAtual = valorProx;
        }
    }
    
    // "consome" próxima linha (vazia ou parada em EOF)
    fgets(ln, sizeof(ln), stdin);

    return newMatrix;
}

/*!
  \brief Imprime matriz

  \param matrix matriz a ser impressa
  \param n tamanho da matriz
*/
void printMatrix(FILE *f_out, float **matrix, int n) {

	fprintf(f_out,"\n");
  for (unsigned int i=0; i<n; ++i) {
		for (unsigned int j=0; j<n; ++j)
			fprintf(f_out,"%-10g ",matrix[i][j]);
		fprintf(f_out,"\n");
	}
  fprintf(f_out,"\n");
}

/*!
  \brief Encontra o maior valor em uma coluna da matriz Mat

  \param Mat struct com a matriz
  \param n dimensao da matriz Mat
  \param i coluna
  \return indice da coluna com max
*/
static unsigned int maxValue (float **Mat, unsigned int n, unsigned int i) {

    unsigned int max = i;
    
    for (int j=i+1; j<n; j++) {
        if (fabs(Mat[j][i]) > fabs(Mat[max][i]))
          max = j;
    }
    
    return max;
}

/*!
  \brief Troca linhas da matriz Mat->A

  \param Mat matriz
  \param i linha a ser trocada com j
  \param j linha a ser trocada com i
*/
static void trocaLinha (float **Mat, unsigned int i, unsigned int j) {

    float *aAux;

    aAux = Mat[i];
    Mat[i] = Mat[j];
    Mat[j] = aAux;
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
        perror("Erro Triangularizacao: falha ao copiar matriz");
        return -1;
    }
    
    *tTotal = timestamp();
    
    if (pivotP) {
      // Transforma a matriz em uma triangular com pivoteamento parcial
      unsigned int pivo;
      
      for (int i=0; i<Mat->n; i++) 
      {
          pivo = maxValue(copia,Mat->n,i);
          if (pivo != i) {
              trocaLinha(copia,i,pivo);
              trocaLinha(Mat->A,i,pivo);
              trocaLinha(Mat->Id,i,pivo);
              trocaLinha(Mat->L,i,pivo);
          }

          Mat->L[i][i] = 1;
          for (int j=i+1; j<Mat->n; j++) {
              double m = copia[j][i] / copia[i][i];
              copia[j][i] = 0.0f;
              Mat->L[j][i] = m;
              for (int k=i+1; k<Mat->n; k++)
                  copia[j][k] -= copia[i][k] * m;
          }
      }
    } else {
      // Transforma a matriz em uma triangular sem pivoteamento
      for (int i=0; i<Mat->n; i++) 
      {
          Mat->L[i][i] = 1;
          for (int j=i+1; j<Mat->n; j++) {
              double m = copia[j][i] / copia[i][i];
              copia[j][i] = 0.0f;
              Mat->L[j][i] = m;
              for (int k=i+1; k<Mat->n; k++)
                  copia[j][k] -= copia[i][k] * m;
          }
      }
    }

    *tTotal = timestamp() - *tTotal;

    /// @todo Mat->U é alocado em alocaStruct, é preciso copiar conteúdo de 'copia' para Mat->U, ou então não alocar memória para Mat->U
    Mat->U = copia;

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
  \param timeLy tempo para calculo de Ly=I
  \param timeUx tempo para calculo de Ux=y
*/
void geraInversa(t_matrix *Mat, double *timeLy, double *timeUx) {
	
  double timeSum = 0.0f;

  // Calcula Ly=I
	
  for (int k=0; k<Mat->n; ++k) 
  {
		*timeLy = timestamp();
    for (int i=0; i<Mat->n; ++i) {
			Mat->Inv[i][k] = Mat->Id[k][i];
			for (int j=i-1; j>=0; --j)
				Mat->Inv[i][k] -= Mat->L[i][j] * Mat->Inv[j][k];
			Mat->Inv[i][k] /= Mat->L[i][i];
		}
	  *timeLy = timestamp() - *timeLy;
    timeSum += *timeLy;
  }
  
  *timeLy = timeSum/Mat->n;
  timeSum = 0.0f;

  // Calcula Ux=y
  
	for (int k=0; k<Mat->n; ++k) 
  {
    *timeUx = timestamp();
		for (int i=Mat->n-1; i>=0; --i) {
			for (int j=i+1; j<Mat->n; ++j)
				Mat->Inv[i][k] -= Mat->U[i][j] * Mat->Inv[j][k];
			Mat->Inv[i][k] /= Mat->U[i][i];
		}
	  *timeUx = timestamp() - *timeUx;
    timeSum += *timeUx;
  }
  
  *timeUx = timeSum/Mat->n;
}
