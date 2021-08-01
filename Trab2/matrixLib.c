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
float* SL_alocaMatrix(unsigned int n) {

  float *newMatrix = calloc(1, n*n*sizeof(float));
  if (!newMatrix) {
    perror("Falha ao alocar matriz\n");
    return NULL;
  }
  return newMatrix;
}

/*!
  \brief Copia sistema linear

  \param SL o sistema linear a ser copiado

  \return cópia da matriz A do sistema. NULL se houve erro de alocação
*/
static float* copiaA(t_sist *SL) {

	float *copia = SL_alocaMatrix(SL->n);
	if (!copia) return NULL;

	for (unsigned int i=0; i<SL->n; ++i)
		for (unsigned int j=0; j<SL->n; ++j)
			copia[SL->n*i+j] = SL->A[SL->n*i+j];
	return copia;
}

/*!
  \brief Aloca memória para o sistema linear
  \todo em caso de falha de alocação, é necessário limpar a memória previamente alocada (se houver)

  \param n tamanho da matriz

  \return ponteiro para t_sist. NULL se houve erro de alocação
*/
t_sist *SL_aloca(unsigned int n) {

	t_sist *newSL = malloc(sizeof(t_sist));
	if (!newSL) return NULL;

	newSL->A = SL_alocaMatrix(n);
	if (!newSL->A) {
    free(newSL);
    return NULL;
  }
	newSL->Inv = SL_alocaMatrix(n);
	if (!newSL->Inv) {
    free(newSL->A);
    free(newSL);
    return NULL;
  }

	newSL->L = SL_alocaMatrix(n);
	if (!newSL->L) {
    free(newSL->Inv);
    free(newSL->A);
    free(newSL);
    return NULL;
  }

	return newSL; 
}

/*!
  \brief Libera recursos alocados por alocaStruct()

  \param SL o sistema linear a ser liberado da memória
*/
void SL_libera(t_sist *SL) {

  free(SL->A);
  free(SL->Inv);
  free(SL->Id);
  free(SL->L);
  free(SL->U);
  free(SL);
}

/*!
  \brief Essa função calcula a norma L2 do resíduo de uma matriz 

  \param SL ponteiro para o sistema linear
  \param matId ponteiro para a matriz identidade
  \param col Coluna de matInv a ser multiplicada

  \return Norma L2 do resíduo.
*/
float SL_normaL2Residuo(t_sist *SL, float *matId, unsigned int col) {

  float sum = 0.0f;
  float res;
    
    for (unsigned int i=0; i<SL->n; ++i) {
        res = 0.0f;
        for (int j=0; j<SL->n; ++j)
            res += SL->A[SL->n*i+j] * SL->Inv[SL->n*j+col];
        res = matId[i] - res;
        sum += res*res;
    }
    return (sqrtf(sum));
}

/*!
  \brief Le valores de stdin para preencher t_sist

  \return ponteiro para t_sist. NULL se houve erro de alocação
*/
t_sist *SL_leitura() {

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
    t_sist *newSL = SL_aloca(n);
    if (!newSL) {
        fputs("Não foi possível alocar 'newSL'\n", stderr);
        return NULL;
    }
    newSL->n = n;

    // extrai as linhas contendo os elementos da matriz
    char *valorAtual, *valorProx;
    for (unsigned int i=0; i < newSL->n; ++i) {
        valorAtual = ln;
        if (!fgets(ln, sizeof(ln), stdin)) {
            perror("Falha de leitura");
            return NULL;
        }

        for (unsigned int j=0; j < newSL->n; ++j) {
            newSL->A[n*i+j] = strtof(valorAtual, &valorProx);
            if (!valorProx) {
                fputs("Não foi possível obter o próximo valor\n", stderr);
                return NULL;
            }
            valorAtual = valorProx;
        }
    }
    
    // "consome" próxima linha (vazia ou parada em EOF)
    fgets(ln, sizeof(ln), stdin);

    return newSL;
}

/*!
  \brief Imprime matriz

  \param matrix matriz a ser impressa
  \param n tamanho da matriz
*/
void SL_printMatrix(FILE *f_out, float *matrix, int n) {

	fprintf(f_out,"\n");
  for (unsigned int i=0; i<n; ++i) {
		for (unsigned int j=0; j<n; ++j)
			fprintf(f_out,"%-10g ",matrix[n*i+j]);
		fprintf(f_out,"\n");
	}
  fprintf(f_out,"\n");
}

/*!
  \brief Encontra o maior valor em uma coluna da matriz

  \param matrix a matriz
  \param n dimensao da matriz
  \param i coluna
  \return indice da coluna com max
*/
static unsigned int maxValue (float *matrix, unsigned int n, unsigned int i) {

    unsigned int max = i;
    
    for (int j=i+1; j<n; j++) {
        if (fabs(matrix[n*j+i]) > fabs(matrix[max*j+i]))
          max = j;
    }
    
    return max;
}

/*!
  \brief Troca linhas de SL->A

  \param matrix a matriz
  \param i linha a ser trocada com j
  \param j linha a ser trocada com i
*/
static void trocaLinha (float *A, unsigned int i, unsigned int j) {

    float aux;

    aux = A[i];
    A[i] = A[j];
    A[j] = aux;
}

/*!
  \brief Calcula determinante a partir da U

  \param U matriz U
  \param n norma da matriz
  \return determinante
*/
static float determinanteU (float *U, unsigned int n) {

    float det = U[0];
    for (int i=1; i < n; ++i)
        det *= U[n*i+i];
    return det;
}

/*!
  \brief Triangulariza a matriz SL->A de norma n
  \note separa SL->A em L e U

  \param SL o sistema linear a ser triangularizado
  \param pivotP pivo parcial
  \param tTotal recebe tempo decorrido para cálculo
  \return 0 se sucesso e -1 em caso de falha
*/
int SL_triangulariza(t_sist *SL, int pivotP, double *tTotal) {
    
    float *copia = copiaA(SL);
    if (!copia) {
        perror("Erro Triangularizacao: falha ao copiar matriz");
        return -1;
    }
    
    *tTotal = timestamp();
    
    if (pivotP) {
      // Transforma a matriz em uma triangular com pivoteamento parcial
      unsigned int pivo;
      
      for (int i=0; i<SL->n; i++) 
      {
          pivo = maxValue(copia,SL->n,i);
          if (pivo != i) {
              trocaLinha(copia,i,pivo);
              trocaLinha(SL->A,i,pivo);
              trocaLinha(SL->Id,i,pivo);
              trocaLinha(SL->L,i,pivo);
          }

          SL->L[SL->n*i+i] = 1;
          for (int j=i+1; j<SL->n; j++) {
              double m = copia[SL->n*j+i] / copia[SL->n*i+i];
              copia[SL->n*j+i] = 0.0f;
              SL->L[SL->n*j+i] = m;
              for (int k=i+1; k<SL->n; k++)
                  copia[SL->n*j+k] -= copia[SL->n*i+k] * m;
          }
      }
    } else {
      // Transforma a matriz em uma triangular sem pivoteamento
      for (int i=0; i<SL->n; i++) 
      {
          SL->L[SL->n*i+i] = 1;
          for (int j=i+1; j<SL->n; j++) {
              double m = copia[SL->n*j+i] / copia[SL->n*i+i];
              copia[SL->n*j+i] = 0.0f;
              SL->L[SL->n*j+i] = m;
              for (int k=i+1; k<SL->n; k++)
                  copia[SL->n*j+k] -= copia[SL->n*i+k] * m;
          }
      }
    }

    *tTotal = timestamp() - *tTotal;

    SL->U = copia;

    // checar se matriz é inversível
    if (0.0f == determinanteU(SL->U, SL->n)) {
      fprintf(stderr, "Erro Triangularizacao: Matriz não é inversível, Det = 0\n");
      return -1;
    }

    return 0;
}

/*!
  \brief Gera uma matriz identidade a partir de uma norma n

  \param n tamanho da matriz
  \return matriz identidade n x n, NULL em caso de falha
*/
float *SL_geraIdentidade(unsigned int n) {

	float *matId = SL_alocaMatrix(n);
	if (!matId) {
		perror("Erro ao alocar matriz identidade");
		return NULL;
	}
	
	for (int i=0; i<n; ++i)
		matId[n*i+i] = 1.0f;

	return matId;
}

/*!
  \brief Gera a inversa de SL->A

  Calcula o sistema Ly=I e Ux=y e armazena em Mat->Inv
  \param Mat matriz original a ser invertida
  \param matId matriz identidade
  \param timeLy tempo para calculo de Ly=I
  \param timeUx tempo para calculo de Ux=y
*/
void SL_geraInversa(t_sist *SL, double *timeLy, double *timeUx) {
	
  double timeSum = 0.0f;

  // Calcula Ly=I
	
  for (int k=0; k<SL->n; ++k) 
  {
		*timeLy = timestamp();
    for (int i=0; i<SL->n; ++i) {
			SL->Inv[SL->n*i+k] = SL->Id[SL->n*k+i];
			for (int j=i-1; j>=0; --j)
				SL->Inv[SL->n*i+k] -= SL->L[SL->n*i+j] * SL->Inv[SL->n*j+k];
			SL->Inv[SL->n*i+k] /= SL->L[SL->n*i+i];
		}
	  *timeLy = timestamp() - *timeLy;
    timeSum += *timeLy;
  }
  
  *timeLy = timeSum/SL->n;
  timeSum = 0.0f;

  // Calcula Ux=y
  
	for (int k=0; k<SL->n; ++k) 
  {
    *timeUx = timestamp();
		for (int i=SL->n-1; i>=0; --i) {
			for (int j=i+1; j<SL->n; ++j)
				SL->Inv[SL->n*i+k] -= SL->U[SL->n*i+j] * SL->Inv[SL->n*j+k];
			SL->Inv[SL->n*i+k] /= SL->U[SL->n*i+i];
		}
	  *timeUx = timestamp() - *timeUx;
    timeSum += *timeUx;
  }
  
  *timeUx = timeSum/SL->n;
}
