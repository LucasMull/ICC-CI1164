/**
 * Luan Machado Bernardt | GRR20190363
 * Lucas Müller          | GRR20197160
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> // memcpy()
#include <math.h>

#include "utils.h"
#include "matrixLib.h"

/*!
  \brief Aloca matriz

  \param n número de valores tabelados
  \param m número de funções tabeladas

  \return ponteiro para matriz. NULL se houve erro de alocação
*/
float* SL_alocaMatrix(unsigned int n, unsigned int m) {

  float *newMatrix = calloc(1, n*m*sizeof(float));
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

	float *copia = SL_alocaMatrix(SL->n, SL->m);
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
t_sist *SL_aloca(unsigned int n, unsigned int m) {

	t_sist *newSL = malloc(sizeof(t_sist));
	if (!newSL) return NULL;

	newSL->A = SL_alocaMatrix(n, m);
	if (!newSL->A) {
    free(newSL);
    return NULL;
  }

	newSL->L = SL_alocaMatrix(n, n);
	if (!newSL->L) {
    free(newSL->A);
    free(newSL);
    return NULL;
  }

  newSL->x = calloc(1, n * sizeof(float));
  if (!newSL->x) {
    free(newSL->L);
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
  free(SL->L);
  if (SL->U)
    free(SL->U);
  if (SL->Int)
    free(SL->Int);
  free(SL->x);
  free(SL);
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

    unsigned int n=0, m=0;
    sscanf(ln, "%d %d", &n, &m);
    if (!n || !m) {
        perror("Não foi possível obter ordem de matriz");
        return NULL;
    }

    // realiza alocação do Sistema Linear
    t_sist *newSL = SL_aloca(n, m);
    if (!newSL) {
        fputs("Não foi possível alocar 'newSL'\n", stderr);
        return NULL;
    }
    newSL->n = n;
    newSL->m = m;

    char *valorAtual, *valorProx;
    if (!fgets(ln, sizeof(ln), stdin)) {
        perror("#1 Falha de leitura");
        return NULL;
    }
    valorAtual = ln;

    for (unsigned int i=0; i < newSL->n; ++i) {
      newSL->x[i] = strtof(valorAtual, &valorProx);
      if (!valorProx) {
          fputs("Não foi possível obter o próximo valor\n", stderr);
          return NULL;
      }
      valorAtual = valorProx;
    }

    // extrai as linhas contendo os elementos da matriz
    for (unsigned int i=0; i < newSL->m; ++i) {
        valorAtual = ln;
        if (!fgets(ln, sizeof(ln), stdin)) {
            perror("#2 Falha de leitura");
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
  \param n número de valores tabelados
  \param m número de funções tabeladas
*/
void SL_printMatrix(FILE *f_out, float *matrix, unsigned int n, unsigned int m) {

	fprintf(f_out,"\n");
  for (unsigned int i=0; i<m; ++i) {
		for (unsigned int j=0; j<n; ++j)
			fprintf(f_out,"%-10g ",matrix[n*i+j]);
		fprintf(f_out,"\n");
	}
  fprintf(f_out,"\n");
}

float *SL_interpolacao(t_sist *SL, unsigned int row) {

  float *mat = SL_alocaMatrix(SL->n, SL->n);
  if (!mat) return NULL;

  for (unsigned int i=0; i<SL->n; ++i)
    for (unsigned int j=0; j<SL->n; ++j)
      mat[SL->n*i+j] = powf(SL->x[i], (float)j);
  SL->Int = mat;
  SL_triangulariza(SL, &(double){0.0});

  float *pol = SL_alocaMatrix(SL->n, 1); // polinomio
  if (!pol) return NULL;

  for (int i=0; i<SL->n; ++i) {
    pol[i] = SL->A[SL->m*row+i];
    for (int j=i-1; j>=0; --j)
      pol[i] -= SL->L[SL->n*i+j] * pol[j];
    pol[i] /= SL->L[SL->n*i+i];
  }

  for (int i=SL->n-1; i>=0; --i) {
    for (int j=i+1; j<SL->n; ++j)
      pol[i] -= SL->U[SL->n*i+j] * pol[j];
    pol[i] /= SL->U[SL->n*i+i];
  }

  return pol;
}

/*!
  \brief Encontra o maior valor em uma coluna da matriz

  \param matrix a matriz
  \param n dimensao da matriz
  \param j coluna
  \return indice da coluna com max
*/
static unsigned int maxValue (float *matrix, unsigned int n, unsigned int j) {

    unsigned int max = j;
    for (int i=max+1; i<n; i++) {
        if (fabs(matrix[n*i+j]) > fabs(matrix[n*max+j]))
          max = i;
    }
    return max;
}

/*!
  \brief Troca elementos

  \param a elemento a
  \param b elemento b
*/
static void trocaElemento (float *a, float *b)
{
  float aux;

  aux = *a;
  *a = *b;
  *b = aux;
}

/*!
  \brief Troca linhas de SL->A

  \param matrix a matriz
  \param i linha a ser trocada com j
  \param j linha a ser trocada com i
*/
static void trocaLinha (float *matrix, unsigned int i, unsigned int j, unsigned int n) {

    float *aux = malloc(n*sizeof(float));;
    if(!aux) {
      perror("Sem memória");
      return;
    }

    memcpy(aux, &matrix[n*i], n*sizeof(float));
    memcpy(&matrix[n*i], &matrix[n*j], n*sizeof(float));
    memcpy(&matrix[n*j], aux, n*sizeof(float));

    free(aux);
}

/*!
  \brief Triangulariza a matriz SL->A de norma n
  \note separa SL->A em L e U

  \param SL o sistema linear a ser triangularizado
  \param tTotal recebe tempo decorrido para cálculo
  \return 0 se sucesso e -1 em caso de falha
*/
int SL_triangulariza(t_sist *SL, double *tTotal) {
    
    float *copia = malloc(SL->n * SL->n * sizeof(float));
    memcpy(copia, SL->Int, SL->n * SL->n * sizeof(float));

    if (!copia) {
        perror("Erro Triangularizacao: falha ao copiar matriz");
        return -1;
    }
    
    *tTotal = timestamp();
    
    // Transforma a matriz em uma triangular com pivoteamento parcial
    unsigned int pivo;
    for (int i=0; i<SL->n; i++) 
    {
        pivo = maxValue(copia,SL->n,i);
        if (pivo != i) {
            trocaElemento(&SL->A[SL->m*i+pivo-1], &SL->A[SL->m*i+i-1]);
            trocaLinha(copia,i,pivo,SL->n);
            trocaLinha(SL->Int,i,pivo,SL->n);
            trocaLinha(SL->L,i,pivo,SL->n);
        }

        SL->L[SL->n*i+i] = 1.0f;
        for (int j=i+1; j<SL->n; j++) {
            double m = copia[SL->n*j+i] / copia[SL->n*i+i];
            copia[SL->n*j+i] = 0.0f;
            SL->L[SL->n*j+i] = m;
            for (int k=i+1; k<SL->n; k++)
                copia[SL->n*j+k] -= copia[SL->n*i+k] * m;
        }
    }

    *tTotal = timestamp() - *tTotal;

    SL->U = copia;
    return 0;
}
