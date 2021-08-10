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

#define FATOR_UNROLL 8

#define L_IDX(n,i,j) (((i)*((i)-1))/2+(j)-1)
#define U_IDX(n,i,j) ((n)*((i)-1)-((i)-2)*((i)-1)/2 + ((j)-(i)))

/*!
 * \brief Retorna base elevada a um expoente inteiro
 *
 * \param base o valor a ser exponencializado
 * \param expoente o expoente da base
 * \return o resultado da potencialização
 */
static inline double pot(double base, unsigned int expoente) {

  if (!expoente) return 1.0f;

  double aux = base;
  for (unsigned int i=1; i < expoente; ++i)
      base *= aux;
  return base;
}

/*!
 * \brief Substituição LU
 *
 * \param SL o sistema linear contendo L e U previamente obtidos
 * \param B termos independentes
 * \return retorna polinômio
 */
static inline double* substituicao(t_sist *SL, double *B) {

  double *pol = SL_alocaMatrix(1, SL->n);
  if (!pol) return NULL;

  for (int i=0; i<SL->n; ++i) {
      pol[i] = B[i];
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
  \brief Aloca matriz

  \param n número de valores tabelados
  \param m número de funções tabeladas

  \return ponteiro para matriz. NULL se houve erro de alocação
*/
double* SL_alocaMatrix(unsigned int n, unsigned int m) {

  double *newMatrix = calloc(1, n*m*sizeof(double));
  if (!newMatrix) {
    perror("Falha ao alocar matriz");
    return NULL;
  }
  return newMatrix;
}

/*!
  \brief Aloca memória para o sistema linear
  \todo em caso de falha de alocação, é necessário limpar a memória previamente alocada (se houver)

  \param n tamanho da matriz

  \return ponteiro para t_sist. NULL se houve erro de alocação
*/
t_sist *SL_aloca(unsigned int n, unsigned int m) {

	t_sist *newSL = calloc(1, sizeof(t_sist));
	if (!newSL) return NULL;

	newSL->A = SL_alocaMatrix(n, m);
	if (!newSL->A) {
      free(newSL);
      return NULL;
  }
  newSL->x = calloc(1, n * sizeof(double));
  if (!newSL->x) {
      free(newSL->L);
      free(newSL->A);
      free(newSL);
      return NULL;
  }

  newSL->n = n;
  newSL->m = m;

	return newSL; 
}

/*!
  \brief Libera recursos alocados por alocaStruct()

  \param SL o sistema linear a ser liberado da memória
*/
void SL_libera(t_sist *SL) {

  free(SL->A);
  free(SL->x);
  if (SL->L)   free(SL->L);
  if (SL->U)   free(SL->U);
  if (SL->Int) free(SL->Int);
  if (SL->Ajc) free(SL->Ajc);
  free(SL);
}

/*!
  \brief Le valores de stdin para preencher t_sist

  \return ponteiro para t_sist. NULL se houve erro de alocação
*/
t_sist *SL_leitura() {

    // extrai a linha contendo a ordem da matriz
    unsigned int n=0, m=0;
    if (EOF == scanf("%d %d%*[^\n]", &n, &m)) return NULL;
    scanf("%*c"); // consome '\n'
    if (!n || !m) {
        fputs("Não foi possível obter ordem de matriz\n", stderr);
        return NULL;
    }

    // realiza alocação do Sistema Linear
    t_sist *newSL = SL_aloca(n, m);
    if (!newSL) {
        fputs("Não foi possível alocar 'newSL'\n", stderr);
        return NULL;
    }

    char c=0;
    for (unsigned int i=0; i < newSL->n; ++i) {
        if (EOF == scanf("%lf%c", &newSL->x[i], &c)) {
            perror("Falha de leitura");
            return NULL;
        }
    }
    if (c != '\n') { // pular espaços extras se houver
        scanf("%*[^\n]");
        scanf("%*c");
    }

    // extrai as linhas contendo os elementos da matriz
    for (unsigned int i=0; i < newSL->m; ++i) {
        for (unsigned int j=0; j < newSL->n; ++j) {
            if (EOF == scanf("%lf%c", &newSL->A[n*i+j], &c)) {
                perror("Falha de leitura");
                return NULL;
            }
        }
        if (c != '\n') { // pular espaços extras se houver
            scanf("%*[^\n]");
            scanf("%*c");
        }
    }
    
    // "consome" próxima linha (vazia ou parada em EOF)
    scanf("%*[^\n]"); 
    scanf("%*c");

    return newSL;
}

/*!
  \brief Imprime matriz

  \param matrix matriz a ser impressa
  \param n número de valores tabelados
  \param m número de funções tabeladas
*/
void SL_printMatrix(FILE *f_out, double *matrix, unsigned int n, unsigned int m) {

  fputc('\n', f_out);
  for (unsigned int i=0; i<m; ++i) {
      for (unsigned int j=0; j<n; ++j)
          fprintf(f_out,"%-1.18g ",matrix[n*i+j]);
      fputc('\n', f_out);
	}
  fputc('\n', f_out);
}

/*!
 * \brief Realiza interpolação na matriz
 *
 * \param SL sistema linear
 * \param row a linha da matriz de entrada
 * \param B vetor de termos independentes do SL->Int
 * \return retorna polinômio gerado pela interpolação
 */
double* SL_interpolacao(t_sist *SL, unsigned int row, double *B) {

  if (!SL->Int) {
      SL->Int = SL_alocaMatrix(SL->n, SL->n);
      if (!SL->Int) return NULL;
      
      for (int i=0; i<(SL->n - (SL->n % FATOR_UNROLL)); i += FATOR_UNROLL)
          for (int j=0; j < SL->n; ++j) {
              SL->Int[SL->n*i+j] = pot(SL->x[i], j);
              SL->Int[SL->n*(i+1)+j] = pot(SL->x[i+1], j);
              SL->Int[SL->n*(i+2)+j] = pot(SL->x[i+2], j);
              SL->Int[SL->n*(i+3)+j] = pot(SL->x[i+3], j);
              SL->Int[SL->n*(i+4)+j] = pot(SL->x[i+4], j);
              SL->Int[SL->n*(i+5)+j] = pot(SL->x[i+5], j);
              SL->Int[SL->n*(i+6)+j] = pot(SL->x[i+6], j);
              SL->Int[SL->n*(i+7)+j] = pot(SL->x[i+7], j);
          }

      for (int i=(SL->n - (SL->n % FATOR_UNROLL)); i < SL->n; ++i)
          for (int j=0; j<SL->n; ++j)
              SL->Int[SL->n*i+j] = pot(SL->x[i], j);
  }

  // copia termos independentes para B
  memcpy(B, SL->A + SL->n*row, SL->n*sizeof(double));

  // separa SL->Int em LU
  if (SL_triangulariza(SL, SL->Int, B)) return NULL;

  // Realiza substituição LU e retorna polinomio resultante
  return substituicao(SL, B);
}

/*!
 * \brief Realiza ajuste de curvas na matriz
 *
 * \param SL sistema linear
 * \param row a linha da matriz de entrada
 * \param B termos independentes
 * \return retorna polinômio gerado pelo ajuste de curvas
 */
double *SL_ajusteDeCurvas(t_sist *SL, unsigned int row, double *B) {

  if (!SL->Ajc) {
      SL->Ajc = SL_alocaMatrix(SL->n, SL->n);
      if (!SL->Ajc) return NULL;
      
    // primeira linha da matriz (j: coluna, k: somatório)
    for (unsigned int j=0; j < SL->n; ++j)
        for (unsigned int k=0; k < SL->n; ++k)
            SL->Ajc[j] += pot(SL->x[k], j);

    // restante das linhas da matriz (i: linhas, j: coluna, k: somatório)
    for (unsigned int i=1; i < SL->n; ++i) {
        for (unsigned int j=0; j < SL->n-1; ++j)
            SL->Ajc[SL->n*i+j] = SL->Ajc[SL->n*(i-1)+(j+1)];
        for (unsigned int k=0; k < SL->n; ++k)
            SL->Ajc[SL->n*i+(SL->n-1)] += pot(SL->x[k], i) * pot(SL->x[k], SL->n-1);
    }
  }

  // zera os valores de B
  memset(B, 0, SL->n*sizeof(double));

  for (unsigned int i=0; i<SL->n; ++i)
      for (unsigned int j=0; j<SL->n; ++j)
          B[i] += SL->A[SL->n*row+j] * pot(SL->x[j], i);

  // separa SL->Ajc em LU
  if (SL_triangulariza(SL, SL->Ajc, B)) return NULL;

  // Realiza substituição LU e retorna polinomio resultante
  return substituicao(SL, B);
}

/*!
  \brief Encontra o maior valor em uma coluna da matriz

  \param matrix a matriz
  \param n dimensao da matriz
  \param j coluna
  \return indice da coluna com max
*/
static unsigned int maxValue(double *matrix, unsigned int n, unsigned int j) {

    unsigned int max = j;

    for (int i=max+1; i<n; i++)
        if (fabs(matrix[n*i+j]) > fabs(matrix[n*max+j]))
            max = i;
    return max;
}

/*!
  \brief Troca elementos

  \param a elemento a
  \param b elemento b
*/
static void trocaElemento(double *a, double *b)
{
  double aux;

  aux = *a;
  *a = *b;
  *b = aux;
}

/*!
  \brief Troca linhas de SL->A

  \param matrix a matriz
  \param i linha a ser trocada com j
  \param j linha a ser trocada com i
  \return 0 para sucesso e -1 para falha
*/
static int trocaLinha(double *mat, unsigned int i, unsigned int j, unsigned int n) {

    double *aux = malloc(n*sizeof(double));

    if(!aux) {
        perror("Sem memória");
        return -1;
    }
    memcpy(aux, &mat[n*i], n*sizeof(double));
    memcpy(&mat[n*i], &mat[n*j], n*sizeof(double));
    memcpy(&mat[n*j], aux, n*sizeof(double));

    free(aux);

    return 0;
}

/*!
  \brief Triangulariza a matriz SL->A de norma n
  \note separa SL->A em L e U

  \param SL o sistema linear
  \param mat a matriz a ser triangularizada
  \param B termos independentes
  \return 0 se sucesso e -1 em caso de falha
*/
int SL_triangulariza(t_sist *SL, double *mat, double *B) {
  
    if (!SL->L) {
        SL->L = SL_alocaMatrix(SL->n, SL->n);
        if (!SL->L) return -1;
    }
    memset(SL->L, 0, SL->n * SL->n * sizeof(double));

    if (!SL->U) {
        SL->U = SL_alocaMatrix(SL->n, SL->n);
        if (!SL->U) return -1;
    }
    memcpy(SL->U, mat, SL->n * SL->n * sizeof(double));

    // Transforma a matriz em uma triangular com pivoteamento parcial
    int pivo;
    for (int i=0; i<SL->n; i++) 
    {
        pivo = maxValue(SL->U,SL->n,i);
        if (pivo != i) {
            trocaElemento(B+i, B+pivo); // troca termo independente
            if (trocaLinha(SL->U, i, pivo, SL->n)) return -1;
            if (trocaLinha(SL->L, i, pivo, SL->n)) return -1;
        }

        SL->L[SL->n*i+i] = 1.0f;
        for (int j=i+1; j<SL->n; j++) {
            double m = SL->U[SL->n*j+i] / SL->U[SL->n*i+i];
            SL->U[SL->n*j+i] = 0.0f;
            SL->L[SL->n*j+i] = m;
            for (int k=i+1; k<SL->n; k++)
                SL->U[SL->n*j+k] -= SL->U[SL->n*i+k] * m;
        }
    }

    return 0;
}

/*!
  \brief Triangulariza a matriz SL->A de norma n
  \note separa SL->A em L e U

  \param SL o sistema linear
  \param mat a matriz a ser triangularizada
  \param B termos independentes
  \return 0 se sucesso e -1 em caso de falha
*/
int SL_triangulariza_otimiz(t_sist *SL, double *mat, double *B) {
  
    if (!SL->L) {
        SL->L = SL_alocaMatrix((SL->n * (SL->n+1))/2, 1);
        if (!SL->L) return -1;
    }
    memset(SL->L, 0, SL->n * SL->n * sizeof(double));

    if (!SL->U) {
        SL->U = SL_alocaMatrix((SL->n * (SL->n+1))/2, 1);
        if (!SL->U) return -1;
    }

    double *copia = SL_alocaMatrix(SL->n, SL->n);
    if (!copia) return -1;

    memcpy(copia, mat, SL->n * SL->n * sizeof(double));
    
    // Transforma a matriz em uma triangular com pivoteamento parcial
    int pivo;
    for (int i=0; i<SL->n; i++) 
    {
        pivo = maxValue(SL->U,SL->n,i);
        if (pivo != i) {
            trocaElemento(B+i, B+pivo); // troca termo independente
            trocaLinha(SL->U, i, pivo, SL->n);
            trocaLinha(SL->L, i, pivo, SL->n);
        }

        SL->L[SL->n*i+i] = 1.0f;
        for (int j=i+1; j<SL->n; j++) {
            double m = SL->U[SL->n*j+i] / SL->U[SL->n*i+i];
            SL->U[SL->n*j+i] = 0.0f;
            SL->L[SL->n*j+i] = m;
            for (int k=i+1; k<SL->n; k++)
                SL->U[SL->n*j+k] -= SL->U[SL->n*i+k] * m;
        }
    }

    free(copia);

    return 0;
}
