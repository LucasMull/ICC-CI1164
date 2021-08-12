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
 * \param pol vetor do polinomio resultante 1xn previamente alocado
 */
void SL_substituicao(t_sist *SL, double *B, double *pol) {

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
      free(newSL->A);
      free(newSL);
      return NULL;
  }
  newSL->L = SL_alocaMatrix(n, n);
  if (!newSL->L) {
    free(newSL->x);
    free(newSL->A);
    free(newSL);
    return NULL;
  }
  newSL->U = SL_alocaMatrix(n, n);
  if (!newSL->U) {
    free(newSL->L);
    free(newSL->x);
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
  free(SL->L);
  free(SL->U);
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
    if (EOF == scanf("%d", &n)) return NULL;
    if (EOF == scanf("%d", &m)) return NULL;
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

    for (unsigned int i=0; i < newSL->n; ++i) {
        if (EOF == scanf("%lf", &newSL->x[i])) {
            perror("Falha de leitura");
            return NULL;
        }
    }

    // extrai as linhas contendo os elementos da matriz
    for (unsigned int i=0; i < newSL->m; ++i) {
        for (unsigned int j=0; j < newSL->n; ++j) {
            if (EOF == scanf("%lf", &newSL->A[n*i+j])) {
                perror("Falha de leitura");
                return NULL;
            }
        }
    }

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
 * \return retorna 0 para sucesso e -1 para falha
 *
 * \note A otimização escolhida foi unroll e jam
 */
int SL_interpolacao(t_sist *SL, unsigned int row, double *B) {

  if (!SL->Int) {
      SL->Int = SL_alocaMatrix(SL->n, SL->n);
      if (!SL->Int) return -1;
      
      for (int i=0; i<(SL->n - (SL->n % 8)); i += 8) {
          // NOTA OTIMIZAÇÃO: como qualquer valor elevado a 0 é possível evitar chamadas desnecessárias a pot()
          SL->Int[SL->n*i] = SL->Int[SL->n*(i+1)] =     \
          SL->Int[SL->n*(i+2)] = SL->Int[SL->n*(i+3)] = \
          SL->Int[SL->n*(i+4)] = SL->Int[SL->n*(i+5)] = \
          SL->Int[SL->n*(i+6)] = SL->Int[SL->n*(i+7)] = 1.0;
          
          // NOTA OTIMIZAÇÃO: como qualquer valor elevado a 1 é ele mesmo é possível evitar chamadas desnecessárias a pot()
          SL->Int[SL->n*i+1] = SL->x[i];
          SL->Int[SL->n*(i+1)+1] = SL->x[i+1];
          SL->Int[SL->n*(i+2)+1] = SL->x[i+2];
          SL->Int[SL->n*(i+3)+1] = SL->x[i+3];
          SL->Int[SL->n*(i+4)+1] = SL->x[i+4];
          SL->Int[SL->n*(i+5)+1] = SL->x[i+5];
          SL->Int[SL->n*(i+6)+1] = SL->x[i+6];
          SL->Int[SL->n*(i+7)+1] = SL->x[i+7];

          // NOTA OTIMIZAÇÃO: é possível evitar chamadas à função pot() se multiplicar o valor da coluna anterior ao da coluna 1
          for (int j=2; j < SL->n; ++j) {
              SL->Int[SL->n*i+j] = SL->Int[SL->n*i+1] * SL->Int[SL->n*i+(j-1)];
              SL->Int[SL->n*(i+1)+j] = SL->Int[SL->n*(i+1)+1] * SL->Int[SL->n*(i+1)+(j-1)];
              SL->Int[SL->n*(i+2)+j] = SL->Int[SL->n*(i+2)+1] * SL->Int[SL->n*(i+2)+(j-1)];
              SL->Int[SL->n*(i+3)+j] = SL->Int[SL->n*(i+3)+1] * SL->Int[SL->n*(i+3)+(j-1)];
              SL->Int[SL->n*(i+4)+j] = SL->Int[SL->n*(i+4)+1] * SL->Int[SL->n*(i+4)+(j-1)];
              SL->Int[SL->n*(i+5)+j] = SL->Int[SL->n*(i+5)+1] * SL->Int[SL->n*(i+5)+(j-1)];
              SL->Int[SL->n*(i+6)+j] = SL->Int[SL->n*(i+6)+1] * SL->Int[SL->n*(i+6)+(j-1)];
              SL->Int[SL->n*(i+7)+j] = SL->Int[SL->n*(i+7)+1] * SL->Int[SL->n*(i+7)+(j-1)];
          }
      }

      // Remainder loop
      for (int i=(SL->n - (SL->n % 8)); i < SL->n; ++i) {
          SL->Int[SL->n*i] = 1.0;
          SL->Int[SL->n*i+1] = SL->x[i];
          for (int j=2; j<SL->n; ++j) {
              SL->Int[SL->n*i+j] = SL->Int[SL->n*i+1] * SL->Int[SL->n*i+(j-1)];
          }
      }
  }

  // copia termos independentes para B
  memcpy(B, SL->A + SL->n*row, SL->n*sizeof(double));

  return 0;
}

/*!
 * \brief Realiza ajuste de curvas na matriz
 *
 * \param SL sistema linear
 * \param row a linha da matriz de entrada
 * \param B termos independentes
 * \return retorna 0 para sucesso e -1 para falha
 */
int SL_ajusteDeCurvas(t_sist *SL, unsigned int row, double *B) {

  if (!SL->Ajc) {
      SL->Ajc = SL_alocaMatrix(SL->n, SL->n);
      if (!SL->Ajc) return -1;
        
      // primeira linha da matriz (j: coluna, k: somatório)
      // zera os valores de B
      memset(B, 0, SL->n*sizeof(double));
      for (unsigned int j=0; j < SL->n; ++j)
          for (unsigned int k=0; k < SL->n; ++k)
              SL->Ajc[j] += pot(SL->x[k], j);

      // restante das linhas da matriz (i: linhas, j: coluna, k: somatório)
      for (unsigned int i=1; i < SL->n; ++i) {
          for (unsigned int j=0; j < SL->n-1; ++j) { // Loop merging
              SL->Ajc[SL->n*i+j] = SL->Ajc[SL->n*(i-1)+(j+1)];
              SL->Ajc[SL->n*i+(SL->n-1)] += pot(SL->x[j], i) * pot (SL->x[j], SL->n-1); // TODO: explicar otimização junta expoentes mesma base
          }
          SL->Ajc[SL->n*i+(SL->n-1)] += pot(SL->x[SL->n-1], i) * pot(SL->x[SL->n-1], SL->n-1); // Última soma
      }
  }

  // zera os valores de B
  memset(B, 0, SL->n*sizeof(double));
  for (unsigned int i=0; i<SL->n; ++i)
      for (unsigned int j=0; j<SL->n; ++j)
          B[i] += SL->A[SL->n*row+j] * pot(SL->x[j], i);

  return 0;
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
static void trocaLinha(double *mat, unsigned int i, unsigned int j, unsigned int n) {

#if 0 // Versão com realocação
    double aux = malloc(n*sizeof(double));
    if(!aux) {
        perror("Sem memória");
        return -1;
    }
    memcpy(aux, &mat[n*i], n*sizeof(double));
    memcpy(&mat[n*i], &mat[n*j], n*sizeof(double));
    memcpy(&mat[n*j], aux, n*sizeof(double));

    free(aux);
#else // Versão sem realocação
    double aux;
    double *a = &mat[n*i], *b = &mat[n*j];

    for (unsigned int k = 0; k<n; ++k) {
        aux = *(a+k);
        *(a+k) = *(b+k);
        *(b+k) = aux;
    }    
#endif
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
  
    memcpy(SL->U, mat, SL->n * SL->n * sizeof(double));

    // Transforma a matriz em uma triangular com pivoteamento parcial
    int pivo;
    double m, divi;

    for (int i=0; i<SL->n; i++) {
        pivo = maxValue(SL->U,SL->n,i);
        if (pivo != i) {
            trocaElemento(B+i, B+pivo); // troca termo independente
            trocaLinha(SL->U, i, pivo, SL->n);
            trocaLinha(SL->L, i, pivo, SL->n);
        }

        divi = SL->U[SL->n*i+i];
        SL->L[SL->n*i+i] = 1.0;
        for (int j=i+1; j<SL->n; j++) {
            m = SL->U[SL->n*j+i] / divi;
            SL->U[SL->n*j+i] = 0.0;
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
int SL_triangulariza(t_sist *SL, double *mat, double *B) {
  
    double *copia = SL_alocaMatrix(SL->n, SL->n);
    if (!copia) return -1;
    memcpy(copia, mat, SL->n * SL->n * sizeof(double));
    
    // Transforma a matriz em uma triangular com pivoteamento parcial
    int pivo;
    for (int i=0; i<SL->n; i++) 
    {
        pivo = maxValue(copia,SL->n,i);
        if (pivo != i) {
            trocaElemento(B+i, B+pivo); // troca termo independente
            trocaLinha(copia, i, pivo, SL->n);
            trocaLinha(SL->L, i, pivo, SL->n);
        }

        SL->L[SL->n*i+i] = 1.0;
        for (int j=i+1; j<SL->n; j++) {
            double m = copia[SL->n*j+i] / copia[SL->n*i+i];
            copia[SL->n*j+i] = 0.0;
            SL->L[SL->n*j+i] = m;
            for (int k=i+1; k<SL->n; k++)
                copia[SL->n*j+k] -= copia[SL->n*i+k] * m;
        }
    }

    if (SL->U) free(SL->U);
    SL->U = copia;

    return 0;
}
