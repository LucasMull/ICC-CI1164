#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "SistemasLineares.h"


/*!
  \brief Retorna o indice da linha pivoLn que possui o maior valor abs
  \param sistema linear SL
  \param linha pivoLn
  \param coluna pivoCol
  \return Índice do maior valor
  */
static unsigned int maxArgIdx(SistLinear_t *SL, unsigned int pivoLn, unsigned int pivoCol) 
{
  unsigned int idxMax = pivoLn;
  for (unsigned int i = pivoLn+1; i < SL->n; ++i)
    if (fabs(SL->A[i][pivoCol] > fabs(SL->A[idxMax][pivoCol])))
      idxMax = i;
  return idxMax;
}

/*!
  \brief Troca linhas pivoLn e idxMax do Sistema Linear
  \param sistema linear SL
  \param linha pivoLn
  \param linha idxMax
 */
static void lnTroca(SistLinear_t *SL, unsigned int pivoLn, unsigned int idxMax)
{
  real_t *auxA = SL->A[pivoLn], auxb = SL->b[pivoLn];
  SL->A[pivoLn] = SL->A[idxMax];
  SL->A[idxMax] = auxA;
  SL->b[pivoLn] = SL->b[idxMax];
  SL->b[idxMax] = auxb;
}

/*!
  \brief Retorna uma cópia do sistema linear
  \param sistema linear SL
  \return Cópia do sistema linear. NULL se houve erro (alocação)
 */
static SistLinear_t* dupSistLinear(SistLinear_t *SL)
{
  SistLinear_t *cpSL = alocaSistLinear(SL->n);
  for (unsigned int i=0; i < SL->n; ++i)
    memcpy(cpSL->A[i], SL->A[i], SL->n*sizeof(real_t));
  memcpy(cpSL->b, SL->b, SL->n*sizeof(real_t));
  cpSL->erro = SL->erro;
  return cpSL;
}

/*!
  \brief Verifica e retorna a maior diferença entre o vetor solução da iteração passada e o da atual
  \param vetor solução atual
  \param vetor solução anterior
  \return Maior diferença encontrada comparando vetor atual com anterior
 */
static real_t maiorDif(real_t *atual, real_t *anterior, unsigned int n)
{
  real_t max=0.0f;
  for (unsigned int i=0; i < n; ++i)
    max = fmax(fabs(anterior[i] - atual[i]), max);
  return max;
}

/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
  \param res Valor do resíduo

  \return Norma L2 do resíduo.
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x, real_t *res)
{
  memset(res, 0, SL->n*sizeof(real_t)); // seta vetor res em 0

  real_t soma=0.0f;
  for (unsigned int i=0; i < SL->n; ++i) {
    for (unsigned int j=0; j < SL->n; ++j)
      res[i] = fmaf(SL->A[i][j], x[j], res[i]);
    res[i] = SL->b[i] - res[i];
    soma += powf(res[i], 2.0f);
  }
  return sqrtf(soma);
}

/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param tTotal tempo gasto pelo método

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss(SistLinear_t *SL, real_t *x, double *tTotal)
{
  SistLinear_t *cpSL = dupSistLinear(SL);
  unsigned int pivoLn=0, pivoCol=0;
  unsigned int idxMax=0;

  *tTotal = timestamp();

  while (pivoLn < cpSL->n && pivoCol < cpSL->n) {
    idxMax = maxArgIdx(cpSL, pivoLn, pivoCol);
    if (0.0f == cpSL->A[idxMax][pivoCol]) {
      ++pivoCol;
      continue;
    }

    if (idxMax != pivoLn)
      lnTroca(cpSL, pivoLn, idxMax);

    for (unsigned int i = pivoLn+1; i < cpSL->n; ++i) {
      const real_t coef = cpSL->A[i][pivoCol] / cpSL->A[pivoLn][pivoCol];
      cpSL->A[i][pivoCol] = 0.0f;
      for (unsigned int j = pivoCol+1; j < cpSL->n; ++j)
        cpSL->A[i][j] = fmaf(-coef, cpSL->A[pivoLn][j], cpSL->A[i][j]);
      cpSL->b[i] = fmaf(-coef, cpSL->b[pivoLn], cpSL->b[i]);
    }
    ++pivoLn;
    ++pivoCol;
  }

  // retrosubstituição
  for (int i = cpSL->n-1; i >= 0; --i) {
    for (int j=i+1; j < cpSL->n; ++j)
      cpSL->b[i] -= cpSL->A[i][j] * cpSL->b[j];
    cpSL->b[i] /= cpSL->A[i][i];
  }

  *tTotal = timestamp() - *tTotal;

  memcpy(x, cpSL->b, SL->n*sizeof(real_t));

  liberaSistLinear(cpSL);

  return 0;
}

/*!
  \brief Método de Jacobi

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
*/
int gaussJacobi(SistLinear_t *SL, real_t *x, double *tTotal)
{
  real_t anterior[SL->n], atual[SL->n], res[SL->n];
  memset(anterior, 0, SL->n*sizeof(real_t));

  *tTotal = timestamp();

  unsigned int iter=0;
  while (iter < MAXIT) {
    ++iter;

    for (unsigned int i=0; i < SL->n; ++i) {
      real_t soma=0.0f;
      for (unsigned int j=0; j < SL->n; ++j) {
        if (j != i) 
          soma += SL->A[i][j] * anterior[j] / SL->A[i][i];
        atual[i] = SL->b[i] / SL->A[i][i] - soma;
      }
    }

    if (SL->erro >= maiorDif(atual, anterior, SL->n))
      break;

    memcpy(anterior, atual, SL->n*sizeof(real_t));
  }

  *tTotal = timestamp() - *tTotal;

  memcpy(x, atual, SL->n*sizeof(real_t));

  return iter;
}

/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int gaussSeidel(SistLinear_t *SL, real_t *x, double *tTotal)
{
  real_t anterior[SL->n], atual[SL->n], res[SL->n];
  memset(anterior, 0, SL->n*sizeof(real_t));

  *tTotal = timestamp();

  unsigned int iter=0;
  while (iter < MAXIT) 
  {
    ++iter;

    for (unsigned int i=0; i < SL->n; ++i) {
      real_t soma=0.0f;
      for (unsigned int j=0; j < i; ++j)
        soma = fmaf(SL->A[i][j], atual[j], soma);
      for (unsigned int j = i+1; j < SL->n; ++j)
        soma = fmaf(SL->A[i][j], anterior[j], soma);
      atual[i] = (SL->b[i] / SL->A[i][i]) - (soma / SL->A[i][i]);
    }

    if (SL->erro >= maiorDif(atual, anterior, SL->n))
      break;

    memcpy(anterior, atual, SL->n*sizeof(real_t));
  }

  *tTotal = timestamp() - *tTotal;

  memcpy(x, atual, SL->n*sizeof(real_t));

  return iter;
}


/*!
  \brief Método de Refinamento

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial para início do refinamento
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int refinamento(SistLinear_t *SL, real_t *x, double *tTotal)
{
  SistLinear_t *cpSL = dupSistLinear(SL);
  real_t norma;
  real_t res[cpSL->n], w[cpSL->n];

  double auxtTotal = timestamp();

  int iter=0;
  while ((iter < MAXIT) && (5.0 < (norma = normaL2Residuo(SL, x, res))) ) 
  {
    memcpy(cpSL->b, res, cpSL->n*sizeof(real_t));
    eliminacaoGauss(cpSL, w, tTotal);
    for (unsigned int i=0; i < cpSL->n; ++i)
      x[i] += w[i];

    if (SL->erro >= maiorDif(x, w, SL->n))
      break;

    memset(w, 0, cpSL->n*sizeof(real_t));

    ++iter;
  }

  *tTotal = timestamp() - auxtTotal;

  liberaSistLinear(cpSL);

  return iter;
}

/*!
  \brief Alocaçao de memória 

  \param n tamanho do SL

  \return ponteiro para SL. NULL se houve erro de alocação
  */
SistLinear_t* alocaSistLinear(unsigned int n) 
{
  if (!n) return NULL;

  SistLinear_t *novoSL = malloc(sizeof(SistLinear_t));
  if (!novoSL) return NULL;

  /* efetua alocação de matriz em 1D para facilitar limpeza */
  novoSL->A = malloc(n*sizeof(real_t*) + n*n*sizeof(real_t));
  if (!novoSL->A) {
    free(novoSL);
    return NULL;
  }
  /* inicializa cada ponteiro para seu bloco de memória
   *        consecutivo alocado */
  real_t *addr = (real_t*)(novoSL->A + n);
  for (unsigned int i=0; i < n; ++i) {
    novoSL->A[i] = addr;
    addr += n;
  }
  /* inicializa termos independentes */
  novoSL->b = malloc(n*sizeof(real_t));
  if (!novoSL->b) {
    free(novoSL);
    free(novoSL->A);
    return NULL;
  }

  novoSL->n = n;

  return novoSL;
}

/*!
  \brief Liberaçao de memória 

  \param sistema linear SL
  */
void liberaSistLinear(SistLinear_t *SL)
{
  if (!SL) return;
  /* por ter sido feita em uma única alocação só é necessário um
   *        único free para SL->A */
  if (SL->A)
    free(SL->A);
  if (SL->b)
    free(SL->b);
  free(SL);
}

/*!
  \brief Leitura de SL a partir de Entrada padrão (stdin).

  \return sistema linear SL. NULL se houve erro (leitura ou alocação)
  */
SistLinear_t *lerSistLinear()
{
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

  // extrai a linha contendo o critério de parada
  if (!fgets(ln, sizeof(ln), stdin)) {
    perror("Falha de leitura");
    return NULL;
  }
  real_t erro = strtof(ln, NULL);
  if (!erro) {
    fputs("Não foi possível obter critério de parada\n", stderr);
    return NULL;
  }

  // realiza alocação do Sistema Linear
  SistLinear_t *novoSL = alocaSistLinear(n);
  if (!novoSL) {
    fputs("Não foi possível alocar 'novo_SL'\n", stderr);
    return NULL;
  }
  novoSL->erro = erro;

  // extrai as linhas contendo os elementos da matriz
  char *valorAtual, *valorProx;
  for (unsigned int i=0; i < novoSL->n; ++i) {
    valorAtual = ln;
    if (!fgets(ln, sizeof(ln), stdin)) {
      perror("Falha de leitura");
      return NULL;
    }

    for (unsigned int j=0; j < novoSL->n; ++j) {
      novoSL->A[i][j] = strtof(valorAtual, &valorProx);
      if (!valorProx) {
        fputs("Não foi possível obter o próximo valor\n", stderr);
        return NULL;
      }
      valorAtual = valorProx;
    }
  }

  // extrai a linha contendo os termos independentes
  valorAtual = ln;
  if (!fgets(ln, sizeof(ln), stdin)) {
    perror("Falha de leitura");
    return NULL;
  }

  for (unsigned int i=0; i < novoSL->n; ++i) {
    novoSL->b[i] = strtof(valorAtual, &valorProx);
    if (!valorProx) {
      fputs("Não foi possível obter o próximo valor\n", stderr);
      return NULL;
    }
    valorAtual = valorProx;
  }

  // "consome" próxima linha (vazia ou parada em EOF)
  fgets(ln, sizeof(ln), stdin);

  return novoSL;
}


// Exibe SL na saída padrão
void prnSistLinear(SistLinear_t *SL)
{
  for (unsigned int i=0; i < SL->n; ++i) {
    prnVetor(SL->A[i], SL->n);
    printf("= %1.7g\n", SL->b[i]);
  }
  putchar('\n');
}

// Exibe um vetor na saída padrão
void prnVetor(real_t *v, unsigned int n)
{
  for (unsigned int i=0; i < n; ++i)
    printf("%1.7g ", v[i]);
  fflush(stdout);
}

