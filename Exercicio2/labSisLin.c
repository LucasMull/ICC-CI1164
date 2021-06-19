#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

void prnSolucao(SistLinear_t *SL, real_t *x)
{
  real_t res[SL->n]; // vetor residuo
  real_t norma = normaL2Residuo(SL, x, res);

  printf("  --> X: ");
  prnVetor(x, SL->n);
  printf("\n  --> Norma L2 do residuo: %1.7g\n", norma);
  if (norma > 5.0) {
    double tTotal;
    int iter = refinamento(SL, x, &tTotal);
    switch (iter) {
    default:
        printf("\n===> Refinamento: %lf ms --> %d iterações\n", tTotal, iter);
        break;
    case -1:
        printf("\n===> Refinamento (Não Converge): %lf ms\n", tTotal);
        break;
    case -2:
        printf("\n===> Refinamento (Sem Solução): %lf ms\n", tTotal);
        break;
    }
    printf("  --> X: ");
    prnVetor(x, SL->n);
    printf("\n  --> Norma L2 do residuo: %1.7g\n", normaL2Residuo(SL, x, res));
  }
  putchar('\n');
}

int main()
{
  SistLinear_t *SL;
  int iter;
  double tTotal;
  unsigned int k=1;

  while (!feof(stdin)) 
  {
    SL = lerSistLinear();
    if (!SL) return -1;

    printf("***** Sistema %u --> n = %u, erro: %g\n", k, SL->n, SL->erro);

    real_t x[SL->n]; // vetor solução

    if (0 == eliminacaoGauss(SL, x, &tTotal))
      printf("===> Eliminação Gauss: %lf ms\n", tTotal);
    else
      printf("===> Eliminação Gauss (Sem Solução): %lf ms\n", tTotal);
    prnSolucao(SL, x);

    switch (iter = gaussJacobi(SL, x, &tTotal)) {
    default:
        printf("===> Jacobi: %lf ms --> %d iterações\n", tTotal, iter);
        break;
    case -1:
        printf("===> Jacobi (Não Converge): %lf ms\n", tTotal);
        break;
    case -2:
        printf("===> Jacobi (Sem Solução): %lf ms\n", tTotal);
        break;
    }
    prnSolucao(SL, x);

    switch (iter = gaussSeidel(SL, x, &tTotal)) {
    default:
        printf("===> Gauss-Seidel: %lf ms --> %d iterações\n", tTotal, iter);
        break;
    case -1:
        printf("===> Gauss-Seidel (Não Converge): %lf ms\n", tTotal);
        break;
    case -2:
        printf("===> Gauss-Seidel (Sem Solução): %lf ms\n", tTotal);
        break;
    }
    prnSolucao(SL, x);

    // libera sistema atual da memória
    liberaSistLinear(SL);

    ++k;
  }

  return 0;
}
