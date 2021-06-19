#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include "utils.h"

/**
 * 1. Para todos os casos, serão gerados resultados para n = 5 e 10 (malha x).
 *      Nos casos (b) e (d) também será considerado m = 3 (malha y).
 * 2. Para cada caso será exibida:
 *  2a. A matriz aumentada do SL resultante (matriz completa OU diagonais e termos independentes apenas),
 *  2b. A solução calculada com o método de Gauss-Seidel
 *  2c. A norma L2 do resíduo e o tempo gasto para a solução
 */
#define N5 5
#define N10 10

/**
 * Equação (a):
 *  y" = 6x - 0.5x², x ∈ (0,12)
 *  y(0) = 0 e y(12) = 0
 */

float eqA_p(float x) { (void)x;return 0.0f; }
float eqA_q(float x) { (void)x;return 0.0f; }
float eqA_r(float x) { return 6.0f*x - 0.5f*powf(x, 2.0f); }

/**
 * Equação (b):
 * Tₓₓ + Tᵧᵧ - T = sin²(x), (x,y) ∈ Ω = (0,L) x (0,W)
 * T(0,y) = 20, T(L,y) = 45, T(x,0) = 0, T(x,W) = 100
 */

/**
 * Equação (c):
 *  y" + y = 0, x ∈ (0,1)
 *  y(0) = 0 e y(1) = 1
 */

float eqC_p(float x) { (void)x;return 0.0f; }
float eqC_q(float x) { (void)x;return 1.0f; }
float eqC_r(float x) { (void)x;return 0.0f; }

/**
 * Equação (d):
 * uₓₓ + uᵧᵧ - T = -cos(x+y) -cos(x-y), (x,y) ∈ (0,π) x (0,π/2)
 * u(0,y) = cos(y), u(π,y) = -cos(y), u(x,0) = cos(x), u(x,π/2) = 0
 */

void prnSolucao(SL_Tridiag *SL, float *Y, double tempo, size_t n)
{
  puts("SL:");
  prnSistLinear(SL, n);
  printf("Y: ");
  prnVetor(Y, n);
  printf("Norma L2: %1.7e, Tempo: %1.7lf\n\n", \
      normaL2Residuo(SL, Y, n), tempo);
}

int main(void) 
{ 
  union { Edo o; Edp p; } edeq={0};
  float Y[N10]={0.0f};
  SL_Tridiag *SL = alocaSL(N10);
  assert(NULL != SL && "Sem memória");

  double tempo;
  int n=N5;
  do {
    float h; // distância malha

    /* Equação (a) */
    edeq.o = (Edo){
      .a = 0, .b = 12,
      .ya = 0, .yb = 0,
      .p = eqA_p, .q = eqA_q, .r = eqA_r 
    };
    h = gaussSeidel(&edeq.o, Y, &tempo, n);
    geraTridiagonal(&edeq.o, SL, n);
    printf("***** item (a): n = %u, H = %1.7f\n", n, h);
    prnSolucao(SL, Y, tempo, n);

    /* Equação (c) */
    edeq.o = (Edo){
      .a = 0, .b = 1,
      .ya = 0, .yb = 1,
      .p = eqC_p, .q = eqC_q, .r = eqC_r 
    };
    gaussSeidel(&edeq.o, Y, &tempo, n);
    geraTridiagonal(&edeq.o, SL, n);
    printf("***** item (c): n = %u, H = %1.7f\n", n, h);
    prnSolucao(SL, Y, tempo, n);

    n *= 2;
  } while (n <= N10);

  free(SL);

  return EXIT_SUCCESS;
}
