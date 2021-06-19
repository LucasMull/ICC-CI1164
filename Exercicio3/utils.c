#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "utils.h"

static double timestamp(void) 
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}

float gaussSeidel(Edo *edoeq, float *Y, double *tempo, size_t n)
{
  int k, i;
  float h, xi, bi, yi, d, di, ds;
  h = (edoeq->b - edoeq->a)/(n+1);

  *tempo = timestamp();
  // Largura do passo da malha
  for (k=0; k < 50; ++k) { // 23 FLOP por iteração do método
    for (i=0; i < n; ++i) { // Para cada equação do SL
      xi = edoeq->a + (i+1)*h; // valor xi da malha:        2 FLOP
      bi = h*h * edoeq->r(xi); // termo independente:       3 FLOP
      di = 1 - h * edoeq->p(xi)/2.0; // diagonal inferior:  3 FLOP
      d = -2 + h*h * edoeq->q(xi); // diagonal principal:   3 FLOP
      ds = 1 + h * edoeq->p(xi)/2.0; // diagonal superior:  3 FLOP
      // 8 FLOP (maximo) ; 4 FLOP (mínimo)
      if (i == 0)
        bi -= ds*Y[i+1] + edoeq->ya * (1 - h*edoeq->p(edoeq->a+h)/2.0);
      else if (i == n-1) 
        bi -= di*Y[i-1] + edoeq->yb * (1 + h*edoeq->p(edoeq->b-h)/2.0);
      else
        bi -= ds*Y[i+1] + di*Y[i-1] ;
      Y[i] = bi / d; // Calcula incógnita: 1 FLOP
    }
  }
  *tempo = timestamp() - *tempo;
  return h;
}

SL_Tridiag *alocaSL(size_t n)
{
  SL_Tridiag *newSL = malloc(sizeof *newSL + 4*n*sizeof(float));
  if (NULL == newSL) return NULL;

  memset(newSL->mem, 0, 4*n*sizeof(float));

  newSL->D = newSL->mem;
  newSL->Di = newSL->mem + n;
  newSL->Ds = newSL->mem + 2*n;
  newSL->B = newSL->mem + 3*n;

  return newSL;
}

void geraTridiagonal(Edo *edoeq, SL_Tridiag *SL, size_t n)
{
  float h = (edoeq->b - edoeq->a) / (n+1.0f);

  float xi;
  for (int i=0; i < n; ++i) {
    xi = edoeq->a + (i+1)*h;                // ponto da malha
    SL->Di[i] = 1 - h * edoeq->p(xi)/2.0f;  // diagonal inferior
    SL->D[i] = -2 + h*h * edoeq->p(xi);     // diagonal principal
    SL->Ds[i] = 1 + h * edoeq->p(xi)/2.0f;  // diagonal superior
    SL->B[i] = h*h * edoeq->r(xi);          // termo independente
  }
  // Condições de contorno subtraídas do 1o e último termos independentes
  SL->B[0] -= edoeq->ya * (1 - h*edoeq->p(edoeq->a+h)/2.0);
  SL->B[n-1] -= edoeq->yb * (1 + h*edoeq->p(edoeq->b-h)/2.0);
}

float normaL2Residuo(SL_Tridiag *SL, float *Y, size_t n)
{
  float res[n];
  memset(res, 0, n*sizeof(float));

  float soma=0.0f;
  for (size_t i=0; i < n; ++i) {
    for (size_t j=0; j < n; ++j) {
      if (j == i)
        res[i] = fmaf(SL->D[i], Y[j], res[i]);
      else if (j == i + 1)
        res[i] = fmaf(SL->Ds[i], Y[j], res[i]);
      else if (j == i - 1)
        res[i] = fmaf(SL->Di[i], Y[j], res[i]);
    }
    res[i] = SL->B[i] - res[i];
    soma += powf(res[i], 2.0f);
  }
  return sqrtf(soma);
}


void prnVetor(float *v, size_t n)
{
  for (size_t i=0; i < n; ++i)
    printf("%1.7g ", v[i]);
  putchar('\n');
}

void prnSistLinear(SL_Tridiag *SL, size_t n)
{
  prnVetor(SL->Ds, n);
  prnVetor(SL->D, n);
  prnVetor(SL->Di, n);
  prnVetor(SL->B, n);
}
