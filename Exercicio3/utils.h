#ifndef __UTILS_H__
#define __UTILS_H__

// Matriz tridiagonal
typedef struct {
  // diagonal principal, inferior, superior e termo independente
  float *D, *Di, *Ds, *B;
  float mem[];
} SL_Tridiag;

// Equação Diferencial Ordinária para MDF
typedef struct {
  float a, b;   // intervalo
  float ya, yb; // condições contorno
  float (*p)(float), (*q)(float), (*r)(float);
} Edo;

// Equação Diferencial Parcial para MDF
typedef struct {
  float a, b;         // intervalo
  float (*u1)(float), // (0, y)
        (*u2)(float), // (a, y)
        (*u3)(float), // (x, 0) 
        (*u4)(float); // (x, b)
  float (*fn)(float, float);
} Edp;

float gaussSeidel(Edo *edoeq, float *Y, double *tempo, size_t n);

SL_Tridiag *alocaSL(size_t n);
void geraTridiagonal(Edo *edoeq, SL_Tridiag *SL, size_t n);
float normaL2Residuo(SL_Tridiag *SL, float *Y, size_t n);
void prnSistLinear(SL_Tridiag *SL, size_t n);
void prnVetor(float *v, size_t n);

#endif /* __UTILS_H__ */
