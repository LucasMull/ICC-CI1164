/**
 * Luan Machado Bernardt | GRR20190363
 * Lucas Müller          | GRR20197160
 */

#ifndef __MATRIXLIB__
#define __MATRIXLIB__

typedef struct {
    unsigned int n, m;
    float *A, *Int;
    float *L, *U;
    float *x;
} t_sist;


float* SL_alocaMatrix(unsigned int n, unsigned int m);
void SL_printMatrix(FILE *f_out, float *matrix, unsigned int n, unsigned int m);

t_sist *SL_aloca(unsigned int n, unsigned int m);
void SL_libera(t_sist *SL);
t_sist *SL_leitura();

float *SL_interpolacao(t_sist *SL, unsigned int row);
float *SL_ajusteDeCurvas(t_sist *SL, unsigned int row);
int SL_triangulariza(t_sist *SL, float *B, double *tTotal);

#endif // __MATRIXLIB__
