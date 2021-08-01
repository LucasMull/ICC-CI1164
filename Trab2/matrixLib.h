/**
 * Luan Machado Bernardt | GRR20190363
 * Lucas Müller          | GRR20197160
 */

#ifndef __MATRIXLIB__
#define __MATRIXLIB__

typedef struct {
    unsigned int n;
    float *A, *Inv, *Id;
    float *L, *U;
} t_sist;


float* SL_alocaMatrix(unsigned int n);
void SL_printMatrix(FILE *f_out, float *matrix, int n);

t_sist *SL_aloca(unsigned int n);
void SL_libera(t_sist *SL);
t_sist *SL_leitura();

float SL_normaL2Residuo(t_sist *SL, float *matId, unsigned int col);
int SL_triangulariza(t_sist *SL, int pivotP, double *tTotal);
float *SL_geraIdentidade(unsigned int n);
void SL_geraInversa(t_sist *SL, double *timeLy, double *timeUx);

#endif // __MATRIXLIB__
