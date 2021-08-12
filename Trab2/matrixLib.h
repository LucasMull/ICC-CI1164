/**
 * Luan Machado Bernardt | GRR20190363
 * Lucas Müller          | GRR20197160
 */

#ifndef __MATRIXLIB__
#define __MATRIXLIB__

typedef struct {
    unsigned int n, m;
    double *A, *Int, *Ajc;
    double *L, *U;
    double *x;
} t_sist;


double* SL_alocaMatrix(unsigned int n, unsigned int m);
void SL_printMatrix(FILE *f_out, double *matrix, unsigned int n, unsigned int m);

t_sist *SL_aloca(unsigned int n, unsigned int m);
void SL_libera(t_sist *SL);
t_sist *SL_leitura();

int SL_interpolacao(t_sist *SL, unsigned int row, double *B);
int SL_ajusteDeCurvas(t_sist *SL, unsigned int row, double *B, double *lookup);
int SL_triangulariza(t_sist *SL, double *mat, double *B);
int SL_triangulariza_otimiz(t_sist *SL, double *mat, double *B);
void SL_substituicao(t_sist *SL, double *B, double *pol);

#endif // __MATRIXLIB__
