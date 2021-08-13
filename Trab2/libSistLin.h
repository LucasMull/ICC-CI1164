/**
 * Luan Machado Bernardt | GRR20190363
 * Lucas Müller          | GRR20197160
 */

#ifndef __LIBSISTLIN__
#define __LIBSISTLIN__

typedef struct {
    unsigned int n, m;
    double *A;
    double *L, *U;
    int *vetTroca;
    union { double *x, *B; };
} t_sist;


double* SL_alocaMatrix(unsigned int n, unsigned int m);
void SL_printMatrix(FILE *f_out, double *matrix, unsigned int n, unsigned int m);

t_sist *SL_aloca(unsigned int n, unsigned int m);
void SL_libera(t_sist *SL);
t_sist *SL_leitura();

int SL_interpolacao(t_sist *SL, t_sist *Int, unsigned int row);
int SL_ajusteDeCurvas(t_sist *SL, t_sist *Ajc, unsigned int row, double *lookup);
int SL_triangulariza(t_sist *SL);
int SL_triangulariza_otimiz(t_sist *SL);
void SL_substituicao(t_sist *SL, double *pol);

#endif // __LIBSISTLIN__
