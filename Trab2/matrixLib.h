/**
 * Luan Machado Bernardt | GRR20190363
 * Lucas Müller          | GRR20197160
 */

#ifndef __MATRIXLIB__
#define __MATRIXLIB__

typedef struct {
    unsigned int n;
    float **A;
    float **Inv;
    float **Id;
#if 0
    float *L, *U;
#endif
    float **L, **U;
} t_matrix;


float** alocaMatrix(unsigned int n);
t_matrix *alocaStruct(unsigned int n);
void limpaStruct(t_matrix *Mat);
float normaL2Residuo(t_matrix *Mat, float *matId, unsigned int col);
t_matrix *readMatrix();
void printMatrix(FILE *f_out, float **matrix, int n);
int triangularizaMatrix(t_matrix *Mat, int pivotP, double *tTotal);
float **geraIdentidade(unsigned int n);
void geraInversa(t_matrix *Mat, double *timeLy, double *timeUx);

#endif // __MATRIXLIB__
