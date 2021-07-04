/**
 * Luan Machado Bernardt | GRR20190363
 * Lucas Müller          | GRR20197160
 */

typedef struct {
    unsigned int n;
    float **A;
    float **Inv;
#if 0
    float *L, *U;
#endif
    float **L, **U;
} t_matrix;


t_matrix *alocaStruct(unsigned int n);
t_matrix *readMatrix();
void printMatrix(float **matrix, int n);
int triangularizaMatrix(t_matrix *Mat, int pivotP, double *tTotal);
float **geraIdentidade(unsigned int n);
int LyI(t_matrix *Mat, float **matId);