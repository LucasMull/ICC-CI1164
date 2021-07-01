/**
 * Luan Machado Bernardt | GRR20190363
 * Lucas Müller          | GRR20197160
 */

typedef struct {
    unsigned int n;
    float **A;
#if 0
    float *L, *U;
#endif
    float **L, **U;
} t_matrix;


t_matrix *alocaStruct(unsigned int n);
t_matrix *readMatrix();
void printMatrix(float **matrix, int n);
