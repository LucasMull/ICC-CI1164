// Luan Machado Bernardt | GRR20190363

typedef struct {
    unsigned int n;
    float **A;
    float *L, *U;
} t_matrix;

t_matrix *alocaStruct(unsigned int n);
t_matrix *readMatrix();
void printMatrix(float **matrix, int n);
