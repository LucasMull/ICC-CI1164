// Luan Machado Bernardt | GRR20190363

typedef struct {
    unsigned int n;
    float **A;
} t_matrix;

t_matrix *alocaMatrix(unsigned int n);
t_matrix *readMatrix();
void printMatrix(float **matrix, int n);
