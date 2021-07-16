DESCRIÇÃO GERAL
    TRABALHO 1
EQUIPE
    Luan Machado Bernardt | GRR20190363
    Lucas Müller          | GRR20197160

NOME
    float** alocaMatrix(unsigned int n);
DESCRIÇÃO
    A função alocaMatrix() faz uma alocação única de tamanho 
    n*sizeof(float*) + n*n*sizeof(float). Esse design foi escolhido
    pois além de otimizar a localidade dos elementos em um único grande
    bloco de memória, também é menos suscetível a erros, uma vez que só
    requer um único free() para liberar memória.
RETORNO
    Retorna uma matriz do tipo float de tamanho nxn, ou NULL se houve
    falha.

NOME
    t_matrix *alocaStruct(unsigned int n);
DESCRIÇÃO
    Aloca na memória uma struct t_matrix contendo os seguintes campos :
        - mat->A :
            A matriz original de norma n.
        - mat->Inv :
            Espaço para a inversa da matriz, obtida em geraInversa().
        - mat->L e mat->U :
            Espaço para a matriz L(lower) e U(upper), obtidas em
            triangularizaMatrix().
RETORNO
    Uma struct t_matrix de norma nxn, ou NULL se houve falha.

NOME
    void limpaStruct(t_matrix *Mat);
DESCRIÇÃO
    Liberação da memória alocada por alocaStruct().

NOME
    float normaL2Residuo(t_matrix *Mat, float *matId, unsigned int col);
DESCRIÇÃO
    É calculada a deviância entre a matriz identidade resultante 
    do produto da matriz original com a inversa obtida por 
    geraInversa(), com a a da matriz identidade esperada.
RETORNO
    Valor do tipo float contendo a norma L2 do residuo.

NOME
    t_matrix *readMatrix();
DESCRIÇÃO
    Efetua leitura da matriz pela entrada padrão (stdin), retorna
    imediatamente após completar a primeira matriz, chamadas
    consecutivas devem ser feitas para conseguir qualquer matriz
    sequêncial.
RETORNO
    Um ponteiro para t_matriz contendo a matriz obtida nessa interação

NOME
    void printMatrix(FILE *f_out, float **matrix, int n);
DESCRIÇÃO
    Imprime a matriz nxn fornecida na saída padrão (stdout).

NOME
    static float determinanteU (float **U, unsigned int n);
DESCRIÇÃO
    Calcula a determinante da matriz já triangularizada,
    multiplicando os termos da diagonal principal de U.
RETORNO
    Determinante (Se != 0, quer dizer que a matriz original é inversível).

NOME
    int triangularizaMatrix(t_matrix *Mat, int pivotP, double *tTotal);
DESCRIÇÃO
    Efetua a triangularização de Mat->A e obtêm os campos Mat->L e 
    Mat->U, que serão utilizados para calcular a inversa em
    geraInversa().
RETORNO
    0 se sucesso e -1 em caso de falha.

NOME
    float **geraIdentidade(unsigned int n);
DESCRIÇÃO
    Gera a matriz identidade de norma n.
RETORNO
    A matriz identidade de norma n, ou NULL em caso de falha.

NOME
    void geraInversa(t_matrix *Mat, float **matId, double *timeLy, double *timeUx);
DESCRIÇÃO
    Gera a inversa da matriz original, apontada por Mat->A, a partir
    dos valores Mat->L e Mat->U encontrados em triangularizaMatrix().
