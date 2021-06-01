/* Lucas Müller GRR20197160 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> /* int32_t */
#include <math.h> /* nextafterf(), INFINITY */
#include <float.h> /* FLT_* */
#include <assert.h>

#define MAX_ULPS 3


/**
 * Retorna o menor valor float do vetor fornecido
 *
 * @param arr array de floats
 * @param size quantidade de elementos
 * @return menor valor float
 */
float faminf(float *arr, size_t size) 
{
  float min = arr[0];
  for (size_t i=1; i < size; ++i)
    min = fminf(min, arr[i]);
  return min;
}

/**
 * Retorna o maior valor float do vetor fornecido
 *
 * @param arr array de floats
 * @param size quantidade de elementos
 * @return maior valor float
 */
float famaxf(float *arr, size_t size)
{
  float max = arr[0];
  for (size_t i=1; i < size; ++i)
    max = fmaxf(max, arr[i]);
  return max;
}

/**
 * Verifica se dois float são aproximadamente iguais
 *
 * @param (a, b) intervalos
 * @param max_ulps distância máxima para ser considerado iguais
 * @return #1 iguais, #0 diferentes
 */
_Bool falmosteq(float a, float b, int max_ulps)
{
  /* Distância 0 para floats iguais */
  if (a == b) return 1;

  /* Distância máxima para NaN ou Inf */
  if (isnan(a) || isnan(b) || isinf(a) || isinf(b))
    return 0;

  union {
    float f;
    int32_t i;
  } ua = { .f = a }, ub = { .f = b };

  /* Não compara floats de sinais diferentes */
  if ((ua.i < 0) != (ub.i < 0)) return 0;

  return abs(ua.i - ub.i) < max_ulps;
}

/**
 * Checa se o intervalo é válido, caso não seja
 *        retorna do programa com -1
 *
 * @param (min, max) delimitação do intervalo
 */
void validate_interval(float min, float max) 
{
  if ((min > max)                           /* ex: [2, 1] */ 
      || (min < -FLT_MAX && max < -FLT_MAX) /* [-inf, -inf] */
      || (min > FLT_MAX && max > FLT_MAX)   /* [inf, inf] */
      || (isnan(min) || isnan(max)))        /* [NaN, NaN] */
  {
    fprintf(stderr, "Nao e um intervalo valido: [%1.8e, %1.8e]\n", min, max);
    exit(-1);
  }
}

int main(void) 
{
  /**
   * @ret_value valor de retorno da saída
   *        #0 SUCESSO, #-1 FALHA
   * @n_values número de valores reais conhecidos 
   * @n_ops operações que devem ser computadas
   */
  size_t n_values=0, n_ops=0; 
  scanf("%zu %zu%*[^\n]", &n_values, &n_ops);
  scanf("%*c"); // consome \n
  assert(n_values != 0 && "Informe a qtd de valores reais");
  assert(n_ops != 0 && "Informe a qtd de operacoes");

  /**
   * @idx percore indice do array @intvs
   * @recv_idx indice extraido da entrada do usuario, para checagem de erros
   * @intvs intervalos min/max extraídos a partir do valor de entrada
   */
  size_t idx=0, recv_idx=0;
  struct {
    float min, max;
  } *intvs = calloc(1, (n_values + n_ops) * sizeof *intvs);
  assert(intvs != NULL && "Sem memória");

  float tmp;
  for (size_t i=0; i < n_values; ++i) {
    scanf("x%zu %f%*[^\n]", &recv_idx, &tmp);
    scanf("%*c"); // consome \n
    assert(recv_idx == idx+1 && "Indice de valores fora de ordem");

    intvs[idx].min = nextafterf(tmp, -INFINITY);
    intvs[idx].max = nextafterf(tmp, INFINITY);
    validate_interval(intvs[idx].min, intvs[idx].max);
    printf("X%zu = [%1.8e, %1.8e]\n", recv_idx, intvs[idx].min, intvs[idx].max);
    
    ++idx;
  }

  /**
   * @op char que indica operação a ser efetuada
   * @(x, y) índice dos operandos @intvs a efetuar operação
   */
  char op=0;
  size_t x=0, y=0;

  puts("nao unitarios:");
  recv_idx=0; // reset
  for (size_t i=0; i < n_ops; ++i) 
  {
    scanf("x%zu = x%zu %c x%zu%*[^\n]", &recv_idx, &x, &op, &y); 
    scanf("%*c"); // consome \n
    assert(recv_idx == idx+1 && "Indice de valores fora de ordem");
    assert(x > 0 && x < (n_values + n_ops) && "'x' fora do intervalo");
    assert(y > 0 && y < (n_values + n_ops) && "'y' fora do intervalo");

    switch (op) {
    case '+': /* SOMA */
        intvs[idx].min = nextafterf(intvs[x-1].min + intvs[y-1].min, -INFINITY);
        intvs[idx].max = nextafterf(intvs[x-1].max + intvs[y-1].max, INFINITY);
        break;
    case '-': /* SUBTRAÇÃO */
        intvs[idx].min = nextafterf(intvs[x-1].min - intvs[y-1].max, -INFINITY);
        intvs[idx].max = nextafterf(intvs[x-1].max - intvs[y-1].min, INFINITY);
        break;
    case '*': { /* MULTIPLICAÇÃO */
        float values[] = {
          intvs[x-1].min * intvs[y-1].min,
          intvs[x-1].min * intvs[y-1].max,
          intvs[x-1].max * intvs[y-1].min,
          intvs[x-1].max * intvs[y-1].max
        };
        size_t size = sizeof(values) / sizeof(float);
        intvs[idx].min = nextafterf(faminf(values, size), -INFINITY);
        intvs[idx].max = nextafterf(famaxf(values, size), INFINITY);
        break; }
    case '/': { /* DIVISÃO */
        float values[] = {
          intvs[x-1].min / intvs[y-1].min,
          intvs[x-1].min / intvs[y-1].max,
          intvs[x-1].max / intvs[y-1].min,
          intvs[x-1].max / intvs[y-1].max
        };
        size_t size = sizeof(values) / sizeof(float);
        intvs[idx].min = nextafterf(faminf(values, size), -INFINITY);
        intvs[idx].max = nextafterf(famaxf(values, size), INFINITY);
        break; }
    default:
        fprintf(stderr, "Simbolo de operacao '%c' nao definido\n", op);
        return -1;
    }
    validate_interval(intvs[idx].min, intvs[idx].max);

    /* Imprime se for não unitário */
    if (!falmosteq(intvs[idx].min, intvs[idx].max, MAX_ULPS)) {
      printf("X%zu = [%1.8e, %1.8e]\n", recv_idx, intvs[idx].min, intvs[idx].max);
    }

    ++idx;
  }

  free(intvs);

  return 0;
}
