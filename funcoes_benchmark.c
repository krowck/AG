#include <math.h>
#define PI 3.14159265359

#define NVARS 6
/*
 * Este arquivo - "funcoes_bechmark.c" - contem a implementacao das funcoes do benchmark.
 *
 *
 * Para cada funcao de benchmark, temos a implementacao de duas funcoes computacionais:
 *  i) uma funcao (em C) que retorna um valor de uma determinada funcao de benchmark, dados x1 e x2;
 *  ii) uma funcao (em C) que "retorna" o limite dos valores de dominio (utilizando passagem de parametro por referencia)
 */

/*
 * Retorna um valor (do tipo float) referente a funcao RASTRIGIN, dado duas entrada, x1 e x2(tambem do tipo float)
 */
float rastrigin(double cromossomo[]){
    double obj = 0;
    int j;
    for(j = 0 ; j < NVARS; j++)
        {
                obj = obj + (pow(cromossomo[j],(double)2)-10*cos(2*PI*cromossomo[j])+10);
        }
    return obj;
}

/*
 * Implementacao da funcao que informa os valores de dominio referente a funcao de benchmark RASTRIGIN
 */
void d_rastrigin(float d[]){
    d[0] = -5.12;
    d[1] = 5.12;
}

