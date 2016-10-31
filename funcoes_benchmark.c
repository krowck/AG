#include <math.h>
#define PI 3.14159265359

#define NVARS 30
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
                obj += (pow(cromossomo[j],(double)2)-cos(2*PI*cromossomo[j]));
        }
    return NVARS + obj;
}

float quadratic(double cromossomo[]){
	double obj = 0;
	int j;
	for (j = 0; j < NVARS; ++j)
	{
		obj += (pow(10, j-1) * pow(cromossomo[j],2));
	}
	return obj;
}

float rosenbrock(double cromossomo[])
{
	double obj = 0;
	int j;
	obj = (100 * pow((cromossomo[1] - pow(cromossomo[0], 2) ),2) ) + pow((1 - cromossomo[0]),2);
	return obj;
}

float schwefel(double cromossomo[])
{
	double obj = 0;
	int j;
	for (j = 0; j < NVARS; ++j)
	{
		obj += (-(cromossomo[j])) * sin(sqrt(fabs(cromossomo[j])));
	}
	return obj;
}
/*
 * Implementacao da funcao que informa os valores de dominio referente a funcao de benchmark RASTRIGIN
 */
void d_rastrigin(float d[])
{
    d[0] = -5.12;
    d[1] = 5.12;
}

void d_quadratic(float d[])
{
	d[0] = -10.0;
	d[1] = 10.0;
}

void d_rosenbrock(float d[])
{
	d[0] = -1.5;
	d[1] = 1.5;
}

void d_schwefel(float d[])
{
	d[0] = -500.0;
	d[1] = 500.0;
}
