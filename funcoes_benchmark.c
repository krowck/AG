#include <math.h>
#include <stdio.h>

#define PI 3.14159265359

#define NVARS 50

double obj;
double aux, aux1;
/*
 * Este arquivo - "funcoes_bechmark.c" - contem a implementacao das funcoes do benchmark.
 *
 *
 * Para cada funcao de benchmark, temos a implementacao de duas funcoes computacionais:
 *  i) uma funcao (em C) que retorna um valor de uma determinada funcao de benchmark
 *  ii) uma funcao (em C) que "retorna" o limite dos valores de dominio (utilizando passagem de parametro por referencia)
 */

float rastrigin(double cromossomo[]){
    obj = 0;
    int j;
    for(j = 0 ; j < NVARS; j++)
        {
                obj += (pow(cromossomo[j],(double)2)-cos(2*PI*cromossomo[j]));
        }
    if (NVARS + obj < 0.01)
    {
    	printf("obj %lf\n", obj);
    }
    return NVARS + obj;
}

float quadratic(double cromossomo[]){
	obj = 0;
	int j;
	for (j = 0; j < NVARS; ++j)
	{
		obj += (pow(10, j-1) * pow(cromossomo[j],2));
	}
	return obj;
}

float rosenbrock(double cromossomo[])
{
	int i;
	obj = 0;
	for (i = 0; i < NVARS-1; i++)
        {
        	obj=obj+100.*pow((cromossomo[i+1] - pow(cromossomo[i],2.)),2) + pow((cromossomo[i] - 1.),2);
        }
	return obj;
}

float schwefel(double cromossomo[])
{
	obj = 0;
	int j;
	for (j = 0; j < NVARS; ++j)
	{
		obj += ((cromossomo[j])) * sin(sqrt(fabs(cromossomo[j])));
	}
	return 418.9829*NVARS - obj;
}

float ackley(double cromossomo[])
{
	obj = 0;
	aux = 0;
	aux1 = 0;
	int i;
	for (i = 0; i < NVARS; ++i)
	{
		aux += cromossomo[i]*cromossomo[i];
	}
	for (i = 0; i < NVARS; ++i)
	{
		aux1 += cos(2.0*PI*cromossomo[i]);
	}

	obj =  (-20.0*(exp(-0.2*sqrt(1.0/(float)NVARS*aux)))-exp(1.0/(float)NVARS*aux1)+20.0+exp(1));

	return obj;
}

float griewank(double cromossomo[])
{	
	aux = 0;
	aux1 = 1;
	int j;
	for(j = 0 ; j < NVARS; j++)
        {
		aux=aux+pow((cromossomo[j]),(double)2);
        	aux1=aux1*cos((((cromossomo[j])/sqrt((double)(j+1)))*PI)/180);
        }
	obj =  (1.0/(double)4000.0)*aux-aux1+1.0;

	return obj;
}

float powell(double cromossomo[])
{
	obj = 0;
	unsigned short int j;
	for (j = 1; j <= (int)NVARS/4; j++) {
		obj +=  pow(cromossomo[4*j-4] + 10 * cromossomo[4*j-3],2.0)
        	  + 5 * pow(cromossomo[4*j-2] - cromossomo[4*j-1],2.0)
        	  + pow(cromossomo[4*j-3] - 2 * cromossomo[4*j-2], 4.0)
		  + 10 * pow(cromossomo[4*j - 4] - cromossomo[4*j-1], 4.0);
	}
	return obj;
}

float zakharov(double cromossomo[])
{
	aux = 0;
	aux1 = 0;
	unsigned short int i;
	for (i = 0; i < NVARS; i++)
        {
        	aux = aux + pow(cromossomo[i],2.0);
	   	aux1 = aux1+0.5*i*cromossomo[i];
        }
	obj =  ( aux+pow(aux1,2.0)+pow(aux1,4.0) );

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
	d[0] = -1.0;
	d[1] = 1.0;
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

void d_ackley(float d[])
{
	d[0] = -32.0;
	d[1] = 32.0;
}

void d_griewank(float d[])
{
	d[0] = -600.0;
	d[1] = 600.0;
}

void d_powell(float d[])
{
	d[0] = -4.0;
	d[1] = 5.0;
}

void d_zakharov(float d[])
{
	d[0] = -5.0;
	d[1] = 10.0;
}	

