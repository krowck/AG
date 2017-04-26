#define PI 3.14159265359

#define NVARS 30

/*
 * Este arquivo - "funcoes_bechmark.h" - contem apenas prototipos de funcoes.
 *
 */

float rastrigin(double cromossomo[]);

float quadratic(double cromossomo[]);

float rosenbrock(double cromossomo[]); 

float schwefel(double cromossomo[]);

float ackley(double cromossomo[]);

float griewank(double cromossomo[]);

float powell(double cromossomo[]);

float zakharov(double cromossomo[]);

void d_quadratic(float d[]);

void d_rastrigin(float d[]);

void d_rosenbrock(float d[]);

void d_schwefel(float d[]);

void d_ackley(float d[]);

void d_griewank(float d[]);

void d_powell(float d[]);

void d_zakharov(float d[]);
