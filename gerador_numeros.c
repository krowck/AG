#include <stdlib.h>
#include <time.h>

/*
 * Este arquivo - "gerador_numeros.c" - contem a implementacao de funcoes que geram aleatoriamente:
 *  - numeros uniforme;
 *  - numeros uniforme continuo;
 *  - numeros uniforme discreto.
 *
*/

 
/*
 * Esta funcao retorna um valor (do tipo float) aleatorio uniforme entre 0 e 1.
 */
float obter_numero_uniforme(){
    int max_mix = rand() % 100; //aumenta a aleatoriedade
    int i;
    for (i = 0; i < max_mix; i++) rand(); //aumenta a aleatoriedade
    
    return (float) rand() / (float) RAND_MAX;
}

/*
 * Esta funcao possui dois parametros de entrada:
 *  - limite inferior
 *  - limite superior
 *
 * Retorna um valor (do tipo float) aleatorio uniforme continuo.
 */
float obter_numero_uniforme_continuo(float l_inf, float l_sup){
    float u = obter_numero_uniforme();
    return l_inf + (l_sup - l_inf) * u;
}

/*
 * Esta funcao possui dois parametros de entrada:
 *  - limite inferior
 *  - limite superior
 *
 * Retorna um valor (do tipo int) aleatorio uniforme discreto.
 */
int obter_numero_uniforme_discreto(float l_inf, float l_sup){
    float u = obter_numero_uniforme();
    
    return l_inf + (int)((l_sup - l_inf + 1.0) * u);
}

int nextInt(int n)
{
	return (int) ((double)0.0 + ((n - 0.0)*rand()/(RAND_MAX+1.0)));
}



double nextDoubleA(double low, double high)
{
	return low + ((double)rand() / RAND_MAX) * ( high - low );
}

double nextDouble()
{
	double result = nextDoubleA(0,1);
	return result;
}