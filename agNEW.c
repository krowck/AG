/******************************************************************************************************************
*                                   Algoritmo desenvolvido para o TCC                                             *
*  para compilar:  gcc agNEW.c -lm -o main principal.c funcoes_benchmark.c gerador_numeros.c ./src-clstr/cluster.c -O3  *
*                                                                                                                 *
*                                                                                                                 *
******************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//#include <windows.h>
#include <assert.h>
#include <string.h>
#include <float.h>
#include <unistd.h>

#include "./src-clstr/cluster.h" 
#include "funcoes_benchmark.h"
#include "gerador_numeros.h"
#include "habitats.h"

#define TAM_TORNEIO 5
#define RUNS 10

double m_nmdf = 0;
double **distances;

/*
 * A estrutura (struct) "t_individuo" representa um unico individuo da populacao
 */
typedef struct {
    double gene[NVARS];
    double fitness;
    int index;
    double distanceFromCenter;
} t_individuo;

void merge(t_individuo vet[], int tam_vet) {

    int mid;
    int i, j, k;
    t_individuo* tmp;
    
    tmp = (t_individuo*) malloc(tam_vet * sizeof(t_individuo));
    
    if (tmp == NULL) {
        exit(1);
    }
    
    mid = tam_vet / 2;
    i = 0;
    j = mid;
    k = 0;
    
    while (i < mid && j < tam_vet) {
        if (vet[i].fitness <= vet[j].fitness) {
            tmp[k] = vet[i++];
        }
        else {
            tmp[k] = vet[j++];
        }
        ++k;
    }
    
    if (i == mid) {
        while (j < tam_vet) {
            tmp[k++] = vet[j++];
        }
    }
    else {
        while (i < mid) {
            tmp[k++] = vet[i++];
        }
    }
    
    for (i = 0; i < tam_vet; ++i) {
        vet[i] = tmp[i];
    }
    
    free(tmp);
}

/*
 * Ordenacao da Populacao pelo Fitness
 * O procedimento abaixo implementa o algoritmo de Ordenacao MERGE SORT.
 */
void mergeSort(t_individuo vet[], int tam_vet) {
    int mid;
    if (tam_vet > 1) {
        mid = tam_vet / 2;
        mergeSort(vet, mid);
        mergeSort(vet+mid, tam_vet-mid);
        merge(vet, tam_vet);
    }
}

float obter_fitness(int funcao, double cromossomo[]){
    float fitness = 0.0;
    switch (funcao)
    {
        case 1:
            fitness = rastrigin(cromossomo);
            break;
        case 2:
            fitness = quadratic(cromossomo);
            break;
        case 3:
            fitness = rosenbrock(cromossomo);
            break;
        case 4:
            fitness = schwefel(cromossomo);
            break;
        case 5:
            fitness = ackley(cromossomo);
            break;
        case 6:
            fitness = griewank(cromossomo);
            break;
        case 7:
            fitness = powell(cromossomo);
            break;
        case 8:
            fitness = zakharov(cromossomo);
            break;
        default:
            printf ("\nERRO!\n");
    }
    
    return fitness;
}

void identificar_dominio(int funcao, float *l_inf, float *l_sup){
    float d[2];
    switch (funcao)
    {
        case 1:
            d_rastrigin(d);
            break;
        case 2:
            d_quadratic(d);
            break;
        case 3:
            d_rosenbrock(d);
            break;
        case 4:
            d_schwefel(d);
            break;
        case 5:
            d_ackley(d);
            break;
        case 6:
            d_griewank(d);
            break;
        case 7:
            d_powell(d);
            break;
        case 8:
            d_zakharov(d);
            break;
        default:
            printf ("\nERRO!\n");
    }
    *l_inf = d[0];
    *l_sup = d[1];
}

/*
 * O procedimento abaixo eh responsavel por encontrar o melhor individuo (o que possui o menor fitness) no vetor "vet".
 */
void encontra_melhor_individuo(t_individuo vet[], int tam_vet, t_individuo *melhor){
    int i;
    *melhor = vet[0];
    for(i = 1; i < tam_vet; i++){
        if (vet[i].fitness < melhor->fitness){
            *melhor = vet[i];
        }            
    }
}


/*
 * Procedimento para imprimir um unico individuo do tipo t_individuo
 */
void imprimir_individuo(t_individuo individuo){
    int i;
    for(i = 0; i < NVARS; i++){
        printf("x[%d] = %f  ",i, individuo.gene[i]);
    }
    printf("fitness = %f\n",individuo.fitness);    
}


/*
 * Procedimento para imprimir um vetor de t_individuo
 */
void imprimir_populacao(t_individuo populacao[], int total_individuos){
    int i, j;
    for(i = 0; i < total_individuos; ++i){
        for (j = 0; j < NVARS; ++j)
        {
            printf("%f  ", populacao[i].gene[j]);
        }
        //printf("(%d) x1 = %f; x2 = %f\n",i,populacao[i].x1,populacao[i].x2);
        printf("\nfitness: %f\n",populacao[i].fitness);
        
        //printf("fitness: %d\n\n",obter_numero_uniforme_discreto(0,4)); //4 inclusive
    }
    printf("\n----------------------------------\n");
}

/*
 * Geracao da Populacao Inicial
 * O procedimento abaixo e responsavel por gerar a populacao inicial do algoritmo genetico.
 * Os parametros de entrada deste procedimento sao:
 *  - a populacao (vetor de "t_individuo")
 *  - o tamanho da populacao ("total_individuos")
 *  - o ID (ou codigo) da funcao a ser otimizada. Este ID eh necessario para se reconhecer a funcao
 *      facilitando a obtencao do dominio e do fitness de tal funcao.
 */
void gerar_populacao_inicial(t_individuo populacao[], int total_individuos, int funcao){

    float l_inf = 0.0;
    float l_sup = 0.0;
    identificar_dominio(funcao,&l_inf,&l_sup);

    int i, j;
    for(i = 0; i < total_individuos; ++i){
        for(j = 0; j < NVARS; j++){
            populacao[i].gene[j] =  obter_numero_uniforme_continuo(l_inf,l_sup);
            //printf("GENE  >>>>> %f  ", populacao[i].gene[j]);    
        }        
        //printf("\n");
        populacao[i].fitness = obter_fitness(funcao, populacao[i].gene);      
    }
    //imprimir_populacao(populacao, total_individuos);
}

void gerar_individuo(t_individuo individuo, int funcao)
{
    float l_inf = 0.0;
    float l_sup = 0.0;
    identificar_dominio(funcao,&l_inf,&l_sup);

    int i, j;
    for(j = 0; j < NVARS; j++)
    {
        individuo.gene[j] =  obter_numero_uniforme_continuo(l_inf,l_sup);
            //printf("GENE  >>>>> %f  ", populacao[i].gene[j]);    
    }        
       //printf("\n");
    individuo.fitness = obter_fitness(funcao, individuo.gene);      
}

/*
 * Operador de Mutacao
 * Este procedimento implementa o operador de mutacao.
 *
 * Os parametros de entrada deste procedimento sao:
 *  - o filho (ou descendente) gerado pela recombinacao de dois pais selecionados
 *  - a probabilidade (baixa) de mutacao informada pelo usuario no Menu (prob_mutacao)
 *  - o ID (ou codigo) da funcao a ser otimizada. Este ID eh necessario para se reconhecer a funcao
 *      facilitando a obtencao do dominio e do fitness de tal funcao.
 * 
 */
void op_mutacao(t_individuo *filho, double prob_mutacao, int funcao){
    int i;    
    for (i = 0; i < NVARS; ++i){   
        float u = obter_numero_uniforme();     
        if (u <= prob_mutacao){        
            float v = obter_numero_uniforme();
            if(v <= 0.5)
            {
                filho->gene[i] += obter_numero_uniforme();
            }
            else
            {
                filho->gene[i] -= obter_numero_uniforme();
            }
        }
    }
    filho->fitness = obter_fitness(funcao, filho->gene);     
}

/*
 * Operador de Recombinacao
 * O procedimento abaixo e responsavel por recombinar os codigos geneticos de ambos os pais
 * 
 * Os parametros de entrada deste procedimento sao:
 *  - os pais selecionados no torneio ("pai" e "mae")
 *  - o novo individuo (filho ou descendente) a ser gerado pela recombinacao dos dois pais selecionados
 *  - o ID (ou codigo) da funcao a ser otimizada. Este ID eh necessario para se reconhecer a funcao
 *      facilitando a obtencao do dominio e do fitness de tal funcao.
 */

void op_crossover(t_individuo *pai, t_individuo *mae, int funcao){
    int point = obter_numero_uniforme_discreto(1, NVARS-2);

    int i;
    double t;

    //printf("%d\n", point);

    for (i = 0; i <= point; i++)
    {
        t = pai->gene[i];
        pai->gene[i] = mae->gene[i];
        mae->gene[i] = t;
    }
    pai->fitness = obter_fitness(funcao, pai->gene);
    mae->fitness = obter_fitness(funcao, mae->gene);
}

/*
 * Operador de Selecao de Pais
 * No procedimento abaixo, candidatos a pais sao sorteados aleatoriamente, sendo que o candidato vencedor do torneio
 * sera aquele que possuir o melhor (menor) fitness.
 *
 * Os parametros de entrada deste procedimento sao:
 *  - a populacao (vetor de "t_individuo")
 *  - o tamanho da populacao ("total_individuos")
 *  - o ID (ou codigo) da funcao a ser otimizada. Este ID eh necessario para se reconhecer a funcao
 *      facilitando a obtencao do dominio e do fitness de tal funcao.
 *  - os pais a serem selecionados no torneio ("pai" e "mae"). [SELECAO POR TORNEIO, onde a letra grega "tau" = 3]
 */
void op_selecao_de_pais(t_individuo populacao[], int total_individuos, t_individuo *pai, t_individuo *mae){
    t_individuo sorteio[TAM_TORNEIO]; //valor eh 3 pois foi definido no enunciado, isto e, letra grega "tau" = 3
    int i;


    for (i = 0; i < TAM_TORNEIO; ++i)
    {
        sorteio[i] = populacao[obter_numero_uniforme_discreto(0,total_individuos-1)];
    }
    encontra_melhor_individuo(sorteio,TAM_TORNEIO,pai);


    for (i = 0; i < TAM_TORNEIO; ++i)
    {
        sorteio[i] = populacao[obter_numero_uniforme_discreto(0,total_individuos-1)];
    }
    encontra_melhor_individuo(sorteio,TAM_TORNEIO,mae);  
}

void selecao_aleatoria(t_individuo populacao[], int total_individuos, t_individuo *sorteio)
{
    *sorteio = populacao[obter_numero_uniforme_discreto(0,total_individuos-1)];
}

/*
apenas atualizar o vetor de população com os novos individuos gerador a partir de crossover e mutação
*/
void op_selecao_de_sobreviventes(t_individuo populacao[], int total_individuos, t_individuo novos_individuos[], int comecar){
    int i = 0;
    int j = 0;
    for(i = comecar; i < total_individuos + comecar; i++){
        populacao[i] = novos_individuos[j];
        j++;
    }

}
/*
double diversity_population(t_individuo pop[], int tamPopulation)
{

    double diversity = 0;
    double aux_1 = 0;
    double aux_2 = 0;
    unsigned short int a = 0;
    unsigned short int b = 0;
    unsigned short int d = 0;
    for(a = 0; a < tamPopulation; a++)
    {
        for(b = (a+1); b < tamPopulation; b++)
        {
            aux_1 = 0;
            for(d = 0; d < NVARS; d++)
            {       
                aux_1 += pow(pop[a].gene[d] - pop[b].gene[d], 2);
            }
            if(b == (a+1) || aux_2 > aux_1)
            {
                aux_2 = aux_1;
            }
        }
        diversity += log((double)1.0 + aux_2);  
    } 
    if(m_nmdf < diversity)
    {
       m_nmdf = diversity;
    }
    return diversity / m_nmdf;
}
*/


double diversity_population(t_individuo pop[], int tamPopulation)
{
    //VARIÁVEIS LOCAIS
    double diversity = 0;
    double aux_1 = 0;
    double aux_2 = 0;
    unsigned short int a = 0;
    unsigned short int b = 0;
    unsigned short int d = 0;
    for (a = 0; a < tamPopulation; a++) 
    {
        for (b = (a + 1); b < tamPopulation; b++)
        {
            aux_1 = 0;
            for (d = 0; d < NVARS; d++) //SOMATÓRIO VARIANDO OS GENES ATÉ O NUMERO DE DIMENSÕES  
            {
                aux_1 += pow(pop[a].gene[d] - pop[b].gene[d], 2);
            }
            aux_1 = sqrt(aux_1); //RAIZ QUADRADA DO RESULTADO DO SOMATÓRIO
            aux_1 = aux_1 / NVARS; //DIVIDE RESULTADO PELO NUMERO DE DIMENSÕES
            if (b == (a + 1) || aux_2 > aux_1) 
            {
                aux_2 = aux_1; //ATRIBUI O MENOR RESULTADO A AUX_2
            }
        }
        diversity += log((double)1.0 + aux_2); //CALCULA O LOGARITMO NATURAL DO RESULTADO SOMADO DE UM
    }
    //ATRIBUI O RESULTADO A N_NMDF 
    if (m_nmdf < diversity)
    {
        m_nmdf = diversity; 
    }
    return diversity / m_nmdf;
}

double sum_array(t_individuo pop[], int num_elements)
{
   int i;
   double sum = 0;
   for (i=0; i<num_elements; i++)
   {
     sum = sum + pop->gene[i];
   }
   return(sum);
}

double euclidiana(double ind1, double ind2)
{
    double soma = 0.0;
    int i;
    for (i = 0; i < NVARS; ++i)
    {
        soma += pow((ind1 - ind2),2);
    }
   return sqrt(soma);
}


double find_minimum_value(double a[], int n) 
{
    int c;
    double min;

    min = a[0];
    int index = 0;
    //printf("%f\n", min);
    for (c = 1; c < n; c++) {
        if (a[c] < min) {
            index = c;
            min = a[c];
        }
    }
    return min;
}

int find_minimum_index(double a[], int n)
{
    int c;
    double min;

    min = a[0];
    int index = 0;
    //printf("%f\n", min);
    for (c = 1; c < n; c++) {
        if (a[c] < min) {
            index = c;
            min = a[c];
        }
    }
    return index; 
}




void shuffle(int *array, int n)
{
    if (n > 1) 
    {
        int i;
        for (i = 0; i < n - 1; i++) 
        {
          int j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}

double find_maximum(double a[], int n)
{
    int c;
    double max, index;

    max = a[0];
    index = 0;

    for (c = 1; c < n; c++) {
        if (a[c] > max) {
            index = c;
            max = a[c];
        }
    }
    return index;
}

double** distancia_euclidiana(t_individuo pop[], int total_individuos)
{
    double media = 0.0;
    double aux = 0;
    double euclid_max = 0;
    double temp = 0;
    int i, j, k;
    unsigned short int a = 0;

    distances = (double**) malloc(total_individuos*sizeof(double*));
    for(a = 0; a < total_individuos; a++)
    {
        distances[a] = (double*) malloc(total_individuos*sizeof(double));
    }

    for (i = 0; i < total_individuos; ++i)
    {
        distances[i][i] = 0;
        for (j = i+1; j < total_individuos; ++j)
        {        
            for (k = 0; k < NVARS; ++k)
            {
                media += pow(fabs(pop[i].gene[k] - pop[j].gene[k]), 2.0);
            }
            distances[i][j] = sqrt(media);
            distances[j][i] = distances[i][j];
            if (euclid_max < distances[i][j])
            {
                euclid_max = distances[i][j];
            }
            media = 0.0;
        }
    }
    for(i = 0; i < total_individuos; i++)
    {
        for(j = i +1; j < total_individuos; j++)
        {
            distances[i][j] = distances[i][j] / euclid_max;
            distances[j][i] = distances[i][j];
        }
    }
    return distances;
}


void clusterAnalysis(t_individuo populacao[], double clusterRadius, int total_individuos, double *clusterID, double *distanceFromCenter, double clusterCenter[total_individuos][NVARS], int numberOfClusters)
{   
    int numberOfChromosomes = NVARS;
    int count, index = 0;
    double EuclideanDistance[numberOfClusters];
    double aux, media, min = 0.0;
    double Y;
    int i, j, k;
    double **distancia;

    initHabitats(total_individuos);
    distancia = distancia_euclidiana(populacao, total_individuos);
    singlelink(total_individuos, NVARS, distancia);
    printDendogram(total_individuos);
    buildHabitats(total_individuos, distancia);

    for (i = 0; i < g_habitatsSize; ++i)
    {
        for (j = 0; j < H[i].h_ind_count; ++j)
        {
            populacao[H[i].h_ind[j]].index = i;
        }       
    }
}

void improveCluster(t_individuo populacao[], double *clusterID, int alvo, double *distanceFromCenter, double killIndex, int funcao, int total_individuos, 
    double prob_mutacao, int *comecar, t_individuo populacao_aux[], t_individuo melhores[])
{
    double distanceMatrix[total_individuos];
    int numberOfChromosomes = 0;
    double u;
    int flag = 0, cont = 0, cont2 = 0;
    t_individuo pai;
    t_individuo mae;
    t_individuo pop_aux[total_individuos];
    t_individuo melhor;

    int i, j, aux = 0;
    for (i = 0; i < total_individuos; ++i)
    {
        if (populacao[i].index == alvo)
        {
            for (j = 0; j < NVARS; ++j)
            {
                pop_aux[aux].gene[j] = populacao[i].gene[j];
                pop_aux[aux].fitness = populacao[i].fitness;
            }
            distanceMatrix[aux] = distanceFromCenter[i];
            aux++;
            numberOfChromosomes++;
        }
    }  

    //printf("%d\n", numberOfChromosomes);

    double percentageChromosomesToCrossover = 0.7;
    int totalChromosomesToCrossover = round(percentageChromosomesToCrossover * numberOfChromosomes);

    int vetor_aux[numberOfChromosomes];
    int vetor_pai[numberOfChromosomes];
    int vetor_naopai[numberOfChromosomes - totalChromosomesToCrossover];
    t_individuo novos_individuos[aux];

    for (i = 0; i < numberOfChromosomes; ++i)
    {
        vetor_aux[i] = i;
    }

    shuffle(vetor_aux, numberOfChromosomes);
    //printf("Total de individuos: %d  & Total de pais: %d\n", numberOfChromosomes, totalChromosomesToCrossover);

    for (i = 0; i < totalChromosomesToCrossover; ++i)
    {
        vetor_pai[i] = vetor_aux[i];
        cont2++;
        //printf("VETOR DE PAIS: %d  %d\n", vetor_pai[i], i);
    }
    //printf("\n\n");
    for (i = cont2; i < numberOfChromosomes; ++i)
    {
        vetor_pai[i] = vetor_aux[i];
        //printf("VETOR DE NAO PAIS: %d  %d\n", vetor_pai[i], i);
    }

    if (totalChromosomesToCrossover%2 == 1 && totalChromosomesToCrossover > 2)
    {
        totalChromosomesToCrossover--;
        flag = 1;
    }
    if (totalChromosomesToCrossover == 1)
    {
        novos_individuos[0] = pop_aux[vetor_pai[totalChromosomesToCrossover-1]];
        //op_selecao_de_sobreviventes(populacao_aux, totalChromosomesToCrossover, novos_individuos, *comecar);
    }
    else
    {
        for (i = 0; i < totalChromosomesToCrossover; ++i)
        {      
            pai = pop_aux[vetor_pai[i]];
            mae = pop_aux[vetor_pai[i+1]];
            //u = obter_numero_uniforme();
            //if (u < 0.9)
            //{
                //printf("%f\n", u);
            op_crossover(&pai, &mae, funcao);
            //}
            op_mutacao(&pai,prob_mutacao,funcao);
            op_mutacao(&mae,prob_mutacao,funcao);
            novos_individuos[i] = pai;
            i++;
            novos_individuos[i] = mae;
            cont += 2;
        }
    }
    
    if (flag == 1)
    {
        novos_individuos[cont] = pop_aux[vetor_pai[totalChromosomesToCrossover]];
        //op_selecao_de_sobreviventes(populacao_aux, totalChromosomesToCrossover+1, novos_individuos, *comecar);
        totalChromosomesToCrossover++;
        //*comecar += totalChromosomesToCrossover;
        //int flag = 0;
    }

    for (i = cont2; i < numberOfChromosomes; ++i)
    {
        novos_individuos[i] = pop_aux[vetor_pai[i]];
    }

    op_selecao_de_sobreviventes(populacao_aux, numberOfChromosomes, novos_individuos, *comecar);

    *comecar += numberOfChromosomes;    

    if(numberOfChromosomes != 0){
        encontra_melhor_individuo(novos_individuos, numberOfChromosomes, &melhor);
        melhores[alvo] = melhor;    
    }/*
    else
    {
        gerar_individuo(melhor, funcao);
        imprimir_individuo(melhor);
        melhores[alvo] = melhor;
    }*/
}

void improveIdeals(t_individuo melhores[], int maxClusters, int funcao, double prob_mutacao)
{
    t_individuo pai;
    t_individuo mae;
    t_individuo novos_individuos[maxClusters];
    int i, cont = 0, flag = 0;

    if (maxClusters % 2 == 1 && maxClusters > 2)
    {
        maxClusters--;
        flag = 1;
    }

    for (i = 0; i < maxClusters; ++i)
    {
        pai = melhores[i];
        mae = melhores[i+1];
        //u = obter_numero_uniforme();
        //if (u < 0.9)
        //{
        op_crossover(&pai, &mae, funcao);
        //}
        op_mutacao(&pai,prob_mutacao,funcao);
        op_mutacao(&mae,prob_mutacao,funcao);
        novos_individuos[i] = pai;
        i++;
        novos_individuos[i] = mae;
        cont += 2;
    }

    if (flag == 1)
    {
        novos_individuos[cont] = melhores[cont];
        maxClusters++;
    }

    //imprimir_populacao(novos_individuos, maxClusters);

    op_selecao_de_sobreviventes(melhores, maxClusters, novos_individuos, 0);
}

void generateNextPopulation(t_individuo populacao[], t_individuo melhores[], int total_individuos, int maxClusters, double prob_mutacao, int funcao)
{
    int i, j;
    t_individuo novos_individuos[total_individuos];
    int vetor_aux[total_individuos + maxClusters];
    int vetor_pai[total_individuos + maxClusters];
    t_individuo mutado;
    int flag = 0;
    float l_inf = 0.0;
    float l_sup = 0.0;
    identificar_dominio(funcao,&l_inf,&l_sup);

    for (i = 0; i < total_individuos + maxClusters; ++i)
    {
        vetor_aux[i] = i;
        //printf("VETOR: %d  I: %d\n", vetor_aux[i], i);
    }
    shuffle(vetor_aux, total_individuos + maxClusters);

    for (i = 0; i < total_individuos; ++i)
    {
        vetor_pai[i] = vetor_aux[i];
        //printf("VETOR: %d  I: %d\n", vetor_pai[i], i);
    }
    for (i = 0; i < total_individuos; ++i)
    {
        if (vetor_pai[i] >= total_individuos)
        {
            for (j = 0; j < NVARS; ++j)
            {
                if (melhores[(total_individuos+maxClusters-1) - vetor_pai[i]].gene[j] < l_inf*2 || melhores[(total_individuos+maxClusters-1) - vetor_pai[i]].gene[j] > l_sup*2)
                {
                    flag = 1;
                }
            }
            if (flag == 0)
            {
                mutado = melhores[(total_individuos+maxClusters-1) - vetor_pai[i]];
                op_mutacao(&mutado, prob_mutacao, funcao);
                novos_individuos[i] = mutado;
            }
            else
            {
                mutado = populacao[vetor_pai[(total_individuos+maxClusters-1) - vetor_pai[i]]];
                op_mutacao(&mutado, prob_mutacao, funcao);
                novos_individuos[i] = mutado;
            }            
        }
        else
        {
            mutado = populacao[vetor_pai[i]];
            op_mutacao(&mutado, prob_mutacao, funcao);
            novos_individuos[i] = mutado;
        }        
    }
    //imprimir_populacao(novos_individuos, total_individuos);
    op_selecao_de_sobreviventes(populacao, total_individuos, novos_individuos, 0);
    //imprimir_populacao(populacao, total_individuos);
}





/*
 * Evolucao da Populacao
 * No procedimento "executar" que e realizada a Evolucao da Populacao.
 
 * Os parametros de entrada (definidos pelo usuario) deste procedimento sao:
 *  - o codigo (ou ID) da funcao a ser otimizada
 *  - o tamanho da populacao (ou total de individuos "total_individuos" da populacao)
 *  - por quantas "geracoes" a populacao inicial sera evoluida
 *  - a probabilidade de mutacao (prob_mutacao)
 */
void executar(int funcao, int total_individuos, int geracoes){
    srand((unsigned)time(NULL)); 

    int run;
    int i;
    int j;
    double vet_melhores[RUNS][geracoes];
    double vet_diversidade[RUNS][geracoes];
    struct timeval timevalA;
    struct timeval timevalB;

    FILE *fpMedia;
    FILE *fpDiversidade;
    FILE *fp;
    FILE *fpNumeroCluster;


    fpMedia = fopen("mediaGeracoes.txt", "w+");
    fpDiversidade = fopen("mediaDiversity.txt", "w+");
    fp = fopen("output.txt", "w+");
    fpNumeroCluster = fopen("NumeroDeCluster.txt", "w+");
    gettimeofday(&timevalA,NULL);
    for (run = 0; run < RUNS; ++run)
    {
        double startingClusterRadius = 0.3;
        double minimumClusterRadiusThreshold = 0.1;
        double clusterRadius = startingClusterRadius;
        double numberOfLessThan2 = 0;
        double numberOfGreaterThan6 = 0;
        double firstnessIndex = 2;
        double killIndex = firstnessIndex + 1;
        int populationOverFlowFlag = 0;
        double stagnationCounter = 0;
        double averageMinimumValue = 0;
        int maxClusters = 0;
        int flagCluster = 0;
        double prob_mutacao = 0.2;
        double minMutationRate = 0.002;
        m_nmdf = 0.0;        
        
        /*
         * A Populacao e representada como um vetor de "t_individuo", cujo o tamanho e "total_individuos" (definido previamente pelo usuario).
         */
        t_individuo populacao[total_individuos];
        t_individuo pop_teste[total_individuos];
        
        gerar_populacao_inicial(populacao, total_individuos, funcao);

        int g = 0; //contador de geracoes
        int j = 0, i, k;
        for(g = 0; g < geracoes; g++){
            double *newClusterData;
            double *newClusterDataFitness;
            double *ideals;
            double clusterID[total_individuos];
            double clusterCenter[total_individuos][NVARS];
            double distanceFromCenter[total_individuos];
            int comecar = 0;

            

            // if (maxClusters < 6 && flagCluster == 0)
            // {
            //     maxClusters++;
            // }
            // if (maxClusters == 6)
            // {
            //     flagCluster = 1;
            // }
            // if (maxClusters > 2 && flagCluster == 1)
            // {
            //     maxClusters--;
            // }
            // if (maxClusters == 2)
            // {
            //     flagCluster = 0;
            // }
            // //printf("%d\n", maxClusters);

            t_individuo populacao_aux[total_individuos];
            t_individuo best;

            encontra_melhor_individuo(populacao, total_individuos, &best);

            /*
            int *c = k_means(populacao, total_individuos, NVARS, 6, 0.0001, 0);
            for (i = 0; i < total_individuos; ++i)
            {
                printf("data point %d is in cluster %d\n", i, c[i]);
            }   
            */

            clusterAnalysis(populacao, clusterRadius, total_individuos, clusterID, distanceFromCenter, clusterCenter, maxClusters);

            maxClusters = g_habitatsSize;
            t_individuo melhores[maxClusters];

            for (i = 0; i < maxClusters; ++i)
            {
                improveCluster(populacao, clusterID, i, distanceFromCenter, killIndex, funcao, total_individuos, prob_mutacao, &comecar, populacao_aux, melhores);
            }

            improveIdeals(melhores, maxClusters, funcao, prob_mutacao);   

            memcpy(populacao, populacao_aux, sizeof(t_individuo)*total_individuos);


            generateNextPopulation(populacao, melhores, total_individuos, maxClusters, prob_mutacao, funcao);

            double diversidade = diversity_population(populacao, total_individuos);

            //imprimir_populacao(populacao, total_individuos);
            //mergeSort(populacao, total_individuos);
            //imprimir_populacao(populacao, total_individuos);

            populacao[total_individuos-1] = best;

            encontra_melhor_individuo(populacao, total_individuos, &best);
            // Saida de Dados
            fprintf(fp, "%d %f %f\n", g, best.fitness, diversidade);
            //printf("%d\t%f  ",g,best.fitness); //saida de dados
            //printf("  %f\n", diversidade);
            //printf("%lf\n", prob_mutacao);
            prob_mutacao -= 0.001;

            if (prob_mutacao <= minMutationRate)
            {
                prob_mutacao = minMutationRate;
            }
            vet_melhores[run][g] = best.fitness;
            vet_diversidade[run][g] = diversidade;
            fprintf(fpNumeroCluster, "%d\n", g_habitatsSize);
        }
        printf("%d\n", run);
        //fprintf(fp, "-------------------------------");    
        //sleep(10);
    }
    gettimeofday(&timevalB,NULL);

    printf("\n\ntempo de execucao: %lf\n",  (timevalB.tv_sec-timevalA.tv_sec+(timevalB.tv_usec-timevalA.tv_usec)/(double)1000000)/10.0);

    double mediabest[geracoes];
    double mediaDiversity[geracoes];


    for (i = 0; i < geracoes; ++i)
    {   
        mediabest[i] = 0.0;
        mediaDiversity[i] = 0.0;
    }

    for (i = 0; i < geracoes; ++i)
    {
        for (j = 0; j < RUNS; ++j)
        {
            mediabest[i] += vet_melhores[j][i];
            mediaDiversity[i] += vet_diversidade[j][i];
        }
        mediaDiversity[i] = mediaDiversity[i]/(double)RUNS;
        mediabest[i] = mediabest[i]/(double)RUNS;
        fprintf(fpMedia, "%f\n", mediabest[i]);
        fprintf(fpDiversidade, "%f\n", mediaDiversity[i]);
    }

    fclose(fp);
    fclose(fpMedia);
    fclose(fpDiversidade);
    fclose(fpNumeroCluster);

}

