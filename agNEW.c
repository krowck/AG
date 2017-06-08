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
#include <windows.h>
#include <assert.h>
#include <string.h>
#include <float.h>
#include <unistd.h>
#include <sys/time.h>

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
    double distance;
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
        case 9:
            fitness = michaelewicz(cromossomo);
            break;
        case 10:
            fitness = dixonprice(cromossomo);
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
        case 9:
            d_michaelewicz(d);
            break;
        case 10:
            d_dixonprice(d);
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

void encontra_pior_individuo(t_individuo vet[], int tam_vet, t_individuo *melhor){
    int i;
    *melhor = vet[0];
    for(i = 1; i < tam_vet; i++){
        if (vet[i].fitness > melhor->fitness){
            *melhor = vet[i];
        }            
    }
}

double encontra_media_populacao(t_individuo pop[], int total_individuos)
{
    int i;
    double sum = 0;
    for (i = 0; i < total_individuos; ++i)
    {
        sum += pop[i].fitness;
    }
    return sum/total_individuos;
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

void gerar_individuo(t_individuo *individuo, int funcao)
{
    float l_inf = 0.0;
    float l_sup = 0.0;
    identificar_dominio(funcao,&l_inf,&l_sup);

    int j;
    for(j = 0; j < NVARS; j++)
    {
        individuo->gene[j] =  obter_numero_uniforme_continuo(l_inf,l_sup); 
    }        
       //printf("\n");
    individuo->fitness = obter_fitness(funcao, individuo->gene);      
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
    float l_inf = 0.0;
    float l_sup = 0.0;
    identificar_dominio(funcao,&l_inf,&l_sup);
    for (i = 0; i < NVARS; ++i){   
        float u = obter_numero_uniforme();     
        if (u <= prob_mutacao){        
            float v = obter_numero_uniforme();
            if(v <= 0.5)
            {
                filho->gene[i] += obter_numero_uniforme();
                if (filho->gene[i] > l_sup)
                {
                    filho->gene[i] = l_sup;
                }
            }
            else
            {
                filho->gene[i] -= obter_numero_uniforme();
                if (filho->gene[i] < l_inf)
                {
                    filho->gene[i] = l_inf;
                }
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
    int point = obter_numero_uniforme_discreto(1, NVARS-1);

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

void op_uniformcrossover(t_individuo *pai, t_individuo *mae, int funcao)
{
    int i;
    double t;

    //printf("%d\n", point);

    for (i = 0; i < NVARS; i++)
    {
        double prob = obter_numero_uniforme();
        if (prob < 0.5)
        {
            t = pai->gene[i];
            pai->gene[i] = mae->gene[i];
            mae->gene[i] = t;
        }
        
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
    t_individuo sorteio[TAM_TORNEIO];
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

double** distancia_euclidiana(t_individuo pop[], int total_individuos, int g)
{
    double media = 0.0;
    double euclid_max = 0;
    int i, j, k;
    unsigned short int a = 0;
    //FILE *fep = fopen("distancias.txt", "a");

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

    // for(i = 0; i < total_individuos; i++)
    // {
    //     for(j = 0; j < total_individuos; j++)
    //     {
    //         if (g == 100)
    //         {
    //             fprintf(fep, "%lf ", distances[i][j]);
    //         }        
    //     }
    //     fprintf(fep, "\n");
    // }
    // fclose(fep);
    return distances;
}


void op_selecao_de_sobreviventes(t_individuo populacao[], int total_individuos, t_individuo novos_individuos[], int comecar){
    int i = 0;
    int j = 0;
    for(i = comecar; i < total_individuos + comecar; i++){
        populacao[i] = novos_individuos[j];
        j++;
    }
}

void freeDistances(int total_individuos)
{    
    unsigned short int a;
    for(a = 0; a < total_individuos; a++)
    {
        free(distances[a]);
    }
    free(distances);
}

void verify_ALL(t_individuo populacao[], int total_individuos, int g, int gg)
{
    int i, j;
    int flag = 0;
    for (i = 0; i < total_individuos; ++i)
    {
        for (j = 0; j < NVARS; ++j)
        {
            if(populacao[i].gene[j] == 0.000000)
            {
                flag = 1;
            }
        }
        if (flag == 1)
        {
            printf("Geracao: %d INDIVIDUO: %d WARNING, NUMBER DA TRETA: %d CLUSTER: %d\n",g, i, gg, populacao[i].index);
            printf("Individuo = %d\n", i);
            imprimir_individuo(populacao[i]);
        }
    }
    if (flag == 1)
    {
        
        Sleep(5000);
    }
}

void clusterAnalysis(t_individuo populacao[], int total_individuos, int geracoes, double u)
{   
    int i, j, k;
    double **distancia;
    double *aux = (double*) malloc(total_individuos*sizeof(double)) ;
    int *index = (int*) malloc(total_individuos*sizeof(int));
    int g_habitatsSizeAntes = 0;
    //printf("%d\n", geracoes);
    
    initHabitats(total_individuos);

    distancia = distancia_euclidiana(populacao, total_individuos, geracoes);

    singlelink(total_individuos, NVARS, distancia, aux);

    buildHabitats(total_individuos, distancia);
    //printDendogram(total_individuos, geracoes);

    for (i = 0; i < total_individuos; ++i)
    {
        populacao[i].distance = aux[i];
        index[i] = 0;
    }

    for (i = 0; i < g_habitatsSize; ++i)
    {
        for (j = 0; j < H[i].h_ind_count; ++j)
        {
            populacao[H[i].h_ind[j]].index = i;
            //printf("%d = %d\n", H[i].h_ind[j], populacao[H[i].h_ind[j]].index);
        }       
    }

    
    for (i = 0; i < g_habitatsSize-1; ++i)
    {        
        for (k = i+1; k < g_habitatsSize; ++k)
        {
            if(fabs(populacao[H[i].h_ind[0]].distance - populacao[H[k].h_ind[0]].distance) < u)
            {
                for (j = 0; j < H[k].h_ind_count; ++j)
                {
                    populacao[H[k].h_ind[j]].index = i;
                    H[i].h_ind_count++;
                    H[i].h_ind[H[k].h_ind[j]] = H[k].h_ind[j];
                    //flag[H[k].h_ind[j]]++;
                    //printf("%d   %d    %d = %d\n",H[i].h_ind[H[k].h_ind[j]], H[k].h_ind[j], H[k].h_ind[j], populacao[H[k].h_ind[j]].index);
                    //H[k].h_ind[j] = -1;
                }
                H[k].h_ind_count = 0;
            }
        }     
    }   

    g_habitatsSizeAntes = g_habitatsSize;
    //printf("Habitats antes = %d\n",g_habitatsSizeAntes);
    for (i = 0; i < total_individuos; ++i)
    {
        index[populacao[i].index]++;
    }
    g_habitatsSize = 0;
    for (i = 0; i < total_individuos; ++i)
    {
        //printf("index = %d\n", index[i]);
        if(index[i] != 0)
        {
            g_habitatsSize++;
        }
    }
    //printf("Habitats depois = %d\n",g_habitatsSize );

    // for (k = 0; k < g_habitatsSizeAntes - g_habitatsSize; ++k)
    // {
    //     for (i = 0; i < g_habitatsSizeAntes; ++i)
    //     {
    //         if(index[i] == 0)
    //         {
    //             for (j = 0; j < H[i+1].h_ind_count; ++j)
    //             {
    //                 H[i].h_ind_count++;
    //                 H[i].h_ind[j] = H[i+1].h_ind[j];
    //             }
    //             H[i+1].h_ind_count=0;
    //             index[i] = H[i].h_ind_count;
    //             index[i+1] = 0;
    //         }       
    //     }
    // }

    for (i = 0; i < total_individuos; ++i)
    {
        if (populacao[i].index > g_habitatsSize-1)
        {
            for (k = 0; k < g_habitatsSize; ++k)
            {
                if(index[k] == 0 || populacao[i].index == index[k])
                {
                    populacao[i].index = index[k];
                    index[k] = populacao[i].index;
                }
            }
        }
    }

    free(index);
    free(aux);
    freeDistances(total_individuos);
    destroyHabitats(total_individuos);
}



void improveCluster(t_individuo populacao[], int alvo, int funcao, int total_individuos, double prob_mutacao, 
    int *comecar, t_individuo populacao_aux[], t_individuo melhores[], int maxClusters, t_individuo best)
{
    int numberOfChromosomes = 0;
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
            aux++;
            numberOfChromosomes++;
        }
    }

    //printf("%d\n", numberOfChromosomes);

    double percentageChromosomesToCrossover = 0.8;
    int totalChromosomesToCrossover = round(percentageChromosomesToCrossover * numberOfChromosomes);

    int vetor_aux[numberOfChromosomes];
    int vetor_pai[numberOfChromosomes];
    t_individuo novos_individuos[aux];

    // if (totalChromosomesToCrossover == 0)
    // {
    //     return;
    // }

    for (i = 0; i < numberOfChromosomes; ++i)
    {
        vetor_aux[i] = i;
    }

    shuffle(vetor_aux, numberOfChromosomes);

    for (i = 0; i < totalChromosomesToCrossover; ++i)
    {
        vetor_pai[i] = vetor_aux[i];
        cont2++;
    }

    for (i = cont2; i < numberOfChromosomes; ++i)
    {
        vetor_pai[i] = vetor_aux[i];
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
            op_uniformcrossover(&pai, &mae, funcao);
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
        totalChromosomesToCrossover++;
    }

    for (i = cont2; i < numberOfChromosomes; ++i)
    {
        novos_individuos[i] = pop_aux[vetor_pai[i]];
    }

    if (numberOfChromosomes == 2)
    {
        for (i = 0; i < numberOfChromosomes; ++i)
        {
            novos_individuos[i] = pop_aux[vetor_pai[i]];
        }
    }

    op_selecao_de_sobreviventes(populacao_aux, numberOfChromosomes, novos_individuos, *comecar);

    *comecar += numberOfChromosomes;

    if(numberOfChromosomes != 0)
    {
        encontra_melhor_individuo(novos_individuos, numberOfChromosomes, &melhor);
        melhores[alvo] = melhor;
    }
    else
    {
        melhores[alvo] = best;
        //gerar_individuo(&melhores[alvo], funcao);
    }
}

void improveIdeals(t_individuo melhores[], int maxClusters, int funcao, double prob_mutacao)
{
    t_individuo pai;
    t_individuo mae;
    t_individuo novos_individuos[maxClusters];
    int vetor_aux[maxClusters];
    int vetor_pai[maxClusters];
    int i, cont = 0, flag = 0;

    for (i = 0; i < maxClusters; ++i)
    {
        vetor_aux[i] = i;
    }

    shuffle(vetor_aux, maxClusters);

    for (i = 0; i < maxClusters; ++i)
    {
        vetor_pai[i] = vetor_aux[i];
    }

    if (maxClusters % 2 == 1 && maxClusters > 2)
    {
        maxClusters--;
        flag = 1;
    }

    if (maxClusters == 1)
    {
        novos_individuos[0] = melhores[0];
    }
    else
    {
        for (i = 0; i < maxClusters; ++i)
        {
            pai = melhores[vetor_pai[i]];
            mae = melhores[vetor_pai[i+1]];
            op_uniformcrossover(&pai, &mae, funcao);
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
        novos_individuos[cont] = melhores[vetor_pai[cont]];
        maxClusters++;
    }
    op_selecao_de_sobreviventes(melhores, maxClusters, novos_individuos, 0);
    //imprimir_populacao(melhores, maxClusters);
}

void generateNextPopulation(t_individuo populacao[], t_individuo melhores[], int total_individuos, int maxClusters, double prob_mutacao, int funcao, t_individuo best_melhores, int g, int geracoes)
{
    int i, j;
    t_individuo novos_individuos[total_individuos];
    int vetor_aux[total_individuos + maxClusters];
    int vetor_pai[total_individuos + maxClusters];
    t_individuo mutado;
    t_individuo new_individuo;
    int flag = 0;
    float l_inf = 0.0;
    float l_sup = 0.0;
    identificar_dominio(funcao,&l_inf,&l_sup);
    double u;

    for (i = 0; i < total_individuos + maxClusters; ++i)
    {
        vetor_aux[i] = i;
    }

    shuffle(vetor_aux, total_individuos + maxClusters);

    for (i = 0; i < total_individuos; ++i)
    {
        vetor_pai[i] = vetor_aux[i];
    }

    for (i = 0; i < total_individuos; ++i)
    {              
        if (g < geracoes*0.15)
        {
            u = obter_numero_uniforme_continuo(0.9, 1.0);
            //printf("Entrou1 %lf\n", u);
        }
        else if (g > geracoes*0.15 && g < geracoes*0.3)
        {
            u = obter_numero_uniforme_continuo(0.7, 1.0);
            //printf("Entrou2 %lf\n", u);
        }
        else if (g > geracoes*0.3 && g < geracoes*0.4)
        {
            u = obter_numero_uniforme_continuo(0.5, 1.0);
            //printf("Entrou3 %lf\n", u);
        }
        else
        {
            u = nextDouble();
            //printf("Entrou4 %lf\n", u);
        }        
        //u = nextDouble();
        if (u < 0.9999)
        {
            if (vetor_pai[i] >= total_individuos)
            {
                for (j = 0; j < NVARS; ++j)
                {
                    if (melhores[(total_individuos+maxClusters-1) - vetor_pai[i]].gene[j] < l_inf || melhores[(total_individuos+maxClusters-1) - vetor_pai[i]].gene[j] > l_sup || melhores[(total_individuos+maxClusters-1) - vetor_pai[i]].gene[j] == 0.0000000)
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
        else
        {
            gerar_individuo(&new_individuo, funcao);
            novos_individuos[i] = new_individuo;
        }     
    }
    encontra_melhor_individuo(melhores, maxClusters, &best_melhores);

    novos_individuos[total_individuos-2] = best_melhores;
    op_selecao_de_sobreviventes(populacao, total_individuos, novos_individuos, 0);
}

void escalonamento_linear(t_individuo pop[], double c, int total_individuos)
{
    t_individuo melhor;
    t_individuo pior;
    double alpha, beta, media;
    int i;
    encontra_melhor_individuo(pop, total_individuos, &melhor);
    encontra_pior_individuo(pop, total_individuos, &pior);
    media = encontra_media_populacao(pop, total_individuos);
    
    if(pior.fitness > (c*media-melhor.fitness)/(c-1.0))
    {
        alpha = (media*(c-1.0)/(melhor.fitness-media));
        beta = (media*(melhor.fitness-c*media)/(melhor.fitness-media));
    }
    else
    {
        alpha = media/(media-pior.fitness);
        beta = -pior.fitness*media/(media-pior.fitness);
    }
        
    for (i = 0; i < total_individuos; ++i)
    {
        pop[i].fitness = pop[i].fitness * alpha + beta;
        if (pop[i].fitness < 0 || pop[i].fitness == 0.000)
        {

            pop[i].fitness = 0.001;
        }
    }

}

double euclidean_distance(t_individuo ind1, t_individuo ind2)
{
    double euclidean = 0.0;
    int i;
    for (i = 0; i < NVARS; ++i)
    {
        euclidean += pow((ind1.gene[i] - ind2.gene[i]),2);
    }
    return sqrt(euclidean);
}

t_individuo encontra_melhor(t_individuo ind1, t_individuo ind2)
{
    if (ind1.fitness < ind2.fitness)
    {
        return ind1;
    }
    else
    {
        return ind2;
    }
}

double encontra_menor_distancia(double *dist, int CROWDING_SIZE)
{
    int i;
    double menor = 9999;
    for (i = 0; i < CROWDING_SIZE; ++i)
    {
        if (dist[i] < menor)
        {
            menor = dist[i];
        }
    }
    return menor;
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


    double **vet_melhores = (double **)malloc(RUNS * sizeof(double*));
        for(int i = 0; i < RUNS; i++) vet_melhores[i] = (double *)malloc(geracoes * sizeof(double));
    double **vet_diversidade = (double **)malloc(RUNS * sizeof(double*));
        for(int i = 0; i < RUNS; i++) vet_diversidade[i] = (double *)malloc(geracoes * sizeof(double));

    struct timeval timevalA;
    struct timeval timevalB;

    FILE *fpMedia;
    FILE *fpDiversidade;
    FILE *fpNumeroCluster;
    FILE *best_geracoes;

    char buf[0x100];
    snprintf(buf, sizeof(buf), "Dados_SINGLELINK/GERACOESSINGLE_FUNCAO%d_%dDIMENSOES.txt", funcao, NVARS);
    char buf1[0x100];
    snprintf(buf1, sizeof(buf1), "Dados_SINGLELINK/DIVERSITYSINGLE_FUNCAO%d_%dDIMENSOES.txt", funcao, NVARS);
    char buf2[0x100];
    snprintf(buf2, sizeof(buf2), "Dados_SINGLELINK/CLUSTER_FUNCAO%d_%dDIMENSOES.txt", funcao, NVARS);
    char buf3[0x100];
    snprintf(buf3, sizeof(buf3), "Dados_SINGLELINK/BEST_GERACOES_FUNCAO%d_%dDIMENSOES.txt", funcao, NVARS);

    fpMedia = fopen(buf, "w+");
    fpDiversidade = fopen(buf1, "w+");
    best_geracoes = fopen(buf3, "w+");
    fpNumeroCluster = fopen(buf2, "w+");
    gettimeofday(&timevalA,NULL);
    for (run = 0; run < RUNS; ++run)
    {
        double aux_inicio = 0.3*geracoes;
        double aux_final = 0.8*geracoes;
        int c = 1.2;
        double firstnessIndex = 2;
        double killIndex = firstnessIndex + 1;
        int maxClusters = 0;
        double prob_mutacao = 0.2;
        double minMutationRate = 0.02;
        m_nmdf = 0.0;        

        t_individuo populacao[total_individuos];
        
        gerar_populacao_inicial(populacao, total_individuos, funcao);

        int g = 0; //contador de geracoes
        int i;
        double u = 0.001;
        for(g = 0; g < geracoes; g++){
            double clusterID[total_individuos];
            double distanceFromCenter[total_individuos];
            int comecar = 0;

            t_individuo populacao_aux[total_individuos];
            t_individuo best;
            t_individuo best_after;
            t_individuo best_melhores;
            int vetor_aux[total_individuos];
            int vetor_pai[total_individuos];
            t_individuo best_ag;
            best_ag.fitness = 99999999999999;
            t_individuo mae;
            t_individuo pai;
            t_individuo *novos_individuos = (t_individuo*) malloc(total_individuos*sizeof(t_individuo));


            encontra_melhor_individuo(populacao, total_individuos, &best);

            if(g < geracoes*0.9)
            {
                clusterAnalysis(populacao, total_individuos, g, u);

                if (u > 0.9)
                {
                    u = 0.9;
                }
                else if (g < geracoes*0.8)
                {
                    u += (0.9/((double)geracoes+geracoes*0.8));
                }
                else
                {
                    u += (0.9/((double)geracoes));
                }

                maxClusters = g_habitatsSize;

                t_individuo melhores[maxClusters];

                for (i = 0; i < maxClusters; ++i)
                {
                    improveCluster(populacao, i, funcao, total_individuos, prob_mutacao, &comecar, populacao_aux, melhores, maxClusters, best);
                    
                    //imprimir_populacao(melhores, i);
                }
                improveIdeals(melhores, maxClusters, funcao, prob_mutacao);   

                memcpy(populacao, populacao_aux, sizeof(t_individuo)*total_individuos);

                generateNextPopulation(populacao, melhores, total_individuos, maxClusters, prob_mutacao, funcao, best_melhores, g, geracoes);

                encontra_melhor_individuo(populacao, total_individuos, &best_after);
            }
            else
            {
                for (i = 0; i < total_individuos; ++i)
                {
                    vetor_aux[i] = i;
                }

                shuffle(vetor_aux, total_individuos);

                for (i = 0; i < total_individuos; ++i)
                {
                    vetor_pai[i] = vetor_aux[i];
                }

                for (i = 0; i < total_individuos; ++i)
                {
                    op_selecao_de_pais(populacao, total_individuos, &pai, &mae);
                    double chance_crossover = nextDouble();
                    op_uniformcrossover(&pai, &mae, funcao);                
                    op_mutacao(&pai,0.02,funcao);
                    op_mutacao(&mae,0.02,funcao);
                    novos_individuos[i] = pai;
                    i++;
                    novos_individuos[i] = mae;
                }                
                op_selecao_de_sobreviventes(populacao, total_individuos, novos_individuos, 0);
                encontra_melhor_individuo(populacao, total_individuos, &best_ag);
            }          

            if (best_after.fitness < best.fitness)
            {
                populacao[total_individuos-1] = best_after;
            }
            else if(best_ag.fitness < best.fitness && best_ag.fitness < best_after.fitness)
            {
                populacao[total_individuos-1] = best_ag;
            }
            else
            {
                populacao[total_individuos-1] = best;
            }

            encontra_melhor_individuo(populacao, total_individuos, &best);

            double diversidade = diversity_population(populacao, total_individuos);
                
            prob_mutacao -= 0.001;

            if (prob_mutacao <= minMutationRate)
            {
                prob_mutacao = minMutationRate;
            }


            vet_melhores[run][g] = best.fitness;
            vet_diversidade[run][g] = diversidade;
            if (run == 0)
            {
                fprintf(fpNumeroCluster, "%d\n", g_habitatsSize);
            }            
        }
        printf("%d\n", run);
        //imprimir_populacao(populacao, total_individuos);
    }
    gettimeofday(&timevalB,NULL);

    printf("\n\ntempo de execucao: %lf\n",  (timevalB.tv_sec-timevalA.tv_sec+(timevalB.tv_usec-timevalA.tv_usec)/(double)1000000)/10.0);

    double *mediabest = (double *)malloc(geracoes * sizeof(double));
    double *mediaDiversity = (double *)malloc(geracoes * sizeof(double));


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
        fprintf(fpMedia, "%.10f\n", mediabest[i]);
        fprintf(fpDiversidade, "%.10f\n", mediaDiversity[i]);
    }

    double sum = 0, sum_squares = 0;

    for (i = 0; i < RUNS; ++i)
    {
        sum_squares += vet_melhores[i][geracoes-1] * vet_melhores[i][geracoes-1];
        sum += vet_melhores[i][geracoes-1];        
        fprintf(best_geracoes, "%lf,\n", vet_melhores[i][geracoes-1]);
    }

    double mean = (double)sum / RUNS;
    double variance = (double)sum_squares / RUNS - (mean * mean);
    double std_dev = sqrt(variance);

    printf("Mean: %.10f +/- %f", mean, std_dev);

    fclose(best_geracoes);
    fclose(fpMedia);
    fclose(fpDiversidade);
    fclose(fpNumeroCluster);
}