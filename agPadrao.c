/*********************************************************************************************************************
*                                   Algoritmo desenvolvido para o TCC                                                *
*  para compilar:  gcc agPadrao.c -lm -o ag principal.c funcoes_benchmark.c gerador_numeros.c -O3              *
*                                                                                                                    *
*                                                                                                                    *
*********************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//#include <windows.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>

#include "funcoes_benchmark.h"
#include "gerador_numeros.h"

#define TAM_TORNEIO 5
#define RUNS 5

double m_nmdf = 0;

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
    printf("\nfitness = %f\n\n",individuo.fitness);    
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
    
}

void gerar_individuo(t_individuo individuo, int funcao)
{
    float l_inf = 0.0;
    float l_sup = 0.0;
    identificar_dominio(funcao,&l_inf,&l_sup);

    int j;
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
    int point = obter_numero_uniforme_discreto(1, NVARS-2);

    int i;
    double t;

    //if (isnan(mae->fitness) || isnan(pai->fitness))
    //{
    //    printf("%lf  %lf\n", mae->fitness, pai->fitness);
    //}

    for (i = 0; i <= point; i++)
    {
        t = pai->gene[i];
        pai->gene[i] = mae->gene[i];
        mae->gene[i] = t;
    }
    pai->fitness = obter_fitness(funcao, pai->gene);
    mae->fitness = obter_fitness(funcao, mae->gene);

    //if (isnan(mae->fitness) || isnan(pai->fitness))
    //{
    //    printf("%lf  %lf\n", mae->fitness, pai->fitness);
    //}

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
 * Evolucao da Populacao
 * No procedimento "executar" que e realizada a Evolucao da Populacao.
 
 * Os parametros de entrada (definidos pelo usuario) deste procedimento sao:
 *  - o codigo (ou ID) da funcao a ser otimizada
 *  - o tamanho da populacao (ou total de individuos "total_individuos" da populacao)
 *  - por quantas "geracoes" a populacao inicial sera evoluida
 */
void executar(int funcao, int total_individuos, int geracoes, double prob_mutacao){
    srand((unsigned)time(NULL)); 
    //double vet_melhores[RUNS][geracoes];
    //double vet_diversidade[RUNS][geracoes];

    double **vet_melhores = (double **)malloc(RUNS * sizeof(double*));
        for(int i = 0; i < RUNS; i++) vet_melhores[i] = (double *)malloc(geracoes * sizeof(double));
    double **vet_diversidade = (double **)malloc(RUNS * sizeof(double*));
        for(int i = 0; i < RUNS; i++) vet_diversidade[i] = (double *)malloc(geracoes * sizeof(double));
    int run;
    int i;
    int j;

    struct timeval timevalA;
    struct timeval timevalB;

    FILE *fpMedia;
    FILE *fpDiversidade;
    FILE *fp;

    char buf[0x100];
    snprintf(buf, sizeof(buf), "DADOS_AGPADRAO/GERACOESPADRAO_FUNCAO%d_%dDIMENSOES.txt", funcao, NVARS);
    char buf1[0x100];
    snprintf(buf1, sizeof(buf1), "DADOS_AGPADRAO/DIVERSITYPADRAO_FUNCAO%d_%dDIMENSOES.txt", funcao, NVARS);

    fpMedia = fopen(buf, "w+");
    fpDiversidade = fopen(buf1, "w+");

    fp = fopen("output_AGPADRAO.txt", "w+");
    gettimeofday(&timevalA,NULL);
    for (run = 0; run < RUNS; ++run)
    {

        int maxClusters = 2;
        int flagCluster = 0;
        double prob_mutacao = 0.2;
        double minMutationRate = 0.002;
        m_nmdf = 0.0;        
        
        /*
         * A Populacao e representada como um vetor de "t_individuo", cujo o tamanho e "total_individuos" (definido previamente pelo usuario).
         */
        t_individuo populacao[total_individuos];
        
        gerar_populacao_inicial(populacao, total_individuos, funcao);

        int g = 0; //contador de geracoes
        int i = 0;
        for(g = 0; g < geracoes; g++){

            t_individuo populacao_aux[total_individuos];
            t_individuo best;
            t_individuo mae;
            t_individuo pai;
            int vetor_aux[total_individuos];
            int vetor_pai[total_individuos];
            t_individuo novos_individuos[total_individuos];

            encontra_melhor_individuo(populacao, total_individuos, &best);

            for (i = 0; i < total_individuos; ++i)
            {
                vetor_aux[i] = i;
            }

            shuffle(vetor_aux, total_individuos);

            for (i = 0; i < total_individuos; ++i)
            {
                vetor_pai[i] = vetor_aux[i];
            }

            for (int i = 0; i < total_individuos; ++i)
            {
                pai = populacao[vetor_pai[i]];
                mae = populacao[vetor_pai[i+1]];
                op_crossover(&pai, &mae, funcao);
                op_mutacao(&pai,prob_mutacao,funcao);
                op_mutacao(&mae,prob_mutacao,funcao);
                novos_individuos[i] = pai;
                i++;
                novos_individuos[i] = mae;
            }
            op_selecao_de_sobreviventes(populacao, total_individuos, novos_individuos, 0);

            double diversidade = diversity_population(populacao, total_individuos);

            populacao[total_individuos-1] = best;

            encontra_melhor_individuo(populacao, total_individuos, &best);

            fprintf(fp, "%d %f %f\n", g, best.fitness, diversidade);

            prob_mutacao -= 0.001;

            if (prob_mutacao <= minMutationRate)
            {
                prob_mutacao = minMutationRate;
            }
            vet_melhores[run][g] = best.fitness;
            vet_diversidade[run][g] = diversidade;
        }

        {
            /* code */
        }
        printf("%d\n", run);
        fprintf(fp, "-------------------------------");    
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

    FILE *in = fopen(buf, "r");
    
    if (in != NULL) {
        double sum = 0, sum_squares = 0, n = 0;
        double val;
        
        while (fscanf(in, "%lf", &val) == 1) {
            sum += val;
            sum_squares += val * val;
            ++n;
        }
        fclose(in);
        
        if (n > 0) {
            double mean = (double)sum / n;
            double variance = (double)sum_squares / n - (mean * mean);
            double std_dev = sqrt(variance);
            
            printf("Mean: %f\nStdDev: %f", mean, std_dev);
        }
    }

}

