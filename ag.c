#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <windows.h>

#include "funcoes_benchmark.h"
#include "gerador_numeros.h"


#define NVARS 6
#define TAM_TORNEIO 5

/*
 * A estrutura (struct) "t_individuo" representa um unico individuo da populacao
 */
typedef struct {
    double gene[NVARS];
    double fitness;
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
    int i;
    switch (funcao)
    {
        case 1:
            fitness = rastrigin(cromossomo);
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
    printf("fitness = %f\n\n",individuo.fitness);    
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
void op_mutacao(t_individuo *filho, float prob_mutacao, int funcao){
    int i;
    float u = obter_numero_uniforme();
    
    if (u <= prob_mutacao){
    
        float l_inf = 0.0;
        float l_sup = 0.0;
        identificar_dominio(funcao,&l_inf,&l_sup);

        for (i = 0; i < NVARS; ++i){        
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
    filho->fitness = obter_fitness(funcao, filho->gene);      
    }
    
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

void op_crossover(t_individuo *pai, t_individuo *mae,  int total_individuos, int funcao){
    int point = obter_numero_uniforme_discreto(0, NVARS-1);

    int i;
    double t;

    //printf("%d\n", point);

    for (i = 0; i < point; i++)
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

/*
apenas atualizar o vetor de população com os novos individuos gerador a partir de crossover e mutação
*/
void op_selecao_de_sobreviventes(t_individuo populacao[], int total_individuos, t_individuo novos_individuos[], int descarte){
    int i = 0;
    int j = total_individuos-1;
    
    for(i = 0; i < descarte; i++){
        populacao[i] = novos_individuos[i];
    }
}

double diversity_population(t_individuo pop[], int tamPopulation)
{
    double m_nmdf = 0;
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
/*
 * Evolucao da Populacao
 * No procedimento "executar" que e realizada a Evolucao da Populacao.
 
 * Os parametros de entrada (definidos pelo usuario) deste procedimento sao:
 *  - o codigo (ou ID) da funcao a ser otimizada
 *  - o tamanho da populacao (ou total de individuos "total_individuos" da populacao)
 *  - por quantas "geracoes" a populacao inicial sera evoluida
 *  - a probabilidade de mutacao (prob_mutacao)
 */
void executar(int funcao, int total_individuos, int geracoes, float prob_mutacao){
    srand((unsigned)time(NULL)); 
    /*
     * A Populacao e representada como um vetor de "t_individuo", cujo o tamanho e "total_individuos" (definido previamente pelo usuario).
     */
    t_individuo populacao[total_individuos]; 
    
    gerar_populacao_inicial(populacao, total_individuos, funcao);
    printf("POPULAÇAO INICIAL\n\n");
    imprimir_populacao(populacao, total_individuos);
    int g = 0; //contador de geracoes
    int j = 0;
    printf("\n#\tx_1\t\tx_2\t\tf(x_1, x_2)                             diversidade\n\n"); //Saida de Dados
    for(g = 0; g < geracoes; g++){
        double diversidade = diversity_population(populacao, total_individuos);
        //imprimir_populacao(populacao, total_individuos);
        t_individuo novos_individuos[total_individuos]; //vetor de novos individuos
        
        int i;

        for(i = 0; i < total_individuos; i++){
        
            t_individuo pai;
            t_individuo mae;

            op_selecao_de_pais(populacao, total_individuos, &pai, &mae);

            op_crossover(&pai, &mae, total_individuos, funcao);

            op_mutacao(&pai,prob_mutacao,funcao);

            op_mutacao(&mae,prob_mutacao,funcao);
            novos_individuos[i] = pai;
            i++;
            novos_individuos[i] = mae;
        }
        mergeSort(populacao, total_individuos);

        op_selecao_de_sobreviventes(populacao,total_individuos,novos_individuos,total_individuos);
        
        mergeSort(populacao, total_individuos);
        
        for (j = 0; j < NVARS; ++j)
        {
            printf("%f ", populacao[0].gene[j]);
        }
        // Saida de Dados
        printf("%d\t%f  ",g,populacao[0].fitness); //saida de dados
        printf("  %f\n", diversidade);
        
    }

}

