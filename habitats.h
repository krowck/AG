#include <math.h>
#include <stdio.h>

#include "gerador_numeros.h"

/**************************************************************************
        ALGORITMO DE CLUSTERIZAÇÃO - SINGLE LINK
**************************************************************************/

#define MINVAL 0.01 
#define MAXVAL 0.99 

struct habitat {
    int h_ind_count;
    int *h_ind;
} *H;

unsigned short int g_habitatsSize = 0;
Node * g_tree;

void initHabitats(int total_individuos)
{
    unsigned short int a = 0;

    H = (struct habitat*) malloc(total_individuos * sizeof(struct habitat));    
    for (a = 0; a < total_individuos; a++)
    {
        H[a].h_ind_count = 0; 
        H[a].h_ind = (int*)malloc (sizeof(int)*total_individuos );
    }
}


void singlelink(int total_individuos, int dimensions, double** distances)
{
    double maxrange;
    double minrange;
    unsigned short int i = 0, j=0;

    //for (i = 0; i < total_individuos; ++i)
    //{
    //    for (j = 0; j < total_individuos; ++j)
    //     {
    //         printf("distances[%d][%d] = %f\n",i, j, distances[i][j]);
    //     }
    //}

    g_tree = treecluster(total_individuos, total_individuos, 0, 0, 0, 0, 'e', 's', distances);

    if(!g_tree)
    {
        printf("Erro ao gerar.\n");
        return;
    }

    maxrange = -1;
    for (i = 0; i < total_individuos-1;i++)
    {
        if (maxrange < g_tree[i].distance)
        {
            maxrange = g_tree[i].distance;
        }
    }
    minrange = maxrange;
    for (i = 0; i < total_individuos-1; i++)
    {
        if (minrange > g_tree[i].distance)
        {
            minrange = g_tree[i].distance;
        }
    }

    for (i = 0 ;i < total_individuos-1;i++)
    {
            g_tree[i].distance=(MAXVAL-MINVAL)/(double)(maxrange-minrange)*g_tree[i].distance-(MAXVAL-MINVAL)/(double)(maxrange-minrange)*minrange+MINVAL;
            //printf("%lf\n", g_tree[i].distance);
    }
    //printf("\n\n\n");   
}

void printDendogram(int total_individuos)
{
    unsigned short int i=0;
    FILE *dendogram;
    dendogram = fopen("dendogram.txt", "w+");
    for(i=0; i< total_individuos-1; i++)
    {
        if (g_tree[i].left >= 0 && g_tree[i].right >= 0)
            fprintf(dendogram, "%d  %d          %f\n", g_tree[i].left+1, g_tree[i].right+1, g_tree[i].distance);
            //printf("%d  %d          %f\n", g_tree[i].left+1, g_tree[i].right+1, g_tree[i].distance);
        if (g_tree[i].left >= 0 && g_tree[i].right < 0)
            fprintf(dendogram, "%d  %d          %f\n", g_tree[i].left+1, (-1*g_tree[i].right)+total_individuos, g_tree[i].distance);
            //printf("%d  %d          %f\n", g_tree[i].left+1, (-1*g_tree[i].right)+total_individuos, g_tree[i].distance);
        if (g_tree[i].left < 0 && g_tree[i].right >= 0)
            fprintf(dendogram, "%d  %d          %f\n", (-1*g_tree[i].left)+total_individuos, g_tree[i].right+1, g_tree[i].distance);
            //printf("%d  %d          %f\n", (-1*g_tree[i].left)+total_individuos, g_tree[i].right+1, g_tree[i].distance);
        if (g_tree[i].left < 0 && g_tree[i].right < 0)
            fprintf(dendogram, "%d  %d          %f\n", (-1*g_tree[i].left)+total_individuos,(-1*g_tree[i].right)+total_individuos, g_tree[i].distance);
            //printf("%d  %d          %f\n", (-1*g_tree[i].left)+total_individuos,(-1*g_tree[i].right)+total_individuos, g_tree[i].distance);
    }
    fclose(dendogram);
}


void buildHabitats(int total_individuos, double** distances)
{
    unsigned short int i = 0;
    unsigned short int j = 0;
    unsigned short int k = 0;
    unsigned short int sum_tabu = 0;
    int cur_node = total_individuos - 1;
    int cur_habitat = 0;
    unsigned short int cur_node_index = 0;
    double sum_left = 0.0;
    double sum_right = 0.0;
    double dist_aux = 0.0;

    g_habitatsSize = 0;

    for (i = 0; i < total_individuos; i++)
    {
        H[i].h_ind_count = 0; 
    }

    while (sum_tabu < total_individuos-1)
    {
        if ( obter_numero_uniforme() >= g_tree[cur_node-1].distance) //link those two itens
        {
            if (cur_node == total_individuos-1) //first iteration.
            {
                H[cur_habitat].h_ind[ H[cur_habitat].h_ind_count ] = g_tree[cur_node-1].left;
                H[cur_habitat].h_ind[ H[cur_habitat].h_ind_count+1 ] = g_tree[cur_node-1].right;
                H[cur_habitat].h_ind_count = 2;
                sum_tabu++;
            }
            else
            {
                if (cur_node == 1)
                {
                    H[cur_habitat].h_ind[ cur_node_index ] = g_tree[0].left;
                    H[cur_habitat].h_ind[ H[cur_habitat].h_ind_count ] = g_tree[0].right;
                    H[cur_habitat].h_ind_count += 1;
                    sum_tabu++;
                }
                else
                {
                    H[cur_habitat].h_ind[ cur_node_index ] = g_tree[cur_node-1].left;
                    H[cur_habitat].h_ind[ H[cur_habitat].h_ind_count ] = g_tree[cur_node-1].right;
                    H[cur_habitat].h_ind_count += 1;
                    sum_tabu++;
                }
            }
        }
        else
        {
            if (cur_node == total_individuos-1)
            {
                H[cur_habitat].h_ind[ H[cur_habitat].h_ind_count ] = g_tree[cur_node-1].left;
                g_habitatsSize++;
                H[g_habitatsSize].h_ind[ H[g_habitatsSize].h_ind_count ] = g_tree[cur_node-1].right;
                H[cur_habitat].h_ind_count += 1;
                H[g_habitatsSize].h_ind_count += 1;
                sum_tabu++;
            }
            else
            {
                if (H[cur_habitat].h_ind_count == 1 && H[cur_habitat].h_ind[0] < 0)
                {
                    H[cur_habitat].h_ind[ cur_node_index ] = g_tree[cur_node-1].left;
                    g_habitatsSize++;
                    H[g_habitatsSize].h_ind[ H[g_habitatsSize].h_ind_count ] = g_tree[cur_node-1].right;
                    H[g_habitatsSize].h_ind_count += 1;
                    sum_tabu++;
                }
                if (g_tree[cur_node-1].left >= 0 && g_tree[cur_node-1].right >= 0 && H[cur_habitat].h_ind_count > 1) 
                {
                    sum_left = sum_right = 0.0;
                    for (i = 0; i < H[cur_habitat].h_ind_count; i++)
                    {
                       if (H[cur_habitat].h_ind[ i ] >= 0)
                       {
                        sum_left  += distances[ H[cur_habitat].h_ind[ i ] ][ g_tree[cur_node-1].left ];
                        sum_right += distances[ H[cur_habitat].h_ind[ i ] ][ g_tree[cur_node-1].right ];
                       }
                    }
                    if (sum_left <= sum_right) 
                    {
                        H[cur_habitat].h_ind[ cur_node_index ] = g_tree[cur_node-1].left;
                        g_habitatsSize++;
                        H[g_habitatsSize].h_ind[ H[g_habitatsSize].h_ind_count ] = g_tree[cur_node-1].right;
                    }
                    else
                    {
                        H[cur_habitat].h_ind[ cur_node_index ] = g_tree[cur_node-1].right;
                        g_habitatsSize++;
                        H[g_habitatsSize].h_ind[ H[g_habitatsSize].h_ind_count ] = g_tree[cur_node-1].left;
                    }
                    H[g_habitatsSize].h_ind_count += 1;
                    sum_tabu++;
                }
                else if (g_tree[cur_node-1].left < 0 && g_tree[cur_node-1].right < 0 && H[cur_habitat].h_ind_count > 1)
                {
                    if ( g_tree[cur_node-1].left < g_tree[cur_node-1].right )
                    {
                        H[cur_habitat].h_ind[ cur_node_index ] = g_tree[cur_node-1].left;
                        g_habitatsSize++;
                        H[g_habitatsSize].h_ind[ H[g_habitatsSize].h_ind_count ] = g_tree[cur_node-1].right;
                    }
                    else
                    {
                        H[cur_habitat].h_ind[ cur_node_index ] = g_tree[cur_node-1].right;
                        g_habitatsSize++;
                        H[g_habitatsSize].h_ind[ H[g_habitatsSize].h_ind_count ] = g_tree[cur_node-1].left;
                    }                   
                    H[g_habitatsSize].h_ind_count += 1;
                    sum_tabu++;
                }
                else if (H[cur_habitat].h_ind_count > 1) 
                {
                    if (g_tree[cur_node-1].left >= 0) 
                    {
                        H[cur_habitat].h_ind[ cur_node_index ] = g_tree[cur_node-1].left;
                        g_habitatsSize++;
                        H[g_habitatsSize].h_ind[ H[g_habitatsSize].h_ind_count ] = g_tree[cur_node-1].right;
                    }else
                    {
                        H[cur_habitat].h_ind[ cur_node_index ] = g_tree[cur_node-1].right;
                        g_habitatsSize++;
                        H[g_habitatsSize].h_ind[ H[g_habitatsSize].h_ind_count ] = g_tree[cur_node-1].left;
                    }
                    H[g_habitatsSize].h_ind_count += 1;
                    sum_tabu++;
                }
            }
        }
        i = 0;
        cur_habitat = -1;
        j = 0;
        while (i <= g_habitatsSize && cur_habitat < 0)
        {
            while (j < H[i].h_ind_count && cur_habitat < 0)
            {
                if (H[i].h_ind[j] < 0)
                {
                    cur_habitat = i;
                }
                else j++;               
            }
            i++;
            j = 0;
        }
        if (cur_habitat >= 0)
        {
            i = 0;
            while (H[cur_habitat].h_ind[i] >= 0)
                i++;
            cur_node = -1*H[cur_habitat].h_ind[i];
            cur_node_index = i;
        }
    }
    g_habitatsSize++;
    
    

}
/**************************************************************************
        FIM ALGORITMO DE CLUSTERIZAÇÃO - SINGLE LINK
**************************************************************************/