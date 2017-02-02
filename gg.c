#include <math.h>
#include <stdio.h>

int main(void)
{
    FILE *in = fopen("mediaGeracoes_ZAKHAROV20D.txt", "r");
    
    if (in != NULL) {
        double sum = 0, sum_squares = 0, n = 0;
        double val;
        
        while (fscanf(in, "%lf", &val) == 1) {
            sum += val;
            sum_squares += val * val;
            ++n;
        }
        printf("%f\n", n);
        fclose(in);
        
        if (n > 0) {
            double mean = (double)sum / n;
            double variance = (double)sum_squares / n - (mean * mean);
            double std_dev = sqrt(variance);
            
            printf("Mean: %f\nStdDev: %f", mean, std_dev);
        }
    }
    
    return 0;
}