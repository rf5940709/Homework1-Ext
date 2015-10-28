#include <stdint.h>
#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

double compute_pi_Baseline(size_t dt)
{
    double pi = 0.0;
    double delta = 1.0 / dt;
    for (size_t i = 0; i < dt; i++) {
        double x = (double) i / dt;
        pi += delta / (1.0 + x * x);
    }
    return pi * 4.0;
}
double compute_pi_Baseline_avx(size_t dt)
{
    double pi = 0.0;
    double delta = 1.0 / dt;
    register __m256d ymm0, ymm1, ymm2, ymm3, ymm4;
    ymm0 = _mm256_set1_pd(1.0);
    ymm1 = _mm256_set1_pd(delta);
    ymm2 = _mm256_set_pd(delta * 3, delta * 2, delta * 1, 0.0);
    ymm4 = _mm256_setzero_pd();

    for (int i = 0; i <= dt - 4; i += 4) {
        ymm3 = _mm256_set1_pd(i * delta);//x = (double) i / dt;
        ymm3 = _mm256_add_pd(ymm3, ymm2);
        ymm3 = _mm256_mul_pd(ymm3, ymm3);//x=x*x
        ymm3 = _mm256_add_pd(ymm0, ymm3);//x=1.0+x
        ymm3 = _mm256_div_pd(ymm1, ymm3);//x=1/x or-1/x
        ymm4 = _mm256_add_pd(ymm4, ymm3);//pi=pi+x
    }
    double tmp[4] __attribute__((aligned(32)));
    _mm256_store_pd(tmp, ymm4);
    pi += tmp[0] + tmp[1] + tmp[2] + tmp[3];

    return pi * 4.0;
}
double compute_pi_Gregory(size_t dt)
{
    double pi = 1;
    for (size_t i = 1; i < dt; i++) {
        if(i%2==1)
            pi -= (1/(double)((i*2)+1));
        else
            pi += (1/(double)((i*2)+1));
    }
    return pi * 4.0;
}
double compute_pi_Euler(size_t dt)
{
    double pi = 0;
    for (size_t i = 1; i < dt; i++)
        pi += (1/((double)i*(double)i));
    return sqrt(pi*6);
}
int main(int argc, char *argv[])
{
    size_t n = atoi(argv[1]);
    clock_t begin, end;
    double time, value;

    begin = clock();
    value = compute_pi_Baseline(n * 1024*1024);
    end = clock();
    time = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("compult_pi loop : %zu MB, vaile : %lf \n", n, value);
    printf("execution time of compute_pi_Baseline()      : %lf sec\n", time);

    begin = clock();
    value = compute_pi_Baseline_avx(n * 1024*1024);
    end = clock();
    time = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("compult_pi loop : %zu MB, vaile : %lf \n", n, value);
    printf("execution time of compute_pi_Baseline_avx()  : %lf sec\n", time);

    begin = clock();
    value = compute_pi_Gregory(n * 1024*1024);
    end = clock();
    time = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("compult_pi loop : %zu MB, vaile : %lf \n", n, value);
    printf("execution time of compute_pi_Gregory()       : %lf sec\n", time);

    begin = clock();
    value = compute_pi_Euler(n * 1024*1024);
    end = clock();
    time = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("compult_pi loop : %zu MB, vaile : %lf \n", n, value);
    printf("execution time of compute_pi_Euler()         : %lf sec\n", time);

    return 0;
}
