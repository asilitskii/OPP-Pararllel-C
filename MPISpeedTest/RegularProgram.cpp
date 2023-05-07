#include <stdio.h>
#include <malloc.h>
#include <time.h>

#define N 200000

double operation(int* a, int* b) {
    double s = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            s += (double)(a[i]) * (double)(b[j]);
    return s;
}

void fillArray(int* array, int value) {
    for (int i = 0; i < N; ++i)
        array[i] = value;
}

int main(int argc, char** argv){
    double timeBefore, timeAfter;
    int* a = (int*)malloc(N * sizeof(int));
    int* b = (int*)malloc(N * sizeof(int));
    fillArray(a, 1);
    fillArray(b, 1);

    timeBefore = clock();
    double result = operation(a, b);
    timeAfter = clock();

    printf("Result: %lf\n Result time: %lf", result, (double)(timeAfter - timeBefore) / CLOCKS_PER_SEC);


    free(a);
    free(b);
    return 0;
}

