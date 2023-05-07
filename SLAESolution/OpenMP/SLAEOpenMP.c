#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <limits.h>
#include <omp.h>

#define MATRIX_SIZE 8000
#define EPSILON 0.000001
#define MAX_ITERATIONS 10000

double* initMatrix(int height, int width) {
    double* matrix = (double*)calloc(height * width, sizeof(double*));
    return matrix;
}

double* initVector(int size) {
    double* vector = (double*)calloc(size, sizeof(double));
    return vector;
}

void fillVectorRandom(double* array, int size) {
    srand(size);
    for (int i = 0; i < size; ++i)
        array[i] = (double)rand() / RAND_MAX * (double)rand();
}

void fillMatrixRandom(double* martix, int size) {
    srand(size);
    for (int i = 0; i < size; ++i)
        martix[i] = (double)rand() / RAND_MAX * 2.0 - 1.0;
}

void addToDiagonal(double* matrix, int size, double value) {
    for (int i = 0; i < size; ++i)
        matrix[i * size + i] += value;
}

void fillVector(double* array, int size, double value) {
    for (int i = 0; i < size; ++i)
        array[i] = value;
}

void fillMatrix(double* matrix, int height, int width, double value) {
    for (int i = 0; i < height; ++i)
        for (int j = 0; j < width; ++j)
            matrix[i * height + j] = value;
}

void freeMatrix(double* matrix) {
    free(matrix);
}

void freeVector(double* vector) {
    free(vector);
}

double* useIterativeSolving(double* A, double* b, int size) {
    double* xCurrent = initVector(size);
    fillVector(xCurrent, size, 0);
    double* resultOfSubtitution = initVector(size);
    double* approximatoinVector = initVector(size);
    double tau = 0;
    double tauNumerator = 0;
    double tauDenominator = 0;
    double currentCondition = 0;
    int currentIteration = 0;
    double exitCondition = 0;

#pragma omp parallel shared(currentIteration, currentCondition, exitCondition)
    {
        //calculating exit condition
#pragma omp for reduction(+:exitCondition) schedule(runtime)
        for (int i = 0; i < size; ++i)
            exitCondition += b[i] * b[i];
#pragma omp single
        exitCondition *= EPSILON * EPSILON;

        //calculating result of subtitution of current solvation vector to system
#pragma omp for schedule(runtime)
        for (int i = 0; i < size; ++i) {
            double dotProduct = 0;
            for (int j = 0; j < size; j++) {
                dotProduct += A[i * size + j] * xCurrent[j];
            }
            resultOfSubtitution[i] = dotProduct - b[i];
        }

        while (currentIteration < MAX_ITERATIONS) {
            //approximatoinVector = A * resultOfSubtitution (Ay)
#pragma omp for schedule(runtime)
            for (int i = 0; i < size; ++i) {
                double dotProduct = 0;
                for (int j = 0; j < size; j++) {
                    dotProduct += A[i * size + j] * resultOfSubtitution[j];
                }
                approximatoinVector[i] = dotProduct;
            }
            //calculate tau
#pragma omp single
            {
                tau = 0;
                tauNumerator = 0;
                tauDenominator = 0;
            }
#pragma omp for reduction(+:tauNumerator, tauDenominator) schedule(runtime)
            for (int i = 0; i < size; ++i) {
                tauNumerator += resultOfSubtitution[i] * approximatoinVector[i];
                tauDenominator += approximatoinVector[i] * approximatoinVector[i];
            }
            tau = tauNumerator / tauDenominator;
            //calculating next solution
#pragma omp for schedule(runtime)
            for (int i = 0; i < size; ++i)
                xCurrent[i] -= tau * resultOfSubtitution[i];
            //calculating next result of subtitution to system
#pragma omp for schedule(runtime)
            for (int i = 0; i < size; ++i) {
                double tmp = 0;
                for (int j = 0; j < size; j++) {
                    tmp += A[i * size + j] * xCurrent[j];
                }
                resultOfSubtitution[i] = tmp - b[i];
            }

#pragma omp single
            {
                ++currentIteration;
                currentCondition = 0;
            }
            //calculating difference between current x and solution and if it fits - finalise
#pragma omp for reduction(+:currentCondition) schedule(runtime)
            for (int i = 0; i < size; ++i)
                currentCondition += resultOfSubtitution[i] * resultOfSubtitution[i];
            if (currentCondition < exitCondition) break;
        }
    }
    //printf("Number of iterations: %d\n", currentIteration);
    free(resultOfSubtitution);
    free(approximatoinVector);
}


int main(int argc, char** argv) {
    double timeBefore, timeAfter;
    double minimalTime = ULLONG_MAX;
    //Ax = b
    double* A = initMatrix(MATRIX_SIZE, MATRIX_SIZE);
    double* b = initVector(MATRIX_SIZE);

    fillMatrixRandom(A, MATRIX_SIZE * MATRIX_SIZE);
    addToDiagonal(A, MATRIX_SIZE, 90);
    fillVectorRandom(b, MATRIX_SIZE);
    for(int i = 0; i < 1; i++){
        timeBefore = omp_get_wtime();
        double* solvation = useIterativeSolving(A, b, MATRIX_SIZE);
        timeAfter = omp_get_wtime();
        if (timeAfter - timeBefore < minimalTime)
            minimalTime = timeAfter - timeBefore;
    }
    //fprintf(stderr, "Result:\n");
    fprintf(stderr, "%lf \n", minimalTime);

    freeMatrix(A);
    freeVector(b);
    return 0;
}

