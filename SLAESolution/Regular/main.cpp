#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include "vectorOperations.h"

#define MATRIX_SIZE 10000
#define EPSILON 0.000001
#define MAX_ITERATIONS 10000
//#define TAU 0.01

double* resultOfSubtitutionToSystem(double* A, double* xCurrent, double* b, int size) {
    double* tmp = multMatrixByVector(A, xCurrent, size);
    double* result = subtractVectors(tmp, b, size);
    freeVector(tmp);
    return result;
}

double calculateTau(double* A, double* resultOfSubtitution, int size) {
    double tau = 0;
    double* dotProducts = (double*)malloc(2 * sizeof(double));
    double* multipliedByMatrix = multMatrixByVector(A, resultOfSubtitution, size);
    dotProducts[1] = dotProduct(multipliedByMatrix, multipliedByMatrix, size);
    if (dotProducts[1] != 0) 
        dotProducts[0] = dotProduct(resultOfSubtitution, multipliedByMatrix, size);
    else
        return 0;
    tau = dotProducts[0] / dotProducts[1];
    free(dotProducts);
    freeVector(multipliedByMatrix);
    return tau;
}

double* calcNextAproximation(double* A, double* xCurrent, double* b, double* resultOfSubtitution, int size) {
    //printf("tau %lf\n", calculateTau(A, resultOfSubtitution, size));
    double* result = multVectByScalar(resultOfSubtitution, size, calculateTau(A, resultOfSubtitution, size));
    return result;
}

double* calcNextSolution(double* A, double* xCurrent, double* b, double* resultOfSubtitution, int size) {
    double* currentApproximation = calcNextAproximation(A, xCurrent, b, resultOfSubtitution, size);
    double* nextSolution = subtractVectors(xCurrent, currentApproximation, size);
    freeVector(currentApproximation);
    return nextSolution;
}

double* useIterativeSolving(double* A, double* b, int size) {
    double* xCurrent = initVector(MATRIX_SIZE);
    fillVector(xCurrent, MATRIX_SIZE, 0);
    double* xNext = NULL;
    double* resultOfSubtitution = NULL;
    double exitCondition = EPSILON * EPSILON * squareVectorNorm(b, MATRIX_SIZE);

    int currentIteration = 0;
    resultOfSubtitution = resultOfSubtitutionToSystem(A, xCurrent, b, MATRIX_SIZE);
    while (currentIteration < MAX_ITERATIONS) {
        if (squareVectorNorm(resultOfSubtitution, MATRIX_SIZE) < exitCondition)
            break;
        currentIteration++;
        if (xNext != NULL) {
            //freeVector(resultOfSubtitution);
            freeVector(xCurrent);
            xCurrent = xNext;
        }
        xNext = calcNextSolution(A, xCurrent, b, resultOfSubtitution, MATRIX_SIZE);
        resultOfSubtitution = resultOfSubtitutionToSystem(A, xNext, b, MATRIX_SIZE);
    }
    printf("Number of iterations: %d\n", currentIteration);
    freeVector(xCurrent); 
    return xNext;
}

int main(int argc, char** argv) {
    double timeBefore, timeAfter;
    double minimalTime = ULLONG_MAX;
    //Ax = b
    double* A = initMatrix(MATRIX_SIZE, MATRIX_SIZE);
    double* b = initVector(MATRIX_SIZE);
  
    fillMatrixRandom(A, MATRIX_SIZE * MATRIX_SIZE);
    addToDiagonal(A, MATRIX_SIZE, 62);
    fillVectorRandom(b, MATRIX_SIZE);

    timeBefore = clock();
    double* solvation = useIterativeSolving(A, b, MATRIX_SIZE);
    timeAfter = clock();
    if (timeAfter - timeBefore < minimalTime)
        minimalTime = timeAfter - timeBefore;

    fprintf(stderr, "Result:\n");
    fprintf(stderr, "%lf \n", minimalTime / CLOCKS_PER_SEC);

    freeMatrix(A);
    freeVector(b);
    return 0;
}

