#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include "vectorOperations.h"

#define MATRIX_SIZE 2
#define EPSILON 0.001
#define TAU 0.01

double* resultOfSubtitutionToSystem(double** A, double* xCurrent, double* b, int size) {
    double* tmp = multMatrixByVector(A, xCurrent, size);
    double* result = subtractVectors(tmp, b, size);
    freeVector(tmp);
    return result;
}

double* calcNextAproximation(double** A, double* xCurrent, double* b, int size) {
    double* tmp = resultOfSubtitutionToSystem(A, xCurrent, b, size);
    double* result = multVectByScalar(tmp, size, TAU);
    freeVector(tmp);
    return result;
}

double* calcNextSolution(double** A, double* xCurrent, double* b, int size) {
    double* currentApproximation = calcNextAproximation(A, xCurrent, b, size);
    double* nextSolution = subtractVectors(xCurrent, currentApproximation, size);
    freeVector(currentApproximation);
    return nextSolution;
}

int main(int argc, char** argv) {
    double timeBefore, timeAfter;
    double* resultVector;
    //Ax = b
    double** A = initMatrix(MATRIX_SIZE);
    double* b = initVector(MATRIX_SIZE);
    double* x = initVector(MATRIX_SIZE);
    
    double* tmp;
    double* tmp2;
    double* xCurrent;
    double* xNext;

    fillMatrix(A, MATRIX_SIZE, 1);
    fillVector(b, MATRIX_SIZE, 2);
    fillVector(x, MATRIX_SIZE, 2);
    
    tmp = multMatrixByVector(A, x, MATRIX_SIZE);
    tmp2 = subtractVectors(tmp, b, MATRIX_SIZE);
    //timeBefore = clock();



    //timeAfter = clock();

    //printf("Result: %lf\n Result time: %lf", resultVector, (double)(timeAfter - timeBefore) / CLOCKS_PER_SEC);


    freeMatrix(A, MATRIX_SIZE);
    freeVector(b);
    freeVector(x);
    freeVector(tmp);
    return 0;
}

