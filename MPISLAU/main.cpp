#include <stdio.h>
#include <time.h>
#include <mpi.h>

#include "vectorOperations.h"

#define MATRIX_SIZE 2
#define EPSILON 0.00001
#define TAU 0.01

double* resultOfSubtitutionToSystem(double** A, double* xCurrent, double* b, int size) {
    double* tmp = multMatrixByVector(A, xCurrent, size);
    double* result = subtractVectors(tmp, b, size);
    freeVector(tmp);
    return result;
}

double calculateTau(double** A, double* resultOfSubstitution, int size) {
    double tau = 0;
    double* multipliedByMatrix = multMatrixByVector(A, resultOfSubstitution, size);
    if (dotProduct(multipliedByMatrix, multipliedByMatrix, size) != 0)
        tau = dotProduct(resultOfSubstitution, multipliedByMatrix, size)
        / dotProduct(multipliedByMatrix, multipliedByMatrix, size);
    freeVector(multipliedByMatrix);
    return tau;
}

double* calcNextAproximation(double** A, double* xCurrent, double* b, int size) {
    double* resultOfSubtitution = resultOfSubtitutionToSystem(A, xCurrent, b, size);
    double* result = multVectByScalar(resultOfSubtitution, size, calculateTau(A, resultOfSubtitution, size));
    freeVector(resultOfSubtitution);
    return result;
}

double* calcNextSolution(double** A, double* xCurrent, double* b, int size) {
    double* currentApproximation = calcNextAproximation(A, xCurrent, b, size);
    double* nextSolution = subtractVectors(xCurrent, currentApproximation, size);
    freeVector(currentApproximation);
    return nextSolution;
}

double* useIterativeSolving(double** A, double* b, int size) {
    double* xCurrent = initVector(MATRIX_SIZE);
    //fillVectorRandom(xCurrent, MATRIX_SIZE);
    fillVector(xCurrent, MATRIX_SIZE, 1000);
    setVectorCell(xCurrent, MATRIX_SIZE, 0, 1000);
    setVectorCell(xCurrent, MATRIX_SIZE, 1, 1001);
    double* xNext = NULL;
    double* resultOfSubtitution = NULL;
    do {
        if (xNext != NULL) {
            freeVector(resultOfSubtitution);
            freeVector(xCurrent);
            xCurrent = xNext;
        }
        xNext = calcNextSolution(A, xCurrent, b, MATRIX_SIZE);
        resultOfSubtitution = resultOfSubtitutionToSystem(A, xNext, b, MATRIX_SIZE);
    } while (vectorNorm(resultOfSubtitution, MATRIX_SIZE)
        / vectorNorm(b, MATRIX_SIZE) >= EPSILON);
    freeVector(xCurrent);
    return xNext;
}

int main(int argc, char** argv) {

    double timeBefore, timeAfter;
    double* resultVector;
    //Ax = b
    double** A = initMatrix(MATRIX_SIZE);
    double* b = initVector(MATRIX_SIZE);

    double* tmp;
    double* tmp2;
    fillMatrix(A, MATRIX_SIZE, 1);

    /*setMatrixCell(A, MATRIX_SIZE, 0, 0, 1);
    setMatrixCell(A, MATRIX_SIZE, 0, 1, 2);
    setMatrixCell(A, MATRIX_SIZE, 1, 0, 1);
    setMatrixCell(A, MATRIX_SIZE, 1, 1, 1);*/

    fillVector(b, MATRIX_SIZE, 2);

    /*setVectorCell(b, MATRIX_SIZE, 0, 7);
    setVectorCell(b, MATRIX_SIZE, 1, 4);*/

    double* solvation = useIterativeSolving(A, b, MATRIX_SIZE);


    /*resultOfSubtitution = multMatrixByVector(A, xCurrent, MATRIX_SIZE);
    tmp2 = subtractVectors(resultOfSubtitution, b, MATRIX_SIZE);*/
    //timeBefore = clock();



    //timeAfter = clock();

    //printf("Result: %lf\n Result time: %lf", resultVector, (double)(timeAfter - timeBefore) / CLOCKS_PER_SEC);


    freeMatrix(A, MATRIX_SIZE);
    freeVector(b);
    return 0;
}

