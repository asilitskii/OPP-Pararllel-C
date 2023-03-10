#include <stdio.h>
#include <malloc.h>
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
    /*setVectorCell(xCurrent, MATRIX_SIZE, 0, 1000);
    setVectorCell(xCurrent, MATRIX_SIZE, 1, 1001);*/
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
    int rankOfProc;
    int numOfProcs;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankOfProc);
    double timeBefore, timeAfter;
    //Ax = b
    double** A;
    double* b;
    int* partSizes;
    int* matrixPartSizes;
    int* matrixPartStarts;
    int* partStarts;
    if (rankOfProc == 0) {
        A = initMatrix(MATRIX_SIZE, MATRIX_SIZE);
        b = initVector(MATRIX_SIZE);
        fillMatrix(A, MATRIX_SIZE, 1);
        fillVector(b, MATRIX_SIZE, 2);
        printVector(b, 1, MATRIX_SIZE);
        printVector(*A, MATRIX_SIZE, MATRIX_SIZE);
        //partSizes = (int*)malloc(numOfProcs * sizeof(int));
        //partStarts = (int*)malloc(numOfProcs * sizeof(int));
        //for (int i = 0; i < numOfProcs; ++i) {
        //    int partStart = ((double)rankOfProc / numOfProcs) * MATRIX_SIZE * MATRIX_SIZE;
        //    int partEnd = ((double)(rankOfProc + 1) / numOfProcs) * MATRIX_SIZE * MATRIX_SIZE;
        //    //int partSize = (MATRIX_SIZE % numOfProcs == 0) ? partEnd - partStart : partEnd - partStart + 1;
        //    partSizes[i] = partEnd - partStart;
        //    partStarts[i] = partStart;
        //}

    }
    partSizes = (int*)malloc(numOfProcs * sizeof(int));
    matrixPartSizes = (int*)malloc(numOfProcs * sizeof(int));

    partStarts = (int*)malloc(numOfProcs * sizeof(int));
    matrixPartStarts = (int*)malloc(numOfProcs * sizeof(int));

    int partStart = ((double)rankOfProc / numOfProcs) * MATRIX_SIZE;
    int partEnd = ((double)(rankOfProc + 1) / numOfProcs) * MATRIX_SIZE;
    int partSize = partEnd - partStart;

    partSizes[rankOfProc] = partSize;
    partStarts[rankOfProc] = partStart;
    matrixPartSizes[rankOfProc] = partSize * MATRIX_SIZE;
    matrixPartStarts[rankOfProc] = partStart * MATRIX_SIZE;
    double** A_part = initMatrix(partSize, MATRIX_SIZE);
    double* b_part = initVector(partSize);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Scatterv(A, matrixPartSizes, matrixPartStarts, MPI_DOUBLE, A_part, partSize * MATRIX_SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(b, partSizes, partStarts, MPI_DOUBLE, b_part, partSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* solvation = useIterativeSolving(A, b, MATRIX_SIZE);
    freeMatrix(A, MATRIX_SIZE);
    freeVector(b);
    MPI_Finalize();
    return 0;
}

