#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <limits.h>
#include <mpi.h>
#include "VectorOperations.h"

#define MATRIX_SIZE 10000
#define EPSILON 0.000001
#define MAX_ITERATIONS 10000

int rankOfProc;
int numOfProcs;

double* resultOfSubtitutionToSystem(double* A, double* xCurrent, double* b, int size, int partSize) {
    double* tmp = multMatrixByVector(A, partSize, xCurrent, size);
    double* result = subtractVectors(tmp, b, partSize);
    freeVector(tmp);
    return result;
}

double calculateTau(double* A_part, double* resultOfSubtitution, int size,
                    int partSize, int* partSizes, int* partStarts) {
    double tau = 0;
    double* dotProducts_part = (double*)malloc(2 * sizeof(double));
    double* dotProducts = (double*)malloc(2 * sizeof(double));

    double* multipliedByMatrix_part = multMatrixByVector(A_part, partSize, resultOfSubtitution, size);
    dotProducts_part[1] = dotProduct(multipliedByMatrix_part, multipliedByMatrix_part, partSize);
    if (dotProducts_part[1] != 0)
        dotProducts_part[0] = dotProduct(&(resultOfSubtitution[partStarts[rankOfProc]]), multipliedByMatrix_part, partSize);
    else 
        return 0; 
    MPI_Allreduce(dotProducts_part, dotProducts, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    tau = dotProducts[0] / dotProducts[1];
    free(dotProducts);
    free(dotProducts_part);
    freeVector(multipliedByMatrix_part);
    return tau;
}

double* calcNextAproximation(double* A_part, double* xCurrent, double* b_part, int size, 
                             double* resultOfSubtitution, int partSize, int* partSizes, int* partStarts) {
    double tau = calculateTau(A_part, resultOfSubtitution, size, partSize, partSizes, partStarts);
    double* result = multVectByScalar(resultOfSubtitution, size, tau);
    return result;
}

double* calcNextSolution(double* A_part, double* xCurrent, double* b_part, int size, 
                         double* resultOfSubtitution, int partSize, int* partSizes, int* partStarts) {

    double* currentApproximation = calcNextAproximation(A_part, xCurrent, b_part, size,
                                                        resultOfSubtitution, partSize, partSizes, partStarts);
    double* nextSolution = subtractVectors(xCurrent, currentApproximation, size);
    freeVector(currentApproximation);
    return nextSolution;
}

double* useIterativeSolving(double* A_part, double* b_part, int size,
                            int partSize, int* partSizes, int* partStarts) {
    double* xCurrent = initVector(size);
    fillVector(xCurrent, size, 0);

    double* xNext = initVector(size);
    
    int currentIteration = 0;

    double exitCondition_part = EPSILON * EPSILON * squareVectorNorm(b_part, partSize);
    double exitCondition;
    MPI_Allreduce(&exitCondition_part, &exitCondition, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double* resultOfSubtitution = initVector(size);
    double* resultOfSubtitution_part = resultOfSubtitutionToSystem(A_part, xCurrent, b_part, size, partSize);

    MPI_Allgatherv(resultOfSubtitution_part, partSizes[rankOfProc], MPI_DOUBLE, resultOfSubtitution, partSizes, partStarts, MPI_DOUBLE, MPI_COMM_WORLD);
    while (currentIteration < MAX_ITERATIONS) {
        if (squareVectorNorm(resultOfSubtitution, size) < exitCondition) 
            break;
        currentIteration++;
        if (xNext != NULL) {
            freeVector(xCurrent);
            xCurrent = xNext;
        }
        xNext = calcNextSolution(A_part, xCurrent, b_part, size, resultOfSubtitution, partSize, partSizes, partStarts);

        freeVector(resultOfSubtitution_part);
        resultOfSubtitution_part = resultOfSubtitutionToSystem(A_part, xNext, b_part, MATRIX_SIZE, partSize);
        MPI_Allgatherv(resultOfSubtitution_part, partSizes[rankOfProc], MPI_DOUBLE, resultOfSubtitution, partSizes, partStarts, MPI_DOUBLE, MPI_COMM_WORLD);

    }

    freeVector(resultOfSubtitution_part);
    freeVector(resultOfSubtitution);
    if(rankOfProc == 0)
        printf("Number of iterations: %d\n", currentIteration);
    return xNext;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankOfProc);
    double timeBefore, timeAfter;
    double minimalTime = ULLONG_MAX;
    //Ax = b
    double* A;
    double* b;
    int* partSizes;
    int* partStarts;
    int* matrixPartSizes;
    int* matrixPartStarts;

    partSizes = (int*)malloc(numOfProcs * sizeof(int));
    partStarts = (int*)malloc(numOfProcs * sizeof(int));
    matrixPartSizes = (int*)malloc(numOfProcs * sizeof(int));
    matrixPartStarts = (int*)malloc(numOfProcs * sizeof(int));

    if (rankOfProc == 0) {
        A = initMatrix(MATRIX_SIZE, MATRIX_SIZE);
        b = initVector(MATRIX_SIZE);
        fillMatrixRandom(A, MATRIX_SIZE * MATRIX_SIZE);
        addToDiagonal(A, MATRIX_SIZE, 85);
        fillVectorRandom(b, MATRIX_SIZE);
    }
    for (int i = 0; i < numOfProcs; ++i) {
        int partStart = ((double)i / numOfProcs) * MATRIX_SIZE;
        int partEnd = ((double)(i + 1) / numOfProcs) * MATRIX_SIZE;
        partSizes[i] = partEnd - partStart;
        partStarts[i] = partStart;
        matrixPartSizes[i] = partSizes[i] * MATRIX_SIZE;
        matrixPartStarts[i] = partStarts[i] * MATRIX_SIZE;
    }
    
    int partStart = ((double)rankOfProc / numOfProcs) * MATRIX_SIZE;
    int partEnd = ((double)(rankOfProc + 1) / numOfProcs) * MATRIX_SIZE;
    int partSize = partEnd - partStart;

    double* A_part = initMatrix(partSize, MATRIX_SIZE);
    double* b_part = initVector(partSize); 

    MPI_Scatterv(b, partSizes, partStarts, MPI_DOUBLE, b_part, partSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Scatterv(A, matrixPartSizes, matrixPartStarts, MPI_DOUBLE, A_part, partSize * MATRIX_SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    timeBefore = MPI_Wtime();
    double* solvation = useIterativeSolving(A_part, b_part, MATRIX_SIZE, partSize, partSizes, partStarts);
    timeAfter = MPI_Wtime();
    if (timeAfter - timeBefore < minimalTime)
        minimalTime = timeAfter - timeBefore;

    if (rankOfProc == 0) {
        fprintf(stderr, "Result:\n");
        fprintf(stderr, "%lf \n", minimalTime);
        freeMatrix(A);
        freeVector(b);
    }
    freeMatrix(A_part);
    freeVector(b_part);
    free(partSizes);
    free(partStarts);
    if (rankOfProc == 0) {
        free(matrixPartSizes);
        free(matrixPartStarts);
    }
    MPI_Finalize();
    return 0;
}