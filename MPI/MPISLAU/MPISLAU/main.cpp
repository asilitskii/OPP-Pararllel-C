#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <mpi.h>
#include "vectorOperations.h"

#define MATRIX_SIZE 10000
#define EPSILON 0.000001
#define MAX_ITERATIONS 10000

int rankOfProc;
int numOfProcs;

double* resultOfSubtitutionToSystem(double* A, double* xCurrent, double* b, int size, int partSize) {
    //fprintf(stderr, "4\n");
    double* tmp = multMatrixByVector(A, partSize, xCurrent, size);
    double* result = subtractVectors(tmp, b, partSize);
    freeVector(tmp);
    /*fprintf(stderr, "sres0 %lf\n", result[0]);
    fprintf(stderr, "sres1 %lf\n", result[1]);*/
    return result;
}

double calculateTau(double* A_part, double* resultOfSubtitution, int size,
                    int partSize, int* partSizes, int* partStarts) {
    double tau = 0;
    double* dotProducts_part = (double*)malloc(2 * sizeof(double));
    double* dotProducts = (double*)malloc(2 * sizeof(double));
    //double* resultOfSubtitution_part = initVector(size);
    //fprintf(stderr, "sres0 %lf\n", resultOfSubtitution_part[0]);
    //MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Allgatherv(resultOfSubtitution_part, partSize, MPI_DOUBLE, resultOfSubtitution_part, partSizes, partStarts, MPI_DOUBLE, MPI_COMM_WORLD);
    //fprintf(stderr, "s0 %lf\n", resultOfSubtitution_part[0]);
    //fprintf(stderr, "s1 %lf\n", resultOfSubtitution_part[1]);

    double* multipliedByMatrix_part = multMatrixByVector(A_part, partSize, resultOfSubtitution, size);
    /*fprintf(stderr, "mbm0 %lf\n", multipliedByMatrix_part[0]);
    fprintf(stderr, "mbm1 %lf\n", multipliedByMatrix_part[1]);*/
    dotProducts_part[1] = dotProduct(multipliedByMatrix_part, multipliedByMatrix_part, partSize);
    /*fprintf(stderr, "dp0 %lf\n", dotProducts_part[0]);
    fprintf(stderr, "dp1 %lf\n", dotProducts_part[1]);*/
    if (dotProducts_part[1] != 0) {
        dotProducts_part[0] = dotProduct(&(resultOfSubtitution[partStarts[rankOfProc]]), multipliedByMatrix_part, partSize);
    }
    else 
        return 0; 

    /*fprintf(stderr, "mm0 %lf\n", multipliedByMatrix_part[0]);
    fprintf(stderr, "dpp0 %lf\n", dotProducts_part[0]);
    fprintf(stderr, "dpp1 %lf\n", dotProducts_part[1]);
    MPI_Barrier(MPI_COMM_WORLD);*/
    MPI_Allreduce(dotProducts_part, dotProducts, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    /*fprintf(stderr, "dp0 %lf\n", dotProducts[0]);
    fprintf(stderr, "dp1 %lf\n", dotProducts[1]);*/
    tau = dotProducts[0] / dotProducts[1];
    //freeVector(resultOfSubtitution_part);
    free(dotProducts);
    free(dotProducts_part);
    freeVector(multipliedByMatrix_part);
    return tau;
}

double* calcNextAproximation(double* A_part, double* xCurrent, double* b_part, int size, 
                             double* resultOfSubtitution, int partSize, int* partSizes, int* partStarts) {
    double tau = calculateTau(A_part, resultOfSubtitution, size, partSize, partSizes, partStarts);
    //fprintf(stderr, "tau %lf\n", tau);
    double* result = multVectByScalar(resultOfSubtitution, size, tau);
    //fprintf(stderr, "res0 %lf\n", result[0]); 
    return result;
}

double* calcNextSolution(double* A_part, double* xCurrent, double* b_part, int size, 
                         double* resultOfSubtitution, int partSize, int* partSizes, int* partStarts) {

    double* currentApproximation = calcNextAproximation(A_part, xCurrent, b_part, size,
                                                        resultOfSubtitution, partSize, partSizes, partStarts);
    //double* currentApproximation = initVector(size);
    /*fprintf(stderr, "app0 %lf\n", currentApproximation[0]);
    fprintf(stderr, "app1 %lf\n", currentApproximation[1]);
    fprintf(stderr, "pss0 %d \n", partSizes[0]);
    fprintf(stderr, "pss1 %d \n", partSizes[1]);
    fprintf(stderr, "pst0 %d \n", partStarts[0]);
    fprintf(stderr, "pst1 %d \n", partStarts[1]);
    fprintf(stderr, "ps %d \n", partSize);*/

    //MPI_Allgatherv(currentApproximation_part, partSize, MPI_DOUBLE, currentApproximation, partSizes, partStarts, MPI_DOUBLE, MPI_COMM_WORLD);
    /*fprintf(stderr, "ap0 %lf\n", currentApproximation[0]);
    fprintf(stderr, "ap1 %lf\n", currentApproximation[1]);*/

    double* nextSolution = subtractVectors(xCurrent, currentApproximation, size);
    /*fprintf(stderr, "sol0 %lf\n", nextSolution[0]);
    fprintf(stderr, "sol1 %lf\n", nextSolution[1]);*/
    //freeVector(currentApproximation_part);
    //freeVector(currentApproximation);

    freeVector(currentApproximation);
    //printVector(nextSolution, 1, MATRIX_SIZE);
    return nextSolution;
}

double* useIterativeSolving(double* A_part, double* b_part, int size,
                            int partSize, int* partSizes, int* partStarts) {
    double* xCurrent = initVector(size);
    //fillVectorRandom(xCurrent, MATRIX_SIZE);
    fillVector(xCurrent, size, 0);
    /*setVectorCell(xCurrent, MATRIX_SIZE, 0, 1000);
    setVectorCell(xCurrent, MATRIX_SIZE, 1, 1001);*/

    double* xNext = initVector(size);
    
    int currentIteration = 0;
    //printVector(b_part, 1, partSize);
    double exitCondition_part = EPSILON * EPSILON * squareVectorNorm(b_part, partSize);
    //fprintf(stderr, "exp %lf\n", exitCondition_part);
    double exitCondition;
    MPI_Allreduce(&exitCondition_part, &exitCondition, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double* resultOfSubtitution = initVector(size);
    double* resultOfSubtitution_part = resultOfSubtitutionToSystem(A_part, xCurrent, b_part, size, partSize);


    /*printVector(A_part, partSize, MATRIX_SIZE);
    fprintf(stderr, "\n");*/

    /*fprintf(stderr, "rs0 %lf\n", resultOfSubtitution_part[0]);
    fprintf(stderr, "rs1 %lf\n", resultOfSubtitution_part[1]);*/
    //fprintf(stderr, "sq %lf\n", squareVectorNorm(resultOfSubtitution, size));
    //fprintf(stderr, "ex %lf\n", exitCondition);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgatherv(resultOfSubtitution_part, partSizes[rankOfProc], MPI_DOUBLE, resultOfSubtitution, partSizes, partStarts, MPI_DOUBLE, MPI_COMM_WORLD);
    //fprintf(stderr, "22 %d\n", rankOfProc);
    while (currentIteration < MAX_ITERATIONS) {
        if (squareVectorNorm(resultOfSubtitution, size) < exitCondition) {
            //fprintf(stderr, "Solvation found\n");
            break;
        }
        //else fprintf(stderr, "Go next\n");
        currentIteration++;
        if (xNext != NULL) {
            freeVector(xCurrent);
            xCurrent = xNext;
        }
        //printVector(xCurrent, 1, MATRIX_SIZE);
        //fprintf(stderr, "NextSol\n");
        xNext = calcNextSolution(A_part, xCurrent, b_part, size, resultOfSubtitution, partSize, partSizes, partStarts);
        //fprintf(stderr, "xn0 %d \n", xNext[0]);
        //fprintf(stderr, "xn1 %d \n", xNext[1]);
        /*fprintf(stderr, "pss0 %d \n", partSizes[0]);
        fprintf(stderr, "pss1 %d \n", partSizes[1]);
        fprintf(stderr, "pst0 %d \n", partStarts[0]);
        fprintf(stderr, "pst1 %d \n", partStarts[1]);
        fprintf(stderr, "ps %d \n", partSize);*/
        //MPI_Barrier(MPI_COMM_WORLD);
        //MPI_Allgatherv(xNext_part, partSize, MPI_DOUBLE, xNext, partSizes, partStarts, MPI_DOUBLE, MPI_COMM_WORLD);
        freeVector(resultOfSubtitution_part);
        //freeVector(resultOfSubtitution);
        //fprintf(stderr, "123\n");
        resultOfSubtitution_part = resultOfSubtitutionToSystem(A_part, xNext, b_part, MATRIX_SIZE, partSize);
        //double* tmp = multMatrixByVector(A_part, partSize, xNext, size);
        ////fprintf(stderr, "4%d\n", rankOfProc);
        //double* tmp2 = subtractVectors(tmp, b_part, partSize);
        ////fprintf(stderr, "5%d\n", rankOfProc);
        //resultOfSubtitution_part = tmp2;
        ////fprintf(stderr, "6%d\n", rankOfProc);
        //free(tmp);

        MPI_Allgatherv(resultOfSubtitution_part, partSizes[rankOfProc], MPI_DOUBLE, resultOfSubtitution, partSizes, partStarts, MPI_DOUBLE, MPI_COMM_WORLD);
        //freeVector(xNext_part);
        //fprintf(stderr, "Found\n");
        //fprintf(stderr, "VecN1%lf\n", vectorNorm(resultOfSubtitution_part, MATRIX_SIZE));
        //fprintf(stderr, "VecN2%lf\n", vectorNorm(b_part, MATRIX_SIZE));
        //fprintf(stderr, "VecN3%lf\n", vectorNorm(resultOfSubtitution_part, MATRIX_SIZE) / vectorNorm(b_part, MATRIX_SIZE));
        //if (vectorNorm(resultOfSubtitution_part, MATRIX_SIZE) / vectorNorm(b_part, MATRIX_SIZE) >= EPSILON)
        //    fprintf(stderr, "Go next\n");


    }
    //vectorNorm(resultOfSubtitution_part, MATRIX_SIZE) / vectorNorm(b_part, MATRIX_SIZE) >= EPSILON
    //freeVector(xCurrent);
    freeVector(resultOfSubtitution_part);
    //freeVector(resultOfSubtitution);
    if(rankOfProc == 0)
        printf("Number of iterations: %d\n", currentIteration);
    return xNext;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankOfProc);
    //fprintf(stderr, "Hello\n");
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
        /*fillMatrix(A, MATRIX_SIZE, MATRIX_SIZE, 1);
        fillVector(b, MATRIX_SIZE, 2);*/
        /*fillMatrix(A, MATRIX_SIZE, MATRIX_SIZE, 1);
        fillDiagonal(A, MATRIX_SIZE, 2);
        fillVector(b, MATRIX_SIZE, MATRIX_SIZE + 1);*/
        fillMatrixRandom(A, MATRIX_SIZE * MATRIX_SIZE);
        addToDiagonal(A, MATRIX_SIZE, 62);
        fillVectorRandom(b, MATRIX_SIZE);
        //printVector(A, MATRIX_SIZE, MATRIX_SIZE);
    }
    //сколько каждому отдать
    for (int i = 0; i < numOfProcs; ++i) {
        int partStart = ((double)i / numOfProcs) * MATRIX_SIZE;
        int partEnd = ((double)(i + 1) / numOfProcs) * MATRIX_SIZE;
        //int partSize = (MATRIX_SIZE % numOfProcs == 0) ? partEnd - partStart : partEnd - partStart + 1;
        partSizes[i] = partEnd - partStart;
        partStarts[i] = partStart;
        matrixPartSizes[i] = partSizes[i] * MATRIX_SIZE;
        matrixPartStarts[i] = partStarts[i] * MATRIX_SIZE;
    }

    //MPI_Bcast(partSizes, numOfProcs, MPI_INT, 0, MPI_COMM_WORLD);
    //MPI_Bcast(partStarts, numOfProcs, MPI_INT, 0, MPI_COMM_WORLD);
    
    /*printVector(b_part, 1, MATRIX_SIZE);
    printVector(A_part, MATRIX_SIZE, MATRIX_SIZE);
    fprintf(stderr, "\n");*/
    
    //сколько каждый примет
    int partStart = ((double)rankOfProc / numOfProcs) * MATRIX_SIZE;
    int partEnd = ((double)(rankOfProc + 1) / numOfProcs) * MATRIX_SIZE;
    int partSize = partEnd - partStart;
    
    /*partSizes = (int*)malloc(numOfProcs * sizeof(int));
    matrixPartSizes = (int*)malloc(numOfProcs * sizeof(int));

    partStarts = (int*)malloc(numOfProcs * sizeof(int));
    matrixPartStarts = (int*)malloc(numOfProcs * sizeof(int));*/

    /*partSizes[rankOfProc] = partSize;
    partStarts[rankOfProc] = partStart;
    matrixPartSizes[rankOfProc] = partSize * MATRIX_SIZE;
    matrixPartStarts[rankOfProc] = partStart * MATRIX_SIZE;*/
    double* A_part = initMatrix(partSize, MATRIX_SIZE);
    double* b_part = initVector(partSize);

    /*printVector(b_part, 1, partSize);
    printVector(A_part, partSize, MATRIX_SIZE);
    fprintf(stderr, "\n");*/
    /*if (rankOfProc == 0) {
        fprintf(stderr, "mpss0 %d \n", matrixPartSizes[0]);
        fprintf(stderr, "mpss1 %d \n", matrixPartSizes[1]);
        fprintf(stderr, "pss0 %d \n", partSizes[0]);
        fprintf(stderr, "pss1 %d \n", partSizes[1]);
    }
    fprintf(stderr, "ps %d \n", partSize);*/

    MPI_Scatterv(b, partSizes, partStarts, MPI_DOUBLE, b_part, partSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //printVector(b_part, 1, partSize);
    //fprintf(stderr, " vec part in %d\n", rankOfProc);
    MPI_Scatterv(A, matrixPartSizes, matrixPartStarts, MPI_DOUBLE, A_part, partSize * MATRIX_SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //printVector(A_part, partSize, MATRIX_SIZE);
    //fprintf(stderr, " mat part in %d\n", rankOfProc);

    /*if (rankOfProc == 0) {
        fprintf(stderr, "mpss0 %d \n", matrixPartSizes[0]);
        fprintf(stderr, "mpss1 %d \n", matrixPartSizes[1]);
        fprintf(stderr, "mpst0 %d \n", matrixPartStarts[0]);
        fprintf(stderr, "mpst1 %d \n", matrixPartStarts[1]);
        fprintf(stderr, "pss0 %d \n", partSizes[0]);
        fprintf(stderr, "pss1 %d \n", partSizes[1]);
        fprintf(stderr, "pst0 %d \n", partStarts[0]);
        fprintf(stderr, "pst1 %d \n", partStarts[1]);
    }
    fprintf(stderr, "ps %d \n", partSize);*/
    /*fprintf(stderr, "mpss0 %d \n", matrixPartSizes[0]);
    fprintf(stderr, "mpss1 %d \n", matrixPartSizes[1]);
    fprintf(stderr, "mpst0 %d \n", matrixPartStarts[0]);
    fprintf(stderr, "mpst1 %d \n", matrixPartStarts[1]);
    fprintf(stderr, "pss0 %d \n", partSizes[0]);
    fprintf(stderr, "pss1 %d \n", partSizes[1]);
    fprintf(stderr, "pst0 %d \n", partStarts[0]);
    fprintf(stderr, "pst1 %d \n", partStarts[1]);*/


    //MPI_Scatterv(b_part, partSizes, partStarts, MPI_DOUBLE, b_part, partSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /*printVector(b_part, 1, partSize);
    fprintf(stderr, " vec part in %d\n", rankOfProc);*/

    /*printVector(b_part, 1, partSize);
    printVector(A_part, partSize, MATRIX_SIZE);
    fprintf(stderr, "\n");*/
    timeBefore = MPI_Wtime();


    double* solvation = useIterativeSolving(A_part, b_part, MATRIX_SIZE, partSize, partSizes, partStarts);

    timeAfter = MPI_Wtime();
    if (timeAfter - timeBefore < minimalTime)
        minimalTime = timeAfter - timeBefore;

    if (rankOfProc == 0) {
        fprintf(stderr, "Result:\n");
        //printVector(solvation, 1, MATRIX_SIZE);
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

