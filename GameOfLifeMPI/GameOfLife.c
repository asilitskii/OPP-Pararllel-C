#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <time.h>
#include <limits.h>
#include <mpi.h>
#include "vectorOperations.h"

#define X 600
#define Y 300
#define MAX_ITERATIONS 24 * 100

void updateMap(char* currentGen, char* nextGen, int length){
    char* currentRow, *upperRow, *lowerRow;
    for(int i = 0; i < length / X; ++i){
        for(int j = 0; j < X; ++j){
            int neighboursNumber = 0;
            currentRow = &currentGen[i * X];
            upperRow = &currentGen[(i - 1) * X];
            lowerRow = &currentGen[(i + 1) * X];

            if(upperRow[(X + j - 1) % X]) neighboursNumber++; 
            if(currentRow[(X + j - 1) % X]) neighboursNumber++; 
            if(lowerRow[(X + j - 1) % X]) neighboursNumber++; 

            if(upperRow[j]) neighboursNumber++;
            if(lowerRow[j]) neighboursNumber++;

            if(upperRow[(j + 1) % X]) neighboursNumber++; 
            if(currentRow[(j + 1) % X]) neighboursNumber++; 
            if(lowerRow[(j + 1) % X]) neighboursNumber++; 
            if (currentRow[j] == 0) 
                if (neighboursNumber == 3) 
                    nextGen[i * X + j] = 1;
                else 
                    nextGen[i * X + j] = 0;

            if (currentRow[j] == 1) 
                if (neighboursNumber == 2 || neighboursNumber == 3)
                    nextGen[i * X + j] = 1;
                else
                    nextGen[i * X + j] = 0;
        }
    }
}

int checkFullHistoryForExit(char* historyCoincidencesFlags, int historyLength, int numOfProcs){
    for (int i = 0; i < historyLength; i++) {
        int matchedCount = 0;
        for (int j = 0; j < numOfProcs; j++) {
            if (historyCoincidencesFlags[i] == historyCoincidencesFlags[i + j * MAX_ITERATIONS] 
                && historyCoincidencesFlags[i] == 1) {
                matchedCount++;
            }
            else break;
        }
        if (matchedCount == numOfProcs) 
            return 1;
    }
    return 0;
}

void findHistoryCoincidence(char* map, int mapLength, char** history, int historyLength, char* coincidenceFlags){
    for(int i = 0; i < historyLength; ++i){
        if(memcmp(map, history[i], mapLength) == 0)
            coincidenceFlags[i] = 1;
        else
            coincidenceFlags[i] = 0;
    }
}

int main(int argc, char** argv) {
    int numOfProcs, rankOfProc;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankOfProc);

    double timeBefore, timeAfter;
    double minimalTime = ULLONG_MAX;

    char* A;
    int* partSizes;
    int* partStarts;

    if (rankOfProc == 0) {
        A = initMatrix(Y, X);
        fillGlider(A, X);
        partSizes = (int*)malloc(numOfProcs * sizeof(int));
        partStarts = (int*)malloc(numOfProcs * sizeof(int));
        //printVector(A, Y, X);
        for (int i = 0; i < numOfProcs; ++i) { //12x3 4pr 3str 
            int partStart = ((double)i / numOfProcs) * Y;
            int partEnd = ((double)(i + 1) / numOfProcs) * Y;
            partSizes[i] = (partEnd - partStart) * X;
            partStarts[i] = partStart * X;
        }
        timeBefore = MPI_Wtime();
    }

    int partStart = ((double)rankOfProc / numOfProcs) * Y;
    int partEnd = ((double)(rankOfProc + 1) / numOfProcs) * Y;
    int partSize = partEnd - partStart;
    char* A_part = initMatrix(partSize + 2, X);//store additional first and last row
    char* nextGen_part = initMatrix(partSize + 2, X);//store additional first and last row
    char** mapPartsHistory = (char**)malloc(MAX_ITERATIONS * sizeof(char*));
    char* historyPartCoincidenceFlags = (char*)malloc(MAX_ITERATIONS * numOfProcs * sizeof(char));
    char* historyCoincidencesFlags = (char*)malloc(MAX_ITERATIONS * numOfProcs * sizeof(char));
    MPI_Scatterv(A, partSizes, partStarts, MPI_CHAR, &(A_part[X]), partSize * X, MPI_CHAR, 0, MPI_COMM_WORLD);
    int upperProcRank =  (numOfProcs + (rankOfProc - 1)) % numOfProcs;
    int lowerProcRank = (rankOfProc + 1) % numOfProcs;
    MPI_Request requests[4];
    int tagForLowerReceiver = 123;
    int tagForUpperReceiver = 321;
    int currentIteration = 0;
    do{
        MPI_Isend(&(A_part[X]), X, MPI_CHAR, upperProcRank, tagForUpperReceiver, MPI_COMM_WORLD, &requests[0]);
        MPI_Isend(&(A_part[partSize * X]), X, MPI_CHAR, lowerProcRank, tagForLowerReceiver, MPI_COMM_WORLD, &requests[2]);
        MPI_Irecv(&(A_part[0]), X, MPI_CHAR, upperProcRank, tagForLowerReceiver, MPI_COMM_WORLD, &requests[1]);
        MPI_Irecv(&(A_part[(partSize + 1) * X]), X, MPI_CHAR, lowerProcRank, tagForUpperReceiver, MPI_COMM_WORLD, &requests[3]);
        
        mapPartsHistory[currentIteration] = (char*)malloc(partSize * X);
        memcpy(mapPartsHistory[currentIteration], &A_part[X], partSize * X * sizeof(char));
        findHistoryCoincidence(&A_part[X], partSize * X, mapPartsHistory, currentIteration, historyPartCoincidenceFlags);
        
        MPI_Request coincidenceRequest;
        MPI_Ialltoall(historyPartCoincidenceFlags, currentIteration, MPI_CHAR, historyCoincidencesFlags, currentIteration, MPI_CHAR, MPI_COMM_WORLD, &coincidenceRequest);
        //MPI_Iallreduce(historyPartCoincidenceFlags, historyCoincidencesFlags, currentIteration, MPI_CHAR, MPI_SUM, MPI_COMM_WORLD, &coincidenceRequest);
        
        updateMap(&(A_part[2 * X]), &(nextGen_part[2 * X]), (partSize - 2) * X);//owned without last and first
        
        MPI_Wait(&requests[0], MPI_STATUS_IGNORE);
        MPI_Wait(&requests[1], MPI_STATUS_IGNORE);
        updateMap(&(A_part[X]), &(nextGen_part[X]), X);//first
       
        MPI_Wait(&requests[2], MPI_STATUS_IGNORE);
        MPI_Wait(&requests[3], MPI_STATUS_IGNORE);
        updateMap(&(A_part[partSize * X]), &(nextGen_part[partSize * X]), X);//last
        
        MPI_Wait(&coincidenceRequest, MPI_STATUS_IGNORE);
        if(checkFullHistoryForExit(historyCoincidencesFlags, currentIteration, numOfProcs)){
            break;
        }
        char* tmp;
        tmp = A_part;
        A_part = nextGen_part;;
        nextGen_part = tmp;
        currentIteration++;
    } while(currentIteration < MAX_ITERATIONS);  
    //fprintf(stderr, "%d\n", currentIteration);
    //printVector(A_part, partSize + 2, X);  

    if (rankOfProc == 0) {
        timeAfter = MPI_Wtime(); 
        if (timeAfter - timeBefore < minimalTime)
            minimalTime = timeAfter - timeBefore;
        fprintf(stderr, "Result:\n");
        fprintf(stderr, "%lf \n", minimalTime);
        free(A);
        free(partSizes);
        free(partStarts);
    }
    for(int i = 0; i < currentIteration; ++i)
        free(mapPartsHistory[i]);
    free(mapPartsHistory);
    free(A_part);
    free(nextGen_part);
    free(historyPartCoincidenceFlags);
    free(historyCoincidencesFlags);
    MPI_Finalize();
    return 0;
}