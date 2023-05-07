#include <mpi.h>
#include <stdio.h>
#include <malloc.h>
#include <time.h>

#define N 200040
#define K 1

double operation(int* a, int* b, int size) {
    double s = 0;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < N; j++)
            s += (double)(a[i]) * (double)(b[j]);
    return s;
}

void fillArray(int* array, int value) {
    for (int i = 0; i < N; ++i)
        array[i] = value;
}

int main(int argc, char** argv) {
    int threadsNumber, threadRank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &threadsNumber);
    MPI_Comm_rank(MPI_COMM_WORLD, &threadRank);
    //printf("Hello, MPI world! I'm number %d from %d and I run.\n", threadRank, threadsNumber);
    double timeBefore = 0, timeAfter = 0;
    double result = 0;
    double semiResult = 0;
    double summaryTime = 0;
    int* a;
    if (threadRank == 0) {
        double semiResult = 0;
        a = (int*)calloc(N, sizeof(int));
        fillArray(a, 1);
        timeBefore = clock();
    }
    int* b = (int*)calloc(N, sizeof(int));
    fillArray(b, 1);

    int pieceStart = ((double)threadRank / threadsNumber) * N;
    int pieceEnd = ((double)(threadRank + 1) / threadsNumber) * N;
    int pieceSize = (N % threadsNumber == 0) ? pieceEnd - pieceStart : pieceEnd - pieceStart + 1;
    int* piece = (int*)calloc((pieceSize), sizeof(int));
    MPI_Scatter(a, pieceSize, MPI_INT, piece, pieceSize, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, N, MPI_INT, 0, MPI_COMM_WORLD);
    semiResult = operation(piece, b, pieceSize);
    MPI_Reduce(&semiResult, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (threadRank == 0) {
        timeAfter = clock();
        //printf("Result time: %lf\n", timeAfter - timeBefore);
        summaryTime = timeAfter - timeBefore;
        printf("Result: %lf\nResult time: %lf\n", result, summaryTime / CLOCKS_PER_SEC);
    }
    MPI_Finalize();
    return 0;
}