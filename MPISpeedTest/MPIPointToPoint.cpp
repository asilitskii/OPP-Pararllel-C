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
    double summaryTime = 0;
    for (int repeats = 0; repeats < K; ++repeats) {
        if (threadRank == 0) {
            double timeBefore, timeAfter;
            double result = 0;
            double semiResult = 0;
            int* a = (int*)calloc(N, sizeof(int));
            int* b = (int*)calloc(N, sizeof(int));
            fillArray(a, 1);
            fillArray(b, 1);
            timeBefore = clock();
            int pieceStart = ((double)0 / threadsNumber) * N;
            int pieceEnd = ((double)(1) / threadsNumber) * N;
            for (int i = 1; i < threadsNumber; ++i) {
                pieceStart = ((double)i / threadsNumber) * N;
                pieceEnd = ((double)(i + 1) / threadsNumber) * N;
                //printf("Sent : %d\n", pieceEnd - pieceStart);
                MPI_Send(&a[pieceStart], pieceEnd - pieceStart, MPI_INT, i, 1, MPI_COMM_WORLD);
                MPI_Send(b, N, MPI_INT, i, 11, MPI_COMM_WORLD);
            }
            result += operation(a, b, pieceEnd - pieceStart);
            for (int i = 1; i < threadsNumber; ++i) {
                MPI_Recv(&semiResult, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                result += semiResult;
            }
            if (threadRank == 0 && repeats == 0)
                printf("Result: %lf\n", result);
            timeAfter = clock();
            summaryTime += timeAfter - timeBefore;
            free(a);
            free(b);
        }
        else {
            int pieceSize = (int)(((double)(threadRank + 1) / threadsNumber) * N) - (int)(((double)threadRank / threadsNumber) * N);
            //printf("Recieved : %d\n", pieceSize);
            int* a = (int*)calloc(pieceSize, sizeof(int));
            int* b = (int*)calloc(N, sizeof(int));
            MPI_Recv(a, pieceSize, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(b, N, MPI_INT, 0, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            double semiResult = operation(a, b, pieceSize);
            //printf("semiResult = %lf\n", semiResult);
            MPI_Send(&semiResult, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
            free(a);
            free(b);
        }
    }
    if(threadRank == 0)
        printf("Result time: %lf\n", (summaryTime / K) / CLOCKS_PER_SEC);
    MPI_Finalize();
    return 0;
}