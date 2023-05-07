#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <limits.h>
#include <mpi.h>
#include "vectorOperations.h"

//#define MATRIX_SIZE 6
#define M 8400
#define N 3600
#define K 8400
#define PROCS_IN_ROW 3
#define PROCS_IN_COLUMN 8
#define DIMS_COUNT 2

int main(int argc, char *argv[]) {
    int size;
    int rank;
    int verticalPartSize = M / PROCS_IN_COLUMN;
    int horisontalPartSize = K / PROCS_IN_ROW;
    double time_start;
    double time_end;

    MPI_Init(&argc, &argv);
    MPI_Comm MPI_COMM_DECART;
    int dims[DIMS_COUNT] = {PROCS_IN_COLUMN, PROCS_IN_ROW};
    int periods[DIMS_COUNT] = {0, 0};
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, DIMS_COUNT, dims, periods, reorder, &MPI_COMM_DECART);

    MPI_Comm MPI_COMM_ROW;
    int rowDims[DIMS_COUNT] = {0, 1};
    MPI_Cart_sub(MPI_COMM_DECART, rowDims, &MPI_COMM_ROW);

    MPI_Comm MPI_COMM_COLUMN;
    int colDims[DIMS_COUNT] = {1, 0};
    MPI_Cart_sub(MPI_COMM_DECART, colDims, &MPI_COMM_COLUMN);

    int coords[DIMS_COUNT] = {0, 0};
    MPI_Comm_size(MPI_COMM_DECART, &size);
    MPI_Comm_rank(MPI_COMM_DECART, &rank);
    MPI_Cart_coords(MPI_COMM_DECART, rank, DIMS_COUNT, coords);

    MPI_Datatype MPI_ROW;
    MPI_Type_contiguous(N * verticalPartSize, MPI_DOUBLE, &MPI_ROW);
    MPI_Type_commit(&MPI_ROW);

    MPI_Datatype MPI_UNSIZED_COLUMN, MPI_COLUMN;
    MPI_Type_vector(N, horisontalPartSize, K, MPI_DOUBLE, &MPI_UNSIZED_COLUMN);
    MPI_Type_commit(&MPI_UNSIZED_COLUMN);
    MPI_Type_create_resized(MPI_UNSIZED_COLUMN, 0, horisontalPartSize * sizeof(double), &MPI_COLUMN);
    MPI_Type_commit(&MPI_COLUMN);

    MPI_Datatype MPI_UNSIZED_BLOCK, MPI_BLOCK;
    MPI_Type_vector(verticalPartSize, horisontalPartSize, K, MPI_DOUBLE, &MPI_UNSIZED_BLOCK);
    MPI_Type_commit(&MPI_UNSIZED_BLOCK);
    MPI_Type_create_resized(MPI_UNSIZED_BLOCK, 0, horisontalPartSize * sizeof(double), &MPI_BLOCK);
    MPI_Type_commit(&MPI_BLOCK);

    double* A;
    double* B;
    double* C;
    if (rank == 0) {
        A = initMatrix(M, N);
        B = initMatrix(N, K);
        fillMatrixRandom(A, M, N);
        //printVector(A, M, N);
        //fillMatrixRandom(B, N, K);
        addToDiagonal(B, N, 1);
        //printVector(B, N, K);
        C = initMatrix(M, K);
        time_start = MPI_Wtime();
    }
    double* A_part = initMatrix(verticalPartSize, N);
    double* B_part = initMatrix(N, horisontalPartSize);
    double* C_part = initMatrix(verticalPartSize, horisontalPartSize);

    if (coords[1] == 0)
        MPI_Scatter(A, 1, MPI_ROW, A_part, 1, MPI_ROW, 0, MPI_COMM_COLUMN);
    MPI_Bcast(A_part, verticalPartSize * N, MPI_DOUBLE, 0, MPI_COMM_ROW);

    if (coords[0] == 0)
        MPI_Scatter(B, 1, MPI_COLUMN,  B_part, horisontalPartSize * N, MPI_DOUBLE, 0, MPI_COMM_ROW);
    MPI_Bcast(B_part, horisontalPartSize * N, MPI_DOUBLE, 0, MPI_COMM_COLUMN);

    for (int i = 0; i < verticalPartSize; ++i) {
        for (int k = 0; k < N; ++k) {
            for (int j = 0; j < horisontalPartSize; ++j)
                C_part[i * horisontalPartSize + j] += A_part[i * N + k] * B_part[k * horisontalPartSize + j];
        }
    }

    int* recvBlockSizes;
    int* recvBlockStarts;
    if(rank == 0){
        recvBlockSizes = (int*)malloc(size * sizeof(int));
        recvBlockStarts = (int*)malloc(size * sizeof(int));
        for(int i = 0; i < size; ++i){
            recvBlockSizes[i] = 1;
            int senderCoords[DIMS_COUNT];
            MPI_Cart_coords(MPI_COMM_DECART, i, DIMS_COUNT, senderCoords);
            recvBlockStarts[i] = senderCoords[0] * verticalPartSize * PROCS_IN_ROW + senderCoords[1];
        }
    }

    MPI_Gatherv(C_part, verticalPartSize * horisontalPartSize, MPI_DOUBLE, C, recvBlockSizes, recvBlockStarts, MPI_BLOCK, 0, MPI_COMM_DECART);

    if (rank == 0) {
        time_end = MPI_Wtime();
        printf("Total time in seconds: %f\n", time_end - time_start);
        //printVector(C, M, K);
        free(A);
        free(B);
        free(C);
        free(recvBlockSizes);
        free(recvBlockStarts);
    }
    MPI_Type_free(&MPI_ROW);
    MPI_Type_free(&MPI_UNSIZED_COLUMN);
    MPI_Type_free(&MPI_COLUMN);
    MPI_Type_free(&MPI_UNSIZED_BLOCK);
    MPI_Type_free(&MPI_BLOCK);
    free(A_part);
    free(B_part);
    free(C_part);

    MPI_Finalize();
    return EXIT_SUCCESS;
}

