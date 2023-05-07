#include <mpi.h>
#include <stdlib.h>

int main(int argc, char** argv){
	int worldRank, worldSize;
	const int dimsCount = 2;
	const int horisontalSize = 3, verticalSize = 2;

	MPI_Comm decartComm;
	int dims[dimsCount] = {horisontalSize, verticalSize};
    int periods[dimsCount] = {0, 0};
    int reorder = 0;
    int coords[dimsCount];
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    MPI_Cart_create(MPI_COMM_WORLD, dimsCount, dims, periods, reorder, &decartComm);
    int rank, size;
    MPI_Comm_size(decartComm, &size);
    MPI_Comm_rank(decartComm, &rank);
    MPI_Cart_coords(decartComm, rank, dimsCount, coords);

    MPI_Comm myColComm;
    MPI_Comm myRowComm;
    int myColCommRank, myColCommSize, myRowCommRank, myRowCommSize;
    int rowDims[dimsCount] = {0, 1};
    int colDims[dimsCount] = {1, 0};
    MPI_Cart_sub(decartComm, rowDims, &myRowComm);
    MPI_Cart_sub(decartComm, colDims, &myColComm);
    MPI_Comm_size(myRowComm, &myRowCommSize);
    MPI_Comm_rank(myRowComm, &myRowCommRank);
    MPI_Comm_size(myColComm, &myColCommSize);
    MPI_Comm_rank(myColComm, &myColCommRank);
    printf("My rank in decart comm is: %d, size: %d\nMy coords: %d, %d\nmyRowCommSize: %d, myRowCommRank: %d, myColCommSize: %d, myColCommRank: %d\n\n", 
        rank, size, coords[0], coords[1], myRowCommSize, myRowCommRank, myColCommSize, myColCommRank);
    
	return 0;
}