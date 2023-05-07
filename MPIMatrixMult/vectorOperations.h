#pragma once

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

double* initMatrix(int height, int width) {
    double* matrix = (double*)calloc(height * width, sizeof(double));
    return matrix;
}

double* initVector(int size) {
    double* vector = (double*)calloc(size, sizeof(double));
    return vector;
}

void fillVectorRandom(double* array, int size) {
    srand(size);
    for (int i = 0; i < size; ++i)
        array[i] = (double)rand() / RAND_MAX * (double)rand();
}

void fillMatrixRandom(double* martix, int height, int width) {
    srand(height);
    for (int i = 0; i < height * width; ++i)
        martix[i] = (double)rand() / RAND_MAX * 2.0 - 1.0;
}

void addToDiagonal(double* matrix, int size, double value) {
    for (int i = 0; i < size; ++i)
        matrix[i * size + i] += value;
}

void fillVector(double* array, int size, double value) {
    for (int i = 0; i < size; ++i)
        array[i] = value;
}

void fillMatrix(double* matrix, int height, int width, double value) {
    for (int i = 0; i < height; ++i)
        for (int j = 0; j < width; ++j)
            matrix[i * height + j] = value;
}

void fillDiagonal(double* matrix, int size, double value) {
    for (int i = 0; i < size; ++i)
        matrix[i * size + i] = value;
}

void setMatrixCell(double* matrix, int height, int width, int i, int j, double value) {
    matrix[i * height + j] = value;
}

void setVectorCell(double* vector, int size, int i, double value) {
    vector[i] = value;
}

void printVector(double* vector, int rows, int cols) {
    /*for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("%lf ", vector[i * rows + j]);
        }
        printf("\n");
    }*/
    for(int i = 0; i < rows * cols; ++i){
        printf("%lf ", vector[i]);
        if(i % cols == cols - 1)
            printf("\n");
    }

}

void freeMatrix(double* matrix) {
    free(matrix);
}

void freeVector(double* vector) {
    free(vector);
}

double dotProduct(double* a, double* b, int size) {
    double result = 0;
    for (int i = 0; i < size; i++)
        result += a[i] * b[i];
    return result;
}

double squareVectorNorm(double* a, int size) {
    double squareNorm = 0;
    for (int i = 0; i < size; ++i)
        squareNorm += a[i] * a[i];
    return squareNorm;
}

double* subtractVectors(double* a, double* b, int size) {
    double* result = (double*)malloc(size * sizeof(double));
    for (int i = 0; i < size; i++)
        result[i] = a[i] - b[i];
    return result;
}

double* addVectors(double* a, double* b, int size) {
    double* result = (double*)malloc(size * sizeof(double));
    for (int i = 0; i < size; i++)
        result[i] = a[i] + b[i];
    return result;
}

double* multVectByScalar(double* a, int size, double scalar) {
    double* result = (double*)malloc(size * sizeof(double));
    for (int i = 0; i < size; i++)
        result[i] = a[i] * scalar;
    return result;
}

double* multMatrixByVector(double* matrix, double* vector, int size) {
    double* result = (double*)malloc(size * sizeof(double));
    for (int i = 0; i < size; ++i)
        result[i] = dotProduct(&matrix[i * size], vector, size);
    return result;
}

double* multMatrixByMatrixTransposed(double* A, int aRows, int aCols, double* B,  int bRows, int bCols){
    if(aCols != bRows){
        fprintf(stderr, "Wrong matrices, unable to multiply\n");
        fflush(stderr);
        return NULL;
    }
    double* resultMatrix = initMatrix(aRows, bCols);
    fillMatrix(resultMatrix, aRows, bCols, 0);
    for(int i = 0; i < aRows; ++i)
        for(int j = 0; j < bCols; ++j)
            resultMatrix[i * bCols + j] += dotProduct(&A[i * aCols], &B[j * bRows], bRows);
    return resultMatrix;
}

double* multMatrixByMatrix(double* A, int aRows, int aCols, double* B,  int bRows, int bCols){
    if(aCols != bRows){
        fprintf(stderr, "Wrong matrices, unable to multiply\n");
        fflush(stderr);
        return NULL;
    }
    double* resultMatrix = initMatrix(aRows, bCols);
    fillMatrix(resultMatrix, aRows, bCols, 0);

    for (int i = 0; i < aRows; ++i) 
        for (int k = 0; k < aCols; ++k) 
            for (int j = 0; j < bCols; ++j) 
                resultMatrix[i * bCols + j] += A[i * aCols + k] * B[k * bCols + j];
    return resultMatrix;
}