#pragma once

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

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

double vectorNorm(double* a, int size) {
    double norm = 0;
    for (int i = 0; i < size; ++i)
        norm += a[i] * a[i];
    return sqrt(norm);
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

double* multMatrixByVector(double* matrix, int matrixRowsCount, double* vector, int vectorSize) {
    double* result = (double*)malloc(matrixRowsCount * sizeof(double));
    for (int i = 0; i < matrixRowsCount; ++i)
        result[i] = dotProduct(&matrix[i * vectorSize], vector, vectorSize);
    return result;
}

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

void fillMatrixRandom(double* martix, int size) {
    srand(size);
    for (int i = 0; i < size; ++i)
        martix[i] = (double)rand() / RAND_MAX * 2.0 - 1.0;
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

void addToDiagonal(double* matrix, int size, double value) {
    for (int i = 0; i < size; ++i)
        matrix[i * size + i] += value;
}

void setMatrixCell(double* matrix, int height, int width, int i, int j, double value) {
    matrix[i * height + j] = value;
}

void setVectorCell(double* vector, int size, int i, double value) {
    vector[i] = value;
}

void printVector(double* vector, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            fprintf(stderr, "%lf ", vector[i * rows + j]);
        }
        fprintf(stderr, "\n");
    }
}

void freeMatrix(double* matrix) {
    free(matrix);
}

void freeVector(double* vector) {
    free(vector);
}
