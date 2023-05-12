#pragma once

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

char* initMatrix(int height, int width) {
    char* matrix = (char*)calloc(height * width, sizeof(char));
    return matrix;
}

void fillGlider(char* map, int length){
    map[1] = map[length + 2] = map[2 * length] = map[2 * length + 1] = map[2 * length + 2] = 1;

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

void fillDiagonal(char* matrix, int size, char value) {
    for (int i = 0; i < size; ++i)
        matrix[i * size + i] = value;
}

void setMatrixCell(double* matrix, int height, int width, int i, int j, double value) {
    matrix[i * height + j] = value;
}

void setVectorCell(double* vector, int size, int i, double value) {
    vector[i] = value;
}

void printVector(char* vector, int rows, int cols) {
    /*for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("%lf ", vector[i * rows + j]);
        }
        printf("\n");
    }*/
    for(int i = 0; i < rows * cols; ++i){
        printf("%d ", vector[i]);
        if(i % cols == cols - 1)
            printf("\n");
    }
    printf("\n");
}

void freeMatrix(char* matrix) {
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

