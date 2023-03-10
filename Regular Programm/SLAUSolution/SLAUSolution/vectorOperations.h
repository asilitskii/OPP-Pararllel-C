#pragma once

double dotProduct(double* a, double* b, int size) {
    double result = 0;
    for (int i = 0; i < size; i++)
        result += a[i] * b[i];
    return result;
}

double* subtractVectors(double* a, double* b, int size) {
    double* result = malloc(size * sizeof(double));
    for (int i = 0; i < size; i++)
        result[i] = a[i] - b[i];
    return result;
}

double* addVectors(double* a, double* b, int size) {
    double* result = malloc(size * sizeof(double));
    for (int i = 0; i < size; i++)
        result[i] = a[i] + b[i];
    return result;
}

double* multVectByScalar(double* a, int size, double scalar) {
    double* result = malloc(size * sizeof(double));
    for (int i = 0; i < size; i++)
        result[i] = a[i] * scalar;
    return result;
}

double* multMatrixByVector(double** matrix, double* vector, int size) {
    double* result = malloc(size * sizeof(double));
    for (int i = 0; i < size; ++i)
        result[i] = dotProduct(matrix[i], vector, size);
    return result;
}

double** initMatrix(int size) {
    double** matrix = calloc(size, sizeof(double));
    for (int i = 0; i < size; ++i)
        matrix[i] = calloc(size, sizeof(double));
    return matrix;
}

double* initVector(int size) {
    double** vector = calloc(size, sizeof(double));
    return vector;
}

void freeMatrix(double** matrix, int size) {
    for (int i = 0; i < size; ++i)
        free(matrix[i]);
    free(matrix);
}

void freeVector(double* vector) {
    free(vector);
}

void fillVector(double* array, int size, double value) {
    for (int i = 0; i < size; ++i)
        array[i] = value;
}

void fillMatrix(double** matrix, int size, double value) {
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            matrix[i][j] = value;
}

void fillDiagonal(double** matrix, int size, double value) {
    for (int i = 0; i < size; ++i)
        matrix[i][i] = value;
}

void setMatrixCell(double** matrix, int size, int i, int j, double value) {
    //if (i >= size || j >= size) printf("ERROR out of bounds\n");
    matrix[i][j] = value;
}
