#include <math.h>
#include <stdio.h>

// The function resets the matrix setting all the cells to 0.0
void matrix_reset(int m, int n, double *c) {
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            c[j * n + i] = 0.0;
        }
    }
}

// This method generates matrix 'A' and 'B' for the matrix-matrix
// multiplication 'C = AB'.
void gendat(int m, int l, int n, double a[][l], double b[][n]){
    int i, j;

    // That's for A
    printf("Generating data for matrix A[%d x %d]\n", m, l);
    for (j = 0; j < m; j++) {
        for (i = 0; i < l; i++) {
            // I just store i+1
            a[j][i] = (double)(i + 1);
        }
    }

    // That's for B
    printf("Generating data for matrix B[%d x %d]\n", l, n);
    for (i = 0; i < n; i++) {
        for (j = 0; j < l; j++) {
            // I store the multiplicative inverse of a[j][i]
            b[j][i] = 1.0 / (double)(j + 1);
        }
    }
}
