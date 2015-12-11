#include <math.h>
#include <stdio.h>
#include "format.h"

// This method generates matrix 'A' and 'B' for the matrix-matrix
// multiplication 'C = AB'.
void gendat(int m, int l, int n, double a[][l], double b[][n]){
    int i, j;

    // That's for A
    printf("Generating data for matrix A[%d x %d]\n", m, l);
    for (i = 0; i < m; i++) {
        for (j = 0; j < l; j++) {
            // I just store j + 1
            a[i][j] = (double)(j + 1);
            debug_printf(__func__, -1, "a[%i][%i] = %i\n", i, j, j + 1);
        }
    }

    // That's for B
    printf("Generating data for matrix B[%d x %d]\n", l, n);
    for (i = 0; i < l; i++) {
        for (j = 0; j < n; j++) {
            // I store the multiplicative inverse of a[j][i]
            b[i][j] = 1.0 / (double)(i + 1);
            debug_printf(__func__, -1, "b[%i][%i] = %G\n", i, j, 1.0 / (double)(i + 1));
        }
    }
}
