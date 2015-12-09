#include <math.h>
#include <stdio.h>
#include "format.h"

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
            verbose_printf(__func__, -1, "a[%i][%i] = %i + 1\n", j, i, i);
        }
    }

    // That's for B
    printf("Generating data for matrix B[%d x %d]\n", l, n);
    for (i = 0; i < n; i++) {
        for (j = 0; j < l; j++) {
            // I store the multiplicative inverse of a[j][i]
            b[j][i] = 1.0 / (double)(j + 1);
            verbose_printf(__func__, -1, "b[%i][%i] = 1.0 / (%i + 1)\n", j, i, j);
        }
    }
}
