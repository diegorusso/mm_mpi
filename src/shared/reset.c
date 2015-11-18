// The function resets the matrix setting all the cells to 0.0
void matrix_reset(int m, int n, double *c) {
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            c[j * n + i] = 0.0;
        }
    }
}
