// Classic (row times column) algorithm to multiply two matrices.
// Complexity is O(n^3)
void mxm(int m, int l, int n, double a[][l], double b[][n],
         double c[][n]) {
    int     i, j, k;

    for (i = 0; i < m; i++) {
        for (j = 0; j < l; j++) {
            for (k = 0; k < n; k++) {
                c[i][k] += a[i][j] * b[j][k];
            }
        }
    }
}
