#if defined CBLAS
#include "cblas.h"
#endif

// The method performs the actual multiplication between A and B storing data
// in C. This method will be executed on every MPI node with local data.
// It has few optimizations with OpenMP and CBlas
void mxm_local(int m, int l, int n, double *a, double *b, double *c) {
    int     i, j, k;
    double  t;

// The golden rule of OpenMP is that all variables defined in an outer scope
// are shared by default in the parallel region.
// All the read-only variables should be share
// The varibla I write should be private

#if defined OPENMP_OUTER
    // The parallel pragma starts a parallel block
    // By default data are shared
    // t, k, j are private, so just i is shared
    // the scheduler is automatic: the compiler will pick up the best one
#pragma omp parallel for default(shared) private(t,k,j) schedule(auto)
    for (i = 0; i < m; i++) {
        for (j = 0; j < l; j++) {
            t = a[i * l + j];
            for (k = 0; k < n; k++) {
                c[i * n + k] = c[i * n + k] + t * b[j * n + k];
            }
        }
    }

#elif defined OPENMP_MIDDLE
    for (i = 0; i < m; i++) {
        // t, j are private, so just i, k are shared
#pragma omp parallel for default(shared) private(t,j) schedule(auto)
        for (k = 0; k < n; k++) {
            t = c[i * n + k];  // t is private and local to every thread
            for (j = 0; j < l; j++) {
                t = t + a[i * l + j] * b[j * n + k];
            }
            c[i * n + k] = t;
        }
    }

#elif defined OPENMP_INNER
    for (i = 0; i < m; i++) {
        for (j = 0; j < l; j++) {
            t = a[i * l + j];
            // t is private, so just i, j, k are shared
            // t is initialized with the value that it has before the parallel
#pragma omp parallel for default(shared) firstprivate(t) schedule(auto)
            for (k = 0; k < n; k++) {
                c[i * n + k] = c[i * n + k] + t * b[j * n + k];
            }
        }
    }

#elif defined OPENMP_NESTED
    // Here I try to combine the two optimizations
#pragma omp parallel for default(shared) private(t,k,j) schedule(auto)
    for (i = 0; i < m; i++) {
#pragma omp parallel for default(shared) private(t,j) schedule(auto)
        for (k = 0; k < n; k++) {
            t = c[i * n + k];  // t is private and local to every thread
            for (j = 0; j < l; j++) {
                t = t + a[i * l + j] * b[j * n + k];
            }
            c[i * n + k] = t;
        }
    }

#elif defined CBLAS
    // Multiply two matrices using cblas dgemm library
    // CblasRowMajor: Indicates that the matrices are stored in row major
    //                order, with the elements of each row of the matrix stored
    //                contiguously
    // CblasNoTrans: Enumeration type indicating that the matrices A and B
    //               should not be transposed or conjugate transposed before
    //               multiplication
    // m, n, l: Integers indicating the size of the matrices:
    //          A[m x l], B[l x n], C[m x n]
    // alpha: Real value used to scale the product of matrices A and B
    // a: matrix A
    // l: number of columns of A
    // b: matrix B
    // n: number of columns of B
    // alpha: Real value used to scale matrix C
    // c: matrix C
    // n: number of columns of C
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, l, 1.0, a, l, b, n, 1.0, c, n);

#else
    // This is the normal matrix multiplication without any optimization
    // Complexity here is O(n^3) because of three for loops
    for (i = 0; i < m; i++)
        for (k = 0; k < n; k++) {
            t = c[i * n + k];
            for (j = 0; j < l; j++)
                t = t + a[i * l + j] * b[j * n + k];
            c[i * n + k] = t;
        }
#endif
}
