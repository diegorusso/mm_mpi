#include <stdio.h>
#include <stdarg.h>
#include "shared.h"

// Function to print debug logs
void verbose_printf(const char* fmt, ...) {
    if (verbose) {
        va_list args;
        va_start(args, fmt);
        vfprintf(stderr, fmt, args);
        va_end(args);
    };
};

// The function prints a header
void header(void) {
    printf("mm: Matrix-matrix multiply test C(m,n) = A(m,l)*B(l,n)\n");
    printf("-------------------------------------------------------\n");
    printf("      Problem size     |            |            |    |\n");
    printf("   m   |   l   |   n   |  Time (s)  |  (Gflop/s) | OK?|\n");
    printf("-------------------------------------------------------\n");
}


// This functions calculates the speed of mm taking in account the data size
// and the time passed the it prints results at screen
void results(int m, int l, int n, double time, int ok) {
    double      numops, gflops;

    // Speed is given by size of the data divided by time passed
    numops = 2 * (double)m * (double)n * (double)l;  // 2 * size^3
    gflops = 1.0e-9 * numops / time;  //1.0e-9 is giga

    printf(" %6d| %6d| %6d| %11.4f| %11.4e| ", m, l, n, time, gflops);

    if (ok == 0)
        printf(" T |\n");
    else
        printf(" F |\n");

    fflush(stdout);
}
