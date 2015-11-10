#include <stdio.h>

#define max(a,b)( ((a) > (b)) ? (a) : (b) )

// This functions caulcates the speed of mm taking in account the data size and
// the time passed
void prtspeed(int m, int l, int n, double time, int ok) {
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
