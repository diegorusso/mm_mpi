#include <stdio.h>
#include <stdlib.h>

int         check(int m, int l, int n, double *c);
void        gendat(int, int, int, double *a, double *b);
void        mxm(int, int, int, double *a, double *b, double *c);
// shared methods
void        matrix_reset(int, int, double *);
double      cclock(void);
void        prthead(void);
void        prtspeed(int, int, int, double, int);

int main(int argc, char **argv) {
    int     m, l, n;
    int     ok = 0;
    int     i, nrep;
    double  *a, *b, *c;
    double  time_start, time_stop, time_compute;
    FILE    *input_file;

    input_file = fopen("mm.in", "r");

    // Read data from the file
    while ((fscanf(input_file, "%d%d%d%d\n", &m, &l, &n, &nrep) != EOF)) {
        // calloc(): zero-initializes the buffer
        a = calloc(m * l, sizeof(double));
        b = calloc(l * n, sizeof(double));
        c = calloc(m * n, sizeof(double));

        // Let's generate some data
        gendat(m, l, n, a, b);

        printf("Matrix C will be [%d x %d]\n", m, n);
        printf("Let's do the math (%d repetitions)\n", nrep);
        for (i = 0; i < nrep; i++) {
            // reset the matrix containing results
            matrix_reset(m, n, c);
            // start the timer
            time_start = cclock();
            // do the actual multiplication
            mxm(m, l, n, a, b, c);
            // stop the timer
            time_stop = cclock();
            // Get the elapsed time
            time_compute += (time_stop - time_start);
            // check if everything was ok (0 means good)
            ok += check(m, l, n, c);
        }

        // Get the average of time spent on computing
        time_compute = time_compute / nrep;
        // Print some nice results
        prthead();
        prtspeed(m, l, n, time_compute, ok);
        fflush(stdout);

        // Free the memory fo the three matrices
        free(a);
        free(b);
        free(c);
    }
}
