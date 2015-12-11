#include <stdio.h>
#include <stdlib.h>
#include "shared.h"
#include "format.h"
#include "reset.h"
#include "utils.h"

int         check(int m, int l, int n, double *c);
void        gendat(int, int, int, double *a, double *b);
void        mxm(int, int, int, double *a, double *b, double *c);

// By default berbose is disabled
// debug is a global variable (in shared.h)
int debug = 0;

int main(int argc, char **argv) {
    int     m, l, n;
    int     ok = 0;
    int     i, num_repititions;
    double  *a, *b, *c;
    double  time_begin, time_end, elapsed_time;
    FILE    *input_file;

    // Parse arguments and return the input file
    parse_arguments(argc, argv, &input_file);

    // Read data from the file
    while ((fscanf(input_file, "%d%d%d%d\n", &m, &l, &n, &num_repititions) != EOF)) {
        // calloc(): zero-initializes the buffer
        a = calloc(m * l, sizeof(double));
        b = calloc(l * n, sizeof(double));
        c = calloc(m * n, sizeof(double));

        // Let's generate some data
        // a, b, c are full matrix
        gendat(m, l, n, a, b);

        printf("Matrix C will be [%d x %d]\n", m, n);
        printf("Let's do the math (%d repetitions)\n", num_repititions);

        // Let's iterate through repititions
        for (i = 0; i < num_repititions; i++) {
            // reset the matrix containing results
            matrix_reset(m, n, c);
            // start the timer
            time_begin = cclock();
            // do the actual multiplication
            mxm(m, l, n, a, b, c);
            // stop the timer
            time_end = cclock();
            // Get the elapsed time
            elapsed_time += (time_end - time_begin);
            // check if everything was ok (0 means good)
            ok += check(m, l, n, c);
        }

        // Get the average of time spent on computing
        elapsed_time = elapsed_time / num_repititions;

        // Print some nice results
        header();
        results(m, l, n, elapsed_time, ok);
        footer();
        fflush(stdout);

        // Free the memory for the three matrices
        free(a);
        free(b);
        free(c);
    }

    return 0;
}
