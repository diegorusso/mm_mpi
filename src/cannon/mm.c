#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "shared.h"
#include "format.h"
#include "reset.h"
#include "utils.h"

int         check(int, int, int, double *);
void        gendat(int, int, int, double *, double *);
void        mxm(int, int, int, double *, double *, double *, MPI_Comm);

// By default berbose is disabled
// verbose is a global variable (in shared.h)
int verbose = 0;

int main(int argc, char **argv) {
    int     m, l, n, M_DBLOCK, N_DBLOCK, L_DBLOCK;
    int     ok, local_check, nrep;
    int     i, myrank, nprocesses, sqrt_nprocesses;
    double  *a, *b, *c;
    double  time_start, time_stop, time_compute;
    int     verbose = 0;
    FILE    *input_file;

    // Initialize the MPI execution environment
    MPI_Init(&argc, &argv);

    // Determines the size of the group associated with a communicator: it
    // returns number of processes in the group MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &nprocesses);

    // Determines the rank of the calling process in the communicator: it
    // returns the rank of the calling process in the group MPI_COMM_WORLD
    // myrank == 0 means it is the first processor
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // Take the square root of number of processes
    sqrt_nprocesses = (int)sqrt((double)nprocesses);

    // Check if 'nprocesses' is a perfect square
    if (nprocesses != (sqrt_nprocesses * sqrt_nprocesses)) {
        if (myrank == 0) {
            printf("ERROR: Number of Processors should be perfect square\n");
            fflush(stdout);
        }

        // Blocks until all processes in the communicator have reached this
        // routine
        MPI_Barrier(MPI_COMM_WORLD);

        // Terminates MPI execution environment
        MPI_Finalize();

        // Let's return something which is not 0
        return 1;
    }

    parse_arguments(argc, argv, &input_file);

    // Read data from the file
    while ((fscanf(input_file, "%d%d%d%d\n", &m, &l, &n, &nrep) != EOF)) {
        // Every block is the one size of the matrix split by the square root
        M_DBLOCK = (int)(m / sqrt_nprocesses);
        L_DBLOCK = (int)(l / sqrt_nprocesses);
        N_DBLOCK = (int)(n / sqrt_nprocesses);

        // This check ensures that every dimension (m, l, n) is full covered by
        // the number of available processors
        if ((m != (M_DBLOCK * sqrt_nprocesses)) ||
            (l != (L_DBLOCK * sqrt_nprocesses)) ||
            (n != (N_DBLOCK * sqrt_nprocesses))) {
            if (myrank == 0) {
                printf("ERROR: dimensions are not covered properly\n");
                printf("m: %i = %i * %i\n", m, M_DBLOCK, sqrt_nprocesses);
                printf("l: %i = %i * %i\n", l, L_DBLOCK, sqrt_nprocesses);
                printf("n: %i = %i * %i\n", n, N_DBLOCK, sqrt_nprocesses);
                fflush(stdout);
            }
        } else {  // every looks good
            // calloc(): zero-initializes the buffer
            a = calloc(M_DBLOCK * L_DBLOCK, sizeof(double));
            b = calloc(L_DBLOCK * N_DBLOCK, sizeof(double));
            c = calloc(M_DBLOCK * N_DBLOCK, sizeof(double));

            // Let's generate some data using MPI. Every node is respondible
            // to generate its own data
            gendat(m, l, n, a, b);

            // Let's wait every node to generate the data
            MPI_Barrier(MPI_COMM_WORLD);

            if (myrank == 0) {
                printf("Matrix C will be [%d x %d]\n", m, n);
                printf("Let's do the math (%d repetitions)\n", nrep);
            }

            // Let's iterate through repititions
            for (i = 0; i < nrep; i++) {
                // reset the matrix containing results
                matrix_reset(M_DBLOCK, N_DBLOCK, c);
                // start the timer
                time_start = MPI_Wtime();
                // do the actual multiplication
                mxm(m, l, n, a, b, c, MPI_COMM_WORLD);
                // stop the timer
                time_stop = MPI_Wtime();
                // Get the elapsed time
                time_compute += (time_stop - time_start);
            }

            // Let's wait all the nodes have computed their calculation
            MPI_Barrier(MPI_COMM_WORLD);

            local_check = check(M_DBLOCK, l, N_DBLOCK, c);

            // Reduces values on all processes to a single value
            // sendbuf: address of send buffer (choice)
            // recvbuf: address of receive buffer (choice)
            // count: number of elements in send buffer (integer)
            // datatype: data type of elements of send buffer (handle)
            // op: reduce operation (handle)
            // root: rank of root process (integer)
            // comm: communicator (handle)
            // I need to sum all the local checks for every processor and see
            // if there was an error somewhere. I store the result to ok,
            // whish hosts the global sum
            MPI_Reduce(&local_check, &ok, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);

            // Get the average of time spent on computing
            time_compute = time_compute / nrep;

            if (myrank == 0) {
                // Print some nice results
                header();
                results(m, l, n, time_compute, ok);
                footer();
                fflush(stdout);
            }

            // Free the memory fo the three matrices
            free(a);
            free(b);
            free(c);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    // Terminates MPI execution environment
    MPI_Finalize();

    return 0;
}
