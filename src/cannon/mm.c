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
// debug is a global variable (in shared.h)
int debug = 0;

int main(int argc, char **argv) {
    int     m, l, n, M_DBLOCK, N_DBLOCK, L_DBLOCK;
    int     ok = 0, local_check, num_repititions;
    int     i, rank, num_processes, sqrt_num_processes;
    double  *a, *b, *c;
    double  time_begin, time_end, elapsed_time;
    FILE    *input_file;

    // Initialize the MPI execution environment
    MPI_Init(&argc, &argv);

    // Determines the size of the group associated with a communicator: it
    // returns number of processes in the group MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

    // Determines the rank of the calling process in the communicator: it
    // returns the rank of the calling process in the group MPI_COMM_WORLD
    // rank == 0 means it is the first processor
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Take the square root of number of processes
    sqrt_num_processes = (int)sqrt((double)num_processes);

    // Check if 'num_processes' is a perfect square
    if (num_processes != (sqrt_num_processes * sqrt_num_processes)) {
        if (rank == 0) {
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

    // Parse arguments and return the input file
    parse_arguments(argc, argv, &input_file);

    if (rank == 0) {
        printf("Executing %s\n", argv[0]);
        printf("Using %i processors\n", num_processes);
        fflush(stdout);
    }

    // Read data from the file
    while ((fscanf(input_file, "%d%d%d%d\n", &m, &l, &n, &num_repititions) \
                != EOF)) {
        // Every block is one size of the matrix split by the square root
        M_DBLOCK = (int)(m / sqrt_num_processes);
        L_DBLOCK = (int)(l / sqrt_num_processes);
        N_DBLOCK = (int)(n / sqrt_num_processes);

        // This check ensures that every dimension (m, l, n) is full covered by
        // the number of available processors
        if ((m != (M_DBLOCK * sqrt_num_processes)) ||
            (l != (L_DBLOCK * sqrt_num_processes)) ||
            (n != (N_DBLOCK * sqrt_num_processes))) {
            if (rank == 0) {
                printf("ERROR: dimensions are not covered properly\n");
                printf("m: %i = %i * %i\n", m, M_DBLOCK, sqrt_num_processes);
                printf("l: %i = %i * %i\n", l, L_DBLOCK, sqrt_num_processes);
                printf("n: %i = %i * %i\n", n, N_DBLOCK, sqrt_num_processes);
                fflush(stdout);
            }
        } else { // Everything looks good, we're ready to rock!
            // calloc(): zero-initializes the buffer
            // a, b, c are submatrix
            a = calloc(M_DBLOCK * L_DBLOCK, sizeof(double));
            b = calloc(L_DBLOCK * N_DBLOCK, sizeof(double));
            c = calloc(M_DBLOCK * N_DBLOCK, sizeof(double));

            // Let's generate some data using MPI. Every node is respondible
            // to generate its own data
            gendat(m, l, n, a, b);

            // Let's wait every node to generate the data
            MPI_Barrier(MPI_COMM_WORLD);

            if (rank == 0) {
                printf("Matrix C will be [%d x %d]\n", m, n);
                printf("Let's do the math (%d repetitions)\n",
                       num_repititions);
            }

            // Let's iterate through repititions
            for (i = 0; i < num_repititions; i++) {
                // reset the matrix containing results
                matrix_reset(M_DBLOCK, N_DBLOCK, c);
                // start the timer
                time_begin = MPI_Wtime();
                // do the actual multiplication
                mxm(m, l, n, a, b, c, MPI_COMM_WORLD);
                // stop the timer
                time_end = MPI_Wtime();
                // Get the elapsed time
                elapsed_time += (time_end - time_begin);
            }

            // Let's wait all the nodes have computed their calculation
            MPI_Barrier(MPI_COMM_WORLD);

            // Once every node has finished, I'll perform a chack of the local
            // data
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
            MPI_Reduce(&local_check, &ok, 1, MPI_INTEGER, MPI_SUM, 0,
                       MPI_COMM_WORLD);

            // Get the average of time spent on computing
            elapsed_time = elapsed_time / num_repititions;

            // Print some nice results
            if (rank == 0) {
                header();
                results(m, l, n, elapsed_time, ok);
                footer();
                fflush(stdout);
            }

            // Free the memory for the three matrices
            free(a);
            free(b);
            free(c);
        }
    }

    // Let's wait any pending job
    MPI_Barrier(MPI_COMM_WORLD);

    // Terminates MPI execution environment
    MPI_Finalize();

    return 0;
}
