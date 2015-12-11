#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "format.h"

// This method generates matrix A and B. Matrices are just vector of doubles
// split in blocks
void gendat(int m, int l, int n, double *a, double *b) {
    int     i, j;
    int     num_processes;
    int     N_DBLOCK, M_DBLOCK, L_DBLOCK;
    int     rank, rank_2d;
    int     periods[2], dimensions[2], coordinates[2];

    // MPI Comminucator
    MPI_Comm    comm_2d;

    // Determines the size of the group associated with a communicator: it
    // returns number of processes in the group MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

    // Determines the rank of the calling process in the communicator: it
    // returns the rank of the calling process in the group MPI_COMM_WORLD
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // I have just two dimensions: they will be passed to build the MPI cart
    // communicator. Every dimension has the squre root of number of processors
    dimensions[0] = dimensions[1] = (int)sqrt(num_processes);

    // I set the periodic attribute to false for every dimension
    // I don't really need it here.
    periods[0] = periods[1] = 0;

    // It creates a new communicator with the following topology information:
    // comm_old: input communicator
    // ndims: number of dimensions of cartesian grid
    // dims: interger array of size ndims specifying the number of process in
    //       each dimension
    // periods: logical array of size ndims specifying whether the grid is
    //          periodic (true) or not (false) in each dimension
    // reorder: ranking may be reordered (true) or not (false)
    // comm_new: output cartesian communicator
    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, periods, 1, &comm_2d);

    // Determines the rank of the calling process in the communicator: it
    // returns the rank of the calling process in the group comm_2d
    MPI_Comm_rank(comm_2d, &rank_2d);

    // Determines process coords in cartesian topology given rank in group
    // comm: communicator with cartesian structure
    // rank: rank of a process within group of comm
    // maxdims: length of vector coords in the calling program
    // coords: integer array (of size ndims ) containing the Cartesian
    //         coordinates of specified process
    MPI_Cart_coords(comm_2d, rank_2d, 2, coordinates);

    debug_printf(__func__, rank_2d, "My coordinates are: %i, %i\n",
                 coordinates[0], coordinates[1]);

    // I split the matrices dimensions in blocks using the square root of the
    // number of processors I have available
    M_DBLOCK = (int)(m / sqrt(num_processes));
    L_DBLOCK = (int)(l / sqrt(num_processes));
    N_DBLOCK = (int)(n / sqrt(num_processes));

    debug_printf(__func__, rank_2d, "BLOCKS - M: %i L: %i N: %i \n",
                 M_DBLOCK, L_DBLOCK, N_DBLOCK);

    // I generate data for matrix A (m x l) represented by an array
    if (rank == 0) printf("Generating data for matrix A[%d x %d]\n", m, l);
    for (i = 0; i < M_DBLOCK; i++) {
        for (j = 0; j < L_DBLOCK; j++) {
            a[i * L_DBLOCK + j] = \
                (double)((coordinates[1] * L_DBLOCK) + j + 1);
            debug_printf(__func__, rank_2d,
                         "A[%i * %i + %i]: %i * %i + %i + 1 = %g\n",
                         i, L_DBLOCK,j, coordinates[1], L_DBLOCK, j,
                         a[i * L_DBLOCK + j] );
        }
    }

    // I generate data for matrix B (n x l) represented by an array
    if (rank == 0) printf("Generating data for matrix B[%d x %d]\n", l, n);
    for (i = 0; i < L_DBLOCK; i++) {
        for (j = 0; j < N_DBLOCK; j++) {
            // Multiplicative inverse
            b[i * N_DBLOCK + j] = \
                1.0 / (double)((coordinates[0] * L_DBLOCK) + i + 1);
            debug_printf(__func__, rank_2d,
                         "B[%i * %i + %i]: 1.0 / ((%i * %i) + %i + 1) = %g\n",
                         i, N_DBLOCK, j,  coordinates[0], L_DBLOCK, i,
                         b[i * N_DBLOCK + j]);
        }
    }

    // I don't need the communicator anymore
    MPI_Comm_free(&comm_2d);
}
