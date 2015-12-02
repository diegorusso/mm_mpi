#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "format.h"

// This method generates matrix A and B. Matrices are just vector of doubles
// split in blocks
void gendat(int m, int l, int n, double *a, double *b) {
    int     i, j, nprocesses, N_DBLOCK, M_DBLOCK, L_DBLOCK;
    int     myrank, my2drank, mycoords[2];
    int     periods[2];
    int     dimensions[2];

    // MPI Comminucator
    MPI_Comm    comm_2d;

    // Determines the size of the group associated with a communicator: it
    // returns number of processes in the group MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &nprocesses);

    // Determines the rank of the calling process in the communicator: it
    // returns the rank of the calling process in the group MPI_COMM_WORLD
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // I have just two dimensions: they will be passed to build the MPI cart
    // communicator. Every dimension has the squre root of number of processors
    dimensions[0] = dimensions[1] = (int)sqrt(nprocesses);
    // I set the periodic attribut to false for every dimension
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
    MPI_Comm_rank(comm_2d, &my2drank);

    // Determines process coords in cartesian topology given rank in group
    // comm: communicator with cartesian structure
    // rank: rank of a process within group of comm
    // maxdims: length of vector coords in the calling program
    // coords: integer array (of size ndims ) containing the Cartesian
    //         coordinates of specified process
    MPI_Cart_coords(comm_2d, my2drank, 2, mycoords);

    verbose_printf("%i x %i, rank %i \n", mycoords[0], mycoords[1], my2drank);

    // I split the matrices dimensions in blocks using the square root of the
    // number of processors I have available
    M_DBLOCK = (int)(m / sqrt(nprocesses));
    L_DBLOCK = (int)(l / sqrt(nprocesses));
    N_DBLOCK = (int)(n / sqrt(nprocesses));

    verbose_printf("BLOCKS - M: %i L: %i N: %i \n", M_DBLOCK, L_DBLOCK, N_DBLOCK);

    // I generate data for matrix A (m x l) represented by an array
    if (myrank == 0) printf("Generating data for matrix A[%d x %d]\n", m, l);
    for (j = 0; j < M_DBLOCK; j++) {
        for (i = 0; i < L_DBLOCK; i++) {
            a[j * L_DBLOCK + i] = (double)((mycoords[1] * L_DBLOCK) + i + 1);
            verbose_printf("A[%i * %i + %i]: %f = %i * %i + %i + 1\n", j, L_DBLOCK, i, a[j * L_DBLOCK +i], mycoords[1], L_DBLOCK, i);
        }
    }

    // I generate data for matrix B (n x l) represented by an array
    if (myrank == 0) printf("Generating data for matrix B[%d x %d]\n", l, n);
    for (i = 0; i < N_DBLOCK; i++) {
        for (j = 0; j < L_DBLOCK; j++) {
            // Multiplicative inverse
            b[j * N_DBLOCK + i] = 1.0 / (double)((mycoords[0] * L_DBLOCK) + j + 1);
            verbose_printf("B[%i * %i + %i]: %f = 1.0 / %i * %i + %i + 1\n", j, N_DBLOCK, i, b[j * N_DBLOCK + i], mycoords[0], L_DBLOCK, j);
        }
    }

    // I don't need the communicator anymore
    MPI_Comm_free(&comm_2d);
}
