#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include "format.h"
#include "mxm-local.h"

// This function performs the actual multiplication between matrix A and B
// storing data in C. The function accepts in input the dimensions, three
// pointers for the matrices and the MPI communicator (MPI_COMM_WORLD)
void mxm(int m, int l, int n, double *a, double *b, double *c, MPI_Comm comm) {
    int     i;
    int     num_processes;
    int     M_DBLOCK, N_DBLOCK, L_DBLOCK, A_DBLOCK, B_DBLOCK;
    int     rank, rank_2d;
    int     source, destination;
    int     up, down, left, right;
    int     periods[2], dimensions[2], coordinates[2];

#if defined NONBLOCKING
    // I need these two buffers. Those are arrays of double pointers
    double  *a_buf[2], *b_buf[2];

    // Array of requests
    MPI_Request request_handles[4];

    // Array of statuses
    MPI_Status  status_handles[4];
#endif

    // That's just a MPI status
    MPI_Status  status;

    // Define 2D communicator. Not to forget there is also comm, which is the
    // world one
    MPI_Comm    comm_2d;

    // Determines the size of the group associated with a communicator: it
    // returns number of processes in the group MPI_COMM_WORLD
    MPI_Comm_size(comm, &num_processes);

    // Determines the rank of the calling process in the communicator: it
    // returns the rank of the calling process in the group MPI_COMM_WORLD
    MPI_Comm_rank(comm, &rank);

    // I have just two dimensions: they will be passed to build the MPI cart
    // communicator. Every dimension has the squre root of number of processors
    dimensions[0] = dimensions[1] = (int)sqrt(num_processes);

    // I set the periodic attribute to true for every dimension
    periods[0] = periods[1] = 1;

    // It creates a new communicator with the following topology information:
    // comm_old: MPI_COMM_WORLD
    // ndims: number of dimensions of cartesian grid
    // dims: interger array of size ndims specifying the number of process in
    //       each dimension
    // periods: logical array of size ndims specifying whether the grid is
    //          periodic (true) or not (false) in each dimension
    // reorder: ranking may be reordered (true) or not (false)
    // comm_new: output cartesian communicator
    MPI_Cart_create(comm, 2, dimensions, periods, 1, &comm_2d);

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

    // Once I have the block dimensions I calculate the blocks for matrix A and B
    A_DBLOCK = M_DBLOCK * L_DBLOCK;
    B_DBLOCK = L_DBLOCK * N_DBLOCK;

#if defined NONBLOCKING
    a_buf[0] = a;
    a_buf[1] = (double *)malloc(A_DBLOCK * sizeof(double));

    b_buf[0] = b;
    b_buf[1] = (double *)malloc(B_DBLOCK * sizeof(double));
#endif

    // Returns the shifted source and destination ranks, given a  shift
    // direction and amount
    // comm: communicator with cartesian structure
    // direction: coordinate dimension of shift. The direction argument is in
    //            the range [0,n-1] for an n-dimensional Cartesian mesh.
    // disp: displacement (> 0: upwards shift, < 0: downwards shift)
    // rank_source: rank of source process
    // rank_dest: rank of destination process
    // With these two commands I get the ranks of the left and up shifts
    MPI_Cart_shift(comm_2d, 1, -1, &right, &left);
    MPI_Cart_shift(comm_2d, 0, -1, &down, &up);
    debug_printf(__func__, rank, "Left rank: %i, Up rank: %i\n", left, up);

    // This is the initial alignment of A
    MPI_Cart_shift(comm_2d, 1, -coordinates[0], &source, &destination);
    debug_printf(__func__, rank,
                 "Initial matrix alignment for A: shit of %i to %i\n",
                 -coordinates[0], destination);

#if defined NONBLOCKING
    // Sends and receives using a single buffer
    // buf: initial address of send and receive buffer
    // count: number of elements in send and receive buffer
    // datatype: type of elements in send and receive buffer
    // dest: rank of destination
    // sendtag: send message tag
    // source: rank of source
    // recvtag: receive message tag
    // comm: communicator
    // status: status object
    MPI_Sendrecv_replace(a_buf[0], A_DBLOCK, MPI_DOUBLE, destination, 0,
                         source, 0, comm_2d, &status);
#else
    MPI_Sendrecv_replace(a, A_DBLOCK, MPI_DOUBLE, destination, 0, source, 0,
                         comm_2d, &status);
#endif

    // This is the initial alignment of B
    MPI_Cart_shift(comm_2d, 0, -coordinates[1], &source, &destination);
    debug_printf(__func__, rank,
                 "Initial matrix alignment for B: shit of %i to %i\n",
                 -coordinates[1], destination);

#if defined NONBLOCKING
    MPI_Sendrecv_replace(b_buf[0], B_DBLOCK, MPI_DOUBLE, destination, 0,
                         source, 0, comm_2d, &status);
#else
    MPI_Sendrecv_replace(b, B_DBLOCK, MPI_DOUBLE, destination, 0, source, 0,
                         comm_2d, &status);
#endif

    debug_printf(__func__, rank,
                 "Iterate through %i dimensions\n", dimensions[0]);
    for (i = 0; i < dimensions[0]; i++) {
#if defined NONBLOCKING
        // Perform the local matrix multiplication
        mxm_local(M_DBLOCK, L_DBLOCK, N_DBLOCK, a_buf[i % 2], b_buf[i % 2], c);

        // Shift matrix A left by one
        // Begins a nonblocking send
        // buf: initial address of send buffer
        // count: number of elements in send buffer
        // datatype: datatype of each send buffer element
        // dest: rank of destination
        // tag: message tag
        // comm: communicator
        // request: communication request
        MPI_Isend(a_buf[i % 2], A_DBLOCK, MPI_DOUBLE, left, 1, comm_2d,
                  &request_handles[0]);
        // Begins a nonblocking receive
        // buf: initial address of receive buffer
        // count: number of elements in receive buffer
        // datatype: datatype of each send buffer element
        // dest: rank of source
        // tag: message tag
        // comm: communicator
        // request: communication request
        MPI_Irecv(a_buf[(i + 1) % 2], A_DBLOCK, MPI_DOUBLE, right, 1, comm_2d,
                  &request_handles[1]);
        debug_printf(__func__, rank,
                     "A: Sending data to %i, Receiving data from %i\n",
                     left, right);

        // Shift matrix B up by one
        MPI_Isend(b_buf[i % 2], B_DBLOCK, MPI_DOUBLE, up, 1, comm_2d,
                  &request_handles[2]);
        MPI_Irecv(b_buf[(i + 1) % 2], B_DBLOCK, MPI_DOUBLE, down, 1, comm_2d,
                  &request_handles[3]);
        debug_printf(__func__, rank,
                     "B: Sending data to %i, Receiving data from %i\n",
                     up, down);

        // Let's wait all the shifts/communications to happen
        // Waits for all given MPI Requests to complete
        // count: list length
        // array_of_requests: array of request handles
        // array_of_statuses: array of status objects
        MPI_Waitall(4, request_handles, status_handles);
#else
        // Perform the local matrix multiplication
        mxm_local(M_DBLOCK, L_DBLOCK, N_DBLOCK, a, b, c);

        // Shift matrix A left by one
        MPI_Sendrecv_replace(a, A_DBLOCK, MPI_DOUBLE, left, 1, right, 1,
                             comm_2d, &status);
        debug_printf(__func__, rank,
                     "A: Sending data to %i, Receiving data from %i\n",
                     left, right);

        // Shift matrix B up by one
        MPI_Sendrecv_replace(b, B_DBLOCK, MPI_DOUBLE, up, 1, down, 1, comm_2d,
                             &status);
        debug_printf(__func__, rank,
                     "B: Sending data to %i, Receiving data from %i\n",
                     up, down);
#endif
    }

    // Final matrix alignement for A
    MPI_Cart_shift(comm_2d, 1, +coordinates[0], &source, &destination);
    debug_printf(__func__, rank,
                 "Final matrix alignment for A: shit of %i to %i\n",
                 +coordinates[0], destination);

#if defined NONBLOCKING
    MPI_Sendrecv_replace(a_buf[0], A_DBLOCK, MPI_DOUBLE, destination, 0,
                         source, 0, comm_2d, &status);
#else
    MPI_Sendrecv_replace(a, A_DBLOCK, MPI_DOUBLE, destination, 0, source, 0,
                         comm_2d, &status);
#endif

    // Final matrix alignement for B
    MPI_Cart_shift(comm_2d, 0, +coordinates[1], &source, &destination);
    debug_printf(__func__, rank,
                 "Final matrix alignment for B: shit of %i to %i\n",
                 +coordinates[1], destination);

#if defined NONBLOCKING
    MPI_Sendrecv_replace(b_buf[0], B_DBLOCK, MPI_DOUBLE, destination, 0,
                         source, 0, comm_2d, &status);
#else
    MPI_Sendrecv_replace(b, B_DBLOCK, MPI_DOUBLE, destination, 0, source, 0,
                         comm_2d, &status);
#endif

#if defined NONBLOCKING
    free(a_buf[1]);
    free(b_buf[1]);
#endif

    // Free up MPI communicator
    MPI_Comm_free(&comm_2d);
}
