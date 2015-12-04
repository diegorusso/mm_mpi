#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include "mxm-local.h"


void mxm(int m, int l, int n, double *a, double *b, double *c, MPI_Comm comm) {
	int		i         , npes, M_DBLOCK, N_DBLOCK, L_DBLOCK, A_DBLOCK,
			B_DBLOCK;
	int		myrank    , my2drank, mycoords[2];
	int		source_a  , destination_a, source_b, destination_b;
	int		Periods    [2], Dimensions[2], Coordinates[2], Remain_dims[2];

	double         *a_buf[2], *b_buf[2];

	MPI_Status	status;
	MPI_Status	req_handlers_status[4];
	MPI_Request	req_handlers[4];

	MPI_Comm	comm_2d, Row_comm, Col_comm;


	MPI_Comm_size(comm, &npes);
	MPI_Comm_rank(comm, &myrank);

	Dimensions[0] = Dimensions[1] = (int)sqrt(npes);
	Periods[0] = Periods[1] = 1;

	MPI_Cart_create(comm, 2, Dimensions, Periods, 1, &comm_2d);
	MPI_Comm_rank(comm_2d, &my2drank);
	MPI_Cart_coords(comm_2d, my2drank, 2, mycoords);

	Remain_dims[0] = 1;
	Remain_dims[1] = 0;
	MPI_Cart_sub(comm_2d, Remain_dims, &Col_comm);

	Remain_dims[0] = 0;
	Remain_dims[1] = 1;
	MPI_Cart_sub(comm_2d, Remain_dims, &Row_comm);

	M_DBLOCK = (int)(m / sqrt(npes));
	L_DBLOCK = (int)(l / sqrt(npes));
	N_DBLOCK = (int)(n / sqrt(npes));

	A_DBLOCK = M_DBLOCK * L_DBLOCK;
	B_DBLOCK = L_DBLOCK * N_DBLOCK;

	a_buf[0] = a;
	b_buf[0] = b;
	a_buf[1] = (double *)malloc(A_DBLOCK * sizeof(double));
	b_buf[1] = (double *)malloc(B_DBLOCK * sizeof(double));

	if (mycoords[0] != 0) {
		source_a = (mycoords[0] + mycoords[1]) % Dimensions[0];
		destination_a = (mycoords[1] + Dimensions[0] - mycoords[0]) % Dimensions[0];
		MPI_Sendrecv_replace(a_buf[0], A_DBLOCK, MPI_DOUBLE, destination_a, 0, source_a, 0, Row_comm, &status);
	}
	if (mycoords[1] != 0) {
		source_b = (mycoords[0] + mycoords[1]) % Dimensions[1];
		destination_b = (mycoords[0] + Dimensions[1] - mycoords[1]) % Dimensions[1];
		MPI_Sendrecv_replace(b_buf[0], B_DBLOCK, MPI_DOUBLE, destination_b, 0, source_b, 0, Col_comm, &status);
	}
	for (i = 0; i < Dimensions[0]; i++) {

		source_a = (mycoords[1] + 1) % Dimensions[0];;
		destination_a = (mycoords[1] + Dimensions[0] - 1) % Dimensions[0];

		source_b = (mycoords[0] + 1) % Dimensions[1];
		destination_b = (mycoords[0] + Dimensions[1] - 1) % Dimensions[1];

#if defined __PREPOSTED_NONBLOCKING
		MPI_Irecv(a_buf[(i + 1) % 2], A_DBLOCK, MPI_DOUBLE, source_a, 303, Row_comm, &req_handlers[2]);
		MPI_Irecv(b_buf[(i + 1) % 2], B_DBLOCK, MPI_DOUBLE, source_b, 303, Col_comm, &req_handlers[3]);

		mxm_local(M_DBLOCK, L_DBLOCK, N_DBLOCK, a_buf[i % 2], b_buf[i % 2], c);

		MPI_Isend(a_buf[i % 2], A_DBLOCK, MPI_DOUBLE, destination_a, 303, Row_comm, &req_handlers[0]);
		MPI_Isend(b_buf[i % 2], B_DBLOCK, MPI_DOUBLE, destination_b, 303, Col_comm, &req_handlers[1]);

#else

		MPI_Isend(a_buf[i % 2], A_DBLOCK, MPI_DOUBLE, destination_a, 303, Row_comm, &req_handlers[0]);
		MPI_Isend(b_buf[i % 2], B_DBLOCK, MPI_DOUBLE, destination_b, 303, Col_comm, &req_handlers[1]);

		MPI_Irecv(a_buf[(i + 1) % 2], A_DBLOCK, MPI_DOUBLE, source_a, 303, Row_comm, &req_handlers[2]);
		MPI_Irecv(b_buf[(i + 1) % 2], B_DBLOCK, MPI_DOUBLE, source_b, 303, Col_comm, &req_handlers[3]);

		mxm_local(M_DBLOCK, L_DBLOCK, N_DBLOCK, a_buf[i % 2], b_buf[i % 2], c);
#endif
		MPI_Waitall(4, req_handlers, req_handlers_status);
	}


	/* Rearange step is MANDATORY */
	if (mycoords[0] != 0) {
		destination_a = (mycoords[0] + mycoords[1]) % Dimensions[0];
		source_a = (mycoords[1] + Dimensions[0] - mycoords[0]) % Dimensions[0];
		MPI_Sendrecv_replace(a_buf[0], A_DBLOCK, MPI_DOUBLE, destination_a, 0, source_a, 0, Row_comm, &status);
	}
	if (mycoords[1] != 0) {
		destination_b = (mycoords[0] + mycoords[1]) % Dimensions[1];
		source_b = (mycoords[0] + Dimensions[1] - mycoords[1]) % Dimensions[1];
		MPI_Sendrecv_replace(b_buf[0], B_DBLOCK, MPI_DOUBLE, destination_b, 0, source_b, 0, Col_comm, &status);
	}
	free(a_buf[1]);
	free(b_buf[1]);

	MPI_Comm_free(&Col_comm);
	MPI_Comm_free(&Row_comm);
	MPI_Comm_free(&comm_2d);
}
