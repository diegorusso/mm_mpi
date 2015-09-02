#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

double		cclock  (void);
int		check      (int, int, int, double *);
void		gendat    (int, int, int, double *, double *);
void		matrix_reset(int, int, double *);
void		MatrixMatrixMultiply(int, int, int, double *, double *, double *, MPI_Comm);
void		prthead   (void);
void		prtspeed  (int, int, int, double, int);
void		state     (char *);

int 
main(int argc, char **argv) {
	int		m         , l, n, M_DBLOCK, N_DBLOCK, L_DBLOCK;
	int		ok        , okl, nrep;
	int		i         , me, nprocs, sqrt_nprocs;
	double         *a, *b, *c;
	double		time_start, time_stop, time_compute;
	FILE           *inl;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &me);

	sqrt_nprocs = (int)sqrt((double)nprocs);

	/* Check if 'nprocs' is a perfect square ... */
	if (nprocs != (sqrt_nprocs * sqrt_nprocs)) {
		if (me == 0) {
			printf("\t >> Number of Processors should be perfect square!!!\n");
			fflush(stdout);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return 0;
	}
	if (me == 0) {
		state("mm");
		prthead();
		fflush(stdout);
	}
	inl = fopen("mm.in", "r");
	while ((fscanf(inl, "%d%d%d%d\n", &m, &l, &n, &nrep) != EOF)) {

		M_DBLOCK = (int)(m / sqrt_nprocs);
		L_DBLOCK = (int)(l / sqrt_nprocs);
		N_DBLOCK = (int)(n / sqrt_nprocs);

		/* Check if '{M,L,N}_DBLOCK' are integer ... */
		if ((m != (M_DBLOCK * sqrt_nprocs)) || (n != (N_DBLOCK * sqrt_nprocs)) || (l != (L_DBLOCK * sqrt_nprocs))) {
			if (me == 0) {
				printf(" %6d| %6d| %6d| Wrong Dimensions (not int)   |\n", m, l, n);
				fflush(stdout);
			}
		} else {

			a = calloc(M_DBLOCK * L_DBLOCK, sizeof(double));
			b = calloc(L_DBLOCK * N_DBLOCK, sizeof(double));
			c = calloc(M_DBLOCK * N_DBLOCK, sizeof(double));

			gendat(m, l, n, a, b);

			MPI_Barrier(MPI_COMM_WORLD);

			for (i = 0; i < nrep; i++) {
				matrix_reset(M_DBLOCK, N_DBLOCK, c);

				time_start = MPI_Wtime();
				MatrixMatrixMultiply(m, l, n, a, b, c, MPI_COMM_WORLD);
				time_stop = MPI_Wtime();
				time_compute += (time_stop - time_start);
			}

			MPI_Barrier(MPI_COMM_WORLD);

			okl = check(M_DBLOCK, l, N_DBLOCK, c);
			MPI_Reduce(&okl, &ok, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);

			time_compute = time_compute / nrep;

			if (me == 0) {
				prtspeed(m, l, n, time_compute, ok);
			}
			free(a);
			free(b);
			free(c);

		}

	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
	return 0;
}
