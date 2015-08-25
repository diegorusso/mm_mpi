#include <math.h>
#include <mpi.h>

void matrix_reset( int m, int n, double *c )
{
   int    i, j;

   for( i = 0; i < n; i++ ){
      for( j = 0; j < m; j++ ){
         c[j*n+i] = 0.0;
      }
   }
}

void gendat( int m, int l, int n, double *a, double *b)
{
        int i, j, npes, N_DBLOCK, M_DBLOCK, L_DBLOCK;
        int myrank, my2drank, mycoords[2];
	int Periods[2];  
	int Dimensions[2];   

        MPI_Status status;
        MPI_Comm comm_2d;

        MPI_Comm_size(MPI_COMM_WORLD, &npes);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        Dimensions[0] = Dimensions[1] = (int)sqrt(npes);
        Periods[0] = Periods[1] = 0;
        
	MPI_Cart_create(MPI_COMM_WORLD, 2, Dimensions, Periods, 1, &comm_2d);
        MPI_Comm_rank(comm_2d, &my2drank);
        MPI_Cart_coords(comm_2d, my2drank, 2, mycoords);

	M_DBLOCK = (int)(m/sqrt(npes));
 	L_DBLOCK = (int)(l/sqrt(npes));
 	N_DBLOCK = (int)(n/sqrt(npes));

	for( j = 0; j < M_DBLOCK; j++ ){
		for( i = 0 ; i < L_DBLOCK; i++ ){
			a[j*L_DBLOCK+i] = (double)( (mycoords[1]*L_DBLOCK) + i + 1 ); 
		}
	}

	for( i = 0; i < N_DBLOCK; i++ ){
		for( j = 0; j < L_DBLOCK; j++ ){
			b[j*N_DBLOCK+i] = 1.0 / (double)( (mycoords[0]*L_DBLOCK ) + j + 1);
		}
	}

	MPI_Comm_free(&comm_2d);
}
