#include <stdio.h>
#include <mpi.h>
#include <math.h>

#if defined __CBLAS
#include "cblas.h"
#endif

void mxms( int m, int l, int n, double *a, double *b, double *c )
{
   int    i, j, k;
   double t;

#if defined __OPENMP_OUTER

#pragma omp parallel for default(shared) private(t,k,j) schedule(runtime)    
   for( i = 0 ; i < m; i++ ){
      for( j = 0; j < l; j++ ){
         t = a[i*l+j];
         for( k = 0; k < n; k++ ){
             c[i*n+k] = c[i*n+k] + t*b[j*n+k];
         }
      }
   }

#elif defined __OPENMP_MIDDLE
   
   for( i = 0 ; i < m; i++ )
   {
#pragma omp parallel for default(shared) private(t,j) schedule(runtime)
      for( k = 0; k < n; k++ ) {
         t = c[i*n+k];
         for( j = 0; j < l; j++ ){
            t = t + a[i*l+j]*b[j*n+k];
         }
         c[i*n+k] = t;
      }
   }


#elif defined __OPENMP_INNER

   for( i = 0 ; i < m; i++ ){
     for( j = 0; j < l; j++ ){
         t = a[i*l+j];
#pragma omp parallel for default(shared) firstprivate(t) schedule(runtime)
         for( k = 0; k < n; k++ ) {
            c[i*n+k] = c[i*n+k] + t*b[j*n+k];
         }
      }
   }

#elif defined __OPENMP_NESTED

#pragma omp parallel for default(shared) private(t,k,j) num_threads(2) schedule(runtime)
   for( i = 0 ; i < m; i++ )
   {
#pragma omp parallel for default(shared) private(t,j) num_threads(2) schedule(runtime)
      for( k = 0; k < n; k++ ) {
         t = c[i*n+k];
         for( j = 0; j < l; j++ ){
            t = t + a[i*l+j]*b[j*n+k];
         }
         c[i*n+k] = t;
      }
   }


#else

   for( i = 0 ; i < m; i++ )
      for( k = 0; k < n; k++ )
      {
         t = c[i*n+k];
         for( j = 0; j < l; j++ )
            t = t + a[i*l+j]*b[j*n+k];
         c[i*n+k] = t;
      }

#endif
}

void MatrixMatrixMultiply(int m, int l, int n, double *a, double *b, double *c, MPI_Comm comm)
{
	int i, npes, M_DBLOCK, N_DBLOCK, L_DBLOCK, A_DBLOCK, B_DBLOCK;
	int myrank, my2drank, mycoords[2];
	int uprank, downrank, leftrank, rightrank, coords[2];
	int shiftsource, shiftdest,source, destination;
	int Periods[2], Dimensions[2], Coordinates[2], Remain_dims[2];

	MPI_Status status;
	MPI_Comm comm_2d, Row_comm, Col_comm;

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

	M_DBLOCK = (int)(m/sqrt(npes));
	L_DBLOCK = (int)(l/sqrt(npes));
	N_DBLOCK = (int)(n/sqrt(npes));

	A_DBLOCK = M_DBLOCK*L_DBLOCK;
	B_DBLOCK = L_DBLOCK*N_DBLOCK; 

	if( mycoords[0] != 0){
		source   = (mycoords[0] + mycoords[1]) % Dimensions[0];
		destination = (mycoords[1] + Dimensions[0] - mycoords[0]) % Dimensions[0]; 
		MPI_Sendrecv_replace(a, A_DBLOCK, MPI_DOUBLE, destination, 0, source, 0, Row_comm, &status);
	}

	if ( mycoords[1] != 0){
		source   = (mycoords[0] + mycoords[1]) % Dimensions[1];
		destination = (mycoords[0] + Dimensions[1] - mycoords[1]) % Dimensions[1]; 
		MPI_Sendrecv_replace(b, B_DBLOCK, MPI_DOUBLE,destination, 0, source, 0, Col_comm, &status);
	}

	for (i=0; i<Dimensions[0]; i++) {

#if defined __CBLAS
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M_DBLOCK, N_DBLOCK, L_DBLOCK, 1.0, a, L_DBLOCK, b, N_DBLOCK, 1.0, c, N_DBLOCK);
#else
		mxms(M_DBLOCK,L_DBLOCK,N_DBLOCK,a,b,c);
#endif

                source   = (mycoords[1] + 1 ) % Dimensions[0];;
                destination = (mycoords[1] + Dimensions[0] - 1) % Dimensions[0];
		MPI_Sendrecv_replace(a, A_DBLOCK, MPI_DOUBLE, destination, 0, source, 0, Row_comm, &status);

                source   = (mycoords[0] + 1) % Dimensions[1];
                destination = (mycoords[0] + Dimensions[1] - 1) % Dimensions[1];
		MPI_Sendrecv_replace(b, B_DBLOCK, MPI_DOUBLE, destination, 0, source, 0, Col_comm, &status);
	}


	/* Rearange step is MANDATORY */
	if( mycoords[0] != 0){
		destination  = (mycoords[0] + mycoords[1]) % Dimensions[0];
		source = (mycoords[1] + Dimensions[0] - mycoords[0]) % Dimensions[0];
		MPI_Sendrecv_replace(a, A_DBLOCK, MPI_DOUBLE, destination, 0, source, 0, Row_comm, &status);
	}

	if( mycoords[1] != 0){
		destination  = (mycoords[0] + mycoords[1]) % Dimensions[1];
		source = (mycoords[0] + Dimensions[1] - mycoords[1]) % Dimensions[1];
		MPI_Sendrecv_replace(b, B_DBLOCK, MPI_DOUBLE, destination, 0, source, 0, Col_comm, &status);
	}

	MPI_Comm_free(&Col_comm);
	MPI_Comm_free(&Row_comm);
	MPI_Comm_free(&comm_2d);
}
