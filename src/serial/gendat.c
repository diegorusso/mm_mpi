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

void gendat_serial( int lda, int m, int l, int n, double a[][lda], double b[][n] )
// ---------------------------------------------------------------------
// --- Routine 'gendat' generates matrix 'A' and 'B' for the
//     matrix-matrix multiplication 'C = AB'.
// ---------------------------------------------------------------------
{
   int    i, j;
// ---------------------------------------------------------------------
   for( j = 0; j < m; j++ ){
      for( i = 0 ; i < l; i++ ){
         a[j][i] = (double)(i + 1);
      }
      a[j][lda-1] = 0.0;
   }
   for( i = 0; i < n; i++ ){
      for( j = 0; j < l; j++ ){
         b[j][i] = 1.0/(double)(j + 1);
      }
   }
}
