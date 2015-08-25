#include <math.h>
#include <float.h>
#include <stdio.h>


int check( int m, int l, int n, double c[][n] )
// ---------------------------------------------------------------------
// --- Routine 'check' does the correctness check for the matrix-matrix
//     multiplication in function 'mxm'. Because of the particular
//     data that are generatedin function 'gendat', this check can be
//     used.
// ---------------------------------------------------------------------
{
   double eps, tvalue = (double)l;
   int    i, j, ok = 0;
// ---------------------------------------------------------------------
   eps = 2.0*l*l*DBL_EPSILON;
   for( i = 0; i < m; i++ ) {
      for( j = 0; j < n; j++ ){
         if ( fabs( tvalue - c[i][j] ) > eps ) ok++;
      }
   }
   return( ok );
}
