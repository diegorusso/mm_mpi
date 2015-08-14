#include <math.h>
#include <float.h>
#include <stdio.h>

int check_cannon( int m, int l, int n, double *c )
// ---------------------------------------------------------------------
// --- Routine 'check' does the correctness check for the matrix-matrix
//     multiplication in function 'mxm'. Because of the particular
//     data that are generatedin function 'gendat', this check can be
//     used.
// ---------------------------------------------------------------------
{
   double eps, tvalue = (double)l;
   int    i, j, ok = 0;
   int    ELEM;

   ELEM = (int)(m*n);

   eps = 2.0*l*l*DBL_EPSILON;
   for( i = 0; i < ELEM; i++ ) {
         if ( fabs( tvalue - c[i] ) > eps ) ok++;
   }
   return( ok );
}
