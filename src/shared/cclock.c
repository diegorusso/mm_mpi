#include <stdlib.h>
#include <sys/types.h>
#include <sys/time.h>

/*  -------------------------------------------------------------------

    This function returns the wall clock time with micro seconds
    accuracy.
    The data type of the returned value is "double".
    The function can be called from a FORTRAN module. The value
    returned by cclock_ and cclock should be of type REAL(Kind = 8).
    Used by mm-serial and mm-2D-cannon only

    -------------------------------------------------------------------
*/

double cclock( void )
{
   const  double  micro = 1.0e-06;    /* Conversion constant */
   static long    start = 0L, startu;
   struct timeval tp; 	              /* Structure used by gettimeofday */
   double         wall_time;          /* To hold the result */


   if ( gettimeofday( &tp, NULL) == -1 )
      wall_time = -1.0e0;
   else if( !start ) {
      start  = tp.tv_sec;
      startu = tp.tv_usec;
      wall_time = 0.0e0;
   }
   else
      wall_time = (double) (tp.tv_sec - start) + micro*(tp.tv_usec - startu);

   return wall_time;
}
