#include <stdio.h>
#include <stdlib.h>

double cclock( void );
int    check_serial( int m, int l, int n, double *c );
void   matrix_reset( int , int , double * );
void   gendat_serial( int, int, int, int, double *a, double *b );
void   mvm( int, int, int, int, double *a, double *b, double *c );
void   prthead( void );
void   prtspeed( int, int, int, double, int, int );
void   state( char * );

int main()
{ 
   int      lda, m, l, n;
   int      ok, nops, nrep;
   int      i;
   double   *a, *b, *c;
   double   time_start, time_stop, time_compute;
   FILE     *inl;
// ------------------------------------------------------------------------
   state( "mm" );
   prthead();
   inl = fopen( "mm.in", "r" );
   while( ( fscanf( inl, "%d%d%d%d\n", &m, &l, &n, &nrep ) != EOF ) ){
      lda = l + 1;
      a = calloc( m*lda, sizeof( double ) );
      b = calloc( l*n, sizeof( double ) );
      c = calloc( m*n, sizeof( double ) );
      gendat_serial( lda, m, l, n, a, b );
      
      for( i = 0; i < nrep; i++ ){
          matrix_reset(m,n,c);
    	  time_start = cclock();
          mxm( lda, m, l, n, a, b, c );
          time_stop = cclock();
          time_compute += (time_stop - time_start);
      }
      
      ok   = check_serial( m, l, n, c );
      time_compute = time_compute/nrep;
      nops  = 2*m*l*n;
      prtspeed( m, l, n, time_compute, ok, nops );fflush(stdout);
      free( a ); free( b ); free( c);
   }
   printf( "-------------------------------------------------------\n" );
}
