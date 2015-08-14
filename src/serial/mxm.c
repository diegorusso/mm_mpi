void mxm( int lda, int m, int l, int n, double a[][lda], double b[][n],
          double c[][n] )
// ---------------------------------------------------------------------
// --- Routine 'mvddot' does a matrix-vector multiplication 'Ab = c'
//     using an dotproduct implementation.
// ---------------------------------------------------------------------
{
   int    i, j, k;
   double t;
// ---------------------------------------------------------------------
   for( i = 0 ; i < m; i++ ){
      for( j = 0; j < l; j++ ){
         t = a[i][j];
         for( k = 0; k < n; k++ ) {
            c[i][k] += t*b[j][k];
         }
      }
   }
}
