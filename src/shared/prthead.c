#include <stdio.h>

void prthead( void )
{
   printf( "mm: Matrix-matrix multiply test C(m,n) = A(m,l)*B(l,n)\n" );
   printf( "-------------------------------------------------------\n" );
   printf( "      Problem size     |            |            |    |\n" );
   printf( "   m   |   l   |   n   |  Time (s)  |  (Mflop/s) | OK?|\n" );
   printf( "-------------------------------------------------------\n" );
}
