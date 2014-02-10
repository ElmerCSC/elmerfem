#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* #include <malloc.h> */

#define MAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MIN(x,y) ( (x) > (y) ? (y) : (x) )
#define ABS(x) ( (x) > 0 ? (x) : ( -(x) ) )

#define AEPS 1.0e-15

#define ALLOCMEM( s ) malloc( s )
#define FREEMEM( p )  free( p )

#define A( i, j ) a[n * ( i ) + ( j )]

static void ludecomp( double *, int, int * );

void lu_mtrinv( double *a, int n )
{
    int i,j,k;
    int *pivot;

    double s;
  
    pivot = (int *)ALLOCMEM(n*sizeof(int));

    /*
     *  AP = LU
     */
    ludecomp( a,n,pivot );

    for( i=0; i<n; i++ )
    {
        if ( ABS(A(i,i))<AEPS )
        {
            fprintf( stderr, "Inv: Matrix is singular.\n" );
            return;
        }
        A(i,i) = 1.0/A(i,i);
    }

    /*  
     *  INV(U)
     */
    for( i=n-2; i>=0; i-- )
        for( j=n-1; j>i; j-- )
        {
            s = -A(i,j);
            for( k=i+1; k<j; k++ ) s -= A(i,k)*A(k,j);
            A(i,j) = s;
        }

    /*
     * INV(L)
     */
    for( i=n-2; i>=0; i-- )
        for( j=n-1; j>i; j-- )
        {
            s = 0.0;
            for( k=i+1; k<=j; k++ ) s -= A(j,k)*A(k,i);
            A(j,i) = A(i,i)*s;
        }
  
    /* 
     * A  = INV(AP)
     */
    for( i=0; i<n; i++ )
        for( j=0; j<n; j++ )
        {
            s = 0.0;
            for( k=MAX(i,j); k<n; k++ )
            {
                if ( k!=i )
                    s += A(i,k)*A(k,j);
                else
                    s += A(k,j);

                A(i,j) = s;
            }
        }

    /*
     * A = INV(A) (at last)
     */
    for( i=n-1; i>=0; i-- )
        if ( pivot[i]!=i )
            for( j=0; j<n; j++ )
            {
                s = A(i,j);
                A(i,j) = A(pivot[i],j);
                A(pivot[i],j) = s;
            }

    FREEMEM((char *)pivot);
}

/*
 * LU- decomposition by gaussian elimination. Row pivoting is used.
 * 
 * result:  AP = L'U;  L' = LD; pviot[i] is the swapped column number
 * for column i.
 *
 * Result is stored in place of original matrix.
 *
 */
static void ludecomp( double *a, int n, int *pivot )
{
    double swap;
    int i,j,k,l;

    for( i=0; i<n-1; i++ )
    {
        j = i;
        for( k=i+1; k<n; k++ ) if ( ABS(A(i,k)) > ABS(A(i,j)) ) j = k;

        if ( ABS(A(i,j)) < AEPS )
        {
            fprintf( stderr, "LUDecomp: Matrix is (at least almost) singular, %d %d.\n", i, j );
            return;
        }

        pivot[i] = j;
    
        if ( j != i )
            for( k=0; k<=i; k++ )
            {
                swap = A(k,j);
                A(k,j) = A(k,i);
                A(k,i) = swap;
            }

        for( k=i+1; k<n; k++ ) A(i,k) /= A(i,i);
    
        for( k=i+1; k<n; k++ )
        {
            if ( j != i )
            {
                swap = A(k,i);
                A(k,i) = A(k,j); 
                A(k,j) = swap;
            }

            for( l=i+1; l<n; l++ ) A(k,l) -= A(k,i) * A(i,l);
        }
    }
  
    pivot[n-1] = n-1;
    if ( ABS(A(n-1,n-1))<AEPS )
    {
        fprintf( stderr, "LUDecomp: Matrix is (at least almost) singular.\n" );
    }
}
