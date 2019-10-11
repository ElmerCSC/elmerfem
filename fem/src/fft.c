/*  
   Elmer, A Finite Element Software for Multiphysical Problems
  
   Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
   
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
  
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with this library (in file ../../LGPL-2.1); if not, write 
   to the Free Software Foundation, Inc., 51 Franklin Street, 
   Fifth Floor, Boston, MA  02110-1301  USA
*/

/*
 * (Power of two) Fast Fourier Transform Subroutine Library
 *
 * Power of two algorithm to compute the transform:
 *
 * the basic idea is to compute partial f(k) sums for which exp(-i*2*pi*n*k/N)
 * is equal and which are also equal to different output frequencies. for
 * example for every (k) and (k+N/2) exp(-i*2*pi*n*k/N) is the same if (n) is
 * even. the algebra is as follows:
 *
 * N even,
 * F(n) = (1/N)*sum(f(k)*exp(-i*2*pi*n*k/N)), k=0..N-1,n=0..N-1
 * take
 * W = exp(-i*2*pi/N)
 * that is
 * F(n) = (1/N)*sum(f(k)*W^(n*k)), k=0..N-1, n=0..N-1
 *      = (1/N)*sum((f(k)*W^(k*n)+f(k+N/2)*W^((k+N/2)*n)), k=0..N/2-1,n=0..N-1
 *      = (1/N)*sum((f(k)+f(k+N/2)*W^(n*N/2))*W^(k*n)), k=0..N/2-1,n=0..N-1
 * now
 * W^(n*N/2) = exp(-i*2*pi*n*N/(2*N))
 *           = exp(-i*pi*n) = { 1: n even, -1: n odd }, n=0..N-1
 * take
 * h(k) = f(k)+f(l+N/2), k=0..N/2-1
 * g(k) = (f(k)-f(l+N/2)))*W^(k), k=0..N/2-1
 * then
 * 2*H(l) = (1/(N/2))*sum(h(k)*exp(-i*2*pi*l*k/(N/2)), l=0..N/2-1, (n=2*l)
 * 2*G(l) = (1/(N/2))*sum(g(k)*exp(-i*2*pi*l*k/(N/2)), l=0..N/2-1, (n=2*l+1)
 * and
 * F(2*l) = H(l), l=0..N/2-1
 * F(2*l+1) = G(l), l=0..N/2-1
 *
 * that is you can compute one size (N) discrete fourier transform as two
 * (M=N/2) size transforms. further, if (M) is even after each successive
 * division (i.e. (N) is power of two) you can continue the division until
 * you feel like writing the sum out explicitly.
 *
 * all the routines below at the end use cfftf (1D complex forward transform)
 * to actually compute the transform and do not in other ways depend on the
 * algorithm (which implies that input sequence lengths must be power of two,
 * serious limitation with storage for multidimensional transforms).so you can
 * substitute it for different algorithm if you like. for rfftf, rfftb, gfftf,
 * and gfftb (the real 1D transforms) the sequence lengths must be even.
 *
 * BTW. zero extend your input to next largest power of two when using the
 *      library. it does not check for it....
 *
 * The following routines are available:
 *
 *------------------------------------------------------------------------------
 *
 * cfftf: forward complex FFT.
 *
 * Parameters:
 *
 * N     : int            / length of sequences. 
 * T     : COMPLEX[N]     / input sequence.
 * F     : COMPLEX[N]     / transform sequence. may point to same memory as T.
 *
 * F(n) = sum(T(k)*exp(-i*2*pi*n*k/N)), k=0..N-1,n=0..N-1
 *
 *------------------------------------------------------------------------------
 *
 * cfftb: inverse complex FFT.
 *
 * Parameters:
 *
 * N     : int          / length of sequences.
 * F     : COMPLEX[N]   / transform sequence.
 * T     : COMPLEX[N]   / output sequence. may point to same memory as F.
 *
 * T(n) = sum(F(k)*exp(i*2*pi*n*k/N)), k=0..N-1,N=0..N-1
 *
 *------------------------------------------------------------------------------
 *
 * rfftf: forward real FFT. First (N/2+1) coefficients returned.
 *        F(n) = Real(F(N-n))-i*Imag(F(N-n)), n=N/2+1..N-1.
 *
 * Parameters:
 *
 * N     : int              / length of input sequence. 
 * T     : double[N]         / input sequence.
 * F     : COMPLEX[N/2+1]   / transform sequence. may point to same memory 
 *                          / as T (which must then be of length (N+2)).
 *
 * F(n) = sum(T(k)*exp(-i*2*pi*n*k/N)), k=0..N-1,n=0..N/2
 *
 * real transform as half size complex transform:
 *
 * f(k) real
 * F(n) = (1/N)*sum(f(k)*exp(-i*2*pi*n*k/N)), k=0..N-1,n=0..N-1
 * Divide to two parts:
 * g(k) = f(2*k),   k=0..N/2-1
 * h(k) = f(2*k+1), k=0..N/2-1
 * then
 * 2*G(n) = (1/(N/2))*sum(g(l)*exp(-i*2*pi*n*l/(N/2)),l=0..N/2-1,(k=2*l)
 * 2*H(n) = (1/(N/2))*sum(h(l)*exp(-i*2*pi*n*l/(N/2)),l=0..N/2-1,(k=2*l+1)
 * and 
 * F(n) = G(n) + exp(-i*2*pi*n/N)*H(n), n=0..N/2-1
 * Now take
 * y(k) = g(k) + i*h(k), k=0..N/2-1
 * then
 * Y(n) = (1/(N/2))*sum(y(k)*exp(-i*2*pi*n*k/(N/2))), k=0..N/2-1, n=0..N/2-1
 *      = G(n) + i*H(n)
 *      = G_even(n) + i*G_odd(n) + i*H_even(n) - H_odd(n)
 * where
 *       G_even is even part of G, G_odd odd part of G,
 *       H_even is even part of H, and H_odd odd part of H.
 * now 
 * 2*G_even(n) = Real( Y(n) + Y(N/2-n))
 * 2*G_odd(n)  = Imag( Y(n) - Y(N/2-n))
 * 2*H_even(n) = Imag( Y(n) + Y(N/2-n))
 * 2*H_odd(n)  = Real(-Y(n) + Y(N/2-n))
 * and (at last)
 * F(n)=G_even(n)+i*G_odd(n)+exp(-i*2*pi*n/N)*(H_even(n)+i*H_odd(n)), n=0..N/2-1
 * F(n)=Real(F(N-n))-i*Imag(F(N-n)), n=N/2+1..N-1
 * 
 *------------------------------------------------------------------------------
 *
 * rfftb: inverse real FFT. 
 *
 * Parameters:
 *
 * N     : int              / length of output sequence.
 * F     : COMPLEX[N/2+1]   / transform sequence.
 * T     : double[N]         / output sequence. may point to same memory as F.
 *
 * T(k) = sum(F(n)*exp(i*2*pi*k*n/N)), n=0..N-1,N=0..N-1
 *
 * Transform of complex input where the other half is complex conjugates of the
 * other (that is the real part is even, complex part odd, and the result is
 * real) as half size transform:
 *
 * F(N-n) = Real(F(n))-i*Imag(F(n))
 * f(k) = sum(F(n)*exp(i*2*pi*k*n/N)), n=0..N-1, k=0..N-1
 * Divide to two parts:
 * G(n) = F(2*n),            n=0..N/2-1
 * H(n) = F(2*n+1)-F(2*n-1), n=0..N/2-1
 * now
 * h(k) = (exp(i*2*pi*k/N)-exp(-i*2*pi*k/N))*x(k)
 *      = i*2*sin(2*pi*k/N)*x(k)
 * where 
 * X(n) = F(2*n+1), k=0..N/2-1
 * x(k) = sum(X(l)*exp(i*2*pi*k*l/(N/2))), l=0..N/2-1,(n=2*l+1)
 * g(k) = sum(G(l)*exp(i*2*pi*k*l/(N/2))), l=0..N/2-1,(n=2*l)
 * take 
 * Y(n) = G(n) + H(n)
 * so
 * y(k) = g(k) + h(k)
 *      = g(k) + x(k) * (i*2*sin(2*pi*k/N)), k=0..N/2-1
 * if
 * 2*g_even(k) = Real( y(k) + y(N/2-k))
 * 2*g_odd(k)  = Real(-y(k) + y(N/2-k))
 * 2*h_odd(k)  = Imag( y(k) - y(N/2-k))
 * 2*h_even(k) = Imag( y(k) + y(N/2-k))
 * k=0..N/2-1
 * then
 * f(k) = g(k) + x(k)
 *      = g_even(k)+g_odd(k) -
 *                    (h_odd(k)+h_evek(k))/(2*sin(2*pi*k/N)), k=0..N/2-1
 *      = g_even(N-k)-g_odd(N-k) -
 *                    (-h_odd(N-k)+h_even(N-k))/(2*sin(2*pi*k/N)), k=N/2+1..N-1
 *
 *------------------------------------------------------------------------------
 *
 * cfftf2D: 2D forward complex FFT.
 *
 * Parameters:
 *
 * M,N   : int              / array dimensions.
 * T     : COMPLEX[M][N]    / input array.
 * F     : COMPLEX[M][N]    / transform array. may point to same memory as T.
 *
 * F(m,n) = sum(T(k,l)*exp(-i*2*pi*(m*k/M+n*l/N))), k,m=0..M-1; l,n=0..N-1
 *
 *------------------------------------------------------------------------------
 *
 * cfftb2D: 2D inverse complex FFT.
 *
 * Parameters:
 *
 * M,N   : int              / array dimensions.
 * F     : COMPLEX[M][N]    / transform array.
 * T     : COMPLEX[M][N]    / output array. may point to same memory as T.
 *
 * T(m,n) = sum(F(k,l)*exp(i*2*pi*(m*k/M+n*l/N))), k,m=0..M-1; l,n=0..N-1
 *
 *------------------------------------------------------------------------------
 *
 * cfftf3D: 3D forward complex FFT.
 *
 * Parameters:
 *
 * L,M,N : int                / array dimensions.
 * T     : COMPLEX[L][M][N]   / input array.
 * F     : COMPLEX[L][M][N]   / transform array. may point to same memory as T.
 *
 * F(l,m,n) = sum(T(i,j,k)*exp(-i*2*pi*(l*i/L+m*j/M+n*k/N))), 
 *                                         l,i=0..L-1; m,j=0..M-1; k,n=0..N-1
 *
 *------------------------------------------------------------------------------
 *
 * cfftb3D: 3D inverse complex FFT.
 *
 * Parameters:
 *
 * L,M,N : int                / array dimensions.
 * F     : COMPLEX[L][M][N]   / transform array.
 * T     : COMPLEX[L][M][N]   / output array. may point to same memory as T.
 *
 * T(l,m,n) = sum(F(i,j,k)*exp(i*2*pi*(l*i/L+m*j/M+n*k/N))), 
 *                                         l,i=0..L-1; m,j=0..M-1; k,n=0..N-1
 *
 *------------------------------------------------------------------------------
 *
 * cfftfND: multidimensional (max 32) forward complex FFT.
 *
 * Parameters:
 *
 * N   : int                                   / number of dimensions.
 * D   : int[N]                                / array dimensions.
 * T   : COMPLEX T[D[N-1]][D[N-2]]...[D[0]]    / input array.
 * F   : COMPLEX F[D[N-1]][D[N-2]]...[D[0]]    / transform array.
 *                                             / may point to same memory as T.
 *
 *------------------------------------------------------------------------------
 *
 * cfftbND: multidimensional inverse complex FFT.
 *
 * Parameters:
 *
 * N   : int                                   / number of dimensions.
 * D   : int[N]                                / array dimensions.
 * F   : COMPLEX T[D[N-1]][D[N-2]]...[D[0]]    / transform array.
 * T   : COMPLEX F[D[N-1]][D[N-2]]...[D[0]]    / output array.
 *                                             / may point to same memory as T.
 *
 *------------------------------------------------------------------------------
 *
 * in all the routines above F and T may point to same memory, if the memory
 * is allocated according to larger of the arrays.
 *
 * The COMPLEX type is defined as:
 *
 * typedef
 * struct {
 *     double Real;
 *     double Imag;
 *     } COMPLEX;
 *
 *------------------------------------------------------------------------------
 *
 * gfftf( nT, time, nF, freq )  - return largest magnitude coeff. (real input)
 *     int nT;
 *     double time[nT + 2];
 *     int nF;
 *     FREQ freq[nF];
 *
 *------------------------------------------------------------------------------
 *
 * gfftb( nF, freq, nT, time )  - given largest magnitude coeff. return
 *     int nF;                    approximation of original sequence.
 *     FREQ freq[nF];
 *     int nT;
 *     double time[nT + 2];
 *
 *------------------------------------------------------------------------------
 *
 * The FREQ type is defined as:
 *
 * typedef
 * struct {
 *     double Real;
 *     double Imag;
 *     double Mag;
 *     int FBin;
 *     } FREQ;
 * 
 * The magnitude entry is returned by gfftf(), not used by gfftb().
 * 
 * Juha Ruokolainen / CSC
 * LAST CHANGE: 19.3. -92
 *
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "../config.h"

#define FALSE 0
#define TRUE  1

typedef
struct {
    double Real;
    double Imag;
    } COMPLEX;

typedef
struct {
    double Real;
    double Imag;
    double Mag;
    int FBin;
    } FREQ;

static double _FFT_PI = 3.14159265358979323844;

#define _FFT_MAX_LEVELS 30

static double _FFT_SIN[_FFT_MAX_LEVELS];
static double _FFT_COS[_FFT_MAX_LEVELS];

static int _FFT_I;
static int _FFT_Level;

/*
 * log2: base two logarithm of an (power of two!!!) integer.
 *       the macro needs and uses variable k.
 *
 * Parameters:
 *
 * R :  int     / put result here.
 * N :  int     / input number.
 *
 */
#define log2( R, N )                   \
   k = 1;                              \
   for( R = 0; R < 32; R++ ) {         \
       if ( k & N ) break;             \
       k <<= 1;                        \
       }

/*
 * BitReverse: reverse the order of given number of lowest bits in an integer.
 *             the macro needs and uses variables i & j.
 *
 * Parameters:
 *
 * R:       int               / put result here.
 * N:       int               / the input number.
 * Nbits:   int               / number of bits to be reversed - 1.
 *
 */
#define BitReverse( R, N, Nbits )      \
    R = 0;                             \
    j = 1;                             \
    for( i = 0; i <= Nbits; i++ ) {    \
        if ( j & N ) {                 \
            R |= 1 << ( Nbits - i );   \
            }                          \
        j <<= 1;                       \
        }

/*
 * BitReverseArray: clear the messed up order in an array,
 *                  the change is done in place.
 *
 * Parameters:
 *
 * N: int                 / number of elements in the array, and
 *                        / k=log(N) (the logarithm is base two),
 *                        / the number of bits.
 * T: COMPLEX[N]          / the array to be operated on.
 */
void BitReverseArray( N, T )
    int N;
    COMPLEX *T;
{
    COMPLEX swap;

    int i;
    int j;

    int k;
    int n;

    int logN;

    log2( logN, N );
    logN--;

    for( k = 0; k < N; k++ ) {
        BitReverse( n, k, logN );
        if ( n > k ) {
            swap = T[k];
            T[k] = T[n];
            T[n] = swap;
            }
       }
}

/*
 * sort: sort an (double) array to ascending order, and move the elements of
 *       another (integer) array accordingly. the latter can be used as track
 *       keeper of where an element in the sorted order at position (k) was in
 *       in the original order (Ord[k]), if it is initialized to contain
 *       numbers (0..N-1) before calling sort. 
 *
 * Parameters:
 *
 * N:      int                  / number of entries in the arrays.
 * Key:    double[N]             / array to be sorted.
 * Ord:    int[N]               / change this accordingly.
 */
void sort_swap(int i,int j,double *Key,int *Ord)
{
   int ival;
   double dval;

   dval=Key[i];
   Key[i]=Key[j];
   Key[j]=dval;

   ival=Ord[i];
   Ord[i]=Ord[j];
   Ord[j]=ival;
}

void sort_shift(int lbeg,int lend,double *Key,int *Ord)
{
   int i,j;

   i = lbeg;
   while( 2*i+1<=lend )
   {
      j=2*i+1;
      if ( j<lend && Key[j]<Key[j+1]) j=j+1;
      if ( Key[i]<Key[j]) {
         sort_swap(i,j,Key,Ord);
         i = j;
      } else break;
   }

}
void sort( N, Key, Ord )
    int N;
    int *Ord;
    double *Key;
{
  int lend,lbeg;

  lend  = N-1;
  lbeg = (lend-1)/2;

  while(lbeg>=0)
  {
     sort_shift( lbeg--,lend,Key,Ord );
  }

  while( lend>0) {
     sort_swap( 0,lend,Key,Ord );
     sort_shift( 0,--lend,Key,Ord );
  }
}

/*
 * FFTInit: initialize internal sin & cos tables for (pi / N)
 *          N being power of two up to _FFT_MAX_LEVELS.
 *
 * Parameters: NONE.
 */
static int FFTInit( )
{
    static int InitDone = FALSE;

    int k;
    int n;

    if ( InitDone ) {
        return 0;
        }

    n = ( 1 << _FFT_MAX_LEVELS );
    for( k = 0; k < _FFT_MAX_LEVELS; k++ ) {
        _FFT_COS[k] =  cos( _FFT_PI / n );    
        _FFT_SIN[k] = -sin( _FFT_PI / n );    
        n /= 2;
        }

    InitDone = TRUE;
}


/*
 * FFTKernel: discrete fourier transform. output is in bitreversed order.
 *
 * Parameters:
 *
 * N     : int          / length of sequences.
 * F     : COMPLEX[N]   / input sequence.
 * T     : COMPLEX[N]   / transform sequence.
 *                      / may point to same memory as F.
 *
 * T(n) = sum(F(k)*exp(-i*2*pi*n*k/N)), k=0..N-1,n=0..N-1
 *
 */
static int FFTKernel( N, F, T )
    int N;
    COMPLEX *F;
    COMPLEX *T;
{
    double ExpR;
    double ExpI;

    double CO;
    double SI;

    double TempR;
    double TempI;

    double t;

    int k;

    if ( N == 4 ) {
        TempR = F[0].Real; TempI = F[0].Imag;

        F[0].Real = TempR + F[2].Real;
        F[0].Imag = TempI + F[2].Imag;

        F[2].Real = TempR - F[2].Real;
        F[2].Imag = TempI - F[2].Imag;

        TempR = F[1].Real; TempI = F[1].Imag;

        F[1].Real = TempR + F[3].Real;
        F[1].Imag = TempI + F[3].Imag;

        TempR -= F[3].Real;
        TempI -= F[3].Imag;

        F[3].Real =  TempI;
        F[3].Imag = -TempR;

        TempR = F[0].Real; TempI = F[0].Imag;

        T[_FFT_I].Real = TempR + F[1].Real;
        T[_FFT_I].Imag = TempI + F[1].Imag;
        _FFT_I++;

        T[_FFT_I].Real = TempR - F[1].Real;
        T[_FFT_I].Imag = TempI - F[1].Imag;
        _FFT_I++;

        TempR = F[2].Real; TempI = F[2].Imag;

        T[_FFT_I].Real = TempR + F[3].Real;
        T[_FFT_I].Imag = TempI + F[3].Imag;
        _FFT_I++;

        T[_FFT_I].Real = TempR - F[3].Real;
        T[_FFT_I].Imag = TempI - F[3].Imag;
        _FFT_I++;

        return 0;
        }
   
    N /= 2;
 
    /*
     *  ExpR + i*ExpI = exp(-i*2*pi*k/N ) = exp(-i*2*pi/N)^(k)
     */
    ExpR = CO = _FFT_COS[_FFT_Level];
    ExpI = SI = _FFT_SIN[_FFT_Level];

    TempR = F[0].Real;
    TempI = F[0].Imag;

    F[0].Real = TempR + F[N].Real;
    F[0].Imag = TempI + F[N].Imag;

    F[N].Real = TempR - F[N].Real;
    F[N].Imag = TempI - F[N].Imag;

    for( k = 1; k < N; k++ ) {
        TempR = F[k].Real - F[k + N].Real;
        TempI = F[k].Imag - F[k + N].Imag;

        F[k].Real += F[k + N].Real;
        F[k].Imag += F[k + N].Imag;

        F[k + N].Real = TempR * ExpR - TempI * ExpI;
        F[k + N].Imag = TempR * ExpI + TempI * ExpR;

        t  = ExpR;
        ExpR = CO * t - SI * ExpI;
        ExpI = CO * ExpI + SI * t;
        }

    _FFT_Level++;

    FFTKernel( N, &F[0], T );
    FFTKernel( N, &F[N], T );

    _FFT_Level--;
}

/*
 * cfftf: forward complex FFT.
 *
 * Parameters:
 *
 * N     : int            / length of sequences. 
 * T     : COMPLEX[N]     / input sequence.
 * F     : COMPLEX[N]     / transform sequence. may point to same memory as T.
 *
 * F(n) = sum(T(k)*exp(-i*2*pi*n*k/N)), k=0..N-1,n=0..N-1
 */
void cfftf( N, T, F )
    int N;
    COMPLEX *T;
    COMPLEX *F;
{
    int logN;
    int k;

    FFTInit( );

    log2( logN, N );
    _FFT_Level = _FFT_MAX_LEVELS - logN + 1;

    _FFT_I = 0;

    if ( F != T ) {
        for( k = 0; k < N; k++ ) F[k] = T[k];
        }

    FFTKernel( N, F, F ); 
    BitReverseArray( N, F );
}

/*
 * cfftb: inverse complex FFT.
 *
 * Parameters:
 *
 * N     : int          / length of sequences.
 * F     : COMPLEX[N]   / transform sequence.
 * T     : COMPLEX[N]   / output sequence. may point to same memory as F.
 *
 * T(n) = sum(F(k)*exp(i*2*pi*n*k/N)), k=0..N-1,N=0..N-1
 *
 */
void cfftb( N, F, T )
    int N;
    COMPLEX *F;
    COMPLEX *T;
{
    int logN;
    int k;

    if ( F != T ) {
        for( k = 0; k < N; k++ ) T[k].Real = F[k].Real;
        }
    for( k = 0; k < N; k++ ) T[k].Imag = -F[k].Imag;
    cfftf( N, T, T );
    for( k = 0; k < N; k++ ) T[k].Imag = -T[k].Imag;
}

#ifdef USE_ISO_C_BINDINGS
void fcfftb( N, F, T )
    int *N;
    COMPLEX *F;
    COMPLEX *T;
{
    cfftb( *N, F, T );
}
#else
void FC_FUNC(fcfftb,FCFFTB)( N, F, T )
    int *N;
    COMPLEX *F;
    COMPLEX *T;
{
    cfftb( *N, F, T );
}
#endif /* USE_ISO_C_BINDINGS */

/*
 * rfftf: forward real FFT. First (N/2+1) coefficients returned.
 *        F(n) = Real(F(N-n))-i*Imag(F(N-n)), n=N/2+1..N-1.
 *
 * parameters:
 *
 * N     : int              / length of input sequence. 
 * T     : double[N]         / input sequence.
 * F     : COMPLEX[N/2+1]   / transform sequence. may point to same memory 
 *                          / as T (which must then be of length (N+2)).
 *
 * F(n) = sum(T(k)*exp(-i*2*pi*n*k/N)), k=0..N-1,n=0..N/2
 *
 */
void rfftf( N, T, F )
    int N;
    double   *T;
    COMPLEX *F;
{
    COMPLEX *W;

    double CO;
    double SI;

    double ExpR;
    double ExpI;

    double pi;
    double t;

    int k;

    N /= 2;

    W = (COMPLEX *)malloc( ( N + 1 ) * sizeof( COMPLEX ) );
/*
 *  if you like (or change the COMPLEX type) uncomment the following
 *
 * for( k = 0; k < N; k++ ) {
 *     W[k].Real = T[2*k];
 *     W[k].Imag = T[2*k+1];
 *     }
 *
 *  and replace 
 *     cfftf( N, (COMPLEX *)T, W );
 *  by
 *     cfftf( N, W, W );
 */
    cfftf( N, (COMPLEX *)T, W );
    W[N] = W[0];

    /*
     *  ExpR + i*ExpI = exp(-i*2*pi*k/N ) = exp(-i*2*pi/N)^(k)
     */
    ExpR = 1.0;
    ExpI = 0.0;

    pi = _FFT_PI / N;
    CO  =  cos( pi );
    SI  = -sin( pi );

    for( k = 0; k <= N; k++ ) {
        F[k].Real =  W[k].Imag + W[N - k].Imag;
        F[k].Imag = -W[k].Real + W[N - k].Real;

        t = F[k].Real;
        F[k].Real = ExpR * t - ExpI * F[k].Imag;
        F[k].Imag = ExpR * F[k].Imag + ExpI* t;

        F[k].Real += W[k].Real + W[N - k].Real;
        F[k].Imag += W[k].Imag - W[N - k].Imag;

        F[k].Real /= 2; F[k].Imag /= 2;

        t = ExpR;                   
        ExpR = CO * t - SI * ExpI;
        ExpI = CO * ExpI + SI * t;
        }

    free( W );
}

/*
 * rfftb: inverse real FFT. 
 *
 * Parameters:
 *
 * N     : int              / length of output sequence.
 * F     : COMPLEX[N/2+1]   / transform sequence.
 * T     : double[N]         / output sequence. may point to same memory as F.
 *
 * T(n) = sum(F(k)*exp(i*2*pi*n*k/N)), k=0..N-1,N=0..N-1
 *
 */
void rfftb( N, F, T )
    int N;
    COMPLEX *F;
    double   *T;
{
    COMPLEX *W;

    double ExpR;
    double ExpI;

    double CO;
    double SI;

    double pi;

    double t;

    double Eves;
    double Odds;

    int k;
    int n;

    N /= 2;

    W = (COMPLEX *)malloc( ( N + 1 ) * sizeof( COMPLEX ) );

    W[0].Imag = F[0].Imag + 2 * F[1].Imag;
    W[0].Real = F[0].Real;

    W[N / 2].Real = F[N].Real;
    W[N / 2].Imag = F[N].Imag - 2 * F[N - 1].Imag;

    for( k = 1; k < N / 2; k++ ) {
        n = 2 * k;
        W[k].Real = F[n].Real + F[n + 1].Real - F[n - 1].Real;
        W[k].Imag = F[n].Imag + F[n + 1].Imag - F[n - 1].Imag;
        }

    for( k = N / 2 + 1; k < N; k++ ) {
        n = 2 * ( N - k );
        W[k].Real =    F[n].Real + F[n - 1].Real - F[n + 1].Real;
        W[k].Imag = -( F[n].Imag + F[n - 1].Imag - F[n + 1].Imag );
        }

    Odds = F[1].Real;
    Eves = 0.0;
    for( k = 1; k < N / 2; k++ ) {
        n = 2 * k;
        Eves += F[n].Real;
        Odds += F[n + 1].Real;
        }
    Odds *= 2;
    Eves *= 2;
    Eves += F[0].Real + F[N].Real;

    cfftb( N, W, W );
    W[N] = W[0];

    /*
     *  ExpR + i*ExpI = exp(i*2*pi*k/N ) = exp(i*2*pi/N)^(k)
     */
    ExpR = 1.0;
    ExpI = 0.0;

    pi = _FFT_PI / N;
    CO  = cos( pi );
    SI  = sin( pi );

    for( k = 1; k < N; k++ ) {
        t = ExpR;                   
        ExpR = CO * t - SI * ExpI;
        ExpI = CO * ExpI + SI * t;

        n = 2 * N - k;
        T[n] = T[k] = 0.5 / ExpI;

        T[k] *= -W[k].Imag;
        T[k] +=  W[k].Real;

        T[n] *= W[N - k].Imag;
        T[n] += W[N - k].Real;
        }
    T[0] = Eves + Odds;
    T[N] = Eves - Odds;

    free( W );
}

#ifdef USE_ISO_C_BINDINGS
void frfftb( N, F, T )
    int *N;
    COMPLEX *F;
    double   *T;
{
   rfftb( *N, F, T );
}
#else
FC_FUNC(frfftb,FRFFTB)( N, F, T )
    int *N;
    COMPLEX *F;
    double   *T;
{
   rfftb( *N, F, T ); 
}
#endif /* USE_ISO_C_BINDINGS */

/*
 *
 * gfftb: given group of frequency coefficients (usually those with largest
 *        magnitude), compute approximation of the original sequence.
 *
 * Parameters:
 *
 * nF   : int              / number of input coefficients.
 * freq : FREQ[nF]         / input coefficients.
 * nT   : int              / length of output sequence.
 * time : time[nT+2]       / output sequence.
 *
 *   FREQ :
 *
 *   struct {
 *       double Real;
 *       double Imag;
 *       double Mag;
 *       int FBin;
 *       }
 *
 */
void gfftb( nF, freq, nT, time )
    int nT;
    int nF;
    double *time;
    FREQ  *freq;
{
    COMPLEX *F;

    int k;

    /* 
     * if you change the COMPLEX type replace 
     *     F = (COMPLEX *)time;
     *     bzero( F, ( nT / 2 + 1 ) * sizeof( COMPLEX ) );
     * by
     *     F = (COMPLEX *)calloc( ( nT / 2 + 1, sizeof( COMPLEX ) );
     * and add 
     *     free( F );
     * to end of the routine.
     */ 
    F = (COMPLEX *)time;
    /* bzero( F, ( nT / 2 + 1 ) * sizeof( COMPLEX ) ); */
    memset( F, 0, ( nT / 2 + 1 ) * sizeof( COMPLEX ) );
        
    for( k = 0; k < nF; k++ ) {
        F[freq[k].FBin].Real = freq[k].Real;
        F[freq[k].FBin].Imag = freq[k].Imag;
        }

    rfftb( nT, F, time );
}


/*
 *
 * gfftf: given real sequence compute nF frequency coefficients
 *        of largest magnitude.
 *
 * Parameters:
 *
 * nT   : int              / length of input sequence.
 * time : double[nT]        / input sequence.
 * nF   : int              / number of coefficients wanted.
 * freq : FREQ[nF]         / output coefficients.
 *
 *   FREQ :
 *
 *   struct {
 *       double Real;
 *       double Imag;
 *       double Mag;
 *       int FBin;
 *       }
 */
void gfftf( nT, time, nF, freq )
    int nT;
    int nF;
    FREQ *freq;
    double *time;
{
    COMPLEX *F;

    double *Magnitude;

    int i;
    int k;
    int *Order;

    nT /= 2;

    F = (COMPLEX *)malloc( ( nT + 1 ) * sizeof( COMPLEX ) );

    rfftf( 2 * nT, time, F );

    /************************************************
     *  search for largest magnitude coefficients.. *
     ************************************************/
    Magnitude = (double *)malloc( ( nT + 1 ) * sizeof( double ) );
    Order = (int *)malloc( ( nT + 1 ) * sizeof( int ) );

    for( k = 0; k <= nT; k++ ) {
        Magnitude[k] = F[k].Real * F[k].Real + F[k].Imag * F[k].Imag;
        Order[k] = k;
        }
    sort( nT + 1, Magnitude, Order );

    k = nT;
    for( i = 0; i < nF; i++, k-- ) {
        freq[i].Real = F[Order[k]].Real;
        freq[i].Imag = F[Order[k]].Imag;
        freq[i].Mag  = Magnitude[k];
        freq[i].FBin = Order[k];
        }
 
    free( F );
    free( Order );
    free( Magnitude );
}

/*
 * cfftf2D: 2D forward complex FFT.
 *
 * Parameters:
 *
 * M,N   : int              / array dimensions.
 * T     : COMPLEX[M][N]    / input array.
 * F     : COMPLEX[M][N]    / transform array. may point to same memory as T.
 *
 * F(m,n) = sum(T(k,l)*exp(-i*2*pi*(m*k/M+n*l/N))), k,m=0..M-1; l,n=0..N-1
 */
void cfftf2D( M, N, T, F ) 
    int M,N;
    COMPLEX *T;
    COMPLEX *F;
{
    COMPLEX *W;

    int l;
    int k;
    int n;

    W = (COMPLEX *)malloc( M * sizeof( COMPLEX ) );

    for( k = 0, n = 0; k < M; k++, n += N )
        cfftf( N, &T[n], &F[n] );

    for( l = 0; l < N; l++ ) {
        for( k = 0, n = l; k < M; k++, n += N ) W[k] = F[n];

        cfftf( M, W, W );

        for( k = 0, n = l; k < M; k++, n += N ) F[n] = W[k];
        }

    free( W );
}

/*
 * cfftb2D: 2D inverse complex FFT.
 *
 * Parameters:
 *
 * M,N   : int              / array dimensions.
 * F     : COMPLEX[M][N]    / transform array.
 * T     : COMPLEX[M][N]    / output array. may point to same memory as T.
 *
 * T(m,n) = sum(F(k,l)*exp(i*2*pi*(m*k/M+n*l/N))), k,m=0..M-1; l,n=0..N-1
 */
void cfftb2D( M, N, F, T ) 
    int M,N;
    COMPLEX *F;
    COMPLEX *T;
{
    COMPLEX *W;

    int k;

    if ( T != F ) {
        for( k = 0; k < M * N; k++ ) T[k].Real = F[k].Real;
        }

    for( k = 0; k < M * N; k++ ) T[k].Imag = -F[k].Imag;

    cfftf2D( M, N, T, T );

    for( k = 0; k < M * N; k++ ) T[k].Imag = -T[k].Imag;
}

/*
 * cfftf3D: 3D forward complex FFT.
 *
 * Parameters:
 *
 * L,M,N : int                / array dimensions.
 * T     : COMPLEX[L][M][N]   / input array.
 * F     : COMPLEX[L][M][N]   / transform array. may point to same memory as T.
 *
 * F(l,m,n) = sum(T(i,j,k)*exp(-i*2*pi*(l*i/L+m*j/M+n*k/N))), 
 *                                         l,i=0..L-1; m,j=0..M-1; k,n=0..N-1
 */
void cfftf3D( L, M, N, T, F )
    int L,M,N;
    COMPLEX *T;
    COMPLEX *F;
{
    COMPLEX *W;

    int j;
    int k;
    int n;

    W = (COMPLEX *)malloc( L * sizeof( COMPLEX ) );

    for( k = 0; k < L; k++ ) cfftf2D( M, N, &T[k*M*N], &F[k*M*N] ); 

    for( j = 0; j < M * N; j++ ) {
        for( k = 0, n = j; k < L; k++, n += M * N ) W[k] = F[n];

        cfftf( L, W, W );

        for( k = 0, n = j; k < L; k++, n += M * N ) F[n] = W[k];
        }

    free( W );
}

/*
 * cfftb3D: 3D inverse complex FFT.
 *
 * Parameters:
 *
 * L,M,N : int                / array dimensions.
 * F     : COMPLEX[L][M][N]   / transform array.
 * T     : COMPLEX[L][M][N]   / output array. may point to same memory as T.
 *
 * T(l,m,n) = sum(F(i,j,k)*exp(i*2*pi*(l*i/L+m*j/M+n*k/N))), 
 *                                         l,i=0..L-1; m,j=0..M-1; k,n=0..N-1
 */
void cfftb3D( L, M, N, F, T )
    int L,M,N;
    COMPLEX *F;
    COMPLEX *T;
{
    int k;

    if ( T != F ) {
        for( k = 0; k < L * M * N; k++ ) T[k].Real = F[k].Real;
        }

    for( k = 0; k < L * M * N; k++ ) T[k].Imag = -F[k].Imag;

    cfftf3D( L, M, N, T, T );

    for( k = 0; k < L * M * N; k++ ) T[k].Imag = -T[k].Imag;
}

/*
 * cfftfND: multidimensional (max 32) forward complex FFT.
 *
 * Parameters:
 *
 * N   : int                                   / number of dimensions.
 * D   : int[N]                                / array dimensions.
 * T   : COMPLEX T[D[N-1]][D[N-2]]...[D[0]]    / input array.
 * F   : COMPLEX F[D[N-1]][D[N-2]]...[D[0]]    / transform array.
 *                                             / may point to same memory as T.
 */
void cfftfND( N, D, T, F )
    int  N;
    int *D;
    COMPLEX *T;
    COMPLEX *F;
{
    COMPLEX *W;

    int TotN;
    int WN;

    int i;
    int j;
    int k;
    int l;
    int m;

    int S[32];
    int P[32];
    
    WN = D[0];
    TotN = 1;
    for( i = 0; i < N; i++ ) {
        if ( D[i] > WN ) WN = D[i];
        S[i]  = TotN;
        TotN *= D[i];
        }

    W = (COMPLEX *)malloc( WN * sizeof( COMPLEX ) );

    if ( F != T )
        for( i = 0; i < TotN; i++ )
            F[i] = T[i];

    for( i = 0; i < N; i++ ) {

        for( j = 0; j < N; j++ ) P[j] = 0;

        k = 0;
        for( j = 0; j < TotN / D[i]; j++ ) {

            if ( j != 0 )
                for( l = 0; l < N; l++ )
                    if ( l != i ) {
                        P[l]++;
                        k += S[l];
                        if ( P[l] == D[l] ) {
                            k -= S[l + 1];
                            P[l] = 0;
                        } else {
                            break;
                            }
                        }

            for( l = 0, m = k; l < D[i]; l++, m += S[i] ) W[l] = F[m];

            cfftf( D[i], W, W );

            for( l = 0, m = k; l < D[i]; l++, m += S[i] ) F[m] = W[l];
            }
        }

    free( W );
}


/*
 * cfftbND: multidimensional inverse complex FFT.
 *
 * Parameters:
 *
 * N   : int                                   / number of dimensions.
 * D   : int[N]                                / array dimensions.
 * F   : COMPLEX T[D[N-1]][D[N-2]]...[D[0]]    / transform array.
 * T   : COMPLEX F[D[N-1]][D[N-2]]...[D[0]]    / output array.
 *                                             / may point to same memory as T.
 */
void cfftbND( N, D, F, T )
    int  N;
    int *D;
    COMPLEX *F;
    COMPLEX *T;
{
    int TotN;
    int k;
    
    TotN = D[0];
    for( k = 1; k < N; k++ ) TotN *= D[k];

    if ( T != F ) {
        for( k = 0; k < TotN; k++ ) T[k].Real = F[k].Real;
        }

    for( k = 0; k < TotN; k++ ) T[k].Imag = -F[k].Imag;

    cfftfND( N, D, T, T );

    for( k = 0; k < TotN; k++ ) T[k].Imag = -T[k].Imag;
}
