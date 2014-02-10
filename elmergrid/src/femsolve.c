/*  
   ElmerGrid - A simple mesh generation and manipulation utility  
   Copyright (C) 1995- , CSC - IT Center for Science Ltd.   

   Author: Peter Råback
   Email: Peter.Raback@csc.fi
   Address: CSC - IT Center for Science Ltd.
            Keilaranta 14
            02101 Espoo, Finland

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/* -------------------------------:  femsolve.c  :----------------------------

   This module is for solving a dense system of linear equations. For sparse systems
   a special sparse matrix library should be used.
*/


#include <math.h>
#include <stdio.h>
#include "common.h"
#include "nrutil.h"
#include "femdef.h"
#include "femsolve.h"



/*******************************************************************************/

void Symmetrize(Real **vf,int sides)
/* Symmetrize the matrix by replacing two elements that should be equal by their
   mean value. If the matrix is not symmetric, Symmetrize should be called before 
   the subroutine Normalize.
*/
{
  int i,j;
  Real temp;

  for (i=1;i<=sides;i++) 
    for (j=i+1;j<=sides;j++) {
      temp = .5*(vf[i][j] + vf[j][i]);
      vf[j][i] = vf[i][j] = temp;
    }
}


void Normalize(Real **vf, const Real *b,int sides)
/* Normalize a symmetric matrix F (vf) by solving a diagonal 
   matrix D (diag) so that the row sums in matrix DFD coincides 
   with those given in vector b. The Newton method 
   is used to solve the diagonal matrix D. 
*/
{
  /* Criteria for ending the iteration */
  int itmax = 10;
  Real eps = 1e-20;
  Real **jac,*rest,*diag,sum;
  int i,j,k,it;
  int *indx;

  printf("Normalization of matrices is currently not supported\n");
  printf("If you want to reactivate it, rewrite the LU decomposition\n");
  return;

  jac = Rmatrix(1,sides,1,sides);
  rest = Rvector(1,sides);
  diag = Rvector(1,sides);
  indx = Ivector(1,sides);

  for (i=1;i<=sides;i++) diag[i] = 1.;


  for(it=1;it<=itmax;it++) {

    /* Calculate the residual sum(DAD)-b. */
    for (i=1; i<=sides; i++) {
      sum = 0.;
      for (j=1; j<=sides; j++) 
        sum += diag[j] * vf[i][j];
      sum *= diag[i];
      rest[i] = sum - b[i];
    }
    sum = 0.;
    for (i=1; i<=sides; i++) 
      sum += rest[i]*rest[i]/b[i];
    sum /= sides;
    
    if(sum < eps) break;

    /* Make the Jacobian matrix. */
    for (i=1; i<=sides; i++) {
      jac[i][i] = 0.;
      for(j=1; j<=sides; j++)
	jac[i][j] = diag[i]*vf[i][j];
      for(k=1; k<=sides;k++)
        jac[i][i] += diag[k]*vf[i][k];
    }

    /* Solve the equation jac * x = rest using Numerical Recipes routines. */     
    /* These have been removed due to copyright violiation */

    /* New approximation. */
    for (i=1; i<=sides; i++)
      diag[i] -= rest[i];    
  } 

  /* Normalize the viewfactors. */
  for (i=1; i<=sides; i++)
    for (j=1; j<=sides; j++)
      vf[i][j] *= diag[i] * diag[j];

  free_Rmatrix(jac,1,sides,1,sides);
  free_Rvector(rest,1,sides);
  free_Rvector(diag,1,sides);
  free_Ivector(indx,1,sides);
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
void SortIndex( int N, double *Key, int *Ord )
{
    double CurrentKey;

    int CurrentOrd;

    int CurLastPos;
    int CurHalfPos;

    int i;
    int j; 
 
    /* Initialize order */
    for(i=1;i<=N;i++)
      Ord[i] = i;

    CurHalfPos = N / 2 + 1;
    CurLastPos = N;
    while( 1 ) {
        if ( CurHalfPos > 1 ) {
            CurHalfPos--;
            CurrentKey = Key[CurHalfPos];
            CurrentOrd = Ord[CurHalfPos];
        } else {
            CurrentKey = Key[CurLastPos];
            CurrentOrd = Ord[CurLastPos];
            Key[CurLastPos] = Key[1];
            Ord[CurLastPos] = Ord[1];
            CurLastPos--;
            if ( CurLastPos == 1 ) {
	      Key[1] = CurrentKey;
	      Ord[1] = CurrentOrd;
	      return;
	      }
            }
        i = CurHalfPos;
        j = 2 * CurHalfPos;
        while( j <= CurLastPos ) {
            if ( j < CurLastPos && Key[j] < Key[j + 1] ) {
                j++;
                }
            if ( CurrentKey < Key[j] ) {
                Key[i] = Key[j];
                Ord[i] = Ord[j];
                i  = j;
                j += i;
            } else {
                j = CurLastPos + 1;
                }
            }
        Key[i] = CurrentKey;
        Ord[i] = CurrentOrd;
        }

}


