/*****************************************************************************
 *
 *  Elmer, A Finite Element Software for Multiphysical Problems
 *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
 * 
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program (in file fem/GPL-2); if not, write to the 
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 *  Boston, MA 02110-1301, USA.
 *
 *****************************************************************************/

/*******************************************************************************
 *
 *     Solve symmetric positive definite eigenvalue problem.
 *
 *******************************************************************************
 *
 *                     Author:       Juha Ruokolainen
 *
 *                    Address: CSC - IT Center for Science Ltd.
 *                                Keilaranta 14, P.O. BOX 405
 *                                  02101 Espoo, Finland
 *                                  Tel. +358 0 457 2723
 *                                Telefax: +358 0 457 2302
 *                              EMail: Juha.Ruokolainen@csc.fi
 *
 *                       Date: 30 May 1996
 *
 *                Modified by:
 *
 *       Date of modification:
 *
 ******************************************************************************/

/*
 * $Id: jacobi.c,v 1.1.1.1 2005/04/14 13:29:14 vierinen Exp $ 
 *
 * $Log: jacobi.c,v $
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:43  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"

VARIABLE *mtr_jacob(a) VARIABLE *a;
{
  VARIABLE *x, *ev;

  double *b, *d, rtol;
  int i, n;

  if (NROW(a) != NCOL(a))
    error("Jacob: Matrix must be square.\n");

  b = MATR(NEXT(a));
  n = NROW(a);

  if (NROW(NEXT(a)) != NCOL(NEXT(a)) || n != NROW(NEXT(a))) 
    error("Jacob: Matrix dimensions incompatible.\n");

  rtol = *MATR(NEXT(NEXT(a)));

  x  = var_new("eigv", TYPE_DOUBLE, NROW(a), NCOL(a));
  d  = (double *)ALLOCMEM(n * sizeof(double));
  ev = var_temp_new(TYPE_DOUBLE, 1, n);

  jacobi(MATR(a), b, MATR(x), MATR(ev), d, n, rtol);
  FREEMEM((char *)d);

  return ev;
};

/************************************************************************
  P R O G R A M 
  To solve the generalized eigenproblem using the
  generalized Jacobi iteration

  INPUT

  a(n,n)    = Stiffness matrix (assumed positive definite)
  b(n,n)    = Mass matrix (assumed positive definite)
  x(n,n)    = Matrix storing eigenvectors on solution exit
  eigv(n)   = Vector storing eigenvalues on solution exit
  d(n)      = Working vector
  n         = Order of the matrices a and b
  rtol      = Convergence tolerance (usually set to 10.**-12)
  nsmax     = Maximum number of sweeps allowed (usually set to 15)

  OUTPUT

  a(n,n)    = Diagonalized stiffness matrix
  b(n,n)    = Diagonalized mass matrix
  x(n,n)    = Eigenvectors stored columnwise
  eigv(n)   = Eigenvalues

*************************************************************************/

int jacobi(a,b,x,eigv,d,n,rtol)
  double a[],b[],x[],eigv[],d[],rtol;
  int n;
{
  register int i,j,k,ii,jj;
  int    nsmax=50,        /* Max number of sweeps */
         nsweep,          /* Current sweep number */
         nr,           
         jp1,jm1,kp1,km1,
         convergence;

  double eps,
         eptola,eptolb,
         akk,ajj,ab,
         check,sqch,
         d1,d2,den,
         ca,cg,
         aj,bj,ak,bk,
         xj,xk,
         tol,dif,
         epsa,epsb,
         bb;              /* Scale mass matrix */

/************************************************************************
  Initialize eigenvalue and eigenvector matrices 
*************************************************************************/
  for( i=0 ; i<n ; i++)
  {
    ii = i*n+i;            /* address of diagonal element */
    if (a[ii]<=0.0 || b[ii]<=0.0) return 0;
    eigv[i] = d[i] = a[ii]/b[ii];
    x[ii] = 1.0;           /* make an unit matrix */
  }
  if(n==1) return 1;      /* Return if single degree of freedom system */

/************************************************************************
  Initialize sweep counter and begin iteration
*************************************************************************/
  nr=n - 1;
  for ( nsweep = 0; nsweep < nsmax; nsweep++)
  {
/************************************************************************
    Check if present off-diagonal element is large enough to require zeroing
*************************************************************************/

    eps = pow(0.01, 2.0 * (nsweep + 1));
    for(j=0 ; j<nr ; j++)
    {
      jj=j+1;
      for(k=jj ; k<n ; k++)
      {
        eptola=( a[j*n+k]*a[j*n+k] / ( a[j*n+j]*a[k*n+k] ) );
        eptolb=( b[j*n+k]*b[j*n+k] / ( b[j*n+j]*b[k*n+k] ) );
        if ( eptola >= eps || eptolb >=eps )
        {

/************************************************************************
          if zeroing is required, calculate the rotation matrix elements
*************************************************************************/

          akk=a[k*n+k]*b[j*n+k] - b[k*n+k]*a[j*n+k];
          ajj=a[j*n+j]*b[j*n+k] - b[j*n+j]*a[j*n+k];
          ab =a[j*n+j]*b[k*n+k] - a[k*n+k]*b[j*n+j];
          check=(ab*ab+4.0*akk*ajj)/4.0;

          if (check<=0.0)
          {
            printf("***Error   solution stop in *jacobi*\n");
            printf("        check = %20.14e\n", check);
            return 1;
          }
          sqch= sqrt(check);
          d1  = ab/2.0+sqch;
          d2  = ab/2.0-sqch;
          den = d1;
          if ( abs(d2) > abs(d1) ) den=d2;

          if (den == 0.0)
          {
            ca=0.0;
            cg= -a[j*n+k] / a[k*n+k];
          }
          else
          {
            ca= akk/den;
            cg= -ajj/den;
          }
/************************************************************************
  Perform the generalized rotation to zero the present off-diagonal element
*************************************************************************/
          if (n != 2)
          {
            jp1=j+1;
            jm1=j-1;
            kp1=k+1;
            km1=k-1;
/**************************************/
            if ( jm1 >= 0)
              for (i=0 ; i<=jm1 ; i++)
              {
                aj = a[i*n+j];
                bj = b[i*n+j];
                ak = a[i*n+k];
                bk = b[i*n+k];
                a[i*n+j] = aj+cg*ak;
                b[i*n+j] = bj+cg*bk;
                a[i*n+k] = ak+ca*aj;
                b[i*n+k] = bk+ca*bj;
              };
/**************************************/
            if ((kp1-n+1) <= 0)
              for (i=kp1 ; i<n ; i++)
              {
                aj = a[j*n+i];
                bj = b[j*n+i];
                ak = a[k*n+i];
                bk = b[k*n+i];
                a[j*n+i] = aj+cg*ak;
                b[j*n+i] = bj+cg*bk;
                a[k*n+i] = ak+ca*aj;
                b[k*n+i] = bk+ca*bj;
              };
/**************************************/
            if ((jp1-km1) <= 0.0)
              for (i=jp1 ; i<=km1 ; i++)
              {
                aj = a[j*n+i];
                bj = b[j*n+i];
                ak = a[i*n+k];
                bk = b[i*n+k];
                a[j*n+i] = aj+cg*ak;
                b[j*n+i] = bj+cg*bk;
                a[i*n+k] = ak+ca*aj;
                b[i*n+k] = bk+ca*bj;
              };
          };
   
          ak = a[k*n+k];
          bk = b[k*n+k];
          a[k*n+k] = ak+2.0*ca*a[j*n+k]+ca*ca*a[j*n+j];
          b[k*n+k] = bk+2.0*ca*b[j*n+k]+ca*ca*b[j*n+j];
          a[j*n+j] = a[j*n+j]+2.0*cg*a[j*n+k]+cg*cg*ak;
          b[j*n+j] = b[j*n+j]+2.0*cg*b[j*n+k]+cg*cg*bk;
          a[j*n+k] = 0.0;
          b[j*n+k] = 0.0;

/************************************************************************
          Update the eigenvector matrix after each rotation
*************************************************************************/
          for (i=0 ; i<n ; i++)
          {
            xj = x[i*n+j];
            xk = x[i*n+k];
            x[i*n+j] = xj+cg*xk;
            x[i*n+k] = xk+ca*xj;
          };
        };
      };
    };
/************************************************************************
    Update the eigenvalues after each sweep
*************************************************************************/
    for (i=0 ; i<n ; i++)
    {
      ii=i*n+i;
      if (a[ii]<=0.0 || b[ii]<=0.0)
      {
         error( "*** Error  solution stop in *jacobi*\n Matrix not positive definite.");
      };
      eigv[i] = a[ii] / b[ii];
    };

/************************************************************************
    check for convergence
*************************************************************************/
    convergence = 1;
    for(i=0 ; i<n ; i++)
    {
      tol = rtol*d[i];
      dif = abs(eigv[i]-d[i]);
      if (dif > tol) convergence = 0;
      if (! convergence) break;
    };
/************************************************************************
    check if all off-diag elements are satisfactorily small
*************************************************************************/
    if (convergence)
    {
      eps=rtol*rtol;
      for (j=0 ; j<nr ; j++)
      {
        jj=j+1;
        for (k=jj ; k<n ; k++)
        {
          epsa = ( a[j*n+k]*a[j*n+k] / ( a[j*n+j]*a[k*n+k] ) );
          epsb = ( b[j*n+k]*b[j*n+k] / ( b[j*n+j]*b[k*n+k] ) );
          if ( epsa >= eps || epsb >=eps ) convergence = 0;
          if (! convergence) break;
        };
        if (! convergence) break;
      };
    };
    if (! convergence)
    {
      for (i=0 ; i<n ; i++)
      d[i] = eigv[i];
    };
    if (convergence) break;
  };
/************************************************************************
  Fill out bottom triangle of resultant matrices and scale eigenvectors
*************************************************************************/
  for (i=0 ; i<n ; i++)
    for (j=i ; j<n ; j++)
    {
      b[j*n+i] = b[i*n+j];
      a[j*n+i] = a[i*n+j];
    };

  for (j=0 ; j<n ; j++)
  {
    bb = sqrt( b[j*n+j] );
    for (k=0 ; k<n ; k++)
      x[k*n+j] /= bb;
  };
  PrintOut( "jacobi: nsweeps %d\n", nsweep );
  return 1 ;
}
