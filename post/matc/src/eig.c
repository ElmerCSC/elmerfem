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
 *     The eigenvalue/eigenvector extraction routines.
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
/***********************************************************************
|
|  EIG.C - Last Edited 8. 8. 1988
|
***********************************************************************/

/*======================================================================
|Syntax of the manual pages:
|
|FUNCTION NAME(...) params ...
|
$  usage of the function and type of the parameters
?  explane the effects of the function
=  return value and the type of value if not of type int
@  globals effected directly by this routine
!  current known bugs or limitations
&  functions called by this function
~  these functions may interest you as an alternative function or
|  because they control this function somehow
^=====================================================================*/


/*
 * $Id: eig.c,v 1.1.1.1 2005/04/14 13:29:14 vierinen Exp $ 
 *
 * $Log: eig.c,v $
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:33  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"

#define A(i,j) a[n * (i) + (j)]

#define MAXITER 1000
#define EPS 1e-16

VARIABLE *mtr_hesse(var)
     VARIABLE *var;
{
  VARIABLE *res;

  double *a;

  int n;
  
  if (NCOL(var) != NROW(var))
  {
     error( "hesse: matrix must be square, current dimensions: [%d,%d]\n", NROW(var), NCOL(var) );
  }
  
  res = var_temp_copy(var);

  a = MATR(res); n = NROW(res);
  if (NROW(res) == 1) return res;
  
  hesse(a, n, n);

  return res;
}

VARIABLE *mtr_eig(var)
     VARIABLE *var;
{
  VARIABLE *ptr, *res;

  int iter, i, j, k, n;
  double *a, b, s, t;
  
  if (NCOL(var) != NROW(var))
  {
    error("eig: matrix must be square, current dimensions: [%d,%d]\n", NROW(var), NCOL(var));
  }
  
  ptr = var_temp_copy(var);

  a = MATR(ptr); n = NROW(ptr);
  
  if (NROW(ptr) == 1) return ptr;
  
  hesse(a, n, n);
  
  for(iter = 0; iter < MAXITER; iter++)
  {
    
    for (i = 0; i < n - 1; i++)
    {
      s = EPS*(abs(A(i,i))+abs(A(i+1,i+1)));
      if (abs(A(i+1,i)) < s) A(i+1,i) = 0.0;
    }
    
    i = 0;
    do
    {
      for(j = i; j < n - 1; j++) 
	if (A(j+1,j) != 0) break;
      
      for(k = j; k < n - 1; k++)
	if (A(k+1,k) == 0) break;
      
      i = k;
      
    } while(i < n - 1 && k - j + 1 < 3);
    
    if (k - j + 1 < 3) break;
    
    francis(&A(j,j), k - j + 1, n);
  }

  res = var_temp_new(TYPE_DOUBLE, n, 2);

  for(i = 0, j = 0; i < n - 1; i++)
    if (A(i+1,i) == 0)
      M(res,j++,0) = A(i,i);
    else
    {
      b = A(i,i) + A(i+1,i+1); s = b * b;
      t = A(i,i) * A(i+1,i+1) - A(i,i+1)*A(i+1,i); 
      s = s - 4 * t;
      if (s < 0)
      {
        M(res,j, 0) = b / 2;
        M(res,j++, 1) = sqrt(-s) / 2;
        M(res,j, 0) = b / 2;
        M(res,j++,1) = -sqrt(-s) / 2;
      }
      else
      {
        M(res,j++, 0) = b / 2 + sqrt(s) / 2;
        M(res,j++, 0) = b / 2 - sqrt(s) / 2;
      }
      i++;
    }
  
  if (A(n-1, n-2) == 0) M(res,j,0) = A(n-1, n-1);
  
  var_delete_temp(ptr);
  
  return res;
}

void vbcalc(x,v,b,beg,end)
     double x[],v[],*b;
     int beg, end;
{
  double alpha,m,mp1;
  int i;

  m = abs(x[beg]);
  for(i = beg + 1; i <= end; i++)
    m = max(m,abs(x[i]));

  if (m < EPS)
  {
/*
 *   for(i = beg; i <= end; i++) v[i] = 0;
 */
    memset(&v[beg], 0, (end-beg+1)<<3);
  }
  else
  {
    alpha = 0;
    mp1 = 1 / m;
    for(i = beg; i <= end; i++)
    {
      v[i] = x[i] * mp1;
      alpha = alpha + v[i]*v[i];
    }
    alpha = sqrt(alpha);
    *b = 1 / (alpha * (alpha + abs(v[beg])));
    v[beg] = v[beg] + sign(v[beg]) * alpha;
  }
  return;
}

#define H(i,j) h[(i) * N + (j)]

void hesse(h, DIM, N)
     int DIM, N;
     double *h;
{
  double *v, *x, b, s;
  int i, j, k;

  x = (double *)ALLOCMEM(DIM * sizeof(double));
  v = (double *)ALLOCMEM(DIM * sizeof(double));

  for (i = 0; i < DIM - 2; i++)
  {

    for(j = i + 1; j < DIM; j++) x[j] = H(j,i);

    vbcalc(x, v, &b, i + 1, DIM - 1);

    if (v[i+1] == 0) break;

    for(j = i + 2; j < DIM; j++)
    {
      x[j] = v[j]/v[i+1];
      v[j] = b * v[i+1] * v[j];
    }
    v[i+1] = b * v[i+1] * v[i+1];

    for(j = 0; j < DIM; j++)
    {
      s = 0.0;
      for(k = i + 1; k < DIM; k++)
	s = s + H(j,k) * v[k];
      H(j,i+1) = H(j,i+1) - s;
      for(k = i + 2; k < DIM; k++)
	H(j,k) = H(j,k) - s * x[k];
    }

    for(j = 0; j < DIM; j++)
    {
      s = H(i+1,j);
      for(k = i + 2; k < DIM; k++)
	s = s + H(k,j) * x[k];
      for(k = i + 1; k < DIM; k++)
	H(k,j) = H(k,j) - s * v[k];
    }

    for(j = i + 2; j < DIM; j++) H(j,i) = 0;
  }
  FREEMEM((char *)x); FREEMEM((char *)v);

  return;
}


void francis(h, DIM, N)
     int DIM, N;
     double *h;
{
  double x[3], v[3], b, s, t, bv, v0i;
  int i, i1, j, k, n, m, end;

  n = DIM - 1; m = n - 1;

  t = H(m,m) * H(n,n) - H(m,n) * H(n,m);
  s = H(m,m) + H(n,n);

  x[0] = H(0,0)*H(0,0)+H(0,1)*H(1,0)-s*H(0,0)+t;
  x[1] = H(1,0) * (H(0,0) + H(1,1) - s);
  x[2] = H(1,0) * H(2,1);

  vbcalc(x, v, &b, 0, 2);

  if (v[0] == 0) return;



/*
 * for(i = 1; i < 3; i++)
 * {
 *   x[i] = v[i]/v[0];
 *   v[i] = b * v[0] * v[i];
 * }
 */
  bv  = b * v[0];

  x[1] = v[1] / v[0];
  v[1] = bv * v[1];

  x[2] = v[2] / v[0];
  v[2] = bv * v[2];



  x[0] = 1; v[0] = b * v[0] * v[0];

  for(i = 0; i < DIM; i++)
  {
/*
 *   s = 0.0;
 *   for(j = 0; j < 3; j++)
 *     s = s + H(i,j) * v[j];
 */
    i1 = i * N;
    s = h[i1] * v[0] + h[i1+1]*v[1] + h[i1+2]*v[2];


    H(i,0) = H(i,0) - s;
/*
 *   for(j = 1; j < 3; j++)
 *     H(i,j) = H(i,j) - s * x[j];
 */
    h[i1+1] = h[i1+1] - s * x[1];
    h[i1+2] = h[i1+2] - s * x[2];
  }

  for(i = 0; i < DIM; i++)
  {
/*
 *   s = H(0,i);
 *   for(j = 1; j < 3; j++)
 *     s = s + H(j,i) * x[j];
 */
    s = h[i] + h[N+i]*x[1] + h[(N<<1)+i]*x[2]; 
/*
 *   for(j = 0; j < 3; j++)
 *     H(j,i) = H(j,i) - s * v[j];
 */
    h[i] = h[i] - s * v[0];
    h[N+i] = h[N+i] - s * v[1];
    h[(N<<1)+i] = h[(N<<1)+i] - s * v[2];
  }

  for (i = 0; i < DIM - 2; i++)
  {

    end = min(2, DIM - i - 2);
    for(j = 0; j <= end; j++) x[j] = H(i+j+1,i);

    vbcalc(x, v, &b, 0, end);

    if (v[0] == 0) break;

    for(j = 1; j <= end; j++)
    {
      x[j] = v[j]/v[0];
      v[j] = b * v[0] * v[j];
    }
    x[0] = 1; v[0] = b * v[0] * v[0];

    for(j = 0; j < DIM; j++)
    {
      s = 0.0;
      for(k = 0; k <= end; k++)
	s = s + H(j,i+k+1) * v[k];
      H(j,i+1) = H(j,i+1) - s;
      for(k = 1; k <= end; k++)
	H(j,i+k+1) = H(j,i+k+1) - s * x[k];
    }

    for(j = 0; j < DIM; j++)
    {
      s = H(i+1,j);
      for(k = 1; k <= end; k++)
	s = s + H(i+k+1,j) * x[k];
      for(k = 0; k <= end; k++)
	H(i+k+1,j) = H(i+k+1,j) - s * v[k];
    }

    for(j = i + 2; j < DIM; j++) H(j,i) = 0;
  }
  return;
}
