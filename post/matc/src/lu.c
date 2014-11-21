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
 *     LU decomposition.
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
|  LU.C - Last Edited 8. 8. 1988
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
 * $Id: lu.c,v 1.1.1.1 2005/04/14 13:29:14 vierinen Exp $ 
 *
 * $Log: lu.c,v $
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:45  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"

#define A(i,j) a[n * (i) + (j)]

VARIABLE *mtr_LUD(var)
     VARIABLE *var;
{
  VARIABLE *res;
  int i, n, *pivot;
  double *a;
  
  if (NCOL(var) != NROW(var))
  {        
    error("LUD: Matrix must be square.\n");
  }
  
  res = var_temp_copy(var); 
  a = MATR(res); n = NROW(res);

  pivot = (int *)ALLOCMEM(n * sizeof(int));
  
  LUDecomp(a, n, pivot);
  
  FREEMEM((char *)pivot);
  
  return res;
}

VARIABLE *mtr_det(var)
     VARIABLE *var;
{
  VARIABLE *tmp, *res;
  int i, n, *pivot;
  double det, *a;
  
  if (NCOL(var) != NROW(var))
  {        
    error("Det: Matrix must be square.\n");
  }
  
  tmp = var_temp_copy(var); 

  a = MATR(tmp); n = NROW(tmp);
  
  pivot = (int *)ALLOCMEM(n * sizeof(int));
  
  LUDecomp(a, n, pivot);
  
  det = 1.0;
  for(i = 0; i < n; i++)
  {
    det *= A(i,i);
    if (pivot[i] != i) det = -det;
  }
  
  FREEMEM((char *)pivot); var_delete_temp(tmp);
  
  res = var_temp_new(TYPE_DOUBLE,1,1);

  M(res,0,0) = det;
  
  return res;
}

VARIABLE *mtr_inv(var)
     VARIABLE *var;
{
  VARIABLE *ptr;

  int i, j , k, n, *pivot;
  double s, *a;
  
  if (NCOL(var) != NROW(var))
  {        
    error("Inv: Matrix must be square.\n");
  }

  ptr = var_temp_copy(var); 

  a = MATR(ptr); n = NROW(ptr);
  
  pivot = (int *)ALLOCMEM(n * sizeof(int));  

  /*
   *  AP = LU
   */
  LUDecomp(a, n, pivot);
  for(i = 0; i < n; i++)
  {
    if (A(i,i) == 0) 
      error("Inv: Matrix is singular.\n");
    A(i,i) = 1.0 / A(i,i);
  }

  /*  
   *  INV(U)
   */
  for(i = n - 2; i >= 0; i--)
    for(j = n - 1; j > i; j--)
    {
      s = 0.0;
      for(k = i+1; k <= j; k++)
	if (k != j)  
	  s = s - A(i,k) * A(k,j);
	else 
	  s = s - A(i,k);
      A(i,j) = s;
    }

  /*
   * INV(L)
   */
  for(i = n - 2; i >= 0; i--)
    for(j = n - 1; j > i; j--)
    {
      s = 0.0;
      for(k = i + 1; k <= j; k++) 
	s = s - A(j,k) * A(k,i);
      A(j,i) = A(i,i) * s;
    }
  
  /* 
   * A  = INV(AP)
   */
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
    {
      s = 0.0;
      for(k = max(i,j); k < n; k++)
	if (k != i)
	  s += A(i,k) * A(k,j);
	else 
	  s += A(k,j);
      A(i,j) = s;
    }

  /*
   * A = INV(A) (at last)
   */
  for(i = 0; i < n; i++)
    if (pivot[i] != i)
      for(j = 0; j < n; j++)
      {
	s = A(i,j);
	A(i,j) = A(pivot[i],j);
	A(pivot[i],j) = s;
      }
  
  FREEMEM((char *)pivot);

  return ptr;
}

/*
 * LU- decomposition by gaussian elimination. Row pivoting is used.
 * 
 * result : AP = L'U ; L' = LD; pivot[i] is the swapped column 
 * number for column i (for pivoting).
 *
 * Result is stored in place of original matrix.
 *
 */
void LUDecomp(a, n, pivot)
   double *a;
   int n, pivot[];
{
  double swap;
  int i, j, k, l;
  
  for (i = 0; i < n - 1; i++)
  {
    j = i;
    for(k = i + 1; k < n; k++)
      if (abs(A(i,k)) > abs(A(j,k))) j = k;
    
    if (A(i,j) == 0.0) 
    {
      error("LUDecomp: Matrix is singular.\n");
    }
    
    pivot[i] = j;
    
    if (j != i)
    {
      swap = A(i,i);
      A(i,i) = A(i,j);
      A(i,j) = swap;
    }
    
    for(k = i + 1; k < n; k++)
      A(i,k) = A(i,k) / A(i,i);
    
    for(k = i + 1; k < n; k++) 
    {
      if (j != i)
      {
	swap = A(k,i);
	A(k,i) = A(k,j); 
	A(k,j) = swap;
      }
      for(l = i + 1; l < n; l++)
	A(k,l) = A(k,l) -  A(i,l) * A(k,i);
    }
  }
  
  pivot[n - 1] = n - 1;
  if ( A(n-1,n-1) == 0.0)
  {
      error("LUDecomp: Matrix is singular.\n");
  }
}
