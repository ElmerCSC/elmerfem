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
 *     MATC operator routines.
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
|   Oper.C - Last Edited 15. 8. 1988
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
~  these functions may intrest you as an alternative function or
|  because they control this function somehow
=====================================================================*/


/*
 * $Id: oper.c,v 1.1.1.1 2005/04/14 13:29:14 vierinen Exp $ 
 *
 * $Log: oper.c,v $
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:51  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"

#ifdef TYPE
#undef TYPE
#endif
#ifdef NROW
#undef NROW
#endif
#ifdef NCOL
#undef NCOL
#endif
#ifdef MATR
#undef MATR
#endif
#ifdef M
#undef M
#endif
#ifdef MATSIZE
#undef MATSIZE
#endif

#define TYPE(mat) (mat)->type
#define NROW(mat) (mat)->nrow
#define NCOL(mat) (mat)->ncol
#define MATR(mat) (mat)->data
#define MATSIZE(mat) (NROW(mat) * NCOL(mat) * sizeof(double))

#define MA(i,j) a[ncola * (i) + (j)]
#define MB(i,j) b[ncolb * (i) + (j)]
#define MC(i,j) c[ncolc * (i) + (j)]

MATRIX *mat_new(type, nrow, ncol) int type, nrow, ncol;
{
   MATRIX *res;

    res = (MATRIX *)ALLOCMEM(MATRIXSIZE);
    TYPE(res) = type;
    NROW(res) = nrow;
    NCOL(res) = ncol;
    MATR(res) = (double *)ALLOCMEM(MATSIZE(res));

    return res;
}

MATRIX *mat_copy(mat) MATRIX *mat;
{
    MATRIX *res;

    if (mat == (MATRIX *)NULL) return NULL;
     
    res = mat_new(TYPE(mat), NROW(mat), NCOL(mat));
    memcpy((char *)MATR(res), (char *)MATR(mat), MATSIZE(mat));

    return res;
}

void mat_free(mat) MATRIX *mat;
{
    if (mat == (MATRIX *)NULL) return;

    FREEMEM((char *)MATR(mat));
    FREEMEM((char *)mat);
}

MATRIX *opr_vector(A, B)
       MATRIX *A, *B;
{
  MATRIX *C;

  int i, n, inc;

  double *a = MATR(A), *b = MATR(B), *c;

  inc = ((int)*b > (int)*a) ? 1:(-1);
  n = abs((int)*b - (int)*a) + 1;

  C = mat_new(TYPE_DOUBLE, 1, n);
  c = MATR(C);

  for(i = 0; i < n; i++) 
    *c++ = (int)*a + i*inc;

  return C;
} 

MATRIX *opr_resize(A, B)
       MATRIX *A, *B;
{
  MATRIX *C;

  double *a = MATR(A), *b = MATR(B), *c;
  int i, j, n, m;
  
  if (NCOL(B) >= 2)
  {
    i = *b++; j = *b++;
  } 
  else
  {
    i = 1; j = *b;
  }
  
  if (i < 1 || j < 1)
    error("resize: invalid size for and array");
  
  C = mat_new(TYPE(A), i, j);
  c = MATR(C);

  n = i * j; m = NROW(A) * NCOL(A);

  for(i = j = 0; i < n; i++)
  {
    *c++ = a[j++];
    if (j == m) j = 0;
  }

  return C;
}

MATRIX *opr_apply(A)
     MATRIX *A;
{
  VARIABLE *store, *ptr;
  MATRIX *C = NULL;
  
  store = var_temp_new(TYPE_STRING, NROW(A), NCOL(A));
  REFCNT(store) = 0;
  mat_free(store->this);
  store->this = A;
  REFCNT(store)++;
  
  ptr = (VARIABLE *)com_apply(store);

  var_delete_temp(store);

  if ( ptr ) C = mat_copy(ptr->this);

  return C;
}

MATRIX *opr_add(A, B)
     MATRIX *A, *B;
{
  MATRIX *C;

  int i;

  int nrowa = NROW(A), ncola = NCOL(A);
  int nrowb = NROW(B), ncolb = NCOL(B);

  double *a = MATR(A), *b = MATR(B), *c;

  double value;
  
  if (nrowa == nrowb && ncola == ncolb)
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C); 
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++)  *c++ = *a++ + *b++;
  }

  else if (nrowa == 1 && ncola == 1)
  {
    C = mat_new(TYPE(B), nrowb, ncolb); c = MATR(C);
    value = *a; nrowb *= ncolb;
    for(i = 0; i < nrowb; i++) *c++ = value + *b++;
  }
  
  else if (nrowb == 1 && ncolb == 1)
  {
    C = mat_new(TYPE(A), nrowa, ncola); c = MATR(C);
    value = *b; nrowa *= ncola;
    for(i = 0; i < nrowa; i++) *c++ = value + *a++;
  }

  else
    error("Add: Incompatible for addition.\n");
  
  return C;
}

MATRIX *opr_minus(A)
     MATRIX *A;
{
  MATRIX *C;
  int i;
  
  int nrowa = NROW(A), ncola = NCOL(A);

  double *a = MATR(A), *c;

  C = mat_new(TYPE(A), nrowa, ncola); c = MATR(C);

  nrowa *= ncola;
  for(i = 0; i < nrowa; i++) *c++ = -*a++;
  
  return C;
}

MATRIX *opr_subs(A, B)
     MATRIX *A, *B;
{
  MATRIX *C;

  double value;

  int i;
  
  int nrowa = NROW(A), ncola = NCOL(A);
  int nrowb = NROW(B), ncolb = NCOL(B);

  double *a = MATR(A), *b = MATR(B), *c;
  
  if (nrowa == nrowb && ncola == ncolb)
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++) *c++ = *a++ - *b++;
  }
  else if (nrowa == 1 && ncola == 1)
  {
    C = mat_new(TYPE(B),nrowb, ncolb); c = MATR(C);
    value = *a; nrowb *= ncolb;
    for(i = 0; i < nrowb; i++) *c++ = value - *b++;
  }
  else if (nrowb == 1 && ncolb == 1)
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    value = *b; nrowa *= ncola;
    for(i = 0; i < nrowa; i++) *c++ = *a++ - value;
  }
  else
    error("Substr: Incompatible for addition.\n");
  
  return C;
}

MATRIX *opr_mul(A, B)
     MATRIX *A, *B;
{
  MATRIX *C;

  double value,s;
  int i, j, k;

  int nrowa = NROW(A), ncola = NCOL(A);
  int nrowb = NROW(B), ncolb = NCOL(B);

  double *a = MATR(A), *b = MATR(B), *c;

  if (nrowa == 1 && ncola == 1)
  {
    C = mat_new(TYPE(B), nrowb, ncolb); c = MATR(C);
    value = *a; nrowb *= ncolb;
    for(i = 0; i < nrowb; i++) *c++ = value * *b++;
  }
  else if (nrowb == 1 && ncolb == 1)
  {
    C = mat_new(TYPE(A), nrowa, ncola); c = MATR(C);
    value = *b; nrowa *= ncola;
    for(i = 0; i < nrowa; i++) *c++ = value * *a++;
  }
  else if (ncola == nrowb)
  {
    C = mat_new(TYPE(A), nrowa,ncolb); 
    c = MATR(C);
    for(i = 0; i < nrowa; i++)
    {
      for(j = 0; j < ncolb; j++)
      {
        s = 0;
	for(k = 0; k < ncola; k++) s += a[k] * MB(k,j);
	*c++ = s;
      }
      a += ncola;
    }
  }
  else if ( ncola == ncolb && nrowa == nrowb )
  {
    C = mat_new(TYPE(A), nrowa,ncolb); 
    c = MATR(C);
    k = 0;
    for( i = 0; i < nrowa; i++ )
        for( j = 0; j < ncolb; j++,k++ ) c[k] = a[k] * b[k];
  }
  else 
  {
     error("Mul: Incompatible for multiplication.\n");
  }
 
  return C;
}

MATRIX *opr_pmul(A, B)
     MATRIX *A, *B;
{
  MATRIX *C;

  double value;
  int i;

  int nrowa = NROW(A), ncola = NCOL(A);
  int nrowb = NROW(B), ncolb = NCOL(B);
  
  double *a = MATR(A), *b = MATR(B), *c;

  if (nrowa == nrowb && ncola == ncolb)
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++) *c++ = *a++ * *b++;
  }
  else if (nrowa == 1 && ncola == 1)
  {
    C = mat_new(TYPE(B), nrowb, ncolb); c = MATR(C);
    value = *a; nrowb *= ncolb;
    for(i = 0; i < nrowb; i++) *c++ = value * *b++;
  }
  else if (nrowb == 1 && ncolb == 1)
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    value = *b; nrowa *= ncola;
    for(i = 0; i < nrowa; i++) *c++ = *a++ * value;
  }
  else 
    error("PMul: Incompatible for pointwise multiplication.\n");
  
  return C;
}

MATRIX *opr_div(A, B)
     MATRIX *A, *B;
{
  MATRIX *C;

  double value;
  int i;

  int nrowa = NROW(A), ncola = NCOL(A);
  int nrowb = NROW(B), ncolb = NCOL(B);
  
  double *a = MATR(A), *b = MATR(B), *c;

  if (nrowa == nrowb && ncola == ncolb)
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++) *c++ = *a++ / *b++;
  }
  else if (nrowa == 1 && ncola == 1)
  {
    C = mat_new(TYPE(B), nrowb, ncolb); c = MATR(C);
    value = *a; nrowb *= ncolb;
    for(i = 0; i < nrowb; i++) *c++ = value / *b++;
  }
  else if (nrowb == 1 && ncolb == 1)
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    value = *b; nrowa *= ncola;
    for(i = 0; i < nrowa; i++) *c++ = *a++ / value;
  }
  else 
    error("Div: Incompatible for division.\n");
  
  return C;
}

MATRIX *opr_pow(A,B)
     MATRIX *A, *B;
{
  MATRIX *C;

  int i, j, k, l, power;

  int nrowa = NROW(A), ncola = NCOL(A);
  int nrowb = NROW(B), ncolb = NCOL(B);
  int ncolc;

  double *a = MATR(A), *b = MATR(B), *c;

  double *v, value;

  if (nrowb != 1 || ncolb != 1) error("Pow: Matrix ^ Matrix ?.\n");

  if (nrowa == 1 || ncola != nrowa)
  {
    C = mat_new(TYPE(A), nrowa, ncola); c = MATR(C);
    value = *b; nrowa *= ncola;
    for(i = 0; i < nrowa; i++) *c++ = pow(*a++, value);

    return C;
  }

  power = (int)*b;

  if (power == 0) {
    C = mat_new(TYPE(A), nrowa, ncola);
    ncolc = ncola; c = MATR(C);
    for(i = 0; i < nrowa; i++) MC(i,i) = 1;
  }
  else if (abs(power) == 1)
  {
    C = mat_copy(A);
  }
  else
  {
    v = (double *)ALLOCMEM(nrowa * sizeof(double));

    C = mat_new(TYPE(A), nrowa, nrowa);
    c = MATR(C); b = MATR(A);

    for(l = 1; l < abs(power); l++) {
      for(i = 0; i < nrowa; i++)
      {
	for(j = 0; j < nrowa; j++)
	{
	  v[j] = 0.0;
	  for(k = 0; k < nrowa; k++) v[j] +=  b[k] * MA(k,j);
	}
	for(j = 0; j < nrowa; j++) *c++ = v[j];
	b += nrowa;
      }
      a = MATR(A); b = c = MATR(C);
    }
    FREEMEM((char *)v);
  }

  if (power < 0)
  {
    VARIABLE *inv, *tmp;
    tmp = (VARIABLE *)ALLOCMEM(VARIABLESIZE);
    tmp -> this = C;
    inv = mtr_inv(tmp);
    mat_free(C);
    FREEMEM((char *)tmp);
    C = inv->this;
    REFCNT(inv)++;
    var_delete_temp(inv);
  }

  return C;
}

MATRIX *opr_trans(A)
     MATRIX *A;
{
  MATRIX *C;

  int i, j;

  int ncola = NCOL(A), nrowa = NROW(A);
  int ncolc;  

  double *a = MATR(A), *c;

  C = mat_new(TYPE(A),ncola,nrowa);
  c = MATR(C); ncolc = nrowa;

  for(i = 0; i < nrowa; i++)
    for(j = 0; j < ncola; j++) MC(j,i) = *a++;
  
  return C;
}

MATRIX *opr_reduction(A, B)
     MATRIX *A, *B;
{
  MATRIX *C;

  int i;

  int nrowa = NROW(A), ncola = NCOL(A);
  int nrowb = NROW(B), ncolb = NCOL(B);

  double *a = MATR(A), *b = MATR(B), *c;

  if ((nrowa == nrowb) && (ncola == ncolb))
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++, a++) *c++ = (*b++) ? (*a) : (0);
  }
  else
    error("Incompatible for reduction.\n");

  return C;
}


MATRIX *opr_lt(A, B)
     MATRIX *A, *B;
{
  MATRIX *C;

  int i;
  
  int nrowa = NROW(A), ncola = NCOL(A);
  int nrowb = NROW(B), ncolb = NCOL(B);

  double *a = MATR(A), *b = MATR(B), *c;

  if ((nrowa == 1) && (ncola == 1))
  {
    C = mat_new(TYPE(B),nrowb, ncolb); c = MATR(C);
    nrowb *= ncolb;
    for(i = 0; i < nrowb; i++,c++) if (*a < b[i]) *c = 1;
  }
  
  else if ((nrowb == 1) && (ncolb == 1))
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++,c++) if (a[i] < *b) *c = 1;
  }
  
  else if ((nrowa == nrowb) && (ncola == ncolb))
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++,c++) if (a[i] < b[i]) *c = 1;
  }
  else
    error("lt: Incompatible for comparison.\n");
  
  return C;
}


MATRIX *opr_le(A, B)
     MATRIX *A, *B;
{
  MATRIX *C;

  int i;
  
  int nrowa = NROW(A), ncola = NCOL(A);
  int nrowb = NROW(B), ncolb = NCOL(B);

  double *a = MATR(A), *b = MATR(B), *c;

  if ((nrowa == 1) && (ncola == 1))
  {
    C = mat_new(TYPE(B),nrowb,ncolb); c = MATR(C);
    nrowb *= ncolb;
    for(i = 0; i < nrowb; i++,c++) if (*a <= b[i]) *c = 1;
  }
  
  else if ((nrowb == 1) && (ncolb == 1))
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++,c++) if (a[i] <= *b) *c = 1;
  }
  
  else if ((nrowa == nrowb) && (ncola == ncolb))
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++,c++) if (a[i] <= b[i]) *c = 1;
  }
  else
    error("le: Incompatible for comparison.\n");
  
  return C;
}


MATRIX *opr_gt(A, B)
     MATRIX *A, *B;
{
  MATRIX *C;
  int i;
  
  int nrowa = NROW(A), ncola = NCOL(A);
  int nrowb = NROW(B), ncolb = NCOL(B);

  double *a = MATR(A), *b = MATR(B), *c;

  if ((nrowa == 1) && (ncola == 1))
  {
    C = mat_new(TYPE(B),nrowb,ncolb); c = MATR(C);
    nrowb *= ncolb;
    for(i = 0; i < nrowb; i++,c++) if (*a > b[i]) *c = 1;
  }
  
  else if ((nrowb == 1) && (ncolb == 1))
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++,c++) if (a[i] > *b) *c = 1;
  }
  
  else if ((nrowa == nrowb) && (ncola == ncolb))
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++,c++) if (a[i] > b[i]) *c = 1;
  }
  
  else
    error("gt: Incompatible for comparison.\n");
  
  return C;
}


MATRIX *opr_ge(A, B)
     MATRIX *A, *B;
{
  MATRIX *C;
  int i;
  
  int nrowa = NROW(A), ncola = NCOL(A);
  int nrowb = NROW(B), ncolb = NCOL(B);

  double *a = MATR(A), *b = MATR(B), *c;

  if ((nrowa == 1) && (ncola == 1))
  {
    C = mat_new(TYPE(B),nrowb,ncolb); c = MATR(C);
    nrowb *= ncolb;
    for(i = 0; i < nrowa; i++,c++) if (*a >= b[i]) *c = 1;
  }
  
  else if ((nrowb == 1) && (ncolb == 1))
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++,c++) if (a[i] >= *b) *c = 1;
  }
  
  else if ((nrowa == nrowb) && (ncola == ncolb))
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++,c++) if (a[i] >= b[i]) *c = 1;
  }
  
  else
    error("ge: Incompatible for comparison.\n");
  
  return C;
}


MATRIX *opr_eq(A, B)
     MATRIX *A, *B;
{
  MATRIX *C;
  int i;
  
  int nrowa = NROW(A), ncola = NCOL(A);
  int nrowb = NROW(B), ncolb = NCOL(B);

  double *a = MATR(A), *b = MATR(B), *c;

  if ((nrowa == 1) && (ncola == 1))
  {
    C = mat_new(TYPE(B),nrowb,ncolb); c = MATR(C);
    nrowb *= ncolb;
    for(i = 0; i < nrowb; i++,c++) if (*a == b[i]) *c = 1;
  }
  
  else if ((nrowb == 1) && (ncolb == 1))
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++,c++) if (a[i] == *b) *c = 1;
  }
  
  else if ((nrowa == nrowb) && (ncola == ncolb))
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++,c++) if (a[i] == b[i]) *c = 1;
  }
  
  else 
    error("eq: Incompatible for comparison.\n");
  
  return C;
}


MATRIX *opr_neq(A, B)
     MATRIX *A, *B;
{
  MATRIX *C;
  int i;
  
  int nrowa = NROW(A), ncola = NCOL(A);
  int nrowb = NROW(B), ncolb = NCOL(B);

  double *a = MATR(A), *b = MATR(B), *c;

  if ((nrowa == 1) && (ncola == 1))
  {
    C = mat_new(TYPE(B),nrowb,ncolb); c = MATR(C);
    nrowb *= ncolb;
    for(i = 0; i < nrowb; i++,c++) if (*a != b[i]) *c = 1;
  }
  
  else if ((nrowb == 1) && (ncolb == 1))
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++,c++) if (a[i] != *b) *c = 1;
  }
  
  else if ((nrowa == nrowb) && (ncola == ncolb))
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++,c++) if (a[i] != b[i]) *c = 1;
  }
  
  else 
    error("neq: Incompatible for comparison.\n");
  
  return C;
}

MATRIX *opr_or(A, B)
     MATRIX *A, *B;
{
  MATRIX *C;
  int i;
  
  int nrowa = NROW(A), ncola = NCOL(A);
  int nrowb = NROW(B), ncolb = NCOL(B);

  double *a = MATR(A), *b = MATR(B), *c;

  if ((nrowa == 1) && (ncola == 1))
  {
    C = mat_new(TYPE(B),nrowb,ncolb); c = MATR(C);
    nrowb *= ncolb;
    for(i = 0; i < nrowb; i++) *c++ = (*a) || (b[i]);
  }
  
  else if ((nrowb == 1) && (ncolb == 1))
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++) *c++ = (a[i]) || (*b);
  }
  
  else if ((nrowa == nrowb) && (ncola == ncolb))
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++) *c++ = (a[i]) || (b[i]);
  }
  
  else 
    error("or: Incompatible for comparison.\n");
  
  return C;
}

MATRIX *opr_and(A, B)
     MATRIX *A, *B;
{
  MATRIX *C;
  int i;
  
  int nrowa = NROW(A), ncola = NCOL(A);
  int nrowb = NROW(B), ncolb = NCOL(B);

  double *a = MATR(A), *b = MATR(B), *c;

  if ((nrowa == 1) && (ncola == 1))
  {
    C = mat_new(TYPE(B),nrowb,ncolb); c = MATR(C);
    nrowb *= ncolb;
    for(i = 0; i < nrowb; i++) *c++ = (*a) && (b[i]);
  }
  
  else if ((nrowb == 1) && (ncolb == 1))
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++) *c++ = (a[i]) && (*b);
  }
  
  else if ((nrowa == nrowb) && (ncola == ncolb))
  {
    C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
    nrowa *= ncola;
    for(i = 0; i < nrowa; i++) *c++ = (a[i]) && (b[i]);
  }
  
  else 
    error("and: Incompatible for comparison.\n");
  
  return C;
}

MATRIX *opr_not(A)
     MATRIX *A;
{
  MATRIX *C;
  int i;

  int nrowa = NROW(A), ncola = NCOL(A);

  double *a = MATR(A), *c;

  C = mat_new(TYPE(A),nrowa,ncola); c = MATR(C);
  nrowa *= ncola;
  for(i = 0; i < nrowa; i++,c++) if (*a == 0) *c = 1;

  return C;
}
