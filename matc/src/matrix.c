/*****************************************************************************
 *
 *  Elmer, A Finite Element Software for Multiphysical Problems
 *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library (in file ../LGPL-2.1); if not, write 
 * to the Free Software Foundation, Inc., 51 Franklin Street, 
 * Fifth Floor, Boston, MA  02110-1301  USA
 *
 *****************************************************************************/

/*******************************************************************************
 *
 *     MATC matrix utilities.
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
|  MATRIX.C - Last Edited 8. 8. 1988
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
 * $Id: matrix.c,v 1.1.1.1 2005/04/14 13:29:14 vierinen Exp $ 
 *
 * $Log: matrix.c,v $
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:50  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"

#define MA(i,j) a[(i) * ncola + (j)]
#define MB(i,j) b[(i) * ncolb + (j)]
#define MC(i,j) c[(i) * ncolc + (j)]

double func_abs(arg) 
     double arg;
{
  return abs(arg);
}

double func_mod(x,y)
     double x,y;
{
  int ix, iy;

  ix = x + 0.5;
  iy = y + 0.5;
  return (double)(ix % iy);
}

VARIABLE *mtr_sum(A) VARIABLE *A;
{
   VARIABLE *C;

   int i, j;

   int nrowa = NROW(A), ncola = NCOL(A);

   double *a = MATR(A), *c;


   if (nrowa == 1 || ncola == 1)
   {
     C = var_temp_new(TYPE_DOUBLE, 1, 1); c = MATR(C);
     nrowa = (nrowa == 1) ? ncola : nrowa;
     for(i = 0; i < nrowa; i++) *c += *a++;
   }
   else
   {
     C = var_temp_new(TYPE_DOUBLE, 1, ncola); c = MATR(C);
     for(i = 0; i < ncola; i++) 
       for(j = 0; j < nrowa; j++) c[i] += MA(j, i);
   }

   return C;
}

VARIABLE *mtr_trace(A) VARIABLE *A;
{
  VARIABLE *C;

  double temp = 0.0;
  int i;

  int nrowa = NROW(A), ncola = NCOL(A);
  double *a = MATR(A);
  
  if (nrowa !=  ncola) error("trace: not square.\n");

  for(i = 0; i < nrowa; i++) temp += MA(i,i);

  C = var_temp_new(TYPE(A), 1, 1); *MATR(C) = temp;

  return C;
}

VARIABLE *mtr_zeros(A) VARIABLE *A;
{
  VARIABLE *C;

  int ind1 = 1, ind2 = 1;
  
  if (NEXT(A) != NULL)
  {
   ind1 = (int)*MATR(A); ind2 = (int)*MATR(NEXT(A));
  } 
  else
  {
    ind2 = (int)*MATR(A);
  }

  if (ind1 < 1 || ind2 < 1) 
    error("Zeros: invalid size for and array");

  C = var_temp_new(TYPE_DOUBLE, ind1, ind2);

  return C;
}

VARIABLE *mtr_ones(A) VARIABLE *A;
{
  VARIABLE *C;
  double *c;
  int i, n;

  C = mtr_zeros(A); c = MATR(C);
  n = NROW(C) * NCOL(C);

  for(i = 0; i < n; i++) *c++ = 1.0; 

  return C;
}

VARIABLE *mtr_rand(A) VARIABLE *A;
{
  VARIABLE *C;

  static int seed = 0;
#pragma omp threadprivate (seed)
  int i, n;

  double *c;

  C = mtr_zeros(A); c = MATR(C);
  n = NROW(C) * NCOL(C);

  if (seed == 0) seed = time(NULL);
  for(i = 0; i < n; i++)  *c++ = urand(&seed);

  return C;
}

VARIABLE *mtr_resize(A) VARIABLE *A;
{
  VARIABLE *C;

  int i, j, n, m, ind1 = 1, ind2;

  double *a = MATR(A), *c;
  
  if (NEXT(NEXT(A)) != NULL)
  {
    ind1 = *MATR(NEXT(A)); ind2 = *MATR(NEXT(NEXT(A)));
  } 
  else
  {
    ind2 = (int)*MATR(NEXT(A));;
  }

  if (ind1 < 1 || ind2 < 1) 
    error("resize: invalid size for and array");
  
  C = var_temp_new(TYPE(A), ind1, ind2); c = MATR(C);
  a = MATR(A); n = ind1 * ind2; m = NROW(A) * NCOL(A);

  for(i = j = 0; i < n; i++)
  {
    *c++ = a[j++]; if (j == m) j = 0;
  }

  return C;
}

VARIABLE *mtr_vector(A) VARIABLE *A;
{
  VARIABLE *C;

  double start, stop, incr, x, *c;

  int i, eval;
  
  start = *MATR(A); stop =  *MATR(NEXT(A));

  if (NEXT(NEXT(A)) != (VARIABLE *)NULL)
    incr = *MATR(NEXT(NEXT(A)));
  else
    incr = (start < stop) ? (1) : (-1);
  
  if (incr == 0)
    incr = (start < stop) ? (1) : (-1);

  eval = (int)(abs(stop-start) / abs(incr)) + 1;
  if (eval < 1) return NULL;

  C = var_temp_new(TYPE_DOUBLE, 1, eval); c = MATR(C);
  
  x = start;
  for(i = 0;  i < eval; i++)
  {
    *c++ = x; x += incr;
  }
  
  return C;
}

VARIABLE *mtr_eye(A) VARIABLE *A;
{
  VARIABLE *C;
  double *c;
  int i, ind, ncolc;
  
  if (*MATR(A) < 1)
  {
    error("eye: Invalid size for an array.\n");
  }

  ind = (int)*MATR(A);
  
  C = var_temp_new(TYPE_DOUBLE, ind, ind); 
  c = MATR(C); ncolc = ind;

  for (i = 0; i < ind; i++) MC(i,i) = 1.0;
  
  return C;
}

VARIABLE *mtr_size(A) VARIABLE *A;
{  
  VARIABLE *C;
  double *c;

  C = var_temp_new(TYPE_DOUBLE, 1, 2); c = MATR(C);
  *c++ = NROW(A); *c = NCOL(A);
  
  return C;
}

VARIABLE *mtr_min(A) VARIABLE *A;
{
   VARIABLE *C;

   double *a = MATR(A), *c;

   int nrowa = NROW(A), ncola = NCOL(A);
   int i, j;

   if (nrowa == 1 || ncola == 1)
   {
     C = var_temp_new(TYPE_DOUBLE, 1, 1); c = MATR(C);
     *c = *a++; nrowa = max(ncola, nrowa);
     for(i = 1; i < nrowa; i++, a++) *c = min(*c, *a);
   }
   else
   {
     C = var_temp_new(TYPE_DOUBLE, 1, ncola); c = MATR(C);
     for(i = 0; i < ncola; i++, c++)
     {
       *c = MA(0, i);
       for(j = 1; j < nrowa; j++) *c = min(*c, MA(j, i));
     }
   }

   return C;
}

VARIABLE *mtr_max(A) VARIABLE *A;
{
   VARIABLE *C;

   double *a = MATR(A), *c;

   int nrowa = NROW(A), ncola = NCOL(A);
   int i, j;

   if (nrowa == 1 || ncola == 1)
   {
     C = var_temp_new(TYPE_DOUBLE, 1, 1); c = MATR(C);
     *c = *a++; nrowa = max(ncola, nrowa);
     for(i = 1; i < nrowa; i++, a++) *c = max(*c, *a);
   }
   else
   {
     C = var_temp_new(TYPE_DOUBLE, 1, ncola); c = MATR(C);
     for(i = 0; i < ncola; i++, c++)
     {
       *c = MA(0, i);
       for(j = 1; j < nrowa; j++) *c = max(*c, MA(j, i));
     }
   }

   return C;
}

VARIABLE *mtr_diag(A) VARIABLE *A;
{
   VARIABLE *C;

   double *a = MATR(A), *c;

   int nrowa = NROW(A), ncola = NCOL(A);
   int ncolc;

   int i;

   if (nrowa == 1 || ncola == 1)
   {
     nrowa = max(nrowa, ncola); ncolc = nrowa;
     C = var_temp_new(TYPE_DOUBLE, nrowa, nrowa); c = MATR(C);
     for(i = 0; i < nrowa; i++) MC(i, i) = *a++;
   }
   else
   {
     C = var_temp_new(TYPE_DOUBLE, 1, nrowa); c = MATR(C);
     for(i = 0; i < min(nrowa,ncola); i++) *c++ = MA(i, i);
   }

   return C;
}

VARIABLE *mtr_pow(A) VARIABLE *A;
{
   VARIABLE *B = NEXT(A), *C;
   double *a = MATR(A), b = M(B,0,0), *c;

   int nrowa = NROW(A), ncola = NCOL(A);

   int i;

   C = var_temp_new(TYPE_DOUBLE, nrowa,ncola );
   c = MATR( C );
   for(i = 0; i < nrowa*ncola; i++) *c++ = pow(*a++,b);

   return C;
}

VARIABLE *mtr_where(A) VARIABLE *A;
{
   VARIABLE *C;

   double *a = MATR(A), *c;

   int nrowa = NROW(A), ncola = NCOL(A);

   int i,n=0;

   for( i=0; i < nrowa*ncola; i++) if ( a[i] ) n++;

   C = var_temp_new( TYPE_DOUBLE,1,n );
   c = MATR(C);

   for( i=0; i < nrowa*ncola; i++ ) if ( a[i] ) { *c++ = i; }

   return C;
}

void mtr_com_init()
{
  static char *minHelp =
  {
     "r = min(matrix)\n"
     "Return value is a vector containing smallest element in columns of given matrix.\n"
     "r=min(min(matrix) gives smallest element of the matrix.\n\n"
  };

  static char *maxHelp =
  {
     "r = max(matrix)\n"
     "Return value is a vector containing largest element in columns of given matrix.\n"
     "r=max(max(matrix)) gives largest element of the matrix.\n\n"
  };

  static char *sumHelp =
  {
     "r = sum(matrix)\n"
     "Return vector is column sums of given matrix. r=sum(sum(matrix)) gives\n"
     "the total sum of elements of the matrix.\n\n"
  };

  static char *traceHelp =
  {
     "r = trace(matrix)\n"
     "Return value is sum of matrix diagonal elements.\n\n"
  };

  static char *detHelp =
  {
     "r = det(matrix)\n"
     "Return value is determinant of given square matrix.\n\n"
  };

  static char *invHelp =
  {
     "r = inv(matrix)\n"
     "Invert given square matrix. Computed also by r=matrix^(-1).\n\n"
  };

  static char *eigHelp =
  {
     "r = eig(matrix)\n"
     "Return eigenvalues of given square matrix. r(n,0) is real part of the\n"
     "n:th eigenvalue, r(n,1) is the imaginary part respectively\n\n"
  };

  static char *jacobHelp =
  {
     "r = jacob(a,b,eps)\n"
     "Solve symmetric positive definite eigenvalue problem by Jacob iteration.\n"
     "Return values are the eigenvalues. Also a variable eigv is created containing\n"
     "eigenvectors.\n\n"
  };

  static char *ludHelp =
  {
     "r = lud(matrix)\n"
     "Return value is lud decomposition of given matrix.\n\n"
  };

  static char *hesseHelp =
  {
     "r = hesse(matrix)\n"
     "Return the upper hessenberg form of given matrix.\n\n"
  };

  static char *eyeHelp =
  {
     "r = eye(n)\n"
     "Return n by n identity matrix.\n\n"
  };

  static char *zerosHelp =
  {
     "r = zeros(n,m)\n"
     "Return n by m matrix with elements initialized to zero.\n"
  };

  static char *onesHelp =
  {
     "r = ones(n,m)\n"
     "Return n by m matrix with elements initialized to one.\n"
  };

  static char *randHelp =
  {
     "r = rand(n,m)\n"
     "Return n by m matrix with elements initialized to with random number from\n"
     "zero to one.\n\n"
  };

  static char *diagHelp =
  {
     "r=diag(matrix) or r=diag(vector)\n"
     "Given matrix return diagonal entries as a vector. Given vector return matrix\n"
     "with diagonal elements from vector. r=diag(diag(a)) gives matrix with diagonal\n"
     "elements from matrix a otherwise elements are zero.\n\n"
  };

  static char *vectorHelp =
  {
     "r=vector(start,end,inc)\n"
     "Return vector of values going from start to end by inc.\n\n"
  };

  static char *sizeHelp =
  {
     "r = size(matrix)\n"
     "Return size of given matrix.\n"
  };

  static char *resizeHelp =
  {
     "r = resize(matrix,n,m)\n"
     "Make a matrix to look as a n by m matrix.\n\n"
  };

  static char *whereHelp =
  {
     "r = where(l)\n"
     "Return linear indexes of where l is true.\n\n"
  };

  com_init( "sin"    , TRUE,  TRUE,  (VARIABLE *(*)())sin    , 1, 1, "r=sin(x)" );
  com_init( "cos"    , TRUE,  TRUE,  (VARIABLE *(*)())cos    , 1, 1, "r=cos(x)" );
  com_init( "tan"    , TRUE,  TRUE,  (VARIABLE *(*)())tan    , 1, 1, "r=tan(x)" );
  com_init( "asin"   , TRUE,  TRUE,  (VARIABLE *(*)())asin   , 1, 1, "r=asin(x)" );
  com_init( "acos"   , TRUE,  TRUE,  (VARIABLE *(*)())acos   , 1, 1, "r=acos(x)" );
  com_init( "atan"   , TRUE,  TRUE,  (VARIABLE *(*)())atan   , 1, 1, "r=atan(x)" );
  com_init( "atan2"  , TRUE,  TRUE,  (VARIABLE *(*)())atan2  , 2, 2, "r=atan2(y,x)" );
  com_init( "sinh"   , TRUE,  TRUE,  (VARIABLE *(*)())sinh   , 1, 1, "r=sinh(x)" );
  com_init( "cosh"   , TRUE,  TRUE,  (VARIABLE *(*)())cosh   , 1, 1, "r=cosh(x)" );
  com_init( "tanh"   , TRUE,  TRUE,  (VARIABLE *(*)())tanh   , 1, 1, "r=tanh(x)" );
  com_init( "exp"    , TRUE,  TRUE,  (VARIABLE *(*)())exp    , 1, 1, "r=exp(x)" );
  com_init( "ln"     , TRUE,  TRUE,  (VARIABLE *(*)())log    , 1, 1, "r=ln(x)\nNatural logarithm." );
  com_init( "log"    , TRUE,  TRUE,  (VARIABLE *(*)())log10  , 1, 1, "r=log(x)\nBase 10 logarithm." );
  com_init( "sqrt"   , TRUE,  TRUE,  (VARIABLE *(*)())sqrt   , 1, 1, "r=sqrt(x)" );
  com_init( "ceil"   , TRUE,  TRUE,  (VARIABLE *(*)())ceil   , 1, 1, "r=ceil(x)\nSmallest integer not less than x." );
  com_init( "floor"  , TRUE,  TRUE,  (VARIABLE *(*)())floor  , 1, 1, "r=floor(x)\nLargest integer not more than x." );
  com_init( "abs"    , TRUE,  TRUE,  (VARIABLE *(*)())func_abs   , 1, 1,"r=abs(x)");
  com_init( "mod"    , TRUE,  TRUE,  (VARIABLE *(*)())func_mod   , 2, 2,"r=mod(x,y)");
  com_init( "pow"    , FALSE, TRUE,  mtr_pow,     2, 2, "r=pow(x,y)" );
  com_init( "min"    , FALSE, TRUE,  mtr_min,     1, 1, minHelp    );
  com_init( "max"    , FALSE, TRUE,  mtr_max,     1, 1, maxHelp    );
  com_init( "sum"    , FALSE, TRUE,  mtr_sum,     1, 1, sumHelp    );
  com_init( "trace"  , FALSE, TRUE,  mtr_trace,   1, 1, traceHelp  );
  com_init( "det"    , FALSE, TRUE,  mtr_det,     1, 1, detHelp    );
  com_init( "inv"    , FALSE, TRUE,  mtr_inv,     1, 1, invHelp    );
  com_init( "eig"    , FALSE, TRUE,  mtr_eig,     1, 1, eigHelp    );
  com_init( "jacob"  , FALSE, TRUE,  mtr_jacob,   3, 3, jacobHelp  );
  com_init( "lud"    , FALSE, TRUE,  mtr_LUD,     1, 1, ludHelp    );
  com_init( "hesse"  , FALSE, TRUE,  mtr_hesse,   1, 1, hesseHelp  );
  com_init( "eye"    , FALSE, TRUE,  mtr_eye,     1, 1, eyeHelp    );
  com_init( "zeros"  , FALSE, TRUE,  mtr_zeros,   1, 2, zerosHelp  );
  com_init( "ones"   , FALSE, TRUE,  mtr_ones,    1, 2, onesHelp   );
  com_init( "rand"   , FALSE, FALSE, mtr_rand,    1, 2, randHelp   );
  com_init( "diag"   , FALSE, TRUE,  mtr_diag,    1, 1, diagHelp   );
  com_init( "vector" , FALSE, TRUE,  mtr_vector,  2, 3, vectorHelp );
  com_init(  "size"  , FALSE, TRUE,  mtr_size,    1, 1, sizeHelp  );
  com_init( "resize" , FALSE, TRUE,  mtr_resize,  2, 3, resizeHelp );
  com_init( "where"  , FALSE, FALSE, mtr_where,   1, 1, whereHelp );
}
