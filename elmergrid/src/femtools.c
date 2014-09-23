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


/* -------------------------------:  FEMTOOLS.C  :----------------------------

   Basic Linear Algebra Subroutines, Gaussian quadratures and such. These subroutines
   should be very general at least in the Euclidian parts of the universe.
*/

#include <stdio.h>
#include <math.h>
#include "common.h"
#include "femdef.h"
#include "femtools.h"

#if 0

void Squad404(Real *xi,Real *eta,Real *sfun, Real *lder)
/*     Computes shape function and its derivatives in local 
       coordinates for a 4-node quadrilatelar element. 
       Input: 
              xi      - x coordinate of the point where functions are required 
              eta     - y coordinate of the point where functions are required 
       Output: 
              sfun    - Shape function values 
              lder    - Partial derivatives of shape functions respect to 
                        xi (first line) and eta (second line). */
{
  sfun[BOTLEFT]  = (1. - (*xi)) * (1. - (*eta)) / 4.;
  sfun[BOTRIGHT] = ((*xi) + 1.) * (1. - (*eta)) / 4.;
  sfun[TOPRIGHT] = ((*xi) + 1.) * ((*eta) + 1.) / 4.;
  sfun[TOPLEFT]  = (1. - (*xi)) * ((*eta) + 1.) / 4.;
  
  lder[BOTLEFT]  = ((*eta) - 1. ) / 4.;
  lder[BOTRIGHT] = (1. - (*eta) ) / 4.;
  lder[TOPRIGHT] = (1. + (*eta) ) / 4.;
  lder[TOPLEFT]  = (-(*eta) - 1.) / 4.;
  
  lder[BOTLEFT+4]  = ((*xi) - 1. ) / 4.;
  lder[BOTRIGHT+4] = (-(*xi) - 1.) / 4.;
  lder[TOPRIGHT+4] = ((*xi) + 1. ) / 4.;
  lder[TOPLEFT+4]  = (1. - (*xi) ) / 4.;
} /* Squad404 */


void Squad303(Real *xi,Real *eta,Real *sfun, Real *lder)
/*     Computes shape function and its derivatives in local 
       coordinates for a 4-node quadrilatelar element. 
       Input: 
              xi      - x coordinate of the point where functions are required 
              eta     - y coordinate of the point where functions are required 
       Output: 
              sfun    - Shape function values 
              lder    - Partial derivatives of shape functions respect to 
                        xi (first line) and eta (second line). */
{
  sfun[BOTLEFT]  = 1.0-(*xi)-(*eta);
  sfun[BOTRIGHT] = (*xi);
  sfun[TOPRIGHT] = (*eta);
  
  lder[BOTLEFT]  = -1.0;
  lder[BOTRIGHT] = 1.0;
  lder[TOPRIGHT] = 0.0;
  
  lder[BOTLEFT+3]  = -1.0;
  lder[BOTRIGHT+3] = 0.0;
  lder[TOPRIGHT+3] = 1.0;
} /* Squad303 */



void Squad408(Real *xi,Real *eta,Real *sfun, Real *lder)
/*     Computes shape function and its derivatives in local 
       coordinates for a 4-node quadrilatelar element. 
       Input: 
              xi      - x coordinate of the point where functions are required 
              eta     - y coordinate of the point where functions are required 
       Output: 
              sfun    - Shape function values 
              lder    - Partial derivatives of shape functions respect to 
                        xi (first line) and eta (second line). */
{
  sfun[0] = (1.-(*xi)) * (1.-(*eta)) * (-(*xi)-(*eta)-1) / 4.;
  sfun[1] = (1.+(*xi)) * (1.-(*eta)) * ((*xi)-(*eta)-1) / 4.;
  sfun[2] = (1.+(*xi)) * (1.+(*eta)) * ((*xi)+(*eta)-1) / 4.;
  sfun[3] = (1.-(*xi)) * (1.+(*eta)) * (-(*xi)+(*eta)-1) / 4.;
  sfun[4] = (1.-(*xi)*(*xi)) * (1.-(*eta)) / 2.;
  sfun[5] = (1.-(*eta)*(*eta)) * (1.+(*xi)) / 2.;
  sfun[6] = (1.-(*xi)*(*xi)) * (1.+(*eta)) / 2.;
  sfun[7] = (1.-(*eta)*(*eta)) * (1.-(*xi)) / 2.;
  
  lder[0] = (1.-(*eta))*(2*(*xi)+(*eta))/4.; 
  lder[1] = (1.-(*eta))*(-2*(*xi)+(*eta))/4.; 
  lder[2] = (1.+(*eta))*(-2*(*xi)-(*eta))/4.; 
  lder[3] = (1.+(*eta))*(2*(*xi)-(*eta))/4.; 
  lder[4] = -(*xi)*(1.-(*eta));
  lder[5] = (1.-(*eta)*(*eta))/2.0;
  lder[6] = -(*xi)*(1.+(*eta));
  lder[7] = -(1.-(*eta)*(*eta))/2.0;
  
  lder[8] = (1.-(*xi))*(2*(*eta)+(*xi))/4.; 
  lder[9] = (1.+(*xi))*(2*(*eta)-(*xi))/4.; 
  lder[10] = (1.+(*xi))*(-2*(*eta)-(*xi))/4.; 
  lder[11] = (1.-(*xi))*(-2*(*eta)+(*xi))/4.; 
  lder[12] = -(1.-(*xi)*(*xi))/2.0;
  lder[13] = -(*eta)*(1.+(*xi));
  lder[14] = -(1.-(*xi)*(*xi))/2.0;
  lder[15] = -(*eta)*(1.-(*xi));
} /* Squad408 */


void Squad409(Real *xi,Real *eta,Real *sfun, Real *lder)
{
  printf("Squad409: Implement this in order to use it!\n");
}



void Squad202(Real *xi, Real *sfun, Real *lder)
/*     Computes shape function and its derivatives in local 
       coordinates for a 2-node linear element. 
       Input: 
              xi      - coordinate of the point where functions are required 
       Output: 
              sfun    - Shape function values 
              lder    - Partial derivatives of shape functions respect to 
                        xi (first line) and eta (second line). */
{
  sfun[0]  = (1. - *xi) / 2.;
  sfun[1] = (*xi + 1.) / 2.;
  
  lder[0]  = -0.5;
  lder[1] =  0.5;
  
} /* Squad2 */


void Squad203(Real *xi, Real *sfun, Real *lder)
/*     Computes shape function and its derivatives in local 
       coordinates for a 3-node quadratic element. */
{
  sfun[0] = -(*xi) * (1. - (*xi)) / 2.;
  sfun[1] = (1 - (*xi)*(*xi)); 
  sfun[2] = (*xi) * ((*xi) + 1.) / 2.;
  
  lder[0] = -0.5 - (*xi);
  lder[1] = -2*(*xi);
  lder[2] = 0.5 + (*xi);
} /* Squad3 */



int LocalToGlobalD2(Real *globalcoord,Real *shapeder,
		    Real *shapefunc,Real *xgauss, Real *ygauss, 
		    Real *det,Real *globalder,int nodesd2)
{
  int i;
  Real jacob00,jacob01,jacob10,jacob11;
  Real x,y;

  /* Calculate the global coordinates corresponding to 
     the local point for gaussian integration */

  (*xgauss) = 0.;  
  (*ygauss) = 0.; 
  
  jacob00 = 0.0;
  jacob01 = 0.0;
  jacob10 = 0.0;
  jacob11 = 0.0;

  /* Calculate the Jacobian matrix */
  for(i = 0; i < nodesd2; i++) {
    x = globalcoord[i];
    y = globalcoord[i+nodesd2];
    
    *xgauss += x * shapefunc[i];
    jacob00 += x * shapeder[i];
    jacob10 += x * shapeder[i+nodesd2];
    
    *ygauss += y * shapefunc[i];
    jacob01 += y * shapeder[i];
    jacob11 += y * shapeder[i+nodesd2];
  }
  
  /* Inverse the jacobian */
  (*det) = jacob11*jacob00 - jacob01*jacob10;
  if((*det) < 0) {
    printf("nodesd2 = %d\n  gauss=[%.4g  %.4g]\n",nodesd2,*xgauss,*ygauss);
    for(i=0;i<2*nodesd2;i++)
      printf("%.5g  ",globalcoord[i]);
    printf("\n");       
    bigerror("Isoparametric transformation is not possible.");
  }

  /* Calculate the derivatives of the global coordinates */
  for(i = 0; i < nodesd2; i++) {
    globalder[i] = ( jacob11 * shapeder[i] - jacob01 * shapeder[i+nodesd2] ) / (*det);  
    globalder[i+nodesd2] = ( jacob00 * shapeder[i+nodesd2] - jacob10 * shapeder[i] ) / (*det);  
  }
  return(0);
}




void LocalToGlobalD1(Real *globalcoord,Real *shapeder,
		     Real *shapefunc,Real *xgauss, Real *ygauss, 
		     Real *ratio,int nodesd1)

/* This subroutine makes the coordinate tranformation required 
   for isoparametric linear elements. This subroutine is used 
   when calculating the boundary conditions for a case with 
   bilinear elements. The number of Gaussian quadrature points 
   is independent of this subroutine.
   Input:   globalcoord  - corner coordinates of the element in real space
            shapefunc    - values of the shape functions at the integration point given
	    shapeder     - values of the derivatives of the shape functions 
   Output:  xgauss       - global x-coordinate corresponding to the integration point     
            ygauss       - global y-coordinate corresponding to the integration point
            ratio        - the ratio of the side lenght in local and global coordinates 
	    */
{
  int i;
  Real jacob0=0.0,jacob1=0.0;

  /* Calculate the global coordinates corresponding to the 
     local point for gaussian integration */
  (*xgauss) = 0.;  (*ygauss) = 0.; 

  for(i = 0; i < nodesd1; i++) {
    *xgauss += globalcoord[i] * shapefunc[i];
    *ygauss += globalcoord[i+nodesd1] * shapefunc[i];
  }

  /* Calculate the Jacobian vector */
  for(i = 0; i < nodesd1; i++) {
    jacob0 += shapeder[i] * globalcoord[i];
    jacob1 += shapeder[i] * globalcoord[i+nodesd1];
  }

  (*ratio) = sqrt( jacob0*jacob0 + jacob1*jacob1 );
}



void SurfaceNormalD1(Real *coord,Real *normal,int nodesd1)
/* Calculates the normal of the surface for 2-point
   elements. The direction is such that it points 
   to the right. */
{
  Real dx,dy,ds;

  dx = coord[1]-coord[0];
  dy = coord[nodesd1+1]-coord[nodesd1];
  ds = sqrt(dx*dx+dy*dy);
  normal[0] = dy/ds;  
  normal[1] = -dx/ds;   
}



int GlobalToLocalD2(Real *coord,Real xglobal,Real yglobal,
		    Real *xlocal,Real *ylocal)
/* Given a global element and a global point calculates the local 
   coordinates of the particular point in the given element. 
   */
{
  int hit;
  int nodesd2;
  Real a0,a1,a2,a3,b0,b1,b2,b3;
  Real ratio=1.0e-2;
 
  nodesd2 = 4;

  a0 = 4*xglobal-coord[0]-coord[1]-coord[2]-coord[3];
  a1 = coord[2]+coord[3]-coord[0]-coord[1];
  a2 = coord[1]+coord[2]-coord[0]-coord[3];
  a3 = coord[0]+coord[2]-coord[1]-coord[3];

  b0 = 4*yglobal-coord[nodesd2]-coord[1+nodesd2]-coord[2+nodesd2]-coord[3+nodesd2];
  b1 = coord[2+nodesd2]+coord[3+nodesd2]-coord[nodesd2]-coord[1+nodesd2];
  b2 = coord[1+nodesd2]+coord[2+nodesd2]-coord[nodesd2]-coord[3+nodesd2];
  b3 = coord[nodesd2]+coord[2+nodesd2]-coord[1+nodesd2]-coord[3+nodesd2];

  if(fabs(a1/a2)<ratio && fabs(a3/a2)<ratio) {
    (*xlocal) = a0/a2;
    (*ylocal) = (b0-b2*(*xlocal))/(b1+b3*(*xlocal));
    hit = 1;
  }
  else if(fabs(b2/b1)<ratio && fabs(b3/b1)<ratio) {
    (*ylocal) = b0/b1;
    (*xlocal) = (a0-a1*(*ylocal))/(a2+a3*(*ylocal));
    hit = 2;
  }
  else {
    hit = FALSE;
    return(hit);
  }

  if(fabs(*xlocal) > 1.0) {
    if(fabs(*xlocal) < 1.01) 
      (*xlocal) /= fabs(*xlocal);
    else 
      hit = FALSE;
  }

  if(fabs(*ylocal) > 1.0) {
    if(fabs(*ylocal) < 1.01) 
      (*ylocal) /= fabs(*ylocal);
    else 
      hit = FALSE;
  }

  if(0) printf("(x,y)=(%.2g,%.2g)\n",(*xlocal),(*ylocal));
  return(hit);
}




int GlobalToLocalD2Complicated(Real *coord,Real xglobal,Real yglobal,
		    Real *xlocal,Real *ylocal)
/* Given a global element and a global point calculates the local 
   coordinates of the particular point in the given element. 
   */
{
  int hit;
  int nodesd2;
  Real a0,a1,a2,a3,b0,b1,b2,b3;
  Real ratio=1.0e-2;
 
  nodesd2 = 4;

  a0 = 4*xglobal-coord[0]-coord[1]-coord[2]-coord[3];
  a1 = coord[2]+coord[3]-coord[0]-coord[1];
  a2 = coord[1]+coord[2]-coord[0]-coord[3];
  a3 = coord[0]+coord[2]-coord[1]-coord[3];

  b0 = 4*yglobal-coord[nodesd2]-coord[1+nodesd2]-coord[2+nodesd2]-coord[3+nodesd2];
  b1 = coord[2+nodesd2]+coord[3+nodesd2]-coord[nodesd2]-coord[1+nodesd2];
  b2 = coord[1+nodesd2]+coord[2+nodesd2]-coord[nodesd2]-coord[3+nodesd2];
  b3 = coord[nodesd2]+coord[2+nodesd2]-coord[1+nodesd2]-coord[3+nodesd2];

  if(fabs(a1/a2)<ratio && fabs(a3/a2)<ratio) {
    (*xlocal) = a0/a2;
    (*ylocal) = (b0-b2*(*xlocal))/(b1+b3*(*xlocal));
    hit = 1;
  }
  else if(fabs(b2/b1)<ratio && fabs(b3/b1)<ratio) {
    (*ylocal) = b0/b1;
    (*xlocal) = (a0-a1*(*ylocal))/(a2+a3*(*ylocal));
    hit = 2;
  }
  else {
    (*xlocal) = 0.0;
    (*ylocal) = 0.0;
    hit = FALSE;
#if 0
    xtest[0] = 0.0;
    xtest[1] = 0.0;
    ytest = 0.0;

    a = a3*b2-a2*b3;
    b = a1*b2-a2*b1+a0*b3-a3*b0;
    c = a0*b1-a1*b0;
    d = b*b-4*a*c;
    if(d < 0) return(FALSE);

    if (fabs(a) < tiny) xtest[0] = -b/c;
    else {
      d = sqrt(d);
      xtest[0] = (-b-d)/(2*a);
      xtest[1] = (-b+d)/(2*a);
    }
    
    for(i=0;i<2;i++) 
      if(fabs(xtest[i]) <= 1.0+tiny) {
	ytest = (b0-b2*xtest[i])/(b1+b3*xtest[i]);
	if(fabs(ytest) <= 1.0+tiny) {
/*	   fabs(a1*ytest+a2*xtest[i]+a3*xtest[i]*ytest-a0) < small  &&
	   fabs(b1*ytest+b2*xtest[i]+b3*xtest[i]*ytest-b0) < small) { */
	  (*xlocal) = xtest[i];
	  (*ylocal) = ytest;
	}
      }
    return(FALSE);
#endif

  }

  if(hit) printf("(x,y)=(%.2g,%.2g)\n",(*xlocal),(*ylocal));
  return(hit);
}

#endif

