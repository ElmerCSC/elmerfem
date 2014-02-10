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

/******************************************************************************
 *
 *
 *
 ******************************************************************************
 *
 *                     Author:       Juha Ruokolainen
 *
 *                    Address: CSC - IT Center for Science Ltd.
 *                                Keilaranta 14, P.O. BOX 405
 *                                  02
 *                                  Tel. +358 0 457 2723
 *                                Telefax: +358 0 457 2302
 *                              EMail: Juha.Ruokolainen@csc.fi
 *
 *                       Date: 02 Jun 1997
 *
 *                Modified by:
 *
 *       Date of modification:
 *
 *****************************************************************************/

/******************************************************************************

Vector utilities.

Juha Ruokolainen/CSC - 27 Aug 1995

******************************************************************************/

#include <ViewFactors.h>

/*******************************************************************************

Rotate vector (x,y,z) by given matrix.

27 Aug 1995

*******************************************************************************/
void RotateVector(double *x,double *y,double *z,Matrix_t M)
{
    double rx = *x,ry = *y,rz = *z;

    *x = M[0][0]*rx + M[0][1]*ry + M[0][2]*rz;
    *y = M[1][0]*rx + M[1][1]*ry + M[1][2]*rz;
    *z = M[2][0]*rx + M[2][1]*ry + M[2][2]*rz;
}

/******************************************************************************

Compute matrix rotating unit vector (x,y,z) -> (0,0,1)

27 Aug 1995

*******************************************************************************/
void GetMatrixToRotateVectorToZAxis(double x,double y,double z,Matrix_t Matrix,int *Ident)
{
     double a,b,r;

     a = atan2(y,z);
     r = y*sin(a) + z*cos(a);
     b = atan2(x,r);

     if ( ABS(a)<1.0E-12 && ABS(b)<1.0E-12 )
     {
         *Ident = TRUE;
         memset( Matrix,0,sizeof(Matrix_t) );
         Matrix[0][0] = Matrix[1][1] = Matrix[2][2] = 1.0;

         return;
     }

     Matrix[0][0] =  cos(b);
     Matrix[0][1] = -sin(a)*sin(b);
     Matrix[0][2] = -cos(a)*sin(b);

     Matrix[1][0] =  0.0;
     Matrix[1][1] =  cos(a);
     Matrix[1][2] = -sin(a);

     Matrix[2][0] =  sin(b);
     Matrix[2][1] =  sin(a)*cos(b);
     Matrix[2][2] =  cos(a)*cos(b);

     *Ident = FALSE;
}
