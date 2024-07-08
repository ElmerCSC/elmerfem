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

#include <stdio.h>
#include <math.h>


void mtrinv(double *,int);

/*
 * sixteen node (cubic) quad surface element
 * 
 *      3---9---8---2
 *      |           |
 *     10  15  14   7
 *      |           |
 *     11  12  13   6 
 *    v |           | 
 *      0---4---5---1
 *       u
 */

static double A[16][16],N[16][16];

static double NodeU[] =
{
       0.0,     1.0,  1.0,  0.0, 1.0/3.0, 2.0/3.0,  1.0,     1.0,
     2.0/3.0, 1.0/3.0, 0.0, 0.0, 1.0/3.0, 2.0/3.0, 2.0/3.0, 1.0/3.0
 };

static double NodeV[] =
{
     0.0, 0.0,   1.0,     1.0,     0.0,    0.0,   1.0/3.0, 2.0/3.0,
     1.0, 1.0, 2.0/3.0, 1.0/3.0, 1.0/3.0,1.0/3.0, 2.0/3.0, 2.0/3.0
};

void elm_16node_quad_add();

void elm_16node_quad_shape_functions(double *B)
{
     double u,v;
     int i,j;

     for( i=0; i<16; i++ )
     {
         u = NodeU[i];
         v = NodeV[i];

         A[i][0]  = 1;
         A[i][1]  = u;
         A[i][2]  = u*u;
         A[i][3]  = u*u*u;
         A[i][4]  = v;
         A[i][5]  = u*v;
         A[i][6]  = u*u*v;
         A[i][7]  = u*u*u*v;
         A[i][8]  = v*v;
         A[i][9]  = u*v*v;
         A[i][10] = u*u*v*v;
         A[i][11] = u*u*u*v*v;
         A[i][12] = v*v*v;
         A[i][13] = u*v*v*v;
         A[i][14] = u*u*v*v*v;
         A[i][15] = u*u*u*v*v*v;
     }

     mtrinv( (double *)A,16 );

     for( i=0; i<16; i++ )
     {
         for( j=0; j<16; j++ ) { N[i][j] = A[j][i]; *B++ = N[i][j]; }
     }
}

void elm_16node_quad_fit(double *X,double *Y,double *Z,double *R)
{
     double u,v;
     int i,j;

     for( i=0; i<16; i++ )
     {
         u = X[i];
         v = Y[i];

         A[i][0]   = 1;
         A[i][1]   = u;
         A[i][2]   = v;
         A[i][3]   = u*u;
         A[i][4]   = u*v;
         A[i][5]   = v*v;
         A[i][6]   = u*u*u;
         A[i][7]   = u*u*v;
         A[i][8]   = u*v*v;
         A[i][9]   = v*v*v;
         A[i][10]  = u*u*u*v;
         A[i][11]  = u*u*v*v;
         A[i][12]  = u*v*v*v;
         A[i][13]  = u*u*u*v*v;
         A[i][14]  = u*u*v*v*v;
         A[i][15]  = u*u*u*v*v*v;
     }

     mtrinv( (double *)A,16 );

     for( i=0; i<16; i++ )
     {
         R[i] = 0.0;
         for( j=0; j<16; j++ ) R[i] += A[i][j]*Z[j];
     }
}

double elm_16node_quad_fvalue(double *F,double u,double v)
{
     double R = 0.0,uv=u*v,uu=u*u,vv=v*v,uuu=u*u*u,vvv=v*v*v,uuv=u*u*v,
                    uvv=u*v*v,uuuv=uuu*v,uvvv=u*vvv,uuvv=uu*vv,
                    uuuvv=uuu*vv,uuvvv=uu*vvv,uuuvvv=uuu*vvv;
     int i;

     for( i=0; i<16; i++ )
     {
         R += F[i]*(
                      N[i][0] +
                      N[i][1]*u +
                      N[i][2]*uu +
                      N[i][3]*uuu + 
                      N[i][4]*v +
                      N[i][5]*uv +
                      N[i][6]*uuv +
                      N[i][7]*uuuv +
                      N[i][8]*vv + 
                      N[i][9]*uvv +
                      N[i][10]*uuvv +
                      N[i][11]*uuuvv +
                      N[i][12]*vvv +
                      N[i][13]*uvvv +
                      N[i][14]*uuvvv +
                      N[i][15]*uuuvvv
                   );
     } 

     return R;
}

void elm_16node_quad_nvalue(double *F,double u,double v)
{
     double R = 0.0,uv=u*v,uu=u*u,vv=v*v,uuu=u*u*u,vvv=v*v*v,uuv=u*u*v,
                    uvv=u*v*v,uuuv=uuu*v,uvvv=u*vvv,uuvv=uu*vv,
                    uuuvv=uuu*vv,uuvvv=uu*vvv,uuuvvv=uuu*vvv;
     int i;

     for( i=0; i<16; i++ )
     {
        F[i] = (
                      N[i][0] +
                      N[i][1]*u +
                      N[i][2]*v +
                      N[i][3]*uu +
                      N[i][4]*uv +
                      N[i][5]*vv + 
                      N[i][6]*uuu + 
                      N[i][7]*uuv +
                      N[i][8]*uvv +
                      N[i][9]*vvv +
                      N[i][10]*uuuv +
                      N[i][11]*uuvv +
                      N[i][12]*uvvv +
                      N[i][13]*uuuvv +
                      N[i][14]*uuvvv +
                      N[i][15]*uuuvvv
                   );
     } 
}

double elm_16node_quad_dndu_fvalue(double *F,double u,double v)
{
     double R=0.0,u2=2*u,uu3=3*u*u,vv=v*v,u2v=2*u*v,uu3v=uu3*v,u2vv=u2*vv,
                  vvv=vv*v,uu3vv=uu3*vv,u2vvv=u2*vvv,uu3vvv=uu3*vvv;
     int i;

     for( i=0; i<16; i++ )
     {
         R += F[i]*(
                      N[i][1] +
                      N[i][3]*u2 +
                      N[i][4]*v +
                      N[i][6]*uu3 + 
                      N[i][7]*u2v +
                      N[i][8]*vv +
                      N[i][10]*uu3v +
                      N[i][11]*u2vv +
                      N[i][12]*vvv +
                      N[i][13]*uu3vv +
                      N[i][14]*u2vvv +
                      N[i][15]*uu3vvv
                   );
     }

     return R;
}

void elm_16node_quad_dndu_nvalue(double *F,double u,double v)
{
     double u2=2*u,uu3=3*u*u,vv=v*v,u2v=2*u*v,uu3v=uu3*v,u2vv=u2*vv,
              vvv=vv*v,uu3vv=uu3*vv,u2vvv=u2*vvv,uu3vvv=uu3*vvv;
     int i;

     for( i=0; i<16; i++ )
     {
         F[i] = (
                     N[i][1] +
                     N[i][3]*u2 +
                     N[i][4]*v +
                     N[i][6]*uu3 + 
                     N[i][7]*u2v +
                     N[i][8]*vv +
                     N[i][10]*uu3v +
                     N[i][11]*u2vv +
                     N[i][12]*vvv +
                     N[i][13]*uu3vv +
                     N[i][14]*u2vvv +
                     N[i][15]*uu3vvv
                 );
     }
}      

double elm_16node_quad_dndv_fvalue(double *F,double u,double v)
{
     int i;

     double R=0.0,v2=2*v,uu=u*u,uv2=u*v2,vv3=3*v*v,uuu=uu*u,uuv2=uu*v2,
                  uvv3=u*vv3,uuuv2=uuu*v2,uuvv3=uu*vv3,uuuvv3=uuu*vv3;

     for( i=0; i<16; i++ )
     {
         R += F[i]*(
                      N[i][2] +
                      N[i][4]*u +
                      N[i][5]*v2 + 
                      N[i][7]*uu +
                      N[i][8]*uv2 +
                      N[i][9]*vv3 +
                      N[i][10]*uuu +
                      N[i][11]*uuv2 +
                      N[i][12]*uvv3 +
                      N[i][13]*uuuv2 +
                      N[i][14]*uuvv3 +
                      N[i][15]*uuuvv3
                   );
     } 
     return R;
}

void elm_16node_quad_dndv_nvalue(double *F,double u,double v)
{
     int i;

     double v2=2*v,uu=u*u,uv2=u*v2,vv3=3*v*v,uuu=uu*u,uuv2=uu*v2,
                  uvv3=u*vv3,uuuv2=uuu*v2,uuvv3=uu*vv3,uuuvv3=uuu*vv3;

     for( i=0; i<16; i++ )
     {
         F[i] = (
                      N[i][2] +
                      N[i][4]*u +
                      N[i][5]*v2 + 
                      N[i][7]*uu +
                      N[i][8]*uv2 +
                      N[i][9]*vv3 +
                      N[i][10]*uuu +
                      N[i][11]*uuv2 +
                      N[i][12]*uvv3 +
                      N[i][13]*uuuv2 +
                      N[i][14]*uuvv3 +
                      N[i][15]*uuuvv3
                   );
     } 
}

double ElementOfArea_16node(double *X,double *Y,double *Z,double u,double v)
{
    double dXdU,dXdV, dYdU,dYdV, dZdU,dZdV;
    double detA, Auu, Auv, Avv;
    int i;


    dXdU = elm_16node_quad_dndu_fvalue(X,u,v);
    dYdU = elm_16node_quad_dndu_fvalue(Y,u,v);
    dZdU = elm_16node_quad_dndu_fvalue(Z,u,v);

    dXdV = elm_16node_quad_dndv_fvalue(X,u,v);
    dYdV = elm_16node_quad_dndv_fvalue(Y,u,v);
    dZdV = elm_16node_quad_dndv_fvalue(Z,u,v);

    Auu = dXdU*dXdU + dYdU*dYdU + dZdU*dZdU; /* surface metric a    */
    Auv = dXdU*dXdV + dYdU*dYdV + dZdU*dZdV; /*                 ij  */
    Avv = dXdV*dXdV + dYdV*dYdV + dZdV*dZdV;

    detA = Auu*Avv - Auv*Auv;

    if ( detA < 1.0E-12 )
    {
        for( i=0; i<16; i++ ) fprintf( stderr, " %8g %8g\n", X[i],Y[i] );

        fprintf( stderr, "area: Element Jacobian singular at: %g %g\n",u,v);
        return 0.0;
    }

    return sqrt(detA);
}
