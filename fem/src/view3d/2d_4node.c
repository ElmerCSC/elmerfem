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
 *
 * shape functions: N0 = (1-u)(1-v)/4       3--------2
 *                  N1 = (1+u)(1-v)/4       |        |
 *                  N2 = (1+u)(1+v)/4     v |        | 
 *                  N3 = (1-u)(1+v)/4       0--------1
 *                                            u
 */
#define AEPS 1.0e-18

/*
extern double NodeArray[][3];
*/

static double N[4][4] = {
     { 1.0/1.0, -1.0/1.0, -1.0/1.0,  1.0/1.0 },  /* 1 - u - v + uv */
     { 1.0/1.0,  1.0/1.0, -1.0/1.0, -1.0/1.0 },
     { 1.0/1.0,  1.0/1.0,  1.0/1.0,  1.0/1.0 },
     { 1.0/1.0, -1.0/1.0,  1.0/1.0, -1.0/1.0 }
  }, A[4][4];

static double NodeU[] = {  0.0,  1.0, 1.0,  0.0 };
static double NodeV[] = {  0.0,  0.0, 1.0,  1.0 };

void elm_4node_quad_shape_functions(double *B)
{
     double u,v;
     int i,j;

     for( i=0; i<4; i++ )
     {
         u = NodeU[i];
         v = NodeV[i];

         A[i][0]  = 1;
         A[i][1]  = u;
         A[i][2]  = v;
         A[i][3]  = u*v;
     }

     mtrinv( (double *)A,4 );

     for( i=0; i<4; i++ )
     {
         for( j=0; j<4; j++ ) { N[i][j] = A[j][i]; *B++ = N[i][j]; }
     }
}

void Nvalue_4node( double *R, double u, double v )
{
     int i;

     for( i=0; i<4; i++ )
     {
          R[i] = N[i][0] + N[i][1]*u + N[i][2]*v + N[i][3]*u*v;
     }
}

double Fvalue_4node( double *F, double u, double v )
{
     double R;
     int i;

     R = 0;
     for( i=0; i<4; i++ )
     {
          R += F[i]*(N[i][0] + N[i][1]*u + N[i][2]*v + N[i][3]*u*v);
     }

     return R;
}

void dNdU_Nvalue_4node( double *R, double u, double v )
{
     int i;

     for( i=0; i<4; i++ )
     {
         R[i] = N[i][1] + N[i][3]*v;
     }
}

void dNdV_Nvalue_4node( double *R, double u, double v )
{
     int i;

     for( i=0; i<4; i++ )
     {
         R[i] = N[i][2] + N[i][3]*u;
     }
}

double dNdU_Fvalue_2node( double *F, double u)
{
     return 0.5*(F[1] - F[0]);
}
     
double dNdU_Fvalue_4node( double *F, double u, double v )
{
     double R;
     int i;

     R = 0.0;
     for( i=0; i<4; i++ )
     {
         R += F[i]*(N[i][1] + N[i][3]*v);
     }

     return R;
}

double dNdV_Fvalue_4node( double *F, double u, double v )
{
     double R;

     int i;

     R = 0.0;
     for( i=0; i<4; i++ )
     {
         R += F[i]*(N[i][2] + N[i][3]*u);
     }

     return R;
}

void ddNddU_Nvalue_4node(double *R,double u, double v)
{
    int i;

    for( i=0; i<4; i++ ) R[i] = 0.0;
}

void ddNddV_Nvalue_4node(double *R,double u, double v)
{
    int i;

    for( i=0; i<4; i++ ) R[i] = 0.0;
}

void ddNdUdV_Nvalue_4node(double *R,double u, double v)
{
    int i;

    for( i=0; i<4; i++ ) R[i] = N[i][3];
}

double ddNddU_Fvalue_node(double *F,double u, double v)
{
    return 0.0;
}

double ddNddV_Fvalue_node(double *F,double u, double v)
{
    return 0.0;
}

double ddNdUdV_Fvalue_4node(double *F,double u,double v)
{
    double G=0.0;
    int i;

    for( i=0; i<4; i++ ) G += F[i]*N[i][3];

    return G;
}

double ElementOfLine_4node(double *X,double *Y,double *Z,double u,double v)
{
   double dXdU,dXdV, dYdU,dYdV, dZdU,dZdV;
   double detA, Auu, Auv, Avv;
   int i;

    dXdU = dNdU_Fvalue_4node(X,u,v);
    dYdU = dNdU_Fvalue_4node(Y,u,v);
    dZdU = dNdU_Fvalue_4node(Z,u,v);

    dXdV = dNdV_Fvalue_4node(X,u,v);
    dYdV = dNdV_Fvalue_4node(Y,u,v);
    dZdV = dNdV_Fvalue_4node(Z,u,v);

    Auu = dXdU*dXdU + dYdU*dYdU + dZdU*dZdU; /* surface metric a    */
    Auv = dXdU*dXdV + dYdU*dYdV + dZdU*dZdV; /*                 ij  */
    Avv = dXdV*dXdV + dYdV*dYdV + dZdV*dZdV;

    detA = Auu*Avv - Auv*Auv;

    if ( detA < AEPS )
    {
        for( i=0; i<4; i++ ) fprintf( stderr, " %8g %8g\n", X[i],Y[i] );

        fprintf( stderr, "area: Element Jacobian singular at: %g %g\n",u,v);
        return 0.0;
    }

    return sqrt(Auu);
}

double ElementOfArea_4node(double *X,double *Y,double *Z,double u,double v)
{
    double detA, Auu, Auv, Avv;
    double dXdU,dXdV, dYdU,dYdV, dZdU,dZdV;

    int i;

    dXdU = dNdU_Fvalue_4node(X,u,v);
    dYdU = dNdU_Fvalue_4node(Y,u,v);
    dZdU = dNdU_Fvalue_4node(Z,u,v);

    dXdV = dNdV_Fvalue_4node(X,u,v);
    dYdV = dNdV_Fvalue_4node(Y,u,v);
    dZdV = dNdV_Fvalue_4node(Z,u,v);

    Auu = dXdU*dXdU + dYdU*dYdU + dZdU*dZdU; /* surface metric a    */
    Auv = dXdU*dXdV + dYdU*dYdV + dZdU*dZdV; /*                 ij  */
    Avv = dXdV*dXdV + dYdV*dYdV + dZdV*dZdV;

    detA = Auu*Avv - Auv*Auv;

    if ( detA < AEPS )
    {
        for( i=0; i<4; i++ ) fprintf( stderr, " %8g %8g\n", X[i],Y[i] );

        fprintf( stderr, "area: Element Jacobian singular at: %g %g\n",u,v);
        return 0.0;
    }

    return sqrt(detA);
}

void derivates_to_global_4node( double *X,double *Y,double *Z,
                                double *dFdX,  /* in @F/@u, out @F/@x */
                                double *dFdY,  /* in @F/@v, out @F/@y */
                                double *dFdZ,  /* out @F/@z */
                                double *u, double *v, int N )
{
    double detA,Auu,Auv,Avv;
    double dXdU,dXdV,dYdU,dYdV,dZdU,dZdV,a,b;
    int i;

    for( i=0; i<N; i++ )
    {
        dXdU = dNdU_Fvalue_4node(X,u[i],v[i]);
        dYdU = dNdU_Fvalue_4node(Y,u[i],v[i]);
        dZdU = dNdU_Fvalue_4node(Z,u[i],v[i]);

        dXdV = dNdV_Fvalue_4node(X,u[i],v[i]);
        dYdV = dNdV_Fvalue_4node(Y,u[i],v[i]);
        dZdV = dNdV_Fvalue_4node(Z,u[i],v[i]);

        Avv = dXdU*dXdU + dYdU*dYdU + dZdU*dZdU; /* surface metric a    */
        Auv = dXdU*dXdV + dYdU*dYdV + dZdU*dZdV; /*                 ij  */
        Auu = dXdV*dXdV + dYdV*dYdV + dZdV*dZdV;

        detA = Auu*Avv - Auv*Auv;
        if (detA < AEPS)
        {
            for( i=0; i<4; i++ ) fprintf( stderr, " %8g %8g\n", X[i],Y[i] );

            fprintf( stderr, "deriv: Element Jacobian singular at: %g %g\n",u[i],v[i]);
            return;
        }

        Auv = -Auv / detA;                       /*   ij  */
        Auu =  Auu / detA;                       /*  a    */
        Avv =  Avv / detA;

        a = dFdX[i];
        b = dFdY[i];
        dFdX[i] = Auu*a + Auv*b;                 /* raise index of the surface */
        dFdY[i] = Auv*a + Avv*b;                 /*   vector (@f/@u,@f/@v)     */

        /*
         * transform to global cartesian frame
         */
        a = dFdX[i];
        b = dFdY[i];
        dFdX[i] = dXdU*a + dXdV*b;
        dFdY[i] = dYdU*a + dYdV*b;
        dFdZ[i] = dZdU*a + dZdV*b;
    } 
}

void dNdXYZ_Nvalue_4node( int *Topology, double *dFdX, double *dFdY, double *dFdZ, double u, double v )
{
    double X[4],Y[4];
    int i;

    for( i=0; i<4; i++ ) { X[i] = u; Y[i] = v; }

    dNdU_Nvalue_4node(dFdX,u,v);
    dNdV_Nvalue_4node(dFdY,u,v);

/*
    derivates_to_global_4node(Topology,dFdX,dFdY,dFdZ,X,Y,4);
*/
}
