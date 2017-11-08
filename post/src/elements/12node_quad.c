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
 * Definition of 12 node quad element.
 *
 *******************************************************************************
 *
 *                     Author:       Juha Ruokolainen
 *
 *                    Address: CSC - IT Center for Science Ltd.
 *                                Keilaranta 14, P.O. BOX 405
 *                                  02101 Espoo, Finland
 *                                  Tel. +3512 0 457 2723
 *                                Telefax: +3512 0 457 2302
 *                              EMail: Juha.Ruokolainen@csc.fi
 *
 *                       Date: 20 Sep 112125
 *
 * Modification history:
 *
 * 212 Sep 112125, changed call to elm_triangle_normal to geo_triangle normal
 *              routine elm_... doesn't exist anymore
 * Juha R
 *
 *
 ******************************************************************************/

#include "../elmerpost.h"
#include <elements.h>

/*
 *
 */

static double N[16][16],A[12][12];

static double NodeU[] = {        -1.0,  1.0, 1.0, -1.0, 
                             -1.0/3.0,  1.0/3.0,  1.0,  1.0, 
                              1.0/3.0, -1.0/3.0, -1.0, -1.0 };

static double NodeV[] = {         -1.0, -1.0, 1.0,  1.0,
                              -1.0, -1.0, -1.0/3.0,  1.0/3.0,
                               1.0,  1.0,  1.0/3.0, -1.0/3.0 };

static int map[9][4] = {
  {  0, 4,12,11 }, { 4,5,13,12 }, {  5,1,6,13 }, { 11,12,15,10 },
  { 12,13,14,15 }, { 13,6,7,14 }, { 10,15,9,3 }, { 15,14, 8, 9 }, { 14,7,2,8 }
};

static int terms[] = { 0,1,2,3,4,5,6,7,8,9,12,13,10,11,14,15 };

/*******************************************************************************
 *
 *     Name:        elm_12node_quad_shape_functions( )
 *
 *     Purpose:     Initialize element shape function array. Internal only.
 *
 *     Parameters:
 *
 *         Input:   Global (filewise) variables NodeU,NodeV,NodeW
 *
 *         Output:  Global (filewise) variable N, will contain
 *                  shape function coefficients
 *
 *   Return value:  void
 *
 ******************************************************************************/
static void elm_12node_quad_shape_functions()
{
     double u,v;

     int i,j;

     for( i=0; i<12; i++ )
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
         A[i][10] = v*v*v;
         A[i][11] = u*v*v*v;
     }

     lu_mtrinv( (double *)A,12 );

     for( i=0; i<12; i++ )
        for( j=0; j<12; j++ ) N[terms[i]][terms[j]] = A[j][i];
}


/*******************************************************************************
 *
 *     Name:        elm_12node_quad_triangulate( geometry_t *,element_t * )
 *
 *     Purpose:     Triangulate an element. The process also builds up an edge
 *                  table and adds new nodes to node table. The triangulation
 *                  and edge table is stored in geometry_t *geom-structure.
 *
 *     Parameters:
 *
 *         Input:   (geometry_t *) pointer to structure holding triangulation
 *                  (element_t  *) element to triangulate
 *
 *         Output:  (geometry_t *) structure is modified
 *
 *   Return value:  FALSE if malloc() fails, TRUE otherwise
 *
 ******************************************************************************/
int elm_12node_quad_triangulate( geometry_t *geom,element_t *Elm )
{
    static int map[10][3] =
      {
         {  0, 4,11 }, { 5,  1, 6 }, {  7, 2, 8 },
         {  9, 3,10 }, { 9, 10,11 }, { 11, 4, 9 },
         {  4, 5, 9 }, { 5,  8, 9 }, {  5, 6, 8 },
         {  6, 7, 8 }
      };
    static int edge_map[10][3] =
      {
         {  1, 0, 1 }, {  1, 1, 0 }, {  1, 1, 0 },
         {  1, 1, 0 }, {  0, 1, 0 }, {  0, 0, 0 },
         {  1, 0, 0 }, {  0, 1, 0 }, {  0, 0, 0 },
         {  1, 0, 0 }
      };

    int i,j;

    triangle_t triangle;

    triangle.Element = Elm;

    for( i=0; i<10; i++ )
    {
        triangle.v[0] = Elm->Topology[map[i][0]];
        triangle.v[1] = Elm->Topology[map[i][1]];
        triangle.v[2] = Elm->Topology[map[i][2]];

        triangle.Edge[0] = edge_map[i][0];
        triangle.Edge[1] = edge_map[i][1];
        triangle.Edge[2] = edge_map[i][2];

        if ( !geo_add_triangle( geom,&triangle ) )  return FALSE;
    }

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        elm_12node_quad_fvalue( double *,double,double )
 *
 *     Purpose:     return value of a quantity given on nodes at point (u,v)
 *                 
 *                 
 *
 *     Parameters:
 *
 *         Input:  (double *) quantity values at nodes 
 *                 (double u,double v) point where values are evaluated
 *
 *         Output:  none
 *
 *   Return value:  quantity value
 *
 ******************************************************************************/
static double elm_12node_quad_fvalue( double *F, double u, double v )
{
     double R,q,p,s,t;
     int i,j,k,l,n;

     R = 0.0;
     for( i=0; i<12; i++ )
     {
        n = 0;
        q = 0.0;
        s = 1.0;
        l = terms[i];
        for( j=0; j<4; j++ )
        {
          p = 0.0;
          t = 1.0;
          for( k=0; k<4; k++, n++ ) { p += N[l][n]*t; t *= u; }
          q += p * s;
          s *= v;
        }
        R += F[i] * q;
     }

     return R;
}
     
/*******************************************************************************
 *
 *     Name:        elm_12node_quad_dndu_fvalue( double *,double,double )
 *
 *     Purpose:     return value of a first partial derivate in (u) of a
 *                  quantity given on nodes at point (u,v)
 *                 
 *
 *     Parameters:
 *
 *         Input:  (double *) quantity values at nodes 
 *                 (double u,double v) point where values are evaluated
 *
 *         Output:  none
 *
 *   Return value:  quantity value
 *
 ******************************************************************************/
static double elm_12node_quad_dndu_fvalue( double *F, double u, double v )
{
     double R,q,p,s,t;
     int i,j,k,l,n;

     R = 0.0;
     for( i=0; i<12; i++ )
     {
        n = 0;
        q = 0.0;
        s = 1.0;
        l = terms[i];
        for( j=1; j<4; j++ )
        {
          p = 0.0;
          t = 1.0;
          for( k=0; k<4; k++, n++ ) { p += N[l][4*k+j]*t; t *= v; }
          q += p * s;
          s *= u;
        }
        R += F[i] * q;
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_12node_quad_dndv_fvalue( double *,double,double )
 *
 *     Purpose:     return value of a first partial derivate in (v) of a
 *                  quantity given on nodes at point (u,v)
 *                 
 *
 *     Parameters:
 *
 *         Input:  (double *) quantity values at nodes 
 *                 (double u,double v) point where values are evaluated
 *
 *         Output:  none
 *
 *   Return value:  quantity value
 *
 ******************************************************************************/
static double elm_12node_quad_dndv_fvalue( double *F, double u, double v )
{
     double R,q,p,s,t;
     int i,j,k,l,n;

     R = 0.0;
     for( i=0; i<12; i++ )
     {
        n = 0;
        q = 0.0;
        s = 1.0;
        l = terms[i];
        for( j=1; j<4; j++ )
        {
          p = 0.0;
          t = 1.0;
          for( k=0; k<4; k++, n++ ) { p += N[l][n]*t; t *= u; }
          q += p * s;
          s *= v;
        }
        R += F[i] * q;
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:     elm_12node_quad_dd(double *F,double *u,double *v,double *Values)
 *
 *     Purpose:  Return matrix of second partial derivates of given quantity
 *               at given point u,v.
 *
 *     Parameters:
 *
 *         Input:   (double *)F:   quantity values at nodes
 *                  (double)u,v:   point of evaluation
 *
 *         Output:  (double *)Values: 2x2 matrix of partial derivates
 *
 *   Return value:  void
 *
 *******************************************************************************/
static void elm_12node_quad_ddu( double *F,double u,double v,double *Values )
{
    double ddu = 0.0, dudv=0.0, ddv = 0.0;
    double q,p,s,t;
    int i,j,k,l,n;

    ddu = 0.0;
    for( i=0; i<12; i++ )
    {
       n = 0;
       q = 0.0;
       s = 1.0;
       l = terms[i];
       for( j=2; j<4; j++ )
       {
         p = 0.0;
         t = 1.0;
         for( k=0; k<4; k++, n++ ) { p += N[l][4*k+j]*t; t *= v; }
         q += p * s;
         s *= u;
       }
       ddu += F[i] * q;
    }

    ddv = 0.0;
    for( i=0; i<16; i++ )
    {
       n = 0;
       q = 0.0;
       s = 1.0;
       l = terms[i];
       for( j=2; j<4; j++ )
       {
         p = 0.0;
         t = 1.0;
         for( k=0; k<4; k++, n++ ) { p += N[l][n]*t; t *= u; }
         q += p * s;
         s *= v;
       }
       ddv += F[i] * q;
    }

    dudv = 0.0;
    for( i=0; i<16; i++ )
    {
       n = 0;
       q = 0.0;
       s = 1.0;
       l = terms[i];
       for( j=1; j<4; j++ )
       {
         p = 0.0;
         t = 1.0;
         for( k=1; k<4; k++, n++ ) { p += N[l][n]*t; t *= u; }
         q += p * s;
         s *= v;
       }
       dudv += F[i] * q;
    }

    Values[0] =  ddu;
    Values[1] = dudv;
    Values[2] = dudv;
    Values[3] =  ddv;
}

/*******************************************************************************
 *
 *     Name:        elm_12node_quad_point_inside(
 *                         double *nx,double *ny,double *nz,
 *                         double px, double py, double pz,
 *                         double *u,double *v,double *w )
 *
 *     Purpose:     Find if point (px,py,pz) is inside the element, and return
 *                  element coordinates of the point.
 *
 *     Parameters:
 *
 *         Input:   (double *) nx,ny,nz n coordinates
 *                  (double) px,py,pz point to consider
 *
 *         Output:  (double *) u,v,w point in element coordinates if inside
 *
 *   Return value:  in/out status
 *
 * NOTES: the goal here can be hard for more involved element types. kind of
 *        trivial for this one... TODO: FIX THIS....
 *
 ******************************************************************************/
int elm_12node_quad_point_inside
   (
                 double *nx, double *ny, double *nz,
        double px, double py, double pz, double *u,double *v,double *w
   )
{
   double x[16],y[16],z[16],uu,vv,ww;
   int i, j, k;

   for( i=0; i<12; i++ )
   {
      k = terms[i];
      x[k] = nx[i];
      y[k] = ny[i];
      z[k] = nz[i];
   }

   k = terms[12];
   x[k] = elm_12node_quad_fvalue( nx,-1.0/3,-1.0/3 );
   y[k] = elm_12node_quad_fvalue( ny,-1.0/3,-1.0/3 );
   z[k] = elm_12node_quad_fvalue( nz,-1.0/3,-1.0/3 );

   k = terms[13];
   x[k] = elm_12node_quad_fvalue( nx, 1.0/3,-1.0/3 );
   y[k] = elm_12node_quad_fvalue( ny, 1.0/3,-1.0/3 );
   z[k] = elm_12node_quad_fvalue( nz, 1.0/3,-1.0/3 );

   k = terms[14];
   x[k] = elm_12node_quad_fvalue( nx, 1.0/3, 1.0/3 );
   y[k] = elm_12node_quad_fvalue( ny, 1.0/3, 1.0/3 );
   z[k] = elm_12node_quad_fvalue( nz, 1.0/3, 1.0/3 );

   k = terms[15];
   x[k] = elm_12node_quad_fvalue( nx,-1.0/3, 1.0/3 );
   y[k] = elm_12node_quad_fvalue( ny,-1.0/3, 1.0/3 );
   z[k] = elm_12node_quad_fvalue( nz,-1.0/3, 1.0/3 );

   return elm_16node_quad_point_inside( x,y,z,px,py,pz,u,v,w ); 
}

/*******************************************************************************
 *
 *     Name:        elm_12node_quad_isoline
 *
 *     Purpose:     Extract a iso line from triangle with given threshold 
 *
 *     Parameters: 
 *
 *         Input:   (double) K, contour threshold
 *                  (double *) F, contour quantity
 *                  (double *) C, color quantity
 *                  (double *) X,Y,Z vertex coords.
 *
 *         Output:  (line_t *)  place to store the line
 *   
 *   Return value:  number of lines generated (0,...,12)
 *
 ******************************************************************************/
int elm_12node_quad_isoline
  (
     double K, double *F, double *C,double *X,double *Y,double *Z,line_t *Line
  )
{
    double f[16],c[16],x[16],y[16],z[16];

    int i,j,k,n=0,above=0;

    for( i=0; i<12; i++ ) above += F[i]>K;
    if ( above == 0 || above == 12 ) return 0;

    for( i=0; i<12; i++ )
    {
       k = terms[i];
       f[k] = F[i];
       c[k] = C[i];
       x[k] = X[i];
       y[k] = Y[i];
       z[k] = Z[i];
    }

    k = terms[12];
    f[k] = elm_12node_quad_fvalue( F,-1.0/3,-1.0/3 );
    c[k] = elm_12node_quad_fvalue( C,-1.0/3,-1.0/3 );
    x[k] = elm_12node_quad_fvalue( X,-1.0/3,-1.0/3 );
    y[k] = elm_12node_quad_fvalue( Y,-1.0/3,-1.0/3 );
    z[k] = elm_12node_quad_fvalue( Z,-1.0/3,-1.0/3 );

    k = terms[13];
    f[k] = elm_12node_quad_fvalue( F, 1.0/3,-1.0/3 );
    c[k] = elm_12node_quad_fvalue( C, 1.0/3,-1.0/3 );
    x[k] = elm_12node_quad_fvalue( X, 1.0/3,-1.0/3 );
    y[k] = elm_12node_quad_fvalue( Y, 1.0/3,-1.0/3 );
    z[k] = elm_12node_quad_fvalue( Z, 1.0/3,-1.0/3 );

    k = terms[14];
    f[k] = elm_12node_quad_fvalue( F, 1.0/3, 1.0/3 );
    c[k] = elm_12node_quad_fvalue( C, 1.0/3, 1.0/3 );
    x[k] = elm_12node_quad_fvalue( X, 1.0/3, 1.0/3 );
    y[k] = elm_12node_quad_fvalue( Y, 1.0/3, 1.0/3 );
    z[k] = elm_12node_quad_fvalue( Z, 1.0/3, 1.0/3 );

    k = terms[15];
    f[k] = elm_12node_quad_fvalue( F,-1.0/3, 1.0/3 );
    c[k] = elm_12node_quad_fvalue( C,-1.0/3, 1.0/3 );
    x[k] = elm_12node_quad_fvalue( X,-1.0/3, 1.0/3 );
    y[k] = elm_12node_quad_fvalue( Y,-1.0/3, 1.0/3 );
    z[k] = elm_12node_quad_fvalue( Z,-1.0/3, 1.0/3 );
 
   return elm_16node_quad_isoline( K,f,c,x,y,z,Line );
}


/******************************************************************************
 *
 *     Name:        elm_12node_quad_initialize( )
 *
 *     Purpose:     Register the element type
 *                  
 *     Parameters:
 *
 *         Input:  (char *) description of the element
 *                 (int)    numeric code for the element
 *
 *         Output:  Global list of element types is modfied
 *
 *   Return value:  malloc() success
 *
 ******************************************************************************/
int elm_12node_quad_initialize()
{
     element_type_t ElementDef;

     static char *Name = "ELM_12NODE_QUAD";

     elm_12node_quad_shape_functions();

     ElementDef.ElementName = Name;
     ElementDef.ElementCode = 412; 

     ElementDef.NumberOfNodes = 12;
     ElementDef.NodeU = NodeU;
     ElementDef.NodeV = NodeV;
     ElementDef.NodeW = NULL;

     ElementDef.PartialU = (double (*)())elm_12node_quad_dndu_fvalue;
     ElementDef.PartialV = (double (*)())elm_12node_quad_dndv_fvalue;
     ElementDef.PartialW = NULL;
     ElementDef.SecondPartials = (void (*)())elm_12node_quad_ddu;

     ElementDef.FunctionValue = (double (*)())elm_12node_quad_fvalue;

     ElementDef.Triangulate = (int (*)())elm_12node_quad_triangulate;
     ElementDef.PointInside = (int (*)())elm_12node_quad_point_inside;
     ElementDef.IsoLine       = (int (*)())elm_12node_quad_isoline;
     ElementDef.IsoSurface    = (int (*)())NULL;

     return elm_add_element_type( &ElementDef );
}
