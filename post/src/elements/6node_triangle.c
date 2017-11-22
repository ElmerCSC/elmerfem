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
 * Definition of six node triangle element.
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
 *                       Date: 20 Sep 1995
 *
 *
 * Modification history:
 *
 ******************************************************************************/

#include "../elmerpost.h"
#include <elements.h>


/*
 * SIX NODE TRIANGLE ELEMENT
 * 
 * shape functions: N0 = (1-u-v)(1-2u-2v)              2  v=1
 *                  N1 = u(2u-1)                     /   \
 *                  N2 = v(2v-1)                    /     \
 *                  N3 = 4u(1-u-v)                 5       4
 *                  N4 = 4uv                      /         \
 *                  N5 = 4v(1-u-v)               0-----3-----1 u=1
 *
 */

static double N[6][6],A[6][6];

static double NodeU[] = { 0.0, 1.0, 0.0, 0.5, 0.5, 0.0 };
static double NodeV[] = { 0.0, 0.0, 1.0, 0.0, 0.5, 0.5 };

/*******************************************************************************
 *
 *     Name:        elm_6node_triangle_shape_functions( )
 *
 *     Purpose:     Initialize element shape function array. Internal only.
 *
 *     Parameters:
 *
 *         Input:   Global (filewise) variables NodeU,NodeV,NodeW
 *
 *         Output:  Global (filewise) variable N[6][6], will contain
 *                  shape function coefficients
 *
 *   Return value:  void
 *
 ******************************************************************************/
static void elm_6node_triangle_shape_functions()
{
     double u,v;

     int i,j;

     for( i=0; i<6; i++ )
     {
         u = NodeU[i];
         v = NodeV[i];

         A[i][0]  = 1;
         A[i][1]  = u;
         A[i][2]  = v;
         A[i][3]  = u*v;
         A[i][4]  = u*u;
         A[i][5]  = v*v;
     }

     lu_mtrinv( (double *)A,6 );

     for( i=0; i<6; i++ )
        for( j=0; j<6; j++ ) N[i][j] = A[j][i];
}

/*******************************************************************************
 *
 *     Name:        elm_6node_triangle_triangulate( geometry_t *,element_t * )
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
int elm_6node_triangle_triangulate( geometry_t *geom, element_t *Elm, element_t *Parent )
{
    double u,v,w,s;

    static int map[4][3] =
      {
         { 0,3,5 }, { 3,1,4 }, { 3,4,5 }, { 4,2,5 }
      };
    static int edge_map[4][3] =
      {
         { 1,0,1 }, { 1,1,0 }, { 0,0,0 }, { 1,1,0 }
      };

    int i,j;

    triangle_t triangle;

    triangle.Element = Parent;
    for( i=0; i<4; i++ )
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
 *     Name:        elm_6node_triangle_point_inside(
 *                         double *nx,double *ny,double *nz,
 *                         double px, double py, double pz,
 *                         double *u,double *v,double *w )
 *
 *     Purpose:     Find if point (px,py,pz) is inside the element, and return
 *                  element coordinates of the point.
 *
 *     Parameters:
 *
 *         Input:   (double *) nx,ny,nz node coordinates
 *                  (double) px,py,pz point to consider
 *
 *         Output:  (double *) u,v,w point in element coordinates if inside
 *
 *   Return value:  in/out status
 *
 * NOTES: the goal here can be hard for more involved element types. kind of
 *        trivial for this one... (TODO: this is for xy-plane tri only)
 *
 ******************************************************************************/
int elm_6node_triangle_point_inside
   ( 
                   double *nx, double *ny, double *nz,
         double px, double py, double pz, double *u,double *v,double *w
   )
{
    return elm_3node_triangle_point_inside( nx,ny,nz,px,py,pz,u,v,w );
}

/*******************************************************************************
 *
 *     Name:        elm_6node_triangle_fvalue( double *,double,double )
 *
 *     Purpose:     return value of a quantity given on nodes at point (u,v)
 *                 
 *
 *     Parameters:
 *
 *         Input:  (double *) quantity values at nodes 
 *                 (double u,double v) point where value is evaluated
 *
 *         Output:  none
 *
 *   Return value:  quantity value
 *
 ******************************************************************************/
static double elm_6node_triangle_fvalue( double *F,double u,double v )
{
     double R=0.0,uv=u*v,u2=u*u,v2=v*v;
     int i;

     for( i=0; i<6; i++ )
     {
         R += F[i]*( N[i][0]    +
                     N[i][1]*u  +
                     N[i][2]*v  +
                     N[i][3]*uv +
                     N[i][4]*u2 +
                     N[i][5]*v2 );
     } 

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_6node_triangle_dndu_fvalue( double *,double,double )
 *
 *     Purpose:     return value of a first partial derivate in (u) of a
 *                  quantity given on nodes at point (u,v)
 *                 
 *
 *     Parameters:
 *
 *         Input:  (double *) quantity values at nodes 
 *                 (double u,double v) point where value is evaluated
 *
 *         Output:  none
 *
 *   Return value:  quantity value
 *
 ******************************************************************************/
static double elm_6node_triangle_dndu_fvalue(double *F,double u,double v)
{
     double R=0.0,u2=2*u;
     int i;

     for( i=0; i<6; i++ )
     {
         R += F[i]*( N[i][1] + N[i][3]*v + N[i][4]*u2 );
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_6node_triangle_dndv_fvalue( double *,double,double )
 *
 *     Purpose:     return value of a first partial derivate in (v) of a
 *                  quantity given on nodes at point (u,v)
 *                 
 *
 *     Parameters:
 *
 *         Input:  (double *) quantity values at nodes 
 *                 (double u,double v) point where value is evaluated
 *
 *         Output:  none
 *
 *   Return value:  quantity value
 *
 ******************************************************************************/
static double elm_6node_triangle_dndv_fvalue(double *F,double u,double v)
{
     double R=0.0,v2=2*v;
     int i;

     for( i=0; i<6; i++ )
     {
         R += F[i]*(N[i][2]+N[i][3]*u+N[i][5]*v2);
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_6node_triangle_isoline
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
 *   Return value:  number of lines generated (0,...,6)
 *
 ******************************************************************************/
int elm_6node_triangle_isoline
  (
     double K, double *F, double *C,double *X,double *Y,double *Z,line_t *Line
  )
{
    double f[3],c[3],x[3],y[3],z[3];

    int i, j, k, n=0, above=0;

    static int map[6][2] =
      {
         { 0,3 }, { 3,1 }, { 1,4 }, { 4,2 }, { 2,5 }, { 5,0 }
      };

    for( i=0; i<6; i++ ) above += F[i]>K;
    if ( above == 0 || above == 6 ) return 0;

    f[0] = elm_6node_triangle_fvalue( F,1.0/3.0,1.0/3.0 );
    c[0] = elm_6node_triangle_fvalue( C,1.0/3.0,1.0/3.0 );
    x[0] = elm_6node_triangle_fvalue( X,1.0/3.0,1.0/3.0 );
    y[0] = elm_6node_triangle_fvalue( Y,1.0/3.0,1.0/3.0 );
    z[0] = elm_6node_triangle_fvalue( Z,1.0/3.0,1.0/3.0 );

    for( i=0; i<6; i++ )
    {
        for( j=1; j<3; j++ )
        {
            k = map[i][j-1];
            f[j] = F[k];
            c[j] = C[k];
            x[j] = X[k];
            y[j] = Y[k];
            z[j] = Z[k];
        }
        n += elm_3node_triangle_isoline( K,f,c,x,y,z,&Line[n] );
    }

    return n;
}


/*******************************************************************************
 *
 *     Name:        elm_6node_triangle_initialize( )
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
int elm_6node_triangle_initialize()
{
     static char *Name = "ELM_6NODE_TRIANGLE";

     element_type_t ElementDef;

     elm_6node_triangle_shape_functions();

     ElementDef.ElementName = Name;
     ElementDef.ElementCode = 306;

     ElementDef.NumberOfNodes = 6;

     ElementDef.NodeU = NodeU;
     ElementDef.NodeV = NodeV;
     ElementDef.NodeW = NULL;

     ElementDef.PartialU = (double (*)())elm_6node_triangle_dndu_fvalue;
     ElementDef.PartialV = (double (*)())elm_6node_triangle_dndv_fvalue;
     ElementDef.PartialW = (double (*)())NULL;

     ElementDef.FunctionValue = (double (*)())elm_6node_triangle_fvalue;
     ElementDef.Triangulate   = (int (*)())elm_6node_triangle_triangulate;
     ElementDef.PointInside   = (int (*)())elm_6node_triangle_point_inside;
     ElementDef.IsoLine       = (int (*)())elm_6node_triangle_isoline;
     ElementDef.IsoSurface    = (int (*)())NULL;

     return elm_add_element_type( &ElementDef ) ;
}
