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
 * Definition of ten node triangle element.
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
 * TEN NODE TRIANGLE ELEMENT
 *
 */

static double N[10][10],A[10][10];

static double NodeU[] = { 0.0, 1.0, 0.0, 1./3., 2./3., 2./3., 1./3., 0.0, 0.0,    1./3. };
static double NodeV[] = { 0.0, 0.0, 1.0, 0.0,   0.0,   1./3., 2./3., 2./3, 1./3., 1./3. };

/*******************************************************************************
 *
 *     Name:        elm_10node_triangle_shape_functions( )
 *
 *     Purpose:     Initialize element shape function array. Internal only.
 *
 *     Parameters:
 *
 *         Input:   Global (filewise) variables NodeU,NodeV,NodeW
 *
 *         Output:  Global (filewise) variable N[10][10], will contain
 *                  shape function coefficients
 *
 *   Return value:  void
 *
 ******************************************************************************/
static void elm_10node_triangle_shape_functions()
{
     double u,v;

     int i,j;

     for( i=0; i<10; i++ )
     {
         u = NodeU[i];
         v = NodeV[i];


         A[i][0]  = 1;
         A[i][1]  = u;
         A[i][2]  = v;
         A[i][3]  = u*v;
         A[i][4]  = u*u;
         A[i][5]  = v*v;
         A[i][6]  = u*u*v;
         A[i][7]  = u*v*v;
         A[i][8]  = u*u*u;
         A[i][9]  = v*v*v;
      
     }

     lu_mtrinv( (double *)A,10 );

     for( i=0; i<10; i++ )
        for( j=0; j<10; j++ ) N[i][j] = A[j][i];
}

/*******************************************************************************
 *
 *     Name:        elm_10node_triangle_triangulate( geometry_t *,element_t * )
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
int elm_10node_triangle_triangulate( geometry_t *geom, element_t *Elm )
{
    double u,v,w,s;

    static int map[9][3] =
      {
         { 0,3,8 }, { 3,4,9 }, { 3,9,8 },
         { 4,5,9 }, { 4,1,5 }, { 5,6,9 },
         { 8,9,7 }, { 9,6,7 }, { 7,6,2 }
      };
    static int edge_map[9][3] =
      {
         { 1,0,1 }, { 1,0,0 }, { 0,0,0 },
         { 0,0,0 }, { 1,1,0 }, { 1,0,0 },
         { 0,0,1 }, { 0,0,0 }, { 0,1,1 }
      };

    int i,j;

    triangle_t triangle;

    triangle.Element = Elm;
    for( i=0; i<9; i++ )
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
 *     Name:        elm_10node_triangle_point_inside(
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
int elm_10node_triangle_point_inside
   ( 
                   double *nx, double *ny, double *nz,
         double px, double py, double pz, double *u,double *v,double *w
   )
{
    return elm_3node_triangle_point_inside( nx,ny,nz,px,py,pz,u,v,w );
}

/*******************************************************************************
 *
 *     Name:        elm_10node_triangle_fvalue( double *,double,double )
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
static double elm_10node_triangle_fvalue( double *F,double u,double v )
{
     double R=0.0,uv=u*v,u2=u*u,v2=v*v;
     int i;


     for( i=0; i<10; i++ )
     {
         R += F[i]*( N[i][0]    +
                     N[i][1]*u  +
                     N[i][2]*v  +
                     N[i][3]*uv +
                     N[i][4]*u2 +
                     N[i][5]*v2 +
                     N[i][6]*u2*v +
                     N[i][7]*u*v2 +
                     N[i][8]*u2*u +
                     N[i][9]*v2*v );
     } 

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_10node_triangle_dndu_fvalue( double *,double,double )
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
static double elm_10node_triangle_dndu_fvalue(double *F,double u,double v)
{
     double R=0.0;
     int i;

     for( i=0; i<10; i++ )
     {
         R += F[i]*( N[i][1] + N[i][3]*v + 2*N[i][4]*u +
               N[i][6]*2*u*v + N[i][7]*v*v + 3*N[i][8]*u*u );
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_10node_triangle_dndv_fvalue( double *,double,double )
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
static double elm_10node_triangle_dndv_fvalue(double *F,double u,double v)
{
     double R=0.0;
     int i;

     for( i=0; i<10; i++ )
     {
         R += F[i]*( N[i][2] + N[i][3]*u + 2*N[i][5]*v + N[i][6]*u*u +
                     N[i][7]*2*u*v + 3*N[i][9]*v*v);
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_10node_triangle_isoline
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
int elm_10node_triangle_isoline
  (
     double K, double *F, double *C,double *X,double *Y,double *Z,line_t *Line
  )
{
    double f[3],c[3],x[3],y[3],z[3];

    int i, j, k, n=0, above=0;

    static int map[9][3] =
      {
         { 0,3,8 }, { 3,4,9 }, { 3,9,8 },
         { 4,5,9 }, { 4,1,5 }, { 5,6,9 },
         { 8,9,7 }, { 9,6,7 }, { 7,6,2 }
      };

    for( i=0; i<10; i++ ) above += F[i]>K;
    if ( above == 0 || above == 10 ) return 0;

    for( i=0; i<9; i++ )
    {
        for( j=0; j<3; j++ )
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
 *     Name:        elm_10node_triangle_initialize( )
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
int elm_10node_triangle_initialize()
{
     static char *Name = "ELM_10NODE_TRIANGLE";

     element_type_t ElementDef;

     elm_10node_triangle_shape_functions();

     ElementDef.ElementName = Name;
     ElementDef.ElementCode = 310;

     ElementDef.NumberOfNodes = 10;

     ElementDef.NodeU = NodeU;
     ElementDef.NodeV = NodeV;
     ElementDef.NodeW = NULL;

     ElementDef.PartialU = (double (*)())elm_10node_triangle_dndu_fvalue;
     ElementDef.PartialV = (double (*)())elm_10node_triangle_dndv_fvalue;
     ElementDef.PartialW = (double (*)())NULL;

     ElementDef.FunctionValue = (double (*)())elm_10node_triangle_fvalue;
     ElementDef.Triangulate   = (int (*)())elm_10node_triangle_triangulate;
     ElementDef.PointInside   = (int (*)())elm_10node_triangle_point_inside;
     ElementDef.IsoLine       = (int (*)())elm_10node_triangle_isoline;
     ElementDef.IsoSurface    = (int (*)())NULL;

     return elm_add_element_type( &ElementDef ) ;
}
