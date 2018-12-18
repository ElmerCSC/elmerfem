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
 * Definition of 5 node quad element.
 *
 *******************************************************************************
 *
 *                     Author:       Juha Ruokolainen
 *
 *                    Address: CSC - IT Center for Science Ltd.
 *                                Keilaranta 14, P.O. BOX 405
 *                                  02101 Espoo, Finland
 *                                  Tel. +355 0 457 2723
 *                                Telefax: +355 0 457 2302
 *                              EMail: Juha.Ruokolainen@csc.fi
 *
 *                       Date: 20 Sep 1995
 *
 * Modification history:
 *
 * 25 Sep 1995, changed call to elm_triangle_normal to geo_triangle normal
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

static double N[5][5],A[5][5];

static double NodeU[] = { -1.0,  1.0, 1.0, -1.0,  0.0 };
static double NodeV[] = { -1.0, -1.0, 1.0,  1.0,  0.0,};

/*******************************************************************************
 *
 *     Name:        elm_5node_quad_shape_functions( )
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
static void elm_5node_quad_shape_functions()
{
     double u,v;

     int i,j;

     for( i=0; i<5; i++ )
     {
         u = NodeU[i];
         v = NodeV[i];

         A[i][0]  = 1;
         A[i][1]  = u;
         A[i][2]  = v;
         A[i][3]  = u*v;
         A[i][4]  = u*u*v*v;
     }

     lu_mtrinv( (double *)A,5 );

     for( i=0; i<5; i++ )
        for( j=0; j<5; j++ ) N[i][j] = A[j][i];
}


/*******************************************************************************
 *
 *     Name:        elm_5node_quad_triangulate( geometry_t *,element_t * )
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
int elm_5node_quad_triangulate( geometry_t *geom,element_t *Elm )
{
    triangle_t triangle;

    triangle.Element = Elm;

    triangle.v[0] = Elm->Topology[0];
    triangle.v[1] = Elm->Topology[1];
    triangle.v[2] = Elm->Topology[4];

    triangle.Edge[0] = TRUE;
    triangle.Edge[1] = FALSE;
    triangle.Edge[2] = FALSE;

    if ( !geo_add_triangle( geom,&triangle ) )  return FALSE;

    triangle.v[0] = Elm->Topology[1];
    triangle.v[1] = Elm->Topology[2];
    triangle.v[2] = Elm->Topology[4];

    triangle.Edge[0] = TRUE; 
    triangle.Edge[1] = FALSE;
    triangle.Edge[2] = FALSE;

    if ( !geo_add_triangle( geom,&triangle ) )  return FALSE;

    triangle.v[0] = Elm->Topology[2];
    triangle.v[1] = Elm->Topology[3];
    triangle.v[2] = Elm->Topology[4];

    triangle.Edge[0] = TRUE; 
    triangle.Edge[1] = FALSE;
    triangle.Edge[2] = FALSE;

    if ( !geo_add_triangle( geom,&triangle ) )  return FALSE;

    triangle.v[0] = Elm->Topology[3];
    triangle.v[1] = Elm->Topology[0];
    triangle.v[2] = Elm->Topology[4];

    triangle.Edge[0] = TRUE; 
    triangle.Edge[1] = FALSE;
    triangle.Edge[2] = FALSE;

    if ( !geo_add_triangle( geom,&triangle ) )  return FALSE;

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        elm_5node_quad_fvalue( double *,double,double )
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
static double elm_5node_quad_fvalue( double *F, double u, double v )
{
     double R,uv=u*v,uu=u*u,vv=v*v,uuv=uu*v,uvv=u*vv;
     int i;

     R = 0.0;
     for( i=0; i<5; i++ )
     {
          R += F[i]*(N[i][0]      +
                     N[i][1]*u    +
                     N[i][2]*v    +
                     N[i][3]*uv   +
                     N[i][4]*uu*vv );
     }

     return R;
}
     
/*******************************************************************************
 *
 *     Name:        elm_5node_quad_dndu_fvalue( double *,double,double )
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
static double elm_5node_quad_dndu_fvalue( double *F, double u, double v )
{
     double R=0.0,uv=u*v,vv=v*v;
     int i;

     for( i=0; i<5; i++ )
     {
         R += F[i]*(N[i][1]+N[i][3]*v+2*N[i][4]*u*v*v );
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_5node_quad_dndv_fvalue( double *,double,double )
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
static double elm_5node_quad_dndv_fvalue( double *F, double u, double v )
{
     double R=0.0,uu=u*u,uv=u*v;

     int i;

     for( i=0; i<5; i++ )
     {
         R += F[i]*(N[i][2]+N[i][3]*u+2*N[i][4]*uu*v );
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:     elm_5node_quad_dd(double *F,double *u,double *v,double *Values)
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
static void elm_5node_quad_ddu( double *F,double u,double v,double *Values )
{
    double ddu = 0.0, dudv=0.0, ddv = 0.0;

    int i;

    for( i=0; i<5; i++ )
    {
        ddu  += F[i]*( 2*N[i][4]*v*v );
        ddv  += F[i]*( 2*N[i][4]*u*u );
        dudv += F[i]*( N[i][3] + 4*N[i][4]*u*v );
    }

    Values[0] =  ddu;
    Values[1] = dudv;
    Values[2] = dudv;
    Values[3] =  ddv;
}

/*******************************************************************************
 *
 *     Name:        elm_5node_quad_point_inside(
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
int elm_5node_quad_point_inside
   (
                 double *nx, double *ny, double *nz,
        double px, double py, double pz, double *u,double *v,double *w
   )
{
    return elm_4node_quad_point_inside( nx,ny,nz,px,py,pz,u,v,w );
}

/*******************************************************************************
 *
 *     Name:        elm_5node_quad_isoline
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
 *   Return value:  number of lines generated (0,...,4)
 *
 ******************************************************************************/
int elm_5node_quad_isoline
  (
     double K, double *F, double *C,double *X,double *Y,double *Z,line_t *Line
  )
{
    double f[4],c[4],x[4],y[4],z[4];

    int i, j, k, n=0, above=0;

    static int map[4][2] =
    {
       { 0,1 }, { 1,2 }, { 2,3 }, { 3,0 }
    };

    for( i=0; i<5; i++ ) above += F[i]>K;
    if ( above == 0 || above == 5 ) return 0;

    f[2] = F[4];
    c[2] = C[4];
    x[2] = X[4];
    y[2] = Y[4];
    z[2] = Z[4];

    for( i=0; i<4; i++ )
    {
        for( j=1; j<2; j++ )
        {
            k = map[i][j];
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


/******************************************************************************
 *
 *     Name:        elm_5node_quad_initialize( )
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
int elm_5node_quad_initialize()
{
     element_type_t ElementDef;

     static char *Name = "ELM_5NODE_QUAD";

     elm_5node_quad_shape_functions();

     ElementDef.ElementName = Name;
     ElementDef.ElementCode = 405;

     ElementDef.NumberOfNodes = 5;
     ElementDef.NodeU = NodeU;
     ElementDef.NodeV = NodeV;
     ElementDef.NodeW = NULL;

     ElementDef.PartialU = (double (*)())elm_5node_quad_dndu_fvalue;
     ElementDef.PartialV = (double (*)())elm_5node_quad_dndv_fvalue;
     ElementDef.PartialW = NULL;
     ElementDef.SecondPartials = (void (*)())elm_5node_quad_ddu;

     ElementDef.FunctionValue = (double (*)())elm_5node_quad_fvalue;

     ElementDef.Triangulate = (int (*)())elm_5node_quad_triangulate;
     ElementDef.PointInside = (int (*)())elm_5node_quad_point_inside;
     ElementDef.IsoLine       = (int (*)())elm_5node_quad_isoline;
     ElementDef.IsoSurface    = (int (*)())NULL;

     return elm_add_element_type( &ElementDef );
}
