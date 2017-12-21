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
 * Definition of 9 node quad element.
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
 * Modification history:
 *
 * 28 Sep 1995, changed call to elm_triangle_normal to geo_triangle normal
 *              routine elm_... doesn't exist anymore
 * Juha R
 *
 *
 ******************************************************************************/

#include "../elmerpost.h"
#include <elements.h>

/*
 *
 * shape functions: N0 = ( uv(1-u-v)+u^2+v^2-1)/4
 *                  N1 = (-uv(1+u-v)+u^2+v^2-1)/4
 *                  N2 = ( uv(1+u+v)+u^2+v^2-1)/4      3---6---2
 *                  N3 = (-uv(1-u+v)+u^2+v^2-1)/4      |       |
 *                  N4 = ((1-v)(1-u^2)/2               7       5
 *                  N5 = ((1+u)(1-v^2)/2             v |       |
 *                  N6 = ((1+v)(1-u^2)/2               0---4---1
 *                  N7 = ((1-u)(1-v^2)/2                 u
 *
 */

static double N[9][9],A[9][9];

static double NodeU[] = { -1.0,  1.0, 1.0, -1.0,  0.0, 1.0, 0.0, -1.0, 0.0 };
static double NodeV[] = { -1.0, -1.0, 1.0,  1.0, -1.0, 0.0, 1.0,  0.0, 0.0 };

/*******************************************************************************
 *
 *     Name:        elm_9node_quad_shape_functions( )
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
static void elm_9node_quad_shape_functions()
{
     double u,v;

     int i,j;

     for( i=0; i<9; i++ )
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
         A[i][8]  = u*u*v*v;
     }

     lu_mtrinv( (double *)A,9 );

     for( i=0; i<9; i++ )
        for( j=0; j<9; j++ ) N[i][j] = A[j][i];
}


/*******************************************************************************
 *
 *     Name:        elm_9node_quad_triangulate( geometry_t *,element_t * )
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
int elm_9node_quad_triangulate( geometry_t *geom,element_t *Elm )
{
    triangle_t triangle;

    triangle.Element = Elm;

    triangle.v[0] = Elm->Topology[0];
    triangle.v[1] = Elm->Topology[4];
    triangle.v[2] = Elm->Topology[8];

    triangle.Edge[0] = TRUE;
    triangle.Edge[1] = FALSE;
    triangle.Edge[2] = FALSE;

    if ( !geo_add_triangle( geom,&triangle ) )  return FALSE;

    triangle.v[0] = Elm->Topology[4];
    triangle.v[1] = Elm->Topology[1];
    triangle.v[2] = Elm->Topology[8];

    triangle.Edge[0] = TRUE; 
    triangle.Edge[1] = FALSE;
    triangle.Edge[2] = FALSE;

    if ( !geo_add_triangle( geom,&triangle ) )  return FALSE;

    triangle.v[0] = Elm->Topology[1];
    triangle.v[1] = Elm->Topology[5];
    triangle.v[2] = Elm->Topology[8];

    triangle.Edge[0] = TRUE; 
    triangle.Edge[1] = FALSE;
    triangle.Edge[2] = FALSE;

    if ( !geo_add_triangle( geom,&triangle ) )  return FALSE;

    triangle.v[0] = Elm->Topology[5];
    triangle.v[1] = Elm->Topology[2];
    triangle.v[2] = Elm->Topology[8];

    triangle.Edge[0] = TRUE; 
    triangle.Edge[1] = FALSE;
    triangle.Edge[2] = FALSE;

    if ( !geo_add_triangle( geom,&triangle ) )  return FALSE;

    triangle.v[0] = Elm->Topology[2];
    triangle.v[1] = Elm->Topology[6];
    triangle.v[2] = Elm->Topology[8];

    triangle.Edge[0] = TRUE; 
    triangle.Edge[1] = FALSE;
    triangle.Edge[2] = FALSE;

    if ( !geo_add_triangle( geom,&triangle ) )  return FALSE;

    triangle.v[0] = Elm->Topology[6];
    triangle.v[1] = Elm->Topology[3];
    triangle.v[2] = Elm->Topology[8];

    triangle.Edge[0] = TRUE; 
    triangle.Edge[1] = FALSE;
    triangle.Edge[2] = FALSE;

    if ( !geo_add_triangle( geom,&triangle ) )  return FALSE;

    triangle.v[0] = Elm->Topology[3];
    triangle.v[1] = Elm->Topology[7];
    triangle.v[2] = Elm->Topology[8];

    triangle.Edge[0] = TRUE; 
    triangle.Edge[1] = FALSE;
    triangle.Edge[2] = FALSE;

    if ( !geo_add_triangle( geom,&triangle ) )  return FALSE;

    triangle.v[0] = Elm->Topology[7];
    triangle.v[1] = Elm->Topology[0];
    triangle.v[2] = Elm->Topology[8];

    triangle.Edge[0] = TRUE; 
    triangle.Edge[1] = FALSE;
    triangle.Edge[2] = FALSE;

    if ( !geo_add_triangle( geom,&triangle ) )  return FALSE;

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        elm_9node_quad_fvalue( double *,double,double )
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
static double elm_9node_quad_fvalue( double *F, double u, double v )
{
     double R,uv=u*v,uu=u*u,vv=v*v,uuv=uu*v,uvv=u*vv;
     int i;

     R = 0.0;
     for( i=0; i<9; i++ )
     {
          R += F[i]*(N[i][0]      +
                     N[i][1]*u    +
                     N[i][2]*v    +
                     N[i][3]*uv   +
                     N[i][4]*uu   +
                     N[i][5]*vv   +
                     N[i][6]*uuv  +
                     N[i][7]*uvv  +
                     N[i][8]*uu*vv);
     }

     return R;
}
     
/*******************************************************************************
 *
 *     Name:        elm_9node_quad_dndu_fvalue( double *,double,double )
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
static double elm_9node_quad_dndu_fvalue( double *F, double u, double v )
{
     double R=0.0,uv=u*v,vv=v*v;
     int i;

     for( i=0; i<9; i++ )
     {
         R += F[i]*(N[i][1]+N[i][3]*v+2*N[i][4]*u+2*N[i][6]*uv+N[i][7]*vv+2*N[i][8]*u*vv);
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_9node_quad_dndv_fvalue( double *,double,double )
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
static double elm_9node_quad_dndv_fvalue( double *F, double u, double v )
{
     double R=0.0,uu=u*u,uv=u*v;

     int i;

     for( i=0; i<9; i++ )
     {
         R += F[i]*(N[i][2]+N[i][3]*u+2*N[i][5]*v+N[i][6]*uu+2*N[i][7]*uv+2*N[i][8]*uu*v);
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:     elm_9node_quad_dd(double *F,double *u,double *v,double *Values)
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
static void elm_9node_quad_ddu( double *F,double u,double v,double *Values )
{
    double ddu = 0.0, dudv=0.0, ddv = 0.0;

    int i;

    for( i=0; i<9; i++ )
    {
        ddu  += F[i]*( 2*N[i][4] + 2*v*N[i][6] + 2*v*v*N[i][8] );
        ddv  += F[i]*( 2*N[i][5] + 2*v*N[i][7] + 2*u*u*N[i][8] );
        dudv += F[i]*( N[i][3] + 2*u*N[i][6] + 2*v*N[i][7] + 4*u*v*N[i][8] );
    }

    Values[0] =  ddu;
    Values[1] = dudv;
    Values[2] = dudv;
    Values[3] =  ddv;
}

/*******************************************************************************
 *
 *     Name:        elm_9node_quad_point_inside(
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
int elm_9node_quad_point_inside
   (
                 double *nx, double *ny, double *nz,
        double px, double py, double pz, double *u,double *v,double *w
   )
{
    double lx[4],ly[4],lz[4],cx,cy,cz;

    int i,j,k;

    static int map[4][3] =
    {
       { 7,0,4 }, { 4,1,5 }, { 5,2,6 }, { 6,3,7 }
    };

    cx = elm_9node_quad_fvalue( nx,0.0,0.0 );
    cy = elm_9node_quad_fvalue( ny,0.0,0.0 );
    cz = elm_9node_quad_fvalue( nz,0.0,0.0 );

    for( i=0; i<4; i++ )
    {

        switch( i )
        {
            case 0:
                lx[0] = nx[0]; ly[0] = ny[0]; lz[0] = nz[0];
                lx[1] = nx[4]; ly[1] = ny[4]; lz[1] = nz[4];
                lx[2] = cx; ly[2] = cy; lz[2] = cz;
                lx[3] = nx[7]; ly[3] = ny[7]; lz[3] = nz[7];
            break;
            case 1:
                lx[0] = nx[4]; ly[0] = ny[4]; lz[0] = nz[4];
                lx[1] = nx[1]; ly[1] = ny[1]; lz[1] = nz[1];
                lx[2] = nx[5]; ly[2] = ny[5]; lz[2] = nz[5];
                lx[3] = cx; ly[3] = cy; lz[3] = cz;
            break;
            case 2:
                lx[0] = cx; ly[0] = cy; lz[0] = cz;
                lx[1] = nx[5]; ly[1] = ny[5]; lz[1] = nz[5];
                lx[2] = nx[2]; ly[2] = ny[2]; lz[2] = nz[2];
                lx[3] = nx[6]; ly[3] = ny[6]; lz[3] = nz[6];
            break;
            case 3:
                lx[0] = nx[7]; ly[0] = ny[7]; lz[0] = nz[7];
                lx[1] = cx; ly[1] = cy; lz[1] = cz;
                lx[2] = nx[6]; ly[2] = ny[6]; lz[2] = nz[6];
                lx[3] = nx[3]; ly[3] = ny[3]; lz[3] = nz[3];
            break;
        }

        if ( elm_4node_quad_point_inside( lx,ly,lz,px,py,pz,u,v,w ) ) 
        {
            switch( i )
            {
                case 0:
                    *u = *u/2 - 0.5;
                    *v = *v/2 - 0.5;
                break;

                case 1:
                    *u = *u/2 + 0.5;
                    *v = *v/2 - 0.5;
                break;

                case 2:
                    *u = *u/2 + 0.5;
                    *v = *v/2 + 0.5;
                break;

                case 3:
                    *u = *u/2 - 0.5;
                    *v = *v/2 + 0.5;
                break;
            }

            return TRUE;
        }
    }

    return FALSE;
}

/*******************************************************************************
 *
 *     Name:        elm_9node_quad_isoline
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
 *   Return value:  number of lines generated (0,...,16)
 *
 ******************************************************************************/
int elm_9node_quad_isoline
  (
     double K, double *F, double *C,double *X,double *Y,double *Z,line_t *Line
  )
{
    double f[4],c[4],x[4],y[4],z[4];

    int i, j, k, n=0, above=0;

    static int map[4][3] =
    {
       { 7,0,4 }, { 4,1,5 }, { 5,2,6 }, { 6,3,7 }
    };

    for( i=0; i<9; i++ ) above += F[i]>K;
    if ( above == 0 || above == 9 ) return 0;

    f[0] = elm_9node_quad_fvalue( F,0.0,0.0 );
    c[0] = elm_9node_quad_fvalue( C,0.0,0.0 );
    x[0] = elm_9node_quad_fvalue( X,0.0,0.0 );
    y[0] = elm_9node_quad_fvalue( Y,0.0,0.0 );
    z[0] = elm_9node_quad_fvalue( Z,0.0,0.0 );

    for( i=0; i<4; i++ )
    {
        for( j=1; j<4; j++ )
        {
            k = map[i][j-1];
            f[j] = F[k];
            c[j] = C[k];
            x[j] = X[k];
            y[j] = Y[k];
            z[j] = Z[k];
        }
        n += elm_4node_quad_isoline( K,f,c,x,y,z,&Line[n] );
    }

    return n;
}


/******************************************************************************
 *
 *     Name:        elm_9node_quad_initialize( )
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
int elm_9node_quad_initialize()
{
     element_type_t ElementDef;

     static char *Name = "ELM_9NODE_QUAD";

     elm_9node_quad_shape_functions();

     ElementDef.ElementName = Name;
     ElementDef.ElementCode = 409;

     ElementDef.NumberOfNodes = 9;
     ElementDef.NodeU = NodeU;
     ElementDef.NodeV = NodeV;
     ElementDef.NodeW = NULL;

     ElementDef.PartialU = (double (*)())elm_9node_quad_dndu_fvalue;
     ElementDef.PartialV = (double (*)())elm_9node_quad_dndv_fvalue;
     ElementDef.PartialW = NULL;
     ElementDef.SecondPartials = (void (*)())elm_9node_quad_ddu;

     ElementDef.FunctionValue = (double (*)())elm_9node_quad_fvalue;

     ElementDef.Triangulate = (int (*)())elm_9node_quad_triangulate;
     ElementDef.PointInside = (int (*)())elm_9node_quad_point_inside;
     ElementDef.IsoLine       = (int (*)())elm_9node_quad_isoline;
     ElementDef.IsoSurface    = (int (*)())NULL;

     return elm_add_element_type( &ElementDef );
}
