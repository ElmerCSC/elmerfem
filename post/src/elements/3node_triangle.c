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
 * Definition of three node triangle element.
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
 * 28 Sep 1995, changed call to elm_triangle_normal to geo_triangle normal
 *              routine elm_... doesn't exist anymore
 *
 ******************************************************************************/

#include "../elmerpost.h"
#include <elements.h>

/*
 * Three node triangular surface element
 * 
 *         2 v=1
 *        / \
 *       /   \        
 *      /     \       
 *     0-------1 u=1
 *
 */
static double N[3][3] =
{
    { 1.0,-1.0,-1.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 }
};

static double NodeU[3] = { 0.0, 1.0, 0.0 };
static double NodeV[3] = { 0.0, 0.0, 1.0 };

/*******************************************************************************
 *
 *     Name:        elm_3node_triangle_triangulate( geometry_t *,element_t * )
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
int elm_3node_triangle_triangulate( geometry_t *geom, element_t *Elm, element_t *Parent)
{
    double u,v,w,s;

    int i,j;

    triangle_t triangle;

    triangle.Element = Parent;
    for( i=0; i<3; i++ )
    {
        triangle.v[i] = Elm->Topology[i];
        triangle.Edge[i] = TRUE;
    }
    geo_triangle_normal( geom,&triangle );

    return geo_add_triangle( geom, &triangle );
}

/*******************************************************************************
 *
 *     Name:        elm_3node_triangle_point_inside(
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
int elm_3node_triangle_point_inside
   ( 
                   double *nx, double *ny, double *nz,
         double px, double py, double pz, double *u,double *v,double *w
   )
{
    double A00,A01,A10,A11,B00,B01,B10,B11,detA;

    if ( px > MAX(MAX( nx[0],nx[1] ),nx[2] ) ) return FALSE;
    if ( px < MIN(MIN( nx[0],nx[1] ),nx[2] ) ) return FALSE;

    if ( py > MAX(MAX( ny[0],ny[1] ),ny[2] ) ) return FALSE;
    if ( py < MIN(MIN( ny[0],ny[1] ),ny[2] ) ) return FALSE;

#if 0
    if ( pz > MAX(MAX( nz[0],nz[1] ),nz[2] ) ) return FALSE;
    if ( pz < MIN(MIN( nz[0],nz[1] ),nz[2] ) ) return FALSE;
#endif

    A00 = nx[1] - nx[0];
    A01 = nx[2] - nx[0];
    A10 = ny[1] - ny[0];
    A11 = ny[2] - ny[0];

    detA = A00*A11 - A01*A10;
    if ( ABS(detA) < AEPS )
    {
        fprintf( stderr, "3node_triangle_inside: SINGULAR, HUH? \n" );
        return FALSE;
    }

    detA = 1/detA;

    B00 =  A11*detA;
    B10 = -A10*detA;
    B01 = -A01*detA;
    B11 =  A00*detA;

    px -= nx[0];
    py -= ny[0];
    *u = *v = *w = 0.0;

    *u = B00*px + B01*py;
    if ( *u < 0.0 || *u > 1.0 ) return FALSE;

    *v = B10*px + B11*py;
    if ( *v < 0.0 || *v > 1.0 ) return FALSE;

    return (*u + *v) <= 1.0;
}

/*******************************************************************************
 *
 *     Name:        elm_3node_triangle_fvalue( double *,double,double )
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
static double elm_3node_triangle_fvalue(double *F,double u,double v)
{
    return F[0]*(1-u-v) + F[1]*u  + F[2]*v;
}

/*******************************************************************************
 *
 *     Name:        elm_3node_triangle_dndu_fvalue( double *,double,double )
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
static double elm_3node_triangle_dndu_fvalue(double *F,double u,double v)
{
    return -F[0] + F[1];
}

/*******************************************************************************
 *
 *     Name:        elm_3node_triangle_dndv_fvalue( double *,double,double )
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
static double elm_3node_triangle_dndv_fvalue(double *F,double u,double v)
{
     return -F[0] + F[2];
}

#define FEPS 1.0E-9

/*******************************************************************************
 *
 *     Name:        elm_3node_triangle_isoline
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
 *   Return value:  number of lines generated (0 or 1)
 *
*******************************************************************************/
int elm_3node_triangle_isoline
  (
     double K, double *F, double *C,double *X,double *Y,double *Z,line_t *Line
  )
{
    double F0=F[0],F1=F[1],F2=F[2],t;
    int i,n=0;

    int S0 = (F0 > K);
    int S1 = (F1 > K);
    int S2 = (F2 > K);
    int S = S0 + S1 + S2;

    if ( S==0 || S==3 ) return 0;

    if (S0 ^ S1)
    {
        if ( ABS(F1-F0) < FEPS ) return 0;
        t = (K-F0)/(F1-F0);

        Line->x[n] = t*(X[1] - X[0]) + X[0];
        Line->y[n] = t*(Y[1] - Y[0]) + Y[0];
        Line->z[n] = t*(Z[1] - Z[0]) + Z[0];
        Line->f[n] = K;
        if ( C ) Line->c[n] = t*(C[1] - C[0]) + C[0];
        n++;
    }

    if (S0 ^ S2)
    {
        if ( ABS(F2-F0) < FEPS ) return 0;
        t = (K-F0)/(F2-F0);

        Line->x[n] = t*(X[2] - X[0]) + X[0];
        Line->y[n] = t*(Y[2] - Y[0]) + Y[0];
        Line->z[n] = t*(Z[2] - Z[0]) + Z[0];
        Line->f[n] = K;
        if ( C ) Line->c[n] = t*(C[2] - C[0]) + C[0];
        n++;
    }

    if (S1 ^ S2)
    {
        if ( ABS(F2-F1) < FEPS ) return 0;
        t = (K-F1)/(F2-F1);

        Line->x[n] = t*(X[2] - X[1]) + X[1];
        Line->y[n] = t*(Y[2] - Y[1]) + Y[1];
        Line->z[n] = t*(Z[2] - Z[1]) + Z[1];
        Line->f[n] = K;
        if ( C ) Line->c[n] = t*(C[2] - C[1]) + C[1];
        n++;
    }

    return (n>=2);
}


/*******************************************************************************
 *
 *     Name:        elm_3node_triangle_initialize()
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
int elm_3node_triangle_initialize()
{
     static char *Name = "ELM_3NODE_TRIANGLE";

     element_type_t ElementDef;

     ElementDef.ElementName = Name;
     ElementDef.ElementCode = 303;

     ElementDef.NumberOfNodes = 3;

     ElementDef.NodeU = NodeU;
     ElementDef.NodeV = NodeV;
     ElementDef.NodeW = NULL;

     ElementDef.PartialU = (double (*)())elm_3node_triangle_dndu_fvalue;
     ElementDef.PartialV = (double (*)())elm_3node_triangle_dndv_fvalue;
     ElementDef.PartialW = (double (*)())NULL;

     ElementDef.FunctionValue = (double (*)())elm_3node_triangle_fvalue;
     ElementDef.Triangulate   = (int (*)())elm_3node_triangle_triangulate;
     ElementDef.IsoLine       = (int (*)())elm_3node_triangle_isoline;
     ElementDef.PointInside   = (int (*)())elm_3node_triangle_point_inside;
     ElementDef.IsoSurface    = (int (*)())NULL;

     return elm_add_element_type( &ElementDef ) ;
}
