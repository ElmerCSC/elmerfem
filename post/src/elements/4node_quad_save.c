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
 * Definition of 4 node quad element.
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
#include "elements.h"

/*
 * Four node quadrilater surface element
 *
 *    3----------2
 *    |          |
 *    |          |
 *  v |          |
 *    0----------1
 *      u
 */

static double N[4][4] =
{
     { 1.0/4.0, -1.0/4.0, -1.0/4.0,  1.0/4.0 },  /* 1 - u - v + uv */
     { 1.0/4.0,  1.0/4.0, -1.0/4.0, -1.0/4.0 },
     { 1.0/4.0,  1.0/4.0,  1.0/4.0,  1.0/4.0 },
     { 1.0/4.0, -1.0/4.0,  1.0/4.0, -1.0/4.0 }
};

static double NodeU[] = { -1.0,  1.0, 1.0, -1.0 };
static double NodeV[] = { -1.0, -1.0, 1.0,  1.0 };

/*******************************************************************************
 *
 *     Name:        elm_4node_quad_triangulate( geometry_t *,element_t * )
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
int elm_4node_quad_triangulate( geometry_t *geom,element_t *Elm )
{
    double x[4],y[4],z[4],f[4][4],elm_4node_quad_fvalue();

    triangle_t triangle;
    vertex_t vertex;

    int i,j,n;

    material_t *material=NULL;

#if 0
    for( i=0; i<4; i++ )
    {
        x[i] = geom->Vertices[Elm->Topology[i]].x[0];
        y[i] = geom->Vertices[Elm->Topology[i]].x[1];
        z[i] = geom->Vertices[Elm->Topology[i]].x[2];

        for( j=0; j<4; j++ )
        {
            f[j][i] = FuncArray[j][Elm->Topology[i]];
        }
    }

    vertex.x[0] = elm_4node_quad_fvalue(x,0.0,0.0);
    vertex.x[1] = elm_4node_quad_fvalue(y,0.0,0.0);
    vertex.x[2] = elm_4node_quad_fvalue(z,0.0,0.0);
    vertex.ElementModelNode = FALSE;

    if ( (n=geo_add_vertex( geom,&vertex )) < 0 ) return FALSE;

    for( j=0; j<4; j++ )
    {
        FuncArray[j][n] = elm_4node_quad_fvalue(f[j],0.0,0.0);
    }

    for( i=0; i<4; i++ )
    {
        triangle.v[0] = Elm->Topology[i];
        if ( i==3 )
            triangle.v[1] = Elm->Topology[0];
        else
            triangle.v[1] = Elm->Topology[i+1];
        triangle.v[2] = n;

        geo_triangle_normal( geom,&triangle );

        triangle.Edge[0] = TRUE;
        triangle.Edge[1] = FALSE;
        triangle.Edge[2] = FALSE;

        if ( !geo_add_triangle( geom,&triangle,material ) )  return FALSE;
    }
#else
    triangle.v[0] = Elm->Topology[0];
    triangle.v[1] = Elm->Topology[1];
    triangle.v[2] = Elm->Topology[2];

    geo_triangle_normal( geom,&triangle );

    triangle.Edge[0] = TRUE;
    triangle.Edge[1] = TRUE;
    triangle.Edge[2] = FALSE;

    if ( !geo_add_triangle( geom,&triangle ) )  return FALSE;

    triangle.v[0] = Elm->Topology[0];
    triangle.v[1] = Elm->Topology[2];
    triangle.v[2] = Elm->Topology[3];

    geo_triangle_normal( geom,&triangle );

    triangle.Edge[0] = FALSE;
    triangle.Edge[1] = TRUE;
    triangle.Edge[2] = TRUE;

    if ( !geo_add_triangle( geom,&triangle ) )  return FALSE;
#endif

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        elm_4node_quad_nvalue( double *,double,double )
 *
 *     Purpose:     return shape function values at given point (u,v)
 *                 
 *                 
 *
 *     Parameters:
 *
 *         Input:   (double u,double v) point where values are evaluated
 *
 *         Output:  (double *) storage for values
 *
 *   Return value:  void
 *
 ******************************************************************************/
static void elm_4node_quad_nvalue( double *R, double u, double v )
{
     int i;

     for( i=0; i<4; i++ )
     {
          R[i] = N[i][0] + N[i][1]*u + N[i][2]*v + N[i][3]*u*v;
     }
}

/*******************************************************************************
 *
 *     Name:        elm_4node_quad_fvalue( double *,double,double )
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
double elm_4node_quad_fvalue( double *F, double u, double v )
{
     double R=0.0,uv=u*v;
     int i;

     for( i=0; i<4; i++ )
     {
          R += F[i]*(N[i][0] + N[i][1]*u + N[i][2]*v + N[i][3]*uv);
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_4node_quad_dndu_fvalue( double *,double,double )
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
static double elm_4node_quad_dndu_fvalue( double *F, double u, double v )
{
     double R=0.0;
     int i;

     for( i=0; i<4; i++ ) R += F[i]*(N[i][1] + N[i][3]*v);

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_4node_quad_dndv_fvalue( double *,double,double )
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
static double elm_4node_quad_dndv_fvalue( double *F, double u, double v )
{
     double R=0.0;
     int i;

     for( i=0; i<4; i++ ) R += F[i]*(N[i][2] + N[i][3]*u);

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_4node_quad_point_inside(
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
 *        trivial for this one...
 *
 ******************************************************************************/
int elm_4node_quad_point_inside
   (
                 double *nx, double *ny, double *nz,
        double px, double py, double pz, double *u,double *v,double *w
   )
{
    double x[3],y[3],z[3],r,s,t;

    if ( px > MAX(MAX(MAX( nx[0],nx[1] ),nx[2] ),nx[3] ) ) return FALSE;
    if ( px < MIN(MIN(MIN( nx[0],nx[1] ),nx[2] ),nx[3] ) ) return FALSE;

    if ( py > MAX(MAX(MAX( ny[0],ny[1] ),ny[2] ),ny[3] ) ) return FALSE;
    if ( py < MIN(MIN(MIN( ny[0],ny[1] ),ny[2] ),ny[3] ) ) return FALSE;

#if 0
    if ( pz > MAX(MAX(MAX( nz[0],nz[1] ),nz[2] ),nz[3] ) ) return FALSE;
    if ( pz < MIN(MIN(MIN( nz[0],nz[1] ),nz[2] ),nz[3] ) ) return FALSE;
#endif


#if 1
{
    double ax,bx,cx,dx,ay,by,cy,dy,a,b,c,d;

    ax = 0.25*( nx[0] + nx[1] + nx[2] + nx[3]);
    bx = 0.25*(-nx[0] + nx[1] + nx[2] - nx[3]);
    cx = 0.25*(-nx[0] - nx[1] + nx[2] + nx[3]);
    dx = 0.25*( nx[0] - nx[1] + nx[2] - nx[3]);

    ay = 0.25*( ny[0] + ny[1] + ny[2] + ny[3]);
    by = 0.25*(-ny[0] + ny[1] + ny[2] - ny[3]);
    cy = 0.25*(-ny[0] - ny[1] + ny[2] + ny[3]);
    dy = 0.25*( ny[0] - ny[1] + ny[2] - ny[3]);

    px -= ax;
    py -= ay;

    a = cy*dx - cx*dy;
    b = bx*cy - by*cx + dy*px - dx*py;
    c = by*px - bx*py;

    if ( ABS(a)<1.0E-16 )
    {
        r = -c/b;
        if ( r <= -1.0 || r >= 1.0 ) return FALSE;

        *v = r;
        *u = (px - cx*r)/(bx + dx*r);
        return (*u >= -1.0 && *u <= 1.0);
    }

    if ( (d=b*b - 4*a*c) < 0.0 ) return FALSE;

    r = (-b + sqrt(d))/(2*a);
    if ( r >= -1.0 && r <= 1.0 )
    {
        *v = r;
        *u = (px - cx*r)/(bx + dx*r);
        if ( *u >= -1.0 && *u <= 1.0 ) return TRUE;
    }

    r = (-b - sqrt(d))/(2*a);
    if ( r >= -1.0 && r <= 1.0 )
    {
        *v = r;
        *u = (px - cx*r)/(bx + dx*r);
        return (*u >= -1.0 && *u <= 1.0 );
    }

    return FALSE;
}
#endif


    x[0] = elm_4node_quad_fvalue( nx,0.0,0.0 ); 
    y[0] = elm_4node_quad_fvalue( ny,0.0,0.0 ); 
    z[0] = elm_4node_quad_fvalue( nz,0.0,0.0 ); 

    x[1] = nx[0];
    y[1] = ny[0];
    z[1] = nz[0];
    x[2] = nx[1];
    y[2] = ny[1];
    z[2] = nz[1];
    if ( elm_3node_triangle_point_inside( x,y,z,px,py,pz,&r,&s,&t ) )
    {
        *u = -r + s;
        *v = -r - s;
        *w = t;
        return TRUE;
    }

    x[1] = nx[1];
    y[1] = ny[1];
    z[1] = nz[1];
    x[2] = nx[2];
    y[2] = ny[2];
    z[2] = nz[2];
    if ( elm_3node_triangle_point_inside( x,y,z,px,py,pz,&r,&s,&t ) )
    {
        *u =  r + s;
        *v = -r + s;
        *w = t;
        return TRUE;
    }

    x[1] = nx[2];
    y[1] = ny[2];
    z[1] = nz[2];
    x[2] = nx[3];
    y[2] = ny[3];
    z[2] = nz[3];
    if ( elm_3node_triangle_point_inside( x,y,z,px,py,pz,&r,&s,&t ) )
    {
        *u = r - s;
        *v = r + s;
        *w = t;
        return TRUE;
    }

    x[1] = nx[3];
    y[1] = ny[3];
    z[1] = nz[3];
    x[2] = nx[0];
    y[2] = ny[0];
    z[2] = nz[0];
    if ( elm_3node_triangle_point_inside( x,y,z,px,py,pz,&r,&s,&t ) )
    {
        *u = -r - s;
        *v =  r - s;
        *w = t;
        return TRUE;
    }

    return FALSE;
}

/******************************************************************************
 *
 *     Name:        elm_4node_quad_initialize( )
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
int elm_4node_quad_initialize()
{
     element_type_t ElementDef;

     static char *Name = "ELM_4NODE_QUAD";

     ElementDef.ElementName = Name;
     ElementDef.ElementCode = 404;

     ElementDef.NumberOfNodes = 4;
     ElementDef.NodeU = NodeU;
     ElementDef.NodeV = NodeV;
     ElementDef.NodeW = NULL;

     ElementDef.PartialU = (double (*)())elm_4node_quad_dndu_fvalue;
     ElementDef.PartialV = (double (*)())elm_4node_quad_dndv_fvalue;
     ElementDef.PartialW = NULL;

     ElementDef.FunctionValue = (double (*)())elm_4node_quad_fvalue;
     ElementDef.Triangulate = (int (*)())elm_4node_quad_triangulate;
     ElementDef.PointInside = (int (*)())elm_4node_quad_point_inside;

     return elm_add_element_type( &ElementDef );
}
