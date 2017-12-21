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
 * Definition of 6 node wedge element.
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
 *                       Date: 4 Oct 1995
 *
 ******************************************************************************/

#include "../elmerpost.h"
#include <elements.h>

/*
 * Six node wedge volume element.
 *
 *                  5-----------4    
 *                 /|          /|          w  v
 *                / |        /  |          | /
 *               /  |      /    |          |/
 *              /   2----/------1          ---u
 *             /   /    /     //
 *            /   /  ///   //
 *           /    //     //
 *          /  //     //
 *         3/  /   //
 *         |  /   //
 *         | / //
 *         |//
 *         0
 */

static double A[6][6],N[6][6];

static double NodeU[] = {  0.0, 1.0, 0.0, 0.0, 1.0, 0.0 };
static double NodeV[] = {  0.0, 0.0, 1.0, 0.0, 0.0, 1.0 };
static double NodeW[] = { -1.0,-1.0,-1.0, 1.0, 1.0, 1.0 };

/*******************************************************************************
 *
 *     Name:        elm_6node_wedge_triangulate
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
static int elm_6node_wedge_triangulate( geometry_t *geom,element_t *wedge )
{
    element_t element;
    int i,j;


    if ( GlobalOptions.VolumeSides )
    {
       int topo[4];

       element.DisplayFlag = TRUE;
       element.Topology = topo;
       for( i=0; i<MAX_GROUP_IDS; i++ ) element.GroupIds[i] = wedge->GroupIds[i];
   
       for( i=0; i<3; i++ )
       {
           for( j=0; j<4; j++ )
           {
               element.Topology[j] = wedge->Topology[ElmWedgeFace[i][j]];
           }
           if ( !elm_4node_quad_triangulate( geom, &element, wedge ) ) return FALSE;
       }

       for( i=3; i<5; i++ )
       {
           for( j=0; j<3; j++ )
           {
               element.Topology[j] = wedge->Topology[ElmWedgeFace[i][j]];
           }
           if ( !elm_3node_triangle_triangulate( geom, &element, wedge ) ) return FALSE;
       }
    } else {
       if ( !geo_add_edge( geom, wedge->Topology[0], wedge->Topology[1],wedge ) ) return FALSE;
       if ( !geo_add_edge( geom, wedge->Topology[0], wedge->Topology[2],wedge ) ) return FALSE;
       if ( !geo_add_edge( geom, wedge->Topology[0], wedge->Topology[3],wedge ) ) return FALSE;
       if ( !geo_add_edge( geom, wedge->Topology[1], wedge->Topology[2],wedge ) ) return FALSE;
       if ( !geo_add_edge( geom, wedge->Topology[1], wedge->Topology[4],wedge ) ) return FALSE;
       if ( !geo_add_edge( geom, wedge->Topology[2], wedge->Topology[5],wedge ) ) return FALSE;
       if ( !geo_add_edge( geom, wedge->Topology[3], wedge->Topology[4],wedge ) ) return FALSE;
       if ( !geo_add_edge( geom, wedge->Topology[3], wedge->Topology[5],wedge ) ) return FALSE;
       if ( !geo_add_edge( geom, wedge->Topology[4], wedge->Topology[5],wedge ) ) return FALSE;
    }

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        elm_6node_wedge_shape_functions
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
 *   Return value:  void
 *
 ******************************************************************************/
static void elm_6node_wedge_shape_functions()
{
     double u,v,w;
     int i,j;

     for( i=0; i<6; i++ )
     {
         u = NodeU[i];
         v = NodeV[i];
         w = NodeW[i];

         A[i][0]   = 1;
         A[i][1]   = u;
         A[i][2]   = v;
         A[i][3]   = w;
         A[i][4]   = u*w;
         A[i][5]   = v*w;
     }

     lu_mtrinv( (double *)A,6 );

     for( i=0; i<6; i++ )
        for( j=0; j<6; j++ ) N[i][j] = A[j][i];
}

/*******************************************************************************
 *
 *     Name:        elm_6node_wedge_fvalue
 *
 *     Purpose:     return value of a quantity given on nodes at point (u,v)
 *                  Use trough (element_type_t *) structure.
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
static double elm_6node_wedge_fvalue(double *F,double u,double v,double w)
{
     double R=0.0,uw=u*w,vw=v*w;
     int i;

     for( i=0; i<6; i++ )
     {
         R += F[i]*(N[i][0]+
                    N[i][1]*u+
                    N[i][2]*v+
                    N[i][3]*w+
                    N[i][4]*uw+
                    N[i][5]*vw);
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_6node_wedge_dndu_fvalue
 *
 *     Purpose:     return value of a first partial derivate in (v) of a
 *                  quantity given on nodes at point (u,v).
 *                  Use trough (element_type_t *) structure.
 *                 
 *
 *     Parameters:
 *
 *         Input:  (double *) quantity values at nodes 
 *                 (double u,double v,double w) point where values are evaluated
 *
 *         Output:  none
 *
 *   Return value:  quantity value
 *
 ******************************************************************************/
static double elm_6node_wedge_dndu_fvalue(double *F,double u,double v,double w)
{
     double R=0.0;
     int i;

     for( i=0; i<6; i++ )
     {
         R += F[i]*( N[i][1] + N[i][4]*w );
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_6node_wedge_dndv_fvalue
 *
 *     Purpose:     return value of a first partial derivate in (v) of a
 *                  quantity given on nodes at point (u,v,w).
 *                  Use trough (element_type_t *) structure.
 *                 
 *
 *     Parameters:
 *
 *         Input:  (double *) quantity values at nodes 
 *                 (double u,double v,double w) point where values are evaluated
 *
 *         Output:  none
 *
 *   Return value:  quantity value
 *
 ******************************************************************************/
static double elm_6node_wedge_dndv_fvalue(double *F,double u,double v,double w)
{
     double R=0.0;
     int i;

     for( i=0; i<6; i++ )
     {
         R += F[i]*( N[i][2] + N[i][5]*w );
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_6node_wedge_dndw_fvalue
 *
 *     Purpose:     return value of a first partial derivate in (w) of a
 *                  quantity given on nodes at point (u,v,w)
 *                  Use trough (element_type_t *) structure.
 *                 
 *
 *     Parameters:
 *
 *         Input:  (double *) quantity values at nodes 
 *                 (double u,double v,double w) point where values are evaluated
 *
 *         Output:  none
 *
 *   Return value:  quantity value
 *
 ******************************************************************************/
static double elm_6node_wedge_dndw_fvalue(double *F,double u,double v,double w)
{
     double R=0.0;
     int i;

     for( i=0; i<6; i++ )
     {
         R += F[i]*( N[i][3] + N[i][4]*u + N[i][5]*v );
     }

     return R;
}



/*******************************************************************************
 *
 *     Name:        elm_6node_wedge_point_inside(
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
int elm_6node_wedge_point_inside
   (
                    double *nx, double *ny, double *nz,
        double px, double py, double pz, double *u,double *v,double *w
   )
{
    double x[4],y[4],z[4],r,s,t,maxx,minx,maxy,miny,maxz,minz;
    int i,j;

    static int map[3][4] =
    { 
        { 4,3,2,0 }, { 4,2,1,0 }, { 4,5,3,2 } 
    };

    maxx = minx = nx[0];
    maxy = miny = ny[0];
    maxz = minz = nz[0];

    for( i=1; i<6; i++ )
    {
        maxx = MAX( nx[i],maxx );
        maxy = MAX( ny[i],maxy );
        maxz = MAX( nz[i],maxz );

        minx = MIN( nx[i],minx );
        miny = MIN( ny[i],miny );
        minz = MIN( nz[i],minz );
    }

    if ( px > maxx || px < minx ) return FALSE;
    if ( py > maxy || py < miny ) return FALSE;
    if ( pz > maxz || pz < minz ) return FALSE;

    for( i=0; i<3; i++ )
    {
        for( j=0; j<4; j++ )
        {
            x[j] = nx[map[i][j]];
            y[j] = ny[map[i][j]];
            z[j] = nz[map[i][j]];
        }

        if ( elm_4node_tetra_point_inside( x,y,z,px,py,pz,&r,&s,&t ) )
        {
            *u = NodeU[map[i][0]]*r + NodeU[map[i][1]]*s + NodeU[map[i][2]]*t;
            *v = NodeV[map[i][0]]*r + NodeV[map[i][1]]*s + NodeV[map[i][2]]*t;
            *w = NodeW[map[i][0]]*r + NodeW[map[i][1]]*s + NodeW[map[i][2]]*t;
            return TRUE;
        }
    }

    return FALSE;
}

/*******************************************************************************
 *
 *     Name:        elm_6node_wedge_isoline
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
 *   Return value:  number of lines generated (0,...,24)
 *
 ******************************************************************************/
int elm_6node_wedge_isoline
  (
     double K, double *F, double *C,double *X,double *Y,double *Z,line_t *Line
  )
{
    double f[4],c[4],x[4],y[4],z[4];

    int i, j, k, n=0, above=0;

    for( i=0; i<6; i++ ) above += F[i]>K;
    if ( above == 0 || above == 6 ) return 0;

    for( i=0; i<3; i++ )
    {
        for( j=0; j<4; j++ )
        {
            k = ElmWedgeFace[i][j];
            f[j] = F[k];
            c[j] = C[k];
            x[j] = X[k];
            y[j] = Y[k];
            z[j] = Z[k];
        }
        n += elm_4node_quad_isoline( K,f,c,x,y,z,&Line[n] );
    }

    for( i=3; i<5; i++ )
    {
        for( j=0; j<3; j++ )
        {
            k = ElmWedgeFace[i][j];
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
 *     Name:        elm_6node_wedge_isosurface
 *
 *     Purpose:     Extract isosurfaces for brick element.
 *
 *     Parameters:
 *
 *         Input:   (double )  K: contour threshold
 *                  (double *) F: contour quantity values at nodes
 *                  (double *) C: color quantity values at nodes
 *                  (double *) X,Y,Z: node coordinates
 *                  (double *) U,V,W: normal vector at nodes
 *
 *         Output:  (polygon_t *)Polygon: output triangles.
 *
 *   Return value:  How many triangles we've got (possible values are 0-48)...
 *
 ******************************************************************************/
int elm_6node_wedge_isosurface
   (
       double  K,double *F,double *C, double *X,double *Y,double *Z,
            double *U,double *V,double *W, polygon_t *Polygon
   )
{
    double tx[4],ty[4],tz[4],tu[4],tv[4],tw[4],tf[4],tc[4];
    int i,j,n;

    int above = 0, below = 0;

     static int map[3][4] =
     { 
         { 4,3,2,0 }, { 4,2,1,0 }, { 4,5,3,2 } 
     };

    for( i=0; i<6; i++ ) above += F[i]>K;
    for( i=0; i<6; i++ ) below += F[i]<K;
    if ( below == 6 || above == 6 ) return 0;

    n = 0;
    for( i=0; i<3; i++ )
    {
        for( j=0; j<4; j++ )
        {
            tx[j] = X[map[i][j]];
            ty[j] = Y[map[i][j]];
            tz[j] = Z[map[i][j]];
            tu[j] = U[map[i][j]];
            tv[j] = V[map[i][j]];
            tw[j] = W[map[i][j]];
            tf[j] = F[map[i][j]];
            tc[j] = C[map[i][j]];
        }
        n += elm_4node_tetra_isosurface( K,tf,tc,tx,ty,tz,tu,tv,tw,&Polygon[n] );
    }
    return n;
}


/******************************************************************************
 *
 *     Name:        elm_6node_wedge_initialize
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
int elm_6node_wedge_initialize()
{
     static char *Name = "ELM_6NODE_WEDGE";

     element_type_t ElementDef;

     elm_6node_wedge_shape_functions();

     ElementDef.ElementName = Name;
     ElementDef.ElementCode = 706;

     ElementDef.NumberOfNodes = 6;

     ElementDef.NodeU = NodeU;
     ElementDef.NodeV = NodeV;
     ElementDef.NodeW = NodeW;

     ElementDef.PartialU = (double (*)())elm_6node_wedge_dndu_fvalue;
     ElementDef.PartialV = (double (*)())elm_6node_wedge_dndv_fvalue;
     ElementDef.PartialW = (double (*)())elm_6node_wedge_dndw_fvalue;

     ElementDef.FunctionValue = (double (*)())elm_6node_wedge_fvalue;
     ElementDef.Triangulate   = (int (*)())elm_6node_wedge_triangulate;
     ElementDef.PointInside   = (int (*)())elm_6node_wedge_point_inside;
     ElementDef.IsoLine       = (int (*)())elm_6node_wedge_isoline;
     ElementDef.IsoSurface    = (int (*)())elm_6node_wedge_isosurface;

     return elm_add_element_type( &ElementDef ) ;
}
