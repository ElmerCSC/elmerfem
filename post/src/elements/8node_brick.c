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
 * Definitions for 8 node brick volume element.
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
 *                       Date: 27 Sep 1995
 *
 *                Modified by:
 *
 *       Date of modification:
 *
 ******************************************************************************/

#include "../elmerpost.h"
#include <elements.h>

/*
 * Eight node brick volume element.
 *
 *       7---------6
 *      /|        /|
 *     4---------5 |
 *     | |       | |
 *     | |       | |
 *     | 3-------|-2     w v
 *    w|/v       |/      |/
 *     0---------1       ---u
 *      u
 */

static double A[8][8],N[8][8];

static double NodeU[] = { -1.0,  1.0,  1.0, -1.0, -1.0,  1.0, 1.0, -1.0 };
static double NodeV[] = { -1.0, -1.0,  1.0,  1.0, -1.0, -1.0, 1.0,  1.0 };
static double NodeW[] = { -1.0, -1.0, -1.0, -1.0,  1.0,  1.0, 1.0,  1.0 };


/*******************************************************************************
 *
 *     Name:        elm_8node_brick_triangulate( geometry_t *,element_t * )
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
static int elm_8node_brick_triangulate( geometry_t *geom,element_t *brick )
{
    element_t quad;
    int i,j;

    if ( GlobalOptions.VolumeSides )
    {
      int topo[4];

      quad.DisplayFlag = TRUE;
      quad.Topology = topo;
      for( i=0; i<MAX_GROUP_IDS; i++ ) quad.GroupIds[i] = brick->GroupIds[i];

      for( i=0; i<6; i++ )
      {
          for( j=0; j<4; j++ ) quad.Topology[j] = brick->Topology[ElmBrickFace[i][j]];
          if ( !elm_4node_quad_triangulate( geom, &quad, brick ) ) return FALSE;
      }
    } else {
      if ( !geo_add_edge( geom, brick->Topology[0], brick->Topology[1],brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[0], brick->Topology[3],brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[0], brick->Topology[4],brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[1], brick->Topology[2],brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[1], brick->Topology[5],brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[2], brick->Topology[3],brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[2], brick->Topology[6],brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[3], brick->Topology[7],brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[4], brick->Topology[5],brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[4], brick->Topology[7],brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[5], brick->Topology[6],brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[6], brick->Topology[7],brick ) ) return FALSE;
    }

    return TRUE;
}




/*******************************************************************************
 *
 *     Name:        elm_8node_brick_shape_functions( )
 *
 *     Purpose:     Initialize element shape function array. Internal only.
 *
 *     Parameters:
 *
 *         Input:   Global (filewise) variables NodeU,NodeV,NodeW
 *
 *         Output:  Global (filewise) variable N[8][8], will contain
 *                  shape function coefficients
 *
 *   Return value:  void
 *
 ******************************************************************************/
static void elm_8node_brick_shape_functions()
{
     double u,v,w;

     int i,j;

     for( i=0; i<8; i++ )
     {
         u = NodeU[i];
         v = NodeV[i];
         w = NodeW[i];

         A[i][0]  = 1;
         A[i][1]  = u;
         A[i][2]  = v;
         A[i][3]  = w;
         A[i][4]  = u*v;
         A[i][5]  = u*w;
         A[i][6]  = v*w;
         A[i][7]  = u*v*w;
     }

     lu_mtrinv( (double *)A,8 );

     for( i=0; i<8; i++ )
        for( j=0; j<8; j++ ) N[i][j] = A[j][i];
}

/*******************************************************************************
 *
 *     Name:        elm_8node_brick_fvalue( double *,double,double,double )
 *
 *     Purpose:     return value of a quantity given on nodes at point (u,v,w)
 *                  Use trough (element_type_t *) structure.
 *
 *     Parameters:
 *
 *         Input:  (double *) quantity values at nodes 
 *                 (double,double,double) point where values are evaluated
 *
 *         Output:  none
 *
 *   Return value:  quantity value
 *
 ******************************************************************************/
static double elm_8node_brick_fvalue(double *F,double u,double v,double w)
{
     double R=0.0,uv=u*v,uw=u*w,vw=v*w,uvw=u*v*w;
     int i;

     for( i=0; i<8; i++ )
     {
         R += F[i]*(N[i][0]+
                    N[i][1]*u+
                    N[i][2]*v+
                    N[i][3]*w+
                    N[i][4]*uv+
                    N[i][5]*uw+
                    N[i][6]*vw+
                    N[i][7]*uvw );
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_8node_brick_dndu_fvalue( double *,double,double,double )
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
static double elm_8node_brick_dndu_fvalue(double *F,double u,double v,double w)
{
     double R=0.0,vw=v*w;
     int i;

     for( i=0; i<8; i++ )
     {
         R += F[i]*(N[i][1] + N[i][4]*v + N[i][5]*w + N[i][7]*vw);
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_8node_brick_dndv_fvalue( double *,double,double,double )
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
static double elm_8node_brick_dndv_fvalue(double *F,double u,double v,double w)
{
     double R=0.0,uw=u*w;
     int i;

     for( i=0; i<8; i++ )
     {
         R += F[i]*(N[i][2] + N[i][4]*u + N[i][6]*w + N[i][7]*uw);
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_8node_brick_dndw_fvalue( double *,double,double,double )
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
static double elm_8node_brick_dndw_fvalue(double *F,double u,double v,double w)
{
     double R=0.0,uv=u*v;
     int i;

     for( i=0; i<8; i++ )
     {
         R += F[i]*(N[i][3] + N[i][5]*u + N[i][6]*v + N[i][7]*uv);
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_8node_brick_point_inside(
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
int elm_8node_brick_point_inside
   (
                    double *nx, double *ny, double *nz,
        double px, double py, double pz, double *u,double *v,double *w
   )
{
    double x[4],y[4],z[4],r,s,t,maxx,minx,maxy,miny,maxz,minz;
    int i,j;

    static int map[12][3] =
    { 
        { 0,1,2 }, { 0,2,3 }, { 4,5,6 }, { 4,6,7 }, { 3,2,6 }, { 3,6,7 },
        { 1,5,6 }, { 1,6,2 }, { 0,4,7 }, { 0,7,3 }, { 0,1,5 }, { 0,5,4 },
    };

    maxx = minx = nx[0];
    maxy = miny = ny[0];
    maxz = minz = nz[0];

    for( i=1; i<8; i++ )
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

    x[0] = 0.125*(nx[0]+nx[1]+nx[2]+nx[3]+nx[4]+nx[5]+nx[6]+nx[7]);
    y[0] = 0.125*(ny[0]+ny[1]+ny[2]+ny[3]+ny[4]+ny[5]+ny[6]+ny[7]);
    z[0] = 0.125*(nz[0]+nz[1]+nz[2]+nz[3]+nz[4]+nz[5]+nz[6]+nz[7]);
    
    for( i=0; i<12; i++ )
    {
        for( j=0; j<3; j++ )
        {
            x[j+1] = nx[map[i][j]];
            y[j+1] = ny[map[i][j]];
            z[j+1] = nz[map[i][j]];
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
 *     Name:        elm_8node_brick_isoline
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
int elm_8node_brick_isoline
  (
     double K, double *F, double *C,double *X,double *Y,double *Z,line_t *Line
  )
{
    double f[4],c[4],x[4],y[4],z[4];

    int i, j, k, n=0, above=0;

    for( i=0; i<8; i++ ) above += F[i]>K;
    if ( above == 0 || above == 8 ) return 0;

    for( i=0; i<6; i++ )
    {
        for( j=0; j<4; j++ )
        {
            k = ElmBrickFace[i][j];
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

/*******************************************************************************
 *
 *     Name:        elm_8node_brick_isosurface
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
int elm_8node_brick_isosurface
   (
       double  K,double *F,double *C, double *X,double *Y,double *Z,
            double *U,double *V,double *W, polygon_t *Polygon
   )
{
    double tx[4],ty[4],tz[4],tu[4],tv[4],tw[4],tf[4],tc[4];
    int i,j,n;

    int above = 0, below = 0;

    for( i=0; i<8; i++ ) above += F[i]>K;
    for( i=0; i<8; i++ ) below += F[i]<K;
    if ( below == 8 || above == 8 ) return 0;

    tx[0] = 0.125 * ( X[0] + X[1] + X[2] + X[3] + X[4] + X[5] + X[6] + X[7] );
    ty[0] = 0.125 * ( Y[0] + Y[1] + Y[2] + Y[3] + Y[4] + Y[5] + Y[6] + Y[7] );
    tz[0] = 0.125 * ( Z[0] + Z[1] + Z[2] + Z[3] + Z[4] + Z[5] + Z[6] + Z[7] );
    tu[0] = 0.125 * ( U[0] + U[1] + U[2] + U[3] + U[4] + U[5] + U[6] + U[7] );
    tv[0] = 0.125 * ( V[0] + V[1] + V[2] + V[3] + V[4] + V[5] + V[6] + V[7] );
    tw[0] = 0.125 * ( W[0] + W[1] + W[2] + W[3] + W[4] + W[5] + W[6] + W[7] );
    tf[0] = 0.125 * ( F[0] + F[1] + F[2] + F[3] + F[4] + F[5] + F[6] + F[7] );
    tc[0] = 0.125 * ( C[0] + C[1] + C[2] + C[3] + C[4] + C[5] + C[6] + C[7] );
    
    n = 0;
    for( i=0; i<6; i++ )
    {
        tx[1] = 0.0;
        ty[1] = 0.0;
        tz[1] = 0.0;
        tu[1] = 0.0;
        tv[1] = 0.0;
        tw[1] = 0.0;
        tf[1] = 0.0;
        tc[1] = 0.0;
        for( j=0; j<4; j++ )
        {
            tx[1] += X[ElmBrickFace[i][j]];
            ty[1] += Y[ElmBrickFace[i][j]];
            tz[1] += Z[ElmBrickFace[i][j]];
            tu[1] += U[ElmBrickFace[i][j]];
            tv[1] += V[ElmBrickFace[i][j]];
            tw[1] += W[ElmBrickFace[i][j]];
            tf[1] += F[ElmBrickFace[i][j]];
            tc[1] += C[ElmBrickFace[i][j]];
        }
        tx[1] /= 4.0;
        ty[1] /= 4.0;
        tz[1] /= 4.0;
        tu[1] /= 4.0;
        tv[1] /= 4.0;
        tw[1] /= 4.0;
        tf[1] /= 4.0;
        tc[1] /= 4.0;
        
        for( j=0; j<3; j++ )
        {
            tx[2] = X[ElmBrickFace[i][j]];
            ty[2] = Y[ElmBrickFace[i][j]];
            tz[2] = Z[ElmBrickFace[i][j]];
            tu[2] = U[ElmBrickFace[i][j]];
            tv[2] = V[ElmBrickFace[i][j]];
            tw[2] = W[ElmBrickFace[i][j]];
            tf[2] = F[ElmBrickFace[i][j]];
            tc[2] = C[ElmBrickFace[i][j]];

            tx[3] = X[ElmBrickFace[i][j+1]];
            ty[3] = Y[ElmBrickFace[i][j+1]];
            tz[3] = Z[ElmBrickFace[i][j+1]];
            tu[3] = U[ElmBrickFace[i][j+1]];
            tv[3] = V[ElmBrickFace[i][j+1]];
            tw[3] = W[ElmBrickFace[i][j+1]];
            tf[3] = F[ElmBrickFace[i][j+1]];
            tc[3] = C[ElmBrickFace[i][j+1]];
            n += elm_4node_tetra_isosurface( K,tf,tc,tx,ty,tz,tu,tv,tw,&Polygon[n] );
        }

        tx[2] = X[ElmBrickFace[i][3]];
        ty[2] = Y[ElmBrickFace[i][3]];
        tz[2] = Z[ElmBrickFace[i][3]];
        tu[2] = U[ElmBrickFace[i][3]];
        tv[2] = V[ElmBrickFace[i][3]];
        tw[2] = W[ElmBrickFace[i][3]];
        tf[2] = F[ElmBrickFace[i][3]];
        tc[2] = C[ElmBrickFace[i][3]];

        tx[3] = X[ElmBrickFace[i][0]];
        ty[3] = Y[ElmBrickFace[i][0]];
        tz[3] = Z[ElmBrickFace[i][0]];
        tu[3] = U[ElmBrickFace[i][0]];
        tv[3] = V[ElmBrickFace[i][0]];
        tw[3] = W[ElmBrickFace[i][0]];
        tf[3] = F[ElmBrickFace[i][0]];
        tc[3] = C[ElmBrickFace[i][0]];
        n += elm_4node_tetra_isosurface( K,tf,tc,tx,ty,tz,tu,tv,tw,&Polygon[n] );
    }
    return n;
}

/*******************************************************************************
 *
 *     Name:        elm_8node_brick_isosurface1
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
int elm_8node_brick_isosurface1
   (
       double  K,double *F,double *C, double *X,double *Y,double *Z,
            double *U,double *V,double *W, polygon_t *Polygon
   )
{
    double tx[4],ty[4],tz[4],tu[4],tv[4],tw[4],tf[4],tc[4];
    int i,j,l,n;

    static int map[12][3] =
    {
        { 0,1,2 }, { 0,2,3 }, { 4,5,6 }, { 4,6,7 }, { 3,2,6 }, { 3,6,7 },
        { 1,5,6 }, { 1,6,2 }, { 0,4,7 }, { 0,7,3 }, { 0,1,5 }, { 0,5,4 },
    };

    int above = 0, below = 0;

    for( i=0; i<8; i++ ) above += F[i]>K;
    for( i=0; i<8; i++ ) below += F[i]<K;
    if ( below == 8 || above == 8 ) return 0;

    tx[0] = 0.125 * ( X[0] + X[1] + X[2] + X[3] + X[4] + X[5] + X[6] + X[7] );
    ty[0] = 0.125 * ( Y[0] + Y[1] + Y[2] + Y[3] + Y[4] + Y[5] + Y[6] + Y[7] );
    tz[0] = 0.125 * ( Z[0] + Z[1] + Z[2] + Z[3] + Z[4] + Z[5] + Z[6] + Z[7] );
    tu[0] = 0.125 * ( U[0] + U[1] + U[2] + U[3] + U[4] + U[5] + U[6] + U[7] );
    tv[0] = 0.125 * ( V[0] + V[1] + V[2] + V[3] + V[4] + V[5] + V[6] + V[7] );
    tw[0] = 0.125 * ( W[0] + W[1] + W[2] + W[3] + W[4] + W[5] + W[6] + W[7] );
    tf[0] = 0.125 * ( F[0] + F[1] + F[2] + F[3] + F[4] + F[5] + F[6] + F[7] );
    tc[0] = 0.125 * ( C[0] + C[1] + C[2] + C[3] + C[4] + C[5] + C[6] + C[7] );

    n = 0;
    for( i=0; i<12; i++ )
    {
        for( j=1; j<4; j++ )
        {
            l = map[i][j-1];
            tx[j] = X[l];
            ty[j] = Y[l];
            tz[j] = Z[l];
            tu[j] = U[l];
            tv[j] = V[l];
            tw[j] = W[l];
            tf[j] = F[l];
            tc[j] = C[l];
        }
        n += elm_4node_tetra_isosurface( K,tf,tc,tx,ty,tz,tu,tv,tw,&Polygon[n] );
    }
    return n;
}


/******************************************************************************
 *
 *     Name:        elm_8node_brick_initialize( )
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
int elm_8node_brick_initialize()
{
     static char *Name = "ELM_8NODE_BRICK";

     element_type_t ElementDef;

     elm_8node_brick_shape_functions();

     ElementDef.ElementName = Name;
     ElementDef.ElementCode = 808;

     ElementDef.NumberOfNodes = 8;

     ElementDef.NodeU = NodeU;
     ElementDef.NodeV = NodeV;
     ElementDef.NodeW = NodeW;

     ElementDef.PartialU = (double (*)())elm_8node_brick_dndu_fvalue;
     ElementDef.PartialV = (double (*)())elm_8node_brick_dndv_fvalue;
     ElementDef.PartialW = (double (*)())elm_8node_brick_dndw_fvalue;

     ElementDef.FunctionValue = (double (*)())elm_8node_brick_fvalue;
     ElementDef.Triangulate   = (int (*)())elm_8node_brick_triangulate;
     ElementDef.PointInside   = (int (*)())elm_8node_brick_point_inside;
     ElementDef.IsoLine       = (int (*)())elm_8node_brick_isoline;
     ElementDef.IsoSurface    = (int (*)())elm_8node_brick_isosurface1;

     return elm_add_element_type( &ElementDef ) ;
}

