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
 * Isoparametric 20 node volume element.
 *
 *                             7------18-------6
 *                           19|             17|
 *                           / |             / |
 *                          4------16-------5  |
 *                          | 15            |  14
 *                          |  |            |  |
 *                         12  |           13  |
 *                          |  3-----10-----|--2
 *                          | 11            | 9        w v
 *                         w|/v             |/         |/
 *                          0-------8-------1           ---u
 *                           u
 *
 * shape functions:  N  = (1+u u)(1+v v)(1+w w)(u u+v v+w w-2)/8   i = 0..7
 *                    i       i      i      i    i   i   i
 *                   N  = (1-u^2)(1+v v)(1+w w)/4     i=8,10,16,18  (u = 0)
 *                    i              i      i                         i
 *                   N  = (1+u u)(1-v^2)(1+w w)/4     i=9,11,17,19  (v = 0)
 *                    i       i             i                         i
 *                   N  = (1+u u)(1+v v)(1-w^2)/4     i=12,13,14,15 (w = 0)
 *                    i       i      i                                i
 *
 * @/@u              N  =  u (1+v v)(1+w w)(2u u+v v+w w-1)/8      i = 0..7
 *                    i     i    i      i     i   i   i
 * @/@u              N  = -2u(1+v v)(1+w w)/4         i=8,10,16,18  (u = 0)
 *                    i          i      i                             i
 * @/@u              N  =  u (1-v^2)(1+w w)/4         i=9,11,17,19  (v = 0)
 *                    i     i           i                             i 
 * @/@u              N  =  u (1+v v)(1-w^2)/4         i=12,13,14,15 (w = 0)
 *                    i     i    i                                    i
 *
 */

static double A[20][20],N[20][20];

static double NodeU[] = {
     -1.0,  1.0,  1.0, -1.0, -1.0,  1.0, 1.0, -1.0, 0.0,  1.0,
      0.0, -1.0, -1.0,  1.0,  1.0, -1.0, 0.0,  1.0, 0.0, -1.0
   };

static double NodeV[] = {
     -1.0, -1.0,  1.0,  1.0, -1.0, -1.0,  1.0, 1.0, -1.0, 0.0,
      1.0,  0.0, -1.0, -1.0,  1.0,  1.0, -1.0, 0.0,  1.0, 0.0
   };

static double NodeW[] = {
     -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0,
     -1.0, -1.0,  0.0,  0.0, 0.0, 0.0, 1.0, 1.0,  1.0,  1.0
   };


/*******************************************************************************
 *
 *     Name:        elm_20node_brick_triangulate( geometry_t *,element_t * )
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
static int elm_20node_brick_triangulate( geometry_t *geom,element_t *brick )
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
          for( j=0; j<4; j++ )
          {
            quad.Topology[j] = brick->Topology[ElmBrickFace[i][j]];
          }
          if ( !elm_4node_quad_triangulate( geom, &quad, brick ) ) return FALSE;
      }
    } else {
      if ( !geo_add_edge( geom, brick->Topology[0],  brick->Topology[8],  brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[8],  brick->Topology[1],  brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[1],  brick->Topology[9],  brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[9],  brick->Topology[2],  brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[2],  brick->Topology[10], brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[10], brick->Topology[3],  brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[3],  brick->Topology[11], brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[11], brick->Topology[0],  brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[0],  brick->Topology[12], brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[12], brick->Topology[4],  brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[1],  brick->Topology[13], brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[13], brick->Topology[5],  brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[2],  brick->Topology[14], brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[14], brick->Topology[6],  brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[3],  brick->Topology[15], brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[15], brick->Topology[7],  brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[15], brick->Topology[7],  brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[4],  brick->Topology[16], brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[16], brick->Topology[5],  brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[5],  brick->Topology[17], brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[17], brick->Topology[6],  brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[6],  brick->Topology[18], brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[18], brick->Topology[7],  brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[7],  brick->Topology[19], brick ) ) return FALSE;
      if ( !geo_add_edge( geom, brick->Topology[19], brick->Topology[4],  brick ) ) return FALSE;
    }

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        elm_20node_brick_shape_functions( )
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
static void elm_20node_brick_shape_functions()
{
     double u,v,w;
     int i,j;

     for( i=0; i<20; i++ )
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
         A[i][8]  = u*u;
         A[i][9]  = v*v;
         A[i][10] = w*w;
         A[i][11] = u*u*v;
         A[i][12] = u*u*w;
         A[i][13] = u*v*v;
         A[i][14] = u*w*w;
         A[i][15] = v*v*w;
         A[i][16] = v*w*w;
         A[i][17] = u*u*v*w;
         A[i][18] = u*v*v*w;
         A[i][19] = u*v*w*w;
     }

     lu_mtrinv( (double *)A,20 );

     for( i=0; i<20; i++ )
        for( j=0; j<20; j++ ) N[i][j] = A[j][i];
}

/*******************************************************************************
 *
 *     Name:        elm_20node_brick_fvalue( double *,double,double,double )
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
static double elm_20node_brick_fvalue(double *F,double u,double v,double w)
{
     double R=0.0,uv=u*v,uw=u*w,vw=v*w,uvw=u*v*w,uu=u*u,vv=v*v,ww=w*w,
                  uuv=uu*v,uuw=uu*w,uvv=u*vv,uww=u*ww,vvw=vv*w,vww=v*ww,
                  uuvw=uu*vw,uvvw=vv*uw,uvww=uv*ww;
     int i;

     for( i=0; i<20; i++ )
     {
         R += F[i]*( N[i][0]       +
                     N[i][1]*u     +
                     N[i][2]*v     +
                     N[i][3]*w     +
                     N[i][4]*uv    +
                     N[i][5]*uw    +
                     N[i][6]*vw    +
                     N[i][7]*uvw   +
                     N[i][8]*uu    +
                     N[i][9]*vv    +
                     N[i][10]*ww   +
                     N[i][11]*uuv  +
                     N[i][12]*uuw  +
                     N[i][13]*uvv  +
                     N[i][14]*uww  +
                     N[i][15]*vvw  +
                     N[i][16]*vww  +
                     N[i][17]*uuvw +
                     N[i][18]*uvvw +
                     N[i][19]*uvww );
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_20node_brick_dndu_fvalue( double *,double,double,double )
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
static double elm_20node_brick_dndu_fvalue(double *F,double u,double v,double w)
{
     double R=0.0,vw=v*w,u2=2*u,u2v=u2*v,u2w=u2*w,vv=v*v,ww=w*w,u2vw=u2*vw,vvw=vv*w,vww=v*ww;
     int i;

     for( i=0; i<20; i++ )
     {
         R += F[i]*( N[i][1]       +
                     N[i][4]*v     +
                     N[i][5]*w     +
                     N[i][7]*vw    +
                     N[i][8]*u2    +
                     N[i][11]*u2v  +
                     N[i][12]*u2w  +
                     N[i][13]*vv   +
                     N[i][14]*ww   +
                     N[i][17]*u2vw +
                     N[i][18]*vvw  +
                     N[i][19]*vww  );
     }

     return R;
}
/*******************************************************************************
 *
 *     Name:        elm_20node_brick_dndv_fvalue( double *,double,double,double )
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
static double elm_20node_brick_dndv_fvalue(double *F,double u,double v,double w)
{
     double R=0.0,uw=u*w,v2=2*v,uu=u*u,uv2=u*v2,v2w=v2*w,ww=w*w,uuw=uu*w,uv2w=uw*v2,uww=u*ww;
     int i;

     for( i=0; i<20; i++ )
     {
         R += F[i]*( N[i][2]       +
                     N[i][4]*u     +
                     N[i][6]*w     +
                     N[i][7]*uw    +
                     N[i][9]*v2    +
                     N[i][11]*uu   +
                     N[i][13]*uv2  +
                     N[i][15]*v2w  +
                     N[i][16]*ww   +
                     N[i][17]*uuw  +
                     N[i][18]*uv2w +
                     N[i][19]*uww  );
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_20node_brick_dndw_fvalue( double *,double,double,double )
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
static double elm_20node_brick_dndw_fvalue(double *F,double u,double v,double w)
{
     double R=0.0,uv=u*v,w2=2*w,uu=u*u,uw2=u*w2,vv=v*v,vw2=v*w2,uuv=uu*v,uvv=u*vv,uvw2=uv*w2;
     int i;

     for( i=0; i<20; i++ )
     {
         R += F[i]*( N[i][3]       +
                     N[i][5]*u     +
                     N[i][6]*v     +
                     N[i][7]*uv    +
                     N[i][10]*w2   +
                     N[i][12]*uu   +
                     N[i][14]*uw2  +
                     N[i][15]*vv   +
                     N[i][16]*vw2  +
                     N[i][17]*uuv  +
                     N[i][18]*uvv  +
                     N[i][19]*uvw2 );
     }

     return R;
}


/*******************************************************************************
 *
 *     Name:        elm_20node_brick_point_inside(
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
 * NOTES: TODO: FIX THIS
 *        trivial for this one...
 *
 ******************************************************************************/
int elm_20node_brick_point_inside
   (
                 double *nx, double *ny, double *nz,
        double px, double py, double pz, double *u,double *v,double *w
   )
{
    return elm_8node_brick_point_inside( nx,ny,nz,px,py,pz,u,v,w );
}

/*******************************************************************************
 *
 *     Name:        elm_20node_brick_isoline
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
 *   Return value:  number of lines generated (0,...,144)
 *
 ******************************************************************************/
static int elm_20node_brick_isoline
  (
     double K, double *F, double *C,double *X,double *Y,double *Z,line_t *Line
  )
{
    double f[8],c[8],x[8],y[8],z[8];

    int i, j, k, n=0, above=0;

    for( i=0; i<20; i++ ) above += F[i]>K;
    if ( above == 0 || above == 20 ) return 0;

    for( i=0; i<6; i++ )
    {
        for( j=0; j<8; j++ )
        {
            k = ElmBrickFace[i][j];
            f[j] = F[k];
            c[j] = C[k];
            x[j] = X[k];
            y[j] = Y[k];
            z[j] = Z[k];
        }
        n += elm_8node_quad_isoline( K,f,c,x,y,z,&Line[n] );
    }

    return n;
}


/*******************************************************************************
 *
 *     Name:        elm_20node_brick_isosurface(
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
 * NOTES: TODO: FIX THIS
 *
 ******************************************************************************/
int elm_20node_brick_isosurface
   (
       double  K,double *F,double *C,
       double *X,double *Y,double *Z,
       double *U,double *V,double *W,
       polygon_t *Polygon
   )
{
    double f[8],c[8],x[8],y[8],z[8],u[8],v[8],w[8];
    int i,j,k, n=0, above = 0;

    for( i=0; i<20; i++ ) above += F[i]>K;
    if ( above == 0 || above == 20 ) return 0;

    n += elm_8node_brick_isosurface( K,F,C,X,Y,Z,U,V,W,Polygon );
    return n;
}


/******************************************************************************
 *
 *     Name:        elm_20node_brick_initialize( )
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
int elm_20node_brick_initialize()
{
     static char *Name = "ELM_20NODE_BRICK";

     element_type_t ElementDef;

     elm_20node_brick_shape_functions();

     ElementDef.ElementName = Name;
     ElementDef.ElementCode = 820;

     ElementDef.NumberOfNodes = 20;

     ElementDef.NodeU = NodeU;
     ElementDef.NodeV = NodeV;
     ElementDef.NodeW = NodeW;

     ElementDef.PartialU = (double (*)())elm_20node_brick_dndu_fvalue;
     ElementDef.PartialV = (double (*)())elm_20node_brick_dndv_fvalue;
     ElementDef.PartialW = (double (*)())elm_20node_brick_dndw_fvalue;

     ElementDef.FunctionValue = (double (*)())elm_20node_brick_fvalue;
     ElementDef.Triangulate   = (int (*)())elm_20node_brick_triangulate;
     ElementDef.PointInside   = (int (*)())elm_20node_brick_point_inside;
     ElementDef.IsoLine       = (int (*)())elm_20node_brick_isoline;
     ElementDef.IsoSurface    = (int (*)())elm_20node_brick_isosurface;

     return elm_add_element_type( &ElementDef ) ;
}
