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
 * Definition of 10 node tetra element
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
#include "elements.h"


/*
 * 10 node tetra volume element.
 *
 */

static double A[10][10],N[10][10];

static double NodeU[] = { 0.0, 1.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.0 };
static double NodeV[] = { 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5 };
static double NodeW[] = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5 };

/*******************************************************************************
 *
 *     Name:        elm_10node_tetra_triangulate
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
static int elm_10node_tetra_triangulate( geometry_t *geom,element_t *tetra )
{
    element_t triangle;
    int i,j;

    if ( GlobalOptions.VolumeSides )
    {
      int topo[7];

      triangle.DisplayFlag = TRUE;
      triangle.Topology = topo;
      for( i=0; i<MAX_GROUP_IDS; i++ ) triangle.GroupIds[i] = tetra->GroupIds[i];

      for( i=0; i<4; i++ )
      {
          for( j=0; j<6; j++ )
          {
              triangle.Topology[j] = tetra->Topology[ElmTetraFace[i][j]];
          }
          if ( !elm_6node_triangle_triangulate( geom, &triangle, tetra ) ) return FALSE;
      }
    } else {
      if ( !geo_add_edge( geom, tetra->Topology[0],tetra->Topology[4],tetra ) ) return FALSE;
      if ( !geo_add_edge( geom, tetra->Topology[4],tetra->Topology[1],tetra ) ) return FALSE;
      if ( !geo_add_edge( geom, tetra->Topology[1],tetra->Topology[5],tetra ) ) return FALSE;
      if ( !geo_add_edge( geom, tetra->Topology[5],tetra->Topology[2],tetra ) ) return FALSE;
      if ( !geo_add_edge( geom, tetra->Topology[2],tetra->Topology[6],tetra ) ) return FALSE;
      if ( !geo_add_edge( geom, tetra->Topology[6],tetra->Topology[0],tetra ) ) return FALSE;
      if ( !geo_add_edge( geom, tetra->Topology[0],tetra->Topology[7],tetra ) ) return FALSE;
      if ( !geo_add_edge( geom, tetra->Topology[7],tetra->Topology[3],tetra ) ) return FALSE;
      if ( !geo_add_edge( geom, tetra->Topology[1],tetra->Topology[8],tetra ) ) return FALSE;
      if ( !geo_add_edge( geom, tetra->Topology[8],tetra->Topology[3],tetra ) ) return FALSE;
      if ( !geo_add_edge( geom, tetra->Topology[2],tetra->Topology[9],tetra ) ) return FALSE;
      if ( !geo_add_edge( geom, tetra->Topology[9],tetra->Topology[3],tetra ) ) return FALSE;
    }

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        elm_10node_tetra_point_inside
 *                     (
 *                         double *nx,double *ny,double *nz,
 *                         double  px,double  py,double  pz,
 *                         double *u, double *v, double *w
 *                      )
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
 * NOTES: the goal here can be hard for more involved element types.
 *         TODO: SUBDIVISION...
 *
 ******************************************************************************/
int elm_10node_tetra_point_inside
   ( 
                   double *nx, double *ny, double *nz,
         double px, double py, double pz, double *u,double *v,double *w
   )
{
    return elm_4node_tetra_point_inside( nx,ny,nz,px,py,pz,u,v,w );
}

/*******************************************************************************
 *
 *     Name:        elm_10node_tetra_isosurface
 *
 *     Purpose:     Extract isosurfaces for element
 *
 *     Parameters:
 *
 *         Input:   (double )  K: contour threshold
 *                  (double *) F: contour quantity values at nodes 
 *                  (double *) C: color quantity values at nodes 
 *                  (double *) X,Y,Z: node coordinates
 *                  (double *) U,V,W: normal vector at nodes
 *
 *         Output:  (polygon_t *)Polygon, output triangles (0,1 or 2) triangles
 *
 *   Return value:  How many triangles we've got..., TODO: SUBDIVISION
 *
 ******************************************************************************/
int elm_10node_tetra_isosurface
   ( 
       double  K,double *F,double *C, double *X,double *Y,double *Z, 
       double *U,double *V,double *W, polygon_t *Polygon
   )
{
     return elm_4node_tetra_isosurface( K,F,C,X,Y,Z,U,V,W,Polygon );
}

/*******************************************************************************
 *
 *     Name:        elm_10node_tetra_shape_functions
 *
 *     Purpose:     Initialize element shape function array. Internal only.
 *
 *     Parameters:
 *
 *         Input:   Global (filewise) variables NodeU,NodeV,NodeW
 *
 *         Output:  Global (filewise) variable N[4][4], will contain
 *                  shape function coefficients
 *
 *   Return value:  void
 *
 ******************************************************************************/
static void elm_10node_tetra_shape_functions()
{
     double u,v,w;
     int i,j;

     for( i=0; i<10; i++ )
     {
         u = NodeU[i];
         v = NodeV[i];
         w = NodeW[i];

         A[i][0]   = 1;
         A[i][1]   = u;
         A[i][2]   = v;
         A[i][3]   = w;
         A[i][4]   = u*u;
         A[i][5]   = v*v;
         A[i][6]   = w*w;
         A[i][7]   = u*v;
         A[i][8]   = u*w;
         A[i][9]   = v*w;
     }

     lu_mtrinv( (double *)A,10 );

     for( i=0; i<10; i++ )
        for( j=0; j<10; j++ ) N[i][j] = A[j][i];
}

/*******************************************************************************
 *
 *     Name:        elm_10node_tetra_fvalue
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
static double elm_10node_tetra_fvalue(double *F,double u,double v,double w)
{
     double R=0.0,uu=u*u,vv=v*v,ww=w*w,uv=u*v,uw=u*w,vw=v*w;
     int i;

     for( i=0; i<10; i++ )
     {
         R += F[i]*( N[i][0]    +
                     N[i][1]*u  +
                     N[i][2]*v  +
                     N[i][3]*w  +
                     N[i][4]*uu +
                     N[i][5]*vv +
                     N[i][6]*ww +
                     N[i][7]*uv +
                     N[i][8]*uw +
                     N[i][9]*vw );
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_10node_tetra_dndu_fvalue
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
static double elm_10node_tetra_dndu_fvalue(double *F,double u,double v,double w)
{
     double R=0.0,u2=2*u;
     int i;

     for( i=0; i<10; i++ )
     {
         R += F[i]*( N[i][1]+
                     N[i][4]*u2 +
                     N[i][7]*v  +
                     N[i][8]*w  );
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_10node_tetra_dndv_fvalue
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
static double elm_10node_tetra_dndv_fvalue(double *F,double u,double v,double w)
{
     double R=0.0,v2=2*v;
     int i;

     for( i=0; i<10; i++ )
     {
         R += F[i]*( N[i][2]    +
                     N[i][5]*v2 +
                     N[i][7]*u  +
                     N[i][9]*w  );
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_10node_tetra_dndw_fvalue
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
static double elm_10node_tetra_dndw_fvalue(double *F,double u,double v,double w)
{
     double R=0.0,w2=2*w;
     int i;

     for( i=0; i<10; i++ )
     {
         R += F[i]*( N[i][3]    +
                     N[i][6]*w2 +
                     N[i][8]*u  +
                     N[i][9]*v  );
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_10node_tetra_isoline
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
static int elm_10node_tetra_isoline
  (
     double K, double *F, double *C,double *X,double *Y,double *Z,line_t *Line
  )
{
    double f[6],c[6],x[6],y[6],z[6];

    int i, j, k, n=0, above=0;

    for( i=0; i<10; i++ ) above += F[i]>K;
    if ( above == 0 || above == 10 ) return 0;

    for( i=0; i<4; i++ )
    {
        for( j=0; j<6; j++ )
        {
            k = ElmTetraFace[i][j];
            f[j] = F[k];
            c[j] = C[k];
            x[j] = X[k];
            y[j] = Y[k];
            z[j] = Z[k];
        }
        n += elm_6node_triangle_isoline( K,f,c,x,y,z,&Line[n] );
    }

    return n;
}


/******************************************************************************
 *
 *     Name:        elm_10node_tetra_initialize
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
int elm_10node_tetra_initialize()
{
     static char *Name = "ELM_10NODE_TETRA";

     element_type_t ElementDef;

     elm_10node_tetra_shape_functions();

     ElementDef.ElementName = Name;
     ElementDef.ElementCode = 510;

     ElementDef.NumberOfNodes = 10;

     ElementDef.NodeU = NodeU;
     ElementDef.NodeV = NodeV;
     ElementDef.NodeW = NodeW;

     ElementDef.PartialU = (double (*)())elm_10node_tetra_dndu_fvalue;
     ElementDef.PartialV = (double (*)())elm_10node_tetra_dndv_fvalue;
     ElementDef.PartialW = (double (*)())elm_10node_tetra_dndw_fvalue;

     ElementDef.FunctionValue = (double (*)())elm_10node_tetra_fvalue;
     ElementDef.Triangulate   = (int (*)())elm_10node_tetra_triangulate;
     ElementDef.PointInside   = (int (*)())elm_10node_tetra_point_inside;
     ElementDef.IsoLine       = (int (*)())elm_10node_tetra_isoline;
     ElementDef.IsoSurface    = (int (*)())elm_10node_tetra_isosurface;

     return elm_add_element_type( &ElementDef );
}
