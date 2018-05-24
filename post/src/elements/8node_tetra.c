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
 * Definition of 4 node tetra element
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
 * Four node tetra volume element.
 *
 *             3
 *            /|\
 *           / | \
 *          /  |  \
 *         /   |   \
 *        /    2    \
 *       /    / \    \
 *      /   /     \   \            w  v
 *     /  /         \  \           | /
 *    //              \ \          |/
 *    0-----------------1          ----u
 *
 *
 */

static double A[4][4],N[4][4];

static double NodeU[] = { 0.0, 1.0, 0.0, 0.0 };
static double NodeV[] = { 0.0, 0.0, 1.0, 0.0 };
static double NodeW[] = { 0.0, 0.0, 0.0, 1.0 };


/*******************************************************************************
 *
 *     Name:        elm_8node_tetra_triangulate
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
static int elm_8node_tetra_triangulate( geometry_t *geom,element_t *tetra )
{
    element_t triangle;
    int i,j;


    if ( GlobalOptions.VolumeSides )
    {
      int topo[4];

      triangle.DisplayFlag = TRUE;
      triangle.Topology = topo;
      for( i=0; i<MAX_GROUP_IDS; i++ ) triangle.GroupIds[i] = tetra->GroupIds[i];

      for( i=0; i<4; i++ )
      {
          for( j=0; j<3; j++ )
          {
              triangle.Topology[j] = tetra->Topology[ElmTetraFace[i][j]];
          }
          if ( !elm_3node_triangle_triangulate( geom, &triangle, tetra ) ) return FALSE;
      }
    } else {
      if ( !geo_add_edge( geom, tetra->Topology[0],tetra->Topology[1],tetra ) ) return FALSE;
      if ( !geo_add_edge( geom, tetra->Topology[0],tetra->Topology[2],tetra ) ) return FALSE;
      if ( !geo_add_edge( geom, tetra->Topology[1],tetra->Topology[2],tetra ) ) return FALSE;
      if ( !geo_add_edge( geom, tetra->Topology[0],tetra->Topology[3],tetra ) ) return FALSE;
      if ( !geo_add_edge( geom, tetra->Topology[1],tetra->Topology[3],tetra ) ) return FALSE;
      if ( !geo_add_edge( geom, tetra->Topology[2],tetra->Topology[3],tetra ) ) return FALSE;
    }

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        elm_8node_tetra_point_inside
 *                     (
 *                         double *nx,double *ny,double *nz,
 *                         double px, double py, double pz,
 *                          double *u,double *v,double *w
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
 * NOTES: the goal here can be hard for more involved element types. kind of
 *        trivial for this one... 
 *
 ******************************************************************************/
int elm_8node_tetra_point_inside
   ( 
                   double *nx, double *ny, double *nz,
         double px, double py, double pz, double *u,double *v,double *w
   )
{
    double B00,B01,B02,B10,B11,B12,B20,B21,B22;
    double A00,A01,A02,A10,A11,A12,A20,A21,A22,detA;

    if ( px > MAX(MAX(MAX(nx[0],nx[1]),nx[2]),nx[3])) return FALSE;
    if ( py > MAX(MAX(MAX(ny[0],ny[1]),ny[2]),ny[3])) return FALSE;
    if ( pz > MAX(MAX(MAX(nz[0],nz[1]),nz[2]),nz[3])) return FALSE;

    if ( px < MIN(MIN(MIN(nx[0],nx[1]),nx[2]),nx[3])) return FALSE;
    if ( py < MIN(MIN(MIN(ny[0],ny[1]),ny[2]),ny[3])) return FALSE;
    if ( pz < MIN(MIN(MIN(nz[0],nz[1]),nz[2]),nz[3])) return FALSE;

    A00 = nx[1] - nx[0];
    A01 = nx[2] - nx[0];
    A02 = nx[3] - nx[0];

    A10 = ny[1] - ny[0];
    A11 = ny[2] - ny[0];
    A12 = ny[3] - ny[0];

    A20 = nz[1] - nz[0];
    A21 = nz[2] - nz[0];
    A22 = nz[3] - nz[0];

    detA  = A00*(A11*A22 - A12*A21);
    detA += A01*(A12*A20 - A10*A22);
    detA += A02*(A10*A21 - A11*A20);
    if ( ABS(detA) < AEPS )
    {
        fprintf( stderr, "8node_tetra_inside: SINGULAR (huh?).\n" );
        return FALSE;
    }

    detA = 1 / detA;

    B00 = (A11*A22 - A12*A21)*detA;
    B10 = (A12*A20 - A10*A22)*detA;
    B20 = (A10*A21 - A11*A20)*detA;

    B01 = (A21*A02 - A01*A22)*detA;
    B11 = (A00*A22 - A20*A02)*detA;
    B21 = (A01*A20 - A00*A21)*detA;

    B02 = (A01*A12 - A11*A02)*detA;
    B12 = (A10*A02 - A00*A12)*detA;
    B22 = (A00*A11 - A10*A01)*detA;

    px -= nx[0];
    py -= ny[0];
    pz -= nz[0];

    *u = B00*px + B01*py + B02*pz;
    if ( *u < 0.0 || *u > 1.0 ) return FALSE;

    *v = B10*px + B11*py + B12*pz;
    if ( *v < 0.0 || *v > 1.0 ) return FALSE;

    *w = B20*px + B21*py + B22*pz;
    if ( *w < 0.0 || *w > 1.0 ) return FALSE;

    return (*u+*v+*w) <= 1.0;
}


/*******************************************************************************
 *
 *     Name:        elm_8node_tetra_isoline
 *
 *     Purpose:     Extract a isoline from triangle with given threshold 
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
static int elm_8node_tetra_isoline
  (
     double K, double *F, double *C,double *X,double *Y,double *Z,line_t *Line
  )
{
    double f[3],c[3],x[3],y[3],z[3];

    int i, j, k, n=0, above=0;

    for( i=0; i<4; i++ ) above += F[i]>K;
    if ( above == 0 || above == 4 ) return 0;

    for( i=0; i<4; i++ )
    {
        for( j=0; j<3; j++ )
        {
            k = ElmTetraFace[i][j];
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
 *     Name:        elm_8node_tetra_isosurface
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
 *   Return value:  How many triangles we've got...
 *
 ******************************************************************************/
int elm_8node_tetra_isosurface
   ( 
       double  K,double *F,double *C, double *X,double *Y,double *Z, 
            double *U,double *V,double *W, polygon_t *Polygon
   )
{
    double t,tx[4],ty[4],tz[4],tu[4],tv[4],tw[4],tf[4],tc[4];
    double ax,ay,az,bx,by,bz,nx,ny,nz;

	int S0 = F[0] > K;
    int S1 = F[1] > K;
    int S2 = F[2] > K;
    int S3 = F[3] > K;

    int S = S0+S1+S2+S3,I[4],j;

    if ( S==0 || S==4 ) return 0;

    if ( S==1 || S==3 )
    {
        if ( (S==1 && S0) || (S==3 && !S0) )
        {
            I[0] = 0;
            I[1] = 1;
            I[2] = 2;
            I[3] = 3;
        } else if ( (S==1 && S1) || (S==3 && !S1) )
        {
            I[0] = 1;
            I[1] = 0;
            I[2] = 2;
            I[3] = 3;
        } else if ( (S==1 && S2) || (S==3 && !S2) )
        {
            I[0] = 2;
            I[1] = 0;
            I[2] = 1;
            I[3] = 3;
        } else if ( (S==1 && S3) || (S==3 && !S3) )
        {
            I[0] = 3;
            I[1] = 0;
            I[2] = 1;
            I[3] = 2;
        } else { return 0; }
   
        for( j=1; j<4; j++ )
        {
            t = (K-F[I[0]]) / (F[I[j]]-F[I[0]]);
            Polygon->x[j-1] = t*(X[I[j]]-X[I[0]]) + X[I[0]];
            Polygon->y[j-1] = t*(Y[I[j]]-Y[I[0]]) + Y[I[0]];
            Polygon->z[j-1] = t*(Z[I[j]]-Z[I[0]]) + Z[I[0]];
            Polygon->u[j-1] = t*(U[I[j]]-U[I[0]]) + U[I[0]];
            Polygon->v[j-1] = t*(V[I[j]]-V[I[0]]) + V[I[0]];
            Polygon->w[j-1] = t*(W[I[j]]-W[I[0]]) + W[I[0]];
            Polygon->c[j-1] = t*(C[I[j]]-C[I[0]]) + C[I[0]];
            Polygon->f[j-1] = K;
        }

        ax = Polygon->x[1] - Polygon->x[0];
        ay = Polygon->y[1] - Polygon->y[0];
        az = Polygon->z[1] - Polygon->z[0];

        bx = Polygon->x[2] - Polygon->x[0];
        by = Polygon->y[2] - Polygon->y[0];
        bz = Polygon->z[2] - Polygon->z[0];

        nx = ay*bz - az*by;
        ny = az*bx - ax*bz;
        nz = ax*by - ay*bx;

        ax = Polygon->u[0] + Polygon->u[1] + Polygon->u[2];
        ay = Polygon->v[0] + Polygon->v[1] + Polygon->v[2];
        az = Polygon->w[0] + Polygon->w[1] + Polygon->w[2];

        if ( nx*ax + ny*ay + nz*az < 0.0 )
        {
            double s;

#define swap( x,y ) { s=x; x=y; y=s; }

            swap( Polygon->x[1], Polygon->x[2] );
            swap( Polygon->y[1], Polygon->y[2] );
            swap( Polygon->z[1], Polygon->z[2] );
            swap( Polygon->u[1], Polygon->u[2] );
            swap( Polygon->v[1], Polygon->v[2] );
            swap( Polygon->w[1], Polygon->w[2] );
            swap( Polygon->f[1], Polygon->f[2] );
            swap( Polygon->c[1], Polygon->c[2] );

#undef swap
        }

        return 1;
    } else
    {
        if ( (S0 && S1) || (!S0 && !S1) )
        {
            t = (K-F[0])/ (F[2]-F[0]);
            tx[0] = t*(X[2]-X[0]) + X[0];
            ty[0] = t*(Y[2]-Y[0]) + Y[0];
            tz[0] = t*(Z[2]-Z[0]) + Z[0];
            tu[0] = t*(U[2]-U[0]) + U[0];
            tv[0] = t*(V[2]-V[0]) + V[0];
            tw[0] = t*(W[2]-W[0]) + W[0];
            tc[0] = t*(C[2]-C[0]) + C[0];
            tf[0] = K;

            t = (K-F[1]) / (F[2]-F[1]);
            tx[1] = t*(X[2]-X[1]) + X[1];
            ty[1] = t*(Y[2]-Y[1]) + Y[1];
            tz[1] = t*(Z[2]-Z[1]) + Z[1];
            tu[1] = t*(U[2]-U[1]) + U[1];
            tv[1] = t*(V[2]-V[1]) + V[1];
            tw[1] = t*(W[2]-W[1]) + W[1];
            tc[1] = t*(C[2]-C[1]) + C[1];
            tf[1] = K;

            t = (K-F[1]) / (F[3]-F[1]);
            tx[2] = t*(X[3]-X[1]) + X[1];
            ty[2] = t*(Y[3]-Y[1]) + Y[1];
            tz[2] = t*(Z[3]-Z[1]) + Z[1];
            tu[2] = t*(U[3]-U[1]) + U[1];
            tv[2] = t*(V[3]-V[1]) + V[1];
            tw[2] = t*(W[3]-W[1]) + W[1];
            tc[2] = t*(C[3]-C[1]) + C[1];
            tf[2] = K;

            t = (K-F[0]) / (F[3]-F[0]);
            tx[3] = t*(X[3]-X[0]) + X[0];
            ty[3] = t*(Y[3]-Y[0]) + Y[0];
            tz[3] = t*(Z[3]-Z[0]) + Z[0];
            tu[3] = t*(U[3]-U[0]) + U[0];
            tv[3] = t*(V[3]-V[0]) + V[0];
            tw[3] = t*(W[3]-W[0]) + W[0];
            tc[3] = t*(C[3]-C[0]) + C[0];
            tf[3] = K;
        }
        else if ( (S0 && S2) || (!S0 && !S2) )
        {
            t = (K-F[0]) / (F[1]-F[0]);
            tx[0] = t*(X[1]-X[0]) + X[0];
            ty[0] = t*(Y[1]-Y[0]) + Y[0];
            tz[0] = t*(Z[1]-Z[0]) + Z[0];
            tu[0] = t*(U[1]-U[0]) + U[0];
            tv[0] = t*(V[1]-V[0]) + V[0];
            tw[0] = t*(W[1]-W[0]) + W[0];
            tc[0] = t*(C[1]-C[0]) + C[0];
            tf[0] = K;

            t = (K-F[2]) / (F[1]-F[2]);
            tx[1] = t*(X[1]-X[2]) + X[2];
            ty[1] = t*(Y[1]-Y[2]) + Y[2];
            tz[1] = t*(Z[1]-Z[2]) + Z[2];
            tu[1] = t*(U[1]-U[2]) + U[2];
            tv[1] = t*(V[1]-V[2]) + V[2];
            tw[1] = t*(W[1]-W[2]) + W[2];
            tc[1] = t*(C[1]-C[2]) + C[2];
            tf[1] = K;

            t = (K-F[2]) / (F[3]-F[2]);
            tx[2] = t*(X[3]-X[2]) + X[2];
            ty[2] = t*(Y[3]-Y[2]) + Y[2];
            tz[2] = t*(Z[3]-Z[2]) + Z[2];
            tu[2] = t*(U[3]-U[2]) + U[2];
            tv[2] = t*(V[3]-V[2]) + V[2];
            tw[2] = t*(W[3]-W[2]) + W[2];
            tc[2] = t*(C[3]-C[2]) + C[2];
            tf[2] = K;

            t = (K-F[0]) / (F[3]-F[0]);
            tx[3] = t*(X[3]-X[0]) + X[0];
            ty[3] = t*(Y[3]-Y[0]) + Y[0];
            tz[3] = t*(Z[3]-Z[0]) + Z[0];
            tu[3] = t*(U[3]-U[0]) + U[0];
            tv[3] = t*(V[3]-V[0]) + V[0];
            tw[3] = t*(W[3]-W[0]) + W[0];
            tc[3] = t*(C[3]-C[0]) + C[0];
            tf[3] = K;
        }
        else if ( (S0 && S3) || (!S0 && !S3) )
        {
            t = (K-F[0]) / (F[1]-F[0]);
            tx[0] = t*(X[1]-X[0]) + X[0];
            ty[0] = t*(Y[1]-Y[0]) + Y[0];
            tz[0] = t*(Z[1]-Z[0]) + Z[0];
            tu[0] = t*(U[1]-U[0]) + U[0];
            tv[0] = t*(V[1]-V[0]) + V[0];
            tw[0] = t*(W[1]-W[0]) + W[0];
            tc[0] = t*(C[1]-C[0]) + C[0];
            tf[0] = K;

            t = (K-F[3]) / (F[1]-F[3]);
            tx[1] = t*(X[1]-X[3]) + X[3];
            ty[1] = t*(Y[1]-Y[3]) + Y[3];
            tz[1] = t*(Z[1]-Z[3]) + Z[3];
            tu[1] = t*(U[1]-U[3]) + U[3];
            tv[1] = t*(V[1]-V[3]) + V[3];
            tw[1] = t*(W[1]-W[3]) + W[3];
            tc[1] = t*(C[1]-C[3]) + C[3];
            tf[1] = K;

            t = (K-F[3]) / (F[2]-F[3]);
            tx[2] = t*(X[2]-X[3]) + X[3];
            ty[2] = t*(Y[2]-Y[3]) + Y[3];
            tz[2] = t*(Z[2]-Z[3]) + Z[3];
            tu[2] = t*(U[2]-U[3]) + U[3];
            tv[2] = t*(V[2]-V[3]) + V[3];
            tw[2] = t*(W[2]-W[3]) + W[3];
            tc[2] = t*(C[2]-C[3]) + C[3];
            tf[2] = K;

            t = (K-F[0]) / (F[2]-F[0]);
            tx[3] = t*(X[2]-X[0]) + X[0];
            ty[3] = t*(Y[2]-Y[0]) + Y[0];
            tz[3] = t*(Z[2]-Z[0]) + Z[0];
            tu[3] = t*(U[2]-U[0]) + U[0];
            tv[3] = t*(V[2]-V[0]) + V[0];
            tw[3] = t*(W[2]-W[0]) + W[0];
            tc[3] = t*(C[2]-C[0]) + C[0];
            tf[3] = K;
        }

        Polygon[0].x[0] = tx[0];
        Polygon[0].y[0] = ty[0];
        Polygon[0].z[0] = tz[0];
        Polygon[0].u[0] = tu[0];
        Polygon[0].v[0] = tv[0];
        Polygon[0].w[0] = tw[0];
        Polygon[0].f[0] = tf[0];
        Polygon[0].c[0] = tc[0];

        ax = tx[1] - tx[0];
        ay = ty[1] - ty[0];
        az = tz[1] - tz[0];

        bx = tx[2] - tx[0];
        by = ty[2] - ty[0];
        bz = tz[2] - tz[0];

        nx = ay*bz - az*by;
        ny = az*bx - ax*bz;
        nz = ax*by - ay*bx;

        ax = tu[0] + tu[1] + tu[2] + tu[3];
        ay = tv[0] + tv[1] + tv[2] + tv[3];
        az = tw[0] + tw[1] + tw[2] + tw[3];

        if ( nx*ax + ny*ay + nz*az >= 0.0 )
        {
            Polygon[0].x[1] = tx[1];
            Polygon[0].y[1] = ty[1];
            Polygon[0].z[1] = tz[1];
            Polygon[0].u[1] = tu[1];
            Polygon[0].v[1] = tv[1];
            Polygon[0].w[1] = tw[1];
            Polygon[0].f[1] = tf[1];
            Polygon[0].c[1] = tc[1];

            Polygon[0].x[2] = tx[2];
            Polygon[0].y[2] = ty[2];
            Polygon[0].z[2] = tz[2];
            Polygon[0].u[2] = tu[2];
            Polygon[0].v[2] = tv[2];
            Polygon[0].w[2] = tw[2];
            Polygon[0].f[2] = tf[2];
            Polygon[0].c[2] = tc[2];
        } else 
        {
            Polygon[0].x[1] = tx[2];
            Polygon[0].y[1] = ty[2];
            Polygon[0].z[1] = tz[2];
            Polygon[0].u[1] = tu[2];
            Polygon[0].v[1] = tv[2];
            Polygon[0].w[1] = tw[2];
            Polygon[0].f[1] = tf[2];
            Polygon[0].c[1] = tc[2];

            Polygon[0].x[2] = tx[1];
            Polygon[0].y[2] = ty[1];
            Polygon[0].z[2] = tz[1];
            Polygon[0].u[2] = tu[1];
            Polygon[0].v[2] = tv[1];
            Polygon[0].w[2] = tw[1];
            Polygon[0].f[2] = tf[1];
            Polygon[0].c[2] = tc[1];
        }

        Polygon[1].x[0] = tx[0];
        Polygon[1].y[0] = ty[0];
        Polygon[1].z[0] = tz[0];
        Polygon[1].u[0] = tu[0];
        Polygon[1].v[0] = tv[0];
        Polygon[1].w[0] = tw[0];
        Polygon[1].f[0] = tf[0];
        Polygon[1].c[0] = tc[0];

        ax = tx[2] - tx[0];
        ay = ty[2] - ty[0];
        az = tz[2] - tz[0];

        bx = tx[3] - tx[0];
        by = ty[3] - ty[0];
        bz = tz[3] - tz[0];

        nx = ay*bz - az*by;
        ny = az*bx - ax*bz;
        nz = ax*by - ay*bx;

        ax = tu[0] + tu[1] + tu[2] + tu[3];
        ay = tv[0] + tv[1] + tv[2] + tv[3];
        az = tw[0] + tw[1] + tw[2] + tw[3];

        if ( nx*ax + ny*ay + nz*az >= 0.0 )
        {
            Polygon[1].x[1] = tx[2];
            Polygon[1].y[1] = ty[2];
            Polygon[1].z[1] = tz[2];
            Polygon[1].u[1] = tu[2];
            Polygon[1].v[1] = tv[2];
            Polygon[1].w[1] = tw[2];
            Polygon[1].f[1] = tf[2];
            Polygon[1].c[1] = tc[2];

            Polygon[1].x[2] = tx[3];
            Polygon[1].y[2] = ty[3];
            Polygon[1].z[2] = tz[3];
            Polygon[1].u[2] = tu[3];
            Polygon[1].v[2] = tv[3];
            Polygon[1].w[2] = tw[3];
            Polygon[1].f[2] = tf[3];
            Polygon[1].c[2] = tc[3];
        } else
        {
            Polygon[1].x[1] = tx[3];
            Polygon[1].y[1] = ty[3];
            Polygon[1].z[1] = tz[3];
            Polygon[1].u[1] = tu[3];
            Polygon[1].v[1] = tv[3];
            Polygon[1].w[1] = tw[3];
            Polygon[1].f[1] = tf[3];
            Polygon[1].c[1] = tc[3];

            Polygon[1].x[2] = tx[2];
            Polygon[1].y[2] = ty[2];
            Polygon[1].z[2] = tz[2];
            Polygon[1].u[2] = tu[2];
            Polygon[1].v[2] = tv[2];
            Polygon[1].w[2] = tw[2];
            Polygon[1].f[2] = tf[2];
            Polygon[1].c[2] = tc[2];
        }

        return 2;
    }

    return 0;
}

/*******************************************************************************
 *
 *     Name:        elm_8node_tetra_shape_functions
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
static void elm_8node_tetra_shape_functions()
{
     double u,v,w;
     int i,j;

     for( i=0; i<4; i++ )
     {
         u = NodeU[i];
         v = NodeV[i];
         w = NodeW[i];

         A[i][0]   = 1;
         A[i][1]   = u;
         A[i][2]   = v;
         A[i][3]   = w;
     }

     lu_mtrinv( (double *)A,4 );

     for( i=0; i<4; i++ )
         for( j=0; j<4; j++ ) N[i][j] = A[j][i];
}

/*******************************************************************************
 *
 *     Name:        elm_8node_tetra_fvalue
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
static double elm_8node_tetra_fvalue(double *F,double u,double v,double w)
{
     double R=0.0;
     int i;

     for( i=0; i<4; i++ )
     {
         R += F[i]*(N[i][0]+N[i][1]*u+N[i][2]*v+N[i][3]*w);
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_8node_tetra_dndu_fvalue
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
static double elm_8node_tetra_dndu_fvalue(double *F,double u,double v,double w)
{
     double R=0.0;
     int i;

     for( i=0; i<4; i++ ) R += F[i]*N[i][1];

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_8node_tetra_dndv_fvalue
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
static double elm_8node_tetra_dndv_fvalue(double *F,double u,double v,double w)
{
     double R=0.0;
     int i;

     for( i=0; i<4; i++ ) R += F[i]*N[i][2];

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_8node_tetra_dndw_fvalue
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
static double elm_8node_tetra_dndw_fvalue(double *F,double u,double v,double w)
{
     double R=0.0;
     int i;

     for( i=0; i<4; i++ ) R += F[i]*N[i][3];

     return R;
}

/******************************************************************************
 *
 *     Name:        elm_8node_tetra_initialize
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
int elm_8node_tetra_initialize()
{
     static char *Name = "ELM_8NODE_TETRA";

     element_type_t ElementDef;

     elm_8node_tetra_shape_functions();

     ElementDef.ElementName = Name;
     ElementDef.ElementCode = 508;

     ElementDef.NumberOfNodes = 8;

     ElementDef.NodeU = NodeU;
     ElementDef.NodeV = NodeV;
     ElementDef.NodeW = NodeW;

     ElementDef.PartialU = (double (*)())elm_8node_tetra_dndu_fvalue;
     ElementDef.PartialV = (double (*)())elm_8node_tetra_dndv_fvalue;
     ElementDef.PartialW = (double (*)())elm_8node_tetra_dndw_fvalue;

     ElementDef.FunctionValue = (double (*)())elm_8node_tetra_fvalue;
     ElementDef.Triangulate   = (int (*)())elm_8node_tetra_triangulate;
     ElementDef.IsoLine       = (int (*)())elm_8node_tetra_isoline;
     ElementDef.IsoSurface    = (int (*)())elm_8node_tetra_isosurface;
     ElementDef.PointInside   = (int (*)())elm_8node_tetra_point_inside;

     return elm_add_element_type( &ElementDef );
}

