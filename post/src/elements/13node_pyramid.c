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
 * Definition of 6 node pyramid element.
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
 * Five node pyramid volume element.
 *
 */

static double A[13][13],N[13][13];

static double NodeU[] = { -1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 1.0, 0.0, -1.0, -0.5,  0.5, 0.5, -0.5 };
static double NodeV[] = { -1.0,-1.0, 1.0,  1.0, 0.0,-1.0, 0.0, 1.0,  0.0, -0.5, -0.5, 0.5,  0.5 };
static double NodeW[] = { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5 };

/*******************************************************************************
 *
 *     Name:        elm_13node_pyramid_triangulate
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
static int elm_13node_pyramid_triangulate( geometry_t *geom,element_t *pyramid )
{
    element_t element;
    int i,j;

    int PrismFace[5][8] = { { 0,1,2,3,5,6,7,8 }, { 0,1,4,5,10, 9, 0, 0 }, 
                                                 { 1,2,4,6,11,10, 0, 0 },
                                                 { 2,3,4,7,12,11, 0, 0 },
                                                 { 3,0,4,8, 9,12, 0, 0 } };
    if ( GlobalOptions.VolumeSides )
    {
       int topo[8];

       element.DisplayFlag = TRUE;
       element.Topology = topo;
       for( i=0; i<MAX_GROUP_IDS; i++ ) element.GroupIds[i] = pyramid->GroupIds[i];
   
        for( j=0; j<8; j++ )
        {
            element.Topology[j] = pyramid->Topology[PrismFace[0][j]];
        }
        if ( !elm_8node_quad_triangulate( geom, &element, pyramid ) ) return FALSE;

        for( i=1; i<5; i++ )
        {
            for( j=0; j<6; j++ )
            {
                element.Topology[j] = pyramid->Topology[PrismFace[i][j]];
            }
            if ( !elm_6node_triangle_triangulate( geom, &element, pyramid ) ) return FALSE;
        }
    } else {
       if ( !geo_add_edge( geom, pyramid->Topology[0], pyramid->Topology[1],pyramid ) ) return FALSE;
       if ( !geo_add_edge( geom, pyramid->Topology[1], pyramid->Topology[2],pyramid ) ) return FALSE;
       if ( !geo_add_edge( geom, pyramid->Topology[2], pyramid->Topology[3],pyramid ) ) return FALSE;
       if ( !geo_add_edge( geom, pyramid->Topology[3], pyramid->Topology[0],pyramid ) ) return FALSE;
       if ( !geo_add_edge( geom, pyramid->Topology[0], pyramid->Topology[4],pyramid ) ) return FALSE;
       if ( !geo_add_edge( geom, pyramid->Topology[1], pyramid->Topology[4],pyramid ) ) return FALSE;
       if ( !geo_add_edge( geom, pyramid->Topology[2], pyramid->Topology[4],pyramid ) ) return FALSE;
       if ( !geo_add_edge( geom, pyramid->Topology[3], pyramid->Topology[4],pyramid ) ) return FALSE;
    }

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        elm_13node_pyramid_shape_functions
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
static void elm_13node_pyramid_shape_functions()
{
     double u,v,w;
     int i,j;

#if 0
     for( i=0; i<13; i++ )
     {
         u = NodeU[i];
         v = NodeV[i];
         w = NodeW[i];

         A[i][0]  = 1;
         A[i][1]  = u;
         A[i][2]  = u*u;
         A[i][3]  = v;
         A[i][4]  = u*v;
         A[i][5]  = u*u*v;
         A[i][6]  = v*v;
         A[i][7]  = u*v*v;
         A[i][8]  = w;
         A[i][9]  = u*w;
         A[i][10] = v*w;
         A[i][11] = u*v*w;
         A[i][12] = w*w;
     }

     lu_mtrinv( (double *)A,13 );

     for( i=0; i<13; i++ )
        for( j=0; j<13; j++ ) N[i][j] = A[j][i];
#endif
}

/*******************************************************************************
 *
 *     Name:        elm_13node_pyramid_fvalue
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
static double elm_13node_pyramid_fvalue(double *F,double u,double v,double w)
{
     double R=0.0,s;
     int i;

     if ( w == 1 ) w = 1.0-1.0e-12;
     s = 1.0 / (1-w);

     R = R + F[0] * (-u-v-1) * ( (1-u) * (1-v) - w + u*v*w * s ) / 4;
     R = R + F[1] * ( u-v-1) * ( (1+u) * (1-v) - w - u*v*w * s ) / 4;
     R = R + F[2] * ( u+v-1) * ( (1+u) * (1+v) - w + u*v*w * s ) / 4;
     R = R + F[3] * (-u+v-1) * ( (1-u) * (1+v) - w - u*v*w * s ) / 4;
     R = R + F[4] * w*(2*w-1);
     R = R + F[5] * (1+u-w)*(1-u-w)*(1-v-w) * s / 2;
     R = R + F[6] * (1+v-w)*(1-v-w)*(1+u-w) * s / 2;
     R = R + F[7] * (1+u-w)*(1-u-w)*(1+v-w) * s / 2;
     R = R + F[8] * (1+v-w)*(1-v-w)*(1-u-w) * s / 2;
     R = R + F[9]  * w * (1-u-w) * (1-v-w) * s;
     R = R + F[10] * w * (1+u-w) * (1-v-w) * s;
     R = R + F[11] * w * (1+u-w) * (1+v-w) * s;
     R = R + F[12] * w * (1-u-w) * (1+v-w) * s;


     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_13node_pyramid_dndu_fvalue
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
static double elm_13node_pyramid_dndu_fvalue(double *F,double u,double v,double w)
{
     double R=0.0,s;
     int i;

     if ( w == 1 ) w = 1.0-1.0e-12;
     s = 1.0 / (1-w);

     R = R + F[0] * ( -( (1-u) * (1-v) - w + u*v*w * s ) +
              (-u-v-1) * ( -(1-v) + v*w * s ) ) / 4;

     R = R + F[1] * (  ( (1+u) * (1-v) - w - u*v*w * s ) +
              ( u-v-1) * (  (1-v) - v*w * s ) ) / 4;

     R = R + F[2] * (  ( (1+u) * (1+v) - w + u*v*w * s ) +
              ( u+v-1) * (  (1+v) + v*w * s ) ) / 4;

     R = R + F[3] * ( -( (1-u) * (1+v) - w - u*v*w * s ) +
              (-u+v-1) * ( -(1+v) - v*w * s ) ) / 4;

     R = R + F[4] * 0.0;

     R = R + F[5] * (  (1-u-w)*(1-v-w) - (1+u-w)*(1-v-w) ) * s / 2;
     R = R + F[6] * (  (1+v-w)*(1-v-w) ) * s / 2;
     R = R + F[7] * (  (1-u-w)*(1+v-w) - (1+u-w)*(1+v-w) ) * s / 2;
     R = R + F[8] * ( -(1+v-w)*(1-v-w) ) * s / 2;

     R = R - F[ 9] * w * (1-v-w) * s;
     R = R + F[10] * w * (1-v-w) * s;
     R = R + F[11] * w * (1+v-w) * s;
     R = R - F[12] * w * (1+v-w) * s;

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_13node_pyramid_dndv_fvalue
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
static double elm_13node_pyramid_dndv_fvalue(double *F,double u,double v,double w)
{
     double R=0.0,s;
     int i;

     if ( w == 1 ) w = 1.0-1.0e-12;
     s = 1.0 / (1-w);

     R = R + F[0] * ( -( (1-u) * (1-v) - w + u*v*w * s ) + 
   	   (-u-v-1) * ( -(1-u) + u*w * s ) ) / 4;

     R = R + F[1] * ( -( (1+u) * (1-v) - w - u*v*w * s ) +
	   ( u-v-1) * ( -(1+u) - u*w * s ) ) / 4;

     R = R + F[2] * (  ( (1+u) * (1+v) - w + u*v*w * s ) +
	   ( u+v-1) * (	 (1+u) + u*w * s ) ) / 4;

     R = R + F[3] * (  ( (1-u) * (1+v) - w - u*v*w * s ) +
	   (-u+v-1) * (	 (1-u) - u*w * s ) ) / 4;

     R = R + F[4] * 0.0;

     R = R - F[5] *  (1+u-w)*(1-u-w) * s / 2;
     R = R + F[6] * ( (1-v-w)*(1+u-w) - (1+v-w)*(1+u-w) ) * s / 2;
     R = R + F[7] *  (1+u-w)*(1-u-w) * s / 2;
     R = R + F[8] * ( (1-v-w)*(1-u-w) - (1+v-w)*(1-u-w) ) * s / 2;

     R = R - F[9]  * w * (1-u-w) * s;
     R = R - F[10] * w * (1+u-w) * s;
     R = R + F[11] * w * (1+u-w) * s;
     R = R + F[12] * w * (1-u-w) * s;

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_13node_pyramid_dndw_fvalue
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
static double elm_13node_pyramid_dndw_fvalue(double *F,double u,double v,double w)
{
     double R=0.0,s;
     int i;

     if ( w == 1 ) w = 1.0-1.0e-12;
     s = 1.0 / (1-w);

     R = R + F[0] * (-u-v-1) * ( -1 + u*v*s*s ) / 4;
     R = R + F[1] * ( u-v-1) * ( -1 - u*v*s*s ) / 4;
     R = R + F[2] * ( u+v-1) * ( -1 + u*v*s*s ) / 4;
     R = R + F[3] * (-u+v-1) * ( -1 - u*v*s*s ) / 4;

     R = R + F[4] * (4*w-1);

     R = R + F[5] * ( ( -(1-u-w)*(1-v-w) - (1+u-w)*(1-v-w) - (1+u-w)*(1-u-w) ) * s +
                       ( 1+u-w)*(1-u-w)*(1-v-w) * s*s ) / 2;

     R = R + F[6] * ( ( -(1-v-w)*(1+u-w) - (1+v-w)*(1+u-w) - (1+v-w)*(1-v-w) ) * s +
                       ( 1+v-w)*(1-v-w)*(1+u-w) * s*s ) / 2;

     R = R + F[7] * ( ( -(1-u-w)*(1+v-w) - (1+u-w)*(1+v-w) - (1+u-w)*(1-u-w) ) * s +
                       ( 1+u-w)*(1-u-w)*(1+v-w) * s*s ) / 2;

     R = R + F[8] * ( ( -(1-v-w)*(1-u-w) - (1+v-w)*(1-u-w) - (1+v-w)*(1-v-w) ) * s +
                       ( 1+v-w)*(1-v-w)*(1-u-w) * s*s ) / 2;

     R = R + F[ 9] * ( ( (1-u-w) * (1-v-w) - w * (1-v-w) - w * (1-u-w) ) * s  +
                      w * (1-u-w) * (1-v-w) * s*s );

     R = R + F[10] * ( ( (1+u-w) * (1-v-w) - w * (1-v-w) - w * (1+u-w) ) * s  +
                      w * (1+u-w) * (1-v-w) * s*s );

     R = R + F[11] * ( ( (1+u-w) * (1+v-w) - w * (1+v-w) - w * (1+u-w) ) * s  +
                      w * (1+u-w) * (1+v-w) * s*s );

     R = R + F[12] * ( ( (1-u-w) * (1+v-w) - w * (1+v-w) - w * (1-u-w) ) * s  +
                      w * (1-u-w) * (1+v-w) * s*s );


     return R;
}

/******************************************************************************
 *
 *     Name:        elm_13node_pyramid_initialize
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
int elm_13node_pyramid_initialize()
{
     static char *Name = "ELM_13NODE_PYRAMID";

     element_type_t ElementDef;

     elm_13node_pyramid_shape_functions();

     ElementDef.ElementName = Name;
     ElementDef.ElementCode = 613;

     ElementDef.NumberOfNodes = 13;

     ElementDef.NodeU = NodeU;
     ElementDef.NodeV = NodeV;
     ElementDef.NodeW = NodeW;

     ElementDef.PartialU = (double (*)())elm_13node_pyramid_dndu_fvalue;
     ElementDef.PartialV = (double (*)())elm_13node_pyramid_dndv_fvalue;
     ElementDef.PartialW = (double (*)())elm_13node_pyramid_dndw_fvalue;

     ElementDef.FunctionValue = (double (*)())elm_13node_pyramid_fvalue;
     ElementDef.Triangulate   = (int (*)())elm_13node_pyramid_triangulate;
     ElementDef.PointInside   = (int (*)())NULL;

     return elm_add_element_type( &ElementDef );
}
