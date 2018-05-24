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
 * Definition of 4 node bar element.
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
 * Two node 1D element
 * 
 *  o---o---o---o u
 *  0    0.5    1
 *
 */


static double NodeU[] = { 0.0, 1.0, 1.0/3.0, 2.0/3.0 };

static double N[4][4],A[4][4];

/*******************************************************************************
 *
 *     Name:        elm_4node_bar_shape_functions( )
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
static void elm_4node_bar_shape_functions()
{
     double u,v;

     int i,j;

     for( i=0; i<4; i++ )
     {
         u = NodeU[i];

         A[i][0]  = 1;
         A[i][1]  = u;
         A[i][2]  = u*u;
         A[i][3]  = u*u*u;
     }

     lu_mtrinv( (double *)A,4 );

     for( i=0; i<4; i++ )
        for( j=0; j<4; j++ ) N[i][j] = A[j][i];
}


/*******************************************************************************
 *
 *     Name:        elm_4node_bar_triangulate( geometry_t *,element_t * )
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
int elm_4node_bar_triangulate( geometry_t *geom, element_t *Elm, element_t *Parent)
{
    geo_add_edge( geom, Elm->Topology[0],Elm->Topology[2],Parent );
    geo_add_edge( geom, Elm->Topology[2],Elm->Topology[3],Parent );
    return geo_add_edge( geom, Elm->Topology[3],Elm->Topology[1],Parent );
}


/*******************************************************************************
 *
 *     Name:        elm_4node_bar_fvalue( double *,double,double )
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
static double elm_4node_bar_fvalue( double *F,double u)
{
     double R=0.0;
     int i;

     for( i=0; i<4; i++ )
     {
         R += F[i]*( N[i][0]    +
                     N[i][1]*u  +
                     N[i][2]*u*u  +
                     N[i][3]*u*u*u );
     } 

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_4node_bar_dndu_fvalue( double *,double,double )
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
static double elm_4node_bar_dndu_fvalue(double *F,double u)
{
     double R=0.0;
     int i;

     for( i=0; i<4; i++ )
     {
         R += F[i]*( N[i][1] + 2*N[i][2]*u + 3*N[i][3]*u*u );
     }

     return R;
}

/*******************************************************************************
 *
 *     Name:        elm_4node_bar_initialize()
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
int elm_4node_bar_initialize()
{
     static char *Name = "ELM_4NODE_LINE";

     element_type_t ElementDef;

     elm_4node_bar_shape_functions();

     ElementDef.ElementName = Name;
     ElementDef.ElementCode = 204;

     ElementDef.NumberOfNodes = 4;

     ElementDef.NodeU = NodeU;
     ElementDef.NodeV = NULL;
     ElementDef.NodeW = NULL;

     ElementDef.PartialU = (double (*)())elm_4node_bar_dndu_fvalue;
     ElementDef.PartialV = (double (*)())NULL;
     ElementDef.PartialW = (double (*)())NULL;

     ElementDef.FunctionValue = (double (*)())elm_4node_bar_fvalue;
     ElementDef.Triangulate   = (int (*)())elm_4node_bar_triangulate;
     ElementDef.IsoLine       = (int (*)())NULL;
     ElementDef.PointInside   = (int (*)())NULL;
     ElementDef.IsoSurface    = (int (*)())NULL;

     return elm_add_element_type( &ElementDef ) ;
}
