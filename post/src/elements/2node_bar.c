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
 * Definition of two node bar 1D element.
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
 *  o---------o u
 *  0         1
 *
 */

static double NodeU[2] = { 0.0, 1.0 };

/*******************************************************************************
 *
 *     Name:        elm_2node_bar_triangulate( geometry_t *,element_t * )
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
int elm_2node_bar_triangulate( geometry_t *geom, element_t *Elm, element_t *Parent)
{

    return geo_add_edge( geom, Elm->Topology[0],Elm->Topology[1],Parent );
}

/*******************************************************************************
 *
 *     Name:        elm_2node_bar_fvalue( double *,double,double )
 *
 *     Purpose:     return value of a quantity given on nodes at point (u)
 *                 
 *
 *     Parameters:
 *
 *         Input:  (double *) quantity values at nodes 
 *                 (double u) point where value is evaluated
 *
 *         Output:  none
 *
 *   Return value:  quantity value
 *
 ******************************************************************************/
static double elm_2node_bar_fvalue(double *F,double u)
{
    return F[0]*(1-u) + F[1]*u;
}

/*******************************************************************************
 *
 *     Name:        elm_2node_bar_dndu_fvalue( double *,double,double )
 *
 *     Purpose:     return value of a first partial derivate in (u) of a
 *                  quantity given on nodes at point (u)
 *                 
 *
 *     Parameters:
 *
 *         Input:  (double *) quantity values at nodes 
 *                 (double u) point where value is evaluated
 *
 *         Output:  none
 *
 *   Return value:  quantity value
 *
 ******************************************************************************/
static double elm_2node_bar_dndu_fvalue(double *F,double u)
{
    return -F[0] + F[1];
}

/*******************************************************************************
 *
 *     Name:        elm_2node_bar_initialize()
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
int elm_2node_bar_initialize()
{
     static char *Name = "ELM_2NODE_LINE";

     element_type_t ElementDef;

     ElementDef.ElementName = Name;
     ElementDef.ElementCode = 202;

     ElementDef.NumberOfNodes = 2;

     ElementDef.NodeU = NodeU;
     ElementDef.NodeV = NULL;
     ElementDef.NodeW = NULL;

     ElementDef.PartialU = (double (*)())elm_2node_bar_dndu_fvalue;
     ElementDef.PartialV = (double (*)())NULL;
     ElementDef.PartialW = (double (*)())NULL;

     ElementDef.FunctionValue = (double (*)())elm_2node_bar_fvalue;
     ElementDef.Triangulate   = (int (*)())elm_2node_bar_triangulate;
     ElementDef.IsoLine       = (int (*)())NULL;
     ElementDef.PointInside   = (int (*)())NULL;
     ElementDef.IsoSurface    = (int (*)())NULL;

     return elm_add_element_type( &ElementDef ) ;
}
