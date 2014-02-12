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
 * Main module for element model descriptions & utility routines.
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
 * 28 Sep 1995, removed elm_triangle_normal 
 * Juha R
 *
 *
 ******************************************************************************/

#define MODULE_ELEMENTS

#include "../elmerpost.h"
#include <elements.h>

/*******************************************************************************
 *
 *     Name:        elm_element_type_add
 *
 *     Purpose:     Add an element to element type list
 *
 *     Parameters:
 *
 *         Input:   (element_type_t *) structure describing the element type
 *
 *         Output:  Global list of element types is modified
 *
 *   Return value:  malloc() success
 *
 ******************************************************************************/
int elm_add_element_type( element_type_t *Def )
{
    element_type_t *ptr;

    ptr = (element_type_t *)calloc( sizeof(element_type_t),1 );

    if ( !ptr )
    {
        fprintf( stderr, "elm_add_element_type: FATAL: can't allocate (a few bytes of) memory.\n" );
        return FALSE;
    }

    *ptr = *Def;

    ptr->Next = ElementDefs.ElementTypes;
    ElementDefs.ElementTypes = ptr;

    return TRUE;
}

int elm_2node_bar_initialize();
int elm_3node_bar_initialize();
int elm_4node_bar_initialize();

int elm_3node_triangle_initialize();
int elm_4node_triangle_initialize();
int elm_6node_triangle_initialize();
int elm_10node_triangle_initialize();

int elm_4node_quad_initialize();
int elm_5node_quad_initialize();
int elm_8node_quad_initialize();
int elm_9node_quad_initialize();
int elm_12node_quad_initialize();
int elm_16node_quad_initialize();
int elm_4node_tetra_initialize();
int elm_8node_tetra_initialize();
int elm_10node_tetra_initialize();

int elm_5node_pyramid_initialize();
int elm_13node_pyramid_initialize();

int elm_6node_wedge_initialize();
int elm_15node_wedge_initialize();

int elm_8node_brick_initialize();
int elm_20node_brick_initialize();
int elm_27node_brick_initialize();

/*******************************************************************************
 *
 *     Name:        elm_initialize_element_types
 *
 *     Purpose:     Initialize all internal element types
 *
 *     Parameters:
 *
 *         Input:   none
 *
 *         Output:  Global list of element types is modified
 *
 *   Return value:  malloc() success
 *
 ******************************************************************************/
int elm_initialize_element_types()
{
    if ( !elm_2node_bar_initialize() )       return FALSE;
    if ( !elm_3node_bar_initialize() )       return FALSE;
    if ( !elm_4node_bar_initialize() )       return FALSE;

    if ( !elm_3node_triangle_initialize() )  return FALSE;
    if ( !elm_4node_triangle_initialize() )  return FALSE;
    if ( !elm_6node_triangle_initialize() )  return FALSE;
    if ( !elm_10node_triangle_initialize() ) return FALSE;

    if ( !elm_4node_quad_initialize() )     return FALSE;
    if ( !elm_5node_quad_initialize() )     return FALSE;
    if ( !elm_8node_quad_initialize() )     return FALSE;
    if ( !elm_9node_quad_initialize() )     return FALSE;
    if ( !elm_12node_quad_initialize() )    return FALSE;
    if ( !elm_16node_quad_initialize() )    return FALSE;

    if ( !elm_4node_tetra_initialize() )    return FALSE;
    if ( !elm_8node_tetra_initialize() )    return FALSE;
    if ( !elm_10node_tetra_initialize() )   return FALSE;

    if ( !elm_5node_pyramid_initialize() )    return FALSE;
    if ( !elm_13node_pyramid_initialize() )    return FALSE;

    if ( !elm_6node_wedge_initialize() )    return FALSE;
    if ( !elm_15node_wedge_initialize() )    return FALSE;

    if ( !elm_8node_brick_initialize() )    return FALSE;
    if ( !elm_20node_brick_initialize() )   return FALSE;
    if ( !elm_27node_brick_initialize() )   return FALSE;

    return TRUE;
}

void elm_force_load(int stat)
{
    if ( stat ) {
      elm_divergence();
      elm_gradient();
      elm_rotor_3D();
    }
}
