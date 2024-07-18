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
 * Visual classes main module + utilities
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
 *                       Date: 26 Sep 1995
 *
 *
 * Modification history:
 *
 * 28 Sep 1995, - added routine vis_get_visual_type, which gives a pointer to 
 *                visual_type_t structure given name of the visual type
 *
 *              - modified vis_add_visual_type according to change in 
 *                visual_type_t and visual_t structures (now holding list of
 *                visual types rather than an array)
 *
 * 29 Sep 1995, - added routines vis_get_visual_type_name and vis_get_visual_name
 *              - added routines vis_set_param and vis_new_visual
 *              - added routines vis_add_visual, vis_delete_visual
 *
 * Juha R.
 *
 ******************************************************************************/

#define MODULE_VISUALS

#include "../elmerpost.h"

/*******************************************************************************
 *
 *     Name:         vis_add_visual_type
 *
 *     Purpose:      Register a visual class
 *
 *     Parameters:
 *
 *         Input:    (visual_t *)  visual class to be added
 *
 *         Output:   Global variable VisualDefs is modified
 *
 *   Return value:   TRUE is success, FALSE if malloc() fails
 *
 ******************************************************************************/
int vis_add_visual_type( visual_type_t *VisualDef )
{
    visual_type_t *ptr;

    ptr = (visual_type_t *)calloc(sizeof(visual_type_t),1);
    if ( !ptr )
    {
        fprintf( stderr, "FATAL: vis_visual_type_add: Can't allocate memory.\n" );
        return FALSE;
    }

    *ptr = *VisualDef;
    ptr->Next = VisualDefs.VisualTypes;

    VisualDefs.VisualTypes = ptr;

    VisualDefs.NumberOfVisualTypes++;

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:         vis_get_visual_type
 *
 *     Purpose:      return visual type pointer given visual name
 *
 *     Parameters:
 *
 *         Input:    (char *) name of visual
 *
 *         Output:   none
 *
 *   Return value:   (visual_type_t *)visual if found, NULL otherwise
 *
 ******************************************************************************/
visual_type_t *vis_get_visual_type(char *name)
{
    visual_type_t *type = VisualDefs.VisualTypes;

    for( type=VisualDefs.VisualTypes; type != NULL; type=type->Next )
    {
        if ( strcmp( name, type->VisualName ) == 0 ) return type;
    }

    fprintf( stderr, "vis_get_visual_type: can't find visual type [%s]\n", name );

    return NULL;
}

/*******************************************************************************
 *
 *     Name:         vis_get_visual_type_name
 *
 *     Purpose:      return pointer to visual type name given visual structure
 *
 *     Parameters:
 *
 *         Input:    (visual_t *) visual structure
 *
 *         Output:   none
 *
 *   Return value:   (char *) name of the visual type
 *
 ******************************************************************************/
char *vis_get_visual_type_name(visual_t *visual)
{
    return visual->VisualType->VisualName;
}

/*******************************************************************************
 *
 *     Name:         vis_get_visual_name
 *
 *     Purpose:      return pointer to visual name given visual structure
 *
 *     Parameters:
 *
 *         Input:    (visual_t *) visual structure
 *
 *         Output:   none
 *
 *   Return value:   (char *) name of the visual
 *
 ******************************************************************************/
char *vis_get_visual_name(visual_t *visual)
{
    return visual->Name;
}

/*******************************************************************************
 *
 *     Name:        vis_set_param
 *
 *     Purpose:     set visual parameters
 *
 *     Parameters: 
 *
 *         Input:   (visual_t *) visual whose parameter is set
 *                  (char *)     parameter name
 *                  (intvalue,doublevalue,void *) pointer to param value 
 *
 *         Output:  graphics
 *   
 *   Return value:  success
 *
 ******************************************************************************/
int vis_set_param
  ( visual_t *visual, char *name,int intvalue,double doublevalue,void *ptrvalue )
{
    visual_type_t  *VisualType   = visual->VisualType;
    visual_param_t *VisualParams = VisualType->VisualParams;

    char *offset = (void *)visual->VisualParams;

    if ( !offset || !name )
    {
       fprintf( stderr, "vis_set_param: argument failure.\n" );
       return FALSE;
    }

    for( ; VisualParams->Name != NULL; VisualParams++ )
    {
        if ( strcmp( name, VisualParams->Name ) == 0 )
        {
            if ( VisualParams->ParamType == VIS_VISUAL_PARAM_LOGICAL )
            {
                *(logical_t *)(offset + VisualParams->Offset) = intvalue;
            }
            else if ( VisualParams->ParamType == VIS_VISUAL_PARAM_INT )
            {
                *(int *)(offset + VisualParams->Offset) = intvalue;
            }
            else if ( VisualParams->ParamType == VIS_VISUAL_PARAM_FLOAT )
            {
                *(double *)(offset + VisualParams->Offset) = doublevalue;
            }
            else if ( VisualParams->ParamType == VIS_VISUAL_PARAM_POINTER )
            {
                *(void **)(offset + VisualParams->Offset) = (void *)ptrvalue;
            }
            break;
        }
    }

    if ( !VisualParams->Name )
    {
        fprintf( stderr, "vis_set_param: no such param [%s,%s]\n", VisualType->VisualName,name );
        return FALSE;
    }

    return TRUE;
}

/*******************************************************************************
 *
 *     Name:        vis_default_params
 *
 *     Purpose:     set default parameter values for a visual
 *
 *     Parameters: 
 *
 *         Input:   (visual_t *) visual structure to modify
 *
 *         Output:  
 *   
 *   Return value:  void
 *
 ******************************************************************************/
void vis_default_params( visual_t *visual )
{
    visual_type_t  *VisualType   = visual->VisualType;
    visual_param_t *VisualParams = VisualType->VisualParams;

    for( ;VisualParams->Name != NULL; VisualParams++ )
    {
        vis_set_param(  visual,VisualParams->Name, VisualParams->IntValue,
             VisualParams->FloatValue, VisualParams->PointerValue );
    }
}

/*******************************************************************************
 *
 *     Name:        vis_new_visual
 *
 *     Purpose:     alloc memory for a new visual structure
 *
 *     Parameters: 
 *
 *         Input:   (char *)  Name of the visual type
 *
 *         Output:  none
 *   
 *   Return value:  pointer to visual_t or NULL if malloc fails
 *
 ******************************************************************************/
visual_t *vis_new_visual(char *name)
{
    visual_t *visual = (visual_t *)calloc( sizeof(visual_t),1 );

    if ( !visual )
    {
       fprintf( stderr, "vis_new_visual: FATAL: can't allocate (a few bytes of) memory\n" );
       return NULL;
    }

    if ( !(visual->VisualType = vis_get_visual_type(name)) )
    {
        free( visual );
        return NULL;
    }

    if ( !(visual->VisualParams = (*visual->VisualType->AllocParams)()) )
    {
        free( visual );
        return NULL;
    }

    vis_default_params( visual );

    return visual;
}

/*******************************************************************************
 *
 *     Name:        vis_link_visual
 *
 *     Purpose:     link a new visual to list of visuals given
 *
 *     Parameters: 
 *
 *         Input:   (visual_t *) list of visuals (can be null)
 *                  (visual_t *) visual to be added    
 *
 *         Output:  (visual_t *) is modified
 *   
 *   Return value:  (visual_t *) head of list
 *
 ******************************************************************************/
visual_t *vis_link_visual( visual_t *list,visual_t *new )
{
    visual_t *ptr = list;

    if ( list )
    {
        while( ptr->Next != NULL ) ptr = ptr->Next; 
        ptr->Next = new;
    } else list = new;

    return list;
}

/*******************************************************************************
 *
 *     Name:        vis_add_visual
 *
 *     Purpose:     add a new visual to list of visuals given
 *
 *     Parameters: 
 *
 *         Input:   (visual_t *) list of visuals (can be null)
 *                  (char *)     name of visual type to add
 *
 *         Output:  (visual_t *) is modified
 *   
 *   Return value:  pointer to new visual_t or NULL if malloc fails
 *
 ******************************************************************************/
visual_t *vis_add_visual(visual_t *visual,char *name)
{
    visual_t *newvisual = (visual_t *)vis_new_visual(name);

    if ( !newvisual ) return NULL;
    
    if ( visual )
    {
        while( visual->Next != NULL ) visual = visual->Next; 
        visual->Next = newvisual;
    }

    return newvisual;
}

/*******************************************************************************
 *
 *     Name:        vis_delete_visual
 *
 *     Purpose:     delete list of visuals given as argument
 *
 *     Parameters: 
 *
 *         Input:   (visual_t *) list to be deleted
 *
 *         Output:  none
 *   
 *   Return value:  void
 *
 ******************************************************************************/
void vis_delete_visual( visual_t *visual )
{
    visual_t *ptr;

    while( visual != NULL )
    {
        (*visual->VisualType->DeleteParams)( visual->VisualParams );
        visual->VisualParams = NULL;

        ptr = visual->Next; 

        free( visual );
        visual = ptr;
    }
}

/*******************************************************************************
 *
 *     Name:        vis_next_visual
 *
 *     Purpose:     return pointer to next visual in a list of given visual
 *
 *     Parameters: 
 *
 *         Input:   (visual_t *)
 *
 *         Output:  none
 *   
 *   Return value:  pointer to visual_t or NULL if no more visuals in the list
 *
 ******************************************************************************/
visual_t *vis_next_visual( visual_t *visual )
{
    if ( visual ) return visual->Next;

    return NULL;
}

/*******************************************************************************
 *
 *     Name:         vis_display_visual
 *
 *     Purpose:      Display one visual_t visual
 *
 *     Parameters:
 *
 *         Input:    (geometry_t *) geometry information
 *                   (visual_t *)   visual
 *
 *         Output:   graphics
 *
 *   Return value:   if mouse interaction is going on and too slow FALSE,
 *                   otherwise true
 *
 ******************************************************************************/
int vis_display_visual( geometry_t *geometry, element_model_t *model, visual_t *VL,double t )
{
    return (*VL->VisualType->RealizeVisual)( geometry, model, VL->VisualParams, t );
}

/*******************************************************************************
 *
 *     Name:         vis_display_list
 *
 *     Purpose:      Display a list of visual_t visuals
 *
 *     Parameters:
 *
 *         Input:    (geometry_t *) geometry information
 *                   (visual_t *)   visuals
 *                   (double t)     real time at invocation
 *
 *         Output:   graphics
 *
 *   Return value:   if mouse interaction is going on and too slow FALSE,
 *                   otherwise true
 *
 ******************************************************************************/
int vis_display_list( geometry_t *geometry, element_model_t *model, visual_t *VL,double t )
{
    for( ; VL != NULL; VL = VL->Next )
    {
#ifdef DEBUG
fprintf( stderr, "DISPLAYING VISUAL TYPE: [%s]\n", VL->VisualType->VisualName );
#endif
        if ( !(*VL->VisualType->RealizeVisual)( geometry, model, VL->VisualParams, t ) )
        {
            return FALSE;
        }

        if ( BreakLoop ) break;
    }

    return TRUE;
}


/*******************************************************************************
 *
 *     Name:         vis_initialize_visual_types
 *
 *     Purpose:      Register all internal visual classes
 *
 *     Parameters:
 *
 *         Input:    none
 *
 *         Output:   Global variable VisualDefs is modified
 *
 *   Return value:   TRUE is success, FALSE if malloc() fails
 *
 ******************************************************************************/
int vis_initialize_visual_types()
{
     if ( !vis_initialize_mesh_visual() ) return FALSE;
     if ( !vis_initialize_arrow_visual() ) return FALSE;
     if ( !vis_initialize_sphere_visual() ) return FALSE;
     if ( !vis_initialize_contour_line_visual() ) return FALSE;
     if ( !vis_initialize_isosurface_visual() ) return FALSE;
     if ( !vis_initialize_particle_visual() ) return FALSE;
     if ( !vis_initialize_colscale_visual() ) return FALSE;

     return TRUE;
}
