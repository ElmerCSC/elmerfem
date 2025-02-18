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
 * Objects main module & utilities
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
 *                       Date: 2 Oct 1995
 *
 *                Modified by:
 *
 *       Date of modification:
 *
 ******************************************************************************/

#define MODULE_OBJECTS

#include "../elmerpost.h"

/*******************************************************************************
 *
 *     Name:          obj_initialize_object
 *
 *     Purpose:       Initialize object to default values
 *
 *     Parameters:
 *
 *         Input:    none
 *
 *         Output:   (object_t *)
 *
 *   Return value:   void
 *
 ******************************************************************************/
void obj_object_initialize( object_t *object )
{
	int i;

    obj_init_transform( &object->Transform );
	for( i=0; i<6;i++ ) object->ClipPlane[i] = -1;
    PiDiv180 = acos(0.0) / 90.0;
}

/*******************************************************************************
 *
 *     Name:          obj_new
 *
 *     Purpose:       Create a new object. Internal only.
 *
 *     Parameters:
 *
 *         Input:    Name of the object to create.
 *
 *         Output:   none
 *
 *   Return value:   (object_t *)
 *
 ******************************************************************************/
object_t *obj_new(char *name)
{
    object_t *object = calloc(1,sizeof(object_t));
    
    if ( !object ) 
    {
       fprintf( stderr, "object_new: FATAL: can't allocate a few bytes of memory.\n" );
       return NULL;
    }

    if ( name )
    {
        if ( !(object->Name = malloc( strlen(name)+1 ) ) )
        {
           fprintf( stderr, "object_new: FATAL: can't allocate a few bytes of memory.\n" );
           free( object );
           return NULL;
        }
 
        strcpy( object->Name,name );
    }

    obj_object_initialize( object );

    return object;
}

/*******************************************************************************
 *
 *     Name:          obj_add_object
 *
 *     Purpose:       Add a new object to list of objects given.
 *
 *     Parameters:
 *
 *         Input:    (object_t *) input list
 *                   (char *) name of the object to create
 *
 *         Output:   (object_t *) is modified
 *
 *   Return value:   pointer to (object_t *), the new entry in the list.
 *
 ******************************************************************************/
object_t *obj_add_object( object_t *object,char *name )
{
     object_t *new = obj_new(name);

     if ( !new ) return NULL;

     if ( object )
     {
         new->Id = 1;
         while( object->Next ) { object = object->Next; new->Id++; }
         object->Next = new;
     }

     return new;
}

/*******************************************************************************
 *
 *     Name:          obj_find
 *
 *     Purpose:       Return pointer to an object with given name if any.
 *
 *     Parameters:
 *
 *         Input:    (object_t *) input list of objects
 *                   (char *) name of the object to find 
 *
 *         Output:   none
 *
 *   Return value:   pointer to (object_t *) if found, NULL otherwise
 *
 ******************************************************************************/
object_t *obj_find( object_t *object,char *name )
{
     while( object )
     {
        if ( strcmp( object->Name, name ) == 0 ) return object;
        object = object->Next;
     }

     return NULL;
}

/*******************************************************************************
 *
 *     Name:          obj_display_list
 *
 *     Purpose:       Display list of objects given.
 *
 *     Parameters:
 *
 *         Input:    (object_t *) input list of objects
 *
 *         Output:   graphics
 *
 *   Return value:   if mouse interaction is going on and too slow FALSE,
 *                   otherwise true
 *
 ******************************************************************************/
int obj_display_list( object_t *object,double t )
{
    extern double XMin,XMax,YMin,YMax,ZMin,ZMax;
    int i;

    for( ; object != NULL; object = object->Next )
    {
        gra_push_matrix();
        obj_set_matrix( object );

        for( i=0; i<6; i++ )
  	   if ( object->ClipPlane[i] >= 0 )
		gra_clip_plane( object->ClipPlane[i],object->ClipEquation[i] );

        if ( user_hook_object_before )
            (*user_hook_object_before)
              (
                 GlobalPass,object->Geometry,object->ElementModel,object->VisualList,t
              );

        if ( epMouseDown && epMouseDownTakesTooLong > 3 )
        {
            gra_bbox( XMin,XMax,YMin,YMax,ZMin,ZMax );
        }
        else if ( !vis_display_list( object->Geometry,object->ElementModel,object->VisualList,t ) )
        {
            gra_pop_matrix();
            return FALSE;
        }

        if ( user_hook_object_after )
            (*user_hook_object_after)
              (
                 GlobalPass,object->Geometry,object->ElementModel,object->VisualList,t
              );

         for( i=0; i<6; i++ )
           if ( object->ClipPlane[i] >= 0 ) gra_disable_clip( object->ClipPlane[i] );


        gra_pop_matrix();
    }

    return TRUE;
}
