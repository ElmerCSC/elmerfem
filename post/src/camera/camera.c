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
 * Camera main module & utilities
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
 *                       Date: 3 Oct 1995
 *
 * Modification history:
 *
 ******************************************************************************/

extern void *TCLInterp;


/*
 * $Id: camera.c,v 1.4 2004/11/29 08:27:09 jpr Exp $ 
 *
 * $Log: camera.c,v $
 * Revision 1.4  2004/11/29 08:27:09  jpr
 * *** empty log message ***
 *
 * Revision 1.3  2003/02/06 09:37:46  jpr
 * *** empty log message ***
 *
 * Revision 1.2  1998/08/01 12:34:09  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#define MODULE_CAMERAS

#include "../elmerpost.h"

#ifdef NOTDEF
void cam_set_frame(int ncam,int status)
{
    camera[--ncam].frame = status;
}

void cam_set_camera_obj_mask(int ncam,int obj,int status)
{
    --ncam;
    camera[ncam].obj_mask[obj] = status;
}
#endif

/*******************************************************************************
 *
 *     Name:        cam_set_viewport( camera_t *,double,double,double,double )
 *
 *     Purpose:     set viewport values given camera structure pointer and
 *                  two points.
 *
 *     Parameters:
 *
 *         Input:   (double,double) lower left hand corner of the viewport
 *                  (double,double) upper right hand corner of the viewport
 *
 *         Output:  (camera_t *) camera strucrure is modified
 *
 *   Return value:  void
 *
 *******************************************************************************/
void cam_set_viewport(camera_t *camera,double lx,double ly,double hx,double hy)
{
    camera->ViewportLowX  = lx;
    camera->ViewportLowY  = ly;
    camera->ViewportHighX = hx;
    camera->ViewportHighY = hy;
}

/*******************************************************************************
 *
 *     Name:        cam_set_projection( camera_t *,camera_proj_t )
 *
 *     Purpose:     Set camera projection type. Internal only.
 *
 *     Parameters:
 *
 *         Input:   (camera_proj_t) projection type to set.
 *
 *         Output:  (camera_t *) camera strucrure is modified
 *
 *   Return value:  void
 *
 *******************************************************************************/
void cam_set_projection( camera_t *camera,camera_proj_t projection )
{
    camera->ProjectionType = projection;
}

/*******************************************************************************
 *
 *     Name:        cam_set_file_angle( camera_t *,float )
 *
 *     Purpose:     Set camera filed angle for perspective projection.
 *                  Internal only.
 *
 *     Parameters:
 *
 *         Input:   (float) intput field angle.
 *
 *         Output:  (camera_t *) camera strucrure is modified
 *
 *   Return value:  void
 *
 *******************************************************************************/
void cam_set_field_angle( camera_t *camera, double angle )
{
    camera->FieldAngle = angle;
}

/*******************************************************************************
 *
 *     Name:        cam_set_look_from( camera_t *,double,double,double,int )
 *
 *     Purpose:     Set camera position to given value.
 *
 *     Parameters:
 *
 *         Input:   (double,double,double) x,y,z coordinates of the camera
 *                  (int) flag saying if inputs should be added to current
 *                        position
 *
 *         Output:  (camera_t *) camera strucrure is modified
 *
 *   Return value:  void
 *
 *******************************************************************************/
void cam_set_look_from(camera_t *camera,double x,double y,double z,int relative)
{
    if ( relative )
    {
        camera->LookFromX += x;
        camera->LookFromY += y;
        camera->LookFromZ += z;
    } else
    {
        camera->LookFromX = x;
        camera->LookFromY = y;
        camera->LookFromZ = z;
    }
}

/*******************************************************************************
 *
 *     Name:        cam_set_up_vector( camera_t *,double,double,double )
 *
 *     Purpose:     Set camera up vector
 *
 *     Parameters:
 *
 *         Input:   (double,double,double) x,y,z direction upwadrds
 *
 *         Output:  (camera_t *) camera strucrure is modified
 *
 *   Return value:  void
 *
 *******************************************************************************/
void cam_set_up_vector(camera_t *camera,double x,double y,double z )
{
    camera->UpX = x;
    camera->UpY = y;
    camera->UpZ = z;
}

/*******************************************************************************
 *
 *     Name:        cam_set_clip( camera_t *,double,double )
 *
 *     Purpose:     Set camera clip viewing planes
 *
 *     Parameters:
 *
 *         Input:   (double,double) near and far plane coordinates
 *
 *         Output:  (camera_t *) camera strucrure is modified
 *
 *   Return value:  void
 *
 *******************************************************************************/
void cam_set_clip(camera_t *camera,double n,double f)
{
    camera->ClipNear = n;
    camera->ClipFar  = f;
}

/*******************************************************************************
 *
 *     Name:        cam_set_look_to( camera_t *,double,double,double,int )
 *
 *     Purpose:     Set point at which the camera is aimed
 *
 *     Parameters:
 *
 *         Input:   (double,double,double) x,y,z coordinates
 *                  (int) flag saying if inputs should be added to current
 *                        position
 *
 *         Output:  (camera_t *) camera strucrure is modified
 *
 *   Return value:  void
 *
 *******************************************************************************/
void cam_set_look_to(camera_t *camera,double x,double y,double z,int relative)
{
    if ( relative )
    {
        camera->LookAtX += x;
        camera->LookAtY += y;
        camera->LookAtZ += z;
    } else
    {
        camera->LookAtX = x;
        camera->LookAtY = y;
        camera->LookAtZ = z;
    }
}

void cam_set_onoff( camera_t *camera,int onoff)
{
    camera->OnOff = onoff;
}

/*******************************************************************************
 *
 *     Name:        cam_add_camera( camera_t *, char * )
 *
 *     Purpose:     Add camera to a given list of camera with given name.
 *                  If one with given name already exists, it is returned.
 *
 *     Parameters:
 *
 *         Input:   (camera_t *) input list of cameras
 *                  (char *)     name to be given to the camera
 *
 *         Output:  (camera_t *) camera strucrure is modified
 *
 *   Return value:  (camera_t *) pointer to the crated camera structure
 *
 *******************************************************************************/
camera_t *cam_add_camera( camera_t *camera,char *name )
{
    camera_t *new_cam = camera;


    while( new_cam )
    {
		if ( strcmp( new_cam->Name,name ) == 0 ) return new_cam;
		new_cam = new_cam->Next;
    }

	new_cam = (camera_t *)calloc(1,sizeof(camera_t) );

    if ( !new_cam )
    {
        fprintf( stderr, "cam_add_camera: FATAL: can't alloc a few bytes of memory\n" );
        return NULL;
    }

    if ( !(new_cam->Name = (char *)malloc( strlen(name)+1 ) ) ) 
    {
        fprintf( stderr, "cam_add_camera: FATAL: can't alloc a few bytes of memory\n" );
        free( new_cam );
        return NULL;
    }

    strcpy( new_cam->Name, name );

    if ( camera )
    {
        while( camera->Next ) camera = camera->Next;
        camera->Next = new_cam;
    }

    return new_cam;
}

/*******************************************************************************
 *
 *     Name:        cam_delete_list( camera_t * )
 *
 *     Purpose:     Delete (free mem ) list of camera definitions
 *
 *     Parameters:
 *
 *         Input:   (camera_t *) input list of cameras
 *
 *         Output:  none
 *
 *   Return value:  void
 *
 *******************************************************************************/
void cam_delete_list( camera_t *camera )
{
    camera_t *ptr;

    while( camera )
    {
        if ( camera->Name ) free( camera->Name );

        ptr    = camera;
        camera = camera->Next;

        free( ptr );
    }
}

static double upx,upy,upz,viewx,viewy,viewz,tox,toy,toz;

#include <tcl.h>
/*******************************************************************************
 *
 *     Name:        cam_display_list( camera_t *, object_t * )
 *
 *     Purpose:     Display list of given objecst is given list of cameras
 *
 *     Parameters:
 *
 *         Input:   (camera_t *) input list of cameras
 *                  (object_t *) input list of objecst
 *
 *         Output:  graphics
 *
 *   Return value:  if mouse interaction is going on and too slow FALSE,
 *                  otherwise TRUE
 *
 *******************************************************************************/
int cam_display_list( camera_t *camera, object_t *object )
{
    double t = RealTime(), ct = CPUTime();

    int FitToPage = 0, nofcameras;
    camera_t *cam;

    if ( GlobalOptions.OutputPS ) {
       initglp( Tcl_GetVar( TCLInterp, "PSFileName", TCL_GLOBAL_ONLY ), 
                GlobalOptions.FitToPagePS );
    }
    if ( user_hook_before_all ) (*user_hook_before_all)( camera,object );

     nofcameras = 0;
     for( cam=camera; cam != NULL; cam = cam->Next, nofcameras++ );

    for( GlobalPass=0; GlobalPass < 2; GlobalPass++ )
    {
        for( cam=camera; cam != NULL; cam = cam->Next )
        {
            if ( !cam->OnOff ) continue;

            gra_set_projection( cam->ProjectionType, cam->FieldAngle,
                                cam->ViewportLowX, cam->ViewportHighX,
                                cam->ViewportLowY, cam->ViewportHighY,
	       	 	        cam->ClipNear, cam->ClipFar, nofcameras>1 );

            gra_push_matrix();

            gra_look_at(
                         cam->LookFromX, cam->LookFromY, cam->LookFromZ,
                            cam->LookAtX, cam->LookAtY, cam->LookAtZ,
                                  cam->UpX, cam->UpY, cam->UpZ
                       );

            if ( user_hook_camera_before ) (*user_hook_camera_before)( GlobalPass,cam,object,t );

            if ( !obj_display_list( object, t ) ) return FALSE;

            if ( user_hook_camera_after ) (*user_hook_camera_after)( GlobalPass,cam,object,t );

            gra_pop_matrix();

             if ( BreakLoop ) break;

        }
        if ( BreakLoop ) break;
    } 

    if ( user_hook_after_all ) (*user_hook_after_all)( camera,object );
    if ( GlobalOptions.OutputPS ) stopglp();

    return TRUE;
}
