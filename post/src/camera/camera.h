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
 * Camera main module includes
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


/*
 * $Id: camera.h,v 1.1.1.1 2005/05/31 06:29:21 vierinen Exp $ 
 *
 * $Log: camera.h,v $
 * Revision 1.1.1.1  2005/05/31 06:29:21  vierinen
 * ads
 *
 * Revision 1.2  1998/08/01 12:34:10  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#ifdef MODULE_CAMERAS
#   define CAM_EXT
#else
#   define CAM_EXT extern
#endif

CAM_EXT int GlobalPass,GlobalMakeSecondPass;

#define CAM_MAX_OBJECTS_MASKED 100

typedef enum
{
   camera_proj_ortho, camera_proj_perspective
} camera_proj_t;

typedef struct camera_s
{
    struct camera_s *Next;
    char *Name;

    float ViewportLowX, ViewportHighX,
          ViewportLowY,ViewportHighY;

    float LookAtX,LookAtY,LookAtZ;
    float UpX,UpY,UpZ;
    float LookFromX,LookFromY,LookFromZ;
    float ClipNear,ClipFar;

    float FieldAngle;
    camera_proj_t ProjectionType;

    object_t *ObjMask[CAM_MAX_OBJECTS_MASKED];

    logical_t OnOff,DrawFrame;
} camera_t;

CAM_EXT camera_t *Camera;


void cam_set_viewport(camera_t *camera,double lx,double ly,double hx,double hy);
void cam_set_projection( camera_t *camera,camera_proj_t projection );
void cam_set_field_angle( camera_t *camera, double angle );
void cam_set_look_from(camera_t *camera,double x,double y,double z,int relative);
void cam_set_up_vector(camera_t *camera,double x,double y,double z );
void cam_set_clip(camera_t *camera,double n,double f);
void cam_set_look_to(camera_t *camera,double x,double y,double z,int relative);
void cam_set_onoff( camera_t *camera,int onoff);
camera_t *cam_add_camera( camera_t *camera,char *name );
void cam_delete_list( camera_t *camera );
int cam_display_list( camera_t *camera, object_t *object );
