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
 *     The character of the routines in this file.
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
 *                       Date: 1 Oct 1995
 *
 * Modification history:
 *
 ******************************************************************************/

#ifdef MODULE_OBJECTS
#   define OBJ_EXT
#else
#   define OBJ_EXT extern
#endif

typedef double matrix_t[4][4];

typedef enum
{
    rot_pri_xyz,
    rot_pri_xzy,
    rot_pri_yxz,
    rot_pri_yzx,
    rot_pri_zxy,
    rot_pri_zyx,
    rot_pri_local,
    rot_pri_parent
} rot_pri_t;

typedef enum
{
    trn_pri_trs,
    trn_pri_tsr,
    trn_pri_rts,
    trn_pri_rst,
    trn_pri_str,
    trn_pri_srt
} trn_pri_t;

typedef struct transform_s
{
    struct transform_list_s *Children;
    struct transform_s *Parent;

    matrix_t Matrix;

    matrix_t RotMatrix;
    matrix_t TrnMatrix;
    matrix_t SclMatrix;

    double RotX,RotY,RotZ;
    double TrnX,TrnY,TrnZ;
    double SclX,SclY,SclZ;

    rot_pri_t RotationPriority;
    trn_pri_t TransformPriority;
} transform_t;

typedef struct transform_list_s
{
    struct transform_list_s *Next;
    transform_t *Entry;
} transform_list_t;

typedef struct object_s
{
    struct object_s *Next;
    int Id;
    char *Name;

    transform_t Transform;

    element_model_t *ElementModel;

    geometry_t *Geometry;

    int ClipPlane[6];
    double ClipEquation[6][4];

    struct visual_s *VisualList;
} object_t;

OBJ_EXT object_t VisualObject,RotObject,*CurrentObject;
OBJ_EXT double PiDiv180;

void obj_object_initialize( object_t *object );
object_t *obj_new(char *name);
object_t *obj_add_object( object_t *object,char *name );
object_t *obj_find( object_t *object,char *name );
int obj_display_list( object_t *object,double t );
void obj_init_transform( transform_t *transform );
void obj_get_matrix(matrix_t matrix,object_t *object);
void obj_get_matrix_transpose( matrix_t matrix,object_t *object );
void obj_set_rotation_priority( object_t *object, rot_pri_t priority);
void obj_set_transform_priority( object_t *object, trn_pri_t priority);
int obj_set_parent( object_t *object,object_t *parent );
void obj_rotate( object_t *object,double x,double y,double z,int which,int relative );
void obj_scale( object_t *object,double x,double y,double z,int which, int relative );
void obj_translate( object_t *object,double x,double y,double z,int which,int relative );
void obj_set_matrix( object_t *object );
