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
 * Main include file of ElmerPost. Mainly includes other include files ;-)
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
 *                Modified by:
 *
 *       Date of modification:
 *
 ******************************************************************************/

#include "../config.h"

#if defined(WIN32) || defined(MINGW32)
#include <windows.h>
#endif

#include <stdlib.h>
#include <stdio.h>

/* #include <malloc.h> */
#include <math.h>

#include <sys/types.h>

#include <signal.h>

#include <elmer/matc.h>


#if defined(MINGW32) || defined(WIN32)

#include "tk/tk.h"

#endif


#ifdef MODULE_MAIN
#define EXT 
#else
#define EXT extern
#endif

#ifndef MIN
#define MIN(x,y) ( (x)>(y) ? (y) : (x) )
#endif

#ifndef MAX
#define MAX(x,y) ( (x)>(y) ? (x) : (y) )
#endif

#ifndef ABS
#define ABS(x) ( (x)>(0) ? (x) : (-(x)) )
#endif

#define FALSE 0
#define TRUE  1

#ifndef DBL_MAX
#define DBL_MAX            1.79769313486231570e+308
#endif

#ifndef M_PI
#define M_PI (3.1415926535897931)
#endif

typedef unsigned char logical_t;

typedef struct
{
    char *name;
    double *f;
    double min,max;
} scalar_t;

typedef struct
{
    char *name;
    double *f;
    double min[3],max[3];
} vector_t;

typedef struct
{
    int VolumeSides;
    int VolumeEdges;
    int SurfaceSides;
    int StereoMode;
    int OutputPS, FitToPagePS;
    double StereoTran,StereoRot;
} global_options_t;

#ifdef MODULE_MAIN
  global_options_t GlobalOptions = { FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,0.03,5.00 };
#else
  extern global_options_t GlobalOptions;
#endif

extern double RealTime(),CPUTime();

#include "geometry.h"
#include "elements/elements.h"
#include "graphics/graphics.h"
#include "visuals/visual.h"
#include "objects/objects.h"
#include "camera/camera.h"

EXT unsigned int epMouseDown,epMouseDownTakesTooLong;
EXT double GraphicsAspect;
EXT unsigned int GraphicsXSize,GraphicsYSize;

#ifndef WIN32
EXT XFontStruct *CurrentXFont;
#endif

EXT int BreakLoop;

#ifdef MODULE_MAIN
  void (*user_hook_before_all)()    = NULL;
  void (*user_hook_after_all)()     = NULL;

  void (*user_hook_camera_before)() = NULL;
  void (*user_hook_camera_after)()  = NULL;

  void (*user_hook_object_before)() = NULL;
  void (*user_hook_object_after)()  = NULL;
#else
  extern void (*user_hook_before_all)();
  extern void (*user_hook_after_all)();

  extern void (*user_hook_camera_before)();
  extern void (*user_hook_camera_after)();

  extern void (*user_hook_object_before)();
  extern void (*user_hook_object_after)();
#endif
