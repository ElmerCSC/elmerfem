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
 *     Object transformation routines.
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
 *                       Date: 6 Jun 1996
 *
 *                Modified by:
 *
 *       Date of modification:
 *
 ******************************************************************************/


/*
 * $Id: transforms.c,v 1.3 1999/06/03 14:12:40 jpr Exp $ 
 *
 * $Log: transforms.c,v $
 * Revision 1.3  1999/06/03 14:12:40  jpr
 * *** empty log message ***
 *
 * Revision 1.2  1998/08/01 12:35:03  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "../elmerpost.h"

#include <tcl.h>
#include <tk.h>


static int TrnPriority(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    static int first=TRUE,n,prior;
    double x=0,y=0,z=0;

    if ( argc != 2 )
    {
        sprintf( interp->result, "cTrnPriority: wrong number of parameters\n" );
        return TCL_ERROR;
    }

    n = sscanf( argv[1], "%d", &prior );

    if ( n <= 0 )
    {
        sprintf( interp->result, "cTrnPriority: wrong number of parameters\n" );
        return TCL_ERROR;
    }

    obj_set_transform_priority( CurrentObject,prior );

    return TCL_OK;
}

static int RotPriority(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    static int first=TRUE,n,prior;
    double x=0,y=0,z=0;

    if ( argc != 2 )
    {
        sprintf( interp->result, "cRotPrioryty: wrong number of parameters\n" );
        return TCL_ERROR;
    }

    n = sscanf( argv[1], "%d", &prior );

    if ( n <= 0 )
    {
        sprintf( interp->result, "cRotPriority: wrong number of parameters\n" );
        return TCL_ERROR;
    }

    if ( prior > 7 )
    {
        if ( CurrentObject == &VisualObject )
            CurrentObject = VisualObject.Next;
        else
            CurrentObject = &VisualObject;

        return TCL_OK;
    }

    obj_set_rotation_priority( CurrentObject,prior );

    return TCL_OK;
}

static int Rotate(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    static int first=TRUE,n;
    double x=0,y=0,z=0;

    if ( argc != 4 )
    {
        sprintf( interp->result, "cRotate: wrong number of parameters\n" );
        return TCL_ERROR;
    }

    n = sscanf( argv[3], "%lf %lf %lf", &x,&y,&z );

    if ( n <= 0 )
    {
        sprintf( interp->result, "cRotate: wrong number of parameters\n" );
        return TCL_ERROR;
    }

    obj_rotate( CurrentObject,x,y,z,argv[1][0],argv[2][0]=='r' );

    opengl_draw();

    return TCL_OK;
}

static int Scale(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    static int n,first=TRUE;
    double x,y,z;

    if ( argc != 4 )
    {
        sprintf( interp->result, "cScale: wrong number of parameters\n" );
        return TCL_ERROR;
    }

    if ( (n=sscanf( argv[3], "%lf %lf %lf", &x,&y,&z )) < 3 ) { y = x; z = x; }

    if ( n <= 0 )
    {
        sprintf( interp->result, "cScale: wrong number of parameters\n" );
        return TCL_ERROR;
    }

    obj_scale( CurrentObject,x,y,z,argv[1][0],argv[2][0]=='r' );

    opengl_draw();

    return TCL_OK;
}

static int Translate(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
    static int n,first=TRUE;
    double x=0,y=0,z=0;

    if ( argc != 4 )
    {
        sprintf( interp->result, "cScale: wrong number of parameters\n" );
        return TCL_ERROR;
    }

    if ( (n = sscanf( argv[3], "%lf %lf %lf", &x,&y,&z ) ) <= 0 )
    {
        sprintf( interp->result, "cScale: wrong number of parameters\n" );
        return TCL_ERROR;
    }

    obj_translate( CurrentObject,x,y,z,argv[1][0],argv[2][0]=='r' );

    opengl_draw();

    return TCL_OK;
}

int Transforms_Init( Tcl_Interp *interp )
{
    Tcl_CreateCommand( interp,"cTranslate",Translate,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand( interp,"cScale",Scale,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand( interp,"cRotate",Rotate,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand( interp,"cRotPriority",RotPriority,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand( interp,"cTrnPriority",TrnPriority,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);

    return TCL_OK;
}
