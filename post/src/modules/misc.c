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
 *     Misc graphics utilities.
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
 * $Id: misc.c,v 1.3 1999/06/04 15:13:20 jim Exp $ 
 *
 * $Log: misc.c,v $
 * Revision 1.3  1999/06/04 15:13:20  jim
 * *** empty log message ***
 *
 * Revision 1.2  1998/08/01 12:35:00  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "../elmerpost.h"

#include <tcl.h>
#include <tk.h>

static int LineSmooth(ClientData cl,Tcl_Interp *interp,int argc,char **argv)
{
   if ( argc>1 )
   {
     if ( strcmp( argv[1], "on" ) == 0 )
     {
        glEnable( GL_LINE_SMOOTH );
        glEnable( GL_BLEND );
     } else {
        glDisable( GL_LINE_SMOOTH );
        glDisable( GL_BLEND );
     }
   } else {
     glEnable( GL_LINE_SMOOTH );
     glEnable( GL_BLEND );
   }

   return TCL_OK;
}

static int GraphicsClear( ClientData cl,Tcl_Interp *interp,int argc,char **argv )
{
   extern double br,bg,bb;

   glClearColor( br,bg,bb,1.0 );
   glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

   return TCL_OK;
}

int Misc_Init( Tcl_Interp *interp )
{
   Tcl_CreateCommand( interp,"linesmooth",LineSmooth,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
   Tcl_CreateCommand( interp,"gclear",GraphicsClear,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);

   return TCL_OK;
}
