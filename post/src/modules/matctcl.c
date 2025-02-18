#/*****************************************************************************
# *
# *  Elmer, A Finite Element Software for Multiphysical Problems
# *
# *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
# * 
# *  This program is free software; you can redistribute it and/or
# *  modify it under the terms of the GNU General Public License
# *  as published by the Free Software Foundation; either version 2
# *  of the License, or (at your option) any later version.
# * 
# *  This program is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with this program (in file fem/GPL-2); if not, write to the 
# *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
# *  Boston, MA 02110-1301, USA.
# *
# *****************************************************************************/

/*******************************************************************************
 *
 *    Elmerpost / MATC utilities.
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
 * $Id: matctcl.c,v 1.3 1999/06/03 14:12:40 jpr Exp $ 
 *
 * $Log: matctcl.c,v $
 * Revision 1.3  1999/06/03 14:12:40  jpr
 * *** empty log message ***
 *
 * Revision 1.2  1998/08/01 12:34:59  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "../elmerpost.h"
#include <tcl.h>
#include <tk.h>

extern Tcl_Interp *TCLInterp;

static VARIABLE *matc_tcl( VARIABLE *ptr )
{
   VARIABLE *res = NULL;
   char *command;

   int i,n;

   command = var_to_string(ptr);

   Tcl_GlobalEval( TCLInterp, command );

   FREEMEM( command );

   if ( TCLInterp->result && (n=strlen(TCLInterp->result))>0 )
   {
       res = var_temp_new( TYPE_STRING,1,n );
       for( i=0; i<n; i++ ) M( res,0,i ) = TCLInterp->result[i];
   }

   return res;
}

static VARIABLE *matc_element( VARIABLE *ptr )
{
    double *num = MATR(ptr);
    char *str = NULL;
   
    int i,j,n,maxn=0;

    VARIABLE *res = NULL;

    element_t *elem = CurrentObject->ElementModel->Elements;

    if ( NEXT(ptr) ) str = var_to_string( NEXT(ptr) );

    if ( CurrentObject->ElementModel->NofElements <= 0 ) error( "element: no elements present.\n" );

    for( i=0; i<NCOL(ptr); i++ )
    {
        n = num[i];
        if ( n <  0 || n >= CurrentObject->ElementModel->NofElements )
        {
            error( "element: Envalid element index: [%d].\n",n );
        }

        maxn = MAX( maxn,elem[n].ElementType->NumberOfNodes );
    }

    res = var_temp_new( TYPE_DOUBLE, NCOL(ptr), maxn );

    for( i=0; i<NCOL(ptr); i++ )
    {
        n = num[i];
        for( j=0; j<elem[n].ElementType->NumberOfNodes; j++ )
        {
            M(res,i,j) = elem[n].Topology[j];
        }
    }

    if ( str ) FREEMEM( str );

    return res;
}

int Matctcl_Init()
{
   extern VARIABLE *elm_gradient(), *elm_divergence(),*elm_rotor_3D(),*elm_rotor_2D();

#if 1
   com_init(
             "grad", FALSE, FALSE, elm_gradient, 1, 1,
             "r = grad(f): compute gradient of a scalar variable f.\n"
           );

   com_init(
             "div", FALSE, FALSE, elm_divergence, 1, 1,
             "r = div(v): compute divergence of a vector variable v.\n"
           );

   com_init(
              "curl3d", FALSE, FALSE, elm_rotor_3D, 1, 1,
              "r = curl3d(v): compute curl of a vector variable v (in 3D).\n"
           );
   com_init(
              "curl2d", FALSE, FALSE, elm_rotor_2D, 1, 1,
              "r = curl2d(v): compute curl of a vector variable v (in 2D).\n"
           );
#endif

   com_init(
             "tcl", FALSE, FALSE, matc_tcl, 1, 1,
             "str = tcl(str): execute a command in tcl part of ElmerPost.\n"
           );

   com_init(
             "element", FALSE, FALSE, matc_element, 1, 2,
             "r = element(v): return element topology for elements given in vector v.\n"
           );

   return TCL_OK;
}
