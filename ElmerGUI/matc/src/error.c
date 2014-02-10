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
 *     Error/signal trapping for MATC.
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
 *                       Date: 30 May 1996
 *
 *                Modified by:
 *
 *       Date of modification:
 *
 ******************************************************************************/

/*
 * $Id: error.c,v 1.3 2007/06/08 08:12:17 jpr Exp $ 
 *
 * $Log: error.c,v $
 * Revision 1.3  2007/06/08 08:12:17  jpr
 * *** empty log message ***
 *
 * Revision 1.2  2005/05/27 12:26:19  vierinen
 * changed header install location
 *
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:34  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"
#include "str.h"

void sig_trap(int sig)
/*======================================================================
?  Interrupt or floating point exeption. Free all memory allocated 
|  after last call to doread, give an error message and
|  longjump back to doread.
&  longjmp, mem_free_all
~  doread
^=====================================================================*/
{
  fprintf( math_out, "^C\n" );

#if 0
  signal(SIGINT, sig_trap);
  signal(SIGFPE, sig_trap);
#endif

  mem_free_all();

  longjmp(*jmpbuf, 2);
}

void error_old(char *fmt,void *p1, void *p2)
/*======================================================================
?  An error is detected by the program. Free all memory allocated 
|  after last call to doread, give an error message and
|  longjump back to doread.
&  longjmp, fprintf, mem_free_all
~  doread
^=====================================================================*/
{
  PrintOut( "MATC ERROR: " );
  PrintOut( fmt,p1,p2 );

  mem_free_all();

  longjmp(*jmpbuf, 2);
}

