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
 *     Random number generator.
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
/***********************************************************************
|
|  URAND.C - Last Edited 6. 8. 1988
|
***********************************************************************/

/*======================================================================
|Syntax of the manual pages:
|
|FUNCTION NAME(...) params ...
|
$  usage of the function and type of the parameters
?  explane the effects of the function
=  return value and the type of value if not of type int
@  globals effected directly by this routine
!  current known bugs or limitations
&  functions called by this function
~  these functions may interest you as an alternative function or
|  because they control this function somehow
^=====================================================================*/


/*
 * $Id: urand.c,v 1.4 2006/02/07 10:24:44 jpr Exp $ 
 *
 * $Log: urand.c,v $
 * Revision 1.4  2006/02/07 10:24:44  jpr
 * *** empty log message ***
 *
 * Revision 1.3  2006/02/07 10:21:42  jpr
 * Changed visibility of some variables to local scope.
 *
 * Revision 1.2  2005/05/27 12:26:22  vierinen
 * changed header install location
 *
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:57  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"

double urand(int *iy)
/*======================================================================
?  urand is a uniform random number generator based  on  theory  and
|  suggestions  given  in  d.e. knuth (1969),  vol  2.   the integer  iy
|  should be initialized to an arbitrary integer prior to the first call
|  to urand.  the calling program should  not  alter  the  value  of  iy
|  between  subsequent calls to urand.  values of urand will be returned
|  in the interval (0,1).
|  see forsythe, malcolm and moler (1977).
^=====================================================================*/
{
   static double s, halfm;
   static int  ia, ic, m, mic;
   static int m2 = 0, itwo = 2;
  
  if (m2 == 0)
  {
    
    /* if first entry, compute machine integer word length */
    
    m = 1;
    do
    {
      m2 = m;
      m = itwo * m2;
    } while(m > m2);
    halfm = m2;
    
    /* compute multiplier and increment for linear congruential method */
    
    ia = 8 * (int)(halfm * atan(1.0) / 8.00) + 5;
    ic = 2*(int)(halfm * (0.50 - sqrt(3.0) / 6.0)) + 1;
    mic = (m2 - ic) + m2;
    
    /* s is the scale factor for converting to floating point */
    
    s = 0.5 / halfm;
    
  }
  
  /* compute next random number */
  
  *iy = *iy * ia;
  
  /*
    the following statement is for computers which do not allow
    integer overflow on addition
    */
  
  if (*iy > mic) *iy = (*iy - m2) - m2;
  
  *iy = *iy + ic;
  
  /*
    the following statement is for computers where the
    word length for addition is greater than for multiplication
    */
  
  if (*iy / 2 > m2) *iy = (*iy - m2) - m2;
  
  /*
    the following statement is for computers where integer
    overflow affects the sign bit
    */
  
  if (*iy < 0) *iy = (*iy + m2) + m2;

  return *iy * s;
}
