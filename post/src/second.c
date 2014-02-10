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
 * Provide system time.
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
#ifndef WIN32

#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>

#include <sys/time.h>
#include <sys/resource.h>

  static struct tms t;
  static struct rusage usage;

  static struct timeval tp;
  static struct timezone tzp;

  double CPUTime( )
  {
    getrusage( RUSAGE_SELF, &usage );
    return (double) usage.ru_utime.tv_sec + usage.ru_utime.tv_usec*1.0e-6;
  }

  double RealTime( )
  {
    gettimeofday( &tp,&tzp );
    return tp.tv_sec + tp.tv_usec*1.0e-6;
  }

#else

  #include <sys/types.h>
  #include <time.h>

  double CPUTime() { return 0.0; }

  double RealTime()
  {
     return clock() / (double)CLOCKS_PER_SEC;
  }

#endif
