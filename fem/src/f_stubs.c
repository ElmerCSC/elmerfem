/*  
   Elmer, A Finite Element Software for Multiphysical Problems
  
   Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
   
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
  
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with this library (in file ../../LGPL-2.1); if not, write 
   to the Free Software Foundation, Inc., 51 Franklin Street, 
   Fifth Floor, Boston, MA  02110-1301  USA
*/


/*
  All kinds of stubs etc that cover up if something is missing from fortran.
 */
#include "../config.h"


#include <sys/types.h>
#include <stdio.h>




#if defined(WIN32) | defined(MINGW32)

#else 

#include <sys/resource.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/param.h>

#endif



#ifndef HAVE_F_ETIME

#if defined(WIN32) | defined(MINGW32)
float STDCALLBULL FC_FUNC(etime,ETIME)(tt)
float tt[2];
{
  return(.0);
}
#else
float FC_FUNC(etime,ETIME)(tt)
float tt[2];
{
   int who;
   struct rusage used;
   who = 0;
   getrusage(who,&used);
   tt[0] = used.ru_utime.tv_sec+((used.ru_utime.tv_usec)/1000000.);
   tt[1] = used.ru_stime.tv_sec+((used.ru_stime.tv_usec)/1000000.);
   return(tt[0]+tt[1]);
}
#endif // win32
#endif // etime_defined

#ifndef HAVE_F_FLUSH
void STDCALLBULL FC_FUNC(flush,FLUSH) (int n)
{
  /*  might as well flush a toilet...? */
}
#endif

void rename_c(const char *old, const char *new)
{
  int rc = rename(old,new);
  (void) rc;
}
