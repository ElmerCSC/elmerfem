/*****************************************************************************
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************
!
! ******************************************************************************
! *
! * Provide system time / memory usage.
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Jun 1997
! *
! *****************************************************************************/

#include "../config.h"

#if defined(MINGW32) || defined(WIN32) 

#include <sys/types.h>
#include <time.h>
#include <windows.h>
#include <psapi.h>

double STDCALLBULL FC_FUNC(realtime,REALTIME) ( )
{
  return clock() / (double)CLOCKS_PER_SEC;
}

double STDCALLBULL FC_FUNC(cputime,CPUTIME) ( )
{
  /* User CPU time consumed by this process, summed over all its threads --
     the analogue of getrusage(RUSAGE_SELF).ru_utime used by the POSIX branch
     below.

     clock() must not be used here.  On Windows it measures elapsed wall time
     rather than consumed CPU time, so it returns the same value as realtime()
     above, and every (CPU,REAL) pair Elmer prints becomes two copies of the
     same number.  That makes a threaded run indistinguishable from a serial
     one in the solver's own output. */
  FILETIME created, exited, kernel, user;
  ULARGE_INTEGER t;

  if ( !GetProcessTimes( GetCurrentProcess(), &created, &exited, &kernel, &user ) )
    return 0.0;

  t.LowPart  = user.dwLowDateTime;
  t.HighPart = user.dwHighDateTime;

  /* FILETIME counts 100-nanosecond intervals. */
  return (double)t.QuadPart * 1.0e-7;
}

double STDCALLBULL  FC_FUNC(cpumemory,CPUMEMORY) ( )
{
  /* Peak working set, the analogue of getrusage(RUSAGE_SELF).ru_maxrss.
     ru_maxrss is reported in kilobytes, so scale to match: the callers of
     this function compare the two platforms' numbers. */
  PROCESS_MEMORY_COUNTERS pmc;

  if ( !GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc) ) )
    return 0.0;

  return (double)pmc.PeakWorkingSetSize / 1024.0;
}

#else

#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>

#include <sys/time.h>
#include <sys/resource.h>

static struct rusage usage;

static struct timeval tp;
static struct timezone tzp;

#ifdef USE_ISO_C_BINDINGS
double cputime ()
{
  getrusage( RUSAGE_SELF, &usage );
  return (double) usage.ru_utime.tv_sec + usage.ru_utime.tv_usec*1.0e-6;
}

double realtime()
{
  gettimeofday( &tp,&tzp );
  return (double) tp.tv_sec + tp.tv_usec*1.0e-6;
}

double cpumemory()
{
  getrusage( RUSAGE_SELF, &usage );
  return (double) 1.0 * usage.ru_maxrss;
}
#else
double FC_FUNC(cputime,CPUTIME) ()
{
  getrusage( RUSAGE_SELF, &usage );
  return (double) usage.ru_utime.tv_sec + usage.ru_utime.tv_usec*1.0e-6;
}

double FC_FUNC(realtime,REALTIME) () 
{
  gettimeofday( &tp,&tzp );
  return (double) tp.tv_sec + tp.tv_usec*1.0e-6;
}

double FC_FUNC(cpumemory,CPUMEMORY) ()
{ 
  getrusage( RUSAGE_SELF, &usage );
  return (double) 1.0 * usage.ru_maxrss;
}
#endif /* USE_ISO_C_BINDINGS*/

#endif // WIN32
