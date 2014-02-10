/*  
   ElmerFront - A graphical user interface of Elmer software
   Copyright (C) 1995- , CSC - IT Center for Science Ltd.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/***********************************************************************
Program:    ELMER Front 
Module:     ecif_timer.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  
 
Abstract:   Implementation. 

************************************************************************/

#if defined(WIN32)
  #include <sys/timeb.h>
#else
  #include <sys/time.h>
  static struct timeval tp;
  static struct timezone tzp;
#endif


#include "ecif_timer.h"

Timer::Timer()
{
  started = 0;
  stopped = 0;
}


double
Timer::getLapTime(enum timeType time_type)
{
  // Not started
  if (!started) {
    return -1;

  //---Process time was asked
  } else if (time_type == PROCESS_TIME) {
    return getProcessTime() - start_process_time;

  //---Wall clock time was asked
  } else {
    return getWallClockTime() - start_wall_time;
  }
}


void
Timer::start()
{
  //---Start process timer
  start_process_time = getProcessTime();

  //---Start wall clock timer
  start_wall_time = getWallClockTime();

  started = 1;
  stopped = 0;
}


void
Timer::stop()
{
  // Not running
  if (!started) {
    return;
  }

  //---Stop process timer 
  end_process_time = getProcessTime();

  //---Stop wall clock timer
  end_wall_time = getWallClockTime();

  stopped = 1;
}


double
Timer::getEndTime(enum timeType time_type)
{
  // Not stopped
  if (!stopped) {
    return -1.0;
  
  //---Process time was asked
  } else if (time_type == PROCESS_TIME) {     
    return  end_process_time - start_process_time;

  //---Wall-clock time was asked
  } else {
    return end_wall_time - start_wall_time;
  }
}


double
Timer::getProcessTime()
{
 return ((double)clock()) / CLOCKS_PER_SEC;
}


double
Timer::getWallClockTime()
{

#if defined(WIN32)
  static struct _timeb ts;
  _ftime(&ts);
  return (ts.time + 1.0e-3 * ts.millitm);
#else
#if 0
  static struct timeb ts;
  ftime(&ts);
  return (ts.time + 1.0e-3 * ts.millitm);
#else
  gettimeofday( &tp,&tzp );
  return tp.tv_sec + tp.tv_usec*1.0e-6;
#endif
#endif

}



// Method creates a time string from 'seconds' variable, 
// result in in the format HH:MM:SS
void
Timer::formTimeString(double seconds, char* buffer)
{
    int hrs =   int( seconds / 3600.0 );
    int mns = int( (seconds - hrs * 3600.0) / 60.0 );
    int scs = int( (seconds - hrs * 3600.0 - mns * 60.0) );
    
    ostrstream ostrm;
    ostrm.setf(ios::right, ios::adjustfield);
    ostrm << setw(4) << hrs << ":";
    ostrm << setfill('0');
    ostrm << setw(2) << mns << ":" << setw(2) << scs;
    int len = ostrm.pcount();
    ostrm << ends;

    for (int i=0; i< len; i++)
        buffer[i] = (ostrm.str())[i];

    buffer[len] = '\0';
}

