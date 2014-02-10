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

/***********************************************************************
Program:    ELMER Front 
Module:     ecif_timer.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Simple timer class for elapsed time.

************************************************************************/

#ifndef _TIMER_
#define _TIMER_

#include "ecif_def.h"

class Timer
{
public:
  Timer();   
  void start();
  void stop();
  double getLapTime(enum timeType time_type = PROCESS_TIME);
  double getEndTime(enum timeType time_type = PROCESS_TIME);
  static void formTimeString(double seconds, char* buffer);

private:
  double getProcessTime();
  double getWallClockTime();

  short started;
  short stopped;
  double start_process_time;
  double end_process_time;
  double start_wall_time;
  double end_wall_time;  
};


#endif
