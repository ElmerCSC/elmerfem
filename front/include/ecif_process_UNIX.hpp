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
Module:     ecif_process_UNIX.hpp
Language:   C++
Date:       20.01.99
Version:    1.00
Author(s):  Martti Verho
Revisions:  
 
Abstract:   Implementation, UNIX specific

************************************************************************/

// Helper routines
// ===============
void
display_system_msg(char* msg)
{
  cerr << msg << endl;
}


// System specific, called in main
// ===============================
void initConsole()
{
}


// =====================
// Process class methods
// =====================

bool
Process::exists()
{
  return true;
}


bool
Process::resume()
{
  return false;
}


Hfile
Process::setLogfile()
{
  if (logfileName == NULL ||
      logfileName[0] == '\0'
     )
    return 0;

  // STDERR
  return 2;
}


bool
Process::start()
{
  strstream strm;
  strm << command << ' ' << arguments << ends;

  system(strm.str());

  return true;
}


bool
Process::stop()
{
  return false;
}


bool
Process::suspend()
{
  return false;
}


void
Process::setPriorityLevel(priorityLevel level)
{ 
  priority = level;
}

