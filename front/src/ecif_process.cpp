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
Module:     ecif_process.cpp
Language:   C++
Date:       20.01.99
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation.

************************************************************************/

#if defined(WIN32)
  // Nothing special
#else
  #include <unistd.h>
#endif

#include "ecif_process.h"
#include "ecif_func.h"

Process::Process()
{
  init();
}


Process::Process(char* process_cmd, char* process_args)
{
  init();
  update_dyna_string(command, process_cmd);
  update_dyna_string(arguments, process_args);
}


Process::Process(char* process_cmd, char* process_args,
                 int process_nbr, char* process_name,
                 enum priorityLevel process_priority,
                 bool show_console, char* logfile_name)
{
  init();
  update_dyna_string(command, process_cmd);
  update_dyna_string(arguments, process_args);
  id = process_nbr;
  update_dyna_string(name, process_name);
  priority = process_priority;
  showConsole = show_console;
  update_dyna_string(logfileName, logfile_name);
}


Process::~Process()
{
  delete[] command;
  delete[] arguments;
  delete[] name;
  delete[] logfileName;
}

void
Process::init()
{
  command = NULL;
  arguments = NULL;
  name = NULL;
  logfileName = NULL;

  priority = ECIF_NO_PRIORITY;
  id = 0;
  processId = 0;
  processHandle = 0;
  showConsole = false;
  logfileHandle = 0;

  threadHandle = NULL;
  threadId = NULL;

  started = 0;
  stopped = 0;
}


void
Process::setLogfile(char* logfile)
{
  if (logfile != NULL ) {
    update_dyna_string(logfileName, logfile);
  }
}


// Platform specific parts

#if defined(WIN32)
  #include "ecif_process_WIN32.hpp"

#else
  #include "ecif_process_UNIX.hpp"
#endif

