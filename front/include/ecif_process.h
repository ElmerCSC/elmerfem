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
Module:     ecif_process.h
Language:   C++
Date:       20.01.99
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Process class for handling external processes.

************************************************************************/

#ifndef _PROCESS_
#define _PROCESS_


#include "ecif_def.h"


class Process
{
public:
  Process();
  Process(char* process_cmd, char* process_args);
  Process(char* process_cmd, char* process_args,
          int process_nbr, char* process_name,
          enum priorityLevel process_priority,
          bool show_console, char* logfile_name);
  ~Process();
  bool exists();
  const char* getName() { return name;}
  Hfile getOutputHandle() { return outputHandle;}
  priorityLevel getPriorityLevel() { return priority; }
  Hprocess getProcessHandle() { return processHandle;}
  ProcessId getProcessId() { return processId;}
  Hprocess getThreadHandle() { return threadHandle;}
  ProcessId getThreadId() { return threadId;}
  int ID() { return id;}
  bool resume();
  void setLogfile(char* logfile);
  void setShowConsole(bool value) { showConsole = value; }
  void setPriorityLevel(priorityLevel level);
  bool start();
  bool stop();
  bool suspend();

private:
  char* arguments;
  char* command;
  int id;
  Hfile logfileHandle;
  char* logfileName;
  char* name;
  Hfile outputHandle;
  enum priorityLevel priority;
  Hprocess processHandle;
  ProcessId processId;
  bool showConsole;
  bool started;
  bool stopped;
  Hprocess threadHandle;
  ProcessId threadId;

  void init();
  Hfile setLogfile();
};


#endif
