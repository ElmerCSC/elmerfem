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
Module:     ecif_process_WIN32.hpp
Language:   C++
Date:       20.01.99
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation, Win32 specific

************************************************************************/


// Helper routines
// ===============


void display_system_msg(char* header)
{
  LPVOID lpMsgBuf;
  FormatMessage( FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM, NULL,
                 GetLastError(),
                 MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
                 (LPTSTR) &lpMsgBuf,    0,    NULL );// Display the string.

  MessageBox( NULL, (const char*)lpMsgBuf, header, MB_OK|MB_ICONINFORMATION );
  // Free the buffer.
  LocalFree( lpMsgBuf );
}


void display_msg(char* header, char* msg)
{
  MessageBox( NULL, (const char*)msg, header, MB_OK|MB_ICONINFORMATION );
}


BOOL CtrlHandler(DWORD fdwCtrlType)
{
  switch (fdwCtrlType) {

  case CTRL_C_EVENT:
    Beep(1000, 500);
    display_system_msg("Ctrl-C: Terminating process!");
    return FALSE;

  case CTRL_BREAK_EVENT:
    Beep(1000, 1000);
    display_system_msg("Ctrl-Break: Terminating process!");
    return FALSE;

  case CTRL_CLOSE_EVENT:
    Beep(1000, 2000);
    display_system_msg("Close: Terminating process!");
    return FALSE;

  case CTRL_LOGOFF_EVENT:
  case CTRL_SHUTDOWN_EVENT:
  default:
    return FALSE;
  }
}



// Console settings, system specific
// Called in main.cpp
//
void initConsole()
{
  BOOL rc = TRUE;

  //rc = SetConsoleCtrlHandler((PHANDLER_ROUTINE) CtrlHandler, TRUE);

  if (rc == FALSE)
    display_system_msg("SetConsoleCtrlHandler");
}


// =====================
// Process class methods
// =====================


bool
Process::exists()
{
  DWORD exit_code;

  BOOL rc =  GetExitCodeProcess(processHandle, &exit_code);

  if ( !rc ) {
    display_system_msg("Checking process existence (GetExitCodeProcess)");
    return false;
  }

  if ( exit_code != STILL_ACTIVE ) {
    return false;

  } else {
    return true;
  }
}


bool
Process::resume()
{
  BOOL rc;

  rc = ResumeThread(threadHandle);

  if (rc == 0xFFFFFFFF) {
    display_system_msg("Resume thread");
    return false;
   }

  return true;
}


Hfile
Process::setLogfile()
{
  if (logfileName == NULL ||
      logfileName[0] == '\0'
     )
    return 0;

  // Create file (for "redirection")
  // ===============================
  // Create inheritable file object handle
  SECURITY_ATTRIBUTES saAttr;
  saAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
  saAttr.bInheritHandle = TRUE;
  saAttr.lpSecurityDescriptor = NULL;

  HANDLE outfile = CreateFile( logfileName,
                               GENERIC_READ | GENERIC_WRITE,
                               FILE_SHARE_READ | FILE_SHARE_WRITE,
                               &saAttr,
                               CREATE_ALWAYS,
                               FILE_ATTRIBUTE_NORMAL,
                               //FILE_FLAG_NO_BUFFERING,
                               //FILE_FLAG_WRITE_THROUGH,
                               NULL);

  if (outfile == INVALID_HANDLE_VALUE) {
    display_system_msg("Creating process logfile");
    return 0;
  }

  return outfile;
}


// ======================
// File based redirection
// ======================
#if 1

bool
Process::start()
{
  BOOL rc;

  bool success_flag = false;

  DWORD winmode; // Window mode for Startupinfo
  DWORD stdmode; // Stdoutput flag use mode

  HANDLE hStdoutRd, hChildStdoutWr;

  hChildStdoutWr = setLogfile();

  // No output,
  if ( !showConsole && hChildStdoutWr == 0) {
    winmode = SW_HIDE;
    stdmode = 0;

  // Output to console
  // Console will be created for the parent, so we can get
  // a handle fore it!
  } else if (showConsole) {
    winmode = SW_SHOW;
    stdmode = 0;

  // Output to file
  // Console will be created for the child
  // and then redirected to a file!
  } else {
    winmode = SW_HIDE;
    stdmode = STARTF_USESTDHANDLES;
  }

  outputHandle = hChildStdoutWr;

  // Set process parameters for the child process
  // ====================================
  STARTUPINFO si;

  GetStartupInfo(&si);
  si.lpTitle = NULL;
  si.dwFlags = STARTF_USESHOWWINDOW | stdmode;
  si.hStdOutput = hChildStdoutWr;
  si.hStdError = hChildStdoutWr;
  si.wShowWindow = winmode;

  DWORD pclass;

  switch (priority) {

  case ECIF_LOW_PRIORITY:
    pclass = IDLE_PRIORITY_CLASS;
    break;

  case ECIF_LOWER_THAN_NORMAL_PRIORITY:
    pclass = IDLE_PRIORITY_CLASS;
    //pclass = BELOW_NORMAL_PRIORITY_CLASS;  // NOTE: This is W2000 class

  case ECIF_NORMAL_PRIORITY:
    pclass = NORMAL_PRIORITY_CLASS;
    break;

  case ECIF_HIGHER_THAN_NORMAL_PRIORITY:
    pclass = NORMAL_PRIORITY_CLASS;
    //pclass = ABOVE_NORMAL_PRIORITY_CLASS;  // NOTE: This is W2000 class
    break;

  case ECIF_HIGH_PRIORITY:
    pclass = HIGH_PRIORITY_CLASS;
    break;

  default:
    pclass = NORMAL_PRIORITY_CLASS;
  }

  PROCESS_INFORMATION pi;
  DWORD flags = CREATE_NEW_CONSOLE | CREATE_NEW_PROCESS_GROUP | pclass;

  // Create process
  // ==============
  strstream strm;
  strm << command << " " << arguments << ends;

  // NOTE: When using this calling format, PATH is used to find
  // the executable, needed for F90 compiler call!
  rc = CreateProcess( NULL, strm.str(),
                      NULL, NULL,
                      TRUE,       // inherit handles
                      flags,      // creation flags
                      NULL, NULL, // use callers environment and cur-dir
                      &si, &pi);

  if (rc == FALSE) {
    strstream strm;

    if ( name != NULL ) {
      strm << "Starting process: " << command << " " << name << ends;
    } else {
      strm << "Starting process: " << command  << ends;
    }

    display_system_msg(strm.str());
    return false;

  } else {
    processId = pi.dwProcessId;
    processHandle = pi.hProcess;

    threadHandle = pi.hThread;
    threadId = pi.dwThreadId;

    success_flag = true;

  } // if process started

  return success_flag;
}
#endif



// ======================
// Pipe based redirection
// ======================
#if 0

bool
Process::start()
{
  char* consoleTitle = "MyTestConsole";

  BOOL rc;
  HANDLE hChildStdoutRd, hChildStdoutWr, hChildStdoutRdDup, hSaveStdout;

  bool success_flag = false;

  DWORD winmode; // Window mode for Startupinfo
  DWORD stdmode; // Stdoutput flag use mode

  SECURITY_ATTRIBUTES saAttr;
  // Set the bInheritHandle flag so pipe handles are inherited.
  saAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
  saAttr.bInheritHandle = TRUE;
  saAttr.lpSecurityDescriptor = NULL;

  // The steps for redirecting child process's STDOUT:
  //  1. Save current STDOUT, to be restored later.
  //  2. Create anonymous pipe to be STDOUT for child process.
  //  3. Set STDOUT of the parent process to be write handle to
  //     the pipe, so it is inherited by the child process.
  //  4. Create a noninheritable duplicate of the read handle and
  //     close the inheritable read handle.

  // Save the handle to the current STDOUT.
  hSaveStdout = GetStdHandle(STD_OUTPUT_HANDLE);

  // Create a pipe for the child process's STDOUT.
  if (! CreatePipe(&hChildStdoutRd, &hChildStdoutWr, &saAttr, 0))
    display_system_msg("Stdout pipe creation");

  // Set a write handle to the pipe to be STDOUT.

  if (! SetStdHandle(STD_OUTPUT_HANDLE, hChildStdoutWr))
    display_system_msg("Redirecting STDOUT");

  // Create noninheritable read handle and close the inheritable read
  // handle.
  rc = DuplicateHandle(GetCurrentProcess(), hChildStdoutRd,
                       GetCurrentProcess(), &hChildStdoutRdDup , 0,
                       FALSE,
                       DUPLICATE_SAME_ACCESS);

  if( !rc )
    display_system_msg("DuplicateHandle");

  CloseHandle(hChildStdoutRd);

  bool silent = false;

  // No output,
  if ( !showConsole && silent) {
    winmode = SW_HIDE;
    stdmode = 0;

  // Output to console
  // Console will be created for the parent, so we can get
  // a handle fore it!
  } else if (showConsole) {
    winmode = SW_SHOW;
    stdmode = 0;

  // Output to file
  // Console will be created for the child
  // and then redirected to a file!
  } else {
    winmode = SW_HIDE;
    stdmode = STARTF_USESTDHANDLES;
  }

  // Set process parameters for the child process
  // ====================================
  STARTUPINFO si;
  //si.cb = sizeof(si);
  GetStartupInfo(&si);
  si.lpTitle = NULL;
  //si.lpTitle = consoleTitle;
  si.dwFlags = STARTF_USESHOWWINDOW | stdmode;
  si.hStdOutput = hChildStdoutWr;
  si.hStdError = hChildStdoutWr;
  si.wShowWindow = winmode;

  outputHandle = hChildStdoutRdDup;

  PROCESS_INFORMATION pi;
  DWORD flags = CREATE_NEW_CONSOLE | CREATE_NEW_PROCESS_GROUP | IDLE_PRIORITY_CLASS;
  //DWORD flags = CREATE_NEW_CONSOLE | IDLE_PRIORITY_CLASS;
  //DWORD flags = DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP | IDLE_PRIORITY_CLASS;

  // Create process
  // ==============

  rc = CreateProcess( command, arguments,
                      NULL, NULL,
                      TRUE,   // inherit handles
                      flags,  // creation flags
                      NULL, NULL,
                      &si, &pi);
  if (rc == FALSE) {
    display_system_msg("Create process");
  }
  else {
    processId = pi.dwProcessId;
    processHandle = pi.hProcess;

    threadHandle = pi.hThread;
    threadId = pi.dwThreadId;

#if 0
    HWND wHandle = FindWindowEx( NULL, NULL, "ConsoleWindowClass",  consoleTitle);
    if (wHandle == FALSE) {
      display_system_msg("Find Window");
    }

    rc = PostMessage( wHandle, WM_CLOSE, 0 ,0 );
    if (rc == FALSE) {
      display_system_msg("Post Msg");
    }

#endif

    //rc = SuspendThread(threadHandle);
    //rc = TerminateThread(threadHandle, 0);

    //if (rc == 0) {
    //  display_system_msg("Terminate thread");
    //}

    //Sleep(500);

    success_flag = true;

  } // if process started


  return success_flag;
}
#endif


bool
Process::stop()
{
  int rc;
  rc = TerminateProcess(processHandle, 0);

  if (rc == 0) {
    return false;
   }

  CloseHandle(threadHandle);
  CloseHandle(processHandle);


  CloseHandle(outputHandle);

  return true;
}


bool
Process::suspend()
{
  BOOL rc;

  rc = SuspendThread(threadHandle);

  if (rc == 0xFFFFFFFF) {
    display_system_msg("Suspend thread");
    return false;
   }

  return true;
}


void
Process::setPriorityLevel(priorityLevel level)
{
  priority = level;

  switch (level) {
  case ECIF_LOW_PRIORITY:
    SetPriorityClass(processHandle, IDLE_PRIORITY_CLASS);
    break;
  case ECIF_NORMAL_PRIORITY:
    SetPriorityClass(processHandle, NORMAL_PRIORITY_CLASS);
    break;
  case ECIF_HIGH_PRIORITY:
    SetPriorityClass(processHandle, HIGH_PRIORITY_CLASS);
    break;
  }
}


