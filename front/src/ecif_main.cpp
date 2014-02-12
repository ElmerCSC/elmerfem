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
Module:     ecif_main.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include <stdio.h>
#include "ecif_userinterface.h"
#include "ecif_userinterface_TCL.h"
#include "ecif_control.h"
#include "ecif_func.h"
#include "ecif_renderer.h"

// Elmer MATC init
extern "C" void mtc_init(FILE* infile, FILE* outfile, FILE* errfile);

// ***** Global variables *****
//Output-file for debug-prints.
ostream* debugFile = new ofstream();

// Set static class variables
enum rendererStatus Renderer::status = HAS_NO_WINDOW;
char* tclMainScript = "ecif_tcl_mainScript.tcl";


// Parse command line (for batch-mode)
//
void
parse_cmd_line(int argc, char** argv, int& cmdc, char**& cmdv)
{
  int i,j;

  cmdc = 0;
  cmdv = NULL;

  bool was_sep = false;

  // Loop arg vectors and change =--= markers to \n
  // for splitting
  //
  for (i = 0; i < argc; i++) {

    int len = strlen(argv[i]);

    int pos = 0;
    for (j = 0; j < len;  j++) {
      char c = argv[i][j];
      if ( c == '-' ) {
        if ( was_sep ) {
          cmdc++;
          argv[i][pos-1] = '\n';
          was_sep = false;
          continue;
        } else {
          was_sep = true;
        }
      } else {
        was_sep = false;
      }

      argv[i][pos++] = c;
    }

    argv[i][pos] = '\0';
  }

  cmdv = new char*[cmdc];

  // Loop arg vectors once more and
  // parse commands
  //
  int idx = 0;
  for (i = 0; i < argc; i++) {

    strstream strm;
    strm << argv[i];

    while ( !strm.eof() && idx < cmdc ) {
      char buffer[1024];

      strm.getline(buffer, 1023);

      if ( strlen(buffer) == 0 ) continue;

      cmdv[idx] = NULL;
      update_dyna_string(cmdv[idx], LibFront::trim(buffer));
      idx++;
    }
  }
}


// Checks is Front is started in batch-mode
//
bool
is_batch_mode(int cmdc, char** cmdv) {

  for ( int i = 0; i < cmdc; i++) {
    if ( 0 == strncmp(LibFront::trimLeft(cmdv[i]), "batch=1", 7) ) {
      return true;
    }
  }

  return false;
}


// Console version (Unix and Win32)
// --------------------------------
int main(int argc, char** argv)
{
  int arg_c = 0;
  char** arg_v = NULL;

  // Program name is dropped
  arg_c = argc - 1;
  arg_v = ++argv;
  Hinst appInstance = (Hinst)argv[0];

  // Note this cal is system specific
  // Ref. ecif_processWIN32/UNIX.hpp
  void initConsole();
  initConsole();

  // Init MATC
  mtc_init(NULL, stdout, stderr);
  LibFront::initMatcFrmt();

  bool is_batch;

  // Copy original command line arguments
  //
  // NOTE: The original argument line is delivered to the
  // gui-mode start()!!!
  //
  int arg_cc = arg_c;
  char** arg_vv = new char*[arg_cc];
  for (int i = 0 ; i < arg_cc; i++) {
    arg_vv[i] = NULL;
    update_dyna_string(arg_vv[i], arg_v[i]);
  }

  // Parse command line into separate arguments and
  // remove leading -- signs from arguments
  //
  int cmdc = 0;
  char** cmdv = NULL;
  parse_cmd_line(arg_cc, arg_vv, cmdc, cmdv);

  // Check if bacth mode flag "batch=1" is given
  //
  if ( is_batch_mode(cmdc, cmdv) ) {
    is_batch = true;
  } else {
    is_batch = false;
  }

  // Select proper ui (batch/normal)
  //
  UserInterface* ui;
  if ( is_batch ) {
    ui = new UserInterface(appInstance, NULL);
  } else {
    ui = new UserInterface_TCL(appInstance, tclMainScript);
  }

  Control* cc = new Control(appInstance, ui);

  // Start processing
  //
  if ( is_batch ) {
    ui->start(cmdc, cmdv);
  } else {
    ui->start(arg_c, arg_v);
  }

  return 0;
}


// -*End of main
