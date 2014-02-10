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
Module:     ecif_datafile.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   A Base class for input and resul file parameter definition

************************************************************************/

#ifndef _ECIF_DATAFILE_
#define _ECIF_DATAFILE_

#include "ecif_parameter.h"


// ****** Datafile parameter class ******
class Datafile : public Parameter{
public:
  Datafile();
  Datafile(int pid);
  Datafile(int pid, char* data_string, char* param_name);
  int getLastId() {return last_id;}
  void setLastId(int lid) {last_id = lid;}
  const char* getGuiName() { return "Datafile"; }
  const char* getArrayName() { return "Datafile"; }
  const char* getEmfName() { return "Datafile"; }
  const char* getSifName() { return "Datafile"; }
  ecif_parameterType getParameterType() { return ECIF_DATAFILE; }
  static void initClass(Model* model);
  virtual ostream& output_sif(ostream& out, short indent_size, short indent_level, SifOutputControl& soc);
  void setName(char* param_name);
protected:
  static int last_id;
  static Model* model;
};


#endif
