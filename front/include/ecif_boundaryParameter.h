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
Module:     ecif_boundaryParameter.h
Language:   C++
Date:       24.01.01
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   A Base class for boundary parameters.

************************************************************************/

#ifndef _ECIF_BOUNDARY_PARAMETER_
#define _ECIF_BOUNDARY_PARAMETER_

#include "ecif_parameter.h"


// ****** BoundaryParameter parameter class ******
class BoundaryParameter : public Parameter{
public:
  BoundaryParameter();
  BoundaryParameter(int pid);
  BoundaryParameter(int pid, int parent_id, char* data_string, char* param_name);
  int getLastId() {return last_id;}
  void setLastId(int lid) {last_id = lid;}
  const char* getGuiName() { return "BoundaryParameter"; }
  const char* getArrayName() { return "BoundaryParameter"; }
  const char* getEmfName() { return "Boundary Parameter"; }
  const char* getSifName() { return "Boundary Parameter"; }
  ecif_parameterType getParameterType() { return ECIF_BOUNDARY_PARAMETER; }
  static void initClass(Model* model);
  void setName(char* param_name);
protected:
  static int last_id;
  static Model* model;
};


#endif
