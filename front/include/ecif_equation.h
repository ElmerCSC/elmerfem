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
Module:     ecif_equation.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   A Base class for equationbody-force.

************************************************************************/

#ifndef _ECIF_EQUATION_
#define _ECIF_EQUATION_

#include "ecif_parameter.h"


// ****** Equation class ******
class Equation : public Parameter {
public:
  Equation();
  Equation(int pid);
  Equation(int pid, char* values, char* param_name);
  int getLastId() {return last_id;}
  void setLastId(int lid) {last_id = lid;}
  const char* getGuiName() { return "Equation"; }
  const char* getArrayName() { return "Equation"; }
  const char* getEmfName() { return "Equation"; }
  const char* getSifName() { return SIF_EQUATION; }
  ecif_parameterType getParameterType() { return ECIF_EQUATION; }
  static void initClass(Model* model);
  virtual ostream& output_sif(ostream& out, short indent_size, short indent_level, SifOutputControl& soc);
  virtual ostream& outputSolverTargetFields_sif(ostream& out, short indent_size, short indent_level, const char* source_eq_name, NameSet& targetFieldNames);
  void setName(char* param_name);
protected:
  static int last_id;
  static Model* model;
  ostream& output_equationWithVariables_sif(ostream& out, short indent_size, short indent_level,
                                            ParameterField* equation_pf,
                                            ParameterField* equation_vars_pf);
};


#endif
