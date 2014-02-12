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
Module:     ecif_simulationParameter.h
Language:   C++
Date:       13.02.01
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   A Base class for (user defined) simulation paremeters.

************************************************************************/

#ifndef _ECIF_SIMULATION_PARAMETER_
#define _ECIF_SIMUALTION_PARAMETER_

#include "ecif_parameter.h"


// ****** SimulationParameter parameter class ******
class SimulationParameter : public Parameter{
public:
  SimulationParameter();
  SimulationParameter(int pid);
  SimulationParameter(int pid, char* data_string, char* param_name);
  int getLastId() {return last_id;}
  void setLastId(int lid) {last_id = lid;}
  const char* getGuiName() { return "Simulation"; }
  const char* getArrayName() { return "SimulationParameter"; }
  const char* getEmfName() { return "Simulation Parameter"; }
  const char* getSifName() { return ""; }
  ecif_parameterType getParameterType() { return ECIF_SIMULATION_PARAMETER; }
  static void initClass(Model* model);
  void setName(char* param_name);
protected:
  static int last_id;
  static Model* model;
};


#endif
