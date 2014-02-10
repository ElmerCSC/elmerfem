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
Module:     ecif_simulationParameter.cpp
Language:   C++
Date:       13.02.01
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_simulationParameter.h"
#include "ecif_model.h"
#include "ecif_parameterField.h"

//Initialize static class variables.
int SimulationParameter::last_id = 0;
Model* SimulationParameter::model = NULL;


// Constructors
SimulationParameter::SimulationParameter()
{
}


SimulationParameter::SimulationParameter(int pid) : Parameter(pid)
{
}


SimulationParameter::SimulationParameter(int pid, char* data_string, char* param_name)
{
  setData(pid, data_string, param_name);
}


void
SimulationParameter::initClass(Model* mdl)
{
  SimulationParameter::model = mdl;
  SimulationParameter::last_id = 0;
}


void
SimulationParameter::setName(char* param_name)
{
  Parameter::setName(param_name, "Simulation");
}

