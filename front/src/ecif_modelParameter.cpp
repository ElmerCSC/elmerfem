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
Module:     ecif_modelParameter.cpp
Language:   C++
Date:       13.02.01
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_modelParameter.h"
#include "ecif_model.h"
#include "ecif_parameterField.h"

//Initialize static class variables.
int ModelParameter::last_id = 0;
Model* ModelParameter::model = NULL;


// Constructors
ModelParameter::ModelParameter()
{
}


ModelParameter::ModelParameter(int pid) : Parameter(pid)
{
}


ModelParameter::ModelParameter(int pid, char* data_string, char* param_name)
{
  setData(pid, data_string, param_name);
}


void
ModelParameter::initClass(Model* mdl)
{
  ModelParameter::model = mdl;
  ModelParameter::last_id = 0;
}


void
ModelParameter::setName(char* param_name)
{
  Parameter::setName(param_name, "ModelParameter");
}

