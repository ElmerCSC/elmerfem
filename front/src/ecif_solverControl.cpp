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
Module:     ecif_solver.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Implementation

************************************************************************/

#include "ecif_model.h"
#include "ecif_parameterField.h"
#include "ecif_solverControl.h"
#include "ecif_userinterface.h"

//Initialize static class variables.
int SolverControl::last_id = 0;
Model* SolverControl::model = NULL;


// Constructors
SolverControl::SolverControl()
{
}


SolverControl::SolverControl(int pid) : Parameter(pid)
{
}


SolverControl::SolverControl(int pid, char* data_string, char* param_name)
{
  setData(pid, data_string, param_name);
}


void
SolverControl::initClass(Model* mdl)
{
  SolverControl::model = mdl;
  SolverControl::last_id = 0;
}


// Outputs all fields with some exceptions (ECHO_ON, CHECK_KEYWORDS)
//
// NOTE: This is used for Simulation block output
//
ostream&
SolverControl::output_sif(ostream& out, short indent_size, short indent_level, SifOutputControl& soc)
{
  char QM = '\"';
 
  for (short i = 0; i < nofFields; i++) {

    ParameterField* pf = fields[i];

    // Check that field exists and it should be output
    if ( !pf->isActiveInstance() ||
         pf->getNofDataStrings() == 0 ||
         ( !soc.outputAll && !pf->isSifOutputField() )
       )
      continue;

    char* fld_name = (char*)pf->getSifName();
    
    // Skip fields which are not output here!!!
    if ( LibFront::in(fld_name, SIF_ECHO_ON) ||
         LibFront::in(fld_name, SIF_CHECK_KEYWORDS) ||
         LibFront::in(fld_name, SIF_RELOAD_INPUT_FILE)
       ) {
      continue;
    }
    
    pf->output_sif(out, indent_size, 1 + indent_level, soc.sectionName);

  }

  return out;
}

void
SolverControl::setName(char* param_name)
{
  Parameter::setName(param_name, "SolverControl");
}




