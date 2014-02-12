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
Module:     ecif_timestep.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Implementation

************************************************************************/

#include "ecif_parameterField.h"
#include "ecif_timestep.h"

//Initialize static class variables.
int Timestep::last_id = 0;
Model* Timestep::model = NULL;


// Constructors
Timestep::Timestep()
{
}


Timestep::Timestep(int pid) : Parameter(pid)
{
  if ( pid > last_id ) last_id = pid;
}


Timestep::Timestep(int pid, char* data_string, char* param_name)
{
  if ( pid > last_id ) last_id = pid;

  setData(pid, data_string, param_name);
}


// Timestep specific output-method for Solver input file
ostream&
Timestep::output_sif(ostream& out, short indent_size, short indent_level, SifOutputControl& soc)
{

  char QM = '\"';


  // Fields NOT proper for steady case
  int nof_non_steady = 6;
  const char* non_steady[] = { 
    "TIMESTEPPING_METHOD",
    "NEWMARK_BETA",
    "BDF_ORDER",
    "TIMESTEP_SIZES",
    "TIMESTEP_INTERVALS",
    "OUTPUT_INTERVALS"
  };

  // Fields NOT proper for transient case
  int nof_non_transient = 1;
  const char* non_transient[] = { 
    "STEADY_STATE_OUTPUT_INTERVAL",
    //"STEADY_STATE_MAX_ITERATIONS"
  };
 

  // Parameter type and id
  if (soc.outputType) {

    LibFront::output_string(out, indent_size, indent_level++, getSifName(), false);

    if (soc.outputId)
      out << ' ' << ID();

    out << endl;
  }

  // Check simulation type
  bool is_steady = true;
  ParameterField* pf = getFieldBySifName("Simulation Type");

  char** data = pf->getDataStrings();
  if ( LibFront::in("transient", data[0]) ) {
    is_steady = false;
  }

  // Output proper fields
  for (short i = 0; i < nofFields; i++) {

    bool is_proper = true;

    ParameterField* pf = fields[i];

    // Check that field exists and it should be output
    if ( !pf->isActiveInstance() ||
         pf->getNofDataStrings() == 0 ||
         ( !soc.outputAll && !pf->isSifOutputField() )
       )
      continue;

    // Check that field is proper for the steady/transient case
    char* fn1 = (char*)pf->getGuiName();

    // Steady state case
    if (is_steady) {

      // Skip non-steady fields
      for (int k = 0; k < nof_non_steady; k++) {
        if ( LibFront::in(fn1, non_steady[k]) ) {
          is_proper = false;
          break;
        }
      }

    // Transient case
    } else {
      // Skip non-transient fields
      for (int k = 0; k < nof_non_transient; k++) {

        if ( LibFront::in(fn1, non_transient[k]) ) {
          is_proper = false;
          break;
        }
      }
    }

    // If field not proper for the case
    if (!is_proper)
      continue;
    
    pf->output_sif(out, indent_size, indent_level, soc.sectionName);

  }

  return out;
}


void
Timestep::initClass(Model* mdl)
{
  Timestep::model = mdl;
  Timestep::last_id = 0;
}


void
Timestep::setName(char* param_name)
{
  Parameter::setName(param_name, "Timestep");
}



