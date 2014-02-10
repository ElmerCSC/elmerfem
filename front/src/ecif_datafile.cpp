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
Module:     ecif_datafile.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Implementation

************************************************************************/

#include "ecif_datafile.h"
#include "ecif_model.h"
#include "ecif_parameterField.h"

//Initialize static class variables.
int Datafile::last_id = 0;
Model* Datafile::model = NULL;


// Constructors
Datafile::Datafile()
{
}


Datafile::Datafile(int pid) : Parameter(pid)
{
}


Datafile::Datafile(int pid, char* data_string, char* param_name)
{
  setData(pid, data_string, param_name);
}


void
Datafile::initClass(Model* mdl)
{
  Datafile::model = mdl;
  Datafile::last_id = 0;
}


// Method prints all datafile parameter's fields into output stream.
// Output of parameter name and id number is controlled by
// output_parameter_type flag
// Check that GebhardtFactors and Viewfactors filenames
// are really needed, ie. mode has DiffuseGray boundary conditions
ostream&
Datafile::output_sif(ostream& out, short indent_size, short indent_level, SifOutputControl& soc)
{
  if (soc.outputType) {

    LibFront::output_string(out, indent_size, indent_level, getSifName());

    if (soc.outputId)
      out << " " << id;

    out << endl;
  } else {
    indent_level++;
  }

  bool has_diff_gray = model->modelHasDiffuseGrayRadiation();

  for (short i = 0; i < nofFields; i++) {

    ParameterField* pf = fields[i];

    if ( LibFront::in(pf->getSifName(), "Gebhardt Factors") ||
         LibFront::in(pf->getSifName(), "View Factors") 
       ) {

      if (!has_diff_gray) {
        continue;
      }
    }

    // Check that field exists and it should be output
    if ( !pf->isActiveInstance() || 
         pf->getNofDataStrings() == 0 ||
         ( !soc.outputAll && !pf->isSifOutputField() )
         ) {
      continue;
    }

    pf->output_sif(out, indent_size, indent_level, soc.sectionName);

  }

  return out;
}


void
Datafile::setName(char* param_name)
{
  Parameter::setName(param_name, "Datafile");
}



