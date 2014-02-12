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
Module:     ecif_equation.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_control.h"
#include "ecif_equation.h"
#include "ecif_model.h"
#include "ecif_parameterField.h"
#include "ecif_userinterface.h"

//Initialize static class variables.
int Equation::last_id = 0;
Model* Equation::model = NULL;


// Constructors
Equation::Equation()
{
}


Equation::Equation(int pid) : Parameter(pid)
{
}


Equation::Equation(int pid, char* values, char* param_name)
{
  setData(pid, values, param_name);
}


void
Equation::initClass(Model* mdl)
{
  Equation::model = mdl;
  Equation::last_id = 0;
}


// Equation specific output-method for Solver input file
ostream&
Equation::output_sif(ostream& out, short indent_size, short indent_level, SifOutputControl& soc)
{
  UserInterface* gui = (UserInterface*)model->getGui();


  char QM = '\"';

  // Parameter type and id  (Equation 1)
  if (soc.outputType) {

    LibFront::output_string(out, indent_size, indent_level++, getSifName(), false);

    if (soc.outputId)
      out << ' ' << ID();

    out << endl;

  } else {
    indent_level++;
  }

  // Output parameter name
  if (soc.outputName) {
    output_sif_name(out, indent_size, indent_level, soc);
    out << endl;
  }

  // Fields
  char* equation_vars_name = NULL;

  for (short i = 0; i < nofFields; i++) {

    ParameterField* pf = fields[i];

    // Check that field exists and it should be output
    if ( !pf->isActiveInstance() ||
         pf->getNofDataStrings() == 0 ||
         ( !soc.outputAll && !pf->isSifOutputField() )
         ) {
      continue;
    }

    // If field should be output to some other sif section (ie. in Solver section)
    if ( gui->getIsSolverTargetField(NULL, pf->getGuiName()) ) {
      continue;
    }

    ParameterField* equation_vars_pf = NULL;

    // Check if variable name should (and could) be added to
    // the equation name
    if ( gui->getUseVariableNameInEquationName(pf->getGuiName()) ) {

      // Check if some subsystem variable values are defined (like AD_VARS=Oxygen;Nitrogen)
      if (gui->getEquationVarsVariable(pf->getGuiName(), equation_vars_name)) {

        if (equation_vars_name != NULL ) {
          equation_vars_pf = getFieldByGuiName(equation_vars_name);
          delete[] equation_vars_name;
          equation_vars_name = NULL;
        }
      }

      // "Indexed" equation, equation_vars defined
      if (equation_vars_pf != NULL ) {
        output_equationWithVariables_sif(out, indent_size, indent_level,
                                         pf, equation_vars_pf);
      }

    // Normal equation, no equation vars
    } else {
      pf->output_sif(out, indent_size, indent_level, soc.sectionName);
    }

  }


  return out;
}


// Equation specific output-method for Solver input file
ostream&
Equation::output_equationWithVariables_sif(ostream& out, short indent_size, short indent_level,
                                           ParameterField* equation_pf,
                                           ParameterField* equation_vars_pf)
{
  const int var_count =  equation_vars_pf->getNofDataStrings();
  char** var_names = equation_vars_pf->getDataStrings();

  for ( int i = 0; i < var_count; i++) {

    LibFront::indent(out, indent_size, indent_level);

    strstream strm;

    strm << equation_pf->getSifName() << " " << var_names[i] << " = " << ends;

    out << strm.str();
    equation_pf->output_sif_type(out);
    equation_pf->output_sif_data(out, indent_size, true, indent_level);

    if ( i < var_count - 1) {
      out << endl;
    }
  }

  return out;
}


// Output those equation fields whose target is actually a solver section
//
ostream&
Equation::outputSolverTargetFields_sif(ostream& out, short indent_size, short indent_level, const char* source_eq_name, NameSet& targetFieldNames)
{
  UserInterface* gui = (UserInterface*)model->getGui();


  // Loop all equation parameter fields and check if tehy belong to the source equation and
  // if they are solver target fields
  //
  for (short i = 0; i < nofFields; i++) {

    ParameterField* pf = fields[i];

    // Check that field is active, contains data and it should be output
    if ( !pf->isActiveInstance() ||
         pf->getNofDataStrings() == 0 ||
         !pf->isSifOutputField()
       )
      continue;

    if ( !gui->getIsSolverTargetField(source_eq_name, pf->getGuiName()) ) {
      continue;
    }

    std::pair<NameSet::iterator, bool> tfPair;

    // Try to insert a new target field, check if it is already inserted
    tfPair = targetFieldNames.insert(pf->getGuiName());

    // Exists already
    if ( !tfPair.second ) continue;

    // Ok, output field
    pf->output_sif(out, indent_size, indent_level, SIF_SOLVER);

  }

  return out;

}




void
Equation::setName(char* param_name)
{
  Parameter::setName(param_name, "Equation");
}


