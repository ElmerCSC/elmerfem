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
Module:     ecif_calculator.cpp
Language:   C++
Date:       24.01.01
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Implementation

************************************************************************/

#include "ecif_model.h"
#include "ecif_parameterField.h"
#include "ecif_calculator.h"
#include "ecif_userinterface.h"

//Initialize static class variables.
int Calculator::last_id = 0;
Model* Calculator::model = NULL;


// Constructors
Calculator::Calculator()
{
}


Calculator::Calculator(int pid) : Parameter(pid)
{
}


Calculator::Calculator(int pid, char* data_string, char* param_name)
{
  setData(pid, data_string, param_name);
}


void
Calculator::initClass(Model* mdl)
{
  Calculator::model = mdl;
  Calculator::last_id = 0;
}


// Calculator section  specific output-method for Solver input file
ostream&
Calculator::output_sif(ostream& out, short indent_size, short indent_level, SifOutputControl& soc)
{
  UserInterface* gui = (UserInterface*)model->getGui();

  static char gui_value_buffer[1025];
  static char value_buffer[1025];
  static char variable_buffer[1025];

  char QM = '\"';

  ParameterField* order_fld = getFieldBySifName("Solving Order");
  
  // Parameter type and order number as the id!
  if (soc.outputType) {

    LibFront::output_string(out, indent_size, indent_level++, getSifName(), false);
    
    out << ' ';
    
    // Solving order as id!
    if (soc.outputId) {
      if (order_fld != NULL)
        out << order_fld->getDataStrings()[0];
      else
        out << id;
    }

    out << endl;

  } else {
    indent_level++;
  }

  char* fld_name;

  // Fields
  for (short i = 0; i < nofFields; i++) {

    ParameterField* pf = fields[i];
    char* fld_name = (char*)fields[i]->getSifName();

    // Check that field exists and it should be output
    if ( !pf->isActiveInstance() ||
         pf->getNofDataStrings() == 0 ||
         ( !soc.outputAll && !pf->isSifOutputField() )
       )
      continue;

    // We don't output SOLVING_ORDER field (it is used as the id!)
    if ( LibFront::in(fld_name, "Solving Order") ) {
      continue;
    }

    pf->output_sif(out, indent_size, indent_level, soc.sectionName);

  }

  return out;
}


void
Calculator::setName(char* param_name)
{
  Parameter::setName(param_name, "Calculator");
}
