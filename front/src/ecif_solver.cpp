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
#include "ecif_solver.h"
#include "ecif_userinterface.h"

//Initialize static class variables.
int Solver::last_id = 0;
Model* Solver::model = NULL;


// Constructors
Solver::Solver()
{
}


Solver::Solver(int pid) : Parameter(pid)
{
}


Solver::Solver(int pid, char* data_string, char* param_name)
{
  setData(pid, data_string, param_name);
}


void
Solver::initClass(Model* mdl)
{
  Solver::model = mdl;
  Solver::last_id = 0;
}


// Solver section  specific output-method for Solver input file
ostream&
Solver::output_sif(ostream& out, short indent_size, short indent_level, SifOutputControl& soc)
{
  UserInterface* gui = (UserInterface*)model->getGui();

  static char value_gui_buffer[1025];
  static char value_sif_buffer[1025];
  static char variable_sif_buffer[1025];

  char QM = '\"';

  int id_shift;

  if ( soc.reloadSolverIsOutput ) {
    id_shift = 1;
  } else  {
    id_shift = 0;
  }

  ParameterField* order_fld = getFieldBySifName("Solving Order");

  // Parameter type and order number as the id!
  if (soc.outputType) {

    LibFront::output_string(out, indent_size, indent_level++, getSifName(), false);

    out << ' ';

    // Solving order as id!
    if (soc.outputId) {
      if (order_fld != NULL)
        out << id_shift + atoi(order_fld->getDataStrings()[0]);
      else
        out << id_shift + id;
    }

    out << endl;

  } else {
    indent_level++;
  }

  char* fld_name;

  // Fields
  // ======

  for (short i = 0; i < nofFields; i++) {

    ParameterField* pf = fields[i];

    char* fld_name = (char*)pf->getSifName();

    // Check that field exists and it should be output
    if ( !pf->isActiveInstance() ||
         pf->getNofDataStrings() == 0 ||
         ( !soc.outputAll && !pf->isSifOutputField() )
       ) {
      continue;
    }

    // We don't output SOLVING_ORDER field (it is used as the id!)
    if ( LibFront::in(fld_name, "Solving Order") ) {
      continue;
    }

    // Meshname
    // ========
    // We output Solver meshname only if is it different from the
    // first active mesh
    if ( LibFront::in(fld_name, "Mesh") ) {

      if ( !getFieldValueBySifName("Mesh", 1024, value_sif_buffer) ) {
        continue;
      }

      const char* amesh1 = model->getActiveMeshName(0);

      // NOTE: Strict comparison for dirnames!
      if ( amesh1 == NULL || 0 == strcmp(value_sif_buffer, (char*)amesh1) ) {
        continue;
      }

      // NOTE: "MESHDIR" directory name must be output before the mesh name!

      pf->output_sif_name(out);
      pf->output_sif_type(out);
      out << " " << QM << model->getMeshDirValue() << QM << " ";
      pf->output_sif_data(out, indent_size, 2 + indent_level, true);

      // Mesh Input File
      // ---------------
      int mesh_index = model->getMeshIndex(value_sif_buffer);

      if ( mesh_index != NO_INDEX ) {

        char* mif_file_name = NULL;
        model->getMeshInputFileName(mif_file_name, mesh_index);

        if ( mif_file_name != NULL ) {
          out << endl;
          LibFront::output_scalar(out, indent_size, indent_level, "Mesh Input File", "File", mif_file_name, true);
          out << endl;

          delete[] mif_file_name;
        }
      }

      if ( i < nofFields - 1 ) out << endl;

      continue;
    }

    // Equation with variable name
    // ===========================
    // Equations with 'real' variable name. like Advection Diffusion, are
    // named with the variable like: "Advection Diffusion Equation Oxygen"
    //
    if ( LibFront::in(fld_name, "Equation") &&
         getFieldValueBySifName("Equation", 1024, value_sif_buffer)
       ) {

      // Change "Heat Equation" --> "HEAT_EQUATION" etc.
      model->fieldNameSifToGui(value_sif_buffer, value_gui_buffer);

      // Check that variable name should (and could) be added to
      // the equation name
      if ( gui->getUseVariableNameInEquationName(value_gui_buffer) &&
           getFieldValueBySifName("Variable", 1024, variable_sif_buffer)
         ) {

        // Form: Equation name + Variable name
        strstream strm;
        strm << value_sif_buffer << " " << variable_sif_buffer << ends;

        // Output Equation-field
        if ( model->getSolverKeywordTypeGiven(SIF_SOLVER, "Equation") ) {
          LibFront::output_scalar(out, indent_size, indent_level, "Equation =", NULL, strm.str(), false);
        } else {
          LibFront::output_scalar(out, indent_size, indent_level, "Equation = String", NULL, strm.str(), false);
        }

        if ( i < nofFields - 1 ) out << endl;

        continue;
      }

    }  // If equation name with varaible name like: "Advection Diffusion Equation Oxygen"


    // Other fields (including "no-variable-named" equation name!)
    // ============
    pf->output_sif(out, indent_size, indent_level, soc.sectionName);

  }

  //--Possible redirected 'TargetFields' from Equation panel
  getFieldValueBySifName("Equation", 1024, value_sif_buffer);
  model->outputSolverTargetFields_sif(out, indent_size, indent_level, value_sif_buffer);

  return out;
}


void
Solver::setName(char* param_name)
{
  Parameter::setName(param_name, "Solver");
}




