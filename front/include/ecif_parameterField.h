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
Module:     ecif_parameterField.h
Language:   C++
Date:       20.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   A class for one parameter field value (like Velocity etc.).

************************************************************************/

#ifndef _ECIF_PARAMETER_FIELD_
#define _ECIF_PARAMETER_FIELD_

#include "ecif_def.h"


class ParameterField {
public:
  ParameterField();
  ParameterField(ParameterFieldInfo* pf_info,
                 char** var_name_buffers,
                 short dim1, short dim2,
                 short nof_variables,
                 short nof_data_strings, char** data_value_strings);
  ~ParameterField();
  char** getDataStrings() { return dataStrings;}
  void getDataDimension(short& dim1, short& dim2, short& nof_variables)
                          { dim1 = dimension1;
                            dim2 = dimension2;
                            nof_variables = nofVariables;
                          }
  const char* getValueType() {return fieldInfo->valueType;}  // Real Integer File etc.
  const char* getSifName() { return fieldInfo->sifName; };   // "Density", "Solving Order" "Navier-Stokes" etc.
  const char* getGuiName() { return fieldInfo->guiName; };   // "DENSITY", "SOLVING_ORDER" "NAVIER-STOKES" etc.
  const char* getGuiIndex() { return fieldInfo->guiIndex; }; // "oxygen", "some_species" etc.
  int getNofDataStrings() { return nofDataStrings;}
  int getNofEntries() { return nofEntries;}
  short getNofVariables() { return nofVariables;}
  const char** getVariableNames() { return (const char**)variableNames;}
  static void initClass(Model* mdl) {ParameterField::model = mdl;};
  bool isActiveInstance() { return isActive; }
  bool isSifOutputField() { return fieldInfo->outputSif; }
  bool isPostIndexed() { return fieldInfo->isPostIndexed; }
  bool isPreIndexed() { return fieldInfo->isPreIndexed; }
  bool isProcedure() { return fieldInfo->isProcName; }
  ostream& output_sif(ostream& out, short indent_size, short indent_level, const char*section_name, bool ouput_equal_sign = true);
  ostream& output_sif_data(ostream& out, short indent_size, short indent_level, bool type_printed,  bool add_eol = true);
  ostream& output_sif_data_as_name(ostream& out, short indent_size, short indent_level, bool type_printed,  bool add_eol = true);
  ostream& output_sif_name(ostream& out, bool add_eol = false);
  ostream& output_sif_size(ostream& out, bool add_eol = false);
  ostream& output_sif_type(ostream& out, bool add_eol = false);
  ostream& output_sif_variableNames(ostream& out, bool add_eol = false);
  void setIsActive(bool value) { isActive = value; }
  void setSifTypeGiven(bool value) { fieldInfo->sifTypeGiven = value; }

protected:
  static Model* model;
  char** dataStrings;       // data rows (originally separated by ';') in character format
  short dimension1;         // row dimension of the data entry (1 for scalar)
  short dimension2;         // column dimension of the data entry ( 1 for scalar)
  ParameterFieldInfo* fieldInfo;
  bool isActive;            // flags if field is active (for cases where we store scalar, table and proc values for the field)
  int nofDataStrings;      // nof data string entries (';' separated strings in the field)
  int nofEntries;          // nof data entries (entry size: dimension1*dimension2)
  // If parameter is dependent on Temperature, Time or any other (scalar) variable
  short nofVariables;       // Normally 1 <--> one argument variable
  char** variableNames;     // "Temperature", "Time" "Coordinate 1" etc.
  //
  void delete_data();
  void init();
};


#endif

