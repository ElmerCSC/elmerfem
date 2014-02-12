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
Module:     ecif_parameterField.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  
 
Abstract:   Implementation. 

************************************************************************/

#include "ecif_func.h"
#include "ecif_model.h"
#include "ecif_parameterField.h"


Model* ParameterField::model = NULL;

// ParameterField class

// Constructors
ParameterField::ParameterField()
{
  fieldInfo = new ParameterFieldInfo();
  init();
}


ParameterField::ParameterField(ParameterFieldInfo* pf_info,
                               char** var_name_buffers,
                               short dim1, short dim2, 
                               short nof_variables,
                               short nof_strings, char** data_strings)
{
  init();
  fieldInfo = pf_info;
  dimension1 = dim1;
  dimension2 = dim2;
  nofVariables = nof_variables;
  nofDataStrings = nof_strings;

  if ( nofDataStrings > 0 ) {

    dataStrings = new char*[nofDataStrings];

    for (short i = 0; i < nofDataStrings; i++) {
      short len = strlen(data_strings[i]);
      dataStrings[i] = new char[1 + len];
      strcpy(dataStrings[i], data_strings[i]);
      dataStrings[i][len] = '\0';
    }
  }

  if (nofVariables > 0) {
    variableNames = new char*[nofVariables];
    for (short i = 0; i < nofVariables; i++) {
      variableNames[i] = new char[1 + strlen(var_name_buffers[i])];
      strcpy(variableNames[i], var_name_buffers[i]);
      variableNames[i][strlen(var_name_buffers[i])] = '\0';
    }
  }

}


ParameterField::~ParameterField() 
{
  delete_data();
  for (short i = 0; i < nofVariables; i++) {
    delete[] variableNames[i];
  }
  delete[] variableNames;

  delete fieldInfo;
}


void
ParameterField::delete_data()
{
  if (dataStrings != NULL) {
    for (short i = 0; i < nofDataStrings; i++) {
      delete[] dataStrings[i];
    }
    delete[] dataStrings;
    dataStrings = NULL;
  }

  nofEntries = 0;
  nofDataStrings = 0;
}


void
ParameterField::init()
{
  dataStrings = NULL;
  nofDataStrings = 0;
  dimension1 = 0;
  dimension2 = 0;
  isActive = true;
  nofEntries = 0;
  nofVariables = 0;
  variableNames = NULL;
}


ostream&
ParameterField::output_sif(ostream& out, short indent_size, short indent_level, 
                           const char* section_name, bool output_equal_sign)
{
  // Empty or non-output field
  if ( nofDataStrings == 0 ||
       !fieldInfo->outputSif
     )
    return out;

  // Logical False is not output
  if ( !fieldInfo->isArray && 
       nofDataStrings == 1 &&
       !fieldInfo->alwaysOutput && 	
       LibFront::ncEqual("Logical", (char*)fieldInfo->valueType) &&
       LibFront::ncEqual("False", dataStrings[0])
     )
    return out;

  // String "None" is not output
  if ( !fieldInfo->isArray && 
       nofDataStrings == 1 &&
       !fieldInfo->alwaysOutput && 	
       LibFront::ncEqual("String", (char*)fieldInfo->valueType) &&
       LibFront::ncEqual("None", dataStrings[0])
     )
    return out;

  short is = indent_size;
  short il = indent_level;

  bool type_printed = false;

  LibFront::indent(out, is, il);

  // Field name (Density etc.)
  // ----------
  output_sif_name(out, false);

  // Dimension and type for arries (Gravity(4) etc.)
  // -----------------------------
  if (fieldInfo->isArray) {
    output_sif_size(out);
  }

  // Output equal sign
  // -----------------
  //
  // NOTE: Inlcude command does not use equal-sign!
  //
  if ( output_equal_sign &&
       !LibFront::ncEqual("Include", fieldInfo->sifName) 
     ) {
    out << " = ";
  } else {
    out << " ";
  }

  // Variable names the parametrer is dependent of (Temperature etc.)
  // ---------------------------------------------
  if (nofVariables > 0) {
    output_sif_variableNames(out);
  }

  // Data type  (Real etc.)
  // ---------
  //
  if ( !model->getSolverKeywordTypeGiven(section_name, getSifName()) &&
       fieldInfo->outputSifType
     ) {

    // Out type on separate line if variables or procedure
    if ( nofVariables > 0 || fieldInfo->isProcName ) {
      out << endl;
      LibFront::indent(out, is, ++il);
    }

    output_sif_type(out);
    type_printed = true;
  }

  // Field data
  // ----------
  output_sif_data(out, is, ++il, type_printed, true);
  
  return out;
}


ostream&
ParameterField::output_sif_data(ostream& out, short indent_size, short indent_level,
                                bool type_printed, bool add_eol)
{
  // Sif continuation mark (\)
  char cc[] = " \\";

  // If data source is a procedure or file name
  if (fieldInfo->isProcName || fieldInfo->isFileName || fieldInfo->isQuoted) {
    output_sif_data_as_name(out, indent_size, indent_level, type_printed);
    return out;
  }

  // If no "real" data
  if ( dataStrings == NULL )
    return out;

  // If this is "Variable" data, start from  a new line after data type name
  //
  // NOTE: Data type is also always used for tables in Sif!
  //
  if (nofVariables > 0) {
    out << endl;
    if ( !type_printed ) {
      LibFront::indent(out, indent_size, indent_level);
      output_sif_type(out, true);
    }

  }

  // Output "normal" data
  for (int i = 0; i < nofDataStrings; i++) {
    
    // If no variable and table, each row at separate line
	  if (nofVariables == 0 && dimension2 > 1) {
		  out << cc << endl;
	  }

    // If this is "Variable" data or table, extra indent before data lines
	  if (nofVariables > 0 || dimension2 > 1) {
        LibFront::indent(out, indent_size, indent_level + 1);
	  }

    char* data_str = dataStrings[i];

    // If this a Matc-definition, remove the extra $-sign from the
    // beginning
    if ( data_str != NULL && data_str[0] == '$' ) {
      data_str++;
    }

    // Replace character '_' (underline) with space in the output
    // Replace character ' (quote) with " (double quote)
    //
    for (short j = 0; j < strlen(data_str); j++) {
      const char c = data_str[j];
      if (c == '_') out << ' ';
      else if ( c == '\'' ) out << '\"';
      else out << c;
    }
      
    // If this is "Variable" data, set each entry to own line
    if (nofVariables > 0)
      out << endl;
  }

  // "Variable" data is limited by "End",
  if (nofVariables > 0) {
    LibFront::indent(out, indent_size, indent_level);
    out << "End";
    out << endl;

  // Otherwise we just put endline
  } else if (add_eol) {
    out << endl;
  }

  return out;
}


// Data is filename or procedure name
ostream&
ParameterField::output_sif_data_as_name(ostream& out, short indent_size, short indent_level,
                                        bool type_printed, bool add_eol)
{
  char QM = '\"'; // quote mark 

  // Procedure 
  if (fieldInfo->isProcName) {

    // If type not given, move to new line before printing
    if ( !type_printed ) {
      out << endl;
      LibFront::indent(out, indent_size, indent_level);
    }
    
    // Check if we have a Matc-proc, it is output in a speliacal way!
    bool is_matc_proc = false;

    if ( dataStrings[1][0] == '$' ) {
      is_matc_proc = true;
    }

    if ( !is_matc_proc ) {
      out << "Procedure ";
      // Module-name, Procedure-name (don't add quotes!)
      if (nofDataStrings > 1) {
        out << QM << dataStrings[0] << QM
            << " "
            << dataStrings[1];
      // Only Procedure-name 
      } else {
        out << QM << QM
            << " "
            << dataStrings[0];
      }
    } else {
      char* tmp = dataStrings[1];
      // Skip firts "$MATC " part, so that we can set quotes around
      // function name
      tmp += 6;
      out << "MATC ";
      out << QM << tmp << QM;
    }
  }

  // File (NOTE: File will get "File" keyword from data type)
  else if (fieldInfo->isFileName) {
    for (int i = 0; i < nofDataStrings; i++) {
      out << QM << dataStrings[i] << QM;
      if (i < nofDataStrings - 1)
        out << " ";
    }
  }

  // Quoted String
  else if (fieldInfo->isQuoted) {
      out << QM << dataStrings[0] << QM;
  }

  // Unknown type
  else {
    return out;
  }

  if (add_eol) {
    out << endl;
  }
 
  return out;

}


ostream&
ParameterField::output_sif_name(ostream& out, bool add_eol)
{
  const char* name = fieldInfo->sifName;

  return LibFront::output_string(out, 0, 0, name, add_eol);
}


ostream&
ParameterField::output_sif_size(ostream& out, bool add_eol)
{
  out << "(" << dimension1;

  if (dimension2 > 1) 
    out << "," << dimension2;

  out << ")";

  if ( add_eol ) {
    //out << endl;
  }

  return out;
}


ostream&
ParameterField::output_sif_type(ostream& out, bool add_eol)
{
  // Data type name (Real etc.)
  out << fieldInfo->valueType << " ";

  if ( add_eol ) {
    out << endl;
  }

  return out;
}



ostream&
ParameterField::output_sif_variableNames(ostream& out, bool add_eol)
{
  int i;

  //--Check that some non-"none" variables defined!
  bool has_variables = false;

  for (i = 0; i < nofVariables; i++) {
    
    if ( !LibFront::ncEqual(variableNames[i], "none") ) {
      has_variables = true;
      break;
    }
  }

  if ( !has_variables) {
    return out;
  }

  //--Print all non-"none" variable names
  for (i = 0; i < nofVariables; i++) {

    if ( LibFront::ncEqual(variableNames[i], "none") ) {
      continue;
    }

    out << "Variable ";
    out << variableNames[i];

    if ( add_eol ) {
      out << endl;
    }
  }

  return out;
}


