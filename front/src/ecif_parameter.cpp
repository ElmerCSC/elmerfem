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
Module:     ecif_parameter.cpp
Language:   C++
Date:       15.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Implementation

************************************************************************/

#include "ecif_body.h"
#include "ecif_control.h"
#include "ecif_model.h"
#include "ecif_parameter.h"
#include "ecif_parameterField.h"
#include "ecif_userinterface.h"

extern char PARAMETER_BUFFER[];

Model* Parameter::model = NULL;


// Constructors
Parameter::Parameter()
{
  applyCount = 0;
  changed = false;
  copied = false;
  dataString = NULL;
  dataString_previous = NULL;
  fields = NULL;
  id = 0;
  name = NULL;
  nofFields = 0;
  nofPreviousFields = 0;
  parentId = NO_INDEX;
  parentEmfTag = NO_INDEX;
  parentEmfType = OT_NONE;
  subParentEmfTag = NO_INDEX;
  subParentEmfType = OT_NONE;
  previousFields = NULL;
  status = STATUS_OK;
  updateFlag = false;
}


Parameter::Parameter(int p_id)
{
  applyCount = 0;
  changed = false;
  copied = false;
  dataString = NULL;
  dataString_previous = NULL;
  fields = NULL;
  id = p_id;
  name = NULL;
  nofFields = 0;
  nofPreviousFields = 0;
  parentId = NO_INDEX;
  parentEmfTag = NO_INDEX;
  parentEmfType = OT_NONE;
  subParentEmfTag = NO_INDEX;
  subParentEmfType = OT_NONE;
  previousFields = NULL;
  status = STATUS_OK;
  updateFlag = false;
}


Parameter::~Parameter() 
{
  delete [] dataString;
  delete [] dataString_previous;
  delete [] name;
  delete_fields();
  delete_previousFields();
}


void
Parameter::checkLastId(int pid)
{
  int lid = getLastId();

  if ( lid  < pid ) {
    setLastId(pid);
  }
}
 

void
Parameter::create_fields()
{
  UserInterface* gui = (UserInterface*)model->getGui();

  char data_buffer[1025];
  short data_buffer_len = 1025;
  ParameterFieldInfo* pf_info;
  short dim1, dim2, nof_variables;
  char data_type;                   // "=" Numeric, ":" File name, "." Proc name
  bool is_file_name, is_proc_name, is_quoted;
  bool is_inactive;

  short i;

  char* name_buffer = new char[MAX_PARAMETER_FIELD_NAME_LENGTH];
  char* field_name_buffer = new char[MAX_PARAMETER_FIELD_NAME_LENGTH];
  char* name_part_buffer = new char[MAX_PARAMETER_FIELD_NAME_LENGTH];
  char* index_part_buffer = new char[MAX_PARAMETER_FIELD_NAME_LENGTH];

  char** var_name_buffers = new char*[MAX_NOF_PARAMETER_VARIABLES];
  for (i = 0; i < MAX_NOF_PARAMETER_VARIABLES; i++) {
    var_name_buffers[i] = new char[MAX_PARAMETER_FIELD_NAME_LENGTH];
  }

  char** data_strings = new char*[MAX_NOF_PARAMETER_DATA_STRINGS];
  for (i = 0; i < MAX_NOF_PARAMETER_DATA_STRINGS; i++) {
    data_strings[i] = NULL;
  }

  // Count nof fields in the parameter set
  istrstream strm1(dataString);

  nofFields = 0;

  while (!strm1.eof()) {
    strm1.getline(PARAMETER_BUFFER, PARAMETER_BUFFER_LEN, PARAMETER_FIELD_SEP);
    nofFields++;
  }

  // Allocate parameter fields
  fields = new ParameterField*[nofFields];

  for (i = 0; i < nofFields; i++) {
    fields[i] = NULL;
  }

  istrstream strm2(dataString);

  // Read the parameter fields
  short counter = 0;
  while (!strm2.eof()) {

    //----Extract one parameter field
    strm2.getline(PARAMETER_BUFFER, PARAMETER_BUFFER_LEN, PARAMETER_FIELD_SEP);
    
    //--This is the stream we use for reading!
    istrstream strm(PARAMETER_BUFFER);

    dim1 = dim2 = 0;
    nof_variables = 0;

    //--Read parameter name, type and possible (d1 d2) size info
    extractFieldInfo(strm, field_name_buffer, var_name_buffers, data_type,
                     dim1, dim2, nof_variables, is_inactive);


    //--Form parameter field info for the field

    // Check if this is an indexed field like: Diffusivity(Oxygen)
    bool has_pre_index = false;
    bool has_post_index = false;

    strstream fn_strm;

    char* field_name = field_name_buffer;;

    char* index_name = NULL;

    if ( field_name_buffer[0] == INDEX_PRE_SEPARATOR ) {
      has_pre_index = true;
      // Jump over first '<' before reading index-name
      char* tmp = field_name_buffer;
      fn_strm << ++tmp;
      fn_strm.getline(index_part_buffer, MAX_PARAMETER_FIELD_NAME_LENGTH, INDEX_POST_SEPARATOR);
      fn_strm.getline(name_part_buffer, MAX_PARAMETER_FIELD_NAME_LENGTH);
      index_name = index_part_buffer;
    }

    else if ( field_name_buffer[strlen(field_name_buffer) - 1] == INDEX_POST_SEPARATOR ) {
      has_post_index = true;
      fn_strm.getline(name_part_buffer, MAX_PARAMETER_FIELD_NAME_LENGTH, INDEX_PRE_SEPARATOR);
      fn_strm.getline(index_part_buffer, MAX_PARAMETER_FIELD_NAME_LENGTH, INDEX_POST_SEPARATOR);
      index_name = index_part_buffer;
    }


    // Get parameter field info by array and field name
    pf_info = model->getParameterFieldInfo(getArrayName(), field_name);

    // If an indexed field was found, store info about it
    if ( pf_info != NULL && (has_pre_index || has_post_index) ) {

      if (has_pre_index) {
        pf_info->isPreIndexed = true;
      } else {
        pf_info->isPostIndexed = true;
      }

       // Index name
       update_dyna_string(pf_info->guiIndex, index_name);
    
    // Normal field
    } else {

      //--ERROR if unknown field
      if (pf_info == NULL ) {
        pf_info = new ParameterFieldInfo();
        update_dyna_string(pf_info->guiName, "UNKNOWN_FIELD");
        update_dyna_string(pf_info->valueType, "UNKNOWN DATATYPE");
      }
    }

    // Check if we have a file-name or procedure name as data
    is_file_name = false;
    is_proc_name = false;

    if (data_type == FILE_NAME_INDICATOR) {
      is_file_name = true;
    } else if (data_type == PROC_NAME_INDICATOR) {
      is_proc_name = true;
    }

    pf_info->isFileName |= is_file_name;
    pf_info->isProcName |= is_proc_name;

    //--Read all DATA values "1.0;;2.0..." or "conduction;convection;..."
    short nof_strings = 0;
    bool append_data = false;

    while (!strm.eof()) {

      data_buffer[0] = '\0';

      strm.getline( data_buffer, data_buffer_len, PARAMETER_DATA_SEP);
      int len = strlen(data_buffer);
      
      // NOTE: Matc-procedures do not have any Library name
      // so the first item can (and must) be empty!!!
      //
      if (!is_proc_name && len == 0)
        continue;
      
      // Data still continues in the stream
      if (append_data) {
        int old_len = strlen(data_strings[nof_strings]);
        char* tmp = new char[old_len + len + 3];
        strcpy(tmp, data_strings[nof_strings]);
        tmp[old_len] = ' ';
        strcpy(tmp + old_len + 1, data_buffer);
        delete data_strings[nof_strings];
        data_strings[nof_strings] = tmp;
        data_strings[old_len + len + 2] = '\0';

      // New data
      } else {
        delete[] data_strings[nof_strings];
        data_strings[nof_strings] = new char[1 + len];
        strcpy(data_strings[nof_strings], data_buffer);
        data_strings[nof_strings][len] = '\0';
      }

      // If current data was not read completely
      if (len == data_buffer_len - 1  &&
          !strm.eof()                 && 
          PARAMETER_DATA_SEP != strm.peek() ) {
        append_data = true;

      // Data was read completely, inrease array index for the next data
      } else {
        append_data = false;

        if ( nof_strings < MAX_NOF_PARAMETER_DATA_STRINGS - 1 ) {
          nof_strings++;
        }
      }

    } // end read values

    // Set field sizes if they were not given in the input

    // Dim1
    // ====
    if ( dim1 == 0 ) {

      //--Array data
      if ( pf_info->isArray ) {
        
        //-String, dim1 alaways 1
        if ( LibFront::in(pf_info->valueType, "string") ) {
          dim1 = 1;

        //-Numeric
        } else {
          dim1 = countNumbersInString(data_strings[0]) - nof_variables;
        }

      //--Scalar
      } else {
        dim1 = 1;
      }
    }

    // Dim2
    // ====
    if ( dim2 == 0 ) {

      //--Array data
      if ( pf_info->isArray ) {
        dim2 = nof_strings;

      //--Scalar
      } else {
        dim2 = 1;
      }
    }

    pf_info->isArray |= (dim1 > 1);

    // Create a new parameter field from the input values
    // ==================================================
    ParameterField* pf;
    pf = new ParameterField(pf_info, var_name_buffers,
                            dim1, dim2, nof_variables,
                            nof_strings, data_strings);

    fields[counter++] = pf;

    // If field is (intentionally!) empty or
    // some alternative (scalar, table, proc) for the field
    // is already stored, make it inactive (we just store it,
    // we do not eg. output it to sif-file!)
    if ( is_inactive || nof_strings == 0 || fieldInstanceExists(pf) ) {
      pf->setIsActive(false);
    }

  } // End read parameter fields

  delete[] name_buffer;
  delete[] field_name_buffer;
  delete[] name_part_buffer;
  delete[] index_part_buffer;

  for (i = 0; i < MAX_NOF_PARAMETER_VARIABLES; i++) {
    delete[] var_name_buffers[i];
  }
  delete[] var_name_buffers;

  for (i = 0; i < MAX_NOF_PARAMETER_DATA_STRINGS; i++) {
    delete[] data_strings[i];
  }
  delete[] data_strings;

}


void
Parameter::extractFieldInfo(istrstream& strm,
                            char* field_name_buffer,
                            char** var_name_buffers,
                            char& data_type,
                            short& dim1, short& dim2,
                            short& nof_variables,
                            bool& is_inactive)
{
  is_inactive = false;

  // Possible inactive marker ("-")
  if ( strm.peek() ==  '-' ) {
    extractInactiveMarker(strm);
    is_inactive = true;
  }

  // Possible description in the beginning (old format!)
  if ( strm.peek() ==  '(' ) {
    extractFieldSizeAndVars(strm, var_name_buffers, dim1, dim2, nof_variables);
  }

  extractFieldNameAndType(strm, field_name_buffer, data_type);

  // Possible description after the name and before the data (new format!)
  if ( strm.peek() ==  '(' ) {
    extractFieldSizeAndVars(strm, var_name_buffers, dim1, dim2, nof_variables);
  }
}


void
Parameter::extractInactiveMarker(istrstream& strm)
{
  // Read inactive marker from the beginning of the field ("-")
  strm.get();
}


void
Parameter::extractFieldNameAndType(istrstream& strm,
                                   char* field_name_buffer,
                                   char& data_type)
{
  // Read field name and type (== =: =.)
  strm.getline(field_name_buffer, MAX_PARAMETER_FIELD_NAME_LENGTH, '=');
  data_type = strm.get();
}


void
Parameter::extractFieldSizeAndVars(istrstream& strm,
                                    char** var_name_buffers,
                                    short& dim1, short& dim2,
                                    short& nof_variables)
{
  dim1 = dim2 = 0;
  nof_variables = 0;
  const int buflen = 81;
  char buffer[1 + buflen];

  // Read dimension desrciption 
  // Format is: "( [variable name1 name2 ; dim1 dim2 )"
  strm.get(); // Pick of "("
  strm.getline(buffer , buflen, ')' );

  bool has_variables = false;

  // Check if there is a variable definition
  for (int i = 0; i < strlen(buffer); i++) {
    if ( isalpha(buffer[i]) ) {
      has_variables = true;
      break;
    }
  }

  strstream strm1;
  strm1 << buffer << ends;
  strm1.getline(buffer, buflen, PARAMETER_DATA_SEP);

  // We have either variable defs or size description
  if ( strlen(buffer) > 0 ) {
    strstream strm2;
    strm2 << buffer << ends;

    if (has_variables) {
      strcpy(var_name_buffers[nof_variables++], strm2.str());
    } else {
      strm2 >> dim1;
      if (!strm2.eof())
        strm2 >> dim2;
    }
  }

  strm1.getline(buffer, 81, PARAMETER_DATA_SEP);

  // If we now have the size description
  if (strlen(buffer) > 0) {
    strstream strm3;
    strm3 << buffer << ends;
    strm3 >> dim1;
    if (!strm3.eof())
      strm3 >> dim2;
  }
}


bool
Parameter::fieldInstanceExists(ParameterField* param_field)
{
  // Check if argument field name matches with any other existing
  for (int i = 0; i < nofFields; i++) {
 
    ParameterField* pf = fields[i];
    
    if ( pf == NULL || pf == param_field ) {
      continue;
    }

    if ( 0 != strcmp(param_field->getGuiName(), pf->getGuiName()) ) {
      continue;
    }

    const char* i1 = param_field->getGuiIndex();
    const char* i2 = pf->getGuiIndex();

    // No indexed fields, must be the same!
    if ( i1 == NULL && i2 == NULL ) {
      return true;
    }
    
    // One, but no both is NULL --> different fields
    if ( i1 == NULL || i2 == NULL ) {
      continue;
    }

    // Different index names --> different fields
    if ( 0 != strcmp(i1, i2) ) {
      continue;
    }

    // Same index names and index types --> must be the same!
    if ( ( param_field->isPreIndexed() && pf->isPreIndexed() ) ||
         ( param_field->isPostIndexed() && pf->isPostIndexed() )
         ) {
      return true;
    }
  }

  // No matching field found
  return false;
}


void
Parameter::delete_fields()
{
  for (short i = 0; i < nofFields; i++) {
    delete fields[i];
  }

  delete[] fields;
  nofFields = 0;
  fields = NULL;
}


void
Parameter::delete_previousFields()
{
  for (short i = 0; i < nofPreviousFields; i++) {
    delete previousFields[i];
  }
  delete[] previousFields;
  nofPreviousFields = 0;
  previousFields = NULL;
}


// If parameter has a "gui-name"-field (like "HEAT_EQUATION") defined,
// return the pointer to it, otherwise return NULL
ParameterField*
Parameter::getFieldByGuiName(const char* name, bool only_active)
{
  for (short i = 0; i < nofFields; i++) {

    ParameterField* pf = fields[i];

    if ( LibFront::in(name, pf->getGuiName()) ) {

      // Check activity (minus mark in the input data like |-TEMPERATUE=1000.0) 
      if ( only_active && !pf->isActiveInstance() ) {
        return NULL;
      } else {
        return pf;
      }
    }
  }

  return NULL;
}


// If parameter has a "gui-name"-field (like "HEAT_EQUATION") defined,
// return the pointer to it, otherwise return NULL
ParameterField*
Parameter::getFieldByGuiName(const char* name, const char* index, bool is_pre_indexed, bool only_active)
{
  if ( index == NULL ) {
    return getFieldByGuiName(name, only_active);
  }

  for (short i = 0; i < nofFields; i++) {

    ParameterField* pf = fields[i];

    if ( LibFront::in(name, pf->getGuiName()) &&
         LibFront::in(index, pf->getGuiIndex()) &&
         ( (is_pre_indexed && pf->isPreIndexed()) ||
           (!is_pre_indexed && pf->isPostIndexed())
         )
         ) {

      // Check activity (minus mark in the input data like |-TEMPERATUE=1000.0) 
      if ( only_active && !pf->isActiveInstance() ) {
        return NULL;
      } else {
        return pf;
      }
    }
  }

  return NULL;
}

// If parameter has a "sif-name"-field (like "Heat Equation") defined,
// return the pointer to it, otherwise return NULL
ParameterField*
Parameter::getFieldBySifName(const char* fname, bool only_active)
{
  for (short i = 0; i < nofFields; i++) {

    ParameterField* parameter_field = fields[i];

    if ( LibFront::in(fname, parameter_field->getSifName()) ) {

      // Check activity (minus mark in the input data like |-TEMPERATUE=1000.0) 
      if ( only_active && !parameter_field->isActiveInstance() ) {
        return NULL;

      } else {
        return parameter_field;
      }
    }
  }
  return NULL;
}




// Check is the parameter has given boolean field 
// and get its value if defined
// NOTE: Data is suppoused to be scalar!
bool
Parameter::getFieldValueBySifName(const char* name, bool& value, bool only_active)
{
  ParameterField* pf = getFieldBySifName(name, only_active);

  if (pf == NULL)
    return false;

  char** data_strings = pf->getDataStrings();

  if ( 0 == strcmp(data_strings[0], "True") )
    value = true;
  else
    value = false;

  // Value found
  return true;
}


// Check is the parameter has given boolean field 
// and get its value if defined
// NOTE: Data can be a vector (but not array)
bool
Parameter::getFieldValueBySifName(const char* name, int& nof_values, bool*& values, bool only_active)
{
  static char str[1024];

  ParameterField* pf = getFieldBySifName(name, only_active);
  
  if (pf == NULL) {
    nof_values = 0;
    values = NULL;
    return false;
  }
  
  short d1, d2, nof_vars;
  pf ->getDataDimension(d1, d2, nof_vars);

  nof_values = d1;

  values = new bool[nof_values];

  char** data_strings = pf->getDataStrings();

  strstream strm;
  strm << data_strings[0] << ends;

  for (int i = 0; i < nof_values; i++) {
  
    strm >> str;

    if ( 0 == strcmp(str, "True") )
      values[i] = true;
    else
      values[i] = false;
  }

  // Value found
  return true;
}


// Check is the parameter has given real field 
// and get its value if defined
// NOTE: Data is suppoused to be scalar!
bool
Parameter::getFieldValueBySifName(const char* name, int max_buffer_len, char* buffer, bool only_active)
{
  ParameterField* pf = getFieldBySifName(name, only_active);

  if (pf == NULL)
    return false;

  char** data = pf->getDataStrings();

  strncpy(buffer, (char*)data[0], max_buffer_len);

  buffer[max_buffer_len] = '\0';

  // Value found
  return true;
}


// Check is the parameter has given integer field 
// and get its value if defined
// NOTE: Data is suppoused to be scalar!
bool
Parameter::getFieldValueBySifName(const char* name, int& value, bool only_active)
{
  ParameterField* pf = getFieldBySifName(name, only_active);

  if (pf == NULL)
    return false;

  char** data = pf->getDataStrings();

  if ( !LibFront::isNumber(data[0]) )
    return false;

  value = atol(data[0]);

  // Value found
  return true;
}


// Check is the parameter has given integer field 
// and get its value if defined
// NOTE: Data can be vector!
bool
Parameter::getFieldValueBySifName(const char* name, int& nof_values, int*& values, bool only_active)
{
  ParameterField* pf = getFieldBySifName(name, only_active);

  if (pf == NULL) {
    nof_values = 0;
    values = NULL;
    return false;
  }
  
  short d1, d2, nof_vars;
  pf ->getDataDimension(d1, d2, nof_vars);

  nof_values = d1;

  values = new int[nof_values];

  char** data_strings = pf->getDataStrings();

  strstream strm;
  strm << data_strings[0];

  for (int i = 0; i < nof_values; i++) {
    strm >> values[i];
  }

  // Value found
  return true;
}


// Check is the parameter has given real field 
// and get its value if defined
// NOTE: Data is suppoused to be scalar!
bool
Parameter::getFieldValueBySifName(const char* name, double& value, bool only_active)
{
  ParameterField* pf = getFieldBySifName(name, only_active);

  if (pf == NULL)
    return false;

  char** data = pf->getDataStrings();

  if ( !LibFront::isNumber(data[0]) )
    return false;

  value = atof(data[0]);

  // Value found
  return true;
}


// Check is the parameter has given real field 
// and get its value if defined
// NOTE: Data is suppoused to be scalar!
bool
Parameter::getFieldValueBySifName(const char* name, int& nof_values, double*& values, bool only_active)
{
  ParameterField* pf = getFieldBySifName(name, only_active);

  if (pf == NULL) {
    nof_values = 0;
    values = NULL;
    return false;
  }
  
  short d1, d2, nof_vars;
  pf ->getDataDimension(d1, d2, nof_vars);

  nof_values = d1;

  values = new double[nof_values];

  char** data_strings = pf->getDataStrings();

  strstream strm;
  strm << data_strings[0];

  for (int i = 0; i < nof_values; i++) {
    strm >> values[i];
  }

  // Value found
  return true;
}



// Check if there is any change in the data
void
Parameter::getValueState(bool& value_has_changed)
{
  if (dataString == NULL && dataString_previous == NULL) {
    value_has_changed = false;
    return;
  }

  if (dataString == NULL || dataString_previous == NULL) {
    value_has_changed = true;
    return;
  }

  if ( 0 == strcmp(dataString, dataString_previous ) ) {
    value_has_changed = false;
    return;
  }

  value_has_changed = true;
}


// Find if the parameter's "field_name"-field was changed from the previous value
void
Parameter::getFieldValueStateBySifName(const char* name, bool& value_is_changed)
{
  // If no change in the wholw parameter
  // we can answer quickly1
  if (!changed) {
    value_is_changed = false;
    return;
  }

  // Ok, we have to really check the field 
  bool has_value;
  getFieldValueStateBySifName(name, has_value, value_is_changed);

  changed = false;  // We have been checked for the change
}


// Find if the parameter has "field_name"-field defined and check if the
// value was changed from the previous value
void
Parameter::getFieldValueStateBySifName(const char* name, bool& has_value, bool& value_is_changed)
{
  ParameterField* pf = getFieldBySifName(name);
  ParameterField* pf_previous = getPreviousFieldBySifName(name);

  // Not at all the parameter field --> it is not "changed"
  if (pf_previous == NULL && pf == NULL ) {
    has_value = false;
    value_is_changed = false;

    return;
  }
  
  // Not more the parametrer field --> it is "changed"
  if (pf == NULL ) {
    has_value = false;
    value_is_changed = true;
    return;
  }

  // New parameter field --> it is "changed"
  if (pf_previous == NULL) {
    has_value = true;
    value_is_changed = true;
    return;
  }
 
  // Ok, we have the field, we test for a change
  has_value = true;
  value_is_changed = true;
  // Compare previous and current
  int nof_ds_previous = pf_previous->getNofDataStrings();
  int nof_ds_current = pf->getNofDataStrings();

  // If nof data strings has changed
  // NOTE: This is certainly not reliable test !!!***!!!
  if (nof_ds_previous != nof_ds_current) {
    return;
  }

  // If string-pairs are not similar
  // NOTE: again not reliable test !!!***!!!
  char** ds_previous = pf_previous->getDataStrings();
  char** ds_current = pf->getDataStrings();
  for (int i = 0; i < nof_ds_current; i++) {
      if ( !strcmp(ds_previous[i], ds_current[i]) )
        return;
  }

  // Value was not (propably :-) changed!
  value_is_changed = false;
}


int
Parameter::getParentEmfTag()
{
  return parentEmfTag;
}


objectType
Parameter::getParentEmfType()
{
  return parentEmfType;
}


int
Parameter::getSubParentEmfTag()
{
  return subParentEmfTag;
}


objectType
Parameter::getSubParentEmfType()
{
  return subParentEmfType;
}



// If parameter previous values has a "gui-name"-field (like "HEAT_EQUATION") defined,
// return the pointer to it, otherwise return NULL
ParameterField*
Parameter::getPreviousFieldByGuiName(const char* fname)
{
  for (short i = 0; i < nofPreviousFields; i++) {

    ParameterField* parameter_field = previousFields[i];

    if ( LibFront::in(fname, parameter_field->getGuiName()) )
      return parameter_field;
  }

  return NULL;
}


// If parameter previous values has a "sif-name"-field (like "Heat Equation") defined,
// return the pointer to it, otherwise return NULL
ParameterField*
Parameter::getPreviousFieldBySifName(const char* fname)
{
  for (short i = 0; i < nofPreviousFields; i++) {

    ParameterField* parameter_field = previousFields[i];

    if ( LibFront::in(fname, parameter_field->getSifName()) )
      return parameter_field;
  }

  return NULL;
}


// Check is the parameter has a "field_name"-field defined
// and that its value is "value"
bool
Parameter::hasFieldValueBySifName(const char* name, char* value)
{
  ParameterField* pf = getFieldBySifName(name);

  if (pf == NULL)
    return false;

  int nof_ds = pf->getNofDataStrings();
  char** data_strings = pf->getDataStrings();

  for (int i = 0; i < nof_ds; i++) {
    if ( 0 == strcmp(value, data_strings[i]) ) {
      return true;
    }
  }

  return false;
}


// Check is the parameter is active
// Note: by default a parameter is active, unless
// it is flagged inactive with the EFN_ACTIVE field
//
bool
Parameter::IsActive() 
{
  ParameterField* pf = getFieldByGuiName("ACTIVE");

  // If no field at all ==> is active!
  if ( pf == NULL ) {
    return true;
  }

  char** data = pf->getDataStrings();

  if ( !LibFront::ncEqual((char*)data[0], "true") ) {
    return false;

  } else {
    return true;
  }
}


// Method prints parameter into output stream in Cadinterface format
// Output of parameter name and id number is controlled by: output_type flag
// The output format of the data is controlled by: data_as_string flag
ostream&
Parameter::output_emf(ostream& out, short indent_size, short indent_level,
                      bool output_type,
                      bool data_as_string,
                      bool data_as_fields)
{
  short& is = indent_size;
  short& il = indent_level;
  char* QM = "\""; // quote-mark
 
  if (output_type) {
    LibFront::output_string(out, indent_size, indent_level++, getEmfName(), false);
    out << ' ' << id << endl;
  }

  // Parent info if applicable
  if ( NO_INDEX != getParentEmfTag() )
    LibFront::output_scalar(out, is, il, "Parent", NULL, getParentEmfTag());

  if ( NO_INDEX != getSubParentEmfTag() )
    LibFront::output_scalar(out, is, il, "Sub Parent", NULL, getSubParentEmfTag());

  if ( OT_NONE != getParentEmfType() )
    LibFront::output_scalar(out, is, il, "Parent Type", NULL, model->objectType2Name(getParentEmfType()), true);


  LibFront::output_scalar(out, is, il, "Name", NULL, name, true);

// These are not in use!
// Mve 06.02.00
#if 0
  LibFront::output_scalar(out, is, il, EFN_APPLY_COUNT, applyCount);

  if (attachMode != ECIF_IGNORE) {
    LibFront::output_scalar(out, is, il, EFN_ATTACH_MODE, attachMode);
  }
#endif

  // Data
  //-as string
  if (data_as_string) { 
    LibFront::output_scalar(out, is, il, "Data", NULL, dataString, true);
  }
  //-as separate fields
  if (data_as_fields) {
    for (short i = 0; i < nofFields; i++) {

      fields[i]->output_sif(out, indent_size, indent_level, NULL);

    }
  }

  return out;
}



// Method prints all parameter's fields into output stream.
// Output of parameter (section) type and id number is
// controlled by output_type flag
ostream&
Parameter::output_sif(ostream& out, short indent_size, short indent_level, SifOutputControl& soc)
{
  // Parameter type and id
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
  for (short i = 0; i < nofFields; i++) {

    ParameterField* pf = fields[i];

    // Check that field is active, contains data and it should be output
    if ( !pf->isActiveInstance() || 
         pf->getNofDataStrings() == 0 ||
         (!soc.outputAll && !pf->isSifOutputField())
       )
      continue;
    
    pf->output_sif(out, indent_size, indent_level, soc.sectionName);
    
  }

  return out;
}


ostream&
Parameter::output_sif_name(ostream& out, short indent_size, short indent_level, SifOutputControl& soc,
                           bool quoted)
{
  if (soc.outputName) {
    if ( model->getSolverKeywordTypeGiven(soc.sectionName, "Name") ) {
      LibFront::output_scalar(out, indent_size, indent_level, "Name =", NULL, name, quoted);
    } else {
      LibFront::output_scalar(out, indent_size, indent_level, "Name = String", NULL, name, quoted);
    }
  }

  return out;
}



bool
Parameter::hasDefaultName()
{
  // If name not given or is in default format
  if ( name == NULL && name[0] == '\0' ) {
    return false;
  }

  const char* base_nm = getGuiName();

  if ( base_nm == NULL ) {
    return false;
  }

  strstream strm;

  strm << base_nm << id << ends;

  if ( 0 == strcmp(name, strm.str()) ) {
    return true;
  } else {
    return false;
  }
}


// Method copies current value in to new values
// and marks parameter unchanged
void
Parameter::resetValue()
{
  // make a copy of the current data for setValue()
  int len = strlen(dataString);
  char* copy_string = new char[1 + len];
  strcpy(copy_string, dataString);
  copy_string[len] = '\0';

  copied = false; // We need this to force copy in setValue() !
  setValue(copy_string);
  changed = false; // Now we can mark data not changed!

  delete[] copy_string;
}


// Method sets object attributes.
void
Parameter::setData(int pid, char* data_string, char* param_name)
{
  id = pid;
  setValue(data_string);
  setName(param_name);
}


// Sets parameter data
void
Parameter::setData(int pid, int parent_id, char* data_string, char* param_name)
{
  id = pid;
  updateParentInfo(parent_id);
  setValue(data_string);
  setName(param_name);
}

 
void
Parameter::setName(char* param_name, char* default_name)
{
  // Set default name 
  if ( param_name == NULL || param_name[0] == '\0' ) {
    strstream strm;
    strm << default_name << id << ends;
    update_dyna_string(name, strm.str());

  // Name given
  } else {
    update_dyna_string(name, param_name);
  }
}


void
Parameter::setValue(char* data_string)
{
  // Copy current values to previous values
  if (dataString != NULL && !copied) {

    // Data string
    update_dyna_string(dataString_previous, dataString);

    // Parameter fields
    delete_previousFields();
    previousFields = fields;
    nofPreviousFields = nofFields;

    fields = NULL;
    nofFields = 0;

    copied = true;
  }

  // Store new value string
  update_dyna_string(dataString, data_string);

  // Create new parameter fields
  // NOTE: don't DELETE parameterFields because previousParameterFields
  // are now pointing to this set!!!***!!!
  if ( dataString != NULL && dataString[0] != '\0' ) {
    create_fields();
  }

  if ( (dataString_previous == NULL) ||
       (0 != strcmp(dataString, dataString_previous))
       ) {
    changed = true;
  }
}


// Set new object parent info
//
void
Parameter::updateParentInfo(int parent_id)
{
  parentId = parent_id;
  parentEmfTag = NO_INDEX;
  parentEmfType = OT_NONE;
  subParentEmfTag = NO_INDEX;
  subParentEmfType = OT_NONE;

  ModelObject* obj = model->getModelObjectById(parent_id);

  if ( obj != NULL ) {
    parentEmfTag = obj->Tag();
    parentEmfType = obj->getObjectType();
  }
}


// Normal case when the emf-parent info defines the actual
// parent object
//
// NOTE: GridParameter (due to Layers) and BoundaryConditions (due to Groups)
// behave  differently
//
void
Parameter::updateParentId()
{
  ModelObject* obj = model->getModelObjectByTag(parentEmfType, parentEmfTag);

  if ( obj != NULL ) {
    parentId = obj->Id();
  }
}


