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
Module:     ecif_parameter.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   An abstract base class for all parameter-types.

************************************************************************/

#ifndef _ECIF_PARAMETER_
#define _ECIF_PARAMETER_

#include "ecif_def.h"
#include "ecif_def_stl.h"


// ****** Parameter class ******
// Abstarct class!
//
class Parameter {
public:
  Parameter();
  Parameter(int p_id);
  virtual ~Parameter();
  virtual void checkLastId(int pid);
  virtual ecif_modelStatus checkStatus() {return status;}
  virtual int getLastId() {return NO_INDEX;}
  int getApplyCount() {return applyCount;}
  virtual const char* getGuiName() { return "UNKNOWN PARAMETER"; }
  virtual const char* getArrayName() { return "UNKNOWN_PARAMETER_ARRAY"; }
  virtual const char* getEmfName() { return "UNKNOWN PARAMETER"; }
  virtual int getParentEmfTag();
  virtual objectType getParentEmfType();
  virtual int getSubParentEmfTag();
  virtual objectType getSubParentEmfType();
  const char* getName() {return name;}
  virtual const char* getSifName() { return "UNKNOWN PARAMETER"; }
  const char* getValue() {return dataString;}
  virtual ParameterField* getFieldByGuiName(const char* name, bool only_active = true);
  virtual ParameterField* getFieldByGuiName(const char* name, const char* index, bool is_pre_indexed, bool only_active = true);
  virtual ParameterField* getFieldBySifName(const char* name, bool only_active = true);
  virtual bool getFieldValueBySifName(const char* name, bool& value, bool only_active = true);
  virtual bool getFieldValueBySifName(const char* name, int& nof_values, bool*& values, bool only_active = true);
  virtual bool getFieldValueBySifName(const char* name, int max_buffer_len, char* buffer, bool only_active = true);
  virtual bool getFieldValueBySifName(const char* name, double& value, bool only_active = true);
  virtual bool getFieldValueBySifName(const char* name, int& nof_values, double*& values, bool only_active = true);
  virtual bool getFieldValueBySifName(const char* name, int& value, bool only_active = true);
  virtual bool getFieldValueBySifName(const char* name, int& nof_values, int*& values, bool only_active = true);
  virtual void getFieldValueStateBySifName(const char* name, bool& has_value, bool& value_is_changed);
  virtual void getFieldValueStateBySifName(const char* name, bool& value_is_changed);
  virtual ecif_parameterType getParameterType() { return ECIF_NOPARAM; }
  virtual int getParentId() { return parentId;};
  virtual ParameterField* getPreviousFieldByGuiName(const char* fname);
  virtual ParameterField* getPreviousFieldBySifName(const char* fname);
  virtual void getValueState(bool& value_has_changed);
  virtual bool getUpdateFlag() { return updateFlag;}
  virtual bool hasDefaultName();
  virtual bool hasFieldValueBySifName(const char* name, char* value);
  int ID() {return id;}
  static void initClass(Model* mdl) {Parameter::model = mdl;};
  bool IsActive();
  virtual bool IsApplied() {return (applyCount > 0);}
  virtual ostream& output_emf(ostream& out, short indent_size, short indent_level,
                              bool output_type,
                              bool data_as_string,
                              bool data_as_fields);
  ostream& output_field1_emf(ostream& out, short indent_size, short indent_level,
                             char* name, char* type, void* data, char* sep = NULL);
  virtual ostream& output_sif(ostream& out, short indent_size, short indent_level, SifOutputControl& soc);
  virtual ostream& output_sif_name(ostream& out, short indent_size, short indent_level, SifOutputControl& soc, bool quoted = true);
  virtual ostream& outputSolverTargetFields_sif(ostream& out, short indent_size, short indent_level, const char* source_eq_name, NameSet& targetFieldNames) { return out; }
  virtual void resetValue();
  virtual void setApplyCount(int count) {applyCount = count;}
  virtual void setChangedState(bool value) {changed = value;}
  virtual void setCopiedState(bool value) {copied = value;}
  void setId(int p_id) {id = p_id;}
  virtual void setLastId(int lid) {}
  virtual void setName(char* param_name) = 0;
  void setName(char* param_name, char* default_name);
  virtual void setParentId(int pr_id) { parentId = pr_id;}
  virtual void setParentEmfTag(int tag) { parentEmfTag = tag;}
  virtual void setParentEmfType(objectType tp) { parentEmfType = tp; }
  virtual void setSubParentEmfTag(int tag) { subParentEmfTag = tag;}
  virtual void setSubParentEmfType(objectType tp) { subParentEmfType = tp;}
  virtual void setUpdateFlag(bool value) {updateFlag = value;}
  virtual void setValue(char* values);
  void updateApplyCount(int increment) {applyCount += increment;}
  virtual void updateParentId();
  virtual void updateParentInfo(int parent_id);
  virtual void updateTargetTags() {}
protected:
  void create_fields();
  //void create_fields(ParameterField**& result_fields, char* source_string);
  void delete_fields();
  void delete_previousFields();
  void extractFieldInfo(istrstream& strm, char* field_name_buffer, char** var_name_buffers,
                        char& data_type, short& dim1, short& dim2, short& nof_variables,
                        bool& is_inactive);
  void extractFieldNameAndType(istrstream& strm, char* field_name_buffer, char& data_type);
  void extractFieldSizeAndVars(istrstream& strm, char** var_name_buffers,
                                short& dim1, short& dim2, short& nof_variables);
  void extractInactiveMarker(istrstream& strm);
  bool fieldInstanceExists(ParameterField* param_field);
  virtual void setData(int id, char* values, char* name);
  virtual void setData(int id, int parent_id, char* values, char* name);

  static Model* model;
  static ecif_parameterType type;
  int applyCount;
  bool changed;
  bool copied;
  char* dataString;
  char* dataString_previous;
  ParameterField** fields;
  int id;
  char* name;
  short nofFields;
  short nofPreviousFields;
  int parentId;                  // Object id for the "creating parent" object.
  int parentEmfTag;                // Parent tag in the emf-file for the "creating parent" (target) object
  enum objectType parentEmfType;   // Parent tag in the emf-file : Body, Face, Edge, Vertex
  int subParentEmfTag;              // Sub parent tag (ex. body layer) in the emf-file for the "creating parent" (target) object
  enum objectType subParentEmfType; // Sub parent type in the emf-file : Body Layer
  ParameterField** previousFields;
  ecif_modelStatus status;
  bool updateFlag;
};



#endif
