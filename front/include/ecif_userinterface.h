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
Module:     ecif_userinterface.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   An absract base class for userinterface.

************************************************************************/

#ifndef _ECIF_USERINTERFACE_
#define _ECIF_USERINTERFACE_

#include "ecif_def.h"


class UserInterface {
public:
  friend class Control;
  UserInterface() {applicationName = NULL; theModel = NULL;}
  UserInterface(Hinst a_n, Model* m) {applicationName = a_n; theModel = m;}
  virtual void acceptEmfParameters(char* msg, int nof_params, char* pdatas, bool* accept_flags) {}
  virtual void checkMeshInfoTs(char* ts) {}
  virtual void colorFileWasRead(char* filename) {}
  virtual void configureButtons(char* buttons, int state) {}
  virtual void configureButtonOption(char* button, char* option, char* value) {}
  virtual void configureMenuButtons(char* menu, char* buttons, int state) {}
  virtual void configureMenuButtonOption(char* menu, char* button, char* option, char* value) {}
  virtual void errMsg(int err_level, char* str1, char* str2 = NULL, char* str3 = NULL, char* str4 = NULL);
  virtual void fieldNameGuiToSif(const char* gui_name, char* sif_name_buffer) {}
  virtual void fieldNameSifToGui(const char* sif_name, char* gui_name_buffer) {}
  virtual void generateEvent() {}
  virtual void getCurrentTimestamp(char* buffer) {};
  virtual bool getEquationVarsVariable(const char* equation_name, char*& equation_vars_name) { return false;}
  virtual bool getIsSolverTargetField(const char* equation_name, const char* field_name) { return false; }
  virtual void getMeshDirectoryInfo(char*& dir, char*& dir_abs) {};
  virtual void getMatcSifDefinitions(int& nof_defs, char**& defs) {};
  virtual bool getMeshInputFileName(char*& mif_file_name) {return false;}
  virtual void getModelDirectoryInfo(char*& dir, char*& dir_abs) {};
  virtual void getModelNameInfo(char*& model_name, char*& problem_name) {};
  virtual bool getParameterFieldInfo(const char* parameter, const char* field, ParameterFieldInfo& finfo) { return false; }
  virtual bool getSolverKeywordTypeGiven(const char* parameter, const char* field) { return false;}
  virtual Renderer* getRenderer() { return NULL; }
  virtual bool getUseModelFileSettings() {return true;}
  virtual bool getUseVariableNameInEquationName(const char* equation_name) { return false;}
  virtual void markSelectedBoundaries() {}
  virtual void matcFileWasRead(char* filename);
  virtual void pause(double seconds = 5.0, bool show=false);
  virtual void saveModelPropertyData(Model* model) {}
  virtual void selectBody(int bd1_id, int lr1_id, int bd2_id, int lr2_id, bool is_selected = true) {}
  virtual void selectBoundary(int elem_id, int bd1_id, int lr1_id, int bd2_id, int lr2_id, bool extend = false) {}
  virtual void selectBoundaries(int nof_elems, int* elem_ids) {}
  virtual int sendCommandToGui(const char* cmd, const char* arg = NULL) {return 0;};
  virtual void setBoundarySelectionMode(int elem_id, bool is_selected = true, bool do_update = true) {}
  virtual void setCurrentMeshH(double mesh_h) {}
  virtual void setExceptionThrown() {}
  virtual void setInitialMeshH(double mesh_h) {}
  virtual void setInitialState() {}
  virtual void setNeedsUpdate(const char* target) {}
  virtual void setModelHasElmerMesh() {}
  virtual void setModelHasMeshParameter() {}
  virtual void setModelHasMatcDefinitions() {}
  virtual void setMeshEdited() {}
  virtual void setMeshExists() {}
  virtual void setMeshInputUnit(double unit) {}
  virtual void setParameterFieldValueState(int parameter_id, const char* field_name,
                                           bool has_value, bool value_has_changed) {}
  virtual void setTimestamp(ecif_parameterType parameter, char* ts) {}
  virtual void setWindowTitle(char* title) {}
  virtual void setWasUpdated(const char* target) {}
  virtual int showMsg(char* messge, short extra_line_feeds = 0, bool append = true);
  virtual void showProgressMsg(Timer& timer, int frequency_nbr,
                                int nbr, int total_nbr,
                                char* text1 = NULL, char* text2 = NULL) {}
  virtual void showUsedTimeMsg(double time, char* text,
                                short extra_line_feeds = 0, bool append = true) {}
  virtual void showUsedTimeMsg(double time, char* text1,int nof_objects,  char* text2,
                                short extra_line_feeds = 0, bool append = true) {}
  virtual void start(int argc, char** argv);
  virtual void update() {}
  virtual void update(int counter, int update_interval) {}
  virtual void updateBodyData(Model* model) {}
  virtual void updateBoundaryData(Model* model) {}
  virtual void updateMeshZeroVelocityElements(int nof_zv_elements) {}
  virtual void updateModelData(Model* model) {}
  virtual void updateModelFlags(Model* model) {}
  virtual void updateModelStatistics(Model* model) {}
  virtual void updateModelStatus(Model* model) {}
  virtual void updateNextActiveSelectionTolerance(double tolerance) {}
  virtual void updateObjectData(Model* model) {}
  virtual void updateParameterDataPre(Model* model) {}
  virtual void updateParameterDataPost(Model* model) {}
  virtual void updateRendererInfo(const RendererInfo& renderer_info) {}
  virtual void variableNameGuiToSif(const char* gui_name, char* sif_name_buffer) {}
  virtual void variableNameSifToGui(const char* sif_name, char* gui_name_buffer) {}
protected:
  Hinst applicationName;
  static Control* theControlCenter;
  Model* theModel;
};


#endif
