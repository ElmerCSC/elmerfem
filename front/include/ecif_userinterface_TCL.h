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
Module:     ecif_userinterface_TCL.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   A class for TCL based userinterface.

************************************************************************/

#ifndef _ECIF_USERINTERFACE_TCL_
#define _ECIF_USERINTERFACE_TCL_

//***NOTE: Order is important!
// Otherwise conflicts between:
// Win32 and X
// Tcl and Stl (list variable argument in tcl.h)
#include "ecif_def.h"
#include <tcl.h>
#include <tk.h>
#include "ecif_def_stl.h"
//***END NOTE

class UserInterface_TCL;

typedef void (UserInterface_TCL::*WIA_CHR) (Tcl_Interp*, char*, int, char*, int, char**);
typedef void (UserInterface_TCL::*WIA_INT) (Tcl_Interp*, char*, int, char*, int, int*);
typedef void (UserInterface_TCL::*WIA_DBL) (Tcl_Interp*, char*, int, char*, int, double*);

struct RendererInfo;

class UserInterface_TCL : public UserInterface {
public:
  UserInterface_TCL(Hinst application, char* script);
  void acceptEmfParameters(char* msg, int nof_params, char* pdatas, bool* accept_flags);
  void checkMeshInfoTs(char* ts);
  void colorFileWasRead(char* filename);
  void configureButtons(char* buttons, int state);
  void configureButtonOption(char* button, char* option, char* value);
  void configureMenuButtons(char* menu, char* buttons, int state);
  void configureMenuButtonOption(char* menu, char* button, char* option, char* value);
  void errMsg(int err_level, char* str1, char* str2 = NULL, char* str3 = NULL, char* str4 = NULL);
  void fieldNameGuiToSif(const char* gui_name, char* sif_name_buffer);
  void fieldNameSifToGui(const char* sif_name, char* gui_name_buffer);
  void generateEvent();
  void getCurrentTimestamp(char* buffer);
  bool getEquationVarsVariable(const char* equation_name, char*& equation_vars_name);
  bool getIsSolverTargetField(const char* equation_name, const char* field_name);
  void getMatcSifDefinitions(int& nof_defs, char**& defs);
  void getMeshDirectoryInfo(char*& dir, char*& dir_abs);
  bool getMeshInputFileName(char*& mif_file_name);
  void getModelDirectoryInfo(char*& dir, char*& dir_abs);
  void getModelNameInfo(char*& model_name, char*& problem_name);
  bool getParameterFieldInfo(const char* parameter, const char* field, ParameterFieldInfo& finfo);
  bool getSolverKeywordTypeGiven(const char* parameter, const char* field);
  Renderer* getRenderer();
  bool getUseModelFileSettings();
  bool getUseVariableNameInEquationName(const char* equation_name);
  void markSelectedBoundaries();
  void matcFileWasRead(char* filename);
  void saveModelPropertyData(Model* model);
  void selectBody(int bd1_id, int lr1_id, int bd2_id, int lr2_id, bool is_selected = true);
  void selectBoundary(int elem_id, int bd1_id, int lr1_id, int bd2_id, int lr2_id, bool extend = false);
  int  sendCommandToGui(const char* cmd, const char* arg = NULL);
  void setBoundarySelectionMode(int elem_id, bool is_selected = true, bool do_update = true);
  void setCurrentMeshH(double mesh_h);
  void setExceptionThrown();
  void setInitialMeshH(double mesh_h);
  void setInitialState();
  void setNeedsUpdate(const char* target);
  void setModelHasElmerMesh();
  void setModelHasMeshParameter();
  void setModelHasMatcDefinitions();
  void setMeshEdited();
  void setMeshExists();
  void setMeshInputUnit(double unit);
  void setParameterFieldValueState(int parameter_id, const char* field_name,
                                    bool has_value, bool value_has_changed);
  void setTimestamp(ecif_parameterType parameter, char* ts);
  void setWindowTitle(char* title);
  void setWasUpdated(const char* target);
  int showMsg(char* messge, short extra_line_feeds = 0, bool append = true);
  void showProgressMsg(Timer& timer, int frequency_nbr,
                       int nbr, int total_nbr,
                       char* text1 = NULL, char* text2 = NULL);
  void showUsedTimeMsg(double time, char* text,
                       short extra_line_feeds = 0, bool append = true);
  void showUsedTimeMsg(double time, char* text1,  int nof_objects, char* text2,
                       short extra_line_feeds = 0, bool append = true);
  void start(int argc, char** argv);
  void update();
  void update(int counter, int update_interval);
  void updateBodyData(Model* model);
  void updateBoundaryData(Model* model);
  void updateMeshZeroVelocityElements(int nof_zv_elements);
  void updateModelData(Model* model);
  void updateModelFlags(Model* model);
  void updateModelStatistics(Model* model);
  void updateModelStatus(Model* model);
  void updateNextActiveSelectionTolerance(double tolerance);
  void updateObjectData(Model* model);
  void updateParameterDataPre(Model* model);
  void updateParameterDataPost(Model* model);
  void updateRendererInfo(const RendererInfo& renderer_info);
  void variableNameGuiToSif(const char* gui_name, char* sif_name_buffer);
  void variableNameSifToGui(const char* sif_name, char* gui_name_buffer);

protected:
  static Tcl_Channel createFileChannel(Tcl_Interp* interp, Hfile native_handle);
  static void createTclCommands( Tcl_Interp* interp);
  static void createTclCommand(Tcl_Interp* interp, char* tcl_command, Tcl_CmdProc* cpp_func);
  Tcl_Interp* createTclEnvironment(Hinst application);
  static int from_tk_AcceptEmfParameters(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_BodySelected(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_BoundarySelected(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_BoundariesSelected(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_CheckMeshCornerElements(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_CheckModelStatus(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ColorHex2Name(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_CombineBoundaries(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_CopyParameters(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_CorrectMeshZeroVelocityElements(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_DoBreak(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_DoMatc(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_Exit(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_LoadMesh(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_OpenCadFile(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_OpenMeshFile(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_OpenModelFile(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ProcessExists(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ProcessResume(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ProcessSetPriorityLevel(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ProcessStart(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ProcessStop(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ProcessSuspend(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_PutModelStatusMessage(ClientData, Tcl_Interp*, int, char**);
  // Read parameters
  static int from_tk_ReadBodyForceData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadBodyParameterData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadBoundariesData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadBoundaryConditionData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadBoundaryConditionsForEdges(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadBoundaryConditionsForFaces(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadBoundaryConditionsForVertices(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadBoundaryParameterData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadCalculatorData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadColorFile(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadConstantData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadCoordinateData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadDatafileData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadEquationData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadEquationVariablesData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadGridHData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadGridParameterData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadInitialConditionData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadMatcFile(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadMaterialData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadMeshDefineData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadModelParameterData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadSimulationParameterData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadSolverData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadSolverControlData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadTimestepData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadUserSettingsData(ClientData, Tcl_Interp*, int, char**);
  // Read converteed older version data
  static int from_tk_ReadConvertedEquationData(ClientData, Tcl_Interp*, int, char**);
  // Read other
  static int from_tk_ReadActiveMeshIndices(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadBodyColors(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadBodyDeleteData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadBodyData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadBodyDisplayData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadBoundaryDisplayData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadBoundaryNames(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadBodyNames(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadDeletedParamIds(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadModelFileCreated(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadModelFileModified(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadModelFileTime(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadModelHasUserDefinitions(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadModelPropertyData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadProcessorData(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadSelectionTolerances(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadTimestamp(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ReadUserSettingFiles(ClientData, Tcl_Interp* interp, int, char**);
  static int from_tk_ReadVertexDisplayData(ClientData, Tcl_Interp*, int, char**);
  //
  static int from_tk_RemoveCadGeometry(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_RendererDisplayModel(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_RendererResetModel(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_RendererRotateModel(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_RendererScaleModel(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_RendererSetEditBoundaries(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_RendererSetRotatePriorities(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_RendererTranslateModel(ClientData, Tcl_Interp*, int, char**);

  static int from_tk_ResetBoundarySelections(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_RestoreBoundaryNames(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_ResetAllBoundarySelections(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SaveElmerMeshFile(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SaveElmerPostMeshFile(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SaveModelFile(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SaveMeshInputFile(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SaveSolverInputFile(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SaveThetisMeshFile(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SaveUserSettingsFile(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SelectMeshBoundaryElements(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SetCurrentDirectory(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SetFlagValue(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SetMatcInputFileEmf(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SetMatcInputFileSif(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SetMeshInputUnit(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SetModelStatus(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SetSelectionsToGui(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SplitBoundary(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SplitCombineBoundariesRedo(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_SplitCombineBoundariesUndo(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_StopEditMeshBoundaries(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_StoreBoundaryNames(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_UnloadMesh(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_UpdateCadGeometry(ClientData, Tcl_Interp*, int, char**);
  static int from_tk_UpdateMatcFile(ClientData, Tcl_Interp*, int, char**);
  static char* getCommandArguments(Tcl_Interp* interp);
  static char* getCommandResults(Tcl_Interp* interp);
  static void initTclVariables(Tcl_Interp* interp, const Model& model);
  static void listInnerIds(Model* model,
              ostrstream& body_ids, ostrstream& elem_ids,
              ostrstream& names, ostrstream& bndr_ids,
              const char obj_sep, const char fld_sep,
              Ids3Set& id_set);
  static void listOuterIds(Model* model,
              ostrstream& body_ids, ostrstream& elem_ids,
              ostrstream& names, ostrstream& bndr_ids,
              const char obj_sep, const char fld_sep,
              Ids2Set& id_set);
  void readUserSettingFiles(Tcl_Interp* interp);
  void readUserSettingsFile(Tcl_Interp* interp, char* filename);
  static int readUserSettingsFileCallBack(void** user_data);

  // Array, single value
  static void readVariable(Tcl_Interp* interp, const char* array, const char* variable,
                            char& value);
  static void readVariable(Tcl_Interp* interp, const char* array, const char* variable,
                            char*& value);
  static void readVariable(Tcl_Interp* interp, const char* array, const char* variable,
                                  int& value);
  static void readVariable(Tcl_Interp* interp, const char* array, const char* variable,
                                  double& value);
  // Array, list of values
  static void readVariable(Tcl_Interp* interp, const char* array, const char* variable,
                                  int& size, char**& values);
  static void readVariable(Tcl_Interp* interp, const char* array, const char* variable,
                                  int& size, int*& values);
  static void readVariable(Tcl_Interp* interp, const char* array, const char* variable,
                                  int& size, double*& values);

  // Id Array, integer id, single value
  static void readIdVariable(Tcl_Interp* interp, const char* array, int id,
                          const char* variable, char& value);
  static void readIdVariable(Tcl_Interp* interp, const char* array, int id,
                          const char* variable, char*& value);
  static void readIdVariable(Tcl_Interp* interp, const char* array, int id,
                          const char* variable, int& value);
  static void readIdVariable(Tcl_Interp* interp, const char* array, int id,
                           const char* variable, double& value);
  // Id array, integer id, list of values
  static void readIdVariable(Tcl_Interp* interp, const char* array, int id,
                           const char* variable, int& size, char**& values);
  static void readIdVariable(Tcl_Interp* interp, const char* array, int id,
                           const char* variable, int& size, int*& values);
  static void readIdVariable(Tcl_Interp* interp, const char* array, int id,
                           const char* variable, int& size, double*& values);

  // Id Array, string id, single value
  static void readIdVariable(Tcl_Interp* interp, const char* array, const char* id,
                          const char* variable, char& value);
  static void readIdVariable(Tcl_Interp* interp, const char* array, const char* id,
                          const char* variable, char*& value);
  static void readIdVariable(Tcl_Interp* interp, const char* array, const char* id,
                          const char* variable, int& value);
  static void readIdVariable(Tcl_Interp* interp, const char* array, const char* id,
                           const char* variable, double& value);
  // Id array, string id, list of values
  static void readIdVariable(Tcl_Interp* interp, const char* array, const char* id,
                           const char* variable, int& size, char**& values);
  static void readIdVariable(Tcl_Interp* interp, const char* array, const char* id,
                           const char* variable, int& size, int*& values);
  static void readIdVariable(Tcl_Interp* interp, const char* array, const char* id,
                           const char* variable, int& size, double*& values);


  // Id2 Array, integer id1,id2, single value
  static void readId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2,
                          const char* variable, char& value);
  static void readId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2,
                          const char* variable, char*& value);
  static void readId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2,
                          const char* variable, int& value);
  static void readId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2,
                           const char* variable, double& value);
  // Id2 array, integer id1,id2, list of values
  static void readId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2,
                           const char* variable, int& size, char**& values);
  static void readId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2,
                           const char* variable, int& size, int*& values);
  static void readId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2,
                           const char* variable, int& size, double*& values);

  // Template based implementations
  template <class T >
  static void readVariable_impl_1(Tcl_Interp* interp,
                                   const char* array, const char* variable,
                                   T& value);
  template <class T >
  static void readVariable_impl_n(Tcl_Interp* interp,
                                   const char* array, const char* variable,
                                   int& size, T*& values);

  static int sendCommandToGui(Tcl_Interp* interp, const char* cmd, const char* arg = NULL);
  static void setBoundaryConditions(Tcl_Interp* interp, Model* model, MultiIdTable& bc_table);
  static void setBoundaryConditionsForEdges(Tcl_Interp* interp, Model* model, MultiIdTable& bc_table);
  static void setBoundaryConditionsForFaces(Tcl_Interp* interp, Model* model, MultiIdTable& bc_table);
  static void setBoundaryConditionsForVertices(Tcl_Interp* interp, Model* model, MultiIdTable& bc_table);
  static int setParameterData(Model* model, ecif_parameterType param_type,
                              Tcl_Interp* interp, const char* array_name);
  static bool start_Tcl_MainLoop();
  static void to_tk_WriteAllParamDataPre(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteAllParamDataPost(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteBodyData(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteBodyLayerData(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteBodyInfoData(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteBodyMeshInfoData(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteBoundaryData(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteBoundaryElements(Tcl_Interp* interp, Model* model, objectType btype);
  static void to_tk_WriteElementGroupData(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteControlParameters(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteModelData(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteModelFlags(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteModelGeometryDimension(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteModelStatus(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteModelStatusMessage(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteStats(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteStatusBodyForces(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteStatusBoundaryConditions(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteStatusEquations(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteStatusInitialConditions(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteStatusMaterials(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteStatusMeshes(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteStatusTimesteps(Tcl_Interp* interp, const Model& model);
  static void to_tk_WriteProcessorData(Tcl_Interp* interp, const Model& model);

  // Array, single value
  static void writeVariable(Tcl_Interp* interp, const char* array, const char* variable,
                            const char* value, bool reset = true);
  static void writeVariable(Tcl_Interp* interp, const char* array, const char* variable,
                            int value, bool reset = true);
  static void writeVariable(Tcl_Interp* interp, const char* array, const char* variable,
                            double value, bool reset = true);
  // Array, list of values
  static void writeVariable(Tcl_Interp* interp, const char* array, const char* variable,
                            int size, const char** values, bool reset = true);
  static void writeVariable(Tcl_Interp* interp, const char* array, const char* variable,
                            int size, const int* values, bool reset = true);
  static void writeVariable(Tcl_Interp* interp, const char* array, const char* variable,
                            int size, const double* values, bool reset = true);

  // Id Array, integer id, single value
  static void writeIdVariable(Tcl_Interp* interp, const char* array, int id, const char* variable,
                              const char* value, bool reset = true);
  static void writeIdVariable(Tcl_Interp* interp, const char* array, int id, const char* variable,
                              int value, bool reset = true);
  static void writeIdVariable(Tcl_Interp* interp, const char* array, int id, const char* variable,
                              ProcessId value, bool reset = true);
  static void writeIdVariable(Tcl_Interp* interp, const char* array, int id, const char* variable,
                              double value, bool reset = true);
  // Id array, integer id, list of values
  static void writeIdVariable(Tcl_Interp* interp, const char* array, int id,
                           const char* variable, int size, const char** values,  bool reset = true);
  static void writeIdVariable(Tcl_Interp* interp, const char* array, int id,
                           const char* variable, int size, const int* values, bool reset = true);
  static void writeIdVariable(Tcl_Interp* interp, const char* array, int id,
                           const char* variable, int size, const double* values, bool reset = true);

    // Id Array, string id, single value
  static void writeIdVariable(Tcl_Interp* interp, const char* array, const char* id, const char* variable,
                              const char* value, bool reset = true);
  static void writeIdVariable(Tcl_Interp* interp, const char* array, const char* id, const char* variable,
                              int value, bool reset = true);
  static void writeIdVariable(Tcl_Interp* interp, const char* array, const char* id, const char* variable,
                              ProcessId value, bool reset = true);
  static void writeIdVariable(Tcl_Interp* interp, const char* array, const char* id, const char* variable,
                              double value, bool reset = true);
  // Id array, string id, list of values
  static void writeIdVariable(Tcl_Interp* interp, const char* array, const char* id,
                           const char* variable, int size, const char** values,  bool reset = true);
  static void writeIdVariable(Tcl_Interp* interp, const char* array, const char* id,
                           const char* variable, int size, const int* values, bool reset = true);
  static void writeIdVariable(Tcl_Interp* interp, const char* array, const char* id,
                           const char* variable, int size, const double* values, bool reset = true);

  // Id2 Array, integer id1, id2, single value
  static void writeId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2, const char* variable,
                              const char* value, bool reset = true);
  static void writeId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2, const char* variable,
                              int value, bool reset = true);
  static void writeId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2, const char* variable,
                              ProcessId value, bool reset = true);
  static void writeId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2, const char* variable,
                              double value, bool reset = true);
  // Id2 array, integer id1, id2, list of values
  static void writeId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2,
                           const char* variable, int size, const char** values,  bool reset = true);
  static void writeId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2,
                           const char* variable, int size, const int* values, bool reset = true);
  static void writeId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2,
                           const char* variable, int size, const double* values, bool reset = true);


  // Template based implementations
  template <class T >
  static void writeVariable_impl_1(Tcl_Interp* interp,
                                   const char* array, const char* variable,
                                   const T value,
                                   bool reset = true);
  template <class T >
  static void writeVariable_impl_n(Tcl_Interp* interp,
                                   const char* array, const char* variable,
                                   int size, const T* values,
                                   bool reset = true);

  static void getTclObjValue(Tcl_Interp* interp, Tcl_Obj* value_obj, bool& value);
  static void getTclObjValue(Tcl_Interp* interp, Tcl_Obj* value_obj, int& value);
  static void getTclObjValue(Tcl_Interp* interp, Tcl_Obj* value_obj, long& value);
  static void getTclObjValue(Tcl_Interp* interp, Tcl_Obj* value_obj, double& value);
  static void getTclObjValue(Tcl_Interp* interp, Tcl_Obj* value_obj, char& value);
  static void getTclObjValue(Tcl_Interp* interp, Tcl_Obj* value_obj, char*& value);

  static void setTclObjValue(Tcl_Obj*& value_obj, const bool value);
  static void setTclObjValue(Tcl_Obj*& value_obj, const int value);
  static void setTclObjValue(Tcl_Obj*& value_obj, const ProcessId value);
  static void setTclObjValue(Tcl_Obj*& value_obj, const long value);
  static void setTclObjValue(Tcl_Obj*& value_obj, const double value);
  static void setTclObjValue(Tcl_Obj*& value_obj, const char value);
  static void setTclObjValue(Tcl_Obj*& value_obj, const char* value);

  static int unknownFieldMsg(emf_ObjectData_X* object_data, bool is_fatal);
  static void update(Tcl_Interp* interp);

  // Class attributes
  static char* controlSideScript;
  static char* tclScriptPath;
  static Tcl_Interp* theInterp;
};

#endif
