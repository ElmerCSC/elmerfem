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
Module:     ecif_userinterface_TCL.cpp
Language:   C++
Date:       21.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "../config.h"

#if defined(WIN32)
  #include <direct.h>
#else
  #include <unistd.h>
#endif

// NOTE: Include order is imporant here!
// Otherwise conflicts with Windows and Stl
#include "ecif_userinterface.h"
#include "ecif_userinterface_TCL.h"
// END NOTE

#include "ecif_body.h"
#include "ecif_bodyLayer.h"
#include "ecif_bodyElement.h"
#include "ecif_bodyElementGroup.h"
#include "ecif_bodyForce.h"
#include "ecif_boundaryCondition.h"
#include "ecif_constant.h"
#include "ecif_coordinate.h"
#include "ecif_control.h"
#include "ecif_datafile.h"
#include "ecif_equation.h"
#include "ecif_geometry.h"
#include "ecif_gridParameter.h"
#include "ecif_initialCondition.h"
#include "ecif_material.h"
#include "ecif_model.h"
#include "ecif_modelObject.h"
#include "ecif_process.h"
#include "ecif_renderer.h"
#include "ecif_renderer_OGL.h"
#include "ecif_solver.h"
#include "ecif_timer.h"
#include "ecif_timestep.h"

#ifdef UNIX
#include "ecif_renderer.h"
#endif


void Tcl_MainLoop();
int My_Tcl_AppInit(Tcl_Interp* interp);
void WishPanic _ANSI_ARGS_(TCL_VARARGS(char *,format));
int TkTest_Init(Tcl_Interp *interp);

void tcl_DisplayIdleProc(ClientData data);
void tcl_InterruptIdleProc(ClientData data);
void tcl_interrupt(ClientData data);

extern void display_system_msg(char* header);
extern void display_msg(char* header, char* msg);

extern UserInterface* TheUI = NULL;

ClientData displayIdleData;
ClientData interruptIdleData;

const int glob_flag = TCL_GLOBAL_ONLY | TCL_LEAVE_ERR_MSG;
const int int_ro_flag = TCL_LINK_INT | TCL_LINK_READ_ONLY;
const int str_ro_flag =   TCL_LINK_STRING | TCL_LINK_READ_ONLY;

const char tclObjectSeparator = '!';
const char tclFieldSeparator =  '&';
const char tclCmdSeparator =    '@';
const char tclArgSeparator =    '^';

const char subIdSep = '.';
const char sis = subIdSep;

const int isOK = 0;



// *****************************************

// =========================================
//          UserInterface base class
// =========================================

void
UserInterface::errMsg(int err_level, char* str1, char* str2, char* str3, char* str4)
{
  strstream strm;

  if ( str1 != NULL ) strm << str1;
  if ( str2 != NULL ) strm << str2;
  if ( str3 != NULL ) strm << str3;
  if ( str4 != NULL ) strm << str4;
  strm << ends;

  if ( strm.str() != NULL ) cerr << strm.str();

}


void
UserInterface::matcFileWasRead(char* filename)
{
  strstream strm;
  strm << "***MATC file loaded: " << filename << ends;

  showMsg(strm.str());
}


// Simulate pause (in seconds)
void
UserInterface::pause(double seconds, bool show)
{
  Timer timer;
  strstream strm;

  if ( show ) {
    strm << "---Pause " << seconds << " seconds..." << ends;
    showMsg(strm.str());
  } else {
    generateEvent();
    update();
  }

  timer.start();

  while ( timer.getLapTime(WALL_TIME) < seconds );

  if ( show ) {
    showMsg("---Pause done!...");
  } else {
    generateEvent();
    update();
  }
}



int
UserInterface::showMsg(char* msg, short extra_line_feeds, bool append)
{
  cerr << msg << endl;

  return 0;
}


// ----------------
// Batch-mode start
// ----------------
//
void
UserInterface::start(int argc, char** argv)
{
  char* model_dir = NULL;
  char* model_name = NULL;
  char* mesh_name = NULL;
  char* mesh_log = NULL;

  for ( int i = 0; i < argc; i++) {

    if ( LibFront::ncEqualPartial(argv[i], "model-directory=") ) {
      char* tmp = LibFront::trimLeft(argv[i]);
      tmp +=16;
      update_dyna_string(model_dir, tmp);

    } else if ( LibFront::ncEqualPartial(argv[i], "model-name=") ) {
      char* tmp = LibFront::trimLeft(argv[i]);
      tmp +=11;
      update_dyna_string(model_name, tmp);

    } else if ( LibFront::ncEqualPartial(argv[i], "mesh-name=") ) {
      char* tmp = LibFront::trimLeft(argv[i]);
      tmp +=10;
      update_dyna_string(mesh_name, tmp);

    } else if ( LibFront::ncEqualPartial(argv[i], "mesh-log=") ) {
      char* tmp = LibFront::trimLeft(argv[i]);
      tmp +=9;
      update_dyna_string(mesh_log, tmp);
    }
  }

  if ( model_name == NULL ||
       mesh_name == NULL
     ) {
    cerr << "Usage: ElmerFront --batch=1 "
         << "[--model-directory=dirname (default: cwd)] "
         << "<--model-name=name> "
         << "<--mesh-name=dirname> "
         << "[--mesh-log=filename (default: cerr)]"
         << endl;
      delete[] model_dir;
      delete[] model_name;
      delete[] mesh_name;
      delete[] mesh_log;
      return;
  }

  if ( model_dir == NULL ) {
    update_dyna_string(model_dir, "./");
  }

  strstream strm, strm1, strm2, strm3, strm4;

  // Check that model file exists
  //
  strm << model_dir << "/" << model_name << ".emf" << ends;

  // NOTE: Possible error message comes from control, no need to output here!
  cerr << endl;
  if ( !theControlCenter->readModelFile(strm.str(), false, true) ) {
    cerr << endl;
    return;
  }

  // Check that mesh directory exists (check by parts!)
  //
  strm1 << model_dir << "/MESHDIR" << ends;  // Model's MESHDIR subdirectory
  if ( !front_mkdir(strm1.str()) ) {
    cerr << "***ERROR: Model's MESHDIR does not exist: " << strm1.str()
         << endl;
   return;
  }

  strm2 << strm1.str() << '/' << mesh_name << ends;  // Mesh-dir

  if ( !front_mkdir(strm2.str()) ) {
    cerr << "***ERROR: Mesh directory does not exist: " << strm2.str()
         << endl;
     return;
  }

  // Create and save mif-file
  //
  strm3 << strm2.str() << '/' << model_name << ".mif" << ends; // Mif-filename
  theControlCenter->saveMeshInputFile(strm3.str());

  char curr_wd[1025];

  front_getcwd(curr_wd, 1024);

  // Go to mesh directory to start ElmerMesh2D
  //
  if ( !front_chdir(strm2.str()) ) {
    cerr << "**ERROR: Cannot change to the mesh directory: " << strm2.str()
         << endl;
     return;
    return;
  }

  // Call: ElmerMesh2D mif-filename
  //
  //strm4 << strm3.str() << ends;
  strm4 << model_name << ".mif" << ends;
  Process* process = new Process("ElmerMesh2D", strm4.str());

  if ( mesh_log != NULL ) {
    process->setLogfile(mesh_log);
  } else {
    process->setShowConsole(true);
  }

  cerr << endl << "Calling: ElmerMesh2D " << strm4.str()
       << endl;

  process->start();

  // Back to current directory
  //
  if ( !front_chdir(curr_wd) ) return;
}


// =========================================
//          UserInterface_TCL class
// =========================================

// Initialize class variables
Control* UserInterface::theControlCenter = NULL;
Tcl_Interp* UserInterface_TCL::theInterp = NULL;
char* UserInterface_TCL::controlSideScript = NULL;
char* UserInterface_TCL::tclScriptPath = NULL;


// Constructor
// ===========
UserInterface_TCL::UserInterface_TCL(Hinst application, char* ctrl_script):
  UserInterface(application, NULL)
{
  controlSideScript = ctrl_script;
  tclScriptPath = NULL;

  theInterp = createTclEnvironment(application);

  if ( theInterp == NULL) {
    exit(1);
  }

  // Hm... this is ugly, but we need
  // a global UI object to give any message
  // for interrupt handling
  TheUI = this;

  // *** Add new commands to the Tcl-enviromnment
  createTclCommands(theInterp);
}


// Open a dialog in gui, where the user can accept list of parmaters
// ( ex. to be copied into model)
void
UserInterface_TCL::acceptEmfParameters(char* msg, int nof_params, char* pdatas, bool* accept_flags)
{
  int i;

  for (i = 0; i < nof_params; i++) {
    accept_flags[i] = false;
  }

  ostrstream strm;
  strm << msg << tclArgSeparator
       << pdatas << ends;

  int rc = sendCommandToGui(theInterp, "Interface::acceptParameters", strm.str());

  char* result = NULL;

  Tcl_CreateTimerHandler(1, tcl_DisplayIdleProc, displayIdleData);

  int event_types = TCL_ALL_EVENTS;

  while ( Tcl_DoOneEvent(event_types) ) {

    result = getCommandResults(theInterp);

    if ( result != NULL && result[0] != '\0' ) break;

    generateEvent();
  }

  if ( result == NULL || result[0] == '\0' ) return;

  istrstream flags(result);
  int flag;

  for (i = 0; i < nof_params; i++) {
    flags >> flag;
    accept_flags[i] = bool(flag);
  }
}


void
UserInterface_TCL::checkMeshInfoTs(char* ts) {
  sendCommandToGui(theInterp,"Interface::checkMeshInfoTs", ts);
}


void
UserInterface_TCL::colorFileWasRead(char* filename)
{
  sendCommandToGui(theInterp, "Interface::colorFileWasRead", filename);
}


// Menu state setting commands
// ===========================

// One option for one button
void
UserInterface_TCL::configureButtonOption(char* button, char* option, char* value)
{
  ostrstream strm;
  strm << button << tclArgSeparator
       << option << tclArgSeparator
       << value << ends;

  sendCommandToGui(theInterp, "Interface::configureButtonOption", strm.str());
}


// State for buttons
void
UserInterface_TCL::configureButtons(char* buttons, int state)
{
  ostrstream strm;
  strm << buttons << tclArgSeparator
       << state << ends;

  sendCommandToGui(theInterp, "Interface::configureButtons", strm.str());
}


// One option for one menu button
void
UserInterface_TCL::configureMenuButtonOption(char* menu, char* button, char* option, char* value)
{
  ostrstream strm;
  strm << menu << tclArgSeparator
       << button << tclArgSeparator
       << option << tclArgSeparator
       << value << ends;

  sendCommandToGui(theInterp, "Interface::configureMenuButtonOption", strm.str());
}


// State for menu buttons
void
UserInterface_TCL::configureMenuButtons(char* menu, char* buttons, int state)
{
  ostrstream strm;
  strm << menu << tclArgSeparator
       << buttons << tclArgSeparator
       << state << ends;

  sendCommandToGui(theInterp, "Interface::configureMenuButtons", strm.str());
}


Tcl_Channel
UserInterface_TCL::createFileChannel(Tcl_Interp* interp, Hfile native_handle)
{
    Tcl_Channel channel = Tcl_MakeFileChannel((ClientData)native_handle, TCL_READABLE);
    Tcl_SetChannelOption(interp, channel, "-blocking", "0");

    Tcl_RegisterChannel(interp, channel);

    return channel;
}


// REGISTRATION of new Tcl commands
// ================================

// Cover function to make writing interface a bit shorter.
void
UserInterface_TCL::createTclCommand(Tcl_Interp* interp, char* tcl_command, Tcl_CmdProc* cpp_func)
{
   Tcl_CreateCommand(interp, tcl_command, cpp_func, (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL);
}

#define TCL_FUNC_PROTO (int (*)(void*,Tcl_Interp*,int,const char**))
// For each command, there is also a cpp-function defined!!!
void
UserInterface_TCL::createTclCommands( Tcl_Interp* interp)
{
  createTclCommand(interp, "cpp_doMatc",  TCL_FUNC_PROTO from_tk_DoMatc);

  // ------ File menu commands ------
  createTclCommand(interp, "cpp_openCadFile", TCL_FUNC_PROTO from_tk_OpenCadFile);
  createTclCommand(interp, "cpp_openModelFile", TCL_FUNC_PROTO from_tk_OpenModelFile);
  createTclCommand(interp, "cpp_openMeshFile", TCL_FUNC_PROTO from_tk_OpenMeshFile);
  createTclCommand(interp, "cpp_loadMesh", TCL_FUNC_PROTO from_tk_LoadMesh);
  createTclCommand(interp, "cpp_unloadMesh", TCL_FUNC_PROTO from_tk_UnloadMesh);
  createTclCommand(interp, "cpp_saveModelFile", TCL_FUNC_PROTO from_tk_SaveModelFile);
  createTclCommand(interp, "cpp_saveMeshInputFile", TCL_FUNC_PROTO from_tk_SaveMeshInputFile);
  createTclCommand(interp, "cpp_saveSolverInputFile", TCL_FUNC_PROTO from_tk_SaveSolverInputFile);
  createTclCommand(interp, "cpp_saveElmerMeshFile", TCL_FUNC_PROTO from_tk_SaveElmerMeshFile);
  createTclCommand(interp, "cpp_saveElmerPostMeshFile", TCL_FUNC_PROTO from_tk_SaveElmerPostMeshFile);
  createTclCommand(interp, "cpp_saveThetisMeshFile", TCL_FUNC_PROTO from_tk_SaveThetisMeshFile);
  createTclCommand(interp, "cpp_copyParameters", TCL_FUNC_PROTO from_tk_CopyParameters);
  // Program exit.
  createTclCommand(interp, "cpp_exit", TCL_FUNC_PROTO from_tk_Exit);

  //------ Edit menu commands -------
  createTclCommand(interp, "cpp_removeCadGeometry", TCL_FUNC_PROTO from_tk_RemoveCadGeometry);
  createTclCommand(interp, "cpp_setMatcInputFileEmf", TCL_FUNC_PROTO from_tk_SetMatcInputFileEmf);
  createTclCommand(interp, "cpp_setMatcInputFileSif", TCL_FUNC_PROTO from_tk_SetMatcInputFileSif);
  createTclCommand(interp, "cpp_saveUserSettingsFile", TCL_FUNC_PROTO from_tk_SaveUserSettingsFile);
  createTclCommand(interp, "cpp_updateCadGeometry", TCL_FUNC_PROTO from_tk_UpdateCadGeometry);
  createTclCommand(interp, "cpp_setMeshInputUnit", TCL_FUNC_PROTO from_tk_SetMeshInputUnit);

  //------ Display menu commands ------
  createTclCommand(interp, "cpp_bodyDisplayPanelOk", TCL_FUNC_PROTO from_tk_ReadBodyDisplayData);
  createTclCommand(interp, "cpp_boundaryDisplayPanelOk", TCL_FUNC_PROTO from_tk_ReadBoundaryDisplayData);
  createTclCommand(interp, "cpp_rendererDisplayModel", TCL_FUNC_PROTO from_tk_RendererDisplayModel);
  createTclCommand(interp, "cpp_rendererResetModel", TCL_FUNC_PROTO from_tk_RendererResetModel);
  createTclCommand(interp, "cpp_rendererRotateModel", TCL_FUNC_PROTO from_tk_RendererRotateModel);
  createTclCommand(interp, "cpp_rendererScaleModel", TCL_FUNC_PROTO from_tk_RendererScaleModel);
  createTclCommand(interp, "cpp_rendererSetEditBoundaries", TCL_FUNC_PROTO from_tk_RendererSetEditBoundaries);
  createTclCommand(interp, "cpp_rendererSetRotatePriorities", TCL_FUNC_PROTO from_tk_RendererSetRotatePriorities);
  createTclCommand(interp, "cpp_rendererTranslateModel", TCL_FUNC_PROTO from_tk_RendererTranslateModel);
  createTclCommand(interp, "cpp_vertexDisplayPanelOk", TCL_FUNC_PROTO from_tk_ReadVertexDisplayData);

  //------ Data transfer and update from GUI ------
  // Panel Ok
  createTclCommand(interp, "cpp_bodyForcePanelOk", TCL_FUNC_PROTO from_tk_ReadBodyForceData);
  createTclCommand(interp, "cpp_bodyParameterPanelOk", TCL_FUNC_PROTO from_tk_ReadBodyParameterData);
  createTclCommand(interp, "cpp_bodyPropertyPanelOk", TCL_FUNC_PROTO from_tk_ReadBodyData);
  createTclCommand(interp, "cpp_bodyPropertyPanelColorsOk", TCL_FUNC_PROTO from_tk_ReadBodyColors);
  createTclCommand(interp, "cpp_bodyPropertyPanelNamesOk", TCL_FUNC_PROTO from_tk_ReadBodyNames);
  createTclCommand(interp, "cpp_bodyPropertyPanelDeleteBodyOk", TCL_FUNC_PROTO from_tk_ReadBodyDeleteData);
  createTclCommand(interp, "cpp_boundariesPanelOk", TCL_FUNC_PROTO from_tk_ReadBoundariesData);
  createTclCommand(interp, "cpp_boundaryConditionPanelOk", TCL_FUNC_PROTO from_tk_ReadBoundaryConditionData);
  createTclCommand(interp, "cpp_boundaryParameterPanelOk", TCL_FUNC_PROTO from_tk_ReadBoundaryParameterData);
  createTclCommand(interp, "cpp_calculatorPanelOk", TCL_FUNC_PROTO from_tk_ReadCalculatorData);
  createTclCommand(interp, "cpp_constantPanelOk", TCL_FUNC_PROTO from_tk_ReadConstantData);
  createTclCommand(interp, "cpp_coordinatePanelOk", TCL_FUNC_PROTO from_tk_ReadCoordinateData);
  createTclCommand(interp, "cpp_datafilePanelOk", TCL_FUNC_PROTO from_tk_ReadDatafileData);
  createTclCommand(interp, "cpp_equationPanelOk", TCL_FUNC_PROTO from_tk_ReadEquationData);
  createTclCommand(interp, "cpp_equationVariablesPanelOk", TCL_FUNC_PROTO from_tk_ReadEquationVariablesData);
  createTclCommand(interp, "cpp_gridHPanelOk", TCL_FUNC_PROTO from_tk_ReadGridHData);
  createTclCommand(interp, "cpp_gridParameterPanelOk", TCL_FUNC_PROTO from_tk_ReadGridParameterData);
  createTclCommand(interp, "cpp_initialConditionPanelOk", TCL_FUNC_PROTO from_tk_ReadInitialConditionData);
  createTclCommand(interp, "cpp_materialPanelOk", TCL_FUNC_PROTO from_tk_ReadMaterialData);
  createTclCommand(interp, "cpp_meshDefinePanelOk", TCL_FUNC_PROTO from_tk_ReadMeshDefineData);
  createTclCommand(interp, "cpp_modelParameterPanelOk", TCL_FUNC_PROTO from_tk_ReadModelParameterData);
  createTclCommand(interp, "cpp_modelPropertyPanelOk", TCL_FUNC_PROTO from_tk_ReadModelPropertyData);
  createTclCommand(interp, "cpp_processorPanelOk", TCL_FUNC_PROTO from_tk_ReadProcessorData);
  createTclCommand(interp, "cpp_simulationParameterPanelOk", TCL_FUNC_PROTO from_tk_ReadSimulationParameterData);
  createTclCommand(interp, "cpp_solverParameterPanelOk", TCL_FUNC_PROTO from_tk_ReadSolverData);
  createTclCommand(interp, "cpp_solverControlPanelOk", TCL_FUNC_PROTO from_tk_ReadSolverControlData);
  createTclCommand(interp, "cpp_timestepPanelOk", TCL_FUNC_PROTO from_tk_ReadTimestepData);
  createTclCommand(interp, "cpp_userSettingPanelOk", TCL_FUNC_PROTO from_tk_ReadUserSettingsData);

  // Older version data conversion save procs
  createTclCommand(interp, "cpp_readConvertedEquationData", TCL_FUNC_PROTO from_tk_ReadConvertedEquationData);

  // Flags and states
  createTclCommand(interp, "cpp_checkMeshCornerElements", TCL_FUNC_PROTO from_tk_CheckMeshCornerElements);
  createTclCommand(interp, "cpp_correctMeshZeroVelocityElements", TCL_FUNC_PROTO from_tk_CorrectMeshZeroVelocityElements);
  createTclCommand(interp, "cpp_checkModelStatus", TCL_FUNC_PROTO from_tk_CheckModelStatus);
  createTclCommand(interp, "cpp_putModelStatusMessage", TCL_FUNC_PROTO from_tk_PutModelStatusMessage);
  createTclCommand(interp, "cpp_setFlagValue", TCL_FUNC_PROTO from_tk_SetFlagValue);
  createTclCommand(interp, "cpp_setModelStatus", TCL_FUNC_PROTO from_tk_SetModelStatus);
  // Other
  createTclCommand(interp, "cpp_deleteParameters", TCL_FUNC_PROTO from_tk_ReadDeletedParamIds);
  createTclCommand(interp, "cpp_readActiveMeshIndices", TCL_FUNC_PROTO from_tk_ReadActiveMeshIndices);
  createTclCommand(interp, "cpp_readModelFileCreated", TCL_FUNC_PROTO from_tk_ReadModelFileCreated);
  createTclCommand(interp, "cpp_readModelFileModified", TCL_FUNC_PROTO from_tk_ReadModelFileModified);
  createTclCommand(interp, "cpp_readModelHasUserDefinitions", TCL_FUNC_PROTO from_tk_ReadModelHasUserDefinitions);
  createTclCommand(interp, "cpp_readModelFileTime", TCL_FUNC_PROTO from_tk_ReadModelFileTime);
  createTclCommand(interp, "cpp_readSelectionTolerances", TCL_FUNC_PROTO from_tk_ReadSelectionTolerances);
  createTclCommand(interp, "cpp_readTimestamp", TCL_FUNC_PROTO from_tk_ReadTimestamp);
  createTclCommand(interp, "cpp_resetAllBoundarySelections", TCL_FUNC_PROTO from_tk_ResetAllBoundarySelections);
  createTclCommand(interp, "cpp_resetBoundarySelections", TCL_FUNC_PROTO from_tk_ResetBoundarySelections);
  createTclCommand(interp, "cpp_setSelectionsToGui", TCL_FUNC_PROTO from_tk_SetSelectionsToGui);
  createTclCommand(interp, "cpp_selectMeshBoundaryElements", TCL_FUNC_PROTO from_tk_SelectMeshBoundaryElements);

  createTclCommand(interp, "cpp_combineBoundaries", TCL_FUNC_PROTO from_tk_CombineBoundaries);
  createTclCommand(interp, "cpp_boundaryPanelNamesOk", TCL_FUNC_PROTO from_tk_ReadBoundaryNames);
  createTclCommand(interp, "cpp_splitCombineBoundariesRedo", TCL_FUNC_PROTO from_tk_SplitCombineBoundariesRedo);
  createTclCommand(interp, "cpp_splitCombineBoundariesUndo", TCL_FUNC_PROTO from_tk_SplitCombineBoundariesUndo);
  createTclCommand(interp, "cpp_splitBoundary", TCL_FUNC_PROTO from_tk_SplitBoundary);
  createTclCommand(interp, "cpp_stopEditMeshBoundaries", TCL_FUNC_PROTO from_tk_StopEditMeshBoundaries);
  createTclCommand(interp, "cpp_restoreBoundaryNames", TCL_FUNC_PROTO from_tk_RestoreBoundaryNames);
  createTclCommand(interp, "cpp_storeBoundaryNames", TCL_FUNC_PROTO from_tk_StoreBoundaryNames);

  //------ List box selection routines ------
  createTclCommand(interp, "cpp_bodySelected", TCL_FUNC_PROTO from_tk_BodySelected);
  createTclCommand(interp, "cpp_boundarySelected", TCL_FUNC_PROTO from_tk_BoundarySelected);
  createTclCommand(interp, "cpp_boundariesSelected", TCL_FUNC_PROTO from_tk_BoundariesSelected);

  //------ Misc ------
  createTclCommand(interp, "cpp_processExists", TCL_FUNC_PROTO from_tk_ProcessExists);
  createTclCommand(interp, "cpp_processResume", TCL_FUNC_PROTO from_tk_ProcessResume);
  createTclCommand(interp, "cpp_processSetPriorityLevel", TCL_FUNC_PROTO from_tk_ProcessSetPriorityLevel);
  createTclCommand(interp, "cpp_processStart", TCL_FUNC_PROTO from_tk_ProcessStart);
  createTclCommand(interp, "cpp_processStop", TCL_FUNC_PROTO from_tk_ProcessStop);
  createTclCommand(interp, "cpp_processSuspend", TCL_FUNC_PROTO from_tk_ProcessSuspend);
  createTclCommand(interp, "cpp_readUserSettingFiles", TCL_FUNC_PROTO from_tk_ReadUserSettingFiles);
  createTclCommand(interp, "cpp_setCurrentDirectory", TCL_FUNC_PROTO from_tk_SetCurrentDirectory);
  createTclCommand(interp, "cpp_readColorFile", TCL_FUNC_PROTO from_tk_ReadColorFile);
  createTclCommand(interp, "cpp_readMatcFile", TCL_FUNC_PROTO from_tk_ReadMatcFile);
  createTclCommand(interp, "cpp_updateMatcFile", TCL_FUNC_PROTO from_tk_UpdateMatcFile);

  // Stop Front processing
  createTclCommand(interp, "cpp_doBreak", TCL_FUNC_PROTO from_tk_DoBreak);

  // Find color name for a hex-value color id
  createTclCommand(interp, "cpp_colorHex2Name", TCL_FUNC_PROTO from_tk_ColorHex2Name);

}


// *** Function creates and initialises the TCL/Tk environment.
Tcl_Interp*
UserInterface_TCL::createTclEnvironment(Hinst application)
{

  Tcl_SetPanicProc((void (*)(const char *, ...)) WishPanic);

  // *** Create Tcl-interpreter and register command procedures
  Tcl_FindExecutable( (char*) application );
  Tcl_Interp* interp = Tcl_CreateInterp();
#ifdef TCL_MEM_DEBUG
  Tcl_InitMemory(interp);
#endif

  // *** Invoke application-specific initialization.
  if (My_Tcl_AppInit(interp) != TCL_OK) {
    WishPanic("My_Tcl_AppInit failed: %s\n", interp->result);
  }

  // Result value is the Tcl interpreter
  return interp;
}


// Show error mesage box in gui
//
// NOTE: str1-str4 are for convenience only, tehy are combined into one string!
// NOTE: Use "////n" in cpp-side to deliver new lines into message-box
//
void
UserInterface_TCL::errMsg(int err_level, char* str1, char* str2, char* str3, char* str4)
{
  ostrstream strm;
  strm.unsetf(ios:: skipws);

  strm << err_level << tclArgSeparator;

  if (str1 != NULL) strm << str1;
  if (str2 != NULL) strm << str2;
  if (str3 != NULL) strm << str3;
  if (str4 != NULL) strm << str4;
  strm << ends;

  sendCommandToGui(theInterp, "Interface::displayErrorMsg", strm.str());
}


// Convert field's gui-name (like "HEAT_EQUATION") to sif-name ("Heat Equation")
void
UserInterface_TCL::fieldNameGuiToSif(const char* gui_name, char* sif_name_buffer)
{
  sendCommandToGui(theInterp, "Interface::fieldNameGuiToSif", gui_name);
  strcpy(sif_name_buffer, getCommandResults(theInterp));
}


// Convert field's sif-name (like "Heat Equation") to gui-name ("HEAT_EQUATION")
void
UserInterface_TCL::fieldNameSifToGui(const char* sif_name, char* gui_name_buffer)
{
  sendCommandToGui(theInterp, "Interface::fieldNameSifToGui", sif_name);
  strcpy(gui_name_buffer, getCommandResults(theInterp));
}


// Call Elmer MATC from Gui
int
UserInterface_TCL::from_tk_DoMatc(ClientData clientData, Tcl_Interp *interp,
                                  int argc, char* argv[])
{
  char* data = getCommandArguments(interp);

  char* result = mtc_domath(data);

  // Write result for Gui
  if ( result != NULL && result[0] != '\0' ) {
    writeVariable(interp, "Info", "MATC_result", result);
  } else {
    writeVariable(interp, "Info", "MATC_result", " ");
  }

  return TCL_OK;
}


// A body was selected in the bodylist.
int
UserInterface_TCL::from_tk_BodySelected(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  istrstream in(getCommandArguments(interp));
  int bd_id, lr_id;
  in >> bd_id >> lr_id;
  theControlCenter->selectBody(bd_id, lr_id);
  return TCL_OK;
}


// One boundary was selected in the listbox
int
UserInterface_TCL::from_tk_BoundarySelected(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  //NOTE: For testing only one element at atime.
  istrstream in(getCommandArguments(interp));

  int bd1_id, lr1_id, bd2_id, lr2_id, bndr_id, accept_body_change, update_gui;

  in >> bndr_id >> bd1_id >> lr1_id >> bd2_id >> lr2_id >> accept_body_change >> update_gui;

  theControlCenter->selectBoundary(bndr_id, bd1_id, lr1_id, bd2_id, lr2_id,
                                  (bool)accept_body_change, (bool)update_gui);
  return TCL_OK;
}


// Group of boundaries were selected in the listbox
// NOTE !!!Works currently only for one boundary!!!
int
UserInterface_TCL::from_tk_BoundariesSelected(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  // NOTE: For testing only one element at a time.
  istrstream in(getCommandArguments(interp));

  int nof_boundaries;

  in >> nof_boundaries;

  int* bndr_ids = new int[nof_boundaries];
  int* bd1_ids = new int[nof_boundaries];
  int* lr1_ids = new int[nof_boundaries];
  int* bd2_ids = new int[nof_boundaries];
  int* lr2_ids = new int[nof_boundaries];

  for (int i = 0; i < nof_boundaries; i++) {
    in >> bndr_ids[i] >> bd1_ids[i] >> lr1_ids[i] >> bd2_ids[i] >> lr2_ids[i];
  }

  int accept_body_change, update_gui;

  in >> accept_body_change >> update_gui;

  theControlCenter->selectBoundaries(nof_boundaries, bndr_ids,
                                     bd1_ids, lr1_ids,
                                     bd2_ids, lr2_ids,
                                     (bool)accept_body_change, (bool)update_gui);
  return TCL_OK;
}


// Function calls the model to check if any mesh corner elements has zero-velocity
// resulting boundary conditions etc.
int
UserInterface_TCL::from_tk_CheckMeshCornerElements(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "CheckMeshCornerElements: Command not legal. No model exists!");
    return TCL_OK;
  }

  model->classifyMeshCornerElements();

  return TCL_OK;
}


// Function ask from the model the status and saves it (also into Gui side)!
int
UserInterface_TCL::from_tk_CheckModelStatus(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "CheckModelStatus: Command not legal. No model exists!");
    return TCL_OK;
  }

  // Get current status
  ecif_modelStatus status = model->checkStatus();

  // Save in model
  model->setModelStatus(status);

  // Save in Gui
  to_tk_WriteModelStatus(interp, *model);

  return TCL_OK;
}


// Function starts boundary splitting
int
UserInterface_TCL::from_tk_CombineBoundaries(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "CombineBoundaries: Command not legal. No model exists!");
    return TCL_OK;
  }

  int body1_id, body2_id;

  istrstream in(getCommandArguments(interp));
  in >> body1_id >> body2_id;

  model->combineBoundaries(body1_id, body2_id);
  return TCL_OK;
}


// Copy parameters from an emf-file.
int
UserInterface_TCL::from_tk_CopyParameters(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  char* emf_filename = getCommandArguments(interp);

#if defined(FRONT_DEBUG)
  theControlCenter->copyParameters(emf_filename);

#else
  try {
    theControlCenter->copyParameters(emf_filename);
  }

  catch (...) {
    TheUI->errMsg(0, "CopyParameters: processing error!");
  }
#endif

  return TCL_OK;
}


// Function calls the model to correct all zero-velocity (corner) elements
int
UserInterface_TCL::from_tk_CorrectMeshZeroVelocityElements(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "CorrectMeshZeroVelocityElements: Command not legal. No model exists!");
    return TCL_OK;
  }

  model->correctMeshZeroVelocityElements();

  return TCL_OK;
}


// Function closes renderer window.
int
UserInterface_TCL::from_tk_Exit(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  theControlCenter->Exit();
  return TCL_OK;
}


// Read mesh for the (CAD) model from Elmer DB
int
UserInterface_TCL::from_tk_LoadMesh(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "LoadMesh: Command not legal. No model exists!");
    return TCL_OK;
  }

  model->loadMesh();

  return TCL_OK;
}



// Read a CAD-file.
int
UserInterface_TCL::from_tk_OpenCadFile(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  // First argument: cad filepath
  // Second argument: cad type (Elmer etc.)
  char* data = getCommandArguments(interp);
  int oc_argc;
  char** oc_argv;
  int code = Tcl_SplitList(interp, data, &oc_argc, (const char ***) &oc_argv);

  if (code != TCL_OK) return TCL_ERROR;

#if defined(FRONT_DEBUG)
  theControlCenter->readCADFile(oc_argv[0], oc_argv[1]);

#else
  try {
    theControlCenter->readCADFile(oc_argv[0], oc_argv[1]);
  }

  catch (...) {
    TheUI->errMsg(0, "Read Cad file: processing error!");
  }
#endif

  return TCL_OK;
}


// Read a mesh file.
int
UserInterface_TCL::from_tk_OpenMeshFile(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{

  // First argument is open mode: 1 = create new model, 0 do not create
  // Second argument: mesh filepath
  // Third argument: mesh type (Default, Ideas etc.)

  char* data = getCommandArguments(interp);
  int om_argc;
  char** om_argv;
  int code = Tcl_SplitList(interp, data, &om_argc, (const char ***) &om_argv);

  if (code != TCL_OK)
    return TCL_ERROR;

  int mode = atoi(om_argv[0]);

#if defined(FRONT_DEBUG)
  // Create new model
  if ( mode == 1 ) {
    theControlCenter->readMeshFile(om_argv[1], om_argv[2], true);

  // Add new mesh (or update existing)
  } else {
    theControlCenter->readMeshFile(om_argv[1], om_argv[2], false);
  }

#else
  try {
    // Create new model
    if ( mode == 1 ) {
      theControlCenter->readMeshFile(om_argv[1], om_argv[2], true);

    // Add new mesh (or update existing)
    } else {
      theControlCenter->readMeshFile(om_argv[1], om_argv[2], false);
    }
  }

  catch (...) {
    TheUI->errMsg(0, "Read mesh file: processing error!");
  }
#endif

  return TCL_OK;
}


// Reads a ecif model file.
int
UserInterface_TCL::from_tk_OpenModelFile(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  char* in_filename = getCommandArguments(interp);

  char* auto_load_mesh = (char *) Tcl_GetVar2(interp, "UserSetting", "AUTO_LOAD_MESH",  glob_flag);

  bool load_mesh = true;

  if (auto_load_mesh != NULL)
    load_mesh = atoi(auto_load_mesh);

  theControlCenter->readModelFile(in_filename, load_mesh);

  return TCL_OK;
}


// Function checks if process  with id 'nbr' exists
int
UserInterface_TCL::from_tk_ProcessExists(ClientData clientData, Tcl_Interp *interp,
                                       int argc, char* argv[])
{
  int nbr = atoi(getCommandArguments(interp));

  bool exists = theControlCenter->processExists(nbr);

  strstream strm;
  strm << nbr << ",exists" << ends;
  writeVariable(interp, "ProcessTable", strm.str(), exists);

  return TCL_OK;
}


// Function restarts a suspended external process with id 'nbr'
int
UserInterface_TCL::from_tk_ProcessResume(ClientData clientData, Tcl_Interp *interp,
                                       int argc, char* argv[])
{
  int nbr = atoi(getCommandArguments(interp));

  theControlCenter->processResume(nbr);

  return TCL_OK;
}


// Function restarts a suspended external process with id 'nbr'
int
UserInterface_TCL::from_tk_ProcessSetPriorityLevel(ClientData clientData, Tcl_Interp *interp,
                                       int argc, char* argv[])
{
  char* data = getCommandArguments(interp);
  int pl_argc;
  char** pl_argv;
  int code = Tcl_SplitList(interp, data, &pl_argc, (const char ***)&pl_argv);

  if (code != TCL_OK)
    return TCL_ERROR;

  int nbr = atoi(pl_argv[0]);

  priorityLevel priority;
  char* plevel = pl_argv[1];
  if ( 0 == strcmp(plevel, "LOW_PRIORITY") )
    priority = ECIF_LOW_PRIORITY;
  else if ( 0 == strcmp(plevel, "LOWER_THAN_NORMAL_PRIORITY") )
    priority = ECIF_LOWER_THAN_NORMAL_PRIORITY;
  else if ( 0 == strcmp(plevel, "NORMAL_PRIORITY") )
    priority = ECIF_NORMAL_PRIORITY;
  else if ( 0 == strcmp(plevel, "HIGER_THAN_NORMAL_PRIORITY") )
    priority = ECIF_HIGHER_THAN_NORMAL_PRIORITY;
  else if ( 0 == strcmp(plevel, "HIGH_PRIORITY") )
    priority = ECIF_HIGH_PRIORITY;

  theControlCenter->processSetPriorityLevel(nbr, priority);

  return TCL_OK;
}

// Function starts external process
int
UserInterface_TCL::from_tk_ProcessStart(ClientData clientData, Tcl_Interp *interp,
                                        int argc, char* argv[])
{
  //---Get process start-data
  char* data = getCommandArguments(interp);
  int pr_argc;
  char** pr_argv;
  int code = Tcl_SplitList(interp, data, &pr_argc, (const char ***) &pr_argv);

  if (code != TCL_OK)
    return TCL_ERROR;

  char* command = pr_argv[0];
  char* args = pr_argv[1];
  int nbr = atoi(pr_argv[2]);
  char* name = pr_argv[3];

  priorityLevel priority;
  char* plevel = pr_argv[4];
  if ( 0 == strcmp(plevel, "LOW_PRIORITY") )
    priority = ECIF_LOW_PRIORITY;
  else if ( 0 == strcmp(plevel, "LOWER_THAN_NORMAL_PRIORITY") )
    priority = ECIF_LOWER_THAN_NORMAL_PRIORITY;
  else if ( 0 == strcmp(plevel, "NORMAL_PRIORITY") )
    priority = ECIF_NORMAL_PRIORITY;
  else if ( 0 == strcmp(plevel, "HIGER_THAN_NORMAL_PRIORITY") )
    priority = ECIF_HIGHER_THAN_NORMAL_PRIORITY;
  else if ( 0 == strcmp(plevel, "HIGH_PRIORITY") )
    priority = ECIF_HIGH_PRIORITY;

  char* logfile = pr_argv[5];
  bool show_console = false;

  if ( 0 == strcmp(logfile, "none") )
    logfile = NULL;
  else if ( 0 == strcmp(logfile, "shell") ) {
    logfile = NULL;
    show_console = true;
  }

  Process* process = new Process(command, args, nbr, name, priority, show_console, logfile);

  // Try to create and start process
  if ( theControlCenter->processStart(process) ) {

    Hfile outputHandle = process->getOutputHandle();
    ProcessId PID = process->getProcessId();

    // Set logfile handle
    // NOTE: use 0 here instead of NULL, because
    // handle is int
    if (outputHandle != 0) {

      Tcl_Channel channel = NULL;
      Timer ch_timer;
      ch_timer.start();

      // Try to open, max 30 seconds
      while( channel == NULL && ch_timer.getLapTime(WALL_TIME) < 30 ) {
        channel = createFileChannel(interp, outputHandle);
      }

      ch_timer.stop();

      if (channel != NULL) {
        char* channel_name = (char *) Tcl_GetChannelName(channel);
        writeIdVariable(interp, "ProcessTable", process->ID(), "channel", channel_name);

      } else {
        strstream strm;
        strm << "Could not open log file for the process " << process->ID() << ends;
        theControlCenter->getGui()->showMsg(strm.str());
      }
    }

    // Set process id (PID)
    if ( PID != 0 ) {
      writeIdVariable(interp, "ProcessTable", process->ID(), "pid", PID);
    }
  }

  return TCL_OK;
}


// Function stops external process with id 'nbr'
int
UserInterface_TCL::from_tk_ProcessStop(ClientData clientData, Tcl_Interp *interp,
                                       int argc, char* argv[])
{
  int nbr = atoi(getCommandArguments(interp));

  theControlCenter->processStop(nbr);

  return TCL_OK;
}


// Function suspends external process with id 'nbr'
int
UserInterface_TCL::from_tk_ProcessSuspend(ClientData clientData, Tcl_Interp *interp,
                                       int argc, char* argv[])
{
  int nbr = atoi(getCommandArguments(interp));

  theControlCenter->processSuspend(nbr);

  return TCL_OK;
}



// Gets model status message from the model
int
UserInterface_TCL::from_tk_PutModelStatusMessage(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "PutModelStatusMessage: Command not legal. No model exists!");
    return TCL_OK;
  }

  to_tk_WriteModelStatusMessage(interp, *model);

  return TCL_OK;
}


// Gets active Solver mesh indices from gui
int
UserInterface_TCL::from_tk_ReadActiveMeshIndices(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadActiveMeshIndices: Command not legal. No model exists!");
    return TCL_OK;
  }

  int nof_meshes;
  int* mesh_indices = NULL;

  readVariable(interp, "Model", "activeMeshIndices", nof_meshes, mesh_indices);

  model->setActiveMeshIndices(nof_meshes, mesh_indices);

  delete[] mesh_indices;

  return TCL_OK;
}


// Function reads and updates body-colors data when OK-button
// is pressed in Tk/GUI body-names dialog.
int
UserInterface_TCL::from_tk_ReadBodyColors(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadBodyColors: Command not legal. No model exists!");
    return TCL_OK;
  }

  Color4 color;
  color[3] = MAX_NOF_COLOR_LEVELS;

  int index = 0;
  while (true) {

    Body* body = model->getBody(index++);

    if (body==NULL) break;

    char* color_buffer = NULL;

    int id = body->Id();

    readIdVariable(interp, "ObjectTable", id, "clr", color_buffer);

    if ( color_buffer == NULL ) continue;

    // Skip leading Tcl # in color hex-values
    char* tmp = color_buffer;
    tmp++;

    model->rgbHex2Color(6, tmp, color);
    body->setColor(color);

    delete[] color_buffer;
  }

  model->refreshRenderer();

  return TCL_OK;
}


// Function reads and updates body data when OK-button
// is pressed in Tk/GUI edit bodies dialog.
int
UserInterface_TCL::from_tk_ReadBodyData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadBodyData: Command not legal. No model exists!");
    return TCL_OK;
  }

  from_tk_ReadBodyColors(clientData, interp, argc, argv);
  from_tk_ReadBodyNames(clientData, interp, argc, argv);

  to_tk_WriteBodyLayerData(interp, *model);

  to_tk_WriteBodyData(interp, *model);

  return TCL_OK;
}


// Function reads body id to be deleted

int
UserInterface_TCL::from_tk_ReadBodyDeleteData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadBodyDeleteData: Command not legal. No model exists!");
    return TCL_OK;
  }

  int body_id = NO_INDEX;

  readVariable(interp, "BodyProperty", "objectId", body_id);

  Body* body = model->getBodyById(body_id);

  if ( body != NULL ) {
    model->removeBody(body);

    to_tk_WriteBodyLayerData(interp, *model);
    to_tk_WriteBodyData(interp, *model);
    to_tk_WriteBoundaryData(interp, *model);
    to_tk_WriteElementGroupData(theInterp, *model);
  }

  return TCL_OK;
}



// Read body display list (from the BodyDisplay panel).
int
UserInterface_TCL::from_tk_ReadBodyDisplayData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadSelectedBodies: Command not legal. No model exists!");
    return TCL_OK;
  }

  int index = 0;
  while (true) {

    Body* body = model->getBody(index++);

    if (body==NULL) break;

    int id = body->Id();
    int display = 1;

    readIdVariable(interp, "ObjectTable", id, "dspl", display);

    if (display == 0) {
      body->setDrawMode(DM_HIDDEN);
    } else {
      body->setDrawMode(DM_NORMAL);
    }
  }

  model->refreshRenderer();

  return TCL_OK;
}


// Function reads and updates body-force data when OK-button
// is pressed in Tk/GUI body-forces dialog.
int
UserInterface_TCL::from_tk_ReadBodyForceData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadBodyForceData: Command not legal. No model exists!");
    return TCL_OK;
  }

  setParameterData(model, ECIF_BODY_FORCE, interp, "BodyForce");

  //----ObjectTable size
  int size;
  int* ids = NULL;

  readVariable(interp, "ObjectTable", "ids", size, ids);

  for (int i = 0; i < size; i++) {

    int id = ids[i];

    char* tmp;
    readIdVariable(interp, "ObjectTable", id, "tp", tmp);

    //--Skip if not a body target group
    if ( 0 != strcmp(tmp, "B") )
      continue;

    int pid;
    readIdVariable(interp, "ObjectTable", id, "bf", pid);

    //--Update body target group's body force-info
    Body* bd = model->getBodyById(id);

    if ( bd != NULL ) {
      bd->setBodyForceId(pid);
    }
  }

  model->updateParametersApplyCounts(ECIF_BODY_FORCE);

  to_tk_WriteStatusBodyForces(interp, *model);

  delete[] ids;

  return TCL_OK;
}


// Function reads and updates body-parameter data when OK-button
// is pressed in Tk/GUI body-parameter dialog.
int
UserInterface_TCL::from_tk_ReadBodyParameterData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadBodyParameterData: Command not legal. No model exists!");
    return TCL_OK;
  }

  setParameterData(model, ECIF_BODY_PARAMETER, interp, "BodyParameter");

  //----ObjectTable size
  int size;
  char** ids = NULL;

  readVariable(interp, "ObjectTable", "ids", size, ids);

  for (int i = 0; i < size; i++) {

    char* oid = ids[i];
    int id = atoi(oid);

    char* tmp;
    readIdVariable(interp, "ObjectTable", oid, "tp", tmp);

    //--Skip if not a body
    if ( 0 != strcmp(tmp, "B") )
      continue;

    int pid;
    readIdVariable(interp, "ObjectTable", oid, "bodyp", pid);

    //--Update body's parameter-info
    Body* body = model->getBodyById(id);

    if ( body != NULL ) {
      body->setBodyParameterId(pid);
    }
  }

  for (int j = 0; j < size; j++) {
    delete[] ids[j];
  }
  delete[] ids;

  return TCL_OK;
}


// Function reads and updates body-names data when OK-button
// is pressed in Tk/GUI body-names dialog.
int
UserInterface_TCL::from_tk_ReadBodyNames(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadBodyNames: Command not legal. No model exists!");
    return TCL_OK;
  }

  int index = 0;
  while (true) {
    Body* body = model->getBody(index++);
    if (body==NULL) break;
    char* name = NULL;
    int id = body->Id();
    readIdVariable(interp, "ObjectTable", id, "nm", name);

    body->setName(name);

    // Update also layer name if it is an implicit (logical) layer
    BodyLayer* bl = (BodyLayer*)body->getLayer(0);
    if ( bl != NULL && bl->getLayerType() == IMPLICIT_LAYER ) {
      bl->setName(name);
      writeIdVariable(interp, "ObjectTable", bl->Id(), "nm", name);
    }

    delete[] name;
  }

  return TCL_OK;
}


// Function updates mesh boundary vertices after boundaries panel split/combine
// operation when OK-button is pressed in Edit/Boundaries panel.
int
UserInterface_TCL::from_tk_ReadBoundariesData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadBoundariesData: Command not legal. No model exists!");
    return TCL_OK;
  }

  model->setFlagValue(GEOMETRY_EDITED, GEOMETRY_EDITED_BOUNDARIES, true);

  if ( !model->getFlagValue(GEOMETRY_TYPE_CAD) ) {
    model->updateBoundaries();
    return TCL_OK;
  }

  // Read discretation settings;
  int cnt = 0;
  linDeltaType* types = NULL;
  double* valuesU = NULL;
  int* useFixed = NULL;
  bool* useFixedN = NULL;

  int nof_vals;
  int nof_fns;

  char* type_str = NULL;

  int index = 0;
  while (true) {
    BodyElement* be = model->getBoundary(index++);
    if (be==NULL) break;
    int id = be->Id();

    readIdVariable(interp, "ObjectTable", id, "nofCmp", cnt);
    readIdVariable(interp, "ObjectTable", id, "dscTp", type_str);
    readIdVariable(interp, "ObjectTable", id, "dscU", nof_vals, valuesU);
    readIdVariable(interp, "ObjectTable", id, "useFN", nof_fns, useFixed);
    
    int nof_types = 0;
    
    if ( type_str != NULL ) {
      nof_types = strlen(type_str);
    }

    // Recode type characters
    //
    if ( nof_types > 0 ) {
      types = new linDeltaType[nof_types];

       for (int i = 0; i < nof_types; i++) {
        if ( type_str[i] == 'H' ) {
          types[i] = LIN_DELTA_H;
        } else if ( type_str[i] == 'N' ) {
          types[i] = LIN_DELTA_N;
        } else if ( type_str[i] == 'U' ) {
          types[i] = LIN_DELTA_U;
        } else {
          types[i] = LIN_DELTA_NONE;
        }
      }
    }
    
    // If data is technically correct, apply it
    //
    if ( (nof_vals == 0  || nof_vals == cnt) &&
         (nof_types == 0 || nof_types == cnt) &&
          nof_fns == cnt
       ) {

      useFixedN = new bool[nof_fns];

      for (int i = 0; i < nof_fns; i++) {
        useFixedN[i] = (bool)useFixed[i];
      }

      be->setDiscretizationData(cnt, types, valuesU, NULL, useFixedN);
    }

    // Delete work arreis for next call
    //
    delete[] types; types = NULL;
    delete[] valuesU; valuesU = NULL;
    delete[] useFixedN; useFixedN = NULL;

    delete[] type_str; type_str = NULL;
    delete[] useFixed; useFixed = NULL;
  }

  model->updateCadGeometry();

  return TCL_OK;
}


// Function reads and updates boundary conditions when OK-button
// is pressed in Tk/GUI boundary-condition dialog.
int
UserInterface_TCL::from_tk_ReadBoundaryConditionData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{

  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadBoundaryConditionData: Command not legal. No model exists!");
    return TCL_OK;
  }

  //---Reset all boundary condition attachments
  model->resetBoundaryConditions();

  //---Boundary condition data
  setParameterData(model, ECIF_BOUNDARY_CONDITION, interp, "BoundaryCondition");

  MultiIdTable bcTable;

  //---Update boundary element's condition ids
  int size;
  int* ids = NULL;

  readVariable(interp, "ObjectTable", "ids", size, ids);

  for (int i = 0; i < size; i++) {

    char* tmp;
    int id = ids[i];

    int pid;
    readIdVariable(interp, "ObjectTable", id, "bc", pid);

    readIdVariable(interp, "ObjectTable", id, "tp", tmp);

    // Read bc from boundary groups
    if ( 0 == strcmp(tmp, "FG") ||
         0 == strcmp(tmp, "EG") ||
         0 == strcmp(tmp, "VG")
       ) {
      bcTable.insert( std::pair<int const, int>(pid, id));
    } else {
      continue;
    }
  }

  setBoundaryConditions(interp, model, bcTable);

  // Update parameter apply count
  model->updateParametersApplyCounts(ECIF_BOUNDARY_CONDITION);

  // Update model's DiffuseGray radiation info
  model->checkDiffuseGrayRadiation();

  to_tk_WriteStatusBoundaryConditions(interp, *model);

  delete[] ids;

  return TCL_OK;
}


// Read boundary display list
int
UserInterface_TCL::from_tk_ReadBoundaryDisplayData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadBoundaryDisplayData: Command not legal. No model exists!");
    return TCL_OK;
  }

  int index = 0;
  while (true) {
    BodyElement* be = model->getBoundary(index++);
    if (be==NULL) break;
    int id = be->Id();
    int display = 1;
    readIdVariable(interp, "ObjectTable", id, "dspl", display);
    if (display == 0) {
      be->setDrawMode(DM_HIDDEN);
    } else {
      be->setDrawMode(DM_NORMAL);
    }
  }

  model->refreshRenderer();

  return TCL_OK;
}


// Function reads and updates boundary name data when OK-button
// is pressed in Tk/GUI edit boundaries dialog.
int
UserInterface_TCL::from_tk_ReadBoundaryNames(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadBoundaryNames: Command not legal. No model exists!");
    return TCL_OK;
  }

  int index = 0;
  while (true) {
    BodyElement* be = model->getBoundary(index++);
    if (be==NULL) break;
    char* name = NULL;
    readIdVariable(interp, "ObjectTable", be->Id(), "nm", name);

    be->setName(name);

    // Update also group name if it is an implicit (logical) group
    BodyElementGroup* beg = (BodyElementGroup*)be->getElementGroup();
    if ( beg != NULL && beg->getGroupType() == IMPLICIT_GROUP ) {
      beg->setName(name);
      writeIdVariable(interp, "ObjectTable", beg->Id(), "nm", name);
    }

    delete[] name;
  }

  return TCL_OK;
}


// Function reads and updates boundary parameters when OK-button
// is pressed in Tk/GUI boundary-parameter dialog.
int
UserInterface_TCL::from_tk_ReadBoundaryParameterData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadBoundaryParameterData: Command not legal. No model exists!");
    return TCL_OK;
  }

  //---Boundary parameter data
  setParameterData(model, ECIF_BOUNDARY_PARAMETER, interp, "BoundaryParameter");

  //---Update boundary element's boundary parameter id
  int size;
  char** ids = NULL;

  readVariable(interp, "ObjectTable", "ids", size, ids);

  for (int i = 0; i < size; i++) {

    char* tmp;
    char* oid = ids[i];
    int id = atoi(oid);

    readIdVariable(interp, "ObjectTable", id, "tp", tmp);

    if ( 0 != strcmp(tmp, "F") &&
         0 != strcmp(tmp, "E") &&
         0 != strcmp(tmp, "V")
       ) {
      continue;
    }

    int pid;
    readIdVariable(interp, "ObjectTable", id, "bndrp", pid);

    //--Update body's parameter-info
    BodyElement* be = model->getBodyElementById(id);

    if ( be != NULL ){
      be->setBoundaryParameterId(pid);
    }
  }

  for (int j = 0; j < size; j++) {
    delete[] ids[j];
  }
  delete[] ids;

  return TCL_OK;
}


// Function reads calulator solvers data from Tk.
int
UserInterface_TCL::from_tk_ReadCalculatorData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadCalculatorData: Command not legal. No model exists!");
    return TCL_OK;
  }

  setParameterData(model, ECIF_CALCULATOR, interp, "Calculator");

  return TCL_OK;
}


// Activate a color file read via Tk.
int
UserInterface_TCL::from_tk_ReadColorFile(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  //---Get tcl-data
  char* fname = getCommandArguments(interp);

  Model::readColorFile(fname);

  return TCL_OK;
}


// Function reads physical constants data from Tk.
int
UserInterface_TCL::from_tk_ReadConstantData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadConstantData: Command not legal. No model exists!");
    return TCL_OK;
  }

  setParameterData(model, ECIF_CONSTANT, interp, "Constant");

  return TCL_OK;
}


// Function reads and updates equation data when OK-button
// is pressed in Tk/GUI equations dialog.
int
UserInterface_TCL::from_tk_ReadConvertedEquationData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadEquationData: Command not legal. No model exists!");
    return TCL_OK;
  }

  setParameterData(model, ECIF_EQUATION, interp, "Equation");
  model->updateParametersApplyCounts(ECIF_EQUATION);

  return TCL_OK;
}


// Function reads coordinate system definition data from Tk.
int
UserInterface_TCL::from_tk_ReadCoordinateData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadCoordinateData: Command not legal. No model exists!");
    return TCL_OK;
  }


  setParameterData(model, ECIF_COORDINATE, interp, "Coordinate");

  return TCL_OK;
}


// Function reads datafile (result and input files) definition data from Tk.
int
UserInterface_TCL::from_tk_ReadDatafileData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadDatafileData: Command not legal. No model exists!");
    return TCL_OK;
  }

  setParameterData(model, ECIF_DATAFILE, interp, "Datafile");

  return TCL_OK;
}


// Function reads ids for deleted parameters when Ok-button is pressed in a panel
// where parameters were deleted.
int
UserInterface_TCL::from_tk_ReadDeletedParamIds(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  //---Get data-string  as a "structured" list from Tcl.
  //--First part of the data is panelt-type
  //--Second part is the list of param-ids which were deleted
  char* deletedList = getCommandArguments(interp);

  //---Break this list into an array of strings, where the strings are
  //the two parts mentioned above.
  int all_argc;
  char** all_argv;
  int code = Tcl_SplitList(interp, deletedList, &all_argc, (const char ***) &all_argv);
  if (code != TCL_OK)
    return TCL_ERROR;

  //---Panel type
   ecif_parameterType param_type;
   char* pt = all_argv[0];
   if (0 == strcmp("IC", pt))
      param_type = ECIF_INITIAL_CONDITION;
   else if (0 == strcmp("BC", pt))
      param_type = ECIF_BOUNDARY_CONDITION;
   else if (0 == strcmp("EQ", pt))
      param_type = ECIF_EQUATION;
   else if (0 == strcmp("BF", pt))
      param_type = ECIF_BODY_FORCE;
   else if (0 == strcmp("MP", pt))
      param_type = ECIF_MATERIAL;
   else if (0 == strcmp("GR", pt))
      param_type = ECIF_GRID_PARAMETER;
   else if (0 == strcmp("GH", pt))
      param_type = ECIF_GRID_H;

   //---Param ids
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadDeletedParamIds: Command not legal. No model exists!");
    return TCL_OK;
  }

  int ids_argc;
  char** ids_argv;
  code = Tcl_SplitList(interp, all_argv[1], &ids_argc, (const char ***) &ids_argv);
  if (code != TCL_OK)
    return TCL_ERROR;
  for (int i=0; i<ids_argc; i++) {
    istrstream in(ids_argv[i]);
    int param_id;
    in >> param_id;
    model->deleteParameter(param_type, param_id);
  }

  // Tcl_SplitList allocates dynamic memory for its argument.
  // We should free this memory, but the following call causes an error !
  //  free((char*) largv);
  return TCL_OK;
}


// Function reads and updates equation data when OK-button
// is pressed in Tk/GUI equations dialog.
int
UserInterface_TCL::from_tk_ReadEquationData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadEquationData: Command not legal. No model exists!");
    return TCL_OK;
  }

  setParameterData(model, ECIF_EQUATION, interp, "Equation");

  //---ObjectTable size
  int size;
  int* ids = NULL;

  readVariable(interp, "ObjectTable", "ids", size, ids);

  for (int i = 0; i < size; i++) {

    int id = ids[i];

    char* tmp;
    readIdVariable(interp, "ObjectTable", id, "tp", tmp);

    //--Skip if not a body
    if ( 0 != strcmp(tmp, "B") )
      continue;

    int pid;

    //--Read variable values
    readIdVariable(interp, "ObjectTable", id, "eq", pid);

    //--Update body's equation-info
    Body* bd = model->getBodyById(id);

    if (bd != NULL ) {
      bd->setEquationId(pid);
    }

  }

  model->updateParametersApplyCounts(ECIF_EQUATION);

  to_tk_WriteStatusEquations(interp, *model);
  to_tk_WriteStatusBodyForces(interp, *model);
  to_tk_WriteStatusBoundaryConditions(interp, *model);
  to_tk_WriteStatusInitialConditions(interp, *model);
  to_tk_WriteStatusMaterials(interp, *model);
  to_tk_WriteStatusTimesteps(interp, *model);
  to_tk_WriteStatusMeshes(interp, *model);

  return TCL_OK;
}


// Function reads equation variables (like Advection Diffusion variables) definition data from Tk.
int
UserInterface_TCL::from_tk_ReadEquationVariablesData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadEquationVariablesData: Command not legal. No model exists!");
    return TCL_OK;
  }


  setParameterData(model, ECIF_EQUATION_VARIABLE, interp, "EquationVariable");

  return TCL_OK;
}


// Function reads and updates grid-param data when OK-button
// is pressed in Tk/GUI grid-control dialog.
int
UserInterface_TCL::from_tk_ReadGridHData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{

  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadGridHData: Command not legal. No model exists!");
    return TCL_OK;
  }

  //----GridH data
  setParameterData(model, ECIF_GRID_H, interp, "GridH");

  //----GridH ids
  int size = 0;
  char** ids = NULL;

  readVariable(interp, "ObjectTable", "ids", size, ids);

  for (int i = 0; i < size; i++) {

    char* oid = ids[i];
    int id = atoi(oid);

    char* tmp;
    readIdVariable(interp, "ObjectTable", oid, "tp", tmp);

    //--Skip if not a bodyelement
    if ( 0 != strcmp(tmp, "F") &&
         0 != strcmp(tmp, "E") &&
         0 != strcmp(tmp, "V")
       ) {
      continue;
    }

    // Read variable values
    int nof_gids = 0;
    int* gids;
    int* mids;
    int pid, bd_lid, bid;

    readIdVariable(interp, "ObjectTable", oid, "ghIds", nof_gids, gids);       // grid paramter ids
    readIdVariable(interp, "ObjectTable", oid, "ghMshIndcs", nof_gids, mids);  // mesh indices

    BodyElement* be = model->getBodyElementById(id);

    if ( be != NULL ) {
      be->setGridHData(nof_gids, gids, mids);
    }

  }

  return TCL_OK;
}


// Function reads and updates grid-param data when OK-button
// is pressed in Tk/GUI grid-param dialog.
int
UserInterface_TCL::from_tk_ReadGridParameterData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadGridParameterData: Command not legal. No model exists!");
    return TCL_OK;
  }

  //----GridParameter data
  setParameterData(model, ECIF_GRID_PARAMETER, interp, "GridParameter");

  //----Body ids for GridParameters
  int size = 0;
  //int* ids = NULL;
  int* ids = NULL;

  readVariable(interp, "ObjectTable", "ids", size, ids);

  for (int i = 0; i < size; i++) {

    int id = ids[i];

    char* tmp;
    readIdVariable(interp, "ObjectTable", id, "tp", tmp);

    //--Skip if not a body layer
    //
    if ( 0 != strcmp(tmp, "BL") )  continue;

    // Read variable values
    int nof_gids, nof_mids, nof_xids;
    int* gids = NULL; // Grid parameter ids
    int* mids = NULL; // Grid mesh indices
    int* xids = NULL; // Exclude mesh indices

    readIdVariable(interp, "ObjectTable", id, "grIds", nof_gids, gids);        // grid paramter ids
    readIdVariable(interp, "ObjectTable", id, "grMshIndcs", nof_mids, mids);  // mesh indices
    readIdVariable(interp, "ObjectTable", id, "excldMshIndcs", nof_xids, xids);  // excluded mesh indices

    BodyLayer* lr = model->getBodyLayerById(id);

    if ( lr != NULL ){
      lr->setGridParameterData(nof_gids, gids, mids);
      lr->setExcludedMeshData(nof_xids, xids);
    }

    delete[] gids;
    delete[] mids;
    delete[] xids;
  }

  //to_tk_WriteStatusMeshParameter(interp, *model);

  delete[] ids;

  return TCL_OK;
}


// Function reads and updates body initial condition when OK-button
// is pressed in Tk/GUI init-condition dialog.
int
UserInterface_TCL::from_tk_ReadInitialConditionData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadInitialConditionData: Command not legal. No model exists!");
    return TCL_OK;
  }

  setParameterData(model, ECIF_INITIAL_CONDITION, interp, "InitialCondition");

  // ObjectTable data
  int size;
  int* ids = NULL;

  readVariable(interp, "ObjectTable", "ids", size, ids);

  for (int i = 0; i < size; i++) {

    int id = ids[i];

    char* tmp;
    readIdVariable(interp, "ObjectTable", id, "tp", tmp);

    //--Skip if not a body
    if ( 0 != strcmp(tmp, "B") )
      continue;

    // Read variable values
    int pid;
    readIdVariable(interp, "ObjectTable", id, "ic", pid);

    //--Update body's condition-info
    Body* bd = model->getBodyById(id);

    if ( bd != NULL ) {
      bd->setInitialConditionId(pid);
    }
  }

  model->updateParametersApplyCounts(ECIF_INITIAL_CONDITION);

  to_tk_WriteStatusInitialConditions(interp, *model);

  return TCL_OK;
}


// Activate a matc variable table file read via Tk.
// NOTE: No model is needed
int
UserInterface_TCL::from_tk_ReadMatcFile(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  //---Get tcl-data
  char* fname = getCommandArguments(interp);

  Model::readMatcFile(fname, NULL, false);

  return TCL_OK;
}


// Function reads and updates body material parameter data when OK-button
// is pressed in Tk/GUI material-parameters dialog.
int
UserInterface_TCL::from_tk_ReadMaterialData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadMaterialData: Command not legal. No model exists!");
    return TCL_OK;
  }

  setParameterData(model, ECIF_MATERIAL, interp, "Material");

  //--ObjectTable size
  int size;
  int* ids = NULL;

  readVariable(interp, "ObjectTable", "ids", size, ids);

  for (int i = 0; i < size; i++) {

    int id = ids[i];

    char* tmp;
    readIdVariable(interp, "ObjectTable", id, "tp", tmp);

    //--Skip if not a body target group
    if ( 0 != strcmp(tmp, "B") )
      continue;

    // Read variable values
    int pid;
    readIdVariable(interp, "ObjectTable", id, "mt", pid);

    //--Update body's material-info
    Body* bd = model->getBodyById(id);

    if ( bd != NULL ) {
      bd->setMaterialId(pid);
    }
  }

  model->updateParametersApplyCounts(ECIF_MATERIAL);

  to_tk_WriteStatusMaterials(interp, *model);

  return TCL_OK;
}


// Function reads model created info from Tk.
int
UserInterface_TCL::from_tk_ReadModelFileCreated(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadModelFileCreated: Command not legal. No model exists!");
    return TCL_OK;
  }

  //---Get tcl-data
  char* data = getCommandArguments(interp);

  //---Model file created info
  model->setModelFileCreated(data);

  return TCL_OK;
}


// Function reads model modified info from Tk.
int
UserInterface_TCL::from_tk_ReadModelFileModified(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadModelFileModified: Command not legal. No model exists!");
    return TCL_OK;
  }

  //---Get tcl-data
  char* data = getCommandArguments(interp);

  //---Model file modified info
  model->setModelFileModified(data);

  return TCL_OK;
}


// Function reads model outfile save timestamp from Tk.
int
UserInterface_TCL::from_tk_ReadModelFileTime(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadModelFileTime: Command not legal. No model exists!");
    return TCL_OK;
  }

  //---Get tcl-data
  char* data = getCommandArguments(interp);

  //---Model outfile timestamp
  model->setModelFileTime(data);
  return TCL_OK;
}


// Function reads model modified info from Tk.
int
UserInterface_TCL::from_tk_ReadModelHasUserDefinitions(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadModelHasUserDefinitions: Command not legal. No model exists!");
    return TCL_OK;
  }

  //---Get tcl-data
  bool has_defs = (bool) atoi(getCommandArguments(interp));

  //---Model file modified info
  model->setModelHasUserDefinitions(has_defs);

  return TCL_OK;
}


// Function reads and updates model level meshH when ok is pressed
// is pressed in Tk/GUI meshdefine dialog.
int
UserInterface_TCL::from_tk_ReadMeshDefineData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  int i;
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadModelMeshHData: Command not legal. No model exists!");
    return TCL_OK;
  }

  // Read mesh names and control
  // ===========================
  int nof_meshes = 0;
  char** mesh_names = NULL;
  double* mesh_hs = NULL;
  double* mesh_fs = NULL;

  readVariable(interp, "Model", "meshNames", nof_meshes, mesh_names);
  readVariable(interp, "Model", "meshHs", nof_meshes, mesh_hs);
  readVariable(interp, "Model", "meshFs", nof_meshes, mesh_fs);

  model->setMeshNames(nof_meshes, mesh_names);
  model->setMeshHs(nof_meshes, mesh_hs);
  model->setMeshFs(nof_meshes, mesh_fs);

  for (i = 0; i < nof_meshes; i++) {
    delete[] mesh_names[i];
  }

  delete[] mesh_names;
  delete[] mesh_hs;
  delete[] mesh_fs;

  // Read background mesh file info
  // ==============================
  //
  // NOTE: Bg-files info is stored (in indexed form) in
  // MeshDefine-array (Model-array contains equivalen "sparse" lists
  // (except indices lists, which is always non-sparse)
  int nof_files = 0;
  int* mesh_bg_file_indices = NULL;
  char** mesh_bg_files = NULL;
  int* mesh_bg_act = NULL;
  int* mesh_bg_ctrl = NULL;

  readVariable(interp, "Model", "meshBgMeshFileIndices", nof_files, mesh_bg_file_indices);

  readVariable(interp, "MeshDefine", "meshBgMeshFiles", nof_files, mesh_bg_files);
  readVariable(interp, "MeshDefine", "meshBgMeshActives", nof_files, mesh_bg_act);
  readVariable(interp, "MeshDefine", "meshBgMeshControls", nof_files, mesh_bg_ctrl);

  if ( mesh_bg_file_indices == NULL ) {
    nof_files = 0;
  }

  model->setMeshBgMeshFileIndices(nof_files, mesh_bg_file_indices);
  model->setMeshBgMeshFiles(nof_files, mesh_bg_files);
  model->setMeshBgMeshActives(nof_files, mesh_bg_act);
  model->setMeshBgMeshControls(nof_files, mesh_bg_ctrl);

  for (i = 0; i < nof_files; i++) {
    delete[] mesh_bg_files[i];
  }

  delete[] mesh_bg_file_indices;
  delete[] mesh_bg_files;
  delete[] mesh_bg_act;
  delete[] mesh_bg_ctrl;

  //--Set current mesh info
  int mesh_index = 0;
  readVariable(interp, "Model", "currentMeshIndex", mesh_index);
  model->setCurrentMeshIndex(mesh_index);

  char* dir = NULL;
  char* dir_abs = NULL;

  //--Mesh directory
  readVariable(interp, "ModelProperty", "CURRENT_MESH_DIRECTORY", dir);
  readVariable(interp, "ModelProperty", "CURRENT_MESH_DIRECTORY,absolute", dir_abs);

  model->setModelMeshDirectory(dir);
  model->setModelMeshDirectoryAbs(dir_abs);

  delete[] dir;
  delete[] dir_abs;

  return TCL_OK;
}


// Function reads user defined model parameters data from Tk.
int
UserInterface_TCL::from_tk_ReadModelParameterData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadModelParameterData: Command not legal. No model exists!");
    return TCL_OK;
  }

  setParameterData(model, ECIF_MODEL_PARAMETER, interp, "ModelParameter");

  return TCL_OK;
}


// Function reads and updates model properties (paths etc.)
int
UserInterface_TCL::from_tk_ReadModelPropertyData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadModelPropertyData: Command not legal. No model exists!");
    return TCL_OK;
  }

  // Init reference buffers and flags
  char* name1 = NULL; // Model name
  char* name2 = NULL; // Problem name

  char* desc1 = NULL;
  char* desc2 = NULL;

  char* dir1 = NULL; // Model Directory
  char* dir2 = NULL; // Include Path
  char* dir3 = NULL; // Results Directory
  char* dir4 = NULL; // Log Directory

  char* dir1_abs = NULL;
  char* dir2_abs = NULL;
  char* dir3_abs = NULL;
  char* dir4_abs = NULL;

  int dir2_save = 0;
  int dir3_save = 0;
  int dir4_save = 0;

  //--Names
  readVariable(interp, "ModelProperty", "MODEL_NAME", name1);
  readVariable(interp, "ModelProperty", "PROBLEM_NAME", name2);
  model->setModelNames(name1, name2);

  //--Model and problem descriptions
  readVariable(interp, "ModelProperty", "MODEL_DESCRIPTION", desc1);
  readVariable(interp, "ModelProperty", "PROBLEM_DESCRIPTION", desc2);
  model->setModelDescriptions(desc1, desc2);

  //--Model directory, Mesh directory, Include path, results directory and temporary files directory info
  readVariable(interp, "ModelProperty", "MODEL_DIRECTORY", dir1);
  readVariable(interp, "ModelProperty", "MODEL_DIRECTORY,absolute", dir1_abs);

  readVariable(interp, "ModelProperty", "INCLUDE_PATH", dir2);
  readVariable(interp, "ModelProperty", "INCLUDE_PATH,absolute", dir2_abs);
  readVariable(interp, "ModelProperty", "INCLUDE_PATH,model,save", dir2_save);

  readVariable(interp, "ModelProperty", "RESULTS_DIRECTORY", dir3);
  readVariable(interp, "ModelProperty", "RESULTS_DIRECTORY,absolute", dir3_abs);
  readVariable(interp, "ModelProperty", "RESULTS_DIRECTORY,model,save", dir3_save);

  readVariable(interp, "ModelProperty", "LOG_DIRECTORY", dir4);
  readVariable(interp, "ModelProperty", "LOG_DIRECTORY,absolute", dir4_abs);
  readVariable(interp, "ModelProperty", "LOG_DIRECTORY,model,save", dir4_save);

  model->setModelFileDirectories(dir1, dir2, dir3, dir4);
  model->setModelFileDirectoriesAbs(dir1_abs, dir2_abs, dir3_abs, dir4_abs);
  model->setModelFileDirectoriesSave(dir2_save, dir3_save, dir4_save);


  //--Delete buffers
  delete[] name1; delete[] name2;
  delete[] desc1; delete[] desc2;
  delete[] dir1; delete[] dir2; delete[] dir3; delete[] dir4;
  delete[] dir1_abs; delete[] dir2_abs; delete[] dir3_abs; delete[] dir4_abs;

  return TCL_OK;
}


// Function reads processor (for parallel processing) settings from Tk.
int
UserInterface_TCL::from_tk_ReadProcessorData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  //---Get tcl-data
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadProcessorData: Command not legal. No model exists!");
    return TCL_OK;
  }

  ParallelInfo pi;
  int all_argc;
  char **all_argv;
  char* allData = getCommandArguments(interp);
  int code = Tcl_SplitList(interp, allData, &all_argc, (const char ***) &all_argv);
  if (code != TCL_OK)
    return TCL_ERROR;
  int index = 0;
  pi.nofProcessors = atoi(all_argv[index++]);
  //
  model->setParallelInfo(pi);

  return TCL_OK;
}


// Function set selection method in the renderer
int
UserInterface_TCL::from_tk_ReadSelectionTolerances(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Renderer* renderer = theControlCenter->getRenderer();

  //---Get tcl-data
  char* data = getCommandArguments(interp);


  // Read toleracnes from gui, check values and
  // store them in the model and update gui
  char* normal_tol = (char *) Tcl_GetVar2(interp, "Info", "normalTolerance", TCL_GLOBAL_ONLY);
  char* distance_tol = (char *) Tcl_GetVar2(interp, "Info", "distanceTolerance", TCL_GLOBAL_ONLY);

  // Normal tolerance is degrees (0-45)
  double normal_tolerance = atof(normal_tol);
  if (normal_tolerance < 0)
    normal_tolerance *= -1;
  if (normal_tolerance > MAX_NORMAL_TOLERANCE)
    normal_tolerance = MAX_NORMAL_TOLERANCE;

  // Distance tolerance is relative (0-1)
  double distance_tolerance = atof(distance_tol);
  if (distance_tolerance < 0)
    distance_tolerance *= -1;
  if (distance_tolerance > MAX_DISTANCE_TOLERANCE)
    distance_tolerance = MAX_DISTANCE_TOLERANCE;

  strstream normal_strm;
  normal_strm << normal_tolerance << ends;

  strstream distance_strm;
  distance_strm << distance_tolerance << ends;

  NORMAL_TOLERANCE = normal_tolerance;
  DISTANCE_TOLERANCE = distance_tolerance;

  Tcl_SetVar2(interp, "Info", "normalTolerance", normal_strm.str(), TCL_GLOBAL_ONLY);
  Tcl_SetVar2(interp, "Info", "distanceTolerance", distance_strm.str(), TCL_GLOBAL_ONLY);

  return TCL_OK;
}


// Function reads user defined model parameters data from Tk.
int
UserInterface_TCL::from_tk_ReadSimulationParameterData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadSimulationParameterData: Command not legal. No model exists!");
    return TCL_OK;
  }

  setParameterData(model, ECIF_SIMULATION_PARAMETER, interp, "SimulationParameter");

  return TCL_OK;
}


// Function reads solver parameter data from Tk.
int
UserInterface_TCL::from_tk_ReadSolverData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadSolverData: Command not legal. No model exists!");
    return TCL_OK;
  }

  setParameterData(model, ECIF_SOLVER, interp, "Solver");

  return TCL_OK;
}


// Function reads solverControl parameter data from Tk.
int
UserInterface_TCL::from_tk_ReadSolverControlData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadSolverControlData: Command not legal. No model exists!");
    return TCL_OK;
  }

  setParameterData(model, ECIF_SOLVER_CONTROL, interp, "SolverControl");

  return TCL_OK;
}


// Function reads timestep parameter settings from Tk.
int
UserInterface_TCL::from_tk_ReadTimestepData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadTimestepData: Command not legal. No model exists!");
    return TCL_OK;
  }


  setParameterData(model, ECIF_TIMESTEP, interp, "Timestep");

  to_tk_WriteStatusTimesteps(interp, *model);

  // NOTE: InitialCOndition status field depends also
  // on Timestep-data (simulation type Steady State or Transient)
  to_tk_WriteStatusInitialConditions(interp, *model);

  return TCL_OK;
}


// Function reads model data related timestamp from Tk.
int
UserInterface_TCL::from_tk_ReadTimestamp(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadTimestamp: Command not legal. No model exists!");
    return TCL_OK;
  }

  //---Get tcl-data
  char* data = getCommandArguments(interp);
  int ts_argc;
  char** ts_argv;
  int code = Tcl_SplitList(interp, data, &ts_argc, (const char ***) &ts_argv);
  if (code != TCL_OK)
    return TCL_ERROR;

  //---Set timestamp
  model->setTimestamp(ts_argv[0], ts_argv[1]);
  return TCL_OK;
}


// Function reads user settings parameter values from Tk.
int
UserInterface_TCL::from_tk_ReadUserSettingsData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadUserSettingsData: Command not legal. No model exists!");
    return TCL_OK;
  }


  setParameterData(model, ECIF_USER_SETTING, interp, "UserSetting");

  return TCL_OK;
}


// Read vertex display list
int
UserInterface_TCL::from_tk_ReadVertexDisplayData(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ReadVertexDisplayData: Command not legal. No model exists!");
    return TCL_OK;
  }

  int index = 0;
  while (true) {
    BodyElement* be = model->getVertex(index++);
    if (be==NULL) break;
    int id = be->Id();
    int display = 1;
    readIdVariable(interp, "ObjectTable", id, "dspl", display);
    if (display == 0) {
      be->setDrawMode(DM_HIDDEN);
    } else {
      be->setDrawMode(DM_NORMAL);
    }
  }

  model->refreshRenderer();

  return TCL_OK;
}


// Function for the user settings file reading
int
UserInterface_TCL::from_tk_ReadUserSettingFiles(ClientData clientData, Tcl_Interp* interp,
                                                int argc, char* argv[])
{
  // NOTE: These are needed for emf_readDataCB callback function
  //
  struct emf_ObjectData_X my_ObjectData;
  emf_ObjectData  = &my_ObjectData;

  char* files = (char *) Tcl_GetVar2(interp, "Info", "userSettingFiles", glob_flag);

  //---Break this list into an array of strings (paths)
  int all_argc;
  char** all_argv;
  int code = Tcl_SplitList(interp, files, &all_argc, (const char ***) &all_argv);

  if (code != TCL_OK) {
    return TCL_ERROR;
  }

  // Read files
  for (int i = 0; i < all_argc; i++) {

    //--Set callback function
    emf_readDataCB = readUserSettingsFileCallBack;

    //--Start parser
    int rc = emf_readData(all_argv[i], (void**)&interp, 0, NULL, 1);

    strstream strm;
    strm << "Reading settings file:  " << all_argv[i] << ends;
    writeVariable(interp, "Info", "userSettingFilesReadInfo", strm.str(), false);
  }

  return TCL_OK;
}



// Function removove Cad geometry from the model
int
UserInterface_TCL::from_tk_RemoveCadGeometry(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "RemoveCadGeometry: Command not legal. No model exists!");
    return TCL_OK;
  }

  model->removeCadGeometry();

  return TCL_OK;
}


// Function draws model when option is selected from menu.
int
UserInterface_TCL::from_tk_RendererDisplayModel(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  theControlCenter->displayModel();
  return TCL_OK;
}


// Function rotates model in the display.
int
UserInterface_TCL::from_tk_RendererRotateModel(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Renderer* renderer = theControlCenter->getRenderer();
  if (renderer == NULL)
    return TCL_OK;

  int all_argc;
  char **all_argv;
  char* allData = getCommandArguments(interp);
  int code = Tcl_SplitList(interp, allData, &all_argc, (const char ***) &all_argv);
  if (code != TCL_OK)
    return TCL_ERROR;

  int axis = atoi(all_argv[0]);
  short direction = atoi(all_argv[1]);
  renderer->rotate(axis, direction);
  return TCL_OK;
}


// Function resets model display.
int
UserInterface_TCL::from_tk_RendererResetModel(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Renderer* renderer = theControlCenter->getRenderer();

  if (renderer == NULL)
    return TCL_OK;

  renderer->reset();

  return TCL_OK;
}


// Function scales model in the display.
int
UserInterface_TCL::from_tk_RendererScaleModel(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Renderer* renderer = theControlCenter->getRenderer();
  if (renderer == NULL)
    return TCL_OK;

  int all_argc;
  char **all_argv;
  char* allData = getCommandArguments(interp);
  int code = Tcl_SplitList(interp, allData, &all_argc, (const char ***) &all_argv);
  if (code != TCL_OK)
    return TCL_ERROR;

  short direction   = atoi(all_argv[0]);
  renderer->scale(direction);

  return TCL_OK;
}


// Set boundary edit mode
int
UserInterface_TCL::from_tk_RendererSetEditBoundaries(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Renderer* renderer = theControlCenter->getRenderer();
  if (renderer == NULL)
    return TCL_OK;

  bool in_edit_mode = bool(atoi(getCommandArguments(interp)));

  renderer->setEditBoundaries(in_edit_mode);

  return TCL_OK;
}


// Set rotation priorities in the renderer
int
UserInterface_TCL::from_tk_RendererSetRotatePriorities(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Renderer* renderer = theControlCenter->getRenderer();
  if (renderer == NULL)
    return TCL_OK;

  int all_argc;
  char **all_argv;
  char* allData = getCommandArguments(interp);
  int code = Tcl_SplitList(interp, allData, &all_argc,(const char ***) &all_argv);
  if (code != TCL_OK)
    return TCL_ERROR;

  int x = atoi(all_argv[0]);
  int y = atoi(all_argv[1]);
  int z = atoi(all_argv[2]);

  renderer->setRotatePriorities(x, y, z);

  return TCL_OK;
}


// Function tranlates model in the display.
int
UserInterface_TCL::from_tk_RendererTranslateModel(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Renderer* renderer = theControlCenter->getRenderer();
  if (renderer == NULL)
    return TCL_OK;

  int all_argc;
  char **all_argv;
  char* allData = getCommandArguments(interp);
  int code = Tcl_SplitList(interp, allData, &all_argc, (const char ***) &all_argv);
  if (code != TCL_OK)
    return TCL_ERROR;

  int coordinate  = atoi(all_argv[0]);
  short direction   = atoi(all_argv[1]);
  renderer->translate(coordinate, direction);

  return TCL_OK;
}


// Reset all boundary selctions
int
UserInterface_TCL::from_tk_ResetAllBoundarySelections(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ResetAllBoundarySelections: Command not legal. No model exists!");
    return TCL_OK;
  }

  model->resetAllBoundarySelections(true);

  return TCL_OK;
}


// Reset boundary selctions when a new boundary is selected
int
UserInterface_TCL::from_tk_ResetBoundarySelections(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ResetBoundarySelections: Command not legal. No model exists!");
    return TCL_OK;
  }

  model->resetBoundarySelections(true, false, 0, (const int*)NULL, false);

  return TCL_OK;
}


// Function sets original names for boundaries
int
UserInterface_TCL::from_tk_RestoreBoundaryNames(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "RestoreBoundaryNames: Command not legal. No model exists!");
    return TCL_OK;
  }

  model->restoreBoundaryNames();

  TheUI->updateBoundaryData(model);

  return TCL_OK;
}


// Function saves external mesh in Elmer (DB) format
int
UserInterface_TCL::from_tk_SaveElmerMeshFile(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{

  char* mesh_dir = getCommandArguments(interp);
  theControlCenter->saveElmerMeshFile(mesh_dir);
  return TCL_OK;
}


// Function saves mesh result file in ElmerPost format
// (ie. activates its construction and saving by the model)
int
UserInterface_TCL::from_tk_SaveElmerPostMeshFile(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  char* out_filename = getCommandArguments(interp);
  theControlCenter->saveElmerPostMeshFile(out_filename);
  return TCL_OK;
}


// Function saves model file (emf-file)
int
UserInterface_TCL::from_tk_SaveModelFile(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "SaveModelFile: Command not legal. No model exists!");
    return TCL_OK;
  }

  // Check if model file matc definitions (dependencies) should be kept
  int drop_matc;
  readVariable(interp, "Info", "dropModelMatcDefinitions", drop_matc);

  if ( 0 == drop_matc ) {
    model->setKeepMatcDefinitions(true);
  } else {
    model->setKeepMatcDefinitions(false);
  }

  char* emf_filename = getCommandArguments(interp);
  theControlCenter->saveFrontModelFile(emf_filename);

  return TCL_OK;
}


// Function saves mesh input file (mif-file)
int
UserInterface_TCL::from_tk_SaveMeshInputFile(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  char* mif_filename = getCommandArguments(interp);
  theControlCenter->saveMeshInputFile(mif_filename);
  return TCL_OK;
}


// Function saves solver input file.(sif-file)
int
UserInterface_TCL::from_tk_SaveSolverInputFile(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  char* sif_filename = getCommandArguments(interp);
  theControlCenter->saveSolverInputFile(sif_filename);
  return TCL_OK;
}


// Function saves mesh result file in Thetis format
// (ie. activates its construction and saving by the model)
int
UserInterface_TCL::from_tk_SaveThetisMeshFile(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  char* out_filename = getCommandArguments(interp);
  theControlCenter->saveThetisMeshFile(out_filename);
  return TCL_OK;
}


// Function saves user settings into (default) file
int
UserInterface_TCL::from_tk_SaveUserSettingsFile(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  char* filename = getCommandArguments(interp);
  theControlCenter->saveUserSettingsFile(filename);
  return TCL_OK;
}


// Function saves model ecf-file.
int
UserInterface_TCL::from_tk_SelectMeshBoundaryElements(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "SelectMeshBoundaryElements: Command not legal. No model exists!");
    return TCL_OK;
  }

  char* type = getCommandArguments(interp);

  if ( 0 == strcmp(type, "do") )
    model->selectMeshBoundaryElements();
  else if ( 0 == strcmp(type, "undo") )
    model->selectMeshBoundaryElementsUndo();
  else if ( 0 == strcmp(type, "redo") )
    model->selectMeshBoundaryElementsRedo();

  return TCL_OK;
}


// Change current working diredctory
int
UserInterface_TCL::from_tk_SetCurrentDirectory(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  char* dir = getCommandArguments(interp);

#if defined(WIN32)
  _chdir(dir);
#else
  chdir(dir);
#endif

  return TCL_OK;
}


// Function sets a (predefined) flag in the model
// In Tk flags are named ModelFlags(drawBodies) (value 0/1)
// Flags are stored in the model in the ModelFlags array where
// indices are enums like DRAW_TARGET_BODIES
// Group name enum for this flag would be DRAW_TARGET
// Tk identifies the group and flag names for the cpp by
// using corresponding strings "DRAW_TARGET" and "DRAW_TARGET_BODIES"
int
UserInterface_TCL::from_tk_SetFlagValue(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "SetFlagValue: Command not legal. No model exists!");
    return TCL_OK;
  }

  //---Get tcl-data
  int all_argc;
  char **all_argv;
  char* allData = getCommandArguments(interp);
  int code = Tcl_SplitList(interp, allData, &all_argc, (const char ***) &all_argv);
  if (code != TCL_OK)
    return TCL_ERROR;

  bool data_ok = true;

  //---First argument: flag group
  flagGroup group;
  if ( 0 == strcmp(all_argv[0], "DRAW_SOURCE") )
    group = DRAW_SOURCE;

  else if ( 0 == strcmp(all_argv[0], "DRAW_TARGET") )
    group = DRAW_TARGET;

  else if ( 0 == strcmp(all_argv[0], "GEOMETRY_TYPE") )
    group = GEOMETRY_TYPE;

  else if ( 0 == strcmp(all_argv[0], "SELECT_METHOD") )
    group = SELECT_METHOD;

  else if ( 0 == strcmp(all_argv[0], "SELECT_MODE") )
    group = SELECT_MODE;

  else if ( 0 == strcmp(all_argv[0], "SELECT_OBJECTS") )
    group = SELECT_OBJECTS;

  else if ( 0 == strcmp(all_argv[0], "LABEL_DISPLAY") )
    group = LABEL_DISPLAY;

  else {
    group = FLAG_GROUP_UNKNOWN;
    data_ok = false;
    strstream strm;
    strm << "Unknown flag group name: " << all_argv[0] << ends;
    TheUI->showMsg(strm.str());
  }

  //---Second argument: flag name
  flagName name;
  if ( 0 == strcmp(all_argv[1], "DRAW_SOURCE_CAD") )
    name = DRAW_SOURCE_CAD;

  else if ( 0 == strcmp(all_argv[1], "DRAW_SOURCE_MESH") )
    name = DRAW_SOURCE_MESH;

  else if ( 0 == strcmp(all_argv[1], "DRAW_TARGET_BODIES") )
    name = DRAW_TARGET_BODIES;

  else if ( 0 == strcmp(all_argv[1], "DRAW_TARGET_SURFACES") )
    name = DRAW_TARGET_SURFACES;

  else if ( 0 == strcmp(all_argv[1], "DRAW_TARGET_EDGES") )
    name = DRAW_TARGET_EDGES;

  else if ( 0 == strcmp(all_argv[1], "SELECT_METHOD_SINGLE") )
    name = SELECT_METHOD_SINGLE;

  else if ( 0 == strcmp(all_argv[1], "SELECT_METHOD_ALL") )
    name = SELECT_METHOD_ALL;

  else if ( 0 == strcmp(all_argv[1], "SELECT_METHOD_BY_NEIGHBOR") )
    name = SELECT_METHOD_BY_NEIGHBOR;

  else if ( 0 == strcmp(all_argv[1], "SELECT_METHOD_BY_NORMAL") )
    name = SELECT_METHOD_BY_NORMAL;

  else if ( 0 == strcmp(all_argv[1], "SELECT_METHOD_BY_PLANE") )
    name = SELECT_METHOD_BY_PLANE;

  else if ( 0 == strcmp(all_argv[1], "SELECT_METHOD_BY_BOX") )
    name = SELECT_METHOD_BY_BOX;

  else if ( 0 == strcmp(all_argv[1], "SELECT_METHOD_BY_RECTANGLE") )
    name = SELECT_METHOD_BY_RECTANGLE;

  else if ( 0 == strcmp(all_argv[1], "SELECT_MODE_TOGGLE") )
    name = SELECT_MODE_TOGGLE;

  else if ( 0 == strcmp(all_argv[1], "SELECT_MODE_EXTEND") )
    name = SELECT_MODE_EXTEND;

  else if ( 0 == strcmp(all_argv[1], "SELECT_MODE_REDUCE") )
    name = SELECT_MODE_REDUCE;

  else if ( 0 == strcmp(all_argv[1], "SELECT_OBJECTS_EXTEND") )
    name = SELECT_OBJECTS_EXTEND;

  else if ( 0 == strcmp(all_argv[1], "SELECT_OBJECTS_TOGGLE") )
    name = SELECT_OBJECTS_TOGGLE;

  else if ( 0 == strcmp(all_argv[1], "LABEL_DISPLAY_NODE") )
    name = LABEL_DISPLAY_NODE;

  else if ( 0 == strcmp(all_argv[1], "LABEL_DISPLAY_ELEMENT") )
    name = LABEL_DISPLAY_ELEMENT;

  else if ( 0 == strcmp(all_argv[1], "LABEL_DISPLAY_VERTEX") )
    name = LABEL_DISPLAY_VERTEX;

  else if ( 0 == strcmp(all_argv[1], "LABEL_DISPLAY_EDGE") )
    name = LABEL_DISPLAY_EDGE;

  else if ( 0 == strcmp(all_argv[1], "LABEL_DISPLAY_FACE") )
    name = LABEL_DISPLAY_FACE;

  else if ( 0 == strcmp(all_argv[1], "LABEL_DISPLAY_BODY") )
    name = LABEL_DISPLAY_BODY;

  else {
    name = FLAG_NAME_UNKNOWN;
    data_ok = false;
    strstream strm;
    strm << "Unknown flag name: " << all_argv[1] << ends;
    TheUI->showMsg(strm.str());
  }

  if (data_ok) {
    bool value = bool(atoi(all_argv[2]));
    model->setFlagValue(group, name, value);
  }

  return TCL_OK;
}


// Set Matc emf-input file name
int
UserInterface_TCL::from_tk_SetMatcInputFileEmf(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "SetMatcInputFileEmf: Command not legal. No model exists!");
    return TCL_OK;
  }

  char* file = getCommandArguments(interp);

  model->setMatcInputFileEmf(file);

  return TCL_OK;
}


// Set Matc sif-input file name
int
UserInterface_TCL::from_tk_SetMatcInputFileSif(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "SetMatcInputFileSif: Command not legal. No model exists!");
    return TCL_OK;
  }

  char* file = getCommandArguments(interp);

  model->setMatcInputFileSif(file);

  return TCL_OK;
}


// Set mesh input unit (for external mesh reading)
// NOTE: uses static model method!
//
int
UserInterface_TCL::from_tk_SetMeshInputUnit(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  char* value = getCommandArguments(interp);

  double unit = atof(value);

  if ( unit > 0 ) {
    Model::setMeshInputUnit(unit);
  } else {
    Model::setMeshInputUnit(1.0);
  }

  return TCL_OK;
}


// Function set model status data from Tk.
int
UserInterface_TCL::from_tk_SetModelStatus(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "SetModelStatus: Command not legal. No model exists!");
    return TCL_OK;
  }

  //---Get tcl-data
  char* data = getCommandArguments(interp);
  int st_argc;
  char** st_argv;
  int code = Tcl_SplitList(interp, data, &st_argc, (const char ***) &st_argv);
  if (code != TCL_OK)
    return TCL_ERROR;
  //---Set model status
  ecif_modelStatus status = (ecif_modelStatus)atol(st_argv[0]);
  model->setModelStatus(status);
  return TCL_OK;
}



// Function reads model outfile save timestamp from Tk.
int
UserInterface_TCL::from_tk_SetSelectionsToGui(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "SetSelectionsToGui: Command not legal. No model exists!");
    return TCL_OK;
  }

  model->setSelectionsToGui();

  return TCL_OK;
}



// Function starts boundary splitting
int
UserInterface_TCL::from_tk_SplitBoundary(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "SplitBoundary: Command not legal. No model exists!");
    return TCL_OK;
  }

  int body1_id, body2_id;

  istrstream in(getCommandArguments(interp));
  in >> body1_id >> body2_id;

  model->splitBoundary(body1_id, body2_id);
  return TCL_OK;
}


// Function redes boundary splitting/combining
int
UserInterface_TCL::from_tk_SplitCombineBoundariesRedo(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "SplitCombineBoundariesRedo: Command not legal. No model exists!");
    return TCL_OK;
  }

  model->splitCombineBoundariesRedo();

  return TCL_OK;
}


// Function undoes boundary splitting/combining
int
UserInterface_TCL::from_tk_SplitCombineBoundariesUndo(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "SplitCombineBoundariesUndo: Command not legal. No model exists!");
    return TCL_OK;
  }

  model->splitCombineBoundariesUndo();

  return TCL_OK;
}


// Function send the stop editing message to the model
int
UserInterface_TCL::from_tk_StopEditMeshBoundaries(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "StopEditMeshBoundaries: Command not legal. No model exists!");
    return TCL_OK;
  }

  char* mode = getCommandArguments(interp);

  bool cancel_edit = false;

  if ( LibFront::ncEqual(mode, "cancel") )
    cancel_edit = true;

  model->stopEditMeshBoundaries(cancel_edit);

  return TCL_OK;
}


// Function accepts current names for boundaries
int
UserInterface_TCL::from_tk_StoreBoundaryNames(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "StoreBoundaryNames: Command not legal. No model exists!");
    return TCL_OK;
  }

  model->storeBoundaryNames();

  return TCL_OK;
}


// Function initiates interrupts processing by setting a proper stop flag on
int
UserInterface_TCL::from_tk_DoBreak(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  theControlCenter->setBreakValue(MESH_INPUT, true);

  return TCL_OK;
}


// Delete mesh from the model data
int
UserInterface_TCL::from_tk_UnloadMesh(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "UnloadMesh: Command not legal. No model exists!");
    return TCL_OK;
  }

  char* msg = getCommandArguments(interp);

  model->unloadMesh(msg);

  return TCL_OK;
}


// Update CAD geometry (based on Matc parameters)
int
UserInterface_TCL::from_tk_UpdateCadGeometry(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "UpdateCadGeometry: Command not legal. No model exists!");
    return TCL_OK;
  }

  model->updateCadGeometry();

  return TCL_OK;
}


// Function updates a matc file with user selcted definitions
int
UserInterface_TCL::from_tk_UpdateMatcFile(ClientData clientData, Tcl_Interp *interp,
                                          int argc, char* argv[])
{
  char* ud_data = getCommandArguments(interp);
  int ud_argc;
  char** ud_argv;
  int code = Tcl_SplitList(interp, ud_data, &ud_argc, (const char ***) &ud_argv);

  if (code != TCL_OK)
    return TCL_ERROR;

  char* filename = ud_argv[0];
  char* mode = ud_argv[1];

  char* mc_data = ud_argv[2];
  int mc_argc;
  char** mc_argv;
  code = Tcl_SplitList(interp, mc_data, &mc_argc, (const char ***) &mc_argv);

  Model::updateMatcFile(filename, mode, mc_argc, mc_argv);

  return TCL_OK;
}

//-----------------------------------

void
UserInterface_TCL::generateEvent()
{
   Tcl_Eval(theInterp, "Util::generateEvent");
}


char*
UserInterface_TCL::getCommandArguments(Tcl_Interp* interp)
{
  return (char *) Tcl_GetVar2(interp, "Info", "arguments", glob_flag);
}


char*
UserInterface_TCL::getCommandResults(Tcl_Interp* interp)
{
  return (char *) Tcl_GetVar2(interp, "Info", "results", glob_flag);
}


Renderer*
UserInterface_TCL::getRenderer()
{
  return theControlCenter->getRenderer();
}


void
UserInterface_TCL::getCurrentTimestamp(char* buffer)
{
  sendCommandToGui(theInterp, "Interface::setCurrentTimestamp");

  char* ts = (char *) Tcl_GetVar2(theInterp, "Info", "currentTimestamp", glob_flag);

  strcpy(buffer, ts);
}


bool
UserInterface_TCL::getEquationVarsVariable(const char* equation_name, char*& equation_vars_name)
{
  // Panel names to check (Equation and generic)
  char* panel_names[] = { "EQ", "#" };
  int panel_count = 2;

  for (int i = 0; i < panel_count; i++) {

    strstream strm;
    strm << panel_names[i] << "," << equation_name << "," << "EquationVars" << ends;
    readVariable(theInterp, "Common", strm.str(), equation_vars_name);

    if ( equation_vars_name != NULL ) {
      return true;
    }
  }

  return false;
}


// Check if an equation panel field (possibly for a specific equation) should
// be output into a solver section
//
bool
UserInterface_TCL::getIsSolverTargetField(const char* equation_name, const char* field_name)
{
  if ( field_name == NULL ) return false;

  ostrstream strm;

  if ( equation_name != NULL ) {
    strm << equation_name << tclArgSeparator;
  } else {
    strm << "" << tclArgSeparator;
  }

  strm << field_name << ends;

  sendCommandToGui(theInterp, "Interface::getIsSolverTargetField", strm.str());

  char* data = getCommandResults(theInterp);
  int tf_argc;
  char** tf_argv;
  int code = Tcl_SplitList(theInterp, data, &tf_argc, (const char ***) &tf_argv);

  if ( tf_argc == 0 ) {
    return false;
  }

  return (bool)atoi(tf_argv[0]);
}


void
UserInterface_TCL::getMatcSifDefinitions(int& nof_defs, char**& defs)
{
  readVariable(theInterp, "MatcDefinitions", "sifDefs", nof_defs, defs);
}


void
UserInterface_TCL::getMeshDirectoryInfo(char*& dir, char*& dir_abs)
{
  readVariable(theInterp, "ModelProperty", "CURRENT_MESH_DIRECTORY", dir);
  readVariable(theInterp, "ModelProperty", "CURRENT_MESH_DIRECTORY,absolute", dir_abs);
}


bool
UserInterface_TCL::getMeshInputFileName(char*& mif_file_name)
{
  sendCommandToGui(theInterp, "Interface::getMeshInputFileName");

  char* mif_fn = getCommandResults(theInterp);

  if ( mif_fn == NULL ) {
    return false;
  }

  mif_file_name = new char[1 + strlen(mif_fn)];

  strcpy(mif_file_name, mif_fn);

  return true;
}


void
UserInterface_TCL::getModelDirectoryInfo(char*& dir, char*& dir_abs)
{
  readVariable(theInterp, "ModelProperty", "MODEL_DIRECTORY", dir);
  readVariable(theInterp, "ModelProperty", "MODEL_DIRECTORY,absolute", dir_abs);
}


void
UserInterface_TCL::getModelNameInfo(char*& model_name, char*& problem_name)
{
  readVariable(theInterp, "ModelProperty", "MODEL_NAME", model_name);
  readVariable(theInterp, "ModelProperty", "PROBLEM_NAME", problem_name);
}


bool
UserInterface_TCL::getParameterFieldInfo(const char* parameter, const char* field, ParameterFieldInfo& finfo)
{
  ostrstream strm;

  strm << parameter << tclArgSeparator << field << ends;

  sendCommandToGui(theInterp, "Interface::getParameterFieldInfo", strm.str());

  char* data = getCommandResults(theInterp);
  int fi_argc;
  char** fi_argv;
  int code = Tcl_SplitList(theInterp, data, &fi_argc, (const char ***) &fi_argv);

  if ( fi_argc == 0 ) {
    return false;
  }

  // NOTE: data order is fixed and data given by Gui should match this!!!
  int idx = 0;
  update_dyna_string(finfo.sifName, fi_argv[idx++]);
  update_dyna_string(finfo.valueType, fi_argv[idx++]);

  finfo.outputSif = (bool)atoi(fi_argv[idx++]);
  finfo.outputSifType  = (bool)atoi(fi_argv[idx++]);
  finfo.alwaysOutput = (bool)atoi(fi_argv[idx++]);
  finfo.isArray  = (bool)atoi(fi_argv[idx++]);
  finfo.isQuoted = (bool)atoi(fi_argv[idx++]);

  // Note: these flags are by default turned off in the data
  finfo.isFileName = (bool)atoi(fi_argv[idx++]);
  finfo.isProcName = (bool)atoi(fi_argv[idx++]);

  return true;
}


bool
UserInterface_TCL::getSolverKeywordTypeGiven(const char* parameter, const char* field)
{
  ostrstream strm;

  strm << parameter << tclArgSeparator << field << ends;

  sendCommandToGui(theInterp, "Interface::getSolverKeywordTypeGiven", strm.str());

  char* data = getCommandResults(theInterp);
  int fi_argc;
  char** fi_argv;
  int code = Tcl_SplitList(theInterp, data, &fi_argc, (const char ***) &fi_argv);

  if ( fi_argc == 0 ) {
    return false;
  }

  return (bool)atoi(fi_argv[0]);
}


bool
UserInterface_TCL::getUseModelFileSettings()
{
  int value;
  readVariable(theInterp, "UserSetting", "DEFAULT_USE_MODEL_SETTINGS", value);

  if (value == 0)
    return false;
  else
    return true;
}


// Check if variable names should be added to equation names in
// Solver sif output
bool
UserInterface_TCL::getUseVariableNameInEquationName(const char* equation_name)
{
  sendCommandToGui(theInterp, "Interface::getUseVariableNameInEquationName", equation_name);

  bool use_var_name = atoi(getCommandResults(theInterp));

  return use_var_name;
}


// *** Set values for Tcl-variables.
void
UserInterface_TCL::initTclVariables(Tcl_Interp* interp, const Model& model)
{
  to_tk_WriteModelStatus(interp, model);
  to_tk_WriteModelGeometryDimension(interp, model);
  to_tk_WriteProcessorData(interp, model);

  to_tk_WriteBodyLayerData(interp, model);
  to_tk_WriteBodyData(interp, model);
  to_tk_WriteBodyInfoData(interp, model);

  // Do this after writing Body level data!!!
  to_tk_WriteBoundaryData(interp, model);
  to_tk_WriteElementGroupData(interp, model);

  to_tk_WriteStats(interp, model);
  to_tk_WriteControlParameters(interp, model);

  // Process object table data in Gui side
  sendCommandToGui(interp, "Interface::applyObjectTableData");

  to_tk_WriteAllParamDataPre(interp, model);
  to_tk_WriteAllParamDataPost(interp, model);
}


// Function makes two "separated" strings from object, elementd ids data
// Into body_ids: (Body-pair-ids, Nof-elements)
// Into elem_ids: (Body1-id, Body2-id, Element-id)
// Into bndr_ids: (Boundary condition id)
// Into grid_ids:  (Grid-param id)
// id_set: the id-data, a set of ids-structures
void
UserInterface_TCL::listInnerIds(Model* model,
              ostrstream& body_ids, ostrstream& elem_ids,
              ostrstream& names,
              ostrstream& bndr_ids,
              const char obj_sep, const char fld_sep,
              Ids3Set& id_set)
{
  Ids3Set::iterator pos = id_set.begin();
  Ids3Set::iterator end = id_set.end();
  int elm_counter;
  int object1_id, object2_id;
  object1_id = -1;

  // We collect unique ids of those bodies which have an
  // inner boundary condition defined into this set
  while (pos != end) {
    Ids3 wid3 = *pos++;

    //New object-pair was encountered
    if (object1_id != wid3.id1 || object2_id != wid3.id2) {

      //Write data after each 'real' new object
      if (object1_id != -1) {
        body_ids << 'I' << fld_sep;
        body_ids << object1_id << fld_sep;
        body_ids << object2_id << fld_sep;
        body_ids << elm_counter;
        body_ids << obj_sep;
      }
      //Update reference variables
      object1_id = wid3.id1;
      object2_id = wid3.id2;
      elm_counter = 0;
    }

    elm_counter++;

    //Write out element-row, boundary condtition and grid-param rows, Body2, Elm, Cond
    BodyElement* be = model->getBoundaryById(wid3.id3);
    int cid = be->getBoundaryConditionId();
    elem_ids << 'I' << fld_sep;
    elem_ids << wid3.id1 << fld_sep << wid3.id2 << fld_sep << wid3.id3;
    names << be->getName();
    bndr_ids << cid;

    //We don't want to add a separator at the end of the string!!!
    if (pos != end) {
      elem_ids << obj_sep;
      names << obj_sep;
      bndr_ids << obj_sep;
    }
  }

  //End of data. Write last object's data.
  if (object1_id != -1) {
    body_ids << 'I' << fld_sep;
    body_ids << object1_id << fld_sep;
    body_ids << object2_id << fld_sep;
    body_ids << elm_counter;
  }

  body_ids << ends;
  elem_ids << ends;
  names << ends;
  bndr_ids << ends;
}


// Function makes two "separated" strings from object, elementd ids data
// Into body_ids: (Body-id, Nof-elements)
// Into elem_ids: (Body-id, Element-id)
// Into bndr_ids: (Boundary condition id)
// id_set: the id-data, a set of ids-structures
void
UserInterface_TCL::listOuterIds(Model* model,
              ostrstream& body_ids, ostrstream& elem_ids,
              ostrstream& names,
              ostrstream& bndr_ids,
              const char obj_sep, const char fld_sep,
              Ids2Set& id_set)
{
  Ids2Set::iterator pos = id_set.begin();
  Ids2Set::iterator end = id_set.end();
  int elm_counter;
  int object_id = -1;
   // We collect unique ids of those bodies which have an
   // outer boundary condition defined into this set
  while (pos != end) {
    Ids2 wid2 = *pos++;

    //New object-pair was encountered
    if (object_id != wid2.id1) {
      //Write data after each 'real' new object
      if (object_id != -1) {
        body_ids << 'O' << fld_sep;
        body_ids << object_id << fld_sep;
        body_ids << elm_counter;
        body_ids << obj_sep;
      }
      //Update reference variables
      object_id = wid2.id1;
      elm_counter = 0;
    }

    elm_counter++;

    //Write out element, boundary condition and grid-param row. Body, Elm, Cond
    BodyElement* be = model->getBoundaryById(wid2.id2);
    int cid = be->getBoundaryConditionId();
    elem_ids << 'O' << fld_sep;
    elem_ids << wid2.id1 << fld_sep << wid2.id2;
    names << be->getName();
    bndr_ids << cid;

    //We don't want to add a separator at the end of the string!!!
    if (pos != end) {
      elem_ids << obj_sep;
      names << obj_sep;
      bndr_ids << obj_sep;
    }
  }

  //End of data. Write last object's data.
  if (object_id != -1) {
    body_ids << 'O' << fld_sep;
    body_ids << object_id << fld_sep;
    body_ids << elm_counter;
  }

  body_ids << ends;
  elem_ids << ends;
  names << ends;
  bndr_ids << ends;
}


void
UserInterface_TCL::markSelectedBoundaries()
{
  sendCommandToGui(theInterp, "Interface::markSelectedBoundaries", NULL);
}


void
UserInterface_TCL::matcFileWasRead(char* filename)
{
  sendCommandToGui(theInterp, "Interface::matcFileWasRead", filename);
}


// Main function for the user settings file reading
void
UserInterface_TCL::readUserSettingsFile(Tcl_Interp* interp, char* filename)
{
  struct emf_ObjectData_X my_ObjectData;
  emf_ObjectData  = &my_ObjectData;

  /*Set callback function*/
  //int my_readDataCB(void**);
  emf_readDataCB = readUserSettingsFileCallBack;

  /*-----Start parser*/
  int rc = emf_readData(filename, (void**)&interp, 0, NULL, 1);

  // No success
  if (rc == 0) {
    return;
  }

  strstream strm;
  strm << "Reading settings file: " << filename << ends;
  writeVariable(interp, "Info", "settingsFileReadInfo", strm.str(), false);

}


// This method gets one "field" at a time from the file
int
UserInterface_TCL::readUserSettingsFileCallBack(void** user_data)
{
  UserInterface* gui = theControlCenter->getGui();

  Tcl_Interp* interp = (Tcl_Interp*)*user_data;

  // Global object pointing to the record, make
  // a short hand for it
  emf_ObjectData_X* od = emf_ObjectData;

  const char* on = od->object_name;
  int onl = od->object_name_length;
  const char* fn = od->field_name;
  int fnl = od->field_name_length;
  const char* dt = od->data_type;
  int dtl = od->data_type_length;


  // Read data according to the object type
  char fld[256];
  bool logical_value;
  bool to_upper = true;

  char* str_val;

  //-Default directory for Elmer models
  if ( LibFront::in("Default Model Directory", fn) ) {
    writeVariable(interp, "UserSetting", "DEFAULT_MODEL_DIRECTORY", (char*)od->data);
    return 1;
  }

  //-Default directory for Elmer Cad files
  if ( LibFront::in("Default Cad Files Directory", fn) ) {
    writeVariable(interp, "UserSetting", "DEFAULT_CAD_FILES_DIRECTORY", (char*)od->data);
    return 1;
  }

  //-Default directory for Elmer external mesh files
  if ( LibFront::in("Default External Mesh Files Directory", fn) ) {
    writeVariable(interp, "UserSetting", "DEFAULT_EXTERNAL_MESH_FILES_DIRECTORY", (char*)od->data);
    return 1;
  }

  //-Default directory for Elmer input files files
  if ( LibFront::in("Default Include Path", fn) ) {
    writeVariable(interp, "UserSetting", "DEFAULT_INCLUDE_PATH", (char*)od->data);
    return 1;
  }

  //-Results directory
  if ( LibFront::in("Default Results Directory", fn) ) {
    writeVariable(interp, "UserSetting", "DEFAULT_RESULTS_DIRECTORY", (char*)od->data);
    return 1;
  }

  //-Default directory for Elmer output files
  if ( LibFront::in("Default Log Directory", fn) ) {
    writeVariable(interp, "UserSetting", "DEFAULT_LOG_DIRECTORY", (char*)od->data);
    return 1;
  }

  //-Use model settings flag
  if ( LibFront::in("Default Use Model Settings", fn) ) {
    LibFront::setNumericData(logical_value, 0);
    writeVariable(interp, "UserSetting", "DEFAULT_USE_MODEL_SETTINGS", logical_value);
    return 1;
  }

  //-Auto load mesh file flag
  if ( LibFront::in("Auto Load Mesh", fn) ) {
    LibFront::setNumericData(logical_value, 0);
    writeVariable(interp, "UserSetting", "AUTO_LOAD_MESH", logical_value);
    return 1;
  }

  //-Auto save model file flag
  if ( LibFront::in("Auto Save Model", fn) ) {
    LibFront::setNumericData(logical_value, 0);
    writeVariable(interp, "UserSetting", "AUTO_SAVE_MODEL", logical_value);
    return 1;
  }

  //-Auto save solver input file flag
  if ( LibFront::in("Auto Save Solver Input", fn) ) {
    LibFront::setNumericData(logical_value, 0);
    writeVariable(interp, "UserSetting", "AUTO_SAVE_SOLVER_INPUT", logical_value);
    return 1;
  }

  //-External browser program command
  if ( LibFront::in("Browser Command", fn) ) {
    LibFront::setNumericData(logical_value, 0);
    writeVariable(interp, "UserSetting", "BROWSER_COMMAND", (char*)od->data);
    return 1;
  }

  //-External editor program command
  if ( LibFront::in("Editor Command", fn) ) {
    LibFront::setNumericData(logical_value, 0);
    writeVariable(interp, "UserSetting", "EDITOR_COMMAND", (char*)od->data);
    return 1;
  }

  //-Gebhardt factors run browse mode
  if ( LibFront::in("Browse Mode Gebhardt Ffactors", fn) ) {
    writeVariable(interp, "UserSetting", "BROWSE_MODE_GEBHARDT_FACTORS", (char*)od->data);
    return 1;
  }

  //-Mesh run browse mode
  if ( LibFront::in("Browse Mode Mesh", fn) ) {
    writeVariable(interp, "UserSetting", "BROWSE_MODE_MESH", (char*)od->data);
    return 1;
  }

  //-Procedure compiler run browse mode
  if ( LibFront::in("Browse Mode Procedure Compiler", fn) ) {
    writeVariable(interp, "UserSetting", "BROWSE_MODE_PROCEDURE_COMPILER", (char*)od->data);
    return 1;
  }

  //-Solver run browse mode
  if ( LibFront::in("Browse Mode Solver", fn) ) {
    writeVariable(interp, "UserSetting", "BROWSE_MODE_SOLVER", (char*)od->data);
    return 1;
  }

  //-View factors run browse mode
  if ( LibFront::in("Browse Mode View Factors", fn) ) {
    writeVariable(interp, "UserSetting", "BROWSE_MODE_VIEW_FACTORS", (char*)od->data);
    return 1;
  }

  //-Font sizes (small, medium, large)
  if ( LibFront::in("Font Sizes", fn) ) {

    int count = od->data_length;
    int* sizes = new int[count];
    for (int i = 0; i < count; i++)
      LibFront::setNumericData(sizes[i], i);

    writeVariable(interp, "UserSetting", "FONT_SIZES", count, sizes);

    delete[] sizes;
    return 1;
  }


  //-UNKNOWN field, collect info in the gui-variable to
  // be displayed later when all the nice gui machinery
  // is loaded!
  else {
    strstream strm;
    strm << "Unknown keyword: "
         << od->field_name
         << ends;

    writeVariable(interp, "Info", "useSettingFilesErrorInfo", strm.str(), false);
  }

  return 1;
}


// ========
// Id Array
// ========

// Id array, integer id, single value
// ----------------------------------
void
UserInterface_TCL::readIdVariable(Tcl_Interp* interp, const char* array, int id, const char* variable,
                               char& value)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  readVariable_impl_1(interp, array, field_name.str(), value);
}

void
UserInterface_TCL::readIdVariable(Tcl_Interp* interp, const char* array, int id, const char* variable,
                               char*& value)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  readVariable_impl_1(interp, array, field_name.str(), value);
}

void
UserInterface_TCL::readIdVariable(Tcl_Interp* interp, const char* array, int id, const char* variable,
                               int& value)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  readVariable_impl_1(interp, array, field_name.str(), value);
}

void
UserInterface_TCL::readIdVariable(Tcl_Interp* interp, const char* array, int id, const char* variable,
                               double& value)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  readVariable_impl_1(interp, array, field_name.str(), value);
}


// Id array, integer id, list of values
// ------------------------------------
void
UserInterface_TCL::readIdVariable(Tcl_Interp* interp, const char* array, int id, const char* variable,
                               int& size, char**& values)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  readVariable_impl_n(interp, array, field_name.str(), size, values);
}

void
UserInterface_TCL::readIdVariable(Tcl_Interp* interp, const char* array, int id, const char* variable,
                               int& size, int*& values)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  readVariable_impl_n(interp, array, field_name.str(), size, values);
}

void
UserInterface_TCL::readIdVariable(Tcl_Interp* interp, const char* array, int id, const char* variable,
                               int& size, double*& values)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  readVariable_impl_n(interp, array, field_name.str(), size, values);
}


// Id array, string id, single value
// ---------------------------------
void
UserInterface_TCL::readIdVariable(Tcl_Interp* interp, const char* array, const char* id, const char* variable,
                               char& value)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  readVariable_impl_1(interp, array, field_name.str(), value);
}

void
UserInterface_TCL::readIdVariable(Tcl_Interp* interp, const char* array, const char* id, const char* variable,
                               char*& value)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  readVariable_impl_1(interp, array, field_name.str(), value);
}

void
UserInterface_TCL::readIdVariable(Tcl_Interp* interp, const char* array, const char* id, const char* variable,
                               int& value)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  readVariable_impl_1(interp, array, field_name.str(), value);
}

void
UserInterface_TCL::readIdVariable(Tcl_Interp* interp, const char* array, const char* id, const char* variable,
                               double& value)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  readVariable_impl_1(interp, array, field_name.str(), value);
}


// Id array, string id, list of values
// -----------------------------------
void
UserInterface_TCL::readIdVariable(Tcl_Interp* interp, const char* array, const char* id, const char* variable,
                               int& size, char**& values)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  readVariable_impl_n(interp, array, field_name.str(), size, values);
}

void
UserInterface_TCL::readIdVariable(Tcl_Interp* interp, const char* array, const char* id, const char* variable,
                               int& size, int*& values)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  readVariable_impl_n(interp, array, field_name.str(), size, values);
}

void
UserInterface_TCL::readIdVariable(Tcl_Interp* interp, const char* array, const char* id, const char* variable,
                               int& size, double*& values)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  readVariable_impl_n(interp, array, field_name.str(), size, values);
}



// Id2 array, integer id1, id2, single value
// -----------------------------------------
void
UserInterface_TCL::readId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2, const char* variable,
                               char& value)
{
  strstream field_name;
  field_name << id1 << "," << id2 << "," << variable << ends;
  readVariable_impl_1(interp, array, field_name.str(), value);
}

void
UserInterface_TCL::readId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2, const char* variable,
                               char*& value)
{
  strstream field_name;
  field_name << id1 << "," << id2 << "," << variable << ends;
  readVariable_impl_1(interp, array, field_name.str(), value);
}

void
UserInterface_TCL::readId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2, const char* variable,
                               int& value)
{
  strstream field_name;
  field_name << id1 << "," << id2 << "," << variable << ends;
  readVariable_impl_1(interp, array, field_name.str(), value);
}

void
UserInterface_TCL::readId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2, const char* variable,
                               double& value)
{
  strstream field_name;
  field_name << id1 << "," << id2 << "," << variable << ends;
  readVariable_impl_1(interp, array, field_name.str(), value);
}


// Id2 array, integer id1, id2, list of values
// -------------------------------------------
void
UserInterface_TCL::readId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2, const char* variable,
                               int& size, char**& values)
{
  strstream field_name;
  field_name << id1 << "," << id2 << "," << variable << ends;
  readVariable_impl_n(interp, array, field_name.str(), size, values);
}

void
UserInterface_TCL::readId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2, const char* variable,
                               int& size, int*& values)
{
  strstream field_name;
  field_name << id1 << "," << id2 << "," << variable << ends;
  readVariable_impl_n(interp, array, field_name.str(), size, values);
}

void
UserInterface_TCL::readId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2, const char* variable,
                               int& size, double*& values)
{
  strstream field_name;
  field_name << id1 << "," << id2 << "," << variable << ends;
  readVariable_impl_n(interp, array, field_name.str(), size, values);
}


// ============
// Normal array
// ============

// Single value
// ------------

void
UserInterface_TCL::readVariable(Tcl_Interp* interp, const char* array, const char* variable,
                                char& value)
{
  readVariable_impl_1(interp, array, variable, value);
}

void
UserInterface_TCL::readVariable(Tcl_Interp* interp, const char* array, const char* variable,
                                char*& value)
{
  readVariable_impl_1(interp, array, variable, value);
}

void
UserInterface_TCL::readVariable(Tcl_Interp* interp, const char* array, const char* variable,
                                int& value)
{
  readVariable_impl_1(interp, array, variable, value);
}


void
UserInterface_TCL::readVariable(Tcl_Interp* interp, const char* array, const char* variable,
                                double& value)
{
  readVariable_impl_1(interp, array, variable, value);
}



// List of values
// --------------

void
UserInterface_TCL::readVariable(Tcl_Interp* interp, const char* array, const char* variable,
                                int& size, char**& values)
{
  readVariable_impl_n(interp, array, variable, size, values);
}

void
UserInterface_TCL::readVariable(Tcl_Interp* interp, const char* array, const char* variable,
                                 int& size, int*& values)
{
  readVariable_impl_n(interp, array, variable, size, values);
}


void
UserInterface_TCL::readVariable(Tcl_Interp* interp, const char* array, const char* variable,
                                int& size, double*& values)
{
  readVariable_impl_n(interp, array, variable, size, values);
}


// *** Function sets single value to array-type Tcl variable.
template <class T > void
UserInterface_TCL::readVariable_impl_1(Tcl_Interp* interp, const char* array, const char* variable,
                                       T& value)
{
  Tcl_Obj* array_obj = Tcl_NewStringObj(NULL, 0);
  Tcl_Obj* value_obj = NULL;

  int tcl_flags = TCL_GLOBAL_ONLY | TCL_LEAVE_ERR_MSG;

  // We have to create a tcl string-object to
  // carry the variable names. Here we use a fully
  // specified array name to create the name-object
  tcl_flags |= TCL_PARSE_PART1;
  strstream strm;
  strm << array << "(" << variable << ")" << ends;
  char* array_name = strm.str();
  Tcl_SetStringObj(array_obj, array_name, strlen(array_name));

  // Try to get existing variable handle
  // NOTE! This does not seem to be working!!! At least not for
  // strings like Default-directories in the settings files.
  // They all get the last string values if this is used
  // Try to get existing variable handle
  ;//value_obj = Tcl_ObjGetVar2(interp, array_obj, NULL, tcl_flags);

  Tcl_ResetResult(interp);

  // Update variable
  value_obj = Tcl_ObjGetVar2(interp, array_obj, NULL, tcl_flags);

  getTclObjValue(interp, value_obj, value);
}


// *** Function sets list of values to array-type Tcl variable.
template <class T > void
UserInterface_TCL::readVariable_impl_n(Tcl_Interp* interp, const char* array, const char* variable,
                                       int& size, T*& values)
{
  size = 0;
  values = NULL;

  Tcl_Obj* array_obj = Tcl_NewStringObj(NULL, 0);
  Tcl_Obj* value_obj = NULL;

  Tcl_Obj** list_val_obj = NULL;
  int list_len = 0;

  int tcl_flags = TCL_GLOBAL_ONLY | TCL_LEAVE_ERR_MSG;

  tcl_flags |= TCL_PARSE_PART1;
  tcl_flags |= TCL_LIST_ELEMENT;
  tcl_flags |= TCL_APPEND_VALUE;

  // In Tcl80 we have to create a tcl string-object to
  // carry the variable names. Here we use a fully
  // specified array name to create the name-object
  strstream strm;
  strm << array << "(" << variable << ")" << ends;
  char* array_name = strm.str();
  Tcl_SetStringObj(array_obj, array_name, strlen(array_name));

  value_obj = Tcl_ObjGetVar2(interp, array_obj, NULL, tcl_flags);


  // If unknown variables
  if (value_obj == NULL) {
    //strstream strm;
    //strm << "ERROR: Cannot read variable: ";
    //strm << array_name;
    //strm << ends;
    //UserInterface* gui = theControlCenter->getGui();
    //gui->showMsg(strm.str());
    return;
  }

  Tcl_ListObjGetElements(interp, value_obj, &list_len, &list_val_obj);

  if (list_len == 0)
    return;

  // Allocate result array
  size = list_len;
  values = new T[size];

  // Insert values
  for (int i = 0; i < size; i++) {
    getTclObjValue(interp, list_val_obj[i], values[i]);

    //if ( list_val_obj[i]->length == 0 ) {
    //  values[i] = (T)0;
    //}
  }

}


// *** Save model property (DB-name, DB-path etc) data
// from the gui-side
void
UserInterface_TCL::saveModelPropertyData(Model* model)
{
  sendCommandToGui(theInterp,"Interface::saveModelPropertyData");
}


void
UserInterface_TCL::selectBody(int bd1_id, int lr1_id, int bd2_id, int lr2_id, bool is_selected)
{
  ostrstream strm;

  strm << bd1_id << " " << lr1_id << " ";
  strm << bd2_id << " " << lr2_id << " ";
  strm << (int)is_selected << ends;

  sendCommandToGui(theInterp, "Interface::selectBody", strm.str());
}


void
UserInterface_TCL::selectBoundary(int bndr_id, int bd1_id, int lr1_id, int bd2_id, int lr2_id, bool extend)
{
  ostrstream strm;
  strm << bndr_id << " ";
  strm << bd1_id << " " << lr1_id << " ";
  strm << bd2_id << " " << lr2_id << " ";
  strm << extend << ends;

  sendCommandToGui(theInterp, "Interface::selectBoundary", strm.str());
}


// Prefixed commands sent by clients
int
UserInterface_TCL::sendCommandToGui(const char* cmd, const char* arg)
{
  if ( 0 == strcmp(cmd, "RESET_RENDERER") )
    return sendCommandToGui(theInterp, "Interface::rendererReset", NULL);

  else if ( 0 == strcmp(cmd, "SELECT_BODIES") )
    return sendCommandToGui(theInterp, "Interface::displayBodySelectPanel");

  else if ( 0 == strcmp(cmd, "SELECT_BOUNDARIES") )
    return sendCommandToGui(theInterp, "Interface::displayBoundarySelectPanel");

  else if ( 0 == strcmp(cmd, "SELECT_LABELS") )
    return sendCommandToGui(theInterp, "Interface::displayLabelSelectPanel");

  else if ( 0 == strcmp(cmd, "SET_ROTATE_PRIORITY_X") )
    return sendCommandToGui(theInterp, "Interface::setRotatePriority", "X");

  else if ( 0 == strcmp(cmd, "SET_ROTATE_PRIORITY_Y") )
    return sendCommandToGui(theInterp, "Interface::setRotatePriority", "Y");

  else if ( 0 == strcmp(cmd, "SET_ROTATE_PRIORITY_Z") )
    return sendCommandToGui(theInterp, "Interface::setRotatePriority", "Z");

  return false;
}


int
UserInterface_TCL::sendCommandToGui(Tcl_Interp* interp, const char* cmd, const char* arg)
{
  int i, len, rc;

  Tcl_ResetResult(interp);

  char* exec = "gui_exec ";

  {
    Tcl_DString dstring;
    char *buf;
    ostrstream strm;

    strm << exec;
    strm << "\"";
    strm << cmd;
    if ( arg != NULL ) {
      if ( arg[0] != '\0' ) {
         strm << tclCmdSeparator << arg;
      }
    }
    strm << "\"";
    strm << ends;

    buf = Tcl_ExternalToUtfDString(NULL, strm.str(), strlen(strm.str()), &dstring );
    rc = Tcl_Eval( interp, buf );
    Tcl_DStringFree( &dstring );
  }

  if (interp->result[0] != '\0') {

    char err_buf[256];
    err_buf[255] = '\0';
    strncpy(err_buf, interp->result, 255);

    char cmd_buf[256];
    cmd_buf[255] = '\0';
    strncpy(cmd_buf, cmd, 255);

    char arg_buf[256];
    arg_buf[0] = '\0';
    if ( arg != NULL && arg[0] != '\0' ) {
      arg_buf[255] = '\0';
      strncpy(arg_buf, arg, 255);
    }

    // Replace quotes with spaces so that we can safetly use showMsg function
    // to display the error message
    len = strlen(err_buf);
    for (i = 0; i < len; i++) {
      if (err_buf[i] == '"')
        err_buf[i] = ' ';
    }

    len = strlen(cmd_buf);
    for (i = 0; i < len; i++) {
      if (cmd_buf[i] == '"')
        cmd_buf[i] = ' ';
    }

    len = strlen(arg_buf);
    for (i = 0; i < len; i++) {
      if (arg_buf[i] == '"')
        arg_buf[i] = ' ';
    }

    UserInterface* gui = theControlCenter->getGui();

    // Command info
    strstream strm;
    strm << "***WARNING INVALID COMMAND:  " << cmd_buf;
    if (arg != NULL)
      strm << "   (arg: " << arg_buf << ")   ";
    strm << "   because:" << ends;
    gui->showMsg(strm.str());

    // Error info from Tcl
    gui->showMsg(err_buf, 1);

    return 0;
  }

  return 1;
}


// Function sets boundary conditions for parent object in model
//
void
UserInterface_TCL::setBoundaryConditions(Tcl_Interp *interp, Model* model,
          MultiIdTable& bc_table)
{

  MultiIdTable::iterator pos = bc_table.begin();
  MultiIdTable::iterator end = bc_table.end();

  while ( pos != end ) {

    int pid = (*pos).first;
    int beg_id = (*pos).second;

    BodyElementGroup* beg = (BodyElementGroup*)model->getModelObjectById(beg_id);

    if ( beg != NULL ) {
      beg->setBoundaryConditionId(pid);
    }

    pos++;
  }
}


#if 0
// Function sets boundary conditions for edges
// Not in use!
void
UserInterface_TCL::setBoundaryConditionsForEdges(Tcl_Interp *interp, Model* model,
          MultiIdTable& bc_table)
{

  MultiIdTable::iterator pos = bc_table.begin();
  MultiIdTable::iterator end = bc_table.end();

  while ( pos != end ) {

    int pid = (*pos).first;
    int be_id = (*pos).second;

    BodyElement* be = model->getEdgeById(be_id);

    if ( be != NULL ) {
      be->setBoundaryConditionId(pid);
    }

    pos++;
  }
}


// Function sets boundary conditions for faces
// Not in use!
void
UserInterface_TCL::setBoundaryConditionsForFaces(Tcl_Interp *interp, Model* model,
          MultiIdTable& bc_table)
{

  MultiIdTable::iterator pos = bc_table.begin();
  MultiIdTable::iterator end = bc_table.end();

  while ( pos != end ) {

    int pid = (*pos).first;
    int be_id = (*pos).second;

    BodyElement* be = model->getFaceById(be_id);

    if ( be != NULL ) {
      be->setBoundaryConditionId(pid);
    }

    pos++;
  }
}


// Function sets boundary conditions for vertices.
// Not in use!
void
UserInterface_TCL::setBoundaryConditionsForVertices(Tcl_Interp *interp, Model* model,
          MultiIdTable& bc_table)
{

  MultiIdTable::iterator pos = bc_table.begin();
  MultiIdTable::iterator end = bc_table.end();

  while ( pos != end ) {

    int pid = (*pos).first;
    int be_id = (*pos).second;

    BodyElement* be = model->getVertexById(be_id);

    if ( be != NULL ) {
      be->setBoundaryConditionId(pid);
    }

    pos++;
  }
}
#endif


void
UserInterface_TCL::setBoundarySelectionMode(int bndr_id, bool is_selected, bool do_update)
{
  ostrstream strm;
  strm << bndr_id << " " << (int)is_selected << " " << (int)do_update << ends;

  sendCommandToGui(theInterp, "Interface::setBoundarySelectionMode", strm.str());
}


void
UserInterface_TCL::setCurrentMeshH(double mesh_h)
{
  strstream strm;
  strm << mesh_h << ends;

  sendCommandToGui(theInterp,"Interface::setCurrentMeshH", strm.str());
}


void
UserInterface_TCL::setExceptionThrown()
{
  sendCommandToGui(theInterp,"Interface::setExceptionThrown");
}


void
UserInterface_TCL::setInitialMeshH(double mesh_h)
{
  strstream strm;
  strm << mesh_h << ends;

  sendCommandToGui(theInterp,"Interface::setInitialMeshH", strm.str());
}


void
UserInterface_TCL::setInitialState()
{
  sendCommandToGui(theInterp,"Interface::setInitialState");
}


void
UserInterface_TCL::setMeshEdited()
{
  sendCommandToGui(theInterp, "Interface::setMeshEdited", "");
}


void
UserInterface_TCL::setMeshExists()
{
  sendCommandToGui(theInterp, "Interface::setMeshExists", "");
}


void
UserInterface_TCL::setMeshInputUnit(double unit)
{
  ostrstream strm;
  strm << unit << ends;

  sendCommandToGui(theInterp, "Interface::setMeshInputUnit", strm.str());
}


void
UserInterface_TCL::setModelHasElmerMesh()
{
  sendCommandToGui(theInterp, "Interface::setModelHasElmerMesh", "");
}


void
UserInterface_TCL::setModelHasMeshParameter()
{
  sendCommandToGui(theInterp, "Interface::setModelHasMeshParameter", "");
}


void
UserInterface_TCL::setModelHasMatcDefinitions()
{
  sendCommandToGui(theInterp, "Interface::setModelHasMatcDefinitions", "");
}


// Target data (database, solver etc.) needs update
void
UserInterface_TCL::setNeedsUpdate(const char* target)
{
  sendCommandToGui(theInterp, "Interface::setNeedsUpdate", (char*)target);
}


int
UserInterface_TCL::setParameterData(Model* model, ecif_parameterType param_type,
                                    Tcl_Interp* interp, const char* array_name)
{
  // Read parameter ids
  int size = 0;
  int* ids = NULL;
  readVariable(interp, array_name, "ids", size, ids);

  // Mark parameter for update, so that non-set parameters
  // can be deleted after update
  model->processParametersBeforeUpdate(param_type);


  char* data_buffer = "";
  char* name_buffer = "";

  int oid = NO_INDEX;
  int prtag = NO_INDEX;
  char* prtp = NULL;
  objectType prtype = OT_NONE;

  int pid;
  int attach_mode;

  for (int i = 0; i < size; i++) {

    pid = ids[i];

    //-Parameter-name
    readIdVariable(interp, array_name,  pid, "name", name_buffer);

    //-Parent object id, tag and type
    readIdVariable(interp, array_name,  pid, "oid", oid);

    if ( oid != NO_INDEX ) {

      readIdVariable(interp, "ObjectTable",  oid, "tg", prtag);
      readIdVariable(interp, "ObjectTable",  oid, "tp", prtp);

      if (prtp != NULL ) {

        if ( 0 == strcmp(prtp, "B") )
          prtype = OT_BODY;
        else if ( 0 == strcmp(prtp, "BL") )
          prtype = OT_BODY_LAYER;
        else if ( 0 == strcmp(prtp, "BP") )
          prtype = OT_BODYPAIR;
        else if ( 0 == strcmp(prtp, "F") )
          prtype = OT_FACE;
        else if ( 0 == strcmp(prtp, "E") )
          prtype = OT_EDGE;
        else if ( 0 == strcmp(prtp, "V") )
          prtype = OT_VERTEX;

        // NOTE: All these are marked as element groups in model
        else if ( 0 == strcmp(prtp, "FG") )
          prtype = OT_ELEMENT_GROUP;
        else if ( 0 == strcmp(prtp, "EG") )
          prtype = OT_ELEMENT_GROUP;
        else if ( 0 == strcmp(prtp, "VG") )
          prtype = OT_ELEMENT_GROUP;

        else
          prtype = OT_NONE;
      }
    }

    //-Parameter data
    readIdVariable(interp, array_name,  pid, "data", data_buffer);

    //---Save data into model
    //model->setParameter(param_type, pid, oid, prtag, prtype, data_buffer, name_buffer);
    model->setParameter(param_type, pid, oid, data_buffer, name_buffer);

    /*    delete[] data_buffer; data_buffer = NULL;
    delete[] name_buffer; name_buffer = NULL;
    delete[] prtp; prtp = NULL; */
  }

  // Delete "old" parameters, ie. those which were not set in this update
  model->processParametersAfterUpdate(param_type);

  delete[] ids;

  return TCL_OK;

}


void
UserInterface_TCL::setParameterFieldValueState(int parameter_id, const char* field_name,
                                               bool has_value, bool value_has_changed)
{
  ostrstream strm;
  strm << parameter_id << ' '
       << (char*)field_name << ' '
       << has_value << ' '
       << value_has_changed
       << ends;
  sendCommandToGui(theInterp, "Interface::setParameterFieldValueState", strm.str());
}


void
UserInterface_TCL::setTimestamp(ecif_parameterType parameter, char* ts)
{
  char* array = NULL;
  char* field = NULL;

  switch (parameter) {
  case ECIF_GRID_PARAMETER:
    array = "Model";
    field = "MeshParameter,timestamp";
    break;
  }

  if (array != NULL && field != NULL) {
    writeVariable(theInterp, array, field, ts);
  }
}


void
UserInterface_TCL::setWindowTitle(char* title)
{
  sendCommandToGui( theInterp, "Util::setMainWindowTitle", title);
}


// Target data (database, solver etc.) was succesfully updated
void
UserInterface_TCL::setWasUpdated(const char* target)
{
  sendCommandToGui(theInterp, "Interface::setWasUpdated", (char*)target);
}


int
UserInterface_TCL::showMsg(char* message, short extra_line_feeds, bool append)
{
  // Replace double quotes with single quotes in order to show the
  // message in Tcl
  int len = strlen(message);
  for (int i = 0; i < len; i++) {
    if ( message[i] == '\"' ) {
      message[i] = '\'';
    }
  }

  ostrstream strm;
  strm << message << tclArgSeparator;
  strm << extra_line_feeds << tclArgSeparator;
  strm << append;
  strm << ends;

  return sendCommandToGui(theInterp, "Interface::showMessage", strm.str());
}


void
UserInterface_TCL::showProgressMsg(Timer& timer, int frequency_nbr,
                                   int nbr, int total_nbr,
                                   char* text1, char* text2)
{
  static char time_string[32];

  // If no need to display progress message
  if ( frequency_nbr > 0 &&
       0 != (nbr % frequency_nbr) &&
       nbr != total_nbr
     )
    return;

  ostrstream strm;

  strm << endl;
  if (text1 != NULL)
    strm << text1;
  strm << total_nbr;
  if (text2 != NULL)
    strm << text2;

  // If at the end (last item beig processed), display total time
  if ( nbr == total_nbr ) {

    double time = timer.getLapTime(WALL_TIME);
    timer.formTimeString(time, time_string);
    strm << " took " << time_string << ends;
    showMsg(strm.str(), 0, false);

  // Otherwise frequency criteria is met, display %-msg
  } else {

    double done = (int((1000.0 * nbr) / total_nbr)) / 10.0;
    strm << "(";
    strm << setiosflags(ios::right | ios::fixed);
    strm << setw(4) << setprecision(1);
    strm << done << "%)";
    strm << ends;
    // If first %-msg, append text
    if (nbr == frequency_nbr) {
      showMsg(strm.str(), 0, true);
    // otherwise overwrite
    } else {
      showMsg(strm.str(), 0, false);
    }
  }
}


void
UserInterface_TCL::showUsedTimeMsg(double time, char* text,
                                   short extra_line_feeds, bool append)
{
  static char time_string[32];
  formTimeString(time, time_string);
  strstream strm;
  strm << text
       << " took "
       << time_string
       << ends;

  showMsg(strm.str(), extra_line_feeds, append);
}


void
UserInterface_TCL::showUsedTimeMsg(double time, char* text1, int nof_objects, char* text2,
                       short extra_line_feeds, bool append)
{
  static char time_string[32];
  formTimeString(time, time_string);
  strstream strm;
  strm << text1
       << " " << nof_objects << " "
       << text2
       << " took "
       << time_string
       << ends;

  showMsg(strm.str(), extra_line_feeds, append);
}


void
UserInterface_TCL::start(int argc, char** argv)
{
  char* elmer_home = (char *) Tcl_GetVar2(theInterp, "env", "ELMER_HOME", glob_flag);
  char* elmer_front_home = (char *) Tcl_GetVar2(theInterp, "env", "ELMER_FRONT_HOME", glob_flag);

  char front_tcl_path[] = "/tcl";

  // ===================================
  // SRCIPT FILE AND SCRIPT LIBRARY PATH
  // ===================================

  // These are buffers for the scriptfile and lib-path
  // picking streams
  char file_buffer[emf_MAX_STRING_LEN];
  char lib_buffer[emf_MAX_STRING_LEN];

  ostrstream file_strm(file_buffer, emf_MAX_STRING_LEN); // Start script name
  ostrstream lib_strm(lib_buffer, emf_MAX_STRING_LEN);   // Front TCL lib

  bool path_found = false;

  //---Find correct path for the start-script

  //-1. try directly under local path:
  if (!path_found) {
    file_strm.seekp(0);
    lib_strm.seekp(0);
    file_strm <<  "./"
              << controlSideScript
              << ends;
    lib_strm  <<  "./"
              << ends;
    // Try to open start-script!
    printf("Trying %s\n",file_buffer);
    if ( NULL != fopen(file_buffer, "r") ) {
      path_found = true;
    }
  }

  //-2. try directly under local TCL-path:
  if (!path_found) {
    file_strm.seekp(0);
    lib_strm.seekp(0);
    file_strm << "./tcl/"
              << controlSideScript
              << ends;
    lib_strm  << "./tcl/"
              << ends;
    // Try to open start-script!
    printf("Trying %s\n",file_buffer);
    if ( NULL != fopen(file_buffer, "r") ) {
      path_found = true;
    }
  }

  //-3. try directly under ELMER_FRONT_HOME environment variables
  if (!path_found) {
    if (elmer_front_home != NULL) {
      file_strm.seekp(0);
      lib_strm.seekp(0);
      file_strm << elmer_front_home
                << "/"
                << controlSideScript
                << ends;
      lib_strm  << elmer_front_home
                << "/"
                << ends;
      printf("Trying %s\n",file_buffer);
      // Try to open start-script!
      if ( NULL != fopen(file_buffer, "r") ) {
        path_found = true;
      }
    }
  }

  //-4. try TCL-path via ELMER_FRONT_HOME environment variables
  if (!path_found) {
    if (elmer_front_home != NULL) {
      file_strm.seekp(0);
      lib_strm.seekp(0);
      file_strm << elmer_front_home
                << front_tcl_path << "/"
                << controlSideScript
                << ends;
      lib_strm  << elmer_front_home << "/"
                << front_tcl_path
                << "/"
                << ends;
      printf("Trying %s\n",file_buffer);
      // Try to open start-script!
      if ( NULL != fopen(file_buffer, "r") ) {
        path_found = true;
      }
    }
  }

  //-5. try TCL-path via installation prefix
  if (!path_found) {
    file_strm.seekp(0);
    lib_strm.seekp(0);
    
    file_strm << ELMER_FRONT_PREFIX 
              << "/share/elmerfront"
	      << front_tcl_path  << "/"
	      << controlSideScript
	      << ends;
    lib_strm  << ELMER_FRONT_PREFIX 
	      << "/share/elmerfront/tcl/"
	      << ends;
    // Try to open start-script!
    printf("Trying %s\n",file_buffer);
    if ( NULL != fopen(file_buffer, "r") ) 
    {
      path_found = true;
    }
  }


  //---Message to the user
  // ERROR
  if (!path_found) {

    // Error message to file
    ((ofstream*)debugFile)->open("ElmerFront.log", ios::out);
    if ( !debugFile->fail() ) {
      *debugFile  << "Can't run Elmer Front program." << endl;
      *debugFile  << "Script " << controlSideScript << " not found!" << endl;
    }

    // Error message to stdout
    cerr        << "Can't run Elmer Front program." << endl;
    cerr        << "Script " << controlSideScript << " not found!" << endl;
    exit(1);
  }

  // OK
  else {
    cerr << "Running Elmer Front program." << endl
         << "Start script=" << file_buffer << endl;
  }

  //---We continue
  //-Copy libpath-name to attribute *tclScriptPath*
  int len = strlen(lib_buffer);
  tclScriptPath = new char[len+1];
  tclScriptPath[len] = '\0';
  strcpy(tclScriptPath, lib_buffer);


  //=====================
  // LOAD GUI MAIN SCRIPT
  //=====================

  //---Load the SCRIPT specified in the fileName argument.
  int code = Tcl_EvalFile(theInterp, file_buffer);

  //--If we can't load the script (= start CONTROL-SIDE interpreter)
  if (code != TCL_OK) {

    char* p = (char *)Tcl_GetVar(theInterp, "errorInfo", glob_flag);

    if ((p == NULL) || (*p == '\0')) {
      p = theInterp->result;
    }

    ((ofstream*)debugFile)->open("ElmerFront.log", ios::out);
    *debugFile  << "Can't run Elmer Front program. " << endl;
    cerr        << "Can't run Elmer Front program. " << endl;
    *debugFile  << "Not able to run Tcl-interpreter: " << p << endl;
    cerr        << "Not able to run Tcl-interpreter: " << p << endl;
    exit(1);

  }

  //---Set path value for in the Tcl-environment (control-side, not GUI-side!)
  writeVariable(theInterp, "Info", "frontScriptPath", tclScriptPath);

  //---Set Elmer Front version numbers for the GUI-side
  writeVariable(theInterp, "Info", "FRONT_PREVIOUS_INPUT_VERSION_NBR", -1);
  writeVariable(theInterp, "Info", "FRONT_INPUT_VERSION_NBR", -1);
  writeVariable(theInterp, "Info", "FRONT_MAIN_VERSION_NBR", ECIF_VERSION_NBR);

  //---Handle posssible command line arguments
  // NOTE: Check ecif_tclMainScript.tcl (proc handleCommandLineArguments)
  // for available options!
  // NOTE: We handle command line args here "quitely" to check if a settings-file name
  // was given as an argument!
  // NOTE: arguments will be reread in the gui-side, when all the nice
  // message machinery has been loaded.
  char* cmd_arguments = Tcl_Concat(argc, argv);

  // Change to behave like Tcl-like dirnames
  // NOTE: We should not use arguments which depend on back-slash!!!
  for (int i = 0; i < strlen(cmd_arguments); i++) {
    if ( cmd_arguments[i] == '\\' )
      cmd_arguments[i] = '/';
  }

  // Transfer args to the gui side and read "quietly"
  Tcl_SetVar2(theInterp, "Info", "commandLineArgs", cmd_arguments, glob_flag);
  sendCommandToGui(theInterp, "handleCommandLineArgs");

  //==========
  // START GUI
  //==========

  //---Init Tk-environment
  if (Tk_Init(theInterp) == TCL_ERROR) {

    ((ofstream*)debugFile)->open("ElmerFront.log", ios::out);
    *debugFile  << "Can't run Elmer Front program. " << endl;
    cerr        << "Can't run Elmer Front program. " << endl;
    *debugFile  << "Not able to run Tk-interpreter" << endl;
    cerr        << "Not able to run Tk-interpreter" << endl;
    exit(1);

  }

  //---Start GUI
  sendCommandToGui(theInterp, "startGUI", NULL);

  Renderer_OGL::setRendererInfo();

  //---Start Tcl-loop
  while ( start_Tcl_MainLoop() );
}


bool
UserInterface_TCL::start_Tcl_MainLoop()
{
  // Start Tcl-loop and don't continue when it is stops normally
  Tcl_MainLoop();

  return false;
}


int
UserInterface_TCL::from_tk_ColorHex2Name(ClientData clientData, Tcl_Interp *interp,
          int argc, char* argv[])
{
  Model* model = theControlCenter->getModel();

  if (model == NULL) {
    TheUI->errMsg(0, "ColorHex2Name: Command not legal. No model exists!");
    return TCL_OK;
  }

  char buffer[64] = "";
  char* hex_value = getCommandArguments(interp);
  int code = model->rgbColor2Id(6, hex_value);
  model->getColorName(code, buffer);
  sendCommandToGui(interp, "Interface::getColorName", buffer);
  return TCL_OK;
}



// Function writes all parameter data into Tk variables.
void
UserInterface_TCL::to_tk_WriteAllParamDataPre(Tcl_Interp* interp, const Model& model)
{
  int nof_params = 20; // Nof parameters types
  int nof_pstrms =  3; // Nof separate streams needed per parameter type

  int ids_buffer[10000];

  // NOTE: Keep these lists consistent!!!
  //
  ecif_parameterType param_types[] = {
    ECIF_BODY_FORCE,
    ECIF_BODY_PARAMETER,
    ECIF_BOUNDARY_CONDITION,
    ECIF_BOUNDARY_PARAMETER,
    ECIF_CALCULATOR,
    ECIF_CONSTANT,
    ECIF_COORDINATE,
    ECIF_DATAFILE,
    ECIF_EQUATION,
    ECIF_EQUATION_VARIABLE,
    ECIF_GRID_H,
    ECIF_GRID_PARAMETER,
    ECIF_INITIAL_CONDITION,
    ECIF_MATERIAL,
    ECIF_MODEL_PARAMETER,
    ECIF_SIMULATION_PARAMETER,
    ECIF_SOLVER,
    ECIF_SOLVER_CONTROL,
    ECIF_TIMESTEP,
    ECIF_USER_SETTING
  };

  char* param_arries[] = {
    "BodyForce",
    "BodyParameter",
    "BoundaryCondition",
    "BoundaryParameter",
    "Calculator",
    "Constant",
    "Coordinate",
    "Datafile",
    "Equation",
    "EquationVariable",
    "GridH",
    "GridParameter",
    "InitialCondition",
    "Material",
    "ModelParameter",
    "SimulationParameter",
    "Solver",
    "SolverControl",
    "Timestep",
    "UserSetting"
  };

  //Ups! An ugly cast needed owing to model's constant type. !!!***!!!
  Model* mdl = (Model*) &model;

  for (int i = 0; i < nof_params; i++) {

    // Initial value
    writeVariable(interp, param_arries[i], "ids", 0, (int*)NULL);
    writeVariable(interp, param_arries[i], "nextNewParameterId", mdl->getNextNewParameterId(param_types[i]));

    int index = 0;
    while (true) {
      Parameter* param = mdl->getParameter(index++, param_types[i]);
      if (param==NULL) break;
      int pid = param->ID();
      ids_buffer[index - 1] = pid;
      writeIdVariable(interp, param_arries[i], pid, "data", (char*)param->getValue());
      writeIdVariable(interp, param_arries[i], pid, "name", (char*)param->getName());
      writeIdVariable(interp, param_arries[i], pid, "oid", param->getParentId());
    }

    if (index > 1 ) {
      writeVariable(interp, param_arries[i], "ids", index-1, ids_buffer);
    }

  }


  // Process parameter data in Gui side
  sendCommandToGui(interp, "Interface::applyAllParameterDataPre");
}


// Parameter mask related stuff
void
UserInterface_TCL::to_tk_WriteAllParamDataPost(Tcl_Interp* interp, const Model& model)
{
  // Process parameter data in Gui side
  sendCommandToGui(interp, "Interface::applyAllParameterDataPost");

  // Write status for the following parameters
  to_tk_WriteStatusEquations(interp, model);
  to_tk_WriteStatusBodyForces(interp, model);
  to_tk_WriteStatusBoundaryConditions(interp, model);
  to_tk_WriteStatusInitialConditions(interp, model);
  to_tk_WriteStatusMaterials(interp, model);
  to_tk_WriteStatusTimesteps(interp, model);
  to_tk_WriteStatusMeshes(interp, model);

}


// Write body layers to be used as Mesh structure panel objects
// NOTE: Layers are read via bodies and only unique layer-tags
// are inserted
//
void
UserInterface_TCL::to_tk_WriteBodyLayerData(Tcl_Interp* interp, const Model& model)
{
  //Ups! An ugly cast needed owing to model's constant type. !!!***!!!
  Model* mdl = (Model*) &model;

  sendCommandToGui(interp, "Interface::resetObjectTableByType", "BL");

  //----Loop all bodies
  //
  int index = 0;
  while (true) {

    Body* bd = mdl->getBody(index++);

    if ( bd == NULL ) break;

    // This is needed to store only unique layer tags
    //
    IdsSet layerTags;
    std::pair<IdsSet::iterator, bool> lrPair;

    //--Loop all layers in the body
    //
    for (int layer = 0; layer < bd->getNofLayers(); layer++) {

      int lr_id = bd->getLayerId(layer);
      BodyLayer* lr = mdl->getBodyLayerById(lr_id);

      if (lr == NULL ) continue;

      // Try to insert a new layer tag, check if it is already inserted
      lrPair = layerTags.insert(lr->Tag());

      // Layer tag exists for the body
      if ( !lrPair.second ) continue;

      int id = lr->Id();

      const char* nm;

      if ( lr->hasName() ) {
        nm = lr->getName();
      } else {

        nm = bd->getName();

        if ( bd->getNofLayers() > 1 ) {
          strstream strm;
          strm << nm << "-Lr" << lr->Tag() << ends;
          nm = strm.str();
        }
      }

      // Add layer id
      writeVariable(interp, "ObjectTable", "ids", 1, &id, false);

      writeIdVariable(interp, "ObjectTable", id, "tg", lr->getBodyTag());
      writeIdVariable(interp, "ObjectTable", id, "tp", "BL");

      if ( lr->isClosed() )
        writeIdVariable(interp, "ObjectTable", id, "cl", "CSD");
      else if ( lr->isOpen() )
        writeIdVariable(interp, "ObjectTable", id, "cl", "OPN");
      else
        writeIdVariable(interp, "ObjectTable", id, "cl", "");

      writeIdVariable(interp, "ObjectTable", id, "nm", nm);
      writeIdVariable(interp, "ObjectTable", id, "bdId", lr->getBodyId());
      writeIdVariable(interp, "ObjectTable", id, "sbIds", "");
      writeIdVariable(interp, "ObjectTable", id, "sbTgs", "");
      writeIdVariable(interp, "ObjectTable", id, "gr", NO_INDEX);
      writeIdVariable(interp, "ObjectTable", id, "grIds", "");
      writeIdVariable(interp, "ObjectTable", id, "grMshIndcs", "");
      writeIdVariable(interp, "ObjectTable", id, "excldMsh", 0);
      writeIdVariable(interp, "ObjectTable", id, "excldMshIndcs", "");
      writeIdVariable(interp, "ObjectTable", id, "accptStrMsh", (int)bd->acceptsStructuredMesh(layer));

      // NOTE: For layers sub ids are normal boundary ids (no groups etc.)
      // because only htese are relevant in mesh structure panel where
      // layer are used!
      //
      int nof_sb_ids;
      int* sb_ids = NULL;

      bd->getElementIds(layer, nof_sb_ids, sb_ids);

      for (int i = 0; i< nof_sb_ids; i++) {
        BodyElement* be = mdl->getBodyElementById(sb_ids[i]);

        if ( be == NULL ) continue;

        int be_id = be->Id();
        int be_tg = be->Tag();

        writeIdVariable(interp, "ObjectTable", id, "sbIds", 1, &be_id, false);
        writeIdVariable(interp, "ObjectTable", id, "sbTgs", 1, &be_tg, false);
      }

      delete[] sb_ids;

      int nof_gids = lr->getNofGridParameterIds();
      if ( nof_gids > 0 ) {
        writeIdVariable(interp, "ObjectTable", id, "grIds", nof_gids, lr->getGridParameterIds());
        writeIdVariable(interp, "ObjectTable", id, "grMshIndcs", nof_gids, lr->getGridParameterMeshIndices());
      }

      int nof_mids = lr->getNofExcludedMeshes();
      if ( nof_mids > 0 ) {
        writeIdVariable(interp, "ObjectTable", id, "excldMshIndcs", nof_mids, lr->getExcludedMeshIndices());
      }
    }
  }
}


// Function writes body-level data into Tk variables.
void
UserInterface_TCL::to_tk_WriteBodyData(Tcl_Interp* interp, const Model& model)
{
  //Ups! An ugly cast needed owing to model's constant type. !!!***!!!
  Model* mdl = (Model*) &model;

  int nof_meshes = mdl->getNofMeshes();

  char color_hex_string[6];
  Color4 color;

  // Clean old body object data
  //
  sendCommandToGui(interp, "Interface::resetObjectTableByType", "B");   // Body, normal close body
  sendCommandToGui(interp, "Interface::resetObjectTableByType", "BP");  // Body, pair of non-virtual bodies

  //----Loop all bodies
  int index = 0;
  while (true) {

    Body* bd = mdl->getBody(index++);

    if (bd==NULL) break;

    int id = bd->Id();

    // Body color
    bd->getColor(color);
    int color_id = mdl->rgbColor2Id(color);
    mdl->rgbColorId2Hex(color_id, color_hex_string, 6);

    // Add body id
    writeVariable(interp, "ObjectTable", "ids", 1, &id, false);

    // Init body variables in OBjectTable
    //
    // NOTE: Body/Bodypair subIds etc. for normal bodies are written
    //       via boundaries/boundary groups!
    //
    writeIdVariable(interp, "ObjectTable", id, "tg", bd->Tag());
    writeIdVariable(interp, "ObjectTable", id, "tp", "B");

    if ( !bd->isVirtual() ) {

      if ( bd->isClosed() )
        writeIdVariable(interp, "ObjectTable", id, "cl", "CSD");
      else if ( bd->isOpen() )
        writeIdVariable(interp, "ObjectTable", id, "cl", "OPN");
      else
        writeIdVariable(interp, "ObjectTable", id, "cl", "");

    } else  {
      writeIdVariable(interp, "ObjectTable", id, "cl", "VIR");
    }

    writeIdVariable(interp, "ObjectTable", id, "nm", bd->getName());
    writeIdVariable(interp, "ObjectTable", id, "blIds", "");  // Grid group ids
    writeIdVariable(interp, "ObjectTable", id, "clr", color_hex_string);
    writeIdVariable(interp, "ObjectTable", id, "msk", "");
    writeIdVariable(interp, "ObjectTable", id, "prId", id);   // Self reference!
    writeIdVariable(interp, "ObjectTable", id, "pr1Id", id);  // Self reference!
    writeIdVariable(interp, "ObjectTable", id, "pr2Id", id);  // Self reference!
    writeIdVariable(interp, "ObjectTable", id, "sbIds", "");
    writeIdVariable(interp, "ObjectTable", id, "sbTgs", "");
    writeIdVariable(interp, "ObjectTable", id, "bodyp", bd->getBodyParameterId());
    writeIdVariable(interp, "ObjectTable", id, "eq", bd->getEquationId());
    writeIdVariable(interp, "ObjectTable", id, "mt", bd->getMaterialId());
    writeIdVariable(interp, "ObjectTable", id, "ic", bd->getInitialConditionId());
    writeIdVariable(interp, "ObjectTable", id, "bf", bd->getBodyForceId());

    // Write layer info and grid info on each GridGroup
    //
    int layer = -1;
    while (true) {
      if (!bd->selectLayer(++layer)) break;
      int layer_id = bd->getLayerId(layer);
      int layer_tag = bd->getLayerTag(layer);
      //--Add layer info for the body
      writeIdVariable(interp, "ObjectTable", id, "blIds", 1, &layer_id, false);
    }

    // Body pairs for this body
    // ========================
    int bp_index = 0;
    while (true) {

      BodyPair* bp = mdl->getBodyPair(bp_index++);

      if (bp==NULL) break;

      const Body* body1;
      const Body* body2;

      bp->getBodies(body1, body2);

      Body* bd1 = (Body*) body1;
      Body* bd2 = (Body*) body2;

      // "bd" must be the first body
      if ( bd1 != bd ) {
        continue;
      }

      int id = bp->Id();

      writeVariable(interp, "ObjectTable", "ids", 1, &id, false);
      writeIdVariable(interp, "ObjectTable", id, "tg", bp->Tag());
      writeIdVariable(interp, "ObjectTable", id, "tp", "BP");
      writeIdVariable(interp, "ObjectTable", id, "cl", "");
      writeIdVariable(interp, "ObjectTable", id, "nm", bp->getName());
      writeIdVariable(interp, "ObjectTable", id, "sbIds", "");
      writeIdVariable(interp, "ObjectTable", id, "sbTgs", "");
      writeIdVariable(interp, "ObjectTable", id, "msk", "");
      writeIdVariable(interp, "ObjectTable", id, "prId", id);  // Self reference!
      writeIdVariable(interp, "ObjectTable", id, "pr1Id", bd1->Id());
      writeIdVariable(interp, "ObjectTable", id, "pr2Id", bd2->Id());

    } // All body-pairs for the body

  } // All bodies

  sendCommandToGui(interp, "Interface::applyBodyData");
}


// *** Write body level info to gui.
void
UserInterface_TCL::to_tk_WriteBodyInfoData(Tcl_Interp* interp, const Model& model)
{
  const char* name;
  RangeVector rv;
  int nof_mesh_elements;

  Model* mdl = (Model*) &model;

  int index = 0;
  while (true) {
    Body* body = mdl->getBody(index++);
    if (body==NULL) break;

    int id = body->Id();
    body->getRangeVector(rv);
    nof_mesh_elements = body->getNofMeshElements();

    writeIdVariable(interp, "ObjectTable", id, "mnX", rv[0]);
    writeIdVariable(interp, "ObjectTable", id, "mxX", rv[1]);
    writeIdVariable(interp, "ObjectTable", id, "mnY", rv[2]);
    writeIdVariable(interp, "ObjectTable", id, "mxY", rv[3]);
    writeIdVariable(interp, "ObjectTable", id, "mnZ", rv[4]);
    writeIdVariable(interp, "ObjectTable", id, "mxZ", rv[5]);
    writeIdVariable(interp, "ObjectTable", id, "nofMshElm", nof_mesh_elements);
  }
}


// *** Write body level mesh info to gui.
void
UserInterface_TCL::to_tk_WriteBodyMeshInfoData(Tcl_Interp* interp, const Model& model)
{
  char* name;
  RangeVector rv;
  int nof_mesh_elements;

  Model* mdl = (Model*) &model;

  int index = 0;
  while (true) {
    Body* body = mdl->getBody(index++);
    if (body==NULL) break;
    int id = body->Id();
    nof_mesh_elements = body->getNofMeshElements();
    writeIdVariable(interp, "ObjectTable", id, "nofMshElm", nof_mesh_elements);
  }
}


// Function writes boundaries related data into Tk variables.
void
UserInterface_TCL::to_tk_WriteBoundaryData(Tcl_Interp* interp, const Model& model)
{
  //Ups! An ugly cast needed owing to model's constant type. !!!***!!!
  Model* mdl = (Model*) &model;

  // Clean old boundary object data
  //
  sendCommandToGui(interp, "Interface::resetObjectTableByType", "V");
  sendCommandToGui(interp, "Interface::resetObjectTableByType", "E");
  sendCommandToGui(interp, "Interface::resetObjectTableByType", "F");

  to_tk_WriteBoundaryElements(interp, mdl, OT_VERTEX);
  to_tk_WriteBoundaryElements(interp, mdl, OT_EDGE);
  to_tk_WriteBoundaryElements(interp, mdl, OT_FACE);

  sendCommandToGui(interp, "Interface::applyBoundaryData");
}


// Function writes boundaries related data into Tk variables.
void
UserInterface_TCL::to_tk_WriteBoundaryElements(Tcl_Interp* interp, Model* model, objectType btype)
{
  int nof_meshes = model->getNofMeshes();

  bool only_active = true;
  BodyElement* be = NULL;

  int index = 0;
  while (true) {

    switch (btype) {
    case OT_FACE:
      be = model->getFace(index++, only_active);
      break;
    case OT_EDGE:
      be = model->getEdge(index++, only_active);
      break;
    case OT_VERTEX:
      be = model->getVertex(index++, only_active);
      break;
    }

    if (be==NULL) break;

    // Find parent for the boundary
    //
    int bd1_id = be->getParentId(1);
    int bd2_id = be->getParentId(2);

    int bd1_lr = be->getParentLayer(1);
    int bd2_lr = be->getParentLayer(2);

    Body* bd1 = model->getBodyById(bd1_id);
    Body* bd2 = model->getBodyById(bd2_id);

    BodyPair* bp = NULL;

    bool is_intra_layer = be->isIntraLayerBoundary();

    if ( bd1 != NULL && bd2 != NULL ) {
      bp = model->getBodyPairById((const Body*)bd1, (const Body*)bd2);
    }

    int parent_oid;

    if ( bp != NULL ) {
      parent_oid = bp->Id();

    } else if ( bd1 != NULL ) {
      parent_oid = bd1->Id();

    } else if ( bd2 != NULL ) {
      parent_oid = bd2->Id();

    } else {
      parent_oid = NO_INDEX;
    }

    if ( parent_oid != NO_INDEX ) {

      // NOTE: here body/body-pair sub-objects are normal boundaries!
      // They are stored for parents in same ObjectTable variables
      // as normal boundary groups (sbIds, sbTgs) !!!
      //
      int sub_oid = be->Id();
      int sub_tag = be->Tag();

#if 0
// Info only to bodies
      if ( bd1 != NULL && bd1_id != NO_INDEX ) {
        writeIdVariable(interp, "ObjectTable", bd1_id, "sbIds", sub_oid, false);
        writeIdVariable(interp, "ObjectTable", bd1_id, "sbTgs", sub_tag, false);
      }

      if ( bd2 != NULL && bd2_id != NO_INDEX ) {
        writeIdVariable(interp, "ObjectTable", bd2_id,"sbIds", sub_oid, false);
        writeIdVariable(interp, "ObjectTable", bd2_id, "sbTgs", sub_tag, false);
      }
#endif

      // Write subelement info into parent
      writeIdVariable(interp, "ObjectTable", parent_oid, "sbIds", sub_oid, false);
      writeIdVariable(interp, "ObjectTable", parent_oid, "sbTgs", sub_tag, false);

    } // if parent-oid != NO_INDEX

    // Store boundary data
    // -------------------
    int id = be->Id();

    writeVariable(interp, "ObjectTable", "ids", 1, &id, false);

    writeIdVariable(interp, "ObjectTable", id, "tg", be->Tag());

    writeIdVariable(interp, "ObjectTable", id, "grpId", be->getElementGroupId());
    writeIdVariable(interp, "ObjectTable", id, "cl", "");
    writeIdVariable(interp, "ObjectTable", id, "msk", "");
    writeIdVariable(interp, "ObjectTable", id, "nm", be->getName());
    writeIdVariable(interp, "ObjectTable", id, "prId", parent_oid);
    writeIdVariable(interp, "ObjectTable", id, "sbTgs", "");
    writeIdVariable(interp, "ObjectTable", id, "sbIds", "");
    writeIdVariable(interp, "ObjectTable", id, "bndrp", be->getBoundaryParameterId());
    //writeIdVariable(interp, "ObjectTable", id, "bc", be->getBoundaryConditionId());
    writeIdVariable(interp, "ObjectTable", id, "bc", NO_INDEX);
    writeIdVariable(interp, "ObjectTable", id, "slctd", 0);
    writeIdVariable(interp, "ObjectTable", id, "gh", NO_INDEX);
    writeIdVariable(interp, "ObjectTable", id, "ghIds", "");
    writeIdVariable(interp, "ObjectTable", id, "ghMshIndcs", "");

    switch (btype) {

    case OT_FACE:
      writeIdVariable(interp, "ObjectTable", id, "tp", "F");
      if ( be->isIntraLayerBoundary() ) {
        writeIdVariable(interp, "ObjectTable", id, "cl", "ILB");
      }
      break;

    case OT_EDGE:
      writeIdVariable(interp, "ObjectTable", id, "tp", "E");

      if ( model->getDimension() == ECIF_3D ) {

        writeIdVariable(interp, "ObjectTable", id, "bc,vtbl", "");

      } else if ( model->getDimension() == ECIF_2D ) {

        // Init discretization data
        writeIdVariable(interp, "ObjectTable", id, "nofCmp", "");
        writeIdVariable(interp, "ObjectTable", id, "dscTp", "");
        writeIdVariable(interp, "ObjectTable", id, "dscU", "");
        writeIdVariable(interp, "ObjectTable", id, "dscV", "");
        writeIdVariable(interp, "ObjectTable", id, "useFN", "");

        if ( be->isIntraLayerBoundary() ) {
          writeIdVariable(interp, "ObjectTable", id, "cl", "ILB");
        }
      }
      break;

    case OT_VERTEX:
      writeIdVariable(interp, "ObjectTable", id, "tp", "V");
      writeIdVariable(interp, "ObjectTable", id, "bc,vtbl", "");
      //writeIdVariable(interp, "ObjectTable", id, "prIds", "");
      break;
    }

    int nof_grid_hs = be->getNofGridHIds();

    if ( nof_grid_hs > 0 ) {
      writeIdVariable(interp, "ObjectTable", id, "ghIds", nof_grid_hs, be->getGridHIds());
      writeIdVariable(interp, "ObjectTable", id, "ghMshIndcs", nof_grid_hs, be->getGridMeshIndices());
    }

    // Pick subelements ids
    for (int i = 0; i < be->getNofSubElements(); i++) {

      BodyElement* se = be->getSubElement(i);

      // NOTE: this should not happend!
      if ( se == NULL ) {
        continue;
      }

      int se_tag = se->Tag();
      int se_oid = se->Id();

      writeIdVariable(interp, "ObjectTable", id, "sbTgs", se_tag, false);
      writeIdVariable(interp, "ObjectTable", id, "sbIds", se_oid, false);

      if ( btype == OT_EDGE ) {
        writeIdVariable(interp, "ObjectTable", se_oid, "prIds", id, false);
      }
    }
    
    // Set discretization data
    //
    if ( model->getDimension() == ECIF_2D && btype == OT_EDGE ) {

      int cnt = 0;
      linDeltaType* types = NULL;
      double* valuesU = NULL;
      double* valuesV = NULL;
      bool* useFixedN = NULL;

      strstream strm;

      be->getDiscretizationData(cnt, types, valuesU, valuesV, useFixedN);
      
      if ( cnt > 0 ) {
        writeIdVariable(interp, "ObjectTable", id, "nofCmp", cnt);

        for (int i = 0; i < cnt; i++) {

          if ( types != NULL ) {
            switch (types[i] ) {
              case LIN_DELTA_NONE: strm << '-'; break;
              case LIN_DELTA_H: strm << 'H'; break;
              case LIN_DELTA_N: strm << 'N'; break;
              case LIN_DELTA_U: strm << 'U'; break;
            }
          }
          
          if ( valuesU != NULL ) {
            writeIdVariable(interp, "ObjectTable", id, "dscU", valuesU[i], false);
          }

          if ( useFixedN != NULL ) {
            writeIdVariable(interp, "ObjectTable", id, "useFN", (int)useFixedN[i], false);
          }
        }

        if ( types != NULL ) {
          strm << ends;
          writeIdVariable(interp, "ObjectTable", id, "dscTp", strm.str());
        }

      }
    } // Set discretization data

  }
}


// Function writes boundary group related data into Tk variables.
void
UserInterface_TCL::to_tk_WriteElementGroupData(Tcl_Interp* interp, const Model& model)
{
  //Ups! An ugly cast needed owing to model's constant type. !!!***!!!
  Model* mdl = (Model*) &model;

  // Clean old boundary object data
  //
  sendCommandToGui(interp, "Interface::resetObjectTableByType", "FG");
  sendCommandToGui(interp, "Interface::resetObjectTableByType", "EG");
  sendCommandToGui(interp, "Interface::resetObjectTableByType", "VG");

  int index = 0;
  while (true) {

    BodyElementGroup* beg = mdl->getBodyElementGroup(index++);

    if (beg==NULL) break;

    if ( beg->getNofElements() == 0 ) continue;

    // Parent info for the boundary group
    //
    int bd1_id = beg->getParentId(1);
    int bd2_id = beg->getParentId(2);

    int bd1_lr = beg->getParentLayer(1);
    int bd2_lr = beg->getParentLayer(2);

    Body* bd1 = mdl->getBodyById(bd1_id);
    Body* bd2 = mdl->getBodyById(bd2_id);

    BodyPair* bp = NULL;

    if ( bd1 != NULL && bd2 != NULL ) {
      bp = mdl->getBodyPairById((const Body*)bd1, (const Body*)bd2);
    }

    int parent_oid;

    if ( bp != NULL ) {
      parent_oid = bp->Id();

    } else if ( bd1 != NULL ) {
      parent_oid = bd1->Id();

    } else if ( bd2 != NULL ) {
      parent_oid = bd2->Id();

    } else {
      parent_oid = NO_INDEX;
    }

    if ( parent_oid != NO_INDEX ) {

      // NOTE: here body/body-pair sub-objects are boundary groups!
      // They are stored for parents in same ObjectTable variables
      // as normal boundaries (sbIds, sbTgs) !!!
      //
      int sub_oid = beg->Id();
      int sub_tag = beg->Tag();

      // Write subelement info into parent
      writeIdVariable(interp, "ObjectTable", parent_oid, "sbIds", sub_oid, false);
      writeIdVariable(interp, "ObjectTable", parent_oid, "sbTgs", sub_tag, false);
    }

    const char* nm;
    int tg;

    // For an implicit group we use boundary's name and tag
    if ( beg->isImplicit() ) {

      BodyElement* be = (BodyElement*)beg->getElement(0);

      if ( be != NULL ) {
        nm = be->getName();
        tg = be->Tag();
      } else {
        nm = "Unknown";
        tg = NO_INDEX;
      }

    // For an explicit and virtual group we use its own name and tag
    } else {
      nm = beg->getName();
      tg = beg->Tag();
    }

    // Store boundary group data
    // -------------------------
    int id = beg->Id();

    writeVariable(interp, "ObjectTable", "ids", 1, &id, false);

    writeIdVariable(interp, "ObjectTable", id, "tg", tg);

    switch (beg->getElementType()) {
    case OT_FACE:
      writeIdVariable(interp, "ObjectTable", id, "tp", "FG");
      break;
    case OT_EDGE:
      writeIdVariable(interp, "ObjectTable", id, "tp", "EG");
      if ( mdl->getDimension() == ECIF_3D )
        writeIdVariable(interp, "ObjectTable", id, "bc,vtbl", "");
      break;
    case OT_VERTEX:
      writeIdVariable(interp, "ObjectTable", id, "tp", "VG");
      writeIdVariable(interp, "ObjectTable", id, "bc,vtbl", "");
      break;
    }

    if ( beg->isImplicit() )
      writeIdVariable(interp, "ObjectTable", id, "cl", "IMP");
    else if ( beg->isExplicit() )
      writeIdVariable(interp, "ObjectTable", id, "cl", "EXP");
    else if ( beg->isVirtual() )
      writeIdVariable(interp, "ObjectTable", id, "cl", "VIR");
    else
      writeIdVariable(interp, "ObjectTable", id, "cl", "");

    writeIdVariable(interp, "ObjectTable", id, "msk", "");
    writeIdVariable(interp, "ObjectTable", id, "nm", nm);
    writeIdVariable(interp, "ObjectTable", id, "prId", parent_oid);
    writeIdVariable(interp, "ObjectTable", id, "mbrIds", "");
    writeIdVariable(interp, "ObjectTable", id, "sbTgs", "");
    writeIdVariable(interp, "ObjectTable", id, "sbIds", "");
    writeIdVariable(interp, "ObjectTable", id, "bndrp", beg->getBoundaryParameterId());
    writeIdVariable(interp, "ObjectTable", id, "bc", beg->getBoundaryConditionId());
    writeIdVariable(interp, "ObjectTable", id, "slctd", 0);

    // NOTE: We do not show subelements for virtual groups
    //
    // If ex. some vertices should also be available in the bc-panel, the user
    // should make a virtual vertex-group and add this to the same body where
    // this group is
    //
    if ( beg->isVirtual() ) continue;


    // Store boundary and subelement (=subelements boundary group) info
    // -----------------------------

    // This is needed to store only unique sub-id entries
    //
    IdsSet subIds;
    std::pair<IdsSet::iterator, bool> subPair;

    for (int i = 0; i < beg->getNofElements(); i++) {

      BodyElement* be = mdl->getBodyElementById(beg->getElementId(i));

      // Store group's member (boundary) id
      //
      writeIdVariable(interp, "ObjectTable", id, "mbrIds", be->Id(), false);

      for (int j = 0; j < be->getNofSubElements(); j++) {

        BodyElement* se = be->getSubElement(j);

        // Get sub-element boundary group ids
        int sg_oid = se->getElementGroupId();
        int sg_tag = se->getElementGroupTag();

        // Try to insert a new sub-id and check if it is already inserted
        subPair = subIds.insert(sg_oid);

        // Sub-id already stored
        if ( !subPair.second ) continue;

        writeIdVariable(interp, "ObjectTable", id, "sbTgs", sg_tag, false);
        writeIdVariable(interp, "ObjectTable", id, "sbIds", sg_oid, false);
      }
    }

  } // All groups
}


// *** Function sets model status flag for gui.
void
UserInterface_TCL::to_tk_WriteControlParameters(Tcl_Interp* interp, const Model& model)
{
  writeVariable(interp, "Info", "maxNormalTolerance", MAX_NORMAL_TOLERANCE);
  writeVariable(interp, "Info", "maxDistanceTolerance", MAX_DISTANCE_TOLERANCE);
  writeVariable(interp, "Info", "normalTolerance", NORMAL_TOLERANCE);
  writeVariable(interp, "Info", "distanceTolerance", DISTANCE_TOLERANCE);

}


// *** Function sets model level objects/'elements' data for gui.
void
UserInterface_TCL::to_tk_WriteModelData(Tcl_Interp* interp, const Model& model)
{
  int i;
  const ModelInfo* mi = model.getModelInfo();

  char* cval = NULL;

  // Model input version numbers. Needed for possible old version
  // conversion etc.
  writeVariable(theInterp, "Info", "FRONT_PREVIOUS_INPUT_VERSION_NBR", mi->frontPreviousInputVersionNbr);
  writeVariable(theInterp, "Info", "FRONT_INPUT_VERSION_NBR", mi->frontInputVersionNbr);

  // Flag for using values in model
  int use_model_settings = 0;
  readVariable(interp, "UserSetting", "DEFAULT_USE_MODEL_SETTINGS", use_model_settings);

  writeVariable(interp, "ModelProperty", "MODEL_NAME", mi->modelName);
  writeVariable(interp, "ModelProperty", "PROBLEM_NAME", mi->problemName);
  writeVariable(interp, "ModelProperty", "MODEL_DESCRIPTION", mi->modelDescription);
  writeVariable(interp, "ModelProperty", "PROBLEM_DESCRIPTION", mi->problemDescription);

  // Read values from model file
  if (use_model_settings ) {

    cval = mi->includePath;
    if ( cval != NULL && cval[0] != '\0' ) {
      writeVariable(interp, "ModelProperty", "INCLUDE_PATH", cval);
      writeVariable(interp, "ModelProperty", "INCLUDE_PATH,model", cval);
      writeVariable(interp, "ModelProperty", "INCLUDE_PATH,model,save", 1);
    }

    cval = mi->resultsDirectory;
    if ( cval != NULL && cval[0] != '\0' ) {
      writeVariable(interp, "ModelProperty", "RESULTS_DIRECTORY", cval);
      writeVariable(interp, "ModelProperty", "RESULTS_DIRECTORY,model", cval);
      writeVariable(interp, "ModelProperty", "RESULTS_DIRECTORY,model,save", 1);
    }

    cval = mi->temporaryFilesDirectory;
    if ( cval != NULL && cval[0] != '\0' ) {
      writeVariable(interp, "ModelProperty", "LOG_DIRECTORY", cval);
      writeVariable(interp, "ModelProperty", "LOG_DIRECTORY.model", cval);
      writeVariable(interp, "ModelProperty", "LOG_DIRECTORY,model,save", 1);
    }
  }

  // History
  writeVariable(interp, "Model", "created", mi->created);
  writeVariable(interp, "Model", "modified", mi->modified);
  writeVariable(interp, "Model", "hasUserDefinitions", mi->hasUserDefinitions);

  // Source files
  writeVariable(interp, "Model", "cadPath", mi->cadSourceFile);
  writeVariable(interp, "Model", "meshPath", mi->meshSourceFile);

  // Matc input files
  writeVariable(interp, "Model", "matcInputFile_emf", mi->matcInputFile_emf);
  writeVariable(interp, "Model", "matcInputFile_sif", mi->matcInputFile_sif);

  // Timestamps
  writeVariable(interp, "Model", "Database,timestamp", mi->databaseTs);
  writeVariable(interp, "Model", "GebhardtFactors,timestamp", mi->gebhardtFactorsTs);
  writeVariable(interp, "Model", "MeshParameter,timestamp", mi->meshParameterTs);
  writeVariable(interp, "Model", "Mesh,timestamp", mi->meshTs);
  writeVariable(interp, "Model", "Solver,timestamp", mi->solverTs);
  writeVariable(interp, "Model", "Viewfactors,timestamp", mi->viewfactorsTs);

  // Mesh names and control
  // ======================
  writeVariable(interp, "Model", "meshNames", "");
  writeVariable(interp, "Model", "meshHs", "");
  writeVariable(interp, "Model", "meshFs", "");
  writeVariable(interp, "Model", "currentMeshIndex", mi->currentMeshIndex);

  int nof_meshes = mi->nofMeshes;
  for (i = 0; i < nof_meshes; i++) {
    writeVariable(interp, "Model", "meshNames", mi->meshNames[i], false);

    if ( mi->meshHs != NULL ) {
      writeVariable(interp, "Model", "meshHs", mi->meshHs[i], false);
      writeVariable(interp, "Model", "meshFs", mi->meshFs[i], false);
    }
  }

  // Bg mesh files
  // =============
  writeVariable(interp, "Model", "meshBgMeshFileIndices", "");
  writeVariable(interp, "Model", "meshBgMeshFiles", "");
  writeVariable(interp, "Model", "meshBgMeshActives", "");
  writeVariable(interp, "Model", "meshBgMeshControls", "");

  int nof_files = mi->nofBgMeshFiles;
  for (i = 0; i < nof_files; i++) {
    writeVariable(interp, "Model", "meshBgMeshFileIndices", mi->meshBgMeshFileIndices[i], false);
    writeVariable(interp, "Model", "meshBgMeshFiles", mi->meshBgMeshFiles[i], false);
    writeVariable(interp, "Model", "meshBgMeshActives", mi->meshBgMeshActives[i], false);
    writeVariable(interp, "Model", "meshBgMeshControls", mi->meshBgMeshControls[i], false);
  }

  // Apply data in Gui
  sendCommandToGui(interp, "Interface::applyModelData");
}


// *** Function sets model status flag for gui.
void
UserInterface_TCL::to_tk_WriteModelStatus(Tcl_Interp* interp, const Model& model)
{
  writeVariable(interp, "Model", "status", model.getModelStatus());
  sendCommandToGui(interp, "Interface::applyModelStatus");
}


// *** Function writes model status message for gui.
void
UserInterface_TCL::to_tk_WriteModelStatusMessage(Tcl_Interp* interp, const Model& model)
{
  ostrstream strm;
  model.getModelStatusMessage(strm);
  strm << ends;

  writeVariable(interp, "Model", "statusMessage", strm.str());
}


// *** Function sets model geometry-dimesnion (2D/3D) variable for gui.
void
UserInterface_TCL::to_tk_WriteModelGeometryDimension(Tcl_Interp* interp, const Model& model)
{
  if ( model.getDimension() == ECIF_2D)
    writeVariable(interp, "Model", "GEOMETRY_DIMENSION", "2D");
  else
    writeVariable(interp, "Model", "GEOMETRY_DIMENSION", "3D");

  sendCommandToGui(interp, "Interface::applyModelGeometryDimension");
}


// *** Function sets model flags for gui.
void
UserInterface_TCL::to_tk_WriteModelFlags(Tcl_Interp* interp, const Model& model)
{
  Model* mdl = (Model*)&model;

  writeVariable(interp, "ModelFlags", "GEOMETRY_TYPE_CAD", mdl->getFlagValue(GEOMETRY_TYPE_CAD));
  writeVariable(interp, "ModelFlags", "GEOMETRY_TYPE_MESH", mdl->getFlagValue(GEOMETRY_TYPE_MESH));

  writeVariable(interp, "ModelFlags", "GEOMETRY_EDITED_BODIES", mdl->getFlagValue(GEOMETRY_EDITED_BODIES));
  writeVariable(interp, "ModelFlags", "GEOMETRY_EDITED_BOUNDARIES", mdl->getFlagValue(GEOMETRY_EDITED_BOUNDARIES));

  writeVariable(interp, "ModelFlags", "DRAW_SOURCE_CAD", mdl->getFlagValue(DRAW_SOURCE_CAD));
  writeVariable(interp, "ModelFlags", "DRAW_SOURCE_MESH", mdl->getFlagValue(DRAW_SOURCE_MESH));
  writeVariable(interp, "ModelFlags", "DRAW_TARGET_BODIES", mdl->getFlagValue(DRAW_TARGET_BODIES));
  writeVariable(interp, "ModelFlags", "DRAW_TARGET_SURFACES", mdl->getFlagValue(DRAW_TARGET_SURFACES));
  writeVariable(interp, "ModelFlags", "DRAW_TARGET_EDGES", mdl->getFlagValue(DRAW_TARGET_EDGES));

  writeVariable(interp, "ModelFlags", "LABEL_DISPLAY_NODE", mdl->getFlagValue(LABEL_DISPLAY_NODE));
  writeVariable(interp, "ModelFlags", "LABEL_DISPLAY_ELEMENT", mdl->getFlagValue(LABEL_DISPLAY_ELEMENT));
  writeVariable(interp, "ModelFlags", "LABEL_DISPLAY_VERTEX", mdl->getFlagValue(LABEL_DISPLAY_VERTEX));
  writeVariable(interp, "ModelFlags", "LABEL_DISPLAY_EDGE", mdl->getFlagValue(LABEL_DISPLAY_EDGE));
  writeVariable(interp, "ModelFlags", "LABEL_DISPLAY_FACE", mdl->getFlagValue(LABEL_DISPLAY_FACE));
  writeVariable(interp, "ModelFlags", "LABEL_DISPLAY_BODY", mdl->getFlagValue(LABEL_DISPLAY_BODY));

  sendCommandToGui(interp, "Interface::applyModelFlags");
}


// *** Function sets values for Tcl-variables in the Statistic-panel.
void
UserInterface_TCL::to_tk_WriteProcessorData(Tcl_Interp* interp, const Model& model)
{
  const ParallelInfo* pi = model.getParallelInfo();

  ostrstream strm;
  strm << pi->nofProcessors;
  strm << ends;
  sendCommandToGui(interp, "Interface::getParallelInfo", strm.str());
}


// *** Function sets values for Tcl-variables in the Statistic-panel.
void
UserInterface_TCL::to_tk_WriteStats(Tcl_Interp* interp, const Model& model)
{
  const ModelStatistics* mstats = model.getModelStatistics();
  const ModelInfo* minfo = model.getModelInfo();
  const MeshInfo* mesh_info = model.getMeshInfo();

  writeVariable(interp, "Model", "nofBodies", mstats->nofBodies);
  writeVariable(interp, "Model", "nofElements", mstats->nofOuterBoundaries + mstats->nofInnerBoundaries);
  writeVariable(interp, "Model", "nofOuterBoundaries", mstats->nofOuterBoundaries);
  writeVariable(interp, "Model", "nofInnerBoundaries", mstats->nofInnerBoundaries);
  writeVariable(interp, "Model", "nofVertices", mstats->nofVertices);
  writeVariable(interp, "Model", "minEdgeSize", minfo->minEdgeSize);

  // Reasonable cad geometry
  if (minfo->minX < minfo->maxX) {
    writeVariable(interp, "Model", "minX", minfo->minX);
    writeVariable(interp, "Model", "maxX", minfo->maxX);
    writeVariable(interp, "Model", "minY", minfo->minY);
    writeVariable(interp, "Model", "maxY", minfo->maxY);
    writeVariable(interp, "Model", "minZ", minfo->minZ);
    writeVariable(interp, "Model", "maxZ", minfo->maxZ);
  // Mesh geometry
  } else {
    writeVariable(interp, "Model", "minX", mesh_info->minX);
    writeVariable(interp, "Model", "maxX", mesh_info->maxX);
    writeVariable(interp, "Model", "minY", mesh_info->minY);
    writeVariable(interp, "Model", "maxY", mesh_info->maxY);
    writeVariable(interp, "Model", "minZ", mesh_info->minZ);
    writeVariable(interp, "Model", "maxZ", mesh_info->maxZ);
  }

  int nof_zv_elements = mesh_info->nofZeroVelocityElements;
  int nof_sp_elements = mesh_info->nofSplittedElements;

  writeVariable(interp, "Model", "nofMeshElements", mesh_info->nofBulkElements);
  writeVariable(interp, "Model", "nofMeshSplittedElements", nof_sp_elements);
  writeVariable(interp, "Model", "nofMeshZeroVelocityElements", nof_zv_elements);
  writeVariable(interp, "Model", "nofMeshBndrElements", mesh_info->nofBoundaryElements);
  writeVariable(interp, "Model", "nofMeshInnerBndrElements", mesh_info->nofInnerBndrElements);
  writeVariable(interp, "Model", "nofMeshOuterBndrElements", mesh_info->nofOuterBndrElements);
  writeVariable(interp, "Model", "nofMeshNodes", mesh_info->nofNodes);

  sendCommandToGui(interp, "Interface::applyStatsData");
}


// *** Function sets model's bodyforce parameters' status field for gui.
void
UserInterface_TCL::to_tk_WriteStatusBodyForces(Tcl_Interp* interp, const Model& model)
{
  const ModelStatistics* st = model.getModelStatistics();

  int count = st->nofBodiesWithBodyForce;
  writeVariable(interp, "Model", "nofBodiesWithBodyForce", count);
  sendCommandToGui(interp, "Interface::setStatusBodyForces");
}


// *** Function sets model's material parameters' status field for gui.
void
UserInterface_TCL::to_tk_WriteStatusBoundaryConditions(Tcl_Interp* interp, const Model& model)
{
  const ModelStatistics* st = model.getModelStatistics();

  int count;

  count = st->nofOuterBoundariesWithCondition;
  writeVariable(interp, "Model", "nofOuterBoundariesWithCondition", count);
  sendCommandToGui(interp, "Interface::setStatusOuterBoundaryConditions");

  count = st->nofInnerBoundariesWithCondition;
  writeVariable(interp, "Model", "nofInnerBoundariesWithCondition", count);
  sendCommandToGui(interp, "Interface::setStatusInnerBoundaryConditions");

}


// *** Function sets model's equation parameters' status field for gui.
void
UserInterface_TCL::to_tk_WriteStatusEquations(Tcl_Interp* interp, const Model& model)
{
  const ModelStatistics* st = model.getModelStatistics();

  int count =  st->nofBodiesWithEquation;
  writeVariable(interp, "Model", "nofBodiesWithEquation", count);
  sendCommandToGui(interp, "Interface::setStatusEquations");
}


// *** Function sets model's initial condition parameters' status field for gui.
void
UserInterface_TCL::to_tk_WriteStatusInitialConditions(Tcl_Interp* interp, const Model& model)
{
  const ModelStatistics* st = model.getModelStatistics();


  int count =  st->nofBodiesWithInitialCondition;
  writeVariable(interp, "Model", "nofBodiesWithInitialCondition", count);
  sendCommandToGui(interp, "Interface::setStatusInitialConditions");
}


// *** Function sets model's material parameters' status field for gui.
void
UserInterface_TCL::to_tk_WriteStatusMaterials(Tcl_Interp* interp, const Model& model)
{
  const ModelStatistics* st = model.getModelStatistics();

  int count = st->nofBodiesWithMaterial;
  writeVariable(interp, "Model", "nofBodiesWithMaterial", count);
  sendCommandToGui(interp, "Interface::setStatusMaterials");
}


// *** Function sets model's meshes status field for gui.
void
UserInterface_TCL::to_tk_WriteStatusMeshes(Tcl_Interp* interp, const Model& model)
{
  sendCommandToGui(interp, "Interface::setStatusMeshes");
}


// *** Function sets model's timestep parameter's status field for gui.
void
UserInterface_TCL::to_tk_WriteStatusTimesteps(Tcl_Interp* interp, const Model& model)
{
  Model* mdl = (Model*) &model;

  int count = mdl->getNofTimestepSteps();
  writeVariable(interp, "Model", "nofTimesteps", count);
  sendCommandToGui(interp, "Interface::setStatusTimestamps");
}


// Convert variable's gui-name (like "oxygen") to sif-name like ("OxyGen")
void
UserInterface_TCL::variableNameGuiToSif(const char* gui_name, char* sif_name_buffer)
{
  sendCommandToGui(theInterp, "Interface::variableNameGuiToSif", gui_name);
  strcpy(sif_name_buffer, getCommandResults(theInterp));
}


// Convert variables's sif-name (like "Some Quantity") to gui-name ("some_quantity")
void
UserInterface_TCL::variableNameSifToGui(const char* sif_name, char* gui_name_buffer)
{
  sendCommandToGui(theInterp, "Interface::variableNameSifToGui", sif_name);
  strcpy(gui_name_buffer, getCommandResults(theInterp));
}


// Id Array, integer id, single value
// ----------------------------------
void
UserInterface_TCL::writeIdVariable(Tcl_Interp* interp, const char* array, int id,
                                       const char* variable,
                                       const char* value,
                                       bool reset)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  writeVariable_impl_1(interp, array, field_name.str(), value, reset);
}

void
UserInterface_TCL::writeIdVariable(Tcl_Interp* interp, const char* array, int id,
                                       const char* variable,
                                       int value,
                                       bool reset)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  writeVariable_impl_1(interp, array, field_name.str(), value, reset);
}

void
UserInterface_TCL::writeIdVariable(Tcl_Interp* interp, const char* array, int id,
                                       const char* variable,
                                       ProcessId value,
                                       bool reset)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  writeVariable_impl_1(interp, array, field_name.str(), value, reset);
}

void
UserInterface_TCL::writeIdVariable(Tcl_Interp* interp, const char* array, int id,
                                       const char* variable,
                                       double value,
                                       bool reset)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  writeVariable_impl_1(interp, array, field_name.str(), value, reset);
}


// Id Array, integer id, list of values
// ------------------------------------

void
UserInterface_TCL::writeIdVariable(Tcl_Interp* interp, const char* array, int id,
                                       const char* variable,
                                       int size, const char** values,
                                       bool reset)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  writeVariable_impl_n(interp, array, field_name.str(), size, values, reset);
}

void
UserInterface_TCL::writeIdVariable(Tcl_Interp* interp, const char* array, int id,
                                       const char* variable,
                                       int size, const int* values,
                                       bool reset)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  writeVariable_impl_n(interp, array, field_name.str(), size, values, reset);
}

void
UserInterface_TCL::writeIdVariable(Tcl_Interp* interp, const char* array, int id,
                                       const char* variable,
                                       int size, const double* values,
                                       bool reset)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  writeVariable_impl_n(interp, array, field_name.str(), size, values, reset);
}


// Id Array, string id, single value
// ---------------------------------
void
UserInterface_TCL::writeIdVariable(Tcl_Interp* interp, const char* array, const char* id,
                                       const char* variable,
                                       const char* value,
                                       bool reset)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  writeVariable_impl_1(interp, array, field_name.str(), value, reset);
}

void
UserInterface_TCL::writeIdVariable(Tcl_Interp* interp, const char* array, const char* id,
                                       const char* variable,
                                       int value,
                                       bool reset)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  writeVariable_impl_1(interp, array, field_name.str(), value, reset);
}

void
UserInterface_TCL::writeIdVariable(Tcl_Interp* interp, const char* array, const char* id,
                                       const char* variable,
                                       ProcessId value,
                                       bool reset)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  writeVariable_impl_1(interp, array, field_name.str(), value, reset);
}

void
UserInterface_TCL::writeIdVariable(Tcl_Interp* interp, const char* array, const char* id,
                                       const char* variable,
                                       double value,
                                       bool reset)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  writeVariable_impl_1(interp, array, field_name.str(), value, reset);
}


// Id Array, string id, list of values
// -----------------------------------

void
UserInterface_TCL::writeIdVariable(Tcl_Interp* interp, const char* array, const char* id,
                                       const char* variable,
                                       int size, const char** values,
                                       bool reset)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  writeVariable_impl_n(interp, array, field_name.str(), size, values, reset);
}

void
UserInterface_TCL::writeIdVariable(Tcl_Interp* interp, const char* array, const char* id,
                                       const char* variable,
                                       int size, const int* values,
                                       bool reset)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  writeVariable_impl_n(interp, array, field_name.str(), size, values, reset);
}

void
UserInterface_TCL::writeIdVariable(Tcl_Interp* interp, const char* array, const char* id,
                                       const char* variable,
                                       int size, const double* values,
                                       bool reset)
{
  strstream field_name;
  field_name << id << "," << variable << ends;
  writeVariable_impl_n(interp, array, field_name.str(), size, values, reset);
}


// Id2 Array, integer id1,id2, single value
// ----------------------------------------
void
UserInterface_TCL::writeId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2,
                                       const char* variable,
                                       const char* value,
                                       bool reset)
{
  strstream field_name;
  field_name << id1 << subIdSep << id2 << "," << variable << ends;
  writeVariable_impl_1(interp, array, field_name.str(), value, reset);
}

void
UserInterface_TCL::writeId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2,
                                       const char* variable,
                                       int value,
                                       bool reset)
{
  strstream field_name;
  field_name << id1 << subIdSep << id2 << "," << variable << ends;
  writeVariable_impl_1(interp, array, field_name.str(), value, reset);
}

void
UserInterface_TCL::writeId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2,
                                       const char* variable,
                                       ProcessId value,
                                       bool reset)
{
  strstream field_name;
  field_name << id1 << subIdSep << id2 << "," << variable << ends;
  writeVariable_impl_1(interp, array, field_name.str(), value, reset);
}

void
UserInterface_TCL::writeId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2,
                                       const char* variable,
                                       double value,
                                       bool reset)
{
  strstream field_name;
  field_name << id1 << subIdSep << id2 << "," << variable << ends;
  writeVariable_impl_1(interp, array, field_name.str(), value, reset);
}


// Id2 Array, integer id1, id2, list of values
// -------------------------------------------
void
UserInterface_TCL::writeId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2,
                                       const char* variable,
                                       int size, const char** values,
                                       bool reset)
{
  strstream field_name;
  field_name << id1 << subIdSep << id2 << "," << variable << ends;
  writeVariable_impl_n(interp, array, field_name.str(), size, values, reset);
}

void
UserInterface_TCL::writeId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2,
                                       const char* variable,
                                       int size, const int* values,
                                       bool reset)
{
  strstream field_name;
  field_name << id1 << subIdSep << id2 << "," << variable << ends;
  writeVariable_impl_n(interp, array, field_name.str(), size, values, reset);
}

void
UserInterface_TCL::writeId2Variable(Tcl_Interp* interp, const char* array, int id1, int id2,
                                       const char* variable,
                                       int size, const double* values,
                                       bool reset)
{
  strstream field_name;
  field_name << id1 << subIdSep << id2 << "," << variable << ends;
  writeVariable_impl_n(interp, array, field_name.str(), size, values, reset);
}


// Array, single value
// -------------------
void
UserInterface_TCL::writeVariable(Tcl_Interp* interp,
                                       const char* array, const char* variable,
                                       const char* value,
                                       bool reset)
{
  writeVariable_impl_1(interp, array, variable, value, reset);
}

void
UserInterface_TCL::writeVariable(Tcl_Interp* interp,
                                       const char* array, const char* variable,
                                       int value,
                                       bool reset)
{
  writeVariable_impl_1(interp, array, variable, value, reset);
}

void
UserInterface_TCL::writeVariable(Tcl_Interp* interp,
                                       const char* array, const char* variable,
                                       double value,
                                       bool reset)
{
  writeVariable_impl_1(interp, array, variable, value, reset);
}



// Array, list of values
// ---------------------
void
UserInterface_TCL::writeVariable(Tcl_Interp* interp,
                                       const char* array, const char* variable,
                                       int size, const char** values,
                                       bool reset)
{
  writeVariable_impl_n(interp, array, variable, size, values, reset);
}

void
UserInterface_TCL::writeVariable(Tcl_Interp* interp,
                                       const char* array, const char* variable,
                                       int size, const int* values,
                                       bool reset)
{
  writeVariable_impl_n(interp, array, variable, size, values, reset);
}

void
UserInterface_TCL::writeVariable(Tcl_Interp* interp,
                                       const char* array, const char* variable,
                                       int size, const double* values,
                                       bool reset)
{
  writeVariable_impl_n(interp, array, variable, size, values, reset);
}


// *** Function sets single value to array-type Tcl variable.
template <class T > void
UserInterface_TCL::writeVariable_impl_1(Tcl_Interp* interp,
                                       const char* array, const char* variable,
                                       const T value, bool reset)
{
  Tcl_Obj* array_obj = Tcl_NewStringObj(NULL, 0);
  Tcl_Obj* value_obj = NULL;

  int tcl_flags = TCL_GLOBAL_ONLY | TCL_LEAVE_ERR_MSG;

  // We have to create a tcl string-object to
  // carry the variable names. Here we use a fully
  // specified array name to create the name-object
  tcl_flags |= TCL_PARSE_PART1;
  strstream strm;
  strm << array << "(" << variable << ")" << ends;
  char* array_name = strm.str();
  Tcl_SetStringObj(array_obj, array_name, strlen(array_name));

  // Try to get existing variable handle
  // NOTE! This does not seem to be working!!! At least not for
  // strings like Default-directories in the settings files.
  // They all get the last string values if this is used
  // Try to get existing variable handle
  ;//value_obj = Tcl_ObjGetVar2(interp, array_obj, NULL, tcl_flags);

  Tcl_ResetResult(interp);

  setTclObjValue(value_obj, value);

  if (!reset) {
    // Turn on the list append mode
    tcl_flags |= TCL_LIST_ELEMENT;
    tcl_flags |= TCL_APPEND_VALUE;
  }

  // Update variable
  Tcl_ObjSetVar2(interp, array_obj, NULL, value_obj, tcl_flags);

}


// *** Function sets list of values to array-type Tcl variable.
template <class T > void
UserInterface_TCL::writeVariable_impl_n(Tcl_Interp* interp,
                                       const char* array, const char* variable,
                                       int size, const T* values,
                                       bool reset)
{
  Tcl_Obj* array_obj = Tcl_NewStringObj(NULL, 0);
  Tcl_Obj* value_obj = NULL;

  int tcl_flags = TCL_GLOBAL_ONLY | TCL_LEAVE_ERR_MSG;

  tcl_flags |= TCL_PARSE_PART1;
  // In Tcl80 we have to create a tcl string-object to
  // carry the variable names. Here we use a fully
  // specified array name to create the name-object
  strstream strm;
  strm << array << "(" << variable << ")" << ends;
  char* array_name = strm.str();
  Tcl_SetStringObj(array_obj, array_name, strlen(array_name));

  // Reset variable
  if (reset) {
    setTclObjValue(value_obj, "");
    Tcl_ObjSetVar2(interp, array_obj, NULL, value_obj, tcl_flags);
  }

  // Turn on the list append mode
  tcl_flags |= TCL_LIST_ELEMENT;
  tcl_flags |= TCL_APPEND_VALUE;

  for (int i = 0; i < size; i++) {

    value_obj = NULL;

    // Try to get existing variable handle
    // NOTE: This apparently does not work!!!
    // Ref. comments for corresponding scalar method "..._impl_1"
    //value_obj = Tcl_ObjGetVar2(interp, array_obj, NULL, tcl_flags);

    Tcl_ResetResult(interp);

    setTclObjValue(value_obj, values[i]);

    // Update variable
    Tcl_ObjSetVar2(interp, array_obj, NULL, value_obj, tcl_flags);

  }
}



void
UserInterface_TCL::getTclObjValue(Tcl_Interp* interp, Tcl_Obj* value_obj, bool& value)
{
  if (value_obj == NULL) {
    return;
  }

  Tcl_GetBooleanFromObj(interp, value_obj, (int*)&value);
}


void
UserInterface_TCL::getTclObjValue(Tcl_Interp* interp, Tcl_Obj* value_obj, int& value)
{
  if (value_obj == NULL) {
    return;
  }

  Tcl_GetIntFromObj(interp, value_obj, &value);
}


void
UserInterface_TCL::getTclObjValue(Tcl_Interp* interp, Tcl_Obj* value_obj, long& value)
{
  if (value_obj == NULL) {
    return;
  }

  Tcl_GetLongFromObj(interp, value_obj, &value);
}


void
UserInterface_TCL::getTclObjValue(Tcl_Interp* interp, Tcl_Obj* value_obj, double& value)
{
  if (value_obj == NULL) {
    return;
  }

  Tcl_GetDoubleFromObj(interp, value_obj, &value);
}


void
UserInterface_TCL::getTclObjValue(Tcl_Interp* interp, Tcl_Obj* value_obj, char& value)
{
  if (value_obj == NULL) {
    return;
  }

  int len = 0;
  char* tmp = Tcl_GetStringFromObj(value_obj, &len);

  if (len > 0) {
    value = tmp[0];
  }

}


void
UserInterface_TCL::getTclObjValue(Tcl_Interp* interp, Tcl_Obj* value_obj, char*& value)
{
  if (value_obj == NULL) {
    return;
  }

  int len = 0;
  char* tmp = Tcl_GetStringFromObj(value_obj, &len);

  if (len > 0) {
    value = new char[1 + len];
    strcpy(value, tmp);
  }
}



void
UserInterface_TCL::setTclObjValue(Tcl_Obj*& value_obj, const bool value)
{
  // If this is new variable, we have to
  // create a Tcl data-object
  if (value_obj == NULL) {
    value_obj = Tcl_NewBooleanObj(0);
  }

  // Set new value into data-object
  Tcl_SetBooleanObj(value_obj, value);
}

void
UserInterface_TCL::setTclObjValue(Tcl_Obj*& value_obj, const int value)
{
  // If this is new variable, we have to
  // create a Tcl data-object
  if (value_obj == NULL) {
    value_obj = Tcl_NewIntObj(0);
  }

  // Set new value into data-object
  Tcl_SetIntObj(value_obj, value);
}

void
UserInterface_TCL::setTclObjValue(Tcl_Obj*& value_obj, const ProcessId value)
{
  // If this is new variable, we have to
  // create a Tcl data-object
  if (value_obj == NULL) {
    value_obj = Tcl_NewIntObj(0);
  }

  // Set new value into data-object
  Tcl_SetIntObj(value_obj, value);
}

void
UserInterface_TCL::setTclObjValue(Tcl_Obj*& value_obj, const long value)
{
  // If this is new variable, we have to
  // create a Tcl data-object
  if (value_obj == NULL) {
    value_obj = Tcl_NewLongObj(0);
  }

  // Set new value into data-object
  Tcl_SetLongObj(value_obj, value);
}

void
UserInterface_TCL::setTclObjValue(Tcl_Obj*& value_obj, const double value)
{
  // If this is new variable, we have to
  // create a Tcl data-object
  if (value_obj == NULL) {
    value_obj = Tcl_NewDoubleObj(0);
  }

  // Set new value into data-object
  Tcl_SetDoubleObj(value_obj, value);
}


void
UserInterface_TCL::setTclObjValue(Tcl_Obj*& value_obj, const char value)
{
  char tmp[1];
  tmp[0] = value;

  // If this is new variable, we have to
  // create a Tcl data-object
  if (value_obj == NULL) {
    value_obj = Tcl_NewStringObj(tmp, 1);
  }

  // Set new value into data-object
  Tcl_SetStringObj(value_obj, (char*)tmp, 1);
}

void
UserInterface_TCL::setTclObjValue(Tcl_Obj*& value_obj, const char* value)
{
  char* val;
  int len;

  if ( value == NULL ) {
    val = "";
    len = 0;
  } else {
    val = (char*)value;
    len = strlen(val);
  }

  // If this is new variable, we have to
  // create a Tcl data-object
  if (value_obj == NULL) {
    value_obj = Tcl_NewStringObj(val, len);
  }

  // Set new value into data-object
  Tcl_SetStringObj(value_obj, val, len);
}




int
UserInterface_TCL::unknownFieldMsg(emf_ObjectData_X* object_data, bool is_fatal)
{
  UserInterface* gui = theControlCenter->getGui();
  strstream strm;
  strm << "WARNING: Unknown field name for object "
       << object_data->object_name
       << ": "
       << object_data->field_name
       << ends;

  gui->showMsg(strm.str());

  if (is_fatal)
    return !isOK;
  else
    return isOK;
}


void
UserInterface_TCL::update()
{
  Tcl_Eval(theInterp, "update idletasks" );
}


// Class method
void
UserInterface_TCL::update(Tcl_Interp* interp)
{
  Tcl_Eval(interp, "update idletasks" );
}


void
UserInterface_TCL::update(int counter, int update_interval)
{
  if ( update_interval > 0 && 0 == counter % update_interval ) {
    update();
  }
}


// Update body data after split/comobine mesh boundaries
void
UserInterface_TCL::updateBodyData(Model* model)
{
  to_tk_WriteBodyData(theInterp, *model);
  to_tk_WriteBodyInfoData(theInterp, *model);

  to_tk_WriteStats(theInterp, *model);
}


void
UserInterface_TCL::updateBoundaryData(Model* model)
{
  to_tk_WriteBoundaryData(theInterp, *model);
  to_tk_WriteElementGroupData(theInterp, *model);

  to_tk_WriteStatusBoundaryConditions(theInterp, *model);

  to_tk_WriteStats(theInterp, *model);

}


void
UserInterface_TCL::updateMeshZeroVelocityElements(int nof_zv_elements)
{
  writeVariable(theInterp, "Model", "nofMeshZeroVelocityElements", nof_zv_elements);
}


// *** Update model data for the gui-side
void
UserInterface_TCL::updateModelData(Model* model)
{
  theModel = model;

  to_tk_WriteModelData(theInterp, *model);

  // Set checked model file dir+name
  sendCommandToGui(theInterp,"Interface::setModelFilePath");
}


void
UserInterface_TCL::updateModelFlags(Model* model)
{
  to_tk_WriteModelFlags(theInterp, *model);
}


// *** Update model statistics for the gui-side
void
UserInterface_TCL::updateModelStatistics(Model* model)
{
  to_tk_WriteStats(theInterp, *model);
  to_tk_WriteBodyMeshInfoData(theInterp, *model);
}


// *** Update model status for the gui-side
void
UserInterface_TCL::updateModelStatus(Model* model)
{
  to_tk_WriteModelStatus(theInterp, *model);
  to_tk_WriteModelGeometryDimension(theInterp, *model);
  to_tk_WriteProcessorData(theInterp, *model);
}


void
UserInterface_TCL::updateNextActiveSelectionTolerance(double tolerance)
{
  ostrstream strm;
  strm <<tolerance << ends;
  sendCommandToGui(theInterp, "Interface::getNextActiveSelectionTolerance", strm.str());
}


// *** Update object data for the gui-side
void
UserInterface_TCL::updateObjectData(Model* model)
{
  theModel = model;

  writeVariable(theInterp, "ObjectTable", "ids", "");
  writeVariable(theInterp, "ObjectTable", "count", 0);

  to_tk_WriteBodyLayerData(theInterp, *model);

  to_tk_WriteBodyData(theInterp, *model);
  to_tk_WriteBodyInfoData(theInterp, *model);

  // Do this after writing Body level data!!!
  to_tk_WriteBoundaryData(theInterp, *model);
  to_tk_WriteElementGroupData(theInterp, *model);

  to_tk_WriteStats(theInterp, *model);

  to_tk_WriteControlParameters(theInterp, *model);

  // Process object table data in Gui side
  sendCommandToGui(theInterp, "Interface::applyObjectTableData");

  sendCommandToGui(theInterp,"Interface::createBodiesMenus");
}


// *** Update parameter data for the gui-side
void
UserInterface_TCL::updateParameterDataPre(Model* model)
{
  to_tk_WriteAllParamDataPre(theInterp, *model);
}


// *** Update object and mask related parameter data for the gui-side
void
UserInterface_TCL::updateParameterDataPost(Model* model)
{
  to_tk_WriteAllParamDataPost(theInterp, *model);
}




void
UserInterface_TCL::updateRendererInfo(const RendererInfo& rinfo)
{
  int* i_val = NULL;
  double* d_val = NULL;
  char* s_val = NULL;

  char* arr = "RendererInfo";

  d_val = (double*) &rinfo.LINE_WIDTH_GRANULARITY;
  writeVariable(theInterp, arr, "LINE_WIDTH_GRANULARITY", 1, d_val);

  d_val = (double*) rinfo.LINE_WIDTH_RANGE;
  writeVariable(theInterp, arr, "LINE_WIDTH_RANGE", 2, d_val);
}


// ********************************************
// ********* Tcl function definitions *********
// ********************************************

// *** Tcl_AppInit
int My_Tcl_AppInit(Tcl_Interp* interp)
{
  if (Tcl_Init(interp) == TCL_ERROR) {
    return TCL_ERROR;
  }

  //NOTE We do not start Tk here. It is started
  // in the start() method just before menus etc. are
  // ready to be built. This way we can (practically)
  // make the tk initial console window invisble

  //int Tklib_Init(Tcl_Interp *interp);
  //Tcl_StaticPackage(interp, "Tklib", Tklib_Init, NULL);

  //if (Tk_Init(interp) == TCL_ERROR) {
  //  return TCL_ERROR;
  //}

#ifdef TK_TEST
  //if (TkTest_Init(interp) == TCL_ERROR) {
  //  return TCL_ERROR;
  //}
#endif /* TK_TEST */

    return TCL_OK;
}


#if 0
// This stuff is a try to controlling
// Xwindows resource usage in the window loop
#ifdef UNIX
#include <signal.h>
#include <sys/time.h>
#include <unistd.h>
// Timer function to control Xwindows looping in Unix
static void timer_dummy() {}
static void timer_sig(int sig)
{
  struct itimerval val;

  int i;

  val.it_interval.tv_sec  = 0;
  val.it_interval.tv_usec = 0;

	val.it_value.tv_sec     = 0;
  val.it_value.tv_usec    = 250000;

  setitimer( ITIMER_REAL, &val, NULL );
  signal( SIGALRM, timer_dummy );
}
#endif
#endif


// *** Main Window-program loop (Tk_MainLoop)
void
Tcl_MainLoop()
{
  // This is needed at least for Unix, to prevent
  // continuos cpu-cycle consumption
  //#ifdef UNIX
  //  Tcl_CreateTimerHandler(1, tcl_DisplayIdleProc, displayIdleData);
  //#else
  //  Tcl_CreateTimerHandler(1, tcl_DisplayIdleProc, displayIdleData);
  //#endif

#if !defined(WIN32)
  Tcl_CreateTimerHandler(1, tcl_DisplayIdleProc, displayIdleData);
#endif

  int event_types = TCL_ALL_EVENTS;
  //int event_types = TCL_WINDOW_EVENTS|TCL_FILE_EVENTS|TCL_TIMER_EVENTS;
  //int event_types = 0;

  // Tcl event loop
  while ( Tcl_DoOneEvent(event_types) );
}


// Function simulates an interrupt with C++ exeception mechanism
// NOTE: This is called eg. from UserInterface_TCL::from_tk_ThrowException
void
tcl_interrupt(ClientData data)
{
#if defined(ELMER_FRONT_EXCEPTIONS)

  // NOTE: You cannot throw execption within Gui, it will
  // crash tcl-interperter sooner or later!!!
  //throw E_ReadE;

#else
  // Do nothing
  TheUI->showMsg("Break not activated, sorry!");
#endif
}


// Tcl proc launched via tcl-timer
void
tcl_DisplayIdleProc(ClientData data)
{
  //NOTE: In WIN32 we do not need a separate call to some window-proc
#if defined(UNIX)
  Renderer* renderer = TheUI->getRenderer();

  if (renderer != NULL) {
    renderer->dummyWindowProc();
  }
#endif

  TheUI->update();

  // Recreate tcl-timer so that I will be called again
  Tcl_CreateTimerHandler(1, tcl_DisplayIdleProc, displayIdleData);
}


// Tcl proc launched via tcl-timer
void
tcl_InterruptIdleProc(ClientData data)
{
}


// *** WishPanic
void
WishPanic TCL_VARARGS_DEF(char *,arg1)
{
  va_list argList;
  char buf[1024];
  char *format;

  format = TCL_VARARGS_START(char *,arg1,argList);
  vsprintf(buf, format, argList);

#if defined(WIN32)
  MessageBeep(MB_ICONEXCLAMATION);
  MessageBox(NULL, buf, "Fatal Tcl Error when running Elmer Front",
             MB_ICONSTOP | MB_OK | MB_TASKMODAL | MB_SETFOREGROUND);
#else
  cerr << "Fatal Tcl error when running Elmer Front: " << endl
       << buf << endl;
#endif

  exit(1);
}
