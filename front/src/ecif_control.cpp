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
Module:     ecif_control.cpp
Language:   C++
Date:       13.11.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_body.h"
#include "ecif_bodyElement.h"
#include "ecif_control.h"
#include "ecif_input.h"
#include "ecif_inputAbaqus.h"
#include "ecif_inputEmf.h"
#include "ecif_inputEgf.h"
#include "ecif_inputElmer.h"
#include "ecif_inputFidap.h"
#include "ecif_inputIges.h"
#include "ecif_inputIdeas.h"
#include "ecif_inputThetis.h"
#include "ecif_mesh.h"
#include "ecif_model.h"
#include "ecif_modelMeshManager.h"
// #include "ecif_modelOutputManager.h"
#include "ecif_process.h"
#include "ecif_renderer_OGL.h"
#include "ecif_userinterface.h"
#include "ecif_timer.h"

void read_Model_File(char*, Model*);

// ***** Constants used within main-function *****
const char* STD_IN  = "stdin";
const char* STD_OUT = "stdout";

Control::Control(Hinst app_name, UserInterface* ui)
{
  appInstance = app_name;

  breakEgfInput = false;
  breakEgfOutput = false;
  breakEmfInput = false;
  breakEmfOutput = false;
  breakMeshInput = false;
  breakMeshOutput = false;

  theModel = NULL;
  theRenderer = NULL;
  theUI = ui;
  Input::theControlCenter = this;
  Model::theControlCenter = this;
  ModelMeshManager::theControlCenter = this;
  ModelOutputManager::theControlCenter = this;
  UserInterface::theControlCenter = this;
  Renderer::theControlCenter = this;

  processTable = new ProcessTable;
}


void
Control::copyParameters(char* emf_filename)
{
  char fname[1 + 1024];

  FILE* emf_file = NULL;

  //---Test if file exists.
  if ( (emf_file = fopen(emf_filename, "r")) == NULL ) {
    theUI->errMsg(0, "Error! Can't open the file: ", emf_filename);
    return;
  }

  //---Read parameters from the emf-file
  ifstream in_file(emf_filename);

  InputEmf* in = new InputEmf(ECIF_ND, in_file, emf_filename);

  in->copyEmfParameters();
}


void
Control::Exit()
{
  delete theRenderer; theRenderer = NULL;
  exit(1);
}


bool
Control::getBreakValue(enum frontProcessType process)
{
  switch (process) {
    case EGF_INPUT: return breakEgfInput;
    case EGF_OUTPUT: return breakEgfOutput;
    case EMF_INPUT: return breakEmfInput;
    case EMF_OUTPUT: return breakEmfOutput;
    case MESH_INPUT: return breakMeshInput;
    case MESH_OUTPUT: return breakMeshOutput;
  }

  return false;
}


void
Control::activateUI()
{
  if (theUI == NULL)
    return;

  enum ecif_modelDimension dim = theModel->getDimension();

  theUI->configureButtons("moveButtons", 1);
  theUI->configureButtons("displayButtons", 1);

  theUI->configureButtons("draw_target_bodies", 1);
  theUI->configureButtons("draw_target_edges", 1);

  if (dim == ECIF_2D) {
    theUI->configureButtons("draw_target_surfaces", 0);
    theUI->configureButtons("rotateButtons2D", 1);
  } else {
    theUI->configureButtons("draw_target_surfaces", 1);
    theUI->configureButtons("rotateButtons", 1);
  }

  if ( theModel->getFlagValue(GEOMETRY_TYPE_CAD) ) {
    theUI->configureButtons("draw_source_cad", 1);
  } else {
    theUI->configureButtons("draw_source_cad", 0);
  }

  if ( theModel->getFlagValue(GEOMETRY_TYPE_MESH) ) {
    theUI->configureMenuButtons("File", "Mesh", 1);
    theUI->configureButtons("draw_source_mesh", 1);
  }
  else {
    theUI->configureMenuButtons("File", "Mesh", 0);
    theUI->configureButtons("draw_source_mesh", 0);
  }

  theUI->configureButtons("saveModel", 1);
  theUI->configureButtons("mesh", 1);
  theUI->configureButtons("solve", 1);
  theUI->configureButtons("results", 1);
  theUI->configureButtons("break", 0);
  theUI->configureButtons("clearMessageArea", 1);
  theUI->configureButtons("info_", 1);
  theUI->configureButtons("process", 1);

  theUI->configureMenuButtons("File", "Save", 1);
  theUI->configureMenuButtons("File", "Browsers", 1);
  theUI->configureMenuButtons("Edit", "All", 1);
  theUI->configureMenuButtons("Edit", "Editors", 1);
  theUI->configureMenuButtons("Display", "All", 1);
  theUI->configureMenuButtons("Problem", "All", 1);
  theUI->configureMenuButtons("Model", "All", 1);
  theUI->configureMenuButtons("Mesh", "All", 1);
  theUI->configureMenuButtons("Solver", "All", 1);
  theUI->configureMenuButtons("Run", "All", 1);
  theUI->setInitialState();
}


void
Control::deactivateUI()
{
  if (theUI == NULL)
    return;

  theUI->configureMenuButtons("File", "Save", 0);
  theUI->configureMenuButtons("File" ,"Browsers", 0);
  theUI->configureMenuButtons("File" , "LoadMesh", 0);
  theUI->configureMenuButtons("File", "Mesh", 0);
  theUI->configureMenuButtons("Edit", "All", 0);
  theUI->configureMenuButtons("Edit" ,"Editors", 0);
  theUI->configureMenuButtons("Edit", "MatcGmtr", 0);
  theUI->configureMenuButtons("Display", "All", 0);
  //theUI->configureMenuButtons("Problem", "All", 0);
  theUI->configureMenuButtons("Problem", "Model", 0);
  //theUI->configureMenuButtons("Model", "All", 0);
  theUI->configureMenuButtons("Model", "Model", 0);
  theUI->configureMenuButtons("Mesh", "All", 0);
  //theUI->configureMenuButtons("Solver", "All", 0);
  theUI->configureMenuButtons("Solver", "Model", 0);
  theUI->configureMenuButtons("Run", "NearlyAll", 0);

  theUI->configureButtons("break", 1);
}


Input*
Control::create_mesh_input(enum ecif_modelDimension m_dim, ifstream& in_file, char* mesh_filename)
{
  Input* in = NULL;

  // *** CHECK INPUT FILE ******
  bool is_abaqus = checkFname(mesh_filename, ".inp");
  bool is_elmer = checkFname(mesh_filename, ".header");
  bool is_ideas = checkFname(mesh_filename, ".unv");
  bool is_thetis = checkFname(mesh_filename, ".tmf");

  //****** KNOWN FILE FORMAT ******
  // Read mesh file and select renderer-type based on the model type.
  if (is_abaqus || is_elmer || is_ideas || is_thetis) {
    //---Select input-type and create model object
    if (is_abaqus) {
      in = new InputAbaqus(m_dim, in_file, mesh_filename);
    } else if (is_elmer) {
      in = new InputElmer(m_dim, in_file, mesh_filename);
    } else if (is_ideas) {
      in = new InputIdeas(m_dim, in_file, mesh_filename);
    } else if (is_thetis) {
      in = new InputThetis(m_dim, in_file, mesh_filename);
    }
  }

  return in;
}


//Method creates a renderer based on model dimension
void
Control::createRenderer(enum ecif_modelDimension m_dim, bool display_gmtr)
{
  // A renderer already exists
  if (theRenderer != NULL) {
    theRenderer->setDimension(m_dim);
    theRenderer->resetData();

  // Create a new renderer
  } else {
    theRenderer = new Renderer_OGL(appInstance, m_dim);
  }

  if ( display_gmtr ) {
    theRenderer->refresh();
  }

  // Renderer has to know about cylindric symmetry axis
  if ( theModel != NULL ) {
    theRenderer->setSimulationDimension(theModel->getSimulationDimension());
  }

}


// Draw the whole model in the Renderer window
void
Control::displayModel()
{
  if (theModel == NULL) {
    return;
  }

  if (theRenderer == NULL) {
    createRenderer(theModel->getDimension());
  }

  theModel->setWindowTitles();
  theRenderer->displayRenderer();
}


void
Control::getCurrentTimestamp(char* buffer)
{
  theUI->getCurrentTimestamp(buffer);
}


Process*
Control::getProcess(int process_nbr)
{
  ProcessTable::iterator itr = processTable->find(process_nbr);

  if ( itr != processTable->end() )
    return (*itr).second;
  else
    return NULL;
}


void
Control::handleKeyAction(enum keyAction action)
{
  switch (action) {
  case KEY_CTRL_B: theUI->sendCommandToGui("SELECT_BODIES"); break;
  case KEY_CTRL_H: theUI->sendCommandToGui("SELECT_BOUNDARIES"); break;
  case KEY_CTRL_L: theUI->sendCommandToGui("SELECT_LABELS"); break;
  case KEY_CTRL_R: theUI->sendCommandToGui("RESET_RENDERER"); break;
  case KEY_CTRL_X: theUI->sendCommandToGui("SET_ROTATE_PRIORITY_X"); break;
  case KEY_CTRL_Y: theUI->sendCommandToGui("SET_ROTATE_PRIORITY_Y"); break;
  case KEY_CTRL_Z: theUI->sendCommandToGui("SET_ROTATE_PRIORITY_Z"); break;
  }
}


// Print %-progress info for "nbr" (eg. element_id)
// when total number is "total_nbr" (eg. nofElements)
// MVe 8/97
iostream&
Control::print_progress_info(iostream& strm, int nbr, int total_nbr, char* text)
{
  double done = (int(1000.0* nbr/total_nbr))/10.0;

  strm << "(";
  strm << setiosflags(ios::right | ios::fixed);
  strm << setw(4) << setprecision(1);
  strm << done << "%)";

  if (text != NULL)
    strm << text;

  return strm;
}


bool
Control::processExists(int process_nbr)
{
  strstream strm;
  Process* process = getProcess(process_nbr);

  if (process == NULL || !process->exists() ) {
    return false;
  }

  return true;
}


bool
Control::processResume(int process_nbr)
{
  strstream strm;
  Process* process = getProcess(process_nbr);

  if (process == NULL || !process->exists() ) {
    strm << "Process: " << process_nbr << " not reachable!" << ends;
    theUI->showMsg(strm.str());
    return false;
  }

  ProcessId pid = process->getProcessId();

  strm << "Process:  " << process_nbr << " (PID " << pid << ") being resumed" << ends;
  theUI->showMsg(strm.str());

  return process->resume();
}


bool
Control::processSetPriorityLevel(int process_nbr, priorityLevel priority)
{
  strstream strm;
  Process* process = getProcess(process_nbr);

  if (process == NULL) {
    strm << "Process: " << process_nbr << " not reachable!" << ends;
    theUI->showMsg(strm.str());
    return false;
  }

  ProcessId pid = process->getProcessId();

  strm << "Process:  " << process_nbr << " (PID " << pid << ") priority being set to " << priority << ends;
  theUI->showMsg(strm.str());

  process->setPriorityLevel(priority);

  return true;
}


bool
Control::processStart(Process* process)
{
  strstream strm;

  int process_nbr = process->ID();

  (*processTable)[process_nbr] = process;


  if ( NULL == getProcess(process_nbr) ) {
    strm << "Process: " << process_nbr << " NOT started!" << ends;
    theUI->showMsg(strm.str());
    return false;
  }

  bool success = process->start();

  ProcessId pid = process->getProcessId();

  if (success) {
    //strm << "Process: " << process_nbr << " (PID " << pid << ") started!" << ends;
    //theUI->showMsg(strm.str());
    return true;
  } else {
    processTable->erase(process_nbr);
    strm << "Process: " << process_nbr << " (PID " << pid << ") NOT started and is removed!" << ends;
    theUI->showMsg(strm.str());
    return false;
  }

  return false;
}


bool
Control::processStop(int process_nbr)
{
  strstream strm;
  Process* process = getProcess(process_nbr);

  if (process == NULL) {
    strm << "Process: " << process_nbr << " not reachable!" << ends;
    theUI->showMsg(strm.str());
    return false;
  }

  ProcessId pid = process->getProcessId();
  const char* nm = process->getName();

  if ( nm == NULL || nm[0] == '\0' ) {
    strm << "Process:  " << process_nbr << " (PID " << pid << ") being stopped!" << ends;
  } else {
    strm << "Process:  " << process_nbr << " (" << nm << ") being stopped!" << ends;
  }

  theUI->showMsg(strm.str());

  bool stopped = process->stop();

  if (stopped) {
    processTable->erase(process_nbr);
    return true;
  }

  return false;
}


bool
Control::processSuspend(int process_nbr)
{
  strstream strm;
  Process* process = getProcess(process_nbr);

  if (process == NULL) {
    strm << "Process: " << process_nbr << " not reachable!" << ends;
    theUI->showMsg(strm.str());
    return false;
  }

  ProcessId pid = process->getProcessId();
  strm << "Process:  " << process_nbr << " (PID " << pid << ") being suspended" << ends;
  theUI->showMsg(strm.str());

  return process->suspend();
}


// Reads in CAD model data file.
//
bool
Control::readCADFile(char* CAD_filename,
                     char* CAD_type,
                     enum ecif_modelDimension m_dim)
{
  char fname[1 + 1024];

  FILE* CAD_file = NULL;

  // *** Check that input-file exists.
  if ( (CAD_file = fopen(CAD_filename, "r")) == NULL ) {
    theUI->errMsg(0, "Error! Can't open the file: ", CAD_filename);
    return false;
  }

  fclose(CAD_file);

  // *** Ok, continue

  // Deactivate first menus
  deactivateUI();

  strncpy(fname, CAD_filename, 1024);
  LibFront::toLower(fname);

  // *** CHECK INPUT FILE ******
  bool is_egf   = false;
  bool is_ideas = false;
  bool is_iges  = false;

  // Quess fron the file extension
  if ( CAD_type == NULL ||
       CAD_type[0] == '\0' ||
       LibFront::in(CAD_type, "default")
     ) {
    is_egf   = checkFname(fname, ".egf");
    is_ideas = checkFname(fname, ".unv");
    is_iges  = checkFname(fname, ".igs");

  // User has specified the type
  } else {
    if ( LibFront::in(CAD_type, "elmer") ) is_egf = true;
    else if ( LibFront::in(CAD_type, "iges") ) is_iges = true;
    else if ( LibFront::in(CAD_type, "ideas") ) is_ideas = true;
  }

  //****** UNKNOWN FILE FORMAT ******
  if ( !(is_egf || is_ideas || is_iges) ) {
    theUI->errMsg(0, "Error! Unknown CAD-file format: ", CAD_filename);
    return false;
  }

  //****** KNOWN FILE FORMAT ******
  // Read CAD file and select renderer-type based on the model type.

  // Clear renderer window
  if ( theRenderer != NULL ) {
    theRenderer->clear();
  }

  // First delete possible old model
  if (theModel != NULL) {
    delete theModel; theModel = NULL;
  }

  Input* in = NULL;

  //ifstream in_file(CAD_filename, ios::nocreate);
  ifstream in_file(CAD_filename);

  theModel = new Model("", ECIF_CAD_FILE, CAD_filename);

  //---Select input-type and create model object
  if (is_egf) {
    // NOTE: Egf-file errors are handled by the input routine, so
    // we can give the dimension already here --> Dimension
    // keyword is not necessary in the egf-file
    // MVe, 12.01.03
    m_dim = ECIF_2D;
    in = new InputEgf(m_dim, in_file, CAD_filename);
  }
  else if (is_ideas) {
    in = new InputIdeas(m_dim, in_file, CAD_filename);
  }
  else if (is_iges) {
    in = new InputIges(m_dim, in_file, CAD_filename);
  }

  //---Read CAD-file, update model-type  and create suitable renderer
  enum ecif_modelDimension m_out_dim;

#if defined(FRONT_DEBUG)

  m_out_dim = in->readCadFile();

#else
  // Try to read cad-geometry file
  // -----------------------------
  try {
      m_out_dim = in->readCadFile();
  }

  // If error
  // --------
  catch (...) {
    m_out_dim = ECIF_ND;
  }
#endif

  in_file.close();

  //---Check that dat was ok.
  if ( m_out_dim == ECIF_ND ) {
    theUI->errMsg(0, "Error! Incorrect data in CAD-inputfile: ", CAD_filename);
    delete in;
    return false;
  }

  //---Ok, process Cad data and create suitable renderer
  bool rc = theModel->processCadFileData();
  if (rc == false) {
    theUI->errMsg(0, "Error! Cad file data not correct!");
    return false;
  }

  //---Ok. Create renderer
  createRenderer(m_out_dim);

  //---Update data and activate GUI menus.
  updateUI();

  delete in;
  return true;

}


// Reads in mesh data file.
//
bool
Control::readMeshFile(char* mesh_filename,
                      char* mesh_type,
                      bool create_new_model,
                      enum ecif_modelDimension m_dim)
{
  char fname[1 + 1024];

  FILE* mesh_file = NULL;

  // *** Check that input-file exists.
  if ( (mesh_file = fopen(mesh_filename, "r")) == NULL ) {
    theUI->errMsg(0, "Error! Can't open the file: ", mesh_filename);
    return false;
  }

  fclose(mesh_file);

  deactivateUI();

  strncpy(fname, mesh_filename, 1024);
  LibFront::toLower(fname);

  //*** CHECK INPUT FILE TYPE ***
  bool is_abaqus = false;
  bool is_elmer = false;
  bool is_fidap = false;
  bool is_ideas = false;
  bool is_thetis = false;

  if ( LibFront::in(mesh_type, "abaqus2d") ) {
    is_abaqus = true;
    m_dim = ECIF_2D;
  } else if ( LibFront::in(mesh_type, "abaqus3d") ) {
    is_abaqus = true;
    m_dim = ECIF_3D;

  } else if ( LibFront::in(mesh_type, "elmer2d") ) {
    is_elmer = true;
    m_dim = ECIF_2D;
  } else if ( LibFront::in(mesh_type, "elmer3d") ) {
    is_elmer = true;
    m_dim = ECIF_3D;

  } else if ( LibFront::in(mesh_type, "fidap2d") ) {
    is_fidap = true;
    m_dim = ECIF_2D;
  } else if ( LibFront::in(mesh_type, "fidap3d") ) {
    is_fidap = true;
    m_dim = ECIF_3D;

  } else if ( LibFront::in(mesh_type, "ideas2d") ) {
    is_ideas = true;
    m_dim = ECIF_2D;
  } else if ( LibFront::in(mesh_type, "ideas3d") ) {
    is_ideas = true;
    m_dim = ECIF_3D;

  } else if ( LibFront::in(mesh_type, "thetis2d") ) {
    is_thetis = true;
    m_dim = ECIF_2D;
  } else if ( LibFront::in(mesh_type, "thetis3d") ) {
    is_thetis = true;
    m_dim = ECIF_3D;
  }

  //****** UNKNOWN FILE FORMAT ******
  if ( !(is_abaqus || is_elmer || is_fidap || is_ideas || is_thetis) ) {
    theUI->errMsg(0, "Error! Unknown mesh-file format: ", mesh_filename);
    return false;
  }

  //*** KNOWN FILE FORMAT ***
  // Read mesh file and select renderer-type based on the model type.

  Input* in = NULL;

  //ifstream in_file(mesh_filename, ios::nocreate);
  ifstream in_file(mesh_filename);

  //---Select input-type and create model object
  if (is_abaqus) {
    in = new InputAbaqus(m_dim, in_file, mesh_filename);

  } else if (is_elmer) {
    in = new InputElmer(m_dim, in_file, mesh_filename);

  } else if (is_fidap) {
    in = new InputFidap(m_dim, in_file, mesh_filename);

  } else if (is_ideas) {
    in = new InputIdeas(m_dim, in_file, mesh_filename);

  } else if (is_thetis) {
    in = new InputThetis(m_dim, in_file, mesh_filename);
  }

  // Clear renderer window
  if ( theRenderer != NULL ) {
    theRenderer->clear();
  }

  // Delete current model and cretate a new model
  if ( create_new_model || theModel == NULL ) {

    // First delete possible old model and create new model
    if (theModel != NULL) {
      delete theModel; theModel = NULL;
    }

    // Create new model
    theModel = new Model("", ECIF_MESH_FILE, mesh_filename);

  // Delete only current mesh data, but use bodies etc.
  // NOTE: Used when remeshing an external mesh!!!
  } else if ( theModel != NULL ) {
    theModel->resetModelData();
    theModel->resetMeshData();
  }

  // Elmer mesh specific
  if (is_elmer) {
    theModel->updateModelNameInfo();
    theModel->updateModelDirectoryInfo();
    theModel->updateMeshDirectoryInfo();
  }

  //---Read meshfile, update model-type  and create suitable renderer
  ecif_modelDimension m_out_dim = ECIF_ND;
  bool mesh_ok = true;

#if defined(FRONT_DEBUG)

  m_out_dim = in->readMeshFile();

#else
  // Try to read mesh-geomtry file
  // -----------------------------
  //
  try {
      m_out_dim = in->readMeshFile();
  }

  // If error
  //
  catch (...) {
    m_out_dim = ECIF_ND;
  }
#endif

  in_file.close();

  if ( m_out_dim == ECIF_ND ) {
    mesh_ok = false;

    if ( !breakMeshInput ) {
      theUI->errMsg(0, "Error! Incorrect data in the mesh inputfile: ", mesh_filename);
    }
  }

  //---If ok so far, process mesh data and create suitable renderer
  if (mesh_ok) {

    mesh_ok = theModel->processMeshFileData(in);

    if (!mesh_ok) {
      theUI->errMsg(0, "Error! Mesh file data not correct, mesh not read!");
    }
  }

  delete in;
  breakMeshInput = false;

  //---If not ok
  if (!mesh_ok) {
    theModel->resetMeshData();
    theUI->configureButtons("break", 0);
    return false;
  }

  //---Ok. Create renderer
  createRenderer(m_out_dim);

  //---Update data and activate GUI menus
  updateUI();

  // NOTE: We reset always the mesh unit factor after a succesful read!!!
  //
  Model::setMeshInputUnit(1.0);

  return true;
}


// Read model (emf) file
//
bool
Control::readModelFile(char* model_filename, bool load_mesh, bool is_batch)
{
  char fname[1 + 1024];

  FILE* model_file = NULL;

  //---Test if file exists.
  if ( model_filename == NULL || model_filename[0] == '\0' ) {
    theUI->errMsg(0, "***ERROR No model file given!");
    return false;
  }

  //---Test if file exists.
  if ( (model_file = fopen(model_filename, "r")) == NULL ) {
    theUI->errMsg(0, "***ERROR Can't open the file: ", model_filename);
    return false;
  }

  fclose(model_file);

  //---Test if we can write to the file
  if ( !fopen(model_filename, "r+") ) {
    theUI->errMsg(0, "***ERROR Could not get write access to the model file:", model_filename, " , file opened in read only mode!");
    fclose(model_file);
  }

  strncpy(fname, model_filename, 1024);
  LibFront::toLower(fname);

  // Check file type
  //bool is_emf  = checkFname(fname, ".emf");
  bool is_emf  = true;

  if (is_emf) {

    // First deactivate menus
    deactivateUI();

    // Clear renderer window
    if ( theRenderer != NULL ) {
      theRenderer->clear();
    }

    // Delete possible old model
    if (theModel != NULL) {
      delete theModel; theModel = NULL;
    }

    //---Create model object
    theModel = new Model();

    //---Update model directory info
    theModel->updateModelDirectoryInfo();

    //---Read model file
    //ifstream in_file(model_filename, ios::nocreate);
    ifstream in_file(model_filename);

    InputEmf* in = new InputEmf(ECIF_ND, in_file, model_filename);

    enum ecif_modelDimension m_out_dim;

    theModel->setReadingModelFile(true);

    // To be sure (for model mesh!)
    Model::setMeshInputUnit(1.0);

#if defined(FRONT_DEBUG)

    m_out_dim = in->readCadFile();

#else
    // Try to read model emf-file
    // --------------------------
    try {
        m_out_dim = in->readCadFile();
    }

    // If error
    // --------
    catch (...) {
      m_out_dim = ECIF_ND;
    }
#endif

    theModel->setReadingModelFile(false);

    in_file.close();
    fclose(model_file);

    //---Check that data was ok.
    if ( m_out_dim == ECIF_ND ) {
      theUI->errMsg(0, "***ERROR Incorrect data in Model-inputfile: ", model_filename);
      return false;
    }

    //---Ok, process model file data and create suitable renderer
    bool rc = theModel->processModelFileData();

    if (rc == false) {
      theUI->errMsg(0, "***ERROR Model file data not correct!");
      return false;
    }

    if (!is_batch) updatePreUI();

    //---Update mesh directory info
    theModel->updateMeshDirectoryInfo();

    //---Mesh loading turned off, make load-button active
    if ( !load_mesh && !is_batch ) {
      theUI->configureMenuButtons("File", "LoadMesh", 1);
    }
    //---Load Elmer mesh
    else {
      theModel->loadMesh();
    }

    theModel->markActiveObjects();

    if (!is_batch) updatePostUI();

    //---Create renderer
    bool display_gmtr = theModel->modelHasCadGeometry() ||
                        (theModel->modelHasMeshGeometry() && load_mesh);
    if (!is_batch) createRenderer(m_out_dim, display_gmtr);

    return true;
  }

  //****** UNKNOWN FILE FORMAT ******
  else {
    theUI->errMsg(0, "***ERROR Unknown Model file format: ", model_filename);
    return false;
  }
}


void
Control::rendererIsClosed()
{
  delete theRenderer;
  theRenderer = NULL;
}


void
Control::setModelDimension(enum ecif_modelDimension model_dim)
{
  theModel->setModelDimension(model_dim);
}


// Saves mesh file in Elmer format
void
Control::saveElmerMeshFile(char* mesh_dir)
{
  theModel->saveElmerMesh(mesh_dir);
}


// Saves mesh result file in ElmerPost format
void
Control::saveElmerPostMeshFile(char* out_filename)
{
  theModel->saveElmerPostMesh(out_filename);
}


// Saves mesh result file in Thetis format
void
Control::saveThetisMeshFile(char* out_filename)
{
  theModel->saveThetisMesh(out_filename);
}


// Saves user settings into a (default) file
void
Control::saveUserSettingsFile(char* out_filename)
{
  theModel->saveUserSettingsFile(out_filename);
}


// Saves model file (.emf)
void
Control::saveFrontModelFile(char* out_filename)
{
  ofstream out_file(out_filename, ios::out);

  if ( !write_ok(out_file, out_filename, "(Writing ELMER Front model file)") )
    return;

#if defined(FRONT_DEBUG)

  theModel->saveFrontModelFile(out_file, out_filename);

#else
  try
  {
    theModel->saveFrontModelFile(out_file, out_filename);
  }

  catch (...)
  {
    theUI->showMsg("ERROR: Unable to save model file (emf-file)!");
  }
#endif

  out_file.close();
}


// Saves Mesh input (.mif) file.
void
Control::saveMeshInputFile(char* out_filename)
{
  ofstream out_file(out_filename, ios::out);

  if ( !write_ok(out_file, out_filename, "(Writing ELMER Mesh input file)") )
    return;

#if defined(FRONT_DEBUG)

  theModel->saveMeshInputFile(out_file, out_filename);

#else
  try
  {
    theModel->saveMeshInputFile(out_file, out_filename);
  }

  catch (...)
  {
    theUI->showMsg("ERROR: Unable to save Mesh2D input file (mif-file)!");
  }
#endif

  out_file.close();
}


// Saves Solver input (.sif) file.
void
Control::saveSolverInputFile(char* out_filename)
{
  ofstream out_file(out_filename, ios::out);

  if ( !write_ok(out_file, out_filename, "(Writing ELMER Solver input file)") )
    return;

#if defined(FRONT_DEBUG)

  theModel->saveSolverInputFile(out_file, out_filename);

#else
  try
  {
    theModel->saveSolverInputFile(out_file, out_filename);
  }

  catch (...)
  {
    theUI->showMsg("ERROR: Unable to save Solver input file (sif-file)!");
  }
#endif

  out_file.close();
}


// Sets drawing mode for a body
void
Control::selectBody(int bd_id, int lr_id, bool update_gui)
{
  return;

  Body* body = theModel->getBodyById(bd_id);

  objectDrawingMode dmode = body->getDrawMode();

  if ( dmode == DM_NORMAL ) {
    body->setDrawMode(DM_HIDDEN);

  } else {
    body->setDrawMode(DM_NORMAL);
  }

  if (theRenderer != NULL) {
    //theRenderer->setParameters();
    theRenderer->drawAllBodies();
  }

}


void
Control::selectBoundary(int elem_id, int body1_id, int layer1_id,
                        int body2_id, int layer2_id,
                        bool accept_body_change, bool update_gui)
{
  if (theRenderer == NULL) {
    createRenderer(theModel->getDimension());
  }

  theRenderer->selectBoundary(elem_id,
                              body1_id, layer1_id,
                              body2_id, layer2_id,
                              accept_body_change, update_gui);
}


void
Control::setBreakValue(enum frontProcessType process, bool value)
{
  switch (process) {
    case EGF_INPUT: breakEgfInput = value;
    case EGF_OUTPUT: breakEgfOutput = value;
    case EMF_INPUT: breakEmfInput = value;
    case EMF_OUTPUT: breakEmfOutput = value;
    case MESH_INPUT: breakMeshInput = value;
    case MESH_OUTPUT: breakMeshOutput = value;
  }
}


void
Control::selectBoundaries(int nof_elems, int* elem_ids,
                          int* body1_ids, int* layer1_ids,
                          int* body2_ids, int* layer2_ids,
                          bool accept_body_change, bool update_gui)
{
  if (theRenderer == NULL) {
    createRenderer(theModel->getDimension());
  }

  theRenderer->selectBoundaries(nof_elems, elem_ids,
                                body1_ids, layer1_ids,
                                body2_ids, layer2_ids,
                                accept_body_change, update_gui);
}

void
Control::setGuiWindowTitle(char* case_name)
{
  if ( theUI == NULL ) {
    return;
  }

  theUI->setWindowTitle(case_name);

}


void
Control::setWindowTitles(char* case_name)
{
  setRendererWindowTitle(case_name);
  setGuiWindowTitle(case_name);
}


void
Control::setRendererWindowTitle(char* case_name)
{
  if ( theModel == NULL || theRenderer == NULL ) {
    return;
  }

  double x1, x2, y1, y2, z1, z2;
  strstream strm;


  theModel->getBoundingBox(x1, x2, y1, y2, z1, z2);

  strm << case_name;

  strm << "  (" << x2 - x1 << " x " << y2 - y1;
  if (theModel->getDimension() == ECIF_3D) {
    strm << " x " << z2 - z1;
  }
  strm << ")";

  const char* mesh_name = theModel->getCurrentMeshName();
  bool drawing_mesh = theModel->getFlagValue(DRAW_SOURCE_MESH);

  if (drawing_mesh && mesh_name != NULL && mesh_name[0] != '\0' ) {
    strm << "    " << mesh_name;
  }

  strm << ends;

  theRenderer->setWindowTitle(strm.str());
}


int
Control::unknownFieldMsg(emf_ObjectData_X* object_data, bool is_fatal)
{
  const int isOK = 0;

  UserInterface* gui = getGui();
  strstream strm;
  strm << "WARNING: Unknown field name ("
       << object_data->field_name
       << ") when reading object: "
       << object_data->object_name
       << ends;

  gui->showMsg(strm.str());

  if (is_fatal)
    return !isOK;
  else
    return isOK;
}


int
Control::unknownObjectMsg(emf_ObjectData_X* object_data, bool is_fatal)
{
  const int isOK = 0;

  UserInterface* gui = getGui();
  strstream strm;
  strm << "WARNING: Unknown object name: "
       << object_data->object_name
       << ends;

  gui->showMsg(strm.str());

  if (is_fatal)
    return !isOK;
  else
    return isOK;
}


void
Control::update(int counter, int update_interval)
{
  theUI->update(counter, update_interval);
}

void
Control::updateUI()
{
  updatePreUI();
  updatePostUI();
}


void
Control::updatePreUI()
{
  //theUI->showMsg("---1. Updating model objects...");

  // Update general model level data in Gui
  theModel->initModelFlags();
  theUI->updateModelFlags(theModel);
  theUI->updateModelStatus(theModel);

  theUI->updateModelData(theModel);

  theUI->updateObjectData(theModel);

  theUI->updateParameterDataPre(theModel);
}

void
Control::updatePostUI()
{
  // Why this was called again here (also in updateObjectData) ???!!!???
  //
  // Update boundary level data
  //theUI->updateBoundaryData(theModel);

  theUI->showMsg("Checking model data...");

  // Set mask related stuff (when all boundaries are now set!)
  theUI->updateParameterDataPost(theModel);

  // Activate menus
  activateUI();

  theModel->setWindowTitles();

  theUI->showMsg("All done");
}


bool
Control::write_ok(ofstream& out_file, char* out_filename, char* msg)
{
  int state = out_file.rdstate();

  if ( state & (ios::failbit | ios::badbit) ) {

    strstream strm;
    strm << "Error! Cannot write data to the file:" << endl;
    strm << out_filename;
    if ( msg != NULL ) {
      strm << endl << msg;
    }
    strm << ends;

    theUI->errMsg(0, strm.str());
    return false;
  }

  return true;
}





