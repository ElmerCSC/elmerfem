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
Module:     ecif_inputEmf.cpp
Language:   C++
Date:       20.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Implementation (read Elmer model files (.emf files)

************************************************************************/
 
#include "ecif_body2D.h"
#include "ecif_body3D.h"
#include "ecif_bodyElement.h"
#include "ecif_control.h"
#include "ecif_geometry.h"
#include "ecif_inputEmf.h" 
#include "ecif_model.h"
#include "ecif_model_aux.h"
#include "ecif_parameter.h" 
#include "ecif_userinterface.h"
#include "ecif_timer.h"

extern char read_buffer[];

extern bool IS_FATAL;
extern bool NOT_FATAL;

const int isOk = emf_OK;
const int notOk = emf_ERROR;

static int int_buffer[10000]; 

template <class T> ostream&
output_array(ostream& out, UserInterface* gui, int length, T* array)
{
  ostrstream strm;
  for (int i = 0; i < length; i++) {
    strm << array[i]; 
  }
  strm << ends;
  gui->showMsg(strm.str());

  return out;
}


struct ParamCopyData {
  bool add_to_model;
  ParameterList paramList;
};

    
InputEmf::InputEmf(enum ecif_modelDimension m_dim,
                       ifstream& in_file, char* in_filename):
InputFront(m_dim, in_file, in_filename)
{
  isChecked = true;
  isEmfInput = true;
}


void
InputEmf::copyEmfParameters()
{ 
  UserInterface* gui = theControlCenter->getGui();
  char objSep = '!';
  char fldSep = '&';

  int i;

  struct emf_ObjectData_X my_ObjectData;
  emf_ObjectData  = &my_ObjectData;

  // Set callback functions
  emf_readDataCB = copyEmfParametersCallBack;
  emf_readMessageCB = readEmfGeometryMsgCallBack;

  /*-----Start parser*/

  ParamCopyData cpd;
  char msg_buffer[256];

  cpd.add_to_model = false;

  int rc = emf_readData(infileName, (void**)&cpd, 255, msg_buffer, 0);
  
  // Parameters which were read from input
  //
  int nof_params = cpd.paramList.size();

  ParameterList::iterator itr = cpd.paramList.begin();

  // Send list of parameters data to Gui to be finally selected
  // by the user
  strstream pdatas;

  for (i = 0; i < nof_params; i++, itr++) {
    Parameter* param = *itr;
    //pdatas << param->getName() << fldSep;
    // Ok, we do not add any extra info
    pdatas << param->getEmfName() << ": ";
    pdatas << param->getName();
    
#if 0    
    // Possible correspondig parent object from the current model
    // NOTE: Matching objects within different models would be so
    // messy that we do not even try it, so no further data is added!!
    //Emf-file object
    objectType etp = param->getParentEmfType();
    const char* enm = model->objectType2Name(etp);
    int etg = param->getParentEmfTag();

    pdatas << enm << fldSep;
    pdatas << etg << fldSep;

    // Current model object
    ModelObject* obj = model->getModelObjectByTag(etp, etg);
    if ( obj != NULL ) {
      objectType otp = obj->getObjectType();
      const char* onm = model->objectType2Name(otp);
      int otg = obj->Tag();
      pdatas << onm << fldSep;
      pdatas << otg << fldSep;

    } else {
      pdatas << "";
      pdatas << NO_INDEX;
    }
#endif

    if ( i < nof_params - 1 ) {
      pdatas << objSep;
    }
  }

  pdatas << ends;

  bool* accept_flags = new bool[nof_params];

  // Ask the user to select parameters to be copied
  //
  gui->acceptEmfParameters("Select parameters to copy", nof_params, pdatas.str(), accept_flags);

  // Copy selected parametrs
  //
  itr = cpd.paramList.begin();

  int count = 0;

  for (i = 0; i < nof_params; i++, itr++) {
    Parameter* param = *itr;

    if ( !accept_flags[i] ) {
      delete param;
      continue;
    }

    count = count + 1;
    
    // Set correct ids for the copied parameter
    //
    int lid = param->getLastId();
    param->setId(1+lid);
    param->setLastId(1+lid);
    param->setParentId(NO_INDEX);
    param->setParentEmfTag(NO_INDEX);
    param->setParentEmfType(OT_NONE);
    param->setSubParentEmfTag(NO_INDEX);
    param->setSubParentEmfType(OT_NONE);
    
    // Select a non-existing name for the paramter to be copied
    // Contrain1 Constrain1-c Contraint1-c-c ...
    //
    while ( model->modelHasParameter(param->getParameterType(), param->getName()) ) {
      strstream strm;
      strm << param->getName() << "-c" << ends;
      param->setName(strm.str());
    }

    // Add parameter
    model->addParameter(param->getParameterType(), param);
  }

  delete[] accept_flags;
  
  // Info message and update the model data in gui
  if ( count > 0 ) {
    strstream strm;
    strm << "Parameters copied: " << count << ends;
    gui->showMsg(strm.str());
    gui->updateParameterDataPre(model);
    gui->updateParameterDataPost(model);
  }

  return;
}


// This method gets one "field" at a time from the file
//
// NOTE: Only those parameter which can have multiple intances
// can be copied
// NOTE: Equations cannot be copied, because they may need 
// an EquationVariables-parameter which a single instance creature
//
int
InputEmf::copyEmfParametersCallBack(void** user_data)
{
  UserInterface* gui = theControlCenter->getGui();

  ParamCopyData* cpd = (ParamCopyData*)user_data;

  // Global object pointing to the record, make
  // a short hand for it
  emf_ObjectData_X* od = emf_ObjectData;

  const char* on = od->object_name;
  int onl = od->object_name_length;
  const char* fn = od->field_name;
  int fnl = od->field_name_length;
  const char* dt = od->data_type;
  int dtl = od->data_type_length;

  const Parameter* param = NULL;
  bool add_to_model = false;

  // Read data according to the object type

  // Parameters
  // ==========
  int rc;

  //-Body Force
  if ( LibFront::in(EMF_BODY_FORCE, on) ) {
    rc = readEmfParameter(od, ECIF_BODY_FORCE, param, add_to_model);
  }

  //-Body Parameter
  else if ( LibFront::in(EMF_BODY_PARAMETER, on) ) {
    rc = readEmfParameter(od, ECIF_BODY_PARAMETER, param, add_to_model);
  }

  //-Boundary condition
  else if ( LibFront::in(EMF_BOUNDARY_CONDITION, on) ) {
    rc = readEmfParameter(od, ECIF_BOUNDARY_CONDITION, param, add_to_model);
  }

  //-Boundary Parameter
  else if ( LibFront::in(EMF_BOUNDARY_PARAMETER, on) ) {
    rc = readEmfParameter(od, ECIF_BOUNDARY_PARAMETER, param, add_to_model);
  }

  //-Calculator solver
  else if ( LibFront::in(EMF_CALCULATOR, on) ) {
    rc = readEmfParameter(od, ECIF_CALCULATOR, param, add_to_model);
  }

  //-Physical Constant
  else if ( LibFront::in(EMF_CONSTANT, on) || LibFront::in("Constants", on) ) {
    //rc = readEmfParameter(od, ECIF_CONSTANT, param, add_to_model);
  }

  //-Coordinate
  else if ( LibFront::in(EMF_COORDINATE, on) || LibFront::in("Coordinates", on) ) {
    //rc = readEmfParameter(od, ECIF_COORDINATE, param, add_to_model);
  }

  //-Equation
  else if ( LibFront::in(EMF_EQUATION, on) ) {
    //rc = readEmfParameter(od, ECIF_EQUATION, param, add_to_model);
  }

  //-EquationVariables
  else if ( LibFront::in(EMF_EQUATION_VARIABLE, on) || LibFront::in("Equation Variables", on) ) {
    //rc = readEmfParameter(od, ECIF_EQUATION_VARIABLE, param, add_to_model);
  }

  //-Grid H parameter
  else if ( LibFront::in(EMF_GRID_H, on) || LibFront::in("GridH", on) ) {
    rc = readEmfParameter(od, ECIF_GRID_H, param, add_to_model);
  }

  //-Grid parameter
  else if ( LibFront::in(EMF_GRID_PARAMETER, on) ||
            LibFront::in("Mesh Parameter", on)
    ) {
    rc = readEmfParameter(od, ECIF_GRID_PARAMETER, param, add_to_model);
  }

  //-Initial condition
  else if ( LibFront::in(EMF_INITIAL_CONDITION, on) ) {
    rc = readEmfParameter(od, ECIF_INITIAL_CONDITION, param, add_to_model);
  }

  //-Material
  else if ( LibFront::in(EMF_MATERIAL, on) ) {
    rc = readEmfParameter(od, ECIF_MATERIAL, param, add_to_model);
  }

  //-Model parameter
  else if ( LibFront::in(EMF_MODEL_PARAMETER, on) ) {
    rc = readEmfParameter(od, ECIF_MODEL_PARAMETER, param, add_to_model);
  }

  //-Simulation parameter
  else if ( LibFront::in(EMF_SIMULATION_PARAMETER, on) ) {
    rc = readEmfParameter(od, ECIF_SIMULATION_PARAMETER, param, add_to_model);
  }

  //-Solver
  else if ( LibFront::in(EMF_SOLVER, on) ) {
    //rc = readEmfParameter(od, ECIF_SOLVER, param, add_to_model);
  }

  //-SolverControl parameter
  else if ( LibFront::in(EMF_SOLVER_CONTROL, on) ) {
    //rc = readEmfParameter(od, ECIF_SOLVER_CONTROL, param, add_to_model);
  }

  //-Timestep
  else if ( LibFront::in(EMF_TIMESTEP, on) ) {
    rc = readEmfParameter(od, ECIF_TIMESTEP, param, add_to_model);
  }

  //-User settings
  else if ( LibFront::in(EMF_USER_SETTING, on) || LibFront::in("User Settings", on) ) {
    //rc = readEmfParameter(od, ECIF_USER_SETTING, param, add_to_model);
  }

  //-UNKNOWN object
  else {
    //return unknownObjectMsg(od, false);
  }

  if ( param != NULL ) {

    const char* pd = ((Parameter*)param)->getValue();

    if ( pd != NULL && pd[0] != '\0' && pd != "" ) {
      cpd->paramList.push_back((Parameter*)param);
    }
  }

  return isOk;
}

 
// Main function for the file reading
bool
InputEmf::readCadGeometry()
{
  struct emf_ObjectData_X my_ObjectData;
  emf_ObjectData  = &my_ObjectData;

  // Set callback functions
  //int my_readDataCB(void**);
  emf_readDataCB = readEmfGeometryCallBack;
  emf_readMessageCB = readEmfGeometryMsgCallBack;

  InputEmf* ci = this;
  char msg_buffer[256];

  // Start parser by calling the interface function
  int rc = emf_readData(infileName, (void**)&ci, 255, msg_buffer, 0);


  if (rc != 0) {
    modelDimension = ECIF_ND;
    return false;
  }

  // Convert tags to ids
  model->convertTags2Ids();

  return true;

}
 

int
InputEmf::readEmfGeometryMsgCallBack(char* msg_buffer)
{
  UserInterface* gui = theControlCenter->getGui();
  return gui->showMsg(msg_buffer, 0, true);
}


// This method gets one "field" at a time from the file
int
InputEmf::readEmfGeometryCallBack(void** user_data)
{
  UserInterface* gui = theControlCenter->getGui();

  InputEmf* emfIn = (InputEmf*)*user_data;
  InputFront* frontIn = (InputFront*)*user_data;

  // Global object pointing to the record, make
  // a short hand for it
  emf_ObjectData_X* od = emf_ObjectData;

  const Parameter* param;

  const char* on = od->object_name;
  int onl = od->object_name_length;
  const char* fn = od->field_name;
  int fnl = od->field_name_length;
  const char* dt = od->data_type;
  int dtl = od->data_type_length;

  // Read data according to the object type

  // Geometry related
  // ================

  //-Body
  if ( LibFront::in(EMF_BODY, on) ) {
    return readBody(od);
  }

  //-Element group
  else if ( LibFront::in(EMF_ELEMENT_GROUP, on) ) {
    return readElementGroup(od);
  }

  //-Body element (inputVersionNbr <= 4 style!)
  else if ( LibFront::in("Element", on) ) {

    if ( model->getDimension() ==ECIF_2D )
      return readEdge(od);
    else
      return readFace(od);
  }

  //-Vertex (inputVersionNbr >= 5 style!)
  else if ( LibFront::in(EMF_VERTEX, on) ) {
    return readVertex(od);
  }

  //-Edge (inputVersionNbr >= 5 style!)
  else if ( LibFront::in(EMF_EDGE, on) ) {
    return readEdge(od);
  }

  //-Face (inputVersionNbr >= 5 style!)
  else if ( LibFront::in(EMF_FACE, on) ) {
    return readFace(od);
  }

  //-Element loop
  else if ( LibFront::in(EMF_ELEMENT_LOOP, on) ) {
    return readElementLoop(od);
  }

  //-Inner boundaries (pre version 5 stuff!)
  else if ( LibFront::in("Inner Boundaries", on) ) {
    return isOk;
  }

  //-Outer boundaries (pre version 5 stuff!)
  else if ( LibFront::in("Outer Boundaries", on) ) {
    return isOk;
  }

  //-Boundary vertices table (pre version 5 stuff!)
  else if ( LibFront::in("Boundary Vertices", on) ) {
    //return readFrontBoundaryVertices(emfIn, od);
  }

  //-Vertex table (version 9 -->)
  else if ( LibFront::in(EMF_VERTICES, on) || 
            LibFront::in(EMF_VERTEX_TABLE, on)  
          ) {
    return readVertexTable(od);
  }


  // Parameters
  // ==========

  //-Body Force
  else if ( LibFront::in(EMF_BODY_FORCE, on) ) {
    return readEmfParameter((emf_ObjectData_X*)od, ECIF_BODY_FORCE, param, true);
  }

  //-Body Parameter
  else if ( LibFront::in(EMF_BODY_PARAMETER, on) ) {
    return readEmfParameter(od, ECIF_BODY_PARAMETER, param, true);
  }

  //-Boundary condition
  else if ( LibFront::in(EMF_BOUNDARY_CONDITION, on) ) {
    return readEmfParameter(od, ECIF_BOUNDARY_CONDITION, param, true);
  }

  //-Boundary Parameter
  else if ( LibFront::in(EMF_BOUNDARY_PARAMETER, on) ) {
    return readEmfParameter(od, ECIF_BOUNDARY_PARAMETER, param, true);
  }

  //-Calculator solver
  else if ( LibFront::in(EMF_CALCULATOR, on) ) {
    return readEmfParameter(od, ECIF_CALCULATOR, param, true);
  }

  //-Physical Constant
  else if ( LibFront::in(EMF_CONSTANT, on) || LibFront::in("Constants", on) ) {
    return readEmfParameter(od, ECIF_CONSTANT, param, true);
  }

  //-Coordinate
  else if ( LibFront::in(EMF_COORDINATE, on) || LibFront::in("Coordinates", on) ) {
    return readEmfParameter(od, ECIF_COORDINATE, param, true);
  }

  //-Datafile
  else if ( LibFront::in(EMF_DATAFILE, on) || LibFront::in("Datafiles", on) ) {
    return readEmfParameter(od, ECIF_DATAFILE, param, true);
  }

  //-Equation
  else if ( LibFront::in(EMF_EQUATION, on) ) {
    return readEmfParameter(od, ECIF_EQUATION, param, true);
  }

  //-EquationVariables
  else if ( LibFront::in(EMF_EQUATION_VARIABLE, on) || LibFront::in("Equation Variables", on) ) {
    return readEmfParameter(od, ECIF_EQUATION_VARIABLE, param, true);
  }

  //-Header
  else if ( LibFront::in(EMF_HEADER, on) ) {
    return readEmfHeader(emfIn, od);
  }

  //-Grid H parameter
  else if ( LibFront::in(EMF_GRID_H, on) || LibFront::in("GridH", on) ) {
    return readEmfParameter(od, ECIF_GRID_H, param, true);
  }

  //-Grid parameter
  else if ( LibFront::in(EMF_GRID_PARAMETER, on) ||
            LibFront::in("Mesh Parameter", on)
    ) {
    return readEmfParameter(od, ECIF_GRID_PARAMETER, param, true);
  }

  //-Initial condition
  else if ( LibFront::in(EMF_INITIAL_CONDITION, on) ) {
    return readEmfParameter(od, ECIF_INITIAL_CONDITION, param, true);
  }

  //-Material
  else if ( LibFront::in(EMF_MATERIAL, on) ) {
    return readEmfParameter(od, ECIF_MATERIAL, param, true);
  }

  //-Model parameter
  else if ( LibFront::in(EMF_MODEL_PARAMETER, on) ) {
    return readEmfParameter(od, ECIF_MODEL_PARAMETER, param, true);
  }

  //-Simulation parameter
  else if ( LibFront::in(EMF_SIMULATION_PARAMETER, on) ) {
    return readEmfParameter(od, ECIF_SIMULATION_PARAMETER, param, true);
  }

  //-Solver
  else if ( LibFront::in(EMF_SOLVER, on) ) {
    return readEmfParameter(od, ECIF_SOLVER, param, true);
  }

  //-SolverControl parameter
  else if ( LibFront::in(EMF_SOLVER_CONTROL, on) ) {
    return readEmfParameter(od, ECIF_SOLVER_CONTROL, param, true);
  }

  //-Model statistics (not read actually!!!)
  else if ( LibFront::in(EMF_STATISTICS, on) ) {
    return isOk;
  }

  //-Update timestamps
  else if ( LibFront::in(EMF_TIMESTAMPS, on) ) {
    return readEmfTimestamps(od);
  }

  //-Timestep
  else if ( LibFront::in(EMF_TIMESTEP, on) ) {
    return readEmfParameter(od, ECIF_TIMESTEP, param, true);
  }

  //-User settings
  else if ( LibFront::in(EMF_USER_SETTING, on) || LibFront::in("User Settings", on) ) {
    return readEmfParameter(od, ECIF_USER_SETTING, param, true);
  }

  //-UNKNOWN object
  else {
    return unknownObjectMsg(od, false);
  }

  return isOk;
}


// This method gets one field at a time for the header
int
InputEmf::readEmfHeader(InputEmf* emfIn, emf_ObjectData_X* od)
{
  UserInterface* gui = theControlCenter->getGui();
  const char* fn = od->field_name;
  ModelInfo* mi = (ModelInfo*) model->getModelInfo();
  ParallelInfo* pi = (ParallelInfo*) model->getParallelInfo();

  //-Geometry edited flag for bodies
  if ( LibFront::in(EMF_BODY_GEOMETRY_EDITED, fn) ) {
    model->setFlagValue(GEOMETRY_EDITED, GEOMETRY_EDITED_BODIES, true);
  }

  //-Geometry edited flag for boundaries
  else if ( LibFront::in(EMF_BOUNDARY_GEOMETRY_EDITED, fn) ) {
    model->setFlagValue(GEOMETRY_EDITED, GEOMETRY_EDITED_BOUNDARIES, true);
  }
 
  //-Created string
  else if ( LibFront::in(EMF_CREATED, fn) ) {
    update_dyna_string(mi->created, (char*)od->data);
  }

  //-Modified string
  else if ( LibFront::in(EMF_MODIFIED, fn) ) {
    update_dyna_string(mi->modified, (char*)od->data);
  }

  //-User defined fields flag
  else if ( LibFront::in(EMF_HAS_USER_DEFINITIONS, fn) ) {
    mi->hasUserDefinitions = true;
  }

  //-Cad source file
  else if ( LibFront::in(EMF_CAD_SOURCE_FILE, fn) ) {
    update_dyna_string(mi->cadSourceFile, (char*)od->data);
  }

  //-Model geometry dimension
  else if ( LibFront::in(EMF_DIMENSION, fn) ) {
    int target;
    LibFront::setNumericData(target, 0);
    mi->dimension = (ecif_modelDimension)target;
    emfIn->modelDimension = mi->dimension;
  }

  //-Elmer Front emf-file format original source version (scalar 1 -->)
  else if ( LibFront::in(EMF_ELMER_FRONT_INPUT_VERSION, fn) ) {
    LibFront::setNumericData(mi->frontPreviousInputVersionNbr, 0);
  }

  //-Elmer Front emf-file format version (scalar 1 -->)
  else if ( LibFront::in(EMF_ELMER_FRONT_VERSION, fn) ) {
    LibFront::setNumericData(mi->frontInputVersionNbr, 0);
  }

  //-Emf Matc input file
  else if ( LibFront::in(EMF_MATC_FILE_EMF, fn) ) {
    char* filename = NULL;
    readName(od, filename);
    if ( filename != NULL && filename[0] != '\0' ) {
      model->readMatcFile(filename, "for emf-file", false);
      model->setMatcInputFileEmf(filename);
    }
    delete[] filename;
  }

  //-Sif Matc input file
  else if ( LibFront::in(EMF_MATC_FILE_SIF, fn) ) {
    char* filename = NULL;
    readName(od, filename);
    if ( filename != NULL && filename[0] != '\0' ) {
      model->readMatcFile(filename, "for sif-file", false);
      model->setMatcInputFileSif(filename);
    }
    delete[] filename;
  }

  //-Include path (if saved in model file)
  else if ( LibFront::in(EMF_INCLUDE_PATH, fn) ) {
    update_dyna_string(mi->includePath, (char*)od->data);
  }

  //-Log files directory (if saved in model file)
  else if ( LibFront::in(EMF_LOG_DIRECTORY, fn) ) {
    update_dyna_string(mi->temporaryFilesDirectory, (char*)od->data);
  }

  //-Mesh name (inputVerionNbr <=4 )
  else if ( LibFront::in("Mesh Name", fn) ) {

    int nof_meshes = 1; // only one string
    char** mesh_names = NULL;

    LibFront::setStringData(nof_meshes, mesh_names);
    model->setMeshNames(nof_meshes, mesh_names);

    for (int i = 0; i < nof_meshes; i++) {
      delete[] mesh_names[i];
    }
    delete[] mesh_names;
  }

  //-Mesh names (inputVerionNbr >= 5 )
  else if ( LibFront::in(EMF_MESH_NAMES, fn) ) {

    int nof_meshes = od->dimension1; // list of strings
    char** mesh_names = NULL;

    LibFront::setStringData(nof_meshes, mesh_names);
    model->setMeshNames(nof_meshes, mesh_names);

    for (int i = 0; i < nof_meshes; i++) {
      delete[] mesh_names[i];
    }
    delete[] mesh_names;
  }

  //-Current mesh index (inputVerionNbr >= 5)
  else if ( LibFront::in(EMF_CURRENT_MESH_INDEX, fn) ) {
    int mesh_index;
    LibFront::setNumericData(mesh_index, 0);
    model->setCurrentMeshIndex(mesh_index);
  }

  //-Model level meshFs
  else if ( LibFront::in(EMF_MESH_F, fn) ) {

    int size = od->data_length;
    double* mesh_f = new double[size];

    for (int i = 0; i < size; i++) 
      LibFront::setNumericData(mesh_f[i], i);

    model->setMeshFs(size, mesh_f);

    delete[] mesh_f;
  }

  //-Model level meshHs
  else if ( LibFront::in(EMF_MESH_H, fn) ) {

    int size = od->data_length;
    double* mesh_h = new double[size];

    for (int i = 0; i < size; i++) 
      LibFront::setNumericData(mesh_h[i], i);

    model->setMeshHs(size, mesh_h);

    delete[] mesh_h;
  }

  else if ( LibFront::in(EMF_MESH_BG_MESH_FILE_INDICES, fn) ) {

    int size = od->data_length;
    int* bg_mesh_file_indices = new int[size];

    for (int i = 0; i < size; i++) 
      LibFront::setNumericData(bg_mesh_file_indices[i], i);

    model->setMeshBgMeshFileIndices(size, bg_mesh_file_indices);

    delete[] bg_mesh_file_indices;
  }

  else if ( LibFront::in(EMF_MESH_BG_MESH_FILES, fn) ) {

    int size = od->dimension1; // list of strings
    char** bg_files = NULL;

    LibFront::setStringData(size, bg_files);
    model->setMeshBgMeshFiles(size, bg_files);

    for (int i = 0; i < size; i++) {
      if ( bg_files [i] != NULL ) {
        delete[] bg_files[i];
      }
    }
    delete[] bg_files;
  }

  //-Bg mesh active indicators
  else if ( LibFront::in(EMF_MESH_BG_MESH_ACTIVES, fn) ) {

    int size = od->data_length;
    int* bg_act = new int[size];

    for (int i = 0; i < size; i++) 
      LibFront::setNumericData(bg_act[i], i);

    model->setMeshBgMeshActives(size, bg_act);

    delete[] bg_act;
  }

  //-Bg mesh control file indicators
  else if ( LibFront::in(EMF_MESH_BG_MESH_CONTROLS, fn) ) {

    int size = od->data_length;
    int* bg_ctrl = new int[size];

    for (int i = 0; i < size; i++)
      LibFront::setNumericData(bg_ctrl[i], i);

    model->setMeshBgMeshControls(size, bg_ctrl);

    delete[] bg_ctrl;
  }

  //-Model level grid parameter id
  else if ( LibFront::in(EMF_GRID_PARAMETER, fn) ) {
    // Do nothing!
    //LibFront::setNumericData(mi->gridParamId, 0);
    //model->setGridParamId(mi->gridParamId);
  }

  //-Mesh result file
  else if ( LibFront::in(EMF_MESH_RESULT_FILE, fn) ) {
    update_dyna_string(mi->meshResultFile, (char*)od->data);
  }

  //-Mesh source file
  else if ( LibFront::in(EMF_MESH_SOURCE_FILE, fn) ) {
    update_dyna_string(mi->meshSourceFile, (char*)od->data);
  }

  //-Minimum edge size
  else if ( LibFront::in(EMF_MINIMUM_EDGE_SIZE, fn) ) {
    LibFront::setNumericData(mi->minEdgeSize, 0);
  }

  //-Model description
  else if ( LibFront::in(EMF_MODEL_DESCRIPTION, fn) ) {
    update_dyna_string(mi->modelDescription, (char*)od->data);
  }

  //-Model name
  else if ( LibFront::in(EMF_MODEL_NAME, fn) ) {
    update_dyna_string(mi->modelName, (char*)od->data);
  }

  //-Model source type (cad/mesh/both)
  else if ( LibFront::in(EMF_MODEL_SOURCE_TYPE, fn) ) {
    int target;
    LibFront::setNumericData(target, 0);
    mi->modelSourceType = (ecif_modelSource)target;
    mi->updateGeometryType();
  }

  //-Model status
  else if ( LibFront::in(EMF_MODEL_STATUS, fn) ) {
    LibFront::setNumericData(mi->modelStatus, 0);
  }

  //-Nof processors (parallel info)
  else if ( LibFront::in(EMF_NOF_PROCESSORS, fn) ) {
    LibFront::setNumericData(pi->nofProcessors, 0);
  }

  //-problem description
  else if ( LibFront::in(EMF_PROBLEM_DESCRIPTION, fn) ) {
    update_dyna_string(mi->problemDescription, (char*)od->data);
  }

  //-Model name
  else if ( LibFront::in(EMF_PROBLEM_NAME, fn) ) {
    update_dyna_string(mi->problemName, (char*)od->data);
  }

  //-Output files directory (if saved in model file)
  else if ( LibFront::in(EMF_RESULTS_DIRECTORY, fn) ) {
    update_dyna_string(mi->resultsDirectory, (char*)od->data);
  }

  //-Modelfile timestamp
  else if ( LibFront::in(EMF_TIMESTAMP, fn) ) {
    update_dyna_string(mi->modelFileTs, (char*)od->data);
  }

  //-Unknown field
  else {
    return unknownFieldMsg(od, false);
  }

  return isOk;
}


// Read one field for a active parameter
int
InputEmf::readEmfParameter(emf_ObjectData_X* od,
                           ecif_parameterType parameter_type,
                           const Parameter*& param,
                           bool add_to_model)
{
  ModelInfo* mi = (ModelInfo*) model->getModelInfo();

  // This stores data during multiple calls
  static Parameter* parameter = NULL;

  UserInterface* gui = theControlCenter->getGui();
  const char* fn = od->field_name;
  
  int parent_tag;
  int sub_parent_tag;
  int parent_type;

  // Nothing for the reference argument yet!
  //
  param = NULL; 

  // New parameter started
  if (od->is_object_start) {
    
    // Create parameter
    // NOTE: NO_INDEX parameter id prevents the last_id class
    // variable being updated, it should be update only when
    // a parameter is really added to the model
    //
    if ( add_to_model ) {
      parameter = model->createNewParameter(parameter_type, od->object_id);
    } else {
      parameter = model->createNewParameter(parameter_type, NO_INDEX);
    }

    if (parameter == NULL)
      return 1;
    
    if ( add_to_model ) {
      model->addParameter(parameter_type, parameter);
    }
  }

  // Parse fields
  // ============

  //-Parent tag for which parameter was created (version 9-)
  // "Parent Object" for versions 5-8
  if ( LibFront::in(EMF_PARENT, fn) ||
       LibFront::in("Parent Object", fn)   
     ) {
    LibFront::setNumericData(parent_tag, 0);
    parameter->setParentEmfTag(parent_tag);    
  }

  //-Sub parent tag for which parameter was created (version 9-)
  // Like body-layer for grid parameter
  if ( LibFront::in(EMF_SUB_PARENT, fn) ) {
    LibFront::setNumericData(sub_parent_tag, 0);
    parameter->setSubParentEmfTag(sub_parent_tag);    
  }

  //-Object type for which parameter was created (versions 9-)
  // "Parent Object Type" for versions 5-8
  else if ( LibFront::in(EMF_PARENT_TYPE, fn) ||
            LibFront::in("Parent Object Type", fn) 
          ) {
    // Object type code (enumrated number like OT_EDGE)
    if ( mi->frontInputVersionNbr < 7 ) {
      LibFront::setNumericData(parent_type, 0);
      parameter->setParentEmfType((objectType)parent_type);

    // Object type name (like "Edge")
    } else {
      parameter->setParentEmfType( model->objectName2Type((char*)od->data) );
    }
  }

  //-Body1 (version 4)
  else if ( LibFront::in("Body1", fn) ) {
    LibFront::setNumericData(parent_tag, 0);
    parameter->setParentEmfTag(parent_tag);    
    parameter->setParentEmfType(OT_BODY);    
  }

  //-Body2 (version 4)
  else if ( LibFront::in("Body2", fn) ) {
  }

  //-Name
  else if ( LibFront::in(EMF_NAME, fn) ) {
    parameter->setName((char*) od->data);    
  }

  //-Apply count (not in use!)
  else if ( LibFront::in("Apply Count", fn) ) {
  }

  //-Attach mode (not in use!)
  else if ( LibFront::in("Attach Mode", fn) ) {
  }

  //-Data (as string!)
  else if ( LibFront::in(EMF_DATA, fn) ) {
    parameter->setValue((char*) od->data);
  }

  //-Mask (is in use!)
  else if ( LibFront::in("Mask", fn) ) {
  }

  //-Unknown field
  else {
    // We accept here anything, because parameters
    // fields could come as fields!!!
    //return unknownFieldMsg(od, false);
  }

  // Parameter at end, store pointer to the reference argument
  //
  if (od->is_object_end) {
    param = parameter;
  }

  return isOk;
}


// This method gets one field for the Tisestamps data
int
InputEmf::readEmfTimestamps(emf_ObjectData_X* od)
{
  UserInterface* gui = theControlCenter->getGui();
  const char* fn = od->field_name;
  ModelInfo* mi = (ModelInfo*) model->getModelInfo();

  //-Elmer Front data itself
  if ( LibFront::in(EMF_TS_FRONT, fn) ) {
    update_dyna_string(mi->databaseTs, (char*) od->data);
  }

  //-Elmer Database
  else if ( LibFront::in(EMF_TS_DATABASE, fn) ) {
    update_dyna_string(mi->databaseTs, (char*) od->data);
  }

  //-Elmer Mesh Parameter 
  else if ( LibFront::in("Grid Parameter", fn) ) {
    update_dyna_string(mi->meshParameterTs, (char*) od->data);
  }

  //-Elmer Mesh generator 
  else if ( LibFront::in(EMF_TS_MESH, fn) ) {
    update_dyna_string(mi->meshTs, (char*) od->data);
  }

  //-Elmer Solver 
  else if ( LibFront::in(EMF_TS_SOLVER, fn) ) {
    update_dyna_string(mi->solverTs, (char*) od->data);
  }

  //-Elmer Gebhardt Factors
  else if ( LibFront::in(EMF_TS_GEBHARDT_FACTORS, fn) ||
            LibFront::in("GebhardtFactors", fn)
          ) {
    update_dyna_string(mi->gebhardtFactorsTs, (char*) od->data);
  }

  //-Elmer View Factors
  else if ( LibFront::in(EMF_TS_VIEW_FACTORS, fn) ||
            LibFront::in("ViewFactors", fn)
          ) {
    update_dyna_string(mi->viewfactorsTs, (char*) od->data);
  }

  //-Unknown field
  else {
    return unknownFieldMsg(od, false);
  }

  return isOk;
}
