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
Module:     ecif_inputEgf.cpp
Language:   C++
Date:       20.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation (read Elmer Geometry Format "cad" file)

************************************************************************/

#include "ecif_body2D.h"
#include "ecif_body3D.h"
#include "ecif_bodyElement.h"
#include "ecif_control.h"
#include "ecif_geometry.h"
#include "ecif_inputEgf.h"
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

InputEgf::InputEgf(enum ecif_modelDimension m_dim,
                       ifstream& in_file, char* in_filename):
InputFront(m_dim, in_file, in_filename)
{
  isEgfInput = true;
  modelDimension = m_dim;
  model->setModelDimension(modelDimension);
}


// Main callable function for the file reading
bool
InputEgf::readCadGeometry()
{
  struct emf_ObjectData_X my_ObjectData;
  emf_ObjectData  = &my_ObjectData;

  // Set callback functions
  //int my_readDataCB(void**);
  emf_readDataCB = readEgfGeometryCB;
  emf_readMessageCB = readEgfGeometryMsgCB;

  InputEgf* ci = this;
  char msg_buffer[256];

  // Start parser by calling the interface function
  int rc = emf_readData(infileName, (void**)&ci, 255, msg_buffer, 0);

  if (rc != 0) {
    modelDimension = ECIF_ND;
    return false;
  }

  // Convert tags to ids
  if ( !model->convertTags2Ids() ) return false;

  return true;
}


int
InputEgf::readEgfGeometryMsgCB(char* msg_buffer)
{
  UserInterface* gui = theControlCenter->getGui();
  return gui->showMsg(msg_buffer);
}


// This method gets one "field" at a time from the file
int
InputEgf::readEgfGeometryCB(void** user_data)
{
  UserInterface* gui = theControlCenter->getGui();

  InputEgf* egfIn = (InputEgf*)*user_data;
  InputFront* frontIn = (InputFront*)*user_data;

  // Global object pointing to the record, make
  // a short hand for it
  emf_ObjectData_X* od = emf_ObjectData;

  const char* on = od->object_name;
  int onl = od->object_name_length;
  const char* fn = od->field_name;
  int fnl = od->field_name_length;
  const char* dt = od->data_type;
  int dtl = od->data_type_length;

  int rc;

  // Read data according to the object type

  //-Header
  if ( LibFront::in("Header", on) ) {
    rc = readEgfHeader(egfIn, od);
  }

  //-Include command
  //else if ( LibFront::in("Include", on) ) {
  //  rc = readIncludeFile(egfIn, od);
  //}

  //-Vertex table
  else if ( LibFront::in(EMF_VERTICES, on) ||
            LibFront::in(EMF_VERTEX_TABLE, on)
          ) {
    rc = readVertexTable(od);
  }

  //-Edge (egf-format)
  else if ( LibFront::in("Edge", on) ) {
    rc = readEdge(od);
  }

  //-Element group
  else if ( LibFront::in("Face Group", on) ||
            LibFront::in("Edge Group", on) ||
            LibFront::in("Vertex Group", on)
          ) {
    rc = readElementGroup(od);
  }

  //-Edge loop
  else if ( LibFront::in("Edge Loop", on) ) {
    rc = readElementLoop(od);
  }

  //-Body
  else if ( LibFront::in("Body", on) ) {
    rc = readBody(od);
  }

  //-UNKNOWN object
  else {
    rc = unknownObjectMsg(od, IS_FATAL);
  }

  return rc;
}


// This method gets one field at a time for the header
int
InputEgf::readEgfHeader(InputEgf* egfIn, emf_ObjectData_X* od)
{
  UserInterface* gui = theControlCenter->getGui();
  const char* fn = od->field_name;
  ModelInfo* mi = (ModelInfo*) model->getModelInfo();

  //-Model geometry dimension
  if ( LibFront::in("Dimension", fn) ) {
    int target;
    LibFront::setNumericData(target, 0);
    mi->dimension = (ecif_modelDimension)target;
    egfIn->modelDimension = mi->dimension;
  }

  //-Geometry input unit  (1.0 (default) <--> m, 0.001 <--> mm etc)
  else if ( LibFront::in("Unit", fn) ) {
    LibFront::setNumericData(inputUnit, 0);
  }

  //-Model name
  else if ( LibFront::in("Model Name", fn) ) {
    readName(od, mi->modelName);
  }

  //-Model DB parent directory
  else if ( LibFront::in("Model Directory", fn) ) {
    //update_dyna_string(mi->modelDirectory, (char*)od->data);
    readName(od, mi->modelDirectory);
  }

  //-Emf Matc input file name
  else if ( LibFront::in(EMF_MATC_FILE, fn) ) {
    char* filename = NULL;
    readName(od, filename);
    if ( filename != NULL && filename[0] != '\0' ) {
      model->readMatcFile(filename, "for egf-file", false);
      model->setMatcInputFileEmf(filename);
    }
    delete[] filename;
  }

  //-Unknown field
  else {
    return unknownFieldMsg(od, IS_FATAL);
  }

  return isOk;
}


int
InputEgf::readIncludeFile(InputEgf* egfIn, emf_ObjectData_X* od)
{
  UserInterface* gui = theControlCenter->getGui();
  char* filename = NULL;

  readName(od, filename);

  if ( filename == NULL ) return isOk;

  // Check that input-file exists.
  if ( NULL == fopen(filename, "r") ) {
    gui->errMsg(0, "***ERROR Can't open the include file: ", filename);
    return false;
  }

  ifstream in_file(filename);

  InputEgf* in = new InputEgf(egfIn->getModelDimension(), in_file, filename);

  // Store curren object data pointer
  struct emf_ObjectData_X* tmp =  emf_ObjectData;

  in->readCadGeometry();

  // Restore object data pointer
  emf_ObjectData = tmp;


  return isOk;
}

