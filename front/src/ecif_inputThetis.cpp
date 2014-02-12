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
Module:     ecif_inputThetis.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_bodyElement.h"
#include "ecif_control.h"
#include "ecif_geometry.h"
#include "ecif_inputThetis.h"
#include "ecif_model.h"
#include "ecif_userinterface.h"
#include "ecif_timer.h"

extern char read_buffer[];

InputThetis::InputThetis(enum ecif_modelDimension m_dim,
                         ifstream& in_file, char* in_filename):
Input(m_dim, in_file, in_filename)
{
}


enum ecif_modelDimension
InputThetis::findModelDimension(ifstream& in_file) {
  return modelDimension;
}


//---Create and update mesh tables
bool
InputThetis::processMeshFileData()
{
  Timer timer;
  double time, time1, time2;

  UserInterface* gui = theControlCenter->getGui();

  MeshElementTable* bt = model->getMeshBulkElements();
  MeshElementTable* bet = model->getMeshBoundaryElements();

  //---Create BULK element connections (neighor ids, sub element indices
  timer.start(); time1 = 0;
  gui->showMsg("---Creating volume element connections ...");

  time2 = timer.getLapTime(WALL_TIME); time = time2 - time1; time1 = time2;
  gui->showUsedTimeMsg(time, "---Creating volume element connections", (short)0, false);

  //---Create bulk element edges
  model->createMeshBulkElementEdges();

  //---Find mesh boundary element neighbors
  model->findMeshElementNeighbors(bet);

  //---Create boundary element edges
  model->createMeshBoundaryElementEdges();

  //---Finish
  timer.stop();
  time = timer.getEndTime(WALL_TIME);
  gui->showUsedTimeMsg(time, "Creating mesh bodies and boundaries", 1,  true);

	return modelDimension;
}


// Read nodes and elements from the Thetis Mesh File (.tmf)
bool
InputThetis::readMeshGeometry()
{
  int i;
  int nofNodes, nofElements;
  int nofBodies;
  int nofBoundaries;
  int nofBoundaryElements; // Fem elements
  Rc rc;

  // Read header information
  readFileLine(infile, read_buffer);
  istrstream strmline(read_buffer);
  if (! (strmline >> nofNodes >> maxExternalNodeId
                  >> nofElements >> maxExternalElementId
                  >> nofBodies
                  >> nofBoundaries
                  >> nofBoundaryElements
        )
     ) {
    cerr << endl << "ERROR: in reading mesh header info in InputThetis::readMeshGeoemtry";
    return false;
  }

  //---Allocate basic tables in the model
  model->allocateMeshBodies(nofBodies);
  model->allocateMeshNodes(nofNodes, maxExternalNodeId);
  model->allocateMeshInputElements(nofElements, maxExternalElementId);

  //---Read nodes from the mesh input file
  if ( ECIF_OK != readMeshNodes(nofNodes) )
    return false;

  //---Read volume elements from the mesh input file
  // NOTE: model dimension known only after reading
  // volume elements
  if ( ECIF_OK != readMeshElements(nofElements) ) return false;

  //---Read all boundary elements from the mesh input file
  // NOTE: We cannot read Thetis boundary elements, because we do not know the
  // material (body id) for them, tehy must be created from bulk faces!!!
  //
  //if ( ECIF_OK != readMeshBoundaryElements(nofBoundaryElements) ) return false;

  model->setModelDimension(modelDimension);

  if (modelDimension == ECIF_ND) return false;

  model->setMeshNodes();

  return true;
}



// Read a boundary element from the Thetis Mesh File and store it into the model.
bool
InputThetis::read_boundary()
{
  UserInterface* gui = theControlCenter->getGui();

  int i;
  bool result_code = true;

  int boundary_tag, external_body1_tag, external_body2_tag;
  int nof_bndr_elements;

  readFileLine(infile, read_buffer);
  istrstream strmline(read_buffer);

  if ( !(strmline >> boundary_tag
                  >> external_body1_tag >> external_body2_tag
                  >> nof_bndr_elements) ) {
    gui->errMsg(1, "could not read boundary data in InputThetis::read_boundaries()");
    return false;
  }

  int body1_tag = model->getBodyTagExt2Int(external_body1_tag);
  int body2_tag = model->getBodyTagExt2Int(external_body2_tag);

  int body1_lr = 0;
  int body2_lr = 0;

  // Check if boundary is already stored (as a normal model boundary!)
  BodyElement* be;

  if ( model->getDimension() == ECIF_2D )
    be = model->getEdgeByTag(boundary_tag);
  else
    be = model->getFaceByTag(boundary_tag);

  if (be != NULL) {
    be->allocateMeshElements(nof_bndr_elements);
  } else {
    model->createBodyElement(boundary_tag, body1_tag, body1_lr, body2_tag, body2_lr, nof_bndr_elements);
  }

  return result_code;
}


// Read a boundary element from the Thetis Mesh File and store it into the model.
Rc
InputThetis::readMeshBoundaryElements(int nof_bndr_elements)
{
  UserInterface* gui = theControlCenter->getGui();

  bool result_code = true;

  int sequence_nbr, bndr_tag;
  int ext_parent1_id, ext_parent2_id;
  int ext_body1_tag, ext_body2_tag;
  int elem_type, nof_nodes;
  meshElementCode elem_code;
  static int ext_node_ids[MAX_NOF_BNDR_NODES];

  char* strm_buffer = new char[BUFFER_LEN];
  strstream strm(strm_buffer, BUFFER_LEN, ios::app);

  short i;
  for (int index = 0; index < nof_bndr_elements; index++) {
    readFileLine(infile, read_buffer);
    reset(strm);
    strm << read_buffer << ends;

    if ( !(strm >> sequence_nbr
                >> bndr_tag
                >> ext_parent1_id
                >> ext_parent2_id
                >> ext_body1_tag
                >> ext_body2_tag
                >> elem_type
           )
       ) {
      gui->errMsg(1, "could not read element data-1 in InputThetis::read_boundary_element()");
      return ECIF_ERROR;
    }

    // Read nodes for the element (external ids!)
    if ( !(strm >> nof_nodes) ) {
      gui->errMsg(1, "could not read nof nodes per element in InputThetis::read_boundary_element()");
      return ECIF_ERROR;
    }

    // NOTE: For boundary elements we have to convert id "manually", because
    // they where not yet available when createMeshBodyTables (which
    // converts bulk elem- and node-ids) was called

    for (i = 0; i < nof_nodes; i++) {
      if (!(strm >> ext_node_ids[i])) {
        gui->errMsg(1, "could not read nodeId in InputThetis::read_boundary_element()");
        return ECIF_ERROR;
      }
    }

    // NOTE: We do not know the material id for the boundaries
    // so they cannot be actually added like this!
    // --> We do not read the boundary elements fromt the Thetis mesh file!!!
    model->addMeshInputElement(elem_type,
                               index, NO_INDEX, bndr_tag,
                               ext_node_ids);
    nofInputBoundaryElements++;

  } // for all boundary elements

  return ECIF_OK;
}


// Read an element from the Thetis Mesh File and store it into the model.
Rc
InputThetis::readMeshElements(int nof_elements)
{
  UserInterface* gui = theControlCenter->getGui();

  int external_id;
  int material_id;
  int elem_type;
  meshElementCode elem_code;
  int nof_nodes, nof_neighbors;
  static int ext_node_ids[27];
  static int ext_neighbor_ids[27];

  char* strm_buffer = new char[BUFFER_LEN];
  strstream strm(strm_buffer, BUFFER_LEN, ios::app);

  short i;
  for (int index = 0; index < nof_elements; index++) {
    readFileLine(infile, read_buffer);
    reset(strm);
    strm << read_buffer << ends;

    // Read element ids and code
    if ( !(strm >> external_id
                >> material_id
                >> elem_type
          )
       ) {
      gui->errMsg(1, "could not read element ids data in InputThetis::read_element()");
      return ECIF_ERROR;
    }

    // Read nodes for the element
    if ( !(strm >> nof_nodes) ) {
      gui->errMsg(1, "could not read nof nodes per element in InputThetis::read_element()");
      return ECIF_ERROR;
    }

    for (i = 0; i < nof_nodes; i++) {
      if ( !(strm >> ext_node_ids[i]) ) {
        gui->errMsg(1, "could not read nodeId in InputThetis::read_element()");
        return ECIF_ERROR;
      }
    }

    // Read neigbour element for the element
    if ( !(strm >> nof_neighbors) ) {
      gui->errMsg(1, "could not read nof-neighbors in InputThetis::read_element()");
      return ECIF_ERROR;
    }

    for (i = 0; i < nof_neighbors; i++) {
      if ( !(strm >> ext_neighbor_ids[i]) ) {
        gui->errMsg(1, "could not read neighborIds in InputThetis::read_element()");
        return ECIF_ERROR;
      }
    }

    model->addMeshInputElement(elem_type,
                               external_id, material_id, NO_INDEX,
                               ext_node_ids);
    nofInputBulkElements++;

  } // for all elements

	//---Model dimension (based on the last element read!)
  if (elem_type > 200 && elem_type < 300)
    modelDimension = ECIF_1D;

  else if (elem_type > 300 && elem_type < 500)
    modelDimension = ECIF_2D;

  else if (elem_type > 500 && elem_type < 900)
    modelDimension = ECIF_3D;

  else
    modelDimension = ECIF_ND;

  return ECIF_OK;
}


// Read a node from the Thetis Mesh File and store it into the model
Rc
InputThetis::readMeshNodes(int nof_nodes)
{
  UserInterface* gui = theControlCenter->getGui();

  int external_id;
  static Point3 p;

  // Init point buffer
  p[0] = p[1] = p[2] = 0.0;

  char* strm_buffer = new char[BUFFER_LEN];
  strstream strm(strm_buffer, BUFFER_LEN, ios::app);

  for (int index = 0; index < nof_nodes; index++) {

    readFileLine(infile, read_buffer);
    reset(strm);
    strm << read_buffer << ends;

    if ( !(strm >> external_id) ) {
      gui->errMsg(1, "could not NodeId in InputThetis::read_node()");
      return ECIF_ERROR;
    }

    for (int i = 0; i < 3; i++) {
      if (!(strm >> p[i])) {
        gui->errMsg(1, "could not read node coordinates in InputThetis::read_node()");
        return ECIF_ERROR;
      }
    }

    model->addMeshNode(index, external_id, p);
  }

  return ECIF_OK;
}


