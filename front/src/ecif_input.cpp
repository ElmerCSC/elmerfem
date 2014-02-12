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
Module:     ecif_input.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_body.h"
#include "ecif_bodyElement.h"
#include "ecif_bodyElement1D.h"
#include "ecif_bodyElement2D.h"
#include "ecif_bodyElement3D.h"
#include "ecif_control.h"
#include "ecif_input.h"
#include "ecif_geometry.h"
#include "ecif_mesh.h"
#include "ecif_model.h"
#include "ecif_timer.h"
#include "ecif_userinterface.h"

Control* Input::theControlCenter = NULL;
Model* Input::model = NULL;


//Constructor
Input::Input(enum ecif_modelDimension m_dim,
             ifstream& in_file, char* filename)
             : inputDimension(m_dim), infile(in_file)
{
  modelDimension = ECIF_ND;

  bulkElementsAllocated = false;
  bndrElementsAllocated = false;
  edgeElementsAllocated = false;
  vrtxElementsAllocated = false;

  maxExternalElementId = -1;
  maxExternalNodeId = -1;

  nofInputBulkElements = 0;
  nofInputBoundaryElements = 0;
  nofInputEdgeElements = 0;
  nofInputVertexElements = 0;

  nofElements = 0;
  nofNodes = 0;

  nofBulkElements = 0;
  nofBoundaryElements = 0;
  nofEdgeElements = 0;
  nofVertexElements = 0;

  create_dyna_string(infileName, filename);
}


Input::~Input()
{
  delete[] infileName;
}


BodyElement*
Input::createBodyElement2D(Body* parent_body, int body_layer, Point3& p1, Point3& p2)
{
  static GcPoint points[2];
  static BodyElement* vertices[2];

  points[0].setPosit(p1[0], p1[1], p1[2]);
  points[1].setPosit(p2[0], p2[1], p2[2]);

  for (int i = 0; i < 2; i++) {

    //-Check if vertex already exists
    vertices[i] = model->findVertex(&points[i]);

    //-Create new vertex
    if ( vertices[i] == NULL ) {
      vertices[i] = new BodyElement1D(&points[i]);
      model->addBodyElement(vertices[i]);
    }
  }

  int v1_id = vertices[0]->Id();
  int v2_id = vertices[1]->Id();

  return createBodyElement2D(parent_body, body_layer, v1_id, v2_id);
}


BodyElement*
Input::createBodyElement2D(Body* parent_body, int body_layer, ecif_EdgeGeometry_X& edge)
{
  int i;
  double w;

  static GcPoint points[2];
  BodyElement* vertices[2];

  if ( edge.start != NULL ) {
    points[0].setPosit(edge.start[0][0], edge.start[0][1], edge.start[0][2]);
  }

  if ( edge.end != NULL ) {
    points[1].setPosit(edge.end[0][0], edge.end[0][1], edge.end[0][2]);
  }

  for (i = 0; i < 2; i++) {

    vertices[i] = model->findVertex(&points[i]);

    if ( vertices[i] == NULL ) {
      vertices[i] = new BodyElement1D(&points[i]);
      model->addBodyElement(vertices[i]);
    }
  }

  int v1_id = vertices[0]->Id();
  int v2_id = vertices[1]->Id();

  // Create a new 2D element into the body
  BodyElement* be = new BodyElement2D(v1_id, v2_id, &edge);

  parent_body->addElement(body_layer, be);
  model->addBodyElement(be);

  return be;
}


BodyElement*
Input::createBodyElement2D(Body* parent_body, int body_layer, int v1_id, int v2_id)
{
  // Create a new 2D line element into the body
  BodyElement* be = new BodyElement2D(v1_id, v2_id);

  parent_body->addElement(body_layer, be);
  model->addBodyElement(be);

  return be;
}


BodyElement*
Input::createBodyElement2D(Body* parent_body, int body_layer, int nof_vertices, int* vertex_ids)
{
  // Create a new 2D line element into the body
  BodyElement* be = new BodyElement2D(nof_vertices, vertex_ids);

  parent_body->addElement(body_layer, be);
  model->addBodyElement(be);

  return be;
}


BodyElement*
Input::createBodyElement2D(Body* parent_body, int body_layer, int nof_points, Point3* points)
{
  GcPoint point;

  int* vertex_ids = new int[nof_points];

  for (int i = 0; i < nof_points; i++) {

    point.setPosit(points[i][0], points[i][1], points[i][2]);

    //-Check if vertex already exists
    BodyElement* v = model->findVertex(&point);

    //-Create new vertex
    if ( v == NULL ) {
      v = new BodyElement1D(&point);
      model->addBodyElement(v);
    }

    vertex_ids[i] = v->Id();
  }

  return createBodyElement2D(parent_body, body_layer, nof_points ,vertex_ids);
}


bool
Input::doBreak()
{
  return theControlCenter->getBreakValue(MESH_INPUT);
}


// Remove leading white space from stream
void
Input::eat_white(istream& strm)
{
  char c = strm.get();

  while ( c == ' ' ||  c == '\t'  ||  c == '\r' || c == '\n' ){
    c = strm.get();
  }

  strm.putback(c);
}


// Helper function for (mesh) models
// Conclude model dimension from geometry (boundíng box) dimension and
// element type dimension
//
ecif_modelDimension
Input::findModelDimension(ecif_modelDimension geom_dim)
{
  // If geometry is 3D, but only boundary elements, then we
  // have a boundary (BEM) element model
  if ( inputDimension == ECIF_3D ||
       geom_dim == ECIF_3D
     ) {
    return ECIF_3D;
  } else {
    return ECIF_2D;
  }
}



void
Input::initClass(Model* mdl)
{
  Input::model = mdl;
}


//---Create mesh tables
bool
Input::processMeshFileData()
{
  Rc rc;
  Timer timer;
  double time, time1, time2;
  UserInterface* gui = theControlCenter->getGui();

  MeshElementTable* bt = NULL;
  MeshElementTable* bet = NULL;

  //---We can allocate and create actual bulk elements
  if ( !bulkElementsAllocated ) {
    model->allocateMeshBulkElements(nofInputBulkElements, maxExternalElementId);
    bulkElementsAllocated = true;
  }

  rc = model->installMeshInputBulkElements();

  if ( rc != ECIF_OK ) return false;

  //---Construct Body2Material table, change body ids to
  // internal ids in the meshBulkElements table
  //
  model->createMeshBodyTables();
  model->createMeshBodies();
  model->convertMeshBulkElementIdsExt2Int();


  // Bulk element stuff
  // ==================
  bt = model->getMeshBulkElements();

  //---Create BULK element connections (neighbor ids)
  timer.start(); time1 = 0;
  gui->showMsg("---Creating volume element connections ...");

  model->findMeshElementNeighbors(bt);

  //---Create bulk element edges
  model->createMeshBulkElementEdges();

  time2 = timer.getLapTime(WALL_TIME); time = time2 - time1; time1 = time2;
  gui->showUsedTimeMsg(time, "---Creating volume element connections", short(0), false);


  // Boundary element stuff
  // ======================

  int nof_bulk_bndr_elems;
  model->findNofBulkBoundaryElements(nof_bulk_bndr_elems);

  // Message
  time2 = timer.getLapTime(WALL_TIME); time = time2 - time1; time1 = time2;
  gui->showMsg("---Creating mesh boundary elements...");

  // A simple model: all boundaries are defined by bulk faces/edges
  //
  if ( nofInputBoundaryElements == 0 ) {

    if ( !bndrElementsAllocated ) {
      model->allocateMeshBoundaryElements(nof_bulk_bndr_elems);
      bndrElementsAllocated = true;
    }

    model->createMeshBoundaries();
    bet = model->getMeshBoundaryElements();

  // Some boundaries were given in the input, so these elements must be first matched
  // with the bulk faces/edges and then the rest of the bulk faces define other
  // boundary elements
  //
  } else {

    if ( !bndrElementsAllocated ) {
      model->allocateMeshBoundaryElements(nofInputBoundaryElements);
      bndrElementsAllocated = true;
    }

    model->installMeshInputBoundaryElements();
    bet = model->getMeshBoundaryElements();

    bool* free_bulk_bndr_elems = new bool[nof_bulk_bndr_elems];
    model->findMeshBoundaryParents(nof_bulk_bndr_elems, free_bulk_bndr_elems);

    // Allocate more space in boundary element table
    int count = 0;
    for (int i = 0; i < nof_bulk_bndr_elems; i++) {
      if ( !free_bulk_bndr_elems[i] ) continue;
      count++;
    }

    if ( count > 0 ) {
      model->reallocateMeshBoundaryElements(nofInputBoundaryElements + count);
    }
    
    model->createMeshBoundaries(nof_bulk_bndr_elems, free_bulk_bndr_elems);
    delete[] free_bulk_bndr_elems;
  }

  // Check if input elements contain some elements which are not installed into any
  // boundary. These will become Bem boundaries
  //
  model->checkMeshInputBoundaryElements();

  //---Create mesh boundary element neighbors
  model->findMeshElementNeighbors(bet);

  time2 = timer.getLapTime(WALL_TIME); time = time2 - time1; time1 = time2;
  gui->showUsedTimeMsg(time, "---Creating", nofBoundaryElements, "mesh boundary elements",
                         0, false);

  //---Create boundary element edges
  model->createMeshBoundaryElementEdges();

  // We can now delete all input elements (they are now all converted to actual
  // elements!)
  model->removeMeshInputElements();

  //---Finish
  timer.stop();
  time = timer.getEndTime(WALL_TIME);
  gui->showUsedTimeMsg(time, "Creating mesh bodies and boundaries", 0,  true);

	return true;
}


ecif_modelDimension
Input::loadMesh()
{

  return  readMeshFile();

#if 0
  try
  {
    return  readMeshFile();
  }

  catch (...)
  {
    UserInterface* gui = theControlCenter->getGui();

    if ( gui != NULL ) {
      strstream strm;
      strm << "***ERROR Unable to read mesh file " << infileName << ends;
      gui->showMsg(strm.str());
    }
  }
#endif

  // No success
  return ECIF_ND;
}


//Open given inputfile and read bodies into *model*.
enum ecif_modelDimension
Input::readCadFile()
{
  // For the very first, check that file makes sense!
  if (!readCadHeader()) {
    modelDimension = ECIF_ND;
    infile.close();
    return modelDimension;
  }

  strstream strm;
  UserInterface* gui = theControlCenter->getGui();
  gui->showMsg("---Reading the geometry ...");
  reset(strm);
  strm << infileName << ends;
  gui->showMsg(strm.str(), 1);

  // Ok, try to read all bodies from the CAD-file
  if (!readCadGeometry()) {
    modelDimension = ECIF_ND;
    infile.close();
    return modelDimension;
  }

  infile.close();
  model->setModelDimension(modelDimension);

  return modelDimension;
}


//Read one line from in_file into buffer.
//Line to be read is in position linepos (so 1 means next line)
//End-of-file condition should be checked before calling this method!!
//Comment lines (<==> first non-space is # in the line) are skipped
void
Input::readFileLine(ifstream& in_file, char* buffer, int linepos)
{
  int i, j, buffer_len;

  for (i = 1; i <= linepos; i++) {

    buffer[0] = '\0';

    if (in_file.eof()) {
      return;
    }

    in_file.getline(buffer, MAXLINE, '\n');

    buffer[MAXLINE] = '\0';

    buffer_len = strlen(buffer);

    // Change tabs and carrage returns to spaces
    for (j = 0; j < buffer_len; j++) {

      if (buffer[i] == '\t' || buffer[i] == '\r' )
        buffer[i] = ' ';
    }

    // Check if this is a comment line
    for (j = 0; j < buffer_len; j++) {

      // Skip leading spaces
      if (buffer[j] == ' ') {
        continue;
      }

      // If first non-space is comment character, skip the line
      if (buffer[j] == '#') {
        i--;
      }
      // Anyway, we must stop checking characters now
      break;
    }
  }
}


ecif_modelDimension
Input::readMeshFile()
{
  UserInterface* gui = theControlCenter->getGui();

  Timer timer;
  strstream strm;
  ecif_modelDimension model_dim;

  timer.start();
  gui->showMsg("Reading the mesh ...");
  reset(strm);
  strm << "***   " << infileName << ends;
  gui->showMsg(strm.str(), 1);

  // No success or interrupted
  //
  if (!readMeshGeometry()) {

    modelDimension = ECIF_ND;

    if ( doBreak() ) {
      gui->showMsg("***NOTE: Mesh reading interrupted!");
    } else {
      gui->showMsg("Could not read the mesh file!");
    }

    infile.close();
    return modelDimension;
  }

  timer.stop();
  double time = timer.getEndTime(WALL_TIME);
  gui->showUsedTimeMsg(time, "Reading the mesh file");

  model->setModelDimension(modelDimension);

  infile.close();

  return modelDimension;
}


void
Input::reset(strstream& strm)
{
  strm.clear();
  strm.seekp(0);
  strm.seekg(0);
}


void
Input::setModelDimension(int dimension)
{
  if (dimension == 1)
    modelDimension = ECIF_1D;
  else if (dimension == 2)
    modelDimension = ECIF_2D;
  else if (dimension == 3)
    modelDimension = ECIF_3D;
}
