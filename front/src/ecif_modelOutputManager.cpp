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
Module:     ecif_modelOutputManager.cpp
Language:   C++
Date:       13.04.00
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include <eio_api.h>
#include "ecif_body.h"
#include "ecif_body2D.h"
#include "ecif_body3D.h"
#include "ecif_bodyLayer.h"
#include "ecif_bodyElement.h"
#include "ecif_bodyElement1D.h"
#include "ecif_bodyElement2D.h"
#include "ecif_bodyElement3D.h"
#include "ecif_bodyElementLoop.h"
#include "ecif_bodyElementGroup.h"
#include "ecif_boundaryCondition.h"
#include "ecif_boundbox.h"
#include "ecif_control.h"
#include "ecif_input.h"
#include "ecif_mesh.h"
// #include "ecif_modelOutputManager.h"
#include "ecif_model_aux.h"
#include "ecif_parameter.h"
#include "ecif_parameterField.h"
#include "ecif_timer.h"
#include "ecif_userinterface.h"


Control* ModelOutputManager::theControlCenter = NULL;
Model* ModelOutputManager::model = NULL;

ModelOutputManager::ModelOutputManager()
{
}


ModelOutputManager::~ModelOutputManager()
{
}



//**************************************************************
//*** EMF output function (Elmer Front's model file format)  ***
//**************************************************************

// Save Elmer Front Model File (emf-file)
//
ostream&
ModelOutputManager::emf_output(ostream& out, char* filename)
{
  const ModelInfo* mi = model->getModelInfo();

  // Info section as comments
  out << "!ElmerFront model file" << endl;
  out << "!Saved        = " << mi->modified << endl;
  out << "!Case         = " << mi->modelName << "  " << mi->problemName << endl;
  out << "!Model dir    = " << mi->modelDirectory_absolute << endl;
  out << "!Include path = " << mi->includePath_absolute << endl;
  out << "!Results dir  = " << mi->resultsDirectory_absolute << endl;
  out << "!Log dir      = " << mi->temporaryFilesDirectory_absolute << endl;

  out << endl;

  emf_outputHeader(out);
  emf_outputTimestamps(out);
  emf_outputStatistics(out);
  emf_outputParameters(out, ECIF_MODEL_PARAMETER);
  emf_outputParameters(out, ECIF_SIMULATION_PARAMETER);
  emf_outputParameters(out, ECIF_SOLVER_CONTROL);
  emf_outputParameters(out, ECIF_CONSTANT);
  emf_outputParameters(out, ECIF_COORDINATE);
  emf_outputParameters(out, ECIF_DATAFILE);
  emf_outputParameters(out, ECIF_EQUATION_VARIABLE);
  emf_outputParameters(out, ECIF_EQUATION);
  emf_outputParameters(out, ECIF_SOLVER);
  emf_outputParameters(out, ECIF_CALCULATOR);
  emf_outputParameters(out, ECIF_TIMESTEP);

  emf_outputVertexTable(out);
  emf_outputElements(out);
  emf_outputElementGroups(out);
  emf_outputElementLoops(out);
  emf_outputBodies(out);
  emf_outputParameters(out, ECIF_BODY_PARAMETER);
  emf_outputParameters(out, ECIF_BOUNDARY_PARAMETER);
  emf_outputParameters(out, ECIF_INITIAL_CONDITION);
  emf_outputBoundaryConditions(out);
  emf_outputParameters(out, ECIF_MATERIAL);
  emf_outputParameters(out, ECIF_BODY_FORCE);
  emf_outputParameters(out, ECIF_GRID_PARAMETER);
  emf_outputParameters(out, ECIF_GRID_H);
  emf_outputParameters(out, ECIF_USER_SETTING);
  emf_outputEof(out);

  UserInterface* gui = theControlCenter->getGui();

  if (mi->gebhardtFactorsNeedsUpdate)
    gui->setNeedsUpdate("Gebhardt Factors");

  if (mi->meshNeedsUpdate)
    gui->setNeedsUpdate("Mesh");

  if (mi->solverNeedsUpdate)
    gui->setNeedsUpdate("Solver");

  return out;
}


// Output bodies.
ostream&
ModelOutputManager::emf_outputBodies(ostream& out)
{
  int index = 0;
  while (true) {
    Body* body = model->getBody(index++);
    if (body==NULL) break;
    emf_outputSectionStart(out);
    body->output_emf(out, ESF_INDENT_SIZE, 0);
    emf_outputSectionEnd(out);
  }

  return out;
}


// Output all boundary conditions (outer, inner).
// NOTE Boundary conditions need extra checking, they
// are not output as standard parameters
ostream&
ModelOutputManager::emf_outputBoundaryConditions(ostream& out)
{
  int index = 0;

  while (true) {

    Parameter* cond = model->getParameter(index++, ECIF_BOUNDARY_CONDITION);

    if (cond==NULL) break;

    // Output
    emf_outputSectionStart(out);
    cond->output_emf(out, ESF_INDENT_SIZE, 0, true, true, false);
    emf_outputSectionEnd(out);

    // Reset
    cond->resetValue();
  }

  return out;
}


// Output boundary groups.
ostream&
ModelOutputManager::emf_outputElementGroups(ostream& out)
{
  int index = 0;
  while (true) {
    BodyElementGroup* beg = model->getBodyElementGroup(index++);
    if (beg==NULL) break;
    if ( IMPLICIT_GROUP == beg->getGroupType() ) continue;
    emf_outputSectionStart(out);
    beg->output_emf(out, ESF_INDENT_SIZE, 0);
    emf_outputSectionEnd(out);
  }

  return out;
}


// Output all model elements
ostream&
ModelOutputManager::emf_outputElements(ostream& out)
{
  emf_outputVertices(out);
  emf_outputEdges(out);
  emf_outputFaces(out);

  return out;
}


// Output all model geometry vertices
ostream&
ModelOutputManager::emf_outputVertices(ostream& out)
{
  int index = 0;
  while (true) {
    BodyElement* v = model->getVertex(index++);
    if (v==NULL) break;
    emf_outputSectionStart(out);

    // Do not output geometry for a vertex given in the model's vertex-table
    if ( model->isInVertexTable(v->Id()) ) {
      v->output_emf(out, ESF_INDENT_SIZE, 0, false);
    } else {
      v->output_emf(out, ESF_INDENT_SIZE, 0, true);
    }
    emf_outputSectionEnd(out);
  }

  return out;
}


// Output all model geometry edges
ostream&
ModelOutputManager::emf_outputEdges(ostream& out)
{
  const ModelInfo* mi = model->getModelInfo();

  double start[3] = {0.0, 0.0, 0.0};
  double end1[3] = {0.0, 0.0, 0.0};
  double end2[3] = {0.0, 0.0, 0.0};

  bool checkSymmetry = model->getSymmetryAxis(start, end1, end2) && (mi->dimension == ECIF_2D);

  GcPoint start_p(start);
  GcPoint end_p1(end1);
  GcPoint end_p2(end1);

  int index = 0;
  while (true) {
    BodyElement* e = model->getEdge(index++);
    if (e==NULL) break;

    //-Check status
    beStatus e_stat = e->getStatus();
    if ( e_stat & (BE_DEVIDED | BE_SWAPPED) ) {
      continue;
    }

    emf_outputSectionStart(out);
    e->output_emf(out, ESF_INDENT_SIZE, 0);
    emf_outputSectionEnd(out);
  }

  return out;
}


// Output all model geometry faces
ostream&
ModelOutputManager::emf_outputFaces(ostream& out)
{
  const ModelInfo* mi = model->getModelInfo();

  double start[3] = {0.0, 0.0, 0.0};
  double end1[3] = {0.0, 0.0, 0.0};
  double end2[3] = {0.0, 0.0, 0.0};

  bool checkSymmetry = model->getSymmetryAxis(start, end1, end2) && (mi->dimension == ECIF_3D);

  GcPoint start_p(start);
  GcPoint end_p1(end1);
  GcPoint end_p2(end1);

  int index = 0;
  while (true) {

    BodyElement* f = model->getFace(index++);
    if (f==NULL) break;

    //-Check status
    beStatus f_stat = f->getStatus();
    if ( f_stat & (BE_DEVIDED | BE_SWAPPED) ) {
      continue;
    }

    emf_outputSectionStart(out);
    f->output_emf(out, ESF_INDENT_SIZE, 0);
    emf_outputSectionEnd(out);
  }

  return out;
}


#if 0
// Output all model elements
ostream&
ModelOutputManager::emf_outputElements(ostream& out)
{
  double start[3] = {0.0, 0.0, 0.0};
  double end1[3] = {0.0, 0.0, 0.0};
  double end2[3] = {0.0, 0.0, 0.0};

  bool checkSymmetry = model->getSymmetryAxis(start, end1, end2);

  GcPoint start_p(start);
  GcPoint end_p1(end1);
  GcPoint end_p2(end1);

  BodyElement* be = model->getFirstBoundary();

  if (be == NULL)
    return out;

  while (be) {

    //-Check status
    beStatus be_stat = be->getStatus();
    if ( be_stat & (BE_DEVIDED | BE_SWAPPED) ) {
      be = getNextBoundary();
      continue;
    }

    // We check if outer boundary is aint possible symmetry axis (2D) or
    // symmetry plane (3D).
    bool isOnSymmetry = false;

    if (checkSymmetry) {
      if (modelInfo->dimension == ECIF_2D)
        isOnSymmetry = be->isOnSameAxis(start_p, end_p1);
      else if (modelInfo->dimension == ECIF_3D)
        isOnSymmetry = be->isOnSamePlane(start_p, end_p1, end_p2);
    }

    emf_outputSectionStart(out);
    be->output_emf(out, ESF_INDENT_SIZE, 0, isOnSymmetry);
    emf_outputSectionEnd(out);

    be = model->getNextBoundary();
  }

  return out;
}

#endif



// Output all model elements loops
ostream&
ModelOutputManager::emf_outputElementLoops(ostream& out)
{
  int index = 0;

  while (true) {

    BodyElementLoop* bel = model->getBodyElementLoop(index++);

    if (bel==NULL) break;

    emf_outputSectionStart(out);
    bel->output_emf(out, ESF_INDENT_SIZE, 0);
    emf_outputSectionEnd(out);
  }

  return out;
}



// Mark End Of File.
ostream&
ModelOutputManager::emf_outputEof(ostream& out)
{
  emf_outputSectionStart(out);

  out << "!End Of File" << endl;

  return out;
}


// Put model header data into ouput file.
ostream&
ModelOutputManager::emf_outputHeader(ostream& out)
{
  char* QM = "\"";     // quote-mark
  char* ES = "\"\"";  // empty string
  short IS = ESF_INDENT_SIZE;

  const ModelInfo* mi = model->getModelInfo();
  const ParallelInfo* pi = model->getParallelInfo();

  emf_outputSectionStart(out);

  // Header
  LibFront::output_string(out, IS, 0, EMF_HEADER, true);

  LibFront::output_scalar(out, IS, 1, EMF_CREATED, NULL, mi->created, true);
  LibFront::output_scalar(out, IS, 1, EMF_MODIFIED, NULL, mi->modified, true);

  if ( mi->hasUserDefinitions ) {
    LibFront::output_scalar(out, IS, 1, EMF_HAS_USER_DEFINITIONS, NULL, mi->hasUserDefinitions);
  }

  // Input version number for version control
  // ========================================
  //--Store input version number if it is different from the current  program version
  if ( mi->frontInputVersionNbr > 0 &&
       mi->frontInputVersionNbr != mi->frontVersionNbr
     ) {
    LibFront::output_scalar(out, IS, 1, EMF_ELMER_FRONT_INPUT_VERSION, NULL, mi->frontInputVersionNbr);
  //--Otherwise store previous input version number if it was given in the input
  } else if (mi->frontPreviousInputVersionNbr > 0 ) {
    LibFront::output_scalar(out, IS, 1, EMF_ELMER_FRONT_INPUT_VERSION, NULL, mi->frontPreviousInputVersionNbr);
  }

  LibFront::output_scalar(out, IS, 1, EMF_ELMER_FRONT_VERSION, NULL, mi->frontVersionNbr);
  LibFront::output_scalar(out, IS, 1, EMF_TIMESTAMP, NULL, mi->modelFileTs, true);
  LibFront::output_scalar(out, IS, 1, EMF_MODEL_STATUS, NULL, mi->modelStatus);
  LibFront::output_scalar(out, IS, 1, EMF_MODEL_SOURCE_TYPE, NULL, mi->modelSourceType);

  bool cad_src = false;
  bool msh_src = false;

  switch (mi->modelSourceType) {
  case ECIF_CAD_FILE:
    cad_src = true; break;
  case ECIF_MESH_FILE:
    msh_src = true; break;
  case ECIF_CAD_AND_MESH_FILE:
    cad_src = msh_src = true; break;
  }

  if (cad_src) {
    LibFront::output_scalar(out, IS, 1, EMF_CAD_SOURCE_FILE, NULL, mi->cadSourceFile, true);
  }
  if (msh_src) {
    LibFront::output_scalar(out, IS, 1, EMF_MESH_SOURCE_FILE, NULL, mi->meshSourceFile, true);
  }
  LibFront::output_scalar(out, IS, 1, EMF_MESH_RESULT_FILE, NULL, mi->meshResultFile, true);
  LibFront::output_scalar(out, IS, 1, EMF_MODEL_NAME, NULL, mi->modelName, true);
  LibFront::output_scalar(out, IS, 1, EMF_PROBLEM_NAME, NULL, mi->problemName, true);
  LibFront::output_scalar(out, IS, 1, EMF_MODEL_DESCRIPTION, NULL, mi->modelDescription, true);
  LibFront::output_scalar(out, IS, 1, EMF_PROBLEM_DESCRIPTION, NULL, mi->problemDescription, true);

  LibFront::output_scalar(out, IS, 1, EMF_MATC_FILE_EMF, NULL, mi->matcInputFile_emf, true);
  LibFront::output_scalar(out, IS, 1, EMF_MATC_FILE_SIF, NULL, mi->matcInputFile_sif, true);

  //--If "save in model file"-flag is turned on (in Gui)
  if ( mi->includePath_save )
    LibFront::output_scalar(out, IS, 1, EMF_INCLUDE_PATH, NULL, mi->includePath, true);
  if ( mi->resultsDirectory_save )
    LibFront::output_scalar(out, IS, 1, EMF_RESULTS_DIRECTORY, NULL, mi->resultsDirectory, true);
  if ( mi->temporaryFilesDirectory_save )
    LibFront::output_scalar(out, IS, 1, EMF_LOG_DIRECTORY, NULL, mi->temporaryFilesDirectory, true);

  LibFront::output_scalar(out, IS, 1, EMF_NOF_PROCESSORS, NULL, pi->nofProcessors);
  LibFront::output_scalar(out, IS, 1, EMF_DIMENSION, NULL, mi->dimension);
  LibFront::output_scalar(out, IS, 1, EMF_MINIMUM_EDGE_SIZE, NULL, mi->minEdgeSize);

  if ( model->getFlagValue(GEOMETRY_EDITED_BODIES) ) {
    LibFront::output_scalar(out, IS, 1, EMF_BODY_GEOMETRY_EDITED, NULL, 1);
  }

  if ( model->getFlagValue(GEOMETRY_EDITED_BOUNDARIES) ) {
    LibFront::output_scalar(out, IS, 1, EMF_BOUNDARY_GEOMETRY_EDITED, NULL, 1);
  }

  //--Mesh names. These define targets for body grid paramters, boundary mesh-values etc.!
  if ( mi->nofMeshes > 0 ) {
    LibFront::output_vector(out, IS, 1, EMF_MESH_NAMES, "String", mi->nofMeshes, (const char**)mi->meshNames, true, true);
  }

  LibFront::output_scalar(out, IS, 1, EMF_CURRENT_MESH_INDEX, NULL, mi->currentMeshIndex);

  //--Mesh control values and background mesh stuff ('mesh lists')
  if ( cad_src && mi->nofMeshes > 0 ) {
    LibFront::output_vector(out, IS, 1, EMF_MESH_H, NULL, mi->nofMeshes, mi->meshHs, false);
    LibFront::output_vector(out, IS, 1, EMF_MESH_F, NULL, mi->nofMeshes, mi->meshFs, false);
  }

   //--Mesh background mesh stuff
  if ( mi->nofBgMeshFiles > 0 ) {
    LibFront::output_vector(out, IS, 1, EMF_MESH_BG_MESH_FILE_INDICES, NULL,  mi->nofBgMeshFiles, mi->meshBgMeshFileIndices, false);
    LibFront::output_vector(out, IS, 1, EMF_MESH_BG_MESH_FILES, "String", mi->nofBgMeshFiles, (const char**)mi->meshBgMeshFiles, true, true);
    LibFront::output_vector(out, IS, 1, EMF_MESH_BG_MESH_ACTIVES, NULL, mi->nofBgMeshFiles, mi->meshBgMeshActives, false);
    LibFront::output_vector(out, IS, 1, EMF_MESH_BG_MESH_CONTROLS, NULL, mi->nofBgMeshFiles, mi->meshBgMeshControls, false);
  }

  emf_outputSectionEnd(out);

  return out;
}


// Output standard parameters.
ostream&
ModelOutputManager::emf_outputParameters(ostream& out, ecif_parameterType param_type)
{
  const ModelInfo* mi = model->getModelInfo();

  int index = 0;

  while (true) {

    Parameter* param = model->getParameter(index++, param_type);

    if (param==NULL) break;

    // Output
    emf_outputSectionStart(out);
    param->output_emf(out, ESF_INDENT_SIZE, 0, true, true, false);
    emf_outputSectionEnd(out);

    // Reset
    param->resetValue();
  }

  return out;
}


#if 0
// Output all model vertex-points.
ostream&
ModelOutputManager::emf_outputVertices(ostream& out)
{

  if (modelStatistics->nofVertices == 0)
    return out;

  short is = ESF_INDENT_SIZE;

  BodyElement* v;

  emf_outputSectionStart(out);

  LibFront::output_string(out, is, 0, "Vertices");

  // Variable: Vertex id
  // Data: Coordinate point
  int dim = 3; // Always 3D
  int data_size = dim + 2; // Two ids; Bndr tag + BC-id

  //LibFront::getFieldInfo(EFN_POINTS_AND_CONSTRAINT_IDS)->output_name(out, is, 1) << endl;
  LibFront::indent(out, is, 1) << "Ids And Points" << endl;
  LibFront::indent(out, is, 2) << "Variable Index" << endl;
  LibFront::indent(out, is, 3) << "Size " << data_size << endl;
  LibFront::indent(out, is, 4) << "Real" << endl;

  v = getFirstVertex();

  while ( v != NULL ) {
    LibFront::indent(out, is, 5);

    // Vertex tag
    out << v->Tag() << "  ";

    // Boundary tag
    out << v->getBoundaryTag() << " ";

    // Boundary Condition Id
    out << v->getBoundaryConditionId() << "  ";

    // Geometric point data
    GcPoint* point = (GcPoint*)v->getGeometry();
    out << *point;

    out << endl;

    v = getNextVertex();
  }

  // Data end marker
  LibFront::indent(out, is, 4) << "End" << endl;

  emf_outputSectionEnd(out);

  return out;
}
#endif


// Data ouput after a section
ostream&
ModelOutputManager::emf_outputSectionEnd(ostream& out)
{
  out << "End" << endl << endl;

  return out;
}


// Data ouput before a section
ostream&
ModelOutputManager::emf_outputSectionStart(ostream& out)
{
  //Currently nothing!
  return out;
}


// Put statistics data into model output file.
ostream&
ModelOutputManager::emf_outputStatistics(ostream& out)
{
  short is = ESF_INDENT_SIZE;

  const ModelStatistics* ms = model->getModelStatistics();

  emf_outputSectionStart(out);

  // Header
  LibFront::output_string(out, is, 0, "Statistics");

  LibFront::output_scalar(out, is, 1, "Nof Bodies", NULL, ms->nofBodies);
  LibFront::output_scalar(out, is, 1, "Nof Loops", NULL, ms->nofElementLoops);

  // Calc nof body_elements to be output
  int nof_elements = 0;

  int index = 0;

  while (true) {

    BodyElement* be = model->getBoundary(index++);

    if (be==NULL) break;

    //-Check status
    beStatus be_stat = be->getStatus();

    if ( !(be_stat & (BE_DEVIDED | BE_SWAPPED)) )
      nof_elements++;
  }

  LibFront::output_scalar(out, is, 1, "Nof Elements", NULL, nof_elements);
  LibFront::output_scalar(out, is, 1, "Nof Outer Boundaries", NULL, ms->nofOuterBoundaries);
  LibFront::output_scalar(out, is, 1, "Nof Inner Boundaries", NULL, ms->nofInnerBoundaries);
  LibFront::output_scalar(out, is, 1, "Nof Vertices", NULL, ms->nofVertices);

  // Calc max loop size
  int max_loop_size = 0;

  index = 0;

  while (true) {

    BodyElementLoop* bel = model->getBodyElementLoop(index++);

    if (bel==NULL) break;

    int size =  bel->getNofElements();

    if (size > max_loop_size)
      max_loop_size = size;

  }

  LibFront::output_scalar(out, is, 1, "Max Loop Count", NULL, max_loop_size);

  emf_outputSectionEnd(out);
  return out;
}


// Put model timestamps data into ouput file.
ostream&
ModelOutputManager::emf_outputTimestamps(ostream& out)
{
  char* QM = "\"";     // quote-mark
  short IS = ESF_INDENT_SIZE;

  const ModelInfo* mi = model->getModelInfo();

  emf_outputSectionStart(out);

  // Timestamps when data was changed in the model file
  // NOTE: Not necessarily the time when the target was updated!

  // Header
  LibFront::output_string(out, IS, 0, "Timestamps");

  LibFront::output_scalar(out, IS, 1, "Front", NULL, mi->frontTs, true);
  LibFront::output_scalar(out, IS, 1, "Database", NULL, mi->databaseTs, true);
  LibFront::output_scalar(out, IS, 1, "Grid Parameter", NULL, mi->meshParameterTs, true);
  LibFront::output_scalar(out, IS, 1, "Mesh", NULL, mi->meshTs, true);
  LibFront::output_scalar(out, IS, 1, "Solver", NULL, mi->solverTs, true);
  LibFront::output_scalar(out, IS, 1, "GebhardtFactors", NULL, mi->gebhardtFactorsTs, true);
  LibFront::output_scalar(out, IS, 1, "ViewFactors", NULL, mi->viewfactorsTs, true);

  emf_outputSectionEnd(out);

  return out;
}


// Output model's vertex-table
//
ostream&
ModelOutputManager::emf_outputVertexTable(ostream& out)
{
  VertexTable* vt = (VertexTable*)model->getVertexTable();

  if ( vt == NULL || vt->dim1 * vt->dim2 == 0 ) return out;

  emf_outputSectionStart(out);

  LibFront::output_string(out, ESF_INDENT_SIZE, 0, EMF_VERTEX_TABLE, true);
  LibFront::output_string(out, ESF_INDENT_SIZE, 1, EMF_POINTS, true);
  LibFront::output_string(out, ESF_INDENT_SIZE, 2, EMF_SIZE, false);
  out << " " << vt->dim1 << " " << vt->dim2 << endl;
  LibFront::output_string(out, ESF_INDENT_SIZE, 3, EMF_REAL, true);

  const char* def = getMatcString(vt->matcTable, EMF_POINTS);

  //--Points defined by a Matc-function/variable
  if ( model->keepMatcDefinitions() && def != NULL ) {
    LibFront::output_matcDef(out, ESF_INDENT_SIZE, 3, NULL, NULL, def);

  //--List of points
  } else {

    for ( int i = 0; i < vt->dim1; i++) {
      BodyElement* v = model->getBodyElementById(vt->vertexIds[i]);
      GcPoint* p = (GcPoint*)v->getGeometry();
      LibFront::output_vector(out, ESF_INDENT_SIZE, 3, NULL, NULL, vt->dim2, (double*)p->getPoint(), false);
    }
  }

  emf_outputSectionEnd(out);

  return out;
}


void
ModelOutputManager::initClass(Model* mdl)
{
  ModelOutputManager::model = mdl;
}



//******************************************
//*** MESH INPUT format output functions ***
//******************************************

// Form Mesh Input File (mif-file)
//
ostream&
ModelOutputManager::mif_output(ostream& out)
{
  const ModelInfo* mi = model->getModelInfo();
  const ModelStatistics* ms = model->getModelStatistics();

  int i;

  //--Set tags for possible linearizing boundary points
  //  NOTE: also mark all points inactive for mesh
  //  Point to be included are marked active below via bodies!
  //--Also set possible boundary point mesh density values
  //
  model->setBoundaryPointData();

  int mesh_index = mi->currentMeshIndex;

  //--Find active mesh objects
  bool* active_objects = new bool[1 + ms->nofModelObjects];

  // Init flags
  for (i = 0; i <= ms->nofModelObjects; i++) {
    active_objects[i] = false;
  }

  int index = 0;

  // Mark active object via bodies!
  while (true) {
    Body* body = model->getBody(index++);
    if ( body == NULL ) break;
    if ( body->isExcludedFromMesh(mesh_index) ) continue;
    body->markActiveMeshObjects(mesh_index, active_objects);
  }

  //--Count active objects by type
  int nof_vertices = 0;
  int nof_edges = 0;
  int nof_bodies = 0;
  int nof_inner_boundaries = 0;
  int nof_outer_boundaries = 0;

  for (i = 0; i <= ms->nofModelObjects; i++) {

    if ( !active_objects[i] )
      continue;

    Body* body = NULL;
    BodyElement* be = NULL;

    ModelObject* obj = model->getModelObjectById(i);
    int bd1_id, bd2_id;

    if ( obj == NULL )
      continue;

    switch ( obj->getObjectType() ) {
    case OT_VERTEX:
      nof_vertices++;
      break;
    case OT_EDGE:
      be = model->getBodyElementById(i);
      nof_edges += be->getNofMifGeometries();
      break;
    case OT_BODY:
      body = model->getBodyById(i);
      if ( !body->isVirtual() && !body->isExcludedFromMesh(mesh_index) ) {
        nof_bodies++;
      }
      break;
    }
  }

  int nof_bpoints = mif_getNofBoundaryPoints();

  //--Output file
  // Comments
  // ========
  char tb[128];
  getCurrentTs(tb, 128);
  out << "!ElmerMesh input file from ElmerFront" << endl;
  out << "!Saved        = " << tb << endl;
  out << "!Case         = " << mi->modelName << "  " << mi->problemName << endl;
  out << "!Model dir    = " << mi->modelDirectory_absolute << endl;
  out << "!Include path = " << mi->includePath_absolute << endl;
  out << "!Results dir  = " << mi->resultsDirectory_absolute << endl;


  double mesh_h = mi->meshHs[mesh_index];
  double mesh_f = mi->meshFs[mesh_index];

  // Header
  // ======
  out << "Geometry2D:" << endl;
  if ( mesh_h > 0 )
    indent(out, 2) << "H: " << mesh_h << endl;

  if ( mesh_f > 0 )
    indent(out, 2) << "MeshScalingFactor: " << mesh_f << endl;

  indent(out, 2) << "Nodes: "  << nof_vertices + nof_bpoints << endl;
  indent(out, 2) << "Edges: "  << nof_edges << endl;
  indent(out, 2) << "Bodies: " << nof_bodies << endl;


  // Vertices + BoundaryPoints
  // =========================
  //out << "!Vertices: id, x, y, meshH" << endl;

  //-Vertices
  index = 0;
  while (true) {
    BodyElement* v = model->getVertex(index++);
    if (v==NULL) break;
    if ( active_objects[v->Id()] ) {
      v->output_mif(out);
    }
  }

  //-BoundaryPoints
  mif_outputBoundaryPoints(out);

  // Edges
  // =====
  //out << "!Edges: id, meshH, meshN, vertex count, vertex ids" << endl;
  index = 0;
  int next_mtag = 1;
  while (true) {
    BodyElement* e = model->getEdge(index++);
    if (e==NULL) break;
    if ( active_objects[e->Id()] ) {
      e->setMifTag(next_mtag);
      e->output_mif(out);
    }
  }

  // Bodies
  // ======
  index = 0;
  while (true) {
    Body* b = model->getBody(index++);
    if ( b == NULL ) break;
    if ( b->isExcludedFromMesh(mesh_index) ) continue;
    if ( !b->isVirtual() && active_objects[b->Id()] ) {
      b->output_mif(out);
    }
  }

  // End of file
  // ===========
  out << "!End" << endl;

  return out;
}


// Get nof active mesh boundary points
//
int
ModelOutputManager::mif_getNofBoundaryPoints()
{
  int i, index;
  int active_count = 0;

  int bp_count;
  BoundaryPoint** bp_points;

  //--First mark all unchecked (to handle copied points)
  //
  index = 0;
  while (true) {
    BodyElement* be = model->getBoundary(index++);
    if (be==NULL) break;

    bp_count = 0;
    bp_points = NULL;

    be->getBoundaryPoints(bp_count, bp_points);

    for (i = 0; i < bp_count; i++) {

      BoundaryPoint* bp = bp_points[i];

      // If this is an original vertex, we do not mark it
      if ( bp->isVertex() ) {
        continue;
      }

      bp->checked = false;
    }

    delete[] bp_points;
    bp_points = NULL;
  }

  //--Next count nof actives in meshing
  //
  index = 0;
  while (true) {
    BodyElement* be = model->getBoundary(index++);
    if (be==NULL) break;

    bp_count = 0;
    bp_points = NULL;

    be->getBoundaryPoints(bp_count, bp_points);

    for (i = 0; i < bp_count; i++) {

      BoundaryPoint* bp = bp_points[i];

      if ( !bp->activeInMeshing ||
           bp->isVertex()       ||
           bp->checked
         ) {
        continue;
      }

      active_count++;

      // Mark as handled
      bp->checked = true;
    }

    delete[] bp_points;
    bp_points = NULL;

  }

  return active_count;
}


// Output linearizing boundary points (as nodes) to
// mesh input file (mif-file)
// NOTE: original vertices are not output here, they
// are output separately (before these points)
//
ostream&
ModelOutputManager::mif_outputBoundaryPoints(ostream& out)
{
  int i, index;

  Point3 p;
  bool only_active = true;

  // For mesh density info
  int mesh_index = model->getCurrentMeshIndex();

  int bp_count;
  BoundaryPoint** bp_points;

  //--First mark all unchecked to handle copied (ex. loop end) points
  //
  index = 0;
  while (true) {
    BodyElement* be = model->getBoundary(index++);
    if (be==NULL) break;

    bp_count = 0;
    bp_points = NULL;

    be->getBoundaryPoints(bp_count, bp_points);

    for (i = 0; i < bp_count; i++) {

      BoundaryPoint* bp = bp_points[i];

      // If this is an original vertex, we do not mark
      if ( bp->isVertex() ) {
        continue;
      }

      bp->checked = false;
    }

    delete[] bp_points;
    bp_points = NULL;
  }

  //---Next output all active
  //
  index = 0;
  while (true) {
    BodyElement* be = model->getBoundary(index++, only_active);
    if (be==NULL) break;

    bp_count = 0;
    bp_points = NULL;

    be->getBoundaryPoints(bp_count, bp_points);

    for (i= 0; i < bp_count; i++) {

      BoundaryPoint* bp = bp_points[i];

      if ( !bp->activeInMeshing ||
           bp->isVertex()       ||
           bp->checked
         ) {
        continue;
      }

      bp->point->getPoint(p);

      bp->checked = true;

      out << "NodeId: " << bp->tag << " -1";
      if ( bp->meshDensityType != ' ' ) {
        out << "  " << bp->meshDensityType << ": " << bp->meshDensityValue;
      }

      // NOTE: Only two coordinates are output for Mesh2D
      //
      out << "  " << p[0] << " " << p[1];
      out << endl;
    }

    delete[] bp_points;
    bp_points = NULL;
  }

  return out;
}


//*******************************
//*** MATC definitions output ***
//*******************************

ostream&
ModelOutputManager::outputMatcDefinitions(ostream& out, int nof_defs, char** defs, bool dsign)
{
  int i;

  // First write variables
  for (i = 0; i < nof_defs; i++) {
    if ( 0 == strncmp("function", LibFront::trimLeft(defs[i]), 8) ) continue;
    if ( dsign ) out << '$';
    out << defs[i] << endl;
  }

  // Then write functions
  for (i = 0; i < nof_defs; i++) {
    if ( 0 != strncmp("function", LibFront::trimLeft(defs[i]), 8) ) continue;
    if ( dsign ) out << '$';
    out << defs[i] << endl;
  }

  return out;
}


//********************************************
//*** SOLVER INPUT format output functions ***
//********************************************

// Form Solver Input File (sif-file)
//
ostream&
ModelOutputManager::sif_output(ostream& out)
{
  UserInterface* gui = theControlCenter->getGui();

  Parameter* p;
  ParameterField* pf;
  int index;

  const ModelInfo* mi = model->getModelInfo();
  const ModelStatistics* ms = model->getModelStatistics();

  //-Calculate nof active solvers
  int nof_solvers = 0;
  index = 0;
  while (true) {
    p = model->getParameter(index++, ECIF_SOLVER);
    if (p==NULL) break;
    pf = p->getFieldBySifName("Active");
    if ( p->IsActive() ) {
        nof_solvers++;
    }
  }

  //-Calculate nof active calculators
  int nof_calculators = 0;
  index = 0;
  while (p != NULL) {
    p = model->getParameter(index++, ECIF_CALCULATOR);
    if (p==NULL) break;
    pf = p->getFieldBySifName("Active");
    if ( p->IsActive() ) {
        nof_calculators++;
    }
  }

  // Info as comments
  out << "!ElmerSolver input file from ElmerFront" << endl;
  out << "!Saved        = " << mi->modified << endl;
  out << "!Case         = " << mi->modelName << "  " << mi->problemName << endl;
  out << "!Model dir    = " << mi->modelDirectory_absolute << endl;
  out << "!Include path = " << mi->includePath_absolute << endl;
  out << "!Results dir  = " << mi->resultsDirectory_absolute << endl;
  out << endl;

  // NOTE: output "max counters" for Solver, ie highest id in use!
  // NOT: These are not used any more by solver (Elmer3.0-->), so printed as comments
  //
  out << "!Bodies "              << ms->nofBodies << endl;
  out << "!Equations "           << ms->nofEquations << endl;
  out << "!Solvers "             << nof_solvers + nof_calculators << endl;
  out << "!Materials "           << ms->nofMaterials << endl;
  out << "!Body Forces "         << ms->nofBodyForces << endl;
  out << "!Initial Conditions "  << ms->nofInitialConditions << endl;
  out << "!Boundary Conditions " << ms->nofBoundaryConditions << endl;
  out << "!Boundaries "          << ms->nofInnerBoundaries + ms->nofOuterBoundaries << endl;
  out << endl;

  // Sif file echo on/off
  p = model->getParameterById(ECIF_SOLVER_CONTROL, 1);
  if ( p != NULL ) {
    bool echo_on = false;
    p->getFieldValueBySifName(SIF_ECHO_ON, echo_on);

    // If not on, make it a comment!
    if ( !echo_on ) {
      out << "!";
    }
    out << "echo on" << endl;
    out << endl;
  }

  sif_outputHeader(out);

  // Possible matc input file
  //
  if ( mi->matcInputFile_sif != NULL &&
       mi->matcInputFile_sif[0] != '\0'
     ) {
    out << "!Matc input file" << endl;
    out << "$source(\"" << mi->matcInputFile_sif << "\")" << endl;
    out << endl;
  }

  // Possible matc definitions
  //
  int nof_defs = 0;
  char** defs = NULL;
  gui->getMatcSifDefinitions(nof_defs, defs);

  // Write matc definitions into Sif-file (true <--> with $-sign)
  if ( nof_defs > 0 ) {
    out << "!Matc definitions" << endl;
    outputMatcDefinitions(out, nof_defs, defs, true);
    out << endl;
  }

  // Possible model parameter
  //
  p = model->getParameterById(ECIF_MODEL_PARAMETER, 1);
  if ( p != NULL ) {

    SifOutputControl soc;
    soc.outputId = false;
    soc.outputName = false;
    soc.outputType = false;
    soc.outputAll = false;
    soc.sectionName = SIF_HEADER;

    p->output_sif(out, ESF_INDENT_SIZE, 0, soc);
    out << endl;
  }

  sif_outputSimulation(out);
  sif_outputConstants(out);
  sif_outputBodies(out);
  sif_outputEquations(out);
  sif_outputSolvers(out);
  sif_outputMaterials(out);
  sif_outputBodyForces(out);
  sif_outputInitialConditions(out);
  sif_outputBoundaryConditions(out);
  sif_outputBoundaries(out);
  sif_outputEof(out);

  return out;
}



//Data after block header
ostream&
ModelOutputManager::sif_outputAfterHeader(ostream& out)
{
  out << endl << "!" << endl;
  return out;
}


// Output bodies.
ostream&
ModelOutputManager::sif_outputBodies(ostream& out)
{
  int index = 0;
  while (true) {
    Body* body = model->getBody(index++);
    if (body==NULL) break;
    sif_outputSectionStart(out);
    body->output_sif(out, ESF_INDENT_SIZE, 0);
    sif_outputSectionEnd(out);
  }
  return out;
}



//----- ESF Output functions (Elemer Solver input format)

// Output body-forces data.
ostream&
ModelOutputManager::sif_outputBodyForces(ostream& out)
{

  SifOutputControl soc;
  soc.sectionName = SIF_BODY_FORCE;

  int index = 0;
  while (true) {
    Parameter* force = model->getParameter(index++, ECIF_BODY_FORCE);
    if (force==NULL) break;
    if ( force->getApplyCount() > 0 ) {
      sif_outputSectionStart(out);
      force->output_sif(out, ESF_INDENT_SIZE, 0, soc);
      sif_outputSectionEnd(out);
    }
  }

  return out;
}


// output All Boundaries.
ostream&
ModelOutputManager::sif_outputBoundaries(ostream& out)
{

  //BodyElement* be = model->getFirstBodyElement();
  int index = 0;

  while (true) {
    BodyElement* be = model->getBoundary(index++);

    if (be==NULL) break;

    //-Check status
    beStatus be_stat = be->getStatus();

    if ( be_stat & (BE_DEVIDED | BE_SWAPPED) ) {
      continue;
    }

    // NOTE: Only those boundaries are output which have
    // a boundary parameter defined (this is the only way
    // to output boundary parameters!)
    //
    if ( NO_INDEX == be->getBoundaryParameterId() ) {
      continue;
    }

    sif_outputSectionStart(out);

    be->output_sif(out, ESF_INDENT_SIZE, 0);

    sif_outputSectionEnd(out);
  }

  return out;
}


#if 0
// output All Boundaries.
ostream&
ModelOutputManager::sif_outputBoundaries(ostream& out)
{
  int count, i;
  int boundary_condition_id;

  // Outer boundaries.
  OuterBoundary* ob;
  count = modelData->outerBoundaries->size();
  for (i = 0; i < count; i++) {
    ob = (*modelData->outerBoundaries)[i];
    if (ob == 0)
      continue;
    boundary_condition_id = ob->outerBoundary->getConditionId();

    sif_outputSectionStart(out);
    out << "Boundary " << ob->outerBoundary->Tag()
        << endl;

    out << LibFront::indent(ESF_INDENT_SIZE, 1) << "Body 1"
        << endl;
    out << LibFront::indent(ESF_INDENT_SIZE, 2) << "Integer "
        << ob->parentBody->Tag()
        << endl;

    // NOTE: bc-cond-id -1 is also output!
    out << endl;
    out << LibFront::indent(ESF_INDENT_SIZE, 1) << "Boundary Condition"
        << endl;
    out << LibFront::indent(ESF_INDENT_SIZE, 2) << "Integer "
        << boundary_condition_id
        << endl;

    sif_outputSectionEnd(out);
  }

  // Inner boundaries
  InnerBoundary* ib;
  count = modelData->innerBoundaries->size();
  for (i = 0; i < count; i++) {
    ib = (*modelData->innerBoundaries)[i];
    if (ib == 0)
      continue;
    boundary_condition_id = ib->innerBoundary->getConditionId();

    sif_outputSectionStart(out);
    out << "Boundary " << ib->innerBoundary->Tag()
        << endl;
    out << LibFront::indent(ESF_INDENT_SIZE, 1) << "Body 1"
        << endl;

    out << LibFront::indent(ESF_INDENT_SIZE, 2) << "Integer "
        << ib->parentPair[0]->body->Tag()
        << endl;

    out << endl;

    out << LibFront::indent(ESF_INDENT_SIZE, 1) << "Body 2"
        << endl;
    out << LibFront::indent(ESF_INDENT_SIZE, 2) << "Integer "
        << ib->parentPair[1]->body->Tag()
        << endl;

    if (boundary_condition_id != NO_INDEX) {
      out << endl;
      out << LibFront::indent(ESF_INDENT_SIZE, 1) << "Boundary Condition"
          << endl;
      out << LibFront::indent(ESF_INDENT_SIZE, 2) << "Integer "
          << boundary_condition_id
          << endl;
    }
    sif_outputSectionEnd(out);
  }

  return out;
}
#endif


// Output all boundary conditions (outer, inner).
ostream&
ModelOutputManager::sif_outputBoundaryConditions(ostream& out)
{

  SifOutputControl soc;
  soc.sectionName = SIF_BOUNDARY_CONDITION;

  int index = 0;
  while (true) {
    Parameter* cond = model->getParameter(index++, ECIF_BOUNDARY_CONDITION);
    if (cond==NULL) break;
    if (cond->getApplyCount() > 0 ) {
      cond->updateTargetTags();
      sif_outputSectionStart(out);
      cond->output_sif(out, ESF_INDENT_SIZE, 0, soc);
      sif_outputSectionEnd(out);
    }
  }

  return out;
}


// Output physial constants.
ostream&
ModelOutputManager::sif_outputConstants(ostream& out)
{
  sif_outputSectionCode(out, SIF_CONSTANT);

  SifOutputControl soc;
  soc.outputType = false;
  soc.outputName = false;
  soc.sectionName = SIF_CONSTANT;

  int index =0;
  while (true) {
    Parameter* constant = model->getParameter(index++, ECIF_CONSTANT);
    if (constant==NULL) break;
    sif_outputSectionStart(out);
    constant->output_sif(out, ESF_INDENT_SIZE, 0, soc);
    sif_outputSectionEnd(out);
  }

  return out;
}


// Output all initial conditions.
ostream&
ModelOutputManager::sif_outputInitialConditions(ostream& out)
{

  SifOutputControl soc;
  soc.sectionName = SIF_INITIAL_CONDITION;

  int index = 0;
  while (true) {
    Parameter* cond = model->getParameter(index++, ECIF_INITIAL_CONDITION);
    if (cond==NULL) break;
    if (cond->getApplyCount() > 0 ) {
      sif_outputSectionStart(out);
      cond->output_sif(out, ESF_INDENT_SIZE, 0, soc);
      sif_outputSectionEnd(out);
    }
  }

  return out;
}


// Mark End Of File.
ostream&
ModelOutputManager::sif_outputEof(ostream& out)
{
  out << "!End Of File" << endl;
  return out;
}


// Output equations data.
ostream&
ModelOutputManager::sif_outputEquations(ostream& out)
{
  SifOutputControl soc;
  soc.sectionName = SIF_EQUATION;

  int index = 0;
  while (true) {
    Parameter* equation = model->getParameter(index++, ECIF_EQUATION);
    if (equation==NULL) break;
    if ( equation->getApplyCount() > 0 ) {
      sif_outputSectionStart(out);
      equation->output_sif(out, ESF_INDENT_SIZE, 0, soc);
      sif_outputSectionEnd(out);
    }
  }

  return out;
}


// Put model header data into ouput file.
ostream&
ModelOutputManager::sif_outputHeader(ostream& out)
{
  Parameter* p;
  ParameterField* pf;

  char key_buffer[1024];
  char value_buffer[1024];

  char QM = '\"';

  const ModelInfo* mi = model->getModelInfo();

  sif_outputSectionStart(out);
  sif_outputSectionCode(out, SIF_HEADER);


  // Keyword checking (ignore, warn, abort)
  p = model->getParameterById(ECIF_SOLVER_CONTROL, 1);
  if ( p != NULL ) {
    LibFront::toUpper(SIF_CHECK_KEYWORDS, key_buffer);
    p->getFieldValueBySifName(SIF_CHECK_KEYWORDS, 1023, value_buffer);
    indent(out, ESF_INDENT_SIZE);
    out << key_buffer << " " << value_buffer;
    out << endl;
  }

  // Meshdir
  LibFront::indent(out, ESF_INDENT_SIZE, 1)
    << "Mesh DB "
    << QM << model->getMeshDirValue() << QM  // "MESHDIR"
    << " ";

    if ( mi->nofActiveMeshes > 0 &&
         mi->meshNames != NULL
       ) {
      out << QM << mi->meshNames[mi->activeMeshIndices[0]] << QM;

    } else {
      out << QM <<  QM;
    }
    out << endl;

  LibFront::indent(out, ESF_INDENT_SIZE, 1) << "Include Path "      << QM << mi->includePath << QM << endl;
  LibFront::indent(out, ESF_INDENT_SIZE, 1) << "Results Directory " << QM << mi->resultsDirectory << QM << endl;

  sif_outputSectionEnd(out);

  return out;
}


// Output materials data.
ostream&
ModelOutputManager::sif_outputMaterials(ostream& out)
{

  SifOutputControl soc;
  soc.sectionName = SIF_MATERIAL;

  int index = 0;
  while (true) {
    Parameter* mater = model->getParameter(index++, ECIF_MATERIAL);
    if (mater==NULL) break;
    if ( mater->getApplyCount() ) {
      sif_outputSectionStart(out);
      mater->output_sif(out, ESF_INDENT_SIZE, 0, soc);
      sif_outputSectionEnd(out);
    }
  }

  return out;
}


// Output code for a new section (SEC_HEADER etc.)
ostream&
ModelOutputManager::sif_outputSectionCode(ostream& out, const char* section_cd)
{
  out << section_cd << endl;
  return out;
}


// Data ouput after a section
ostream&
ModelOutputManager::sif_outputSectionEnd(ostream& out)
{
  out << "End" << endl << endl;
  return out;
}


// Data ouput before a section
ostream&
ModelOutputManager::sif_outputSectionStart(ostream& out)
{
  //Currently nothing!
  return out;
}


// Put problem data into ouput file.
ostream&
ModelOutputManager::sif_outputSimulation(ostream& out)
{
  Parameter* p;

  const ModelInfo* mi = model->getModelInfo();

  sif_outputSectionStart(out);
  sif_outputSectionCode(out, SIF_SIMULATION);

  SifOutputControl soc;
  soc.outputType = false;
  soc.outputName = false;
  soc.sectionName = SIF_SIMULATION;

  // Solver control
  p = model->getParameterById(ECIF_SOLVER_CONTROL, 1);

  if ( p != NULL ) {
    SifOutputControl soc;
    soc.outputId = false;
    soc.outputName = false;
    soc.outputType = false;
    soc.outputAll = false;
    soc.sectionName = SIF_SIMULATION;

    p->output_sif(out, ESF_INDENT_SIZE, 0, soc);
  }

  // Simulation parameter
  p = model->getParameterById(ECIF_SIMULATION_PARAMETER, 1);

  if ( p != NULL ) {

    SifOutputControl soc;
    soc.outputId = false;
    soc.outputName = false;
    soc.outputType = false;
    soc.outputAll = false;
    soc.sectionName = SIF_SIMULATION;

    p->output_sif(out, ESF_INDENT_SIZE, 0, soc);
    out << endl;
  }

  // Coordinate stuff
  Parameter* coordinate = model->getParameter(0, ECIF_COORDINATE);
  if (coordinate != NULL ) {
    sif_outputSectionStart(out);
    coordinate->output_sif(out, ESF_INDENT_SIZE, 0, soc);
    out << endl;
  }

  // Timestep stuff (first active!)
  int index = 0;
  while (true) {
    Parameter* timestep = model->getParameter(index++, ECIF_TIMESTEP);
    if (timestep==NULL) break;
    ParameterField* pf = timestep->getFieldBySifName("Active");

    // If active, output and stop looping
    if ( timestep->IsActive()) {
      sif_outputSectionStart(out);
      timestep->output_sif(out, ESF_INDENT_SIZE, 1, soc);
      break;
    }
  }
  out << endl;

  // Result and inputfile stuff
  Parameter* datafile = model->getParameter(0, ECIF_DATAFILE);
  if ( datafile != NULL ) {
    sif_outputSectionStart(out);
    datafile->output_sif(out, ESF_INDENT_SIZE, 0, soc);
  }

  // Mesh Input File name
  char* mif_file_name = NULL;
  model->getMeshInputFileName(mif_file_name, 0);

  if ( mif_file_name != NULL ) {
    LibFront::output_scalar(out, ESF_INDENT_SIZE, 1, "Mesh Input File", "File", mif_file_name, true);
    out << endl;

    delete[] mif_file_name;
  }

  sif_outputSectionEnd(out);

  return out;
}


// Output solvers data.
ostream&
ModelOutputManager::sif_outputSolvers(ostream& out)
{
  int i;
  int index;
  Parameter* solver;
  Parameter* calculator;

  const ModelStatistics* ms = model->getModelStatistics();

  int count = ms->nofSolvers + ms->nofCalculators;

  Parameter** solvers = new Parameter*[count];

  for (i = 0; i < count; i++) {
    solvers[i] = NULL;
  }

  // Set active solvers and calculators according to the SOLVING_ORDER field
  // into the solvers table for output
  int counter = 0;

  //-Solvers
  index = 0;
  while (true) {
    solver = model->getParameter(index++, ECIF_SOLVER);
    if (solver==NULL) break;
    ParameterField* order_fld = solver->getFieldBySifName("Solving Order");

    if ( solver->IsActive() ) {
      counter++;
      int index;
      if (order_fld != NULL)
        index = atol(order_fld->getDataStrings()[0]);
      else
        index = counter;
      solvers[index - 1] = solver;
    }
  }

  //-Calculators
  index = 0;
  while (true) {
    calculator = model->getParameter(index++, ECIF_CALCULATOR);
    if (calculator==NULL) break;
    if ( calculator->IsActive() ) {
      counter++;

      int index;
      ParameterField* order_fld = calculator->getFieldBySifName("Solving Order");

      if (order_fld != NULL)
        index = atol(order_fld->getDataStrings()[0]);
      else
        index = counter;

      solvers[index - 1] = calculator;
    }
  }

  SifOutputControl soc;
  soc.outputName = false;
  soc.sectionName = SIF_SOLVER;

  Parameter* sc = model->getParameter(0, ECIF_SOLVER_CONTROL);
  ParameterField* pf = NULL;
  if ( sc != NULL ) {
    pf = sc->getFieldByGuiName("RELOAD_INPUT_FILE", true);
  }

  // If Reload-file field active, output reload info as a Solver
  //
  if ( pf != NULL ) {

    soc.reloadSolverIsOutput = true;

    LibFront::output_string(out, ESF_INDENT_SIZE, 0, "Solver 1", true);

    pf = sc->getFieldByGuiName("EXEC_SOLVER", true);
    if ( pf != NULL ) {
      pf->output_sif(out, ESF_INDENT_SIZE, 1, "Solver", true);
    }

    // Equation for the 'solver'
    if ( model->getSolverKeywordTypeGiven("Solver", "Equation" ) )
      LibFront::output_string(out, ESF_INDENT_SIZE, 1, "Equation =  \"ReloadInput\"", true);
    else
      LibFront::output_string(out, ESF_INDENT_SIZE, 1, "Equation =  String \"ReloadInput\"", true);

    // Procedure for the 'solver'
    if ( model->getSolverKeywordTypeGiven("Solver", "Procedure" ) )
      LibFront::output_string(out, ESF_INDENT_SIZE, 1, "Procedure = \"ReloadInput\" \"ReloadInput\"", true);
    else
      LibFront::output_string(out, ESF_INDENT_SIZE, 1, "Procedure = File \"ReloadInput\" \"ReloadInput\"", true);

    LibFront::output_string(out, ESF_INDENT_SIZE, 0, "End", true);
    out << endl;
  }

  // Output active solvers in solving order
  for (i = 0; i < count; i++) {
    if ( solvers[i] == NULL )
      continue;
    sif_outputSectionStart(out);
    solvers[i]->output_sif(out, ESF_INDENT_SIZE, 0, soc);
    sif_outputSectionEnd(out);
  }

  delete[] solvers;

  return out;
}


// Output those equation fields in whose target is actually a solver section
//
ostream&
ModelOutputManager::sif_outputSolverTargetFields(ostream& out, short indent_size, short indent_level, const char* source_eq_name)
{
  NameSet targetFieldNames;

  int nof_values;
  bool* values = NULL;

  int index = 0;
  while (true) {

    delete[] values;
    values = NULL;

    Parameter* equation = model->getParameter(index++, ECIF_EQUATION);
    if (equation==NULL) break;
    if ( equation->getApplyCount() == 0 ) continue;

    if ( !equation->getFieldValueBySifName(source_eq_name, nof_values, values) ) {
      continue;
    }

    if ( !values[0] ) continue;

    equation->outputSolverTargetFields_sif(out, indent_size, indent_level, source_eq_name, targetFieldNames);
  }

  return out;
}



//***********************************
//*** MESH FILES output functions ***
//***********************************

// Elmer format file
// =================

void
ModelOutputManager::write_Elmer_mesh(char* mesh_dir)
{
  const ModelInfo* mi = model->getModelInfo();
  const MeshInfo* mei = model->getMeshInfo();
  const ModelStatistics* ms = model->getModelStatistics();


  UserInterface* gui = theControlCenter->getGui();

  int info;
  int nodeC = mei->nofNodes;
  int elemC_net = mei->nofBulkElements - mei->nofSplittedElements; // To header!
  int belemC = mei->nofBoundaryElements;

  int elemC = mei->nofBulkElements; // For looping bulk table!
  int counter;
  int i, index;

  // These are 101-type boundary elements!
  mei->nofUsedElementTypes[0] = ms->nofVertices;

  // Count nof element types used
  int typeC = 0;
  for (i = 0; i < MAX_NOF_ELEM_CODES; i++) {
    if ( mei->nofUsedElementTypes[i] == 0 )
      continue;
    typeC++;
  }

  int* types = new int[typeC];
  int* typesC = new int[typeC];

  // Count nof elements by used types
  counter = 0;
  for (i = 0; i < MAX_NOF_ELEM_CODES; i++) {

    if ( mei->nofUsedElementTypes[i] == 0 )
      continue;

    types[counter] = MeshElementDesc[i][DESC_ELEM_TYPE];
    typesC[counter] = mei->nofUsedElementTypes[i];

    // Update boundary element counter
    if ( mi->dimension == ECIF_3D ) {

      // Add possible edges and vertices
      if ( types[counter] < 300 ) {
        belemC += typesC[counter];
      }

    } else {

      // Add possible vertices
      if ( types[counter] < 200 ) {
        belemC += typesC[counter];
      }
    }

    counter++;
  }


  // Init DB
  eio_init(info);

  // EIO (create mesh files)
  // -----------------------
  eio_create_mesh(mesh_dir, info);

  // Ouch, no success with the mesh!
  if (info == -1) {
    gui->errMsg(0, "Cannot create mesh in database using directory and name:", mesh_dir);
    return;
  }

  // EIO (mesh.header)
  // -----------------
  eio_set_mesh_description(nodeC, elemC_net, belemC, typeC, types, typesC, info);

  // Nodes
  // =====
  Point3* nodes = model->getMeshNodeData();

  double coordinates[3];
  for (i = 0; i < nodeC; i++) {
    for (short j = 0; j < 3; j++) {
      coordinates[j] = nodes[i][j];
    }

    // Elmer DB ids
    int node_id = i + 1;
    int cond_id = 0;

    // EIO (mesh.nodes)
    // ----------------
    eio_set_mesh_node(node_id, cond_id, coordinates, info);
  }

  int node_ids[32];

  // Elements
  // ========

  // Bulk elements
  // -------------
  MeshElementTable* bt = model->getMeshBulkElements();
  counter = 0;

  for (i = 0; i < elemC; i++) {

    if (bt->splitted[i]) {
      continue;
    }

    int elem_code = bt->getElementCode(i);
    int elem_type = MeshElementDesc[elem_code][DESC_ELEM_TYPE];
    int body_id = model->getModelObjectTagById(bt->parentIds[i][0]);
    int nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];
    const int* elem_node_ids = bt->getNodeIds(i);

    // Elmer DB ids
    int elem_ext_id = 1 + counter++;

    for (short j = 0; j < nof_nodes; j++)  {
      node_ids[j] = 1 + elem_node_ids[j];
    }

    // EIO (mesh.elements)
    // -------------------
    eio_set_mesh_element_conns(elem_ext_id, body_id, elem_type, node_ids, info);
  }

  // Boundary elements
  // -----------------
  counter = 0;
  BodyElement* be;

  //-Boundary elements proper (faces or edges)
  MeshElementTable* bet = model->getMeshBoundaryElements();

  index = 0;
  while (true) {
    be = model->getBoundary(index++);

    if (be==NULL) break;

    int nof_elements = be->getNofMeshElements();

    for (i = 0; i < nof_elements; i++) {

      int index = be->getMeshElementId(i);

      // parent element ids, Elmer DB numbering
      int p_elem1_int_id = bet->parentIds[index][0]; // Int parent elem1
      int p_elem2_int_id = bet->parentIds[index][1]; // Int parent elem2

      if ( model->modelHasBulkRenumbering() ) {
        p_elem1_int_id = model->getRenumberedMeshBulkElementId(p_elem1_int_id);
        p_elem2_int_id = model->getRenumberedMeshBulkElementId(p_elem2_int_id);
      }

      int p_elem1_ext_id = 1 + p_elem1_int_id;  // Ext parent elem1 id
      int p_elem2_ext_id = 1 + p_elem2_int_id;  // Ext parent elem2 id

      // Element code for the boundary element
      meshElementCode elem_code = bet->getElementCode(index);
      int elem_type = MeshElementDesc[elem_code][DESC_ELEM_TYPE];
      int nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];

      // Boundary tag for the element
      int bndr_tag = be->getBoundaryTag(); // Internal boundary tag
      const int* elem_node_ids = bet->getNodeIds(index);

      // Elmer DB ids (DB index origo is 1)
      for (int k = 0; k < nof_nodes; k++) {
        node_ids[k] = 1 + (int)elem_node_ids[k];
      }
      //int bndr_elem_id = index + 1;
      int elem_id = ++counter;

      info = 0;

      // EIO (mesh.boundary)
      // -------------------
      eio_set_mesh_bndry_element(elem_id, bndr_tag,
                                 p_elem1_ext_id, p_elem2_ext_id,
                                 elem_type, node_ids,
                                 info);

    } // for all mesh elements in the boundary

  } // All boundaries


  //-Boundary elements edges (3D models)
  //
  MeshElementTable* eet = model->getMeshBoundaryElementEdges();

  index = 0;
  while ( mi->dimension == ECIF_3D) {

    be = model->getEdge(index++);

    if (be==NULL) break;

    int nof_elements = be->getNofMeshElements();

    for (i = 0; i < nof_elements; i++) {

      int index = be->getMeshElementId(i);

      // No parents for edges in 3D!
      int p_elem1_ext_id = -1;
      int p_elem2_ext_id = -1;

       // Element code for the boundary element
      meshElementCode elem_code = eet->getElementCode(index);
      int elem_type = MeshElementDesc[elem_code][DESC_ELEM_TYPE];
      int nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];

      // Boundary tag for the element
      int bndr_tag = be->getBoundaryTag(); // Internal boundary tag
      const int* elem_node_ids = eet->getNodeIds(index);

      // Elmer DB ids (DB index origo is 1)
      for (int k = 0; k < nof_nodes; k++) {
        node_ids[k] = 1 + elem_node_ids[k];
      }
      //int bndr_elem_id = index + 1;
      int elem_id = ++counter;

      info = 0;

      // EIO (mesh.boundary)
      // -------------------
      eio_set_mesh_bndry_element(elem_id, bndr_tag,
                                 p_elem1_ext_id, p_elem2_ext_id,
                                 elem_type, node_ids,
                                 info);

    } // for all meshelements in the edge

  } // All 3D edges


  // Vertices
  // =================
  index = 0;
  while (true) {

    be = model->getVertex(index++);

    if (be==NULL) break;

    int elem_id = ++counter;
    int bndr_tag = be->getBoundaryTag();
    int p_elem1_ext_id = -1;
    int p_elem2_ext_id = -1;
    int elem_type = 101;
    node_ids[0] = 1 + be->getMeshElementId(0);

    // EIO (mesh.boundary)
    // -------------------
    eio_set_mesh_bndry_element(elem_id, bndr_tag,
                               p_elem1_ext_id, p_elem2_ext_id,
                               elem_type, node_ids,
                               info);

  }

  // Close EIO
  eio_close_mesh(info);
  eio_close(info);

  delete[] types;
  delete[] typesC;

  gui->showMsg("Mesh saved in Elmer format!");

}


// Elmer Post format file
// ======================

void
ModelOutputManager::write_ElmerPost_mesh(ostream& outfile)
{
  int i, index;
  const ModelInfo* mi = model->getModelInfo();
  const MeshInfo* mei = model->getMeshInfo();
  const ModelStatistics* ms = model->getModelStatistics();

  UserInterface* gui = theControlCenter->getGui();

  int nofNodes = mei->nofNodes;
  int nofElements = mei->nofBulkElements - mei->nofSplittedElements;
  int nofBoundaryElements = mei->nofBoundaryElements;

  // Header
  // =====
  outfile << nofNodes << "  " << nofElements + nofBoundaryElements;
  outfile << endl;

  // Node coordinates
  // ================
  Point3* nodes = model->getMeshNodeData();

  for (i = 0; i < nofNodes; i++) {

    Point3& p = nodes[i];
    outfile << "  " << p[0] << " " << p[1] << " " << p[2] << endl;
  }

  // Bulk elements from bodies
  // =========================
  MeshElementTable* bt = model->getMeshBulkElements();

  outfile << "#group all-bodies" << endl;

  while (true) {

    Body* body = model->getBody(index++);

    if (body==NULL) break;

    int nof_elements = body->getNofMeshElements();

    if (nof_elements == 0) {
      continue;
    }

    //--Body name as the group name
    outfile << "#group " << body->getName() << endl;

    // All body mesh elements
    for (i = 0; i < nof_elements; i++) {

      int eid = body->getMeshElementId(i);

      const int* nodes = bt->getNodeIds(eid);
      meshElementCode elem_code = bt->getElementCode(eid);
      int body_id = bt->parentIds[eid][0];

      outfile << "all ";

      // Elm type for the bulk element
      outfile << MeshElementDesc[elem_code][DESC_ELEM_TYPE];

      // Internal node ids
      outfile << " ";
      int nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];
      for (int j = 0; j < nof_nodes; j++) {
        outfile << " " << nodes[j];
      }
      outfile << endl;
    }

    //--End for this body group
    outfile << "#endgroup " << body->getName() << endl;
  }

  //---End for all bodies
  outfile << "#endgroup all-bodies" << endl;


  // Boundary elements from boundaries
  // =================================
  MeshElementTable* bet = model->getMeshBoundaryElements();

  outfile << "#group all-boundaries" << endl;

  // Loop all boundaries
  index = 0;
  while (true) {
    BodyElement* be = model->getBoundary(index++);
    if (be==NULL) break;
    int nof_elements = be->getNofMeshElements();
    if (nof_elements == 0) {
      continue;
    }

    const beStatus be_status = be->getStatus();
    if ( be_status & (BE_DEVIDED | BE_SWAPPED) ) {
      continue;
    }

    //--Boundary name as the group name
    outfile << "#group " << be->getName() << endl;

    for (i = 0; i < nof_elements; i++) {

      int eid = be->getMeshElementId(i);

      const int* nodes = bet->getNodeIds(eid);
      int parent1_id = bet->parentIds[eid][0];
      int body1_id = bt->parentIds[parent1_id][0];
      meshElementCode elem_code = bet->getElementCode(eid);

      outfile << "all ";

      // Elm type for the boundary element
      outfile << MeshElementDesc[elem_code][DESC_ELEM_TYPE];

      // Internal node ids for the boundary element
      int nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];
      outfile << " ";
      for (int k = 0; k < nof_nodes; k++) {
        outfile << " " << nodes[k];
      }
      outfile << endl;
    }

    //--End for this boundary
    outfile << "#endgroup " << be->getName() << endl;
  }

  //---End for all boundaries
  outfile << "#endgroup all-boundaries" << endl;

  //---End for all elements
  outfile << "#endgroup all" << endl;


  // All done
  // ========
  gui->showMsg("Mesh saved in Elmer Post format!");

}



// Thetis format file
// ==================

void
ModelOutputManager::write_Thetis_mesh(ostream& outfile)
{
  const ModelInfo* mi = model->getModelInfo();
  const MeshInfo* mei = model->getMeshInfo();
  const ModelStatistics* ms = model->getModelStatistics();

  UserInterface* gui = theControlCenter->getGui();

  int i,j, index;
  int nofNodes = mei->nofNodes;
  int nofBulkElements = mei->nofBulkElements;
  int nofBulkElements_net = mei->nofBulkElements - mei->nofSplittedElements;
  int nofBodies = mei->nofBodies;
  int nofBoundaries = ms->nofInnerBoundaries + ms->nofOuterBoundaries;
  int nofBndrElements = mei->nofBoundaryElements;

  //---Header
  outfile << "# nof: Nodes, Max external-id, Elements, Max external-id, Bodies, Nof boundaries, Nof boundary elements\n";
  outfile <<         nofNodes    << " " << mei->maxExternalNodeId;
  outfile << "  " << nofBulkElements_net << " " << mei->maxExternalElementId;
  outfile << "  " << nofBodies;
  outfile << "  " << nofBoundaries;
  outfile << "  " << nofBndrElements;
  outfile << endl;

  //---Node coordinates
  outfile << "##" << endl;
  outfile << "#Nodes:" << endl;
  //outfile << "External-id, Coordinates, Nof bodies, External body ids" << endl;
  outfile << "# External-id, Coordinates" << endl;

  Point3* nodes = model->getMeshNodeData();

  for (i = 0; i < nofNodes; i++) {

    Point3& p = nodes[i];
    outfile << model->getMeshNodeIdInt2Ext(i); // External id
    outfile << "  " << p[0] << " " << p[1] << " " << p[2];
    // nof node bodies and external body_ids
    //short nof_node_bodies = meshNode2Bodies[i][0];
    //outfile << "  " << nof_node_bodies;
    //for (short j = 0; j < nof_node_bodies; j++) {
      //int body_id = meshNode2Bodies[i][1 + j];
      //outfile << " " << meshBodyInt2Ext[body_id];
    //}
    outfile << endl;
  }


  MeshElementTable* bt = model->getMeshBulkElements();
  MeshElementTable* bet = model->getMeshBoundaryElements();

  //---Elements
  outfile << "##" << endl;
  outfile << "#Elements:" << endl;
  outfile << "# External-id, External body-id, Element-code" << endl;
  outfile << "# Nof nodes, External node-ids" << endl;
  outfile << "# Nof neighbour elements , External neighbour ids" << endl;

  for (i = 0; i < nofBulkElements; i++) {

    if ( bt->splitted[i] ) {
      continue;
    }

    const int* node_ids = bt->getNodeIds(i);
    meshElementCode elem_code = bt->getElementCode(i);

    int body_id = bt->parentIds[i][0];
    int body_tag = model->getModelObjectTagById(body_id);

    outfile << model->getMeshBulkElementIdInt2Ext(i);
    outfile << " " << model->getBodyTagInt2Ext(body_tag);
    outfile << " " << MeshElementDesc[elem_code][DESC_ELEM_TYPE];

    // External node ids
    outfile << "   ";
    int nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];
    outfile << nof_nodes << " ";

    for (j = 0; j < nof_nodes; j++) {
      outfile << " " << model->getMeshNodeIdInt2Ext(node_ids[j]);
    }

    // External neigbour element ids
    outfile << "    ";
    int nof_bndr_elems = MeshElementDesc[elem_code][DESC_NOF_BNDR_ELEMS];
    outfile << nof_bndr_elems << " ";

    for (j = 0; j < nof_bndr_elems; j++) {
      int elem_id = bt->neighborIds[i][j];
      outfile << " " << model->getMeshBulkElementIdInt2Ext(elem_id);
    }
    outfile << endl;
  }

  //---Boundaries
  outfile << "##" << endl;
  outfile << "#Boundaries:" << endl;
  outfile << "# Boundary id, "
          << "External body1-id, body2-id, Nof bndr elements" << endl;

  index = 0;
  while (true) {

    BodyElement* be = model->getBoundary(index++);

    if (be==NULL) break;

    // Name as comment line
    outfile << "#" << be->getName() << endl;

    // Parent tags
    int body1_tag = be->getParentTag(1);
    int body2_tag = be->getParentTag(2);

    outfile << be->Tag();

    outfile << "  " << model->getBodyTagInt2Ext(body1_tag);

    if (body2_tag != NO_INDEX)
      outfile << "  " << model->getBodyTagInt2Ext(body2_tag);
    else
      outfile << "  " << NO_INDEX;

    // Nof elements
    outfile << "  " << be->getNofMeshElements();
    outfile << endl;
  }

  //---Boundary elements
  outfile << "##" << endl;
  outfile << "#Boundary elements:" << endl;
  outfile << "# Sequence nbr, Boundary id" << endl;
  outfile << "# External parent-elem1-id, parent-elem2-id" << endl;
  outfile << "# External body1-id, body2-id, Element-code" << endl;
  outfile << "# Nof nodes, External node-ids" << endl;

  // Loop all boundary elements
  index = 0;
  while (true) {

    BodyElement* be = model->getBoundary(index++);

    if (be==NULL) break;

    int nof_elements = be->getNofMeshElements();

    for (i = 0; i < nof_elements; i++) {

      int index = be->getMeshElementId(i);
      int bndr_tag = be->Tag();

      const int* node_ids = bet->getNodeIds(index);
      meshElementCode elem_code = bet->getElementCode(index);

      outfile << i + 1                                 // Sequence nbr
              << " "
              << bndr_tag;                             // Boundary tag

      int p1_id = bet->parentIds[index][0];            // Int parent elem1
      int p2_id = bet->parentIds[index][1];            // Int parent elem2

      // parent elements
      outfile << "  " << model->getMeshBulkElementIdInt2Ext(p1_id);    // -->ext parent elem1
      outfile << " "  << model->getMeshBulkElementIdInt2Ext(p2_id);    // -->ext parent elem2

      // parent bodies
      int b1_id = bt->parentIds[p1_id][0];  // Int parent body1
      int b1_tag = model->getModelObjectTagById(b1_id);

      outfile << "  " << model->getBodyTagInt2Ext(b1_tag);            // Ext parent body1

      if (p2_id != NO_INDEX) {

        int b2_id = bt->parentIds[p2_id][0];// Int parent body2
        int b2_tag = model->getModelObjectTagById(b2_id);

        outfile << " " << model->getBodyTagInt2Ext(b2_tag);           // Ext parent body2
      } else {
        outfile << " " << NO_INDEX;                         // No parent elem2
      }

      // Elm type for the boundary element
      outfile << "  " << MeshElementDesc[elem_code][DESC_ELEM_TYPE];

      // Nodes in the boundary element
      int nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];
      outfile << "   ";
      outfile << nof_nodes << " ";

      // External node ids for the boundary element
      for (int k = 0; k < nof_nodes; k++) {
        outfile << " " << model->getMeshNodeIdInt2Ext(node_ids[k]);
      }
      outfile << endl;
    }
  }

  gui->showMsg("Mesh saved in Thetis format!");
}

