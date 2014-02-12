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
Module:     ecif_model_aux.cpp
Language:   C++
Date:       20.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_bodyElement.h"
#include "ecif_geometry.h"
#include "ecif_modelObject.h"
#include "ecif_model_aux.h"
#include "ecif_parameter.h"
#include "ecif_mesh.h"


// MeshConnectionTable
// ===================

// Init static class variables
int MeshConnectionTable::allocSize = 0;
int MeshConnectionTable::reallocSize = 0;
int MeshConnectionTable::totalCount = 0;

MeshConnectionTable::MeshConnectionTable(int count, int alloc_size, int realloc_size,
                                         bool enumerate_connections)
{
  allocSize = alloc_size;
  reallocSize = realloc_size;
  totalCount = count;
  enumerateConnections = enumerate_connections;

  nofConnections = new int[totalCount];
  maxNofConnections = new int[totalCount];
  lastConnectionNbr = 0;

  parentIds = new int*[totalCount];

  if (enumerateConnections) {
    connectionNbrs = new int*[totalCount];
  } else {
    connectionNbrs = NULL;
  }

  for (int i = 0; i < totalCount; i++) {

    nofConnections[i] = 0;
    maxNofConnections[i] = allocSize;
    parentIds[i] = new int[allocSize];

    if (enumerateConnections) {
      connectionNbrs[i] = new int[allocSize];
    }

    for (int j = 0; j < allocSize; j++) {
      parentIds[i][j] = NO_INDEX;

      if (enumerateConnections) {
        connectionNbrs[i][j] = NO_INDEX;
      }
    }

  }
}


MeshConnectionTable::~MeshConnectionTable()
{
  for (int i = 0; i < totalCount; i++) {
    delete[] parentIds[i];

    if (enumerateConnections) {
      delete[] connectionNbrs[i];
    }
  }

  delete[] nofConnections;
  delete[] maxNofConnections;
  delete[] parentIds;
  delete[] connectionNbrs;
}


bool
MeshConnectionTable::addConnection(int index, int parent_id, bool only_unique)
{
  int connection_nbr;

  return addConnection(index, parent_id, connection_nbr, only_unique);
}

bool
MeshConnectionTable::addConnection(int index, int parent_id, int& connection_nbr,
                                   bool only_unique)
{
  connection_nbr = NO_INDEX;

  if ( index < 0 || index >= totalCount ) {
    return false;
  }

  int pos = nofConnections[index];

  // If only unique occurances and id already exist, do
  // not store it twice!
  if (only_unique) {

    for (int i = 0; i < pos; i++) {

      if ( parentIds[index][i] == parent_id ) {

        if (enumerateConnections) {
          connection_nbr = connectionNbrs[index][i];
        }

        return false;
      }
    }
  }

  // Allocate new space if needed
  if ( pos >= maxNofConnections[index] ) {

    int nsize = reallocSize + pos;

    int* tmp = new int[nsize];
    int* tmp2 = NULL;

    if (enumerateConnections) {
      tmp2 = new int[nsize];
    }

    for (int i = 0; i < nsize; i++) {

      if ( i < pos) {
        tmp[i] = parentIds[index][i];

        if (enumerateConnections) {
          tmp2[i] = connectionNbrs[index][i];
        }

      } else {
        tmp[i] = NO_INDEX;

        if (enumerateConnections) {
          tmp2[i] = NO_INDEX;
        }
      }
    }

    delete[] parentIds[index];
    parentIds[index] = tmp;

    if (enumerateConnections) {
      delete[] connectionNbrs[index];
      connectionNbrs[index] = tmp2;
    }

    maxNofConnections[index] = nsize;
  }

  // Add new index
  parentIds[index][pos] = parent_id;

  if (enumerateConnections) {
    connectionNbrs[index][pos] = ++lastConnectionNbr;
    connection_nbr = lastConnectionNbr;
  }

  nofConnections[index] = ++pos;

  return true;
}


bool
MeshConnectionTable::addConnections(int nof_indices, int* indices, int parent_id, bool only_unique)
{
  bool rc = false;
  int connection_nbr;

  for (int i = 0; i < nof_indices; i++) {

    if ( addConnection(indices[i], parent_id, connection_nbr, only_unique) ) {
      rc = true;
    }
  }

  return rc;
}


int
MeshConnectionTable::getConnectionNbr(int index, int parent_index)
{
  if (!enumerateConnections) {
    return NO_INDEX;
  }

  if ( index < 0 || index >= totalCount )
    return NO_INDEX;

  if ( parent_index < 0 || parent_index >= nofConnections[index] )
    return NO_INDEX;

  return connectionNbrs[index][parent_index];
}


const int*
MeshConnectionTable::getConnectionNbrs(int index)
{
  if ( !enumerateConnections || index < 0 || index >= totalCount )
    return NULL;

  return connectionNbrs[index];
}


int
MeshConnectionTable::getNofConnections(int index)
{
  if ( index < 0 || index >= totalCount )
    return 0;

  return nofConnections[index];
}


int
MeshConnectionTable::getParentId(int index, int parent_index)
{
  if ( index < 0 || index >= totalCount )
    return NO_INDEX;

  if ( parent_index < 0 || parent_index >= nofConnections[index] )
    return NO_INDEX;

  return parentIds[index][parent_index];
}


const int*
MeshConnectionTable::getParentIds(int index)
{
  if ( index < 0 || index >= totalCount )
    return NULL;

  return parentIds[index];
}


bool
MeshConnectionTable::getParentIds(int nof_indices, int* indices,
                               int max_buffer_size, int* parent_id_buffer,
                               int minimum_connection_count, int& nof_ids)
{
  static int slot_indices[64];
  static short slot_positions[64];
  int limit = 0x7fffffff; // 'Largest' as a comparison seed

  int i;
  bool rc = false;
  nof_ids = 0;

  //--Count nof active index slots
  int nof_slots = 0;
  for (i = 0; i < nof_indices; i++) {

    int index = indices[i];

    if ( index < 0 || index >= totalCount ) {
      continue;
    }

    slot_indices[i] = index;
    slot_positions[i] = 0;
    nof_slots++;
  }

  //--Pick those parent ids from the id-sets which meet the the connection
  //  count criterium
  while (1) {

    int smallest = limit;

    //-Compare head items in each id-set slot and select the smallest
    for (i = 0; i < nof_slots; i++) {

      int index = slot_indices[i];
      int pos = slot_positions[i];

      if ( pos < nofConnections[index] &&
           smallest > parentIds[index][pos]
         ) {
        smallest = parentIds[index][pos];
      }
    }

    //-If all slots are at the end, stop
    if (smallest == limit) {
      break;
    }

    //-Calc nof "smallest"
    int counter = 0;
    for (i = 0; i < nof_slots; i++) {

      int index = slot_indices[i];
      int pos = slot_positions[i];

      if ( pos < nofConnections[index] &&
           smallest == parentIds[index][pos]
         ) {
        counter++;
        slot_positions[i]++;
      }
    }

    //-Add it to the the result if the connection  count criterium is fulfilled
    if ( counter >= minimum_connection_count ) {

      if ( nof_ids < max_buffer_size ) {
        parent_id_buffer[nof_ids] = smallest;
      }

      nof_ids++;
    }

  }

  //--If something was found!
  if ( nof_ids > 0 ) {
    rc = true;
  } else {
    rc = false;
  }

  return rc;
}


// Sort ids in each parent-ids slot
// NOTE: Simple insert sort!
void
MeshConnectionTable::sortParentIds()
{
  for (int index = 0; index < totalCount; index++) {

    for (int i = 0; i < nofConnections[index] - 1; i++) {

      int min_pos = i;

      //-Find if there is smaller than i-th element
      for (int j = i + 1; j < nofConnections[index] - 1; j++) {

        if ( parentIds[index][j] < parentIds[index][min_pos] ) {
          min_pos = j;
        }
      }

      //-If smaller was found, swap i-th and smaller
      if (min_pos != i) {

        int tmp;

        tmp = parentIds[index][min_pos];
        parentIds[index][min_pos] = parentIds[index][i];
        parentIds[index][i] = tmp;

        if (enumerateConnections) {
          tmp = connectionNbrs[index][min_pos];
          connectionNbrs[index][min_pos] = connectionNbrs[index][i];
          connectionNbrs[index][i] = tmp;
        }
      }
    }

  } // All slots
}


// MeshInputElement
// =================
MeshInputElement::MeshInputElement()
{
  type = ' ';
  elementCode = MEC_000;
  elementId = NO_INDEX;
  extElementId = NO_INDEX;
  extParentTag = NO_INDEX;
  extNodeIds = NULL;
  parentTag = NO_INDEX;
  isAdded = '0';
}

MeshInputElement::~MeshInputElement()
{
  delete[] extNodeIds;
}


// MeshCornerElement
// =================
MeshCornerElement::MeshCornerElement()
{
  elementId = NO_INDEX;
  elementCode = 0;
  bodyId = NO_INDEX;
  nodeIds = NULL;
  hasZeroVelocity = false;
  corrected = false;
  splitted = false;
  nofSubElements = 0;
  subElementIds = NULL;
  boundaryElementIds = NULL;
  boundaryIds = NULL;
  newNodeId = NO_INDEX;
  nofNewElements = 0;
  newElementIds = NULL;
}


MeshCornerElement::~MeshCornerElement()
{
  delete[] nodeIds;
  delete[] subElementIds;
  delete[] boundaryElementIds;
  delete[] boundaryIds;
  delete[] newElementIds;
}


//************
// MESH DATA
// *********

MeshData::MeshData()
{
  nofNodes = 0;
  nofAllocNodes = 0;

  bulkEdgeRendered = NULL;
  bulkEdges = NULL;
  bulkElements = NULL;
  bulkRenumbering = NULL;

  boundaryElements = NULL;
  boundaryEdges = NULL;
  boundaryVertices = NULL;
  boundaryElementBoundaryIds = NULL;

  nodes = NULL;

  bulkElementExt2Int = NULL;
  bulkElementInt2Ext = NULL;
  nodeExt2Int = NULL;
  nodeInt2Ext = NULL;
  bodyExt2Int = NULL;
  bodyInt2Ext = NULL;

  nofBoundaryNodeParents = NULL;
  boundaryNodeParentIds = NULL;

  lastBoundaryElementIndex = 0;
  lastBulkElementIndex = 0;

  // Allocate body<-->matrial id tables
  bodyExt2Int = new int[1 + MAX_NOF_BODIES];
  bodyInt2Ext = new int[1 + MAX_NOF_BODIES];
  for (int i = 0; i < 1 + MAX_NOF_BODIES; i++) {
    bodyExt2Int[i] = NO_INDEX;
    bodyInt2Ext[i] = NO_INDEX;
  }
}


MeshData::~MeshData()
{
  // Elements related data
  delete bulkElements;
  delete bulkEdges;
  delete[] bulkEdgeRendered;
  delete[] bulkRenumbering;

  delete boundaryElements;
  delete boundaryEdges;
  delete boundaryVertices;
  delete[] boundaryElementBoundaryIds;

  delete[] nofBoundaryNodeParents;

  for (int i = 0; i < nofNodes; i++) {
    if ( boundaryNodeParentIds != NULL &&
         boundaryNodeParentIds[i] != NULL
       ) {
      delete[] boundaryNodeParentIds[i];
    }
  }

  delete[] boundaryNodeParentIds;

  delete[] bulkElementInt2Ext;
  delete[] bulkElementExt2Int;

  // Node related data
  delete[] nodes;
  delete[] nodeExt2Int;
  delete[] nodeInt2Ext;

  // Body related
  delete[] bodyInt2Ext;
  delete[] bodyExt2Int;
}


void
MeshData::resizeFlagTable(bool*& table, int old_size, int new_size, bool* copy_flags)
{
  int i;
  bool* new_table = new bool[new_size];

  int copy_size = (new_size > old_size) ? old_size: new_size;

  for (i = 0; i < copy_size; i++) {

    if ( copy_flags == NULL || copy_flags[i] ) {
      new_table[i] = table[i];
    }
  }

  for (i = old_size; i < new_size; i++) {
      new_table[i] = false;
  }

  delete[] table;

  table = new_table;
}


void
MeshData::resizeIdTable(int*& table, int old_size, int new_size, bool* copy_flags)
{
  int i;
  int* new_table = new int[new_size];

  int copy_size = (new_size > old_size) ? old_size: new_size;

  for (i = 0; i < copy_size; i++) {

    if ( copy_flags == NULL || copy_flags[i] ) {
      new_table[i] = table[i];
    }
  }

  for (i = old_size; i < new_size; i++) {
      new_table[i] = NO_INDEX;
  }

  delete[] table;

  table = new_table;
}


void
MeshData::resizeNodeTable(int old_size, int new_size, bool* copy_flags)
{
  int i;
  Point3* new_table = new Point3[new_size];

  int copy_size = (new_size > old_size) ? old_size: new_size;

  for (i = 0; i < copy_size; i++) {

    if ( copy_flags == NULL || copy_flags[i] ) {
      new_table[i][0] = nodes[i][0];
      new_table[i][1] = nodes[i][1];
      new_table[i][2] = nodes[i][2];
    }
  }

  for (i = old_size; i < new_size; i++) {
      new_table[i][0] = 0.0;
      new_table[i][1] = 0.0;
      new_table[i][2] = 0.0;
  }

  delete[] nodes;

  nodes = new_table;
}


// *********
// MESH INFO
// *********

// Method initializes model's mesh info structure
MeshInfo::MeshInfo()
{
  selectedBndrElementId = NO_INDEX;
  selectedBulkElementId = NO_INDEX;
  nofAllocElements = 0;
  nofCornerElements = 0;
  nofZeroVelocityElements = 0;
  nofSplittedElements = 0;
  nofBulkEdges = 0;
  nofBulkElements = 0;
  nofVertices = 0;
  maxExternalElementId = NO_INDEX;
  nofBodies = 0;
  nofAllocNodes = 0;
  nofNodes = 0;
  maxExternalNodeId = NO_INDEX;
  nofBoundaryEdges = 0;
  nofBoundaryElements = 0;
  nofBoundaryVertices = 0;
  nofInputBoundaryElements = 0;
  nofInnerBndrElements = 0;
  nofOuterBndrElements = 0;
  minX = minY = minZ = MAX_RANGE;
  maxX = maxY = maxZ = MIN_RANGE;

  nofUsedElementTypes = new int[MAX_NOF_ELEM_CODES];
  for (int i = 0; i < MAX_NOF_ELEM_CODES; i++) {
    nofUsedElementTypes[i] = 0;
  }
}


MeshInfo::~MeshInfo()
{
  delete[] nofUsedElementTypes;
}


// **********
// MODEL DATA
// **********

VertexTable::VertexTable()
{
  dim1 = 0;
  dim2 = 0;
  vertexIds = NULL;
}

VertexTable::~VertexTable()
{
  delete[]vertexIds;
  purgeMatcValueTable(matcTable);
}


ModelData::ModelData()
{
  modelObjects = new ModelObjectArray;
  modelObjects->push_back(NULL); // Seed inserted because object ids start from 1

  modelPoints = new PointHashTable;
  modelPoint2Vertices = new Point2VertexTable;
  lastPointId = 0;

  coveringElementTable = new IdListTable;
  swapElementTable = new IdNumberTable;

  meshCornerElements = new MeshCornerElementList;

  createdModelElements = new BodyElementTable;
  removedModelElements = new BodyElementTable;

  splitCombineInfos = NULL;
  splitCombineInfoIndex = NO_INDEX;

  vertexTable = new VertexTable;

  // Parameter tables
  bodyParameters = new ParameterTable;
  boundaryParameters = new ParameterTable;
  bodyForces = new ParameterTable;
  boundaryConditions = new ParameterTable;
  calculators = new ParameterTable;
  constants = new ParameterTable;
  coordinates = new ParameterTable;
  datafiles = new ParameterTable;
  equations = new ParameterTable;
  equationVariables = new ParameterTable;
  gridHs = new ParameterTable;
  gridParameters = new ParameterTable;
  initialConditions = new ParameterTable;
  materials = new ParameterTable;
  modelParameters = new ParameterTable;
  simulationParameters = new ParameterTable;
  solvers = new ParameterTable;
  solverControls = new ParameterTable;
  timesteps = new ParameterTable;
  userSettings = new ParameterTable;
}


ModelData::~ModelData()
{
  purgeBodyElementTable(removedModelElements);
  purgeModelPoints();

  // Parameter tables
  purgeParameterTable(bodyParameters);
  purgeParameterTable(boundaryParameters);
  purgeParameterTable(bodyForces);
  purgeParameterTable(boundaryConditions);
  purgeParameterTable(calculators);
  purgeParameterTable(constants);
  purgeParameterTable(coordinates);
  purgeParameterTable(datafiles);
  purgeParameterTable(equations);
  purgeParameterTable(equationVariables);
  purgeParameterTable(gridHs);
  purgeParameterTable(gridParameters);
  purgeParameterTable(initialConditions);
  purgeParameterTable(materials);
  purgeParameterTable(modelParameters);
  purgeParameterTable(simulationParameters);
  purgeParameterTable(solvers);
  purgeParameterTable(solverControls);
  purgeParameterTable(timesteps);
  purgeParameterTable(userSettings);

  purgeModelObjects();
  purgeMeshCornerElements();

  delete vertexTable;
}


void
ModelData::purgeBodyElementTable(BodyElementTable*& table)
{
  BodyElementTable::iterator end_pos = table->end();
  BodyElementTable::iterator pos = table->begin();

  while (pos != end_pos) {
    delete (*pos).second;
    pos++;
  }

  delete table;
  table = NULL;
}


void
ModelData::purgeCreatedElements()
{
  purgeBodyElementTable(createdModelElements);
}


void
ModelData::purgeRemovedElements()
{
  purgeBodyElementTable(removedModelElements);
}


void
ModelData::purgeMeshCornerElements()
{
  MeshCornerElementList::iterator end_pos = meshCornerElements->end();
  MeshCornerElementList::iterator pos = meshCornerElements->begin();

  while (pos != end_pos) {
    delete (*pos);
    pos++;
  }

  delete meshCornerElements;
  meshCornerElements = NULL;
}


void
ModelData::purgeModelObjects()
{
  ModelObjectArray::iterator end_pos = modelObjects->end();
  ModelObjectArray::iterator pos = modelObjects->begin();

  while (pos != end_pos) {
    delete (*pos);
    pos++;
  }

  delete modelObjects;
  modelObjects = NULL;
}


void
ModelData::purgeModelPoints()
{
  PointHashTable::iterator end_pos = modelPoints->end();
  PointHashTable::iterator pos = modelPoints->begin();

  while (pos != end_pos) {
    delete (*pos).second;
    //delete (*pos);
    pos++;
  }

  delete modelPoints;
  modelPoints = NULL;
}


void
ModelData::purgeParameterTable(ParameterTable*& table)
{
  ParameterTable::iterator end_pos = table->end();
  ParameterTable::iterator pos = table->begin();

  while (pos != end_pos) {
    delete (*pos).second;
    pos++;
  }

  delete table;
  table = NULL;
}


void
ModelData::reset()
{
  purgeBodyElementTable(createdModelElements);
  createdModelElements = new BodyElementTable;

  purgeBodyElementTable(removedModelElements);
  removedModelElements = new BodyElementTable;
}


// Store vertex table in model data (from egf/emf files)
void
ModelData::setVertexTable(int dim1, int dim2, int* vertex_ids, MatcValueTable& matc_table)
{
  vertexTable->dim1 = dim1;
  vertexTable->dim2 = dim2;
  delete[] vertexTable->vertexIds;
  vertexTable->vertexIds = NULL;

  if ( dim1 > 0 ) {
    vertexTable->vertexIds = new int[dim1];
    for (int i = 0; i < dim1; i++) {
      vertexTable->vertexIds[i] = vertex_ids[i];
    }
  }

  purgeMatcValueTable(vertexTable->matcTable);
  copyMatcValueTable(matc_table, vertexTable->matcTable);
}


#if 0
template <class T> void
ModelData::purgeContainer(T*& container)
{
  T::iterator end_pos = container->end();
  T::iterator pos = container->begin();

  while (pos != end_pos) {
    delete (*pos);
    pos++;
  }

  delete container;
  container = NULL;
}


template <class T> void
ModelData::purgeIdContainer(T*& container)
{
  T::iterator end_pos = container->end();
  T::iterator pos = container->begin();

  while (pos != end_pos) {
    delete (*pos).second;
    pos++;
  }

  delete container;
  container = NULL;
}
#endif


// **********
// MODEL INFO
// **********

ModelInfo::ModelInfo(char* model_name, ecif_modelSource source, char* in_file_name)
{
  create_dyna_string(modelFileTs, NULL);

  frontPreviousInputVersionNbr = -1;
  frontInputVersionNbr = -1;
  frontVersionNbr = ECIF_VERSION_NBR;

  modelStatus = STATUS_OK;
  modelSourceType = source;
  updateGeometryType();

  updateCoordinateMapping("1 2 3");
  updateCoordinateType("Cartesian");
  updateSimulationDimension("Cartesian");

  dimension = ECIF_ND;
  simulationDimension = ECIF_ND;

  hasUserDefinitions = false;
  hasMatcDefinitions = true;
  hasGeometry = false;

  keepMatcDefinitions = true;

  editingMeshBoundaries = false;
  selectedBodyElementId = NO_INDEX;
  selectedBody1Id = NO_INDEX;
  selectedBody2Id = NO_INDEX;
  selectedLayer1Id = NO_INDEX;
  selectedLayer2Id = NO_INDEX;
  minEdgeSize = MAX_RANGE;
  dimAvg = dimMax = dimMin = 0.0;
  minX = minY = minZ = MAX_RANGE;
  maxX = maxY = maxZ = MIN_RANGE;

  //--Database stuff
  create_dyna_string(meshDirectory, NULL);
  create_dyna_string(meshDirectory_absolute, NULL);

  create_dyna_string(modelDirectory, NULL);
  create_dyna_string(modelDirectory_absolute, NULL);

  create_dyna_string(modelName, model_name);
  create_dyna_string(modelFileName, in_file_name);

  //--History stuff
  create_dyna_string(created, NULL);
  create_dyna_string(modified, NULL);

  //--Model cad-sourcefile, mesh-sourcefile, mesh-resultfile
  create_dyna_string(cadSourceFile, NULL);
  create_dyna_string(meshSourceFile, NULL);
  create_dyna_string(meshResultFile, NULL);
  if (in_file_name != NULL) {
    if (modelSourceType == ECIF_CAD_FILE)
      update_dyna_string(cadSourceFile, in_file_name);
    if (modelSourceType == ECIF_MESH_FILE)
      update_dyna_string(meshSourceFile, in_file_name);
  }

  //--Model description
  create_dyna_string(modelDescription, NULL);

  //--Matc input file names
  create_dyna_string(matcInputFile_emf, NULL);
  create_dyna_string(matcInputFile_sif, NULL);

  // Mesh stuff
  nofMeshes = 0;
  meshNames = NULL;
  meshFs = NULL;
  meshHs = NULL;
  nofBgMeshFiles = 0;
  meshBgMeshFileIndices = NULL;
  meshBgMeshFiles = NULL;
  meshBgMeshActives = NULL;
  meshBgMeshControls = NULL;
  nofActiveMeshes = 0;
  activeMeshIndices = NULL;
  currentMeshIndex = NO_INDEX;

  //--Problem stuff
  create_dyna_string(problemDescription, NULL);
  create_dyna_string(problemName, NULL);

  //--File path stuff
  create_dyna_string(includePath, NULL);
  create_dyna_string(includePath_absolute, NULL);
  includePath_save = false;

  create_dyna_string(resultsDirectory, NULL);
  create_dyna_string(resultsDirectory_absolute, NULL);
  resultsDirectory_save = false;

  create_dyna_string(temporaryFilesDirectory, NULL);
  create_dyna_string(temporaryFilesDirectory_absolute, NULL);
  temporaryFilesDirectory_save = false;

  //--Model file timestamps
  create_dyna_string(databaseTs, NULL);
  create_dyna_string(frontTs, NULL);
  create_dyna_string(gebhardtFactorsTs, NULL);
  create_dyna_string(meshParameterTs, NULL);
  create_dyna_string(meshTs, NULL);
  create_dyna_string(solverTs, NULL);
  create_dyna_string(viewfactorsTs, NULL);

  //--Session flags
  databaseNeedsUpdate = false;
  hasDiffuseGrayRadiation = false;
  gebhardtFactorsNeedsUpdate = false;
  meshNeedsUpdate = false;
  solverNeedsUpdate = false;
  viewfactorsNeedsUpdate = false;
  readingModelFile = false;
}


ModelInfo::~ModelInfo()
{
  delete[] modelFileTs;

  //--Datbase stuff
  delete[] meshDirectory;
  delete[] meshDirectory_absolute;
  delete[] modelDirectory;
  delete[] modelDirectory_absolute;
  delete[] modelName;
  delete[] modelFileName;

  //--History stuff
  delete[] created;
  delete[] modified;

  //--Model cad-sourcefile, mesh-sourcefile, mesh-resultfile
  delete[] cadSourceFile;
  delete[] meshSourceFile;
  delete[] meshResultFile;

  //--Model description
  delete[] modelDescription;

  //--Matc input file name
  delete[] matcInputFile_emf;
  delete[] matcInputFile_sif;

  //--Mesh stuff
  for (int i = 0; i < nofMeshes; i++) {
    delete[] meshNames[i];
  }
  delete[] meshNames;
  delete[] meshFs;
  delete[] meshHs;
  delete[] activeMeshIndices;

  delete[] meshBgMeshFileIndices;
  delete[] meshBgMeshFiles;
  delete[] meshBgMeshActives;
  delete[] meshBgMeshControls;

  //--Problem stuff
  delete[] problemDescription;
  delete[] problemName;

  //--File directories
  delete[] includePath;
  delete[] includePath_absolute;
  delete[] resultsDirectory;
  delete[] resultsDirectory_absolute;
  delete[] temporaryFilesDirectory;
  delete[] temporaryFilesDirectory_absolute;

  //--Model file timestamps
  delete[] databaseTs;
  delete[] frontTs;
  delete[] gebhardtFactorsTs;
  delete[] meshParameterTs;
  delete[] meshTs;
  delete[] solverTs;
  delete[] viewfactorsTs;
}


void
ModelInfo::updateCoordinateMapping(const char* coordinate_mapping_str)
{
  coordinateMapping[0] = 1;
  coordinateMapping[1] = 2;
  coordinateMapping[2] = 3;

  if ( coordinate_mapping_str == NULL )
    return;

  strstream strm;
  strm << coordinate_mapping_str;

  int index = 0;

  while ( !strm.eof() ) {
    strm >> coordinateMapping[index++];
  }

}


modelCoordinateType
ModelInfo::updateCoordinateType(const char* coordinate_system)
{
  // Cartesian systems
  if ( 0 == strcmp(coordinate_system, "Cartesian")    ||
       0 == strcmp(coordinate_system, "Cartesian 2D") ||
       0 == strcmp(coordinate_system, "Cartesian 3D")
     ) {
    coordinateType = COORD_CARTESIAN;

  // Cylindric systems
  } else if ( 0 == strcmp(coordinate_system, "Axi Symmetric")       ||
              0 == strcmp(coordinate_system, "Cylindric Symmetric") ||
              0 == strcmp(coordinate_system, "Cylindrical")
            ) {
    coordinateType = COORD_CYLINDRIC;

  // Polar (spherical) systems
  } else {
    coordinateType = COORD_POLAR;
  }

  return coordinateType;
}


modelGeometryType
ModelInfo::updateGeometryType()
{
  // Mesh geometry type
  switch (modelSourceType) {
  case ECIF_CAD_OR_MESH_FILE:
    geometryType = GEOM_CAD_OR_MESH;
    break;
  case ECIF_CAD_FILE:
    geometryType = GEOM_CAD;
    break;
  case ECIF_MESH_FILE:
    geometryType = GEOM_MESH;
    break;
  case ECIF_CAD_AND_MESH_FILE:
    geometryType = GEOM_CAD_AND_MESH;
    break;
  }

  return geometryType;
}


ecif_modelDimension
ModelInfo::updateSimulationDimension(const char* coordinate_system)
{
  // 2D simulation
  if ( 0 == strcmp(coordinate_system, "Cartesian")     ||
       0 == strcmp(coordinate_system, "Cartesian 2D")  ||
       0 == strcmp(coordinate_system, "Axi Symmetric") ||
       0 == strcmp(coordinate_system, "Polar 2D")
     ) {
    simulationDimension = ECIF_2D;

  // 3D simulation
  } else {
    simulationDimension = ECIF_3D;
  }

  return simulationDimension;
}


// ****************
// MODEL STATISTICS
// ****************

ModelStatistics::ModelStatistics()
{
  nofBodies = 0;
  nofBodiesWithBodyForce = 0;
  nofBodiesWithEquation = 0;
  nofBodiesWithInitialCondition = 0;
  nofBodiesWithMaterial = 0;
  nofBodyPairs = 0;
  nofBoundaryPoints = 0;
  nofEdgesWithCondition = 0;
  nofElements = 0;
  nofElementGroups = 0;
  nofElementLoops = 0;
  nofInnerBoundaries = 0;
  nofInnerBoundariesWithCondition = 0;
  nofMeshBodies = 0;
  nofMeshBodyElements = 0;
  nofMeshBoundaryVertices = 0;
  nofMaterials = 0;
  nofModelObjects = 0;
  nofOuterBoundaries = 0;
  nofOuterBoundariesWithCondition = 0;
  nofPatterns = 0;
  nofPoints = 0;
  nofVertices = 0;
  nofVerticesWithCondition = 0;
  maxLoopSize =0;
  // Parameter objects in use (by bodies etc.)
  nofBodyForces = 0;
  nofBodyParameters = 0;
  nofBoundaryConditions = 0;
  nofBoundaryParameters = 0;
  nofCalculators = 0;
  nofConstants = 0;
  nofCoordinates = 0;
  nofDatafiles = 0;
  nofEquations = 0;
  nofEquationVariables = 0;
  nofGridHs = 0;
  nofGridParameters = 0;
  nofInitialConditions = 0;
  nofModelParameters = 0;
  nofSimulationParameters = 0;
  nofSolvers = 0;
  nofSolverControls = 0;
  nofTimesteps = 0;
  nofUserSettings = 0;
}


ModelStatistics::~ModelStatistics()
{
}


// *************
// PARALLEL INFO
// *************

ParallelInfo::ParallelInfo()
{
  nofProcessors = 1;
}

ParallelInfo::~ParallelInfo()
{
}


// ******************
// SPLIT/COMBINE INFO
// ******************

SplitCombineInfo::SplitCombineInfo()
{
  canRedo = false;
  canUndo = false;
  splitSourceId = NO_INDEX;
  splitTargetId = NO_INDEX;
  nofSplitTargetMeshElements = 0;
  splitTargetMeshElementSourceIndices = NULL;
  combineTargetId  = NO_INDEX;
  nofCombineSourceIds = 0;
  combineSourceIds = NULL;
}

SplitCombineInfo::~SplitCombineInfo()
{
  delete[] splitTargetMeshElementSourceIndices;
  delete[] combineSourceIds;
}
