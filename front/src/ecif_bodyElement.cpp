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
Module:     ecif_bodyelement.cpp
Language:   C++
Date:       15.10.98
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
#include "ecif_bodyElementLoop.h"
#include "ecif_bodyElementGroup.h"
#include "ecif_boundaryCondition.h"
#include "ecif_geometry.h"
#include "ecif_mesh.h"
#include "ecif_model.h"
#include "ecif_parameter.h"
#include "ecif_renderer.h"
#include "ecif_userinterface.h"

//Initialize static class variables.
int BodyElement::last_boundaryTag = 0;


// BoundaryPoint struct
// ===================

BoundaryPoint::BoundaryPoint() {
  activeInMeshing = false;
  checked = false;
  tag = NO_INDEX;
  //boundaryTag = NO_INDEX;
  point = NULL;
  vertexTag = NO_INDEX;
  meshDensityType = ' ';
  meshDensityValue = 0.0;
}

BoundaryPoint::~BoundaryPoint() {};

void
BoundaryPoint::copy(BoundaryPoint& other) {
  activeInMeshing = other.activeInMeshing;
  checked = other.checked;
  meshDensityType = other.meshDensityType;
  meshDensityValue = other.meshDensityValue;
  tag = other.tag;
  vertexTag = other.vertexTag;
  point->setPos(X, other.point->Pos(X));
  point->setPos(Y, other.point->Pos(Y));
  point->setPos(Z, other.point->Pos(Z));
}


// GridHData struct
// ==============

BodyElementGridHData::BodyElementGridHData()
{
  nofIds = 0;
  gridHIds = NULL;
  meshIndices = NULL;

  meshIndex = -1;
  hValueType = ' ';
  hValue = 0.0;
}


BodyElementGridHData::BodyElementGridHData(int nof_ids, int* grid_h_ids, int* mesh_indices)
{
  nofIds = nof_ids;

  meshIndex = -1;
  hValueType = ' ';
  hValue = 0.0;

  if ( nofIds <= 0)
    return;

  gridHIds = new int[nofIds];
  meshIndices = new int[nofIds];

  for (int i = 0; i < nofIds; i++) {
    gridHIds[i] = grid_h_ids[i];
    meshIndices[i] = mesh_indices[i];
  }
}


BodyElementGridHData::~BodyElementGridHData()
{
  delete[] gridHIds;
  delete[] meshIndices;
}


void
BodyElementGridHData::setData(int nof_ids, int* grid_h_ids, int* mesh_indices, bool force)
{
  int i;

  nofIds = nof_ids;

  // If directly set density values for the mesh already exist and parameters
  // are to be forced, delete direct definitions!
  for ( i = 0; i < nofIds; i++) {
    if ( mesh_indices[i] == meshIndex && force ) {
      meshIndex = -1;
      hValueType = ' ';
      hValue = 0.0;
      break;
    }
  }

  delete[] gridHIds;
  delete[] meshIndices;

  gridHIds = NULL;
  meshIndices = NULL;

  if ( nofIds <= 0)
    return;

  gridHIds = new int[nofIds];
  meshIndices = new int[nofIds];

  for (i = 0; i < nofIds; i++) {
    gridHIds[i] = grid_h_ids[i];
    meshIndices[i] = mesh_indices[i];
  }
}


void
BodyElementGridHData::setData(int mesh_index, char value_type, double value, bool force)
{
  int i;

  // If grid parameter for the mesh already exista and values
  // not to forced (overwrite), do nothing!
  for ( i = 0; i < nofIds; i++) {
    if ( meshIndices[i] == mesh_index && !force ) return;
  }

  meshIndex = mesh_index;
  hValueType = value_type;
  hValue = value;
}


// QuadGridData struct
// ===================

BodyElementQuadGridData::BodyElementQuadGridData()
{
  nofIds = 0;
  meshIndices = NULL;
  nValues = NULL;
}


BodyElementQuadGridData::BodyElementQuadGridData(int nof_ids, int* mesh_indices, int* n_values)
{
  nofIds = nof_ids;

  if ( nofIds <= 0)
    return;

  meshIndices = new int[nofIds];
  nValues = new int[nofIds];

  for (int i = 0; i < nofIds; i++) {
    meshIndices[i] = mesh_indices[i];
    nValues[i] = n_values[i];
  }
}


BodyElementQuadGridData::~BodyElementQuadGridData()
{
  delete[] meshIndices;
  delete[] nValues;
}


void
BodyElementQuadGridData::setData(int nof_ids, int* mesh_indices, int* n_values)
{
  nofIds = nof_ids;

  delete[] meshIndices;
  delete[] nValues;

  meshIndices = NULL;
  nValues = NULL;

  if ( nofIds <= 0)
    return;

  meshIndices = new int[nofIds];
  nValues = new int[nofIds];

  for (int i = 0; i < nofIds; i++) {
    meshIndices[i] = mesh_indices[i];
    nValues[i] = n_values[i];
  }
}


// ModelData struct
// ================

BodyElementModelData::BodyElementModelData()
{
  status = 0;

  nofBoundaryPoints = 0;
  nofBoundaryPointsU = 0;
  nofBoundaryPointsV = 0;
  boundaryPoints = NULL;

  nofSubElements = 0;
  subElementIds = new IdArray;

  nofVertices = 0;
  vertexIds = new IdArray;

  code = 0;

  isClosedU = false;
  isClosedV = false;
  isIntraLayerBoundary = false;

  parent1Id = NO_INDEX;
  parent2Id = NO_INDEX;

  parent1Tag = NO_INDEX;
  parent2Tag = NO_INDEX;

  parent1Layer = NO_INDEX;
  parent2Layer = NO_INDEX;
}


BodyElementModelData::~BodyElementModelData()
{
  deleteBoundaryPoints();
  delete subElementIds;
  delete vertexIds;
}


void
BodyElementModelData::deleteBoundaryPoints()
{
  // In closed boundaries, start and end points
  // are same pointers ==> we cannot delete them twice!!!

  int nof_points_u = nofBoundaryPointsU;
  int nof_points_v = nofBoundaryPointsV;

  if ( isClosedU ) {
    nof_points_u--;
  }

  if ( isClosedV ) {
    nof_points_v--;
  }

  // Find nof-points to be deleted
  // =============================
  int nof_points;

  if ( nof_points_v > 0 )
    nof_points = nof_points_u * nof_points_v;
  else
    nof_points = nof_points_u;

  // Delete point coordinates (but not if original vertex!)
  // ==============================================
  for (int i = 0; i < nof_points; i++) {

    BoundaryPoint* bp = boundaryPoints[i];

    // We can delete "internal" points, but NOT the original
    // vertex points!!!
    if ( !bp->isVertex() ) {
      delete bp->point;
    }

    delete bp;
  }

  // Delete table itself
  // ===================
  nofBoundaryPoints = 0;
  nofBoundaryPointsU = 0;
  nofBoundaryPointsV = 0;

  delete[] boundaryPoints;
  boundaryPoints = NULL;
}


// Return subelement tags int the buffer "element_tags"
// NOTE: Buffer is allocated here, but each client should delete it!
bool
BodyElementModelData::getSubElementTags(Model* model, objectType sub_type, int& count,
                                        int*& element_tags)

{
  count = 0;
  element_tags = NULL;

  if ( nofSubElements == 0 ) return false;

  int i;

  // Count sub-elements of proper type
  for (i = 0; i < nofSubElements; i++) {
    int element_id = (*subElementIds)[i];
    BodyElement* se = model->getBodyElementById(element_id);
    if ( se != NULL && sub_type == se->getObjectType() ) {
      count++;
    }
  }

  element_tags = new int[count];

  // Collect sub-element tags of proper type
  for (i = 0; i < nofSubElements; i++) {
    int element_id = (*subElementIds)[i];
    BodyElement* se = model->getBodyElementById(element_id);
    if ( se != NULL && sub_type == se->getObjectType() ) {
      element_tags[i] = se->Tag();
    }
  }

  return true;
}


// Return vertex tags int the buffer "vertex_tags"
// NOTE: Buffer is allocated here, but each client should delete it!
bool
BodyElementModelData::getVertexTags(Model* model, int& count, int*& vertex_tags)
{
  count = nofVertices;
  vertex_tags = NULL;

  if ( count == 0 )
    return false;

  vertex_tags = new int[count];

  for (int i = 0; i < count; i++) {

    int vertex_id = (*vertexIds)[i];

    BodyElement* v = model->getVertexById(vertex_id);

    if ( v != NULL )
      vertex_tags[i] = v->Tag();
    else
      vertex_tags[i] = NO_INDEX;
  }

  return true;

}


// MeshData struct
// ===============

BodyElementMeshData::BodyElementMeshData()
{
  maxNofMeshElements = 0;
  nofMeshElements = 0;
  meshElementIds  = NULL;
  meshElementDirs = NULL;

  maxNofPendingMeshElements = 0;
  nofPendingMeshElements = 0;
  pendingMeshElementInfos = NULL;

  nofMeshSelectedElements = 0;

  nofMeshBorderElements = 0;
  meshBorderElementIds  = NULL;
  meshBorderElementDirs = NULL;

  backupData = NULL;
}


BodyElementMeshData::~BodyElementMeshData()
{
  delete[] meshElementIds;
  delete[] meshElementDirs;

  delete[] meshBorderElementIds;
  delete[] meshBorderElementDirs;

  for ( int i = 0; i < maxNofPendingMeshElements; i++) {
    delete[] pendingMeshElementInfos[i];
  }
  delete[] pendingMeshElementInfos;

  delete backupData;
}


// MeshDataBackup struct
// =====================

BodyElementMeshDataBackup::BodyElementMeshDataBackup()
{
  nofMeshElements = 0;
  meshElementIds = NULL;
  meshElementDirs = NULL;

  nofMeshBorderElements = 0;
  meshBorderElementIds = NULL;
  meshBorderElementDirs = NULL;
};


BodyElementMeshDataBackup::~BodyElementMeshDataBackup()
{
  delete[] meshElementIds;
  delete[] meshElementDirs;

  delete[] meshBorderElementIds;
  delete[] meshBorderElementDirs;
};



// =================
// BodyElement class
// =================

BodyElement::BodyElement()
{
  init();
}


BodyElement::BodyElement(int int_tag)
{
  tag = int_tag;

  init();
}


BodyElement::BodyElement(int cd, char* nm)
{
  init();

  if ( modelData != NULL )
    modelData->code = cd;

  update_dyna_string(name, nm);
  update_dyna_string(name_old, nm);
}


//-Create a bodyelement from a modelfile bodyelement
BodyElement::BodyElement(ecif_Element_X& tx)
{
  tag = tx.tag;

  init();


  boundaryTag = tx.bndr_tag;
  checkLastBoundaryTag();

  // NOTE: Boundary conditions are stored
  // finally in element groups
  //
  txBndrConditionId = tx.bndr_cond_id;

  boundaryParameterId = tx.bndr_param_id;

  elementGroupTag = tx.bndr_group_tag;

  if ( tx.nof_gridh_ids > 0 ) {
    gridHData = new BodyElementGridHData(tx.nof_gridh_ids, tx.gridh_ids, tx.gridh_mesh_indices);
  }

}


//-Create a bodyelement from a modelfile vertex
BodyElement::BodyElement(ecif_Vertex_X& tx)
{
  tag = tx.tag;
  init();

  boundaryTag = tx.bndr_tag;
  elementGroupTag = tx.bndr_group_tag;

  checkLastBoundaryTag();

  txBndrConditionId = tx.bndr_cond_id;

  if ( tx.nof_grid_h_ids > 0 ) {
    gridHData = new BodyElementGridHData(tx.nof_grid_h_ids, tx.grid_h_ids, tx.mesh_indices);
  }

}


// NOTE: ptrGmtr is deleted in the implementation
// classes BodyElement3D ..., cause they know when to
// delete and when not (those vertices...!!!)
//
BodyElement::~BodyElement()
{
  delete[] name_old;
}


int
BodyElement::addAllPendingMeshElements()
{
  MeshInfo* mi = (MeshInfo*)model->getMeshInfo();

  int element_count = 0;

  if ( getObjectType() != OT_EDGE ||
       meshData == NULL ||
       meshData->nofPendingMeshElements == 0
     ) {
    return 0;
  }

  MeshElementTable* source = model->getMeshBoundaryElementEdges();

  for ( int i = 0; i < meshData->nofPendingMeshElements; i++) {

    int* info = meshData->pendingMeshElementInfos[i];

    int elem_code = info[0];

    info++;

    int nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];

    short dir, start;

    for ( int j = 0; j < source->NofElements(); j++) {

      if ( source->isSameElement(j, dir, start, elem_code, nof_nodes, info) ) {

        addMeshElement(j, dir);

        mi->nofUsedElementTypes[elem_code]++;

        element_count++;

        if ( element_count == meshData->nofPendingMeshElements ) {
          break;
        }

      }
    }

    if ( element_count == meshData->nofPendingMeshElements ) {
          break;
    }
  }


  return element_count;
}


// Add one mesh (fem) element id
// Returns the internal index of the new mesh element
int
BodyElement::addMeshElement(int elem_id, short direction)
{
  if ( meshData == NULL )
    return 0;

  if (meshData->meshElementIds == NULL)
    allocateMeshElements(meshData->maxNofMeshElements);

  if (meshData->nofMeshElements == meshData->maxNofMeshElements)
    reallocateMeshElements();

  meshData->meshElementDirs[meshData->nofMeshElements] = direction;

  meshData->meshElementIds[meshData->nofMeshElements] = elem_id;

  return meshData->nofMeshElements++;
}


int
BodyElement::addPendingMeshElementAsNodes(int element_type, int* ext_node_ids)
{
  if ( meshData == NULL )
    return 0;

  if ( meshData->pendingMeshElementInfos == NULL )
    allocatePendingMeshElementInfos(meshData->maxNofPendingMeshElements);

  if (meshData->nofPendingMeshElements == meshData->maxNofPendingMeshElements)
    reallocatePendingMeshElementInfos();

  meshElementCode el_code = model->convertElementType(element_type);
  int nof_nodes = MeshElementDesc[el_code][DESC_NOF_NODES];

  int* node_info = new int[1 + nof_nodes];

  node_info[0] = el_code;

  for (int i = 0; i < nof_nodes; i++) {
    node_info[1 + i] = model->getMeshNodeIdExt2Int(ext_node_ids[i]);
  }

  int index = meshData->nofPendingMeshElements;

  meshData->pendingMeshElementInfos[index] = node_info;

  return ++meshData->nofPendingMeshElements;
}


beStatus
BodyElement::addStatus(beStatus new_status)
{
  if ( modelData == NULL )
    return BE_NONE;

  else
    return modelData->status |= new_status;
}


// Add sublement ( Edge for a Face, Veretx for an Edge)
int
BodyElement::addSubElement(BodyElement* be)
{
  return addSubElement(be->Id());
}


// Add sublement ( Edge for a Face, Vertex for an Edge)
int
BodyElement::addSubElement(int sub_id)
{
  if ( modelData == NULL )
    return NO_INDEX;

  if ( modelData->subElementIds->size() == 0 ||
       (int)NULL == (*(modelData->subElementIds))[sub_id]
     ) {
    modelData->subElementIds->push_back(sub_id);
    modelData->nofSubElements++;
  }

  // Return index for the element, so that
  // clients can use it
  return modelData->nofSubElements - 1;
}


// Allocate space for fem elements
void
BodyElement::allocateMeshElements(int nof_elements)
{
  delete[] meshData->meshElementDirs; meshData->meshElementDirs = NULL;
  delete[] meshData->meshElementIds; meshData->meshElementIds = NULL;

  if (nof_elements > 0) {
    meshData->meshElementDirs = new short[nof_elements];
    meshData->meshElementIds = new int[nof_elements];
  }

  meshData->maxNofMeshElements = nof_elements;

  // No elements stored yet!
  meshData->nofMeshElements = 0;
}


// Allocate space for pending fem element infos
void
BodyElement::allocatePendingMeshElementInfos(int nof_infos)
{
  delete[] meshData->pendingMeshElementInfos;
  meshData->pendingMeshElementInfos = NULL;

  if (nof_infos > 0) {
    meshData->pendingMeshElementInfos = new int*[nof_infos];

    for (int i = 0; i < nof_infos; i++) {
      meshData->pendingMeshElementInfos[i] = NULL;
    }
  }

  meshData->maxNofPendingMeshElements = nof_infos;

  // No infos stored yet!
  meshData->nofPendingMeshElements = 0;
}



void
BodyElement::backupMeshData()
{
  if ( meshData == NULL )
    return;

  meshData->backupData = new BodyElementMeshDataBackup;

  int i;

  //--Mesh elements
  meshData->backupData->nofMeshElements = meshData->nofMeshElements;
  meshData->backupData->meshElementIds = new int[meshData->nofMeshElements];
  meshData->backupData->meshElementDirs = new short[meshData->nofMeshElements];

  for (i = 0; i < meshData->nofMeshElements; i++) {
    meshData->backupData->meshElementIds[i] = meshData->meshElementIds[i];
    meshData->backupData->meshElementDirs[i] = meshData->meshElementDirs[i];
  }

  //--Mesh border elements
  meshData->backupData->nofMeshBorderElements = meshData->nofMeshBorderElements;
  meshData->backupData->meshBorderElementIds = new int[meshData->nofMeshBorderElements];
  meshData->backupData->meshBorderElementDirs = new short[meshData->nofMeshBorderElements];

  for (i = 0; i < meshData->nofMeshBorderElements; i++) {
    meshData->backupData->meshBorderElementIds[i] = meshData->meshBorderElementIds[i];
    meshData->backupData->meshBorderElementDirs[i] = meshData->meshBorderElementDirs[i];
  }

}


// Linearize geometry and store linearizing points in the 'bpoints' buffer
// NOTE: Buffer is allocated here!!!
void
BodyElement::calcBoundaryPoints(double delta_u, double delta_v)
{
  if ( modelData == NULL ) {
    return;
  }


  modelData->deleteBoundaryPoints();

  if ( ptrGmtr == NULL ) return;

  ptrGmtr->calcBoundaryPoints();

  return;


  // Old version!!!
  //
  int nof_points_u;
  int nof_points_v;
  GcPoint** points;

  ptrGmtr->calcLinearizingPoints(nof_points_u, nof_points_v, points, delta_u, delta_v);

  int nof_points = 0;

  if ( nof_points_v > 0 )
    nof_points = nof_points_u * nof_points_v;
  else
    nof_points = nof_points_u;

  if (nof_points == 0)
    return;

  modelData->boundaryPoints = new BoundaryPoint*[nof_points];

  for (int i = 0; i < nof_points; i++) {

    BoundaryPoint* bp = new BoundaryPoint;

    modelData->boundaryPoints[i] = bp;

    // Check if bp is an existing vertex!
    BodyElement* v = model->getVertex(points[i]);

    if ( v != NULL ) {
      bp->tag = v->Tag();
      bp->vertexTag = v->Tag();

    }

    //bp->boundaryTag = boundaryTag;
    bp->point = points[i];

  }

  // Handle possibly closed geometry
  // -------------------------------
  int nof_u = nof_points_u;
  int nof_v = (nof_points_v > 0)?nof_points_v:1;

  //-Copy start points to end points if U-closed
  //
  if (ptrGmtr->isClosedU()) {

    for (int i = 1; i <= nof_v; i++) {
      int pos0 = (i - 1) * nof_u;
      int pos1 = pos0 + nof_u - 1;

      BoundaryPoint* tmp = modelData->boundaryPoints[pos1];
      modelData->boundaryPoints[pos1] = modelData->boundaryPoints[pos0];

      delete tmp;
    }
  }

  //-Copy start points to end points if V-closed
  //
  if (ptrGmtr->isClosedV()) {

    for (int i = 1; i <= nof_u; i++) {
      int pos0 = (i - 1) * nof_v;
      int pos1 = pos0 + nof_v - 1;

      BoundaryPoint* tmp = modelData->boundaryPoints[pos1];
      modelData->boundaryPoints[pos1] = modelData->boundaryPoints[pos0];

      delete tmp;
    }
  }

  modelData->nofBoundaryPointsU = nof_points_u;
  modelData->nofBoundaryPointsV = nof_points_v;

  if ( nof_points_v == 0 )
    modelData->nofBoundaryPoints = nof_points_u;
  else
    modelData->nofBoundaryPoints = nof_points_u * nof_points_v;

  delete[] points;

}


// Check boundary group data
//
// NOTE: Each element must belong to some boundary group!!!
// If boundary group data is not complete, make sure that
// element is inserted to the group (the grop will be created
// if it does not exist)
//
bool
BodyElement::checkElementGroupData()
{
  // NOTE: No group is created (explicitely) for an intra layer boundary!
  if ( isIntraLayerBoundary() ) return true;

  if ( elementGroupTag == NO_INDEX ||
       elementGroupId == NO_INDEX
     ) {

    // Get group id, possibly create a new group if tag is no-index
    elementGroupId = model->addToElementGroup(this, elementGroupTag);

    // Pick group
    BodyElementGroup* beg = model->getBodyElementGroupById(elementGroupId);

    if ( beg == NULL ) return false;

    // If a new group was created for this boundary, get the group tag and set also
    // boundaryCondition-id for the group
    // NOTE: This is a single boundary group, only for this element!
    //
    if ( elementGroupTag == NO_INDEX) {
      elementGroupTag = beg->Tag();
      beg->setType(IMPLICIT_GROUP);
      beg->setBoundaryConditionId(txBndrConditionId);

      // Set parent info into group
      beg->setParentId(1,getParentId(1));
      beg->setParentId(2,getParentId(2));
      beg->setParentLayer(1,getParentLayer(1));
      beg->setParentLayer(2,getParentLayer(2));
      beg->setParentTag(1,getParentTag(1));
      beg->setParentTag(2,getParentTag(2));
    }
 }

  return true;
}


// If body element is added or deleted from the model
// we have to check that class attribute "lastBoundaryTag"
// remains correct.
// NOTE: Each subclasses (1D, 2D, 3D) have a common last_boundaryTag
void
BodyElement::checkLastBoundaryTag(bool being_deleted)
{
  // Nothing to do
  if ( boundaryTag == NO_INDEX )
    return;

  // If the last element is being deleted
  if (being_deleted) {
    if (last_boundaryTag == boundaryTag) {
      last_boundaryTag--;
    }
  }

  // If this element has larger id than what
  // corresponds to the current lastId
  else {
    if (last_boundaryTag < boundaryTag) {
      last_boundaryTag = boundaryTag;
    }
  }

}


// If body element is added or deleted from the model
// we have to check that class attribute "last_tag"
// remains correct.
// NOTE: Each subclass (1D, 2D, 3D) has its own tags which
// starts from 1
void
BodyElement::checkLastTag(bool being_deleted)
{
  int last_tag = getLastTag();

  // Nothing to do (no last tag defined!)
  if ( last_tag == NO_INDEX )
    return;

  // If the last element is being deleted
  if (being_deleted) {
    if (last_tag == tag) {
      last_tag--;
    }
  }

  // If this element has larger id than what
  // corresponds to the current lastId
  else {
    if (last_tag < tag) {
      last_tag = tag;
    }
  }

  setLastTag(last_tag);
}

// Check if boundary is an intra-layer boundary
//
// NOTE: Such boundaries are displayed only on GridParameter
// panel for grid controlling purposes
//
// NOTE: They will be however output to Elmer mesh, because Mesh2D
// cannot drop them away!
//
bool
BodyElement::isIntraLayerBoundary()
{
  if ( modelData == NULL )
    return false;
  else
    return modelData->isIntraLayerBoundary;
}


// Check if we have to create at least one vertex for the
// body element (this is mainly needed for mesh geometries
// where we want to have at least one vertex at the boundaries
// for Dirichlet boundary conditions
//
void
BodyElement::createExtraVertices()
{
  // No need to create an extra vertex!
  if ( modelData == NULL ||
       modelData->nofVertices > 0 ||
       modelData->nofSubElements > 0
     ) {
    return;
  }

  // No way to create an extra vertex!
  if ( meshData == NULL ||
       meshData->nofMeshElements == 0
       ) {
    return;
  }

  MeshElementTable* source;
  Point3* nodes = model->getMeshNodeData();
  const int* node_ids;

  objectType otype = getObjectType();
  ecif_modelDimension dim = model->getDimension();

  // 3D face or a normal 2D boundary
  if ( (dim == ECIF_3D && otype == OT_FACE) ||
       (dim == ECIF_2D && otype == OT_EDGE)
     ) {
    source = model->getMeshBoundaryElements();

  // 3D edge boundary or 2D vertex boundary
  } else if ( (dim == ECIF_3D && otype == OT_EDGE) ||
              (dim == ECIF_2D && otype == OT_VERTEX)
     ) {
    source = model->getMeshBoundaryElementEdges();

  } else {
    return;
  }

  // Pick first node from the first element!
  node_ids = source->getNodeIds(meshData->meshElementIds[0]);

  // Something went wrong!
  if ( node_ids == NULL ) {
    return;
  }

  // Check if vertex already exists, and if not
  // create a new vertex
  BodyElement* v = NULL;

  GcPoint gp(nodes[node_ids[0]]);

  // Find model point
  GcPoint* mp = model->getPoint(&gp);

  if ( mp != NULL ) {
    v = model->getVertex(mp);
  }

  if ( v == NULL ) {

    v = new BodyElement1D(&gp);

    v->addMeshElement(node_ids[0], 1);

    model->addBodyElement(v);
  }

  // Add vertex to the body element
  modelData->nofSubElements = 1;
  modelData->subElementIds->push_back(v->Id());
}


bool
BodyElement::convertTags2Ids()
{
  UserInterface* gui = (UserInterface*)model->getGui();
  strstream strm1, strm2;

  if ( modelData == NULL ) {
    strm1 << "***ERROR in definition for the edge " << tag << ":" << ends;
    strm2 << "---No data defined in the input!" << ends;
    gui->showMsg(strm1.str());
    gui->showMsg(strm2.str(), 1);
    return false;
  }

  int i;

  // Subelement tag -> ids
  for (i = 0; i < modelData->nofSubElements; i++) {

    int se_tag = (*modelData->subElementIds)[i];

    BodyElement* se;

    if ( model->getDimension() == ECIF_3D )
      se = model->getEdgeByTag(se_tag);
    else
      se = model->getVertexByTag(se_tag);

    if ( se != NULL ) {
      (*modelData->subElementIds)[i] = se->Id();
    } else {
      (*modelData->subElementIds)[i] = NO_INDEX;
      strm1 << "***ERROR in geometry for the edge " << tag << ":" << ends;
      strm2 << "---Sub element " << se_tag << " not found!" << ends;
      gui->showMsg(strm1.str());
      gui->showMsg(strm2.str(), 1);
      return false;
    }
  }

  // Vertex tags -> ids
  for (i = 0; i < modelData->nofVertices; i++) {

    int v_tag = (*modelData->vertexIds)[i];

    BodyElement* v = model->getVertexByTag(v_tag);

    if ( v != NULL ) {
      (*modelData->vertexIds)[i] = v->Id();
    } else {
      (*modelData->vertexIds)[i] = NO_INDEX;
      strm1 << "***ERROR in geometry for the edge " << tag << ":" << ends;
      strm2 << "---Vertex " << v_tag << " not defined in the input!" << ends;
      gui->showMsg(strm1.str());
      gui->showMsg(strm2.str(), 1);
      return false;
    }
  }

  return true;
}


void
BodyElement::deleteCadData()
{
  delete ptrGmtr; ptrGmtr = NULL;
}


void
BodyElement::deleteMeshDataBackup()
{
  if (meshData == NULL || meshData->backupData == NULL )
    return;

  delete[] meshData->backupData->meshElementIds;
  delete[] meshData->backupData->meshElementDirs;
  delete[] meshData->backupData->meshBorderElementIds;
  delete[] meshData->backupData->meshBorderElementDirs;

  delete meshData->backupData;
  meshData->backupData = NULL;
}


void
BodyElement::deleteMeshElements()
{
  if ( meshData == NULL ) {
    return;
  }

  delete meshData;
  meshData = new BodyElementMeshData;
}


void
BodyElement::deleteMeshElements(bool* delete_flags)
{
  if ( meshData == NULL ) {
    return;
  }

  int i;

  // Count nof elements left
  int nof_elements = 0;

  for (i = 0; i < meshData->nofMeshElements; i++) {

    if ( !delete_flags[i] ) {
      nof_elements++;
    }
  }

  // If no elements left!!!
  if (nof_elements == 0) {
    deleteMeshElements();
    return;
  }

  int* tmp_ids = new int[nof_elements];
  short* tmp_dirs = new short[nof_elements];

  int counter = 0;

  for (i = 0; i < meshData->nofMeshElements; i++) {

    if ( delete_flags[i] ) {
      continue;
    }

    int index = meshData->meshElementIds[i];
    short dir = meshData->meshElementDirs[i];

    tmp_ids[counter] = index;
    tmp_dirs[counter] = dir;
    counter++;
  }

  delete[] meshData->meshElementIds;
  delete[] meshData->meshElementDirs;

  meshData->meshElementIds = tmp_ids;
  meshData->meshElementDirs = tmp_dirs;
  meshData->nofMeshElements = nof_elements;
  meshData->maxNofMeshElements = nof_elements;

  // For sake of sure!
  meshData->nofMeshSelectedElements = 0;
}


// Draw element by geometry component index
//
void BodyElement::draw(int gmtr_index, Renderer* renderer, flagName geometry_type, int body_id, int direction, bool is_first_loop)
{
  if ( isIntraLayerBoundary() ) {
    //drawMode = DM_INTRA_LAYER;
  }

  Body* bd = model->getBodyById(body_id);

  int layer;
  if ( body_id == modelData->parent1Id )
    layer = modelData->parent1Layer;
  else if ( body_id == modelData->parent2Id )
    layer = modelData->parent2Layer;
  else
    layer = NO_INDEX;

  // If this is a hidden boundary or an exluded mesh boundary
  if ( (drawMode == DM_HIDDEN && drawState != DS_SELECTED ) ||
       ( geometry_type == DRAW_SOURCE_MESH &&
         bd->isExcludedFromCurrentMesh(layer)
       )
     ) {
    return;
  }

  // Store element name and draw
  renderer->name_save(id);

  if ( geometry_type == DRAW_SOURCE_MESH ) {
    drawMesh(renderer, body_id, direction);

  } else {
    ptrGmtr->draw(gmtr_index, renderer, drawMode, drawState, direction, is_first_loop);
  }

  // Draw sub elements (edges/vertices)
  if ( gmtr_index == 0 ) {
    drawSubElements(renderer, geometry_type, body_id, direction, is_first_loop);
  }

  renderer->name_delete(id);
}


// Draw element
//
void BodyElement::draw(Renderer* renderer, flagName geometry_type, int body_id, int direction, bool is_first_loop)
{
  if ( isIntraLayerBoundary() ) {
    //drawMode = DM_INTRA_LAYER;
  }

  Body* bd = model->getBodyById(body_id);

  int layer;
  if ( body_id == modelData->parent1Id )
    layer = modelData->parent1Layer;
  else if ( body_id == modelData->parent2Id )
    layer = modelData->parent2Layer;
  else
    layer = NO_INDEX;

  // If this is a hidden boundary or an exluded mesh boundary
  if ( (drawMode == DM_HIDDEN && drawState != DS_SELECTED) ||
       ( geometry_type == DRAW_SOURCE_MESH &&
         bd->isExcludedFromCurrentMesh(layer)
       )
     ) {
    return;
  }

  // Store element name and draw
  renderer->name_save(id);

  if ( geometry_type == DRAW_SOURCE_MESH ) {
    drawMesh(renderer, body_id, direction);

  } else {
    ptrGmtr->draw(renderer, drawMode, drawState, direction, is_first_loop);
  }

  // Draw sub elements (edges/vertices)
  drawSubElements(renderer, geometry_type, body_id, direction, is_first_loop);

  renderer->name_delete(id);
}


void
BodyElement::drawLabel(Renderer* renderer, bool draw_sub_labels)
{
  if ( isIntraLayerBoundary() ) return;

  if (labelData != NULL && model->getLabelDisplayFlagValue(this) ) {

    renderer->drawElementLabel(labelData->label, labelData->position);
  }

  if ( !draw_sub_labels || modelData == NULL ) {
    return;
  }

  for (int i = 0; i < modelData->nofSubElements; i++ ) {

    BodyElement* se = model->getBodyElementById((*modelData->subElementIds)[i]);

    if ( se != NULL && (se->getDrawMode() == DM_NORMAL || se->getDrawState() == DS_SELECTED) ) {
      se->drawLabel(renderer, draw_sub_labels);
    }
  }
}


void
BodyElement::drawSubElements(Renderer* renderer, flagName geometry_type, int body_id, int direction, bool is_first_loop)
{
  int nof_sub_elements = getNofSubElements();

  for (int i = 0; i < nof_sub_elements; i++) {

    BodyElement* se = getSubElement(i);

    if ( se != NULL ) {
      se->draw(renderer, geometry_type, body_id, direction, is_first_loop);
    }
  }
}


void
BodyElement::drawMesh(Renderer* renderer, int body_id, int direction)
{
  if ( meshData == NULL || meshData->nofMeshElements == 0) {
    return;
  }

  if ( renderer->isPickingMesh() && drawState != DS_SELECTED ) {
    return;
  }

  //=============================================

  if ( renderer->isEditingMesh() ) {

    //---Draw border for the element
    if ( ECIF_3D == model->getDimension() ) {
      drawMeshBorder(renderer, true);

    } else {
      drawMeshBorder(renderer, false);
    }
  }

  MeshElementTable* elements;
  MeshElementTable* parents;
  MeshElementTable* et;

  bool isLevel2 = false;

  if ( model->getDimension() == ECIF_3D &&
       getObjectType() == OT_EDGE
     ) {
    isLevel2 = true;
  }

  if ( isLevel2) {
    elements = model->getMeshBoundaryElementEdges();
    parents = NULL;
    et = model->getMeshBoundaryElementVertices();

  } else {
    elements = model->getMeshBoundaryElements();
    parents = model->getMeshBulkElements();
    et = model->getMeshBoundaryElementEdges();
  }

  Point3* node_data = model->getMeshNodeData();

  int i, j, index;
  short dir;


  //=============================================
  //---Surface elements drawing loop (draw "boundaries" filled)
  renderer->startDrawingMeshSurface();

  meshData->nofMeshSelectedElements = 0;

  // Set state flags
  // ===============
  bool filled = false;
  bool selected = false;
  bool two_sided = false;

  //--3D surface mode
  if ( model->getDimension() == ECIF_3D ) {

    two_sided = true;

    if ( model->getFlagValue(DRAW_TARGET_BOUNDARIES) || drawState == DS_SELECTED ) {
      filled = true;
    }

  //--2D "surface" mode
  } else {
    if ( model->getFlagValue(DRAW_TARGET_EDGES) ) {
      filled = true;
    }
  }

  if ( !model->getFlagValue(DRAW_TARGET_BOUNDARIES) && drawState == DS_SELECTED ) {
    two_sided = true;
  }

  if ( two_sided ) {
    renderer->drawTwoSidedPolygons(true);
  }

  if ( !filled) {
    renderer->setLightUndirected();
  }

  // Loop mesh elements in the boundary
  // ==================================
  for (i = 0; i < meshData->nofMeshElements; i++) {

    selected = false;

    index = meshData->meshElementIds[i];
    dir = meshData->meshElementDirs[i];

    if ( parents != NULL ) {
      int parent1_elem_id = elements->parentIds[index][0];
      int parent1_body_id = parents->getParentId(parent1_elem_id, 0);

      if ( parent1_body_id != body_id ) {
        dir = -dir;
      }
    }

    // Mesh element's selected state is taken into account only when
    // boundary itself is selected
    if ( drawState == DS_SELECTED ) {

      if ( elements->getSelected(index) ) {

        selected = true;
        meshData->nofMeshSelectedElements++;
      }
    }

    // Draw
    // ====
    if (filled) {
      elements->drawElement(renderer, index, dir, selected);

    } else {
      elements->drawEdges(renderer, et, index, selected);
    }

    elements->resetRendered(index);
  }

  // Reset state flags
  // =================
  if ( two_sided ) {
    renderer->drawTwoSidedPolygons(false);
  }

  if ( !filled) {
    renderer->setLightDirected();
  }

  renderer->stopDrawingMeshSurface();


  //=============================================
  // If boundary is not selected, we can stop here
  if (drawState != DS_SELECTED) {
    return;
  }


  //=============================================
  //---Mark selected element's border edges as
  //   selected (to be drawn thicker!)

  // Set first all edges as non-selected
  for (i = 0; i < meshData->nofMeshElements; i++) {

    index = meshData->meshElementIds[i];
    meshElementCode bndr_code = elements->getElementCode(index);
    int nof_edges = MeshElementDesc[bndr_code][DESC_NOF_EDGES];

    const int* edge_ids = elements->getEdgeIds(index);

    for (j = 0; j < nof_edges; j++) {
      et->setSelected(edge_ids[j], false, false);
    }
  }

  // Then set the edges of the selected elements selected so
  // that the interior edges are left non-selected!
  for (i = 0; i < meshData->nofMeshElements; i++) {

    index = meshData->meshElementIds[i];

    if ( !elements->selected[index] ) {
      continue;
    }

    meshElementCode bndr_code = elements->getElementCode(index);
    int nof_edges = MeshElementDesc[bndr_code][DESC_NOF_EDGES];

    const int* edge_ids = elements->getEdgeIds(index);

    for (j = 0; j < nof_edges; j++) {
      et->setSelected(edge_ids[j], true, true);
    }
  }


  //=============================================
  //---Surface mesh drawing loop  (element edges)
  renderer->startDrawingMeshSurfaceEdges();

  for (i = 0; i < meshData->nofMeshElements; i++) {

    index = meshData->meshElementIds[i];

    bool as_selected = false;

    if (isLevel2) {
      as_selected = true;
    }

    elements->drawEdges(renderer, et, index, as_selected);
  }

  renderer->stopDrawingMeshSurfaceEdges();
}


// Draw border elements (edges (3D) or vertices (2D))
void
BodyElement::drawMeshBorder(Renderer* renderer, bool selected)
{
  if ( model->modelHasCadGeometry() ||
       meshData == NULL ||
       meshData->nofMeshBorderElements == 0
     ) {
    return;
  }

  MeshElementTable* et = model->getMeshBoundaryElementEdges();
  Point3* node_data = model->getMeshNodeData();

  int i, index;

  for (i = 0; i < meshData->nofMeshBorderElements; i++) {

    index = meshData->meshBorderElementIds[i];
    et->drawElement(renderer, index, selected);
  }
}


// Checks if bodyelement has a common boundary with an other bodyelement.
// Returns the matching type
matchType
BodyElement::findCommonBoundary(BodyElement* other, BodyElement*& common)
{
  matchType match_type = MATCH_NONE;
  common = NULL;

  if ( ptrGmtr == NULL)
    return MATCH_NONE;

  // Method to be called depends on geometry.
  ecif_geometryType tg1 = this->ptrGmtr->getType();
  ecif_geometryType tg2 = other->ptrGmtr->getType();

  // We have edges which are both linear:
  if ( (tg1 == ECIF_LINE || tg1 == ECIF_POLYLINE) &&
       (tg2 == ECIF_LINE || tg2 == ECIF_POLYLINE)
     ) {
    match_type = this->matchToLinear(other, common);
  }

  else if (tg1 == ECIF_NURBS && tg2 == ECIF_NURBS) {
    match_type = this->matchToNurbs(other, common);
  }

  return match_type;
}


int
BodyElement::findMeshBorder()
{
  if ( meshData == NULL || meshData->nofMeshElements == 0)
    return 0;

  MeshElementTable* bet = model->getMeshBoundaryElements();

  MeshElementTable* et = model->getMeshBoundaryElementEdges();

  int i, index;

  // Set first all subelements (edges or vertices) as non-selected
  for (i = 0; i < meshData->nofMeshElements; i++) {

    index = meshData->meshElementIds[i];
    meshElementCode bndr_code = bet->getElementCode(index);
    int nof_edges = MeshElementDesc[bndr_code][DESC_NOF_EDGES];

    const int* edge_ids = bet->getEdgeIds(index);

    for (int j = 0; j < nof_edges; j++) {
      et->setSelected(edge_ids[j], false, false);
    }
  }

  // Then set the sublements selected so
  // that the interior ones are left non-selected!
  for (i = 0; i < meshData->nofMeshElements; i++) {

    index = meshData->meshElementIds[i];
    meshElementCode bndr_code = bet->getElementCode(index);
    int nof_edges = MeshElementDesc[bndr_code][DESC_NOF_EDGES];

    const int* edge_ids = bet->getEdgeIds(index);

    for (int j = 0; j < nof_edges; j++) {
      et->setSelected(edge_ids[j], true, true);
    }
  }

  meshData->nofMeshBorderElements = 0;

  // Finally calculate number of selected sub_elements
  for (i = 0; i < meshData->nofMeshElements; i++) {

    index = meshData->meshElementIds[i];
    meshElementCode bndr_code = bet->getElementCode(index);
    int nof_edges = MeshElementDesc[bndr_code][DESC_NOF_EDGES];

    const int* edge_ids = bet->getEdgeIds(index);

    for (int j = 0; j < nof_edges; j++) {

      if ( et->getSelected(edge_ids[j]) ) {
        meshData->nofMeshBorderElements++;
      }
    }
  }

  // Create ids and direction indicator tables
  delete[] meshData->meshBorderElementIds;
  meshData->meshBorderElementIds = new int[meshData->nofMeshBorderElements];

  delete[] meshData->meshBorderElementDirs;
  meshData->meshBorderElementDirs = new short[meshData->nofMeshBorderElements];

  // Pick mesh border element ids
  int counter = 0;
  for (i = 0; i < meshData->nofMeshElements; i++) {

    index = meshData->meshElementIds[i];

    meshElementCode bndr_code = bet->getElementCode(index);
    int nof_edges = MeshElementDesc[bndr_code][DESC_NOF_EDGES];

    const int* edge_ids = bet->getEdgeIds(index);

    for (int j = 0; j < nof_edges; j++) {

      int edge_id = edge_ids[j];

      // If subelement is selected
      if ( et->getSelected(edge_ids[j]) ) {

        meshData->meshBorderElementIds[counter] = edge_id;
        meshData->meshBorderElementDirs[counter] = 1; // Currently always 1 !!!###!!!
        counter++;
      }
    }
  }

  return meshData->nofMeshBorderElements;
}


int
BodyElement::findSelectedMeshElement(Renderer* renderer,
                                     Point3& ray_start, Point3& ray_dir,
                                     Point3& isec_point,
                                     int& bndr_id,
                                     int& body1_id, int& layer1_id,
                                     int& body2_id, int& layer2_id)
{
  bndr_id = id;
  body1_id = NO_INDEX;
  body2_id = NO_INDEX;
  layer1_id = NO_INDEX;
  layer2_id = NO_INDEX;

  if ( meshData == NULL || meshData->nofMeshElements == 0 ||
       drawMode == DM_HIDDEN
     ) {
    return NO_INDEX;
  }

  int bd1_id = modelData->parent1Id;
  int bd2_id = modelData->parent2Id;

  Body* bd1 = model->getBodyById(bd1_id);
  Body* bd2 = model->getBodyById(bd2_id);

  int lr1_id = NO_INDEX;
  int lr2_id = NO_INDEX;

  if ( bd1 != NULL ) {
    lr1_id = bd1->getLayerId(modelData->parent1Layer);
  }

  if ( bd2 != NULL ) {
    lr2_id = bd2->getLayerId(modelData->parent2Layer);
  }

  // Filter: If both parent bodies are hidden
  if ( (bd1 == NULL || bd1->getDrawMode() == DM_HIDDEN) &&
       (bd2 == NULL || bd2->getDrawMode() == DM_HIDDEN)
     ) {
    return NO_INDEX;
  }

  // Filter: In 3D, if both parent bodies are visible, we cannot select this boundary!
  //
  // NOTE: If bd1_id==bd2_id this is an outer boundary and then no problems!
  //
  if ( ECIF_3D == model->getDimension() &&
       (bd1_id != bd2_id) &&
       (bd1 != NULL && bd1->getDrawMode() != DM_HIDDEN) &&
       (bd2 != NULL && bd2->getDrawMode() != DM_HIDDEN)
     ) {
    return NO_INDEX;
  }

  short direction;

  if ( bd1->getDrawMode() == DM_HIDDEN ) {
    body1_id = bd2_id;
    body2_id = bd1_id;
    layer1_id = lr2_id;
    layer2_id = lr1_id;
    direction = -1;

  } else {
    body1_id = bd1_id;
    body2_id = bd2_id;
    layer1_id = lr1_id;
    layer2_id = lr2_id;
    direction = 1;
  }

  // A Bem boundary (body) can be selected
  // form both sides!
  if ( isBemBoundary() ) {
    direction = 0;
  }

  flagName select_mode = model->getSelectionMode();

  MeshElementTable* elements = model->getMeshBoundaryElements();
  MeshElementTable* parents = model->getMeshBulkElements();

  double min_distance = 1.0e100;
  int isec_elem_index = NO_INDEX;

  for (int i = 0; i < meshData->nofMeshElements; i++) {

    int elem_index = meshData->meshElementIds[i];

    int nof_isections = elements->calcLineIntersections(elem_index, direction,
                                                        ray_start, ray_dir,
                                                        &isec_point);
#if 0
    if ( nof_isections > 0 ) {
      return elem_index;
    }
#endif

    if ( nof_isections == 0 ) {
      continue;
    }

    double distance = ray_start[2] - isec_point[2];

    distance = (distance < 0 ) ? -distance : distance;

    if ( distance < min_distance ) {
      min_distance = distance;
      isec_elem_index = elem_index;
    }

  } // for all mesh elements

  //return NO_INDEX;

  return isec_elem_index;
}


// Linearize original geometry and replace it with the polyline geometry
//
// NOTE: This should be used only  for internal purposes in the cases it is easier
// to handle geometry in discretized form (intersection calculations etc.)
// Geometry imported from Cad files should be stored in the original format!!!
//
void
BodyElement::formLinearizedGeometry()
{
  if ( modelData == NULL ) {
    return;
  }

  int nof_points_u;
  int nof_points_v;
  GcPoint** points;

  ptrGmtr->calcLinearizingPoints(nof_points_u, nof_points_v, points);

  if ( points == NULL ) return;

  int nof_points = 0;

  if ( nof_points_v > 0 ) {
    nof_points = nof_points_u * nof_points_v;

  } else {
    nof_points = nof_points_u;
  }

  // Replace the original geometry
  //
  delete ptrGmtr;

  ptrGmtr = new GcPolyLine(nof_points, points);

  delete[] points;
}


int
BodyElement::getBoundaryConditionId() const
{
  BodyElementGroup* beg = model->getBodyElementGroupById(elementGroupId);

  if ( beg == NULL ) return NO_INDEX;

  return beg->getBoundaryConditionId();
}


int
BodyElement::getBoundaryParameterId() const
{
  return boundaryParameterId;
}


void
BodyElement::getBoundaryPoints(int& count, BoundaryPoint**& points)
{
  if ( ptrGmtr == NULL ) {
    count = 0;
    points = NULL;

  } else {
    ptrGmtr->getBoundaryPoints(count, points);
  }
}


void
BodyElement::getDiscretizationData(int& nof_components, linDeltaType*& types, double*& valuesU, double*& valuesV, bool*& useFixedN)
{
  nof_components = 0;
  types = NULL;
  valuesU = NULL;
  valuesV = NULL;
  useFixedN = NULL;

  if ( ptrGmtr == NULL ) return;

  ptrGmtr->getDiscretizationData(nof_components, types, valuesU, valuesV, useFixedN);
}


const BodyElementGroup*
BodyElement::getElementGroup() const
{
  return model->getBodyElementGroupById(elementGroupId);
}


enum elementGroupType
BodyElement::getElementGroupType() const
{
  BodyElementGroup* beg = model->getBodyElementGroupById(elementGroupId);

  if ( beg == NULL ) return NONE_GROUP;

  return beg->getGroupType();
}


BodyElement*
BodyElement::getFirstSubElement()
{
  if ( modelData == NULL || modelData->nofSubElements == 0 )
    return NULL;

  return getSubElement(0);
}


const int*
BodyElement::getGridHIds()
{
  if ( gridHData != NULL )
    return gridHData->gridHIds;
  else
    return NULL;
}


const int*
BodyElement::getGridMeshIndices()
{
  if ( gridHData != NULL )
    return gridHData->meshIndices;
  else
    return NULL;
}


BodyElement*
BodyElement::getLastSubElement()
{
  if ( modelData == NULL || modelData->nofSubElements == 0 )
    return NULL;

  return getSubElement(modelData->nofSubElements - 1);
}


bool
BodyElement::getMeshDensityValues(int mesh_index, char& type, int& nof_values, double values[4])
{
  if ( gridHData == NULL && modelData == NULL )
    return false;

  if ( mesh_index < 0 || mesh_index >= model->getNofMeshes() )
    return false;

  // If directly set value available for the mesh
  // ============================================
  if ( gridHData != NULL &&
       gridHData->meshIndex == mesh_index &&
       gridHData->hValueType != ' ' &&
       gridHData->hValue > 0
     ) {
    type = gridHData->hValueType;
    nof_values = 1;
    values[0] = gridHData->hValue;
    return true;
  }

  // Check quadgrid mesh density
  // ===========================

  Body* bd1 = model->getBodyById(getParentId(1));
  Body* bd2 = model->getBodyById(getParentId(2));

  if ( bd1 != NULL ) {
    double n = bd1->getMeshQuadGridN(modelData->parent1Layer, mesh_index, id);
    if ( n > 0 ) {
      type = 'N';
      nof_values = 1;
      values[0] = n;
      return true;
    }
  }

  if ( bd2 != NULL ) {
    double n = bd2->getMeshQuadGridN(modelData->parent2Layer, mesh_index, id);
    if ( n > 0 ) {
      type = 'N';
      nof_values = 1;
      values[0] = n;
      return true;
    }
  }

  // Check normal mesh density
  // =========================

  if ( gridHData == NULL ) return false;

  type = ' ';
  nof_values = 0;
  values[0] = values[1] = values[2] = values[3] = 0;

  char dtype[2] = " ";
  double* dvalues = NULL;

  // Check if some parameter defined for the mesh
  for (int i = 0; i < gridHData->nofIds; i++) {

    if ( mesh_index == gridHData->meshIndices[i] ) {

      Parameter* param = model->getParameterById(ECIF_GRID_H, gridHData->gridHIds[i]);


      if ( param == NULL ) return false;

      // If value type given
      if ( param->getFieldValueBySifName("Mesh Density Type", 1, (char*)dtype) ) {

        // MeshH
        if ( dtype[0] == 'H' && param->getFieldValueBySifName("Mesh H", nof_values, dvalues) ) {
          type = 'H';

        // Mesh R
        } else if ( dtype[0] == 'R' && param->getFieldValueBySifName("Mesh R", nof_values, dvalues) ) {
          type = 'R';

        // Mesh N
        } else if ( dtype[0] == 'N' && param->getFieldValueBySifName("Mesh N", nof_values, dvalues) ) {
          type = 'N';
        }
      }

      break;
    }
  }

  if ( nof_values == 0 ) return false;

  for (int j = 0; j < nof_values; j++) {
    values[j] = dvalues[j];
  }
  delete[] dvalues;

  return true;
}


short
BodyElement::getMeshElementDir(int index)
{
  if ( meshData == NULL || meshData->nofMeshElements == 0 )
    return 0;

  if ( index < 0 || index >= meshData->nofMeshElements )
    return 0;

  return meshData->meshElementDirs[index];
}


int
BodyElement::getMeshElementId(int index)
{
  if ( meshData == NULL || index < 0 || index >= meshData->nofMeshElements)
    return NO_INDEX;

  return meshData->meshElementIds[index];
}


const char*
BodyElement::getName()
{
  if ( name == NULL || name[0] == '\0' ) {
    initName();
  }

  return name;
}


int
BodyElement::getNofBoundaryPoints() const
{
  if ( ptrGmtr == NULL ) {
    return 0;
  } else {
    return ptrGmtr->getNofBoundaryPoints();
  }
}


int
BodyElement::getNofComponents(bool only_emf_components) const
{
  if ( ptrGmtr != NULL )
    return ptrGmtr->getNofComponents(only_emf_components);
  else if ( meshData != NULL )
    return 1;
  else
    return 0;
}


int
BodyElement::getNofGridHIds() const
{
  if ( gridHData != NULL )
    return gridHData->nofIds;
  else
    return 0;
}


int
BodyElement::getNofMeshElements() const
{
  if ( meshData == NULL )
    return 0;
  else
    return meshData->nofMeshElements;
}


int
BodyElement::getNofMeshSelectedElements() const
{
  if ( meshData == NULL )
    return 0;
  else
    return meshData->nofMeshSelectedElements;
}


int
BodyElement::getNofSubElements() const
{
  if ( modelData == NULL )
    return 0;
  else
    return modelData->nofSubElements;
}


int
BodyElement::getNofVertices() const
{
  if ( modelData == NULL )
    return 0;

  return modelData->nofVertices;
}


// Method returns a list of bodyelements which together form
// the outer boundary-part of the element.
BodyElementList*
BodyElement::getOuterBoundary()
{
  if ( modelData == NULL || !(modelData->status & BE_OUTER_CHECKED) )
    return NULL;

  return findOuterBoundary();
}


int
BodyElement::getParentId(short parent) const
{
  if ( modelData == NULL )
    return NO_INDEX;

  if (parent == 1)
    return modelData->parent1Id;

  else if (parent == 2)
    return modelData->parent2Id;

  else
    return NO_INDEX;
}


int
BodyElement::getParentLayer(short parent) const
{
  if ( modelData == NULL )
    return NO_INDEX;

  if (parent == 1)
    return modelData->parent1Layer;

  else if (parent == 2)
    return modelData->parent2Layer;

  else
    return NO_INDEX;
}


int
BodyElement::getParentTag(short parent) const
{
  if ( modelData == NULL )
    return NO_INDEX;

  if (parent == 1)
    return modelData->parent1Tag;

  else if (parent == 2)
    return modelData->parent2Tag;

  else
    return NO_INDEX;
}


bool
BodyElement::getRangeVector(RangeVector rv)
{
  if ( ptrGmtr == NULL ) {
    return false;
  }

  ptrGmtr->getRangeVector(rv);

  return true;
}


void
BodyElement::getSelectedMeshElementIds(int& nof_ids, int*& ids)
{
  nof_ids = 0;

  if ( meshData == NULL || meshData->nofMeshElements == 0)
    return;

  MeshElementTable* bet = model->getMeshBoundaryElements();

  int counter, i;

  counter = 0;
  for (i = 0; i < meshData->nofMeshElements; i++) {

    int index = meshData->meshElementIds[i];

    if ( bet->selected[index] ) {
      counter++;
    }
  }

  nof_ids = counter;

  if (nof_ids == 0)
    return;

  ids = new int[nof_ids];

  counter = 0;
  for (i = 0; i < meshData->nofMeshElements; i++) {

    int index = meshData->meshElementIds[i];

    if ( bet->selected[index] ) {
      ids[counter++] = index;
    }
  }
}


int
BodyElement::getSubElementTag(int index)
{
  if ( modelData == NULL || modelData->nofSubElements == 0 )
    return NO_INDEX;

  if ( index < 0 || index >= modelData->nofSubElements )
    return NO_INDEX;

  return model->getModelObjectTagById((*modelData->subElementIds)[index]);
}


beStatus
BodyElement::getStatus()
{
  if ( modelData == NULL )
    return BE_NONE;

  else
    return modelData->status;
}



// Checks if bodyelement contains completely an other bodyelement.
bool
BodyElement::hasInside(BodyElement* other_element)
{
  if (ptrGmtr == NULL)
    return false;

  BoundBox* box1 = ptrGmtr->boundbox;
  BoundBox* box2 = other_element->ptrGmtr->boundbox;
  return box1->contains(box2);
}


/*
// Method tells if bodyelement has an outer boundary.
bool
BodyElement::hasOuterBoundary(int dir_in_body)
{
  if ( !(status & BE_OUTER_CHECKED) )
    checkOuterBoundary(dir_in_body);
  //-Return true if element is an outer boundary
  // or includes some as sub-elements.
  return (status & BE_INCLUDES_OUTER);
}
*/


bool
BodyElement::hasZeroVelocityBC()
{
  int bc_id = getBoundaryConditionId();

  if ( bc_id == NO_INDEX ) {
    return false;
  }

  BoundaryCondition* bc;

  bc = (BoundaryCondition*) model->getParameterById(ECIF_BOUNDARY_CONDITION,
                                                    bc_id);

  if ( bc != NULL &&  bc->hasZeroVelocity() ) {
    return true;

  } else {
    return false;
  }

}


void
BodyElement::init()
{
  gridHData = NULL;
  modelData = NULL;
  meshData = NULL;
  labelData = NULL;
  ptrGmtr = NULL;

  boundaryTag = NO_INDEX;
  elementGroupTag = NO_INDEX;
  elementGroupId = NO_INDEX;

  txBndrConditionId = NO_INDEX;

  boundaryParameterId = NO_INDEX;

  create_dyna_string(name, NULL);
  create_dyna_string(name_old, NULL);

  drawMode = DM_NORMAL;
  drawState = DS_NORMAL;
}


void
BodyElement::initLabelData()
{
  if ( ptrGmtr == NULL )
    return;

  if ( labelData == NULL ) {
    labelData = new LabelData;
    labelData->label = NULL;
  }

  // Make an id-string.
  ostrstream id_str;
  id_str << Tag() << ends;
  update_dyna_string(labelData->label, id_str.str());

  //---Get 'nice' position for the id-label
  ptrGmtr->getLabelPoint(labelData->position);
}


void
BodyElement::initName(char* be_name)
{
  // If name already set
  if (name != NULL && name[0] != '\0')
    return;

  // If name given as argument
  if ( be_name != NULL && be_name[0] !='\0' ) {
    update_dyna_string(name, be_name);
    update_dyna_string(name_old, name);
    return;
  }

  // Set default name
  strstream strm;

  objectType otp = getObjectType();

  switch (otp) {

  case OT_FACE:
    strm << "Boundary";
    break;
  case OT_EDGE:
    if ( model->getDimension() == ECIF_3D )
      strm << "Edge";
    else
      strm << "Boundary";
    break;
  case OT_VERTEX:
    strm << "Vertex";
    break;
  }

  if (tag != NO_INDEX) {
    strm << tag;
  }

  strm << ends;

  update_dyna_string(name, strm.str());
  update_dyna_string(name_old, name);

}


void
BodyElement::initClass(Model* mdl)
{
  BodyElement::last_boundaryTag = 0;

  // Init also sub-classe (for setting last_tag!)
  BodyElement1D::initClass(mdl);
  BodyElement2D::initClass(mdl);
  BodyElement3D::initClass(mdl);
}


bool
BodyElement::isBemBoundary()
{
  int bd_id = getParentId(1);

  Body* bd = model->getBodyById(bd_id);

  if ( bd != NULL && BEM_BODY == bd->getType() )
    return true;
  else
    return false;
}


// Nof of intersections with a line which is parallel to x-axis.
RayHit*
BodyElement::isectionsWithXdir(GcPoint* startp, bool& negative_on_left)
{
  negative_on_left = false;

  if (ptrGmtr == NULL ) {
    return NULL;

  } else {
    return ptrGmtr->isectionsWithXdir(startp, negative_on_left);
  }
}


bool
BodyElement::isInnerBoundary()
{
  if ( ( (getObjectType() == OT_EDGE && model->getDimension() == ECIF_2D) ||
         (getObjectType() == OT_FACE && model->getDimension() == ECIF_3D)
       ) &&
       modelData != NULL &&
       modelData->parent1Id != NO_INDEX &&
       modelData->parent2Id != NO_INDEX &&
       !isIntraLayerBoundary()
     ) {
    return true;

  } else {
    return false;
  }
}


void
BodyElement::markActiveMeshObjects(bool*& active_object_flags)
{
  // Subelements
  int nof_subs = getNofSubElements();
  int i;

  for (i = 0; i < nof_subs; i++) {

    BodyElement* se = getSubElement(i);

    if ( se == NULL || !se->isActive() )
      continue;

    active_object_flags[se->Id()] = true;

    se->markActiveMeshObjects(active_object_flags);
  }

  // Boundary points
  int bp_count;
  BoundaryPoint** bp_points;

  getBoundaryPoints(bp_count, bp_points);

  for (i = 0; i < bp_count; i++ ) {

    BoundaryPoint* bp = bp_points[i];

    bp->activeInMeshing = true;
  }

}


void
BodyElement::markActiveObjects()
{
  if ( modelData == NULL )
    return;

  int i;

  // Subelements
  for (i = 0; i < modelData->nofSubElements; i++) {
    int se_id = (*modelData->subElementIds)[i];
    model->markObjectActive(se_id);

    BodyElement* sbe = model->getBodyElementById(se_id);
    sbe->markActiveObjects();
  }

  // Vertices
  for (i = 0; i < modelData->nofVertices; i++) {
    int v_id = (*modelData->vertexIds)[i];
    model->markObjectActive(v_id);
  }

}


ostream&
BodyElement::output_emf(ostream& out, short indent_size, short indent_level,
                        bool output_geometry)
{
  char* QM = "\"";
  short is = indent_size;
  short il = indent_level;

  char* fn = NULL;

  objectType otp = getObjectType();

  //--Element type and tag
  if ( otp == OT_VERTEX ) {
    LibFront::output_scalar(out, is, il, EMF_VERTEX, NULL, tag);
    fn = "Vertex";
  } else if (otp == OT_EDGE ) {
    LibFront::output_scalar(out, is, il, EMF_EDGE, NULL, tag);
  } else if (otp == OT_FACE ) {
    LibFront::output_scalar(out, is, il, EMF_FACE, NULL, tag);
  } else {
    return out;
  }

  //--Boundary tag
  LibFront::output_scalar(out, is, 1 + il, "Boundary Tag", NULL, boundaryTag);

  //--Name
  if ( name != NULL && name[0] != '\0' ) {
    LibFront::output_scalar(out, is, 1 + il, "Name", NULL, name);
  }

  //--Subelements (edges or vertices)
  //
  // NOTE: No subelements for Edges. Their vertices are ouput via Components!!!
  //
  if ( modelData != NULL && modelData->nofSubElements > 0 ) {
    int count = 0;
    int* element_tags;

    if ( otp == OT_FACE ) {
      // Try edges
      modelData->getSubElementTags(model, OT_EDGE, count, element_tags);
      if ( count > 0 ) {
        LibFront::output_vector(out, is, 1 + il, EMF_EDGES, NULL, count, element_tags, false);
        delete[] element_tags;
      }

      // Try vertices
      modelData->getSubElementTags(model, OT_VERTEX, count, element_tags);
      if ( count > 0 ) {
        LibFront::output_vector(out, is, 1 + il, EMF_VERTICES, NULL, count, element_tags, false);
        delete[] element_tags;
      }

    // NOTE: No actual subelements for Edges.
    // Their vertices are ouput via Components!
    //
    } else if ( otp == OT_EDGE &&
                ( ptrGmtr == NULL ||
                  0 == ptrGmtr->getNofVertices()
                )
              ) {
      // Try (extra) vertices
      modelData->getSubElementTags(model, OT_VERTEX, count, element_tags);
      if ( count > 0 ) {
        LibFront::output_vector(out, is, 1 + il, EMF_EXTRA_VERTICES, NULL, count, element_tags, false);
        delete[] element_tags;
      }
    }
  }

  //--Geometry
  if ( output_geometry && ptrGmtr != NULL ) {
    ptrGmtr->output_emf(out, is, 1 + il);
  }

  // If boundary forms a single (logical) group), output
  // boundary condition and parameter info (because no group
  // is output!)
  //
  if ( EXPLICIT_GROUP != getElementGroupType() ) {

    //--Boundary Parameter Id
    int pid = getBoundaryParameterId();
    if ( pid != NO_INDEX ) {
      LibFront::output_scalar(out, is, 1 + il, EMF_BOUNDARY_PARAMETER, NULL, pid);
    }

    //--Boundary Condition Id
    int cid = getBoundaryConditionId();
    if ( cid != NO_INDEX ) {
      LibFront::output_scalar(out, is, 1 + il, EMF_BOUNDARY_CONDITION, NULL, cid);
    }
  }

  //--GridH
  if ( gridHData != NULL && gridHData->nofIds > 0 ) {
    int nof_ids = gridHData->nofIds;
    int* gids = gridHData->gridHIds;
    int* mids = gridHData->meshIndices;

    LibFront::output_vector(out, is, 1 + il, "Gridh Ids", NULL, nof_ids, gids, false);
    LibFront::output_vector(out, is, 1 + il, "Gridh Mesh Indices", NULL, nof_ids, mids, false);
  }

  return out;
}


// Output Boundary data into sif-file
//
// NOTE: Only boundary parameter data is output!
//
ostream&
BodyElement::output_sif(ostream& out, short indent_size, short indent_level)
{
  char* QM = "\"";
  short is = indent_size;
  short il = indent_level;

  // Boundary header and id
  LibFront::output_scalar(out, is, il, SIF_BOUNDARY, NULL, boundaryTag);

  // Name (in quotes)
  if ( model->getSolverKeywordTypeGiven(SIF_BOUNDARY, "Name") ) {
    LibFront::output_scalar(out, is, il + 1, "Name =", NULL, (char*)getName(), true);
  } else {
    LibFront::output_scalar(out, is, il + 1, "Name = String", NULL, (char*)getName(), true);
  }
  out << endl;

#if 0
  // Body1 tag
  if ( model->getSolverKeywordTypeGiven(SIF_BOUNDARY, "Body 1") ) {
    LibFront::output_scalar(out, is, il + 1, "Body 1 =", NULL, getParentTag(1));
  } else {
    LibFront::output_scalar(out, is, il + 1, "Body 1 = Integer", NULL, getParentTag(1));
  }

  // Body2 tag
  if ( model->getSolverKeywordTypeGiven(SIF_BOUNDARY, "Body 2") ) {
    LibFront::output_scalar(out, is, il + 1, "Body 2 =", NULL, getParentTag(2));
  } else {
    LibFront::output_scalar(out, is, il + 1, "Body 2 = Integer", NULL, getParentTag(2));
  }
#endif

  // Boundary parameter
  if ( boundaryParameterId != NO_INDEX ) {

    Parameter* p = model->getParameterById(ECIF_BOUNDARY_PARAMETER, boundaryParameterId);

    if ( p != NULL ) {
      SifOutputControl soc;
      soc.outputId = false;
      soc.outputName = false;
      soc.outputType = false;
      soc.outputAll = false;
      soc.sectionName = SIF_BOUNDARY;

      p->output_sif(out, is, il, soc);
    }
  }

  return out;
}


// Output boundary point tags (for Elmer Mesh input file)
//
ostream&
BodyElement::outputBoundaryPointTags(ostream& out, int indent_size)
{
  if ( modelData == NULL ) {
    return out;
  }

  int max_per_line = 20;

  //--Normal vertices
  if ( modelData->nofBoundaryPoints == 0 ) {

    int count;
    int* vertex_tags;

    modelData->getVertexTags(model, count, vertex_tags);

    for (int i = 0; i < count;  i++) {
      out << vertex_tags[i] << " ";
    }

    delete[] vertex_tags;

  //--Discretizating points
  } else {

    int counter = 0;

    for (int i = 0; i < modelData->nofBoundaryPoints;  i++) {

      BoundaryPoint* bp = modelData->boundaryPoints[i];

      if ( bp != NULL ) {
        out << bp->tag << " ";
        counter++;
      }

      // Continue from next line
      if ( counter >= max_per_line ) {
        out << endl;
        indent(out, indent_size);
        counter = 0;
      }
    }
  }

  return out;
}


// Output bodylements as a list of sub-element ids.
// Direction tells bodyelement's direction in the parent body.
//
ostream&
BodyElement::outputCoveringElementTags(ostream& out, int direction, int& nof_output_tags)
{
  nof_output_tags = 0;

  if ( modelData == NULL || modelData->nofSubElements == 0 )
    return out;

  if ( !(modelData->status & BE_DEVIDED) ) {
    out << direction * Tag();
    return out;
  }

  int se_tag;
  IdList* se_list = model->getCoveringElementList(this->tag);

  // Element is in positive loop-order (ccw)
  // Start reading from head forwards.
  if (direction == 1) {
    IdList::iterator pos = se_list->begin();
    while (pos != se_list->end()) {
      se_tag = *pos++;
      out << direction * se_tag << " ";
      nof_output_tags++;
    }
  }
  // Element is in negative loop-order (cw)
  // Start reading from tail backwards.
  else {
    IdList::reverse_iterator pos = se_list->rbegin();
    while (pos != se_list->rend()) {
      se_tag = *pos++;
      out << direction * se_tag << " ";
      nof_output_tags++;
    }
  }

  return out;
}


// Get point corresponding parametric values.
GcPoint*
BodyElement::param2Point(double u_p, double v_p)
{
  if ( ptrGmtr != NULL ) {
    return ptrGmtr->param2Point(u_p, v_p);
  } else {
    return NULL;
  }
}


// Get parametric values corresponding a point.
ParamPair*
BodyElement::point2Param(GcPoint* p)
{
  if ( ptrGmtr != NULL ) {
    return ptrGmtr->point2Param(p);
  } else {
    return NULL;
  }
}


// Increase mesh element tables
void
BodyElement::reallocateMeshElements()
{
  if ( meshData == NULL )
    return;

  meshData->maxNofMeshElements += 50;

  int size = meshData->maxNofMeshElements;

  int* element_ids = new int[size];
  short* element_dirs = new short[size];

  for (int i = 0; i < meshData->nofMeshElements; i++) {
    element_ids[i] = meshData->meshElementIds[i];
    element_dirs[i] = meshData->meshElementDirs[i];
  }

  delete[] meshData->meshElementIds;
  meshData->meshElementIds = element_ids;

  delete[] meshData->meshElementDirs;
  meshData->meshElementDirs = element_dirs;
}


// Increase pending mesh element info table
void
BodyElement::reallocatePendingMeshElementInfos()
{
  if ( meshData == NULL )
    return;

  meshData->maxNofPendingMeshElements += 50;

  int size = meshData->maxNofPendingMeshElements;

  int** infos = new int*[size];

  for (int i = 0; i < size; i++) {

    if ( i < meshData->nofPendingMeshElements ) {
      infos[i] = meshData->pendingMeshElementInfos[i];

    } else {
      infos[i] = NULL;
    }
  }

  delete[] meshData->pendingMeshElementInfos;
  meshData->pendingMeshElementInfos = infos;

}


// Delete given (with index!) mesh elements
// from the boundary
void
BodyElement::removeMeshElements(int nof_elements, int* element_indices)
{
  if ( meshData == NULL )
    return;

  if (nof_elements == 0 || meshData->nofMeshElements == 0)
    return;

  // construct a flag array and init it
  bool* delete_flags = new bool[meshData->nofMeshElements];

  for (int i = 0; i < meshData->nofMeshElements; i++) {
    delete_flags[i] = false;
  }

  // Set flag active for the indices to be removed
  for (int j = 0; j < nof_elements; j++) {
    delete_flags[element_indices[j]] = true;
  }

  deleteMeshElements(delete_flags);

  delete[] delete_flags;
}


// Delete selected (in the display meaning!) mesh elements
// from the boundary
void
BodyElement::removeSelectedMeshElements(int* removed_elements_storage)
{
  if ( meshData == NULL || meshData->nofMeshSelectedElements == 0)
    return;

  bool* delete_flags = new bool[meshData->nofMeshElements];

  MeshElementTable* elements = model->getMeshBoundaryElements();

  int remove_counter = 0;
  for (int i = 0; i < meshData->nofMeshElements; i++) {

    int index = meshData->meshElementIds[i];

    //--Selected mesh element
    if ( elements->selected[index] ) {
      delete_flags[i] = true;

      // If we should store removed indices (for split redo!)
      if (removed_elements_storage != NULL)
        removed_elements_storage[remove_counter++] = i;
    }

    //--Non selected mesh element
    else {
      delete_flags[i] = false;
    }
  }

  deleteMeshElements(delete_flags);

  meshData->nofMeshSelectedElements = 0;

  delete[] delete_flags;
}


beStatus
BodyElement::removeStatus(beStatus old_status)
{
  if ( modelData == NULL )
    return BE_NONE;

  else
    return modelData->status &= ~old_status;
}


// Copy data from backup
// NOTE: "pointer copy" DO NOT delete
// meshDataBackup arries, just repoint
// them and then make "unvisible"
void
BodyElement::restoreMeshDataBackup()
{
  if ( meshData == NULL || meshData->backupData == NULL)
    return;

  //--Mesh elements
  meshData->maxNofMeshElements = meshData->backupData->nofMeshElements;
  meshData->nofMeshElements = meshData->backupData->nofMeshElements;

  delete[] meshData->meshElementIds;
  meshData->meshElementIds = meshData->backupData->meshElementIds;

  delete[] meshData->meshElementDirs;
  meshData->meshElementDirs = meshData->backupData->meshElementDirs;

  //--Mesh border elements
  meshData->nofMeshBorderElements = meshData->backupData->nofMeshBorderElements;

  delete[] meshData->meshBorderElementIds;
  meshData->meshBorderElementIds = meshData->backupData->meshBorderElementIds;

  delete[] meshData->meshBorderElementDirs;
  meshData->meshBorderElementDirs = meshData->backupData->meshBorderElementDirs;

  //--Delete backup struct (but NOT the arries inside!!!)
  delete meshData->backupData;
  meshData->backupData = NULL;
}


void
BodyElement::restoreName()
{
  create_dyna_string(name, name_old);
}


void
BodyElement::setBoundaryParameterId(int pid)
{
  boundaryParameterId = pid;
}


void
BodyElement::setBoundaryPointsMeshDensityValue(int mesh_index)
{
  char dtype;
  int nof_dvalues;
  double dvalues[4];

  BodyElement* v;
  Point3 p1;
  Point3 p2;
  double delta;
  double tot_length, cum_length;
  int i;

  if ( modelData->nofBoundaryPoints == 0 ) {
    return;
  }

  BoundaryPoint* bp;

  // Reset old values
  for (i = 1; i < modelData->nofBoundaryPoints; i++) {
    bp = modelData->boundaryPoints[i];
    bp->meshDensityType = ' ';
    bp->meshDensityValue = 0;

    // Delete density values also from any vertex
    if ( bp->vertexTag != NO_INDEX ) {
      v = model->getVertexByTag(bp->vertexTag);
      if ( v != NULL ) {
        v->setGridHValues(NO_INDEX, ' ', 0, true);
      }
    }
  }

  getMeshDensityValues(mesh_index, dtype, nof_dvalues, dvalues);

  if ( nof_dvalues < 2 || dtype == 'N' ) return;

  delta = dvalues[1] - dvalues[0];

  BoundaryPoint* bp1;
  BoundaryPoint* bp2;

  // Calc total boundary length
  tot_length = 0.0;
  for (i = 1; i < modelData->nofBoundaryPoints; i++) {
    bp1 = modelData->boundaryPoints[i-1];
    bp2 = modelData->boundaryPoints[i];
    bp1->point->getPoint(p1);
    bp2->point->getPoint(p2);
    tot_length += dist3(p1, p2);
  }

  // Set values for the first point
  bp = modelData->boundaryPoints[0];
  bp->meshDensityType = dtype;
  bp->meshDensityValue = dvalues[0];

  if ( bp->vertexTag != NO_INDEX ) {
    v = model->getVertexByTag(bp->vertexTag);
    if ( v != NULL ) {
      v->setGridHValues(mesh_index, dtype, dvalues[0]);
    }
  }

  int last;

  // If closed, do not set value again for the first point!
  if (modelData->isClosedU) {
    last = modelData->nofBoundaryPoints - 1;
  } else {
    last = modelData->nofBoundaryPoints;
  }

  cum_length = 0.0;
  for (i = 1; i < last; i++) {
    bp1 = modelData->boundaryPoints[i-1];
    bp2 = modelData->boundaryPoints[i];
    bp1->point->getPoint(p1);
    bp2->point->getPoint(p2);
    cum_length += dist3(p1, p2);

    double dvalue = dvalues[0] + (cum_length / tot_length) * delta;

    bp2->meshDensityType = dtype;
    bp2->meshDensityValue = dvalue;

    if ( bp2->vertexTag != NO_INDEX ) {
      v = model->getVertexByTag(bp2->vertexTag);
      if ( v != NULL ) {
        v->setGridHValues(mesh_index, dtype, dvalue);
      }
    }
  }

}


void
BodyElement::setElementGroupId(int gr_id)
{
  elementGroupId = gr_id;
}


void
BodyElement::setElementGroupTag(int gr_tag)
{
  elementGroupTag = gr_tag;
}


void
BodyElement::setBoundaryTag(int tag)
{
  boundaryTag = tag;
}


void
BodyElement::setDiscretizationData(int nof_components, linDeltaType* types, double* valuesU, double* valuesV, bool* useFixedN)
{
  if ( ptrGmtr == NULL ) return;

  ptrGmtr->setDiscretizationData(nof_components, types, valuesU, valuesV, useFixedN);
}


void
BodyElement::setMeshElementDir(int index, short direction)
{
  if ( meshData == NULL || meshData->nofMeshElements == 0 )
    return;

  if ( index < 0 || index >= meshData->nofMeshElements)
    return;

  meshData->meshElementDirs[index] = direction;
}


void
BodyElement::setGridHData(int nof_ids, int* grid_h_ids, int* mesh_indices, bool force)
{
  if ( gridHData == NULL ) {

    if ( nof_ids > 0 ) {
      gridHData = new BodyElementGridHData();

    } else {
      return;
    }
  }

  gridHData->setData(nof_ids, grid_h_ids, mesh_indices);
}


void
BodyElement::setGridHValues(int mesh_index, char value_type, double value, bool force)
{
  if ( gridHData == NULL ) {
    gridHData = new BodyElementGridHData();
  }

  gridHData->setData(mesh_index, value_type, value, force);

#if 0
  if ( value_type == 'N' ) {
    ptrGmtr->setDeltaU(LIN_DELTA_N, value);
  } else if ( value_type == 'H' ) {
    ptrGmtr->setDeltaU(LIN_DELTA_H, value);
  } else if ( value_type == 'R' ) {
    ptrGmtr->setDeltaU(LIN_DELTA_H, value * model->getMeshH(mesh_index));
  }

  ptrGmtr->calcBoundaryPoints();
#endif

}


void
BodyElement::setIsIntraLayerBoundary(bool value)
{
  if ( modelData != NULL ) {
    modelData->isIntraLayerBoundary = value;
  }
}


void
BodyElement::setName(char* new_name)
{
  update_dyna_string(name, new_name);
}


void
BodyElement::setParentId(short parent_nbr, int parent_id)
{
  if ( modelData == NULL )
    return;

  if ( parent_nbr == 1 ) {
    modelData->parent1Id = parent_id;
  } else  {
    modelData->parent2Id = parent_id;
  }
}


void
BodyElement::setParentIds(int parent1_id, int parent2_id)
{
  if ( modelData == NULL ) return;

  modelData->parent1Id = parent1_id;
  modelData->parent2Id = parent2_id;
}


void
BodyElement::setParentLayer(short parent_nbr, int layer)
{
  if ( modelData == NULL ) return;

  if ( parent_nbr == 1 ) {
    modelData->parent1Layer = layer;
  } else  {
    modelData->parent2Layer = layer;
  }
}


void
BodyElement::setParentLayers(int parent1_layer, int parent2_layer)
{
  if ( modelData == NULL ) return;

  modelData->parent1Layer = parent1_layer;
  modelData->parent2Layer = parent2_layer;
}


void
BodyElement::setParentTag(short parent_nbr, int parent_tag)
{
  if ( modelData == NULL )
    return;

  if ( parent_nbr == 1 ) {
    modelData->parent1Tag = parent_tag;
  } else  {
    modelData->parent2Tag = parent_tag;
  }
}


void
BodyElement::setParentTags(int parent1_tag, int parent2_tag)
{
  if ( modelData == NULL ) return;

  modelData->parent1Tag = parent1_tag;
  modelData->parent2Tag = parent2_tag;
}


beStatus
BodyElement::setStatus(beStatus new_status)
{
  if ( modelData == NULL )
    return BE_NONE;
  else
    return modelData->status = new_status;
}


void
BodyElement::storeName()
{
  create_dyna_string(name_old, name);
}


// Add new data into element
//
// NOTE: Only non-geometry related data is update here!
//
void
BodyElement::update(ecif_Element_X& tx)
{
  if ( tx.name != NULL ) {
    update_dyna_string(name, tx.name);
  }

  if ( tx.bndr_cond_id != NO_INDEX ) {
    txBndrConditionId = tx.bndr_cond_id;
  }

  if ( tx.bndr_param_id != NO_INDEX ) {
    boundaryParameterId = tx.bndr_param_id;
  }

  if ( tx.bndr_tag != NO_INDEX ) {
    boundaryTag = tx.bndr_tag;
    checkLastBoundaryTag();
  }

  if ( tx.bndr_group_tag != NO_INDEX ) {
    elementGroupTag = tx.bndr_group_tag;
  }

  if ( tx.nof_gridh_ids > 0 ) {
    gridHData = new BodyElementGridHData(tx.nof_gridh_ids, tx.gridh_ids, tx.gridh_mesh_indices);
  }

}


bool
BodyElement::updateGeometry()
{
  IdList vertex_tags;

  if ( ptrGmtr != NULL ) {
    if ( !ptrGmtr->updateGeometry(tag, vertex_tags) ) return false;
  }

  update();

  return true;
}


void
BodyElement::update()
{
  initLabelData();

  if ( modelData != NULL && ptrGmtr != NULL ) {

    modelData->isClosedU = ptrGmtr->isClosedU();
    modelData->isClosedV = ptrGmtr->isClosedV();
  }
}
