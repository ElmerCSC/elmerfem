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
Module:     ecif_modelMeshManager.cpp
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
#include "ecif_bodyElement.h"
#include "ecif_bodyElement1D.h"
#include "ecif_bodyElement2D.h"
#include "ecif_bodyElement3D.h"
#include "ecif_boundaryCondition.h"
#include "ecif_boundbox.h"
#include "ecif_control.h"
#include "ecif_input.h"
#include "ecif_mesh.h"
#include "ecif_modelMeshManager.h"
#include "ecif_model_aux.h"
#include "ecif_timer.h"
#include "ecif_userinterface.h"

void pressAnyKey();
extern char read_buffer[];

Control* ModelMeshManager::theControlCenter = NULL;
Model* ModelMeshManager::model = NULL;


ModelMeshManager::ModelMeshManager()
{
  // NOTE: These are mesh-manager's 'own' data
  // They are deallocated when input elements are removed!
  //
  meshInputElements = NULL;
  meshInputElementsExt2Int = NULL;

  meshInputElementsMaxExtId = 0;

  nofMeshInputBulkElements = 0;
  nofMeshInputBoundaryElements = 0;
  nofMeshInputEdgeElements = 0;
  nofMeshInputVertexElements = 0;

  nofMeshInputElements = 0;
  nofAllocMeshInputElements = 0;

  // NOTE: These pointers are set by the Model
  // Do NOT delete them in this class!!!
  meshData = NULL;
  meshInfo = NULL;
  meshBox = NULL;

  modelData = NULL;
  modelInfo = NULL;

}


ModelMeshManager::~ModelMeshManager()
{
}


int
ModelMeshManager::addMeshBoundaryElement(
                              int bndr_tag, meshElementCode el_code,
                              int ext_parent1_id, int ext_parent2_id,
                              const int* ext_node_ids,
                              bool& is_added_to_bndr)
{
  BodyElement* be = model->getBoundaryByTag(bndr_tag);

  // Create boundary if all necesary info is available
  if (be == NULL && meshData->bulkElements != NULL) {

    // Convert parent element ids external --> internal
    int elem1_id = getMeshBulkElementIdExt2Int(ext_parent1_id);
    int elem2_id = getMeshBulkElementIdExt2Int(ext_parent2_id);

    int body1_tag = meshData->bulkElements->getParentId(elem1_id, 0);
    int body2_tag = meshData->bulkElements->getParentId(elem2_id, 0);

    int body1_lr = 0;
    int body2_lr = 0;

    if (body1_tag != NO_INDEX) {
      be = model->createBodyElement(bndr_tag, body1_tag, body1_lr, body2_tag, body2_lr, 0);
    }
  }
 
  return addMeshBoundaryElement(be, el_code,
                                ext_parent1_id, ext_parent2_id,
                                ext_node_ids,
                                is_added_to_bndr);
}


// Add element with external node ids
int
ModelMeshManager::addMeshBoundaryElement(
                              BodyElement* boundary, meshElementCode el_code,
                              int ext_parent1_id, int ext_parent2_id,
                              const int* ext_node_ids,
                              bool& is_added_to_bndr)
{
  static int parent_pair[2];
  static int int_node_ids[MAX_NOF_NODES];

  // Convert parent element ids external --> internal
  parent_pair[0] = getMeshBulkElementIdExt2Int(ext_parent1_id);
  parent_pair[1] = getMeshBulkElementIdExt2Int(ext_parent2_id);

  meshInfo->nofUsedElementTypes[el_code]++;
  int nof_nodes = MeshElementDesc[el_code][DESC_NOF_NODES];

  // Convert node ids external --> internal
  for (short i = 0; i < nof_nodes; i++) {
      int_node_ids[i] = getMeshNodeIdExt2Int(ext_node_ids[i]);
  }

  int elem_id = meshData->lastBoundaryElementIndex++;

  meshData->boundaryElements->createTableEntry(elem_id, el_code, parent_pair, nof_nodes, int_node_ids);

  short direction = 1;

  is_added_to_bndr = false;

  // Add element to the boundary
  if (boundary != NULL) {
    meshData->boundaryElementBoundaryIds[elem_id] = boundary->Id();
    boundary->addMeshElement(elem_id, direction);
    is_added_to_bndr = true;
  }

  return elem_id;
}


// Add element with internal node ids
int
ModelMeshManager::addMeshBoundaryElement(BodyElement* boundary, meshElementCode elem_code,
                              int int_parent1_id, int int_parent2_id,
                              int nof_nodes, const int* int_node_ids,
                              bool& is_added_to_bndr)
{
  static int parent_pair[2];

  parent_pair[0] = int_parent1_id;
  parent_pair[1] = int_parent2_id;

  UserInterface* gui = theControlCenter->getGui();

  meshInfo->nofUsedElementTypes[elem_code]++;

  int elem_id = meshData->lastBoundaryElementIndex++;

  //-Check element table overflow!!!
	if ( elem_id  >= meshData->boundaryElements->nofElements ) {
    gui->errMsg(1, "Fatal Error: Mesh boundary element internal id error");
		return NO_INDEX;
	}

  meshData->boundaryElements->createTableEntry(elem_id, elem_code, parent_pair, nof_nodes, int_node_ids);

  short direction = 1;

  is_added_to_bndr = false;

  // Add element to the boundary and store boundary id in the ids-array
  if (boundary != NULL) {
    meshData->boundaryElementBoundaryIds[elem_id] = boundary->Id();
    boundary->addMeshElement(elem_id, direction);
    is_added_to_bndr = true;
  }

  return elem_id;
}




Rc
ModelMeshManager::addMeshBulkElement(int ext_id, int material_id,
                          meshElementCode el_code, const int* ext_node_ids,
                          int nof_ngbrs, const int* ext_ngbr_ids)
{
  static Timer timer;

  if (meshInfo->nofBulkElements == 0)
    timer.start();

  UserInterface* gui = theControlCenter->getGui();

  static int parent_pair_ids[2];
  static int int_node_ids[MAX_NOF_NODES];

  int int_id = meshData->lastBulkElementIndex++;

  //-Check element table overflow!!!
	if ( int_id  >= meshInfo->nofAllocElements ) {
    gui->errMsg(1, "Fatal Error: Mesh bulk element internal id error");
		return ECIF_MESH_ELEMENT_INDEX_ERROR;
	}

  //-Check element-dictionary table overflow!!!
	if ( ext_id < 0 || ext_id  > meshInfo->maxExternalElementId ) {
    gui->errMsg(1, "Fatal Error: Mesh bulk element external id error");
		return ECIF_MESH_ELEMENT_INDEX_ERROR;
	}

  //-Check material id and update counter
  if ( material_id > MAX_NOF_BODIES ) {
    return ECIF_INCORRECT_MATERIAL_ID;
  }

  //-Mark the existence of the material
  // NOTE: we will change values later in this table
  // into body ids (when we have read all elements)
  // NOTE: NO_INDEX (-1) <--> body/material not yet known
  if (material_id >= 0) {
    meshData->bodyExt2Int[material_id] = 1;
  }

	//-Update external <--> internal element-id dictionaries
	meshData->bulkElementExt2Int[ext_id] = int_id;
	meshData->bulkElementInt2Ext[int_id] = ext_id;
	meshInfo->nofBulkElements++;

  meshInfo->nofUsedElementTypes[el_code]++;
  int nof_nodes = MeshElementDesc[el_code][DESC_NOF_NODES];

  // Convert node ids external --> internal
  for (short i = 0; i < nof_nodes; i++) {
      int_node_ids[i] = getMeshNodeIdExt2Int(ext_node_ids[i]);
  }

	// Allocate space for the element's data: code, material, + nodes
  // NOTE: nof_nodes for the element can be concluded from the element code,
  // there is no need to store nof-nodes.
  // NOTE: material-id will be changed later into body-id ! 
  // NOTE: external node ids will be changed into internal ids later !
  parent_pair_ids[0] = material_id;
  parent_pair_ids[1] = NO_INDEX;

  meshData->bulkElements->createTableEntry(int_id, el_code, parent_pair_ids,
                                           nof_nodes, int_node_ids);

  //-If neighbour data is available 
  // External element ids in this phase!
  if (nof_ngbrs > 0) {
    meshData->bulkElements->neighborIds[int_id] = new int[nof_ngbrs];
    for (int j = 0; j < nof_ngbrs; j++)
      meshData->bulkElements->neighborIds[int_id][j] = ext_ngbr_ids[j];
  }

  gui->showProgressMsg(timer, 2000,
                       meshInfo->nofBulkElements,
                       meshInfo->nofAllocElements,
                       "Reading ", " elements ");

	return ECIF_OK;
}


Rc
ModelMeshManager::addMeshNode(int int_id, int ext_id, Point3& point)
{
  static Timer timer;

  if (meshInfo->nofNodes == 0) {
    timer.start();
  }

  UserInterface* gui = theControlCenter->getGui();

	//-Check node table overflow!!!
	if ( int_id < 0 || int_id  >= meshInfo->nofAllocNodes ) {
    gui->errMsg(1, "Fatal Error: Mesh: node internal id index error");
		return ECIF_MESH_NODE_INDEX_ERROR;
	}
	//-Check node-dictionary-table overflow!!!
	if ( ext_id < 0 || ext_id  > meshInfo->maxExternalNodeId ) {
    gui->errMsg(1, "Fatal Error: Mesh: node external id index error");
		return ECIF_MESH_NODE_INDEX_ERROR;
	}

	//-Update external <--> internal node-id dictionaries
	meshData->nodeExt2Int[ext_id] = int_id;
	meshData->nodeInt2Ext[int_id] = ext_id;
  meshInfo->nofNodes++;

  meshData->nofNodes++;

  // NOTE: Mesh input units are applied here!!!
  //
  // meshInfo->unit = 0.001;
  double unit = Model::getMeshInputUnit();

  meshData->nodes[int_id][0] = unit * point[0];
  meshData->nodes[int_id][1] = unit * point[1];
  meshData->nodes[int_id][2] = unit * point[2];

  double x = meshData->nodes[int_id][0];
  double y = meshData->nodes[int_id][1];
  double z = meshData->nodes[int_id][2];

  if (x < meshInfo->minX)
    meshInfo->minX = x;
  if (x > meshInfo->maxX)
    meshInfo->maxX = x;
  if (y < meshInfo->minY)
    meshInfo->minY = y;
  if (y > meshInfo->maxY)
    meshInfo->maxY = y;
  if (z < meshInfo->minZ)
    meshInfo->minZ = z;
  if (z > meshInfo->maxZ)
    meshInfo->maxZ = z;


  gui->showProgressMsg(timer, 2000,
                       meshInfo->nofNodes,
                       meshInfo->nofAllocNodes,
                       "Reading ", " nodes ");
 
	return ECIF_OK;
}


// Add a temporary element (bulk or boundary )
//
// NOTE: This can be used only when the total number elements is
// known and we do not want to read the mesh file twice
// to know separately the nofBulkElements and nofBoundaryElements
//
Rc
ModelMeshManager::addMeshInputElement(int el_type,
                                      int ext_elem_id, int parent_tag, int ext_parent_tag,
                                      int* ext_node_ids)
{
  static Timer timer;

  if (nofMeshInputElements == 0)
    timer.start();

  UserInterface* gui = theControlCenter->getGui();

  //-Check table overflow!!!
	if ( nofMeshInputElements >= nofAllocMeshInputElements ) {
    gui->errMsg(1, "Fatal Error: Mesh: nof temporary elements incorrect");
    removeMeshInputElements();
		return ECIF_MESH_ELEMENT_INDEX_ERROR;
	}

  // Find nof-nodes
  meshElementCode el_code = model->convertElementType(el_type);
  int nof_nodes = MeshElementDesc[el_code][DESC_NOF_NODES];

  MeshInputElement& te = meshInputElements[nofMeshInputElements];
  
  bool is_bulk, is_bndr, is_edge;
  checkMeshElementType(el_type, is_bulk, is_bndr, is_edge);
  
  if ( is_bulk) 
    te.type = 'U';
  else if ( is_bndr ) 
    te.type = 'B';
  else if ( is_edge ) 
    te.type = 'E';
  else
    te.type = 'V';
  
  te.elementCode = el_code;
  te.extElementId = ext_elem_id;
  te.elementId = NO_INDEX; // Internal id is here unknown
  te.parentTag = parent_tag;
  te.extParentTag = ext_parent_tag;
  te.extNodeIds = new int[nof_nodes];

	// Copy node ids
  for (int i = 0; i < nof_nodes; i++) {
    te.extNodeIds[i] = ext_node_ids[i];
  }
  
  meshInputElementsExt2Int[ext_elem_id] = nofMeshInputElements;

  nofMeshInputElements++;

  if ( te.type == 'U' ) {
    nofMeshInputBulkElements++;
  } else if ( te.type == 'B' ) {
    nofMeshInputBoundaryElements++;
  } else if ( te.type == 'E' ) {
    nofMeshInputEdgeElements++;
  } else if ( te.type == 'V' ) {
    nofMeshInputVertexElements++;
  }


  gui->showProgressMsg(timer, 2000,
                       nofMeshInputElements,
                       nofAllocMeshInputElements,
                       "Reading ", " elements ");

	return ECIF_OK;
}


void
ModelMeshManager::allocateMeshBodies(int nof_bodies)
{
  // Currently nothing else is done!
  meshInfo->nofBodies = nof_bodies;
}


// Allocates a table of mesh boundary elements
// NOTE: memory for each item is allocated when data is
// actually read, because we dont know yet the nof-node per each
// element.
void
ModelMeshManager::allocateMeshBoundaryElements(int nof_elements)
{
  UserInterface* gui = theControlCenter->getGui();

  strstream strm;
  strm << "Nof boundary elements: " << nof_elements << ends;
  gui->showMsg(strm.str(), 1);

  meshInfo->nofBoundaryElements = nof_elements;

  // 2D
  if (modelInfo->dimension == ECIF_2D) {
    meshData->boundaryElements = new MeshEdgeElementTable(nof_elements, true);
  // 3D
  } else {
    meshData->boundaryElements = new MeshFaceElementTable(nof_elements, true);
  }

  // Allocate data arries which are needed for boundary elements
  meshData->boundaryElements->parentIds = new int*[meshInfo->nofBoundaryElements];
  meshData->boundaryElements->actionLevels = new char[meshInfo->nofBoundaryElements];
  meshData->boundaryElements->centers = new Point3[meshInfo->nofBoundaryElements];
  meshData->boundaryElements->dirInParents = new short*[meshInfo->nofBoundaryElements];
  meshData->boundaryElements->edgeIds = new int*[meshInfo->nofBoundaryElements];
  meshData->boundaryElements->normalDistances = new double[meshInfo->nofBoundaryElements];
  meshData->boundaryElements->normals = new Point3[meshInfo->nofBoundaryElements];
  meshData->boundaryElements->nodeNormals = new Point3*[meshInfo->nofBoundaryElements];
  meshData->boundaryElements->rSquares = new double[meshInfo->nofBoundaryElements];

  meshData->boundaryElementBoundaryIds = new int[meshInfo->nofBoundaryElements];

  // Init data
  for (int i = 0; i < meshInfo->nofBoundaryElements; i++) {
    meshData->boundaryElements->parentIds[i] = new int[2];
    meshData->boundaryElements->actionLevels[i] = 0;
    meshData->boundaryElements->centers[i][0] = 0.0;
    meshData->boundaryElements->centers[i][1] = 0.0;
    meshData->boundaryElements->centers[i][2] = 0.0;

    meshData->boundaryElements->dirInParents[i] = new short[2];
    meshData->boundaryElements->dirInParents[i][0] = 0;
    meshData->boundaryElements->dirInParents[i][1] = 0;

    meshData->boundaryElements->edgeIds[i] = NULL;

    meshData->boundaryElements->normalDistances[i] = 0.0;

    meshData->boundaryElements->normals[i][0] = 0.0;
    meshData->boundaryElements->normals[i][1] = 0.0;
    meshData->boundaryElements->normals[i][2] = 0.0;

    meshData->boundaryElements->nodeNormals[i] = NULL;

    meshData->boundaryElements->rSquares[i] = 0.0;

    meshData->boundaryElementBoundaryIds[i] = NO_INDEX;

  }
}


// Allocates a table of mesh bulk elements
// NOTE: memory for each item is allocated when data is
// actually read, because we dont know yet the nof-nodes for each
// element.
void
ModelMeshManager::allocateMeshBulkElements(int nof_elements, int max_external_element_id)
{
  UserInterface* gui = theControlCenter->getGui();

  strstream strm;
  strm << "Nof elements: " << nof_elements << ends;
  gui->showMsg(strm.str());

  meshData->bulkElements = new MeshBulkElementTable(nof_elements);

  meshInfo->nofAllocElements = nof_elements;

  if (max_external_element_id == NO_INDEX)
    return;

  int i;

  meshData->bulkElementExt2Int = new int[max_external_element_id + 1];
  for (i = 0; i <= max_external_element_id; i++)
    meshData->bulkElementExt2Int[i] = NO_INDEX;

  meshData->bulkElementInt2Ext = new int[nof_elements];
  for (i = 0; i < nof_elements; i++)
    meshData->bulkElementInt2Ext[i] = NO_INDEX;

  meshInfo->maxExternalElementId = max_external_element_id;
}


// Allocates space for the mesh-node table
void
ModelMeshManager::allocateMeshNodes(int nof_nodes, int max_external_node_id)
{
  UserInterface* gui = theControlCenter->getGui();

  strstream strm;
  strm << "Nof nodes: " << nof_nodes << ends;
  gui->showMsg(strm.str());

  meshData->nodes = new Point3[nof_nodes];

  meshInfo->nofAllocNodes = nof_nodes;
  meshData->nofAllocNodes = nof_nodes;

  if (max_external_node_id == NO_INDEX)
    return;

  meshData->nodeExt2Int = new int[max_external_node_id + 1];
  meshData->nodeInt2Ext = new int[nof_nodes];
  
  meshData->nofBoundaryNodeParents = new int[nof_nodes];
  meshData->boundaryNodeParentIds = new int*[nof_nodes];

  for (int i = 0; i < nof_nodes; i++) {
    meshData->boundaryNodeParentIds[i] = NULL;
  }

  meshInfo->maxExternalNodeId = max_external_node_id;

}


// Allocates space for the mesh (temporary) input element table
void
ModelMeshManager::allocateMeshInputElements(int nof_elements, int max_ext_elem_id)
{
  UserInterface* gui = theControlCenter->getGui();

  strstream strm;
  strm << "Nof bulk and boundary elements: " << nof_elements << ends;
  gui->showMsg(strm.str());

  meshInputElements = new MeshInputElement[nof_elements];

  meshInputElementsExt2Int = new int[1 + max_ext_elem_id];
  meshInputElementsMaxExtId = max_ext_elem_id;

  for (int i = 0; i < max_ext_elem_id; i++) {
    meshInputElementsExt2Int[i] = NO_INDEX;
  }

  nofMeshInputElements = 0;
  nofAllocMeshInputElements = nof_elements;
}

// Calculate average normals for mesh boundary nodes
//
void 
ModelMeshManager::calcMeshBoundaryNodeNormals()
{
  int i, j, k, n;

  MeshElementTable* bet = meshData->boundaryElements;
  
  // Init element node normals
  for (i = 0; i < meshInfo->nofBoundaryElements; i++) {

    meshElementCode elem_code = bet->getElementCode(i);
    int nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];

    delete[] bet->nodeNormals[i];
    bet->nodeNormals[i] = new Point3[nof_nodes];

    for (j = 0; j < nof_nodes; j++) {
      for (k = 0; k < 3; k++) {
        bet->nodeNormals[i][j][k] = bet->normals[i][k];
      }
    }
  }

  // Calc average normals by connected elements for each node in the element
  //
  for (i = 0; i < meshInfo->nofBoundaryElements; i++) {

    meshElementCode elem_code = bet->getElementCode(i);
    int nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];

    const int* bndr_node_ids = bet->nodeIds[i];

    for (j = 0; j < nof_nodes; j++) {
      
      int nid = bndr_node_ids[j];

      for (n = 0; n < meshData->nofBoundaryNodeParents[nid]; n++) {

        int nbr_id = meshData->boundaryNodeParentIds[nid][n];
        
        // Skip element itself!
        if ( nbr_id == i ) continue;

        // Check that the parent element does not differ too much
        // by orientation (ie. if we have a 'real' corner, skip)
        // Limit is pi/8 = 30.0 degrees
        double dp = dot3(bet->normals[nbr_id], bet->normals[i]);

        if ( acos(fabs(dp)) > (PI / 6) ) {
          continue;
        }
        
        // If element normals happen to be dirrently oriented
        int sign = 1;

        if ( dp < 0.0 ) {
          sign = -1;
        }
        
        // Add parent's contribution
        for (k = 0; k < 3; k++) {
          bet->nodeNormals[i][j][k] += sign * bet->normals[nbr_id][k];
        }

      } // for each parent for the node

      normalize(bet->nodeNormals[i][j]);

    } // for each node in elment
  } // for each element
}


bool
ModelMeshManager::checkMeshElement(meshElementCode element_code, int* node_ids, bool swap_to_ccw)
{
  bool result; 
  switch (element_code) {
  case MEC_303:
  case MEC_304:
    result = checkMeshElement303(node_ids, swap_to_ccw);
    break;
  case MEC_306:
  case MEC_307:
    result = checkMeshElement306(node_ids, swap_to_ccw);
    break;
  case MEC_404:
    result = checkMeshElement404(node_ids, swap_to_ccw);
    break;
  case MEC_408:
  case MEC_409:
    result = checkMeshElement408(node_ids, swap_to_ccw);
    break;
  case MEC_504:
    result = checkMeshElement504(node_ids, swap_to_ccw);
    break;
  case MEC_508:
    result = checkMeshElement508(node_ids, swap_to_ccw);
    break;
  case MEC_510:
    result = checkMeshElement510(node_ids, swap_to_ccw);
    break;
  case MEC_808:
    result = checkMeshElement808(node_ids, swap_to_ccw);
    break;
  case MEC_820:
    result = checkMeshElement820(node_ids, swap_to_ccw);
    break;
  default:
    break;
  }

  return result;
}


// Note this is called from input-objects, because
// they know when the mesh data must be checked!
//
// NOTE: Not in use!!!
//
bool
ModelMeshManager::checkMeshElements(MeshElementTable* table, bool swap_to_ccw)
{
  bool result = true;

  for (int i = 0; i < table->NofElements(); i++) {

    meshElementCode code = table->getElementCode(i);
    int* nodes = table->nodeIds[i];
    result = checkMeshElement(code, nodes, swap_to_ccw);
  }

  return result;
}


bool
ModelMeshManager::checkMeshElement303(int* nodes, bool swap_to_ccw)
{
  double* a = (double*) &meshData->nodes[nodes[0]];
	double*	b = (double*) &meshData->nodes[nodes[1]];
	double*	c = (double*) &meshData->nodes[nodes[2]];


  bool is_ccw = isCcwTriangle(a, b, c);

	//If triangle is not in CCW-order, swap nodes 1,2
	if (swap_to_ccw && !is_ccw) {
		int tmp;
    tmp = nodes[1]; nodes[1] = nodes[2]; nodes[2] = tmp;
    is_ccw = true;
  }

	return is_ccw;
}


bool
ModelMeshManager::checkMeshElement306(int* nodes, bool swap_to_ccw)
{
  double* a = (double*) &meshData->nodes[nodes[0]];
	double*	b = (double*) &meshData->nodes[nodes[1]];
	double*	c = (double*) &meshData->nodes[nodes[2]];


  bool is_ccw = isCcwTriangle(a, b, c);

	//If triangle is not in CCW-order, swap nodes 1,2 and 3,5
	if (swap_to_ccw && !is_ccw) {
		int tmp;
    tmp = nodes[1]; nodes[1] = nodes[2]; nodes[2] = tmp;
    tmp = nodes[3]; nodes[3] = nodes[5]; nodes[5] = tmp;
    is_ccw = true;
  }

	return is_ccw;
}



bool
ModelMeshManager::checkMeshElement404(int* nodes, bool swap_to_ccw)
{
  double* a = (double*) &meshData->nodes[nodes[0]];
	double*	b = (double*) &meshData->nodes[nodes[1]];
	double*	c = (double*) &meshData->nodes[nodes[2]];


  bool is_ccw = isCcwTriangle(a, b, c);

	//If rectangle is not in CCW-order, swap nodes 1,3
	if (swap_to_ccw && !is_ccw) {
		int tmp;
    tmp = nodes[1]; nodes[1] = nodes[3]; nodes[3] = tmp;
    is_ccw = true;
  }

	return is_ccw;
}


bool
ModelMeshManager::checkMeshElement408(int* nodes, bool swap_to_ccw)
{
  double* a = (double*) &meshData->nodes[nodes[0]];
	double*	b = (double*) &meshData->nodes[nodes[2]];
	double*	c = (double*) &meshData->nodes[nodes[4]];


  bool is_ccw = isCcwTriangle(a, b, c);

	//If rectangle is not in CCW-order, swap nodes: 1,3 4,7 and 5,6.
	if (swap_to_ccw && !is_ccw) {
    int tmp;
		tmp = nodes[1]; nodes[1] = nodes[3]; nodes[3] = tmp;
		tmp = nodes[4]; nodes[4] = nodes[7]; nodes[7] = tmp;
		tmp = nodes[5]; nodes[5] = nodes[6]; nodes[6] = tmp;
    is_ccw = true;
  }

	return is_ccw;
}


bool
ModelMeshManager::checkMeshElement504(int* nodes, bool swap_to_ccw)
{
  Point3& a = meshData->nodes[nodes[0]];
	Point3&	b = meshData->nodes[nodes[1]];
	Point3&	c = meshData->nodes[nodes[2]];
	Point3&	d = meshData->nodes[nodes[3]];

	//Calculate tetrahedron volume from determinant
	// O'Rourke, Computational geometry, p. 26)
	// |a0 a1 a2 1|
	// |b0 b1 b2 1|
	// |c0 c1 c2 1|
	// |d0 d1 d2 1|
	double vol = 0 
		-b[0]*c[1]*d[2] +a[0]*c[1]*d[2] +b[1]*c[0]*d[2] -a[1]*c[0]*d[2] -a[0]*b[1]*d[2] 
		+a[1]*b[0]*d[2] +b[0]*c[2]*d[1] -a[0]*c[2]*d[1] -b[2]*c[0]*d[1] +a[2]*c[0]*d[1] 
		+a[0]*b[2]*d[1] -a[2]*b[0]*d[1] -b[1]*c[2]*d[0] +a[1]*c[2]*d[0] +b[2]*c[1]*d[0] 
		-a[2]*c[1]*d[0] -a[1]*b[2]*d[0] +a[2]*b[1]*d[0] +a[0]*b[1]*c[2] -a[1]*b[0]*c[2] 
		-a[0]*b[2]*c[1] +a[2]*b[0]*c[1] +a[1]*b[2]*c[0] -a[2]*b[1]*c[0];

	//Final orientation of the bottom-triangle is concluded from sign 
	//of the volume. We want the bottom-triangle to be in CCW-order looked
	//from outside (ie. CW-ordered when looked from the top node).
	bool is_ccw = (vol < 0.0);

	//If bottom-triangle is not in CCW-order, swap nodes 1,2
	if (swap_to_ccw && !is_ccw) {
		int tmp;
    tmp = nodes[1]; nodes[1] = nodes[2]; nodes[2] = tmp;
    is_ccw = true;
  }

	return is_ccw;
}


bool
ModelMeshManager::checkMeshElement508(int* nodes, bool swap_to_ccw)
{
  Point3& a = meshData->nodes[nodes[0]];
	Point3&	b = meshData->nodes[nodes[1]];
	Point3&	c = meshData->nodes[nodes[2]];
	Point3&	d = meshData->nodes[nodes[3]];

	//Calculate tetrahedron volume from determinant
	// O'Rourke, Computational geometry, p. 26)
	// |a0 a1 a2 1|
	// |b0 b1 b2 1|
	// |c0 c1 c2 1|
	// |d0 d1 d2 1|
	double vol = 0 
		-b[0]*c[1]*d[2] +a[0]*c[1]*d[2] +b[1]*c[0]*d[2] -a[1]*c[0]*d[2] -a[0]*b[1]*d[2] 
		+a[1]*b[0]*d[2] +b[0]*c[2]*d[1] -a[0]*c[2]*d[1] -b[2]*c[0]*d[1] +a[2]*c[0]*d[1] 
		+a[0]*b[2]*d[1] -a[2]*b[0]*d[1] -b[1]*c[2]*d[0] +a[1]*c[2]*d[0] +b[2]*c[1]*d[0] 
		-a[2]*c[1]*d[0] -a[1]*b[2]*d[0] +a[2]*b[1]*d[0] +a[0]*b[1]*c[2] -a[1]*b[0]*c[2] 
		-a[0]*b[2]*c[1] +a[2]*b[0]*c[1] +a[1]*b[2]*c[0] -a[2]*b[1]*c[0];

	//Final orientation of the bottom-triangle is concluded from sign 
	//of the volume. We want the bottom-triangle to be in CCW-order looked
	//from outside (ie. CW-ordered when looked from the top node).
	bool is_ccw = (vol < 0.0);

	//If bottom-triangle is not in CCW-order, swap nodes 1,2  5,7  is this ok???
	if (swap_to_ccw && !is_ccw) {
		int tmp;
    tmp = nodes[1]; nodes[1] = nodes[2]; nodes[2] = tmp;
    tmp = nodes[5]; nodes[5] = nodes[7]; nodes[7] = tmp;
    is_ccw = true;
  }

	return is_ccw;
}


bool
ModelMeshManager::checkMeshElement510(int* nodes, bool swap_to_ccw)
{
  Point3& a = meshData->nodes[nodes[0]];
	Point3&	b = meshData->nodes[nodes[1]];
	Point3&	c = meshData->nodes[nodes[2]];
	Point3&	d = meshData->nodes[nodes[3]];

	//Calculate tetrahedron volume from determinant
	// O'Rourke, Computational geometry, p. 26)
	// |a0 a1 a2 1|
	// |b0 b1 b2 1|
	// |c0 c1 c2 1|
	// |d0 d1 d2 1|
	double vol = 0 
		-b[0]*c[1]*d[2] +a[0]*c[1]*d[2] +b[1]*c[0]*d[2] -a[1]*c[0]*d[2] -a[0]*b[1]*d[2] 
		+a[1]*b[0]*d[2] +b[0]*c[2]*d[1] -a[0]*c[2]*d[1] -b[2]*c[0]*d[1] +a[2]*c[0]*d[1] 
		+a[0]*b[2]*d[1] -a[2]*b[0]*d[1] -b[1]*c[2]*d[0] +a[1]*c[2]*d[0] +b[2]*c[1]*d[0] 
		-a[2]*c[1]*d[0] -a[1]*b[2]*d[0] +a[2]*b[1]*d[0] +a[0]*b[1]*c[2] -a[1]*b[0]*c[2] 
		-a[0]*b[2]*c[1] +a[2]*b[0]*c[1] +a[1]*b[2]*c[0] -a[2]*b[1]*c[0];

	//Final orientation of the bottom-triangle is concluded from sign 
	//of the volume. We want the bottom-triangle to be in CCW-order looked
	//from outside (ie. CW-ordered when looked from the top node).
	bool is_ccw = (vol < 0.0);

	//If bottom-triangle is not in CCW-order, swap nodes 1,2  4,5, 6,7  is this ok???
	if (swap_to_ccw && !is_ccw) {
		int tmp;
    tmp = nodes[1]; nodes[1] = nodes[2]; nodes[2] = tmp;
    tmp = nodes[4]; nodes[4] = nodes[5]; nodes[5] = tmp;
    tmp = nodes[6]; nodes[6] = nodes[7]; nodes[7] = tmp;
    is_ccw = true;
  }

	return is_ccw;
}


bool
ModelMeshManager::checkMeshElement808(int* nodes, bool swap_to_ccw)
{
	return true;
}


bool
ModelMeshManager::checkMeshElement820(int* nodes, bool swap_to_ccw)
{
	return true;
}


void
ModelMeshManager::checkMeshElementType(int elem_type,
                                       bool& is_bulk_element,
                                       bool& is_bndr_element,
                                       bool& is_edge_element)
{
  is_bulk_element = false;
  is_bndr_element = false;
  is_edge_element = false;

  if (modelInfo->dimension == ECIF_ND) return;

  // 2D model
  if ( modelInfo->dimension == ECIF_2D ) {
    if ( elem_type < 200 ) {
      is_edge_element = true;
    } else if ( elem_type < 300 ) {
      is_bndr_element = true;
    } else {
      is_bulk_element = true;
    }

  // 3D model
  } else {
    if ( elem_type < 200 ) {
      return;
    } else if ( elem_type < 300 ) {
      is_edge_element = true;
    } else if ( elem_type < 500 ) {
      is_bndr_element = true;
    } else {
      is_bulk_element = true;
    }
  }
}


// Check if some of the boundary elements are not installed into any boundaries
// These elements will create Bem boundaries !!!
//
void
ModelMeshManager::checkMeshInputBoundaryElements()
{
  Rc rc;
  int i;
  int elem_index;
  Body* bd;
  BodyElement* be;
  
  IdTable bndrTagsTbl;
  IdTable nofElementsTbl;

  // Count nof elements per boundary
  //
  for (i = 0; i < nofMeshInputElements; i++) {

    MeshInputElement& te = meshInputElements[i];

    // Skip if not a boundary element or boundary given or element already added
    if ( te.type != 'B' ||
         te.parentTag != NO_INDEX ||
         te.isAdded == '1'
       ) {
      continue;
    }
    
    int count = nofElementsTbl[te.extParentTag];
    nofElementsTbl[te.extParentTag] = ++count;
  }

  // Create boundaries and add elements
  //
  // NOTE: Create also parent bodies if needed (relevant for BEM models!)
  //
  for (i = 0; i < nofMeshInputElements; i++) {

    MeshInputElement& te = meshInputElements[i];

    // Skip if not a boundary element
    if ( te.type != 'B' ) {
      continue;
    }

    be = NULL;

    // Boundary tag given
    // ------------------
    //
    if ( te.parentTag != NO_INDEX ) {
      
      // Try to get the boundary
      //
      if ( modelInfo->dimension == ECIF_2D ) {
        be = model->getBodyElementByTag(OT_EDGE, te.parentTag);
      } else {
        be = model->getBodyElementByTag(OT_FACE, te.parentTag);
      }

      // Check that parent body exists, create if necessary
      // 
      if ( be != NULL ) {

        int bd_tg = be->getParentTag(1);
        Body* bd = model->getBodyByTag(bd_tg);

        if ( bd == NULL ) {
          if ( modelInfo->dimension == ECIF_2D ) {
            bd = new Body2D();
          } else {
            bd = new Body3D();
          }
    
          bd->setType(BEM_BODY);

          be->setParentId(1, bd->Id());
          be->setParentTag(1, bd->Tag());
          be->setParentLayer(1,0);

          model->addBody(bd);
          model->addBodyElement(be, true);
        }
      }

    // Only external boundary tag given
    // --------------------------------
    //
    } else {

      // Check if boundary already created
      int bndr_tg = bndrTagsTbl[te.extParentTag];

      // Create a new boundary and the parent body for the boundary
      if ( bndr_tg == 0 ) {
      
        int count = nofElementsTbl[te.extParentTag];

        if ( modelInfo->dimension == ECIF_2D ) {
          bd = new Body2D();
          be = new BodyElement2D(bd->Tag(), NO_INDEX, count);
        } else {
          bd = new Body3D();
          be = new BodyElement3D(bd->Tag(), NO_INDEX, count);
        }

        bd->setType(BEM_BODY);
        be->setParentLayer(1,0);

        model->addBody(bd);
        model->addBodyElement(be, true);

        bndrTagsTbl[te.extParentTag] = bd->Tag();

      } else {
        if ( modelInfo->dimension == ECIF_2D ) {
          be = model->getBodyElementByTag(OT_EDGE, bndr_tg);
        } else {
          be = model->getBodyElementByTag(OT_FACE, bndr_tg);
        }
      }
    }

    // Add mesh element into boundary
    //
    if ( be != NULL &&
         te.elementId != NO_INDEX  &&
         te.isAdded == '0'
       ) {
      be->addMeshElement(te.elementId, 1);
    }
  }
}


// Classify "corner" bulk elements, currently:
// No-slip bc --> zero velocity element
void
ModelMeshManager::classifyMeshCornerElements()
{
  if ( meshInfo->nofCornerElements == 0 ) {
    meshInfo->nofZeroVelocityElements = 0;
    return;
  }

  MeshElementTable* bt = meshData->bulkElements;

  UserInterface* gui = theControlCenter->getGui();

  short dimension =  modelInfo->dimension;

  int zero_velocity_counter = 0;
 
  int index = 0;
  while (true) {

    MeshCornerElement* mce = getMeshCornerElement(index++);

    if (mce==NULL) break;

    // If already corrected
    if (mce->corrected) {
      continue;
    }

    mce->hasZeroVelocity = false;

    // Try to negate this
    bool has_zero_velocity = true;

    // If elements' parent is not a flow body!
    int obj_id = bt->parentIds[mce->elementId][0];
    Body* body = (Body*)model->getModelObjectById(obj_id);

    if ( !body->canHaveZeroVelocityElements() ) {
      continue;
    }

    meshElementCode bulk_code = bt->getElementCode(mce->elementId);
    int nof_bulk_nodes = MeshElementDesc[bulk_code][DESC_NOF_NODES];

    for (int i = 0; i < nof_bulk_nodes; i++) {

      int nd_id = mce->nodeIds[i];
      int be_id = mce->boundaryIds[i];

      // Check first if node is also a vertex
      BodyElement* be = model->getVertexByNodeId(nd_id);

      // Otherwise use the boundary id
      if ( be == NULL ) {
        be = model->getBoundaryById(be_id);
      }

      if ( be == NULL || !be->hasZeroVelocityBC() ) {
        has_zero_velocity = false;
        break;
      }

    }

    // If all boundary elements really had a zero-velocity
    // boundary condition
    // Test!
    //if (true || has_zero_velocity) {
    if (has_zero_velocity) {
      
      mce->hasZeroVelocity = true;

      zero_velocity_counter++;

#if 0
      strstream strm;
      strm << "***WARNING Corner element "
           << 1 + mce->elementId
           << " has zero-velocity boundary condition!" 
           << ends;

      gui->showMsg(strm.str());
#endif
    }
  }

  // Update info in GUI
  gui->updateMeshZeroVelocityElements(zero_velocity_counter);

  meshInfo->nofZeroVelocityElements = zero_velocity_counter;

}


// External ids --> internal ids
void
ModelMeshManager::convertMeshBulkElementIdsExt2Int()
{
  // If no bulk elements (pure boundary model)
  //
  if ( meshInfo->nofBulkElements == 0 ) return;

  // If neighbor ids are not yet set, don't try 
  // to convert them!
  bool convert_neighbor_ids = true;

  if ( meshData->bulkElements->neighborIds[0][0] == UNSET_INDEX ) {
    convert_neighbor_ids = false;
  }

  short j;
  for (int i = 0; i < meshInfo->nofBulkElements; i++) {

    meshElementCode elm_code = meshData->bulkElements->getElementCode(i);

    // Body tags (material ids) --> body object ids in the element data
    int material_id = meshData->bulkElements->parentIds[i][0];
    int body_tag = meshData->bodyExt2Int[material_id];
    Body* body = model->getBodyByTag(body_tag);
    meshData->bulkElements->parentIds[i][0] = body->Id();

#if 0
    // NOTE: Nodes are already converted via addMeshElement... methods
    // External node ids --> internal node ids
    int elm_nof_nodes = MeshElementDesc[elm_code][DESC_NOF_NODES];
    for (j = 0; j < elm_nof_nodes; j++) {
      int ext_node_id = meshData->bulkElements->node_ids[i][j];
      int int_node_id = meshData->nodeExt2Int[ext_node_id];
      meshData->bulkElements->node_ids[i][j] = int_node_id;
    }
#endif

    if ( !convert_neighbor_ids ||
         meshData->bulkElements->neighborIds[i] == NULL
       )
      continue;

    // External neighbor element ids --> internal element ids
    int elm_nof_ngbrs = MeshElementDesc[elm_code][DESC_NOF_BNDR_ELEMS];
    for (j = 0; j < elm_nof_ngbrs; j++) {
      int ext_elem_id = meshData->bulkElements->neighborIds[i][j];

      int int_elem_id;
      if (ext_elem_id != NO_INDEX)
        int_elem_id = meshData->bulkElementExt2Int[ext_elem_id];
      else
        int_elem_id = NO_INDEX;

      meshData->bulkElements->neighborIds[i][j] = int_elem_id;
    }

  } // for all mesh elements

}


int
ModelMeshManager::convertElementCode(meshElementCode element_code)
{
  switch (element_code) {

  case MEC_101: //Vertex
    return 101; 
  case MEC_202: //Linear Beam
    return 202; 
  case MEC_203: //Parabolic Beam
    return 203; 
  case MEC_303: //Linear Triangle
    return 303;
  case MEC_304: //Linear Triangle, one in middle (not any more in Solver!)
    return 304;
  case MEC_306: //Parabolic Triangle
    return 306;
  case MEC_307: //Parabolic Triangle, one in middle
    return 307;
  case MEC_404: //Linear Quadrilateral
    return 404;
  case MEC_408: //Parabolic Quadrilateral
    return 408;
  case MEC_409: //Parabolic Quadrilateral, one in middle
    return 409;
  case MEC_504: //Linear Tetra
    return 504;
  case MEC_508: //'SemiLinear' Tetra (center nodes at faces)
    return 508;
  case MEC_510: //Parabolic Tetra
    return 510;
  case MEC_605: //Linear Prism
    return 605;
  case MEC_613: //Parabolic Prism
    return 613;
  case MEC_706: //Linear Wedge
    return 706;
  case MEC_715: //Parabolic Wedge
    return 715;
  case MEC_808: //Linear Brick
    return 808;
  case MEC_820: //Parabolic Brick
    return 820;
  case MEC_827: //Parabolic Brick, one in middle of faces, one in the middle
    return 827;

  default:
    strstream strm;
    strm << "Unknown element code: "
         << element_code
         << ". Stopping reading the mesh!"
         << endl;
    theControlCenter->getGui()->showMsg(strm.str());
    theControlCenter->setBreakValue(MESH_INPUT, true);
    return 0;
  }
}



meshElementCode
ModelMeshManager::convertElementType(int element_type)
{
  switch (element_type) {

  case 101:   //Vertex
    return MEC_101; 
  case 202:   //Linear Beam
    return MEC_202; 
  case 203:   //Parabolic Beam
    return MEC_203; 
  case 303: //Linear Triangle
    return MEC_303;
  case 304: //Linear Triangle, one in middle (not any more in Solver!)
    return MEC_304;
  case 306: //Parabolic Triangle
    return MEC_306;
  case 307: //Parabolic Triangle, one in middle
    return MEC_307;
  case 404: //Linear Quadrilateral
    return MEC_404;
  case 408: //Parabolic Quadrilateral
    return MEC_408;
  case 409: //Parabolic Quadrilateral, one in middle
    return MEC_409;
  case 504: //Linear Tetra
    return MEC_504;
  case 508: //'SemiLinear' Tetra (center nodes at faces)
    return MEC_508;
  case 510: //Parabolic Tetra
    return MEC_510;
  case 605: //Linear Prism
    return MEC_605;
  case 613: //Parabolic Prism
    return MEC_613;
  case 706: //Linear Wedge
    return MEC_706;
  case 715: //Parabolic Wedge
    return MEC_715;
  case 808: //Linear Brick
    return MEC_808;
  case 820: //Parabolic Brick
    return MEC_820;
  case 827: //Parabolic Brick, one in middle of faces, one in the middle
    return MEC_827;

  default:
    strstream strm;
    strm << "Unknown element type: "
         << element_type
         << ". Stopping reading the mesh!"
         << endl;
    theControlCenter->getGui()->showMsg(strm.str());
    theControlCenter->setBreakValue(MESH_INPUT, true);
    return MEC_000;
  }
}


void
ModelMeshManager::correctMeshZeroVelocityElements()
{
  // Nothing to do!
  if ( meshInfo->nofZeroVelocityElements == 0 ) {
    return;
  }

  // Message for GUI
  // ===============
  UserInterface* gui = theControlCenter->getGui();

  strstream strm;
  strm << "Correcting " << meshInfo->nofZeroVelocityElements
       << " zero-velocity element(s)!" << ends;

  gui->showMsg(strm.str());

  // Start processing
  // ================
  int nof_zv_elements = meshInfo->nofZeroVelocityElements;

  MeshCornerElement* mce;
  int i,j, index;

  int* body_ids = new int[nof_zv_elements];

  for (i = 0; i < nof_zv_elements; i++) {
    body_ids[i] = NO_INDEX;
  }

  int old_bulk_counter = 0; // Nof removed (splitted) bulk corners
  int new_bulk_counter = 0; // Nof new bulk elements from splitting
  int new_edge_counter = 0; // Nof new edges from splitting
  int new_node_counter = 0; // Nof new nodes from splitting

  // Create bulk remove flags table
  // We need this to mark original elemetns as NOT
  // to be copied!
  int nof_elements = meshInfo->nofBulkElements;

  bool* bulk_remove_flags = new bool[nof_elements];
  for (i = 0; i < nof_elements; i++) {
    bulk_remove_flags[i] = false;
  }

  // Count nof new elements
  // ======================
  index = 0;
  while (true) {
    
    mce = getMeshCornerElement(index++);

    if (mce==NULL) break;

    if (!mce->hasZeroVelocity) {
      continue;
    }

    // Store corner element body id (if not yet stored)
    for (i = 0; i < nof_zv_elements; i++) {

      if ( mce->bodyId == body_ids[i] ) {
        break;
      }

      if ( body_ids[i] == NO_INDEX ) {
        body_ids[i] = mce->bodyId;
        break;
      }
    }

    // We will remove the original elements from
    // bodies!
    bulk_remove_flags[mce->elementId] = true;

    if ( mce->elementCode == MEC_303 ) {
      old_bulk_counter += 1;
      new_bulk_counter += 3;
      new_edge_counter += 3;
    }

    if ( mce->elementCode == MEC_504 ) {
      old_bulk_counter += 1;
      new_bulk_counter += 4;
      new_edge_counter += 4;
    }

    new_node_counter += 1;

  }

  // Reallocate mesh tables
  // ======================
  MeshElementTable* table;
  int old_size;
  int new_size;

  // Bulk related
  // ------------
  table = meshData->bulkElements;
  new_size = table->NofElements() + new_bulk_counter;
  table->resize(new_size, NULL);
  meshInfo->nofAllocElements = new_size;
  meshInfo->nofSplittedElements = old_bulk_counter;

  //-Reallocate Ext2Int ids
  old_size = 1 + meshInfo->maxExternalElementId;
  new_size = 1 + meshInfo->maxExternalElementId + new_bulk_counter;
  meshData->resizeIdTable(meshData->bulkElementExt2Int, old_size, new_size, NULL); 

  //-Reallocate Int2Ext ids
  old_size = meshInfo->nofBulkElements;
  new_size = meshInfo->nofBulkElements + new_bulk_counter;
  meshData->resizeIdTable(meshData->bulkElementInt2Ext, old_size, new_size, NULL); 

  //-Reallocate bulk edge rendered flags
  old_size = meshInfo->nofBulkEdges;
  new_size = meshInfo->nofBulkEdges + new_edge_counter;
  meshData->resizeFlagTable(meshData->bulkEdgeRendered, old_size, new_size, NULL); 

  //-Remove splitted elements fromt the bodies
  for (i = 0; i < nof_zv_elements; i++) {

    if ( body_ids[i] == NO_INDEX )
      break;

    Body* body = model->getBodyById(body_ids[i]);

    if ( body != NULL ) {
      body->removeMeshElements(bulk_remove_flags);
    }
  }

  // Node related
  // ------------
  old_size = meshInfo->nofNodes;
  new_size = new_node_counter + table->NofElements();
  meshData->resizeNodeTable(old_size, new_size, NULL);
  setMeshNodes(); // NOTE: Rememeber to do this!!!
  meshInfo->nofAllocNodes = new_size;

  //-Resize Ext2Int ids
  old_size = 1 + meshInfo->maxExternalNodeId;
  new_size = 1 + meshInfo->maxExternalNodeId + new_node_counter;
  meshData->resizeIdTable(meshData->nodeExt2Int, old_size, new_size, NULL); 

  //-Resize Int2Ext ids
  old_size = meshInfo->nofNodes;
  new_size = meshInfo->nofNodes + new_node_counter;
  meshData->resizeIdTable(meshData->nodeInt2Ext, old_size, new_size, NULL); 

  // Set correct starting values for bulk element counters before
  // splitting the elements
  // ======================
  meshData->bulkElements->nofElements = meshInfo->nofBulkElements;
  int nof_old_elements = meshInfo->nofBulkElements;

  // Split corner elements
  // =====================
  index = 0;
  while (true) {

    mce = getMeshCornerElement(index++);
    
    if (mce==NULL) break;

    if (mce->hasZeroVelocity) {
    
      splitMeshCornerElement(mce);

      mce->corrected = true;
      mce->splitted = true;
      meshInfo->nofZeroVelocityElements--;
    }
  }

  // Bulk renumbering table
  // ======================
  delete[] meshData->bulkRenumbering;
  meshData->bulkRenumbering = new int[meshInfo->nofBulkElements];

  index = 0;
  for (i = 0; i < meshInfo->nofBulkElements; i++) {

    if ( !meshData->bulkElements->splitted[i] ) {
      meshData->bulkRenumbering[i] = index++;

    } else {
      meshData->bulkRenumbering[i] = NO_INDEX;
    }
  }

  // Add new elements to the bodies
  // ==============================
  // NOTE: At most nof_zero_velocity_elements different
  // bodies were encountered!
  int* nof_new_bbulks = new int[nof_zv_elements];
  int** new_bbulk_ids = new int*[nof_zv_elements];

  for (i = 0; i < nof_zv_elements; i++) {
    nof_new_bbulks[i] = 0;
    new_bbulk_ids[i] = NULL;
  }

  //-Count nof new mesh bulk elements per body
  for (i = 0; i < new_bulk_counter; i++) {

    int body_id = meshData->bulkElements->parentIds[nof_old_elements + i][0];

    for (j = 0; j < nof_zv_elements; j++) {

      if ( body_id == body_ids[j] ) {
        nof_new_bbulks[j]++;
        break;
      }
    }
  }

  //-Allocate ids tables per body to store new bulk ids
  for (i = 0; i < nof_zv_elements; i++) {

    if ( body_ids[i] == NO_INDEX )
      break;

    new_bbulk_ids[i] = new int[nof_new_bbulks[i]];

    // Reset counter for the next step!
    nof_new_bbulks[i] = 0;
  }

  //-Copy bulk ids to tables per body
  for (i = 0; i < new_bulk_counter; i++) {

    int body_id = meshData->bulkElements->parentIds[nof_old_elements + i][0];
    int elem_id = nof_old_elements + i;

    for (j = 0; j < nof_zv_elements; j++) {

      if ( body_id == body_ids[j] ) {
        new_bbulk_ids[j][nof_new_bbulks[j]] = elem_id;
        nof_new_bbulks[j]++;
        break;
      }
    }
  }

  //-Add elements to the bodies
  for (i = 0; i < nof_zv_elements; i++) {

    int body_id = body_ids[i];

    if ( body_id == NO_INDEX )
      break;
    
    Body* body = model->getBodyById(body_id);

    if ( body != NULL ) {
      body->addMeshElements(nof_new_bbulks[i], new_bbulk_ids[i]);
    }
  }

  // Delete buffers
  // ==============
  for (i = 0; i < nof_zv_elements; i++) {
    delete[] new_bbulk_ids[i];
  }

  delete[] bulk_remove_flags;
  delete[] body_ids;
  delete[] nof_new_bbulks;
  delete[] new_bbulk_ids;


  // Update environment 
  // ==================
  gui->updateModelStatistics(model);

  model->refreshRenderer();
}



void
ModelMeshManager::createMeshBodies()
{
  int i;
  int nof_bodies = meshInfo->nofBodies;

  //-----Find fem elements per body
  //-We store this data in temporary tables
  int* nof_body_fem_elements = new int[1 + MAX_NOF_BODIES];
  int** body_fem_element_ids = new int*[1 + MAX_NOF_BODIES];

  for (i = 0; i < 1 + MAX_NOF_BODIES; i++) {
    nof_body_fem_elements[i] = 0;
    body_fem_element_ids[i] = NULL;
  }

  //-Calculate first fem elements per body (for alloc size)
  for (i = 0; i < meshInfo->nofBulkElements; i++) {
    int ext_bd_tg = meshData->bulkElements->parentIds[i][0];
    int int_bd_tg = meshData->bodyExt2Int[ext_bd_tg];
    nof_body_fem_elements[int_bd_tg]++;
  }

  //-Allocate space for fem element ids
  for (i = 0; i < 1+ MAX_NOF_BODIES; i++) {

    int count = nof_body_fem_elements[i];
    nof_body_fem_elements[i] = 0; // we use this for position counter below!

    if (count > 0) {
      body_fem_element_ids[i] = new int[count];
    }
  }

  //-Create internal bulk element indices for bodies
  // (use sequential numbering implied by the bulk-table!)
  int position;
  for (i = 0; i < meshInfo->nofBulkElements; i++) {

    int ext_bd_tg = meshData->bulkElements->parentIds[i][0];
    int int_bd_tg = meshData->bodyExt2Int[ext_bd_tg];
    position = nof_body_fem_elements[int_bd_tg];
    body_fem_element_ids[int_bd_tg][position] = i;
    nof_body_fem_elements[int_bd_tg] = ++position;
  }

  //-----Create now all bodies
  for (i = 0; i < nof_bodies; i++) {

    int int_bd_tag = i + 1;  // internal tag are: 1...

    int ext_bd_tag = meshData->bodyInt2Ext[int_bd_tag];

    int nof_fem_elems = nof_body_fem_elements[int_bd_tag];
    int* fem_elem_ids = body_fem_element_ids[int_bd_tag];
    Body* body = NULL;

    // Check if this body is already created
    body = model->getBodyByTag(int_bd_tag);

    if (body != NULL) {
      body->setMeshElements(nof_fem_elems, fem_elem_ids);
      body->setExternalTag(ext_bd_tag);
      continue;
    }

    // Create a new body
    switch (modelInfo->dimension) {
    case 2:
      body = new Body2D(MESH_BODY, int_bd_tag, ext_bd_tag, nof_fem_elems, fem_elem_ids);
      break;
    case 3:
      body = new Body3D(MESH_BODY, int_bd_tag, ext_bd_tag, nof_fem_elems, fem_elem_ids);
      break;
    }

    model->addBody(body);

  } // for nof-mesh-bodies

  // NOTE: delete only tables, NOT the items in the tables!
  delete[] nof_body_fem_elements;
  delete[] body_fem_element_ids;
}




// Construct body <--> material table
//
void 
ModelMeshManager::createMeshBodyTables()
{
  int i;

  // Calc nof bodies from material counter table if
  // it is not yet set
  if (meshInfo->nofBodies == 0) {
    for(i= 0; i < MAX_NOF_BODIES; i++) {
      if (meshData->bodyExt2Int[i] != NO_INDEX) {
        meshInfo->nofBodies++;
      }
    }
  }

  // Construct internal body id <--> external body id (= material id) -tables
  // Condition: external ids >= 0 !!!***!!!
  int counter = 1;
  for(i= 0; i < MAX_NOF_BODIES; i++) {
    if (meshData->bodyExt2Int[i] != NO_INDEX) {
      meshData->bodyExt2Int[i] = counter;
      meshData->bodyInt2Ext[counter] = i;
      counter++;
    }
  }
}


// Create mesh geometry boundaries
// Add mesh boundary elements to the boundaries
//
// Argument: free_bulk_bndr_flags (possibly) still free mesh bulk boundary 
// element are flagged here.
// If NULL, all mesh boundary elements are created here!
//
void
ModelMeshManager::createMeshBoundaries(int nof_bulk_bndr_elems, bool* free_bulk_bndr_elems)
{
  int i, j, counter;

  int nof_bodies = meshInfo->nofBodies;

  //-Create a table to store boundary ids
  // stored in an upper triangle matrix using body1,body2
  // indices as keys
  //
  int** bodyPairTable = new int*[nof_bodies];
  int** boundaryIdTable = new int*[nof_bodies];
  BodyElement*** boundaryTable = new BodyElement**[nof_bodies];

  for (i = 0; i < nof_bodies; i++) {

    bodyPairTable[i] = new int[nof_bodies];
    boundaryIdTable[i] = new int[nof_bodies];
    boundaryTable[i] = new BodyElement*[nof_bodies];

    for (int j = 0; j < nof_bodies; j++) {
      bodyPairTable[i][j] = 0;
      boundaryIdTable[i][j] = NO_INDEX;
      boundaryTable[i][j] = NULL;
    }
  }

  MeshElementTable* bt = meshData->bulkElements;

  // Calculate nof boundary elments per boundary
  // ============================================
  counter = -1;
  bool is_first_parent;
  int nof_bndr_elems;
  meshElementCode bulk_code;
  

  for (int bulk_id = 0; bulk_id < bt->NofElements(); bulk_id++) {

    bulk_code = bt->getElementCode(bulk_id);
    
    nof_bndr_elems = MeshElementDesc[bulk_code][DESC_NOF_BNDR_ELEMS];

    // Pick each face and check if it is a boundary element
    for (int face = 0; face < nof_bndr_elems; face++) {

      // NOTE: counter is correct if we take only the first parent, because
      // this was the parent which created the boundary element!
      //
      if ( !bt->isBoundaryFace(bulk_id, face, is_first_parent) ||
           !is_first_parent
         ) {
        continue;
      }

      counter++;

      // If boundary element is already flagged (used), skip it
      if ( free_bulk_bndr_elems != NULL && 
           !free_bulk_bndr_elems[counter]
           ) {
        continue;
      }

      //---Check if face was a boundary element
      int body1_id = bt->parentIds[bulk_id][0];
      int body2_id = NO_INDEX;

      int neighbor_id = bt->neighborIds[bulk_id][face];

      if ( neighbor_id != NO_INDEX) {
        body2_id = bt->parentIds[neighbor_id][0];
      }

      int body1_tag, body2_tag;

      model->getBodyTags(body1_id, body2_id, body1_tag, body2_tag);

      if ( body2_tag == NO_INDEX ) {
        body2_tag = body1_tag;
      }
        
      // Update body-pair element counter
      if (body1_tag <= body2_tag) {
        bodyPairTable[body1_tag - 1][body2_tag - 1]++;
      } else {
        bodyPairTable[body2_tag - 1][body1_tag - 1]++;
      }

    } // for all bulk element faces      

  } // for all bulk elements


  // Create outer/inner boundaries
  // =============================
  //-Take a body corresponding the row in the meshBodyPairTable
  BodyElement* be;

  for (i = 0; i < nof_bodies; i++) {

    for (j = i; j < nof_bodies; j++) {

      int count = bodyPairTable[i][j];

      // If no common element
      if (count == 0) {
        continue;
      }

      int body1_tag = i + 1;
      int body2_tag = j + 1;

      int body1_lr = 0;
      int body2_lr = 0;

      // Outer boundary
      if (i == j)
        body2_tag = NO_INDEX;

      // Try get get possible existing boundary
      be = model->getBoundaryByTags(body1_tag, body2_tag);

      //--Use existing element, IF it does not have any
      //  mesh elements yet!
      if ( be != NULL && be->getNofMeshElements() == 0 ) {
        be->allocateMeshElements(count);

      //--Otherwise, create new element
      } else {
        be = model->createBodyElement(body1_tag, body1_lr, body2_tag, body2_lr, count);
      }

      if ( be != NULL ) {
        // Store boundary and its tag
        boundaryTable[i][j] = be;
        boundaryIdTable[i][j] = be->Tag();
      }

    } // for columns (first body)
  } // for rows (second body)

 
  // Create mesh boundary elements
  // =============================
  createMeshBoundaryElements(nof_bulk_bndr_elems, free_bulk_bndr_elems, boundaryTable);

  // Delete work tables
  for (i = 0; i < nof_bodies; i++) {
    delete[] bodyPairTable[i];
    delete[] boundaryIdTable[i];
    delete[] boundaryTable[i];
  }

  delete[] bodyPairTable;
  delete[] boundaryIdTable;
  delete[] boundaryTable;
}


void
ModelMeshManager::createMeshBoundaryElements(int nof_bulk_bndr_elems, 
                                             bool* free_bulk_bndr_elems,
                                             BodyElement*** boundary_table)
{
  int node_ids[27];

  MeshElementTable* bt = meshData->bulkElements;

  int counter = -1;
  bool is_first_parent;
  bool is_added_to_bndr;
  meshElementCode bulk_code;
  meshElementCode bndr_code;
  int nof_bndr_elems;
  const int* bulk_bndr_node_ids;

  for (int bulk_id = 0; bulk_id < bt->NofElements(); bulk_id++) {

    bulk_code = bt->getElementCode(bulk_id);

    nof_bndr_elems = MeshElementDesc[bulk_code][DESC_NOF_BNDR_ELEMS];

    // Pick each face and check if it is a boundary element
    for (int face = 0; face < nof_bndr_elems; face++) {

      // NOTE: counter is correct if we take only the first parent, because
      // this was the parent which created the boundary element!
      //
      if ( !bt->isBoundaryFace(bulk_id, face, is_first_parent) ||
          !is_first_parent
         ) {
        continue;
      }

      counter++;

      // If boundary element are already flagged for us, skip non-free (used)
      // elements
      if ( free_bulk_bndr_elems != NULL && 
           !free_bulk_bndr_elems[counter]
         ) {
        continue;
      }
     
      int body1_id = bt->parentIds[bulk_id][0];
      int body2_id = NO_INDEX;

      int neighbor_id = bt->neighborIds[bulk_id][face];

      if ( neighbor_id != NO_INDEX) {
        body2_id = bt->parentIds[neighbor_id][0];
      }
 
      int body1_tag, body2_tag;

      model->getBodyTags(body1_id, body2_id, body1_tag, body2_tag);
      
      bndr_code = MeshElementBndrCodes[bulk_code][face];

      int nof_nodes = MeshElementDesc[bndr_code][DESC_NOF_NODES];

      // Pick bulk nodes ids
      const int* bulk_node_ids = bt->getNodeIds(bulk_id);

      // Pick local bndr node ids
      bulk_bndr_node_ids = MeshElementBndrNodes[bulk_code][face];

      for (int i = 0; i < nof_nodes; i++) {
        node_ids[i] = bulk_node_ids[bulk_bndr_node_ids[i]];
      }

      BodyElement* be;

      // From upper triangle, bodies are numbered 1,...
      if (body2_tag == NO_INDEX) {
        body2_tag = body1_tag;
      }

      if (body1_tag <= body2_tag) {
        be = boundary_table[body1_tag - 1][body2_tag - 1];
      } else {
        be = boundary_table[body2_tag - 1][body1_tag - 1];
      }

      // Create new mesh bndr element entry
      int elem_index = addMeshBoundaryElement(be, bndr_code, bulk_id, neighbor_id,
                                              nof_nodes, node_ids,
                                              is_added_to_bndr);

      // Dir in parents
      short dir_in_parent1, dir_in_parent2;
 
      dir_in_parent1 = 1;
    
      if (body2_id != NO_INDEX) {
        dir_in_parent2 = -1;
      } else {
        dir_in_parent2 = 0;
      }

      meshData->boundaryElements->setDirInParents(elem_index, dir_in_parent1, dir_in_parent2);

    } // for all faces in the bulk element

  } // for all bulk elements
}


//---Create mesh element edges (edges in 3D, edges or vertices! in 2D)
void
ModelMeshManager::createMeshElementEdges(MeshElementTable* source, MeshElementTable* target)
{
  int node_ids[27];

  //--Loop all source ELEMENTS
  for (int elem_id = 0; elem_id < source->NofElements(); elem_id++) {

    meshElementCode elem_code = source->getElementCode(elem_id);

    // Source edge info
    int nof_edges = MeshElementDesc[elem_code][DESC_NOF_EDGES];
    meshElementCode edge_code = (meshElementCode)MeshElementDesc[elem_code][DESC_EDGE_ELEM_CODE];

    // Get elem nodes
    const int* elem_nodes = source->getNodeIds(elem_id);
 
    //--Loop each EDGE in the source element
    for (int edge = 0; edge < nof_edges; edge++) {

      int edge_id = source->edgeIds[elem_id][edge];

      // If edge already created
      if ( target->getElementCode(edge_id) != MEC_000 ) {
        continue;
      }

      // Pick edge nodes from the element
      int nof_edge_nodes = MeshElementDesc[edge_code][DESC_NOF_NODES];
      int* elem_edge_node_ids = (int*)MeshElementEdgeNodes[elem_code][edge];

      for (int i = 0; i < nof_edge_nodes; i++) {
        node_ids[i] = elem_nodes[elem_edge_node_ids[i]];
      }

      target->createTableEntry(edge_id, edge_code, NULL, nof_edge_nodes, node_ids);

    } // for each edge in the source element

  } // for all source elements
}


//---Create mesh geometry subelements
// Argument hlevel (hierarchy level):
//   hlevel = 3: edges    (in 3D)
//   hlevel = 2: vertices (in 3D or 2D)
//
// NOTE: Here we are checking mesh element connectivity to body elements
// (face,edges,vertices), not to other mesh elements! So, in 3D meshelements
// which belong to two boundaries form an edge element etc.
//
void
ModelMeshManager::createMeshSubElements(int hlevel, int min_nof_parents)
{
  MeshElementTable* sub_source;
  MeshElementTable* sub_table;

  // 3D model
  // ========
  if ( modelInfo->dimension == ECIF_3D ) {

    if ( hlevel == 3 ) {
      sub_source = meshData->boundaryElements;
      sub_table = meshData->boundaryEdges;

    } else if ( hlevel == 2 ) {
      sub_source = meshData->boundaryEdges;
      sub_table = meshData->boundaryVertices;

    } else {
      return;
    }
  
  // 2D model
  // ========
  } else {
    if ( hlevel = 2 ) {
      sub_source = meshData->boundaryElements;
      sub_table = meshData->boundaryVertices;

    } else {
      return;
    }
  }

  int count = sub_table->NofElements();
  MeshConnectionTable* connInfo = new MeshConnectionTable(count, 5, 3);

  // Loop all higer level elements
  // =============================
  BodyElement* be;
  int index = 0;

  while (true) {

    if ( hlevel == 3 ) {
      be = model->getFace(index++);
    } else {
      be = model->getEdge(index++);
    }
    
    if (be==NULL) break;
  
    // Mesh boundary border is subelement for a Bem boundary!
    //
    // NOTE: Currently ElmerSolver's BEM solver does not actually support
    // any 'sublements', so we do not use this. In any this is very coarse way
    // to 'see' bondary edges, we should add all necessary sub edges and not only the
    // whole border!!!###!!! MVe 19.05.03
    //
#if 0
    if ( be->isBemBoundary() ) {
      be->addMeshBorderAsSubElement();
    }
#endif

    int nof_elems = be->getNofMeshElements();

    // Loop all mesh elements in the parent
    for (int i = 0; i < nof_elems; i++) {

      int elem_id = be->getMeshElementId(i);

      meshElementCode elem_code = sub_source->getElementCode(elem_id);
      short nof_sub_elems = MeshElementDesc[elem_code][DESC_NOF_BNDR_ELEMS];

      // Add parent info for each mesh subelement
      for (short l = 0; l < nof_sub_elems; l++) {

        int sub_index = sub_source->edgeIds[elem_id][l];

        // NOTE: We store parent object ids, only unique occurances!
        connInfo->addConnection(sub_index, be->Id(), true);
      }
    }
  }

  bool* checked = new bool[count];  // To flag out unsuitable or already handled mesh subelements
  bool* is_seed = new bool[count];  // True if this mesh subelement originated the new body subelement
  int* subTags = new int[count];    // To store the body subelement tag where this mesh subelement belongs

  int i;

  // Loop all mesh subelements and check if it
  // fullfills the minimum connection criterium
  // ===========================================
  for (i = 0; i < count;  i++) {

    int nof_conn = connInfo->getNofConnections(i);

    // Flagged away (is unusable!)
    if ( nof_conn < min_nof_parents ) {
      checked[i] = true;

    // Still usable
    } else {
      checked[i] = false;
    }

    is_seed[i] = false;
    subTags[i] = NO_INDEX;
  }


  int sub_counter = 0;

  // Loop all mesh subelements to find the body subelement where it
  // ==============================================================
  // possibly belongs
  // ================
  // NOTE: Now all nonchecked belong to some body subeleemnt!
  for (i = 0; i < count;  i++) {

    if (checked[i])
      continue;
    
    // Pick first nonchecked and make it a new
    // edge seed
    if ( subTags[i] == NO_INDEX ) {

      checked[i] = true;
      is_seed[i] = true;
      subTags[i] = ++sub_counter;

      if ( hlevel == 2 ) {
        continue;
      }

      // In level 3 (for face edges) check against all the rest (nonchecked)
      // if they have same connections and then also belong
      // to the same body subelement!
      int nof_conn_i = connInfo->getNofConnections(i);
      const int* pr_ids_i = connInfo->getParentIds(i);

      for (int j = i + 1; j < count; j++) {

        if (checked[j]) {
          continue;
        }

        if ( nof_conn_i != connInfo->getNofConnections(j) ) {
          continue;
        }

        const int* pr_ids_j = connInfo->getParentIds(j);

        bool is_same_se = true;

        for (int k = 0; k < nof_conn_i; k++) {

          if ( pr_ids_i[k] != pr_ids_j[k] ) {
            is_same_se = false;
            break;
          }
        }

        if (is_same_se) {
          checked[j] = true;
          subTags[j] = subTags[i];
        }

      } // For all rest from i

    } // If new boundry seed found

  } // For all from beginning


  // Count nof meshelements in each body subelement
  // ==============================================
  int* elem_counts = new int[1 + sub_counter];
  int* sub_obj_ids = new int[1 + sub_counter];

  for (i = 0; i <= sub_counter; i++) {
    elem_counts[i] = 0;
    sub_obj_ids[i] = NO_INDEX;
  }

  for (i = 0; i < count; i++) {
    if ( subTags[i] != NO_INDEX ) {
      elem_counts[subTags[i]]++;
    }
  }

  // Create bodyelements and allocate space for meshelements
  //========================================================
  static GcPoint p;
  double* np;

  for (i = 0; i < count; i++) {

    int sub_tag = subTags[i];

    if ( sub_tag == NO_INDEX || !is_seed[i] )
      continue;

    // In 3D highest level subelement is an edge
    if ( hlevel == 3 ) {
      be = new BodyElement2D(sub_tag);
      be->allocateMeshElements(elem_counts[sub_tag]);

    // Otherwise subelement is a vertex, so we need the node-id
    // in the vertex
    } else {
      const int* nd_id = sub_table->getNodeIds(i);
      np = meshData->nodes[nd_id[0]];
      p.setPosit(np[0], np[1], np[2]);
      be = new BodyElement1D(sub_tag, &p);
    }

    sub_obj_ids[sub_tag] = be->Id();

    model->addBodyElement(be);
  }


  // Add mesh elements to body subelements and add body subelements
  // ===============================================================
  // to parents
  // ==========
  for (i = 0; i < count; i++) {

    if ( subTags[i] == NO_INDEX ) {
      continue;
    }

    int pe_id = (sub_obj_ids[subTags[i]]);

    if ( hlevel == 3 ) {
      be = model->getEdgeById(pe_id);
    } else {
      be = model->getVertexById(pe_id);
    }

    // Update element type counters for Elmer Mesh header!
    meshElementCode elem_code = sub_table->getElementCode(i);
    meshInfo->nofUsedElementTypes[elem_code]++;
    
    be->addMeshElement(i, 0);

    // If first occurance of the subelement, add it
    // to the parent element
    if ( is_seed[i] ) {

      for (int j = 0; j < connInfo->getNofConnections(i); j++) {

        int pe_id = connInfo->getParentId(i, j);

        BodyElement* pe;
        
        if ( hlevel == 3 ) {
          pe = model->getFaceById(pe_id);
        } else {
          pe = model->getEdgeById(pe_id);
        }

        be->setParentId(1, pe->getParentId(1));
        be->setParentLayer(1, pe->getParentLayer(1));

        pe->addSubElement(be);
      }
    }
  }

}



void
ModelMeshManager::findMeshBoundaryBorders()
{
  int index = 0;
  while (true) {
    BodyElement* be = model->getBoundary(index++);
    if (be==NULL) break;
    be->findMeshBorder();
  }
}


//---Find nof bndr elements defined by bulk elements
//
// NOTE: Actual boundary elements entries will be created later!!!
//
void
ModelMeshManager::findNofBulkBoundaryElements(int& nof_bndr_elements)
{
  MeshElementTable* source = meshData->bulkElements;
  
  nof_bndr_elements = 0;

  // Loop all bulk elements
  //
  for (int elem_id = 0; elem_id < source->NofElements(); elem_id++) {
    int elem_code = source->getElementCode(elem_id);
    int nof_bndr_elems = MeshElementDesc[elem_code][DESC_NOF_BNDR_ELEMS];

    // Loop all faces in the bulk element
    //
    for (int face = 0; face < nof_bndr_elems; face++) {
    
      bool is_first_parent;

      if ( !source->isBoundaryFace(elem_id, face, is_first_parent) ||
           !is_first_parent
           ) {
        continue;
      }

      nof_bndr_elements++;

    } // for each face 
  } // for each bulk

}



// Methods finds parent bodies for each mesh boundary.
// Matching is done by finding the parent bulk element and
// corresponding subelement index for the first element
// in the boundary.
//
// Reference argument: free_bndr_element_indices will flag all those bulk boundary elements
// (faces/edges) which are not matched to any existing boundary, ie. those boundary elements
// which were not members of any boundary list in the input.
//
// NOTE: Those boundary elements in the input which are not bulk faces/edges must be
// 'BEM' boundaries and they define a boundary and a body at the same time!!!
//
void
ModelMeshManager::findMeshBoundaryParents(int nof_bulk_bndr_elems, bool* free_bulk_bndr_flags)
{
  int node_ids[27];

  const ModelStatistics* ms = model->getModelStatistics();

  MeshElementTable* bt = meshData->bulkElements;
  MeshElementTable* bet = meshData->boundaryElements;

  int i, counter, index;

  MeshConnectionTable* bulk_bndr2nodes = new MeshConnectionTable(meshInfo->nofNodes, 5, 3);

  // Allocate flag tables
  int* bulk_bndr_parents1 = new int[nof_bulk_bndr_elems];
  int* bulk_bndr_parents2 = new int[nof_bulk_bndr_elems];

  // Initialize flags
  //
  for (i = 0; i < nof_bulk_bndr_elems; i++) {
    free_bulk_bndr_flags[i] = true;
  }

  meshElementCode bulk_code;
  meshElementCode bndr_code;
  int nof_bndr_elems;
  const int* bulk_bndr_node_ids;

  //--1. Find all bndr elements in the bulk elements table
  //     Fill also bulk sub elements connection table
  //
  counter = 0;
  for (int bulk_id = 0; bulk_id < bt->NofElements(); bulk_id++) {

    bulk_code = bt->getElementCode(bulk_id);

    nof_bndr_elems = MeshElementDesc[bulk_code][DESC_NOF_BNDR_ELEMS];

    // Pick each face and check if it is a boundary element
    for (int face = 0; face < nof_bndr_elems; face++) {
      
      bool is_first_parent;

      // NOTE: counter is correct if we take only the first parent, because
      // this was the parent which created the boundary element!
      //
      if ( !bt->isBoundaryFace(bulk_id, face, is_first_parent) ||
          !is_first_parent
         ) {
        continue;
      }

      bulk_bndr_parents1[counter] = bulk_id;

      bulk_bndr_parents2[counter] = bt->neighborIds[bulk_id][face];

      const int* bulk_node_ids = bt->getNodeIds(bulk_id);

      bndr_code = MeshElementBndrCodes[bulk_code][face];

      int nof_match_nodes = MeshElementDesc[bndr_code][DESC_NOF_MATCH_NODES];

      // Pick local bndr node ids
      bulk_bndr_node_ids = MeshElementBndrNodes[bulk_code][face];

      for (int i = 0; i < nof_match_nodes; i++) {
        node_ids[i] = bulk_node_ids[bulk_bndr_node_ids[i]];
      }

      bulk_bndr2nodes->addConnections(nof_match_nodes, (int*)node_ids, counter, true);

      counter++;

    } // all faces in the bulk

  } // all bulk elements

  bulk_bndr2nodes->sortParentIds();

  // Table of original boundary sets 
  // We have to read boundary sets from this table
  // because we may create some new boundaries when processing
  // old boundaries and model's getFirst/getNextBoundary paradigm would not
  // then be useable!
  //
  int nof_boundary_sets = ms->nofElements;
  BodyElement** boundary_sets;
  bool** delete_flags;
  bool* delete_indicators;

  boundary_sets = new BodyElement*[nof_boundary_sets];
  delete_flags = new bool*[nof_boundary_sets];
  delete_indicators = new bool[nof_boundary_sets];

  counter = 0;
  
  index =0;
  while (true) {
    BodyElement* be = model->getBoundary(index++);
    if (be==NULL) break;
    boundary_sets[counter] = be;
    counter++;
  }

  Ids3 ids3;


  //--2. Match mesh elements in the boundaries with the bulk bndr-elements

  // Loop all boundaries
  // ===================
  for (int set_index = 0; set_index < nof_boundary_sets; set_index++) {

    // In this set of ids-triplets (body1-id,body2-id,bndr-id)
    // we store info on all boundaries which possible are
    // created from THIS boundary element set (including itself)
    //
    Ids3Set* body_bndr_ids = new Ids3Set;
    int nof_created_boundaries = 0;

    BodyElement* be = boundary_sets[set_index];

    int nof_elements = be->getNofMeshElements();

    delete_flags[set_index] = new bool[nof_elements];
    delete_indicators[set_index] = false;

    int body1_tag, body2_tag;

    int buffer_size = 1024;
    int ids_buffer[1024];
    int nof_ids;

    bool match_found;

    bool is_first_mesh_element = true;

    // Loop all mesh elements in the boundary
    // ======================================
    for (int set_elem_index = 0; set_elem_index < nof_elements; set_elem_index++) {

      int bndr_elem_id = be->getMeshElementId(set_elem_index);

      delete_flags[set_index][set_elem_index] = false;

      // Final bodyelement for the mesh element!
      BodyElement* be_for_mesh_element = NULL;

      meshElementCode elem_code = bet->getElementCode(bndr_elem_id);

      int nof_match_nodes = MeshElementDesc[elem_code][DESC_NOF_MATCH_NODES];
      int* node_ids = bet->nodeIds[bndr_elem_id];

      // No match, the boundary is a 'BEM boundary'!
      // -------------------------------------------
      if ( !bulk_bndr2nodes->getParentIds(nof_match_nodes, node_ids,
                                          buffer_size, ids_buffer,
                                          nof_match_nodes, nof_ids)
         ) {
        match_found = false;
        break;
      }
      
      // Set parent for the possible input element
      setMeshInputElementParentTag(bndr_elem_id, be->Tag());
      setMeshInputElementIsAdded(bndr_elem_id, true);

      match_found = true;

      // Set mesh boundary element data
      // ------------------------------
      for (int i = 0; i < nof_ids; i++) {

        int bulk_bndr_id = ids_buffer[i]; // This is the "counter" above

        // This bulk bndr element is not free!
        //
        free_bulk_bndr_flags[bulk_bndr_id] = false;

        // Pick parent info from the matching bulk sub element
        int bulk_elem1_id = bulk_bndr_parents1[bulk_bndr_id];
        int bulk_elem2_id = bulk_bndr_parents2[bulk_bndr_id];

        bet->parentIds[bndr_elem_id][0] = bulk_elem1_id;
        bet->parentIds[bndr_elem_id][1] = bulk_elem2_id;

        // Find parent body info for the bndr element and
        // add it to the corresponding boundary 
        // (Parent info for the boundary is not yet known!)
        //
        int body1_id, body2_id;

        body1_id = bt->parentIds[bulk_elem1_id][0];

        if (bulk_elem2_id != NO_INDEX)
          body2_id = bt->parentIds[bulk_elem2_id][0];
        else
          body2_id = NO_INDEX;
 
        model->getBodyTags(body1_id, body2_id, body1_tag, body2_tag);

        // First body-id should be the smaller
        if ( body2_tag != NO_INDEX && body2_tag < body1_tag) {
          int tmp = body1_tag;
          body1_tag = body2_tag;
          body2_tag = tmp;
        }
        
        int body1_lr = 0;
        int body2_lr = 0;
       
        // If original set does not yet have parent body ids,
        // update parent ids and add boundary to the model
        //
        if ( is_first_mesh_element ) {

          if ( !model->modelHasCadGeometry() ) {
            be->setParentTags(body1_tag, body2_tag);
            model->addBodyElement(be, body1_tag, body1_lr, body2_tag, body2_lr);
          }

          // Store info in the bndr ids set
          ids3.id1 = body1_tag;
          ids3.id2 = body2_tag;
          ids3.id3 = be->Tag();
          body_bndr_ids->insert(ids3);

          is_first_mesh_element = false;
        }

        // Now, all necessary info should be available: find boundary
        // which matches the body-pair ids
        //
        int be_tag = find3(*body_bndr_ids, body1_tag, body2_tag);

        // OLD boundary
        if (be_tag != NO_INDEX) {
          be_for_mesh_element = model->getBoundaryByTag(be_tag);

        // NEW boundary 
        } else {

          be_for_mesh_element = model->createBodyElement(body1_tag, body1_lr, body2_tag, body2_lr, 0);
          
          if ( be_for_mesh_element == NULL ) {
            match_found = false;
            break;
          }

          nof_created_boundaries++;

          // Store info on the new boundary
          ids3.id1 = body1_tag;
          ids3.id2 = body2_tag;
          ids3.id3 = be_for_mesh_element->Tag();
          body_bndr_ids->insert(ids3);

          // Set name for the new boundary
          strstream strm;
          strm << be->getName() << nof_created_boundaries << ends;
          be_for_mesh_element->setName(strm.str());
        }

        // Set dir-in-parents for the boundary element
        //
        short dir_in_parent1, dir_in_parent2;
 
        dir_in_parent1 = 1;
    
        if (bulk_elem2_id != NO_INDEX) {
          dir_in_parent2 = -1;
        } else {
          dir_in_parent2 = 0;
        }

        // If the boundary is not the same as the original set
        // add mesh element to the new boundary and mark it to
        // be deleted form the boundary set
        //
        if ( be_for_mesh_element != be) {

          be_for_mesh_element->addMeshElement(bndr_elem_id, dir_in_parent1);
          
          delete_flags[set_index][set_elem_index] = true;
          delete_indicators[set_index] = true;

        // Otherwise just update current element dir in the set
        } else {
          be_for_mesh_element->setMeshElementDir(set_elem_index, dir_in_parent1);
        }

      } // For all sub-ids found in connections

    } // For all mesh elements in the boundary

    // Delete mesh elements moved from the set into new boundaries
    //
    if ( match_found && delete_indicators[set_index] ) {

      be->deleteMeshElements(delete_flags[set_index]);

      // If no mesh elements left in the original set, remove
      // it from the model
      if ( be->getNofMeshElements() == 0 ) {
        Body* body1 = model->getBodyByTag(body1_tag);
        Body* body2 = model->getBodyByTag(body2_tag);
        model->removeBodyElement(be, body1, body2, true);
        delete be;
      }
    }

    delete body_bndr_ids;

  } // For all boundaries


  // --3. Clean data
  // 

  delete bulk_bndr2nodes;

  delete[] bulk_bndr_parents1;
  delete[] bulk_bndr_parents2;

  for (i = 0; i < nof_boundary_sets; i++) {
    delete[] delete_flags[i];
  }
  delete[] delete_flags;
  delete[] delete_indicators;
}


//---Find mesh corner element: bulk elements whose all nodes are at boundaries
void
ModelMeshManager::findMeshCornerElements()
{
  static bool node_flags[27];
  int i;

  modelData->purgeMeshCornerElements();
  modelData->meshCornerElements = new MeshCornerElementList;

  MeshElementTable* bulk_table = meshData->bulkElements;
  MeshElementTable* bndr_table = meshData->boundaryElements;

  int corner_count = 0;

  // Loop all bulk elements
  // ======================
  for (int bulk_id = 0; bulk_id < bulk_table->NofElements(); bulk_id++) {
  
    meshElementCode bulk_code = bulk_table->getElementCode(bulk_id);
    int nof_nodes = MeshElementDesc[bulk_code][DESC_NOF_NODES];
    int nof_bndr_elems = MeshElementDesc[bulk_code][DESC_NOF_BNDR_ELEMS];

    // Init flags
    for (i = 0; i < nof_nodes; i++) {
      node_flags[i] = false;
    }

    bool is_first_parent; // This is not actually used here!

    // Pick each face and check if it is a boundary element
    // Mark all nodes at boundary face as "boundary nodes"
    for (int face = 0; face < nof_bndr_elems; face++) {

      if ( bulk_table->isBoundaryFace(bulk_id, face, is_first_parent) ) {

        meshElementCode bndr_code = MeshElementBndrCodes[bulk_code][face];
        int nof_bndr_nodes = MeshElementDesc[bndr_code][DESC_NOF_NODES];
        const int* bndr_nodes = MeshElementBndrNodes[bulk_code][face];

        for (i = 0; i < nof_bndr_nodes; i++) {
          node_flags[bndr_nodes[i]] = true;
        }
      }
    }
  
    // If all nodes are really at boundary, we
    // have a "corner" element
    bool at_corner = true;
    for (i = 0; i < nof_nodes; i++) {

      if ( !node_flags[i] ) {
        at_corner = false;
        break;
      }
    }

    if (!at_corner) {
      continue;
    }
    
    // Create new corner element info
    // ==============================
    MeshCornerElement* corner = new MeshCornerElement();

    corner->elementId = bulk_id;
    corner->elementCode = bulk_code;
    corner->bodyId = bulk_table->parentIds[bulk_id][0];
    corner->nodeIds = new int[nof_nodes];

    bulk_table->centerPoint(bulk_id, corner->centerPoint);

    // Copy node ids
    const int* node_ids = bulk_table->getNodeIds(bulk_id);

    for (i = 0; i < nof_nodes; i++) {
      corner->nodeIds[i] = node_ids[i];
    }

    corner->nofSubElements = nof_bndr_elems;
    corner->subElementIds = new int[nof_bndr_elems];
    corner->boundaryIds = new int[nof_nodes];
    corner->boundaryElementIds = new int[nof_bndr_elems];

    // Init data
    for (i = 0; i < nof_nodes; i++) {
      corner->boundaryIds[i] = NO_INDEX;
    }

    for (i = 0; i < nof_bndr_elems; i++) {
      corner->boundaryElementIds[i] = NO_INDEX;
    }

    modelData->meshCornerElements->push_back(corner);
    corner_count++;

#if 0
    strstream strm;
    strm << corner_count << ". **Corner element: " << elem_id << endl << ends;
    theControlCenter->getGui()->showMsg(strm.str());
#endif

    // Find boundary ids (loop all boundary elements!)
    // =================
    for (int bndr_elem_id = 0; bndr_elem_id < bndr_table->NofElements(); bndr_elem_id++) {

      int bulk1_id = bndr_table->parentIds[bndr_elem_id][0];
      int bulk2_id = bndr_table->parentIds[bndr_elem_id][1];

      // Face found (this bndr element is a face in the current bulk)
      // ----------
      if ( bulk1_id == bulk_id || bulk2_id == bulk_id ) {

        //--Find bndr element face index in the parent
        meshElementCode bndr_elem_code = bndr_table->getElementCode(bndr_elem_id);
        int nof_bndr_nodes = MeshElementDesc[bndr_elem_code][DESC_NOF_NODES];
        int nof_match_nodes = MeshElementDesc[bndr_elem_code][DESC_NOF_MATCH_NODES];
        const int* bndr_elem_node_ids = bndr_table->getNodeIds(bndr_elem_id);
 
        int sub_dir, sub_start_pos;
        int face = bulk_table->findSubElementIndex(bulk_id, sub_dir, sub_start_pos,
                                                   nof_match_nodes, bndr_elem_node_ids);

        //--This should not happen!
        if ( face == NO_INDEX ) {
          strstream strm;
          strm << "Corner element: " << bulk_id
               << "   Face index not found for the boundary element " << bndr_elem_id
               << endl << ends;
          theControlCenter->getGui()->showMsg(strm.str());

          continue;
        }

        int bndr_id = meshData->boundaryElementBoundaryIds[bndr_elem_id];


        // Store boundary info for bulk nodes
        // ----------------------------------
        const int* bndr_nodes = MeshElementBndrNodes[bulk_code][face];
        for (i = 0; i < nof_bndr_nodes;  i++) {

          int nd_index = bndr_nodes[i];

          if ( corner->boundaryIds[nd_index] < bndr_id ) {
            corner->boundaryIds[nd_index] = bndr_id;
          }
        }

        // Store bndr element id for bulk face
        corner->boundaryElementIds[face] = bndr_elem_id;

      }

    } // For all boundary elements
  
  } // For all bulk elements

  meshInfo->nofCornerElements = corner_count;
}


//---Find global mesh edge element ids for edges in the source mesh element table
//
void
ModelMeshManager::findMeshElementEdges(MeshElementTable* source, int& nof_edge_elements)
{
  nof_edge_elements = 0;

  int nof_nodes = meshInfo->nofNodes;
 
  MeshConnectionTable* edgeInfos = new MeshConnectionTable(nof_nodes, 5, 3, true);

  //--Loop all source ELEMENTS
  for (int elem_id = 0; elem_id < source->NofElements(); elem_id++) {

    meshElementCode elem_code = source->getElementCode(elem_id);

    // Edge info
    int nof_edges = MeshElementDesc[elem_code][DESC_NOF_EDGES];
    meshElementCode edge_code = (meshElementCode)MeshElementDesc[elem_code][DESC_EDGE_ELEM_CODE];

    // Get element nodes
    const int* elem_nodes = source->getNodeIds(elem_id);
 
    //--Loop each EDGE in the element
    for (int edge = 0; edge < nof_edges; edge++) {

      int edge_id = source->edgeIds[elem_id][edge];

      // Already set!
      if ( edge_id != NO_INDEX ) {
        continue;
      }

      // Pick edge end nodes from the element
      int nof_edge_nodes = MeshElementDesc[edge_code][DESC_NOF_NODES];
      int nof_match_nodes = 2;

      if ( nof_edge_nodes == 1 ) {
        nof_match_nodes = 1;
      }

      int* elem_edge_node_ids = (int*)MeshElementEdgeNodes[elem_code][edge];

      int node_id1, node_id2;

      node_id1 = elem_nodes[elem_edge_node_ids[0]];

      if ( nof_match_nodes == 2 ) {
        node_id2 = elem_nodes[elem_edge_node_ids[1]];
      } else {
        node_id2 = NO_INDEX;
      }

      // Make node_id1 the smaller one (pairs stored only once!)
      if ( node_id2 != NO_INDEX && node_id1 > node_id2) {
        int tmp = node_id1;
        node_id1 = node_id2;
        node_id2 = tmp;
      }

      int edge_nbr;

      bool new_edge = edgeInfos->addConnection(node_id1, node_id2, edge_nbr, true);
      
      edge_id = edge_nbr - 1;

      //--Add edge to the element
      if ( source->edgeIds != NULL ) {
        source->edgeIds[elem_id][edge] = edge_id;
      }

      // If we had an existing edge
      if ( !new_edge ) {
        continue;
      }
        
      // A new edge was found
      nof_edge_elements++;

    } // for each edge in the element

  } // for all source elements

  delete edgeInfos;

}


//---Find element neighbor ids
// NOTE! Neighbor ids simply indices of the table elements (0... nofElements)
//       themselves, so neigbors are from the same table!!!      
// NOTE! Storing neighbor ids works only from uppermost connectivity level: 
// ==> for BULKS via FACES (ie. bulk element neigbors)
// ==> for FACES via EDGES (ie. bulk element neighbors (2D), bndr element neighbors (3D))
// ==> doesn't work for BULKS via EDGES in 3D (connec. can be any finite number!)
// ==> doesn't work for FACES via VERTICES in 2D (connec. can be any finite number!)
void
ModelMeshManager::findMeshElementNeighbors(MeshElementTable* source)
{
  int i, j, k;

  // This structure stores info on node's connections to table elements
  struct nodeInfo {
    int node_id;            // Node id
    short nof_elem;         // Current nof node's parent elements
    short max_nof_elem;     // Max nof of parents elements currenly allocated
    int* elems;             // List of ids to parent elements
  };
  
  //==============================
  //Create Node --> Elements connections
  int nofElements = source->NofElements();
  int nofNodes = meshInfo->nofNodes;

  short alloc_size = 5;
  short realloc_size = 5;

  // Initialize node-info table (Loop all NODES)
  nodeInfo* node_table = new nodeInfo[nofNodes];

  for (i = 0; i < nofNodes; i++) {
    node_table[i].node_id = i;
    node_table[i].nof_elem = 0;
    node_table[i].max_nof_elem = alloc_size;
    node_table[i].elems = new int[alloc_size];
  }
  
  //---Store info on each node's elements (loop all elements)
  for (i = 0; i < nofElements; i++) {

    meshElementCode elem_code = source->getElementCode(i);
    int elem_nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];

    for (int j = 0; j < elem_nof_nodes; j++) {

      // Pick node info entry
      int node_id = source->nodeIds[i][j];
      nodeInfo* ni = &node_table[node_id];

      //--If needed, ALLOCATE more size for elems-table
      if (ni->nof_elem == ni->max_nof_elem) {

        short new_size = ni->nof_elem + realloc_size;
        int* new_elems = new int[new_size];
        //-Copy old info
        for (short k = 0; k < ni->nof_elem; k++)
          new_elems[k] = ni->elems[k];
        //-Update max-size and elems-table pointer
        ni->max_nof_elem = new_size;
        delete[] ni->elems;
        ni->elems = new_elems;
      }

      //--Add new INTERNAL element-id into node elems-table
      ni->elems[ni->nof_elem++] = i;
    }
  }

  //---Sort the elems-lists in each node_table entry 
  for (int nd_index = 0; nd_index < nofNodes; nd_index++) {

    nodeInfo* ni = &node_table[nd_index];
    short last = ni->nof_elem;

    for (i = 0; i < last - 1; i++) {
      short min_pos = i;

      //-Find if there is smaller than i-th element
      for (j = i + 1; j <  last; j++) {

        if (ni->elems[j] < ni->elems[min_pos]) {
          min_pos = j;
        }
      }

      //-If smaller was found, swap i-th and smaller
      if (min_pos != i) {

        int tmp = ni->elems[min_pos];
        ni->elems[min_pos] = ni->elems[i];
         ni->elems[i] = tmp;
      }
    }
  }


  //=================================================
  // Find neighbors for the elements (Loop all ELEMENTS)
  //--Helper tables and variables
  int bndr_node_counter = 0;
  nodeInfo* infos[64];
  short positions[64];
  int check_table[8192];
  int limit = 0x7fffffff; // 'Largest' as a comparison seed
 
  //--Check all sub-elements for each element
  // Basic idea is to pick node-info entries for each node in the
  // boundary and to check if the boundary is referencing also to an other
  // element. In that case the boundary cannot be an outer boundary

  int* sub_node_ids = new int[MAX_NOF_BNDR_NODES];

  // For debugging
  int nof_unknown_ngbr_faces = 0;

  for (int elem_id = 0; elem_id < nofElements; elem_id++) {

    int elem_code = source->getElementCode(elem_id);
    int nof_bndr_elems = MeshElementDesc[elem_code][DESC_NOF_BNDR_ELEMS];
    //int bndr_elem_code = MeshElementDesc[elem_code][DESC_BNDR_ELEM_CODE];

    // Pick each sub-element (face/edge) 
    for (int sub = 0; sub < nof_bndr_elems; sub++) {

      // If neighbor info is already stored (ie. via the neighbor!)
      if ( source->neighborIds[elem_id][sub] != UNSET_INDEX)
        continue;

      int bndr_elem_code = MeshElementBndrCodes[elem_code][sub];

      // Make a table of the nodeInfos for the nodes in the face
      // by copying the entries from *node_table*
      int nof_match_nodes = MeshElementDesc[bndr_elem_code][DESC_NOF_MATCH_NODES];

      for (k = 0; k < nof_match_nodes; k++) {

        int local_id = MeshElementBndrNodes[elem_code][sub][k];
        int node_id = source->nodeIds[elem_id][local_id];
        infos[k] = &node_table[node_id];
        positions[k] = 0;
        sub_node_ids[k] = node_id;
      }
 
      // Sort the element-ids in these infos (kind of merge sort)
      // If we drop the id for the current element from that sorted list
      // and still find a sequence of same element numbers, we know that
      // current boundary beints also to an other element!
      int smallest = limit;
      int box, pos;
      int counter = 0;

      while (1) {
        //-Compare head items in each node's element reference list
        // and select the smallest
        for (k = 0; k < nof_match_nodes; k++) {

          if ( positions[k] < infos[k]->nof_elem &&
                smallest > infos[k]->elems[positions[k]]
             ) {
            smallest = infos[k]->elems[positions[k]];
            box = k;
            pos = positions[k];
          }
        }

        //-If all list are at the end, stop
        if (smallest == limit) {
          break;
        }

        //-Add the smallest into check table if it is not current element!
        int e_id = infos[box]->elems[pos];

        if (e_id != elem_id) {
          check_table[counter++] = e_id;
        }

        //-Update the list position of the smallest and start again
        positions[box] = pos + 1 ;
        smallest = limit;
      } // while

      //-Check if a *continuos* group of same element-ids was found.
      // If group's length is the same as the number of bndr nodes in element
      // these nodes must be also a boundary in that element
      int neighbor_id = NO_INDEX;
      pos = 0;
      short sames = 1; // At least one neighbor (so must be: nof_elems > 1)!

      for (k = 1; k < counter; k++) {
        // Check if still are in the group
        if (check_table[pos] == check_table[k]) {
          sames++;
        } else {
          sames = 1;
          pos = k;
        }
        //-If a group was found break right away
        //NOTE: criterium for a group depends on the boundary element type!
        if (sames == nof_match_nodes) {
          break;
        }
      }
 
      // If group was found, we will get the
      // neighbor id from the check-table
      if (sames == nof_match_nodes) {
        neighbor_id = check_table[pos];
      }
 
      // Set info into the current element
      source->neighborIds[elem_id][sub] = neighbor_id;

      //---Set info in the neighbor element
      if ( neighbor_id != NO_INDEX ) {

        int direction, start_position;
        int ngbr_sub = source->findSubElementIndex(neighbor_id,
                                                   direction, start_position,
                                                   nof_match_nodes, sub_node_ids);
        // Found, ok 
        if (ngbr_sub != NO_INDEX) {
          source->neighborIds[neighbor_id][ngbr_sub] = elem_id;

        // Error!
        } else {
          nof_unknown_ngbr_faces++;
        }

      }

    } // for each sub-element in the element

  } // for elements

#if defined(FRONT_DEBUG)
  // ERROR message
  if (nof_unknown_ngbr_faces > 0) {
    UserInterface* gui = theControlCenter->getGui();
    strstream strm;
    strm << "ERROR: Nof unknown neighbor face indices = "
         << nof_unknown_ngbr_faces
         << ends;

    cerr << "ERROR: Nof unknown neighbor face indices = "
         << nof_unknown_ngbr_faces
         << endl;
    gui->showMsg(strm.str());
  }
#endif

  // Delete work tables
  for (i = 0; i < nofNodes; i++) {
    delete[] node_table[i].elems;
  }
  delete[] node_table;
}


//---Find parent elements for the source table nodes. Set ids in the result
//   tables
// NOTE: nof_nodes should be at least as large as the largest node id in the
// source table!
// This info is needed to calculate average normals for boundary element nodes!
//
void
ModelMeshManager::findMeshElementNodeParents(MeshElementTable* source, int nofNodes, int*& nofNodeParents, int**& nodeParentIds)
{
  MeshConnectionTable* nodeInfos = new MeshConnectionTable(nofNodes, 5, 3, true);

  //--Loop all source ELEMENTS
  for (int elem_id = 0; elem_id < source->NofElements(); elem_id++) {

    meshElementCode elem_code = source->getElementCode(elem_id);

    int nof_elem_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];

    // Get element nodes
    const int* elem_nodes = source->getNodeIds(elem_id);
 
    //--Loop each NODE in the surce element
    for (int nd = 0; nd < nof_elem_nodes; nd++) {

      int node_id = source->nodeIds[elem_id][nd];

      nodeInfos->addConnection(node_id, elem_id, true);

    } // for each node in the element

  } // for all source elements

  // Store parent ids 
  for (int i = 0; i < nofNodes; i++) {
    int nof_parents = nodeInfos->getNofConnections(i);
    nofNodeParents[i] = nof_parents;
    nodeParentIds[i] = new int[nof_parents];

    for (int j = 0; j < nof_parents; j++) {
      nodeParentIds[i][j] = nodeInfos->getParentId(i, j);
    }
  }
  
  delete nodeInfos;

}


// Find selected mesh boundary element index using ray-casting
int
ModelMeshManager::findSelectedMeshBoundaryElement(Renderer* renderer, Point3& ray_start, Point3& ray_dir,
                                       bool try_current_bndr, int& bndr_id,
                                       int& body1_id, int& layer1_id,
                                       int& body2_id, int& layer2_id)
{
  bndr_id = NO_INDEX;
  body1_id = NO_INDEX;
  body2_id = NO_INDEX;
  layer1_id = NO_INDEX;
  layer2_id = NO_INDEX;

  Point3 isec_point;
  BodyElement* be;

  // Selected boundary already given
  // ===============================
  if ( try_current_bndr && modelInfo->selectedBodyElementId != NO_INDEX ) {
    
    be = model->getBoundaryById(modelInfo->selectedBodyElementId);
    
    // Ouch, some error!
    //
    if ( be == NULL ) return NO_INDEX;

    return be->findSelectedMeshElement(renderer, ray_start, ray_dir, isec_point,
                                       bndr_id,
                                       body1_id, layer1_id,
                                       body2_id, layer2_id);
  }

  // Otherwise loop all boundaries
  // =============================

  double min_distance = 1.0e100;
  int max_fem_id = NO_INDEX;
  int be_id;
  int bd1_id, lr1_id;
  int bd2_id, lr2_id;
  
  int index = 0;
  while (true) {
    
    be = model->getBoundary(index++);

    if (be==NULL) break;

    int fem_id = be->findSelectedMeshElement(renderer, ray_start, ray_dir, isec_point,
                                             be_id, bd1_id, lr1_id, bd2_id, lr2_id);

    if ( fem_id != NO_INDEX ) {

      double distance = ray_start[2] - isec_point[2];
      
      distance = (distance < 0 ) ? -distance : distance;

      if ( distance < min_distance ) {

        min_distance = distance;
        max_fem_id = fem_id;

        bndr_id = be_id;
        body1_id = bd1_id;
        body2_id = bd2_id;
        layer1_id = lr1_id;
        layer2_id = lr2_id;
      }
    }
  }

  return max_fem_id;
}

MeshCornerElement*
ModelMeshManager::getMeshCornerElement(int index)
{
  if ( index < 0 ) return NULL;

  MeshCornerElementList::iterator pos = modelData->meshCornerElements->begin();
  MeshCornerElementList::iterator end_pos = modelData->meshCornerElements->end();

  for (int i = 0; i < index && pos != end_pos; pos++);

  if ( pos == end_pos || ++pos == end_pos ) {
    return NULL;
  } else {
    return *pos;
  }
}


// Get all bodies bounding box for the model;
void
ModelMeshManager::getMeshBoundingBox(RangeVector rv) const
{

  meshBox->getRangeVector(rv);
}


int
ModelMeshManager::getMeshBulkElementIdExt2Int(int ext_id)
{
  if (ext_id != NO_INDEX)
    return meshData->bulkElementExt2Int[ext_id];
  else
    return NO_INDEX;
}


int
ModelMeshManager::getMeshInputElementIdExt2Int(int ext_id)
{
  if (ext_id >= 0 && ext_id <= meshInputElementsMaxExtId ) 
    return meshInputElementsExt2Int[ext_id];
  else
    return NO_INDEX;
}


int
ModelMeshManager::getMeshInputElementType(int elem_id)
{
  if ( elem_id < 0 || elem_id >= nofMeshInputElements ) return 0;

  meshElementCode elem_code  = meshInputElements[elem_id].elementCode;

  int elem_type = convertElementCode(elem_code);

  return elem_type;
}


int
ModelMeshManager::getMeshBulkElementIdInt2Ext(int int_id)
{
  if (int_id != NO_INDEX)
    return meshData->bulkElementInt2Ext[int_id];
  else
    return NO_INDEX;
}


int
ModelMeshManager::getMeshNodeIdExt2Int(int ext_id)
{
  if (ext_id != NO_INDEX)
    return meshData->nodeExt2Int[ext_id];
  else
    return NO_INDEX;
}


int
ModelMeshManager::getMeshNodeIdInt2Ext(int int_id)
{
  if (int_id != NO_INDEX)
    return meshData->nodeInt2Ext[int_id];
  else
    return NO_INDEX;
}


// Range vector from the model's mesh bounding box
void
ModelMeshManager::getMeshRangeVector(RangeVector rv) const
{
  meshBox->getRangeVector(rv);
}



void
ModelMeshManager::initClass(Model* mdl)
{
  ModelMeshManager::model = mdl;
}


// Convert input boundary elements to actual boundary elements
Rc
ModelMeshManager::installMeshInputBoundaryElements(bool clear_nodes)
{
  Rc rc;
  int elem_index;

  for (int i = 0; i < nofMeshInputElements; i++) {

    MeshInputElement& te = meshInputElements[i];

    // Add if boundary element
    if ( te.type != 'B' ) continue;
    
    bool is_added;
    elem_index = addMeshBoundaryElement(te.parentTag, te.elementCode, NO_INDEX, NO_INDEX, te.extNodeIds, is_added);
    
    if ( is_added ) 
      te.isAdded = '1';
    else
      te.isAdded = '0';

    if ( elem_index == NO_INDEX ) {
      return ECIF_ERROR;
    }
    
    te.elementId = elem_index;

    if ( clear_nodes ) {
      delete[] te.extNodeIds;
      te.extNodeIds = NULL;
    }
  }

  return ECIF_OK;
}


// Convert input bulk elements to actual bulk elements
Rc
ModelMeshManager::installMeshInputBulkElements(bool clear_nodes)
{
  Rc rc;

  for (int i = 0; i < nofMeshInputElements; i++) {

    MeshInputElement& te = meshInputElements[i];

    // Add if bulk element
    if ( te.type != 'U' ) continue;

    rc = addMeshBulkElement(te.extElementId, te.extParentTag,
                            te.elementCode, te.extNodeIds);

    if (rc != ECIF_OK)
      return rc;
 
    if ( clear_nodes ) {
      delete[] te.extNodeIds;
      te.extNodeIds = NULL;
    }

  }

  return ECIF_OK;
}


// Converts temporary elements (bulk+boundary) to bulk
// and boundary elements
Rc
ModelMeshManager::installMeshInputElements(bool clear_nodes)
{
  Rc rc;
  bool is_added_to_bndr;

  for (int i = 0; i < nofMeshInputElements; i++) {

    MeshInputElement& te = meshInputElements[i];

    // Add element
    if ( te.type == 'U' ) {
      rc = addMeshBulkElement(te.extElementId, te.extParentTag,
                              te.elementCode, te.extNodeIds);
    } else if ( te.type == 'B' )  {
      if ( NO_INDEX == addMeshBoundaryElement(te.parentTag, te.elementCode,
                                              NO_INDEX, NO_INDEX, te.extNodeIds,
                                              is_added_to_bndr)
         ) {
        rc = ECIF_ERROR;
      }
    }

    if (rc != ECIF_OK)
      return rc;

    if ( clear_nodes ) {
      delete[] te.extNodeIds;
      te.extNodeIds = NULL;
    }

  }

  return ECIF_OK;
}


bool
ModelMeshManager::isMeshInputBulkElement(int elem_id)
{
  if ( elem_id < 0 || elem_id >= nofMeshInputElements ) return false;

  return  meshInputElements[elem_id].type == 'U';
}



#if 0
void
ModelMeshManager::normalizeMeshPoints()
{
  return;

  RangeVector rv;
  getRangeVector(rv);
  // Normalize factors;
  double norm[3], shift[3];
  norm[0]   = +1 * (rv[1] - rv[0]) / 2;
  norm[1]   = +1 * (rv[3] - rv[2]) / 2;
  norm[2]   = +1 * (rv[5] - rv[4]) / 2;
  shift[0]  = -1 * (rv[1] + rv[0]) / 2;
  shift[1]  = -1 * (rv[3] + rv[2]) / 2;
  shift[2]  = -1 * (rv[5] + rv[4]) / 2;

  for (int i = 0; i < meshInfo->nofNodes; i++) {
    GcPoint* vp = meshNodeData[i];
    vp->normalize(norm, shift);
  }

}
#endif



// Allocates a table of mesh boundary elements
//
// NOTE: memory for each item is allocated when data is
// actually read, because we dont know yet the nof-node per each
// element.
//
void
ModelMeshManager::reallocateMeshBoundaryElements(int new_size)
{
  // Resize standard mesh table components
  // =====================================
  meshData->boundaryElements->resize(new_size, NULL);

  // Mesh boundary element specific
  // ==============================
  int old_size = meshData->boundaryElements->NofElements();
  int copy_size = (new_size >= old_size)? old_size: new_size;
  int i;

  //--Allocate new tables
  int* tmp1 = new int[new_size]; //  Mesh boundary element boundary ids


  //--Copy old data
  for (i = 0; i < copy_size; i++) {
    tmp1[i] = meshData->boundaryElementBoundaryIds[i];
  }

  //--Init new data
  for(i = old_size; i < new_size; i++) {
    tmp1[i] = NO_INDEX;
  }

  //--Update data
  delete[] meshData->boundaryElementBoundaryIds;
  meshData->boundaryElementBoundaryIds = tmp1;
  meshInfo->nofBoundaryElements = new_size;
}


void
ModelMeshManager::removeMeshInputElements()
{
  for (int i = 0; i < nofMeshInputElements; i++) {
    delete[] meshInputElements[i].extNodeIds;
    meshInputElements[i].extNodeIds = NULL;
  }

  delete[] meshInputElements;
  meshInputElements = NULL;

  delete[] meshInputElementsExt2Int;
  meshInputElementsMaxExtId = 0;

  nofMeshInputElements = 0;
  nofAllocMeshInputElements = 0;
}


void
ModelMeshManager::resetMeshEdgesSelected()
{
#if 0
  if (meshData->bulkEdges != NULL)
    meshData->bulkEdges->resetSelected();
#endif

  if (meshData->boundaryEdges!= NULL)
    meshData->boundaryEdges->resetSelected();

  if (meshData->boundaryVertices!= NULL)
    meshData->boundaryVertices->resetSelected();

}


void
ModelMeshManager::resetMeshRendered()
{
  if (meshData->bulkElements != NULL)
    meshData->bulkElements->resetRendered();

#if 0
  if (meshData->bulkEdges != NULL)
    meshData->bulkEdges->resetRendered();
#endif

  if (meshData->boundaryElements != NULL)
    meshData->boundaryElements->resetRendered();

  if (meshData->boundaryEdges != NULL)
    meshData->boundaryEdges->resetRendered();

  if (meshData->boundaryVertices != NULL)
    meshData->boundaryVertices->resetRendered();

}


void
ModelMeshManager::resetMeshSelected()
{
  if (meshData->bulkElements != NULL)
    meshData->bulkElements->resetSelected();

#if 0
  if (meshData->bulkEdges != NULL)
    meshData->bulkEdges->resetSelected();
#endif

  if (meshData->boundaryElements != NULL)
    meshData->boundaryElements->resetSelected();

  if (meshData->boundaryEdges!= NULL)
    meshData->boundaryEdges->resetSelected();

  if (meshData->boundaryVertices!= NULL)
    meshData->boundaryVertices->resetSelected();

}


void
ModelMeshManager::selectMeshBoundaryElement(int fem_id)
{
  short current_level = meshData->boundaryElements->currentActionLevel;

  if ( model->getFlagValue(SELECT_MODE_TOGGLE) ) {
    meshData->boundaryElements->selected[fem_id] = !meshData->boundaryElements->selected[fem_id];
  }

  else if ( model->getFlagValue(SELECT_MODE_EXTEND) ) {
    if ( !meshData->boundaryElements->selected[fem_id] ) {
      meshData->boundaryElements->selected[fem_id] = true;
      meshData->boundaryElements->actionLevels[fem_id] = current_level;
    }
  }

  else if ( model->getFlagValue(SELECT_MODE_REDUCE) ) {
    if ( meshData->boundaryElements->selected[fem_id] ) {
      meshData->boundaryElements->selected[fem_id] = false;
      meshData->boundaryElements->actionLevels[fem_id] = -current_level;
    }
  }

}


void
ModelMeshManager::unselectMeshBoundaryElement(int fem_id)
{
  short current_level = meshData->boundaryElements->currentActionLevel;

  if ( model->getFlagValue(SELECT_MODE_TOGGLE) ) {
    meshData->boundaryElements->selected[fem_id] = !meshData->boundaryElements->selected[fem_id];
  }

  else if ( model->getFlagValue(SELECT_MODE_EXTEND) ) {
    if ( meshData->boundaryElements->selected[fem_id] ) {
      meshData->boundaryElements->selected[fem_id] = false;
      meshData->boundaryElements->actionLevels[fem_id] = -current_level;
    }
  }

  else if ( model->getFlagValue(SELECT_MODE_REDUCE) ) {
    if ( !meshData->boundaryElements->selected[fem_id] ) {
      meshData->boundaryElements->selected[fem_id] = true;
      meshData->boundaryElements->actionLevels[fem_id] = current_level;
    }
  }

}


void
ModelMeshManager::selectMeshBoundaryElements()
{
  MeshElementTable* table = meshData->boundaryElements;
  table->updateActionLevel(1);

  if ( model->getFlagValue(SELECT_METHOD_ALL) ) {
    selectMeshBoundaryElementsAll();
  }

  else if ( model->getFlagValue(SELECT_METHOD_BY_NEIGHBOR) ) {
    selectMeshBoundaryElementsByNeighbor();
  }

  else if ( model->getFlagValue(SELECT_METHOD_BY_NORMAL) ) {
    selectMeshBoundaryElementsByNormal();
  }

  else if ( model->getFlagValue(SELECT_METHOD_BY_PLANE) ) {
    selectMeshBoundaryElementsByPlane();
  }
}


void
ModelMeshManager::selectMeshBoundaryElementsAll()
{
  BodyElement* be = model->getBoundaryById(modelInfo->selectedBodyElementId);

  if (be == NULL)
    return;

  int nof_elements = be->getNofMeshElements();

  for (int i = 0; i < nof_elements; i++) {

    int index = be->getMeshElementId(i);
    selectMeshBoundaryElement(index);
  }

  model->refreshRenderer();
}


void
ModelMeshManager::selectMeshBoundaryElementsByNeighbor()
{
  if (modelInfo->selectedBodyElementId == NO_INDEX)
    return;

  if (meshInfo->selectedBndrElementId == NO_INDEX)
    return;
  

  BodyElement* be = model->getBoundaryById(modelInfo->selectedBodyElementId);
  if (be == NULL)
    return;

  MeshElementTable* table = meshData->boundaryElements;

  int ref_id = meshInfo->selectedBndrElementId;
  int ref_bndr_id = meshData->boundaryElementBoundaryIds[ref_id];

  double normal_tol = cos(TWO_PI * NORMAL_TOLERANCE / 360);

  double next_active_dgr;
  double next_active_tol = 0.0;
  int next_active_id = NO_INDEX;

  Ids2 id_pair;

  //Ids2Stack pair_stack;
  std::stack<Ids2> pair_stack;

  id_pair.id1 = ref_id;
  id_pair.id2 = NO_INDEX;

  // Init stack
  pair_stack.push(id_pair);

  while (!pair_stack.empty() ) {
    
    // Pick topmost from the stack
    Ids2 ids = pair_stack.top();
    pair_stack.pop();

    int id1 = ids.id1;
    int id2 = ids.id2;

    // Check neighbors if second id not
    // yet set
    if ( id2 == NO_INDEX ) {

      int elem_code = table->getElementCode(id1);
      int nof_sub_elems = MeshElementDesc[elem_code][DESC_NOF_BNDR_ELEMS];

      for (short i = 0; i < nof_sub_elems; i++) {

        int ngbr_id = table->neighborIds[id1][i];
        int ngbr_bndr_id = meshData->boundaryElementBoundaryIds[ngbr_id];

        // If neighbor beints to the boundary and
        // is not yet checked, accept it to the pair
        if ( !table->checked[ngbr_id] &&
             ngbr_bndr_id == ref_bndr_id
           ) {
          id_pair.id2 = ngbr_id;
          pair_stack.push(id_pair);
        }
      } // for sub-elements

      table->checked[id1] = true;
      continue;
    }

    double* ref_normal = table->normals[id1];
    double* normal = table->normals[id2];
    double normal_diff =  dot3(ref_normal, normal);

    // If change between the normals is small enough,
    // select the second and make it a new "point"
    // for neighbor selecting
    if ( normal_diff >= normal_tol ) {
      selectMeshBoundaryElement(id2);
      id_pair.id1 = id2;
      id_pair.id2 = NO_INDEX;
      pair_stack.push(id_pair);
    }
    // Otherwise the second is not useful
    else {
      if ( normal_diff > next_active_tol ) {
        next_active_tol = normal_diff;
        next_active_dgr = ( 360 * acos(next_active_tol) ) / TWO_PI;
        next_active_id = id2;
      }
      table->checked[id2] = true;
    }

  } // while stack not empty

  table->resetChecked();

  UserInterface* gui = theControlCenter->getGui();
  gui->updateNextActiveSelectionTolerance(.001 * int(1000 * next_active_dgr));

  model->refreshRenderer();
}


void
ModelMeshManager::selectMeshBoundaryElementsByNormal()
{
  if (modelInfo->selectedBodyElementId == NO_INDEX)
    return;

  if (meshInfo->selectedBndrElementId == NO_INDEX)
    return;
  

  BodyElement* be = model->getBoundaryById(modelInfo->selectedBodyElementId);
  if (be == NULL)
    return;

  MeshElementTable* table = meshData->boundaryElements;

  int ref_id = meshInfo->selectedBndrElementId;
  double* ref_normal = table->normals[ref_id];

  double normal_tol = cos(TWO_PI * NORMAL_TOLERANCE / 360);

  double next_active_dgr;
  double next_active_tol = 0.0;
  int next_active_id = NO_INDEX;

  // Loop all elements in the boundary
  int nof_elements = be->getNofMeshElements();
  for (int i = 0; i < nof_elements; i++) {

    int index = be->getMeshElementId(i);

    // Don't compare to myself!
    if (index == ref_id)
      continue;

    // Pick the the normal to compare and
    // calculate dot-product
    double* normal = table->normals[index];

    double normal_diff =  dot3(ref_normal, normal);

    bool same_direction = false;

    if ( normal_diff > normal_tol ) {
        same_direction = true;
    } else {
      if ( normal_diff > next_active_tol ) {
        next_active_tol = normal_diff;
        next_active_dgr = ( 360 * acos(next_active_tol) ) / TWO_PI;
        next_active_id = index;
      }
    }

    if (same_direction)
      selectMeshBoundaryElement(index);

  }

  UserInterface* gui = theControlCenter->getGui();
  gui->updateNextActiveSelectionTolerance(next_active_dgr);

  model->refreshRenderer();
}


void
ModelMeshManager::selectMeshBoundaryElementsByPlane()
{
  if (modelInfo->selectedBodyElementId == NO_INDEX)
    return;

  if (meshInfo->selectedBndrElementId == NO_INDEX)
    return;
  

  BodyElement* be = model->getBoundaryById(modelInfo->selectedBodyElementId);
  if (be == NULL)
    return;

  MeshElementTable* table = meshData->boundaryElements;

  int ref_id = meshInfo->selectedBndrElementId;
  double ref_distance = table->normalDistances[ref_id];
  double* ref_normal = table->normals[ref_id];

  double normal_tol = cos(NORMAL_TOLERANCE);
  double distance_tol = meshInfo->dimAvg * DISTANCE_TOLERANCE;

  double next_active_dst;
  double next_active_tol = 0.0;
  int next_active_id = NO_INDEX;

  // Loop all elements in the boundary
  int nof_elements = be->getNofMeshElements();
  for (int i = 0; i < nof_elements; i++) {

    int index = be->getMeshElementId(i);

    // Don't compare to myself!
    if (index == ref_id)
      continue;

    // Pick values to compare
    double distance = table->normalDistances[index];
    double* normal = table->normals[index];

    double distance_diff = ref_distance - distance;
    double normal_diff =  dot3(ref_normal, normal);

    if (distance_diff < 0)
      distance_diff *= -1;

    bool same_direction = false;

    // Compare normals and distances
    if ( normal_diff > normal_tol  &&
         distance_diff <= distance_tol
       ) {
     same_direction = true;
    } else {
      if (distance_diff < next_active_tol) {
        next_active_tol = distance_diff;

        if ( meshInfo->dimAvg > 0 )
          next_active_dst = next_active_tol / meshInfo->dimAvg;
        else
          next_active_dst = 0.0;

        next_active_id = index;
      }
    }

    if (same_direction)
      selectMeshBoundaryElement(index);

  }

  UserInterface* gui = theControlCenter->getGui();
  gui->updateNextActiveSelectionTolerance(next_active_dst);

  model->refreshRenderer();
}


void
ModelMeshManager::selectMeshBoundaryElementsRedo()
{
  BodyElement* be = model->getBoundaryById(modelInfo->selectedBodyElementId);

  if (be == NULL)
    return;

  int nof_elements = be->getNofMeshElements();

  MeshElementTable* table = meshData->boundaryElements;

  short current_level = 0;

  if ( table->currentActionLevel > 0 &&
       table->currentActionLevel < table->maxActionLevel
     )
    current_level = ++table->currentActionLevel;

  for (int i = 0; i < nof_elements; i++) {
    int index = be->getMeshElementId(i);
    short level = table->actionLevels[index];

    if (level == 0 ) 
      continue;

    bool is_neg = (level < 0);

    if (is_neg)
      level *= -1;

    if (level == current_level) {

      if (is_neg)
        table->selected[index] = false;
      else
        table->selected[index] = true;
    }
  }

  model->refreshRenderer();
}


void
ModelMeshManager::selectMeshBoundaryElementsUndo()
{
  BodyElement* be = model->getBoundaryById(modelInfo->selectedBodyElementId);

  if (be == NULL)
    return;

  int nof_elements = be->getNofMeshElements();

  MeshElementTable* table = meshData->boundaryElements;

  short current_level = 0;

  if (table->currentActionLevel > 0)
    current_level = table->currentActionLevel--;

  for (int i = 0; i < nof_elements; i++) {
    int index = be->getMeshElementId(i);
    short level = table->actionLevels[index];

    if (level == 0 ) 
      continue;

    bool is_neg = (level < 0);

    if (is_neg)
      level *= -1;

    if (level == current_level) {

      if (is_neg)
        table->selected[index] = true;
      else
        table->selected[index] = false;
    }
  }

  model->refreshRenderer();
}


// Set an existence flag for a mesh body (material)
// NOTE: This is needed mainly when element's material
// is not known when element is added to model
// (with NO_INDEX material-id). Later body-id is read
// from element set etc.
void
ModelMeshManager::setMeshBodyExt2IntFlag(int external_id)
{
  if (external_id < 0 || external_id >  MAX_NOF_BODIES)
    return;

  meshData->bodyExt2Int[external_id] = 1;
}


void
ModelMeshManager::setMeshBulkElementParentId(int elem_id, int parent_id)
{
  if ( elem_id < 0 || elem_id >= meshInfo->nofBulkElements ) return;

  meshData->bulkElements->parentIds[elem_id][0] = parent_id;
}


void
ModelMeshManager::setMeshData(MeshData* mesh_data, MeshInfo* mesh_info, BoundBox* mesh_box)
{
  meshData = mesh_data;
  meshInfo = mesh_info;
  meshBox = mesh_box;
}


void
ModelMeshManager::setMeshInputElementIsAdded(int elem_id, bool value)
{
  if ( elem_id < 0 || elem_id >= nofMeshInputElements ) return;

  meshInputElements[elem_id].isAdded = value;
}


void
ModelMeshManager::setMeshInputElementParentTag(int elem_id, int parent_tag)
{
  if ( elem_id < 0 || elem_id >= nofMeshInputElements ) return;

  meshInputElements[elem_id].parentTag = parent_tag;
}


void
ModelMeshManager::setMeshInputElementExtParentTag(int elem_id, int ext_parent_tag)
{
  if ( elem_id < 0 || elem_id >= nofMeshInputElements ) return;

  meshInputElements[elem_id].extParentTag = ext_parent_tag;
}


void
ModelMeshManager::setMeshNodes()
{
  MeshElementTable::setMeshNodes(meshData->nodes);
}


void
ModelMeshManager::setModelData(ModelData* model_data, ModelInfo* model_info)
{
  modelData = model_data;
  modelInfo = model_info;
}


// Method splits one mesh corner element of type MEC_303 or MEC_504
// by adding a new node in the center and then forming
// new elements by connecting original edges/faces with the new node.
//
void
ModelMeshManager::splitMeshCornerElement(MeshCornerElement* mce)
{
  if ( mce == NULL ) {
    return;
  }

  MeshElementTable* bt = meshData->bulkElements;

  meshElementCode bulk_code = bt->getElementCode(mce->elementId);

  if ( bulk_code != MEC_303 && bulk_code != MEC_504 ) {
    return;
  }

  int nof_bulk_nodes = MeshElementDesc[bulk_code][DESC_NOF_NODES];
  int nof_bulk_edges = MeshElementDesc[bulk_code][DESC_NOF_EDGES];
  int nof_bulk_faces = MeshElementDesc[bulk_code][DESC_NOF_BNDR_ELEMS];

  int i,j;
 
  // These are for handling new objects
  int int_id, ext_id;
  int new_bulk_id, new_edge_id, new_node_id;

  // Work tables (max size for MEC_504)
  int bulk_edge_node_ids[6][2];
  int parent_ids[2];

  int new_bulk_node_ids[4];
  int new_bulk_edge_ids[6];


  // Create a new center node
  // ========================
  int_id = meshInfo->nofNodes++;
  ext_id = ++meshInfo->maxExternalNodeId;

  meshData->nodeExt2Int[ext_id] = int_id;
  meshData->nodeInt2Ext[int_id] = ext_id;

  new_node_id = int_id;
  mce->newNodeId = int_id;
  meshData->nodes[int_id][0] = mce->centerPoint[0];
  meshData->nodes[int_id][1] = mce->centerPoint[1];
  meshData->nodes[int_id][2] = mce->centerPoint[2];

  const int* bulk_node_ids = bt->getNodeIds(mce->elementId);
  const int* bulk_edge_ids = bt->getEdgeIds(mce->elementId);

  // Collect bulk edge node ids for identifying the new bulk edges
  for (int edge = 0; edge < nof_bulk_edges; edge++) {

    const int* edge_nodes = MeshElementEdgeNodes[bulk_code][edge];

    bulk_edge_node_ids[edge][0] = bulk_node_ids[edge_nodes[0]];
    bulk_edge_node_ids[edge][1] = bulk_node_ids[edge_nodes[1]];
  }

  int first_new_edge_id = meshInfo->nofBulkEdges;

  meshInfo->nofBulkEdges += nof_bulk_faces;

  // Create new bulk element at each face
  // ====================================
  for (int face = 0; face < nof_bulk_faces; face++) {

    parent_ids[0] = mce->bodyId;
    parent_ids[1] = NO_INDEX;

    const int* face_nodes = MeshElementBndrNodes[bulk_code][face];

    // Store new bulk nodes
    for (j = 0 ; j < (nof_bulk_nodes - 1); j++) {
      new_bulk_node_ids[j] = bulk_node_ids[face_nodes[j]];
    }
    new_bulk_node_ids[nof_bulk_nodes - 1] = new_node_id;

    int_id = meshInfo->nofBulkElements++;
    ext_id = ++meshInfo->maxExternalElementId;

    new_bulk_id = int_id;
    meshData->bulkElementExt2Int[ext_id] = int_id;
    meshData->bulkElementInt2Ext[int_id] = ext_id;

    // Add bulk element
    bt->createTableEntry(int_id, bulk_code, parent_ids, nof_bulk_nodes, new_bulk_node_ids);
    bt->nofElements++;

    // Find edge ids for the new bulk element
    for (int new_edge = 0; new_edge < nof_bulk_edges; new_edge++) {

      new_bulk_edge_ids[new_edge] = NO_INDEX;

      const int* edge_nodes = MeshElementEdgeNodes[bulk_code][new_edge];

      int n_nd1_id = new_bulk_node_ids[edge_nodes[0]];
      int n_nd2_id = new_bulk_node_ids[edge_nodes[1]];

      // Try first if new-bulk-edge matches with
      // some of the old bulk edges
      for ( int old_edge = 0; old_edge < nof_bulk_edges; old_edge++) {

        int o_nd1_id = bulk_edge_node_ids[old_edge][0];
        int o_nd2_id = bulk_edge_node_ids[old_edge][1];

        if ( (n_nd1_id == o_nd1_id && n_nd2_id == o_nd2_id) ||
             (n_nd1_id == o_nd2_id && n_nd2_id == o_nd1_id)
           ) {
          new_bulk_edge_ids[new_edge] = bulk_edge_ids[old_edge];
          break;
        }
      }

      if ( new_bulk_edge_ids[new_edge] != NO_INDEX ) {
        continue;
      }

      // Next try if new-bulk-edge matches with
      // some of the new edgesnew edges (created by bulk-nodes with new center node)
      for (int old_node = 0; old_node < nof_bulk_nodes; old_node++) {

        int o_nd1_id = bulk_node_ids[old_node];
        int o_nd2_id = new_node_id;

        if ( (n_nd1_id == o_nd1_id && n_nd2_id == o_nd2_id) ||
             (n_nd1_id == o_nd2_id && n_nd2_id == o_nd1_id)
           ) {
          // NOTE: Node index is implicitely used for indexing new edges
          new_bulk_edge_ids[new_edge] = first_new_edge_id + old_node;
          break;
        }
      }

    } // for each new edge

    // Check if this is a boundary face, then we must
    // update the parent in the boundary mesh element!
    int bndr_elem_id = mce->boundaryElementIds[face];

    if ( bndr_elem_id != NO_INDEX ) {

      int parent_index = meshData->boundaryElements->updateParentId(bndr_elem_id, mce->elementId, int_id);

      int direction = 1;
      meshData->boundaryElements->dirInParents[bndr_elem_id][parent_index] = direction;
    }

    // Set edge ids for the new
    bt->setEdgeIds(new_bulk_id, nof_bulk_edges, new_bulk_edge_ids);

  } // for each bulk face, ie. for each new bulk


  // Net nof new bulk elements created
  // ==================================
  bt->splitted[mce->elementId] = true;
  meshInfo->nofUsedElementTypes[bulk_code] += nof_bulk_faces - 1;
}


