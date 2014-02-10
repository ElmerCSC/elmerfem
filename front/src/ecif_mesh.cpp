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
Module:     ecif_mesh.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  
    
Abstract:   Implementation

************************************************************************/

#include "ecif_mesh.h"
#include "ecif_model.h"
#include "ecif_renderer.h"

//Initialize static class variables.
Model* MeshElementTable::model = NULL;
Point3* MeshElementTable::meshNodes = NULL;
const MeshData* MeshElementTable::meshData = NULL;
const MeshInfo* MeshElementTable::meshInfo = NULL;

int MeshElementTable::elementNodesBuffer[27];
int MeshBulkElementTable::elementNodesBuffer[27];
int MeshFaceElementTable::elementNodesBuffer[27];
int MeshEdgeElementTable::elementNodesBuffer[27];
int MeshVertexElementTable::elementNodesBuffer[27];
double MeshEdgeElementTable::pickingTolerance = 0.0;

const short MAX_MESH_ACTION_LEVEL = 127;

// *** Mesh general element methods *****

MeshElementTable::MeshElementTable(int size)
{
  nofElements = size;

  currentActionLevel = 0;
  maxActionLevel = 0;
  actionLevels = NULL;

  centers = NULL;
  dirInParents = NULL;
  edgeIds = NULL;
  elementCodes = NULL;
  nodeIds = NULL;
  neighborIds = NULL;
  normalDistances = NULL;
  normals = NULL;
  nodeNormals = NULL;
  parentIds = NULL;
  rSquares = NULL;

  checked = NULL;
  rendered = NULL;
  selected = NULL;
  splitted = NULL;

  if (nofElements <= 0) {
    return;
  }

  elementCodes = new meshElementCode[nofElements];
  nodeIds = new int*[nofElements];
  edgeIds = new int*[nofElements];

  // Init data
  for (int i = 0; i < nofElements; i++) {
    elementCodes[i] = MEC_000;
    nodeIds[i] = NULL;
    edgeIds[i] = NULL;
  }
}


MeshElementTable::~MeshElementTable()
{

  for (int i = 0; i < nofElements; i++) {

    // These arries may have subcomponents allocated
    // dynamically with new-operator
    //
    if (edgeIds != NULL && edgeIds[i] != NULL)
      delete[] edgeIds[i];

    if (neighborIds != NULL && neighborIds[i] != NULL)
      delete[] neighborIds[i];

    if (nodeIds != NULL && nodeIds[i] != NULL)
      delete[] nodeIds[i];

    if (nodeNormals != NULL && nodeNormals[i] != NULL)
      delete[] nodeNormals[i];

    if (parentIds != NULL) 
      delete[] parentIds[i];

  }

  // Table itself
  delete[] edgeIds;
  delete[] neighborIds;
  delete[] nodeIds;
  delete[] nodeNormals;
  delete[] parentIds;

  // These arrays do not have dynamically allocated
  // components, so it is enough to delete only the table
  // itself
  delete[] actionLevels;
  delete[] centers;
  delete[] checked;
  delete[] dirInParents;
  delete[] elementCodes;
  delete[] normalDistances;
  delete[] normals;
  delete[] rendered;
  delete[] rSquares;
  delete[] selected;
  delete[] splitted;
}
 

void
MeshElementTable::calcCenter(int index)
{
  if (centers == NULL)
    return;

  const int* nodeIds = getNodeIds(index);

  // Create center point
  Point3& center = centers[index];

  centerPoint(index, center);

  // Calculate square of the distance of the first node
  // from the center
  Point3 r;
  diff3(meshNodes[nodeIds[0]], center, r);
  
  double rsquared = dot3(r, r);
  rSquares[index] = rsquared;
}


void
MeshElementTable::calcNormalDistance(int index)
{
  if ( normalDistances == NULL ||
       centers == NULL          ||
       normals == NULL
     )
    return;

  normalDistances[index] = dot3(normals[index], centers[index]);
}


void
MeshElementTable::centerPoint(int index, Point3& center)
{
  const int* nodeIds = getNodeIds(index);

  center[0] = 0;
  center[1] = 0;
  center[2] = 0;

  meshElementCode elem_code = getElementCode(index);
  short nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];

  for (short i = 0; i < nof_nodes; i++) {
    Point3& p = meshNodes[nodeIds[i]];
    center[0] += p[0];
    center[1] += p[1];
    center[2] += p[2];
  }

  center[0] /= nof_nodes;
  center[1] /= nof_nodes;
  center[2] /= nof_nodes;
}


void
MeshElementTable::createTableEntry(int int_id, meshElementCode elem_code,
                                    const int* parent_pair_ids,
                                    int nof_nodes, const int* entry_node_ids)
{
  int i;
  elementCodes[int_id] = elem_code;

  // Edge ids
  if ( edgeIds != NULL ) {

    int nof_edges = MeshElementDesc[elem_code][DESC_NOF_EDGES];
    edgeIds[int_id] = new int[nof_edges];

    for (i = 0; i < nof_edges; i++) {
      edgeIds[int_id][i] = NO_INDEX;
    }
  }

  // Parent ids
  if  (parentIds != NULL && parent_pair_ids != NULL ) {

    if ( parentIds[int_id] == NULL ) {
      parentIds[int_id] = new int[2];
    }

    for (i = 0; i < 2; i++) {
      parentIds[int_id][i]  = parent_pair_ids[i];
    }
  }

  // Node ids
  if ( nodeIds != NULL && nof_nodes > 0 && entry_node_ids != NULL) {

    nodeIds[int_id] = new int[nof_nodes];

    for (i = 0; i < nof_nodes; i++) {
      nodeIds[int_id][i] = entry_node_ids[i];
    }
  }

  // Neighbor element info
  if (neighborIds != NULL) {

    int nof_sub_elems = MeshElementDesc[elem_code][DESC_NOF_BNDR_ELEMS];
    neighborIds[int_id] = new int[nof_sub_elems];

    for (i = 0; i < nof_sub_elems; i++) {
      neighborIds[int_id][i] = UNSET_INDEX;
    }
  }
 
  // Calc center and max radius
  calcCenter(int_id);
}


void
MeshElementTable::calcNormalsAndDistances()
{
  for (int int_id = 0; int_id < nofElements; int_id++) {
    // Calc normal vector (if needed)
    calcNormalVector(int_id);
    // NOTE: This can be calculated only when center and normal
    // has been calculated!
    calcNormalDistance(int_id);
  }
}


void
MeshElementTable::drawEdges(Renderer* renderer, MeshElementTable* edge_table,
                            int index, bool selected)
{
  if ( index < 0 || index >= nofElements ) {
    return;
  }

  if ( edgeIds == NULL || edgeIds[index] == NULL ) {
    return;
  }

  meshElementCode elem_code = getElementCode(index);

  int nof_edges = MeshElementDesc[elem_code][DESC_NOF_EDGES];

  for (int i = 0; i < nof_edges; i++) {

    int edge_id = edgeIds[index][i];
    bool as_selected = selected || edge_table->selected[edge_id];

    edge_table->drawElement(renderer, edge_id, as_selected);
  }

}


// Method finds subelement ("face") index in the parent by
// comparing subelement's node-ids to that of the potential
// parent.
// Returns NO_INDEX if subelement doesn't belong to the parent
int
MeshElementTable::findSubElementIndex(int elem_index,
                                      int& direction, int& start_position,
                                      int nof_sub_nodes, const int* sub_node_ids)
{
  static int my_sub_node_ids[64];

  meshElementCode elem_code = getElementCode(elem_index);
  int my_nof_sub_elems = MeshElementDesc[elem_code][DESC_NOF_BNDR_ELEMS];

  int nof_nodes_to_check = nof_sub_nodes;

  int i, my_pos, my_pos2;

  const int* my_node_ids = getNodeIds(elem_index);

  //---Loop all faces in the element and compare node ids
  for (int my_sub_index = 0; my_sub_index < my_nof_sub_elems; my_sub_index++) {

    // Collect face node ids and check if this COULD be the face
    for (i = 0; i < nof_sub_nodes; i++) {
      int local_id = MeshElementBndrNodes[elem_code][my_sub_index][i];
      my_sub_node_ids[i] = my_node_ids[local_id];

    }

    int start2;

    if ( !idLoopsAreMatching(nof_sub_nodes, nof_nodes_to_check,
                             my_sub_node_ids, sub_node_ids,
                             direction, start_position, start2) ) {
      continue;

    // Match found
    } else {
      return my_sub_index;
    }

  }

  // No match
  return (int) NO_INDEX;
}


// Get edge ids for the mesh element (if defined)
const int*
MeshElementTable::getEdgeIds(int elem_id)
{
  if ( elem_id < 0 || elem_id >= nofElements) {
    return NULL;
  }

  if ( edgeIds == NULL ) {
    return NULL;
  }

  return edgeIds[elem_id];
}


// Get element code for the mesh element (if defined)
meshElementCode
MeshElementTable::getElementCode(int elem_id)
{
  if ( elem_id < 0 || elem_id >= nofElements) {
    return MEC_000;
  }

  if ( elementCodes == NULL ) {
    return MEC_000;
  }

  return elementCodes[elem_id];
}




// Get nodes for the mesh element
const int*
MeshElementTable::getNodeIds(int elem_id)
{
  if ( nodeIds != NULL ) {
    return nodeIds[elem_id];

  } else {
    return NULL;
  }

}


// Method return a table where proper indices of the sub element nodes are
// given in the parent element
const int*
MeshElementTable::getNodeIdsReversed(int elem_index)
{

  meshElementCode elem_code = getElementCode(elem_index);
  int nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];
  int* reversed = (int*)MeshElementReversedNodeIndices[elem_code];

  int* node_ids = nodeIds[elem_index];
  int* node_ids_reversed = getElementNodesBuffer();

  for (int i = 0; i < nof_nodes; i++) {
    node_ids_reversed[i] = node_ids[reversed[i]];
  }

  return node_ids_reversed;
}


int
MeshElementTable::getParentId(int elem_index, short parent_index)
{
  if (elem_index < 0 || elem_index >= nofElements) {
    return NO_INDEX;
  }

  if ( parent_index >= 2 ) {
    return NO_INDEX;
  }

  return parentIds[elem_index][parent_index];
}


// Getproper indices of the sub element nodes in the parent element
// NOTE: subelement's "face" position and its corner's start position in  the
// parent is taken into account!
void
MeshElementTable::getSubElementNodeIndices(int elem_code,
                                           short sub_elem_index, short sub_start_pos,
                                           int* node_index_buffer)
{
  int sub_elem_code = MeshElementBndrCodes[elem_code][sub_elem_index];
  int nof_sub_nodes = MeshElementDesc[sub_elem_code][DESC_NOF_NODES];
  const int* sub_node_indices = MeshElementBndrNodes[elem_code][sub_elem_index];

  int nof1, nof2, start2;
  nof1 = nof2 = 0;
  start2 = 0;

  switch (sub_elem_code) {

  case MEC_203: nof1 = 2;
    break;
  case MEC_304: nof1 = 3;
    break;
  case MEC_306:
  case MEC_307: nof1 = 3; nof2 = 3; start2 = 3;
    break;
  case MEC_405: nof1 = 4;
    break;
  case MEC_408:
  case MEC_409: nof1 = 4; nof2 = 4; start2 = 4;
    break;
  case MEC_820:
  case MEC_827: nof1 = 8; nof2 = 8; start2 = 8;
    break;
  default:
    nof1 = nof_sub_nodes;
    break;
  }

  // First part (corner nodes)
  if (nof1 > 0) {
    for (int i = 0; i < nof1; i++) {
      int pos = (sub_start_pos + i) % nof1;
      node_index_buffer[i] = sub_node_indices[pos];
    }
  }

  // Second part (middle nodes)
  if (nof2 > 0) {
    for (int i = 0; i < nof2; i++) {
      int pos = (sub_start_pos + i) % nof2;
      node_index_buffer[start2 + i] = sub_node_indices[start2 + pos];
    }
  }

  // Possible inner node
  if (nof_sub_nodes > (nof1 + nof2) )
    node_index_buffer[nof_sub_nodes - 1] = sub_node_indices[nof_sub_nodes - 1];
}


void
MeshElementTable::initClass(Model* mdl)
{
  MeshElementTable::model = mdl;
}


bool
MeshElementTable::isSameElement(int index,
                                short& direction, short& start_position,
                                int other_elem_code, int other_nof_nodes, int* other_node_ids)
{
  meshElementCode elem_code = getElementCode(index);

  if (elem_code != other_elem_code)
    return false;

  int nof_nodes = other_nof_nodes;
  const int* my_node_ids = getNodeIds(index);

  int i;

  // First trivial sum check
  int my_sum = 0;
  int other_sum = 0;
  for (i = 0; i < nof_nodes; i++) {
    my_sum += my_node_ids[i];
    other_sum += other_node_ids[i];
  }

  if (my_sum != other_sum)
    return false;
 
  // Now simply compare in order if node ids matches
  start_position = -1;
  int other_start_position;
  for (i = 0; i < nof_nodes; i++) {
    for (int j = i; j < nof_nodes; j++) {
      if (my_node_ids[i] == other_node_ids[j]) {
        start_position = i;
        other_start_position = j;
        break;
      }
    }
    // If common node was found
    if (start_position != -1)
      break;
  }

  if (start_position == -1)
    return false;

  int my_pos, other_pos;

  // Find direction for my nodes
  my_pos = start_position + 1;
  if (my_pos == nof_nodes)
    my_pos = 0;

  other_pos = other_start_position + 1;
  if (other_pos == nof_nodes)
    other_pos = 0;

  // Try first the positive direction
  if ( my_node_ids[my_pos] == other_node_ids[other_pos]) {
    direction = 1;
  // Try next the negative direction
  } else {
    my_pos = start_position - 1;
    if (my_pos < 0)
      my_pos = nof_nodes - 1;

    if ( my_node_ids[my_pos] == other_node_ids[other_pos]) {
        direction = -1;
    // Ok, no match in either direction
    } else {
      return false;
    }
  }

  // Continue checking the match
  my_pos += direction;
  other_pos +=1;

  // Loop enough nodes the match
  for (i = 2; i < nof_nodes; i++) {

    //-check first array edges
    // "dropping" from left
    if (my_pos < 0)
      my_pos = nof_nodes - 1;

    // "dropping" from right
    else if (my_pos == nof_nodes)
      my_pos = 0;

    if (other_pos == nof_nodes)
      other_pos = 0;

    //-if nodes are different, we can stop
    if (my_node_ids[my_pos] != other_node_ids[other_pos]) {
      return false;
    }
    //-next check positions
    my_pos += direction;
    other_pos += 1;
  }

  return true;
}


void
MeshElementTable::resetActionLevels()
{
  currentActionLevel = 0;
  maxActionLevel = 0;

  if ( actionLevels == NULL ) {
    return;
  }

  for (int i = 0; i < nofElements; i++) {
    actionLevels[i] = 0;
  }
}


void
MeshElementTable::resetState(bool* target)
{
  if ( target == NULL ) {
    return;
  }

  for (int i = 0; i < nofElements; i++) {
    target[i] = false;
  }
}


void
MeshElementTable::resetChecked()
{
  resetState(checked);
}


void
MeshElementTable::resetRendered()
{
  resetState(rendered);
}


void
MeshElementTable::resetSelected()
{
  resetState(selected);
  resetActionLevels();
}


void
MeshElementTable::reverseElementNodes(int index)
{
  if ( nodeIds == NULL ||
       nodeIds[index] == NULL
     )
    return;

  meshElementCode elem_code = getElementCode(index);
  int nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];
  int* nodes = nodeIds[index];
  
  for (int i = 1; i < nof_nodes / 2; i++) {
    long tmp = nodes[i];
    nodes[i] = nodes[nof_nodes - i];
    nodes[nof_nodes - i] = tmp;
  }
}


bool
MeshElementTable::resize(int new_size, bool* copy_flags)
{
  resize_component_impl1((int*&)elementCodes, nofElements, new_size, (int)0, copy_flags);

  resize_component_impl3(edgeIds, nofElements, new_size, 0, (int)0, copy_flags);
  resize_component_impl3(neighborIds, nofElements, new_size, 0, (int)0, copy_flags);
  resize_component_impl3(nodeIds, nofElements, new_size, 0, (int)0, copy_flags);
  resize_component_impl3(parentIds, nofElements, new_size, 0, (int)0, copy_flags);
  resize_component_impl3(dirInParents, nofElements, new_size, 2, (short)0, copy_flags);

  resize_component_impl2((Point3*&)centers, nofElements, new_size, 3, 0.0, copy_flags);
  resize_component_impl2((Point3*&)normals, nofElements, new_size, 3, 0.0, copy_flags);
  resize_component_impl1(normalDistances, nofElements, new_size, 0.0, copy_flags);
  resize_component_impl1(rSquares, nofElements, new_size, 0.0, copy_flags);

  resize_component_impl1(actionLevels, nofElements, new_size, '\0', copy_flags);

  resize_component_impl1(checked, nofElements, new_size, false, copy_flags);
  resize_component_impl1(rendered, nofElements, new_size, false, copy_flags);
  resize_component_impl1(selected, nofElements, new_size, false, copy_flags);
  resize_component_impl1(splitted, nofElements, new_size, false, copy_flags);

  // Finally resize 'nodeNormals' explicitely (sorry, but do not know how to handle
  // this array with template call as above :-)
  //
  if ( nodeNormals != NULL ) {
    Point3** tmp = new Point3*[new_size];
    for (int i = 0; i < new_size; i++) {
      if ( i < nofElements ) {
        tmp[i] = nodeNormals[i];
      } else {
        tmp[i] = NULL;
      }
    }
    delete[] nodeNormals;
    nodeNormals = tmp;
  }

  nofElements = new_size;

  return true;
}


void 
MeshElementTable::resize_component(char*& component, int old_size, int new_size,
                                   char init_value, bool* copy_flags)
{
  resize_component_impl1(component, old_size, new_size, init_value, copy_flags);
}

void 
MeshElementTable::resize_component(bool*& component, int old_size, int new_size,
                                   bool init_value, bool* copy_flags)
{
  resize_component_impl1(component, old_size, new_size, init_value, copy_flags);
}

void 
MeshElementTable::resize_component(short*& component, int old_size, int new_size,
                                   short init_value, bool* copy_flags)
{
  resize_component_impl1(component, old_size, new_size, init_value, copy_flags);
}

void 
MeshElementTable::resize_component(int*& component, int old_size, int new_size,
                                   int init_value, bool* copy_flags)
{
  resize_component_impl1(component, old_size, new_size, init_value, copy_flags);
}

void 
MeshElementTable::resize_component(double*& component, int old_size, int new_size,
                                   double init_value, bool* copy_flags)
{
  resize_component_impl1(component, old_size, new_size, init_value, copy_flags);
}


void 
MeshElementTable::resize_component(char**& component, int old_size, int new_size,
                                   int nof_values, char init_value, bool* copy_flags)
{
  resize_component_impl3(component, old_size, new_size, nof_values, init_value, copy_flags);
}

void 
MeshElementTable::resize_component(bool**& component, int old_size, int new_size,
                                   int nof_values, bool init_value, bool* copy_flags)
{
  resize_component_impl3(component, old_size, new_size, nof_values, init_value, copy_flags);
}

void 
MeshElementTable::resize_component(short**& component, int old_size, int new_size,
                                   int nof_values, short init_value, bool* copy_flags)
{
  resize_component_impl3(component, old_size, new_size, nof_values, init_value, copy_flags);
}

void 
MeshElementTable::resize_component(int**& component, int old_size, int new_size,
                                   int nof_values, int init_value, bool* copy_flags)
{
  resize_component_impl3(component, old_size, new_size, nof_values, init_value, copy_flags);
}

void 
MeshElementTable::resize_component(double**& component, int old_size, int new_size,
                                   int nof_values, double init_value, bool* copy_flags)
{
  resize_component_impl3(component, old_size, new_size, nof_values, init_value, copy_flags);
}


template <class T> void 
MeshElementTable::resize_component_impl1(T*& component, int old_size, int new_size,
                                        T init_value, bool* copy_flags) {

  int i;

  // Nothing to resize!
  if (component == NULL ) {
    return;
  }

  // Allocate new table
  T* component_tmp = new T[new_size];

  int copy_size = (new_size > old_size) ? old_size: new_size;

  // Copy old component's items, use copy_flags to check
  // item should be copied
  int new_index = 0;
  for (i = 0; i < copy_size; i++) {

    // Check if copied
    if ( copy_flags != NULL && !copy_flags[i]) {
      continue;
    }

    component_tmp[new_index] = component[i];
    new_index++;
  }

  // Init possible new items
  for (i = new_index; i < new_size; i++) {
    component_tmp[i] = init_value;
  }

  // Remove old component
  delete[] component;

  // update pointer
  component = component_tmp;
}


// NOTE: This does not allocate component slots, type-T is suppoused to
// be of proper type!
template <class T, class T1> void 
MeshElementTable::resize_component_impl2(T*& component, int old_size, int new_size,
                                         int nof_values, T1 init_value, bool* copy_flags) {

  int i;

  if ( component == NULL ) {
    return;
  }

  // Allocate new table
  T* component_tmp = new T[new_size];

  int copy_size = (new_size > old_size) ? old_size: new_size;

  // Copy old component's items, use copy_flags to check
  // if item should be copied
  int new_index = 0;
  for (i = 0; i < copy_size; i++) {

    // Check if copied
    if ( copy_flags != NULL && !copy_flags[i]) {
      continue;
    }

    // Copy values
    for (int j = 0; j < nof_values; j++) {
      component_tmp[new_index][j] = component[i][j];
    }

    new_index++;
  }

  // Increase size
  if ( new_size > new_index ) {

    for (i = new_index; i < new_size; i++) {

      // Copy values
      // NOTE: No slot allocation!!!
      for (int j = 0; j < nof_values; j++) {
        component_tmp[i][j] = init_value;
      }
    }

  }

  // Delete old component table (but NOT
  // its items!)
  delete[] component;

  // Update pointer
  component = component_tmp;
}


template <class T> void 
MeshElementTable::resize_component_impl3(T**& component, int old_size, int new_size,
                                         int nof_values, T init_value, bool* copy_flags) {

  int i;

  if ( component == NULL ) {
    return;
  }

  // Allocate new table
  T** component_tmp = new T*[new_size];

  int copy_size = (new_size > old_size) ? old_size: new_size;

  // Copy old component's items, use copy_flags to check
  // if item should be copied
  int new_index = 0;
  for (i = 0; i < copy_size; i++) {

#if 0
    if ( copy_flags != NULL && !copy_flags[i]) {
      delete[] component[i];
      continue;
    }
#endif

    component_tmp[new_index] = component[i];
    new_index++;
  }

  // Increase size
  if ( new_size > new_index ) {

    for (i = new_index; i < new_size; i++) {

      if ( nof_values > 0 ) {
        component_tmp[i] = new T[nof_values];
      } else {
        component_tmp[i] = NULL;
      }

      // Copy values
      for (int j = 0; j < nof_values; j++) {
        component_tmp[i][j] = init_value;
      }
    }

  // Decrease size
  } else {

    for (i = new_index; i < old_size; i++) {
      delete[] component[i];
    }
  }

  // Delete old component table (but NOT
  // its items!)
  delete[] component;

  // Update pointer
  component = component_tmp;
}


void
MeshElementTable::setDirInParents(int index, short dir_in_parent1, short dir_in_parent2)
{
  if ( index < 0 || index >= nofElements ) {
    return;
  }

  if ( dirInParents == NULL || dirInParents[index] == NULL ) {
    return;
  }

  dirInParents[index][0] = dir_in_parent1;
  dirInParents[index][1] = dir_in_parent2;
}


bool
MeshElementTable::setEdgeIds(int index, int nof_ids, int* edge_ids)
{
  if ( index < 0 || index >= nofElements ) {
    return false;
  }

  if ( edgeIds == NULL || edgeIds[index] == NULL ) {
    return false;
  }

  meshElementCode elem_code = getElementCode(index);

  int nof_edges = MeshElementDesc[elem_code][DESC_NOF_EDGES];

  if ( nof_edges != nof_ids ) {
    return false;
  }

  for (int i = 0; i < nof_edges; i++) {

    edgeIds[index][i] = edge_ids[i];
  }

  return true;
}


void
MeshElementTable::setSelected(int index, bool value, bool toggle)
{
  if ( selected == NULL ) {
    return;
  }

  if ( index < 0 || index >= nofElements ) {
    return;
  }

  // If toggle-flag is on, we just convert the current value
  if ( toggle ) {
    selected[index] = !selected[index];

  // Otherwise we use the value given in the argument
  } else {
    selected[index] = value;
  }

}



bool
MeshElementTable::setEntry(int index, meshElementCode elem_code, int nof_ids, int* nodes)
{
  if (index < 0 || index >= nofElements) {
    return false;
  }

  delete[] nodeIds[index];
  nodeIds[index] = new int[nof_ids];

  elementCodes[index] = elem_code;

  for (int i = 0; i < nof_ids; i++) {
    nodeIds[index][i] = nodes[i];
  }

  return true;
}


void
MeshElementTable::setMeshData(MeshData* mesh_data, MeshInfo* mesh_info)
{
  meshData = mesh_data;
  meshInfo = mesh_info;
}



void
MeshElementTable::setParentIds(int index, int* parent_ids)
{
  if (index < 0 || index >= nofElements) {
    return;
  }

  parentIds[index][0] = parent_ids[0];
  parentIds[index][1] = parent_ids[1];
}


short
MeshElementTable::updateActionLevel(short change_in_value)
{
  currentActionLevel += change_in_value;

  if (currentActionLevel < 0)
    currentActionLevel = 0;

  if (currentActionLevel > MAX_MESH_ACTION_LEVEL)
    currentActionLevel = MAX_MESH_ACTION_LEVEL;

  if (currentActionLevel > maxActionLevel)
    maxActionLevel = currentActionLevel;

  return currentActionLevel;
}


int
MeshElementTable::updateParentId(int index, int old_id, int new_id)
{
  if (parentIds == NULL || parentIds[index] == NULL)
    return NO_INDEX;


  if ( parentIds[index][0] == old_id ) {
    parentIds[index][0] = new_id;
    return 0;
  }

  if ( parentIds[index][1] == old_id ) {
    parentIds[index][1] = new_id;
    return 1;
  }

  return NO_INDEX;
}

// *********************************
// *** Mesh BULK element methods ***
// *********************************

MeshBulkElementTable::MeshBulkElementTable(int nof_elems)
:MeshElementTable(nof_elems)

{
  parentIds = new int*[nofElements];
  neighborIds = new int*[nofElements];
  splitted = new bool[nofElements];

  // Init data
  for (int i = 0; i < nofElements; i++) {
    neighborIds[i] = NULL;    
    parentIds[i] = new int[2];
    splitted[i] = false;
  }
}


MeshBulkElementTable::~MeshBulkElementTable()
{
}


void
MeshBulkElementTable::drawEdges(Renderer* renderer, const Point3* node_data,
                                int index, bool selected)
{
  static int node_ids[27];

  if ( index < 0 || index >= nofElements ) {
    return;
  }

  if ( edgeIds == NULL || edgeIds[index] == NULL ) {
    return;
  }

  meshElementCode elem_code = getElementCode(index);
  int nof_edges = MeshElementDesc[elem_code][DESC_NOF_EDGES];

  meshElementCode edge_code = (meshElementCode)MeshElementDesc[elem_code][DESC_EDGE_ELEM_CODE];
  int edge_type = MeshElementDesc[edge_code][DESC_ELEM_TYPE];
  int nof_edge_nodes = MeshElementDesc[edge_code][DESC_NOF_NODES];

  for (int edge = 0; edge < nof_edges; edge++) {

    int edge_id = edgeIds[index][edge];

    if ( meshData->bulkEdgeRendered[edge_id] ) {
      continue;
    }

    meshData->bulkEdgeRendered[edge_id] = true;

    for ( int i = 0; i < nof_edge_nodes; i++) {

      const int* edge_nodes = (int*)MeshElementEdgeNodes[elem_code][edge];
      node_ids[i] = nodeIds[index][edge_nodes[i]];
    }

    renderer->drawMeshElement(edge_type, node_ids, NULL, node_data, 1, selected);
  }

}



#if 0
int
MeshBulkElementTable::findFaceIndex(int elem_index, short direction,
                                    int nof_face_nodes, int* face_node_ids)
{
  static int* my_face_node_ids = new int[MAX_NOF_BNDR_NODES];

  int elem_code = getElementCode(index);
  int bndr_elem_code = MeshElementDesc[elem_code][DESC_BNDR_ELEM_CODE];
  int elem_type = MeshElementDesc[elem_code][DESC_ELEM_CODE];
  int nof_faces = MeshElementDesc[elem_code][DESC_NOF_BNDR_ELEMS];
  int nof_nodes = MeshElementDesc[bndr_elem_code][DESC_NOF_NODES];

  // This should not be possible, but we check it anyway!
  if ( nof_face_nodes != nof_nodes)
    return (int)NO_INDEX;

  //int nof_nodes_to_check = nof_face_nodes;
  int nof_nodes_to_check = 2; // This should be enough for our element types!!!

  int i, my_pos;
  // Loop all faces in the element and compare node ids
  for (int face_index = 0; face_index < nof_faces; face_index++) {
    bool is_same_face = true;
    // Collect face node ids and check if this COULD the face
    my_pos = -1;
    for (i = 0; i < nof_nodes; i++) {
      int local_id = MeshElementBndrDesc[elem_code][face_index][i];
      my_face_node_ids[i] = nodeIds[elem_index][local_id];
      if (my_face_node_ids[i] == face_node_ids[0])
        my_pos = i;
    }

    // If no common node was found
    if (my_pos == -1)
      continue;

    // Loop enough nodes in the faces to check the match
    // external face is looped in order: 0 1 2 ...
    // my face is looped: my_pos +- 1 depending of the direction variable
    for (i = 0; i < nof_nodes_to_check; i++) {
      // if nodes are different, we can stop
      if (my_face_node_ids[my_pos] != face_node_ids[i]) {
        is_same_face = false;
        break;
      }
      // my next check position
      my_pos = (my_pos + direction);
      // "dropping" from left
      if (my_pos < 0)
        my_pos = nof_nodes - 1;
      // "dropping" from right
      else if (my_pos == nof_nodes)
        my_pos = 0;
    }

    // If same face was found
    if (is_same_face) {
      return face_index;
    }
  }

  return (int) NO_INDEX;
}
#endif


meshElementCode
MeshBulkElementTable::getElementCode(int elem_index)
{
  return elementCodes[elem_index];
}


bool
MeshBulkElementTable::isBoundaryFace(int index, int face, bool& is_first_parent) {

  is_first_parent = true;

  if ( index < 0 || index >= nofElements ) {
    return false;
  }
  
  if ( neighborIds == NULL || parentIds == NULL) {
    return false;
  }

  int neighbor_id = neighborIds[index][face];

  //---Check if face is a boundary element
  int parent1_id = parentIds[index][0];
  int parent2_id = NO_INDEX;

  if ( neighbor_id != NO_INDEX) {
    parent2_id = parentIds[neighbor_id][0];
  }

  if ( parent2_id != NO_INDEX && parent1_id > parent2_id ) {
    is_first_parent = false;
  }

  //---If we had a boundary element, update counter
  if ( parent1_id == parent2_id ) {
    return false;
  } else {
    return true;
  }
}


void
MeshBulkElementTable::resetRendered()
{
  resetState(rendered);

  for (int i = 0; i < meshInfo->nofBulkEdges; i++) {
    meshData->bulkEdgeRendered[i] = false;
  }
}



// *** Mesh FACE element methods *****

MeshFaceElementTable::MeshFaceElementTable(int nof_elems, bool store_neighbors)
:MeshElementTable(nof_elems)
{
  if (store_neighbors) {
    neighborIds = new int*[nofElements];
  }

  checked = new bool[nofElements];
  rendered = new bool[nofElements];
  selected = new bool[nofElements];

  // init data
  for (int i = 0; i < nofElements; i++) {

    if (store_neighbors) {
      neighborIds[i] = NULL;
    }

    checked[i] = false;
    rendered[i] = false;
    selected[i] = false;
  }
}


MeshFaceElementTable::~MeshFaceElementTable()
{
}


/*
  Graph.Alg.Faq:
Subject 5.05: How do I find the intersection of a line and a plane?
    If the plane is defined as:
        a*x + b*y + c*z + d = 0
    and the line is defined as:
        x = x1 + (x2 - x1)*t = x1 + i*t
        y = y1 + (y2 - y1)*t = y1 + j*t
        z = z1 + (z2 - z1)*t = z1 + k*t
    Then just substitute these into the plane equation. You end up
    with:
        t = - (a*x1 + b*y1 + c*z1 + d)/(a*i + b*j + c*k)
          When the denominator is zero: 
            -the line is contained in the plane 
             if the numerator is also zero (the point at t=0 satisfies the
             plane equation),
            -otherwise the line is parallel to the plane.
*/


// NOTE: This should work for all face elements (303 - 409), but only
// if mid-edge and middle nodes are after "corner" nodes in the node list!
//
// Return nof intersections
//
int
MeshFaceElementTable::calcLineIntersections(int elem_index, short direction,
                                            Point3& lstart, Point3& ldir, Point3* isec_points)
{

  meshElementCode elem_code = getElementCode(elem_index);

  if ( elem_code >= MEC_404 && elem_code < MEC_504 ) {
    return calcQuadriLineIntersections(elem_index, direction, lstart, ldir, isec_points);

  } else {
    return calcTriangleLineIntersections(elem_index, direction, lstart, ldir, isec_points);
  }
}



// NOTE: This should work for all quadri elments (404 - 409), but only
// if mid-edge and middle nodes are after "corner" nodes in the node list!
// Return nof intersections
//
int
MeshFaceElementTable::calcQuadriLineIntersections(int elem_index, short direction,
                                                  Point3& lstart, Point3& ldir, Point3* isec_points)
{
  // Ccw ordered nodes (if direction is -1, it is from the parent2
  // and nodes are cw oriented and must me reordered
  Point3& p0 = meshNodes[nodeIds[elem_index][0]];
  Point3& p1 = meshNodes[nodeIds[elem_index][1]];
  Point3& p2 = meshNodes[nodeIds[elem_index][2]];
  Point3& p3 = meshNodes[nodeIds[elem_index][3]];
  Point3& normal = normals[elem_index];

  static Point3* points[4];

  points[0] = &p0;
  points[1] = (direction >= 0)? &p1 : &p2;
  points[2] = (direction >= 0)? &p2 : &p1;
  points[3] = &p3;

#if 0
  // If element is looking into wrong direction
  //
  // NOTE: This makes picking much faster, so wew have to use it (unless some better
  // method is found), although it means that elements cannot be selected from 'inside',
  // which would be quite convenient in some cases!
  //
  if ( direction == 1 ) {

    if ( !isLess(dot3(normal, ldir), 0) ) {
      return 0;
    }

  } else if ( direction == -1 ){

    if ( !isGreater(dot3(normal, ldir), 0) ) {
      return 0;
    }
  }
#endif

  // Plane equation for the normal and a point in the plane is;
  // (r) dot (normal) = (r0) dot (normal) = d
  // So for the form Ax + By + Cz + D = 0, we have
  // A = normal[0], B = normal[1], C = normal[2], D = -d

  double D = -1 * dot3(p0, normal);

  double numer = dot3(normal, lstart) + D;
  double denom = dot3(normal, ldir);

  double t;

  // Intersection
  if (denom != 0) {
    t = - numer / denom;

  // Line is on the plane
  } else if (numer == 0) {
    t = 0.0;

  // Line is parallel,but not in the plane
  } else {
    return 0;
  }

  //-Calc intersection point from the line equation
  Point3 tmp;
  scalarmult(t, ldir, tmp);
  Point3 isec_point;
  add3(lstart, tmp, isec_point);

  // Finally check if intersection point
  // is inside the element (quadrilateral)
  if ( pointInsideRectangle(isec_point, points, centers[elem_index], rSquares[elem_index]) ) {

    copy3(isec_point, *isec_points);
    return 1;

  } else {
    return 0;
  }

}


// NOTE: This should work for all triangle elments (303 - 306), but only
// if mid-edge and middle nodes are after "corner" nodes in the node list!
// Return nof intersections
//
int
MeshFaceElementTable::calcTriangleLineIntersections(int elem_index, short direction,
                                                    Point3& lstart, Point3& ldir, Point3* isec_points)
{
  // Ccw ordered nodes (if elem_dir is -1, it is from the parent2
  // and nodes are cw oriented and must me reordered
  Point3& p0 = meshNodes[nodeIds[elem_index][0]];
  Point3& p1 = meshNodes[nodeIds[elem_index][1]];
  Point3& p2 = meshNodes[nodeIds[elem_index][2]];
  Point3& normal = normals[elem_index];

  static Point3* points[3];

  points[0] = (direction >= 0)? &p0 : &p1;
  points[1] = (direction >= 0)? &p1 : &p0;
  points[2] = &p2;

#if 0
  // If element is looking into wrong direction
  //
  // NOTE: This makes picking much faster, so wew have to use it (unless some better
  // method is found), although it means that elements cannot be selected from 'inside',
  // which would be quite convenient in some cases!
  //
  if ( direction == 1 ) {
    if ( !isLess(dot3(normal, ldir), 0) ) {
      return 0;
    }

  } else if ( direction == -1 ) {
    if ( !isGreater(dot3(normal, ldir), 0) ) {
      return 0;
    }
  }
#endif

  // Plane equation for the normal and a point in the plane is;
  // (r) dot (normal) = (r0) dot (normal) = d
  // So for the form Ax + By + Cz + D = 0, we have
  // A = normal[0], B = normal[1], C = normal[2], D = -d

  double D = -1 * dot3(p0, normal);

  double numer = dot3(normal, lstart) + D;
  double denom = dot3(normal, ldir);

  double t;

  // Intersection
  if (denom != 0) {
    t = - numer / denom;

  // Line is on the plane
  } else if (numer == 0) {
    t = 0.0;

  // Line is parallel,but not in the plane
  } else {
    return 0;
  }


  //-Calc intersection point from the line equation
  Point3 tmp;
  scalarmult(t, ldir, tmp);
  Point3 isec_point;
  add3(lstart, tmp, isec_point);

  // Finally check if intersection point
  // is inside the element (triangle)
  if ( pointInsideTriangle(isec_point, points, centers[elem_index], rSquares[elem_index]) ) {

    copy3(isec_point, *isec_points);
    return 1;

  } else {
    return 0;
  }

}



void
MeshFaceElementTable::calcNormalVector(int index)
{
  if (normals == NULL)
    return;

  meshElementCode elem_code = getElementCode(index);

  // Pick normal
  Point3& normal = normals[index];
  Point3 *p1, *p2, *p3;

  switch (elem_code) {
  case MEC_303:
  case MEC_304:
  case MEC_306:
  case MEC_307:
    p1 = &meshNodes[nodeIds[index][0]];
    p2 = &meshNodes[nodeIds[index][1]];
    p3 = &meshNodes[nodeIds[index][2]];
    break;
  case MEC_404:
  case MEC_405:
  case MEC_408:
  case MEC_409:
    p1 = &meshNodes[nodeIds[index][0]];
    p2 = &meshNodes[nodeIds[index][1]];
    p3 = &meshNodes[nodeIds[index][2]];
    break;
  default:
    return;
  }

  Point3 a, b;
  diff3(*p3, *p1, a);
  diff3(*p2, *p1, b);
  cross3(b, a, normal);
  normalize(normal);

  for (int i = 0; i < 3; i++) {
    if ( isZero(normal[i]) ) {
      normal[i] = 0.0;
    }
  }
}


void
MeshFaceElementTable::drawElement(Renderer* renderer, int index,
                                  short direction, bool selected)
{
  if ( rendered[index] ) {
    return;
  }

  rendered[index] = true;

  MeshElementTable* edges = model->getMeshBoundaryElementEdges();

  meshElementCode elem_code = getElementCode(index);
  int elem_type = MeshElementDesc[elem_code][DESC_ELEM_TYPE];

  // Call triangle or some other area drawing functions
  renderer->drawMeshElement(elem_type, nodeIds[index],
                            nodeNormals[index], meshNodes,
                            direction, selected);
}


// ***** Mesh EDGE element methods *****

MeshEdgeElementTable::MeshEdgeElementTable(int nof_elems, bool store_neighbors)
:MeshElementTable(nof_elems)
{
  if (store_neighbors) {
    neighborIds = new int*[nofElements];
  }

  checked = new bool[nofElements];
  rendered = new bool[nofElements];
  selected = new bool[nofElements];

  // init data
  for (int i = 0; i < nofElements; i++) {

    if (store_neighbors) {
      neighborIds[i] = NULL;
    }

    checked[i] = false;
    rendered[i] = false;
    selected[i] = false;
  }
}


MeshEdgeElementTable::~MeshEdgeElementTable()
{
}


int
MeshEdgeElementTable::calcLineIntersections(int elem_index, short elem_dir,
                                            Point3& lstart, Point3& ldir, Point3* isec_points)
{
  meshElementCode elem_code = getElementCode(elem_index);

  if (elem_code <  MEC_202 || elem_code >= 303)
    return 0;

  const int* nodeIds = getNodeIds(elem_index, elem_dir);

  Point3& p0 = meshNodes[nodeIds[0]];
  Point3& p1 = meshNodes[nodeIds[1]];
  
  Point3& normal = normals[elem_index];
  if (elem_dir == -1)
    scalarmult(-1, normal, normal);

  // Check end-point cases
  if ( samepoint(p0, lstart) ){
    copy3(p0, *isec_points);
    return 1;
  }
  if ( samepoint(p1, lstart) ){
    copy3(p1, *isec_points);
    return 1;
  }


  Point3 edge_dir, l_delta0, tmp;
 
  // Edge direction vector (normalized)
  edge_dir[0] = normal[1];
  edge_dir[1] = -1 * normal[0];
  edge_dir[2] = 0.0;

  // Edge length
  diff3(p1, p0, tmp);
  double edge_len = dot3(edge_dir, tmp);

  // Vector l_delta0 = lstart - p0
  diff3(lstart, p0, l_delta0);

  // Check that intersection is "within" the edge
  // project the intersection point to the edge
  double t = dot3(edge_dir, l_delta0);
  if ( isLess(t, 0.0) ||
       isGreater(t, edge_len)
     )
    return 0;

  // Check that intersection distance from the edge is ok
  // project intersection point to the edge normal
  double d = dot3(normal, l_delta0);
  if (d < 0)
    d *= -1;
  if ( isGreater(d, MeshEdgeElementTable::pickingTolerance) )
    return 0;

  // Intersection point is: p0 + t * (p1 - p0)
  scalarmult(t, edge_dir, tmp);
  add3(p0, tmp, *isec_points);

  return 1;
}


void
MeshEdgeElementTable::calcNormalVector(int index)
{
  if (normals == NULL)
    return;

  // Create normal
  const int* nodeIds = getNodeIds(index);
  Point3& normal = normals[index];

  Point3& p1 = meshNodes[nodeIds[0]];
  Point3& p2 = meshNodes[nodeIds[1]];

  Point3 a;
  diff3(p2, p1, a);
  normal[0] = -a[1];
  normal[1] = a[0];
  normal[2] = 0.0;

  normalize(normal);

  for (int i = 0; i < 3; i++) {
    if ( isZero(normal[i]) ) {
      normal[i] = 0.0;
    }
  }
}


void
MeshEdgeElementTable::draw(Renderer* renderer, bool el_selected)
{
  meshElementCode elem_code;
  int elem_type;

  for (int i = 0; i < nofElements; i++) {

    elem_code = getElementCode(i);
    elem_type = MeshElementDesc[elem_code][DESC_ELEM_TYPE];

    renderer->drawMeshElement(elem_type, getNodeIds(i), NULL, meshNodes,
                              1, el_selected);
  }
}


void
MeshEdgeElementTable::drawElement(Renderer* renderer, int index, bool selected)
{
  drawElement(renderer, index, 1, selected);
}


void
MeshEdgeElementTable::drawElement(Renderer* renderer, int index,
                                  short direction, bool selected)
{
  if ( rendered[index] ) {
    return;
  }
    
  rendered[index] = true;

  meshElementCode elem_code = getElementCode(index);
  int elem_type = MeshElementDesc[elem_code][DESC_ELEM_TYPE];

  renderer->drawMeshElement(elem_type, nodeIds[index], NULL, meshNodes,
                            direction, selected);
}

 
const int*
MeshEdgeElementTable::getNodeIds(int elem_index, short direction)
{
  if (nodeIds != NULL) {

    if (direction == 1)
      return nodeIds[elem_index];
    else
      return getNodeIdsReversed(elem_index);

  } else {
    return NULL;
  }
}



// *** Mesh VERTEX element methods *****

MeshVertexElementTable::MeshVertexElementTable(int nof_elems)
:MeshElementTable(nof_elems)
{
  checked = new bool[nofElements];
  rendered = new bool[nofElements];
  selected = new bool[nofElements];

  // init data
  for (int i = 0; i < nofElements; i++) {
    checked[i] = false;
    rendered[i] = false;
    selected[i] = false;
  }

}


MeshVertexElementTable::~MeshVertexElementTable()
{
}


void
MeshVertexElementTable::drawElement(Renderer* renderer, int index, bool selected)
{
  if ( rendered[index] ) {
    return;
  }
  
  rendered[index] = true;
  
  renderer->drawMeshElement(101, nodeIds[index], NULL, meshNodes, 1, selected);
}


const int*
MeshVertexElementTable::getNodeIds(int elem_index, short direction)
{
  static int nodeIds[1];
  nodeIds[0] = elem_index;

  return nodeIds;
}


