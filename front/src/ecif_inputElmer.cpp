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
Module:     ecif_inputElmer.cpp
Language:   C++
Date:       17.11.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation (read ELmer mesh from DB)

************************************************************************/

#include "eio_api.h"

#include "ecif_body2D.h"
#include "ecif_body3D.h"
#include "ecif_bodyElement.h"
#include "ecif_control.h"
#include "ecif_geometry.h"
#include "ecif_inputElmer.h"
#include "ecif_model.h"
#include "ecif_userinterface.h"
#include "ecif_timer.h"

extern char read_buffer[];


InputElmer::InputElmer(enum ecif_modelDimension m_dim,
                       ifstream& in_file, char* in_filename):
Input(m_dim, in_file, in_filename)
{
  int len = strlen(infileName);

  int pos = len - 1;

  // Remove possible trailing path separator
  while (pos > 0 ) {

    if ( infileName[pos] == '/' ||
         infileName[pos] == '\\'
         ) {

      infileName[pos] = '\0';
      break;
    }

    pos--;
  }
}


#if 0
enum ecif_modelDimension
InputElmer::findMeshModelDimension() {

  int info = 0;

  UserInterface* gui = theControlCenter->getGui();

  ModelInfo* minfo = (ModelInfo*) model->getModelInfo();

  enum ecif_modelDimension dimension = minfo->dimension;

  //---Init EIO
  eio_init(info);

  //---Find model dimension from the mesh element types
  if (dimension == ECIF_ND) {

    int nodeC, elementC, bndrElementC;
    nodeC = elementC = bndrElementC = 0;

    int nof_used_element_types = 0;
    int* used_element_types = new int[MAX_NOF_ELEM_CODES];
    int* used_element_type_counts = new int[MAX_NOF_ELEM_CODES];

    // Read mesh info and close
    eio_open_mesh(infileName, info);

    if (info == -1) {
      gui->errMsg(0, "Cannot read Elmer mesh header in the directory:", infileName);
      return ECIF_ND;
    }

    eio_get_mesh_description(nodeC, elementC, bndrElementC,
                             nof_used_element_types,
                             used_element_types, used_element_type_counts, info);
    eio_close_mesh(info);

    // Problems for reading Elmer mesh!
    if (info == -1 || nodeC <= 0 || elementC <= 0 || bndrElementC <= 0) {

      gui->errMsg(0, "Cannot read Elmer mesh header description in the directory:", infileName);

      return ECIF_ND;
    }

    // Find dimension from element types
    int max_type = 0;
    for (int i = 0; i < nof_used_element_types; i++) {

      if (used_element_types[i] > max_type)
        max_type = used_element_types[i];
    }

    // 2D
    if (max_type >= 200 && max_type <= 500)
      dimension = ECIF_2D;
    // 3D
    else if (max_type >= 500)
      dimension = ECIF_3D;

    delete[] used_element_types;
    delete[] used_element_type_counts;
  }

  //---Close EIO
  eio_close(info);

  return dimension;
}
#endif


enum ecif_modelDimension
InputElmer::findMeshModelDimension() {

  ecif_modelDimension dim = ECIF_ND;
  ecif_modelDimension box_dim = ECIF_ND;
  ecif_modelDimension elem_dim = ECIF_ND;

  RangeVector rv;
  bool x_is_const = false;
  bool y_is_const = false;
  bool z_is_const = false;

  meshBoundingBox.getRangeVector(rv);

  // Find model dimension by boundary box
  if ( isEqual(rv[0], rv[1]) ) {
    x_is_const = true;
  }

  if ( isEqual(rv[2], rv[3]) ) {
    y_is_const = true;
  }

  if ( isEqual(rv[4], rv[5]) ) {
    z_is_const = true;
  }

  if ( !(x_is_const || y_is_const || z_is_const) ) {
    box_dim = ECIF_3D;
  } else if ( !( x_is_const || y_is_const) && z_is_const ) {
    box_dim = ECIF_2D;
  } else if ( !(x_is_const && y_is_const && z_is_const) ) {
    box_dim = ECIF_1D;
  } else {
    box_dim = ECIF_ND;
  }

  // Find actual dimension
  dim = findModelDimension(box_dim);

  return dim;
}


//---Create and update mesh tables
//
// NOTE: Elmer mesh is processed differently compared to other
// 'external' meshes. It has its own 'processMeshFileData'-method!
//
bool
InputElmer::processMeshFileData()
{
  Rc rc;

  const ModelStatistics* ms = model->getModelStatistics();

  //---Elmer mesh as an external mesh!
  //
  if ( ms->nofBodies == 0 ) {
    return Input::processMeshFileData();
  }

  //---Elmer mesh as a model mesh
  //
  Timer timer;
  double time, time1, time2;

  //---We can allocate and create actual bulk and boundary elements elements
  //model->allocateMeshBulkElements(nofInputBulkElements, nofInputBulkElements);
  //rc = model->installMeshInputBulkElements();

  //model->createMeshBodyTables();
  //model->createMeshBodies();
  //model->convertMeshBulkElementIdsExt2Int();

  //model->allocateMeshBoundaryElements(nofInputBoundaryElements);
  //rc = model->installMeshInputBoundaryElements();

  //model->removeMeshInputElements();

  UserInterface* gui = theControlCenter->getGui();

  // Bulk element stuff
  // ==================

  MeshElementTable* bt = model->getMeshBulkElements();

  //---Create BULK element connections (neighor ids, sub element indices
  timer.start(); time1 = 0;
  gui->showMsg("---Creating volume element connections ...");

  model->findMeshElementNeighbors(bt);

  time2 = timer.getLapTime(WALL_TIME); time = time2 - time1; time1 = time2;
  gui->showUsedTimeMsg(time, "---Creating volume element connections", short(0), false);

  //---Create bulk elmement edges
  model->createMeshBulkElementEdges();

  // Boundary element stuff
  // ======================

  MeshElementTable* bet = model->getMeshBoundaryElements();

  //---Create mesh boundary element neighbors
  model->findMeshElementNeighbors(bet);

  //---Create boundary element edges
  model->createMeshBoundaryElementEdges();

  //---Finish
  timer.stop();
  time = timer.getEndTime(WALL_TIME);
  gui->showUsedTimeMsg(time, "Creating mesh bodies and boundaries", 0,  true);

	return modelDimension;
}


// Read mesh (nodes, bulk and bndr elements) from the Elmer DB
//
bool
InputElmer::readMeshGeometry()
{
  int info;
  int bulk_element_count;
  int bndr_element_count;
  int node_count;
  int nof_used_element_types;
  int used_element_types[64];
  int used_element_type_counts[64];

  int input_elem_count = 0;
  int actual_bulk_elem_count = 0;
  int actual_bndr_elem_count = 0;

  int edge_elem_count = 0;
  int vrtx_elem_count = 0;

  UserInterface* gui = theControlCenter->getGui();

  const ModelInfo* mi = model->getModelInfo();
  const ModelStatistics* ms = model->getModelStatistics();

  bool model_has_geometry = (ms->nofElements > 0);

  //---Init DB
  eio_init(info);

  if (info == -1) {
    gui->errMsg(0, "WARNING: Cannot open Elmer mesh files using directory:\n\n", infileName);
    return resetMeshData();
  }

  //---Open DB mesh data
  eio_open_mesh(infileName, info);

  if (info == -1) {
    gui->errMsg(0, "WARNING: Cannot open Elmer mesh files using directory:\n\n", infileName);
    return resetMeshData();
  }

  // Read mesh info
  eio_get_mesh_description(node_count, bulk_element_count, bndr_element_count,
                           nof_used_element_types,
                           used_element_types, used_element_type_counts, info);

  maxExternalElementId = max(bulk_element_count, bndr_element_count);

  //---Delete possible old mesh data and init data
  model->resetMeshData();

  //---Allocate mesh tables
  model->allocateMeshBodies(ms->nofBodies);
  model->allocateMeshNodes(node_count, node_count);

  //---Read nodes from the mesh input file
  if ( ECIF_OK != readMeshNodes(node_count) ) {
    return resetMeshData();
  }
  
  // Set model dimension
  //
  if ( !model_has_geometry ) {
    modelDimension = findMeshModelDimension();
    model->setModelDimension(modelDimension);

  } else {
    modelDimension = inputDimension;
  }

  //---Calculate nof "real" bulk/boundary elements
  short table_index = DESC_ELEM_TYPE; // Description table column index

  for (short i = 0; i < nof_used_element_types; i++) {

    int elem_code = model->convertElementType(used_element_types[i]);

    if ( elem_code == MEC_000 ) {
      return false;
    }

    // Subtract elements which are not actual bulk/boundary elements
    //
    // NOTE: Normally we just drop vertices(2D)/edges(3D) from boundary elements
    // but especially a 2D model read as 3D-BEM model needs special treatment
    //

    // These are bulk elements
    if (
        (modelDimension == ECIF_3D && MeshElementDesc[elem_code][table_index] > 500) ||
        (modelDimension == ECIF_2D && MeshElementDesc[elem_code][table_index] > 300)
       ) {
      actual_bulk_elem_count += used_element_type_counts[i];

    // These are boundary elements
    } else if (
        (modelDimension == ECIF_3D && MeshElementDesc[elem_code][table_index] > 300) ||
        (modelDimension == ECIF_2D && MeshElementDesc[elem_code][table_index] > 200)
       ) {
      actual_bndr_elem_count += used_element_type_counts[i];

    // These are edge elements
    } else if ( MeshElementDesc[elem_code][table_index] > 200 ) {
      edge_elem_count += used_element_type_counts[i];

      // These are vertex elements
    } else {
      vrtx_elem_count += used_element_type_counts[i];
    }
  }

  int nof_input_elements = 0;

  if ( !model_has_geometry ) {
    nof_input_elements += edge_elem_count + vrtx_elem_count;
  }

  if ( bulk_element_count > actual_bulk_elem_count ) {
    nof_input_elements += bulk_element_count - actual_bulk_elem_count;
  }

  if ( nof_input_elements > 0 ) {
    model->allocateMeshInputElements(nof_input_elements, maxExternalElementId);
  }

  // Bulk element table
  // ------------------

  model->allocateMeshBulkElements(actual_bulk_elem_count, maxExternalElementId);
  bulkElementsAllocated = true;

   //---Read volume elements from the mesh input file
  if ( ECIF_OK != readMeshBulkElements() ) {
    return resetMeshData();
  }

  model->setMeshNodes();

  // NOTE: We have to construct mesh bodies already here, because boundary
  // elements need to know the parent body for their parent boundary when
  // they are installed into boundary-element-table
  //
  // NOTE: Do not convert here the bulk element ids into internal ids!
  // (external bulk element ids are stored in the input in the boundary elements!)
  //
  model->createMeshBodyTables();
  model->createMeshBodies();


  // Boundary element table
  // ----------------------
  model->allocateMeshBoundaryElements(actual_bndr_elem_count);
  bndrElementsAllocated = true;
 
  //---Read all boundary elements from the mesh input file
  if ( ECIF_OK != readMeshBoundaryElements() ) {
    return resetMeshData();
  }

  //---Close DB connection
  eio_close_mesh(info);
  eio_close(info);

  //---Update Gui
  gui->updateModelStatistics(model);
  gui->setWasUpdated("Mesh From Db");
  gui->setModelHasElmerMesh();

  return true;
}


bool
InputElmer::readMeshHeader()
{
  return true;
}


// Read boundary elements from the Elmer DB Mesh File and store them into the model.
Rc
InputElmer::readMeshBoundaryElements()
{
  UserInterface* gui = theControlCenter->getGui();
  strstream strm1, strm2;

  const ModelStatistics* ms = model->getModelStatistics();
  bool model_has_geometry = (ms->nofElements > 0);

  Rc rc;
  int info;
  int boundary_tag, last_boundary_tag;
  int element_id;
  int elem_type;
  meshElementCode elem_code;
  int ext_node_ids[MAX_NOF_NODES];
  int ext_parent1_id, ext_parent2_id;
  int int_parent1_id, int_parent2_id;
  double coord[3 * MAX_NOF_NODES];

  last_boundary_tag = NO_INDEX;

  while (1) {

    eio_get_mesh_bndry_element(element_id, boundary_tag,
                               ext_parent1_id, ext_parent2_id,
				                       elem_type, ext_node_ids,
                               coord, info);
    if (info == -1)
      break;

    bool add_as_boundary = true;
    bool add_as_sub_element = false;
    bool add_as_edge_element = false;
    bool add_as_vrtx_element = false;

    // These are not "real" boundary elements but
    // lower level in hierachy like vertex in 2D and edge in 3D
    if (modelDimension == ECIF_2D) {
      switch (elem_type) {
      case 101:
        add_as_boundary = false;
        add_as_sub_element = true;
        add_as_vrtx_element = true;
      default:
        break;
      }

    } else if (modelDimension == ECIF_3D) {
      switch (elem_type) {
      case 101:
        add_as_boundary = false;
        add_as_sub_element = true;
        add_as_vrtx_element = true;
        break;
      case 202:
      case 203:
        add_as_boundary = false;
        add_as_sub_element = true;
        add_as_edge_element = true;
       break;
      default:
        break;
      }
    }

    // Possibly the second parent is the correct one
    if ( ext_parent1_id <= 0 ) {
      ext_parent1_id = ext_parent2_id;
      ext_parent2_id = NO_INDEX;
    }

    // We prefer NO-INDEX as the no-parent indicator
    if (ext_parent1_id <= 0) {
      ext_parent1_id = NO_INDEX;
    }

    if (ext_parent2_id <= 0) {
      ext_parent2_id = NO_INDEX;
    }

    elem_code = model->convertElementType(elem_type);
    bool is_added;

    if ( add_as_boundary ) {
        

      rc = model->addMeshBoundaryElement(boundary_tag,
                                         elem_code,
                                         ext_parent1_id, ext_parent2_id,
                                         ext_node_ids,
                                         is_added);

      if ( rc != ECIF_OK ) {
        strm1 << "***ERROR when reading mesh boundary elements!" << ends;
        strm2 << "Could not insert boundary element " << element_id << ends;
        gui->showMsg(strm1.str());
        gui->showMsg(strm2.str());
        return rc;
      }
    } // is bndr_element

    // We must add the mesh element pending, defined just by node ids, and
    // then add the nodes later to proper mesh element (edge or vertex)
    // NOTE: This is used only if we have model geometry available, so that
    // we can pick the correct body element
    // We cannot *create* the bodyelement here, because we do not know whose
    // subelement it is!
    //
    if ( add_as_sub_element ) {
      
      if ( model_has_geometry ) {

        static BodyElement* se = NULL;
        if ( last_boundary_tag != boundary_tag ) {
          se = model->getBodyElementByBoundaryTag(boundary_tag);
        }
        if ( se != NULL ) {
          se->addPendingMeshElementAsNodes(elem_type, (int*)ext_node_ids);
        }
 
      } else {
        rc = model->addMeshInputElement(elem_type, element_id,
                                        NO_INDEX, boundary_tag,
                                        ext_node_ids);

        if ( rc != ECIF_OK ) {
          strm1 << "***ERROR when reading mesh boundary elements!" << ends;
          strm2 << "Could not add mesh edge/vetrex element " << element_id << ends;
          gui->showMsg(strm1.str());
          gui->showMsg(strm2.str());
          return rc;
        }

        if ( add_as_edge_element ) {
          nofInputEdgeElements++;
        } else {
          nofInputVertexElements++;
        }
      }
    }

    // Update last boundary tag
    if ( last_boundary_tag != boundary_tag ) {
      last_boundary_tag = boundary_tag;
    }

  }

  return ECIF_OK;
}


// Read an element from the Elmer Mesh File and store it into the model.
Rc
InputElmer::readMeshBulkElements()
{
  UserInterface* gui = theControlCenter->getGui();
  strstream strm1, strm2;

  Rc rc;
  int info;
  int body_tag;
  int ext_element_id;
  int ext_node_ids[MAX_NOF_NODES];
  int pdofs[10];
  int elem_type;
 
  const ModelStatistics* ms = model->getModelStatistics();
  bool model_has_geometry = (ms->nofElements > 0);

  while (1) {

    eio_get_mesh_element_conns(ext_element_id, body_tag, elem_type, pdofs, ext_node_ids, info);

    if (info == -1) break;

    meshElementCode elem_code = model->convertElementType(elem_type);
    
    if ( (modelDimension == ECIF_3D && elem_type < 500) ||
         (modelDimension == ECIF_2D && elem_type < 300)
       ) {
      
      // This will become a boundary element and will cause a body to be created!
      rc = model->addMeshInputElement(elem_type, ext_element_id,
                                      NO_INDEX, body_tag,
                                      ext_node_ids);
      nofInputBoundaryElements++;
  
    } else {
      rc = model->addMeshBulkElement(ext_element_id, body_tag,
                                     elem_code, ext_node_ids);
    }

    if (rc != ECIF_OK) {
      strm1 << "***ERROR when reading bulk elements!" << ends;
      strm2 << "Could not insert element " << ext_element_id << ends;
      gui->showMsg(strm1.str());
      gui->showMsg(strm2.str());
      return rc;
    }

  }

  return ECIF_OK;
}


// Read all nodes from the Elmer Mesh File and store them into the model
Rc
InputElmer::readMeshNodes(int nof_nodes)
{
  Rc rc;

  int* node_ids = new int[nof_nodes];
  double* coordinates = new double[3 * nof_nodes];

  int info;
  eio_get_mesh_nodes(node_ids, coordinates, info);

  Point3 p;

  int point_position = 0;

  for (int i = 0; i < nof_nodes; i++) {

    double* c = coordinates + point_position;

    p[0] = c[0];
    p[1] = c[1];
    p[2] = c[2];

    rc = model->addMeshNode( node_ids[i] - 1, node_ids[i], p);

    meshBoundingBox.extendByPoint(p);

    if (rc != ECIF_OK)
      return rc;

   point_position += 3;
  }

  delete[] node_ids;
  delete[] coordinates;

  return ECIF_OK;
}


// Clean table after an error
bool
InputElmer::resetMeshData()
{
  int info;

  //---Close EIO
  eio_close_mesh(info);
  eio_close(info);

  return false;
}


