/*****************************************************************************
 *
 *  Elmer, A Finite Element Software for Multiphysical Problems
 *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
 * 
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program (in file fem/GPL-2); if not, write to the 
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 *  Boston, MA 02110-1301, USA.
 *
 *****************************************************************************/

/***********************************************************************
Program:    ELMER Front
Module:     ecif_modelMeshManager.h
Language:   C++
Date:       13.04.00
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Model mesh manager class. Helper class for hanfling mesh related stuff

************************************************************************/

#ifndef _ECIF_MODEL_MESH_MGR_
#define _ECIF_MODEL_MESH_MGR_

#include "ecif_model.h"


class ModelMeshManager {
friend class Control;
friend class Model;

public:
  static void initClass(Model* model);
protected:
  ModelMeshManager();
  ~ModelMeshManager();

  int addMeshBoundaryElement(int bndr_tag, meshElementCode elem_code,
                             int paren1_tag, int parent2_tag, const int* node_ids,
                             bool& is_added_to_bndr);
  int addMeshBoundaryElement(BodyElement* boundary, meshElementCode elem_code,
                             int ext_parent1_tag, int ext_parent2_tag, const int* node_ids,
                             bool& is_added_to_bndr);
  int addMeshBoundaryElement(BodyElement* boundary, meshElementCode elem_code,
                             int int_parent1_tag, int int_parent2_tag,
                             int nof_nodes, const int* int_node_ids,
                             bool& is_added_to_bndr);
  Rc addMeshBulkElement(int ext_id, int material_id,
                        meshElementCode elem_code, const int* node_ids,
                        int nof_ngbrs = 0, const int* ngbr_ids = NULL);
  Rc addMeshNode(int int_id, int ext_id, Point3& point);
  Rc addMeshInputElement(int elem_type,
                         int ext_elem_id, int parent_tag, int ext_parent_tag,
                         int* node_ids);
  void allocateMeshBodies(int nof_bodies);
  void allocateMeshBoundaryElements(int nof_elements);
  void allocateMeshBulkElements(int nof_elements, int max_element_id);
  void allocateMeshNodes(int nof_nodes, int max_node_id);
  void allocateMeshInputElements(int nof_elements, int max_ext_elem_id);
  void calcMeshBoundaryNodeNormals();
  bool checkMeshElement(meshElementCode elem_code, int* node_ids, bool swap_to_ccw);
  bool checkMeshElements(MeshElementTable* table, bool swap_to_ccw);
  bool checkMeshElement303(int* node_ids, bool swap_to_ccw);
  bool checkMeshElement306(int* node_ids, bool swap_to_ccw);
  bool checkMeshElement404(int* node_ids, bool swap_to_ccw);
  bool checkMeshElement408(int* node_ids, bool swap_to_ccw);
  bool checkMeshElement504(int* node_ids, bool swap_to_ccw);
  bool checkMeshElement508(int* node_ids, bool swap_to_ccw);
  bool checkMeshElement510(int* node_ids, bool swap_to_ccw);
  bool checkMeshElement808(int* node_ids, bool swap_to_ccw);
  bool checkMeshElement820(int* node_ids, bool swap_to_ccw);
  void checkMeshElementType(int elem_type, bool& is_bulk, bool& is_bndr, bool& is_edge);
  void checkMeshInputBoundaryElements();
  void classifyMeshCornerElements();
  int convertElementCode(meshElementCode element_code);
  meshElementCode convertElementType(int element_type);
  void convertMeshBulkElementIdsExt2Int();
  void correctMeshZeroVelocityElements();
  void createMeshBodies();
  void createMeshBoundaries(int nof_bulk_bndr_elems, bool* free_bulk_bndr_elems);
  void createMeshBoundaryElements(int nof_bulk_bndr_elems,
                                  bool* free_bulk_bndr_elems,
                                  BodyElement*** boundary_table);
  void createMeshBodyTables();
  //void createMeshBulkElementEdges();
  void createMeshElementEdges(MeshElementTable* source, MeshElementTable* target);
  void createMeshSubElements(int hlevel, int min_nof_parents);
  void findMeshBoundaryBorders();
  void findMeshBoundaryParents(int nof_bulk_bndr_elems, bool* free_bulk_bndr_flags);
  void findMeshCornerElements();
  void findMeshElementEdges(MeshElementTable* source, int& nof_edge_elements);
  void findMeshElementNeighbors(MeshElementTable* source);
  void findMeshElementNodeParents(MeshElementTable* source, int nof_nodes, int*& nofNodeParents, int**& nodeParentIds);
  void findNofBulkBoundaryElements(int& nof_boundary_elements);
  int findSelectedMeshBoundaryElement(Renderer* renderer, Point3& ray_start, Point3& ray_dir,
                                      bool try_current_bndr,
                                      int& bndr_id,
                                      int& bd1_id, int& layer1_id,
                                      int& bd2_id, int& layer2_id);
  MeshCornerElement* getMeshCornerElement(int index);
  void getMeshBoundingBox(RangeVector rv) const;
  int getMeshBulkElementIdExt2Int(int ext_id);
  int getMeshBulkElementIdInt2Ext(int int_id);
  int getMeshInputElementIdExt2Int(int ext_id);
  int getMeshInputElementType(int elem_id);
  int getMeshNodeIdExt2Int(int ext_id);
  int getMeshNodeIdInt2Ext(int int_id);
  int getNofActiveMeshingBoundaryPoints();
  Rc installMeshInputBulkElements(bool clear_nodes);
  Rc installMeshInputBoundaryElements(bool clear_nodes);
  Rc installMeshInputElements(bool clear_nodes);
  bool isMeshInputBulkElement(int elem_id);
  //bool meshBoundaryElementSelected(Renderer* renderer, int fem_id);
  //bool meshBoundaryElementSelected(Renderer* renderer, Point3& ray_start, Point3& ray_dir);
  //bool meshBulkElementSelected(Renderer* renderer, int fem_id);

  void reallocateMeshBoundaryElements(int new_size);

  void removeMeshGeometry();
  void removeMeshInputElements();

  void resetMeshEdgesSelected();
  void resetMeshRendered();
  void resetMeshSelected();

  void selectMeshBoundaryElement(int fem_id);
  void selectMeshBoundaryElements();
  void selectMeshBoundaryElementsUndo();
  void selectMeshBoundaryElementsRedo();
  void selectMeshBoundaryElementsAll();
  void selectMeshBoundaryElementsByNeighbor();
  void selectMeshBoundaryElementsByNormal();
  void selectMeshBoundaryElementsByPlane();
  void setMeshBodyExt2IntFlag(int external_id);
  void setMeshBulkElementParentId(int elem_id, int parent_id);
  void setMeshInputElementIsAdded(int elem_id, bool value);
  void setMeshInputElementParentTag(int elem_id, int parent_tag);
  void setMeshInputElementExtParentTag(int elem_id, int ext_parent_tag);
  void setMeshNodes();
  void setMeshData(MeshData* mesh_data, MeshInfo* mesh_info, BoundBox* mesh_box);
  void setModelData(ModelData* model_data, ModelInfo* model_info);
  void unselectMeshBoundaryElement(int fem_id);

protected:
  void getMeshRangeVector(RangeVector rv) const;
  void normalizeMeshPoints();
  void sortMeshBoundaryIndices();
  void splitMeshCornerElement(MeshCornerElement* mce);

  static Control* theControlCenter;
  static Model* model;
  MeshInputElement* meshInputElements;
  int* meshInputElementsExt2Int;
  int meshInputElementsMaxExtId;

  int nofMeshInputElements;
  int nofMeshInputBoundaryElements;
  int nofMeshInputBulkElements;
  int nofMeshInputEdgeElements;
  int nofMeshInputVertexElements;
  int nofAllocMeshInputElements;

  // NOTE: These pointers are set by the Model
  // Do NOT delete them in this class!!!
  //
  MeshData* meshData;
  MeshInfo* meshInfo;
  BoundBox* meshBox;

  ModelData* modelData;
  ModelInfo* modelInfo;
};


#endif
