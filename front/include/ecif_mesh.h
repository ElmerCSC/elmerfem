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
Module:     ecif_mesh.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Classes for mesh data.

************************************************************************/

#ifndef _ECIF_MESH_
#define _ECIF_MESH_

#include "ecif_def.h"


// Table of mesh elements (abstract base class)
class MeshElementTable {
  friend class Model;
  friend class ModelMeshManager;
  friend class ModelOutputManager;
  friend class Body;
  friend class BodyElement;
public:
  MeshElementTable(int size);
  virtual ~MeshElementTable();
  virtual void calcCenter(int elem_index);
  virtual int calcLineIntersections(int elem_index, short direction,
                                    Point3& lstart, Point3& ldir,
                                    Point3* isec_points) { return 0;}
  virtual void calcNormalDistance(int elem_index);
  void calcNormalsAndDistances();
  virtual void calcNormalVector(int elem_index) {}
  virtual void centerPoint(int elem_index, Point3& center);
  virtual void createTableEntry(int int_id, meshElementCode elem_code,
                                const int* parent_pair_ids,
                                int nof_nodes, const int* entry_node_ids);
  virtual void draw(Renderer* renderer, bool filled, bool selected) {};
  virtual void drawEdges(Renderer* renderer, MeshElementTable* edge_table, int index, bool selected);
  virtual void drawEdges(Renderer* renderer, const Point3* node_data, int index, bool selected) {}
  virtual void drawElement(Renderer* renderer, int index, bool selected) {};
  virtual void drawElement(Renderer* renderer, int index, short direction, bool selected) {};
  int findSubElementIndex(int elem_index, int& direction, int& start_position,
                          int nof_sub_nodes, const int* sub_node_ids);
  bool getChecked(int index) {return checked[index];}
  virtual const int* getEdgeIds(int elem_index);
  virtual meshElementCode getElementCode(int elem_index);
  virtual int* getElementNodesBuffer() {return MeshElementTable::elementNodesBuffer;}
  virtual const int* getNodeIds(int elem_index);
  virtual const int* getNodeIdsReversed(int elem_index);
  int getParentId(int elem_index, short parent_index);
  bool getSelected(int index) {return selected[index];}
  void getSubElementNodeIndices(int elem_code, short sub_elem_index, short sub_start_pos,
                                int* node_index_buffer);
  static void initClass(Model* mdl);
  virtual bool isBoundaryFace(int index, int face_index, bool& is_first_parent) { return false;}
  bool isSameElement(int my_index, short& direction, short& start_position,
                     int other_elem_code, int other_nof_nodes, int* other_node_ids);
  int NofElements() { return nofElements;}
  virtual void resetActionLevels();
  virtual void resetChecked();
  virtual void resetChecked(int index) { checked[index] = false; }
  virtual void resetRendered();
  virtual void resetRendered(int index) { rendered[index] = false; }
  virtual void resetSelected();
  virtual void resetSelected(int index) { selected[index] = false; }
  void reverseElementNodes(int index);
  bool resize(int new_size, bool* copy_flags);
  void setDirInParents(int index, short dir_in_parent1, short dir_in_parent2);
  virtual void setSelected(int index, bool value, bool toggle);
  virtual bool setEdgeIds(int index, int nof_ids, int* edge_ids);
  virtual bool setEntry(int index, meshElementCode elem_code, int nof_ids, int* node_ids);
  static void setMeshData(MeshData* mesh_data, MeshInfo* mesh_info);
  static void setMeshNodes(Point3* mesh_nodes) { meshNodes = mesh_nodes;}
  void setParentIds(int index, int* parent_ids);
  short updateActionLevel(short change_in_level);
  int updateParentId(int index, int old_id, int new_id);

protected:
  void resize_component(char*& component,int old_size, int new_size, char init_value, bool* copy_flags);
  void resize_component(bool*& component,int old_size, int new_size, bool init_value, bool* copy_flags);
  void resize_component(short*& component,int old_size, int new_size, short init_value, bool* copy_flags);
  void resize_component(int*& component,int old_size, int new_size, int init_value, bool* copy_flags);
  void resize_component(double*& component,int old_size, int new_size, double init_value, bool* copy_flags);

  void resize_component(char**& component,int old_size, int new_size, int nof_values, char init_value, bool* copy_flags);
  void resize_component(bool**& component,int old_size, int new_size, int nof_values, bool init_value, bool* copy_flags);
  void resize_component(int**& component,int old_size, int new_size, int nof_values, int init_value, bool* copy_flags);
  void resize_component(short**& component,int old_size, int new_size, int nof_values, short init_value, bool* copy_flags);
  void resize_component(double**& component,int old_size, int new_size, int nof_values, double init_value, bool* copy_flags);

  template <class T> void resize_component_impl1(T*& component,int old_size, int new_size,
                                                T init_value, bool* copy_flags);

  template <class T, class T1> void resize_component_impl2(T*& component,int old_size, int new_size,
                                                           int nof_values, T1 init_value, bool* copy_flags);


  template <class T> void resize_component_impl3(T**& component,int old_size, int new_size,
                                                           int nof_values, T init_value, bool* copy_flags);

  char* actionLevels;       //  flag to handle selection levels etc.
  static int elementNodesBuffer[27];
  Point3* centers;          // center point of the element
  bool* checked;            // helper flag for misc processing
  short currentActionLevel; // a flag value for currently active level
  short** dirInParents;     // Direction of the sub element in the parents
  int** edgeIds;            // Edge element ids
  meshElementCode* elementCodes; // Element code like: MEC_303 for linear triangle
  short maxActionLevel;     // a flag value for maximum action level applied in elements
  static Point3* meshNodes; // Class variable (for fast access to points)
  static const MeshData* meshData;// Class variable
  static const MeshInfo* meshInfo;// Class variable
  static Model* model;      // Class variable
  int** neighborIds;        // Neighbor element ids
  int** nodeIds;            // Element node ids
  int nofElements;          // Total number of enries in the table
  double* normalDistances;  // center distance from a plane through origin (orthogonal to the normal)
  Point3* normals;          // normal vectors (when meaningful!)
  Point3** nodeNormals;     // normal vectors at element nodes (when meaningful, normally used only for boundary nodes!)
  int** parentIds;          // Parent ids (NOTE: nof parents = 2, if parents defined!)
  bool* rendered;           // true if boundary element is already drawn (for bndr elemetns)
  double* rSquares;         // R2-values for the elements (when sensible)
  bool* selected;           // true if selected
  bool* splitted;           // element is splitted into sub-elements (like zero corners)
  void resetState(bool* target);

};


class MeshBulkElementTable : public MeshElementTable {
  friend class Body;
  friend class Model;
public:
  MeshBulkElementTable(int size);
  ~MeshBulkElementTable();
  void drawEdges(Renderer* renderer, const Point3* node_data, int index, bool selected);
  virtual meshElementCode getElementCode(int elem_index);
  int* getElementNodesBuffer() {return MeshBulkElementTable::elementNodesBuffer;}
  bool isBoundaryFace(int index, int face_index, bool& is_first_parent);
  void resetRendered();
protected:
  static int elementNodesBuffer[27];
};


class MeshFaceElementTable : public MeshElementTable {
  friend class Model;
public:
  MeshFaceElementTable(int size, bool store_neighbors);
  ~MeshFaceElementTable();
  void calcNormalVector(int elem_index);
  int calcLineIntersections(int elem_index, short direction,
                            Point3& lstart, Point3& ldir, Point3* isec_points);
  void drawElement(Renderer* renderer, int index, short direction, bool selected);
  int* getElementNodesBuffer() {return MeshFaceElementTable::elementNodesBuffer;}

protected:
  int calcTriangleLineIntersections(int elem_index, short direction,
                                    Point3& lstart, Point3& ldir, Point3* isec_points);
  int calcQuadriLineIntersections(int elem_index, short direction,
                                  Point3& lstart, Point3& ldir, Point3* isec_points);
  static int elementNodesBuffer[27];
};


class MeshEdgeElementTable : public MeshElementTable {
  friend class Model;
public:
  MeshEdgeElementTable(int size, bool store_neighbors);
  ~MeshEdgeElementTable();
  int calcLineIntersections(int elem_index, short elem_dir,
                             Point3& lstart, Point3& ldir, Point3* isec_points);
  void calcNormalVector(int elem_index);
  void draw(Renderer* renderer, bool selected);
  void drawElement(Renderer* renderer, int index, bool selected);
  void drawElement(Renderer* renderer, int index, short direction, bool selected);
  int* getElementNodesBuffer() {return MeshEdgeElementTable::elementNodesBuffer;}
  const int* getNodeIds(int elem_index, short direction = 1);
  static void setPickingTolerance( double value) {pickingTolerance = value;}
protected:
  static int elementNodesBuffer[27];
  static double pickingTolerance;
};


class MeshVertexElementTable : public MeshElementTable {
  friend class Model;
public:
  MeshVertexElementTable(int size);
  ~MeshVertexElementTable();
  void drawElement(Renderer* renderer, int index, bool selected);
  int* getElementNodesBuffer() {return MeshVertexElementTable::elementNodesBuffer;}
  const int* getNodeIds(int elem_index, short direction);
protected:
  static int elementNodesBuffer[27];
};

#endif