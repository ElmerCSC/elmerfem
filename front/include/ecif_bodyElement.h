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
Module:     ecif_bodyelement.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   A Base class for bodyelemets.
            A bodyelement is a separable entity in a body.
            They can be edges, surfaces etc. depending on the dimension.

************************************************************************/

#ifndef _ECIF_BODYELEMENT_
#define _ECIF_BODYELEMENT_

#include "ecif_def.h"
#include "ecif_def_stl.h"
#include "ecif_boundbox.h"
#include "ecif_geometry.h"
#include "ecif_modelObject.h"


// Used when linearizating boundaries
struct BoundaryPoint {
  
  BoundaryPoint();
  ~BoundaryPoint();
  bool isActiveInMif() { return activeInMeshing;}
  bool isVertex() { return vertexTag != NO_INDEX;}
  void copy(BoundaryPoint& other);

  bool activeInMeshing;    // A flag to control output for mif-files
  bool checked;            // A helper flag to handle fex. closed boundary copied points
  int tag;                 // Identifier (used also as node = vertex tag for mesh)
  GcPoint* point;          // Coordinate point, 3D
  int vertexTag;           // Vertex identifier if point is also a vertex
  char meshDensityType;    // H or R (points cannot have N!)
  double meshDensityValue; // Mesh density value
};


struct BodyElementGridHData {
  BodyElementGridHData();
  BodyElementGridHData(int nof_ids, int* grid_h_ids, int* mesh_indices);
  ~BodyElementGridHData();

  void setData(int nof_ids, int* grid_h_ids, int* mesh_indices, bool force = true);
  void setData(int mesh_index, char value_type, double value, bool force = false);
  int nofIds;
  int* gridHIds;
  int* meshIndices;

  // Possible directly inserted values (mainly for vertices)
  int meshIndex;
  int hValueType;
  double hValue;
};


struct BodyElementQuadGridData {
  BodyElementQuadGridData();
  BodyElementQuadGridData(int nof_ids, int* mesh_indices, int* n_values);
  ~BodyElementQuadGridData();

  void setData(int nof_ids, int* mesh_indices, int* n_values);
  int nofIds;
  int* meshIndices;
  int* nValues;
};


struct BodyElementModelData {
  BodyElementModelData();
  ~BodyElementModelData();
  void deleteBoundaryPoints();
  bool getSubElementTags(Model* model, objectType sub_type, int& count, int*& element_tags);
  bool getVertexTags(Model* model, int& count, int*& vertex_tags);
  int code;
  bool isClosedU;
  bool isClosedV;
  bool isIntraLayerBoundary;
  int nofBoundaryPoints;
  int nofBoundaryPointsU;
  int nofBoundaryPointsV;
  BoundaryPoint** boundaryPoints;
  int nofSubElements;
  int nofVertices;
  int parent1Id;
  int parent2Id;
  int parent1Tag;
  int parent2Tag;
  int parent1Layer;
  int parent2Layer;

  beStatus status;
  IdArray* subElementIds;
  IdArray* vertexIds;
};


struct BodyElementMeshDataBackup {
  BodyElementMeshDataBackup();
  ~BodyElementMeshDataBackup();

  int nofMeshElements;
  int* meshElementIds;
  short* meshElementDirs;
  int nofMeshBorderElements;
  int* meshBorderElementIds;
  short* meshBorderElementDirs;
};


struct BodyElementMeshData {
  BodyElementMeshData();
  ~BodyElementMeshData();

  int maxNofMeshElements;
  int nofMeshElements;
  int* meshElementIds;
  short* meshElementDirs;

  int maxNofPendingMeshElements;
  int nofPendingMeshElements;
  int** pendingMeshElementInfos;

  int nofMeshSelectedElements;
  int nofMeshBorderElements;
  int* meshBorderElementIds;
  short* meshBorderElementDirs;
  BodyElementMeshDataBackup* backupData;
};


struct LabelData {
  char* label;
  Point3 position;
};



class BodyElement : public ModelObject
{
friend class Control;
friend class Model;

public:
  BodyElement();
  BodyElement(int tag);
  BodyElement(int code, char* name);
  BodyElement(ecif_Element_X& trx_element);
  BodyElement(ecif_Vertex_X& trx_vertex);
  virtual ~BodyElement();

  virtual int addAllPendingMeshElements();
  virtual int addAllPendingSubElements() { return 0; }
  virtual void addCoveringElement(BodyElement* se, beStatus se_stat) {}
  virtual void addMeshBorderAsSubElement() {}
  virtual int addMeshElement(int elem_id, short direction);
  virtual int addPendingMeshElementAsNodes(int element_type, int* ext_node_ids);
  virtual int addPendingEdge(int edge_tag) { return 0; }
  virtual int addPendingVertex(int vertex_tag) { return 0; }
  beStatus addStatus(beStatus new_status);
  int addSubElement(BodyElement* be);
  int addSubElement(int sub_id);
  virtual void allocateMeshElements(int nof_mesh_elements);
  virtual void allocatePendingMeshElementInfos(int nof_infos);
  void backupMeshData();
  virtual void calcBoundaryPoints(double delta_u = 0.0, double delta_v = 0.0);
  virtual void checkBoundaryDiscretization(int mesh_index) {}
  virtual bool checkElementGroupData();
  void checkLastBoundaryTag(bool being_deleted = false);
  virtual void checkLastTag(bool being_deleted = false);
  virtual bool checkOuterBoundaries() = 0;
  virtual BodyElement* createElement(int nof_vertices, int* vertex_ids, ecif_geometryType gt) { return NULL;}
  virtual void createExtraVertices();
  virtual int compareOrientation(BodyElement* oe) = 0;
  virtual bool convertTags2Ids();
  void deleteCadData();
  void deleteMeshElements();
  void deleteMeshElements(bool* delete_flags);
  void deleteMeshDataBackup();
  virtual void draw(Renderer* renderer, flagName geometry_type, int body_id, int direction, bool is_first_loop = true);
  virtual void draw(int gmtr_index, Renderer* renderer, flagName geometry_type, int body_id, int direction, bool is_first_loop = true);
  void drawMesh(Renderer* renderer, int body_id, int direction);
  void drawMeshBorder(Renderer* renderer, bool selected);
  void drawLabel(Renderer* renderer, bool draw_sub_labels = true);
  virtual BodyElement* duplicate(bool* object_duplicated_flags) { return NULL; }
  matchType findCommonBoundary(BodyElement* other_element, BodyElement*& common );
  int findMeshBorder();
  virtual int findMeshBorderNodes(int buf_size, int* ids_buffer) = 0;
  int findSelectedMeshElement( Renderer* renderer,
                              Point3& ray_start, Point3& ray_dir, Point3& isec_point,
                              int& bndr_id,
                              int& body1_id, int& layer1_id,
                              int& body2_id, int& layer2_id);
  BoundBox* getBoundBox() {return ptrGmtr->getBoundBox();}
  int getBoundaryConditionId() const;
  int getBoundaryParameterId() const;
  int getBoundaryTag() const { return boundaryTag;}
  void getBoundaryPoints(int& count, BoundaryPoint**& points);
  double getDeltaU() { if ( ptrGmtr != NULL) return ptrGmtr->getDeltaU(); else return -1.0; }
  double getDeltaV() { if ( ptrGmtr != NULL) return ptrGmtr->getDeltaV(); else return -1.0; }
  void getDiscretizationData(int& nof_components, linDeltaType*& types, double*& valuesU, double*& valuesV, bool*& useFixedN);
  enum objectDrawingMode getDrawMode() { return drawMode; }
  enum objectDrawingState getDrawState() { return drawState; }
  const BodyElementGroup* getElementGroup() const;
  int getElementGroupId() const { return elementGroupId;}
  int getElementGroupTag() const { return elementGroupTag;}
  enum elementGroupType getElementGroupType() const;
  BodyElement* getFirstSubElement();
  Geometry* getGeometry() { return ptrGmtr;}
  const int* getGridHIds();
  const int* getGridMeshIndices();
  BodyElement* getLastSubElement();
  short getMeshElementDir(int index);
  virtual int getMeshElementId(int index);
  bool getMeshDensityValues(int mesh_index, char& type, int& nof_values, double values[4]);
  virtual int getMifGeometryTag(int index) const { return NO_INDEX; } 
  virtual const char* getName();
  int getNofBoundaryPoints() const;
  int getNofComponents(bool only_emf_components = false) const;
  int getNofMeshElements() const;
  virtual int getNofMifGeometries() const { return 0; };
  int getNofGridHIds() const;
  int getNofMeshSelectedElements() const;
  int getNofSubElements() const;
  virtual int getNofVertices() const;
  BodyElementList* getOuterBoundary();
  virtual double getParamArea(Geometry* gp) { return 0;}
  virtual ParamValues* getParamValues(Geometry* gp) {return NULL;}
  int getParentId(short parent) const;
  int getParentLayer(short parent) const;
  int getParentTag(short parent) const;
  bool getRangeVector(RangeVector rv);
  void getSelectedMeshElementIds(int& nof_ids, int*& ids);
  beStatus getStatus();
  int getSubElementTag(int index);
  virtual BodyElement* getSubElement(int index) { return NULL;}
  bool hasInside(BodyElement* other_element);
  bool hasZeroVelocityBC();
  int incrMeshElementCount();
  static void initClass(Model* model);
  virtual bool isBemBoundary();
  virtual bool isClosedU() { if ( ptrGmtr != NULL) return ptrGmtr->isClosedU(); else return false; }
  virtual bool isClosedV() { if ( ptrGmtr != NULL) return ptrGmtr->isClosedV(); else return false; }
  virtual bool isDiscretized() { if ( ptrGmtr != NULL) return ptrGmtr->isDiscretized(); else return false; }
  virtual bool isInnerBoundary();
  virtual bool isIntraLayerBoundary();
  virtual RayHit* isectionsWithXdir(GcPoint* startp, bool& negative_on_left);
  virtual bool isOk() {return 0;}
  virtual bool isOnSameAxis(GcPoint& p1, GcPoint& p2) {return false;}
  virtual bool isOnSamePlane(GcPoint& p1, GcPoint& p2, GcPoint& p3) {return false;}
  virtual void markActiveMeshObjects(bool*& active_object_flags);
  virtual void markActiveObjects();
  virtual matchType matchToLinear(BodyElement* be2, BodyElement*& common) { return MATCH_NONE;};
  virtual matchType matchToNurbs(BodyElement* be2, BodyElement*& common) { return MATCH_NONE;};
  virtual ostream& output_emf(ostream& out, short indent_size, short indent_level, bool output_geometry = true);
  virtual ostream& output_mif(ostream& out) { return  out;}
  virtual ostream& output_sif(ostream& out, short indent_size, short indent_level);
  virtual ostream& outputCoveringElementTags(ostream& out, int direction, int& nof_output_tag);
  GcPoint* param2Point(double u_p, double v_p = 0);
  ParamPair* point2Param(GcPoint* p);
  beStatus removeStatus(beStatus old_status);
  void removeMeshElements(int nof_elements, int* element_indices);
  void removeSelectedMeshElements(int* removed_indices_storage = NULL);
  void restoreMeshDataBackup();
  void restoreName();
  void setBoundaryParameterId(int pid);
  void setBoundaryPointsMeshDensityValue(int mesh_index);
  void setBoundaryTag(int tag);
  void setDiscretizationData(int nof_components, linDeltaType* types, double* valuesU, double* valuesV, bool* useFixedN);
  void setDrawMode(enum objectDrawingMode mode) { drawMode = mode; }
  void setDrawState(enum objectDrawingState state) { drawState = state; }
  void setElementGroupId(int id);
  void setElementGroupTag(int tag);
  void setGridHData(int nof_ids, int* grid_h_ids, int* mesh_indices, bool force = true);
  void setGridHValues(int mesh_index, char value_type, double value, bool force = false);
  void setIsIntraLayerBoundary(bool value);
  int setMeshElementCount(int count);
  void setMeshElementDir(int index, short direction);
  virtual void setMifTag(int& next_tag) {}
  void setName(char* new_name);
  void setParentId(short parent_nbr, int parent_id);
  void setParentIds(int parent1_id, int parent2_id);
  void setParentLayer(short parent_nbr, int parent_layer);
  void setParentLayers(int parent1_layer, int parent2_layer);
  void setParentTag(short parent_nbr, int parent_tag);
  void setParentTags(int parent1_tag, int parent2_tag);
  beStatus setStatus(beStatus new_status);
  void storeName();
  void update(ecif_Element_X& trx_element);
  bool updateGeometry(); // Update Matc-related geometry
  //beStatus updateStatus(beStatus st) {return status |= st;}

protected:
  virtual int calcDirection(BodyElement* sub_elm) { return 0;}
  virtual void drawSubElements(Renderer* renderer, flagName geometry_type, int body_id, int direction, bool is_first_loop);
  static int getLastBoundaryTag() { return last_boundaryTag; }
  virtual int getLastTag() { return NO_INDEX; }
  virtual BodyElementList* findOuterBoundary() { return NULL;}
  virtual void formLinearizedGeometry();
  virtual void init();
  virtual void initLabelData();
  void initName(char* be_name = NULL);
  virtual int newTag() { return NO_INDEX;}
  virtual ostream& outputBoundaryPointTags(ostream& out, int indent_size);
  void reallocateMeshElements();
  void reallocatePendingMeshElementInfos();
  static void setLastBoundaryTag(int lbtag) { last_boundaryTag = lbtag; }
  virtual void setLastTag(int ltag) {}
  virtual void update();

  int elementGroupId;
  int elementGroupTag;
  int boundaryTag;
  static int last_boundaryTag;
  enum objectDrawingMode drawMode;
  enum objectDrawingState drawState;

  int txBndrConditionId;

  int boundaryParameterId;

  char* name_old;
  Geometry* ptrGmtr;
  BodyElementModelData* modelData;
  BodyElementGridHData* gridHData;
  BodyElementMeshData* meshData;
  LabelData* labelData;
};


#endif

