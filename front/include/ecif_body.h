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
Module:     ecif_body.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   A Base class for model's bodies.
  A body is a separable entity in the model geometry.
  Body consist of body-elements, which can be edges, surfaces
  etc. depending on the topology dimension.
  Often also called a part in CAD programs.

************************************************************************/

#ifndef _ECIF_BODY_
#define _ECIF_BODY_

#include "ecif_def.h"
#include "ecif_def_stl.h"
#include "ecif_def_trx.h"
#include "ecif_boundbox.h"
#include "ecif_modelObject.h"


class BodyPair : public ModelObject
{
public:
  BodyPair();
  BodyPair(const Body* bd1, const Body* bd2);
  BodyPair(int tag, const Body* bd1, const Body* bd2);
  ~BodyPair();
  static void initClass(Model* model);
  void getBodies(const Body*& bd1, const Body*& bd2) { bd1 = body1; bd2 = body2;}
  const char* getName();
  int setTag();
  int setId();

protected:
  void init();
  const Body* body1;
  const Body* body2;
  static int last_tag;
};


class Body : public ModelObject
{
friend class Control;
friend class Model;
public:
  Body();
  Body(bodyGmtrType type, int ext_tag, char* name, colorIndices color_index = DEFAULT_COLOR_INDEX);
  Body(bodyGmtrType type, int int_tag, int ext_tag, char* name, colorIndices color_index = DEFAULT_COLOR_INDEX);
  Body(ecif_Body_X& trx_body, bool add_default_layer = false);
  virtual ~Body();
  virtual bool addAllLoopElements(int layer);
  virtual int addAllPendingElements(int layer);
  virtual int addAllPendingVertices(int layer) = 0;
  bool addElementLoopTags(int layer, int nof_loops, int* loop_tags);
  bool addElementGroupTags(int nof_groups, int* group_tags);
  bool addElementGroupTag(int group_tag);
  int addElement(int layer, BodyElement* be);
  int addLayer(ecif_BodyLayer_X& trx_layer);
  int addMeshElements(int nof_new_elements, int* new_elment_ids);
  int addPendingElement(int layer,int element_tag);
  int addPendingVertex(int layer,int vertex_group_id, int vertex_tag);
  virtual bool acceptsStructuredMesh(int layer) { return false; }
  void calcMeshBox();
  bool canHaveZeroVelocityElements();
  bool checkElements();
  bool checkLayerIndex(int layer);
  ecif_modelStatus checkStatus();
  virtual bool convertTags2Ids();
  void deleteMeshElements();
  virtual void draw(Renderer* renderer, flagName geometry_type, objectDrawingMode dmode);
  BoundBox* getBoundBox() {return boundbox;}
  BoundBox* getBoundBox(int layer);
  void getColor(Color4& color);
  void getColor(Color4f& color);
  colorIndices getColorIndex() { return colorIndex;}
  int getBodyForceId();
  BodyLayer* getBodyLayer(int layer);
  int getBodyParameterId() { return bodyParameterId;}
  int getDirectedElementLoopId(int layer, int index);
  int getDirectedElementLoopTag(int layer, int index);
  enum objectDrawingMode getDrawMode() { return drawMode;}
  enum objectDrawingState getDrawState() { return drawState;}
  int getElementCount();
  int getElementCount(int layer);
  BodyElement** getElementArray();
  BodyElement** getElementArray(int layer);
  BodyElement* getElement(int index);
  BodyElement* getElement(int layer, int index);
  BodyElement* getElementById(int id);
  BodyElement* getElementById(int layer, int id);
  int getElementGroupId(int index);
  void getElementGroupIds(int& nof_groups, int*& group_ids);
  int getElementGroupTag(int index);
  void getElementGroupTags(int& nof_groups, int*& group_tags);
  BodyElementLoop* getElementLoop(int index);
  BodyElementLoop* getElementLoop(int layer, int index);
  void getElementIds(int& nof_elements, int*& element_ids);
  void getElementIds(int layer, int& nof_elements, int*& element_ids);
  void getElementTags(int& nof_elements, int*& element_tags);
  void getElementTags(int layer, int& nof_elements, int*& element_tags);
  int getElementLoopId(int index);
  int getElementLoopTag(int index);
  int getElementLoopId(int layer, int index);
  int getElementLoopTag(int layer, int index);
  int getEquationId();
  const int* getExcludedMeshIndices(int layer);
  int getExternalTag() {return externalTag;}
  const int* getGridParameterMeshIndices(int layer);
  const int* getGridParameterIds(int layer);
  int getInitialConditionId();
  BodyLayer* getLayer(int layer);
  int getLayerId(int layer);
  int getLayerIndexById(int layer_id);
  int getLayerIndexByTag(int layer_tag);
  int getLayerTag(int layer);
  int getLayerTagById(int layer_id);
  int getMaterialId();
  int getMeshElementId(int index);
  bool getMeshDensityValue(int layer, int mesh_index, char& type, double& value);
  int getMeshQuadGridN(int layer, int mesh_index, int element_id);
  const char* getName();
  const char* getName(int layer);
  int getNofBoundaries(int layer, elementType type);
  int getNofElementGroups();
  int getNofElementLoops(int layer);
  int getNofLayers() {return nofLayers;}
  int getNofMeshElements() {return nofMeshElements;}
  int getNofMifLayers();
  void getRangeVector(RangeVector rv) {boundbox->getRangeVector(rv);}
  bool getSeparateFlag() { return separateFlag; }
  int getGmtrType() { return gmtrType;}
  int getTplgType() { return tplgType;}
  int getType() { return type;}
  bool hasInside(int layer, Body* other_body, int other_body_layer);
  void incrBndrCount(int layer, elementType type);
  static void initClass(Model* model);
  void initName();
  virtual bool isClosed() { return tplgType == CLOSED_BODY; }
  bool isExcludedFromCurrentMesh();
  bool isExcludedFromCurrentMesh(int layer);
  bool isExcludedFromMesh(int mesh_index);
  bool isExcludedFromMesh(int layer, int mesh_index);
  virtual bool isCadBody() { return gmtrType == GEOM_CAD; }
  virtual bool isBemBody() { return type == BEM_BODY; }
  virtual bool isMeshBody() { return gmtrType == GEOM_MESH; }
  virtual bool isOk();
  virtual bool isOpen() { return tplgType == OPEN_BODY; }
  virtual bool isSeparable(int layer);
  virtual bool isVirtual() { return type == VIRTUAL_BODY; }
  virtual void markActiveObjects();
  virtual void markActiveMeshObjects(int mesh_index, bool*& active_object_flags);
  virtual ostream& output_emf(ostream& out, short indent_size, short indent_level);
  virtual ostream& output_sif(ostream& out, short indent_size, short indent_level);
  virtual ostream& output_mif(ostream& out);
  int removeMeshElements(bool* remove_flags);
  bool selectLayer(int layer);
  virtual void separate(bool* object_duplicate_flags);
  void setSeparateFlag(bool value) { separateFlag = value; }
  void setBodyForceId(int bf_id) { bodyForceId = bf_id; }
  void setBodyParameterId(int pid);
  void setColor(Color4& color);
  void setColorIndex(colorIndices color_index);
  void setDrawMode(enum objectDrawingMode mode) { drawMode = mode; }
  void setDrawState(enum objectDrawingState state) { drawState = state; }
  void setEquationId(int eq_id) { equationId = eq_id; }
  void setElementLoops(int layer, int nof_loops, int* loop_ids, int* loop_tags);
  void setExternalTag(int ext_tag) {externalTag = ext_tag;}
  void setInitialConditionId(int ic_id) { initialConditionId = ic_id; }
  void setLayerTag(int layer, int tag);
  void setMaterialId(int mt_id) { materialId = mt_id; }
  void setMeshElements(int nof_elements, int* element_ids);
  void setGmtrType(bodyGmtrType value) { gmtrType = value; };
  void setTplgType(bodyTplgType value) { tplgType = value; };
  void setType(bodyType value) { type = value; };
  virtual void swapBodyElements(int layer, BodyElement* b1, BodyElement* b2);
  int removeElement(int layer, BodyElement* be);
  virtual void updateBoundBox();
  void updateBoundBox(int layer, RangeVector rv);
  void updateMinimumBox(RangeVector rv);

protected:
  virtual bool check();
  virtual bool checkVirtual();
  virtual void drawMesh(Renderer* renderer);
  void init(bool add_default_layer = true);
  void removeElementLoops(int layer);
  void removeElementGroups();
  void removeLayers();
  void updateBodyElementMeshH(int layer, int updateMode);

  // Body level stuff
  //
  BoundBox* boundboxMesh;
  BoundBox* boundbox;
  bool checked;
  Color4 color;
  colorIndices colorIndex;
  enum objectDrawingMode drawMode;
  enum objectDrawingState drawState;
  bool separateFlag;
  int* meshElementIds;
  int externalTag;
  static int last_tag;
  int maxNofMeshElements;
  int nofMeshElements;

  bool status;
  enum bodyType type;
  enum bodyGmtrType gmtrType;
  enum bodyTplgType tplgType;

  int bodyParameterId;
  int initialConditionId;
  int bodyForceId;
  int equationId;
  int materialId;

  // Body Layer level stuff
  //
  int nofLayers;
  int* layerIds;
  BodyElementTable** belements;
  BoundBox** boundboxes;
  IdList** elementLoopIds;
  IdList** elementLoopTags;
  int* nofElements;
  int* nofElementLoops;
  IdArray** pendingElementTags;
  IdArray** pendingVertexGroups;
  IdArray** pendingVertexTags;
  int* nofInnerBoundaries;
  int* nofOuterBoundaries;

  // For virtual bodies
  int nofElementGroups;
  IdList* elementGroupIds;
  IdList* elementGroupTags;
};

#endif
