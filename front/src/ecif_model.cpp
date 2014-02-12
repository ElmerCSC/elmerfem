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
Module:     ecif_model.cpp
Language:   C++
Date:       12.11.98
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
#include "ecif_bodyElementGroup.h"
#include "ecif_bodyElementLoop.h"
#include "ecif_bodyForce.h"
#include "ecif_bodyParameter.h"
#include "ecif_boundbox.h"
#include "ecif_boundaryCondition.h"
#include "ecif_boundaryParameter.h"
#include "ecif_control.h"
#include "ecif_calculator.h"
#include "ecif_constant.h"
#include "ecif_coordinate.h"
#include "ecif_datafile.h"
#include "ecif_equation.h"
#include "ecif_equationVariables.h"
#include "ecif_geometry.h"
#include "ecif_gridH.h"
#include "ecif_gridParameter.h"
#include "ecif_initialCondition.h"
#include "ecif_inputIdeas.h"
#include "ecif_inputElmer.h"
#include "ecif_inputThetis.h"
#include "ecif_material.h"
#include "ecif_mesh.h"
#include "ecif_model.h"
#include "ecif_modelMeshManager.h"
// #include "ecif_modelOutputManager.h"
#include "ecif_modelParameter.h"
#include "ecif_model_aux.h"
#include "ecif_modelObject.h"
#include "ecif_parameter.h"
#include "ecif_parameterField.h"
#include "ecif_renderer.h"
#include "ecif_simulationParameter.h"
#include "ecif_solver.h"
#include "ecif_solverControl.h"
#include "ecif_timer.h"
#include "ecif_timestep.h"
#include "ecif_userinterface.h"
#include "ecif_userSettings.h"

void pressAnyKey();
extern char read_buffer[];

Control* Model::theControlCenter = NULL;
NameTable* Model::colorNameTable = new NameTable;
RGBColorTable* Model::colorValueTable = new RGBColorTable;
double Model::meshInputUnit = 1.0;

Model::Model(char* model_name, ecif_modelSource source, char* in_file_name)
{
  id = 1;
  lastObjectId = 0; // System-wide unique object id (for anyone who needs it)

  Body::initClass(this);
  BodyElement::initClass(this);
  BodyElementGroup::initClass(this);
  BodyElementLoop::initClass(this);
  BodyForce::initClass(this);
  BodyLayer::initClass(this);
  BodyPair::initClass(this);
  BoundaryCondition::initClass(this);
  Constant::initClass(this);
  Coordinate::initClass(this);
  Datafile::initClass(this);
  Equation::initClass(this);
  Geometry::initClass(this);
  GridParameter::initClass(this);
  InitCondition::initClass(this);
  Input::initClass(this);
  Material::initClass(this);
  ModelMeshManager::initClass(this);
  ModelOutputManager::initClass(this);
  MeshElementTable::initClass(this);
  ModelObject::initClass(this);
  Parameter::initClass(this);
  ParameterField::initClass(this);
  Renderer::initClass(this);
  Solver::initClass(this);
  Timestep::initClass(this);

  modelData = new ModelData;
  modelFlags = new bool[MAX_NOF_FLAG_NAMES];
  for (short i = 0; i < MAX_NOF_FLAG_NAMES; i++) {
    modelFlags[i] = false;
  }
  modelInfo = new ModelInfo(model_name, source, in_file_name);
  modelStatistics = new ModelStatistics;
  parallelInfo = new ParallelInfo;
  pickInfo = new PickInfo;

  modelBox = new BoundBox;

  initMeshData();

  initSplitCombineInfos();

  meshManager = new ModelMeshManager();
  meshManager->setMeshData(meshData, meshInfo, meshBox);
  meshManager->setModelData(modelData, modelInfo);

  outputManager = new ModelOutputManager();

  MeshElementTable::setMeshData(meshData, meshInfo);

}


Model::~Model()
{
  modelInfo->dimension = ECIF_ND;

  deleteMeshData();

  delete modelBox;
  delete modelData;
  delete[] modelFlags;

  delete outputManager;
  delete meshManager;
}


// Method adds a new body to the model.
bool
Model::addBody(Body* body)
{
  if ( !body->objectIsOk() ) return false;

  GBODY = body;
  modelStatistics->nofBodies++;
  return true;
}


// Method adds a new transfer-format body to the model.
bool
Model::addBody(ecif_Body_X& trx_body, bool add_default_layer)
{
  int body_tag = trx_body.tag;
  Body* body = this->getBodyByTag(body_tag);

  //-Possible old body with the same id is deleted first
  if (NULL != body) {
    removeBody(body);
    delete body;
  }

  //-New body is created with preset id!
  switch (modelInfo->dimension) {
  case 2:
    body = new Body2D(trx_body, add_default_layer);
    break;
  case 3:
    body = new Body3D(trx_body, add_default_layer);
    break;
  default:
    UserInterface* gui = theControlCenter->getGui();
    strstream strm;
    strm << "***ERROR Model dimension is zero, cannot add Body " << body_tag << ends;
    gui->showMsg(strm.str());
    return NO_INDEX;
  }

  return this->addBody(body);
}


// Method adds a new bodyelement into the model.
bool
Model::addBodyElement(BodyElement* be)
{
  if ( !be->objectIsOk() ) {
    return false;
  }

  be->checkLastTag();
  be->checkLastBoundaryTag();

  if ( OT_VERTEX == be->getObjectType() ) {
    addVertex(be);
  } else {
    modelStatistics->nofElements++;
  }

  return true;
}


// Method adds a new bodyelement into the model.
bool
Model::addBodyElement(BodyElement* be, bool add_to_bodies)
{
  if (add_to_bodies) {
    int body1_tag = be->getParentTag(1);
    int body1_layer = be->getParentLayer(1);
    int body2_tag = be->getParentTag(2);
    int body2_layer = be->getParentLayer(2);
    return addBodyElement(be, body1_tag, body1_layer, body2_tag, body2_layer);

  } else {
    return addBodyElement(be);
  }
}


bool
Model::addBodyElement(BodyElement* be,
                      int body1_tag, int body1_layer,
                      int body2_tag, int body2_layer)
{
  Body* body1 = getBodyByTag(body1_tag);
  Body* body2 = getBodyByTag(body2_tag);

  // No success!
  if ( body1 == NULL ) {
    return false;
  }

  // Add to relevant bodies and also create corresponding
  // inner/outer boundary entry

  // Outer boundary
  if (body2 == NULL) {
    body1->addElement(body1_layer, be);
    be->addStatus(BE_OUTER);

  // Inner boundary
  } else {
    body1->addElement(body1_layer, be);
    body2->addElement(body2_layer, be);
    be->addStatus(BE_INNER);
  }

  // Add the body element itself
  return addBodyElement(be);
}


// Method adds a new transfer-format bodyelement (edge or face) into the model.
bool
Model::addBodyElement(ecif_Element_X& tx)
{
  int be_tag = tx.tag;

  BodyElement* be = NULL;

  if ( tx.tplg_type == ECIF_VERTEX )
    be = getVertexByTag(be_tag);

  else if ( tx.tplg_type == ECIF_EDGE )
    be = getEdgeByTag(be_tag);

  else if ( tx.tplg_type == ECIF_FACE )
    be = getFaceByTag(be_tag);

  //-Possible old bodyelement with the same tag
  // is updated!
  if (NULL != be) {
    //delete be;
    be->update(tx);
    return true;
  }

  //-New bodyelement is created
  switch (tx.tplg_type) {

  case ECIF_VERTEX:
    be = new BodyElement1D(tx);
    break;
  case ECIF_EDGE:
    be = new BodyElement2D(tx);
    break;
  case ECIF_FACE:
    be = new BodyElement3D(tx);
    break;
  default:
    return NO_INDEX;
  }

  return this->addBodyElement(be);
}


// Method adds a new transfer-format bodyelement (vertex) to the model.
bool
Model::addBodyElement(ecif_Vertex_X& tx)
{
  BodyElement* be = new BodyElement1D(tx);

  return addBodyElement(be);
}


// Method adds a new bodyelement into the model.
bool
Model::addBodyElementGroup(BodyElementGroup* beg)
{
  if ( !beg->objectIsOk() ) return false;

  modelStatistics->nofElementGroups++;
  return true;
}


// Method adds a new bodyelement-loop into the model.
bool
Model::addBodyElementLoop(ecif_ElementLoop_X& tx)
{
  int bel_tag = tx.tag;

  BodyElementLoop* bel = this->getBodyElementLoopByTag(bel_tag);

  //-Possible old bodyelement-loop with the same id is deleted first
  if (NULL != bel) {
    removeBodyElementLoop(bel->Id());
    delete bel;
  }


  //-New bodyelement-loop is created with preset id!
  //
  // NOTE: Loop element type is given as a generic boundary!
  //
  bel = new BodyElementLoop(bel_tag,
                            tx.nof_elements,
                            tx.element_tags,
                            tx.is_open,
                            OT_BOUNDARY);

  return this->addBodyElementLoop(bel);
}


// Method adds a new bodyelement into the model.
bool
Model::addBodyElementLoop(BodyElementLoop* bel)
{
  if ( !bel->objectIsOk() ) return false;

  modelStatistics->nofElementLoops++;
  return true;
}


bool
Model::addBodyPair(BodyPair* bp)
{
  if ( !bp->objectIsOk() ) return false;

  modelStatistics->nofBodyPairs++;
  return true;

}


void
Model::addColorName(char* name, Color4 color)
{
  static char buf[1024];
  int i;

  Color4* c = new Color4[1];

  int len = strlen(name);

  for (i = 0; i < len; i++) {
    buf[i] = tolower(name[i]);
  }
  buf[len] = '\0';

  string nm(buf);

  for (i = 0; i < 4; i++) {
    (*c)[i] = color[i];
  }

  (*colorValueTable)[nm] = c;
}


void
Model::addColorValue(int color_value, char* name)
{
  (*colorNameTable)[color_value] = name;
}


void
Model::addCoveringElementList(int elem_id, IdList* se_list)
{
  IdList* old_list = (*modelData->coveringElementTable)[elem_id];

  if ( old_list != NULL ) {
    delete old_list;
  }

  // Temprary debugging
  IdList::iterator pos = se_list->begin();

  while (pos != se_list->end()) {
    int id = *pos++;
  }

  (*modelData->coveringElementTable)[elem_id] = se_list;
}


void
Model::addMatcDefinition(char* def)
{
  char* res = mtc_domath(def);

  if ( res != NULL && res[0] != '\0' ) {
    UserInterface* gui = theControlCenter->getGui();
    strstream strm;
    int len = strlen(def);

    for (int i = 0; i < len; i++ ) {
      if ( def[i] < 0 ) {
        def[i] = '\0';
        break;
      }
    }
    strm << "***" << res << " in: " << def << ends;
    gui->showMsg(strm.str());
  }
}


Rc
Model::addMeshBoundaryElement(int bndr_tag, meshElementCode elem_code,
                              int ext_parent1_id, int ext_parent2_id,
                              const int* ext_node_ids,
                              bool& is_added_to_bndr)
{
  int elem_id = meshManager->addMeshBoundaryElement(bndr_tag, elem_code, ext_parent1_id, ext_parent2_id,
                                                    ext_node_ids, is_added_to_bndr);

  if ( elem_id == NO_INDEX )
    return ECIF_ERROR;
  else
    return ECIF_OK;
}


Rc
Model::addMeshBulkElement(int ext_id, int material_id,
                          meshElementCode elem_code, const int* ext_node_ids,
                          int nof_ngbrs, const int* ext_ngbr_ids)
{
  return meshManager->addMeshBulkElement(ext_id, material_id, elem_code, ext_node_ids, nof_ngbrs, ext_ngbr_ids);
}


// Add a temporary element (bulk or boundary )
// This is used when only the total number elements is
// known and we do not want to read the mesh file twice
// to know separately the nofBulkElements and nofBoundaryElements
Rc
Model::addMeshInputElement(int elem_type,
                           int ext_id, int body_tag, int bndr_tag,
                           int* ext_node_ids)
{
  return meshManager->addMeshInputElement(elem_type, 
                                          ext_id, body_tag, bndr_tag,
                                          ext_node_ids);
}



Rc
Model::addMeshNode(int int_id, int ext_id, Point3& point)
{
  return meshManager->addMeshNode(int_id, ext_id, point);
}


// Add new model object. Increase id-cunter and set new counter as the object id
//
bool
Model::addModelObject(ModelObject* object, enum objectType otype)
{
  object->setObjectType(otype);

  object->setId(++lastObjectId);

  modelData->modelObjects->push_back(object);

  modelStatistics->nofModelObjects++;

  return true;
}


// Method adds a parameter of the type "parameter-type" to the model.
bool
Model::addParameter(ecif_parameterType parameter_type, Parameter* parameter)
{
  ParameterTable* table;

  // NOTE: This will be scalar pointer to a modelStatistics-counter,
  // do NOT delete it!
  int* use_counter;

  if ( !selectParameterTable(parameter_type, table, use_counter) )
    return false;

  int pid = parameter->ID();

  (*table)[pid] = parameter;

  int count = table->size();

  // NOTE!! Success of addition should be checked first! ***!!!***
  (*use_counter)++;

  return true;
}


void
Model::addPendingMeshElements()
{
  int index = 0;
  while (true) {
    BodyElement* be = getBodyElement(index++);
    if (be==NULL) break;
    be->addAllPendingMeshElements();
    be->createExtraVertices();
  }
  
  markActiveObjects();
  updateModelStatistics();
}


void
Model::addSwapElement(int orig_elem_id, int swap_elem_id)
{
  (*modelData->swapElementTable)[orig_elem_id] = swap_elem_id;
}


// Add body element to proper boundary group. Create the group
// if it does not exist.
// Return group's object id
//
int
Model::addToElementGroup(BodyElement* be, int group_tag)
{
  objectType ot = OT_ELEMENT_GROUP;
  objectType bt = be->getObjectType();

  BodyElementGroup* beg = (BodyElementGroup*)getModelObjectByTag(ot, group_tag);

  if ( beg == NULL ) {
    beg = new BodyElementGroup();
    beg->setObjectType(ot);
    beg->setElementType(bt);
  }

  beg->addElement(be->Id());

  return beg->Id();
}


// Method adds a new point into model.
// Idea is to store points only once.
// This concerns points which are connected to vertices.
// Points used in the representation of curves and surfaces are not
// stored in this table.
// Returns true if new point was added.
bool
Model::addToHashPoints(GcPoint* point)
{
  bool new_point = false;
  bool hash_key_found = false;

  // Test if point is already in the table.
  int key = point->hashKey();
  PointHashTable::iterator itr = modelData->modelPoints->find(key);
  PointList* plist;

  if (itr != modelData->modelPoints->end()) {
    plist = (*itr).second;
    hash_key_found = true;
  }

  new_point = new_point || !hash_key_found;

  // Hash-key not in the table, a completely new point-list is added.
  if (!hash_key_found) {

    plist = new PointList;
    plist->push_front(point);
    modelStatistics->nofPoints++;
    (*modelData->modelPoints)[key] = plist;
  }

  // Hash-key was found, but we have to check if point really exist.
  else {

    bool point_exists = ( NULL != findPoint(plist, point) );

    new_point = new_point || !point_exists;

    if (!point_exists) {

      plist->push_front(point);
      modelStatistics->nofPoints++;
    }
  }

  return new_point;
}


// Method adds a new vertex to the model.
// Point2Vertex-table and model's point hashtable is also updated
int
Model::addVertex(BodyElement* vertex)
{
  GcPoint* vpoint = (GcPoint*)vertex->getGeometry();

  (*modelData->modelPoint2Vertices)[vpoint] = vertex;

  addToHashPoints(vpoint);

  return ++modelStatistics->nofVertices;
}


void
Model::bodySelected(Renderer* renderer, int body_id, int lr_id)
{
  Body* body = getBodyById(body_id);

  if (body == NULL) {
    return;
  }

  modelInfo->selectedBody1Id = body_id;
  modelInfo->selectedLayer1Id = lr_id;
  modelInfo->selectedBody2Id = NO_INDEX;
  modelInfo->selectedLayer2Id = NO_INDEX;
  modelInfo->selectedBodyElementId = NO_INDEX;

  objectDrawingMode dmode = body->getDrawMode();

#if 0
  if (dmode == SELECTED) {
    body->setDrawMode(NORMAL);
  } else {
    body->setDrawMode(SELECTED);
  }
#endif

  UserInterface* gui = theControlCenter->getGui();

  gui->selectBody(body_id, lr_id, NO_INDEX, NO_INDEX, true);
}


void
Model::boundarySelected(Renderer* renderer, int bndr_id,
                        int body1_id, int layer1_id,
                        int body2_id, int layer2_id,
                        bool accept_body_change, bool update_gui)
{
  // Selected body element or boundary group
  ModelObject* obj = getModelObjectById(bndr_id);

  const int* bndr_ids = NULL;
  int nof_bndrs = 0;

  int bndr_id1[1];

  BodyElement* be = NULL;
  BodyElementGroup* beg = NULL;

  switch (obj->getObjectType()) {

  // Objects are normal boundaries
  case OT_FACE:
  case OT_EDGE:
  case OT_VERTEX:
    bndr_id1[0] = bndr_id;
    bndr_ids = bndr_id1;
    nof_bndrs = 1;
    be = getBodyElementById(bndr_id);
    break;

  // Objects are element groups (BC panel)
  case OT_ELEMENT_GROUP:
    beg = (BodyElementGroup*) obj;
    nof_bndrs = beg->getNofElements();
    bndr_ids = beg->getElementIds();

    int be_id = beg->getElementId(0);
    be = getBodyElementById(be_id);
    break;
  }

#if 0
  if ( modelInfo->dimension == ECIF_2D )
    be = getEdge(bndr_id);
  else
    be = getFace(bndr_id);
#endif

  // If a non existing id (this should not really happen!)
  //
  if (be == NULL) return;

  int current_body_ids[2] = {modelInfo->selectedBody1Id, modelInfo->selectedBody2Id};

  int body_ids[2] = {NO_INDEX, NO_INDEX};

  int count, i;
  bool is_selected;

  UserInterface* gui = theControlCenter->getGui();

  int bd1_id = be->getParentId(1);
  int bd2_id = be->getParentId(2);

  if ( body2_id == NO_INDEX ) {

    // Element is an inner boundary
    if ( be->isInnerBoundary() ) {

      // inner boundaries are stored the smaller body id
      // always as the first, so we have to check which id
      // is for the body2
      if (body1_id == bd1_id)
        body2_id = bd2_id;
      else
        body2_id = bd1_id;
    }

    // Outer boundary (we make body2_id = body1_id)
    else {
      body2_id = body1_id;
    }
  }

  body_ids[0] = body1_id;
  body_ids[1] = body2_id;

  bool same_bodies = ( current_body_ids[0] == body_ids[0] &&
                       current_body_ids[1] == body_ids[1]
                     );

  // Reset current bodies (body-pair) if we have changed
  // the body
  // NOTE: We don't reset current body if it remains the
  // same, because we can have multiple elements selected
  // in a body
  if ( accept_body_change ||
       same_bodies ||
       current_body_ids[0] == NO_INDEX
     ) {
    count = 0;

  } else if ( current_body_ids[0] == current_body_ids[1] ||
              current_body_ids[1] == NO_INDEX
            ) {
    count = 1;

  } else {
    count = 2;
  }

  // Reset current bodies if changed
  for (i = 0; i < count; i++) {

    int current_body_id = current_body_ids[i];
    Body* current_body = getBodyById(current_body_id);

    if ( current_body == NULL ) continue;

    int index = 0;
    while (true) {
      BodyElement* cbe = current_body->getElement(index++);
      if (cbe==NULL) break;
      gui->setBoundarySelectionMode(cbe->Id(), false);
      cbe->setDrawState(DS_NORMAL);
    }
  }

  // Reset current boundary selections
  resetBoundarySelections(update_gui, beg != NULL, nof_bndrs, bndr_ids, false);

  // Update current ids
  current_body_ids[0] = body1_id;
  current_body_ids[1] = body2_id;

  modelInfo->selectedBody1Id = body1_id;
  modelInfo->selectedBody2Id = body2_id;

  modelInfo->selectedBodyElementId = bndr_id;

  for (i = 0; i < nof_bndrs; i++) {

    be = getBodyElementById(bndr_ids[i]);

    // Pick selected bodyelement state
    objectDrawingState dstate = be->getDrawState();

    // Set state for the current
    if (dstate == DS_SELECTED) {

      be->setDrawState(DS_NORMAL);
      is_selected = false;

    } else {

      be->setDrawState(DS_SELECTED);
      is_selected = true;
    }
  }

  // Set element's selection state in the Gui
  if (update_gui) {

    gui->setBoundarySelectionMode(bndr_id, is_selected, true);

    if ( body1_id == body2_id ) {
      body2_id = NO_INDEX;
    }

    gui->selectBoundary(bndr_id, body1_id, layer1_id, body2_id, layer2_id, modelFlags[SELECT_OBJECTS_EXTEND]);
  }

  gui->markSelectedBoundaries();
}


double
Model::calcInitialMeshH()
{
  double initial_mesh_h = MAX_RANGE;

  if (modelInfo->dimMin > 0 && modelInfo->dimMin < initial_mesh_h)
    initial_mesh_h = MESH_H_INIT_FACTOR * modelInfo->dimMin;

  if (modelInfo->minEdgeSize > 0 && modelInfo->minEdgeSize < initial_mesh_h)
    initial_mesh_h = modelInfo->minEdgeSize;

  return initial_mesh_h;

}


// Method checks that all model's bodies are ok.
// Return true if there are no bodies or if all are ok.
//
bool
Model::checkBodies()
{
  // A table (a stack) for available default colors.
  ColorIndexArray  color_table;

  for(int i = MAX_NOF_BODIES-1; i >= 0; i--)
    color_table.push_back(defaultColorIndices[i]);

  RangeVector rv;
  colorIndices body_color;
  ColorIndexArray::iterator pos;

  int index = 0;
  while (true) {

    Body* body = getBody(index++);

    if (body==NULL) break;

    body->check();

    if ( !body->isOk() ) {
      return false;
    }

    // If body is ok, update model's bounding box.
    body->getRangeVector(rv);
    modelBox->extendByRange(rv);

    // Remove used non-default color from unused-colors.
    body_color = body->getColorIndex();

    if (body_color != DEFAULT_COLOR_INDEX) {

		  pos = std::find(color_table.begin(), color_table.end(), body_color);

      if (pos != color_table.end())
        color_table.erase(pos);
    }
  }

  // Replace non-defined (=DEFAULT_COLOR_INDEX) with
  // colors from default-indices list.
  setBodyColors(color_table);

  // Set bounding box info
  modelBox->getRangeVector(rv);

  modelInfo->minX = rv[0];
  modelInfo->maxX = rv[1];
  modelInfo->minY = rv[2];
  modelInfo->maxY = rv[3];
  modelInfo->minZ = rv[4];
  modelInfo->maxZ = rv[5];

  // Model dimensions
  modelInfo->dimX = modelInfo->maxX - modelInfo->minX;
  modelInfo->dimY = modelInfo->maxY - modelInfo->minY;
  modelInfo->dimZ = modelInfo->maxZ - modelInfo->minZ;

  modelInfo->dimAvg = modelInfo->dimX + modelInfo->dimY + modelInfo->dimZ;
  modelInfo->dimAvg /= modelInfo->dimension;

  modelInfo->dimMax = modelInfo->dimX;
  if (modelInfo->dimY > modelInfo->dimMax)
    modelInfo->dimMax = modelInfo->dimY;
  if (modelInfo->dimension == ECIF_3D && modelInfo->dimZ > modelInfo->dimMax )
    modelInfo->dimMax = modelInfo->dimZ;

  modelInfo->dimMin = modelInfo->dimX;
  if (modelInfo->dimY < modelInfo->dimMin)
    modelInfo->dimMin = modelInfo->dimY;
  if (modelInfo->dimension == ECIF_3D && modelInfo->dimZ < modelInfo->dimMin )
    modelInfo->dimMin = modelInfo->dimZ;

  if ( modelInfo->modelSourceType != ECIF_MESH_FILE &&
       modelInfo->dimAvg > 0.0
     )
    setFlagValue(GEOMETRY_TYPE, GEOMETRY_TYPE_CAD, true);
  else
    setFlagValue(GEOMETRY_TYPE, GEOMETRY_TYPE_CAD, false);

  return true;
}


// Method checks that all model's bodyelement loops are correct.
// 2D: ccw-oriented
// 3D: nothing special currently
// Return true if no loops or if all could be coreclty oriented
//
// NOTE: An error (in 2D) means that loop is not continuous!
//
bool
Model::checkBodyElementLoops()
{
  int index = 0;
  while (true) {
    BodyElementLoop* bel = getBodyElementLoop(index++);
    if (bel==NULL) break;
    if ( !bel->check() ) {
      return false;
    }
  }

  return true;
}


bool
Model::checkBodyElements()
{
  int index = 0;
  while (true) {
    BodyElement* be = getBodyElement(index++);
    if (be==NULL) break;
    be->addAllPendingSubElements();
  }

  return true;
}


bool
Model::checkBoundaries()
{
  // Check outer boundaries which possibly have "hanging"
  // components to be created
  int index = 0;
  while (true) {
    BodyElement* be = getBoundary(index++);
    if (be==NULL) break;
    if ( !be->checkOuterBoundaries() ) {
      return false;
    }
  }

  // Next remove splitted boundaries from bodies
  // and replace them with corresponding covering
  // subelements
  index = 0;
  while (true) {
    Body* body = getBody(index++);
    if (body==NULL) break;
    if ( !body->checkElements() ) {
      return false;
    }
  }

  // Next remove splitted boundaries from element loops
  // and replace them with corresponding covering
  // subelements
  index = 0;
  while (true) {
    BodyElementLoop* bel = getBodyElementLoop(index++);
    if (bel==NULL) break;
    if ( !bel->checkElements() ) {
      return false;
    }
  }

  return true;
}


bool
Model::checkElementGroupData()
{
  int index, be_id;
  BodyElementGroup* beg;

  // Check boundary  groups
  // ----------------------
  index = 0;
  while (true) {
    BodyElementGroup* beg = getBodyElementGroup(index++);
    if (beg==NULL) break;
    if ( !beg->check() ) return false;
  }

  // Check all boundaries
  // --------------------
  index = 0;
  while (true) {
    BodyElement* be = getFace(index++);
    if (be==NULL) break;
    if (!be->checkElementGroupData()) return false;
  }

  index = 0;
  while (true) {
    BodyElement* be = getEdge(index++);
    if (be==NULL) break;
    if (!be->checkElementGroupData()) return false;
  }

  index = 0;
  while (true) {
    BodyElement* be = getVertex(index++);
    if (be==NULL) break;
    if (!be->checkElementGroupData()) return false;
  }

  return true;
}


void
Model::checkDiffuseGrayRadiation()
{
  modelInfo->hasDiffuseGrayRadiation = false;

  int index = 0;
  while (true) {
    Parameter* bc = getParameter(index++, ECIF_BOUNDARY_CONDITION);
    if (bc==NULL) break;
    if ( bc->IsApplied() &&
         bc->hasFieldValueBySifName("Radiation", "Diffuse gray")
       ) {
      modelInfo->hasDiffuseGrayRadiation = true;
      return;
    }
  }
}


// Method checks that all model's mesh bodies are ok.
// Return true if there are no bodies or if all are ok.
bool
Model::checkMeshBodies()
{
  int index = 0;
  while (true) {
    Body* bd = getBody(index++);
    if (bd==NULL) break;
    bd->calcMeshBox();
  }

  MeshInfo* mi = meshInfo;

  //-----Update model's mesh bounding box
  RangeVector rv;
  rv[0] = mi->minX;
  rv[1] = mi->maxX;
  rv[2] = mi->minY;
  rv[3] = mi->maxY;
  rv[4] = mi->minZ;
  rv[5] = mi->maxZ;

  meshBox->extendByRange(rv);

  mi->dimX = mi->maxX - mi->minX;
  mi->dimY = mi->maxY - mi->minY;
  mi->dimZ = mi->maxZ - mi->minZ;

  mi->dimAvg = mi->dimX + mi->dimY + mi->dimZ;
  mi->dimAvg /= modelInfo->dimension;

  double pick_tol = 0.05 * mi->dimAvg;
  MeshEdgeElementTable::setPickingTolerance(pick_tol);

  if (meshInfo->dimAvg < 0.0)
    setFlagValue(GEOMETRY_TYPE, GEOMETRY_TYPE_MESH, false);
  else
    setFlagValue(GEOMETRY_TYPE, GEOMETRY_TYPE_MESH, true);

  return true;
}



void
Model::checkVertexExistence()
{
  int index = 0;
  while (true) {
    BodyElement* be = getBodyElement(index++);
    if (be==NULL) break;
    be->createExtraVertices();
  }

  markActiveObjects();
  updateModelStatistics();
}

// Check if a verex exists 
bool
Model::checkVertexExistence(int tag)
{
  if ( NULL == getVertexByTag(tag) ) 
    return false;
  else
    return true;
}


ecif_modelStatus
Model::checkStatus()
{
  ecif_modelStatus modelStatus = STATUS_OK;

  int index = 0;
  while (true) {
    Body* body = getBody(index++);
    if (body==NULL) break;
    ecif_modelStatus body_status = body->checkStatus();
    modelStatus |= body_status;
  }

  return modelStatus;
}


// Combine mesh boundaries (in one body or body-pair)
//
void
Model::combineBoundaries(int body1_id, int body2_id)
{
  int i,j;
  UserInterface* gui = theControlCenter->getGui();

  Body* body1 = getBodyById(body1_id);
  Body* body2 = getBodyById(body2_id);

  int body1_tag = NO_INDEX;
  int body2_tag = NO_INDEX;

  int body1_layer = NO_INDEX;
  int body2_layer = NO_INDEX;

  if ( body1 != NULL ) {
    body1_tag = body1->Tag();
  }

  if ( body2 != NULL ) {
    body2_tag = body2->Tag();
  }

  int* source_ids = new int[modelStatistics->nofElements];
  int nof_source_ids = 0;

  //---Count nof mesh elements in the new body element
  //   and the number of contributing boundaries
  int nof_mesh_elements = 0;

  int index = 0;
  while (true) {
    BodyElement* be = getBoundary(index++);
    if (be==NULL) break;
    if ( DS_SELECTED == be->getDrawState() ) {
      source_ids[nof_source_ids++] = be->Id();
      nof_mesh_elements += be->getNofMeshElements();
      body1_layer = be->getParentLayer(1);
      body2_layer = be->getParentLayer(2);
    }
  }

  //---If nothing or only one element was selected
  if (nof_source_ids <= 1) {
    if (nof_source_ids == 0)
      gui->showMsg("No boundaries selected for combine!");
    else
      gui->showMsg("Only one boundary selected. Nothing to combine!");

    delete[] source_ids;

    return;
  }

  // If this is the first edit operation
  if ( modelInfo->editingMeshBoundaries == false ) {
    startEditMeshBoundaries();
  }

  // Delete old info
  short idx = modelData->splitCombineInfoIndex;
  delete (*modelData->splitCombineInfos)[idx];

  // Create new info entry and store it
  SplitCombineInfo* sc_info = new SplitCombineInfo;

  (*modelData->splitCombineInfos)[idx] = sc_info;

  // Update split/combine index for next insert
  setNewSplitCombineIndex(+1);

  //---Create new bodyelement, add to the model and to the bodies
  BodyElement* new_be = createBodyElement(body1_tag, body1_layer, body2_tag, body2_layer, nof_mesh_elements);

  new_be->setParentIds(body1_id, body2_id);
  new_be->setParentTags(body1_tag, body2_tag);
  new_be->checkElementGroupData();

  (*modelData->createdModelElements)[new_be->Id()] = new_be;

  // Set info data
  sc_info->canRedo = false;
  sc_info->canUndo = true;
  sc_info->combineTargetId = new_be->Id();
  sc_info->nofCombineSourceIds = nof_source_ids;
  sc_info->combineSourceIds = new int[nof_source_ids];

  //---Copy selected mesh elements to the new boundary
  //   and remove the original bodyelements from model
  for (i = 0; i < nof_source_ids; i++) {

    sc_info->combineSourceIds[i] = source_ids[i];

    BodyElement* be = getBoundaryById(source_ids[i]);

    int nof_elements = be->getNofMeshElements();

    //-Copy all subelements
    for (j = 0; j < be->getNofSubElements(); j++) {
      BodyElement* se = be->getSubElement(j);
      new_be->addSubElement(se);
    }

    //-Copy all mesh elements
    for (int j = 0; j < nof_elements; j++) {
      int elem_id = be->getMeshElementId(j);
      short elem_dir = be->getMeshElementDir(j);

      // Add mesh element to the new boundary
      new_be->addMeshElement(elem_id, elem_dir);
      meshData->boundaryElementBoundaryIds[elem_id] = sc_info->combineTargetId;
    }

    //-Remove boundary element from the model, but not the
    // subelements!
    //
    removeBodyElement(be, body1, body2, false);

  } // for source ids

  delete[] source_ids;

  body1->check();

  if (body2 != NULL) {
    body2->check();
  }

  // Update model
  new_be->findMeshBorder();
  updateModelStatistics();

  gui->updateBodyData(this);
  gui->updateBoundaryData(this);

  gui->setMeshEdited();

  refreshRenderer();
}



bool
Model::convertTags2Ids()
{
  int index;

  // Bodies
  // ======
  index = 0;
  while (true) {
    Body* body = getBody(index++);
    if (body==NULL) break;
    if ( !body->convertTags2Ids() ) return false;
  }

  // Boundary loop
  // =============
  index = 0;
  while (true) {
    BodyElementLoop* bel = getBodyElementLoop(index++);
    if (bel==NULL) break;
    if ( !bel->convertTags2Ids() ) return false;
  }

  // Boundaries
  // =============
  index = 0;
  while (true) {
    BodyElement* be = getBoundary(index++);
    if (be==NULL) break;
    if ( !be->convertTags2Ids() ) return false;
  }

  return true;
}


// Creates one mesh geometry body element with given element-id
BodyElement*
Model::createBodyElement(int element_tag,
                         int body1_tag, int body1_layer,
                         int body2_tag, int body2_layer,
                         int nof_mesh_elements)
{
  Body* body1 = getBodyByTag(body1_tag);
  Body* body2 = getBodyByTag(body2_tag);

  if (body1 == NULL ) {
    return NULL;
  }

  BodyElement* be = NULL;

  switch (modelInfo->dimension) {
  case ECIF_2D:
    be = new BodyElement2D(element_tag, body1_tag, body2_tag, nof_mesh_elements);
    break;
  case ECIF_3D:
    be = new BodyElement3D(element_tag, body1_tag, body2_tag, nof_mesh_elements);
    break;
  }

  if (be != NULL) {
    be->setParentLayers(body1_layer, body2_layer);
  }

  // No success (nofElements still zero!)
  if ( 0 == addBodyElement(be, body1_tag, body1_layer, body2_tag, body2_layer) ) {
    delete be;
    be = NULL;
  }

  return be;
}


// Creates one mesh geometry body element
BodyElement*
Model::createBodyElement(int body1_tag, int body1_layer, int body2_tag, int body2_layer, int nof_mesh_elements)
{
  Body* body1 = getBodyByTag(body1_tag);
  Body* body2 = getBodyByTag(body2_tag);

  if (body1 == NULL ) {
    return NULL;
  }

  BodyElement* be = NULL;

  switch (modelInfo->dimension) {
  case ECIF_2D:
     be = new BodyElement2D(body1_tag, body2_tag, nof_mesh_elements);
    break;
  case ECIF_3D:
     be = new BodyElement3D(body1_tag, body2_tag, nof_mesh_elements);
    break;
  }

  if (be != NULL) {
    be->setParentLayers(body1_layer, body2_layer);
  }

  if ( 0 == addBodyElement(be, body1_tag, body1_layer, body2_tag, body2_layer) ) {
    delete be;
    be = NULL;
  }

  return be;
}


void
Model::createMeshBoundaryElementEdges()
{
  //---3D
  if ( modelInfo->dimension == ECIF_3D ) {

    // Boundary edges
    meshManager->findMeshElementEdges(meshData->boundaryElements, meshInfo->nofBoundaryEdges);

    delete meshData->boundaryEdges;
    meshData->boundaryEdges = new MeshEdgeElementTable(meshInfo->nofBoundaryEdges, false);

    meshManager->createMeshElementEdges(meshData->boundaryElements, meshData->boundaryEdges);

#if 1
    // Boundary vertices
    meshManager->findMeshElementEdges(meshData->boundaryEdges, meshInfo->nofBoundaryVertices);

    delete meshData->boundaryVertices;
    meshData->boundaryVertices = new MeshVertexElementTable(meshInfo->nofBoundaryVertices);

    meshManager->createMeshElementEdges(meshData->boundaryEdges, meshData->boundaryVertices);
#endif

  //---2D
  } else {

    // Boundary vertices
    meshManager->findMeshElementEdges(meshData->boundaryElements, meshInfo->nofBoundaryVertices);

    delete meshData->boundaryVertices;
    meshData->boundaryVertices = new MeshVertexElementTable(meshInfo->nofBoundaryVertices);

    meshManager->createMeshElementEdges(meshData->boundaryElements, meshData->boundaryVertices);
  }

}


void
Model::createMeshBulkElementEdges()
{
  meshManager->findMeshElementEdges(meshData->bulkElements, meshInfo->nofBulkEdges);

  delete[] meshData->bulkEdgeRendered;
  meshData->bulkEdgeRendered = new bool[meshInfo->nofBulkEdges];

  UserInterface* gui = theControlCenter->getGui();

  strstream strm;
  strm << "Nof bulk edges: " << meshInfo->nofBulkEdges << ends;
  gui->showMsg(strm.str());
}


const char*
Model::objectType2Name(objectType type, int max_buf_len, char* name_buffer)
{
  char* name;

  switch (type) {
  case OT_NONE:
    name = "None"; break;
  case OT_BODY:
    name = "Body"; break;
  case OT_BODY_LAYER:
    name = "Body Layer"; break;
  case OT_BODYPAIR:
    name = "Body Pair"; break;
  case OT_ELEMENT_LOOP:
    name = "Element Loop"; break;
  case OT_BOUNDARY:
    name = "Boundary"; break;
  case OT_FACE:
    name = "Face"; break;
  case OT_EDGE:
    name = "Edge"; break;
  case OT_VERTEX:
    name = "Vertex"; break;
  case OT_ELEMENT_GROUP:
    name = "Element Group"; break;
  default:
    name = "None"; break;
  }

  if ( max_buf_len > 0  && name_buffer != NULL ) {
    strncpy(name_buffer, name, max_buf_len);
  }

  return name;
}


objectType
Model::objectName2Type(const char* name)
{
  if ( name == NULL || LibFront::in(name, "None") )
    return OT_NONE;
  else if (LibFront::in(name, "Body") )
    return OT_BODY;
  else if (LibFront::in(name, "Body Layer") )
    return OT_BODY_LAYER;
  else if (LibFront::in(name, "Body Pair") )
    return OT_BODYPAIR;
  else if (LibFront::in(name, "Element Loop") )
    return OT_ELEMENT_LOOP;
  else if (LibFront::in(name, "Boundary") )
    return OT_BOUNDARY;
  else if (LibFront::in(name, "Face") )
    return OT_FACE;
  else if (LibFront::in(name, "Edge") )
    return OT_EDGE;
  else if (LibFront::in(name, "Vertex") )
    return OT_VERTEX;
  else if (LibFront::in(name, "Element Group") )
    return OT_ELEMENT_GROUP;
  else
    return OT_NONE;
}


ostream&
Model::outputSolverTargetFields_sif(ostream& out, short indent_size, short indent_level, const char* source_eq_name)
{
  return outputManager->sif_outputSolverTargetFields(out, indent_size, indent_level, source_eq_name);

}



Parameter*
Model::createNewParameter(ecif_parameterType param_type,
                          int pid, int parent_id,
                          char* param_value, char* param_name)
{

  switch (param_type) {
  case ECIF_BODY_FORCE:
    return new BodyForce(pid, parent_id, param_value, param_name);
  case ECIF_BODY_PARAMETER:
    return new BodyParameter(pid, parent_id, param_value, param_name);
  case ECIF_BOUNDARY_CONDITION:
    return new BoundaryCondition(pid, parent_id, param_value, param_name);
  case ECIF_BOUNDARY_PARAMETER:
    return new BoundaryParameter(pid, parent_id, param_value, param_name);
  case ECIF_CALCULATOR:
    return new Calculator(pid, param_value, param_name);
  case ECIF_CONSTANT:
    return new Constant(pid, param_value, param_name);
  case ECIF_COORDINATE:
    return new Coordinate(pid, param_value, param_name);
  case ECIF_DATAFILE:
    return new Datafile(pid, param_value, param_name);
  case ECIF_EQUATION:
    return new Equation(pid, param_value, param_name);
  case ECIF_EQUATION_VARIABLE:
    return new EquationVariable(pid, param_value, param_name);
  case ECIF_GRID_H:
    return new GridH(pid, parent_id, param_value, param_name);
  case ECIF_GRID_PARAMETER:
    return new GridParameter(pid, parent_id, param_value, param_name);
  case ECIF_INITIAL_CONDITION:
    return new InitCondition(pid, parent_id, param_value, param_name);
  case ECIF_MATERIAL:
    return new Material(pid, parent_id, param_value, param_name);
  case ECIF_MODEL_PARAMETER:
    return new ModelParameter(pid, param_value, param_name);
  case ECIF_SIMULATION_PARAMETER:
    return new SimulationParameter(pid, param_value, param_name);
  case ECIF_SOLVER:
    return new Solver(pid, param_value, param_name);
  case ECIF_SOLVER_CONTROL:
    return new SolverControl(pid, param_value, param_name);
  case ECIF_TIMESTEP:
    return new Timestep(pid, param_value, param_name);
  case ECIF_USER_SETTING:
    return new UserSetting(pid, param_value, param_name);
  default:
    return NULL;
  }
}


Parameter*
Model::createNewParameter(ecif_parameterType parameter_type, int p_id)
{
  Parameter* parameter;

  switch (parameter_type) {
  case ECIF_BODY_FORCE:
    parameter = new BodyForce(p_id);
    break;
  case ECIF_BODY_PARAMETER:
    parameter = new BodyParameter(p_id);
    break;
  case ECIF_BOUNDARY_CONDITION:
    parameter = new BoundaryCondition(p_id);
    break;
  case ECIF_BOUNDARY_PARAMETER:
    parameter = new BoundaryParameter(p_id);
    break;
  case ECIF_CALCULATOR:
    parameter = new Calculator(p_id);
    break;
  case ECIF_CONSTANT:
    parameter = new Constant(p_id);
    break;
  case ECIF_COORDINATE:
    parameter = new Coordinate(p_id);
    break;
  case ECIF_DATAFILE:
    parameter = new Datafile(p_id);
    break;
  case ECIF_EQUATION:
    parameter = new Equation(p_id);
    break;
  case ECIF_EQUATION_VARIABLE:
    parameter = new EquationVariable(p_id);
    break;
  case ECIF_GRID_H:
    parameter = new GridH(p_id);
    break;
  case ECIF_GRID_PARAMETER:
    parameter = new GridParameter(p_id);
    break;
  case ECIF_INITIAL_CONDITION:
    parameter = new InitCondition(p_id);
    break;
  case ECIF_MATERIAL:
    parameter = new Material(p_id);
    break;
  case ECIF_MODEL_PARAMETER:
    parameter = new ModelParameter(p_id);
    break;
  case ECIF_SIMULATION_PARAMETER:
    parameter = new SimulationParameter(p_id);
    break;
  case ECIF_SOLVER:
    parameter = new Solver(p_id);
    break;
  case ECIF_SOLVER_CONTROL:
    parameter = new SolverControl(p_id);
    break;
  case ECIF_TIMESTEP:
    parameter = new Timestep(p_id);
    break;
  case ECIF_USER_SETTING:
    parameter = new UserSetting(p_id);
    break;
  default:
    parameter = NULL;
    break;
  }

  parameter->checkLastId(p_id);

  parameter->setChangedState(false);

  return parameter;
}


void
Model::deleteMeshData()
{
  int index = 0;

  setFlagValue(GEOMETRY_TYPE, GEOMETRY_TYPE_MESH, false);
  setFlagValue(DRAW_SOURCE, DRAW_SOURCE_MESH, false);

  theControlCenter->getGui()->updateModelFlags(this);
  theControlCenter->getGui()->configureButtons("draw_source_mesh", 0);

  delete meshData; meshData = NULL;
  delete meshInfo; meshInfo = NULL;
  delete meshBox;  meshBox  = NULL;

  index = 0;
  while(true) {
    Body* body = getBody(index++);
    if (body==NULL) break;
    body->deleteMeshElements();
  }

  index = 0;
  while(true) {
    BodyElement* be = getBoundary(index++);
    if (be==NULL) break;
    be->deleteMeshElements();
  }
}


// Delete all parameters from the table defined
// by the parameter_type
void
Model::deleteParameters(ecif_parameterType parameter_type)
{
  ParameterTable* table;

  // NOTE: This will be scalar pointer to a modelStatistics-counter,
  // do NOT delete it!
  int* use_counter;

  if ( !selectParameterTable(parameter_type, table, use_counter) ) {
    return;
  }

  table->clear();
  (*use_counter) = 0;

  // Reset objects' parameter counters
  switch (parameter_type) {
  case ECIF_EQUATION:
    modelStatistics->nofBodiesWithEquation = 0;
    break;
  case ECIF_MATERIAL:
    modelStatistics->nofBodiesWithMaterial = 0;
    break;
  case ECIF_BODY_FORCE:
    modelStatistics->nofBodiesWithBodyForce = 0;
    break;
  case ECIF_INITIAL_CONDITION:
    modelStatistics->nofBodiesWithInitialCondition  = 0;
    break;
  case ECIF_BOUNDARY_CONDITION:
    modelStatistics->nofInnerBoundariesWithCondition = 0;
    modelStatistics->nofOuterBoundariesWithCondition = 0;
    break;
  }

}


// Delete one parameter by id from the table defined
// by the parameter_type
int
Model::deleteParameter(ecif_parameterType parameter_type, int pid)
{
  ParameterTable* table;

  // NOTE: This will be scalar pointer to a modelStatistics-counter,
  // do NOT delete it!
  int* use_counter;

  if ( !selectParameterTable(parameter_type, table, use_counter) ) {
    return NO_INDEX;
  }

  if ( 0 == table->erase(pid))
    return NO_INDEX;

  return --(*use_counter);
}


void
Model::deleteCadData()
{
  UserInterface* gui = theControlCenter->getGui();

  setFlagValue(GEOMETRY_TYPE, GEOMETRY_TYPE_CAD, false);
  setFlagValue(DRAW_SOURCE, DRAW_SOURCE_CAD, false);

  gui->updateModelFlags(this);
  gui->configureButtons("draw_source_cad", 0);

  int index = 0;
  while(true) {
    BodyElement* be = getBoundary(index++);
    if (be==NULL) break;
    be->deleteCadData();
  }

}


void
Model::fieldNameGuiToSif(const char* gui_name, char* sif_name_buffer)
{
  theControlCenter->getUI()->fieldNameGuiToSif(gui_name, sif_name_buffer);
}


void
Model::fieldNameSifToGui(const char* sif_name, char* gui_name_buffer)
{
  theControlCenter->getUI()->fieldNameSifToGui(sif_name, gui_name_buffer);
}


// Find body pairs from body element parents
//
void
Model::findBodyPairs()
{

  int be_id, bd1_id, bd2_id;

  // Loop all elements
  int index = 0;
  while (true) {

    BodyElement* be = getBoundary(index++);

    if (be==NULL) break;

    int bd1_id = be->getParentId(1);
    int bd2_id = be->getParentId(2);

    if ( bd1_id == NO_INDEX || bd2_id == NO_INDEX ) {
      continue;
    }

    // This is a boundary between connected body layers!
    // It does not create a body pair!!!
    //
    if ( bd1_id == bd2_id ) {
      continue;
    }

    Body* body1 = getBodyById(bd1_id);
    Body* body2 = getBodyById(bd2_id);

    // Virtual bodies cannot be a pair member
    //
    if ( body1->isVirtual() || body2->isVirtual() ) {
      continue;
    }

    BodyPair* bp = getBodyPairById(body1, body2);

    // Create new bodypair
    //
    if ( bp == NULL ) {
      bp = new BodyPair(body1, body2);
      bp->setTag();
      bp->setId();
      addBodyPair(bp);
    }

  }
}


// Method finds neigbouring bodies and innerBoundaries between these
// bodies.
void
Model::findBoundaries()
{
  // Find all candiate neigbours (body pairs).
  modelData->neighbourCandidates = findNeighbourCandidates();

  // And find all candidate innerBoundaries beteween candidate bodies.
  modelData->boundaryCandidates = findBoundaryCandidates();

  // And find real neigobour bodies and adjacent innerBoundaries.
  findInnerBoundaries();

  checkBoundaries();
}


// Method finds all potential adjacent innerBoundaries beteween potential
// neighbour bodies.
AdjacentPairArray*
Model::findBoundaryCandidates()
{
  // This method creates a table where we have pairs of (Body,Element)-pairs.
  // Each pair is a candidate for an adjancent surface because bodyelement
  // pairs are overlapping (this is a necessary but not a sufficient
  // condition for surfaces being adjacent).
  AdjacentPairArray* cp_array = new AdjacentPairArray;
  BodyElement* be1;
  BoundBox* be1_box;
  BodyElement* be2;
  BoundBox* be2_box;
  BoundBox* body2_box;
  BodyElement** be_array1;
  BodyElement** be_array2;
  int be_count1, be_count2;

  int bp_count = modelData->neighbourCandidates->size();
  BodyPair* bp;

  //We read each pair from *bp_array* and check if the elements in the bodies
  //overlap.
  for (int i = 0; i < bp_count; i++) {
    bp = (*modelData->neighbourCandidates)[i];

    if (bp == NULL)
      continue;

    const Body* bd1;
    const Body* bd2;

    bp->getBodies(bd1, bd2);

    Body* body1 = (Body*) bd1;
    Body* body2 = (Body*) bd2;

    int layer1 = -1;
    while (true) {

      if (!body1->selectLayer(++layer1)) break;

      // Elements from the body1-layer into a table.
      be_array1 = body1->getElementArray(layer1);
      be_count1 = body1->getElementCount(layer1);

      int layer2 = -1;
      while (true) {

        if (!body2->selectLayer(++layer2)) break;

        // Elements from the body2-layer into a table.
        be_array2 = body2->getElementArray();
        be_count2 = body2->getElementCount(layer2);

        //Next we loop all elements of the body1 and check if they overlap
        //with elements in the body2.
        body2_box = body2->getBoundBox(layer2); //each be1 must overlap this!
        for (int j = 0; j < be_count1; j++) {

          be1 = be_array1[j];
          be1_box = be1->getBoundBox();

          if (! be1_box->overlap(body2_box)) continue;

          for (int k = 0; k < be_count2; k++) {

            be2 = be_array2[k];
            be2_box = be2->getBoundBox();

            if (! be1_box->overlap(be2_box)) continue;

            AdjacentHalf** adj_pair = new AdjacentHalf*[2];
            adj_pair[0] = new AdjacentHalf;
            adj_pair[1] = new AdjacentHalf;
            adj_pair[0]->body = body1;
            adj_pair[0]->bodyLayer = layer1;
            adj_pair[0]->element = be1;
            adj_pair[1]->body = body2;
            adj_pair[1]->bodyLayer = layer2;
            adj_pair[1]->element = be2;
            cp_array->push_back(adj_pair);
          }
        }

        delete[] be_array2;
      }

      delete[] be_array1;
    }
  }

  //delete[] b_array;
  return cp_array;
}

// Method finds body layers which are completely inside an other body layer.
// The outer boundary of the contained layer becomes an inner boundary!
void
Model::findContainedBodies()
{
  //Read bodies into an array.
  int b_count = modelStatistics->nofBodies;

  Body** b_array = new Body*[b_count];

  int row = 0;
  int index = 0;
  while (true) {

    Body* body = getBody(index++);
    if (body==NULL) break;

    if ( body->isVirtual() ) {
      b_count--;
      continue;
    }

    b_array[row++] = body;
  }

  // Loop each body1
  for (int i = 0; i < b_count; i++) {

    // Loop each body2
    for (int j = i+1; j < b_count; j++) {

      Body* body1 = b_array[i];
      Body* body2 = b_array[j];

      // If this is already an existing body-pair, do nothing!
      if ( NULL != getBodyPairById(body1, body2) ) {
        continue;
      }

      int layer1 = -1;
      int layer2 = -1;

      // Each body1 layer
      //
      while (true) {

        if (!body1->selectLayer(++layer1)) break;

        // Each body2 layer
        while (true) {

          if (!body2->selectLayer(++layer2)) break;

          Body* bd1 = NULL;
          Body* bd2 = NULL;
          int lr1 = NO_INDEX;
          int lr2 = NO_INDEX;

          // Select "outer" and "inner" bodies (layers)
          if ( body1->hasInside(layer1, body2, layer2) ) {
            bd1 = body1;
            lr1 = layer1;
            bd2 = body2;
            lr2 = layer2;
          } else if ( body2->hasInside(layer2, body1, layer1) ) {
            bd1 = body2;
            lr1 = layer2;
            bd2 = body1;
            lr2 = layer1;
          }

          // Layers do not form any inner boundaries
          if ( bd1 == NULL || bd2 == NULL ) continue;

          BodyElementLoop* bel = bd2->getElementLoop(lr2);

          // Insert bd2 outer boundary as an inner loop
          // into the bd1
          int nof_loop_ids1 = bd1->getNofElementLoops(lr1);
          int* new_loop_ids = new int[1 + nof_loop_ids1];
          int* new_loop_tags = new int[1 + nof_loop_ids1];

          new_loop_ids[0] = bd1->getElementLoopId(lr1, 0);
          new_loop_ids[1] = -1 * bel->Id();

          new_loop_tags[0] = bd1->getElementLoopTag(lr1, 0);
          new_loop_tags[1] = bel->Tag();

          for(int i = 1; i < nof_loop_ids1; i++) {
            new_loop_ids[1 + i] = bd1->getElementLoopId(lr1, i);
            new_loop_tags[1 + i] = bd1->getElementLoopTag(lr1, i);
          }

          bd1->setElementLoops(lr1, 1 + nof_loop_ids1, new_loop_ids, new_loop_tags);

          // Add all elements in the loop to the bd1 and marks it
          // as the second parent
          int index = 0;
          while (true) {
            BodyElement* be = bel->getElement(index++);
            if (be==NULL) break;
            int bd1_id = bd1->Id();
            int bd2_id = bd2->Id();

            if ( bd1_id < bd2_id ) {
              be->setParentIds(bd1_id, bd2_id);
            } else {
              be->setParentIds(bd2_id, bd1_id);
            }

            bd1->addElement(lr1, be);
          }

        } // Each body2 layer
      } // Each body1 layer
    } // Body2 loop
  } // Body1 loop
}


// Method checks which of the candidate pair of bodyelelements
// really are adjacent and creates bodypairs
// NOTE: body-pairs must be 'lexicographically' ordered in
// *boundaryCandidates*, otherwise *neigbours* is not correctly created.
void
Model::findInnerBoundaries()
{
  IdArray swap1_ids;
  IdArray swap2_ids;
  IdArray relative_dirs;

  int count = modelData->boundaryCandidates->size();

  for (int i = 0; i < count; i++) {

    AdjacentHalf** adj_pair = (*modelData->boundaryCandidates)[i];

    if (adj_pair == NULL)
      continue;

    BodyElement* be1 = adj_pair[0]->element;
    BodyElement* be2 = adj_pair[1]->element;
    BodyElement* common = NULL;

    matchType match_type = be1->findCommonBoundary(be2, common);

    if (match_type != MATCH_NONE) {
      //new body-pair and boundary structures are created.

      // Bodies and layers:
      Body* b1 = adj_pair[0]->body;
      Body* b2 = adj_pair[1]->body;
      int lr1  = adj_pair[0]->bodyLayer;
      int lr2  = adj_pair[1]->bodyLayer;

      b1->incrBndrCount(lr1, INNER_ELEMENT);
      b2->incrBndrCount(lr2, INNER_ELEMENT);

      // info concerning the common boundary is added to both bodyelements.
      switch (match_type) {
      case MATCH_1_INSIDE:
        be2->addCoveringElement(common, BE_INNER);
        break;
      case MATCH_2_INSIDE:
        be1->addCoveringElement(common, BE_INNER);
        break;
      case MATCH_OVERLAP:
        be1->addCoveringElement(common, BE_INNER);
        be2->addCoveringElement(common, BE_INNER);
        break;
      case MATCH_EXACT:
        be1->addStatus(BE_INNER);
        be2->addStatus(BE_INNER);
        b1->swapBodyElements(lr1, be1, common);
        b2->swapBodyElements(lr2, be2, common);

        if ( be1 == common) {
          swap1_ids.push_back(be2->Id());
        } else {
          swap1_ids.push_back(be1->Id());
        }

        swap2_ids.push_back(common->Id());

        relative_dirs.push_back(getRelativeOrientation(be1, be2));
        break;
      }

    } // if match-type != none

  } // for all boundaries candidates

  if ( swap1_ids.size() > 0 ) {
    swapBodyElements(&swap1_ids, &swap2_ids, &relative_dirs);
  }

}


// Method finds all body-pairs which potentially are adjacent neighbour bodies.
BodyPairArray*
Model::findNeighbourCandidates()
{
  //Method finds bodies which are 'neighbour bodies', ie their bounding
  // boxes are overlapping (this is a necessary but not sufficient
  // condition for common innerBoundaries between bodies).

  //Read bodies into an array.
  int b_count = modelStatistics->nofBodies;

  Body** b_array = new Body*[b_count];

  int row = 0;
  int index = 0;
  while (true) {
    Body* body = getBody(index++);
    if (body==NULL) break;
    b_array[row++] = body;
  }

  //Create a table of body-pairs which are overlapping.
  BodyPairArray* bp_array = new BodyPairArray;

  BoundBox* box1;
  BoundBox* box2;
  BodyPair* bp;

  int bp_count = 0;

  for (int i = 0; i < b_count; i++) {
    box1 = b_array[i]->getBoundBox();

    for (int j = i+1; j < b_count; j++) {
      box2 = b_array[j]->getBoundBox();

      if (box1->overlap(box2)) {
        bp = new BodyPair(b_array[i], b_array[j]);
        bp_array->push_back(bp);
        bp_count++;
      }
    }
  }

  return bp_array;
}


// Methods checks if a point is stored in the point-list *plist*
GcPoint*
Model::findPoint(PointList* plist, GcPoint* point)
{
  // Note: Point comparision must be based on objects, not pointers, because same
  // point-geometry may have been read into a new point-object from the input file.
  PointList::iterator i = plist->begin();
  while ( i != plist->end()) {

    // NOTE: GcPont class comparison operator is
    // epsilon sensitive!
    if ( **i == *point)
      return *i;
    i++;
  }

  // Point was not found.
  return NULL;
}


// Find selected mesh boundary element index using ray-casting
int
Model::findSelectedMeshBoundaryElement(Renderer* renderer, Point3& ray_start, Point3& ray_dir,
                                       bool try_current_bndr, int& bndr_id,
                                       int& body1_id, int& layer1_id,
                                       int& body2_id, int& layer2_id)
{
  return meshManager->findSelectedMeshBoundaryElement(renderer, ray_start, ray_dir,
                                                      try_current_bndr, bndr_id,
                                                      body1_id, layer1_id,
                                                      body2_id,layer2_id);
}


// Find vertex by geometry point
//
BodyElement*
Model::findVertex(GcPoint* point)
{
  BodyElement* vertex = NULL;

  // Try if point stored in the model
  GcPoint* stored_point = this->getPoint(point);

  // If point already exist in the model, check if a vertex also exists
  if ( stored_point != NULL ) {
    vertex = getVertex(stored_point);
  }

  return vertex;
}


const char*
Model::getActiveMeshName(int active_mesh_index)
{
  if ( active_mesh_index < 0 ||
       active_mesh_index >= modelInfo->nofActiveMeshes
     ) {
    return NULL;
  }

  int mindex = modelInfo->activeMeshIndices[active_mesh_index];

  if ( mindex < 0 || mindex >= modelInfo->nofMeshes ) {
    return NULL;
  }

  return modelInfo->meshNames[mindex];
}


// Method finds a body by index
Body*
Model::getBody(int index, bool only_active)
{
  return (Body*)getModelObject(index, OT_BODY, only_active);
}


Body*
Model::getBodyById(int id) const
{
  ModelObject* obj = getModelObjectById(id);

  if ( obj == NULL || obj->getObjectType() != OT_BODY)
    return NULL;
  else
    return (Body*)obj;
}


// Method checks if a body with the id-key exists.
// Returns pointer to the body or NULL if not found.
Body*
Model::getBodyByTag(int tag) const
{
  return (Body*)getModelObjectByTag(OT_BODY, tag);
}


// Method finds a body element by index.
// NOTE: order of Faces/Edges/Vertices is not specified, so they are
// delivered in what order they are encountered
//
BodyElement*
Model::getBodyElement(int index, bool only_active)
{
  int count = modelData->modelObjects->size();

  if ( index < 0 || index >= count ) {
    return NULL;
  }

  int obj_index = 0;
  int counter = 0;

  while (true) {

    ModelObject*obj = getModelObject(obj_index++);

    if ( obj == NULL ) {
      return NULL;
    }

    if ( only_active && !obj->isActive() ) {
      continue;
    }

    objectType otp = obj->getObjectType();
    if ( otp != OT_FACE && otp != OT_EDGE && otp != OT_VERTEX ) {
      continue;
    }

    if ( counter++ == index ) {
      return (BodyElement*)obj;
    }
  }

  return NULL;
}


BodyElement*
Model::getBodyElementById(int id) const
{
  ModelObject* obj = getModelObjectById(id);

  if ( obj == NULL )
    return NULL;

  objectType otp = obj->getObjectType();

  if ( otp != OT_FACE && otp != OT_EDGE && otp != OT_VERTEX)
    return NULL;

  return (BodyElement*)obj;
}


BodyElement*
Model::getBodyElementByBoundaryTag(int btag)
{
  int index = 0;
  while (true) {
    BodyElement* be = getBodyElement(index++);
    if (be==NULL) break;
    if ( btag == be->getBoundaryTag() ) {
      return be;
    }
  }

  return NULL;
}


BodyElement*
Model::getBodyElementByTag(objectType type, int tag)
{
  int index = 0;
  while (true) {
    BodyElement* be = getBodyElement(index++);
    if (be==NULL) break;
    if ( type == be->getObjectType() && tag == be->Tag() ) {
      return be;
    }
  }

  return NULL;
}


// Method finds a bodyelementloop by index.
BodyElementLoop*
Model::getBodyElementLoop(int index, bool only_active)
{
  return (BodyElementLoop*)getModelObject(index, OT_ELEMENT_LOOP, only_active);
}


BodyElementLoop*
Model::getBodyElementLoopById(int id) const
{
  ModelObject* obj = getModelObjectById(id);

  if ( obj == NULL || obj->getObjectType() != OT_ELEMENT_LOOP)
    return NULL;
  else
    return (BodyElementLoop*)obj;
}


// Method checks if a bodyElement-loop with the id-key exists.
// Returns pointer to the bodyElementLoop or NULL if not found.
BodyElementLoop*
Model::getBodyElementLoopByTag(int tag) const
{
  return (BodyElementLoop*)getModelObjectByTag(OT_ELEMENT_LOOP, tag);
}


// Find element-loop by elements in the loop.
//
// Return true if a matching loop found
// Return loop id in the reference argument 'id'
bool
Model::getBodyElementLoopId(IdList* loop_ids, int& id, int& direction)
{
  // Collect element ids into an array
  //
  int nof_ids = loop_ids->size();

  int* ids = new int[nof_ids];

  IdList::iterator pos = loop_ids->begin();

  int index = 0;
  while ( pos != loop_ids->end() ) {
    ids[index++] = *pos++;
  }

  return getBodyElementLoopId(nof_ids, ids, id, direction);

  delete[] ids;
}


// Find element-loop by elements in the loop.
//
// Return true if a matching loop found
// Return loop id in the reference argument 'id'
//
// A helper function
bool
Model::getBodyElementLoopId(int nof_ids, const int* loop_ids, int& id, int& direction)
{
  const int* ids;
  int start1, start2;

  int index = 0;

  while (true) {
    BodyElementLoop* bel = getBodyElementLoop(index++);

    if (bel==NULL) break;

    if ( bel->getNofElements() != nof_ids ) continue;

    ids = bel->getElementIds();

    // Matching loops!
    if ( idLoopsAreMatching(nof_ids, nof_ids,
                            ids, loop_ids,
                            direction, start1, start2)
       ) {
      id = bel->Id();
      return true;
    }
  }

  id = NO_INDEX;
  direction = 0;

  return false;
}


BodyLayer*
Model::getBodyLayer(int index, bool only_active)
{
  return (BodyLayer*)getModelObject(index, OT_BODY_LAYER, only_active);
}


BodyLayer*
Model::getBodyLayerById(int id)
{
  ModelObject* obj = getModelObjectById(id);

  if ( obj == NULL || obj->getObjectType() != OT_BODY_LAYER)
    return NULL;
  else
    return (BodyLayer*)obj;
}


BodyLayer*
Model::getBodyLayerByTag(int tag)
{
  int index = 0;
  while (true) {
    BodyLayer* bl = getBodyLayer(index++);
    if (bl==NULL) break;
    if ( OT_BODY_LAYER == bl->getObjectType() && tag == bl->Tag() ) {
      return bl;
    }
  }

  return NULL;
}


// Find layer for a body by body-tag
// Return layer-id or NO_INDEX if not found
//
int
Model::getBodyLayerByBodyId(int bd_id)
{
  int index = 0;

  while(true) {
    BodyLayer* bl = getBodyLayer(index++);
    if (bl==NULL) break;
    if (bd_id == bl->getBodyId() ) {
      return bl->Id();
    }
  }

  return NO_INDEX;
}


// Method finds a bodypair by index.
BodyPair*
Model::getBodyPair(int index, bool only_active)
{
  return (BodyPair*)getModelObject(index, OT_BODYPAIR, only_active);
}


BodyPair*
Model::getBodyPairById(const Body* bd1, const Body* bd2)
{
  int index = 0;
  while (true) {

    BodyPair* bp = getBodyPair(index++);

    if (bp==NULL) break;

    const Body* b1;
    const Body* b2;

    bp->getBodies(b1, b2);

    if ( b1 == bd1 && b2 == bd2 ||
         b1 == bd2 && b2 == bd1
         ) {
      return bp;
    }
  }

  return NULL;
}


int
Model::getBodyTagExt2Int(int external_id)
{
  if (external_id < 0 || external_id > MAX_NOF_BODIES ) {
    return NO_INDEX;
  }

  return meshData->bodyExt2Int[external_id];
}


int
Model::getBodyTagInt2Ext(int internal_id)
{
  if (internal_id < 0 || internal_id > MAX_NOF_BODIES ) {
    return NO_INDEX;
  }

  return meshData->bodyInt2Ext[internal_id];
}


void
Model::getBodyTags(int bd1_id, int bd2_id, int& bd1_tag, int& bd2_tag)
{
  //bd1_tag = NO_INDEX;
  //bd2_tag = NO_INDEX;

  bd1_tag = bd1_id;
  bd2_tag = bd2_id;

  if (!modelInfo->hasGeometry) {

    //bd1_tag = bd1_id;
    //bd2_tag = bd2_id;

    //return;
  }

  Body* bd1 = getBodyById(bd1_id);
  Body* bd2 = getBodyById(bd2_id);

  if ( bd1 != NULL ) {
    bd1_tag = bd1->Tag();
  }

  if ( bd2 != NULL ) {
    bd2_tag = bd2->Tag();
  }
}


// Method finds the next boundary (edge or face)
//
BodyElement*
Model::getBoundary(int index, bool only_active)
{
  if ( modelInfo->dimension == ECIF_2D )
    return (BodyElement*)getModelObject(index, OT_EDGE, only_active);
  else
    return (BodyElement*)getModelObject(index, OT_FACE, only_active);

}


// Method checks if a bodyElement with the id-key exists.
// Returns pointer to the bodyElement or NULL if not found.
BodyElement*
Model::getBoundaryById(int id) const
{
  if ( modelInfo->dimension == ECIF_2D )
    return getEdgeById(id);
  else
    return getFaceById(id);
}


// Method checks if a bodyElement with the id-key exists.
// Returns pointer to the bodyElement or NULL if not found.
BodyElement*
Model::getBoundaryByTag(int tag) const
{
  if ( modelInfo->dimension == ECIF_2D )
    return getEdgeByTag(tag);
  else
    return getFaceByTag(tag);
}


// Methods finds (the FIRSTt) body element for the body
// pair (body1,body2). This is meant mainly when there
// are unique boundaries per body/body-pair like in Ideas
// mesh files where boundaries are implied by material ids
BodyElement*
Model::getBoundaryByTags(int body1_tag, int body2_tag)
{
  int index = 0;
  while (true) {
    BodyElement* be = getBoundary(index++);
    if (be==NULL) break;
    int bd1_tag = be->getParentTag(1);
    int bd2_tag = be->getParentTag(2);

    if ( bd1_tag == body1_tag && bd2_tag == body2_tag ) {
      return be;
    }
  }

  // Not found
  return NULL;
}


void
Model::getBoundaries(int body1_id, int body2_id,
                     int& nof_boundaries,
                     BodyElement**& boundaries)
{
  nof_boundaries = 0;
  boundaries = NULL;

  if (body1_id == NO_INDEX) {
    return;
  }

  // Count boundaries
  int index = 0;
  while (true) {
    BodyElement* be = getBoundary(index++);
    if (be==NULL) break;
    int bd1_id = be->getParentId(1);
    int bd2_id = be->getParentId(2);

    if ( bd1_id == body2_id && bd2_id == body2_id ) {
      nof_boundaries++;
    }
  }

  if ( nof_boundaries == 0 ) return;

  boundaries = new BodyElement*[nof_boundaries];
  nof_boundaries = 0;

  // Store boundaries
  index = 0;
  while (true) {
    BodyElement* be = getBoundary(index++);
    if (be==NULL) break;
    int bd1_id = be->getParentId(1);
    int bd2_id = be->getParentId(2);

    if ( bd1_id == body2_id && bd2_id == body2_id ) {
      boundaries[nof_boundaries++] = be;
    }
  }

}


// Method finds the next boundary group
//
BodyElementGroup*
Model::getBodyElementGroup(int index, bool only_active)
{
  return (BodyElementGroup*)getModelObject(index, OT_ELEMENT_GROUP, only_active);
}


// Method checks if a boundary group with the id-key exists.
// Returns pointer to the boundary group or NULL if not found.
BodyElementGroup*
Model::getBodyElementGroupById(int id) const
{
  ModelObject* obj = getModelObjectById(id);

  if ( obj == NULL ) return NULL;

  objectType ot = obj->getObjectType();

  if ( OT_ELEMENT_GROUP != obj->getObjectType() ) return NULL;

  return (BodyElementGroup*)obj;
}


// Method finds a boundary group by tag
BodyElementGroup*
Model::getBodyElementGroupByTag(int tag) const
{
  return (BodyElementGroup*)getModelObjectByTag(OT_ELEMENT_GROUP, tag);
}


void
Model::getBoundingBox(double& x1, double& x2, double& y1, double& y2, double& z1, double&z2) const
{
  RangeVector rv;
  modelBox->getRangeVector(rv);

  x1 = rv[0]; x2 = rv[1];
  y1 = rv[2]; y2 = rv[3];
  z1 = rv[4]; z2 = rv[5];
}


// Get all bodies bounding box for the model;
void
Model::getBoundingBox(RangeVector rv) const
{
  modelBox->getRangeVector(rv);
}


void
Model::getCoordinateLabels(int max_len, char* label_x, char* label_y, char* label_z)
{
  static int len = 1;
  static char cartesian[] = "XYZ";
  static char cylindric[] = "RZP";
  static char polar[]     = "RTP";

  int i1 = modelInfo->coordinateMapping[0] - 1;
  int i2 = modelInfo->coordinateMapping[1] - 1;
  int i3 = modelInfo->coordinateMapping[2] - 1;

  label_x[max_len] = '\0';
  label_y[max_len] = '\0';
  label_z[max_len] = '\0';

  if (max_len < len) return;

#if 0
  // NOTE: Fixed labels!
  label_x[0] = cartesian[0];
  label_y[0] = cartesian[1];
  label_z[0] = cartesian[2];
  return;
#endif

  // Set labels according to selected
  // coordinate system
  switch (modelInfo->coordinateType) {

  case COORD_CARTESIAN:
    label_x[0] = cartesian[i1];
    label_y[0] = cartesian[i2];
    label_z[0] = cartesian[i3];
    break;

  case COORD_CYLINDRIC:
    label_x[0] = cylindric[i1];
    label_y[0] = cylindric[i2];
    label_z[0] = cylindric[i3];
    break;

  case COORD_POLAR:
    label_x[0] = polar[i1];
    label_y[0] = polar[i2];
    label_z[0] = polar[i3];
    break;
  }

  label_x[len] = '\0';
  label_y[len] = '\0';
  label_z[len] = '\0';

}


// Get active bodies bounding box for the model;
void
Model::getCurrentBoundingBox(RangeVector rv)
{
  BoundBox current_box(rv);
  RangeVector body_rv;

  int index = 0;
  while (true) {
    Body* body = getBody(index++);
    if (body==NULL) break;
    if ( !DM_HIDDEN == body->getDrawMode() ) {
      body->getRangeVector(body_rv);
      current_box.extendByRange(body_rv);
    }
  }

  current_box.getRangeVector(rv);
}


// Get active bodies bounding box for the model;
void
Model::getCurrentMeshBoundingBox(RangeVector rv)
{
  BoundBox current_box(rv);
  RangeVector body_rv;

  int index = 0;
  while (true) {
    Body* body = getBody(index++);
    if (body==NULL) break;
    if ( !DM_HIDDEN == body->getDrawMode() ) {
      body->getRangeVector(body_rv);
      current_box.extendByRange(body_rv);
    }
  }

  current_box.getRangeVector(rv);
}


const char*
Model::getCurrentMeshName()
{
  if ( modelInfo->nofMeshes == 0 || modelInfo->currentMeshIndex == NO_INDEX )
    return NULL;
  else
    return modelInfo->meshNames[modelInfo->currentMeshIndex];
}


void
Model::getCurrentTime(char* buffer)
{
  time_t t = time(NULL);
  struct tm* mt = localtime(&t);
  strftime(buffer, 64, "%d.%m.%Y  %H:%M:%S.", mt);
}


bool
Model::getColorValue(char* name, Color4 color)
{
  static char buf[1024];

  int i;
  int len = strlen(name);

  for (i = 0; i < len; i++) {
    buf[i] = tolower(name[i]);
  }
  buf[len] = '\0';

  string nm(buf);

  //Color4* c = colorTable[buf];
  Color4* c = (*colorValueTable)[nm];

  if ( c == NULL ) return false;

  for (i = 0; i < 4; i++) {
    color[i] = (*c)[i] ;
  }

  return true;
}


bool
Model::getColorName(int cn_id, char* buffer)
{
  buffer[0] = '\0';

  char* name = (*colorNameTable)[cn_id];

  if ( name == NULL ) return false;

#if 0
  if (itr != modelData->colorNames_RGB->end()) {
    char* name = (*itr).second;
    int len = strlen(name);
    strcpy(buffer, name);
    buffer[len] = '\0';
    return true;
  }
#endif

    int len = strlen(name);
    strcpy(buffer, name);
    buffer[len] = '\0';

    return true;
}


IdList*
Model::getCoveringElementList(int elem_id)
{
  return (*modelData->coveringElementTable)[elem_id];
}


bool
Model::getCoveringElementList(int elem_id, IdList& covering_elem_ids)
{
  IdList* ids = (*modelData->coveringElementTable)[elem_id];

  if ( ids == NULL ) {
    return false;
  }

  IdList::iterator pos = ids->begin();

  while (pos != ids->end()) {
    int id = *pos++;
    covering_elem_ids.push_back(id);
  }

  return true;
}


// Method finds a edge by index.
BodyElement*
Model::getEdge(int index, bool only_active)
{
  return (BodyElement*)getModelObject(index, OT_EDGE, only_active);
}


BodyElement*
Model::getEdgeById(int id) const
{
  ModelObject* obj = getModelObjectById(id);

  if ( obj == NULL || obj->getObjectType() != OT_EDGE)
    return NULL;
  else
    return (BodyElement*)obj;
}


BodyElement*
Model::getEdgeByTag(int tag) const
{
  return (BodyElement*)getModelObjectByTag(OT_EDGE, tag);
}


// Method finds a face by index.
BodyElement*
Model::getFace(int index, bool only_active)
{
  return (BodyElement*)getModelObject(index, OT_FACE, only_active);
}


BodyElement*
Model::getFaceById(int id) const
{
  ModelObject* obj = getModelObjectById(id);

  if ( obj == NULL || obj->getObjectType() != OT_FACE)
    return NULL;
  else
    return (BodyElement*)obj;
}



BodyElement*
Model::getFaceByTag(int tag) const
{
  return (BodyElement*)getModelObjectByTag(OT_FACE, tag);
}



bool
Model::getFlagValue(flagName name)
{
  bool value = modelFlags[name];

  return value;
}


const UserInterface*
Model::getGui()
{
  return theControlCenter->getGui();
}


bool
Model::getLabelDisplayFlagValue(BodyElement* be)
{
  objectType ot = be->getObjectType();

  switch (ot) {
  case OT_VERTEX:
    return modelFlags[LABEL_DISPLAY_VERTEX];
  case OT_EDGE:
    return modelFlags[LABEL_DISPLAY_EDGE];
  case OT_FACE:
    return modelFlags[LABEL_DISPLAY_FACE];
  }

  return false;
}


// Get all bodies bounding box for the model;
void
Model::getMeshBoundingBox(RangeVector rv) const
{
  meshBox->getRangeVector(rv);
}


MeshElementTable*
Model::getMeshBoundaryElementEdges()
{
  if ( modelInfo->dimension == ECIF_3D ) {
    return meshData->boundaryEdges;

  } else {
    return meshData->boundaryVertices;
  }
}


double
Model::getMeshF(int mesh_index)
{
  if ( mesh_index  < 0 || mesh_index >= modelInfo->nofMeshes )
    return -1;

  if (modelInfo->meshHs == NULL)
    return 0;

  return modelInfo->meshFs[mesh_index];
}


double
Model::getMeshH(int mesh_index)
{
  if ( mesh_index  < 0 || mesh_index >= modelInfo->nofMeshes )
    return -1;

  if (modelInfo->meshHs == NULL)
    return 0;

  return modelInfo->meshHs[mesh_index];
}


int
Model::getMeshIndex(const char* mesh_name)
{
  if ( modelInfo->nofMeshes == 0 ) {
    return NO_INDEX;
  }

  for (int i = 0; i < modelInfo->nofMeshes; i++) {

    if ( 0 == strcmp(mesh_name, modelInfo->meshNames[i]) ) {
      return i;
    }
  }

  return NO_INDEX;
}


bool
Model::getMeshInputFileName(char*& mif_file_name, int mesh_index)
{
  if ( modelInfo->nofActiveMeshes == 0 ||
       modelInfo->meshNames == NULL
     ) {
    return false;
  }

  double h = getMeshH(mesh_index);

  if ( h < 0 ) {
    return false;
  }

  UserInterface* gui = theControlCenter->getGui();

  // Ok, ask the mif-file-name
  return gui->getMeshInputFileName(mif_file_name);
}


void
Model::getMeshNames(int& nof_meshes, char**& mesh_names)
{
  nof_meshes = modelInfo->nofMeshes;
  mesh_names = NULL;

  if ( nof_meshes == 0 )
    return;

  mesh_names = new char*[nof_meshes];

  for (int i = 0; i < nof_meshes; i++) {

    int len = strlen(modelInfo->meshNames[i]);

    // Alloc space for each mesh name
    mesh_names[i] = new char[1 + len];

    strcpy(mesh_names[i], modelInfo->meshNames[i]);
  }

}


int
Model::getMeshNodeIdExt2Int(int ext_id)
{
  if (ext_id != NO_INDEX)
    return meshData->nodeExt2Int[ext_id];
  else
    return NO_INDEX;
}


int
Model::getMeshNodeIdInt2Ext(int int_id)
{
  if (int_id != NO_INDEX)
    return meshData->nodeInt2Ext[int_id];
  else
    return NO_INDEX;
}


// Range vector from the model's mesh bounding box
void
Model::getMeshRangeVector(RangeVector rv) const
{
  meshBox->getRangeVector(rv);
}


ostream&
Model::getModelStatusMessage(ostream& out) const
{
  ecif_modelStatus status = modelInfo->modelStatus;

  if (status == STATUS_OK) {
    out << "Model is ok" << endl;
    return out;
  }

  if (status & BODY_EQUATION_MISSING ) {
    out << "All bodies do not have an EQUATION definition." << endl;
  }

  else if (status & BODY_MATERIAL_MISSING ) {
    out << "All bodies do not have a MATERIAL definition." << endl;
  }

  return out;
}


// Get model object by type and index
ModelObject*
Model::getModelObject(int index) const
{
  int count = modelData->modelObjects->size();

  if (index < 0 || index >= count ) return NULL;

  int counter = -1;
  for (int i = 0; i < count && counter < index; i++) {
    ModelObject* obj = (*modelData->modelObjects)[i];
    if ( obj != NULL && ++counter == index) {
      return obj;
    }
  }

  return NULL;
}


// Get model object by index, type and activity flag
ModelObject*
Model::getModelObject(int index, objectType type, bool only_active) const
{
  int count = modelData->modelObjects->size();

  if (index < 0 || index >= count ) return NULL;

  int counter = -1;
  for (int i = 0; i < count; i++) {

    ModelObject* obj = (*modelData->modelObjects)[i];

    if ( obj != NULL && obj->getObjectType() == type ) {
      if ( (!only_active || obj->isActive()) &&
           (++counter == index)
         ) {
        return obj;
      }
    }
  }

  return NULL;
}


// Get model object by object id
ModelObject*
Model::getModelObjectById(int object_id) const
{
  int count = modelData->modelObjects->size();

  for (int i = 0; i < count; i++) {

    ModelObject* obj = (*modelData->modelObjects)[i];

    if ( obj != NULL &&
         obj->Id() == object_id
       ) {
      return obj;
    }
  }

  return NULL;
}


// Get model object by type and id tag
ModelObject*
Model::getModelObjectByTag(objectType type, int tag) const
{
  int count = modelData->modelObjects->size();

  for (int i = 0; i < count; i++) {

    ModelObject* obj = (*modelData->modelObjects)[i];

    if ( obj != NULL &&
         obj->getObjectType() == type &&
         obj->Tag() == tag
       ) {
      return obj;
    }
  }

  return NULL;
}


const char*
Model::getModelObjectNameById(int oid)
{
  ModelObject* obj = getModelObjectById(oid);

  if (obj == NULL ) {
    return NULL;
  } else {
    return obj->getName();
  }
}


// Get model object tag by object id
int
Model::getModelObjectTagById(int object_id) const
{
  ModelObject* obj = getModelObjectById(object_id);

  if ( obj == NULL )
    return NO_INDEX;

  else
    return obj->Tag();
}


int
Model::getNewObjectId()
{
  return ++lastObjectId;
}


MeshCornerElement*
Model::getMeshCornerElement(int index)
{
  if ( index < 0 ) return NULL;

  MeshCornerElementList::iterator pos = modelData->meshCornerElements->begin();
  MeshCornerElementList::iterator end_pos = modelData->meshCornerElements->end();

  for (int i = 0; i < index && pos != end_pos; i++,pos++);

  if ( pos == end_pos ) {
    return NULL;
  } else {
    return *pos;
  }
}


//Method finds the next available (1+previous largest) parameter id for
// parameter type "parameter_type"
int
Model::getNextNewParameterId(ecif_parameterType parameter_type)
{
  int largest_pid = 0;
  int index = 0;
  while (true) {
    Parameter* param = getParameter(index++, parameter_type);
    if ( param == NULL ) break;
    int pid = param->ID();
    if (largest_pid < pid) {
      largest_pid = pid;
    }
  }

  return ++largest_pid;
}


int
Model::getNofTimestepSteps()
{
  int index = 0;

  while (true) {

    Parameter* ts = getParameter(index++, ECIF_TIMESTEP);

    if (ts == NULL) break;

    // Take first active!
    if ( !ts->IsActive() ) {
      continue;
    }

    ParameterField* ts_intervals = NULL;

    if (modelHasSteadyStateProblem()) {
      ts_intervals = ts->getFieldBySifName("Steady State Max Iterations");
    } else {
      ts_intervals = ts->getFieldBySifName("Timestep Intervals");
    }

    if (ts_intervals == NULL) {
      return -1;
    }

    int ts_count = 0;

    char** data = ts_intervals->getDataStrings();
    int nof_entries = ts_intervals->getNofDataStrings();

    if (nof_entries > 0) {
      strstream strm;
      strm << data[0];

      int steps = 0;

      while (!strm.eof()) {
        strm >> steps;
        ts_count += steps;
        steps = 0; //to nullify nasty effects of eol in the stream
      }
    }

    return ts_count;
  }

  // No active found!
  return -1;
}


Parameter*
Model::getParameter(int index, ecif_parameterType parameter_type)
{
  ParameterTable* table;
  Parameter* param = NULL;

  // NOTE: This will be scalar pointer to a modelStatistics-counter,
  // do NOT delete it!
  int* use_counter;

  if ( !selectParameterTable(parameter_type, table, use_counter) ) {
    return NULL;
  }

  ParameterTable::iterator pos = table->begin();
  ParameterTable::iterator end_pos = table->end();

  param = (*pos).second;

  int count = table->size();

  for (int i = 0; i < index && pos != end_pos; i++,pos++);

  if (pos == end_pos) {
    return NULL;
  } else {
    param = (*pos).second;
    return param;
  }

}


// Method checks if a paramter with the type paramter-type and id-key pid exists.
// Returns pointer to the parameter or NULL if nothing found.
Parameter*
Model::getParameterById(ecif_parameterType parameter_type, int pid)
{
  if ( pid == NO_INDEX ) return NULL;

  ParameterTable* table;

  // NOTE: This will be a scalar pointer to a
  // modelStatistic counter, so shoud not be deleted!
  int* use_counter;

  if ( !selectParameterTable(parameter_type, table, use_counter) ) {
    return NULL;
  }

  ParameterTable::iterator pos = table->find(pid);

  if (pos != table->end()) {
    return (*pos).second;
  } else {
    return NULL;
  }
}



ParameterFieldInfo*
Model::getParameterFieldInfo(const char* parameter, const char* field)
{
  static char name_buffer[1 + 1024];
  name_buffer[0] = name_buffer[1024] = '\0';

  ParameterFieldInfo* fi = new ParameterFieldInfo;

  //--Try to get parameter field info from Gui
  if ( !theControlCenter->getUI()->getParameterFieldInfo(parameter, field, *fi) ) {
    delete fi;
    return NULL;
  }

  //--Add some more info

  //-Store field's Gui-name
  update_dyna_string(fi->guiName, (char*)field);

  //-Create field's Sif-name, if not delivered from Gui
  if ( fi->sifName == NULL || fi->sifName[0] == '\0' ) {
    fieldNameGuiToSif(field, name_buffer);
    update_dyna_string(fi->sifName, name_buffer);
  }

  return fi;
}


bool
Model::getSolverKeywordTypeGiven(const char* parameter, const char* field)
{
  if ( parameter == NULL || field == NULL ) {
    return false;
  } else {
    return theControlCenter->getUI()->getSolverKeywordTypeGiven(parameter, field);
  }
}


// Returns pointer to an existing point if the geometric-point
// is already in the model or returns null-pointer.
//
GcPoint*
Model::getPoint(GcPoint* point)
{
  int key = point->hashKey();

  PointHashTable::iterator itr = modelData->modelPoints->find(key);

  // If nothing found by the key return null.
  if (itr == modelData->modelPoints->end())
    return NULL;

  // Find the point from the key's bucket.
  else {
    PointList* plist = (*itr).second;
    return findPoint(plist, point);
  }
}


// Range vector from the model's bounding box
void
Model::getRangeVector(RangeVector rv) const
{
  modelBox->getRangeVector(rv);
}


// Method checks relative orientation of two matching edges
// NOTE: We suppouse that elements are exactly matching, but
// possibly differently oriented
// This should be used only with swapped elements!!!
int
Model::getRelativeOrientation(BodyElement* be1, BodyElement* be2)
{
  int dir12 = 1;

  // We have exact match. Mutual orientation can
  // be concluded from starting vertices only !!!###!!!
  BodyElement* v11 = be1->getFirstSubElement();
  BodyElement *v21 = be2->getFirstSubElement();

  if (v11 != v21)
    dir12 = -1;

  return dir12;
}


// Method checks if a bodyElement with the id-key exists in teh removed list.
// Returns pointer to the bodyElement or NULL if not found.
BodyElement*
Model::getRemovedBodyElement(int be_id) const
{
  BodyElementTable::iterator itr = modelData->removedModelElements->find(be_id);

  if (itr != modelData->removedModelElements->end()) {
    return (*itr).second;
  }
  else {
    return NULL;
  }
}


int
Model::getRenumberedMeshBulkElementId(int original_id)
{
  if ( meshData->bulkRenumbering == NULL ) {
    return original_id;

  } else {
    return meshData->bulkRenumbering[original_id];
  }
}


flagName
Model::getSelectionMethod()
{
  short nof_flags = 7;
  flagName flags[] = {
    SELECT_METHOD_SINGLE,
    SELECT_METHOD_ALL,
    SELECT_METHOD_BY_NEIGHBOR,
    SELECT_METHOD_BY_NORMAL,
    SELECT_METHOD_BY_PLANE,
    SELECT_METHOD_BY_BOX,
    SELECT_METHOD_BY_RECTANGLE
  };

  for (short i = 0; i < nof_flags; i++) {
    if ( modelFlags[flags[i]] )
      return flags[i];
  }

  return SELECT_METHOD_SINGLE;
}


flagName
Model::getSelectionMode()
{
  short nof_flags = 3;
  flagName flags[] = {
    SELECT_MODE_TOGGLE, SELECT_MODE_EXTEND, SELECT_MODE_REDUCE
  };

  for (short i = 0; i < nof_flags; i++) {
    if ( modelFlags[flags[i]] )
      return flags[i];
  }

  return SELECT_MODE_EXTEND;
}


int
Model::getSwapElementId(int orig_elem_id)
{
  IdNumberTable::iterator
  itr = modelData->swapElementTable->find(orig_elem_id);

  if (itr != modelData->swapElementTable->end())
    return (*itr).second;
  else
    return NO_INDEX;
}


bool
Model::getSymmetryAxis(double start[3], double end1[3], double end2[3] )
{
  //-Initialize
  for (int i=0; i < 3; i++) {
    start[i] = end1[i] = end2[i] = 0;
  }

  //-If there is no symmetry-axis defined
  Parameter* simulation = getParameter(0, ECIF_COORDINATE);

  if (simulation == NULL )
    return false;

  ParameterField* param_field = simulation->getFieldBySifName("Symmetry Plane");

  if (param_field == NULL)
    return false;

  char** symmetry_data = param_field->getDataStrings();
  char* symmetry = (char*)symmetry_data[0];

  //-Ok, find symmetry axis based on model box
  //-NOTE: A nurbs near an edge could confuse things here
  // Check how nurbs-bounding box is calulated !!!###!!!
  int crn1, crn2;
  // Bounding box corners are numbered ccw.
  // z=0: 1,2,3,4; z=1: 5,6,7,8

  if ( LibFront::ncEqual(symmetry, "None") )
    return false;
  else if ( LibFront::ncEqual(symmetry, "Y Left") ) {
    crn1 = 1; crn2 = 2;
  }
  else if ( LibFront::ncEqual(symmetry, "Y Right") ) {
    crn1 = 4; crn2 = 3;
  }
  else if ( LibFront::ncEqual(symmetry, "X Bottom") ) {
    crn1 = 1; crn2 = 4;
  }
  else if ( LibFront::ncEqual(symmetry, "X Top") ) {
    crn1 = 2; crn2 = 3;
  }

  modelBox->getBoxAxis(crn1, crn2, start, end1, end2);

  return true;
}

// Method finds a vertex by index.
BodyElement*
Model::getVertex(int index, bool only_active)
{
  return (BodyElement*)getModelObject(index, OT_VERTEX, only_active);
}


BodyElement*
Model::getVertex(GcPoint* point)
{
  Point2VertexTable::iterator itr = modelData->modelPoint2Vertices->find(point);

  if (itr != modelData->modelPoint2Vertices->end())
    return (*itr).second;
  else
    return NULL;
}


BodyElement*
Model::getVertexById(int id)
{
  ModelObject* obj = getModelObjectById(id);

  if ( obj == NULL || obj->getObjectType() != OT_VERTEX)
    return NULL;
  else
    return (BodyElement*)obj;
}


BodyElement*
Model::getVertexByNodeId(int node_id)
{
  bool only_active = false;

  int index = 0;
  while (true) {
    BodyElement* v = getVertex(index++, only_active);
    if (v==NULL) break;
    int nd_id = v->getMeshElementId(0);
    if ( nd_id == node_id ) {
      return v;
    }
  }
  return NULL;
}


BodyElement*
Model::getVertexByTag(int tag)
{
  return (BodyElement*)getModelObjectByTag(OT_VERTEX, tag);
}


// Init boundaries so that vertices get boundary ids as their
// parent-ids
void
Model::initBoundaries()
{
  int index = 0;
  while (true) {
    BodyElement* be = getBoundary(index++);
    if (be==NULL) break;
    //-Check status
    beStatus be_stat = be->getStatus();
  }
}


void
Model::initMeshData()
{
  meshBox = new BoundBox;
  meshData = new MeshData;
  meshInfo = new MeshInfo;
}


// NOTE: This is called after reading the model file!
// So, do not init any flags here which are set in the model file!
void
Model::initModelFlags()
{
  // Cad geometry is shown if it exists
  if (modelFlags[GEOMETRY_TYPE_CAD]) {
    modelFlags[DRAW_SOURCE_CAD] = true;
  }

  // Mesh geometry is also shown if it exists
  if (modelFlags[GEOMETRY_TYPE_MESH]) {
    modelFlags[DRAW_SOURCE_MESH] = true;
  }


  // Draw mode depends on geometry type and dimension
  //--2D
  if (modelInfo->dimension == ECIF_2D) {
    // If only mesh geometry
    if (!modelFlags[GEOMETRY_TYPE_CAD]) {
      modelFlags[DRAW_TARGET_EDGES] = true;
      modelFlags[DRAW_TARGET_BOUNDARIES] = true;
      modelFlags[DRAW_TARGET_BODIES] = false;
    // If at least Cad geometry
    } else {
      modelFlags[DRAW_TARGET_BODIES] = true;
      modelFlags[DRAW_TARGET_EDGES] = false;
      modelFlags[DRAW_TARGET_BOUNDARIES] = false;
    }

    if (modelFlags[GEOMETRY_TYPE_CAD]) {
      modelFlags[LABEL_DISPLAY_EDGE] = true;
    }

  //--3D
  } else if (modelInfo->dimension == ECIF_3D) {
    modelFlags[DRAW_TARGET_SURFACES] = true;
    modelFlags[DRAW_TARGET_BOUNDARIES] = true;

    if (modelFlags[GEOMETRY_TYPE_CAD]) {
      modelFlags[LABEL_DISPLAY_FACE] = true;
    }
  }

  modelFlags[SELECT_MODE_TOGGLE] = true;
  modelFlags[SELECT_METHOD_SINGLE] = true;
  modelFlags[SELECT_OBJECTS_TOGGLE] = true;
}


bool
Model::modelHasCadGeometry()
{
  if ( modelInfo->geometryType == GEOM_CAD ||
       modelInfo->geometryType == GEOM_CAD_AND_MESH
       ) {
    return true;
  } else {
    return false;
  }
}


bool
Model::modelHasDiffuseGrayRadiation()
{
  return modelInfo->hasDiffuseGrayRadiation;
}


bool
Model::modelHasEquation(const char* equation_name)
{
  int index = 0;
  while (true) {
    Parameter* eq = getParameter(index++, ECIF_EQUATION);
    if (eq==NULL) break;
    if ( eq->hasFieldValueBySifName(equation_name, "True") ) {
      return true;
    }
  }

  return false;
}


bool
Model::modelHasMeshGeometry()
{
  if ( modelInfo->geometryType == GEOM_MESH ||
       modelInfo->geometryType == GEOM_CAD_AND_MESH
       ) {
    return true;
  } else {
    return false;
  }

}


bool
Model::modelHasParameter(ecif_parameterType parameter_type, const char* param_name)
{
  int index = 0;
  while (true) {
    Parameter* p = getParameter(index++, parameter_type);
    if (p==NULL) break;
    if ( 0 == strcmp(p->getName(), param_name) ) {
      return true;
    }
  }

  return false;
}


// Check if the model is a steady-state model
//
bool
Model::modelHasSteadyStateProblem()
{
  int index = 0;
  while (true) {
    Parameter* ts = getParameter(index++, ECIF_TIMESTEP);
    if (ts==NULL) break;
    if ( ts->IsActive() && ts->hasFieldValueBySifName("Simulation Type", "Steady State") ) {
      return true;
    }
  }

  return false;
}


void
Model::initSplitCombineInfos()
{
  modelData->splitCombineInfos = new SplitCombineInfoArray(MAX_NOF_SPLIT_COMBINE_INFOS);
  for (short i = 0; i < MAX_NOF_SPLIT_COMBINE_INFOS; i++)
    (*modelData->splitCombineInfos)[i] = NULL;

  modelData->splitCombineInfoIndex = 0;
}



#if 0
void
Model::linearizeBoundaries(double delta_u, double delta_v)
{
  modelData->purgeBoundaryVertices();
  modelData->boundaryVertices = new BoundaryVertexList;

  modelStatistics->nofBoundaryVertices = 0;
  modelData->lastBoundaryVertexId = 1 + modelStatistics->nofVertices;

  // Loop all bodies in the model
  // ============================
  Body* body = getFirstBody();
  while (body != NULL ) {

    // Loop all elements in the body
    // =============================
    BodyElement* be = body->getFirstElement();
    while ( be != NULL ) {

      //---Calculate linearizing vertices
      be->calcBoundaryVertices(delta_u, delta_v);

      be = body->getNextElement();
    } // all elements in the body

    body = getNextBody();
  } // all bodies in the model

  modelStatistics->nofBoundaryVertices =  modelData->boundaryVertices->size();

  sortBoundaryVertices(GEOM_CAD);
}
#endif

#if 0
void
Model::linearizeBoundaries(double delta_u, double delta_v)
{
  modelData->purgeBoundaryVertices();
  modelData->boundaryVertices = new BoundaryVertexList;

  modelStatistics->nofBoundaryVertices = 0;
  modelData->lastBoundaryVertexId = 1 + modelStatistics->nofVertices;

  int nof_bpoints = 0;
  BoundaryVertex** bpoints;
  BodyElement** mapping_vertices;
  int* local_vertex_indices;
  int* bp_ids;

  // Table to handle current (real) vertices
  // First column is a direct key to check if the vertex
  // is already added to the list of linearized vertices
  // Second column is the new renumbered index for the vertex
  int vertex_count = modelData->modelVertices->size();
  bool* vertex_done = new bool[1 + vertex_count];

  for (int i = 0; i <= vertex_count; i++) {
    vertex_done[i] = false;
  }

  // Loop all bodies in the model
  // ============================
  Body* body = getFirstBody();
  while (body != NULL ) {

    bp_ids = NULL;
    mapping_vertices = NULL;
    local_vertex_indices = NULL;

    // Loop all elements in the body
    // =============================
    BodyElement* be = body->getFirstElement();
    while ( be != NULL ) {

      //---Calculate linearizing vertices
      be->calcBoundaryVertices(nof_bpoints,
                               bpoints,
                               mapping_vertices,
                               local_vertex_indices,
                               delta_u, delta_v);

      // Allocate table for element's linearizing vertex ids
      // NOTE: DO NOT DELETE this!!!
      // This will be given the boundary element!!!
      bp_ids = new int[nof_bpoints];
      int bp_id;

      //---Add (new) linearizing vertices to the model
      //   and store vertex ids (renumbered) in the element
      for (int i = 0; i < nof_bpoints; i++) {

        BoundaryVertex* bp = bpoints[i];

        BodyElememt* v = mapping_vertices[i];

        bool add_point;

        //-Check possible mapping vertex
        if (v != NULL) {

          int v_id = v->Tag();

          bp_id = v_id;

          // If already handled
          if ( vertex_done[v_id] ) {
            delete bp;
            add_point = false;            // do not add

          // Otherwise add
          } else {
            add_point = true;
            bp->id = bp_id;
            bp->original = true;
          }

          // Mark this vertex visited
          vertex_done[v_id] = true;
        }

        //-Check if possible closed curve/surface has marked vertices the same
        // as some previous local vertex
        else if (local_vertex_indices[i] != i) {

          int l_id = local_vertex_indices[i];

          // Pick id from already stored ids
          bp_id = bp_ids[l_id];
          bp->id = bp_id;
          add_point = false;
        }

        //-Append as a new boundary vertex
        else {
          bp_id = modelData->lastBoundaryVertexId;
          bp->id = bp_id;
          add_point = true;
          modelData->lastBoundaryVertexId++;
        }

        //-Add point to the model if it is new
        if (add_point) {
          bp->refId = bp->id;
          modelData->boundaryVertices->push_back(bp);
        }

        //-Set final id also for the element
        bp_ids[i] = bp_id;
      }

      //---Store ids in the element
      be->getTopology()->setBoundaryVertices(nof_bpoints, bp_ids);

      delete[] mapping_vertices;
      delete[] local_vertex_indices;

      be = body->getNextElement();
    } // all elements in the body

    body = getNextBody();
  } // all bodies in the model

  modelStatistics->nofBoundaryVertices =  modelData->boundaryVertices->size();

  sortBoundaryVertices(GEOM_CAD);

  delete[] vertex_done;
}
#endif


// Check if a vertex is in model's vertex-table
//
bool
Model::isInVertexTable(int vertex_id)
{
  VertexTable* vt = modelData->vertexTable;

  if ( vt == NULL ) return false;

  for (int i = 0; i < vt->dim1; i++) {
    if ( vt->vertexIds[i] == vertex_id ) {
      return true;
    }
  }

  return false;
}


bool
Model::loadDBMesh(char* disp_msg)
{
  if ( modelInfo->meshDirectory_absolute[0] == '\0' ||
       modelInfo->meshNames == NULL ||
       modelInfo->currentMeshIndex == NO_INDEX
     ) {
    return false;
  }

  Timer timer;
  timer.start();

  UserInterface* gui = theControlCenter->getGui();

  if (disp_msg != NULL) {
    gui->showMsg(disp_msg);
  }

  // NOTE: infile-stream (dummy_strm) is not really used by InputElmer, cause
  // eio-library does not use directly the given meshfile name, but the parent
  // directory
  // NOTE: It is important anyway here to add the filename itself to the path
  // (mesh.header) because it is suppoused to be there in InputElmer!
  // ( this in turn is because complete filename is delivered from Gui via
  //   Control's readMeshFile method to InputElmer!)
  strstream strm;
  strm << modelInfo->meshDirectory_absolute << "/" << "mesh.header" << ends;

  ifstream dummy_strm;
  InputElmer in(modelInfo->dimension, dummy_strm, strm.str());

  if ( ECIF_ND == in.loadMesh() ) {
    gui->showMsg("NOTE: No mesh file found for the model!");
    gui->configureMenuButtons("File", "LoadMesh", 0);
    gui->configureMenuButtons("File", "Mesh", 0);
    resetMeshData();
    return false;
  }

  //gui->showUsedTimeMsg(timer.stop(), "Reading the DB mesh file");

  if ( !processMeshFileData(&in) )
    return false;

  addPendingMeshElements();

  return true;
}


bool
Model::loadMesh()
{
  bool has_mesh_geometry = false;

  //---Read mesh geometry
  // NOTE: missing mesh geometry is not a fatal error,
  // so we return true although it is missing
  char* file1 = modelInfo->meshResultFile;
  char* file2 = modelInfo->meshSourceFile;

  Input* input = NULL;
  ifstream* in_file = NULL;

  FILE* mesh_file = NULL;

  UserInterface* ui = theControlCenter->getGui();

  // Try first to read Elmer DB mesh file
  // ====================================
  if ( !has_mesh_geometry ) {

    if ( loadDBMesh() ) {
      has_mesh_geometry = true;

    } else {
      ui->showMsg("***WARNING Could not read Elmer mesh!");
    }
  }

  // Try next to read mesh result file if it exists
  // ==============================================
  if ( !has_mesh_geometry && file1[0] != '\0' ) {

    ui->showMsg("***NOTE Trying to read mesh result file");

    // *** Check that input-file exists.
    mesh_file = fopen(file1, "r");

    if ( !mesh_file ) {
      fclose(mesh_file);
      ui->errMsg(0, "Error! Can't open mesh result file: ", file1);
    }
    else {

      //in_file = new ifstream(file1, ios::nocreate);
      in_file = new ifstream(file1);
      input = theControlCenter->create_mesh_input(modelInfo->dimension,
                                                  *in_file, file1);

      if (input == NULL) {
        fclose(mesh_file);
        ui->errMsg(0, "Error! Can't open mesh result file: ", file1);
      }

      else if ( ECIF_ND != input->readMeshFile() ) {
          has_mesh_geometry = true;
      }
    }
  }

  // Try finally to read the original mesh file
  // ==========================================
  if ( !has_mesh_geometry && file2[0] != '\0' ) {

    ui->showMsg("***NOTE Trying to read mesh source file");

    // *** Check that input-file exists.
    mesh_file = fopen(file2, "r");

    if ( mesh_file == NULL) {
      ui->errMsg(0, "Error! Can't find mesh source file: ", file2);

    } else {

      //in_file = new ifstream(file2, ios::nocreate);
      in_file = new ifstream(file2);

      input = theControlCenter->create_mesh_input(modelInfo->dimension,
                                                  *in_file, file2);

      if (input == NULL) {
        fclose(mesh_file);
        ui->errMsg(0, "Error! Can't open mesh source file: ", file2);
        has_mesh_geometry = false;

      } else if ( ECIF_ND != input->readMeshFile() ) {
          has_mesh_geometry = true;
      }
    }

  }

  if (mesh_file != NULL)
    fclose(mesh_file);

  // Process mesh file data
  // ======================
  if (has_mesh_geometry && input != NULL) {
    if ( !processMeshFileData(input) ) {
      has_mesh_geometry = false;
    }
  }

  // If success
  // ==========
  if (has_mesh_geometry) {

    UserInterface* gui = theControlCenter->getGui();

    setFlagValue(GEOMETRY_TYPE, GEOMETRY_TYPE_MESH, true);
    setFlagValue(DRAW_SOURCE, DRAW_SOURCE_MESH, true);
    gui->updateModelFlags(this);

    modelInfo->meshNeedsUpdate = false;
    gui->setWasUpdated("Mesh");

    gui->configureButtons("draw_source_mesh", 1);
    gui->configureMenuButtons("File", "LoadMesh", -1);

    setWindowTitles();
  }


  delete input;
  delete in_file;

  return has_mesh_geometry;
}


// Set active flags in objects
// NOTE: This must be done hierachically: bodies -> loops -> body elements
//
void
Model::markActiveObjects()
{
  int i, index;

  // Mark first all object inactive
  // ==============================
  int count = modelData->modelObjects->size();

  for (i = 0; i < count; i++) {

    ModelObject* obj = getModelObject(i);

    if ( obj != NULL ) {
      obj->setActive(false);
    }
  }

  // Mark all body layers active!
  // ============================
  index = 0;
  while (true) {
    BodyLayer* bl = getBodyLayer(index++, false);
    if (bl==NULL) break;
    bl->setActive(true);
  }

  // Mark first all bodies active!
  // ==============================
  index = 0;
  while (true) {
    Body* body = getBody(index++, false);
    if (body==NULL) break;
    body->setActive(true);
    body->markActiveObjects();
  }

  // Then bodypairs from bodies active!
  // ==================================
  index = 0;
  while (true) {
    BodyPair* bp = getBodyPair(index++, false);
    if (bp==NULL) break;
    const Body* body1;
    const Body* body2;
    bp->getBodies(body1, body2);
    Body* bd1 = (Body*) body1;
    Body* bd2 = (Body*) body2;

    if ( bd1 != NULL && bd1->isActive() &&
         bd2 != NULL && bd2->isActive()
       ) {
      bp->setActive(true);
    }
  }

  // Mark all boundary groups active!
  // ================================
  index = 0;
  while (true) {
    BodyElementGroup* beg = getBodyElementGroup(index++, false);
    if (beg==NULL) break;
    beg->setActive(true);
 }

  // Then read only active (sub) objects!
  // ====================================

  // Body element loop
  // =================
  index = 0;
  while (true) {
    BodyElementLoop* bel = getBodyElementLoop(index++);
    if (bel==NULL) break;
    bel->markActiveObjects();
  }

  // Body elements
  // =============
  index = 0;
  while (true) {
    BodyElement* be = getBoundary(index++);
    if (be==NULL) break;
    be->markActiveObjects();
 }

}


void
Model::updateModelStatistics()
{
  int index;

  // Bodies
  // ======
  modelStatistics->nofBodies = 0;
  modelStatistics->nofElementLoops = 0;
  modelStatistics->nofElements = 0;
  modelStatistics->nofInnerBoundaries = 0;
  modelStatistics->nofOuterBoundaries = 0;
  modelStatistics->nofVertices = 0;

  modelStatistics->maxLoopSize = 0;

  index = 0;
  while (true) {
    Body* body = getBody(index++);
    if (body==NULL) break;
    if ( body->isActive() ) {
      modelStatistics->nofBodies++;
    }
  }

  // Body element loops
  // ==================
  index = 0;
  while (true) {
    BodyElementLoop* bel = getBodyElementLoop(index++);
    if (bel==NULL) break;
    if ( bel->isActive() ) {
      modelStatistics->nofElementLoops++;

      int elem_c = bel->getNofElements();
      if ( elem_c > modelStatistics->maxLoopSize ) {
        modelStatistics->maxLoopSize = elem_c;
      }
    }
  }

  // Body elements
  // =============
  index = 0;
  while (true) {
    BodyElement* be = getBoundary(index++);
    if (be==NULL) break;
    if ( be->isActive() ) {
      modelStatistics->nofElements++;

      if ( NO_INDEX == be->getParentId(2) )
        modelStatistics->nofOuterBoundaries++;
      else
        modelStatistics->nofInnerBoundaries++;
    }
  }

  // Vertices
  // ========
  index = 0;
  while (true) {
    BodyElement* v = getVertex(index++);
    if (v==NULL) break;
    if ( v->isActive() ) {
      modelStatistics->nofVertices++;
    }
  }

}


bool
Model::markObjectActive(int id)
{
  ModelObject* obj = getModelObjectById(id);

  if (obj == NULL ) {
//theControlCenter->getGui()->showMsg("Model: not marking NULL object!", 1);
    return false;
  }

  obj->setActive(true);

  return true;
}


// Mesh boundary element selected by  renderer selection hit method
bool
Model::meshBoundaryElementSelectionHit(Renderer* renderer, int fem_id)
{
  if ( fem_id < 0 ||
       fem_id >= meshData->boundaryElements->NofElements()
     ) {
    return false;
  }

  meshInfo->selectedBndrElementId = fem_id;

  selectMeshBoundaryElement(fem_id);

  return true;
}


bool
Model::meshBoundaryElementSelected(Renderer* renderer, int fem_id)
{
  static int prev_fem_id = NO_INDEX;
  static flagName prev_select_mode = SELECT_MODE_EXTEND;

  flagName select_mode = getSelectionMode();

  int mk_state = renderer->getMouseKeyboardState();

  bool paint_action = ( 0 != (mk_state & MK_SHIFT & (MK_LBUTTON | MK_RBUTTON) ) );

  // When painting we want to avoid all
  // unnecesary rendering
  // If same element was selected and
  // selection mode remained the same
  // (extend or reduce), skip redisplay
  if ( paint_action &&
       prev_fem_id == fem_id &&
       select_mode != SELECT_MODE_TOGGLE &&
       prev_select_mode == select_mode
     ) {
    return false;
  }

  prev_fem_id = fem_id;
  prev_select_mode = select_mode;

  meshInfo->selectedBndrElementId = fem_id;

  // Left/right buttons do the opposite action!
  if (mk_state & MK_RBUTTON)
    unselectMeshBoundaryElement(fem_id);
  else
    selectMeshBoundaryElement(fem_id);

  return true;
}


// Mesh bulk element selected by  renderer selection hit method
bool
Model::meshBulkElementSelectionHit(Renderer* renderer, int fem_id)
{
  if ( fem_id < 0 ||
       fem_id >= meshData->bulkElements->NofElements()
     )
    return false;

  meshInfo->selectedBulkElementId = fem_id;


  if ( getFlagValue(SELECT_MODE_TOGGLE) )
    meshData->bulkElements->selected[fem_id] = !meshData->bulkElements->selected[fem_id];
  else if ( getFlagValue(SELECT_MODE_EXTEND) )
    meshData->bulkElements->selected[fem_id] = true;
  else if ( getFlagValue(SELECT_MODE_REDUCE) )
    meshData->bulkElements->selected[fem_id] = false;

  return true;
}


#if 0
void
Model::normalizeMeshPoints()
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

void
Model::normalizeVertices()
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

  int index = 0;
  while (true) {
    BodyElement* v = getVertex(index++);
    if (v==NULL) break;
    GcPoint* vp = (GcPoint*)v->getGeometry();
    vp->normalize(norm, shift);
  }
}


bool
Model::processCadFileData()
{
  UserInterface* gui = theControlCenter->getGui();

  strstream strm;
  strm << "Nof bodies found: " << modelStatistics->nofBodies << ends;
  gui->showMsg(strm.str());

  if ( !checkBodyElementLoops() ) {
    return false;
  }

  //--Update bodies' boundboxes
  int index = 0;
  while (true) {
    Body* body = getBody(index++);
    if (body==NULL) break;
    body->updateBoundBox();
  }

  if ( !checkBodies() ) {
    return false;
  }

  //--Ok, now we have decent geometry
  modelInfo->hasGeometry = true;

  findBoundaries();
  initBoundaries();
  setBoundaryParentIdsAndTags();
  findBodyPairs();
  findContainedBodies();
  //normalizeVertices();

  //if ( !checkBodyElementLoops() ) {
  //  return false;
  //}

  // Check boundary group data
  if ( !checkElementGroupData() ) {
    return false;
  }

  markActiveObjects();

  setBoundaryTags();

  updateModelStatistics();

  setFlagValue(GEOMETRY_TYPE, GEOMETRY_TYPE_CAD, true);

  // Set default mesh-h
  setInitialMeshH();

  return true;
}


// NOTE: The order of processing is important here!
bool
Model::processMeshFileData(Input* input)
{
  if ( !input->processMeshFileData() ) {
    return false;
  }

  if (!modelInfo->hasGeometry) {
    setBoundaryParentIdsAndTags();
    findBodyPairs();
    removeEmptyBoundaries();
  }

  meshManager->findMeshElementNodeParents(meshData->boundaryElements,
                                          meshData->nofAllocNodes,
                                          meshData->nofBoundaryNodeParents,
                                          meshData->boundaryNodeParentIds);

  meshManager->findMeshBoundaryBorders();

  meshManager->findMeshCornerElements();

  theControlCenter->getGui()->showMsg("Reading mesh done.", 1);

  meshData->boundaryElements->calcNormalsAndDistances();

  meshManager->calcMeshBoundaryNodeNormals();

  if ( !checkMeshBodies() ) {
    return false;
  }

  // If we have only mesh source
  if ( !modelInfo->hasGeometry ) {

    checkBodies();

    // Create subelements
    //--3D: Create Edges and Vertices
    if ( modelInfo->dimension == ECIF_3D ) {
      meshManager->createMeshSubElements(3, 2);
      //meshManager->createMeshSubElements(2, 3);

    //--2D: Create Vertices
    } else {
      meshManager->createMeshSubElements(2, 2);
    }

    checkVertexExistence();

    //--Ok, now we have decent geometry
    modelInfo->hasGeometry = true;
    
    // This should not be applied here, it would set invalid values for
    // boundaries' 'isIntraLayer' flags!!! MVe 02.06.2004
    //
    //setBoundaryParentIdsAndTags();

    // Check boundary group data
    if ( !checkElementGroupData() ) {
      return false;
    }

    markActiveObjects();

    setBoundaryTags();

    updateModelStatistics();
  }

  RangeVector rv;
  meshBox->getRangeVector(rv);
  modelBox->extendByRange(rv);

  setInitialMeshH();

  // Dircet call for test, remove this!
  //classifyMeshCornerElements();
  //correctMeshZeroVelocityElements();

  return true;
}


bool
Model::processModelFileData()
{
  UserInterface* ui = theControlCenter->getGui();

  //--Update bodies' boundboxes
  int index = 0;
  while (true) {
    Body* body = getBody(index++);
    if (body==NULL) break;
    body->updateBoundBox();
  }

  //--Check that bodies' geometry makes sense
  if (!checkBodies()) {
    return false;
  }

  if (!checkBodyElements()) {
    return false;
  }

  //--Ok, now we have decent geometry
  modelInfo->hasGeometry = true;

  setBoundaryParentIdsAndTags();
  findBodyPairs();

  // Check boundary group data
  if ( !checkElementGroupData() ) {
    return false;
  }

  markActiveObjects();

  // NOTE: Bondary tags are now read from the model file!
  // They must the same as in the elmer mesh boundary file, so
  // we cannot renumber them afterwards here (or anyway not without
  // rewriting the mesh!)
  //setBoundaryTags();

  updateModelStatistics();

  if ( modelInfo->frontInputVersionNbr <= 4 ) {
    updateParametersParentTags();
    setMeshFs();
    setMeshHs();
    setCurrentMeshIndex(0);
  }

  updateParametersParentId();

  updateParametersApplyCounts(ECIF_BODY_FORCE);
  updateParametersApplyCounts(ECIF_BOUNDARY_CONDITION);
  updateParametersApplyCounts(ECIF_EQUATION);
  updateParametersApplyCounts(ECIF_INITIAL_CONDITION);
  updateParametersApplyCounts(ECIF_MATERIAL);

  setInitialMeshH();

  //--Check if model has DiffuseGray radiation boundary condition
  //  This info is needed when checking if Viewfactors filename
  //  should be output etc.
  checkDiffuseGrayRadiation();

  modelInfo->meshNeedsUpdate = true;
  ui->setNeedsUpdate("Mesh");

  return true;
}


// Delete all un-marked parameters (they are considered as inactive!)
// This is done after updating all active parameters
void
Model::processParametersAfterUpdate(ecif_parameterType parameter_type)
{
#if 0
  // All un-marked parameters will be deleted as
  // inactive
  //
  Parameter* param = getFirstParameter(parameter_type);

  while (param != NULL) {

    if ( !param->getUpdateFlag() ) {
      deleteParameter(parameter_type, param->ID());
    }

    param = getNextParameter(parameter_type);
  }
#endif

}


// Mark all parameters as inactive
// This is used before updating all active parameters, so
// that all inactive can be easily deleted by the update-flag
void
Model::processParametersBeforeUpdate(ecif_parameterType parameter_type)
{
  deleteParameters(parameter_type);

#if 0
  // Mark all parameter as non-touched, so that
  // those which are non-touched after mass update,
  // can be deleted as inactive parametrs
  //
  Parameter* param = getFirstParameter(parameter_type);
  while (param != NULL) {
    param->setUpdateFlag(false);
    param->setApplyCount(0);
    param = getNextParameter(parameter_type);
  }
#endif

#if 0
// Old version of: setParameter(...)

  ParameterTable* table;
  ParameterTable::iterator* table_position;
  int* use_counter;

  if ( !selectParameterTable(parameter_type, table, table_position, use_counter) )
    return;

  *use_counter = 0;

  // Reset objects' parameter counters
  switch (parameter_type) {
  case ECIF_EQUATION:
    modelStatistics->nofBodiesWithEquation = 0;
    break;
  case ECIF_MATERIAL:
    modelStatistics->nofBodiesWithMaterial = 0;
    break;
  case ECIF_BODY_FORCE:
    modelStatistics->nofBodiesWithBodyForce = 0;
    break;
  case ECIF_INITIAL_CONDITION:
    modelStatistics->nofBodiesWithInitialCondition  = 0;
    break;
  case ECIF_BOUNDARY_CONDITION:
    modelStatistics->nofInnerBoundariesWithCondition = 0;
    modelStatistics->nofOuterBoundariesWithCondition = 0;
    break;
  }
#endif

}

// Read RGB colors into model colorTable from a file
//
int
Model::readColorFile(char* in_file)
{
  UserInterface* gui = theControlCenter->getGui();

  char buffer[256];
  Color4 color;
  color[3] = MAX_NOF_COLOR_LEVELS;
  int parts;
  char name_parts[3][32];
  //
  ifstream infile(in_file);
  infile.flags(ios::skipws);

  int counter = 0;
  while (!infile.eof()) {

    // Read one line into buffer
    infile.getline(buffer, 256);

    if ( strlen(buffer) == 0 ) break;

    // Skip comment rows
    if (buffer[0] == '!') continue;

    counter++;

    strstream strm;
    strm << buffer << ends;

    // Red, Green, Blue values (0-255)
    strm >> color[0] >> color[1] >> color[2];

    // Read name part(s)
    parts = -1;
    while (!strm.eof()) {
      strm >> name_parts[++parts];
    }
    // Combine all parts into the name string
    int name_len = 0;
    for (int i= 0; i <= parts; i++) {
      int len = strlen(name_parts[i]);
      strcpy(&buffer[name_len], name_parts[i]);
      name_len += len;
      if (i < parts)
        buffer[name_len++] = ' ';
    }

    char* name = NULL;
    update_dyna_string(name, buffer);

    // Store name in colorName2ValueTable (name --> value table)
    addColorName(name, color);

    // Store name in colorValue2Name table (value --> name table)
    int color_value = rgbColor2Id(color);
    addColorValue(color_value, name);
  }

  // Inform Gui about the file
  gui->colorFileWasRead(in_file);

  return counter;
}


// Read Matc definitions and load them into Matc
//
int
Model::readMatcFile(char* in_file, char* info, bool must_exists)
{
  UserInterface* gui = theControlCenter->getGui();

  char* source = NULL;
  char matc_buffer[10001];
  char name_buffer[256];
  bool is_var, is_func;

  int rc = LibFront::readFile(in_file, source);

  // If file not found or it is empty
  //
  if ( rc <= 0 ) {
    strstream strm;
    if ( rc == -1 ) {
      strm << "***WARNING: Matc file " << in_file << " not found!" << ends;
    } else {
      strm << "***WARNING: Matc file " << in_file << " is empty!" << ends;
    }

    gui->showMsg(strm.str());
    delete[] source;
    return 0;
  }

  strstream strm;
  strm << "format(8,\"rowform\")" << ends;
  mtc_domath( strm.str() );

  int source_pos = 0;
  int source_len = strlen(source);

  // Read source buffer till end
  //
  while (source_pos < source_len ) {

    // Read expression and set new source-pos just after the expression
    //
    source_pos = LibFront::readMatcExpression(source, source_pos, source_len,
                                              matc_buffer, 10000);

    LibFront::trim(matc_buffer);

    if ( matc_buffer[0] == '\0' ) continue;

    LibFront::getMatcName(matc_buffer, name_buffer, 255, is_var, is_func);

    // Evaluate the Matc-expression
    //
    char* matc_result = mtc_domath(matc_buffer);

    // An empty result
    if ( matc_result == NULL || matc_result[0] == '\0' ) continue;

    // Matc Error
    if ( LibFront::isMatcError(matc_result) ) {
      strstream strm, strm2;
      strm << "***WARNING: Matc evaluation error (" << matc_buffer << ")!" << ends;
      gui->showMsg(strm.str());
      strm2 << matc_result<< ends;
      gui->showMsg(strm2.str());
      continue;
    }
  }

  delete[] source;

  // Inform Gui about the file
  gui->matcFileWasRead(in_file);

  return 1;
}


// Makes an integer valued color-id from rgb vector
// rgb values: 0 - MAX_NOF_COLOR_LEVES
int
Model::rgbColor2Id(Color4& color)
{
  // create combined code
  int code = 0;
  code |= (color[0] << 16); // red
  code |= (color[1] << 8);  // green
  code |= color[2];       // blue
  return code;
}


// Makes an integer id from hex-string color value
int
Model::rgbColor2Id(int len, char* hex_value)
{
  // create combined code from hex string
  int code = 0;
  for (int i = 0; i < len; i++) {
    char c = hex_value[i];
    if ( c >= '0' && c <= '9')
      code |= (hex_value[i] - '0') << (4 * (len - i - 1));
    else
      code |= (10 + hex_value[i] - 'a') << (4 * (len - i - 1));
  }
  return code;
}


// Makes a hex string from integer valued color value
void
Model::rgbColorId2Hex(int color_id, char* hex_value, int len)
{
  int mask = 0x000000f << 4*(len - 1);
  for (int i = 0; i < len; i++) {
    int nbr = (color_id & mask) >> 4*(len - i - 1);
    if (nbr > 9)
      hex_value[i] = char('a' + nbr - 10);
    else
      hex_value[i] = char('0' + nbr);
    mask >>= 4;
  }
  hex_value[len] = '\0';
}


// Makes a color vector from hex string color
void
Model::rgbHex2Color(int len, char* hex_value, Color4& color)
{
  int nof_hex_nbrs = len / 3;
  int mask;
  int grp_len = nof_hex_nbrs * 4;
  switch (nof_hex_nbrs) {
  case 1: mask = 0x0000000f << 1 * 8; break;
  case 2: mask = 0x000000ff << 2 * 8; break;
  }
  int color_id = rgbColor2Id(len, hex_value);
  for (int i = 0; i < 3; i++) {
    int color_level = (color_id & mask) >> (3 - i - 1) * grp_len;
    color[i] = color_level;
    mask >>= grp_len;
  }
}


void
Model::refreshRenderer()
{
  if ( modelInfo->dimension == ECIF_ND ) {
    return;
  }

  Renderer* renderer = theControlCenter->getRenderer();

  if ( renderer != NULL && renderer->isVisible() )
    renderer->refresh();
}


// Methods removes a body from the model.
// What about element loops in the body??? Mve 12.11.02
int
Model::removeBody(Body* body)
{
  int index = 0;
  while ( true ) {

    BodyElement*be = body->getElement(index++);

    if (be==NULL) break;

    if ( be->isInnerBoundary() ) {

      beStatus be_st = be->getStatus();
      be->setStatus( (be_st & ~BE_INNER) | BE_OUTER );

      int bd1_id = be->getParentId(1);
      int bd2_id = be->getParentId(2);
      int bd1_tag = be->getParentTag(1);
      int bd2_tag = be->getParentTag(2);
      int bd1_lr = be->getParentLayer(1);
      int bd2_lr = be->getParentLayer(2);

      Body* body2;

      if ( bd1_id == body->Id() ) {
        body2 = getBodyById(bd2_id);
        be->setParentIds(bd2_id, NO_INDEX);
        be->setParentTags(bd2_tag, NO_INDEX);
        be->setParentLayers(bd2_lr, NO_INDEX);

      } else {
        body2 = getBodyById(bd1_id);
        be->setParentIds(bd1_id, NO_INDEX);
        be->setParentTags(bd1_tag, NO_INDEX);
        be->setParentLayers(bd1_lr, NO_INDEX);
      }

      removeBodyPair(body, body2);

      modelStatistics->nofInnerBoundaries--;
      modelStatistics->nofOuterBoundaries++;

    } else {
      removeBodyElement(be, true);
    }

  }

  removeModelObject(body->Id());

  return --modelStatistics->nofBodies;
}


// Methods removes a bodyelement from the model.
int
Model::removeBodyElement(BodyElement* be, bool remove_subs)
{
  be->checkLastTag(true);
  be->checkLastBoundaryTag(true);

  // Remove all sub elements
  if ( remove_subs ) {
    int nof_subs = be->getNofSubElements();
    for (int i = 0; i < nof_subs; i++) {

      BodyElement* se = be->getSubElement(i);

      if ( se != NULL ) {
        removeBodyElement(se, true);
      }
    }
  }

  // Update also point2Vertex info
  if ( OT_VERTEX == be->getObjectType() ) {
    removeVertex(be);
  } else {
    modelStatistics->nofElements--;
  }

  // Remove element from the group
  BodyElementGroup* beg = (BodyElementGroup*)be->getElementGroup();
  if ( beg != NULL ) {
    beg->removeElement(be->Id());
  }

  // Remove from active elements
  removeModelObject(be->Id());

  // Add to removed elements
  (*modelData->removedModelElements)[be->Id()] = be;

  return modelStatistics->nofElements;
}


// Methods removes a bodyelement from the model.
int
Model::removeBodyElement(BodyElement* be, bool remove_from_bodies, bool remove_subs)
{
  if (remove_from_bodies) {

    Body* body1 = getBodyById(be->getParentId(1));
    Body* body2 = getBodyById(be->getParentId(2));

    return removeBodyElement(be, body1, body2, remove_subs);
  }

  else {
    return removeBodyElement(be, remove_subs);
  }
}


// Methods removes a bodyelement from the model.
int
Model::removeBodyElement(BodyElement* be, Body* body1, Body* body2, bool remove_subs)
{
  // Remove element from bodies and remove also
  // the corresponding inner/outer boundary entry

  //-Orphant boundary!!!
  if (body1 == NULL) {
    this->removeBodyElement(be, remove_subs);

  //-Outer boundary
  } else if (body2 == NULL) {
    body1->removeElement(be->getParentLayer(1), be);

  //-Inner boundary
  } else {
    body1->removeElement(be->getParentLayer(1), be);
    body2->removeElement(be->getParentLayer(2), be);
  }

  // Remove the element itself
  return removeBodyElement(be, remove_subs);
}


// Methods removes a bodyelement group from the model.
int
Model::removeBodyElementGroup(int beg_id)
{
  BodyElementGroup* beg = getBodyElementGroupById(beg_id);

  if ( beg == NULL )
    return modelStatistics->nofElementGroups;

  removeModelObject(beg->Id());

  return --modelStatistics->nofElementGroups;
}


// Methods removes a bodyelement loop from the model.
int
Model::removeBodyElementLoop(int bel_id)
{
  BodyElementLoop* bel = getBodyElementLoopById(bel_id);

  if ( bel == NULL )
    return modelStatistics->nofElementLoops;

  removeModelObject(bel->Id());

  return --modelStatistics->nofElementLoops;
}


// Methods removes a bodylayer from the model.
int
Model::removeBodyLayer(BodyLayer* lr)
{
  if ( lr == NULL )
    return 0;

  removeModelObject(lr->Id());

  return 1;
}


int
Model::removeBodyPair(const Body* bd1, const Body* bd2)
{
  BodyPair* bp = NULL;

  int index = 0;
  while (true) {
    bp = getBodyPair(index++);
    if (bp==NULL) break;
    const Body* b1;
    const Body* b2;

    bp->getBodies(b1, b2);

    if ( b1 == bd1 && b2 == bd2 ||
         b1 == bd2 && b2 == bd1
       )
      break;
  }

  if ( bp != NULL ) {
    removeModelObject(bp->Id());
    modelStatistics->nofBodyPairs--;
  }

  return modelStatistics->nofBodyPairs;
}


void
Model::removeCadGeometry()
{
  deleteCadData();

  int index = 0;
  while (true) {
    Body* body = getBody(index++);
    if (body==NULL) break;
    body->setGmtrType(MESH_BODY);
  }
}


void
Model::removeEmptyBoundaries()
{
  BodyElement** bels = new BodyElement*[modelStatistics->nofElements];

  int counter = 0;

  int index = 0;
  while (true) {
    BodyElement* be = getBoundary(index++);
    if (be==NULL) break;
    if ( 0 == be->getNofMeshElements() ) {
      bels[counter++] = be;
    }
  }

  for (int i = 0; i < counter; i++) {
    removeBodyElement(bels[i], true);
  }
}



void
Model::removeMeshGeometry()
{
  resetMeshData();
}


int
Model::removeModelObject(int obj_id)
{
  (*(modelData->modelObjects))[obj_id] = NULL;

  return --modelStatistics->nofModelObjects;
}


// Method reomves a vertex from the model.
// Point2Vertex-table is also updated
int
Model::removeVertex(BodyElement* vertex)
{
  removeModelObject(vertex->Id());

  modelData->modelPoint2Vertices->erase( (GcPoint*)vertex->getGeometry() );

  return --modelStatistics->nofVertices;
}


// Methods stores a bodyelement back into the model.
int
Model::restoreBodyElement(BodyElement* be, bool add_to_bodies)
{
  int nof_model_objects = restoreModelObject(be);

  addBodyElement(be, true);

  BodyElementGroup* beg = (BodyElementGroup*)be->getElementGroup();

  if ( beg != NULL ) {
    beg->addElement(be->Id());
  } else {
    be->checkElementGroupData();
  }

  return nof_model_objects;
}


void
Model::restoreBoundaryNames()
{
  int index = 0;
  while(true) {
    BodyElement* be = getBoundary(index++);
    if (be==NULL) break;
    be->restoreName();
  }
}


// Methods stores a model object back into the model.
int
Model::restoreModelObject(ModelObject* obj)
{
  (*(modelData->modelObjects))[obj->Id()] = obj;

  return ++modelStatistics->nofModelObjects;
}


void
Model::resetAllBoundarySelections(bool update_gui)
{
  UserInterface* gui = theControlCenter->getUI();

  // We loop all boundaries, edges and vertices!
  int index = 0;
  while (true) {

    BodyElement* be = getBodyElement(index++);
    if (be==NULL) break;

    be->setDrawState(DS_NORMAL);

    // Set element's selection state in the Gui
    if (update_gui) {
      gui->setBoundarySelectionMode(be->Id(), false);
    }
  }

  // If we have mesh geometry, deselect all
  // mesh elements
  if ( modelFlags[GEOMETRY_TYPE_MESH] ) {
    meshData->boundaryElements->resetSelected();
  }

  refreshRenderer();
}


void
Model::resetBoundarySelections(bool update_gui, bool use_boundary_groups,
                               int nof_skip_ids, const int* skip_bndr_ids,
                               bool call_update)
{
  // Don't remove any current selections if we are in
  // extended object selection mode
  if ( modelFlags[SELECT_OBJECTS_EXTEND] ) return;

  UserInterface* gui = theControlCenter->getUI();

  int index = 0;
  while (true) {

    BodyElement* be  = getBodyElement(index++);

    if ( be == NULL ) break;

    bool skip = false;

    for (int i = 0; i < nof_skip_ids; i++) {
      if ( be->Id() == skip_bndr_ids[i] ) {
        skip = true;
        break;
      }
    }

    if ( skip ) continue;

    be->setDrawState(DS_NORMAL);

    if (update_gui) {

      if (false && use_boundary_groups) {
        gui->setBoundarySelectionMode(be->getElementGroupId(), false, call_update);
      } else {
        gui->setBoundarySelectionMode(be->Id(), false, call_update);
      }
    }

  }
}


// **** Parameter id attachment reset methods.

void
Model::resetBoundaryConditions()
{
  int index = 0;
  while (true) {
    BodyElementGroup* beg = getBodyElementGroup(index++);
    if (beg==NULL) break;
    beg->setBoundaryConditionId(NO_INDEX);
  }
}


void
Model::resetInitialConditions()
{
  int index = 0;
  while (true) {
    Body* bd = getBody(index++);
    if (bd==NULL) break;
    bd->setInitialConditionId(NO_INDEX);
  }
}


void
Model::resetMeshData()
{
  deleteMeshData();
  initMeshData();

  meshManager->setMeshData(meshData, meshInfo, meshBox);
  MeshElementTable::setMeshData(meshData, meshInfo);

  setFlagValue(GEOMETRY_TYPE, GEOMETRY_TYPE_MESH, false);
  setFlagValue(DRAW_SOURCE, DRAW_SOURCE_MESH, false);
}


void
Model::resetModelData()
{
  modelData->reset();
  initSplitCombineInfos();
}


void
Model::saveElmerMesh(char* mesh_dir)
{
  outputManager->write_Elmer_mesh(mesh_dir);
}

void
Model::saveElmerPostMesh(char* filename)
{
  // Create outputfile
  ofstream out_file(filename, ios::out);
  outputManager->write_ElmerPost_mesh(out_file);
  out_file.close();
}


void
Model::saveThetisMesh(char* filename)
{
  // Update mesh result filename
  update_dyna_string(modelInfo->meshResultFile, filename);

  // Create outputfile
  ofstream out_file(filename, ios::out);
  outputManager->write_Thetis_mesh(out_file);
  out_file.close();
}


// Save Elmer Front Model File (emf-file)
//
ostream&
Model::saveFrontModelFile(ostream& out, char* filename)
{
  modelInfo->gebhardtFactorsNeedsUpdate = false;
  modelInfo->meshNeedsUpdate = false;
  modelInfo->solverNeedsUpdate = false;
  update_dyna_string(modelInfo->modelFileName, filename);

  return outputManager->emf_output(out, filename);
}


// Save Mesh Input File (mif-file)
//
ostream&
Model::saveMeshInputFile(ostream& out, char* filename)
{
  return outputManager->mif_output(out);
}


// Save Solver Input File (sif-file)
//
ostream&
Model::saveSolverInputFile(ostream& out, char* filename)
{
  return outputManager->sif_output(out);
}



// Save user settings into a (default) file
void
Model::saveUserSettingsFile(char* filename)
{
  Parameter* us = getParameter(0, ECIF_USER_SETTING);

  if (us == NULL) {
    return;
  }

  ofstream out(filename, ios::out);

  char time_buffer[80];
  getCurrentTime(time_buffer);

  out << "!Elmer Cadi default user settings" << endl;
  out << "!" << time_buffer << endl;
  out << endl;

  //-write parameter type "User Settings"
  //-no id
  //-output all
  SifOutputControl soc;
  soc.outputId = false;
  soc.outputAll = true;

  us->output_sif(out, 2, 0, soc);

  out << "End" << endl;

  out.close();
}


void
Model::selectMeshBoundaryElement(int fem_id)
{
  short current_level = meshData->boundaryElements->currentActionLevel;

  if ( getFlagValue(SELECT_MODE_TOGGLE) ) {
    meshData->boundaryElements->selected[fem_id] = !meshData->boundaryElements->selected[fem_id];
  }

  else if ( getFlagValue(SELECT_MODE_EXTEND) ) {
    if ( !meshData->boundaryElements->selected[fem_id] ) {
      meshData->boundaryElements->selected[fem_id] = true;
      meshData->boundaryElements->actionLevels[fem_id] = current_level;
    }
  }

  else if ( getFlagValue(SELECT_MODE_REDUCE) ) {
    if ( meshData->boundaryElements->selected[fem_id] ) {
      meshData->boundaryElements->selected[fem_id] = false;
      meshData->boundaryElements->actionLevels[fem_id] = -current_level;
    }
  }

}


void
Model::unselectMeshBoundaryElement(int fem_id)
{
  short current_level = meshData->boundaryElements->currentActionLevel;

  if ( getFlagValue(SELECT_MODE_TOGGLE) ) {
    meshData->boundaryElements->selected[fem_id] = !meshData->boundaryElements->selected[fem_id];
  }

  else if ( getFlagValue(SELECT_MODE_EXTEND) ) {
    if ( meshData->boundaryElements->selected[fem_id] ) {
      meshData->boundaryElements->selected[fem_id] = false;
      meshData->boundaryElements->actionLevels[fem_id] = -current_level;
    }
  }

  else if ( getFlagValue(SELECT_MODE_REDUCE) ) {
    if ( !meshData->boundaryElements->selected[fem_id] ) {
      meshData->boundaryElements->selected[fem_id] = true;
      meshData->boundaryElements->actionLevels[fem_id] = current_level;
    }
  }

}


void
Model::selectMeshBoundaryElements()
{
  MeshElementTable* table = meshData->boundaryElements;
  table->updateActionLevel(1);

  if ( getFlagValue(SELECT_METHOD_ALL) ) {
    selectMeshBoundaryElementsAll();
  }

  else if ( getFlagValue(SELECT_METHOD_BY_NEIGHBOR) ) {
    selectMeshBoundaryElementsByNeighbor();
  }

  else if ( getFlagValue(SELECT_METHOD_BY_NORMAL) ) {
    selectMeshBoundaryElementsByNormal();
  }

  else if ( getFlagValue(SELECT_METHOD_BY_PLANE) ) {
    selectMeshBoundaryElementsByPlane();
  }
}


void
Model::selectMeshBoundaryElementsAll()
{
  BodyElement* be = getBoundaryById(modelInfo->selectedBodyElementId);

  if (be == NULL) {
    return;
  }

  int nof_elements = be->getNofMeshElements();

  for (int i = 0; i < nof_elements; i++) {
    int index = be->getMeshElementId(i);
    selectMeshBoundaryElement(index);
  }

  refreshRenderer();
}


void
Model::selectMeshBoundaryElementsByNeighbor()
{
  if (modelInfo->selectedBodyElementId == NO_INDEX) {
    return;
  }

  if (meshInfo->selectedBndrElementId == NO_INDEX) {
    return;
  }


  BodyElement* be = getBoundaryById(modelInfo->selectedBodyElementId);

  if (be == NULL) {
    return;
  }

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

  refreshRenderer();
}


void
Model::selectMeshBoundaryElementsByNormal()
{
  if (modelInfo->selectedBodyElementId == NO_INDEX)
    return;

  if (meshInfo->selectedBndrElementId == NO_INDEX)
    return;


  BodyElement* be = getBoundaryById(modelInfo->selectedBodyElementId);
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

  refreshRenderer();
}


void
Model::selectMeshBoundaryElementsByPlane()
{
  if (modelInfo->selectedBodyElementId == NO_INDEX)
    return;

  if (meshInfo->selectedBndrElementId == NO_INDEX)
    return;


  BodyElement* be = getBoundaryById(modelInfo->selectedBodyElementId);
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

  refreshRenderer();
}


void
Model::selectMeshBoundaryElementsRedo()
{
  BodyElement* be = getBoundaryById(modelInfo->selectedBodyElementId);

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

  refreshRenderer();
}


void
Model::selectMeshBoundaryElementsUndo()
{
  BodyElement* be = getBoundaryById(modelInfo->selectedBodyElementId);

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

  refreshRenderer();
}

// Select correct parameter table for the parameter type
//
// NOTE: use_counter is a pointer, because callres may want to
// upadte the corresponding modelStatistic-counter
//
bool
Model::selectParameterTable(ecif_parameterType parameter_type,
                            ParameterTable*& table,
                            int*& use_counter)
{
  switch (parameter_type) {
  case ECIF_BODY_FORCE:
    table = modelData->bodyForces;
    use_counter = &modelStatistics->nofBodyForces;
    break;
  case ECIF_BODY_PARAMETER:
    table = modelData->bodyParameters;
    use_counter = &modelStatistics->nofBodyParameters;
    break;
  case ECIF_BOUNDARY_CONDITION:
    table = modelData->boundaryConditions;
    use_counter = &modelStatistics->nofBoundaryConditions;
    break;
  case ECIF_BOUNDARY_PARAMETER:
    table = modelData->boundaryParameters;
    use_counter = &modelStatistics->nofBoundaryParameters;
    break;
  case ECIF_CALCULATOR:
    table = modelData->calculators;
    use_counter = &modelStatistics->nofCalculators;
    break;
  case ECIF_CONSTANT:
    table = modelData->constants;
    use_counter = &modelStatistics->nofConstants;
    break;
  case ECIF_COORDINATE:
    table = modelData->coordinates;
    use_counter = &modelStatistics->nofCoordinates;
    break;
  case ECIF_DATAFILE:
    table = modelData->datafiles;
    use_counter = &modelStatistics->nofDatafiles;
    break;
  case ECIF_EQUATION:
    table = modelData->equations;
    use_counter = &modelStatistics->nofEquations;
    break;
  case ECIF_EQUATION_VARIABLE:
    table = modelData->equationVariables;
    use_counter = &modelStatistics->nofEquationVariables;
    break;
   case ECIF_GRID_H:
    table = modelData->gridHs;
    use_counter = &modelStatistics->nofGridHs;
    break;
  case ECIF_GRID_PARAMETER:
    table = modelData->gridParameters;
    use_counter = &modelStatistics->nofGridParameters;
    break;
  case ECIF_INITIAL_CONDITION:
    table = modelData->initialConditions;
    use_counter = &modelStatistics->nofInitialConditions;
    break;
  case ECIF_MATERIAL:
    table = modelData->materials;
    use_counter = &modelStatistics->nofMaterials;
    break;
  case ECIF_MODEL_PARAMETER:
    table = modelData->modelParameters;
    use_counter = &modelStatistics->nofModelParameters;
    break;
  case ECIF_SIMULATION_PARAMETER:
    table = modelData->simulationParameters;
    use_counter = &modelStatistics->nofSimulationParameters;
    break;
  case ECIF_SOLVER:
    table = modelData->solvers;
    use_counter = &modelStatistics->nofSolvers;
    break;
  case ECIF_SOLVER_CONTROL:
    table = modelData->solverControls;
    use_counter = &modelStatistics->nofSolverControls;
    break;
  case ECIF_TIMESTEP:
    table = modelData->timesteps;
    use_counter = &modelStatistics->nofTimesteps;
    break;
  case ECIF_USER_SETTING:
    table = modelData->userSettings;
    use_counter = &modelStatistics->nofUserSettings;
    break;
  default:
    return false;
  }
  return true;
}


void
Model::separateBodies()
{
  bool* duplicated_object_flags = new bool[modelStatistics->nofModelObjects];

  for (int i = 0; i < modelStatistics->nofModelObjects; i++) {
    duplicated_object_flags[i] = false;
  }

  int index = 0;
  while (true) {
    Body* body = getBody(index++);
    if (body==NULL) break;
    if ( body->getSeparateFlag() ) {
      body->separate(duplicated_object_flags);
    }
  }
}


void
Model::setActiveMeshIndices(int nof_meshes,  int* mesh_indices)
{
  delete[] modelInfo->activeMeshIndices;
  modelInfo->activeMeshIndices = NULL;

  modelInfo->nofActiveMeshes = nof_meshes;

  if ( nof_meshes == 0 )
    return;

  modelInfo->activeMeshIndices = new int[nof_meshes];

  for (int i = 0; i < nof_meshes; i++) {
    modelInfo->activeMeshIndices[i] = mesh_indices[i];
  }
}


void
Model::setBodyColors(ColorIndexArray& unused_colors)
{
  int index = 0;
  while (unused_colors.size() != 0) {
    Body* body = getBody(index++);
    if (body==NULL) break;
    if (DEFAULT_COLOR_INDEX == body->getColorIndex()) {
      body->setColorIndex(unused_colors.back());
      unused_colors.pop_back();
    }
  }
}


// Set parent tags and layer indices for boundaries
//
void
Model::setBoundaryParentIdsAndTags()
{
  // Loop all bodies
  int index = 0;
  while (true) {

    Body* body = getBody(index++);

    if (body==NULL) break;

    int layer = -1;

    // Loop all layer and their elements in the body
    //
    while (true) {

      if (!body->selectLayer(++layer)) break;

      int be_index = 0;

      while (true) {

        BodyElement* be = body->getElement(layer, be_index++);

        if (be == NULL) break;

        int bd_id = body->Id();
        int bd_tg = body->Tag();

        int bd1_id = be->getParentId(1);
        int bd2_id = be->getParentId(2);

        // If body is already a parent this must be an intra layer boundary!
        if ( bd1_id == bd_id || bd2_id == bd_id ) {
          be->setIsIntraLayerBoundary(true);
          continue;
        }

        short body_nbr;

        // If body1-id not yet set this is set
        // as the first body
        if ( bd1_id == NO_INDEX ) {
          body_nbr = 1;

        // If body1-id already set (via an other body)
        // this is the second body
        } else {
          body_nbr = 2;
        }

        be->setParentId(body_nbr, bd_id);
        be->setParentTag(body_nbr, bd_tg);
        be->setParentLayer(body_nbr, layer);

      } // Each body element in the layer
    } // Each layer
  }
}

#if 0
// Set boundary element boundary condtions from boundary groups
// NOTE: Boundary condityions are set in hierarchical order!
void
Model::setBoundaryConditions()
{
  ElementGroup* eg;
  BodyElement* be;
  int i, index;

  // Face groups
  index = 0;
  while (true) {
    eg = getElementGroup(index++);
    if (eg==NULL) break;
    if ( OT_FACE != eg->getElementType() ) continue;
    int bc_pid = eg->getBoundaryConditionId();
    for (i = 0; i < eg->getNofElements(); i++) {
      be = (BodyElement*)eg->getElement(i);
      if ( be == NULL ) continue;
      be->setBoundaryConditionId(bc_pid);
    }
  }

  // Edge groups
  index = 0;
  while (true) {
    eg = getElementGroup(index++);
    if (eg==NULL) break;
    if ( OT_EDGE != eg->getElementType() ) continue;
    int bc_pid = eg->getBoundaryConditionId();
    for (i = 0; i < eg->getNofElements(); i++) {
      be = (BodyElement*)eg->getElement(i);
      if ( be == NULL ) continue;
      be->setBoundaryConditionId(bc_pid);
    }
  }

  // Vertex groups
  index = 0;
  while (true) {
    eg = getElementGroup(index++);
    if (eg==NULL) break;
    if ( OT_VERTEX != eg->getElementType() ) continue;
    int bc_pid = eg->getBoundaryConditionId();
    for (i = 0; i < eg->getNofElements(); i++) {
      be = (BodyElement*)eg->getElement(i);
      if ( be == NULL ) continue;
      be->setBoundaryConditionId(bc_pid);
    }
  }

}
#endif


// Set parent tags for boundaries
void
Model::setBoundaryTags()
{
  // Boundary tag counter
  int max_btag = 0;

  bool only_active = true;
  BodyElement* be;
  int index;

  // Normal boundaries use their own tags!
  //
  index = 0;
  while (true) {
    be = getBoundary(index++, only_active);
    if (be==NULL) break;
    int tag = be->Tag();

    be->setBoundaryTag(tag);

    if (max_btag < tag)
      max_btag = tag;
  }

  // Edges in 3D will be numbered simply after boundaries
  //
  if ( modelInfo->dimension == ECIF_3D ) {
    index = 0;
    while (true) {
      be = getEdge(index++, only_active);
      if (be==NULL) break;
      be->setBoundaryTag(++max_btag);
    }
  }

  // Vertices are numbered simply after boundaries or  edges
  //
  index = 0;
  while (true) {
    be = getVertex(index++, only_active);
    if (be==NULL) break;
    be->setBoundaryTag(++max_btag);
  }

  // Update class-variable
  BodyElement::setLastBoundaryTag(max_btag);
}


// Set boundary point tags and mesh density values
void
Model::setBoundaryPointData()
{

  int bp_counter = 0;
  int next_bp_tag = 1 + modelStatistics->nofVertices;

  int index = 0;
  while (true) {

    BodyElement* be = getBoundary(index++);

    if (be==NULL) break;
    
    be->checkBoundaryDiscretization(modelInfo->currentMeshIndex);

    int bp_count;
    BoundaryPoint** bp_points;

    be->getBoundaryPoints(bp_count, bp_points);

    int i;

    // Reset tags (this way we can control that
    // copied point (for closed boundaried) are
    // not numbered twice
    for (i = 0; i < bp_count; i++) {

      BoundaryPoint* bp = bp_points[i];

      if ( bp->isVertex() ) {
        continue;
      }

      bp->checked = false;
      bp->activeInMeshing = false;
      bp->tag = NO_INDEX;
    }

    // Now number points (but skip those which already have a tag!)
    for (i = 0;  i < bp_count; i++) {

      BoundaryPoint* bp = bp_points[i];

      if ( bp->isVertex() || bp->checked ) {
        continue;
      }

      bp_counter++;
      bp->tag = next_bp_tag++;
      bp->checked = true;
    }

    // Finally set possible density values for boundary points
    be->setBoundaryPointsMeshDensityValue(getCurrentMeshIndex());

  }

  modelStatistics->nofBoundaryPoints = bp_counter;
}


void
Model::setCurrentAnchorPoint(int vertex_id)
{
  BodyElement* vertex = getVertexById(vertex_id);
  GcPoint* point = (GcPoint*)vertex->getGeometry();

  Point3* p = point->getPoint();

  for (int i = 0; i < 3; i++) {
    pickInfo->anchorPoint[i] = (*p)[i];
  }
}


void
Model::setCurrentAnchorPoint(Point3 point)
{
  for (int i = 0; i < 3; i++) {
    pickInfo->anchorPoint[i] = point[i];
  }
}


void
Model::setCurrentPickInfo(Point3 point, Point3 dir)
{
  for (int i = 0; i < 2; i++) {
    pickInfo->pickPoint[i] = point[i];
    pickInfo->pickDir[i] = dir[i];
  }

  pickInfo->pickPoint[2] = 0.0;
  pickInfo->pickDir[2] = 0.0;

  refreshRenderer();
}


void
Model::drawCurrentPickVector()
{
  Renderer* renderer = theControlCenter->getRenderer();

  renderer->drawLine(DM_NORMAL, DS_NORMAL, 1, &(pickInfo->anchorPoint), &(pickInfo->pickPoint));
}


void
Model::setCurrentMeshIndex(int index)
{
  if ( index < 0 || index >= modelInfo->nofMeshes )
    modelInfo->currentMeshIndex = NO_INDEX;

  else
    modelInfo->currentMeshIndex = index;
}


void
Model::setInitialMeshH()
{
  UserInterface* ui = theControlCenter->getGui();

  double initial_mesh_h = calcInitialMeshH();

  if (initial_mesh_h < MAX_RANGE)
    ui->setInitialMeshH(initial_mesh_h);
  else
    ui->setInitialMeshH(-1.0);
}


void
Model::setFlagValue(flagGroup group, flagName name, bool value)
{
  bool refresh_renderer = false;

  // Reset exlusive flag sets
  switch (group) {

  case DRAW_SOURCE:
    refresh_renderer = true;
    break;

  case DRAW_TARGET:
    modelFlags[DRAW_TARGET_BODIES]      = false;
    modelFlags[DRAW_TARGET_SURFACES]    = false;
    modelFlags[DRAW_TARGET_EDGES]       = false;
    modelFlags[DRAW_TARGET_BOUNDARIES]  = false;
    refresh_renderer = true;
    break;

  case SELECT_METHOD:
    modelFlags[SELECT_METHOD_SINGLE]        = false;
    modelFlags[SELECT_METHOD_ALL]           = false;
    modelFlags[SELECT_METHOD_BY_NEIGHBOR]   = false;
    modelFlags[SELECT_METHOD_BY_NORMAL]     = false;
    modelFlags[SELECT_METHOD_BY_PLANE]      = false;
    modelFlags[SELECT_METHOD_BY_BOX]        = false;
    modelFlags[SELECT_METHOD_BY_RECTANGLE]  = false;
    break;

  case SELECT_MODE:
    modelFlags[SELECT_MODE_TOGGLE] = false;
    modelFlags[SELECT_MODE_EXTEND] = false;
    modelFlags[SELECT_MODE_REDUCE] = false;
    //refresh_renderer = true;
    break;

  case SELECT_OBJECTS:
    modelFlags[SELECT_OBJECTS_TOGGLE] = false;
    modelFlags[SELECT_OBJECTS_EXTEND] = false;
    //refresh_renderer = true;
    break;

  case LABEL_DISPLAY:
    // No exclusive values here!
    //modelFlags[LABEL_DISPLAY_NODE]     = false;
    //modelFlags[LABEL_DISPLAY_ELEMENT]  = false;
    //modelFlags[LABEL_DISPLAY_VERTEX]   = false;
    //modelFlags[LABEL_DISPLAY_EDGE]     = false;
    //modelFlags[LABEL_DISPLAY_FACE]     = false;
    //modelFlags[LABEL_DISPLAY_BODY]     = false;
    //modelFlags[LABEL_DISPLAY_BODYPAIR] = false;
    refresh_renderer = true;
    break;
  }

  // Set flag value
  modelFlags[name] = value;

  // Postprocess flags
  switch (group) {

  case DRAW_TARGET: // Check if we are drawing boundaries
    if ( modelInfo->dimension == ECIF_2D ) {
      if ( name == DRAW_TARGET_EDGES ) {
        modelFlags[DRAW_TARGET_BOUNDARIES] = true;
      }
    }
    if ( modelInfo->dimension == ECIF_3D ) {
      if ( name == DRAW_TARGET_SURFACES ) {
        modelFlags[DRAW_TARGET_BOUNDARIES] = true;
      }
    }
    break;

  case SELECT_METHOD:
    //selectMeshBoundaryElements();
    break;

  case SELECT_OBJECTS:
    // Set deafult value if no value set
    if ( modelFlags[SELECT_OBJECTS_TOGGLE] == false &&
         modelFlags[SELECT_OBJECTS_EXTEND] == false
       )
      modelFlags[SELECT_OBJECTS_TOGGLE] = true;
    break;
  }

  if ( refresh_renderer ) {
    refreshRenderer();
  }

}


void
Model::setMatcInputFileEmf(char* matc_input_file)
{
  update_dyna_string(modelInfo->matcInputFile_emf, matc_input_file);
}


void
Model::setMatcInputFileSif(char* matc_input_file)
{
  update_dyna_string(modelInfo->matcInputFile_sif, matc_input_file);
}


void
Model::setMeshBgMeshActives(int nof_files, int* bg_mesh_actives)
{
  delete[] modelInfo->meshBgMeshActives;
  modelInfo->meshBgMeshActives = NULL;

  if ( nof_files > 0 ) {

    // Alloc table
    modelInfo->meshBgMeshActives = new bool[nof_files];

    for (int i = 0; i < nof_files; i++) {
      modelInfo->meshBgMeshActives[i] = (bool) bg_mesh_actives[i];
    }
  }
}


void
Model::setMeshBgMeshControls(int nof_files, int* bg_mesh_controls)
{
  delete[] modelInfo->meshBgMeshControls;
  modelInfo->meshBgMeshControls = NULL;

  if ( nof_files > 0 ) {

    // Alloc table
    modelInfo->meshBgMeshControls = new bool[nof_files];

    for (int i = 0; i < nof_files; i++) {
      modelInfo->meshBgMeshControls[i] = (bool) bg_mesh_controls[i];
    }
  }
}


void
Model::setMeshBgMeshFiles(int nof_files, char** bg_mesh_files)
{
  delete[] modelInfo->meshBgMeshFiles;
  modelInfo->meshBgMeshFiles = NULL;

  if ( nof_files > 0 ) {

    // Alloc table
    modelInfo->meshBgMeshFiles = new char*[nof_files];

    for (int i = 0; i < nof_files; i++) {

      modelInfo->meshBgMeshFiles[i] = NULL;

      if ( bg_mesh_files[i] == NULL ) {
        continue;
      }

      int len = strlen(bg_mesh_files[i]);

      if ( len == 0 ) {
        continue;
      }

      // Alloc space for each mesh name
      modelInfo->meshBgMeshFiles[i] = new char[1 + len];

      strcpy(modelInfo->meshBgMeshFiles[i], bg_mesh_files[i]);
    }
  }

}


void
Model::setMeshBgMeshFileIndices(int nof_files, int* bg_mesh_file_indices)
{
  int i;

  delete[] modelInfo->meshBgMeshFileIndices;

  modelInfo->nofBgMeshFiles = nof_files;
  modelInfo->meshBgMeshFileIndices = NULL;

  if ( nof_files > 0 ) {

    // Alloc table
    modelInfo->meshBgMeshFileIndices = new int[nof_files];

    for (i = 0; i < nof_files; i++) {
      modelInfo->meshBgMeshFileIndices[i] = bg_mesh_file_indices[i];
    }
  }
}


void
Model::setMeshFs(int nof_values, double* mesh_fs)
{
  delete[] modelInfo->meshFs;
  modelInfo->meshFs = NULL;

  if ( nof_values == 0 )
    return;

  modelInfo->meshFs = new double[nof_values];

  for (int i = 0; i < nof_values; i++) {
    modelInfo->meshFs[i] = mesh_fs[i];
  }
}


void
Model::setMeshHs(int nof_values, double* mesh_hs)
{
  delete[] modelInfo->meshHs;
  modelInfo->meshHs = NULL;

  if ( nof_values == 0 )
    return;

  modelInfo->meshHs = new double[nof_values];

  for (int i = 0; i < nof_values; i++) {
    modelInfo->meshHs[i] = mesh_hs[i];
  }
}


void
Model::setMeshInputUnit(double unit)
{
  meshInputUnit = unit;
}


void
Model::setMeshNames(int nof_meshes, char** mesh_names)
{
  int i;

  for (i = 0; i < modelInfo->nofMeshes; i++) {
    delete[] modelInfo->meshNames[i];
  }
  delete[] modelInfo->meshNames;

  modelInfo->meshNames = NULL;
  modelInfo->nofMeshes = nof_meshes;

  if ( nof_meshes > 0 ) {

    // Alloc table
    modelInfo->meshNames = new char*[nof_meshes];

    for (i = 0; i < nof_meshes; i++) {

      int len = strlen(mesh_names[i]);

      // Alloc space for each mesh name
      modelInfo->meshNames[i] = new char[1 + len];

      strcpy(modelInfo->meshNames[i], mesh_names[i]);
    }
  }

}


void
Model::setMeshNodes()
{
  MeshElementTable::setMeshNodes(meshData->nodes);
}


void
Model::setModelDescriptions(char* model_desc, char* problem_desc)
{
  update_dyna_string(modelInfo->modelDescription, model_desc);
  update_dyna_string(modelInfo->problemDescription, problem_desc);
}


void
Model::setModelDimension(ecif_modelDimension dim)
{
  modelInfo->dimension = dim;
}


void
Model::setModelFileDirectories(char* model_dir, char* include_path, char* results_dir, char* temporary_dir)
{
  update_dyna_string(modelInfo->modelDirectory, model_dir);
  update_dyna_string(modelInfo->includePath, include_path);
  update_dyna_string(modelInfo->resultsDirectory, results_dir);
  update_dyna_string(modelInfo->temporaryFilesDirectory, temporary_dir);
}


void
Model::setModelFileDirectoriesAbs(char* model_dir, char* include_path, char* results_dir, char* temporary_dir)
{
  update_dyna_string(modelInfo->modelDirectory_absolute, model_dir);
  update_dyna_string(modelInfo->includePath_absolute, include_path);
  update_dyna_string(modelInfo->resultsDirectory_absolute, results_dir);
  update_dyna_string(modelInfo->temporaryFilesDirectory_absolute, temporary_dir);
}

void
Model::setModelFileDirectoriesSave(bool include_path, bool results_dir, bool temporary_dir)
{
  modelInfo->includePath_save = include_path;
  modelInfo->resultsDirectory_save = results_dir;
  modelInfo->temporaryFilesDirectory_save = temporary_dir;
}


void
Model::setModelFileCreated(char* created_str)
{
  update_dyna_string(modelInfo->created, created_str);
}


void
Model::setModelFileModified(char* modified_str)
{
  update_dyna_string(modelInfo->modified, modified_str);
}


void
Model::setModelMeshDirectory(char* mesh_dir)
{
  update_dyna_string(modelInfo->meshDirectory, mesh_dir);
}


void
Model::setModelMeshDirectoryAbs(char* mesh_dir)
{
  update_dyna_string(modelInfo->meshDirectory_absolute, mesh_dir);
}


void
Model::setModelNameAndDirectory(char* file_path)
{
  char name_buffer[1024];
  char dir_buffer[1024];

  int len = strlen(file_path);

  for (int i = 0; i< len; i++) {
    if ( file_path[i] == '\\' )
      file_path[i] = '/';
  }

  int counter;
  int pos = len -1;

  // Move at the end of name
  while ( pos >= 0 && file_path[pos] != '/' )
    pos--;

  // Read name (in reversed order!)
  if ( pos >= 0 && file_path[pos] == '/' )
    pos--;
  counter = 0;
  while ( pos >= 0 && file_path[pos] != '/' ) {
    name_buffer[counter++] = file_path[pos--];
  }
  name_buffer[counter] = '\0';

  // Read dir (in reversed order!)
  if ( pos >= 0 && file_path[pos] == '/' )
    pos--;
  counter = 0;
  while ( pos >= 0 ) {
    dir_buffer[counter++] = file_path[pos--];
  }
  dir_buffer[counter] = '\0';


  //- Reverse name buffer
  len = strlen(name_buffer);
  counter = 0;
  while (counter < len / 2) {
    char tmp = name_buffer[counter];
    name_buffer[counter] = name_buffer[len - 1 - counter];
    name_buffer[len - 1 - counter] = tmp;
    counter++;
  }

  //- Reverse dir buffer
  len = strlen(dir_buffer);
  counter = 0;
  while (counter < len / 2) {
    char tmp = dir_buffer[counter];
    dir_buffer[counter] = dir_buffer[len - 1 - counter];
    dir_buffer[len - 1 - counter] = tmp;
    counter++;
  }

  update_dyna_string(modelInfo->modelName, name_buffer);
  update_dyna_string(modelInfo->modelDirectory, dir_buffer);

  setWindowTitles();
}


void
Model::setModelNames(char* model_name, char* problem_name)
{
  update_dyna_string(modelInfo->modelName, model_name);
  update_dyna_string(modelInfo->problemName, problem_name);

  setWindowTitles();
}


bool
Model::setModelObjectNameById(int oid, char*name)
{
  ModelObject* obj = getModelObjectById(oid);

  if (obj == NULL ) return false;

  obj->setName(name);

  return true;
}


bool
Model::setModelObjectTagById(int oid, int tag)
{
  ModelObject* obj = getModelObjectById(oid);

  if (obj == NULL ) return false;

  obj->setTag(tag);

  return true;
}


// Method updates current model outfile time
void
Model::setModelFileTime(char* outfile_time)
{
  update_dyna_string(modelInfo->modelFileTs, outfile_time);
}


void
Model::setModelHasDiffuseGrayRadiation(bool value)
{
  modelInfo->hasDiffuseGrayRadiation = value;
}


void
Model::setSelectionsToGui()
{
  bool has_selections = false;

  UserInterface* gui = theControlCenter->getGui();

  int index = 0;
  while (true) {

    BodyElement* be = getBoundary(index++);

    if (be==NULL) break;

    if ( DS_SELECTED == be->getDrawState() ) {

      int bndr_id = be->Id();

      // If first selected bodyelement was found
      if (!has_selections) {
        has_selections = true;

        int bd1_id = modelInfo->selectedBody1Id;
        int lr1_id = modelInfo->selectedLayer1Id;

        int bd2_id = modelInfo->selectedBody2Id;
        int lr2_id = modelInfo->selectedLayer2Id;

        gui->selectBoundary(bndr_id, bd1_id, lr1_id, bd2_id, lr2_id, true);
      }

      gui->setBoundarySelectionMode(bndr_id, true);
    }

  }
}


// Method updates timestamp for the "target_name" related data.
void
Model::setTimestamp(char* target_name, char* ts)
{
  // Predefined target names
  const char* target_names[7]= {
    "Database",
    "Front",
    "Gebhardt Factors",
    "Grid Parameter",
    "Mesh",
    "Solver",
    "View Factors"
  };

  // Predefined target variables
  char** targets[7] = {
    &modelInfo->databaseTs,
    &modelInfo->frontTs,
    &modelInfo->gebhardtFactorsTs,
    &modelInfo->meshParameterTs,
    &modelInfo->meshTs,
    &modelInfo->solverTs,
    &modelInfo->viewfactorsTs
  };

  // Check which of the targets (if any) should
  // be updated
  for (int i = 0; i < 7; i++) {
    if ( LibFront::ncEqual(target_name, (char*)target_names[i]) ) {
      char*& target = *(targets[i]);
      update_dyna_string(target, ts);
      break;
    }
  }
}


void
Model::setVertexTable(int dim1, int dim2, int* vertex_ids, MatcValueTable& matc_table)
{
  modelData->setVertexTable(dim1, dim2, vertex_ids, matc_table);
}


void
Model::setWindowTitles()
{
  strstream strm;

  strm << "ELMER Front - ";

  // If model exists
  if ( modelInfo->modelName != NULL && modelInfo->modelName[0] != '\0' ) {
    strm << modelInfo->modelName;

    // if problem name defined
    if ( modelInfo->problemName != NULL && modelInfo->problemName[0] != '\0' ) {
      strm << "." << modelInfo->problemName;
    }

  }
  // No model exists
  else {
    strm << "No model name";
  }

  strm << ends;

  theControlCenter->setWindowTitles(strm.str());
}


void
Model::setParameter(ecif_parameterType parameter_type,
                    int pid, int parent_id,
                    char* data_string, char* param_name)
{
  Parameter* parameter = getParameterById(parameter_type, pid);

  // Create a new parameter entry and add to model
  if (parameter == NULL) {
    parameter = createNewParameter( parameter_type, pid, parent_id,
                                    data_string, param_name);
    addParameter(parameter_type, parameter);
  }

  // Update an existing parameter
  else {
    parameter->setName(param_name);
    parameter->setValue(data_string);
  }

  // Mark this parameter updated
  parameter->setUpdateFlag(true);

  UserInterface* gui = theControlCenter->getGui();

  char ts_buffer[64];

  theControlCenter->getCurrentTimestamp(ts_buffer);

  // Update coordinate system info
  if (parameter_type == ECIF_COORDINATE) {

    ParameterField* pf;
    char** data = NULL;

    pf = parameter->getFieldBySifName("Coordinate System");

    if(pf != NULL) 
    {
      data = pf->getDataStrings();

      modelInfo->updateCoordinateType(data[0]);
      modelInfo->updateSimulationDimension(data[0]);

    } else {
      printf("Parameter field is null, something went wrong. Contact <elmerdiscussion@postit.csc.fi>\n");
    }


    pf = parameter->getFieldBySifName("Coordinate Mapping");

    if(pf != NULL) 
    {
      data = pf->getDataStrings();
      modelInfo->updateCoordinateMapping(data[0]);
    } else {
      printf("Parameter field is null, something went wrong. Contact <elmerdiscussion@postit.csc.fi>\n");
    }

    Renderer* renderer = theControlCenter->getRenderer();

    if ( renderer != NULL ) {
      renderer->setSimulationDimension(modelInfo->simulationDimension);
      refreshRenderer();
    }

  }

  // Write mesh parameter ts
  if (parameter_type == ECIF_GRID_PARAMETER) {
    setTimestamp("Mesh_Parameter", ts_buffer);
    gui->setTimestamp(ECIF_GRID_PARAMETER, ts_buffer);
  }

  // NOTE: Change checks do not work anymore, because
  // all parameter are written anew each time parameters
  // are updated (old are deleted and existing are added as new!)

#if 0
  // Check change status
  bool value_changed;
  bool emissivity_changed;

  parameter->getParameterValueState(value_changed);

  // If any parameter's value has changed (potential change for Solver)
  if (value_changed) {
    modelInfo->solverNeedsUpdate = true;
    gui->setNeedsUpdate(EFN_PROGRAM_SOLVER);
  }

  // Check if mesh parameter has changed
  if (parameter_type == ECIF_MESH && value_changed) {
    setTimestamp("Mesh_Parameter", ts_buffer);
    gui->setTimestamp(ECIF_MESH, ts_buffer);
    //modelInfo->meshNeedsUpdate = true;
    //gui->setNeedsUpdate(EFN_PROGRAM_MESH);
  }

  // Check if emissivity boundary condition has changed
  if (parameter_type == ECIF_BOUNDARY_CONDITION) {

    parameter->getParameterFieldValueState(EFN_EMISSIVITY, emissivity_changed);

    if (emissivity_changed) {
      modelInfo->gebhardtFactorsNeedsUpdate = true;
      gui->setNeedsUpdate(EFN_PROGRAM_GEBHARDT_FACTORS);
    }
  }
#endif


}


int
Model::setNewSplitCombineIndex(short shift)
{
  if (shift == 0)
    return modelData->splitCombineInfoIndex;

  int index = modelData->splitCombineInfoIndex + shift;

  if ( index < 0 )
    return modelData->splitCombineInfoIndex;

  // Allocate more space if needed
  if ( index >= MAX_NOF_SPLIT_COMBINE_INFOS - 1 ) {

    int current_size = MAX_NOF_SPLIT_COMBINE_INFOS;

    MAX_NOF_SPLIT_COMBINE_INFOS += 16;

    SplitCombineInfoArray* tmp = new SplitCombineInfoArray(MAX_NOF_SPLIT_COMBINE_INFOS);

    for (int i = 0; i < MAX_NOF_SPLIT_COMBINE_INFOS; i++) {
      if ( i < current_size)
        (*tmp)[i] = (*modelData->splitCombineInfos)[i];
      else
        (*tmp)[i] = NULL;
    }

    delete modelData->splitCombineInfos;

    modelData->splitCombineInfos = tmp;
  }

  modelData->splitCombineInfoIndex = index;

  return modelData->splitCombineInfoIndex;
}


void
Model::setParallelInfo(ParallelInfo& pi)
{
  parallelInfo->nofProcessors = pi.nofProcessors;
}


#if 0
void
Model::sortBoundaryVertices(enum modelGeometryType gtype)
{

  BoundaryVertexList* bvlist;

  // Select vertex list to use
  if (gtype == GEOM_CAD) {
    bvlist = modelData->boundaryVertices;
  } else {
    bvlist = modelData->meshBoundaryVertices;
  }

  int count = 1 + bvlist->size();

  BoundaryVertex** bpoints = new BoundaryVertex*[count];

  for (int i = 0; i < count; i++) {
    bpoints[i] = NULL;
  }

  BoundaryVertex* bp = getFirstBoundaryVertex(gtype);

  while (bp != NULL ) {

    bpoints[bp->id] = bp;

    bp = getNextBoundaryVertex(gtype);
  }

  bvlist->clear();

  for (i = 0; i < count; i++) {
    if ( bpoints[i] != NULL ) {
      bvlist->push_back(bpoints[i]);
    }
  }

  if (gtype == GEOM_CAD) {
    modelStatistics->nofBoundaryVertices = bvlist->size();
  } else {
    modelStatistics->nofMeshBoundaryVertices = bvlist->size();
  }

}
#endif

#if 0
void
Model::sortMeshBoundaryIndices()
{
  int i, j;
  int nof_belems = meshInfo->nofInnerBndrElements + meshInfo->nofOuterBndrElements;
  if (nof_belems == 0)
    return;
  // Sorted indices are stored here
  int** new_table = new int*[nof_belems];
  //
  // Position counter for each (body1,body2)-pair
  // are stored in this table
  int start = 0;
  int** pair_counters = new int*[meshInfo->nofBodies];
  for (i = 0; i < meshInfo->nofBodies; i++) {
    pair_counters[i] = new int[meshInfo->nofBodies];
    for (j = 0; j < meshInfo->nofBodies; j++) {
      pair_counters[i][j] = start;
      start += meshBodyPairTable[i][j];
    }
  }
  // Sort data
  int M = meshInfo->nofBodies; // pair table size
  for (i=0; i < nof_belems; i++) {
    int p1_id = meshBoundaryIndices[i][BELM_PARENT1_ID];
    int p2_id = meshBoundaryIndices[i][BELM_PARENT2_ID];
    int m1 = meshBulkElements[p1_id][ELM_BODY];
    int m2 = meshBulkElements[p2_id][ELM_BODY];
    if (m2 == NO_INDEX)
      m2 = m1;
    // Boundary element ID:
    // pair-table-index + 1(index origo = 1) - lower diag indices
    int bndr_elem_id = (m1 * M + m2) + 1 - ( (m1+1) * ((m1+0)/2) );
    //meshBoundaryIndices[i][BELM_ID] = m1 * M + m2 + 1;
    meshBoundaryIndices[i][BELM_ID] = bndr_elem_id;
    int index = pair_counters[m1][m2];
    new_table[index] = meshBoundaryIndices[i];
    pair_counters[m1][m2]++;
  }

  for (i = 0; i < meshInfo->nofBodies; i++)
    delete[] pair_counters[i];
  delete[] pair_counters;

  delete[] meshBoundaryIndices;
  meshBoundaryIndices = new_table;
}
#endif


// NOTE: Only one boundary element at a time
// can be under splitting !!!
void
Model::splitBoundary(int body1_id, int body2_id)
{
  UserInterface* gui = theControlCenter->getGui();

  Body* body1 = getBodyById(body1_id);
  Body* body2 = getBodyById(body2_id);

  int body1_tag = NO_INDEX;
  int body2_tag = NO_INDEX;

  if ( body1 != NULL ) {
    body1_tag = body1->Tag();
  }

  if ( body2 != NULL ) {
    body2_tag = body2->Tag();
  }

  int nof_selected;
  int source_id;

  // Split source element (boundary)
  BodyElement* split_source = NULL;

  //---Count nof mesh elements in the new body element
  //   Only one boundary at a time!!!
  int index = 0;
  while (true) {
    BodyElement* be = getBoundary(index++);
    if (be==NULL) break;
    nof_selected = be->getNofMeshSelectedElements();
    if (nof_selected > 0) {
      split_source = be;
      source_id = be->Id();
      break;
    }
  }

  //---If nothing was selected
  if (split_source == NULL) {
    gui->showMsg("Nothing selected for splitting!");
    return;
  }

  int nof_source_elements = split_source->getNofMeshElements();

  // If all source elements were selected, split
  // does not make sense!
  if (nof_selected == nof_source_elements) {
    gui->showMsg("All elements were selected. Nothing to split!");
    return;
  }

  if (modelInfo->editingMeshBoundaries == false) {
    startEditMeshBoundaries();
  }

  // Delete possible old info entry and replace with a new
  short idx = modelData->splitCombineInfoIndex;
  delete (*modelData->splitCombineInfos)[idx];
  SplitCombineInfo* sc_info = new SplitCombineInfo;
  (*modelData->splitCombineInfos)[idx] = sc_info;

  // Update index for the next entry
  setNewSplitCombineIndex(+1);

  int body1_lr = split_source->getParentLayer(1);
  int body2_lr = split_source->getParentLayer(2);

  //---Create new bodyelement, add to the model and to the bodies
  BodyElement* new_be = createBodyElement(body1_tag, body1_lr, body2_tag, body2_lr, nof_selected);

  new_be->setParentIds(body1_id, body2_id);
  new_be->checkElementGroupData();

  (*modelData->createdModelElements)[new_be->Id()] = new_be;

  //-Store split info
  sc_info->canRedo = false;
  sc_info->canUndo = true;
  sc_info->splitSourceId = source_id;
  sc_info->splitTargetId = new_be->Id();
  sc_info->nofSplitTargetMeshElements = nof_selected;
  sc_info->splitTargetMeshElementSourceIndices = new int[nof_selected];

  body1->check();

  if (body2 != NULL) {
    body2->check();
  }

  //---Copy selected mesh elements to the new boundary
  //   and remove them from the original boundary
  split_source->setDrawState(DS_NORMAL);
  gui->setBoundarySelectionMode(source_id, false);

  //-Copy all selected
  for (int i = 0; i < split_source->getNofMeshElements(); i++) {

    int elem_id = split_source->getMeshElementId(i);
    short elem_dir = split_source->getMeshElementDir(i);

    if ( !meshData->boundaryElements->selected[elem_id] )
      continue;

    // Set mesh element to the new boundary
    new_be->addMeshElement(elem_id, elem_dir);
    meshData->boundaryElementBoundaryIds[elem_id] = sc_info->splitTargetId;
  }

  //-Remove (selected) elements from the original
  // Store also original source ids for possible redo!!!
  int* removed = sc_info->splitTargetMeshElementSourceIndices;
  split_source->removeSelectedMeshElements(removed);

  //-Update model
  split_source->findMeshBorder();
  new_be->findMeshBorder();

  updateModelStatistics();

  gui->updateBodyData(this);
  gui->updateBoundaryData(this);

  gui->setMeshEdited();

  meshData->boundaryElements->resetSelected();

  refreshRenderer();
}


void
Model::splitCombineBoundariesRedo()
{
  UserInterface* gui = theControlCenter->getGui();

  short cur_index = modelData->splitCombineInfoIndex;

  SplitCombineInfo* sc_info = (*modelData->splitCombineInfos)[cur_index];

  // If there is nothing to undo
  if (sc_info == NULL || !sc_info->canRedo) {
    gui->showMsg("Nothing for redo!");
    return;
  }

  // Update index for the next insert
  int new_index = setNewSplitCombineIndex(+1);

  sc_info->canRedo = false;
  sc_info->canUndo = true;

  //---If we have to REDO a SPLIT
  if (sc_info->splitTargetId != NO_INDEX) {

    BodyElement* be_target = getRemovedBodyElement(sc_info->splitTargetId);

    // This is not a "removed" element anymore
    modelData->removedModelElements->erase(be_target->Id());

    // Add target element (again) to the model
    restoreBodyElement(be_target, true);

    BodyElement* be_source = getBoundaryById(sc_info->splitSourceId);

    // Remove target mesh elements from the source
    // NOTE: They are still stored in the target, no need to
    // add them there again!
    be_source->removeMeshElements(sc_info->nofSplitTargetMeshElements,
                                  sc_info->splitTargetMeshElementSourceIndices);

    // Update moved mesh element boundary ids
    for (int i = 0; i < sc_info->nofSplitTargetMeshElements; i++) {

      int elem_id = be_target->getMeshElementId(i);
      meshData->boundaryElementBoundaryIds[elem_id] = sc_info->splitTargetId;
    }

    be_source->findMeshBorder();
    be_target->findMeshBorder();
  }

  //---If we have to REDO a COMBINE
  else if (sc_info->nofCombineSourceIds > 0) {

    BodyElement* be_target = getRemovedBodyElement(sc_info->combineTargetId);

    // This should not be a "removed" element anymore
    modelData->removedModelElements->erase(be_target->Id());

    // Add target element (again) to the model
    restoreBodyElement(be_target, true);

    // Remove source elements from the model
    for (int i = 0; i < sc_info->nofCombineSourceIds; i++) {
      BodyElement* be_source = getBoundaryById(sc_info->combineSourceIds[i]);
      removeBodyElement(be_source, true, false);
    }

    be_target->findMeshBorder();
  }

  //-Check body data (rebuild element loops etc.)
  checkBodies();

  //-Update model
  updateModelStatistics();

  // Update gui and graphics screen
  gui->updateBodyData(this);
  gui->updateBoundaryData(this);

  refreshRenderer();

#if 0
  // For debugging
  strstream cur_strm;
  strstream new_strm;
  cur_strm << "Cur SplitCombineIndex = " << cur_index << ends;
  new_strm << "New SplitCombineIndex = " << new_index << ends;
  gui->showMsg(cur_strm.str());
  gui->showMsg(new_strm.str());
#endif

}


void
Model::splitCombineBoundariesUndo()
{
  UserInterface* gui = theControlCenter->getGui();

  // For debugging only!
  short cur_index = modelData->splitCombineInfoIndex;

  // Move backwards to undo
  short new_index = setNewSplitCombineIndex(-1);

  SplitCombineInfo* sc_info = (*modelData->splitCombineInfos)[new_index];

  // If there is nothing to undo!
  if (sc_info == NULL || !sc_info->canUndo) {
    setNewSplitCombineIndex(+1);
    gui->showMsg("Nothing for undo!");
    return;
  }

  sc_info->canRedo = true;
  sc_info->canUndo = false;

  //---If we have to UNDO a SPLIT
  if (sc_info->splitTargetId != NO_INDEX) {

    BodyElement* be_target = getBoundaryById(sc_info->splitTargetId);
    BodyElement* be_source = getBoundaryById(sc_info->splitSourceId);

    // Remove target element from the model
    removeBodyElement(be_target, true, false);

    // Add target mesh elements (back) to the source
    for (int i = 0; i < sc_info->nofSplitTargetMeshElements; i++) {

      int elem_id = be_target->getMeshElementId(i);
      short elem_dir = be_target->getMeshElementDir(i);
      int int_index = be_source->addMeshElement(elem_id, elem_dir);
      sc_info->splitTargetMeshElementSourceIndices[i] = int_index;

      meshData->boundaryElementBoundaryIds[elem_id] = sc_info->splitSourceId;
    }

    be_source->findMeshBorder();
  }

  //---If we have to UNDO a COMBINE
  else if (sc_info->nofCombineSourceIds > 0) {

    BodyElement* be_target = getBoundaryById(sc_info->combineTargetId);

    // Remove target element from the model
    removeBodyElement(be_target, true, false);

    // Add source elements (back) to the model
    for (int i = 0; i < sc_info->nofCombineSourceIds; i++) {
      BodyElement* be_source = getRemovedBodyElement(sc_info->combineSourceIds[i]);
      restoreBodyElement(be_source, true);
      be_source->findMeshBorder();
    }

  }

  //-Check body data (rebuild element loops etc.)
  checkBodies();

  //-Update model
  updateModelStatistics();

  // Update gui and graphics screen
  gui->updateBodyData(this);
  gui->updateBoundaryData(this);

  refreshRenderer();

#if 0
  // For debugging
  strstream cur_strm;
  strstream new_strm;
  cur_strm << "cur SplitCombineIndex = " << cur_index << ends;
  new_strm << "New SplitCombineIndex = " << new_index << ends;
  gui->showMsg(cur_strm.str());
  gui->showMsg(new_strm.str());
#endif

}


void
Model::startEditMeshBoundaries()
{
  modelInfo->editingMeshBoundaries = true;

  int body1_id = modelInfo->selectedBody1Id;
  int body2_id = modelInfo->selectedBody2Id;

  if (body1_id == body2_id)
    body2_id = NO_INDEX;

  // Backup all boundaries under (possible) edit
  // We will copy these back in user quit the boundaries
  // edit panel wit CANCEL
  int& nof_boundaries = modelInfo->nofEditableMeshBoundaries;
  BodyElement**& boundaries = modelInfo->editableMeshBoundaries;

  // NOTE: We allocate here editableMeshBoundaries container!!!
  getBoundaries(body1_id, body2_id, nof_boundaries, boundaries);

  for (int i = 0; i < nof_boundaries; i++) {
    boundaries[i]->backupMeshData();
  }

}


void
Model::stopEditMeshBoundaries(bool cancel_edit)
{
  // If nothing was really edited
  if (modelInfo->editingMeshBoundaries == false)
    return;

  modelInfo->editingMeshBoundaries = false;

  modelData->splitCombineInfoIndex = 0;

  // Reset split/combine info data
  for (int i = 0; i < MAX_NOF_SPLIT_COMBINE_INFOS; i++) {
    delete (*modelData->splitCombineInfos)[i];
    (*modelData->splitCombineInfos)[i] = NULL;
  }

  // If editing was cancelled, delete current boundaries
  // and restore the original from the backup
  if (cancel_edit) {
    int i;
    int nof_cur_boundaries;
    BodyElement** cur_boundaries;

    int body1_id = modelInfo->selectedBody1Id;
    int body2_id = modelInfo->selectedBody2Id;

    if (body1_id == body2_id)
      body2_id = NO_INDEX;

    getBoundaries(body1_id, body2_id, nof_cur_boundaries, cur_boundaries);

    bool update_bodies = true;

    for (i = 0; i < nof_cur_boundaries; i++) {
      removeBodyElement(cur_boundaries[i], update_bodies, true);
    }

    // Add backup elements again into the model and restore
    // their mesh data
    for (i = 0; i < modelInfo->nofEditableMeshBoundaries; i++) {

      BodyElement* be = modelInfo->editableMeshBoundaries[i];

      addBodyElement(be, update_bodies);

      // This should not be a "removed" element anymore
      modelData->removedModelElements->erase(be->Id());

      be->restoreMeshDataBackup();
    }

    Body* body1 = getBodyById(body1_id);
    Body* body2 = getBodyById(body2_id);

    body1->check();

    if (body2 != NULL)
      body2->check();

    UserInterface* gui = theControlCenter->getGui();

    gui->updateBodyData(this);
    gui->updateBoundaryData(this);

    meshData->boundaryElements->resetSelected();
    refreshRenderer();
  } // end cancel-edit mode

  // Now it is safe to purge all "removed" elements from the
  // memory
  modelData->purgeRemovedElements();
  modelData->removedModelElements = new BodyElementTable;

  delete[] modelInfo->editableMeshBoundaries;
  modelInfo->editableMeshBoundaries = NULL;
  modelInfo->nofEditableMeshBoundaries = 0;

}


void
Model::storeBoundaryNames()
{
  int index = 0;
  while(true) {
    BodyElement* be = getBoundary(index++);
    if (be==NULL) break;
    be->storeName();
  }
}


void
Model::swapBodyElements(IdArray* ids1, IdArray* ids2, IdArray* relative_dirs)
{
  int index = 0;
  while (true) {
    BodyElementLoop* bel = getBodyElementLoop(index++);
    if (bel==NULL) break;
    bel->swapElements(ids1, ids2, relative_dirs);
  }

}


void
Model::unloadMesh(char* msg)
{
  UserInterface* gui = theControlCenter->getGui();

  if (msg == NULL || msg[0] == '\0')
    gui->showMsg("Removing mesh data from the memory");
  else
    gui->showMsg(msg);

  resetMeshData();

  gui->configureMenuButtons("File", "LoadMesh", 1);
  gui->configureMenuButtons("File", "Mesh", 1);

  setWindowTitles();
}



void
Model::updateMeshDirectoryInfo()
{
  char* dir = NULL;
  char* dir_abs = NULL;
  char* dir_head = NULL;
  char* dir_head_abs = NULL;

  UserInterface* gui = theControlCenter->getGui();

  gui->getMeshDirectoryInfo(dir, dir_abs);

  update_dyna_string(modelInfo->meshDirectory, dir);
  update_dyna_string(modelInfo->meshDirectory_absolute, dir_abs);

  delete[] dir; delete[] dir_abs;
}


void
Model::setMeshFs()
{
  delete[] modelInfo->meshFs;
  modelInfo->meshFs = NULL;

  if ( modelInfo->nofMeshes == 0 )
    return;

  modelInfo->meshFs = new double[modelInfo->nofMeshes];

  // Init
  for (int i = 0; i < modelInfo->nofMeshes; i++ )
    modelInfo->meshFs[i] = 1.0;
}


void
Model::setMeshHs()
{
  delete[] modelInfo->meshHs;
  modelInfo->meshHs = NULL;

  if ( modelInfo->nofMeshes == 0 )
    return;

  modelInfo->meshHs = new double[modelInfo->nofMeshes];


  // For version 4 files we have normally only one
  // model level parameter, take mesh-h from it!
  double mesh_h;
  Parameter* param = getParameter(0, ECIF_GRID_PARAMETER);

  if ( param != NULL && param->getFieldValueBySifName("Mesh H", mesh_h) ) {

  } else {
    mesh_h = calcInitialMeshH();
  }

  // Init with same value
  for (int i = 0; i < modelInfo->nofMeshes; i++ )
    modelInfo->meshHs[i] = mesh_h;

}


void
Model::updateModelDirectoryInfo()
{
  char* dir = NULL;
  char* dir_abs = NULL;

  UserInterface* gui = theControlCenter->getGui();

  gui->getModelDirectoryInfo(dir, dir_abs);

  update_dyna_string(modelInfo->modelDirectory, dir);
  update_dyna_string(modelInfo->modelDirectory_absolute, dir_abs);

  delete[] dir; delete[] dir_abs;
}


void
Model::updateModelNameInfo()
{
  char* model_name = NULL;
  char* problem_name = NULL;

  UserInterface* gui = theControlCenter->getGui();

  gui->getModelNameInfo(model_name, problem_name);

  update_dyna_string(modelInfo->modelName, model_name);
  update_dyna_string(modelInfo->problemName, problem_name);

  delete[] model_name;
  delete[] problem_name;
}


void
Model::updateBodyForceApplyCounts()
{
  modelStatistics->nofBodiesWithBodyForce = 0;

  int index = 0;
  while (true) {
    Body* body = getBody(index++);
    if (body==NULL) break;
    int param_id = body->getBodyForceId();
    if ( param_id != NO_INDEX ) {
      Parameter* param = getParameterById(ECIF_BODY_FORCE, param_id);
      if ( param != NULL ) {
        modelStatistics->nofBodiesWithBodyForce++;
        param->updateApplyCount(1);
      } else {
        body->setBodyForceId(NO_INDEX);
      }
    }
  }
}


void
Model::updateBoundaryConditionApplyCounts()
{
  modelStatistics->nofInnerBoundariesWithCondition = 0;
  modelStatistics->nofOuterBoundariesWithCondition = 0;
  modelStatistics->nofEdgesWithCondition = 0;
  modelStatistics->nofVerticesWithCondition = 0;

  bool param_found;
  int index = 0;
  while (true) {

    param_found = true;
    BodyElementGroup* beg = getBodyElementGroup(index++);

    if (beg==NULL) break;

    int param_id = beg->getBoundaryConditionId();

    Parameter* param = getParameterById(ECIF_BOUNDARY_CONDITION, param_id);

    if ( param == NULL ) {
      beg->setBoundaryConditionId(NO_INDEX);
      continue;
    }

    for (int i = 0; i < beg->getNofElements(); i++) {

      const BodyElement* be = beg->getElement(i);

      if ( be == NULL ) break;

      // For a virtual group we increse boundary counter
      // only if the element doe not have already an other
      // boundary condition
      //
      if ( beg->isVirtual() ) {

        BodyElementGroup* be_eg = (BodyElementGroup*)be->getElementGroup();

        if ( NO_INDEX != be_eg->getBoundaryConditionId() ) {
          param->updateApplyCount(1);
          continue;
        }
      }

      int bd1_tag = be->getParentTag(1);
      int bd2_tag = be->getParentTag(2);

      // Inner boundary
      if ( bd1_tag != NO_INDEX && bd2_tag != NO_INDEX ) {
        modelStatistics->nofInnerBoundariesWithCondition++;

      // Outer boundary
      } else if ( bd1_tag != NO_INDEX || bd2_tag != NO_INDEX ) {
        modelStatistics->nofOuterBoundariesWithCondition++;

      // Edge in 3D
      } else if ( OT_EDGE == be->getObjectType() ) {
        modelStatistics->nofEdgesWithCondition++;

      // Vertex
      } else  {
        modelStatistics->nofVerticesWithCondition++;
      }

      param->updateApplyCount(1);

    } // each group element

  } // each group
}


void
Model::updateBoundaries()
{
  setBoundaryTags();
}


void
Model::updateEquationApplyCounts()
{
  modelStatistics->nofBodiesWithEquation = 0;

  int index = 0;
  while (true) {
    Body* body = getBody(index++);
    if (body==NULL) break;
    int eq_id = body->getEquationId();

    if ( eq_id != NO_INDEX ) {
      Parameter* param = getParameterById(ECIF_EQUATION, eq_id);
      if ( param != NULL ) {
        modelStatistics->nofBodiesWithEquation++;
        param->updateApplyCount(1);
      } else {
        body->setEquationId(NO_INDEX);
      }
    }
  }
}


// Update (Matc-dependent) geometry
//
void
Model::updateCadGeometry()
{
  // Update vertex-table points
  //
  VertexTable* vt = modelData->vertexTable;
  const char* def = getMatcString(vt->matcTable, EMF_POINTS);

  if ( vt->dim1 != 0 && def != NULL ) {

    // Evaluate Matc-string for points
    strstream strm;
    char* res = LibFront::evalMatcString(def);

    if ( res != NULL ) {
      strm << res << ends;

      // Loop points
      for (int i = 0; i < vt->dim1; i++) {

        BodyElement* v = getBodyElementById(vt->vertexIds[i]);
        GcPoint* p = (GcPoint*)v->getGeometry();
        double val;

        // Loop point coordinates
        for (int j = 0; j < vt->dim2; j++) {
          strm >> val;
          p->setPos(j, val);
        }
      }
    }
  }

  int index = 0;
  while (true) {
    BodyElement* be = getBodyElement(index++, true);
    if (be==NULL) break;

    // Error!
    if ( !be->updateGeometry() ) {
      return;
    }
  }

  // Check body data (bound boxes etc.
  checkBodies();

  // Update model statistics
  updateModelStatistics();

  // Update gui
  UserInterface* gui = theControlCenter->getGui();
  gui->updateBodyData(this);
  gui->updateBoundaryData(this);

  // Remove old display lists and update geoemtry display
  //
  Renderer* renderer = theControlCenter->getRenderer();
  if ( renderer != NULL ) {
    renderer->removeDisplayLists();
    refreshRenderer();
  }
}


void
Model::updateInitialConditionApplyCounts()
{
  modelStatistics->nofBodiesWithInitialCondition = 0;

  int index = 0;
  while (true) {
    Body* body = getBody(index++);
    if (body==NULL) break;
    int param_id = body->getInitialConditionId();

    if ( param_id != NO_INDEX ) {
      Parameter* param = getParameterById(ECIF_INITIAL_CONDITION, param_id);
      if ( param != NULL ) {
        modelStatistics->nofBodiesWithInitialCondition++;
        param->updateApplyCount(1);
      } else {
        body->setInitialConditionId(NO_INDEX);
      }
    }
  }
}


bool
Model::updateMatcFile(char* filename, char* mode, int def_count, char** defs)
{
  UserInterface* gui = theControlCenter->getGui();

  // Overwrite mode
  // ==============

  // Overwrite mode, just write user definition into the file
  if ( LibFront::ncEqualPartial("w", mode ) ) {
    fstream matc_file(filename, ios::out | ios::trunc);

    for (int i = 0; i < def_count; i++) {
      matc_file << defs[i] << endl;
    }
    matc_file.close();

    return true;
  }

  // Update/append mode
  // ==================
  MatcValueTable all_defs, old_vars, old_funcs;
  NameList old_names, new_names;

  // First read all old definitions and store also variables and
  // functions in separate tables
  int comment_count = 0;
  char* source = NULL;
  char matc_buffer[10001];
  char name_buffer[256];
  bool is_var, is_func;
  LibFront::readFile(filename, source, false);

  int source_pos = 0;
  int source_len = strlen(source);

  while (source_pos < source_len ) {
    source_pos = LibFront::readMatcExpression(source, source_pos, source_len,
                                              matc_buffer, 10000,
                                              false);
    LibFront::trim(matc_buffer);

    if ( matc_buffer[0] == '\0' ) continue;

    LibFront::getMatcName(matc_buffer, name_buffer, 255, is_var, is_func);

    // Store old name
    // NOTE: Allocate new array for the name!
    char* name;

    if ( name_buffer != NULL && name_buffer[0] != '\0' ) {
      create_dyna_string(name, name_buffer);
    } else {
      strstream strm;
      strm << "comment " << ++comment_count << ends;
      create_dyna_string(name, strm.str());
    }
    old_names.push_back(name);

    // Store value
    // NOTE: Allocate new array for the value!
    char* value;
    create_dyna_string(value, matc_buffer);
    all_defs[name] =value;

    if ( is_var ) {
      old_vars[name] = "1";
    } else if ( is_var )  {
      old_funcs[name] = "1";
    }

  }

  // Next update old-defs by new defs
  for (int i = 0; i < def_count; i++) {
    LibFront::getMatcName(defs[i], name_buffer, 255, is_var, is_func);

    // A new variable cannot be an old function!
    if ( is_var && NULL != old_funcs[name_buffer] ) {
      strstream strm;
      strm << "Cannot update the Matc file!" << "\\\\n"
           << "variable: " << name_buffer << "\\\\n"
           << "already defined as a function in the file!" << ends;
      gui->errMsg(0, strm.str());
      return false;

    // A new function cannot be an old variable!
    } else if ( is_func && NULL != old_vars[name_buffer] ) {
      strstream strm;
      strm << "Cannot update the Matc file!" << "\\\\n"
           << "function: " << name_buffer << "\\\\n"
           << "already defined as a variable in the file!" << ends;
      gui->errMsg(0, strm.str());
      return false;
    }

    // NOTE: Allocate new array for the name!
    char* name;
    create_dyna_string(name, name_buffer);

    // New object: store name
    if ( NULL == all_defs[name_buffer] ) {
      new_names.push_back(name);

    // Old object: delete old value
    } else {
      delete[] all_defs[name];
    }

    // Ok, update old or add new
    // NOTE: Allocate new array for the value
    char* value;
    create_dyna_string(value, defs[i]);
    all_defs[name] = value;
  }

  // Then write all definitions to the file
  // NOTE: If a 'definition' is just the linefeed, it is
  // not actually written to file,  just lf!
  //
  fstream matc_file(filename, ios::out | ios::trunc);

  NameList::iterator itr;

  itr = old_names.begin();
  while ( itr != old_names.end() ) {
    char* nm = (*itr);
    char* def = all_defs[nm];
    if ( strlen(def) == 1 && def[0] == '\n' ) {
      matc_file << def;
    } else {
      matc_file << def << endl;
    }
    itr++;
  }

  itr = new_names.begin();
  while ( itr != new_names.end() ) {
    char* nm = (*itr);
    char* def = all_defs[nm];
    if ( strlen(def) == 1 && def[0] == '\n' ) {
      matc_file << def;
    } else {
      matc_file << def << endl;
    }
    itr++;
  }

  matc_file.close();

  purgeMatcValueTable(all_defs);
  purgeNameList(old_names);
  purgeNameList(new_names);

  return true;
}


void
Model::updateMaterialApplyCounts()
{
  modelStatistics->nofBodiesWithMaterial = 0;

  int index = 0;
  while (true) {
    Body* body = getBody(index++);
    if (body==NULL) break;
    int param_id = body->getMaterialId();

    if ( param_id != NO_INDEX ) {
      Parameter* param = getParameterById(ECIF_MATERIAL, param_id);
      if ( param != NULL ) {
        modelStatistics->nofBodiesWithMaterial++;
        param->updateApplyCount(1);
      } else {
        body->setMaterialId(NO_INDEX);
      }
    }
  }
}


void
Model::updateMinimumEdgeSize(int nof_points, GcPoint** points)
{
  for (int i = 1 ; i < nof_points; i++) {

    Point3* p1 = points[i - 1]->getPoint();
    Point3* p2 = points[i]->getPoint();

    double dist = dist3(*p1, *p2);

    if (dist < modelInfo->minEdgeSize)
      modelInfo->minEdgeSize = dist;
  }
}


void
Model::updateParametersApplyCounts(ecif_parameterType parameter_type)
{
  switch (parameter_type) {
  case ECIF_BODY_FORCE:
    updateBodyForceApplyCounts();
    break;
  case ECIF_BOUNDARY_CONDITION:
    updateBoundaryConditionApplyCounts();
    break;
  case ECIF_EQUATION:
    updateEquationApplyCounts();
    break;
  case ECIF_INITIAL_CONDITION:
    updateInitialConditionApplyCounts();
    break;
  case ECIF_MATERIAL:
    updateMaterialApplyCounts();
    break;
  }
}


// Update parent object info which is read from a emf-file
// Needed ex. for BoundaryConditions, whose parents are mainly
// normal boundaries, but in the model parents are bounadry groups
// (logical or real)
//
void
Model::updateParametersParentId()
{
  updateParametersParentId(ECIF_INITIAL_CONDITION);
  updateParametersParentId(ECIF_BOUNDARY_CONDITION);
  updateParametersParentId(ECIF_MATERIAL);
  updateParametersParentId(ECIF_BODY_FORCE);
  updateParametersParentId(ECIF_GRID_PARAMETER);
  updateParametersParentId(ECIF_GRID_H);
}


void
Model::updateParametersParentId(ecif_parameterType param_type)
{
  int index = 0;
  while (true) {
    Parameter* param = getParameter(index++, param_type);
    if (param==NULL) break;
    param->updateParentId();
  }

}


// NOTE: Used for <=4 version emf-files
//
void
Model::updateParametersParentTags()
{
  Parameter* param;
  int pid;
  int ctag;

  int index = 0;
  while (true) {

    Body* body = getBody(index++);
    if (body==NULL) break;

    pid = body->getBodyForceId();
    param = getParameterById(ECIF_BODY_FORCE, pid);
    if ( param != NULL && NO_INDEX == param->getParentEmfTag() ) {
      param->setParentEmfTag(body->Tag());
      param->setParentEmfType(OT_BODY);
    }

    pid = body->getInitialConditionId();
    param = getParameterById(ECIF_INITIAL_CONDITION, pid);
    if ( param != NULL && NO_INDEX == param->getParentEmfTag() ) {
      param->setParentEmfTag(body->Tag());
      param->setParentEmfType(OT_BODY);
    }

    pid = body->getMaterialId();
    param = getParameterById(ECIF_MATERIAL, pid);
    if ( param != NULL && NO_INDEX == param->getParentEmfTag() ) {
      param->setParentEmfTag(body->Tag());
      param->setParentEmfType(OT_BODY);
    }
  }

  index = 0;
  while (true) {

    BodyElementGroup* beg = getBodyElementGroup(index++);
    if (beg==NULL) break;

    if ( EXPLICIT_GROUP != beg->getGroupType() ) continue;

    pid = beg->getBoundaryConditionId();
    param = getParameterById(ECIF_BOUNDARY_CONDITION, pid);
    if ( param != NULL &&
         ( NO_INDEX == param->getParentEmfTag() ||
           OT_BODY == param->getParentEmfType()
         )
         ) {
      param->setParentEmfTag(beg->Tag());
      param->setParentEmfType(beg->getObjectType());
    }
  }

  index = 0;
  while (true) {

    BodyElement* be = getBodyElement(index++);
    if (be==NULL) break;

    pid = be->getBoundaryConditionId();
    param = getParameterById(ECIF_BOUNDARY_CONDITION, pid);
    if ( param != NULL &&
         ( NO_INDEX == param->getParentEmfTag() ||
           OT_BODY == param->getParentEmfType()
         )
         ) {
      param->setParentEmfTag(be->Tag());
      param->setParentEmfType(be->getObjectType());
    }
  }
}


void
Model::variableNameGuiToSif(const char* gui_name, char* sif_name_buffer)
{
  theControlCenter->getUI()->variableNameGuiToSif(gui_name, sif_name_buffer);
}


void
Model::variableNameSifToGui(const char* sif_name, char* gui_name_buffer)
{
  theControlCenter->getUI()->variableNameSifToGui(sif_name, gui_name_buffer);
}



// **********************
//   STL etc. utilities
// **********************

int
find2(Ids2Set& id_set, int key1)
{
  Ids2Set::iterator pos;
  pos = id_set.begin();
  while (pos != id_set.end()) {
    Ids2 ids = *pos++;
    if (ids.id1 < key1)
      continue;
    else if (ids.id1 == key1)
      return ids.id2;
    else
      return NO_INDEX;
  }

  return NO_INDEX;
}


int
find3(Ids3Set& id_set, int key1, int key2)
{
  Ids3Set::iterator pos;
  pos = id_set.begin();
  while (pos != id_set.end()) {
    Ids3 ids = *pos++;
    if (ids.id1 < key1)
      continue;
    else if (ids.id1 == key1 && ids.id2 < key2 )
      continue;
    else if (ids.id1 == key1 && ids.id2 == key2)
      return ids.id3;
    else
      return NO_INDEX;
  }

  return NO_INDEX;
}


// Get Matc-expression string using a name-key
// NOTE: An empty definition is returned as a NULL-string
//
const char*
getMatcString(MatcValueTable& table, const char* key)
{
  char* val = table[key];

  if ( val == NULL || val[0] == '\0' ) {
    return NULL;
  } else {
    return val;
  }
}


// Store Matc-expression string using a name-key
// NOTE: An empty definition will nullify a possible older
// definition
//
void
storeMatcString(MatcValueTable& table, const char* key, const char* value)
{
  // Delete a possible old definition
  char* old_val = table[key];
  if ( old_val != NULL ) {
    delete[] old_val;
    table.erase(key);
  }

  // Create new definition string if a non-empty
  // matc-definition given
  //
  char* str = NULL;

  // If non-empty value given
  if ( value != NULL && strlen(value) > 1 ) {

    char* new_val;
    const char* work = value;
    
    // Remove possible leading $-sign
    //
    if ( work[0] == '$' ) {
      work++;
    }

    // Allocate str and copy tmp
    create_dyna_string(new_val, work);

    // Store key/value pair
    table[key] = new_val;
  }
}


// Copy any Matc-expression string from one table to an other
//
void
copyMatcValueTable(MatcValueTable& source, MatcValueTable& target)
{
  MatcValueTable::iterator pos = source.begin();
  MatcValueTable::iterator end_pos = source.end();

  while ( pos != end_pos ) {

    // Pick value
    const char* key = (*pos).first.c_str();
    char* val = (*pos).second;

    // Delete possible old value
    char* old_val = target[key];
    if ( old_val != NULL ) {
      delete[] old_val;
      target.erase(key);
    }

    // Create new value
    char* new_val;
    create_dyna_string(new_val, val);

    // Store new value using the key
    target[key] = new_val;

    pos++;
  }
}


// Delete values from the Matc value table and clear it
//
// NOTE: Use this only if values (char*) are allocated with new!!!
//
void
purgeMatcValueTable(MatcValueTable& table)
{
  MatcValueTable::iterator pos = table.begin();
  MatcValueTable::iterator end_pos = table.end();

  while (pos != end_pos) {
    delete[] (*pos).second;
    pos++;
  }
  table.clear();
}


// Delete values from a name list and clear it
//
// NOTE: Use this only if values (char*) are allocated with new!!!
//
void
purgeNameList(NameList& list)
{
  NameList::iterator pos = list.begin();
  NameList::iterator end_pos = list.end();

  while (pos != end_pos) {
    delete[] (*pos);
    pos++;
  }
  list.clear();
}


// Delete values from a name table and clear it
//
// NOTE: Use this only if values (char*) are allocated with new!!!
//
void
purgeNameTable(NameTable& table)
{
  NameTable::iterator pos = table.begin();
  NameTable::iterator end_pos = table.end();

  while (pos != end_pos) {
    delete[] (*pos).second;
    pos++;
  }
  table.clear();
}



// Array reallocate functions
// ==========================

void
reallocate_array(int old_size, int new_size, int*& array, int default_value)
{
  reallocate_array_impl1(old_size, new_size, array, default_value);
}

void
reallocate_array(int old_size, int new_size, char**& array, char* default_value)
{
  reallocate_array_impl2(old_size, new_size, array, default_value);
}

void
reallocate_array(int old_size, int new_size, int**& array, int* default_value)
{
  reallocate_array_impl2(old_size, new_size, array, default_value);
}

void
reallocate_array(int old_size, int new_size, BoundBox**& array)
{
  reallocate_array_impl3(old_size, new_size, array);
}

void
reallocate_array(int old_size, int new_size, BodyElementTable**& array)
{
  reallocate_array_impl3(old_size, new_size, array);
}

void
reallocate_array(int old_size, int new_size, IdList**& array)
{
  reallocate_array_impl3(old_size, new_size, array);
}

void
reallocate_array(int old_size, int new_size, IdArray**& array)
{
  reallocate_array_impl3(old_size, new_size, array);
}


// Implmentation of the template functions
// ---------------------------------------


// Function reallocates an object array and sets a default value for
// new entries
//
template <class T > void
reallocate_array_impl1(int old_size, int new_size, T*& array, T default_value)
{
  int i;
  T* new_array = new T[new_size];

  for (i = 0; i < old_size; i++) {
    new_array[i] = array[i];
  }

  for (i = old_size; i < new_size; i++) {
    new_array[i] = default_value;
  }

  delete[] array;

  array = new_array;
}


// Function reallocates an object array and sets a pointer default value for
// new entries
//
template <class T > void
reallocate_array_impl2(int old_size, int new_size, T**& array, T* default_value)
{
  int i;
  T** new_array = new T*[new_size];

  for (i = 0; i < old_size; i++) {
    new_array[i] = array[i];
  }

  for (i = old_size; i < new_size; i++) {
    new_array[i] = default_value;
  }

  delete[] array;

  array = new_array;
}


// Function reallocates an object array and allocates objects for
// new entries
//
template <class T > void
reallocate_array_impl3(int old_size, int new_size, T**& array)
{
  int i;
  T** new_array = new T*[new_size];

  for (i = 0; i < old_size; i++) {
    new_array[i] = array[i];
  }

  for (i = old_size; i < new_size; i++) {
    new_array[i] = new T;
  }

  delete[] array;

  array = new_array;
}
