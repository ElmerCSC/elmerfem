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
Module:     ecif_body.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_body.h"
#include "ecif_bodyLayer.h"
#include "ecif_bodyElement.h"
#include "ecif_bodyElement2D.h"
#include "ecif_bodyElement3D.h"
#include "ecif_bodyElementGroup.h"
#include "ecif_bodyElementLoop.h"
#include "ecif_mesh.h"
#include "ecif_model.h"
#include "ecif_parameter.h"
#include "ecif_parameterField.h"
#include "ecif_renderer.h"
#include "ecif_userinterface.h"

//Initialize static class variables.
int BodyPair::last_tag = 0;
int Body::last_tag = 0;


// BodyPair class
// ==============

BodyPair::BodyPair()
{
  id = NO_INDEX;
  tag = NO_INDEX;
  init();
}


BodyPair::BodyPair( const Body* bd1, const Body* bd2)
{
  id = NO_INDEX;
  tag = NO_INDEX;
  init();

  body1 = bd1;
  body2 = bd2;
}



BodyPair::BodyPair(int int_tag, const Body* bd1, const Body* bd2)
{
  tag = int_tag;

  init();

  body1 = bd1;
  body2 = bd2;

  if (last_tag < tag)
    last_tag = tag;
}


BodyPair::~BodyPair()
{
}


const char*
BodyPair::getName()
{
  strstream strm;
  strm << ((Body*)body1)->getName() << "-" << ((Body*)body2)->getName() << ends;
  return strm.str();
}


void
BodyPair::init()
{
  body1 = NULL;
  body2 = NULL;
}


void
BodyPair::initClass(Model* mdl)
{
  BodyPair::model = mdl;
  BodyPair::last_tag = 0;
}


int
BodyPair::setTag()
{
  tag = ++last_tag;

  return tag;
}


int
BodyPair::setId()
{
  model->addModelObject(this, OT_BODYPAIR);

  return Id();
}




// Body class
// ==========

Body::Body()
{
  tag = ++last_tag;
  init();
}


Body::Body(bodyGmtrType gmtr_type, int ext_tag, char* body_name, colorIndices color_index)
{
  tag = ++last_tag;

  init();

  colorIndex = color_index;
  for (int i = 0; i < 4; i++)
    color[i] = colorValues[colorIndex][i];

  externalTag = ext_tag;

  update_dyna_string(name, body_name);

  gmtrType = gmtr_type;
}



Body::Body(bodyGmtrType gmtr_type, int int_tag, int ext_tag, char* body_name, colorIndices color_index)
{
  tag = int_tag;

  if (last_tag < tag) {
    last_tag = tag;
  }

  init();

  colorIndex = color_index;
  for (int i = 0; i < 4; i++)
    color[i] = colorValues[colorIndex][i];

  externalTag = ext_tag;

  update_dyna_string(name, body_name);

  gmtrType = gmtr_type;
}


Body::Body(ecif_Body_X& tx, bool add_default_layer)
{
  tag = tx.tag;

  if (last_tag < tag) {
    last_tag = tag;
  }

  init(add_default_layer);

  int i;

  if ( tx.has_color ) {
    colorIndex = ef_nodefault;
    for (i = 0; i < 4; i++)
      color[i] = tx.color[i];

  } else {
    colorIndex = DEFAULT_COLOR_INDEX;
  }

  update_dyna_string(name, tx.name);

  bodyForceId = tx.body_force_id;
  bodyParameterId = tx.body_param_id;
  equationId = tx.equation_id;
  initialConditionId = tx.init_cond_id;
  materialId = tx.material_id;

  if (tx.is_checked) {
    status = true;
    checked = true;
  }


  if ( tx.is_open ) {
    tplgType = OPEN_BODY;
  }

  if ( tx.is_bem ) {
    type = BEM_BODY;

  } else if ( tx.is_virtual ) {
    type = VIRTUAL_BODY;
    drawMode = DM_HIDDEN;
  }


}


Body::~Body()
{
  removeLayers();

  removeElementGroups();

  delete boundbox;
  delete boundboxMesh;

  delete[] meshElementIds;
}


// Basic version which just creates the loop of elements without
// checking their orientation or connectivity.
// Used mainly for mesh-based boundaries for which checking is not relevant
// Return true if body is ok, otherwise false.
//
bool
Body::check()
{
  // NOTE: Virtual bodies are checked separately, because their
  // possible pending edges must be inserted into a virtual group!!
  //
  if ( isVirtual() ) {
    return checkVirtual();
  }

  if (nofLayers == 0) return false;

  initName();

  int layer = -1;
  while (true) {

    if (!selectLayer(++layer)) break;

    addAllPendingVertices(layer);
    addAllPendingElements(layer);

    checked = true;

    nofElements[layer] = belements[layer]->size();

    // We do not have yet any elements!
    if (nofElements[layer] == 0 ) {

      // If even no loops, something is really wrong!
      if (nofElementLoops[layer] == 0) {
        status = false;
        return false;

      // Element ids should be in the loops
      } else if ( !addAllLoopElements(layer) ) {
        status = false;
        return false;
      }
    }

    nofElements[layer] = belements[layer]->size();

    // If still no elements, something is wrong!!!
    if (nofElements[layer] == 0 ) {
      status = false;
      return false;
    }

    // NOTE: Currently we don't check the orientation of real loops (like
    // a cube inside a cube !!!***!!!
    // we just put all existing elements to ONE SINGLE LOOP !!!
    //
    // So this is meant for: mesh bodies, open bodies

    // Create (one) new loop from all elements
    // =======================================
    status = false;
    removeElementLoops(layer);

    int nofFreeElements = nofElements[layer];

    while (nofFreeElements > 0) {

      // Create one element loop
      int be_index = 0;

      IdList elementLoop;
      BodyElement* be = getElement(layer, be_index++);

      while (be) {
        nofFreeElements--;
        elementLoop.push_back(be->Id());
        be = getElement(layer, be_index++);
      }

      BodyElementLoop* bel = new BodyElementLoop(&elementLoop, false, OT_BOUNDARY);

      // Store loop in the model
      model->addBodyElementLoop(bel);

      // Store loop info also in the body
      int bid = bel->Id();
      elementLoopIds[layer]->push_back(bel->Id());
      elementLoopTags[layer]->push_back(bel->Tag());
      nofElementLoops[layer]++;
    }

  } // All layers


  status = true;
  return true;
}


bool
Body::checkLayerIndex(int layer)
{
  if ( layer < 0 || layer >= nofLayers ) return false;
  return true;
}


// Basic version which just creates the loop of elements without
// checking their orientation or connectivity.
// Used mainly for mesh-based boundaries for which checking is not relevant
// Return true if body is ok, otherwise false.
//
bool
Body::checkVirtual()
{
  if (nofLayers == 0) return false;

  initName();

  int layer = -1;
  while (true) {

    if (!selectLayer(++layer)) break;

    checked = true;

    addAllPendingVertices(layer);
    addAllPendingElements(layer);

    nofElements[layer] = belements[layer]->size();

    // We do not have yet any elements!
    if (nofElements[layer] == 0 ) {

      // If even no loops, something is really wrong!
      if (nofElementLoops[layer] == 0) {
        status = false;
        return false;

      // Element ids should be in the loops
      } else if ( !addAllLoopElements(layer) ) {
        status = false;
        return false;
      }
    }

    nofElements[layer] = belements[layer]->size();

    // If still no elements, something is wrong!!!
    if (nofElements[layer] == 0 ) {
      status = false;
      return false;
    }

    int be_count = 0;
    int* be_ids = NULL;
    int* be_tgs = NULL;

    getElementIds(layer, be_count, be_ids);
    getElementTags(layer, be_count, be_tgs);

    // Group eleemnt type by th efirst element!
    BodyElement* be = getElement(layer, 0);

    BodyElementGroup* beg = new BodyElementGroup();
    beg->setElementIdsAndTags(be_count, be_ids, be_tgs);
    beg->setType(VIRTUAL_GROUP);
    beg->setElementType(be->getObjectType());
    beg->setParentId(1, id);
    beg->setParentTag(1, tag);
    beg->setName((char*)getName(layer));

    delete[] be_ids;
    delete[] be_tgs;

  } // All layers


  status = true;
  return true;
}


// NOTE: For a virtual body only!
//
bool
Body:: addElementGroupTag(int group_tag)
{
  nofElementGroups += 1;

  elementGroupTags->push_back(group_tag);
  elementGroupIds->push_back(NO_INDEX);

  return true;
}


// NOTE: For a virtual body only!
//
bool
Body:: addElementGroupTags(int nof_groups, int* group_tags)
{
  nofElementGroups += nof_groups;

  for(int i = 0; i < nof_groups; i++) {
    elementGroupTags->push_back(group_tags[i]);
    elementGroupIds->push_back(NO_INDEX);
  }

  return true;
}


bool
Body:: addElementLoopTags(int layer, int nof_loops, int* loop_tags)
{
  if (!checkLayerIndex(layer)) return false;

  nofElementLoops[layer] += nof_loops;

  for(int i = 0; i < nof_loops; i++) {
    elementLoopTags[layer]->push_back(loop_tags[i]);
    elementLoopIds[layer]->push_back(NO_INDEX);
  }

  return true;
}



// Add all elements from all bodyelement loops in the body
// NOTE: This should be used only when elements are not
// added directly to the body.
// Returns true if success
bool
Body::addAllLoopElements(int layer)
{
  if (!checkLayerIndex(layer)) return false;

  UserInterface* gui = (UserInterface*)model->getGui();

  IdList::iterator id_beg = elementLoopIds[layer]->begin();
  IdList::iterator id_end = elementLoopIds[layer]->end();

  IdList::iterator tg_beg = elementLoopTags[layer]->begin();
  IdList::iterator tg_end = elementLoopTags[layer]->end();

  IdList::iterator id_itr;
  IdList::iterator tg_itr;

  DirectedBodyElement d_be;

  bool is_first_loop = true;

  // All body element loops
  for (id_itr=id_beg,tg_itr=tg_beg; id_itr != id_end; id_itr++,tg_itr++) {

    int bel_id = *id_itr;
    int bel_tg = *tg_itr;

    bel_id = (bel_id < 0) ? -bel_id : bel_id;

    BodyElementLoop* bel = model->getBodyElementLoopById(bel_id);

    if ( bel == NULL ) {

      strstream strm1, strm2;
      strm1 << "***ERROR in geometry for the body " << tag << ":" << ends;
      strm2 << "---Element loop " << bel_tg << " not defined in the input!" << ends;
      gui->showMsg(strm1.str());
      gui->showMsg(strm2.str(), 1);

      return false;
    }


    // All elements in the loop
    int index =0;
    while (true) {

      bel->getElement(index++, d_be);
      BodyElement* be = d_be.element;

      if (be==NULL) break;

      bool add = true;

      // If element already in the body
      if ( NULL != getElementById(layer, be->Id()) ) {

        // Outermost loop cannot contain multiple elements
        if ( is_first_loop ) {
          strstream strm1, strm2;
          strm1 << "***ERROR in geometry for the body " << tag << ":" << ends;
          strm2 << "---Element " << be->Tag() << " included more than once!" << ends;
          gui->showMsg(strm1.str());
          gui->showMsg(strm2.str(), 1);
          return false;

        // Inner loops could introduce same elements many times, but these
        // must be inner boundaries and they are not accepted
        } else {
          removeElement(layer, be);
          add = false;
        }
      }

      if ( add ) {
        addElement(layer, be);
      }

    } // all elements in the loop

    is_first_loop = false;

  } // all loops

  return true;
}


// Add all pending element ids
// Returns nof elements added
int
Body::addAllPendingElements(int layer)
{
  if (!checkLayerIndex(layer)) return false;

  int nof_elements = 0;

  for (int i = 0; i < pendingElementTags[layer]->size(); i++) {

    // Read and remove last element
    int be_tag = (*pendingElementTags[layer])[i];

    BodyElement* be = model->getBoundaryByTag(be_tag);

    if (be != NULL) {
      addElement(layer, be);
      nof_elements++;
    }
  }

  pendingElementTags[layer]->clear();

  return nof_elements;
}


// NOTE: Take care that bodies are
// checked (method check()) after calling
// this method (before body is really used) !!!
//
// Otherwise body's state is not ok. Especially
// element is NOT yet inserted to any bodyElement
// loop here!!!
//
int
Body::addElement(int layer, BodyElement* be)
{
  if (!checkLayerIndex(layer)) return 0;

  BodyElementTable* elements = belements[layer];

  RangeVector rv;

  (*elements)[be->Id()] = be;

  nofElements[layer]++;

  if ( !be->getRangeVector(rv) )
    return 0;

  updateBoundBox(layer, rv);

  return 0;
}


// "Add" layer by allocating more space to layer-related arries
int
Body::addLayer(ecif_BodyLayer_X& tx)
{
  int old_count = nofLayers++;

  reallocate_array(old_count, nofLayers, layerIds, NO_INDEX);
  reallocate_array(old_count, nofLayers, boundboxes);
  reallocate_array(old_count, nofLayers, belements);
  reallocate_array(old_count, nofLayers, pendingElementTags);
  reallocate_array(old_count, nofLayers, pendingVertexGroups);
  reallocate_array(old_count, nofLayers, pendingVertexTags);
  reallocate_array(old_count, nofLayers, elementLoopIds);
  reallocate_array(old_count, nofLayers, elementLoopTags);
  reallocate_array(old_count, nofLayers, nofElements, 0);
  reallocate_array(old_count, nofLayers, nofElementLoops, 0);
  reallocate_array(old_count, nofLayers, nofInnerBoundaries, 0);
  reallocate_array(old_count, nofLayers, nofOuterBoundaries, 0);

  // Create new model object layer
  BodyLayer* lr = new BodyLayer(tx);

  // If body is declared open, all layers must
  // be open!
  //
  if ( isOpen() ) {
    lr->setTplgType(OPEN_LAYER);
  }

  layerIds[nofLayers-1] = lr->Id();

  return nofLayers;
}


int
Body::addMeshElements(int nof_new_elements, int* new_element_ids)
{
  int new_size = nofMeshElements + nof_new_elements;

  // Allocate new table
  int* new_ids = new int[new_size];

  // Copy old ids
  for (int i = 0; i < nofMeshElements; i++)
    new_ids[i] = meshElementIds[i];

  // Add new ids
  for (int j = 0; j < nof_new_elements; j++)
    new_ids[nofMeshElements + j] = new_element_ids[j];

  // Update info
  nofMeshElements = new_size;
  delete[] meshElementIds;
  meshElementIds = new_ids;

  return nofMeshElements;
}


// NOTE: Take care that pending elements tags are
// finally added as element ids to the body
//
int
Body::addPendingElement(int layer, int element_tag)
{
  if (!checkLayerIndex(layer)) return 0;

  pendingElementTags[layer]->push_back(element_tag);

  return 0;
}


// NOTE: Take care that pending vertices are
// finally added as an element to the body
//
int
Body::addPendingVertex(int layer, int layer_id, int vertex_tag)
{
  if (!checkLayerIndex(layer)) return 0;

  pendingVertexGroups[layer]->push_back(layer_id);
  pendingVertexTags[layer]->push_back(vertex_tag);

  return 0;
}


ecif_modelStatus
Body::checkStatus()
{
  ecif_modelStatus body_status = STATUS_OK;

  int equationId = getEquationId();
  int materialId = getMaterialId();

  if (equationId == NO_INDEX) {
    body_status |= BODY_EQUATION_MISSING;
  }

  if (materialId == NO_INDEX) {
    body_status |= BODY_MATERIAL_MISSING;
  } else {
    Parameter* material = model->getParameterById(ECIF_MATERIAL, materialId);
    if (material != NULL) {
      body_status |= material->checkStatus();
    }
  }
  return body_status;
}


bool
Body::convertTags2Ids()
{
  UserInterface* gui = (UserInterface*)model->getGui();
  strstream strm1, strm2;

  int i,j;

  IdList::iterator id_itr;
  IdList::iterator tg_itr;

  // Layer level data
  //
  for (i = 0; i < nofLayers; i++) {

    // Element loop tags to ids
    id_itr = (elementLoopIds[i])->begin();
    tg_itr = (elementLoopTags[i])->begin();

    for (j = 0; j < nofElementLoops[i]; j++, id_itr++,tg_itr++) {

      int bel_tag = *tg_itr;

      int sign = (bel_tag < 0)?-1:1;

      BodyElementLoop* bel = model->getBodyElementLoopByTag(sign * bel_tag);

      // Tags are stored unsigned!
      *tg_itr = sign * bel_tag;

      // Ids are stored signed!
      if ( bel != NULL ) {
        *id_itr = sign * bel->Id();
      } else {
        *id_itr = NO_INDEX;
        strm1 << "***ERROR in geometry for the body " << tag << ":" << ends;
        strm2 << "---Edge loop " << bel_tag << " not defined in the input!" << ends;
        gui->showMsg(strm1.str());
        gui->showMsg(strm2.str(), 1);
        return false;
      }
    }

  } // all layers


  // Element grops (virtual body stuff)
  //
  id_itr = elementGroupIds->begin();
  tg_itr = elementGroupTags->begin();

  for (i = 0; i < nofElementGroups; i++, id_itr++,tg_itr++) {

    int eg_tag = *tg_itr;

    BodyElementGroup* eg = model->getBodyElementGroupByTag(eg_tag);

    if ( eg != NULL ) {
      *id_itr = eg->Id();

      // ERROR if this is not a virtual group
      if ( !eg->isVirtual() ) {
        strm1 << "***ERROR in geometry for the virtual body " << tag << ":" << ends;
        strm2 << "---Edge group " << eg_tag << " should be virtual (Type=Virtual)!" << ends;
        gui->showMsg(strm1.str());
        gui->showMsg(strm2.str(), 1);
        return false;
      }

      // Set parent body info
      eg->setParentId(1, id);
      eg->setParentTag(1, tag);

    } else {
      *id_itr = NO_INDEX;
      strm1 << "***ERROR in geometry for Body " << tag << ":" << ends;
      strm2 << "---Edge group " << eg_tag << " not defined in the input!" << ends;
      gui->showMsg(strm1.str());
      gui->showMsg(strm2.str(), 1);
      return false;
    }
  } // all element groups

  return true;
}


void
Body::deleteMeshElements()
{
  delete[] meshElementIds;
  meshElementIds = NULL;
  nofMeshElements = 0;
}


void
Body::calcMeshBox()
{
  MeshElementTable* bulks = model->getMeshBulkElements();
  Point3* nodes = model->getMeshNodeData();

  RangeVector rv;
  rv[0] = rv[2] = rv[4] = MAX_RANGE;
  rv[1] = rv[3] = rv[5] = MIN_RANGE;

  // Loop all mesh elements
  for (int i = 0; i < nofMeshElements; i++) {
    int elem_id = meshElementIds[i];

    // Pick mesh element nodes
    meshElementCode elem_code = bulks->getElementCode(elem_id);
    const int* node_ids = bulks->getNodeIds(elem_id);
    short nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];

    // Loop all nodes in the mesh element
    for (short j = 0; j < nof_nodes; j++) {
      int node_id = node_ids[j];
      double* point = nodes[node_id];
      double x = point[0];
      double y = point[1];
      double z = point[2];
      if (x < rv[0])
        rv[0] = x;
      if (x > rv[1])
        rv[1] = x;
      if (y < rv[2])
        rv[2] = y;
      if (y > rv[3])
        rv[3] = y;
      if (z < rv[4])
        rv[4] = z;
      if (z > rv[5])
        rv[5] = z;
    }
  }

  // If no cad geometry, use mesh dimensions
  // body's boundbox
  if (gmtrType == MESH_BODY) {
    boundbox->extendByRange(rv);
  }

  boundboxMesh->extendByRange(rv);
}


bool
Body::canHaveZeroVelocityElements()
{
  int equationId = getEquationId();

  if ( equationId == NO_INDEX ) {
    return false;
  }

  Parameter* eq = model->getParameterById(ECIF_EQUATION, equationId);

  if ( eq == NULL ) {
    return false;
  }

  ParameterField* pf = eq->getFieldBySifName("Navier-Stokes");

  if ( pf != NULL ) {
    char** data = pf->getDataStrings();

    if ( LibFront::ncEqual((char*)data[0], "true") ) {
      return true;
    }
  }

  return false;
}


bool
Body::checkElements()
{
  BodyElementList devided_elements;
  BodyElementList covering_elements;
  IdList covering_ids;

  int layer = 0;

  if (!selectLayer(layer)) return false;

  while (true) {

    // First check outer boundaries which
    // have "hanging" components to be created
    int be_index = 0;
    BodyElement* be = getElement(layer, be_index++);

    while (be) {

      if ( !(be->getStatus() & BE_DEVIDED) ) {
        be = getElement(layer, be_index++ );
        continue;
      }

      devided_elements.push_back(be);

      covering_ids.clear();

      if ( !model->getCoveringElementList(be->Id(), covering_ids) ) {
        be = getElement(layer, be_index++);
        continue;
      }

      while (!covering_ids.empty()) {

        int ce_id = covering_ids.front();
        covering_ids.pop_front();

        ce_id = ( ce_id < 0)?-ce_id:ce_id;

        BodyElement* ce = model->getBodyElementById(ce_id);

        if ( ce != NULL ) {
          covering_elements.push_back(ce);
        }
      }

      be = getElement(layer, be_index++);
    }

    while (!devided_elements.empty()) {
      removeElement(layer, devided_elements.front() );
      devided_elements.pop_front();
    }

    while ( !covering_elements.empty() ) {
      addElement(layer, covering_elements.front() );
      covering_elements.pop_front();
    }

    if (!selectLayer(++layer)) break;
  }

  return true;
}


void
Body::draw(Renderer* renderer, flagName geometry_type, objectDrawingMode dmode)
{
  if ( isVirtual() ) return;

  int i, idx;

  // Check if this body is possibly excluded from meshing
  int layer = 0;

  if (!selectLayer(layer)) return;

  bool mesh_excluded = isExcludedFromCurrentMesh(layer);

  // Draw bulk elements' mesh
  if ( geometry_type == DRAW_SOURCE_MESH &&
       model->getFlagValue(DRAW_TARGET_BODIES) &&
       !mesh_excluded
     ) {
    drawMesh(renderer);
    return;
  }

  renderer->name_save(id);

  // Draw all bodyelements of the body in counter-clock-wise order.
  // This is how OpenGL thinks.
  // NOTE: originally elements are thought to be in ccw
  // If body is stored in ccw-order, then *loopDirection* = 1 (cw = -1).
  // So the way we read *edgeLoop* depends on attribute *loopDirection*.
  // If it is negative, then bodyelements are in cw-direction in the elementloop
  // and we must start from the end of loop.    !!!***!!!
  int loopDirection;
  DirectedBodyElement d_be;

  for (layer = 0; layer <  nofLayers; layer++) {

    IdList::iterator beg = (elementLoopIds[layer])->begin();
    IdList::iterator end = (elementLoopIds[layer])->end();

    BodyLayer* lr = getLayer(layer);

    if ( lr == NULL ) continue;

    int lid = lr->Id();

    renderer->name_save(lid);

    int disp_list_id = OBJECT_DISPLAY_LIST_BASE + lid;

    bool creating_display_list = false;

    bool drawing_cad_body = geometry_type == DRAW_SOURCE_CAD &&
                            !lr->isOpen() &&
                            dmode == DM_NORMAL;

    // Set display list stuff
    //
    if ( drawing_cad_body ) {

      if ( renderer->hasDisplayList(disp_list_id) ) {
        renderer->useDisplayList(disp_list_id);
        //renderer->name_delete(id);
        //return;

      } else {
        renderer->startDisplayList(disp_list_id);
        creating_display_list = true;
      }
    }

    // If body starts with a multi-geometry element, we have to draw each component
    // as separate 'body'!
    //
    BodyElementLoop* bel = getElementLoop(layer, 0);

    if ( bel == NULL )continue;

    BodyElement* be = bel->getElement(0);

    if ( be == NULL ) continue;

    int cmp_count = be->getNofComponents();

    for ( idx = 0; idx < cmp_count; idx++) {

      if ( drawing_cad_body ) {
        renderer->startDrawingCadBody();
      }

      // Loop all body element loops in the body layer
      //
      bool is_first_loop = true;

      for (IdList::iterator itr = beg; itr != end; itr++) {

        int bel_id = *itr;

        if ( bel_id < 0 )
          loopDirection = -1;
        else
          loopDirection = 1;

        BodyElementLoop* bel = model->getBodyElementLoopById(loopDirection * bel_id);

        // Pre-loop processing
        //
        if ( drawing_cad_body ) {
          renderer->startDrawingCadBodyElementLoop(is_first_loop);
        }

        // Draw all elements in the loop
        //
        int be_index;
        if ( loopDirection == 1 ) {
          be_index = 0;
        } else {
          be_index = bel->getNofElements() - 1;
        }

        while (true) {

          if ( loopDirection == 1 ) {
            bel->getElement(be_index++, d_be);
          } else {
            bel->getElement(be_index--, d_be);
          }

          BodyElement* be = d_be.element;
          bool mode_changed = false;

          if (be==NULL) break;

          // Let's make the final decision concerning bodyelement's ccw-direction
          // If it is in ccw-direction, we set: sign = 1.
          int direction = loopDirection * d_be.direction;

          objectDrawingMode be_dm = be->getDrawMode();

          if ( be_dm == DM_HIDDEN ) continue;

          // If wireframe mode is needed for this layer
          if ( !drawing_cad_body ||
               lr->isOpen() ||
               (geometry_type == DRAW_SOURCE_MESH && mesh_excluded)
             ) {

            objectDrawingState be_ds = be->getDrawState();

            if ( be_ds != DS_SELECTED ) {
              mode_changed = true;
              be->setDrawMode(DM_WIREFRAME);
            }
          }

          // Draw one body element
          //
          // If a 'mid' element has multi-geometry, we havet to draw them all
          // here (this is possible because their 'inner' components are not
          // visible at this level (they would be visible only ig tje element
          // were the first in the body)!
          //
          if ( cmp_count > 1 ) {
            be->draw(idx, renderer, geometry_type, id, direction, is_first_loop);

          } else {
            int  ccount = be->getNofComponents();
            for (i = 0; i < ccount; i++) {
              be->draw(i, renderer, geometry_type, id, direction, is_first_loop);
              if ( drawing_cad_body && i < ccount - 1 ) {
                renderer->stopDrawingCadBodyElementLoop();
                renderer->startDrawingCadBodyElementLoop(is_first_loop);
              }
            }
          }

          // Restore the original mode if changed for wireframe
          if ( mode_changed ) {
            be->setDrawMode(be_dm);
          }

        } //-All elements in the loop

        is_first_loop = false;

        // Post-loop processing
        if ( drawing_cad_body ) {
          renderer->stopDrawingCadBodyElementLoop();
        }

      } //-All loops in the layer

      if ( drawing_cad_body ) {
        renderer->stopDrawingCadBody();
      }

    } // By geometry component

    if ( creating_display_list ) {
     renderer->stopDisplayList(disp_list_id);
    }

    renderer->name_delete(lid);

  } //-All layer in the body

  renderer->name_delete(id);
}


void
Body::drawMesh(Renderer* renderer)
{
  if (nofMeshElements == 0) return;

  renderer->setLightUndirected();

  renderer->name_save(id);
  renderer->name_save(0); // Space for possible fem element id

  const Point3* node_data = model->getMeshNodeData();
  MeshElementTable* bt = model->getMeshBulkElements();
  //MeshElementTable* et = model->getMeshBulkElementEdges();

  int index;
  bool selected;

  for (int i = 0; i < nofMeshElements; i++) {

    index = meshElementIds[i];

    selected = false;

    if (drawState == DS_SELECTED && bt->selected[index]) {
      selected = true;
    }

    // Bulk elements are always drawn as wireframe!
    bt->drawEdges(renderer, node_data, index, selected);
  }

  renderer->name_delete(0);
  renderer->name_delete(id);

  renderer->setLightDirected();
}


#if 0
void
Body::drawMesh(Renderer* renderer)
{
  if (nofFemElements == 0)
    return;

  renderer->name_save(id);
  renderer->name_save(0);

  MeshBulkElementTable* elements = model->getMeshBulkElements();
  Point3* node_data = model->getMeshNodeData();

  for (int i = 0; i < nofFemElements; i++) {
    int index = femElementIds[i];
    int elem_code = elements->element_codes[index];
    int elem_type = MeshElementDesc[elem_code][DESC_ELEM_CODE];
    int* elem_nodes = elements->node_ids[index];

    // If body is selected, we store also
    // fem element ids
    if (drawMode == SELECTED) {
      renderer->name_replace(index);
    }

    objectDisplayMode display_mode;
    if ( elements->selected[index] )
      display_mode = HIGHLITED;
    else
      display_mode = NORMAL_DISPLAY;

    renderer->drawMeshBulkElement(elem_type, elem_nodes,
                                  display_mode, node_data);
  }

  renderer->name_delete(0);
  renderer->name_delete(id);
}
#endif


int
Body::getBodyForceId()
{
  return bodyForceId;
}


BodyLayer*
Body::getBodyLayer(int layer)
{
  if (!checkLayerIndex(layer)) return NULL;

  return model->getBodyLayerById(layerIds[layer]);
}


BoundBox*
Body::getBoundBox(int layer)
{
  if (!checkLayerIndex(layer)) return NULL;

  return boundboxes[layer];
}

// Int color coding
void
Body::getColor(Color4& clr)
{
  for (int i = 0; i < 4; i++)
    clr[i] = color[i];
}

// Float color coding
void
Body::getColor(Color4f& clr)
{
  for (int i = 0; i < 4; i++)
    clr[i] = float(color[i])/MAX_NOF_COLOR_LEVELS;
}


int
Body::getDirectedElementLoopId(int layer, int index)
{
  if (!checkLayerIndex(layer)) return 0;

  if ( index >= elementLoopIds[layer]->size() ) return 0;

  IdList::iterator itr = elementLoopIds[layer]->begin();

  for (int i = 0; i < index; i++, itr++);

  return *itr;
}


int
Body::getDirectedElementLoopTag(int layer, int index)
{
  if (!checkLayerIndex(layer)) return 0;

  int directed_id = getDirectedElementLoopId(layer, index);

  if ( directed_id == 0 ) return 0;

  int sign = (directed_id < 0)?-1:1;

  ModelObject* obj = model->getModelObjectById(sign * directed_id);

  if ( obj != NULL )
    return sign * obj->Tag();
  else
    return 0;
}


BodyElement**
Body::getElementArray()
{
  int tot_count = getElementCount();

  if (tot_count == 0) return NULL;

  BodyElement** be_array = new BodyElement*[tot_count];

  int index = 0;
  int pos = 0;
  while (true) {
    BodyElement* be = getElement(index++);
    if (be==NULL) break;
    be_array[pos++] = be;
  }

  return be_array;
}


BodyElement**
Body::getElementArray(int layer)
{
  if (!checkLayerIndex(layer)) return NULL;

  int row;

  int be_count = getElementCount(layer);

  if ( be_count == 0 ) return NULL;

  BodyElement** be_array = new BodyElement*[be_count];

  int be_index = 0;

  BodyElement* be = getElement(layer, be_index++);
  row = 0;
  while (be != NULL) {
    be_array[row++] = be;
    be = getElement(layer, be_index++);
  }

  return be_array;
}


int
Body::getElementCount()
{
  int count = 0;

  for (int i = 0; i < nofLayers; i++)
    count += nofElements[i];

  return count;
}


int
Body::getElementCount(int layer)
{
  if (!checkLayerIndex(layer)) return 0;

  return nofElements[layer];
}


// Get body element group id by index
//
int
Body::getElementGroupId(int index)
{
  if (index >= elementGroupIds->size()) return NO_INDEX;

  IdList::iterator itr = elementGroupIds->begin();

  for (int i = 0; i < index; i++, itr++);

  int eg_id = *itr;

  return eg_id;
}


// Get all body element groups ids
//
void
Body::getElementGroupIds(int& nof_groups, int*& group_ids)
{
  nof_groups = getNofElementGroups();
  group_ids = NULL;

  if ( nof_groups == 0 ) return;

  group_ids = new int[nof_groups];

  int counter = 0;

  for (int i = 0; i < nof_groups; i++) {
    group_ids[counter++] = getElementGroupId(i);
  }
}


// Get body element group tag by layer and index
//
int
Body::getElementGroupTag(int index)
{
  if (index >= getNofElementGroups() ) return NO_INDEX;

  IdList::iterator itr = elementGroupTags->begin();

  for (int i = 0; i < index; i++, itr++);

  int eg_tag = *itr;

  return eg_tag;
}


// Get all body element groups tags
//
void
Body::getElementGroupTags(int& nof_groups, int*& group_tags)
{
  nof_groups = getNofElementGroups();
  group_tags = NULL;

  if ( nof_groups == 0 ) return;

  group_tags = new int[nof_groups];

  int counter = 0;

  for (int i = 0; i < nof_groups; i++) {
    group_tags[counter++] = getElementGroupTag(i);
  }

}


// Get all (unsigned) body element ids in loop order
//
void
Body::getElementIds(int& nof_elements, int*& element_ids)
{
  nof_elements = getElementCount();
  element_ids = new int[nof_elements];

  int index = 0;

  int layer = -1;
  while (true ) {

    if (!checkLayerIndex(++layer)) break;

    // If loops created (normal bodies)
    //
    if ( getNofElementLoops(layer) > 0 ) {
      int bel_index = 0;
      BodyElementLoop* bel = getElementLoop(layer, bel_index++);

      while ( bel != NULL ) {

        for (int i = 0; i < bel->getNofElements(); i++) {
          element_ids[index++] = bel->getElementId(i);
        }

        bel = getElementLoop(layer, bel_index++);
      }

    // If no loops created (virtual bodies)
    //
    } else {
      int be_index = 0;
      while (true) {
        BodyElement* be = getElement(layer, be_index++);
        if (be == NULL ) break;
        element_ids[index++] = be->Id();
      }
    }
  }

}


// Get all (unsigned) layer element ids in loop order
//
void
Body::getElementIds(int layer, int& nof_elements, int*& element_ids)
{
  nof_elements = 0;
  element_ids = NULL;

  if (!checkLayerIndex(layer)) return;

  nof_elements = nofElements[layer];

  if ( nof_elements == 0 ) return;

  element_ids = new int[nof_elements];

  // If loops created (normal bodies)
  //
  if ( getNofElementLoops(layer) > 0 ) {

    int bel_index = 0;
    BodyElementLoop* bel = getElementLoop(layer, bel_index++);

    int index = 0;

    while ( bel != NULL ) {

      for (int i = 0; i < bel->getNofElements(); i++) {
        element_ids[index++] = bel->getElementId(i);
      }

      bel = getElementLoop(layer, bel_index++);
    }

  // If no loops created (virtual bodies)
  //
  } else {
    int index = 0;
    int be_index = 0;
    while (true) {
      BodyElement* be = getElement(layer, be_index++);
      if (be == NULL ) break;
      element_ids[index++] = be->Id();
    }
  }
}


// Get all (unsigned) body element tags in loop order
//
void
Body::getElementTags(int& nof_elements, int*& element_tags)
{
  nof_elements = getElementCount();
  element_tags = new int[nof_elements];

  int index = 0;

  int layer = -1;

  while (true ) {

    if (!checkLayerIndex(++layer)) break;

    // If loops created (normal bodies)
    //
    if ( getNofElementLoops(layer) > 0 ) {
      int bel_index = 0;
      BodyElementLoop* bel = getElementLoop(layer, bel_index++);

      while ( bel != NULL ) {

        for (int i = 0; i < bel->getNofElements(); i++) {
          element_tags[index++] = bel->getElementTag(i);
        }

        bel = getElementLoop(layer, bel_index++);
      }

    // If no loops created (virtual bodies)
    //
    } else {
      int be_index = 0;
      while (true) {
        BodyElement* be = getElement(layer, be_index++);
        if (be == NULL ) break;
        element_tags[index++] = be->Tag();
      }
    }
  }

}


// Get all (unsigned) layer element tags in loop order
//
void
Body::getElementTags(int layer, int& nof_elements, int*& element_tags)
{
  nof_elements = 0;
  element_tags = NULL;

  if (!checkLayerIndex(layer)) return;

  nof_elements = nofElements[layer];

  if ( nof_elements == 0 ) return;

  element_tags = new int[nof_elements];

  // If loops created (normal bodies)
  //
  if ( getNofElementLoops(layer) > 0 ) {

    int bel_index = 0;
    BodyElementLoop* bel = getElementLoop(layer, bel_index++);

    int index = 0;

    while ( bel != NULL ) {

      for (int i = 0; i < bel->getNofElements(); i++) {
        element_tags[index++] = bel->getElementTag(i);
      }

      bel = getElementLoop(layer, bel_index++);
    }

  // If no loops created (virtual bodies)
  //
  } else {
    int be_index = 0;
    int index = 0;
    while (true) {
      BodyElement* be = getElement(layer, be_index++);
      if (be == NULL ) break;
      element_tags[index++] = be->Tag();
    }
  }

}


int
Body::getElementLoopId(int index)
{
  int count = 0;
  int old_count;

  for (int i = 0; i < nofLayers; i++) {
    old_count = count;
    count += elementLoopIds[i]->size();

    if ( index < count )
      return getElementLoopId(i, index-old_count);
  }

  return NO_INDEX;
}


int
Body::getElementLoopId(int layer, int index)
{
  if (!checkLayerIndex(layer)) return NO_INDEX;

  if (index >= elementLoopIds[layer]->size()) return NO_INDEX;

  IdList::iterator itr = elementLoopIds[layer]->begin();

  for (int i = 0; i < index; i++, itr++);

  int bel_id = *itr;

  return ( bel_id < 0 )?-bel_id:bel_id;
}


int
Body::getElementLoopTag(int index)
{
  int count = 0;
  int old_count;

  for (int i = 0; i < nofLayers; i++) {
    old_count = count;
    count += elementLoopTags[i]->size();

    if ( index < count ) {
      return getElementLoopTag(i, index-old_count);
    }
  }

  return NO_INDEX;
}


int
Body::getElementLoopTag(int layer, int index)
{
  if (!checkLayerIndex(layer)) return NO_INDEX;

  if (index >= elementLoopTags[layer]->size()) return NO_INDEX;

  IdList::iterator itr = elementLoopTags[layer]->begin();

  for (int i = 0; i < index; i++, itr++);

  int bel_tag = *itr;

  return bel_tag;
}



int
Body::getEquationId()
{
  return equationId;
}


BodyElement*
Body::getElement(int index)
{
  int count = 0;
  int old_count;

  for (int i = 0; i < nofLayers; i++) {
    old_count = count;
    count += belements[i]->size();

    if ( index < count )
      return getElement(i, index-old_count);
  }

  return NULL;
}


BodyElement*
Body::getElement(int layer, int index)
{
  if (!checkLayerIndex(layer)) return NULL;

  if (index >= belements[layer]->size()) return NULL;

  BodyElementTable::iterator itr = belements[layer]->begin();

  for (int i = 0; i < index; i++, itr++);

  return (*itr).second;
}


BodyElement*
Body::getElementById(int be_id)
{
  BodyElement* be;

  for (int i = 0; i < nofLayers; i++) {

    be = getElementById(i, be_id);

    if ( be != NULL )
      return be;
  }

  return NULL;
}


BodyElement*
Body::getElementById(int layer, int be_id)
{
  if (!checkLayerIndex(layer)) return NULL;

  BodyElementTable* elements = belements[layer];

  BodyElement* be = NULL;

  if ( elements != NULL ) {
    be = (*elements)[be_id];
  }

  return be;
}


BodyElementLoop*
Body::getElementLoop(int index)
{
  int count = 0;
  int old_count;

  for (int i = 0; i < nofLayers; i++) {
    old_count = count;
    count += nofElementLoops[i];

    if ( index < count )
      return getElementLoop(i, index-old_count);
  }

  return NULL;
}


BodyElementLoop*
Body::getElementLoop(int layer, int index)
{
  if (!checkLayerIndex(layer)) return NULL;

  if (index >= nofElementLoops[layer]) return NULL;

  int bel_id = getElementLoopId(layer, index);

  return model->getBodyElementLoopById(bel_id);
}


BodyLayer*
Body::getLayer(int layer)
{
  if (!checkLayerIndex(layer)) return NULL;

  return model->getBodyLayerById(layerIds[layer]);
}


int
Body::getLayerId(int layer)
{
  if (!checkLayerIndex(layer)) return NO_INDEX;

  return layerIds[layer];
}


int
Body::getLayerIndexById(int layer_id)
{
  for (int i = 0; i < nofLayers; i++) {
    if (layerIds[i] == layer_id) {
      return i;
    }
  }

  return NO_INDEX;
}


int
Body::getLayerIndexByTag(int layer_tag)
{
  for (int i = 0; i < nofLayers; i++) {
    if ( layer_tag == getLayerTag(i) ) {
      return i;
    }
  }

  return NO_INDEX;
}


int
Body::getLayerTag(int layer)
{
  if (!checkLayerIndex(layer)) return NO_INDEX;

  return model->getModelObjectTagById(layerIds[layer]);
}


int
Body::getLayerTagById(int layer_id)
{
  return model->getModelObjectTagById(layer_id);
}


int
Body::getInitialConditionId()
{
  return initialConditionId;
}


int
Body::getMaterialId()
{
  return materialId;
}


int
Body::getMeshElementId(int index)
{
  if (index < 0 || index >= nofMeshElements) {
    return NO_INDEX;
  }

  return meshElementIds[index];
}

// Get mesh density (H) value for a layer
bool
Body::getMeshDensityValue(int layer, int mesh_index, char& dtype, double& dvalue)
{
  if (!checkLayerIndex(layer)) return false;

  if ( mesh_index < 0 || mesh_index >= model->getNofMeshes() )
    return false;

  dtype = ' ';
  dvalue = 0;

  BodyLayer* bl = model->getBodyLayerById(layerIds[layer]);

  if (bl == NULL ) return false;

  return bl->getMeshDensityValue(mesh_index, dtype, dvalue);
}


// Get boundary meshing N for a body layer
// If layer not 'griddable' return 0
int
Body::getMeshQuadGridN(int layer, int mesh_index, int element_id)
{
  if (!checkLayerIndex(layer)) return 0;

  BodyLayer* bl = model->getBodyLayerById(layerIds[layer]);

  if (bl == NULL ) return 0;

  return bl->getMeshQuadGridN(mesh_index, element_id);
}


// Body level name
const char*
Body::getName()
{
  if ( name == NULL || name[0] == '\0' ) {
    strstream strm;
    strm << "Body" << tag << ends;
    return strm.str();

  } else {
    return name;
  }
}


// Body layer name
const char*
Body::getName(int layer)
{
  if (!checkLayerIndex(layer)) return NULL;

  const char* nm = model->getModelObjectNameById(layerIds[layer]);

  if ( nm == NULL || nm[0] == '\0' ) {
    strstream strm;
    strm << "Body" << tag << "-L" << 1 + layer << ends;
    return strm.str();

  } else {
    return nm;
  }
}


int
Body::getNofBoundaries(int layer, elementType type)
{
  if (!checkLayerIndex(layer)) return 0;

  switch (type) {
  case INNER_ELEMENT:
    return nofInnerBoundaries[layer];
    break;
  case OUTER_ELEMENT:
    return nofOuterBoundaries[layer];
    break;
  default:
    return nofInnerBoundaries[layer] + nofOuterBoundaries[layer];
    break;
  }
}


int
Body::getNofElementGroups()
{
  return nofElementGroups;
}


int
Body::getNofElementLoops(int layer)
{
  if (!checkLayerIndex(layer)) return 0;

  return nofElementLoops[layer];
}


int
Body::getNofMifLayers()
{
  int count = 0;

  for (int i = 0; i < nofLayers; i++) {
    BodyLayer* bl = model->getBodyLayerById(layerIds[i]);
    if ( bl == NULL ) break;
    count += bl->getNofMifLayers(elementLoopIds[i]);
  }

  return count;
}


// Check if the body layer is inside a layer of an other body
bool
Body::hasInside(int layer, Body* other_body, int other_body_layer)
{
  BoundBox* box1 = getBoundBox(layer);
  BoundBox* box2 = other_body->getBoundBox(other_body_layer);

  if ( !box1->contains(box2) ) {
    return false;
  }

  int o_bel_index = 0;
  BodyElementLoop* o_bel = other_body->getElementLoop(other_body_layer, o_bel_index);

  BodyElement* o_be = o_bel->getElement(0);

  if ( o_be == NULL ) {
    return false;
  }

  BodyElement* o_se = o_be->getSubElement(0);

  if ( o_se == NULL ) {
    return false;
  }

  // Pick one vertex from other body element
  BodyElement* o_v;

  if ( model->getDimension() == ECIF_2D ) {
    o_v = o_se;
  } else {
    o_v = o_se->getSubElement(0);
  }

  if ( o_v == NULL ) {
    return false;
  }

  //---Calculate intersection from the vertex
  //  with all elements of this body layer
  GcPoint* vp = (GcPoint*)o_v->getGeometry();
  double vp_x = vp->Pos(X);

  //-Pick first (outermost is enough!) loop from this body layer
  int bel_index = 0;
  BodyElementLoop* bel = getElementLoop(layer, bel_index++);

  if ( bel == NULL ) {
    return false;
  }

  bool negative_on_left;

  Point3List hit_points;

  Point3* p1;
  Point3* p2;

  Point3List::iterator itr1;
  Point3List::iterator itr2;

  //-Add intersections with this element to the total count
  int index = 0;
  while (true) {

    BodyElement* be = bel->getElement(index++);

    if (be==NULL) break;

    RayHit* be_hits =  be->isectionsWithXdir(vp, negative_on_left);

    if ( be_hits != NULL && !isEqual(vp_x, be_hits->min_value) ) {

      for (itr1 = hit_points.begin(); itr1 != hit_points.end(); itr1++) {
          p1 = *itr1;
        for (itr2 = be_hits->points.begin(); itr2 != be_hits->points.end(); itr2++) {
          p2 = *itr2;

          if ( !isZero(dist3(*p1, *p2)) ) {
            hit_points.push_front(p2);
          }
        }
      }
    }
  }

  //--Non-zero and even number of intersections --> point is inside
  if ( 1 == hit_points.size() % 2 ) {
    return true;

  //--Otherwise it is outside
  } else {
    return false;
  }
}


void
Body::initClass(Model* mdl)
{
  Body::last_tag = 0;
}


void
Body::init(bool add_default_layer)
{
  int i;

  model->addModelObject(this, OT_BODY);

  name = NULL;

  bodyParameterId = NO_INDEX;
  boundbox = new BoundBox;
  boundboxMesh = new BoundBox;
  checked = false;

  colorIndex = DEFAULT_COLOR_INDEX;

  for (i = 0; i< 4; i++) {
    color[i] = colorValues[DEFAULT_COLOR_INDEX][i];
  }

  drawMode = DM_NORMAL;
  drawState = DS_NORMAL;
  externalTag = NO_INDEX;
  maxNofMeshElements = 0;

  nofMeshElements = 0;
  meshElementIds = NULL;

  nofLayers = 0;
  layerIds = NULL;

  bodyParameterId = NO_INDEX;
  initialConditionId = NO_INDEX;
  bodyForceId = NO_INDEX;
  equationId = NO_INDEX;
  materialId = NO_INDEX;

  boundboxes = new BoundBox*[nofLayers];
  belements = new BodyElementTable*[nofLayers];
  pendingElementTags = new IdArray*[nofLayers];
  pendingVertexGroups = new IdArray*[nofLayers];
  pendingVertexTags = new IdArray*[nofLayers];
  elementLoopIds = new IdList*[nofLayers];
  elementLoopTags = new IdList*[nofLayers];

  nofElements = new int[nofLayers];
  nofElementLoops = new int[nofLayers];
  nofInnerBoundaries = new int[nofLayers];
  nofOuterBoundaries = new int[nofLayers];

  for (i=0; i < nofLayers; i++) {

    layerIds[i] = NO_INDEX;

    boundboxes[i] = new BoundBox;
    belements[i] = new BodyElementTable;
    pendingElementTags[i] = new IdArray;
    pendingVertexGroups[i] = new IdArray;
    pendingVertexTags[i] = new IdArray;
    elementLoopIds[i] = new IdList;
    elementLoopTags[i] = new IdList;

    nofElements[i] = 0;
    nofElementLoops[i] = 0;
    nofInnerBoundaries[i] = 0;
    nofOuterBoundaries[i] = 0;
  }

  nofElementGroups = 0;
  elementGroupIds = new IdList;
  elementGroupTags = new IdList;

  separateFlag = false;
  status = false;

  type = NORMAL_BODY;
  if (model->modelHasCadGeometry()) {
    gmtrType = GEOM_BODY;
  } else {
    gmtrType = MESH_BODY;
  }

  tplgType = CLOSED_BODY;

  // Add default layer
  if (add_default_layer) {
    ecif_BodyLayer_X bl;
    init_trx_data(bl);
    addLayer(bl);
  }

}


void
Body::initName()
{

  if ( name == NULL || name[0] == '\0' ) {
    strstream strm;
    strm << "Body" << tag << ends;

    update_dyna_string(name, strm.str());
  }
}


void
Body::incrBndrCount(int layer, elementType type)
{
  if (!checkLayerIndex(layer)) return;

  switch (type) {
  case INNER_ELEMENT:
    nofInnerBoundaries[layer]++;
    break;
  case OUTER_ELEMENT:
    nofOuterBoundaries[layer]++;
    break;
  }
}


// If body is excluded from the currently drawn mesh
bool
Body::isExcludedFromCurrentMesh()
{
  int mesh_index = model->getCurrentMeshIndex();

  return isExcludedFromMesh(mesh_index);
}



// If body layer is excluded from the currently drawn mesh
bool
Body::isExcludedFromCurrentMesh(int layer)
{
  int mesh_index = model->getCurrentMeshIndex();

  return isExcludedFromMesh(layer, mesh_index);
}


// If body is excluded from the given mesh
bool
Body::isExcludedFromMesh(int mesh_index)
{
  if (mesh_index == NO_INDEX) return false;

  int layer = -1;
  while (true) {

    if (!selectLayer(++layer)) break;

    BodyLayer* bl = model->getBodyLayerById(layerIds[layer]);

    if (bl == NULL ) continue;

    // Body is 'included' if even one of the layers is active!!!
    if ( !bl->isExcludedFromMesh(mesh_index) ) {
      return false;
    }
  }

  return true;
}


// If body layer is excluded from the currently drawn mesh
bool
Body::isExcludedFromMesh(int layer, int mesh_index)
{
  if (!checkLayerIndex(layer)) return true;

  if (mesh_index == NO_INDEX) return false;

  BodyLayer* bl = model->getBodyLayerById(layerIds[layer]);

  if (bl == NULL ) return false;

  return bl->isExcludedFromMesh(mesh_index);
}


// Is this a 'decent' body?
bool
Body::isOk()
{
  if (!checked)
    status = check();
  return status;
}


bool
Body::isSeparable(int layer)
{
  if (!checkLayerIndex(layer)) return false;

  int be_index = 0;
  BodyElement* be = getElement(layer, be_index++);

  while ( be != NULL ) {

    if ( be->isInnerBoundary() ) {
      return true;
    }

    be = getElement(layer, be_index++);
  }

  return false;
}


void
Body::markActiveMeshObjects(int mesh_index, bool*& active_object_flags)
{

  // Mark body active
  active_object_flags[id] = true;

  int layer = -1;
  while (true) {

    if (!selectLayer(++layer)) break;

    BodyLayer* bl = model->getBodyLayerById(layerIds[layer]);

    if (bl == NULL ) continue;;

    if ( !bl->isExcludedFromMesh(mesh_index) ) {

      // Mark all body related elements active
      int be_index = 0;
      BodyElement* be = getElement(layer, be_index++);

      while ( be != NULL ) {

        if ( !be->isActive() ) {
          be = getElement(layer, be_index++);
          continue;
        }

        active_object_flags[be->Id()] = true;

        be->markActiveMeshObjects(active_object_flags);

        be = getElement(layer, be_index++);
      }
    }
  }
}


void
Body::markActiveObjects()
{
  int layer = 0;
  if (!selectLayer(layer)) return;

  while (true) {

    int bel_index = 0;
    BodyElementLoop* bel = getElementLoop(layer, bel_index++);
    while (bel != NULL) {
      int bel_id = bel->Id();
      model->markObjectActive(bel_id);
      bel = getElementLoop(layer, bel_index++);
    }

    if (!selectLayer(++layer)) break;
  }
}


ostream&
Body::output_emf(ostream& out, short indent_size, short indent_level)
{
  char* QM = "\"";
  short is = indent_size;
  short il = indent_level;

  int i;

  const short nof_parameter_ids = 5;

  const char* param_names[nof_parameter_ids] = {
    "Body Parameter", "Body Force", "Equation", "Initial Condition", "Material"
  };

  int param_ids[nof_parameter_ids] = {
    bodyParameterId, getBodyForceId(), getEquationId(), getInitialConditionId(), getMaterialId()
  };

  // Header
  LibFront::output_scalar(out, is, il, EMF_BODY, NULL, tag);

  // Type
  if ( isOpen() ) {
    LibFront::output_scalar(out, is, il + 1, EMF_TYPE, NULL, "Open");
  } else if ( isBemBody() ) {
    LibFront::output_scalar(out, is, il + 1, EMF_TYPE, NULL, "Bem");
  } else if ( isVirtual() ) {
    LibFront::output_scalar(out, is, il + 1, EMF_TYPE, NULL, "Virtual");
  }

  // Body name
  LibFront::output_scalar(out, is, 1 + il, EMF_NAME, NULL, name);

  // We do not print size-info for the color vector
  LibFront::output_vector(out, is, 1 + il, EMF_COLOR, NULL, 4, color, false);

  // Parameter (equation,body force etc.)ids
  for (i = 0; i < nof_parameter_ids; i++) {
    if (param_ids[i] != NO_INDEX) {
      LibFront::output_scalar(out, is, 1 + il, param_names[i], NULL, param_ids[i]);
    }
  }

  if ( isVirtual() ) {
    int eg_count = getNofElementGroups();

    if ( eg_count > 0 ) {

      int* eg_tags = new int[eg_count];
      for (int i = 0; i < eg_count; i++) {
        eg_tags[i] = getElementGroupTag(i);
      }

      LibFront::output_vector(out, is, il + 1, "Element Groups", NULL, eg_count, eg_tags, false);

      delete[] eg_tags;
    }

    return out;
  }


  int layer, nof_ids;
  int* mids;
  int* gids;

  // Loop all layers
  // ===============
  layer = -1;
  int lr_tag_prev = NO_INDEX;

  while (true) {

    if (!selectLayer(++layer)) break;

    int lr_id = layerIds[layer];

    BodyLayer* bl = model->getBodyLayerById(lr_id);

    if ( bl == NULL ) continue;

    int lr_tag = bl->Tag();

    // Store directed element loop and element group tags into an arries
    //
    int bel_count = elementLoopTags[layer]->size();
    int* bel_tags = new int[bel_count];
    for (i = 0; i < bel_count; i++) {
      bel_tags[i] = getDirectedElementLoopTag(layer, i);
    }

    //--Header and tag
    if ( EXPLICIT_LAYER == bl->getLayerType() && lr_tag != lr_tag_prev ) {
      lr_tag_prev = lr_tag;
      LibFront::output_scalar(out, is, 1 + il, EMF_LAYER, NULL, lr_tag);

      // If name given explicitely
      if ( bl->hasName() ) {
        LibFront::output_scalar(out, is, 1 + il, EMF_LAYER_NAME, NULL, bl->getName());
      }
    }

    //--Directed element loop tags and gridparameter ids etc. from layer.
    //
    int ill = il;

    if ( EXPLICIT_LAYER == bl->getLayerType() )
      ill += 2;
    else
      ill += 1;

    if ( bel_count > 0 ) {
      LibFront::output_vector(out, is, ill, "Element Loops", NULL, bel_count, bel_tags, false);
    }

    bl->output_emf(out, is, ill);

    delete[] bel_tags;

  } // All Layers

  return out;
}


ostream&
Body::output_mif(ostream& out)
{
  char* QM = "\"";

  int i;
  char value_buffer[1 + 128];

  // Read bodylevel ghrid info from the first layer
  //
  BodyLayer* bl = getBodyLayer(0);
  if ( bl == NULL ) return out;
  int mesh_index = model->getCurrentMeshIndex();
  int pid = bl->getGridParameterId(mesh_index);
  Parameter* grid_param = model->getParameterById(ECIF_GRID_PARAMETER, pid);

  // Body id
  // -------
  indent(out, 0) << "BodyId: " << tag << endl;

  // Element order (common for all layers)
  // -------------
  indent(out, 2) << "ElementOrder: ";

  if ( grid_param != NULL &&
       grid_param->getFieldValueBySifName("Mesh Element Order", 128, value_buffer)
     )
    out << value_buffer;
  // Default
  else
    out << "Linear";
  out << endl;

  // Nof layers
  // ----------
  indent(out, 2) << "Layers: " << getNofMifLayers() << endl;

  int next_layer_id = 1;

  // Output data for each Layer
  // ==========================
  int layer = -1;
  while (true) {

    if (!selectLayer(++layer)) break;

    BodyLayer* bl = getBodyLayer(layer);

    if ( bl == NULL ) break;

    bl->output_mif(out, next_layer_id, elementLoopIds[layer]);

  } // All body Layers

  return out;
}


#if 0
ostream&
Body::output_mif(ostream& out)
{
  char* QM = "\"";

  int i;
  char value_buffer[1 + 128];

  int mesh_index = model->getCurrentMeshIndex();

  int layer = -1;

  // Output data for each Layer
  // ==========================
  while (true) {

    if (!selectLayer(++layer)) break;

    Parameter* grid_param = NULL;

    BodyLayer* bl = getBodyLayer(layer);

    if ( bl == NULL ) continue;

    int pid = bl->getGridParameterId(mesh_index);

    grid_param = model->getParameterById(ECIF_GRID_PARAMETER, pid);

    // ----------------------------------------
    // Output body level data based on 1. layer
    // ----------------------------------------

    if ( layer == 0 ) {

      // Body id
      // -------
      indent(out, 0) << "BodyId: " << tag << endl;

      // Element order (common for all layers)
      // -------------
      indent(out, 2) << "ElementOrder: ";

      if ( grid_param != NULL &&
           grid_param->getFieldValueBySifName("Mesh Element Order", 128, value_buffer)
         )
        out << value_buffer;
      // Default
      else
        out << "Linear";
      out << endl;

      // Nof layers
      // ----------
      indent(out, 2) << "Layers: " << nofLayers << endl;

    }

    char dtype;
    double dvalue;
    bool dgiven = getMeshDensityValue(layer, mesh_index, dtype, dvalue);

    bool is_quadGrid;
    // Check if Quadrilateral elements
    if ( grid_param != NULL &&
         grid_param->getFieldValueBySifName("Mesh Element Type", 128, value_buffer) &&
         LibFront::ncEqual(value_buffer, "Quad")
       ) {
      is_quadGrid = true;
    } else {
      is_quadGrid = false;
    }

    // Layer id
    // --------
    indent(out, 4) << "LayerId: " << layer + 1 << endl;

    // Layer mesh density value
    // ------------------------
    if ( dgiven && dvalue > 0 )
      indent(out, 4) << dtype << ": " << dvalue << endl;

    // Layer type (meshing method)
    // ---------------------------
    indent(out, 6) << "LayerType: ";

    //-Structured grid
    if ( is_quadGrid ) {
      out << "QuadGrid";

      indent(out, 6) << "GridSize: ";
      grid_param->getFieldValueBySifName("Mesh Quadgrid N1", 128, value_buffer);
      out << value_buffer << " ";
      grid_param->getFieldValueBySifName("Mesh Quadgrid N2", 128, value_buffer);
      out << value_buffer << endl;

    //-Triangles
    } else if ( grid_param != NULL &&
                grid_param->getFieldValueBySifName("Mesh Layer Type", 128, value_buffer)
              ) {
      out << value_buffer;

    //-Default
    } else {
      out << "MovingFront";
    }

    out << endl;

    // Background mesh stuff
    // ---------------------
    if ( !is_quadGrid ) {

      //-Fixed nodes
      indent(out, 6) << "FixedNodes: " << 0 << endl;  // NOTE: Currently always 0!, MVe 16.02.00

      //-Background mesh type
      indent(out, 6) << "BGMesh: ";

      //--Background mesh method specified
      if ( grid_param != NULL &&
           grid_param->getFieldValueBySifName("Mesh Bg Mesh", 128, value_buffer)
         ) {

        //-External bg-mesh file
        if ( LibFront::in(value_buffer, "External") ) {

          //-File name given
          if ( grid_param->getFieldValueBySifName("Mesh Bg Mesh File", 128, value_buffer) &&
               value_buffer != NULL && value_buffer[0] != '\0'
             ) {
            out << "External " << QM << value_buffer << QM;

          //-File name missing, use default method!
          } else {
            out << "Delaunay";
          }

        //-Some of the predefined methods selected by the user
        } else {
          out << value_buffer;
        }

      //--Use default method (no bg mesh)
      } else {
        out << "Delaunay";
      }

      out << endl;

      // Mesh seed (optional)
      // --------------------
      if ( grid_param != NULL &&
           grid_param->getFieldValueBySifName("Mesh Seed Type", 128, value_buffer)
         ) {

        // Only Implicit seed currently supported, MVe 16.02.99
        if ( LibFront::ncEqual(value_buffer, "Implicit") &&
             grid_param->getFieldValueBySifName("Mesh Seed Edge", 128, value_buffer)
           ) {
          indent(out,6) << "Seed: Implicit" << endl;
          indent(out,8) << "Edge: " << value_buffer << endl;
        }
      }
    }

    // Layer element loops
    // -------------------
    indent(out, 6) << "Loops: " << nofElementLoops[layer] << endl;
    for (i = 0; i < nofElementLoops[layer]; i++) {

      int did = getDirectedElementLoopId(layer, i);
      int direction = (did < 0 )?-1:1;

      BodyElementLoop* bel = model->getBodyElementLoopById(direction * did);

      if ( bel == NULL )
        continue;

      // Print positions before and after printing fixed data
      // (for getting the indent size for ids)
      int pos1, pos2;

      indent(out, 8) << "LoopId: " << bel->Tag() << endl;
      indent(out, 10) << "Direction: " << direction << endl;

      pos1 = out.tellp();
      indent(out, 10) << "Edges: " << bel->getNofElements() << "  ";
      pos2 = out.tellp();

      bel->outputDirectedElementTags(out, pos2 - pos1);
      out << endl;
    }

  } // All body Layers

  return out;
}
#endif


ostream&
Body::output_sif(ostream& out, short indent_size, short indent_level)
{
  // NOTE: If no equation defined, the body is not output to sif!!!
  //
  if ( NO_INDEX == getEquationId() ) {
    return out;
  }

  char* QM = "\"";
  short is = indent_size;
  short il = indent_level;

  int nof_parameter_ids = 4;
  const char* names[4] = { SIF_BODY_FORCE, SIF_EQUATION, SIF_INITIAL_CONDITION, SIF_MATERIAL};
  int ids[4] = { getBodyForceId(), getEquationId(), getInitialConditionId(), getMaterialId()};

  // Body header and id
  LibFront::output_scalar(out, is, il, SIF_BODY, NULL, tag);

  // Name (typed and in quotes)
  if ( model->getSolverKeywordTypeGiven(SIF_BODY, "Name") ) {
    LibFront::output_scalar(out, is, il + 1, "Name =", NULL, name, true);
  } else {
    LibFront::output_scalar(out, is, il + 1, "Name = String", NULL, name, true);
  }
  out << endl;

  // Body parameters
  if ( bodyParameterId != NO_INDEX ) {
    Parameter* p = model->getParameterById(ECIF_BODY_PARAMETER, bodyParameterId);

    if ( p != NULL ) {

      SifOutputControl soc;
      soc.outputId = false;
      soc.outputName = false;
      soc.outputType = false;
      soc.outputAll = false;
      soc.sectionName = SIF_BODY;

      p->output_sif(out, indent_size, indent_level, soc);
    }
  }

  // Parameter ids
  for (short i = 0; i < nof_parameter_ids; i++) {

    strstream strm;

    // NOTE: Equation is an exception: id = -1 is output!
    // This mainly for convenience, so that the euquation can be
    // nicely turned on/off when using sif-file directly
    //
    //
    if ( 0 == strcmp(SIF_EQUATION, names[i]) || ids[i] != NO_INDEX ) {

      strm << names[i] << " =";
      if ( !model->getSolverKeywordTypeGiven(SIF_BODY, names[i]) ) {
        strm << " Integer" << ends;
      }
      strm << ends;
      LibFront::output_scalar(out, is, 1 + il, strm.str(), NULL, ids[i]);
    }
  }

  return out;
}


bool
Body::selectLayer(int layer)
{
  return checkLayerIndex(layer);
}


// NOTE: Take care that bodies are
// checked (method check()) before they are
// used, when any elements have been removed!!!
int
Body::removeElement(int layer, BodyElement* be)
{
  if (!checkLayerIndex(layer)) return 0;

  belements[layer]->erase(be->Id());
  nofElements[layer]--;

  return 0;
}


// NOTE: We cannot remove the group object themselves, because
// they could be in use somewhere else!
//
void
Body::removeElementGroups()
{
  nofElementGroups = 0;

  delete elementGroupIds;
  delete elementGroupTags;
}


void
Body::removeElementLoops(int layer)
{
  if (!checkLayerIndex(layer)) return;

  if (nofElementLoops[layer] > 0 ) {

    IdList::iterator beg = elementLoopIds[layer]->begin();
    IdList::iterator end = elementLoopIds[layer]->end();

    for (IdList::iterator itr = beg; itr != end; itr++) {

      int bel_id = *itr;

      if (bel_id < 0)
        bel_id = -1 * bel_id;

      //model->removeBodyElementLoop(bel_id);
    }

    elementLoopIds[layer]->clear();
    elementLoopTags[layer]->clear();
    nofElementLoops[layer] = 0;
  }
}


// Remove all layer data from body
//
// NOTE: Layers themselves are removed as model object by the model!
//
//
void
Body::removeLayers()
{
  for (int i = 0; i < nofLayers; i++) {

    delete belements[i];
    delete boundboxes[i];
    delete elementLoopIds[i];
    delete elementLoopTags[i];
    delete pendingElementTags[i];
    delete pendingVertexGroups[i];
    delete pendingVertexTags[i];
  }

  nofLayers = 0;

  delete[] belements;
  delete[] elementLoopIds;
  delete[] elementLoopTags;
  delete[] layerIds;
  delete[] nofElements;
  delete[] nofElementLoops;
  delete[] nofInnerBoundaries;
  delete[] nofOuterBoundaries;
}


int
Body::removeMeshElements(bool* remove_flags)
{
  int counter, i;

  //--Count nof remaining ids
  counter = 0;
  for (i = 0; i < nofMeshElements; i++) {

    if (remove_flags[meshElementIds[i]])
      continue;

    counter++;
  }

  //--Allocate new table
  int* new_ids = new int[counter];

  //--Copy new ids
  counter = 0;
  for (i = 0; i < nofMeshElements; i++) {

    if (remove_flags[meshElementIds[i]])
      continue;

    new_ids[counter] = meshElementIds[i];
    counter++;
  }

  //--Update vars
  nofMeshElements = counter;
  delete[] meshElementIds;
  meshElementIds = new_ids;

  return nofMeshElements;
}


// fix this!!!***!!!
void
Body::separate(bool* duplicated_object_flags)
{
  int layer = 0;

  if (!selectLayer(layer)) return;

  IdArray ids1, ids2;

  BodyElementLoop* bel = model->getBodyElementLoopById( getElementLoopId(layer, 0) );

  int index = 0;
  while (true) {

    BodyElement* be = bel->getElement(index++);

    if (be==NULL) break;

    if ( !be->isInnerBoundary() ) {
      continue;
    }

    int neighbor_id;

    if ( id == be->getParentId(1) ) {
      neighbor_id = be->getParentId(2);

    } else {
      neighbor_id = be->getParentId(1);
    }

    Body* neighbor = model->getBodyById(neighbor_id);

    if ( neighbor->getSeparateFlag() ) {
      continue;
    }

    int neighbor_tag = neighbor->Tag();

    BodyElement* new_be = be->duplicate(duplicated_object_flags);

    if ( new_be != NULL ) {

      ids1.push_back(be->Id());
      ids2.push_back(new_be->Id());

      new_be->setParentIds(id, NO_INDEX);
      new_be->setParentTags(tag, NO_INDEX);

      be->setParentIds(neighbor_id, NO_INDEX);
      be->setParentTags(neighbor_tag, NO_INDEX);

      removeElement(layer, be);
      addElement(layer, new_be);
    }
  }

  if ( 0 < ids1.size() ) {
    bel->swapElements(&ids1, &ids2);
  }
}


void
Body::setBodyParameterId(int pid)
{
  //model->updateParameterApplyCounts(A_BODY, id, ECIF_BODY_PARAMETER, bodyParameterId, pid);
  bodyParameterId = pid;
}


void
Body::setColor(Color4& bd_color)
{
  for (int i = 0; i< 4; i++)
    color[i] = bd_color[i];

  colorIndex = ef_nodefault;
}


void
Body::setColorIndex(colorIndices color_index)
{
  // Set color index
  colorIndex = color_index;
  // Set actual color
  for (int i = 0; i< 4; i++) {
    color[i] = colorValues[color_index][i];
  }
}


void
Body::setElementLoops(int layer, int nof_loops, int* loop_ids, int* loop_tags)
{
  if (!checkLayerIndex(layer)) return;

  nofElementLoops[layer] = nof_loops;

  delete elementLoopIds[layer];
  delete elementLoopTags[layer];

  elementLoopIds[layer] = new IdList;
  elementLoopTags[layer] = new IdList;

  for(int i=0; i < nofElementLoops[layer]; i++) {
    elementLoopIds[layer]->push_back(loop_ids[i]);
    elementLoopTags[layer]->push_back(loop_tags[i]);
  }

}



void
Body::setMeshElements(int nof_elements, int* element_ids)
{
  delete[] meshElementIds; meshElementIds = NULL;

  nofMeshElements = nof_elements;

  meshElementIds = element_ids;
}


void
Body::setLayerTag(int layer, int tag)
{
  if (!checkLayerIndex(layer)) return;

  model->setModelObjectTagById(layerIds[layer], tag);
}


void
Body::swapBodyElements(int layer, BodyElement* be1, BodyElement* be2)
{
  if (!checkLayerIndex(layer)) return;

  if (be1 == be2) {
    return;
  }

  be1->addStatus(BE_SWAPPED);

  removeElement(layer, be1);
  addElement(layer, be2);
}


// Method updates body's:
//  -bounding box

// Update must be done only after when all elements
// have been read into the model
void
Body::updateBoundBox()
{
  int i;
  int direction;
  RangeVector rv;

  int layer = -1;

  // Each Layer
  while (true) {

    if (!selectLayer(++layer)) break;

    IdList::iterator end = elementLoopIds[layer]->end();

    for (IdList::iterator itr = elementLoopIds[layer]->begin(); itr != end; itr++) {

      int bel_id = *itr;

      direction = (bel_id < 0)?-1:1;

      BodyElementLoop* bel = model->getBodyElementLoopById(direction * bel_id);

      if (bel == NULL) continue;

      int index = 0;
      while (true) {
        BodyElement* be = bel->getElement(index++);
        if (be==NULL) break;
        if ( be->getRangeVector(rv) ) {
          updateBoundBox(layer, rv);
        }
      }
    }

  } // Each layer
}


void
Body::updateBoundBox(int layer, RangeVector rv)
{
  boundbox->extendByRange(rv);

  if (!checkLayerIndex(layer)) return;

  boundboxes[layer]->extendByRange(rv);
}


void
Body::updateBodyElementMeshH(int layer, int updateMode)
{
  if (!checkLayerIndex(layer)) return;

  int nof_meshes = model->getNofMeshes();

  int be_index = 0;
  BodyElement* be = getElement(layer, be_index++);

  while(be != NULL) {

    for (int i = 0; i < nof_meshes; i++) {

      //double mesh_h = getMeshH(i);

      //be->setMeshH(i, mesh_h, updateMode);
    }

    be = getElement(layer, be_index++);
  }
}
