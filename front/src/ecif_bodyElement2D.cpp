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
Module:     ecif_bodyelement2D.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_body.h"
#include "ecif_bodyElement1D.h"
#include "ecif_bodyElement2D.h"
#include "ecif_control.h"
#include "ecif_geometry.h"
#include "ecif_parameter.h"
#include "ecif_renderer.h"

extern Control* theControlCenter;


//Initialize static class variables.
int BodyElement2D::last_tag = 0;


// Constructors.

BodyElement2D::BodyElement2D()
  : BodyElement()
{
  tag = newTag();
  init();

  update();
}


//-Boundary with given id
BodyElement2D::BodyElement2D(int int_tag)
  : BodyElement(int_tag)
{
  checkLastTag();
  init();

  update();
}

 
//-Topology and geometry directly defined
BodyElement2D::BodyElement2D(int v1_id, int v2_id, Geometry* pG, int cd, char* nm)
  : BodyElement(cd, nm)
{
  tag = newTag();
  init();

  if ( v1_id != v2_id ) {
    modelData->nofSubElements = 2;
    modelData->subElementIds->push_back(v1_id);
    modelData->subElementIds->push_back(v2_id);

    modelData->nofVertices = 2;
    modelData->vertexIds->push_back(v1_id);
    modelData->vertexIds->push_back(v2_id);

  } else {
    modelData->nofSubElements = 1;
    modelData->subElementIds->push_back(v1_id);

    modelData->nofVertices = 1;
    modelData->vertexIds->push_back(v1_id);
  }

  ptrGmtr = pG;
  
  calcBoundaryPoints();

  update();
}


//-Create a line or multiline element
BodyElement2D::BodyElement2D(int nof_vertices, int* vertex_ids)
: BodyElement()
{
  tag = newTag();

  if (nof_vertices <= 0 )
    return;

  init();

  modelData->nofSubElements = nof_vertices;
  modelData->nofVertices = nof_vertices;

  BodyElement** vertices = new BodyElement*[nof_vertices];

  for (int i = 0; i < nof_vertices; i++) {

    modelData->subElementIds->push_back(vertex_ids[i]);
    modelData->vertexIds->push_back(vertex_ids[i]);

    vertices[i] = model->getVertexById(vertex_ids[i]);
  }

  if ( nof_vertices == 2 ) {
    ptrGmtr = new GcLine(vertices[0], vertices[1]);
  } else {
    ptrGmtr = new GcPolyLine(nof_vertices, vertices);
  }

  delete[] vertices;
  
  calcBoundaryPoints();

  update();
}


//-Create a line-element from vertices
BodyElement2D::BodyElement2D(int v1_id, int v2_id)
  : BodyElement()
{
  tag = newTag();
  init();

  modelData->nofSubElements = 2;
  modelData->subElementIds->push_back(v1_id);
  modelData->subElementIds->push_back(v2_id);

  modelData->nofVertices = 2;
  modelData->vertexIds->push_back(v1_id);
  modelData->vertexIds->push_back(v2_id);

  BodyElement* v1 = model->getVertexById(v1_id);
  BodyElement* v2 = model->getVertexById(v2_id);

  // Geometry, line segment
  ptrGmtr = new GcLine(v1, v2);
  
  calcBoundaryPoints();

  update();
}


//-Create an edge element
BodyElement2D::BodyElement2D(int v1_id, int v2_id, ecif_EdgeGeometry_X* params)
  : BodyElement()
{
  tag = newTag();
  init();

  modelData->nofSubElements = 2;
  modelData->subElementIds->push_back(v1_id);
  modelData->subElementIds->push_back(v2_id);

  modelData->nofVertices = 2;
  modelData->vertexIds->push_back(v1_id);
  modelData->vertexIds->push_back(v2_id);

  BodyElement* v1 = model->getVertexById(v1_id);
  BodyElement* v2 = model->getVertexById(v2_id);

  // Geometry by type
  switch (params->type) {

  case ECIF_CIRCLE:
    ptrGmtr = new GcCircle(v1, v2, params);
    break;
  case ECIF_NURBS:
    ptrGmtr = new GcNurbsCurve(v1, v2, params);
    break;
  default:
    ptrGmtr = new GcLine(v1, v2);
    break;
  }

  calcBoundaryPoints();

  update();
}


//-Create a bodyelement from a model-file structure.
BodyElement2D::BodyElement2D(ecif_Element_X& tx)
  : BodyElement(tx)
{
  UserInterface* gui = (UserInterface*)model->getGui();
  strstream strm;

  IdList vertex_tags;

  if ( tag == NO_INDEX ) {
    tag = newTag();
  }

  checkLastTag();

  init(tx.name);

  // Store extra vertices as subelements
  //
  for (int i = 0; i < tx.nof_extra_vertices; i++) {
    modelData->subElementIds->push_back(tx.extra_vertex_tags[i]);
    modelData->nofSubElements++;
  }

  // Create element geometry
  // =======================

  // NOTE: For 3D-model edges, no actual geometry
  ptrGmtr = NULL;

  if ( tx.nof_components > 0 ) {

    //---Single component geometry which is not function
    //
    if ( tx.nof_components == 1 && !tx.components[0]->isFunction ) {

      switch ( tx.components[0]->gmtr_type ) {

      case ECIF_CIRCLE:
        ptrGmtr = new GcCircle(*tx.components[0], vertex_tags);
        break;

      case ECIF_LINE:
        ptrGmtr = new GcLine(*tx.components[0], vertex_tags);
        break;

      case ECIF_POLYLINE:
        ptrGmtr = new GcPolyLine(*tx.components[0], vertex_tags);
        break;

      case ECIF_NURBS:
        ptrGmtr = new GcNurbsCurve(*tx.components[0], vertex_tags);

        // NOTE: Nurbs are replace with linearized geometry!!!
        formLinearizedGeometry();

        break;

      default:
        break;
      }

    //---Multi component geometry or a function
    //
    } else {
        ptrGmtr = new GcMulti2D(tx, vertex_tags);
    }

  } // If geometry components

  if ( ptrGmtr != NULL && !ptrGmtr->geometryIsOk() ) {
    objectOk = false;
    strm << "***Error in Edge " << tag << ends;
    gui->showMsg(strm.str());
    return;
  }

  // Store Vertices
  // ==============
  modelData->nofSubElements += vertex_tags.size();
  modelData->nofVertices += vertex_tags.size();;

  IdList::iterator itr = vertex_tags.begin();

  // NOTE: Tags will be converted later to ids!!!
  while ( itr != vertex_tags.end() ) {

    int vtag = *itr++;

    modelData->subElementIds->push_back(vtag);
    modelData->vertexIds->push_back(vtag);
  }

  // Linearize geometry for drawing and for Mesh2D if needed
  if ( ptrGmtr != NULL ) {
    ptrGmtr->calcBoundaryPoints();
  }

  update();
}


//-Empty mesh boundary with given id
BodyElement2D::BodyElement2D(int int_tag, int parent1_tag, int parent2_tag, int nof_mesh_elements)
  : BodyElement(int_tag)
{
  checkLastTag();
  init();

  modelData->parent1Tag = parent1_tag;
  modelData->parent2Tag = parent2_tag;

  allocateMeshElements(nof_mesh_elements);

  update();
}


//-Empty mesh boundary
BodyElement2D::BodyElement2D(int parent1_tag, int parent2_tag, int nof_mesh_elements)
  : BodyElement()
{
  tag = newTag();
  init();

  modelData->parent1Tag = parent1_tag;
  modelData->parent2Tag = parent2_tag;

  allocateMeshElements(nof_mesh_elements);

  update();
}


BodyElement2D::~BodyElement2D()
{
  delete ptrGmtr;
  delete gridHData;
  delete modelData;
  delete meshData;
}


// Method updates the list which contains info on all bodyelements
// into which this bodyelement is divided (ie. inner & outer boundaries).
// Normally we first add inner boundaries into this list and outer boundaries
// are added as a 'residual'. List is kept in 'consecutive' order.
// This method is called (at least) from the MODEL when
// inner boundaries are searched and created.
// Element is marked 'unchecked' if a new common boundary is added
// because we have to check outer-boundary etc. after that.
//
void
BodyElement2D::addCoveringElement(BodyElement* covering_elem, beStatus se_stat)
{
  BodyElement2D* se = (BodyElement2D*) covering_elem;

  IdList* se_list = model->getCoveringElementList(this->id);
  
  bool is_new_list = false;

  // If this is first sub-element.
  if (se_list == NULL) {
    se_list = new IdList;
    is_new_list = true;
    //model->addCoveringElementList(this->id, se_list);
    modelData->status |= BE_DEVIDED;
  }

  // Sub-element status is updated.
  se->addStatus(se_stat);

  // sub-element's direction compared to owner element (*this*)
  int se_direction = calcDirection(se);

  // New element is added in 'ascending' oder.
  IdList::iterator pos = se_list->begin();
  bool inserted = false;

  while (pos != se_list->end()) {

    // Remeber the current position and get item at it.
    IdList::iterator cur_pos = pos++;

    int oe_id = *cur_pos;

    if ( oe_id == 0 ) {
      continue;
    }

    int oe_sign = (oe_id < 0)?-1:1;
    oe_id = oe_sign * oe_id;
    BodyElement* oe = model->getEdgeById(oe_id);

    // If new item is 'smaller' than the retrieved item, we put it
    // before the retrieved and stop the loop.
    // NOTE: We add directed id!
    //
    if (isBefore(se, se_direction, oe, oe_sign, ptrGmtr)) {
      // STL's insert is "before"?
      se_list->insert(cur_pos, se_direction * se->Id());
      inserted = true;
      break;
    }
  }

  // If we didn't insert above, we have to add to tail.
  if (!inserted) {
    se_list->push_back(se_direction * se->Id());
  }

  // Element must be checked after an add-operation!
  modelData->status = modelData->status & ~BE_OUTER_CHECKED;

  if ( is_new_list ) {
    model->addCoveringElementList(this->id, se_list);
  }
  
}


// Method calculates the direction of the covering element
// compared to *this* element.
// Conclusion is done by calculating the parametric-values for covering element's
// start- and end-points in the geometry of the parent element.
int
BodyElement2D::calcDirection(BodyElement* covering_elm)
{
  int direction = 0;

  // Note: parent's geometry id given as the parmeter.
  ParamValues* pv = covering_elm->getParamValues(ptrGmtr);

  // If we don't have an edge, or edges are not
  // adjacent, we can't say anything definite
  if (pv->count != 2 ||
     pv->values[0] == NULL ||
     pv->values[1] == NULL)
     return 0;

  // Ok, direction are comparable
  double u1 = pv->values[0]->u;
  double u2 = pv->values[1]->u;

  direction = (u1 <= u2)?1:-1;

  return direction;
}


// Check if boundaries should be linearized with new mesh-h/mesh-n
//
void
BodyElement2D::checkBoundaryDiscretization(int mesh_index)
{
  // NOT IN USE!!!
  // *************
  return;

  char type;
  int nof_values;
  double values[4];

  // If Edge has not a gridH parameter attached, set default discretization!
  //
  if ( !( getMeshDensityValues(mesh_index, type, nof_values, values) &&
          nof_values == 1 &&
          values[0] > 0
        )
     ) {
    ptrGmtr->setDeltaU(LIN_DELTA_NONE, 0.0);

  // Otherwise, apply the mesh density value
  //
  } else {
    double meshH = model->getMeshH(mesh_index);
    double meshF = model->getMeshF(mesh_index);
    double value = values[0];
  
    if ( type == 'N' ) {
      ptrGmtr->setDeltaU(LIN_DELTA_N, value);

    } else if ( type == 'H' ) {
      ptrGmtr->setDeltaU(LIN_DELTA_H, value * meshF);

    } else if ( type == 'R' ) {
      ptrGmtr->setDeltaU(LIN_DELTA_H, value * meshH * meshF);
    }
  }

  // If geometry component should now use fixed discretization
  // for ELmerMesh2D  (N: n style in mif)
  //
  if ( ptrGmtr->useFixedMeshN() ) {
    ptrGmtr->calcBoundaryPoints();
  }
}


// Method checks if element is Ok and if there is some outer boundary.
// Outer boundary elements are created if necessary between two
// inner boundary sub-elements.!!!
// Outer-boundary info is updated into status-attribute
// Returns true if outer-boundies exist.
bool
BodyElement2D::checkOuterBoundaries()
{
  if ( modelData->status & (BE_OUTER_CHECKED | !BE_INCLUDES_OUTER) )
    return true;

  // Mark that I'm outer-boundary checked!
  modelData->status |= BE_OUTER_CHECKED;

  // If this is an inner boundary, it is a single
  // element and cannot include (or be itself such)
  // outer-boundaries
  if ( modelData->status & BE_INNER )
    return true;

  // No sub-elements, no inner boundary:
  // --> element itself is a single outer boundary;
  if ( !(modelData->status & (BE_DEVIDED | BE_INNER)) ) {
    modelData->status |= (BE_INCLUDES_OUTER | BE_OUTER);
    return true;
  }

  // Now we loop through sub-elements-list and check if something
  // is missing between two consecutive inner-boundary-elements
  // NOTE: this works in 2D, but in 3D!!!
  IdList* se_list = model->getCoveringElementList(this->id);

  if ( se_list == NULL ) {
    return false;
  }

  int se_id;
  int se_sign;
  BodyElement2D *be, *se;

  // We start from body-element's starting vertex.
  // Casting casting!, well, but we know that we have edges, don't we ....
  BodyElement* v1 = getFirstSubElement();
  BodyElement* v2 = NULL;

  bool isAtParentEnd = false;

  IdList::iterator pos = se_list->begin();

  while (true) {

    // If there still is a sub-element, we compare *v1* to its
    // starting vertex ...
    IdList::iterator cur_pos;
    if (pos != se_list->end()) {

      // Current position is remembered
      cur_pos = pos;

      // Current sub-element is copied and pos -iterator is increased
      se_id = *pos++;

      if ( se_id == 0 ) {
        continue;
      }

      se_sign = (se_id < 0)?-1:1;
      se_id = se_sign * se_id;
      se = (BodyElement2D*)model->getEdgeById(se_id);

      if (se_sign == 1)
        v2 = se->getFirstSubElement();
      else
        v2 = se->getLastSubElement();
    }

    // ... otherwise we compare bodyelement's ending vertex.
    // NOTE: because element is BE_DEVIDED, *v1* cannot any more
    // be bodyelements's starting vertex, so there is no danger
    // to create this element again!
    else {
      isAtParentEnd = true;
      v2 = getLastSubElement();
    }

    // If vertices are not the same, a new sub-element is created.
    // ... and this an outer-boundary!!!
    if (v1 != v2) {

      modelData->status |= BE_INCLUDES_OUTER;

      //--New body-element
      be = (BodyElement2D*) createElement(v1->Id(), v2->Id(), ptrGmtr->getType());

      model->addBodyElement(be);

      //--It will be a SUB-element of an existing element,
      //  and we don't add it explicitely to the parent body.
      be->addStatus(BE_OUTER);

      // If we were able to read the 'next' sub-element for comparison
      // new element must be inserted before it, otherwise it is added
      // at the tail.
      if (!isAtParentEnd) {
        // STL's insert is "before"?
        se_list->insert(cur_pos, be->Id());

      } else {
        se_list->push_back(be->Id());
      }

    }

    // We break if the last sub-element is handled and we at the end
    // of the parent element!!!
    // NOTE: pos -iterator is the correct position to compare
    // because cur_pos was possibly decreased and we possibly pushed
    // something at the end of the list
    if (isAtParentEnd)
      break;

    // Otherwise we continue and change the comparison vertex to the
    // end of the current sub-element.
    if (se_sign == 1)
      v1 = se->getLastSubElement();
    else
      v1 = se->getFirstSubElement();
  }

  return true;
}


// This doesn't do anything useful currently !!!###!!!
// It should check if two body-elements are similarly
// oriented in the space (normals are in same side etc.)
int
BodyElement2D::compareOrientation(BodyElement* other)
{
  return 0;
}


#if 0
void BodyElement2D::draw(Renderer* renderer, flagName geometry_type, int body_id, int direction)
{
  Body* bd = model->getBody(body_id);

  // If this is a mesh boundary, no name, just draw
  if ( geometry_type == DRAW_SOURCE_MESH &&
       !bd->isExcludedFromCurrentMesh()
     ) {

    drawMesh(renderer, body_id, direction);
  }

  //--Element itself is and name id saved.
  else {

    // Store element name and draw
    renderer->name_save(id);
    ptrGmtr->draw(renderer, drawMode, direction);

    // Draw sub elements (vertices)
    drawSubElements(renderer, geometry_type, body_id, direction);

    renderer->name_delete(id);
  }

}
#endif


// Method creates a new 2D edge-form bodyelement from vertices.
BodyElement*
BodyElement2D::createElement(int v1_id, int v2_id, ecif_geometryType gt)
{
  BodyElement* v1 = model->getVertexById(v1_id);
  BodyElement* v2 = model->getVertexById(v2_id);

  Geometry* g;

  switch (gt) {

  case ECIF_LINE:
      g = new GcLine(v1,v2);
      break;

    default:
      return NULL;
      break;
  }
 
  BodyElement* be = new BodyElement2D(v1_id, v2_id, g);

  return be;
}


int
BodyElement2D::findMeshBorderNodes(int buf_size, int* ids_buffer)
{
  // In 2D this is easy!
  for (int i = 0; i < meshData->nofMeshBorderElements; i++) {

    int node_id = meshData->meshBorderElementIds[i];

    ids_buffer[node_id] = node_id;
  }

  return meshData->nofMeshBorderElements;
}


// Method picks outer-boundary (sub)elements within a bodyelement.
// Returns a list of bodyelements.
BodyElementList*
BodyElement2D::findOuterBoundary()
{
  if ( !(modelData->status & (BE_INCLUDES_OUTER)) ||
       modelData->status & BE_INNER
    )
    return 0;

  BodyElementList* oe_list = new BodyElementList;

  // If bodyelement itself is an outer boundary ...
  if ( modelData->status & BE_OUTER ) {
    oe_list->push_front(this);
    return oe_list;
  }

  // ... otherwise we select outer-boundaries from coveringElements-list.
  IdList* se_list = model->getCoveringElementList(this->tag);
  IdList::iterator pos = se_list->begin();
  while (pos != se_list->end()) {
    int se_id = *pos++;
    if (se_id < 0)
      se_id = -1 * se_id;

    BodyElement2D* se = (BodyElement2D*)model->getEdgeById(se_id);

    if ( se->getStatus() & BE_OUTER )
      oe_list->push_back(se);
  }
  return oe_list;
}


int
BodyElement2D::getMifGeometryTag(int index) const
{
  if ( index < 0 || index > nofMifTags - 1 ) return NO_INDEX;

  return mifTags[index];
}


int
BodyElement2D::getNofMifGeometries() const
{
  if ( ptrGmtr == NULL ) {
    return 0;
  } else {
    return ptrGmtr->getNofComponents();
  }
}

double
BodyElement2D::getParamArea(Geometry* gp)
{
  // Calculate the 'arc length' of the bodyelement (edge)
  // in the argument's parametric space;

  GcPoint* p1 = (GcPoint*)getFirstSubElement()->getGeometry();
  GcPoint* p2 = (GcPoint*)getLastSubElement()->getGeometry();

  ParamPair* pp0 = gp->point2Param(p1);
  ParamPair* pp1 = gp->point2Param(p2);

  // For a edge we need only the u-value.
  double u0 = pp0->u;
  double u1 = pp1->u;

  // Result should be non-negative.
  if (u0 < u1)
    return (u1 - u0);
  else
    return (u0 - u1);
}


ParamValues*
BodyElement2D::getParamValues(Geometry* gp)
{
  ParamValues* pv = new ParamValues;

  pv->count = 2;
  pv->values = new ParamPair*[pv->count];

  pv->t_type = ECIF_EDGE;
  pv->g_type = gp->getType();

  GcPoint* p1 = (GcPoint*)getFirstSubElement()->getGeometry();
  GcPoint* p2 = (GcPoint*)getLastSubElement()->getGeometry();

  pv->values[0] = gp->point2Param(p1);
  pv->values[1] = gp->point2Param(p2);

  return pv;
}


BodyElement*
BodyElement2D::getSubElement(int index)
{
  if ( modelData == NULL || modelData->nofSubElements == 0 )
    return NULL;

  if ( index < 0 || index >= modelData->nofSubElements )
    return NULL;

  return model->getVertexById((*modelData->subElementIds)[index]);
}


void
BodyElement2D::init(char* be_name)
{
  model->addModelObject(this, OT_EDGE);

  initName(be_name);

  modelData = new BodyElementModelData;
  meshData = new BodyElementMeshData;

  nofMifTags = 0;
  mifTags = NULL;
}


void
BodyElement2D::initClass(Model* model)
{
  BodyElement2D::last_tag = 0;
}


// Method checks if (sub)element *se1* is 'before' (sub)element *se2*.
// Conclusion is based on parametric values calculated in the
// geometry *gp*, which is given as a parameter. Normally this refers
// to parent element's geometry.
bool
BodyElement2D::isBefore(BodyElement* se1, int dir1,
                BodyElement* se2, int dir2,
                Geometry* gp)
{
  bool result = false;

  bool canCalculate = true;

  ParamValues* pv1 = se1->getParamValues(gp);
  ParamValues* pv2 = se2->getParamValues(gp);

  if (pv1->count != 2 || pv2->count != 2 ||
     pv1->t_type != ECIF_EDGE || pv2->t_type != ECIF_EDGE
     || (dir1 == 0)
     || (dir2 == 0)
    )
    canCalculate = false;

  if (canCalculate) {
    int indx1 = (dir1 == 1)?0:1;
    int indx2 = (dir2 == 1)?0:1;
    if (pv1->values[indx1]->u < pv2->values[indx2]->u)
      result = true;
  }

  delete pv1;
  delete pv2;

  return result;
}

#if 0
bool
BodyElement2D::isBemBoundary()
{
  if ( model->getDimension() == ECIF_3D )
    return true;
  else
    return BodyElement::isBemBoundary();
}
#endif

// Method checks if bodyelement is Ok.
bool
BodyElement2D::isOk()
{
//  if ( ! (status & BE_OUTER_CHECKED) )
//    checkOuterBoundary();
  return ( 0 != modelData->status & BE_OUTER_CHECKED );
}


bool
BodyElement2D::isOnSameAxis(GcPoint& p1, GcPoint& p2)
{
  return ptrGmtr->isOnSameAxis(p1, p2);
}


// Method checks if two linear edges are adjacent.
// Returns the adjacent section as a new bodyelement or 0.
matchType
BodyElement2D::matchToLinear(BodyElement* be2, BodyElement*& common )
{

  // No subelements, no hope!
  if ( getNofSubElements() == 0 || be2->getNofSubElements() == 0 ) {
    return MATCH_NONE;
  }

  // A polyline is currently nerver matching to any other line!
  if ( getNofVertices() > 2 || be2->getNofVertices() > 2 ) {
    return MATCH_NONE;
  }

  // Some nice casts :-), well both are lines!
  GcLine* g1 = (GcLine*) ptrGmtr;
  GcLine* g2 = (GcLine*) be2->getGeometry();

  BodyElement* v11 = getFirstSubElement();
  BodyElement* v12 = getLastSubElement();
  BodyElement* v21 = be2->getFirstSubElement();
  BodyElement* v22 = be2->getLastSubElement();

 
  GcPoint *p11, *p12, *p21, *p22; //, *dir1, *dir2;
  GcPoint dir1, dir2; // NOTE: for SGI, pointers won't work for it!

  p11 = (GcPoint*)v11->getGeometry();
  p12 = (GcPoint*)v12->getGeometry();
  p21 = (GcPoint*)v21->getGeometry();
  p22 = (GcPoint*)v22->getGeometry();

  //dir1 = &(*p12 - *p11);
  dir1 = *p12 - *p11;
  //dir2 = &(*p22 - *p21);
  dir2 = *p22 - *p21;

  common = NULL;

  //------Negative cases are checked first
  //---Lines must be parallel at least.
  if (! (&dir1)->isParallel(&dir2))
    return MATCH_NONE;

  //---Their start-points must lay on the same line.
  //if (! dir1->isParallel(&(*p11 - *p21)))
  GcPoint tmp = *p11 - *p21;
  if (! (&dir1)->isParallel(&tmp))
    return MATCH_NONE;


  //------Positive cases:
  //---Case1: Element be2 is inside be1.
  //   Result is be2 or smaller of the ids if match is exact
  if (this->hasInside(be2)) {
    // Exact match
    if (be2->hasInside(this)) {
      if (tag < be2->Tag() )
        common = this;
      else
        common = be2;
      return  MATCH_EXACT;
    }
    // be2 is inside this
    else {
      common = be2;
      return  MATCH_2_INSIDE;
    }
  }

  //---Case2: Element be1 is inside be2. Result is be1.
  if (be2->hasInside(this)) {
    common = this;
    return  MATCH_1_INSIDE;
  }

  //---Case3: Lines are overlapping. A new element as a result!
  BodyElement *v1, *v2;

  // One of the be1's vertices must be within the line of be2
  if (p11->isBetween(p21,p22))
    v1 = v11;
  else
    v1 = v12;
  // and one of the be'2s vertices must be within the line of be1.
  if (p21->isBetween(p11,p12))
    v2 = v21;
  else
    v2 = v22;

  // As a result we return a completely new 2D bodyelement.
  Geometry* g  = new GcLine(v1, v2);
  BodyElement* be = new BodyElement2D(v1->Id(), v2->Id(), g);
  common = be;

  return MATCH_OVERLAP;
}


// Method checks if two nurbs-curve edges are adjacent.
// Returns match type, currently only MATCH_NONE or MATCH_EXACT
// ie. we cannot split nurb-curves!!!
matchType
BodyElement2D::matchToNurbs(BodyElement* be2, BodyElement*& common)
{
  common = NULL;
  bool sameDirection = true;

  // Currently we check if curves are completely adjacent (exact match)
  // ie. they start and end at the same points.
  // This test is naturally too strict and
  // we should accept partial adjacency !!!***!!!
  GcNurbsCurve* g1 = (GcNurbsCurve*) this->ptrGmtr;
  GcNurbsCurve* g2 = (GcNurbsCurve*) be2->getGeometry();

  BodyElement* v11 = getFirstSubElement();
  BodyElement* v12 = getLastSubElement();
  BodyElement* v21 = be2->getFirstSubElement();
  BodyElement* v22 = be2->getLastSubElement();

  GcPoint *p11, *p12, *p21, *p22;

  p11 = (GcPoint*)v11->getGeometry();
  p12 = (GcPoint*)v12->getGeometry();
  p21 = (GcPoint*)v21->getGeometry();
  p22 = (GcPoint*)v22->getGeometry();

  //---End points must at least be the same
  //-At first we swap other element's points if
  // index-1 points are different
  if ((const GcPoint&)(*p11) != (const GcPoint&)(*p21)) {
    GcPoint* tmp = p21;
    p21 = p22;
    p22 = tmp;
    sameDirection = false;
  }

  //-Then we continue to check in the normal way;
  if ((const GcPoint&)(*p11) != (const GcPoint&)(*p21) ||
     (const GcPoint&)(*p12) != (const GcPoint&)(*p22)
    )
    return MATCH_NONE;

  //---Next we calculate a sample of points on the curves and
  //   test that they are the same
  int smpl_size = 10;
  GcPoint* testp1;
  ParamPair* param2;
  double start2_u = 0.0;

  double delta;

  if ( smpl_size > 0 )
    delta = 1.0 / smpl_size;
  else
    delta = 1.0;

  int direction = 1;

  if (!sameDirection) {
    start2_u = 1.0;
    direction = -1;
  }

  for (int i = 1; i < smpl_size; i++) {

    testp1 = g1->param2Point(i * delta, 0.0);
    param2 = g2->point2Param(testp1, direction, start2_u);
    start2_u = param2->u;

    delete param2;
    delete testp1;

    if (start2_u == NSVD)
      return MATCH_NONE;
  }

  //-----Well, we really have same nurbs curve!!!
  common = be2;

  return MATCH_EXACT;
}


// Output edge for Elmer Mesh input file (mif-file)
//
ostream&
BodyElement2D::output_mif(ostream& out)
{
  if ( ptrGmtr == NULL )return out;

  strstream strm;

  //--Boundary tag is output first into the suffix string
  //
  strm << boundaryTag;

  // If geometry component uses fixed Mesh-N, it will handle
  // all the rest of the output!
  //
  if ( ptrGmtr->useFixedMeshN() ) {
    strm << ends;
    return ptrGmtr->output_mif(out, strm.str(), true);
  }

  //--Check possible mesh-density value
  //
  char type;
  int nof_values;
  double values[4];

  int mesh_index = model->getCurrentMeshIndex();

  // If mesh density value set, output value here (like H: 0.2, N:50 etc.)
  //
  if ( getMeshDensityValues(mesh_index, type, nof_values, values) &&
       nof_values == 1 &&
       values[0] > 0
     ) {

    double value = values[0];

    // Force user defined mesh-h to obey boundary discretization
    //
    // NOTE: NOT IN USE!!!
    // *******************
    //
    if ( type == 'H' &&
         isDiscretized() &&
         ptrGmtr->getMeshH() > 0
       ) {
      //value = ptrGmtr->getMeshH();
    }

    strm << "  " << type << ": " << value;

  // Force a mesh-h which obeys boundary discretization
  //
  // NOTE: NOT IN USE!!!
  // *******************
  //
  } else if ( isDiscretized() && ptrGmtr->getMeshH() > 0 ) {
    double value = ptrGmtr->getMeshH();
    //strm << "  " << 'H' << ": " << value;
  }

  // NOTE: Multi2D geometry may output multiple instances into mif-file
  // using same boundary tag!
  //
  strm << ends;
  return ptrGmtr->output_mif(out, strm.str(), false);
}


// Set edge mif-tag into elements geometry
//
void
BodyElement2D::setMifTag(int& next_tag)
{
  nofMifTags = 0;
  delete[] mifTags;
  mifTags = NULL;

  if ( ptrGmtr != NULL ) {
    ptrGmtr->setMifTag(next_tag);
    ptrGmtr->getMifTags(nofMifTags, mifTags);
  }

}


