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
Module:     ecif_bodyelement3D.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_body.h"
#include "ecif_bodyElement2D.h"
#include "ecif_bodyElement3D.h"
#include "ecif_control.h"
#include "ecif_geometry.h"
#include "ecif_mesh.h"
#include "ecif_renderer.h"

extern Control* theControlCenter;


//Initialize static class variables.
int BodyElement3D::last_tag = 0;

// Constructors.

BodyElement3D::BodyElement3D()
: BodyElement()
{
  tag = newTag();
  init();
  update();
}

//-Geometry directly defined
BodyElement3D::BodyElement3D(int nof_vertices, int* vertex_ids, Geometry* pG, int cd, char* nm)
  : BodyElement(cd, nm)
{
  tag = newTag();
  init();

  ptrGmtr = pG;

  calcBoundaryPoints();

  update();
}


//-Create a plane-element
BodyElement3D::BodyElement3D(int nof_vertices, int* vertex_ids)
: BodyElement()
{
  tag = newTag();

  if (nof_vertices != 4 )
    return;

  init();

  GcPoint* points[4];

  modelData->nofSubElements = 4;
  modelData->nofVertices = 4;

  for (int i = 1; i < 4; i++) {
    modelData->subElementIds->push_back(vertex_ids[i]);
    modelData->vertexIds->push_back(vertex_ids[i]);
    BodyElement* be = model->getVertexById(vertex_ids[i]);
    points[i] = (GcPoint*)be->getGeometry();
  }

  // Geometry, plane
  Geometry* pg = new GcPlane(4, points);

  ptrGmtr = pg;

  calcBoundaryPoints();

  update();
}


//-Create a nurbs-srf element
BodyElement3D::BodyElement3D(int nof_vertices, int* vertex_ids, ecif_FaceGeometry_X* params)
  : BodyElement()
{
  tag = newTag();

  if (nof_vertices != 4 )
    return;

  init();

  GcPoint* points[4];

  modelData->nofSubElements = 4;
  modelData->nofVertices = 4;

  for (int i = 1; i < 4; i++) {
    modelData->subElementIds->push_back(vertex_ids[i]);
    modelData->vertexIds->push_back(vertex_ids[i]);
    BodyElement* be = model->getVertexById(vertex_ids[i]);
    points[i] = (GcPoint*)be->getGeometry();
  }

  // Geometry, nurbs-curve
  Geometry* pg = new GcNurbsSurface(params);

  ptrGmtr = pg;

  calcBoundaryPoints();

  update();
}


//-Create a bodyelement from model-file structure.
BodyElement3D::BodyElement3D(ecif_Element_X& tx)
  : BodyElement(tx)
{
  int i;

  IdList all_vertex_tags;

  if ( tag == NO_INDEX ) {
    tag = newTag();
  }

  checkLastTag();

  init(tx.name);

  ptrGmtr = NULL;

  // If we have some eal geometry
  //
  // NOTE: Currently we cannot have any 3D Cad-geometry!
  //
  if ( tx.nof_components > 0 ) {

    switch (tx.components[0]->gmtr_type) {

    case ECIF_LINE:
      ptrGmtr = new GcPlane(*tx.components[0], all_vertex_tags);
      break;
    case ECIF_NURBS:
      ptrGmtr = new GcNurbsSurface(*tx.components[0], all_vertex_tags);
      break;
    default:
      ptrGmtr = NULL; // Mesh "geometry"
      break;
    }

  }


  // NOTE: This is actually incorrect!!!
  // Subelements should be edges, but currently we cannot have 3D cad-geometry, so this
  // does not matter!

  // Vertices
  // ========
  modelData->nofSubElements = all_vertex_tags.size();
  modelData->nofVertices = all_vertex_tags.size();;

  IdList::iterator itr = all_vertex_tags.begin();

  // NOTE: Tags will be converted later to ids!
  while ( itr != all_vertex_tags.end() ) {
    modelData->subElementIds->push_back(*itr);
    modelData->vertexIds->push_back(*itr);
    itr++;
  }

  //formLinearizedGeometry();
  calcBoundaryPoints();

  update();
}


//-Empty mesh boundary with given id
BodyElement3D::BodyElement3D(int int_tag, int parent1_tag, int parent2_tag,int nof_mesh_elements)
  : BodyElement(int_tag)
{
  checkLastTag(tag);

  init();

  modelData->parent1Tag = parent1_tag;
  modelData->parent2Tag = parent2_tag;

  allocateMeshElements(nof_mesh_elements);

  update();
}


//-Mesh boundary (not fem elements added yet)
BodyElement3D::BodyElement3D(int parent1_tag, int parent2_tag, int nof_mesh_elements)
  : BodyElement()
{
  tag = newTag();
  init();

  modelData->parent1Tag = parent1_tag;
  modelData->parent2Tag = parent2_tag;

  allocateMeshElements(nof_mesh_elements);

  update();
}


BodyElement3D::~BodyElement3D()
{
  delete ptrGmtr;
  delete gridHData;
  delete modelData;
  delete meshData;
}


// Add all pending element ids
// Returns nof elements added
int
BodyElement3D::addAllPendingSubElements()
{
  int i;
  int nof_elements = 0;

  for (i = 0; i < pendingEdgeTags.size(); i++) {
    int be_tag = pendingEdgeTags[i];
    BodyElement* be = model->getBodyElementByTag(OT_EDGE, be_tag);
    if (be != NULL) {
      addSubElement(be->Id());
      nof_elements++;
    }
  }
  pendingEdgeTags.clear();

  for (i = 0; i < pendingVertexTags.size(); i++) {
    int be_tag = pendingVertexTags[i];
    BodyElement* be = model->getBodyElementByTag(OT_VERTEX, be_tag);
    if (be != NULL) {
      addSubElement(be->Id());
      nof_elements++;
    }
  }
  pendingVertexTags.clear();

  return nof_elements;
}


// Method updates the list which contains info on all bodyelements
// into which this bodyelement is divided (ie. inner & outer boundaries).
// Normally we first add inner boundaries into this list and outer boundaries
// are added as a 'residual'. List is kept in 'consecutive' order.
// This method is called (at least) from the MODEL when
// inner boundaries are searched and created.
// Element is marked 'unchecked' if a new common boundary is added
// because we have to check outer-boundary etc. after that.
void
BodyElement3D::addCoveringElement(BodyElement* covering_elem, beStatus se_stat)
{
  BodyElement3D* se = (BodyElement3D*) covering_elem;
  IdList* se_list = model->getCoveringElementList(this->id);

  // If this is first sub-element.
  if (se_list == 0) {
    se_list = new IdList;
    model->addCoveringElementList(this->tag, se_list);
    modelData->status |= BE_DEVIDED;
  }

  // Sub-element status is updated.
  se->addStatus(se_stat);

  // sub-element's direction compared to owner element (*this*)
  int se_direction = calcDirection(se);
  // New element is added in 'ascending' oder.
  IdList::iterator pos = se_list->begin();
  bool inserted = false;

  // If we didn't insert above, we have to add to tail.
  if (!inserted)
    se_list->push_back(se_direction * se->Id());

  // Element must be checked after an add-operation!
  modelData->status = modelData->status & ~BE_OUTER_CHECKED;
}

void
BodyElement3D::addMeshBorderAsSubElement()
{
  if ( meshData == NULL ||
       meshData->nofMeshBorderElements == 0
     ) {
    return;
  }

  // Create new body element, set this element's parent info also
  // inot new element!
  //
  int nof_elems = meshData->nofMeshBorderElements;

  BodyElement* se = new BodyElement2D(getParentTag(1), NO_INDEX, nof_elems);
  
  se->setParentId(1, getParentId(1));
  se->setParentLayer(1, getParentLayer(1));

  // Add into this element and into the model
  //
  addSubElement(se);
  model->addBodyElement(se);

  // Add border elements as mesh elements into new sub element
  //
  for (int i = 0; i < nof_elems; i++) {
    se->addMeshElement(meshData->meshBorderElementIds[i], 
                       meshData->meshBorderElementDirs[i]);
  }
}

 
// NOTE: Take care that pending edges tags are
// finally added as edge ids to the element
//
int
BodyElement3D::addPendingEdge(int edge_tag)
{
  pendingEdgeTags.push_back(edge_tag);

  return 0;
}


// NOTE: Take care that pending vertex tags are
// finally added as vertex ids to the element
//
int
BodyElement3D::addPendingVertex(int vertex_tag)
{
  pendingVertexTags.push_back(vertex_tag);

  return 0;
}



// Method calculates the direction of the sub-element
// compared to *this* element.
// Conclusion is done by calculating the parametric-values for subelement's
// start- and end-points in the geometry of the parent element.
int
BodyElement3D::calcDirection(BodyElement* sub_elm)
{
  int direction = 1;

  // Note: parent's geometry id given as the parmeter.
  ParamValues* pv = sub_elm->getParamValues(ptrGmtr);

  // If we don't have an edge, or edges are not
  // adjacent, we can't say anything definite
  if (pv->count != 2 ||
     pv->values[0] == NULL ||
     pv->values[1] == NULL)
     return 0;

  //Ok, direction are comparable
  //double u1 = pv->values[0]->u;
  //double u2 = pv->values[1]->u;
  //direction = (u1 <= u2)?1:-1;

  return direction;
}


// Method checks if element is Ok and if there is some outer boundary.
// Outer boundary elements are created if necessary between two
// inner boundary sub-elements.!!!
// Outer-boundary info is updated into status-attribute
// Returns true if outer-boundies exist.
bool
BodyElement3D::checkOuterBoundaries()
{
  if ( modelData->status & (BE_OUTER_CHECKED | !BE_INCLUDES_OUTER) )
    return false;

  // Mark that I'm outer-boundary checked!
  modelData->status |= BE_OUTER_CHECKED;

  // If this is an inner boundary, it is a single
  // element and cannot include (or be itself such)
  // outer-boundaries
  if ( modelData->status & BE_INNER )
    return false;

  // No sub-elements, no inner boundary:
  // --> element itself is a single outer boundary;
  if ( !(modelData->status & (BE_DEVIDED | BE_INNER)) ) {
    modelData->status |= (BE_INCLUDES_OUTER | BE_OUTER);
    return true;
  }

  // Return value based on checked, current outer-boundary status
  return ( 0 != modelData->status & BE_INCLUDES_OUTER );
}


// This doesn't do anything useful currently !!!###!!!
// It should check if two body-elements are similarly
// oriented in the space (normals are in same side etc.)
int
BodyElement3D::compareOrientation(BodyElement* other)
{
  return 0;
}


// Method creates a new 3D face from vertices.
BodyElement*
BodyElement3D::createElement(int nof_vertices, int* vertex_ids, ecif_geometryType gt)
{
  if (nof_vertices < 3)
    return NULL;

  GcPoint** points = new GcPoint*[nof_vertices];

  for (int i = 0; i < nof_vertices; i++) {
    BodyElement* v = model->getVertexById(vertex_ids[i]);
    points[i] = (GcPoint*)v->getGeometry();
  }

  BodyElement* be = NULL;
  Geometry* gmtr = NULL;

  switch (nof_vertices) {

  // Triangular
  case 3:
    break;

  // Quadrilateral
  case 4:

    // Here we should create first edges from vertices
    //tplg = new TcFace(nof_vertices, vertices);

    switch (gt) {

    case ECIF_LINE:
      gmtr = new GcPlane(4, points);
      break;

    default:
      gmtr = new GcPlane(4, points);
      break;
    }

  default:
    break;
  }

  if ( gmtr != NULL )
    be = new BodyElement3D(nof_vertices, vertex_ids, gmtr);

  return be;
}


#if 0
void BodyElement3D::draw(Renderer* renderer, flagName geometry_type, int body_id, int direction)
{
  // If this is a mesh boundary
  if (geometry_type == DRAW_SOURCE_MESH) {
    drawMesh(renderer, body_id, direction);
    return;
  }

  // Topology is here always supposed to be face-topology!!!
  int* vertices = 0;
  // There are no sub-elements. Only element itself is drawn.
  if ( !( modelData->status & BE_DEVIDED ) ) {
    ptrGmtr->draw(renderer, drawMode, direction, this->id);
    return;
  }

  // There are covering elements. Each covering element must be drawn separately
  // so that its geometry and direction is taken into account.
  else {
    IdList* se_list = model->getCoveringElementList(this->id);
    Geometry* pg;
    beStatus sub_direction;
    BodyElement* se;
    int se_id;
    int se_sign;

    // Because covering elements are stored according to parent-element's
    // direction, we have to take into account parent-element's direction
    // in the body's edege-loop when reading covering elements.
    IdList::iterator pos = se_list->begin();
    IdList::reverse_iterator rpos = se_list->rbegin();

    while (pos != se_list->end() &&
         rpos != se_list->rend()) {

      // Reading forwards
      if (direction == 1)
        se_id = *pos++;

      // Reading backwards
      else
        se_id = *rpos++;

      // Check sign i.e. direction compared to parent element.
      se_sign = (se_id < 0)?-1:1;

      // Find sub-element
      se = model->getFace(se_sign * se_id);

      // Sub-element is handled
      pg = se->getGeometry();

      // store element name
      renderer->name_replace(se->Id());

      // NOTE positive direction is originally ccw-direction and this
      // is how OpenGL wants to see edge-loop data.
      sub_direction = direction * se_sign;
      pg->draw(renderer, drawMode, sub_direction, this->id);
    }
  }
}
#endif


int
BodyElement3D::findMeshBorderNodes(int buf_size, int* ids_buffer)
{
  // In 3D border elements are edges!
  MeshElementTable* elements = model->getMeshBoundaryElementEdges();

  int node_count = 0;

  // Loop all border elements
  for (int i = 0; i < meshData->nofMeshBorderElements; i++) {

    int elem_id = meshData->meshBorderElementIds[i];
    meshElementCode elem_code = elements->getElementCode(elem_id);
    int nof_nodes = MeshElementDesc[elem_code][DESC_NOF_NODES];

    const int* node_ids = elements->getNodeIds(elem_id);

    // Loop all nodes in each element
    for (int j = 0; j < nof_nodes; j++) {
      int node_id = node_ids[j];

      // If already stored
      if ( ids_buffer[node_id] != NO_INDEX ) {
        continue;
      }

      node_count++;
      ids_buffer[node_id] = node_id;
    }
  }

  return node_count;
}


// Method picks outer-boundary (sub)elements within a bodyelement.
// Returns a list of bodyelements.
//
BodyElementList*
BodyElement3D::findOuterBoundary()
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

  // ... otherwise we select outer-boundaries from subElements-list.
  IdList* se_list = model->getCoveringElementList(this->tag);
  IdList::iterator pos = se_list->begin();
  while (pos != se_list->end()) {
    int se_id = *pos++;
    if (se_id < 0)
      se_id = -1 * se_id;

    BodyElement3D* se = (BodyElement3D*)model->getFaceById(se_id);

    if ( se->getStatus() & BE_OUTER )
      oe_list->push_back(se);
  }
  return oe_list;
}


BodyElement*
BodyElement3D::getSubElement(int index)
{
  if ( modelData == NULL || modelData->nofSubElements == 0 )
    return NULL;

  if ( index < 0 || index >= modelData->nofSubElements )
    return NULL;

  return model->getBodyElementById((*modelData->subElementIds)[index]);
}


void
BodyElement3D::init(char* be_name)
{
  model->addModelObject(this, OT_FACE);

  initName(be_name);

  modelData = new BodyElementModelData;
  meshData = new BodyElementMeshData;
}


void
BodyElement3D::initClass(Model* model)
{
  BodyElement3D::last_tag = 0;
}


// Method checks if bodyelement is Ok.
bool
BodyElement3D::isOk()
{
  //  if ( ! (status & BE_OUTER_CHECKED) )
  //    checkOuterBoundary();
  return ( 0 != modelData->status & BE_OUTER_CHECKED );
}


// NOTE: Not working!!!
matchType
BodyElement3D::matchToLinear(BodyElement* be2, BodyElement*& common)
{
  return MATCH_NONE;
}


// NOTE: Not working!!!
matchType
BodyElement3D::matchToNurbs(BodyElement* be2, BodyElement*& common)
{
  return MATCH_NONE;
}

