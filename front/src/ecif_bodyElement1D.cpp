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
Module:     ecif_bodyelement1D.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_body.h"
#include "ecif_bodyElement1D.h"
#include "ecif_control.h"
#include "ecif_geometry.h"
#include "ecif_parameter.h"
#include "ecif_renderer.h"

extern Control* theControlCenter;


//Initialize static class variables.
int BodyElement1D::last_tag = 0;

// Constructors.

BodyElement1D::BodyElement1D()
  : BodyElement()
{
  init();

  update();
}


BodyElement1D::BodyElement1D(int int_tag)
  : BodyElement(int_tag)
{
  checkLastTag();
  init();

  update();
}


//-Create element from point
BodyElement1D::BodyElement1D(int int_tag, GcPoint* p)
  : BodyElement(int_tag)
{
  init();

  ptrGmtr = new GcPoint(p);

  update();
}


//-Create element from point
BodyElement1D::BodyElement1D(GcPoint* p)
  : BodyElement()
{
  tag = newTag();
  init();

  ptrGmtr = new GcPoint(p);

  update();
}


BodyElement1D::BodyElement1D(ecif_Element_X& tx)
  : BodyElement(tx)
{

  if ( tag == NO_INDEX ) {
    tag = newTag();
  }

  checkLastTag();
  init();

  ptrGmtr = new GcPoint(tx.components[0]->geometry.vertex->point);

  update();
}


BodyElement1D::BodyElement1D(ecif_Vertex_X& tx)
  : BodyElement(tx)
{
  checkLastTag();
  init();

  ptrGmtr = new GcPoint(tx.point);

  update();
}


BodyElement1D::~BodyElement1D()
{
  delete ptrGmtr;
  delete gridHData;
}


void BodyElement1D::draw(Renderer* renderer, flagName geometry_type, int body_id)
{
  if ( drawMode == DM_HIDDEN && drawState != DS_SELECTED ) return;

  // We draw vertices in 3D only when they are selected
  if ( model->getDimension() == ECIF_3D &&
       drawState != DS_SELECTED
     ) {
    return;
  }

  renderer->name_save(id);
  ptrGmtr->draw(renderer, drawMode, drawState);
  renderer->name_delete(id);
}


// Method creates a new 1D vertex-type from a vertex.
BodyElement*
BodyElement1D::createElement(int nof_vertices, int* vertex_ids, ecif_geometryType gt)
{
  return NULL;
}


int
BodyElement1D::findMeshBorderNodes(int buf_size, int* ids_buffer)
{
  if ( nodeId != NO_INDEX ) {
    // In 1D this is really easy!
    ids_buffer[0] = nodeId;
    return 1;

  } else {
    return 0;
  }
}


int
BodyElement1D::getMeshElementId(int index)
{
  if ( index != 0 ) {
    return NO_INDEX;

  } else {
    return nodeId;
  }
}


void
BodyElement1D::init(char* be_name)
{
  model->addModelObject(this, OT_VERTEX);

  nodeId = NO_INDEX;

  if (be_name != NULL) {
    initName(be_name);
  }

  nofCurrentMeshHSources = 0;
}


void
BodyElement1D::initClass(Model* model)
{
  BodyElement1D::last_tag = 0;
}


void
BodyElement1D::initLabelData()
{
  if ( ptrGmtr == NULL )
    return;

  if ( labelData == NULL ) {
    labelData = new LabelData;
    labelData->label = NULL;
  }

  // Make an id-string.
  ostrstream id_str;
  id_str << 'V' << Tag() << ends;
  update_dyna_string(labelData->label, id_str.str());

  //---Get 'nice' position for the id-label
  ptrGmtr->getLabelPoint(labelData->position);
}

#if 0
// Do nothing! (vertices are output as table)
ostream&
BodyElement1D::output_emf(ostream& out, short indent_size, short indent_level, bool isOnSymmAxis)
{
  return out;
}
#endif


// Output vertex for Elmer Mesh input file (mif-file)
//
ostream&
BodyElement1D::output_mif(ostream& out)
{
  Point3 p;

  ((GcPoint*)ptrGmtr)->getPoint(p);

  int mesh_index = model->getCurrentMeshIndex();

  //---Tag and boundary tag
  out << "NodeId: " << tag << " " << boundaryTag;

  //---Possible mesh-density value
  char type;
  int nof_values;
  double values[4];

  if ( getMeshDensityValues(mesh_index, type, nof_values, values) &&
       nof_values == 1 &&
       values[0] > 0
     ) {
    out << "  " << type << ": " << values[0];
  }

  //---Coordinates (x,y)
  out << "  " << p[0] << " " << p[1];

  out << endl;

  return out;
}

