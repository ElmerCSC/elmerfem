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
Module:     ecif_bodyElementLoop.cpp
Language:   C++
Date:       01.01.038
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_bodyElementGroup.h"
#include "ecif_bodyElement.h"
#include "ecif_model.h"
#include "ecif_userinterface.h"

//Initialize static class variables.
int BodyElementGroup::last_tag = 0;


// BodyElementGroup  class
// ====================

BodyElementGroup::BodyElementGroup()
{
  tag = ++last_tag;
  init();
  //initName();
}


BodyElementGroup::BodyElementGroup(ecif_ElementGroup_X& tx, enum objectType elem_tp)
{
  if ( tx.tag == NO_INDEX ) {
    tag = ++last_tag;
  } else {
    tag = tx.tag;
  }

  if (last_tag < tag) {
    last_tag = tag;
  }

  init();

  elemType = elem_tp;

  if ( tx.is_virtual ) {
    groupType = VIRTUAL_GROUP;
  }

  update_dyna_string(name, tx.name);

  // If name not defined
  if ( name == NULL || name[0] == '\0' ) {
    strstream strm;
    strm << " Bg-" << tag << ends;
    update_dyna_string(name, strm.str());
  }

  nofElements = tx.nof_elements;

  elementTags = new int[nofElements];
  elementIds = new int[nofElements];

  for (int i = 0; i < nofElements; i++) {
    elementTags[i] = tx.element_tags[i];
    elementIds[i] = NO_INDEX;
  }

  boundaryConditionId = tx.boundary_cond_id;
  boundaryParameterId = tx.boundary_param_id;

}


BodyElementGroup::~BodyElementGroup()
{
  delete[] elementIds;
  delete[] elementTags;
}


int
BodyElementGroup::addElement(int be_id)
{
  // Element exists
  if ( hasElement(be_id) ) return nofElements;

  BodyElement* be = model->getBodyElementById(be_id);

  // Add element

  // Allocate new arries
  // NOTE: Do not delete these!!!
  int* tmp_ids = new int[nofElements+1];
  int* tmp_tags = new int[nofElements+1];

  // Copy old data
  for (int i = 0; i < nofElements; i++) {
    tmp_ids[i] = elementIds[i];
    tmp_tags[i] = elementTags[i];
  }

  // Delete old data
  delete[] elementIds;
  delete[] elementTags;

  // Add new info
  tmp_ids[nofElements] = be_id;
  tmp_tags[nofElements] = be->Tag();

  // Update object data
  nofElements++;
  elementIds = tmp_ids;
  elementTags = tmp_tags;

  return nofElements;
}


bool
BodyElementGroup::check()
{
  UserInterface* gui = (UserInterface*)model->getGui();
  strstream strm1, strm2;

  BodyElement* be = NULL;
  bool is_first_bndr = true;

  for (int i = 0; i < nofElements; i++) {

    int elem_tag = elementTags[i];

    switch (elemType) {
    case OT_VERTEX:
      be = model->getBodyElementByTag(OT_VERTEX, elem_tag);
      break;
    case OT_EDGE:
      be = model->getBodyElementByTag(OT_EDGE, elem_tag);
      break;
    case OT_FACE:
      be = model->getBodyElementByTag(OT_FACE, elem_tag);
      break;
    }

    if ( be != NULL ) {

      elementIds[i] = be->Id();

      // For non-virtual groups:
      // -set group info into member boundaries
      // -set parent body info
      // -check that elements are homogenous
      if ( !isVirtual() ) {

        be->setElementGroupTag(tag);
        be->setElementGroupId(id);

        // NOTE: Body info from the first element
        // (boundaries must be homogenous in this respect!)
        //
        if ( is_first_bndr ) {
          is_first_bndr = false;
          parent1Id = be->getParentId(1);
          parent2Id = be->getParentId(2);
          parent1Layer = be->getParentLayer(1);
          parent2Layer = be->getParentLayer(2);
          parent1Tag = be->getParentTag(1);
          parent2Tag = be->getParentTag(2);
        } else if ( parent1Id != be->getParentId(1) ||
                    parent2Id != be->getParentId(2)
                  ) {
            strm1 << "***ERROR in Edge Group " << tag << ends;
            strm2 << "---Edge " << be->Tag() << " is not similar to previous edges!" << ends;
            gui->showMsg(strm1.str());
            gui->showMsg(strm2.str());
            return false;
        }
      }

    // Error!
    } else {
      strm1 << "***ERROR in Edge Group " << tag << ends;
      strm2 << "---Cannot find edge " << elem_tag << ends;
      gui->showMsg(strm1.str());
      gui->showMsg(strm2.str());
      return false;
    }
  }

  return true;
}


const BodyElement*
BodyElementGroup::getElement(int index)
{
  if ( index < 0 || index >= nofElements ) return NULL;

  int be_id = elementIds[index];

  return model->getBodyElementById(be_id);
}


int
BodyElementGroup::getElementId(int index)
{
  if ( index < 0 || index >= nofElements ) return NO_INDEX;

  return elementIds[index];
}


int
BodyElementGroup::getElementTag(int index)
{
  if ( index < 0 || index >= nofElements ) return NO_INDEX;

  return elementTags[index];
}


int
BodyElementGroup::getParentId(short parent) const
{
  if (parent == 1)
    return parent1Id;

  else if (parent == 2)
    return parent2Id;

  else
    return NO_INDEX;
}


int
BodyElementGroup::getParentLayer(short parent) const
{
  if (parent == 1)
    return parent1Layer;

  else if (parent == 2)
    return parent2Layer;

  else
    return NO_INDEX;
}


int
BodyElementGroup::getParentTag(short parent) const
{
  if (parent == 1)
    return parent1Tag;

  else if (parent == 2)
    return parent2Tag;

  else
    return NO_INDEX;
}



bool
BodyElementGroup::hasElement(int bn_id)
{
  for (int i = 0; i< nofElements; i++) {

    if ( elementIds[i] == bn_id ) {
      return true;
    }
  }

  return false;
}


void
BodyElementGroup::initClass(Model* mdl)
{
  BodyElementGroup::last_tag = 0;
}


void
BodyElementGroup::init()
{
  int i;

  model->addModelObject(this, OT_ELEMENT_GROUP);

  elemType = OT_NONE;
  groupType = EXPLICIT_GROUP;

  parent1Id = NO_INDEX;
  parent2Id = NO_INDEX;
  parent1Layer = NO_INDEX;
  parent2Layer = NO_INDEX;
  parent1Tag = NO_INDEX;
  parent2Tag = NO_INDEX;

  name = NULL;

  nofElements = 0;
  elementIds = NULL;
  elementTags = NULL;

  boundaryConditionId = NO_INDEX;
  boundaryParameterId = NO_INDEX;
}


void
BodyElementGroup::initName()
{
  if ( name == NULL || name[0] == '\0' ) {
    strstream strm;
    strm << "Bg-" << tag << ends;

    update_dyna_string(name, strm.str());
  }
}



ostream&
BodyElementGroup::output_emf(ostream& out, short indent_size, short indent_level)
{
  char* QM = "\"";
  short is = indent_size;
  short il = indent_level;


  // Header: Element Group
  LibFront::output_scalar(out, is, il, EMF_ELEMENT_GROUP, NULL, tag);

  // Type (output only if virtual)
  if ( isVirtual() ) {
    LibFront::output_scalar(out, is, 1 + il, EMF_TYPE, NULL, "Virtual");
  }

  // Name
  LibFront::output_scalar(out, is, 1 + il, EMF_NAME, NULL, name);

  int count, i, nof_ids;
  int* bids;
  int* mids;
  int* gids;

  //--BodyElementGroup group type and element tags
  if ( nofElements > 0 ) {
    if ( elemType == OT_FACE ) {
      LibFront::output_vector(out, is, 1 + il, "Faces", NULL, nofElements, elementTags, false);
    } else if ( elemType == OT_EDGE ) {
      LibFront::output_vector(out, is, 1 + il, "Edges", NULL, nofElements, elementTags, false);
    } else if ( elemType == OT_VERTEX ) {
      LibFront::output_vector(out, is, 1 + il, "Vertices", NULL, nofElements, elementTags, false);
    } else {
      LibFront::output_vector(out, is, 1 + il, EMF_ELEMENTS, NULL, nofElements, elementTags, false);
    }

  }

  //--Boundary parameter id
  if ( boundaryParameterId != NO_INDEX ) {
    LibFront::output_scalar(out, is, 1 + il, EMF_BOUNDARY_PARAMETER, NULL, boundaryParameterId);
  }

  //--Boundary condition id
  if ( boundaryConditionId != NO_INDEX ) {
    LibFront::output_scalar(out, is, 1 + il, EMF_BOUNDARY_CONDITION, NULL, boundaryConditionId);
  }

  return out;
}


int
BodyElementGroup::removeElement(int be_id)
{
  // Element does not exist
  if ( !hasElement(be_id) ) return nofElements;

  // Remove one element

  int* tmp_ids = NULL;
  int* tmp_tags = NULL;

  // Allocate new arries
  // NOTE: Do not delete these!!!
  if ( nofElements > 1 ) {
    tmp_ids = new int[nofElements-1];
    tmp_tags = new int[nofElements-1];

    // Copy old data, but skip be_id
    for (int i = 0; i < nofElements; i++) {
      if ( be_id == elementIds[i] ) continue;
      tmp_ids[i] = elementIds[i];
      tmp_tags[i] = elementTags[i];
    }
  }

  // Delete old data
  delete[] elementIds;
  delete[] elementTags;

  // Update object data
  nofElements--;
  elementIds = tmp_ids;
  elementTags = tmp_tags;

  return nofElements;
}


void
BodyElementGroup::setBoundaryConditionId(int bc_id)
{
  boundaryConditionId = bc_id;
}


void
BodyElementGroup::setElementIdsAndTags(int nof_ids, int* bndr_ids, int* bndr_tags)
{
  delete[] elementIds;
  delete[] elementTags;

  nofElements = nof_ids;

  elementIds = new int[nofElements];
  elementTags = new int[nofElements];

  for (int i = 0; i < nofElements; i++) {
    elementIds[i] = bndr_ids[i];
    elementTags[i] = bndr_tags[i];
  }
}


void
BodyElementGroup::setBoundaryParameterId(int pid)
{
  boundaryParameterId = pid;
}


void
BodyElementGroup::setElementType(enum objectType et)
{
  elemType = et;
}


void
BodyElementGroup::setType(enum elementGroupType gt)
{
  groupType = gt;
}


void
BodyElementGroup::setParentId(short parent, int parent_id)
{
  if (parent == 1)
    parent1Id = parent_id;
  else if (parent == 2)
    parent2Id = parent_id;
}


void
BodyElementGroup::setParentLayer(short parent, int layer)
{
  if (parent == 1)
    parent1Layer = layer;
  else if (parent == 2)
    parent2Layer = layer;
}


void
BodyElementGroup::setParentTag(short parent, int parent_tag)
{
  if (parent == 1)
    parent1Tag = parent_tag;
  else if (parent == 2)
    parent2Tag = parent_tag;
}
