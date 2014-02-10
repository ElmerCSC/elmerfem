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
Module:     ecif_bodyelement_loop.cpp
Language:   C++
Date:       15.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_bodyElementLoop.h"
#include "ecif_bodyElement.h"
#include "ecif_boundbox.h"
#include "ecif_model.h"
#include "ecif_userinterface.h"

//Initialize static class variables.
int BodyElementLoop::last_tag = 0;


BodyElementLoop::BodyElementLoop()
{
  tag = ++last_tag;

  init();

  nofElements = 0;
  elementIds = NULL;
  elementTags = NULL;
  tplgType = CLOSED_LOOP;
  elementType = OT_BOUNDARY;
}


#if 0
BodyElementLoop::BodyElementLoop(int nof_elems, int* elem_tags, bool is_open, objectType etype)
{
  tag = ++last_tag;

  init();

  nofElements = nof_elems;
  elementIds = new int[nof_elems];
  elementTags = new int[nof_elems];

  for (int i = 0; i < nofElements; i++) {
    elementTags[i] = elem_tags[i];
    elementIds[i] = NO_INDEX;
  }

  if ( is_open ) {
    tplgType = OPEN_LOOP;
  }

  elementType = etype;
}
#endif


// Construct with predefined tag
//
BodyElementLoop::BodyElementLoop(int ext_tag, int nof_elems, int* elem_tags, bool is_open, objectType etype)
{
  tag = ext_tag;

  if (last_tag < tag)
    last_tag = tag;

  init();

  nofElements = nof_elems;
  elementIds = new int[nofElements];
  elementTags = new int[nofElements];

  for (int i = 0; i < nofElements; i++) {
    elementIds[i] = NO_INDEX;
    elementTags[i] = elem_tags[i];
  }

  if ( is_open ) {
    tplgType = OPEN_LOOP;
  }

  elementType = etype;
}


BodyElementLoop::BodyElementLoop(IdList* elem_ids, bool is_open, objectType etype)
{
  tag = ++last_tag;

  init();

  nofElements = elem_ids->size();

  elementIds = new int[nofElements];
  elementTags = new int[nofElements];

  int pos =0;
  for (IdList::iterator itr = elem_ids->begin(); itr != elem_ids->end(); itr++) {
    elementIds[pos] = *itr;

    BodyElement* be = model->getBodyElementById(*itr);

    if ( be != NULL ) {
      elementTags[pos] = be->Tag();
    } else {
      elementTags[pos] = NO_INDEX;
    }

    pos++;
  }

  if ( is_open ) {
    tplgType = OPEN_LOOP;
  }

  elementType = etype;
}


BodyElementLoop::~BodyElementLoop()
{
  delete[] elementIds;
  delete[] elementTags;
  delete boundbox;
  delete minimumbox;
}


bool
BodyElementLoop::check()
{
  enum geomError ge_rc = GE_NONE;

  // No need to check this, it is just a technical
  // container for elements
  //
  if ( type == LOGICAL_LOOP ) return true;

  update();

  // If ok
  if ( check2D(ge_rc) ) return true;

  UserInterface* gui = (UserInterface*)model->getGui();
  strstream strm1, strm2;

  switch (ge_rc) {
  case GE_NOT_CLOSED:
      strm1 << "***ERROR in geometry for the edge loop " << tag << ":" << ends;
      strm2 << "---Loop not closed!" << ends;
      gui->showMsg(strm1.str());
      gui->showMsg(strm2.str(), 1);
  case GE_NOT_CONTINUOUS:
      strm1 << "***ERROR in geometry for the edge loop " << tag << ":" << ends;
      strm2 << "---Loop not continuous!" << ends;
      gui->showMsg(strm1.str());
      gui->showMsg(strm2.str(), 1);
  }

  return false;
}


// Method checks that edge-loop is ccw oriented.
bool
BodyElementLoop::check2D(enum geomError& rc)
{
  rc = GE_NONE;

  if ( !checkLoopSequence2D(rc) ) {
    return false;
  }

  // If loop is open,  no need to check furtheer!
  //
  if ( tplgType == OPEN_LOOP ) return true;

  // We select a point which is outside the loop and from where
  // a ray which is parallel to x-axis won't hit any vertices
  // (and is not along some edge) of the loop.
  // Idea is to find the first bodyelement which intersect
  // with the ray, and from that info to conclude what is the
  // counter-clock-wise order for the loop.
  RangeVector rv;
  minimumbox->getRangeVector(rv);

  // Select a starting point whose: x-value is certainly smaller
  // than body's min-x-value, y-value is such that no vertices are hit.
  GcPoint* startp = new GcPoint(rv[0] -1, (rv[2] + rv[3]) / 2, 0);

  // Find the element which is hit first (looked along x-axis).
  double min_xval = NSVD;
  RayHit* isection;
  BodyElement* min_elm = NULL;
  bool negative_on_left = false;

  int index = 0;
  while (true) {

    BodyElement* be = getElement(index++);

    if (be==NULL) break;

    bool neg_on_left;

    isection = be->isectionsWithXdir(startp, neg_on_left);

    if (isection) {

      // Did we find a smaller x-value.
      if ( min_xval == NSVD ||
           1 == compare(min_xval, isection->min_value)
         ) {
        min_xval = isection->min_value;
        negative_on_left = neg_on_left;
        min_elm = be;
      }
    }
  }

  // We still have to check how element is oriented in the loop.
  // Without that information we can't decide which is cw- or ccw-direction.
  bool along = inLoopDirection(min_elm);

  // Default is ccw.
  int result = 1;

  // We have cw-order if:
  // element is in 'cw-order' (negative_on_left) and is along the loop order  OR
  // element is in 'ccw-order' but is not along loop order.
  if ( ( negative_on_left && along ) || ( !negative_on_left && !along) ) {
    result = -1;
  }

  if ( result == -1 ) {
    reverseLoopOrder2D();
  }

  return true;
}


// Check that elements are in correct sequence
// Correct sequence if needed
// Correct also element id signs if elements are in reverse
// order with regard to the loop order, ie, if sign is inittially
// positive althoug the second vertex should be the first
//
bool
BodyElementLoop::checkLoopSequence2D(enum geomError& rc)
{
  int i;

  rc = GE_NONE;

  if ( nofElements == 0 ) {
    rc = GE_EMPTY;
    return false;
  }

  if ( nofElements == 1 && tplgType != OPEN_LOOP ) {
    BodyElement* be = getElement(0);
    if ( !(be->isClosedU()) ) {
      rc = GE_NOT_CLOSED;
      return false;
    } else {
      return true;
    }
  }

  BodyElement* be;
  BodyElement* be2;
  BodyElement* vrtx11;
  BodyElement* vrtx12;
  BodyElement* vrtx21;
  BodyElement* vrtx22;

  int eid = elementIds[0];
  int sign = (eid < 0)?-1:1;
  be = model->getBodyElementById(sign * eid);

  // Pick start vertex = end of first
  if ( sign == 1 ) {
    vrtx11 = be->getFirstSubElement();
    vrtx12 = be->getLastSubElement();
  } else {
    vrtx11 = be->getLastSubElement();
    vrtx12 = be->getFirstSubElement();
  }

  int start = 1;

  // Continue with other elements
  while ( start < nofElements ) {

    bool next_found = false;

    for (i = start; i < nofElements; i++) {

      int eid2= elementIds[i];
      int sign2 = (eid2 < 0)?-1:1;
      be2 = model->getBodyElementById(sign2 * eid2);

      // Pick vertex to check = first of next
      if ( sign2 == 1 ) {
        vrtx21 = be2->getFirstSubElement();
        vrtx22 = be2->getLastSubElement();
      } else {
        vrtx21 = be2->getLastSubElement();
        vrtx22 = be2->getFirstSubElement();
      }

      // Next found
      if ( vrtx12 == vrtx21 || vrtx12 == vrtx22 ) {

        // If next is not in the correct palce, swap
        if ( start != i ) {
          int tmp = elementIds[start];
          elementIds[start] = elementIds[i];
          elementIds[i] = tmp;
        }

        // Check element2 sign and pick new start
        //
        if ( vrtx12 == vrtx21 ) {
          if ( sign2 == -1 ) {
            elementIds[start] = -1 * eid2 * sign2;
          }
          vrtx12 = vrtx22;

        } else if ( vrtx12 == vrtx22 ) {
          if ( sign2 == 1 ) {
            elementIds[start] = -1 * eid2 * sign2;
          }
          vrtx12 = vrtx21;
        }

        next_found = true;
        start++;
        break;
      }
    }

    if ( !next_found ) {
      rc = GE_NOT_CONTINUOUS;
      return false;
    }

  }

  // Check that loop is closed if it should be closed!
  //
  if ( tplgType == CLOSED_LOOP && vrtx11 != vrtx12 ) {
    rc = GE_NOT_CLOSED;
    return false;
  }

  return  true;
}



bool
BodyElementLoop::check3D(enum geomError& rc)
{
  rc = GE_NONE;
  return true;
}


// Replace devided elements with their covering
// elements
bool
BodyElementLoop::checkElements()
{
  IdList ce_list;
  IdArray new_ids;


  int nof_elements = 0;
  int pos = -1;
  int index = 0;

  while (true) {

    BodyElement* be = getElement(index++);

    if (be==NULL) break;

    pos++;

    // Not devided, keep the element as it is
    if ( !(BE_DEVIDED & be->getStatus()) ) {
      new_ids.push_back(elementIds[pos]);
      nof_elements++;
      continue;
    }

    // Devided, skip the element and add the covering
    // elements instead
    // NOTE: Add them in right order!
    ce_list.clear();

    model->getCoveringElementList(be->Id(), ce_list);

    int direction = (elementIds[pos] < 0)?-1:1;

    while ( !ce_list.empty() ) {
      int ce_id;

      // Same order, same signs
      if ( direction == 1 ) {
        ce_id = ce_list.front();
        ce_list.pop_front();

      // Reverse order and opposite signs
      } else {
        ce_id = -1 * ce_list.back();
        ce_list.pop_back();
      }

      new_ids.push_back(ce_id);
      nof_elements++;
    }
  }

  nofElements = nof_elements;

  // Find minimum tag (we start loops from minimum tag)
  int i;
  int min_tag = NO_INDEX;
  int min_pos = 0;
  for (i = 0; i < nofElements; i++) {

    int el_id = new_ids[i];
    int be_id = (el_id < 0)?-el_id:el_id;

    int be_tag = NO_INDEX;

    BodyElement* be = model->getBodyElementById(be_id);

    if ( be != NULL ) {
      be_tag = be->Tag();
    }

    if ( min_tag == NO_INDEX  || be_tag < min_tag ) {
      min_tag = be_tag;
      min_pos = i;
    }
  }

  // Copy new ids
  delete[] elementIds;
  elementIds = new int[nofElements];

  for (i = 0; i < nofElements; i++) {
    elementIds[i] = new_ids[(min_pos + i) % nofElements];
  }

  return true;
}


bool
BodyElementLoop::convertTags2Ids()
{
  UserInterface* gui = (UserInterface*)model->getGui();
  strstream strm1, strm2;

  for (int i = 0; i < nofElements; i++) {

    int be_tag = elementTags[i];

    int sign = (be_tag < 0)?-1:1;

    BodyElement* be = NULL;

    switch (elementType) {
    case OT_BOUNDARY:
      be = model->getBoundaryByTag(sign * be_tag);
      break;
    case OT_FACE:
      be = model->getFaceByTag(sign * be_tag);
      break;
    case OT_EDGE:
      be = model->getEdgeByTag(sign * be_tag);
      break;
    case OT_VERTEX:
      be = model->getVertexByTag(sign * be_tag);
      break;
    }

    if ( be != NULL ) {

      // NOTE: Be ids and tags are stored unsigned
      //
      elementIds[i] = sign * be->Id();
      elementTags[i] = sign * be_tag;

    } else {
      elementIds[i] = NO_INDEX;
      strm1 << "***ERROR in definition for the edge loop " << tag << ":" << ends;
      strm2 << "---Edge " << sign * be_tag << " not defined in the input!" << ends;
      gui->showMsg(strm1.str());
      gui->showMsg(strm2.str(), 1);
      return false;
    }

  }

  return true;

}


BodyElement*
BodyElementLoop::getBodyElementById(int be_id)
{
  return model->getBodyElementById(be_id);
}


// Return signed) element id
//
int
BodyElementLoop::getDirectedElementId(int index)
{
  if ( index < 0 || index >= nofElements )
    return 0;

  return elementIds[index];
}


// Get body element by index
//
BodyElement*
BodyElementLoop::getElement(int index)
{
  if ( index < 0 || index >= nofElements ) {
    return NULL;
  }

  int elem_id = elementIds[index];

  elem_id = (elem_id < 0)?-elem_id:elem_id;

  return getBodyElementById(elem_id);
}


// Get directed body element by index
//
void
BodyElementLoop::getElement(int index, DirectedBodyElement& dbe)
{
  int direction = 1;

  if ( index < 0 || index >= nofElements ) {
    dbe.element = NULL;

  } else {

    int elem_id = elementIds[index];

    if (elem_id < 0) {
      direction = -1;
      elem_id = -1 * elem_id;
    }

    dbe.element = getBodyElementById(elem_id);
  }

  dbe.direction = direction;
}


// Return (real ie. unsigned) element id
//
int
BodyElementLoop::getElementId(int index)
{
  if ( index < 0 || index >= nofElements )
    return NO_INDEX;

  int eid = elementIds[index];

  return (eid < 0 )?-eid:eid;
}

// Return (real ie. unsigned) element tag
//
int
BodyElementLoop::getElementTag(int index)
{
  if ( index < 0 || index >= nofElements )
    return NO_INDEX;

  int etg = elementTags[index];

  return (etg < 0 )?-etg:etg;
}




BodyElement*
BodyElementLoop::getLastElement()
{
  return getElement(nofElements-1);
}


void
BodyElementLoop::getLastElement(DirectedBodyElement& dbe)
{
  getElement(nofElements-1, dbe);
}


// Get nof elements to be output into mif-file for
// the body element loop
//
// NOTE: Normally this the number of elements in the bel, but
// in multigeometry element case some bel elements may 'run out'
// geometries before the first element (which defines the maximum
// number of geometries to be taken into account) and the count
// would be zero!
//
int
BodyElementLoop::getNofElementMifTags(int gmtr_index)
{
  if ( gmtr_index == NO_INDEX ) {
    return nofElements;
  }

  int count = 0;

  for (int i = 0; i < nofElements; i++) {

    BodyElement* be = getElement(i);

    if ( be == NULL ) return 0;

    if ( NO_INDEX != be->getMifGeometryTag(gmtr_index) ) {
      count += 1;
    }
  }

  return count;
}


// Total number of mif-file layer loops which must be created
// for one body element loop
//
// NOTE: If bel starts with a single-geometry element (when we have one
// layer and gmtr_index=NO_INDEX) a following muli-geometry elements creat
// one mif-loop per each geoemtry component they have
//
// index: the multi-geometry (mif-layer index) index
//
int
BodyElementLoop::getNofMifLoops(int gmtr_index)
{
  int count = 0;

  BodyElement* be = getElement(0);

  if ( be == NULL ) return 0;

  if ( gmtr_index == NO_INDEX ) {
    count += be->getNofMifGeometries();
  } else {
    if ( NO_INDEX != be->getMifGeometryTag(gmtr_index) ) {
      count = 1;
    }
  }

  return count;
}


void
BodyElementLoop::init()
{
  model->addModelObject(this, OT_ELEMENT_LOOP);

  boundbox = new BoundBox;
  minimumbox = new BoundBox(MAX_RANGE, MAX_RANGE);
  parentIds[0] = parentIds[1] = NO_INDEX;

  type = NORMAL_LOOP;
}


// Method checks if bodyelement is directed postively or negatively
// in the elementloop.
// True means according to loop-table order.
bool
BodyElementLoop::inLoopDirection(BodyElement* be)
{
  for (int i = 0; i < nofElements; i++) {

    int eid = elementIds[i];
    int sign = (eid < 0)?-1:1;

    BodyElement* my_be = model->getBodyElementById(sign * eid);

    if ( my_be == be ) {

      if ( sign < 0)
        return false;
      else
        return true;
    }

  }

  return false;
}


void
BodyElementLoop::initClass(Model* mdl)
{
  BodyElementLoop::last_tag = 0;
}


void
BodyElementLoop::markActiveObjects()
{
  for (int i = 0; i < nofElements; i++) {
    int be_id = elementIds[i];
    int sign = (be_id < 0)?-1:1;
    model->markObjectActive(sign * be_id);
  }

}


// Get all element mif-tags for element loop
//
// NOTE: For a multi-geometry case it is possible that some elements
// 'run outof' geometry components before the first multi-geometry
// element and then we do not ouput any tags!
//
ostream&
BodyElementLoop::outputDirectedElementMifTags(ostream& out, bool as_open,
                                              int indent_size, int gmtr_index, int max_per_line)
{
  int nof_elements = getNofElementMifTags(gmtr_index);
  int counter = 0;
  int i;

  if ( nof_elements == 0 ) return out;

  for (i = 0; i < nof_elements; i++) {

    int eid = elementIds[i];
    int sign = (eid < 0)?-1:1;

    BodyElement* be = model->getBodyElementById(sign * eid);
    int mif_tag = NO_INDEX;

    // Normal single-geometry call
    if ( gmtr_index == NO_INDEX ) {
      mif_tag = be->getMifGeometryTag(0);
    // Multi-geometry call, geometry index given
    } else {
      mif_tag = be->getMifGeometryTag(gmtr_index);
    }

    // If no geometry (component) any more available
    //
    if ( mif_tag == NO_INDEX ) continue;

    out << sign * mif_tag;
    counter++;

    if ( i < (nofElements - 1) ) {
      if ( indent_size > 0 && counter >= max_per_line ) {
        out << endl;
        indent(out, indent_size);
        counter = 0;
      } else {
        out << " ";
      }
    }
  }


  // NOTE: Open loops are printed like:
  // Edges: 6   1 2 3 -3 -2 -1
  // so this will output the extra '-3 -2 -1'
  //
  if ( as_open ) {

    out << ' ';

    for (i = nof_elements-1; i >= 0; i--) {

      int eid = elementIds[i];
      int sign = (eid < 0)?1:-1;

      BodyElement* be = model->getBodyElementById(-1 * sign * eid);
      int mif_tag = NO_INDEX;

      // Normal single-geometry call
      if ( gmtr_index == NO_INDEX ) {
        mif_tag = be->getMifGeometryTag(0);
      // Multi-geometry call, geometry index given
      } else {
        mif_tag = be->getMifGeometryTag(gmtr_index);
      }

      // If no geometry (component) any more available
      //
      if ( mif_tag == NO_INDEX ) continue;

      out << sign * mif_tag;
      counter++;

      if ( i > 0 ) {
        if ( indent_size > 0 && counter >= max_per_line ) {
          out << endl;
          indent(out, indent_size);
          counter = 0;
        } else {
          out << " ";
        }
      }
    }
  }

  return out;
}


ostream&
BodyElementLoop::outputDirectedElementTags(ostream& out, int indent_size, int max_per_line)
{
  int counter = 0;

  for (int i = 0; i < nofElements; i++) {

    int eid = elementIds[i];

    int sign = (eid < 0)?-1:1;

    BodyElement* be = model->getBodyElementById(sign * eid);

    out << sign * be->Tag();

    counter++;

    if ( i < (nofElements - 1) ) {

      if ( indent_size > 0 && counter >= max_per_line ) {
        out << endl;
        indent(out, indent_size);
        counter = 0;

      } else {
        out << " ";
      }
    }
  }

  return out;
}


ostream&
BodyElementLoop::output_emf(ostream& out, short indent_size, short indent_level)
{
  // NOTE: This is just a technical 'container' for
  // body elements, it is not output to emf
  //
  if ( type == LOGICAL_LOOP ) return out;

  short is = indent_size;
  short il = indent_level;

  // Header
  LibFront::output_scalar(out, is, il, "Element Loop", NULL, tag);

  LibFront::output_string(out, is, 1 + il, "Elements", false) << ' ';

  outputDirectedElementTags(out) << endl;

  return out;
}


void
BodyElementLoop::reverseLoopOrder2D()
{
  int i, j;

  for (i = 0, j = nofElements -1; i <= j; i++, j--) {

    int tmp = elementIds[i];

    elementIds[i] = -1 * elementIds[j];
    elementIds[j] = -1 * tmp;
  }
}


void
BodyElementLoop::setParentId(int parent_index, int parent_id)
{
  if ( parent_index < 0 | parent_index > 1 ) {
    return;
  }

  parentIds[parent_index] = parent_id;
}


void
BodyElementLoop::swapElements(IdArray* ids1, IdArray* ids2, IdArray* relative_dirs)
{
  int sign, eid, id1, id2, dir;

  int count = ids1->size();

  for (int i = 0; i < nofElements; i++) {

    eid = elementIds[i];

    sign = (eid < 0)?-1:1;

    for (int j = 0; j < count; j++) {

      id1 = (*ids1)[j];
      id2 = (*ids2)[j];

      if ( relative_dirs != NULL ) {
        dir = (*relative_dirs)[j];

      } else {
        dir = 1;
      }

      if ( id1 == (sign * eid) ) {
        elementIds[i] = sign * dir * id2;
      }
    }
  }
}


void
BodyElementLoop::update()
{
  RangeVector rv;

  int index = 0;
  while (true) {
    BodyElement* be = getElement(index++);
    if (be==NULL) break;
    be->getRangeVector(rv);
    updateBoundBox(rv);
    updateMinimumBox(rv);
  }
}


void
BodyElementLoop::updateBoundBox(RangeVector rv)
{
  boundbox->extendByRange(rv);
}


void
BodyElementLoop::updateMinimumBox(RangeVector rv)
{
  minimumbox->restrictByRange(rv);
}

