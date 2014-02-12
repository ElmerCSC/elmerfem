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
Module:     ecif_bodyElementLoop.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   A Base class for bodyelement loops.
            They define a body by curves (2D) or surfaces (3D)

************************************************************************/

#ifndef _ECIF_BODYELEMENT_LOOP_
#define _ECIF_BODYELEMENT_LOOP_

#include "ecif_def.h"
#include "ecif_def_stl.h"
#include "ecif_modelObject.h"


class BodyElementLoop : public ModelObject
{

public:
  BodyElementLoop();
  BodyElementLoop(int bel_tag, int nof_elements, int* element_ids, bool is_open, objectType etype);
  //BodyElementLoop(int nof_elements, int* element_ids, bool is_open, objectType etype);
  BodyElementLoop(IdList* element_ids, bool is_open, objectType etype);
  ~BodyElementLoop();
  bool check();
  bool checkElements();
  bool convertTags2Ids();
  int getDirectedElementId(int index);
  int getElementId(int index);
  const int* getElementIds() { return elementIds; }
  int getElementTag(int index);
  const int* getElementTags() { return elementTags; }
  BodyElement* getElement(int index);
  void getElement(int index, DirectedBodyElement& dbe);
  BodyElement* getLastElement();
  void getLastElement(DirectedBodyElement& dbe);
  void getParentIds(int& parent1_id, int& parent2_id) { parent1_id = parentIds[0]; parent2_id = parentIds[1]; }
  int getNofElements() {return nofElements;}
  int getNofElementMifTags(int gmtr_index = NO_INDEX);
  int getNofMifLoops(int gmtr_index);
  elementLoopTplgType getTplgType() { return tplgType; }
  elementLoopType getType() { return type; }
  static void initClass(Model* mdl);
  bool inLoopDirection(BodyElement* be);
  void markActiveObjects();
  ostream& outputDirectedElementTags(ostream& out, int indent_size = 0, int max_per_line = 20);
  ostream& outputDirectedElementMifTags(ostream& out, bool as_open, int indent_size = 0, int gmtr_index = NO_INDEX, int max_per_line = 20);
  ostream& output_emf(ostream& out, short indent_size, short indent_level);
  void setParentId(int parent_index, int parent_id);
  void setParentIds(int parent1_id, int parent2_id) { parentIds[0] = parent1_id; parentIds[1] = parent2_id; }
  void setType(elementLoopType value) { type = value; }
  void setTplgType(elementLoopTplgType value) { tplgType = value; }
  void swapElements(IdArray* ids1, IdArray* ids2, IdArray* relative_dirs = NULL);

protected:
  bool check2D(enum geomError& rc);
  bool check3D(enum geomError& rc);
  bool checkLoopSequence2D(enum geomError& rc);
  BodyElement* getBodyElementById(int be_id);
  void init();
  void reverseLoopOrder2D();
  void update();
  void updateBoundBox(RangeVector rv);
  void updateMinimumBox(RangeVector rv);

  static int last_tag;
  int nofElements;
  int* elementIds;
  int* elementTags;
  elementLoopTplgType tplgType;
  elementLoopType type;
  objectType elementType;
  BoundBox* boundbox;
  BoundBox* minimumbox;
  int parentIds[2];  // Parent body ids
};

#endif

