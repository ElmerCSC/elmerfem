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
Module:     ecif_bodyElementGroup.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   A Base class for bodyelemets.
            A bodyelement is a separable entity in a body.
            They can be edges, surfaces etc. depending on the dimension.

************************************************************************/

#ifndef _ECIF_BODYELEMENT_GROUP_
#define _ECIF_BODYELEMENT_GROUP_

#include "ecif_def.h"
#include "ecif_def_stl.h"
#include "ecif_boundbox.h"
#include "ecif_geometry.h"
#include "ecif_modelObject.h"


// Element group for boundary conditions etc.
//
class BodyElementGroup : public ModelObject
{
friend class Control;
friend class Model;
public:
  BodyElementGroup();
  BodyElementGroup(ecif_ElementGroup_X& trx_bg, enum objectType bndr_type);
  virtual ~BodyElementGroup();
  int addElement(int bn_id);
  bool check();
  const BodyElement* getElement(int index);
  int getElementId(int index);
  const int* getElementIds() { return elementIds; }
  int getElementTag(int index);
  const int* getElementTags() { return elementTags; }
  int getBoundaryConditionId() { return boundaryConditionId;}
  int getBoundaryParameterId() { return boundaryParameterId;}
  enum objectType getElementType() const { return elemType;}
  enum elementGroupType getGroupType() const { return groupType;}
  int getNofElements() { return nofElements; }
  int getParentId(short parent) const;
  int getParentLayer(short parent) const;
  int getParentTag(short parent) const;
  bool hasElement(int be_id);
  void initName();
  static void initClass(Model* model);
  bool isExplicit() { return groupType == EXPLICIT_GROUP; }
  bool isImplicit() { return groupType == IMPLICIT_GROUP; }
  bool isVirtual() { return groupType == VIRTUAL_GROUP; }
  virtual ostream& output_emf(ostream& out, short indent_size, short indent_level);
  int removeElement(int bn_id);
  void setElementIdsAndTags(int nof_ids, int* be_ids, int* be_tgs);
  void setBoundaryConditionId(int bc_id);
  void setBoundaryParameterId(int pid);
  void setElementType(enum objectType et);
  void setType(enum elementGroupType gt);
  void setParentId(short parent_nbr, int parent_id);
  void setParentLayer(short parent_nbr, int layer);
  void setParentTag(short parent_nbr, int parent_tag);

protected:
  void init();

  static int last_tag;
  int* elementIds;
  int* elementTags;
  int boundaryConditionId;
  int boundaryParameterId;
  int nofElements;
  enum elementGroupType groupType;  // REAL_GROUP etc.
  enum objectType elemType;  // OT_EDGE etc.
  int parent1Id;
  int parent2Id;
  int parent1Layer;
  int parent2Layer;
  int parent1Tag;
  int parent2Tag;

} ;

#endif

