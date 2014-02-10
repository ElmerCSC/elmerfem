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
Module:     ecif_boundaryCondition.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   A Base class for boundary conditions.

************************************************************************/

#ifndef _ECIF_BOUNDARY_CONDITION_
#define _ECIF_BOUNDARY_CONDITION_

#include "ecif_parameter.h"


// ****** Boundary Condition class ******
class BoundaryCondition : public Parameter{
public:
  BoundaryCondition();
  BoundaryCondition(int p_id);
  BoundaryCondition(int cid, int parent_id, char* data_string, char* param_name);
  ~BoundaryCondition();
  int getLastId() {return last_id;}
  void setLastId(int lid) {last_id = lid;}
  const char* getGuiName() { return "Constraint"; }
  const char* getArrayName() { return "BoundaryCondition"; }
  const char* getEmfName() { return "Boundary Condition"; }
  const char* getSifName() { return SIF_BOUNDARY_CONDITION; }
  ecif_parameterType getParameterType() { return ECIF_BOUNDARY_CONDITION; }
  bool hasZeroVelocity();
  static void initClass(Model* model);
  virtual ostream& output_sif(ostream& out, short indent_size, short indent_level, SifOutputControl& soc);
  void setName(char* param_name);
  void updateParentId();
  void updateParentInfo(int parent_id);
  void updateTargetTags();

protected:
  static int last_id;
  static Model* model;
  int nofTargetBoundaries;
  int* targetBoundaryTags;
  int targetBodyTag;
  const ModelObject* getParentEmfObject();
};


#endif
