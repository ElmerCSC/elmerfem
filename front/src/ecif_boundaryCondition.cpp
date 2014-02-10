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
Module:     ecif_boundaryCondition.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Implementation
  
************************************************************************/

#include "ecif_bodyElement.h"
#include "ecif_bodyElementGroup.h"
#include "ecif_boundaryCondition.h"
#include "ecif_model.h"
#include "ecif_parameterField.h"

//Initialize static class variables.
int BoundaryCondition::last_id = 0;
Model* BoundaryCondition::model = NULL;


// Constructors
BoundaryCondition::BoundaryCondition()
{
  nofTargetBoundaries = 0;
  targetBoundaryTags = NULL;
  targetBodyTag = NO_INDEX;
}


BoundaryCondition::BoundaryCondition(int p_id) : Parameter(p_id)
{
  nofTargetBoundaries = 0;
  targetBoundaryTags = NULL;
  targetBodyTag = NO_INDEX;
}


BoundaryCondition::BoundaryCondition(int cid, int parent_id, char* data_string, char* param_name)
{
  setData(cid, parent_id, data_string, param_name);
  nofTargetBoundaries = 0;
  targetBoundaryTags = NULL;
  targetBodyTag = NO_INDEX;
}


BoundaryCondition::~BoundaryCondition()
{
  delete[] targetBoundaryTags;
}


// Parent object or emf-file info
//
const ModelObject*
BoundaryCondition::getParentEmfObject()
{
  ModelObject* obj = model->getModelObjectById(parentId);

  if ( obj == NULL ) return NULL;

  BodyElementGroup* beg = (BodyElementGroup*)obj;
  
  // For implicit groups we use the first element as the parent
  //
  if ( IMPLICIT_GROUP == beg->getGroupType() ) {
    int be_id = beg->getElementId(0);
    obj = model->getModelObjectById(be_id);
  }

  return obj;
}


bool
BoundaryCondition::hasZeroVelocity()
{
  ParameterField* pf;
  short dim1, dim2, nofVars;

  short dimension = (short)model->getDimension();
  
  const char* fields[] = {"Velocity 1", "Velocity 2", "Velocity 3"};

  // Loop all relevant velocity components
  for (short d = 0; d < dimension; d++) {

    pf = getFieldBySifName(fields[d]);

    // If component is not set 
    if ( pf == NULL || 0 == pf->getNofDataStrings() )
      return false;

    // If component is a procedure
    // NOTE: We do not trust procedures at all
    if ( pf->isProcedure() )
      return true;

    // Get constraint data dims
    pf->getDataDimension(dim1, dim2, nofVars);

    // Get constraint data package
    char** data = pf->getDataStrings();

    double value;

    // Check if velocity constraint value is zero
    // NOTE: One velocity component is scalar -->
    // no need to dim2 looping
    for (int i = 0; i < dim1; i++ ) {
      strstream strm;
      strm << data[i];
      strm >> value;
      if ( value != 0 )
        return false;
    }
  }

  // Ok, we really are a zero-velocity condition!
  return true;
}


void
BoundaryCondition::initClass(Model* mdl)
{
  BoundaryCondition::model = mdl;
  BoundaryCondition::last_id = 0;
}


// Boundary condition specific output-method for Solver input file
ostream&
BoundaryCondition::output_sif(ostream& out, short indent_size, short indent_level, SifOutputControl& soc)
{
  char QM = '\"';

  // Parameter type and id
  if (soc.outputType) {

    LibFront::output_string(out, indent_size, indent_level++, getSifName(), false);

    if (soc.outputId)
      out << ' ' << ID();

    out << endl;

  } else {
    indent_level++;
  }

  // Output parameter name
  if (soc.outputName) {
    output_sif_name(out, indent_size, indent_level, soc);
  }

  // Target body tag
  // NOTE: This is only for virtual boundary group bc:s
  //
  if ( targetBodyTag != NO_INDEX ) {
    if ( model->getSolverKeywordTypeGiven(SIF_BOUNDARY_CONDITION, "Body Id") ) {
      LibFront::output_scalar(out, indent_size, indent_level, "Body Id = ", NULL, targetBodyTag);
    } else {
      LibFront::output_scalar(out, indent_size, indent_level, "Body Id = Integer ", NULL, targetBodyTag);
    }
  }
    
  // Target boundary tags
  if ( nofTargetBoundaries > 0 ) {
    strstream strm;

    if ( model->getSolverKeywordTypeGiven(SIF_BOUNDARY_CONDITION, "Target Boundaries") ) {
      strm << "Target Boundaries(" << nofTargetBoundaries << ") =" << ends;
      LibFront::output_vector(out, indent_size, indent_level,
                              strm.str(), NULL,
                              nofTargetBoundaries, targetBoundaryTags,
                              false);

    } else {
      strm << "Target Boundaries(" << nofTargetBoundaries << ") = Integer" << ends;
      LibFront::output_vector(out, indent_size, indent_level,
                              strm.str(), NULL,
                              nofTargetBoundaries, targetBoundaryTags,
                              false);
    }
  }

  out << endl;

  // Fields
  for (short i = 0; i < nofFields; i++) {

    ParameterField* pf = fields[i];

    // Check that field exists and it should be output
    if ( !pf->isActiveInstance() ||
         pf->getNofDataStrings() == 0 ||
         ( !soc.outputAll && !pf->isSifOutputField() )
       )
      continue;

   pf->output_sif(out, indent_size, indent_level, soc.sectionName);

  }

  return out;
}


void
BoundaryCondition::setName(char* param_name)
{
  Parameter::setName(param_name, "Constraint");
}


// Set parent object id (after reading from model emf-file)
void
BoundaryCondition::updateParentId()
{
  // If parent object type is a normal boundary, replace the
  // parent object with the boundary-group object
  //
  if ( parentEmfType != OT_ELEMENT_GROUP ) {

    ModelObject* obj = model->getModelObjectByTag(parentEmfType, parentEmfTag);
    if ( obj != NULL ) {
      BodyElement* be = (BodyElement*)obj;
      parentEmfTag = be->getElementGroupTag();
      parentEmfType = OT_ELEMENT_GROUP;
    } else {
      parentEmfTag = NO_INDEX;
      parentEmfType = OT_NONE;
    }
  }

  ModelObject* obj = model->getModelObjectByTag(parentEmfType, parentEmfTag);

  if ( obj != NULL ) {
    parentId = obj->Id();
  }
}


// Update parent info from Gui
//
void
BoundaryCondition::updateParentInfo(int parent_id)
{
  // Set parent object's id
  //
  parentId = parent_id;

  BodyElementGroup* beg = model->getBodyElementGroupById(parent_id);

  if ( beg != NULL ) {

    // If parent object type is an implicit group use
    // the boundary element's tag
    if ( IMPLICIT_GROUP == beg->getGroupType() ) {

        const BodyElement* be = beg->getElement(0);
        parentEmfTag = be->Tag();
        parentEmfType = be->getObjectType();

    // Otherwise use group's own tag
    } else {
        parentEmfTag = beg->Tag();
        parentEmfType = beg->getObjectType();
    }

  } else {
    parentEmfTag = NO_INDEX;
    parentEmfType = OT_NONE;
  }

}


// Update target boundary tags for the boundary condition
//
void
BoundaryCondition::updateTargetTags()
{
  nofTargetBoundaries = 0;
  delete[] targetBoundaryTags;
  targetBoundaryTags = NULL;
  targetBodyTag = NO_INDEX;

  IdsSet be_tags;

  // Body elements
  // =============
  int index = 0;
  while (true) {

    BodyElement* be = model->getBodyElement(index++);

    if ( be == NULL ) break;

    if ( id == be->getBoundaryConditionId() ) {
      be_tags.insert(be->getBoundaryTag());

      if ( be->isBemBoundary() ) {
        int bd_tg = be->getParentTag(1);
        if ( bd_tg != NO_INDEX ) {
          targetBodyTag = bd_tg;
        }
      }
    } // Boundary's bc-id matches

  } // Loop all boundaries

  // Element groups
  // =============
  index = 0;
  while (true) {

    BodyElementGroup* beg = model->getBodyElementGroup(index++);

    if (beg==NULL) break;

    if ( id == beg->getBoundaryConditionId() ) {

      for (int i = 0; i < beg->getNofElements(); i++) {
        const BodyElement* be  = beg->getElement(i);

        if ( be != NULL ) {
          be_tags.insert(be->getBoundaryTag());
        }
      }

      if ( beg->getGroupType() == VIRTUAL_GROUP ) {
        int bd_tg = beg->getParentTag(1);
        if ( bd_tg != NO_INDEX ) {
          targetBodyTag = bd_tg;
        }
      }
    } // Group's bc-id matched

  } // Loop boundary groups

  // Store data
  // ==========
  nofTargetBoundaries = be_tags.size();
  targetBoundaryTags = new int[nofTargetBoundaries];

  IdsSet::iterator itr = be_tags.begin();

  for (int i = 0; i < nofTargetBoundaries; i++, itr++) {
    targetBoundaryTags[i] = *itr;
  }

}

