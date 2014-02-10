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
Module:     ecif_gridparameter.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Implementation

************************************************************************/
 
#include "ecif_gridParameter.h"
#include "ecif_body.h"
#include "ecif_bodyLayer.h"
#include "ecif_model.h"

//Initialize static class variables.
int GridParameter::last_id = 0;
Model* GridParameter::model = NULL;


// Constructors
GridParameter::GridParameter()
{
}


GridParameter::GridParameter(int pid) : Parameter(pid)
{
  if ( pid > last_id ) last_id = pid;
}


GridParameter::GridParameter(int pid, int parent_id, char* data_string, char* param_name)
{
  if ( pid > last_id ) last_id = pid;
  setData(pid, parent_id, data_string, param_name);
}


// Parent object tag for emf-file
int
GridParameter::getParentEmfTag()
{
  // Parent object is a body layer, but we need the body as the parent and
  // the layer as the subParent in emf-file
  //
  // NOTE: This is because we do not store Layers as objects in Emf-file, but
  // they are (if shown explicitely) always under bodies!
  //
  BodyLayer* bl = model->getBodyLayerById(parentId);

  if ( bl == NULL ) return NO_INDEX;

  Body* body = model->getBodyById(bl->getBodyId());

  if ( body == NULL ) return NO_INDEX;

  return body->Tag();
}


objectType
GridParameter::getParentEmfType()
{
  // Parent object is a body layer, but we need the body as the parent and
  // the layer as the subParent in emf-file
  //
  BodyLayer* bl = model->getBodyLayerById(parentId);

  if ( bl == NULL ) return OT_NONE;

  Body* body = model->getBodyById(bl->getBodyId());

  if ( body == NULL ) return OT_NONE;

  return body->getObjectType();
}


int
GridParameter::getSubParentEmfTag()
{
  // Parent object is a body layer, but we need the body as the parent and
  // the layer as the subParent in emf-file
  //
  BodyLayer* bl = model->getBodyLayerById(parentId);

  if ( bl == NULL ) return NO_INDEX;

  return bl->Tag();
}



void
GridParameter::initClass(Model* mdl)
{
  GridParameter::model = mdl;
  GridParameter::last_id = 0;
}


void
GridParameter::setName(char* param_name)
{
  Parameter::setName(param_name, "MeshStructure");
}


void
GridParameter::updateParentId()
{
  // We first read from emf-file the body and the body layer info
  // Then get body-layer object id with this info
  
  parentId = NO_INDEX;

  Body* body = model->getBodyByTag(parentEmfTag);

  if ( body == NULL ) return;

  int layer = body->getLayerIndexByTag(subParentEmfTag);

  parentId = body->getLayerId(layer);
}


void
GridParameter::updateParentInfo(int parent_id)
{
  // Update parent object's id
  //
  parentId = parent_id;

  BodyLayer* lr = model->getBodyLayerById(parent_id);

  if ( lr != NULL ) {

    const Body* bd = lr->getBody();
    
    // Emf parent is always the body
    if ( bd != NULL ) {
      parentEmfTag = bd->Tag();
      parentEmfType = bd->getObjectType();

      // For 'technical' layers we don use emf sub-parents
      //
      if ( IMPLICIT_LAYER == lr->getLayerType() ) {
        subParentEmfTag = NO_INDEX;
        subParentEmfType = OT_NONE;

      // For 'real'layer we the layer's tag
      } else {
        subParentEmfTag = lr->Tag();
        subParentEmfType = lr->getObjectType();
      }

    } else {
      parentEmfTag = NO_INDEX;
      parentEmfType = OT_NONE;
      subParentEmfTag = NO_INDEX;
      subParentEmfType = OT_NONE;
    }
  }
}



