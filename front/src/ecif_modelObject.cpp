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
Module:     ecif_modelObject.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Implementation

************************************************************************/

#include "ecif_modelObject.h"
#include "ecif_func.h"

// Init static class variables
Model* ModelObject::model = NULL;


// =================
// ModelObject class
// =================

ModelObject::ModelObject()
{
  id = NO_INDEX;
  objectOk = true;
  otype = OT_NONE;

  tag = NO_INDEX;
  name = NULL;
  active = true;
}

ModelObject::ModelObject(int oid, enum objectType tp, int tg, char* nm)
{
  id = oid;
  objectOk = true;
  otype = tp;

  tag = tg;
  name = NULL; update_dyna_string(name, nm);
  active = true;
}


ModelObject::~ModelObject()
{
  delete[] name;
}


// Check if name defined
//
bool
ModelObject::hasName() const
{
  if ( name == NULL || name[0] == '\0' ) {
    return false;
  } else {
    return true;
  }
}


void
ModelObject::initClass(Model* mdl)
{
  ModelObject::model = mdl;
}


