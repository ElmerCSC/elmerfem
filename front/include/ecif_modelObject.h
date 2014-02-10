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
Module:     ecif_modelObject.h
Language:   C++
Date:       20.12.99
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Common base classes for model level objects

************************************************************************/

#ifndef _ECIF_MODELOBJECT_
#define _ECIF_MODELOBJECT_

#include "ecif_def.h"


// Absract base class for all model objects
class ModelObject {
public:
  ModelObject();
  ModelObject(int oid, enum objectType tp, int tg, char* nm);
  virtual ~ModelObject();
  virtual int Id() const { return id; }
  bool isActive() { return active; }
  bool objectIsOk() { return objectOk; }
  virtual enum objectType getObjectType() const { return otype; }
  virtual const char* getName() const { return name; }
  virtual bool hasName() const;
  virtual int Tag() const { return tag; }
  static void initClass(Model* mdl);
  virtual void setId(int oid) { id = oid; }
  virtual void setObjectType(enum objectType tp) { otype = tp; }
  virtual void setTag(int tg) { tag = tg; }
  virtual void setName(char* nm) { update_dyna_string(name, nm); }
  virtual void setActive(bool value) { active = value; }

protected:
  static Model* model;
  int id;
  bool objectOk;
  enum objectType otype;
  int tag;
  char* name;
  bool active;
};

#endif
