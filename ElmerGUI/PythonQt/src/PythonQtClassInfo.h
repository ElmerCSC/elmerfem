#ifndef _PYTHONQTCLASSINFO_H
#define _PYTHONQTCLASSINFO_H

/*
 *
 *  Copyright (C) 2006 MeVis Research GmbH All Rights Reserved.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  Further, this software is distributed without any warranty that it is
 *  free of the rightful claim of any third person regarding infringement
 *  or the like.  Any license provided herein, whether implied or
 *  otherwise, applies only to this software file.  Patent licenses, if
 *  any, provided herein do not apply to combinations of this program with
 *  other software, or any other product whatsoever.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *  Contact information: MeVis Research GmbH, Universitaetsallee 29,
 *  28359 Bremen, Germany or:
 *
 *  http://www.mevis.de
 *
 */

#include <QMetaObject>
#include <QMetaMethod>
#include <QHash>
#include <QByteArray>
#include <QList>

class PythonQtSlotInfo;

struct PythonQtMemberInfo {
  enum Type {
    Invalid, Slot, EnumValue, Property 
  };

  PythonQtMemberInfo() { _type = Invalid; }
  
  PythonQtMemberInfo(PythonQtSlotInfo* info) {
    _type = Slot;
    _slot = info;
    _enumValue = 0;
  }

  PythonQtMemberInfo(unsigned int enumValue) {
    _type = EnumValue;
    _slot = NULL;
    _enumValue = enumValue;
  }

  PythonQtMemberInfo(const QMetaProperty& prop) {
    _type = Property;
    _slot = NULL;
    _enumValue = 0;
    _property = prop;
  }

  PythonQtSlotInfo* _slot;
  unsigned int      _enumValue;
  QMetaProperty     _property;
  Type              _type;
};

//! a class that stores all required information about a Qt object (and an optional associated C++ class name)
/*! for fast lookup of slots when calling the object from Python
*/
class PythonQtClassInfo {

public:
  PythonQtClassInfo(const QMetaObject* meta, const QByteArray& wrappedClassName = QByteArray());

  ~PythonQtClassInfo();

  //! get the Python method definition for a given slot name (without return type and signature)
  PythonQtMemberInfo member(const char* member);

  PythonQtSlotInfo* constructors();
  
  //! get the Qt classname
  const char* className();

  //! get the classname of the wrapped C++ class
  const QByteArray& wrappedCPPClassName();

  //! returns if the object is a CPP wrapper
  bool isCPPWrapper() { return !_wrappedClassName.isEmpty(); }

  //! get the meta object
  const QMetaObject* metaObject() { return _meta; }

  //! set the meta object, this will reset the caching
  void setMetaObject(const QMetaObject* meta);

  //! returns if the meta object inherits the given classname
  bool inherits(const char* name);
  
  //! get help string for the metaobject
  QString help();

  //! get list of all members
  QStringList memberList(bool metaOnly = false);

private:
  PythonQtSlotInfo* findDecoratorSlots(const char* classname, const char* memberName, int memberNameLen, PythonQtSlotInfo* tail, bool &found);
  int findCharOffset(const char* sigStart, char someChar);
  QHash<QByteArray, PythonQtMemberInfo> _cachedMembers;
  PythonQtSlotInfo*                    _constructors;
  const QMetaObject*                   _meta;
  QByteArray                           _wrappedClassName;
};

//---------------------------------------------------------------


#endif
