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

//----------------------------------------------------------------------------------
/*!
// \file    PythonQt.cpp
// \author  Florian Link
// \author  Last changed by $Author: florian $
// \date    2006-05
*/
//----------------------------------------------------------------------------------

#include "PythonQtClassInfo.h"
#include "PythonQtMethodInfo.h"
#include "PythonQt.h"
#include <QMetaMethod>

QHash<QByteArray, int> PythonQtMethodInfo::_parameterTypeDict;

PythonQtClassInfo::PythonQtClassInfo(const QMetaObject* meta, const QByteArray& wrappedClassName) {
  _meta = meta;
  _wrappedClassName = wrappedClassName;
  _constructors = NULL;
}

PythonQtClassInfo::~PythonQtClassInfo()
{
  QHashIterator<QByteArray, PythonQtMemberInfo> i(_cachedMembers);
  while (i.hasNext()) {
    PythonQtMemberInfo member = i.next().value();
    if (member._type== PythonQtMemberInfo::Slot) {
      PythonQtSlotInfo* info = member._slot;
      while (info) {
        PythonQtSlotInfo* next = info->nextInfo();
        delete info;
        info = next;
      }
    }
  }
}

int PythonQtClassInfo::findCharOffset(const char* sigStart, char someChar)
{
  const char* sigEnd = sigStart;
  char c;
  do {
    c = *sigEnd++;
  } while (c!=someChar && c!=0);
  return sigEnd-sigStart-1;
}
          
PythonQtMemberInfo PythonQtClassInfo::member(const char* memberName)
{
  PythonQtMemberInfo info = _cachedMembers.value(memberName);
  if (info._type != PythonQtMemberInfo::Invalid) {
    return info;
  } else {
    bool found = false;
    const char* attributeName = memberName;
    bool nameMapped = false;
    // look for properties
    int i = _meta->indexOfProperty(attributeName);
    if (i==-1) {
      // try to map name to objectName
      if (qstrcmp(attributeName, "name")==0) {
        attributeName = "objectName";
        nameMapped = true;
        i = _meta->indexOfProperty(attributeName);
      }
    }
    if (i!=-1) {
      PythonQtMemberInfo newInfo(_meta->property(i));
      _cachedMembers.insert(attributeName, newInfo);
      if (nameMapped) {
        _cachedMembers.insert(memberName, newInfo);
      }
#ifdef PYTHONQT_DEBUG
      std::cout << "caching property " << memberName << " on " << _meta->className() << std::endl;
#endif
      found = true;
    } else {
      int memberNameLen = strlen(memberName);
      // if it is not a property, look for slots
      PythonQtSlotInfo* tail = NULL;
      int numMethods = _meta->methodCount();
      for (int i = 0; i < numMethods; i++) {
        QMetaMethod m = _meta->method(i);
        if ((m.methodType() == QMetaMethod::Method ||
          m.methodType() == QMetaMethod::Slot) && m.access() == QMetaMethod::Public) {
          
          const char* sigStart = m.signature();
          // find the first '('
          int offset = findCharOffset(sigStart, '(');
          
          // check if same length and same name
          if (memberNameLen == offset && qstrncmp(memberName, sigStart, offset)==0) {
            found = true;
            PythonQtSlotInfo* info = new PythonQtSlotInfo(m, i);
            if (tail) {
              tail->setNextInfo(info);
            } else {
              PythonQtMemberInfo newInfo(info);
              _cachedMembers.insert(memberName, newInfo);
            }
            tail = info;
          }
        }
      }
      
      // look for decorators
      if (!_wrappedClassName.isEmpty()) {
        tail = findDecoratorSlots(_wrappedClassName.constData(), memberName, memberNameLen, tail, found);
      }
      const QMetaObject* meta = _meta;
      while (meta) {
        tail = findDecoratorSlots(meta->className(), memberName, memberNameLen, tail, found);
        meta = meta->superClass();
      }

    }
    if (!found) {
      // look for enum values
      int enumCount = _meta->enumeratorCount();
      for (i=0;i<enumCount; i++) {
        QMetaEnum e = _meta->enumerator(i);
        for (int j=0; j < e.keyCount(); j++) {
          if (qstrcmp(e.key(j), attributeName)==0) {
            PythonQtMemberInfo newInfo(e.value(j));
            _cachedMembers.insert(memberName, newInfo);
#ifdef PYTHONQT_DEBUG
            std::cout << "caching enum " << memberName << " on " << _meta->className() << std::endl;
#endif
            found = true;
          }
        }
      }
    }
    return _cachedMembers.value(memberName);
  }
}

PythonQtSlotInfo* PythonQtClassInfo::findDecoratorSlots(const char* classname, const char* memberName, int memberNameLen, PythonQtSlotInfo* tail, bool& found)
{
  QListIterator<PythonQtSlotInfo*> it(PythonQt::priv()->getDecoratorSlots(classname));
  while (it.hasNext()) {

    PythonQtSlotInfo* infoOrig = it.next();
          
    const char* sigStart = infoOrig->metaMethod()->signature();
    if (qstrncmp("static_", sigStart, 7)==0) {
      sigStart += 7;
      sigStart += findCharOffset(sigStart, '_')+1;
    }
    int offset = findCharOffset(sigStart, '(');
    if (memberNameLen == offset && qstrncmp(memberName, sigStart, offset)==0) {
      //make a copy, otherwise we will have trouble on overloads!
      PythonQtSlotInfo* info = new PythonQtSlotInfo(*infoOrig);
      found = true;
      if (tail) {
        tail->setNextInfo(info);
      } else {
        PythonQtMemberInfo newInfo(info);
        _cachedMembers.insert(memberName, newInfo);
      }
      tail = info;
    }
  }
  return tail;
}

 
QStringList PythonQtClassInfo::memberList(bool metaOnly)
{
  QStringList l;
  QString h;
  if (_wrappedClassName.isEmpty()) {
    int i;
    int numProperties = _meta->propertyCount();
    for (i = 0; i < numProperties; i++) {
      QMetaProperty p = _meta->property(i);
      l << QString(p.name());
    }
  }
  
  if (!metaOnly) {
    int numMethods = _meta->methodCount();
    bool skipQObj = !_wrappedClassName.isEmpty();
    for (int i = skipQObj?QObject::staticMetaObject.methodCount():0; i < numMethods; i++) {
      QMetaMethod m = _meta->method(i);
      if ((m.methodType() == QMetaMethod::Method ||
        m.methodType() == QMetaMethod::Slot) && m.access() == QMetaMethod::Public) {
        QByteArray signa(m.signature()); 
        if (signa.startsWith("new_")) continue;
        if (signa.startsWith("delete_")) continue;
        if (signa.startsWith("static_")) continue;
        PythonQtSlotInfo slot(m, i);
        l << slot.slotName();
      }
    }
  }
  // look for decorators
  QList<const char*> names;
  if (!_wrappedClassName.isEmpty()) {
    names << _wrappedClassName.constData();
  }
  const QMetaObject* meta = _meta;
  while (meta) {
    if (meta==&QObject::staticMetaObject && !_wrappedClassName.isEmpty()) break;
    names << meta->className();
    meta = meta->superClass();
  }

  QListIterator<const char*> nameIt(names);
  while (nameIt.hasNext()) {
    QListIterator<PythonQtSlotInfo*> it(PythonQt::priv()->getDecoratorSlots(nameIt.next()));
    while (it.hasNext()) {
      PythonQtSlotInfo* slot = it.next();
      if (metaOnly) {
        if (slot->isClassDecorator()) {
          QByteArray first = slot->slotName();
          if (first.startsWith("static_")) {
            int idx = first.indexOf('_');
            idx = first.indexOf('_', idx+1);
            first = first.mid(idx+1);
          }
          l << first;
        }
      } else {
        l << slot->slotName();
      }
    }
  }
   
  if (_meta->enumeratorCount()) {
    for (int i = 0; i<_meta->enumeratorCount(); i++) {
      QMetaEnum e = _meta->enumerator(i);
      for (int j=0; j < e.keyCount(); j++) {
        l << QString(e.key(j));
      }
    }
  }
  return l;
}

const char* PythonQtClassInfo::className()
{
  if (!_wrappedClassName.isEmpty()) {
    return _wrappedClassName.constData();
  } else {
    return _meta->className();
  }
}

bool PythonQtClassInfo::inherits(const char* name)
{
  const QMetaObject* m = _meta;
  while (m) {
    if (strcmp(name, m->className())==0) {
      return true;
    }
    m = m->superClass();
  }
  return false;
}

const QByteArray& PythonQtClassInfo::wrappedCPPClassName()
{
  return _wrappedClassName;
}

QString PythonQtClassInfo::help()
{
  QString h;
  bool isVariant = QMetaType::type(className())!=QMetaType::Void;
  h += QString("--- ") + QString(className()) + QString(" ---\n");
  
  if (_wrappedClassName.isEmpty()) {
    h += "Properties:\n";
  
    int i;
    int numProperties = _meta->propertyCount();
    for (i = 0; i < numProperties; i++) {
      QMetaProperty p = _meta->property(i);
      h += QString(p.name()) + " (" + QString(p.typeName()) + " )\n";
    }
  }
  
  if (constructors()) {
    h += "Constructors:\n";
    PythonQtSlotInfo* constr = constructors();
    while (constr) {
      h += constr->fullSignature(false) + "\n";
      constr = constr->nextInfo();
    }
  }

  h += "Slots:\n";
  h += "QString help()\n";
  h += "QString className()\n";

  int numMethods = _meta->methodCount();
  for (int i = isVariant?QObject::staticMetaObject.methodCount():0; i < numMethods; i++) {
    QMetaMethod m = _meta->method(i);
    if ((m.methodType() == QMetaMethod::Method ||
      m.methodType() == QMetaMethod::Slot) && m.access() == QMetaMethod::Public) {
      QByteArray signa(m.signature()); 
      if (signa.startsWith("new_")) continue;
      if (signa.startsWith("delete_")) continue;
      if (signa.startsWith("static_")) continue;
      PythonQtSlotInfo slot(m, i);
      h += slot.fullSignature(isVariant)+ "\n";
    }
  }
  // look for decorators
  QList<const char*> names;
  if (!_wrappedClassName.isEmpty()) {
    names << _wrappedClassName.constData();
  }
  const QMetaObject* meta = _meta;
  while (meta) {
    names << meta->className();
    meta = meta->superClass();
    if (isVariant && meta==&QObject::staticMetaObject) break;
  }

  QListIterator<const char*> nameIt(names);
  while (nameIt.hasNext()) {
    QListIterator<PythonQtSlotInfo*> it(PythonQt::priv()->getDecoratorSlots(nameIt.next()));
    while (it.hasNext()) {
      PythonQtSlotInfo* slot = it.next();
      h += slot->fullSignature(slot->isInstanceDecorator()) + "\n";
    }
  }
   
  if (_meta->enumeratorCount()) {
    h += "Enums:\n";
    for (int i = 0; i<_meta->enumeratorCount(); i++) {
      QMetaEnum e = _meta->enumerator(i);
      h += QString(e.name()) + " {";
      for (int j=0; j < e.keyCount(); j++) {
        if (j) { h+= ", "; }
        h += e.key(j);
      }
      h += " }\n";
    }
  }

  if (_wrappedClassName.isEmpty()) {
    int numMethods = _meta->methodCount();
    if (numMethods>0) {
      h += "Signals:\n";
      for (int i = isVariant?QObject::staticMetaObject.methodCount():0; i < numMethods; i++) {
        QMetaMethod m = _meta->method(i);
        if (m.methodType() == QMetaMethod::Signal) {
          h += QString(m.signature()) + "\n";
        }
      }
    }
  }
  return h;
}

PythonQtSlotInfo* PythonQtClassInfo::constructors()
{
  if (!_constructors) {
    _constructors = PythonQt::priv()->getConstructorSlot(!_wrappedClassName.isEmpty()?_wrappedClassName:QByteArray(_meta->className()));
  }
  return _constructors;
}

void PythonQtClassInfo::setMetaObject(const QMetaObject* meta)
{
  _meta = meta;
  _cachedMembers.clear();
}
