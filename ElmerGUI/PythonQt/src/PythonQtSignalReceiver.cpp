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
// \file    PythonQtSignalReceiver.cpp
// \author  Florian Link
// \author  Last changed by $Author: florian $
// \date    2006-05
*/
//----------------------------------------------------------------------------------

#include "PythonQtSignalReceiver.h"
#include "PythonQtClassInfo.h"
#include "PythonQtMethodInfo.h"
#include "PythonQtConversion.h"
#include <QMetaObject>
#include <QMetaMethod>

void PythonQtSignalTarget::call(void **arguments) const
{
  const PythonQtMethodInfo* m = methodInfo();
  // paramterCount includes return value:
  int count = m->parameterCount();

  PyObject* pargs = NULL;
  if (count>1) {
    pargs = PyTuple_New(count-1);
  }
  bool err = false;
  // transform Qt values to Python
  const QList<PythonQtMethodInfo::ParameterInfo>& params = m->parameters();
  for (int i = 1; i < count; i++) {
    const PythonQtMethodInfo::ParameterInfo& param = params.at(i);
    PyObject* arg = PythonQtConv::ConvertQtValueToPython(param, arguments[i]);
    if (arg) {
      // steals reference, no unref
      PyTuple_SetItem(pargs, i-1,arg);
    } else {
      err = true;
      break;
    }
  }

  if (!err) {
    PyErr_Clear();
    PyObject* result = PyObject_CallObject(_callable, pargs);
    if (result) {
      // ok
      Py_DECREF(result);
    } else {
      PythonQt::self()->handleError();
    }
  }
  if (pargs) {
    // free the arguments again
    Py_DECREF(pargs);
  }
}

//------------------------------------------------------------------------------

PythonQtSignalReceiver::PythonQtSignalReceiver(QObject* obj):PythonQtSignalReceiverBase(obj)
{
  _obj = obj;
  _slotCount = staticMetaObject.methodOffset();
}

PythonQtSignalReceiver::~PythonQtSignalReceiver()
{
}


bool PythonQtSignalReceiver::addSignalHandler(const char* signal, PyObject* callable)
{
  bool flag = false;
  int sigId = getSignalIndex(signal);
  if (sigId>=0) {
    // create PythonQtMethodInfo from signal
    QMetaMethod meta = _obj->metaObject()->method(sigId);
    const PythonQtMethodInfo* signalInfo = PythonQtMethodInfo::getCachedMethodInfo(meta);
    PythonQtSignalTarget t(sigId, signalInfo, _slotCount, callable);
    _targets.append(t);
    // now connect to ourselves with the new slot id
    QMetaObject::connect(_obj, sigId, this, _slotCount, Qt::AutoConnection, 0);

    _slotCount++;
    flag = true;
  }
  return flag;
}

bool PythonQtSignalReceiver::removeSignalHandler(const char* signal, PyObject* callable)
{
  bool found = false;
  int sigId = getSignalIndex(signal);
  if (sigId>=0) {
    QMutableListIterator<PythonQtSignalTarget> i(_targets);
    while (i.hasNext()) {
      if (i.next().isSame(sigId, callable)) {
        i.remove();
        found = true;
        break;
      }
    }
  }
  return found;
}

void PythonQtSignalReceiver::removeSignalHandlers()
{
  _targets.clear();
}

int PythonQtSignalReceiver::getSignalIndex(const char* signal)
{
  int sigId = _obj->metaObject()->indexOfSignal(signal+1);
  if (sigId<0) {
    QByteArray tmpSig = QMetaObject::normalizedSignature(signal+1);
    sigId = _obj->metaObject()->indexOfSignal(tmpSig);
  }
  return sigId;
}

int PythonQtSignalReceiver::qt_metacall(QMetaObject::Call c, int id, void **arguments)
{
//  mlabDebugConst("PythonQt", "PythonQtSignalReceiver invoke " << _obj->className() << " " << _obj->name() << " " << id);
  if (c != QMetaObject::InvokeMetaMethod) {
    QObject::qt_metacall(c, id, arguments);
  }

  bool found = false;
  foreach(const PythonQtSignalTarget& t, _targets) {
    if (t.slotId() == id) {
      found = true;
      t.call(arguments);
      break;
    }
  }
  return 0;
}

