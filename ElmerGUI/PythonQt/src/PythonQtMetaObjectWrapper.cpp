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
// \file    PythonQtMetaObjectWrapper.cpp
// \author  Florian Link
// \author  Last changed by $Author: florian $
// \date    2006-05
*/
//----------------------------------------------------------------------------------

#include "PythonQtMetaObjectWrapper.h"
#include <QObject>

#include "PythonQt.h"
#include "PythonQtSlot.h"
#include "PythonQtClassInfo.h"
#include "PythonQtConversion.h"

static void PythonQtMetaObjectWrapper_dealloc(PythonQtMetaObjectWrapper* self)
{
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* PythonQtMetaObjectWrapper_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PythonQtMetaObjectWrapper *self;

  self = (PythonQtMetaObjectWrapper *)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->_info = NULL;
  }
  return (PyObject *)self;
}

static int PythonQtMetaObjectWrapper_init(PythonQtMetaObjectWrapper *self, PyObject *args, PyObject *kwds)
{
  return 0;
}

PyObject *PythonQtMetaObjectWrapper_call(PyObject *func, PyObject *args, PyObject *kw) {
  PythonQtMetaObjectWrapper* wrapper = (PythonQtMetaObjectWrapper*)func;
  PyObject* result = NULL;
  QString error;
  PyObject* err = NULL;
  if (wrapper->_info->constructors()) {
    result = PythonQtSlotFunction_CallImpl(NULL, wrapper->_info->constructors(), args, kw);
    err = PyErr_Occurred();
  }
  if (!result) {
    QObject* v = NULL;
    QListIterator<PythonQtConstructorHandler*> it(PythonQt::self()->constructorHandlers());
    while (!v && it.hasNext()) {
      v = it.next()->create(wrapper->_info->metaObject(), args, kw, error);
    }
    if (v) {
      result = PythonQt::priv()->wrapQObject(v);
    }
  }
  if (result) {
    // change ownershipflag to be owned by PythonQt
    if (result->ob_type == &PythonQtWrapper_Type) {
      ((PythonQtWrapper*)result)->_ownedByPythonQt = true;
    }
  } else {
    if (!wrapper->_info->constructors()) {
      if (!err) {
        if (error.isEmpty()) {
          error = QString("No constructors available for ") + wrapper->_info->className();
        }
        PyErr_SetString(PyExc_ValueError, error.toLatin1().data());
      }
    }
  }
  return result;
}

static PyObject *PythonQtMetaObjectWrapper_classname(PythonQtMetaObjectWrapper* type)
{
  return PyString_FromString((QString("Meta_") + type->_info->className()).toLatin1().data());
}

static PyObject *PythonQtMetaObjectWrapper_help(PythonQtMetaObjectWrapper* type)
{
  return PythonQt::self()->helpCalled(type->_info);
}


static PyMethodDef PythonQtMetaObjectWrapper_methods[] = {
    {"className", (PyCFunction)PythonQtMetaObjectWrapper_classname, METH_NOARGS,
     "Return the classname of the object"
    },
    {"help", (PyCFunction)PythonQtMetaObjectWrapper_help, METH_NOARGS,
    "Shows the help of available methods for this class"
    },
    {NULL}  /* Sentinel */
};


static PyObject *PythonQtMetaObjectWrapper_getattro(PyObject *obj,PyObject *name)
{
  const char *attributeName;
  PythonQtMetaObjectWrapper *wt = (PythonQtMetaObjectWrapper *)obj;

  if ((attributeName = PyString_AsString(name)) == NULL) {
    return NULL;
  }

  PythonQtMemberInfo member = wt->_info->member(attributeName);
  if (member._type == PythonQtMemberInfo::EnumValue) {
    return PyInt_FromLong(member._enumValue);
  }
  if (member._type == PythonQtMemberInfo::Slot && member._slot->isClassDecorator()) {
    return PythonQtSlotFunction_New(member._slot, obj, NULL);
  }
  
  // look for the interal methods (className(), help())
  PyObject* internalMethod = Py_FindMethod( PythonQtMetaObjectWrapper_methods, obj, (char*)attributeName);
  if (internalMethod) {
    return internalMethod;
  }
  PyErr_Clear();

  if (qstrcmp(attributeName, "__dict__")==0) {
    QStringList l = wt->_info->memberList(true);
    PyObject* dict = PyDict_New();
    foreach (QString name, l) {
      //PyObject* o = PyObject_GetAttrString(obj, name.toLatin1().data());
      PyDict_SetItemString(dict, name.toLatin1().data(), Py_None);
      //Py_DECREF(o);
    }
    return dict;
  }

  QString error = QString(wt->_info->className()) + " has no attribute named '" + QString(attributeName) + "'";
  PyErr_SetString(PyExc_AttributeError, error.toLatin1().data());
  return NULL;
}

static PyObject * PythonQtMetaObjectWrapper_repr(PyObject * obj)
{
  PythonQtMetaObjectWrapper* wt = (PythonQtMetaObjectWrapper*)obj;
  if (wt->_info->isCPPWrapper()) {
    return PyString_FromFormat("%s Class (C++ wrapped by %s)", wt->_info->className(), wt->_info->metaObject()->className());
  } else {
    return PyString_FromFormat("%s Class", wt->_info->className());
  }
}

static int PythonQtMetaObjectWrapper_compare(PyObject * obj1, PyObject * obj2)
{
  if (obj1->ob_type == &PythonQtMetaObjectWrapper_Type &&
    obj2->ob_type == &PythonQtMetaObjectWrapper_Type) {

    PythonQtMetaObjectWrapper* w1 = (PythonQtMetaObjectWrapper*)obj1;
    PythonQtMetaObjectWrapper* w2 = (PythonQtMetaObjectWrapper*)obj2;
    if (w1->_info == w2->_info) {
      return 0;
    } else {
      return -1;
    }
  } else {
    return -1;
  }
}

PyTypeObject PythonQtMetaObjectWrapper_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "PythonQt.PythonQtMetaObjectWrapper",             /*tp_name*/
    sizeof(PythonQtMetaObjectWrapper),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)PythonQtMetaObjectWrapper_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    PythonQtMetaObjectWrapper_compare,         /*tp_compare*/
    PythonQtMetaObjectWrapper_repr,            /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    PythonQtMetaObjectWrapper_call,                         /*tp_call*/
    0,                         /*tp_str*/
    PythonQtMetaObjectWrapper_getattro,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "PythonQtMetaObjectWrapper object",           /* tp_doc */
    0,                   /* tp_traverse */
    0,                   /* tp_clear */
    0,                   /* tp_richcompare */
    0,                   /* tp_weaklistoffset */
    0,                   /* tp_iter */
    0,                   /* tp_iternext */
    0,             /* tp_methods */
    0,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)PythonQtMetaObjectWrapper_init,      /* tp_init */
    0,                         /* tp_alloc */
    PythonQtMetaObjectWrapper_new,                 /* tp_new */
};

//-------------------------------------------------------

