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
// \file    PythonQtWrapper.cpp
// \author  Florian Link
// \author  Last changed by $Author: florian $
// \date    2006-05
*/
//----------------------------------------------------------------------------------

#include "PythonQtWrapper.h"
#include <QObject>
#include "PythonQt.h"
#include "PythonQtSlot.h"
#include "PythonQtClassInfo.h"
#include "PythonQtConversion.h"

static void PythonQtWrapper_dealloc(PythonQtWrapper* self)
{
  if (self->_wrappedPtr) {
    
    //mlabDebugConst("Python","c++ wrapper removed " << self->_wrappedPtr << " " << self->_obj->className() << " " << self->_info->wrappedClassName().latin1());
    
    PythonQt::priv()->removeWrapperPointer(self->_wrappedPtr);
    // we own our qobject, so we delete it now:
    delete self->_obj;
    self->_obj = NULL;
    if (self->_ownedByPythonQt) {
      PythonQtSlotInfo* slot = PythonQt::priv()->getDestructorSlot(self->_info->wrappedCPPClassName());
      if (slot) {
        void* args[2];
        args[0] = NULL;
        args[1] = &self->_wrappedPtr;
        slot->decorator()->qt_metacall(QMetaObject::InvokeMetaMethod, slot->slotIndex(), args);
        self->_wrappedPtr = NULL;
      } else {
        // TODO: print a warning? we can not destroy that object
      }
    }
  } else if (self->_obj) {
    //mlabDebugConst("Python","qobject wrapper removed " << self->_obj->className() << " " << self->_info->wrappedClassName().latin1());
    PythonQt::priv()->removeWrapperPointer(self->_obj);
    if (self->_ownedByPythonQt) {
      if (!self->_obj->parent()) {
        delete self->_obj;
        self->_obj = NULL;
      }
    }
  }
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* PythonQtWrapper_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PythonQtWrapper *self;
  
  self = (PythonQtWrapper *)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->_info = NULL;
    self->_obj = NULL;
    self->_wrappedPtr = NULL;
    self->_ownedByPythonQt = false;
  }
  return (PyObject *)self;
}

static int PythonQtWrapper_init(PythonQtWrapper *self, PyObject *args, PyObject *kwds)
{
  return 0;
}

static PyObject *PythonQtWrapper_classname(PythonQtWrapper* type)
{
  return PyString_FromString(type->_info->className());
}

static PyObject *PythonQtWrapper_help(PythonQtWrapper* type)
{
  return PythonQt::self()->helpCalled(type->_info);
}


static PyMethodDef PythonQtWrapper_methods[] = {
    {"className", (PyCFunction)PythonQtWrapper_classname, METH_NOARGS,
     "Return the classname of the object"
    },
    {"help", (PyCFunction)PythonQtWrapper_help, METH_NOARGS,
    "Shows the help of available methods for this class"
    },
    {NULL}  /* Sentinel */
};


static PyObject *PythonQtWrapper_getattro(PyObject *obj,PyObject *name)
{
  const char *attributeName;
  PythonQtWrapper *wt = (PythonQtWrapper *)obj;
  
  if ((attributeName = PyString_AsString(name)) == NULL) {
    return NULL;
  }
  
  if (!wt->_obj && !wt->_wrappedPtr) {
    QString error = QString("Trying to read attribute '") + attributeName + "' from a destroyed " + wt->_info->className() + " object";
    PyErr_SetString(PyExc_ValueError, error.toLatin1().data());
    return NULL;
  }
  
  //  mlabDebugConst("Python","get " << attributeName);
  
  // TODO: dynamic properties are missing
  
  PythonQtMemberInfo member = wt->_info->member(attributeName);
  switch (member._type) {
  case PythonQtMemberInfo::Property:
    if (wt->_obj) {
      return PythonQtConv::QVariantToPyObject(member._property.read(wt->_obj));
    }
    break;
  case PythonQtMemberInfo::Slot:
    return PythonQtSlotFunction_New(member._slot, obj, NULL);
    break;
  case PythonQtMemberInfo::EnumValue:
    return PyInt_FromLong(member._enumValue);
    break;
  }
  
  // look for the interal methods (className(), help())
  PyObject* internalMethod = Py_FindMethod( PythonQtWrapper_methods, obj, (char*)attributeName);
  if (internalMethod) {
    return internalMethod;
  }
  PyErr_Clear();

  if (wt->_obj) {
    // look for a child
    QObjectList children = wt->_obj->children();
    for (int i = 0; i < children.count(); i++) {
      QObject *child = children.at(i);
      if (child->objectName() == attributeName) {
        return PythonQt::self()->priv()->wrapQObject(child);
      }
    }
  }

  if (qstrcmp(attributeName, "__dict__")==0) {
    QStringList l = wt->_info->memberList(false);
    PyObject* dict = PyDict_New();
    foreach (QString name, l) {
      //PyObject* o = PyObject_GetAttrString(obj, name.toLatin1().data());
      PyDict_SetItemString(dict, name.toLatin1().data(), Py_None);
      //Py_DECREF(o);
    }
    // Note: we do not put children into the dict, is would look confusing?!
    return dict;
  }

  
  QString error = QString(wt->_info->className()) + " has no attribute named '" + QString(attributeName) + "'";
  PyErr_SetString(PyExc_AttributeError, error.toLatin1().data());
  return NULL;
}

static int PythonQtWrapper_setattro(PyObject *obj,PyObject *name,PyObject *value)
{
  QString error;
  char *attributeName;
  PythonQtWrapper *wt = (PythonQtWrapper *)obj;
  
  if ((attributeName = PyString_AsString(name)) == NULL)
    return -1;
  
  if (!wt->_obj) {
    error = QString("Trying to set attribute '") + attributeName + "' on a destroyed " + wt->_info->className() + " object";
    PyErr_SetString(PyExc_AttributeError, error.toLatin1().data());
    return -1;
  }
  
  PythonQtMemberInfo member = wt->_info->member(attributeName);
  if (member._type == PythonQtMemberInfo::Property) {
    QMetaProperty prop = member._property;
    if (prop.isWritable()) {
      QVariant v;
      if (prop.isEnumType()) {
        // this will give us either a string or an int, everything else will probably be an error
        v = PythonQtConv::PyObjToQVariant(value);
      } else {
        int t = prop.userType();
        v = PythonQtConv::PyObjToQVariant(value, t);
      }
      bool success = false;
      if (v.isValid()) {
        success = prop.write(wt->_obj, v);
      }
      if (success) {
        return 0;
      } else {
        error = QString("Property '") + attributeName + "' of type '" +
          prop.typeName() + "' does not accept an object of type "
          + QString(value->ob_type->tp_name) + " (" + PythonQtConv::PyObjGetRepresentation(value) + ")";
      }
    } else {
      error = QString("Property '") + attributeName + "' of " + wt->_info->className() + " object is not writable";
    }
  } else {
    if (member._type == PythonQtMemberInfo::Slot) {
      error = QString("Slot '") + attributeName + "' can not be overwritten on " + wt->_info->className() + " object";
    } else if (member._type == PythonQtMemberInfo::EnumValue) {
      error = QString("EnumValue '") + attributeName + "' can not be overwritten on " + wt->_info->className() + " object";
    }
  }
  
  PyErr_SetString(PyExc_AttributeError, error.toLatin1().data());
  return -1;
}

static PyObject * PythonQtWrapper_repr(PyObject * obj)
{
  PythonQtWrapper* wt = (PythonQtWrapper*)obj;
  if (wt->_wrappedPtr) {
    if (wt->_obj) {
      return PyString_FromFormat("%s (C++ Object 0x%x wrapped by %s 0x%x))", wt->_info->className(), wt->_wrappedPtr, wt->_obj->metaObject()->className(), wt->_obj);
    } else {
      return PyString_FromFormat("%s (C++ Object 0x%x unwrapped)", wt->_info->className(), wt->_wrappedPtr);
    }
  } else {
    return PyString_FromFormat("%s (QObject 0x%x)", wt->_info->className(), wt->_obj, wt->_wrappedPtr);
  }
}

static int PythonQtWrapper_compare(PyObject * obj1, PyObject * obj2)
{
  if (obj1->ob_type == &PythonQtWrapper_Type &&
    obj2->ob_type == &PythonQtWrapper_Type) {
    
    PythonQtWrapper* w1 = (PythonQtWrapper*)obj1;
    PythonQtWrapper* w2 = (PythonQtWrapper*)obj2;
    if (w1->_wrappedPtr != NULL) {
      if (w1->_wrappedPtr == w1->_wrappedPtr) {
        return 0;
      } else {
        return -1;
      }
    } else if (w1->_obj == w2->_obj) {
      return 0;
    } else {
      return -1;
    }
  } else {
    return -1;
  }
}

static int PythonQtWrapper_nonzero(PyObject *obj)
{
  PythonQtWrapper* wt = (PythonQtWrapper*)obj;
  return (wt->_wrappedPtr == NULL && wt->_obj == NULL)?0:1;
}

// we override nb_nonzero, so that one can do 'if' expressions to test for a NULL ptr
static PyNumberMethods PythonQtWrapper_as_number = {
  0,      /* nb_add */
    0,      /* nb_subtract */
    0,      /* nb_multiply */
    0,      /* nb_divide */
    0,      /* nb_remainder */
    0,      /* nb_divmod */
    0,      /* nb_power */
    0,      /* nb_negative */
    0,      /* nb_positive */
    0,      /* nb_absolute */
    PythonQtWrapper_nonzero,      /* nb_nonzero */
    0,      /* nb_invert */
    0,      /* nb_lshift */
    0,      /* nb_rshift */
    0,    /* nb_and */
    0,    /* nb_xor */
    0,    /* nb_or */
    0,      /* nb_coerce */
    0,      /* nb_int */
    0,      /* nb_long */
    0,      /* nb_float */
    0,      /* nb_oct */
    0,      /* nb_hex */
    0,      /* nb_inplace_add */
    0,      /* nb_inplace_subtract */
    0,      /* nb_inplace_multiply */
    0,      /* nb_inplace_divide */
    0,      /* nb_inplace_remainder */
    0,      /* nb_inplace_power */
    0,      /* nb_inplace_lshift */
    0,      /* nb_inplace_rshift */
    0,      /* nb_inplace_and */
    0,      /* nb_inplace_xor */
    0,      /* nb_inplace_or */
    0,      /* nb_floor_divide */
    0,      /* nb_true_divide */
    0,      /* nb_inplace_floor_divide */
    0,      /* nb_inplace_true_divide */
};

PyTypeObject PythonQtWrapper_Type = {
  PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "PythonQt.PythonQtWrapper",             /*tp_name*/
    sizeof(PythonQtWrapper),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)PythonQtWrapper_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    PythonQtWrapper_compare,         /*tp_compare*/
    PythonQtWrapper_repr,            /*tp_repr*/
    &PythonQtWrapper_as_number,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    PythonQtWrapper_getattro,                         /*tp_getattro*/
    PythonQtWrapper_setattro,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "PythonQtWrapper object",           /* tp_doc */
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
    (initproc)PythonQtWrapper_init,      /* tp_init */
    0,                         /* tp_alloc */
    PythonQtWrapper_new,                 /* tp_new */
};

//-------------------------------------------------------

