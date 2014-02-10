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
// \file    PythonQtVariantWrapper.cpp
// \author  Florian Link
// \author  Last changed by $Author: florian $
// \date    2006-05
*/
//----------------------------------------------------------------------------------

#include "PythonQtVariantWrapper.h"
#include <QObject>
#include <QDate>
#include <QDateTime>
#include <QTime>
#include "PythonQt.h"
#include "PythonQtSlot.h"
#include "PythonQtClassInfo.h"
#include "PythonQtConversion.h"

static void PythonQtVariantWrapper_dealloc(PythonQtVariantWrapper* self)
{
  if (self->_variant) {
    delete self->_variant;
    self->_variant = NULL;
  }
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* PythonQtVariantWrapper_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PythonQtVariantWrapper *self;
  
  self = (PythonQtVariantWrapper *)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->_variant = new QVariant();
    self->_info = NULL;
  }
  return (PyObject *)self;
}

static int PythonQtVariantWrapper_init(PythonQtVariantWrapper *self, PyObject *args, PyObject *kwds)
{
  return 0;
}

static PyObject *PythonQtVariantWrapper_classname(PythonQtVariantWrapper* type)
{
  return PyString_FromString(type->_info->className());
}

static PyObject *PythonQtVariantWrapper_help(PythonQtVariantWrapper* type)
{
  return PythonQt::self()->helpCalled(type->_info);
}


static PyMethodDef PythonQtVariantWrapper_methods[] = {
    {"className", (PyCFunction)PythonQtVariantWrapper_classname, METH_NOARGS,
     "Return the classname of the object"
    },
    {"help", (PyCFunction)PythonQtVariantWrapper_help, METH_NOARGS,
    "Shows the help of available methods for this class"
    },
    {NULL}  /* Sentinel */
};


static PyObject *PythonQtVariantWrapper_getattro(PyObject *obj,PyObject *name)
{
  const char *attributeName;
  PythonQtVariantWrapper *wt = (PythonQtVariantWrapper *)obj;
  
  if ((attributeName = PyString_AsString(name)) == NULL) {
    return NULL;
  }

  if (wt->_wrapper && wt->_info) {
    PythonQtMemberInfo member = wt->_info->member(attributeName);
    if (member._type == PythonQtMemberInfo::Slot) {
      return PythonQtSlotFunction_New(member._slot, obj, NULL);
    } else if (member._type == PythonQtMemberInfo::EnumValue) {
      return PyInt_FromLong(member._enumValue);
    }
  }

  // look for the interal methods (className(), help())
  PyObject* internalMethod = Py_FindMethod( PythonQtVariantWrapper_methods, obj, (char*)attributeName);
  if (internalMethod) {
    return internalMethod;
  }
  PyErr_Clear();

  if (qstrcmp(attributeName, "__dict__")==0) {
    QStringList l = wt->_info->memberList(false);
    PyObject* dict = PyDict_New();
    foreach (QString name, l) {
      //PyObject* o = PyObject_GetAttrString(obj, name.toLatin1().data());
      PyDict_SetItemString(dict, name.toLatin1().data(), Py_None);
      //Py_DECREF(o);
    }
    return dict;
  }

  QString error = QString(wt->_variant->typeName()) + " has no attribute named '" + QString(attributeName) + "'";
  PyErr_SetString(PyExc_AttributeError, error.toLatin1().data());
    
  return NULL;
}

QString qVariantToString(const QVariant& v) {
  QString r;
  switch (v.type()) {
  case QVariant::Size:
    r = QString::number(v.toSize().width()) + ", " + QString::number(v.toSize().height());
    break;
  case QVariant::SizeF:
    r = QString::number(v.toSizeF().width()) + ", " + QString::number(v.toSizeF().height());
    break;
  case QVariant::Point:
    r = QString::number(v.toPoint().x()) + ", " + QString::number(v.toPoint().y());
    break;
  case QVariant::PointF:
    r = QString::number(v.toPointF().x()) + ", " + QString::number(v.toPointF().y());
    break;
  case QVariant::Rect:
    r = QString::number(v.toRect().x()) + ", " + QString::number(v.toRect().y());
    r += ", " + QString::number(v.toRect().width()) + ", " + QString::number(v.toRect().height());
    break;
  case QVariant::RectF:
    r = QString::number(v.toRectF().x()) + ", " + QString::number(v.toRectF().y());
    r += ", " + QString::number(v.toRectF().width()) + ", " + QString::number(v.toRectF().height());
    break;
  case QVariant::Date:
    r = v.toDate().toString("ddMMyyyy");
    break;
  case QVariant::DateTime:
    r = v.toDateTime().toString("ddMMyyyy,hh:mm:ss");
    break;
  case QVariant::Time:
    r = v.toTime().toString("hh:mm:ss");
    break;
    //TODO: add more printing for other variant types
  default:
    r = v.toString();
  }
  return r;
}

static PyObject * PythonQtVariantWrapper_str(PyObject * obj)
{
  PythonQtVariantWrapper* wt = (PythonQtVariantWrapper*)obj;
  QString val = qVariantToString(*wt->_variant);
  return PyString_FromFormat("(%s)", val.toLatin1().constData());
}

static PyObject * PythonQtVariantWrapper_repr(PyObject * obj)
{
  PythonQtVariantWrapper* wt = (PythonQtVariantWrapper*)obj;
  QString val = qVariantToString(*wt->_variant);
  return PyString_FromFormat("%s(%s)", wt->_variant->typeName(), val.toLatin1().constData());
}

static int PythonQtVariantWrapper_compare(PyObject * obj1, PyObject * obj2)
{
  if (obj1->ob_type == &PythonQtVariantWrapper_Type &&
    obj2->ob_type == &PythonQtVariantWrapper_Type) {
    
    PythonQtVariantWrapper* w1 = (PythonQtVariantWrapper*)obj1;
    PythonQtVariantWrapper* w2 = (PythonQtVariantWrapper*)obj2;
    if (*w1->_variant == *w2->_variant) {
      return 0;
    } else {
      return -1;
    }
  } else {
    return -1;
  }
}


PyTypeObject PythonQtVariantWrapper_Type = {
  PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "PythonQt.PythonQtVariantWrapper",             /*tp_name*/
    sizeof(PythonQtVariantWrapper),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)PythonQtVariantWrapper_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    PythonQtVariantWrapper_compare,         /*tp_compare*/
    PythonQtVariantWrapper_repr,            /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    PythonQtVariantWrapper_str,                         /*tp_str*/
    PythonQtVariantWrapper_getattro,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "PythonQtVariantWrapper object",           /* tp_doc */
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
    (initproc)PythonQtVariantWrapper_init,      /* tp_init */
    0,                         /* tp_alloc */
    PythonQtVariantWrapper_new,                 /* tp_new */
};

//-------------------------------------------------------

