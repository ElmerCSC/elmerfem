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
// \file    PythonQtSlot.cpp
// \author  Florian Link
// \author  Last changed by $Author: florian $
// \date    2006-05
*/
//----------------------------------------------------------------------------------

#include "PythonQt.h"
#include "PythonQtSlot.h"
#include "PythonQtWrapper.h"
#include "PythonQtClassInfo.h"
#include "PythonQtMisc.h"
#include "PythonQtConversion.h"
#include <iostream>

#define PYTHONQT_MAX_ARGS 32

// Correction by A. Powell:
#if PY_VERSION_HEX >= 0x02060000 || defined EG_MODS
#define WRITE_RESTRICTED PY_WRITE_RESTRICTED
#endif

PyObject* PythonQtCallSlot(QObject* objectToCall, PyObject* args, bool strict, PythonQtSlotInfo* info, bool isVariantCall, void* firstArgument)
{
  if (isVariantCall && info->isInstanceDecorator()) return NULL;

  static unsigned int recursiveEntry = 0;
  
  // store the current storage position, so that we can get back to this state after a slot is called
  // (do this locally, so that we have all positions on the stack
  PythonQtValueStoragePosition globalValueStoragePos;
  PythonQtValueStoragePosition globalPtrStoragePos;
  PythonQtValueStoragePosition globalVariantStoragePos;
  PythonQtConv::global_valueStorage.getPos(globalValueStoragePos);
  PythonQtConv::global_ptrStorage.getPos(globalPtrStoragePos);
  PythonQtConv::global_variantStorage.getPos(globalVariantStoragePos);
  
  recursiveEntry++;
  
  // the arguments that are passed to qt_metacall
  void* argList[PYTHONQT_MAX_ARGS];
  PyObject* result = NULL;
  int argc = info->parameterCount();
  const QList<PythonQtSlotInfo::ParameterInfo>& params = info->parameters();
  
  bool returnValueIsEnum = false;
  const PythonQtSlotInfo::ParameterInfo& returnValueParam = params.at(0);
  if (returnValueParam.typeId != QMetaType::Void) {
    // extra handling of enum return value
    if (!returnValueParam.isPointer && returnValueParam.typeId == PythonQtMethodInfo::Unknown) {
      returnValueIsEnum = PythonQt::priv()->isEnumType(objectToCall->metaObject(), returnValueParam.name);
      if (returnValueIsEnum) {
        PythonQtValueStorage_ADD_VALUE(PythonQtConv::global_valueStorage, long, 0, argList[0]);
      }
    } else {
      // create empty default value for the return value
      argList[0] = PythonQtConv::CreateQtReturnValue(returnValueParam);
    }
  } else {
    argList[0] = NULL;
  }
  
  const QMetaObject* meta = objectToCall?objectToCall->metaObject():NULL;
  bool ok = true;
  if (info->isInstanceDecorator() || isVariantCall) {
    if (!firstArgument) {
      argList[1] = &objectToCall;
    } else {
      // for the variant call we take the ptr to the variant data, for decorators on CPP objects, we take the cpp ptr
      argList[1] = &firstArgument;
    }
    if (ok) {
      for (int i = 2; i<argc && ok; i++) {
        const PythonQtSlotInfo::ParameterInfo& param = params.at(i);
        //std::cout << param.name.data() << " " << param.typeId << (param.isPointer?"*":"") << (param.isConst?" const":"") << std::endl;
        argList[i] = PythonQtConv::ConvertPythonToQt(param, PyTuple_GET_ITEM(args, i-2), strict, meta);
        if (argList[i]==NULL) {
          ok = false;
          break;
        }
      }
    }
  } else {
    for (int i = 1; i<argc && ok; i++) {
      const PythonQtSlotInfo::ParameterInfo& param = params.at(i);
      //std::cout << param.name.data() << " " << param.typeId << (param.isPointer?"*":"") << (param.isConst?" const":"") << std::endl;
      argList[i] = PythonQtConv::ConvertPythonToQt(param, PyTuple_GET_ITEM(args, i-1), strict, meta);
      if (argList[i]==NULL) {
        ok = false;
        break;
      }
    }
  }
  
  if (ok) {
    (info->decorator()?info->decorator():objectToCall)->qt_metacall(QMetaObject::InvokeMetaMethod, info->slotIndex(), argList);
    
    if (!returnValueIsEnum) {
      result = PythonQtConv::ConvertQtValueToPython(returnValueParam, argList[0]);
    } else {
      result = PyInt_FromLong(*((unsigned int*)argList[0]));
    }
  }
  recursiveEntry--;
  
  // reset the parameter storage position to the stored pos to "pop" the parameter stack
  PythonQtConv::global_valueStorage.setPos(globalValueStoragePos);
  PythonQtConv::global_ptrStorage.setPos(globalPtrStoragePos);
  PythonQtConv::global_variantStorage.setPos(globalVariantStoragePos);
  
  // NOTE: it is important to only return here, otherwise the stack will not be popped!!!
  return result;
}

//-----------------------------------------------------------------------------------

static PythonQtSlotFunctionObject *pythonqtslot_free_list = NULL;

PyObject *PythonQtSlotFunction_Call(PyObject *func, PyObject *args, PyObject *kw)
{
  PythonQtSlotFunctionObject* f = (PythonQtSlotFunctionObject*)func;
  PythonQtSlotInfo*    info = f->m_ml;
  if (f->m_self->ob_type == &PythonQtWrapper_Type) {
    PythonQtWrapper* self = (PythonQtWrapper*) f->m_self;
    return PythonQtSlotFunction_CallImpl(self->_obj, info, args, kw, false, self->_wrappedPtr);
  } else if (f->m_self->ob_type == &PythonQtVariantWrapper_Type) {
    PythonQtVariantWrapper* self = (PythonQtVariantWrapper*) f->m_self;
    if (!info->isClassDecorator()) {
      return PythonQtSlotFunction_CallImpl(self->_wrapper, info, args, kw, true, (void*)self->_variant->constData());
    } else {
      return PythonQtSlotFunction_CallImpl(NULL, info, args, kw);
    }
  } else if (f->m_self->ob_type == &PythonQtMetaObjectWrapper_Type) {
    return PythonQtSlotFunction_CallImpl(NULL, info, args, kw);
  } else {
    return NULL;
  }
}

PyObject *PythonQtSlotFunction_CallImpl(QObject* objectToCall, PythonQtSlotInfo* info, PyObject *args, PyObject *kw, bool isVariantCall, void* firstArg)
{
  int argc = PyTuple_Size(args);
  
#ifdef PYTHONQT_DEBUG
  std::cout << "called " << info->metaMethod()->typeName() << " " << info->metaMethod()->signature() << std::endl;
#endif

  PyObject* r = NULL;
  
  if (info->nextInfo()) {
    // overloaded slot call, try on all slots with strict conversion first
    PythonQtSlotInfo* i = info;
    while (i && r==NULL) {
      bool skipFirst = (i->isInstanceDecorator() || isVariantCall);
      if (i->parameterCount()-1-(skipFirst?1:0) == argc) {
        r = PythonQtCallSlot(objectToCall, args, true, i, isVariantCall, firstArg);
      }
      i = i->nextInfo();
    }
    if (!r) {
      // try on all slots with non-strict conversion
      i = info;
      while (i && r==NULL) {
        bool skipFirst = (i->isInstanceDecorator() || isVariantCall);
        if (i->parameterCount()-1-(skipFirst?1:0) == argc) {
          r = PythonQtCallSlot(objectToCall, args, false, i, isVariantCall, firstArg);
        }
        i = i->nextInfo();
      }
    }
    if (r==0) {
      QString e = QString("Could not find matching overload for given arguments:\n" + PythonQtConv::PyObjGetString(args) + "\n The following slots are available:\n");
      PythonQtSlotInfo* i = info;
      while (i) {
        bool skipFirst = (i->isInstanceDecorator() || isVariantCall);
        e += QString(i->fullSignature(skipFirst)) + "\n";
        i = i->nextInfo();
      }
      PyErr_SetString(PyExc_ValueError, e.toLatin1().data());
    }
  } else {
    // simple (non-overloaded) slot call
    bool skipFirst = (info->isInstanceDecorator() || isVariantCall);
    if (info->parameterCount()-1-(skipFirst?1:0) == argc) {
      r = PythonQtCallSlot(objectToCall, args, false, info, isVariantCall, firstArg);
      if (!r) {
        QString e = QString("Called ") + info->fullSignature(skipFirst) + " with wrong arguments: " + PythonQtConv::PyObjGetString(args);
        PyErr_SetString(PyExc_ValueError, e.toLatin1().data());
      }
    } else {
      QString e = QString("Called ") + info->fullSignature(skipFirst) + " with wrong number of arguments: " + PythonQtConv::PyObjGetString(args);
      PyErr_SetString(PyExc_ValueError, e.toLatin1().data());
    }
  }
  
  return r;
}

PyObject *
PythonQtSlotFunction_New(PythonQtSlotInfo *ml, PyObject *self, PyObject *module)
{
  PythonQtSlotFunctionObject *op;
  op = pythonqtslot_free_list;
  if (op != NULL) {
    pythonqtslot_free_list = (PythonQtSlotFunctionObject *)(op->m_self);
    PyObject_INIT(op, &PythonQtSlotFunction_Type);
  }
  else {
    op = PyObject_GC_New(PythonQtSlotFunctionObject, &PythonQtSlotFunction_Type);
    if (op == NULL)
      return NULL;
  }
  op->m_ml = ml;
  Py_XINCREF(self);
  op->m_self = self;
  Py_XINCREF(module);
  op->m_module = module;
  PyObject_GC_Track(op);
  return (PyObject *)op;
}

PythonQtSlotInfo*
PythonQtSlotFunction_GetSlotInfo(PyObject *op)
{
  if (!PythonQtSlotFunction_Check(op)) {
    PyErr_BadInternalCall();
    return NULL;
  }
  return ((PythonQtSlotFunctionObject *)op) -> m_ml;
}

PyObject *
PythonQtSlotFunction_GetSelf(PyObject *op)
{
  if (!PythonQtSlotFunction_Check(op)) {
    PyErr_BadInternalCall();
    return NULL;
  }
  return ((PythonQtSlotFunctionObject *)op) -> m_self;
}

/* Methods (the standard built-in methods, that is) */

static void
meth_dealloc(PythonQtSlotFunctionObject *m)
{
  PyObject_GC_UnTrack(m);
  Py_XDECREF(m->m_self);
  Py_XDECREF(m->m_module);
  m->m_self = (PyObject *)pythonqtslot_free_list;
  pythonqtslot_free_list = m;
}

static PyObject *
meth_get__doc__(PythonQtSlotFunctionObject *m, void *closure)
{
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
meth_get__name__(PythonQtSlotFunctionObject *m, void *closure)
{
  return PyString_FromString(m->m_ml->metaMethod()->signature());
}

static int
meth_traverse(PythonQtSlotFunctionObject *m, visitproc visit, void *arg)
{
  int err;
  if (m->m_self != NULL) {
    err = visit(m->m_self, arg);
    if (err)
      return err;
  }
  if (m->m_module != NULL) {
    err = visit(m->m_module, arg);
    if (err)
      return err;
  }
  return 0;
}

static PyObject *
meth_get__self__(PythonQtSlotFunctionObject *m, void *closure)
{
  PyObject *self;
  if (PyEval_GetRestricted()) {
    PyErr_SetString(PyExc_RuntimeError,
      "method.__self__ not accessible in restricted mode");
    return NULL;
  }
  self = m->m_self;
  if (self == NULL)
    self = Py_None;
  Py_INCREF(self);
  return self;
}

static PyGetSetDef meth_getsets [] = {
  {"__doc__",  (getter)meth_get__doc__,  NULL, NULL},
  {"__name__", (getter)meth_get__name__, NULL, NULL},
  {"__self__", (getter)meth_get__self__, NULL, NULL},
  {0}
};

#define OFF(x) offsetof(PythonQtSlotFunctionObject, x)

static PyMemberDef meth_members[] = {
  {"__module__",    T_OBJECT,     OFF(m_module), WRITE_RESTRICTED},
  {NULL}
};

static PyObject *
meth_repr(PythonQtSlotFunctionObject *m)
{
  return PyString_FromFormat("<built-in qt slot %s of %s object at %p>",
    m->m_ml->metaMethod()->signature(),
    m->m_self->ob_type->tp_name,
    m->m_self);
}

static int
meth_compare(PythonQtSlotFunctionObject *a, PythonQtSlotFunctionObject *b)
{
  if (a->m_self != b->m_self)
    return (a->m_self < b->m_self) ? -1 : 1;
  if (a->m_ml == b->m_ml)
    return 0;
  if (strcmp(a->m_ml->metaMethod()->signature(), b->m_ml->metaMethod()->signature()) < 0)
    return -1;
  else
    return 1;
}

static long
meth_hash(PythonQtSlotFunctionObject *a)
{
  long x,y;
  if (a->m_self == NULL)
    x = 0;
  else {
    x = PyObject_Hash(a->m_self);
    if (x == -1)
      return -1;
  }
  y = _Py_HashPointer((void*)(a->m_ml));
  if (y == -1)
    return -1;
  x ^= y;
  if (x == -1)
    x = -2;
  return x;
}


PyTypeObject PythonQtSlotFunction_Type = {
  PyObject_HEAD_INIT(&PyType_Type)
    0,
    "builtin_qt_slot",
    sizeof(PythonQtSlotFunctionObject),
    0,
    (destructor)meth_dealloc,     /* tp_dealloc */
    0,          /* tp_print */
    0,          /* tp_getattr */
    0,          /* tp_setattr */
    (cmpfunc)meth_compare,      /* tp_compare */
    (reprfunc)meth_repr,      /* tp_repr */
    0,          /* tp_as_number */
    0,          /* tp_as_sequence */
    0,          /* tp_as_mapping */
    (hashfunc)meth_hash,      /* tp_hash */
    PythonQtSlotFunction_Call,      /* tp_call */
    0,          /* tp_str */
    PyObject_GenericGetAttr,    /* tp_getattro */
    0,          /* tp_setattro */
    0,          /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,/* tp_flags */
    0,          /* tp_doc */
    (traverseproc)meth_traverse,    /* tp_traverse */
    0,          /* tp_clear */
    0,          /* tp_richcompare */
    0,          /* tp_weaklistoffset */
    0,          /* tp_iter */
    0,          /* tp_iternext */
    0,          /* tp_methods */
    meth_members,       /* tp_members */
    meth_getsets,       /* tp_getset */
    0,          /* tp_base */
    0,          /* tp_dict */
};

/* Clear out the free list */

void
PythonQtSlotFunction_Fini(void)
{
  while (pythonqtslot_free_list) {
    PythonQtSlotFunctionObject *v = pythonqtslot_free_list;
    pythonqtslot_free_list = (PythonQtSlotFunctionObject *)(v->m_self);
    PyObject_GC_Del(v);
  }
}

