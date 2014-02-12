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

#include "PythonQt.h"
#include "PythonQtImporter.h"
#include "PythonQtClassInfo.h"
#include "PythonQtMethodInfo.h"
#include "PythonQtSignalReceiver.h"
#include "PythonQtConversion.h"
#include "PythonQtStdOut.h"
#include "PythonQtCppWrapperFactory.h"
#include "PythonQtVariants.h"
#include "PythonQtStdDecorators.h"
#include <pydebug.h>

PythonQt* PythonQt::_self = NULL;


void PythonQt::init(int flags)
{
  if (!_self) {
    _self = new PythonQt(flags);
  }

  PythonQt::self()->addDecorators(new PythonQtStdDecorators());
  
  PythonQt::priv()->addVariantWrapper("QBitArray", new PythonQtQBitArrayWrapper);
  PythonQt::priv()->addVariantWrapper("QDate", new PythonQtQDateWrapper);
  PythonQt::priv()->addVariantWrapper("QTime", new PythonQtQTimeWrapper);
  PythonQt::priv()->addVariantWrapper("QDateTime", new PythonQtQDateTimeWrapper);
  PythonQt::priv()->addVariantWrapper("QUrl", new PythonQtQUrlWrapper);
  PythonQt::priv()->addVariantWrapper("QLocale", new PythonQtQLocaleWrapper);
  PythonQt::priv()->addVariantWrapper("QRect", new PythonQtQRectWrapper);
  PythonQt::priv()->addVariantWrapper("QRectF", new PythonQtQRectFWrapper);
  PythonQt::priv()->addVariantWrapper("QSize", new PythonQtQSizeWrapper);
  PythonQt::priv()->addVariantWrapper("QSizeF", new PythonQtQSizeFWrapper);
  PythonQt::priv()->addVariantWrapper("QLine", new PythonQtQLineWrapper);
  PythonQt::priv()->addVariantWrapper("QLineF", new PythonQtQLineFWrapper);
  PythonQt::priv()->addVariantWrapper("QPoint", new PythonQtQPointWrapper);
  PythonQt::priv()->addVariantWrapper("QPointF", new PythonQtQPointFWrapper);
  PythonQt::priv()->addVariantWrapper("QRegExp", new PythonQtQRegExpWrapper);
  PythonQt::priv()->addVariantWrapper("QFont", new PythonQtQFontWrapper);
  PythonQt::priv()->addVariantWrapper("QPixmap", new PythonQtQPixmapWrapper);
  PythonQt::priv()->addVariantWrapper("QBrush", new PythonQtQBrushWrapper);
  PythonQt::priv()->addVariantWrapper("QColor", new PythonQtQColorWrapper);
  PythonQt::priv()->addVariantWrapper("QPalette", new PythonQtQPaletteWrapper);
  PythonQt::priv()->addVariantWrapper("QIcon", new PythonQtQIconWrapper);
  PythonQt::priv()->addVariantWrapper("QImage", new PythonQtQImageWrapper);
  PythonQt::priv()->addVariantWrapper("QPolygon", new PythonQtQPolygonWrapper);
  PythonQt::priv()->addVariantWrapper("QRegion", new PythonQtQRegionWrapper);
  PythonQt::priv()->addVariantWrapper("QBitmap", new PythonQtQBitmapWrapper);
  PythonQt::priv()->addVariantWrapper("QCursor", new PythonQtQCursorWrapper);
  PythonQt::priv()->addVariantWrapper("QSizePolicy", new PythonQtQSizePolicyWrapper);
  PythonQt::priv()->addVariantWrapper("QKeySequence", new PythonQtQKeySequenceWrapper);
  PythonQt::priv()->addVariantWrapper("QPen", new PythonQtQPenWrapper);
  PythonQt::priv()->addVariantWrapper("QTextLength", new PythonQtQTextLengthWrapper);
  PythonQt::priv()->addVariantWrapper("QTextFormat", new PythonQtQTextFormatWrapper);
  PythonQt::priv()->addVariantWrapper("QMatrix", new PythonQtQMatrixWrapper);
  
}

void PythonQt::cleanup()
{
  if (_self) {
    delete _self;
    _self = NULL;
  }
}

PythonQt::PythonQt(int flags)
{
  _p = new PythonQtPrivate;
  _p->_initFlags = flags;

  _p->_PythonQtObjectPtr_metaId = qRegisterMetaType<PythonQtObjectPtr>("PythonQtObjectPtr");
  
  Py_SetProgramName("PythonQt");
  if (flags & IgnoreSiteModule) {
    // this prevents the automatic importing of Python site files
    Py_NoSiteFlag = 1;
  }
  Py_Initialize();
  
  // add our own python object types for qt object slots
  if (PyType_Ready(&PythonQtSlotFunction_Type) < 0) {
    std::cerr << "could not initialize PythonQtSlotFunction_Type" << ", in " << __FILE__ << ":" << __LINE__ << std::endl;
  }
  Py_INCREF(&PythonQtSlotFunction_Type);
  
  // add our own python object types for qt objects
  if (PyType_Ready(&PythonQtWrapper_Type) < 0) {
    std::cerr << "could not initialize PythonQtWrapper_Type" << ", in " << __FILE__ << ":" << __LINE__ << std::endl;
  }
  Py_INCREF(&PythonQtWrapper_Type);
  
  // add our own python object types for qt objects
  if (PyType_Ready(&PythonQtVariantWrapper_Type) < 0) {
    std::cerr << "could not initialize PythonQtVariantWrapper_Type" << ", in " << __FILE__ << ":" << __LINE__ << std::endl;
  }
  Py_INCREF(&PythonQtVariantWrapper_Type);

  // add our own python object types for qt objects
  if (PyType_Ready(&PythonQtMetaObjectWrapper_Type) < 0) {
    std::cerr << "could not initialize PythonQtMetaObjectWrapper_Type" << ", in " << __FILE__ << ":" << __LINE__ << std::endl;
  }
  Py_INCREF(&PythonQtMetaObjectWrapper_Type);
  
  // add our own python object types for redirection of stdout
  if (PyType_Ready(&PythonQtStdOutRedirectType) < 0) {
    std::cerr << "could not initialize PythonQtStdOutRedirectType" << ", in " << __FILE__ << ":" << __LINE__ << std::endl;
  }
  Py_INCREF(&PythonQtStdOutRedirectType);
  
  initPythonQtModule(flags & RedirectStdOut);

}

PythonQt::~PythonQt() {
  delete _p;
  _p = NULL;
}

PythonQtPrivate::~PythonQtPrivate() {
  {
    QHashIterator<QByteArray, PythonQtClassInfo *> i(_knownQtClasses);
    while (i.hasNext()) {
      delete i.next().value();
    }
  }
  {
    QHashIterator<QByteArray, PythonQtClassInfo *> i(_knownQtWrapperClasses);
    while (i.hasNext()) {
      delete i.next().value();
    }
  }
  {
    QHashIterator<int , QPair<PythonQtClassInfo*, QObject*> > i(_knownVariantWrappers);
    while (i.hasNext()) {
      delete i.next().value().first;
    }
  }
  {
    QHashIterator<QByteArray, PythonQtSlotInfo *> i(_constructorSlots);
    while (i.hasNext()) {
      delete i.next().value();
    }
  }
  {
    QHashIterator<QByteArray, PythonQtSlotInfo *> i(_destructorSlots);
    while (i.hasNext()) {
      delete i.next().value();
    }
  }
  PythonQtConv::global_valueStorage.clear();
  PythonQtConv::global_ptrStorage.clear();
  PythonQtConv::global_variantStorage.clear();

  PythonQtMethodInfo::cleanupCachedMethodInfos();

  delete _qtNamespace;
}

PythonQtImportFileInterface* PythonQt::importInterface()
{
  return _self->_p->_importInterface;
}

void PythonQt::registerClass(const QMetaObject* metaobject)
{
  _p->registerClass(metaobject);
}

void PythonQtPrivate::registerClass(const QMetaObject* metaobject)
{
  // we register all classes in the hierarchy
  const QMetaObject* m = metaobject;
  while (m) {
    PythonQtClassInfo* info = _knownQtClasses.value(m->className());
    if (!info) {
      info = new PythonQtClassInfo(m);
      _knownQtClasses.insert(m->className(), info);
      PyModule_AddObject(_pythonQtModule, m->className(), (PyObject*)createNewPythonQtMetaObjectWrapper(info));
    }
    m = m->superClass();
  }
}

bool PythonQtPrivate::isEnumType(const QMetaObject* meta, const QByteArray& name) {
  int i = meta?meta->indexOfEnumerator(name.constData()):-1;
  if (i!=-1) {
    return true;
  } else {
    // look for scope
    int scopePos = name.indexOf("::");
    if (scopePos != -1) {
      // slit into scope and enum name
      QByteArray enumScope = name.mid(0,scopePos);
      QByteArray enumName = name.mid(scopePos+2);
      if (enumScope == "Qt") {
        // special qt namespace case
        return isEnumType(&staticQtMetaObject, enumName);
      } else {
        // look for known classes as scope
        // TODO: Q_GADGETS are not yet handled
        PythonQtClassInfo* info = _knownQtClasses.value(enumScope);
        if (info) {
          return isEnumType(info->metaObject(), enumName);
        }
      }
    }
  }
  return false;
}

PyObject* PythonQtPrivate::wrapQObject(QObject* obj)
{
  if (!obj) {
    Py_INCREF(Py_None);
    return Py_None;
  }
  PythonQtWrapper* wrap = _wrappedObjects.value(obj);
  if (!wrap) {
    // smuggling it in...
    PythonQtClassInfo* classInfo = _knownQtClasses.value(obj->metaObject()->className());
    if (!classInfo) {
      registerClass(obj->metaObject());
      classInfo = _knownQtClasses.value(obj->metaObject()->className());
    }
    wrap = createNewPythonQtWrapper(obj, classInfo);
    // insert destroyed handler
    connect(obj, SIGNAL(destroyed(QObject*)), this, SLOT(wrappedObjectDestroyed(QObject*)));
    //    mlabDebugConst("MLABPython","new qobject wrapper added " << " " << wrap->_obj->className() << " " << wrap->_info->wrappedClassName().latin1());
  } else {
    Py_INCREF(wrap);
    //    mlabDebugConst("MLABPython","qobject wrapper reused " << wrap->_obj->className() << " " << wrap->_info->wrappedClassName().latin1());
  }
  return (PyObject*)wrap;
}

PyObject* PythonQtPrivate::wrapPtr(void* ptr, const QByteArray& name)
{
  if (!ptr) {
    Py_INCREF(Py_None);
    return Py_None;
  }
  PythonQtWrapper* wrap = _wrappedObjects.value(ptr);
  if (!wrap) {
    PythonQtClassInfo* info = _knownQtClasses.value(name);
    if (!info) {
      // we do not know the metaobject yet, but we might know it by it's name:
      if (_knownQObjectClassNames.find(name)!=_knownQObjectClassNames.end()) {
        // yes, we know it, so we can convert to QObject
        QObject* qptr = (QObject*)ptr;
        registerClass(qptr->metaObject());
        info = _knownQtClasses.value(qptr->metaObject()->className());
      }
    }
    if (info) {
      QObject* qptr = (QObject*)ptr;
      // if the object is a derived object, we want to switch the class info to the one of the derived class:
      if (name!=(qptr->metaObject()->className())) {
        registerClass(qptr->metaObject());
        info = _knownQtClasses.value(qptr->metaObject()->className());
      }
      wrap = createNewPythonQtWrapper(qptr, info);
      // insert destroyed handler
      connect(qptr, SIGNAL(destroyed(QObject*)), this, SLOT(wrappedObjectDestroyed(QObject*)));
      //    mlabDebugConst("MLABPython","new qobject wrapper added " << " " << wrap->_obj->className() << " " << wrap->_info->wrappedClassName().latin1());
    } else {
      // maybe it is a PyObject, which we can return directly
      if (name == "PyObject") {
        PyObject* p = (PyObject*)ptr;
        Py_INCREF(p);
        return p;
      }
        // not a known QObject, so try our wrapper factory:
      QObject* wrapper = NULL;
      for (int i=0; i<_cppWrapperFactories.size(); i++) {
        wrapper = _cppWrapperFactories.at(i)->create(name, ptr);
        if (wrapper) {
          break;
        }
      }
      PythonQtClassInfo* info = _knownQtWrapperClasses.value(name);
      if (!info) {
        info = new PythonQtClassInfo(wrapper?wrapper->metaObject():&QObject::staticQtMetaObject, name);
        _knownQtWrapperClasses.insert(name, info);
        PyModule_AddObject(_pythonQtModule, name, (PyObject*)createNewPythonQtMetaObjectWrapper(info));
      } else {
        if (wrapper && (info->metaObject() != wrapper->metaObject())) {
          info->setMetaObject(wrapper->metaObject());
        }
      }
      wrap = createNewPythonQtWrapper(wrapper, info, ptr);
      //          mlabDebugConst("MLABPython","new c++ wrapper added " << wrap->_wrappedPtr << " " << wrap->_obj->className() << " " << wrap->_info->wrappedClassName().latin1());
    }
  } else {
    Py_INCREF(wrap);
    //mlabDebugConst("MLABPython","c++ wrapper reused " << wrap->_wrappedPtr << " " << wrap->_obj->className() << " " << wrap->_info->wrappedClassName().latin1());
  }
  return (PyObject*)wrap;
}

void PythonQt::registerCPPClassNames(const QStringList& names)
{
  foreach ( QString n, names) {
    QByteArray name = n.toLatin1();
    PythonQtClassInfo* info = _p->_knownQtWrapperClasses.value(name);
    if (!info) {
      info = new PythonQtClassInfo(&QObject::staticMetaObject, name);
      _p->_knownQtWrapperClasses.insert(name, info);
      PyModule_AddObject(_p->_pythonQtModule, name.data(), (PyObject*)_p->createNewPythonQtMetaObjectWrapper(info));
    }
  }
}

PythonQtWrapper* PythonQtPrivate::createNewPythonQtWrapper(QObject* obj, PythonQtClassInfo* info, void* wrappedPtr) {
  PythonQtWrapper* result;
  result = (PythonQtWrapper *)PythonQtWrapper_Type.tp_new(&PythonQtWrapper_Type,
    NULL, NULL);

  result->_obj = obj;
  result->_info = info;
  result->_wrappedPtr = wrappedPtr;
  result->_ownedByPythonQt = false;

  if (wrappedPtr) {
    _wrappedObjects.insert(wrappedPtr, result);
  } else {
    _wrappedObjects.insert(obj, result);
  }
  return result;
}

PythonQtVariantWrapper* PythonQtPrivate::createNewPythonQtVariantWrapper(const QVariant& variant) {
  PythonQtVariantWrapper* result;
  result = (PythonQtVariantWrapper *)PythonQtVariantWrapper_Type.tp_new(&PythonQtVariantWrapper_Type,
    NULL, NULL);

  *result->_variant = variant;
  QPair<PythonQtClassInfo*, QObject*> pair = _knownVariantWrappers.value(variant.userType());
  result->_wrapper = pair.second;
  result->_info = pair.first;
  return result;
}

PythonQtMetaObjectWrapper* PythonQtPrivate::createNewPythonQtMetaObjectWrapper(PythonQtClassInfo* info) {
  PythonQtMetaObjectWrapper* result;
  result = (PythonQtMetaObjectWrapper *)PythonQtMetaObjectWrapper_Type.tp_new(&PythonQtMetaObjectWrapper_Type,
    NULL, NULL);
  result->_info = info;
  return result;
}


PythonQtSignalReceiver* PythonQt::getSignalReceiver(QObject* obj)
{
  PythonQtSignalReceiver* r = _p->_signalReceivers[obj];
  if (!r) {
    r = new PythonQtSignalReceiver(obj);
    _p->_signalReceivers.insert(obj, r);
    // insert destroyed handler
    connect(obj, SIGNAL(destroyed(QObject*)), _p ,SLOT(destroyedSignalEmitter(QObject*)));
  }
  return r;
}

bool PythonQt::addSignalHandler(QObject* obj, const char* signal, PyObject* module, const QString& objectname)
{
  bool flag = false;
  PythonQtObjectPtr callable = lookupCallable(module, objectname);
  if (callable) {
    PythonQtSignalReceiver* r = getSignalReceiver(obj);
    flag = r->addSignalHandler(signal, callable);
    if (!flag) {
      // signal not found
    }
  } else {
    // callable not found
  }
  return flag;
}

bool PythonQt::addSignalHandler(QObject* obj, const char* signal, PyObject* receiver)
{
  bool flag = false;
  PythonQtSignalReceiver* r = getSignalReceiver(obj);
  if (r) {
    flag = r->addSignalHandler(signal, receiver);
  }
  return flag;
}

bool PythonQt::removeSignalHandler(QObject* obj, const char* signal, PyObject* module, const QString& objectname)
{
  bool flag = false;
  PythonQtObjectPtr callable = lookupCallable(module, objectname);
  if (callable) {
    PythonQtSignalReceiver* r = _p->_signalReceivers[obj];
    if (r) {
      flag = r->removeSignalHandler(signal, callable);
    }
  } else {
    // callable not found
  }
  return flag;
}

bool PythonQt::removeSignalHandler(QObject* obj, const char* signal, PyObject* receiver)
{
  bool flag = false;
  PythonQtSignalReceiver* r = _p->_signalReceivers[obj];
  if (r) {
    flag = r->removeSignalHandler(signal, receiver);
  }
  return flag;
}

PythonQtObjectPtr PythonQt::lookupCallable(PyObject* module, const QString& name)
{
  PythonQtObjectPtr p = lookupObject(module, name);
  if (p) {
    if (PyCallable_Check(p)) {
      return p;
    }
  }
  PyErr_Clear();
  return NULL;
}

PythonQtObjectPtr PythonQt::lookupObject(PyObject* module, const QString& name)
{
  QStringList l = name.split('.');
  PythonQtObjectPtr p = module;
  PythonQtObjectPtr prev;
  QString s;
  QByteArray b;
  for (QStringList::ConstIterator i = l.begin(); i!=l.end() && p; ++i) {
    prev = p;
    b = (*i).toLatin1();
    p.setNewRef(PyObject_GetAttrString(p, b.data()));
  }
  PyErr_Clear();
  return p;
}

PythonQtObjectPtr PythonQt::getMainModule() {
  //both borrowed
  PythonQtObjectPtr dict = PyImport_GetModuleDict();
  return PyDict_GetItemString(dict, "__main__");
}

QVariant PythonQt::evalCode(PyObject* module, PyObject* pycode) {
  QVariant result;
  if (pycode) {
    PyObject* r = PyEval_EvalCode((PyCodeObject*)pycode, PyModule_GetDict((PyObject*)module) , PyModule_GetDict((PyObject*)module));
    if (r) {
      result = PythonQtConv::PyObjToQVariant(r);
      Py_DECREF(r);
    } else {
      handleError();
    }
  } else {
    handleError();
  }
  return result;
}

QVariant PythonQt::evalScript(PyObject* module, const QString& script, int start)
{
  QVariant result;
  PythonQtObjectPtr p;
  p.setNewRef(PyRun_String(script.toLatin1().data(), start, PyModule_GetDict(module), PyModule_GetDict(module)));
  if (p) {
    result = PythonQtConv::PyObjToQVariant(p);
  } else {
    handleError();
  }
  return result;
}

void PythonQt::evalFile(PyObject* module, const QString& filename)
{
  PythonQtObjectPtr code = parseFile(filename);
  if (code) {
    evalCode(module, code);
  } else {
    handleError();
  }
}

PythonQtObjectPtr PythonQt::parseFile(const QString& filename)
{
  PythonQtObjectPtr p;
  p.setNewRef(PythonQtImport::getCodeFromPyc(filename));
  if (!p) {
    handleError();
  }
  return p;
}

void PythonQt::addObject(PyObject* module, const QString& name, QObject* object)
{
  PyModule_AddObject(module, name.toLatin1().data(), _p->wrapQObject(object));
}

void PythonQt::addVariable(PyObject* module, const QString& name, const QVariant& v)
{
  PyModule_AddObject(module, name.toLatin1().data(), PythonQtConv::QVariantToPyObject(v));
}

void PythonQt::removeVariable(PyObject* module, const QString& name)
{
  PyObject_DelAttrString(module, name.toLatin1().data());
}

QVariant PythonQt::getVariable(PyObject* module, const QString& objectname)
{
  QVariant result;
  PythonQtObjectPtr obj = lookupObject(module, objectname);
  if (obj) {
    result = PythonQtConv::PyObjToQVariant(obj);
  }
  return result;
}

QStringList PythonQt::introspection(PyObject* module, const QString& objectname, PythonQt::ObjectType type)
{
  QStringList results;
  
  PythonQtObjectPtr object;
  if (objectname.isEmpty()) {
    object = module;
  } else {
    object = lookupObject(module, objectname);
  }

  if (object) {
    if (type == CallOverloads) {
      if (PythonQtSlotFunction_Check(object)) {
        PythonQtSlotFunctionObject* o = (PythonQtSlotFunctionObject*)object.object();
        PythonQtSlotInfo* info = o->m_ml;
        
        while (info) {
          results << info->fullSignature(info->isInstanceDecorator() || o->m_self->ob_type == &PythonQtVariantWrapper_Type);
          info = info->nextInfo();
        }
      } else if (object->ob_type == &PythonQtMetaObjectWrapper_Type) {
        PythonQtMetaObjectWrapper* o = (PythonQtMetaObjectWrapper*)object.object();
        PythonQtSlotInfo* info = o->_info->constructors();
    
        while (info) {
          results << info->fullSignature(false);
          info = info->nextInfo();
        }
      } else {
        PyObject* doc = PyObject_GetAttrString(object, "__doc__");
        if (doc) {
          results << PyString_AsString(doc);
          Py_DECREF(doc);
        }
      }
    } else {
      PyObject* keys = PyObject_Dir(object);
      if (keys) {
        int count = PyList_Size(keys);
        PyObject* key;
        PyObject* value;
        QString keystr;
        for (int i = 0;i<count;i++) {
          key = PyList_GetItem(keys,i);
          value = PyObject_GetAttr(object, key);
          if (!value) continue;
          keystr = PyString_AsString(key);
          static const QString underscoreStr("__");
          if (!keystr.startsWith(underscoreStr)) {
            switch (type) {
            case Anything:
              results << keystr;
              break;
            case Class:
              if (value->ob_type == &PyClass_Type) {
                results << keystr;
              }
              break;
            case Variable:
              if (value->ob_type != &PyClass_Type
                && value->ob_type != &PyCFunction_Type
                && value->ob_type != &PyFunction_Type
                && value->ob_type != &PyModule_Type
                ) {
                results << keystr;
              }
              break;
            case Function:
              if (value->ob_type == &PyFunction_Type ||
                value->ob_type == &PyMethod_Type
                ) {
                results << keystr;
              }
              break;
            case Module:
              if (value->ob_type == &PyModule_Type) {
                results << keystr;
              }
              break;
            default:
              std::cerr << "PythonQt: introspection: unknown case" << ", in " << __FILE__ << ":" << __LINE__ << std::endl;
            }
          }
          Py_DECREF(value);
        }
        Py_DECREF(keys);
      }
    }
  }
  return results;
}

QVariant PythonQt::call(PyObject* module, const QString& name, const QVariantList& args)
{
  QVariant r;
  
  PythonQtObjectPtr callable = lookupCallable(module, name);
  if (callable) {
    PythonQtObjectPtr pargs;
    int count = args.size();
    if (count>0) {
      pargs.setNewRef(PyTuple_New(count));
    }
    bool err = false;
    // transform QVariants to Python
    for (int i = 0; i < count; i++) {
      PyObject* arg = PythonQtConv::QVariantToPyObject(args.at(i));
      if (arg) {
        // steals reference, no unref
        PyTuple_SetItem(pargs, i,arg);
      } else {
        err = true;
        break;
      }
    }
    
    if (!err) {
      PyErr_Clear();
      PythonQtObjectPtr result;
      result.setNewRef(PyObject_CallObject(callable, pargs));
      if (result) {
        // ok
        r = PythonQtConv::PyObjToQVariant(result);
      } else {
        PythonQt::self()->handleError();
      }
    }
  }
  return r;
}

void PythonQt::addInstanceDecorators(QObject* o)
{
  _p->addDecorators(o, true, false);
}

void PythonQt::addClassDecorators(QObject* o)
{
  _p->addDecorators(o, false, true);
}

void PythonQt::addDecorators(QObject* o)
{
  _p->addDecorators(o, true, true);
}

void PythonQt::registerQObjectClassNames(const QStringList& names)
{
  _p->registerQObjectClassNames(names);
}

void PythonQt::setImporter(PythonQtImportFileInterface* importInterface)
{
  static bool first = true;
  if (first) {
    first = false;
    _p->_importInterface = importInterface;
    PythonQtImport::init();
  }
}

void PythonQt::setImporterIgnorePaths(const QStringList& paths)
{
  _p->_importIgnorePaths = paths;
}

const QStringList& PythonQt::getImporterIgnorePaths()
{
  return _p->_importIgnorePaths;
}

void PythonQt::addWrapperFactory(PythonQtCppWrapperFactory* factory)
{
  _p->_cppWrapperFactories.append(factory);
}

void PythonQt::addConstructorHandler(PythonQtConstructorHandler* factory)
{
  _p->_constructorHandlers.append(factory);
}

const QList<PythonQtConstructorHandler*>& PythonQt::constructorHandlers()
{ 
  return _p->_constructorHandlers;
};

//---------------------------------------------------------------------------------------------------
PythonQtPrivate::PythonQtPrivate()
{
  _importInterface = NULL;
}

void PythonQtPrivate::addDecorators(QObject* o, bool instanceDeco, bool classDeco)
{
  o->setParent(this);
  int numMethods = o->metaObject()->methodCount();
  for (int i = 0; i < numMethods; i++) {
    QMetaMethod m = o->metaObject()->method(i);
    if ((m.methodType() == QMetaMethod::Method ||
      m.methodType() == QMetaMethod::Slot) && m.access() == QMetaMethod::Public) {
     const PythonQtMethodInfo* info = PythonQtMethodInfo::getCachedMethodInfo(m);
      if (qstrncmp(m.signature(), "new_", 4)==0) {
        if (!classDeco) continue; 
        // either it returns a * or a QVariant and the name starts with "new_"
        bool isVariantReturn = info->parameters().at(0).typeId == PythonQtMethodInfo::Variant;
        if ((info->parameters().at(0).isPointer || isVariantReturn)) {
          QByteArray signature = m.signature();
          QByteArray nameOfClass = signature.mid(4, signature.indexOf('(')-4);
          PythonQtSlotInfo* prev = _constructorSlots.value(nameOfClass);
          PythonQtSlotInfo* newSlot = new PythonQtSlotInfo(m, i, o, PythonQtSlotInfo::ClassDecorator);
          if (prev) {
            newSlot->setNextInfo(prev->nextInfo());
            prev->setNextInfo(newSlot);
          } else {
            _constructorSlots.insert(nameOfClass, newSlot);
          }
        }
      } else if (qstrncmp(m.signature(), "delete_", 7)==0) {
        if (!classDeco) continue; 
        QByteArray signature = m.signature();
        QByteArray nameOfClass = signature.mid(7, signature.indexOf('(')-7);
        PythonQtSlotInfo* newSlot = new PythonQtSlotInfo(m, i, o, PythonQtSlotInfo::ClassDecorator);
        _destructorSlots.insert(nameOfClass, newSlot);
      } else if (qstrncmp(m.signature(), "static_", 7)==0) {
        if (!classDeco) continue; 
        QByteArray signature = m.signature();
        QByteArray nameOfClass = signature.mid(signature.indexOf('_')+1);
        nameOfClass = nameOfClass.mid(0, nameOfClass.indexOf('_'));
        PythonQtSlotInfo* slotCopy = new PythonQtSlotInfo(m, i, o, PythonQtSlotInfo::ClassDecorator);
        _knownQtDecoratorSlots.insert(nameOfClass, slotCopy);
      } else {
        if (!instanceDeco) continue; 
        if (info->parameters().count()>1) {
          PythonQtMethodInfo::ParameterInfo p = info->parameters().at(1);
          if (p.isPointer) {
            PythonQtSlotInfo* slotCopy = new PythonQtSlotInfo(m, i, o, PythonQtSlotInfo::InstanceDecorator);
            _knownQtDecoratorSlots.insert(p.name, slotCopy);
          }
        }
      }
    }
  }
}

void PythonQtPrivate::registerQObjectClassNames(const QStringList& names)
{
  foreach(QString name, names) {
    _knownQObjectClassNames.insert(name.toLatin1(), true);
  }
}

QList<PythonQtSlotInfo*> PythonQtPrivate::getDecoratorSlots(const QByteArray& className)
{
  return _knownQtDecoratorSlots.values(className);
}

void PythonQtPrivate::wrappedObjectDestroyed(QObject* obj)
{
  // mlabDebugConst("MLABPython","PyWrapper QObject destroyed " << o << " " << o->name() << " " << o->className());
  PythonQtWrapper* wrap = _wrappedObjects[obj];
  if (wrap) {
    _wrappedObjects.remove(obj);
    // remove the pointer but keep the wrapper alive in python
    wrap->_obj = NULL;
  }
}

void PythonQtPrivate::destroyedSignalEmitter(QObject* obj)
{
  _signalReceivers.take(obj);
}

bool PythonQt::handleError()
{
  bool flag = false;
  if (PyErr_Occurred()) {
    
    // currently we just print the error and the stderr handler parses the errors
    PyErr_Print();
    
    /*
    // EXTRA: the format of the ptype and ptraceback is not really documented, so I use PyErr_Print() above
    PyObject *ptype;
    PyObject *pvalue;
    PyObject *ptraceback;
    PyErr_Fetch( &ptype, &pvalue, &ptraceback);
    
      Py_XDECREF(ptype);
      Py_XDECREF(pvalue);
      Py_XDECREF(ptraceback);
    */
    PyErr_Clear();
    flag = true;
  }
  return flag;
}

void PythonQt::overwriteSysPath(const QStringList& paths)
{
  PythonQtObjectPtr sys;
  sys.setNewRef(PyImport_ImportModule("sys"));
  PyModule_AddObject(sys, "path", PythonQtConv::QStringListToPyList(paths));
}

void PythonQt::setModuleImportPath(PyObject* module, const QStringList& paths)
{
  PyModule_AddObject(module, "__path__", PythonQtConv::QStringListToPyList(paths));
}

void PythonQt::stdOutRedirectCB(const QString& str)
{
  emit PythonQt::self()->pythonStdOut(str);
}

void PythonQt::stdErrRedirectCB(const QString& str)
{
  emit PythonQt::self()->pythonStdErr(str);
}



static PyMethodDef PythonQtMethods[] = {
  {NULL, NULL, 0, NULL}
};

void PythonQt::initPythonQtModule(bool redirectStdOut)
{
  _p->_pythonQtModule.setNewRef(Py_InitModule("PythonQt", PythonQtMethods));
  _p->_qtNamespace = new PythonQtClassInfo(&staticQtMetaObject);
  PyModule_AddObject(_p->_pythonQtModule, "Qt", (PyObject*)_p->createNewPythonQtMetaObjectWrapper(_p->_qtNamespace));
  
  if (redirectStdOut) {
    PythonQtObjectPtr sys;
    PythonQtObjectPtr out;
    PythonQtObjectPtr err;
    sys.setNewRef(PyImport_ImportModule("sys"));
    // create a redirection object for stdout and stderr
    out = PythonQtStdOutRedirectType.tp_new(&PythonQtStdOutRedirectType,NULL, NULL);
    ((PythonQtStdOutRedirect*)out.object())->_cb = stdOutRedirectCB;
    err = PythonQtStdOutRedirectType.tp_new(&PythonQtStdOutRedirectType,NULL, NULL);
    ((PythonQtStdOutRedirect*)err.object())->_cb = stdErrRedirectCB;
    // replace the built in file objects with our own objects
    PyModule_AddObject(sys, "stdout", out);
    PyModule_AddObject(sys, "stderr", err);
  }
}

void PythonQt::addVariantWrapper(const char* typeName, QObject* wrapper)
{
  _p->addVariantWrapper(typeName, wrapper);
}


void PythonQtPrivate::addVariantWrapper(const char* typeName, QObject* wrapper)
{
  int type = QMetaType::type(typeName);
  PythonQtClassInfo* info = new PythonQtClassInfo(wrapper->metaObject(), typeName);
  _knownVariantWrappers.insert(type, qMakePair(info, wrapper));
  addDecorators(wrapper, false, true);
  PyModule_AddObject(_pythonQtModule, typeName, (PyObject*)createNewPythonQtMetaObjectWrapper(info));
}

PyObject* PythonQt::helpCalled(PythonQtClassInfo* info)
{ 
  if (_p->_initFlags & ExternalHelp) {
    emit pythonHelpRequest(QByteArray(info->className()));
    return Py_BuildValue("");
  } else {
    return PyString_FromString(info->help().toLatin1().data());
  }
}
