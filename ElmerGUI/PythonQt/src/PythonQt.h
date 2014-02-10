#ifndef _PYTHONQT_H
#define _PYTHONQT_H

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
// \file    PythonQt.h
// \author  Florian Link
// \author  Last changed by $Author: florian $
// \date    2006-05
*/
//----------------------------------------------------------------------------------

#include "PythonQtSystem.h"
#include "PythonQtWrapper.h"
#include "PythonQtVariantWrapper.h"
#include "PythonQtMetaObjectWrapper.h"
#include "PythonQtSlot.h"
#include "PythonQtObjectPtr.h"
#include <QObject>
#include <QVariant>
#include <QList>
#include <QHash>
#include <QByteArray>
#include <QStringList>
#include <QtDebug>
#include <iostream>


class PythonQtClassInfo;
class PythonQtPrivate;
class PythonQtMethodInfo;
class PythonQtSignalReceiver;
class PythonQtImportFileInterface;
class PythonQtCppWrapperFactory;
class PythonQtConstructorHandler;

//! the main interface to the Python Qt binding, realized as a singleton
class PYTHONQT_EXPORT PythonQt : public QObject {

  Q_OBJECT

public:
  enum InitFlags {
    RedirectStdOut = 1,   //!<< sets if the std out/err is redirected to pythonStdOut() and pythonStdErr() signals
    IgnoreSiteModule = 2, //!<< sets if Python should ignore the site module
    ExternalHelp = 4      //!<< sets if help() calls on PythonQt modules are forwarded to the pythonHelpRequest() signal
  };

  //! initialize the python qt binding (flags are a or combination of InitFlags)
  static void init(int flags = IgnoreSiteModule | RedirectStdOut);

  //! cleanup
  static void cleanup();

  //! get the singleton instance
  static PythonQt* self() { return _self; }

  //-----------------------------------------------------------------------------
  // Public API:

  //! defines the object types for introspection
  enum ObjectType {
    Class,
    Function,
    Variable,
    Module,
    Anything,
    CallOverloads
  };

  //! overwrite the python sys path (call this directly after PythonQt::init() if you want to change the std python sys path)
  void overwriteSysPath(const QStringList& paths);

  //! sets the __path__ list of a module to the given list (important for local imports)
  void setModuleImportPath(PyObject* module, const QStringList& paths);

  //! get the __main__ module of python
  PythonQtObjectPtr getMainModule();

  //! registers a QObject derived class to PythonQt (this is implicitly called by addObject as well)
  //! All added metaobjects will be visible under the className in the PythonQt module as MetaObjectWrappers and the enums
  //! and constructors (added by addConstructors) will be available.
  /* Since Qt4 does not offer a way to detect if a given classname is derived from QObject and thus has a QMetaObject,
     you MUST register all your QObject derived classes here when you want them to be detected in signal and slot calls */
  void registerClass(const QMetaObject* metaobject);

  //! as an alternative to registerClass, you can tell PythonQt the names of QObject derived classes
  //! and it will register the classes when it first sees a pointer to such a derived class
  void registerQObjectClassNames(const QStringList& names);

  //! this will register CPP classnames as known CPP classes (NOT QObjects) and make their MetaObjectWrapper available in
  //! the PythonQt module. In combination with addConstuctors(), this can be used to create CPP objects from PythonQt
  void registerCPPClassNames(const QStringList& names);

  //! parses the given file and returns the python code object, this can then be used to call evalCode()
  PythonQtObjectPtr parseFile(const QString& filename);

  //! evaluates the given code and returns the result value (use Py_Compile etc. to create pycode from string)
  //! If pycode is NULL, a python error is printed.
  QVariant evalCode(PyObject* module, PyObject* pycode);

  //! evaluates the given script code and returns the result value
  QVariant evalScript(PyObject* module, const QString& script, int start = Py_file_input);

  //! evaluates the given script code from file
  void evalFile(PyObject* module, const QString& filename);

  //@{ Signal handlers

  //! add a signal handler to the given \c signal of \c obj  and connect it to a callable \c objectname in module
  bool addSignalHandler(QObject* obj, const char* signal, PyObject* module, const QString& objectname);

  //! remove a signal handler from the given \c signal of \c obj
  bool removeSignalHandler(QObject* obj, const char* signal, PyObject* module, const QString& objectname);

  //! add a signal handler to the given \c signal of \c obj  and connect it to a callable \c receiver
  bool addSignalHandler(QObject* obj, const char* signal, PyObject* receiver);

  //! remove a signal handler from the given \c signal of \c obj
  bool removeSignalHandler(QObject* obj, const char* signal, PyObject* receiver);

  //@}

  //@{ Variable access

  //! add the given \c object to the \c module as a variable with \c name (it can be removed via clearVariable)
  void addObject(PyObject* module, const QString& name, QObject* object);

  //! add the given variable to the module
  void addVariable(PyObject* module, const QString& name, const QVariant& v);

  //! remove the given variable
  void removeVariable(PyObject* module, const QString& name);

  //! get the variable with the \c name of the \c module, returns an invalid QVariant on error
  QVariant getVariable(PyObject* module, const QString& name);

  //! read vars etc. in scope of a module, optional looking inside of an object \c objectname
  QStringList introspection(PyObject* module, const QString& objectname, ObjectType type);

  //! returns the found callable object or NULL
  //! @return new reference
  PythonQtObjectPtr lookupCallable(PyObject* module, const QString& name);

  //@}

  //@{ Calling of python callables

  //! call the given python method, returns the result converted to a QVariant
  QVariant call(PyObject* module, const QString& callable, const QVariantList& args);

  //@}

  //@{ Decorations, constructors, wrappers...


  //! add an object whose slots will be used as decorator slots for
  //! other QObjects or CPP classes. The slots need to follow the
  //! convention that the first argument is a pointer to the wrapped object.
  //! (ownership is passed to PythonQt)
  /*!
  Example:

  A slot with the signature

  \code
  bool doSomething(QWidget* w, int a)
  \endcode

  will extend QWidget instances (and derived classes) with a "bool doSomething(int a)" slot
  that will be called with the concrete instance as first argument.
  So in Python you can now e.g. call

  \code
  someWidget.doSomething(12)
  \endcode

  without QWidget really having this method. This allows to easily make normal methods
  of Qt classes callable by forwarding them with such decorator slots
  or to make CPP classes (which are not derived from QObject) callable from Python.
  */
  void addInstanceDecorators(QObject* o);

  //! add an object whose slots will be used as decorator slots for
  //! class objects (ownership is passed to PythonQt)
  /*!
  The slots need to follow the following convention:
  - SomeClass* new_SomeClass(...)
  - QVariant new_SomeClass(...)
  - void delete_SomeClass(SomeClass*)
  - ... static_SomeClass_someName(...)

  This will add:
  - a constructor
  - a constructor which generates a QVariant
  - a destructor (only useful for CPP objects)
  - a static decorator slot which will be available on the MetaObject (visible in PythonQt module)

  */
  void addClassDecorators(QObject* o);

  //! this will add the object both as class and instance decorator (ownership is passed to PythonQt)
  void addDecorators(QObject* o);

  //! add a wrapper object for the given QMetaType typeName, also does an addClassDecorators() to add constructors for variants
  //! (ownership of wrapper is passed to PythonQt)
  /*! Make sure that you have done a qRegisterMetaType first, if typeName is a user type!

  This will add a wrapper object that is used to make calls to the given classname \c typeName.
  All slots that take a pointer to typeName as the first argument will be callable from Python on
  a variant object that contains such a type.
  */
  void addVariantWrapper(const char* typeName, QObject* wrapper);

  //! add the given factory to PythonQt (ownership stays with caller)
  void addWrapperFactory(PythonQtCppWrapperFactory* factory);

  //! add the given constructor handler to PythonQt (ownership stays with caller)
  void addConstructorHandler(PythonQtConstructorHandler* handler);

  //! get list of constructor handlers
  const QList<PythonQtConstructorHandler*>& constructorHandlers();

  //@}

  //@{ Custom importer (to replace internal import implementation of python)

  //! replace the internal import implementation and use the supplied interface to load files (both py and pyc files)
  //! (this method should be called directly after initialization of init() and before calling overwriteSysPath().
  //! It can only be called once, further calls will be ignored silently. (ownership stays with caller)
  void setImporter(PythonQtImportFileInterface* importInterface);

  //! set paths that the importer should ignore
  void setImporterIgnorePaths(const QStringList& paths);

  //! get paths that the importer should ignore
  const QStringList& getImporterIgnorePaths();

  //@}

  //! get access to internal data (should not be used on the public API, but is used by some C functions)
  static PythonQtPrivate* priv() { return _self->_p; }

  //! get access to the file importer (if set)
  static PythonQtImportFileInterface* importInterface();

  //! handle a python error, call this when a python function fails. If no error occurred, it returns false.
  //! The error is currently just output to the python stderr, future version might implement better trace printing
  bool handleError();

signals:
  //! emitted when python outputs something to stdout (and redirection is turned on)
  void pythonStdOut(const QString& str);
  //! emitted when python outputs something to stderr (and redirection is turned on)
  void pythonStdErr(const QString& str);

  //! emitted when help() is called on a PythonQt object and \c ExternalHelp is enabled
  void pythonHelpRequest(const QByteArray& cppClassName);


public:
  //! called by internal help methods
  PyObject* helpCalled(PythonQtClassInfo* info);

private:
  void initPythonQtModule(bool redirectStdOut);

  //! callback for stdout redirection, emits pythonStdOut signal
  static void stdOutRedirectCB(const QString& str);
  //! callback for stderr redirection, emits pythonStdErr signal
  static void stdErrRedirectCB(const QString& str);

  //! returns the found object or NULL
  //! @return new reference
  PythonQtObjectPtr lookupObject(PyObject* module, const QString& name);

  //! get (and create if not available) the signal receiver of that QObject, signal receiver is made child of the passed \c obj
  PythonQtSignalReceiver* getSignalReceiver(QObject* obj);

  PythonQt(int flags);
  ~PythonQt();

  static PythonQt* _self;

  PythonQtPrivate* _p;

};

//! internal PythonQt details
class PythonQtPrivate : public QObject {

  Q_OBJECT

public:
  PythonQtPrivate();
  ~PythonQtPrivate();

  //! returns if the id is the id for PythonQtObjectPtr
  bool isPythonQtObjectPtrMetaId(int id) { return _PythonQtObjectPtr_metaId == id; }

  //! remove the wrapper ptr again
  void removeWrapperPointer(void* obj) { _wrappedObjects.take(obj); }

  //! wrap the given QObject into a Python object (or return existing wrapper!)
  PyObject* wrapQObject(QObject* obj);

  //! wrap the given ptr into a Python object (or return existing wrapper!) if there is a known QObject of that name or a known wrapper in the factory
  PyObject* wrapPtr(void* ptr, const QByteArray& name);

  //! registers a QObject derived class to PythonQt (this is implicitly called by addObject as well)
  /* Since Qt4 does not offer a way to detect if a given classname is derived from QObject and thus has a QMetaObject,
     you MUST register all your QObject derived classes here when you want them to be detected in signal and slot calls */
  void registerClass(const QMetaObject* metaobject);

  //! as an alternative to registerClass, you can tell PythonQt the names of QObject derived classes
  //! and it will register the classes when it first sees a pointer to such a derived class
  void registerQObjectClassNames(const QStringList& names);

  //! add a decorator object
  void addDecorators(QObject* o, bool instanceDeco, bool classDeco);

  //! add a wrapper object for the given qvariant, also does an addConstructors() to add constructors for variants
  void addVariantWrapper(const char* typeName, QObject* wrapper);

  //! get list of all slots that are available as decorator slots
  QList<PythonQtSlotInfo*> getDecoratorSlots(const QByteArray& className);

  //! check if the enum is either part of the \c meta class or contains a scope and is
  //! an enum of another known metaobject (and as last resort, of the Qt namespace)
  bool isEnumType(const QMetaObject* meta, const QByteArray& name);

  //! helper method that creates a PythonQtMetaObjectWrapper object
  PythonQtMetaObjectWrapper* createNewPythonQtMetaObjectWrapper(PythonQtClassInfo* info);

  //! helper method that creates a PythonQtWrapper object and registers it in the object map
  PythonQtWrapper* createNewPythonQtWrapper(QObject* obj, PythonQtClassInfo* info, void* wrappedPtr = NULL);

  //! helper method that creates a PythonQtVariantWrapper object
  PythonQtVariantWrapper* createNewPythonQtVariantWrapper(const QVariant& variant);

  //! get the class info for a meta object (if available)
  PythonQtClassInfo* getClassInfo(const QMetaObject* meta) { return _knownQtClasses.value(meta->className()); }

  //! get the constructor slot for the given classname
  PythonQtSlotInfo* getConstructorSlot(const QByteArray& className) { return _constructorSlots.value(className); }

  //! get the destructor slot for the given classname
  PythonQtSlotInfo* getDestructorSlot(const QByteArray& className) { return _destructorSlots.value(className); }

protected slots:
  //! called when a wrapped QObject is destroyed
  void wrappedObjectDestroyed(QObject* obj);

  //! called when a signal emitting QObject is destroyed to remove the signal handler from the hash map
  void destroyedSignalEmitter(QObject* obj);

private:

  //! stores pointer to PyObject mapping of wrapped QObjects AND C++ objects
  QHash<void* , PythonQtWrapper *>       _wrappedObjects;

  //! stores the meta info of known Qt classes
  QHash<QByteArray, PythonQtClassInfo *> _knownQtClasses;

  //! stores the meta info of known Qt classes
  QHash<QByteArray, PythonQtClassInfo *> _knownQtWrapperClasses;

  //! stores the meta info of known Qt C++ wrapper classes
  QMultiHash<QByteArray, PythonQtSlotInfo *> _knownQtDecoratorSlots;

  //! names of qobject derived classes that can be casted to qobject savely
  QHash<QByteArray, bool> _knownQObjectClassNames;

  //! stores signal receivers for QObjects
  QHash<QObject* , PythonQtSignalReceiver *> _signalReceivers;

  //! the PythonQt python module
  PythonQtObjectPtr _pythonQtModule;

  //! the importer interface (if set)
  PythonQtImportFileInterface* _importInterface;

  QStringList _importIgnorePaths;

  //! the cpp object wrapper factories
  QList<PythonQtCppWrapperFactory*> _cppWrapperFactories;

  //! the cpp object wrapper factories
  QList<PythonQtConstructorHandler*> _constructorHandlers;

  QHash<QByteArray , PythonQtSlotInfo *> _constructorSlots;
  QHash<QByteArray , PythonQtSlotInfo *> _destructorSlots;

  QHash<int , QPair<PythonQtClassInfo*, QObject*> > _knownVariantWrappers;

  PythonQtClassInfo* _qtNamespace;

  int _initFlags;
  int _PythonQtObjectPtr_metaId;

  friend class PythonQt;
};

#endif
