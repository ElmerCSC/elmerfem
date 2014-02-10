#ifndef _PYTHONQTDOC_H
#define _PYTHONQTDOC_H

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
// \file    PythonQtDoc.h
// \author  Florian Link
// \author  Last changed by $Author: florian $
// \date    2006-10
*/
//----------------------------------------------------------------------------------

/*!
\mainpage

 \section Introduction

 \b PythonQt is a dynamic Python (http://www.python.org) binding for Qt (http://www.trolltech.com).
 It offers an easy way to embedd the Python scripting language into
 your Qt applications. It makes heavy use of the QMetaObject system and thus requires Qt4.x.

 In contrast to <a href="http://www.riverbankcomputing.co.uk/pyqt/">PyQt</a> , PythonQt is \b not a complete
 Python wrapper around the complete Qt functionality. So if you are looking for a way to
 write complete applications in Python using the Qt GUI, you should use PyQt.

 If you are looking for a simple way to embed the Python language into your Qt Application
 and to script parts of your application via Python, PythonQt is the way to go!

 PythonQt is a stable library that was developed to make the Image Processing and Visualization platform MeVisLab (http://www.mevislab.de)
 scriptable from Python.

 \section Licensing

 PythonQt is distributed under the LGPL license.

 \section Download

 PythonQt is hosted on SourceForge at http://sourceforge.net/projects/pythonqt , you can access it via SVN
 or download a tarball.

 \section Features

 - Access all \b slots, \b properties, children and registered enums of any QObject derived class from Python
 - Connecting Qt Signals to Python functions (both from within Python and from C++)
 - Wrapping of C++ objects (which are not derived from QObject) via PythonQtCPPWrapperFactory
 - Extending C++ and QObject derived classes with additional slots, static methods and constructors (see Decorators)
 - StdOut/Err redirection to Qt signals instead of cout
 - Interface for creating your own \c import replacement, so that Python scripts can be e.g. signed/verified before they are executed (PythonQtImportInterface)
 - Mapping of plain-old-datatypes and ALL QVariant types to and from Python
 - Support for wrapping of user QVariant types which are registerd via QMetaType
 - Support for Qt namespace (with all enumerators)
 - All PythonQt wrapped objects support the dir() statement, so that you can see easily which attributes a QObject, CPP object or QVariant has
 - No preprocessing/wrapping tool needs to be started, PythonQt can script any QObject without prior knowledge about it (except for the MetaObject information from the \b moc)

 \section Non-Features

 Features that PythonQt does NOT support (and will not support):

 - you can not derive from QObjects inside of Python, this would require wrapper generation like PyQt does
 - you can only script QObject derived classes, for normal C++ classes you need to create a PythonQtCPPWrapperFactory and adequate wrapper classes or add decorator slots
 - you can not access normal member functions of QObjects, only slots and properties, because the \b moc does not store normal member functions in the MetaObject system

 \section Interface

 The main interface to PythonQt is the PythonQt singleton.
 PythonQt needs to be initialized via PythonQt::init() once.
 Afterwards you communicate with the singleton via PythonQt::self().
 PythonQt offers a default binding for the complete QWidget set, which
 needs to be enabled via PythonQtGui::init().


 \section Datatype Datatype Mapping

  The following table shows the mapping between Python and Qt objects:
  <table>
  <tr><th>Qt/C++</th><th>Python</th></tr>
  <tr><td>bool</td><td>bool</td></tr>
  <tr><td>double</td><td>float</td></tr>
  <tr><td>float</td><td>float</td></tr>
  <tr><td>char/uchar,int/uint,short,ushort,QChar</td><td>integer</td></tr>
  <tr><td>long</td><td>integer</td></tr>
  <tr><td>ulong,longlong,ulonglong</td><td>long</td></tr>
  <tr><td>QString</td><td>unicode string</td></tr>
  <tr><td>QByteArray</td><td>str</td></tr>
  <tr><td>char*</td><td>str</td></tr>
  <tr><td>QStringList</td><td>tuple of unicode strings</td></tr>
  <tr><td>QVariantList</td><td>tuple of objects</td></tr>
  <tr><td>QVariantMap</td><td>dict of objects</td></tr>
  <tr><td>QVariant</td><td>depends on type, see below</td></tr>
  <tr><td>QSize, QRect and all other standard Qt QVariants</td><td>variant wrapper that supports complete API of the respective Qt classes</td></tr>
  <tr><td>OwnRegisteredMetaType</td><td>variant wrapper, optionally with a wrapper provided by addVariantWrapper()</td></tr>
  <tr><td>EnumType</td><td>integer (all enums that are known via the moc and the Qt namespace are supported)</td></tr>
  <tr><td>QObject (and derived classes)</td><td>QObject wrapper</td></tr>
  <tr><td>C++ object</td><td>CPP wrapper, either wrapped via PythonQtCPPWrapperFactory or just decorated with decorators</td></tr>
  <tr><td>PyObject</td><td>PyObject</td></tr>
  </table>

  PyObject is passed as simple pointer, which allows to pass/return any Python Object directly to/from
  a Qt slot.
  QVariants are mapped recursively as given above, e.g. a dictionary can
  contain lists of dictionaries of doubles.
  For example a QVariant of type "String" is mapped to a python unicode string.
  All Qt QVariant types are implemented, PythonQt supports the complete Qt API for these object.

 \section QObject QObject Wrapping

 All classes derived from QObject are automatically wrapped with a python wrapper class
 when they become visible to the Python interpreter. This can happen via
 - the PythonQt::addObject() method
 - when a Qt \b slot returns a QObject derived object to python
 - when a Qt \b signal contains a QObject and is connected to a python function

 It is important that you call PythonQt::registerClass() for any QObject derived class
 that may become visible to Python, except when you add it via PythonQt::addObject().
 This will register the complete parent hierachy of the registered class, so that
 when you register e.g. a QPushButton, QWidget will be registered as well (and all intermediate
 parents).

 From Python, you can talk to the returned QObjects in a natural way by calling
 their slots and receiving the return values. You can also read/write all
 properties of the objects as if they where normal python properties.

 In addition to this, the wrapped objects support
 - className() - returns a string that reprents the classname of the QObject
 - help() - shows all properties, slots, enums, decorator slots and constructors of the object, in a printable form
 - connect(signal, function) - connect the signal of the given object to a python function
 - connect(signal, qobject, slot) - connect the signal of the given object to a slot of another QObject
 - disconnect(signal, function) - disconnect the signal of the given object from a python function
 - disconnect(signal, qobject, slot) - disconnect the signal of the given object from a slot of another QObject
 - children() - returns the children of the object
 - setParent(QObject) - set the parent
 - QObject* parent() - get the parent

 The below example shows how to connect signals in Python:

 \code
 # define a signal handler function
 def someFunction(flag):
   print flag

 # button1 is a QPushButton that has been added to Python via addObject()
 # connect the clicked signal to a python function:
 button1.connect("clicked(bool)", someFunction)

 \endcode

\section CPP CPP Wrapping

You can create dedicated wrapper QObject for any C++ class. This is done by deriving from PythonQtCPPWrapperFactory
and adding your factory via addWrapperFactory(). Whenever PythonQt encounters a CPP pointer (e.g. on a slot or signal)
and it does not known it as a QObject derived class, it will create a generic CPP wrapper. So even unknown C++ objects
can be passed through Python. If the wrapper factory supports the CPP class, a QObject wrapper will be created for each
instance that enters Python. An alternative to a complete wrapper via the wrapper factory are decorators, see \ref Decorators

\section MetaObject Meta Object/Class access

For each known CPP class, QObject derived class and QVariant type, PythonQt provides a Meta class. These meta classes are visible
inside of the "PythonQt" python module.

A Meta class supports:

- access to all declared enum values
- constructors
- static decorator slots
- help() and className()

From within Python, you can import the module "PythonQt" to access these meta objects and the Qt namespace.

\code
from PythonQt import *

# namespace access:
print Qt.AlignLeft

# constructors
a = QSize(12,13)
b = QFont()

# static method
QDate.currentDate()

# enum value
QFont.UltraCondensed

\endcode

\section Decorators Decorator slots

PythonQt introduces a new generic approach to extend any wrapped QObject or CPP object with

- constructors
- destructors (for CPP objects)
- additional slots
- static slots (callable on both the Meta object and the instances)

The idea behind decorators is that we wanted to make it as easy as possible to extend
wrapped objects. Since we already have an implementation for invoking any Qt Slot from
Python, it looked promising to use this approach for the extension of wrapped objects as well.
This avoids that the PythonQt user needs to care about how Python arguments are mapped from/to
Qt when he wants to create static methods, constructors and additional member functions.

The basic idea about decorators is to create a QObject derived class that implements slots
which take one of the above roles (e.g. constructor, destructor etc.) via a naming convention.
These slots are then assigned to other classes via the naming convention.

- QVariant new_SomeClassName(...) - defines a constructor for "SomeClassName" that returns a QVariant
- SomeClassName* new_SomeClassName(...) - defines a constructor for "SomeClassName" that returns a new object of type SomeClassName (where SomeClassName can be any CPP class, not just QObject classes)
- void delete_SomeClassName(SomeClassName* o) - defines a destructor, which should delete the passed in object o
- anything static_SomeClassName_someMethodName(...) - defines a static method that is callable on instances and the meta class
- anything someMethodName(SomeClassName* o, ...) - defines a slot that will be available on SomeClassName instances (and derived instances). When such a slot is called the first argument is the pointer to the instance and the rest of the arguments can be used to make a call on the instance.

The below example shows all kinds of decorators in action:

\code

// an example CPP object
class YourCPPObject {
public:
  YourCPPObject(int arg1, float arg2) { a = arg1; b = arg2; }

  float doSomething(int arg1) { return arg1*a*b; };

  private:

  int a;
  float b;
};

// an example decorator
class ExampleDecorator : public QObject
{
  Q_OBJECT

public slots:
  // add a constructor to QSize variant that takes a QPoint
  QVariant new_QSize(const QPoint& p) { return QSize(p.x(), p.y()); }

  // add a constructor for QPushButton that takes a text and a parent widget
  QPushButton* new_QPushButton(const QString& text, QWidget* parent=NULL) { return new QPushButton(text, parent); }

  // add a constructor for a CPP object
  YourCPPObject* new_YourCPPObject(int arg1, float arg2) { return new YourCPPObject(arg1, arg2); }

  // add a destructor for a CPP object
  void delete_YourCPPObject(YourCPPObject* obj) { delete obj; }

  // add a static method to QWidget
  QWidget* static_QWidget_mouseGrabber() { return QWidget::mouseGrabber(); }

  // add an additional slot to QWidget (make move() callable, which is not declared as a slot in QWidget)
  void move(QWidget* w, const QPoint& p) { w->move(p); }

  // add an additional slot to QWidget, overloading the above move method
  void move(QWidget* w, int x, int y) { w->move(x,y); }

  // add a method to your own CPP object
  int doSomething(YourCPPObject* obj, int arg1) { return obj->doSomething(arg1); }
};

...

PythonQt::self()->addDecorators(new ExampleDecorator());
PythonQt::self()->registerClass(&QPushButton::staticMetaObject);
PythonQt::self()->registerCPPClassNames(QStringList() << "YourCPPObject");

\endcode

After you have registered an instance of the above ExampleDecorator, you can do the following from Python
(all these calls are mapped to the above decorator slots):

\code
from PythonQt import *

# call our new constructor of QSize
size = QSize(QPoint(1,2));

# call our new QPushButton constructor
button = QPushButton("sometext");

# call the move slot (overload1)
button.move(QPoint(0,0))

# call the move slot (overload2)
button.move(0,0)

# call the static method
grabber = QWidget.mouseWrapper();

# create a CPP object via constructor
yourCpp = YourCPPObject(1,11.5)

# call the wrapped method on CPP object
print yourCpp.doSomething(1);

# destructor will be called:
yourCpp = None

\endcode

 \section Building

 PythonQt requires at least Qt 4.2.2 (or higher) and Python 2.3, 2.4 and 2.5 on Windows, Linux and MacOS X.
 To compile PythonQt, you will need a python developer installation which includes Python's header files and
 the python2x.[lib | dll | so | dynlib].

 For building PythonQt, you will need to set some environment variables:
 \b PYTHON_ROOT should point to the Python sources/headers.
 \b PYTHON_LIB should point to where the Python library files are located.
 \b PYTHONQT_ROOT should point to the root directory of PythonQt.

 Run qmake on PythonQt.pro to generate a project file for your system and then build it.

 \section Tests

 There is a unit test that tests most features of PythonQt, see the \b tests subdirectory for details.

 \section Examples

 Examples are available in the \b examples directory. The PyScriptingConsole implements a simple
 interactive scripting console that shows how to script a simple application.

 The following shows how to integrate PythonQt into you Qt application:

 \code
 #include "PythonQt.h"
 #include <QApplication>
 ...

 int main( int argc, char **argv )
 {

  QApplication qapp(argc, argv);

  // init PythonQt and Python itself
  PythonQt::init(PythonQt::IgnoreSiteModule | PythonQt::RedirectStdOut);

  // get a smart pointer to the __main__ module of the Python interpreter
  PythonQtObjectPtr mainContext = PythonQt::self()->getMainModule();

  // add a QObject as variable of name "example" to the namespace of the __main__ module
  PyExampleObject example;
  PythonQt::self()->addObject(mainContext, "example", &example);

  // register all other QObjects that you want to script and that are returned by your API
  PythonQt::self()->registerClass(&QMainWindow::staticMetaObject);
  PythonQt::self()->registerClass(&QPushButton::staticMetaObject);
  ...

  // do something
  PythonQt::self()->runScript(mainContext, "print example\n");
  PythonQt::self()->runScript(mainContext, "def multiply(a,b):\n  return a*b;\n");
  QVariantList args;
  args << 42 << 47;
  QVariant result = PythonQt::self()->call(mainContext,"multiply", args);
  ...
 \endcode


  \section TODOs

  - improve qmake profiles for non mevis users
  - add more information on how to distribute an application that uses PythonQt, including the Python distribution

*/
