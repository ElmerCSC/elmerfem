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
// \file    PythonQtImporter.h
// \author  Florian Link
// \author  Last changed by $Author: florian $
// \date    2006-05
*/
// This module was inspired by the zipimport.c module of the original
// Python distribution. Most of the functions are identical or slightly
// modified to do all the loading of Python files via an external file interface.
// In contrast to zipimport.c, this module also writes *.pyc files
// automatically if it has write access/is not inside of a zip file.
//----------------------------------------------------------------------------------

#include "PythonQtImporter.h"
#include "PythonQtImportFileInterface.h"
#include "PythonQt.h"
#include <QFile>
#include <QFileInfo>

#define IS_SOURCE   0x0
#define IS_BYTECODE 0x1
#define IS_PACKAGE  0x2

struct st_mlab_searchorder {
  char suffix[14];
  int type;
};

/* mlab_searchorder defines how we search for a module in the Zip
   archive: we first search for a package __init__, then for
   non-package .pyc, .pyo and .py entries. The .pyc and .pyo entries
   are swapped by initmlabimport() if we run in optimized mode. Also,
   '/' is replaced by SEP there. */
 struct st_mlab_searchorder mlab_searchorder[] = {
  {"/__init__.pyc", IS_PACKAGE | IS_BYTECODE},
  {"/__init__.pyo", IS_PACKAGE | IS_BYTECODE},
  {"/__init__.py", IS_PACKAGE | IS_SOURCE},
  {".pyc", IS_BYTECODE},
  {".pyo", IS_BYTECODE},
  {".py", IS_SOURCE},
  {"", 0}
};

extern PyTypeObject PythonQtImporter_Type;
PyObject *PythonQtImportError;

QString PythonQtImport::getSubName(const QString& str)
{
  int idx = str.lastIndexOf('.');
  if (idx!=-1) {
    return str.mid(idx+1);
  } else {
    return str;
  }
}

PythonQtImport::module_info PythonQtImport::getModuleInfo(PythonQtImporter* self, const QString& fullname)
{
  QString subname;
  struct st_mlab_searchorder *zso;

  subname = getSubName(fullname);
  QString path = *self->_path + "/" + subname;

  QString test;
  for (zso = mlab_searchorder; *zso->suffix; zso++) {
    test = path + zso->suffix;
    if (PythonQt::importInterface()->exists(test)) {
      if (zso->type & IS_PACKAGE)
        return MI_PACKAGE;
      else
        return MI_MODULE;
    }
  }
  return MI_NOT_FOUND;
}


/* PythonQtImporter.__init__
  Just store the path argument
*/
int PythonQtImporter_init(PythonQtImporter *self, PyObject *args, PyObject *kwds)
{
  self->_path = NULL;

  const char* path;
  if (!PyArg_ParseTuple(args, "s",
    &path))
    return -1;

  if (PythonQt::importInterface()->exists(path)) {
    //qDebug("path %s", path);
    QString p(path);
    const QStringList& ignorePaths = PythonQt::self()->getImporterIgnorePaths();
    foreach(QString a, ignorePaths) {
      if (a==p) {
        PyErr_SetString(PythonQtImportError,
          "path ignored");
        return -1;
      }
    }

    self->_path = new QString(p);

    //mlabDebugConst("MLABPython", "PythonQtImporter init: " << *self->_path);

    return 0;
  } else {
    PyErr_SetString(PythonQtImportError,
        "path does not exist error");
    return -1;
  }
}

void
PythonQtImporter_dealloc(PythonQtImporter *self)
{
  // free the stored path
  if (self->_path) delete self->_path;
  // free ourself
  self->ob_type->tp_free((PyObject *)self);
}


/* Check whether we can satisfy the import of the module named by
   'fullname'. Return self if we can, None if we can't. */
PyObject *
PythonQtImporter_find_module(PyObject *obj, PyObject *args)
{
  PythonQtImporter *self = (PythonQtImporter *)obj;
  PyObject *path = NULL;
  char *fullname;

  if (!PyArg_ParseTuple(args, "s|O:PythonQtImporter.find_module",
            &fullname, &path))
    return NULL;

//  mlabDebugConst("MLABPython", "FindModule " << fullname << " in " << *self->_path);

  PythonQtImport::module_info info = PythonQtImport::getModuleInfo(self, fullname);
  if (info == PythonQtImport::MI_MODULE || info == PythonQtImport::MI_PACKAGE) {
    Py_INCREF(self);
    return (PyObject *)self;
  } else {
    Py_INCREF(Py_None);
    return Py_None;
  }
}

/* Load and return the module named by 'fullname'. */
PyObject *
PythonQtImporter_load_module(PyObject *obj, PyObject *args)
{
  PythonQtImporter *self = (PythonQtImporter *)obj;
  PyObject *code, *mod, *dict;
  char *fullname;
  QString modpath;
  int ispackage;

  if (!PyArg_ParseTuple(args, "s:PythonQtImporter.load_module",
            &fullname))
    return NULL;

  code = PythonQtImport::getModuleCode(self, fullname, &ispackage, modpath);
  if (code == NULL)
    return NULL;

  mod = PyImport_AddModule(fullname);
  if (mod == NULL) {
    Py_DECREF(code);
    return NULL;
  }
  dict = PyModule_GetDict(mod);

  if (PyDict_SetItemString(dict, "__loader__", (PyObject *)self) != 0)
    goto error;

  if (ispackage) {
    PyObject *pkgpath, *fullpath;
    QString subname = PythonQtImport::getSubName(fullname);
    int err;

    fullpath = PyString_FromFormat("%s%c%s",
          self->_path->toLatin1().constData(),
          SEP,
          subname.toLatin1().constData());
    if (fullpath == NULL)
      goto error;

    pkgpath = Py_BuildValue("[O]", fullpath);
    Py_DECREF(fullpath);
    if (pkgpath == NULL)
      goto error;
    err = PyDict_SetItemString(dict, "__path__", pkgpath);
    Py_DECREF(pkgpath);
    if (err != 0)
      goto error;
  }
  mod = PyImport_ExecCodeModuleEx(fullname, code, (char*)modpath.toLatin1().data());
  Py_DECREF(code);
  if (Py_VerboseFlag)
    PySys_WriteStderr("import %s # loaded from %s\n",
          fullname, (char*)modpath.toLatin1().data());
  return mod;
error:
  Py_DECREF(code);
  Py_DECREF(mod);
  return NULL;
}


PyObject *
PythonQtImporter_get_data(PyObject *obj, PyObject *args)
{
  // EXTRA, NOT YET IMPLEMENTED
  return NULL;
}

PyObject *
PythonQtImporter_get_code(PyObject *obj, PyObject *args)
{
  PythonQtImporter *self = (PythonQtImporter *)obj;
  char *fullname;

  if (!PyArg_ParseTuple(args, "s:PythonQtImporter.get_code", &fullname))
    return NULL;

  QString notused;
  return PythonQtImport::getModuleCode(self, fullname, NULL, notused);
}

PyObject *
PythonQtImporter_get_source(PyObject *obj, PyObject *args)
{
  // EXTRA, NOT YET IMPLEMENTED
/*
  PythonQtImporter *self = (PythonQtImporter *)obj;
  PyObject *toc_entry;
  char *fullname, *subname, path[MAXPATHLEN+1];
  int len;
  enum module_info mi;

  if (!PyArg_ParseTuple(args, "s:PythonQtImporter.get_source", &fullname))
    return NULL;

  mi = get_module_info(self, fullname);
  if (mi == MI_ERROR)
    return NULL;
  if (mi == MI_NOT_FOUND) {
    PyErr_Format(PythonQtImportError, "can't find module '%.200s'",
           fullname);
    return NULL;
  }
  subname = get_subname(fullname);

  len = make_filename(PyString_AsString(self->prefix), subname, path);
  if (len < 0)
    return NULL;

  if (mi == MI_PACKAGE) {
    path[len] = SEP;
    strcpy(path + len + 1, "__init__.py");
  }
  else
    strcpy(path + len, ".py");

  toc_entry = PyDict_GetItemString(self->files, path);
  if (toc_entry != NULL)
    return get_data(PyString_AsString(self->archive), toc_entry);

  Py_INCREF(Py_None);
  return Py_None;
*/
  return NULL;
}

PyDoc_STRVAR(doc_find_module,
"find_module(fullname, path=None) -> self or None.\n\
\n\
Search for a module specified by 'fullname'. 'fullname' must be the\n\
fully qualified (dotted) module name. It returns the PythonQtImporter\n\
instance itself if the module was found, or None if it wasn't.\n\
The optional 'path' argument is ignored -- it's there for compatibility\n\
with the importer protocol.");

PyDoc_STRVAR(doc_load_module,
"load_module(fullname) -> module.\n\
\n\
Load the module specified by 'fullname'. 'fullname' must be the\n\
fully qualified (dotted) module name. It returns the imported\n\
module, or raises PythonQtImportError if it wasn't found.");

PyDoc_STRVAR(doc_get_data,
"get_data(pathname) -> string with file data.\n\
\n\
Return the data associated with 'pathname'. Raise IOError if\n\
the file wasn't found.");

PyDoc_STRVAR(doc_get_code,
"get_code(fullname) -> code object.\n\
\n\
Return the code object for the specified module. Raise PythonQtImportError\n\
is the module couldn't be found.");

PyDoc_STRVAR(doc_get_source,
"get_source(fullname) -> source string.\n\
\n\
Return the source code for the specified module. Raise PythonQtImportError\n\
is the module couldn't be found, return None if the archive does\n\
contain the module, but has no source for it.");

PyMethodDef PythonQtImporter_methods[] = {
  {"find_module", PythonQtImporter_find_module, METH_VARARGS,
   doc_find_module},
  {"load_module", PythonQtImporter_load_module, METH_VARARGS,
   doc_load_module},
  {"get_data", PythonQtImporter_get_data, METH_VARARGS,
   doc_get_data},
  {"get_code", PythonQtImporter_get_code, METH_VARARGS,
   doc_get_code},
  {"get_source", PythonQtImporter_get_source, METH_VARARGS,
   doc_get_source},
  {NULL,    NULL} /* sentinel */
};


PyDoc_STRVAR(PythonQtImporter_doc,
"PythonQtImporter(path) -> PythonQtImporter object\n\
\n\
Create a new PythonQtImporter instance. 'path' must be a valid path on disk/or inside of a zip file known to MeVisLab\n\
. Every path is accepted.");

#define DEFERRED_ADDRESS(ADDR) 0

PyTypeObject PythonQtImporter_Type = {
  PyObject_HEAD_INIT(DEFERRED_ADDRESS(&PyType_Type))
  0,
  "PythonQtImport.PythonQtImporter",
  sizeof(PythonQtImporter),
  0,          /* tp_itemsize */
  (destructor)PythonQtImporter_dealloc, /* tp_dealloc */
  0,          /* tp_print */
  0,          /* tp_getattr */
  0,          /* tp_setattr */
  0,          /* tp_compare */
  0,    /* tp_repr */
  0,          /* tp_as_number */
  0,          /* tp_as_sequence */
  0,          /* tp_as_mapping */
  0,          /* tp_hash */
  0,          /* tp_call */
  0,          /* tp_str */
  PyObject_GenericGetAttr,    /* tp_getattro */
  0,          /* tp_setattro */
  0,          /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE ,    /* tp_flags */
  PythonQtImporter_doc,     /* tp_doc */
  0,      /* tp_traverse */
  0,          /* tp_clear */
  0,          /* tp_richcompare */
  0,          /* tp_weaklistoffset */
  0,          /* tp_iter */
  0,          /* tp_iternext */
  PythonQtImporter_methods,     /* tp_methods */
  0,          /* tp_members */
  0,          /* tp_getset */
  0,          /* tp_base */
  0,          /* tp_dict */
  0,          /* tp_descr_get */
  0,          /* tp_descr_set */
  0,          /* tp_dictoffset */
  (initproc)PythonQtImporter_init,    /* tp_init */
  PyType_GenericAlloc,      /* tp_alloc */
  PyType_GenericNew,      /* tp_new */
  PyObject_Del,     /* tp_free */
};


/* Given a buffer, return the long that is represented by the first
   4 bytes, encoded as little endian. This partially reimplements
   marshal.c:r_long() */
long
PythonQtImport::getLong(unsigned char *buf)
{
  long x;
  x =  buf[0];
  x |= (long)buf[1] <<  8;
  x |= (long)buf[2] << 16;
  x |= (long)buf[3] << 24;
#if SIZEOF_LONG > 4
  /* Sign extension for 64-bit machines */
  x |= -(x & 0x80000000L);
#endif
  return x;
}

FILE *
open_exclusive(const QString& filename)
{
#if defined(O_EXCL)&&defined(O_CREAT)&&defined(O_WRONLY)&&defined(O_TRUNC)
  /* Use O_EXCL to avoid a race condition when another process tries to
     write the same file.  When that happens, our open() call fails,
     which is just fine (since it's only a cache).
     XXX If the file exists and is writable but the directory is not
     writable, the file will never be written.  Oh well.
  */
  QFile::remove(filename);

  int fd;
  int flags = O_EXCL|O_CREAT|O_WRONLY|O_TRUNC;
#ifdef O_BINARY
    flags |= O_BINARY;   /* necessary for Windows */
#endif
#ifdef WIN32
  fd = _wopen(filename.ucs2(), flags, 0666);
#else
  fd = open(filename.local8Bit(), flags, 0666);
#endif
  if (fd < 0)
    return NULL;
  return fdopen(fd, "wb");
#else
  /* Best we can do -- on Windows this can't happen anyway */
  return fopen(filename.toLocal8Bit().constData(), "wb");
#endif
}


void PythonQtImport::writeCompiledModule(PyCodeObject *co, const QString& filename, long mtime)
{
  FILE *fp;

  fp = open_exclusive(filename);
  if (fp == NULL) {
    if (Py_VerboseFlag)
      PySys_WriteStderr(
      "# can't create %s\n", filename.toLatin1().constData());
    return;
  }
#if PY_VERSION_HEX < 0x02040000
  PyMarshal_WriteLongToFile(PyImport_GetMagicNumber(), fp);
#else
  PyMarshal_WriteLongToFile(PyImport_GetMagicNumber(), fp, Py_MARSHAL_VERSION);
#endif
  /* First write a 0 for mtime */
#if PY_VERSION_HEX < 0x02040000
  PyMarshal_WriteLongToFile(0L, fp);
#else
  PyMarshal_WriteLongToFile(0L, fp, Py_MARSHAL_VERSION);
#endif
#if PY_VERSION_HEX < 0x02040000
  PyMarshal_WriteObjectToFile((PyObject *)co, fp);
#else
  PyMarshal_WriteObjectToFile((PyObject *)co, fp, Py_MARSHAL_VERSION);
#endif
  if (ferror(fp)) {
    if (Py_VerboseFlag)
      PySys_WriteStderr("# can't write %s\n", filename.toLatin1().constData());
    /* Don't keep partial file */
    fclose(fp);
    QFile::remove(filename);
    return;
  }
  /* Now write the true mtime */
  fseek(fp, 4L, 0);
#if PY_VERSION_HEX < 0x02040000
  PyMarshal_WriteLongToFile(mtime, fp);
#else
  PyMarshal_WriteLongToFile(mtime, fp, Py_MARSHAL_VERSION);
#endif
  fflush(fp);
  fclose(fp);
  if (Py_VerboseFlag)
    PySys_WriteStderr("# wrote %s\n", filename.toLatin1().constData());
//#ifdef macintosh
//  PyMac_setfiletype(cpathname, 'Pyth', 'PYC ');
//#endif
}

/* Given the contents of a .py[co] file in a buffer, unmarshal the data
   and return the code object. Return None if it the magic word doesn't
   match (we do this instead of raising an exception as we fall back
   to .py if available and we don't want to mask other errors).
   Returns a new reference. */
PyObject *
PythonQtImport::unmarshalCode(const QString& path, const QByteArray& data, time_t mtime)
{
  PyObject *code;
  // ugly cast, but Python API is not const safe
  char *buf = (char*) data.constData();
  int size = data.size();

  if (size <= 9) {
    PySys_WriteStderr("# %s has bad pyc data\n",
            path.toLatin1().constData());
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (getLong((unsigned char *)buf) != PyImport_GetMagicNumber()) {
    if (Py_VerboseFlag)
      PySys_WriteStderr("# %s has bad magic\n",
            path.toLatin1().constData());
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (mtime != 0 && !(getLong((unsigned char *)buf + 4) == mtime)) {
    if (Py_VerboseFlag)
      PySys_WriteStderr("# %s has bad mtime\n",
            path.toLatin1().constData());
    Py_INCREF(Py_None);
    return Py_None;
  }

  code = PyMarshal_ReadObjectFromString(buf + 8, size - 8);
  if (code == NULL)
    return NULL;
  if (!PyCode_Check(code)) {
    Py_DECREF(code);
    PyErr_Format(PyExc_TypeError,
         "compiled module %.200s is not a code object",
         path.toLatin1().constData());
    return NULL;
  }
  return code;
}


/* Given a string buffer containing Python source code, compile it
   return and return a code object as a new reference. */
PyObject *
PythonQtImport::compileSource(const QString& path, const QByteArray& data)
{
  PyObject *code;
  QByteArray data1 = data;
// in qt4, data is null terminated
//  data1.resize(data.size()+1);
//  data1.data()[data.size()-1] = 0;
  code = Py_CompileString(data.data(), path.toLatin1().constData(),
        Py_file_input);
  return code;
}


/* Return the code object for the module named by 'fullname' from the
   Zip archive as a new reference. */
PyObject *
PythonQtImport::getCodeFromData(const QString& path, int isbytecode,int ispackage, time_t mtime)
{
  bool hasImporter = PythonQt::importInterface()!=NULL;

  PyObject *code;

  QByteArray qdata;
  if (!hasImporter) {
    QFile file(path);
    QIODevice::OpenMode flags = QIODevice::ReadOnly;
    if (!isbytecode) {
      flags |= QIODevice::Text;
    }
    if (!file.open(flags)) {
      return NULL;
    } 
    qdata = file.readAll();
  } else {
    if (!isbytecode) {
      //    mlabDebugConst("MLABPython", "reading source " << path);
      bool ok;
      qdata = PythonQt::importInterface()->readSourceFile(path, ok);
      if (!ok) {
        //    mlabErrorConst("PythonQtImporter","File could not be verified" << path);
        return NULL;
      }
      if (qdata == " ") {
        qdata.clear();
      }
    } else {
      qdata = PythonQt::importInterface()->readFileAsBytes(path);
    }
  }

  if (isbytecode) {
//    mlabDebugConst("MLABPython", "reading bytecode " << path);
    code = unmarshalCode(path, qdata, mtime);
  }
  else {
  //  mlabDebugConst("MLABPython", "compiling source " << path);
    code = compileSource(path, qdata);
    // save a pyc file if possible
    QDateTime time;
    time = hasImporter?PythonQt::importInterface()->lastModifiedDate(path):QFileInfo(path).lastModified();
    writeCompiledModule((PyCodeObject*)code, path+"c", time.toTime_t());
  }
  return code;
}

time_t
PythonQtImport::getMTimeOfSource(const QString& path)
{
  time_t mtime = 0;
  QString path2 = path;
  path2.truncate(path.length()-1);
  if (PythonQt::importInterface()->exists(path2)) {
    mtime = PythonQt::importInterface()->lastModifiedDate(path2).toTime_t();
  }
  return mtime;
}

/* Get the code object associated with the module specified by
   'fullname'. */
PyObject *
PythonQtImport::getModuleCode(PythonQtImporter *self, char *fullname,
    int *p_ispackage, QString& modpath)
{
  QString subname;
  struct st_mlab_searchorder *zso;

  subname = getSubName(fullname);
  QString path = *self->_path + "/" + subname;

  QString test;
  for (zso = mlab_searchorder; *zso->suffix; zso++) {
    PyObject *code = NULL;
    test = path + zso->suffix;

    if (Py_VerboseFlag > 1)
      PySys_WriteStderr("# trying %s\n",
            test.toLatin1().constData());
    if (PythonQt::importInterface()->exists(test)) {
      time_t mtime = 0;
      int ispackage = zso->type & IS_PACKAGE;
      int isbytecode = zso->type & IS_BYTECODE;

      if (isbytecode)
        mtime = getMTimeOfSource(test);
      if (p_ispackage != NULL)
        *p_ispackage = ispackage;
      code = getCodeFromData(test, isbytecode, ispackage, mtime);
      if (code == Py_None) {
        Py_DECREF(code);
        continue;
      }
      if (code != NULL)
        modpath = test;
      return code;
    }
  }
  PyErr_Format(PythonQtImportError, "can't find module '%.200s'", fullname);

  return NULL;
}

QString PythonQtImport::replaceExtension(const QString& str, const QString& ext)
{
 QString r;
 int i = str.lastIndexOf('.');
 if (i!=-1) {
   r = str.mid(0,i) + "." + ext;
 } else {
   r = str + "." + ext;
 }
 return r;
}

PyObject* PythonQtImport::getCodeFromPyc(const QString& file)
{
  bool hasImporter = PythonQt::importInterface()!=NULL;

  PyObject* code;
  const static QString pycStr("pyc");
  QString pyc = replaceExtension(file, pycStr);
  if ((hasImporter && PythonQt::importInterface()->exists(pyc)) ||
    (!hasImporter && QFile::exists(pyc))) {
    time_t mtime = 0;
    mtime = getMTimeOfSource(pyc);
    code = getCodeFromData(pyc, true, false, mtime);
    if (code != Py_None && code != NULL) {
      return code;
    }
    if (code) {
      Py_DECREF(code);
    }
  }
  code = getCodeFromData(file,false,false,0);
  return code;
}

/* Module init */

PyDoc_STRVAR(mlabimport_doc,
"Imports python files into MeVisLab, completely replaces internal python import");

void PythonQtImport::init()
{
  PyObject *mod;

  if (PyType_Ready(&PythonQtImporter_Type) < 0)
    return;

  /* Correct directory separator */
  mlab_searchorder[0].suffix[0] = SEP;
  mlab_searchorder[1].suffix[0] = SEP;
  mlab_searchorder[2].suffix[0] = SEP;
  if (Py_OptimizeFlag) {
    /* Reverse *.pyc and *.pyo */
    struct st_mlab_searchorder tmp;
    tmp = mlab_searchorder[0];
    mlab_searchorder[0] = mlab_searchorder[1];
    mlab_searchorder[1] = tmp;
    tmp = mlab_searchorder[3];
    mlab_searchorder[3] = mlab_searchorder[4];
    mlab_searchorder[4] = tmp;
  }

  mod = Py_InitModule4("PythonQtImport", NULL, mlabimport_doc,
           NULL, PYTHON_API_VERSION);

  PythonQtImportError = PyErr_NewException("PythonQtImport.PythonQtImportError",
              PyExc_ImportError, NULL);
  if (PythonQtImportError == NULL)
    return;

  Py_INCREF(PythonQtImportError);
  if (PyModule_AddObject(mod, "PythonQtImportError",
             PythonQtImportError) < 0)
    return;

  Py_INCREF(&PythonQtImporter_Type);
  if (PyModule_AddObject(mod, "PythonQtImporter",
             (PyObject *)&PythonQtImporter_Type) < 0)
    return;

  // set our importer into the path_hooks to handle all path on sys.path
  PyObject* classobj = PyDict_GetItemString(PyModule_GetDict(mod), "PythonQtImporter");
  PyObject* path_hooks = PySys_GetObject("path_hooks");
  PyList_Append(path_hooks, classobj);
}
