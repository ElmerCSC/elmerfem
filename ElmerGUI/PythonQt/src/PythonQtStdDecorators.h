#ifndef _PYTHONQTSTDDECORATORS_H
#define _PYTHONQTSTDDECORATORS_H

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
// \file    PythonQtStdDecorators.h
// \author  Florian Link
// \author  Last changed by $Author: florian $
// \date    2007-04
*/
//----------------------------------------------------------------------------------

#include "PythonQtSystem.h"
#include <Python.h>
#include <QObject>
#include <QVariantList>
#include <QTextDocument>
#include <QColor>

class PYTHONQT_EXPORT PythonQtStdDecorators : public QObject
{
  Q_OBJECT

public slots:
  // additional constructors
  QVariant new_QSize(const QSize& o) { QSize a = o; return a; }
  QVariant new_QSizeF(const QSizeF& o) { QSizeF a = o; return a; }
  QVariant new_QPoint(const QPoint& o) { QPoint a = o; return a; }
  QVariant new_QPointF(const QPointF& o) { QPointF a = o; return a; }
  QVariant new_QRect(const QRect& o) { QRect a = o; return a; }
  QVariant new_QRectF(const QRectF& o) { QRectF a = o; return a; }

  bool connect(QObject* sender, const QByteArray& signal, PyObject* callable);
  bool connect(QObject* sender, const QByteArray& signal, QObject* receiver, const QByteArray& slot);
  bool disconnect(QObject* sender, const QByteArray& signal, PyObject* callable);
  bool disconnect(QObject* sender, const QByteArray& signal, QObject* receiver, const QByteArray& slot);

  QObject* parent(QObject* o);
  void setParent(QObject* o, QObject* parent);

  QVariantList children(QObject* o);

  QString static_Qt_escape(const QString& s) { return Qt::escape(s); }

  //TODO: add findChild/findChildren/children/...
};


#endif
