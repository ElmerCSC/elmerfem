/*****************************************************************************
 *                                                                           *
 *  Elmer, A Finite Element Software for Multiphysical Problems              *
 *                                                                           *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland    *
 *                                                                           *
 *  This program is free software; you can redistribute it and/or            *
 *  modify it under the terms of the GNU General Public License              *
 *  as published by the Free Software Foundation; either version 2           *
 *  of the License, or (at your option) any later version.                   *
 *                                                                           *
 *  This program is distributed in the hope that it will be useful,          *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU General Public License for more details.                             *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License        *
 *  along with this program (in file fem/GPL-2); if not, write to the        *
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,         *
 *  Boston, MA 02110-1301, USA.                                              *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 *  ElmerGUI booundarypropertyeditor                                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Authors: Mikko Lyly, Juha Ruokolainen and Peter Råback                   *
 *  Email:   Juha.Ruokolainen@csc.fi                                         *
 *  Web:     http://www.csc.fi/elmer                                         *
 *  Address: CSC - IT Center for Science Ltd.                                 *
 *           Keilaranta 14                                                   *
 *           02101 Espoo, Finland                                            *
 *                                                                           *
 *  Original Date: 15 Mar 2008                                               *
 *                                                                           *
 *****************************************************************************/

#ifndef BOUNDARYPROPERTYEDITOR_H
#define BOUNDARYPROPERTYEDITOR_H

#include <QWidget>
#include <QDomDocument>
#include "projectio.h"
#include "dynamiceditor.h"
#include "bodypropertyeditor.h"
#include "ui_boundarypropertyeditor.h"

class BoundaryPropertyEditor : public QDialog
{
  Q_OBJECT

signals:
  void BoundaryAsABodyChanged(BoundaryPropertyEditor *, int);
  void BoundaryComboChanged(BoundaryPropertyEditor *, QString);

public:
  BoundaryPropertyEditor(QWidget *parent = 0);
  ~BoundaryPropertyEditor();

  void appendToProject(QDomDocument*, QDomElement*);
  void readFromProject(QDomDocument*, QDomElement*);

  Ui::boundaryPropertyDialog ui;
  DynamicEditor *condition;

  int bodyID;
  BodyPropertyEditor *bodyProperties;

  bool touched;

public slots:
  void boundaryAsABodyChanged(int);
  void boundaryComboChanged(QString);

private slots:
  void applySlot();
  void discardSlot();

private:
  ProjectIO projectIO;

};

#endif // BOUNDARYPROPERTYEDITOR_H
