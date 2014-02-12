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
 *  ElmerGUI edfeditor                                                       *
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

#ifndef EDFEDITOR_H
#define EDFEDITOR_H

#include <QWidget>
#include <QDomDocument>
#include <QIcon>
#include <QTreeWidget>
#include <QHash>
#include "dynamiceditor.h"

class QPushButton;

class EdfEditor : public QWidget
{
  Q_OBJECT

public:
  EdfEditor(QWidget *parent = 0);
  ~EdfEditor();

  QSize minimumSizeHint() const;
  QSize sizeHint() const;

  void setupEditor(QDomDocument *elmerDefs);

signals:

protected:
  void keyPressEvent(QKeyEvent*);
  void keyReleaseEvent(QKeyEvent*);
  
private slots:
  void addButtonClicked();
  void removeButtonClicked();
  void expandCollapseAllButtonClicked();
  void openButtonClicked();
  void appendButtonClicked();
  void saveAsButtonClicked();
  void applyButtonClicked();
  void previewButtonClicked();

  void treeItemClicked(QTreeWidgetItem*, int);
  void updateElement(QTreeWidgetItem*, int);

private:
  QIcon addIcon;
  QIcon removeIcon;
  QIcon collapseIcon;
  QIcon expandIcon;
  QIcon openIcon;
  QIcon appendIcon;
  QIcon saveAsIcon;
  QIcon applyIcon;
  QIcon previewIcon;

  QTreeWidget *edfTree;
  QPushButton *addButton;
  QPushButton *removeButton;
  QPushButton *expandCollapseAllButton;
  QPushButton *openButton;
  QPushButton *appendButton;
  QPushButton *saveAsButton;
  QPushButton *applyButton;
  QPushButton *previewButton;

  DynamicEditor *dynamicEditorSimulation;
  DynamicEditor *dynamicEditorConstants;
  DynamicEditor *dynamicEditorMaterial;
  DynamicEditor *dynamicEditorSolver;
  DynamicEditor *dynamicEditorBC;
  DynamicEditor *dynamicEditorIC;
  DynamicEditor *dynamicEditorEquation;
  DynamicEditor *dynamicEditorBodyForce;

  QDomDocument *elmerDefs;
  QHash<QTreeWidgetItem*, QDomElement> elementForItem;
  void insertItemForElement(QDomElement, QTreeWidgetItem*);

  bool expandCollapseAll;
  QTreeWidgetItem *lastActiveItem;
  bool ctrlPressed;
  bool altPressed;
};

#endif // EDFEDITOR_H
