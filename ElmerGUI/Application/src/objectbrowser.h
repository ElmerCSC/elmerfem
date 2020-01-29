/*****************************************************************************
 *                                                                           *
 *  Elmer, A Finite Element Software for Multiphysical Problems              *
 *                                                                           *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland   *
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
 *  ElmerGUI objectbrowser.h                                                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Author: Saeki Takayuki                                                   *
 *  Original Date: 11 Jan 2020                                               *
 *                                                                           *
 *****************************************************************************/

#ifndef OBJECTBROWSER_H
#define OBJECTBROWSER_H

#include <QWidget>
#include "projectio.h"
#include "dynamiceditor.h"
#include "glwidget.h"
#include "boundarypropertyeditor.h"


class ObjectBrowser : public QDockWidget
{
  Q_OBJECT

//signals:

public:
  ObjectBrowser(QMainWindow *parent, Qt::WindowFlags flags=0);
  ~ObjectBrowser();

//public slots:

private slots:
  void modelSetupSlot();
  void addEquationSlot();
  void addMaterialSlot();
  void addBodyForceSlot();
  void addInitialConditionSlot();
  void addBoundaryConditionSlot();
  void equationSelectedSlot(QAction*);
  void materialSelectedSlot(QAction*); 
  void bodyForceSelectedSlot(QAction*);
  void initialConditionSelectedSlot(QAction*);
  void boundaryConditionSelectedSlot(QAction*);  
  void equationEditorFinishedSlot(int signal, int id);
  void materialEditorFinishedSlot(int signal, int id);
  void bodyForceEditorFinishedSlot(int signal, int id);
  void initialConditionEditorFinishedSlot(int signal, int id);
  void boundaryConditionEditorFinishedSlot(int signal, int id);
  void treeItemClickedSlot(QTreeWidgetItem *item, int column);
  void treeItemDoubleClickedSlot(QTreeWidgetItem *item, int column);
  void treeItemExpandedSlot(QTreeWidgetItem* item);
  void treeItemSelectionChangedSlot();
   
  void openSlot();
  void loadSlot();
  void loadProjectSlot();
  void saveProjectSlot();
  void modelClearSlot();

  void viewFullScreenSlot();
  void viewNormalModeSlot();  
 
  void boundaryDividedSlot(double);
  void boundaryUnifiedSlot();
 
  void boundarySelectedSlot(list_t*);
  
  void boundaryComboChanged(BoundaryPropertyEditor *,QString);
  void bodyComboChanged(BodyPropertyEditor *,QString);  
  void bodyPropertyEditorAccepted(bool); 
  void bodyPropertyEditorDiscarded(bool);
  void boundaryPropertyEditorAccepted(bool);
  void boundaryPropertyEditorDiscarded(bool);
  void bodyCheckBoxChangedSlot(int);
  void boundaryCheckBoxChangedSlot(int);
  
  void focusChangedSlot(QWidget*, QWidget*);

  void meshingStartedSlot();
  void meshingTerminatedSlot();
  void meshingFinishedSlot();
      
private:
  QMainWindow *mainWindow;
  QTreeWidget *tree;

  QTreeWidgetItem *geometryTopLevelTreeItem;
  QTreeWidgetItem *modelTopLevelTreeItem;
  QTreeWidgetItem *bodyPropertyParentTreeItem;
  QTreeWidgetItem *boundaryPropertyParentTreeItem;  
  QTreeWidgetItem *geometryParentTreeItem;
  QTreeWidgetItem *setupParentTreeItem;
  QTreeWidgetItem *equationParentTreeItem;
  QTreeWidgetItem *materialParentTreeItem;
  QTreeWidgetItem *bodyForceParentTreeItem;
  QTreeWidgetItem *initialConditionParentTreeItem;
  QTreeWidgetItem *boundaryConditionParentTreeItem;  
  void addToTree(DynamicEditor*, bool select = false);
  void updateBoundaryProperties(BoundaryPropertyEditor* selectThis = NULL);
  void updateBodyProperties(BodyPropertyEditor* selectThis = NULL);
  void updateEquation();
  void updateMaterial();
  void updateBodyForce();
  void updateInitialCondition();
  void updateBoundaryCondition();
  void snap(QWidget*);
  void selectTreeItemByID(QTreeWidgetItem* parent, int id);
  
  int boundaryListToBodyIndex(list_t*);
  list_t* selectBoundary(BoundaryPropertyEditor *pe, bool append = false, bool select = true);
  list_t* selectBody(BodyPropertyEditor *pe, bool append = false, bool select = true); 
  list_t* boundaryList(int index);
  list_t* bodyList(int index);
  
  bool connect1(const QObject * sender, const char * signal, const QObject * receiver, const char * method, Qt::ConnectionType type = Qt::AutoConnection);
};

#endif // OBJECTBROWSER_H
