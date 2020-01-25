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
 *  ElmerGUI objectbrowser.cpp                                               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Author: Saeki Takayuki                                                   *
 *  Original Date: 11 Jan 2020                                               *
 *                                                                           *
 *****************************************************************************/

#include <QtGui>
#include <iostream>
#include <QSettings>
#include "objectbrowser.h"
#include "mainwindow.h"

using namespace std;

ObjectBrowser::ObjectBrowser(QMainWindow *parent, Qt::WindowFlags flags) : QDockWidget("Object Browser", parent, flags)
{
  setTitleBarWidget(new QWidget());
  setFeatures(QDockWidget::NoDockWidgetFeatures /*QDockWidget::DockWidgetFloatable | QDockWidget::DockWidgetMovable*/);
  setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea );
  if(parent) parent->addDockWidget(Qt::LeftDockWidgetArea, this); 
  
  tree = new QTreeWidget(this); 
  tree->setColumnCount(2);
  QStringList headerLabels;
  headerLabels << "Object" << "Value";
  tree->setHeaderLabels(headerLabels);

  geometryTopLevelTreeItem = new QTreeWidgetItem();
  modelTopLevelTreeItem = new QTreeWidgetItem();
  geometryTopLevelTreeItem->setText(0, "Geometry");
  modelTopLevelTreeItem->setText(0, "Model"); 
  
 
  bodyPropertyParentTreeItem = new QTreeWidgetItem();
  boundaryPropertyParentTreeItem = new QTreeWidgetItem();
  geometryParentTreeItem = new QTreeWidgetItem();
  setupParentTreeItem = new QTreeWidgetItem();
  equationParentTreeItem = new QTreeWidgetItem();
  materialParentTreeItem = new QTreeWidgetItem();
  bodyForceParentTreeItem = new QTreeWidgetItem();
  initialConditionParentTreeItem = new QTreeWidgetItem();
  boundaryConditionParentTreeItem = new QTreeWidgetItem();

  bodyPropertyParentTreeItem->setText(0, "Body");
  boundaryPropertyParentTreeItem->setText(0, "Boundary");
  geometryParentTreeItem->setText(0, "Input file");
  setupParentTreeItem->setText(0, "Setup");
  equationParentTreeItem->setText(0, "Equation");
  equationParentTreeItem->setText(1, "[Add...]"); 
  materialParentTreeItem->setText(0, "Material");
  materialParentTreeItem->setText(1, "[Add...]");
  bodyForceParentTreeItem->setText(0, "Body force");
  bodyForceParentTreeItem->setText(1, "[Add...]");
  initialConditionParentTreeItem->setText(0, "Initial condition");
  initialConditionParentTreeItem->setText(1, "[Add...]");
  boundaryConditionParentTreeItem->setText(0, "Boundary condition");
  boundaryConditionParentTreeItem->setText(1, "[Add...]");

  tree->addTopLevelItem(geometryTopLevelTreeItem);
  tree->addTopLevelItem(modelTopLevelTreeItem);

  
  geometryTopLevelTreeItem->addChild(geometryParentTreeItem);  
  geometryTopLevelTreeItem->addChild(bodyPropertyParentTreeItem);    
  geometryTopLevelTreeItem->addChild(boundaryPropertyParentTreeItem);    
  //modelTopLevelTreeItem->addChild(setupParentTreeItem); Setup will not be used so frequently...
  modelTopLevelTreeItem->addChild(equationParentTreeItem);
  modelTopLevelTreeItem->addChild(materialParentTreeItem);
  modelTopLevelTreeItem->addChild(bodyForceParentTreeItem);
  modelTopLevelTreeItem->addChild(initialConditionParentTreeItem);
  modelTopLevelTreeItem->addChild(boundaryConditionParentTreeItem);  
  
  geometryTopLevelTreeItem->setExpanded(true);
  modelTopLevelTreeItem->setExpanded(true);  
  tree->resizeColumnToContents(0);
  tree->resizeColumnToContents(1); 

  setWidget(tree); 
  
  mainWindow = parent;
  MainWindow* mainwindow = (MainWindow*) mainWindow;
  
  connect(mainwindow->modelSetupAct, SIGNAL(triggered()), this, SLOT(modelSetupSlot()));
  connect(mainwindow->addEquationAct, SIGNAL(triggered()), this, SLOT(addEquationSlot()));
  connect(mainwindow->addMaterialAct, SIGNAL(triggered()), this, SLOT(addMaterialSlot()));
  connect(mainwindow->addBodyForceAct, SIGNAL(triggered()), this, SLOT(addBodyForceSlot()));
  connect(mainwindow->addInitialConditionAct, SIGNAL(triggered()), this, SLOT(addInitialConditionSlot()));
  connect(mainwindow->addBoundaryConditionAct, SIGNAL(triggered()), this, SLOT(addBoundaryConditionSlot()));
  connect(mainwindow->equationMenu, SIGNAL(triggered(QAction*)), this, SLOT(equationSelectedSlot(QAction*)));
  connect(mainwindow->materialMenu, SIGNAL(triggered(QAction*)), this, SLOT(materialSelectedSlot(QAction*)));
  connect(mainwindow->bodyForceMenu, SIGNAL(triggered(QAction*)), this, SLOT(bodyForceSelectedSlot(QAction*)));
  connect(mainwindow->initialConditionMenu, SIGNAL(triggered(QAction*)), this, SLOT(initialConditionSelectedSlot(QAction*)));
  connect(mainwindow->boundaryConditionMenu, SIGNAL(triggered(QAction*)), this, SLOT(boundaryConditionSelectedSlot(QAction*)));    
  connect(tree, SIGNAL(itemClicked(QTreeWidgetItem*, int)), this, SLOT(treeItemClickedSlot(QTreeWidgetItem*, int)));     
  connect(tree, SIGNAL(itemDoubleClicked(QTreeWidgetItem*, int)), this, SLOT(treeItemDoubleClickedSlot(QTreeWidgetItem*, int)));     
  connect(tree, SIGNAL(itemExpanded(QTreeWidgetItem*)), this, SLOT(treeItemExpandedSlot(QTreeWidgetItem*)));   //used to update body and boundary parent
  connect(tree, SIGNAL(itemSelectionChanged()), this, SLOT(treeItemSelectionChangedSlot()));
  connect(mainwindow->openAct, SIGNAL(triggered()), this, SLOT(openSlot()));
  connect(mainwindow->loadAct, SIGNAL(triggered()), this, SLOT(loadSlot()));
  connect(mainwindow->loadProjectAct, SIGNAL(triggered()), this, SLOT(loadProjectSlot()));  
  connect(mainwindow->recentProject0Act, SIGNAL(triggered()), this, SLOT(loadProjectSlot()));  
  connect(mainwindow->recentProject1Act, SIGNAL(triggered()), this, SLOT(loadProjectSlot()));  
  connect(mainwindow->recentProject2Act, SIGNAL(triggered()), this, SLOT(loadProjectSlot()));  
  connect(mainwindow->recentProject3Act, SIGNAL(triggered()), this, SLOT(loadProjectSlot()));  
  connect(mainwindow->recentProject4Act, SIGNAL(triggered()), this, SLOT(loadProjectSlot()));
  connect(mainwindow->saveProjectAct, SIGNAL(triggered()), this, SLOT(saveProjectSlot()));
  connect(mainwindow->modelClearAct, SIGNAL(triggered()), this, SLOT(modelClearSlot()));
  
  connect(mainwindow->viewFullScreenAct, SIGNAL(triggered()), this, SLOT(viewFullScreenSlot()));
  connect(mainwindow->glWidget, SIGNAL(escPressed()), this, SLOT(viewNormalModeSlot()));  
 
  connect(mainwindow->boundaryDivide, SIGNAL(signalDoDivideSurface(double)), this, SLOT(boundaryDividedSlot(double)));  
  connect(mainwindow->boundaryDivide, SIGNAL(signalDoDivideEdge(double)), this, SLOT(boundaryDividedSlot(double)));  
  connect(mainwindow->surfaceUnifyAct, SIGNAL(triggered()), this, SLOT(boundaryUnifiedSlot()));
  connect(mainwindow->edgeUnifyAct, SIGNAL(triggered()), this, SLOT(boundaryUnifiedSlot()));
  
  connect(mainwindow->glWidget, SIGNAL(signalBoundarySelected(list_t*)), this, SLOT(boundarySelectedSlot(list_t*)));
  
  connect(mainwindow->meshingThread, SIGNAL(started()), this, SLOT(meshingStartedSlot()));
  connect(mainwindow->meshingThread, SIGNAL(finished()), this, SLOT(meshingFinishedSlot()));
  connect(mainwindow->meshingThread, SIGNAL(terminated()), this, SLOT(meshingTerminatedSlot())); 
    
  QApplication* q = (QApplication*) QCoreApplication::instance();
  connect(q, SIGNAL(focusChanged(QWidget*, QWidget*)), this, SLOT(focusChangedSlot(QWidget*, QWidget*)));
  
  loadProjectSlot();
}

ObjectBrowser::~ObjectBrowser(){
}

void ObjectBrowser::focusChangedSlot(QWidget* old, QWidget* now){

  if(now == NULL) return;

  QWidget* p = now->parentWidget();
  while(p != NULL){
    now = p;
    p = now->parentWidget();
  }
 
  MainWindow* mainwindow = (MainWindow*) mainWindow;
  
  if(now->windowTitle() == "Equation"){
    selectTreeItemByID( equationParentTreeItem, ((DynamicEditor*)now)->ID);
    return;
  }
  if(now->windowTitle() == "Material"){
    selectTreeItemByID( materialParentTreeItem, ((DynamicEditor*)now)->ID);
    return;
  }
  if(now->windowTitle() == "BodyForce"){
    selectTreeItemByID( bodyForceParentTreeItem, ((DynamicEditor*)now)->ID);
    return;
  }
  if(now->windowTitle() == "InitialCondition"){
    selectTreeItemByID( initialConditionParentTreeItem, ((DynamicEditor*)now)->ID);
    return;
  }
  if(now->windowTitle() == "BoundaryCondition"){
    selectTreeItemByID( boundaryConditionParentTreeItem, ((DynamicEditor*)now)->ID);
    return;
  }
  if(now->windowTitle().indexOf(QString("Properties for body")) == 0){
    int n = bodyPropertyParentTreeItem->childCount();
    for(int i = 0; i < n; i++){
      QTreeWidgetItem* item = bodyPropertyParentTreeItem->child(i);
      if(item->data(0, Qt::UserRole).value<qulonglong>() == (qulonglong)now){
        tree->setCurrentItem(item);
        return;
      }
    }
    return;
  }
  if(now->windowTitle().indexOf(QString("Properties for boundary")) == 0){
    int n = boundaryPropertyParentTreeItem->childCount();
    for(int i = 0; i < n; i++){
      QTreeWidgetItem* item = boundaryPropertyParentTreeItem->child(i);
      if(item->data(0, Qt::UserRole).value<qulonglong>() == (qulonglong)now){
        tree->setCurrentItem(item);
        return;
      }
    }
    return;
  }
}

void ObjectBrowser::treeItemClickedSlot(QTreeWidgetItem *item, int column){
/*  Moved to treeItemDoubleClickedSlot(QTreeWidgetItem *item, int column)
  MainWindow* mainwindow = (MainWindow*) mainWindow;
  if(item == geometryParentTreeItem) return;

  if(item == setupParentTreeItem){
    mainwindow->modelSetupSlot();
    return;
  }
  if(item == equationParentTreeItem){
	if(column == 1){
		mainwindow->addEquationSlot();
		addEquationSlot();
	}
	return;
  }
  if(item == materialParentTreeItem){
	if(column == 1){
		mainwindow->addMaterialSlot();
		addMaterialSlot();
	}
	return;
  }
  if(item == bodyForceParentTreeItem){
	if(column == 1){
		mainwindow->addBodyForceSlot();
		addBodyForceSlot();
	}
	return;
  }
  if(item == initialConditionParentTreeItem){
	if(column == 1){
		mainwindow->addInitialConditionSlot();
		addInitialConditionSlot();
	}
	return;
  }
  if(item == boundaryConditionParentTreeItem){
	if(column == 1){
		mainwindow->addBoundaryConditionSlot();
		addBoundaryConditionSlot();
	}
	return;
  } 
*/
}

void ObjectBrowser::treeItemDoubleClickedSlot(QTreeWidgetItem *item, int column){

  MainWindow* mainwindow = (MainWindow*) mainWindow;
  
    // Show selected dynamic editor  
  if(item->parent() == equationParentTreeItem){
    DynamicEditor *pe = (DynamicEditor*) item->data(0, Qt::UserRole).value<qulonglong>();
    //snap(pe); 
    mainwindow->equationSelectedSlot(pe->menuAction);

    QGroupBox *box = (QGroupBox*) pe->spareScroll->widget();
    if(!mainwindow->glWidget->hasMesh() || box == NULL) return;
    for(int i=2; i< box->children().size(); i++){
      QCheckBox* cb = (QCheckBox*) box->children()[i];
      connect1(cb, SIGNAL(stateChanged(int)), this, SLOT(bodyCheckBoxChangedSlot(int)));
    }

    return;
  } 
  if(item->parent() == materialParentTreeItem){
  	DynamicEditor *pe = (DynamicEditor*) item->data(0, Qt::UserRole).value<qulonglong>();
    //snap(pe);
    mainwindow->materialSelectedSlot(pe->menuAction);

    QGroupBox *box = (QGroupBox*) pe->spareScroll->widget();
    if(!mainwindow->glWidget->hasMesh() || box == NULL) return;
    for(int i=2; i< box->children().size(); i++){
      QCheckBox* cb = (QCheckBox*) box->children()[i];
      connect1(cb, SIGNAL(stateChanged(int)), this, SLOT(bodyCheckBoxChangedSlot(int)));
    }

	return;
  }
  if(item->parent() == bodyForceParentTreeItem){
  	DynamicEditor *pe = (DynamicEditor*) item->data(0, Qt::UserRole).value<qulonglong>();  
    //snap(pe);
    mainwindow->bodyForceSelectedSlot(pe->menuAction);

    QGroupBox *box = (QGroupBox*) pe->spareScroll->widget();
    if(!mainwindow->glWidget->hasMesh() || box == NULL) return;
    for(int i=2; i< box->children().size(); i++){
      QCheckBox* cb = (QCheckBox*) box->children()[i];
      connect1(cb, SIGNAL(stateChanged(int)), this, SLOT(bodyCheckBoxChangedSlot(int)));
    }

    return;
  }
  if(item->parent() == initialConditionParentTreeItem){
  	DynamicEditor *pe = (DynamicEditor*) item->data(0, Qt::UserRole).value<qulonglong>();   
    //snap(pe);
    mainwindow->initialConditionSelectedSlot(pe->menuAction);

    QGroupBox *box = (QGroupBox*) pe->spareScroll->widget();
    if(!mainwindow->glWidget->hasMesh() || box == NULL) return;
    for(int i=2; i< box->children().size(); i++){
      QCheckBox* cb = (QCheckBox*) box->children()[i];
      connect1(cb, SIGNAL(stateChanged(int)), this, SLOT(bodyCheckBoxChangedSlot(int)));
    }

    return;
  }
  if(item->parent() == boundaryConditionParentTreeItem){
  	DynamicEditor *pe = (DynamicEditor*) item->data(0, Qt::UserRole).value<qulonglong>();   
    //snap(pe);
    mainwindow->boundaryConditionSelectedSlot(pe->menuAction);

    //connect checkbox
    QGroupBox *box = (QGroupBox*) pe->spareScroll->widget();
    if(!mainwindow->glWidget->hasMesh() || box == NULL) return;
    for(int i=2; i< box->children().size(); i++){
      QCheckBox* cb = (QCheckBox*) box->children()[i];
      connect1(cb, SIGNAL(stateChanged(int)), this, SLOT(boundaryCheckBoxChangedSlot(int)));
    }

    return;
  }
  
  
  //show property editor
  //bool bodyEditActive = mainwindow->bodyEditActive;
  //bool bcEditActive = mainwindow->bcEditActive;
  if(item->parent() == boundaryPropertyParentTreeItem ){
    BoundaryPropertyEditor* pe = (BoundaryPropertyEditor*) item->data(0, Qt::UserRole).value<qulonglong>();
    if(pe == NULL ){ cout << " BoundaryPropertyEditor NULL" << endl; return;}
    list_t* l = selectBoundary(pe);
    mainwindow->bcEditActive = true;
    mainwindow->bodyEditActive = false;
    mainwindow->glWidget->bodyEditActive = false;
    mainwindow->bcEditAct->setChecked(true);
    mainwindow->bodyEditAct->setChecked(false);
    mainwindow->boundarySelectedSlot(l);
    //mainwindow->bcEditActive = bcEditActive;
    //mainwindow->bodyEditActive = bodyEditActive;
		connect1( pe,SIGNAL(BoundaryComboChanged(BoundaryPropertyEditor *,QString)), this, SLOT(boundaryComboChanged(BoundaryPropertyEditor *,QString)) );
		connect1( pe->ui.applyButton,SIGNAL(clicked(bool)), this, SLOT(boundaryPropertyEditorAccepted(bool)));
		connect1( pe->ui.discardButton,SIGNAL(clicked(bool)), this, SLOT(boundaryPropertyEditorDiscarded(bool)));	    
    snap(pe); pe->show(); pe->raise();
    return;
  }
  if(item->parent() != NULL && item->parent()->parent() == boundaryPropertyParentTreeItem){
    BoundaryPropertyEditor* pe = (BoundaryPropertyEditor*) item->parent()->data(0, Qt::UserRole).value<qulonglong>();
    if(pe == NULL){ cout << " BoundaryPropertyEditor NULL" << endl; return;}
    list_t* l = selectBoundary(pe);
    mainwindow->bcEditActive = true;	
    mainwindow->bodyEditActive = false;
    mainwindow->glWidget->bodyEditActive = false;
    mainwindow->bcEditAct->setChecked(true);
    mainwindow->bodyEditAct->setChecked(false);
    mainwindow->boundarySelectedSlot(l);
    //mainwindow->bcEditActive = bcEditActive;
    //mainwindow->bodyEditActive = bodyEditActive;
		connect1( pe,SIGNAL(BoundaryComboChanged(BoundaryPropertyEditor *,QString)), this, SLOT(boundaryComboChanged(BoundaryPropertyEditor *,QString)) );
		connect1( pe->ui.applyButton,SIGNAL(clicked(bool)), this, SLOT(boundaryPropertyEditorAccepted(bool)));
		connect1( pe->ui.discardButton,SIGNAL(clicked(bool)), this, SLOT(boundaryPropertyEditorDiscarded(bool)));
	  snap(pe); pe->show(); pe->raise();
    return;
  }
  if(item->parent() == bodyPropertyParentTreeItem){
    BodyPropertyEditor* pe = (BodyPropertyEditor*) item->data(0, Qt::UserRole).value<qulonglong>();
    if(pe == NULL){ cout << " BodyPropertyEditor NULL" << endl; return;}
    list_t* l = selectBody(pe);
    mainwindow->bcEditActive = false;
    mainwindow->bodyEditActive = true;
    mainwindow->glWidget->bodyEditActive = true;
    mainwindow->bcEditAct->setChecked(false);
    mainwindow->bodyEditAct->setChecked(true);
    mainwindow->boundarySelectedSlot(l);
    //mainwindow->bcEditActive = bcEditActive;
    //mainwindow->bodyEditActive = bodyEditActive;
    connect1( pe,SIGNAL(BodyMaterialComboChanged(BodyPropertyEditor *,QString)), this, SLOT(bodyComboChanged(BodyPropertyEditor *,QString)) );
    connect1( pe,SIGNAL(BodyInitialComboChanged(BodyPropertyEditor *,QString)), this, SLOT(bodyComboChanged(BodyPropertyEditor *,QString)) );
    connect1( pe,SIGNAL(BodyForceComboChanged(BodyPropertyEditor *,QString)), this, SLOT(bodyComboChanged(BodyPropertyEditor *,QString)) );
    connect1( pe,SIGNAL(BodyEquationComboChanged(BodyPropertyEditor *,QString)), this, SLOT(bodyComboChanged(BodyPropertyEditor *,QString)) );
    connect1( pe->ui.applyButton,SIGNAL(clicked(bool)), this, SLOT(bodyPropertyEditorAccepted(bool)));
    connect1( pe->ui.discardButton,SIGNAL(clicked(bool)), this, SLOT(bodyPropertyEditorDiscarded(bool)));	
    snap(pe); pe->show(); pe->raise();
    tree->collapseItem(item);
    return;
  }
  if(item->parent() != NULL && item->parent()->parent() == bodyPropertyParentTreeItem){
    BodyPropertyEditor* pe = (BodyPropertyEditor*) item->parent()->data(0, Qt::UserRole).value<qulonglong>();
    if(pe == NULL){ cout << " BodyPropertyEditor NULL" << endl; return;}
    list_t* l = selectBody(pe);
    mainwindow->bcEditActive = false;
    mainwindow->bodyEditActive = true;
    mainwindow->glWidget->bodyEditActive = true;
    mainwindow->bcEditAct->setChecked(false);
    mainwindow->bodyEditAct->setChecked(true);	
    mainwindow->boundarySelectedSlot(l);
    //mainwindow->bcEditActive = bcEditActive;
    //mainwindow->bodyEditActive = bodyEditActive;
    connect1( pe,SIGNAL(BodyMaterialComboChanged(BodyPropertyEditor *,QString)), this, SLOT(bodyComboChanged(BodyPropertyEditor *,QString)) );
    connect1( pe,SIGNAL(BodyInitialComboChanged(BodyPropertyEditor *,QString)), this, SLOT(bodyComboChanged(BodyPropertyEditor *,QString)) );
    connect1( pe,SIGNAL(BodyForceComboChanged(BodyPropertyEditor *,QString)), this, SLOT(bodyComboChanged(BodyPropertyEditor *,QString)) );
    connect1( pe,SIGNAL(BodyEquationComboChanged(BodyPropertyEditor *,QString)), this, SLOT(bodyComboChanged(BodyPropertyEditor *,QString)) );
    connect1( pe->ui.applyButton,SIGNAL(clicked(bool)), this, SLOT(bodyPropertyEditorAccepted(bool)));
    connect1( pe->ui.discardButton,SIGNAL(clicked(bool)), this, SLOT(bodyPropertyEditorDiscarded(bool)));	
    snap(pe); pe->show(); pe->raise();
    return;
  }
  
  // double clicking [Add...] button
  if(item == equationParentTreeItem){
	if(column == 1){
		mainwindow->addEquationSlot();
		addEquationSlot();
		tree->collapseItem(equationParentTreeItem);
	}
	return;
  }
  if(item == materialParentTreeItem){
	if(column == 1){
		mainwindow->addMaterialSlot();
		addMaterialSlot();	
		tree->collapseItem(materialParentTreeItem);
	}
	return;
  }
  if(item == bodyForceParentTreeItem){
	if(column == 1){
		mainwindow->addBodyForceSlot();
		addBodyForceSlot();
		tree->collapseItem(bodyForceParentTreeItem);
	}
	return;
  }
  if(item == initialConditionParentTreeItem){
	if(column == 1){
		mainwindow->addInitialConditionSlot();
		addInitialConditionSlot();
		tree->collapseItem(initialConditionParentTreeItem);
	}
	return;
  }
  if(item == boundaryConditionParentTreeItem){
	if(column == 1){
		mainwindow->addBoundaryConditionSlot();
		addBoundaryConditionSlot();
		tree->collapseItem(boundaryConditionParentTreeItem);
	}
	return;
  } 

}

void ObjectBrowser::addToTree(DynamicEditor *de, bool select /*= false*/){
	if(de == NULL) return;
	if(de->ID < 0) return; //removed item
  QTreeWidgetItem *treeItem = new QTreeWidgetItem();
  treeItem->setText(0, de->nameEdit->text().trimmed());// +" ["  + QString::number(de->ID) + "]"); 
  QString title = de->windowTitle();
  if(title == "Equation"){
    equationParentTreeItem->addChild(treeItem);
    connect1(de, SIGNAL(dynamicEditorReady(int,int)), this, SLOT(equationEditorFinishedSlot(int,int)));
  
    QGroupBox *box = (QGroupBox*) de->spareScroll->widget();
    bool used = false;
    if(box != NULL){
      for(int i=2; i< box->children().size(); i++){
        QCheckBox* cb = (QCheckBox*) box->children()[i];
        connect1(cb, SIGNAL(stateChanged(int)), this, SLOT(bodyCheckBoxChangedSlot(int)));
        used |= cb->isChecked(); 
      }
    }
    //if(!used) treeItem->setText(0, "*" + de->nameEdit->text().trimmed());
  }
  else if(title == "Material"){
    materialParentTreeItem->addChild(treeItem);
    connect1(de, SIGNAL(dynamicEditorReady(int,int)), this, SLOT(materialEditorFinishedSlot(int,int)));	
    
    QGroupBox *box = (QGroupBox*) de->spareScroll->widget();
    bool used = false;
    if(box != NULL){
      for(int i=2; i< box->children().size(); i++){
        QCheckBox* cb = (QCheckBox*) box->children()[i];
        connect1(cb, SIGNAL(stateChanged(int)), this, SLOT(bodyCheckBoxChangedSlot(int)));
        used |= cb->isChecked(); 
      }
    }
    //if(!used) treeItem->setText(0, "*" + de->nameEdit->text().trimmed());
  }
  else if(title == "BodyForce"){
    bodyForceParentTreeItem->addChild(treeItem);
    connect1(de, SIGNAL(dynamicEditorReady(int,int)), this, SLOT(bodyForceEditorFinishedSlot(int,int)));

    QGroupBox *box = (QGroupBox*) de->spareScroll->widget();
    bool used = false;
    if(box != NULL){
      for(int i=2; i< box->children().size(); i++){
        QCheckBox* cb = (QCheckBox*) box->children()[i];
        connect1(cb, SIGNAL(stateChanged(int)), this, SLOT(bodyCheckBoxChangedSlot(int)));
        used |= cb->isChecked(); 
      }
    }
    //if(!used) treeItem->setText(0, "*" + de->nameEdit->text().trimmed());
  }
  else if(title == "InitialCondition"){
    initialConditionParentTreeItem->addChild(treeItem);
    connect1(de, SIGNAL(dynamicEditorReady(int,int)), this, SLOT(initialConditionEditorFinishedSlot(int,int)));

    QGroupBox *box = (QGroupBox*) de->spareScroll->widget();
    bool used = false;
    if(box != NULL){
      for(int i=2; i< box->children().size(); i++){
        QCheckBox* cb = (QCheckBox*) box->children()[i];
        connect1(cb, SIGNAL(stateChanged(int)), this, SLOT(bodyCheckBoxChangedSlot(int)));
        used |= cb->isChecked(); 
      }
    }
    //if(!used) treeItem->setText(0, "*" + de->nameEdit->text().trimmed());
  }
  else if(title == "BoundaryCondition"){
    boundaryConditionParentTreeItem->addChild(treeItem);
    connect1(de, SIGNAL(dynamicEditorReady(int,int)), this, SLOT(boundaryConditionEditorFinishedSlot(int,int)));

    QGroupBox *box = (QGroupBox*) de->spareScroll->widget();
    bool used = false;
    if(box != NULL){
      for(int i=2; i< box->children().size(); i++){
        QCheckBox* cb = (QCheckBox*) box->children()[i];
        connect1(cb, SIGNAL(stateChanged(int)), this, SLOT(boundaryCheckBoxChangedSlot(int)));
        used |= cb->isChecked(); 
      }
    }
    //if(!used) treeItem->setText(0, "*" + de->nameEdit->text().trimmed());
  }
  else {
    cout << " ***Warning*** In ObjectBrowser::addTotree(), the DynamicEditor did not match to any of categories."; 	
  }
  
  treeItem->setData(0, Qt::UserRole, QVariant::fromValue( (qulonglong)de ));
  if(select) tree->setCurrentItem(treeItem);
  else treeItem->parent()->setExpanded(true);
  
//  tree->resizeColumnToContents(0);
//  tree->resizeColumnToContents(1);   

}

void ObjectBrowser::modelSetupSlot(){}
void ObjectBrowser::addEquationSlot(){
	MainWindow* mainwindow = (MainWindow*) mainWindow;
	DynamicEditor* de = (DynamicEditor*) mainwindow->equationEditor.last();
	addToTree(de, true);
	snap(de);
}
void ObjectBrowser::addMaterialSlot(){
	MainWindow* mainwindow = (MainWindow*) mainWindow;
	DynamicEditor* de = (DynamicEditor*) mainwindow->materialEditor.last();
	addToTree(de, true);
	snap(de);
}
void ObjectBrowser::addBodyForceSlot(){
	MainWindow* mainwindow = (MainWindow*) mainWindow;
	DynamicEditor* de = (DynamicEditor*) mainwindow->bodyForceEditor.last();
	addToTree(de, true);
	snap(de);
}
void ObjectBrowser::addInitialConditionSlot(){
	MainWindow* mainwindow = (MainWindow*) mainWindow;
	DynamicEditor* de = (DynamicEditor*) mainwindow->initialConditionEditor.last();
	addToTree(de, true);
	snap(de);
}
void ObjectBrowser::addBoundaryConditionSlot(){
	MainWindow* mainwindow = (MainWindow*) mainWindow;
	DynamicEditor* de = (DynamicEditor*) mainwindow->boundaryConditionEditor.last();
	addToTree(de, true);
	snap(de);
}
void ObjectBrowser::equationSelectedSlot(QAction* action){
	int n = equationParentTreeItem->childCount();
	for(int i= 0; i < n; i++){
		QTreeWidgetItem* child = equationParentTreeItem->child(i);
		DynamicEditor* de = (DynamicEditor*) child->data(0, Qt::UserRole).value<qulonglong>();
		if(de->menuAction == action){
			tree->setCurrentItem(child);
			return;
		}
	}
}
void ObjectBrowser::materialSelectedSlot(QAction* action){
	int n = materialParentTreeItem->childCount();
	for(int i= 0; i < n; i++){
		QTreeWidgetItem* child = materialParentTreeItem->child(i);
		DynamicEditor* de = (DynamicEditor*) child->data(0, Qt::UserRole).value<qulonglong>();
		if(de->menuAction == action){
			tree->setCurrentItem(child);
			return;
		}
	}	
}
void ObjectBrowser::bodyForceSelectedSlot(QAction* action){
	int n = bodyForceParentTreeItem->childCount();
	for(int i= 0; i < n; i++){
		QTreeWidgetItem* child = bodyForceParentTreeItem->child(i);
		DynamicEditor* de = (DynamicEditor*) child->data(0, Qt::UserRole).value<qulonglong>();
		if(de->menuAction == action){
			tree->setCurrentItem(child);
			return;
		}
	}	
}
void ObjectBrowser::initialConditionSelectedSlot(QAction* action){
	int n = initialConditionParentTreeItem->childCount();
	for(int i= 0; i < n; i++){
		QTreeWidgetItem* child = initialConditionParentTreeItem->child(i);
		DynamicEditor* de = (DynamicEditor*) child->data(0, Qt::UserRole).value<qulonglong>();
		if(de->menuAction == action){
			tree->setCurrentItem(child);
			return;
		}
	}	
}
void ObjectBrowser::boundaryConditionSelectedSlot(QAction* action){
	int n = boundaryConditionParentTreeItem->childCount();
	for(int i= 0; i < n; i++){
		QTreeWidgetItem* child = boundaryConditionParentTreeItem->child(i);
		DynamicEditor* de = (DynamicEditor*) child->data(0, Qt::UserRole).value<qulonglong>();
		if(de->menuAction == action){
			tree->setCurrentItem(child);
			return;
		}
	}	
}  

void ObjectBrowser::equationEditorFinishedSlot(int signal, int id)
{
	updateEquation();	
	updateBodyProperties();
	updateBoundaryProperties();
  if ( signal != MAT_NEW ) selectTreeItemByID(equationParentTreeItem, id);
  else{
    MainWindow* mainwindow = (MainWindow*) mainWindow;
    DynamicEditor*de = mainwindow->equationEditor.last();
    if(de != NULL) selectTreeItemByID(equationParentTreeItem, de->ID);
    else selectTreeItemByID(equationParentTreeItem, -1);
  }
}

void ObjectBrowser::materialEditorFinishedSlot(int signal, int id )
{
	updateMaterial();
	updateBodyProperties();
	updateBoundaryProperties();
  if ( signal != MAT_NEW ) selectTreeItemByID(materialParentTreeItem, id);
  else{
    MainWindow* mainwindow = (MainWindow*) mainWindow;
    DynamicEditor*de = mainwindow->materialEditor.last();
    if(de != NULL) selectTreeItemByID(materialParentTreeItem, de->ID);
    else selectTreeItemByID(materialParentTreeItem, -1);
  }
}

void ObjectBrowser::bodyForceEditorFinishedSlot(int signal, int id)
{
	updateBodyForce();
	updateBodyProperties();
	updateBoundaryProperties();
  if ( signal != MAT_NEW ) selectTreeItemByID(bodyForceParentTreeItem, id);
  else{
    MainWindow* mainwindow = (MainWindow*) mainWindow;
    DynamicEditor*de = mainwindow->bodyForceEditor.last();
    if(de != NULL) selectTreeItemByID(bodyForceParentTreeItem, de->ID);
    else selectTreeItemByID(bodyForceParentTreeItem, -1);
  }
}

void ObjectBrowser::initialConditionEditorFinishedSlot(int signal, int id)
{
	updateInitialCondition();
	updateBodyProperties();
	updateBoundaryProperties();
  if ( signal != MAT_NEW ) selectTreeItemByID(initialConditionParentTreeItem, id);
  else{
    MainWindow* mainwindow = (MainWindow*) mainWindow;
    DynamicEditor*de = mainwindow->initialConditionEditor.last();
    if(de != NULL) selectTreeItemByID(initialConditionParentTreeItem, de->ID);
    else selectTreeItemByID(initialConditionParentTreeItem, -1);
  }
}

void ObjectBrowser::boundaryConditionEditorFinishedSlot(int signal, int id)
{
	updateBoundaryCondition();
	updateBodyProperties();
	updateBoundaryProperties();
  if ( signal != MAT_NEW ) selectTreeItemByID(boundaryConditionParentTreeItem, id);
  else{
    MainWindow* mainwindow = (MainWindow*) mainWindow;
    DynamicEditor*de = mainwindow->boundaryConditionEditor.last();
    if(de != NULL) selectTreeItemByID(boundaryConditionParentTreeItem, de->ID);
    else selectTreeItemByID(boundaryConditionParentTreeItem, -1);
  }
}


void ObjectBrowser::openSlot()
{
	MainWindow* mainwindow = (MainWindow*) mainWindow;
	QFileInfo fi(mainwindow->geometryInputFileName);
	geometryParentTreeItem->setText(0, "Input file");
	geometryParentTreeItem->setText(1, fi.fileName());
	
	bodyPropertyParentTreeItem->setExpanded(false);
	bodyPropertyParentTreeItem->addChild(new QTreeWidgetItem()); //dummy
	boundaryPropertyParentTreeItem->setExpanded(false);
	boundaryPropertyParentTreeItem->addChild(new QTreeWidgetItem()); //dummy
}

void ObjectBrowser::loadSlot()
{
	MainWindow* mainwindow = (MainWindow*) mainWindow;
	QFileInfo fi(mainwindow->saveDirName);
	geometryParentTreeItem->setText(0, "Elmer mesh dir");
	geometryParentTreeItem->setText(1, fi.fileName());
	
	bodyPropertyParentTreeItem->setExpanded(false);
	bodyPropertyParentTreeItem->addChild(new QTreeWidgetItem()); //dummy
	boundaryPropertyParentTreeItem->setExpanded(false);
	boundaryPropertyParentTreeItem->addChild(new QTreeWidgetItem()); //dummy
}


void ObjectBrowser::loadProjectSlot()
{
  MainWindow* mainwindow = (MainWindow*) mainWindow;
	if(mainwindow->currentProjectDirName.isEmpty()) return;
	
	modelClearSlot();

	int n = mainwindow->equationEditor.size();
	for(int i= 0; i < n; i++){
		addToTree(mainwindow->equationEditor.at(i));
	}
	n = mainwindow->materialEditor.size();
	for(int i= 0; i < n; i++){
		addToTree(mainwindow->materialEditor.at(i));
	}
	n = mainwindow->bodyForceEditor.size();
	for(int i= 0; i < n; i++){
		addToTree(mainwindow->bodyForceEditor.at(i));
	}
	n = mainwindow->initialConditionEditor.size();
	for(int i= 0; i < n; i++){
		addToTree(mainwindow->initialConditionEditor.at(i));
	}
	n = mainwindow->boundaryConditionEditor.size();
	for(int i= 0; i < n; i++){
		addToTree(mainwindow->boundaryConditionEditor.at(i));
	}
	
	if(!mainwindow->geometryInputFileName.isEmpty()){
		QFileInfo fi(mainwindow->geometryInputFileName);
		geometryParentTreeItem->setText(0, "Input file");	
		geometryParentTreeItem->setText(1, fi.fileName());
	}else if(!mainwindow->saveDirName.isEmpty()){
		QFileInfo fi(mainwindow->saveDirName);
		geometryParentTreeItem->setText(0, "Elmer mesh dir");
		geometryParentTreeItem->setText(1, fi.fileName());
	}else{
		geometryParentTreeItem->setText(0, "Input file");	
		geometryParentTreeItem->setText(1, "");	
	}
	
	bodyPropertyParentTreeItem->setExpanded(false);
	bodyPropertyParentTreeItem->addChild(new QTreeWidgetItem()); //dummy
	boundaryPropertyParentTreeItem->setExpanded(false);
	boundaryPropertyParentTreeItem->addChild(new QTreeWidgetItem()); //dummy
	
	tree->resizeColumnToContents(0);
  tree->resizeColumnToContents(1);
}



void ObjectBrowser::saveProjectSlot()
{	
  MainWindow* mainwindow = (MainWindow*) mainWindow;
	if(mainwindow->currentProjectDirName.isEmpty()) return;
	
	if( !mainwindow->geometryInputFileName.isEmpty()){
		QFileInfo fi(mainwindow->geometryInputFileName);
		geometryParentTreeItem->setText(0, "Input file");	
		geometryParentTreeItem->setText(1, fi.fileName());	
	}
	else if( !mainwindow->saveDirName.isEmpty()){
		QFileInfo fi(mainwindow->saveDirName);
		geometryParentTreeItem->setText(0, "Elmer mesh dir");
		geometryParentTreeItem->setText(1, fi.fileName());
	}else{
		geometryParentTreeItem->setText(0, "Input file");	
		geometryParentTreeItem->setText(1, "");	
	}
}

void ObjectBrowser::modelClearSlot(){
  geometryParentTreeItem->setText(1,"");
  QTreeWidgetItem* treeItem = 0;
  tree->setCurrentItem(geometryTopLevelTreeItem);
  while(treeItem = bodyPropertyParentTreeItem->takeChild(0)) delete treeItem;
  while(treeItem = boundaryPropertyParentTreeItem->takeChild(0)) delete treeItem;
  while(treeItem = equationParentTreeItem->takeChild(0)) delete treeItem;
  while(treeItem = materialParentTreeItem->takeChild(0)) delete treeItem;
  while(treeItem = bodyForceParentTreeItem->takeChild(0)) delete treeItem;
  while(treeItem = initialConditionParentTreeItem->takeChild(0)) delete treeItem;
  while(treeItem = boundaryConditionParentTreeItem->takeChild(0)) delete treeItem;

  openSlot();
}

void ObjectBrowser::treeItemExpandedSlot(QTreeWidgetItem* item){ //used to update body and boundary parent

	if(item == bodyPropertyParentTreeItem ){
		updateBodyProperties();	
		for(int i = 0; i < bodyPropertyParentTreeItem->childCount(); i++) bodyPropertyParentTreeItem->child(i)->setExpanded(true);
		updateBoundaryProperties(); //This was added to avoid crushing when chosing body property just after loading 3D mesh without expanding boundary condition parent item.
	}
	
	if(item == boundaryPropertyParentTreeItem ) updateBoundaryProperties();
}

void ObjectBrowser::updateBoundaryProperties(BoundaryPropertyEditor* selectThis /*=NULL*/){

	MainWindow* mainwindow = (MainWindow*) mainWindow;
	if(tree->currentItem() != NULL && tree->currentItem()->parent() == boundaryPropertyParentTreeItem){
    tree->setCurrentItem(boundaryPropertyParentTreeItem);// to improve speed
  }
	QTreeWidgetItem* treeItem = 0;
	while(treeItem = boundaryPropertyParentTreeItem->takeChild(0)) delete treeItem;

	BoundaryPropertyEditor* pe;
	
  for( int i=0; i< mainwindow->glWidget->boundaryMap.count(); i++ )
  {
    int n = mainwindow->glWidget->boundaryMap.key(i);
    if ( n >= 0 ) {
      int m = mainwindow->glWidget->boundaryMap.value(n);

      if(m >= mainwindow->boundaryPropertyEditor.size()) mainwindow->boundaryPropertyEditor.resize(m + 1);

      if(!mainwindow->boundaryPropertyEditor[m])  mainwindow->boundaryPropertyEditor[m] = new BoundaryPropertyEditor;

      BoundaryPropertyEditor *pe = mainwindow->boundaryPropertyEditor[m];		
    
      treeItem = new QTreeWidgetItem();
      boundaryPropertyParentTreeItem->addChild(treeItem);			
      
      treeItem->setText(0, "Boundary " + QString::number(n));
      //if(pe->touched && pe->condition != 0) treeItem->setText(0, "Boundary " + QString::number(n));
      //else treeItem->setText(0, "*Boundary " + QString::number(n));

      treeItem->setData(0, Qt::UserRole, QVariant::fromValue( (qulonglong)pe));		
      if(pe->condition != 0){			
        if(pe->touched) treeItem->setText(1, pe->condition->nameEdit->text().trimmed());
        else treeItem->setText(1, "<Canceled>" + pe->condition->nameEdit->text().trimmed());
        treeItem->setData(1, Qt::UserRole, QVariant::fromValue( (qulonglong)pe->condition));	
      }		
    }
  }	 
  
  if(selectThis != NULL){  
	QTreeWidgetItem* item;
	for(int i = 0; i < boundaryPropertyParentTreeItem->childCount(); i++){
		item = boundaryPropertyParentTreeItem->child(i);
		if( item->data(0, Qt::UserRole).value<qulonglong>() == (qulonglong)selectThis ){
			tree->setCurrentItem(item);
			return;
		}
	}
  }
}

void ObjectBrowser::boundaryUnifiedSlot(){
  updateBodyProperties();  
  updateBoundaryProperties();
}
void ObjectBrowser::boundaryDividedSlot(double d){
	bodyPropertyParentTreeItem->setExpanded(false);
	boundaryPropertyParentTreeItem->setExpanded(false);
}

void ObjectBrowser::boundarySelectedSlot(list_t* l){

	bodyPropertyParentTreeItem->setExpanded(true);
	boundaryPropertyParentTreeItem->setExpanded(true);


	MainWindow* mainwindow = (MainWindow*) mainWindow; 
	QDialog* dialog = NULL;



  if(mainwindow->bodyEditActive  && mainwindow->glWidget->currentlySelectedBody != -1)
  {
    int m = mainwindow->glWidget->bodyMap.value(mainwindow->glWidget->currentlySelectedBody);
    dialog = (QDialog*) mainwindow->bodyPropertyEditor[m];
  }else{
  
    for( int i=0; i< mainwindow->glWidget->boundaryMap.count(); i++ ) {
      int n = mainwindow->glWidget->boundaryMap.key(i);
      if ( n >= 0 ) {
        int m = mainwindow->glWidget->boundaryMap.value(n);     
        if(m >= 0 && m < mainwindow->boundaryPropertyEditor.size() && mainwindow->boundaryPropertyEditor[m] != NULL &&  boundaryList(n) == l){
          dialog = mainwindow->boundaryPropertyEditor[m];	
        }
      }
    }
    
    if(dialog == NULL){
      for( int i=0; i< mainwindow->glWidget->bodyMap.count(); i++ ) {
        int n = mainwindow->glWidget->bodyMap.key(i);
        if ( n >= 0 ) {
          int m = mainwindow->glWidget->bodyMap.value(n);     
          if(m >= 0 && m < mainwindow->bodyPropertyEditor.size() && mainwindow->bodyPropertyEditor[m] != NULL &&  bodyList(n) == l){
            dialog = mainwindow->bodyPropertyEditor[m];	
          }
        }
      }
    }
  }
	
	if(dialog == NULL){
		cout << " could not find selected body/boundary in ObjectBrowser::boundarySelectedSlot(list_t*).  list_t:" << (qulonglong)l<<endl;
		return;
	}

	for(int i = 0; i < boundaryPropertyParentTreeItem->childCount(); i++){
		QTreeWidgetItem* child = boundaryPropertyParentTreeItem->child(i) ;
		if( child->data(0, Qt::UserRole) == (qulonglong) dialog){
			tree->setCurrentItem(child);
			BoundaryPropertyEditor* pe = (BoundaryPropertyEditor*) dialog;
			connect1( pe,SIGNAL(BoundaryComboChanged(BoundaryPropertyEditor *,QString)), this, SLOT(boundaryComboChanged(BoundaryPropertyEditor *,QString)) );
      connect1( pe->ui.applyButton,SIGNAL(clicked(bool)), this, SLOT(boundaryPropertyEditorAccepted(bool)));
      connect1( pe->ui.discardButton,SIGNAL(clicked(bool)), this, SLOT(boundaryPropertyEditorDiscarded(bool)));				
			return;
		}
	}

	for(int i = 0; i < bodyPropertyParentTreeItem->childCount(); i++){
		QTreeWidgetItem* child = bodyPropertyParentTreeItem->child(i) ;
		if( child->data(0, Qt::UserRole) == (qulonglong) dialog){
			tree->setCurrentItem(child);
			BodyPropertyEditor* pe = (BodyPropertyEditor*) dialog;
			connect1( pe,SIGNAL(BodyMaterialComboChanged(BodyPropertyEditor *,QString)), this, SLOT(bodyComboChanged(BodyPropertyEditor *,QString)) );
			connect1( pe,SIGNAL(BodyInitialComboChanged(BodyPropertyEditor *,QString)), this, SLOT(bodyComboChanged(BodyPropertyEditor *,QString)) );
			connect1( pe,SIGNAL(BodyForceComboChanged(BodyPropertyEditor *,QString)), this, SLOT(bodyComboChanged(BodyPropertyEditor *,QString)) );
			connect1( pe,SIGNAL(BodyEquationComboChanged(BodyPropertyEditor *,QString)), this, SLOT(bodyComboChanged(BodyPropertyEditor *,QString)) );
      connect1( pe->ui.applyButton,SIGNAL(clicked(bool)), this, SLOT(bodyPropertyEditorAccepted(bool)));
      connect1( pe->ui.discardButton,SIGNAL(clicked(bool)), this, SLOT(bodyPropertyEditorDiscarded(bool)));	
			return;
		}
		
	}

}



list_t* ObjectBrowser::boundaryList(int index){

  MainWindow* mainwindow = (MainWindow*) mainWindow; 
  
	int n = mainwindow->glWidget->getLists();
	int boundaryCount = 0;
	for(int i = 0; i < n; i++){
		list_t* l = mainwindow->glWidget->getList(i);
		if(l->getNature() == PDE_BOUNDARY && l->getIndex() == index) return l;		
	}
	return NULL;
}

list_t* ObjectBrowser::selectBoundary(BoundaryPropertyEditor *pe_in_tree, bool append/*=false*/, bool select /*=true*/){

  MainWindow* mainwindow = (MainWindow*) mainWindow;
  if(mainwindow->glWidget->getMesh() == NULL) return NULL;

   // unselect all the bodies
  for( int i=0; i< mainwindow->glWidget->bodyMap.count(); i++ ) {
    int n = mainwindow->glWidget->bodyMap.key(i);
    if(n >= 0 && bodyList(n) != NULL) bodyList(n)->setSelected(false);
		if(bodyList(n) == NULL){
//			cout << " *** bodyList(" << n << ") is NULL" << endl;
//			cout << " *** body index is " << mainwindow->glWidget->bodyMap.value(i) << endl;			
		}
  }

  list_t* l = NULL; 
  list_t* l_selected = NULL;
  //select boundary
  for( int i=0; i< mainwindow->glWidget->boundaryMap.count(); i++ ) {
    int n = mainwindow->glWidget->boundaryMap.key(i);
    if ( n >= 0 ) {
      int m = mainwindow->glWidget->boundaryMap.value(n);     
      if(m >= 0 && m < mainwindow->boundaryPropertyEditor.size() && mainwindow->boundaryPropertyEditor[m] != NULL) {
        BoundaryPropertyEditor *pe = mainwindow->boundaryPropertyEditor[m];	
        l = boundaryList(n);
        if(l != NULL /*&& pe_in_tree != NULL*/){
          if(pe == pe_in_tree) {l->setSelected(select); l_selected = l;}
          if(pe != pe_in_tree && append == false){l->setSelected(!select);} 
        }
      }
    }
  }
  
  if(!append){
    mainwindow->glWidget->rebuildEdgeLists();
    mainwindow->glWidget->rebuildSurfaceLists();
    mainwindow->glWidget->updateGL();  
  }
  return l_selected;  
}

list_t* ObjectBrowser::selectBody(BodyPropertyEditor *pe_in_tree, bool append /*=false*/, bool select /*=true*/){

  MainWindow* mainwindow = (MainWindow*) mainWindow;
  if(mainwindow->glWidget->getMesh() == NULL) return NULL;  

  //unselect all the boundaries
  if(!append){
    for( int i=0; i< mainwindow->glWidget->boundaryMap.count(); i++ ) {
      int n = mainwindow->glWidget->boundaryMap.key(i);
      if ( n >= 0 && boundaryList(n) != NULL) boundaryList(n)->setSelected(false);
      if(boundaryList(n) == NULL) cout << " *** boundaryList(" << n << ") is NULL" << endl;	
    }
  }
  mainwindow->glWidget->currentlySelectedBody = -1;
   
  list_t* l = NULL; 
  list_t* l_selected = NULL;
  
  //select body
  for( int i=0; i< mainwindow->glWidget->bodyMap.count(); i++ ) {
    int n = mainwindow->glWidget->bodyMap.key(i);
    if ( n >= 0 ) {
      int m = mainwindow->glWidget->bodyMap.value(n);
      if(m >= 0 && m < mainwindow->bodyPropertyEditor.size() && mainwindow->bodyPropertyEditor[m] != NULL){
        BodyPropertyEditor *pe = mainwindow->bodyPropertyEditor[m];
        l = bodyList(n);
		
        //new
        if(l != NULL /*&& pe_in_tree != NULL*/){
          if(pe == pe_in_tree) {l->setSelected(select); l_selected = l; mainwindow->glWidget->currentlySelectedBody = select ? n : -1;}
          if(pe != pe_in_tree && append == false){l->setSelected(!select);} 
        }
      
      
        //oled  
        //if(l != NULL)l->setSelected(/*pe_in_tree != NULL &&*/ pe == pe_in_tree); 
        //if(l != NULL /*&& pe_in_tree != NULL*/ && pe == pe_in_tree){
        //l_selected = l;
			
      }
	  }
	}

  if(l_selected == NULL){
  //For 3D geometry
    mainwindow->glWidget->currentlySelectedBody = -1;
    int index_selected= -1;
    if(l_selected == NULL){
      for( int i=0; i< mainwindow->glWidget->boundaryMap.count(); i++ ) {
        int n = mainwindow->glWidget->boundaryMap.key(i);
        if ( n >= 0 ) {
          int m = mainwindow->glWidget->boundaryMap.value(n);     
          if(m >= 0 && m < mainwindow->boundaryPropertyEditor.size() && mainwindow->boundaryPropertyEditor[m] != NULL){
            BoundaryPropertyEditor *pe = mainwindow->boundaryPropertyEditor[m];
            l = boundaryList(n);
            if(l != NULL /*&& pe_in_tree != NULL*/){
              int index = boundaryListToBodyIndex(l);
              if( index >= 0){
                int m2 = mainwindow->glWidget->bodyMap.value(index);
                if(m2 >= 0 && m2 < mainwindow->bodyPropertyEditor.size() && mainwindow->bodyPropertyEditor[m2] != NULL){
                  BodyPropertyEditor *bpe = mainwindow->bodyPropertyEditor[m2];
                  if(bpe == pe_in_tree){
                    l_selected = l;
                    l = selectBoundary(pe, true, select);
                    index_selected = index;
                  }
                }
              }
            }
          }
        }
      }
      mainwindow->glWidget->currentlySelectedBody = select ? index_selected : -1;
    }
   }
   
  if(!append){
    mainwindow->glWidget->rebuildEdgeLists();
    mainwindow->glWidget->rebuildSurfaceLists();
    mainwindow->glWidget->updateGL();
  }
  return l_selected;  
}

int ObjectBrowser::boundaryListToBodyIndex(list_t* l){

	MainWindow* mainwindow = (MainWindow*) mainWindow;  
	mesh_t* mesh = mainwindow->glWidget->getMesh();
	int ret = -1;

	if(true){
      // determine the max bulk index
      int MAX_BULK_INDEX = -1;
		

      for(int i = 0; i < mesh->getElements(); i++) {
	element_t *elem = mesh->getElement(i);
	if(elem->getNature() != PDE_BULK)
	  break;
	if(elem->getIndex() > MAX_BULK_INDEX)
	  MAX_BULK_INDEX = elem->getIndex();
      }

      for(int i = 0; i < mesh->getSurfaces(); i++) {
	surface_t *surf = mesh->getSurface(i);
	if(surf->getNature() != PDE_BULK)
	  break;
	if(surf->getIndex() > MAX_BULK_INDEX)
	  MAX_BULK_INDEX = surf->getIndex();
      }

      for(int i = 0; i < mesh->getEdges(); i++) {
	edge_t *edge = mesh->getEdge(i);
	if(edge->getNature() != PDE_BULK)
	  break;
	if(edge->getIndex() > MAX_BULK_INDEX)
	  MAX_BULK_INDEX = edge->getIndex();
      }
      
      MAX_BULK_INDEX++;
      if(MAX_BULK_INDEX == 0) {
	cout << "Error in body selection: "
	  "There are no legal body indiced from which to choose" << endl;
	cout.flush();
	goto body_selection_finished;
      }

      // allocate temp arrays:
      bool *tmp1 = new bool[MAX_BULK_INDEX];
      bool *tmp2 = new bool[MAX_BULK_INDEX];

      for(int i = 0; i < MAX_BULK_INDEX; i++) {
	tmp1[i] = true;
	tmp2[i] = false;
      }
      
      // check if the selected lists uniquely determine a bulk body:
      for(int i = 0; i < mainwindow->glWidget->getLists(); i++) {
	list_t *l2 = mainwindow->glWidget->getList(i);

	if(l2 == l && (l2->getNature() == PDE_BULK)) { //:-)	if(l2->isSelected() && (l2->getNature() == PDE_BULK)) {
	  for(int j = 0; j < MAX_BULK_INDEX; j++) {
	    if(j != l2->getIndex())
	      tmp1[j] = false;
	  }
	}
	
	if(l2 == l && (l2->getNature() == PDE_BOUNDARY) &&//:-)	if(l2->isSelected() && (l2->getNature() == PDE_BOUNDARY) && 
	   (l2->getType() == SURFACELIST)) {	  
	  for(int j = 0; j < mesh->getSurfaces(); j++) {
	    surface_t *surf = mesh->getSurface(j);
	    if(surf->getIndex() == l2->getIndex()) {
	      for(int k = 0; k < surf->getElements(); k++) {
		int l = surf->getElementIndex(k);
		if(l < 0) 
		  break;
		element_t *elem = mesh->getElement(l);
		if((elem->getIndex() < 0) || (elem->getIndex() >= MAX_BULK_INDEX))
		  break;
		tmp2[elem->getIndex()] = true;
	      }
	      for(int k = 0; k < MAX_BULK_INDEX; k++) {
		tmp1[k] &= tmp2[k];
		tmp2[k] = false;
	      }
	    }
	  }
	}
      }

      // array "tmp1" should contain only one entry with value "true"
      int count = 0;
      int found = -1;
      for(int i = 0; i < MAX_BULK_INDEX; i++) {
	if( tmp1[i] ) {
	  count++;
	  found = i;
	}
      }
      
      if((count == 1) && (found >= 0)){
		mainwindow->glWidget->currentlySelectedBody = found;
		ret = found; 
		//cout << "*****boundary found: " << found << endl;
	}
      delete [] tmp1;
      delete [] tmp2;
   }
  body_selection_finished:
    
    // Emit result to mainwindow:
    //emit(signalBoundarySelected(l));
	
	return ret;

}

void ObjectBrowser::updateBodyProperties(BodyPropertyEditor* selectThis /*=NULL*/){

	MainWindow* mainwindow = (MainWindow*) mainWindow;
	if(tree->currentItem() != NULL && tree->currentItem()->parent() == bodyPropertyParentTreeItem){
    tree->setCurrentItem(bodyPropertyParentTreeItem);// to improve speed
  }
	QTreeWidgetItem* treeItem = 0;
	while(treeItem = bodyPropertyParentTreeItem->takeChild(0)) delete treeItem;

	BodyPropertyEditor* pe;
	int count = 1;
	
  for( int i=0; i< mainwindow->glWidget->bodyMap.count(); i++ )
  {
    int n = mainwindow->glWidget->bodyMap.key(i);
    if ( n >= 0 ) {
      int m = mainwindow->glWidget->bodyMap.value(n);

      if(m >= mainwindow->bodyPropertyEditor.size()) mainwindow->bodyPropertyEditor.resize(m + 1);

      if(!mainwindow->bodyPropertyEditor[m])  mainwindow->bodyPropertyEditor[m] = new BodyPropertyEditor;

      BodyPropertyEditor *pe = mainwindow->bodyPropertyEditor[m];

      treeItem = new QTreeWidgetItem();
      bodyPropertyParentTreeItem->addChild(treeItem);
      QString title = pe->ui.nameEdit->text().trimmed();
      if(title.isEmpty()) {
        title = "Body Property " + QString::number(n);
      }
      
      if(!pe->touched){
        treeItem->setText(0, title);
        if(pe->equation || pe->material || pe->force || pe->initial) treeItem->setText(1, "<Canceled>");
      }else{
        QString sifbody = "Body " + QString::number(count++);
        QString sifname = title;
        if(pe->ui.nameEdit->text().trimmed().isEmpty()){
          sifname = sifbody;
        }
        treeItem->setText(0, title);
        //treeItem->setText(1, "<" + sifbody + " \"" + sifname +"\">");
        treeItem->setText(1,  sifbody + " in sif");
      }

      treeItem->setData(0, Qt::UserRole, QVariant::fromValue( (qulonglong)pe));

      if(pe->equation != 0){
        QTreeWidgetItem* child = new QTreeWidgetItem();
        treeItem->addChild(child);
        child->setText(0, "Equation");			
        child->setText(1, pe->equation->nameEdit->text().trimmed());
        child->setData(1, Qt::UserRole, QVariant::fromValue( (qulonglong)pe->equation));	
      }
      if(pe->material != 0){
        QTreeWidgetItem* child = new QTreeWidgetItem();
        treeItem->addChild(child);
        child->setText(0, "Material");	
        child->setText(1, pe->material->nameEdit->text().trimmed());
        child->setData(1, Qt::UserRole, QVariant::fromValue( (qulonglong)pe->material));	
      }
      if(pe->force != 0){
        QTreeWidgetItem* child = new QTreeWidgetItem();
        treeItem->addChild(child);
        child->setText(0, "Body force");	
        child->setText(1, pe->force->nameEdit->text().trimmed());
        child->setData(1, Qt::UserRole, QVariant::fromValue( (qulonglong)pe->force));	
      }
      if(pe->initial != 0){
        QTreeWidgetItem* child = new QTreeWidgetItem();
        treeItem->addChild(child);
        child->setText(0, "Initial condition");	
        child->setText(1, pe->initial->nameEdit->text().trimmed());
        child->setData(1, Qt::UserRole, QVariant::fromValue( (qulonglong)pe->initial));	
      }
      treeItem->setExpanded(true);
    }
  }	 
  
  if(selectThis != NULL){  
	QTreeWidgetItem* item;
	for(int i = 0; i < bodyPropertyParentTreeItem->childCount(); i++){
		item = bodyPropertyParentTreeItem->child(i);
		if( item->data(0, Qt::UserRole).value<qulonglong>() == (qulonglong)selectThis ){
			tree->setCurrentItem(item);
			return;
		}
	}
  }

  /* 
  cout << "solverparametereditor:" << mainwindow->solverParameterEditor.size() << endl;
  for(int i= 0; i < mainwindow->solverParameterEditor.size(); i++){
    if(mainwindow->solverParameterEditor.at(i)) mainwindow->logMessage(  mainwindow->solverParameterEditor.at(i)->solverName );
    else cout << "null" << endl;
  }
  */
}


list_t* ObjectBrowser::bodyList(int index){

  MainWindow* mainwindow = (MainWindow*) mainWindow; 
  
	int n = mainwindow->glWidget->getLists();
	int bodyCount = 0;
	for(int i = 0; i < n; i++){
		list_t* l = mainwindow->glWidget->getList(i);
		if(l->getNature() == PDE_BULK && l->getIndex() == index) return l;
	}
	return NULL;
}

void ObjectBrowser::snap(QWidget* widget){
	//widget->move(mapToGlobal(QPoint(width(),0)));
}


void ObjectBrowser::updateEquation(){
  tree->setCurrentItem(equationParentTreeItem);
	QTreeWidgetItem* treeItem;
	while(treeItem = equationParentTreeItem->takeChild(0)) delete treeItem;
	MainWindow* mainwindow = (MainWindow*) mainWindow;
	int n = mainwindow->equationEditor.size();
	for(int i= 0; i < n; i++){
		addToTree(mainwindow->equationEditor.at(i));
	}
}
void ObjectBrowser::updateMaterial(){
  tree->setCurrentItem(materialParentTreeItem);
	QTreeWidgetItem* treeItem;
	while(treeItem = materialParentTreeItem->takeChild(0)) delete treeItem;
	MainWindow* mainwindow = (MainWindow*) mainWindow;
	int n = mainwindow->materialEditor.size();
	for(int i= 0; i < n; i++){
		addToTree(mainwindow->materialEditor.at(i));
	}
}
void ObjectBrowser::updateBodyForce(){
  tree->setCurrentItem(bodyForceParentTreeItem);
	QTreeWidgetItem* treeItem;
	while(treeItem = bodyForceParentTreeItem->takeChild(0)) delete treeItem;
	MainWindow* mainwindow = (MainWindow*) mainWindow;
	int n = mainwindow->bodyForceEditor.size();
	for(int i= 0; i < n; i++){
		addToTree(mainwindow->bodyForceEditor.at(i));
	}
}
void ObjectBrowser::updateInitialCondition(){
  tree->setCurrentItem(initialConditionParentTreeItem);
	QTreeWidgetItem* treeItem;
	while(treeItem = initialConditionParentTreeItem->takeChild(0)) delete treeItem;
	MainWindow* mainwindow = (MainWindow*) mainWindow;
	int n = mainwindow->initialConditionEditor.size();
	for(int i= 0; i < n; i++){
		addToTree(mainwindow->initialConditionEditor.at(i));
	}
}
void ObjectBrowser::updateBoundaryCondition(){
  tree->setCurrentItem(boundaryConditionParentTreeItem);
	QTreeWidgetItem* treeItem;
	while(treeItem = boundaryConditionParentTreeItem->takeChild(0)) delete treeItem;
	MainWindow* mainwindow = (MainWindow*) mainWindow;
	int n = mainwindow->boundaryConditionEditor.size();
	for(int i= 0; i < n; i++){
		addToTree(mainwindow->boundaryConditionEditor.at(i));
	}
}
void ObjectBrowser::boundaryComboChanged(BoundaryPropertyEditor *pe,QString text){
	QTreeWidgetItem* treeItem;
	QTreeWidgetItem* child;
	
	int n = boundaryPropertyParentTreeItem->childCount();
	for(int i= 0; i < n; i++){
		treeItem = boundaryPropertyParentTreeItem->child(i);
		if( treeItem->data(0, Qt::UserRole) == (qulonglong) pe){
			if(pe->condition != 0){
				treeItem->setText(1, pe->condition->nameEdit->text().trimmed());
			}else{
				treeItem->setText(1, "");
			}
			treeItem->setData(1, Qt::UserRole, QVariant::fromValue( (qulonglong)pe->condition));
			return;
		}	
	}
}
	
void ObjectBrowser::bodyComboChanged(BodyPropertyEditor *pe,QString text){
	MainWindow* mainwindow = (MainWindow*) mainWindow;
	QTreeWidgetItem* treeItem;
	QTreeWidgetItem* child;
	int n = bodyPropertyParentTreeItem->childCount();
	for(int i= 0; i < n; i++){
	  treeItem = bodyPropertyParentTreeItem->child(i);
	  if( treeItem->data(0, Qt::UserRole) == (qulonglong) pe){
	    while(child = treeItem->takeChild(0)) delete child;
		
		if(pe->equation != 0){
			child = new QTreeWidgetItem();
			treeItem->addChild(child);
			child->setText(0, "Equation");
			child->setText(1, pe->equation->nameEdit->text().trimmed());
			child->setData(1, Qt::UserRole, QVariant::fromValue( (qulonglong)pe->equation));
		}			

		if(pe->material != 0){
			child = new QTreeWidgetItem();
			treeItem->addChild(child);
			child->setText(0, "Material");
			child->setText(1, pe->material->nameEdit->text().trimmed());
			child->setData(1, Qt::UserRole, QVariant::fromValue( (qulonglong)pe->material));	
		}

		if(pe->force != 0){
			child = new QTreeWidgetItem();
			treeItem->addChild(child);
			child->setText(0, "Body force");
			child->setText(1, pe->force->nameEdit->text().trimmed());
			child->setData(1, Qt::UserRole, QVariant::fromValue( (qulonglong)pe->force));	
		}

		if(pe->initial != 0) {
			child = new QTreeWidgetItem();
			treeItem->addChild(child);
			child->setText(0, "Initial condition");
			child->setText(1, pe->initial->nameEdit->text().trimmed());
			child->setData(1, Qt::UserRole, QVariant::fromValue( (qulonglong)pe->initial));	
		}
		
		return;
	  }
	}
}

void ObjectBrowser::selectTreeItemByID(QTreeWidgetItem* parent, int id){
	DynamicEditor* de;
	if(id >= 0){
		for(int i = 0; i < parent->childCount(); i++){
			de = (DynamicEditor*) parent->child(i)->data(0, Qt::UserRole).value<qulonglong>();
			if(de == NULL) cout << "*** de==NULL" << endl;
			if( de != NULL && de->ID == id){
				tree->setCurrentItem(parent->child(i));
				return;
			}
		}
	}
	tree->setCurrentItem(parent);			
}

void ObjectBrowser::bodyPropertyEditorAccepted(bool checked){
	QWidget *a = (QWidget *)QObject::sender();
	updateBodyProperties((BodyPropertyEditor*)a->parent());
}

void ObjectBrowser::bodyPropertyEditorDiscarded(bool checked){
	QWidget *a = (QWidget *)QObject::sender();
	updateBodyProperties((BodyPropertyEditor*)a->parent());
}
void ObjectBrowser::boundaryPropertyEditorAccepted(bool checked){
	QWidget *a = (QWidget *)QObject::sender();
	updateBoundaryProperties((BoundaryPropertyEditor*)a->parent());
}
void ObjectBrowser::boundaryPropertyEditorDiscarded(bool checked){
	QWidget *a = (QWidget *)QObject::sender();
	updateBoundaryProperties((BoundaryPropertyEditor*)a->parent());
}

void ObjectBrowser::viewFullScreenSlot(){
	if(isVisible()) hide();
	else show();
}

void ObjectBrowser::viewNormalModeSlot(){
	show();
} 

bool ObjectBrowser::connect1(const QObject * sender, const char * signal, const QObject * receiver, const char * method, Qt::ConnectionType type /*= Qt::AutoConnection*/){
  disconnect(sender, signal, receiver, method);
  return connect(sender, signal, receiver, method, type);
}

void ObjectBrowser::bodyCheckBoxChangedSlot(int state){
  updateBodyProperties();

  QWidget *a = (QWidget *)QObject::sender();
  if(a == NULL) return;
  BodyPropertyEditor* pe = (BodyPropertyEditor*) a->property("body").toULongLong();
  if(pe == NULL) return;
  
  selectBody(pe, true, state == Qt::Checked);
	MainWindow* mainwindow = (MainWindow*) mainWindow;  
  mainwindow->glWidget->rebuildEdgeLists();
  mainwindow->glWidget->rebuildSurfaceLists();
  mainwindow->glWidget->updateGL();  
}

void ObjectBrowser::boundaryCheckBoxChangedSlot(int state){
  updateBoundaryProperties();

  QWidget *a = (QWidget *)QObject::sender();
  if(a == NULL) return;
  BoundaryPropertyEditor *pe = (BoundaryPropertyEditor*) a->property("boundary").toULongLong();
  if(pe == NULL) return;

  selectBoundary(pe, true, state == Qt::Checked);
	MainWindow* mainwindow = (MainWindow*) mainWindow;
  mainwindow->glWidget->rebuildEdgeLists();
  mainwindow->glWidget->rebuildSurfaceLists();
  mainwindow->glWidget->updateGL();  
}

void ObjectBrowser::treeItemSelectionChangedSlot(){
/*
  if(item == bodyPropertyParentTreeItem){
    updateBodyProperties();
    return;
  }
  if(item == boundaryPropertyParentTreeItem){
    updateBoundaryProperties();
    return;
  }
  */

	MainWindow* mainwindow = (MainWindow*) mainWindow;
	if(mainwindow->glWidget->getMesh() == NULL) return;
	if(mainwindow->glWidget->ctrlPressed) return; // ignore if selecting multiple surfaces/edges in glWidget to maintain the selection at glWidget
	QTreeWidgetItem* item = tree->currentItem();
  
    //select boundaries/bodies checked in DyanmicEditor
  if(item->parent() == equationParentTreeItem){
    updateBodyProperties();
    updateBoundaryProperties();  
    DynamicEditor *pe = (DynamicEditor*) item->data(0, Qt::UserRole).value<qulonglong>(); 
    mainwindow->createBodyCheckBoxes(BODY_EQUATION, pe);
    QGroupBox *box = (QGroupBox*) pe->spareScroll->widget();
    selectBody(NULL);
    for(int i=2; i< box->children().size(); i++){
      QCheckBox* cb = (QCheckBox*) box->children()[i];
      connect1(cb, SIGNAL(stateChanged(int)), this, SLOT(bodyCheckBoxChangedSlot(int)));
      if(cb->isChecked() && i-2 < mainwindow->bodyPropertyEditor.size()){
        selectBody(mainwindow->bodyPropertyEditor[i-2], true);
      }
    }
    mainwindow->glWidget->rebuildEdgeLists();
    mainwindow->glWidget->rebuildSurfaceLists();
    mainwindow->glWidget->updateGL();
    return;
  } 
  if(item->parent() == materialParentTreeItem){
    updateBodyProperties();
    updateBoundaryProperties();  
  	DynamicEditor *pe = (DynamicEditor*) item->data(0, Qt::UserRole).value<qulonglong>();
  	mainwindow->createBodyCheckBoxes(BODY_MATERIAL, pe);
    QGroupBox *box = (QGroupBox*) pe->spareScroll->widget();
    selectBody(NULL);
    for(int i=2; i< box->children().size(); i++){
      QCheckBox* cb = (QCheckBox*) box->children()[i];
      connect1(cb, SIGNAL(stateChanged(int)), this, SLOT(bodyCheckBoxChangedSlot(int)));
      if(cb->isChecked()) selectBody((BodyPropertyEditor*) bodyPropertyParentTreeItem->child(i-2)->data(0, Qt::UserRole).value<qulonglong>(), true);
    }
    mainwindow->glWidget->rebuildEdgeLists();
    mainwindow->glWidget->rebuildSurfaceLists();
    mainwindow->glWidget->updateGL();
    return;
  }
  if(item->parent() == bodyForceParentTreeItem){
    updateBodyProperties();
    updateBoundaryProperties();  
  	DynamicEditor *pe = (DynamicEditor*) item->data(0, Qt::UserRole).value<qulonglong>();
    mainwindow->createBodyCheckBoxes(BODY_FORCE, pe);
    QGroupBox *box = (QGroupBox*) pe->spareScroll->widget();
    selectBody(NULL);
    for(int i=2; i< box->children().size(); i++){
      QCheckBox* cb = (QCheckBox*) box->children()[i];
      connect1(cb, SIGNAL(stateChanged(int)), this, SLOT(bodyCheckBoxChangedSlot(int)));
      if(cb->isChecked()) selectBody((BodyPropertyEditor*) bodyPropertyParentTreeItem->child(i-2)->data(0, Qt::UserRole).value<qulonglong>(), true);
    }
    mainwindow->glWidget->rebuildEdgeLists();
    mainwindow->glWidget->rebuildSurfaceLists();
    mainwindow->glWidget->updateGL();         
    return;
  }
  if(item->parent() == initialConditionParentTreeItem){
    updateBodyProperties();
    updateBoundaryProperties();  
  	DynamicEditor *pe = (DynamicEditor*) item->data(0, Qt::UserRole).value<qulonglong>();
    mainwindow->createBodyCheckBoxes(BODY_INITIAL, pe);
    QGroupBox *box = (QGroupBox*) pe->spareScroll->widget();
    selectBody(NULL);    
    for(int i=2; i< box->children().size(); i++){
      QCheckBox* cb = (QCheckBox*) box->children()[i];
      connect1(cb, SIGNAL(stateChanged(int)), this, SLOT(bodyCheckBoxChangedSlot(int)));
      if(cb->isChecked()) selectBody((BodyPropertyEditor*) bodyPropertyParentTreeItem->child(i-2)->data(0, Qt::UserRole).value<qulonglong>(), true);
    }
    mainwindow->glWidget->rebuildEdgeLists();
    mainwindow->glWidget->rebuildSurfaceLists();
    mainwindow->glWidget->updateGL();     
    return;
  }
  if(item->parent() == boundaryConditionParentTreeItem){
    updateBodyProperties();
    updateBoundaryProperties();  
  	DynamicEditor *pe = (DynamicEditor*) item->data(0, Qt::UserRole).value<qulonglong>();
    mainwindow->createBoundaryCheckBoxes(pe);
    QGroupBox *box = (QGroupBox*) pe->spareScroll->widget();
    selectBoundary(NULL);
    for(int i=2; i< box->children().size(); i++){
      QCheckBox* cb = (QCheckBox*) box->children()[i];
      connect1(cb, SIGNAL(stateChanged(int)), this, SLOT(boundaryCheckBoxChangedSlot(int)));
      if(cb->isChecked()) selectBoundary((BoundaryPropertyEditor*) boundaryPropertyParentTreeItem->child(i-2)->data(0, Qt::UserRole).value<qulonglong>(), true);
    }
    mainwindow->glWidget->rebuildEdgeLists();
    mainwindow->glWidget->rebuildSurfaceLists();
    mainwindow->glWidget->updateGL(); 
    return;
  }

  //select boundary/body
  if(item->parent() == boundaryPropertyParentTreeItem ){
 	BoundaryPropertyEditor* pe = (BoundaryPropertyEditor*) item->data(0, Qt::UserRole).value<qulonglong>();
	if(pe == NULL ){ cout << " BoundaryPropertyEditor NULL" << endl; return;}
	selectBoundary(pe);
	return;
  }
  if(item->parent() != NULL && item->parent()->parent() == boundaryPropertyParentTreeItem){
 	BoundaryPropertyEditor* pe = (BoundaryPropertyEditor*) item->parent()->data(0, Qt::UserRole).value<qulonglong>();
	if(pe == NULL){ cout << " BoundaryPropertyEditor NULL" << endl; return;}
	selectBoundary(pe);
	return;
  }
  if(item->parent() == bodyPropertyParentTreeItem){
 	BodyPropertyEditor* pe = (BodyPropertyEditor*) item->data(0, Qt::UserRole).value<qulonglong>();
	if(pe == NULL){ cout << " BodyPropertyEditor NULL" << endl; return;}
	selectBody(pe);
	return;
  }
  if(item->parent() != NULL && item->parent()->parent() == bodyPropertyParentTreeItem){
 	BodyPropertyEditor* pe = (BodyPropertyEditor*) item->parent()->data(0, Qt::UserRole).value<qulonglong>();
	if(pe == NULL){ cout << " BodyPropertyEditor NULL" << endl; return;}
	selectBody(pe);
	return;
  }

  selectBoundary(NULL);
  selectBody(NULL);
 
}

void ObjectBrowser::meshingStartedSlot()
{
}
void ObjectBrowser::meshingTerminatedSlot()
{
}
void ObjectBrowser::meshingFinishedSlot()
{
	bodyPropertyParentTreeItem->setExpanded(false);
	boundaryPropertyParentTreeItem->setExpanded(false);
}
