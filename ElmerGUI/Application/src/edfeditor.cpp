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

#include <QtGui>
#include <iostream>
#include "edfeditor.h"

using namespace std;

// ctor...
//----------------------------------------------------------------------------
EdfEditor::EdfEditor(QWidget *parent)
  : QWidget(parent)
{
  addIcon = QIcon(":/icons/list-add.png");
  removeIcon = QIcon(":/icons/list-remove.png");
  collapseIcon = QIcon(":/icons/arrow-up.png");
  expandIcon = QIcon(":/icons/arrow-down.png");
  openIcon = QIcon(":/icons/document-open.png");
  appendIcon = QIcon(":/icons/tab-new-background.png"); // todo
  saveAsIcon = QIcon(":/icons/document-save.png");
  applyIcon = QIcon(":/icons/dialog-close.png");
  previewIcon = QIcon(":/icons/document-preview.png");

  lastActiveItem = NULL;
  ctrlPressed = false;

  // Set up tree widget:
  //--------------------
  edfTree = new QTreeWidget;

  connect(edfTree, SIGNAL(itemClicked(QTreeWidgetItem*,int)),
	  this, SLOT(treeItemClicked(QTreeWidgetItem*,int)));

  edfTree->setColumnCount(3);
  edfTree->setColumnWidth(0,200);
  edfTree->setColumnWidth(1,200);
  edfTree->setColumnWidth(2,200);

  edfTree->setAnimated(true);

  // Set internal drag'n drop mode on:
  //----------------------------------
  edfTree->setDragEnabled(true);
  edfTree->setDragDropMode(QAbstractItemView::InternalMove);
  edfTree->setDropIndicatorShown(true);
  edfTree->setDragDropOverwriteMode(false);

  QStringList qsl;
  qsl << "Tag" << "Attributes" << "Value";
  edfTree->setHeaderLabels(qsl);
  edfTree->setAlternatingRowColors(true);

  // Buttons:
  //---------
  addButton = new QPushButton(tr("&Add child"));
  addButton->setIcon(addIcon);
  connect(addButton, SIGNAL(clicked()), this, SLOT(addButtonClicked()));
  
  removeButton = new QPushButton(tr("&Remove item"));
  removeButton->setIcon(removeIcon);
  connect(removeButton, SIGNAL(clicked()), this, SLOT(removeButtonClicked()));

  expandCollapseAllButton = new QPushButton(tr("Collapse all"));
  expandCollapseAllButton->setIcon(collapseIcon);
  connect(expandCollapseAllButton, SIGNAL(clicked()),
	  this, SLOT(expandCollapseAllButtonClicked()));

  openButton = new QPushButton(tr("&Open"));
  openButton->setIcon(openIcon);
  connect(openButton, SIGNAL(clicked()), this, SLOT(openButtonClicked()));

  appendButton = new QPushButton(tr("&Append"));
  appendButton->setIcon(appendIcon);
  connect(appendButton, SIGNAL(clicked()), this, SLOT(appendButtonClicked()));

  previewButton = new QPushButton(tr("&Preview"));
  previewButton->setIcon(previewIcon);
  connect(previewButton, SIGNAL(clicked()), this, SLOT(previewButtonClicked()));

  saveAsButton = new QPushButton(tr("&Save as"));
  saveAsButton->setIcon(saveAsIcon);
  connect(saveAsButton, SIGNAL(clicked()), this, SLOT(saveAsButtonClicked()));

  applyButton = new QPushButton(tr("&Close"));
  applyButton->setIcon(applyIcon);
  connect(applyButton, SIGNAL(clicked()), this, SLOT(applyButtonClicked()));

  QHBoxLayout *buttonLayout = new QHBoxLayout;  
  buttonLayout->addWidget(addButton);
  buttonLayout->addWidget(removeButton);
  buttonLayout->addWidget(expandCollapseAllButton);
  buttonLayout->addWidget(openButton);
  buttonLayout->addWidget(appendButton);
  buttonLayout->addWidget(previewButton);
  buttonLayout->addWidget(saveAsButton);
  buttonLayout->addWidget(applyButton);

  // Main layout:
  //-------------
  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(edfTree);
  mainLayout->addLayout(buttonLayout);
  setLayout(mainLayout);

  setWindowTitle("Elmer Definitions File editor");

  setFocusPolicy(Qt::ClickFocus);

  expandCollapseAll = false;

  dynamicEditorSimulation = new DynamicEditor;
  dynamicEditorConstants = new DynamicEditor;
  dynamicEditorEquation = new DynamicEditor;
  dynamicEditorSolver = new DynamicEditor;
  dynamicEditorMaterial = new DynamicEditor;
  dynamicEditorBodyForce = new DynamicEditor;
  dynamicEditorBC = new DynamicEditor;
  dynamicEditorIC = new DynamicEditor;

  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
}

// dtor...
//----------------------------------------------------------------------------
EdfEditor::~EdfEditor()
{
}

// Min window size...
//----------------------------------------------------------------------------
QSize EdfEditor::minimumSizeHint() const
{
  return QSize(64, 64);
}

// Default window size...
//----------------------------------------------------------------------------
QSize EdfEditor::sizeHint() const
{
  return QSize(720, 480);
}

// preview panels
//-----------------------------------------------------------------------------
void EdfEditor::previewButtonClicked()
{
  if(elmerDefs == NULL)
    return;
  
  // always create a new instance:
  QMainWindow *sandBox = new QMainWindow;
  QMdiArea *mdiArea = new QMdiArea;

  sandBox->setCentralWidget(mdiArea);
  sandBox->setWindowTitle(tr("Preview definitions"));
  sandBox->show();

  if(dynamicEditorSimulation) delete dynamicEditorSimulation;
  dynamicEditorSimulation = new DynamicEditor;
  dynamicEditorSimulation->setupTabs(elmerDefs, "Simulation",1);
  QMdiSubWindow *simulationSubWindow = new QMdiSubWindow;
  simulationSubWindow->setWidget(dynamicEditorSimulation);
  mdiArea->addSubWindow(simulationSubWindow);
  simulationSubWindow->show();

  if(dynamicEditorConstants) delete dynamicEditorConstants;
  dynamicEditorConstants = new DynamicEditor;
  dynamicEditorConstants->setupTabs(elmerDefs, "Constants",1);
  QMdiSubWindow *constantsSubWindow = new QMdiSubWindow;
  constantsSubWindow->setWidget(dynamicEditorConstants);
  mdiArea->addSubWindow(constantsSubWindow);
  constantsSubWindow->show();

  if(dynamicEditorEquation) delete dynamicEditorEquation;
  dynamicEditorEquation = new DynamicEditor;
  dynamicEditorEquation->setupTabs(elmerDefs, "Equation",1);
  QMdiSubWindow *equationSubWindow = new QMdiSubWindow;
  equationSubWindow->setWidget(dynamicEditorEquation);
  mdiArea->addSubWindow(equationSubWindow);
  equationSubWindow->show();

  if(dynamicEditorSolver) delete dynamicEditorSolver;
  dynamicEditorSolver = new DynamicEditor;
  dynamicEditorSolver->setupTabs(elmerDefs, "Solver",1 );
  QMdiSubWindow *solverSubWindow = new QMdiSubWindow;
  solverSubWindow->setWidget(dynamicEditorSolver);
  mdiArea->addSubWindow(solverSubWindow);
  solverSubWindow->show();

  if(dynamicEditorMaterial) delete dynamicEditorMaterial;
  dynamicEditorMaterial = new DynamicEditor;
  dynamicEditorMaterial->setupTabs(elmerDefs, "Material",1 );
  QMdiSubWindow *materialSubWindow = new QMdiSubWindow;
  materialSubWindow->setWidget(dynamicEditorMaterial);
  mdiArea->addSubWindow(materialSubWindow);
  materialSubWindow->show();

  if(dynamicEditorBodyForce) delete dynamicEditorBodyForce;
  dynamicEditorBodyForce = new DynamicEditor;
  dynamicEditorBodyForce->setupTabs(elmerDefs, "BodyForce",1 );
  QMdiSubWindow *bodyForceSubWindow = new QMdiSubWindow;
  bodyForceSubWindow->setWidget(dynamicEditorBodyForce);
  mdiArea->addSubWindow(bodyForceSubWindow);
  bodyForceSubWindow->show();

  if(dynamicEditorIC) delete dynamicEditorIC;
  dynamicEditorIC = new DynamicEditor;
  dynamicEditorIC->setupTabs(elmerDefs, "InitialCondition",1 );
  QMdiSubWindow *initialConditionSubWindow = new QMdiSubWindow;
  initialConditionSubWindow->setWidget(dynamicEditorIC);
  mdiArea->addSubWindow(initialConditionSubWindow);
  initialConditionSubWindow->show();

  if(dynamicEditorBC) delete dynamicEditorBC;
  dynamicEditorBC = new DynamicEditor;
  dynamicEditorBC->setupTabs(elmerDefs, "BoundaryCondition",1 );
  QMdiSubWindow *bcSubWindow = new QMdiSubWindow;
  bcSubWindow->setWidget(dynamicEditorBC);
  mdiArea->addSubWindow(bcSubWindow);
  bcSubWindow->show();

  mdiArea->tileSubWindows();
  //mdiArea->cascadeSubWindows();

}



// Add items from document to tree view...
//----------------------------------------------------------------------------
void EdfEditor::insertItemForElement(QDomElement element,
				     QTreeWidgetItem *parentItem)
{
  if(element.isNull())
    return;

  // set expanded
//  if(parentItem != NULL)
//    parentItem->setExpanded(true);

  // create new tree item
  QTreeWidgetItem *newItem = new QTreeWidgetItem(parentItem);

  newItem->setText(0, element.tagName().trimmed());
  newItem->setFlags(newItem->flags() | Qt::ItemIsEditable);

  // display attributes
  QStringList list;
  QDomNamedNodeMap attributeMap = element.attributes();
  for(int index = 0; index < (int)attributeMap.length(); index++) {
    QDomNode attribute = attributeMap.item(index);
    list << attribute.nodeName() + "=\"" + attribute.nodeValue() + "\"";
  }
  newItem->setText(1, list.join(" "));
  
  // display value
  if(element.firstChildElement().isNull()) 
    newItem->setText(2, element.text().split("\n").join(" ").trimmed());
  
  // update hash
  elementForItem.insert(newItem, element);

  // add item
  edfTree->addTopLevelItem(newItem);
  
  if(!element.firstChildElement().isNull()) 
    insertItemForElement(element.firstChildElement(), newItem);
  
  insertItemForElement(element.nextSiblingElement(), parentItem);      
}

// Construct tree view...
//----------------------------------------------------------------------------
void EdfEditor::setupEditor(QDomDocument *elmerDefs)
{
  this->elmerDefs = elmerDefs;

  disconnect(edfTree, SIGNAL(itemChanged(QTreeWidgetItem*, int)),
	     this, SLOT(updateElement(QTreeWidgetItem*, int)));

  // clear tree view & hash
  edfTree->clear();
  elementForItem.clear(); 

  // get root entry & recursively add items to the tree:
  QDomElement root = elmerDefs->documentElement();
  insertItemForElement(root, NULL);
  edfTree->setCurrentItem(NULL);

  connect(edfTree, SIGNAL(itemChanged(QTreeWidgetItem*, int)),
	  this, SLOT(updateElement(QTreeWidgetItem*, int)));

  edfTree->expandAll();
  expandCollapseAllButton->setText("Collapse all");
  expandCollapseAllButton->setIcon(collapseIcon);
  expandCollapseAll = false;

  edfTree->setCurrentItem(NULL);
}

// Tree view item has been edited: update document accordingly...
//----------------------------------------------------------------------------
void EdfEditor::updateElement(QTreeWidgetItem *item, int column)
{
  // get element from hash
  QDomElement element = elementForItem.value(item);

  if(element.isNull())
    return;

  // set new tag
  element.setTagName(item->text(0).trimmed());

  // delete old attributes
  QDomNamedNodeMap oldAttributes = element.attributes();
  for(int i = 0; i<(int)oldAttributes.length(); i++) {
    QDomNode node = oldAttributes.item(i);
    QString name = node.nodeName();
    element.removeAttribute(name);
  }
  
  // parse and set new attributes
  QString pattern = "([a-zA-Z0-9]+)[ \t]*=[ \t]*[\"]([^\"]+)[\"]";
  QRegExp expression(pattern);
  QString qs = item->text(1).trimmed();

  int index = qs.indexOf(expression);

  QString parsedString = "";
  if(index < 0)
    parsedString = qs;

  while(index >= 0) {
    int length = expression.matchedLength();
    QString currentMatch = qs.mid(index, length);
    QStringList currentList = currentMatch.split("=");
    QString name = currentList.at(0);
    QString value = currentList.at(1);

    int firstPar = value.indexOf("\"", 0);
    int secondPar = value.indexOf("\"", firstPar+1);
    value = value.mid(firstPar+1, secondPar-firstPar-1);

    parsedString.append(name + "=\"" + value + "\" ");

    element.setAttribute(name.trimmed(), value.trimmed());
    index = qs.indexOf(expression, index + length);
  }

  // update display with parsed attributes
  item->setText(1, parsedString);

  // set new text (iff old element has no children)
  if(element.firstChildElement().isNull()) {

    // remove old text node
    QDomNodeList children = element.childNodes();
    for(int i=0;  i<(int)children.length(); i++) {
      QDomNode node = children.at(i);
      if(node.isText())
	element.removeChild(node);
    }
    
    // new text node
    QDomText text = elmerDefs->createTextNode(item->text(2));
    element.appendChild(text);
    
  } else {
    
    // clear value from tree view to avoid confusions:
    item->setText(2, "");
  }

  // no need to update hash
}

// Add tree view item & document element...
//----------------------------------------------------------------------------
void EdfEditor::addButtonClicked()
{
  QTreeWidgetItem *current = edfTree->currentItem();
  
  if(current == NULL)
    return;

  QString newTag = "empty";
  QString newAttribute = "attribute";
  QString newAttributeValue = "empty";
  QString newValue = "empty";
  
  // add item to tree view:
  QTreeWidgetItem *newItem = new QTreeWidgetItem(current);

  newItem->setFlags(newItem->flags() | Qt::ItemIsEditable);

  newItem->setText(0, newTag);
  newItem->setText(1, newAttribute + "=\"" + newAttributeValue + "\"");
  newItem->setText(2, newValue);
  current->addChild(newItem);
  newItem->parent()->setExpanded(true);

  // clear the value field for current item (as it just became parent)
  current->setText(2, "");

  // add to document
  QDomElement newElement = elmerDefs->createElement(newTag);
  newElement.setAttribute(newAttribute, newAttributeValue);

  QDomText newText = elmerDefs->createTextNode(newValue);
  newElement.appendChild(newText);

  QDomElement parent = elementForItem.value(newItem->parent());
  parent.appendChild(newElement);

  // update hash
  elementForItem.insert(newItem, newElement);

  edfTree->setCurrentItem(newItem);
}

// Remove item from tree view item & element from document...
//----------------------------------------------------------------------------
void EdfEditor::removeButtonClicked()
{
  QTreeWidgetItem *currentItem = edfTree->currentItem();

  if(currentItem == NULL)
    return;

  QTreeWidgetItem *parentItem = currentItem->parent();
  QDomElement element = elementForItem.value(currentItem);
  QDomElement parentElement = elementForItem.value(parentItem);

  parentItem->removeChild(currentItem);
  parentElement.removeChild(element);

  // update hash
  elementForItem.remove(currentItem);

  edfTree->setCurrentItem(NULL);
}

// Save as...
//----------------------------------------------------------------------------
void EdfEditor::saveAsButtonClicked()
{
  QString fileName;

  fileName = QFileDialog::getSaveFileName(this,
                 tr("Save definitions"), "", tr("EDF (*.xml)") );

  if(fileName.isEmpty())
    return;

  const int indent = 3;
  
  QFile file;
  file.setFileName(fileName);
  file.open(QIODevice::WriteOnly);
  QTextStream out(&file);
  elmerDefs->save(out, indent);
}


// Expand/collapse tree view...
//----------------------------------------------------------------------------
void EdfEditor::expandCollapseAllButtonClicked()
{
  if(expandCollapseAll) {
    edfTree->expandAll();
    expandCollapseAllButton->setText("Collapse all");
    expandCollapseAllButton->setIcon(collapseIcon);
    expandCollapseAll = false;
  } else {
    edfTree->collapseAll();
    expandCollapseAllButton->setText("Expand all");
    expandCollapseAllButton->setIcon(expandIcon);
    expandCollapseAll = true;
  }
}

// Open...
//----------------------------------------------------------------------------
void EdfEditor::openButtonClicked()
{
  QString fileName;

  fileName = QFileDialog::getOpenFileName(this, 
	      tr("Open definitions"), "", tr("EDF (*.xml)") );

  if(fileName.isEmpty())
    return;

  QFile file;
  file.setFileName(fileName);
  file.open(QIODevice::ReadOnly);

  QString errStr;
  int errRow;
  int errCol;

  if(!elmerDefs->setContent(&file, true, &errStr, &errRow, &errCol)) {
    QMessageBox::information(window(), tr("Elmer definitions file"),
			     tr("Parse error at line %1, col %2:\n%3")
			     .arg(errRow).arg(errCol).arg(errStr));
    file.close();
    return;

  } else {
      
    if(elmerDefs->documentElement().tagName() != "edf") {
      QMessageBox::information(window(), tr("Elmer definitions file"),
			       tr("This is not an edf file"));
      delete elmerDefs;
      file.close();
      return;
      
    }
  }
  
  setupEditor(elmerDefs);

  edfTree->setCurrentItem(NULL);
}


// Append...
//----------------------------------------------------------------------------
void EdfEditor::appendButtonClicked()
{
  QString fileName;

  fileName = QFileDialog::getOpenFileName(this, 
	      tr("Open definitions"), "", tr("EDF (*.xml)") );

  if(fileName.isEmpty())
    return;

  QFile file;
  file.setFileName(fileName);
  file.open(QIODevice::ReadOnly);

  QDomDocument tmpDoc;

  QString errStr;
  int errRow;
  int errCol;

  if(!tmpDoc.setContent(&file, true, &errStr, &errRow, &errCol)) {
    QMessageBox::information(window(), tr("Elmer definitions file"),
			     tr("Parse error at line %1, col %2:\n%3")
			     .arg(errRow).arg(errCol).arg(errStr));
    file.close();
    return;

  } else {
      
    if(tmpDoc.documentElement().tagName() != "edf") {
      QMessageBox::information(window(), tr("Elmer definitions file"),
			       tr("This is not an edf file"));
      file.close();
      return;      
    }
  }

  // add new elements to the document
  QDomElement root = elmerDefs->documentElement();
  QDomElement tmpRoot = tmpDoc.documentElement();

  QDomElement element = tmpRoot.firstChildElement();
  while(!element.isNull()) {
    root.appendChild(element);
    element = tmpRoot.firstChildElement();
  }
  
  setupEditor(elmerDefs);

  edfTree->setCurrentItem(NULL);
}


// Close...
//----------------------------------------------------------------------------
void EdfEditor::applyButtonClicked()
{
  // rebuild document from tree view:
  //---------------------------------
  this->close();
}

// Change the place of two items and elements...
//----------------------------------------------------------------------------
void EdfEditor::treeItemClicked(QTreeWidgetItem *item, int column)
{
  if(item == lastActiveItem)
    return;

  if(lastActiveItem == NULL) {
    lastActiveItem = item;
    return;
  }

  if(!ctrlPressed)
    return;

  // items must have the same parent:
  if(item->parent() != lastActiveItem->parent()) {
    // cout << "Items have different parent - unable to swap" << endl;
    // cout.flush();
    lastActiveItem = item;
    return;
  }

  // get elements:
  QDomElement element = elementForItem.value(item);  
  QDomElement lastActiveItemElement = elementForItem.value(lastActiveItem);  

  // elements must have the same parent (should always be true):
  if(element.parentNode() != lastActiveItemElement.parentNode()) {
    // cout << "Parent element mismatch - unable to swap items" << endl;
    // cout.flush();
    lastActiveItem = item;
    return;
  }
  
  // clone elements:
  QDomNode clone = element.cloneNode(true);
  QDomNode lastActiveItemClone = lastActiveItemElement.cloneNode(true);

  // replace elements with their clones:
  element.parentNode().replaceChild(lastActiveItemClone, element);
  lastActiveItemElement.parentNode().replaceChild(clone, lastActiveItemElement);

  // remove old elements from the document:
  element.parentNode().removeChild(element);
  lastActiveItemElement.parentNode().removeChild(lastActiveItemElement);

  // make sure that old elements are cleared (they should be already):
  element.clear();
  lastActiveItemElement.clear();

  // rebuild tree & hash:
  setupEditor(elmerDefs);

  // set focus back to the last selected item:
  lastActiveItem = NULL;
  for(int i = 0; i < elementForItem.count(); i++) {
    if(elementForItem.values().at(i) == lastActiveItemClone) {
      edfTree->setCurrentItem(elementForItem.keys().at(i));
      edfTree->scrollToItem(elementForItem.keys().at(i), 
			    QAbstractItemView::PositionAtCenter);
    }
  }

  return;
}

// Key pressed...
//-----------------------------------------------------------------------------
void EdfEditor::keyPressEvent(QKeyEvent *event)
{
  if(event->key() == Qt::Key_Control)
    ctrlPressed = true;

  if(event->key() == Qt::Key_Alt)
    altPressed = true;
}


// Key released...
//-----------------------------------------------------------------------------
void EdfEditor::keyReleaseEvent(QKeyEvent *event)
{
  if(event->key() == Qt::Key_Control)
    ctrlPressed = false;

  if(event->key() == Qt::Key_Alt)
    altPressed = false;
}

