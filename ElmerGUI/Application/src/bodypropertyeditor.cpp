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
 *  ElmerGUI bodypropertyeditor                                              *
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
#include "mainwindow.h"
#include "bodypropertyeditor.h"

using namespace std;

BodyPropertyEditor::BodyPropertyEditor(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  touched = false;

  material = NULL;
  initial  = NULL;
  force    = NULL;
  equation = NULL;

  connect(ui.applyButton, SIGNAL(clicked()), this, SLOT(applySlot()));

  connect(ui.discardButton, SIGNAL(clicked()), this, SLOT(discardSlot()));

  connect(ui.materialCombo, SIGNAL(currentIndexChanged(QString)), this, 
           SLOT(materialComboChanged(QString)));

  connect(ui.initialConditionCombo, SIGNAL(currentIndexChanged(QString)), this, 
           SLOT(initialComboChanged(QString)));

  connect(ui.bodyForceCombo, SIGNAL(currentIndexChanged(QString)), this, 
           SLOT(forceComboChanged(QString)));

  connect(ui.equationCombo, SIGNAL(currentIndexChanged(QString)), this, 
           SLOT(equationComboChanged(QString)));
}

BodyPropertyEditor::~BodyPropertyEditor()
{
}

void BodyPropertyEditor::applySlot()
{
  touched = true;
  this->close();
}

void BodyPropertyEditor::discardSlot()
{
  touched = false;
  this->close();
}

void BodyPropertyEditor::materialComboChanged(QString text)
{
   emit( BodyMaterialComboChanged(this,text) );
}

void BodyPropertyEditor::initialComboChanged(QString text)
{
   emit( BodyInitialComboChanged(this,text) );
}

void BodyPropertyEditor::forceComboChanged(QString text)
{
   emit( BodyForceComboChanged(this,text) );
}

void BodyPropertyEditor::equationComboChanged(QString text)
{
   emit( BodyEquationComboChanged(this,text) );
}

void BodyPropertyEditor::appendToProject(QDomDocument *projectDoc, QDomElement *item)
{
  projectIO.parentWidget = this;
  projectIO.appendToProject(projectDoc, item);
}

void BodyPropertyEditor::readFromProject(QDomDocument *projectDoc, QDomElement *item)
{
  projectIO.parentWidget = this;
  projectIO.readFromProject(projectDoc, item);
}
