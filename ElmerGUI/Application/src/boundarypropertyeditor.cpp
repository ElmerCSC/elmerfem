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

#include <QtGui>
#include <iostream>
#include "boundarypropertyeditor.h"

using namespace std;

BoundaryPropertyEditor::BoundaryPropertyEditor(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  touched = false;
  condition = NULL;
  bodyProperties = NULL;

  connect(ui.applyButton, SIGNAL(clicked()), this, SLOT(applySlot()));
  connect(ui.discardButton, SIGNAL(clicked()), this, SLOT(discardSlot()));

  connect(ui.boundaryAsABody, SIGNAL(stateChanged(int)), this, 
           SLOT(boundaryAsABodyChanged(int)));

  connect(ui.boundaryConditionCombo, SIGNAL(currentIndexChanged(QString)), this, 
           SLOT(boundaryComboChanged(QString)));
}

BoundaryPropertyEditor::~BoundaryPropertyEditor()
{
}

void BoundaryPropertyEditor::applySlot()
{
  touched = true;
  this->close();
}

void BoundaryPropertyEditor::discardSlot()
{
  touched = false;
  this->close();
}

void BoundaryPropertyEditor::boundaryComboChanged(QString text)
{
   emit( BoundaryComboChanged(this,text) );
}

void BoundaryPropertyEditor::boundaryAsABodyChanged(int status)
{
   emit( BoundaryAsABodyChanged(this,status) );
}

void BoundaryPropertyEditor::appendToProject(QDomDocument *projectDoc, QDomElement *item)
{
  projectIO.parentWidget = this;
  projectIO.appendToProject(projectDoc, item);
}

void BoundaryPropertyEditor::readFromProject(QDomDocument *projectDoc, QDomElement *item)
{
  projectIO.parentWidget = this;
  projectIO.readFromProject(projectDoc, item);
}
