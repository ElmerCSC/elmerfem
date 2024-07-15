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
 *  ElmerGUI solverparameters                                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Authors: Mikko Lyly, Juha Ruokolainen and Peter R�back                   *
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
#include "solverparameters.h"

using namespace std;

SolverParameterEditor::SolverParameterEditor(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);
  generalOptions = NULL;
  solverName = "";

  this->setWindowIcon(QIcon(":/icons/Mesh3D.png"));

  connect(ui.applyButton, SIGNAL(clicked()), this, SLOT(close()));

  connect(ui.useHypre, SIGNAL(stateChanged(int)), this, SLOT(hypreStateChanged(int)));
  connect(ui.useParasails, SIGNAL(stateChanged(int)), this, SLOT(parasailsStateChanged(int)));
  connect(ui.useBoomerAMG, SIGNAL(stateChanged(int)), this, SLOT(boomerAMGStateChanged(int)));

  hypreStateChanged(0);
  
  ui.applyButton->setIcon(QIcon::fromTheme("dialog-accept"));

  whatsThisButton = new QPushButton(tr("Whatis"));
  whatsThisButton->setIcon(QIcon::fromTheme("text-questionmark"));
  connect(whatsThisButton, SIGNAL(clicked()), this, SLOT(whatsThisButtonClicked()));
  whatsThisButton->setWhatsThis("Press this button, then click the widget to be explained.");
  ui.buttonLayout->addWidget(whatsThisButton);  
}

SolverParameterEditor::~SolverParameterEditor()
{
}

void SolverParameterEditor::appendToProject(QDomDocument *projectDoc, QDomElement *item)
{
  projectIO.parentWidget = this;
  projectIO.appendToProject(projectDoc, item);
}

void SolverParameterEditor::readFromProject(QDomDocument *projectDoc, QDomElement *item)
{
  projectIO.parentWidget = this;
  projectIO.readFromProject(projectDoc, item);
}

void SolverParameterEditor::hypreStateChanged(int)
{
  if(ui.useHypre->isChecked()) {
    ui.parasailsGroup->setEnabled(true);
    ui.boomerAMGGroup->setEnabled(true);
  } else {
    ui.parasailsGroup->setEnabled(false);
    ui.boomerAMGGroup->setEnabled(false);
  }

}

void SolverParameterEditor::parasailsStateChanged(int)
{
  if(ui.useParasails->isChecked()) {
    ui.boomerAMGGroup->setEnabled(false);
  } else {
    ui.boomerAMGGroup->setEnabled(true);
  }
}

void SolverParameterEditor::boomerAMGStateChanged(int)
{
  if(ui.useBoomerAMG->isChecked()) {
    ui.parasailsGroup->setEnabled(false);
  } else {
    ui.parasailsGroup->setEnabled(true);
  }
}

void SolverParameterEditor::whatsThisButtonClicked()
{
  QWhatsThis::enterWhatsThisMode();
}
