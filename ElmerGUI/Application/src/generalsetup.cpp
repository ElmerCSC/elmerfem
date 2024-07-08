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
 *  ElmerGUI generalsetup                                                    *
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
#include "generalsetup.h"
#include <QWhatsThis>

using namespace std;

GeneralSetup::GeneralSetup(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  connect(ui.acceptButton, SIGNAL(clicked()), 
	  this, SLOT(acceptButtonClicked()));

  setWindowIcon(QIcon(":/icons/Mesh3D.png"));

  // Set minimum height for free text editors:
  QTextEdit *te = ui.headerFreeTextEdit;
  QFont currentFont = te->currentFont();
  QFontMetrics fontMetrics(currentFont);
  int fontHeight = fontMetrics.height();

  ui.headerFreeTextEdit->setMinimumHeight(3*fontHeight);
  ui.simulationFreeTextEdit->setMinimumHeight(3*fontHeight);
  ui.constantsFreeTextEdit->setMinimumHeight(3*fontHeight);  

  ui.acceptButton->setIcon(QIcon::fromTheme("dialog-accept"));

  whatsThisButton = new QPushButton(tr("Whatis"));
  whatsThisButton->setIcon(QIcon::fromTheme("text-questionmark"));
  connect(whatsThisButton, SIGNAL(clicked()), this, SLOT(whatsThisButtonClicked()));
  whatsThisButton->setWhatsThis("Press this button, then click the widget to be explained.");
  ui.buttonLayout->addWidget(whatsThisButton);   
}

GeneralSetup::~GeneralSetup()
{
}

void GeneralSetup::acceptButtonClicked()
{
  this->close();
}

void GeneralSetup::appendToProject(QDomDocument *projectDoc, QDomElement *item)
{
  projectIO.parentWidget = this;
  projectIO.appendToProject(projectDoc, item);
}

void GeneralSetup::readFromProject(QDomDocument *projectDoc, QDomElement *item)
{
  projectIO.parentWidget = this;
  projectIO.readFromProject(projectDoc, item);
}

void GeneralSetup::whatsThisButtonClicked()
{
  QWhatsThis::enterWhatsThisMode();
}