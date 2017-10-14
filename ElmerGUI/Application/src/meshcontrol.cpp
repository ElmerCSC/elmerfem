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
 *  ElmerGUI meshcontrol                                                     *
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
#include "meshcontrol.h"

#include <stdio.h>

using namespace std;

MeshControl::MeshControl(QWidget *parent)
  : QDialog(parent)
{
  tetlibPresent = true;
  nglibPresent = true;

  ui.setupUi(this);

  connect(ui.tetlibRadioButton, SIGNAL(clicked()), this, SLOT(tetlibClicked()));
  connect(ui.nglibRadioButton, SIGNAL(clicked()), this, SLOT(nglibClicked()));
  connect(ui.elmerGridRadioButton, SIGNAL(clicked()), this, SLOT(elmerGridClicked()));

  connect(ui.tetlibStringEdit, SIGNAL(textChanged(const QString&)), this, SLOT(defineTetlibControlString(const QString&)));

  connect(ui.nglibMaxHEdit, SIGNAL(textChanged(const QString&)), this, SLOT(defineNglibMaxH(const QString&)));
  connect(ui.nglibFinenessEdit, SIGNAL(textChanged(const QString&)), this, SLOT(defineNglibFineness(const QString&)));
  connect(ui.nglibBgmeshEdit, SIGNAL(textChanged(const QString&)), this, SLOT(defineNglibBackgroundmesh(const QString&)));

  connect(ui.defaultsButton, SIGNAL(clicked()), this, SLOT(defaultControls()));
  connect(ui.closeButton, SIGNAL(clicked()), this, SLOT(close()));

  connect(ui.elmerGridStringEdit, SIGNAL(textChanged(const QString&)), this, SLOT(defineElmerGridControlString(const QString&)));

  connect(ui.elmerGridStringEdit, SIGNAL(textChanged(const QString&)), this, SLOT(defineElmerGridControlString(const QString&)));

  connect(ui.elementCodesStringEdit, SIGNAL(textChanged(const QString&)), this, SLOT(defineElementCodesString(const QString&)));

  //  defaultControls(); // Note - already called from mainwindow.cpp
}

MeshControl::~MeshControl()
{
}

void MeshControl::tetlibClicked()
{
  generatorType = GEN_TETLIB;
}

void MeshControl::nglibClicked()
{
  generatorType = GEN_NGLIB;
}

void MeshControl::elmerGridClicked()
{
  generatorType = GEN_ELMERGRID;
}

void MeshControl::defineElementCodesString(const QString &qs)
{
  elementCodesString = qs;
}

void MeshControl::defineTetlibControlString(const QString &qs)
{
  tetlibControlString = qs;
}

void MeshControl::defineNglibMaxH(const QString &qs)
{
  nglibMaxH = qs;
}

void MeshControl::defineNglibFineness(const QString &qs)
{
  nglibFineness = qs;
}

void MeshControl::defineNglibBackgroundmesh(const QString &qs)
{
  nglibBackgroundmesh = qs;
}

void MeshControl::defineElmerGridControlString(const QString &qs)
{
  elmerGridControlString = qs;
}

void MeshControl::defaultControls()
{
  // generatorType = GEN_TETLIB;
  // ui.tetlibRadioButton->setChecked(true);

  generatorType = GEN_NGLIB;
  ui.nglibRadioButton->setChecked(true);

  // if(!tetlibPresent) {
  // generatorType = GEN_NGLIB;
  // ui.nglibRadioButton->setChecked(true);
  // }

  // if(!tetlibPresent && !nglibPresent) {
  // generatorType = GEN_ELMERGRID;
  // ui.elmerGridRadioButton->setChecked(true);
  // }

  // Check for valid tetgen quality in cmd line arguments:
  //---------------------------
  QStringList args = QCoreApplication::arguments();
  int tq_ind = args.indexOf("-tq");
  bool tq_isnum=0;
  QString tetlib_options = "nnJApVq";
  if (tq_ind > 0) args.at(tq_ind + 1).toFloat(&tq_isnum);
#if WITH_QT5
  if (tq_ind > 0 && !tq_isnum) cout << "Ignoring -tq option: " << args.at(tq_ind + 1).toLatin1().data() << endl;
#else
  if (tq_ind > 0 && !tq_isnum) cout << "Ignoring -tq option: " << args.at(tq_ind + 1).toAscii().data() << endl;
#endif
  if (tq_isnum) tetlib_options.append(args.at(tq_ind + 1));
  else tetlib_options.append("1.414");
 
  ui.tetlibStringEdit->setText(tetlib_options);
  ui.nglibMaxHEdit->setText("1000000");
  ui.nglibFinenessEdit->setText("0.5");
  ui.nglibBgmeshEdit->setText("");
  ui.elmerGridStringEdit->setText("-relh 1.0");
  ui.elementCodesStringEdit->setText("");
}

void MeshControl::appendToProject(QDomDocument *projectDoc, QDomElement *item)
{
  projectIO.parentWidget = this;
  projectIO.appendToProject(projectDoc, item);
}

void MeshControl::readFromProject(QDomDocument *projectDoc, QDomElement *item)
{
  projectIO.parentWidget = this;
  projectIO.readFromProject(projectDoc, item);
}
