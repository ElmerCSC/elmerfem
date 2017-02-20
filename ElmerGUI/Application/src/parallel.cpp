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
 *  ElmerGUI parallel control                                                *
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
#include "parallel.h"

#if WITH_QT5
#include <QtWidgets>
#endif

using namespace std;

Parallel::Parallel(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  connect(ui.browseButton, SIGNAL(clicked()), 
	  this, SLOT(browseButtonClicked()));

  connect(ui.defaultsButton, SIGNAL(clicked()), 
	  this, SLOT(defaultsButtonClicked()));

  connect(ui.okButton, SIGNAL(clicked()), 
	  this, SLOT(okButtonClicked()));

  defaultsButtonClicked();

  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
}

Parallel::~Parallel()
{
}

void Parallel::browseButtonClicked()
{
  QString fileName = QFileDialog::getOpenFileName(this);

  if(fileName.isEmpty())
    return;

  ui.parallelExecLineEdit->setText(fileName);
}

void Parallel::defaultsButtonClicked()
{
  ui.parallelActiveCheckBox->setChecked(false);
  ui.skipPartitioningCheckBox->setChecked(false);
  ui.nofProcessorsSpinBox->setValue(2);
  
#ifdef WIN32
  ui.parallelExecLineEdit->setText("C:/Program Files/MPICH2/bin/mpiexec.exe");
  ui.parallelArgsLineEdit->setText("-localonly %n -genvlist PATH,ELMER_HOME ElmerSolver_mpi.exe");
#else
  ui.parallelExecLineEdit->setText("mpirun");
  ui.parallelArgsLineEdit->setText("-np %n ElmerSolver_mpi");
#endif

  ui.divideLineEdit->setText("ElmerGrid 2 2 %msh -metis %n");
  ui.mergeLineEdit->setText("ElmerGrid 15 3 %ep -partjoin %n");
}

void Parallel::okButtonClicked()
{
  this->close();
}

void Parallel::appendToProject(QDomDocument *projectDoc, QDomElement *item)
{
  projectIO.parentWidget = this;
  projectIO.appendToProject(projectDoc, item);
}

void Parallel::readFromProject(QDomDocument *projectDoc, QDomElement *item)
{
  projectIO.parentWidget = this;
  projectIO.readFromProject(projectDoc, item);
}
