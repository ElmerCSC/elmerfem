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
 *  ElmerGUI newproject                                                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Author: Saeki Takayuki                                                   *
 *  Original Date: 15 Feb 2020                                               *
 *                                                                           *
 *****************************************************************************/

#include <QtGui>
#include <QFileDialog>
#include <iostream>
#include "newprojectdialog.h"

using namespace std;

NewProjectDialog::NewProjectDialog(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
   
  connect(ui.buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()), this, SLOT(reject()));
  connect(ui.radioButton_elmerMesh, SIGNAL(toggled(bool)), this, SLOT(elmerMeshToggled(bool)));
  connect(ui.radioButton_geometryFile, SIGNAL(toggled(bool)), this, SLOT(geometryFileToggled(bool)));  
  connect(ui.radioButton_later, SIGNAL(toggled(bool)), this, SLOT(laterToggled(bool)));
  connect(ui.pushButton_projectDir, SIGNAL(clicked(bool)), this, SLOT(projectDirClicked(bool)));
  connect(ui.pushButton_meshDir, SIGNAL(clicked(bool)), this, SLOT(meshDirClicked(bool)));
  connect(ui.pushButton_geometryFile, SIGNAL(clicked(bool)), this, SLOT(geometryFileClicked(bool)));
  connect(ui.pushButton_addSolver, SIGNAL(clicked(bool)), this, SLOT(addSolverClicked(bool)));
  connect(ui.pushButton_removeSolver, SIGNAL(clicked(bool)), this, SLOT(removeSolverClicked(bool)));
  connect(ui.listWidget_selectedSolvers, SIGNAL(currentRowChanged(int)), this, SLOT(selectedSolverChanged(int)));
  connect(ui.listWidget_unselectedSolvers, SIGNAL(currentRowChanged(int)), this, SLOT(unselectedSolverChanged(int)));
  ui.pushButton_addSolver->setEnabled(false);
  ui.pushButton_removeSolver->setEnabled(false);

  ui.radioButton_later->setChecked(true);
  geometryFileToggled(false);
  elmerMeshToggled(false);
  
  ui.buttonBox->button(QDialogButtonBox::Ok)->setEnabled(false);
}

NewProjectDialog::~NewProjectDialog()
{
}

void NewProjectDialog::elmerMeshToggled(bool b){
  ui.pushButton_meshDir->setEnabled(b);
  ui.label_meshDir->setEnabled(b);  
}

void NewProjectDialog::geometryFileToggled(bool b){
  ui.pushButton_geometryFile->setEnabled(b);
  ui.label_geometryFile->setEnabled(b);  
}

void NewProjectDialog::laterToggled(bool b){
}

void NewProjectDialog::projectDirClicked(bool b){
  QString dirName = QFileDialog::getExistingDirectory(this, tr("Select directory to store the new project"), defaultDirName);
  if (!dirName.isEmpty()){
    ui.label_projectDir->setText(dirName);
    ui.buttonBox->button(QDialogButtonBox::Ok)->setEnabled(true);
  }
}

void NewProjectDialog::meshDirClicked(bool b){
  QString dirName = QFileDialog::getExistingDirectory(this, tr("Select Elmer mesh directory"), defaultDirName);
  if (!dirName.isEmpty()){
    ui.label_meshDir->setText(dirName);
  }
}

void NewProjectDialog::geometryFileClicked(bool b){
  QString fileName = QFileDialog::getOpenFileName(this, tr("Select geometry input file"), defaultDirName);
  if (!fileName.isEmpty()) {
     ui.label_geometryFile->setText(fileName);
  }
}

void NewProjectDialog::setDirectories(QString& defaultDir, QString& extraDirName){
  defaultDirName = defaultDir;
  extraDirPath = extraDirName;

  QString name;
  QDir extraDir(extraDirPath);
  QStringList nameFilters;
  nameFilters << "*.xml";
  QStringList fileNameList = extraDir.entryList(nameFilters, QDir::Files | QDir::Readable);
  
  ui.listWidget_unselectedSolvers->addItems(fileNameList);
  
  QString labelString;
  QDir edfDir(extraDirPath + "/../edf");
  fileNameList = edfDir.entryList(nameFilters, QDir::Files | QDir::Readable);
  for(int i=0; i<fileNameList.size(); i++){
    if( fileNameList[i] != "edf.xml" && fileNameList[i] != "egini.xml" && fileNameList[i] != "egmaterials.xml")
    {
      labelString +=  " " + fileNameList[i] + "\n";
    }
  }
  ui.label_defaultSolvers->setText(labelString);
}

void NewProjectDialog::addSolverClicked(bool b){
  int i = ui.listWidget_unselectedSolvers->currentRow();
  if(i < 0) return;
  QListWidgetItem* item = ui.listWidget_unselectedSolvers->takeItem(i);
  if(item != NULL) ui.listWidget_selectedSolvers->addItem(item->text());
}

void NewProjectDialog::removeSolverClicked(bool b){
  int i = ui.listWidget_selectedSolvers->currentRow();
  if(i < 0) return;
  QListWidgetItem* item = ui.listWidget_selectedSolvers->takeItem(i);
  if(item != NULL) ui.listWidget_unselectedSolvers->addItem(item->text());
}

void NewProjectDialog::selectedSolverChanged(int i){
  ui.pushButton_removeSolver->setEnabled(i >= 0);
}

void NewProjectDialog::unselectedSolverChanged(int i){
  ui.pushButton_addSolver->setEnabled(i >= 0); 
}