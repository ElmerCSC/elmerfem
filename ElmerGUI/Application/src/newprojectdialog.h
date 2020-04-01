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
 *  ElmerGUI newproject                                                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Author: Saeki Takayuki                                                   *
 *  Original Date: 15 Feb 2020                                               *
 *                                                                           *
 *****************************************************************************/

#ifndef NEWPROJECTDIALOG_H
#define NEWPROJECTDIALOG_H

#include <QWidget>
#include "ui_newproject.h"

class NewProjectDialog : public QDialog
{
  Q_OBJECT

public:
  NewProjectDialog(QWidget *parent = 0);
  ~NewProjectDialog();

  Ui::newProjectDialog ui;

  void setDirectories(QString& deaultDir, QString& extraDir);
  
signals:

private slots:
  void elmerMeshToggled(bool);
  void geometryFileToggled(bool);
  void laterToggled(bool);
  void projectDirClicked(bool);
  void meshDirClicked(bool);
  void geometryFileClicked(bool);
  void addSolverClicked(bool);
  void removeSolverClicked(bool);
  void selectedSolverChanged(int);
  void unselectedSolverChanged(int);
    
private:
  QString defaultDirName;
  QString extraDirPath;
};

#endif // NEWPROJECTDIALOG_H
