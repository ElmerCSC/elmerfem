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
 *  ElmerGUIlogger                                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Authors: Mikko Lyly, Juha Ruokolainen and Peter Råback                   *
 *  Small fixes: Boris Pek                                                   *
 *  Email:   Juha.Ruokolainen@csc.fi                                         *
 *  Web:     http://www.csc.fi/elmer                                         *
 *  Address: CSC - IT Center for Science Ltd.                                 *
 *           Keilaranta 14                                                   *
 *           02101 Espoo, Finland                                            *
 *                                                                           *
 *  Original Date: 24 Aug 2009                                               *
 *                                                                           *
 *****************************************************************************/

#ifndef MAINWINDOW_H
#define MAINWONDOW_H

#include <QMainWindow>
#include <QProcess>

class QTextEdit;
class QMenu;
class QAction;

class MainWindow : public QMainWindow 
{
  Q_OBJECT

 public:
  MainWindow();
  ~MainWindow();

 private slots:
  void startElmerGUISlot();
  void saveAsSlot();
  void printSlot();
  void exitSlot();
  void errorSlot(QProcess::ProcessError);
  void finishedSlot(int, QProcess::ExitStatus);
  void stdoutSlot();
  void stderrSlot();
  void startedSlot();
  void stateChangedSlot(QProcess::ProcessState);

 private:
  void createActions();
  void createMenus();
  
  QMenu* fileMenu;
  QAction* startElmerGUIAct;
  QAction* saveAsAct;
  QAction* printAct;
  QAction* exitAct;
  QTextEdit* textEdit;
  QProcess* elmerGUI;

};

#endif // MAINWINDOW_H
