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
 *  ElmerGUI timestep                                                        *
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

#ifndef TIMESTEP_H
#define TIMESTEP_H

#include <QWidget>
#include "ui_timestep.h"

class vtkRenderWindow;
class QIcon;

class TimeStep : public QDialog
{
  Q_OBJECT

public:
  TimeStep(QWidget *parent = 0);
  ~TimeStep();

  Ui::timeStepDialog ui;

  int maxSteps;

signals:
  void timeStepChangedSignal();

public slots:
  void SetCurrent(int);                             // set current step
  void SetStart(int);                               // set first step
  void SetStop(int);                                // set last step
  void SetIncrement(int);                           // set increment
  void SetMatcCmd(QString);                         // set MATC cmd
  void RegenerateBeforeDrawing(bool);               // regenerate actors
  void SaveFrames(bool);                            // save frames
  void SetSaveDirectory(QString);                   // set save dir
  void Loop();                                      // toggle looping
  bool IsLooping();                                 // is loop on?
  void DrawCurrent();                               // draw current frame

  void canProceedWithNextSlot(vtkRenderWindow*);

private slots:
  void cancelButtonClicked();
  void applyButtonClicked();
  void okButtonClicked();
  void loopButtonClicked();
  void browseButtonClicked();

private:
  bool loopOn;
  QString saveDir;
};

#endif // TIMESTEP_H
