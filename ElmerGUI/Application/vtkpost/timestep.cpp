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

#include <QtGui>
#include <iostream>
#include <vtkRenderWindow.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include "timestep.h"

using namespace std;

TimeStep::TimeStep(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  connect(ui.loopButton, SIGNAL(clicked()), this, SLOT(loopButtonClicked()));
  connect(ui.cancelButton, SIGNAL(clicked()), this, SLOT(cancelButtonClicked()));
  connect(ui.applyButton, SIGNAL(clicked()), this, SLOT(applyButtonClicked()));
  connect(ui.okButton, SIGNAL(clicked()), this, SLOT(okButtonClicked()));
  connect(ui.browseButton, SIGNAL(clicked()), this, SLOT(browseButtonClicked()));

  maxSteps = 0;
  loopOn = false;
  saveDir = "";

  setWindowTitle("Time step control");
  setWindowIcon(QIcon(":/icons/Mesh3D.png"));

  ui.applyButton->setIcon(QIcon(":/icons/dialog-ok.png"));
  ui.loopButton->setIcon(QIcon(":/icons/dialog-ok.png"));
}

TimeStep::~TimeStep()
{
}

void TimeStep::browseButtonClicked()
{
  QString saveDir = QFileDialog::getExistingDirectory(this,
	 tr("Save directory"), "", QFileDialog::ShowDirsOnly
	       | QFileDialog::DontResolveSymlinks);

  ui.saveDirectory->setText(saveDir);
}

void TimeStep::cancelButtonClicked()
{
  close();
}

void TimeStep::okButtonClicked()
{
  applyButtonClicked();
  cancelButtonClicked();
}

void TimeStep::applyButtonClicked()
{
  int current = ui.timeStep->value();

  if(current < 1) {
    current = 1;
    ui.timeStep->setValue(current);
  }

  if(current > maxSteps) {
    current = maxSteps;
    ui.timeStep->setValue(current);
  }

  emit(timeStepChangedSignal());
}

void TimeStep::canProceedWithNextSlot(vtkRenderWindow *renderWindow)
{
  if(!loopOn) return;

  bool saveFrames = ui.saveFrames->isChecked();
  int current = ui.timeStep->value();
  int stop = ui.stop->value();
  int increment = ui.increment->value();

  if(saveFrames) {
    QString saveDir = ui.saveDirectory->text().trimmed();
    if(saveDir.isEmpty()) saveDir = ".";
    QString frameName = "frame" + QString::number(current) + ".png";
    QString fileName = saveDir + "/" + frameName;

    vtkWindowToImageFilter *image = vtkWindowToImageFilter::New();
    image->SetInput(renderWindow);
    image->Update();
    
    vtkPNGWriter *writer =  vtkPNGWriter::New();
    writer->SetInputConnection(image->GetOutputPort());
    writer->SetFileName(fileName.toAscii().data());

    renderWindow->Render();
    writer->Write();

    image->Delete();
    writer->Delete();
  }

  if(increment < 1) {
    increment = 1;
    ui.increment->setValue(increment);
  }

  if(stop > maxSteps) {
    stop = maxSteps;
    ui.stop->setValue(stop);
  }

  if(current > stop) {
    loopOn = false;
    ui.loopButton->setText("Loop");
    ui.loopButton->setIcon(QIcon(":/icons/dialog-ok.png"));
    this->repaint();

  } else {
    ui.loopButton->setText("Stop");
    ui.loopButton->setIcon(QIcon(":/icons/dialog-close.png"));
    this->repaint();
    current += increment;

    if(current > stop) {
      loopOn = false;
      ui.loopButton->setText("Loop");
      ui.loopButton->setIcon(QIcon(":/icons/dialog-ok.png"));
      this->repaint();
      return;
    }

    ui.timeStep->setValue(current);
    this->repaint();
    applyButtonClicked();
  }
}

void TimeStep::loopButtonClicked()
{
  if(loopOn) {
    loopOn = false;
    ui.loopButton->setText("Loop");
    ui.loopButton->setIcon(QIcon(":/icons/dialog-ok.png"));
    this->repaint();
 
  } else {
    loopOn = true;
    ui.loopButton->setText("Stop");
    ui.loopButton->setIcon(QIcon(":/icons/dialog-close.png"));
    this->repaint();
    int start = ui.start->value();
    int stop = ui.stop->value();
    int increment = ui.increment->value();

    if(start < 1) {
      start = 1;
      ui.start->setValue(start);
    }

    if(start > maxSteps) {
      start = maxSteps;
      ui.start->setValue(start);
    }

    if(stop < 1) {
      stop = 1;
      ui.stop->setValue(stop);
    }

    if(stop > maxSteps) {
      stop = maxSteps;
      ui.stop->setValue(stop);
    }

    if(stop < start) {
      stop = start;
      ui.stop->setValue(stop);
    }

    if(increment < 1) {
      increment = 1;
      ui.increment->setValue(increment);
    }

    ui.timeStep->setValue(start);
    this->repaint();
    applyButtonClicked();
  }
}

void TimeStep::SetCurrent(int n)
{
  ui.timeStep->setValue(n);
}

void TimeStep::SetStart(int n)
{
  ui.start->setValue(n);
}

void TimeStep::SetStop(int n)
{
  ui.stop->setValue(n);
}

void TimeStep::SetIncrement(int n)
{
  ui.increment->setValue(n);
}

void TimeStep::SetMatcCmd(QString cmd)
{
  ui.doBefore->setText(cmd);
}

void TimeStep::RegenerateBeforeDrawing(bool b)
{
  ui.regenerateBeforeDrawing->setChecked(b);
}

void TimeStep::SaveFrames(bool b)
{
  ui.saveFrames->setChecked(b);
}

void TimeStep::SetSaveDirectory(QString dirName)
{
  ui.saveDirectory->setText(dirName);
}

void TimeStep::Loop()
{
  loopButtonClicked();
}

bool TimeStep::IsLooping()
{
  return this->loopOn;
}

void TimeStep::DrawCurrent()
{
  if(loopOn) return;
  this->applyButtonClicked();
}
