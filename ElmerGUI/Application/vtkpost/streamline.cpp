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
 *  ElmerGUI streamline                                                      *
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
#include "epmesh.h"
#include "vtkpost.h"
#include "streamline.h"
#include "timestep.h"

#include <vtkUnstructuredGrid.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPointSource.h>
#include <vtkLineSource.h>
#include <vtkStreamLine.h>
#include <vtkRungeKutta4.h>
#include <vtkRibbonFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkLookupTable.h>
#include <vtkActor.h>

using namespace std;

StreamLine::StreamLine(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  connect(ui.cancelButton, SIGNAL(clicked()), this, SLOT(cancelButtonClicked()));
  connect(ui.applyButton, SIGNAL(clicked()), this, SLOT(applyButtonClicked()));
  connect(ui.okButton, SIGNAL(clicked()), this, SLOT(okButtonClicked()));
  connect(ui.colorCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(colorSelectionChanged(int)));
  connect(ui.keepLimits, SIGNAL(stateChanged(int)), this, SLOT(keepLimitsSlot(int)));

  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
}

StreamLine::~StreamLine()
{
}

void StreamLine::applyButtonClicked()
{
  emit(drawStreamLineSignal());
}

void StreamLine::cancelButtonClicked()
{
  emit(hideStreamLineSignal());
  close();
}

void StreamLine::okButtonClicked()
{
  applyButtonClicked();
  close();
}

void StreamLine::populateWidgets(VtkPost* vtkPost)
{
  this->scalarField = vtkPost->GetScalarField();
  this->scalarFields = vtkPost->GetScalarFields();

  ui.vectorCombo->clear();

  int index = -1;
  for(int i = 0; i < scalarFields; i++) {
    ScalarField *sf = &scalarField[i];
    QString name = sf->name;
    if((index = name.indexOf("_x")) >= 0) {
      ui.vectorCombo->addItem(name.mid(0, index));
    }
  }

  QString name = ui.colorCombo->currentText();

  ui.colorCombo->clear();

  for(int i = 0; i < scalarFields; i++) {
    ScalarField *sf = &scalarField[i];
    ui.colorCombo->addItem(sf->name);
  }

  this->SetColorName(name);

  colorSelectionChanged(ui.colorCombo->currentIndex());
}

void StreamLine::colorSelectionChanged(int newIndex)
{
  ScalarField *sf = &this->scalarField[newIndex];
  if(!ui.keepLimits->isChecked()) {
    ui.minVal->setText(QString::number(sf->minVal));
    ui.maxVal->setText(QString::number(sf->maxVal));
  }
}

void StreamLine::keepLimitsSlot(int state)
{
  if(state == 0)
    colorSelectionChanged(ui.colorCombo->currentIndex());
}

void StreamLine::draw(VtkPost* vtkPost, TimeStep* timeStep)
{

  QString vectorName = ui.vectorCombo->currentText();

  if(vectorName.isEmpty()) return;

  int i, j, index = -1;
  for(i = 0; i < scalarFields; i++) {
    ScalarField* sf = &scalarField[i];
    QString name = sf->name;
    if((j = name.indexOf("_x")) >= 0) {
      if(vectorName == name.mid(0, j)) {
	index = i;
	break;
      }
    }
  }

  if(index < 0) return;

  // Controls:
  double propagationTime = ui.propagationTime->text().toDouble();
  double stepLength = ui.stepLength->text().toDouble();
  double integStepLength = ui.integStepLength->text().toDouble();
  int threads = ui.threads->value();
  bool useSurfaceGrid = ui.useSurfaceGrid->isChecked();
  bool lineSource = ui.lineSource->isChecked();
  bool sphereSource = ui.sphereSource->isChecked();
  bool pickSource = ui.pickSource->isChecked();
  bool forward = ui.forward->isChecked();
  bool backward = ui.backward->isChecked();

  if(!(forward || backward)) return;

  // Color:
  int colorIndex = ui.colorCombo->currentIndex();
  QString colorName = ui.colorCombo->currentText();
  double minVal = ui.minVal->text().toDouble();
  double maxVal = ui.maxVal->text().toDouble();

  // Appearance:
  bool drawRibbon = ui.drawRibbon->isChecked();
  int ribbonWidth = ui.ribbonWidth->value();
  int lineWidth = ui.lineWidth->text().toInt();

  // Line source:
  double startX = ui.startX->text().toDouble();
  double startY = ui.startY->text().toDouble();
  double startZ = ui.startZ->text().toDouble();
  double endX = ui.endX->text().toDouble();
  double endY = ui.endY->text().toDouble();
  double endZ = ui.endZ->text().toDouble();
  int lines = ui.lines->value();
  bool drawRake = ui.rake->isChecked();
  int rakeWidth = ui.rakeWidth->value();

  ScalarField* sf_x = &scalarField[index + 0];
  ScalarField* sf_y = &scalarField[index + 1];
  ScalarField* sf_z = &scalarField[index + 2];
  int maxDataStepVector = sf_x->values / vtkPost->NofNodes();
  int step = timeStep->ui.timeStep->value();
  if(step > maxDataStepVector) step = maxDataStepVector;
  if(step > timeStep->maxSteps) step = timeStep->maxSteps;
  int vectorOffset = vtkPost->NofNodes() * (step-1);

  ScalarField* sf = &scalarField[colorIndex];
  int maxDataStepColor = sf->values / vtkPost->NofNodes();
  step = timeStep->ui.timeStep->value();
  if(step > maxDataStepColor) step = maxDataStepColor;
  if(step > timeStep->maxSteps) step = timeStep->maxSteps;
  int colorOffset = vtkPost->NofNodes() * (step-1);

  // Sphere source:
  double centerX = ui.centerX->text().toDouble();
  double centerY = ui.centerY->text().toDouble();
  double centerZ = ui.centerZ->text().toDouble();
  double radius = ui.radius->text().toDouble();
  int points = ui.points->value();

  // Pick source:
  double* currentPickPosition = vtkPost->GetCurrentPickPosition();
  double pickX = currentPickPosition[0];
  double pickY = currentPickPosition[1];
  double pickZ = currentPickPosition[2];

  // Choose the grid:
  //------------------
  vtkUnstructuredGrid* grid = NULL;
  if(useSurfaceGrid)
    grid = vtkPost->GetSurfaceGrid();
  else
    grid = vtkPost->GetVolumeGrid();

  if(!grid) return;

  if(grid->GetNumberOfCells() < 1) return;

  // Vector data:
  //-------------
  grid->GetPointData()->RemoveArray("VectorData");
  vtkFloatArray* vectorData = vtkFloatArray::New();
  vectorData->SetNumberOfComponents(3);
  vectorData->SetNumberOfTuples(vtkPost->NofNodes());
  vectorData->SetName("VectorData");
  for(int i = 0; i < vtkPost->NofNodes(); i++) {
    double val_x  = sf_x->value[i + vectorOffset];
    double val_y  = sf_y->value[i + vectorOffset];
    double val_z  = sf_z->value[i + vectorOffset];
    vectorData->SetComponent(i,0,val_x); 
    vectorData->SetComponent(i,1,val_y); 
    vectorData->SetComponent(i,2,val_z); 
  }
  grid->GetPointData()->AddArray(vectorData);

  // Color data:
  //-------------
  grid->GetPointData()->RemoveArray("StreamLineColor");
  vtkFloatArray* vectorColor = vtkFloatArray::New();
  vectorColor->SetNumberOfComponents(1);
  vectorColor->SetNumberOfTuples(vtkPost->NofNodes());
  vectorColor->SetName("StreamLineColor");
  for(int i = 0; i < vtkPost->NofNodes(); i++) 
    vectorColor->SetComponent(i, 0, sf->value[i + colorOffset]); 
  grid->GetPointData()->AddArray(vectorColor);

  // Generate stream lines:
  //-----------------------
  grid->GetPointData()->SetActiveVectors("VectorData");
  grid->GetPointData()->SetActiveScalars("StreamLineColor");
  vtkPointSource* point = vtkPointSource::New();
  vtkLineSource* line = vtkLineSource::New();
  if(lineSource) {
    line->SetPoint1(startX, startY, startZ);
    line->SetPoint2(endX, endY, endZ);
    line->SetResolution(lines);
  } else {
    if(sphereSource) {
      point->SetCenter(centerX, centerY, centerZ);
      point->SetRadius(radius);
      point->SetNumberOfPoints(points);
      point->SetDistributionToUniform();
    } else {
      point->SetCenter(pickX, pickY, pickZ);
      point->SetRadius(0.0);
      point->SetNumberOfPoints(1);
    }
  }

  vtkStreamLine* streamer = vtkStreamLine::New();
  vtkRungeKutta4* integrator = vtkRungeKutta4::New();

  streamer->SetInput(grid);

  if(lineSource) {
    streamer->SetSource(line->GetOutput());
  } else {
    streamer->SetSource(point->GetOutput());
  }

  streamer->SetIntegrator(integrator);
  streamer->SetMaximumPropagationTime(propagationTime);
  streamer->SetIntegrationStepLength(integStepLength);
  if(forward && backward) {
    streamer->SetIntegrationDirectionToIntegrateBothDirections();
  } else if(forward) {
    streamer->SetIntegrationDirectionToForward();    
  } else {
    streamer->SetIntegrationDirectionToBackward();
  }
  streamer->SetStepLength(stepLength);
  streamer->SetNumberOfThreads(threads);
  
  vtkRibbonFilter* ribbon = vtkRibbonFilter::New();

  if(drawRibbon) {
    double length = grid->GetLength();
    ribbon->SetInputConnection(streamer->GetOutputPort());
    ribbon->SetWidth(ribbonWidth * length / 1000.0);
    ribbon->SetWidthFactor(5);
  }

  vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();

  mapper->ScalarVisibilityOn();
  mapper->SetScalarRange(minVal, maxVal);

  if(drawRibbon) {
    mapper->SetInputConnection(ribbon->GetOutputPort());
  } else {
    mapper->SetInputConnection(streamer->GetOutputPort());
  }

  mapper->SetColorModeToMapScalars();
  mapper->SetLookupTable(vtkPost->GetCurrentLut());

  vtkPost->GetStreamLineActor()->SetMapper(mapper);

  if(!drawRibbon)
    vtkPost->GetStreamLineActor()->GetProperty()->SetLineWidth(lineWidth);

  vtkPost->SetCurrentStreamLineName(colorName);

  // Clean up:
  //----------
  line->Delete();
  point->Delete();
  vectorData->Delete();
  vectorColor->Delete();
  integrator->Delete();
  streamer->Delete();
  ribbon->Delete();
  mapper->Delete();
}

// Public slots:
//---------------
QString StreamLine::GetFieldName()
{
  return ui.vectorCombo->currentText();
}

QString StreamLine::GetColorName()
{
  return ui.colorCombo->currentText();
}

bool StreamLine::SetFieldName(QString name)
{
  for(int i = 0; i < ui.vectorCombo->count(); i++) {
    if(ui.vectorCombo->itemText(i) == name) {
      ui.vectorCombo->setCurrentIndex(i);
      return true;
    }
  }
  return false;
}

bool StreamLine::SetColorName(QString name)
{
  for(int i = 0; i < ui.colorCombo->count(); i++) {
    if(ui.colorCombo->itemText(i) == name) {
      ui.colorCombo->setCurrentIndex(i);
      return true;
    }
  }
  return false;
}

void StreamLine::SetMaxTime(double f)
{
  ui.propagationTime->setText(QString::number(f));
}

void StreamLine::SetStepLength(double f)
{
  ui.stepLength->setText(QString::number(f));
}

void StreamLine::SetThreads(int n)
{
  ui.threads->setValue(n);
}

void StreamLine::UseSurfaceMesh(bool b)
{
  ui.useSurfaceGrid->setChecked(b);
}

void StreamLine::UseVolumeMesh(bool b)
{
  ui.useVolumeGrid->setChecked(b);
}

void StreamLine::IntegrateForwards(bool b)
{
  ui.forward->setChecked(b);
}

void StreamLine::IntegrateBackwards(bool b)
{
  ui.backward->setChecked(b);
}

void StreamLine::SetIntegStepLength(double f)
{
  ui.integStepLength->setText(QString::number(f));
}

void StreamLine::SetMinColorVal(double f)
{
  ui.minVal->setText(QString::number(f));
}

void StreamLine::SetMaxColorVal(double f)
{
  ui.maxVal->setText(QString::number(f));
}

void StreamLine::KeepColorLimits(bool b)
{
  ui.keepLimits->setChecked(b);
}

void StreamLine::DrawLines(bool b)
{
  ui.drawLines->setChecked(b);
}

void StreamLine::DrawRibbons(bool b)
{
  ui.drawRibbon->setChecked(b);
}

void StreamLine::SetLineWidth(int n)
{
  ui.lineWidth->setValue(n);
}

void StreamLine::SetRibbonWidth(int n)
{
  ui.ribbonWidth->setValue(n);
}

void StreamLine::UseSphereSource(bool b)
{
  ui.sphereSource->setChecked(b);
}

void StreamLine::UseLineSource(bool b)
{
  ui.lineSource->setChecked(b);
}

void StreamLine::UsePointSource(bool b)
{
  ui.pickSource->setChecked(b);
}

void StreamLine::SetSphereSourceX(double f)
{
  ui.centerX->setText(QString::number(f));
}

void StreamLine::SetSphereSourceY(double f)
{
  ui.centerY->setText(QString::number(f));
}

void StreamLine::SetSphereSourceZ(double f)
{
  ui.centerZ->setText(QString::number(f));
}

void StreamLine::SetSphereSourceRadius(double f)
{
  ui.radius->setText(QString::number(f));
}

void StreamLine::SetSphereSourcePoints(int n)
{
  ui.points->setValue(n);
}

void StreamLine::SetLineSourceStartX(double f)
{
  ui.startX->setText(QString::number(f));
}

void StreamLine::SetLineSourceStartY(double f)
{
  ui.startY->setText(QString::number(f));
}

void StreamLine::SetLineSourceStartZ(double f)
{
  ui.startZ->setText(QString::number(f));
}

void StreamLine::SetLineSourceEndX(double f)
{
  ui.endX->setText(QString::number(f));
}

void StreamLine::SetLineSourceEndY(double f)
{
  ui.endY->setText(QString::number(f));
}

void StreamLine::SetLineSourceEndZ(double f)
{
  ui.endZ->setText(QString::number(f));
}

void StreamLine::SetLineSourcePoints(int n)
{
  ui.lines->setValue(n);
}
