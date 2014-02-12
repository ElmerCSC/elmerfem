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
 *  ElmerGUI isocontour                                                      *
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
#include "isocontour.h"
#include "timestep.h"

#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkContourFilter.h>
#include <vtkDataSetMapper.h>
#include <vtkLookupTable.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkTubeFilter.h>
#include <vtkClipPolyData.h>
#include <vtkPlane.h>

using namespace std;

IsoContour::IsoContour(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  connect(ui.cancelButton, SIGNAL(clicked()), this, SLOT(cancelButtonClicked()));
  connect(ui.applyButton, SIGNAL(clicked()), this, SLOT(applyButtonClicked()));
  connect(ui.okButton, SIGNAL(clicked()), this, SLOT(okButtonClicked()));
  connect(ui.contoursCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(contoursSelectionChanged(int)));
  connect(ui.colorCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(colorSelectionChanged(int)));
  connect(ui.keepContourLimits, SIGNAL(stateChanged(int)), this, SLOT(keepContourLimitsSlot(int)));
  connect(ui.keepColorLimits, SIGNAL(stateChanged(int)), this, SLOT(keepColorLimitsSlot(int)));

  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
}

IsoContour::~IsoContour()
{
}

void IsoContour::applyButtonClicked()
{
  emit(drawIsoContourSignal());
}

void IsoContour::cancelButtonClicked()
{
  emit(hideIsoContourSignal());
  close();
}

void IsoContour::okButtonClicked()
{
  emit(drawIsoContourSignal());
  close();
}

void IsoContour::populateWidgets(VtkPost* vtkPost)
{
  this->scalarField = vtkPost->GetScalarField();
  this->scalarFields = vtkPost->GetScalarFields();

  QString contoursName = ui.contoursCombo->currentText();
  QString colorName = ui.colorCombo->currentText();

  ui.contoursCombo->clear();
  ui.colorCombo->clear();

  for(int i = 0; i < scalarFields; i++) {
    ScalarField *sf = &scalarField[i];
    ui.contoursCombo->addItem(sf->name);
    ui.colorCombo->addItem(sf->name);
  }
  
  this->SetFieldName(contoursName);
  this->SetColorName(colorName);

  contoursSelectionChanged(ui.contoursCombo->currentIndex());
  colorSelectionChanged(ui.colorCombo->currentIndex());

  ui.contourList->clear();
}

void IsoContour::contoursSelectionChanged(int newIndex)
{
  ScalarField *sf = &this->scalarField[newIndex];
  if(!ui.keepContourLimits->isChecked()) {
    ui.contoursMinEdit->setText(QString::number(sf->minVal));
    ui.contoursMaxEdit->setText(QString::number(sf->maxVal));
  }
}

void IsoContour::colorSelectionChanged(int newIndex)
{
  ScalarField *sf = &this->scalarField[newIndex];
  if(!ui.keepColorLimits->isChecked()) {
    ui.colorMinEdit->setText(QString::number(sf->minVal));
    ui.colorMaxEdit->setText(QString::number(sf->maxVal));
  }
}

void IsoContour::keepContourLimitsSlot(int state)
{
  if(state == 0)
    contoursSelectionChanged(ui.contoursCombo->currentIndex());
}

void IsoContour::keepColorLimitsSlot(int state)
{
  if(state == 0)
    colorSelectionChanged(ui.colorCombo->currentIndex());
}

void IsoContour::draw(VtkPost* vtkPost, TimeStep* timeStep)
{
  int contourIndex = ui.contoursCombo->currentIndex();
  QString contourName = ui.contoursCombo->currentText();
  int contours = ui.contoursSpin->value() + 1;
  int lineWidth = ui.lineWidthSpin->value();
  double contourMinVal = ui.contoursMinEdit->text().toDouble();
  double contourMaxVal = ui.contoursMaxEdit->text().toDouble();
  int colorIndex = ui.colorCombo->currentIndex();
  QString colorName = ui.colorCombo->currentText();
  double colorMinVal = ui.colorMinEdit->text().toDouble();
  double colorMaxVal = ui.colorMaxEdit->text().toDouble();
  bool useTubeFilter = ui.useTubeFilter->isChecked();
  int tubeRadius = ui.tubeRadius->value();
  int tubeQuality = ui.tubeQuality->value();
  bool useClip = ui.useClip->isChecked();
  useClip |= vtkPost->GetClipAll();

  // contour list:
  QString contourListText = ui.contourList->text().trimmed();
  QStringList contourList = contourListText.split(";");
  int contourValues = contourList.count();

  QVector<double> contourValue(contourValues);
  for(int i = 0; i < contourValues; i++)
    contourValue[i] = contourList.at(i).toDouble();  
  qSort(contourValue);

  bool useListValues = false;
  if(!contourListText.isEmpty())
    useListValues = true;

  ScalarField* sf = &scalarField[contourIndex];
  int maxDataStepsContour = sf->values / vtkPost->NofNodes();
  int step = timeStep->ui.timeStep->value();
  if(step > maxDataStepsContour) step = maxDataStepsContour;
  if(step > timeStep->maxSteps) step = timeStep->maxSteps;
  int contourOffset = vtkPost->NofNodes() * (step-1);

  sf = &scalarField[colorIndex];
  int maxDataStepsColor = sf->values / vtkPost->NofNodes();
  step = timeStep->ui.timeStep->value();
  if(step > maxDataStepsColor) step = maxDataStepsColor;
  if(step > timeStep->maxSteps) step = timeStep->maxSteps;
  int colorOffset = vtkPost->NofNodes() * (step-1);

  if(contourName == "Null") return;

  // Scalars:
  //----------
  vtkPost->GetSurfaceGrid()->GetPointData()->RemoveArray("IsoContour");
  vtkFloatArray* contourArray = vtkFloatArray::New();
  sf = &scalarField[contourIndex];
  contourArray->SetNumberOfComponents(1);
  contourArray->SetNumberOfTuples(vtkPost->NofNodes());
  contourArray->SetName("IsoContour");
  for(int i = 0; i < vtkPost->NofNodes(); i++)
    contourArray->SetComponent(i, 0, sf->value[i + contourOffset]);
  vtkPost->GetSurfaceGrid()->GetPointData()->AddArray(contourArray);

  vtkPost->GetSurfaceGrid()->GetPointData()->RemoveArray("IsoContourColor");
  vtkFloatArray* colorArray = vtkFloatArray::New();
  sf = &scalarField[colorIndex];
  colorArray->SetName("IsoContourColor");
  colorArray->SetNumberOfComponents(1);
  colorArray->SetNumberOfTuples(vtkPost->NofNodes());
  for(int i = 0; i < vtkPost->NofNodes(); i++)
    colorArray->SetComponent(i, 0, sf->value[i + colorOffset]);
  vtkPost->GetSurfaceGrid()->GetPointData()->AddArray(colorArray);

  // Isocontours:
  //--------------
  vtkContourFilter* iso = vtkContourFilter::New();
  vtkPost->GetSurfaceGrid()->GetPointData()->SetActiveScalars("IsoContour");
  iso->SetInput(vtkPost->GetSurfaceGrid());
  iso->ComputeScalarsOn();
  if(useListValues) {
    iso->SetNumberOfContours(contourValues);
    for(int i = 0; i < contourValues; i++)
      iso->SetValue(i, contourValue[i]);    
  } else {
    iso->GenerateValues(contours, contourMinVal, contourMaxVal);
  }

  // Tube filter:
  //-------------
  vtkTubeFilter* tubes = vtkTubeFilter::New();
  if(useTubeFilter) {
    double r = vtkPost->GetLength() * tubeRadius / 2000.0;
    tubes->SetInputConnection(iso->GetOutputPort());
    tubes->SetNumberOfSides(tubeQuality);
    tubes->SetRadius(r);
  }

  // Apply clip plane:
  //-------------------
  vtkClipPolyData* clipper = vtkClipPolyData::New();
  if(useClip) {
    if(useTubeFilter) {
      clipper->SetInputConnection(tubes->GetOutputPort());
    } else {
      clipper->SetInputConnection(iso->GetOutputPort());
    }
    clipper->SetClipFunction(vtkPost->GetClipPlane());
    clipper->GenerateClipScalarsOn();
    clipper->GenerateClippedOutputOn();
  }

  // Mapper:
  //---------
  vtkDataSetMapper* mapper = vtkDataSetMapper::New();
  if(useClip) {
    mapper->SetInputConnection(clipper->GetOutputPort());
  } else {
    if(useTubeFilter) {
      mapper->SetInputConnection(tubes->GetOutputPort());
    } else {
      mapper->SetInputConnection(iso->GetOutputPort());
    }
  }
  mapper->ScalarVisibilityOn();
  mapper->SelectColorArray("IsoContourColor");
  mapper->SetScalarModeToUsePointFieldData();
  mapper->SetScalarRange(colorMinVal, colorMaxVal);
  mapper->SetLookupTable(vtkPost->GetCurrentLut());
  // mapper->ImmediateModeRenderingOn();

  // Actor & renderer:
  //-------------------
  vtkPost->GetIsoContourActor()->SetMapper(mapper);
  vtkPost->GetIsoContourActor()->GetProperty()->SetLineWidth(lineWidth);
  vtkPost->SetCurrentIsoContourName(colorName);

  // Clean up:
  //----------
  mapper->Delete();
  clipper->Delete();
  tubes->Delete();
  iso->Delete();
  colorArray->Delete();
  contourArray->Delete();
}


// Public slots:
//---------------
QString IsoContour::GetFieldName()
{
  return ui.contoursCombo->currentText();
}

QString IsoContour::GetColorName()
{
  return ui.colorCombo->currentText();
}

bool IsoContour::SetFieldName(QString name)
{
  for(int i = 0; i < ui.contoursCombo->count(); i++) {
    if(ui.contoursCombo->itemText(i) == name) {
      ui.contoursCombo->setCurrentIndex(i);
      return true;
    }
  }
  return false;
}

bool IsoContour::SetColorName(QString name)
{
  for(int i = 0; i < ui.colorCombo->count(); i++) {
    if(ui.colorCombo->itemText(i) == name) {
      ui.colorCombo->setCurrentIndex(i);
      return true;
    }
  }
  return false;
}

void IsoContour::SetMinFieldVal(double f)
{
  ui.contoursMinEdit->setText(QString::number(f));
}

void IsoContour::SetMaxFieldVal(double f)
{
  ui.contoursMaxEdit->setText(QString::number(f));
}

void IsoContour::SetContours(int n)
{
  ui.contoursSpin->setValue(n);
}

void IsoContour::SetContourValues(QString values)
{
  ui.contourList->setText(values);
}

void IsoContour::KeepFieldLimits(bool b)
{
  ui.keepContourLimits->setChecked(b);
}

void IsoContour::SetMinColorVal(double f)
{
  ui.colorMinEdit->setText(QString::number(f));
}

void IsoContour::SetMaxColorVal(double f)
{
  ui.colorMaxEdit->setText(QString::number(f));
}

void IsoContour::KeepColorLimits(bool b)
{
  ui.keepColorLimits->setChecked(b);
}

void IsoContour::UseTubeFilter(bool b)
{
  ui.useTubeFilter->setChecked(b);
}

void IsoContour::UseClipPlane(bool b)
{
  ui.useClip->setChecked(b);
}

void IsoContour::SetLineWidth(int n)
{
  ui.lineWidthSpin->setValue(n);
}

void IsoContour::SetTubeQuality(int n)
{
  ui.tubeQuality->setValue(n);
}

void IsoContour::SetTubeRadius(int n)
{
  ui.tubeRadius->setValue(n);
}
