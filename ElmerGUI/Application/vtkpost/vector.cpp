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
 *  ElmerGUI vector                                                          *
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
#include "vector.h"
#include "timestep.h"

#include <vtkPolyDataNormals.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkGlyph3D.h>
#include <vtkArrowSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkLookupTable.h>
#include <vtkActor.h>
#include <vtkClipPolyData.h>
#include <vtkPlane.h>
#include <vtkMaskPoints.h>

using namespace std;

Vector::Vector(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  connect(ui.cancelButton, SIGNAL(clicked()), this, SLOT(cancelButtonClicked()));
  connect(ui.applyButton, SIGNAL(clicked()), this, SLOT(applyButtonClicked()));
  connect(ui.okButton, SIGNAL(clicked()), this, SLOT(okButtonClicked()));
  connect(ui.colorCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(colorSelectionChanged(int)));
  connect(ui.thresholdCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(thresholdSelectionChanged(int)));
  connect(ui.keepLimits, SIGNAL(stateChanged(int)), this, SLOT(keepLimitsSlot(int)));
  connect(ui.keepThresholdLimits, SIGNAL(stateChanged(int)), this, SLOT(keepThresholdLimitsSlot(int)));

  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
}

Vector::~Vector()
{
}

void Vector::applyButtonClicked()
{
  emit(drawVectorSignal());
}

void Vector::cancelButtonClicked()
{
  emit(hideVectorSignal());
  close();
}

void Vector::okButtonClicked()
{
  emit(drawVectorSignal());
  close();
}

void Vector::populateWidgets(VtkPost *vtkPost)
{
  this->scalarField = vtkPost->GetScalarField();
  this->scalarFields = vtkPost->GetScalarFields();

  // Arrow:
  //-------
  ui.vectorCombo->clear();

  int index = -1;
  for(int i = 0; i < scalarFields; i++) {
    ScalarField *sf = &scalarField[i];
    QString name = sf->name;
    if((index = name.indexOf("_x")) >= 0) {
      ui.vectorCombo->addItem(name.mid(0, index));
    }
  }

  // Color:
  //-------
  QString name = ui.colorCombo->currentText();

  ui.colorCombo->clear();

  for(int i = 0; i < scalarFields; i++) {
    ScalarField *sf = &scalarField[i];
    ui.colorCombo->addItem(sf->name);
  }

  this->SetColorName(name);

  colorSelectionChanged(ui.colorCombo->currentIndex());

  // Threshold:
  //-----------
  name = ui.thresholdCombo->currentText();

  ui.thresholdCombo->clear();
  
  for(int i = 0; i < scalarFields; i++) {
    ScalarField *sf = &scalarField[i];
    ui.thresholdCombo->addItem(sf->name);
  }

  this->SetThresholdName(name);

  thresholdSelectionChanged(ui.thresholdCombo->currentIndex());
}

void Vector::colorSelectionChanged(int newIndex)
{
  ScalarField *sf = &this->scalarField[newIndex];
  if(!ui.keepLimits->isChecked()) {
    ui.minVal->setText(QString::number(sf->minVal));
    ui.maxVal->setText(QString::number(sf->maxVal));
  }
}

void Vector::thresholdSelectionChanged(int newIndex)
{
  ScalarField *sf = &this->scalarField[newIndex];
  if(!ui.keepThresholdLimits->isChecked()) {
    ui.thresholdMin->setText(QString::number(sf->minVal));
    ui.thresholdMax->setText(QString::number(sf->maxVal));
  }
}

void Vector::keepLimitsSlot(int state)
{
  if(state == 0)
    colorSelectionChanged(ui.colorCombo->currentIndex());
}

void Vector::keepThresholdLimitsSlot(int state)
{
  if(state == 0)
    thresholdSelectionChanged(ui.thresholdCombo->currentIndex());
}

void Vector::draw(VtkPost* vtkPost, TimeStep* timeStep)
{
  QString vectorName = ui.vectorCombo->currentText();

  if(vectorName.isEmpty()) return;

  int i, j, index = -1;
  for(i = 0; i < scalarFields; i++) {
    ScalarField *sf = &scalarField[i];
    QString name = sf->name;
    if((j = name.indexOf("_x")) >= 0) {
      if(vectorName == name.mid(0, j)) {
	index = i;
	break;
      }
    }
  }

  if(index < 0) return;

  int colorIndex = ui.colorCombo->currentIndex();
  QString colorName = ui.colorCombo->currentText();
  double minVal = ui.minVal->text().toDouble();
  double maxVal = ui.maxVal->text().toDouble();
  int quality = ui.qualitySpin->value();
  double scaleMultiplier = ui.scaleSpin->value() / 100.0;
  bool scaleByMagnitude = ui.scaleByMagnitude->isChecked();
  bool useClip = ui.useClip->isChecked();
  useClip |= vtkPost->GetClipAll();
  bool useNormals = ui.useNormals->isChecked();
  int everyNth = ui.everyNth->value();
  bool randomMode = ui.randomMode->isChecked();
  bool useThreshold = ui.useThreshold->isChecked();
  int thresholdIndex = ui.thresholdCombo->currentIndex();
  double thresholdMin = ui.thresholdMin->text().toDouble();
  double thresholdMax = ui.thresholdMax->text().toDouble();

  ScalarField* sf_x = &scalarField[index + 0];
  ScalarField* sf_y = &scalarField[index + 1];
  ScalarField* sf_z = &scalarField[index + 2];
  int maxDataStepsVector = sf_x->values / vtkPost->NofNodes();
  int step = timeStep->ui.timeStep->value();
  if(step > maxDataStepsVector) step = maxDataStepsVector;
  if(step > timeStep->maxSteps) step = timeStep->maxSteps;
  int vectorOffset = vtkPost->NofNodes() * (step-1);

  ScalarField* sf = &scalarField[colorIndex];
  int maxDataStepsColor = sf->values / vtkPost->NofNodes();
  step = timeStep->ui.timeStep->value();
  if(step > maxDataStepsColor) step = maxDataStepsColor;
  if(step > timeStep->maxSteps) step = timeStep->maxSteps;
  int colorOffset = vtkPost->NofNodes() * (step-1);

  ScalarField* sf_threshold = &scalarField[thresholdIndex];
  int maxDataStepsThreshold = sf_threshold->values / vtkPost->NofNodes();
  step = timeStep->ui.timeStep->value();
  if(step > maxDataStepsThreshold) step = maxDataStepsThreshold;
  if(step > timeStep->maxSteps) step = timeStep->maxSteps;
  int thresholdOffset = vtkPost->NofNodes() * (step-1);

  // Vector data:
  //-------------
  vtkPost->GetVolumeGrid()->GetPointData()->RemoveArray("VectorData");
  vtkFloatArray* vectorData = vtkFloatArray::New();
  vectorData->SetNumberOfComponents(3);
  vectorData->SetNumberOfTuples(vtkPost->NofNodes());
  vectorData->SetName("VectorData");
  double scaleFactor = 0.0;
  for(int i = 0; i < vtkPost->NofNodes(); i++) {
    double val_x  = sf_x->value[i + vectorOffset];
    double val_y  = sf_y->value[i + vectorOffset];
    double val_z  = sf_z->value[i + vectorOffset];

    if(useThreshold) {
      double thresholdVal = sf_threshold->value[i + thresholdOffset];
      if((thresholdVal < thresholdMin) || (thresholdVal > thresholdMax)) {
	val_x = 0; val_y = 0; val_z = 0;
      }
    }

    double absval = sqrt(val_x*val_x + val_y*val_y + val_z*val_z);
    if(absval > scaleFactor) scaleFactor = absval;

    vectorData->SetComponent(i, 0, val_x);
    vectorData->SetComponent(i, 1, val_y); 
    vectorData->SetComponent(i, 2, val_z); 
  }
  vtkPost->GetVolumeGrid()->GetPointData()->AddArray(vectorData);

  // Size of volume grid:
  //---------------------
  double length = vtkPost->GetVolumeGrid()->GetLength();
  if ( scaleByMagnitude ) 
    scaleFactor = scaleFactor * 100.0 / length;
  else
    scaleFactor = 100.0 / length;

  // Color data:
  //-------------
  vtkPost->GetVolumeGrid()->GetPointData()->RemoveArray("VectorColor");
  vtkFloatArray* vectorColor = vtkFloatArray::New();
  vectorColor->SetNumberOfComponents(1);
  vectorColor->SetNumberOfTuples(vtkPost->NofNodes());
  vectorColor->SetName("VectorColor");
  for(int i = 0; i < vtkPost->NofNodes(); i++) 
    vectorColor->SetComponent(i, 0, sf->value[i + colorOffset]); 
  vtkPost->GetVolumeGrid()->GetPointData()->AddArray(vectorColor);

  // Mask points:
  //-------------
  vtkMaskPoints* maskPoints = vtkMaskPoints::New();
  maskPoints->SetInput(vtkPost->GetVolumeGrid());
  if(randomMode) {
    maskPoints->RandomModeOn();
  } else {
    maskPoints->RandomModeOff();
  }
  maskPoints->SetOnRatio(everyNth);

  // Glyphs:
  //---------
  vtkPost->GetVolumeGrid()->GetPointData()->SetActiveVectors("VectorData"); 
  vtkGlyph3D* glyph = vtkGlyph3D::New();
  vtkArrowSource* arrow = vtkArrowSource::New();
  arrow->SetTipResolution(quality);
  arrow->SetShaftResolution(quality);
  glyph->SetInputConnection(maskPoints->GetOutputPort());
  glyph->SetSourceConnection(arrow->GetOutputPort());
  glyph->SetVectorModeToUseVector();

  glyph->SetScaleFactor(scaleMultiplier/scaleFactor);
  if ( scaleByMagnitude )
    glyph->SetScaleModeToScaleByVector();
  else
    glyph->SetScaleModeToDataScalingOff();

  glyph->SetColorModeToColorByScale();
  
  vtkClipPolyData* clipper = vtkClipPolyData::New();

  if(useClip) {
    clipper->SetInputConnection(glyph->GetOutputPort());
    clipper->SetClipFunction(vtkPost->GetClipPlane());
    clipper->GenerateClipScalarsOn();
    clipper->GenerateClippedOutputOn();
  }

  vtkPolyDataNormals *normals = vtkPolyDataNormals::New();

  if(useNormals) {
    if(useClip) {
      normals->SetInputConnection(clipper->GetOutputPort());
    } else {
      normals->SetInputConnection(glyph->GetOutputPort());
    }
    normals->SetFeatureAngle(80.0);
  }


  vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();

  if ( useNormals ) {
      mapper->SetInputConnection(normals->GetOutputPort());    
  } else {
    if(useClip) {
      mapper->SetInputConnection(clipper->GetOutputPort());    
    } else {
      mapper->SetInputConnection(glyph->GetOutputPort());
    }
  }
  mapper->SetScalarModeToUsePointFieldData();
  mapper->ScalarVisibilityOn();
  mapper->SetScalarRange(minVal, maxVal);
  mapper->SelectColorArray("VectorColor");
  mapper->SetLookupTable(vtkPost->GetCurrentLut());
  // mapper->ImmediateModeRenderingOn();

  vtkPost->GetVectorActor()->SetMapper(mapper);
  vtkPost->SetCurrentVectorName(colorName);

  maskPoints->Delete();
  normals->Delete();
  mapper->Delete();
  clipper->Delete();
  arrow->Delete();
  glyph->Delete();
  vectorColor->Delete();
  vectorData->Delete();
}

// Public slots:
//---------------
QString Vector::GetFieldName()
{
  return ui.vectorCombo->currentText();
}

QString Vector::GetColorName()
{
  return ui.colorCombo->currentText();
}

QString Vector::GetThresholdName()
{
  return ui.thresholdCombo->currentText();
}

bool Vector::SetFieldName(QString name)
{
  for(int i = 0; i < ui.vectorCombo->count(); i++) {
    if(ui.vectorCombo->itemText(i) == name) {
      ui.vectorCombo->setCurrentIndex(i); 
      return true;
    }
  }
  return false;
}

bool Vector::SetColorName(QString name)
{
  for(int i = 0; i < ui.colorCombo->count(); i++) {
    if(ui.colorCombo->itemText(i) == name) {
      ui.colorCombo->setCurrentIndex(i); 
      return true;
    }
  }
  return false;
}

bool Vector::SetThresholdName(QString name)
{
  for(int i = 0; i < ui.thresholdCombo->count(); i++) {
    if(ui.thresholdCombo->itemText(i) == name) {
      ui.thresholdCombo->setCurrentIndex(i);
      return true;
    }
  }
  return false;
}

void Vector::SetLength(int n)
{
  ui.scaleSpin->setValue(n);
}

void Vector::SetQuality(int n)
{
  ui.qualitySpin->setValue(n);
}

void Vector::SetEveryNth(int n)
{
  ui.everyNth->setValue(n);
}

void Vector::SetClipPlane(bool b)
{
  ui.useClip->setChecked(b);
}

void Vector::ComputeNormals(bool b)
{
  ui.useNormals->setChecked(b);
}

void Vector::SetRandomMode(bool b)
{
  ui.randomMode->setChecked(b);
}

void Vector::ScaleByMagnitude(bool b)
{
  ui.scaleByMagnitude->setChecked(b);
}

void Vector::SetMinColorVal(double f)
{
  ui.minVal->setText(QString::number(f));
}

void Vector::SetMaxColorVal(double f)
{
  ui.maxVal->setText(QString::number(f));
}

void Vector::KeepColorLimits(bool b)
{
  ui.keepLimits->setChecked(b);
}

void Vector::SetMinThresholdVal(double f)
{
  ui.thresholdMin->setText(QString::number(f));
}

void Vector::SetMaxThresholdVal(double f)
{
  ui.thresholdMax->setText(QString::number(f));
}

void Vector::UseThreshold(bool b)
{
  ui.useThreshold->setChecked(b);
}

void Vector::KeepThresholdLimits(bool b)
{
  ui.keepThresholdLimits->setChecked(b);
}

