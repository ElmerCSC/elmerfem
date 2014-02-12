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
 *  ElmerGUI isosurface                                                      *
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
#include <QVector>
#include <iostream>
#include "epmesh.h"
#include "vtkpost.h"
#include "isosurface.h"
#include "timestep.h"

#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkContourFilter.h>
#include <vtkClipPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkDataSetMapper.h>
#include <vtkPlane.h>
#include <vtkLookupTable.h>
#include <vtkProperty.h>
#include <vtkActor.h>

using namespace std;

IsoSurface::IsoSurface(QWidget *parent)
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

IsoSurface::~IsoSurface()
{
}

void IsoSurface::applyButtonClicked()
{
  emit(drawIsoSurfaceSignal());
}

void IsoSurface::cancelButtonClicked()
{
  emit(hideIsoSurfaceSignal());
  close();
}

void IsoSurface::okButtonClicked()
{
  emit(drawIsoSurfaceSignal());
  close();
}

void IsoSurface::populateWidgets(VtkPost* vtkPost)
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

void IsoSurface::contoursSelectionChanged(int newIndex)
{
  ScalarField *sf = &this->scalarField[newIndex];
  if(!ui.keepContourLimits->isChecked()) {
    ui.contoursMinEdit->setText(QString::number(sf->minVal));
    ui.contoursMaxEdit->setText(QString::number(sf->maxVal));
  }
}

void IsoSurface::colorSelectionChanged(int newIndex)
{
  ScalarField *sf = &this->scalarField[newIndex];
  if(!ui.keepColorLimits->isChecked()) {
    ui.colorMinEdit->setText(QString::number(sf->minVal));
    ui.colorMaxEdit->setText(QString::number(sf->maxVal));
  }
}

void IsoSurface::keepContourLimitsSlot(int state)
{
  if(state == 0)
    contoursSelectionChanged(ui.contoursCombo->currentIndex());
}

void IsoSurface::keepColorLimitsSlot(int state)
{
  if(state == 0)
    colorSelectionChanged(ui.colorCombo->currentIndex());
}

void IsoSurface::draw(VtkPost* vtkPost, TimeStep* timeStep)
{
  int contourIndex = ui.contoursCombo->currentIndex();
  QString contourName = ui.contoursCombo->currentText();
  int contours = ui.contoursSpin->value() + 1;
  double contourMinVal = ui.contoursMinEdit->text().toDouble();
  double contourMaxVal = ui.contoursMaxEdit->text().toDouble();
  bool useNormals = ui.normalsCheck->isChecked();
  int colorIndex = ui.colorCombo->currentIndex();
  QString colorName = ui.colorCombo->currentText();
  double colorMinVal = ui.colorMinEdit->text().toDouble();
  double colorMaxVal = ui.colorMaxEdit->text().toDouble();
  int featureAngle = ui.featureAngle->value();
  double opacity = ui.opacitySpin->value() / 100.0;
  bool useClip = ui.clipPlane->isChecked();
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
  int maxDataStepContour = sf->values / vtkPost->NofNodes();
  int step = timeStep->ui.timeStep->value();
  if(step > maxDataStepContour) step = maxDataStepContour;
  if(step > timeStep->maxSteps) step = timeStep->maxSteps;
  int contourOffset = vtkPost->NofNodes() * (step-1);

  sf = &scalarField[colorIndex];
  int maxDataStepColor = sf->values / vtkPost->NofNodes();
  step = timeStep->ui.timeStep->value();
  if(step > maxDataStepColor) step = maxDataStepColor;
  if(step > timeStep->maxSteps) step = timeStep->maxSteps;
  int colorOffset = vtkPost->NofNodes() * (step-1);

  if(contourName == "Null") return;

  // Scalars:
  //----------
  vtkPost->GetVolumeGrid()->GetPointData()->RemoveArray("IsoSurface");
  vtkFloatArray* contourArray = vtkFloatArray::New();
  sf = &scalarField[contourIndex];
  contourArray->SetNumberOfComponents(1);
  contourArray->SetNumberOfTuples(vtkPost->NofNodes());
  contourArray->SetName("IsoSurface");
  for(int i = 0; i < vtkPost->NofNodes(); i++)
    contourArray->SetComponent(i, 0, sf->value[i + contourOffset]);
  vtkPost->GetVolumeGrid()->GetPointData()->AddArray(contourArray);

  vtkPost->GetVolumeGrid()->GetPointData()->RemoveArray("IsoSurfaceColor");
  vtkFloatArray* colorArray = vtkFloatArray::New();
  sf = &scalarField[colorIndex];
  colorArray->SetName("IsoSurfaceColor");
  colorArray->SetNumberOfComponents(1);
  colorArray->SetNumberOfTuples(vtkPost->NofNodes());
  for(int i = 0; i < vtkPost->NofNodes(); i++)
    colorArray->SetComponent(i, 0, sf->value[i + colorOffset]);
  vtkPost->GetVolumeGrid()->GetPointData()->AddArray(colorArray);

  // Isosurfaces:
  //--------------
  vtkContourFilter* iso = vtkContourFilter::New();
  vtkPost->GetVolumeGrid()->GetPointData()->SetActiveScalars("IsoSurface");
  iso->SetInput(vtkPost->GetVolumeGrid());
  iso->ComputeScalarsOn();
  if(useListValues) {
    iso->SetNumberOfContours(contourValues);
    for(int i = 0; i < contourValues; i++)
      iso->SetValue(i, contourValue[i]);
  } else {
    iso->GenerateValues(contours, contourMinVal, contourMaxVal);
  }

  // Apply the clip plane:
  //-----------------------
  vtkClipPolyData* clipper = vtkClipPolyData::New();

  if(useClip) {
    clipper->SetInputConnection(iso->GetOutputPort());
    clipper->SetClipFunction(vtkPost->GetClipPlane());
    clipper->GenerateClipScalarsOn();
    clipper->GenerateClippedOutputOn();
  }

  // Normals:
  //---------
  vtkPolyDataNormals* normals = vtkPolyDataNormals::New();

  if(useNormals) {
    if(useClip) {
      normals->SetInputConnection(clipper->GetOutputPort());
    } else {
      normals->SetInputConnection(iso->GetOutputPort());
    }
    normals->SetFeatureAngle(featureAngle);
  }

  // Mapper:
  //--------
  vtkDataSetMapper* mapper = vtkDataSetMapper::New();
  
  if(useNormals) {
    mapper->SetInputConnection(normals->GetOutputPort());
  } else {
    if(useClip) {
      mapper->SetInputConnection(clipper->GetOutputPort());      
    } else {
      mapper->SetInputConnection(iso->GetOutputPort());
    }
  }

  mapper->ScalarVisibilityOn();
  mapper->SelectColorArray("IsoSurfaceColor");
  mapper->SetScalarModeToUsePointFieldData();
  mapper->SetScalarRange(colorMinVal, colorMaxVal);
  mapper->SetLookupTable(vtkPost->GetCurrentLut());
  // mapper->ImmediateModeRenderingOn();

  // Actor:
  //--------
  vtkPost->GetIsoSurfaceActor()->SetMapper(mapper);
  vtkPost->GetIsoSurfaceActor()->GetProperty()->SetOpacity(opacity);
  vtkPost->SetCurrentIsoSurfaceName(colorName);

  // Clean up:
  //----------
  mapper->Delete();
  normals->Delete();
  clipper->Delete();
  iso->Delete();
  colorArray->Delete();
  contourArray->Delete();
}

// Public slots:
//---------------
QString IsoSurface::GetFieldName()
{
  return ui.contoursCombo->currentText();
}

QString IsoSurface::GetColorName()
{
  return ui.colorCombo->currentText();
}

bool IsoSurface::SetFieldName(QString name)
{
  for(int i = 0; i < ui.contoursCombo->count(); i++) {
    if(ui.contoursCombo->itemText(i) == name) {
      ui.contoursCombo->setCurrentIndex(i);
      return true;
    }
  }
  return false;
}

bool IsoSurface::SetColorName(QString name)
{
  for(int i = 0; i < ui.colorCombo->count(); i++) {
    if(ui.colorCombo->itemText(i) == name) {
      ui.colorCombo->setCurrentIndex(i);
      return true;
    }
  }
  return false;
}

void IsoSurface::SetMinFieldVal(double f)
{
  ui.contoursMinEdit->setText(QString::number(f));
}

void IsoSurface::SetMaxFieldVal(double f)
{
  ui.contoursMaxEdit->setText(QString::number(f));
}

void IsoSurface::SetContours(int n)
{
  ui.contoursSpin->setValue(n);
}

void IsoSurface::SetContourValues(QString values)
{
  ui.contourList->setText(values);
}

void IsoSurface::KeepFieldLimits(bool b)
{
  ui.keepContourLimits->setChecked(b);
}

void IsoSurface::SetMinColorVal(double f)
{
  ui.colorMinEdit->setText(QString::number(f));
}

void IsoSurface::SetMaxColorVal(double f)
{
  ui.colorMaxEdit->setText(QString::number(f));
}

void IsoSurface::KeepColorLimits(bool b)
{
  ui.keepColorLimits->setChecked(b);
}

void IsoSurface::ComputeNormals(bool b)
{
  ui.normalsCheck->setChecked(b);
}

void IsoSurface::UseClipPlane(bool b)
{
  ui.clipPlane->setChecked(b);
}

void IsoSurface::SetFeatureAngle(int n)
{
  ui.featureAngle->setValue(n);
}

void IsoSurface::SetOpacity(int n)
{
  ui.opacitySpin->setValue(n);
}

