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
 *  ElmerGUI surface                                                         *
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
#include "surface.h"
#include "timestep.h"

#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkGeometryFilter.h>
#include <vtkClipPolyData.h>
#include <vtkPlane.h>
#include <vtkPolyDataNormals.h>
#include <vtkDataSetMapper.h>
#include <vtkLookupTable.h>
#include <vtkProperty.h>
#include <vtkActor.h>

using namespace std;

Surface::Surface(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  connect(ui.cancelButton, SIGNAL(clicked()), this, SLOT(cancelButtonClicked()));
  connect(ui.applyButton, SIGNAL(clicked()), this, SLOT(applyButtonClicked()));
  connect(ui.okButton, SIGNAL(clicked()), this, SLOT(okButtonClicked()));
  connect(ui.surfaceCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(surfaceSelectionChanged(int)));
  connect(ui.keepLimits, SIGNAL(stateChanged(int)), this, SLOT(keepLimitsSlot(int)));

  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
}

Surface::~Surface()
{
}

void Surface::cancelButtonClicked()
{
  emit(hideSurfaceSignal());
  close();
}

void Surface::applyButtonClicked()
{
  emit(drawSurfaceSignal());
}

void Surface::okButtonClicked()
{
  applyButtonClicked();
  close();
}

void Surface::populateWidgets(VtkPost* vtkPost)
{
  this->scalarField = vtkPost->GetScalarField();
  this->scalarFields = vtkPost->GetScalarFields();

  QString name = ui.surfaceCombo->currentText();

  ui.surfaceCombo->clear();

  for(int i = 0; i < scalarFields; i++) {
    ScalarField *sf = &scalarField[i];
    QString name = sf->name;
    ui.surfaceCombo->addItem(sf->name);
  }

  this->SetFieldName(name);

  surfaceSelectionChanged(ui.surfaceCombo->currentIndex());
}

void Surface::surfaceSelectionChanged(int newIndex)
{
  ScalarField *sf = &this->scalarField[newIndex];
  if(!ui.keepLimits->isChecked()) {
    ui.minEdit->setText(QString::number(sf->minVal));
    ui.maxEdit->setText(QString::number(sf->maxVal));
  }
}

void Surface::keepLimitsSlot(int state)
{
  if(state == 0)
    surfaceSelectionChanged(ui.surfaceCombo->currentIndex());
}

void Surface::draw(VtkPost* vtkPost, TimeStep* timeStep)
{
  int surfaceIndex = ui.surfaceCombo->currentIndex();
  QString surfaceName = ui.surfaceCombo->currentText();
  double minVal = ui.minEdit->text().toDouble();
  double maxVal = ui.maxEdit->text().toDouble();
  bool useNormals = ui.useNormals->isChecked();
  int featureAngle = ui.featureAngle->value();
  double opacity = ui.opacitySpin->value() / 100.0;
  bool useClip = ui.clipPlane->isChecked();
  useClip |= vtkPost->GetClipAll();

  ScalarField* sf = &scalarField[surfaceIndex];
  int maxDataSteps = sf->values / vtkPost->NofNodes();
  int step = timeStep->ui.timeStep->value();
  if(step > maxDataSteps) step = maxDataSteps;
  if(step > timeStep->maxSteps) step = timeStep->maxSteps;
  int offset = vtkPost->NofNodes() * (step-1);

  // Scalars:
  //---------
  vtkPost->GetSurfaceGrid()->GetPointData()->RemoveArray("Surface");
  vtkFloatArray* scalars = vtkFloatArray::New();
  scalars->SetNumberOfComponents(1);
  scalars->SetNumberOfTuples(vtkPost->NofNodes());
  scalars->SetName("Surface");
  for(int i = 0; i < vtkPost->NofNodes(); i++)
    scalars->SetComponent(i, 0, sf->value[i + offset]);  
  vtkPost->GetSurfaceGrid()->GetPointData()->AddArray(scalars);

  // Convert from vtkUnstructuredGrid to vtkPolyData:
  //-------------------------------------------------
  vtkGeometryFilter* filter = vtkGeometryFilter::New();

  filter->SetInput(vtkPost->GetSurfaceGrid());
  filter->GetOutput()->ReleaseDataFlagOn();

  // Apply the clip plane:
  //-----------------------
  vtkClipPolyData *clipper = vtkClipPolyData::New();

  if(useClip) {
    clipper->SetInputConnection(filter->GetOutputPort());
    clipper->SetClipFunction(vtkPost->GetClipPlane());
    clipper->GenerateClipScalarsOn();
    clipper->GenerateClippedOutputOn();
  }

  // Normals:
  //---------
  vtkPolyDataNormals *normals = vtkPolyDataNormals::New();
  
  if(useNormals) {
    if(useClip) {
      normals->SetInputConnection(clipper->GetOutputPort());
    } else {
      normals->SetInputConnection(filter->GetOutputPort());
    }
    normals->SetFeatureAngle(featureAngle);
  }

  // Mapper:
  //--------
  vtkDataSetMapper *mapper = vtkDataSetMapper::New();

  if(useNormals) {
    mapper->SetInputConnection(normals->GetOutputPort());
  } else {
    if(useClip) {
      mapper->SetInputConnection(clipper->GetOutputPort());
    } else {
      mapper->SetInput(vtkPost->GetSurfaceGrid());
    }
  }

  mapper->SetScalarModeToUsePointFieldData();
  mapper->SelectColorArray("Surface");
  mapper->ScalarVisibilityOn();
  mapper->SetScalarRange(minVal, maxVal);
  mapper->SetResolveCoincidentTopologyToPolygonOffset();
  mapper->SetLookupTable(vtkPost->GetCurrentLut());
  // mapper->ImmediateModeRenderingOn();

  // Actor:
  //--------
  vtkPost->GetSurfaceActor()->SetMapper(mapper);
  vtkPost->GetSurfaceActor()->GetProperty()->SetOpacity(opacity);
  vtkPost->SetCurrentSurfaceName(sf->name);

  // Clean up:
  //-----------
  mapper->Delete();
  normals->Delete();
  clipper->Delete();
  filter->Delete();
  scalars->Delete();
}

// Public slots:
//---------------
QString Surface::GetFieldName()
{
  return ui.surfaceCombo->currentText();
}

bool Surface::SetFieldName(QString name)
{
  for(int i = 0; i < ui.surfaceCombo->count(); i++) {
    if(ui.surfaceCombo->itemText(i) == name) {
      ui.surfaceCombo->setCurrentIndex(i);
      return true;
    }
  }
  return false;
}

void Surface::SetMinVal(double f)
{
  ui.minEdit->setText(QString::number(f));
}

void Surface::SetMaxVal(double f)
{
  ui.maxEdit->setText(QString::number(f));
}

void Surface::KeepLimits(bool b)
{
  ui.keepLimits->setChecked(b);
}

void Surface::SetComputeNormals(bool b)
{
  ui.useNormals->setChecked(b);
}

void Surface::SetFeatureAngle(int n)
{
  ui.featureAngle->setValue(n);
}

void Surface::SetOpacity(int n)
{
  ui.opacitySpin->setValue(n);
}

void Surface::SetClipPlane(bool b)
{
  ui.clipPlane->setChecked(b);
}
