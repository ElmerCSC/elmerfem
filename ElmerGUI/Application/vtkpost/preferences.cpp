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
 *  ElmerGUI preferences                                                     *
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
#include "vtkpost.h"
#include "preferences.h"

using namespace std;

Preferences::Preferences(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  connect(ui.cancelButton, SIGNAL(clicked()), this, SLOT(cancelButtonClicked()));
  connect(ui.okButton, SIGNAL(clicked()), this, SLOT(okButtonClicked()));
  connect(ui.applyButton, SIGNAL(clicked()), this, SLOT(applyButtonClicked()));

  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
}

Preferences::~Preferences()
{
}

void Preferences::okButtonClicked()
{
  emit(redrawSignal());
  close();
}

void Preferences::applyButtonClicked()
{
  emit(redrawSignal());  
}

void Preferences::cancelButtonClicked()
{
  close();
}


// Public slots:
//---------------
int Preferences::GetFeatureAngle()
{
  return ui.angleSpin->value();
}

void Preferences::UseSurfaceMeshForPoints(bool b)
{
  ui.meshPointsSurface->setChecked(b);
}
void Preferences::UseVolumeMeshForPoints(bool b)
{
  ui.meshPointsVolume->setChecked(b);
}
void Preferences::SetPointSize(int n)
{
  ui.pointSize->setValue(n);
}
void Preferences::SetPointQuality(int n)
{
  ui.pointQuality->setValue(n);
}
void Preferences::UseClipPlaneForPoints(bool b)
{
  ui.meshPointsClip->setChecked(b);
}

void Preferences::UseSurfaceMeshForEdges(bool b)
{
  ui.meshEdgesSurface->setChecked(b);
}

void Preferences::UseVolumeMeshForEdges(bool b)
{
  ui.meshEdgesVolume->setChecked(b);
}

void Preferences::UseTubeFilterForEdges(bool b)
{
  ui.meshEdgeTubes->setChecked(b);
}

void Preferences::UseClipPlaneForEdges(bool b)
{
  ui.meshEdgesClip->setChecked(b);
}

void Preferences::SetLineWidthForEdges(int n)
{
  ui.meshLineWidth->setValue(n);
}

void Preferences::SetTubeQualityForEdges(int n)
{
  ui.meshEdgeTubeQuality->setValue(n);
}

void Preferences::SetTubeRadiusForEdges(int n)
{
  ui.meshEdgeTubeRadius->setValue(n);
}

void Preferences::UseSurfaceMeshForFeatureEdges(bool b)
{
  ui.surfaceRButton->setChecked(b);
}

void Preferences::UseVolumeMeshForFeatureEdges(bool b)
{
  ui.volumeRButton->setChecked(b);
}

void Preferences::UseTubeFilterForFeatureEdges(bool b)
{
  ui.featureEdgeTubes->setChecked(b);
}

void Preferences::UseClipPlaneForFeatureEdges(bool b)
{
  ui.featureEdgesClip->setChecked(b);
}

void Preferences::DrawBoundaryEdges(bool b)
{
  ui.drawBoundaryEdges->setChecked(b);
}

void Preferences::SetFeatureAngle(int n)
{
  ui.angleSpin->setValue(n);
}

void Preferences::SetLineWidthForFeatureEdges(int n)
{
  ui.lineWidthSpin->setValue(n);
}

void Preferences::SetTubeQualityForFeatureEdges(int n)
{
  ui.featureEdgeTubeQuality->setValue(n);
}

void Preferences::SetTubeRadiusForFeatureEdges(int n)
{
  ui.featureEdgeTubeRadius->setValue(n);
}

void Preferences::SetClipPlaneOx(double f)
{
  ui.clipPointX->setText(QString::number(f));
}

void Preferences::SetClipPlaneOy(double f)
{
  ui.clipPointY->setText(QString::number(f));
}

void Preferences::SetClipPlaneOz(double f)
{
  ui.clipPointZ->setText(QString::number(f));
}

void Preferences::SetClipPlaneNx(double f)
{
  ui.clipNormalX->setText(QString::number(f));
}

void Preferences::SetClipPlaneNy(double f)
{
  ui.clipNormalY->setText(QString::number(f));
}

void Preferences::SetClipPlaneNz(double f)
{
  ui.clipNormalZ->setText(QString::number(f));
}
