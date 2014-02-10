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
 *  ElmerGUI meshpoint                                                       *
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
#include "meshpoint.h"
#include "preferences.h"

#include <vtkUnstructuredGrid.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkClipPolyData.h>
#include <vtkPlane.h>

using namespace std;

MeshPoint::MeshPoint(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  setWindowTitle("Mesh points");
  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
}

MeshPoint::~MeshPoint()
{
}

void MeshPoint::draw(VtkPost* vtkPost, Preferences* preferences)
{
  double length = vtkPost->GetLength();
  int pointQuality = preferences->ui.pointQuality->value();
  int pointSize = preferences->ui.pointSize->value();
  bool useSurfaceGrid = preferences->ui.meshPointsSurface->isChecked();
  bool useClip = preferences->ui.meshPointsClip->isChecked();
  useClip |= vtkPost->GetClipAll();

  vtkSphereSource* sphere = vtkSphereSource::New();
  sphere->SetRadius((double)pointSize * length / 2000.0);
  sphere->SetThetaResolution(pointQuality);
  sphere->SetPhiResolution(pointQuality);

  vtkUnstructuredGrid* grid = NULL;

  if(useSurfaceGrid) {
    grid = vtkPost->GetSurfaceGrid();
  } else {
    grid = vtkPost->GetVolumeGrid();
  }

  if(!grid) return;

  if(grid->GetNumberOfPoints() < 1) return;

  vtkGlyph3D* glyph = vtkGlyph3D::New();
  glyph->SetInput(grid);
  glyph->SetSourceConnection(sphere->GetOutputPort());

  vtkClipPolyData* clipper = vtkClipPolyData::New();
  if(useClip) {
    clipper->SetInputConnection(glyph->GetOutputPort());
    clipper->SetClipFunction(vtkPost->GetClipPlane());
    clipper->GenerateClippedOutputOn();
  }

  vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
  if(useClip) {
    mapper->SetInputConnection(clipper->GetOutputPort());
  } else {
    mapper->SetInputConnection(glyph->GetOutputPort());
  }
  mapper->ScalarVisibilityOff();

  vtkPost->GetMeshPointActor()->SetMapper(mapper);
  vtkPost->GetMeshPointActor()->GetProperty()->SetColor(0.5, 0.5, 0.5);

  mapper->Delete();
  clipper->Delete();
  glyph->Delete();
  sphere->Delete();
}
