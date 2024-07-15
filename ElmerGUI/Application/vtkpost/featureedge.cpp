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
 *  ElmerGUI featureedge                                                     *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Authors: Mikko Lyly, Juha Ruokolainen and Peter R�back                   *
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
#include "featureedge.h"
#include "preferences.h"

#include <vtkGeometryFilter.h>
#include <vtkFeatureEdges.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>
#include <vtkClipPolyData.h>
#include <vtkPlane.h>

using namespace std;

FeatureEdge::FeatureEdge(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  setWindowTitle("Feature edges");
  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
  
  actor = vtkActor::New();
}

FeatureEdge::~FeatureEdge()
{
  actor->Delete();
}

/*
// The original draw() function which draws all the groups by one FeatureEdge instance
// using one vtkUnstructuredGrid. This ends up with boundary of two groups not drawn. 
void FeatureEdge::draw(VtkPost* vtkPost, Preferences* preferences)
{ 
  bool useSurfaceGrid = preferences->ui.surfaceRButton->isChecked();
  int featureAngle = preferences->ui.angleSpin->value();
  int lineWidth = preferences->ui.lineWidthSpin->value();
  bool useTubeFilter = preferences->ui.featureEdgeTubes->isChecked();
  int tubeQuality = preferences->ui.featureEdgeTubeQuality->value();
  int radius = preferences->ui.featureEdgeTubeRadius->value();  
  bool useClip = preferences->ui.featureEdgesClip->isChecked();
  bool boundaryEdges = preferences->ui.drawBoundaryEdges->isChecked();
  useClip |= vtkPost->GetClipAll();

  vtkUnstructuredGrid* grid = NULL;

  if(useSurfaceGrid) {
    grid = vtkPost->GetSurfaceGrid();
  } else {
    grid = vtkPost->GetVolumeGrid();
  }

  if(!grid) return;

  if(grid->GetNumberOfCells() < 1) return;

  // Convert from vtkUnstructuredGrid to vtkPolyData:
  vtkGeometryFilter* filter = vtkGeometryFilter::New();
#if VTK_MAJOR_VERSION <= 5
  filter->SetInput(grid);
#else
  filter->SetInputData(grid);
#endif
  // filter->GetOutput()->ReleaseDataFlagOn();

  vtkFeatureEdges* edges = vtkFeatureEdges::New();
  edges->SetInputConnection(filter->GetOutputPort());
  edges->SetFeatureAngle(featureAngle);
  edges->NonManifoldEdgesOn();
  edges->ManifoldEdgesOn();
  if(boundaryEdges) {
    edges->BoundaryEdgesOn();
  } else {
    edges->BoundaryEdgesOff();
  }

  vtkTubeFilter* tubes = vtkTubeFilter::New();
  if(useTubeFilter) {
    double r = vtkPost->GetLength() * radius / 2000.0;
    tubes->SetInputConnection(edges->GetOutputPort());
    tubes->SetNumberOfSides(tubeQuality);
    tubes->SetRadius(r);
  }

  vtkClipPolyData* clipper = vtkClipPolyData::New();
  if(useClip) {
    if(useTubeFilter) {
      clipper->SetInputConnection(tubes->GetOutputPort());
    } else {
      clipper->SetInputConnection(edges->GetOutputPort());
    }
    clipper->SetClipFunction(vtkPost->GetClipPlane());
    clipper->GenerateClippedOutputOn();    
  }

  vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
  if(useClip) {
    mapper->SetInputConnection(clipper->GetOutputPort());
  } else {
    if(useTubeFilter) {
      mapper->SetInputConnection(tubes->GetOutputPort());
    } else {
      mapper->SetInputConnection(edges->GetOutputPort());
    }
  }
  mapper->ScalarVisibilityOff();
  mapper->SetResolveCoincidentTopologyToPolygonOffset();
  // mapper->ImmediateModeRenderingOn();

  vtkPost->GetFeatureEdgeActor()->GetProperty()->SetLineWidth(lineWidth);
  vtkPost->GetFeatureEdgeActor()->GetProperty()->SetColor(0, 0, 0);
  vtkPost->GetFeatureEdgeActor()->SetMapper(mapper);


  mapper->Delete();
  clipper->Delete();
  tubes->Delete();
  edges->Delete();
  filter->Delete();
}
*/


// The new draw() function which draws one group by one FeatureEdge instance 
// using the specified vtkUnstructuredGrid to draw boundary of two groups.
void FeatureEdge::draw(VtkPost* vtkPost, Preferences* preferences, vtkUnstructuredGrid* grid)
{ 
  bool useSurfaceGrid = preferences->ui.surfaceRButton->isChecked();
  int featureAngle = preferences->ui.angleSpin->value();
  int lineWidth = preferences->ui.lineWidthSpin->value();
  bool useTubeFilter = preferences->ui.featureEdgeTubes->isChecked();
  int tubeQuality = preferences->ui.featureEdgeTubeQuality->value();
  int radius = preferences->ui.featureEdgeTubeRadius->value();  
  bool useClip = preferences->ui.featureEdgesClip->isChecked();
  bool boundaryEdges = preferences->ui.drawBoundaryEdges->isChecked();
  useClip |= vtkPost->GetClipAll();
  QColor color = preferences->getFeatureEdgeColor();
  

  if(!grid) return;

  if(grid->GetNumberOfCells() < 1) return;

  // Convert from vtkUnstructuredGrid to vtkPolyData:
  vtkGeometryFilter* filter = vtkGeometryFilter::New();
#if VTK_MAJOR_VERSION <= 5
  filter->SetInput(grid);
#else
  filter->SetInputData(grid);
#endif
  // filter->GetOutput()->ReleaseDataFlagOn();

  vtkFeatureEdges* edges = vtkFeatureEdges::New();
  edges->SetInputConnection(filter->GetOutputPort());
  edges->SetFeatureAngle(featureAngle);
  edges->NonManifoldEdgesOn();
  edges->ManifoldEdgesOn();
  if(boundaryEdges) {
    edges->BoundaryEdgesOn();
  } else {
    edges->BoundaryEdgesOff();
  }

  vtkTubeFilter* tubes = vtkTubeFilter::New();
  if(useTubeFilter) {
    double r = vtkPost->GetLength() * radius / 2000.0;
    tubes->SetInputConnection(edges->GetOutputPort());
    tubes->SetNumberOfSides(tubeQuality);
    tubes->SetRadius(r);
  }

  vtkClipPolyData* clipper = vtkClipPolyData::New();
  if(useClip) {
    if(useTubeFilter) {
      clipper->SetInputConnection(tubes->GetOutputPort());
    } else {
      clipper->SetInputConnection(edges->GetOutputPort());
    }
    clipper->SetClipFunction(vtkPost->GetClipPlane());
    clipper->GenerateClippedOutputOn();    
  }

  vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
  if(useClip) {
    mapper->SetInputConnection(clipper->GetOutputPort());
  } else {
    if(useTubeFilter) {
      mapper->SetInputConnection(tubes->GetOutputPort());
    } else {
      mapper->SetInputConnection(edges->GetOutputPort());
    }
  }
  mapper->ScalarVisibilityOff();
  mapper->SetResolveCoincidentTopologyToPolygonOffset();
  // mapper->ImmediateModeRenderingOn();

  actor->GetProperty()->SetLineWidth(lineWidth);
  actor->GetProperty()->SetColor(color.redF(), color.greenF(), color.blueF());
  actor->SetMapper(mapper);

  mapper->Delete();
  clipper->Delete();
  tubes->Delete();
  edges->Delete();
  filter->Delete();
}

void FeatureEdge::removeActorFrom(vtkRenderer* renderer){
  renderer->RemoveActor(actor);
}
void FeatureEdge::addActorTo(vtkRenderer* renderer){
  renderer->AddActor(actor);
}