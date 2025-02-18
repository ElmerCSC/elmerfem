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
 *  ElmerGUI axes                                                            *
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
#include "axes.h"

#include <vtkAxes.h>
#include <vtkTubeFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkFollower.h>
#include <vtkVectorText.h>
#include <vtkProperty.h>

using namespace std;

Axes::Axes(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  setWindowTitle("Coordinate axes");
  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
}

Axes::~Axes()
{
}

void Axes::draw(VtkPost* vtkPost)
{
  double scl = vtkPost->GetLength() / 8.0;

  vtkAxes* axes = vtkAxes::New();
  axes->SetOrigin(0, 0, 0);
  axes->SetScaleFactor(scl);

  vtkTubeFilter* axesTubes = vtkTubeFilter::New();
  axesTubes->SetInputConnection(axes->GetOutputPort());
  axesTubes->SetRadius(axes->GetScaleFactor() / 33.0);
  axesTubes->SetNumberOfSides(20);

  vtkPolyDataMapper* axesMapper = vtkPolyDataMapper::New();
  axesMapper->SetInputConnection(axesTubes->GetOutputPort());

  vtkPost->GetAxesActor()->SetMapper(axesMapper);

  vtkVectorText* XText = vtkVectorText::New();
  vtkVectorText* YText = vtkVectorText::New();
  vtkVectorText* ZText = vtkVectorText::New();

  XText->SetText("X");
  YText->SetText("Y");
  ZText->SetText("Z");
  
  vtkPolyDataMapper* XTextPolyDataMapper = vtkPolyDataMapper::New();
  vtkPolyDataMapper* YTextPolyDataMapper = vtkPolyDataMapper::New();
  vtkPolyDataMapper* ZTextPolyDataMapper = vtkPolyDataMapper::New();

  XTextPolyDataMapper->SetInputConnection(XText->GetOutputPort());
  YTextPolyDataMapper->SetInputConnection(YText->GetOutputPort());
  ZTextPolyDataMapper->SetInputConnection(ZText->GetOutputPort());

  vtkPost->GetAxesXTextActor()->SetMapper(XTextPolyDataMapper);
  vtkPost->GetAxesYTextActor()->SetMapper(YTextPolyDataMapper);
  vtkPost->GetAxesZTextActor()->SetMapper(ZTextPolyDataMapper);

  scl = axes->GetScaleFactor() / 5.0;

  vtkPost->GetAxesXTextActor()->SetScale(scl, scl, scl);
  vtkPost->GetAxesYTextActor()->SetScale(scl, scl, scl);
  vtkPost->GetAxesZTextActor()->SetScale(scl, scl, scl);

  scl = axes->GetScaleFactor();

  vtkPost->GetAxesXTextActor()->SetPosition(scl, 0.0, 0.0);
  vtkPost->GetAxesYTextActor()->SetPosition(0.0, scl, 0.0);
  vtkPost->GetAxesZTextActor()->SetPosition(0.0, 0.0, scl);
  
  vtkPost->GetAxesXTextActor()->GetProperty()->SetColor(0, 0, 0);
  vtkPost->GetAxesYTextActor()->GetProperty()->SetColor(0, 0, 0);
  vtkPost->GetAxesZTextActor()->GetProperty()->SetColor(0, 0, 0);

  // Clean up:
  //----------
  XTextPolyDataMapper->Delete();
  YTextPolyDataMapper->Delete();
  ZTextPolyDataMapper->Delete();
  XText->Delete();
  YText->Delete();
  ZText->Delete();
  axesTubes->Delete();
  axesMapper->Delete();
  axes->Delete();  
}
