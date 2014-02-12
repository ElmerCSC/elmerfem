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
 *  ElmerGUI boundarydivision                                                *
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
#include "boundarydivision.h"

BoundaryDivide::BoundaryDivide(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  connect(ui.angleDegreeEdit, SIGNAL(textChanged(const QString&)), this, SLOT(defineAngle(const QString&)));
  connect(ui.divideButton, SIGNAL(clicked()), this, SLOT(divideBoundary()));
  connect(ui.closeButton, SIGNAL(clicked()), this, SLOT(close()));

  ui.angleDegreeEdit->setText("20.0");

  target = TARGET_UNKNOWN;
}

BoundaryDivide::~BoundaryDivide()
{
}

void BoundaryDivide::defineAngle(const QString &qs)
{
  angleDegree = qs;
}

void BoundaryDivide::divideBoundary()
{
  double angle = angleDegree.toDouble();

  if(target == TARGET_SURFACES)
    emit(signalDoDivideSurface(angle));

  if(target == TARGET_EDGES)
    emit(signalDoDivideEdge(angle));
  
  close();
}
