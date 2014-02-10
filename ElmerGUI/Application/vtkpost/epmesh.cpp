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
 *  ElmerGUI epmesh                                                          *
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

#include <QString>
#include "epmesh.h"

// EpNode:
//========
EpNode::EpNode()
{
  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;
}

EpNode::~EpNode()
{
}

// EpElement:
//===========
EpElement::EpElement()
{
  groupName = "";
  code = 0;
  indexes = 0;
  index = NULL;
}

EpElement::~EpElement()
{
  delete [] index;
}

// EpMesh:
//=========
EpMesh::EpMesh()
{
  epNodes = 0;
  epNode = NULL;

  epElements = 0;
  epElement = NULL;
}

EpMesh::~EpMesh()
{
  delete [] epNode;
  delete [] epElement;
}

// ScalarField:
//==============
ScalarField::ScalarField()
{
  name = "";
  values = 0;
  value = NULL;
  minVal = +9.9e99;
  maxVal = -9.9e99;
}

ScalarField::~ScalarField()
{
  delete [] value;
}
