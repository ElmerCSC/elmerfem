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
 *  ElmerGUI helpers                                                         *
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

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "helpers.h"
#include <QMatrix4x4>

Helpers :: Helpers()
{
}

Helpers :: ~Helpers()
{
}

//====================================================================
//                             Normalize
//====================================================================

void Helpers::normalize(qreal *a)
{
  double b;

  b = vlen(a);
  a[0] /= b;
  a[1] /= b;
  a[2] /= b;
}

//====================================================================
//                              Length
//====================================================================

qreal Helpers::vlen(qreal *a)
{
  return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

//====================================================================
//                           Cross product
//====================================================================

void Helpers::crossProduct(qreal *a, qreal *b, qreal *c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

//====================================================================
//            Invert 4x4 matrix (for visualiztion only)
//====================================================================
void Helpers::invertMatrix(const qreal *a, qreal *inva)
{
  QMatrix4x4 matrix(a);

  bool ok(true);

  QMatrix4x4 inverse(matrix.transposed().inverted(&ok));

  if(!ok) fprintf(stderr, "Error: Singular 4x4 matrix\n");

  for(int i = 0; i < 16; ++i)
    inva[i] = double(inverse.data()[i]);
}
