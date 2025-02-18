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
 *  ElmerGUI meshutils                                                       *
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

#ifndef MESHUTILS_H
#define MESHUTILS_H

#include <math.h>
#include "meshtype.h"
#include "helpers.h"

#define SHARPEDGE     0
#define UNKNOWN      -1
#define MYPI 3.14159265

class Meshutils
{
 public:
  Meshutils();
  ~Meshutils();

  void clearMesh(mesh_t*);
  void findEdgeElementPoints(mesh_t*);
  void findSurfaceElements(mesh_t*);
  void findSurfaceElementEdges(mesh_t*);
  void findSurfaceElementParents(mesh_t*);
  void findEdgeElementParents(mesh_t*);
  void findSurfaceElementNormals(mesh_t*);
  void findSharpEdges(mesh_t*, double);
  void findSharpPoints(mesh_t*, double);
  int divideEdgeBySharpPoints(mesh_t*);
  int divideSurfaceBySharpEdges(mesh_t*);
  void sort_index(int n, double *a, int *b);
  void increaseElementOrder(mesh_t*);
  void decreaseElementOrder(mesh_t*);
  int cleanHangingSharpEdges(mesh_t*);
};
#endif // #ifndef MESHUTILS_H
