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
 *  ELMER/Mesh3D tetlib_api                                                  *
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
#include "tetlib_api.h"
#include <iostream>
using namespace std;

TetlibAPI::TetlibAPI()
{
}

TetlibAPI::~TetlibAPI()
{
}

bool TetlibAPI::loadTetlib()
{
  cout << "Load tetgen plugin... ";

  libtet = new QLibrary("tetplugin");

  if(!libtet->load()) {
    cout << "not found" << endl;
    cout << "Tetgen functionality unavailable" << endl;
    cout.flush();
    return false;
  }

  cout << "done" << endl;
  cout.flush();

  ptetgenio = (tetgenio_t)libtet->resolve("CreateObjectOfTetgenio");

  if(!ptetgenio) {
    cout << "Unable to resolve 'CreateObjectOfTetgenio'" << endl;
    cout.flush();
    return false;
  }

  in = (ptetgenio)();
  out = (ptetgenio)(); 

  delegate_tetrahedralize = (delegate_tetrahedralize_t)libtet->resolve("delegate_tetrahedralize");

  if(!delegate_tetrahedralize) {
    cout << "Unable to resolve 'delegate_tetrahedralize'" << endl;
    cout.flush();
    return false;
  }

  return true;
}

// Populate Elmer's mesh structure:
//-----------------------------------------------------------------------------
mesh_t *TetlibAPI::createElmerMeshStructure()
{
  Helpers helpers;
  Meshutils meshutils;

  // Create new mesh structure:
  mesh_t *mesh = new mesh_t;

  mesh->setNodes(0);
  mesh->setPoints(0);
  mesh->setEdges(0);
  mesh->setSurfaces(0);
  mesh->setElements(0);
  
  // Nodes:
  mesh->setNodes(out->numberofpoints);
  mesh->newNodeArray(mesh->getNodes());

  REAL *pointlist = out->pointlist;

  for(int i=0; i < mesh->getNodes(); i++) {
    node_t *node = mesh->getNode(i);
    
    node->setX(0, *pointlist++);
    node->setX(1, *pointlist++);
    node->setX(2, *pointlist++);

    node->setIndex(-1); // default
  }

  // Elements:
  mesh->setElements(out->numberoftetrahedra);
  mesh->newElementArray(mesh->getElements());

  int *tetrahedronlist = out->tetrahedronlist;
  REAL *attribute = out->tetrahedronattributelist;
  int na = out->numberoftetrahedronattributes;
  
  for(int i=0; i< mesh->getElements(); i++) {
    element_t *element = mesh->getElement(i);

    element->setNature(PDE_BULK);
    element->setCode(504);
    element->setNodes(4);
    element->newNodeIndexes(4);
    
    element->setNodeIndex(0, (*tetrahedronlist++) - out->firstnumber);
    element->setNodeIndex(1, (*tetrahedronlist++) - out->firstnumber);
    element->setNodeIndex(2, (*tetrahedronlist++) - out->firstnumber);
    element->setNodeIndex(3, (*tetrahedronlist++) - out->firstnumber);
    
    element->setIndex(1); // default
    // must have "A" in control string:
    if(out->tetrahedronattributelist != (REAL*)NULL) 
      element->setIndex((int)attribute[na*(i+1)-1]);
  }
  
  // Boundary elements:
  mesh->setSurfaces(out->numberoftrifaces);
  mesh->newSurfaceArray(mesh->getSurfaces());

  int *trifacelist = out->trifacelist;
  int *face2tetlist = out->face2tetlist;

  for(int i=0; i < mesh->getSurfaces(); i++) {
    surface_t *surface = mesh->getSurface(i);

    surface->setNature(PDE_BOUNDARY);
    surface->setCode(303);
    surface->setNodes(3);
    surface->newNodeIndexes(3);
    surface->setEdges(3);
    surface->newEdgeIndexes(3);

    surface->setElements(2);
    surface->newElementIndexes(2);

    surface->setIndex(1); // default
    if(out->trifacemarkerlist != (int*)NULL)
      surface->setIndex(out->trifacemarkerlist[i]);

    surface->setEdgeIndex(0, -1);
    surface->setEdgeIndex(1, -1);
    surface->setEdgeIndex(2, -1);
    
    surface->setElementIndex(0, -1);
    surface->setElementIndex(1, -1);

    // must have "nn" in control string:
    if(out->face2tetlist != (int*)NULL) {
      surface->setElementIndex(0, (*face2tetlist++) - out->firstnumber);
      surface->setElementIndex(1, (*face2tetlist++) - out->firstnumber);
    }

    int u = (*trifacelist++) - out->firstnumber;
    int v = (*trifacelist++) - out->firstnumber;
    int w = (*trifacelist++) - out->firstnumber;

    surface->setNodeIndex(0, u);
    surface->setNodeIndex(1, v);
    surface->setNodeIndex(2, w);
  }

  // Edges:
  meshutils.findSurfaceElementEdges(mesh);
  meshutils.findSurfaceElementNormals(mesh);

  // Points:
  mesh->setPoints(0);
  // mesh->point == NULL;
  
  mesh->setDim(3);
  mesh->setCdim(3);

  return mesh;
}
