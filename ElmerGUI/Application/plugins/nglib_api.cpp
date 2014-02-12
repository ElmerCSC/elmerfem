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
 *  ELMER/Mesh3D nglib_api                                                   *
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
#include "nglib_api.h"

using namespace std;

NglibAPI::NglibAPI()
{
}


NglibAPI::~NglibAPI()
{
}

void NglibAPI::setDim(int ngDim)
{
  this->ngDim = ngDim;
}

int NglibAPI::getDim() const
{
  return this->ngDim;
}

void NglibAPI::setNgmesh(nglib::Ng_Mesh *ngmesh)
{
  this->ngmesh = ngmesh;
}

void NglibAPI::setNggeom2D(nglib::Ng_Geometry_2D *geom2d)
{
  this->geom2d = geom2d;
}

// Populate elmer's mesh structure:
//-----------------------------------------------------------------------------
mesh_t* NglibAPI::createElmerMeshStructure()
{
  mesh_t *mesh = new mesh_t;

  mesh->setNodes(0);
  mesh->setPoints(0);
  mesh->setEdges(0);
  mesh->setSurfaces(0);
  mesh->setElements(0);

  bool twod = (ngDim == 2) ? true : false;
  
  if(twod) {
    create2D(mesh);
  } else {
    create3D(mesh);
  }

  return mesh;
}

void NglibAPI::create2D(mesh_t *mesh)
{
  Meshutils meshutils;
  
  // Node points:
  //--------------
  mesh->setNodes(nglib::Ng_GetNP_2D(ngmesh));
  mesh->newNodeArray(mesh->getNodes());

  for(int i = 0; i < mesh->getNodes(); i++) {
    node_t *node = mesh->getNode(i);

    double *x = node->getXvec();
    x[0] = 0; x[1] = 0; x[2] = 0;
    nglib::Ng_GetPoint_2D(ngmesh, i + 1, x);
    
    node->setIndex(-1); // default
  }

  // Boundary elements:
  //--------------------
  mesh->setEdges(nglib::Ng_GetNSeg_2D(ngmesh));
  mesh->newEdgeArray(mesh->getEdges());
  
  for(int i = 0; i < mesh->getEdges(); i++) {
    edge_t *edge = mesh->getEdge(i);
    
    edge->setNature(PDE_BOUNDARY);
    edge->setCode(202);
    edge->setNodes(2);
    edge->newNodeIndexes(2);
    edge->setPoints(2);
    edge->newPointIndexes(2);
    
    edge->setPointIndex(0, -1);
    edge->setPointIndex(1, -1);
    
    int matIdx;
    nglib::Ng_GetSegment_2D(ngmesh, i + 1, edge->getNodeIndexes(), &matIdx);

    int bcIdx;
    nglib::EG_GetSegmentBCProperty(ngmesh, geom2d, matIdx-1, &bcIdx);

    edge->setIndex(bcIdx);

    edge->setNodeIndex(0, edge->getNodeIndex(0) - 1);
    edge->setNodeIndex(1, edge->getNodeIndex(1) - 1);
    
    // swap orientation:
    //------------------
    int tmp = edge->getNodeIndex(0);
    edge->setNodeIndex(0, edge->getNodeIndex(1));
    edge->setNodeIndex(1, tmp);
  }
  
  // Elements:
  //-----------
  mesh->setSurfaces(nglib::Ng_GetNE_2D(ngmesh));
  mesh->newSurfaceArray(mesh->getSurfaces()); 
  
  double n[3];
  n[0] = 0; n[1] = 0; n[2] = -1;
  
  for(int i = 0; i < mesh->getSurfaces(); i++) {
    surface_t *surface = mesh->getSurface(i);
    
    surface->setNature(PDE_BULK);
    surface->setCode(303);
    surface->setNodes(3);
    surface->newNodeIndexes(3);
    
    int matIdx;

    nglib::Ng_GetElement_2D(ngmesh, i+1, surface->getNodeIndexes(), &matIdx);

    surface->setIndex(matIdx);
    
    surface->setNodeIndex(0, surface->getNodeIndex(0) - 1);
    surface->setNodeIndex(1, surface->getNodeIndex(1) - 1);
    surface->setNodeIndex(2, surface->getNodeIndex(2) - 1);
    
    surface->setNormalVec(n);    
  }
  
  // Find parents for edge elements:
  //---------------------------------
  meshutils.findEdgeElementParents(mesh);
  
  mesh->setDim(ngDim);
  mesh->setCdim(ngDim);
}

void NglibAPI::create3D(mesh_t *mesh)
{  
  Meshutils meshutils;

  // Node points:
  //--------------
  mesh->setNodes(nglib::Ng_GetNP(ngmesh));
  mesh->newNodeArray(mesh->getNodes());

  for(int i = 0; i < mesh->getNodes(); i++) {
    node_t *node = mesh->getNode(i);
    nglib::Ng_GetPoint(ngmesh, i+1, node->getXvec());
    node->setIndex(-1); // default
  }

  // Boundary elements:
  //--------------------
  mesh->setSurfaces(nglib::Ng_GetNSE(ngmesh));
  mesh->newSurfaceArray(mesh->getSurfaces());
  
  for(int i = 0; i < mesh->getSurfaces(); i++) {
    surface_t *surface = mesh->getSurface(i);
    
    surface->setNature(PDE_BOUNDARY);
    surface->setCode(303);
    surface->setNodes(3);
    surface->newNodeIndexes(3);
    surface->setEdges(3);
    surface->newEdgeIndexes(3);
    
    int face = nglib::EG_GetSurfaceElementBCProperty(ngmesh, i+1);
    
    surface->setIndex(face);
    
    surface->setEdgeIndex(0, -1);
    surface->setEdgeIndex(1, -1);
    surface->setEdgeIndex(2, -1);
    
    nglib::Ng_GetSurfaceElement(ngmesh, i+1, surface->getNodeIndexes());
    
    surface->setNodeIndex(0, surface->getNodeIndex(0) - 1);
    surface->setNodeIndex(1, surface->getNodeIndex(1) - 1);
    surface->setNodeIndex(2, surface->getNodeIndex(2) - 1);
    
    // swap orientation:
    //------------------
    int tmp = surface->getNodeIndex(1);
    surface->setNodeIndex(1, surface->getNodeIndex(2));
    surface->setNodeIndex(2, tmp);
  }

  // Elements:
  //-----------
  mesh->setElements(nglib::Ng_GetNE(ngmesh));
  mesh->newElementArray(mesh->getElements()); 
  
  for(int i = 0; i < mesh->getElements(); i++) {
    element_t *element = mesh->getElement(i);
    
    element->setNature(PDE_BULK);
    element->setCode(504);
    element->setNodes(4);
    element->newNodeIndexes(4);
    
    nglib::Ng_GetVolumeElement(ngmesh, i+1, element->getNodeIndexes());
    
    element->setNodeIndex(0, element->getNodeIndex(0) - 1);
    element->setNodeIndex(1, element->getNodeIndex(1) - 1);
    element->setNodeIndex(2, element->getNodeIndex(2) - 1);
    element->setNodeIndex(3, element->getNodeIndex(3) - 1);
    
    element->setIndex(1); // default (no multibody meshing atm)
  }
  
  // Find parents for surface elements:
  //------------------------------------
  meshutils.findSurfaceElementParents(mesh);
  
  // Find edges for surface elements:
  //----------------------------------
  meshutils.findSurfaceElementEdges(mesh);
  
  // Compute normals for boundary elements:
  //---------------------------------------
  meshutils.findSurfaceElementNormals(mesh);

  mesh->setDim(ngDim);
  mesh->setCdim(ngDim);
}
