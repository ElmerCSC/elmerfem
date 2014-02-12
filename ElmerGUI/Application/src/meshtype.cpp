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
 *  ElmerGUI mesh_t (Elmer mesh structure)                                   *
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
#include <fstream>
#include "meshtype.h"
using namespace std;

// node_t
//-----------------------------------------------------------------------------
node_t::node_t()
{
}

node_t::~node_t()
{
}

void node_t::setXvec(double* y)
{
  this->x[0] = y[0];
  this->x[1] = y[1];
  this->x[2] = y[2];
}

double* node_t::getXvec()
{
  return &this->x[0];
}

void node_t::setX(int n, double y)
{
  this->x[n] = y;
}

double node_t::getX(int n) const
{
  return this->x[n];
}

void node_t::setIndex(int n)
{
  this->index = n;

}

int node_t::getIndex() const
{
  return this->index;
}

// element_t
//-----------------------------------------------------------------------------
element_t::element_t()
{
}

element_t::~element_t()
{
}

void element_t::setNature(int n)
{
  this->nature = n;
}

int element_t::getNature() const
{
  return this->nature;
}

void element_t::setCode(int n)
{
  this->code = n;
}

int element_t::getCode() const
{
  return this->code;
}

void element_t::setNodes(int n)
{
  this->nodes = n;
}

int element_t::getNodes() const
{
  return this->nodes;
}

void element_t::setIndex(int n)
{
  this->index = n;
}

int element_t::getIndex() const
{
  return this->index;
}

void element_t::setSelected(int n)
{
  this->selected = n;
}

int element_t::getSelected() const
{
  return this->selected;
}

int element_t::getNodeIndex(int n) const
{
  return this->node[n];
}

void element_t::setNodeIndex(int m, int n)
{
  this->node[m] = n;
}

void element_t::newNodeIndexes(int n)
{
  this->node = new int[n];
}

void element_t::deleteNodeIndexes()
{
  delete [] this->node;
}

int* element_t::getNodeIndexes() const
{
  return this->node;
}

// point_t
//-----------------------------------------------------------------------------
point_t::point_t()
{
}

point_t::~point_t()
{
}

void point_t::setSharp(bool b)
{
  this->sharp_point = b;
}

bool point_t::isSharp() const
{
  return this->sharp_point;
}

void point_t::setEdges(int n)
{
  this->edges = n;
}

int point_t::getEdges() const
{
  return this->edges;
}

int point_t::getEdgeIndex(int n) const
{
  return this->edge[n];
}

void point_t::setEdgeIndex(int m, int n)
{
  this->edge[m] = n;
}

void point_t::newEdgeIndexes(int n)
{
  this->edge = new int[n];
}

void point_t::deleteEdgeIndexes()
{
  delete [] this->edge;
}

// edge_t
//-----------------------------------------------------------------------------
edge_t::edge_t()
{
}

edge_t::~edge_t()
{
}

void edge_t::setSharp(bool b)
{
  this->sharp_edge = b;
}

bool edge_t::isSharp() const
{
  return this->sharp_edge;
}

void edge_t::setPoints(int n)
{
  this->points = n;
}

int edge_t::getPoints() const
{
  return this->points;
}

void edge_t::setPointIndex(int m, int n)
{
  this->point[m] = n;
}

int edge_t::getPointIndex(int n) const
{
  return this->point[n];
}

void edge_t::newPointIndexes(int n)
{
  this->point = new int[n];
}

void edge_t::deletePointIndexes()
{
  delete [] this->point;
}

void edge_t::setSurfaces(int n)
{
  this->surfaces = n;
}

int edge_t::getSurfaces() const
{
  return this->surfaces;
}

void edge_t::setSurfaceIndex(int m, int n)
{
  this->surface[m] = n;
}

int edge_t::getSurfaceIndex(int n) const
{
  return this->surface[n];
}

void edge_t::newSurfaceIndexes(int n)
{
  this->surface = new int[n];
}

void edge_t::deleteSurfaceIndexes()
{
  delete [] this->surface;
}

// surface_t
//-----------------------------------------------------------------------------
surface_t::surface_t()
{
}

surface_t::~surface_t()
{
}

void surface_t::setEdges(int n)
{
  this->edges = n;
}

int surface_t::getEdges() const
{
  return this->edges;
}

void surface_t::setEdgeIndex(int m, int n)
{
  this->edge[m] = n;
}

int surface_t::getEdgeIndex(int n) const
{
  return this->edge[n];
}

void surface_t::newEdgeIndexes(int n)
{
  this->edge = new int[n];
}

void surface_t::deleteEdgeIndexes()
{
  delete [] this->edge;
}

void surface_t::setElements(int n)
{
  this->elements = n;
}

int surface_t::getElements() const
{
  return this->elements;
}

void surface_t::setElementIndex(int m, int n)
{
  this->element[m] = n;
}

int surface_t::getElementIndex(int n) const
{
  return this->element[n];
}

void surface_t::newElementIndexes(int n)
{
  this->element = new int[n];
}

void surface_t::deleteElementIndexes()
{
  delete [] this->element;
}

double* surface_t::getNormalVec()
{
  return &this->normal[0];
}

void surface_t::setNormalVec(double* d)
{
  this->normal[0] = d[0];
  this->normal[1] = d[1];
  this->normal[2] = d[2];
}

double surface_t::getNormal(int n) const
{
  return this->normal[n];
}

void surface_t::setNormal(int n, double d)
{
  this->normal[n] = d;
}

void surface_t::setVertexNormalVec(int n, double* d)
{
  this->vertex_normals[n][0] = d[0];
  this->vertex_normals[n][1] = d[1];
  this->vertex_normals[n][2] = d[2];
}

void surface_t::addVertexNormalVec(int n, double* d)
{
  this->vertex_normals[n][0] += d[0];
  this->vertex_normals[n][1] += d[1];
  this->vertex_normals[n][2] += d[2];
}

void surface_t::subVertexNormalVec(int n, double* d)
{
  this->vertex_normals[n][0] -= d[0];
  this->vertex_normals[n][1] -= d[1];
  this->vertex_normals[n][2] -= d[2];
}

double* surface_t::getVertexNormalVec(int n)
{
  return &this->vertex_normals[n][0];
}

// mesh_t
//-----------------------------------------------------------------------------
mesh_t::mesh_t()
{
  this->setDefaults();
}

mesh_t::~mesh_t()
{
}

bool mesh_t::isUndefined() const
{
  if((cdim < 0) || (dim < 0) || (nodes < 1))
    return true;

  return false;
}

void mesh_t::clear()
{
  delete [] element;
  delete [] surface;
  delete [] edge;
  delete [] point;
  delete [] node;
  
  setDefaults();
}

void mesh_t::setDefaults()
{
  cdim = -1;
  dim = -1;
  nodes = 0;
  node = 0;
  points = 0;
  point = 0;
  edges = 0;
  edge = 0;
  surfaces = 0;
  surface = 0;
  elements = 0;
  element = 0;
}

// Load Elmer mesh files and populate mesh structures
//---------------------------------------------------------------------------
bool mesh_t::load(char* dirName)
{
  char fileName[1024];
  ifstream mesh_header;
  ifstream mesh_nodes;
  ifstream mesh_elements;
  ifstream mesh_boundary;

  // Header:
  //--------
  sprintf(fileName, "%s/mesh.header", dirName);
  mesh_header.open(fileName);

  if(!mesh_header.is_open()) {
    cout << "Mesh: load: unable to open " << fileName << endl;
    return false;
  }

  int nodes, elements, boundaryelements, types, type, ntype;

  mesh_header >> nodes >> elements >> boundaryelements;
  mesh_header >> types;

  int elements_zero_d = 0;
  int elements_one_d = 0;
  int elements_two_d = 0;
  int elements_three_d = 0;
  
  for(int i = 0; i < types; i++) {
    mesh_header >> type >> ntype;
    
    switch(type/100) {
    case 1:
      elements_zero_d += ntype;
      break;
    case 2:
      elements_one_d += ntype;
      break;
    case 3:
    case 4:
      elements_two_d += ntype;
      break;
    case 5:
    case 6:
    case 7:
    case 8:
      elements_three_d += ntype;
      break;
    default:
      cout << "Unknown element family (possibly not implamented)" << endl;
      cout.flush();
      return false;
      // exit(0);
    }
  }
  
  cout << "Summary:" << endl;
  cout << "Nodes: " << nodes << endl;
  cout << "Point elements: " << elements_zero_d << endl;
  cout << "Edge elements: " << elements_one_d << endl;
  cout << "Surface elements: " << elements_two_d << endl;
  cout << "Volume elements: " << elements_three_d << endl;
  cout.flush();

  // Set mesh dimension:
  this->dim = 0;

  if(elements_one_d > 0)
    this->dim = 1;

  if(elements_two_d > 0)
    this->dim = 2;

  if(elements_three_d > 0)
    this->dim = 3;

  this->nodes = nodes;
  node = new node_t[nodes];

  this->points = elements_zero_d;
  point = new point_t[this->points];

  this->edges = elements_one_d;
  edge = new edge_t[this->edges];

  this->surfaces = elements_two_d;
  surface = new surface_t[this->surfaces];

  this->elements = elements_three_d;
  element = new element_t[this->elements];

  mesh_header.close();

  // Nodes:
  //-------
  sprintf(fileName, "%s/mesh.nodes", dirName);
  mesh_nodes.open(fileName);

  if(!mesh_nodes.is_open()) {
    cout << "Mesh: load: unable to open " << fileName << endl;
    return false;
  }

  int number, index;
  double x, y, z;

  for(int i = 0; i < nodes; i++) {
    node_t *node = &this->node[i];
    mesh_nodes >> number >> index >> x >> y >> z;
    node->setX(0, x);
    node->setX(1, y);
    node->setX(2, z);
    node->setIndex(index);
  }

  mesh_nodes.close();  

  // Elements:
  //----------
  sprintf(fileName, "%s/mesh.elements", dirName);
  mesh_elements.open(fileName);

  if(!mesh_elements.is_open()) {
    cout << "Mesh: load: unable to open " << fileName << endl;
    return false;
  }

  int current_point = 0;
  int current_edge = 0;
  int current_surface = 0;
  int current_element = 0;

  point_t *point = NULL;
  edge_t *edge = NULL;
  surface_t *surface = NULL;
  element_t *element = NULL;

  for(int i = 0; i < elements; i++) {
    mesh_elements >> number >> index >> type;

    switch(type/100) {
    case 1:
      point = &this->point[current_point++];
      point->setNature(PDE_BULK);
      point->setIndex(index);
      point->setCode(type);
      point->setNodes(point->getCode() % 100);
      point->newNodeIndexes(point->getNodes());
      for(int j = 0; j < point->getNodes(); j++) {
	int k;
	mesh_elements >> k;
	point->setNodeIndex(j, k-1);
      }
      point->setEdges(2);
      point->newEdgeIndexes(point->getEdges());
      point->setEdgeIndex(0, -1);
      point->setEdgeIndex(1, -1);
      break;

    case 2:
      edge = &this->edge[current_edge++];
      edge->setNature(PDE_BULK);
      edge->setIndex(index);
      edge->setCode(type);
      edge->setNodes(edge->getCode() % 100);
      edge->newNodeIndexes(edge->getNodes());
      for(int j = 0; j < edge->getNodes(); j++) {
	int k;	
	mesh_elements >> k;
	edge->setNodeIndex(j, k-1);
      }
      edge->setSurfaces(0);
      edge->newSurfaceIndexes(edge->getSurfaces());
      edge->setSurfaceIndex(0, -1);
      edge->setSurfaceIndex(1, -1);

      break;

    case 3:
    case 4:
      surface = &this->surface[current_surface++];
      surface->setNature(PDE_BULK);
      surface->setIndex(index);
      surface->setCode(type);
      surface->setNodes(surface->getCode() % 100);
      surface->newNodeIndexes(surface->getNodes());
      for(int j = 0; j < surface->getNodes(); j++) {
	int k;
	mesh_elements >> k;
	surface->setNodeIndex(j, k-1);
      }      
      surface->setEdges((int)(surface->getCode() / 100));
      surface->newEdgeIndexes(surface->getEdges());
      for(int j = 0; j < surface->getEdges(); j++)
	surface->setEdgeIndex(j, -1);
      surface->setElements(2);
      surface->newElementIndexes(surface->getElements());
      surface->setElementIndex(0, -1);
      surface->setElementIndex(1, -1);

      break;

    case 5:
    case 6:
    case 7:
    case 8:
      element = &this->element[current_element++];
      element->setNature(PDE_BULK);
      element->setIndex(index);
      element->setCode(type);
      element->setNodes(element->getCode() % 100);
      element->newNodeIndexes(element->getNodes());
      for(int j = 0; j < element->getNodes(); j++) {
	int k;
	mesh_elements >> k;
	element->setNodeIndex(j, k-1);
      }
      break;

    default:
      cout << "Unknown element type (possibly not implemented" << endl;
      cout.flush();
      return false;
      // exit(0);
    }
  }

  mesh_elements.close();

  // Boundary elements:
  //-------------------
  sprintf(fileName, "%s/mesh.boundary", dirName);
  mesh_boundary.open(fileName);

  if(!mesh_boundary.is_open()) {
    cout << "Mesh: load: unable to open " << fileName << endl;
    return false;
  }

  int parent0, parent1;

  for(int i = 0; i < boundaryelements; i++) {
    mesh_boundary >> number >> index >> parent0 >> parent1 >> type;

    switch(type/100) {
    case 1:
      point = &this->point[current_point++];
      point->setNature(PDE_BOUNDARY);
      point->setIndex(index);
      point->setEdges(2);
      point->newEdgeIndexes(point->getEdges());
      point->setEdgeIndex(0, parent0 - 1);
      point->setEdgeIndex(1, parent0 - 1);
      point->setCode(type);
      point->setNodes(point->getCode() % 100);
      point->newNodeIndexes(point->getNodes());
      for(int j = 0; j < point->getNodes(); j++) {
	int k;
	mesh_boundary >> k;
	point->setNodeIndex(j, k-1);
      }
      break;

    case 2:
      edge = &this->edge[current_edge++];
      edge->setNature(PDE_BOUNDARY);
      edge->setIndex(index);
      edge->setSurfaces(2);
      edge->newSurfaceIndexes(edge->getSurfaces());
      edge->setSurfaceIndex(0, parent0 - 1);
      edge->setSurfaceIndex(1, parent1 - 1);
      edge->setCode(type);
      edge->setNodes(edge->getCode() % 100);
      edge->newNodeIndexes(edge->getNodes());
      for(int j = 0; j < edge->getNodes(); j++) {
	int k;
	mesh_boundary >> k;
	edge->setNodeIndex(j, k-1);
      }

      break;

    case 3:
    case 4:
      surface = &this->surface[current_surface++];
      surface->setNature(PDE_BOUNDARY);
      surface->setIndex(index);
      surface->setElements(2);
      surface->newElementIndexes(surface->getElements());
      surface->setElementIndex(0, parent0 - 1);
      surface->setElementIndex(1, parent1 - 1);
      surface->setCode(type);
      surface->setNodes(surface->getCode() % 100);
      surface->newNodeIndexes(surface->getNodes());
      for(int j = 0; j < surface->getNodes(); j++) {
	int k;
	mesh_boundary >> k;
	surface->setNodeIndex(j, k-1);
      }
      surface->setEdges((int)(surface->getCode() / 100));
      surface->newEdgeIndexes(surface->getEdges());
      for(int j = 0; j < surface->getEdges(); j++)
	surface->setEdgeIndex(j, -1);
      
      break;

    case 5:
    case 6:
    case 7:
    case 8:
      // these can't be boundary elements
      break;

    default:
      break;
    }
  }

  mesh_boundary.close();

  this->boundingBox();

  return true;
}

// Save Elmer mesh files and populate mesh structures
//---------------------------------------------------------------------------
bool mesh_t::save(char *dirName)
{
  char fileName[1024];
  ofstream mesh_header;
  ofstream mesh_nodes;
  ofstream mesh_elements;
  ofstream mesh_boundary;

  // Elmer's elements codes are smaller than 1000
  int maxcode = 1000;
  int *bulk_by_type = new int[maxcode];
  int *boundary_by_type = new int[maxcode];
  
  for(int i = 0; i < maxcode; i++) {
    bulk_by_type[i] = 0;
    boundary_by_type[i] = 0;
  }
  
  for(int i = 0; i < elements; i++) {
    element_t *e = &element[i];
      
    if(e->getNature() == PDE_BULK) 
      bulk_by_type[e->getCode()]++;
    
    if(e->getNature() == PDE_BOUNDARY)
      boundary_by_type[e->getCode()]++;
  }
	
  for(int i = 0; i < surfaces; i++) {
    surface_t *s = &surface[i];

    if(s->getNature() == PDE_BULK)
      bulk_by_type[s->getCode()]++;

    if(s->getNature() == PDE_BOUNDARY)
      boundary_by_type[s->getCode()]++;
  }

  for(int i = 0; i < edges; i++) {
    edge_t *e = &edge[i];

    if(e->getNature() == PDE_BULK)
      bulk_by_type[e->getCode()]++;

    if(e->getNature() == PDE_BOUNDARY)
      boundary_by_type[e->getCode()]++;
  }

  for(int i = 0; i < points; i++) {
    point_t *p = &point[i];

    if(p->getNature() == PDE_BULK)
      bulk_by_type[p->getCode()]++;

    if(p->getNature() == PDE_BOUNDARY)
      boundary_by_type[p->getCode()]++;
  }

  int bulk_elements = 0;
  int boundary_elements = 0;
  int element_types = 0;
  
  for(int i = 0; i < maxcode; i++) {
    bulk_elements += bulk_by_type[i];
    boundary_elements += boundary_by_type[i];
    
    if((bulk_by_type[i] > 0) || (boundary_by_type[i] > 0))
      element_types++;
  }
  
  // Header:
  //---------
  sprintf(fileName, "%s/mesh.header", dirName);
  
  mesh_header.open(fileName);
  
  if(!mesh_header.is_open()) {
    cout << "Unable to open " << fileName << endl;
    return false;
  }
  
  cout << "Saving " << nodes << " nodes" << endl;
  cout << "Saving " << bulk_elements << " elements" << endl;
  cout << "Saving " << boundary_elements << " boundary elements" << endl;
  cout.flush();
  
  mesh_header << nodes << " ";
  mesh_header << bulk_elements << " ";
  mesh_header << boundary_elements << endl;
  mesh_header << element_types << endl;
  
  for(int i = 0; i < maxcode; i++) {
    int j = bulk_by_type[i] + boundary_by_type[i];
    if(j > 0) 
      mesh_header << i << " " << j << endl;
  }
  
  mesh_header.close();

  // Nodes:
  //--------
  sprintf(fileName, "%s/mesh.nodes", dirName);

  mesh_nodes.open(fileName);
  
  if(!mesh_nodes.is_open()) {
    cout << "Unable to open " << fileName << endl;
    return false;
  }

  for(int i = 0; i < this->nodes; i++) {
    node_t *node = &this->node[i];

    int ind = node->getIndex();

    mesh_nodes << i+1 << " " << ind << " ";
    mesh_nodes << node->getX(0) << " ";
    mesh_nodes << node->getX(1) << " ";
    mesh_nodes << node->getX(2) << endl;
  }

  mesh_nodes.close();

  
  // Elements:
  //----------
  sprintf(fileName, "%s/mesh.elements", dirName);

  mesh_elements.open(fileName);

  if(!mesh_elements.is_open()) {
    cout << "Unable to open " << fileName << endl;
    return false;
  }

  int current = 0;

  for(int i = 0; i < this->elements; i++) {
    element_t *e = &this->element[i];

    int ind = e->getIndex();

    if(ind < 1)
      ind = 1;

    if(e->getNature() == PDE_BULK) {
      mesh_elements << ++current << " ";
      mesh_elements << ind << " ";
      mesh_elements << e->getCode() << " ";

      for(int j = 0; j < e->getNodes(); j++) 
	mesh_elements << e->getNodeIndex(j) + 1 << " ";

      mesh_elements << endl;
    }
  }

  for(int i = 0; i < this->surfaces; i++) {
    surface_t *s = &this->surface[i];

    int ind = s->getIndex();

    if(ind < 1)
      ind = 1;

    if(s->getNature() == PDE_BULK) {
      mesh_elements << ++current << " ";
      mesh_elements << ind << " ";
      mesh_elements << s->getCode() << " ";

      for(int j = 0; j < s->getNodes(); j++) 
	mesh_elements << s->getNodeIndex(j) + 1 << " ";

      mesh_elements << endl;
    }
  }

  for(int i = 0; i < this->edges; i++) {
    edge_t *e = &this->edge[i];

    int ind = e->getIndex();

    if(ind < 1)
      ind = 1;

    if(e->getNature() == PDE_BULK) {
      mesh_elements << ++current << " ";
      mesh_elements << ind << " ";
      mesh_elements << e->getCode() << " ";

      for(int j = 0; j < e->getNodes(); j++)
	mesh_elements << e->getNodeIndex(j) + 1 << " ";

      mesh_elements << endl;
    }
  }

  for(int i = 0; i < this->points; i++) {
    point_t *p = &this->point[i];

    int ind = p->getIndex();

    if(ind < 1)
      ind = 1;

    if(p->getNature() == PDE_BULK) {
      mesh_elements << ++current << " ";
      mesh_elements << ind << " ";
      mesh_elements << p->getCode() << " ";

      for(int j = 0; j < p->getNodes(); j++)
	mesh_elements << p->getNodeIndex(j) + 1 << " ";

      mesh_elements << endl;
    }
  }

  mesh_elements.close();

  // Boundary elements:
  //-------------------
  sprintf(fileName, "%s/mesh.boundary", dirName);

  mesh_boundary.open(fileName);

  if(!mesh_boundary.is_open()) {
    cout << "Unable to open " << fileName << endl;
    return false;
  }

  current = 0;

  for(int i = 0; i < this->surfaces; i++) {
    surface_t *s = &this->surface[i];

    if(s->getNature() == PDE_BULK)
      continue;

    int e0 = s->getElementIndex(0) + 1;
    int e1 = s->getElementIndex(1) + 1;

    if(e0 < 0)
      e0 = 0;

    if(e1 < 0)
      e1 = 0;

    int ind = s->getIndex();

    if(ind < 1)
      ind = 1;

    if(s->getNature() == PDE_BOUNDARY) {
      mesh_boundary << ++current << " ";
      mesh_boundary << ind << " ";
      mesh_boundary << e0 << " " << e1 << " ";
      mesh_boundary << s->getCode() << " ";

      for(int j = 0; j < s->getNodes(); j++) 
	mesh_boundary << s->getNodeIndex(j) + 1 << " ";

      mesh_boundary << endl;
    }
  }


  for(int i = 0; i < this->edges; i++) {
    edge_t *e = &this->edge[i];

    if(e->getNature() == PDE_BULK)
      continue;

    int s0 = e->getSurfaceIndex(0) + 1;
    int s1 = e->getSurfaceIndex(1) + 1;

    if(s0 < 0)
      s0 = 0;

    if(s1 < 0)
      s1 = 0;

    int ind = e->getIndex();

    if(ind < 1)
      ind = 1;

    if(e->getNature() == PDE_BOUNDARY) {
      mesh_boundary << ++current << " ";
      mesh_boundary << ind << " ";
      mesh_boundary << s0 << " " << s1 << " ";
      mesh_boundary << e->getCode() << " ";

      for(int j = 0; j < e->getNodes(); j++) 
	mesh_boundary << e->getNodeIndex(j) + 1 << " ";

      mesh_boundary << endl;
    }
  }

  for(int i = 0; i < this->points; i++) {
    point_t *p = &this->point[i];

    int e0 = p->getEdgeIndex(0) + 1;
    int e1 = p->getEdgeIndex(1) + 1;

    if(e0 < 0)
      e0 = 0;

    if(e1 < 0)
      e1 = 0;

    int ind = p->getIndex();

    if(ind < 1)
      ind = 1;

    if(p->getNature() == PDE_BOUNDARY) {
      mesh_boundary << ++current << " ";
      mesh_boundary << ind << " ";
      mesh_boundary << e0 << " " << e1 << " ";
      mesh_boundary << p->getCode() << " ";

      for(int j = 0; j < p->getNodes(); j++) 
	mesh_boundary << p->getNodeIndex(j) + 1 << " ";

      mesh_boundary << endl;
    }
  }

  mesh_boundary.close();

  delete [] bulk_by_type;
  delete [] boundary_by_type;

  return true;
}

// Bounding box...
//-----------------------------------------------------------------------------
double* mesh_t::boundingBox()
{
  double *result = new double[10];

  double xmin = +9e9;
  double xmax = -9e9;

  double ymin = +9e9;
  double ymax = -9e9;

  double zmin = +9e9;
  double zmax = -9e9;

  for(int i=0; i < this->nodes; i++) {
    node_t *node = &this->node[i];
    
    if(node->getX(0) > xmax) 
      xmax = node->getX(0);
    
    if(node->getX(0) < xmin) 
      xmin = node->getX(0);
    
    if(node->getX(1) > ymax) 
      ymax = node->getX(1);

    if(node->getX(1) < ymin) 
      ymin = node->getX(1);

    if(node->getX(2) > zmax) 
      zmax = node->getX(2);

    if(node->getX(2) < zmin) 
      zmin = node->getX(2);
  }
  
  double xmid = (xmin + xmax)/2.0;
  double ymid = (ymin + ymax)/2.0;
  double zmid = (zmin + zmax)/2.0;

  double xlen = (xmax - xmin)/2.0;
  double ylen = (ymax - ymin)/2.0;
  double zlen = (zmax - zmin)/2.0;

  double s = xlen;

  if(ylen > s)
    s = ylen;
  
  if(zlen > s)
    s = zlen;
  
  s *= 1.1;

  bool sx = xmin==xmax;
  bool sy = ymin==ymax;
  bool sz = zmin==zmax;

  this->cdim = 3;

  if((sz && sy) || (sz && sx) || (sx && sy))
    this->cdim = 1;

  else if(sz || sy || sx)
    this->cdim = 2;

  result[0] = xmin;
  result[1] = xmax;
  result[2] = ymin;
  result[3] = ymax;
  result[4] = zmin;
  result[5] = zmax;

  result[6] = xmid;
  result[7] = ymid;
  result[8] = zmid;

  result[9] = s;

  return result;
}

void mesh_t::setCdim(int n)
{
  this->cdim = n;
}

int mesh_t::getCdim() const
{
  return this->cdim;
}

void mesh_t::setDim(int n)
{
  this->dim = n;
}

int mesh_t::getDim() const
{
  return this->dim;
}

void mesh_t::setNodes(int n)
{
  this->nodes = n;
}

int mesh_t::getNodes() const
{
  return this->nodes;
}

void mesh_t::setPoints(int n)
{
  this->points = n;
}

int mesh_t::getPoints() const
{
  return this->points;
}

void mesh_t::setEdges(int n)
{
  this->edges = n;
}

int mesh_t::getEdges() const
{
  return this->edges;
}

void mesh_t::setSurfaces(int n)
{
  this->surfaces = n;
}

int mesh_t::getSurfaces() const
{
  return this->surfaces;
}

void mesh_t::setElements(int n)
{
  this->elements = n;
}

int mesh_t::getElements() const
{
  return this->elements;
}

node_t* mesh_t::getNode(int n)
{
  return &this->node[n];
}

void mesh_t::setNodeArray(node_t* n)
{
  this->node = n;
}

void mesh_t::newNodeArray(int n)
{
  this->node = new node_t[n];
}

void mesh_t::deleteNodeArray()
{
  delete [] this->node;
}

point_t* mesh_t::getPoint(int n)
{
  return &this->point[n];
}

void mesh_t::setPointArray(point_t* p)
{
  this->point = p;
}

void mesh_t::newPointArray(int n)
{
  this->point = new point_t[n];
}

void mesh_t::deletePointArray()
{
  delete [] this->point;
}

edge_t* mesh_t::getEdge(int n)
{
  return &this->edge[n];
}

void mesh_t::setEdgeArray(edge_t* e)
{
  this->edge = e;
}

void mesh_t::newEdgeArray(int n)
{
  this->edge = new edge_t[n];
}

void mesh_t::deleteEdgeArray()
{
  delete [] this->edge;
}

surface_t* mesh_t::getSurface(int n)
{
  return &this->surface[n];
}

void mesh_t::setSurfaceArray(surface_t* s)
{
  this->surface = s;
}

void mesh_t::newSurfaceArray(int n)
{
  this->surface = new surface_t[n];
}

void mesh_t::deleteSurfaceArray()
{
  delete [] this->surface;
}

element_t* mesh_t::getElement(int n)
{
  return &this->element[n];
}

void mesh_t::setElementArray(element_t* e)
{
  this->element = e;
}

void mesh_t::newElementArray(int n)
{
  this->element = new element_t[n];
}

void mesh_t::deleteElementArray()
{
  delete [] this->element;
}
