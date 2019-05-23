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

#include <QtCore>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include "meshutils.h"
using namespace std;

template <typename T> 
T **AllocateDynamicArray(int nRows, int nCols) {
      T **dynamicArray;
      dynamicArray = new T*[nRows];
      for( int i = 0 ; i < nRows ; i++ )
      dynamicArray[i] = new T [nCols];
      return dynamicArray;
}

template <typename T>
void FreeDynamicArray(T** dArray) {
      delete [] *dArray;
      delete [] dArray;
}

template <typename T> 
T *AllocateDynamicVector(int nRows) {
      T *dynamicArray;
      dynamicArray = new T[nRows];
      return dynamicArray;
}

template <typename T>
void FreeDynamicVector(T* dArray) {
      delete [] dArray;
}

Meshutils::Meshutils()
{
}


Meshutils::~Meshutils()
{
}


// Clear mesh...
//-----------------------------------------------------------------------------
void Meshutils::clearMesh(mesh_t *mesh)
{
  if(mesh == NULL)
    return;
  
  mesh->clear();
  delete mesh;

  cout << "Old mesh cleared" << endl;
  cout.flush();
}



// Find surface elements for 3D elements when they are not provided 
//----------------------------------------------------------------------------
void Meshutils::findSurfaceElements(mesh_t *mesh)
{
#define RESETENTRY1             \
  h->node[0] = UNKNOWN;		\
  h->node[1] = UNKNOWN;  	\
  h->element[0] = UNKNOWN;	\
  h->element[1] = UNKNOWN;	\
  h->index[0] = UNKNOWN;        \
  h->index[1] = UNKNOWN;        \
  h->next = NULL;               

  class hashEntry {
  public:
    int node[2];
    int element[2];
    int index[2];
    int face;
    hashEntry *next;
  };


  if(mesh->getElements() == 0) return;

  int keys = mesh->getNodes();  
  hashEntry *hash = new hashEntry[keys];

  bool found;
  hashEntry *h;

  for(int i=0; i<keys; i++) {
    h = &hash[i];
    RESETENTRY1;
  }

  // TODO: only 1st and 2nd order elements

  static int familyfaces[9] = {0, 0, 0, 0, 0, 4, 5, 5, 6};

  static int faceedges8[] = {4, 4, 4, 4, 4, 4};
  static int facemap8[][8] = {{0,1,2,3,8,9,10,11}, {4,5,6,7,16,17,18,19}, {0,1,5,4,8,13,16,12}, {1,2,6,5,9,14,17,13}, {2,3,7,6,10,15,18,14}, {3,0,4,7,11,12,19,15}};

  static int faceedges7[] = {3, 3, 4, 4, 4};
  static int facemap7[][8] = {{0,1,2,6,7,8}, {3,4,5,12,13,14}, {0,1,4,3,6,10,12,9}, {1,2,5,4,7,11,13,10}, {2,0,3,5,8,9,14,11}};

  static int faceedges6[] = {4, 3, 3, 3, 3};
  static int facemap6[][8] = {{0,1,2,3,5,6,7,8}, {0,1,4,5,10,9}, {1,2,4,6,11,10}, {2,3,4,7,12,11}, {3,0,4,8,9,12}};

  static int faceedges5[4] = {3, 3, 3, 3};
  static int facemap5[][6] = {{0,1,2,4,5,6}, {0,1,3,4,8,7}, {1,2,3,5,9,8}, {2,0,3,6,7,9}};


  for(int i=0; i < mesh->getElements(); i++) {
    element_t *e = mesh->getElement(i);

    int facenodes = 0;
    int *facemap = NULL;
    int n[4];

    int family = e->getCode() / 100;
    int faces = familyfaces[family];

    for(int f=0; f<faces; f++) {
      if(family == 5) {
	facenodes = faceedges5[f];
	facemap = &facemap5[f][0];
      }
      else if(family == 6) {
	facenodes = faceedges6[f];
	facemap = &facemap6[f][0];
      }
      else if(family == 7) {
	facenodes = faceedges7[f];
	facemap = &facemap7[f][0];
      }
      else if(family == 8) {
	facenodes = faceedges8[f];
	facemap = &facemap8[f][0];
      }

      for(int j=0; j < facenodes; j++) 
	n[j] = e->getNodeIndex(facemap[j]);

      // Order indexes in an increasing order 
      for(int k=facenodes-1;k>0;k--) {
	for(int j=0;j<k;j++) {
	  if(n[j] > n[j+1]) {
	    int tmp = n[j+1];
	    n[j+1] = n[j];
	    n[j] = tmp;
	  } 
	}
      }
      
      // three nodes define a face uniquely also for rectangles
      h = &hash[n[0]];
      found = false;
      while(h->next) {                                       
	if((h->node[0] == n[1]) && (h->node[1] == n[2])) {
	  found = true;
	  break;
	}
	h = h->next;
      }                                                      
      
      if(!found) {
	h->node[0] = n[1];
	h->node[1] = n[2];
	h->element[0] = i;
	h->index[0] = e->getIndex();
	h->face = f;
	
	h->next = new hashEntry;
	h = h->next;
	RESETENTRY1;
      } else {
	h->index[1] = e->getIndex();
	h->element[1] = i;
      }
    }
  }

  // count faces that have different materials at either sides:
  int allsurfaces = 0;
  int surfaces = 0;
  int maxindex1 = 0;
  int maxindex2 = 0;
  for(int i=0; i<keys; i++) {
    h = &hash[i];
    while(h->next){
      if(h->element[0] > UNKNOWN) allsurfaces++;
      if(h->index[0] != h->index[1]) surfaces++;
      if(h->index[0] > maxindex1) maxindex1 = h->index[0];
      if(h->index[1] > maxindex2) maxindex2 = h->index[1];
      h = h->next;
    }
  }
  cout << "Found " << surfaces << " interface faces of " << allsurfaces << endl;


  // Create a index table such that all combinations of materials 
  // get different BC index  
  //int indextable[maxindex1+1][maxindex2+1];
  int **indextable = AllocateDynamicArray<int>(maxindex1+1, maxindex2+1);
  int index1,index2;
  for(int i=0;i<=maxindex1;i++)
    for(int j=0;j<=maxindex2;j++)
      indextable[i][j] = 0;

  for(int i=0; i<keys; i++) {
    h = &hash[i];
    while(h->next){
      if(h->index[0] != h->index[1]) {
	index1 = h->index[0];
	index2 = h->index[1];
	if(index2 == -1) index2 = 0;
	indextable[index1][index2] = 1;
      }
      h = h->next;
    }
  }
  index1=0;
  for(int i=0;i<=maxindex1;i++)
    for(int j=0;j<=maxindex2;j++)
      if(indextable[i][j]) indextable[i][j] = ++index1;

  cout << "Boundaries were numbered up to index " << index1 << endl;

  // Finally set the surfaces:
  mesh->setSurfaces(surfaces);
  mesh->newSurfaceArray(mesh->getSurfaces());

  surfaces = 0;
  for(int i=0; i<keys; i++) {
    h = &hash[i];

    while(h->next) {
      if(h->index[0] != h->index[1]) {
	surface_t *s = mesh->getSurface(surfaces);

	s->setElements(1);
	s->newElementIndexes(2); 
	s->setElementIndex(0, h->element[0]);
	s->setElementIndex(1, h->element[1]);
	if(s->getElementIndex(1) >= 0) s->setElements(2); // ?????
	
	element_t *e = mesh->getElement(h->element[0]);
	int code = e->getCode();
	int family = code / 100;
	int f = h->face;

	int faceedges = 0;
	int *facemap = NULL;
	int degree = 1;

	if(family == 5) {
	  faceedges = faceedges5[f];
	  if( code == 510) degree = 2;
	  facemap = &facemap5[f][0];
	}
	else if(family == 6) {
	  faceedges = faceedges6[f];
	  if( code == 613) degree = 2;
	  facemap = &facemap6[f][0];
	}
	else if(family == 7) {
	  faceedges = faceedges7[f];
	  if( code == 715) degree = 2;
	  facemap = &facemap7[f][0];
	}
	else if(family == 8) {
	  faceedges = faceedges8[f];
	  if( code == 820 ) degree = 2;
	  facemap = &facemap8[f][0];
	}
	
	int facenodes = degree * faceedges;

	s->setNodes(facenodes);
	s->setCode(100 * faceedges + facenodes);
	s->newNodeIndexes(s->getNodes());
	for(int j=0; j < s->getNodes(); j++) 
	  s->setNodeIndex(j, e->getNodeIndex(facemap[j]));

	index1 = h->index[0];
	index2 = h->index[1];
	if(index2 < 0) index2 = 0;
	
	s->setIndex(indextable[index1][index2]);

	s->setEdges(s->getNodes());
	s->newEdgeIndexes(s->getEdges());
	for(int j=0; j < s->getEdges(); j++)
	  s->setEdgeIndex(j, UNKNOWN);

	s->setNature(PDE_BOUNDARY);
	surfaces++;
      }
      h = h->next;
    }
  }

  delete [] hash;
}



// Find parent elements for existing surfaces...
//----------------------------------------------------------------------------
void Meshutils::findEdgeElementParents(mesh_t *mesh)
{
#define RESETENTRY0             \
  h->node[0] = UNKNOWN;		\
  h->node[1] = UNKNOWN;  	\
  h->element[0] = UNKNOWN;	\
  h->element[1] = UNKNOWN;	\
  h->next = NULL;

  class hashEntry {
  public:
    int node[2];
    int element[2];
    hashEntry *next;
  };

  int keys = mesh->getNodes();
  hashEntry *hash = new hashEntry[keys];

  bool found;
  hashEntry *h;

  for(int i=0; i<keys; i++) {
    h = &hash[i];
    RESETENTRY0;
  }

  // TODO: only tetrahedron at the moment

  static int edgemap[][2] = {{0,1}, {1,2}, {2,0}};
  
  for(int i = 0; i < mesh->getSurfaces(); i++) {
    surface_t *s = mesh->getSurface(i);

    for(int f = 0; f < 3; f++) {
      int n0 = s->getNodeIndex(edgemap[f][0]);
      int n1 = s->getNodeIndex(edgemap[f][1]);

      if(n1 < n0) {
	int tmp = n1;
	n1 = n0;
	n0 = tmp;
      }
      
      h = &hash[n0];
      found = false;
      while(h->next) {                                       
	if(h->node[0] == n1) {
	  found = true;
	  break;
	}
	h = h->next;
      }                                                      
      
      if(!found) {
	h->node[0] = n1;
	h->element[0] = i;
	h->next = new hashEntry;
	h = h->next;
	RESETENTRY0;
      } else {
	h->element[1] = i;
      }      
    }
  }

  // count faces:
  int edges = 0;
  for(int i = 0; i < keys; i++) {
    h = &hash[i];
    while((h = h->next) != NULL) 
      edges++;
  }
  
  cout << "Found total of " << edges << " edges" << endl;

  // Finally find parents:
  for(int i = 0; i < mesh->getEdges(); i++) {
    edge_t *e = mesh->getEdge(i);
    
    int n0 = e->getNodeIndex(0);
    int n1 = e->getNodeIndex(1);
    
    if(n1 < n0) {
      int tmp = n1;
      n1 = n0;
      n0 = tmp;
    }
    
    h = &hash[n0];
    while(h->next) {
      if(h->node[0] == n1) {

	// should we deallocate s->element if it exists?
	e->setSurfaces(2);
	e->newSurfaceIndexes(2);

	e->setSurfaceIndex(0, h->element[0]);
	e->setSurfaceIndex(1, h->element[1]);
      }
      h = h->next;
    }
  }

  delete [] hash;
}


// Find parent elements for existing surfaces...
//----------------------------------------------------------------------------
void Meshutils::findSurfaceElementParents(mesh_t *mesh)
{
#define RESETENTRY0             \
  h->node[0] = UNKNOWN;		\
  h->node[1] = UNKNOWN;  	\
  h->element[0] = UNKNOWN;	\
  h->element[1] = UNKNOWN;	\
  h->next = NULL;

  class hashEntry {
  public:
    int node[2];
    int element[2];
    hashEntry *next;
  };

  int keys = mesh->getNodes();
  hashEntry *hash = new hashEntry[keys];

  bool found;
  hashEntry *h;

  for(int i=0; i<keys; i++) {
    h = &hash[i];
    RESETENTRY0;
  }

  // TODO: only tetrahedron at the moment

  static int facemap[][3] = {{0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}};
  
  for(int i=0; i < mesh->getElements(); i++) {
    element_t *e = mesh->getElement(i);

    for(int f=0; f<4; f++) {
      int n0 = e->getNodeIndex(facemap[f][0]);
      int n1 = e->getNodeIndex(facemap[f][1]);
      int n2 = e->getNodeIndex(facemap[f][2]);

      if(n2 < n1) {
	int tmp = n2;
	n2 = n1;
	n1 = tmp;
      }
      
      if(n2 < n0) {
	int tmp = n2;
	n2 = n0;
	n0 = tmp;
      }
      
      if(n1 < n0) {
	int tmp = n1;
	n1 = n0;
	n0 = tmp;
      }
      
      h = &hash[n0];
      found = false;
      while(h->next) {                                       
	if((h->node[0] == n1) && (h->node[1] == n2)) {
	  found = true;
	  break;
	}
	h = h->next;
      }                                                      
      
      if(!found) {
	h->node[0] = n1;
	h->node[1] = n2;
	h->element[0] = i;
	h->next = new hashEntry;
	h = h->next;
	RESETENTRY0;
      } else {
	h->element[1] = i;
      }      
    }
  }

  // count faces:
  int faces = 0;
  for(int i=0; i<keys; i++) {
    h = &hash[i];
    while((h = h->next) != NULL) 
      faces++;
  }
  
  cout << "Found total of " << faces << " faces" << endl;

  // Finally find parents:
  for(int i=0; i < mesh->getSurfaces(); i++) {
    surface_t *s = mesh->getSurface(i);
    
    int n0 = s->getNodeIndex(0);
    int n1 = s->getNodeIndex(1);
    int n2 = s->getNodeIndex(2);
    
    if(n2 < n1) {
      int tmp = n2;
      n2 = n1;
      n1 = tmp;
    }
    
    if(n2 < n0) {
      int tmp = n2;
      n2 = n0;
      n0 = tmp;
    }
    
    if(n1 < n0) {
      int tmp = n1;
      n1 = n0;
      n0 = tmp;
    }
    
    h = &hash[n0];
    while(h->next) {
      if((h->node[0] == n1) && (h->node[1] == n2)) {

	// should we deallocate s->element if it exists?
	s->setElements(2);
	s->newElementIndexes(2);

	s->setElementIndex(0, h->element[0]);
	s->setElementIndex(1, h->element[1]);
      }
      h = h->next;
    }
  }

  delete [] hash;
}



// Find points for edge elements...
//-----------------------------------------------------------------------------
void Meshutils::findEdgeElementPoints(mesh_t *mesh)
{
  class hashEntry {
  public:
    int edges;
    int *edge;
  };

  int keys = mesh->getNodes();
  
  hashEntry *hash = new hashEntry[keys];

  for(int i = 0; i < keys; i++) {
    hashEntry *h = &hash[i];
    h->edges = 0;
    h->edge = NULL;
  }

  for(int i = 0; i < mesh->getEdges(); i++) {
    edge_t *e = mesh->getEdge(i);

    if(e->getNature() == PDE_BOUNDARY) {      
      for(int k = 0; k < 2; k++) {
	int n = e->getNodeIndex(k);
	
	hashEntry *h = &hash[n];
	
	bool found = false;
	for(int j = 0; j < h->edges; j++) {
	  if(h->edge[j] == i) {
	    found = true;
	    break;
	  }
	}
	
	if(!found) {
	  int *tmp = new int[h->edges+1];
	  for(int j = 0; j < h->edges; j++)
	    tmp[j] = h->edge[j];
	  tmp[h->edges] = i;
	  delete [] h->edge;
	  
	  h->edges++;
	  h->edge = tmp;	
	}
      }
    }
  }
  
  // count points:
  int count = 0;
  for(int i = 0; i < keys; i++) {
    hashEntry *h = &hash[i];
    if(h->edges > 0) 
      count++;
  }

  cout << "Found " << count << " points on boundary edges" << endl;
  cout.flush();

  // delete old points, if any:
  if(mesh->getPoints() > 0) {
    cout << "Deleteing old points and creating new" << endl;
    cout.flush();
    for(int i = 0; i < mesh->getPoints(); i++) {
      mesh->getPoint(i)->deleteNodeIndexes();
      mesh->getPoint(i)->deleteEdgeIndexes();
    }
    mesh->deletePointArray();
  }

  mesh->setPoints(count);
  mesh->newPointArray(mesh->getPoints());

  count = 0;
  for(int i = 0; i < keys; i++) {
    hashEntry *h = &hash[i];
    
    if(h->edges > 0) {
      point_t *p = mesh->getPoint(count++);
      p->setNodes(1);
      p->newNodeIndexes(1);
      p->setNodeIndex(0, i);
      p->setEdges(h->edges);
      p->newEdgeIndexes(p->getEdges());
      for(int j = 0; j < p->getEdges(); j++) {
	p->setEdgeIndex(j, h->edge[j]);
      }
      p->setSharp(false);
    }
  }

  // delete temp stuff
  for(int i = 0; i < keys; i++) {
    hashEntry *h = &hash[i];
    if(h->edges > 0)
      delete [] h->edge;
  }

  delete [] hash;

  // Inverse map
  cout << "Constructing inverse map from edges to points" << endl;
  cout.flush();

  for(int i=0; i < mesh->getPoints(); i++) {
    point_t *p = mesh->getPoint(i);

    for(int j=0; j < p->getEdges(); j++) {
      int k = p->getEdgeIndex(j);
      if ( k<0 ) continue;

      edge_t *e = mesh->getEdge(k);

      // allocate space for two points, if not yet done:
      if(e->getPoints() < 2) {
	e->setPoints(2);
	e->newPointIndexes(2);
	e->setPointIndex(0, -1);
	e->setPointIndex(1, -1);
      }
            
      for(int r=0; r < e->getPoints(); r++) {
	if(e->getPointIndex(r) < 0) {
	  e->setPointIndex(r, i); 
	  break;
	}
      }
    }
  }  
}



// Find edges for surface elements...
//-----------------------------------------------------------------------------
void Meshutils::findSurfaceElementEdges(mesh_t *mesh)
{
#define RESETENTRY              \
    h->surfaces = 0;            \
    h->surface = NULL;          \
    h->parentindex[0] = UNKNOWN; \
    h->parentindex[1] = UNKNOWN; \
    h->next = NULL;

  int keys = mesh->getNodes();

  class hashEntry {
  public:
    int nodes;
    int *node;
    int nature;
    int index;
    int surfaces;
    int *surface;
    int parentindex[2];
    hashEntry *next;
  };
  
  hashEntry *hash = new hashEntry[keys];

  bool found;
  hashEntry *h;

  for(int i=0; i<keys; i++) {
    h = &hash[i];
    RESETENTRY;
  }

  bool createindexes = ((!mesh->getEdges()) && (!mesh->getElements()));

  // ????? should we test for mesh->edges ?????
  // if ( mesh->edge && (mesh->getEdges() > 0) ) {  
  if ( mesh->getEdges() > 0 ) {  
    // add existing edges first:
    for( int i=0; i<mesh->getEdges(); i++ )
    {
      edge_t *edge = mesh->getEdge(i);
      int n0 = edge->getNodeIndex(0);
      int n1 = edge->getNodeIndex(1);

      int m = (n0<n1) ? n0 : n1;
      int n = (n0<n1) ? n1 : n0;

      h = &hash[m];
      found = false;
      while(h->next) {                                       
        if(h->node[0] == n) {
          found = true;
          break;
        }
        h = h->next;
      }                                                      
      
      if(!found) {
        h->nodes = edge->getNodes() - 1;
        h->node = new int[h->nodes];
        h->node[0] = n;
        for( int j=1; j<h->nodes; j++ )
        {
          h->node[j] = edge->getNodeIndex(j+1);
        }
        h->surfaces = edge->getSurfaces();
        h->surface = new int[edge->getSurfaces()];
        for( int j=0; j < edge->getSurfaces(); j++ ) {
          h->surface[j] = edge->getSurfaceIndex(j);
        }
        h->index  = edge->getIndex();
        h->nature = edge->getNature();
        h->next = new hashEntry;
        h = h->next;
        RESETENTRY;
      }
    }

    mesh->setEdges(0);
    mesh->deleteEdgeArray(); // delete [] mesh->edge;
  }


  static int triedgemap[][4] = { {0,1,3,6}, {1,2,4,7}, {2,0,5,8} };
  static int quadedgemap[][4] = {{0,1,4,8}, {1,2,5,9}, {2,3,6,10}, {3,0,7,11}};

  
  for(int i=0; i < mesh->getSurfaces(); i++) {
    surface_t *s = mesh->getSurface(i);
    
    // loop over edges
    for(int e=0; e < s->getEdges(); e++) {
      int n0, n1,nrest[2] = { -1,-1 };
      if((int)(s->getCode() / 100) == 3) {
	n0 = s->getNodeIndex(triedgemap[e][0]);
	n1 = s->getNodeIndex(triedgemap[e][1]);
	if ( s->getCode() >= 306) nrest[0] = s->getNodeIndex(triedgemap[e][2]);
	if ( s->getCode() >= 309) nrest[1] = s->getNodeIndex(triedgemap[e][3]);
      } else if((int)(s->getCode() / 100) == 4) {
	n0 = s->getNodeIndex(quadedgemap[e][0]);
	n1 = s->getNodeIndex(quadedgemap[e][1]);
	// corrected 4.11.2008:
	if ( s->getCode() >= 408) nrest[0] = s->getNodeIndex(quadedgemap[e][2]);
	if ( s->getCode() >= 412) nrest[1] = s->getNodeIndex(quadedgemap[e][3]);
      } else {
	cout << "findBoundaryElementEdges: error: unknown element code" << endl;
	exit(0);
      }

      int m = (n0<n1) ? n0 : n1;
      int n = (n0<n1) ? n1 : n0;

      h = &hash[m];
      found = false;
      while(h->next) {                                       
	if(h->node[0] == n) {
	  found = true;
	  break;
	}
	h = h->next;
      }                                                      
      
      if(!found) {
        h->nodes = 1;
        h->nodes += nrest[0] >=0 ? 1:0;
        h->nodes += nrest[1] >=0 ? 1:0;
        h->node = new int[h->nodes];
        h->node[0] = n;
        for( int j=1; j<h->nodes; j++ )
          h->node[j] = nrest[j-1];

	h->surfaces = 1;
	h->surface = new int[1];
	h->surface[0] = i;
	h->parentindex[0] = s->getIndex();
        h->index = UNKNOWN;
        h->nature = PDE_UNKNOWN;
	h->next = new hashEntry;
	h = h->next;
	RESETENTRY;
      } else {
        int *tmp  = new int[h->surfaces];

        found = false;
        for(int j=0; j<h->surfaces; j++)
        {
          tmp[j] = h->surface[j];
          if ( tmp[j] == i ) found = true; 
        }
        if ( found ) {
          delete [] tmp;
        } else {
          delete [] h->surface;
          h->surface = new int[h->surfaces+1];
          for(int j=0; j<h->surfaces; j++)
             h->surface[j] = tmp[j];
          h->surface[h->surfaces++] = i;
	  if( s->getIndex() < h->parentindex[0]) {	    
	    h->parentindex[1] = h->parentindex[0];
	    h->parentindex[0] = s->getIndex();
	  }
	  else {
	    h->parentindex[1] = s->getIndex();
	  }	    
          delete [] tmp;
        }
      }
    }
  }

  // count edges:
  int edges = 0;
  for(int i=0; i<keys; i++) {
    h = &hash[i];
    while((h = h->next) != NULL) 
      edges++;
  }

    

  cout << "Found " << edges << " edges on boundary" << endl;

  mesh->setEdges(edges);
  mesh->newEdgeArray(edges); // edge = new edge_t[edges];

  if(createindexes) {
    cout << "Creating edge indexes " << endl;

    int maxindex1 = UNKNOWN;
    int maxindex2 = UNKNOWN;
    
    for(int i=0; i<keys; i++) {
      h = &hash[i];
      while(h->next){
	if(h->parentindex[0] > maxindex1) maxindex1 = h->parentindex[0];
	if(h->parentindex[1] > maxindex2) maxindex2 = h->parentindex[1];
	h = h->next;
      }
    }
 
    // Create a index table such that all combinations of materials 
    // get different BC index  
    //int indextable[maxindex1+2][maxindex2+2];
    int **indextable = AllocateDynamicArray<int>(maxindex1+2, maxindex2+2);
    int index1,index2;
    for(int i=-1;i<=maxindex1;i++)
      for(int j=-1;j<=maxindex2;j++)
	indextable[i+1][j+1] = 0;

    int edgebcs=0;
    for(int i=0; i<keys; i++) {
      h = &hash[i];
      while(h->next){
	if(h->parentindex[0] != h->parentindex[1]) {
	  index1 = h->parentindex[0];
	  index2 = h->parentindex[1];
	  indextable[index1+1][index2+1] = 1;
	  edgebcs += 1;
	}
	h = h->next;
      }
    }

    index1=0;
    for(int i=-1;i<=maxindex1;i++)
      for(int j=-1;j<=maxindex2;j++) 
	if(indextable[i+1][j+1]) indextable[i+1][j+1] = ++index1;

    cout << edgebcs << " boundary edges were numbered up to index " << index1 << endl;
    
    for(int i=0; i<keys; i++) {
      h = &hash[i];
      while(h->next){
	if(h->parentindex[0] != h->parentindex[1]) {
	  index1 = h->parentindex[0];
	  index2 = h->parentindex[1];
	  h->index = indextable[index1+1][index2+1];
	  h->nature = PDE_BOUNDARY;
	}
	h = h->next;
      }
    }
  }


  // Create edges:
  edges = 0;
  for(int i=0; i<keys; i++) {
    h = &hash[i];
    while(h->next) {
      edge_t *e = mesh->getEdge(edges++);
      
      e->setNature(h->nature);
      e->setNodes(h->nodes + 1);
      e->newNodeIndexes(e->getNodes());
      e->setNodeIndex(0, i); // ?????
      for( int j=1; j < e->getNodes(); j++ )
        e->setNodeIndex(j, h->node[j-1]);

      e->setCode(200 + e->getNodes());

      e->setSurfaces(h->surfaces);
      e->newSurfaceIndexes(max(e->getSurfaces(), 2));
      e->setSurfaceIndex(0, -1);
      e->setSurfaceIndex(1, -1);

      for(int j=0; j < e->getSurfaces(); j++)
	e->setSurfaceIndex(j, h->surface[j]);

      e->setSharp(false);

      e->setIndex(h->index);
      e->setPoints(0);
      h = h->next;
    }
  }

  delete [] hash;

  // Inverse map
  for(int i=0; i < mesh->getEdges(); i++) {
    edge_t *e = mesh->getEdge(i);

    for(int j=0; j < e->getSurfaces(); j++) {
      int k = e->getSurfaceIndex(j);
      if ( k< 0 ) continue;

      surface_t *s = mesh->getSurface(k);
      
      for(int r=0; r < s->getEdges(); r++) {
	if(s->getEdgeIndex(r) < 0) {
	  s->setEdgeIndex(r, i);
	  break;
	}
      }
    }
  }  

#if 0
  cout << "*********************" << endl;
  for(int i=0; i<mesh->getEdges(); i++)
    cout << "Edge " << i << " nodes " << mesh->edge[i].node[0] << " "<< mesh->edge[i].node[0] << endl;

  for(int i=0; i < mesh->getSurfaces(); i++)
    cout << "Surface " << i << " nodes " 
	 << mesh->surface[i].node[0] << " " 
	 << mesh->surface[i].node[1] << " "
	 << mesh->surface[i].node[2] << " "
	 << " Edges " 
	 << mesh->surface[i].edge[0] << " " 
	 << mesh->surface[i].edge[1] << " "
	 << mesh->surface[i].edge[2] << " "
	 << " Parents " 
	 << mesh->surface[i].element[0] << " " 
	 << mesh->surface[i].element[1] << " "
	 << endl;

  cout.flush();
#endif

}




// Find sharp points for edge elements...
//-----------------------------------------------------------------------------
void Meshutils::findSharpPoints(mesh_t *mesh, double limit)
{
  qreal t0[3], t1[3];

  cout << "Limit: " << limit << " degrees" << endl;
  cout.flush();
  
  double angle;
  int count = 0;
  point_t *point = NULL;
  edge_t *edge = NULL;
  Helpers *helpers = new Helpers;
  
  for(int i=0; i < mesh->getPoints(); i++) {
    point = mesh->getPoint(i);

    if(point->getEdges() == 2) {
      int n = point->getNodeIndex(0);

      int e0 = point->getEdgeIndex(0);
      int e1 = point->getEdgeIndex(1);

      edge = mesh->getEdge(e0);
      int n0 = edge->getNodeIndex(0);
      if(edge->getNodeIndex(1) != n)
	n0 = edge->getNodeIndex(1);

      edge = mesh->getEdge(e1);
      int n1 = edge->getNodeIndex(0);
      if(edge->getNodeIndex(1) != n)
	n1 = edge->getNodeIndex(1);

      // unit tangent from node to node0
      t0[0] = mesh->getNode(n0)->getX(0) - mesh->getNode(n)->getX(0);
      t0[1] = mesh->getNode(n0)->getX(1) - mesh->getNode(n)->getX(1);
      t0[2] = mesh->getNode(n0)->getX(2) - mesh->getNode(n)->getX(2);
      
      // unit tangent from node to node1
      t1[0] = mesh->getNode(n1)->getX(0) - mesh->getNode(n)->getX(0);
      t1[1] = mesh->getNode(n1)->getX(1) - mesh->getNode(n)->getX(1);
      t1[2] = mesh->getNode(n1)->getX(2) - mesh->getNode(n)->getX(2);
      
      helpers->normalize(t0);
      helpers->normalize(t1);

      double cosofangle = t0[0]*t1[0] + t0[1]*t1[1] + t0[2]*t1[2];
      angle = acos(cosofangle) / MYPI * 180.0;
    } else {
      angle = 0.0;
    }    
    
    point->setSharp(false);
    if(sqrt(angle*angle) < (180.0-limit) ) {
      point->setSharp(true);
      count++;
    }
  }

  cout << "Found " << count << " sharp points" << endl;
  delete helpers;
}



// Find sharp edges for surface elements...
//-----------------------------------------------------------------------------
void Meshutils::findSharpEdges(mesh_t *mesh, double limit)
{

  cout << "Limit: " << limit << " degrees" << endl;
  cout.flush();

  double angle;
  int count = 0;
  
  for(int i=0; i<mesh->getEdges(); i++) {
    edge_t *edge = mesh->getEdge(i);

    if(edge->getSurfaces() == 2) {
      int s0 = edge->getSurfaceIndex(0);
      int s1 = edge->getSurfaceIndex(1);
      double *n0 = mesh->getSurface(s0)->getNormalVec();
      double *n1 = mesh->getSurface(s1)->getNormalVec();
      double cosofangle = n0[0]*n1[0] + n0[1]*n1[1] + n0[2]*n1[2];
	  cosofangle = std::abs(cosofangle);
      angle = acos(cosofangle) / MYPI * 180.0;
    } else {
      angle = 180.0;
    }    
    
    edge->setSharp(false);
    if(sqrt(angle*angle) > limit) {
      edge->setSharp(true);
      count++;
    }
  }

  cout << "Found " << count << " sharp edges" << endl;
}



// Divide edge by sharp points...
//-----------------------------------------------------------------------------
int Meshutils::divideEdgeBySharpPoints(mesh_t *mesh)
{
  class Bc {
  public:
    void propagateIndex(mesh_t* mesh, int index, int i) {
      edge_t *edge = mesh->getEdge(i);

      // index is ok
      if(!edge->getSelected() || (edge->getIndex() != UNKNOWN)
	 || (edge->getNature() != PDE_BOUNDARY)) return;
      
      // set index
      edge->setIndex(index);

      // propagate index
      for(int j=0; j < edge->getPoints(); j++) {
	int k = edge->getPointIndex(j);
	point_t *point = mesh->getPoint(k);

	// skip sharp points
	if(!point->isSharp()) {
	  for(int m = 0; m < point->getEdges(); m++) {
	    int n = point->getEdgeIndex(m);
	    propagateIndex(mesh, index, n);
	  }
	}

      }
    }
  };
  
  
  // reset bc-indices on edges:
  int index = 0;
  int count = 0;

  for(int i=0; i < mesh->getEdges(); i++)
  {
    edge_t *edge = mesh->getEdge(i);
    if((edge->getNature() == PDE_BOUNDARY) && !edge->getSelected())
      index = max(index, edge->getIndex());
  }
  index++;

  // int edge_index[index];
  int *edge_index = new int[index];

  for( int i=0; i<index; i++ )
    edge_index[i] = UNKNOWN;

  index = 0;
  for(int i=0; i < mesh->getEdges(); i++)
  {
    edge_t *edge = mesh->getEdge(i);
    if (edge->getNature() == PDE_BOUNDARY)
    {
      if ( edge->getSelected() ) {
        count++;
        edge->setIndex(UNKNOWN);
      } else {
        if ( edge_index[edge->getIndex()] == UNKNOWN )
          edge_index[edge->getIndex()] = ++index;
        edge->setIndex(edge_index[edge->getIndex()]);
      }
    }
  }

  if ( count==0 ) {
    cout << "No boundary edges to divde." << endl;
    return 0;
  }

  Bc *bc = new Bc;

  // recursively determine boundary parts:
  for(int i=0; i < mesh->getEdges(); i++)
  {
    edge_t *edge = mesh->getEdge(i);
    if(edge->getSelected() && (edge->getIndex() == UNKNOWN) && (edge->getNature() == PDE_BOUNDARY))
      bc->propagateIndex(mesh, ++index, i);
  }
  index++;
  
  // Create a hopefully mesh indepedent indexing of groupings to enable
  // reapplying merge/division operations after remeshing. The indices
  // are given based on group bounding box corner distances from a given
  // point. Fails if two groups have same bbox, which should not happen 
  // often though (?)

  //double xmin[index], ymin[index], zmin[index];
  //double xmax[index], ymax[index], zmax[index];
  //double xc = 0, yc = 0, zc = 0, dist[index];
  //int cc[index], order[index], sorder[index];
  //double g_xmin,g_xmax,g_ymin,g_ymax,g_zmin,g_zmax;

  double *xmin = new double[index];
  double *ymin = new double[index];
  double *zmin = new double[index];
  double *xmax = new double[index];
  double *ymax = new double[index];
  double *zmax = new double[index];
  double xc = 0, yc = 0, zc = 0;
  double *dist = new double[index];
  int *cc = new int[index];
  int *order = new int[index];
  int *sorder = new int[index];
  double g_xmin,g_xmax,g_ymin,g_ymax,g_zmin,g_zmax;

  for( int i=0; i<index; i++ )
  {
    cc[i] = 0;
    order[i] = i;
    xmin[i] = ymin[i] = zmin[i] =  1e20;
    xmax[i] = ymax[i] = zmax[i] = -1e20;
  }

  for( int i=0; i<mesh->getEdges(); i++ )
  {
    edge_t *edge = mesh->getEdge(i);
    if (edge->getNature() == PDE_BOUNDARY) {
      int k = edge->getIndex();
      for( int j=0; j < edge->getNodes(); j++ ) {
        int n = edge->getNodeIndex(j);
        cc[k]++;
        xmin[k] = min(xmin[k], mesh->getNode(n)->getX(0));
        ymin[k] = min(ymin[k], mesh->getNode(n)->getX(1));
        zmin[k] = min(zmin[k], mesh->getNode(n)->getX(2));
 
        xmax[k] = max(xmax[k], mesh->getNode(n)->getX(0));
        ymax[k] = max(ymax[k], mesh->getNode(n)->getX(1));
        zmax[k] = max(zmax[k], mesh->getNode(n)->getX(2));
       }
    }
  }

  g_xmin = g_ymin = g_zmin =  1e20;
  g_xmax = g_ymax = g_zmax = -1e20;
  for( int i=0; i<index; i++)
  {
    g_xmin = min(xmin[i],g_xmin);
    g_ymin = min(ymin[i],g_ymin);
    g_zmin = min(zmin[i],g_zmin);

    g_xmax = max(xmax[i],g_xmax);
    g_ymax = max(ymax[i],g_ymax);
    g_zmax = max(zmax[i],g_zmax);
  }

  double g_scale = max(max(g_xmax-g_xmin,g_ymax-g_ymin),g_zmax-g_zmin);
  double g_xp = g_xmax + 32.1345  * g_scale;
  double g_yp = g_ymin - 5.3*MYPI * g_scale;
  double g_zp = g_zmax + 8.1234   * g_scale;

  for( int i=0; i<index; i++ )
  {
    dist[i] = 0;
    if ( cc[i]>0 ) {
      for( int j=0; j<8; j++ ) {
        switch(j) {
          case 0: xc=xmin[i]; yc=ymin[i]; zc=zmin[i]; break;
          case 1: xc=xmax[i]; break;
          case 2: yc=xmax[i]; break;
          case 3: xc=xmin[i]; break;
          case 4: zc=zmax[i]; break;
          case 5: yc=ymin[i]; break;
          case 6: xc=xmax[i]; break;
          case 7: yc=ymax[i]; break;
        }
        dist[i] += (xc-g_xp)*(xc-g_xp);
        dist[i] += (yc-g_yp)*(yc-g_yp);
        dist[i] += (zc-g_zp)*(zc-g_zp);
      }
    }
  }

  sort_index( index, dist, order );
  for( int i=0; i<index; i++ )
    sorder[order[i]] = i;

  for( int i=0; i<mesh->getEdges(); i++ )
    if ( mesh->getEdge(i)->getNature() == PDE_BOUNDARY )
      mesh->getEdge(i)->setIndex(sorder[mesh->getEdge(i)->getIndex()]);

  --index;
  cout << "Edge divided into " << index << " parts" << endl;
  
  delete bc;

  return index;
}


void Meshutils::sort_index(int n, double *a, int *b)
{
  QList<QPair<double, int> > list;

  for(int i = 0; i < n; i++)
    list << qMakePair(a[i], b[i]);

  qSort(list);

  for(int i = 0; i < n; i++) {
    a[i] = list[i].first;
    b[i] = list[i].second;
  }
}


// Divide surface by sharp edges...
//-----------------------------------------------------------------------------
int Meshutils::divideSurfaceBySharpEdges(mesh_t *mesh)
{
  class Bc {
  public:
    void propagateIndex(mesh_t* mesh, int index, int i) {
      surface_t *surf = mesh->getSurface(i);

      // index is ok
      if(!surf->getSelected() || (surf->getIndex() != UNKNOWN) 
	 || (surf->getNature() != PDE_BOUNDARY) ) return;

      // set index
      surf->setIndex(index);

      // propagate index
      for(int j=0; j < surf->getEdges(); j++) {
	int k = surf->getEdgeIndex(j);
	edge_t *edge = mesh->getEdge(k);

	// skip sharp edges
	if(!edge->isSharp()) {
	  for(int m=0; m < edge->getSurfaces(); m++) {
	    int n = edge->getSurfaceIndex(m);
	    propagateIndex(mesh, index, n);
	  }
	}
      }
    }
  };
  
  // reset bc-indices:
  int count = 0;
  int index = 0;

  for( int i=0; i<mesh->getSurfaces(); i++ )
  {
    surface_t *surf = mesh->getSurface(i);
    if( (surf->getNature() == PDE_BOUNDARY) && !surf->getSelected() )
      index = max(index, surf->getIndex());
  }
  index++;

  // int surf_index[index];
  int *surf_index = new int[index];

  for( int i=0; i<index; i++ )
    surf_index[i] = UNKNOWN;

  index = 0;
  for(int i=0; i < mesh->getSurfaces(); i++)
  {
    surface_t *surf = mesh->getSurface(i);
    if (surf->getNature() == PDE_BOUNDARY) {
      if ( surf->getSelected() ) {
        count++;
        surf->setIndex(UNKNOWN);
      } else {
        if ( surf_index[surf->getIndex()] == UNKNOWN )
          surf_index[surf->getIndex()] = ++index;
        surf->setIndex(surf_index[surf->getIndex()]);
      }
    }
  }

  if ( count==0 ) {
    cout << "No boundary surfaces to divde." << endl;
    return 0;
  }


  // recursively determine boundary parts:
  Bc *bc = new Bc;

  for(int i=0; i < mesh->getSurfaces(); i++)
  {
    surface_t *surf = mesh->getSurface(i);
    if(surf->getSelected() && (surf->getIndex() == UNKNOWN) && (surf->getNature() == PDE_BOUNDARY))
      bc->propagateIndex(mesh, ++index, i);
  }
  index++;

  // Create a hopefully mesh indepedent indexing of groupings to enable
  // reapplying merge/division operations after remeshing. The indices
  // are given based on group bounding box corner distances from a given
  // point. Fails if two groups have same bbox, which should not happen 
  // often though (?)

  //double xmin[index], ymin[index], zmin[index];
  //double xmax[index], ymax[index], zmax[index];
  //double xc = 0, yc = 0, zc = 0, dist[index];
  //int cc[index], order[index], sorder[index];
  //double g_xmin,g_xmax,g_ymin,g_ymax,g_zmin,g_zmax;

  double *xmin = new double[index];
  double *ymin = new double[index];
  double *zmin = new double[index];
  double *xmax = new double[index];
  double *ymax = new double[index];
  double *zmax = new double[index];
  double xc = 0, yc = 0, zc = 0;
  double *dist = new double[index];
  int *cc = new int[index];
  int *order = new int[index];
  int *sorder = new int[index];
  double g_xmin,g_xmax,g_ymin,g_ymax,g_zmin,g_zmax;

  for( int i=0; i<index; i++ )
  {
    cc[i] = 0;
    order[i] = i;
    xmin[i] = ymin[i] = zmin[i] =  1e20;
    xmax[i] = ymax[i] = zmax[i] = -1e20;
  }

  for( int i=0; i < mesh->getSurfaces(); i++ )
  {
    surface_t *surf = mesh->getSurface(i);
    if ( mesh->getSurface(i)->getNature() == PDE_BOUNDARY ) {
      int k = surf->getIndex();
      for( int j=0; j < surf->getNodes(); j++ ) {
        int n = surf->getNodeIndex(j);
        cc[k]++;
        xmin[k] = min( xmin[k], mesh->getNode(n)->getX(0) );
        ymin[k] = min( ymin[k], mesh->getNode(n)->getX(1) );
        zmin[k] = min( zmin[k], mesh->getNode(n)->getX(2) );
 
        xmax[k] = max( xmax[k], mesh->getNode(n)->getX(0) );
        ymax[k] = max( ymax[k], mesh->getNode(n)->getX(1) );
        zmax[k] = max( zmax[k], mesh->getNode(n)->getX(2) );
      }
    }
  }

  g_xmin = g_ymin = g_zmin =  1e20;
  g_xmax = g_ymax = g_zmax = -1e20;
  for( int i=0; i<index; i++)
  {
    g_xmin = min(xmin[i],g_xmin);
    g_ymin = min(ymin[i],g_ymin);
    g_zmin = min(zmin[i],g_zmin);

    g_xmax = max(xmax[i],g_xmax);
    g_ymax = max(ymax[i],g_ymax);
    g_zmax = max(zmax[i],g_zmax);
  }

  double g_scale = max(max(g_xmax-g_xmin,g_ymax-g_ymin),g_zmax-g_zmin);
  double g_xp = g_xmax + 32.1345  * g_scale;
  double g_yp = g_ymin - 5.3*MYPI * g_scale;
  double g_zp = g_zmax + 8.1234   * g_scale;

  for( int i=0; i<index; i++ )
  {
    dist[i] = 0;
    if ( cc[i]>0 ) {
      for( int j=0; j<8; j++ ) {
        switch(j) {
          case 0: xc=xmin[i]; yc=ymin[i]; zc=zmin[i]; break;
          case 1: xc=xmax[i]; break;
          case 2: yc=xmax[i]; break;
          case 3: xc=xmin[i]; break;
          case 4: zc=zmax[i]; break;
          case 5: yc=ymin[i]; break;
          case 6: xc=xmax[i]; break;
          case 7: yc=ymax[i]; break;
        }
        dist[i] += (xc-g_xp)*(xc-g_xp);
        dist[i] += (yc-g_yp)*(yc-g_yp);
        dist[i] += (zc-g_zp)*(zc-g_zp);
      }
    }
  }

  sort_index( index, dist, order );
  for( int i=0; i<index; i++ )
    sorder[order[i]] = i;

  for( int i=0; i < mesh->getSurfaces(); i++ )
    if ( mesh->getSurface(i)->getNature() == PDE_BOUNDARY )
      mesh->getSurface(i)->setIndex(sorder[mesh->getSurface(i)->getIndex()]);
  --index;

  cout << "Surface divided into " << index << " parts" << endl;

  delete bc;

  return index;
}



// Find surface element normals...
//-----------------------------------------------------------------------------
void Meshutils::findSurfaceElementNormals(mesh_t *mesh)
{
  static qreal a[3], b[3], c[3];
  double center_surface[3], center_element[3], center_difference[3];
  Helpers *helpers = new Helpers;
  int u, v, w, e0, e1, i0, i1, bigger;

  for(int i=0; i < mesh->getSurfaces(); i++) {
    surface_t *surface = mesh->getSurface(i);
    
    u = surface->getNodeIndex(0);
    v = surface->getNodeIndex(1);

    if((int)(surface->getCode() / 100) == 3) {
      w = surface->getNodeIndex(2);
    } else if((int)(surface->getCode() / 100) == 4) {
      w = surface->getNodeIndex(3);
    } else {
      cout << "findBoundaryElementNormals: error: unknown code" << endl;
      cout.flush();
      exit(0);
    }

    // Calculate normal (modulo sign):
    a[0] = mesh->getNode(v)->getX(0) - mesh->getNode(u)->getX(0);
    a[1] = mesh->getNode(v)->getX(1) - mesh->getNode(u)->getX(1);
    a[2] = mesh->getNode(v)->getX(2) - mesh->getNode(u)->getX(2);
    
    b[0] = mesh->getNode(w)->getX(0) - mesh->getNode(u)->getX(0);
    b[1] = mesh->getNode(w)->getX(1) - mesh->getNode(u)->getX(1);
    b[2] = mesh->getNode(w)->getX(2) - mesh->getNode(u)->getX(2);
    
    helpers->crossProduct(a,b,c);
    helpers->normalize(c);
    
    surface->setNormal(0, -c[0]);
    surface->setNormal(1, -c[1]);
    surface->setNormal(2, -c[2]);
    
    // Determine sign:
    //----------------

    // a) which parent element has bigger index?

    e0 = surface->getElementIndex(0);
    e1 = surface->getElementIndex(1);
    
    if( (e0<0) && (e1<0) ) {
      // both parents unknown
      bigger = -1;
    } else if(e1<0) {
      // e0 known, e1 unknown
      bigger = e0;
    } else {
      // both parents known
      bigger = e0;
      i0 = mesh->getElement(e0)->getIndex();
      i1 = mesh->getElement(e1)->getIndex();
      if(i1 > i0)
	bigger = e1;
    }
    
    // b) normal should point to the parent with smaller index:

    if(bigger > -1) {

      // Compute center point of the surface element:
      center_surface[0] = 0.0;
      center_surface[1] = 0.0;
      center_surface[2] = 0.0;

      for(int i=0; i < surface->getNodes(); i++) {
	int j = surface->getNodeIndex(i);
	node_t *n = mesh->getNode(j);
	center_surface[0] += n->getX(0);
	center_surface[1] += n->getX(1);
	center_surface[2] += n->getX(2);
      }

      center_surface[0] /= (double)(surface->getNodes());
      center_surface[1] /= (double)(surface->getNodes());
      center_surface[2] /= (double)(surface->getNodes());
      
      element_t *e = mesh->getElement(bigger);

      // compute center point of the parent element:
      center_element[0] = 0.0;
      center_element[1] = 0.0;
      center_element[2] = 0.0;

      for(int i=0; i < e->getNodes(); i++) {
	int j = e->getNodeIndex(i);
	node_t *n = mesh->getNode(j);
	center_element[0] += n->getX(0);
	center_element[1] += n->getX(1);
	center_element[2] += n->getX(2);
      }

      center_element[0] /= (double)(e->getNodes());
      center_element[1] /= (double)(e->getNodes());
      center_element[2] /= (double)(e->getNodes());

      // difference of the centers:
      center_difference[0] = center_element[0] - center_surface[0];
      center_difference[1] = center_element[1] - center_surface[1];
      center_difference[2] = center_element[2] - center_surface[2];
      
      // dot product must be negative
      double dp = center_difference[0]*c[0]
                + center_difference[1]*c[1] 
                + center_difference[2]*c[2];
      
      if(dp > 0.0) {
	surface->setNormal(0, -surface->getNormal(0));
	surface->setNormal(1, -surface->getNormal(1));
	surface->setNormal(2, -surface->getNormal(2));

	// change orientation of the surface element:
	if(surface->getCode() == 303) {
	  int tmp = surface->getNodeIndex(1);
	  surface->setNodeIndex(1, surface->getNodeIndex(2));
	  surface->setNodeIndex(2, tmp);

	} else if(surface->getCode() == 404) {
	  int tmp = surface->getNodeIndex(1);
	  surface->setNodeIndex(1, surface->getNodeIndex(3));
	  surface->setNodeIndex(3, tmp);

	} else if(surface->getCode() == 306) {
	  int tmp = surface->getNodeIndex(1);
	  surface->setNodeIndex(1, surface->getNodeIndex(2));
	  surface->setNodeIndex(2, tmp);
	  tmp = surface->getNodeIndex(3);
	  surface->setNodeIndex(3, surface->getNodeIndex(5));
	  surface->setNodeIndex(5, tmp);

	} else if(surface->getCode() == 408) {
	  int tmp = surface->getNodeIndex(1);
	  surface->setNodeIndex(1, surface->getNodeIndex(3));
	  surface->setNodeIndex(3, tmp);
	  tmp = surface->getNodeIndex(4);
	  surface->setNodeIndex(4, surface->getNodeIndex(7));
	  surface->setNodeIndex(7, tmp);
	  tmp = surface->getNodeIndex(5);
	  surface->setNodeIndex(5, surface->getNodeIndex(6));
	  surface->setNodeIndex(6, tmp);

	} else {
	  cout << "findSurfaceElementNormals: error: unable to change element orientation" << endl;
	  cout.flush();
	  exit(0);
	}
      }
    }
  }

  for( int i=0; i < mesh->getSurfaces(); i++ )
  {
    surface_t *surface = mesh->getSurface(i);
    int n = surface->getCode() / 100;
    for(int j=0; j<n; j++ )
    {
      surface->setVertexNormalVec(j, surface->getNormalVec());
    }
  }


  // List surfaces connected to nodes
  //class n_s_t {
  //public:
  //int index;
  //n_s_t *next;
  //} n_s[mesh->getNodes()];

  class n_s_t {
  public:
    int index;
    n_s_t *next;
  };

  n_s_t *n_s = new n_s_t[mesh->getNodes()];

  for( int i=0; i<mesh->getNodes(); i++ )
  {
    n_s[i].index = -1;
     n_s[i].next = NULL;
  }

  for( int i=0; i < mesh->getSurfaces(); i++ ) 
  {
    surface_t *surface = mesh->getSurface(i);
    int n = surface->getCode() / 100;

    for( int j=0; j<n; j++ )
      {
	n_s_t *p = &n_s[surface->getNodeIndex(j)];
        if ( p->index >= 0 ) {
          n_s_t *q = new n_s_t;
          q->next = p->next;
          p->next = q;
          q->index = i;
        } else p->index = i;
      }
  }

  // average normals over surfaces connected to vertices if
  // normals within the limit_angle:
  double limit_angle = cos(50.*3.14159/180.);

  for( int i=0; i < mesh->getSurfaces(); i++ )
  {
    surface_t *surf1 = mesh->getSurface(i);
    int n = surf1->getCode() / 100;

    for( int j=0; j<n; j++ )
    {
      n_s_t *p = &n_s[surf1->getNodeIndex(j)];
      for( ; p && p->index>=0; p=p->next )
      {
        if ( p->index == i ) continue;

        surface_t *surf2 = mesh->getSurface(p->index);
        double s = 0.;

        s += surf1->getNormal(0) * surf2->getNormal(0);
        s += surf1->getNormal(1) * surf2->getNormal(1);
        s += surf1->getNormal(2) * surf2->getNormal(2);
        if ( fabs(s) > limit_angle )
        {
           if ( s > 0 ) {
	     surf1->addVertexNormalVec(j, surf2->getNormalVec());
           } else {
	     surf1->subVertexNormalVec(j, surf2->getNormalVec());
           }
        }
      }
    }
  }

  // delete lists:
  for( int i=0; i<mesh->getNodes(); i++ )
  {
     n_s_t *p=&n_s[i], *q;
     p = p->next;
     while( p )
     {
        q = p->next;
        delete p;
        p = q;
     }
  }

  // And finally normalize:
  for( int i=0; i < mesh->getSurfaces(); i++ )
  {
    surface_t *surface = mesh->getSurface(i);
    int n = surface->getCode() / 100;

    for(int j=0; j<n; j++ )
    {
       double* v = surface->getVertexNormalVec(j);
       double s = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
       if ( s != 0 ) {
         s = sqrt(s);
	 v[0] /= s;
	 v[1] /= s;
	 v[2] /= s;
       }
    }
  }


  delete helpers;
}



// Increase elementorder from 1 to 2
//----------------------------------------------------------------------------
void Meshutils::increaseElementOrder(mesh_t *mesh)
{
  class hashEntry {
  public:
    int node1;
    int noedge;
    hashEntry *next;
  };

  if((mesh->getElements() == 0) && (mesh->getSurfaces() == 0)) return;

  int keys = mesh->getNodes();  
  hashEntry *hash = new hashEntry[keys];
  hashEntry *h;

  for(int i=0; i<keys; i++) {
    h = &hash[i];
    h->node1 = UNKNOWN;
    h->noedge = UNKNOWN;
    h->next = NULL;
  }

  static int familylinnodes[9] = {0, 1, 2, 3, 4, 4, 5, 6, 8};
  static int familyedges[9]    = {0, 0, 1, 3, 4, 6, 8, 9, 12};

  static int edgemap8[][2] = {{0,1}, {1,2}, {2,3}, {3,0}, {0,4}, {1,5}, {2,6}, {3,7}, {4,5}, {5,6}, {6,7}, {7,4}};
  static int edgemap7[][2] = {{0,1}, {1,2}, {2,0}, {0,3}, {1,4}, {2,5}, {3,4}, {4,5}, {5,6}};
  static int edgemap6[][2] = {{0,1}, {1,2}, {2,3}, {3,0}, {0,4}, {1,4}, {2,4}, {3,4}};
  static int edgemap5[][2] = {{0,1}, {1,2}, {2,0}, {0,3}, {1,3}, {2,3}};
  static int edgemap4[][2] = {{0,1}, {1,2}, {2,3}, {3,0}};
  static int edgemap3[][2] = {{0,1}, {1,2}, {2,0}};
  static int edgemap2[][2] = {{0,1}};


  // First go through elements and find all new edges introduced by the 2nd order 
  int noedges = 0;
  for(int i=0; i < mesh->getElements() + mesh->getSurfaces(); i++) {

    element_t *e;
    if(i < mesh->getElements()) 
      e = mesh->getElement(i);
    else 
      e = mesh->getSurface(i - mesh->getElements());

    int *edgemap = NULL;
    int family = e->getCode() / 100;
    int edges = familyedges[family];

    // Skip undtermined and nonlinear element 
    if((e->getNodes() < 2) || (e->getNodes() > familylinnodes[family])) continue;

    for(int f=0; f<edges; f++) {
      if(family == 3) 
	edgemap = &edgemap3[f][0];
      else if(family == 4) 
	edgemap = &edgemap4[f][0];
      if(family == 5) 
	edgemap = &edgemap5[f][0];
      else if(family == 6) 
	edgemap = &edgemap6[f][0];
      else if(family == 7) 
	edgemap = &edgemap7[f][0];
      else if(family == 8) 
	edgemap = &edgemap8[f][0];

      int n0 = e->getNodeIndex(edgemap[0]);
      int n1 = e->getNodeIndex(edgemap[1]);

      if(n0 < 0 || n1 < 0) continue;

      // Order indexes in an increasing order 
      if(n0 > n1) {
	int tmp = n0;
	n0 = n1;
	n1 = tmp;
      }

      h = &hash[n0];
      bool found = false;
      while(h->next) {                                       
	if(h->node1 == n1) {
	  found = true;
	  break;
	}
	h = h->next;
      }                                                      
      
      if(!found) {
	h->node1 = n1;
	h->noedge = noedges;

	h->next = new hashEntry;
	h = h->next;
	h->node1 = UNKNOWN;
	h->noedge = UNKNOWN;
	h->next = NULL;
	noedges++;
      }
    }
  }
  
  cout << "Found " << noedges << " edges in the mesh " << endl;
  
  if(noedges == 0) {
    return;
    delete [] hash;
  }

  // Then redifine the mesh using the additional nodes
  int quadnodes = mesh->getNodes() + noedges;  
  node_t *quadnode = new node_t[quadnodes];
  
  // Copy linear nodes
  for(int i = 0; i < mesh->getNodes(); i++) {
    quadnode[i].setX(0, mesh->getNode(i)->getX(0));
    quadnode[i].setX(1, mesh->getNode(i)->getX(1));
    quadnode[i].setX(2, mesh->getNode(i)->getX(2));
    quadnode[i].setIndex(mesh->getNode(i)->getIndex());
  }

  // Quadratic nodes are taken to be means of the linear nodes
  for(int i=0; i<keys; i++) {
    int j = 0;
    h = &hash[i];
    while(h->next){
      if(h->noedge > UNKNOWN) {
	j++;
	int n0 = i;
	int n1 = h->node1;
	int j = mesh->getNodes() + h->noedge;
	quadnode[j].setX(0, 0.5*(quadnode[n0].getX(0) + quadnode[n1].getX(0)));
	quadnode[j].setX(1, 0.5*(quadnode[n0].getX(1) + quadnode[n1].getX(1)));
	quadnode[j].setX(2, 0.5*(quadnode[n0].getX(2) + quadnode[n1].getX(2)));
	quadnode[j].setIndex(UNKNOWN);
      }
      h = h->next;     
    }
  }
  mesh->deleteNodeArray();
  mesh->setNodeArray(quadnode);

  cout << "Set " << quadnodes << " additional nodes" << endl;


  int tmpnode[27];

  for(int i=0; i < mesh->getElements() + mesh->getSurfaces() + mesh->getEdges(); i++) {
    
    element_t *e;
    if(i < mesh->getElements()) 
      e = mesh->getElement(i);
    else if (i < mesh->getElements() + mesh->getSurfaces()) 
      e = mesh->getSurface(i - mesh->getElements());
    else 
      e = mesh->getEdge(i - mesh->getElements() - mesh->getSurfaces());

    int *edgemap = NULL;
    int family = e->getCode() / 100;
    int edges = familyedges[family];

    /* Not a linear element */
    if((e->getNodes() < 2) || (e->getNodes() > familylinnodes[family])) continue;

    int linnodes = e->getNodes();
    for(int j=0;j<linnodes;j++)
      tmpnode[j] = e->getNodeIndex(j);

    e->deleteNodeIndexes();
    e->setCode(e->getCode() + edges);
    e->setNodes(e->getNodes() + edges);
    e->newNodeIndexes(e->getNodes());
    for(int j=0;j<linnodes;j++)
      e->setNodeIndex(j, tmpnode[j]);


    for(int f=0; f<edges; f++) {

      if(family == 2) 
	edgemap = &edgemap2[f][0];
      if(family == 3) 
	edgemap = &edgemap3[f][0];
      else if(family == 4) 
	edgemap = &edgemap4[f][0];
      if(family == 5) 
	edgemap = &edgemap5[f][0];
      else if(family == 6) 
	edgemap = &edgemap6[f][0];
      else if(family == 7) 
	edgemap = &edgemap7[f][0];
      else if(family == 8) 
	edgemap = &edgemap8[f][0];
      
      int n0 = e->getNodeIndex(edgemap[0]);
      int n1 = e->getNodeIndex(edgemap[1]);
      if((n0 < 0) || (n1 < 0)) continue;

      // Order indexes in an increasing order 
      if(n0 > n1) {
	int tmp = n0;
	n0 = n1;
	n1 = tmp;
      }
      
      h = &hash[n0];
      bool found = false;
      while(h->next) {                                       
	if(h->node1 == n1) {
	  found = true;
	  e->setNodeIndex(linnodes+f, h->noedge + mesh->getNodes());
	  break;
	}
	h = h->next;
      }                                                            
      if(!found) {
	cout << "All edges should be found!? " << endl;
      }
    }    
  }
  mesh->setNodes(quadnodes);

  delete [] hash;
}



// Decrease elementorder from >1 to 1
//----------------------------------------------------------------------------
void Meshutils::decreaseElementOrder(mesh_t *mesh)
{
  if((mesh->getElements() == 0) && (mesh->getSurfaces() == 0)) return;

  static int familylinnodes[9] = {0, 1, 2, 3, 4, 4, 5, 6, 8};

  int *activenodes = new int[mesh->getNodes()];
  for(int i=0; i < mesh->getNodes(); i++)
    activenodes[i] = -1;

  // int noedges = 0;
  for(int i=0; i < mesh->getElements() + mesh->getSurfaces(); i++) {

    element_t *e;
    if(i < mesh->getElements()) 
      e = mesh->getElement(i);
    else 
      e = mesh->getSurface(i - mesh->getElements());

    int family = e->getCode() / 100;
    int linnodes = familylinnodes[family];

    for(int j=0;j<linnodes;j++)
      activenodes[e->getNodeIndex(j)] = 0;
  }

  int linnodes = 0;
  for(int i=0; i < mesh->getNodes(); i++)
    if(activenodes[i] > -1) activenodes[i] = linnodes++;

  cout << "Setting " << linnodes << " linear nodes of the total " << mesh->getNodes() << " nodes" << endl;
  if(linnodes == mesh->getNodes()) return;


  // Copy linear nodes
  node_t *linnode = new node_t[linnodes];
  for(int i=0; i < mesh->getNodes(); i++) {
    int j = activenodes[i];
    if(j > -1) {   
      linnode[j].setX(0, mesh->getNode(i)->getX(0));
      linnode[j].setX(1, mesh->getNode(i)->getX(1));
      linnode[j].setX(2, mesh->getNode(i)->getX(2));
    }
  }
  mesh->deleteNodeArray();
  mesh->setNodeArray(linnode);
  mesh->setNodes(linnodes);

  int tmpnode[8];

  for(int i=0; i < mesh->getElements() + mesh->getSurfaces() + mesh->getEdges(); i++) {
    
    element_t *e;
    if(i < mesh->getElements()) 
      e = mesh->getElement(i);
    else if (i < mesh->getElements() + mesh->getSurfaces()) 
      e = mesh->getSurface(i - mesh->getElements());
    else 
      e = mesh->getEdge(i - mesh->getElements() - mesh->getSurfaces());

    int family = e->getCode() / 100;
    int linnodes = familylinnodes[family];

    for(int j=0;j<linnodes;j++)
      tmpnode[j] = e->getNodeIndex(j);

    e->deleteNodeIndexes();
    e->setCode(100*family + linnodes);
    e->setNodes(linnodes);
    e->newNodeIndexes(e->getNodes());

    for(int j=0;j<linnodes;j++)
      e->setNodeIndex(j, activenodes[tmpnode[j]]);
  }
}

// Clean up hanging sharp edges (for visualization)...
//----------------------------------------------------------------------------
int Meshutils::cleanHangingSharpEdges(mesh_t *mesh)
{
  int count_total = 0;
  int count_removed = 0;
  
  if(mesh->getEdges() == 0)
    return 0;

  for(int i = 0; i < mesh->getEdges(); i++) {
    edge_t *edge = mesh->getEdge(i);
    if(edge->isSharp()) {
      count_total++;
      if(edge->getSurfaces() == 2) {
	// has two parents
	int i0 = edge->getSurfaceIndex(0);
	int i1 = edge->getSurfaceIndex(1);
	surface_t *s0 = NULL;
	surface_t *s1 = NULL;
	if(i0 >= 0)
	  s0 = mesh->getSurface(i0);
	if(i1 >= 0)
	  s1 = mesh->getSurface(i1);
	if((s0 != NULL) && (s1 != NULL)) {
	  int index0 = s0->getIndex();
	  int index1 = s1->getIndex();
	  if(index0 == index1) {
	    // remove
	    count_removed++;
	    edge->setSharp(false);
	  }
	}
      }
    }
  }

  return count_removed;
}
