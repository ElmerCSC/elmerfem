/*  
   Elmer, A Finite Element Software for Multiphysical Problems
   Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/***********************************************************************
Program:    ELMER Data base interface (EIO)
Author(s):  Harri Hakula 10.03.98
************************************************************************/

#ifndef EIOMESHAGENT_H
#define EIOMESHAGENT_H

#include "EIOModelManager.h"

struct cache_node
{
  int tag;
  int constraint;
  double x,y,z;
};


class EIOMeshAgent
{
public:
  EIOMeshAgent(EIOModelManager *mm, int split=0, int part=0);
  ~EIOMeshAgent();

  int createMesh(const char *dir);
  int openMesh(const char *dir);
  int closeMesh();
  
  // READ
  int read_descriptor(int& nodeC, int& elementC, int& boundaryElementC, 
		      int& usedElementTypes, int* usedElementTypeTags,
		      int* usedElementTypeCount);
  int read_nextElementConnections(int& tag, int& part, int& body, int& type, int *pdofs, int* nodes);
  int read_nextElementCoordinates(int& tag, int& body, int& type, int* nodes,
			     double *coord);
  int read_nextBoundaryElement(int& tag, int& part, int& boundary,
                               int& leftElement, int& rightElement,
                               int& type, int* nodes, double* coord);
  int read_allNodes(int *tags,double* coord);

  // WRITE
  int write_descriptor(int& nodeC, int& elementC, int& boundaryElementC, 
		       int& usedElementTypes,
		       int* elementTypeTags,
		       int* elementCountByType);
  int write_node(int& tag, int& type, double* coord);
  int write_elementConnections(int& tag, int& body, int& type, int* nodes);
  int write_boundaryElement(int& tag, int& boundary, 
			    int& leftElement, int& rightElement, 
			    int& type, int* nodes);
  // NEW
  int read_partDescriptor(int& shared);
  int read_sharedNode(int& tag, 
			 int& constraint,      
			 double *coord, 
			 int& partcount, 
			 int *parts);

private:
  // We "use" ModelManager in every Agent.
  EIOModelManager *manager;

  // All streams
  fstream *meshFileStream;
  // Sizes
  // Sizes
  char newdir[PATH_MAX];

  int parts;
  int me;

  int nodeCount;
  int elementCount;
  int boundaryElementCount;
  int elementTypes;
  int *elementTypeTags;
  int *elementTypeCount;

  int sharedNodeCount;
  int borderElementCount;

  // Storage
  cache_node *clist;

  int dim;

  // Setting
  int parallel;
  int meshFiles;

  void cache_nodes();

  int copy_coords(double *target, const int address);
  cache_node * search_node(const int address);
};

#endif /* EIOMESHAGENT_H */
