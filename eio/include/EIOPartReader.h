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

#ifndef EIOPARTREADER_H
#define EIOPARTREADER_H

#include "EIOModelManager.h"
using namespace std;
#include <fstream>

const int partReaderFiles = 5;

struct cache_node
{
  int tag;
  int type;
  double x,y;
};

class EIOPartReader
{
public:
  EIOPartReader(int& partCount, EIOModelManager *mm);
  ~EIOPartReader();

  int openPartitioning(int& part);
  int closePartitioning(int& part);

  // READ
  int read_descriptor(int& nodeC,
		      int& sharedC,
		      int& elementC, 
		      int& borderC, 
		      int& usedElementTypes,
		      int* elementTypeTags,
		      int* elementCountByType);

  // Reading nodes and elements is more compilated than writing them, since
  // we must cache first.

  int read_nextSharedNode(int& tag, int& partC, int *parts);
  int read_borderElements(int& len, int* tags);

  int 
  read_nextElementConnections(int& tag, int& body, int& type, int* nodes);
  int read_nextElementCoordinates(
       int& tag, int& body, int& type, int* nodes, double *coord);

private:
  // We "use" ModelManager in every Agent.
  EIOModelManager *manager;

  // All streams
  fstream meshFileStream[partReaderFiles];

  // Sizes
  char newdir[PATH_MAX];

  int parts;
  int me;

  int nodeCount;
  int sharedNodeCount;
  int elementCount;
  int boundaryElementCount;
  int elementTypes;
  int *elementTypeTags;
  int *elementTypeCount;

  // Storage
  int dim;
  cache_node *clist;

  //
  void openStreams();
  void closeStreams();
};

#endif /* EIOPARTREADER_H */
