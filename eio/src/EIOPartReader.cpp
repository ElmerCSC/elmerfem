/*  
   Elmer, A Finite Element Software for Multiphysical Problems
  
   Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
   
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
  
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with this library (in file ../../LGPL-2.1); if not, write 
   to the Free Software Foundation, Inc., 51 Franklin Street, 
   Fifth Floor, Boston, MA  02110-1301  USA
*/

/***********************************************************************
Program:    ELMER Data base interface (EIO)
Author(s):  Harri Hakula 10.03.98
************************************************************************/

#include "EIOPartReader.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "eio-config.h"
 
extern void make_filename(char *buf, const char *model, const char *suffix);

int CDECL nodecomp(const void *a, const void *b);

static char *extension[] = {
  "%s/part.%d.header",
  "%s/part.%d.nodes",
  "%s/part.%d.shared",
  "%s/part.%d.elements",
  "%s/part.%d.border"
};

enum { HEADER = 0, NODES, SHARED, ELEMENTS, BORDER};

EIOPartReader::EIOPartReader(int& partCount, EIOModelManager *mm)
{
  parts = partCount;
  me = -1;
  manager = mm;
  dim = 2;
  clist = (cache_node *)0;
}

EIOPartReader::~EIOPartReader()
{
}

int EIOPartReader::openPartitioning(int& part)
{
  int i;
  me = part;
  sprintf(newdir, "%s/partitioning.%d", meshdir, parts);

  openStreams();

  // Read header
  fstream& str = meshFileStream[HEADER];
  str >> nodeCount;
  str >> sharedNodeCount;
  str >> elementCount;
  str >> boundaryElementCount;
  str >> elementTypes;
  
  elementTypeTags = new int[elementTypes];
  elementTypeCount = new int[elementTypes];

  for(i = 0; i < elementTypes; ++i)
    {
      int etype, ecount;
      str >> etype >> ecount;
      elementTypeTags[i] = etype;
      elementTypeCount[i] = ecount;
    }

  return 0;
}

int EIOPartReader::closePartitioning(int& part)
{
  me = -1;
  closeStreams();
  return 0;
}

int EIOPartReader::read_descriptor(int& nodeC,                 /* nodes */
		 int& sharedC,               /* shared nodes */
		 int& elementC,              /* elements (inner) */
		 int& borderC,               /* elements bordering the part */
		 int& usedElementTypes, 
		 int* usedElementTypeTags,
		 int* usedElementCountByType)
{
  int i;
  nodeC = nodeCount;
  sharedC = sharedNodeCount;
  elementC = elementCount;
  borderC = boundaryElementCount;
  usedElementTypes = elementTypes;

  for(i = 0; i < elementTypes; ++i)
    {
      usedElementTypeTags[i] = elementTypeTags[i];
      usedElementCountByType[i] = elementTypeCount[i];
    }
  return 0;
}

int EIOPartReader::read_borderElements(int& len, int* tags)
{
  int i;
  fstream& str = meshFileStream[SHARED];
  
  for(i = 0; i < boundaryElementCount; ++i)
    {
      str >> tags[i];
    }

  streampos pos = 0;
  filebuf *fbuf = str.rdbuf();
  fbuf->seekpos(pos, ios::in);
  return 0;
}

static int step = 0;
int EIOPartReader::read_nextSharedNode(int& tag, int& partC, int *parts)
{
  int i;
  fstream& str = meshFileStream[SHARED];
  if(step == sharedNodeCount)
    {
      streampos pos = 0;
      filebuf *fbuf = str.rdbuf();
      fbuf->seekpos(pos, ios::in);
      step = 0;
      return -1;
    }
  str >> tag;
  str >> partC;
  for(i = 0; i < partC; ++i)
    {
      str >> parts[i];
    }
  return 0;
}

int EIOPartReader::read_nextElementConnections(int& tag, int& body, int& type, int* nodes)
{
  int i;
  fstream& str = meshFileStream[ELEMENTS];
  if(step == elementCount)
    {
      streampos pos = 0;
      filebuf *fbuf = str.rdbuf();
      fbuf->seekpos(pos, ios::in);
      step = 0;
      return -1;
    }
  str >> tag >> body >> type;
  switch(type)
    {
    default:
      for(i = 0; i < 3; ++i)
	{
	  str >> nodes[i];
	}
    }
  ++step;
  return 0;
}

int EIOPartReader::read_nextElementCoordinates(
       int& tag, int& body, int& type, int* nodes, double *coord)
{
  int i;
  fstream& str = meshFileStream[ELEMENTS];
  if(step == 0)
    {
      if(!clist)
	{
	  clist = new cache_node[nodeCount];
	  fstream& nstr = meshFileStream[NODES];
	  for(i = 0; i < nodeCount; ++i)
	    {
	      nstr >> clist[i].tag >> clist[i].type 
		   >> clist[i].x >> clist[i].y;
	    }
	  streampos pos = 0;
	  filebuf *fbuf = nstr.rdbuf();
	  fbuf->seekpos(pos, ios::in);
	}
    }
  else if(step == elementCount)
    {
      streampos pos = 0;
      step = 0;
      filebuf *fbuf = str.rdbuf();
      fbuf->seekpos(pos, ios::in);
      return -1;
    }

  str >> tag >> body >> type;
  switch(type)
    {
    default:
      for(i = 0; i < 3; ++i)
	{
	  cache_node entry;
	  cache_node *retval;

	  str >> nodes[i];
	  entry.tag = nodes[i];
	  retval = (cache_node *) bsearch((void *)&entry, (void *)clist, 
					  nodeCount, 
					  sizeof(cache_node),
					  nodecomp);
	  if(retval == NULL) 
	  {
	      cout << "PANIC!" << endl;
	      exit(14);
	  }
	  coord[i*dim] = retval->x;
	  coord[i*dim+1] = retval->y;
	}
    }
  ++step;
  return 0;  
}

void EIOPartReader::openStreams()
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < partReaderFiles; ++i)
    {
      sprintf(filename, extension[i], newdir, me);
      manager->openStream(meshFileStream[i], filename, ios::in);
    }
}

void EIOPartReader::closeStreams()
{
  int i;
  for(i = 0; i < partReaderFiles; ++i)
    {
      manager->closeStream(meshFileStream[i]);
    }
}
