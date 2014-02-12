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

#include "EIODualMeshAgent.h"
#include <string.h>

extern void make_filename(char *buf, const char *model, const char *suffix);

static char *extension[] = {
  "/mesh.header",
  "/dual.elements"
};

enum { HEADER = 0, ELEMENTS};

EIODualMeshAgent::EIODualMeshAgent(EIOModelManager *mm)
{
  manager = mm;
}

EIODualMeshAgent::~EIODualMeshAgent()
{
}

void EIODualMeshAgent::
readHeader()
{
  int i;
  // Read header
  fstream& str = meshFileStream[HEADER];
  str >> nodeCount;
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
}

int EIODualMeshAgent::
createMesh(const char *dir)
{
  char filename[PATH_MAX];

  make_filename(filename, dir, extension[HEADER]);
  manager->openStream(meshFileStream[HEADER], filename, std::ios::in);

  make_filename(filename, dir, extension[ELEMENTS]);
  manager->openStream(meshFileStream[ELEMENTS], filename, std::ios::out);

  readHeader();

  return 0;
}

int EIODualMeshAgent::
openMesh(const char *dir)
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < dualMeshFiles; ++i)
    {
      make_filename(filename, dir, extension[i]);
      manager->openStream(meshFileStream[i], filename, std::ios::in);
    }

  readHeader();

  return 0;
}

int EIODualMeshAgent::
closeMesh()
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < dualMeshFiles; ++i)
    {
      manager->closeStream(meshFileStream[i]);
    } 
  return 0;
}

static int step = 0;
int EIODualMeshAgent::
read_nextElementConnections(int& tag, int& type, int* nodes)
{
  int i;
  fstream& str = meshFileStream[ELEMENTS];
  if(step == elementCount)
    {
      streampos pos = 0;
      filebuf *fbuf = str.rdbuf();
      fbuf->pubseekpos(pos, std::ios::in);
      step = 0;
      return -1;
    }
  str >> tag >> type;
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

int EIODualMeshAgent::
write_elementConnections(int& tag, int& type, int* nodes)
{
  int i;
  fstream& str = meshFileStream[ELEMENTS];
  str << tag << ' ' << type << ' ';  
  switch(type)
    {
    default:
      for(i = 0; i < 3; ++i)
	{
	  str << nodes[i] << ' ';
	}
    }
  str << std::endl;
  return 0;
}









