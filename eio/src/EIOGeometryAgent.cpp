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

#include "EIOGeometryAgent.h"

extern void make_filename(char *buf, const char *model, const char *suffix);

static const char *extension[] = {
  "geometry.header",
  "geometry.nodes",
  "geometry.elements",
  "geometry.bodies",
  "geometry.loops",
  "geometry.boundaries"
};

enum { HEADER = 0, NODES, ELEMENTS, BODIES, LOOPS, BOUNDARIES };

EIOGeometryAgent::EIOGeometryAgent(EIOModelManager *mm)
{
  manager = mm;
}

EIOGeometryAgent::~EIOGeometryAgent()
{
}

int EIOGeometryAgent::
createGeometry()
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < geometryFiles; ++i)
    {
      //      make_filename(filename, manager->name(), extension[i]);
      manager->openStream(geometryFileStream[i], extension[i], std::ios::out);
    }
  return 0;
}

int EIOGeometryAgent::
openGeometry()
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < geometryFiles; ++i)
    {
      //      make_filename(filename, manager->name(), extension[i]);
      manager->openStream(geometryFileStream[i], extension[i], std::ios::in);
    }

  // Read header
  int ftypes;
  fstream& str = geometryFileStream[HEADER];
  str >> bodies;
  str >> boundaries;
  str >> outer;
  str >> inner;
  str >> vertices;
  str >> loops;
  str >> maxloop;
  return 0;
}

int EIOGeometryAgent::
closeGeometry()
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < geometryFiles; ++i)
    {
      manager->closeStream(geometryFileStream[i]);
    } 
  return 0;
}

int EIOGeometryAgent::
descriptor(int& bodyC, int& boundaryC, int& outerC, int& innerC,
	   int& vertexC, int& maxLooplen, int& loopC)
{
  bodyC = bodies;
  boundaryC = boundaries;
  outerC = outer;
  innerC = inner;
  vertexC = vertices;
  loopC = loops;
  maxLooplen = maxloop;
  return 0;
}

int EIOGeometryAgent::
setDescriptor(int& bodyC, int& boundaryC, int& outerC, int& innerC,
	   int& vertexC, int& maxLooplen, int& loopC)
{
  bodies = bodyC;
  boundaries = boundaryC;
  outer = outerC;
  inner = innerC;
  vertices = vertexC;
  loops = loopC; 
  maxloop = maxLooplen;

  fstream& str = geometryFileStream[HEADER];
  str << bodies << ' '
      << boundaries << ' '
      << outer << ' '
      << inner << ' '
      << vertices << ' '
      << loops << ' '
      << maxloop << ' '
      << std::endl;
  return 0;
}

int EIOGeometryAgent::
writeNode(int& tag, int& cTag, double *coord)
{
  int i;
  fstream& str = geometryFileStream[NODES];
  str << tag << ' '
      << cTag << ' ';
  for(i = 0; i < 3; ++i)
    {
      str << coord[i] << ' ';
    }
  str << std::endl;
  return 0;
}

int EIOGeometryAgent::
writeBody(int& tag, int& meshControl, int& loopC, int *loopv) 
{ 
  int i;
  fstream& str = geometryFileStream[BODIES];
  str << tag << ' '
      << meshControl << ' '
      << loopC << '\n';

  for(i = 0; i < loopC; ++i)
    {
      str << loopv[i] << ' ';
    }
  str << std::endl;
  return 0;
}

static int step = 0;
int EIOGeometryAgent::
nextBody(int& tag, int& meshControl, int& loopC, int *loopv) 
{ 
  int i;
  fstream& str = geometryFileStream[BODIES];
  if(step == bodies)
    {
      step = 0;
      return -1;
    }
  str >> tag >> meshControl >> loopC;
  for(i = 0; i < loopC; ++i)
    {
      str >> loopv[i];
    }
  ++step;
  return 0;
}

int EIOGeometryAgent::
writeLoop(int& tag, int& field, int *nodes)
{
  int i;
  fstream& str = geometryFileStream[LOOPS];
  str << tag << ' '
      << field << ' ';
  for(i = 0; i < field; ++i)
    {
      str << nodes[i] << ' ';
    }
  str << std::endl;  
  return 0;
}

static int lstep = 0;
int EIOGeometryAgent::
nextLoop(int& tag, int& field, int *nodes)
{
  int i;
  fstream& str = geometryFileStream[LOOPS];

  if(lstep == maxloop)
    {      
      streampos pos = 0;
      filebuf *fbuf = str.rdbuf();
      fbuf->pubseekpos(pos, std::ios::in);
      lstep = 0;
      return -1;
    }
  str >> tag >> field;
  for(i = 0; i < field; ++i)
    {
      str >> nodes[i];
    }
  ++lstep;
  return 0;
}

// Modified to use nodeC: Martti Verho, 17.03.99
int EIOGeometryAgent::
writeElement(int& tag, int& cTag, int& meshControl, int& type,
	     int& nodeC, int *nodes) 
{ 
  int i;
  fstream& str = geometryFileStream[ELEMENTS];
  str << tag << ' '
      << cTag << ' '
      << meshControl << ' '
      << type << ' ';

  switch(type)
    {
    // No nodeC for 2-vertex edges!
    case 101:
      break;
    // Multi vertex edge, nodeC needed
    default:
      str << nodeC << ' ';
      break;
    }

  for(i = 0; i < nodeC; ++i)
    {
      str << nodes[i] << ' ';
    }
  str << std::endl;

  return 0;
}

// Modified to use nodeC: Martti Verho, 17.03.99
int EIOGeometryAgent::
nextElement(int& tag, int& cTag, int& meshControl, int& type,
	    int& nodeC, int *nodes) 
{ 
  fstream& str = geometryFileStream[ELEMENTS];
  if(step == boundaries)
    {
      step = 0;
      return -1;
    }

  str >> tag >> cTag >> meshControl >> type;

  switch(type)
    {
    // No nodeC stored in db for a 2-vertex edge
    case 101:
      nodeC = 2;	      
      break;
    // Read nodeC from db for a multi verex edge
    default:
      str >> nodeC;
      break;
    }

  int junk;
  // Deliver node tags only if nodes-buffer really given
  for(int i = 0; i < nodeC; i++)
    {
      if ( nodes != NULL )
	str >> nodes[i];
      else
        str >> junk;
    }

  ++step;
  return 0;
}

int EIOGeometryAgent::
nextNode(int& tag, int& cTag, double *coord) 
{ 
  fstream& str = geometryFileStream[NODES];
  if(step == vertices)
    {
      step = 0;
      return -1;
    }
  str >> tag >> cTag >> coord[0] >> coord[1] >> coord[2];
  ++step;
  return 0;
}

int EIOGeometryAgent::
writeBoundary(int& tag, int& left, int& right) 
{ 
  fstream& str = geometryFileStream[BOUNDARIES];
  str << tag << ' '
      << left << ' '
      << right << std::endl;
  return 0;
}

int EIOGeometryAgent::
nextBoundary(int& tag, int& left, int& right) 
{ 
  fstream& str = geometryFileStream[BOUNDARIES];
  if(step == (inner + outer))
    {
      step = 0;
      return -1;
    }
  str >> tag >> left >> right;
  ++step;
  return 0;
}
