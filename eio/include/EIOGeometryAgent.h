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

#ifndef EIOGEOMETRYAGENT_H
#define EIOGEOMETRYAGENT_H

#include "EIOModelManager.h"

const int geometryFiles = 6;

class EIOGeometryAgent
{
public:
  EIOGeometryAgent(EIOModelManager *mm);
  ~EIOGeometryAgent();
  int createGeometry();
  int openGeometry();
  int closeGeometry();

  int setDescriptor(int& bodyC, int& boundaryC, int& outerC, int& innerC,
		    int& vertexC, int& maxLooplen, int& loopC);
  int descriptor(int& bodyC, int& boundaryC, int& outerC, int& innerC,
		 int& vertexC, int& maxLooplen, int& loopC);

  int writeBody(int& tag, int& meshControl, int& loopC, int *loopv);
  int nextBody(int& tag, int& meshControl, int& loopC, int *loopv);

  int writeLoop(int& tag, int& filed, int *nodes);
  int nextLoop(int& tag, int& filed, int *nodes);

  int writeElement(int& tag, int& cTag, int& meshControl, int& type,
		   int& nodeC, int *nodes);
  int nextElement(int& tag, int& cTag, int& meshControl, int& type,
		  int& nodeC, int *nodes);

  int writeNode(int& tag, int& cTag, double *coord);
  int nextNode(int& tag, int& cTag, double *coord);

  int writeBoundary(int& tag, int& left, int& right);
  int nextBoundary(int& tag, int& left, int& right);
protected:
  // We "use" ModelManager in every Agent.
  EIOModelManager *manager;

  // All streams
  fstream geometryFileStream[geometryFiles];

  // Sizes
  int bodies;
  int boundaries;
  int outer;
  int inner;
  int vertices;
  int loops;
  int maxloop;
};

#endif /* EIOGEOMETRYAGENT_H */
