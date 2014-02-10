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

#ifndef EIODUALMESHAGENT_H
#define EIODUALMESHAGENT_H

#include "EIOModelManager.h"

const int dualMeshFiles = 2;

class EIODualMeshAgent
{
public:
  EIODualMeshAgent(EIOModelManager *mm);
  ~EIODualMeshAgent();

  int createMesh(const char *dir);
  int openMesh(const char *dir);
  int closeMesh();
  int read_nextElementConnections(int& tag, int& type, int* nodes);
  int write_elementConnections(int& tag, int& type, int* nodes);

private:
  // We "use" ModelManager in every Agent.
  EIOModelManager *manager;

  // All streams
  fstream meshFileStream[dualMeshFiles];

  // Sizes
  int nodeCount;
  int elementCount;
  int boundaryElementCount;
  int elementTypes;
  int *elementTypeTags;
  int *elementTypeCount;

  void readHeader();
};

#endif /* EIODUALMESHAGENT_H */
