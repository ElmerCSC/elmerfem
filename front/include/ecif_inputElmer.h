/*****************************************************************************
 *
 *  Elmer, A Finite Element Software for Multiphysical Problems
 *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
 * 
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program (in file fem/GPL-2); if not, write to the 
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 *  Boston, MA 02110-1301, USA.
 *
 *****************************************************************************/

/***********************************************************************
Program:    ELMER Front 
Module:     ecif_inputElmer.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Elmer DB reader

************************************************************************/

#ifndef _ECIF_INPUT_ELMER 
#define _ECIF_INPUT_ELMER

#include "ecif_def.h"
#include "ecif_input.h"


struct ecf_ObjectData_X;

class InputElmer : public Input
{       
public:        
  InputElmer(enum ecif_modelDimension m_dim,
             ifstream& infile, char* filename);   
  bool readMeshGeometry();
  bool processMeshFileData();

protected:
  enum ecif_modelDimension findMeshModelDimension();
  bool readMeshHeader();
  Rc readMeshBoundaryElements();
  Rc readMeshBulkElements();
  Rc readMeshNodes(int nof_nodes);
  bool resetMeshData();
};
 
#endif

