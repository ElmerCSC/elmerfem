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
Module:     ecif_inputThetis.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   A base class for input classes. 
  These classes are used for reading CAD-specific files.

************************************************************************/

#ifndef _ECIF_INPUT_THETIS 
#define _ECIF_INPUT_THETIS

#include "ecif_def.h"
#include "ecif_input.h"


// Virtual base class for input - actual class depends on type of the file to read in.

class InputThetis : public Input
{       
public:        
  InputThetis(enum ecif_modelDimension m_dim,
              ifstream& infile, char* filename);   
  virtual bool processMeshFileData();
  bool readMeshGeometry();
protected:
  enum ecif_modelDimension findModelDimension(ifstream& infile);
  bool read_boundary();
  Rc readMeshBoundaryElements(int nof_bndr_elements);
  Rc readMeshElements(int nof_elements);
  Rc readMeshNodes(int nof_nodes);
} ;
 
#endif
