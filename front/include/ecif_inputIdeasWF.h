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
Module:     ecif_inputIdeasWF.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Class reads I-deas wireframe-type input files. 

************************************************************************/

#ifndef _ECIF_INPUT_IDEAS_WF_ 
#define _ECIF_INPUT_IDEAS_WF_

#include "ecif_inputIdeas.h"

//*****
class InputIdeasWF : public InputIdeas
{
public:        
  InputIdeasWF(enum ecif_modelDimension m_dim,
               ifstream& in_file, char* in_filename);  
protected:
  IdNumberTable bodyNumbers;
  bool readCadGeometry(); 
  char* readBodyName(char* fileline); 
  int readBodyNbr(char* fileline); 
  int readWireFrameBody();
  ecif_geometryType readGeomType(char* s); 
  bool readLine(Body* body, char* buffer); 
  bool readNurbs(Body* body, char* buffer); 
}; 

#endif
