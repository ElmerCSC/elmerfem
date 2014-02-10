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
Module:     ecif_inputEmf.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Elmer model file reader (.emf files)

************************************************************************/

#ifndef _ECIF_INPUT_EMF 
#define _ECIF_INPUT_EMF

#include "ecif_def.h"
#include "ecif_input.h"
#include "ecif_inputFront.h"
 
struct emf_ObjectData_X;

class InputEmf : public InputFront
{       
public:        
  InputEmf(enum ecif_modelDimension m_dim, ifstream& infile, char* filename);   
  void copyEmfParameters();
protected:
  static int copyEmfParametersCallBack(void** user_data);
  bool readCadGeometry();
  static int readEmfGeometryCallBack(void** user_data);
  static int readEmfGeometryMsgCallBack(char* msg_buffer);
  //static int readEmfElement(InputEmf* emf_input, emf_ObjectData_X* object_data, ecif_topologyType element_type);
  //static int readEmfElementLoop(InputEmf* emf_input, emf_ObjectData_X* object_data);
  static int readEmfHeader(InputEmf* emf_input, emf_ObjectData_X* object_data);
  static int readEmfParameter(emf_ObjectData_X* object_data,
                              ecif_parameterType parameter_type,
                              const Parameter*& param,
                              bool add_to_model);
  static int readEmfTimestamps(emf_ObjectData_X* object_data);

};
 
#endif
