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
Module:     ecif_inputIdeas.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   A base class for input of I-deas CAD input files. 

************************************************************************/

#ifndef _ECIF_INPUT_IDEAS 
#define _ECIF_INPUT_IDEAS

#include "ecif_input.h" 

const int NEW_DATASET = -1;
const int DS_UNITS = 164;
const int DS_COLOR = 436;
const int DS_WIRE_FRAME_CURVES = 801;
const int DS_BODIES_AND_BOUNDARIES_2430 = 2430;  // Obsolate Ideas format!
const int DS_BODIES_AND_BOUNDARIES_2432 = 2432;  // Obsolate Ideas format!
const int DS_BODIES_AND_BOUNDARIES_2435 = 2435;  // Current Ideas format!
const int DS_ELEMENTS_780 = 780;
const int DS_ELEMENTS_2412 = 2412;
const int DS_NODES_781 = 781;
const int DS_NODES_2411 = 2411;

const short MAX_IDEAS_ELEMENT_CODE = 200;
const short MAX_IDEAS_NOF_ELEMENT_SETS = 9999;


const int MESH_INDEX = -1;
const int NO_SUCH_MATERIAL  = -1;	// A 'reset' value
const int NO_SUCH_INDEX = -1;	// A 'reset' value


struct IdeasElementInfo {
  IdeasElementInfo();
  int elementSetId;
};


inline
IdeasElementInfo::IdeasElementInfo()
{
  elementSetId = -1;
}


struct IdeasElementSetInfo {
  IdeasElementSetInfo();
  ~IdeasElementSetInfo();
  int internalId;
  bool isBoundarySet;
  char* name;
  int externalBodyId;
  int externalBoundaryId;
  int nofElements;
  int nofMatchedElements;
};


inline
IdeasElementSetInfo::IdeasElementSetInfo()
{
  name = NULL;
  internalId = NO_INDEX;
  isBoundarySet = false;
  externalBodyId = NO_INDEX;
  externalBoundaryId = NO_INDEX;
  nofElements = 0;
  nofMatchedElements = 0;
}


inline
IdeasElementSetInfo::~IdeasElementSetInfo()
{
  delete[] name;
}

const int FORMWIDTH = 10;
const int DATASETLEN = 7;   // lenght of dataset number when read in from file
 
class GcPoint;
class Model;

//*****
class InputIdeas : public Input
{
public:        
  InputIdeas(enum ecif_modelDimension m_dim,
             ifstream& infile, char* filename);
  ~InputIdeas();
  bool readMeshGeometry();

protected:
  int* elementCodeCounters;
  IdeasElementInfo* elementInfos;
  IdeasElementSetInfo  elementSetInfos[MAX_IDEAS_NOF_ELEMENT_SETS];
  int maxExternalBulkElemId;
  int maxExternalBndrElemId;
  int maxMaterialId;

  //
  void countNofBoundaryAndBulkElements();
  bool endofDataset();
  bool endofDataset(char* fileline);
  int findDataset(int new_flag = 1);
  enum ecif_modelDimension findCadModelDimension();
  enum ecif_modelDimension findMeshModelDimension();
  virtual bool readCadGeometry(); 
  bool readCadHeader();
  int readColor();
  int readDatasetNbr(int newflag); 
  virtual ecif_geometryType readGeomType(char* s); 
  //bool processMeshFileData();
  bool readMeshHeader();
  virtual bool readLine(Body* body, char* buffer); 
  Rc readMeshBodiesAndBoundaries(int data_set);
  Rc readMeshElements(int data_set);
  Rc readMeshNodes(int data_set);
  Rc readNofMeshElements(int data_set, int& nof_elements, int& max_element_id, int& max_material_id);
  Rc readNofMeshNodes(int data_set, int& nof_elements, int& max_element_id);
  virtual bool readNurbs(Body* body, char* buffer); 
  virtual GcPoint* readPoint(char* fileline);
  virtual bool readPoint(char* fileline, Point3& point);
  virtual BodyElement* readVertex(char* fileline);
} ; 

#endif
