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
Module:     ecif_modelOutputManager.h
Language:   C++
Date:       20.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Model output manager class. A helper class
************************************************************************/

#ifndef _ECIF_MODEL_OUTPUT_MGR_
#define _ECIF_MODEL_OUTPUT_MGR

// #include "ecif_model.h"


//enum ecifFieldInfo;

class ModelOutputManager
{     
public:
  /*  friend class Control;
  friend class Model; */

  ModelOutputManager();
  ~ModelOutputManager();

  static void initClass(Model* model);
  /* protected: */

  ostream& emf_output(ostream& out, char* filename);
  ostream& emf_outputBodies(ostream& out);
  ostream& emf_outputBoundaryConditions(ostream& out);
  ostream& emf_outputElementGroups(ostream& out);
  ostream& emf_outputEdges(ostream& out);
  ostream& emf_outputElements(ostream& out);
  ostream& emf_outputElementLoops(ostream& out);
  ostream& emf_outputEof(ostream& out);
  ostream& emf_outputFaces(ostream& out);
  ostream& emf_outputHeader(ostream& out);
  ostream& emf_outputNeighbours(ostream& out);
  ostream& emf_outputParameters(ostream& out, ecif_parameterType param_type);
  ostream& emf_outputSectionEnd(ostream& out);
  ostream& emf_outputSectionStart(ostream& out);
  ostream& emf_outputStatistics(ostream& out);
  ostream& emf_outputTimestamps(ostream& out);
  ostream& emf_outputVertexTable(ostream& out);
  ostream& emf_outputVertices(ostream& out);

  // Mesh input file
  int mif_getNofBoundaryPoints();
  ostream& mif_output(ostream& out);
  ostream& mif_outputBoundaryPoints(ostream& out);

  ostream& outputMatcDefinitions(ostream& out, int nof_defs, char** defs, bool dsign);

  // Solver input file
  ostream& sif_output(ostream& out);
  ostream& sif_outputAfterHeader(ostream& out);
  ostream& sif_outputBodies(ostream& out);
  ostream& sif_outputBodyForces(ostream& out);
  ostream& sif_outputBoundaryConditions(ostream& out);
  ostream& sif_outputBoundaries(ostream& out);
  ostream& sif_outputConstants(ostream& out);
  ostream& sif_outputEquations(ostream& out);
  ostream& sif_outputEof(ostream& out);
  ostream& sif_outputHeader(ostream& out);
  ostream& sif_outputInitialConditions(ostream& out);
  ostream& sif_outputMaterials(ostream& out);
  ostream& sif_outputSectionCode(ostream& out, const char* section_cd);
  ostream& sif_outputSectionEnd(ostream& out);
  ostream& sif_outputSectionStart(ostream& out);
  ostream& sif_outputSolvers(ostream& out);
  ostream& sif_outputSimulation(ostream& out);

  ostream& sif_outputSolverTargetFields(ostream& out, short indent_size, short indent_level, const char* source_eq_name);
  
  void write_Elmer_mesh(char* mesh_dir);
  void write_ElmerPost_mesh(ostream& outfile);
  void write_Thetis_mesh(ostream& outfile);

  static Control* theControlCenter;
  static Model* model;

};



#endif
