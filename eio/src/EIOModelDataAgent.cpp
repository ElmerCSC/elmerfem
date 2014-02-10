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

#include "EIOModelDataAgent.h"
#include <string.h>

extern
void make_filename(char *buf, const char *model, const char *suffix);

void eio_output_head(fstream& str, int tag, int fields)
{
  str << tag << ' ' << fields << '\n';
}

void eio_input_head(fstream& str, int& tag, int& fields)
{
  str >> tag >> fields;
}

void eio_output_field(fstream& str, 
		      int& name, int& type, 
		      int& len, int* fields, double* values)
{
  int i = 0;
  str << name << ' '
      << type << ' '
      << len << ' ';
  for(i = 0; i < len; ++i)
    {
      str << fields[i] << ' ';
    }
  for(i = 0; i < len; ++i)
    {
      str << values[i] << ' ';
    }
  str << '\n';
}

void eio_input_field(fstream& str, 
		     int& name, int& type, 
		     int& len, int* fields, double* values)
{
  int i;
  str >> name >> type >> len;
  for(i = 0; i < len; ++i)
    {
      str >> fields[i];
    }
  for(i = 0; i < len; ++i)
    {
      str >> values[i];
    }
}


// File list changed
// Martti Verho 19.10.98
#if 0
static char *extension[] = {
  "modeldata.header",
  "modeldata.materials",
  "modeldata.boundary_conditions",
  "modeldata.initial_conditions",
  "modeldata.body_equations",
  "modeldata.body_forces",
  "modeldata.mesh_parameters",
  "modeldata.constants",
  "modeldata.coordinates",
  "modeldata.bodies"
};
enum { HEADER = 0, MATERIALS, BOUNDARYCONDITIONS,
INITIALCONDITIONS, BODYEQUATIONS, BODYFORCES, MESHPARAMETERS,
CONSTANTS, COORDINATES, BODIES};
#endif

// Only first three are needed currently
// Ref. the value of modelDataFiles in EIOModelDataAgent.h
// Martti Verho 19.10.98
static char *extension[] = {
  "modeldata.header",
  "modeldata.coordinates",
  "modeldata.mesh_parameters",
  "modeldata.materials",
  "modeldata.boundary_conditions",
  "modeldata.initial_conditions",
  "modeldata.body_equations",
  "modeldata.body_forces",
  "modeldata.constants",
  "modeldata.bodies"
};

enum { HEADER = 0, COORDINATES,  MESHPARAMETERS,
MATERIALS, BOUNDARYCONDITIONS,
INITIALCONDITIONS, BODYEQUATIONS, BODYFORCES,
CONSTANTS, BODIES};

EIOModelDataAgent::EIOModelDataAgent(EIOModelManager *mm)
{
  manager = mm;
}

EIOModelDataAgent::~EIOModelDataAgent()
{
}

int EIOModelDataAgent::
createModelData()
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < modelDataFiles; ++i)
    {
      //      make_filename(filename, manager->name(), extension[i]);
      manager->openStream(modelDataFileStream[i], extension[i], std::ios::out);
    }

  return 0;
}

int EIOModelDataAgent::
openModelData()
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < modelDataFiles; ++i)
    {
      //      make_filename(filename, manager->name(), extension[i]);
      manager->openStream(modelDataFileStream[i], extension[i], std::ios::in);
    }

  return 0;
}

int EIOModelDataAgent::
closeModelData()
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < modelDataFiles; ++i)
    {
      manager->closeStream(modelDataFileStream[i]);
    } 
  return 0;
}
int EIOModelDataAgent::
writeDescription(int& bodies,
		 int& body_forces,
		 int& body_equations,
		 int& materials,
		 int& boundary_conditions,
		 int& initial_conditions,
		 int& mesh_parameters)
{
  fstream& str = modelDataFileStream[HEADER];
  str << bodies << ' ' 
      << body_forces << ' '
      << body_equations << ' '
      << materials << ' '
      << boundary_conditions << ' '
      << initial_conditions << ' '
      << mesh_parameters << std::endl;
  return 0;
}

int EIOModelDataAgent::
readDescription(int& bodies,
		int& body_forces,
		int& body_equations,
		int& materials,
		int& boundary_conditions,
		int& initial_conditions,
		int& mesh_parameters)
{
  fstream& str = modelDataFileStream[HEADER];
  str >> bodies
      >> body_forces
      >> body_equations
      >> materials
      >> boundary_conditions
      >> initial_conditions
      >> mesh_parameters;
  return 0;
}

int EIOModelDataAgent::
writeBodyRecord(int& tag, int& body_force_id,
			  int& equation_id, int& init_cond_id,
			  int& material_id, int& mesh_param_id)
{
  fstream& str = modelDataFileStream[BODIES];
  str << tag << ' '
      << body_force_id << ' '
      << equation_id << ' '
      << init_cond_id << ' '
      << material_id << ' '
      << mesh_param_id << std::endl;
  return 0;
}

int EIOModelDataAgent::
readBodyRecord(int& tag, int& body_force_id,
			  int& equation_id, int& init_cond_id,
			  int& material_id, int& mesh_param_id)
{
  fstream& str = modelDataFileStream[BODIES];
  str >> tag
      >> body_force_id
      >> equation_id
      >> init_cond_id
      >> material_id
      >> mesh_param_id;
  return 0;
}

int EIOModelDataAgent::
writeConstants(double* gravity, double& boltz)
{
  int i;
  fstream& str = modelDataFileStream[CONSTANTS];
  for(i = 0; i < 4; ++i)
    {
      str << gravity[i] << std::endl;
    }
  str << boltz << std::endl;
  return 0;
}

int EIOModelDataAgent::
readConstants(double* gravity, double& boltz)
{
  int i;
  fstream& str = modelDataFileStream[CONSTANTS];
  for(i = 0; i < 4; ++i)
    {
      str >> gravity[i];
    }
  str >> boltz;
  return 0;
}

int EIOModelDataAgent::
writeCoordinates(int& dim, int& coordsys, int *mapping,
				  int& symmetry,
				  double *start,
				  double *end1, double* end2)
{
  int i;
  fstream& str = modelDataFileStream[COORDINATES];
  str << dim << ' ' << coordsys << ' ';
  for(i = 0; i < 3; ++i)
    {
      str << mapping[i] << ' ';
    }
  str << std::endl;
  str << symmetry << std::endl;
  for(i = 0; i < 3; ++i)
    {
      str << start[i] << ' ';
    }
  str << std::endl;
   for(i = 0; i < 3; ++i)
    {
      str << end1[i] << ' ';
    }
  str << std::endl;
  for(i = 0; i < 3; ++i)
    {
      str << end1[i] << ' ';
    }
  str << std::endl;
 
  return 0;
}

int EIOModelDataAgent::
readCoordinates(int& dim, int& coordsys, int *mapping,
				  int& symmetry,
				  double *start,
				  double *end1, double* end2)
{
  int i;
  fstream& str = modelDataFileStream[COORDINATES];
  str >> dim >> coordsys;
  for(i = 0; i < 3; ++i)
    {
      str >> mapping[i];
    }
  str >> symmetry;
  for(i = 0; i < 3; ++i)
    {
      str >> start[i];
    }
   for(i = 0; i < 3; ++i)
    {
      str >> end1[i];
    }
  for(i = 0; i < 3; ++i)
    {
      str >> end1[i];
    }
  return 0;
}

int EIOModelDataAgent::
writeMaterialHead(int& tag, int& fields)
{
  fstream& str = modelDataFileStream[MATERIALS];
  eio_output_head(str, tag, fields);
  return 0;
}

int EIOModelDataAgent::
writeMaterialField(int& name, int& type, int& len, int* fields, double* values)
{
  fstream& str = modelDataFileStream[MATERIALS];
  eio_output_field(str, name, type, len, fields, values);
  return 0;
}

int EIOModelDataAgent::
readMaterialHead(int& tag, int& fields)
{
  fstream& str = modelDataFileStream[MATERIALS];
  eio_input_head(str, tag, fields);
  return 0;
}

int EIOModelDataAgent::
readMaterialField(int& name, int& type, int& len, int* fields, double* values)
{
  fstream& str = modelDataFileStream[MATERIALS];
  eio_input_field(str, name, type, len, fields, values);
  return 0;
}

int EIOModelDataAgent::
writeBoundaryConditionHead(int& tag, int& fields)
{
  fstream& str = modelDataFileStream[BOUNDARYCONDITIONS];
  eio_output_head(str, tag, fields);
  return 0;
}

int EIOModelDataAgent::
writeBoundaryConditionField(int& name, int& type, int& len,
			    int* fields, double* values)
{
  fstream& str = modelDataFileStream[BOUNDARYCONDITIONS];
  eio_output_field(str, name, type, len, fields, values);
  return 0;
}

int EIOModelDataAgent::
readBoundaryConditionHead(int& tag, int& fields)
{
  fstream& str = modelDataFileStream[BOUNDARYCONDITIONS];
  eio_input_head(str, tag, fields);
  return 0;
}

int EIOModelDataAgent::
readBoundaryConditionField(int& name, int& type, int& len, int* fields, double* values)
{
  fstream& str = modelDataFileStream[BOUNDARYCONDITIONS];
  eio_input_field(str, name, type, len, fields, values);
  return 0;
}

int EIOModelDataAgent::
writeInitialConditionHead(int& tag, int& fields)
{
  fstream& str = modelDataFileStream[INITIALCONDITIONS];
  eio_output_head(str, tag, fields);
  return 0;
}

int EIOModelDataAgent::
writeInitialConditionField(int& name, int& type, int& len,
			    int* fields, double* values)
{
  fstream& str = modelDataFileStream[INITIALCONDITIONS];
  eio_output_field(str, name, type, len, fields, values);
  return 0;
}

int EIOModelDataAgent::
readInitialConditionHead(int& tag, int& fields)
{
  fstream& str = modelDataFileStream[INITIALCONDITIONS];
  eio_input_head(str, tag, fields);
  return 0;
}

int EIOModelDataAgent::
readInitialConditionField(int& name, int& type, int& len, int* fields, double* values)
{
  fstream& str = modelDataFileStream[INITIALCONDITIONS];
  eio_input_field(str, name, type, len, fields, values);
  return 0;
}

int EIOModelDataAgent::
writeBodyEquationHead(int& tag, int& fields)
{
  fstream& str = modelDataFileStream[BODYEQUATIONS];
  eio_output_head(str, tag, fields);
  return 0;
}

int EIOModelDataAgent::
writeBodyEquationField(int& name, int& type, int& len,
			    int* fields, double* values)
{
  fstream& str = modelDataFileStream[BODYEQUATIONS];
  eio_output_field(str, name, type, len, fields, values);
  return 0;
}

int EIOModelDataAgent::
readBodyEquationHead(int& tag, int& fields)
{
  fstream& str = modelDataFileStream[BODYEQUATIONS];
  eio_input_head(str, tag, fields);
  return 0;
}

int EIOModelDataAgent::
readBodyEquationField(int& name, int& type, int& len, int* fields, double* values)
{
  fstream& str = modelDataFileStream[BODYEQUATIONS];
  eio_input_field(str, name, type, len, fields, values);
  return 0;
}

int EIOModelDataAgent::
writeBodyForceHead(int& tag, int& fields)
{
  fstream& str = modelDataFileStream[BODYFORCES];
  eio_output_head(str, tag, fields);
  return 0;
}

int EIOModelDataAgent::
writeBodyForceField(int& name, int& type, int& len,
			    int* fields, double* values)
{
  fstream& str = modelDataFileStream[BODYFORCES];
  eio_output_field(str, name, type, len, fields, values);
  return 0;
}

int EIOModelDataAgent::
readBodyForceHead(int& tag, int& fields)
{
  fstream& str = modelDataFileStream[BODYFORCES];
  eio_input_head(str, tag, fields);
  return 0;
}

int EIOModelDataAgent::
readBodyForceField(int& name, int& type, int& len, int* fields, double* values)
{
  fstream& str = modelDataFileStream[BODYFORCES];
  eio_input_field(str, name, type, len, fields, values);
  return 0;
}

int EIOModelDataAgent::
writeMeshParameterHead(int& tag, int& fields)
{
  fstream& str = modelDataFileStream[MESHPARAMETERS];
  eio_output_head(str, tag, fields);
  return 0;
}

int EIOModelDataAgent::
writeMeshParameterField(int& name, int& type, int& len,
			    int* fields, double* values)
{
  fstream& str = modelDataFileStream[MESHPARAMETERS];
  eio_output_field(str, name, type, len, fields, values);
  return 0;
}

int EIOModelDataAgent::
readMeshParameterHead(int& tag, int& fields)
{
  fstream& str = modelDataFileStream[MESHPARAMETERS];
  eio_input_head(str, tag, fields);
  return 0;
}

int EIOModelDataAgent::
readMeshParameterField(int& name, int& type, int& len, int* fields, double* values)
{
  fstream& str = modelDataFileStream[MESHPARAMETERS];
  eio_input_field(str, name, type, len, fields, values);
  return 0;
}


