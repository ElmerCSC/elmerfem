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

#include "EIOModelManager.h"
#include "EIOGeometryAgent.h"
#include "EIOMeshAgent.h"
#include "EIOSolverAgent.h"
#include "EIOModelDataAgent.h"
#include "EIODualMeshAgent.h"
#include "EIOPartWriter.h"

#include "../config.h"

struct parallel
{
  int isParallel;
  int numProc;
  int myProc;
};

static struct parallel paraState = { 0, 0, 0};
static EIOModelManager *modelManager = 0;
static EIOGeometryAgent *geometryAgent = 0;
static EIOMeshAgent *meshAgent = 0;
static EIOSolverAgent *solverAgent = 0;
static EIOModelDataAgent *modelDataAgent = 0;
static EIODualMeshAgent *dualMeshAgent = 0;
static EIOPartWriter *partitioningWriter = 0;

//#define EIOFC(funname) extern "C" void WITH_BINDING_EXTENSION(funname)

extern "C" void STDCALLBULL FC_FUNC_(eio_init,EIO_INIT)
  (int& info)
{
  paraState.isParallel = 0;
  paraState.numProc    = 1;
  paraState.myProc     = 0;

  if(modelManager = new EIOModelManager)
    {
      info = 0;
    }
  else
    {
      info = -1;
    }
}

extern "C" void STDCALLBULL FC_FUNC_(eio_init_parallel,EIO_INIT_PARALLEL)
  (int& procs, int& me, int& info)
{
  paraState.isParallel = 1;
  paraState.numProc    = procs;
  paraState.myProc     = me;

  if(modelManager = new EIOModelManager)
    {
      info = 0;
    }
  else
    {
      info = -1;
    }
}

extern "C" void STDCALLBULL FC_FUNC_(eio_close,EIO_CLOSE)
  (int& info)
{
  delete modelManager;
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_create_model,EIO_CREATE_MODEL)
  (const char *directory, int& info)
{
  info = modelManager->createModel(directory);
}

extern "C" void STDCALLBULL FC_FUNC_(eio_open_model,EIO_OPEN_MODEL)
  (const char *directory, int& info)
{
  info = modelManager->openModel(directory);
}

extern "C" void STDCALLBULL FC_FUNC_(eio_close_model,EIO_CLOSE_MODEL)
  (int& info)
{
  info = modelManager->closeModel();
}

extern "C" void STDCALLBULL FC_FUNC_(eio_create_geometry,EIO_CREATE_GEOMETRY)
  (int& info)
{
  if(geometryAgent = new EIOGeometryAgent(modelManager))
    {
      info = geometryAgent->createGeometry();
    }
  else
    {
      info = -1;
    }
}

extern "C" void STDCALLBULL FC_FUNC_(eio_open_geometry,EIO_OPEN_GEOMETRY)
  (int& info)
{
  if(geometryAgent = new EIOGeometryAgent(modelManager))
    {
      info = geometryAgent->openGeometry();
    }
  else
    {
      info = -1;
    }
}

extern "C" void STDCALLBULL FC_FUNC_(eio_close_geometry,EIO_CLOSE_GEOMETRY)
  (int& info)
{
  geometryAgent->closeGeometry();
  delete geometryAgent;
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_geometry_description,EIO_SET_GEOMETRY_DESCRIPTION)
  (int& bodyC, int& boundaryC, int& outerC, 
   int& innerC, int& vertexC, 
   int& loopC, int& maxLooplen, int& info)
{
  geometryAgent->setDescriptor(bodyC, boundaryC, outerC, innerC, vertexC,
			       loopC, maxLooplen);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_geometry_description,EIO_GET_GEOMETRY_DESCRIPTION)
  (int& bodyC, int& boundaryC, int& outerC, 
   int& innerC, int& vertexC, 
   int& loopC, int& maxLooplen, int& info)
{
  geometryAgent->descriptor(bodyC, boundaryC, outerC, innerC, vertexC,
			    loopC, maxLooplen);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_geometry_body,EIO_SET_GEOMETRY_BODY)
  (int& tag, int& meshControl, int& loopC,
   int *loops,
   int& info)
{
  geometryAgent->writeBody(tag, meshControl, loopC, loops);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_geometry_body,EIO_GET_GEOMETRY_BODY)
  (int& tag, int& meshControl, int& loopC,
   int *loops,
   int& info)
{
  if(geometryAgent->nextBody(tag, meshControl, loopC, loops) != -1)
    info = 0;
  else
    info = -1;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_geometry_body_loop,EIO_SET_GEOMETRY_BODY_LOOP)
  (int& tag, int& field, int *nodes, int& info)
{
  geometryAgent->writeLoop(tag, field, nodes);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_geometry_body_loop,EIO_GET_GEOMETRY_BODY_LOOP)
  (int& tag, int& field, int *nodes, int& info)
{
  if(geometryAgent->nextLoop(tag, field, nodes) != -1)
    info = 0;
  else
    info = -1;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_geometry_element,EIO_SET_GEOMETRY_ELEMENT)
  (int& tag, int& cTag, int& meshControl,
   int& type, int& nodeC, int *nodes, int& info)
{
  // Added nodeC argument: Martti Verho 1703.99
  geometryAgent->writeElement(tag, cTag, meshControl, type, nodeC, nodes);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_geometry_element,EIO_GET_GEOMETRY_ELEMENT)
  (int& tag, int& cTag, int& meshControl,
   int& type, int& nodeC, int *nodes, int& info)
{
  // Added nodeC argument: Martti Verho 1703.99
  if(geometryAgent->nextElement(tag, cTag, meshControl, type, nodeC, nodes) != -1)
    info = 0;
  else
    info = -1;
}

// Added: Martti Verho, 17.03.99
extern "C" void STDCALLBULL FC_FUNC_(eio_get_geometry_element_description,EIO_GET_GEOMETRY_ELEMENT_DESCRIPTION)
  (int& tag, int& cTag, int& meshControl,
   int& type, int& nodeC, int& info)
{
  //NOTE: nodes argument is NULL <--> no node tags returned!
  if(geometryAgent->nextElement(tag, cTag, meshControl, type, nodeC, NULL) != -1)
    info = 0;
  else
    info = -1;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_geometry_node,EIO_SET_GEOMETRY_NODE)
  (int& tag, int& cTag, double *coord, int& info)
{
  geometryAgent->writeNode(tag, cTag, coord);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_geometry_node,EIO_GET_GEOMETRY_NODE)
  (int& tag, int& cTag, double *coord, int& info)
{
  if(geometryAgent->nextNode(tag, cTag, coord) != -1)
    info = 0;
  else
    info = -1;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_geometry_boundary,EIO_SET_GEOMETRY_BOUNDARY)
  (int& tag, int& left, int& right, int& info)
{
  geometryAgent->writeBoundary(tag, left, right);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_geometry_boundary,EIO_GET_GEOMETRY_BOUNDARY)
  (int& tag, int& left, int& right, int& info)
{
  if(geometryAgent->nextBoundary(tag, left, right) != -1)
    info = 0;
  else
    info = -1;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_create_mesh,EIO_CREATE_MESH)
  (const char *directory, int& info)
{
  if(meshAgent = new EIOMeshAgent(modelManager))
    {
      info = meshAgent->createMesh(directory);
    }
  else
    {
      info = -1;
    }  
}

extern "C" void STDCALLBULL FC_FUNC_(eio_open_mesh,EIO_OPEN_MESH)
  (const char *directory, int& info)
{
  if(paraState.isParallel)
    meshAgent = new EIOMeshAgent(modelManager, paraState.numProc,
				 paraState.myProc);
  else
    meshAgent = new EIOMeshAgent(modelManager);

  if(meshAgent)
    {
      info = meshAgent->openMesh(directory);
    }
  else
    {
      info = -1;
    }
}

extern "C" void STDCALLBULL FC_FUNC_(eio_close_mesh,EIO_CLOSE_MESH)
  (int& info)
{
  meshAgent->closeMesh();
  delete meshAgent;
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_mesh_description,EIO_SET_MESH_DESCRIPTION)
  (int& nodeCount, int& elementCount, 
   int& boundaryElementCount, 
   int& usedElementTypes, int* elementTypeTags,
   int *elementCountByType,
   int& info)
{
  meshAgent->write_descriptor(nodeCount, elementCount, boundaryElementCount,
			      usedElementTypes, elementTypeTags,
			      elementCountByType);    
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_mesh_node,EIO_SET_MESH_NODE)
  (int &tag, int& type, double *coord, int& info)
{
  meshAgent->write_node(tag, type, coord);
  info = 0;  
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_mesh_element_conns,EIO_SET_MESH_ELEMENT_CONNS)
  (int& tag, int& body, 
   int& type, int *nodes, 
   int& info)
{
  meshAgent->write_elementConnections(tag, body, type, nodes);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_mesh_bndry_element,EIO_SET_MESH_BNDRY_ELEMENT)
  (int& tag, int& boundary,
   int& leftElement, int& rightElement,
   int& type, int* nodes,
   int& info)
{
  meshAgent->write_boundaryElement(tag, boundary,
				   leftElement, rightElement, 
				   type,
				   nodes);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_mesh_bndry_element,EIO_GET_MESH_BNDRY_ELEMENT)
  (int& tag, int& part, int& boundary, 
   int& leftElement, int& rightElement,
   int& type, int* nodes, double *coord,
   int& info)
{
  if(meshAgent->read_nextBoundaryElement(tag, part, boundary, 
					 leftElement, rightElement, 
					 type, nodes, coord) != -1)
    info = 0;
  else
    info = -1;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_mesh_description,EIO_GET_MESH_DESCRIPTION)
  (int& nodeCount, int& elementCount, 
   int& boundaryElementCount, 
   int& usedElementTypes, int* elementTypeTags,
   int *elementCountByType, int& info)
{
  meshAgent->read_descriptor(nodeCount, elementCount, boundaryElementCount, 
			     usedElementTypes, elementTypeTags, 
			     elementCountByType);
  info = 0;
}
extern "C" void STDCALLBULL FC_FUNC_(eio_get_mesh_element_conns,EIO_GET_MESH_ELEMENT_CONNS)
  (int& tag, int& part, int& body, int& type, int *pdofs, int *nodes, 
   int& info)
{
  if(meshAgent->read_nextElementConnections(tag, part, body, type, pdofs, nodes) != -1)
    info = 0;
  else info = -1;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_mesh_element_coords,EIO_GET_MESH_ELEMENT_COORDS)
  (int& tag, int& body, int& type, int *nodes, 
   double *coord, int& info)
{
  if(meshAgent->read_nextElementCoordinates(tag, body, type, nodes, coord) != -1)
    info = 0;
  else info = -1;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_mesh_nodes,EIO_GET_MESH_NODES)
  (int *tags,double *coord, int& info)
{
  meshAgent->read_allNodes(tags,coord);
  info = 0;
}


extern "C" void STDCALLBULL FC_FUNC_(eio_create_dual_mesh,EIO_CREATE_DUAL_MESH)
  (const char *dir, int& info)
{
  if(dualMeshAgent = new EIODualMeshAgent(modelManager))
    {
      info = dualMeshAgent->createMesh(dir);
    }
  else
    {
      info = -1;
    }    
}

extern "C" void STDCALLBULL FC_FUNC_(eio_open_dual_mesh,EIO_OPEN_DUAL_MESH)
  (const char *dir, int& info)
{
  if(dualMeshAgent = new EIODualMeshAgent(modelManager))
    {
      info = dualMeshAgent->openMesh(dir);
    }
  else
    {
      info = -1;
    }  
}

extern "C" void STDCALLBULL FC_FUNC_(eio_close_dual_mesh,EIO_CLOSE_DUAL_MESH)
  (int& info)
{
  dualMeshAgent->closeMesh();
  delete dualMeshAgent;
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_dual_mesh_element_conns,EIO_SET_DUAL_MESH_ELEMENT_CONNS)
  (int& tag, int& type, int *nodes, int& info)
{
  if(dualMeshAgent->write_elementConnections(tag, type, nodes) != -1)
    info = 0;
  else info = -1;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_dual_mesh_element_conns,EIO_GET_DUAL_MESH_ELEMENT_CONNS)
  (int& tag, int& type, int *nodes, int& info)
{
  if(dualMeshAgent->read_nextElementConnections(tag, type, nodes) != -1)
    info = 0;
  else info = -1;
}


extern "C" void STDCALLBULL FC_FUNC_(eio_create_part,EIO_CREATE_PART)
  (const char *dir, int& parts, int& info)
{
  if(partitioningWriter = new EIOPartWriter(parts, modelManager))
    {
      info = partitioningWriter->createPartitioning(dir);
    }
  else
    {
      info = -1;
    }   
}

extern "C" void STDCALLBULL FC_FUNC_(eio_close_part,EIO_CLOSE_PART)
  (int& info)
{
  if(paraState.isParallel == 0)
    {
      partitioningWriter->closePartitioning();
      delete partitioningWriter;
    }
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_part_description,EIO_SET_PART_DESCRIPTION)
  (int& nodeCount, 
   int& sharedNodeCount,
   int& elementCount, 
   int& borderElementCount, 
   int& boundaryElementCount, 
   int& usedElementTypes, 
   int* elementTypeTags,
   int *elementCountByType,
   int& info)
{
  if(partitioningWriter->write_descriptor(nodeCount, 
					  sharedNodeCount,
					  elementCount, 
					  borderElementCount, 
					  boundaryElementCount, 
					  usedElementTypes, 
					  elementTypeTags,
					  elementCountByType) != -1)
    info = 0;
  else
    info = -1;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_part_description,EIO_GET_PART_DESCRIPTION)
  (int& sharedNodeCount,
   int& info)
{
  meshAgent->read_partDescriptor(sharedNodeCount);
  info = 0;
}


extern "C" void STDCALLBULL FC_FUNC_(eio_activate_part_part,EIO_ACTIVATE_PART_PART)
  (int& part, int& info)
{
  if(partitioningWriter->activatePart(part) != -1)
    info = 0;
  else
    info = -1;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_deactivate_part_part,EIO_DEACTIVATE_PART_PART)
  (int& info)
{
  if(partitioningWriter->deactivatePart() != -1)
    info = 0;
  else
    info = -1;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_part_node,EIO_SET_PART_NODE)
  (int& tag, 
   int& type,      
   double *coord, 
   int& partcount, 
   int *parts,     
   int& info)
{
  if(partitioningWriter->write_node(tag, type, coord, partcount, parts) != -1)
    info = 0;
  else
    info = -1;  
}


extern "C" void STDCALLBULL FC_FUNC_(eio_get_part_node,EIO_GET_PART_NODE)
  (int& tag, 
   int& constraint,      
   double *coord, 
   int& partcount, 
   int *parts,     
   int& info)
{
  if(meshAgent->read_sharedNode(tag, 
				constraint, coord, partcount, parts) != -1)
    info = 0;
  else
    info = -1;  
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_part_element,EIO_SET_PART_ELEMENT)
  (int& tag, 
   int& body, 
   int& type, 
   int *nodes,
   int& border, 
   int& info)
{
  if(partitioningWriter->write_element(tag, body, type, 
				       nodes, border) != -1)
    info = 0;
  else
    info = -1;  
}


extern "C" void STDCALLBULL FC_FUNC_(eio_create_modeldata,EIO_CREATE_MODELDATA)
  (int& info)
{
  if(modelDataAgent = new EIOModelDataAgent(modelManager))
    {
      info = modelDataAgent->createModelData();
    }
  else
    {
      info = -1;
    }
}

extern "C" void STDCALLBULL FC_FUNC_(eio_open_modeldata,EIO_OPEN_MODELDATA)
  (int& info)
{
  if(modelDataAgent = new EIOModelDataAgent(modelManager))
    {
      info = modelDataAgent->openModelData();
    }
  else
    {
      info = -1;
    }
}

extern "C" void STDCALLBULL FC_FUNC_(eio_close_modeldata,EIO_CLOSE_MODELDATA)
  (int& info)
{
  modelDataAgent->closeModelData();
  delete modelDataAgent;
  info = 0;
}


extern "C" void STDCALLBULL FC_FUNC_(eio_set_modeldata_description,EIO_SET_MODELDATA_DESCRIPTION)
  (int& bodies,
   int& body_forces,
   int& body_equations,
   int& materials,
   int& boundary_conditions,
   int& initial_conditions,
   int& mesh_parameters,
   int& info)
{
  modelDataAgent->writeDescription(bodies,
				   body_forces,
				   body_equations,
				   materials,
				   boundary_conditions,
				   initial_conditions,
				   mesh_parameters);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_modeldata_description,EIO_GET_MODELDATA_DESCRIPTION)
  (int& bodies,
   int& body_forces,
   int& body_equations,
   int& materials,
   int& boundary_conditions,
   int& initial_conditions,
   int& mesh_parameters,
   int& info)
{
  modelDataAgent->readDescription(bodies,
				  body_forces,
				  body_equations,
				  materials,
				  boundary_conditions,
				  initial_conditions,
				  mesh_parameters);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_body,EIO_SET_BODY)
  (int& tag, int& body_force_id,
   int& equation_id, int& init_cond_id,
   int& material_id, int& mesh_param_id, int& info)
{
  modelDataAgent->writeBodyRecord(tag, body_force_id,
				  equation_id, init_cond_id,
				  material_id, mesh_param_id);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_body,EIO_GET_BODY)
  (int& tag, int& body_force_id,
   int& equation_id, int& init_cond_id,
   int& material_id, int& mesh_param_id, int& info)
{
  modelDataAgent->readBodyRecord(tag, body_force_id,
				 equation_id, init_cond_id,
				 material_id, mesh_param_id);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_constants,EIO_SET_CONSTANTS)
  (double* gravity, double& boltz, int& info)
{
  modelDataAgent->writeConstants(gravity, boltz);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_constants,EIO_GET_CONSTANTS)
  (double* gravity, double& boltz, int& info)
{
  modelDataAgent->readConstants(gravity, boltz);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_coords,EIO_SET_COORDS)
  (int& dim, int& coordsys, int *mapping,
   int& symmetry,
   double *start,
   double *end1, double* end2, int& info)
{
  modelDataAgent->writeCoordinates( dim,  coordsys, mapping,
				    symmetry,
				    start,
				    end1, end2);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_coords,EIO_GET_COORDS)
  (int& dim, int& coordsys, int *mapping,
   int& symmetry,
   double *start,
   double *end1, double* end2, int& info)
{
  modelDataAgent->readCoordinates(dim,  coordsys, mapping,
				  symmetry,
				  start,
				  end1, end2);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_material_head,EIO_SET_MATERIAL_HEAD)
  (int& tag, int& fields, int& info)
{
  modelDataAgent->writeMaterialHead(tag, fields);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_material_field,EIO_SET_MATERIAL_FIELD)
  (int& name,
   int& type,
   int& len,
   int *fields,
   double *values,
   int& info)
{
  modelDataAgent->writeMaterialField(name, type, len, fields, values);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_bndry_condition_head,EIO_SET_BNDRY_CONDITION_HEAD)
  (int& tag, int& fields, int& info)
{
  modelDataAgent->writeBoundaryConditionHead(tag, fields);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_bndry_condition_field,EIO_SET_BNDRY_CONDITION_FIELD)
  (int& name,
   int& type,
   int& len,
   int *fields,
   double *values,
   int& info)
{
  modelDataAgent->writeBoundaryConditionField(name, type, len, fields, values);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_material_head,EIO_GET_MATERIAL_HEAD)
  (int& tag, int& fields, int& info)
{
  modelDataAgent->readMaterialHead(tag, fields);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_material_field,EIO_GET_MATERIAL_FIELD)
  (int& name,
   int& type,
   int& len,
   int *fields,
   double *values,
   int& info)
{
  modelDataAgent->readMaterialField(name, type, len, fields, values);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_bndry_condition_head,EIO_GET_BNDRY_CONDITION_HEAD)
  (int& tag, int& fields, int& info)
{
  modelDataAgent->readBoundaryConditionHead(tag, fields);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_bndry_condition_field,EIO_GET_BNDRY_CONDITION_FIELD)
  (int& name,
   int& type,
   int& len,
   int *fields,
   double *values,
   int& info)
{
  modelDataAgent->readBoundaryConditionField(name, type, len, fields, values);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_initial_condition_head,EIO_SET_INITIAL_CONDITION_HEAD)
  (int& tag, int& fields, int& info)
{
  modelDataAgent->writeInitialConditionHead(tag, fields);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_initial_condition_field,EIO_SET_INITIAL_CONDITION_FIELD)
  (int& name,
   int& type,
   int& len,
   int *fields,
   double *values,
   int& info)
{
  modelDataAgent->writeInitialConditionField(name, type, len, fields, values);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_initial_condition_head,EIO_GET_INITIAL_CONDITION_HEAD)
  (int& tag, int& fields, int& info)
{
  modelDataAgent->readInitialConditionHead(tag, fields);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_initial_condition_field,EIO_GET_INITIAL_CONDITION_FIELD)
  (int& name,
   int& type,
   int& len,
   int *fields,
   double *values,
   int& info)
{
  modelDataAgent->readInitialConditionField(name, type, len, fields, values);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_body_equation_head,EIO_SET_BODY_EQUATION_HEAD)
  (int& tag, int& fields, int& info)
{
  modelDataAgent->writeBodyEquationHead(tag, fields);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_body_equation_field,EIO_SET_BODY_EQUATION_FIELD)
  (int& name,
   int& type,
   int& len,
   int *fields,
   double *values,
   int& info)
{
  modelDataAgent->writeBodyEquationField(name, type, len, fields, values);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_body_force_head,EIO_SET_BODY_FORCE_HEAD)
  (int& tag, int& fields, int& info)
{
  modelDataAgent->writeBodyForceHead(tag, fields);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_body_force_field,EIO_SET_BODY_FORCE_FIELD)
  (int& name,
   int& type,
   int& len,
   int *fields,
   double *values,
   int& info)

{
  modelDataAgent->writeBodyForceField(name, type, len, fields, values);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_mesh_parameter_head,EIO_SET_MESH_PARAMETER_HEAD)
  (int& tag, int& fields, int& info)
{
  modelDataAgent->writeMeshParameterHead(tag, fields);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_mesh_parameter_field,EIO_SET_MESH_PARAMETER_FIELD)
  (int& name,
   int& type,
   int& len,
   int *fields,
   double *values,
   int& info)
{
  modelDataAgent->writeMeshParameterField(name, type, len, fields, values);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_mesh_parameter_head,EIO_GET_MESH_PARAMETER_HEAD)
  (int& tag, int& fields, int& info)
{
  modelDataAgent->readMeshParameterHead(tag, fields);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_mesh_parameter_field,EIO_GET_MESH_PARAMETER_FIELD)
  (int& name,
   int& type,
   int& len,
   int *fields,
   double *values,
   int& info)
{
  modelDataAgent->readMeshParameterField(name, type, len, fields, values);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_body_equation_head,EIO_GET_BODY_EQUATION_HEAD)
  (int& tag, int& fields, int& info)
{
  modelDataAgent->readBodyEquationHead(tag, fields);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_body_equation_field,EIO_GET_BODY_EQUATION_FIELD)
  (int& name,
   int& type,
   int& len,
   int *fields,
   double *values,
   int& info)
{
  modelDataAgent->readBodyEquationField(name, type, len, fields, values);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_body_force_head,EIO_GET_BODY_FORCE_HEAD)
  (int& tag, int& fields, int& info)
{
  modelDataAgent->readBodyForceHead(tag, fields);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_body_force_field,EIO_GET_BODY_FORCE_FIELD)
  (int& name,
   int& type,
   int& len,
   int *fields,
   double *values,
   int& info)
{
  modelDataAgent->readBodyForceField(name, type, len, fields, values);
  info = 0;
}


extern "C" void STDCALLBULL FC_FUNC_(eio_create_solver,EIO_CREATE_SOLVER)
  (int& info)
{
  if(solverAgent = new EIOSolverAgent(modelManager))
    {
      info = solverAgent->createSolver();
    }
  else
    {
      info = -1;
    }
}

extern "C" void STDCALLBULL FC_FUNC_(eio_open_solver,EIO_OPEN_SOLVER)
  (int& info)
{
  if(solverAgent = new EIOSolverAgent(modelManager))
    {
      info = solverAgent->openSolver();
    }
  else
    {
      info = -1;
    }
}

extern "C" void STDCALLBULL FC_FUNC_(eio_close_solver,EIO_CLOSE_SOLVER)
  (int& info)
{
  solverAgent->closeSolver();
  delete solverAgent;
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_solver_description,EIO_SET_SOLVER_DESCRIPTION)
  (int& linsys, int& procs, int& info)
{
  solverAgent->writeDescription(linsys, procs);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_solver_description,EIO_GET_SOLVER_DESCRIPTION)
  (int& linsys, int& procs, int& info)
{
  solverAgent->readDescription(linsys, procs);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_solver,EIO_SET_SOLVER)
  (int& equation,
   int& main_type,
   int& sub_type,
   int& precond_type,
   int& stabilization,
   int& max_iter,
   double& stop_tol,
   double& steady_stop_tol,
   int& linearization,
   int& lin_max_iter,
   double& lin_stop_tol,
   int& lin_use_picard,
   int& lin_use_newton,
   int& newton_after_iter,
   double& newton_after_tol,
   int& info)
{
  solverAgent->writeSolverRecord(equation,
				 main_type,
				 sub_type,
				 precond_type,
				 stabilization,
				 max_iter,
				 stop_tol,
				 steady_stop_tol,
				 linearization,
				 lin_max_iter,
				 lin_stop_tol,
				 lin_use_picard,
				 lin_use_newton,
				 newton_after_iter,
				 newton_after_tol);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_solver,EIO_GET_SOLVER)
  (int& equation,
   int& main_type,
   int& sub_type,
   int& precond_type,
   int& stabilization,
   int& max_iter,
   double& stop_tol,
   double& steady_stop_tol,
   int& linearization,
   int& lin_max_iter,
   double& lin_stop_tol,
   int& lin_use_picard,
   int& lin_use_newton,
   int& newton_after_iter,
   double& newton_after_tol,
   int& info)
{
  solverAgent->readSolverRecord(equation,
				main_type,
				sub_type,
				precond_type,
				stabilization,
				max_iter,
				stop_tol,
				steady_stop_tol,
				linearization,
				lin_max_iter,
				lin_stop_tol,
				lin_use_picard,
				lin_use_newton,
				newton_after_iter,
				newton_after_tol);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_timestep_head,EIO_SET_TIMESTEP_HEAD)
  (int& dependence, int& len, int& info)
{
  solverAgent->writeTimestepDescription(dependence, len);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_timestep_head,EIO_GET_TIMESTEP_HEAD)
  (int& dependence, int& len, int& info)
{
  solverAgent->readTimestepDescription(dependence, len);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_set_timestep_field,EIO_SET_TIMESTEP_FIELD)
  (int& type,
   int *nof_timesteps,
   double *timestep_sizes,
   int *output_intervals,
   int& steady_max_iter,
   int& info)
{
  solverAgent->writeTimestepRecord(type,
				   nof_timesteps,
				   timestep_sizes,
				   output_intervals,
				   steady_max_iter);
  info = 0;
}

extern "C" void STDCALLBULL FC_FUNC_(eio_get_timestep_field,EIO_GET_TIMESTEP_FIELD)
  (int& type,
   int *nof_timesteps,
   double *timestep_sizes,
   int *output_intervals,
   int& steady_max_iter,
   int& info)
{
  solverAgent->readTimestepRecord(type,
				  nof_timesteps,
				  timestep_sizes,
				  output_intervals,
				  steady_max_iter);
  info = 0;
}



