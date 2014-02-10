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
 
struct parallel
{
  int isParallel;
  int numProc;
  int myProc;
};
#if 0
static struct parallel paraState = { 0, 0, 0};
static EIOModelManager *modelManager = 0;
static EIOGeometryAgent *geometryAgent = 0;
static EIOMeshAgent *meshAgent = 0;
static EIOSolverAgent *solverAgent = 0;
static EIOModelDataAgent *modelDataAgent = 0;
static EIODualMeshAgent *dualMeshAgent = 0;
static EIOPartWriter *partitioningWriter = 0;
#else
struct parallel paraState = { 0, 0, 0};
EIOModelManager *modelManager = 0;
EIOGeometryAgent *geometryAgent = 0;
EIOMeshAgent *meshAgent = 0;
EIOSolverAgent *solverAgent = 0;
EIOModelDataAgent *modelDataAgent = 0;
EIODualMeshAgent *dualMeshAgent = 0;
EIOPartWriter *partitioningWriter = 0;
#endif

#if defined(CBINDING)
#  define WITH_BINDING_EXTENSION(x) x
#else
#  if defined(FBINDING)
#     if !defined(AIX)
#        if defined (_UNICOS) | defined(WIN32)
#          define WITH_BINDING_EXTENSION(x) x
#        else
#          define WITH_BINDING_EXTENSION(x) x##_
#        endif
#     else
#        define WITH_BINDING_EXTENSION(x) x
#     endif
#     if defined(_UNICOS) | defined(WIN32)
#        define UCFUNNAME 1
#     endif
#     if defined(_UNICOS) | defined(WIN32)
#        define STRLENINT 1
#     endif
#  endif
#endif

#define EIOFC(funname) extern "C" void WITH_BINDING_EXTENSION(funname)

#if defined(UCFUNNAME)
EIOFC(EIO_INIT)
#else
  EIOFC(eio_init)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_INIT_PARALLEL)
#else
  EIOFC(eio_init_parallel)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_CLOSE)
#else
  EIOFC(eio_close)
#endif
  (int& info)
{
  delete modelManager;
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_CREATE_MODEL)
#else
  EIOFC(eio_create_model)
#endif
#if defined(STRLENINT)
  (const char *directory, int l1, int& info)
#else
  (const char *directory, int& info)
#endif
{
  info = modelManager->createModel(directory);
}

#if defined(UCFUNNAME)
EIOFC(EIO_OPEN_MODEL)
#else
  EIOFC(eio_open_model)
#endif
#if defined(STRLENINT)
  (const char *directory, int l1, int& info)
#else
  (const char *directory, int& info)
#endif
{
  info = modelManager->openModel(directory);
}

#if defined(UCFUNNAME)
EIOFC(EIO_CLOSE_MODEL)
#else
  EIOFC(eio_close_model)
#endif
  (int& info)
{
  info = modelManager->closeModel();
}

#if defined(UCFUNNAME)
EIOFC(EIO_CREATE_GEOMETRY)
#else
  EIOFC(eio_create_geometry)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_OPEN_GEOMETRY)
#else
  EIOFC(eio_open_geometry)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_CLOSE_GEOMETRY)
#else
  EIOFC(eio_close_geometry)
#endif
  (int& info)
{
  geometryAgent->closeGeometry();
  delete geometryAgent;
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_GEOMETRY_DESCRIPTION)
#else
  EIOFC(eio_set_geometry_description)
#endif
  (int& bodyC, int& boundaryC, int& outerC, 
   int& innerC, int& vertexC, 
   int& loopC, int& maxLooplen, int& info)
{
  geometryAgent->setDescriptor(bodyC, boundaryC, outerC, innerC, vertexC,
			       loopC, maxLooplen);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_GEOMETRY_DESCRIPTION)
#else
  EIOFC(eio_get_geometry_description)
#endif
  (int& bodyC, int& boundaryC, int& outerC, 
   int& innerC, int& vertexC, 
   int& loopC, int& maxLooplen, int& info)
{
  geometryAgent->descriptor(bodyC, boundaryC, outerC, innerC, vertexC,
			    loopC, maxLooplen);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_GEOMETRY_BODY)
#else
  EIOFC(eio_set_geometry_body)
#endif
  (int& tag, int& meshControl, int& loopC,
   int *loops,
   int& info)
{
  geometryAgent->writeBody(tag, meshControl, loopC, loops);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_GEOMETRY_BODY)
#else
  EIOFC(eio_get_geometry_body)
#endif
  (int& tag, int& meshControl, int& loopC,
   int *loops,
   int& info)
{
  if(geometryAgent->nextBody(tag, meshControl, loopC, loops) != -1)
    info = 0;
  else
    info = -1;
}

#if defined(UCFUNNAME)
EIOFC (EIO_SET_GEOMETRY_BODY_LOOP)
#else
  EIOFC (eio_set_geometry_body_loop)
#endif
  (int& tag, int& field, int *nodes, int& info)
{
  geometryAgent->writeLoop(tag, field, nodes);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC (EIO_GET_GEOMETRY_BODY_LOOP)
#else
  EIOFC (eio_get_geometry_body_loop)
#endif
  (int& tag, int& field, int *nodes, int& info)
{
  if(geometryAgent->nextLoop(tag, field, nodes) != -1)
    info = 0;
  else
    info = -1;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_GEOMETRY_ELEMENT)
#else
  EIOFC(eio_set_geometry_element)
#endif
  (int& tag, int& cTag, int& meshControl,
   int& type, int& nodeC, int *nodes, int& info)
{
  // Added nodeC argument: Martti Verho 1703.99
  geometryAgent->writeElement(tag, cTag, meshControl, type, nodeC, nodes);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_GEOMETRY_ELEMENT)
#else
  EIOFC(eio_get_geometry_element)
#endif
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
#if defined(UCFUNNAME)
EIOFC(EIO_GET_GEOMETRY_ELEMENT_DESCRIPTION)
#else
  EIOFC(eio_get_geometry_element_description)
#endif
  (int& tag, int& cTag, int& meshControl,
   int& type, int& nodeC, int& info)
{
  //NOTE: nodes argument is NULL <--> no node tags returned!
  if(geometryAgent->nextElement(tag, cTag, meshControl, type, nodeC, NULL) != -1)
    info = 0;
  else
    info = -1;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_GEOMETRY_NODE)
#else
  EIOFC(eio_set_geometry_node)
#endif
  (int& tag, int& cTag, double *coord, int& info)
{
  geometryAgent->writeNode(tag, cTag, coord);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_GEOMETRY_NODE)
#else
  EIOFC(eio_get_geometry_node)
#endif
  (int& tag, int& cTag, double *coord, int& info)
{
  if(geometryAgent->nextNode(tag, cTag, coord) != -1)
    info = 0;
  else
    info = -1;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_GEOMETRY_BOUNDARY)
#else
  EIOFC(eio_set_geometry_boundary)
#endif
  (int& tag, int& left, int& right, int& info)
{
  geometryAgent->writeBoundary(tag, left, right);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_GEOMETRY_BOUNDARY)
#else
  EIOFC(eio_get_geometry_boundary)
#endif
  (int& tag, int& left, int& right, int& info)
{
  if(geometryAgent->nextBoundary(tag, left, right) != -1)
    info = 0;
  else
    info = -1;
}

#if defined(UCFUNNAME)
EIOFC(EIO_CREATE_MESH)
#else
  EIOFC(eio_create_mesh)
#endif
  // Added: directory argument, Jouni Malinen 15-Oct-99
#if defined(STRLENINT)
  (const char *directory, int l1, int& info)
#else
  (const char *directory, int& info)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_OPEN_MESH)
#else
  EIOFC(eio_open_mesh)
#endif
  // Added: directory argument, Jouni Malinen 15-Oct-99
#if defined(STRLENINT)
  (const char *directory, int l1, int& info)
#else
  (const char *directory, int& info)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_CLOSE_MESH)
#else
  EIOFC(eio_close_mesh)
#endif
  (int& info)
{
  meshAgent->closeMesh();
  delete meshAgent;
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_MESH_DESCRIPTION)
#else
  EIOFC(eio_set_mesh_description)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_SET_MESH_NODE)
#else
  EIOFC(eio_set_mesh_node)
#endif
  (int &tag, int& type, double *coord, int& info)
{
  meshAgent->write_node(tag, type, coord);
  info = 0;  
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_MESH_ELEMENT_CONNS)
#else
  EIOFC(eio_set_mesh_element_conns)
#endif
  (int& tag, int& body, 
   int& type, int *nodes, 
   int& info)
{
  meshAgent->write_elementConnections(tag, body, type, nodes);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_MESH_BNDRY_ELEMENT)
#else
  EIOFC(eio_set_mesh_bndry_element)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_GET_MESH_BNDRY_ELEMENT)
#else
  EIOFC(eio_get_mesh_bndry_element)
#endif
  (int& tag, int& boundary, 
   int& leftElement, int& rightElement,
   int& type, int* nodes, double *coord,
   int& info)
{
  int part;
  if(meshAgent->read_nextBoundaryElement(tag, &part, boundary, 
					 leftElement, rightElement, 
					 type, nodes, coord) != -1)
    info = 0;
  else
    info = -1;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_MESH_DESCRIPTION)
#else
  EIOFC(eio_get_mesh_description)
#endif
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
#if defined(UCFUNNAME)
EIOFC(EIO_GET_MESH_ELEMENT_CONNS)
#else
  EIOFC(eio_get_mesh_element_conns)
#endif
  (int& tag, int& body, int& type, int *pdofs, int *nodes, 
   int& info)
{
  int part;
  if(meshAgent->read_nextElementConnections(tag, &part, body, type, pdofs, nodes) != -1)
    info = 0;
  else info = -1;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_MESH_ELEMENT_COORDS)
#else
  EIOFC(eio_get_mesh_element_coords)
#endif
  (int& tag, int& body, int& type, int *nodes, 
   double *coord, int& info)
{
  if(meshAgent->read_nextElementCoordinates(tag, body, type, nodes, coord) != -1)
    info = 0;
  else info = -1;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_MESH_NODES)
#else
  EIOFC(eio_get_mesh_nodes)
#endif
  (int *tags,double *coord, int& info)
{
  meshAgent->read_allNodes(tags,coord);
  info = 0;
}


#if defined(UCFUNNAME)
EIOFC(EIO_CREATE_DUAL_MESH)
#else
  EIOFC(eio_create_dual_mesh)
#endif
#if defined(STRLENINT)
  (const char *dir, int l1, int& info)
#else
  (const char *dir, int& info)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_OPEN_DUAL_MESH)
#else
  EIOFC(eio_open_dual_mesh)
#endif
#if defined(STRLENINT)
  (const char *dir, int l1, int& info)
#else
  (const char *dir, int& info)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_CLOSE_DUAL_MESH)
#else
  EIOFC(eio_close_dual_mesh)
#endif
  (int& info)
{
  dualMeshAgent->closeMesh();
  delete dualMeshAgent;
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_DUAL_MESH_ELEMENT_CONNS)
#else
  EIOFC(eio_set_dual_mesh_element_conns)
#endif
  (int& tag, int& type, int *nodes, int& info)
{
  if(dualMeshAgent->write_elementConnections(tag, type, nodes) != -1)
    info = 0;
  else info = -1;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_DUAL_MESH_ELEMENT_CONNS)
#else
  EIOFC(eio_get_dual_mesh_element_conns)
#endif
  (int& tag, int& type, int *nodes, int& info)
{
  if(dualMeshAgent->read_nextElementConnections(tag, type, nodes) != -1)
    info = 0;
  else info = -1;
}


#if defined(UCFUNNAME)
EIOFC(EIO_CREATE_PART)
#else
  EIOFC(eio_create_part)
#endif
#if defined(STRLENINT)
  (const char *dir, int l1, int& parts, int& info)
#else
  (const char *dir, int& parts, int& info)
#endif

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

#if defined(UCFUNNAME)
EIOFC(EIO_CLOSE_PART)
#else
  EIOFC(eio_close_part)
#endif
  (int& info)
{
  if(paraState.isParallel == 0)
    {
      partitioningWriter->closePartitioning();
      delete partitioningWriter;
    }
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_PART_DESCRIPTION)
#else
  EIOFC(eio_set_part_description)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_GET_PART_DESCRIPTION)
#else
  EIOFC(eio_get_part_description)
#endif
  (int& sharedNodeCount,
   int& info)
{
  meshAgent->read_partDescriptor(sharedNodeCount);
  info = 0;
}


#if defined(UCFUNNAME)
EIOFC(EIO_ACTIVATE_PART_PART)
#else
  EIOFC(eio_activate_part_part)
#endif
  (int& part, int& info)
{
  if(partitioningWriter->activatePart(part) != -1)
    info = 0;
  else
    info = -1;
}

#if defined(UCFUNNAME)
EIOFC(EIO_DEACTIVATE_PART_PART)
#else
  EIOFC(eio_deactivate_part_part)
#endif
  (int& info)
{
  if(partitioningWriter->deactivatePart() != -1)
    info = 0;
  else
    info = -1;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_PART_NODE)
#else
  EIOFC(eio_set_part_node)
#endif
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


#if defined(UCFUNNAME)
EIOFC(EIO_GET_PART_NODE)
#else
  EIOFC(eio_get_part_node)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_SET_PART_ELEMENT)
#else
  EIOFC(eio_set_part_element)
#endif
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


#if defined(UCFUNNAME)
EIOFC(EIO_CREATE_MODELDATA)
#else
  EIOFC(eio_create_modeldata)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_OPEN_MODELDATA)
#else
  EIOFC(eio_open_modeldata)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_CLOSE_MODELDATA)
#else
  EIOFC(eio_close_modeldata)
#endif
  (int& info)
{
  modelDataAgent->closeModelData();
  delete modelDataAgent;
  info = 0;
}


#if defined(UCFUNNAME)
EIOFC(EIO_SET_MODELDATA_DESCRIPTION)
#else
  EIOFC(eio_set_modeldata_description)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_GET_MODELDATA_DESCRIPTION)
#else
  EIOFC(eio_get_modeldata_description)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_SET_BODY)
#else
  EIOFC(eio_set_body)
#endif
  (int& tag, int& body_force_id,
   int& equation_id, int& init_cond_id,
   int& material_id, int& mesh_param_id, int& info)
{
  modelDataAgent->writeBodyRecord(tag, body_force_id,
				  equation_id, init_cond_id,
				  material_id, mesh_param_id);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_BODY)
#else
  EIOFC(eio_get_body)
#endif
  (int& tag, int& body_force_id,
   int& equation_id, int& init_cond_id,
   int& material_id, int& mesh_param_id, int& info)
{
  modelDataAgent->readBodyRecord(tag, body_force_id,
				 equation_id, init_cond_id,
				 material_id, mesh_param_id);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_CONSTANTS)
#else
  EIOFC(eio_set_constants)
#endif
  (double* gravity, double& boltz, int& info)
{
  modelDataAgent->writeConstants(gravity, boltz);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_CONSTANTS)
#else
  EIOFC(eio_get_constants)
#endif
  (double* gravity, double& boltz, int& info)
{
  modelDataAgent->readConstants(gravity, boltz);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_COORDS)
#else
  EIOFC(eio_set_coords)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_GET_COORDS)
#else
  EIOFC(eio_get_coords)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_SET_MATERIAL_HEAD)
#else
  EIOFC(eio_set_material_head)
#endif
  (int& tag, int& fields, int& info)
{
  modelDataAgent->writeMaterialHead(tag, fields);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_MATERIAL_FIELD)
#else
  EIOFC(eio_set_material_field)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_SET_BNDRY_CONDITION_HEAD)
#else
  EIOFC(eio_set_bndry_condition_head)
#endif
  (int& tag, int& fields, int& info)
{
  modelDataAgent->writeBoundaryConditionHead(tag, fields);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_BNDRY_CONDITION_FIELD)
#else
  EIOFC(eio_set_bndry_condition_field)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_GET_MATERIAL_HEAD)
#else
  EIOFC(eio_get_material_head)
#endif
  (int& tag, int& fields, int& info)
{
  modelDataAgent->readMaterialHead(tag, fields);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_MATERIAL_FIELD)
#else
  EIOFC(eio_get_material_field)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_GET_BNDRY_CONDITION_HEAD)
#else
  EIOFC(eio_get_bndry_condition_head)
#endif
  (int& tag, int& fields, int& info)
{
  modelDataAgent->readBoundaryConditionHead(tag, fields);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_BNDRY_CONDITION_FIELD)
#else
  EIOFC(eio_get_bndry_condition_field)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_SET_INITIAL_CONDITION_HEAD)
#else
  EIOFC(eio_set_initial_condition_head)
#endif
  (int& tag, int& fields, int& info)
{
  modelDataAgent->writeInitialConditionHead(tag, fields);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_INITIAL_CONDITION_FIELD)
#else
  EIOFC(eio_set_initial_condition_field)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_GET_INITIAL_CONDITION_HEAD)
#else
  EIOFC(eio_get_initial_condition_head)
#endif
  (int& tag, int& fields, int& info)
{
  modelDataAgent->readInitialConditionHead(tag, fields);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_INITIAL_CONDITION_FIELD)
#else
  EIOFC(eio_get_initial_condition_field)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_SET_BODY_EQUATION_HEAD)
#else
  EIOFC(eio_set_body_equation_head)
#endif
  (int& tag, int& fields, int& info)
{
  modelDataAgent->writeBodyEquationHead(tag, fields);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_BODY_EQUATION_FIELD)
#else
  EIOFC(eio_set_body_equation_field)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_SET_BODY_FORCE_HEAD)
#else
  EIOFC(eio_set_body_force_head)
#endif
  (int& tag, int& fields, int& info)
{
  modelDataAgent->writeBodyForceHead(tag, fields);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_BODY_FORCE_FIELD)
#else
  EIOFC(eio_set_body_force_field)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_SET_MESH_PARAMETER_HEAD)
#else
  EIOFC(eio_set_mesh_parameter_head)
#endif
  (int& tag, int& fields, int& info)
{
  modelDataAgent->writeMeshParameterHead(tag, fields);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_MESH_PARAMETER_FIELD)
#else
  EIOFC(eio_set_mesh_parameter_field)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_GET_MESH_PARAMETER_HEAD)
#else
  EIOFC(eio_get_mesh_parameter_head)
#endif
  (int& tag, int& fields, int& info)
{
  modelDataAgent->readMeshParameterHead(tag, fields);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_MESH_PARAMETER_FIELD)
#else
  EIOFC(eio_get_mesh_parameter_field)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_GET_BODY_EQUATION_HEAD)
#else
  EIOFC(eio_get_body_equation_head)
#endif
  (int& tag, int& fields, int& info)
{
  modelDataAgent->readBodyEquationHead(tag, fields);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_BODY_EQUATION_FIELD)
#else
  EIOFC(eio_get_body_equation_field)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_GET_BODY_FORCE_HEAD)
#else
  EIOFC(eio_get_body_force_head)
#endif
  (int& tag, int& fields, int& info)
{
  modelDataAgent->readBodyForceHead(tag, fields);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_BODY_FORCE_FIELD)
#else
  EIOFC(eio_get_body_force_field)
#endif
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


#if defined(UCFUNNAME)
EIOFC(EIO_CREATE_SOLVER)
#else
  EIOFC(eio_create_solver)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_OPEN_SOLVER)
#else
  EIOFC(eio_open_solver)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_CLOSE_SOLVER)
#else
  EIOFC(eio_close_solver)
#endif
  (int& info)
{
  solverAgent->closeSolver();
  delete solverAgent;
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_SOLVER_DESCRIPTION)
#else
  EIOFC(eio_set_solver_description)
#endif
  (int& linsys, int& procs, int& info)
{
  solverAgent->writeDescription(linsys, procs);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_SOLVER_DESCRIPTION)
#else
  EIOFC(eio_get_solver_description)
#endif
  (int& linsys, int& procs, int& info)
{
  solverAgent->readDescription(linsys, procs);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_SOLVER)
#else
  EIOFC(eio_set_solver)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_GET_SOLVER)
#else
  EIOFC(eio_get_solver)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_SET_TIMESTEP_HEAD)
#else
  EIOFC(eio_set_timestep_head)
#endif
  (int& dependence, int& len, int& info)
{
  solverAgent->writeTimestepDescription(dependence, len);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_GET_TIMESTEP_HEAD)
#else
  EIOFC(eio_get_timestep_head)
#endif
  (int& dependence, int& len, int& info)
{
  solverAgent->readTimestepDescription(dependence, len);
  info = 0;
}

#if defined(UCFUNNAME)
EIOFC(EIO_SET_TIMESTEP_FIELD)
#else
  EIOFC(eio_set_timestep_field)
#endif
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

#if defined(UCFUNNAME)
EIOFC(EIO_GET_TIMESTEP_FIELD)
#else
  EIOFC(eio_get_timestep_field)
#endif
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






