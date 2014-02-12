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
Module:     ecif_model_aux.h
Language:   C++
Date:       20.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Model helper structures

************************************************************************/

#ifndef _ECIF_MODEL_AUX_
#define _ECIF_MODEL_AUX_

#include "ecif_def.h"
#include "ecif_def_stl.h"

class MeshElementTable;
class MeshBulkElementTable;
class MeshFaceElementTable;
class MeshEdgeElementTable;


// This class stores info on mesh element connections to parent
// elements, like edge elements to boundaries, node to edges etc.
// If connections are numerated, each new connection stores also
// numerating number for the connections
//
class MeshConnectionTable {
public:
  MeshConnectionTable(int count, int alloc_size, int realloc_size,
                      bool enumerate_connections = false);
  ~MeshConnectionTable();
  bool addConnection(int index, int parent_id, bool only_unique = false);
  bool addConnection(int index, int parent_id, int& connection_nbr, bool only_unique = false);
  bool addConnections(int nof_indices, int* indices, int parent_id, bool only_unique = false);
  int getConnectionNbr(int index, int parent_index);
  const int* getConnectionNbrs(int index);
  int getNofConnections(int index);
  int getParentId(int index, int parent_index);
  const int* getParentIds(int index);
  bool getParentIds(int nof_indices, int* indices, int max_buffer_size, int* parent_id_buffer,
                 int minimum_connection_count, int& nof_ids);
  void sortParentIds();
protected:
  bool enumerateConnections;
  int lastConnectionNbr;
  int* nofConnections;    // Current nof connections
  int* maxNofConnections; // Max nof of connection
  int** parentIds;        // List of ids to connected parents
  int** connectionNbrs;   // List of enumerating numbers for the connections
  static int allocSize;   // Initial connection allocation size
  static int reallocSize; // Connection reallocation size
  static int totalCount;  // Total nof items
};


// This structure is used to store info for mesh elements during
// mesh input file reading
//
struct MeshInputElement {
  MeshInputElement();
  ~MeshInputElement();
  char type; // U=bulk, B=bndr, E=edge, V=vertex
  meshElementCode elementCode;  // Elmer elemen code (MEC_303 etc.)
  int elementId;    // Internal element id (both types have own numbering!)
  int extElementId; // External element id
  int extParentTag; // External body/boundary tag (material-id)
  int* extNodeIds;  // External node ids
  int parentTag;    // Internal body/boundary tag (counter)
  char isAdded;     // 0/1
};


// Used  to store info on mesh corner elements, ie elements
// whose n-1 faces are at boundary
struct MeshCornerElement {
  MeshCornerElement();
  ~MeshCornerElement();
  int elementId;              // Original element id
  int elementCode;            // Original element code
  int bodyId;                 // Element bodyi id
  int* nodeIds;               // Current node ids
  bool hasZeroVelocity;       // Zero-velocity bc flag
  bool corrected;             // True if corrected for zero-velocity
  bool splitted;              // True if corrected by splitting
  int nofSubElements;         // Original element nof sub elems
  int* subElementIds;         // Original element subelem ids
  int* boundaryElementIds;    // Original element subelem ids in boundarElements table
  int* boundaryIds;           // Original element subelem bounary ids
  int newNodeId;              // New (center) node id for splitted elemtent
  int nofNewElements;         // Nof of new elements for splitted element
  int* newElementIds;         // New elment ids for splitted element
  Point3 centerPoint;         // Element center point (location of the new node!)
};


struct MeshData {
  MeshData();
  ~MeshData();
  void resizeFlagTable(bool*& id_table, int old_size, int new_size, bool* copy_flags);
  void resizeIdTable(int*& id_table, int old_size, int new_size, bool* copy_flags);
  void resizeNodeTable(int old_size, int new_size, bool* copy_flags);

  MeshElementTable* bulkElements;
  MeshElementTable* bulkEdges;
  bool* bulkEdgeRendered;
  int* bulkRenumbering; // This is needed if splitted ("holes") elements

  MeshElementTable* boundaryElements;
  MeshElementTable* boundaryEdges;
  MeshElementTable* boundaryVertices;
  int* boundaryElementBoundaryIds;

  int* bulkElementExt2Int;
  int* bulkElementInt2Ext;
  int* bodyExt2Int;
  int* bodyInt2Ext;
  int nofAllocNodes;
  int nofNodes;
  int* nofBoundaryNodeParents;  // Nof boundary element parents for each boundary node (NOTE: a sparse table, size is nofNodes)
  int** boundaryNodeParentIds;  // Boundary element parent ids for each boundary node (NOTE: a sparse table, size is nofNodes)
  int lastBulkElementIndex;
  int lastBoundaryElementIndex;
  Point3* nodes;
  int* nodeExt2Int;
  int* nodeInt2Ext;
};


struct MeshInfo {
  MeshInfo();
  ~MeshInfo();
  int selectedBndrElementId;
  int selectedBulkElementId;
  int nofCornerElements;
  int nofZeroVelocityElements;
  int nofSplittedElements;
  int nofBulkEdges;
  int nofBulkElements;
  int nofVertices;
  int nofAllocElements;
  int maxExternalElementId;
  int nofNodes;
  int nofAllocNodes;
  int maxExternalNodeId;
  int nofBodies;
  //int nofBoundaries;
  int nofBoundaryEdges;
  int nofBoundaryElements;
  int nofBoundaryVertices;
  int nofInputBoundaryElements;
  //int nofInnerBoundaries;
  int nofInnerBndrElements;
  //int nofOuterBoundaries;
  int nofOuterBndrElements;
  int* nofUsedElementTypes;
  double minX, maxX, minY, maxY, minZ, maxZ;
  double dimX, dimY, dimZ, dimAvg;
};


// This is used ex. to store vertices to be printed as a vertex table
// in emf-file
//
struct VertexTable {
  VertexTable();
  ~VertexTable();

  int dim1;
  int dim2;
  int* vertexIds;
  MatcValueTable matcTable;
};


struct ModelData {
  ModelData();
  ~ModelData();

  void purgeBodyElementTable(BodyElementTable*& table);
  void purgeMeshCornerElements();
  void purgeCreatedElements();
  void purgeModelObjects();
  void purgeModelPoints();
  void purgeParameterTable(ParameterTable*& table);
  void purgeRemovedElements();
  void setVertexTable(int dim1, int dim2, int* vertex_ids, MatcValueTable& matc_table);
  void reset();

  AdjacentPairArray* boundaryCandidates;

  int lastPointId;

  MeshCornerElementList* meshCornerElements;

  ModelObjectArray* modelObjects;

  VertexTable* vertexTable;

  PointHashTable* modelPoints;
  Point2VertexTable* modelPoint2Vertices;
  BodyPairArray* neighbourCandidates;

  BodyElementTable* createdModelElements;
  BodyElementTable* removedModelElements;

  SplitCombineInfoArray* splitCombineInfos;
  int splitCombineInfoIndex;

  IdListTable* coveringElementTable;
  IdNumberTable* swapElementTable;

  ParameterTable* bodyParameters;
  ParameterTable* boundaryParameters;
  ParameterTable* boundaryConditions;
  ParameterTable* bodyForces;
  ParameterTable* calculators;
  ParameterTable* constants;
  ParameterTable* coordinates;
  ParameterTable* datafiles;
  ParameterTable* equations;
  ParameterTable* equationVariables;
  ParameterTable* gridParameters;
  ParameterTable* gridHs;
  ParameterTable* initialConditions;
  ParameterTable* materials;
  ParameterTable* modelParameters;
  ParameterTable* simulationParameters;
  ParameterTable* solvers;
  ParameterTable* solverControls;
  ParameterTable* timesteps;
  ParameterTable* userSettings;
};


// A structure to collect general model info
struct ModelInfo {
  ModelInfo(char* model_name, ecif_modelSource source, char* in_file_name);
  ~ModelInfo();
  modelGeometryType updateGeometryType();
  void updateCoordinateMapping(const char* coordinate_mapping_str);
  modelCoordinateType updateCoordinateType(const char* coordinate_system);
  ecif_modelDimension updateSimulationDimension(const char* coordinate_system);

  // NOTE: We want to track one level back with the input version number!
  int frontPreviousInputVersionNbr; // Front previous version number for the input model file (emf-file)
  int frontInputVersionNbr; // Front version number for the input model file (emf-file)
  int frontVersionNbr;      // Front version number for the output model file (emf-file)

  ecif_modelStatus modelStatus;
  enum ecif_modelSource modelSourceType;
  enum ecif_modelDimension dimension;
  enum ecif_modelDimension simulationDimension;
  enum modelGeometryType geometryType;
  enum modelCoordinateType coordinateType;
  int coordinateMapping[3];
  bool hasGeometry;

  int* activeMeshIndices;
  char* cadSourceFile;
  char* created;
  bool hasUserDefinitions;
  bool hasMatcDefinitions;
  char* includePath;
  char* includePath_absolute;
  bool  includePath_save;
  bool keepMatcDefinitions;
  char* matcInputFile_emf;
  char* matcInputFile_sif;
  int currentMeshIndex;
  int nofActiveMeshes;
  int nofMeshes;
  int nofBgMeshFiles;
  int* meshBgMeshFileIndices;
  char** meshBgMeshFiles;
  bool* meshBgMeshActives;
  bool* meshBgMeshControls;
  char** meshNames;
  double* meshFs;
  double* meshHs;
  char* meshDirectory;
  char* meshDirectory_absolute;
  char* meshSourceFile;
  char* meshResultFile;
  char* modelName;
  char* modelDescription;
  char* modelDirectory;
  char* modelDirectory_absolute;
  char* modelFileName;
  char* modelFileTs;
  char* modified;
  char* problemDescription;
  char* problemName;
  char* resultsDirectory;
  char* resultsDirectory_absolute;
  bool  resultsDirectory_save;
  char* temporaryFilesDirectory;
  char* temporaryFilesDirectory_absolute;
  bool  temporaryFilesDirectory_save;

  int selectedBodyElementId;
  int selectedBody1Id;
  int selectedBody2Id;
  int selectedLayer1Id;
  int selectedLayer2Id;
  bool editingMeshBoundaries;
  int nofEditableMeshBoundaries;
  BodyElement** editableMeshBoundaries;
  // Model (Cad) dimensions
  double minEdgeSize;
  double minX, maxX, minY, maxY, minZ, maxZ;
  double dimX, dimY, dimZ, dimAvg, dimMin, dimMax;

  // Timestamps when stored data was last change in the model file
  // Note! These are not (necesarily) times when the target was updated!
  // This time is stored in the target iself (Gebhard factors file etc.)
  char* databaseTs;
  char* frontTs;
  char* gebhardtFactorsTs;
  char* meshParameterTs;
  char* meshTs;
  char* solverTs;
  char* viewfactorsTs;

  // Flags if data was changed during the current session
  bool databaseNeedsUpdate;
  bool gebhardtFactorsNeedsUpdate;
  bool hasDiffuseGrayRadiation;
  bool meshNeedsUpdate;
  bool solverNeedsUpdate;
  bool viewfactorsNeedsUpdate;
  // Other flags
  bool readingModelFile;          // Model is reading ecf-model input file

};


struct ModelStatistics
{
  ModelStatistics();
  ~ModelStatistics();
  int nofBodies;
  int nofBodiesWithBodyForce;
  int nofBodiesWithEquation;
  int nofBodiesWithInitialCondition;
  int nofBodiesWithMaterial;
  int nofBodyForces;
  int nofBodyPairs;
  int nofBodyParameters;
  int nofBoundaryConditions;
  int nofBoundaryParameters;
  int nofBoundaryPoints;
  int nofCalculators;
  int nofConstants;
  int nofCoordinates;
  int nofDatafiles;
  int nofEdgesWithCondition;
  int nofElements;
  int nofElementLoops;
  int nofElementGroups;
  int nofEquations;
  int nofEquationVariables;
  int nofGridHs;
  int nofGridParameters;
  int nofInitialConditions;
  int nofInnerBoundaries;
  int nofInnerBoundariesWithCondition;
  int nofMaterials;
  int nofMeshBodies;
  int nofMeshBodyElements;
  int nofMeshBoundaryVertices;
  int nofModelObjects;
  int nofModelParameters;
  int nofOuterBoundaries;
  int nofOuterBoundariesWithCondition;
  int nofPatterns;
  int nofPoints;
  int nofSimulationParameters;
  int nofSolvers;
  int nofSolverControls;
  int nofTimesteps;  // NOTE: nof timestepping schemes!
  int nofUserSettings;
  int nofVertices;
  int nofVerticesWithCondition;
  int maxLoopSize;
};


struct ParallelInfo {
  ParallelInfo();
  ~ParallelInfo();
  int nofProcessors;
};


struct PickInfo {
  PickInfo() {
    anchorPoint[0] = anchorPoint[1] = anchorPoint[2] = NSVD;
    pickPoint[0] = pickPoint[1] = pickPoint[2] = NSVD;
    pickDir[0] = pickDir[1] = pickDir[2] = NSVD;
  }

  ~PickInfo() {};

  Point3 anchorPoint;
  Point3 pickPoint;
  Point3 pickDir;
};


#endif
