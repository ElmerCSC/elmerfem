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
Module:     ecif_model.h
Language:   C++
Date:       20.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   A Base class for model.
  Model stores general information concerning the CAD model.
  Datastructures storing bodies, elements and boundaries are
  also stored here.

************************************************************************/

#ifndef _ECIF_MODEL_
#define _ECIF_MODEL_

#include "ecif_def.h"
#include "ecif_def_trx.h"
#include "ecif_def_stl.h"
#include "ecif_model_aux.h"
#include "ecif_modelMeshManager.h"
#include "ecif_modelObject.h"
#include "ecif_modelOutputManager.h"

//enum ecifFieldInfo;

class Model
{
public:
  /*  friend class Control;
  friend class ModelMeshManager;
  friend class ModelOutputManager;
  friend class Renderer_OGL; */
  Model(char* name = NULL, ecif_modelSource = ECIF_CAD_OR_MESH_FILE, char* in_filename = NULL);
  ~Model();
  struct Units unit; // Muuta tämä suojatuksi !!
  //
  static void addColorName(char* name, Color4 color);
  static void addColorValue(int color_value, char* name);
  static void addMatcDefinition(char* def);
  static bool getColorValue(char* name, Color4 color);
  static bool getColorName(int id, char* buffer);
  static double getMeshInputUnit() { return meshInputUnit;};
  static int readColorFile(char* infile);
  static int readMatcFile(char* infile, char* info, bool must_exists);
  static int rgbColor2Id(Color4& color);
  static int rgbColor2Id(int len, char* hex_value);
  static void rgbColorId2Hex(int color_id, char* hex_value, int len);
  static void rgbHex2Color(int len , char* hex_value, Color4& color);
  static void setMeshInputUnit(double unit);
  static bool updateMatcFile(char* infile, char* mode, int def_count, char** defs);

  bool addBody(Body* body);
  bool addBody(ecif_Body_X& trx_body, bool add_default_layer = false);
  bool addBodyElement(BodyElement* be);
  bool addBodyElement(BodyElement* be, bool add_to_bodies);
  bool addBodyElement(BodyElement* be, int body1_tag, int body1_layer, int body2_tag, int body2_layer);
  bool addBodyElement(ecif_Element_X& trx_element);
  bool addBodyElement(ecif_Vertex_X& trx_vertex);
  bool addBodyElementGroup(BodyElementGroup* beg);
  bool addBodyElementLoop(BodyElementLoop* bel);
  bool addBodyElementLoop(ecif_ElementLoop_X& trx_element_loop);
  bool addBodyPair(BodyPair* bp);
  int addToElementGroup(BodyElement* be, int group_tag);
  void addCoveringElementList(int elem_id, IdList* se_list);
  Rc addMeshBoundaryElement(int bndr_tag, meshElementCode elem_code,
                            int paren1_tag, int parent2_tag, const int* node_ids,
                            bool& is_added_to_bndr);
  Rc addMeshBulkElement(int ext_id, int material_id,
                        meshElementCode elem_code, const int* node_ids,
                        int nof_ngbrs = 0, const int* ngbr_ids = NULL);
  Rc addMeshNode(int int_id, int ext_id, Point3& point);
  Rc addMeshInputElement(int elem_typr,
                         int ext_elem_id, int parent_tag, int ext_parent_tag,
                         int* node_ids);
  bool addModelObject(ModelObject* object, enum objectType otype);
  bool addParameter(ecif_parameterType parameter_type, Parameter* parameter);
  void addPendingMeshElements();
  void addSwapElement(int orig_elem_id, int swap_ele_id);
  bool addToHashPoints(GcPoint* point);
  void allocateMeshBodies(int nof_bodies) {meshManager->allocateMeshBodies(nof_bodies);}
  void allocateMeshBoundaryElements(int nof_elements) {meshManager->allocateMeshBoundaryElements(nof_elements);}
  void allocateMeshBulkElements(int nof_elements, int max_element_id) {meshManager->allocateMeshBulkElements(nof_elements, max_element_id);}
  void allocateMeshNodes(int nof_nodes, int max_node_id) {meshManager->allocateMeshNodes(nof_nodes, max_node_id);}
  void allocateMeshInputElements(int nof_elements, int max_ext_elem_id) {meshManager->allocateMeshInputElements(nof_elements, max_ext_elem_id);}
  void bodySelected(Renderer* renderer, int body_id, int lr_id);
  void boundarySelected(Renderer* renderer, int bndr_id,
                        int body1_id, int lr1_id, int body2_id, int lr2_id,
                        bool accept_body_change = false, bool update_gui = true);
	double calcInitialMeshH();
  void calcMeshBoundaryNodeNormals() { meshManager->calcMeshBoundaryNodeNormals(); }
  bool checkBodies();
  bool checkBodyElementLoops();
  bool checkBodyElements();
  bool checkBoundaries();
  void checkDiffuseGrayRadiation();
  bool checkMeshBodies();
  bool checkMeshElements(MeshElementTable* table, bool swap_to_ccw) {return meshManager->checkMeshElements(table, swap_to_ccw);}
  void checkMeshElementType(int elem_type, bool& is_bulk, bool& is_bndr, bool& is_edge) {meshManager->checkMeshElementType(elem_type, is_bulk, is_bndr, is_edge);}
  ecif_modelStatus checkStatus();
  void checkVertexExistence();
  bool checkVertexExistence(int tag);
  void checkMeshInputBoundaryElements() { meshManager->checkMeshInputBoundaryElements(); }
  void classifyMeshCornerElements() {meshManager->classifyMeshCornerElements();}
  void combineBoundaries(int body1_id, int body2_id);
  int convertElementCode(meshElementCode element_code) {return meshManager->convertElementCode(element_code);}
  meshElementCode convertElementType(int element_type) {return meshManager->convertElementType(element_type);}
  void convertMeshBulkElementIdsExt2Int() {meshManager->convertMeshBulkElementIdsExt2Int();}
  bool convertTags2Ids();
  void correctMeshZeroVelocityElements() {meshManager->correctMeshZeroVelocityElements();}
  BodyElement* createBodyElement(int body1_tag, int body1_layer, int body2_tag, int body2_layer, int nof_fem_elements);
  BodyElement* createBodyElement(int element_tag, int body1_tag, int body1_layer, int body2_tag, int body2_layer, int nof_fem_elements);
  void createMeshBodies() {meshManager->createMeshBodies();}
  void createMeshBodyTables() {meshManager->createMeshBodyTables();}
  void createMeshBoundaries(int nof_bulk_bndr_elems = 0, bool* free_bulk_bndr_elems = NULL) {meshManager->createMeshBoundaries(nof_bulk_bndr_elems, free_bulk_bndr_elems);}
  void createMeshBoundaryElementEdges();
  void createMeshBulkElementEdges();
  Parameter* createNewParameter(ecif_parameterType param_type,
                                int pid, int parent_id,
                                char* param_value, char* param_name);
  Parameter* createNewParameter(ecif_parameterType parameter_type, int id);
  void deleteParameters(ecif_parameterType param_type);
  int deleteParameter(ecif_parameterType param_type, int pid);
  void drawCurrentPickVector();
  void fieldNameGuiToSif(const char* gui_name, char* sif_name_buffer);
  void fieldNameSifToGui(const char* sif_name, char* gui_name_buffer);
  void findBoundaries();
  void findMeshBoundaryParents(int nof_bulk_bndr_elems, bool* free_bulk_bndr_flags) {meshManager->findMeshBoundaryParents(nof_bulk_bndr_elems, free_bulk_bndr_flags);}
  void findMeshElementNodeParents(MeshElementTable* source, int nof_nodes, int*& nofNodeParents, int**& nodeParentIds) {
                                      meshManager->findMeshElementNodeParents(source, nof_nodes, nofNodeParents, nodeParentIds);}
  void findNofBulkBoundaryElements(int& nof_bndr_elements) {meshManager->findNofBulkBoundaryElements(nof_bndr_elements);}
  void findMeshElementNeighbors(MeshElementTable* source) {meshManager->findMeshElementNeighbors(source);}
  int findSelectedMeshBoundaryElement(Renderer* renderer, Point3& ray_start, Point3& ray_dir,
                                      bool try_current_bndr, int& bndr_id,
                                      int& bd1_id, int& layer1_id,
                                      int& bd2_id, int& layer2_id);

  BodyElement* findVertex(GcPoint* point);
  const char* getActiveMeshName(int active_mesh_index);
  Body* getBody(int index, bool only_active = true);
  Body* getBodyById(int id) const;
  Body* getBodyByTag(int tag) const;
  BodyElement* getBodyElement(int index, bool only_active = true);
  BodyElement* getBodyElementById(int id) const;
  BodyElement* getBodyElementByBoundaryTag(int btag);
  BodyElement* getBodyElementByTag(objectType type, int tag);
  BodyElementGroup* getBodyElementGroup(int index, bool only_active = true);
  BodyElementGroup* getBodyElementGroupById(int id) const;
  BodyElementGroup* getBodyElementGroupByTag(int tag) const;
  BodyElementLoop* getBodyElementLoop(int index, bool only_active = true);
  BodyElementLoop* getBodyElementLoopById(int id) const;
  BodyElementLoop* getBodyElementLoopByTag(int tag) const;
  bool getBodyElementLoopId(IdList* loopd_ids, int& id, int& direction);
  bool getBodyElementLoopId(int nof_ids, const int* loopd_ids, int& id, int& direction);
  BodyLayer* getBodyLayer(int index, bool only_active = true);
  BodyLayer* getBodyLayerById(int id);
  BodyLayer* getBodyLayerByTag(int tag);
  int getBodyLayerByBodyId(int bd_id);
  BodyPair* getBodyPair(int index, bool only_active = true);
  BodyPair* getBodyPairById(const Body* bd1, const Body* bd2);
  int getBodyTagExt2Int(int external_tag);
  int getBodyTagInt2Ext(int internal_tag);
  BodyElement* getBoundary(int index, bool only_active = true);
  BodyElement* getBoundaryById(int id) const;
  BodyElement* getBoundaryByTag(int tag) const;
  BodyElement* getBoundaryByTags(int body1_tag, int body2_tag);
  void getBoundaries(int body1_id, int body2_id,
                     int& nof_boundaries,
                     BodyElement**& boundaries);
  void getBoundingBox(RangeVector rv) const;
  void getBoundingBox(double& x1, double& x2, double& y1, double& y2, double& z1, double&z2) const;
  const Control* getControlCenter() { return theControlCenter;}
  void getCoordinateLabels(int max_len, char* label_x, char* label_y, char* label_z);
  IdList* getCoveringElementList(int elem_id);
  bool getCoveringElementList(int master_elem_id, IdList& covering_elem_ids);
  void getCurrentBoundingBox(RangeVector rv);
  void getCurrentMeshBoundingBox(RangeVector rv);
  int getCurrentMeshIndex() { return modelInfo->currentMeshIndex; }
  const char* getCurrentMeshName();
  void getCurrentTime(char* buffer);
  ecif_modelDimension getDimension() const {return modelInfo->dimension;}

  BodyElement* getEdge(int index, bool only_active = true);
  BodyElement* getEdgeById(int id) const;
  BodyElement* getEdgeByTag(int tag) const;

  BodyElement* getFace(int index, bool only_active = true);
  BodyElement* getFaceById(int id) const;
  BodyElement* getFaceByTag(int tag) const;
  MeshCornerElement* getMeshCornerElement(int index);

  bool getFlagValue(flagName name);
  modelGeometryType getGeometryType() const {return modelInfo->geometryType;}
  const UserInterface* getGui();
  bool getLabelDisplayFlagValue(BodyElement* be);
  MeshElementTable* getMeshBoundaryElementEdges();
  MeshElementTable* getMeshBoundaryElements() {return meshData->boundaryElements;}
  MeshElementTable* getMeshBoundaryElementVertices() {return meshData->boundaryVertices;}
  void getMeshBoundingBox(RangeVector rv) const;
  MeshElementTable* getMeshBulkElements() {return meshData->bulkElements;}
  const char* getMeshDirValue() { return MESH_DIRECTORY_NAME; }
  MeshElementTable* getMeshBulkElementEdges() {return meshData->bulkEdges;}
  const MeshData* getMeshData() const { return meshData;}
  int getMeshBulkElementIdExt2Int(int ext_id) { return meshManager->getMeshBulkElementIdExt2Int(ext_id); }
  int getMeshBulkElementIdInt2Ext(int int_id) { return meshManager->getMeshBulkElementIdInt2Ext(int_id); }
  int getMeshInputElementIdExt2Int(int ext_id) { return meshManager->getMeshInputElementIdExt2Int(ext_id); }
  int getMeshInputElementType(int elem_id) { return meshManager->getMeshInputElementType(elem_id); }
  double getMeshF(int mesh_index);
  double getMeshH(int mesh_index);
  int getMeshIndex(const char* mesh_name);
  const MeshInfo* getMeshInfo() const { return meshInfo;}
  bool getMeshInputFileName(char*& mif_file_name, int mesh_index);
  int getMeshNodeIdExt2Int(int ext_id);
  int getMeshNodeIdInt2Ext(int int_id);
  void getMeshNames(int& nof_meshes, char**& mesh_names);
  Point3* getMeshNodeData() {return meshData->nodes;}

  const ModelInfo* getModelInfo() const { return modelInfo;}
  ecif_modelStatus getModelStatus() const {return modelInfo->modelStatus;}
  ostream& getModelStatusMessage(ostream& out) const;
  ModelObject* getModelObject(int index) const;
  ModelObject* getModelObject(int index, objectType type, bool only_active) const;
  ModelObject* getModelObjectById(int object_id) const;
  ModelObject* getModelObjectByTag(enum objectType type, int tag) const;
  const char* getModelObjectNameById(int oid);
  int getModelObjectTagById(int object_id) const;
  const ModelStatistics* getModelStatistics() const { return modelStatistics;}
  int getNewObjectId();

  int getNextNewParameterId(ecif_parameterType param_type);
  int getNofMeshes() { return modelInfo->nofMeshes; };
  int getNofMeshInputBoundaryElements() { return meshInfo->nofInputBoundaryElements; }
  int getNofTimestepSteps();
  Parameter* getParameter(int index, ecif_parameterType param_type);
  Parameter* getParameterById(ecif_parameterType parameter_type, int pid) ;
  ParameterFieldInfo* getParameterFieldInfo(const char* parameter, const char* field);
  GcPoint* getPoint(GcPoint* point);
  const ParallelInfo* getParallelInfo() const {return parallelInfo;}
  int getRelativeOrientation(BodyElement* be1, BodyElement* be2);
  BodyElement* getRemovedBodyElement(int id) const;
  int getRenumberedMeshBulkElementId(int original_id);
  flagName getSelectionMode();
  flagName getSelectionMethod();
  ecif_modelDimension getSimulationDimension() const {return modelInfo->simulationDimension;}
  bool getSolverKeywordTypeGiven(const char* parameter, const char* field);
  bool getSymmetryAxis(double start[3], double end1[3], double end2[3]);
  int getSwapElementId(int orig_elem_id);
  BodyElement* getVertex(GcPoint* point);
  BodyElement* getVertex(int index, bool only_active = true);
  BodyElement* getVertexById(int id);
  BodyElement* getVertexByNodeId(int node_id);
  BodyElement* getVertexByTag(int tag);
  const VertexTable* getVertexTable() { return modelData->vertexTable; }
  void initSplitCombineInfos();
  Rc installMeshInputBoundaryElements(bool clear_nodes = true) {return meshManager->installMeshInputBoundaryElements(clear_nodes);}
  Rc installMeshInputBulkElements(bool clear_nodes = true) {return meshManager->installMeshInputBulkElements(clear_nodes);}
  Rc installMeshInputElements(bool clear_nodes = true) {return meshManager->installMeshInputElements(clear_nodes);}
  bool isMeshInputBulkElement(int elem_id) {return meshManager->isMeshInputBulkElement(elem_id);}
  bool isInVertexTable(int vertex_id);
  bool keepMatcDefinitions() { return modelInfo->keepMatcDefinitions;}
  bool loadDBMesh(char* display_msg = NULL);
  bool loadMesh();
  void markActiveObjects();
  bool markObjectActive(int id);
  bool modelHasBulkRenumbering() { return meshData->bulkRenumbering != NULL; }
  bool modelHasCadGeometry();
  bool modelHasDiffuseGrayRadiation();
  bool modelHasEquation(const char* equation_name);
  bool modelHasMeshGeometry();
  bool modelHasParameter(ecif_parameterType parameter_type, const char* param_name);
  bool modelHasSteadyStateProblem();
  bool meshBoundaryElementSelected(Renderer* renderer, int fem_id);
  bool meshBoundaryElementSelectionHit(Renderer* renderer, int fem_id);
  bool meshBulkElementSelectionHit(Renderer* renderer, int fem_id);

  const char* objectType2Name(objectType type, int max_buf_len = 0, char* name_buffer = NULL);
  objectType objectName2Type(const char* name);

  ostream& outputSolverTargetFields_sif(ostream& out, short indent_size, short indent_level, const char* source_eq_name);

  bool processCadFileData();
  bool processMeshFileData(Input* input);
  bool processModelFileData();
  void processParametersAfterUpdate(ecif_parameterType parameter_type);
  void processParametersBeforeUpdate(ecif_parameterType parameter_type);

  void reallocateMeshBoundaryElements(int new_size) { meshManager->reallocateMeshBoundaryElements(new_size); }

  void refreshRenderer();

  int removeBody(Body* body);
  int removeBodyPair(const Body* body1, const Body* body2);
  int removeBodyElement(BodyElement* be, Body* body1, Body* body2, bool remove_subs);
  int removeBodyElement(BodyElement* be, bool remove_subs);
  int removeBodyElement(BodyElement* be, bool remove_from_bodies, bool remove_subs);
  int removeBodyElementGroup(int beg_id);
  int removeBodyElementLoop(int bel_id);
  int removeBodyLayer(BodyLayer* lr);
  void removeCadGeometry();
  void removeEmptyBoundaries();
  void removeMeshGeometry();
  void removeMeshInputElements() { meshManager->removeMeshInputElements(); }
  int removeVertex(BodyElement* vertex);

  void resetAllBoundarySelections(bool update_gui);
  void resetBoundaryConditions();
  void resetBoundarySelections(bool update_gui, bool use_boundary_groups, int nof_skip_ids = 0, const int* skip_bndr_ids = NULL, bool call_update = true);
  void resetInitialConditions();
  void resetMeshData();
  void resetMeshEdgesSelected() {meshManager->resetMeshEdgesSelected();}
  void resetMeshRendered() {meshManager->resetMeshRendered();}
  void resetMeshSelected() {meshManager->resetMeshSelected();}
  void resetModelData();

  int restoreBodyElement(BodyElement* be, bool remove_from_bodies);
  void restoreBoundaryNames();
  ostream& saveFrontModelFile(ostream& out, char* filename);
  ostream& saveSolverInputFile(ostream& out, char* filename);
  ostream& saveMeshInputFile(ostream& out, char* filename);
  void saveElmerMesh(char* mesh_dir);
  void saveElmerPostMesh(char* filename);
  void saveThetisMesh(char* filename);
  void saveUserSettingsFile(char* filename);
  void selectMeshBoundaryElement(int fem_id);
  void selectMeshBoundaryElements();
  void selectMeshBoundaryElementsUndo();
  void selectMeshBoundaryElementsRedo();
  void selectMeshBoundaryElementsAll();
  void selectMeshBoundaryElementsByNeighbor();
  void selectMeshBoundaryElementsByNormal();
  void selectMeshBoundaryElementsByPlane();
  void separateBodies();
  void setActiveMeshIndices(int nof_meshes,  int* mesh_indices);
  void setBoundaryConditions();
  void setCurrentAnchorPoint(int vertex_id);
  void setCurrentAnchorPoint(Point3 point);
  void setCurrentMeshIndex(int index);
  void setCurrentPickInfo(Point3 point, Point3 dir);
  void setInitialMeshH();
  void setFlagValue(flagGroup group, flagName name, bool value);
  void setKeepMatcDefinitions(bool value) {modelInfo->keepMatcDefinitions = value;}
  void setMatcInputFileEmf(char* matc_input_file);
  void setMatcInputFileSif(char* matc_input_file);
  void setMeshBodyExt2IntFlag(int external_id) { meshManager->setMeshBodyExt2IntFlag(external_id); }
  void setMeshBgMeshFiles(int nof_files, char** bg_mesh_file);
  void setMeshBgMeshFileIndices(int nof_files, int* bg_mesh_file_indices);
  void setMeshBgMeshActives(int nof_files, int* bg_mesh_actives);
  void setMeshBgMeshControls(int nof_files, int* bg_mesh_controls);
  void setMeshBulkElementParentId(int elem_id, int parent_id) { meshManager->setMeshBulkElementParentId(elem_id, parent_id); }
  void setMeshFs(int nof_values, double* mesh_fs);
  void setMeshHs(int nof_values, double* mesh_hs);
  void setMeshInputElementIsAdded(int elem_id, bool value) { meshManager->setMeshInputElementIsAdded(elem_id, value); }
  void setMeshInputElementParentTag(int elem_id, int parent_tag) { meshManager->setMeshInputElementParentTag(elem_id, parent_tag); }
  void setMeshInputElementExtParentTag(int elem_id, int ext_parent_tag) { meshManager->setMeshInputElementExtParentTag(elem_id, ext_parent_tag); }
  void setMeshNames(int nof_meshes, char** mesh_names);
  void setMeshNodes();
  void setModelCadSource(char* filename) {modelInfo->cadSourceFile = filename;}
  void setModelDescriptions(char* model_desc, char* problem_desc);
  void setModelDimension( ecif_modelDimension dim);
  void setModelFileCreated(char* created_str);
  void setModelFileDirectories(char* model_dir, char* include_path, char* results_dir, char* temp_dir);
  void setModelFileDirectoriesAbs(char* model_dir, char* include_path, char* results_dir, char* temp_dir);
  void setModelFileDirectoriesSave(bool include_path, bool results_dir, bool temp_dir);
  void setModelFileModified(char* modified_str);
  void setModelFileTime(char* time_str);
  void setModelHasUserDefinitions(bool value) { modelInfo->hasUserDefinitions = value;}
  void setModelHasDiffuseGrayRadiation(bool value);
  void setModelMeshDirectory(char* mesh_dir);
  void setModelMeshDirectoryAbs(char* mesh_dir);
  void setModelNameAndDirectory(char* file_path);
  void setModelNames(char* model_name, char* problem_name);
  bool setModelObjectNameById(int oid, char* name);
  bool setModelObjectTagById(int oid, int tag);
  void setModelStatus(ecif_modelStatus status) {modelInfo->modelStatus = status;}
  void setNofMeshInputBoundaryElements(int count) { meshInfo->nofInputBoundaryElements = count; }
  void setReadingModelFile(bool value) {modelInfo->readingModelFile = value;}
  void setParameter(ecif_parameterType param_type, int pid, int parent_id,
                    char* param_value, char* param_name);
  void setParallelInfo(ParallelInfo& pi);
  void setSelectionsToGui();
  void setTimestamp(char* target_name, char* time_str);
  void setVertexTable(int dim1, int dim2, int* vertex_ids, MatcValueTable& matc_table);
  void setWindowTitles();
  void sortBoundaryVertices(enum modelGeometryType gtype);
  void splitBoundary(int body1_id, int body2_id);
  void splitCombineBoundariesRedo();
  void splitCombineBoundariesUndo();
  void startEditMeshBoundaries();
  void stopEditMeshBoundaries(bool cancel_edit);
  void storeBoundaryNames();
  void swapBodyElements(IdArray* ids1, IdArray* ids2, IdArray* relative_dirs);
  void unloadMesh(char* msg = NULL);
  void unselectMeshBoundaryElement(int fem_id);
  void updateBodyForceApplyCounts();
  void updateBoundaryConditionApplyCounts();
  void updateBoundaries();
  void updateEquationApplyCounts();
  void updateCadGeometry();
  void updateInitialConditionApplyCounts();
  void updateMaterialApplyCounts();
  void updateMeshDirectoryInfo();
  void updateMinimumEdgeSize(int nof_points, GcPoint** boundary_points);
  void updateModelDirectoryInfo();
  void updateModelNameInfo();
  void updateModelStatistics();
  void updateParametersApplyCounts(ecif_parameterType parameter_type);
  void variableNameGuiToSif(const char* gui_name, char* sif_name_buffer);
  void variableNameSifToGui(const char* sif_name, char* gui_name_buffer);

  /* protected: */
  static Control* theControlCenter;
  static NameTable* colorNameTable;
  static RGBColorTable* colorValueTable;
  static double meshInputUnit; // Mesh scaling unit (0.001 <--> mm etc)

  Body* GBODY;
  MeshData* meshData;
  MeshInfo* meshInfo;
  ModelMeshManager* meshManager;
  ModelOutputManager* outputManager;
  ModelData* modelData;
  bool* modelFlags;
  ModelInfo* modelInfo;
  ModelStatistics* modelStatistics;
  ParallelInfo* parallelInfo;
  PickInfo* pickInfo;
  int id;
  int lastObjectId;
  BoundBox* modelBox;
  BoundBox* meshBox;

  int addVertex(BodyElement* vertex); // Adds new vertex2point entry
  bool checkElementGroupData();
  // Delete Cad /Mesh data (geometry)
  void deleteCadData();
  void deleteMeshData();

  void findContainedBodies();
  BodyPairArray* findNeighbourCandidates();
  void findBodyPairs();
  AdjacentPairArray* findBoundaryCandidates();
  void findInnerBoundaries();
  GcPoint* findPoint(PointList* plist, GcPoint* point);
  void getBodyTags(int bd1_id, int bd2_id, int& bd1_tag, int& bd2_tag);
  void getMeshRangeVector(RangeVector rv) const;
  void getRangeVector(RangeVector rv) const;
  void initBoundaries();
  void initMeshData();
  void initModelFlags();
  void normalizeMeshPoints();
  void normalizeVertices();
  int removeModelObject(int object_id);
  int restoreModelObject(ModelObject* obj);
  bool selectParameterTable(ecif_parameterType param_type,
                            ParameterTable*& table,
                            int*& use_counter);
  void setBodyColors(ColorIndexArray& unused_colors);
  void setBoundaryParentIdsAndTags();
  void setBoundaryPointData();
  void setBoundaryTags();
  void setMeshFs();
  void setMeshHs();
  int setNewSplitCombineIndex(short shift);
  void sortMeshBoundaryIndices();
  void updateParametersParentId();  // Update emf-file info (boundary groups etc.)
  void updateParametersParentId(ecif_parameterType param_type);
  void updateParametersParentTags();  // For pre version-5 files

};



#endif
