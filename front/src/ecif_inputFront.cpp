/*  
   ElmerFront - A graphical user interface of Elmer software
   Copyright (C) 1995- , CSC - IT Center for Science Ltd.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/***********************************************************************
Program:    ELMER Front
Module:     ecif_inputCadi.cpp
Language:   C++
Date:       20.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation (read Elmer Cadinterface model file)

************************************************************************/

#include "ecif_body.h"
#include "ecif_bodyLayer.h"
#include "ecif_bodyElement.h"
#include "ecif_bodyElementGroup.h"
#include "ecif_bodyElement1D.h"
#include "ecif_control.h"
#include "ecif_geometry.h"
#include "ecif_inputFront.h"
#include "ecif_model.h"
#include "ecif_model_aux.h"
#include "ecif_parameter.h"
#include "ecif_userinterface.h"
#include "ecif_timer.h"

extern char read_buffer[];

bool IS_FATAL = true;
bool NOT_FATAL = false;

const int isOk = emf_OK;
const int notOk = emf_ERROR;

// Init static class variables
//
bool InputFront::isChecked = false;
bool InputFront::isEgfInput = false;
bool InputFront::isEmfInput = false;
double InputFront::inputUnit = 1.0;  // By default input is in meters

template <class T> ostream&
output_array(ostream& out, UserInterface* gui, int length, T* array)
{
  ostrstream strm;
  for (int i = 0; i < length; i++) {
    strm << array[i];
  }
  strm << ends;
  gui->showMsg(strm.str());

  return out;
}


InputFront::InputFront(enum ecif_modelDimension m_dim,
                       ifstream& in_file, char* in_filename):
Input(m_dim, in_file, in_filename)
{
  isChecked = false;
  isEgfInput = false;
  isEmfInput = false;
  inputUnit = 1.0; // Default input unit is 1m
}


// Add normal (single) component
//
bool
InputFront::addElementComponent(ecif_Element_X& tx, ecif_geometryType gtype)
{
  // Current and new nof components
  int nof_components = tx.nof_components;
  tx.nof_components++;

  // Allocate new array and copy old components
  ecif_ElementComponent_X** components = new ecif_ElementComponent_X*[tx.nof_components];

  for (int i = 0; i < nof_components; i++) {
    components[i] = tx.components[i];
  }

  // Update tx
  delete[] tx.components;
  tx.components = components;

  // Allocate and store one new component
  ecif_ElementComponent_X* txc = new ecif_ElementComponent_X;
  init_trx_data(*txc);
  txc->gmtr_type = gtype;

  tx.components[tx.nof_components - 1] = txc;

  // Set geometry (union member) by element type
  if ( tx.tplg_type == ECIF_VERTEX) {
    txc->geometry.vertex = new ecif_VertexGeometry_X;
    init_trx_data(*txc->geometry.vertex);

  } else if ( tx.tplg_type == ECIF_EDGE) {
    txc->geometry.edge = new ecif_EdgeGeometry_X;
    init_trx_data(*txc->geometry.edge);

  } else if ( tx.tplg_type == ECIF_FACE ) {
    txc->geometry.face = new ecif_FaceGeometry_X;
    init_trx_data(*txc->geometry.face);
  }

  return true;
}


bool
InputFront::addFunctionComponentData(ecif_DllArg& da ,ecif_ElementComponent_X& tx)
{
  tx.isFunction = true;
  tx.isCpp = da.is_cpp;
  tx.isF95 = da.is_f95;
  tx.isMatc = da.is_matc;

  update_dyna_string(tx.functionName, da.func);
  update_dyna_string(tx.libraryName, da.lib);

  tx.startVertex = da.start_vertex;
  tx.endVertex = da.end_vertex;

  tx.argc = da.argc;

  if ( tx.argc > 0 ) {
    tx.argv = new double[tx.argc];
    for (int i = 0; i < tx.argc; i++) {
      tx.argv[i] = da.argv[i];
    }
  }

  if ( da.has_start_point ) {
    tx.startPoint = new Point3[1];
    (*tx.startPoint)[0] = da.start_point[0];
    (*tx.startPoint)[1] = da.start_point[1];
    (*tx.startPoint)[2] = da.start_point[2];
  }

  if ( da.has_end_point ) {
    tx.endPoint = new Point3[1];
    (*tx.endPoint)[0] = da.end_point[0];
    (*tx.endPoint)[1] = da.end_point[1];
    (*tx.endPoint)[2] = da.end_point[2];
  }

  // Copy any Matc-expression strings into component
  //
  copyMatcValueTable(da.matcTable, tx.matcTable);

  return true;
}


// Find geometry type code
//
ecif_geometryType
InputFront::getElementGeometryType(const char* type_name , ecif_Element_X& tx)
{
  ecif_geometryType gmtr_type;
  char* type = (char*) type_name;

  if ( LibFront::ncEqual(type, "Circle") ) {
    gmtr_type = ECIF_CIRCLE;
  }
  else if ( LibFront::ncEqual(type, "Ellipse") ) {
    gmtr_type = ECIF_ELLIPSE;
  }
  else if ( LibFront::ncEqual(type, "Hyperbola") ) {
    gmtr_type = ECIF_HYPERBOLA;
  }
  else if ( LibFront::ncEqual(type, "Line") ) {
    gmtr_type = ECIF_LINE;
  }
  else if ( LibFront::ncEqual(type, "Parabola") ) {
    gmtr_type = ECIF_PARABOLA;
  }
  else if ( LibFront::ncEqual(type, "PolyLine") ) {
    gmtr_type = ECIF_POLYLINE;
  }
  else if ( LibFront::ncEqual(type, "Spline") ) {
    gmtr_type = ECIF_SPLINE;
  }
  else if ( LibFront::ncEqual(type, "Nurbs") ) {
    gmtr_type = ECIF_NURBS;
  }
  else {
    gmtr_type = ECIF_NODIM;
  }

  return gmtr_type;
}


// Read one field for a body
int
InputFront::readBody(emf_ObjectData_X* od)
{
  static UserInterface* gui = theControlCenter->getGui();
  strstream strm1, strm2;

  static ecif_Body_X tx;
  static ecif_BodyLayer_X tx_layer;

  static int current_layer_tag = NO_INDEX;
  static char* current_layer_name = NULL;

  static IdArray edge_tags;
  static IdArray edge_group_tags;
  static IdArray vertex_groups;
  static IdArray vertex_tags;
  static int vertex_group_id = 0;
  static int layer = 0;
  static bool new_layer = false;

  static Body* body = NULL;
  static BodyLayer* bl = NULL;

  const char* fn = od->field_name;

  // This is the field id (when it make sense!)
  int id;
  LibFront::setNumericData(id, 0);

  bool isEgf = isEgfInput;
  bool isEmf = isEmfInput;

  // New body object started
  // =======================
  if (od->is_object_start) {

     // Init tx data
    init_trx_data(tx);
    init_trx_data(tx_layer);

    current_layer_tag = NO_INDEX;
    delete[] current_layer_name;
    current_layer_name = NULL;
    bl = NULL;

    tx.tag = od->object_id;
    tx.is_checked = geometryIsChecked();
  }


  // ---------------
  // Body level data
  // ---------------

  //-Body name
  if ( LibFront::in(EMF_NAME, fn) ) {
    //update_dyna_string(tx.name, (char*) od->data);
    readName(od, tx.name);
  }

  //-Body type
  else if ( LibFront::in(EMF_TYPE, fn) ) {
    if ( LibFront::in("Open", (char*) od->data) ) {
      tx.is_open = true;
    } else if ( LibFront::in("Bem", (char*) od->data) ) {
      tx.is_checked = true;
      tx.is_bem = true;
    } else if ( LibFront::in("Virtual", (char*) od->data) ) {
      tx.is_checked = true;
      tx.is_virtual = true;
    }
   }

  //-Body color
  else if ( LibFront::in(EMF_COLOR, fn) ) {
    if ( readColor(od, tx.color) ) {
      tx.has_color = true;
    }
  }

  //-Body Force id
  else if ( LibFront::in(EMF_BODY_FORCE, fn) ) {
    LibFront::setNumericData(tx.body_force_id, 0);
  }

  //-Body Parameter id
  else if ( LibFront::in(EMF_BODY_PARAMETER, fn) ) {
    LibFront::setNumericData(tx.body_param_id, 0);
  }

   //-Equation id
  else if ( LibFront::in(EMF_EQUATION, fn) ) {
    LibFront::setNumericData(tx.equation_id, 0);
  }

  //-Initial condition
  else if ( LibFront::in(EMF_INITIAL_CONDITION, fn) ) {
    LibFront::setNumericData(tx.init_cond_id, 0);
  }

  //-Material
  else if ( LibFront::in(EMF_MATERIAL, fn) ) {
    LibFront::setNumericData(tx.material_id, 0);
  }

  //-Edge groups (for virtual body)
  //
  else if ( LibFront::in("Edge Groups", fn) ||
            LibFront::in("Element Groups", fn)
          ) {

    int tag;
    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(tag, i);
      edge_group_tags.push_back(tag);
    }
  }


  // ----------------
  // Layer level data
  // ----------------

  //-Layer tag
  else if ( LibFront::in(EMF_LAYER, fn) ||
            LibFront::in(EMF_LAYER_TAG, fn)
          ) {
    LibFront::setNumericData(current_layer_tag, 0);
  }

  //-Layer type
  else if ( LibFront::in(EMF_LAYER_TYPE, fn) ) {
    if ( LibFront::in("Open", (char*) od->data) ) {
      tx_layer.is_open = true;
    }
  }

   //-Layer name
  else if ( LibFront::in(EMF_LAYER_NAME, fn) ) {
    tx_layer.has_name = true;
    update_dyna_string(current_layer_name, (char*) od->data);
  }

  //-Layer color
  else if ( LibFront::in(EMF_LAYER_COLOR, fn) ) {
    if ( readColor(od, tx_layer.color) ) {
      tx_layer.has_color = true;
    }
  }

  //-Element loop tags
  // EMF_EDGE_LOOPS (egf-style)
  // EMF_ELEMENT_LOOPS (emf-style)
  //
  else if ( LibFront::in(EMF_EDGE_LOOPS, fn) ||
            LibFront::in(EMF_ELEMENT_LOOPS, fn)
          ) {
    new_layer = true;

    tx_layer.nof_elem_loops = od->data_length;
    delete[] tx_layer.elem_loop_tags;
    tx_layer.elem_loop_tags = NULL;

    if ( tx_layer.nof_elem_loops > 0 ) {
      tx_layer.elem_loop_tags = new int[tx_layer.nof_elem_loops];
    }

    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(tx_layer.elem_loop_tags[i], i);
    }
  }

  //-Body vertices as tags (short egf-style)
  //
  else if ( LibFront::in(EMF_VERTICES, fn) ||
            LibFront::in(EMF_POLYGON, fn)
          ) {

    new_layer = true;

    vertex_group_id++;

    int tag, first_tag;

    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(tag, i);
      vertex_groups.push_back(vertex_group_id);
      vertex_tags.push_back(tag);

      if ( i == 0) {
        first_tag = tag;
      }
    }

    // Close the loop if given as open!
    if ( tag != first_tag ) {
      vertex_tags.push_back(first_tag);
      vertex_groups.push_back(vertex_group_id);
    }
  }

  //-Edge tags
  //
  else if ( LibFront::in(EMF_EDGES, fn) ) {

    new_layer = true;

    int tag;
    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(tag, i);
      edge_tags.push_back(tag);
    }
  }

  //-Grid parameter ids
  //
  else if ( LibFront::in(EMF_GRID_PARAMETER_IDS, fn) ) {
    tx_layer.nof_grid_param_ids = od->data_length;
    tx_layer.grid_param_ids = new int[od->data_length];

    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(tx_layer.grid_param_ids[i], i);
    }

    if ( bl != NULL ) {
      bl->setGridParameterIds(od->data_length,
                              tx_layer.grid_param_ids);
    }
  }

  //-Excluded mesh indices
  //
  else if ( LibFront::in(EMF_EXCLUDED_MESH_INDICES, fn) ) {
    tx_layer.nof_excluded_meshes = od->data_length;
    tx_layer.excluded_mesh_indices = new int[od->data_length];

    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(tx_layer.excluded_mesh_indices[i], i);
    }

    if ( bl != NULL ) {
      bl->setExcludedMeshData(od->data_length,
                              tx_layer.excluded_mesh_indices);
    }

  }

  //-Grid parameter mesh indices
  //
  else if ( LibFront::in(EMF_GRID_PARAMETER_MESH_INDICES, fn) ) {
    tx_layer.grid_param_mesh_indices = new int[od->data_length];

    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(tx_layer.grid_param_mesh_indices[i], i);
    }

    if ( bl != NULL ) {
      bl->setGridParameterMeshIndices(od->data_length,
                                      tx_layer.grid_param_mesh_indices);
    }

  }

  //-Unknown field
  //
  else {
    unknownFieldMsg(od, IS_FATAL);
  }

  //----Body and the first layer or some additional layer is ready
  //
  if ( new_layer ) {

    // First layer --> a new body
    if ( layer == 0 ) {
      model->addBody(tx);
      body = model->getBodyByTag(tx.tag);
    }

    // Add new layer to the body
    // =========================
    if ( body != NULL ) {

      tx_layer.tag = current_layer_tag;
      update_dyna_string(tx_layer.name, current_layer_name);

      tx_layer.body_id = body->Id();
      tx_layer.body_tag = body->Tag();

      body->addLayer(tx_layer);
      body->addElementLoopTags(layer, tx_layer.nof_elem_loops, tx_layer.elem_loop_tags);

      int nof_edges = edge_tags.size();
      int nof_vertices = vertex_tags.size();

      // Store edges as pending elements (will be resolved later when
      // the whole file is read
      if ( nof_edges > 0 ) {
        for (int i = 0; i < nof_edges; i++) {
          body->addPendingElement(layer, edge_tags[i]);
        }
        edge_tags.clear();
      }

      // Store vertices as pending start-end points (will be resolved later when
      // the whole file is read
      if ( nof_vertices > 0 ) {
        for (int i = 0; i < nof_vertices; i++) {
          body->addPendingVertex(layer, vertex_groups[i], vertex_tags[i]);
        }
        vertex_groups.clear();
        vertex_tags.clear();
        vertex_group_id = 0;
      }

      bl = body->getLayer(layer);

      new_layer = false;
      layer++;
      init_trx_data(tx_layer);
    }
  }

  // Reset for the next body
  //
  if ( od->is_object_end ) {

    int nof_edge_groups = edge_group_tags.size();

    if ( !tx.is_virtual && nof_edge_groups > 0 ) {
      strm1 << "***ERROR in definition for Body " << tx.tag << ends;
      strm2 << "Edge Groups are valid only for virtual bodies (Type=virtual)" << ends;
      gui->showMsg(strm1.str());
      gui->showMsg(strm2.str());

      goto error;
    }

    if ( !tx.is_virtual && body == NULL ) {
      strm1 << "***ERROR in definition for Body " << tx.tag << ends;
      strm2 << "No geometry defined for the body!" << ends;
      gui->showMsg(strm1.str());
      gui->showMsg(strm2.str());

      goto error;
    }

    if ( tx.is_virtual && body == NULL && nof_edge_groups == 0) {
      strm1 << "***ERROR in definition for Body " << tx.tag << ends;
      strm2 << "No edge groups defined for a virtual body!" << ends;
      gui->showMsg(strm1.str());
      gui->showMsg(strm2.str());

      goto error;
    }

    // If body not yet created (a virtual body)
    //
    if ( body == NULL ) {
      model->addBody(tx);
      body = model->getBodyByTag(tx.tag);
    }

    // Add possible pending group tags
    //
    if ( body != NULL && nof_edge_groups > 0 ) {
      for (int i = 0; i < nof_edge_groups; i++) {
        body->addElementGroupTag(edge_group_tags[i]);
      }
      edge_group_tags.clear();
    }

    body = NULL;
    bl = NULL;
    layer = 0;
    new_layer = false;

    reset_trx_data(tx);
  }

  return isOk;

  // Error
  error:
  reset_trx_data(tx);
  return notOk;
}


// Read boundary group
int
InputFront::readElementGroup(emf_ObjectData_X* od)
{
  UserInterface* gui = theControlCenter->getGui();
  const char* fn = od->field_name;

  static ecif_ElementGroup_X tx;
  static enum objectType group_type = OT_ELEMENT_GROUP;
  static enum objectType bndr_type = OT_NONE;

  bool isEgf = isEgfInput;
  bool isEmf = isEmfInput;

  //---New loop is started
  if (od->is_object_start) {

    // Init tx data
    init_trx_data(tx);
    tx.tag = od->object_id;

  }

  // Parse fields
  // ============

  //-Element ids
  if ( LibFront::in(EMF_ELEMENTS, fn) ||
       LibFront::in("Faces", fn) ||
       LibFront::in("Edges", fn) ||
       LibFront::in("Vertices", fn)
     ) {
    int size = od->data_length;
    tx.element_tags = new int[size];
    tx.nof_elements = size;
    for (int i = 0; i < size; i++)
      LibFront::setNumericData(tx.element_tags[i], i);

    // Set boundary group type
    if ( LibFront::in("Faces", fn) ) {
      bndr_type = OT_FACE;
    } else if ( LibFront::in("Edges", fn) ) {
      bndr_type = OT_EDGE;
    } else if ( LibFront::in("Vertices", fn) ) {
      bndr_type = OT_VERTEX;
    }
  }

  //-Name
  else if ( LibFront::in(EMF_NAME, fn) ) {
    //update_dyna_string(tx.name, (char*) od->data);
    readName(od, tx.name);
  }

  //-Group type
  else if ( LibFront::in(EMF_TYPE, fn) ) {
    if ( LibFront::in("Virtual", (char*) od->data) ) {
      tx.is_virtual = true;
    }
  }

  //-Boundary Condition id
  else if ( LibFront::in(EMF_BOUNDARY_CONDITION, fn) ) {
    LibFront::setNumericData(tx.boundary_cond_id, 0);
  }

  //-Boundary Parameter id
  else if ( LibFront::in(EMF_BOUNDARY_PARAMETER, fn) ) {
    LibFront::setNumericData(tx.boundary_param_id, 0);
  }

  //-Unknown field
  else {
    return unknownFieldMsg(od, IS_FATAL);
  }

  //----Group is ready, create it (as a model object)
  if (od->is_object_end) {
    BodyElementGroup* bg = new BodyElementGroup(tx, bndr_type);
    reset_trx_data(tx);
  }

  return isOk;
}


bool
InputFront::readColor(emf_ObjectData_X* od, Color4 color)
{
  static UserInterface* gui = theControlCenter->getGui();
  strstream strm;

  if ( LibFront::in(od->data_type, "string") ) {

    char* nm = NULL;
    readName(od, nm);
    bool rc = Model::getColorValue(nm, color);

    if (rc) {
      delete[] nm;
      return true;
    } else {

      strm << "***WARNING Invalid color name: " << nm << ends;
      gui->showMsg(strm.str());
      delete[] nm;
      return false;
    }

  } else {

    int nof_cindices = 4;

    if ( od->data_length < nof_cindices )
      nof_cindices = od->data_length;

    for (short i = 0; i < nof_cindices; i++) {
      LibFront::setNumericData(color[i], i);
    }

    return true;
  }

  return false;
}


// Read one field for an edge
int
InputFront::readEdge(emf_ObjectData_X* od)
{
  static UserInterface* gui = theControlCenter->getGui();
  strstream strm1;
  strstream strm2;

  int rc = isOk;

  const char* fn = od->field_name;

  static ecif_Element_X tx;
  static ecif_geometryType gtype;
  static ecif_DllArg da;

  static bool component_added;
  static bool is_by_geometry;

  static bool is_cpp = true;  // Default!
  static bool is_f95 = false;
  static bool is_matc = false;

  bool isEgf = isEgfInput;
  bool isEmf = isEmfInput;

  //---New element is started
  if (od->is_object_start) {

    // Init tx data
    init_trx_data(tx);
    tx.tag = od->object_id;
    tx.bndr_tag = NO_INDEX;
    tx.tplg_type = ECIF_EDGE;

    gtype = ECIF_NODIM;
    component_added = false;
    is_by_geometry = false;

    is_cpp = false;
    is_f95 = false;
    is_matc = false;
  }

  // Parse fields
  // ============

  //-No data for the element
  if ( od->field_name_length == 0 ) {
    ; // Do nothing!
  }

  //-User function
  else if ( LibFront::in(EMF_FUNCTION, fn) ||
            LibFront::in("FunctionC", fn) ||
            LibFront::in("FunctionF", fn)
          ) {

    char* data = (char*) od->data;
    char* dll_nms[2];

    int pos = 0;
    for (int i = 0; i < od->nof_entries; i++) {
      int len = od->data_lengths[i];
      dll_nms[i] = new char[1+len];
      for (int j = 0; j < len; j++) {
        dll_nms[i][j] = data[pos++];
      }
      dll_nms[i][len] = '\0';
    }

    bool must_have_lib = false;

    if ( LibFront::in("FunctionC", fn) ) {
      must_have_lib = true;
      is_cpp = true;
    }

    else if ( LibFront::in("FunctionF", fn) ) {
      must_have_lib = true;
      is_f95 = true;
    }

    else if ( LibFront::in("FunctionM", fn) ) {
      is_matc = true;
    }

    if ( !is_by_geometry ) {
      strm1 << "***ERROR Invalid function for Edge " << tx.tag << ends;
      strm2 << "No Geometry keyword given (Geometry = Circle/Curve)!" << ends;
      gui->showMsg(strm1.str());
      gui->showMsg(strm2.str());
      goto error;
    }

    if ( od->nof_entries == 1 || dll_nms[1] == NULL || dll_nms[1][0] == '\0' ) {
      if ( must_have_lib ) {
        strm1 << "***ERROR Invalid function definition for Edge " << tx.tag << ends;
        strm2 << "Syntax is: Function Func-name Lib-name!" << ends;
        gui->showMsg(strm1.str());
        gui->showMsg(strm2.str());
        goto error;
      } else {
        is_matc = true;
      }
    }

    if ( !(is_cpp || is_f95 || is_matc)) {
      is_cpp = true;
    }

    // Create new component
    if ( !component_added && !addElementComponent(tx, gtype) ) {
      goto error;
    }

    da.is_cpp = is_cpp;
    da.is_f95 = is_f95;
    da.is_matc = is_matc;

    update_dyna_string(da.func, dll_nms[0]);

    if ( is_cpp || is_f95 ) {
      update_dyna_string(da.lib, dll_nms[1]);
    }

    storeMatcData(da.matcTable, EMF_FUNCTION, od);

    // Pick component and copy dll-argument etc. to the component
    ecif_ElementComponent_X* txc = tx.components[tx.nof_components-1];
    if ( txc != NULL ) {
      txc->isFunction = true;
      addFunctionComponentData(da, *txc);
    } else {
      goto error;
    }

    // Init for the possible next component
    component_added = false;
    da.init();
  }

  //-Function arguments
  else if ( LibFront::in(EMF_ARGUMENTS, (char*)fn) ) {
    da.argc = od->data_length;
    da.argv = new double[od->data_length];
    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(da.argv[i], i);
    }
    storeMatcData(da.matcTable, EMF_ARGUMENTS, od);
  }

  //-Start vertex for a function geometry
  else if ( LibFront::in(EMF_START_VERTEX, (char*)fn) ) {
   LibFront::setNumericData(da.start_vertex, 0);
   storeMatcData(da.matcTable, EMF_START_VERTEX, od);
  }

  //-End vertex for a function geometry
  else if ( LibFront::in(EMF_END_VERTEX, (char*)fn) ) {
   LibFront::setNumericData(da.end_vertex, 0);
   storeMatcData(da.matcTable, EMF_END_VERTEX, od);
  }

  //-Start point for a function geometry
  else if ( LibFront::in(EMF_START_POINT, (char*)fn) ) {
    da.has_start_point = true;
    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(da.start_point[i], i);
      da.start_point[i] *= inputUnit;
    }
    storeMatcData(da.matcTable, EMF_START_POINT, od);
  }

  //-End point for a function geometry
  else if ( LibFront::in(EMF_END_POINT, (char*)fn) ) {
    da.has_end_point = true;
    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(da.end_point[i], i);
      da.end_point[i] *= inputUnit;
    }
    storeMatcData(da.matcTable, EMF_END_POINT, od);
  }

  //-Boundary tag
  else if ( LibFront::in(EMF_BOUNDARY_TAG, fn) ) {
    LibFront::setNumericData(tx.bndr_tag, 0);
  }

  //-Boundary Condition id
  else if ( LibFront::in(EMF_BOUNDARY_CONDITION, fn) ) {
    LibFront::setNumericData(tx.bndr_cond_id, 0);
  }

  //-Boundary Parameter id
  else if ( LibFront::in(EMF_BOUNDARY_PARAMETER, fn) ) {
    LibFront::setNumericData(tx.bndr_param_id, 0);
  }

  //-Name
  else if ( LibFront::in(EMF_NAME, fn) ) {
    //update_dyna_string(tx.name, (char*) od->data);
    readName(od, tx.name);
  }

  //-Boundary group tag
  else if ( LibFront::in(EMF_GROUP, fn) ) {
    LibFront::setNumericData(tx.bndr_group_tag, 0);
  }

  //-Symmetry axis flag
  // NOTE: Not stored in the element!
  else if ( LibFront::in("On Symmetry Axis", fn) ) {
  }

  //-Boundary vertices (= boundary points)
  // NOTE: Not read any more from emf-file (ie. is pre version 5 stuff!)
  else if ( LibFront::in("Boundary Vertices", fn) ) {
  }

  //-Grid parameter ids
  else if ( LibFront::in("Gridh Ids", fn) ) {
    tx.nof_gridh_ids = od->data_length;
    tx.gridh_ids = new int[od->data_length];

    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(tx.gridh_ids[i], i);
    }
  }

  //-GridH mesh indices
  else if ( LibFront::in("Gridh Mesh Indices", fn) ) {
    tx.gridh_mesh_indices = new int[od->data_length];

    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(tx.gridh_mesh_indices[i], i);
    }
  }

  //-Extra vertices (neede mainly for closed geometry and mesh geometry)
  else if ( LibFront::in(EMF_EXTRA_VERTICES, fn) ) {
    tx.nof_extra_vertices = od->data_length;
    tx.extra_vertex_tags = new int[od->data_length];
    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(tx.extra_vertex_tags[i], i);
    }
  }

  //-Linear geometry by vertex tags
  // NOTE: Only linear geometry can be defined without type,
  //       by just giving the vertices!
  else if ( LibFront::in(EMF_VERTICES, fn) ) {

    // If no type given --> Linear, add new linear component
    if ( !is_by_geometry ) {

      if ( od->data_length == 2 ) {
        addElementComponent(tx, ECIF_LINE);

      } else if ( od->data_length > 2 ) {
        addElementComponent(tx, ECIF_POLYLINE);

      } else {
        strm1 << "***ERROR Invalid definition for Edge " << tx.tag << ends;
        strm2 << "At least two vertices must be given!" << ends;
        gui->showMsg(strm1.str());
        gui->showMsg(strm2.str());
        goto error;
      }
    }

    ecif_ElementComponent_X* txc = tx.components[tx.nof_components-1];

    txc->nof_vertices = od->data_length;
    txc->vertex_tags = new int[txc->nof_vertices];
    for (int i = 0; i < od->data_length; i++) {
      int tag;
      LibFront::setNumericData(tag, i);

      if ( !model->checkVertexExistence(tag) ) {
        strm1 << "***ERROR Invalid definition for Edge " << tx.tag << ends;
        strm2 << "Illegal vertex tag: " << tag << ends;
        gui->showMsg(strm1.str());
        gui->showMsg(strm2.str());
        goto error;
      }

      txc->vertex_tags[i] = tag;
    }
  }

  //-Geometry type
  // NOTE: This type must be given before any parametric data
  // NOTE: Only linear geometry can be defined without type,
  // by just giving the vertices!
  else if ( LibFront::in(EMF_GEOMETRY, fn) ) {
    char* tp_nm = NULL;
    readName(od, tp_nm);
    gtype = getElementGeometryType(tp_nm, tx);

    if ( gtype == ECIF_NODIM ) {
      strm1 << "***ERROR Invalid geometry definition for Edge " << tx.tag << ends;
      strm2 << "Unknown geometry type: " << tp_nm << ends;
      gui->showMsg(strm1.str());
      gui->showMsg(strm2.str());
      goto error;
    }

    delete[] tp_nm;

    is_by_geometry = true;

    if ( !addElementComponent(tx, gtype) ) {
      goto error;
    }

    component_added = true;

#if 0
    if ( !component_added ) {
      if ( !addElementComponent(tx, gtype) ) {
        goto error;
      }
      component_added = true;
    }
#endif

#if 0
    component_added = false;
#endif

    storeMatcData(da.matcTable, EMF_GEOMETRY, od);
  }

  // Some edge geometry component
  else {

    if ( !is_by_geometry ) {
      strm1 << "***ERROR Invalid definition for Edge " << tx.tag << ends;
      strm2 << "A geometry component given, but geometry type not defined!" << ends;
      gui->showMsg(strm1.str());
      gui->showMsg(strm2.str());
      goto error;
    }

    if ( !component_added ) {
      if ( !addElementComponent(tx, gtype) ) {
        goto error;
      }
      component_added = true;
    }

    ecif_ElementComponent_X* txc = tx.components[tx.nof_components-1];

    rc = readElementGeometry2D(od, *txc);
  }

  if (rc != isOk) {
    goto error;
  }

  //----Element is ready
  if ( od->is_object_end ) {

    if ( tx.nof_components <= 0 ) {
      strm1 << "***ERROR Invalid definition for Edge " << tx.tag << ends;
      strm2 << "Edge at end, but no actual geometry defined!" << ends;
      //gui->showMsg(strm1.str());
      //gui->showMsg(strm2.str());
      //goto error;

    }

    if ( !model->addBodyElement(tx) ) {
      goto error;
    }

    reset_trx_data(tx);
  }

  // Ok
  return isOk;

  // Error
  error:
  reset_trx_data(tx);
  return notOk;
}


// Read 1D element geometry (vertex)
int
InputFront::readElementGeometry1D(emf_ObjectData_X* od, ecif_ElementComponent_X& tx)
{
  const char* fn = od->field_name;

  bool isEgf = isEgfInput;
  bool isEmf = isEmfInput;

  //-Vertex point (always 3D coordinates)
  if ( LibFront::in("Point", fn) ) {
    for (short j = 0; j < 3; j++) {
      LibFront::setNumericData(tx.geometry.vertex->point[j], j);
      tx.geometry.vertex->point[j] *= inputUnit;
    }
  }

  //-ERROR: Unknown field
  else {
    return unknownFieldMsg(od, IS_FATAL);
  }

  return isOk;
}


// Read 2D element geometry (edge)
int
InputFront::readElementGeometry2D(emf_ObjectData_X* od, ecif_ElementComponent_X& tx)
{
  const char* fn = od->field_name;

  bool isEgf = isEgfInput;
  bool isEmf = isEmfInput;

  //-Radius
  if ( LibFront::in("Radius", (char*)fn) ) {
    LibFront::setNumericData(tx.geometry.edge->radius1, 0);
    tx.geometry.edge->radius1 *= inputUnit;

    storeMatcData(tx.geometry.edge->matcTable, EMF_RADIUS, od);
  }

  //-Center
  else if ( LibFront::in("Center", (char*)fn) ) {
    tx.geometry.edge->location = new Point3[1];
    initPoint3(*tx.geometry.edge->location, 0.0);
    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData((*tx.geometry.edge->location)[i], i);
      (*tx.geometry.edge->location)[i] *= inputUnit;
    }

    storeMatcData(tx.geometry.edge->matcTable, EMF_CENTER, od);
  }

  //-Defining points
  else if ( LibFront::in("Defining Points", (char*)fn) ) {
    tx.geometry.edge->nofDefiningPoints = od->dimension1;
    tx.geometry.edge->location = new Point3[od->dimension1];
    for (int i = 0; i < od->dimension1; i++) {
      for (int j = 0; j < od->dimension2; j) {
        int pos = i * od->dimension1;
        LibFront::setNumericData(tx.geometry.edge->definingPoints[i][j], pos + j);
        tx.geometry.edge->definingPoints[i][j] *= inputUnit;
      }
    }
  }

  //-Delta-h for linearizing (interval length space)
  else if ( LibFront::in(EMF_DELTA_H, fn) ||
            LibFront::in(EMF_MESH_H, fn)
          ) {
    tx.lin_delta_type = LIN_DELTA_H;
    for (int i = 0; i < 1 && i < od->data_length; i++) {
      LibFront::setNumericData(tx.lin_delta[i], i);
    }

    storeMatcData(tx.geometry.edge->matcTable, EMF_MESH_H, od);
  }

  //-Delta-n for linearizing (interval count space)
  else if ( LibFront::in(EMF_DELTA_N, fn) ||
            LibFront::in(EMF_MESH_N, fn)
          ) {
    tx.lin_delta_type = LIN_DELTA_N;
    for (int i = 0; i < 1 && i < od->data_length; i++) {
      LibFront::setNumericData(tx.lin_delta[i], i);
    }

    storeMatcData(tx.geometry.edge->matcTable, EMF_MESH_N, od);
  }

#if 0
  //-Delta-u for linearizing (parameter space)
  else if ( LibFront::in(EMF_DELTA_U, fn) ||
            LibFront::in(EMF_MESH_U, fn)
          ) {
    tx.lin_delta_type = LIN_DELTA_U;
    for (int i = 0; i < 1 && i < od->data_length; i++) {
      LibFront::setNumericData(tx.lin_delta[i], i);
    }
  }
#endif

  // Use (force)  N: setting in meshing
  else if ( LibFront::in(EMF_USE_MESH_N, fn) ) {
    LibFront::setNumericData(tx.use_fixed_mesh_n, 0);
  }

  //-ERROR: Unknown field
  else {
    return unknownFieldMsg(od, IS_FATAL);
  }

  return isOk;
}


// Read element loop
int
InputFront::readElementLoop(emf_ObjectData_X* od)
{
  UserInterface* gui = theControlCenter->getGui();
  const char* fn = od->field_name;

  static ecif_ElementLoop_X tx;

  bool isEgf = isEgfInput;
  bool isEmf = isEmfInput;

  //---New loop is started
  if (od->is_object_start) {

    // Init tx data
    init_trx_data(tx);
    tx.is_checked = geometryIsChecked();
    tx.tag = od->object_id;
  }

  // Parse fields
  // ============

  //-Element tags
  // NOTE: EMF_ELEMENT_IDS: pre version 9 emf-format
  if ( LibFront::in(EMF_ELEMENTS, fn) ||
       LibFront::in(EMF_ELEMENT_IDS, fn) ||
       LibFront::in("Faces", fn) ||
       LibFront::in("Edges", fn)
     ) {
    int size = od->data_length;
    tx.element_tags = new int[size];
    tx.nof_elements = size;
    for (int i = 0; i < size; i++)
      LibFront::setNumericData(tx.element_tags[i], i);
  }

  //-Boundary group tag
  else if ( LibFront::in(EMF_GROUP, fn) ) {
    LibFront::setNumericData(tx.bndr_group_tag, 0);
  }

  //-Loop type
  else if ( LibFront::in(EMF_TYPE, fn) ) {
    if ( LibFront::in("Open", (char*) od->data) ) {
      tx.is_open = true;
    }
  }

  //-Unknown field
  else {
    return unknownFieldMsg(od, IS_FATAL);
  }

  //----ElementLoop is ready
  if (od->is_object_end) {
    model->addBodyElementLoop(tx);
    reset_trx_data(tx);
  }

  return isOk;
}


// Read one field for a face
int
InputFront::readFace(emf_ObjectData_X* od)
{
  static UserInterface* gui = theControlCenter->getGui();
  strstream strm1;
  strstream strm2;
  int i;

  int rc = isOk;

  static IdArray edge_tags;
  static IdArray vertex_tags;

  const char* fn = od->field_name;

  static ecif_Element_X tx;

  //---New element is started
  if (od->is_object_start) {

    // Init tx data
    init_trx_data(tx);
    tx.tag = od->object_id;
    tx.bndr_tag = NO_INDEX;
    tx.tplg_type = ECIF_FACE;
  }

  // Parse fields
  // ============

  //-No data for the element (but id)
  if ( od->field_name_length == 0) {
    ; // Do nothing!
  }

  //-Boundary tag
  else if ( LibFront::in(EMF_BOUNDARY_TAG, fn) ) {
    LibFront::setNumericData(tx.bndr_tag, 0);
  }

  //-Boundary Condition id
  else if ( LibFront::in(EMF_BOUNDARY_CONDITION, fn) ) {
    LibFront::setNumericData(tx.bndr_cond_id, 0);
  }

  //-Boundary Parameter id
  else if ( LibFront::in(EMF_BOUNDARY_PARAMETER, fn) ) {
    LibFront::setNumericData(tx.bndr_param_id, 0);
  }

  //-Boundary name
  else if ( LibFront::in(EMF_NAME, fn) ) {
    //update_dyna_string(tx.name, (char*) od->data);
    readName(od, tx.name);
  }

  //-Edge tags
  else if ( LibFront::in(EMF_EDGES, fn) ) {
    int tag;
    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(tag, i);
      edge_tags.push_back(tag);
    }
  }

  //-Vertex tags
  // NOTE: This is for mesh bodies which may not have any edges
  // but we still want to store some vertices in the model
  //
  else if ( LibFront::in(EMF_VERTICES, fn) ) {
    int tag;
    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(tag, i);
      vertex_tags.push_back(tag);
    }
  }

  //-Grid parameter ids
  else if ( LibFront::in("Gridh Ids", fn) ) {
    tx.nof_gridh_ids = od->data_length;
    tx.gridh_ids = new int[od->data_length];

    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(tx.gridh_ids[i], i);
    }
  }

  //-GridH mesh indices
  else if ( LibFront::in("Gridh Mesh Indices", fn) ) {
    tx.gridh_mesh_indices = new int[od->data_length];

    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(tx.gridh_mesh_indices[i], i);
    }
  }

  if (rc != isOk)
    return !isOk;

  //----Element is ready
  if ( od->is_object_end ) {

    model->addBodyElement(tx);

    // Add possible pending edge definitions
    // into the created face bodyelement
    BodyElement* be = model->getBodyElementByTag(OT_FACE, tx.tag);

    if ( be != NULL ) {

      // Store edges as pending elements (will be resolved later when
      // the whole file is read
      int nof_edges = edge_tags.size();
      if ( nof_edges > 0 ) {
        for (i = 0; i < nof_edges; i++)
          be->addPendingEdge(edge_tags[i]);
      }
      edge_tags.clear();

      // Store vertices as pending elements (will be resolved later when
      // the whole file is read
      int nof_vertices = vertex_tags.size();
      if ( nof_vertices > 0 ) {
        for (i = 0; i < nof_vertices; i++)
          be->addPendingVertex(vertex_tags[i]);
      }
      vertex_tags.clear();
    }

    reset_trx_data(tx);
  }

  return isOk;
}


int
InputFront::readIncludeFile(emf_ObjectData_X* od)
{
  return isOk;
}


bool
InputFront::readName(emf_ObjectData_X* od, char*& name)
{
  name = NULL;

  char* data = (char*) od->data;
  strstream strm;

  int pos = 0;
  for (int i = 0; i < od->nof_entries; i++) {
    int len = od->data_lengths[i];
    for (int j = 0; j < len; j++) {
      strm << data[pos++];
    }

    if ( i < od->nof_entries - 1 ) {
      strm << ' ';
    }
  }

  strm << ends;

  update_dyna_string(name, strm.str());

  return true;
}


void
InputFront::readNumericVector(int read_count, int*& target_vector,
                              int init_count, int init_value)
{
  int i;

  for (i = 0; i < init_count; i++) {
    target_vector[i] = init_value;
  }

  for (i = 0; i < read_count; i++) {
    LibFront::setNumericData(target_vector[i], i);
  }
}


void
InputFront::readNumericVector(int read_count, double*& target_vector,
                             int init_count, double init_value)
{
  int i;

  for (i = 0; i < init_count; i++) {
    target_vector[i] = init_value;
  }

  for (i = 0; i < read_count; i++) {
    LibFront::setNumericData(target_vector[i], i);
  }
}


void
InputFront::readPoint3(int read_count, Point3& point, double init_value)
{
  int i;

  for (i = 0; i < 3; i++) {
    point[i] = init_value;
  }

  for (i = 0; i < read_count; i++) {
    LibFront::setNumericData(point[i], i);
  }
}


// Read one field for a vertex
int
InputFront::readVertex(emf_ObjectData_X* od)
{
  static UserInterface* gui = theControlCenter->getGui();
  strstream strm1;
  strstream strm2;

  int rc = isOk;

  const char* fn = od->field_name;

  static ecif_Element_X tx;

  //---New element is started
  if (od->is_object_start) {

    // Init tx data
    init_trx_data(tx);
    tx.tag = od->object_id;
    tx.bndr_tag = NO_INDEX;
    tx.tplg_type = ECIF_VERTEX;
  }

  // Parse fields
  // ============

  //-No data for the element
  if ( od->field_name_length == 0 ) {
    ; // Do nothing!
  }

  //-Boundary tag
  else if ( LibFront::in(EMF_BOUNDARY_TAG, fn) ) {
    LibFront::setNumericData(tx.bndr_tag, 0);
  }

  //-Boundary Condition id
  else if ( LibFront::in(EMF_BOUNDARY_CONDITION, fn) ) {
    LibFront::setNumericData(tx.bndr_cond_id, 0);
  }

  //-Boundary Parameter id
  else if ( LibFront::in(EMF_BOUNDARY_PARAMETER, fn) ) {
    LibFront::setNumericData(tx.bndr_param_id, 0);
  }

  //-Name
  else if ( LibFront::in(EMF_NAME, fn) ) {
    update_dyna_string(tx.name, (char*) od->data);
  }

  //-Boundary group tag
  else if ( LibFront::in(EMF_GROUP, fn) ) {
    LibFront::setNumericData(tx.bndr_group_tag, 0);
  }

  //-Grid parameter ids
  else if ( LibFront::in("Gridh Ids", fn) ) {
    tx.nof_gridh_ids = od->data_length;
    tx.gridh_ids = new int[od->data_length];

    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(tx.gridh_ids[i], i);
    }
  }

  //-GridH mesh indices
  else if ( LibFront::in("Gridh Mesh Indices", fn) ) {
    tx.gridh_mesh_indices = new int[od->data_length];

    for (int i = 0; i < od->data_length; i++) {
      LibFront::setNumericData(tx.gridh_mesh_indices[i], i);
    }
  }

  // Some vertex geometry component
  else {

    if ( tx.nof_components == 0 ) {
      addElementComponent(tx, ECIF_POINT);
    }

    ecif_ElementComponent_X* txc = tx.components[tx.nof_components-1];

    rc = readElementGeometry1D(od, *txc);
  }

  if (rc != isOk) {
    goto error;
  }

  //----Element is ready
  if ( od->is_object_end ) {

    // NOTE: Check NOT in use!
    if ( false && tx.nof_components <= 0 ) {
      strm1 << "***ERROR Invalid definition for Vertex " << tx.tag << ends;
      strm2 << "Vertex at end, but no geometry point given!" << ends;
      gui->showMsg(strm1.str());
      gui->showMsg(strm2.str());
      goto error;
    }

    model->addBodyElement(tx);
  }

  return isOk;

  // Error
  error:
  return notOk;
}


// Methods read a vertex table
int
InputFront::readVertexTable(emf_ObjectData_X* od)
{
  const ModelInfo* mi = model->getModelInfo();
  int inputVersionNbr = mi->frontInputVersionNbr;

  UserInterface* gui = theControlCenter->getGui();
  const char* fn = od->field_name;

  static ecif_Vertex_X tx;
  static GcPoint point;

  int nofVertices = 0;
  int* vertexIds = NULL;
  MatcValueTable matcTable;

  bool isEgf = isEgfInput;
  bool isEmf = isEmfInput;

  reset_trx_data(tx);

  // Point table (typical Egf format)
  // ================================
  if ( LibFront::in("Points", fn) ) {

    if ( !LibFront::ncEqual(od->data_type, "real") ) {
      gui->showMsg("Incorrect data format for vertex points (should be Real!)");
      return 1;
    }

    int tag, dim1, dim2;
    bool indexed_table;
    double coord[3];
    BodyElement* vertex;

    if ( od->nof_variables == 1 &&
         LibFront::ncEqual(od->variable_names[0], "index")
       ) {
      indexed_table = true;
      dim1 = od->nof_entries;
      dim2 = od->dimension1;
    } else {
      indexed_table = false;
      dim1 = od->dimension1;
      dim2 = od->dimension2;
    }

    vertexIds = new int[dim1];

    storeMatcData(matcTable, EMF_POINTS, od);

    // Read vertex table
    for (int i = 0; i < dim1; i++) {

      coord[0] = coord[1] = coord[2] = 0.0;

      for (short j = 0; j < dim2; j++) {
        LibFront::setNumericData(coord[j], j + i * dim2);
        coord[j] *= inputUnit;
      }

      point.setPosit(coord[0], coord[1], coord[2]);

      if ( indexed_table ) {
        LibFront::setNumericVariable(tag, i);
        vertex = new BodyElement1D(tag, &point);

      } else {
        vertex = new BodyElement1D(&point);
      }

      vertexIds[i] = vertex->Id();

      model->addBodyElement(vertex);
    }

    model->setVertexTable(dim1, dim2, vertexIds, matcTable);
  }

  // Ids and Point table (typical Emf format)
  // ========================================
  // NOTE: New format
  // ----------------
  else if ( LibFront::in("Ids And Points", fn) ) {

    if ( !LibFront::ncEqual(od->data_type, "real") ) {
      gui->showMsg("Incorrect data format for vertex points (should be Real!)");
      return 1;
    }

    int data_size = od->dimension1;

    int var_size = 1; // Id field (index)

    int nof_coord;

    if ( inputVersionNbr < 5 ) {
      nof_coord = data_size;     // Nothing but index variable (id)  before coordinates
    } else {
      nof_coord = data_size - 2; // Boundary tag and boundary condition id added before coordinates
    }

     // Read table
    for (int i = 0; i < od->nof_entries; i++) {

      int var_pos = i * var_size;
      int data_pos = i * data_size;

      // Read id (as index variable!)
      LibFront::setNumericVariable(tx.tag, var_pos);

      if ( inputVersionNbr >= 5 ) {
        LibFront::setNumericData(tx.bndr_tag, data_pos++);
        LibFront::setNumericData(tx.bndr_cond_id, data_pos++);
      }

      // Read coordinates (3D)
      for (short j = 0; j < nof_coord; j++) {
        LibFront::setNumericData(tx.point[j], data_pos++);
        tx.point[j] *= inputUnit;
      }

      model->addBodyElement(tx);
    }
  }

  // Point and Constraint ids table (typical Emf format)
  // ===================================================
  // NOTE: Old format!
  // -----------------
  else if ( LibFront::in("Points And Constraint Ids", fn) ) {
    if ( !LibFront::ncEqual(od->data_type, "real") ) {
      gui->showMsg("Incorrect data format for vertex points (should be Real!)");
      return 1;
    }

    int data_size = od->dimension1;

    int nof_coord = data_size - 2;  // two ids after coordinates
    int var_size = 1;

     // Read vertex table
    for (int i = 0; i < od->nof_entries; i++) {

      int var_pos = i * var_size;
      int data_pos = i * data_size;

      // Read id
      LibFront::setNumericVariable(tx.tag, var_pos);

      // Read coordinates (3D)
      for (short j = 0; j < nof_coord; j++) {
        LibFront::setNumericData(tx.point[j], data_pos++);
        tx.point[j] *= inputUnit;
      }

      // Read boundary condition id
      LibFront::setNumericData(tx.bndr_cond_id, data_pos++);

      // Read boundary id
      LibFront::setNumericData(tx.bndr_tag, data_pos++);

      model->addBodyElement(tx);
    }
  }

  //-Unknown field
  else {
    return unknownFieldMsg(od, IS_FATAL);
  }

  delete[] vertexIds;
  purgeMatcValueTable(matcTable);

  return isOk;
}


// Store field data in Matc-expression table if it contains any Matc-variables
//
bool
InputFront::storeMatcData(MatcValueTable& matcTable, const char* key, emf_ObjectData_X* od)
{
  if ( !od->dataHasMatcVars ) {
    return false;
  }

  // Call external store function
  storeMatcString(matcTable, key, od->dataAsString);

  // Inform gui!
  UserInterface* gui = theControlCenter->getGui();
  gui->setModelHasMatcDefinitions();

  return true;
}


int
InputFront::unknownFieldMsg(emf_ObjectData_X* object_data, bool is_fatal)
{
  UserInterface* gui = theControlCenter->getGui();
  strstream strm;

  if ( is_fatal ) {
    strm << "***ERROR: Unknown field name (";
  } else {
    strm << "***WARNING: Unknown field name (";
  }

  strm << object_data->field_name
       << ") when reading object: "
       << object_data->object_name;

  if ( object_data->object_id != NO_INDEX ) {
    strm << " " << object_data->object_id;
  }

  strm << ends;

  gui->showMsg(strm.str());

  if (is_fatal)
    return notOk;
  else
    return isOk;
}


int
InputFront::unknownObjectMsg(emf_ObjectData_X* object_data, bool is_fatal)
{
  UserInterface* gui = theControlCenter->getGui();
  strstream strm;

  if ( is_fatal ) {
    strm << "***ERROR: Unknown object name:";
  } else {
    strm << "***WARNING: Unknown object name:";
  }

  strm << object_data->object_name
       << ends;

  gui->showMsg(strm.str());

  if (is_fatal)
    return notOk;
  else
    return isOk;
}
