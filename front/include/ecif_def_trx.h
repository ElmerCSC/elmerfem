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
Module:     ecif_def_trx.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Header file for the interface for reading ecif-output file.
  All necessary data structures, functions and global variables needed
  for the interface are defined/declared here.


************************************************************************/

#ifndef _ECIF_DEF_TRX_
#define _ECIF_DEF_TRX_

#include "ecif_def_stl.h"

// ********
// TYPEDEFS
// ********
typedef double ecif_Point3[3];
typedef double ecif_Point4[4];

typedef double egf_Point3[3];
typedef double egf_Point4[4];

typedef double Point3[3];
typedef double Point4[4];

#define ecif_MAX_NAME_LEN 128

// Max number of vertices per face/edge
#define ecif_MAX_NOF_VERTICES 1024

// Max number of control/knots points per nurbs face/edge
#define ecif_MAX_NOF_CPOINTS 8192
#define ecif_MAX_NOF_KNOTS 8200



// This struct stores info for a egf user function call
struct ecif_DllArg
{
  ecif_DllArg() {
    init();
  }

  ~ecif_DllArg() {
    delete[] func;
    delete[] lib;
    delete[] argv;
    purgeMatcValueTable(matcTable);
  }

  void init() {
    is_cpp = true;
    is_f95 = false;
    is_matc = false;
    gtype = ECIF_NODIM;
    func = NULL;
    lib = NULL;
    argc = 0;
    argv = NULL;
    start_vertex = end_vertex = NO_INDEX;
    has_start_point = has_end_point = false;
    start_point[0] = start_point[1] = start_point[2] = 0.0;
    end_point[0] = end_point[1] = end_point[2] = 0.0;
    purgeMatcValueTable(matcTable);
  }
  enum ecif_geometryType gtype;
  bool is_cpp;
  bool is_f95;
  bool is_matc;
  char* func;
  char* lib;
  int argc;
  double* argv;
  int start_vertex;
  int end_vertex;
  bool has_start_point;
  bool has_end_point;
  egf_Point3 start_point;
  egf_Point3 end_point;
  MatcValueTable matcTable;
};


// This struct transfers body related data
struct ecif_Body_X
{
  int tag;
  char* name;
  int color[4];
  bool has_color;
  bool is_checked;
  bool is_bem;
  bool is_open;
  bool is_virtual;

  int body_force_id;
  int body_param_id;
  int equation_id;
  int init_cond_id;
  int material_id;

  int nof_layers;
  int* layer_tags;
};


// This struct transfers body grid group data
struct ecif_BodyLayer_X
{
  int tag;
  char* name;
  int color[4];
  bool has_name;
  bool has_color;
  bool is_open;

  int body_id;
  int body_tag;

  int nof_elem_loops;
  int* elem_loop_tags;

  int nof_grid_param_ids;
  int* grid_param_ids;
  int* grid_param_mesh_indices;

  int nof_excluded_meshes;
  int* excluded_mesh_indices;
  
  // For virtual bodies actually
  int nof_elem_groups;
  int* elem_group_tags;
};


// This is a generic structure for edge-type geometry
struct ecif_EdgeGeometry_X
{
  // Generic for all
  enum ecif_geometryType type;
  Point3* start;
  Point3* end;
  bool isClosed;
  int nofDefiningPoints;
  Point3* definingPoints;
  bool* pointVertexFlags; // Which of the def. points should become a vertex

  // Line specific
  bool onSymmAxis;

  // Conic specific
  Point3* location;
  Point3 direction;
  double radius1;
  double radius2;
  double apex;
  double focalLength;

  // Nurbs/Spline specific
  bool isRational;
  int degree;

  int nofKnots;
  double* knots;

  int nofCpoints;
  Point4* cpoints; //NOTE: is (nofCpoints * 4)-array

  MatcValueTable matcTable;
};


// This is a generic structure for face-type geometry
struct ecif_FaceGeometry_X
{
  // Generic
  enum ecif_geometryType type;
  int nofDefiningPoints;
  Point3* definingPoints;
  bool isClosed;

  // Plane specific
  bool onSymmPlane;

  // Conic specific
  Point3* location;
  Point3 direction;
  double radius1;
  double radius2;
  double radius3;

  // Nurbs/Spline specific
  bool isRational;
  int degree_u;
  int degree_v;

  int nofKnots_u;
  int nofKnots_v;
  double* knots_u;
  double* knots_v;

  int nofCpoints_u;
  int nofCpoints_v;
  int nofCpoints;
  Point4* cpoints; //NOTE: is (nofCpoints * 4)-array

  MatcValueTable matcTable;
};


union ecif_Geometry_X
{
  // 1D
  struct ecif_VertexGeometry_X* vertex;
  // 2D
  struct ecif_EdgeGeometry_X* edge;
  // 3D/
  struct ecif_FaceGeometry_X* face;
};


// This structure transfers element-loop info
// A loop defines a boundary for a body
// It can be a loop of edges or faces
struct ecif_ElementLoop_X
{
  int tag;
  bool is_checked;
  bool is_open;
  int nof_elements;
  int* element_tags;
  int bndr_group_tag;
};


// This struct transfers bodyelement (boundary) info
struct ecif_ElementComponent_X
{
  double lin_delta[2];
  enum linDeltaType lin_delta_type;
  int use_fixed_mesh_n;

  int nof_vertices;
  int* vertex_tags;

  enum ecif_geometryType gmtr_type;
  union ecif_Geometry_X geometry;

  // Function related stuff
  bool isFunction;
  bool isCpp;
  bool isF95;
  bool isMatc;
  int argc;
  double* argv;
  Point3* startPoint;
  Point3* endPoint;
  int startVertex;
  int endVertex;
  char* functionName;
  char* libraryName;
  MatcValueTable matcTable;
};


// This struct transfers bodyelement (boundary) info
struct ecif_Element_X
{
  int tag;
  int bndr_tag;
  int bndr_cond_id;
  int bndr_param_id;
  char* name;
  
  ecif_topologyType tplg_type;
  int nof_components;
  ecif_ElementComponent_X** components;
  
  int nof_extra_vertices;
  int* extra_vertex_tags;

  int nof_gridh_ids;
  int* gridh_ids;
  int* gridh_mesh_indices;

  int bndr_group_tag;
};


// This struct transfers info for a 3D-point
// NOTE This is for: inputVersionNbr <= 4
struct ecif_Vertex_X
{
  int tag;
  int bndr_tag;
  int bndr_cond_id;
  int* grid_h_ids;
  int nof_grid_h_ids;
  int* mesh_indices;
  Point3 point;

  int bndr_group_tag;
};


// This struct transfers 1D geoemtry (point)
// NOTE This is for: inputVersionNbr >= 5
struct ecif_VertexGeometry_X
{
  Point3 point;
};


// This struct transfers element group data
struct ecif_ElementGroup_X
{
  int tag;
  bool is_virtual;
  char* name;
  bool has_name;

  int nof_elements;
  int* element_tags;

  int boundary_cond_id;
  int boundary_param_id;
};


// This struct transfers info for a outer boundary
struct ecif_OuterBoundary_X
{
  int body_tag;
  int elem_tag;
};


// This struct transfers info for a inner boundary
struct ecif_InnerBoundary_X
{
  int body1_tag;
  int body2_tag;
  int elem_tag;
};



// This struct transfers mesh_h parameter data
struct ecif_MeshParameter_X
{
  int id;
  double mesh_h;
};



// Ecif-trx struc handling functions
// =================================

// Initialize
extern void init_trx_data(ecif_Body_X& tx);
extern void init_trx_data(ecif_BodyLayer_X& tx);
extern void init_trx_data(ecif_Element_X& tx);
extern void init_trx_data(ecif_ElementComponent_X& tx);
extern void init_trx_data(ecif_ElementLoop_X& tx);
extern void init_trx_data(ecif_Geometry_X& tx);
extern void init_trx_data(ecif_VertexGeometry_X& tx);
extern void init_trx_data(ecif_EdgeGeometry_X& tx);
extern void init_trx_data(ecif_FaceGeometry_X& tx);
extern void init_trx_data(ecif_Vertex_X& tx);
extern void init_trx_data(ecif_ElementGroup_X& tx);

// Reset (also delete any dynamically allocated data)
extern void reset_trx_data(ecif_Body_X& tx);
extern void reset_trx_data(ecif_BodyLayer_X& tx);
extern void reset_trx_data(ecif_Element_X& tx);
extern void reset_trx_data(ecif_ElementComponent_X& tx, ecif_topologyType tplg_type);
extern void reset_trx_data(ecif_ElementLoop_X& tx);
extern void reset_trx_data(ecif_Geometry_X& tx, ecif_topologyType tplg_type);
extern void reset_trx_data(ecif_VertexGeometry_X& tx);
extern void reset_trx_data(ecif_EdgeGeometry_X& tx);
extern void reset_trx_data(ecif_FaceGeometry_X& tx);
extern void reset_trx_data(ecif_Vertex_X& tx);
extern void reset_trx_data(ecif_ElementGroup_X& tx);


#endif

