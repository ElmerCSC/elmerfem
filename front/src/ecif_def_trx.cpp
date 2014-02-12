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
Module:     ecif_def_trx.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation of ecif_***_X struct handling routines

************************************************************************/


#include "ecif_def.h"
#include "ecif_def_stl.h"
#include "ecif_def_trx.h"

// Body
// ====
void init_trx_data(ecif_Body_X& tx)
{
  tx.tag = NO_INDEX;
  tx.name = NULL;
  tx.color[0] = tx.color[1] = tx.color[2] = 0;
  tx.color[3] = 255;
  tx.has_color = false;

  tx.is_checked = false;
  tx.is_bem = false;
  tx.is_open = false;
  tx.is_virtual = false;

  tx.body_param_id = NO_INDEX;
  tx.equation_id = NO_INDEX;
  tx.body_force_id = NO_INDEX;
  tx.init_cond_id = NO_INDEX;
  tx.material_id = NO_INDEX;

  tx.nof_layers = 0;
  tx.layer_tags = NULL;
};


void reset_trx_data(ecif_Body_X& tx)
{
  delete[] tx.name;
  delete[] tx.layer_tags;

  init_trx_data(tx);
};


// Body Layer
// ==========
void init_trx_data(ecif_BodyLayer_X& tx)
{
  tx.tag = NO_INDEX;
  tx.name = NULL;
  tx.color[0] = tx.color[1] = tx.color[2] = 0;
  tx.color[3] = 255;
  tx.has_color = false;
  tx.is_open = false;

  tx.body_id = NO_INDEX;
  tx.body_tag = NO_INDEX;

  tx.nof_elem_loops = 0;
  tx.elem_loop_tags = NULL;

  tx.nof_grid_param_ids = 0;
  tx.grid_param_ids = NULL;
  tx.grid_param_mesh_indices = NULL;

  tx.nof_excluded_meshes = 0;
  tx.excluded_mesh_indices = NULL;

  tx.nof_elem_groups = 0;
  tx.elem_group_tags = NULL;

};


void reset_trx_data(ecif_BodyLayer_X& tx)
{
  delete[] tx.name;
  delete[] tx.elem_loop_tags;
  delete[] tx.grid_param_ids;
  delete[] tx.grid_param_mesh_indices;
  delete[] tx.excluded_mesh_indices;

  delete[] tx.elem_group_tags;

  init_trx_data(tx);
};


// Element Group
// =============
void init_trx_data(ecif_ElementGroup_X& tx)
{
  tx.tag = NO_INDEX;
  tx.is_virtual = false;
  tx.nof_elements = 0;
  tx.element_tags = NULL;
  tx.name = NULL;
  tx.has_name = false;
  tx.boundary_cond_id = NO_INDEX;
  tx.boundary_param_id = NO_INDEX;
};


void reset_trx_data(ecif_ElementGroup_X& tx)
{
  delete[] tx.element_tags;
  delete[] tx.name;
  init_trx_data(tx);
};


// Geometry
// ========
void init_trx_data(ecif_Geometry_X& tx)
{
  tx.vertex = NULL;
  tx.edge = NULL;
  tx.face = NULL;
};


void reset_trx_data(ecif_Geometry_X& tx, ecif_topologyType tplg_type)
{
  if ( tplg_type == ECIF_VERTEX && tx.vertex != NULL) {
    reset_trx_data(*tx.vertex);

  } else if ( tplg_type == ECIF_EDGE && tx.edge != NULL ) {
    reset_trx_data(*tx.edge);

  } else if ( tplg_type == ECIF_FACE && tx.face != NULL ) {
    reset_trx_data(*tx.face);
  }

  init_trx_data(tx);
};


// Vertex geometry
void init_trx_data(ecif_VertexGeometry_X& tx)
{
  tx.point[0] = tx.point[1] = tx.point[2] = 0.0;
};


void reset_trx_data(ecif_VertexGeometry_X& tx)
{
  init_trx_data(tx);
};


// Edge geometry
void init_trx_data(ecif_EdgeGeometry_X& tx)
{
  // Generic for all
  tx.type = ECIF_NODIM;
  tx.start = NULL;
  tx.end = NULL;
  tx.isClosed = true;

  tx.nofDefiningPoints = 0;
  tx.definingPoints = NULL;
  tx.pointVertexFlags = NULL;

  tx.onSymmAxis = false;

  tx.location = NULL;

  // Default direction towards positive Y-axis
  tx.direction[0] = 0;
  tx.direction[1] = 1;
  tx.direction[2] = 0;

  tx.radius1 = 0;
  tx.radius2 = 0;
  tx.apex = 0;
  tx.focalLength = 0;

  // Nurbs/Spline specific
  tx.isRational = true;
  tx.degree = 0;

  tx.nofKnots = 0;
  tx.knots = NULL;

  tx.nofCpoints = 0;
  tx.cpoints = NULL; //NOTE: is (nofCpoints,4)-array

 };

void reset_trx_data(ecif_EdgeGeometry_X& tx)
{
  delete[] tx.start;
  delete[] tx.end;
  delete[] tx.definingPoints;
  delete[] tx.pointVertexFlags;
  delete[] tx.location;
  delete[] tx.knots;
  delete[] tx.cpoints;
  purgeMatcValueTable(tx.matcTable);
  init_trx_data(tx);
}

// Face geometry (not ready yet!)
void init_trx_data(ecif_FaceGeometry_X& tx)
{
  // Generic
  tx.type = ECIF_NODIM;
  tx.nofDefiningPoints = 0;
  tx.definingPoints = NULL;
  tx.isClosed = true;

  // Plane specific
  tx.onSymmPlane = false;

  tx.radius1 = 0;
  tx.radius2 = 0;
  tx.radius3 = 0;

  tx.location = NULL;
  // Default direction towards positive Y-plane
  tx.direction[0] = 0;
  tx.direction[1] = 1;
  tx.direction[2] = 0;

  // Nurbs/Spline specific
  tx.isRational = true;
  tx.degree_u = 0;
  tx.degree_v = 0;

  tx.nofKnots_u = 0;
  tx.nofKnots_v = 0;
  tx.knots_u = NULL;
  tx.knots_v = NULL;

  tx.nofCpoints_u = 0;
  tx.nofCpoints_v = 0;
  tx.nofCpoints = 0;
  tx.cpoints = NULL; //NOTE: is (nofCpoints,4)-array

};

void reset_trx_data(ecif_FaceGeometry_X& tx)
{
  delete[] tx.definingPoints;
  delete[] tx.location;
  delete[] tx.knots_u;
  delete[] tx.knots_v;
  delete[] tx.cpoints;
  purgeMatcValueTable(tx.matcTable);
  init_trx_data(tx);
}



// ElementLoop
// ===========
void init_trx_data(ecif_ElementLoop_X& tx)
{
  tx.tag = NO_INDEX;
  tx.is_checked = false;
  tx.is_open = false;
  tx.nof_elements = 0;
  tx.element_tags = NULL;
  tx.bndr_group_tag = NO_INDEX;
};


void reset_trx_data(ecif_ElementLoop_X& tx)
{
  delete[] tx.element_tags;
  init_trx_data(tx);
};


// Element
// =======
void init_trx_data(ecif_Element_X& tx)
{
  tx.tag = NO_INDEX;
  tx.bndr_cond_id = NO_INDEX;
  tx.bndr_param_id = NO_INDEX;
  tx.name = NULL;
  tx.bndr_group_tag = NO_INDEX;

  tx.nof_extra_vertices = 0;
  tx.extra_vertex_tags = 0;

  tx.nof_gridh_ids = 0;
  tx.gridh_ids = NULL;
  tx.gridh_mesh_indices = NULL;

  tx.nof_components = 0;
  tx.components = NULL;
};


void reset_trx_data(ecif_Element_X& tx)
{
  delete[] tx.extra_vertex_tags;

  delete[] tx.gridh_ids;
  delete[] tx.gridh_mesh_indices;
  delete[] tx.name;

  for (int i = 0; i < tx.nof_components; i++) {
    reset_trx_data(*tx.components[i], tx.tplg_type);
    delete tx.components[i];
  }
  delete[] tx.components;
  init_trx_data(tx);
};


// Element component
// =================
void init_trx_data(ecif_ElementComponent_X& tx)
{
  tx.gmtr_type = ECIF_NODIM;
  tx.nof_vertices = 0;
  tx.vertex_tags = NULL;

  tx.lin_delta[0] = tx.lin_delta[1] = -1.0 ;
  tx.lin_delta_type = LIN_DELTA_NONE;
  tx.use_fixed_mesh_n = 0;

  init_trx_data(tx.geometry);

 // Function related
  tx.isFunction = false;
  tx.isCpp = true;
  tx.isF95 = false;
  tx.isMatc = false;
  tx.argc = 0;
  tx.argv = NULL;
  tx.startPoint = NULL;
  tx.endPoint = NULL;
  tx.startVertex = NO_INDEX;
  tx.endVertex = NO_INDEX;
  tx.functionName = NULL;
  tx.libraryName = NULL;
};


void reset_trx_data(ecif_ElementComponent_X& tx, ecif_topologyType tplg_type)
{
  delete[] tx.vertex_tags;
  delete[] tx.argv;
  delete tx.startPoint;
  delete tx.endPoint;
  delete[] tx.functionName;
  delete[] tx.libraryName;
  reset_trx_data(tx.geometry, tplg_type);
  purgeMatcValueTable(tx.matcTable);
  init_trx_data(tx);
};


// Vertex
// ======
void init_trx_data(ecif_Vertex_X& tx)
{
  tx.tag = NO_INDEX;
  tx.bndr_cond_id = NO_INDEX;
  tx.grid_h_ids = NULL;
  tx.nof_grid_h_ids = 0;
  tx.mesh_indices = NULL;
  tx.point[0] = tx.point[1] = tx.point[2] = 0.0;
  tx.bndr_group_tag = NO_INDEX;
};


void reset_trx_data(ecif_Vertex_X& tx)
{
  delete[] tx.mesh_indices;
  delete[] tx.grid_h_ids;
  init_trx_data(tx);
};
