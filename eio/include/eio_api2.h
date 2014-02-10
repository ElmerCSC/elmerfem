/*  
   Elmer, A Finite Element Software for Multiphysical Problems
   Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland

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
Program:    ELMER Data base interface (EIO)
Author(s):  Harri Hakula 10.03.98
************************************************************************/

#ifndef EIOAPI_H
#define EIOAPI_H

#include "../config.h"

#ifdef __cplusplus
extern "C"
{
#endif		/* __cplusplus */
  

//#if defined(__cplusplus)
//#define EIOFUN extern "C" void
#define IREF int&
#define DREF double&
  //#else
  //#define EIOFUN void
  //#define IREF int*
  //#define DREF double*
  //#endif

void FC_FUNC_(eio_init,eio_init) (IREF info);
void FC_FUNC_(eio_init_parallel,eio_init_parallel) (IREF procs, IREF me, IREF info);
void FC_FUNC_(eio_close,eio_close) (IREF info);

void FC_FUNC_(eio_create_model,eio_create_model) (const char *directory, IREF info);
void FC_FUNC_(eio_open_model,eio_open_model) (const char *directory, IREF info);
void FC_FUNC_(eio_close_model,eio_close_model) (IREF info);


void FC_FUNC_(eio_create_geometry,eio_create_geometry) (IREF info);
void FC_FUNC_(eio_open_geometry,eio_open_geometry) (IREF info);
void FC_FUNC_(eio_close_geometry,eio_close_geometry) (IREF info);
void FC_FUNC_(eio_set_geometry_description,eio_set_geometry_description) (IREF bodyC, IREF boundaryC, IREF outerC, 
				    IREF innerC, IREF vertexC, 
				    IREF loopC, IREF maxLooplen, IREF info);
void FC_FUNC_(eio_get_geometry_description,eio_get_geometry_description) (IREF bodyC, IREF boundaryC, IREF outerC, 
				    IREF innerC, IREF vertexC, 
				    IREF loopC, IREF maxLooplen, IREF info);
void FC_FUNC_(eio_set_geometry_body,eio_set_geometry_body) (IREF tag, IREF meshControl, IREF loopC,
			     int *loops,
			     IREF info);
void FC_FUNC_(eio_get_geometry_body,eio_get_geometry_body) (IREF tag, IREF meshControl, IREF loopC, 
			     int *loops,
			     IREF info);
void FC_FUNC_(eio_set_geometry_body_loop,eio_set_geometry_body_loop) (IREF tag, IREF field, int *nodes, IREF info);
void FC_FUNC_(eio_get_geometry_body_loop,eio_get_geometry_body_loop) (IREF tag, IREF field, int *nodes, IREF info);

// Added nodeC argument: Martti Verho 17.03.99
void FC_FUNC_(eio_set_geometry_element,eio_set_geometry_element) (IREF tag, IREF cTag, IREF meshControl,
				IREF type, IREF nodeC, int *nodes, IREF info);

// Added nodeC argument: Martti Verho 17.03.99
void FC_FUNC_(eio_get_geometry_element,eio_get_geometry_element) (IREF tag, IREF cTag, IREF meshControl,
				IREF type, IREF nodeC, int *nodes, IREF info);

// Added: Martti Verho 17.03.99
void FC_FUNC_(eio_get_geometry_element_description,eio_get_geometry_element_description) (IREF tag, IREF cTag, IREF meshControl,
				     IREF type, IREF nodeC, IREF info);
void FC_FUNC_(eio_set_geometry_node,eio_set_geometry_node) (IREF tag, IREF cTag, double *coord, IREF info);
void FC_FUNC_(eio_get_geometry_node,eio_get_geometry_node) (IREF tag, IREF cTag, double *coord, IREF info);
void FC_FUNC_(eio_set_geometry_boundary,eio_set_geometry_boundary) (IREF tag, IREF left, IREF right, IREF info);
void FC_FUNC_(eio_get_geometry_boundary,eio_get_geometry_boundary) (IREF tag, IREF left, IREF right, IREF info);

void FC_FUNC_(eio_create_mesh,eio_create_mesh) (const char *directory, IREF info);
void FC_FUNC_(eio_open_mesh,eio_open_mesh) (const char *directory, IREF info);
void FC_FUNC_(eio_close_mesh,eio_close_mesh) (IREF info);

void FC_FUNC_(eio_set_mesh_description,eio_set_mesh_description) (IREF nodeCount, IREF elementCount, 
				IREF boundaryElementCount,
				IREF usedElementTypes, int* elementTypeTags,
				int *elementCountByType, 
				IREF info);
void FC_FUNC_(eio_set_mesh_node,eio_set_mesh_node) (IREF tag, IREF constraint, double *coord, IREF info);

void FC_FUNC_(eio_set_mesh_element_conns,eio_set_mesh_element_conns) (IREF tag, IREF body, 
    IREF type, int *nodes, 
				   IREF info);
/*
 *
 * Obsolete versions of these routines
void FC_FUNC_(eio_set_mesh_bndry_element,eio_set_mesh_bndry_element) (IREF tag, IREF constraint, 
				     IREF leftbody, IREF rightbody,
				     IREF left, IREF right,
			   double *auxCoord,
				     IREF type, int *nodes,
				     IREF info);
void FC_FUNC_(eio_get_mesh_bndry_element,eio_get_mesh_bndry_element) (IREF tag, IREF constraint, 
				     IREF leftbody, IREF rightbody,
				     IREF left, IREF right,
				     double *auxCoord, IREF type, int *nodes,
				     double *coord, IREF info);
				     */
/*
 *  Modified routines
 */
void FC_FUNC_(eio_set_mesh_bndry_element,eio_set_mesh_bndry_element) (IREF tag, IREF boundary,
				     IREF leftElement, IREF rightElement,
				     IREF type, int *nodes,
				     IREF info);
void FC_FUNC_(eio_get_mesh_bndry_element,eio_get_mesh_bndry_element) (IREF tag, IREF boundary,
				     IREF leftElement, IREF rightElement,
				     IREF type, int *nodes,
				     double *coord, IREF info);
void FC_FUNC_(eio_get_mesh_description,eio_get_mesh_description) (IREF nodeCount, IREF elementCount, 
				IREF boundaryElementCount,
				IREF usedElementTypes, int* elementTypeTags,
				int *elementCountByType, 
				IREF info);
void FC_FUNC_(eio_get_mesh_element_conns,eio_get_mesh_element_conns) (IREF tag, IREF body, IREF type, 
		int *pdofs, int *nodes, IREF info);
void FC_FUNC_(eio_get_mesh_element_coords,eio_get_mesh_element_coords) (IREF tag, IREF body, IREF type, 
					int *nodes, 
					double *coord, IREF info);
void FC_FUNC_(eio_get_mesh_nodes,eio_get_mesh_nodes) (int *tags, double *coord, IREF info);

void FC_FUNC_(eio_create_dual_mesh,eio_create_dual_mesh) (const char *dir, IREF info);
void FC_FUNC_(eio_open_dual_mesh,eio_open_dual_mesh) (const char *dir, IREF info);
void FC_FUNC_(eio_close_dual_mesh,eio_close_dual_mesh) (IREF info);

void FC_FUNC_(eio_set_dual_mesh_element_conns,eio_set_dual_mesh_element_conns) (IREF tag, 
					     IREF type, int *nodes, 
					     IREF info);
void FC_FUNC_(eio_get_dual_mesh_element_conns,eio_get_dual_mesh_element_conns) (IREF tag, 
					     IREF type, int *nodes, 
					     IREF info);

void FC_FUNC_(eio_create_part,eio_create_part) (const char *dir, IREF parts, IREF info);
void FC_FUNC_(eio_open_part,eio_open_part) (IREF info);
void FC_FUNC_(eio_close_part,eio_close_part) (IREF info);
void FC_FUNC_(eio_set_part_description,eio_set_part_description) (IREF nodeCount, 
			 IREF sharedNodeCount,
			 IREF elementCount, 
			 IREF borderElementCount, 
			 IREF boundaryElementCount, 
			 IREF usedElementTypes, 
			 int* elementTypeTags,
			 int *elementCountByType,
			 IREF info);
EIOFUN
eio_get_part_description(IREF sharedNodeCount,
			 IREF info);
void FC_FUNC_(eio_activate_part_part,eio_activate_part_part) (IREF part, IREF info);
void FC_FUNC_(eio_deactivate_part_part,eio_deactivate_part_part) (IREF info);
void FC_FUNC_(eio_set_part_node,eio_set_part_node) ( IREF tag, 
		   IREF constraint,      
		   double *coord, 
		   IREF partcount, 
		   int *parts,     
		   IREF info);
void FC_FUNC_(eio_get_part_node,eio_get_part_node) ( IREF tag, 
		   IREF constraint,      
		   double *coord, 
		   IREF partcount, 
		   int *parts,     
		   IREF info);
void FC_FUNC_(eio_set_part_element,eio_set_part_element) (IREF tag, 
				    IREF body, 
				    IREF type, 
				    int *nodes,
				    IREF border,
				    IREF info);
/* MODEL DATA */

void FC_FUNC_(eio_create_modeldata,eio_create_modeldata) (IREF info);
void FC_FUNC_(eio_open_modeldata,eio_open_modeldata) (IREF info);
void FC_FUNC_(eio_close_modeldata,eio_close_modeldata) (IREF info);
void FC_FUNC_(eio_set_modeldata_description,eio_set_modeldata_description) (IREF bodies,
				     IREF body_forces,
				     IREF body_equations,
				     IREF materials,
				     IREF boundary_conditions,
				     IREF initial_conditions,
				     IREF mesh_parameters,
				     IREF info);

void FC_FUNC_(eio_get_modeldata_description,eio_get_modeldata_description) (IREF bodies,
				     IREF body_forces,
				     IREF body_equations,
				     IREF materials,
				     IREF boundary_conditions,
				     IREF initial_conditions,
				     IREF mesh_parameters,
				     IREF info);
void FC_FUNC_(eio_set_body,eio_set_body) (IREF tag, IREF body_force_id,
			  IREF equation_id, IREF init_cond_id,
			  IREF material_id, IREF mesh_param_id, IREF info);
void FC_FUNC_(eio_get_body,eio_get_body) (IREF tag, IREF body_force_id,
			  IREF equation_id, IREF init_cond_id,
			  IREF material_id, IREF mesh_param_id, IREF info);
void FC_FUNC_(eio_set_constants,eio_set_constants) (double* gravity, DREF boltz, IREF info);
void FC_FUNC_(eio_get_constants,eio_get_constants) (double* gravity, DREF boltz, IREF info);
void FC_FUNC_(eio_set_coords,eio_set_coords) (IREF dim, IREF coordsys, int *mapping,
				  IREF symmetry,
				  double *start,
				  double *end1, double* end2, IREF info);
void FC_FUNC_(eio_get_coords,eio_get_coords) (IREF dim, IREF coordsys, int *mapping,
				  IREF symmetry,
				  double *start,
				  double *end1, double* end2, IREF info);
void FC_FUNC_(eio_set_material_head,eio_set_material_head) (IREF tag, IREF fields, IREF info);
void FC_FUNC_(eio_set_material_field,eio_set_material_field) (IREF name,
				     IREF type, IREF len,
				     int *fields,
				     double *values,
				     IREF info);
void FC_FUNC_(eio_set_bndry_condition_head,eio_set_bndry_condition_head) (IREF tag, 
					      IREF fields, 
					      IREF info);
void FC_FUNC_(eio_set_bndry_condition_field,eio_set_bndry_condition_field) (IREF name,
				     IREF type, IREF len,
				     int *fields,
				     double *values,
				     IREF info);
void FC_FUNC_(eio_set_initial_condition_head,eio_set_initial_condition_head) (IREF tag, 
					      IREF fields, 
					      IREF info);
void FC_FUNC_(eio_set_initial_condition_field,eio_set_initial_condition_field) (IREF name,
				     IREF type, IREF len,
				     int *fields,
				     double *values,
				     IREF info);
void FC_FUNC_(eio_set_body_equation_head,eio_set_body_equation_head) (IREF tag, 
					      IREF fields, 
					      IREF info);
void FC_FUNC_(eio_set_body_equation_field,eio_set_body_equation_field) (IREF name,
				     IREF type, IREF len,
				     int *fields,
				     double *values,
				     IREF info);
void FC_FUNC_(eio_set_body_force_head,eio_set_body_force_head) (IREF tag, 
					      IREF fields, 
					      IREF info);
void FC_FUNC_(eio_set_body_force_field,eio_set_body_force_field) (IREF name,
				     IREF type, IREF len,
				     int *fields,
				     double *values,
				     IREF info);
void FC_FUNC_(eio_set_mesh_parameter_head,eio_set_mesh_parameter_head) (IREF tag, 
					      IREF fields, 
					      IREF info);
void FC_FUNC_(eio_set_mesh_parameter_field,eio_set_mesh_parameter_field) (IREF name,
				     IREF type, IREF len,
				     int *fields,
				     double *values,
				     IREF info);

void FC_FUNC_(eio_get_material_head,eio_get_material_head) (IREF tag, IREF fields, IREF info);
void FC_FUNC_(eio_get_material_field,eio_get_material_field) (IREF name,
				     IREF type, IREF len,
				     int *fields,
				     double *values,
				     IREF info);
void FC_FUNC_(eio_get_bndry_condition_head,eio_get_bndry_condition_head) (IREF tag, 
					      IREF fields, 
					      IREF info);
void FC_FUNC_(eio_get_bndry_condition_field,eio_get_bndry_condition_field) (IREF name,
				     IREF type, IREF len,
				     int *fields,
				     double *values,
				     IREF info);
void FC_FUNC_(eio_get_initial_condition_head,eio_get_initial_condition_head) (IREF tag, 
					      IREF fields, 
					      IREF info);
void FC_FUNC_(eio_get_initial_condition_field,eio_get_initial_condition_field) (IREF name,
				     IREF type, IREF len,
				     int *fields,
				     double *values,
				     IREF info);
void FC_FUNC_(eio_get_body_equation_head,eio_get_body_equation_head) (IREF tag, 
					      IREF fields, 
					      IREF info);
void FC_FUNC_(eio_get_body_equation_field,eio_get_body_equation_field) (IREF name,
				     IREF type, IREF len,
				     int *fields,
				     double *values,
				     IREF info);
void FC_FUNC_(eio_get_body_force_head,eio_get_body_force_head) (IREF tag, 
					      IREF fields, 
					      IREF info);
void FC_FUNC_(eio_get_body_force_field,eio_get_body_force_field) (IREF name,
				     IREF type, IREF len,
				     int *fields,
				     double *values,
				     IREF info);
void FC_FUNC_(eio_get_mesh_parameter_head,eio_get_mesh_parameter_head) (IREF tag, 
					      IREF fields, 
					      IREF info);
void FC_FUNC_(eio_get_mesh_parameter_field,eio_get_mesh_parameter_field) (IREF name,
				     IREF type, IREF len,
				     int *fields,
				     double *values,
				     IREF info);

void FC_FUNC_(eio_create_solver,eio_create_solver) (IREF info);
void FC_FUNC_(eio_open_solver,eio_open_solver) (IREF info);
void FC_FUNC_(eio_close_solver,eio_close_solver) (IREF info);
void FC_FUNC_(eio_set_solver_description,eio_set_solver_description) (IREF linsys, IREF procs, IREF info);
void FC_FUNC_(eio_get_solver_description,eio_get_solver_description) (IREF linsys, IREF procs, IREF info);
void FC_FUNC_(eio_set_solver,eio_set_solver) (IREF equation,
			     IREF main_type,
			     IREF sub_type,
			     IREF precond_type,
			     IREF stabilization,
			     IREF max_iter,
			     DREF stop_tol,
			     DREF steady_stop_tol,
			     IREF linearization,
			     IREF lin_max_iter,
			     DREF lin_stop_tol,
			     IREF lin_use_picard,
			     IREF lin_use_newton,
			     IREF newton_after_iter,
			     DREF newton_after_tol,
			     IREF info);
void FC_FUNC_(eio_get_solver,eio_get_solver) (IREF equation,
			     IREF main_type,
			     IREF sub_type,
			     IREF precond_type,
			     IREF stabilization,
			     IREF max_iter,
			     DREF stop_tol,
			     DREF steady_stop_tol,
			     IREF linearization,
			     IREF lin_max_iter,
			     DREF lin_stop_tol,
			     IREF lin_use_picard,
			     IREF lin_use_newton,
			     IREF newton_after_iter,
			     DREF newton_after_tol,
			     IREF info);
void FC_FUNC_(eio_set_timestep_head,eio_set_timestep_head) (IREF dependence, IREF len, IREF info);
void FC_FUNC_(eio_get_timestep_head,eio_get_timestep_head) (IREF dependence, IREF len, IREF info);
void FC_FUNC_(eio_set_timestep_field,eio_set_timestep_field) (IREF type,
			       int *nof_timesteps,
			       double *timestep_sizes,
			       int *output_intervals,
			       IREF steady_max_iter,
			       IREF info);
void FC_FUNC_(eio_get_timestep_field,eio_get_timestep_field) (IREF type,
			       int *nof_timesteps,
			       double *timestep_sizes,
			       int *output_intervals,
			       IREF steady_max_iter,
			       IREF info);




#ifdef __cplusplus
}
#endif		/* __cplusplus */


#endif

















