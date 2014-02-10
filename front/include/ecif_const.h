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
Module:     ecif_const.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Common constant values for all modules.

************************************************************************/

#ifndef _ECIF_CONST_
#define _ECIF_CONST_

#include "ecif_def.h"

#define MAX_DIMENSION 3
#define MAXDEPTH 4
#define MAXNAME 256
#define MINDEPTH 0
#define X 0
#define Y 1
#define Z 2

/* Stream etc. formatting */
#define ACCURACY 15
#define WIDTH  8

// 'No Such Value (Double)' <--> a double value was not found
extern double NSVD;

// Max initial value for boundbox coord.
extern double MAX_RANGE;

// Min initial value for boundbox coord.
extern double MIN_RANGE;

// Zero tolerances for point etc. comparitions
extern double EPSILON;  // Generic
extern double POINT_EPSILON;  // For vertex-points comparision

// Gap tolerance for adjagency etc. comparitions.
extern double GAP_TOLERANCE;


#define MAX_NOF_ELEM_CODES 32
#define MAX_NOF_BNDR_ELEMS 6
#define MAX_NOF_BNDR_NODES 12
#define MAX_NOF_EDGE_ELEMS 12
#define MAX_NOF_EDGE_NODES 3
#define MAX_NOF_NODES 32

//Cadi ecf-file format version number
extern const int ECIF_VERSION_NBR;

//Environment variables to be tested when
//loading tcl-sripts
extern const char* ELMER_TCL_LIB;
extern const char* ELMER_ECIF_TCL_LIB;


extern const int MAXBODYNAME;
extern const int MAXLINE;           // Maximum input file line length

extern const char NL1;              // Enter+Lf
extern const char NL2;              // Not in use

extern const int BUFFER_LEN;        // common read-in buffer length

extern const int NO_INDEX;          // to indicate a nonexisting tableindex
extern const int UNSET_INDEX;       // to indicate a uninitialized tableindex

extern const int MAX_FILE_NAME_LEN;  // Max size for a path,file,procedure name.
extern const int TIME_STR_SIZE;     // Size for the model time string.

extern const int MAX_NOF_BODIES;    // Max nof bodies and colors
extern const int MAX_NOF_SOLVERS;   // Max nof solvers (different linear systems)

extern const int MAX_NOF_COLOR_LEVELS; //Maximum R or G or B levels (fex. 255)

extern const int MAX_NOF_FLAG_NAMES;    // Max nof different named flags in the model

extern const short ESF_INDENT_SIZE;    // indent size in Elemer Solver Format file

extern const double PI;
extern const double HALF_PI;
extern const double TWO_PI;

extern const char* MESH_DIRECTORY_NAME; // "MESHDIR" value
extern const double MESH_H_INIT_FACTOR; // Factor * AvgDim is initial mesh size

extern double MAX_NORMAL_TOLERANCE;   // max comparison tol. for normal-vectors (deg. 0-180)
extern double MAX_DISTANCE_TOLERANCE; // max comparison tol. for planes etc. (relative 0-1)
extern double NORMAL_TOLERANCE;       // comparison tol. for normal-vectors (deg. 0-180)
extern double DISTANCE_TOLERANCE;     // comparison tol. for planes etc. (relative 0-1)

extern int OBJECT_DISPLAY_LIST_BASE;  // display-list id for a geometric object is: _BASE + objectId

extern int MAX_NOF_SPLIT_COMBINE_INFOS; // max array size <--> max undo level for mesh boundary splitting/combining

// Character separators in the parameter data string
extern const short MAX_NOF_PARAMETER_DATA_STRINGS;    // Max nof separate ';' separated "string-data" slots
extern const short MAX_NOF_PARAMETER_VARIABLES; // Max nof indipendent variables in parameters slots
extern const short MAX_PARAMETER_FIELD_NAME_LENGTH; // Max length for a field/variable name
extern const int PARAMETER_BUFFER_LEN;
extern char PARAMETER_BUFFER[];              // Buffer for parameter string
extern const char PARAMETER_PART_SEP;        // Paret separator in a parameter set
extern const char PARAMETER_FIELD_SEP;       // Parameter field separator in a parameter set
extern const char PARAMETER_DATA_SEP;        // Data value separator in a parameter field
extern const char FILE_NAME_INDICATOR;             // File name indicator in a parameter field like =:file.txt
extern const char PROC_NAME_INDICATOR;             // Procedure name indicator in a parameter field like =.libray_file;function_name

extern const int GUI_UPDATE_INTERVAL;       // Nof iterations after Gui is given an update message
/*
//MV 07.08.96
//This declaration (with const)works in VC++ 4.0, but not in SGI compilers.
//So, it is transformed into a non-const definition
extern const enum colorIndices;
*/
// Color indexes (as Tk names).
enum colorIndices  { ef_nodefault = -1,
  ef_black, ef_blue, ef_DodgerBlue, ef_DeepSkyBlue,
  ef_cyan, ef_DarkGreen, ef_LimeGreen, ef_green,
  ef_yellow, ef_orange, ef_OrangeRed, ef_red,
  ef_magenta, ef_violet, ef_pink, ef_white
};
extern const colorIndices defaultColorIndices[];
extern const colorIndices DEFAULT_COLOR_INDEX;
extern const char* colorNames[];
extern const int colorValues[][4];

// Constants used in output-file.
extern const char* topologyNames[];
extern const char* geometryNames[];

extern const char LB[]; //Left-brace for grouping-start
extern const char RB[]; //Right-brace for grouping-end

extern char read_buffer[];  // Buffer variable for general use.

typedef int ecif_modelStatus;
// Model status constants
// If code > 0 ==> error in the named area
extern const ecif_modelStatus STATUS_OK;
extern const ecif_modelStatus BODY_EQUATION_MISSING;
extern const ecif_modelStatus BODY_MATERIAL_MISSING;

typedef unsigned short beStatus;
// Bodyelement status constants
extern const beStatus BE_NONE;
extern const beStatus BE_DEVIDED;
extern const beStatus BE_SWAPPED;
extern const beStatus BE_OUTER;
extern const beStatus BE_INCLUDES_OUTER;
extern const beStatus BE_OUTER_CHECKED;
extern const beStatus BE_INNER_CANDIDATE;
extern const beStatus BE_INNER;
extern const beStatus BE_INNER_CHECKED;

extern const int DESC_ELEM_TYPE;            // Pos. of the elem-type field in elem description
extern const int DESC_NOF_NODES;            // Pos. of the nof nodes  field in elem description
extern const int DESC_NOF_MATCH_NODES;      // Pos. of the nof nodes neede to identify trhe element (corner nodes!)
extern const int DESC_HAS_INNER_NODE;       // Pos. of the inner node flag  field in elem description
extern const int DESC_NOF_BNDR_ELEMS;       // Pos. of the nof boundary elements field in elem description
extern const int DESC_NOF_NEEDED_BNDR_ELEMS;// Pos. of the nof neededboundary elements for defining all nodes
extern const int DESC_IS_1D_BNDR_ELEM;      // Pos. of the 1D boundary element flag
extern const int DESC_IS_2D_BNDR_ELEM;      // Pos. of the 2D boundary element flag
extern const int DESC_IS_3D_BNDR_ELEM;      // Pos. of the 3D boundary element flag
extern const int DESC_NOF_EDGES;            // Pos. of the nof edges
extern const int DESC_EDGE_ELEM_CODE;       // Pos. of the edge code (MEC_202 etc.)

extern const int MeshElementDesc[][12];
extern const int MeshElementReversedNodeIndices[][MAX_NOF_NODES];
extern const meshElementCode MeshElementBndrCodes[][MAX_NOF_BNDR_ELEMS];
extern const int MeshElementBndrNodes[][MAX_NOF_BNDR_ELEMS][MAX_NOF_BNDR_NODES];
extern const int MeshElementEdgeNodes[][MAX_NOF_EDGE_ELEMS][MAX_NOF_EDGE_NODES];


// Indexing variable separators in field names
extern const char INDEX_PRE_SEPARATOR;
extern const char INDEX_POST_SEPARATOR;

// ===============
// Emf field names
// ===============

//--Header
extern const char* EMF_HEADER;
extern const char* EMF_CREATED;
extern const char* EMF_MODIFIED;
extern const char* EMF_HAS_USER_DEFINITIONS;
extern const char* EMF_ELMER_FRONT_VERSION;
extern const char* EMF_ELMER_FRONT_INPUT_VERSION;
extern const char* EMF_TIMESTAMP;
extern const char* EMF_MODEL_STATUS;
extern const char* EMF_MODEL_SOURCE_TYPE;
extern const char* EMF_CAD_SOURCE_FILE;
extern const char* EMF_MESH_SOURCE_FILE;
extern const char* EMF_MESH_RESULT_FILE;
extern const char* EMF_MODEL_NAME;
extern const char* EMF_PROBLEM_NAME;
extern const char* EMF_MODEL_DESCRIPTION;
extern const char* EMF_PROBLEM_DESCRIPTION;
extern const char* EMF_MATC_FILE;
extern const char* EMF_MATC_FILE_EMF;
extern const char* EMF_MATC_FILE_SIF;
extern const char* EMF_INCLUDE_PATH;
extern const char* EMF_LOG_DIRECTORY;
extern const char* EMF_RESULTS_DIRECTORY;
extern const char* EMF_NOF_PROCESSORS;
extern const char* EMF_DIMENSION;
extern const char* EMF_MINIMUM_EDGE_SIZE;
extern const char* EMF_MESH_NAMES;
extern const char* EMF_CURRENT_MESH_INDEX;
extern const char* EMF_MESH_BG_MESH_FILE_INDICES;
extern const char* EMF_MESH_BG_MESH_FILES;
extern const char* EMF_MESH_BG_MESH_ACTIVES;
extern const char* EMF_MESH_BG_MESH_CONTROLS;

//--Timestamps
extern const char* EMF_TIMESTAMPS;
extern const char* EMF_TS_FRONT;
extern const char* EMF_TS_DATABASE;
extern const char* EMF_TS_GEBHARDT_FACTORS;
extern const char* EMF_TS_MESH;
extern const char* EMF_TS_SOLVER;
extern const char* EMF_TS_VIEW_FACTORS;

//--Statistics
extern const char* EMF_STATISTICS;
extern const char* EMF_NOF_BODIES;
extern const char* EMF_NOF_LOOPS;
extern const char* EMF_NOF_ELEMENTS;
extern const char* EMF_NOF_OUTER_BOUNDARIES;
extern const char* EMF_NOF_INNER_BOUNDARIES;
extern const char* EMF_NOF_VERTICES;
extern const char* EMF_MAX_LOOP_COUNT;

//--Geometry
extern const char* EMF_BODY;
extern const char* EMF_BODY1;
extern const char* EMF_BODY2;
extern const char* EMF_LAYER;
extern const char* EMF_LAYER_TAG;
extern const char* EMF_LAYER_TYPE;
extern const char* EMF_LAYER_NAME;
extern const char* EMF_LAYER_COLOR;
extern const char* EMF_CENTER;
extern const char* EMF_DEFINING_POINTS;
extern const char* EMF_DELTA_H;
extern const char* EMF_DELTA_N;
extern const char* EMF_DELTA_U;
extern const char* EMF_DIRECTION;
extern const char* EMF_EDGE;
extern const char* EMF_EDGES;
extern const char* EMF_EDGE_GROUP;
extern const char* EMF_EDGE_GROUPS;
extern const char* EMF_EDGE_LOOP;
extern const char* EMF_EDGE_LOOPS;
extern const char* EMF_ELEMENT;
extern const char* EMF_ELEMENTS;
extern const char* EMF_ELEMENT_GROUP;
extern const char* EMF_ELEMENT_GROUPS;
extern const char* EMF_ELEMENT_ID;
extern const char* EMF_ELEMENT_IDS;
extern const char* EMF_ELEMENT_TAG;
extern const char* EMF_ELEMENT_TAGS;
extern const char* EMF_ELEMENT_LOOP;
extern const char* EMF_ELEMENT_LOOPS;
extern const char* EMF_END_POINT;
extern const char* EMF_END_VERTEX;
extern const char* EMF_EXTRA_EDGES;
extern const char* EMF_EXTRA_VERTICES;
extern const char* EMF_FACE;
extern const char* EMF_FACES;
extern const char* EMF_FACE_GROUP;
extern const char* EMF_FACE_LOOP;
extern const char* EMF_FACE_LOOPS;
extern const char* EMF_GEOMETRY;
extern const char* EMF_GEOMETRY_SEGMENT;
extern const char* EMF_IDS_AND_POINTS;
extern const char* EMF_INNER_BOUNDARY;
extern const char* EMF_INNER_BOUNDARIES;
extern const char* EMF_MESH_H;
extern const char* EMF_MESH_F;
extern const char* EMF_MESH_N;
extern const char* EMF_MESH_U;
extern const char* EMF_OUTER_BOUNDARY;
extern const char* EMF_OUTER_BOUNDARIES;
extern const char* EMF_POINT;
extern const char* EMF_POINTS;
extern const char* EMF_POLYGON;
extern const char* EMF_RADIUS;
extern const char* EMF_START_POINT;
extern const char* EMF_START_VERTEX;
extern const char* EMF_USE_MESH_N;
extern const char* EMF_VERTEX;
extern const char* EMF_VERTEX_GROUP;
extern const char* EMF_VERTEX_TABLE;
extern const char* EMF_VERTICES;

//--General
extern const char* EMF_ACTIVE;
extern const char* EMF_ARGUMENT;
extern const char* EMF_ARGUMENTS;
extern const char* EMF_BOUNDARY_TAG;
extern const char* EMF_BODY_GEOMETRY_EDITED;
extern const char* EMF_BOUNDARY_GEOMETRY_EDITED;
extern const char* EMF_COLOR;
extern const char* EMF_DATA;
extern const char* EMF_FUNCTION;
extern const char* EMF_GROUP;
extern const char* EMF_ID_TABLE;
extern const char* EMF_INCLUDE;
extern const char* EMF_LIBRARY;
extern const char* EMF_NAME;
extern const char* EMF_OBJECT;
extern const char* EMF_OBJECT_TYPE;
extern const char* EMF_PARENT;
extern const char* EMF_SUB_PARENT;
extern const char* EMF_PARENT_TYPE;
extern const char* EMF_SIZE;
extern const char* EMF_TYPE;
extern const char* EMF_UNIT;

extern const char* EMF_INTEGER;
extern const char* EMF_LOGICAL;
extern const char* EMF_PROCEDURE;
extern const char* EMF_REAL;
extern const char* EMF_STRING;

//--Parameters
extern const char* EMF_BODY_FORCE;
extern const char* EMF_BODY_PARAMETER;
extern const char* EMF_BOUNDARY_CONDITION;
extern const char* EMF_BOUNDARY_PARAMETER;
extern const char* EMF_CALCULATOR;
extern const char* EMF_CONSTANT;
extern const char* EMF_COORDINATE;
extern const char* EMF_DATAFILE;
extern const char* EMF_EQUATION;
extern const char* EMF_EQUATION_VARIABLE;
extern const char* EMF_GRID_H;
extern const char* EMF_GRID_PARAMETER;
extern const char* EMF_INITIAL_CONDITION;
extern const char* EMF_MATERIAL;
extern const char* EMF_MODEL_PARAMETER;
extern const char* EMF_SIMULATION_PARAMETER;
extern const char* EMF_SOLVER;
extern const char* EMF_SOLVER_CONTROL;
extern const char* EMF_TIMESTEP;
extern const char* EMF_USER_SETTING;
//
extern const char* EMF_EXCLUDED_MESH_INDICES;
extern const char* EMF_INCLUDED_MESH_INDICES;
extern const char* EMF_MESH_INDICES;
extern const char* EMF_GRID_H_IDS;
extern const char* EMF_GRID_H_MESH_INDICES;
extern const char* EMF_GRID_PARAMETER_IDS;
extern const char* EMF_GRID_PARAMETER_MESH_INDICES;
extern const char* EMF_GRID_PARAMETER_PARENT;
extern const char* EMF_QUADGRID_MESH_INDICES;

//--User setting fields
extern const char* EMF_DEFAULT_MODEL_DIRECTORY;
extern const char* EMF_DEFAULT_MODEL_NAME;
extern const char* EMF_DEFAULT_CAD_FILES_DIRECTORY;
extern const char* EMF_DEFAULT_EXTERNAL_MESH_FILES_DIRECTORY;
extern const char* EMF_DEFAULT_INCLUDE_PATH;
extern const char* EMF_DEFAULT_RESULTS_DIRECTORY;
extern const char* EMF_DEFAULT_LOG_DIRECTORY;
extern const char* EMF_DEFAULT_USE_MODEL_SETTINGS;
extern const char* EMF_DEFAULT_AUTO_SAVE_EXTERNAL_MESH;
extern const char* EMF_AUTO_LOAD_DEFINITION_FILE;
extern const char* EMF_AUTO_LOAD_MESH;
extern const char* EMF_AUTO_SAVE_MODEL;
extern const char* EMF_AUTO_SAVE_SOLVER_INPUT;
extern const char* EMF_BROWSER_COMMAND;
extern const char* EMF_EDITOR_COMMAND;
extern const char* EMF_BROWSE_MODE_GEBHARDT_FACTORS;
extern const char* EMF_BROWSE_MODE_MESH;
extern const char* EMF_BROWSE_MODE_PROCEDURE_COMPILER;
extern const char* EMF_BROWSE_MODE_SOLVER;
extern const char* EMF_BROWSE_MODE_VIEW_FACTORS;
extern const char* EMF_FONT_SIZES;

// Sif names
// =========
// Section names
extern const char* SIF_HEADER;
extern const char* SIF_SIMULATION;
extern const char* SIF_CONSTANT;
extern const char* SIF_EQUATION;
extern const char* SIF_BODY;
extern const char* SIF_BOUNDARY;
extern const char* SIF_BOUNDARY_CONDITION;
extern const char* SIF_INITIAL_CONDITION;
extern const char* SIF_BODY_FORCE;
extern const char* SIF_MATERIAL;
extern const char* SIF_SOLVER;

// Generic field names
extern const char* SIF_ECHO_ON;
extern const char* SIF_CHECK_KEYWORDS;
extern const char* SIF_RELOAD_INPUT_FILE;


#endif
