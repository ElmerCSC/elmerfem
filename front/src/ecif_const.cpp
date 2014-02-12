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
Module:     ecif_const.cpp
Language:   C++
Date:       01.11.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_def.h"

/*********************************************************************
Cadi ecf-file format version number (this can be used to control
model file conversion etc.

----------------------------------------------------------------------
001:-First release
002:-Boundary (linearizing) vertices added
003:
004:-Parameter masks removed from model file, vertices changed
005:-Boundary tags added, vertices changed, boundary vertices removed
006:-Equation parameter:
       AD_VARS added
       DIFFUSION_VARIABLES changed to AD_VAR

007: 1/2001
  -User defined equations and keywords added (wide changes!!!)
  -Enumerated field name constants and fieldName containers dropped,
     replaced with char constants
  -Equation parameter:
     AD_VAR dropped
     AD_VARS->ADVECTION_DIFFUSION_EQUATION_vars
     In general: *_vars added for each equation (like HEAT_EQUATION_vars)
  -Solver parameter:
     EQUATION_NAME->EQUATION, SOLVER_MESH->MESH

008: 12/2001
  -Blank fields are written into model file to block non-blank initial
     values to become in force (fields are now set first to init-values
     and after that replaced with model file value!)
  -'Size'-field option added to edf-file like (DIM*DIM) 2
  -'Always Output'-field option added (also to edf-file)
  -Errors in mesh-based geometry handling (model size etc.) fixed
  -Errors in system created vertex/edge handling fixed
  -Field help system added (this does not affect to the model file!)

009: 01/2003
  -New object types (OT_BODY_LAYER, OT_BOUNDARY_GROUP) added
  -Matc support
  -New sif format supported
----------------------------------------------------------------------

NOTE: From version 6 --> Field name updates are done in Tcl procs:
                         Interface::pre/postConverseInputData
**********************************************************************/

//*****************************
//
const int ECIF_VERSION_NBR = 9;
//
//*****************************

// 'No Such Value (Double)' <--> a double value was not found
double NSVD = 1.0e15;

// Max initial value for boundbox coord.
double MAX_RANGE = 1.0e15;

// Min initial value for boundbox coord.
double MIN_RANGE = -1.0e15;

// Default zero tolerance for point etc. comparitions
extern double EPSILON = 1.0e-08;
extern double POINT_EPSILON = 1.0e-08;

// Default gap tolerance for adjagency etc. comparitions.
extern double GAP_TOLERANCE = 1.0e-08;

//Environment variables to be tested when
//loading tcl-sripts
const char* ELMER_TCL_LIB = "ELMER_TCL_LIB";
const char* ELMER_ECIF_TCL_LIB = "ELMER_ECIF_TCL_LIB";

const int MAXBODYNAME = 20;
const int MAXLINE = 160 ;     // Maximum input file line length

#ifdef WIN32
  const char NL1 = '\n';      // Enter+Lf
  const char NL2 = '\0';      // Not in use
#else
  const char NL1 = '\n';      // Lf
  const char NL2 = '\r';      // Return
#endif

const int BUFFER_LEN = 1024;        // common read-in buffer length

const int NO_INDEX = -1;            // to indicate a nonexisting tableindex
const int UNSET_INDEX = -2;         // to indicate an uninitialized tableindex

const int MAX_FILE_NAME_LEN = 512;  // Max size for a path,file,procedure name.
const int TIME_STR_SIZE = 80;       // Size for the model time string.

const int MAX_NOF_BODIES = 999;     // Max nof bodies and colors (and mesh materials!)

const int MAX_NOF_SOLVERS = 12;     // Max nof solvers (different linear systems)


const int MAX_NOF_FLAG_NAMES =  256; // max nof different named flags in the model

const short ESF_INDENT_SIZE = 2; // indent size in Elemer Solver Format file


const double PI       = 3.1415926535;
const double HALF_PI  = PI / 2;
const double TWO_PI   = PI * 2;

const char* MESH_DIRECTORY_NAME = "MESHDIR";
const double MESH_H_INIT_FACTOR = 0.1; // Factor * AvgMin is initial mesh size


double MAX_NORMAL_TOLERANCE     = 180.0;   // max comparison tol. for normal-vectors (deg. 0-180)
double MAX_DISTANCE_TOLERANCE   = 1.0;  // max comparison tol. for planes etc. (relative 0-1)

// NOTE These are not constant, they are updated from Gui!
double NORMAL_TOLERANCE     = 0.5;    // comparison tol. for normal-vectors (deg. 0-180)
double DISTANCE_TOLERANCE   = 0.01;   // comparison tol. for planes etc. (relative 0-1)

int OBJECT_DISPLAY_LIST_BASE = 1000;  // display-list id for a geometric object is: _BASE + objectId

int MAX_NOF_SPLIT_COMBINE_INFOS = 16; // max array size <--> max undo level for mesh boundary splitting/combining

// Character separators in the parameter data string
const short MAX_NOF_PARAMETER_DATA_STRINGS = 1024;  // Max nof separate ';' separeated "string-data" slots
const short MAX_NOF_PARAMETER_VARIABLES = 4;      // Max nof indipendent variables in parameters slots
const short MAX_PARAMETER_FIELD_NAME_LENGTH = 81; // Max length for a field/variable name
const int PARAMETER_BUFFER_LEN = 2048;
char PARAMETER_BUFFER[1 + PARAMETER_BUFFER_LEN];    // Buffer for parameter string

const char PARAMETER_PART_SEP = '|';    // Part separator (|) in a parameter set
const char PARAMETER_FIELD_SEP = '|';   // Field separator (|) in a parameter set
const char PARAMETER_DATA_SEP = ';';    // Data separator (;) in a parameter field
const char FILE_NAME_INDICATOR = ':';   // File name indicator in a parameter field like =:file.txt
const char PROC_NAME_INDICATOR = '.';   // Procedure name indicator in a parameter field like =.libray_file;function_name


// Constants used in output-file.
// NOTE: Keep these congruent to the enum-definitions
// in the ecif_def.h and cadi_CPPdefs.h files !!!
const char* modeltypeNames[3] = {
  "NO_MODEL" "MODEL_2D", "MODEL_3D",
};
const char* topologyNames[4] = {
  "VERTEX", "EDGE", "FACE", "SHELL",
};
const char* geometryNames[4] = {
  "NODIM", "LINEAR", "CUBIC", "NURBS",
};


// Sif section keywords (names)
//
const char* SIF_HEADER = "Header";
const char* SIF_SIMULATION = "Simulation";
const char* SIF_CONSTANT = "Constants";
const char* SIF_EQUATION = "Equation";
const char* SIF_BODY = "Body";
const char* SIF_BOUNDARY = "Boundary";
const char* SIF_BOUNDARY_CONDITION = "Boundary Condition";
const char* SIF_INITIAL_CONDITION = "Initial Condition";
const char* SIF_BODY_FORCE = "Body Force";
const char* SIF_MATERIAL = "Material";
const char* SIF_SOLVER = "Solver";

const char LB[] = "{";
const char RB[] = "}";

// Model status constants (powers of 2)
// If code > 0 ==> error in the named area
const ecif_modelStatus STATUS_OK                = 0;
const ecif_modelStatus BODY_EQUATION_MISSING    = 1;
const ecif_modelStatus BODY_MATERIAL_MISSING    = 2;
const ecif_modelStatus MATERIAL_DENSITY_MISSING = 4;


// Bodyelement status constansts (powers of 2)
const beStatus BE_NONE            = 0;
const beStatus BE_DEVIDED         = 1;
const beStatus BE_SWAPPED         = 2;
const beStatus BE_INCLUDES_OUTER  = 4;
const beStatus BE_OUTER           = 8;
const beStatus BE_OUTER_CHECKED   = 16;
const beStatus BE_INNER_CANDIDATE = 32;
const beStatus BE_INNER           = 64;
const beStatus BE_INNER_CHECKED   = 128;


const int GUI_UPDATE_INTERVAL = 1000; // Nof iterations after Gui is given an update message

const int MAX_NOF_COLOR_LEVELS = 255; //maximum R or G or B levels

// The order of default colors. To be used if they
// are not defined in the Cad-input file.
const colorIndices defaultColorIndices[MAX_NOF_BODIES] = {
  ef_blue, ef_cyan, ef_green,
  ef_yellow, ef_orange, ef_red, ef_magenta,
  ef_violet, ef_pink, ef_DodgerBlue, ef_DarkGreen,
  ef_OrangeRed, ef_DeepSkyBlue, ef_LimeGreen, ef_white, ef_black,
};

const colorIndices DEFAULT_COLOR_INDEX = ef_black;

// Tcl/Tk names (Xlib)
const char* colorNames[MAX_NOF_BODIES] = {
  "black", "blue", "DodgerBlue", "DeepSkyBlue",
  "cyan", "DarkGreen", "LimeGreen", "green",
  "yellow", "orange", "OrangeRed", "red",
  "magenta", "violet", "pink", "white",
};

// Define a few colors. In OpenGL, RGB values are specified as
// floating-point numbers between 0.0 and 1.0.
// RGB values (and original I-deas names)
const int colorValues[MAX_NOF_BODIES][4] = {
  {0,   0,   0,   255}, // black
  {0,   0,   255, 255}, // blue
  {0,   84,  255, 255}, // gray_blue
  {0,   168, 255, 255}, // light_blue
  {0,   255, 255, 255}, // cyan
  {0,   84,  0,   255}, // dark_olive
  {0,   168, 0,   255}, // dark_green
  {0,   255, 0,   255}, // green
  {255, 255, 0,   255}, //  yellow
  {255, 168, 0,   255}, // golden_orange
  {255, 84,  0,   255}, // orange
  {255, 0,   0,   255}, // red
  {255, 0,   255, 255}, // magenta
  {255, 84,  255, 255}, // light_magenta
  {255, 168, 255, 255}, // pink
  {255, 255, 255, 255}, // white
};

// Buffers
char read_buffer[BUFFER_LEN+1]; // Buffer variable for general use.

const int DESC_ELEM_TYPE              = 0;  // Pos. of the elment type field in elem-desc
const int DESC_NOF_NODES              = 1;  // Pos. of the nof nodes field in elem-desc
const int DESC_NOF_MATCH_NODES        = 2;  // Pos. of the nof nodes needed to identify the elemetns (the corner nodes!)
const int DESC_HAS_INNER_NODE         = 3;  // Pos. of the inner node flag field in elem-desc
const int DESC_NOF_BNDR_ELEMS         = 4;  // Pos. of the nof boundary elemments field in elem-desc
const int DESC_NOF_NEEDED_BNDR_ELEMS  = 5;  // Pos. of the nof needed boundary elements for defining all nodes
const int DESC_IS_1D_BNDR_ELEM        = 6;  // Pos. of the bndr elem flag in 1D
const int DESC_IS_2D_BNDR_ELEM        = 7;  // Pos. of the bndr elem flag in 2D
const int DESC_IS_3D_BNDR_ELEM        = 8;  // Pos. of the bndr elem flag in 3D
const int DESC_NOF_EDGES              = 9;  // Pos. of the nof edges
const int DESC_EDGE_ELEM_CODE         = 10; // Pos. of the edge code (MEC_202 etc.)


// NOTE: VERY IMPORTANT!!!
// ***********************
// Keep the following arries consintent with the number and order of the element types!!!
// ***********************

// NOTE: DESC_NOF_BNDR_ELEMS info is also used when we have to match a bulk's subelement with neigbor's
// subelements. Eg. for all tetras we need only three triangle "boundary" nodes to match terta's
// subelements (faces) with neighborin faces.
// NOTE: NOF_MATCH_NODES depends strongly on the fact that mid-edge nodes an internal nodes are stored
// after corner nodes in all elements!!!

// 1.  Elem type
//
// 2.  Nof nodes
// 3.  Nof match nodes
// 4.  Inner node flag
//
// 5.  Nof boundary elems
// 6.  Nof needed bndr elems (to get all bndr nodes!)
//
// 7.  is 1D bndr elem flag
// 8.  is 2D bndr elem flag
// 9.  is 3D bndr elem flag
//
// 10. Nof edges
// 11. Edge code
//
const int MeshElementDesc[][12] =  {
  {101,  1, 1, 0,   0, 1,   1, 0, 0,   0, MEC_000},   // Vertex
  {202,  2, 2, 0,   2, 2,   0, 1, 0,   2, MEC_101},   // Line 202 (2 node)
  {203,  3, 2, 1,   2, 2,   0, 1, 0,   2, MEC_101},   // Line 203 (3 node)
  {303,  3, 3, 0,   3, 2,   0, 0, 1,   3, MEC_202},   // Triangle 303 (3 node)
  {304,  4, 3, 1,   3, 2,   0, 0, 1,   3, MEC_202},   // Triangle 304 (4 node, 1 in center)
  {306,  6, 3, 0,   3, 3,   0, 0, 1,   3, MEC_203},   // Triangle 306 (6 node)
  {307,  7, 3, 1,   3, 3,   0, 0, 1,   3, MEC_203},   // Triangle 307 (7 node, 1 in center)
  {404,  4, 4, 0,   4, 2,   0, 0, 1,   4, MEC_202},   // Quadrilateral 404 (4 node)
  {405,  5, 4, 1,   4, 2,   0, 0, 1,   4, MEC_202},   // Quadrilateral 405 (4 node, 1 in center)
  {408,  8, 4, 0,   4, 4,   0, 0, 1,   4, MEC_203},   // Quadrilateral 408 (8 node)
  {409,  9, 4, 1,   4, 4,   0, 0, 1,   4, MEC_203},   // Quadrilateral 409 (9 node, 1 in center)
  {504,  4, 4, 0,   4, 2,   0, 0, 0,   6, MEC_202},   // Tetra 504 (4 node)
  {508,  8, 4, 0,   4, 4,   0, 0, 0,   6, MEC_202},   // Tetra 508 (8 node)
  {510, 10, 4, 0,   4, 3,   0, 0, 0,   6, MEC_203},   // Tetra 510 (10 node)
  {605,  5, 5, 0,   5, 2,   0, 0, 0,   8, MEC_202},   // Prism 605 (5 node)  // NOTE Boundaries: 4 tris, 1 quadri
  {613, 13, 5, 0,   5, 3,   0, 0, 0,   8, MEC_203},   // Prims 613 (13 node) // NOTE Boundaries: 4 tris, 1 quadri
  {706,  6, 6, 0,   5, 2,   0, 0, 0,   9, MEC_202},   // Wedge 706 (6 node)  // NOTE Boundaries: 2 tris, 3 quadris
  {715, 15, 6, 0,   5, 3,   0, 0, 0,   9, MEC_203},   // Wedge 715 (15 node) // NOTE Boundaries: 2 tris, 3 quadris
  {808,  8, 8, 0,   6, 2,   0, 0, 0,  12, MEC_202},   // Brick 808 (8 node)
  {820, 20, 8, 0,   6, 4,   0, 0, 0,  12, MEC_203},   // Brick 820 (20 node)
  {827, 27, 8, 1,   6, 6,   0, 0, 0,  12, MEC_203}    // Brick 827 (27 node, 1 in center)
};


const int MeshElementReversedNodeIndices[][MAX_NOF_NODES] = {
  {0},                  // Vertex
  {1,0},                // Line 202 (2 node)
  {1,0,2},              // Line 203 (3 node)
  {0,2,1},              // Triangle 303 (3 node)
  {0,2,1,3},            // Triangle 304 (4 node, 1 in center)
  {0,2,1,5,4,3},        // Triangle 306 (6 node)
  {0,2,1,5,4,3,6},      // Triangle 307 (7 node, 1 in center)
  {0,3,2,1},            // Quadrilateral 404 (4 node)
  {0,3,2,1,4},          // Quadrilateral 405 (4 node, 1 in center)
  {0,3,2,1,7,6,5,4},    // Quadrilateral 408 (8 node)
  {0,3,2,1,7,6,5,4,8},  // Quadrilateral 409 (9 node, 1 in center)
  // NOTE: Reverse for bulk elements is actually meaningless!
  {0,1,2,3},                                            // Tetra 504 (4 node)
  {0,1,2,3,4,5,6,7},                                    // Tetra 508 (8 node)
  {0,1,2,3,4,5,6,7,8,9},                                // Tetra 510 (10 node)
  {0,1,2,3,4},                                          // Prism 605 (5 node)
  {0,1,2,3,4,5,6,7,8,9,10,11,12},                       // Prism 613 (13 node)
  {0,1,2,3,4,5},                                        // Wedge 706 (6 node)
  {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14},                 // Wedge 715 (15 node)
  {0,1,2,3,4,5,6,7},                                    // Brick 808 (8 node)
  {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19},  // Brick 820 (20 node)
  {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,232,24,25,26} // Brick 827 (27 node, 1 in center)
};


//--Boundary element codes for bulk elements
const meshElementCode MeshElementBndrCodes[][MAX_NOF_BNDR_ELEMS] = {
  {MEC_000},   // Vertex
  {MEC_101, MEC_101},   // Line 202 (2 node)
  {MEC_101, MEC_101},   // Line 203 (3 node)
  {MEC_202, MEC_202, MEC_202},   // Triangle 303 (3 node)
  {MEC_202, MEC_202, MEC_202},   // Triangle 304 (4 node, 1 in center)
  {MEC_203, MEC_203, MEC_203},   // Triangle 306 (6 node)
  {MEC_203, MEC_203, MEC_203},   // Triangle 307 (7 node, 1 in center)
  {MEC_202, MEC_202, MEC_202, MEC_202},   // Quadri 404 (4 node)
  {MEC_202, MEC_202, MEC_202, MEC_202},   // Quadri 405 (4 node, 1 in center)
  {MEC_203, MEC_203, MEC_203, MEC_203},   // Quadri 408 (8 node)
  {MEC_203, MEC_203, MEC_203, MEC_203},   // Quadri 409 (9 node, 1 in center)
  {MEC_303, MEC_303, MEC_303, MEC_303},   // Tetra 504 (4 node)
  {MEC_304, MEC_304, MEC_304, MEC_304},   // Tetra 508 (8 node)
  {MEC_306, MEC_306, MEC_306, MEC_306},   // Tetra 510 (10 node)
  {MEC_404, MEC_303, MEC_303, MEC_303, MEC_303},   // Prism 605 (5 node)
  {MEC_408, MEC_306, MEC_306, MEC_306, MEC_306},   // Prism 613 (13 node)
  {MEC_303, MEC_404, MEC_404, MEC_404, MEC_303},   // Wedge 706 (6 node)
  {MEC_306, MEC_408, MEC_408, MEC_408, MEC_306},   // Wedge 715 (15 node)
  {MEC_404, MEC_404, MEC_404, MEC_404, MEC_404, MEC_404}, // Brick 808 (8 node)
  {MEC_408, MEC_408, MEC_408, MEC_408, MEC_408, MEC_408}, // Brick 820 (20 node)
  {MEC_409, MEC_409, MEC_409, MEC_409, MEC_409, MEC_409}  // Brick 827 (27 node, 1 in center)
};

// Boundary elements as parent element faces defined by parent element nodes
// Direct orientation (faces ccw from outside when using Elmerpost local numbering)
const int MeshElementBndrNodes[][MAX_NOF_BNDR_ELEMS][MAX_NOF_BNDR_NODES] = {
  { {0} },                                  // Vertex (101)
  { {0},{1} },                              // Line (202)
  { {0},{1}, },                             // Line (203)
  { {0,1},{1,2},{2,0} },                    // Triangle (303)
  { {0,1},{1,2},{2,0} },                    // Triangle (304)
  { {0,1,3},{1,2,4},{2,0,5} },              // Triangle (306)
  { {0,1,3},{1,2,4},{2,0,5} },              // Triangle (307)
  { {0,1},{1,2},{2,3},{3,0} },              // Quadrilateral (404)
  { {0,1},{1,2},{2,3},{3,0} },              // Quadrilateral (405)
  { {0,1,4},{1,2,5},{2,3,6},{3,0,7} },      // Quadrilateral (408)
  { {0,1,4},{1,2,5},{2,3,6},{3,0,7} },      // Quadrilateral (409)
  { {0,2,1},{2,3,1},{0,3,2},{0,1,3} },      // Tetra (504)
  { {0,2,1,4},{2,3,1,6},{0,3,2,7},{0,1,3,5} }, // Tetra (508)
  { {0,2,1, 6,5,4},{2,3,1, 9,8,5},{0,3,2, 7,9,6},{0,1,3, 4,8,7} }, // Tetra (510)
  { {0,3,2,1},{0,1,4},{1,2,4},{2,3,4},{0,4,3} }, // Prism (605)
  { {0,3,2,1, 8,7,6,5},
    {0,1,4, 5,10,9},
    {1,2,4, 6,11,10},
    {2,3,4, 7,12,11},
    {0,4,3, 9,12,8} }, // Prism (613)
  { {0,1,2},{1,4,5,2},{0,2,5,3},{0,3,4,1},{3,5,4} }, // Wedge (706)
  { {0,1,2, 6,7,8},
    {1,4,5,2, 12,10,13,7},
    {0,2,5,3, 8,13,11,14},
    {0,3,4,1, 14,9,12,6},
    {3,5,4, 11,10,9} }, // Wedge (715)
  { {0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},{3,2,1,0},{4,5,6,7} }, // Brick (808)
  { {0,1,5,4,  8,13,16,12},
    {1,2,6,5,  9,14,17,13},
    {2,3,7,6, 10,15,18,14},
    {3,0,4,7, 11,12,19,15},
    {3,2,1,0, 10, 9, 8,11},
    {4,5,6,7, 16,17,18,19} }, // Brick (820)
  { {0,1,5,4,  8,13,16,12,20},
    {1,2,6,5,  9,14,17,13,21},
    {2,3,7,6, 10,15,18,14,22},
    {3,0,4,7, 11,12,19,15,23},
    {3,2,1,0, 10, 9, 8,11,24},
    {4,5,6,7, 16,17,18,19,25} } // Brick (827)
};


// Bulk element edges defined by bulk element nodes
// NOTE: No orientation, just all edges, start from smallest node id!
// (NOTE: Possible middle nodes always last in the list)
//
const int MeshElementEdgeNodes[][MAX_NOF_EDGE_ELEMS][MAX_NOF_EDGE_NODES] = {
  { {0} },                                  // Vertex (101)
  { {0},{1} },                              // Line (202)
  { {0},{1},},                              // Line (203)
  { {0,1},{0,2},{1,2} },                    // Triangle (303)
  { {0,1},{0,2},{1,2} },                    // Triangle (304)
  { {0,1,3},{0,2,5},{1,2,4} },              // Triangle (306)
  { {0,1,3},{0,2,5},{1,2,4} },              // Triangle (307)
  { {0,1},{0,3},{1,2},{2,3} },              // Quadrilateral (404)
  { {0,1},{0,3},{1,2},{2,3} },              // Quadrilateral (405)
  { {0,1,4},{0,3,7},{1,2,5},{2,3,6} },      // Quadrilateral (408)
  { {0,1,4},{0,3,7},{1,2,5},{2,3,6} },      // Quadrilateral (409)
  { {0,1},{0,2},{0,3},{1,2},{1,3},{2,3} },      // Tetra (504)
  { {0,1},{0,2},{0,3},{1,2},{1,3},{2,3} },      // Tetra (508)
  { {0,1,4},{0,2,6},{0,3,7},{1,2,5},{1,3,8},{2,3,9} },      // Tetra (510)
  { {0,1},{0,3},{0,4},{1,2},{1,4},{2,3},{2,4},{3,4} }, // Prism (605)
  { {0,1,5},{0,3,8},{0,4,9},{1,2,6},{1,4,10},{2,3,7},{2,4,11},{3,4,12} }, // Prism (613)
  { {0,1},{0,2},{0,3},{1,2},{1,4},{2,5},{3,4},{3,5},{4,5} }, // Wedge (706)
  { {0,1,6},{0,2,8},{0,3,14},{1,2,7},{1,4,12},{2,3,13},{3,4,9},{3,5,11},{4,5,10} }, // Wedge (715)
  { {0,1},{0,3},{0,4},{1,2},{1,5},{2,3},{2,6},{3,7},{4,5},{4,7},{5,6},{6,7} }, // Brick (808)
  { {0,1,8},{0,3,11},{0,4,12},{1,2,9},{1,5,13},{2,3,10},{2,6,14},{3,7,15},
    {4,5,16},{4,7,19},{5,6,17},{6,7,18} }, // Brick (820)
  { {0,1,8},{0,3,11},{0,4,12},{1,2,9},{1,5,13},{2,3,10},{2,6,14},{3,7,15},
    {4,5,16},{4,7,19},{5,6,17},{6,7,18} }, // Brick (827)
};



// Indexing variable separators in field names
const char INDEX_PRE_SEPARATOR = '<';
const char INDEX_POST_SEPARATOR = '>';

// ===============
// Emf field names
// ===============

//--Header
const char* EMF_HEADER = "Header";
const char* EMF_CREATED = "Created";
const char* EMF_MODIFIED = "Modified";
const char* EMF_HAS_USER_DEFINITIONS = "Has Definitions";
const char* EMF_ELMER_FRONT_VERSION = "Elmer Front Version";
const char* EMF_ELMER_FRONT_INPUT_VERSION = "Elmer Front Input Version";
const char* EMF_TIMESTAMP = "Timestamp";
const char* EMF_MODEL_STATUS = "Model Status";
const char* EMF_MODEL_SOURCE_TYPE = "Model Source Type";
const char* EMF_CAD_SOURCE_FILE = "Cad Source File";
const char* EMF_MESH_SOURCE_FILE = "Mesh Source File";
const char* EMF_MESH_RESULT_FILE = "Mesh Result File";
const char* EMF_MODEL_NAME = "Model Name";
const char* EMF_PROBLEM_NAME = "Problem Name";
const char* EMF_MODEL_DESCRIPTION = "Model Description";
const char* EMF_PROBLEM_DESCRIPTION = "Problem Description";
const char* EMF_MATC_FILE = "Matc File";
const char* EMF_MATC_FILE_EMF = "Matc File Emf";
const char* EMF_MATC_FILE_SIF = "Matc File Sif";
const char* EMF_INCLUDE_PATH = "Include Path";
const char* EMF_LOG_DIRECTORY = "Log Directory";
const char* EMF_RESULTS_DIRECTORY = "Results Directory";
const char* EMF_NOF_PROCESSORS = "Nof Processors";
const char* EMF_DIMENSION = "Dimension";
const char* EMF_MINIMUM_EDGE_SIZE = "Minimum Edge Size";
const char* EMF_MESH_NAMES = "Mesh Names";
const char* EMF_CURRENT_MESH_INDEX = "Current Mesh Index";
const char* EMF_MESH_BG_MESH_FILE_INDICES = "Mesh Bg Mesh File Indices";
const char* EMF_MESH_BG_MESH_FILES = "Mesh Bg Mesh Files";
const char* EMF_MESH_BG_MESH_ACTIVES = "Mesh Bg Mesh Actives";
const char* EMF_MESH_BG_MESH_CONTROLS = "Mesh Bg Mesh Controls";

//--Timestamps
const char* EMF_TIMESTAMPS = "Timestamps";
const char* EMF_TS_FRONT = "Front";
const char* EMF_TS_DATABASE = "Database";
const char* EMF_TS_GEBHARDT_FACTORS = "Gebhardt Factors";
const char* EMF_TS_MESH = "Mesh";
const char* EMF_TS_SOLVER = "Solver";
const char* EMF_TS_VIEW_FACTORS = "View Factors";

//--Statistics
const char* EMF_STATISTICS = "Statistics";
const char* EMF_NOF_BODIES = "Nof Bodies";
const char* EMF_NOF_LOOPS = "Nof Loops";
const char* EMF_NOF_ELEMENTS = "Nof Elements";
const char* EMF_NOF_OUTER_BOUNDARIES = "Nof Outer Boundaries";
const char* EMF_NOF_INNER_BOUNDARIES = "Nof Inner Boundaries";
const char* EMF_NOF_VERTICES = "Nof Vertices";
const char* EMF_MAX_LOOP_COUNT = "Max Loop Count";

//--Geometry
const char* EMF_BODY = "Body";
const char* EMF_BODY1 = "Body1";
const char* EMF_BODY2 = "Body2";
const char* EMF_LAYER = "Layer";
const char* EMF_LAYER_TAG = "Layer Tag";
const char* EMF_LAYER_TYPE = "Layer Type";
const char* EMF_LAYER_NAME = "Layer Name";
const char* EMF_LAYER_COLOR = "Layer Color";
const char* EMF_BODY_PAIR = "Body Pair";
const char* EMF_CENTER = "Center";
const char* EMF_DEFINING_POINTS = "Defining Points";
const char* EMF_DELTA_H = "Delta H";
const char* EMF_DELTA_N = "Delta N";
const char* EMF_DELTA_U = "Delta U";
const char* EMF_DIRECTION = "Direction";
const char* EMF_EDGE = "Edge";
const char* EMF_EDGES = "Edges";
const char* EMF_EDGE_GROUP = "Edge Group";
const char* EMF_EDGE_GROUPS = "Edge Groups";
const char* EMF_EDGE_LOOP = "Edge Loop";
const char* EMF_EDGE_LOOPS = "Edge Loops";
const char* EMF_ELEMENT = "Element";
const char* EMF_ELEMENTS = "Elements";
const char* EMF_ELEMENT_GROUP = "Element Group";
const char* EMF_ELEMENT_GROUPS = "Element Groups";
const char* EMF_ELEMENT_ID = "Element Id";
const char* EMF_ELEMENT_IDS = "Element Ids";
const char* EMF_ELEMENT_TAG = "Element Tag";
const char* EMF_ELEMENT_TAGS = "Element Tags";
const char* EMF_ELEMENT_LOOP = "Element Loop";
const char* EMF_ELEMENT_LOOPS = "Element Loops";
const char* EMF_END_POINT = "End Point";
const char* EMF_END_VERTEX = "End Vertex";
const char* EMF_EXTRA_EDGES = "Extra Edges";
const char* EMF_EXTRA_VERTICES = "Extra Vertices";
const char* EMF_FACE = "Face";
const char* EMF_FACES = "Faces";
const char* EMF_FACE_GROUP = "Face Group";
const char* EMF_FACE_LOOP = "Face Loop";
const char* EMF_FACE_LOOPS = "Face Loops";
const char* EMF_GEOMETRY = "Geometry";
const char* EMF_GEOMETRY_SEGMENT = "Geometry Segment";
const char* EMF_IDS_AND_POINTS = "Ids And Points";
const char* EMF_INNER_BOUNDARY = "Inner Boundary";
const char* EMF_INNER_BOUNDARIES = "Inner Boundaries";
const char* EMF_MESH_H = "Mesh H";
const char* EMF_MESH_F = "Mesh F";
const char* EMF_MESH_N = "Mesh N";
const char* EMF_MESH_U = "Mesh U";
const char* EMF_OUTER_BOUNDARY = "Outer Boundary";
const char* EMF_OUTER_BOUNDARIES = "Outer Boundaries";
const char* EMF_POINT = "Point";
const char* EMF_POINTS = "Points";
const char* EMF_POLYGON = "Polygon";
const char* EMF_RADIUS = "Radius";
const char* EMF_START_POINT = "Start Point";
const char* EMF_START_VERTEX = "Start Vertex";
const char* EMF_USE_MESH_N = "Use Mesh N";
const char* EMF_VERTEX = "Vertex";
const char* EMF_VERTEX_GROUP = "Vertex Group";
const char* EMF_VERTEX_TABLE = "Vertex Table";
const char* EMF_VERTICES = "Vertices";

//--General
const char* EMF_ACTIVE = "Active";
const char* EMF_ARGUMENT = "Argument";
const char* EMF_ARGUMENTS = "Arguments";
const char* EMF_BOUNDARY_TAG = "Boundary Tag";
const char* EMF_BODY_GEOMETRY_EDITED = "Body Geometry Edited";
const char* EMF_BOUNDARY_GEOMETRY_EDITED = "Boundary Geometry Edited";
const char* EMF_COLOR = "Color";
const char* EMF_DATA = "Data";
const char* EMF_FUNCTION = "Function";
const char* EMF_GROUP = "Group";
const char* EMF_ID_TABLE = "Id Table";
const char* EMF_INCLUDE = "Include";
const char* EMF_LIBRARY = "Library";
const char* EMF_NAME = "Name";
const char* EMF_OBJECT = "Object";
const char* EMF_OBJECT_TYPE = "Object Type";
const char* EMF_PARENT = "Parent";
const char* EMF_SUB_PARENT = "Sub Parent";
const char* EMF_PARENT_TYPE = "Parent Type";
const char* EMF_SIZE = "Size";
const char* EMF_TYPE = "Type";
const char* EMF_UNIT = "Unit";
const char* EMF_FILE = "File";
const char* EMF_INTEGER = "Integer";
const char* EMF_LOGICAL = "Logical";
const char* EMF_PROCEDURE = "Procedure";
const char* EMF_REAL = "Real";
const char* EMF_STRING = "String";

//--Parameters
const char* EMF_BODY_FORCE = "Body Force";
const char* EMF_BODY_PARAMETER = "Body Parameter";
const char* EMF_BOUNDARY_CONDITION = "Boundary Condition";
const char* EMF_BOUNDARY_PARAMETER = "Boundary Parameter";
const char* EMF_CALCULATOR = "Calculator";
const char* EMF_CONSTANT = "Constant";
const char* EMF_COORDINATE = "Coordinate";
const char* EMF_DATAFILE = "Datafile";
const char* EMF_EQUATION = "Equation";
const char* EMF_EQUATION_VARIABLE = "Equation Variable";
const char* EMF_GRID_H = "Grid H";
const char* EMF_GRID_PARAMETER = "Grid Parameter";
const char* EMF_INITIAL_CONDITION = "Initial Condition";
const char* EMF_MATERIAL = "Material";
const char* EMF_MODEL_PARAMETER = "Model Parameter";
const char* EMF_SIMULATION_PARAMETER = "Simulation Parameter";
const char* EMF_SOLVER = "Solver";
const char* EMF_SOLVER_CONTROL = "SolverControl";
const char* EMF_TIMESTEP = "Timestep";
const char* EMF_USER_SETTING = "User Setting";
//
const char* EMF_EXCLUDED_MESH_INDICES = "Excluded Mesh Indices";
const char* EMF_INCLUDED_MESH_INDICES = "Included Mesh Indices";
const char* EMF_MESH_INDICES = "Mesh Indices";
const char* EMF_GRID_H_IDS = "Grid H Ids";
const char* EMF_GRID_H_MESH_INDICES = "Grid H Mesh Indices";
const char* EMF_GRID_PARAMETER_IDS = "Grid Parameter Ids";
const char* EMF_GRID_PARAMETER_MESH_INDICES = "Grid Parameter Mesh Indices";
const char* EMF_GRID_PARAMETER_PARENT = "Grid Parameter Parent";
const char* EMF_QUADGRID_MESH_INDICES = "QuadGrid Mesh Indices";


//--User setting fields
const char* EMF_DEFAULT_MODEL_DIRECTORY = "Default Model Directory";
const char* EMF_DEFAULT_MODEL_NAME = "Default Model Name";
const char* EMF_DEFAULT_CAD_FILES_DIRECTORY = "Default Cad Files Directory";
const char* EMF_DEFAULT_EXTERNAL_MESH_FILES_DIRECTORY = "Default External Mesh Files Directory";
const char* EMF_DEFAULT_INCLUDE_PATH = "Default Include Path";
const char* EMF_DEFAULT_RESULTS_DIRECTORY = "Default Results Directory";
const char* EMF_DEFAULT_LOG_DIRECTORY = "Default Log Directory";
const char* EMF_DEFAULT_USE_MODEL_SETTINGS = "Default Use Model Settings";
const char* EMF_DEFAULT_AUTO_SAVE_EXTERNAL_MESH = "Default Auto Save External Mesh";
const char* EMF_AUTO_LOAD_DEFINITION_FILE = "Auto Load Definition File";
const char* EMF_AUTO_LOAD_MESH = "Auto Load Mesh";
const char* EMF_AUTO_SAVE_MODEL = "Auto Save Model";
const char* EMF_AUTO_SAVE_SOLVER_INPUT = "Auto Save Solver Input";
const char* EMF_BROWSER_COMMAND = "Browser Command";
const char* EMF_EDITOR_COMMAND = "Editor Command";
const char* EMF_BROWSE_MODE_GEBHARDT_FACTORS = "Browse Mode Gebhardt Factors";
const char* EMF_BROWSE_MODE_MESH = "Browse Mode Mesh";
const char* EMF_BROWSE_MODE_PROCEDURE_COMPILER = "Browse Mode Procedure Compiler";
const char* EMF_BROWSE_MODE_SOLVER = "Browse Mode Solver";
const char* EMF_BROWSE_MODE_VIEW_FACTORS = "Browse Mode View Factors";
const char* EMF_FONT_SIZES = "Font Sizes";


// Sif file keywords
// =================
const char* SIF_ECHO_ON = "Echo On";
const char* SIF_CHECK_KEYWORDS = "Check Keywords";
const char* SIF_RELOAD_INPUT_FILE = "Reload Input File";
