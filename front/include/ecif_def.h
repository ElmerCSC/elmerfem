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
Module:     ecif_def.h
Language:   C++
Date:       15.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   File includes all common definition for all C++ modules.

************************************************************************/

#ifndef _ECIF_DEF_
#define _ECIF_DEF_

// Debuggin flag (turns off try-catch block etc)
//
#undef FRONT_DEBUG
#define FRONT_DEBUG

// By default we use exceptions
#undef ELMER_FRONT_EXCEPTIONS
#if !defined(ELMER_FRONT_NOEXCEPTIONS)
  #define ELMER_FRONT_EXCEPTIONS
#endif

#undef UNIX
#if !defined(WIN32)
  #define UNIX
#endif

// Predefined names
extern char *progname;

// Elmer MATC do-math call
extern "C" char* mtc_domath(char* s);

// ***********************************
// Window-system dependent definitions
#if defined(WIN32)

  // For MS-WIN32
  #define Hdisplay int
  #define Hwindow  HWND
  #define Hinst HINSTANCE
  #define Hdll HINSTANCE
  #define Hfunc void*
  #define Hfile HANDLE
  #define Hprocess HANDLE
  #define ProcessId DWORD
  #define Callback CALLBACK
  #define Callbackp CALLBACK*
  #define LresCallback LRESULT CALLBACK
  #define Event UINT
  #define Wparam WPARAM
  #define Lpfunc FARPROC
  #define Lparam LPARAM
  #define Lpstr LPSTR
  #define Hdc HDC
  #define Hglrc HGLRC
  #define BIN binary
  #define NSTDint

#else

  // For Unix ( and X-Windows)
  #define Hdisplay Display
  #define Hwindow Window
  #define Hinst char*
  #define Hdll void*
  #define Hfunc void*
  #define Hfile int
  #define Hprocess int
  #define ProcessId unsigned int
  #define Callback
  #define Callbackp *
  #define LresCallback int
  #define Event XEvent
  #define Wparam void*
  #define Lpfunc int
  #define Lparam void*
  #define Lpstr char*
  #define Hdc GC
  #define Hglrc GLXContext
  #define BIN bin
  #define UINT  unsigned int
  #define UCHAR char
  #define MK_LBUTTON  Button1Mask
  #define MK_MBUTTON  Button2Mask
  #define MK_RBUTTON  Button3Mask
  #define MK_SHIFT    ShiftMask
  #define MK_CONTROL  ControlMask

#endif
// ***********************************



// ***********************************
// Enum defs
// ***********************************

// ======================
// Geometry related types
// ======================

// Geometry type
enum modelGeometryType {
  GEOM_CAD,
  GEOM_MESH,
  GEOM_CAD_AND_MESH,
  GEOM_CAD_OR_MESH
};

// Coordinate type
enum modelCoordinateType {
  COORD_CARTESIAN,
  COORD_CYLINDRIC,
  COORD_POLAR
};


// Body type
enum bodyType {
  NORMAL_BODY,
  BEM_BODY,
  VIRTUAL_BODY
};

// Body geometry type
enum bodyGmtrType {
  GEOM_BODY,
  MESH_BODY
};


// Body topology type
enum bodyTplgType {
  CLOSED_BODY,
  OPEN_BODY
};


// Body layer type
enum bodyLayerType {
  NONE_LAYER,
  EXPLICIT_LAYER,
  IMPLICIT_LAYER
};

// Body layer topology type
enum bodyLayerTplgType {
  CLOSED_LAYER,
  OPEN_LAYER
};


// Element loop type
enum elementLoopType {
  NORMAL_LOOP,
  LOGICAL_LOOP
};


// Element loop topology type
enum elementLoopTplgType {
  CLOSED_LOOP,
  OPEN_LOOP
};

// Element group type
enum elementGroupType {
  NONE_GROUP,
  EXPLICIT_GROUP,
  IMPLICIT_GROUP,
  VIRTUAL_GROUP
};


// Sub-element type
enum elementType {
  OUTER_ELEMENT,
  INNER_ELEMENT,
  INTERNAL_ELEMENT,
  BEM_ELEMENT,
  ANY_ELEMENT
};




// Object types for handlind Tcl ObjectTable types in cpp
// ***********************************
// ***********************************
//
// NOTE:
// --Do NOT change or delete these numbers!
// --Add new types always at the end!
//   This is because old input file (pre 7) parameters have parent object types
//   stored as (these) numbers!!!
//
// NOTE: Remember to add new type in functions:
//  objectType Model::objectName2Type(...)
//  const char* Model::objectType2Name(...)

enum objectType {
  OT_NONE =         0,
  OT_BODY =         1,
  OT_BODYPAIR =     2,
  OT_NOT_IN_USE  =  3,
  OT_ELEMENT_LOOP = 4,
  OT_BOUNDARY =     5,
  OT_FACE =         6,
  OT_EDGE =         7,
  OT_VERTEX =       8,
  OT_BODY_LAYER =   9,
  OT_ELEMENT_GROUP = 10,
};
//
// ***********************************
// ***********************************


// Geometry linearization factor type
enum linDeltaType {
  LIN_DELTA_NONE = 0,
  LIN_DELTA_H =    1,
  LIN_DELTA_N =    2,
  LIN_DELTA_U =    3
};


// Element by element intersection codes
enum matchType {
  MATCH_NONE,
  MATCH_1_INSIDE,
  MATCH_2_INSIDE,
  MATCH_EXACT,
  MATCH_OVERLAP
};

// Geometry error codes
enum geomError {
  GE_NONE =  0,
  GE_EMPTY =  1,
  GE_NOT_IN_ORDER = 2,
  GE_NOT_CLOSED = 3,
  GE_NOT_CONTINUOUS = 4,
};

// ========================
// Module and process types
// ========================

enum moduleType {
  CONTROL,
  GUI,
  RENDERER
};

enum frontProcessType {
  EGF_INPUT,
  EGF_OUTPUT,
  EMF_INPUT,
  EMF_OUTPUT,
  MESH_INPUT,
  MESH_OUTPUT
};


// Process priority level
enum priorityLevel {
  ECIF_NO_PRIORITY,
  ECIF_LOW_PRIORITY,
  ECIF_LOWER_THAN_NORMAL_PRIORITY,
  ECIF_NORMAL_PRIORITY,
  ECIF_HIGHER_THAN_NORMAL_PRIORITY,
  ECIF_HIGH_PRIORITY
};

// Miscellanous for info messages etc.
enum timeType {
  PROCESS_TIME,
  WALL_TIME
};

// =======================
// Rendering related types
// =======================

enum objectDrawingMode {
  DM_NORMAL,
  DM_WIREFRAME,
  DM_INTRA_LAYER,
  DM_HIDDEN
};

enum objectDrawingState {
  DS_NORMAL,
  DS_SELECTED
};

enum flagGroup {
  FLAG_GROUP_UNKNOWN = -1,
  GEOMETRY_TYPE = 0,
  GEOMETRY_EDITED,
  DRAW_SOURCE,
  DRAW_TARGET,
  SELECT_METHOD,
  SELECT_MODE,
  SELECT_OBJECTS,
  LABEL_DISPLAY
};

// codes for different flags which are
// set from/to Gui
enum flagName {
  FLAG_NAME_UNKNOWN = -1,
  // geometry type
  GEOMETRY_TYPE_CAD = 0,
  GEOMETRY_TYPE_MESH,
  // edited mode
  GEOMETRY_EDITED_BODIES,
  GEOMETRY_EDITED_BOUNDARIES,
  // geometry to be drawn
  DRAW_SOURCE_CAD,
  DRAW_SOURCE_MESH,
  // object to be drawn
  DRAW_TARGET_BODIES,
  DRAW_TARGET_SURFACES,
  DRAW_TARGET_EDGES,
  DRAW_TARGET_BOUNDARIES,
  DRAW_TARGET_BOUNDARY_EDGES,
  DRAW_TARGET_VERTICES,
  DRAW_TARGET_BOUNDARY_VERTICES,
  // element selection method
  SELECT_METHOD_SINGLE,
  SELECT_METHOD_ALL,
  SELECT_METHOD_BY_NEIGHBOR,
  SELECT_METHOD_BY_NORMAL,
  SELECT_METHOD_BY_PLANE,
  SELECT_METHOD_BY_RECTANGLE,
  SELECT_METHOD_BY_BOX,
  // element selection mode
  SELECT_MODE_TOGGLE,
  SELECT_MODE_EXTEND,
  SELECT_MODE_REDUCE,
  // object (boundaries, bodies) selection mode
  SELECT_OBJECTS_TOGGLE,
  SELECT_OBJECTS_EXTEND,
  // Label display control
  LABEL_DISPLAY_NODE,
  LABEL_DISPLAY_ELEMENT,
  LABEL_DISPLAY_VERTEX,
  LABEL_DISPLAY_EDGE,
  LABEL_DISPLAY_FACE,
  LABEL_DISPLAY_BODY,
};

enum renderingMode {
  RENDER_NORMAL,
  RENDER_SELECTION,
  RENDER_CALLBACK
};

enum mouseAction {
  MOUSE_DBL_CLCK,
  MOUSE_ALT_DBL_CLCK,
  MOUSE_CTRL_DBL_CLCK,
  MOUSE_SHFT_DBL_CLCK,
  MOUSE_CLCK,
  MOUSE_ALT_CLCK,
  MOUSE_CTRL_CLCK,
  MOUSE_SHFT_CLCK,
  MOUSE_MOVE,
  MOUSE_ALT_MOVE,
  MOUSE_CTRL_MOVE,
  MOUSE_SHFT_MOVE
};

enum keyAction {
  KEY_CTRL_B,
  KEY_CTRL_H,
  KEY_CTRL_L,
  KEY_CTRL_R,
  KEY_CTRL_X,
  KEY_CTRL_Y,
  KEY_CTRL_Z
};

// Renderer status codes
enum rendererStatus {
  HAS_NO_WINDOW,
  HAS_WINDOW,
  SHOW
};


// ==================
// Mesh related types
// ==================

// !!!===VERY IMPORTANT===!!!
// NOTE: Be careful that these ids match exactly the postions in mesh description
// tables in ecif_const.cpp, because they are used as indices
// into those tables
// NOTE: MEC_202S etc.are 'shell' elements, used as BEM bulk elements!
//
enum meshElementCode {
  MEC_000 = -1,
  MEC_101 = 0,

  // Lines
  MEC_202,
  MEC_203,

  // Triangles
  MEC_303,
  MEC_304,
  MEC_306,
  MEC_307,

  // Quadrilaterals
  MEC_404,
  MEC_405,
  MEC_408,
  MEC_409,

  // Tetras
  MEC_504,
  MEC_508,
  MEC_510,

  // Prism
  MEC_605,
  MEC_613,

  // Wedges
  MEC_706,
  MEC_715,

  // Bricks
  MEC_808,
  MEC_820,
  MEC_827
};


// Return codes for input etc.
enum returnCode {
	ECIF_OK,
	ECIF_STOP_READING,
	ECIF_ERROR,
	ECIF_MESH_FILE_NOT_FOUND,
	ECIF_MESH_ELEMENT_INDEX_ERROR,
	ECIF_MESH_NODE_INDEX_ERROR,
	ECIF_MESH_ERROR,
  ECIF_INCORRECT_MATERIAL_ID,
	ECIF_TOO_MANY_MATERIALS_PER_NODE
};
typedef returnCode Rc;


// ******************
// ENUM TYPES for trx
// ******************

// Model source (origin) type (CAD=cad file, MESH=mesh file, CAD_MESH=both
enum ecif_modelSource {
  ECIF_CAD_FILE,
  ECIF_MESH_FILE,
  ECIF_CAD_AND_MESH_FILE,
  ECIF_CAD_OR_MESH_FILE
};


// Model geometry dimension
enum ecif_modelDimension {
  ECIF_ND = 0,
  ECIF_1D = 1,
  ECIF_2D = 2,
  ECIF_3D = 3
};


// Model's symmetry axis position in screen coordinates
enum ecif_axisType {
  ECIF_NOAXIS,
  ECIF_YLEFT,
  ECIF_YRIGHT,
  ECIF_XBOT,
  ECIF_XTOP
};


// Geometry type of a boundary
enum ecif_geometryType {
  ECIF_NODIM,
  ECIF_GEOMETRY_MESH,
  ECIF_POINT,
  ECIF_LINE,
  ECIF_POLYLINE,
  ECIF_CIRCLE,
  ECIF_ELLIPSE,
  ECIF_PARABOLA,
  ECIF_HYPERBOLA,
  ECIF_SPLINE,
  ECIF_NURBS,
  ECIF_MULTI2D
};


// Geometry function type
enum ecif_functionType {
  ECIF_NOFUNCTION,
  ECIF_CPP,
  ECIF_F95,
  ECIF_MATC,
};

// Topology type of a boundary
enum ecif_topologyType {
  ECIF_NOFORM,
  ECIF_TOPOLOGY_MESH,
  ECIF_VERTEX,
  ECIF_EDGE,
  ECIF_FACE,
  ECIF_SHELL
};


// Type of parameter-data (initial cond., boundary cond. , material param. etc)
enum ecif_parameterType {
  ECIF_NOPARAM,
  ECIF_BODY_PARAMETER,
  ECIF_BOUNDARY_PARAMETER,
  ECIF_BODY_FORCE,
  ECIF_BOUNDARY_CONDITION,
  ECIF_CALCULATOR,
  ECIF_CONSTANT,
  ECIF_COORDINATE,
  ECIF_DATAFILE,
  ECIF_EQUATION,
  ECIF_EQUATION_VARIABLE,
  ECIF_GRID_PARAMETER,
  ECIF_GRID_H,
  ECIF_INITIAL_CONDITION,
  ECIF_MATERIAL,
  ECIF_MODEL_PARAMETER,
  ECIF_SIMULATION_PARAMETER,
  ECIF_SOLVER,
  ECIF_SOLVER_CONTROL,
  ECIF_TIMESTEP,
  ECIF_USER_SETTING
};



// *************************************
// Common types and common include files
// *************************************

// ================
// General typedefs
// ================
typedef float Color3f[3];
typedef int Color3[3];
typedef int Color4[4];
typedef float Color4f[4];
typedef double Point3[3];
typedef float Point3f[3];
typedef double Point4[4];
typedef float Point4f[4];


// By default we do not use the new stdlib streams
// (problems in Unix: Sgi,...)
  #include <iostream>
  #include <fstream>
  #include <strstream>
  #include <iomanip>

  #include <cctype>
  #include <cmath>
  #include <cstdio>
  #include <cstdlib>
  #include <cstring>
  #include <ctime>

  using namespace std;

#if 0

  #include <ctype.h>
  #include <math.h>
  #include <stdio.h>
  #include <stdlib.h>
  #include <string.h>
  #include <time.h>

#endif

#ifdef WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <windows.h>
  #undef WIN32_LEAN_AND_MEAN
  #undef NOMINMAX
#else
  #include <memory>
  #include <X11/Xlib.h>
#endif

// Define min/max macros if missing
// #ifndef max
// #define max(a,b)            (((a) > (b)) ? (a) : (b))
// #endif
// #ifndef min
// #define min(a,b)            (((a) < (b)) ? (a) : (b))
// #endif

#include <frontlib.h>
#include "ecif_const.h"
#include "ecif_func.h"

// ***********************************
// Misc typedefs and structs
// ***********************************

// Common declarations.
class Body;
class BodyLayer;
class BodyElement;
class BodyElementGroup;
class BodyElementLoop;
class BodyForce;
class BodyPair;
class BodyParameter;
class BoundaryCondition;
class BoundaryParameter;
class BoundBox;
class Control;
class Equation;
class EquationVariable;
class GcCircle;
class GcLine;
class GcMulti2D;
class GcPlane;
class GcPoint;
class GcPolyLine;
class GcNurbsCurve;
class GcNurbsSurface;
class Geometry;
class GridParameter;
class InitCondition;
class Input;
class Material;
class Model;
class ModelObject;
class Renderer;
class Parameter;
class ParameterField;
class Process;
class Solver;
class GcPoint;
class Timer;
class Timestep;
class UserInterface;
class UserSetting;

struct BoundaryPoint;
struct Compare_Ctring;
struct ParameterFieldInfo;
struct MeshCornerElement;
struct MeshData;
struct MeshInfo;
struct ModelData;
struct ModelInfo;
struct ModelStatistics;
struct ParallelInfo;
struct ParameterFieldInfo;
struct RendererInfo;

extern ostream* debugFile;
extern ostream* ecif_df; // debug-file!!



// Point in a paramtric form on a 3D surface;
struct ParamPair {
  double u;
  double v;
};

// Structure to parametric values for a body element;
struct ParamValues {
  int count; // number of parameter-pairs.
  ecif_topologyType t_type;
  ecif_geometryType g_type;
  ParamPair** values; // table of param-pairs.
};

// Color-type for OpenGL
typedef float RGBfloat;
const RGBfloat defaultColor[] = {0.0f, 0.0f, 0.0f, 1.0f};

// Used in bounds (bounding box, min-value-pairs etc)
// -- (min-x, max-x, min-y, max-y, min-z, max-z)
typedef double RangeVector[2 * MAX_DIMENSION];

// Used in bodyelement when checking outer-boundaries
typedef double ParamVector[4];


// Used for element loops
struct DirectedBodyElement {
  int direction;
  BodyElement* element;
};

// A helper structure. Needed when matching object boundaries.
struct AdjacentHalf {
  Body* body;
  int bodyLayer;
  BodyElement* element;
};

struct SplitCombineInfo {
  SplitCombineInfo();
  ~SplitCombineInfo();
  bool canRedo;
  bool canUndo;
  int splitSourceId;
  int splitTargetId;
  int nofSplitTargetMeshElements;
  int* splitTargetMeshElementSourceIndices;
  int combineTargetId;
  int nofCombineSourceIds;
  int* combineSourceIds;
};


// ====================================================
// Id-structures for Body-initial, Inner/Outer-boundary
// elements and their possible condition ids.
// ====================================================
struct Ids1 {
  int id1; //Body id
};

struct Ids2 {
  int id1; //Body id
  int id2; //Outer boundary element id
  int find2(int key1);
};

struct Ids3 {
  int id1; //Body1 id
  int id2; //Body2 id
  int id3; //Inner boundary element id
  int find3(int key1, int key2);
};


typedef int meshDictionaryEntry;

typedef int MeshEdge[2];


// This structure stores info on node's connection to other nodes
struct node2EdgeInfo {
  meshElementCode* edge_codes;  // Edge element codes
  int* edge_ids;                // Edge ids starting from this node
  int** edge_nodes;       // List of tables of nodes ids for each edge (except the first parent node!)
  short nof_conn;         // Current nof node's connections
  short max_nof_conn;     // Max nof of connections currenly allocated
};


// ============
// Miscellanous
// ============
struct ParameterFieldInfo {
  char* guiName;   // "HEAT_EQUATION"
  char* guiIndex;  // "oxygen"
  char* sifName;   // "Heat Equation"
  char* valueType; // "File", "Integer", "Logical", "Real", "String"
  bool alwaysOutput;
  bool isArray;
  bool isQuoted;
  bool isFileName;
  bool isPostIndexed;
  bool isPreIndexed;
  bool isProcName;
  bool sifTypeGiven;
  bool outputSif;
  bool outputSifType;

  ParameterFieldInfo() {
    guiName = NULL;
    guiIndex = NULL;
    sifName = NULL;
    valueType = NULL;
    alwaysOutput = true;
    isArray = false;
    isQuoted = false;
    isFileName = false;
    isPostIndexed = false;
    isPreIndexed = false;
    isProcName = false;
    sifTypeGiven = false;
    outputSif = true;
    outputSifType = true;
  }

  ~ParameterFieldInfo() {
    delete[] guiName;
    delete[] guiIndex;
    delete[] sifName;
    delete[] valueType;
  }
};


// A structure for controlling output to
// Solver input file (sif-file)
struct SifOutputControl {
  bool outputId;
  bool outputName;
  bool outputType;
  bool outputAll;
  bool reloadSolverIsOutput;
  const char* sectionName;

  SifOutputControl() {
    outputId = true;
    outputName = true;
    outputType = true;
    outputAll = false;
    sectionName = NULL;
    reloadSolverIsOutput = false;
  }
};

// A structure for storing information on model's units.
struct Units {

  Units() {
    conv[0] = conv[1] = conv[2] = 1.0;
    conv[3] = 1.0;
  }

  short code;
  char str[32];
  int descr;
  double conv[4];  // x,y,z and w scaling factors

};

#endif
