#ifndef NGLIB
#define NGLIB

/**************************************************************************/
/* File:   nglib.hh                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   7. May. 2000                                                   */
/**************************************************************************/

/*
  
Interface to the netgen meshing kernel
  
*/

/// Data type for NETGEN mesh
typedef void * Ng_Mesh;

/// Data type for NETGEN CSG geomty
typedef void * Ng_CSG_Geometry;

/// Data type for NETGEN 2D geomty
typedef void * Ng_Geometry_2D;

/// Data type for NETGEN STL geomty
typedef void * Ng_STL_Geometry;



// max number of nodes per element
#define NG_VOLUME_ELEMENT_MAXPOINTS 10

// implemented element types:
enum Ng_Volume_Element_Type { NG_TET = 1, NG_PYRAMID = 2, NG_PRISM = 3,
                              NG_TET10 = 4 };

// max number of nodes per surface element
#define NG_SURFACE_ELEMENT_MAXPOINTS 6

// implemented element types:
enum Ng_Surface_Element_Type { NG_TRIG = 1, NG_QUAD = 2, 
			       NG_TRIG6 = 3 };



class Ng_Meshing_Parameters 
{
 public:
  
  double maxh;
  double fineness;   // 0 .. coarse, 1 .. fine
  int secondorder;
  char * meshsize_filename;
  int quad_dominated;

  Ng_Meshing_Parameters();
};


enum Ng_Result { NG_OK = 0, 
		 NG_SURFACE_INPUT_ERROR = 1,
		 NG_VOLUME_FAILURE = 2, 
		 NG_STL_INPUT_ERROR = 3,
		 NG_SURFACE_FAILURE = 4,
		 NG_FILE_NOT_FOUND = 5 };




  
// initialize, deconstruct Netgen library:
void Ng_Init ();
void Ng_Exit ();
  
// ===== Elmer related additions =====

// return surface and volume element in pi + underlying stl number
Ng_Surface_Element_Type 
EG_GetSurfaceElement (Ng_Mesh * mesh, int num, int * pi, int * ptrignum);

// return surface element bc property:
int
EG_GetSurfaceElementBCProperty(Ng_Mesh * mesh, int i);

// return segment bcnum
void
EG_GetSegmentBCProperty (Ng_Mesh *mesh, Ng_Geometry_2D *geom, int num, int * bcnum);


// Generates new mesh structure
Ng_Mesh * Ng_NewMesh ();
void Ng_DeleteMesh (Ng_Mesh * mesh);
  
// feeds points, surface elements and volume elements to the mesh
void Ng_AddPoint (Ng_Mesh * mesh, double * x);
void Ng_AddSurfaceElement (Ng_Mesh * mesh, Ng_Surface_Element_Type et,
                           int * pi);
void Ng_AddVolumeElement (Ng_Mesh * mesh, Ng_Volume_Element_Type et,
                          int * pi);
  
// ask for number of points, surface and volume elements
int Ng_GetNP (Ng_Mesh * mesh);
int Ng_GetNSE (Ng_Mesh * mesh);
int Ng_GetNE (Ng_Mesh * mesh);
  
//  return point coordinates
void Ng_GetPoint (Ng_Mesh * mesh, int num, double * x);

// return surface and volume element in pi
Ng_Surface_Element_Type 
Ng_GetSurfaceElement (Ng_Mesh * mesh, int num, int * pi);

Ng_Volume_Element_Type
Ng_GetVolumeElement (Ng_Mesh * mesh, int num, int * pi);


// Defines MeshSize Functions
void Ng_RestrictMeshSizeGlobal (Ng_Mesh * mesh, double h);
void Ng_RestrictMeshSizePoint (Ng_Mesh * mesh, double * p, double h);
void Ng_RestrictMeshSizeBox (Ng_Mesh * mesh, double * pmin, double * pmax, double h);
  
// generates volume mesh from surface mesh
Ng_Result Ng_GenerateVolumeMesh (Ng_Mesh * mesh, Ng_Meshing_Parameters * mp);

void Ng_SaveMesh(Ng_Mesh * mesh, const char* filename);
Ng_Mesh * Ng_LoadMesh(const char* filename);





// **********************************************************
// **   2D Meshing                                         **
// **********************************************************


// feeds points and boundary to mesh

void Ng_AddPoint_2D (Ng_Mesh * mesh, double * x);
void Ng_AddBoundarySeg_2D (Ng_Mesh * mesh, int pi1, int pi2);
  
// ask for number of points, elements and boundary segments
int Ng_GetNP_2D (Ng_Mesh * mesh);
int Ng_GetNE_2D (Ng_Mesh * mesh);
int Ng_GetNSeg_2D (Ng_Mesh * mesh);
  
//  return point coordinates
void Ng_GetPoint_2D (Ng_Mesh * mesh, int num, double * x);

// return 2d triangles
void Ng_GetElement_2D (Ng_Mesh * mesh, int num, int * pi, int * matnum = NULL);

// return 2d boundary segment
void Ng_GetSegment_2D (Ng_Mesh * mesh, int num, int * pi, int * matnum = NULL);


// load 2d netgen spline geometry
Ng_Geometry_2D * Ng_LoadGeometry_2D (const char * filename);

// generate 2d mesh, mesh is allocated by function
Ng_Result Ng_GenerateMesh_2D (Ng_Geometry_2D * geom,
                              Ng_Mesh ** mesh,
                              Ng_Meshing_Parameters * mp);
  
void Ng_HP_Refinement (Ng_Geometry_2D * geom,
                       Ng_Mesh * mesh,
                       int levels);
  




// **********************************************************
// **   STL Meshing                                        **
// **********************************************************


// loads geometry from STL file
Ng_STL_Geometry * Ng_STL_LoadGeometry (const char * filename, int binary = 0);


// generate new STL Geometry
Ng_STL_Geometry * Ng_STL_NewGeometry ();
  

// fills STL Geometry
// positive orientation
// normal vector may be null-pointer
void Ng_STL_AddTriangle (Ng_STL_Geometry * geom, 
                         double * p1, double * p2, double * p3, 
                         double * nv = NULL);

// add (optional) edges :
void Ng_STL_AddEdge (Ng_STL_Geometry * geom, 
                     double * p1, double * p2);

// after adding triangles (and edges) initialize
Ng_Result Ng_STL_InitSTLGeometry (Ng_STL_Geometry * geom);

// automatically generates edges:
Ng_Result Ng_STL_MakeEdges (Ng_STL_Geometry * geom,
                            Ng_Mesh* mesh,
                            Ng_Meshing_Parameters * mp);


// generates mesh, empty mesh must be already created.
Ng_Result Ng_STL_GenerateSurfaceMesh (Ng_STL_Geometry * geom,
                                      Ng_Mesh * mesh,
                                      Ng_Meshing_Parameters * mp);


#ifdef ACIS

// **********************************************************
// **   ACIS Meshing                                       **
// **********************************************************

/// Data type for NETGEN STL geomty
typedef void * Ng_ACIS_Geometry;

// loads geometry from STL file
Ng_ACIS_Geometry * Ng_ACIS_LoadGeometry (const char * filename);
  
// generates mesh, empty mesh must be already created.
Ng_Result Ng_ACIS_GenerateSurfaceMesh (Ng_ACIS_Geometry * geom,
                                       Ng_Mesh * mesh,
                                       Ng_Meshing_Parameters * mp);


#endif


#endif
