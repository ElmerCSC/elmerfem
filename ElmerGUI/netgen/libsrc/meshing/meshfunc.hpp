#ifndef FILE_MESHFUNC
#define FILE_MESHFUNC

/**************************************************************************/
/* File:   meshfunc.hh                                                    */
/* Author: Johannes Gerstmayr                                             */
/* Date:   26. Jan. 98                                                    */
/**************************************************************************/


/*
  Functions for mesh-generations strategies
 */

class Mesh;
// class CSGeometry;

/// Build tet-mesh
MESHING3_RESULT MeshVolume(MeshingParameters & mp, Mesh& mesh3d);

/// Build mixed-element mesh
MESHING3_RESULT MeshMixedVolume(MeshingParameters & mp, Mesh& mesh3d);

/// Optimize tet-mesh
MESHING3_RESULT OptimizeVolume(MeshingParameters & mp, Mesh& mesh3d);
//			       const CSGeometry * geometry = NULL);

void RemoveIllegalElements (Mesh & mesh3d);


enum MESHING_STEP { 
  MESHCONST_ANALYSE = 1,
  MESHCONST_MESHEDGES = 2,
  MESHCONST_MESHSURFACE = 3,
  MESHCONST_OPTSURFACE = 4,
  MESHCONST_MESHVOLUME = 5,
  MESHCONST_OPTVOLUME = 6
};


#endif
