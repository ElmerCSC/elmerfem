#ifndef TOPOLOGY
#define TOPOLOGY

/**************************************************************************/
/* File:   topology.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   27. Apr. 01                                                    */
/**************************************************************************/

/*
    Mesh topology
    (Elements, Faces, Edges, Vertices
*/


class MeshTopology
{
  const Mesh & mesh;
  bool buildedges;
  bool buildfaces;

  MoveableArray<INDEX_2> edge2vert;
  MoveableArray<INDEX_4> face2vert;
  MoveableArray<int[12]> edges;
  MoveableArray<int[6]> faces;
  MoveableArray<int[4]> surfedges;
  MoveableArray<int> segedges;
  MoveableArray<int> surffaces;
  MoveableArray<INDEX_2> surf2volelement;
  MoveableArray<int> face2surfel;
  TABLE<int,PointIndex::BASE> *vert2element;
  TABLE<int,PointIndex::BASE> *vert2surfelement;
  TABLE<int,PointIndex::BASE> *vert2segment;
  int timestamp;
public:
  int GetNSurfedges() const {return surfedges.Size();}

  MeshTopology (const Mesh & amesh);
  ~MeshTopology ();

  void SetBuildEdges (bool be)
  { buildedges = be; }
  void SetBuildFaces (bool bf)
  { buildfaces = bf; }

  bool HasEdges () const
  { return buildedges; }
  bool HasFaces () const
  { return buildfaces; }

  void Update();


  int GetNEdges () const
  { return edge2vert.Size(); }
  int GetNFaces () const
  { return face2vert.Size(); }

  static int GetNVertices (ELEMENT_TYPE et);
  static int GetNEdges (ELEMENT_TYPE et);
  static int GetNFaces (ELEMENT_TYPE et);

  static const Point3d * GetVertices (ELEMENT_TYPE et);
  inline static const ELEMENT_EDGE * GetEdges (ELEMENT_TYPE et);
  inline static const ELEMENT_FACE * GetFaces (ELEMENT_TYPE et);

  
  int GetSegmentEdge (int segnr) const { return abs(segedges[segnr-1]); }
  int GetSegmentEdgeOrientation (int segnr) const { return sgn(segedges[segnr-1]); }

  void GetSegmentEdge (int segnr, int & enr, int & orient) const
  {
    enr = abs(segedges.Get(segnr));
    orient = segedges.Get(segnr) > 0 ? 1 : -1;
  }

  void GetElementEdges (int elnr, ARRAY<int> & edges) const;
  void GetElementFaces (int elnr, ARRAY<int> & faces, bool withorientation = false) const;
  void GetElementEdgeOrientations (int elnr, ARRAY<int> & eorient) const;
  void GetElementFaceOrientations (int elnr, ARRAY<int> & forient) const;

  int GetElementEdges (int elnr, int * edges, int * orient) const;
  int GetElementFaces (int elnr, int * faces, int * orient) const;

  void GetFaceVertices (int fnr, ARRAY<int> & vertices) const;
  void GetFaceVertices (int fnr, int * vertices) const;
  void GetEdgeVertices (int fnr, int & v1, int & v2) const;
  void GetFaceEdges (int fnr, ARRAY<int> & edges, bool withorientation = false) const;

  ELEMENT_TYPE GetFaceType (int fnr) const;

  void GetSurfaceElementEdges (int elnr, ARRAY<int> & edges) const;
  int GetSurfaceElementFace (int elnr) const;
  void GetSurfaceElementEdgeOrientations (int elnr, ARRAY<int> & eorient) const;
  int GetSurfaceElementFaceOrientation (int elnr) const;

  int GetSurfaceElementEdges (int elnr, int * edges, int * orient) const;

  void GetSurface2VolumeElement (int selnr, int & elnr1, int & elnr2) const
  { 
    elnr1 = surf2volelement.Get(selnr)[0];
    elnr2 = surf2volelement.Get(selnr)[1];
  }

  int GetFace2SurfaceElement (int fnr) const { return face2surfel[fnr-1]; }
  
  void GetVertexElements (int vnr, ARRAY<int> & elements) const;
  FlatArray<int> GetVertexElements (int vnr) const;

  void GetVertexSurfaceElements( int vnr, ARRAY<int>& elements ) const;
  FlatArray<int> GetVertexSurfaceElements (int vnr) const;


  
  int GetVerticesEdge ( int v1, int v2) const;
  void GetSegmentVolumeElements ( int segnr, ARRAY<int> & surfels ) const;
};













const ELEMENT_EDGE * MeshTopology :: GetEdges (ELEMENT_TYPE et)
{
  static int segm_edges[1][2] =
    { { 1, 2 }};

  static int trig_edges[3][2] =
    { { 3, 1 },
      { 2, 3 },        
      { 1, 2 }};

  static int quad_edges[4][2] =
    { { 1, 2 },
      { 3, 4 },
      { 4, 1 },
      { 2, 3 }};


  static int tet_edges[6][2] =
    { { 4, 1 },
      { 4, 2 },
      { 4, 3 }, 
      { 1, 2 },
      { 1, 3 },
      { 2, 3 }};

  static int prism_edges[9][2] =
    { { 3, 1 },
      { 1, 2 },
      { 3, 2 },
      { 6, 4 },
      { 4, 5 },
      { 6, 5 },
      { 3, 6 },
      { 1, 4 },
      { 2, 5 }};

  static int pyramid_edges[8][2] =
    { { 1, 2 },
      { 2, 3 },
      { 1, 4 },
      { 4, 3 },
      { 1, 5 },
      { 2, 5 },
      { 3, 5 },
      { 4, 5 }};

  static int hex_edges[12][2] =
    {
      { 1, 2 },
      { 3, 4 },
      { 4, 1 },
      { 2, 3 },
      { 5, 6 },
      { 7, 8 },
      { 8, 5 },
      { 6, 7 },
      { 1, 5 },
      { 2, 6 },
      { 3, 7 },
      { 4, 8 },
    };

  switch (et)
    {
    case SEGMENT:
    case SEGMENT3:
      return segm_edges;

    case TRIG:
    case TRIG6:
      return trig_edges;

    case QUAD:
    case QUAD6:
    case QUAD8:
      return quad_edges;

    case TET:
    case TET10:
      return tet_edges;

    case PYRAMID:
      return pyramid_edges;

    case PRISM:
    case PRISM12:
      return prism_edges;

    case HEX:
      return hex_edges;
    default:
      cerr << "Ng_ME_GetEdges, illegal element type " << et << endl;
    }
   return 0;  
}


const ELEMENT_FACE * MeshTopology :: GetFaces (ELEMENT_TYPE et)
{
  static const int trig_faces[1][4] = 
    { { 1, 2, 3, 0 } };
  static const int quad_faces[1][4] = 
    { { 1, 2, 3, 4 } };

  static const int tet_faces[4][4] =
    { { 4, 2, 3, 0 },
      { 4, 3, 1, 0 },
      { 4, 1, 2, 0 },
      { 1, 3, 2, 0 } };
  
  static const int prism_faces[5][4] =
    {
      { 1, 3, 2, 0 },
      { 4, 5, 6, 0 },
      { 3, 1, 4, 6 },
      { 1, 2, 5, 4 },
      { 2, 3, 6, 5 } 
    };

  static const int pyramid_faces[5][4] =
    {
      { 1, 2, 5, 0 },
      { 2, 3, 5, 0 },
      { 3, 4, 5, 0 },
      { 4, 1, 5, 0 },
      { 1, 4, 3, 2 } 
    };

  static const int hex_faces[6][4] =
    {
      { 1, 4, 3, 2 },
      { 5, 6, 7, 8 },
      { 1, 2, 6, 5 },
      { 2, 3, 7, 6 },
      { 3, 4, 8, 7 },
      { 4, 1, 5, 8 }
    };


  
  switch (et)
    {
    case TRIG:
    case TRIG6:
      return trig_faces;

    case QUAD:
    case QUAD6:
    case QUAD8:
      return quad_faces;


    case TET:
    case TET10:
      return tet_faces;

    case PRISM:
    case PRISM12:
      return prism_faces;

    case PYRAMID:
      return pyramid_faces;

    case SEGMENT:
    case SEGMENT3:

    case HEX:
      return hex_faces;

    default:
      cerr << "Ng_ME_GetVertices, illegal element type " << et << endl;
    }
  return 0;
}




















#endif
