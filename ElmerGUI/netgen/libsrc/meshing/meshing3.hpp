#ifndef FILE_MESHING3
#define FILE_MESHING3




enum MESHING3_RESULT
{
  MESHING3_OK = 0,
  MESHING3_GIVEUP = 1,
  MESHING3_NEGVOL = 2,
  MESHING3_OUTERSTEPSEXCEEDED = 3,
  MESHING3_TERMINATE = 4,
  MESHING3_BADSURFACEMESH = 5
};


/// 3d volume mesh generation
class Meshing3
{
  /// current state of front
  AdFront3 * adfront;
  /// 3d generation rules
  ARRAY<vnetrule*> rules;
  /// counts how often a rule is used
  ARRAY<int> ruleused, canuse, foundmap;
  /// describes, why a rule is not applied
  ARRAY<char*> problems;
  /// tolerance criterion
  double tolfak;
public:
  /// 
  Meshing3 (const string & rulefilename); 
  /// 
  Meshing3 (const char ** rulep);
  ///
  virtual ~Meshing3 ();
  
  ///
  void LoadRules (const char * filename, const char ** prules);
  ///
  MESHING3_RESULT GenerateMesh (Mesh & mesh, const MeshingParameters & mp);
  
  ///
  int ApplyRules (ARRAY<Point3d> & lpoints, ARRAY<int> & allowpoint,
		  ARRAY<MiniElement2d> & lfaces, INDEX lfacesplit,
		  INDEX_2_HASHTABLE<int> & connectedpairs,
		  ARRAY<Element> & elements,
		  ARRAY<INDEX> & delfaces, int tolerance, 
		  double sloppy, int rotind1,
		  float & retminerr);
  
  ///
  PointIndex AddPoint (const Point3d & p, PointIndex globind);
  ///
  void AddBoundaryElement (const Element2d & elem);
  ///
  void AddBoundaryElement (const MiniElement2d & elem);
  ///
  int AddConnectedPair (const INDEX_2 & pair);
  
  ///
  void BlockFill (Mesh & mesh, double gh);
  ///
  void BlockFillLocalH (Mesh & mesh, const MeshingParameters & mp);

  /// uses points of adfront, and puts new elements into mesh
  void Delaunay (Mesh & mesh, int domainnr, const MeshingParameters & mp);
  ///
  friend class PlotVolMesh;
  ///
  friend void TestRules ();
};




/// status of mesh generation
class MeshingStat3d
{
public:
  ///
  MeshingStat3d ();
  ///
  int cntsucc;
  ///
  int cnttrials;
  ///
  int cntelem;
  ///
  int nff;
  ///
  int qualclass;
  ///
  double vol0;
  ///
  double vol;
  ///
  double h;
  ///
  int problemindex;
};





/*
template <typename POINTARRAY, typename FACEARRAY>
extern int FindInnerPoint (POINTARRAY & grouppoints,
			   FACEARRAY & groupfaces,
			   Point3d & p);

*/





#endif










