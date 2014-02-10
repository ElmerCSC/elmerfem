#ifndef FILE_STLTOPOLOGY
#define FILE_STLTOPOLOGY

/**************************************************************************/
/* File:   stltopology.hpp                                                */
/* Author: Joachim Schoeberl                                              */
/* Author2: Johannes Gerstmayr                                            */
/* Date:   26. Jul. 99                                                    */
/**************************************************************************/

/*
  The STLTopology contains topologic information as
  triangle->point, point->triangles, triangle->edge, 2-points->edge,...
*/


class STLGeometry;

#define STLBASE 1

class STLPointIndex
{
  int i;
public:
  STLPointIndex () { ; }
  STLPointIndex (int ai) : i(ai) { ; }
  STLPointIndex & operator= (const STLPointIndex & ai) { i = ai.i; return *this; }
  STLPointIndex & operator= (int ai) { i = ai; return *this; }
  operator int () const { return i; }
  STLPointIndex operator++ (int) { return i++; }
  STLPointIndex operator-- (int) { return i--; }
};



class STLTrigIndex
{
  int i;
public:
  STLTrigIndex () { ; }
  STLTrigIndex (int ai) : i(ai) { ; }
  STLTrigIndex & operator= (const STLTrigIndex & ai) { i = ai.i; return *this; }
  STLTrigIndex & operator= (int ai) { i = ai; return *this; }
  operator int () const { return i; }
  STLTrigIndex operator++ (int) { return i++; }
  STLTrigIndex operator-- (int) { return i--; }
};





// triangle structure for loading stl files
class STLReadTriangle
{
  Vec<3> normal;
  Point<3> pts[3];
public:
  STLReadTriangle (const Point<3> * apts, const Vec<3> & anormal);
  STLReadTriangle () {};
  const Point<3> & operator[] (int i) const { return pts[i]; }
  const Vec<3> & Normal() const { return normal; }
};



class STLTriangle
{
  // topology edges of triangle, edge[i] opposite to point[i]
  int topedges[3];
  // neighbour triangles, trig[i] opposite to point[i]
  int nbtrigs[2][3]; 
  // normalized stored normal vector ??
  Vec<3> normal;
  // point numbers of triangle
  int pts[3];
  // front-side and back-side domains
  int domains[2];


public:

  Box<3> box;
  Point<3> center;
  double rad;
  int facenum;

  struct 
  {
    unsigned int toperror : 1;
  } flags;




  STLTriangle (const int * apts);
  STLTriangle () {pts[0]=0;pts[1]=0;pts[2]=0;}

  int operator[] (int i) const { return pts[i]; }
  int & operator[] (int i) { return pts[i]; }

  int EdgeNum(int i) const { return topedges[(i-1)]; }
  int & EdgeNum(int i) { return topedges[(i-1)]; }

  int NBTrig (bool side, int i) const { return nbtrigs[side][i]; }
  int & NBTrig (bool side, int i) { return nbtrigs[side][i]; }

  
  int Domain (bool side) const { return domains[side]; }
  int & Domain (bool side) { return domains[side]; }



  // obsolete:
  int PNum(int i) const { return pts[(i-1)]; }
  int & PNum(int i) { return pts[(i-1)]; }
  int PNumMod(int i) const { return pts[(i-1)%3]; }
  int & PNumMod(int i)  { return pts[(i-1)%3]; }

  int EdgeNumMod(int i) const { return topedges[(i-1)%3]; }
  int & EdgeNumMod(int i)  { return topedges[(i-1)%3]; }

  int NBTrigNum(int i) const { return nbtrigs[0][(i-1)]; }
  int & NBTrigNum(int i) { return nbtrigs[0][(i-1)]; }
  int NBTrigNumMod(int i) const { return nbtrigs[0][(i-1)%3]; }
  int & NBTrigNumMod(int i)  { return nbtrigs[0][(i-1)%3]; }
  

  // consistently oriented neighbour:
  int IsNeighbourFrom(const STLTriangle& t) const;
  // opposite to consistently oriented neighbour:
  int IsWrongNeighbourFrom(const STLTriangle& t) const;

  ///Get the two points of neighbour-Triangles in orientation of this-Triangle
  void GetNeighbourPoints(const STLTriangle& t, int& p1, int& p2) const;
  int GetNeighbourPointsAndOpposite(const STLTriangle& t, int& p1, int& p2, int& po) const;



  // NON-normalized geometry - normal vector
  Vec<3> GeomNormal(const ARRAY<Point<3> >& ap) const;
  
  // Stored normal vector, normalized
  void SetNormal (const Vec<3> & n);
  const Vec<3> & Normal () const { return normal; }


  void ChangeOrientation(); 

  //project with a certain normal vector in plane
  void ProjectInPlain(const ARRAY<Point<3> >& ap, 
		      const Vec<3> & n, Point<3> & pp) const;
  //project with the triangle's normal vector in plane
  void ProjectInPlain(const ARRAY<Point<3> > & ap, Point<3> & pp) const;


  /*
    Project the point pp along the nproj into the plane of
    the triangle. The triangle normal is given by ntrig to 
    avoid numerical instabilities.
    The local coordinates lam are defined by

    pp(input) = P1 + lam1 v1 + lam2 v2 + lam3 n

    the result is
    
    pp(output) = P1 + lam1 v1 + lam2 v2
  */
  int ProjectInPlain (const ARRAY<Point<3> >& ap, 
		      const Vec<3> & nproj, 
		      Point<3> & pp, Vec<3> & lam) const;

  int PointInside(const ARRAY<Point<3> >& ap, const Point<3> & pp) const;

  //get nearest point on triangle and distance to it
  double GetNearestPoint(const ARRAY<Point<3> >& ap, 
			 Point<3> & p3d) const;

  double Area(const ARRAY<Point<3> >& ap) const;

  double MinHeight(const ARRAY<Point<3> >& ap) const;
  double MaxLength(const ARRAY<Point<3> >& ap) const; 
  //max length of a side of triangle

  int GetFaceNum() {return facenum;}
  void SetFaceNum(int i) {facenum = i;}

  int HasEdge(int p1, int p2) const;
};


/**
   Topology Edge:
   Useful unside a face.
   A edges sharing more than 2 faces: trigs are undefined 
 */
class STLTopEdge 
{
  int pts[2];  
  int trigs[2];  
  double cosangle;
  int status;  // excluded, confirmed, candidate, undefined
public:
  STLTopEdge ();
  STLTopEdge (int p1, int p2, int trig1, int trig2);

  int operator[] (int i) const { return pts[i]; }
  int & operator[] (int i) { return pts[i]; }


  int PNum(int i) const { return pts[(i-1)]; }
  int & PNum(int i) { return pts[(i-1)]; }
  int PNumMod(int i) const { return pts[(i-1)%2]; }
  int & PNumMod(int i)  { return pts[(i-1)%2]; }

  int TrigNum(int i) const { return trigs[(i-1)]; }
  int & TrigNum(int i) { return trigs[(i-1)]; }
  int TrigNumMod(int i) const { return trigs[(i-1)%2]; }
  int & TrigNumMod(int i)  { return trigs[(i-1)%2]; }

  void SetCosAngle (double ca) { cosangle = ca; }
  double CosAngle () const { return cosangle; }
  double Angle () const { return acos (cosangle); }

  void SetStatus (int stat) { status = stat; }
  int GetStatus () const { return status; }
};



ostream& operator<<(ostream& os, const STLTriangle& t);







class STLTopology
{
protected:
  ARRAY<STLTriangle> trias;
  ARRAY<STLTopEdge> topedges;
  ARRAY<Point<3> > points;

  // mapping of sorted pair of points to topedge
  INDEX_2_HASHTABLE<int> * ht_topedges;
  // mapping of node to trigs
  TABLE<int> trigsperpoint; 
  // mapping of node to edges
  TABLE<int> topedgesperpoint; 
  
  // searchtree for trigs and points

  Box3dTree * searchtree; // ADT
  Point3dTree * pointtree;

  Box<3> boundingbox;
  double pointtol;

public:
  enum STL_GEOM_STATUS { STL_GOOD, STL_WARNING, STL_ERROR };

protected:
  STL_GEOM_STATUS status;
  string statustext;
  
  bool topology_ok;
  bool orientation_ok;

public:
  STLTopology();
  virtual ~STLTopology();

  static STLGeometry * LoadNaomi (istream & ist);
  static STLGeometry * Load (istream & ist);
  static STLGeometry * LoadBinary (istream & ist);

  void Save (const char* filename);
  void SaveBinary (const char* filename, const char* aname);
  void SaveSTLE (const char * filename); // stores trigs and edges
  
  virtual void InitSTLGeometry (const ARRAY<STLReadTriangle> & readtrigs);

  virtual void TopologyChanged() {}; //do some things, if topology changed!

  /// Generate topology tables
  void FindNeighbourTrigs();

  
  void GetTrianglesInBox (const Box<3> & box,
			  ARRAY<int> & trias) const;


  int GetNP() const { return points.Size(); }
  int AddPoint(const Point<3> & p) { return points.Append(p); }
  const Point<3> & GetPoint(int nr) const { return points.Get(nr); }
  int GetPointNum (const Point<3> & p);
  void SetPoint(int nr, const Point<3> & p) { points.Elem(nr) = p; }
  const ARRAY<Point<3> >& GetPoints() const { return points; }

  const Point<3> & operator[] (STLPointIndex i) const { return points[i]; }
  Point<3> & operator[] (STLPointIndex i) { return points[i]; }




  int GetNT() const { return trias.Size(); }
  void AddTriangle(const STLTriangle& t);
  const STLTriangle & GetTriangle (int nr) const { return trias.Get(nr); }
  STLTriangle & GetTriangle (int nr) { return trias.Elem(nr); }
  
  const STLTriangle & operator[] (STLTrigIndex i) const { return trias[i]; }
  STLTriangle & operator[] (STLTrigIndex i) { return trias[i]; }


  int GetNTE() const { return topedges.Size(); }
  const STLTopEdge & GetTopEdge (int nr) const { return topedges.Get(nr); }
  STLTopEdge & GetTopEdge (int nr)  { return topedges.Elem(nr); }
  int GetTopEdgeNum (int pi1, int pi2) const;


  int NOTrigsPerPoint(int pn) { return trigsperpoint.EntrySize(pn); }
  int TrigPerPoint(int pn, int i) { return trigsperpoint.Get(pn, i); }


  int NTopEdgesPerPoint (int pn) const { return topedgesperpoint.EntrySize(pn); }
  int TopEdgePerPoint (int pn, int ei) const { return topedgesperpoint.Get(pn, ei); }

  
  bool Topology_Ok() const { return topology_ok; }
  bool Orientation_Ok() const { return orientation_ok; }

  STL_GEOM_STATUS GetStatus () const { return status; }
  const string & GetStatusText () const { return statustext; }

  void InvertTrig (int trig);
  void DeleteTrig (int trig);
  void OrientAfterTrig (int trig);


  // Table will be constructed, if topology is not ok
  /// neighbourtrigs for surfacetrigs
  TABLE<int> neighbourtrigs;

  /// get nr-th neighbour Triangle for triangle trig
  int NONeighbourTrigs(int trig) const { return neighbourtrigs.EntrySize(trig); }
  int NeighbourTrig(int trig, int nr) const { return neighbourtrigs.Get(trig,nr); }
  int NeighbourTrigSorted(int trig, int nr) const;
  void AddNeighbourTrig(int i, int nt) { neighbourtrigs.Add1(i, nt); }




  int GetLeftTrig (int p1, int p2) const;
  int GetRightTrig (int p1, int p2) const;

  const Box<3> & GetBoundingBox () const { return boundingbox; }
};


#endif
