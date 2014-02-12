#ifndef FILE_RULER3
#define FILE_RULER3


/**
  3D element generation rule.
 */
class vnetrule
{
private:
  /// rule is applicable for quality classes above this value
  int quality;
  /// name of rule
  char * name;
  /// point coordinates in reference position
  ARRAY<Point3d> points;
  /// old and new faces in reference numbering
  ARRAY<Element2d> faces;
  /// additional edges of rule
  ARRAY<twoint> edges;

  /// points of freezone in reference coordinates
  ARRAY<Point3d> freezone;
  /// points of freezone in reference coordinates if tolcalss to infty
  ARRAY<Point3d> freezonelimit;
  /// point index, if point equal to mappoint, otherwise 0
  ARRAY<int> freezonepi;
  /// faces of each convex part of freezone
  ARRAY<ARRAY<threeint>*> freefaces;
  /// set of points of each convex part of freezone
  ARRAY<ARRAY<int>*> freesets;
  /// points of transformed freezone
  ARRAY<Point3d> transfreezone;
  /// edges of each convex part of freezone
  ARRAY<ARRAY<twoint>*> freeedges;

  /// face numbers to be deleted
  ARRAY<int> delfaces;
  /// elements to be generated
  ARRAY<Element> elements;
  /// tolerances for points and faces (used ??)
  ARRAY<double> tolerances, linetolerances;
  /// transformation matrix 
  DenseMatrix oldutonewu;
  /// transformation matrix: deviation old point to dev. freezone
  DenseMatrix * oldutofreezone;
  /** transformation matrix: deviation old point to dev. freezone, 
    quality class to infinity */
  DenseMatrix * oldutofreezonelimit;

  // can be deleted:
  // BaseMatrix *outf, *outfl;

  /**
    a point is outside of convex part of freezone, 
    iff mat * (point, 1) >= 0 for each component (correct ?)
    */
  ARRAY<DenseMatrix*> freefaceinequ;
  /// 
  ARRAY<fourint> orientations;
  /**
    flags specified in rule-description file:
    t .. test rule
    */
  ARRAY<char> flags;

  /**
    topological distance of face to base element
    non-connected: > 100  (??) 
    */
  ARRAY<int> fnearness;
  ARRAY<int> pnearness;
  int maxpnearness;

  /// number of old points in rule
  int noldp;
  /// number of new poitns in rule
  int noldf;
  /// box containing free-zone
public:  
  // double fzminx, fzmaxx, fzminy, fzmaxy, fzminz, fzmaxz;
  Box3d fzbox;

public:
  
  ///
  vnetrule ();
  ///
  ~vnetrule ();
  ///
  int GetNP () const { return points.Size(); }
  ///
  int GetNF () const { return faces.Size(); }
  ///
  int GetNE () const { return elements.Size(); }
  ///
  int GetNO () const { return orientations.Size(); }
  ///
  int GetNEd () const { return edges.Size(); }
  ///
  int GetNOldP () const { return noldp; }
  ///
  int GetNOldF () const { return noldf; }
  ///
  int GetNDelF () const { return delfaces.Size(); }
  ///
  int GetQuality () const { return quality; }
  ///
  int GetFNearness (int fi) const { return fnearness.Get(fi); }
  ///
  int GetPNearness (int pi) const { return pnearness.Get(pi); }
  ///
  int GetMaxPNearness () const { return maxpnearness; }


  ///
  const Point3d & GetPoint (int i) const { return points.Get(i); }
  ///
  const Element2d & GetFace (int i) const { return faces.Get(i); }
  ///
  const Element & GetElement (int i) const { return elements.Get(i); }
  ///
  const twoint & GetEdge (int i) const { return edges.Get(i); }
  ///
  int GetDelFace (int i) const { return delfaces.Get(i); }
  ///
  int IsDelFace (int fn) const;
  
  ///
  float CalcPointDist (int pi, const Point3d & p) const;
  ///
  double PointDistFactor (int pi) const
    {
      return tolerances.Get(pi);
    }
  ///
  void SetFreeZoneTransformation (const Vector & allp,
				  int tolclass);
  ///
  int IsInFreeZone (const Point3d & p) const;
  /**
    0 not in free-zone
    1 in free-zone
    -1 maybe 
   */
  int IsTriangleInFreeZone (const Point3d & p1, const Point3d & p2,
                            const Point3d & p3, const ARRAY<int> & pi, int newone);
  ///
  int IsQuadInFreeZone (const Point3d & p1, const Point3d & p2,
			const Point3d & p3, const Point3d & p4,
			const ARRAY<int> & pi, int newone);
  ///
  int IsTriangleInFreeSet (const Point3d & p1, const Point3d & p2,
                           const Point3d & p3, int fs, const ARRAY<int> & pi, int newone);

  ///
  int IsQuadInFreeSet (const Point3d & p1, const Point3d & p2,
		       const Point3d & p3, const Point3d & p4,
		       int fs, const ARRAY<int> & pi, int newone);
  
  ///
  int ConvexFreeZone () const;
  
  /// if t1 and t2 are neighbourtriangles, NTP returns the opposite Point of t1 in t2
  int NeighbourTrianglePoint (const threeint & t1, const threeint & t2) const;
  ///
  const Point3d & GetTransFreeZone (int i) { return transfreezone.Get(i); }

  ///
  int GetNP (int fn) const
  { return faces.Get(fn).GetNP(); }
  ///
  int GetPointNr (int fn, int endp) const
  { return faces.Get(fn).PNum(endp); }
  ///
  int GetPointNrMod (int fn, int endp) const
  { return faces.Get(fn).PNumMod(endp); }
  ///
  const fourint & GetOrientation (int i) { return orientations.Get(i); }

  ///
  int TestFlag (char flag) const;

  ///
  const DenseMatrix & GetOldUToNewU () const { return oldutonewu; }
  //
  //  const DenseMatrix & GetOldUToFreeZone () const { return oldutofreezone; }
  //
  //  const DenseMatrix & GetOldUToFreeZoneLimit () const 
  //    { return oldutofreezonelimit; }
  ///
  const char * Name () const { return name; }
  ///
  void LoadRule (istream & ist);

  ///
  const ARRAY<Point3d> & GetTransFreeZone () { return transfreezone; }
  ///
  int TestOk () const;

  ///
  friend void TestRules ();
  ///
  //  friend void Plot3DRule (const ROT3D & r, char key);
};



#endif

