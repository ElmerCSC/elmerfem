#ifndef FILE_ADFRONT3
#define FILE_ADFRONT3

/**************************************************************************/
/* File:   adfront3.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

/*
  Advancing front class for volume meshing
*/



/// Point in advancing front
class FrontPoint3
{
  /// coordinates
Point<3> p;           
  /// global node index
PointIndex globalindex;   
  /// number of faces connected to point 
int nfacetopoint;    
  /// distance to original boundary
int frontnr;
  /// 
int cluster;
public:
  ///
  FrontPoint3 ();
  ///
  FrontPoint3 (const Point<3> & ap, PointIndex agi);
  
  ///
  const Point<3> & P () const
  { return p; }
  ///
  PointIndex GlobalIndex () const
  { return globalindex; }
  
  ///
  void AddFace ()
  { nfacetopoint++; }

  ///
  void RemoveFace()
  { 
    nfacetopoint--;
    if (nfacetopoint == 0) nfacetopoint = -1;
  }
  
  ///
  int Valid () const
  { return nfacetopoint >= 0; }

  ///
  void DecFrontNr (int afrontnr)
  {
    if (frontnr > afrontnr) frontnr = afrontnr;
  }
  
  ///
  int FrontNr () const
  { return frontnr; }

  ///
  friend class AdFront3;
};



class MiniElement2d
{
protected:
  int np;
  PointIndex pnum[4];
  bool deleted;
public:
  MiniElement2d ()
  { np = 3; deleted = 0; }
  MiniElement2d (int anp)
  { np = anp; deleted = 0; }

  int GetNP() const { return np; }
  PointIndex & operator[] (int i) { return pnum[i]; }
  const PointIndex operator[] (int i) const { return pnum[i]; }

  const PointIndex PNum (int i) const { return pnum[i-1]; }
  PointIndex & PNum (int i) { return pnum[i-1]; }
  const PointIndex PNumMod (int i) const { return pnum[(i-1)%np]; }

  void Delete () { deleted = 1; pnum[0] = pnum[1] = pnum[2] = pnum[3] = PointIndex::BASE-1; }
  bool IsDeleted () const { return deleted; }
};


inline ostream & operator<<(ostream  & s, const MiniElement2d & el)
{
  s << "np = " << el.GetNP();
  for (int j = 0; j < el.GetNP(); j++)
    s << " " << el[j];
  return s;
}




/// Face in advancing front
class FrontFace
{
private:
  ///
  MiniElement2d f;
  ///
  int qualclass;
  ///
  char oldfront;
  ///
  int hashvalue;
  ///
  int cluster;

public:
  ///
  FrontFace ();
  ///
  FrontFace (const MiniElement2d & af);
  ///
  const MiniElement2d & Face () const
  { return f; }
  
  ///
  int QualClass () const
  { return qualclass; }

  ///
  void IncrementQualClass ()
  { qualclass++; }

  ///
  void ResetQualClass ()
  {
    if (qualclass > 1)
      {
	qualclass = 1;
	oldfront = 0;
      }
  }
  
  ///
  bool Valid () const
  { return !f.IsDeleted(); }

  ///
  void Invalidate ();

  ///
  int HashValue() const 
  { return hashvalue; }

  ///
  void SetHashValue(int hv) 
  { hashvalue = hv; }

  ///
  friend class AdFront3;

  int Cluster () const { return cluster; }
};  




/// Advancing front, 3D.
class AdFront3
{
  ///
ARRAY<FrontPoint3, PointIndex::BASE> points;
  ///
ARRAY<FrontFace> faces;
  ///
ARRAY<PointIndex> delpointl;

  /// which points are connected to pi ?
TABLE<int, PointIndex::BASE> * connectedpairs;
  
  /// number of total front faces;
int nff;
  /// number of quads in front
int nff4; 
  
  ///
double vol;

  ///
GeomSearch3d hashtable;

  /// 
int hashon;

  ///
int hashcreated;

  /// counter for rebuilding internal tables
int rebuildcounter;
  /// last base element
int lasti;
  /// minimal selection-value of baseelements
int minval;

  ///
class Box3dTree * facetree;
public:

  ///
  AdFront3 ();
  ///
  ~AdFront3 ();
  ///
  void GetPoints (ARRAY<Point<3> > & apoints) const;
  ///
  int GetNP() const 
  { return points.Size(); }
  ///
  const Point<3> & GetPoint (PointIndex pi) const
  { return points[pi].P(); }
  ///
  int GetNF() const
  { return nff; }
  ///
  const MiniElement2d & GetFace (int i) const
  { return faces.Get(i).Face(); }
  ///
  void Print () const;
  ///
  bool Empty () const
  { return nff == 0; }
  ///
  bool Empty (int elnp) const
  {
    if (elnp == 4)
      return (nff4 == 0);
    return (nff - nff4 == 0);
  }
  ///
  int SelectBaseElement ();

  ///
  void CreateTrees ();

  ///
  void GetIntersectingFaces (const Point<3> & pmin, const Point<3> & pmax, 
			     ARRAY<int> & ifaces) const;

  ///
  void GetFaceBoundingBox (int i, Box3d & box) const;

  ///
  int GetLocals (int baseelement,
		 ARRAY<Point3d > & locpoints,
                 ARRAY<MiniElement2d> & locfaces,   // local index
                 ARRAY<PointIndex> & pindex,
                 ARRAY<INDEX> & findex,
		 INDEX_2_HASHTABLE<int> & connectedpairs,
                 float xh,
		 float relh,
		 INDEX& facesplit);
  
  ///
  void GetGroup (int fi,
                 ARRAY<MeshPoint> & grouppoints,
                 ARRAY<MiniElement2d> & groupelements,
                 ARRAY<PointIndex> & pindex,
                 ARRAY<INDEX> & findex
                 ) const;

  ///
  void DeleteFace (INDEX fi);
  ///
  PointIndex AddPoint (const Point<3> & p, PointIndex globind);
  ///
  INDEX AddFace (const MiniElement2d & e);
  ///
  INDEX AddConnectedPair (const INDEX_2 & pair);
  ///
  void IncrementClass (INDEX fi)
  { faces.Elem(fi).IncrementQualClass(); }

  ///
  void ResetClass (INDEX fi)
  { faces.Elem(fi).ResetQualClass(); }

  ///
  void SetStartFront (int baseelnp = 0);

  /// is Point p inside Surface ?
  bool Inside (const Point<3> & p) const;
  /// both points on same side ?
  int SameSide (const Point<3> & lp1, const Point<3> & lp2, 
		const ARRAY<int> * testfaces = NULL) const;


  ///
  PointIndex GetGlobalIndex (PointIndex pi) const
  { return points[pi].GlobalIndex(); }
  ///
  double Volume () const
  { return vol; }


private:
  void RebuildInternalTables();
};




#endif
