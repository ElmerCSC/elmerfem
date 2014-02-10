
#ifndef FILE_IDENTIFY
#define FILE_IDENTIFY

/**************************************************************************/
/* File:   identify.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   1. Aug. 99                                                    */
 /**************************************************************************/

/**
   Identify surfaces for periodic b.c. or
   thin domains
*/


class SpecialPoint;
class Identification
{
protected:
  const CSGeometry & geom;
  // identified faces, index sorted
  INDEX_2_HASHTABLE<int> identfaces;
  int nr;

public:
  Identification (int anr, const CSGeometry & ageom);
  virtual ~Identification ();
  virtual void Print (ostream & ost) const = 0;
  virtual void GetData (ostream & ost) const = 0;

  /// obsolete
  //  virtual void IdentifySpecialPoints (ARRAY<class SpecialPoint> & points);

  /// can identify both special points (fixed direction)
  /// (identified points, same tangent)
  virtual int Identifyable (const SpecialPoint & sp1, const SpecialPoint & sp2,
			    const TABLE<int> & specpoint2solid,			  
			    const TABLE<int> & specpoint2surface) const;
  ///
  virtual int Identifyable (const Point<3> & p1, const Point<3> & sp2) const;
  /// is it possible to identify sp1 with some other ?
  virtual int IdentifyableCandidate (const SpecialPoint & sp1) const;
  
  /// are points (if connected) by a short edge (direction anyhow) ?
  virtual int ShortEdge (const SpecialPoint & sp1, const SpecialPoint & sp2) const;

  /// add entries in mesh identification tables
  virtual void IdentifyPoints (class Mesh & mesh);

  /// add entries to identified faces (based on segment infos)
  virtual void IdentifyFaces (class Mesh & mesh);

  /// get point on other surface, add entry in mesh identifications
  virtual int GetIdentifiedPoint (class Mesh & mesh, int pi1);

  /// copy surfaces, or fill rectangles
  virtual void BuildSurfaceElements (ARRAY<class Segment> & segs,
				     class Mesh & mesh,
				     const Surface * surf);

  /// insert volume elements in thin layers
  virtual void BuildVolumeElements (ARRAY<class Element2d> & surfels,
				    class Mesh & mesh);

  /// get list of identified faces
  virtual void GetIdentifiedFaces (ARRAY<INDEX_2> & idfaces) const;

  friend ostream & operator<< (ostream & ost, Identification & ident);
};


class PeriodicIdentification : public Identification
{
  const Surface * s1;
  const Surface * s2;
public:
  PeriodicIdentification (int anr,
			  const CSGeometry & ageom,
			  const Surface * as1,
			  const Surface * as2);
  virtual ~PeriodicIdentification ();
  virtual void Print (ostream & ost) const;
  virtual void GetData (ostream & ost) const;


  //  virtual void IdentifySpecialPoints (ARRAY<class SpecialPoint> & points);
  virtual int Identifyable (const SpecialPoint & sp1, const SpecialPoint & sp2,
			    const TABLE<int> & specpoint2solid,
			    const TABLE<int> & specpoint2surface) const;

  virtual int Identifyable (const Point<3> & p1, const Point<3> & sp2) const;
  virtual int GetIdentifiedPoint (class Mesh & mesh, int pi1);
  virtual void IdentifyPoints (class Mesh & mesh);
  virtual void IdentifyFaces (class Mesh & mesh);
  virtual void BuildSurfaceElements (ARRAY<class Segment> & segs,
				     class Mesh & mesh,
				     const Surface * surf);
};


///
class TopLevelObject;
class CloseSurfaceIdentification : public Identification
{
  const Surface * s1;
  const Surface * s2;
  const TopLevelObject * domain;
  ///
  int dom_nr;
  /// number of refinement levels (in Z-refinement)
  int ref_levels;
  /// number of refinement levels for layer next to s1 (in Z-refinement)
  int ref_levels_s1;
  /// number of refinement levels for layer next to s2 (in Z-refinement)
  int ref_levels_s2;
  ///
  double eps_n;
  ARRAY<double> slices;
  /// used only for domain-local identification:
  ARRAY<int> domain_surfaces;
  ///
  bool dom_surf_valid;

  ///
  Vec<3> direction;
  ///
  bool usedirection;
public:
  CloseSurfaceIdentification (int anr, 
			      const CSGeometry & ageom,
			      const Surface * as1,
			      const Surface * as2,
			      const TopLevelObject * adomain,
			      const Flags & flags);
  virtual ~CloseSurfaceIdentification ();

  virtual void Print (ostream & ost) const;
  virtual void GetData (ostream & ost) const;


  //  virtual void IdentifySpecialPoints (ARRAY<class SpecialPoint> & points);
  virtual int Identifyable (const SpecialPoint & sp1, const SpecialPoint & sp2,
			    const TABLE<int> & specpoint2solid,
			    const TABLE<int> & specpoint2surface) const;
  virtual int Identifyable (const Point<3> & p1, const Point<3> & sp2) const;
  virtual int IdentifyableCandidate (const SpecialPoint & sp1) const;
  virtual int ShortEdge (const SpecialPoint & sp1, const SpecialPoint & sp2) const;
  virtual int GetIdentifiedPoint (class Mesh & mesh, int pi1);
  const ARRAY<double> & GetSlices () const { return slices; }
  virtual void IdentifyPoints (class Mesh & mesh);
  virtual void IdentifyFaces (class Mesh & mesh);
  virtual void BuildSurfaceElements (ARRAY<class Segment> & segs,
				     class Mesh & mesh,
				     const Surface * surf);
  void BuildSurfaceElements2 (ARRAY<class Segment> & segs,
			      class Mesh & mesh,
			      const Surface * surf);

  virtual void BuildVolumeElements (ARRAY<class Element2d> & surfels,
				    class Mesh & mesh);

  int RefLevels () const { return ref_levels; }
  int RefLevels1 () const { return ref_levels_s1; }
  int RefLevels2 () const { return ref_levels_s2; }

  bool IsSkewIdentification(void) const {return usedirection;}
  const Vec<3> & GetDirection(void) const {return direction;}

  const Surface & GetSurface1(void) const
  { return *s1;}
  const Surface & GetSurface2(void) const
  { return *s2;}
};


class CloseEdgesIdentification : public Identification
{
  const Surface * facet;
  const Surface * s1;
  const Surface * s2;
public:
  CloseEdgesIdentification (int anr,
			    const CSGeometry & ageom,
			    const Surface * afacet,
			    const Surface * as1,
			    const Surface * as2);
  virtual ~CloseEdgesIdentification ();
  virtual void Print (ostream & ost) const;
  virtual void GetData (ostream & ost) const;

  //  virtual void IdentifySpecialPoints (ARRAY<class SpecialPoint> & points);
  virtual int Identifyable (const SpecialPoint & sp1, const SpecialPoint & sp2,
			    const TABLE<int> & specpoint2solid,
			    const TABLE<int> & specpoint2surface) const;


  virtual void IdentifyPoints (class Mesh & mesh);
  virtual void BuildSurfaceElements (ARRAY<class Segment> & segs,
				     class Mesh & mesh,
				     const Surface * surf);
};

#endif
