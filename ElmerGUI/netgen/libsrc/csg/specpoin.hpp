#ifndef FILE_SPECPOIN
#define FILE_SPECPOIN


/**************************************************************************/
/* File:   specpoin.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

/*

Special Point Calculation
  
*/

class Surface;
class Solid;

/// Special point.
class SpecialPoint
{
public:
  /// coordinates
  Point<3> p;
  /// tangential to edge
  Vec<3> v;
  ///
  int layer;
  /// point must be used in mesh 
  bool unconditional; 

  /// surfaces defining edge 
  int s1, s2;
  /// if s1 and s2 are only representatives, then these are the original indices
  int s1_orig, s2_orig;
  int nr;
  ///
  SpecialPoint () : p(0,0,0), v(0,0,0), layer(0), unconditional(0), s1(0), s2(0), s1_orig(0), s2_orig(0)
  { ; }

  ///
  SpecialPoint (const SpecialPoint & sp2);

  ///
  SpecialPoint & operator= (const SpecialPoint & sp2);
  
  ///
  void Print (ostream & str) const;


  int GetLayer() const { return layer; }

  ///
  bool HasSurfaces (int as1, int as2) const
  {
    return (s1 == as1 && s2 == as2 || s1 == as2 && s2 == as1);
  }
};

inline ostream & operator<< (ostream & ost, const SpecialPoint & sp)
{
  sp.Print (ost);
  return ost;
}




///
class SpecialPointCalculation
{
private:
  ///
  const CSGeometry * geometry;
  ///
  ARRAY<MeshPoint> * points;
  ///
  ARRAY<long int> boxesinlevel;

  ///
  double size;
  ///
  double relydegtest;   // maximal dimension of bisection intervall for
                        /// test of degeneration parameters
  double cpeps1, epeps1, epeps2, epspointdist2;

  double ideps;

public: 

  ///
  SpecialPointCalculation (); 
  
  ///
  void SetIdEps(const double epsin) {ideps = epsin;}

  ///
  void CalcSpecialPoints (const CSGeometry & ageometry, 
			  ARRAY<MeshPoint> & points);
  ///
  void AnalyzeSpecialPoints (const CSGeometry & geometry, 
			     ARRAY<MeshPoint> & points, 
			     ARRAY<SpecialPoint> & specpoints);

protected:
  ///
  void CalcSpecialPointsRec (const Solid * sol, int layer,
			     const BoxSphere<3> & box, 
			     int level, 
			     bool calccp, bool calcep);


  ///
  bool CrossPointNewtonConvergence (const Surface * f1, const Surface * f2, 
				    const Surface * f3, const BoxSphere<3> & box);  
  ///
  bool CrossPointDegenerated (const Surface * f1, const Surface * f2,
			      const Surface * f3, const BoxSphere<3> & box) const;
  ///
  void CrossPointNewton (const Surface * f1, const Surface * f2, 
			 const Surface * f3, Point<3> & p);
  
  bool EdgeNewtonConvergence (const Surface * f1, const Surface * f2, 
			      const Point<3> & p);  
  ///
  bool EdgeDegenerated (const Surface * f1, const Surface * f2,
			const BoxSphere<3> & box) const;
  ///
  void EdgeNewton (const Surface * f1, const Surface * f2, 
		   Point<3> & p);
  ///
  bool IsEdgeExtremalPoint (const Surface * f1, const Surface * f2, 
			    const Point<3> & p, Point<3> & pp, double rad);



  /*
  ///
  bool ExtremalPointPossible (const Surface * f1, const Surface * f2, 
			      int dir, const BoxSphere<3> & box);
  ///
  bool ExtremalPointDegenerated (const Surface * f1, const Surface * f2, 
				 int dir, const BoxSphere<3> & box);
  ///
  bool ExtremalPointNewtonConvergence (const Surface * f1, const Surface * f2, 
				       int dir, const BoxSphere<3> & box);
  */
  ///
  void ExtremalPointNewton (const Surface * f1, const Surface * f2, 
			    int dir, Point<3> & p);


  ///
  bool AddPoint (const Point<3> & p, int layer);

  void ComputeExtremalPoints (const Plane * plane, 
			      const QuadraticSurface * quadric, 
			      ARRAY<Point<3> > & pts);

  void ComputeCrossPoints (const Plane * plane1, 
			   const Plane * plane2, 
			   const Plane * plane3, 
			   ARRAY<Point<3> > & pts);

  void ComputeCrossPoints (const Plane * plane1, 
			   const Plane * plane2, 
			   const QuadraticSurface * quadratic, 
			   ARRAY<Point<3> > & pts);
};

#endif


