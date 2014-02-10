

OLD IMPLEMENTATION, NOT USED ANYMORE!




#ifndef FILE_SPLINEGEOMETRY2
#define FILE_SPLINEGEOMETRY2

/**************************************************************************/
/* File:   splinegeometry2.hpp                                            */
/* Author: Joachim Schoeberl                                              */
/* Date:   24. Jul. 96                                                    */
/**************************************************************************/

#include "splinegeometry.hpp"

#ifdef OLDSPLINEGEOMETRY


/// 
extern void LoadBoundarySplines (const char * filename,
				 ARRAY<GeomPoint2d> & geompoints,
				 ARRAY<SplineSegment*> & splines, 
				 double & elto0);
///
extern void PartitionBoundary (const ARRAY<SplineSegment*> & splines,
			       double h, double elto0,
			       Mesh & mesh2d);


class CSGScanner;

class SplineGeometry2d
{
  ARRAY<GeomPoint2d> geompoints;
  ARRAY<SplineSegment*> splines;
  double elto0;


private:
  void AppendSegment(SplineSegment * spline, const int leftdomain, const int rightdomain,
		     const int bc,
		     const double reffac, const bool hprefleft, const bool hprefright,
		     const int copyfrom);

public:
  ~SplineGeometry2d();

  void Load (const char * filename);
  void CSGLoad (CSGScanner & scan);
  void PartitionBoundary (double h, Mesh & mesh2d);

  void CopyEdgeMesh (int from, int to, Mesh & mesh2d, Point3dTree & searchtree);

  const ARRAY<SplineSegment*> & GetSplines () const
  { return splines; }

  int GetNSplines (void) const { return splines.Size(); }
  string GetSplineType (const int i) const { return splines[i]->GetType(); }
  SplineSegment & GetSpline (const int i) {return *splines[i];}
  const SplineSegment & GetSpline (const int i) const {return *splines[i];}

  void GetBoundingBox (Box<2> & box) const;

  int GetNP () const { return geompoints.Size(); }
  const GeomPoint2d & GetPoint(int i) const { return geompoints[i]; }

  void SetGrading (const double grading);
  void AppendPoint (const double x, const double y, const double reffac = 1., const bool hpref = false);
  
  void AppendLineSegment (const int n1, const int n2,
			  const int leftdomain, const int rightdomain, const int bc = -1,
			  const double reffac = 1.,
			  const bool hprefleft = false, const bool hprefright = false,
			  const int copyfrom = -1);
  void AppendSplineSegment (const int n1, const int n2, const int n3,
			    const int leftdomain, const int rightdomain, const int bc = -1,
			    const double reffac = 1.,
			    const bool hprefleft = false, const bool hprefright = false,
			    const int copyfrom = -1);
  void AppendCircleSegment (const int n1, const int n2, const int n3,
			    const int leftdomain, const int rightdomain, const int bc = -1,
			    const double reffac = 1.,
			    const bool hprefleft = false, const bool hprefright = false,
			    const int copyfrom = -1);
  void AppendDiscretePointsSegment (const ARRAY< Point<2> > & points, 
				    const int leftdomain, const int rightdomain, const int bc = -1,
				    const double reffac = 1.,
				    const bool hprefleft = false, const bool hprefright = false,
				    const int copyfrom = -1);
  
};


void MeshFromSpline2D (SplineGeometry2d & geometry,
		       Mesh *& mesh, 
		       MeshingParameters & mp);

#endif

#endif
