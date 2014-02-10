



das File sollte nicht mehr verwendet werden ---> spline.hpp







#ifndef FILE_SPLINE2D
#define FILE_SPLINE2D

/**************************************************************************/
/* File:   spline2d.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   24. Jul. 96                                                    */
/**************************************************************************/


/*
  Spline curves for 2D mesh generation
  */

#include "spline.hpp"


//#define OLDSPLINEVERSION
#ifdef OLDSPLINEVERSION

/// Geometry point
class GeomPoint2d : public Point<2>
{
public:
  /// refinement to point
  double refatpoint;
  bool hpref;

  GeomPoint2d ()
  { ; }

  ///
  GeomPoint2d (double ax, double ay, double aref = 1)
    : Point<2> (ax, ay), refatpoint(aref) { ; }
};



/// base class for 2d - segment
class SplineSegment
{
public:
  /// left domain
  int leftdom;
  /// right domain
  int rightdom;
  /// refinement at line
  double reffak;
  /// boundary condition number
  int bc;
  /// copy spline mesh from other spline (-1.. do not copy)
  int copyfrom;
  /// perfrom anisotropic refinement (hp-refinement) to edge
  bool hpref_left;
  bool hpref_right;
  /// calculates length of curve
  virtual double Length () const;
  /// returns point at curve, 0 <= t <= 1
  virtual Point<2> GetPoint (double t) const = 0;
  /// partitionizes curve
  void Partition (double h, double elto0,
		  Mesh & mesh, Point3dTree & searchtree, int segnr) const;
  /// returns initial point on curve
  virtual const GeomPoint2d & StartPI () const = 0;
  /// returns terminal point on curve
  virtual const GeomPoint2d & EndPI () const = 0;
  /** writes curve description for fepp:
    for implicitly given quadratic curves, the 6 coefficients of
    the polynomial
    $$ a x^2 + b y^2 + c x y + d x + e y + f = 0 $$
    are written to ost */
  void PrintCoeff (ostream & ost) const;

  virtual void GetCoeff (Vector & coeffs) const = 0;

  virtual void GetPoints (int n, ARRAY<Point<2> > & points);

  /** calculates lineintersections:
      for lines $$ a x + b y + c = 0 $$ the interecting points are calculated
      and stored in points */
  virtual void LineIntersections (const double a, const double b, const double c,
				  ARRAY < Point<2> > & points, const double eps) const
  {points.SetSize(0);}

  virtual double MaxCurvature(void) const = 0;

  virtual string GetType(void) const {return "splinebase";}
};


/// Straight line form p1 to p2
class LineSegment : public SplineSegment
{
  ///
  const GeomPoint2d &p1, &p2;
public:
  ///
  LineSegment (const GeomPoint2d & ap1, const GeomPoint2d & ap2);
  ///
  virtual double Length () const;
  ///
  virtual Point<2> GetPoint (double t) const;
  ///
  virtual const GeomPoint2d & StartPI () const { return p1; };
  ///
  virtual const GeomPoint2d & EndPI () const { return p2; }
  ///
  //virtual void PrintCoeff (ostream & ost) const;
  virtual void GetCoeff (Vector & coeffs) const;

  virtual string GetType(void) const {return "line";}

  virtual void LineIntersections (const double a, const double b, const double c,
				  ARRAY < Point<2> > & points, const double eps) const;

  virtual double MaxCurvature(void) const {return 0;}
};


/// curve given by a rational, quadratic spline (including ellipses)
class SplineSegment3 : public SplineSegment
{
  ///
  const GeomPoint2d &p1, &p2, &p3;
public:
  ///
  SplineSegment3 (const GeomPoint2d & ap1, 
		  const GeomPoint2d & ap2, 
		  const GeomPoint2d & ap3);
  ///
  virtual Point<2> GetPoint (double t) const;
  ///
  virtual const GeomPoint2d & StartPI () const { return p1; };
  ///
  virtual const GeomPoint2d & EndPI () const { return p3; }
  ///
  //virtual void PrintCoeff (ostream & ost) const;
  virtual void GetCoeff (Vector & coeffs) const;

  virtual string GetType(void) const {return "spline3";}

  const GeomPoint2d & TangentPoint (void) const { return p2; }

  virtual void LineIntersections (const double a, const double b, const double c,
				  ARRAY < Point<2> > & points, const double eps) const;

  virtual double MaxCurvature(void) const;
};


// Gundolf Haase  8/26/97
/// A circle
class CircleSegment : public SplineSegment
{
  ///
private:
  const GeomPoint2d	&p1, &p2, &p3;
  Point<2>		pm;
  double		radius, w1,w3;
public:
  ///
  CircleSegment (const GeomPoint2d & ap1, 
		 const GeomPoint2d & ap2, 
		 const GeomPoint2d & ap3);
  ///
  virtual Point<2> GetPoint (double t) const;
  ///
  virtual const GeomPoint2d & StartPI () const { return p1; }
  ///
  virtual const GeomPoint2d & EndPI () const { return p3; }
  ///
  //virtual void PrintCoeff (ostream & ost) const;
  virtual void GetCoeff (Vector & coeffs) const;
  ///
  double Radius() const { return radius; }
  ///
  double StartAngle() const { return w1; }
  ///
  double EndAngle() const { return w3; }
  ///
  const Point<2> & MidPoint(void) const {return pm; }

  virtual string GetType(void) const {return "circle";}

  virtual void LineIntersections (const double a, const double b, const double c,
				  ARRAY < Point<2> > & points, const double eps) const;

  virtual double MaxCurvature(void) const {return 1./radius;}
};






/// 
class DiscretePointsSegment : public SplineSegment
{
  ARRAY<Point<2> > pts;
  GeomPoint2d p1, p2;
public:
  ///
  DiscretePointsSegment (const ARRAY<Point<2> > & apts);
  ///
  virtual ~DiscretePointsSegment ();
  ///
  virtual Point<2> GetPoint (double t) const;
  ///
  virtual const GeomPoint2d & StartPI () const { return p1; };
  ///
  virtual const GeomPoint2d & EndPI () const { return p2; }
  ///
  //virtual void PrintCoeff (ostream & /* ost */) const { ; }
  virtual void GetCoeff (Vector & coeffs) const {;}

  virtual double MaxCurvature(void) const {return 1;}
};


#endif

#endif
