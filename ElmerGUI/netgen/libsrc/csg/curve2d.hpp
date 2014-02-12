#ifndef FILE_CURVE2D
#define FILE_CURVE2D

/**************************************************************************/
/* File:   curve2d.hh                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   24. Jul. 96                                                    */
/**************************************************************************/

/*

  2D Curve repesentation

*/



///
class Curve2d : public Manifold
  {
  public:
  ///
  virtual void Project (Point<2> & p) const = 0;
  ///
  virtual void NormalVector (const Point<2> & p, Vec<2> & n) const = 0;
  };
  
///
class CircleCurve2d : public Curve2d
  {
  ///
  Point<2> center;
  ///
  double rad;
  public:
  ///
  CircleCurve2d (const Point<2> & acenter, double arad);
  ///
  virtual void Project (Point<2> & p) const;
  ///
  virtual void NormalVector (const Point<2> & p, Vec<2> & n) const;
  };
  
///
class QuadraticCurve2d : public Curve2d
{
  ///
  double cxx, cyy, cxy, cx, cy, c;
public:
  ///
  QuadraticCurve2d ();
  ///
  void Read (istream & ist);
  ///
  virtual void Project (Point<2> & p) const;
  ///
  virtual void NormalVector (const Point<2> & p, Vec<2> & n) const;
};
#endif
