#ifndef FILE_EXPLICITCURVE2D
#define FILE_EXPLICITCURVE2D

/**************************************************************************/
/* File:   explicitcurve2d.hh                                             */
/* Author: Joachim Schoeberl                                              */
/* Date:   14. Oct. 96                                                    */
/**************************************************************************/

/*

  Explicit 2D Curve repesentation

*/



///
class ExplicitCurve2d : public Curve2d
{
public:
  ///
  ExplicitCurve2d ();

  ///
  virtual void Project (Point<2> & p) const;
  ///
  virtual double ProjectParam (const Point<2> & p) const = 0;
  ///
  virtual double NumericalProjectParam (const Point<2> & p, double lb, double ub) const;
  ///
  virtual double MinParam () const = 0;
  ///
  virtual double MaxParam () const = 0;
  ///
  virtual Point<2> Eval (double t) const = 0;
  ///
  virtual Vec<2> EvalPrime (double t) const = 0;
  ///
  virtual Vec<2> Normal (double t) const;
  ///
  virtual void NormalVector (const Point<2> & p, Vec<2> & n) const;
  ///
  virtual Vec<2> EvalPrimePrime (double t) const = 0;

  ///
  virtual double MaxCurvature () const;
  ///
  virtual double MaxCurvatureLoc (const Point<2> & p, double rad) const;

  ///
  virtual Point<2> CurvCircle (double t) const;
  ///
  virtual void Print (ostream & /* str */) const { };
  
  ///
  virtual int SectionUsed (double /* t */) const { return 1; }
  ///
  virtual void Reduce (const Point<2> & /* p */, double /* rad */) { };
  ///
  virtual void UnReduce () { };
}; 
  
  
///
class BSplineCurve2d : public ExplicitCurve2d
{
  ///
  ARRAY<Point<2> > points;
  ///
  ARRAY<int> intervallused;
  ///
  int redlevel;
  
public:
  ///
  BSplineCurve2d ();
  ///
  void AddPoint (const Point<2> & apoint);

  bool Inside (const Point<2> & p, double & dist) const;
  
  ///
  virtual double ProjectParam (const Point<2> & p) const;
  ///
  virtual double MinParam () const { return 0; }
  ///
  virtual double MaxParam () const { return points.Size(); }
  ///
  virtual Point<2> Eval (double t) const;
  ///
  virtual Vec<2> EvalPrime (double t) const;  
  ///
  virtual Vec<2> EvalPrimePrime (double t) const;
  ///
  virtual void Print (ostream & str) const;

  ///
  virtual int SectionUsed (double t) const;
  ///
  virtual void Reduce (const Point<2> & p, double rad);
  ///
  virtual void UnReduce ();
};  




#endif
