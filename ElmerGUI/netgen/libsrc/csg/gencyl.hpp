#ifndef FILE_GENCYL
#define FILE_GENCYL

/**************************************************************************/
/* File:   gencyl.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   14. Oct. 96                                                    */
/**************************************************************************/

/*
  
  Generalized Cylinder
  
*/


///
class GeneralizedCylinder : public Surface
{
  ///
  ExplicitCurve2d & crosssection;
  ///
  Point<3> planep;
  ///
  Vec<3> planee1, planee2, planee3;
  
  ///  Vec<3> ex, ey, ez;
  Vec2d e2x, e2y;
    ///
  Point<3> cp;
  
public:
  ///
  GeneralizedCylinder (ExplicitCurve2d & acrosssection,
		       Point<3> ap, Vec<3> ae1, Vec<3> ae2);
  
  ///
  virtual void Project (Point<3> & p) const;
  
  ///
  virtual int BoxInSolid (const BoxSphere<3> & box) const;
  /// 0 .. no, 1 .. yes, 2 .. maybe
  
  virtual double CalcFunctionValue (const Point<3> & point) const;
  ///
  virtual void CalcGradient (const Point<3> & point, Vec<3> & grad) const;
  ///
  virtual void CalcHesse (const Point<3> & point, Mat<3> & hesse) const;
  ///
  virtual double HesseNorm () const;
  ///
  virtual double MaxCurvatureLoc (const Point<3> & c, double rad) const;
  ///
  virtual Point<3> GetSurfacePoint () const;
  ///
  virtual void Print (ostream & str) const;
  
  ///
  virtual void Reduce (const BoxSphere<3> & box);
  ///
  virtual void UnReduce ();
};  

#endif
