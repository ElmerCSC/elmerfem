///
class splinesegment3d
  {
  ///
  Point<3> p1, p2, p3;
  
  public:
  ///
  splinesegment3d (const Point<3> & ap1, const Point<3> & ap2, 
  	const Point<3> & ap3);
  ///
  void Evaluate (double t, Point<3> & p) const;
  ///
  void EvaluateTangent (double t, Vec<3> & tang) const;
  ///
  const Point<3> & P1() const { return p1; }
  ///
  const Point<3> & P2() const { return p2; }
  ///
  const Point<3> & P3() const { return p3; }
  };

///
class spline3d
  {
  ///
  ARRAY<splinesegment3d *> segments;
  
  public:
  ///
  spline3d () { };
  ///
  void AddSegment (const Point<3> & ap1, const Point<3> & ap2, const Point<3> & ap3);
  ///
  int GetNumSegments () const { return segments.Size(); }
  ///
  double ProjectToSpline (Point<3> & p) const;
  ///
  double ProjectToSpline (Point<3> & p, double t) const;
  ///
  void Evaluate (double t, Point<3> & p) const;
  ///
  void EvaluateTangent (double t, Vec<3> & tang) const;
  ///
  const Point<3> & P1(int i) const { return segments.Get(i)->P1(); }
  ///
  const Point<3> & P2(int i) const { return segments.Get(i)->P2(); }
  ///
  const Point<3> & P3(int i) const { return segments.Get(i)->P3(); }
  };
  
///
class splinetube : public Surface
  {
  ///
  const spline3d & middlecurve;
  ///
  double r;
///  Vec<3> ex, ey, ez;
  Vec<2> e2x, e2y;
  ///
  Point<3> cp;
  
  public:
  ///
  splinetube (const spline3d & amiddlecurve, double ar);
  
  ///
  virtual void DefineTangentialPlane (const Point<3> & ap1, const Point<3> & ap2);
  ///
  virtual void ToPlane (const Point<3> & p, Point<2> & pplain, double h, int & zone) const;
  ///
  virtual void FromPlane (const Point<2> & pplain, Point<3> & p, double h) const;
  ///
  virtual void Project (Point<3> & p) const;

//  virtual int RootInBox (const box3d & box) const { return 0; }
    /// 0 .. no, 1 .. yes, 2 .. maybe

  virtual int BoxInSolid (const BoxSphere<3> & box) const;
    /// 0 .. no, 1 .. yes, 2 .. maybe

  virtual double CalcFunctionValue (const Point<3> & point) const;
  ///
  virtual void CalcGradient (const Point<3> & point, Vec<3> & grad) const;
  ///
  virtual double HesseNorm () const { return 0.5 / r; }
  ///
  virtual Point<3> GetSurfacePoint () const;
  ///
  virtual void Print (ostream & str) const;
  };  
