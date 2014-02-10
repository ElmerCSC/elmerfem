#ifndef FILE_GEOM2D
#define FILE_GEOM2D

/* *************************************************************************/
/* File:   geom2d.hh                                                       */
/* Author: Joachim Schoeberl                                               */
/* Date:   5. Aug. 95                                                      */
/* *************************************************************************/



/* Geometric Algorithms */

#define EPSGEOM 1E-5


// extern void MyError (const char * ch);

class Point2d;
class Vec2d;

class LINE2D;
class Line2d;
class PLine2d;
class TRIANGLE2D;
class PTRIANGLE2D;


inline Vec2d operator- (const Point2d & p1, const Point2d & p2);
inline Point2d operator- (const Point2d & p1, const Vec2d & v);
inline Point2d operator+ (const Point2d & p1, const Vec2d & v);
inline Point2d Center (const Point2d & p1, const Point2d & p2);

inline void PpSmV (const Point2d & p1, double s, const Vec2d & v, Point2d & p2);
inline void PmP (const Point2d & p1, const Point2d & p2, Vec2d & v);
ostream & operator<<(ostream  & s, const Point2d & p);
inline Vec2d operator- (const Point2d & p1, const Point2d & p2);
inline Point2d operator- (const Point2d & p1, const Vec2d & v);
inline Point2d operator+ (const Point2d & p1, const Vec2d & v);
inline Vec2d operator- (const Vec2d & p1, const Vec2d & v);
inline Vec2d operator+ (const Vec2d & p1, const Vec2d & v);
inline Vec2d operator* (double scal, const Vec2d & v);
double Angle (const Vec2d & v);
double FastAngle (const Vec2d & v);
double Angle (const Vec2d & v1, const Vec2d & v2);
double FastAngle (const Vec2d & v1, const Vec2d & v2);
ostream & operator<<(ostream  & s, const Vec2d & v);
double Dist2(const Line2d & g, const Line2d & h );		// GH
int Near (const Point2d & p1, const Point2d & p2, const double eps);

int Parallel (const Line2d & l1, const Line2d & l2, double peps = EPSGEOM);
int IsOnLine (const Line2d & l, const Point2d & p, double heps = EPSGEOM);
int IsOnLongLine (const Line2d & l, const Point2d & p);
int Hit (const Line2d & l1, const Line2d & l2, double heps = EPSGEOM);
ostream & operator<<(ostream  & s, const Line2d & l);
Point2d CrossPoint (const PLine2d & l1, const PLine2d & l2);
int Parallel (const PLine2d & l1, const PLine2d & l2, double peps = EPSGEOM);
int IsOnLine (const PLine2d & l, const Point2d & p, double heps = EPSGEOM);
int IsOnLongLine (const PLine2d & l, const Point2d & p);
int Hit (const PLine2d & l1, const Line2d & l2, double heps = EPSGEOM);
ostream & operator<<(ostream  & s, const Line2d & l);
ostream & operator<<(ostream  & s, const TRIANGLE2D & t); 
ostream & operator<<(ostream & s, const PTRIANGLE2D & t);

///
class Point2d
{
  ///
  friend class Vec2d;

protected:
  ///
  double px, py;

public:
  ///
  Point2d() { /* px = py = 0; */ }
  ///
  Point2d(double ax, double ay) { px = ax; py = ay; }
  ///
  Point2d(const Point2d & p2) { px = p2.px; py = p2.py; }

  Point2d (const Point<2> & p2)
  {
    px = p2(0);
    py = p2(1);
  }
  ///
  Point2d & operator= (const Point2d & p2)
    { px = p2.px; py = p2.py; return *this; }
    
  ///
  int operator== (const Point2d & p2) const			// GH
    { return (px == p2.px  &&  py == p2.py) ; }

  ///
  double & X() { return px; }
  ///
  double & Y() { return py; }
  ///
  double X() const { return px; }
  ///
  double Y() const { return py; }

  operator Point<2> () const
  {
    return Point<2> (px, py);
  }


  ///
  friend inline Vec2d operator- (const Point2d & p1, const Point2d & p2);
  ///
  friend inline Point2d operator- (const Point2d & p1, const Vec2d & v);
  ///
  friend inline Point2d operator+ (const Point2d & p1, const Vec2d & v);

  ///
  friend inline Point2d Center (const Point2d & p1, const Point2d & p2);

  const Point2d & SetToMin (const Point2d & p2)
    {
      if (p2.px < px) px = p2.px;
      if (p2.py < py) py = p2.py;
      return *this;
    }


  ///
  const Point2d & SetToMax (const Point2d & p2)
    {
      if (p2.px > px) px = p2.px;
      if (p2.py > py) py = p2.py;
      return *this;
    }

  ///
  friend double Dist (const Point2d & p1, const Point2d & p2)
    { return sqrt ( (p1.px - p2.px) * (p1.px - p2.px) +
		    (p1.py - p2.py) * (p1.py - p2.py) ); }
  //    { return sqrt ( sqr (p1.X()-p2.X()) + sqr (p1.Y()-p2.Y()) ); }

  ///
  friend double Dist2 (const Point2d & p1, const Point2d & p2)
    { return ( (p1.px - p2.px) * (p1.px - p2.px) +
	       (p1.py - p2.py) * (p1.py - p2.py) ); }
  //    { return sqr (p1.X()-p2.X()) + sqr (p1.Y()-p2.Y()) ; }


  /**
    Points clock-wise ?
    Are the points (p1, p2, p3) clock-wise ?
    */
  friend inline int CW (const Point2d & p1, const Point2d & p2, const Point2d & p3)
    {
      //      return Cross (p2 - p1, p3 - p2) < 0;      
      return
	(p2.px - p1.px) * (p3.py - p2.py) - 
	(p2.py - p1.py) * (p3.px - p2.px) < 0;
    }
  /**
    Points counter-clock-wise ?
    Are the points (p1, p2, p3) counter-clock-wise ?
    */
  friend inline bool CCW (const Point2d & p1, const Point2d & p2, const Point2d & p3)
    {
      //      return Cross (p2 - p1, p3 - p2) > 0;
      return
	(p2.px - p1.px) * (p3.py - p2.py) - 
	(p2.py - p1.py) * (p3.px - p2.px) > 0;
    }  /**
    Points counter-clock-wise ?
    Are the points (p1, p2, p3) counter-clock-wise ?
    */
  friend inline bool CCW (const Point2d & p1, const Point2d & p2, const Point2d & p3, double eps)
    {
      //      return Cross (p2 - p1, p3 - p2) > 0;
      double ax = p2.px - p1.px;
      double ay = p2.py - p1.py;
      double bx = p3.px - p2.px;
      double by = p3.py - p2.py;

      return ax*by - ay*bx > eps*eps*max2(ax*ax+ay*ay,bx*bx+by*by);
    }

  ///
  friend inline void PpSmV (const Point2d & p1, double s, const Vec2d & v, Point2d & p2);
  ///
  friend inline void PmP (const Point2d & p1, const Point2d & p2, Vec2d & v);

  ///
  friend ostream & operator<<(ostream  & s, const Point2d & p);
};


inline int Near (const Point2d & p1, const Point2d & p2, 
	  const double eps = 1e-4 )
{ 
  return  Dist2(p1,p2) <= eps*eps; 
}






///
class Vec2d
  {
protected:
  ///
  double vx, vy;

public:
  ///
    Vec2d() { /* vx = vy = 0; */ }
    ///
    Vec2d(double ax, double ay)
    { vx = ax; vy = ay; }
    ///
    Vec2d(const Vec2d & v2) { vx = v2.vx; vy = v2.vy; }

    ///
    explicit Vec2d(const Vec<2> & v2) { vx = v2(0); vy = v2(1); }

    ///
    Vec2d(const Point2d & p1, const Point2d & p2)
    { vx = p2.px - p1.px; vy = p2.py - p1.py; }
    
  ///
  Vec2d & operator= (const Vec2d & p2)
    { vx = p2.vx; vy = p2.vy; return *this; }

  ///
  double & X() { return vx; }
  ///
  double & Y() { return vy; }
  ///
  double X() const { return vx; }
  ///
  double Y() const { return vy; }

  ///
  double Length() const { return sqrt (vx * vx + vy * vy); }
  ///
  double Length2() const { return vx * vx + vy * vy; }

  void GetNormal (Vec2d & n) const { n.vx=-vy; n.vy=vx; }		// GH

  ///
  inline Vec2d & operator+= (const Vec2d & v2);
  ///
  inline Vec2d & operator-= (const Vec2d & v2);
  ///
  inline Vec2d & operator*= (double s);
  ///
  inline Vec2d & operator/= (double s);

  ///
  friend inline Vec2d operator- (const Point2d & p1, const Point2d & p2);
  ///
  friend inline Point2d operator- (const Point2d & p1, const Vec2d & v);
  ///
  friend inline Point2d operator+ (const Point2d & p1, const Vec2d & v);
  ///
  friend inline Vec2d operator- (const Vec2d & p1, const Vec2d & v);
  ///
  friend inline Vec2d operator+ (const Vec2d & p1, const Vec2d & v);
  ///
  friend inline Vec2d operator* (double scal, const Vec2d & v);

  ///
  friend double operator* (const Vec2d & v1, const Vec2d & v2)
    { return v1.X() * v2.X() + v1.Y() * v2.Y(); }


  ///
  friend double Cross (const Vec2d & v1, const Vec2d & v2)
    { return double(v1.X()) * double(v2.Y()) -
             double(v1.Y()) * double(v2.X()); }

  ///
  friend inline void PpSmV (const Point2d & p1, double s, const Vec2d & v, Point2d & p2);
  ///
  friend inline void PmP (const Point2d & p1, const Point2d & p2, Vec2d & v);

///						Angle in [0,2*PI)

  ///
  friend double Angle (const Vec2d & v);
  ///
  friend double FastAngle (const Vec2d & v);
  ///
  friend double Angle (const Vec2d & v1, const Vec2d & v2);
  ///
  friend double FastAngle (const Vec2d & v1, const Vec2d & v2);

  ///
  friend ostream & operator<<(ostream  & s, const Vec2d & v);
  };



///
class Line2d
  {
protected:
  ///
  Point2d p1, p2;

public:
  ///
  Line2d() : p1(), p2() { };
  ///
  Line2d(const Point2d & ap1, const Point2d & ap2)
    { p1 = ap1; p2 = ap2; }

  ///
  Line2d & operator= (const Line2d & l2)
    { p1 = l2.p1; p2 = l2.p2; return *this;}

  ///
  Point2d & P1() { return p1; }
  ///
  Point2d & P2() { return p2; }
  ///
  const Point2d & P1() const { return p1; }
  ///
  const Point2d & P2() const { return p2; }

  ///
  double XMax() const { return max2 (p1.X(), p2.X()); }
  ///
  double YMax() const { return max2 (p1.Y(), p2.Y()); }
  ///
  double XMin() const { return min2 (p1.X(), p2.X()); }
  ///
  double YMin() const { return min2 (p1.Y(), p2.Y()); }

  ///
  Vec2d Delta () const { return Vec2d (p2.X()-p1.X(), p2.Y()-p1.Y()); }
  ///
  double Length () const { return Delta().Length(); }
  ///
  double Length2 () const
        { return sqr (p1.X() - p2.X()) +
                 sqr (p1.Y() - p2.Y()); }

  void GetNormal (Line2d & n) const;					// GH
  Vec2d NormalDelta () const;						// GH

  /// square of the distance between two 2d-lines.
  friend double Dist2(const Line2d & g, const Line2d & h );		// GH

  ///
  friend Point2d CrossPoint (const Line2d & l1, const Line2d & l2);
    /// returns 1 iff parallel
    friend int CrossPointBarycentric (const Line2d & l1, const Line2d & l2,
				      double & lam1, double & lam2);

  ///
  friend int Parallel (const Line2d & l1, const Line2d & l2, double peps);
  ///
  friend int IsOnLine (const Line2d & l, const Point2d & p, double heps);
  ///
  friend int IsOnLongLine (const Line2d & l, const Point2d & p);
  ///
  friend int Hit (const Line2d & l1, const Line2d & l2, double heps);

  ///
  friend ostream & operator<<(ostream  & s, const Line2d & l);
  };


#ifdef NONE
///
class PLine2d
  {
protected:
  ///
  Point2d const * p1, *p2;

public:
  ///
  PLine2d() { };
  ///
  PLine2d(Point2d const * ap1, Point2d const * ap2)
    { p1 = ap1; p2 = ap2; }

  ///
  PLine2d & operator= (const PLine2d & l2)
    { p1 = l2.p1; p2 = l2.p2; return *this;}

  ///
  const Point2d *& P1() { return p1; }
  ///
  const Point2d *& P2() { return p2; }
  ///
  const Point2d & P1() const { return *p1; }
  ///
  const Point2d & P2() const { return *p2; }

  ///
  double XMax() const { return max2 (p1->X(), p2->X()); }
  ///
  double YMax() const { return max2 (p1->Y(), p2->Y()); }
  ///
  double XMin() const { return min2 (p1->X(), p2->X()); }
  ///
  double YMin() const { return min2 (p1->Y(), p2->Y()); }


  ///
  Vec2d Delta () const { return Vec2d (p2->X()-p1->X(), p2->Y()-p1->Y()); }
  ///
  double Length () const { return Delta().Length(); }
  ///
  double Length2 () const
        { return sqr (p1->X() - p2->X()) +
                 sqr (p1->Y() - p2->Y()); }


    
  ///
  friend Point2d CrossPoint (const PLine2d & l1, const PLine2d & l2);
  ///
  friend int Parallel (const PLine2d & l1, const PLine2d & l2, double peps);
  ///
  friend int IsOnLine (const PLine2d & l, const Point2d & p, double heps);
  ///
  friend int IsOnLongLine (const PLine2d & l, const Point2d & p);
  ///
  friend int Hit (const PLine2d & l1, const Line2d & l2, double heps);

  ///
  friend ostream & operator<<(ostream  & s, const Line2d & l);
  };



///
class ILINE
  {
  ///
  INDEX i[2];

  public:
  ///
  ILINE() {};
  ///
  ILINE(INDEX i1, INDEX i2) { i[0] = i1; i[1] = i2; }
  ///
  ILINE(const ILINE & l) { i[0] = l.i[0]; i[1] = l.i[1]; }

  ///
  ILINE & operator= (const ILINE & l)
    { i[0] = l.i[0]; i[1] = l.i[1]; return *this; }

  ///
  const INDEX & I(int ai) const { return i[ai-1]; }
  ///
  const INDEX & X() const { return i[0]; }
  ///
  const INDEX & Y() const { return i[1]; }
  ///
  const INDEX & I1() const { return i[0]; }
  ///
  const INDEX & I2() const { return i[1]; }

  ///
  INDEX & I(int ai) { return i[ai-1]; }
  ///
  INDEX & X() { return i[0]; }
  ///
  INDEX & Y() { return i[1]; }
  ///
  INDEX & I1() { return i[0]; }
  ///
  INDEX & I2() { return i[1]; }
  };




///
class TRIANGLE2D
  {
private:
  ///
  Point2d p1, p2, p3;

public:
  ///
  TRIANGLE2D() { };
  ///
  TRIANGLE2D (const Point2d & ap1, const Point2d & ap2,
              const Point2d & ap3)
    { p1 = ap1; p2 = ap2; p3 = ap3;}

  ///
  TRIANGLE2D & operator= (const TRIANGLE2D & t2)
    { p1 = t2.p1; p2 = t2.p2; p3 = t2.p3; return *this; }

  ///
  Point2d & P1() { return p1; }
  ///
  Point2d & P2() { return p2; }
  ///
  Point2d & P3() { return p3; }
  ///
  const Point2d & P1() const { return p1; }
  ///
  const Point2d & P2() const { return p2; }
  ///
  const Point2d & P3() const { return p3; }

  ///
  double XMax() const { return max3 (p1.X(), p2.X(), p3.X()); }
  ///
  double YMax() const { return max3 (p1.Y(), p2.Y(), p3.Y()); }
  ///
  double XMin() const { return min3 (p1.X(), p2.X(), p3.X()); }
  ///
  double YMin() const { return min3 (p1.Y(), p2.Y(), p3.Y()); }

  ///
  inline Point2d Center () const
   { return Point2d( (p1.X()+p2.X()+p3.X())/3, (p1.Y()+p2.Y()+p3.Y())/3); }

  ///
  int Regular() const;
  /// 
  int CW () const;
  ///
  int CCW () const;

  ///
  int IsOn (const Point2d & p) const;
  ///
  int IsIn (const Point2d & p) const;
  ///
  friend ostream & operator<<(ostream  & s, const TRIANGLE2D & t);
  };


///
class PTRIANGLE2D
  {
private:
  ///
  Point2d const *p1, *p2, *p3;

public:
  ///
  PTRIANGLE2D() { };
  ///
  PTRIANGLE2D (const Point2d * ap1, const Point2d * ap2,
              const Point2d * ap3)
    { p1 = ap1; p2 = ap2; p3 = ap3;}

  ///
  PTRIANGLE2D & operator= (const PTRIANGLE2D & t2)
    { p1 = t2.p1; p2 = t2.p2; p3 = t2.p3; return *this; }

  ///
  const Point2d *& P1() { return p1; }
  ///
  const Point2d *& P2() { return p2; }
  ///
  const Point2d *& P3() { return p3; }
  ///
  const Point2d * P1() const { return p1; }
  ///
  const Point2d * P2() const { return p2; }
  ///
  const Point2d * P3() const { return p3; }

  ///
  double XMax() const { return max3 (p1->X(), p2->X(), p3->X()); }
  ///
  double YMax() const { return max3 (p1->Y(), p2->Y(), p3->Y()); }
  ///
  double XMin() const { return min3 (p1->X(), p2->X(), p3->X()); }
  ///
  double YMin() const { return min3 (p1->Y(), p2->Y(), p3->Y()); }

  ///
  Point2d Center () const
   { return Point2d( (p1->X()+p2->X()+p3->X())/3, (p1->Y()+p2->Y()+p3->Y())/3); }


  ///
  int Regular() const;
  ///
  int CW () const;
  ///
  int CCW () const;

  ///
  int IsOn (const Point2d & p) const;
  ///
  int IsIn (const Point2d & p) const;
  ///
  friend ostream & operator<<(ostream & s, const PTRIANGLE2D & t);
  };
#endif



class Polygon2d
{
protected:
  ARRAY<Point2d> points;
  
public:
  Polygon2d ();
  ~Polygon2d ();
  
  void AddPoint (const Point2d & p);
  int GetNP() const { return points.Size(); }
  void GetPoint (int i, Point2d & p) const
  { p = points.Get(i); }
  void GetLine (int i, Point2d & p1, Point2d & p2) const
  { p1 = points.Get(i); p2 = points.Get(i%points.Size()+1); }

  double Area () const { return fabs (HArea()); }
  int CW () const { return HArea() > 0; }
  int CCW () const { return HArea() < 0; }
  
  int IsOn (const Point2d & p) const;
  int IsIn (const Point2d & p) const;

  int IsConvex () const;

  int IsStarPoint (const Point2d & p) const;
  Point2d Center() const;
  Point2d EqualAreaPoint () const;
private:
  double HArea () const;
};


/** Cheap approximation to atan2.
  A monotone function of atan2(x,y) is computed.
 */
extern double Fastatan2 (double x, double y);


inline Vec2d & Vec2d :: operator+= (const Vec2d & v2)
  {
  vx += v2.vx;
  vy += v2.vy;
  return *this;
  }

inline Vec2d & Vec2d :: operator-= (const Vec2d & v2)
  {
  vx -= v2.vx;
  vy -= v2.vy;
  return *this;
  }

inline Vec2d & Vec2d :: operator*= (double s)
  {
  vx *= s;
  vy *= s;
  return *this;
  }


inline Vec2d & Vec2d :: operator/= (double s)
{
  if (s != 0)
    {
      vx /= s;
      vy /= s;
    }
  else
    {
      MyError ("Vec2d::operator /=: Division by zero");
    }
  return *this;
}



inline Vec2d operator- (const Point2d & p1, const Point2d & p2)
  {
  return Vec2d (p1.X() - p2.X(), p1.Y() - p2.Y());
  }


inline Point2d operator- (const Point2d & p1, const Vec2d & v)
  {
  return Point2d (p1.X() - v.X(), p1.Y() - v.Y());
  }


inline Point2d operator+ (const Point2d & p1, const Vec2d & v)
  {
  return Point2d (p1.X() + v.X(), p1.Y() + v.Y());
  }


inline Point2d Center (const Point2d & p1, const Point2d & p2)
  {
  return Point2d ((p1.X() + p2.X()) / 2, (p1.Y() + p2.Y()) / 2);
  }


inline Vec2d operator- (const Vec2d & v1, const Vec2d & v2)
  {
  return Vec2d (v1.X() - v2.X(), v1.Y() - v2.Y());
  }


inline Vec2d operator+ (const Vec2d & v1, const Vec2d & v2)
  {
  return Vec2d (v1.X() + v2.X(), v1.Y() + v2.Y());
  }


inline Vec2d operator* (double scal, const Vec2d & v)
  {
  return Vec2d (scal * v.X(), scal * v.Y());
  }


inline void PpSmV (const Point2d & p1, double s,
                   const Vec2d & v, Point2d & p2)
  {
  p2.X() = p1.X() + s * v.X();
  p2.Y() = p1.Y() + s * v.Y();
  }


inline void PmP (const Point2d & p1, const Point2d & p2, Vec2d & v)
  {
  v.X() = p1.X() - p2.X();
  v.Y() = p1.Y() - p2.Y();
  }





#ifdef none
inline int TRIANGLE2D :: Regular() const
    {
    return fabs(Cross ( p2 - p1, p3 - p2)) > EPSGEOM;
    }


inline int TRIANGLE2D :: CW () const
    {
    return Cross ( p2 - p1, p3 - p2) < 0;
    }


inline int TRIANGLE2D :: CCW () const
    {
    return Cross ( p2 - p1, p3 - p2) > 0;
    }




inline int PTRIANGLE2D :: Regular() const
    {
    return fabs(Cross ( *p2 - *p1, *p3 - *p2)) > EPSGEOM;
    }


inline int PTRIANGLE2D :: CW () const
    {
    return Cross ( *p2 - *p1, *p3 - *p2) < 0;
    }


inline int PTRIANGLE2D :: CCW () const
    {
    return Cross ( *p2 - *p1, *p3 - *p2) > 0;
    }


#endif


///
class Mat2d
{
protected:
  ///
  double coeff[4];

public:
  ///
  Mat2d() { coeff[0] = coeff[1] = coeff[2] = coeff[3] = 0; }
  ///
  Mat2d(double a11, double a12, double a21, double a22)
    { coeff[0] = a11; coeff[1] = a12; coeff[2] = a21; coeff[3] = a22; }
  ///
  Mat2d(const Mat2d & m2)
    { for (int i = 0; i < 4; i++) coeff[i] = m2.Get(i); }

  ///
  double & Elem (INDEX i, INDEX j) { return coeff[2*(i-1)+j-1]; }
  ///
  double & Elem (INDEX i) {return coeff[i]; }
  ///
  double Get (INDEX i, INDEX j) const { return coeff[2*(i-1)+j-1]; }
  ///
  double Get (INDEX i) const {return coeff[i]; }

  ///  
  double Det () const { return coeff[0] * coeff[3] - coeff[1] * coeff[2]; }

  ///
  void Mult (const Vec2d & v, Vec2d & prod) const;
  ///
  void MultTrans (const Vec2d & v , Vec2d & prod) const;
  ///
  void Solve (const Vec2d & rhs, Vec2d & x) const;
  /// Solves mat * x = rhs, but using a positive definite matrix instead of mat
  void SolvePositiveDefinite (const Vec2d & rhs, Vec2d & x) const;
  /// add a term \alpha * v * v^T
  void AddDiadicProduct (double alpha, Vec2d & v);
};



inline void Mat2d :: Mult (const Vec2d & v, Vec2d & prod) const
{
  prod.X() = coeff[0] * v.X() + coeff[1] * v.Y();
  prod.Y() = coeff[2] * v.X() + coeff[3] * v.Y();
}


inline  void Mat2d :: MultTrans (const Vec2d & v, Vec2d & prod) const
{
  prod.X() = coeff[0] * v.X() + coeff[2] * v.Y();
  prod.Y() = coeff[1] * v.X() + coeff[3] * v.Y();
}



inline void Mat2d :: Solve (const Vec2d & rhs, Vec2d & x) const
{
  double det = Det();
  
  if (det == 0)
    MyError ("Mat2d::Solve: zero determinant");
  else
    {
      x.X() = (coeff[3] * rhs.X() - coeff[1] * rhs.Y()) / det;
      x.Y() = (-coeff[2] * rhs.X() + coeff[0] * rhs.Y()) / det;
    }
}


inline void Mat2d :: SolvePositiveDefinite (const Vec2d & rhs, Vec2d & x) const
{
  double a = max2(coeff[0], 1e-8);
  double b = coeff[1] / a;
  double c = coeff[2] / a;
  double d = max2(coeff[3] - a *b * c, 1e-8);

  x.X() = (rhs.X() - b * rhs.Y()) / a;
  x.Y() = rhs.Y() / d - c * x.X();
}


inline void Mat2d :: AddDiadicProduct (double alpha, Vec2d & v)
{
  coeff[0] += alpha * v.X() * v.X();
  coeff[1] += alpha * v.X() * v.Y();
  coeff[2] += alpha * v.Y() * v.X();
  coeff[3] += alpha * v.Y() * v.Y();
}



#endif
