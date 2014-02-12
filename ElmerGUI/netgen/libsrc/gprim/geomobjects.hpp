#ifndef FILE_OBJECTS
#define FILE_OBJECTS

/* *************************************************************************/
/* File:   geomobjects.hpp                                                 */
/* Author: Joachim Schoeberl                                               */
/* Date:   20. Jul. 02                                                     */
/* *************************************************************************/



template <int D> class Vec;
template <int D> class Point;


template <int D>
class Point
{

protected:
  double x[D];

public:
  Point () { ; }
  Point (double ax) { x[0] = ax; }
  Point (double ax, double ay) { x[0] = ax; x[1] = ay; }
  Point (double ax, double ay, double az)
  { x[0] = ax; x[1] = ay; x[2] = az; }
  Point (double ax, double ay, double az, double au)
  { x[0] = ax; x[1] = ay; x[2] = az; x[3] = au;}

  Point (const Point<D> & p2)
  { for (int i = 0; i < D; i++) x[i] = p2.x[i]; }

  explicit Point (const Vec<D> & v)
  { for (int i = 0; i < D; i++) x[i] = v(i); }


  Point & operator= (const Point<D> & p2)
  {
    for (int i = 0; i < D; i++) x[i] = p2.x[i]; 
    return *this;
  }

  Point & operator= (double val)
  {
    for (int i = 0; i < D; i++) x[i] = val;
    return *this;
  }

  double & operator() (int i) { return x[i]; }
  const double & operator() (int i) const { return x[i]; }

  operator const double* () const { return x; }
};





template <int D>
class Vec
{

protected:
  double x[D];

public:
  Vec () { ; } // for (int i = 0; i < D; i++) x[i] = 0; }
  Vec (double ax) { for (int i = 0; i < D; i++) x[i] = ax; }
  Vec (double ax, double ay) { x[0] = ax; x[1] = ay; }
  Vec (double ax, double ay, double az)
  { x[0] = ax; x[1] = ay; x[2] = az; }
  Vec (double ax, double ay, double az, double au)
  { x[0] = ax; x[1] = ay; x[2] = az; x[3] = au; }

  Vec (const Vec<D> & p2)
  { for (int i = 0; i < D; i++) x[i] = p2.x[i]; }

  explicit Vec (const Point<D> & p)
  { for (int i = 0; i < D; i++) x[i] = p(i); }

  Vec (const Vec<D> & p1, const Vec<D> & p2)
  { for(int i=0; i<D; i++) x[i] = p2(i)-p1(1); }
  


  Vec & operator= (const Vec<D> & p2)
  {
    for (int i = 0; i < D; i++) x[i] = p2.x[i]; 
    return *this;
  }

  Vec & operator= (double s)
  {
    for (int i = 0; i < D; i++) x[i] = s;
    return *this;
  }

  double & operator() (int i) { return x[i]; }
  const double & operator() (int i) const { return x[i]; }

  operator const double* () const { return x; }

  double Length () const
  {
    double l = 0;
    for (int i = 0; i < D; i++)
      l += x[i] * x[i];
    return sqrt (l);
  }

  double Length2 () const
  {
    double l = 0;
    for (int i = 0; i < D; i++)
      l += x[i] * x[i];
    return l;
  }

  const Vec<D> & Normalize ()
  {
    double l = Length();
    if (l != 0)
      for (int i = 0; i < D; i++)
	x[i] /= l;
    return *this;
  }

  Vec<D> GetNormal () const;
};





template <int H, int W=H>
class Mat
{

protected:
  double x[H*W];

public:
  Mat () { ; }
  Mat (const Mat & b)
  { for (int i = 0; i < H*W; i++) x[i] = b.x[i]; }
  
  Mat & operator= (double s)
  {
    for (int i = 0; i < H*W; i++) x[i] = s;
    return *this;
  }

  Mat & operator= (const Mat & b)
  {
    for (int i = 0; i < H*W; i++) x[i] = b.x[i]; 
    return *this;
  }

  double & operator() (int i, int j) { return x[i*W+j]; }
  const double & operator() (int i, int j) const { return x[i*W+j]; }
  double & operator() (int i) { return x[i]; }
  const double & operator() (int i) const { return x[i]; }

  Vec<H> Col (int i) const
  {
    Vec<H> hv; 
    for (int j = 0; j < H; j++)
      hv(j) = x[j*W+i];
    return hv; 
  }

  Vec<W> Row (int i) const
  {
    Vec<W> hv; 
    for (int j = 0; j < W; j++)
      hv(j) = x[i*W+j];
    return hv; 
  }

  void Solve (const Vec<H> & rhs, Vec<W> & sol) const
  {
    Mat<W,H> inv;
    CalcInverse (*this, inv);
    sol = inv * rhs;
  }
};




template <int D>
class Box
{
protected:
  Point<D> pmin, pmax;
public:
  Box () { ; }
  Box ( const Point<D> & p1, const Point<D> & p2)
  {
    for (int i = 0; i < D; i++)
      {
	pmin(i) = min2(p1(i), p2(i));
	pmax(i) = max2(p1(i), p2(i));
      }
  }

  enum EB_TYPE { EMPTY_BOX = 1 };
  Box ( EB_TYPE et ) 
  {
    pmin = Point<3> (1e99, 1e99, 1e99);
    pmax = Point<3> (-1e99, -1e99, -1e99);
  }

  const Point<D> & PMin () const { return pmin; }
  const Point<D> & PMax () const { return pmax; }
  
  void Set (const Point<D> & p)
  { pmin = pmax = p; }

  void Add (const Point<D> & p)
  { 
    for (int i = 0; i < D; i++)
      {
	if (p(i) < pmin(i)) pmin(i) = p(i);
	else if (p(i) > pmax(i)) pmax(i) = p(i);
      }
  }

  Point<D> Center () const 
  { 
    Point<D> c;
    for (int i = 0; i < D; i++)
      c(i) = 0.5 * (pmin(i)+pmax(i)); 
    return c;
  }
  double Diam () const { return Abs (pmax-pmin); }

  Point<D> GetPointNr (int nr) const
  {
    Point<D> p;
    for (int i = 0; i < D; i++)
      {
	p(i) = (nr & 1) ? pmax(i) : pmin(i);
	nr >>= 1;
      }
    return p;
  }


  bool Intersect (const Box<D> & box2) const
  {
    for (int i = 0; i < D; i++)
      if (pmin(i) > box2.pmax(i) ||
	  pmax(i) < box2.pmin(i)) return 0;
    return 1;
  }


  bool IsIn (const Point<D> & p) const
  {
    for (int i = 0; i < D; i++)
      if (p(i) < pmin(i) || p(i) > pmax(i)) return 0;
    return 1;
  }


  void Increase (double dist)
  {
    for (int i = 0; i < D; i++)
      {
	pmin(i) -= dist;
	pmax(i) += dist;
      }
  }
};




template <int D>
class BoxSphere : public Box<D>
{
protected:
  ///
  Point<D> c;
  ///
  double diam;
  ///
  double inner;
public:
  ///
  BoxSphere () { };
 ///
  BoxSphere (const Box<D> & box) 
  : Box<D> (box) 
  { 
    CalcDiamCenter();
  };

  ///
  BoxSphere ( Point<D> apmin, Point<D> apmax )
    : Box<D> (apmin, apmax)
  {
    CalcDiamCenter();
  }

  ///
  const Point<D> & Center () const { return c; }
  ///
  double Diam () const { return diam; }
  ///
  double Inner () const { return inner; }


  ///
  void GetSubBox (int nr, BoxSphere & sbox) const
  {
    for (int i = 0; i < D; i++)
      {
	if (nr & 1)
	  {
	    sbox.pmin(i) = c(i);
	    sbox.pmax(i) = this->pmax(i);
	  }
	else
	  {
	    sbox.pmin(i) = this->pmin(i);
	    sbox.pmax(i) = c(i);
	  }
	sbox.c(i) = 0.5 * (sbox.pmin(i) + sbox.pmax(i));
	nr >>= 1;
      }
    sbox.diam = 0.5 * diam;
    sbox.inner = 0.5 * inner;
  }


  ///
  void CalcDiamCenter ()
  {
    c = Box<D>::Center ();
    diam = Dist (this->pmin, this->pmax);

    inner = this->pmax(0) - this->pmin(0);
    for (int i = 1; i < D; i++)
      if (this->pmax(i) - this->pmin(i) < inner)
	inner = this->pmax(i) - this->pmin(i);
  }

};






#endif
