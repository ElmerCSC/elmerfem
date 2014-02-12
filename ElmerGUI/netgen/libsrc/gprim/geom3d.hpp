#ifndef FILE_GEOM3D
#define FILE_GEOM3D

/* *************************************************************************/
/* File:   geom3d.hh                                                       */
/* Author: Joachim Schoeberl                                               */
/* Date:   5. Aug. 95                                                      */
/* *************************************************************************/




extern void MyError (const char * ch);

class Point3d;
class Vec3d;

inline Vec3d operator- (const Point3d & p1, const Point3d & p2);
inline Point3d operator- (const Point3d & p1, const Vec3d & v);
inline Point3d operator+ (const Point3d & p1, const Vec3d & v);
Point3d & Add (double d, const Vec3d & v);
Point3d & Add2 (double d, const Vec3d & v,
			 double d2, const Vec3d & v2);
inline Point3d Center (const Point3d & p1, const Point3d & p2);
inline Point3d Center (const Point3d & p1, const Point3d & p2, const Point3d & p3);
inline Point3d Center (const Point3d & p1, const Point3d & p2, 
				const Point3d & p3, const Point3d & p4);
ostream & operator<<(ostream  & s, const Point3d & p);
inline Vec3d operator- (const Vec3d & p1, const Vec3d & v);
inline Vec3d operator+ (const Vec3d & p1, const Vec3d & v);
inline Vec3d operator* (double scal, const Vec3d & v);
inline double operator* (const Vec3d & v1, const Vec3d & v2);
inline Vec3d Cross (const Vec3d & v1, const Vec3d & v2);
inline void Cross (const Vec3d & v1, const Vec3d & v2, Vec3d & prod);
double Angle (const Vec3d & v);
double FastAngle (const Vec3d & v);
double Angle (const Vec3d & v1, const Vec3d & v2);
double FastAngle (const Vec3d & v1, const Vec3d & v2);
ostream & operator<<(ostream  & s, const Vec3d & v);
void Transpose (Vec3d & v1, Vec3d & v2, Vec3d & v3);
int SolveLinearSystem (const Vec3d & col1,
		       const Vec3d & col2,
		       const Vec3d & col3,
		       const Vec3d & rhs,
		       Vec3d & sol);
int SolveLinearSystemLS (const Vec3d & col1,
			 const Vec3d & col2,
			 const Vec2d & rhs,
			 Vec3d & sol);
int SolveLinearSystemLS2 (const Vec3d & col1,
			  const Vec3d & col2,
			  const Vec2d & rhs, 
			  Vec3d & sol,
			  double & x, double & y);
int PseudoInverse (const Vec3d & col1,
		   const Vec3d & col2,
		   Vec3d & inv1,
		   Vec3d & inv2);
double Determinant (const Vec3d & col1,
		    const Vec3d & col2,
		    const Vec3d & col3);

/// Point in R3
class Point3d
{
protected:
  ///
  double x[3];
  
public:
  ///
  Point3d () { x[0] = x[1] = x[2] = 0; }
  ///
  Point3d(double ax, double ay, double az) 
    { x[0] = ax; x[1] = ay; x[2] = az; }
  ///
  Point3d(double ax[3])
    { x[0] = ax[0]; x[1] = ax[1]; x[2] = ax[2]; }

  ///
  Point3d(const Point3d & p2) 
    { x[0] = p2.x[0]; x[1] = p2.x[1]; x[2] = p2.x[2]; }

  Point3d (const Point<3> & p2)
  {
    for (int i = 0; i < 3; i++)
      x[i] = p2(i);
  }
  
  ///
  Point3d & operator= (const Point3d & p2)
    { x[0] = p2.x[0]; x[1] = p2.x[1]; x[2] = p2.x[2]; return *this; }
  
  ///
  int operator== (const Point3d& p) const
    { return (x[0] == p.x[0] && x[1] == p.x[1] && x[2] == p.x[2]); }
  
  ///
  double & X() { return x[0]; }
  ///
  double & Y() { return x[1]; }
  ///
  double & Z() { return x[2]; }
  ///
  double X() const { return x[0]; }
  ///
  double Y() const { return x[1]; }
  ///
  double Z() const { return x[2]; }
  ///
  double & X(int i) { return x[i-1]; }
  ///
  double X(int i) const { return x[i-1]; }
  ///
  const Point3d & SetToMin (const Point3d & p2)
    {
      if (p2.x[0] < x[0]) x[0] = p2.x[0];
      if (p2.x[1] < x[1]) x[1] = p2.x[1];
      if (p2.x[2] < x[2]) x[2] = p2.x[2];
      return *this;
    }

  ///
  const Point3d & SetToMax (const Point3d & p2)
    {
      if (p2.x[0] > x[0]) x[0] = p2.x[0];
      if (p2.x[1] > x[1]) x[1] = p2.x[1];
      if (p2.x[2] > x[2]) x[2] = p2.x[2];
      return *this;
    }

  ///
  friend inline Vec3d operator- (const Point3d & p1, const Point3d & p2);
  ///
  friend inline Point3d operator- (const Point3d & p1, const Vec3d & v);
  ///
  friend inline Point3d operator+ (const Point3d & p1, const Vec3d & v);
  ///
  inline Point3d & operator+= (const Vec3d & v);
  inline Point3d & operator-= (const Vec3d & v);
  ///
  inline Point3d & Add (double d, const Vec3d & v);
  ///
  inline Point3d & Add2 (double d, const Vec3d & v,
			 double d2, const Vec3d & v2);
  ///
  friend inline double Dist (const Point3d & p1, const Point3d & p2)
    { return sqrt (  (p1.x[0]-p2.x[0]) * (p1.x[0]-p2.x[0]) + 
		     (p1.x[1]-p2.x[1]) * (p1.x[1]-p2.x[1]) +
                     (p1.x[2]-p2.x[2]) * (p1.x[2]-p2.x[2])); }
  ///
  inline friend double Dist2 (const Point3d & p1, const Point3d & p2)
    { return  (  (p1.x[0]-p2.x[0]) * (p1.x[0]-p2.x[0]) + 
		 (p1.x[1]-p2.x[1]) * (p1.x[1]-p2.x[1]) +
		 (p1.x[2]-p2.x[2]) * (p1.x[2]-p2.x[2])); }

  ///
  friend inline Point3d Center (const Point3d & p1, const Point3d & p2);
  ///
  friend inline Point3d Center (const Point3d & p1, const Point3d & p2, const Point3d & p3);
  ///
  friend inline Point3d Center (const Point3d & p1, const Point3d & p2, 
				const Point3d & p3, const Point3d & p4);
  ///
  friend ostream & operator<<(ostream  & s, const Point3d & p);
  
  ///
  friend class Vec3d;
  ///
  friend class Box3d;


  operator Point<3> () const
  {
    return Point<3> (x[0], x[1], x[2]);
  }
};


///
class Vec3d
{
protected:
  ///
  double x[3];

public:
  ///
  inline Vec3d() { x[0] = x[1] = x[2] = 0; }
  ///
  inline Vec3d(double ax, double ay, double az)
    { x[0] = ax; x[1] = ay; x[2] = az; }
  ///
  Vec3d(double ax[3])
    { x[0] = ax[0]; x[1] = ax[1]; x[2] = ax[2]; }
  ///
  inline Vec3d(const Vec3d & v2) 
    { x[0] = v2.x[0]; x[1] = v2.x[1]; x[2] = v2.x[2]; }
  ///
  inline Vec3d(const Point3d & p1, const Point3d & p2)
    { 
      x[0] = p2.x[0] - p1.x[0];
      x[1] = p2.x[1] - p1.x[1];
      x[2] = p2.x[2] - p1.x[2]; 
    }
  ///
  inline Vec3d(const Point3d & p1)
    { 
      x[0] = p1.x[0];
      x[1] = p1.x[1];
      x[2] = p1.x[2];
    }
  
  Vec3d (const Vec<3> & v2)
  {
    for (int i = 0; i < 3; i++)
      x[i] = v2(i);
  }

  operator Vec<3> () const
  {
    return Vec<3> (x[0], x[1], x[2]);
  }


  Vec3d & operator= (const Vec3d & v2)
    { x[0] = v2.x[0]; x[1] = v2.x[1]; x[2] = v2.x[2]; return *this; }
  ///
  Vec3d & operator= (double val)
    { x[0] = x[1] = x[2] = val; return *this; }
  ///
  double & X() { return x[0]; }
  ///
  double & Y() { return x[1]; }
  ///
  double & Z() { return x[2]; }
  ///
  double & X(int i) { return x[i-1]; }

  ///
  double X() const { return x[0]; }
  ///
  double Y() const { return x[1]; }
  ///
  double Z() const { return x[2]; }
  ///
  double X(int i) const { return x[i-1]; }

  ///
  double Length() const 
    { return sqrt (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]); }
  ///
  double Length2() const 
    { return x[0] * x[0] + x[1] * x[1] + x[2] * x[2]; }

  ///
  Vec3d & operator+= (const Vec3d & v2);
  ///
  Vec3d & operator-= (const Vec3d & v2);
  ///
  Vec3d & operator*= (double s);
  ///
  Vec3d & operator/= (double s);
  ///
  inline Vec3d & Add (double d, const Vec3d & v);
  ///
  inline Vec3d & Add2 (double d, const Vec3d & v,
			 double d2, const Vec3d & v2);

  ///
  friend inline Vec3d operator- (const Point3d & p1, const Point3d & p2);
  ///
  friend inline Point3d operator- (const Point3d & p1, const Vec3d & v);
  ///
  friend inline Point3d operator+ (const Point3d & p1, const Vec3d & v);
  ///
  friend inline Vec3d operator- (const Vec3d & p1, const Vec3d & v);
  ///
  friend inline Vec3d operator+ (const Vec3d & p1, const Vec3d & v);
  ///
  friend inline Vec3d operator* (double scal, const Vec3d & v);

  ///
  friend inline double operator* (const Vec3d & v1, const Vec3d & v2);
  ///
  friend inline Vec3d Cross (const Vec3d & v1, const Vec3d & v2);
  ///
  friend inline void Cross (const Vec3d & v1, const Vec3d & v2, Vec3d & prod);

  /// Returns one normal-vector to n
  void GetNormal (Vec3d & n) const;
  ///
  friend double Angle (const Vec3d & v);
  ///
  friend double FastAngle (const Vec3d & v);
  ///
  friend double Angle (const Vec3d & v1, const Vec3d & v2);
  ///
  friend double FastAngle (const Vec3d & v1, const Vec3d & v2);

  void Normalize() 
  {
    double len = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    if (len == 0) return;
    len = sqrt (len);
    x[0] /= len; x[1] /= len; x[2] /= len;
  }

  ///
  friend ostream & operator<<(ostream  & s, const Vec3d & v);

  ///
  friend class Point3d;
  friend void Transpose (Vec3d & v1, Vec3d & v2, Vec3d & v3);
  friend int SolveLinearSystem (const Vec3d & col1,
				const Vec3d & col2,
				const Vec3d & col3,
				const Vec3d & rhs,
				Vec3d & sol);
  friend int SolveLinearSystemLS (const Vec3d & col1,
				  const Vec3d & col2,
				  const Vec2d & rhs,
				  Vec3d & sol);
  friend int SolveLinearSystemLS2 (const Vec3d & col1,
				   const Vec3d & col2,
				   const Vec2d & rhs, 
				   Vec3d & sol,
				   double & x, double & y);
  friend int PseudoInverse (const Vec3d & col1,
			    const Vec3d & col2,
			    Vec3d & inv1,
			    Vec3d & inv2);
  friend double Determinant (const Vec3d & col1,
			     const Vec3d & col2,
			     const Vec3d & col3);
};



class QuadraticFunction3d
{
  double c0, cx, cy, cz;
  double cxx, cyy, czz, cxy, cxz, cyz;

public:
  QuadraticFunction3d (const Point3d & p, const Vec3d & v);
  double Eval (const Point3d & p)
    {
      return 
	c0 
	+ p.X() * (cx + cxx * p.X() + cxy * p.Y() + cxz * p.Z())
	+ p.Y() * (cy + cyy * p.Y() + cyz * p.Z())
	+ p.Z() * (cz + czz * p.Z());
    }
};



inline Point3d Center (const Point3d & p1, const Point3d & p2)
{
  return Point3d (0.5 * (p1.x[0] + p2.x[0]),
		  0.5 * (p1.x[1] + p2.x[1]),
		  0.5 * (p1.x[2] + p2.x[2]));
}


inline Point3d Center (const Point3d & p1, const Point3d & p2,
                       const Point3d & p3)
{
  return Point3d (1.0/3.0 * (p1.x[0] + p2.x[0] + p3.x[0]),
		  1.0/3.0 * (p1.x[1] + p2.x[1] + p3.x[1]),
		  1.0/3.0 * (p1.x[2] + p2.x[2] + p3.x[2]));
}

inline Point3d Center (const Point3d & p1, const Point3d & p2,
                       const Point3d & p3, const Point3d & p4)
{
  return Point3d (0.25 * (p1.x[0] + p2.x[0] + p3.x[0] + p4.x[0]),
		  0.25 * (p1.x[1] + p2.x[1] + p3.x[1] + p4.x[1]),
		  0.25 * (p1.x[2] + p2.x[2] + p3.x[2] + p4.x[2]));
}



inline Vec3d & Vec3d :: operator+= (const Vec3d & v2)
{
  x[0] += v2.x[0];
  x[1] += v2.x[1];
  x[2] += v2.x[2];
  return *this;
}

inline Vec3d & Vec3d :: operator-= (const Vec3d & v2)
{
  x[0] -= v2.x[0];
  x[1] -= v2.x[1];
  x[2] -= v2.x[2];
  return *this;
}


inline Vec3d & Vec3d :: operator*= (double s)
{
  x[0] *= s;
  x[1] *= s;
  x[2] *= s;
  return *this;
}


inline Vec3d & Vec3d :: operator/= (double s)
{
  if (s != 0)
    {
      x[0] /= s;
      x[1] /= s;
      x[2] /= s;
    }
#ifdef DEBUG
  else
    {
      cerr << "Vec div by 0, v = " << (*this) << endl;
      //      MyError ("Vec3d::operator /=: Divisioin by zero");
    }
#endif
  return *this;
}

inline Vec3d & Vec3d::Add (double d, const Vec3d & v)
{
  x[0] += d * v.x[0]; 
  x[1] += d * v.x[1]; 
  x[2] += d * v.x[2];
  return *this;
}

inline Vec3d & Vec3d::Add2 (double d, const Vec3d & v,
			    double d2, const Vec3d & v2)
{
  x[0] += d * v.x[0] + d2 * v2.x[0]; 
  x[1] += d * v.x[1] + d2 * v2.x[1]; 
  x[2] += d * v.x[2] + d2 * v2.x[2]; 
  return *this;
}








inline Vec3d operator- (const Point3d & p1, const Point3d & p2)
{
  return Vec3d (p1.x[0] - p2.x[0], p1.x[1] - p2.x[1],p1.x[2] - p2.x[2]);
}


inline Point3d operator- (const Point3d & p1, const Vec3d & v)
{
  return Point3d (p1.x[0] - v.x[0], p1.x[1] - v.x[1],p1.x[2] - v.x[2]);
}


inline Point3d operator+ (const Point3d & p1, const Vec3d & v)
{
  return Point3d (p1.x[0] + v.x[0], p1.x[1] + v.x[1],p1.x[2] + v.x[2]);
}

inline Point3d & Point3d::operator+= (const Vec3d & v) 
{
  x[0] += v.x[0]; 
  x[1] += v.x[1]; 
  x[2] += v.x[2];
  return *this;
}

inline Point3d & Point3d::operator-= (const Vec3d & v) 
{
  x[0] -= v.x[0]; 
  x[1] -= v.x[1]; 
  x[2] -= v.x[2];
  return *this;
}

inline Point3d & Point3d::Add (double d, const Vec3d & v)
{
  x[0] += d * v.x[0]; 
  x[1] += d * v.x[1]; 
  x[2] += d * v.x[2];
  return *this;
}

inline Point3d & Point3d::Add2 (double d, const Vec3d & v,
				double d2, const Vec3d & v2)
{
  x[0] += d * v.x[0] + d2 * v2.x[0]; 
  x[1] += d * v.x[1] + d2 * v2.x[1]; 
  x[2] += d * v.x[2] + d2 * v2.x[2]; 
  return *this;
}


inline Vec3d operator- (const Vec3d & v1, const Vec3d & v2)
{
  return Vec3d (v1.x[0] - v2.x[0], v1.x[1] - v2.x[1],v1.x[2] - v2.x[2]);
}


inline Vec3d operator+ (const Vec3d & v1, const Vec3d & v2)
{
  return Vec3d (v1.x[0] + v2.x[0], v1.x[1] + v2.x[1],v1.x[2] + v2.x[2]);
}


inline Vec3d operator* (double scal, const Vec3d & v)
{
  return Vec3d (scal * v.x[0], scal * v.x[1], scal * v.x[2]);
}



inline double operator* (const Vec3d & v1, const Vec3d & v2)
{
  return v1.x[0] * v2.x[0] + v1.x[1] * v2.x[1] + v1.x[2] * v2.x[2];
}



inline Vec3d Cross (const Vec3d & v1, const Vec3d & v2)
{
  return Vec3d
    ( v1.x[1] * v2.x[2] - v1.x[2] * v2.x[1],
      v1.x[2] * v2.x[0] - v1.x[0] * v2.x[2],
      v1.x[0] * v2.x[1] - v1.x[1] * v2.x[0]);
}

inline void Cross (const Vec3d & v1, const Vec3d & v2, Vec3d & prod)
{
  prod.x[0] = v1.x[1] * v2.x[2] - v1.x[2] * v2.x[1];
  prod.x[1] = v1.x[2] * v2.x[0] - v1.x[0] * v2.x[2];
  prod.x[2] = v1.x[0] * v2.x[1] - v1.x[1] * v2.x[0];
}

inline double Determinant (const Vec3d & col1,
			   const Vec3d & col2,
			   const Vec3d & col3)
{
  return
    col1.x[0] * ( col2.x[1] * col3.x[2] - col2.x[2] * col3.x[1]) +
    col1.x[1] * ( col2.x[2] * col3.x[0] - col2.x[0] * col3.x[2]) +
    col1.x[2] * ( col2.x[0] * col3.x[1] - col2.x[1] * col3.x[0]);
}


///
class Box3d
{
protected:
  ///
  double minx[3], maxx[3];

public:
  ///
  Box3d () { };
  ///
  Box3d ( double aminx, double amaxx,
          double aminy, double amaxy,
          double aminz, double amaxz );
  ///
  Box3d ( const Box3d & b2 );
  ///
  Box3d (const Point3d& p1, const Point3d& p2);
  ///
  Box3d (const Box<3> & b2);
  ///
  double MinX () const { return minx[0]; }
  ///
  double MaxX () const { return maxx[0]; }
  ///
  double MinY () const { return minx[1]; }
  ///
  double MaxY () const { return maxx[1]; }
  ///
  double MinZ () const { return minx[2]; }
  ///
  double MaxZ () const { return maxx[2]; }

  ///
  double Mini (int i) const { return minx[i-1]; }
  ///
  double Maxi (int i) const { return maxx[i-1]; }

  ///
  Point3d PMin () const { return Point3d(minx[0], minx[1], minx[2]); }
  ///
  Point3d PMax () const { return Point3d(maxx[0], maxx[1], maxx[2]); }

  ///
  void GetPointNr (int i, Point3d & point) const;
  /// increase Box at each side with dist 
  void Increase (double dist);
  /// increase Box by factor rel
  void IncreaseRel (double rel);
  /// return 1 if closures are intersecting
  int Intersect (const Box3d & box2) const
  {
    if (minx[0] > box2.maxx[0] || maxx[0] < box2.minx[0] ||
	minx[1] > box2.maxx[1] || maxx[1] < box2.minx[1] ||
	minx[2] > box2.maxx[2] || maxx[2] < box2.minx[2])
      return 0;
    return 1;
  }
  /// return 1 if point p in closure
  int IsIn (const Point3d & p) const
    {
      if (minx[0] <= p.x[0] && maxx[0] >= p.x[0] &&
	  minx[1] <= p.x[1] && maxx[1] >= p.x[1] &&
	  minx[2] <= p.x[2] && maxx[2] >= p.x[2])
	return 1;
      return 0;
    }
  ///
  inline void SetPoint (const Point3d & p)
    {
      minx[0] = maxx[0] = p.X();
      minx[1] = maxx[1] = p.Y();
      minx[2] = maxx[2] = p.Z();    
    }

  ///
  inline void AddPoint (const Point3d & p)
    {
      if (p.x[0] < minx[0]) minx[0] = p.x[0];
      if (p.x[0] > maxx[0]) maxx[0] = p.x[0];
      if (p.x[1] < minx[1]) minx[1] = p.x[1];
      if (p.x[1] > maxx[1]) maxx[1] = p.x[1];
      if (p.x[2] < minx[2]) minx[2] = p.x[2];
      if (p.x[2] > maxx[2]) maxx[2] = p.x[2];
    }

  ///
  const Box3d& operator+=(const Box3d& b);

  ///
  Point3d MaxCoords() const;
  ///
  Point3d MinCoords() const;

  /// Make a negative sized box;
  //  void CreateNegMinMaxBox();
  
  ///
  Point3d CalcCenter () const { return Point3d(0.5*(minx[0] + maxx[0]),
					       0.5*(minx[1] + maxx[1]),
					       0.5*(minx[2] + maxx[2])); }
  ///
  double CalcDiam () const { return sqrt(sqr(maxx[0]-minx[0])+
					 sqr(maxx[1]-minx[1])+
					 sqr(maxx[2]-minx[2])); }

  ///
  void WriteData(ofstream& fout) const;
  ///
  void ReadData(ifstream& fin);
};


class Box3dSphere : public Box3d
{
protected:
  ///
  double diam, inner;
  ///
  Point3d c;
public:
  ///
  Box3dSphere () { };
  ///
  Box3dSphere ( double aminx, double amaxx,
		double aminy, double amaxy,
		double aminz, double amaxz);
  ///
  const Point3d & Center () const { return c; }

  ///
  double Diam () const { return diam; }
  ///
  double Inner () const { return inner; }
  ///
  void GetSubBox (int i, Box3dSphere & sbox) const;

  // private:
  ///
  void CalcDiamCenter ();
};




///
class referencetransform
{
  ///
  Vec3d ex, ey, ez;
  ///
  Vec3d exh, eyh, ezh;
  ///
  Vec3d ex_h, ey_h, ez_h;
  ///
  Point3d rp;
  ///
  double h;

public:

  ///
  void Set (const Point3d & p1, const Point3d & p2,
            const Point3d & p3, double ah);

  ///
  void ToPlain (const Point3d & p, Point3d & pp) const;
  ///
  void ToPlain (const ARRAY<Point3d> & p, ARRAY<Point3d> & pp) const;
  ///
  void FromPlain (const Point3d & pp, Point3d & p) const;
};



#endif
