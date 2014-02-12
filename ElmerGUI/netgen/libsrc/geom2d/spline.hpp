#ifndef FILE_SPLINE_HPP
#define FILE_SPLINE_HPP

/**************************************************************************/
/* File:   spline.hpp                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   24. Jul. 96                                                    */
/**************************************************************************/


void CalcPartition (double l, double h, double r1, double r2,
		    double ra, double elto0, ARRAY<double> & points);

/*
  Spline curves for 2D mesh generation
*/


/// Geometry point
template < int D >
class GeomPoint : public Point<D>
{
public:
  /// refinement to point
  double refatpoint;
  bool hpref;

  GeomPoint ()
  { ; }

  ///
  GeomPoint (double ax, double ay, double aref = 1, bool ahpref=false)
  : Point<D> (ax, ay), refatpoint(aref), hpref(ahpref) { ; }
  GeomPoint (double ax, double ay, double az, double aref, bool ahpref=false)
  : Point<D> (ax, ay, az), refatpoint(aref), hpref(ahpref) { ; }
  GeomPoint (const Point<D> & ap, double aref = 1, bool ahpref=false)
  : Point<D>(ap), refatpoint(aref), hpref(ahpref) { ; }
};



/// base class for 2d - segment
template < int D >
class SplineSeg
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
  virtual Point<D> GetPoint (double t) const = 0;
  /// returns a (not necessarily uniform) tangent vector for 0 <= t <= 1
  virtual Vec<D> GetTangent (const double t) const
  { cerr << "GetTangent not implemented for spline base-class"  << endl; Vec<D> dummy; return dummy;}
  virtual void GetDerivatives (const double t, 
			       Point<D> & point,
			       Vec<D> & first,
			       Vec<D> & second) const {;}
  /// partitionizes curve
  void Partition (double h, double elto0,
		  Mesh & mesh, Point3dTree & searchtree, int segnr) const;
  /// returns initial point on curve
  virtual const GeomPoint<D> & StartPI () const = 0;
  /// returns terminal point on curve
  virtual const GeomPoint<D> & EndPI () const = 0;
  /** writes curve description for fepp:
      for implicitly given quadratic curves, the 6 coefficients of
      the polynomial
      $$ a x^2 + b y^2 + c x y + d x + e y + f = 0 $$
      are written to ost */
  void PrintCoeff (ostream & ost) const;

  virtual void GetCoeff (Vector & coeffs) const = 0;

  virtual void GetPoints (int n, ARRAY<Point<D> > & points);

  /** calculates (2D) lineintersections:
      for lines $$ a x + b y + c = 0 $$ the interecting points are calculated
      and stored in points */
  virtual void LineIntersections (const double a, const double b, const double c,
				  ARRAY < Point<D> > & points, const double eps) const
  {points.SetSize(0);}

  virtual double MaxCurvature(void) const = 0;

  virtual string GetType(void) const {return "splinebase";}

  virtual void Project (const Point<D> point, Point<D> & point_on_curve, double & t) const
  { cerr << "Project not implemented for spline base-class" << endl;}

  virtual void GetRawData (ARRAY<double> & data) const
  { cerr << "GetRawData not implemented for spline base-class" << endl;}

};


/// Straight line form p1 to p2
template< int D >
class LineSeg : public SplineSeg<D>
{
  ///
  GeomPoint<D> p1, p2;
  //const GeomPoint<D> &p1, &p2;
public:
  ///
  LineSeg (const GeomPoint<D> & ap1, const GeomPoint<D> & ap2);
  ///
  virtual double Length () const;
  ///
  virtual Point<D> GetPoint (double t) const;
  ///
  virtual Vec<D> GetTangent (const double t) const;

  
  virtual void GetDerivatives (const double t, 
			       Point<D> & point,
			       Vec<D> & first,
			       Vec<D> & second) const;
  ///
  virtual const GeomPoint<D> & StartPI () const { return p1; };
  ///
  virtual const GeomPoint<D> & EndPI () const { return p2; }
  ///
  virtual void GetCoeff (Vector & coeffs) const;

  virtual string GetType(void) const {return "line";}

  virtual void LineIntersections (const double a, const double b, const double c,
				  ARRAY < Point<D> > & points, const double eps) const;

  virtual double MaxCurvature(void) const {return 0;}

  virtual void Project (const Point<D> point, Point<D> & point_on_curve, double & t) const;

  virtual void GetRawData (ARRAY<double> & data) const;
};


/// curve given by a rational, quadratic spline (including ellipses)
template< int D >
class SplineSeg3 : public SplineSeg<D>
{
  ///
  GeomPoint<D> p1, p2, p3;
  //const GeomPoint<D> &p1, &p2, &p3;

  mutable double proj_latest_t;
public:
  ///
  SplineSeg3 (const GeomPoint<D> & ap1, 
	      const GeomPoint<D> & ap2, 
	      const GeomPoint<D> & ap3);
  ///
  virtual Point<D> GetPoint (double t) const;
  ///
  virtual Vec<D> GetTangent (const double t) const;

  
  virtual void GetDerivatives (const double t, 
			       Point<D> & point,
			       Vec<D> & first,
			       Vec<D> & second) const;
  ///
  virtual const GeomPoint<D> & StartPI () const { return p1; };
  ///
  virtual const GeomPoint<D> & EndPI () const { return p3; }
  ///
  virtual void GetCoeff (Vector & coeffs) const;

  virtual string GetType(void) const {return "spline3";}

  const GeomPoint<D> & TangentPoint (void) const { return p2; }

  virtual void LineIntersections (const double a, const double b, const double c,
				  ARRAY < Point<D> > & points, const double eps) const;

  virtual double MaxCurvature(void) const;

  virtual void Project (const Point<D> point, Point<D> & point_on_curve, double & t) const;

  virtual void GetRawData (ARRAY<double> & data) const;
};


// Gundolf Haase  8/26/97
/// A circle
template < int D >
class CircleSeg : public SplineSeg<D>
{
  ///
private:
  GeomPoint<D>	p1, p2, p3;
  //const GeomPoint<D>	&p1, &p2, &p3;
  Point<D>		pm;
  double		radius, w1,w3;
public:
  ///
  CircleSeg (const GeomPoint<D> & ap1, 
	     const GeomPoint<D> & ap2, 
	     const GeomPoint<D> & ap3);
  ///
  virtual Point<D> GetPoint (double t) const;
  ///
  virtual const GeomPoint<D> & StartPI () const { return p1; }
  ///
  virtual const GeomPoint<D> & EndPI () const { return p3; }
  ///
  virtual void GetCoeff (Vector & coeffs) const;
  ///
  double Radius() const { return radius; }
  ///
  double StartAngle() const { return w1; }
  ///
  double EndAngle() const { return w3; }
  ///
  const Point<D> & MidPoint(void) const {return pm; }

  virtual string GetType(void) const {return "circle";}

  virtual void LineIntersections (const double a, const double b, const double c,
				  ARRAY < Point<D> > & points, const double eps) const;

  virtual double MaxCurvature(void) const {return 1./radius;}
};






/// 
template<int D>
class DiscretePointsSeg : public SplineSeg<D>
{
  ARRAY<Point<D> > pts;
  GeomPoint<D> p1, p2;
public:
  ///
  DiscretePointsSeg (const ARRAY<Point<D> > & apts);
  ///
  virtual ~DiscretePointsSeg ();
  ///
  virtual Point<D> GetPoint (double t) const;
  ///
  virtual const GeomPoint<D> & StartPI () const { return p1; };
  ///
  virtual const GeomPoint<D> & EndPI () const { return p2; }
  ///
  virtual void GetCoeff (Vector & coeffs) const {;}

  virtual double MaxCurvature(void) const {return 1;}
};







// calculates length of spline-curve
template<int D>
double SplineSeg<D> :: Length () const
{
  Point<D> p, pold;

  int i, n = 100;
  double dt = 1.0 / n;

  pold = GetPoint (0);

  double l = 0;
  for (i = 1; i <= n; i++)
    {
      p = GetPoint (i * dt);
      l += Dist (p, pold);
      pold = p;
    }
  return l;
}



// partitionizes spline curve
template<int D>
void SplineSeg<D> :: Partition (double h, double elto0,
				Mesh & mesh, Point3dTree & searchtree, int segnr) const
{
  int i, j;
  double l, r1, r2, ra;
  double lold, dt, frac;
  int n = 100;
  Point<D> p, pold, mark, oldmark;
  ARRAY<double> curvepoints;
  double edgelength, edgelengthold;
  l = Length();

  r1 = StartPI().refatpoint;
  r2 = EndPI().refatpoint;
  ra = reffak;

  //  cout << "Partition, l = " << l << ", h = " << h << endl;
  CalcPartition (l, h, r1, r2, ra, elto0, curvepoints);
  //  cout << "curvepoints = " << curvepoints << endl;

  dt = 1.0 / n;

  l = 0;
  j = 1;

  pold = GetPoint (0);
  lold = 0;
  oldmark = pold;
  edgelengthold = 0;
  ARRAY<int> locsearch;

  for (i = 1; i <= n; i++)
    {
      p = GetPoint (i*dt);
      l = lold + Dist (p, pold);
      while (j < curvepoints.Size() && (l >= curvepoints[j] || i == n))
	{
	  frac = (curvepoints[j]-lold) / (l-lold);
	  mark = pold + frac * (p-pold);
	  edgelength = i*dt + (frac-1)*dt;
	  {
	    PointIndex pi1 = -1, pi2 = -1;
	  
	    Point3d mark3(mark(0), mark(1), 0);
	    Point3d oldmark3(oldmark(0), oldmark(1), 0);

	    Vec<3> v (1e-4*h, 1e-4*h, 1e-4*h);
	    searchtree.GetIntersecting (oldmark3 - v, oldmark3 + v, locsearch);
	    if (locsearch.Size()) pi1 = locsearch[0];
	      
	    searchtree.GetIntersecting (mark3 - v, mark3 + v, locsearch);
	    if (locsearch.Size()) pi2 = locsearch[0];
	    /*	    
	      for (PointIndex pk = PointIndex::BASE; 
	      pk < mesh.GetNP()+PointIndex::BASE; pk++)
	      {
	      if (Dist (mesh[pk], oldmark3) < 1e-4 * h) pi1 = pk;
	      if (Dist (mesh[pk], mark3) < 1e-4 * h) pi2 = pk;
	      }
	    */
	    

	    //	    cout << "pi1 = " << pi1 << endl;
	    //	    cout << "pi2 = " << pi2 << endl;
	    
	    if (pi1 == -1)
	      {
		pi1 = mesh.AddPoint(oldmark3);
		searchtree.Insert (oldmark3, pi1);
	      }
	    if (pi2 == -1)
	      {
		pi2 = mesh.AddPoint(mark3);
		searchtree.Insert (mark3, pi2);
	      }

	    // cout << "pi1 = " << pi1 << endl;
	    // cout << "pi2 = " << pi2 << endl;
	  
	    Segment seg;
	    seg.edgenr = segnr;
	    seg.si = bc; // segnr;
	    seg.p1 = pi1;
	    seg.p2 = pi2;
	    seg.domin = leftdom;
	    seg.domout = rightdom;
	    seg.epgeominfo[0].edgenr = segnr;
	    seg.epgeominfo[0].dist = edgelengthold;
	    seg.epgeominfo[1].edgenr = segnr;
	    seg.epgeominfo[1].dist = edgelength;
	    seg.singedge_left = hpref_left;
	    seg.singedge_right = hpref_right;
	    mesh.AddSegment (seg);
	  }
	
	  oldmark = mark;
	  edgelengthold = edgelength;
	  j++;
	}
    
      pold = p;
      lold = l;
    }
}


template<int D>
void SplineSeg<D> :: GetPoints (int n, ARRAY<Point<D> > & points)
{
  points.SetSize (n);
  if (n >= 2)
    for (int i = 0; i < n; i++)
      points[i] = GetPoint(double(i) / (n-1));
}

template<int D>
void SplineSeg<D> :: PrintCoeff (ostream & ost) const
{
  Vector u(6);

  GetCoeff(u);

  for ( int i=0; i<6; i++)
    ost << u[i] << "  ";
  ost << endl;
}



/* 
   Implementation of line-segment from p1 to p2
*/


template<int D>
LineSeg<D> :: LineSeg (const GeomPoint<D> & ap1, 
		       const GeomPoint<D> & ap2)
  : p1(ap1), p2(ap2)
{
  ;
}


template<int D>
Point<D> LineSeg<D> :: GetPoint (double t) const
{
  return p1 + t * (p2 - p1);
}

template<int D>
Vec<D> LineSeg<D> :: GetTangent (const double t) const
{
  return p2-p1;
}

template<int D>
void LineSeg<D> :: GetDerivatives (const double t, 
				   Point<D> & point,
				   Vec<D> & first,
				   Vec<D> & second) const
{
  first = p2 - p1;
  point = p1 + t * first;
  second = 0;
}


template<int D>
double LineSeg<D> :: Length () const
{
  return Dist (p1, p2);
}


template<int D>
void LineSeg<D> :: GetCoeff (Vector & coeffs) const
{
  coeffs.SetSize(6);

  double dx = p2(0) - p1(0);
  double dy = p2(1) - p1(1);

  coeffs[0] = coeffs[1] = coeffs[2] = 0;
  coeffs[3] = -dy;
  coeffs[4] = dx;
  coeffs[5] = -dx * p1(1) + dy * p1(0);
}



template<int D>
void LineSeg<D> :: LineIntersections (const double a, const double b, const double c,
				      ARRAY < Point<D> > & points, const double eps) const
{
  points.SetSize(0);

  double denom = -a*p2(0)+a*p1(0)-b*p2(1)+b*p1(1);
  if(fabs(denom) < 1e-20)
    return;

  double t = (a*p1(0)+b*p1(1)+c)/denom;
  if((t > -eps) && (t <  1.+eps))
    points.Append(GetPoint(t));
}



template<int D>
void LineSeg<D> :: Project (const Point<D> point, Point<D> & point_on_curve, double & t) const
{
  Vec<D> v = p2-p1;
  double l = v.Length();
  v *= 1./l;
  t = (point-p1)*v;

  if(t<0) t = 0;
  if(t>l) t = l;

  point_on_curve = p1+t*v;

  t *= 1./l;
}


template<int D>
void LineSeg<D> :: GetRawData (ARRAY<double> & data) const
{
  data.Append(2);
  for(int i=0; i<D; i++)
    data.Append(p1[i]);
  for(int i=0; i<D; i++)
    data.Append(p2[i]);
}


template<int D>
void SplineSeg3<D> :: Project (const Point<D> point, Point<D> & point_on_curve, double & t) const
{
  double t_old = -1;

  if(proj_latest_t > 0. && proj_latest_t < 1.)
    t = proj_latest_t;
  else
    t = 0.5;
	
  Point<D> phi;
  Vec<D> phip,phipp,phimp;
    
  int i=0;

  while(t > -0.5 && t < 1.5 && i<20 && fabs(t-t_old) > 1e-15 )
    {
      GetDerivatives(t,phi,phip,phipp);
	
      t_old = t;

      phimp = phi-point;

      //t = min2(max2(t-(phip*phimp)/(phipp*phimp + phip*phip),0.),1.);
      t -= (phip*phimp)/(phipp*phimp + phip*phip);

      i++;
    }
    
  //if(i<10 && t > 0. && t < 1.)
  if(i<20 && t > -0.4 && t < 1.4)
    {
      if(t < 0)
	{
	  t = 0.;
	}
      if(t > 1)
	{
	  t = 1.;
	}

      point_on_curve = GetPoint(t);
	
      double dist = Dist(point,point_on_curve);
	
      phi = GetPoint(0);
      double auxdist = Dist(phi,point);
      if(auxdist < dist)
	{
	  t = 0.;
	  point_on_curve = phi;
	  dist = auxdist;
	}
      phi = GetPoint(1);
      auxdist = Dist(phi,point);
      if(auxdist < dist)
	{
	  t = 1.;
	  point_on_curve = phi;
	  dist = auxdist;
	}
    }
  else
    {
      double t0 = 0;
      double t1 = 0.5;
      double t2 = 1.;

      double d0,d1,d2;

	
      //(*testout) << "newtonersatz" << endl;
      while(t2-t0 > 1e-8)
	{
	    
	  phi = GetPoint(t0); d0 = Dist(phi,point);
	  phi = GetPoint(t1); d1 = Dist(phi,point);
	  phi = GetPoint(t2); d2 = Dist(phi,point);

	  double a = (2.*d0 - 4.*d1 +2.*d2)/pow(t2-t0,2);

	  if(a <= 0)
	    {
	      if(d0 < d2)
		t2 -= 0.3*(t2-t0);
	      else
		t0 += 0.3*(t2-t0);

	      t1 = 0.5*(t2+t0);
	    }
	  else
	    {
	      double b = (d1-d0-a*(t1*t1-t0*t0))/(t1-t0);

	      double auxt1 = -0.5*b/a;

	      if(auxt1 < t0)
		{
		  t2 -= 0.4*(t2-t0);
		  t0 = max2(0.,t0-0.1*(t2-t0));
		}
	      else if (auxt1 > t2)
		{
		  t0 += 0.4*(t2-t0);
		  t2 = min2(1.,t2+0.1*(t2-t0));
		}
	      else
		{
		  t1 = auxt1;
		  auxt1 = 0.25*(t2-t0);
		  t0 = max2(0.,t1-auxt1);
		  t2 = min2(1.,t1+auxt1);
		}
		
	      t1 = 0.5*(t2+t0);
	    }  

	}

	
      phi = GetPoint(t0); d0 = Dist(phi,point);
      phi = GetPoint(t1); d1 = Dist(phi,point);
      phi = GetPoint(t2); d2 = Dist(phi,point);

      double mind = d0;
      t = t0;
      if(d1 < mind)
	{
	  t = t1;
	  mind = d1;
	}
      if(d2 < mind)
	{
	  t = t2;
	  mind = d2;
	}

      point_on_curve = GetPoint(t);
    }
  //(*testout) << " latest_t " << proj_latest_t << " t " << t << endl;

  proj_latest_t = t;
}




template<int D>
SplineSeg3<D> :: SplineSeg3 (const GeomPoint<D> & ap1, 
			     const GeomPoint<D> & ap2,
			     const GeomPoint<D> & ap3)
  : p1(ap1), p2(ap2), p3(ap3)
{
  proj_latest_t = 0.5;
}

template<int D>
Point<D> SplineSeg3<D> :: GetPoint (double t) const
{
  double x, y, w;
  double b1, b2, b3;

  b1 = (1-t)*(1-t);
  b2 = sqrt(2.0) * t * (1-t);
  b3 = t * t;

  x = p1(0) * b1 + p2(0) * b2 + p3(0) * b3;
  y = p1(1) * b1 + p2(1) * b2 + p3(1) * b3;
  w = b1 + b2 + b3;

  if(D==3)
    {
      double z = p1(2) * b1 + p2(2) * b2 + p3(2) * b3;
      return Point<D> (x/w, y/w, z/w);
    }
  else
    return Point<D> (x/w, y/w);
}



template<int D>
void SplineSeg3<D> :: GetDerivatives (const double t, 
				      Point<D> & point,
				      Vec<D> & first,
				      Vec<D> & second) const
{
  Vec<D> v1(p1), v2(p2), v3(p3);

  double b1 = (1.-t)*(1.-t);
  double b2 = sqrt(2.)*t*(1.-t);
  double b3 = t*t;
  double w = b1+b2+b3;
  b1 *= 1./w; b2 *= 1./w; b3 *= 1./w;

  double b1p = 2.*(t-1.);
  double b2p = sqrt(2.)*(1.-2.*t);
  double b3p = 2.*t;
  const double wp = b1p+b2p+b3p;
  const double fac1 = wp/w;
  b1p *= 1./w; b2p *= 1./w; b3p *= 1./w;

  const double b1pp = 2.;
  const double b2pp = -2.*sqrt(2.);
  const double b3pp = 2.;
  const double wpp = b1pp+b2pp+b3pp;
  const double fac2 = (wpp*w-2.*wp*wp)/(w*w);

  for(int i=0; i<D; i++)
    point(i) = b1*p1(i) + b2*p2(i) + b3*p3(i);
    
 
  first = (b1p - b1*fac1) * v1 +
    (b2p - b2*fac1) * v2 +
    (b3p - b3*fac1) * v3;

  second = (b1pp/w - b1p*fac1 - b1*fac2) * v1 +
    (b2pp/w - b2p*fac1 - b2*fac2) * v2 +
    (b3pp/w - b3p*fac1 - b3*fac2) * v3;
}



template<int D>
Vec<D> SplineSeg3<D> :: GetTangent (const double t) const
{
  const double b1 = (1.-t)*((sqrt(2.)-2.)*t-sqrt(2.));
  const double b2 = sqrt(2.)*(1.-2.*t);
  const double b3 = t*((sqrt(2.)-2)*t+2.);


  Vec<D> retval;
  for(int i=0; i<D; i++) 
    retval(i) = b1*p1(i) + b2*p2(i) + b3*p3(i);

  return retval;

}


template<int D>
void SplineSeg3<D> :: GetCoeff (Vector & u) const
{
  double t;
  int i;
  Point<D> p;
  DenseMatrix a(6, 6);
  DenseMatrix ata(6, 6);
  Vector f(6);

  u.SetSize(6);

  //  ata.SetSymmetric(1);

  t = 0;
  for (i = 1; i <= 5; i++, t += 0.25)
    {
      p = GetPoint (t);
      a.Elem(i, 1) = p(0) * p(0);
      a.Elem(i, 2) = p(1) * p(1);
      a.Elem(i, 3) = p(0) * p(1);
      a.Elem(i, 4) = p(0);
      a.Elem(i, 5) = p(1);
      a.Elem(i, 6) = 1;
    }
  a.Elem(6, 1) = 1;

  CalcAtA (a, ata);

  u = 0;
  u.Elem(6) = 1;
  a.MultTrans (u, f);
  ata.Solve (f, u);
}

/*
template<int D>
double SplineSeg3<D> :: MaxCurvature(void) const
{
  Vec<D> v1 = p1-p2;
  Vec<D> v2 = p3-p2;
  double l1 = v1.Length();
  double l2 = v2.Length();
  (*testout) << "v1 " << v1 << " v2 " << v2 << endl;

  double cosalpha = v1*v2/(l1*l2);

  (*testout) << "cosalpha " << cosalpha << endl;

  return sqrt(cosalpha + 1.)/(min2(l1,l2)*(1.-cosalpha));
}
*/
  

template<int D>
void SplineSeg3<D> :: LineIntersections (const double a, const double b, const double c,
					 ARRAY < Point<D> > & points, const double eps) const
{
  points.SetSize(0);

  double t;

  const double c1 = a*p1(0) - sqrt(2.)*a*p2(0) + a*p3(0) 
    + b*p1(1) - sqrt(2.)*b*p2(1) + b*p3(1) 
    + (2.-sqrt(2.))*c;
  const double c2 = -2.*a*p1(0) + sqrt(2.)*a*p2(0) -2.*b*p1(1) + sqrt(2.)*b*p2(1) + (sqrt(2.)-2.)*c;
  const double c3 = a*p1(0) + b*p1(1) + c;

  if(fabs(c1) < 1e-20)
    {
      if(fabs(c2) < 1e-20)
	return;

      t = -c3/c2;
      if((t > -eps) && (t < 1.+eps))
	points.Append(GetPoint(t));
      return;
    }

  const double discr = c2*c2-4.*c1*c3;

  if(discr < 0)
    return;

  if(fabs(discr/(c1*c1)) < 1e-14)
    {
      t = -0.5*c2/c1;
      if((t > -eps) && (t < 1.+eps))
	points.Append(GetPoint(t));
      return;
    }

  t = (-c2 + sqrt(discr))/(2.*c1);
  if((t > -eps) && (t < 1.+eps))
    points.Append(GetPoint(t));

  t = (-c2 - sqrt(discr))/(2.*c1);
  if((t > -eps) && (t < 1.+eps))
    points.Append(GetPoint(t));
}


template < int D >
void SplineSeg3<D> :: GetRawData (ARRAY<double> & data) const
{
  data.Append(3);
  for(int i=0; i<D; i++)
    data.Append(p1[i]);
  for(int i=0; i<D; i++)
    data.Append(p2[i]);
  for(int i=0; i<D; i++)
    data.Append(p3[i]);
}


//########################################################################
//		circlesegment

template<int D>
CircleSeg<D> :: CircleSeg (const GeomPoint<D> & ap1, 
			   const GeomPoint<D> & ap2,
			   const GeomPoint<D> & ap3)
  : p1(ap1), p2(ap2), p3(ap3)
{
  Vec<D> v1,v2;
  
  v1 = p1 - p2;
  v2 = p3 - p2;
  
  Point<D> p1t(p1+v1);
  Point<D> p2t(p3+v2);

  // works only in 2D!!!!!!!!!
    
  Line2d g1t,g2t;

  g1t.P1() = Point<2>(p1(0),p1(1));
  g1t.P2() = Point<2>(p1t(0),p1t(1));
  g2t.P1() = Point<2>(p3(0),p3(1));
  g2t.P2() = Point<2>(p2t(0),p2t(1));

  Point<2> mp = CrossPoint (g1t,g2t);

  pm(0) = mp(0); pm(1) = mp(1);
  radius  = Dist(pm,StartPI());
  Vec2d auxv;
  auxv.X() = p1(0)-pm(0); auxv.Y() = p1(1)-pm(1);
  w1      = Angle(auxv);
  auxv.X() = p3(0)-pm(0); auxv.Y() = p3(1)-pm(1);
  w3      = Angle(auxv);
  if ( fabs(w3-w1) > M_PI )
    {  
      if ( w3>M_PI )   w3 -= 2*M_PI;
      if ( w1>M_PI )   w1 -= 2*M_PI;
    }
}
 

template<int D>
Point<D> CircleSeg<D> :: GetPoint (double t) const
{
  if (t >= 1.0)  { return p3; }
     
  double phi = StartAngle() + t*(EndAngle()-StartAngle());
  Vec<D>  tmp(cos(phi),sin(phi));
     
  return pm + Radius()*tmp;
}
  
template<int D>
void CircleSeg<D> :: GetCoeff (Vector & coeff) const
{ 
  coeff[0] = coeff[1] = 1.0;
  coeff[2] = 0.0;
  coeff[3] = -2.0 * pm[0];
  coeff[4] = -2.0 * pm[1];
  coeff[5] = sqr(pm[0]) + sqr(pm[1]) - sqr(Radius());
}

  
template<int D>
void CircleSeg<D> :: LineIntersections (const double a, const double b, const double c,
					ARRAY < Point<D> > & points, const double eps) const
{
  points.SetSize(0);

  double px=0,py=0;

  if(fabs(b) > 1e-20)
    py = -c/b;
  else
    px = -c/a;

  const double c1 = a*a + b*b;
  const double c2 = 2. * ( a*(py-pm(1)) - b*(px-pm(0)));
  const double c3 = pow(px-pm(0),2) + pow(py-pm(1),2) - pow(Radius(),2);
    
  const double discr = c2*c2 - 4*c1*c3;

  if(discr < 0)
    return;

  ARRAY<double> t;

  if(fabs(discr) < 1e-20)
    t.Append(-0.5*c2/c1);
  else
    {
      t.Append((-c2+sqrt(discr))/(2.*c1));
      t.Append((-c2-sqrt(discr))/(2.*c1));
    }

  for(int i=0; i<t.Size(); i++)
    {
      Point<D> p (px-t[i]*b,py+t[i]*a);

      double angle = atan2(p(1),p(0))+M_PI;

      if(angle > StartAngle()-eps && angle < EndAngle()+eps)
	points.Append(p);
    }
}




template<int D>
DiscretePointsSeg<D> ::   DiscretePointsSeg (const ARRAY<Point<D> > & apts)
  : pts (apts)
{ 
  for(int i=0; i<D; i++)
    {
      p1(i) = apts[0](i);
      p2(i) = apts.Last()(i);
    }
  p1.refatpoint = true;
  p2.refatpoint = true;
}


template<int D>
DiscretePointsSeg<D> :: ~DiscretePointsSeg ()
{ ; }

template<int D>
Point<D> DiscretePointsSeg<D> :: GetPoint (double t) const
{
  double t1 = t * (pts.Size()-1);
  int segnr = int(t1);
  if (segnr < 0) segnr = 0;
  if (segnr >= pts.Size()) segnr = pts.Size()-1;

  double rest = t1 - segnr;
    
  return pts[segnr] + rest*Vec<D>(pts[segnr+1]-pts[segnr]);
}

  

typedef GeomPoint<2> GeomPoint2d;
typedef SplineSeg<2> SplineSegment;
typedef LineSeg<2> LineSegment;
typedef SplineSeg3<2> SplineSegment3;
typedef CircleSeg<2> CircleSegment;
typedef DiscretePointsSeg<2> DiscretePointsSegment;




#endif
