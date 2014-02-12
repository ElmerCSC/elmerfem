
2d Spline curve for Mesh generator


OLD IMPLEMENTATION, NOT USED ANYMORE!



#include <mystdlib.h>
#include <csg.hpp>
#include <linalg.hpp>
#include <meshing.hpp>

namespace netgen
{
#include "spline2d.hpp"


  void CalcPartition (double l, double h, double r1, double r2,
		      double ra, double elto0, ARRAY<double> & points);



  // calculates length of spline-curve
  double SplineSegment :: Length () const
  {
    Point<2> p, pold;

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
  void SplineSegment :: Partition (double h, double elto0,
				   Mesh & mesh, Point3dTree & searchtree, int segnr) const
  {
    int i, j;
    double l, r1, r2, ra;
    double lold, dt, frac;
    int n = 100;
    Point<2> p, pold, mark, oldmark;
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


  void SplineSegment :: GetPoints (int n, ARRAY<Point<2> > & points)
  {
    points.SetSize (n);
    if (n >= 2)
      for (int i = 0; i < n; i++)
	points[i] = GetPoint(double(i) / (n-1));
  }

  void SplineSegment :: PrintCoeff (ostream & ost) const
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


  LineSegment :: LineSegment (const GeomPoint2d & ap1, 
			      const GeomPoint2d & ap2)
    : p1(ap1), p2(ap2)
  {
    ;
  }


  Point<2> LineSegment :: GetPoint (double t) const
  {
    return p1 + t * (p2 - p1);
  }

  double LineSegment :: Length () const
  {
    return Dist (p1, p2);
  }


  void LineSegment :: GetCoeff (Vector & coeffs) const
  {
    coeffs.SetSize(6);

    double dx = p2(0) - p1(0);
    double dy = p2(1) - p1(1);

    coeffs[0] = coeffs[1] = coeffs[2] = 0;
    coeffs[3] = dy;
    coeffs[4] = -dx;
    coeffs[5] = dx * p1(1) - dy * p1(0);
  }



  void LineSegment :: LineIntersections (const double a, const double b, const double c,
					 ARRAY < Point<2> > & points, const double eps) const
  {
    points.SetSize(0);

    double denom = -a*p2(0)+a*p1(0)-b*p2(1)+b*p1(1);
    if(fabs(denom) < 1e-20)
      return;

    double t = (a*p1(0)+b*p1(1)+c)/denom;
    if((t > -eps) && (t < 1.+eps))
      points.Append(GetPoint(t));
  }


  SplineSegment3 :: SplineSegment3 (const GeomPoint2d & ap1, 
				    const GeomPoint2d & ap2,
				    const GeomPoint2d & ap3)
    : p1(ap1), p2(ap2), p3(ap3)
  {
    ;
  }

  Point<2> SplineSegment3 :: GetPoint (double t) const
  {
    double x, y, w;
    double b1, b2, b3;

    b1 = (1-t)*(1-t);
    b2 = sqrt(2.0) * t * (1-t);
    b3 = t * t;

    x = p1(0) * b1 + p2(0) * b2 + p3(0) * b3;
    y = p1(1) * b1 + p2(1) * b2 + p3(1) * b3;
    w = b1 + b2 + b3;

    return Point<2> (x/w, y/w);
  }


  void SplineSegment3 :: GetCoeff (Vector & u) const
  {
    double t;
    int i;
    Point<2> p;
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

  double SplineSegment3 :: MaxCurvature(void) const
  {
    Vec<2> v1 = p1-p2;
    Vec<2> v2 = p3-p2;
    double l1 = v1.Length();
    double l2 = v2.Length();

    double cosalpha = v1*v2/(l1*l2);


    return sqrt(cosalpha + 1.)/(min2(l1,l2)*(1.-cosalpha));
  }
  

  void SplineSegment3 :: LineIntersections (const double a, const double b, const double c,
					    ARRAY < Point<2> > & points, const double eps) const
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

    const double D = c2*c2-4.*c1*c3;

    if(D < 0)
      return;

    if(fabs(D/(c1*c1)) < 1e-14)
      {
	t = -0.5*c2/c1;
	if((t > -eps) && (t < 1.+eps))
	  points.Append(GetPoint(t));
	return;
      }

    t = (-c2 + sqrt(D))/(2.*c1);
    if((t > -eps) && (t < 1.+eps))
      points.Append(GetPoint(t));

    t = (-c2 - sqrt(D))/(2.*c1);
    if((t > -eps) && (t < 1.+eps))
      points.Append(GetPoint(t));
  }



  //########################################################################
  //		circlesegment

  CircleSegment :: CircleSegment (const GeomPoint2d & ap1, 
				  const GeomPoint2d & ap2,
				  const GeomPoint2d & ap3)
    : p1(ap1), p2(ap2), p3(ap3)
  {
    Vec<2> v1,v2;
  
    v1 = p1 - p2;
    v2 = p3 - p2;
  
    Point<2> p1t(p1(0)+v1[1], p1(1)-v1[0]);
    Point<2> p2t(p3(0)+v2[1], p3(1)-v2[0]);
    Line2d g1t(p1, p1t), g2t(p3, p2t);
  
    pm 	  = CrossPoint (g1t,g2t);
    radius  = Dist(pm,StartPI());
    w1      = Angle(Vec2d (p1 - pm));
    w3      = Angle(Vec2d (p3 - pm));
    if ( fabs(w3-w1) > M_PI )
      {  
	if ( w3>M_PI )   w3 -= 2*M_PI;
	if ( w1>M_PI )   w1 -= 2*M_PI;
      }
  }
 
  Point<2>  CircleSegment :: GetPoint (double t) const
  {
    if (t >= 1.0)  { return p3; }
     
    double phi = StartAngle() + t*(EndAngle()-StartAngle());
    Vec2d  tmp(cos(phi),sin(phi));
     
    return pm + Radius()*tmp;
  }
  
  void CircleSegment :: GetCoeff (Vector & coeff) const
  { 
    coeff[0] = coeff[1] = 1.0;
    coeff[2] = 0.0;
    coeff[3] = -2.0 * pm[0];
    coeff[4] = -2.0 * pm[1];
    coeff[5] = sqr(pm[0]) + sqr(pm[1]) - sqr(Radius());
  }

  
  void CircleSegment :: LineIntersections (const double a, const double b, const double c,
					   ARRAY < Point<2> > & points, const double eps) const
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
    
    const double D = c2*c2 - 4*c1*c3;

    if(D < 0)
      return;

    ARRAY<double> t;

    if(fabs(D) < 1e-20)
      t.Append(-0.5*c2/c1);
    else
      {
	t.Append((-c2+sqrt(D))/(2.*c1));
	t.Append((-c2-sqrt(D))/(2.*c1));
      }

    for(int i=0; i<t.Size(); i++)
      {
	Point<2> p (px-t[i]*b,py+t[i]*a);

	double angle = atan2(p(1),p(0))+M_PI;

	if(angle > StartAngle()-eps && angle < EndAngle()+eps)
	  points.Append(p);
      }
  }



  DiscretePointsSegment ::   DiscretePointsSegment (const ARRAY<Point<2> > & apts)
    : pts (apts),
      p1 (apts[0](0), apts[0](1), 1),
      p2 (apts.Last()(0), apts.Last()(1), 1)

  { ; }


  DiscretePointsSegment :: ~DiscretePointsSegment ()
  { ; }

  Point<2> DiscretePointsSegment :: GetPoint (double t) const
  {
    double t1 = t * (pts.Size()-1);
    int segnr = int(t1);
    if (segnr < 0) segnr = 0;
    if (segnr >= pts.Size()) segnr = pts.Size()-1;

    double rest = t1 - segnr;
    
    return Point<2> ((1-rest)*pts[segnr](0) + rest*pts[segnr+1](0),
		     (1-rest)*pts[segnr](1) + rest*pts[segnr+1](1));
  }

  




  //########################################################################




  void CalcPartition (double l, double h, double r1, double r2,
		      double ra, double elto0, ARRAY<double> & points)
  {
    int i, j, n, nel;
    double sum, t, dt, fun, fperel, oldf, f;

    n = 1000;

    points.SetSize (0);

    sum = 0;
    dt = l / n;
    t = 0.5 * dt;
    for (i = 1; i <= n; i++)
      {
	fun = min3 (h/ra, t/elto0 + h/r1, (l-t)/elto0 + h/r2);
	sum += dt / fun;
	t += dt;
      }

    nel = int (sum+1);
    fperel = sum / nel;

    points.Append (0);

    i = 1;
    oldf = 0;
    t = 0.5 * dt;
    for (j = 1; j <= n && i < nel; j++)
      {
	fun = min3 (h/ra, t/elto0 + h/r1, (l-t)/elto0 + h/r2);

	f = oldf + dt / fun;

	while (f > i * fperel && i < nel)
	  {
	    points.Append ( (l/n) * (j-1 +  (i * fperel - oldf) / (f - oldf)) );
	    i++;
	  }
	oldf = f;
	t += dt;
      }
    points.Append (l);
  }


}
