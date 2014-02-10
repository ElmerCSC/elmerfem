/*

Spline curve for Mesh generator

*/

#include <mystdlib.h>
#include <csg.hpp>
#include <linalg.hpp>
#include <meshing.hpp>

namespace netgen
{
#include "spline.hpp"

  /*
  template<> void SplineSeg3<2> :: Project (const Point<2> point, Point<2> & point_on_curve, double & t) const
  {
    double t_old = 0;
    t = 0.5;
    
    Point<2> phi;
    Vec<2> phip,phipp,phimp;
    
    int i=0;

    while(fabs(t-t_old) > 1e-8 && i<10)
      {
	GetDerivatives(t,phi,phip,phipp);
	
	t_old = t;

	phimp = phi-point;

	t = min2(max2(t-(phip*phimp)/(phipp*phimp + phip*phip),0.),1.);

	i++;
      }

    if(i<10)
      {
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

	
	//(*testout) << "2d newtonersatz" << endl;
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
  }

  template<> void SplineSeg3<3> :: Project (const Point<3> point, Point<3> & point_on_curve, double & t) const
  {
    double t_old = -1;

    if(proj_latest_t > 0. && proj_latest_t < 1.)
      t = proj_latest_t;
    else
      t = 0.5;
	
    Point<3> phi;
    Vec<3> phip,phipp,phimp;
    
    int i=0;

    while(fabs(t-t_old) > 1e-8 && t > -0.5 && t < 1.5 && i<10)
      {
	GetDerivatives(t,phi,phip,phipp);
	
	t_old = t;

	phimp = phi-point;

	//t = min2(max2(t-(phip*phimp)/(phipp*phimp + phip*phip),0.),1.);
	t -= (phip*phimp)/(phipp*phimp + phip*phip);

	i++;
      }
    
    //if(i<10 && t > 0. && t < 1.)
    if(i<10 && t > -0.4 && t < 1.4)
      {
	if(t < 0)
	  t = 0.;
	if(t > 1)
	  t = 1.;

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
  */

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

  template<>
  double SplineSeg3<2> :: MaxCurvature(void) const
  {
    Vec<2> v1 = p1-p2;
    Vec<2> v2 = p3-p2;
    double l1 = v1.Length();
    double l2 = v2.Length();
        
    double cosalpha = (v1*v2)/(l1*l2);
    
            
    return sqrt(cosalpha + 1.)/(min2(l1,l2)*(1.-cosalpha));
  }

  template<>
  double SplineSeg3<3> :: MaxCurvature(void) const
  {
    Vec<3> v1 = p1-p2;
    Vec<3> v2 = p3-p2;
    double l1 = v1.Length();
    double l2 = v2.Length();
        
    double cosalpha = v1*v2/(l1*l2);
    
        
    return sqrt(cosalpha + 1.)/(min2(l1,l2)*(1.-cosalpha));
  }




}
