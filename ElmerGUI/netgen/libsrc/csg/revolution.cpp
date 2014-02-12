#include <mystdlib.h>

#include <linalg.hpp>
#include <csg.hpp>

namespace netgen
{
  void RevolutionFace :: Init(void)
  {
    const LineSeg<2> * line = dynamic_cast<const LineSeg<2>*>(spline);
    const SplineSeg3<2> * spline3 = dynamic_cast<const SplineSeg3<2>*>(spline);

    if(line)
      {
	checklines_start.Append(new Point<2>(line->StartPI()));
	checklines_vec.Append(new Vec<2>(line->EndPI() - line->StartPI()));
	(*checklines_vec.Last()) *= 1./pow(checklines_vec.Last()->Length(),2); //!!
      }
    else if (spline3)
      {
	checklines_start.Append(new Point<2>(spline3->EndPI()));
	checklines_start.Append(new Point<2>(spline3->TangentPoint()));
	checklines_start.Append(new Point<2>(spline3->StartPI()));
	checklines_vec.Append(new Vec<2>(spline3->StartPI() - spline3->EndPI()));
	(*checklines_vec.Last()) *= 1./pow(checklines_vec.Last()->Length(),2); //!!
	checklines_vec.Append(new Vec<2>(spline3->EndPI() - spline3->TangentPoint()));
	(*checklines_vec.Last()) *= 1./pow(checklines_vec.Last()->Length(),2); //!!
	checklines_vec.Append(new Vec<2>(spline3->TangentPoint() - spline3->StartPI()));
	(*checklines_vec.Last()) *= 1./pow(checklines_vec.Last()->Length(),2); //!!
	
      }
    
    for(int i=0; i<checklines_vec.Size(); i++)
      {
	checklines_normal.Append(new Vec<2>);
	(*checklines_normal.Last())(0) = - (*checklines_vec[i])(1);
	(*checklines_normal.Last())(1) = (*checklines_vec[i])(0);
	checklines_normal.Last()->Normalize();
      }
  }

  RevolutionFace :: RevolutionFace(const SplineSeg<2> & spline_in,
				   const Point<3> & p,
				   const Vec<3> & vec,
				   bool first,
				   bool last,
				   const int id_in) :
    spline(&spline_in), p0(p), v_axis(vec), isfirst(first), islast(last), id(id_in)
  {    
    deletable = false;
    Init();
  }


  RevolutionFace :: RevolutionFace(const ARRAY<double> & raw_data)
  {
    deletable = true;
    
    int pos = 0;

    ARRAY< Point<2> > p(3);

    int stype = int(raw_data[pos]); pos++;

    for(int i=0; i<stype; i++)
      {
	p[i](0) = raw_data[pos]; pos++;
	p[i](1) = raw_data[pos]; pos++;
      }

    if(stype == 2)
      {
	spline = new LineSeg<2>(GeomPoint<2>(p[0],1),
				GeomPoint<2>(p[1],1));
	//(*testout) << "appending LineSeg<2> " << p[0] 
	//	   << " to " << p[1] << endl;
      }
    else if(stype == 3)
      {
	spline = new SplineSeg3<2>(GeomPoint<2>(p[0],1),
				   GeomPoint<2>(p[1],1),
				   GeomPoint<2>(p[2],1));
	//(*testout) << "appending SplineSeg<3> "
	//	   << p[0] << " -- " << p[1] << " -- " << p[2] << endl;
      }

    for(int i=0; i<3; i++)
      {
	p0(i) = raw_data[pos];
	pos++;
      }
    for(int i=0; i<3; i++)
      {
	v_axis(i) = raw_data[pos];
	pos++;
      }
    isfirst = (raw_data[pos] > 0.9);
    pos++;
    islast = (raw_data[pos] < 0.1);
    pos++;
    

  }

  
  RevolutionFace :: ~RevolutionFace()
  {
    for(int i=0; i<checklines_start.Size(); i++)
      {
	delete checklines_start[i];
	delete checklines_vec[i];
	delete checklines_normal[i];
      }

    if(deletable)
      delete spline;
  }
  
  void RevolutionFace :: CalcProj(const Point<3> & point3d, Point<2> & point2d,
				  const Vec<3> & vector3d, Vec<2> & vector2d) const
  {
    Vec<3> pmp0 = point3d-p0;
    CalcProj0(pmp0,point2d);
    Vec<3> y=pmp0-point2d(0)*v_axis; y.Normalize();
    vector2d(0) = vector3d*v_axis;
    vector2d(1) = vector3d*y;
  }
  

  void RevolutionFace :: CalcProj(const Point<3> & point3d, Point<2> & point2d) const
  {
    Vec<3> pmp0 = point3d-p0;
    CalcProj0(pmp0,point2d);
  }
  
  void RevolutionFace :: CalcProj0(const Vec<3> & point3d_minus_p0, Point<2> & point2d) const
  {
    point2d(0) = point3d_minus_p0 * v_axis;
    point2d(1) = sqrt( point3d_minus_p0 * point3d_minus_p0 - point2d(0)*point2d(0) );
  }
  

  int  RevolutionFace ::IsIdentic (const Surface & s2, int & inv, double eps) const
  {
    const RevolutionFace * rev2 = dynamic_cast<const RevolutionFace*>(&s2);
    
    if(!rev2) return 0;
    
    if(rev2 == this)
      return 1;
        
    return 0;
  }

  double RevolutionFace :: CalcFunctionValue (const Point<3> & point) const
  {
    if(spline_coefficient.Size() == 0)
      spline->GetCoeff(spline_coefficient);

    Point<2> p;
    CalcProj(point,p);

    return spline_coefficient(0)*p(0)*p(0) + spline_coefficient(1)*p(1)*p(1)
      + spline_coefficient(2)*p(0)*p(1) + spline_coefficient(3)*p(0)
      + spline_coefficient(4)*p(1) + spline_coefficient(5);
  }

  void RevolutionFace :: CalcGradient (const Point<3> & point, Vec<3> & grad) const
  {
    if(spline_coefficient.Size() == 0)
      spline->GetCoeff(spline_coefficient);

    Vec<3> point_minus_p0 = point-p0;

    Point<2> p;
    CalcProj0(point_minus_p0,p);

    const double dFdxbar = 2.*spline_coefficient(0)*p(0) + spline_coefficient(2)*p(1) + spline_coefficient(3);

    if(fabs(p(1)) > 1e-10)
      {
	const double dFdybar = 2.*spline_coefficient(1)*p(1) + spline_coefficient(2)*p(0) + spline_coefficient(4);

	grad(0) = dFdxbar*v_axis(0) + dFdybar * ( point_minus_p0(0)-v_axis(0)*p(0) )/p(1);
	grad(1) = dFdxbar*v_axis(1) + dFdybar * ( point_minus_p0(1)-v_axis(1)*p(0) )/p(1);
	grad(2) = dFdxbar*v_axis(2) + dFdybar * ( point_minus_p0(2)-v_axis(2)*p(0) )/p(1);
	//(*testout) << "grad1("<<point<<") = " << grad << endl;
      }
    else
      {
	grad(0) = dFdxbar*v_axis(0);
	grad(1) = dFdxbar*v_axis(1);
	grad(2) = dFdxbar*v_axis(2);
	//(*testout) << "grad2("<<point<<") = " << grad << endl;
      }
  }

  
  void RevolutionFace :: CalcHesse (const Point<3> & point, Mat<3> & hesse) const
  {
    if(spline_coefficient.Size() == 0)
      spline->GetCoeff(spline_coefficient);

    Vec<3> point_minus_p0 = point-p0;

    Point<2> p;
    CalcProj0(point_minus_p0,p);


    if(fabs(p(1)) > 1e-10)
      {
	const double dFdybar = 2.*spline_coefficient(1)*p(1) + spline_coefficient(2)*p(0) + spline_coefficient(4);
	
	const double aux = -pow(p(1),-3);
	const double aux0 = point_minus_p0(0) - v_axis(0)*p(0);
	const double aux1 = point_minus_p0(1) - v_axis(1)*p(0);
	const double aux2 = point_minus_p0(2) - v_axis(2)*p(0);
	

	const double dybardx = aux0/p(1);
	const double dybardy = aux1/p(1);
	const double dybardz = aux2/p(1);
    
	const double dybardxx = aux*aux0*aux0 + (1.-v_axis(0)*v_axis(0))/p(1);
	const double dybardyy = aux*aux1*aux1 + (1.-v_axis(1)*v_axis(1))/p(1);
	const double dybardzz = aux*aux2*aux2 + (1.-v_axis(2)*v_axis(2))/p(1);
	const double dybardxy = aux*aux0*aux1 - v_axis(0)*v_axis(1)/p(1);
	const double dybardxz = aux*aux0*aux2 - v_axis(0)*v_axis(2)/p(1);
	const double dybardyz = aux*aux1*aux2 - v_axis(1)*v_axis(2)/p(1);
	
	hesse(0,0) = 2.*spline_coefficient(0)*v_axis(0)*v_axis(0) + 2.*spline_coefficient(2)*v_axis(0)*dybardx + 2.*spline_coefficient(1)*dybardx*dybardx
	  + dFdybar*dybardxx;
	hesse(1,1) = 2.*spline_coefficient(0)*v_axis(1)*v_axis(1) + 2.*spline_coefficient(2)*v_axis(1)*dybardy + 2.*spline_coefficient(1)*dybardy*dybardy
	  + dFdybar*dybardyy;
	hesse(2,2) = 2.*spline_coefficient(0)*v_axis(2)*v_axis(2) + 2.*spline_coefficient(2)*v_axis(2)*dybardz + 2.*spline_coefficient(1)*dybardz*dybardz
	  + dFdybar*dybardzz;
	
	hesse(0,1) = hesse(1,0) = 2.*spline_coefficient(0)*v_axis(0)*v_axis(1) + spline_coefficient(2)*v_axis(0)*dybardy + spline_coefficient(2)*dybardx*v_axis(1)
	  + 2.*spline_coefficient(2)*dybardx*dybardy + dFdybar*dybardxy;
	hesse(0,2) = hesse(2,0) = 2.*spline_coefficient(0)*v_axis(0)*v_axis(2) + spline_coefficient(2)*v_axis(0)*dybardz + spline_coefficient(2)*dybardx*v_axis(2)
	  + 2.*spline_coefficient(2)*dybardx*dybardz + dFdybar*dybardxz;
	hesse(1,2) = hesse(2,1) = 2.*spline_coefficient(0)*v_axis(1)*v_axis(2) + spline_coefficient(2)*v_axis(1)*dybardz + spline_coefficient(2)*dybardy*v_axis(2)
	  + 2.*spline_coefficient(2)*dybardy*dybardz + dFdybar*dybardyz;

	//(*testout) << "hesse1: " << hesse <<endl;
      }
    else if (fabs(spline_coefficient(2)) + fabs(spline_coefficient(4)) < 1.e-9 &&
	     fabs(spline_coefficient(0)) > 1e-10)
      {
	double aux = spline_coefficient(0)-spline_coefficient(1);
	
	hesse(0,0) = aux*v_axis(0)*v_axis(0) + spline_coefficient(1);
	hesse(0,0) = aux*v_axis(1)*v_axis(1) + spline_coefficient(1);
	hesse(0,0) = aux*v_axis(2)*v_axis(2) + spline_coefficient(1);

	hesse(0,1) = hesse(1,0) = aux*v_axis(0)*v_axis(1);
	hesse(0,2) = hesse(2,0) = aux*v_axis(0)*v_axis(2);
	hesse(1,2) = hesse(2,1) = aux*v_axis(1)*v_axis(2);
	//(*testout) << "hesse2: " << hesse <<endl;
	
      }
    else if (fabs(spline_coefficient(1)) + fabs(spline_coefficient(3)) + fabs(spline_coefficient(4)) + fabs(spline_coefficient(5)) < 1.e-9) // line
      {
	hesse = 0;
	//(*testout) << "hesse3: " << hesse <<endl;
      }
    else
      {
	(*testout) << "hesse4: " << hesse <<endl;
      }
  }

  

  double RevolutionFace ::HesseNorm () const
  {
    if (fabs(spline_coefficient(1)) + fabs(spline_coefficient(3)) + fabs(spline_coefficient(4)) + fabs(spline_coefficient(5)) < 1.e-9) // line
      return 0;
      
    if (fabs(spline_coefficient(2)) + fabs(spline_coefficient(4)) < 1.e-9 &&
	fabs(spline_coefficient(0)) > 1e-10)
      return 2.*max2(fabs(spline_coefficient(0)),fabs(spline_coefficient(1)));


    double alpha = fabs(spline_coefficient(2)*(spline->StartPI()(0)-spline->EndPI()(0))) /
      max2(fabs(spline->StartPI()(1)),fabs(spline->EndPI()(1)));

    return max2(2.*fabs(spline_coefficient(0))+sqrt(2.)*fabs(spline_coefficient(2)),
		2.*fabs(spline_coefficient(1))+spline_coefficient(2)+1.5*alpha);
  }

  double  RevolutionFace :: MaxCurvature() const
  {
    double retval = spline->MaxCurvature();

    ARRAY < Point<2> > checkpoints;

    const SplineSeg3<2> * ss3 = dynamic_cast<const SplineSeg3<2> *>(spline);
    const LineSeg<2> * ls = dynamic_cast<const LineSeg<2> *>(spline);
    
    if(ss3)
      {
	checkpoints.Append(ss3->StartPI());
	checkpoints.Append(ss3->TangentPoint());
	checkpoints.Append(ss3->TangentPoint());
	checkpoints.Append(ss3->EndPI());
      }
    else if(ls)
      {
	checkpoints.Append(ls->StartPI());
	checkpoints.Append(ls->EndPI());
      }

    for(int i=0; i<checkpoints.Size(); i+=2)
      {
	Vec<2> v = checkpoints[i+1]-checkpoints[i];
	Vec<2> n(v(1),-v(0)); n.Normalize();

	//if(ss3)
	//  (*testout) << "n " << n << endl;

	if(fabs(n(1)) < 1e-15)
	  continue;

	double t1 = -checkpoints[i](1)/n(1);
	double t2 = -checkpoints[i+1](1)/n(1);
	
	double c1 = (t1 > 0) ? (1./t1) : -1;
	double c2 = (t2 > 0) ? (1./t2) : -1;
	
	//if(ss3)
	//  (*testout) << "t1 " << t1 << " t2 " << t2 << " c1 " << c1 << " c2 " << c2 << endl;

	if(c1 > retval)
	  retval = c1;
	if(c2 > retval)
	  retval = c2;
      }
	
    //if(ss3)
    //  (*testout) << "curvature " << retval << endl;

    return retval;

    /*


    // find smallest y value of spline:
    ARRAY<double> testt;

    if(!isfirst)
      testt.Append(0);
    if(!islast)
      testt.Append(1);
    
    const SplineSegment3 * s3 = dynamic_cast<const SplineSegment3 *>(&spline);

    if(s3)
      {
	double denom = (2.-sqrt(2.))*(s3->EndPI()(1) - s3->StartPI()(1));
	
	if(fabs(denom) < 1e-20)
	  testt.Append(0.5);
	else
	  {
	    double sD = sqrt(pow(s3->TangentPoint()(1) - s3->StartPI()(1),2)+
			     pow(s3->TangentPoint()(1) - s3->EndPI()(1),2));
	    testt.Append((s3->StartPI()(1)*(sqrt(2.)-1.) - sqrt(2.)*s3->TangentPoint()(1) + s3->EndPI()(1) + sD)/denom);
	    testt.Append((s3->StartPI()(1)*(sqrt(2.)-1.) - sqrt(2.)*s3->TangentPoint()(1) + s3->EndPI()(1) - sD)/denom);
	  }	
      }

    double miny = fabs(spline.GetPoint(testt[0])(1));
    for(int i=1; i<testt.Size(); i++)
      {
	double thisy = fabs(spline.GetPoint(testt[i])(1));
	if(thisy < miny)
	  miny = thisy;
      }

    return max2(splinecurvature,1./miny);
    */
  }

  void RevolutionFace :: Project (Point<3> & p) const
  {
    Point<2> p2d;

    CalcProj(p,p2d);

    const Vec<3> y = (p-p0)-p2d(0)*v_axis;
    const double yl = y.Length();

    double dummy;

    spline->Project(p2d,p2d,dummy);

    p = p0 + p2d(0)*v_axis;

    if(yl > 1e-20*Dist(spline->StartPI(),spline->EndPI()))
      p+= (p2d(1)/yl)*y;
  }



  
  Point<3>  RevolutionFace :: GetSurfacePoint () const
  {
    Vec<3> random_vec(0.760320,-0.241175,0.60311534);

    Vec<3> n = Cross(v_axis,random_vec); n.Normalize();

    Point<2> sp = spline->GetPoint(0.5);

    Point<3> retval = p0 + sp(0)*v_axis + sp(1)*n;

    return retval;
  }


  void RevolutionFace :: Print (ostream & str) const
  {
    if(spline_coefficient.Size() == 0)
      spline->GetCoeff(spline_coefficient);

    str << p0(0) << " " << p0(1) << " " << p0(2) << " "
	<< v_axis(0) << " " << v_axis(1) << " " << v_axis(2) << " ";
    for(int i=0; i<6; i++) str << spline_coefficient(i) << " ";
    str << endl;
  }


  void RevolutionFace :: GetTriangleApproximation (TriangleApproximation & tas, 
						   const Box<3> & boundingbox, 
						   double facets) const
  {
    Vec<3> random_vec(0.760320,-0.241175,0.60311534);

    Vec<3> v1 = Cross(v_axis,random_vec); v1.Normalize();
    Vec<3> v2 = Cross(v1,v_axis); v2.Normalize();

    int n = int(2.*facets) + 1;

    int i,j;
    double phi;
    Point<3> p;

    for(i=0; i<=n; i++)
      {
	Point<2> sp = spline->GetPoint(double(i)/double(n));
	for(j=0; j<=n; j++)
	  {
	    phi = 2.*M_PI*double(j)/double(n);
	    
	    p = p0 + sp(0)*v_axis + sp(1)*cos(phi)*v1 + sp(1)*sin(phi)*v2;
	    tas.AddPoint(p);   
	  }
      }
    
    for(i=0; i<n; i++)
      for(j=0; j<n; j++)
	{
	  int pi = (n+1)*i+j;

	  tas.AddTriangle( TATriangle (id, pi,pi+1,pi+n+1));
	  tas.AddTriangle( TATriangle (id, pi+1,pi+n+1,pi+n+2));
	}
  }
  

  bool RevolutionFace :: BoxIntersectsFace(const Box<3> & box) const
  {
    Point<3> center = box.Center();

    Project(center);

    return (Dist(box.Center(),center) < 0.5*box.Diam());
  }


  /*
  bool RevolutionFace :: BoxIntersectsFace (const BoxSphere<3> & box, bool & uncertain) const
  {
    Point<2> c,pmin,pmax;
    CalcProj(box.Center(),c);
    double aux = box.Diam()/sqrt(8.);
    pmin(0) = c(0)-aux; pmin(1) = c(1)-aux;
    pmax(0) = c(0)+aux; pmax(1) = c(1)+aux;
    
    BoxSphere<2> box2d(pmin,pmax);
    return BoxIntersectsFace(box2d, uncertain);
  }

  bool RevolutionFace :: BoxIntersectsFace (const BoxSphere<2> & box, bool & uncertain) const
  {
    const LineSegment * line = dynamic_cast<const LineSegment*>(&spline);
    const SplineSegment3 * spline3 = dynamic_cast<const SplineSegment3*>(&spline);

    bool always_right = true, always_left = true;

    bool retval = false;
    bool thisint;
    bool intdirect = false;
    bool inttangent = false;
    uncertain = false;

    if(line) inttangent = true;
  
    for(int i=0; i<checklines_start.Size(); i++)
      {
	Vec<2> b = box.Center()- (*checklines_start[i]);

	double d;

	double checkdist = b * (*checklines_vec[i]);
	double ncomp = b * (*checklines_normal[i]);

	if(checkdist < 0)
	  d = b.Length();
	else if (checkdist > 1)
	  {
	    if(spline3)
	      d = Dist(box.Center(),*checklines_start[(i+1)%3]);
	    else
	      d = Dist(box.Center(),(*checklines_start[i]) 
		       + pow(checklines_vec[i]->Length(),2)*(*checklines_vec[i]));
	  }
	else 
	  d = fabs(ncomp);
	  
	thisint = (box.Diam() >= 2.*d);
	retval = retval || thisint;
	if(thisint)
	  {
	    if(i==0)
	      intdirect = true;
	    else
	      inttangent = true;
	  }

	if(ncomp > 0) always_right = false;
	else if(ncomp < 0) always_left = false;
      }

    if(retval && !(intdirect && inttangent))
      uncertain = true;

    if(!retval && spline3 && (always_right || always_left))
      {
	retval = true;
	uncertain = true;
      }
    
    return retval;	
  }  
  */
  

  INSOLID_TYPE RevolutionFace :: PointInFace (const Point<3> & p, const double eps) const
  {
    Point<2> p2d;
    CalcProj(p,p2d);

    double val = spline_coefficient(0)*p2d(0)*p2d(0) + spline_coefficient(1)*p2d(1)*p2d(1) + spline_coefficient(2)*p2d(0)*p2d(1) +
      spline_coefficient(3)*p2d(0) + spline_coefficient(4)*p2d(1) + spline_coefficient(5);

    if(val > eps)
      return IS_OUTSIDE;
    if(val < -eps)
      return IS_INSIDE;
     
    return DOES_INTERSECT;
  }

  

  void RevolutionFace :: GetRawData(ARRAY<double> & data) const
  {
    data.DeleteAll();
    spline->GetRawData(data);
    for(int i=0; i<3; i++)
      data.Append(p0(i));
    for(int i=0; i<3; i++)
      data.Append(v_axis(i));
    data.Append((isfirst) ? 1. : 0.);
    data.Append((islast) ? 1. : 0.);
  }



  Revolution :: Revolution(const Point<3> & p0_in,
			   const Point<3> & p1_in,
			   const SplineGeometry2d & spline_in) :
    p0(p0_in), p1(p1_in), splinecurve(spline_in),
    nsplines(spline_in.GetNSplines())
  {
    surfaceactive.SetSize(0);
    surfaceids.SetSize(0);

    v_axis = p1-p0;

    v_axis.Normalize();

    if(splinecurve.GetSpline(0).StartPI()(1) <= 0. &&
       splinecurve.GetSpline(nsplines-1).EndPI()(1) <= 0.)
      type = 2;
    else if (Dist(splinecurve.GetSpline(0).StartPI(),
		  splinecurve.GetSpline(nsplines-1).EndPI()) < 1e-7)
      type = 1;
    else
      cerr << "Surface of revolution cannot be constructed" << endl;

    for(int i=0; i<splinecurve.GetNSplines(); i++)
      {
	RevolutionFace * face = new RevolutionFace(splinecurve.GetSpline(i),
						   p0,v_axis,
						   type==2 && i==0,
						   type==2 && i==splinecurve.GetNSplines()-1);
	faces.Append(face);
	surfaceactive.Append(1);
	surfaceids.Append(0);
      }
  }
  
  Revolution::~Revolution()
  {
    for(int i=0; i<faces.Size(); i++)
      delete faces[i];
  }


  INSOLID_TYPE Revolution :: BoxInSolid (const BoxSphere<3> & box) const
  {
    for(int i=0; i<faces.Size(); i++)
      if(faces[i]->BoxIntersectsFace(box))
	return DOES_INTERSECT;
    
    
    return PointInSolid(box.Center(),0);
	 

    /*
    Point<2> c,pmin,pmax;
    faces[0]->CalcProj(box.Center(),c);
    double aux = box.Diam()/sqrt(8.);
    pmin(0) = c(0)-aux; pmin(1) = c(1)-aux;
    pmax(0) = c(0)+aux; pmax(1) = c(1)+aux;
    

    BoxSphere<2> box2d(pmin,pmax);

    bool intersection = false;
    bool uncertain = true;

    for(int i=0; !(intersection && !uncertain) && i<faces.Size(); i++)
      {
	bool thisintersects;
	bool thisuncertain;
	thisintersects = faces[i]->BoxIntersectsFace(box2d,thisuncertain);
	intersection = intersection || thisintersects;
	if(thisintersects && !thisuncertain)
	  uncertain = false;
      }

    if(intersection)
      {
	if(!uncertain)
	  return DOES_INTERSECT;
	else
	  {
	    ARRAY < Point<3> > pext(2);
	    Point<3> p;

	    pext[0] = box.PMin();
	    pext[1] = box.PMax();

	    INSOLID_TYPE position;
	    bool firsttime = true;

	    for(int i=0; i<2; i++)
	      for(int j=0; j<2; j++)
		for(int k=0; k<2; k++)
		  {
		    p(0) = pext[i](0);
		    p(1) = pext[j](1);
		    p(2) = pext[k](2);
		    INSOLID_TYPE ppos = PointInSolid(p,0);
		    if(ppos == DOES_INTERSECT)
		      return DOES_INTERSECT;
		    
		    if(firsttime)
		      {
			firsttime = false;
			position = ppos;
		      }
		    if(position != ppos)
		      return DOES_INTERSECT;	    
		  }
	    return position;

	  }
      }

    return PointInSolid(box.Center(),0);
    */ 
  }

  INSOLID_TYPE Revolution :: PointInSolid (const Point<3> & p,
					   double eps) const
  {
    Point<2> p2d;
    faces[0]->CalcProj(p,p2d);

    int intersections_before(0), intersections_after(0);
    double randomx = 7.42357;
    double randomy = 1.814756;
    randomx *= 1./sqrt(randomx*randomx+randomy*randomy);
    randomy *= 1./sqrt(randomx*randomx+randomy*randomy);
    

    const double a = randomy;
    const double b = -randomx;
    const double c = -a*p2d(0)-b*p2d(1);

    ARRAY < Point<2> > points;

    //(*testout) << "face intersections at: " << endl;
    for(int i=0; i<faces.Size(); i++)
      {
	faces[i]->GetSpline().LineIntersections(a,b,c,points,eps);
	
	for(int j=0; j<points.Size(); j++)
	  {
	    double t = (points[j](0)-p2d(0))/randomx;

	    //(*testout) << t << endl;
	    if ( t < -eps )
	      intersections_before++;
	    else if ( t > eps )
	      intersections_after++;
	    else
	      {
		intersecting_face = i;
		return DOES_INTERSECT;
	      }
	  }
      }

    if(intersections_before % 2 == 0)
      return IS_OUTSIDE;
    else
      return IS_INSIDE;
  }

  INSOLID_TYPE Revolution :: VecInSolid (const Point<3> & p,
					 const Vec<3> & v,
					 double eps) const
  {
    INSOLID_TYPE pInSolid = PointInSolid(p,eps);

    if(pInSolid != DOES_INTERSECT)
      {
	//(*testout) << "pInSolid" << endl;
	return pInSolid;
      }

    ARRAY<int> intersecting_faces;

    for(int i=0; i<faces.Size(); i++)
      if(faces[i]->PointInFace(p,eps) == DOES_INTERSECT)
	intersecting_faces.Append(i);

     Vec<3> hv;

    if(intersecting_faces.Size() == 1)
      {
	faces[intersecting_faces[0]]->CalcGradient(p,hv);

	double hv1;
	hv1 = v * hv;
	
	if (hv1 <= -eps)
	  return IS_INSIDE;
	if (hv1 >= eps)
	  return IS_OUTSIDE;
	
	return DOES_INTERSECT; 
      }
    else if(intersecting_faces.Size() == 2)
      {
	Point<2> p2d;
	Vec<2> v2d;
	faces[intersecting_faces[0]]->CalcProj(p,p2d,v,v2d);

	if(Dist(faces[intersecting_faces[0]]->GetSpline().StartPI(),p2d) <
	   Dist(faces[intersecting_faces[0]]->GetSpline().EndPI(),p2d))
	  {
	    int aux = intersecting_faces[0];
	    intersecting_faces[0] = intersecting_faces[1];
	    intersecting_faces[1] = aux;
	  }
	
	const SplineSeg3<2> * splinesegment3 = 
	  dynamic_cast<const SplineSeg3<2> *>(&faces[intersecting_faces[0]]->GetSpline());
	const LineSeg<2> * linesegment = 
	  dynamic_cast<const LineSeg<2> *>(&faces[intersecting_faces[0]]->GetSpline());
		
	Vec<2> t1,t2;

	if(linesegment)
	  t1 = linesegment->StartPI() - linesegment->EndPI();
	else if(splinesegment3)
	  t1 = splinesegment3->TangentPoint() - splinesegment3->EndPI();

	linesegment = 
	  dynamic_cast<const LineSeg<2> *>(&faces[intersecting_faces[1]]->GetSpline());
	splinesegment3 = 
	  dynamic_cast<const SplineSeg3<2> *>(&faces[intersecting_faces[1]]->GetSpline());
	
	if(linesegment)
	  t2 = linesegment->EndPI() - linesegment->StartPI();
	else if(splinesegment3)
	  t2 = splinesegment3->TangentPoint() - splinesegment3->StartPI();

	t1.Normalize();
	t2.Normalize();

	double d1 = v2d*t1;
	double d2 = v2d*t2;

	Vec<2> n;

	if(d1 > d2)
	  {
	    n(0) = t1(1);
	    n(1) = -t1(0);
	  }
	else
	  {
	    n(0) = -t2(1);
	    n(1) = t2(0);
	  }

	double d = v2d*n;

	if(d > eps)
	  return IS_OUTSIDE;
	else if (d < -eps)
	  return IS_INSIDE;
	else
	  return DOES_INTERSECT;


      }
    else
      {
	cerr << "Jo gibt's denn des?" << endl;
      }

    return DOES_INTERSECT;    
  }

  INSOLID_TYPE Revolution :: VecInSolid2 (const Point<3> & p,
					  const Vec<3> & v1,
					  const Vec<3> & v2,
					  double eps) const
  {
    INSOLID_TYPE ret1 = VecInSolid(p,v1,eps);
    if(ret1 != DOES_INTERSECT)
      return ret1;

    return VecInSolid(p,v2,eps);
  }
  
  int Revolution :: GetNSurfaces() const
  {
    return faces.Size();
  }

  Surface & Revolution :: GetSurface (int i)
  {
    return *faces[i];
  }

  const Surface & Revolution :: GetSurface (int i) const
  {
    return *faces[i];
  }


  void Revolution :: Reduce (const BoxSphere<3> & box)
  { 
    //bool dummy;
    for(int i=0; i<faces.Size(); i++)
      surfaceactive[i] = (faces[i]->BoxIntersectsFace(box));
    //surfaceactive[i] = (faces[i]->BoxIntersectsFace(box,dummy));
  }

  void Revolution :: UnReduce ()
  {
    for(int i=0; i<faces.Size(); i++)
      surfaceactive[i] = true;
  }
}
