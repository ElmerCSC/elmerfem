#include <mystdlib.h>


#include <linalg.hpp>
#include <csg.hpp>


namespace netgen
{

double 
QuadraticSurface :: CalcFunctionValue (const Point<3> & p) const
{
  return p(0) * (cxx * p(0) + cxy * p(1) + cxz * p(2) + cx) +
    p(1) * (cyy * p(1) + cyz * p(2) + cy) +
    p(2) * (czz * p(2) + cz) + c1;
}

void 
QuadraticSurface :: CalcGradient (const Point<3> & p, Vec<3> & grad) const
{
  grad(0) = 2 * cxx * p(0) + cxy * p(1) + cxz * p(2) + cx;
  grad(1) = 2 * cyy * p(1) + cxy * p(0) + cyz * p(2) + cy;
  grad(2) = 2 * czz * p(2) + cxz * p(0) + cyz * p(1) + cz;
}

void 
QuadraticSurface :: CalcHesse (const Point<3> & /* p */, Mat<3> & hesse) const
{
  hesse(0,0) = 2 * cxx;
  hesse(1,1) = 2 * cyy;
  hesse(2,2) = 2 * czz;
  hesse(0,1) = hesse(1,0) = cxy;
  hesse(0,2) = hesse(2,0) = cxz;
  hesse(1,2) = hesse(2,1) = cyz;
}


void QuadraticSurface :: Read (istream & ist)
{
  ist >> cxx >> cyy >> czz >> cxy >> cxz >> cyz >> cx >> cy >> cz >> c1;
}

void QuadraticSurface :: Print (ostream & ost) const
{
  ost << cxx << "  " << cyy << "  " << czz << "  "
      << cxy << "  " << cxz << "  " << cyz << "  "
      << cx << "  " << cy << "  " << cz << "  "
      << c1;
}


void QuadraticSurface :: PrintCoeff (ostream & ost) const
{
  ost << " cxx = " << cxx
      << " cyy = " << cyy
      << " czz = " << czz
      << " cxy = " << cxy
      << " cxz = " << cxz
      << " cyz = " << cyz
      << " cx = " << cx
      << " cy = " << cy
      << " cz = " << cz
      << " c1 = " << c1 << endl;
}



Point<3> QuadraticSurface :: GetSurfacePoint () const
{
  MyError ("GetSurfacePoint called for QuadraticSurface");
  return Point<3> (0, 0, 0);
}


Plane :: Plane (const Point<3> & ap, Vec<3> an)
{
  eps_base = 1e-8;

  p = ap;
  n = an;
  n.Normalize();

  cxx = cyy = czz = cxy = cxz = cyz = 0;
  cx = n(0); cy = n(1); cz = n(2);
  c1 = - (cx * p(0) + cy * p(1) + cz * p(2));
}

Primitive * Plane :: Copy () const
{
  return new Plane (p, n);
}

void Plane :: Transform (Transformation<3> & trans)
{
  Point<3> hp;
  Vec<3> hn;
  trans.Transform (p, hp);
  trans.Transform (n, hn);
  p = hp;
  n = hn;

  cxx = cyy = czz = cxy = cxz = cyz = 0;
  cx = n(0); cy = n(1); cz = n(2);
  c1 = - (cx * p(0) + cy * p(1) + cz * p(2));
}



void Plane :: GetPrimitiveData (const char *& classname, 
				ARRAY<double> & coeffs) const
{
  classname = "plane";
  coeffs.SetSize (6);
  coeffs.Elem(1) = p(0);
  coeffs.Elem(2) = p(1);
  coeffs.Elem(3) = p(2);
  coeffs.Elem(4) = n(0);
  coeffs.Elem(5) = n(1);
  coeffs.Elem(6) = n(2);
}

void Plane :: SetPrimitiveData (ARRAY<double> & coeffs)
{
  p(0) = coeffs.Elem(1);
  p(1) = coeffs.Elem(2);
  p(2) = coeffs.Elem(3);
  n(0) = coeffs.Elem(4);
  n(1) = coeffs.Elem(5);
  n(2) = coeffs.Elem(6);

  n.Normalize();

  cxx = cyy = czz = cxy = cxz = cyz = 0;
  cx = n(0); cy = n(1); cz = n(2);
  c1 = - (cx * p(0) + cy * p(1) + cz * p(2));
}

Primitive * Plane :: CreateDefault ()
{
  return new Plane (Point<3> (0,0,0), Vec<3> (0,0,1));
}


int Plane :: IsIdentic (const Surface & s2, int & inv, double eps) const
{
  const Plane * ps2 = dynamic_cast<const Plane*>(&s2);

  if(ps2)
    {
      Point<3> pp2 = ps2->GetSurfacePoint();
      Vec<3> n2 = s2.GetNormalVector(pp2);

      if(fabs(n*n2) < 1.-eps_base)
	return 0;

      if (fabs (s2.CalcFunctionValue(p)) > eps) return 0;
    }
  else
    {
      //(*testout) << "pos1" << endl;
      if (fabs (s2.CalcFunctionValue(p)) > eps) return 0;
      Vec<3> hv1, hv2;
      hv1 = n.GetNormal ();
      hv2 = Cross (n, hv1);
      
      Point<3> hp = p + hv1;
      //(*testout) << "pos2" << endl;
      //(*testout) << "eps " << eps << " s2.CalcFunctionValue(hp) " << s2.CalcFunctionValue(hp) << endl;
      if (fabs (s2.CalcFunctionValue(hp)) > eps) return 0;
      //(*testout) << "pos3" << endl;
      hp = p + hv2;
      if (fabs (s2.CalcFunctionValue(hp)) > eps) return 0;
    }

  //(*testout) << "pos4" << endl;
  Vec<3> n1, n2;
  n1 = GetNormalVector (p);
  n2 = s2.GetNormalVector (p);
  inv = (n1 * n2 < 0);
  //(*testout) << "inv " << inv << endl;
  return 1;
}



void Plane :: DefineTangentialPlane (const Point<3> & ap1, const Point<3> & ap2)
{
  Surface::DefineTangentialPlane (ap1, ap2);
}


void Plane :: ToPlane (const Point<3> & p3d, 
		       Point<2> & pplane, 
		       double h, int & zone) const
{
  Vec<3> p1p;

  p1p = p3d - p1;
  p1p /= h;
  pplane(0) = p1p * ex;
  pplane(1) = p1p * ey;
  zone = 0;
}

void Plane :: FromPlane (const Point<2> & pplane, Point<3> & p3d, double h) const
{
  /*
  Vec<3> p1p;
  Point<2> pplane2 = pplane;
  
  pplane2 *= h;
  p1p = pplane2(0) * ex + pplane2(1) * ey;
  p3d = p1 + p1p;
  */
  p3d = p1 + (h * pplane(0)) * ex + (h * pplane(1)) * ey;
}


void Plane :: Project (Point<3> & p3d) const
{
  double val = Plane::CalcFunctionValue (p3d);
  p3d -= val * n;
}

INSOLID_TYPE Plane :: BoxInSolid (const BoxSphere<3> & box) const
{
  int i;
  double val;
  Point<3> pp;

  val = Plane::CalcFunctionValue (box.Center());
  if (val > box.Diam() / 2) return IS_OUTSIDE;
  if (val < -box.Diam() / 2) return IS_INSIDE;

  if (val > 0)
    {
      /*
      double modify = 
	((box.MaxX()-box.MinX()) * fabs (cx) + 
	 (box.MaxY()-box.MinY()) * fabs (cy) + 
	 (box.MaxZ()-box.MinZ()) * fabs (cz)) / 2;
      */
      Vec<3> vdiag = box.PMax() - box.PMin();
      double modify = (vdiag(0) * fabs (cx) + 
		       vdiag(1) * fabs (cy) + 
		       vdiag(2) * fabs (cz) ) / 2;

      if (val - modify < 0)
	return DOES_INTERSECT;
      return IS_OUTSIDE;

      // only outside or intersect possible
      for (i = 0; i < 8; i++)
	{
	  pp = box.GetPointNr (i);
	  val = Plane::CalcFunctionValue (pp);
	  if (val < 0) 
	    return DOES_INTERSECT;
	}
      return IS_OUTSIDE;
    }
  else
    {
      /*
	double modify = 
	((box.MaxX()-box.MinX()) * fabs (cx) + 
	(box.MaxY()-box.MinY()) * fabs (cy) + 
	(box.MaxZ()-box.MinZ()) * fabs (cz)) / 2;
      */
      Vec<3> vdiag = box.PMax() - box.PMin();
      double modify =  (vdiag(0) * fabs (cx) + 
			vdiag(1) * fabs (cy) + 
			vdiag(2) * fabs (cz) ) / 2;
      if (val + modify > 0)
	return DOES_INTERSECT;
      return IS_INSIDE;


      // only inside or intersect possible
      for (i = 0; i < 8; i++)
	{
	  pp = box.GetPointNr (i);
	  val = Plane::CalcFunctionValue (pp);
	  if (val > 0) 
	    return DOES_INTERSECT;
	}
      return IS_INSIDE;
    }



  /*
  for (i = 1; i <= 8; i++)
    {
      box.GetPointNr (i, p);
      val = CalcFunctionValue (p);
      if (val > 0) inside = 0;
      if (val < 0) outside = 0;
    }

  if (inside) return IS_INSIDE;
  if (outside) return IS_OUTSIDE;
  return DOES_INTERSECT;
  */
}



// double Plane :: CalcFunctionValue (const Point<3> & p3d) const
// {
//   return cx * p3d(0) + cy * p3d(1) + cz * p3d(2) + c1;
// }

void Plane :: CalcGradient (const Point<3> & /* p */, Vec<3> & grad) const
{
  grad(0) = cx;
  grad(1) = cy;
  grad(2) = cz;
}

void Plane :: CalcHesse (const Point<3> & /* p */, Mat<3> & hesse) const
{
  hesse = 0;
}

double Plane :: HesseNorm () const
{
  return 0;
}


Point<3> Plane :: GetSurfacePoint () const
{
  return p;
}


void Plane :: GetTriangleApproximation 
(TriangleApproximation & tas, 
 const Box<3> & boundingbox, double facets) const
{
  // find triangle, such that
  // boundingbox /cap plane is contained in it

  Point<3> c = boundingbox.Center();
  double r = boundingbox.Diam();

  Project (c);
  Vec<3> t1 = n.GetNormal();
  Vec<3> t2 = Cross (n, t1);

  t1.Normalize();
  t2.Normalize();

  tas.AddPoint (c + (-0.5 * r) * t2 + (sqrt(0.75) * r) * t1);
  tas.AddPoint (c + (-0.5 * r) * t2 + (-sqrt(0.75) * r) * t1);
  tas.AddPoint (c +  r * t2);

  tas.AddTriangle (TATriangle (0, 0, 1, 2));
}




Sphere :: Sphere (const Point<3> & ac, double ar)
{
  c = ac;
  r = ar;
  
  cxx = cyy = czz = 0.5 / r;
  cxy = cxz = cyz = 0;
  cx = - c(0) / r;
  cy = - c(1) / r;
  cz = - c(2) / r;
  c1 = (c(0) * c(0) + c(1) * c(1) + c(2) * c(2)) / (2 * r) - r / 2;
}

void Sphere :: GetPrimitiveData (const char *& classname, ARRAY<double> & coeffs) const
{
  classname = "sphere";
  coeffs.SetSize (4);
  coeffs.Elem(1) = c(0);
  coeffs.Elem(2) = c(1);
  coeffs.Elem(3) = c(2);
  coeffs.Elem(4) = r;
}

void Sphere :: SetPrimitiveData (ARRAY<double> & coeffs)
{
  c(0) = coeffs.Elem(1);
  c(1) = coeffs.Elem(2);
  c(2) = coeffs.Elem(3);
  r = coeffs.Elem(4);

  cxx = cyy = czz = 0.5 / r;
  cxy = cxz = cyz = 0;
  cx = - c(0) / r;
  cy = - c(1) / r;
  cz = - c(2) / r;
  c1 = (c(0) * c(0) + c(1) * c(1) + c(2) * c(2)) / (2 * r) - r / 2;
}

Primitive * Sphere :: CreateDefault ()
{
  return new Sphere (Point<3> (0,0,0), 1);
}



Primitive * Sphere :: Copy () const
{
  return new Sphere (c, r);
}

void Sphere :: Transform (Transformation<3> & trans)
{
  Point<3> hp;
  trans.Transform (c, hp);
  c = hp;

  cxx = cyy = czz = 0.5 / r;
  cxy = cxz = cyz = 0;
  cx = - c(0) / r;
  cy = - c(1) / r;
  cz = - c(2) / r;
  c1 = (c(0) * c(0) + c(1) * c(1) + c(2) * c(2)) / (2 * r) - r / 2;
}




int Sphere :: IsIdentic (const Surface & s2, int & inv, double eps) const
{
  const Sphere * sp2 = dynamic_cast<const Sphere*>  (&s2);

  if (!sp2) return 0;

  if (Dist (sp2->c, c) > eps) return 0;
  if (fabs (sp2->r - r) > eps) return 0;

  inv = 0;

  return 1;
}


void Sphere :: DefineTangentialPlane (const Point<3> & ap1, const Point<3> & ap2)
{
  Surface::DefineTangentialPlane (ap1, ap2);

  ez = p1 - c;
  ez /= ez.Length();

  ex = p2 - p1;
  ex -= (ex * ez) * ez;
  ex /= ex.Length();

  ey = Cross (ez, ex);
}


void Sphere :: ToPlane (const Point<3> & p, Point<2> & pplane, double h, int & zone) const
{
  Vec<3> p1p;
  
  p1p = p - p1;
  
  /*
  if (p1p * ez < -r)
    {
      zone = -1;
      pplane = Point<2> (1E8, 1E8);
    }
  else
    { 
      zone = 0;
      p1p /= h;
      pplane(0) = p1p * ex;
      pplane(1) = p1p * ey;
    }
  */

  Point<3> p1top = c + (c - p1);

  Vec<3> p1topp = p - p1top;
  Vec<3> p1topp1 = p1 - p1top;
  Vec<3> lam;
  //  SolveLinearSystem (ex, ey, p1topp, p1topp1, lam);

  Mat<3> m;
  for (int i = 0; i < 3; i++)
    {
      m(i, 0) = ex(i);
      m(i, 1) = ey(i);
      m(i, 2) = p1topp(i);
    }
  m.Solve (p1topp1, lam);

  pplane(0) = -lam(0) / h;
  pplane(1) = -lam(1) / h;
  
  if (lam(2) > 2)
    zone = -1;
  else 
    zone = 0;
}

void Sphere :: FromPlane (const Point<2> & pplane, Point<3> & p, double h) const
{
  /*
    //  Vec<3> p1p;
    double z;
    Point<2> pplane2 (pplane);

    pplane2(0) *= h;
    pplane2(1) *= h;
    z = -r + sqrt (sqr (r) - sqr (pplane2(0)) - sqr (pplane2(1)));
    //  p = p1;
    p(0) = p1(0) + pplane2(0) * ex(0) + pplane2(1) * ey(0) + z * ez(0);
    p(1) = p1(1) + pplane2(0) * ex(1) + pplane2(1) * ey(1) + z * ez(1);
    p(2) = p1(2) + pplane2(0) * ex(2) + pplane2(1) * ey(2) + z * ez(2);
    */

  Point<2> pplane2 (pplane);

  pplane2(0) *= h;
  pplane2(1) *= h;

  p(0) = p1(0) + pplane2(0) * ex(0) + pplane2(1) * ey(0);
  p(1) = p1(1) + pplane2(0) * ex(1) + pplane2(1) * ey(1);
  p(2) = p1(2) + pplane2(0) * ex(2) + pplane2(1) * ey(2);
  Project (p);
}


void Sphere :: Project (Point<3> & p) const
{
  Vec<3> v;
  v = p - c;
  v *= (r / v.Length());
  p = c + v;
}


INSOLID_TYPE Sphere :: BoxInSolid (const BoxSphere<3> & box) const
{
  double dist;
  dist = Dist (box.Center(), c);

  if (dist - box.Diam()/2 > r) return IS_OUTSIDE;
  if (dist + box.Diam()/2 < r) return IS_INSIDE;
  return DOES_INTERSECT;
}

double Sphere :: HesseNorm () const
{
  return 2 / r;
}


Point<3> Sphere :: GetSurfacePoint () const
{
  return c + Vec<3> (r, 0, 0);
}


void Sphere :: GetTriangleApproximation 
(TriangleApproximation & tas, 
 const Box<3> & boundingbox, double facets) const
{
  int i, j;
  double lg, bg;
  int n = int(facets) + 1;  

  for (j = 0; j <= n; j++)
    for (i = 0; i <= n; i++)
      {
	lg = 2 * M_PI * double (i) / n;
	bg = M_PI * (double(j) / n - 0.5);

	Point<3> p(c(0) + r * cos(bg) * sin (lg),
		  c(1) + r * cos(bg) * cos (lg),
		  c(2) + r * sin(bg));
	tas.AddPoint (p);
      }

  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++)
      {
	int pi = i + (n+1) * j;
	tas.AddTriangle (TATriangle (0, pi, pi+1, pi+n+2));
	tas.AddTriangle (TATriangle (0, pi, pi+n+2, pi+n+1));
      }
}





Ellipsoid :: 
Ellipsoid (const Point<3> & aa,
	   const Vec<3> & av1, const Vec<3> & av2, const Vec<3> & av3)
{
  a = aa;
  v1 = av1;
  v2 = av2;
  v3 = av3;

  CalcData();
}


void Ellipsoid :: CalcData ()
{
  // f = (x-a, vl)^2 / |vl|^2 + (x-a, vs)^2 / |vs|^2 -1
  // f = sum_{i=1}^3  (x-a,v_i)^2 / |vi|^4 - 1   =  sum (x-a,hv_i)^2
  
  Vec<3> hv1, hv2, hv3;
  double lv1 = v1.Length2 ();
  if (lv1 < 1e-32) lv1 = 1;
  double lv2 = v2.Length2 ();
  if (lv2 < 1e-32) lv2 = 1;
  double lv3 = v3.Length2 ();
  if (lv3 < 1e-32) lv3 = 1;

  rmin = sqrt (min3 (lv1, lv2, lv3));

  hv1 = (1.0 / lv1) * v1;
  hv2 = (1.0 / lv2) * v2;
  hv3 = (1.0 / lv3) * v3;

  cxx = hv1(0) * hv1(0) + hv2(0) * hv2(0) + hv3(0) * hv3(0);
  cyy = hv1(1) * hv1(1) + hv2(1) * hv2(1) + hv3(1) * hv3(1);
  czz = hv1(2) * hv1(2) + hv2(2) * hv2(2) + hv3(2) * hv3(2);

  cxy = 2 * (hv1(0) * hv1(1) + hv2(0) * hv2(1) + hv3(0) * hv3(1));
  cxz = 2 * (hv1(0) * hv1(2) + hv2(0) * hv2(2) + hv3(0) * hv3(2));
  cyz = 2 * (hv1(1) * hv1(2) + hv2(1) * hv2(2) + hv3(1) * hv3(2));

  Vec<3> va (a);
  c1 = sqr(va * hv1) + sqr(va * hv2) + sqr(va * hv3) - 1;
  
  Vec<3> v = -2 * (va * hv1) * hv1 - 2 * (va * hv2) * hv2  - 2 * (va * hv3) * hv3;
  cx = v(0);
  cy = v(1);
  cz = v(2);
}


INSOLID_TYPE Ellipsoid :: BoxInSolid (const BoxSphere<3> & box) const
{
  // double grad = 2.0 / rmin;
  // double grad = 3*(box.Center()-a).Length() / (rmin*rmin*rmin);

  double ggrad = 1.0 / (rmin*rmin);
  Vec<3> g;
  double val = CalcFunctionValue (box.Center());
  CalcGradient (box.Center(), g);
  double grad = g.Length();

  double r = box.Diam() / 2;
  double maxval = grad * r + ggrad * r * r;

  //  (*testout) << "box = " << box << ", val = " << val << ", maxval = " << maxval << endl;

  if (val > maxval) return IS_OUTSIDE;
  if (val < -maxval) return IS_INSIDE;
  return DOES_INTERSECT;
}


double Ellipsoid :: HesseNorm () const
{
  return 1.0/ (rmin * rmin);
}

double Ellipsoid :: MaxCurvature () const
{
  const double a2 = v1.Length2();
  const double b2 = v2.Length2();
  const double c2 = v3.Length2();

  return max3 ( sqrt(a2)/min2(b2,c2), sqrt(b2)/min2(a2,c2), sqrt(c2)/min2(a2,b2) );
}

Point<3> Ellipsoid :: GetSurfacePoint () const
{
  return a + v1;
}



void Ellipsoid :: GetTriangleApproximation 
(TriangleApproximation & tas, 
 const Box<3> & boundingbox, double facets) const
{
  int i, j;
  double lg, bg;
  int n = int(facets) + 1;  

  for (j = 0; j <= n; j++)
    for (i = 0; i <= n; i++)
      {
	lg = 2 * M_PI * double (i) / n;
	bg = M_PI * (double(j) / n - 0.5);


	Point<3> p(a + 
		   sin (bg) * v1 + 
		   cos (bg) * sin (lg) * v2 +
		   cos (bg) * cos (lg) * v3);

	tas.AddPoint (p);
      }

  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++)
      {
	int pi = i + (n+1) * j;
	tas.AddTriangle (TATriangle (0, pi, pi+1, pi+n+2));
	tas.AddTriangle (TATriangle (0, pi, pi+n+2, pi+n+1));
      }
}

















Cylinder :: Cylinder (ARRAY<double> & coeffs)
{
  SetPrimitiveData(coeffs);
}

Cylinder :: Cylinder (const Point<3> & aa, const Point<3> & ab, double ar)
{
  a = aa;
  b = ab;
  vab = (b - a);
  vab /= vab.Length();
  r = ar;

  // ( <x,x> - 2 <x,a> + <a,a>
  //   - <x,vab>^2 + 2 <x,vab> <a, vab> - <a, vab>^2
  //   - r^2) / (2r) = 0

  double hv;
  cxx = cyy = czz = 0.5 / r;
  cxy = cxz = cyz = 0;
  cx = - a(0) / r;
  cy = - a(1) / r;
  cz = - a(2) / r;
  c1 = (a(0) * a(0) + a(1) * a(1) + a(2) * a(2)) / (2 * r);
  hv = a(0) * vab(0) + a(1) * vab(1) + a(2) * vab(2);
  cxx -= vab(0) * vab(0) / (2 * r);
  cyy -= vab(1) * vab(1) / (2 * r);
  czz -= vab(2) * vab(2) / (2 * r);
  cxy -= vab(0) * vab(1) / r;
  cxz -= vab(0) * vab(2) / r;
  cyz -= vab(1) * vab(2) / r;
  cx += vab(0) * hv / r;
  cy += vab(1) * hv / r;
  cz += vab(2) * hv / r;
  c1 -= hv * hv / (2 * r);
  c1 -= r / 2;
  //  PrintCoeff ();
}




void Cylinder :: GetPrimitiveData (const char *& classname, ARRAY<double> & coeffs) const
{
  classname = "cylinder";
  coeffs.SetSize (7);
  coeffs.Elem(1) = a(0);
  coeffs.Elem(2) = a(1);
  coeffs.Elem(3) = a(2);
  coeffs.Elem(4) = b(0);
  coeffs.Elem(5) = b(1);
  coeffs.Elem(6) = b(2);
  coeffs.Elem(7) = r;
}

void Cylinder :: SetPrimitiveData (ARRAY<double> & coeffs)
{
  a(0) = coeffs.Elem(1);
  a(1) = coeffs.Elem(2);
  a(2) = coeffs.Elem(3);
  b(0) = coeffs.Elem(4);
  b(1) = coeffs.Elem(5);
  b(2) = coeffs.Elem(6);
  r = coeffs.Elem(7);


  vab = (b - a);
  vab /= vab.Length();


  double hv;
  cxx = cyy = czz = 0.5 / r;
  cxy = cxz = cyz = 0;
  cx = - a(0) / r;
  cy = - a(1) / r;
  cz = - a(2) / r;
  c1 = (a(0) * a(0) + a(1) * a(1) + a(2) * a(2)) / (2 * r);
  hv = a(0) * vab(0) + a(1) * vab(1) + a(2) * vab(2);
  cxx -= vab(0) * vab(0) / (2 * r);
  cyy -= vab(1) * vab(1) / (2 * r);
  czz -= vab(2) * vab(2) / (2 * r);
  cxy -= vab(0) * vab(1) / r;
  cxz -= vab(0) * vab(2) / r;
  cyz -= vab(1) * vab(2) / r;
  cx += vab(0) * hv / r;
  cy += vab(1) * hv / r;
  cz += vab(2) * hv / r;
  c1 -= hv * hv / (2 * r);
  c1 -= r / 2;
}

Primitive * Cylinder :: CreateDefault ()
{
  return new Cylinder (Point<3> (0,0,0), Point<3> (1,0,0), 1);
}




Primitive * Cylinder :: Copy () const
{
  return new Cylinder (a, b, r);
}


int Cylinder :: IsIdentic (const Surface & s2, int & inv, double eps) const
{
  const Cylinder * cyl2 = dynamic_cast<const Cylinder*>  (&s2);

  if (!cyl2) return 0;

  if (fabs (cyl2->r - r) > eps) return 0;

  Vec<3> v1 = b - a;
  Vec<3> v2 = cyl2->a - a;

  if ( fabs (v1 * v2) < (1-eps) * v1.Length() * v2.Length()) return 0;
  v2 = cyl2->b - a;
  if ( fabs (v1 * v2) < (1-eps) * v1.Length() * v2.Length()) return 0;

  inv = 0;
  return 1;
}



void Cylinder :: Transform (Transformation<3> & trans)
{
  Point<3> hp;
  trans.Transform (a, hp);
  a = hp;
  trans.Transform (b, hp);
  b = hp;

  vab = (b - a);
  vab /= vab.Length();

  // ( <x,x> - 2 <x,a> + <a,a>
  //   - <x,vab>^2 + 2 <x,vab> <a, vab> - <a, vab>^2
  //   - r^2) / (2r) = 0

  double hv;
  cxx = cyy = czz = 0.5 / r;
  cxy = cxz = cyz = 0;
  cx = - a(0) / r;
  cy = - a(1) / r;
  cz = - a(2) / r;
  c1 = (a(0) * a(0) + a(1) * a(1) + a(2) * a(2)) / (2 * r);
  hv = a(0) * vab(0) + a(1) * vab(1) + a(2) * vab(2);
  cxx -= vab(0) * vab(0) / (2 * r);
  cyy -= vab(1) * vab(1) / (2 * r);
  czz -= vab(2) * vab(2) / (2 * r);
  cxy -= vab(0) * vab(1) / r;
  cxz -= vab(0) * vab(2) / r;
  cyz -= vab(1) * vab(2) / r;
  cx += vab(0) * hv / r;
  cy += vab(1) * hv / r;
  cz += vab(2) * hv / r;
  c1 -= hv * hv / (2 * r);
  c1 -= r / 2;
  //  PrintCoeff ();
}









void Cylinder :: DefineTangentialPlane (const Point<3> & ap1, const Point<3> & ap2)
{
  Surface::DefineTangentialPlane (ap1, ap2);

  ez = Center (p1, p2) - a;
  ez -= (ez * vab) * vab;
  ez /= ez.Length();

  ex = p2 - p1;
  ex -= (ex * ez) * ez;
  ex /= ex.Length();

  ey = Cross (ez, ex);
}


void Cylinder :: ToPlane (const Point<3> & p, 
			  Point<2> & pplane, 
			  double h, int & zone) const
{
  Point<3> cp1p2 = Center (p1, p2);
  Project (cp1p2);
  
  Point<3> ccp1p2 = a + ( (cp1p2 - a) * vab ) * vab;

  Vec<3> er = cp1p2 - ccp1p2;
  er.Normalize();
  Vec<3> ephi = Cross (vab, er);

  double co, si;
  Point<2> p1p, p2p, pp;

  co = er * (p1 - ccp1p2);
  si = ephi * (p1 - ccp1p2);
  p1p(0) = r * atan2 (si, co);
  p1p(1) = vab * (p1 - ccp1p2);

  co = er * (p2 - ccp1p2);
  si = ephi * (p2 - ccp1p2);
  p2p(0) = r * atan2 (si, co);
  p2p(1) = vab * (p2 - ccp1p2);
  
  co = er * (p - ccp1p2);
  si = ephi * (p - ccp1p2);

  double phi = atan2 (si, co);
  pp(0) = r * phi;
  pp(1) = vab * (p - ccp1p2);
  
  zone = 0;
  if (phi > 1.57) zone = 1;
  if (phi < -1.57) zone = 2;



  Vec<2> e2x = p2p - p1p;
  e2x /= e2x.Length();

  Vec<2> e2y (-e2x(1), e2x(0));

  Vec<2> p1pp = pp - p1p;


  pplane(0) = (p1pp * e2x) / h;
  pplane(1) = (p1pp * e2y) / h;

  /*
  (*testout) << "p1 = " << p1 << ",  p2 = " << p2 << endl;
  (*testout) << "p = " << p << ",  pp = " << pp << ",  pplane = " << pplane << endl;
  */

  /*
  Vec<3> p1p;

  p1p = p - p1;

  if (p1p * ez < -1 * r)
    {
      zone = -1;
      pplane(0) = 1e8;
      pplane(1) = 1e8;
    }
  else
    {
      zone = 0;
      p1p /= h;
      pplane(0) = p1p * ex;
      pplane(1) = p1p * ey;
    }
    */
}

void Cylinder :: FromPlane (const Point<2> & pplane, Point<3> & p, double h) const
{
  Point<2> pplane2 (pplane);

  pplane2(0) *= h;
  pplane2(1) *= h;

  p(0) = p1(0) + pplane2(0) * ex(0) + pplane2(1) * ey(0);
  p(1) = p1(1) + pplane2(0) * ex(1) + pplane2(1) * ey(1);
  p(2) = p1(2) + pplane2(0) * ex(2) + pplane2(1) * ey(2);
  Project (p);
}


void Cylinder :: Project (Point<3> & p) const
{
  Vec<3> v;
  Point<3> c;

  c = a + ((p - a) * vab) * vab;
  v = p - c;
  v *= (r / v.Length());
  p = c + v;
}
/*
int Cylinder :: RootInBox (const BoxSphere<3> & box) const
  {
  double dist;
  dist = sqrt (2 * CalcFunctionValue(box.Center()) * r + r * r);
  if (fabs (dist - r) > box.Diam()/2) return 0;
  return 2;
  }
*/

INSOLID_TYPE Cylinder :: BoxInSolid (const BoxSphere<3> & box) const
{
  double dist;
  //  dist = sqrt (2 * CalcFunctionValue(box.Center()) * r + r * r);

  dist =  (2 * CalcFunctionValue(box.Center()) * r + r * r);
  if (dist <= 0) dist = 0;
  else dist = sqrt (dist + 1e-16);

  if (dist - box.Diam()/2 > r) return IS_OUTSIDE;
  if (dist + box.Diam()/2 < r) return IS_INSIDE;
  return DOES_INTERSECT;
}


double Cylinder :: HesseNorm () const
{
  return 2 / r;
}

Point<3> Cylinder :: GetSurfacePoint () const
{
  Vec<3> vr;
  if (fabs (vab(0)) > fabs(vab(2)))
    vr = Vec<3> (vab(1), -vab(0), 0);
  else
    vr = Vec<3> (0, -vab(2), vab(1));
    
  vr *= (r / vr.Length());
  return a + vr;
}

void Cylinder :: GetTriangleApproximation 
(TriangleApproximation & tas, 
 const Box<3> & boundingbox, double facets) const
{
  int i, j;
  double lg, bg;
  int n = int(facets) + 1;  

  Vec<3> lvab = b - a;
  Vec<3> n1 = lvab.GetNormal();
  Vec<3> n2 = Cross (lvab, n1);
  
  n1.Normalize();
  n2.Normalize();


  for (j = 0; j <= n; j++)
    for (i = 0; i <= n; i++)
      {
	lg = 2 * M_PI * double (i) / n;
	bg = double(j) / n;

	Point<3> p = a + (bg * lvab) 
	  + ((r * cos(lg)) * n1) 
	  + ((r * sin(lg)) * n2);

	tas.AddPoint (p);
      }

  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++)
      {
	int pi = i + (n+1) * j;
	tas.AddTriangle (TATriangle (0, pi, pi+1, pi+n+2));
	tas.AddTriangle (TATriangle (0, pi, pi+n+2, pi+n+1));
      }
}









EllipticCylinder :: 
EllipticCylinder (const Point<3> & aa,
		  const Vec<3> & avl, const Vec<3> & avs)
{
  a = aa;
  if(avl.Length2() > avs.Length2())
    {
      vl = avl;
      vs = avs;
    }
  else
    {
      vl = avs;
      vs = avl;
    }

  CalcData();
  // Print (cout);
}


void EllipticCylinder :: CalcData ()
{
  // f = (x-a, vl)^2 / |vl|^2 + (x-a, vs)^2 / |vs|^2 -1

  Vec<3> hvl, hvs;
  double lvl = vl.Length2 ();
  if (lvl < 1e-32) lvl = 1;
  double lvs = vs.Length2 ();
  if (lvs < 1e-32) lvs = 1;

  hvl = (1.0 / lvl) * vl;
  hvs = (1.0 / lvs) * vs;

  cxx = hvl(0) * hvl(0) + hvs(0) * hvs(0);
  cyy = hvl(1) * hvl(1) + hvs(1) * hvs(1);
  czz = hvl(2) * hvl(2) + hvs(2) * hvs(2);

  cxy = 2 * (hvl(0) * hvl(1) + hvs(0) * hvs(1));
  cxz = 2 * (hvl(0) * hvl(2) + hvs(0) * hvs(2));
  cyz = 2 * (hvl(1) * hvl(2) + hvs(1) * hvs(2));

  Vec<3> va (a);
  c1 = pow(va * hvl,2) + pow(va * hvs,2) - 1;
  
  Vec<3> v = -2 * (va * hvl) * hvl - 2 * (va * hvs) * hvs;
  cx = v(0);
  cy = v(1);
  cz = v(2);
}


INSOLID_TYPE EllipticCylinder :: BoxInSolid (const BoxSphere<3> & box) const
{
  double grad = 2.0 / vs.Length ();
  double ggrad = 1.0 / vs.Length2 ();

  double val = CalcFunctionValue (box.Center());
  double r = box.Diam() / 2;
  double maxval = grad * r + ggrad * r * r;

  // (*testout) << "box = " << box << ", val = " << val << ", maxval = " << maxval << endl;

  if (val > maxval) return IS_OUTSIDE;
  if (val < -maxval) return IS_INSIDE;
  return DOES_INTERSECT;
}


double EllipticCylinder :: HesseNorm () const
{
  return 1.0/min(vs.Length2 (),vl.Length2());
}

double EllipticCylinder :: MaxCurvature () const
{
  double aa = vs.Length();
  double bb = vl.Length();

  return max2(bb/(aa*aa),aa/(bb*bb));
}

double EllipticCylinder :: MaxCurvatureLoc (const Point<3> & c, 
                                            double rad) const
{
#ifdef JOACHIMxxx
  cout << "max curv local" << endl;
  return 0.02;
#endif

  // saubere Loesung wird noch notwendig !!!
  double aa = vs.Length();
  double bb = vl.Length();
  return max2(bb/(aa*aa),aa/(bb*bb));
}



Point<3> EllipticCylinder :: GetSurfacePoint () const
{
  return a + vl;
}



void EllipticCylinder :: GetTriangleApproximation 
(TriangleApproximation & tas, 
 const Box<3> & boundingbox, double facets) const
{
  int i, j;
  double lg, bg;
  int n = int(facets) + 1;  

  Vec<3> axis = Cross (vl, vs);

  for (j = 0; j <= n; j++)
    for (i = 0; i <= n; i++)
      {
	lg = 2 * M_PI * double (i) / n;
	bg = double(j) / n;

	Point<3> p = a + (bg * axis)
	  + cos(lg) * vl + sin(lg) * vs;

	tas.AddPoint (p);
      }

  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++)
      {
	int pi = i + (n+1) * j;
	tas.AddTriangle (TATriangle (0, pi, pi+1, pi+n+2));
	tas.AddTriangle (TATriangle (0, pi, pi+n+2, pi+n+1));
      }
}










Cone :: Cone (const Point<3> & aa, const Point<3> & ab, 
	      double ara, double arb)
{
  a = aa;
  b = ab;
  ra = ara;
  rb = arb;

  CalcData();
  // Print (cout);
}


Primitive * Cone :: CreateDefault ()
{
  return new Cone (Point<3> (0,0,0), Point<3> (1,0,0), 0.5, 0.2);
}




void Cone :: GetPrimitiveData (const char *& classname, ARRAY<double> & coeffs) const
{
  classname = "cone";
  coeffs.SetSize (8);
  coeffs.Elem(1) = a(0);
  coeffs.Elem(2) = a(1);
  coeffs.Elem(3) = a(2);
  coeffs.Elem(4) = b(0);
  coeffs.Elem(5) = b(1);
  coeffs.Elem(6) = b(2);
  coeffs.Elem(7) = ra;
  coeffs.Elem(8) = rb;
}

void Cone :: SetPrimitiveData (ARRAY<double> & coeffs)
{
  a(0) = coeffs.Elem(1);
  a(1) = coeffs.Elem(2);
  a(2) = coeffs.Elem(3);
  b(0) = coeffs.Elem(4);
  b(1) = coeffs.Elem(5);
  b(2) = coeffs.Elem(6);
  ra = coeffs.Elem(7);
  rb = coeffs.Elem(8);

  CalcData();
}

void Cone :: CalcData ()
{

  minr = (ra < rb) ? ra : rb;

  vab = b - a;
  vabl = vab.Length();

  Vec<3> va (a);

  //
  //   f = r(P)^2 - R(z(P))^2
  //
  //   z(P) = t0vec * P + t0 = (P-a, b-a)/(b-a,b-a)
  //   R(z(P)) = t1vec * P + t1 = rb * z + ra * (1-z)
  //   r(P)^2 =||P-a||^2 - ||a-b||^2 z^2k


  t0vec = vab;
  t0vec /= (vabl * vabl);
  t0 = -(va * vab) / (vabl * vabl);

  t1vec = t0vec;
  t1vec *= (rb - ra);
  t1 = ra + (rb - ra) * t0; 

  cxx = cyy = czz = 1;
  cxy = cxz = cyz = 0;

  cxx = 1 - (vab*vab) * t0vec(0) * t0vec(0) - t1vec(0) * t1vec(0);
  cyy = 1 - (vab*vab) * t0vec(1) * t0vec(1) - t1vec(1) * t1vec(1);
  czz = 1 - (vab*vab) * t0vec(2) * t0vec(2) - t1vec(2) * t1vec(2);
  
  cxy = -2 * (vab * vab) * t0vec(0) * t0vec(1) - 2 * t1vec(0) * t1vec(1);
  cxz = -2 * (vab * vab) * t0vec(0) * t0vec(2) - 2 * t1vec(0) * t1vec(2);
  cyz = -2 * (vab * vab) * t0vec(1) * t0vec(2) - 2 * t1vec(1) * t1vec(2);

  cx = -2 * a(0) - 2 * (vab * vab) * t0 * t0vec(0) - 2 * t1 * t1vec(0);
  cy = -2 * a(1) - 2 * (vab * vab) * t0 * t0vec(1) - 2 * t1 * t1vec(1);
  cz = -2 * a(2) - 2 * (vab * vab) * t0 * t0vec(2) - 2 * t1 * t1vec(2);

  c1 = va.Length2() - (vab * vab) * t0 * t0 - t1 * t1;


  double maxr = max2(ra,rb);
  cxx /= maxr; cyy /= maxr; czz /= maxr;
  cxy /= maxr; cxz /= maxr; cyz /= maxr;
  cx /= maxr; cy /= maxr; cz /= maxr;
  c1 /= maxr;


  // (*testout) << "t0vec = " << t0vec << " t0 = " << t0 << endl;
  // (*testout) << "t1vec = " << t1vec << " t1 = " << t1 << endl;
  // PrintCoeff (*testout);
}


INSOLID_TYPE Cone :: BoxInSolid (const BoxSphere<3> & box) const
{
  double rp, dist;

  Vec<3> cv(box.Center());

  rp = cv * t1vec + t1;
  dist = sqrt (CalcFunctionValue(box.Center()) *max2(ra,rb) + rp * rp) - rp;

  if (dist - box.Diam() > 0) return IS_OUTSIDE;
  if (dist + box.Diam() < 0) return IS_INSIDE;
  return DOES_INTERSECT;
}


double Cone :: HesseNorm () const
{
  return 2 / minr;
}


double Cone ::  LocH (const Point<3> & p, double x, 
				  double c, double hmax) const
{
  //double bloch = Surface::LocH (p, x, c, hmax);
  Vec<3> g;
  CalcGradient (p, g);

  double lam = Abs(g);
  double meancurv = 
    -( 2  * g(0)*g(1)*cxy - 2 * czz * (g(0)*g(0)+g(1)*g(1))
       +2 * g(1)*g(2)*cyz - 2 * cxx * (g(1)*g(1)+g(2)*g(2))
       +2 * g(0)*g(2)*cxz - 2 * cyy * (g(0)*g(0)+g(2)*g(2))) / (3*lam*lam*lam);

  // cout << "type = " << typeid(*this).name() << ", baseh = " << bloch << ", meancurv = " << meancurv << endl;
  // return bloch;
  
  meancurv = fabs (meancurv);
  if (meancurv < 1e-20) meancurv = 1e-20;

  // cout << "c = " << c << ", safety = " << mparam.curvaturesafety << endl;
  double hcurv = 1.0/(4*meancurv*mparam.curvaturesafety);

  return min2 (hmax, hcurv);
}


Point<3> Cone :: GetSurfacePoint () const
{
  Vec<3> vr = vab.GetNormal ();
  
  vr *= (ra / vr.Length());
  return a + vr;
}





void Cone :: GetTriangleApproximation 
(TriangleApproximation & tas, 
 const Box<3> & boundingbox, double facets) const
{
  int i, j;
  double lg, bg;
  int n = int(facets) + 1;  

  Vec<3> lvab = b - a;
  Vec<3> n1 = lvab.GetNormal();
  Vec<3> n2 = Cross (lvab, n1);
  
  n1.Normalize();
  n2.Normalize();


  for (j = 0; j <= n; j++)
    for (i = 0; i <= n; i++)
      {
	lg = 2 * M_PI * double (i) / n;
	bg = double(j) / n;

	Point<3> p = a + (bg * lvab) 
	  + (( (ra+(rb-ra)*bg)  * cos(lg)) * n1) 
	  + (( (ra+(rb-ra)*bg)  * sin(lg)) * n2);

	tas.AddPoint (p);
      }

  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++)
      {
	int pi = i + (n+1) * j;
	tas.AddTriangle (TATriangle (0, pi, pi+1, pi+n+2));
	tas.AddTriangle (TATriangle (0, pi, pi+n+2, pi+n+1));
      }
}









 /// Torus 
 /// Lorenzo Codecasa (codecasa@elet.polimi.it)
 /// April 27th, 2005 
 ///
 /// begin...
Torus :: Torus (const Point<3> & ac, const Vec<3> & an, double aR, double ar)
{
  c = ac;
  n = an;
  R = aR;
  r = ar;
}

void Torus :: GetPrimitiveData (const char *& classname, ARRAY<double> & coeffs) const
{
  classname = "torus";
  coeffs.SetSize (8);
  coeffs.Elem(1) = c(0);
  coeffs.Elem(2) = c(1);
  coeffs.Elem(3) = c(2);
  coeffs.Elem(4) = n(0);
  coeffs.Elem(5) = n(1);
  coeffs.Elem(6) = n(2);
  coeffs.Elem(7) = R;
  coeffs.Elem(8) = r;
}

void Torus :: SetPrimitiveData (ARRAY<double> & coeffs)
{
  c(0) = coeffs.Elem(1);
  c(1) = coeffs.Elem(2);
  c(2) = coeffs.Elem(3);
  n(0) = coeffs.Elem(4);
  n(1) = coeffs.Elem(5);
  n(2) = coeffs.Elem(6);
  R = coeffs.Elem(7);
  r = coeffs.Elem(8);
}

Primitive * Torus :: CreateDefault ()
{
  return new Torus (Point<3> (0,0,0), Vec<3> (0,0,1), 2, 1);
}

Primitive * Torus :: Copy () const
{
  return new Torus (c, n, R, r);
}

void Torus :: Transform (Transformation<3> & trans)
{
  Point<3> hc;
  trans.Transform (c, hc);
  c = hc;
  
  Vec<3> hn;
  trans.Transform (n, hn);
  n = hn;
}

int Torus :: IsIdentic (const Surface & s2, int & inv, double eps) const
{
  const Torus * torus2 = dynamic_cast<const Torus*>  (&s2);

  if (!torus2) return 0;

  if (fabs (torus2->R - R) > eps) return 0;
  
  if (fabs (torus2->r - r) > eps) return 0;

  Vec<3> v2 = torus2->n - n;
  if ( v2 * v2 > eps ) return 0;
  
  v2 = torus2->c - c;
  if ( v2 * v2 > eps ) return 0;

  inv = 0;
  return 1;
}

double Torus :: CalcFunctionValue (const Point<3> & point) const
{
	Vec<3> v1 = point - c;
	double a1 = v1(0) * v1(0) + v1(1) * v1(1) + v1(2) * v1(2);
	double a2 = n(0) * v1(0) + n(1) * v1(1) + n(2) * v1(2);
	double a3 = a1 + R * R - r * r;
	double a4 = n(0) * n(0) + n(1) * n(1) + n(2) * n(2);
	return ( a3 * a3 -4 * R * R * ( a1 - a2 * a2 / a4 ) ) / ( R * R * R );
}

void Torus :: CalcGradient (const Point<3> & point, Vec<3> & grad) const
{
	Vec<3> v1 = point - c;
	double a1 = v1(0) * v1(0) + v1(1) * v1(1) + v1(2) * v1(2);
	double a2 = n(0) * v1(0) + n(1) * v1(1) + n(2) * v1(2);
	double a3 = a1 - R * R - r * r;
	double a4 = n(0) * n(0) + n(1) * n(1) + n(2) * n(2);
	grad(0) = ( 4 * a3 * v1(0) + 8 * R * R * a2 / a4 * n(0) ) / ( R * R * R );
	grad(1) = ( 4 * a3 * v1(1) + 8 * R * R * a2 / a4 * n(1) ) / ( R * R * R );
	grad(2) = ( 4 * a3 * v1(2) + 8 * R * R * a2 / a4 * n(2) ) / ( R * R * R );
}

void Torus :: CalcHesse (const Point<3> & point, Mat<3> & hesse) const
{
	Vec<3> v1 = point - c;
	double a1 = v1(0) * v1(0) + v1(1) * v1(1) + v1(2) * v1(2);
	double a3 = a1 - R * R - r * r;
	double a4 = n(0) * n(0) + n(1) * n(1) + n(2) * n(2);
	hesse(0,0) = ( 4 * a3 + 8 * (v1(0) * v1(0) + (R * n(0)) * (R * n(0)) / a4 ) ) / ( R * R * R );
	hesse(1,1) = ( 4 * a3 + 8 * (v1(1) * v1(1) + (R * n(1)) * (R * n(1)) / a4 ) ) / ( R * R * R );
	hesse(2,2) = ( 4 * a3 + 8 * (v1(2) * v1(2) + (R * n(2)) * (R * n(2)) / a4 ) ) / ( R * R * R );
	hesse(0,1) = hesse(1,0) = 8 * (v1(0) * v1(1) + (R * n(0)) * (R * n(1)) / a4 ) / ( R * R * R );
	hesse(1,2) = hesse(2,1) = 8 * (v1(1) * v1(2) + (R * n(1)) * (R * n(2)) / a4) / ( R * R * R );
	hesse(0,2) = hesse(2,0) = 8 * (v1(0) * v1(2) + (R * n(0)) * (R * n(2)) / a4) / ( R * R * R );
}

double Torus :: HesseNorm () const
{	
	return  ( 2 / r + 2 / ( R - r ) );
}

Point<3> Torus :: GetSurfacePoint () const
{
  Vec<3> vn = n.GetNormal();
  return c + ( R + r ) * vn.Normalize();
}

/// void Torus :: DefineTangentialPlane (const Point<3> & ap1, const Point<3> & ap2)
/// {
/// }

/// void Torus :: ToPlane (const Point<3> & p, 
///			  Point<2> & pplane, 
///			  double h, int & zone) const
/// {
/// }

/// void Torus :: FromPlane (const Point<2> & pplane, Point<3> & p, double h) const
/// {
/// }

/// void Torus :: Project (Point<3> & p) const
/// {
/// }

INSOLID_TYPE Torus :: BoxInSolid (const BoxSphere<3> & box) const
{
  Vec<3> v1 = box.Center() - c;
  double a1 = v1(0) * v1(0) + v1(1) * v1(1) + v1(2) * v1(2);
  double a2 = n(0) * v1(0) + n(1) * v1(1) + n(2) * v1(2);
  double a4 = n(0) * n(0) + n(1) * n(1) + n(2) * n(2);
 
  double dist = sqrt( a1 + R * R - 2 * R * sqrt( a1 - a2 * a2 / a4) );

  if (dist - box.Diam()/2 > r) return IS_OUTSIDE;
  if (dist + box.Diam()/2 < r) return IS_INSIDE;
  return DOES_INTERSECT;
}

void Torus :: GetTriangleApproximation (TriangleApproximation & tas, 
                                        const Box<3> & boundingbox, double facets) const
{
	int i, j;
	double lg, bg;
	int N = int(facets) + 1;  
	
	Vec<3> lvab = n ;
	lvab.Normalize();
	
	Vec<3> n1 = lvab.GetNormal();
	n1.Normalize();
	
	Vec<3> n2 = Cross(lvab, n1);
	n2.Normalize();
	
	for (j = 0; j <= N; j++)
	for (i = 0; i <= N; i++)
	{
		lg = 2 * M_PI * double (i) / N;
		bg = 2 * M_PI * double(j) / N;
	
		Point<3> p = c + ( R + r * cos(lg) ) * ( cos(bg) * n1 + sin(bg) * n2 ) + r * sin(lg) * n;
		tas.AddPoint (p);
	}
	
	for (j = 0; j < N; j++)
	for (i = 0; i < N; i++)
	{
		int pi = i + (N+1) * j;
		tas.AddTriangle (TATriangle (0, pi, pi+1, pi+N+2));
		tas.AddTriangle (TATriangle (0, pi, pi+N+2, pi+N+1));
	}
} 
  
void Torus :: Read (istream & ist)
{
  ist >> c(0) >> c(1) >> c(2) >> n(0) >> n(1) >> n(2) >> R >> r;
}

void Torus :: Print (ostream & ost) const
{
  ost << c(0) << "  " << c(1) << "  " << c(2) << "  "
      << n(0) << "  " << n(1) << "  " << n(2) << "  "
      << R    << "  " << r    << endl;
}

/// end...
















}
