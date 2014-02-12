#include <mystdlib.h>

#include <myadt.hpp>
#include <csg.hpp>

#include <linalg.hpp>
#include <meshing.hpp>


namespace netgen
{
Surface :: Surface ()
{
  maxh = 1e10;
  name = new char[7];
  strcpy (name, "noname");
  bcprop = -1;
  bcname = "default";
}

Surface :: ~Surface()
{
  delete [] name;
}


void Surface :: SetName (const char * aname)
{
  delete [] name;
  name = new char[strlen (aname)+1];
  strcpy (name, aname);
}


int Surface :: PointOnSurface (const Point<3> & p,
			       double eps) const
{
  double val = CalcFunctionValue (p);
  return fabs (val) < eps;
}


void Surface :: CalcHesse (const Point<3> & point, Mat<3> & hesse) const
{
  double dx = 1e-5;
  Point<3> hp1, hp2;
  Vec<3> g1, g2;

  for (int i = 0; i < 3; i++)
    {
      hp1 = point;
      hp2 = point;

      hp1(i) += dx;
      hp2(i) -= dx;

      CalcGradient (hp1, g1);
      CalcGradient (hp2, g2);
      	
      for (int j = 0; j < 3; j++)
	hesse(i, j) = (g1(j) - g2(j)) / (2 * dx);
    }
}
  
/*
void Surface :: GetNormalVector (const Point<3> & p, Vec<3> & n) const
{
  CalcGradient (p, n);
  n.Normalize();
}
*/
Vec<3> Surface :: GetNormalVector (const Point<3> & p) const
{
  Vec<3> n;
  CalcGradient (p, n);
  n.Normalize();
  return n;
}

void Surface :: DefineTangentialPlane (const Point<3> & ap1, 
				       const Point<3> & ap2)
{
  p1 = ap1;
  p2 = ap2;
  
  ez = GetNormalVector (p1);
  ex = p2 - p1;
  ex -= (ex * ez) * ez;
  ex.Normalize();
  ey = Cross (ez, ex);  
}

void Surface :: ToPlane (const Point<3> & p3d, Point<2> & pplane, 
			 double h, int & zone) const
{
  Vec<3> p1p, n;

  n = GetNormalVector (p3d);
  if (n * ez < 0)
    {
      zone = -1;
      pplane(0) = 1e8;
      pplane(1) = 1e9;
      return;
    }
  
  p1p = p3d - p1;
  pplane(0) = (p1p * ex) / h;
  pplane(1) = (p1p * ey) / h;
  zone = 0;
}	

void Surface :: FromPlane (const Point<2> & pplane, 
			   Point<3> & p3d, double h) const 
{ 
  p3d = p1 
    + (h * pplane(0)) * ex 
    + (h * pplane(1)) * ey;
  
  Project (p3d);
}

void Surface :: Project (Point<3> & p) const
{
  Vec<3> n;
  double val;

  for (int i = 1; i <= 10; i++)
    {
      val = CalcFunctionValue (p);
      if (fabs (val) < 1e-12) return;
	
      CalcGradient (p, n);
      p -= (val / Abs2 (n)) * n;
    }
}

void Surface :: SkewProject (Point<3> & p, const Vec<3> & direction) const
{
  Point<3> startp(p);
  double t_old(0),t_new(1);
  Vec<3> grad;
  for(int i=0; fabs(t_old-t_new) > 1e-20 && i<15; i++)
    {
      t_old = t_new;
      CalcGradient(p,grad);
      t_new = t_old - CalcFunctionValue(p)/(grad*direction);
      p = startp + t_new*direction;
    }
}


double Surface :: MaxCurvature () const
{ 
  return 0.5 * HesseNorm (); 
}

double Surface :: 
MaxCurvatureLoc (const Point<3> & /* c */ , double /* rad */) const
{ 
  return MaxCurvature (); 
}
              


double Surface :: LocH (const Point<3> & p, double x, 
			double c, double hmax) const
  // finds h <= hmax, s.t.  h * \kappa_x*h < c
{
  /*
    double h, hmin, kappa;
    hmin = 0;
  
    while (hmin < 0.9 * hmax)
    {
    h = 0.5 * (hmin + hmax);
    kappa = 2 * MaxCurvatureLoc (p, x * h);
      
    if (kappa * h >= c)
    hmax = h;
    else
    hmin = h;
    }
    return h;
  */

  double hret;
  double kappa = MaxCurvatureLoc (p, x*hmax);

  kappa *= c *  mparam.curvaturesafety;
  
  if (hmax * kappa < 1)
    hret = hmax;
  else
    hret = 1 / kappa;

  if (maxh < hret)
    hret = maxh;

  return hret;
}




Primitive :: Primitive ()
{
  surfaceids.SetSize (1);
  surfaceactive.SetSize (1);
  surfaceactive[0] = 1;
}

Primitive :: ~Primitive()
{
  ;
}

int Primitive :: GetSurfaceId (int i) const
{
  return surfaceids[i];
}

void Primitive :: SetSurfaceId (int i, int id) 
{
  surfaceids[i] = id;
}




void Primitive :: GetPrimitiveData (const char *& classname, 
				    ARRAY<double> & coeffs) const
{
  classname = "undef";
  coeffs.SetSize (0);
}

void Primitive :: SetPrimitiveData (ARRAY<double> & coeffs)
{
  ;
}

Primitive * Primitive :: CreatePrimitive (const char * classname)
{
  if (strcmp (classname, "sphere") == 0)
    return Sphere::CreateDefault();
  if (strcmp (classname, "plane") == 0)
    return Plane::CreateDefault();
  if (strcmp (classname, "cylinder") == 0)
    return Cylinder::CreateDefault();
  if (strcmp (classname, "cone") == 0)
    return Cone::CreateDefault();
  if (strcmp (classname, "brick") == 0)
    return Brick::CreateDefault();


  stringstream ost;
  ost << "Primitve::CreatePrimitive not implemented for " << classname << endl;
  throw NgException (ost.str());
}


Primitive * Primitive :: Copy () const
{
  stringstream ost;
  ost << "Primitve::Copy not implemented for " << typeid(*this).name() << endl;
  throw NgException (ost.str());
}


void Primitive :: Transform (Transformation<3> & trans)
{
  stringstream ost;
  ost << "Primitve::Transform not implemented for " << typeid(*this).name() << endl;
  throw NgException (ost.str());
}

void Primitive :: GetTangentialSurfaceIndices (const Point<3> & p, 
					       ARRAY<int> & surfind, double eps) const
{
  for (int j = 0; j < GetNSurfaces(); j++)
    if (fabs (GetSurface(j).CalcFunctionValue (p)) < eps)
      if (!surfind.Contains (GetSurfaceId(j)))
	surfind.Append (GetSurfaceId(j));
}


void Primitive :: 
GetTangentialVecSurfaceIndices (const Point<3> & p, const Vec<3> & v,
				ARRAY<int> & surfind, double eps) const
{
  cout << "get tangvecsurfind not implemented" << endl;
  surfind.SetSize (0);
}

void Primitive :: 
GetTangentialVecSurfaceIndices2 (const Point<3> & p, const Vec<3> & v1, const Vec<3> & v2,
				 ARRAY<int> & surfind, double eps) const
{
  for (int j = 0; j < GetNSurfaces(); j++)
    {
      if (fabs (GetSurface(j).CalcFunctionValue (p)) < eps)
	{
	  Vec<3> grad;
	  GetSurface(j).CalcGradient (p, grad);
	  if (sqr (grad * v1) < 1e-6 * v1.Length2() * grad.Length2()  && 
	      sqr (grad * v2) < 1e-6 * v2.Length2() * grad.Length2() )   // new, 18032006 JS
	    {
	      if (!surfind.Contains (GetSurfaceId(j)))
		surfind.Append (GetSurfaceId(j));
	    }
	}
    }
}




INSOLID_TYPE Primitive :: 
VecInSolid2 (const Point<3> & p,
	     const Vec<3> & v1,
	     const Vec<3> & v2,
	     double eps) const
{
  //(*testout) << "Primitive::VecInSolid2" << endl;
  Point<3> hp = p + 1e-3 * v1 + 1e-5 * v2;

  INSOLID_TYPE res = PointInSolid (hp, eps);
  //  (*testout) << "vectorin2, type = " << typeid(*this).name() << ", res = " << res << endl;

  return res;
}

INSOLID_TYPE Primitive :: 
VecInSolid3 (const Point<3> & p,
	     const Vec<3> & v1,
	     const Vec<3> & v2,
	     double eps) const
{
  //(*testout) << "Primitive::VecInSolid3" << endl;
  return VecInSolid (p, v1, eps);
}

INSOLID_TYPE Primitive :: 
VecInSolid4 (const Point<3> & p,
	     const Vec<3> & v,
	     const Vec<3> & v2,
	     const Vec<3> & m,
	     double eps) const
{
  return VecInSolid2 (p, v, m, eps);
}





OneSurfacePrimitive :: OneSurfacePrimitive()
{
  ;
}

OneSurfacePrimitive :: ~OneSurfacePrimitive()
{
  ;
}


INSOLID_TYPE OneSurfacePrimitive :: 
PointInSolid (const Point<3> & p,
	      double eps) const
{
  double hv1 = (GetSurface(0).CalcFunctionValue(p));
  if (hv1 <= -eps)
    return IS_INSIDE;
  if (hv1 >= eps)
    return IS_OUTSIDE;
  return DOES_INTERSECT;
}
 

INSOLID_TYPE OneSurfacePrimitive :: 
VecInSolid (const Point<3> & p, const Vec<3> & v,
	    double eps) const
{
  double hv1 = (GetSurface(0).CalcFunctionValue(p));
  if (hv1 <= -eps)
    return IS_INSIDE;
  if (hv1 >= eps)
    return IS_OUTSIDE;


  Vec<3> hv;
  GetSurface(0).CalcGradient (p, hv);

  hv1 = v * hv;

  if (hv1 <= -eps)
    return IS_INSIDE;
  if (hv1 >= eps)
    return IS_OUTSIDE;

  return DOES_INTERSECT;
}


INSOLID_TYPE OneSurfacePrimitive :: 
VecInSolid2 (const Point<3> & p,
	     const Vec<3> & v1,
	     const Vec<3> & v2,
	     double eps) const
{
  double hv1 = (GetSurface(0).CalcFunctionValue(p));
  if (hv1 <= -eps)
    return IS_INSIDE;
  if (hv1 >= eps)
    return IS_OUTSIDE;

  Vec<3> hv;

  GetSurface(0).CalcGradient (p, hv);

  hv1 = v1 * hv;
  if (hv1 <= -eps)
    return IS_INSIDE;
  if (hv1 >= eps)
    return IS_OUTSIDE;

  double hv2 = v2 * hv;
  if (hv2 <= 0)
    return IS_INSIDE;
  else
    return IS_OUTSIDE;
}
  


INSOLID_TYPE OneSurfacePrimitive :: 
VecInSolid3 (const Point<3> & p, const Vec<3> & v, const Vec<3> & v2,
	     double eps) const
{
  //(*testout) << "OneSurfacePrimitive::VecInSolid3" << endl;
  double hv1 = (GetSurface(0).CalcFunctionValue(p));
  if (hv1 <= -eps)
    return IS_INSIDE;
  if (hv1 >= eps)
    return IS_OUTSIDE;

  Vec<3> grad;
  GetSurface(0).CalcGradient (p, grad);

  hv1 = v * grad;
  if (hv1 <= -eps) return IS_INSIDE;
  if (hv1 >= eps) return IS_OUTSIDE;

  Mat<3> hesse;
  GetSurface(0).CalcHesse (p, hesse);

  double hv2 = v2 * grad + v * (hesse * v);

  if (hv2 <= -eps) return IS_INSIDE;
  if (hv2 >= eps) return IS_OUTSIDE;

  return DOES_INTERSECT;
}




INSOLID_TYPE OneSurfacePrimitive :: 
VecInSolid4 (const Point<3> & p, const Vec<3> & v, const Vec<3> & v2,
	     const Vec<3> & m,
	     double eps) const
{
  double hv1 = (GetSurface(0).CalcFunctionValue(p));
  if (hv1 <= -eps)
    return IS_INSIDE;
  if (hv1 >= eps)
    return IS_OUTSIDE;

  Vec<3> grad;
  GetSurface(0).CalcGradient (p, grad);

  hv1 = v * grad;
  if (hv1 <= -eps) return IS_INSIDE;
  if (hv1 >= eps) return IS_OUTSIDE;

  Mat<3> hesse;
  GetSurface(0).CalcHesse (p, hesse);

  double hv2 = v2 * grad + v * (hesse * v);

  if (hv2 <= -eps) return IS_INSIDE;
  if (hv2 >= eps) return IS_OUTSIDE;


  double hv3 = m * grad;
  if (hv3 <= -eps) return IS_INSIDE;
  if (hv3 >= eps) return IS_OUTSIDE;

  return DOES_INTERSECT;
}







int OneSurfacePrimitive :: GetNSurfaces() const
{
  return 1;
}

Surface & OneSurfacePrimitive :: GetSurface (int i)
{
  return *this;
}

const Surface & OneSurfacePrimitive :: GetSurface (int i) const
{
  return *this;
}






void ProjectToEdge (const Surface * f1, const Surface * f2, Point<3> & hp)
{
  Vec<2> rs, lam;
  Vec<3> a1, a2;
  Mat<2> a;

  int i = 10;
  while (i > 0)
    {
      i--;
      rs(0) = f1 -> CalcFunctionValue (hp);
      rs(1) = f2 -> CalcFunctionValue (hp);
      f1->CalcGradient (hp, a1);
      f2->CalcGradient (hp, a2);

      double alpha = fabs(a1*a2)/sqrt(a1.Length2()*a2.Length2());
      if(fabs(1.-alpha) < 1e-6)
	{
	  if(fabs(rs(0)) >= fabs(rs(1)))
	    f1 -> Project(hp);
	  else
	    f2 -> Project(hp);
	}
      else
	{

	  a(0,0) = a1 * a1;
	  a(0,1) = a(1,0) = a1 * a2;
	  a(1,1) = a2 * a2;
	  
	  a.Solve (rs, lam);

	  hp -= lam(0) * a1 + lam(1) * a2;
	}

      if (Abs2 (rs) < 1e-24 && i > 1) i = 1;
    }
}
}
