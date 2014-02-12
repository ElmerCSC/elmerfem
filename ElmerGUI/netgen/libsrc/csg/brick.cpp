#include <mystdlib.h>

#include <linalg.hpp>
#include <csg.hpp>

namespace netgen
{

Parallelogram3d :: Parallelogram3d (Point<3> ap1, Point<3> ap2, Point<3> ap3)
{
  p1 = ap1;
  p2 = ap2;
  p3 = ap3;

  CalcData();
}

Parallelogram3d ::~Parallelogram3d ()
{
  ;
}

void Parallelogram3d :: SetPoints (Point<3> ap1, 
				   Point<3> ap2, 
				   Point<3> ap3)
{
  p1 = ap1;
  p2 = ap2;
  p3 = ap3;

  CalcData();
}

void Parallelogram3d :: CalcData()
{
  v12 = p2 - p1;
  v13 = p3 - p1;
  p4 = p2 + v13;

  n = Cross (v12, v13);
  n.Normalize();
}

int Parallelogram3d :: 
IsIdentic (const Surface & s2, int & inv, double eps) const
{
  int id = 
    (fabs (s2.CalcFunctionValue (p1)) <= eps) &&
    (fabs (s2.CalcFunctionValue (p2)) <= eps) &&
    (fabs (s2.CalcFunctionValue (p3)) <= eps);

  if (id)
    {
      Vec<3> n2;
      n2 = s2.GetNormalVector(p1);
      inv = (n * n2) < 0;
    }
  return id;
}


double Parallelogram3d :: CalcFunctionValue (const Point<3> & point) const
{
  return n * (point - p1);
}

void Parallelogram3d :: CalcGradient (const Point<3> & point, 
				      Vec<3> & grad) const
{
  grad = n;
}

void Parallelogram3d :: CalcHesse (const Point<3> & point, Mat<3> & hesse) const
{
  hesse = 0;
}

double Parallelogram3d :: HesseNorm () const
{
  return 0;
}

Point<3> Parallelogram3d :: GetSurfacePoint () const
{
  return p1;
}

void Parallelogram3d :: Print (ostream & str) const
{
  str << "Parallelogram3d " << p1 << " - " << p2 << " - " << p3 << endl;
}

  
void Parallelogram3d :: 
GetTriangleApproximation (TriangleApproximation & tas, 
			  const Box<3> & bbox, 
			  double facets) const
{
  tas.AddPoint (p1);
  tas.AddPoint (p2);
  tas.AddPoint (p3);
  tas.AddPoint (p4);
  tas.AddTriangle (TATriangle (0, 0, 1, 2));
  tas.AddTriangle (TATriangle (0, 2, 1, 3));
}










Brick :: Brick (Point<3> ap1, Point<3> ap2, 
		Point<3> ap3, Point<3> ap4)
{
  faces.SetSize (6);
  surfaceids.SetSize (6);
  surfaceactive.SetSize(6);

  p1 = ap1; p2 = ap2;
  p3 = ap3; p4 = ap4;

  for (int i = 0; i < 6; i++)
    {
      faces[i] = new Plane (Point<3>(0,0,0), Vec<3> (0,0,1));
      surfaceactive[i] = 1;
    }

  CalcData();
}

Brick :: ~Brick ()
{
  for (int i = 0; i < 6; i++)
    delete faces[i];
}

Primitive * Brick :: CreateDefault ()
{
  return new Brick (Point<3> (0,0,0),
		    Point<3> (1,0,0),
		    Point<3> (0,1,0),
		    Point<3> (0,0,1));
}



Primitive * Brick :: Copy () const
{
  return new Brick (p1, p2, p3, p4);
}

void  Brick :: Transform (Transformation<3> & trans)
{
  trans.Transform (p1);
  trans.Transform (p2);
  trans.Transform (p3);
  trans.Transform (p4);

  CalcData();
}









INSOLID_TYPE Brick :: BoxInSolid (const BoxSphere<3> & box) const
{
  /*
  int i;
  double maxval;
  for (i = 1; i <= 6; i++)
    {
      double val = faces.Get(i)->CalcFunctionValue (box.Center());
      if (i == 1 || val > maxval)
	maxval = val;
    }
  
  if (maxval > box.Diam()) return IS_OUTSIDE;
  if (maxval < -box.Diam()) return IS_INSIDE;
  return DOES_INTERSECT;
  */

  bool inside = 1;
  bool outside = 0;

  Point<3> p[8];
  for (int j = 0; j < 8; j++)
    p[j] = box.GetPointNr(j);

  for (int i = 0; i < 6; i++)
    {
      bool outsidei = 1;
      for (int j = 0; j < 8; j++)
	{
	  // Point<3> p = box.GetPointNr (j);
	  double val = faces[i]->Plane::CalcFunctionValue (p[j]);

	  if (val > 0)  inside = 0;
	  if (val < 0)  outsidei = 0;
	}
      if (outsidei) outside = 1;
    }

  if (outside) return IS_OUTSIDE;
  if (inside) return IS_INSIDE;
  return DOES_INTERSECT;
}

INSOLID_TYPE Brick :: PointInSolid (const Point<3> & p,
			   double eps) const
{
  double maxval = faces[0] -> Plane::CalcFunctionValue (p);
  for (int i = 1; i < 6; i++)
    {
      double val = faces[i] -> Plane::CalcFunctionValue (p);
      if (val > maxval) maxval = val;
    }

  if (maxval > eps) return IS_OUTSIDE;
  if (maxval < -eps) return IS_INSIDE;
  return DOES_INTERSECT;
}


INSOLID_TYPE Brick :: VecInSolid (const Point<3> & p,
				  const Vec<3> & v,
				  double eps) const
{
  INSOLID_TYPE result = IS_INSIDE;
  for (int i = 0; i < faces.Size(); i++)
    {
      INSOLID_TYPE hres = faces[i]->VecInSolid(p, v, eps);
      if (hres == IS_OUTSIDE || result == IS_OUTSIDE) result = IS_OUTSIDE;
      else if (hres == DOES_INTERSECT || result == DOES_INTERSECT) result = DOES_INTERSECT;
      else result = IS_INSIDE;
    }
  return result;

  /*
  INSOLID_TYPE is = IS_INSIDE;
  Vec<3> grad;
  double scal;

  for (int i = 0; i < faces.Size(); i++)
    {
      if (faces[i] -> PointOnSurface (p, eps))
	{
	  GetSurface(i).CalcGradient (p, grad);
	  scal = v * grad;
	  
	  if (scal >= eps) 
	    is = IS_OUTSIDE;
	  if (scal >= -eps && is == IS_INSIDE)
	    is = DOES_INTERSECT;
	}
    }
  return is;
  */

  /*
  Point<3> p2 = p + 1e-2 * v;
  return PointInSolid (p2, eps);
  */
}





INSOLID_TYPE Brick :: VecInSolid2 (const Point<3> & p,
				    const Vec<3> & v1,
				    const Vec<3> & v2,
				    double eps) const
{
  INSOLID_TYPE result = IS_INSIDE;
  for (int i = 0; i < faces.Size(); i++)
    {
      INSOLID_TYPE hres = faces[i]->VecInSolid2(p, v1, v2, eps);
      if (hres == IS_OUTSIDE || result == IS_OUTSIDE) result = IS_OUTSIDE;
      else if (hres == DOES_INTERSECT || result == DOES_INTERSECT) result = DOES_INTERSECT;
      else result = IS_INSIDE;
    }
  return result;
}

INSOLID_TYPE Brick :: VecInSolid3 (const Point<3> & p,
				    const Vec<3> & v1,
				    const Vec<3> & v2,
				    double eps) const
{
  INSOLID_TYPE result = IS_INSIDE;
  for (int i = 0; i < faces.Size(); i++)
    {
      INSOLID_TYPE hres = faces[i]->VecInSolid3(p, v1, v2, eps);
      if (hres == IS_OUTSIDE || result == IS_OUTSIDE) result = IS_OUTSIDE;
      else if (hres == DOES_INTERSECT || result == DOES_INTERSECT) result = DOES_INTERSECT;
      else result = IS_INSIDE;
    }
  return result;
}

INSOLID_TYPE Brick :: VecInSolid4 (const Point<3> & p,
				    const Vec<3> & v,
				    const Vec<3> & v2,
				    const Vec<3> & m,
				    double eps) const
{
  INSOLID_TYPE result = IS_INSIDE;
  for (int i = 0; i < faces.Size(); i++)
    {
      INSOLID_TYPE hres = faces[i]->VecInSolid4(p, v, v2, m, eps);
      if (hres == IS_OUTSIDE || result == IS_OUTSIDE) result = IS_OUTSIDE;
      else if (hres == DOES_INTERSECT || result == DOES_INTERSECT) result = DOES_INTERSECT;
      else result = IS_INSIDE;
    }
  return result;
}



















void Brick :: 
GetPrimitiveData (const char *& classname, ARRAY<double> & coeffs) const
{
  classname = "brick";
  coeffs.SetSize(12);
  coeffs.Elem(1) = p1(0);
  coeffs.Elem(2) = p1(1);
  coeffs.Elem(3) = p1(2);

  coeffs.Elem(4) = p2(0);
  coeffs.Elem(5) = p2(1);
  coeffs.Elem(6) = p2(2);

  coeffs.Elem(7) = p3(0);
  coeffs.Elem(8) = p3(1);
  coeffs.Elem(9) = p3(2);

  coeffs.Elem(10) = p4(0);
  coeffs.Elem(11) = p4(1);
  coeffs.Elem(12) = p4(2);
}

void Brick :: SetPrimitiveData (ARRAY<double> & coeffs)
{
  p1(0) = coeffs.Elem(1);
  p1(1) = coeffs.Elem(2);
  p1(2) = coeffs.Elem(3);

  p2(0) = coeffs.Elem(4);
  p2(1) = coeffs.Elem(5);
  p2(2) = coeffs.Elem(6);

  p3(0) = coeffs.Elem(7);
  p3(1) = coeffs.Elem(8);
  p3(2) = coeffs.Elem(9);

  p4(0) = coeffs.Elem(10);
  p4(1) = coeffs.Elem(11);
  p4(2) = coeffs.Elem(12);

  CalcData();
}



void Brick :: CalcData()
{
  v12 = p2 - p1;
  v13 = p3 - p1;
  v14 = p4 - p1;

  Point<3> pi[8];
  int i1, i2, i3;
  int i, j;
  
  i = 0;
  for (i3 = 0; i3 <= 1; i3++)
    for (i2 = 0; i2 <= 1; i2++)
      for (i1 = 0; i1 <= 1; i1++)
	{
	  pi[i] = p1 + i1 * v12 + i2 * v13 + i3 * v14;
	  i++;
	}

  static int lface[6][4] =
  { { 1, 3, 2, 4 },
    { 5, 6, 7, 8 },
    { 1, 2, 5, 6 },
    { 3, 7, 4, 8 },
    { 1, 5, 3, 7 },
    { 2, 4, 6, 8 } };
  
  ARRAY<double> data(6);
  for (i = 0; i < 6; i++)
    {
      const Point<3> lp1 = pi[lface[i][0]-1];
      const Point<3> lp2 = pi[lface[i][1]-1];
      const Point<3> lp3 = pi[lface[i][2]-1];

      Vec<3> n = Cross ((lp2-lp1), (lp3-lp1));
      n.Normalize();
      
      for (j = 0; j < 3; j++)
	{
	  data[j] = lp1(j);
	  data[j+3] = n(j);
	}
      faces[i] -> SetPrimitiveData (data);
      /* 
	 {
	 faces.Elem(i+1) -> SetPoints
	 (pi[lface[i][0]-1],
	 pi[lface[i][1]-1],
	 pi[lface[i][2]-1]);
	 }
      */
    }
}


void Brick :: Reduce (const BoxSphere<3> & box)
{
  double val;
  // Point<3> p;
  Point<3> p[8];
  for(int j=0;j<8;j++)
    p[j]=box.GetPointNr(j);

  for (int i = 0; i < 6; i++)
    {
      bool hasout = 0;
      bool hasin = 0;
      for (int j = 0; j < 8; j++)
	{
	  // p = box.GetPointNr (j);
	  val = faces[i]->Plane::CalcFunctionValue (p[j]);
	  if (val > 0)  hasout = 1;
	  else if (val < 0)  hasin = 1;
	  if (hasout && hasin) break;
	}
      surfaceactive[i] =  hasout && hasin;
    }
}

void Brick :: UnReduce ()
{ 
  for (int i = 0; i < 6; i++)
    surfaceactive[i] = 1;
}



OrthoBrick :: OrthoBrick (const Point<3> & ap1, const Point<3> & ap2)
  : Brick (ap1, 
	   Point<3> (ap2(0), ap1(1), ap1(2)),
	   Point<3> (ap1(0), ap2(1), ap1(2)),
	   Point<3> (ap1(0), ap1(1), ap2(2)))
{
  pmin = ap1;
  pmax = ap2;
}
	 
INSOLID_TYPE OrthoBrick :: BoxInSolid (const BoxSphere<3> & box) const
{
  if (pmin(0) > box.PMax()(0) ||
      pmin(1) > box.PMax()(1) ||
      pmin(2) > box.PMax()(2) ||
      pmax(0) < box.PMin()(0) ||
      pmax(1) < box.PMin()(1) ||
      pmax(2) < box.PMin()(2))
    return IS_OUTSIDE;

  if (pmin(0) < box.PMin()(0) &&
      pmin(1) < box.PMin()(1) &&
      pmin(2) < box.PMin()(2) &&
      pmax(0) > box.PMax()(0) &&
      pmax(1) > box.PMax()(1) &&
      pmax(2) > box.PMax()(2))
    return IS_INSIDE;

  return DOES_INTERSECT;
}


void OrthoBrick :: Reduce (const BoxSphere<3> & box)
{
  surfaceactive.Elem(1) =
    (box.PMin()(2) < pmin(2)) && (pmin(2) < box.PMax()(2));
  surfaceactive.Elem(2) =
    (box.PMin()(2) < pmax(2)) && (pmax(2) < box.PMax()(2));

  surfaceactive.Elem(3) =
    (box.PMin()(1) < pmin(1)) && (pmin(1) < box.PMax()(1));
  surfaceactive.Elem(4) =
    (box.PMin()(1) < pmax(1)) && (pmax(1) < box.PMax()(1));

  surfaceactive.Elem(5) =
    (box.PMin()(0) < pmin(0)) && (pmin(0) < box.PMax()(0));
  surfaceactive.Elem(6) =
    (box.PMin()(0) < pmax(0)) && (pmax(0) < box.PMax()(0));
}
}
