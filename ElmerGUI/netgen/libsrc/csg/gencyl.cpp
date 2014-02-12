#include <linalg.hpp>
#include <csg.hpp>


namespace netgen
{

GeneralizedCylinder :: GeneralizedCylinder (ExplicitCurve2d & acrosssection,
					    Point<3> ap, Vec<3> ae1, Vec<3> ae2)
  : crosssection(acrosssection)
{
  planep = ap;
  planee1 = ae1;
  planee2 = ae2;
  planee3 = Cross (planee1, planee2);
  (*testout) << "Vecs = " << planee1 << " " << planee2 << " " << planee3 << endl;
};
  

void GeneralizedCylinder :: Project (Point<3> & p) const
{
  Point<2> p2d;
  double z;
  
  p2d = Point<2> (planee1 * (p - planep), planee2 * (p - planep));
  z = planee3 * (p - planep);

  crosssection.Project (p2d);
  
  p = planep + p2d(0) * planee1 + p2d(1) * planee2 + z * planee3;
}

int GeneralizedCylinder ::BoxInSolid (const BoxSphere<3> & box) const
{
  Point<3> p3d;
  Point<2> p2d, projp;
  double t;
  Vec<2> tan, n;
  
  p3d = box.Center();
  
  p2d = Point<2> (planee1 * (p3d - planep), planee2 * (p3d - planep));
  t = crosssection.ProjectParam (p2d);
  
  projp = crosssection.Eval (t);
  tan = crosssection.EvalPrime (t);
  n(0) = tan(1);
  n(1) = -tan(0);
    
  if (Dist (p2d, projp) < box.Diam()/2)
    return 2;
    
  if (n * (p2d - projp) > 0) 
    {
      return 0;   
    }
    
  return 1;
}

double GeneralizedCylinder :: CalcFunctionValue (const Point<3> & point) const
{
  Point<2> p2d, projp;
  double t;
  Vec<2> tan, n;
  
  
  p2d = Point<2> (planee1 * (point - planep), planee2 * (point - planep));
  t = crosssection.ProjectParam (p2d);
  
  projp = crosssection.Eval (t);
  tan = crosssection.EvalPrime (t);
  n(0) = tan(1);
  n(1) = -tan(0);
    
  n /= n.Length();
  return n * (p2d - projp);
}
  
void GeneralizedCylinder :: CalcGradient (const Point<3> & point, Vec<3> & grad) const
{
  Point<2> p2d, projp;
  double t;
  Vec<2> tan, n;
  
  
  p2d = Point<2> (planee1 * (point - planep), planee2 * (point - planep));
  t = crosssection.ProjectParam (p2d);
  
  projp = crosssection.Eval (t);
  tan = crosssection.EvalPrime (t);
  n(0) = tan(1);
  n(1) = -tan(0);
    
  n /= n.Length();
  grad = n(0) * planee1 + n(1) * planee2;
}
  
  
void GeneralizedCylinder :: CalcHesse (const Point<3> & point, Mat<3> & hesse) const
{
  Point<2> p2d, projp;
  double t, dist, val;
  Point<2> curvp;
  Vec<2> curvpp;
  Mat<2> h2d;
  Mat<3,2> vmat;
  int i, j, k, l;
  
  p2d = Point<2> (planee1 * (point - planep), planee2 * (point - planep));
  t = crosssection.ProjectParam (p2d);

  curvp = crosssection.CurvCircle (t);
  curvpp = p2d-curvp;
  dist = curvpp.Length();
  curvpp /= dist;
    
  h2d(1, 1) = (1 - curvpp(0) * curvpp(0) ) / dist;  
  h2d(1, 2) = h2d(2, 1) = (- curvpp(0) * curvpp(1) ) / dist;  
  h2d(2, 2) = (1 - curvpp(1) * curvpp(1) ) / dist;  
  
  vmat(0,0) = planee1(0);
  vmat(1,0) = planee1(1);
  vmat(2,0) = planee1(2);
  vmat(0,1) = planee2(0);
  vmat(1,1) = planee2(1);
  vmat(2,1) = planee2(2);
  
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      {
	val = 0;
	for (k = 0; k < 2; k++)
	  for (l = 0; l < 2; l++)
	    val += vmat(i,k) * h2d(k,l) * vmat(j,l);
	hesse(i,j) = val;
      }
}


double GeneralizedCylinder :: HesseNorm () const
{
  return crosssection.MaxCurvature();
}

double GeneralizedCylinder :: MaxCurvatureLoc (const Point<3> & c, double rad) const
{
  Point<2> c2d = Point<2> (planee1 * (c - planep), planee2 * (c - planep));
  return crosssection.MaxCurvatureLoc(c2d, rad);
}
  

  
Point<3> GeneralizedCylinder :: GetSurfacePoint () const
{
  Point<2> p2d; 
  p2d = crosssection.Eval(0);
  return planep + p2d(0) * planee1 + p2d(1) * planee2;
}

void GeneralizedCylinder :: Reduce (const BoxSphere<3> & box)
{
  Point<2> c2d = Point<2> (planee1 * (box.Center() - planep), 
			   planee2 * (box.Center() - planep));
  crosssection.Reduce (c2d, box.Diam()/2);
}

void GeneralizedCylinder :: UnReduce ()
{
  crosssection.UnReduce ();
}

void GeneralizedCylinder :: Print (ostream & str) const
{
  str << "Generalized Cylinder" << endl;
  crosssection.Print (str);
}
  
#ifdef MYGRAPH  
void GeneralizedCylinder :: Plot (const class ROT3D & rot) const
{
  Point<2> p2d;
  Point<3> p, oldp;
  double t, tmin, tmax, dt;
  
  tmin = crosssection.MinParam();
  tmax = crosssection.MaxParam();
  dt = (tmax - tmin)/ 500;
  
  p2d = crosssection.Eval(tmin);
  p = planep + p2d(0) * planee1 + p2d(1) * planee2;
  
  for (t = tmin; t <= tmax+dt; t += dt)
    {
      if (crosssection.SectionUsed (t))
	MySetColor (RED);
      else
	MySetColor (BLUE);
      
      oldp = p;
      p2d = crosssection.Eval(t);
      p = planep + p2d(0) * planee1 + p2d(1) * planee2;
      MyLine3D (p, oldp, rot);
    }

}

#endif  
}
