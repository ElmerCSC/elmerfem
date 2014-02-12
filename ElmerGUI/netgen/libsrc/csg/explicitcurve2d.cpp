#include <mystdlib.h>
#include <csg.hpp>

namespace netgen
{
ExplicitCurve2d :: ExplicitCurve2d ()
  {
    ;
  }
  
  
void ExplicitCurve2d :: Project (Point<2> & p) const
  {
  double t;
  t = ProjectParam (p);
  p = Eval (t);
  }

double ExplicitCurve2d :: NumericalProjectParam (const Point<2> & p, double lb, double ub) const
  {
  double t(-1);
  Vec<2> tan;
  Vec<2> curv;
  Point<2> cp;
  double f, fl, fu;
  int cnt;
  
  tan = EvalPrime (lb);
  cp = Eval (lb);
  fl = tan * (cp - p);
  if (fl > 0)			// changed by wmf, originally fl >= 0
    {
      //      cerr << "tan = " << tan << " cp - p = " << (cp - p) << endl;
      //      cerr << "ExplicitCurve2d::NumericalProject: lb wrong" << endl;
      return 0;
    }
  
  tan = EvalPrime (ub);
  cp = Eval (ub);
  fu = tan * (cp - p);
  if (fu < 0)			// changed by wmf, originally fu <= 0
    {
      //    cerr << "tan = " << tan << " cp - p = " << (cp - p) << endl;
      //    cerr << "ExplicitCurve2d::NumericalProject: ub wrong" << endl;
    return 0;
    }
    
  cnt = 0;
  while (ub - lb > 1e-12 && fu - fl > 1e-12)
    {
    cnt++;
    if (cnt > 50)
      {
      (*testout) << "Num Proj, cnt = " << cnt << endl;
      }
     
    t = (lb * fu - ub * fl) / (fu - fl);
    if (t > 0.9 * ub + 0.1 * lb) t = 0.9 * ub + 0.1 * lb;
    if (t < 0.1 * ub + 0.9 * lb) t = 0.1 * ub + 0.9 * lb;
    
    tan = EvalPrime (t);
    cp = Eval (t);
    f = tan * (cp - p);
    
    if (f >= 0)
      {
      ub = t;
      fu = f;
      }
    else
      {
      lb = t;
      fl = f;
      }
    }
    
  return t;
  }


Vec<2> ExplicitCurve2d :: Normal (double t) const
{
  Vec<2> tan = EvalPrime (t);
  tan.Normalize();
  return Vec<2> (tan(1), -tan(0));
}


void ExplicitCurve2d :: NormalVector (const Point<2> & p, Vec<2> & n) const
  {
  double t = ProjectParam (p);
  n = Normal (t);
  }


Point<2> ExplicitCurve2d :: CurvCircle (double t) const
  {
  Point<2> cp;
  Vec<2> tan, n, curv;
  double den;
  
  cp = Eval (t);
  tan = EvalPrime (t);
  n = Normal (t);
  curv = EvalPrimePrime (t);
  
  den = n * curv;
  if (fabs (den) < 1e-12)
    return cp + 1e12 * n;  
    
  return cp + (tan.Length2() / den) * n;  
  }


double ExplicitCurve2d :: MaxCurvature () const
  {
  double t, tmin, tmax, dt;
  double curv;
  Vec<2> tan;
  double maxcurv;

  maxcurv = 0;  
  
  tmin = MinParam ();
  tmax = MaxParam ();
  dt = (tmax - tmin) / 1000;
  for (t = tmin; t <= tmax+dt; t += dt)
    if (SectionUsed (t))
      {
      tan = EvalPrime (t);
      curv = fabs ( (Normal(t) * EvalPrimePrime(t)) / tan.Length2());
      if (curv > maxcurv) maxcurv = curv; 
      }
  return maxcurv;
  }  
  
double ExplicitCurve2d :: MaxCurvatureLoc (const Point<2> & p, double rad) const
  {
  double t, tmin, tmax, dt;
  double curv;
  Vec<2> tan;
  double maxcurv;

  maxcurv = 0;  
  
  tmin = MinParam ();
  tmax = MaxParam ();
  dt = (tmax - tmin) / 1000;
  for (t = tmin; t <= tmax+dt; t += dt)
    if (Dist (Eval(t), p) < rad)
      {
      tan = EvalPrime (t);
      curv = fabs ( (Normal(t) * EvalPrimePrime(t)) / tan.Length2());
      if (curv > maxcurv) maxcurv = curv; 
      }
    
  return maxcurv;
  }  
  
}
