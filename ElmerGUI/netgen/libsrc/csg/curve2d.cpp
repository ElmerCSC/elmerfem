#include <mystdlib.h>

#include <myadt.hpp>
#include <csg.hpp>

namespace netgen
{
CircleCurve2d :: CircleCurve2d (const Point<2> & acenter, double arad)
  {
  center = acenter;
  rad = arad;
  }
  
void CircleCurve2d :: Project (Point<2> & p) const
  {
  Vec<2> v = p - center;
  v *= rad/v.Length();
  p = center + v;
  }

void CircleCurve2d :: NormalVector (const Point<2> & p, Vec<2> & n) const
  {
  n = p - center;
  n /= n.Length();
  }






QuadraticCurve2d ::  QuadraticCurve2d ()
{
  cxx = cyy = cxy = cx = cy = c = 0;
}

void QuadraticCurve2d :: Read (istream & ist)
{
  ist >> cxx >> cyy >> cxy >> cx >> cy >> c;
}


void QuadraticCurve2d :: Project (Point<2> & p) const
{
  double f, x, y, gradx, grady, grad2;
  int its = 0;

  x = p(0);
  y = p(1);

  do
    {
      f = cxx * x * x + cyy * y * y + cxy * x * y + cx * x + cy * y + c;
      gradx = 2 * cxx * x + cxy * y + cx;
      grady = 2 * cyy * y + cxy * x + cy;
      grad2 = gradx * gradx + grady * grady;
      
      x -= f * gradx / grad2;
      y -= f * grady / grad2;

      //      (*mycout) << "x = " << x << " y = " << y << " f = " << f << endl;
      its++;
    }
  while (fabs (f) > 1e-8 && its < 20);
  if (its >= 20)
    cerr << "QuadraticCurve2d::Project:  many iterations, f = " << f << endl;
  p(0) = x;
  p(1) = y;
}


void QuadraticCurve2d :: NormalVector (const Point<2> & p, Vec<2> & n) const
{
  n(0) = 2 * cxx * p(0) + cxy * p(1) + cx;
  n(1) = 2 * cyy * p(1) + cxy * p(0) + cy;
  n.Normalize();
}
}
