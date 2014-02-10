#include <mystdlib.h>

#include <myadt.hpp>
#include <gprim.hpp>

#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif

namespace netgen
{

ostream & operator<<(ostream  & s, const Point2d & p)
{
  return s << "(" << p.px << ", " << p.py << ")";
}

ostream & operator<<(ostream  & s, const Vec2d & v)
{
  return s << "(" << v.vx << ", " << v.vy << ")";
}

#ifdef none
ostream & operator<<(ostream  & s, const Line2d & l)
  {
  return s << l.p1 << "-" << l.p2;
}

ostream & operator<<(ostream  & s, const TRIANGLE2D & t)
{
  return s << t.p1 << "-" << t.p2 << "-" << t.p3;
}
#endif


double Fastatan2 (double x, double y)
{
  if (y > 0)
    {
      if (x > 0)
	return y / (x+y);
      else
	return 1 - x / (y-x);
    }
  else if (y < 0)
    {
      if (x < 0)
	return 2 + y / (x+y);
      else
	return 3 - x / (y-x);
    }
  else 
    {
      if (x >= 0)
	return 0;
      else
	return 2;
    }
}


double Angle (const Vec2d & v)
{
  if (v.X() == 0 && v.Y() == 0)
    return 0;
    
  double ang = atan2 (v.Y(), v.X());
  if (ang < 0) ang+= 2 * M_PI;
  return ang;
}

double FastAngle (const Vec2d & v)
{
  return Fastatan2 (v.X(), v.Y());
}

double Angle (const Vec2d & v1, const Vec2d & v2)
{
  double ang = Angle(v2) - Angle(v1);
  if (ang < 0) ang += 2 * M_PI;
  return ang;
}

double FastAngle (const Vec2d & v1, const Vec2d & v2)
{
  double ang = FastAngle(v2) - FastAngle(v1);
  if (ang < 0) ang += 4;
  return ang;
}

/*
int CW (const Point2d & p1,const Point2d & p2,const Point2d & p3)
{
  return Cross (p2 - p1, p3 - p2) < 0;
}

int CCW (const Point2d & p1,const Point2d & p2,const Point2d & p3)
{
  return Cross (p2 - p1, p3 - p2) > 0;
}
*/

double  Dist2(const Line2d & g, const Line2d & h )
  {
  double   dd = 0.0, d1,d2,d3,d4;
  Point2d  cp = CrossPoint(g,h);
  
  if ( Parallel(g,h) || !IsOnLine(g,cp) || !IsOnLine(h,cp) )
    {
      d1 = Dist2(g.P1(),h.P1());
      d2 = Dist2(g.P1(),h.P2());
      d3 = Dist2(g.P2(),h.P1());
      d4 = Dist2(g.P2(),h.P2());
      if (d1<d2)  d2 = d1;
      if (d3<d4)  d4 = d3;
      dd = ( d2 < d4 ) ? d2 : d4;
    }
  return dd;
}


Point2d CrossPoint (const Line2d & l1, const Line2d & l2)
  {
  double den = Cross (l1.Delta(), l2.Delta());
  double num = Cross ( (l2.P1() - l1.P1()), l2.Delta());

  if (den == 0)
    return l1.P1();
  else
    return l1.P1() + (num/den) * l1.Delta();
}


int CrossPointBarycentric (const Line2d & l1, const Line2d & l2,
			   double & lam1, double & lam2)
{
  // p = l1.1 + lam1 (l1.2-l1.1) = l2.1 + lam2 (l2.2-l2.1)
  double a11 = l1.p2.X() - l1.p1.X();
  double a21 = l1.p2.Y() - l1.p1.Y();
  double a12 = -(l2.p2.X() - l2.p1.X());
  double a22 = -(l2.p2.Y() - l2.p1.Y());

  double b1 = l2.p1.X() - l1.p1.X();
  double b2 = l2.p1.Y() - l1.p1.Y();
  
  double det = a11*a22 - a12 * a21;
  if (det == 0)
    return 1;
    
  lam1 = (a22 * b1 - a12 * b2) / det;
  lam2 = (a11 * b2 - a21 * b1) / det;
  return 0;
}




int Parallel (const Line2d & l1, const Line2d & l2, double peps)
  {
  double p = fabs (Cross (l1.Delta(), l2.Delta()));
  //  (*mycout) << endl << p << "  " <<  l1.Length() << "  " << l2.Length() << endl;
  return p <= peps * l1.Length() * l2.Length();
}

int IsOnLine (const Line2d & l, const Point2d & p, double heps)
  {
  double c1 = (p - l.P1()) * l.Delta();
  double c2 = (p - l.P2()) * l.Delta();
  double d = fabs (Cross ( (p - l.P1()), l.Delta()));
  double len2 = l.Length2();

  return c1 >= -heps * len2 && c2 <= heps * len2 && d <= heps * len2;
}

#ifdef none
int IsOnLine (const PLine2d & l, const Point2d & p, double heps)
  {
  double c1 = (p - l.P1()) * l.Delta();
  double c2 = (p - l.P2()) * l.Delta();
  double d = fabs (Cross ( (p - l.P1()), l.Delta()));
  double len2 = l.Length2();

  return c1 >= -heps * len2 && c2 <= heps * len2 && d <= heps * len2;
}

int IsOnLongLine (const Line2d & l, const Point2d & p)
  {
  double d = fabs (Cross ( (p - l.P1()), l.Delta()));
  return d <= EPSGEOM * l.Length();
}

int Hit (const Line2d & l1, const Line2d & l2, double heps)
  {
  double den =  Cross ( (l1.P2() - l1.P1()), (l2.P1() - l2.P2()));
  double num1 = Cross ( (l2.P1() - l1.P1()), (l2.P1() - l2.P2()));
  double num2 = Cross ( (l1.P2() - l1.P1()), (l2.P1() - l1.P1()));
  num1 *= sgn (den);
  num2 *= sgn (den);
  den = fabs (den);

  int ch = (-den * heps <= num1 && num1 <= den * (1 + heps) &&
	    -den * heps <= num2 && num2 <= den * (1 + heps));
  return ch;
}


void Line2d :: GetNormal (Line2d & n) const
{
  double 	ax  = P2().X()-P1().X(),
    ay  = P2().Y()-P1().Y();
  Point2d 	mid(P1().X()+.5*ax, P1().Y()+.5*ay);
 
 n=Line2d(mid,Point2d(mid.X()+ay,mid.Y()-ax)) ;
}

Vec2d Line2d :: NormalDelta () const
{
 Line2d tmp;
 GetNormal(tmp);
 return tmp.Delta();
}

int TRIANGLE2D :: IsOn (const Point2d & p) const
  {
  return IsOnLine (Line2d (p1, p2), p) ||
         IsOnLine (Line2d (p1, p3), p) ||
         IsOnLine (Line2d (p2, p3), p);
  }


int TRIANGLE2D :: IsIn (const Point2d & p) const
{
  return ::CW(p, p1, p2) == ::CW(p, p2, p3) &&
         ::CW(p, p1, p2) == ::CW(p, p3, p1);
}



int PTRIANGLE2D :: IsOn (const Point2d & p) const
{
  return IsOnLine (Line2d (*p1, *p2), p) ||
         IsOnLine (Line2d (*p1, *p3), p) ||
         IsOnLine (Line2d (*p2, *p3), p);
}


int PTRIANGLE2D :: IsIn (const Point2d & p) const
{
  return ::CW(p, *p1, *p2) == ::CW(p, *p2, *p3) &&
         ::CW(p, *p1, *p2) == ::CW(p, *p3, *p1);
}

#endif






Polygon2d :: Polygon2d ()
{
  ;
}

Polygon2d :: ~Polygon2d ()
{
  ;
}

void Polygon2d :: AddPoint (const Point2d & p)
{ 
  points.Append(p); 
}


double Polygon2d :: HArea () const
{
  int i;
  double ar = 0;
  for (i = 1; i <= points.Size(); i++)
    {
      const Point2d & p1 = points.Get(i);
      const Point2d & p2 = points.Get(i%points.Size()+1);
      ar += 
	(p2.X()-p1.X()) * p1.Y() -
	(p2.Y()-p1.Y()) * p1.X();
    }
  return ar/2;
  /*
  CURSOR c;
  double ar = 0;
  Point2d * p1, * p2, p0 = Point2d(0, 0);
  Vec2d v1, v2 = Vec2d(1, 0);

  p2 = points[points.Last()];
  for (c = points.First(); c != points.Head(); c++)
    {
    p1 = p2;
    p2 = points[c];
    ar += Cross ( (*p2-*p1), (*p1 - p0));
    }
  return ar / 2;
  */
}


int Polygon2d :: IsOn (const Point2d & p) const
{
  int i;
  for (i = 1; i <= points.Size(); i++)
    {
      const Point2d & p1 = points.Get(i);
      const Point2d & p2 = points.Get(i%points.Size()+1);
      if (IsOnLine (Line2d(p1, p2), p)) return 1;
    }
  return 0;
  /*
  CURSOR c;
  Point2d * p1, * p2;
  
  p2 = points[points.Last()];
  for (c = points.First(); c != points.Head(); c++)
    {
      p1 = p2;
      p2 = points[c];
      if (IsOnLine (Line2d(*p1, *p2), p)) return 1;
    }
  return 0;
  */
}


int Polygon2d :: IsIn (const Point2d & p) const
{
  int i;
  double sum = 0, ang;
  for (i = 1; i <= points.Size(); i++)
    {
      const Point2d & p1 = points.Get(i);
      const Point2d & p2 = points.Get(i%points.Size()+1);
      ang = Angle ( (p1 - p), (p2 - p) );
      if (ang > M_PI) ang -= 2 * M_PI;
      sum += ang;
    }
  return fabs(sum) > M_PI;
  /*
  CURSOR c;
  Point2d * p1, * p2;
  double sum = 0, ang;

  p2 = points[points.Last()];
  for (c = points.First(); c != points.Head(); c++)
    {
    p1 = p2;
    p2 = points[c];
    ang = Angle ( (*p1 - p), (*p2 - p) );
    if (ang > M_PI) ang -= 2 * M_PI;
    sum += ang;
    }

  return fabs(sum) > M_PI;
  */
}

int Polygon2d :: IsConvex () const
  {
    /*
  Point2d *p, *pold, *pnew;
  char cw;
  CURSOR c;

  if (points.Length() < 3) return 0;

  c = points.Last();
  p = points[c];
  c--;
  pold = points[c];
  pnew = points[points.First()];
  cw = ::CW (*pold, *p, *pnew);

  for (c = points.First(); c != points.Head(); c++)
    {
    pnew = points[c];
    if (cw != ::CW (*pold, *p, *pnew))
      return 0;
    pold = p;
    p = pnew;
    }
    */
    return 0;
  }


int Polygon2d :: IsStarPoint (const Point2d & p) const
  {
    /*
  Point2d *pnew, *pold;
  char cw;
  CURSOR c;

  if (points.Length() < 3) return 0;

  pold = points[points.Last()];
  pnew = points[points.First()];

  cw = ::CW (p, *pold, *pnew);

  for (c = points.First(); c != points.Head(); c++)
    {
    pnew = points[c];
    if (cw != ::CW (p, *pold, *pnew))
      return 0;
    pold = pnew;
    }
  return 1;
    */
    return 0;
  }


Point2d Polygon2d :: Center () const
  {
    /*
  double ai, a = 0, x = 0, y = 0;
  Point2d * p, *p2;
  Point2d p0 = Point2d(0, 0);
  CURSOR c;

  p2 = points[points.Last()];

  for (c = points.First(); c != points.Head(); c++)
    {
    p = points[c];
    ai = Cross (*p2 - p0, *p - p0);
    x += ai / 3 * (p2->X() + p->X());
    y += ai / 3 * (p2->Y() + p->Y());
    a+= ai;
    p2 = p;
    }
  if (a != 0)
    return Point2d (x / a, y / a);
  else
    return Point2d (0, 0);
    */
    return Point2d (0, 0);
  }



Point2d Polygon2d :: EqualAreaPoint () const
  {
    /*
  double a11 = 0, a12 = 0, a21= 0, a22 = 0;
  double b1 = 0, b2 = 0, dx, dy;
  double det;
  Point2d * p, *p2;
  CURSOR c;

  p = points[points.Last()];

  for (c = points.First(); c != points.Head(); c++)
    {
    p2 = p;
    p = points[c];

    dx = p->X() - p2->X();
    dy = p->Y() - p2->Y();

    a11 += sqr (dy);
    a12 -= dx * dy;
    a21 -= dx * dy;
    a22 += sqr (dx);
    b1 -= dy * (p->X() * p2->Y() - p2->X() * p->Y());
    b2 -= dx * (p->Y() * p2->X() - p2->Y() * p->X());
    }

  det = a11 * a22 - a21 * a12;

  if (det != 0)
    return Point2d ( (b1 * a22 - b2 * a12) / det,
                     (a11 * b2 - a21 * b1) / det);
  else
    return Point2d (0, 0);
*/
    return Point2d (0, 0);
  }


}
