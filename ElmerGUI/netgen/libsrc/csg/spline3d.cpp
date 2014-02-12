#include <mystdlib.h>

#include <myadt.hpp>

#include <linalg.hpp>
#include <csg.hpp>


namespace netgen
{
splinesegment3d :: splinesegment3d (const Point<3> & ap1, const Point<3> & ap2, 
				    const Point<3> & ap3)
{
  p1 = ap1;
  p2 = ap2;
  p3 = ap3;
}


/*
  todo
  Tip von Joerg Stiller:
  setzt Du in 
  void splinesegment3d :: Evaluate
  Zeilen 54 und 56
  b2 = 2 * t * (1-t);
  b2 /= sqrt(2);
  Das heisst, Du wichtest das zweite Bersteinpolynom mit 
  w2 = 1 / sqrt(2);
  Das ist aber nur fuer 45-Grad-Segmente korrekt. Fuer den
  allgemeinen Fall funktioniert
  w2 = ( e(p3 - p1), e(p2 - p1) );  // also cos(winkel(p3-p1, p2-p1))
  bzw. schoen symmetrisch
  w2 = ( e(p3 - p1), e(p2 - p1) )/2 + ( e(p1 - p3), e(p2 - p3) )/2;
  Das ist natuerlich kein C++ Code sondern symbolisch, wobei
  e(p3 - p1)    ist der von p1 zu p3 zeigende Einheitsvektor und
  (a, b)        steht fuer das Skalarprodukt zweier Vektoren etc.

  Eine vergleichbare Information steht auch irgendwo im Hoscheck & Lasser.
  Ich habe das Buch aber eben nicht zur Hand.
*/

void splinesegment3d :: Evaluate (double t, Point<3> & p) const
{
  double x, y, z, w;
  double b1, b2, b3;

  b1 = (1-t)*(1-t);
  b2 = 2 * t * (1-t);
  b3 = t * t;

  b2 /= sqrt(double(2));

  x = p1(0) * b1 + p2(0) * b2 + p3(0) * b3;
  y = p1(1) * b1 + p2(1) * b2 + p3(1) * b3;
  z = p1(2) * b1 + p2(2) * b2 + p3(2) * b3;
  w = b1 + b2 + b3;

  p(0) = x / w;
  p(1) = y / w;
  p(2) = z / w;
}

void splinesegment3d :: EvaluateTangent (double t, Vec<3> & tang) const
{
  double x, y, z, w, xprime, yprime, zprime, wprime;
  double b1, b2, b3, b1prime, b2prime, b3prime;

  b1 = (1-t)*(1-t);
  b2 = 2 * t * (1-t);
  b3 = t * t;
  b2 /= sqrt(double(2));

  b1prime = 2 * t - 2;
  b2prime = - 4 * t + 2;
  b3prime = 2 * t;
  b2prime /= sqrt(double(2));

 
  x = p1(0) * b1 + p2(0) * b2 + p3(0) * b3;
  y = p1(1) * b1 + p2(1) * b2 + p3(1) * b3;
  z = p1(2) * b1 + p2(2) * b2 + p3(2) * b3;
  w = b1 + b2 + b3;

  xprime = p1(0) * b1prime + p2(0) * b2prime + p3(0) * b3prime;
  yprime = p1(1) * b1prime + p2(1) * b2prime + p3(1) * b3prime;
  zprime = p1(2) * b1prime + p2(2) * b2prime + p3(2) * b3prime;
  wprime = b1prime + b2prime + b3prime;

  tang(0) = (w * xprime - x * wprime) / (w * w);
  tang(1) = (w * yprime - y * wprime) / (w * w);
  tang(2) = (w * zprime - z * wprime) / (w * w);
}
 

void spline3d :: AddSegment (const Point<3> & ap1, const Point<3> & ap2, 
			     const Point<3> & ap3)
{
  segments.Append (new splinesegment3d (ap1, ap2, ap3));
}

void spline3d :: Evaluate (double t, Point<3> & p) const
{
  int nr;
  double loct;
  static int cnt = 0;
  
  cnt++;
  if (cnt % 10000 == 0) (*mycout) << "Evaluate calls: " << cnt << endl;

  while (t < 0) t += GetNumSegments();
  while (t >= GetNumSegments()) t -= GetNumSegments();
  nr = 1 + int (t);
  loct = t - nr + 1;
  segments.Get(nr)->Evaluate (loct, p);
}
  
void spline3d :: EvaluateTangent (double t, Vec<3> & tang) const
{
  int nr;
  double loct;

  while (t < 0) t += GetNumSegments();
  while (t >= GetNumSegments()) t -= GetNumSegments();
  nr = 1 + int (t);
  loct = t - nr + 1;
  segments.Get(nr)->EvaluateTangent (loct, tang);
}


double spline3d :: ProjectToSpline (Point<3> & p) const
{
  double t, tl, tu, dt, dist, mindist, optt(0);
  Point<3> hp;
  Vec<3> tanx, px;
  
  dt = 0.01;
  mindist = 0;
  for (t = 0; t <= GetNumSegments() + dt/2; t += dt)
    {
      Evaluate (t, hp);
      dist = Dist (hp, p);
      if (t == 0 || dist < mindist)
	{
	  optt = t;
	  mindist = dist;
	} 
    }

  
  tu = optt + dt;
  tl = optt - dt;
  while (tu - tl > 1e-2)
    {
      optt = 0.5 * (tu + tl);
      Evaluate (optt, hp);
      EvaluateTangent (optt, tanx);
      if (tanx * (hp - p) > 0)
	tu = optt;
      else
	tl = optt;
    } 

  optt = 0.5 * (tu + tl);

  optt = ProjectToSpline (p, optt);
  return optt;
}
 
 
double spline3d :: ProjectToSpline (Point<3> & p, double optt) const
{ 
  double tl, tu, dt, val, dval, valu, vall;
  Point<3> hp;
  Vec<3> tanx, px;
  int its = 0;
  int cnt = 1000;
  do
    {
      dt = 1e-8;
      tl = optt - dt;
      tu = optt + dt;
    
      EvaluateTangent (optt, tanx); 
      Evaluate (optt, hp);
      px = hp - p;
      val =  px * tanx;
    
      EvaluateTangent (tl, tanx); 
      Evaluate (tl, hp);
      px = hp - p;
      vall =  px * tanx;
    
      EvaluateTangent (tu, tanx); 
      Evaluate (tu, hp);
      px = hp - p;
      valu =  px * tanx;
    
      dval = (valu - vall) / (2 * dt);

      if (its % 100 == 99)    
	(*testout) << "optt = " << optt 
		   << " val = " << val 
		   << " dval = " << dval << endl;
      optt -= val / dval;
      its++;
      if (fabs(val) < 1e-8 && cnt > 5) cnt = 5;
      cnt--;
    }
  while (cnt > 0);
        
  Evaluate (optt, p);
  return optt;
}
  
  
splinetube :: splinetube (const spline3d & amiddlecurve, double ar)
  : Surface(), middlecurve (amiddlecurve), r(ar)
{
  (*mycout) << "Splinetube Allocated, r = " << r << endl;

}
  
void splinetube :: DefineTangentialPlane (const Point<3> & ap1, 
					  const Point<3> & ap2)
{
  double t;
  double phi, z;
  
  p1 = ap1;
  p2 = ap2;
  cp = p1;
  t = middlecurve.ProjectToSpline (cp);
  ex = p1 - cp;
  middlecurve.EvaluateTangent (t, ez); 
  ex.Normalize();
  ez.Normalize();
  ey = Cross (ez, ex);
  
  phi = r * atan2 (ey * (p2-cp), ex * (p2-cp));
  z = ez * (p2 - cp); 
  e2x(0) = phi;
  e2x(1) = z;
  e2x.Normalize();
  e2y(1) = e2x(0);
  e2y(0) = -e2x(1);
  
  //  (*testout) << "Defineplane: " << endl
  //  	<< "p1 = " << p1 << "   p2 = " << p2 << endl
  //  	<< "pc = " << cp << endl
  //  	<< "ex = " << ex << " ey = " << ey << " ez = " << ez << endl
  //  	<< "phi = " << phi << "  z = " << z << endl
  //  	<< "e2x = " << e2x << " e2y = " << e2y << endl;
}
  
void splinetube :: ToPlane (const Point<3> & p3d, Point<2> & pplain, double h, 
			    int & zone) const
{
  Vec<2> v;
  v(0) = r * atan2 (ey * (p3d-cp), ex * (p3d-cp));
  v(1) = ez * (p3d - cp); 
  zone = 0;
  if (v(0) > r * 2) zone = 1;
  if (v(0) < r * 2) zone = 2;
  
  pplain(0) = (v * e2x) / h;
  pplain(1) = (v * e2y) / h;
}
  
void splinetube :: FromPlane (const Point<2> & pplain, Point<3> & p3d, double h) const
{
  Vec<2> v;
  
  v(0) = pplain(0) * h * e2x(0) + pplain(1) * h * e2y(0);
  v(1) = pplain(0) * h * e2x(1) + pplain(1) * h * e2y(1);
  
  p3d = p1 + v(0) * ey + v(1) * ez;

  Project (p3d);
}
  
void splinetube :: Project (Point<3> & p3d) const
{
  Point<3> hp;
  
  hp = p3d;
  middlecurve.ProjectToSpline (hp);
  
  p3d = hp + (r / Dist(p3d, hp)) * (p3d - hp); 
}



double splinetube :: CalcFunctionValue (const Point<3> & point) const
{
  Point<3> hcp;
  double rad;

  hcp = point;
  middlecurve.ProjectToSpline (hcp);
  rad = Dist (hcp, point);
  return 0.5 * (rad * rad / r - r);
}
  
void splinetube :: CalcGradient (const Point<3> & point, Vec<3> & grad) const
{
  Point<3> hcp;

  hcp = point;
  middlecurve.ProjectToSpline (hcp);

  grad = point - hcp;
  grad /= r;
}
  
  


Point<3> splinetube :: GetSurfacePoint () const
{
  Point<3> p;
  Vec<3> t, n;
  
  middlecurve.Evaluate (0, p);
  middlecurve.EvaluateTangent (0, t);
  n = t.GetNormal ();
  n *= r;
  (*mycout) << "p = " << p << " t = " << t << "  n = " << n << endl;
  return p + n;
}

void splinetube :: Print (ostream & str) const
{
  int i;
  str << "SplineTube, " 
      << middlecurve.GetNumSegments () << " segments, r = " << r << endl;
  for (i = 1; i <= middlecurve.GetNumSegments(); i++)
    str << middlecurve.P1(i) << " - " 
	<< middlecurve.P2(i) << " - " 
	<< middlecurve.P3(i) << endl;
}


int splinetube :: BoxInSolid (const BoxSphere<3> & box) const
  // 0 .. no, 1 .. yes, 2 .. maybe
{
  Point<3> pc = box.Center();
  middlecurve.ProjectToSpline (pc);
  double d = Dist (pc, box.Center());
  
  if (d < r - box.Diam()/2) return 1;
  if (d > r + box.Diam()/2) return 0;
  return 2;
}
}
