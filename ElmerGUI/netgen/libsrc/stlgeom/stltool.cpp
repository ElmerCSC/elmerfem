#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <gprim.hpp>

#include <meshing.hpp>

#include "stlgeom.hpp"

namespace netgen
{


//add a point into a pointlist, return pointnumber
int AddPointIfNotExists(ARRAY<Point3d>& ap, const Point3d& p, double eps)
{
  int i;
  for (i = 1; i <= ap.Size(); i++)
    {
      if (Dist(ap.Get(i),p) <= eps ) {return i;}
    }
  return ap.Append(p);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double GetDistFromLine(const Point<3> & lp1, const Point<3> & lp2, 
		       Point<3> & p)
{
  Vec3d vn = lp2 - lp1;
  Vec3d v1 = p - lp1;
  Vec3d v2 = lp2 - p;

  Point3d pold = p;

  if (v2 * vn <= 0) {p = lp2; return (pold - p).Length();}
  if (v1 * vn <= 0) {p = lp1; return (pold - p).Length();}
    
  double vnl = vn.Length();
  if (vnl == 0) {return Dist(lp1,p);}

  vn /= vnl;
  p = lp1 + (v1 * vn) * vn;
  return (pold - p).Length();
};

double GetDistFromInfiniteLine(const Point<3>& lp1, const Point<3>& lp2, const Point<3>& p)
{
  Vec3d vn(lp1, lp2);
  Vec3d v1(lp1, p);

  double vnl = vn.Length();

  if (vnl == 0)
    {
      return Dist (lp1, p);
    }
  else
    {
      return Cross (vn, v1).Length() / vnl;
    }
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Binary IO-Manipulation



void FIOReadInt(istream& ios, int& i)
{
  const int ilen = sizeof(int);
  
  char buf[ilen];
  int j;
  for (j = 0; j < ilen; j++)
    {
      ios.get(buf[j]);
    }
  memcpy(&i, &buf, ilen);
}

void FIOWriteInt(ostream& ios, const int& i)
{
  const int ilen = sizeof(int);
  
  char buf[ilen];
  memcpy(&buf, &i, ilen);

  int j;
  for (j = 0; j < ilen; j++)
    {
      ios << buf[j];
    }
}

void FIOReadDouble(istream& ios, double& i)
{
  const int ilen = sizeof(double);
  
  char buf[ilen];
  int j;
  for (j = 0; j < ilen; j++)
    {
      ios.get(buf[j]);
    }
  memcpy(&i, &buf, ilen);
}

void FIOWriteDouble(ostream& ios, const double& i)
{
  const int ilen = sizeof(double);
  
  char buf[ilen];
  memcpy(&buf, &i, ilen);

  int j;
  for (j = 0; j < ilen; j++)
    {
      ios << buf[j];
    }
}

void FIOReadFloat(istream& ios, float& i)
{
  const int ilen = sizeof(float);
  
  char buf[ilen];
  int j;
  for (j = 0; j < ilen; j++)
    {
      ios.get(buf[j]);
    }
  memcpy(&i, &buf, ilen);
}

void FIOWriteFloat(ostream& ios, const float& i)
{
  const int ilen = sizeof(float);
  
  char buf[ilen];
  memcpy(&buf, &i, ilen);

  int j;
  for (j = 0; j < ilen; j++)
    {
      ios << buf[j];
     }
}

void FIOReadString(istream& ios, char* str, int len)
{
  int j;
  for (j = 0; j < len; j++)
    {
      ios.get(str[j]);
    }
}

//read string and add terminating 0
void FIOReadStringE(istream& ios, char* str, int len)
{
  int j;
  for (j = 0; j < len; j++)
    {
      ios.get(str[j]);
    }
  str[len] = 0;
}

void FIOWriteString(ostream& ios, char* str, int len)
{
  int j;
  for (j = 0; j < len; j++)
    {
      ios << str[j];
    }
}


/*
void FIOReadInt(istream& ios, int& i)
{
  const int ilen = sizeof(int);
  
  char buf[ilen];
  int j;
  for (j = 0; j < ilen; j++)
    {
      ios.get(buf[ilen-j-1]);
    }
  memcpy(&i, &buf, ilen);
}

void FIOWriteInt(ostream& ios, const int& i)
{
  const int ilen = sizeof(int);
  
  char buf[ilen];
  memcpy(&buf, &i, ilen);

  int j;
  for (j = 0; j < ilen; j++)
    {
      ios << buf[ilen-j-1];
    }
}

void FIOReadDouble(istream& ios, double& i)
{
  const int ilen = sizeof(double);
  
  char buf[ilen];
  int j;
  for (j = 0; j < ilen; j++)
    {
      ios.get(buf[ilen-j-1]);
    }
  memcpy(&i, &buf, ilen);
}

void FIOWriteDouble(ostream& ios, const double& i)
{
  const int ilen = sizeof(double);
  
  char buf[ilen];
  memcpy(&buf, &i, ilen);

  int j;
  for (j = 0; j < ilen; j++)
    {
      ios << buf[ilen-j-1];
    }
}

void FIOReadFloat(istream& ios, float& i)
{
  const int ilen = sizeof(float);
  
  char buf[ilen];
  int j;
  for (j = 0; j < ilen; j++)
    {
      ios.get(buf[ilen-j-1]);
    }
  memcpy(&i, &buf, ilen);
}

void FIOWriteFloat(ostream& ios, const float& i)
{
  const int ilen = sizeof(float);
  
  char buf[ilen];
  memcpy(&buf, &i, ilen);

  int j;
  for (j = 0; j < ilen; j++)
    {
      ios << buf[ilen-j-1];
    }
}

void FIOReadString(istream& ios, char* str, int len)
{
  int j;
  for (j = 0; j < len; j++)
    {
      ios.get(str[j]);
    }
}

//read string and add terminating 0
void FIOReadStringE(istream& ios, char* str, int len)
{
  int j;
  for (j = 0; j < len; j++)
    {
      ios.get(str[j]);
    }
  str[len] = 0;
}

void FIOWriteString(ostream& ios, char* str, int len)
{
  int j;
  for (j = 0; j < len; j++)
    {
      ios << str[j];
    }
}
*/

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

STLReadTriangle :: STLReadTriangle (const Point<3> * apts,
				    const Vec<3> & anormal)
{
  pts[0] = apts[0];
  pts[1] = apts[1];
  pts[2] = apts[2]; 
  normal = anormal;
}



STLTriangle :: STLTriangle(const int * apts)
{
  pts[0] = apts[0];
  pts[1] = apts[1];
  pts[2] = apts[2];

  facenum = 0;
}

int STLTriangle :: IsNeighbourFrom(const STLTriangle& t) const
{
  //triangles must have same orientation!!!
  int i, j;
  for(i = 0; i <= 2; i++)
    {
      for(j = 0; j <= 2; j++)
	{
	  if (t.pts[(i+1)%3] == pts[j] &&
	      t.pts[i] == pts[(j+1)%3])
	    {return 1;}
	}
    }
  return 0;      
}

int STLTriangle :: IsWrongNeighbourFrom(const STLTriangle& t) const
{
  //triangles have not same orientation!!!
  int i, j;
  for(i = 0; i <= 2; i++)
    {
      for(j = 0; j <= 2; j++)
	{
	  if (t.pts[(i+1)%3] == pts[(j+1)%3] &&
	      t.pts[i] == pts[j])
	    {return 1;}
	}
    }
  return 0;      
}

void STLTriangle :: GetNeighbourPoints(const STLTriangle& t, int& p1, int& p2) const
{
  int i, j;
  for(i = 1; i <= 3; i++)
    {
      for(j = 1; j <= 3; j++)
	{
	  if (t.PNumMod(i+1) == PNumMod(j) &&
	      t.PNumMod(i) == PNumMod(j+1))
	    {p1 = PNumMod(j); p2 = PNumMod(j+1); return;}
	}
    }
  PrintSysError("Get neighbourpoints failed!");
}

int STLTriangle :: GetNeighbourPointsAndOpposite(const STLTriangle& t, int& p1, int& p2, int& po) const
{
  int i, j;
  for(i = 1; i <= 3; i++)
    {
      for(j = 1; j <= 3; j++)
	{
	  if (t.PNumMod(i+1) == PNumMod(j) &&
	      t.PNumMod(i) == PNumMod(j+1))
	    {p1 = PNumMod(j); p2 = PNumMod(j+1); po = PNumMod(j+2); return 1;}
	}
    }
  return 0;
}

Vec<3> STLTriangle :: GeomNormal(const ARRAY<Point<3> >& ap) const
{
  const Point<3> & p1 = ap.Get(PNum(1));
  const Point<3> & p2 = ap.Get(PNum(2));
  const Point<3> & p3 = ap.Get(PNum(3));
  
  return Cross(p2-p1, p3-p1);
}


void STLTriangle :: SetNormal (const Vec<3> & n)
{
  double len = n.Length();
  if (len > 0)
    {
      normal = n;
      normal.Normalize();
    }
  else
    {
      normal = Vec<3> (1, 0, 0);
    }
}


void STLTriangle :: ChangeOrientation()
{ 
  normal *= -1;
  Swap(pts[0],pts[1]); 
}



double STLTriangle :: Area(const ARRAY<Point<3> >& ap) const
{
  return 0.5 * Cross(ap.Get(PNum(2))-ap.Get(PNum(1)), 
		     ap.Get(PNum(3))-ap.Get(PNum(1))).Length();
}

double STLTriangle :: MinHeight(const ARRAY<Point<3> >& ap) const
{
  double ml = MaxLength(ap);
  if (ml != 0) {return 2.*Area(ap)/ml;}
  PrintWarning("max Side Length of a triangle = 0!!!");
  return 0;
}

double STLTriangle :: MaxLength(const ARRAY<Point<3> >& ap) const
{
  return max3(Dist(ap.Get(PNum(1)),ap.Get(PNum(2))),
	      Dist(ap.Get(PNum(2)),ap.Get(PNum(3))),
	      Dist(ap.Get(PNum(3)),ap.Get(PNum(1))));
}

void STLTriangle :: ProjectInPlain(const ARRAY<Point<3> >& ap, 
				   const Vec<3> & n, Point<3> & pp) const
{
  const Point<3> & p1 = ap.Get(PNum(1));
  const Point<3> & p2 = ap.Get(PNum(2));
  const Point<3> & p3 = ap.Get(PNum(3));
  
  Vec<3> v1 = p2 - p1;
  Vec<3> v2 = p3 - p1;
  Vec<3> nt = Cross(v1, v2);

  double c = - (p1(0)*nt(0) + p1(1)*nt(1) + p1(2)*nt(2));

  double prod = n * nt;  

  if (fabs(prod) == 0) 
    {
      pp = Point<3>(1.E20,1.E20,1.E20); 
      return; 
    }

  double nfact = -(pp(0)*nt(0) + pp(1)*nt(1) + pp(2)*nt(2) + c) / (prod);
  pp = pp + (nfact) * n;

}


int STLTriangle :: ProjectInPlain (const ARRAY<Point<3> >& ap, 
				   const Vec<3> & nproj, 
				   Point<3> & pp, Vec<3> & lam) const
{
  const Point<3> & p1 = ap.Get(PNum(1));
  const Point<3> & p2 = ap.Get(PNum(2));
  const Point<3> & p3 = ap.Get(PNum(3));
  
  Vec<3> v1 = p2-p1;
  Vec<3> v2 = p3-p1;

  Mat<3> mat;
  for (int i = 0; i < 3; i++)
    {
      mat(i,0) = v1(i);
      mat(i,1) = v2(i);
      mat(i,2) = nproj(i);
    }

  int err = 0;
  mat.Solve (pp-p1, lam);
  //  int err = SolveLinearSystem (v1, v2, nproj, pp-p1, lam);

  if (!err)
    {
      //      pp = p1 + lam(0) * v1 + lam(1) * v2;

      pp(0) = p1(0) + lam(0) * v1(0) + lam(1) * v2(0);
      pp(1) = p1(1) + lam(0) * v1(1) + lam(1) * v2(1);
      pp(2) = p1(2) + lam(0) * v1(2) + lam(1) * v2(2);
    }
  return err;
}





void STLTriangle :: ProjectInPlain(const ARRAY<Point<3> >& ap, 
				   Point<3> & pp) const
{
  const Point<3> & p1 = ap.Get(PNum(1));
  const Point<3> & p2 = ap.Get(PNum(2));
  const Point<3> & p3 = ap.Get(PNum(3));
  
  Vec<3> v1 = p2 - p1;
  Vec<3> v2 = p3 - p1;
  Vec<3> nt = Cross(v1, v2);

  double c = - (p1(0)*nt(0) + p1(1)*nt(1) + p1(2)*nt(2));
  
  double prod = nt * nt;  

  double nfact = -(pp(0)*nt(0) + pp(1)*nt(1) + pp(2)*nt(2) + c) / (prod);

  pp = pp + (nfact) * nt;
}

int STLTriangle :: PointInside(const ARRAY<Point<3> > & ap, 
			       const Point<3> & pp) const
{
  const Point<3> & p1 = ap.Get(PNum(1));
  const Point<3> & p2 = ap.Get(PNum(2));
  const Point<3> & p3 = ap.Get(PNum(3));
  
  Vec<3> v1 = p2 - p1;
  Vec<3> v2 = p3 - p1;
  Vec<3> v  = pp - p1;
  double det, l1, l2;
  Vec<3> ex, ey, ez;


  ez = GeomNormal(ap);
  ez /= ez.Length();
  ex = v1;
  ex /= ex.Length();
  ey = Cross (ez, ex);
  
  Vec<2> v1p(v1*ex, v1*ey);
  Vec<2> v2p(v2*ex, v2*ey);
  Vec<2> vp(v*ex, v*ey);

  det = v2p(1) * v1p(0) - v2p(0) * v1p(1);

  if (fabs(det) == 0) {return 0;}
  
  l2 = (vp(1) * v1p(0) - vp(0) * v1p(1)) / det;
  
  if (v1p(0) != 0.)
    {
      l1 = (vp(0) - l2 * v2p(0)) / v1p(0);
    }
  else if (v1p(1) != 0.)
    {
      l1 = (vp(1) - l2 * v2p(1)) / v1p(1);
    }
  else {return 0;}
  
  if (l1 >= -1E-10 && l2 >= -1E-10 && l1 + l2 <= 1.+1E-10) {return 1;}
  return 0; 
}

double STLTriangle :: GetNearestPoint(const ARRAY<Point<3> >& ap, 
				      Point<3> & p3d) const
{
  Point<3> p = p3d;
  ProjectInPlain(ap, p);
  double dist = (p - p3d).Length();

  if (PointInside(ap, p)) {p3d = p; return dist;}
  else
    {
      Point<3> pf;
      double nearest = 1E50;
      //int fi = 0;
      int j;
      
      for (j = 1; j <= 3; j++)
	{
	  p = p3d;
	  dist = GetDistFromLine(ap.Get(PNum(j)), ap.Get(PNumMod(j+1)), p);
	  if (dist < nearest)
	    {
	      nearest = dist; 
	      pf = p;
	    }
	}
      p3d = pf;
      return nearest;
    }
}

int STLTriangle :: HasEdge(int p1, int p2) const
{
  int i;
  for (i = 1; i <= 3; i++)
    {
      if (p1 == PNum(i) && p2 == PNumMod(i+1)) {return 1;}
    }
  return 0;
}

ostream& operator<<(ostream& os, const STLTriangle& t)
{
  os << "[";
  os << t[0] << ",";
  os << t[1] << ",";
  os << t[2] << "]";

  return os;
};



STLTopEdge :: STLTopEdge ()
{
  pts[0] = pts[1] = 0;
  trigs[0] = trigs[1] = 0;
  cosangle = 1;
  status = ED_UNDEFINED;
}

STLTopEdge :: STLTopEdge (int p1, int p2, int trig1, int trig2)
{ 
  pts[0] = p1; 
  pts[1] = p2; 
  trigs[0] = trig1; 
  trigs[1] = trig2; 
  cosangle = 1;
  status = ED_UNDEFINED;
}




//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++   STL CHART   +++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

STLChart :: STLChart(STLGeometry * ageometry)
{
  charttrigs = new ARRAY<int> (0,0);
  outertrigs = new ARRAY<int> (0,0);
  ilimit = new ARRAY<twoint> (0,0);
  olimit = new ARRAY<twoint> (0,0);

  geometry = ageometry;

  if ( (stlparam.usesearchtree == 1))
    searchtree = new Box3dTree (geometry->GetBoundingBox().PMin() - Vec3d(1,1,1),
				geometry->GetBoundingBox().PMax() + Vec3d(1,1,1));
  else
    searchtree = NULL;
}

void STLChart :: AddChartTrig(int i)
{
  charttrigs->Append(i);
  
  const STLTriangle & trig = geometry->GetTriangle(i);
  const Point3d & p1 = geometry->GetPoint (trig.PNum(1));
  const Point3d & p2 = geometry->GetPoint (trig.PNum(2));
  const Point3d & p3 = geometry->GetPoint (trig.PNum(3));

  Point3d pmin(p1), pmax(p1);
  pmin.SetToMin (p2);
  pmin.SetToMin (p3);
  pmax.SetToMax (p2);
  pmax.SetToMax (p3);
  
  if (!geomsearchtreeon && (stlparam.usesearchtree == 1))
    {searchtree->Insert (pmin, pmax, i);}
}

void STLChart :: AddOuterTrig(int i)
{
  outertrigs->Append(i);

  const STLTriangle & trig = geometry->GetTriangle(i);
  const Point3d & p1 = geometry->GetPoint (trig.PNum(1));
  const Point3d & p2 = geometry->GetPoint (trig.PNum(2));
  const Point3d & p3 = geometry->GetPoint (trig.PNum(3));

  Point3d pmin(p1), pmax(p1);
  pmin.SetToMin (p2);
  pmin.SetToMin (p3);
  pmax.SetToMax (p2);
  pmax.SetToMax (p3);
  
  if (!geomsearchtreeon && (stlparam.usesearchtree==1))
    {searchtree->Insert (pmin, pmax, i);}
}

int STLChart :: IsInWholeChart(int nr) const
{
  int i;
  for (i = 1; i <= charttrigs->Size(); i++)
    {
      if (charttrigs->Get(i) == nr) {return 1;}
    }
  for (i = 1; i <= outertrigs->Size(); i++)
    {
      if (outertrigs->Get(i) == nr) {return 1;}
    }
  return 0;
}

void STLChart :: GetTrianglesInBox (const Point3d & pmin,
				    const Point3d & pmax,
				    ARRAY<int> & trias) const
{
  if (geomsearchtreeon) {PrintMessage(5,"geomsearchtreeon is set!!!");}

  if (searchtree)
    searchtree -> GetIntersecting (pmin, pmax, trias);
  else
    {
      int i;
      Box3d box1(pmin, pmax);
      box1.Increase (1e-4);
      Box3d box2;

      trias.SetSize(0);
      
      int nt = GetNT();
      for (i = 1; i <= nt; i++)
	{

	  int trignum = GetTrig(i);
	  const STLTriangle & trig = geometry->GetTriangle(trignum);
	  box2.SetPoint (geometry->GetPoint (trig.PNum(1)));
	  box2.AddPoint (geometry->GetPoint (trig.PNum(2)));
	  box2.AddPoint (geometry->GetPoint (trig.PNum(3)));
	  
	  if (box1.Intersect (box2))
	    {
	      trias.Append (trignum);
	    }
	}
    }
}

//trigs may contain the same triangle double
void STLChart :: MoveToOuterChart(const ARRAY<int>& trigs)
{
  if (!trigs.Size()) {return;}
  int i;
  for (i = 1; i <= trigs.Size(); i++)
    {
      if (charttrigs->Get(trigs.Get(i)) != -1) 
	{AddOuterTrig(charttrigs->Get(trigs.Get(i)));}
      charttrigs->Elem(trigs.Get(i)) = -1;
    }
  DelChartTrigs(trigs);
}

//trigs may contain the same triangle double
void STLChart :: DelChartTrigs(const ARRAY<int>& trigs)
{
  if (!trigs.Size()) {return;}

  int i;
  for (i = 1; i <= trigs.Size(); i++)
    {
      charttrigs->Elem(trigs.Get(i)) = -1;
    }

  int cnt = 0;
  for (i = 1; i <= charttrigs->Size(); i++)
    {
      if (charttrigs->Elem(i) == -1)
	{
	  cnt++;
	}
      if (cnt != 0 && i < charttrigs->Size())
	{
	  charttrigs->Elem(i-cnt+1) = charttrigs->Get(i+1);
	}
    }
  i = charttrigs->Size() - trigs.Size();
  charttrigs->SetSize(i);

  if (!geomsearchtreeon && stlparam.usesearchtree == 1)
    {
      PrintMessage(7, "Warning: unsecure routine due to first use of searchtrees!!!");
      //bould new searchtree!!!
      searchtree = new Box3dTree (geometry->GetBoundingBox().PMin() - Vec3d(1,1,1),
				  geometry->GetBoundingBox().PMax() + Vec3d(1,1,1));

      for (i = 1; i <= charttrigs->Size(); i++)
	{
	  const STLTriangle & trig = geometry->GetTriangle(i);
	  const Point3d & p1 = geometry->GetPoint (trig.PNum(1));
	  const Point3d & p2 = geometry->GetPoint (trig.PNum(2));
	  const Point3d & p3 = geometry->GetPoint (trig.PNum(3));
	  
	  Point3d pmin(p1), pmax(p1);
	  pmin.SetToMin (p2);
	  pmin.SetToMin (p3);
	  pmax.SetToMax (p2);
	  pmax.SetToMax (p3);
	  
	  searchtree->Insert (pmin, pmax, i);	  
	}
    }
}


void STLChart :: SetNormal (const Point<3> & apref, const Vec<3> & anormal)
{
  pref = apref;
  normal = anormal;
  double len = normal.Length();
  if (len) normal /= len;
  else normal = Vec<3> (1, 0, 0);

  t1 = normal.GetNormal ();
  t2 = Cross (normal, t1);
}

Point<2> STLChart :: Project2d (const Point<3> & p3d) const
{
  Vec<3> v = p3d-pref;
  return Point<2> (t1 * v, t2 * v);
}



/*
  Point3d p1, p2, center;
  double rad;
  int i1, i2;
public:
*/
STLBoundarySeg :: 
STLBoundarySeg (int ai1, int ai2, const ARRAY<Point<3> > & points,
		const STLChart * chart)
{
  i1 = ai1;
  i2 = ai2; 
  p1 = points.Get(i1);
  p2 = points.Get(i2);
  center = ::netgen::Center (p1, p2);
  rad = Dist (p1, center);

  p2d1 = chart->Project2d (p1);
  p2d2 = chart->Project2d (p2);
  
  boundingbox.Set (p2d1);
  boundingbox.Add (p2d2);
}

void STLBoundarySeg :: Swap ()
{
  ::netgen::Swap (i1, i2);
  ::netgen::Swap (p1, p2);
}



STLBoundary :: STLBoundary (STLGeometry * ageometry)
  : boundary(), geometry(ageometry)
{
  ;
}


void STLBoundary :: AddOrDelSegment(const STLBoundarySeg & seg)
{
  int i;
  int found = 0;
  for (i = 1; i <= boundary.Size(); i++)
    {
      if (found) {boundary.Elem(i-1) = boundary.Get(i);}
      if (boundary.Get(i) == seg) {found = 1;}
    }
  if (!found) 
    {
      boundary.Append(seg);
    }
  else 
    {
      boundary.SetSize(boundary.Size()-1);
    }
}

void STLBoundary ::AddTriangle(const STLTriangle & t)
{
  int i;
  int found1 = 0;
  int found2 = 0;
  int found3 = 0;
  //int offset = 0;
  

  STLBoundarySeg seg1(t[0],t[1], geometry->GetPoints(), chart);
  STLBoundarySeg seg2(t[1],t[2], geometry->GetPoints(), chart);
  STLBoundarySeg seg3(t[2],t[0], geometry->GetPoints(), chart);

  seg1.SetSmoothEdge (geometry->IsSmoothEdge (seg1.I1(), seg1.I2()));
  seg2.SetSmoothEdge (geometry->IsSmoothEdge (seg2.I1(), seg2.I2()));
  seg3.SetSmoothEdge (geometry->IsSmoothEdge (seg3.I1(), seg3.I2()));

  /*
  for (i = 1; i <= boundary.Size(); i++)
    {
      if (offset) {boundary.Elem(i-offset) = boundary.Get(i);}
      if (boundary.Get(i) == seg1) {found1 = 1; offset++;}
      if (boundary.Get(i) == seg2) {found2 = 1; offset++;}
      if (boundary.Get(i) == seg3) {found3 = 1; offset++;}
    }

  if (offset)
    {
      boundary.SetSize(boundary.Size()-offset);
    }    
  */
  for (i = boundary.Size(); i >= 1; i--)
    {
      if (boundary.Get(i) == seg1) 
	{ boundary.DeleteElement (i); found1 = 1; } 
      else if (boundary.Get(i) == seg2) 
	{ boundary.DeleteElement (i); found2 = 1; } 
      else if (boundary.Get(i) == seg3) 
	{ boundary.DeleteElement (i); found3 = 1; } 
    }

  if (!found1) {seg1.Swap(); boundary.Append(seg1);}
  if (!found2) {seg2.Swap(); boundary.Append(seg2);}
  if (!found3) {seg3.Swap(); boundary.Append(seg3);}
}

int STLBoundary :: TestSeg(const Point<3>& p1, const Point<3> & p2, const Vec<3> & sn, 
			   double sinchartangle, int divisions, ARRAY<Point<3> >& points, double eps)
{

  if (usechartnormal)
    return TestSegChartNV (p1, p2, sn);

  // for statistics
  {
    int i;
    static ARRAY<int> cntclass;
    static int cnt = 0;
    static int cnti = 0, cnto = 0;
    static long int cntsegs = 0;
    if (cntclass.Size() == 0)
      {
	cntclass.SetSize (20);
	for (i = 1; i <= cntclass.Size(); i++)
	  cntclass.Elem(i) = 0;
      }
    
    cntsegs += NOSegments();
    int cla = int (log (double(NOSegments()+1)) / log(2.0));
    if (cla < 1) cla = 1;
    if (cla > cntclass.Size()) cla = cntclass.Size();
    cntclass.Elem(cla)++;
    cnt++;
    if (divisions)
      cnti++;
    else
      cnto++;
    if (cnt > 100000) 
      {
	cnt = 0;
	/*
	(*testout) << "TestSeg-calls for classes:" << endl;
	(*testout) << cnti << " inner calls, " << cnto << " outercalls" << endl;
	(*testout) << "total testes segments: " << cntsegs << endl;
	for (i = 1; i <= cntclass.Size(); i++)
	  {
	    (*testout) << int (exp (i * log(2.0))) << " bnd segs: " << cntclass.Get(i) << endl;
	  }
	*/
      }
  }


  int i,j,k;
  Point<3> seg1p/*, seg2p*/;
  Point<3> sp1,sp2;
  double lambda1, lambda2, vlen2;
  Vec<3> vptpl;
  double sinchartangle2 = sqr(sinchartangle);
  double scal;
  int possible;

  //double maxval = -1;
  //double maxvalnew = -1;



  double scalp1 = p1(0) * sn(0) + p1(1) * sn(1) + p1(2) * sn(2);
  double scalp2 = p2(0) * sn(0) + p2(1) * sn(1) + p2(2) * sn(2);
  double minl = min2(scalp1, scalp2);
  double maxl = max2(scalp1, scalp2);
  Point<3> c = Center (p1, p2);
  double dist1 = Dist (c, p1);
 
  int nseg = NOSegments();
  for (j = 1; j <= nseg; j++)
    {
      const STLBoundarySeg & seg = GetSegment(j);


      if (seg.IsSmoothEdge())
	continue;


      sp1 = seg.P1();
      sp2 = seg.P2();

      // Test, ob Spiral Konfikt moeglich
      
      possible = 1;

      double scalsp1 = sp1(0) * sn(0) + sp1(1) * sn(1) + sp1(2) * sn(2);
      double scalsp2 = sp2(0) * sn(0) + sp2(1) * sn(1) + sp2(2) * sn(2);

      double minsl = min2(scalsp1, scalsp2);
      double maxsl = max2(scalsp1, scalsp2);
      
      double maxdiff = max2 (maxsl - minl, maxl - minsl);
      
      /*
      Point3d sc = Center (sp1, sp2);
      double mindist = Dist(c, sc) - dist1 - GetSegment(j).Radius();
      if (maxdiff < sinchartangle * mindist)
	{
	  possible = 0;
	}
      */
       
      double hscal = maxdiff + sinchartangle * (dist1 + seg.Radius());
      if (hscal * hscal < sinchartangle * Dist2(c, seg.center ))
	possible = 0;


      /*      
      if (possible)
	{
	  double mindist2ex = MinDistLL2 (p1, p2, sp1, sp2);
	  if (maxdiff * maxdiff < sinchartangle2 * mindist2ex)
	    possible = 0;
	}
      */

      if (possible)
      	{
	  LinearPolynomial2V lp (scalp1 - scalsp1,
				 scalp2 - scalp1,
				 -(scalsp2 - scalsp1));
	  QuadraticPolynomial2V slp;
	  slp.Square (lp);
	  
      
	  Vec3d v (p1, sp1);
	  Vec3d vl (p1, p2);
	  Vec3d vsl (sp1, sp2);
      
	  QuadraticPolynomial2V qp (v.Length2(),
				    -2 * (v * vl),
				    2 * (v * vsl),
				    vl.Length2(),
				    -2 * (vl * vsl),
				    vsl.Length2());
	  
	  slp.Add (-sinchartangle2, qp);

	  double hv = slp.MaxUnitSquare();

	  if (hv > eps) return 0;
	  /*
	  if (hv > maxvalnew)
	    maxvalnew = hv;
	  */
	}
      

      if (possible && 0)

	for (i = 0; i <= divisions; i++)
	  {
	    
	    lambda1 = (double)i/(double)divisions;
	    seg1p = Point3d(p1(0)*lambda1+p2(0)*(1.-lambda1),
			    p1(1)*lambda1+p2(1)*(1.-lambda1),
			    p1(2)*lambda1+p2(2)*(1.-lambda1));
	    

	    
	    for (k = 0; k <= divisions; k++)
	      {
		lambda2 = (double)k/(double)divisions;
		vptpl = Vec3d(sp1(0)*lambda2+sp2(0)*(1.-lambda2)-seg1p(0),
			      sp1(1)*lambda2+sp2(1)*(1.-lambda2)-seg1p(1),
			      sp1(2)*lambda2+sp2(2)*(1.-lambda2)-seg1p(2));
		
		vlen2 = vptpl.Length2();

		//		if (vlen2 > 0)
		  {
		    scal = vptpl * sn;
		    double hv = scal*scal - sinchartangle2*vlen2;



		    /*
		    if (hv > maxval)
		      maxval = hv;
		    */
		    if (hv > eps) return 0;
		  }
	      } 
	  }
    }
  
  return 1;
  //  return (maxvalnew < eps);
}



// checks, whether 2d projection intersects
int STLBoundary :: TestSegChartNV(const Point3d & p1, const Point3d& p2, 
				  const Vec3d& sn)
{
  int j;
  int nseg = NOSegments();

  Point<2> p2d1 = chart->Project2d (p1);
  Point<2> p2d2 = chart->Project2d (p2);

  Box<2> box2d;
  box2d.Set (p2d1);
  box2d.Add (p2d2);
  /*
  Point2d pmin(p2d1);
  pmin.SetToMin (p2d2);
  Point2d pmax(p2d1);
  pmax.SetToMax (p2d2);
  */

  Line2d l1 (p2d1, p2d2);

  double lam1, lam2;
  double eps = 1e-3;
  
  for (j = 1; j <= nseg; j++)
    {
      const STLBoundarySeg & seg = GetSegment(j);

      if (!box2d.Intersect (seg.BoundingBox()))
	continue;
      /*
      if (seg.P2DMin()(0) > pmax(0)) continue;
      if (seg.P2DMin()(1) > pmax(1)) continue;
      if (seg.P2DMax()(0) < pmin(0)) continue;
      if (seg.P2DMax()(1) < pmin(1)) continue;
      */

      if (seg.IsSmoothEdge()) continue;

      const Point<2> & sp1 = seg.P2D1();
      const Point<2> & sp2 = seg.P2D2();
      

      Line2d l2 (sp1, sp2);
      
      int err =
	CrossPointBarycentric (l1, l2, lam1, lam2);
      /*
      if (chartdebug)
	{
	  
	  (*testout) << "lam1 = " << lam1 << ", lam2 = " << lam2 << endl;
	  (*testout) << "p2d = " << p2d1 << ", " << p2d2 << endl;
	  (*testout) << "sp2d = " << sp1 << ", " << sp2 << endl;
	  (*testout) << "i1,2 = " << seg.I1() << ", " << seg.I2() << endl;
	  
	}
      */
      if (!err && lam1 > eps && lam1 < 1-eps &&
	  lam2 > eps && lam2 < 1-eps)
	return 0;
    }
  return 1;
}



STLDoctorParams :: STLDoctorParams()
{
  drawmeshededges = 1;
  geom_tol_fact = 1E-6;
  longlinefact = 0;
  showexcluded = 1;

  selectmode = 0;
  edgeselectmode = 0;
  useexternaledges = 0;
  showfaces = 0;
  showtouchedtrigchart = 1;
  showedgecornerpoints = 1;
  conecheck = 1;
  spiralcheck = 1;
  selecttrig = 0;
  nodeofseltrig = 1;
  selectwithmouse = 1;
  showmarkedtrigs = 1;
  dirtytrigfact = 0.001;
  smoothangle = 90;
  smoothnormalsweight = 0.2;
  vicinity = 0;
  showvicinity = 0;
}



STLDoctorParams stldoctor;

void STLDoctorParams :: Print (ostream & ost) const
{
  ost << "STL doctor parameters:" << endl
      << "selecttrig = " << selecttrig << endl
      << "selectlocalpoint = " << nodeofseltrig << endl
      << "selectwithmouse = " << selectwithmouse << endl
      << "showmarkedtrigs = " << showmarkedtrigs << endl
      << "dirtytrigfact = " << dirtytrigfact << endl
      << "smoothangle = " << smoothangle << endl;
}


STLParameters ::   STLParameters()
{
  yangle = 30;
  contyangle = 20;
  edgecornerangle = 60;
  chartangle = 15;
  outerchartangle = 70;
     
  usesearchtree = 0;
  atlasminh = 1E-4;
  resthsurfcurvfac = 2;
  resthsurfcurvenable = 0;
  resthatlasfac = 2;
  resthatlasenable = 1;
  resthchartdistfac = 1.2;
  resthchartdistenable = 1;
  resthlinelengthfac = 0.5;
  resthlinelengthenable = 1;
  resthcloseedgefac = 1;
  resthcloseedgeenable = 1;
  resthedgeanglefac = 1;
  resthedgeangleenable = 0;
  resthsurfmeshcurvfac = 1;
  resthsurfmeshcurvenable = 0;
  recalc_h_opt = 1;
}

void STLParameters :: Print (ostream & ost) const
{
  ost << "STL parameters:" << endl
      << "yellow angle = " << yangle << endl
      << "continued yellow angle = " << contyangle << endl
      << "edgecornerangle = " << edgecornerangle << endl
      << "chartangle = " << chartangle << endl
      << "outerchartangle = " << outerchartangle << endl
      << "restrict h due to ..., enable and safety factor: " << endl
      << "surface curvature: " << resthsurfcurvenable
      << ", fac = " << resthsurfcurvfac << endl
      << "atlas surface curvature: " << resthatlasenable
      << ", fac = " << resthatlasfac << endl
      << "chart distance: " << resthchartdistenable
      << ", fac = " << resthchartdistfac << endl
      << "line length: " << resthlinelengthenable
      << ", fac = " << resthlinelengthfac << endl
      << "close edges: " << resthcloseedgeenable
      << ", fac = " << resthcloseedgefac << endl
      << "edge angle: " << resthedgeangleenable
      << ", fac = " << resthedgeanglefac << endl;
}


STLParameters stlparam;


}
