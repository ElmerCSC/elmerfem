#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <gprim.hpp>

#include <meshing.hpp>

#include "stlgeom.hpp"

namespace netgen
{

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++  EDGE DATA     ++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/*
void STLEdgeData :: Write(ofstream& of) const
{
  of // << angle << " "
     << p1 << " "
     << p2 << " "
     << lt << " "
     << rt << " "
    //     << status
     << endl;
}

void STLEdgeData :: Read(ifstream& ifs)
{
  // ifs >> angle;
  ifs >> p1;
  ifs >> p2;
  ifs >> lt;
  ifs >> rt;
  //  ifs >> status;
}


int STLEdgeData :: GetStatus () const
{
  if (topedgenr <= 0 || topedgenr > top->GetNTE()) return 0;
  return top->GetTopEdge (topedgenr).GetStatus(); 
}

void STLEdgeData ::SetStatus (int stat)
{
  if (topedgenr >= 1 && topedgenr <= top->GetNTE())
    top->GetTopEdge (topedgenr).SetStatus(stat); 
}


float STLEdgeData :: CosAngle() const
{
  return top->GetTopEdge (topedgenr).CosAngle(); 
}



void STLEdgeDataList :: ResetAll()
{
  int i;
  for (i = 1; i <= edgedata.Size(); i++)
    {
      edgedata.Elem(i).SetUndefined();
    }
}

void STLEdgeDataList :: ResetCandidates()
{
  int i;
  for (i = 1; i <= edgedata.Size(); i++)
    {
      if (edgedata.Get(i).Candidate())
	{edgedata.Elem(i).SetUndefined();}
    }
}

int STLEdgeDataList :: GetNConfEdges() const
{
  int i;
  int cnt = 0;
  for (i = 1; i <= edgedata.Size(); i++)
    {
      if (edgedata.Get(i).Confirmed()) {cnt++;}
    }
  return cnt;
}

void STLEdgeDataList :: ConfirmCandidates()
{
  int i;
  for (i = 1; i <= edgedata.Size(); i++)
    {
      if (edgedata.Get(i).Candidate())
	{edgedata.Elem(i).SetConfirmed();}
    }
}

int STLEdgeDataList :: GetEdgeNum(int np1, int np2) const
{
  INDEX_2 ed(np1,np2);
  ed.Sort();
  if (hashtab.Used(ed))
    {
      return hashtab.Get(ed);
    }

//   int i;
//   for (i = 1; i <= Size(); i++)
//     {
//       if ((Get(i).p1 == np1 && Get(i).p2 == np2) ||
// 	  (Get(i).p2 == np1 && Get(i).p1 == np2))
// 	{
// 	  return i;
// 	}
//     }

  return 0;
}

const STLEdgeDataList& STLEdgeDataList :: operator=(const STLEdgeDataList& edl)
{
  int i;
  SetSize(edl.Size());
  for (i = 1; i <= Size(); i++)
    {
      Add(edl.Get(i), i);
    }
  return *this;
} 

void STLEdgeDataList :: Add(const STLEdgeData& ed, int i)
{
  INDEX_2 edge(ed.p1,ed.p2);
  edge.Sort();
  hashtab.Set(edge, i);
  Elem(i) = ed;
  AddEdgePP(ed.p1,i);
  AddEdgePP(ed.p2,i);
}

void STLEdgeDataList :: Write(ofstream& of) const
{
  of.precision(16);
  int i;
  of << Size() << endl;
  
  for (i = 1; i <= Size(); i++)
    {
      Get(i).Write(of);
    }
}

void STLEdgeDataList :: Read(ifstream& ifs)
{
  int i,n;
  ifs >> n;

  SetSize(n);
  STLEdgeData ed;
  for (i = 1; i <= n; i++)
    {
      ed.Read(ifs);
      Add(ed,i);
    }
}

int STLEdgeDataList :: GetNEPPStat(int p, int status) const
{
  int i;
  int cnt = 0;
  for (i = 1; i <= GetNEPP(p); i++)
    {
      if (Get(GetEdgePP(p,i)).GetStatus() == status)
	{
	  cnt++;
	}
    }
  return cnt;
}

int STLEdgeDataList :: GetNConfCandEPP(int p) const
{
  int i;
  int cnt = 0;
  for (i = 1; i <= GetNEPP(p); i++)
    {
      if (Get(GetEdgePP(p,i)).ConfCand())
	{
	  cnt++;
	}
    }
  return cnt;
}


void STLEdgeDataList :: BuildLineWithEdge(int ep1, int ep2, ARRAY<twoint>& line)
{
  int status = Get(GetEdgeNum(ep1,ep2)).GetStatus();

  int found, pstart, p, en, pnew, ennew;
  int closed = 0;
  int j, i;
  for (j = 1; j <= 2; j++)
    {
      if (j == 1) {p = ep1;}
      if (j == 2) {p = ep2;}

      pstart = p;
      en = GetEdgeNum(ep1,ep2);

      found = 1;
      while (found && !closed)
	{
	  found = 0;
	  
	  if (GetNEPPStat(p,status) == 2)
	    {
	      for (i = 1; i <= GetNEPP(p); i++)
		{		
		  const STLEdgeData& e = Get(GetEdgePP(p,i));
		  if (GetEdgePP(p,i) != en && e.GetStatus() == status) 
		    {
		      if (e.p1 == p) 
			{pnew = e.p2;}
		      else 
			{pnew = e.p1;}

		      ennew = GetEdgePP(p,i);
		    }
		}
	      if (pnew == pstart) {closed = 1;}
	      else
		{
		  line.Append(twoint(p,pnew));
		  p = pnew;
		  en = ennew;
		  found = 1;
		}
	    }
	}
    }
  
}
*/




STLEdgeDataList :: STLEdgeDataList (STLTopology & ageom)
  : geom(ageom)
{
  ;
}

STLEdgeDataList :: ~STLEdgeDataList()
{
  ;
}


void STLEdgeDataList :: Store ()
{
  int i, ne = geom.GetNTE();
  storedstatus.SetSize(ne);
  for (i = 1; i <= ne; i++)
    {
      storedstatus.Elem(i) = Get(i).GetStatus();
    }
}

void STLEdgeDataList :: Restore ()
{
  int i, ne = geom.GetNTE();
  if (storedstatus.Size() == ne)
    for (i = 1; i <= ne; i++)
      geom.GetTopEdge(i).SetStatus (storedstatus.Elem(i));
}


void STLEdgeDataList :: ResetAll()
{
  int i, ne = geom.GetNTE();
  for (i = 1; i <= ne; i++)
    geom.GetTopEdge (i).SetStatus (ED_UNDEFINED);
}

int STLEdgeDataList :: GetNConfEdges() const
{
  int i, ne = geom.GetNTE();
  int cnt = 0;
  for (i = 1; i <= ne; i++)
    if (geom.GetTopEdge (i).GetStatus() == ED_CONFIRMED)
      cnt++;
  return cnt; 
}

void STLEdgeDataList :: ChangeStatus(int status1, int status2)
{
  int i, ne = geom.GetNTE();
  for (i = 1; i <= ne; i++)
    if (geom.GetTopEdge (i).GetStatus() == status1)
      geom.GetTopEdge (i).SetStatus (status2);
}

/*
void STLEdgeDataList :: Add(const STLEdgeData& ed, int i)
{
  INDEX_2 edge(ed.p1,ed.p2);
  edge.Sort();
  hashtab.Set(edge, i);
  Elem(i) = ed;
  AddEdgePP(ed.p1,i);
  AddEdgePP(ed.p2,i);
}
*/

void STLEdgeDataList :: Write(ofstream& of) const
{
  
  /*
  of.precision(16);
  int i;
  of << Size() << endl;
  
  for (i = 1; i <= Size(); i++)
    {
      Get(i).Write(of);
    }

  */
  of.precision(16);
  int i, ne = geom.GetNTE();
  //of << GetNConfEdges() << endl;
  of << geom.GetNTE() << endl;

  for (i = 1; i <= ne; i++)
    {
      const STLTopEdge & edge = geom.GetTopEdge(i);
      //if (edge.GetStatus() == ED_CONFIRMED)
      of << edge.GetStatus() << " ";

      const Point3d & p1 = geom.GetPoint (edge.PNum(1));
      const Point3d & p2 = geom.GetPoint (edge.PNum(2));
      of << p1.X() << " "
	 << p1.Y() << " "
	 << p1.Z() << " "
	 << p2.X() << " "
	 << p2.Y() << " "
	 << p2.Z() << endl;
    }
  
}

void STLEdgeDataList :: Read(ifstream& ifs)
{
  int i, nce;
  Point3d p1, p2;
  int pi1, pi2;
  int status, ednum;

  ifs >> nce;
  for (i = 1; i <= nce; i++)
    {
      ifs >> status;
      ifs >> p1.X() >> p1.Y() >> p1.Z();
      ifs >> p2.X() >> p2.Y() >> p2.Z();

      pi1 = geom.GetPointNum (p1);
      pi2 = geom.GetPointNum (p2);
      ednum = geom.GetTopEdgeNum (pi1, pi2);


      if (ednum)
	{ 
	  geom.GetTopEdge(ednum).SetStatus (status);
	//	geom.GetTopEdge (ednum).SetStatus (ED_CONFIRMED);
	}
    }
    /*
  int i,n;
  ifs >> n;

  SetSize(n);
  STLEdgeData ed;
  for (i = 1; i <= n; i++)
    {
      ed.Read(ifs);
      Add(ed,i);
    }
  */
}

int STLEdgeDataList :: GetNEPPStat(int p, int status) const
{
  int i;
  int cnt = 0;
  for (i = 1; i <= GetNEPP(p); i++)
    {
      if (Get(GetEdgePP(p,i)).GetStatus() == status)
	{
	  cnt++;
	}
    }
  return cnt;
}

int STLEdgeDataList :: GetNConfCandEPP(int p) const
{
  int i;
  int cnt = 0;
  for (i = 1; i <= GetNEPP(p); i++)
    {
      if (Get(GetEdgePP(p,i)).GetStatus() == ED_CANDIDATE || 
	  Get(GetEdgePP(p,i)).GetStatus() == ED_CONFIRMED)
	{
	  cnt++;
	}
    }
  return cnt;
}


void STLEdgeDataList :: BuildLineWithEdge(int ep1, int ep2, ARRAY<twoint>& line)
{
  int status = Get(GetEdgeNum(ep1,ep2)).GetStatus();

  int found, pstart, p(0), en, pnew(0), ennew(0);
  int closed = 0;
  int j, i;
  for (j = 1; j <= 2; j++)
    {
      if (j == 1) {p = ep1;}
      if (j == 2) {p = ep2;}

      pstart = p;
      en = GetEdgeNum(ep1,ep2);

      found = 1;
      while (found && !closed)
	{
	  found = 0;
	  
	  if (GetNEPPStat(p,status) == 2)
	    {
	      for (i = 1; i <= GetNEPP(p); i++)
		{		
		  const STLTopEdge & e = Get(GetEdgePP(p,i));
		  if (GetEdgePP(p,i) != en && e.GetStatus() == status) 
		    {
		      if (e.PNum(1) == p) 
			{pnew = e.PNum(2);}
		      else 
			{pnew = e.PNum(1);}

		      ennew = GetEdgePP(p,i);
		    }
		}
	      if (pnew == pstart) {closed = 1;}
	      else
		{
		  line.Append(twoint(p,pnew));
		  p = pnew;
		  en = ennew;
		  found = 1;
		}
	    }
	}
    }
  
}

int Exists(int p1, int p2, const ARRAY<twoint>& line)
{
  int i;
  for (i = 1; i <= line.Size(); i++)
    {
      if (line.Get(i).i1 == p1 && line.Get(i).i2 == p2 ||
	  line.Get(i).i1 == p2 && line.Get(i).i2 == p1) {return 1;}
    }
  return 0;
}

void STLEdgeDataList :: BuildClusterWithEdge(int ep1, int ep2, ARRAY<twoint>& line)
{
  int status = Get(GetEdgeNum(ep1,ep2)).GetStatus();

  int p(0), en;
  int j, i, k;
  int oldend;
  int newend = 1;
  int pnew, ennew(0);

  int changed = 1;
  while (changed)
    {
      changed = 0;
      for (j = 1; j <= 2; j++)
	{
	  oldend = newend;
	  newend = line.Size();
	  for (k = oldend; k <= line.Size(); k++)
	    {
	      if (j == 1) p = line.Get(k).i1;
	      if (j == 2) p = line.Get(k).i2;
	      en = GetEdgeNum(line.Get(k).i1, line.Get(k).i2);

	      for (i = 1; i <= GetNEPP(p); i++)
		{		
		  pnew = 0;
		  const STLTopEdge & e = Get(GetEdgePP(p,i));
		  if (GetEdgePP(p,i) != en && e.GetStatus() == status) 
		    {
		      if (e.PNum(1) == p) 
			{pnew = e.PNum(2);}
		      else 
			{pnew = e.PNum(1);}

		      ennew = GetEdgePP(p,i);
		    }
		  if (pnew && !Exists(p,pnew,line))
		    {
		      changed = 1;
		      line.Append(twoint(p,pnew));
		      p = pnew;
		      en = ennew;
		    }
		}
	      
	    }
	}

    }

}










//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++   STL LINE    +++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

STLLine :: STLLine(const STLGeometry * ageometry)
  : pts(), lefttrigs(), righttrigs()
{
  geometry = ageometry;
  split = 0;
};

int STLLine :: GetNS() const
{
  if (pts.Size() <= 1) {return 0;}
  return pts.Size()-1;
}
void STLLine :: GetSeg(int nr, int& p1, int& p2) const
{
  p1 = pts.Get(nr);
  p2 = pts.Get(nr+1);
}

int STLLine :: GetLeftTrig(int nr) const 
{
  if (nr > lefttrigs.Size()) {PrintSysError("In STLLine::GetLeftTrig!!!"); return 0;}
  return lefttrigs.Get(nr);
};

int STLLine :: GetRightTrig(int nr) const 
{
  if (nr > righttrigs.Size()) {PrintSysError("In STLLine::GetRightTrig!!!"); return 0;}
  return righttrigs.Get(nr);
};

double STLLine :: GetSegLen(const ARRAY<Point<3> >& ap, int nr) const
{
  return Dist(ap.Get(PNum(nr)),ap.Get(PNum(nr+1)));
}

double STLLine :: GetLength(const ARRAY<Point<3> >& ap) const
{
  double len = 0;
  for (int i = 2; i <= pts.Size(); i++)
    {
      len += (ap.Get(pts.Get(i)) - ap.Get(pts.Get(i-1))).Length();
    }
  return len;
}

void STLLine :: GetBoundingBox (const ARRAY<Point<3> > & ap, Box<3> & box) const
{
  box.Set (ap.Get (pts[0]));
  for (int i = 1; i < pts.Size(); i++)
    box.Add (ap.Get(pts[i]));
}



Point<3> STLLine :: 
GetPointInDist(const ARRAY<Point<3> >& ap, double dist, int& index) const
{
  if (dist <= 0)
    {
      index = 1;
      return ap.Get(StartP());
    }
  
  double len = 0;
  int i;
  for (i = 1; i < pts.Size(); i++)
    {
      double seglen = Dist (ap.Get(pts.Get(i)),
			    ap.Get(pts.Get(i+1)));

      if (len + seglen > dist)
	{
	  index = i;
	  double relval = (dist - len) / (seglen + 1e-16);
	  Vec3d v (ap.Get(pts.Get(i)), ap.Get(pts.Get(i+1)));
	  return ap.Get(pts.Get(i)) + relval * v;
	}

      len += seglen;
    }

  index = pts.Size() - 1;
  return ap.Get(EndP());
}


/*
double stlgh;
double GetH(const Point3d& p, double x) 
{
  return stlgh;//+0.5)*(x+0.5);
}
*/
STLLine* STLLine :: Mesh(const ARRAY<Point<3> >& ap, 
			 ARRAY<Point3d>& mp, double ghi,
			 class Mesh& mesh) const
{
  STLLine* line = new STLLine(geometry);

  //stlgh = ghi; //uebergangsloesung!!!!
  
  double len = GetLength(ap);
  double inthl = 0; //integral of 1/h
  double dist = 0;
  double h;
  int ind;
  Point3d p;

  int i, j;

  Box<3> bbox;
  GetBoundingBox (ap, bbox);
  double diam = bbox.Diam();

  double minh = mesh.LocalHFunction().GetMinH (bbox.PMin(), bbox.PMax());

  double maxseglen = 0;
  for (i = 1; i <= GetNS(); i++)
    maxseglen = max2 (maxseglen, GetSegLen (ap, i));
  
  int nph = 10+int(maxseglen / minh); //anzahl der integralauswertungen pro segment

  ARRAY<double> inthi(GetNS()*nph);
  ARRAY<double> curvelen(GetNS()*nph);


  for (i = 1; i <= GetNS(); i++)
    {
      //double seglen = GetSegLen(ap,i);
      for (j = 1; j <= nph; j++)
	{
	  p = GetPointInDist(ap,dist,ind);
	  //h = GetH(p,dist/len);
	  h = mesh.GetH(p);

	  
	  dist += GetSegLen(ap,i)/(double)nph;
	  
	  inthl += GetSegLen(ap,i)/nph/(h);
	  inthi.Elem((i-1)*nph+j) = GetSegLen(ap,i)/nph/h;
	  curvelen.Elem((i-1)*nph+j) = GetSegLen(ap,i)/nph;
	}
    }


  int inthlint = int(inthl+1);

  if ( (inthlint < 3) && (StartP() == EndP()))
    {
      inthlint = 3;
    }
  if ( (inthlint == 1) && ShouldSplit())
    {
      inthlint = 2; 
    }
     
  double fact = inthl/(double)inthlint;
  dist = 0;
  j = 1;


  p = ap.Get(StartP());
  int pn = AddPointIfNotExists(mp, p, 1e-10*diam);

  int segn = 1;
  line->AddPoint(pn);
  line->AddLeftTrig(GetLeftTrig(segn));
  line->AddRightTrig(GetRightTrig(segn));
  line->AddDist(dist);

  inthl = 0; //restart each meshseg
  for (i = 1; i <= inthlint; i++)
    {
      while (inthl < 1.000000001 && j <= inthi.Size())
      //      while (inthl-1. < 1e-9) && j <= inthi.Size())
	{
	  inthl += inthi.Get(j)/fact;
	  dist += curvelen.Get(j);
	  j++;
	}

      //went to far:
      j--;
      double tofar = (inthl - 1)/inthi.Get(j);
      inthl -= tofar*inthi.Get(j);
      dist -= tofar*curvelen.Get(j)*fact;

      if (i == inthlint && fabs(dist - len) >= 1E-8) 
	{
	  PrintSysError("meshline failed!!!"); 
	}

      if (i != inthlint) 
	{
	  p = GetPointInDist(ap,dist,ind);
	  pn = AddPointIfNotExists(mp, p, 1e-10*diam);
	  segn = ind;
	  line->AddPoint(pn);
	  line->AddLeftTrig(GetLeftTrig(segn));
	  line->AddRightTrig(GetRightTrig(segn));
	  line->AddDist(dist);
	}

      inthl = tofar*inthi.Get(j);
      dist += tofar*curvelen.Get(j)*fact;
      j++;
    }

  p = ap.Get(EndP());
  pn = AddPointIfNotExists(mp, p, 1e-10*diam);
  segn = GetNS();
  line->AddPoint(pn);
  line->AddLeftTrig(GetLeftTrig(segn));
  line->AddRightTrig(GetRightTrig(segn));
  line->AddDist(dist);
  
  for (int ii = 1; ii <= line->GetNS(); ii++)
    {
      int p1, p2;
      line->GetSeg(ii,p1,p2);
    }
  /*  
  (*testout) << "line, " << ap.Get(StartP()) << "-" << ap.Get(EndP())
	     << " len = " << Dist (ap.Get(StartP()), ap.Get(EndP())) << endl;
  */
  return line;
}
}
