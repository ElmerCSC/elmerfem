//20.11.1999 second part of stlgeom.cc, mainly mesh functions

#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <gprim.hpp>

#include <meshing.hpp>

#include "stlgeom.hpp"

namespace netgen
{
int EdgeUsed(int p1, int p2, ARRAY<INDEX_2>& edges, INDEX_2_HASHTABLE<int>& hashtab)
{
  if (p1 > p2) {swap (p1,p2);}

  if (hashtab.Used(INDEX_2(p1,p2))) 
    {return hashtab.Get(INDEX_2(p1,p2));}

  return 0;
}

Point<3> STLGeometry :: PointBetween(const Point<3> & ap1, int t1, 
				     const Point<3> & ap2, int t2)
{
  //funktioniert nicht in allen Fällen!

  PrintWarning("Point between");


  ClearMarkedSegs();

  InitMarkedTrigs();
  SetMarkedTrig(t1,1);
  SetMarkedTrig(t2,1);

  TABLE<Point3d> edgepoints;
  TABLE<double> edgepointdists;
  TABLE<int> edgepointorigines;
  TABLE<int> edgepointoriginps;

  ARRAY<int> edgetrigs;
  ARRAY<INDEX_2> edgepointnums;
  ARRAY<int> edgetriglocinds;

  int size = 3*GetNT();
  INDEX_2_HASHTABLE<int> hashtab(size);

  int divisions = 10;

  edgepoints.SetSize(size);
  edgepointdists.SetSize(size);
  edgepointorigines.SetSize(size);
  edgepointoriginps.SetSize(size);

  edgetrigs.SetSize(size);
  edgepointnums.SetSize(size);
  edgetriglocinds.SetSize(size);

  ARRAY<int> edgelist1;
  ARRAY<int> edgelist2;

  edgelist1.SetSize(0);
  edgelist2.SetSize(0);


  int i, j, k, l, m;
  int edgecnt = 0;

  //first triangle:
  for (i = 1; i <= 3; i++)
    {
      int ptn1 = GetTriangle(t1).PNum(i);
      int ptn2 = GetTriangle(t1).PNumMod(i+1);

      if (ptn1 > ptn2) {swap(ptn1,ptn2);}

      Point3d pt1 = GetPoint(ptn1);
      Point3d pt2 = GetPoint(ptn2);

      edgecnt++;
      edgetrigs.Elem(edgecnt) = t1;
      edgepointnums.Elem(edgecnt) = INDEX_2(ptn1,ptn2);
      hashtab.Set(edgepointnums.Get(edgecnt),edgecnt);

      edgetriglocinds.Elem(edgecnt) = i;
      edgelist1.Append(edgecnt);

      for (j = 1; j <= divisions; j++)
	{
	  double lfact = (double)j/(double)divisions;
	  Point3d pbtw(lfact*pt1.X()+(1.-lfact)*pt2.X(),
		       lfact*pt1.Y()+(1.-lfact)*pt2.Y(),
		       lfact*pt1.Z()+(1.-lfact)*pt2.Z());

	  //AddMarkedSeg(ap1,pbtw);
	
	  edgepoints.Add1(edgecnt,pbtw);
	  edgepointdists.Add1(edgecnt,Dist(pbtw,ap1));
	  edgepointorigines.Add1(edgecnt,0);
	  edgepointoriginps.Add1(edgecnt,0);
	}
    }

  int finished = 0;
  int endpointorigine = 0;
  int endpointoriginp = 0;
  double endpointmindist = 1E50;

  int maxsize = 0;
  while (!finished)
    {
      finished = 1;
      
      if (edgelist1.Size() > maxsize) {maxsize = edgelist1.Size();}

      for (i = 1; i <= edgelist1.Size(); i++)
	{
	  int en = edgelist1.Get(i);
	  int trig = edgetrigs.Get(en);
	  int edgenum = edgetriglocinds.Get(en);
	  int tn = NeighbourTrigSorted(trig,edgenum);

	  if (tn != t2)
	    {
	      for (k = 1; k <= 3; k++)
		{
		  int pnt1 = GetTriangle(tn).PNum(k);
		  int pnt2 = GetTriangle(tn).PNumMod(k+1);
		      
		  if (pnt1 > pnt2) {swap(pnt1,pnt2);}

		  Point3d pt1 = GetPoint(pnt1);
		  Point3d pt2 = GetPoint(pnt2);
		      
		  //AddMarkedSeg(pt1,pt2);
		  
		  //if (!(pnt1 == ep1 && pnt2 == ep2))
		  //  {
		  int edgeused = 0;
		  edgenum = EdgeUsed(pnt1, pnt2, edgepointnums, hashtab);
		  if (edgenum != en)
		    {
		      if (edgenum != 0) 
			{edgeused = 1;}
		      else 
			{
			  edgecnt++; 
			  edgenum = edgecnt;
			  
			  edgetrigs.Elem(edgenum) = tn;
			  edgepointnums.Elem(edgenum) = INDEX_2(pnt1,pnt2);
			  hashtab.Set(edgepointnums.Get(edgenum),edgenum);
			  edgetriglocinds.Elem(edgenum) = k;
			}
		      
		      if (edgenum > size || edgenum == 0) {PrintSysError("edgenum = ", edgenum);}
			  
		      double minofmindist = 1E50;
		      int changed = 0;
		      
		      for (l = 1; l <= divisions; l++)
			{
			  double lfact = (double)l/(double)divisions;
			  Point3d pbtw(lfact*pt1.X()+(1.-lfact)*pt2.X(),
				       lfact*pt1.Y()+(1.-lfact)*pt2.Y(),
				       lfact*pt1.Z()+(1.-lfact)*pt2.Z());
			  
			  double mindist = 1E50;
			  int index=0;
			  
			  for (m = 1; m <= divisions; m++)
			    {
			      const Point3d& p = edgepoints.Get(en,m);
			      if (Dist(pbtw,p) + edgepointdists.Get(en,m) < mindist)
				{mindist = Dist(pbtw,p) + edgepointdists.Get(en,m); index = m;}
			    }
			  
			  //if (mindist < endpointmindist) {finished = 0;}
			  if (mindist < minofmindist) {minofmindist = mindist;}
			  
			  
			  if (!edgeused)
			    {
			      //AddMarkedSeg(pbtw,edgepoints.Get(en,index));

			      edgepoints.Add1(edgenum,pbtw);
			      edgepointdists.Add1(edgenum,mindist);
			      edgepointorigines.Add1(edgenum,en);
			      edgepointoriginps.Add1(edgenum,index);
			      changed = 1;
			    }
			  else
			    {
			      if (mindist < edgepointdists.Get(edgenum,l))
				{
				  edgepointdists.Set(edgenum,l,mindist);
				  edgepointorigines.Set(edgenum,l,en);
				  edgepointoriginps.Set(edgenum,l,index);
				  changed = 1;
				}			      
			    }
			}
		      if (minofmindist < endpointmindist-1E-10 && changed)
			{
			  finished = 0;
			  edgelist2.Append(edgenum);
			}
		    }
		}
	    }
	  else
	    {
	      double mindist = 1E50;
	      int index(0);
	      for (m = 1; m <= divisions; m++)
		{
		  const Point3d& p = edgepoints.Get(en,m);
		  if (Dist(ap2,p) + edgepointdists.Get(en,m) < mindist)
		    {mindist = Dist(ap2,p) + edgepointdists.Get(en,m); index = m;}
		}
	      if (mindist < endpointmindist)
		{
		  endpointorigine = en;
		  endpointoriginp = index;
		  endpointmindist = mindist;
		}
	    }
	}
      edgelist1.SetSize(0);
      for (i = 1; i <= edgelist2.Size(); i++)
	{
	  edgelist1.Append(edgelist2.Get(i));
	}
    }

  if (!endpointorigine) {PrintSysError("No connection found!");}

  ARRAY<Point3d> plist;

  plist.Append(ap2);
  int laste = endpointorigine;
  int lastp = endpointoriginp;
  int lle, llp;


  while (laste)
    {
      plist.Append(edgepoints.Get(laste,lastp));

      lle = laste;
      llp = lastp; 
      laste = edgepointorigines.Get(lle,llp);
      lastp = edgepointoriginps.Get(lle,llp);
    }

  plist.Append(ap1);

  for (i = 1; i <= plist.Size()-1; i++)
    {
      AddMarkedSeg(plist.Get(i),plist.Get(i+1));
    }

  PrintMessage(5,"PointBetween: complexity=", maxsize);


  Point3d pm;
  double dist = 0;
  int found = 0;
  
  for (i = 1; i <= plist.Size()-1; i++)
    {
      dist += Dist(plist.Get(i),plist.Get(i+1));
      if (dist > endpointmindist*0.5) 
	{
	  double segl = Dist(plist.Get(i), plist.Get(i+1));
	  double d = dist - endpointmindist * 0.5;
	  pm = Point3d(d/segl*plist.Get(i).X() + (1.-d/segl)*plist.Get(i+1).X(),
		       d/segl*plist.Get(i).Y() + (1.-d/segl)*plist.Get(i+1).Y(),
		       d/segl*plist.Get(i).Z() + (1.-d/segl)*plist.Get(i+1).Z());
	  found = 1;
	  break;
	}
    }
  if (!found) {PrintWarning("Problem in PointBetween"); pm = Center(ap1,ap2);}

  AddMarkedSeg(pm, Point3d(0.,0.,0.));
  
  return pm;
  
}


void STLGeometry :: PrepareSurfaceMeshing()
{
  meshchart = -1; //clear no old chart
  meshcharttrigs.SetSize(GetNT());
  int i;
  for (i = 1; i <= GetNT(); i++) 
    {meshcharttrigs.Elem(i) = 0;}
}

void STLGeometry::GetMeshChartBoundary (ARRAY<Point2d > & apoints,
					ARRAY<Point3d > & points3d,
					ARRAY<INDEX_2> & alines, double h)
{
  int i, j;
  twoint seg, newseg;
  int zone;
  Point<2> p2;

  const STLChart& chart = GetChart(meshchart);


  for (i = 1; i <= chart.GetNOLimit(); i++)
    {
      seg = chart.GetOLimit(i);
      INDEX_2 i2;
      for (j = 1; j <= 2; j++)
	{
	  int pi = (j == 1) ? seg.i1 : seg.i2;
	  int lpi;
	  if (ha_points.Get(pi) == 0)
	    {
	      const Point<3> & p3d = GetPoint (pi);
	      Point<2> p2d;

	      points3d.Append (p3d);
	      ToPlane(p3d, 0, p2d, h, zone, 0);
	      apoints.Append (p2d);
	      
	      lpi = apoints.Size();
	      ha_points.Elem(pi) = lpi;
	    }
	  else
	    lpi = ha_points.Get(pi);

	  i2.I(j) = lpi;
	}
      alines.Append (i2);

      /*
      seg = chart.GetOLimit(i);
      psize = points.Size();

      newseg.i1 = psize+1;
      newseg.i2 = psize+2;

      ToPlane(GetPoint(seg.i1), 0, p2, h, zone, 0);
      points.Append(p2);
      points3d.Append (GetPoint(seg.i1));
      ToPlane(GetPoint(seg.i2), 0, p2, h, zone, 0);
      points.Append(p2);
      points3d.Append (GetPoint(seg.i2));
      lines.Append (INDEX_2 (points.Size()-1, points.Size()));
      */
    }

  for (i = 1; i <= chart.GetNOLimit(); i++)
    {
      seg = chart.GetOLimit(i);
      ha_points.Elem(seg.i1) = 0;
      ha_points.Elem(seg.i2) = 0;
    }
}

void STLGeometry :: DefineTangentialPlane (const Point<3> & ap1, const Point<3> & ap2, int trig)
{
  p1 = ap1; //save for ToPlane, in data of STLGeometry class
  Point<3> p2 = ap2; //only locally used

  meshchart = GetChartNr(trig);

  if (usechartnormal)
    meshtrignv = GetChart(meshchart).GetNormal();
  else
    meshtrignv = GetTriangle(trig).Normal();

  //meshtrignv = GetTriangle(trig).Normal(points);

  meshtrignv /= meshtrignv.Length();

  GetTriangle(trig).ProjectInPlain(points, meshtrignv, p2);


  ez = meshtrignv;
  ez /= ez.Length();
  ex = p2 - p1;
  ex -= (ex * ez) * ez;
  ex /= ex.Length();
  ey = Cross (ez, ex);

}


void STLGeometry :: SelectChartOfTriangle (int trignum)
{
  meshchart = GetChartNr(trignum);
  meshtrignv = GetTriangle(trignum).Normal();	
}


void STLGeometry :: SelectChartOfPoint (const Point<3> & p)
{
  int i, ii;

  ARRAY<int> trigsinbox;
  
  Box<3> box(p,p);
  box.Increase (1e-6);
  GetTrianglesInBox (box, trigsinbox);
  

  //  for (i = 1; i <= GetNT(); i++)
  for (ii = 1; ii <= trigsinbox.Size(); ii++)
    {
      i = trigsinbox.Get(ii);
      Point<3> hp = p;
      if (GetTriangle(i).GetNearestPoint(points, hp) <= 1E-8)
	{
	  SelectChartOfTriangle (i);
	  break;
      }
    }
  return;
}



void STLGeometry :: ToPlane (const Point<3> & locpoint, int * trigs,
			     Point<2> & plainpoint, double h, int& zone,
			     int checkchart)
{
  if (checkchart)
    {

      //check if locpoint lies on actual chart:
      zone = 0;
      
      
      //  Point3d p;
      int i = 1;
      const STLChart& chart = GetChart(meshchart);
      int foundinchart = 0;
      const double range = 1e-6; //1e-4 old
      
      
      
      
      if (trigs)
	{
	  int * htrigs = trigs;
	  while (*htrigs)
	    {
	      if (TrigIsInOC (*htrigs, meshchart))
		{
		  foundinchart = 1;
		  break;
		}
	      htrigs++;
	    }
	}
      
      else
	{
	  ARRAY<int> trigsinbox;

	  if (!geomsearchtreeon)
	    {
	      //alter chart-tree
	      Box<3> box(locpoint, locpoint);
	      box.Increase (range);
	      chart.GetTrianglesInBox (box.PMin(), box.PMax(), trigsinbox);
	    }
	  else
	    {
	      ARRAY<int> trigsinbox2;
	      Box<3> box(locpoint, locpoint);
	      box.Increase (range);
	      GetTrianglesInBox (box, trigsinbox2);
	      for (i = 1; i <= trigsinbox2.Size(); i++)
		{
		  if (TrigIsInOC(trigsinbox2.Get(i),meshchart)) {trigsinbox.Append(trigsinbox2.Get(i));}
		}
	      
	    }
	  
	  
	  for (i = 1; i <= trigsinbox.Size(); i++)
	    {
	      Point<3> p = locpoint;
	      if (GetTriangle(trigsinbox.Get(i)).GetNearestPoint(points, p) 
		  <= 1E-8)
		{
		  foundinchart = 1;
		  break;
		}
	      
	    }
	}
      
  //do not use this point (but do correct projection (joachim)
      if (!foundinchart) 
	{
	  zone = -1; // plainpoint.X() = 11111; plainpoint.Y() = 11111; return; 
	}
    }
  
  else
    {
      zone = 0;
    }
  
  //transform in plane
  Vec<3> p1p = locpoint - p1;
  plainpoint(0) = (p1p * ex) / h;
  plainpoint(1) = (p1p * ey) / h;

}

int STLGeometry :: FromPlane (const Point<2> & plainpoint, 
			      Point<3> & locpoint, double h)
{
  Point2d plainpoint2 (plainpoint);

  plainpoint2.X() *= h;
  plainpoint2.Y() *= h;
  Vec3d p1p = plainpoint2.X() * ex + plainpoint2.Y() * ey;
  locpoint = p1 + p1p;


  int rv = Project(locpoint);
  if (!rv) {return 1;} //project nicht gegangen
  return 0;
}

int lasttrig;
int STLGeometry :: LastTrig() const {return lasttrig;};

//project normal to tangential plane
int STLGeometry :: Project(Point<3> & p3d) const
{
  Point<3> p, pf;

  int i, j;
  int fi = 0;
  int cnt = 0;
  int different = 0;
  const double lamtol = 1e-6;

  const STLChart& chart = GetChart(meshchart);

  int nt = chart.GetNT();

   QuadraticFunction3d quadfun(p3d, meshtrignv);
 
   /*
     Vec3d hv = meshtrignv;
     hv /= hv.Length();
     Vec3d t1, t2;
     hv.GetNormal (t1);
     Cross (hv, t1, t2);
   */
  
  for (j = 1; j <= nt; j++)
    {
      i = chart.GetTrig(j);

      const Point<3> & c = GetTriangle(i).center;
      /*
      double d1 = t1 * (c-p3d);
      double d2 = t2 * (c-p3d);
      */
      /*
      if (d1 * d1 + d2 * d2 > sqr (GetTriangle(i).rad))
	continue;
      */
      if (quadfun.Eval(c) > sqr (GetTriangle(i).rad))
	continue;

      p = p3d;
      Vec<3> lam;
      int err = GetTriangle(i).ProjectInPlain(points, meshtrignv, p, lam);      
      int inside = (err == 0 && lam(0) > -lamtol && 
		    lam(1) > -lamtol && (1-lam(0)-lam(1)) > -lamtol);


      /*
      p = p3d;
      GetTriangle(i).ProjectInPlain(points, meshtrignv, p);
      if (GetTriangle(i).PointInside(points, p)) 
      */
      if (inside)
	{
	  if (cnt != 0) 
	    {
	      if (Dist2(p,pf)>=1E-16) 
		{
		  //		  (*testout) << "ERROR: found two points to project which are different" << endl;
		  //(*testout) << "p=" << p << ", pf=" << pf << endl;
		  different = 1;
		}
	    }
	  pf = p; fi = i; cnt++;
	}

      if (inside)
	break;

    }

  //  if (cnt == 2) {(*testout) << "WARNING: found 2 triangles to project" << endl;}
  //if (cnt == 3) {(*testout) << "WARNING: found 3 triangles to project" << endl;}
  //if (cnt > 3) {(*testout) << "WARNING: found more than 3 triangles to project" << endl;}

  if (fi != 0) {lasttrig = fi;}
  if (fi != 0 && !different) {p3d = pf; return fi;}

  //  (*testout) << "WARNING: Project failed" << endl;
  return 0;
  
}

//project normal to tangential plane
int STLGeometry :: ProjectOnWholeSurface(Point<3> & p3d) const
{
  Point<3> p, pf;

  int i;
  int fi = 0;
  int cnt = 0;
  int different = 0;
  const double lamtol = 1e-6;

  for (i = 1; i <= GetNT(); i++)
    {
      p = p3d;
      Vec<3> lam;
      int err =
	GetTriangle(i).ProjectInPlain(points, meshtrignv, p, lam);      
      int inside = (err == 0 && lam(0) > -lamtol && 
		    lam(1) > -lamtol && (1-lam(0)-lam(1)) > -lamtol);

      /*
      p = p3d;
      GetTriangle(i).ProjectInPlain(points, meshtrignv, p);
      if (GetTriangle(i).PointInside(points, p)) 
      */
      if (inside)
	{
	  if (cnt != 0) 
	    {
	      if (Dist2(p,pf)>=1E-16) 
		{
		  //		  (*testout) << "ERROR: found two points to project which are different" << endl;
		  //		  (*testout) << "p=" << p << ", pf=" << pf << endl;
		  different = 1;
		}
	    }
	  pf = p; fi = i; cnt++;
	}
    }
  /*
  if (cnt == 2) {(*testout) << "WARNING: found 2 triangles to project" << endl;}
  if (cnt == 3) {(*testout) << "WARNING: found 3 triangles to project" << endl;}
  if (cnt > 3) {(*testout) << "WARNING: found more than 3 triangles to project" << endl;}
  */
  if (fi != 0) {lasttrig = fi;}
  if (fi != 0 && !different) {p3d = pf; return fi;}

  //  (*testout) << "WARNING: Project failed" << endl;
  return 0;
  
}


int STLGeometry :: ProjectNearest(Point<3> & p3d) const
{
  Point<3> p, pf;

  //set new chart
  const STLChart& chart = GetChart(meshchart);
  int i;
  double nearest = 1E50;
  double dist;
  int ft = 0;

  for (i = 1; i <= chart.GetNT(); i++)
    {
      p = p3d;
      dist  = GetTriangle(chart.GetTrig(i)).GetNearestPoint(points, p);
      if (dist < nearest)
	{
	  pf = p;
	  nearest = dist;
	  ft = chart.GetTrig(i);
	}      
    }
  p3d = pf;
  //if (!ft) {(*testout) << "ERROR: ProjectNearest failed" << endl;}
  
  return ft;
}



	
//Restrict local h due to curvature for make atlas
void STLGeometry :: RestrictLocalHCurv(class Mesh & mesh, double gh)
{
  PushStatusF("Restrict H due to surface curvature");

  //bei jedem Dreieck alle Nachbardreiecke vergleichen, und, fallskein Kante dazwischen,
  //die Meshsize auf ein bestimmtes Mass limitieren
  int i,j;

  int ap1,ap2,p3,p4;
  Point<3> p1p, p2p, p3p, p4p;
  Vec<3> n, ntn;
  double rzyl, localh;

  //  double localhfact = 0.5;
  // double geometryignorelength = 1E-4;
  double minlocalh = stlparam.atlasminh;

  Box<3> bb = GetBoundingBox();
  //  mesh.SetLocalH(bb.PMin() - Vec3d(10, 10, 10),bb.PMax() + Vec3d(10, 10, 10),
  //		 mparam.grading);

  //  mesh.SetGlobalH(gh);

  double mincalch = 1E10;
  double maxcalch = -1E10;

  double objectsize = bb.Diam();
  double geometryignoreedgelength = objectsize * 1e-5;


  if (stlparam.resthatlasenable)
    {
      ARRAY<double> minh; //minimales h pro punkt
      minh.SetSize(GetNP());
      for (i = 1; i <= GetNP(); i++)
	{
	  minh.Elem(i) = gh;
	}
      
      for (i = 1; i <= GetNT(); i++)
	{
	  SetThreadPercent((double)i/(double)GetNT()*100.);

	  if (multithread.terminate)
	    {PopStatus(); return;}

	  const STLTriangle& trig = GetTriangle(i);
	  n = GetTriangle(i).Normal();
	  for (j = 1; j <= 3; j++)
	    {
	      const STLTriangle& nt = GetTriangle(NeighbourTrig(i,j));
	      
	      trig.GetNeighbourPointsAndOpposite(nt,ap1,ap2,p3);	    	    
	      
	      //checken, ob ap1-ap2 eine Kante sind
	      if (IsEdge(ap1,ap2)) continue;
	      
	      p4 = trig.PNum(1) + trig.PNum(2) + trig.PNum(3) - ap1 - ap2;
	      
	      p1p = GetPoint(ap1); p2p = GetPoint(ap2); 
	      p3p = GetPoint(p3); p4p = GetPoint(p4);
	      
	      double h1 = GetDistFromInfiniteLine(p1p,p2p, p4p);
	      double h2 = GetDistFromInfiniteLine(p1p,p2p, p3p);
	      double diaglen = Dist (p1p, p2p);
	      
	      if (diaglen < geometryignoreedgelength)
		continue;
	      rzyl = ComputeCylinderRadius 
		(n, GetTriangle(NeighbourTrig(i,j)).Normal(), 
		 h1, h2);
	      
	      
	      if (h1 < 1e-3 * diaglen && h2 < 1e-3 * diaglen)
		continue;
	      if (h1 < 1e-5 * objectsize && h2 < 1e-5 * objectsize)
		continue;
	      
	      
	      //	      rzyl = mindist/(2*sinang);
	      localh = 10.*rzyl / stlparam.resthatlasfac;
	      if (localh < mincalch) {mincalch = localh;}
	      if (localh > maxcalch) {maxcalch = localh;}

	      if (localh < minlocalh) {localh = minlocalh;}
	      if (localh < gh)
		{
		  minh.Elem(ap1) = min2(minh.Elem(ap1),localh);
		  minh.Elem(ap2) = min2(minh.Elem(ap2),localh);
		}
	      
	      mesh.RestrictLocalHLine(p1p, p2p, localh);
	    }
	  
	}
    }
  PrintMessage(5, "done\nATLAS H: nmin local h=", mincalch);
  PrintMessage(5, "ATLAS H: max local h=", maxcalch);
  PrintMessage(5, "Local h tree has ", mesh.LocalHFunction().GetNBoxes(), " boxes of size ",
	       (int)sizeof(GradingBox));

  PopStatus();

}
  //restrict local h due to near edges and due to outer chart distance
void STLGeometry :: RestrictLocalH(class Mesh & mesh, double gh)
{
  
  //bei jedem Dreieck alle Nachbardreiecke vergleichen, und, fallskein Kante dazwischen,
  //die Meshsize auf ein bestimmtes Mass limitieren
  int i,j;

  int ap1,ap2,p3,p4;
  Point3d p1p, p2p, p3p, p4p;
  Vec3d n, ntn;
  double rzyl, localh;

  //  double localhfact = 0.5;
  // double geometryignorelength = 1E-4;

  Box<3> bb = GetBoundingBox();
  //mesh.SetLocalH(bb.PMin() - Vec3d(10, 10, 10),bb.PMax() + Vec3d(10, 10, 10),
  //		 mparam.grading);

  //mesh.SetGlobalH(gh);

  double mincalch = 1E10;
  double maxcalch = -1E10;

  double objectsize = bb.Diam();
  double geometryignoreedgelength = objectsize * 1e-5;

  if (stlparam.resthsurfcurvenable)
    {
      PushStatusF("Restrict H due to surface curvature");

      ARRAY<double> minh; //minimales h pro punkt
      minh.SetSize(GetNP());
      for (i = 1; i <= GetNP(); i++)
	{
	  minh.Elem(i) = gh;
	}

      for (i = 1; i <= GetNT(); i++)
	{
	  SetThreadPercent((double)i/(double)GetNT()*100.);
	  if (i%20000==19999) {PrintMessage(7, (double)i/(double)GetNT()*100. , "%");}

	  if (multithread.terminate)
	    {PopStatus(); return;}
	  
	  const STLTriangle& trig = GetTriangle(i);
	  n = GetTriangle(i).Normal();
	  for (j = 1; j <= 3; j++)
	    {
	      const STLTriangle& nt = GetTriangle(NeighbourTrig(i,j));
	      
	      trig.GetNeighbourPointsAndOpposite(nt,ap1,ap2,p3);	    	    
	      
	      //checken, ob ap1-ap2 eine Kante sind
	      if (IsEdge(ap1,ap2)) continue;
	      
	      p4 = trig.PNum(1) + trig.PNum(2) + trig.PNum(3) - ap1 - ap2;
	      
	      p1p = GetPoint(ap1); p2p = GetPoint(ap2); 
	      p3p = GetPoint(p3); p4p = GetPoint(p4);
	      
	      double h1 = GetDistFromInfiniteLine(p1p,p2p, p4p);
	      double h2 = GetDistFromInfiniteLine(p1p,p2p, p3p);
	      double diaglen = Dist (p1p, p2p);
	      
	      if (diaglen < geometryignoreedgelength)
		continue;
	      rzyl = ComputeCylinderRadius 
		(n, GetTriangle (NeighbourTrig(i,j)).Normal(), 
		 h1, h2);
	      
	      
	      if (h1 < 1e-3 * diaglen && h2 < 1e-3 * diaglen)
		continue;
	      
	      if (h1 < 1e-5 * objectsize && h2 < 1e-5 * objectsize)
		continue;
	      
	      
	      //	      rzyl = mindist/(2*sinang);
	      localh = rzyl / stlparam.resthsurfcurvfac;
	      if (localh < mincalch) {mincalch = localh;}
	      if (localh > maxcalch) {maxcalch = localh;}
	      if (localh < gh) 
		{
		  minh.Elem(ap1) = min2(minh.Elem(ap1),localh);
		  minh.Elem(ap2) = min2(minh.Elem(ap2),localh);
		}
	      
	      //if (localh < 0.2) {localh = 0.2;}

	      if(localh < objectsize)
		mesh.RestrictLocalHLine(p1p, p2p, localh);
	      (*testout) << "restrict h along " << p1p << " - " << p2p << " to " << localh << endl;
	      
	      if (localh < 0.1)
		{
		  localh = 0.1;
		}
	      
	    }
	}
      PrintMessage(7, "done\nmin local h=", mincalch, "\nmax local h=", maxcalch);
      PopStatus();
    }

  if (stlparam.resthcloseedgeenable)
    {
      PushStatusF("Restrict H due to close edges");
      //geht nicht für spiralen!!!!!!!!!!!!!!!!!!
      
      double disttohfact = sqr(10.0 / stlparam.resthcloseedgefac);
      int k,l;
      double h1, h2, dist;
      int rc = 0;
      Point3d p3p1;
      double mindist = 1E50;
      
      PrintMessage(7,"build search tree...");
      Box3dTree* lsearchtree = new Box3dTree (GetBoundingBox().PMin() - Vec3d(1,1,1),
					     GetBoundingBox().PMax() + Vec3d(1,1,1));
      
      ARRAY<Point3d> pmins(GetNLines());
      ARRAY<Point3d> pmaxs(GetNLines());

      double maxhline;
      for (i = 1; i <= GetNLines(); i++)
	{
	  maxhline = 0;
	  STLLine* l1 = GetLine(i);
	  Point3d pmin(GetPoint(l1->StartP())), pmax(GetPoint(l1->StartP())), px;

	  for (j = 2; j <= l1->NP(); j++)
	    {
	      px = GetPoint(l1->PNum(j));
	      maxhline = max2(maxhline,mesh.GetH(px));
	      pmin.SetToMin (px);
	      pmax.SetToMax (px);
	    }
	  Box3d box(pmin,pmax);
	  box.Increase(maxhline);

	  lsearchtree->Insert (box.PMin(), box.PMax(), i);
	  pmins.Elem(i) = box.PMin();
	  pmaxs.Elem(i) = box.PMax();
	}

      ARRAY<int> linenums;
      int k2;

      for (i = 1; i <= GetNLines(); i++)
	{
	  SetThreadPercent((double)i/(double)GetNLines()*100.);
	  if (multithread.terminate)
	    {PopStatus(); return;}

	  linenums.SetSize(0);
	  lsearchtree->GetIntersecting(pmins.Get(i),pmaxs.Get(i),linenums);
	      
	  STLLine* l1 = GetLine(i);
	  for (j = 1; j <= l1->NP(); j++)
	    {
	      p3p1 = GetPoint(l1->PNum(j));
	      h1 = sqr(mesh.GetH(p3p1));
	      
	      for (k2 = 1; k2 <= linenums.Size(); k2++)
		{
		  k = linenums.Get(k2);
		  if (k <= i) {continue;} 
		  /*  
		   //old, without searchtrees
		     for (k = i+1; k <= GetNLines(); k++)
		     {
		  */
		  STLLine* l2 = GetLine(k);
		  for (l = 1; l <= l2->NP(); l++)
		    {
		      const Point3d& p3p2 = GetPoint(l2->PNum(l));
		      h2 = sqr(mesh.GetH(p3p2));
		      dist = Dist2(p3p1,p3p2)*disttohfact;		  
		      if (dist > 1E-12)
			{
			  if (dist < h1) 
			    {
			      mesh.RestrictLocalH(p3p1,sqrt(dist)); 
			      rc++;
			      mindist = min2(mindist,sqrt(dist));
			    }
			  if (dist < h2) 
			    {
			      mesh.RestrictLocalH(p3p2,sqrt(dist)); 
			      rc++;
			      mindist = min2(mindist,sqrt(dist));
			    }
			}
		    }
		}	  
	    }
	}
      PrintMessage(5, "done\n Restricted h in ", rc, " points due to near edges!");
      PopStatus(); 
    }

  if (stlparam.resthedgeangleenable)
    {
      PushStatusF("Restrict h due to close edges");

      int lp1, lp2;
      Vec3d v1,v2;
      mincalch = 1E50;
      maxcalch = -1E50;

      for (i = 1; i <= GetNP(); i++)
	{
	  SetThreadPercent((double)i/(double)GetNP()*100.);
	  if (multithread.terminate)
	    {PopStatus(); return;}

	  if (GetNEPP(i) == 2 && !IsLineEndPoint(i))
	    {
	      if (GetEdge(GetEdgePP(i,1)).PNum(2) == GetEdge(GetEdgePP(i,2)).PNum(1) ||
		  GetEdge(GetEdgePP(i,1)).PNum(1) == GetEdge(GetEdgePP(i,2)).PNum(2))
		{
		  lp1 = 1; lp2 = 2;
		}
	      else
		{
		  lp1 = 2; lp2 = 1;
		}

	      v1 = Vec3d(GetPoint(GetEdge(GetEdgePP(i,1)).PNum(1)),
			 GetPoint(GetEdge(GetEdgePP(i,1)).PNum(2)));
	      v2 = Vec3d(GetPoint(GetEdge(GetEdgePP(i,2)).PNum(lp1)),
			 GetPoint(GetEdge(GetEdgePP(i,2)).PNum(lp2)));

	      rzyl = ComputeCylinderRadius(v1, v2, v1.Length(), v2.Length());
	      	      
	      localh = rzyl / stlparam.resthedgeanglefac;
	      if (localh < mincalch) {mincalch = localh;}
	      if (localh > maxcalch) {maxcalch = localh;}
	      
	      if (localh != 0)
		mesh.RestrictLocalH(GetPoint(i), localh);
	    }	  
	}
      PrintMessage(7,"edge-angle min local h=", mincalch, "\nedge-angle max local h=", maxcalch);
      PopStatus();
    }

  if (stlparam.resthchartdistenable)
    {
      PushStatusF("Restrict H due to outer chart distance");
      
      // mesh.LocalHFunction().Delete();

      //berechne minimale distanz von chart zu einem nicht-outerchart-punkt in jedem randpunkt einer chart
      
      ARRAY<int> acttrigs; //outercharttrigs
      acttrigs.SetSize(GetNT());
      for (i = 1; i <= GetNT(); i++)
	{
	  acttrigs.Elem(i) = 0;
	}
      for (i = 1; i <= GetNOCharts(); i++)
	{
	  SetThreadPercent((double)i/(double)GetNOCharts()*100.);
	  if (multithread.terminate)
	    {PopStatus(); return;}

	  RestrictHChartDistOneChart(i, acttrigs, mesh, gh, 1., 0.);
	}
      
      PopStatus();
    }

  if (stlparam.resthlinelengthenable)
    {
      //restrict h due to short lines
      PushStatusF("Restrict H due to line-length");
      
      double minhl = 1E50;
      double linefact = 1./stlparam.resthlinelengthfac;
      double l;
      for (i = 1; i <= GetNLines(); i++)
	{
	  SetThreadPercent((double)i/(double)GetNLines()*100.);
	  if (multithread.terminate)
	    {PopStatus(); return;}
	  
	  l = GetLine(i)->GetLength(points);
	  
	  const Point3d& pp1 = GetPoint(GetLine(i)->StartP());
	  const Point3d& pp2 = GetPoint(GetLine(i)->EndP());
	  
	  if (l != 0)
	    {
	      minhl = min2(minhl,l*linefact);
	      
	      mesh.RestrictLocalH(pp1, l*linefact);
	      mesh.RestrictLocalH(pp2, l*linefact);      
	    }
	}
      PopStatus();
      PrintMessage(5, "minh due to line length=", minhl);
  }
}

void STLGeometry :: RestrictHChartDistOneChart(int chartnum, ARRAY<int>& acttrigs, 
					       class Mesh & mesh, double gh, double fact, double minh)
{
  int i = chartnum;
  int j;

  double limessafety = stlparam.resthchartdistfac*fact;  // original: 2
  double localh;

  double f1,f2;
  //  mincalch = 1E10;
  //maxcalch = -1E10;  
  ARRAY<int> limes1;
  ARRAY<int> limes2;
	  
  ARRAY<Point3d> plimes1;
  ARRAY<Point3d> plimes2;
	  
  ARRAY<int> plimes1trigs; //check from wich trig the points come
  ARRAY<int> plimes2trigs;
	  
  ARRAY<int> plimes1origin; //either the original pointnumber or zero, if new point

  int divisions = 10;
	  
  int k, t, nt, np1, np2;
  Point3d p3p1, p3p2;
  STLTriangle tt;
      
  limes1.SetSize(0);
  limes2.SetSize(0);
  plimes1.SetSize(0);
  plimes2.SetSize(0);
  plimes1trigs.SetSize(0);
  plimes2trigs.SetSize(0);
  plimes1origin.SetSize(0);

  STLChart& chart = GetChart(i);
  chart.ClearOLimit();
  chart.ClearILimit();

  for (j = 1; j <= chart.GetNChartT(); j++)
    {
      t = chart.GetChartTrig(j); 
      tt = GetTriangle(t);
      for (k = 1; k <= 3; k++)
	{
	  nt = NeighbourTrig(t,k); 
	  if (GetChartNr(nt) != i)
	    {	      
	      tt.GetNeighbourPoints(GetTriangle(nt),np1,np2);
	      if (!IsEdge(np1,np2) && !GetSpiralPoint(np1) && !GetSpiralPoint(np2))
		{
		  p3p1 = GetPoint(np1);
		  p3p2 = GetPoint(np2);
		  if (AddIfNotExists(limes1,np1)) 
		    {
		      plimes1.Append(p3p1); 
		      plimes1trigs.Append(t);
		      plimes1origin.Append(np1); 			      
		    }
		  if (AddIfNotExists(limes1,np2)) 
		    {
		      plimes1.Append(p3p2); 
		      plimes1trigs.Append(t);
		      plimes1origin.Append(np2); 			      
		    }
		  chart.AddILimit(twoint(np1,np2));

		  for (int di = 1; di <= divisions; di++)
		    {
		      f1 = (double)di/(double)(divisions+1.);
		      f2 = (divisions+1.-(double)di)/(double)(divisions+1.);
			      
		      plimes1.Append(Point3d(p3p1.X()*f1+p3p2.X()*f2,
					     p3p1.Y()*f1+p3p2.Y()*f2,
					     p3p1.Z()*f1+p3p2.Z()*f2));
		      plimes1trigs.Append(t);
		      plimes1origin.Append(0); 			      
		    }
		}
	    }
	}
    }
	  
	 
  for (j = 1; j <= chart.GetNT(); j++)
    {
      acttrigs.Elem(chart.GetTrig(j)) = i;
    }
	  
  for (j = 1; j <= chart.GetNOuterT(); j++)
    {
      t = chart.GetOuterTrig(j); 
      tt = GetTriangle(t);
      for (k = 1; k <= 3; k++)
	{
	  nt = NeighbourTrig(t,k);

	  if (acttrigs.Get(nt) != i)
	    {
	      tt.GetNeighbourPoints(GetTriangle(nt),np1,np2);
		      
	      if (!IsEdge(np1,np2))
		{
		  p3p1 = GetPoint(np1);
		  p3p2 = GetPoint(np2);
			  
		  if (AddIfNotExists(limes2,np1)) {plimes2.Append(p3p1); plimes2trigs.Append(t);}
		  if (AddIfNotExists(limes2,np2)) {plimes2.Append(p3p2); plimes2trigs.Append(t);}
		  chart.AddOLimit(twoint(np1,np2));

		  for (int di = 1; di <= divisions; di++)
		    {
		      f1 = (double)di/(double)(divisions+1.);
		      f2 = (divisions+1.-(double)di)/(double)(divisions+1.);
			      
		      plimes2.Append(Point3d(p3p1.X()*f1+p3p2.X()*f2,
					     p3p1.Y()*f1+p3p2.Y()*f2,
					     p3p1.Z()*f1+p3p2.Z()*f2));
		      plimes2trigs.Append(t);
		    }
		}
	    }
	}
    }
	  
	  
  double chartmindist = 1E50;

  if (plimes2.Size())
    {
      Box3d bbox;
      bbox.SetPoint (plimes2.Get(1));
      for (j = 2; j <= plimes2.Size(); j++)
	bbox.AddPoint (plimes2.Get(j));
      Point3dTree stree(bbox.PMin(), bbox.PMax());
      for (j = 1; j <= plimes2.Size(); j++)
	stree.Insert (plimes2.Get(j), j);
      ARRAY<int> foundpts;
	  
      for (j = 1; j <= plimes1.Size(); j++)
	{
	  double mindist = 1E50;
	  double dist;

	  const Point3d & ap1 = plimes1.Get(j);
	  double boxs = mesh.GetH (plimes1.Get(j)) * limessafety;

	  Point3d pmin = ap1 - Vec3d (boxs, boxs, boxs);
	  Point3d pmax = ap1 + Vec3d (boxs, boxs, boxs);

	  stree.GetIntersecting (pmin, pmax, foundpts);


	  for (int kk = 1; kk <= foundpts.Size(); kk++)
	    {
	      k = foundpts.Get(kk);
	      dist = Dist2(plimes1.Get(j),plimes2.Get(k));
	      if (dist < mindist) 
		{
		  mindist = dist;
		}
	    }

	  /*
	    const Point3d & ap1 = plimes1.Get(j);
	    double his = mesh.GetH (plimes1.Get(j));

	    double xmin = ap1.X() - his * limessafety;
	    double xmax = ap1.X() + his * limessafety;	      
	    double ymin = ap1.Y() - his * limessafety;
	    double ymax = ap1.Y() + his * limessafety;	      
	    double zmin = ap1.Z() - his * limessafety;
	    double zmax = ap1.Z() + his * limessafety;	      

	    for (k = 1; k <= plimes2.Size(); k++)
	    {
	    const Point3d & ap2 = plimes2.Get(k);
	    if (ap2.X() >= xmin && ap2.X() <= xmax &&
	    ap2.Y() >= ymin && ap2.Y() <= ymax &&
	    ap2.Z() >= zmin && ap2.Z() <= zmax)
	    {
	    dist = Dist2(plimes1.Get(j),plimes2.Get(k));
	    if (dist < mindist) 
	    {
	    mindist = dist;
	    }
	    }
	    }
	  */
	  mindist = sqrt(mindist);
	  localh = mindist/limessafety;

	  if (localh < minh && localh != 0) {localh = minh;} //minh is generally 0! (except make atlas)
	  if (localh < gh && localh > 0)
	    {
	      mesh.RestrictLocalH(plimes1.Get(j), localh);
	      //	      if (mindist < mincalch) {mincalch = mindist;}
	      //	      if (mindist > maxcalch) {maxcalch = mindist;}
	      if (mindist < chartmindist) {chartmindist = mindist;}
	    }
	}
    }

}


//void * STLMeshingDummy (void *)
int STLMeshingDummy (STLGeometry* stlgeometry, Mesh*& mesh,
			    int perfstepsstart, int perfstepsend, char* optstring)
{
  if (perfstepsstart > perfstepsend) return 0;

  multithread.terminate = 0;
  int success = 1;
  //int trialcntouter = 0;

  if (perfstepsstart <= MESHCONST_MESHEDGES)
    {

      mesh = new Mesh();
      mesh -> SetGlobalH (mparam.maxh);
      mesh -> SetLocalH (stlgeometry->GetBoundingBox().PMin() - Vec3d(10, 10, 10),
			 stlgeometry->GetBoundingBox().PMax() + Vec3d(10, 10, 10),
			 mparam.grading);
      mesh -> LoadLocalMeshSize (mparam.meshsizefilename);
      
      success = 0;
  
      //mesh->DeleteMesh();
 
      STLMeshing (*stlgeometry, *mesh);

      stlgeometry->edgesfound = 1;
      stlgeometry->surfacemeshed = 0;
      stlgeometry->surfaceoptimized = 0;
      stlgeometry->volumemeshed = 0;
    }

  if (multithread.terminate)
    return 0;

  if (perfstepsstart <= MESHCONST_MESHSURFACE && 
      perfstepsend >= MESHCONST_MESHSURFACE)
    {

      if (!stlgeometry->edgesfound) 
	{
	  PrintUserError("You have to do 'analyse geometry' first!!!");
	  return 0; 
	}
      if (stlgeometry->surfacemeshed || stlgeometry->surfacemeshed) 
	{
	  PrintUserError("Already meshed. Please start again with 'Analyse Geometry'!!!"); 
	  return 0; 
	}

      success = 0;
      int retval = STLSurfaceMeshing (*stlgeometry, *mesh);
      if (retval == MESHING3_OK)
	{
	  PrintMessage(3,"Success !!!!");
	  stlgeometry->surfacemeshed = 1;
	  stlgeometry->surfaceoptimized = 0;
	  stlgeometry->volumemeshed = 0;
	  success = 1;
	} 
      else if (retval == MESHING3_OUTERSTEPSEXCEEDED)
	{
	  PrintError("Give up because of too many trials. Meshing aborted!");
	}
      else if (retval == MESHING3_TERMINATE)
	{
	  PrintWarning("Meshing Stopped by user!");
	}
      else
	{
	  PrintError("Surface meshing not successful. Meshing aborted!");
	}
      
#ifdef STAT_STREAM
      (*statout) << mesh->GetNSeg() << " & " << endl
		 << mesh->GetNSE() << " & " << endl
		 << GetTime() << " & ";
#endif
    }
  if (multithread.terminate)
    return 0;

  if (success)
    {
      if (perfstepsstart <= MESHCONST_OPTSURFACE && 
	  perfstepsend >= MESHCONST_OPTSURFACE)
	{
	  if (!stlgeometry->edgesfound) 
	    {
	      PrintUserError("You have to do 'meshing->analyse geometry' first!!!"); 
	      return 0; 
	    }
	  if (!stlgeometry->surfacemeshed) 
	    {
	      PrintUserError("You have to do 'meshing->mesh surface' first!!!"); 
	      return 0; 
	    }
	  if (stlgeometry->volumemeshed) 
	    {
	      PrintWarning("Surface optimization with meshed volume is dangerous!!!"); 
	    }

	  if (!optstring || strlen(optstring) == 0)
	    {
	      mparam.optimize2d = "smcm";
	    }
	  else
	    {
	      mparam.optimize2d = optstring;
	    }

	  STLSurfaceOptimization (*stlgeometry, *mesh, mparam);
	  
	  if (stlparam.recalc_h_opt)
	    {
	      mesh -> SetLocalH (stlgeometry->GetBoundingBox().PMin() - Vec3d(10, 10, 10),
				 stlgeometry->GetBoundingBox().PMax() + Vec3d(10, 10, 10),
				 mparam.grading);
	      mesh -> LoadLocalMeshSize (mparam.meshsizefilename);	      
	      mesh -> CalcLocalHFromSurfaceCurvature (stlparam.resthsurfmeshcurvfac);
	      mparam.optimize2d = "cmsmSm";
	      STLSurfaceOptimization (*stlgeometry, *mesh, mparam);
#ifdef STAT_STREAM
	      (*statout) << GetTime() << " & ";
#endif

#ifdef OPENGL
	      extern void Render();
	      Render();
#endif	      
	    }
	  stlgeometry->surfaceoptimized = 1;
	}
      if (multithread.terminate)
	return 0;

      if (perfstepsstart <= MESHCONST_MESHVOLUME && 
	  perfstepsend >= MESHCONST_MESHVOLUME)
	{
	  if (stlgeometry->volumemeshed) 
	    {
	      PrintUserError("Volume already meshed!"); return 0;
	    }

	  if (!stlgeometry->edgesfound) 
	    {
	      PrintUserError("You have to do 'meshing->analyse geometry' first!!!"); 
	      return 0; 
	    }
	  if (!stlgeometry->surfacemeshed) 
	    {
	      PrintUserError("You have to do 'meshing->mesh surface' first!!!"); 
	      return 0; 
	    }
	  if (!stlgeometry->surfaceoptimized) 
	    {
	      PrintWarning("You should do 'meshing->optimize surface' first!!!"); 
	    }


	  PrintMessage(5,"Check Overlapping boundary: ");
	  mesh->FindOpenElements();
	  mesh->CheckOverlappingBoundary();
	  PrintMessage(5,"");


	  if (stlparam.recalc_h_opt)
	    {
	      mesh -> SetLocalH (stlgeometry->GetBoundingBox().PMin() - Vec3d(10, 10, 10),
				 stlgeometry->GetBoundingBox().PMax() + Vec3d(10, 10, 10),
				 mparam.grading);	  
	      mesh -> LoadLocalMeshSize (mparam.meshsizefilename);
	      mesh -> CalcLocalH ();
	    }
	  
	  
	  PrintMessage(5,"Volume meshing");
	  int retval = MeshVolume (mparam, *mesh);
	  if (retval == MESHING3_OK)
	    {
	      RemoveIllegalElements(*mesh);
	      stlgeometry->volumemeshed = 1;
	    } 
	  else if (retval == MESHING3_OUTERSTEPSEXCEEDED)
	    {
	      PrintError("Give up because of too many trials. Meshing aborted!");
	      return 0;
	    }
	  else if (retval == MESHING3_TERMINATE)
	    {
	      PrintWarning("Meshing Stopped by user!");
	    }
	  else
	    {
	      PrintError("Volume meshing not successful. Meshing aborted!");
	      return 0;
	    }

#ifdef STAT_STREAM
	  (*statout) << GetTime() << " & " << endl;
#endif
	  MeshQuality3d (*mesh);
	}

      if (multithread.terminate)
	return 0;

      if (perfstepsstart <= MESHCONST_OPTVOLUME && 
	  perfstepsend >= MESHCONST_OPTVOLUME)
	{
	  if (!stlgeometry->edgesfound) 
	    {
	      PrintUserError("You have to do 'meshing->analyse geometry' first!!!"); 
	      return 0; 
	    }
	  if (!stlgeometry->surfacemeshed) 
	    {
	      PrintUserError("You have to do 'meshing->mesh surface' first!!!"); 
	      return 0; 
	    }
	  if (!stlgeometry->volumemeshed) 
	    {
	      PrintUserError("You have to do 'meshing->mesh volume' first!!!"); 
	      return 0; 
	    }

	  if (!optstring || strlen(optstring) == 0)
	    {
	      mparam.optimize3d = "cmdmstm";
	    }
	  else
	    {
	      mparam.optimize3d = optstring;
	    }


	  OptimizeVolume (mparam, *mesh);
	  
#ifdef STAT_STREAM
	  (*statout) << GetTime() << " & " << endl;
	  (*statout) << mesh->GetNE() << " & " << endl
		     << mesh->GetNP() << " " << '\\' << '\\' << " \\" << "hline" << endl;
#endif

#ifdef OPENGL
	  extern void Render();
	  Render();
#endif	      

	}
    }
  

  return 0;
}



}
