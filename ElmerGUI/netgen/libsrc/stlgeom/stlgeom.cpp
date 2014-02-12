#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <gprim.hpp>

#include <meshing.hpp>

#include "stlgeom.hpp"


namespace netgen
{

//globalen searchtree fuer gesamte geometry aktivieren
int geomsearchtreeon = 0;

int usechartnormal = 1;  

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void STLMeshing (STLGeometry & geom,
		 Mesh & mesh)
{
  geom.Clear();
  geom.BuildEdges();
  geom.MakeAtlas(mesh);
  geom.CalcFaceNums();
  geom.AddFaceEdges();
  geom.LinkEdges();

  mesh.ClearFaceDescriptors();
  for (int i = 1; i <= geom.GetNOFaces(); i++)
    mesh.AddFaceDescriptor (FaceDescriptor (i, 1, 0, 0));
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++   STL GEOMETRY   ++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


STLGeometry :: STLGeometry()
  : edges(), edgesperpoint(),
    normals(),  externaledges(),
    atlas(), chartmark(), 
    lines(), outerchartspertrig(), vicinity(), markedtrigs(), markedsegs(),
    lineendpoints(), spiralpoints(), selectedmultiedge()
{
	edgedata = new STLEdgeDataList(*this);
  externaledges.SetSize(0);
  Clear();
  meshchart = 0; // initialize all ?? JS

  if (geomsearchtreeon)
    searchtree = new Box3dTree (GetBoundingBox().PMin() - Vec3d(1,1,1),
				GetBoundingBox().PMax() + Vec3d(1,1,1));
  else
    searchtree = NULL;

  status = STL_GOOD;
  statustext = "Good Geometry";
  smoothedges = NULL;
}

STLGeometry :: ~STLGeometry()
{
  delete edgedata;
}

void STLGeometry :: STLInfo(double* data)
{
  data[0] = GetNT();

  Box<3> b = GetBoundingBox();
  data[1] = b.PMin()(0);
  data[2] = b.PMax()(0);
  data[3] = b.PMin()(1);
  data[4] = b.PMax()(1);
  data[5] = b.PMin()(2);
  data[6] = b.PMax()(2);

  int i;
 
  int cons = 1;
  for (i = 1; i <= GetNT(); i++)
    {
      if (NONeighbourTrigs(i) != 3) {cons = 0;}
    }
  data[7] = cons;
}

void STLGeometry :: MarkNonSmoothNormals()
{

  PrintFnStart("Mark Non-Smooth Normals");

  int i,j;

  markedtrigs.SetSize(GetNT());

  for (i = 1; i <= GetNT(); i++)
    {
      SetMarkedTrig(i, 0);
    }

  double dirtyangle = stlparam.yangle/180.*M_PI;

  int cnt = 0;
  int lp1,lp2;
  for (i = 1; i <= GetNT(); i++)
    {
      for (j = 1; j <= NONeighbourTrigs(i); j++)
	{
	  if (GetAngle(i, NeighbourTrig(i,j)) > dirtyangle)
	    {
	      GetTriangle(i).GetNeighbourPoints(GetTriangle(NeighbourTrig(i,j)), lp1, lp2);
	      if (!IsEdge(lp1,lp2))
		{
		  if (!IsMarkedTrig(i)) {SetMarkedTrig(i,1); cnt++;}
		}
	    }
	}
    }

  PrintMessage(5,"marked ",cnt," non-smooth trig-normals");

}

void STLGeometry :: SmoothNormals()
{
  multithread.terminate = 0;

  //  UseExternalEdges();

  BuildEdges();


  DenseMatrix m(3), hm(3);
  Vector rhs(3), sol(3), hv(3), hv2(3);

  Vec<3> ri;

  double wnb = stldoctor.smoothnormalsweight;   // neigbour normal weight
  double wgeom = 1-wnb;   // geometry normal weight


  // minimize 
  //  wgeom sum_T  \sum ri  \| ri^T (n - n_geom) \|^2  
  //  + wnb sum_SE  \| ri x (n - n_nb) \|^2
  
  int i, j, k, l;
  int nt = GetNT();
  
  PushStatusF("Smooth Normals");
    
  //int testmode;

  for (i = 1; i <= nt; i++)
    {

      SetThreadPercent( 100.0 * (double)i / (double)nt);

      const STLTriangle & trig = GetTriangle (i);
      
      m = 0;
      rhs = 0;

      // normal of geometry:
      Vec<3> ngeom = trig.GeomNormal(points);
      ngeom.Normalize();

      for (j = 1; j <= 3; j++)
	{ 
	  int pi1 = trig.PNumMod (j);
	  int pi2 = trig.PNumMod (j+1);

	  // edge vector
	  ri = GetPoint (pi2) - GetPoint (pi1);
	  
	  for (k = 0; k < 3; k++)
	    for (l = 0; l < 3; l++)
	      hm.Elem(k+1, l+1) = wgeom * ri(k) * ri(l);
	  
	  
	  for (k = 0; k < 3; k++)
	    hv.Elem(k+1) = ngeom(k);
	  
	  hm.Mult (hv, hv2);
	  /*
	  if (testmode)
	    (*testout) << "add vec " << hv2 << endl 
		       << " add m " << hm << endl;
	  */
	  rhs.Add (1, hv2);
	  m += hm;


	  int nbt = 0;
	  int fp1,fp2;
	  for (k = 1; k <= NONeighbourTrigs(i); k++)
	    {
	      trig.GetNeighbourPoints(GetTriangle(NeighbourTrig(i, k)),fp1,fp2);
	      if (fp1 == pi1 && fp2 == pi2)
		{
		  nbt = NeighbourTrig(i, k);
		}
	    }

	  if (!nbt)
	    {
	      cerr << "ERROR: stlgeom::Smoothnormals, nbt = 0" << endl;
	    }

	  // smoothed normal
	  Vec<3> nnb = GetTriangle(nbt).Normal();   // neighbour normal
	  nnb.Normalize();

	  if (!IsEdge(pi1,pi2)) 
	    {
	      double lr2 = ri * ri;
	      for (k = 0; k < 3; k++)
		{
		  for (l = 0; l < k; l++)
		    {
		      hm.Elem(k+1, l+1) = -wnb * ri(k) * ri(l);
		      hm.Elem(l+1, k+1) = -wnb * ri(k) * ri(l);
		    }
		  
		  hm.Elem(k+1, k+1) = wnb * (lr2 - ri(k) * ri(k));
		}
	      
	      for (k = 0; k < 3; k++)
		hv.Elem(k+1) = nnb(k);
	      
	      hm.Mult (hv, hv2);
	      /*
	      if (testmode)
		(*testout) << "add nb vec " << hv2 << endl 
			   << " add nb m " << hm << endl;
	      */

	      rhs.Add (1, hv2);
	      m += hm;
	    }
	}

      m.Solve (rhs, sol);
      Vec3d newn(sol.Get(1), sol.Get(2), sol.Get(3));
      newn /= (newn.Length() + 1e-24);      

      GetTriangle(i).SetNormal(newn);
      // setnormal (sol);
    }

  /*
  for (i = 1; i <= nt; i++)
    SetMarkedTrig(i, 0);      		



  int crloop;
  for (crloop = 1; crloop <= 3; crloop++)
    {

  // find critical:

  ARRAY<INDEX_2> critpairs;
  for (i = 1; i <= nt; i++)
    {
      const STLTriangle & trig = GetTriangle (i);
      
      Vec3d ngeom = GetTriangleNormal (i); // trig.Normal(points);
      ngeom /= (ngeom.Length() + 1e-24);

      for (j = 1; j <= 3; j++)
	{ 
	  int pi1 = trig.PNumMod (j);
	  int pi2 = trig.PNumMod (j+1);

	  int nbt = 0;
	  int fp1,fp2;
	  for (k = 1; k <= NONeighbourTrigs(i); k++)
	    {
	      trig.GetNeighbourPoints(GetTriangle(NeighbourTrig(i, k)),fp1,fp2);
	      if (fp1 == pi1 && fp2 == pi2)
		{
		  nbt = NeighbourTrig(i, k);
		}
	    }
	  
	  if (!nbt)
	    {
	      cerr << "ERROR: stlgeom::Smoothnormals, nbt = 0" << endl;
	    }

	  Vec3d nnb = GetTriangleNormal(nbt);   // neighbour normal
	  nnb /= (nnb.Length() + 1e-24);

	  if (!IsEdge(pi1,pi2)) 
	    {
	      if (Angle (nnb, ngeom) > 150 * M_PI/180)
		{
		  SetMarkedTrig(i, 1);      		
		  SetMarkedTrig(nbt, 1);      		
		  critpairs.Append (INDEX_2 (i, nbt));
		}
	    }

	}
    }

  if (!critpairs.Size())
    {
      break;
    }

  if (critpairs.Size())
    {

      ARRAY<int> friends;
      double area1 = 0, area2 = 0;

      for (i = 1; i <= critpairs.Size(); i++)
	{
	  int tnr1 = critpairs.Get(i).I1();
	  int tnr2 = critpairs.Get(i).I2();
	  (*testout) << "t1 = " << tnr1 << ", t2 = " << tnr2
		     << " angle = " << Angle (GetTriangleNormal (tnr1),
					      GetTriangleNormal (tnr2))
		     << endl;

	  // who has more friends ?
	  int side;
	  area1 = 0;
	  area2 = 0;
	  for (side = 1; side <= 2; side++)
	    {
	      friends.SetSize (0);
	      friends.Append ( (side == 1) ? tnr1 : tnr2);

	      for (j = 1; j <= 3; j++)
		{
		  int fsize = friends.Size();
		  for (k = 1; k <= fsize; k++)
		    {
		      int testtnr = friends.Get(k);
		      Vec3d ntt = GetTriangleNormal(testtnr);
		      ntt /= (ntt.Length() + 1e-24);
		      
		      for (l = 1; l <= NONeighbourTrigs(testtnr); l++)
			{
			  int testnbnr = NeighbourTrig(testtnr, l);
			  Vec3d nbt = GetTriangleNormal(testnbnr);
			  nbt /= (nbt.Length() + 1e-24);

			  if (Angle (nbt, ntt) < 15 * M_PI/180)
			    {
			      int ii;
			      int found = 0;
			      for (ii = 1; ii <= friends.Size(); ii++)
				{
				  if (friends.Get(ii) == testnbnr)
				    {
				      found = 1;
				      break;
				    }
				}
			      if (!found)
				friends.Append (testnbnr);
			    }
			}
		    }
		}

	      // compute area:
	      for (k = 1; k <= friends.Size(); k++)
		{
		  double area = 
		    GetTriangle (friends.Get(k)).Area(points);

		  if (side == 1)
		    area1 += area;
		  else
		    area2 += area;
		}
	      
	    }

	  (*testout) << "area1 = " << area1 << " area2 = " << area2 << endl;
	  if (area1 < 0.1 * area2)
	    {
	      Vec3d n = GetTriangleNormal (tnr1);
	      n *= -1;
	      SetTriangleNormal(tnr1, n);
	    }
	  if (area2 < 0.1 * area1)
	    {
	      Vec3d n = GetTriangleNormal (tnr2);
	      n *= -1;
	      SetTriangleNormal(tnr2, n);
	    }
	}
    }
    }
  */

  calcedgedataanglesnew = 1;
  PopStatus();
}


int STLGeometry :: AddEdge(int ap1, int ap2)
{
  STLEdge e(ap1,ap2);
  e.SetLeftTrig(GetLeftTrig(ap1,ap2));
  e.SetRightTrig(GetRightTrig(ap1,ap2));
  return edges.Append(e);
}

void STLGeometry :: STLDoctorConfirmEdge()
{
  StoreEdgeData();
  if (GetSelectTrig() >= 1 && GetSelectTrig() <= GetNT() && GetNodeOfSelTrig())
    {
      if (stldoctor.selectmode == 1)
	{
	  int ap1 = GetTriangle(GetSelectTrig()).PNum(GetNodeOfSelTrig());
	  int ap2 = GetTriangle(GetSelectTrig()).PNumMod(GetNodeOfSelTrig()+1);
	  edgedata->Elem(edgedata->GetEdgeNum(ap1,ap2)).SetStatus (ED_CONFIRMED);
	}
      else if (stldoctor.selectmode == 3 || stldoctor.selectmode == 4)
	{
	  int i;
	  for (i = 1; i <= selectedmultiedge.Size(); i++)
	    {
	      int ap1 = selectedmultiedge.Get(i).i1;
	      int ap2 = selectedmultiedge.Get(i).i2;
	      edgedata->Elem(edgedata->GetEdgeNum(ap1,ap2)).SetStatus (ED_CONFIRMED);
	    }
	}
    }
}

void STLGeometry :: STLDoctorCandidateEdge()
{
  StoreEdgeData();
  if (GetSelectTrig() >= 1 && GetSelectTrig() <= GetNT() && GetNodeOfSelTrig())
    {
      if (stldoctor.selectmode == 1)
	{
	  int ap1 = GetTriangle(GetSelectTrig()).PNum(GetNodeOfSelTrig());
	  int ap2 = GetTriangle(GetSelectTrig()).PNumMod(GetNodeOfSelTrig()+1);
	  edgedata->Elem(edgedata->GetEdgeNum(ap1,ap2)).SetStatus (ED_CANDIDATE);
	}
      else if (stldoctor.selectmode == 3 || stldoctor.selectmode == 4)
	{
	  int i;
	  for (i = 1; i <= selectedmultiedge.Size(); i++)
	    {
	      int ap1 = selectedmultiedge.Get(i).i1;
	      int ap2 = selectedmultiedge.Get(i).i2;
	      edgedata->Elem(edgedata->GetEdgeNum(ap1,ap2)).SetStatus (ED_CANDIDATE);
	    }
	}
    }
}

void STLGeometry :: STLDoctorExcludeEdge()
{
  StoreEdgeData();
  if (GetSelectTrig() >= 1 && GetSelectTrig() <= GetNT() && GetNodeOfSelTrig())
    {
      if (stldoctor.selectmode == 1)
	{
	  int ap1 = GetTriangle(GetSelectTrig()).PNum(GetNodeOfSelTrig());
	  int ap2 = GetTriangle(GetSelectTrig()).PNumMod(GetNodeOfSelTrig()+1);
	  edgedata->Elem(edgedata->GetEdgeNum(ap1,ap2)).SetStatus(ED_EXCLUDED);
	}
      else if (stldoctor.selectmode == 3 || stldoctor.selectmode == 4)
	{
	  int i;
	  for (i = 1; i <= selectedmultiedge.Size(); i++)
	    {
	      int ap1 = selectedmultiedge.Get(i).i1;
	      int ap2 = selectedmultiedge.Get(i).i2;
	      edgedata->Elem(edgedata->GetEdgeNum(ap1,ap2)).SetStatus(ED_EXCLUDED);
	    }
	}
    }
}

void STLGeometry :: STLDoctorUndefinedEdge()
{
  StoreEdgeData();
  if (GetSelectTrig() >= 1 && GetSelectTrig() <= GetNT() && GetNodeOfSelTrig())
    {
      if (stldoctor.selectmode == 1)
	{
	  int ap1 = GetTriangle(GetSelectTrig()).PNum(GetNodeOfSelTrig());
	  int ap2 = GetTriangle(GetSelectTrig()).PNumMod(GetNodeOfSelTrig()+1);
	  edgedata->Elem(edgedata->GetEdgeNum(ap1,ap2)).SetStatus(ED_UNDEFINED);
	}
      else if (stldoctor.selectmode == 3 || stldoctor.selectmode == 4)
	{
	  int i;
	  for (i = 1; i <= selectedmultiedge.Size(); i++)
	    {
	      int ap1 = selectedmultiedge.Get(i).i1;
	      int ap2 = selectedmultiedge.Get(i).i2;
	      edgedata->Elem(edgedata->GetEdgeNum(ap1,ap2)).SetStatus(ED_UNDEFINED);
	    }
	}
    }
}

void STLGeometry :: STLDoctorSetAllUndefinedEdges()
{
  edgedata->ResetAll();
}

void STLGeometry :: STLDoctorEraseCandidateEdges()
{
  StoreEdgeData();
  edgedata->ChangeStatus(ED_CANDIDATE, ED_UNDEFINED);
}

void STLGeometry :: STLDoctorConfirmCandidateEdges()
{
  StoreEdgeData();
  edgedata->ChangeStatus(ED_CANDIDATE, ED_CONFIRMED);
}

void STLGeometry :: STLDoctorConfirmedToCandidateEdges()
{
  StoreEdgeData();
  edgedata->ChangeStatus(ED_CONFIRMED, ED_CANDIDATE);
}

void STLGeometry :: STLDoctorDirtyEdgesToCandidates()
{
  StoreEdgeData();
}

void STLGeometry :: STLDoctorLongLinesToCandidates()
{
  StoreEdgeData();
}

twoint STLGeometry :: GetNearestSelectedDefinedEdge()
{
  Point<3> pestimate = Center(GetTriangle(GetSelectTrig()).center,
  			     GetPoint(GetTriangle(GetSelectTrig()).PNum(GetNodeOfSelTrig())));
    //Point3d pestimate = GetTriangle(GetSelectTrig()).center;

  int i, j, en;
  ARRAY<int> vic;
  GetVicinity(GetSelectTrig(),4,vic);
  

  twoint fedg;
  fedg.i1 = 0;
  fedg.i2 = 0;
  double mindist = 1E50;
  double dist;
  Point<3> p;

  for (i = 1; i <= vic.Size(); i++)
  {
    const STLTriangle& t = GetTriangle(vic.Get(i));
    for (j = 1; j <= 3; j++)
      {
	en = edgedata->GetEdgeNum(t.PNum(j),t.PNumMod(j+1));
	if (edgedata->Get(en).GetStatus() != ED_UNDEFINED)
	  {
	    p = pestimate;
	    dist = GetDistFromLine(GetPoint(t.PNum(j)),GetPoint(t.PNumMod(j+1)),p);
	    if (dist < mindist)
	      {
		mindist = dist;
		fedg.i1 = t.PNum(j);
		fedg.i2 = t.PNumMod(j+1);
	      }
	  }
      }
  }
  return fedg;
}
 
void STLGeometry :: BuildSelectedMultiEdge(twoint ep)
{
  if (edgedata->Size() == 0 || 
      !GetEPPSize()) 
    {
      return; 
    }

  selectedmultiedge.SetSize(0);
  int tenum = GetTopEdgeNum (ep.i1, ep.i2);

  if (edgedata->Get(tenum).GetStatus() == ED_UNDEFINED)
    {
      twoint epnew = GetNearestSelectedDefinedEdge();
      if (epnew.i1) 
	{
	  ep = epnew;
	  tenum = GetTopEdgeNum (ep.i1, ep.i2);
	}
    }

  selectedmultiedge.Append(twoint(ep));

  if (edgedata->Get(tenum).GetStatus() == ED_UNDEFINED)
    {
      return;
    }

  edgedata->BuildLineWithEdge(ep.i1,ep.i2,selectedmultiedge);
}

void STLGeometry :: BuildSelectedEdge(twoint ep)
{
  if (edgedata->Size() == 0 || 
      !GetEPPSize()) 
    {
      return; 
    }

  selectedmultiedge.SetSize(0);

  selectedmultiedge.Append(twoint(ep));
}

void STLGeometry :: BuildSelectedCluster(twoint ep)
{
  if (edgedata->Size() == 0 || 
      !GetEPPSize()) 
    {
      return; 
    }

  selectedmultiedge.SetSize(0);

  int tenum = GetTopEdgeNum (ep.i1, ep.i2);

  if (edgedata->Get(tenum).GetStatus() == ED_UNDEFINED)
    {
      twoint epnew = GetNearestSelectedDefinedEdge();
      if (epnew.i1) 
	{
	  ep = epnew;
	  tenum = GetTopEdgeNum (ep.i1, ep.i2);
	}
    }

  selectedmultiedge.Append(twoint(ep));

  if (edgedata->Get(tenum).GetStatus() == ED_UNDEFINED)
    {
      return;
    }

  edgedata->BuildClusterWithEdge(ep.i1,ep.i2,selectedmultiedge);
}

void STLGeometry :: ImportEdges()
{
  StoreEdgeData();

  PrintMessage(5, "import edges from file 'edges.ng'");
  ifstream fin("edges.ng");

  int ne;
  fin >> ne;

  ARRAY<Point<3> > eps;

  int i;
  Point<3> p;
  for (i = 1; i <= 2*ne; i++)
    {
      fin >> p(0); 
      fin >> p(1); 
      fin >> p(2);
      eps.Append(p);
    }
  AddEdges(eps);
}

void STLGeometry :: AddEdges(const ARRAY<Point<3> >& eps)
{
  int i;
  int ne = eps.Size()/2;
  
  ARRAY<int> epsi;
  Box<3> bb = GetBoundingBox();
  bb.Increase(1);

  Point3dTree ptree (bb.PMin(), 
			 bb.PMax());
  ARRAY<int> pintersect;

  double gtol = GetBoundingBox().Diam()/1.E10;
  Point<3> p;

  for (i = 1; i <= GetNP(); i++)
    {
      p = GetPoint(i);
      ptree.Insert (p, i);
    }
  
  int error = 0;
  for (i = 1; i <= 2*ne; i++)
    {
      p = eps.Get(i);
      Point3d pmin = p - Vec3d (gtol, gtol, gtol);
      Point3d pmax = p + Vec3d (gtol, gtol, gtol);
	  
      ptree.GetIntersecting (pmin, pmax, pintersect);
      if (pintersect.Size() > 1)
	{
	  PrintError("Found too much points in epsilon-dist");
	  error = 1;
	}
      else if (pintersect.Size() == 0)
	{
	  error = 1;
	  PrintError("edgepoint does not exist!");
	  PrintMessage(5,"p=",Point3d(eps.Get(i)));
	}
      else
	{
	  epsi.Append(pintersect.Get(1));
	}
    }

  if (error) return;

  int en;
  for (i = 1; i <= ne; i++)
    {
      if (epsi.Get(2*i-1) == epsi.Get(2*i)) {PrintError("Edge with zero length!");}
      else 
	{
	  en = edgedata->GetEdgeNum(epsi.Get(2*i-1),epsi.Get(2*i));
	  edgedata->Elem(en).SetStatus (ED_CONFIRMED);
	}
    }

}



void STLGeometry :: ImportExternalEdges(const char * filename)
{
  //AVL edges!!!!!!

  ifstream inf (filename);
  char ch;
  //int cnt = 0;
  int records, units, i, j;
  PrintFnStart("Import edges from ",filename);
  
  const int flen=30;
  char filter[flen+1];
  filter[flen] = 0;
  char buf[20];

  ARRAY<Point3d> importpoints;
  ARRAY<int> importlines;
  ARRAY<int> importpnums;

  while (inf.good())
    {
      inf.get(ch);
      //      (*testout) << cnt << ": " << ch << endl;
      
      for (i = 0; i < flen; i++)
	filter[i] = filter[i+1];
      filter[flen-1] = ch;
      //      (*testout) << filter << endl;

      if (strcmp (filter+flen-7, "RECORDS") == 0)
	{
	  inf.get(ch);  // '='
	  inf >> records;
	}
      if (strcmp (filter+flen-5, "UNITS") == 0)
	{
	  inf.get(ch);  // '='
	  inf >> units;
	}

      if (strcmp (filter+flen-17, "EDGE NODE NUMBERS") == 0)
	{
	  int nodenr;
	  importlines.SetSize (units);
	  for (i = 1; i <= units; i++)
	    {
	      inf >> nodenr;
	      importlines.Elem(i) = nodenr;
	      //	      (*testout) << nodenr << endl;
	    }
	}

      if (strcmp (filter+flen-23, "EDGE POINT COORD IN DIR") == 0)
	{
	  int coord;

	  inf >> coord;
	  
	  importpoints.SetSize (units);

	  inf >> ch;
	  inf.putback (ch);

	  for (i = 1; i <= units; i++)
	    {
	      for (j = 0; j < 12; j++)
		inf.get (buf[j]);
	      buf[12] = 0;

	      importpoints.Elem(i).X(coord) = 1000 * atof (buf);
	    }
	}
    }

  /*
  (*testout) << "lines: " << endl;
  for (i = 1; i <= importlines.Size(); i++)
    (*testout) << importlines.Get(i) << endl;
  (*testout) << "points: " << endl;
  for (i = 1; i <= importpoints.Size(); i++)
    (*testout) << importpoints.Get(i) << endl;
  */



  importpnums.SetSize (importpoints.Size());
  

  Box3d bb (GetBoundingBox().PMin() + Vec3d (-1,-1,-1),
	    GetBoundingBox().PMax() + Vec3d (1, 1, 1));

  Point3dTree ptree (bb.PMin(), 
			 bb.PMax());


  PrintMessage(7,"stl - bb: ",bb.PMin(), " - ", bb.PMax());
  
  Box3d ebb;
  ebb.SetPoint (importpoints.Get(1));
  for (i = 1; i <= importpoints.Size(); i++)
    ebb.AddPoint (importpoints.Get(i));
  PrintMessage(7,"edgep - bb: ", ebb.PMin(), " - ", ebb.PMax());

  ARRAY<int> pintersect;

  double gtol = GetBoundingBox().Diam()/1.E6;

  for (i = 1; i <= GetNP(); i++)
    {
      Point3d p = GetPoint(i);
      //      (*testout) << "stlpt: " << p << endl;
      ptree.Insert (p, i);
    }
  

  for (i = 1; i <= importpoints.Size(); i++)
    {
      Point3d p = importpoints.Get(i);
      Point3d pmin = p - Vec3d (gtol, gtol, gtol);
      Point3d pmax = p + Vec3d (gtol, gtol, gtol);
	  
      ptree.GetIntersecting (pmin, pmax, pintersect);
      if (pintersect.Size() > 1)
	{
	  importpnums.Elem(i) = 0;
	  PrintError("Found too many points in epsilon-dist");
	}
      else if (pintersect.Size() == 0)
	{
	  importpnums.Elem(i) = 0;
	  PrintError("Edgepoint does not exist!");
	}
      else
	{
	  importpnums.Elem(i) = pintersect.Get(1);
	}
    }

  //  if (!error) 
    {
      PrintMessage(7,"found all edge points in stl file");


      StoreEdgeData();

      int oldp = 0;

      for (i = 1; i <= importlines.Size(); i++)
	{
	  int newp = importlines.Get(i);
	  if (!importpnums.Get(abs(newp)))
	    newp = 0;

	  if (oldp && newp)
	    {
	      int en = edgedata->GetEdgeNum(importpnums.Get(oldp), 
					   importpnums.Get(abs(newp)));
	      edgedata->Elem(en).SetStatus (ED_CONFIRMED);
	    }
	  
	  if (newp < 0)
	    oldp = 0;
	  else
	    oldp = newp;
	}
    }


}



void STLGeometry :: ExportEdges()
{
  PrintFnStart("Save edges to file 'edges.ng'");

  ofstream fout("edges.ng");
  fout.precision(16);

  int n = edgedata->GetNConfEdges();
  
  fout << n << endl;

  int i;
  for (i = 1; i <= edgedata->Size(); i++)
    {
      if (edgedata->Get(i).GetStatus() == ED_CONFIRMED)
	{
	  const STLTopEdge & e = edgedata->Get(i);
	  fout << GetPoint(e.PNum(1))(0) << " " << GetPoint(e.PNum(1))(1) << " " << GetPoint(e.PNum(1))(2) << endl;
	  fout << GetPoint(e.PNum(2))(0) << " " << GetPoint(e.PNum(2))(1) << " " << GetPoint(e.PNum(2))(2) << endl;
	}
    }

}

void STLGeometry :: LoadEdgeData(const char* file)
{
  StoreEdgeData();

  PrintFnStart("Load edges from file '", file, "'");
  ifstream fin(file);

  edgedata->Read(fin);

  //  calcedgedataanglesnew = 1;
}

void STLGeometry :: SaveEdgeData(const char* file)
{
  PrintFnStart("save edges to file '", file, "'");
  ofstream fout(file);

  edgedata->Write(fout);
}







/*
void STLGeometry :: SaveExternalEdges()
{
  ofstream fout("externaledgesp3.ng");
  fout.precision(16);

  int n = NOExternalEdges();
  fout << n << endl;

  int i;
  for (i = 1; i <= n; i++)
    {
      twoint e = GetExternalEdge(i);
      fout << GetPoint(e.i1)(0) << " " << GetPoint(e.i1)(1) << " " << GetPoint(e.i1)(2) << endl;
      fout << GetPoint(e.i2)(0) << " " << GetPoint(e.i2)(1) << " " << GetPoint(e.i2)(2) << endl;
    }

}
*/
void STLGeometry :: StoreExternalEdges()
{
  storedexternaledges.SetSize(0);
  undoexternaledges = 1;
  int i;
  for (i = 1; i <= externaledges.Size(); i++)
    {
      storedexternaledges.Append(externaledges.Get(i));      
    }

}

void STLGeometry :: UndoExternalEdges()
{
  if (!undoexternaledges) 
    {
      PrintMessage(1, "undo not further possible!");
      return;
    }
  RestoreExternalEdges();
  undoexternaledges = 0;
}

void STLGeometry :: RestoreExternalEdges()
{
  externaledges.SetSize(0);
  int i;
  for (i = 1; i <= storedexternaledges.Size(); i++)
    {
      externaledges.Append(storedexternaledges.Get(i));      
    }

}


void STLGeometry :: AddExternalEdgeAtSelected()
{
  StoreExternalEdges();
  if (GetSelectTrig() >= 1 && GetSelectTrig() <= GetNT())
    {
      int ap1 = GetTriangle(GetSelectTrig()).PNum(GetNodeOfSelTrig());
      int ap2 = GetTriangle(GetSelectTrig()).PNumMod(GetNodeOfSelTrig()+1);
      if (!IsExternalEdge(ap1,ap2)) {AddExternalEdge(ap1,ap2);}
    }
}

void STLGeometry :: AddClosedLinesToExternalEdges()
{
  StoreExternalEdges();

  int i, j;
  for (i = 1; i <= GetNLines(); i++)
    {
      STLLine* l = GetLine(i);
      if (l->StartP() == l->EndP()) 
	{
	  for (j = 1; j < l->NP(); j++)
	    {
	      int ap1 = l->PNum(j);
	      int ap2 = l->PNum(j+1);

	      if (!IsExternalEdge(ap1,ap2)) {AddExternalEdge(ap1,ap2);}	      
	    }
	}
    }
}

void STLGeometry :: AddLongLinesToExternalEdges()
{
  StoreExternalEdges();

  double diamfact = stldoctor.dirtytrigfact;
  double diam = GetBoundingBox().Diam();

  int i, j;
  for (i = 1; i <= GetNLines(); i++)
    {
      STLLine* l = GetLine(i);
      if (l->GetLength(points) >= diamfact*diam) 
	{
	  for (j = 1; j < l->NP(); j++)
	    {
	      int ap1 = l->PNum(j);
	      int ap2 = l->PNum(j+1);

	      if (!IsExternalEdge(ap1,ap2)) {AddExternalEdge(ap1,ap2);}	      
	    }
	}
    }
}

void STLGeometry :: AddAllNotSingleLinesToExternalEdges()
{
  StoreExternalEdges();

  int i, j;
  for (i = 1; i <= GetNLines(); i++)
    {
      STLLine* l = GetLine(i);
      if (GetNEPP(l->StartP()) > 1 || GetNEPP(l->EndP()) > 1) 
	{
	  for (j = 1; j < l->NP(); j++)
	    {
	      int ap1 = l->PNum(j);
	      int ap2 = l->PNum(j+1);

	      if (!IsExternalEdge(ap1,ap2)) {AddExternalEdge(ap1,ap2);}	      
	    }
	}
    }
}

void STLGeometry :: DeleteDirtyExternalEdges()
{
  //delete single triangle edges and single edge-lines in clusters"
  StoreExternalEdges();

  int i, j;
  for (i = 1; i <= GetNLines(); i++)
    {
      STLLine* l = GetLine(i);
      if (l->NP() <= 3 || (l->StartP() == l->EndP() && l->NP() == 4))
	{
	  for (j = 1; j < l->NP(); j++)
	    {
	      int ap1 = l->PNum(j);
	      int ap2 = l->PNum(j+1);

	      if (IsExternalEdge(ap1,ap2)) {DeleteExternalEdge(ap1,ap2);}	      
	    }
	}
    }
}

void STLGeometry :: AddExternalEdgesFromGeomLine()
{
  StoreExternalEdges();
  if (GetSelectTrig() >= 1 && GetSelectTrig() <= GetNT())
    {
      int ap1 = GetTriangle(GetSelectTrig()).PNum(GetNodeOfSelTrig());
      int ap2 = GetTriangle(GetSelectTrig()).PNumMod(GetNodeOfSelTrig()+1);

      if (IsEdge(ap1,ap2))
	{
	  int edgenum = IsEdgeNum(ap1,ap2);
	  if (!IsExternalEdge(ap1,ap2)) {AddExternalEdge(ap1,ap2);}
	  
	  int noend = 1;
	  int startp = ap1;
	  int laste = edgenum;
	  int np1, np2;
	  while (noend)
	    {
	      if (GetNEPP(startp) == 2)
		{
		  if (GetEdgePP(startp,1) != laste) {laste = GetEdgePP(startp,1);}
		  else {laste = GetEdgePP(startp,2);}
		  np1 = GetEdge(laste).PNum(1);
		  np2 = GetEdge(laste).PNum(2);
		  
		  if (!IsExternalEdge(np1, np2)) {AddExternalEdge(np1, np2);}
		  else {noend = 0;}
		  if (np1 != startp) {startp = np1;}
		  else {startp = np2;}
		}
	      else {noend = 0;}
	    }

	  startp = ap2;
	  laste = edgenum;
	  noend = 1;
	  while (noend)
	    {
	      if (GetNEPP(startp) == 2)
		{
		  if (GetEdgePP(startp,1) != laste) {laste = GetEdgePP(startp,1);}
		  else {laste = GetEdgePP(startp,2);}
		  np1 = GetEdge(laste).PNum(1);
		  np2 = GetEdge(laste).PNum(2);
		  
		  if (!IsExternalEdge(np1, np2)) {AddExternalEdge(np1, np2);}
		  else {noend = 0;}
		  if (np1 != startp) {startp = np1;}
		  else {startp = np2;}
		}
	      else {noend = 0;}
	    }
	  
	}

    }
  
}

void STLGeometry :: ClearEdges()
{
  edgesfound = 0;
  edges.SetSize(0);
  //edgedata->SetSize(0);
  // externaledges.SetSize(0);
  edgesperpoint.SetSize(0);
  undoexternaledges = 0;

}

void STLGeometry :: STLDoctorBuildEdges()
{
  //  if (!trigsconverted) {return;}
  ClearEdges();

  meshlines.SetSize(0);
  FindEdgesFromAngles();
}

void STLGeometry :: DeleteExternalEdgeAtSelected()
{
  StoreExternalEdges();
  if (GetSelectTrig() >= 1 && GetSelectTrig() <= GetNT())
    {
      int ap1 = GetTriangle(GetSelectTrig()).PNum(GetNodeOfSelTrig());
      int ap2 = GetTriangle(GetSelectTrig()).PNumMod(GetNodeOfSelTrig()+1);
      if (IsExternalEdge(ap1,ap2)) {DeleteExternalEdge(ap1,ap2);}
    }
}

void STLGeometry :: DeleteExternalEdgeInVicinity()
{
  StoreExternalEdges();
  if (!stldoctor.showvicinity || vicinity.Size() != GetNT()) {return;}

  int i, j, ap1, ap2;
  
  for (i = 1; i <= GetNT(); i++)
    {
      if (vicinity.Elem(i))
	{
	  for (j = 1; j <= 3; j++)
	    {
	      ap1 = GetTriangle(i).PNum(j);
	      ap2 = GetTriangle(i).PNumMod(j+1);

	      if (IsExternalEdge(ap1,ap2))
		{
		  DeleteExternalEdge(ap1,ap2);
		}
	    }
	}
    }
}

void STLGeometry :: BuildExternalEdgesFromEdges()
{
  StoreExternalEdges();

  if (GetNE() == 0) {PrintWarning("Edges possibly not generated!");}

  int i;
  externaledges.SetSize(0);

  for (i = 1; i <= GetNE(); i++)
    {
      STLEdge e = GetEdge(i);
      AddExternalEdge(e.PNum(1), e.PNum(2));
    }

}


void STLGeometry :: AddExternalEdge(int ap1, int ap2)
{
  externaledges.Append(twoint(ap1,ap2));
}

void STLGeometry :: DeleteExternalEdge(int ap1, int ap2)
{

  int i;
  int found = 0;
  for (i = 1; i <= NOExternalEdges(); i++)
    {
      if ((GetExternalEdge(i).i1 == ap1 && GetExternalEdge(i).i2 == ap2) ||
	  (GetExternalEdge(i).i1 == ap2 && GetExternalEdge(i).i2 == ap1)) {found = 1;};
      if (found && i < NOExternalEdges())
	{
	  externaledges.Elem(i) = externaledges.Get(i+1);
	}
    }
  if (!found) {PrintWarning("edge not found");}
  else
    {
      externaledges.SetSize(externaledges.Size()-1);
    }

}

int STLGeometry :: IsExternalEdge(int ap1, int ap2)
{
  int i;
  for (i = 1; i <= NOExternalEdges(); i++)
    {
      if ((GetExternalEdge(i).i1 == ap1 && GetExternalEdge(i).i2 == ap2) ||
	  (GetExternalEdge(i).i1 == ap2 && GetExternalEdge(i).i2 == ap1)) {return 1;};
    }
  return 0;
}

void STLGeometry :: DestroyDirtyTrigs()
{

  PrintFnStart("Destroy dirty triangles");
  PrintMessage(5,"original number of triangles=", GetNT());

  //destroy every triangle with other than 3 neighbours;
  int changed = 1;
  int i, j, k;
  while (changed)
    {
      changed = 0;
      Clear();

      for (i = 1; i <= GetNT(); i++)
	{
	  int dirty = NONeighbourTrigs(i) < 3;

	  for (j = 1; j <= 3; j++)
	    {
	      int pnum = GetTriangle(i).PNum(j);
	      /*
	      if (pnum == 1546)
		{
		// for (k = 1; k <=  NOTrigsPerPoint(pnum); k++)
		}
	      */
	      if (NOTrigsPerPoint(pnum) <= 2) 
		dirty = 1;
	    }
	  
	  int pi1 = GetTriangle(i).PNum(1);
	  int pi2 = GetTriangle(i).PNum(2);
	  int pi3 = GetTriangle(i).PNum(3);
	  if (pi1 == pi2 || pi1 == pi3 || pi2 == pi3)
	    {
	      PrintMessage(5,"triangle with Volume 0: ", i, "  nodes: ", pi1, ", ", pi2, ", ", pi3);
	      dirty = 1;
	    }

	  if (dirty)
	    {
	      for (k = i+1; k <= GetNT(); k++)
		{
		  trias.Elem(k-1) = trias.Get(k);
		  // readtrias: not longer permanent, JS
		  //		  readtrias.Elem(k-1) = readtrias.Get(k); 
		}
	      int size = GetNT();
	      trias.SetSize(size-1);
	      //	      readtrias.SetSize(size-1);
	      changed = 1;
	      break;
	    }
	}
    }  

  FindNeighbourTrigs();
  PrintMessage(5,"final number of triangles=", GetNT());
}

void STLGeometry :: CalcNormalsFromGeometry()
{
  int i;
  for (i = 1; i <= GetNT(); i++)
    {
      const STLTriangle & tr = GetTriangle(i);
      const Point3d& ap1 = GetPoint(tr.PNum(1));
      const Point3d& ap2 = GetPoint(tr.PNum(2));
      const Point3d& ap3 = GetPoint(tr.PNum(3));

      Vec3d normal = Cross (ap2-ap1, ap3-ap1);
      
      if (normal.Length() != 0)
	{
	  normal /= (normal.Length());		  
	}
      GetTriangle(i).SetNormal(normal);
    }
  PrintMessage(5,"Normals calculated from geometry!!!");

  calcedgedataanglesnew = 1;
}

void STLGeometry :: SetSelectTrig(int trig)
{
  stldoctor.selecttrig = trig;
}

int STLGeometry :: GetSelectTrig() const
{
  return stldoctor.selecttrig;
}

void STLGeometry :: SetNodeOfSelTrig(int n)
{
  stldoctor.nodeofseltrig = n;
}

int STLGeometry :: GetNodeOfSelTrig() const
{
  return stldoctor.nodeofseltrig;
}

void STLGeometry :: MoveSelectedPointToMiddle()
{
  if (GetSelectTrig() >= 1 && GetSelectTrig() <= GetNT())
    {
      int p = GetTriangle(GetSelectTrig()).PNum(GetNodeOfSelTrig());
      Point<3> pm(0.,0.,0.); //Middlevector;
      Point<3> p0(0.,0.,0.);
      PrintMessage(5,"original point=", Point3d(GetPoint(p)));

      int i;
      int cnt = 0;
      for (i = 1; i <= trigsperpoint.EntrySize(p); i++)
	{
	  const STLTriangle& tr = GetTriangle(trigsperpoint.Get(p,i));
	  int j;
	  for (j = 1; j <= 3; j++)
	    {
	      if (tr.PNum(j) != p)
		{
		  cnt++;
		  pm(0) += GetPoint(tr.PNum(j))(0);
		  pm(1) += GetPoint(tr.PNum(j))(1);
		  pm(2) += GetPoint(tr.PNum(j))(2);
		}
	    }
	}

      Point<3> origp = GetPoint(p);
      double fact = 0.2;

      SetPoint(p, p0 + fact*(1./(double)cnt)*(pm-p0)+(1.-fact)*(origp-p0));

      PrintMessage(5,"middle point=", Point3d (GetPoint(p)));
      
      PrintMessage(5,"moved point ", Point3d (p));

    }
}

void STLGeometry :: PrintSelectInfo()
{

  //int trig = GetSelectTrig();
  //int p = GetTriangle(trig).PNum(GetNodeOfSelTrig());
  
  PrintMessage(1,"touch triangle ", GetSelectTrig()
       , ", local node ", GetNodeOfSelTrig()
       , " (=", GetTriangle(GetSelectTrig()).PNum(GetNodeOfSelTrig()), ")");
  if (AtlasMade() && GetSelectTrig() >= 1 && GetSelectTrig() <= GetNT())
    {
      PrintMessage(1,"           chartnum=",GetChartNr(GetSelectTrig()));
      /*      
      PointBetween(Center(Center(GetPoint(GetTriangle(270).PNum(1)),
				 GetPoint(GetTriangle(270).PNum(2))),
			  GetPoint(GetTriangle(270).PNum(3))),270,
		   Center(Center(GetPoint(GetTriangle(trig).PNum(1)),
				 GetPoint(GetTriangle(trig).PNum(2))),
			  GetPoint(GetTriangle(trig).PNum(3))),trig);
      */
      //PointBetween(Point3d(5.7818, 7.52768, 4.14879),260,Point3d(6.80292, 6.55392, 4.70184),233);
    }
}

void STLGeometry :: ShowSelectedTrigChartnum()
{
  int st = GetSelectTrig();

  if (st >= 1 && st <= GetNT() && AtlasMade())
    PrintMessage(1,"selected trig ", st, " has chartnumber ", GetChartNr(st));
}

void STLGeometry :: ShowSelectedTrigCoords()
{
  int st = GetSelectTrig();

  /*
  //testing!!!!
  ARRAY<int> trigs;
  GetSortedTrianglesAroundPoint(GetTriangle(st).PNum(GetNodeOfSelTrig()),st,trigs);
  */

  if (st >= 1 && st <= GetNT())
    {
      PrintMessage(1, "coordinates of selected trig ", st, ":");
      PrintMessage(1, "   p1 = ", GetTriangle(st).PNum(1), " = ", 
		   Point3d (GetPoint(GetTriangle(st).PNum(1))));
      PrintMessage(1, "   p2 = ", GetTriangle(st).PNum(2), " = ", 
		   Point3d (GetPoint(GetTriangle(st).PNum(2))));
      PrintMessage(1, "   p3 = ", GetTriangle(st).PNum(3), " = ", 
		   Point3d (GetPoint(GetTriangle(st).PNum(3))));
    }
}

void STLGeometry :: LoadMarkedTrigs()
{
  PrintFnStart("load marked trigs from file 'markedtrigs.ng'");
  ifstream fin("markedtrigs.ng");

  int n;
  fin >> n;
  if (n != GetNT() || n == 0) {PrintError("Not a suitable marked-trig-file!"); return;}

  int i, m;
  for (i = 1; i <= n; i++)
    {
      fin >> m;
      SetMarkedTrig(i, m);      
    }

  fin >> n;
  if (n != 0) 
    {
 
      Point<3> ap1, ap2;
      for (i = 1; i <= n; i++)
	{
	  fin >> ap1(0); fin >> ap1(1); fin >> ap1(2);
	  fin >> ap2(0); fin >> ap2(1); fin >> ap2(2);
	  AddMarkedSeg(ap1,ap2);      
	}
    }
}

void STLGeometry :: SaveMarkedTrigs()
{
  PrintFnStart("save marked trigs to file 'markedtrigs.ng'");
  ofstream fout("markedtrigs.ng");

  int n = GetNT();
  fout << n << endl;

  int i;
  for (i = 1; i <= n; i++)
    {
      fout << IsMarkedTrig(i) << "\n";
    }

  n = GetNMarkedSegs();
  fout << n << endl;

  Point<3> ap1,ap2;
  for (i = 1; i <= n; i++)
    {
      GetMarkedSeg(i,ap1,ap2);
      fout << ap1(0) << " " << ap1(1) << " " << ap1(2) << "  ";
      fout << ap2(0) << " " << ap2(1) << " " << ap2(2) << " " << "\n";
    }

}

void STLGeometry :: NeighbourAnglesOfSelectedTrig()
{
  int st = GetSelectTrig();

  if (st >= 1 && st <= GetNT())
    {
      int i;
      PrintMessage(1,"Angle to triangle ", st, ":");
      for (i = 1; i <= NONeighbourTrigs(st); i++)
	{
	  PrintMessage(1,"   triangle ", NeighbourTrig(st,i), ": angle = ", 
		       180./M_PI*GetAngle(st, NeighbourTrig(st,i)), "°",
		       ", calculated = ", 180./M_PI*Angle(GetTriangle(st).GeomNormal(points), 
							  GetTriangle(NeighbourTrig(st,i)).GeomNormal(points)), "°");
	}
    }
}

void STLGeometry :: GetVicinity(int starttrig, int size, ARRAY<int>& vic)
{
  if (starttrig == 0 || starttrig > GetNT()) {return;} 

  ARRAY<int> vicarray;
  vicarray.SetSize(GetNT());

  int i;
  for (i = 1; i <= vicarray.Size(); i++)
    {
      vicarray.Elem(i) = 0;
    }
 
  vicarray.Elem(starttrig) = 1;
  
  int j = 0,k;

  ARRAY <int> list1;
  list1.SetSize(0);
  ARRAY <int> list2;
  list2.SetSize(0);
  list1.Append(starttrig);

  while (j < size)
    {
      j++;
      for (i = 1; i <= list1.Size(); i++)
	{
	  for (k = 1; k <= NONeighbourTrigs(i); k++)
	    {
	      int nbtrig = NeighbourTrig(list1.Get(i),k);
	      if (nbtrig && vicarray.Get(nbtrig) == 0)
		{
		  list2.Append(nbtrig);
		  vicarray.Elem(nbtrig) = 1;
		}
	    }
	}
      list1.SetSize(0);
      for (i = 1; i <= list2.Size(); i++)
	{
	  list1.Append(list2.Get(i));
	}
      list2.SetSize(0);
    }

  vic.SetSize(0);
  for (i = 1; i <= vicarray.Size(); i++)
    {
      if (vicarray.Get(i)) {vic.Append(i);}
    }
}

void STLGeometry :: CalcVicinity(int starttrig)
{
  if (starttrig == 0 || starttrig > GetNT()) {return;} 

  vicinity.SetSize(GetNT());

  if (!stldoctor.showvicinity) {return;}

  int i;
  for (i = 1; i <= vicinity.Size(); i++)
    {
      vicinity.Elem(i) = 0;
    }
 
  vicinity.Elem(starttrig) = 1;
  
  int j = 0,k;

  ARRAY <int> list1;
  list1.SetSize(0);
  ARRAY <int> list2;
  list2.SetSize(0);
  list1.Append(starttrig);

  //  int cnt = 1;
  while (j < stldoctor.vicinity)
    {
      j++;
      for (i = 1; i <= list1.Size(); i++)
	{
	  for (k = 1; k <= NONeighbourTrigs(i); k++)
	    {
	      int nbtrig = NeighbourTrig(list1.Get(i),k);
	      if (nbtrig && vicinity.Get(nbtrig) == 0)
		{
		  list2.Append(nbtrig);
		  vicinity.Elem(nbtrig) = 1;
		  //cnt++;
		}
	    }
	}
      list1.SetSize(0);
      for (i = 1; i <= list2.Size(); i++)
	{
	  list1.Append(list2.Get(i));
	}
      list2.SetSize(0);
    }

}

int STLGeometry :: Vicinity(int trig) const 
{
  if (trig <= vicinity.Size() && trig >=1)
    {
      return vicinity.Get(trig);
    }
  else {PrintSysError("In STLGeometry::Vicinity");}
  return 0;
}

void STLGeometry :: InitMarkedTrigs()
{
  markedtrigs.SetSize(GetNT());
  int i;
  for (i = 1; i <= GetNT(); i++)
    {
      SetMarkedTrig(i, 0);
    }
}

void STLGeometry :: MarkDirtyTrigs()
{
  PrintFnStart("mark dirty trigs");
  int i,j;

  markedtrigs.SetSize(GetNT());

  for (i = 1; i <= GetNT(); i++)
    {
      SetMarkedTrig(i, 0);
    }

  int found;
  double dirtyangle = stlparam.yangle/2./180.*M_PI;
  int cnt = 0;
  for (i = 1; i <= GetNT(); i++)
    {
      found = 0;
      for (j = 1; j <= NONeighbourTrigs(i); j++)
	{
	  if (GetAngle(i, NeighbourTrig(i,j)) > dirtyangle)
	    {
	      found++;
	    }
	}
      if (found && GetTriangle(i).MinHeight(points) < 
	  stldoctor.dirtytrigfact*GetTriangle(i).MaxLength(points))
	{
	  SetMarkedTrig(i, 1); cnt++;
	}
      /*
      else if (found == 3)
	{
	  SetMarkedTrig(i, 1); cnt++;	  
	}
      */
    }

  PrintMessage(1, "marked ", cnt, " dirty trigs");
}


void STLGeometry :: MarkTopErrorTrigs()
{
  int cnt = 0;
  markedtrigs.SetSize(GetNT());
  for (int i = 1; i <= GetNT(); i++)
    {
      const STLTriangle & trig = GetTriangle(i);

      SetMarkedTrig(i, trig.flags.toperror);
      if (trig.flags.toperror) cnt++;
    }
  PrintMessage(1,"marked ", cnt, " inconsistent triangles");
}



double STLGeometry :: CalcTrigBadness(int i)
{
  int j;
  double maxbadness = 0;
  int ap1, ap2;
  for (j = 1; j <= NONeighbourTrigs(i); j++)
    {
      GetTriangle(i).GetNeighbourPoints(GetTriangle(NeighbourTrig(i,j)), ap1, ap2);
      
      if (!IsEdge(ap1,ap2) && GetGeomAngle(i, NeighbourTrig(i,j)) > maxbadness)
	{
	  maxbadness = GetGeomAngle(i, NeighbourTrig(i,j));
	}
    }
  return maxbadness;

}

void STLGeometry :: GeomSmoothRevertedTrigs()
{
  //double revertedangle = stldoctor.smoothangle/180.*M_PI;
  double fact = stldoctor.dirtytrigfact;

  MarkRevertedTrigs();

  int i, j, k, l, p;

  for (i = 1; i <= GetNT(); i++)
    {
      if (IsMarkedTrig(i)) 
	{
	  for (j = 1; j <= 3; j++)
	    {
	      double origbadness = CalcTrigBadness(i);

	      p = GetTriangle(i).PNum(j);
	      Point<3> pm(0.,0.,0.); //Middlevector;
	      Point<3> p0(0.,0.,0.);

	      int cnt = 0;

	      for (k = 1; k <= trigsperpoint.EntrySize(p); k++)
		{
		  const STLTriangle& tr = GetTriangle(trigsperpoint.Get(p,k));
		  for (l = 1; l <= 3; l++)
		    {
		      if (tr.PNum(l) != p)
			{
			  cnt++;
			  pm(0) += GetPoint(tr.PNum(l))(0);
			  pm(1) += GetPoint(tr.PNum(l))(1);
			  pm(2) += GetPoint(tr.PNum(l))(2);
			}
		    }
		}
	      Point3d origp = GetPoint(p);
	      Point3d newp = p0 + fact*(1./(double)cnt)*(pm-p0)+(1.-fact)*(origp-p0);

	      SetPoint(p, newp);

	      if (CalcTrigBadness(i) > 0.9*origbadness) {SetPoint(p,origp); PrintDot('f');}
	      else {PrintDot('s');}
	    }
	}
    }
  MarkRevertedTrigs();
}

void STLGeometry :: MarkRevertedTrigs()
{
  int i,j;
  if (edgesperpoint.Size() != GetNP()) {BuildEdges();}

  PrintFnStart("mark reverted trigs");

  InitMarkedTrigs();

  int found;
  double revertedangle = stldoctor.smoothangle/180.*M_PI;

  int cnt = 0;
  int ap1, ap2;
  for (i = 1; i <= GetNT(); i++)
    {
      found = 0;
      for (j = 1; j <= NONeighbourTrigs(i); j++)
	{
	  GetTriangle(i).GetNeighbourPoints(GetTriangle(NeighbourTrig(i,j)), ap1, ap2);

	  if (!IsEdge(ap1,ap2))
	    {
              if (GetGeomAngle(i, NeighbourTrig(i,j)) > revertedangle)
		{
		  found = 1;
		  break;
		}
	    }
	}
      
      if (found)
	{
	  SetMarkedTrig(i, 1); cnt++;
	}
      
    }

  PrintMessage(5, "found ", cnt, " reverted trigs");


}

void STLGeometry :: SmoothDirtyTrigs()
{
  PrintFnStart("smooth dirty trigs");

  MarkDirtyTrigs();

  int i,j;
  int changed = 1;
  int ap1, ap2;
  
  while (changed)
    {
      changed = 0;
      for (i = 1; i <= GetNT(); i++)
	{
	  if (IsMarkedTrig(i))
	    {
	      int foundtrig = 0;
	      double maxlen = 0;
	      // JS: darf normalvector nicht ueber kurze Seite erben
	      maxlen = GetTriangle(i).MaxLength(GetPoints()) / 2.1; //JG: bei flachem dreieck auch kurze Seite

	      for (j = 1; j <= NONeighbourTrigs(i); j++)
		{
		  if (!IsMarkedTrig(NeighbourTrig(i,j)))
		    {
		      GetTriangle(i).GetNeighbourPoints(GetTriangle(NeighbourTrig(i,j)),ap1,ap2);
		      if (Dist(GetPoint(ap1),GetPoint(ap2)) >= maxlen)
			{
			  foundtrig = NeighbourTrig(i,j);
			  maxlen = Dist(GetPoint(ap1),GetPoint(ap2));
			}
		    }
		}
	      if (foundtrig)
		{
		  GetTriangle(i).SetNormal(GetTriangle(foundtrig).Normal());
		  changed = 1;
		  SetMarkedTrig(i,0);
		}
	    }
	}
    }

  calcedgedataanglesnew = 1;


  MarkDirtyTrigs();

  int cnt = 0;
  for (i = 1; i <= GetNT(); i++)
    {
      if (IsMarkedTrig(i)) {cnt++;}
    }

  PrintMessage(5,"NO marked dirty trigs=", cnt);

}

int STLGeometry :: IsMarkedTrig(int trig) const 
{
  if (trig <= markedtrigs.Size() && trig >=1)
    {
      return markedtrigs.Get(trig);
    }
  else {PrintSysError("In STLGeometry::IsMarkedTrig");}

  return 0;  
}

void STLGeometry :: SetMarkedTrig(int trig, int num)
{
  if (trig <= markedtrigs.Size() && trig >=1)
    {
      markedtrigs.Elem(trig) = num;
    }
  else {PrintSysError("In STLGeometry::SetMarkedTrig");}
}

void STLGeometry :: Clear()
{
  PrintFnStart("Clear");

  surfacemeshed = 0;
  surfaceoptimized = 0;
  volumemeshed = 0;

  selectedmultiedge.SetSize(0);
  meshlines.SetSize(0);
  // neighbourtrigs.SetSize(0);
  outerchartspertrig.SetSize(0);
  atlas.SetSize(0);
  ClearMarkedSegs();
  ClearSpiralPoints();
  ClearLineEndPoints();

  SetSelectTrig(0);
  SetNodeOfSelTrig(1);
  facecnt = 0;

  SetThreadPercent(100.);

  ClearEdges();
}

double STLGeometry :: Area()
{
  double ar = 0;
  int i;
  for (i = 1; i <= GetNT(); i++)
    {
      ar += GetTriangle(i).Area(points);
    }
  return ar;
}

double STLGeometry :: GetAngle(int t1, int t2)
{
  return Angle(GetTriangle(t1).Normal(),GetTriangle(t2).Normal());
}

double STLGeometry :: GetGeomAngle(int t1, int t2)
{
  Vec3d n1 = GetTriangle(t1).GeomNormal(points);
  Vec3d n2 = GetTriangle(t2).GeomNormal(points);
  return Angle(n1,n2);
}


void STLGeometry :: InitSTLGeometry(const ARRAY<STLReadTriangle> & readtrias)
{
  PrintFnStart("Init STL Geometry");
  STLTopology::InitSTLGeometry(readtrias);

  int i, k;

  //const double geometry_tol_fact = 1E8; //distances lower than max_box_size/tol are ignored

  int np = GetNP();
  PrintMessage(5,"NO points= ", GetNP());
  normals.SetSize(GetNP());
  ARRAY<int> normal_cnt(GetNP()); // counts number of added normals in a point

  for (i = 1; i <= np; i++)
    {
      normal_cnt.Elem(i) = 0;
      normals.Elem(i) = Vec3d (0,0,0);
    }

  for(i = 1; i <= GetNT(); i++)
    {
      //      STLReadTriangle t = GetReadTriangle(i);
      //      STLTriangle st;

      Vec<3> n = GetTriangle(i).Normal ();

      for (k = 1; k <= 3; k++)
	{
	  int pi = GetTriangle(i).PNum(k);
	  
	  normal_cnt.Elem(pi)++;
	  SetNormal(pi, GetNormal(pi) + n);
	}
    } 

  //normalize the normals
  for (i = 1; i <= GetNP(); i++)
    {
      SetNormal(i,1./(double)normal_cnt.Get(i)*GetNormal(i));
    }

  trigsconverted = 1;

  vicinity.SetSize(GetNT());
  markedtrigs.SetSize(GetNT());
  for (i = 1; i <= GetNT(); i++)
    {
      markedtrigs.Elem(i) = 0;
      vicinity.Elem(i) = 1;
    }

  ha_points.SetSize(GetNP());
  for (i = 1; i <= GetNP(); i++)
    ha_points.Elem(i) = 0;

  calcedgedataanglesnew = 0;
  edgedatastored = 0;
  edgedata->Clear();


  if (GetStatus() == STL_ERROR) return;

  CalcEdgeData();
  CalcEdgeDataAngles();

  ClearLineEndPoints();

  CheckGeometryOverlapping();
}

void STLGeometry :: TopologyChanged()
{
  calcedgedataanglesnew = 1;
}

int STLGeometry :: CheckGeometryOverlapping()
{
  int i, j, k;

  Box<3> geombox = GetBoundingBox();
  Point<3> pmin = geombox.PMin();
  Point<3> pmax = geombox.PMax();

  Box3dTree setree(pmin, pmax);
  ARRAY<int> inters;

  int oltrigs = 0;
  markedtrigs.SetSize(GetNT());

  for (i = 1; i <= GetNT(); i++)
    SetMarkedTrig(i, 0);

  for (i = 1; i <= GetNT(); i++)
    {
      const STLTriangle & tri = GetTriangle(i);
      
      Point<3> tpmin = tri.box.PMin();
      Point<3> tpmax = tri.box.PMax();
      Vec<3> diag = tpmax - tpmin;

      tpmax = tpmax + 0.001 * diag;
      tpmin = tpmin - 0.001 * diag;

      setree.Insert (tpmin, tpmax, i);
    }

  for (i = 1; i <= GetNT(); i++)
    {
      const STLTriangle & tri = GetTriangle(i);
      
      Point<3> tpmin = tri.box.PMin();
      Point<3> tpmax = tri.box.PMax();

      setree.GetIntersecting (tpmin, tpmax, inters);

      for (j = 1; j <= inters.Size(); j++)
	{
	  const STLTriangle & tri2 = GetTriangle(inters.Get(j));

	  const Point<3> *trip1[3], *trip2[3];	
	  Point<3> hptri1[3], hptri2[3];
	  /*
	  for (k = 1; k <= 3; k++)
	    {
	      trip1[k-1] = &GetPoint (tri.PNum(k));
	      trip2[k-1] = &GetPoint (tri2.PNum(k));
	    }
	  */

	  for (k = 0; k < 3; k++)
	    {
	      hptri1[k] = GetPoint (tri[k]);
	      hptri2[k] = GetPoint (tri2[k]);
	      trip1[k] = &hptri1[k];
	      trip2[k] = &hptri2[k];
	    }

	  if (IntersectTriangleTriangle (&trip1[0], &trip2[0]))
	    {
	      oltrigs++;
	      PrintMessage(5,"Intersecting Triangles: trig ",i," with ",inters.Get(j),"!");
	      SetMarkedTrig(i, 1);
	      SetMarkedTrig(inters.Get(j), 1);
	    }
	}
    }

  PrintMessage(3,"Check Geometry Overlapping: overlapping triangles = ",oltrigs);
  return oltrigs;
}

/*
void STLGeometry :: InitSTLGeometry()
{
  STLTopology::InitSTLGeometry();

  int i, j, k;

  const double geometry_tol_fact = 1E8; //distances lower than max_box_size/tol are ignored


  trias.SetSize(0);
  points.SetSize(0);
  normals.SetSize(0);

  ARRAY<int> normal_cnt; // counts number of added normals in a point

  Box3d bb (GetBoundingBox().PMin() + Vec3d (-1,-1,-1),
	    GetBoundingBox().PMax() + Vec3d (1, 1, 1));

  Point3dTree pointtree (bb.PMin(), 
			 bb.PMax());
  ARRAY<int> pintersect;

  double gtol = GetBoundingBox().CalcDiam()/geometry_tol_fact;

  for(i = 1; i <= GetReadNT(); i++)
    {
      //if (i%500==499) {(*mycout) << (double)i/(double)GetReadNT()*100. << "%" << endl;}

      STLReadTriangle t = GetReadTriangle(i);
      STLTriangle st;
      Vec3d n = t.normal;

      for (k = 0; k < 3; k++)
	{
	  Point3d p = t.pts[k];

	  Point3d pmin = p - Vec3d (gtol, gtol, gtol);
	  Point3d pmax = p + Vec3d (gtol, gtol, gtol);
	  
	  pointtree.GetIntersecting (pmin, pmax, pintersect);
	  
	  if (pintersect.Size() > 1)
	    (*mycout) << "found too much  " << char(7) << endl;
	  int foundpos = 0;
	  if (pintersect.Size())
	    foundpos = pintersect.Get(1);

	  if (foundpos) 
	    {
	      normal_cnt[foundpos]++;
	      SetNormal(foundpos,GetNormal(foundpos)+n);
	      //	      (*testout) << "found p " << p << endl;
	    }
	  else
	    {
	      foundpos = AddPoint(p);
	      AddNormal(n);
	      normal_cnt.Append(1);

	      pointtree.Insert (p, foundpos);
	    }
	  //(*mycout) << "foundpos=" << foundpos << endl;
	  st.pts[k] = foundpos;
	}

      if ( (st.pts[0] == st.pts[1]) || 
	   (st.pts[0] == st.pts[2]) || 
	   (st.pts[1] == st.pts[2]) )
	{
	  (*mycout) << "ERROR: STL Triangle degenerated" << endl;
	}
      else
	{
	  // do not add ? js
	  AddTriangle(st);
	}
      //(*mycout) << "TRIG" << i << " = " << st << endl;
      
    } 
  //normal the normals
  for (i = 1; i <= GetNP(); i++)
    {
      SetNormal(i,1./(double)normal_cnt[i]*GetNormal(i));
    }

  trigsconverted = 1;

  vicinity.SetSize(GetNT());
  markedtrigs.SetSize(GetNT());
  for (i = 1; i <= GetNT(); i++)
    {
      markedtrigs.Elem(i) = 0;
      vicinity.Elem(i) = 1;
    }

  ha_points.SetSize(GetNP());
  for (i = 1; i <= GetNP(); i++)
    ha_points.Elem(i) = 0;

  calcedgedataanglesnew = 0;
  edgedatastored = 0;
  edgedata->Clear();

  CalcEdgeData();
  CalcEdgeDataAngles();

  ClearLineEndPoints();

  (*mycout) << "done" << endl;
}
*/



void STLGeometry :: SetLineEndPoint(int pn) 
{
  if (pn <1 || pn > lineendpoints.Size()) {PrintSysError("Illegal pnum in SetLineEndPoint!!!"); return; }
  lineendpoints.Elem(pn) = 1;
}

int STLGeometry :: IsLineEndPoint(int pn) 
{
  //  return 0;
  if (pn <1 || pn > lineendpoints.Size()) 
    {PrintSysError("Illegal pnum in IsLineEndPoint!!!"); return 0;}
  return lineendpoints.Get(pn);
}

void STLGeometry :: ClearLineEndPoints()
{
  lineendpoints.SetSize(GetNP());
  int i;
  for (i = 1; i <= GetNP(); i++)
    {
      lineendpoints.Elem(i) = 0;
    }
}

int STLGeometry :: IsEdge(int ap1, int ap2)
{
  int i,j;
  for (i = 1; i <= GetNEPP(ap1); i++)
    {
      for (j = 1; j <= GetNEPP(ap2); j++)
	{
	  if (GetEdgePP(ap1,i) == GetEdgePP(ap2,j)) {return 1;}
	}
    }
  return 0;
}

int STLGeometry :: IsEdgeNum(int ap1, int ap2)
{
  int i,j;
  for (i = 1; i <= GetNEPP(ap1); i++)
    {
      for (j = 1; j <= GetNEPP(ap2); j++)
	{
	  if (GetEdgePP(ap1,i) == GetEdgePP(ap2,j)) {return GetEdgePP(ap1,i);}
	}
    }
  return 0;
}


void STLGeometry :: BuildEdges()
{
  //PrintFnStart("build edges");
  edges.SetSize(0);
  meshlines.SetSize(0);
  FindEdgesFromAngles();
}

void STLGeometry :: UseExternalEdges()
{
  int i;
  for (i = 1; i <= NOExternalEdges(); i++)
    {
      AddEdge(GetExternalEdge(i).i1,GetExternalEdge(i).i2);
    }
  //BuildEdgesPerPointy();
}

void STLGeometry :: UndoEdgeChange()
{
  if (edgedatastored) 
    {
      RestoreEdgeData();
    }
  else
    {
      PrintWarning("no edge undo possible");
    }
}


void STLGeometry :: StoreEdgeData()
{
  //  edgedata_store = *edgedata;
  
  edgedata->Store();
  edgedatastored = 1;

  // put stlgeom-edgedata to stltopology edgedata 
  /*
  int i;
  for (i = 1; i <= GetNTE(); i++)
    {
      const STLTopEdge & topedge = GetTopEdge (i);
      int ednum = edgedata->GetEdgeNum (topedge.PNum(1),
				       topedge.PNum(2));
      topedges.Elem(i).SetStatus (edgedata->Get (ednum).status);
    }
  */
}

void STLGeometry :: RestoreEdgeData()
{
  //  *edgedata = edgedata_store;
  edgedata->Restore();
  edgedatastored=0;
}


void STLGeometry :: CalcEdgeData()
{
  PushStatus("Calc Edge Data");

  int np1, np2;
  int i;
  
  int ecnt = 0;
  edgedata->SetSize(GetNT()/2*3);

  for (i = 1; i <= GetNT(); i++)
    {
      SetThreadPercent((double)i/(double)GetNT()*100.);
      
      const STLTriangle & t1 = GetTriangle(i);

      for (int j = 1; j <= NONeighbourTrigs(i); j++)
	{
	  int nbti = NeighbourTrig(i,j);
	  if (nbti > i)
	    {
	      const STLTriangle & t2 = GetTriangle(nbti);

	      if (t1.IsNeighbourFrom(t2))
		{
		  ecnt++; if (ecnt > edgedata->Size()) {PrintError("In Calc edge data, illegal geometry");}

		  t1.GetNeighbourPoints(t2,np1,np2);

		  /* ang = GetAngle(i,nbti);
		     if (ang < -M_PI) {ang += 2*M_PI;}*/


		  // edgedata->Add(STLEdgeData(0, np1, np2, i, nbti),ecnt);
		  edgedata->Elem(ecnt).SetStatus(ED_UNDEFINED);

		  // edgedata->Elem(ecnt).top = this;
		  // edgedata->Elem(ecnt).topedgenr = GetTopEdgeNum (np1, np2);
		}
	    }
	}      
    }
  
  //BuildEdgesPerPoint();
  PopStatus();  
}

void STLGeometry :: CalcEdgeDataAngles()
{
  PrintMessage(5,"calc edge data angles");

  int i;

  for (i = 1; i <= GetNTE(); i++)
    {
      STLTopEdge & edge = GetTopEdge (i);
      double cosang = 
	GetTriangle(edge.TrigNum(1)).Normal() *
	GetTriangle(edge.TrigNum(2)).Normal();
      edge.SetCosAngle (cosang);
    }

  for (i = 1; i <= edgedata->Size(); i++)
    {
      /*
      const STLEdgeData& e = edgedata->Get(i);
      ang = GetAngle(e.lt,e.rt);
      if (ang < -M_PI) {ang += 2*M_PI;}
      edgedata->Elem(i).angle = fabs(ang);
      */
    }
  
}

void STLGeometry :: FindEdgesFromAngles()
{
  //  PrintFnStart("find edges from angles");

  double min_edge_angle = stlparam.yangle/180.*M_PI;
  double cont_min_edge_angle = stlparam.contyangle/180.*M_PI;

  double cos_min_edge_angle = cos (min_edge_angle);
  double cos_cont_min_edge_angle = cos (cont_min_edge_angle);

  if (calcedgedataanglesnew) {CalcEdgeDataAngles(); calcedgedataanglesnew = 0;}

  int i;
  for (i = 1; i <= edgedata->Size(); i++)
    {
      STLTopEdge & sed = edgedata->Elem(i);
      if (sed.GetStatus() == ED_CANDIDATE || 
	  sed.GetStatus() == ED_UNDEFINED)
	{
	  if (sed.CosAngle() <= cos_min_edge_angle)
	    {
	      sed.SetStatus (ED_CANDIDATE);
	    }
	  else
	    {
	      sed.SetStatus(ED_UNDEFINED);
	    }
	} 
    }

  if (stlparam.contyangle < stlparam.yangle)
    {
      int changed = 1;
      int its = 0;
      while (changed && stlparam.contyangle < stlparam.yangle)
	{
	  its++;
	  //(*mycout) << "." << flush;
	  changed = 0;
	  for (i = 1; i <= edgedata->Size(); i++)
	    {
	      STLTopEdge & sed = edgedata->Elem(i);
	      if (sed.CosAngle() <= cos_cont_min_edge_angle 
		  && sed.GetStatus() == ED_UNDEFINED && 
		  (edgedata->GetNConfCandEPP(sed.PNum(1)) == 1 || 
		   edgedata->GetNConfCandEPP(sed.PNum(2)) == 1))
		{
		  changed = 1;
		  sed.SetStatus (ED_CANDIDATE);
		}
	    }
	}
    }
  
  int confcand = 0;
  if (edgedata->GetNConfEdges() == 0) 
    {
      confcand = 1;
    }
  
  for (i = 1; i <= edgedata->Size(); i++)
    {
      STLTopEdge & sed = edgedata->Elem(i);
      if (sed.GetStatus() == ED_CONFIRMED || 
	  (sed.GetStatus() == ED_CANDIDATE && confcand))
	{
	  STLEdge se(sed.PNum(1),sed.PNum(2));
	  se.SetLeftTrig(sed.TrigNum(1));
	  se.SetRightTrig(sed.TrigNum(2));
	  AddEdge(se);
	}
    }
  BuildEdgesPerPoint();

  

  //(*mycout) << "its for continued angle = " << its << endl;
  PrintMessage(5,"built ", GetNE(), " edges with yellow angle = ", stlparam.yangle, " degree");
  
}

/*
void STLGeometry :: FindEdgesFromAngles()
{
  double yangle = stlparam.yangle;
  char * savetask = multithread.task;
  multithread.task = "find edges";

  const double min_edge_angle = yangle/180.*M_PI;

  int np1, np2;
  double ang;
  int i;

  //(*mycout) << "area=" << Area() << endl;

  for (i = 1; i <= GetNT(); i++)
    {
      multithread.percent = (double)i/(double)GetReadNT()*100.;
      
      const STLTriangle & t1 = GetTriangle(i);
      //NeighbourTrigs(nt,i);

      for (int j = 1; j <= NONeighbourTrigs(i); j++)
	{
	  int nbti = NeighbourTrig(i,j);
	  if (nbti > i)
	    {
	      const STLTriangle & t2 = GetTriangle(nbti);

	      if (t1.IsNeighbourFrom(t2))
		{
		  ang = GetAngle(i,nbti);
		  if (ang < -M_PI*0.5) {ang += 2*M_PI;}

		  t1.GetNeighbourPoints(t2,np1,np2);
		  
		  if (fabs(ang) >= min_edge_angle)
		    {
		      STLEdge se(np1,np2);
		      se.SetLeftTrig(i);
		      se.SetRightTrig(nbti);
		      AddEdge(se);
		    }
		}
	    }
	}      
    }
  
  (*mycout) << "added " << GetNE() << " edges" << endl;

  //BuildEdgesPerPoint();

  multithread.percent = 100.;
  multithread.task = savetask;
  
}
*/
void STLGeometry :: BuildEdgesPerPoint()
{
  //cout << "*** build edges per point" << endl;
  edgesperpoint.SetSize(GetNP());

  //add edges to points
  int i;
  for (i = 1; i <= GetNE(); i++)
    {
      //(*mycout) << "EDGE " << GetEdge(i).PNum(1) << " - " << GetEdge(i).PNum(2) << endl;
      for (int j = 1; j <= 2; j++)
	{
	  AddEdgePP(GetEdge(i).PNum(j),i);
	}
    }
}

void STLGeometry :: AddFaceEdges()
{
  PrintFnStart("Add starting edges for faces");

  //für Kugel eine STLLine hinzufügen (Vorteil: verfeinerbar, unabhängig von Auflösung der Geometrie!!!):
  //Grenze von 1. gefundener chart

  ARRAY<int> edgecnt;
  ARRAY<int> chartindex;
  edgecnt.SetSize(GetNOFaces());
  chartindex.SetSize(GetNOFaces());

  int i,j;
  for (i = 1; i <= GetNOFaces(); i++)
    {
      edgecnt.Elem(i) = 0;
      chartindex.Elem(i) = 0;
    }

  for (i = 1; i <= GetNT(); i++)
    {
      int fn = GetTriangle(i).GetFaceNum();
      if (!chartindex.Get(fn)) {chartindex.Elem(fn) = GetChartNr(i);}
      for (j = 1; j <= 3; j++)
	{
	  edgecnt.Elem(fn) += GetNEPP(GetTriangle(i).PNum(j));
	}
    }

  for (i = 1; i <= GetNOFaces(); i++)
    {
      if (!edgecnt.Get(i)) {PrintMessage(5,"Face", i, " has no edge!");}
    }
  
  int changed = 0;
  int k, ap1, ap2;
  for (i = 1; i <= GetNOFaces(); i++)
    {
      if (!edgecnt.Get(i))
      {
	const STLChart& c = GetChart(chartindex.Get(i));
	for (j = 1; j <= c.GetNChartT(); j++)
	  {
	    const STLTriangle& t1 = GetTriangle(c.GetChartTrig(j));
	    for (k = 1; k <= 3; k++)
	      {
		int nt = NeighbourTrig(c.GetChartTrig(j),k);
		if (GetChartNr(nt) != chartindex.Get(i))
		  {
		    t1.GetNeighbourPoints(GetTriangle(nt),ap1,ap2);
		    AddEdge(ap1,ap2);
		    changed = 1;
		  }
	      }
	  }
      }
      
    }
  
  if (changed) BuildEdgesPerPoint();
  
}

void STLGeometry :: LinkEdges()
{
  PushStatusF("Link Edges");
  PrintMessage(5,"have now ", GetNE(), " edges with yellow angle = ", stlparam.yangle, " degree");

  int i;

  lines.SetSize(0);
  int starte(0);
  int edgecnt = 0;
  int found;
  int rev(0); //indicates, that edge is inserted reverse

  //worked edges
  ARRAY<int> we(GetNE());

  //setlineendpoints; wenn 180°, dann keine endpunkte
  //nur punkte mit 2 edges kommen in frage, da bei mehr oder weniger punkten ohnehin ein meshpoint hinkommt

  Vec3d v1,v2;
  double cos_eca = cos(stlparam.edgecornerangle/180.*M_PI);
  int ecnt = 0;
  int lp1, lp2;
  if (stlparam.edgecornerangle < 180)
    {
      for (i = 1; i <= GetNP(); i++)
	{
	  if (GetNEPP(i) == 2)
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
	      if ((v1*v2)/sqrt(v1.Length2()*v2.Length2()) < cos_eca) 
		{
		  //(*testout) << "add edgepoint " << i << endl;
		  SetLineEndPoint(i);
		  ecnt++;
		}
	    }	  
	}
    }
  PrintMessage(5, "added ", ecnt, " mesh_points due to edge corner angle (", 
	       stlparam.edgecornerangle, " degree)");

  for (i = 1; i <= GetNE(); i++) {we.Elem(i) = 0;}

  while(edgecnt < GetNE())
    {
      SetThreadPercent((double)edgecnt/(double)GetNE()*100.);

      STLLine* line = new STLLine(this);

      //find start edge
      int j = 1;
      found = 0;
      //try second time, if only rings are left!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      int second = 0;

      //find a starting edge at point with 1 or more than 2 edges or at lineendpoint
      while (!found && j<=GetNE())
	{
	  if (!we.Get(j))
	    {
	      if (GetNEPP(GetEdge(j).PNum(1)) != 2 || IsLineEndPoint(GetEdge(j).PNum(1)))
		{
		  starte = j;
		  found = 1;
		  rev = 0;
		}
	      else 
	      if (GetNEPP(GetEdge(j).PNum(2)) != 2 || IsLineEndPoint(GetEdge(j).PNum(2)))
		{
		  starte = j;
		  found = 1;
		  rev = 1;
		}
	      else if (second)
		{
		  starte = j;
		  found = 1;
		  rev = 0; //0 or 1 are possible
		}
	    }
	  j++;
	  if (!second && j == GetNE()) {second = 1; j = 1;}
	}

      if (!found) {PrintSysError("No starting edge found, edgecnt=", edgecnt, ", GETNE=", GetNE());}

      line->AddPoint(GetEdge(starte).PNum(1+rev));
      line->AddPoint(GetEdge(starte).PNum(2-rev));
      if (!rev)
	{
	  line->AddLeftTrig(GetEdge(starte).LeftTrig());
	  line->AddRightTrig(GetEdge(starte).RightTrig());
	}
      else
	{
	  line->AddLeftTrig(GetEdge(starte).RightTrig());
	  line->AddRightTrig(GetEdge(starte).LeftTrig());
	}
      edgecnt++; we.Elem(starte) = 1;

      //add segments to line as long as segments other than starting edge are found or lineendpoint is reached 
      found = 1;
      int other;
      while(found)
	{
	  found = 0;
	  int fp = GetEdge(starte).PNum(2-rev);
	  if (GetNEPP(fp) == 2 && !IsLineEndPoint(fp))
	    {
	      //find the "other" edge of point fp
	      other = 0;
	      if (GetEdgePP(fp,1) == starte) {other = 1;}

	      starte = GetEdgePP(fp,1+other);

	      //falls ring -> aufhoeren !!!!!!!!!!!
	      if (!we.Elem(starte))
		{
		  found = 1;
		  rev = 0;
		  if (GetEdge(starte).PNum(2) == fp) {rev = 1;}
		  else if (GetEdge(starte).PNum(1) != fp) {PrintSysError("In Link Edges!");}

		  line->AddPoint(GetEdge(starte).PNum(2-rev));
		  if (!rev) 
		    {
		      line->AddLeftTrig(GetEdge(starte).LeftTrig());
		      line->AddRightTrig(GetEdge(starte).RightTrig());
		    }
		  else
		    {
		      line->AddLeftTrig(GetEdge(starte).RightTrig());
		      line->AddRightTrig(GetEdge(starte).LeftTrig());
		    }
		  edgecnt++; we.Elem(starte) = 1;
		}
	    }     
	}
      AddLine(line);      
    }
  PrintMessage(5,"number of lines generated = ", GetNLines());

  //check, which lines must have at least one midpoint
  INDEX_2_HASHTABLE<int> lineht(GetNLines()+1);

  for (i = 1; i <= GetNLines(); i++)
    {
      if (GetLine(i)->StartP() == GetLine(i)->EndP())
	{
	  GetLine(i)->DoSplit();	  
	}
    }

  for (i = 1; i <= GetNLines(); i++)
    {
      INDEX_2 lineep (GetLine(i)->StartP(),GetLine(i)->EndP());
      lineep.Sort();

      if (lineht.Used (lineep))
	{
	  GetLine(i)->DoSplit();
	  int other = lineht.Get(lineep);
	  GetLine(other)->DoSplit();
	}
      else
	{
	  lineht.Set (lineep, i);
	}
    }

  for (i = 1; i <= GetNLines(); i++)
    {
      STLLine* line = GetLine(i);
      for (int ii = 1; ii <= line->GetNS(); ii++)
	{
	  int ap1, ap2;
	  line->GetSeg(ii,ap1,ap2);
	  //	  (*mycout) << "SEG " << p1 << " - " << p2 << endl;
	}
    }

  PopStatus();
}

int STLGeometry :: GetNOBodys()
{
  int markedtrigs1 = 0;
  int starttrig = 1;
  int i, k, nnt;
  int bodycnt = 0;

  ARRAY<int> bodynum(GetNT());

  for (i = 1; i <= GetNT(); i++)
    bodynum.Elem(i)=0;


  while (markedtrigs1 < GetNT())
    {
      for (i = starttrig; i <= GetNT(); i++)
	{
	  if (!bodynum.Get(i))
	    {
	      starttrig = i;
	      break;
	    }
	} 
      //add all triangles around starttriangle, which is reachable without going over an edge
      ARRAY<int> todolist;
      ARRAY<int> nextlist;
      bodycnt++;
      markedtrigs1++;
      bodynum.Elem(starttrig) = bodycnt;
      todolist.Append(starttrig);

      while(todolist.Size())
	{
	  for (i = 1; i <= todolist.Size(); i++)
	    {
	      //const STLTriangle& tt = GetTriangle(todolist.Get(i));
	      for (k = 1; k <= NONeighbourTrigs(todolist.Get(i)); k++)
		{
		  nnt = NeighbourTrig(todolist.Get(i),k);
		  if (!bodynum.Get(nnt))
		    {
		      nextlist.Append(nnt);
		      bodynum.Elem(nnt) = bodycnt;
		      markedtrigs1++;
		    }
		}
	    }
	  
	  todolist.SetSize(0);
	  for (i = 1; i <= nextlist.Size(); i++)
	    {
	      todolist.Append(nextlist.Get(i));
	    }
	  nextlist.SetSize(0);	  
	}
    }
  PrintMessage(3, "Geometry has ", bodycnt, " separated bodys");

  return bodycnt;
}

void STLGeometry :: CalcFaceNums()
{
  int markedtrigs1 = 0;
  int starttrig(0);
  int laststarttrig = 1;
  int i, k, nnt;
  facecnt = 0;


  for (i = 1; i <= GetNT(); i++)
    GetTriangle(i).SetFaceNum(0);


  while (markedtrigs1 < GetNT())
    {
      for (i = laststarttrig; i <= GetNT(); i++)
	{
	  if (!GetTriangle(i).GetFaceNum()) 
	    {
	      starttrig = i;
	      laststarttrig = i;
	      break;
	    }
	} 
      //add all triangles around starttriangle, which is reachable without going over an edge
      ARRAY<int> todolist;
      ARRAY<int> nextlist;
      facecnt++;
      markedtrigs1++;
      GetTriangle(starttrig).SetFaceNum(facecnt);
      todolist.Append(starttrig);
      int ap1, ap2;

      while(todolist.Size())
	{
	  for (i = 1; i <= todolist.Size(); i++)
	    {
	      const STLTriangle& tt = GetTriangle(todolist.Get(i));
	      for (k = 1; k <= NONeighbourTrigs(todolist.Get(i)); k++)
		{
		  nnt = NeighbourTrig(todolist.Get(i),k);
		  STLTriangle& nt = GetTriangle(nnt);
		  if (!nt.GetFaceNum())
		    {
		      tt.GetNeighbourPoints(nt,ap1,ap2);
		      if (!IsEdge(ap1,ap2))
			{
			  nextlist.Append(nnt);
			  nt.SetFaceNum(facecnt);
			  markedtrigs1++;
			}
		    }
		}
	    }
	  
	  todolist.SetSize(0);
	  for (i = 1; i <= nextlist.Size(); i++)
	    {
	      todolist.Append(nextlist.Get(i));
	    }
	  nextlist.SetSize(0);	  
	}
    }
  GetNOBodys();
  PrintMessage(3,"generated ", facecnt, " faces");
}
 
void STLGeometry :: ClearSpiralPoints()
{
  spiralpoints.SetSize(GetNP());
  int i;
  for (i = 1; i <= spiralpoints.Size(); i++)
    {
      spiralpoints.Elem(i) = 0;
    }
}


void STLGeometry :: BuildSmoothEdges ()
{
  if (smoothedges) delete smoothedges;

  smoothedges = new INDEX_2_HASHTABLE<int> (GetNE()/10 + 1);


  // Jack: Ok ?
  //  UseExternalEdges();

  PushStatusF("Build Smooth Edges");

  int i, j;//, k, l;
  int nt = GetNT();
  Vec3d ng1, ng2;

  for (i = 1; i <= nt; i++)
    {
      if (multithread.terminate)
	{PopStatus();return;}

      SetThreadPercent(100.0 * (double)i / (double)nt);

      const STLTriangle & trig = GetTriangle (i);
      
      ng1 = trig.GeomNormal(points);
      ng1 /= (ng1.Length() + 1e-24);

      for (j = 1; j <= 3; j++)
	{ 
	  int nbt = NeighbourTrig (i, j);
	  
	  ng2 = GetTriangle(nbt).GeomNormal(points);
	  ng2 /= (ng2.Length() + 1e-24);
	  
	  
	  int pi1, pi2;

	  trig.GetNeighbourPoints(GetTriangle(nbt), pi1, pi2);

	  if (!IsEdge(pi1,pi2)) 
	    {
	      if (ng1 * ng2 < 0)
		{
		  PrintMessage(7,"smoothedge found");
		  INDEX_2 i2(pi1, pi2);
		  i2.Sort();
		  smoothedges->Set (i2, 1);
		}
	    }
	}
    }

  PopStatus();
}





int STLGeometry :: IsSmoothEdge (int pi1, int pi2) const
{
  if (!smoothedges)
    return 0;
  INDEX_2 i2(pi1, pi2);
  i2.Sort();
  return smoothedges->Used (i2);
}




//function is not used now
int IsInArray(int n, const ARRAY<int>& ia)
{
  int i;
  for (i = 1; i <= ia.Size(); i++)
    {
      if (ia.Get(i) == n) {return 1;}
    }
  return 0;
}

void STLGeometry :: AddConeAndSpiralEdges()
{
  PrintMessage(5,"have now ", GetNE(), " edges with yellow angle = ", stlparam.yangle, " degree");

  PrintFnStart("AddConeAndSpiralEdges");

  int i,j,k,n;
  //  int changed = 0;

  //check edges, where inner chart and no outer chart come together without an edge
  int np1, np2, nt;
  int cnt = 0;

  for (i = 1; i <= GetNOCharts(); i++)
    {
      STLChart& chart = GetChart(i);
      for (j = 1; j <= chart.GetNChartT(); j++)
	{
	  int t = chart.GetChartTrig(j); 
	  const STLTriangle& tt = GetTriangle(t);

	  for (k = 1; k <= 3; k++)
	    {
	      nt = NeighbourTrig(t,k); 
	      if (GetChartNr(nt) != i && !TrigIsInOC(nt,i))
		{	      
		  tt.GetNeighbourPoints(GetTriangle(nt),np1,np2);
		  if (!IsEdge(np1,np2))
		    {
		      STLEdge se(np1,np2);
		      se.SetLeftTrig(t);
		      se.SetRightTrig(nt);
		      int edgenum = AddEdge(se);
		      AddEdgePP(np1,edgenum);
		      AddEdgePP(np2,edgenum);
		      //changed = 1;
		      PrintWarning("Found a spiral like structure: chart=", i,
				   ", trig=", t, ", p1=", np1, ", p2=", np2);
		      cnt++;
		    }
		}
	    }
	}
	  
    }

  PrintMessage(5, "found ", cnt, " spiral like structures");
  PrintMessage(5, "added ", cnt, " edges due to spiral like structures");
  
  cnt = 0;
  int edgecnt = 0;

  ARRAY<int> trigsaroundp;
  ARRAY<int> chartpointchecked; //gets number of chart, if in this chart already checked
  chartpointchecked.SetSize(GetNP());

  for (i = 1; i <= GetNP(); i++)
    {
      chartpointchecked.Elem(i) = 0;
    }

  int onoc, notonoc, tpp, pn;
  int ap1, ap2, tn1, tn2, l, problem;

  if (!stldoctor.conecheck) {PrintWarning("++++++++++++ \ncone checking deactivated by user!!!!!\n+++++++++++++++"); return ;}

  PushStatus("Find Critical Points");

  int addedges = 0;

  for (i = 1; i <= GetNOCharts(); i++)
    {
      SetThreadPercent((double)i/(double)GetNOCharts()*100.);
      if (multithread.terminate)
	{PopStatus();return;}

      STLChart& chart = GetChart(i);
      for (j = 1; j <= chart.GetNChartT(); j++)
	{
	  int t = chart.GetChartTrig(j); 
	  const STLTriangle& tt = GetTriangle(t);

	  for (k = 1; k <= 3; k++)
	    {
	      pn = tt.PNum(k);
	      if (chartpointchecked.Get(pn) == i)
		{continue;}
	      
	      int checkpoint = 0;
	      for (n = 1; n <= trigsperpoint.EntrySize(pn); n++)
		{
		  if (trigsperpoint.Get(pn,n) != t && 
		      GetChartNr(trigsperpoint.Get(pn,n)) != i &&
		      !TrigIsInOC(trigsperpoint.Get(pn,n),i)) {checkpoint = 1;};
		}
	      if (checkpoint)
		{
		  chartpointchecked.Elem(pn) = i;

		  int worked = 0;
		  int spworked = 0;
		  GetSortedTrianglesAroundPoint(pn,t,trigsaroundp);
		  trigsaroundp.Append(t);
		      
		  problem = 0;
		  for (l = 2; l <= trigsaroundp.Size()-1; l++)
		    {
		      tn1 = trigsaroundp.Get(l-1);
		      tn2 = trigsaroundp.Get(l);
		      const STLTriangle& t1 = GetTriangle(tn1);
		      const STLTriangle& t2 = GetTriangle(tn2);
		      t1.GetNeighbourPoints(t2, ap1, ap2);
		      if (IsEdge(ap1,ap2)) break;
		      
		      if (GetChartNr(tn2) != i && !TrigIsInOC(tn2,i)) {problem = 1;}
		    }

		  if (problem)
		    {
		      for (l = 2; l <= trigsaroundp.Size()-1; l++)
			{
			  tn1 = trigsaroundp.Get(l-1);
			  tn2 = trigsaroundp.Get(l);
			  const STLTriangle& t1 = GetTriangle(tn1);
			  const STLTriangle& t2 = GetTriangle(tn2);
			  t1.GetNeighbourPoints(t2, ap1, ap2);
			  if (IsEdge(ap1,ap2)) break;
			  
			  if ((GetChartNr(tn1) == i && GetChartNr(tn2) != i && TrigIsInOC(tn2,i)) ||
			      (GetChartNr(tn2) == i && GetChartNr(tn1) != i && TrigIsInOC(tn1,i))) 				 
			    {
			      if (addedges || !GetNEPP(pn))
				{
				  STLEdge se(ap1,ap2);
				  se.SetLeftTrig(tn1);
				  se.SetRightTrig(tn2);
				  int edgenum = AddEdge(se);
				  AddEdgePP(ap1,edgenum);
				  AddEdgePP(ap2,edgenum);
				  edgecnt++;
				}
			      if (!addedges && !GetSpiralPoint(pn))
				{
				  SetSpiralPoint(pn);
				  spworked = 1;
				}
			      worked = 1;
			    }
			}
		    }
		  //backwards:
		  problem = 0;
		  for (l = trigsaroundp.Size()-1; l >= 2; l--)
		    {
		      tn1 = trigsaroundp.Get(l+1);
		      tn2 = trigsaroundp.Get(l);
		      const STLTriangle& t1 = GetTriangle(tn1);
		      const STLTriangle& t2 = GetTriangle(tn2);
		      t1.GetNeighbourPoints(t2, ap1, ap2);
		      if (IsEdge(ap1,ap2)) break;
		      
		      if (GetChartNr(tn2) != i && !TrigIsInOC(tn2,i)) {problem = 1;}
		    }
		  if (problem)
		    for (l = trigsaroundp.Size()-1; l >= 2; l--)
		      {
			tn1 = trigsaroundp.Get(l+1);
			tn2 = trigsaroundp.Get(l);
			const STLTriangle& t1 = GetTriangle(tn1);
			const STLTriangle& t2 = GetTriangle(tn2);
			t1.GetNeighbourPoints(t2, ap1, ap2);
			if (IsEdge(ap1,ap2)) break;
			
			if ((GetChartNr(tn1) == i && GetChartNr(tn2) != i && TrigIsInOC(tn2,i)) ||
			    (GetChartNr(tn2) == i && GetChartNr(tn1) != i && TrigIsInOC(tn1,i))) 				 
			  {
			    if (addedges || !GetNEPP(pn))
			      {
				STLEdge se(ap1,ap2);
				se.SetLeftTrig(tn1);
				se.SetRightTrig(tn2);
				int edgenum = AddEdge(se);
				AddEdgePP(ap1,edgenum);
				AddEdgePP(ap2,edgenum);
				edgecnt++;
			      }
			    if (!addedges && !GetSpiralPoint(pn))
			      {
				SetSpiralPoint(pn);
				spworked = 1;
				//if (GetNEPP(pn) == 0) {(*mycout) << "ERROR: spiralpoint with no edge found!" << endl;}
			      }
			    worked = 1;
			  }
		      }

		  if (worked)
		    {		      
		      //(*testout) << "set edgepoint due to spirals: pn=" << i << endl;
		      SetLineEndPoint(pn);
		    }
		  if (spworked)
		    {		
		      /*      
		      (*mycout) << "Warning: Critical Point " << tt.PNum(k) 
			   << "( chart " << i << ", trig " << t
			   << ") has been neutralized!!!" << endl;
		      */
		      cnt++;
		    }
		  //		  markedpoints.Elem(tt.PNum(k)) = 1;
		}
	    }
	}
    }
  PrintMessage(5, "found ", cnt, " critical points!");
  PrintMessage(5, "added ", edgecnt, " edges due to critical points!");

  PopStatus();

  //search points where inner chart and outer chart and "no chart" trig come together at edge-point

  PrintMessage(7,"search for special chart points");
  for (i = 1; i <= GetNOCharts(); i++)
    {
      STLChart& chart = GetChart(i);
      for (j = 1; j <= chart.GetNChartT(); j++)
	{
	  int t = chart.GetChartTrig(j); 
	  const STLTriangle& tt = GetTriangle(t);

	  for (k = 1; k <= 3; k++)
	    {
	      pn = tt.PNum(k);
	      if (GetNEPP(pn) == 2)
		{
		  onoc = 0;
		  notonoc = 0;
		  for (n = 1; n <= trigsperpoint.EntrySize(pn); n++)
		    {
		      tpp = trigsperpoint.Get(pn,n);
		      if (tpp != t && GetChartNr(tpp) != i)
			{
			  if (TrigIsInOC(tpp,i)) {onoc = 1;}
			  if (!TrigIsInOC(tpp,i)) {notonoc = 1;}
			}
		    }
		  if (onoc && notonoc && !IsLineEndPoint(pn)) 
		    {
		      GetSortedTrianglesAroundPoint(pn,t,trigsaroundp);
		      int here = 1; //we start on this side of edge, !here = there
		      int thereOC = 0;
		      int thereNotOC = 0;
		      for (l = 2; l <= trigsaroundp.Size(); l++)
			{
			  GetTriangle(trigsaroundp.Get(l-1)).
			    GetNeighbourPoints(GetTriangle(trigsaroundp.Get(l)), ap1, ap2);
			  if (IsEdge(ap1,ap2)) {here = (here+1)%2;}
			  if (!here && TrigIsInOC(trigsaroundp.Get(l),i)) {thereOC = 1;}
			  if (!here && !TrigIsInOC(trigsaroundp.Get(l),i)) {thereNotOC = 1;}
			}
		      if (thereOC && thereNotOC)
			{
			  //(*mycout) << "Special OCICnotC - point " << pn << " found!" << endl;
			  //(*testout) << "set edgepoint due to spirals: pn=" << i << endl;
			  SetLineEndPoint(pn);
			}
		    }
		}
	    }
	}
    }
  PrintMessage(5,"have now ", GetNE(), " edges with yellow angle = ", stlparam.yangle, " degree");
}

//get trigs at a point, started with starttrig, then every left
void STLGeometry :: GetSortedTrianglesAroundPoint(int p, int starttrig, ARRAY<int>& trigs)
{
  int acttrig = starttrig;
  trigs.SetAllocSize(trigsperpoint.EntrySize(p));
  trigs.SetSize(0);
  trigs.Append(acttrig);
  int i, j, t, ap1, ap2, locindex1(0), locindex2(0);

  //(*mycout) << "trigs around point " << p << endl;

  int end = 0;
  while (!end)
    {
      const STLTriangle& at = GetTriangle(acttrig);
      for (i = 1; i <= trigsperpoint.EntrySize(p); i++)
	{
	  t = trigsperpoint.Get(p,i);
	  const STLTriangle& nt = GetTriangle(t);
	  if (at.IsNeighbourFrom(nt))
	    {
	      at.GetNeighbourPoints(nt, ap1, ap2);
	      if (ap2 == p) {Swap(ap1,ap2);}
	      if (ap1 != p) {PrintSysError("In GetSortedTrianglesAroundPoint!!!");}
	      
	      for (j = 1; j <= 3; j++) 
		{
		  if (at.PNum(j) == ap1) {locindex1 = j;};
		  if (at.PNum(j) == ap2) {locindex2 = j;};
		}
	      if ((locindex2+1)%3+1 == locindex1) 
		{
		  if (t != starttrig)
		    {
		      trigs.Append(t);
		      //		      (*mycout) << "trig " << t << endl;
		      acttrig = t;
		    }
		  else
		    {
		      end = 1;
		    }
		  break;
		}
	    }
	}
    }
  
}

/*
int STLGeometry :: NeighbourTrig(int trig, int nr) const
{
  return neighbourtrigs.Get(trig,nr);
}
*/



void STLGeometry :: SmoothGeometry ()
{
  int i, j, k;
  
  double maxerr0, maxerr;

  for (i = 1; i <= GetNP(); i++)
    {
      if (GetNEPP(i)) continue;
      
      maxerr0 = 0;
      for (j = 1; j <= NOTrigsPerPoint(i); j++)
	{
	  int tnum = TrigPerPoint(i, j);
	  double err = Angle (GetTriangle(tnum).Normal (), 
			      GetTriangle(tnum).GeomNormal(GetPoints()));
	  if (err > maxerr0)
	    maxerr0 = err;
	}

      Point3d pi = GetPoint (i);
      if (maxerr0 < 1.1) continue;    // about 60 degree

      maxerr0 /= 2;  // should be at least halfen
      
      for (k = 1; k <= NOTrigsPerPoint(i); k++)
	{
	  const STLTriangle & trig = GetTriangle (TrigPerPoint (i, k));
	  Point3d c = Center(GetPoint (trig.PNum(1)),
			     GetPoint (trig.PNum(2)),
			     GetPoint (trig.PNum(3)));

	  Point3d np = pi + 0.1 * (c - pi);
	  SetPoint (i, np);
	  
	  maxerr = 0;
	  for (j = 1; j <= NOTrigsPerPoint(i); j++)
	    {
	      int tnum = TrigPerPoint(i, j);
	      double err = Angle (GetTriangle(tnum).Normal (), 
				  GetTriangle(tnum).GeomNormal(GetPoints()));
	      if (err > maxerr)
		maxerr = err;
	    }
	  
	  if (maxerr < maxerr0)
	    {
	      pi = np;
	    }
	}

      SetPoint (i, pi);
    }
}
}
