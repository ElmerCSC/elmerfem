#include <mystdlib.h>

#include "meshing.hpp"
#include <opti.hpp>

#ifndef SMALLLIB
//#ifndef NOTCL
//#include <visual.hpp>
//#endif
#endif

namespace netgen
{

class Neighbour
{
  int nr[3];
  int orient[3];

public:
  Neighbour () { nr[0] = nr[1] = nr[2] = -1; orient[0] = orient[1] = orient[2] = 0; }

  void SetNr (int side, int anr) { nr[side-1] = anr; }
  int GetNr (int side) { return nr[side-1]; }

  void SetOrientation (int side, int aorient) { orient[side-1] = aorient; }
  int GetOrientation (int side) { return orient[side-1]; }
};




class trionedge
{
public:
  int tnr;
  int sidenr;

  trionedge () { tnr = 0; sidenr = 0; }
  trionedge (int atnr, int asidenr)
    { tnr = atnr; sidenr = asidenr; }
};



 
void MeshOptimize2d :: EdgeSwapping (Mesh & mesh, int usemetric)
{
  //  return; 

  if (!faceindex)
    {
      if (usemetric)
	PrintMessage (3, "Edgeswapping, metric");
      else
	PrintMessage (3, "Edgeswapping, topological");

      for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++)
	{
	  EdgeSwapping (mesh, usemetric);

	  if (multithread.terminate)
	    throw NgException ("Meshing stopped");
	}

      faceindex = 0;
      mesh.CalcSurfacesOfNode();
      return;
    }


  static int timer = NgProfiler::CreateTimer ("EdgeSwapping 2D");
  NgProfiler::RegionTimer reg1 (timer);


  int i, i2, j, j2;
  bool should;
  PointIndex pi;

  ARRAY<SurfaceElementIndex> seia;
  mesh.GetSurfaceElementsOfFace (faceindex, seia);

  for (i = 0; i < seia.Size(); i++)
    if (mesh[seia[i]].GetNP() != 3)
      {
	GenericImprove (mesh);
	return;
      }

  int surfnr = mesh.GetFaceDescriptor (faceindex).SurfNr();

  ARRAY<Neighbour> neighbors(mesh.GetNSE());
  INDEX_2_HASHTABLE<trionedge> other(seia.Size() + 2);


  ARRAY<char> swapped(mesh.GetNSE());
  ARRAY<int,PointIndex::BASE> pdef(mesh.GetNP());
  ARRAY<double,PointIndex::BASE> pangle(mesh.GetNP());

  SurfaceElementIndex t1, t2;
  int o1, o2;

  PointIndex pi1, pi2, pi3, pi4;
  PointGeomInfo gi1, gi2, gi3, gi4;


  int nswaps = 0;
  int e, done;
  double d;
  Vec3d nv1, nv2;
  double horder;
  double loch(-1);
  static const double minangle[] = { 0, 1.481, 2.565, 3.627, 4.683, 5.736, 7, 9 };

  pangle = 0;

  for (i = 0; i < seia.Size(); i++)
    {
      const Element2d & sel = mesh[seia[i]];
      for (j = 0; j < 3; j++)
	{
	  POINTTYPE typ = mesh[sel[j]].Type();
	  if (typ == FIXEDPOINT || typ == EDGEPOINT)
	    {
	      pangle[sel[j]] +=
		Angle (mesh[sel[(j+1)%3]] - mesh[sel[j]],
		       mesh[sel[(j+2)%3]] - mesh[sel[j]]);
	    }
	}
    }

  for (pi = PointIndex::BASE; 
       pi < mesh.GetNP()+PointIndex::BASE; pi++)
    {
      if (mesh[pi].Type() == INNERPOINT || mesh[pi].Type() == SURFACEPOINT)
	pdef[pi] = -6;
      else
	for (j = 0; j < 8; j++)
	  if (pangle[pi] >= minangle[j])
	    pdef[pi] = -1-j;
    }

  for (i = 0; i < seia.Size(); i++)
    {
      const Element2d & sel = mesh[seia[i]];
      for (j = 0; j < 3; j++)
	pdef[sel[j]]++;
    }

  for (i = 0; i < seia.Size(); i++)
    {
      //const Element2d & sel = mesh[seia[i]];
      for (j = 0; j < 3; j++)
	{
	  neighbors[seia[i]].SetNr (j+1, -1);
	  neighbors[seia[i]].SetOrientation (j+1, 0);
	}
    }

  /*
  ARRAY<Vec3d> normals(mesh.GetNP());
  for (i = 1; i <= mesh.GetNSE(); i++)
    {
      Element2d & hel = mesh.SurfaceElement(i);
      if (hel.GetIndex() == faceindex)
	for (k = 1; k <= 3; k++)
	  {
	    int pi = hel.PNum(k);
	    SelectSurfaceOfPoint (mesh.Point(pi), hel.GeomInfoPi(k));
	    int surfi = mesh.GetFaceDescriptor(faceindex).SurfNr();
	    GetNormalVector (surfi, mesh.Point(pi), normals.Elem(pi));
	    normals.Elem(pi) /= normals.Elem(pi).Length();
	  }
    }
  */	    

  
  for (i = 0; i < seia.Size(); i++)
    {
      const Element2d & sel = mesh[seia[i]];

      for (j = 1; j <= 3; j++)
	{
	  pi1 = sel.PNumMod(j+1);
	  pi2 = sel.PNumMod(j+2);
	  
	  loch = mesh.GetH(mesh[pi1]);
	    
	  INDEX_2 edge(pi1, pi2);
	  edge.Sort();
	  
	  if (mesh.IsSegment (pi1, pi2))
	    continue;
	  
	  /*
	    if (segments.Used (edge))
	    continue;
	  */
	  INDEX_2 ii2 (pi1, pi2);
	  if (other.Used (ii2))
	    {
	      // INDEX_2 i2s(ii2);
	      // i2s.Sort();
	      
	      i2 = other.Get(ii2).tnr;
	      j2 = other.Get(ii2).sidenr;
		
	      neighbors[seia[i]].SetNr (j, i2);
	      neighbors[seia[i]].SetOrientation (j, j2);
	      neighbors[i2].SetNr (j2, seia[i]);
	      neighbors[i2].SetOrientation (j2, j);
	    }
	  else
	    {
	      other.Set (INDEX_2 (pi2, pi1), trionedge (seia[i], j));
	    }
	}
    }

  for (i = 0; i < seia.Size(); i++)
    swapped[seia[i]] = 0;


  int t = 4;
  done = 0;
  while (!done && t >= 2)
    {
      for (i = 0; i < seia.Size(); i++)
	{
	  t1 = seia[i];

	  if (mesh[t1].IsDeleted())
	    continue;

	  if (mesh[t1].GetIndex() != faceindex)
	    continue;

	  if (multithread.terminate)
	    throw NgException ("Meshing stopped");

	  for (o1 = 1; o1 <= 3; o1++)
	    {
	      t2 = neighbors[t1].GetNr (o1);
	      o2 = neighbors[t1].GetOrientation (o1);

	      if (t2 == -1) continue;
	      if (swapped[t1] || swapped[t2]) continue;
	      

	      pi1 = mesh[t1].PNumMod(o1+1);
	      pi2 = mesh[t1].PNumMod(o1+2);
	      pi3 = mesh[t1].PNumMod(o1);
	      pi4 = mesh[t2].PNumMod(o2);
	      
	      gi1 = mesh[t1].GeomInfoPiMod(o1+1);
	      gi2 = mesh[t1].GeomInfoPiMod(o1+2);
	      gi3 = mesh[t1].GeomInfoPiMod(o1);
	      gi4 = mesh[t2].GeomInfoPiMod(o2);
	    
	      bool allowswap = true;

	      Vec3d auxvec1,auxvec2;

	      auxvec1 = mesh.Point(pi3)-mesh.Point(pi4);
	      auxvec2 = mesh.Point(pi1)-mesh.Point(pi4);
	      allowswap = allowswap && fabs(1.-(auxvec1*auxvec2)/(auxvec1.Length()*auxvec2.Length())) > 1e-4;

	      if(!allowswap)
		continue;

	      // normal of new
	      nv1 = Cross (auxvec1, 
			   auxvec2);

	      auxvec1 = mesh.Point(pi4)-mesh.Point(pi3);
	      auxvec2 = mesh.Point(pi2)-mesh.Point(pi3);
	      allowswap = allowswap && fabs(1.-(auxvec1*auxvec2)/(auxvec1.Length()*auxvec2.Length())) > 1e-4;


	      if(!allowswap)
		continue;

	      nv2 = Cross (auxvec1, 
			   auxvec2);

	      
	      // normals of original
	      Vec3d nv3, nv4;
	      nv3 = Cross (mesh.Point(pi1)-mesh.Point(pi4), 
			   mesh.Point(pi2)-mesh.Point(pi4));
	      nv4 = Cross (mesh.Point(pi2)-mesh.Point(pi3), 
			   mesh.Point(pi1)-mesh.Point(pi3));
	      
	      nv3 *= -1;
	      nv4 *= -1;
	      nv3.Normalize();
	      nv4.Normalize();

	      nv1.Normalize();
	      nv2.Normalize();
	    
	      Vec<3> nvp3, nvp4;
	      SelectSurfaceOfPoint (mesh.Point(pi3), gi3);
	      GetNormalVector (surfnr, mesh.Point(pi3), gi3, nvp3);

	      nvp3.Normalize();

	      SelectSurfaceOfPoint (mesh.Point(pi4), gi4);
	      GetNormalVector (surfnr, mesh.Point(pi4), gi4, nvp4);
	    
	      nvp4.Normalize();
	      
	      
	      
	      double critval = cos (M_PI / 6);  // 30 degree
	      allowswap = allowswap &&
		(nv1 * nvp3 > critval) && 
		(nv1 * nvp4 > critval) && 
		(nv2 * nvp3 > critval) && 
		(nv2 * nvp4 > critval) &&
		(nvp3 * nv3 > critval) && 
		(nvp4 * nv4 > critval);
	      

	      horder = Dist (mesh.Point(pi1), mesh.Point(pi2));

	      if ( // nv1 * nv2 >= 0 &&
		  nv1.Length() > 1e-3 * horder * horder &&
		  nv2.Length() > 1e-3 * horder * horder &&
		  allowswap )
		{
		  if (!usemetric)
		    {
		      e = pdef[pi1] + pdef[pi2] - pdef[pi3] - pdef[pi4];
		      d = 
			Dist2 (mesh.Point(pi1), mesh.Point(pi2)) - 
			Dist2 (mesh.Point(pi3), mesh.Point(pi4));
		      
		      should = e >= t && (e > 2 || d > 0);
		    }
		  else
		    {
		      should = 
			CalcTriangleBadness (mesh.Point(pi4), mesh.Point(pi3), mesh.Point(pi1), 
					     metricweight, loch) +
			CalcTriangleBadness (mesh.Point(pi3), mesh.Point(pi4), mesh.Point(pi2), 
					     metricweight, loch) <
			CalcTriangleBadness (mesh.Point(pi1), mesh.Point(pi2), mesh.Point(pi3), 
					     metricweight, loch) +
			CalcTriangleBadness (mesh.Point(pi2), mesh.Point(pi1), mesh.Point(pi4), 
					     metricweight, loch);

		    }
		  
		  if (allowswap)
		    {
		      Element2d sw1 (pi4, pi3, pi1);
		      Element2d sw2 (pi3, pi4, pi2);

		      int legal1 = 
			mesh.LegalTrig (mesh.SurfaceElement (t1)) + 
			mesh.LegalTrig (mesh.SurfaceElement (t2));
		      int legal2 = 
			mesh.LegalTrig (sw1) + mesh.LegalTrig (sw2);

		      if (legal1 < legal2) should = 1;
		      if (legal2 < legal1) should = 0;
		    }
		  
		  if (should)
		    {
		      // do swapping !
		      
		      // cout << "swap " << endl;

		      nswaps ++;
		      
		      // testout << "nv1 = " << nv1 << "   nv2 = " << nv2 << endl;
		     

		      done = 1;
		      
		      mesh[t1].PNum(1) = pi1;
		      mesh[t1].PNum(2) = pi4;
		      mesh[t1].PNum(3) = pi3;
		      
		      mesh[t2].PNum(1) = pi2;
		      mesh[t2].PNum(2) = pi3;
		      mesh[t2].PNum(3) = pi4;
		      
		      mesh[t1].GeomInfoPi(1) = gi1;
		      mesh[t1].GeomInfoPi(2) = gi4;
		      mesh[t1].GeomInfoPi(3) = gi3;
		      
		      mesh[t2].GeomInfoPi(1) = gi2;
		      mesh[t2].GeomInfoPi(2) = gi3;
		      mesh[t2].GeomInfoPi(3) = gi4;
		      
		      pdef[pi1]--;
		      pdef[pi2]--;
		      pdef[pi3]++;
		      pdef[pi4]++;
		      
		      swapped[t1] = 1;
		      swapped[t2] = 1;
		    }
		}
	    }
	}
      t--;
    }

  mesh.SetNextTimeStamp();
}







 
void MeshOptimize2d :: CombineImprove (Mesh & mesh)
{
  if (!faceindex)
    {
      PrintMessage (3, "Combine improve");

      for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++)
	{
	  CombineImprove (mesh);

	  if (multithread.terminate)
	    throw NgException ("Meshing stopped");
	}
      faceindex = 0;
      return;
    }


  static int timer = NgProfiler::CreateTimer ("Combineimprove 2D");
  NgProfiler::RegionTimer reg (timer);



  int i, j, k, l;
  PointIndex pi;
  SurfaceElementIndex sei;


  ARRAY<SurfaceElementIndex> seia;
  mesh.GetSurfaceElementsOfFace (faceindex, seia);


  for (i = 0; i < seia.Size(); i++)
    if (mesh[seia[i]].GetNP() != 3)
      return;



  int surfnr = 0;
  if (faceindex)
    surfnr = mesh.GetFaceDescriptor (faceindex).SurfNr();


  PointIndex pi1, pi2;
  MeshPoint p1, p2, pnew;
  double bad1, bad2;
  Vec<3> nv;

  int np = mesh.GetNP();
  //int nse = mesh.GetNSE();

  TABLE<SurfaceElementIndex,PointIndex::BASE> elementsonnode(np); 
  ARRAY<SurfaceElementIndex> hasonepi, hasbothpi;

  for (i = 0; i < seia.Size(); i++)
    {
      Element2d & el = mesh[seia[i]];
      for (j = 0; j < el.GetNP(); j++)
	{
	  elementsonnode.Add (el[j], seia[i]);
	}      
    }


  ARRAY<bool,PointIndex::BASE> fixed(np);
  fixed = false;

  SegmentIndex si;
  for (si = 0; si < mesh.GetNSeg(); si++)
    {
      INDEX_2 i2(mesh[si].p1, mesh[si].p2);
      fixed[i2.I1()] = true;
      fixed[i2.I2()] = true;
    }

  for(i = 0; i<mesh.LockedPoints().Size(); i++)
    fixed[mesh.LockedPoints()[i]] = true;


  ARRAY<Vec<3>,PointIndex::BASE> normals(np);

  for (pi = PointIndex::BASE; 
       pi < np + PointIndex::BASE; pi++)
    {
      if (elementsonnode[pi].Size())
	{
	  Element2d & hel = mesh[elementsonnode[pi][0]];
	  for (k = 0; k < 3; k++)
	    if (hel[k] == pi)
	      {
		SelectSurfaceOfPoint (mesh[pi], hel.GeomInfoPi(k+1));
		GetNormalVector (surfnr, mesh[pi], hel.GeomInfoPi(k+1), normals[pi]);
		break;
	      }
	  if (k == 3)
	    {
	      cerr << "Neuer Fehler von Joachim, code 17121" << endl;
	    }
	}
    }


  for (i = 0; i < seia.Size(); i++)
    {

      sei = seia[i];
      Element2d & elem = mesh[sei];
      if (elem.IsDeleted()) continue;

      for (j = 0; j < 3; j++)
	{
	  pi1 = elem[j];
	  pi2 = elem[(j+1) % 3];

	  if (pi1 < PointIndex::BASE || 
	      pi2 < PointIndex::BASE)
	    continue;

	  /*
	  INDEX_2 i2(pi1, pi2);
	  i2.Sort();
	  if (segmentht.Used(i2))
	    continue;
	  */

	  bool debugflag = 0;

	  if (debugflag)
	    {
	      (*testout) << "Combineimprove, face = " << faceindex 
			 << "pi1 = " << pi1 << " pi2 = " << pi2 << endl;
	    }

	  /*
	  // save version:
	  if (fixed.Get(pi1) || fixed.Get(pi2)) 
	    continue;
	  if (pi2 < pi1) swap (pi1, pi2);
	  */

	  // more general 
	  if (fixed[pi2]) 
	    Swap (pi1, pi2);

	  if (fixed[pi2])  
	    continue;

	  double loch = mesh.GetH (mesh[pi1]);

	  INDEX_2 si2 (pi1, pi2);
	  si2.Sort();

	  /*	  
	  if (edgetested.Used (si2))
	    continue;
	  edgetested.Set (si2, 1);
	  */

	  hasonepi.SetSize(0);
	  hasbothpi.SetSize(0);

	  for (k = 0; k < elementsonnode[pi1].Size(); k++)
	    {
	      const Element2d & el2 = mesh[elementsonnode[pi1][k]];

	      if (el2.IsDeleted()) continue;

	      if (el2[0] == pi2 || el2[1] == pi2 || el2[2] == pi2)
		{
		  hasbothpi.Append (elementsonnode[pi1][k]);
		  nv = Cross (Vec3d (mesh[el2[0]], mesh[el2[1]]),
			      Vec3d (mesh[el2[0]], mesh[el2[2]]));
		}
	      else
		{
		  hasonepi.Append (elementsonnode[pi1][k]);
		}
	    } 


	  Element2d & hel = mesh[hasbothpi[0]];
	  for (k = 0; k < 3; k++)
	    if (hel[k] == pi1)
	      {
		SelectSurfaceOfPoint (mesh[pi1],
				      hel.GeomInfoPi(k+1));
		GetNormalVector (surfnr, mesh[pi1], hel.GeomInfoPi(k+1), nv);
		break;
	      }
	  if (k == 3)
	    {
	      cerr << "Neuer Fehler von Joachim, code 32434" << endl;
	    }


	  //	  nv = normals.Get(pi1);



	  for (k = 0; k < elementsonnode[pi2].Size(); k++)
	    {
	      const Element2d & el2 = mesh[elementsonnode[pi2][k]];
	      if (el2.IsDeleted()) continue;

	      if (el2[0] == pi1 || el2[1] == pi1 || el2[2] == pi1)
		;
	      else
		hasonepi.Append (elementsonnode[pi2][k]);
	    } 

	  bad1 = 0;
	  int illegal1 = 0, illegal2 = 0;
	  for (k = 0; k < hasonepi.Size(); k++)
	    {
	      const Element2d & el = mesh[hasonepi[k]];
	      bad1 += CalcTriangleBadness (mesh[el[0]], mesh[el[1]], mesh[el[2]],
					   nv, -1, loch);
	      illegal1 += 1-mesh.LegalTrig(el);
	    }
	  
	  for (k = 0; k < hasbothpi.Size(); k++)
	    {
	      const Element2d & el = mesh[hasbothpi[k]];
	      bad1 += CalcTriangleBadness (mesh[el[0]], mesh[el[1]], mesh[el[2]],
					   nv, -1, loch);
	      illegal1 += 1-mesh.LegalTrig(el);
	    }
	  bad1 /= (hasonepi.Size()+hasbothpi.Size());

	  p1 = mesh[pi1];
	  p2 = mesh[pi2];

	  pnew = p1;
	  mesh[pi1] = pnew;
	  mesh[pi2] = pnew;

	  bad2 = 0;
	  for (k = 0; k < hasonepi.Size(); k++)
	    {
	      Element2d & el = mesh[hasonepi[k]];
	      double err = 
		CalcTriangleBadness (mesh[el[0]], mesh[el[1]], mesh[el[2]],
				     nv, -1, loch);
	      bad2 += err;

	      Vec<3> hnv = Cross (Vec3d (mesh[el[0]],
                                         mesh[el[1]]),
                                  Vec3d (mesh[el[0]],
                                         mesh[el[2]]));
	      if (hnv * nv < 0)
		bad2 += 1e10;
              
	      for (l = 0; l < 3; l++)
		if ( (normals[el[l]] * nv) < 0.5)
		  bad2 += 1e10;

	      illegal2 += 1-mesh.LegalTrig(el);
	    }
	  bad2 /= hasonepi.Size();

	  mesh[pi1] = p1;
	  mesh[pi2] = p2;
	  
       
	  if (debugflag)
	    {
	      (*testout) << "bad1 = " << bad1 << ", bad2 = " << bad2 << endl;
	    }


	  bool should = (bad2 < bad1 && bad2 < 1e4);
	  if (bad2 < 1e4)
	    {
	      if (illegal1 > illegal2) should = 1;
	      if (illegal2 > illegal1) should = 0;
	    }
	  

	  if (should)
	    {
	      // (*testout) << "combine !" << endl;
	      // (*testout) << "bad1 = " << bad1 << ", bad2 = " << bad2 << endl;


	      mesh[pi1] = pnew;
	      PointGeomInfo gi;
	      bool gi_set(false);
	      
	      
	      Element2d *el1p(NULL);
	      l=0;
	      while(mesh[elementsonnode[pi1][l]].IsDeleted() && l<elementsonnode.EntrySize(pi1)) l++;
	      if(l<elementsonnode.EntrySize(pi1))
		el1p = &mesh[elementsonnode[pi1][l]];
	      else
		cerr << "OOPS!" << endl;

	      for (l = 0; l < el1p->GetNP(); l++)
		if ((*el1p)[l] == pi1)
		  {
		    gi = el1p->GeomInfoPi (l+1);
		    gi_set = true;
		  }

	      // (*testout) << "Connect point " << pi2 << " to " << pi1 << "\n";
	      for (k = 0; k < elementsonnode[pi2].Size(); k++)
		{
		  Element2d & el = mesh[elementsonnode[pi2][k]];
		  if(el.IsDeleted()) continue;
		  elementsonnode.Add (pi1, elementsonnode[pi2][k]);

		  bool haspi1 = 0;
		  for (l = 0; l < el.GetNP(); l++)
		    if (el[l] == pi1)
		      haspi1 = 1;
		  if (haspi1) continue;

		  for (l = 0; l < el.GetNP(); l++)
		    {
		      if (el[l] == pi2)
			{
			  el[l] = pi1;
			  el.GeomInfoPi (l+1) = gi;
			}

		      fixed[el[l]] = true;
		    }
		}

	      /*
	      for (k = 0; k < hasbothpi.Size(); k++)
		{
		  cout << mesh[hasbothpi[k]] << endl;
		  for (l = 0; l < 3; l++)
		    cout << mesh[mesh[hasbothpi[k]][l]] << " ";
		  cout << endl;
		}
	      */

	      for (k = 0; k < hasbothpi.Size(); k++)
		{
		  mesh[hasbothpi[k]].Delete();
		  /*
		  for (l = 0; l < 4; l++)
		    mesh[hasbothpi[k]][l] = PointIndex::BASE-1;
		  */
		}

	    }
	}
    }

  //  mesh.Compress();
  mesh.SetNextTimeStamp();
}


void MeshOptimize2d :: CheckMeshApproximation (Mesh & mesh)
{
  // Check angles between elements and normals at corners
  /*
  
  int i, j;
  int ne = mesh.GetNSE();
  int surfnr;
  
  Vec3d n, ng;
  ARRAY<Vec3d> ngs(3);

  (*mycout) << "Check Surface Approxiamtion" << endl;
  (*testout) << "Check Surface Approxiamtion" << endl;

  for (i = 1; i <= ne; i++)
    {
      const Element2d & el = mesh.SurfaceElement(i);
      surfnr = mesh.GetFaceDescriptor (el.GetIndex()).SurfNr();
      Vec3d n = Cross (mesh.Point (el.PNum(1)) - mesh.Point (el.PNum(2)),
		       mesh.Point (el.PNum(1)) - mesh.Point (el.PNum(3)));
      n /= n.Length();

      for (j = 1; j <= el.GetNP(); j++)
	{
	  SelectSurfaceOfPoint (mesh.Point(el.PNum(j)), el.GeomInfoPi(j));
	  GetNormalVector (surfnr, mesh.Point(el.PNum(j)), ng);
	  ng /= ng.Length();
	  ngs.Elem(j) = ng;

	  double angle =  (180.0 / M_PI) * Angle (n, ng);
	  if (angle > 60)
	    {
	      (*testout) << "el " << i << " node " << el.PNum(j)
			 << "has angle = " << angle << endl;
	    }
	}	

      for (j = 1; j <= 3; j++)
	{
	  double angle =  (180.0 / M_PI) * Angle (ngs.Get(j), ngs.Get(j%3+1));
	  if (angle > 60)
	    {
	      (*testout) << "el " << i << " node-node " 
			 << ngs.Get(j) << " - " << ngs.Get(j%3+1)
			 << " has angle = " << angle << endl;
	    }
	}
    }
  */
}
}
