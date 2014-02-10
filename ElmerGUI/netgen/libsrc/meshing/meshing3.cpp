#include <mystdlib.h>
#include "meshing.hpp"

namespace netgen
{

double minother;
double minwithoutother;





MeshingStat3d :: MeshingStat3d ()
{
  cntsucc = cnttrials = cntelem = qualclass = 0;
  vol0 = h = 1;
  problemindex = 1;
}  
  

Meshing3 :: Meshing3 (const string & rulefilename) 
{
  tolfak = 1;

  LoadRules (rulefilename.c_str(), NULL);
  adfront = new AdFront3;

  problems.SetSize (rules.Size());
  foundmap.SetSize (rules.Size());
  canuse.SetSize (rules.Size());
  ruleused.SetSize (rules.Size());

  for (int i = 1; i <= rules.Size(); i++)
    {
      problems.Elem(i) = new char[255];
      foundmap.Elem(i) = 0;
      canuse.Elem(i) = 0;
      ruleused.Elem(i) = 0;
    }
}


Meshing3 :: Meshing3 (const char ** rulep)
{
  tolfak = 1;

  LoadRules (NULL, rulep);
  adfront = new AdFront3;

  problems.SetSize (rules.Size());
  foundmap.SetSize (rules.Size());
  canuse.SetSize (rules.Size());
  ruleused.SetSize (rules.Size());

  for (int i = 0; i < rules.Size(); i++)
    {
      problems[i] = new char[255];
      foundmap[i] = 0;
      canuse[i]   = 0;
      ruleused[i] = 0;
    }
}

Meshing3 :: ~Meshing3 ()
{
  delete adfront;
  for (int i = 0; i < rules.Size(); i++)
    {
      delete [] problems[i];
      delete rules[i];
    }
}



static double CalcLocH (const ARRAY<Point3d> & locpoints,
			const ARRAY<MiniElement2d> & locfaces,
			double h)
{
  return h;

  // was war das ????
  
  int i, j;
  double hi, h1, d, dn, sum, weight, wi;
  Point3d p0, pc;
  Vec3d n, v1, v2;

  p0.X() = p0.Y() = p0.Z() = 0;
  for (j = 1; j <= 3; j++)
    {
      p0.X() += locpoints.Get(locfaces.Get(1).PNum(j)).X();
      p0.Y() += locpoints.Get(locfaces.Get(1).PNum(j)).Y();
      p0.Z() += locpoints.Get(locfaces.Get(1).PNum(j)).Z();
    }
  p0.X() /= 3; p0.Y() /= 3; p0.Z() /= 3;
  
  v1 = locpoints.Get(locfaces.Get(1).PNum(2)) -
    locpoints.Get(locfaces.Get(1).PNum(1));
  v2 = locpoints.Get(locfaces.Get(1).PNum(3)) -
    locpoints.Get(locfaces.Get(1).PNum(1));

  h1 = v1.Length();
  n = Cross (v2, v1);
  n /= n.Length();

  sum = 0;
  weight = 0;

  for (i = 1; i <= locfaces.Size(); i++)
    {
      pc.X() = pc.Y() = pc.Z() = 0;
      for (j = 1; j <= 3; j++)
	{
	  pc.X() += locpoints.Get(locfaces.Get(i).PNum(j)).X();
	  pc.Y() += locpoints.Get(locfaces.Get(i).PNum(j)).Y();
	  pc.Z() += locpoints.Get(locfaces.Get(i).PNum(j)).Z();
	}
      pc.X() /= 3; pc.Y() /= 3; pc.Z() /= 3;

      d = Dist (p0, pc);
      dn = n * (pc - p0);
      hi = Dist (locpoints.Get(locfaces.Get(i).PNum(1)),
		 locpoints.Get(locfaces.Get(i).PNum(2)));
		 
      if (dn > -0.2 * h1)
	{
	  wi = 1 / (h1 + d);
	  wi *= wi;
	}
      else
	wi = 0;

      sum += hi * wi;
      weight += wi;
    }

  return sum/weight;
}


PointIndex Meshing3 :: AddPoint (const Point3d & p, PointIndex globind)
{
  return adfront -> AddPoint (p, globind);  
}  

void Meshing3 :: AddBoundaryElement (const Element2d & elem)
{
  MiniElement2d mini(elem.GetNP());
  for (int j = 0; j < elem.GetNP(); j++)
    mini[j] = elem[j];
  adfront -> AddFace(mini);
}  


void Meshing3 :: AddBoundaryElement (const MiniElement2d & elem)
{
  adfront -> AddFace(elem);
}

int Meshing3 :: AddConnectedPair (const INDEX_2 & apair)
{
  return adfront -> AddConnectedPair (apair);
}

MESHING3_RESULT Meshing3 :: 
GenerateMesh (Mesh & mesh, const MeshingParameters & mp)
{
  static int meshing3_timer = NgProfiler::CreateTimer ("Meshing3::GenerateMesh");
  static int meshing3_timer_a = NgProfiler::CreateTimer ("Meshing3::GenerateMesh a");
  static int meshing3_timer_b = NgProfiler::CreateTimer ("Meshing3::GenerateMesh b");
  static int meshing3_timer_c = NgProfiler::CreateTimer ("Meshing3::GenerateMesh c");
  static int meshing3_timer_d = NgProfiler::CreateTimer ("Meshing3::GenerateMesh d");
  NgProfiler::RegionTimer reg (meshing3_timer);


  ARRAY<Point3d > locpoints;      // local points
  ARRAY<MiniElement2d> locfaces;     // local faces
  ARRAY<PointIndex> pindex;      // mapping from local to front point numbering
  ARRAY<int> allowpoint;         // point is allowd ?
  ARRAY<INDEX> findex;           // mapping from local to front face numbering
  //INDEX_2_HASHTABLE<int> connectedpairs(100);  // connecgted pairs for prism meshing

  ARRAY<Point3d > plainpoints;       // points in reference coordinates
  ARRAY<int> delpoints, delfaces;   // points and lines to be deleted
  ARRAY<Element> locelements;       // new generated elements

  int i, j, oldnp, oldnf;
  int found;
  referencetransform trans;
  int rotind;
  INDEX globind;
  Point3d inp;
  float err;

  INDEX locfacesplit;             //index for faces in outer area
  
  bool loktestmode = false;

  int uselocalh = mp.uselocalh;

  // int giveuptol = mp.giveuptol; // 
  MeshingStat3d stat;      // statistics
  int plotstat_oldne = -1;

  
  // for star-shaped domain meshing
  ARRAY<MeshPoint> grouppoints;      
  ARRAY<MiniElement2d> groupfaces;
  ARRAY<PointIndex> grouppindex;
  ARRAY<INDEX> groupfindex;
  
  
  float minerr;
  int hasfound;
  double tetvol;
  // int giveup = 0;

  
  ARRAY<Point3d> tempnewpoints;
  ARRAY<MiniElement2d> tempnewfaces;
  ARRAY<int> tempdelfaces;
  ARRAY<Element> templocelements;


  stat.h = mp.maxh;

  adfront->SetStartFront (mp.baseelnp);


  found = 0;
  stat.vol0 = adfront -> Volume();
  tetvol = 0;

  stat.qualclass = 1;

  while (1)
    {
      if (multithread.terminate)
	throw NgException ("Meshing stopped");

      // break if advancing front is empty
      if (!mp.baseelnp && adfront->Empty())
	break;

      // break, if advancing front has no elements with
      // mp.baseelnp nodes  
      if (mp.baseelnp && adfront->Empty (mp.baseelnp))
	break;

      locpoints.SetSize(0);
      locfaces.SetSize(0);
      locelements.SetSize(0);
      pindex.SetSize(0);
      findex.SetSize(0);

      INDEX_2_HASHTABLE<int> connectedpairs(100);  // connected pairs for prism meshing
      
      // select base-element (will be locface[1])
      // and get local environment of radius (safety * h)


      int baseelem = adfront -> SelectBaseElement ();
      if (mp.baseelnp && adfront->GetFace (baseelem).GetNP() != mp.baseelnp)
	{
	  adfront->IncrementClass (baseelem);	  
	  continue;
	}

      const MiniElement2d & bel = adfront->GetFace (baseelem);
      const Point3d & p1 = adfront->GetPoint (bel.PNum(1));
      const Point3d & p2 = adfront->GetPoint (bel.PNum(2));
      const Point3d & p3 = adfront->GetPoint (bel.PNum(3));

      // (*testout) << endl << "base = " << bel << endl;


      Point3d pmid = Center (p1, p2, p3);

      double his = (Dist (p1, p2) + Dist(p1, p3) + Dist(p2, p3)) / 3;
      double hshould;

      hshould = mesh.GetH (pmid);

      if (adfront->GetFace (baseelem).GetNP() == 4)
	hshould = max2 (his, hshould);

      double hmax = (his > hshould) ? his : hshould;
      
      // qualclass should come from baseelem !!!!!
      double hinner = hmax * (1 + stat.qualclass);
      double houter = hmax * (1 + 2 * stat.qualclass);

      NgProfiler::StartTimer (meshing3_timer_a);
      stat.qualclass =
        adfront -> GetLocals (baseelem, locpoints, locfaces, 
			      pindex, findex, connectedpairs,
			      houter, hinner,
			      locfacesplit);
      NgProfiler::StopTimer (meshing3_timer_a);

      // (*testout) << "locfaces = " << endl << locfaces << endl;

      int pi1 = pindex.Get(locfaces[0].PNum(1));
      int pi2 = pindex.Get(locfaces[0].PNum(2));
      int pi3 = pindex.Get(locfaces[0].PNum(3));

      //loktestmode = 1;
      testmode = loktestmode;  //changed 
      // loktestmode = testmode =  (adfront->GetFace (baseelem).GetNP() == 4) && (rules.Size() == 5);

      loktestmode = stat.qualclass > 5;
      

      if (loktestmode)
	{
	  (*testout) << "baseel = " << baseelem << ", ind = " << findex.Get(1) << endl;
	  (*testout) << "pi = " << pi1 << ", " << pi2 << ", " << pi3 << endl;
	}





      if (testmode)
	{
	  (*testout) << "baseelem = " << baseelem << " qualclass = " << stat.qualclass << endl;
	  (*testout) << "locpoints = " << endl << locpoints << endl;
	  (*testout) << "connected = " << endl << connectedpairs << endl;
	}



      // loch = CalcLocH (locpoints, locfaces, h);
      
      stat.nff = adfront->GetNF();
      stat.vol = adfront->Volume();
      if (stat.vol < 0) break;

      oldnp = locpoints.Size();
      oldnf = locfaces.Size();


      allowpoint.SetSize(locpoints.Size());
      if (uselocalh && stat.qualclass <= 3)
	for (i = 1; i <= allowpoint.Size(); i++)
	  {
	    allowpoint.Elem(i) =
	      (mesh.GetH (locpoints.Get(i)) > 0.4 * hshould / mp.sloppy) ? 2 : 1;
	  }
      else
	allowpoint = 2;


      
      if (stat.qualclass >= mp.starshapeclass &&
	  mp.baseelnp != 4)   
	{
	  NgProfiler::RegionTimer reg1 (meshing3_timer_b);
	  // star-shaped domain removing

	  grouppoints.SetSize (0);
	  groupfaces.SetSize (0);
	  grouppindex.SetSize (0);
	  groupfindex.SetSize (0);
	  
	  adfront -> GetGroup (findex[0], grouppoints, groupfaces, 
			       grouppindex, groupfindex);

	  bool onlytri = 1;
	  for (i = 0; i < groupfaces.Size(); i++)
	    if (groupfaces[i].GetNP() != 3) 
	      onlytri = 0;
	  
	  if (onlytri && groupfaces.Size() <= 20 + 2*stat.qualclass &&
	      FindInnerPoint (grouppoints, groupfaces, inp))
	    {
	      (*testout) << "inner point found" << endl;

	      for (i = 1; i <= groupfaces.Size(); i++)
		adfront -> DeleteFace (groupfindex.Get(i));
	      
	      for (i = 1; i <= groupfaces.Size(); i++)
		for (j = 1; j <= locfaces.Size(); j++)
		  if (findex.Get(j) == groupfindex.Get(i))
		    delfaces.Append (j);
	      
	      
	      delfaces.SetSize (0);
	      
	      INDEX npi;
	      Element newel;
	      
	      npi = mesh.AddPoint (inp);
	      newel.SetNP(4);
	      newel.PNum(4) = npi;
	      
	      for (i = 1; i <= groupfaces.Size(); i++)
		{
		  for (j = 1; j <= 3; j++)
		    {
		      newel.PNum(j) = 
			adfront->GetGlobalIndex 
			(grouppindex.Get(groupfaces.Get(i).PNum(j)));
		    }
		  mesh.AddVolumeElement (newel);
		}
	      continue;
	    }
	}
      
      found = 0;
      hasfound = 0;
      minerr = 1e6;

      //      int optother = 0;

      /*
      for (i = 1; i <= locfaces.Size(); i++)
	{
	  (*testout) << "Face " << i << ": ";
	  for (j = 1; j <= locfaces.Get(i).GetNP(); j++)
	    (*testout) << pindex.Get(locfaces.Get(i).PNum(j)) << " ";
	  (*testout) << endl;
	}
      for (i = 1; i <= locpoints.Size(); i++)
	{
	  (*testout) << "p" << i 
		     << ", gi = " << pindex.Get(i) 
		     << " = " << locpoints.Get(i) << endl;
	}
	*/

      minother = 1e10;
      minwithoutother = 1e10;

      bool impossible = 1;

      for (rotind = 1; rotind <= locfaces[0].GetNP(); rotind++)
	{
	  // set transformatino to reference coordinates

	  if (locfaces.Get(1).GetNP() == 3)
	    {
	      trans.Set (locpoints.Get(locfaces.Get(1).PNumMod(1+rotind)),
			 locpoints.Get(locfaces.Get(1).PNumMod(2+rotind)),
			 locpoints.Get(locfaces.Get(1).PNumMod(3+rotind)), hshould);
	    }
	  else
	    {
	      trans.Set (locpoints.Get(locfaces.Get(1).PNumMod(1+rotind)),
			 locpoints.Get(locfaces.Get(1).PNumMod(2+rotind)),
			 locpoints.Get(locfaces.Get(1).PNumMod(4+rotind)), hshould);
	    }

	  trans.ToPlain (locpoints, plainpoints);


	  for (i = 1; i <= allowpoint.Size(); i++)
	    {
	      if (plainpoints.Get(i).Z() > 0)
		{
		  //if(loktestmode)
		  //  (*testout) << "plainpoints.Get(i).Z() = " << plainpoints.Get(i).Z() << " > 0" << endl;
		  allowpoint.Elem(i) = 0;
		}
	    }

	  stat.cnttrials++;


	  if (stat.cnttrials % 100 == 0)
	    {
	      (*testout) << "\n";
	      for (i = 1; i <= canuse.Size(); i++)
	      {
		(*testout) << foundmap.Get(i) << "/" 
			   << canuse.Get(i) << "/"
			   << ruleused.Get(i) << " map/can/use rule " << rules.Get(i)->Name() << "\n";
	      }
	      (*testout) << endl;
	    }

	  NgProfiler::StartTimer (meshing3_timer_c);	  

	  found = ApplyRules (plainpoints, allowpoint, 
			      locfaces, locfacesplit, connectedpairs,
			      locelements, delfaces, 
			      stat.qualclass, mp.sloppy, rotind, err);

	  if (found >= 0) impossible = 0;
	  if (found < 0) found = 0;


	  NgProfiler::StopTimer (meshing3_timer_c);	  

	  if (!found) loktestmode = 0;

	  NgProfiler::RegionTimer reg2 (meshing3_timer_d);	  
	  
	  if (loktestmode)
	    {
	      (*testout) << "plainpoints = " << endl << plainpoints << endl;
	      (*testout) << "Applyrules found " << found << endl;
	    }

	  if (found) stat.cntsucc++;

	  locpoints.SetSize (plainpoints.Size());
	  for (i = oldnp+1; i <= plainpoints.Size(); i++)
	    trans.FromPlain (plainpoints.Elem(i), locpoints.Elem(i));
	  


	  // avoid meshing from large to small mesh-size
	  if (uselocalh && found && stat.qualclass <= 3)
	    {
	      for (i = 1; i <= locelements.Size(); i++)
		{
		  Point3d pmin = locpoints.Get(locelements.Get(i).PNum(1));
		  Point3d pmax = pmin;
		  for (j = 2; j <= 4; j++)
		    {
		      const Point3d & hp = locpoints.Get(locelements.Get(i).PNum(j));
		      pmin.SetToMin (hp);
		      pmax.SetToMax (hp);
		    }

		  if (mesh.GetMinH (pmin, pmax) < 0.4 * hshould / mp.sloppy)
		    found = 0;
		}
	    }
	  if (found)
	    {
	      for (i = 1; i <= locelements.Size(); i++)
		for (j = 1; j <= 4; j++)
		  {
		    const Point3d & hp = locpoints.Get(locelements.Get(i).PNum(j));
		    if (Dist (hp, pmid) > hinner)
		      found = 0;
		  }
	    }


	  if (found)
	    ruleused.Elem(found)++;


	  // plotstat->Plot(stat);
	  
	  if (stat.cntelem != plotstat_oldne)
	    {
	      plotstat_oldne = stat.cntelem;

	      PrintMessageCR (5, "El: ", stat.cntelem,
			      //	    << " trials: " << stat.cnttrials
			      " faces: ", stat.nff,
			      " vol = ", float(100 * stat.vol / stat.vol0));
  
	      multithread.percent = 100 -  100.0 * stat.vol / stat.vol0;
	    }


	  if (found && (!hasfound || err < minerr) )
	    {
	      
	      if (testmode)
		{
		  (*testout) << "found is active, 3" << endl;
		  for (i = 1; i <= plainpoints.Size(); i++)
		    {
		      (*testout) << "p";
		      if (i <= pindex.Size())
			(*testout) << pindex.Get(i) << ": ";
		      else
			(*testout) << "new: ";
		      (*testout) << plainpoints.Get(i) << endl;
		    }
		}
	      
	      
	      
	      hasfound = found;
	      minerr = err;
	      
	      tempnewpoints.SetSize (0);
	      for (i = oldnp+1; i <= locpoints.Size(); i++)
		tempnewpoints.Append (locpoints.Get(i));
	      
	      tempnewfaces.SetSize (0);
	      for (i = oldnf+1; i <= locfaces.Size(); i++)
		tempnewfaces.Append (locfaces.Get(i));
	      
	      tempdelfaces.SetSize (0);
	      for (i = 1; i <= delfaces.Size(); i++)
		tempdelfaces.Append (delfaces.Get(i));
	      
	      templocelements.SetSize (0);
	      for (i = 1; i <= locelements.Size(); i++)
		templocelements.Append (locelements.Get(i));

	      /*
	      optother =
		strcmp (problems[found], "other") == 0;
	      */
	    }
	  
	  locpoints.SetSize (oldnp);
	  locfaces.SetSize (oldnf);
	  delfaces.SetSize (0);
	  locelements.SetSize (0);
	}
      
      

      if (hasfound)
	{

	  /*
	  if (optother)
	    (*testout) << "Other is optimal" << endl;

	  if (minother < minwithoutother)
	    {
	      (*testout) << "Other is better, " << minother << " less " << minwithoutother << endl;
	    }
	    */

	  for (i = 1; i <= tempnewpoints.Size(); i++)
	    locpoints.Append (tempnewpoints.Get(i));
	  for (i = 1; i <= tempnewfaces.Size(); i++)
	    locfaces.Append (tempnewfaces.Get(i));
	  for (i = 1; i <= tempdelfaces.Size(); i++)
	    delfaces.Append (tempdelfaces.Get(i));
	  for (i = 1; i <= templocelements.Size(); i++)
	    locelements.Append (templocelements.Get(i));


	  if (loktestmode)
	    {
	      (*testout) << "apply rule" << endl;
	      for (i = 1; i <= locpoints.Size(); i++)
		{
		  (*testout) << "p";
		  if (i <= pindex.Size())
		    (*testout) << pindex.Get(i) << ": ";
		  else
		    (*testout) << "new: ";
		  (*testout) << locpoints.Get(i) << endl;
		}
	    }



	  pindex.SetSize(locpoints.Size());

	  for (i = oldnp+1; i <= locpoints.Size(); i++)
	    {
	      globind = mesh.AddPoint (locpoints.Get(i));
	      pindex.Elem(i) = adfront -> AddPoint (locpoints.Get(i), globind);
	    }

	  for (i = 1; i <= locelements.Size(); i++)
	    {
	      Point3d * hp1, * hp2, * hp3, * hp4;
	      hp1 = &locpoints.Elem(locelements.Get(i).PNum(1));
	      hp2 = &locpoints.Elem(locelements.Get(i).PNum(2));
	      hp3 = &locpoints.Elem(locelements.Get(i).PNum(3));
	      hp4 = &locpoints.Elem(locelements.Get(i).PNum(4));
	      
	      tetvol += (1.0 / 6.0) * ( Cross ( *hp2 - *hp1, *hp3 - *hp1) * (*hp4 - *hp1) );

	      for (j = 1; j <= locelements.Get(i).NP(); j++)
		locelements.Elem(i).PNum(j) =
		  adfront -> GetGlobalIndex (pindex.Get(locelements.Get(i).PNum(j)));

	      mesh.AddVolumeElement (locelements.Get(i));
	      stat.cntelem++;
	    }

	  for (i = oldnf+1; i <= locfaces.Size(); i++)
	    {
	      for (j = 1; j <= locfaces.Get(i).GetNP(); j++)
		locfaces.Elem(i).PNum(j) = 
		  pindex.Get(locfaces.Get(i).PNum(j));
	      // (*testout) << "add face " << locfaces.Get(i) << endl;
	      adfront->AddFace (locfaces.Get(i));
	    }
	  
	  for (i = 1; i <= delfaces.Size(); i++)
	    adfront->DeleteFace (findex.Get(delfaces.Get(i)));
	}
      else
	{
	  adfront->IncrementClass (findex.Get(1));
	  if (impossible && mp.check_impossible)
	    {
	      (*testout) << "skip face since it is impossible" << endl;
	      for (j = 0; j < 100; j++)
		adfront->IncrementClass (findex.Get(1));
	    }
	}

      locelements.SetSize (0);
      delpoints.SetSize(0);
      delfaces.SetSize(0);

      if (stat.qualclass >= mp.giveuptol)
	break;
    }
  
  PrintMessage (5, "");  // line feed after statistics

  for (i = 1; i <= ruleused.Size(); i++)
    (*testout) << setw(4) << ruleused.Get(i)
	       << " times used rule " << rules.Get(i) -> Name() << endl;


  if (!mp.baseelnp && adfront->Empty())
    return MESHING3_OK;

  if (mp.baseelnp && adfront->Empty (mp.baseelnp))
    return MESHING3_OK;

  if (stat.vol < -1e-15)
    return MESHING3_NEGVOL;

  return MESHING3_NEGVOL;
}




enum blocktyp { BLOCKUNDEF, BLOCKINNER, BLOCKBOUND, BLOCKOUTER };

void Meshing3 :: BlockFill (Mesh & mesh, double gh)
{
  PrintMessage (3, "Block-filling called (obsolete) ");

  int i, j(0), i1, i2, i3, j1, j2, j3;
  int n1, n2, n3, n, min1, min2, min3, max1, max2, max3;
  int changed, filled;
  double xmin(0), xmax(0), ymin(0), ymax(0), zmin(0), zmax(0);
  double xminb, xmaxb, yminb, ymaxb, zminb, zmaxb;
  //double rad = 0.7 * gh;
  
  for (i = 1; i <= adfront->GetNP(); i++)
    {
      const Point3d & p = adfront->GetPoint(i);
      if (i == 1)
	{
	  xmin = xmax = p.X();
	  ymin = ymax = p.Y();
	  zmin = zmax = p.Z();
	}
      else
	{
	  if (p.X() < xmin) xmin = p.X();
	  if (p.X() > xmax) xmax = p.X();
	  if (p.Y() < ymin) ymin = p.Y();
	  if (p.Y() > ymax) ymax = p.Y();
	  if (p.Z() < zmin) zmin = p.Z();
	  if (p.Z() > zmax) zmax = p.Z();
	}
    }
  
  xmin -= 5 * gh;
  ymin -= 5 * gh;
  zmin -= 5 * gh;
  
  n1 = int ((xmax-xmin) / gh + 5);
  n2 = int ((ymax-ymin) / gh + 5);
  n3 = int ((zmax-zmin) / gh + 5);
  n = n1 * n2 * n3;
  
  PrintMessage (5, "n1 = ", n1, " n2 = ", n2, " n3 = ", n3);

  ARRAY<blocktyp> inner(n);
  ARRAY<int> pointnr(n), frontpointnr(n);


  // initialize inner to 1

  for (i = 1; i <= n; i++)
    inner.Elem(i) = BLOCKUNDEF;


  // set blocks cutting surfaces to 0

  for (i = 1; i <= adfront->GetNF(); i++)
    {
      const MiniElement2d & el = adfront->GetFace(i);
      xminb = xmax; xmaxb = xmin;
      yminb = ymax; ymaxb = ymin;
      zminb = zmax; zmaxb = zmin;

      for (j = 1; j <= 3; j++)
	{
	  const Point3d & p = adfront->GetPoint (el.PNum(j));
	  if (p.X() < xminb) xminb = p.X();
	  if (p.X() > xmaxb) xmaxb = p.X();
	  if (p.Y() < yminb) yminb = p.Y();
	  if (p.Y() > ymaxb) ymaxb = p.Y();
	  if (p.Z() < zminb) zminb = p.Z();
	  if (p.Z() > zmaxb) zmaxb = p.Z();
	}

	

      double filldist = 0.2; // globflags.GetNumFlag ("filldist", 0.4);
      xminb -= filldist * gh;
      xmaxb += filldist * gh;
      yminb -= filldist * gh;
      ymaxb += filldist * gh;
      zminb -= filldist * gh;
      zmaxb += filldist * gh;

      min1 = int ((xminb - xmin) / gh) + 1;
      max1 = int ((xmaxb - xmin) / gh) + 1;
      min2 = int ((yminb - ymin) / gh) + 1;
      max2 = int ((ymaxb - ymin) / gh) + 1;
      min3 = int ((zminb - zmin) / gh) + 1;
      max3 = int ((zmaxb - zmin) / gh) + 1;


      for (i1 = min1; i1 <= max1; i1++)
	for (i2 = min2; i2 <= max2; i2++)
	  for (i3 = min3; i3 <= max3; i3++)
	    inner.Elem(i3 + (i2-1) * n3 + (i1-1) * n2 * n3) = BLOCKBOUND;      
    }

  


  while (1)
    {
      int undefi = 0;
      Point3d undefp;

      for (i1 = 1; i1 <= n1 && !undefi; i1++)
	for (i2 = 1; i2 <= n2 && !undefi; i2++)
	  for (i3 = 1; i3 <= n3 && !undefi; i3++)
	    {
	      i = i3 + (i2-1) * n3 + (i1-1) * n2 * n3;
	      if (inner.Elem(i) == BLOCKUNDEF)
		{
		  undefi = i;
		  undefp.X() = xmin + (i1-0.5) * gh;
		  undefp.Y() = ymin + (i2-0.5) * gh;
		  undefp.Z() = zmin + (i3-0.5) * gh;
		}
	    }
	      
      if (!undefi)
	break;

      //      PrintMessage (5, "Test point: ", undefp);
      
      if (adfront -> Inside (undefp))
	{
	  //	  (*mycout) << "inner" << endl;
	  inner.Elem(undefi) = BLOCKINNER;
	}
      else
	{
	  //	  (*mycout) << "outer" << endl;
	  inner.Elem(undefi) = BLOCKOUTER;
	}

      do
	{
	  changed = 0;
	  for (i1 = 1; i1 <= n1; i1++)
	    for (i2 = 1; i2 <= n2; i2++)
	      for (i3 = 1; i3 <= n3; i3++)
		{
		  i = i3 + (i2-1) * n3 + (i1-1) * n2 * n3;

		  for (int k = 1; k <= 3; k++)
		    {
		      switch (k)
			{
			case 1: j = i + n2 * n3; break;
			case 2: j = i + n3; break;
			case 3: j = i + 1; break;
			}
		  
		      if (j > n1 * n2 * n3) continue;

		      if (inner.Elem(i) == BLOCKOUTER && inner.Elem(j) == BLOCKUNDEF)
			{
			  changed = 1;
			  inner.Elem(j) = BLOCKOUTER;
			}
		      if (inner.Elem(j) == BLOCKOUTER && inner.Elem(i) == BLOCKUNDEF)
			{
			  changed = 1;
			  inner.Elem(i) = BLOCKOUTER;
			}
		      if (inner.Elem(i) == BLOCKINNER && inner.Elem(j) == BLOCKUNDEF)
			{
			  changed = 1;
			  inner.Elem(j) = BLOCKINNER;
			}
		      if (inner.Elem(j) == BLOCKINNER && inner.Elem(i) == BLOCKUNDEF)
			{
			  changed = 1;
			  inner.Elem(i) = BLOCKINNER;
			}
		    }
		}
	}
      while (changed); 

    }



  filled = 0;
  for (i = 1; i <= n; i++)
    if (inner.Elem(i) == BLOCKINNER)
      {
	filled++;
      }
  PrintMessage (5, "Filled blocks: ", filled);

  for (i = 1; i <= n; i++)
    {
      pointnr.Elem(i) = 0;
      frontpointnr.Elem(i) = 0;
    }
  
  for (i1 = 1; i1 <= n1-1; i1++)
    for (i2 = 1; i2 <= n2-1; i2++)
      for (i3 = 1; i3 <= n3-1; i3++)
	{
	  i = i3 + (i2-1) * n3 + (i1-1) * n2 * n3;
	  if (inner.Elem(i) == BLOCKINNER)
	    {
	      for (j1 = i1; j1 <= i1+1; j1++)
		for (j2 = i2; j2 <= i2+1; j2++)
		  for (j3 = i3; j3 <= i3+1; j3++)
		    {
		      j = j3 + (j2-1) * n3 + (j1-1) * n2 * n3;
		      if (pointnr.Get(j) == 0)
			{
			  Point3d hp(xmin + (j1-1) * gh, 
				     ymin + (j2-1) * gh, 
				     zmin + (j3-1) * gh);
			  pointnr.Elem(j) = mesh.AddPoint (hp);
			  frontpointnr.Elem(j) =
			    AddPoint (hp, pointnr.Elem(j));

			}
		    }
	    }
	}


  for (i1 = 2; i1 <= n1-1; i1++)
    for (i2 = 2; i2 <= n2-1; i2++)
      for (i3 = 2; i3 <= n3-1; i3++)
	{
	  i = i3 + (i2-1) * n3 + (i1-1) * n2 * n3;
	  if (inner.Elem(i) == BLOCKINNER)
	    {
	      int pn[9];
	      pn[1] = pointnr.Get(i);
	      pn[2] = pointnr.Get(i+1);
	      pn[3] = pointnr.Get(i+n3);
	      pn[4] = pointnr.Get(i+n3+1);
	      pn[5] = pointnr.Get(i+n2*n3);
	      pn[6] = pointnr.Get(i+n2*n3+1);
	      pn[7] = pointnr.Get(i+n2*n3+n3);
	      pn[8] = pointnr.Get(i+n2*n3+n3+1);
	      static int elind[][4] =
	      {
		{ 1, 8, 2, 4 },
		{ 1, 8, 4, 3 },
		{ 1, 8, 3, 7 },
		{ 1, 8, 7, 5 },
		{ 1, 8, 5, 6 },
		{ 1, 8, 6, 2 }
	      };
	      for (j = 1; j <= 6; j++)
		{
		  Element el(4);
		  for (int k = 1; k <= 4;  k++)
		    el.PNum(k) = pn[elind[j-1][k-1]];

		  mesh.AddVolumeElement (el);
		}
	    }
	}



  for (i1 = 2; i1 <= n1-1; i1++)
    for (i2 = 2; i2 <= n2-1; i2++)
      for (i3 = 2; i3 <= n3-1; i3++)
	{
	  i = i3 + (i2-1) * n3 + (i1-1) * n2 * n3;
	  if (inner.Elem(i) == BLOCKINNER)
	    {    
	      int pi1(0), pi2(0), pi3(0), pi4(0);

	      int pn1 = frontpointnr.Get(i);
	      int pn2 = frontpointnr.Get(i+1);
	      int pn3 = frontpointnr.Get(i+n3);
	      int pn4 = frontpointnr.Get(i+n3+1);
	      int pn5 = frontpointnr.Get(i+n2*n3);
	      int pn6 = frontpointnr.Get(i+n2*n3+1);
	      int pn7 = frontpointnr.Get(i+n2*n3+n3);
	      int pn8 = frontpointnr.Get(i+n2*n3+n3+1);

	      for (int k = 1; k <= 6; k++)
		{
		  switch (k)
		    {
		    case 1: // j3 = i3+1
		      j = i + 1;
		      pi1 = pn2;
		      pi2 = pn6;
		      pi3 = pn4;
		      pi4 = pn8;
		      break;
		    case 2: // j3 = i3-1
		      j = i - 1;
		      pi1 = pn1;
		      pi2 = pn3;
		      pi3 = pn5;
		      pi4 = pn7;
		      break;
		    case 3: // j2 = i2+1
		      j = i + n3;
		      pi1 = pn3;
		      pi2 = pn4;
		      pi3 = pn7;
		      pi4 = pn8;
		      break;
		    case 4: // j2 = i2-1
		      j = i - n3;
		      pi1 = pn1;
		      pi2 = pn5;
		      pi3 = pn2;
		      pi4 = pn6;
		      break;
		    case 5: // j1 = i1+1
		      j = i + n3*n2;
		      pi1 = pn5;
		      pi2 = pn7;
		      pi3 = pn6;
		      pi4 = pn8;
		      break;
		    case 6: // j1 = i1-1
		      j = i - n3*n2;
		      pi1 = pn1;
		      pi2 = pn2;
		      pi3 = pn3;
		      pi4 = pn4;
		      break;
		    }

		  if (inner.Get(j) == BLOCKBOUND)
		    {
		      MiniElement2d face;
		      face.PNum(1) = pi4;
		      face.PNum(2) = pi1;
		      face.PNum(3) = pi3;
		      AddBoundaryElement (face);

		      face.PNum(1) = pi1;
		      face.PNum(2) = pi4;
		      face.PNum(3) = pi2;
		      AddBoundaryElement (face);

		    }
		}
	    }
	}
}



static const AdFront3 * locadfront;
static int TestInner (const Point3d & p)
{
  return locadfront->Inside (p);
}
static int TestSameSide (const Point3d & p1, const Point3d & p2)
{
  return locadfront->SameSide (p1, p2);
}




void Meshing3 :: BlockFillLocalH (Mesh & mesh, 
				  const MeshingParameters & mp)
{
  int i, j;
  
  double filldist = mp.filldist;

  (*testout) << "blockfill local h" << endl;
  (*testout) << "rel filldist = " << filldist << endl;
  PrintMessage (3, "blockfill local h");

  /*  
  (*mycout) << "boxes: " << mesh.LocalHFunction().GetNBoxes() << endl
	    << "filldist = " << filldist << endl;
  */
  ARRAY<Point3d> npoints;
  
  adfront -> CreateTrees();

  Point3d mpmin, mpmax;
  // mesh.GetBox (mpmin, mpmax);
  bool firstp = 1;

  double maxh = 0;
  for (i = 1; i <= adfront->GetNF(); i++)
    {
      const MiniElement2d & el = adfront->GetFace(i);
      for (j = 1; j <= 3; j++)
	{
	  const Point3d & p1 = adfront->GetPoint (el.PNumMod(j));
	  const Point3d & p2 = adfront->GetPoint (el.PNumMod(j+1));
	  double hi = Dist (p1, p2);
	  if (hi > maxh)
	    {
	      maxh = hi;
	      //(*testout) << "reducing maxh to " << maxh << " because of " << p1 << " and " << p2 << endl;
	    }

	  if (firstp)
	    {
	      mpmin = p1;
	      mpmax = p1;
	      firstp = 0;
	    }
	  else
	    {
	      mpmin.SetToMin  (p1);
	      mpmax.SetToMax  (p1);
	    }
	}
    }

  Point3d mpc = Center (mpmin, mpmax);
  double d = max3(mpmax.X()-mpmin.X(), 
		  mpmax.Y()-mpmin.Y(), 
		  mpmax.Z()-mpmin.Z()) / 2;
  mpmin = mpc - Vec3d (d, d, d);
  mpmax = mpc + Vec3d (d, d, d);
  Box3d meshbox (mpmin, mpmax);

  LocalH loch2 (mpmin, mpmax, 1);


  if (mp.maxh < maxh)
    {
      maxh = mp.maxh;
      //(*testout) << "reducing maxh to " << maxh << " because of mp.maxh" << endl;
    }

  int changed;
  do 
    {
      mesh.LocalHFunction().ClearFlags();

      for (i = 1; i <= adfront->GetNF(); i++)
	{
	  const MiniElement2d & el = adfront->GetFace(i);
	  Point3d pmin = adfront->GetPoint (el.PNum(1));
	  Point3d pmax = pmin;
	  
	  for (j = 2; j <= 3; j++)
	    {
	      const Point3d & p = adfront->GetPoint (el.PNum(j));
	      pmin.SetToMin (p);
	      pmax.SetToMax (p);
	    }
	  

	  double filld = filldist * Dist (pmin, pmax);
	  
	  pmin = pmin - Vec3d (filld, filld, filld);
	  pmax = pmax + Vec3d (filld, filld, filld);
	  //	  (*testout) << "cut : " << pmin << " - " << pmax << endl;
	  mesh.LocalHFunction().CutBoundary (pmin, pmax);
	}

      locadfront = adfront;
      mesh.LocalHFunction().FindInnerBoxes (adfront, NULL);

      npoints.SetSize(0);
      mesh.LocalHFunction().GetInnerPoints (npoints);

      changed = 0;
      for (i = 1; i <= npoints.Size(); i++)
	{
	  if (mesh.LocalHFunction().GetH(npoints.Get(i)) > 1.5 * maxh)
	    {
	      mesh.LocalHFunction().SetH (npoints.Get(i), maxh);
	      changed = 1;
	    }
	}
    }
  while (changed);

  if (debugparam.slowchecks)
    (*testout) << "Blockfill with points: " << endl;
  for (i = 1; i <= npoints.Size(); i++)
    {
      if (meshbox.IsIn (npoints.Get(i)))
	{
	  int gpnum = mesh.AddPoint (npoints.Get(i));
	  adfront->AddPoint (npoints.Get(i), gpnum);

	  if (debugparam.slowchecks)
	    {
	      (*testout) << npoints.Get(i) << endl;
	      if (!adfront->Inside(npoints.Get(i)))
		{
		  cout << "add outside point" << endl;
		  (*testout) << "outside" << endl;
		}
	    }

	}
    }

  

  // find outer points
  
  loch2.ClearFlags();

  for (i = 1; i <= adfront->GetNF(); i++)
    {
      const MiniElement2d & el = adfront->GetFace(i);
      Point3d pmin = adfront->GetPoint (el.PNum(1));
      Point3d pmax = pmin;
      
      for (j = 2; j <= 3; j++)
	{
	  const Point3d & p = adfront->GetPoint (el.PNum(j));
	  pmin.SetToMin (p);
	  pmax.SetToMax (p);
	}
      
      loch2.SetH (Center (pmin, pmax), Dist (pmin, pmax));
    }

  for (i = 1; i <= adfront->GetNF(); i++)
    {
      const MiniElement2d & el = adfront->GetFace(i);
      Point3d pmin = adfront->GetPoint (el.PNum(1));
      Point3d pmax = pmin;
      
      for (j = 2; j <= 3; j++)
	{
	  const Point3d & p = adfront->GetPoint (el.PNum(j));
	  pmin.SetToMin (p);
	  pmax.SetToMax (p);
	}
      
      double filld = filldist * Dist (pmin, pmax);
      pmin = pmin - Vec3d (filld, filld, filld);
      pmax = pmax + Vec3d (filld, filld, filld);
      loch2.CutBoundary (pmin, pmax);
    }

  locadfront = adfront;
  loch2.FindInnerBoxes (adfront, NULL);

  npoints.SetSize(0);
  loch2.GetOuterPoints (npoints);
  
  for (i = 1; i <= npoints.Size(); i++)
    {
      if (meshbox.IsIn (npoints.Get(i)))
	{
	  int gpnum = mesh.AddPoint (npoints.Get(i));
	  adfront->AddPoint (npoints.Get(i), gpnum);
	}
    }  
}

}
