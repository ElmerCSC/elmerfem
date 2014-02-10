#include <mystdlib.h>
#include <myadt.hpp>

#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>


namespace netgen
{
Identification :: Identification (int anr, const CSGeometry & ageom)
  : geom(ageom), identfaces(10)
{
  nr = anr;
}

Identification :: ~Identification ()
{
  ;
}


ostream & operator<< (ostream & ost, Identification & ident)
{
  ident.Print (ost);
  return ost;
}


/*
void Identification :: IdentifySpecialPoints (ARRAY<class SpecialPoint> & points)
{
  ;
}
*/


int Identification :: 
Identifyable (const SpecialPoint & sp1, const SpecialPoint & sp2,
	      const TABLE<int> & specpoint2solid,
	      const TABLE<int> & specpoint2surface) const
{
  cout << "Identification::Identifyable called for base-class" << endl;
  return 0;
}

int Identification :: 
Identifyable (const Point<3> & p1, const Point<3> & sp2) const
{
  cout << "Identification::Identifyable called for base-class" << endl;
  return 0;
}


int Identification :: 
IdentifyableCandidate (const SpecialPoint & sp1) const
{
  return 1;
}


int Identification :: 
ShortEdge (const SpecialPoint & sp1, const SpecialPoint & sp2) const
{
  return 0;
}

int Identification :: GetIdentifiedPoint (class Mesh & mesh, int pi)
{
  cout << "Identification::GetIdentifiedPoint called for base-class" << endl;
  return -1;
}

void Identification :: IdentifyPoints (Mesh & mesh)
{
  cout << "Identification::IdentifyPoints called for base-class" << endl;
  ;
}

void Identification :: IdentifyFaces (class Mesh & mesh)
{
  cout << "Identification::IdentifyFaces called for base-class" << endl;
  ;
}

void Identification :: 
BuildSurfaceElements (ARRAY<Segment> & segs,
		      Mesh & mesh, const Surface * surf)
{
  cout << "Identification::BuildSurfaceElements called for base-class" << endl;
  ;
}


void Identification :: 
BuildVolumeElements (ARRAY<class Element2d> & surfels,
			  class Mesh & mesh)
{
  ;
}

void Identification :: 
GetIdentifiedFaces (ARRAY<INDEX_2> & idfaces) const
{
  idfaces.SetSize(0);
  for (int i = 1; i <= identfaces.GetNBags(); i++)
    for (int j = 1; j <= identfaces.GetBagSize(i); j++)
      {
	INDEX_2 i2;
	int val;
	identfaces.GetData (i, j, i2, val);
	idfaces.Append (i2);
      }
}




PeriodicIdentification ::
PeriodicIdentification (int anr,
			const CSGeometry & ageom,
			const Surface * as1,
			const Surface * as2)
  : Identification(anr, ageom)
{
  s1 = as1;
  s2 = as2;
}

PeriodicIdentification :: ~PeriodicIdentification ()
{
  ;
}

/*
void PeriodicIdentification :: IdentifySpecialPoints 
(ARRAY<class SpecialPoint> & points)
{
  int i, j;
  int bestj;
  double bestval, val;

  for (i = 1; i <= points.Size(); i++)
    {
      Point<3> p1 = points.Get(i).p;
      Point<3> hp1 = p1;
      s1->Project (hp1);
      if (Dist (p1, hp1) > 1e-6) continue;

      Vec<3> n1;
      s1->GetNormalVector (p1, n1);
      n1 /= n1.Length();
      if ( fabs(n1 * points.Get(i).v) > 1e-3)
	continue;

      bestval = 1e8;
      bestj = 1;
      for (j = 1; j <= points.Size(); j++)
	{
	  Point<3> p2= points.Get(j).p;
	  Point<3> hp2 = p2;
	  s2->Project (hp2);
	  if (Dist (p2, hp2) > 1e-6) continue;
	  
	  Vec<3> n2;
	  s2->GetNormalVector (p2, n2);
	  n2 /= n2.Length();
	  if ( fabs(n2 * points.Get(j).v) > 1e-3)
	    continue;


	  Vec<3> v(p1, p2);
	  double vl = v.Length();
	  double cl = fabs (v*n1);

	  val = 1 - cl*cl/(vl*vl);

	  val += (points.Get(i).v - points.Get(j).v).Length();

	  if (val < bestval)
	    {
	      bestj = j;
	      bestval = val;
	    }
	}

      (*testout) << "Identify Periodic special points: pi = " 
		 << points.Get(i).p << ", vi = " << points.Get(i).v 
		 << " pj = " << points.Get(bestj).p 
		 << ", vj = " << points.Get(bestj).v 
		 << " bestval = " << bestval << endl;
    }
}
*/

int PeriodicIdentification :: 
Identifyable (const SpecialPoint & sp1, const SpecialPoint & sp2,
	      const TABLE<int> & specpoint2solid,
	      const TABLE<int> & specpoint2surface) const
{
  int i;
  double val;
  
  SpecialPoint hsp1 = sp1;
  SpecialPoint hsp2 = sp2;

  for (i = 1; i <= 1; i++)
    {
      //      Swap (hsp1, hsp2);

      if (!s1->PointOnSurface (hsp1.p))
	continue;

      Vec<3> n1;
      n1 = s1->GetNormalVector (hsp1.p);
      n1 /= n1.Length();
      if ( fabs(n1 * hsp1.v) > 1e-3)
	continue;


      if (!s2->PointOnSurface(hsp2.p))
	continue;

      Vec<3> n2;
      n2 = s2->GetNormalVector (hsp2.p);
      n2 /= n2.Length();
      if ( fabs(n2 * hsp2.v) > 1e-3)
	continue;


      Vec<3> v = hsp2.p - hsp1.p;
      double vl = v.Length();
      double cl = fabs (v*n1);
      

      val = 1 - cl*cl/(vl*vl);
      val += (hsp1.v - hsp2.v).Length();
    
      if (val < 1e-6)
        return 1;
    }

  return 0;
}

int PeriodicIdentification :: 
Identifyable (const Point<3> & p1, const Point<3> & p2) const
{
  return (s1->PointOnSurface (p1) &&
	  s2->PointOnSurface (p2));
}
  



int PeriodicIdentification :: 
GetIdentifiedPoint (class Mesh & mesh,  int pi)
{
  const Surface * sold, *snew;
  const Point<3> & p = mesh.Point (pi);

  if (s1->PointOnSurface (p))
    {
      snew = s2;
    }
  else
    {
      if (s2->PointOnSurface (p))
	{
	  snew = s1;
	}
      else
	{
	  cerr << "GetIdenfifiedPoint: Not possible" << endl;
	  exit (1);
	}    
    }
  
  // project to other surface
  Point<3> hp = p;
  snew->Project (hp);

  int i;
  int newpi = 0;
  for (i = 1; i <= mesh.GetNP(); i++)
    if (Dist2 (mesh.Point(i), hp) < 1e-12)
      {
	newpi = i;
	break;
      }
  if (!newpi)
    newpi = mesh.AddPoint (hp);

  if (snew == s2)
    mesh.GetIdentifications().Add (pi, newpi, nr);
  else
    mesh.GetIdentifications().Add (newpi, pi, nr);

  mesh.GetIdentifications().SetType(nr,Identifications::PERIODIC);
	   
  /* 
  (*testout) << "Identify points(periodic), nr = " << nr << ": " << mesh.Point(pi)
	     << " and " << mesh.Point(newpi) 
	     << ((snew == s2) ? "" : " inverse")
	     << endl;
  */
  return newpi;
}


void PeriodicIdentification :: IdentifyPoints (class Mesh & mesh)
{
  int i, j;
  for (i = 1; i <= mesh.GetNP(); i++)
    {
      Point<3> p = mesh.Point(i);
      if (s1->PointOnSurface (p))
	{
	  Point<3> pp = p;
	  s2->Project (pp);
	  for (j = 1; j <= mesh.GetNP(); j++)
	    if (Dist2(mesh.Point(j), pp) < 1e-6)
	      {
		mesh.GetIdentifications().Add (i, j, nr);
		/*
		(*testout) << "Identify points(periodic:), nr = " << nr << ": "
			   << mesh.Point(i) << " - " << mesh.Point(j) << endl;
		*/
	      }
	}
    }

  mesh.GetIdentifications().SetType(nr,Identifications::PERIODIC);
}


void PeriodicIdentification :: IdentifyFaces (class Mesh & mesh)
{
  int i, j, k, l;
  int fi1, fi2, side;
  for (i = 1; i <= mesh.GetNFD(); i++)
    for (j = 1; j <= mesh.GetNFD(); j++)
      {
	int surfi = mesh.GetFaceDescriptor(i).SurfNr();
	int surfj = mesh.GetFaceDescriptor(j).SurfNr();
	if (surfi == surfj)
	  continue;
	
	if (geom.GetSurface (surfi) != s1 ||
	    geom.GetSurface (surfj) != s2)
	  continue;
	    
	int idok = 1;


	//	(*testout) << "check faces " << i << " and " << j << endl;
	for (side = 1; side <= 2 && idok; side++)
	  {
	    if (side == 1)
	      {
		fi1 = i; 
		fi2 = j;
	      }
	    else
	      {
		fi1 = j;
		fi2 = i;
	      }

	    for (k = 1; k <= mesh.GetNSeg(); k++)
	      {
		const Segment & seg1 = mesh.LineSegment(k);
		if (seg1.si != fi1)
		  continue;

		int foundother = 0;
		for (l = 1; l <= mesh.GetNSeg(); l++)
		  {
		    const Segment & seg2 = mesh.LineSegment(l);
		    if (seg2.si != fi2)
		      continue;
		    
		    //		    (*testout) << "seg1 = " << seg1.p1 << "-" << seg1.p2 << ", seg2 = " << seg2.p1 << "-" << seg2.p2;

		    if (side == 1)
		      {
			if (mesh.GetIdentifications().Get (seg1.p1, seg2.p1) &&
			    mesh.GetIdentifications().Get (seg1.p2, seg2.p2))
			  {
			    foundother = 1;
			    break;
			  }
			
			if (mesh.GetIdentifications().Get (seg1.p1, seg2.p2) &&
			    mesh.GetIdentifications().Get (seg1.p2, seg2.p1))
			  {
			    foundother = 1;
			    break;
			  }
		      }
		    else
		      {
			if (mesh.GetIdentifications().Get (seg2.p1, seg1.p1) &&
			    mesh.GetIdentifications().Get (seg2.p2, seg1.p2))
			  {
			    foundother = 1;
			    break;
			  }
			
			if (mesh.GetIdentifications().Get (seg2.p1, seg1.p2) &&
			    mesh.GetIdentifications().Get (seg2.p2, seg1.p1))
			  {
			    foundother = 1;
			    break;
			  }
		      }
		  }

		if (!foundother)
		  {
		    idok = 0;
		    break;
		  }
	      }
	  }


	if (idok)
	  {
	    // (*testout) << "Identify faces " << i << " and " << j << endl;
	    INDEX_2 fpair(i,j);
	    fpair.Sort();
	    identfaces.Set (fpair, 1);
	  }
      }
}



void PeriodicIdentification :: 
BuildSurfaceElements (ARRAY<Segment> & segs,
		      Mesh & mesh, const Surface * surf)
{
  int found = 0;
  int fother;

  int facei = segs.Get(1).si;
  int surfnr = mesh.GetFaceDescriptor(facei).SurfNr();

  if (geom.GetSurface(surfnr) == s1 ||
      geom.GetSurface(surfnr) == s2)
    {
      ARRAY<int> copy_points;

      for (int i = 1; i <= mesh.GetNSE(); i++)
	{
	  const Element2d & sel = mesh.SurfaceElement(i);
	  INDEX_2 fpair (facei, sel.GetIndex());
	  fpair.Sort();
	  if (identfaces.Used (fpair))
            {
	      for (int k = 0; k < sel.GetNP(); k++)
                if (!copy_points.Contains (sel[k]))
                  copy_points.Append (sel[k]);
            }      
        }
      BubbleSort (copy_points);
      for (int k = 0; k < copy_points.Size(); k++)
        GetIdentifiedPoint (mesh, copy_points[k]);



      for (int i = 1; i <= mesh.GetNSE(); i++)
	{
	  const Element2d & sel = mesh.SurfaceElement(i);
	  INDEX_2 fpair (facei, sel.GetIndex());
	  fpair.Sort();
	  if (identfaces.Used (fpair))
	    {
	      found = 1;
	      fother = sel.GetIndex();

	      // copy element
	      Element2d newel(sel.GetType());
	      newel.SetIndex (facei);
	      for (int k = 0; k < sel.GetNP(); k++)
                newel[k] = GetIdentifiedPoint (mesh, sel[k]);

	      Vec<3> nt = Cross (Point<3> (mesh[newel[1]])- Point<3> (mesh[newel[0]]),
				 Point<3> (mesh[newel[2]])- Point<3> (mesh[newel[0]]));
	      
	      Vec<3> nsurf = geom.GetSurface (surfnr)->GetNormalVector (mesh[newel[0]]);
	      if (nsurf * nt < 0)
                Swap (newel[0], newel[2]);
				
	      mesh.AddSurfaceElement (newel);
	    }
	}
    }
  
  if (found)
    {
      // (*mycout) << " copy face " << facei << " from face " << fother;
      PrintMessage (4, " copy face ", facei, " from face ", fother);
      
      segs.SetSize(0);
    }
}








void PeriodicIdentification :: Print (ostream & ost) const
{
  ost << "Periodic Identifiaction, surfaces: " 
      << s1->Name() << " - " << s2->Name() << endl;
  s1->Print (ost);
  ost << " - ";
  s2->Print (ost);
  ost << endl;
}


void PeriodicIdentification :: GetData (ostream & ost) const
{
  ost << "periodic " << s1->Name() << " " << s2->Name();
}







CloseSurfaceIdentification ::
CloseSurfaceIdentification (int anr,
			    const CSGeometry & ageom,
			    const Surface * as1,
			    const Surface * as2,
			    const TopLevelObject * adomain,
			    const Flags & flags)
  : Identification(anr, ageom)
{
  s1 = as1;
  s2 = as2;
  domain = adomain;
  ref_levels = int (flags.GetNumFlag ("reflevels", 2));
  ref_levels_s1 = int (flags.GetNumFlag ("reflevels1", 0));
  ref_levels_s2 = int (flags.GetNumFlag ("reflevels2", 0));
  slices = flags.GetNumListFlag ("slices");
  for(int i=0; i<slices.Size(); i++)
    if((i==0 && slices[i] <= 0) ||
       (i>0 && slices[i] <= slices[i-1]) ||
       (slices[i] >= 1))
      throw NgException ("slices have to be in ascending order, between 0 and 1");

  // eps_n = 1e-3;
  eps_n = 1e-6;

  dom_surf_valid = 0;

  if (domain)
    for (int i = 0; i < geom.GetNTopLevelObjects(); i++)
      if (domain == geom.GetTopLevelObject(i))
	dom_nr = i;

  usedirection = flags.NumListFlagDefined("direction");
  if(usedirection)
    {
      for(int i=0; i<3; i++)
	direction(i) = flags.GetNumListFlag("direction")[i];

      direction.Normalize();
    }
}

CloseSurfaceIdentification :: ~CloseSurfaceIdentification ()
{
  ;
}

void CloseSurfaceIdentification :: Print (ostream & ost) const
{
  ost << "CloseSurface Identifiaction, surfaces: " 
      << s1->Name() << " - " << s2->Name() << endl;
  s1->Print (ost);
  s2->Print (ost);
  ost << endl;
}


void CloseSurfaceIdentification :: GetData (ostream & ost) const
{
  ost << "close surface " << s1->Name() << " " << s2->Name();
}


/*
void CloseSurfaceIdentification :: IdentifySpecialPoints 
(ARRAY<class SpecialPoint> & points)
{
  int i, j;
  int bestj;
  double bestval, val;

  for (i = 1; i <= points.Size(); i++)
    {
      Point<3> p1 = points.Get(i).p;
      Vec<3> n1;

      if (!s1->PointOnSurface (p1))
	continue;

	s1->GetNormalVector (p1, n1);
      n1 /= n1.Length();
      if ( fabs(n1 * points.Get(i).v) > 1e-3)
	continue;

      bestval = 1e8;
      bestj = 1;
      for (j = 1; j <= points.Size(); j++)
	{
	  Point<3> p2= points.Get(j).p;
	  if (!s2->PointOnSurface (p2))
	    continue;
	  
	  Vec<3> n2;
	  s2->GetNormalVector (p2, n2);
	  n2 /= n2.Length();
	  if ( fabs(n2 * points.Get(j).v) > 1e-3)
	    continue;


	  Vec<3> v(p1, p2);
	  double vl = v.Length();
	  double cl = fabs (v*n1);

	  val = 1 - cl*cl/(vl*vl);

	  val += (points.Get(i).v - points.Get(j).v).Length();

	  if (val < bestval)
	    {
	      bestj = j;
	      bestval = val;
	    }
	}

      (*testout) << "Identify close surfaces special points: pi = " 
		 << points.Get(i).p << ", vi = " << points.Get(i).v 
		 << " pj = " << points.Get(bestj).p 
		 << ", vj = " << points.Get(bestj).v 
		 << " bestval = " << bestval << endl;
    }
}
*/

int CloseSurfaceIdentification :: 
Identifyable (const SpecialPoint & sp1, const SpecialPoint & sp2,
	      const TABLE<int> & specpoint2solid,
	      const TABLE<int> & specpoint2surface) const
{
  //(*testout) << "identcheck: " << sp1.p << "; " << sp2.p << endl;

  if (!dom_surf_valid)
    {
      const_cast<bool&> (dom_surf_valid) = 1;
      ARRAY<int> & hsurf = const_cast<ARRAY<int>&> (domain_surfaces);

      if (domain)
	{
	  BoxSphere<3> hbox (geom.BoundingBox());
	  geom.GetIndependentSurfaceIndices (domain->GetSolid(), hbox, hsurf);
	  //(*testout) << "surfs of identification " << nr << ": " << endl << hsurf << endl;
	}
      else
	{
	  hsurf.SetSize (geom.GetNSurf());
	  for (int j = 0; j < hsurf.Size(); j++)
	    hsurf[j] = j;
	}
    }

  if (domain)
    {
      bool has1 = 0, has2 = 0;
      for (int i = 0; i < specpoint2solid[sp1.nr].Size(); i++)
	if (specpoint2solid[sp1.nr][i] == dom_nr)
	  { has1 = 1; break; }
      for (int i = 0; i < specpoint2solid[sp2.nr].Size(); i++)
	if (specpoint2solid[sp2.nr][i] == dom_nr)
	  { has2 = 1; break; }

      if (!has1 || !has2) 
	{
	  //(*testout) << "failed at pos1" << endl;
	  return 0;
	}
    }

  if (!s1->PointOnSurface (sp1.p))
    {
      //(*testout) << "failed at pos2" << endl;
      return 0;
    }

//   (*testout) << "sp1 " << sp1.p << " sp2 " << sp2.p << endl
// 	     << "specpoint2solid[sp1.nr] " << specpoint2solid[sp1.nr] << endl
// 	     << "specpoint2solid[sp2.nr] " << specpoint2solid[sp2.nr] << endl;


  Vec<3> n1 = s1->GetNormalVector (sp1.p);
  n1.Normalize();
  if ( fabs(n1 * sp1.v) > eps_n)
    {
      //(*testout) << "failed at pos3" << endl;
      return 0;
    }

  if (!s2->PointOnSurface(sp2.p))
    {
      //(*testout) << "failed at pos4" << endl;
      return 0;
    }


  Vec<3> n2 = s2->GetNormalVector (sp2.p);
  n2.Normalize();
  if ( fabs(n2 * sp2.v) > eps_n)
    {
      //(*testout) << "failed at pos5" << endl;
      return 0;
    }
  
  // must have joint surface 
  bool joint = 0;

  int j = 0, k = 0;
  while (1)
    {
      int snr1 = specpoint2surface[sp1.nr][j];
      int snr2 = specpoint2surface[sp2.nr][k];
      if (snr1 < snr2) 
	{
	  j++;
	  if (j == specpoint2surface[sp1.nr].Size()) break;
	}
      else if (snr2 < snr1) 
	{
	  k++;
	  if (k == specpoint2surface[sp2.nr].Size()) break;
	}
      else
	{
	  bool dom_surf = 0;
	  for (int l = 0; l < domain_surfaces.Size(); l++)
	    if (domain_surfaces[l] == snr1)
	      dom_surf = 1;

	  if (dom_surf)
	    {
	      Vec<3> hn1 = geom.GetSurface(snr1)->GetNormalVector (sp1.p);
	      Vec<3> hn2 = geom.GetSurface(snr1)->GetNormalVector (sp2.p);
	      
	      if (hn1 * hn2 > 0)
		{
		  joint = 1;
		  break;
		}
	    }

	  j++;
	  if (j == specpoint2surface[sp1.nr].Size()) break;
	  k++;
	  if (k == specpoint2surface[sp2.nr].Size()) break;
	}
    }

  if (!joint)
    {
      //(*testout) << "failed at pos6" << endl;
      return 0;
    }

  Vec<3> v = sp2.p - sp1.p;
  double vl = v.Length();
  double cl = (usedirection) ? fabs(v*direction) : fabs (v*n1);


  if(cl <= (1-eps_n*eps_n) * vl)
    {
      //(*testout) << "failed at pos7" << endl;
      return 0;
    }
  
  double dl;

  if(usedirection)
    {
      Vec<3> v1 = sp1.v - (sp1.v*direction)*direction; v1.Normalize();
      Vec<3> v2 = sp2.v - (sp2.v*direction)*direction; v2.Normalize();
      
      dl = (v1 - v2).Length();
    }
  else
    dl = (sp1.v - sp2.v).Length();

  if (dl < 0.1)
    return 1;
   

  //(*testout) << "failed at pos8" << endl;
  return 0;
}

int CloseSurfaceIdentification :: 
Identifyable (const Point<3> & p1, const Point<3> & p2) const
{  
//   if (domain)
//     if (!domain->GetSolid()->IsIn (p1) || !domain->GetSolid()->IsIn (p2))
//       return 0;
  return (s1->PointOnSurface (p1) && s2->PointOnSurface (p2));
}
  



int CloseSurfaceIdentification :: 
IdentifyableCandidate (const SpecialPoint & sp1) const
{
  if (domain)
    if (!domain->GetSolid()->IsIn (sp1.p))
      return 0;

  if (s1->PointOnSurface (sp1.p))
    {
      Vec<3> n1;
      n1 = s1->GetNormalVector (sp1.p);
      n1.Normalize();
      if ( fabs(n1 * sp1.v) > eps_n)
	return 0;
      return 1;
    }

  if (s2->PointOnSurface (sp1.p))
    {
      Vec<3> n1;
      n1 = s2->GetNormalVector (sp1.p);
      n1.Normalize();
      if ( fabs(n1 * sp1.v) > eps_n)
	return 0;
      return 1;
    }
  return 0;
}



int CloseSurfaceIdentification :: 
ShortEdge (const SpecialPoint & sp1, const SpecialPoint & sp2) const
{  
  if ( (s1->PointOnSurface (sp1.p) && s2->PointOnSurface (sp2.p)) ||
       (s1->PointOnSurface (sp2.p) && s2->PointOnSurface (sp1.p)) )
    {
      return 1;
    }
  return 0;
}



int CloseSurfaceIdentification :: 
GetIdentifiedPoint (class Mesh & mesh,  int pi)
{
  const Surface * sold, *snew;
  const Point<3> & p = mesh.Point (pi);

  ARRAY<int,PointIndex::BASE> identmap(mesh.GetNP());
  mesh.GetIdentifications().GetMap (nr, identmap);
  if (identmap.Get(pi))
    return identmap.Get(pi);

  
  if (s1->PointOnSurface (p))
    snew = s2;
  else if (s2->PointOnSurface (p))
    snew = s1;
  else
    {
      (*testout)  << "GetIdenfifiedPoint: Not possible" << endl;
      (*testout) << "p = " << p << endl;
      (*testout) << "surf1: " << (*s1) << endl
		 << "surf2: " << (*s2) << endl;
      
      cerr << "GetIdenfifiedPoint: Not possible" << endl;
      throw NgException ("GetIdenfifiedPoint: Not possible");
    }    

  // project to other surface
  Point<3> hp = p;
  if(usedirection)
    snew->SkewProject(hp,direction);
  else
    snew->Project (hp);

  //(*testout) << "projecting " << p << " to " << hp << endl;

  int newpi = 0;
  for (int i = 1; i <= mesh.GetNP(); i++)
    if (Dist2 (mesh.Point(i), hp) < 1e-12)
      //    if (Dist2 (mesh.Point(i), hp) < 1 * Dist2 (hp, p))
      {
	newpi = i;
	break;
      }
  if (!newpi)
    newpi = mesh.AddPoint (hp);

  if (snew == s2)
    {
      mesh.GetIdentifications().Add (pi, newpi, nr);
      //(*testout) << "add identification(1) " << pi << " - " << newpi << ", " << nr << endl;
    }
  else
    {
      mesh.GetIdentifications().Add (newpi, pi, nr);
      //(*testout) << "add identification(2) " << newpi << " - " << pi << ", " << nr << endl;
    }
  mesh.GetIdentifications().SetType(nr,Identifications::CLOSESURFACES);
	   

  /*
  (*testout) << "Identify points(closesurface), nr = " << nr << ": " << mesh.Point(pi)
	     << " and " << mesh.Point(newpi) 
	     << ((snew == s2) ? "" : " inverse")
	     << endl;
  */
  return newpi;
}





void CloseSurfaceIdentification :: IdentifyPoints (Mesh & mesh)
{
  int np = mesh.GetNP();

  ARRAY<int> points_on_surf2;

  for (int i2 = 1; i2 <= np; i2++)
    if (s2->PointOnSurface (mesh.Point(i2)))
      points_on_surf2.Append (i2);
    
  ARRAY<int> surfs_of_p1;

  for (int i1 = 1; i1 <= np; i1++)
    {
      Point<3> p1 = mesh.Point(i1);
      //      (*testout) << "p1 = " << i1 << " = " << p1 << endl;
      if (domain && !domain->GetSolid()->IsIn (p1))
	continue;
      
      //if(domain) (*testout) << "p1 is in " << domain->GetSolid()->Name() << endl;

      if (s1->PointOnSurface (p1))
	{
	  int candi2 = 0;
	  double mindist = 1e10;

	  Vec<3> n1;
	  n1 = s1->GetNormalVector (p1);
	  n1.Normalize();
	   
	  surfs_of_p1.SetSize(0);
	  for (int jj = 0; jj < domain_surfaces.Size(); jj++)
	    {
	      int j = domain_surfaces[jj];
	      if (geom.GetSurface(j) -> PointOnSurface(p1))
		surfs_of_p1.Append (j);
	    }
	  //(*testout) << " surfs of p1 = " << endl << surfs_of_p1 << endl;

	  for (int ii2 = 0; ii2 < points_on_surf2.Size(); ii2++)
	    {
	      int i2 = points_on_surf2[ii2];
	      if (i2 == i1) continue;
	      const Point<3> p2 = mesh.Point(i2);
	      
	      Vec<3> n = p2 - p1;
	      n.Normalize();
	      
	      bool joint = 0;
	      for (int jj = 0; jj < surfs_of_p1.Size(); jj++)
		{
		  int j = surfs_of_p1[jj];
		  if (geom.GetSurface(j) -> PointOnSurface(p2))
		    {
		      Vec<3> hn1 = geom.GetSurface(j)->GetNormalVector (p1);
		      Vec<3> hn2 = geom.GetSurface(j)->GetNormalVector (p2);
		      
		      if (hn1 * hn2 > 0)
			{
			  joint = 1;
			  break;
			}
		    }
		}

	      if (!joint) continue;
	      
	      if(usedirection)
		{
		  if (fabs (n*direction) > 0.9)
		    {
		      Vec<3> p1p2 = p2-p1;
		      double ndist = p1p2.Length2() - pow(p1p2*direction,2);
		      if(ndist < mindist)
			{
			  candi2 = i2;
			  mindist = ndist;
			}
		    }
		      
		}
	      else
		{
		  if (fabs (n * n1) > 0.9 &&
		      Dist (p1, p2) < mindist)
		    {
		      candi2 = i2;
		      mindist = Dist (p1, p2);
		    }
		}
	    
	    }

	  if (candi2)
	    {
	      //(*testout) << "identify points " << p1 << " - " << mesh.Point(candi2) << endl;

	      /*
	      (*testout) << "Add Identification from CSI2, nr = " << nr << ", p1 = " 
			 << i1 << " = " 
			 << mesh[PointIndex(i1)] << ", p2 = " << candi2 << " = " 
			 << mesh[PointIndex(candi2)] << endl;
	      */
	      mesh.GetIdentifications().Add (i1, candi2, nr);
	      mesh.GetIdentifications().SetType(nr,Identifications::CLOSESURFACES);
	      //(*testout) << "add identification " << i1 << " - " << candi2 << ", " << nr << endl;
	    }
	}
    }
}



void CloseSurfaceIdentification :: IdentifyFaces (class Mesh & mesh)
{
  int fi1, fi2, side;
  int s1rep, s2rep;

  for (int i = 0; i < geom.GetNSurf(); i++)
    {
      if (geom.GetSurface (i) == s1) 
	s1rep = geom.GetSurfaceClassRepresentant(i);
      if (geom.GetSurface (i) == s2) 
	s2rep = geom.GetSurfaceClassRepresentant(i);
    }

  ARRAY<int> segs_on_face1, segs_on_face2;

  identfaces.DeleteData();

  //(*testout) << "identify faces, nr = " << nr << endl;
  
  for (int i = 1; i <= mesh.GetNFD(); i++)
    {
      int surfi = mesh.GetFaceDescriptor(i).SurfNr();
      if (s1rep != surfi) continue;


      if (domain &&
	  domain != geom.GetTopLevelObject (mesh.GetFaceDescriptor(i).DomainIn()-1) &&
	  domain != geom.GetTopLevelObject (mesh.GetFaceDescriptor(i).DomainOut()-1))
	continue;

      for (int j = 1; j <= mesh.GetNFD(); j++)
	{
	  int surfj = mesh.GetFaceDescriptor(j).SurfNr();

	  if (surfi == surfj) continue;
	  if (s2rep != surfj) continue;

  
	  int idok = 1;
	  
	  for (side = 1; side <= 2 && idok; side++)
	    {
	      if (side == 1)
		{
		  fi1 = i; 
		  fi2 = j;
		}
	      else
		{
		  fi1 = j;
		  fi2 = i;
		}
	      

	      segs_on_face1.SetSize(0);
	      segs_on_face2.SetSize(0);

	      for (int k = 1; k <= mesh.GetNSeg(); k++)
		{
		  if (mesh.LineSegment(k).si == fi1)
		    segs_on_face1.Append (k);
		  if (mesh.LineSegment(k).si == fi2)
		    segs_on_face2.Append (k);
		}


	      for (int k = 1; k <= mesh.GetNSeg(); k++)
		{
		  const Segment & seg1 = mesh.LineSegment(k);
		  if (seg1.si != fi1)
		    continue;
		  
		  int foundother = 0;
		  /*
		  for (int l = 1; l <= mesh.GetNSeg(); l++)
		    {
		      const Segment & seg2 = mesh.LineSegment(l);
		      if (seg2.si != fi2)
			continue;
		  */
		  for (int ll = 0; ll < segs_on_face2.Size(); ll++)
		    {
		      int l = segs_on_face2[ll];
		      const Segment & seg2 = mesh.LineSegment(l);
		      
		      if (side == 1)
			{
			  if (mesh.GetIdentifications().Get (seg1.p1, seg2.p1) &&
			      mesh.GetIdentifications().Get (seg1.p2, seg2.p2))
			    {
			      foundother = 1;
			      break;
			    }
			  
			  if (mesh.GetIdentifications().Get (seg1.p1, seg2.p2) &&
			      mesh.GetIdentifications().Get (seg1.p2, seg2.p1))
			    {
			      foundother = 1;
			      break;
			    }
			}
		      else
			{
			  if (mesh.GetIdentifications().Get (seg2.p1, seg1.p1) &&
			      mesh.GetIdentifications().Get (seg2.p2, seg1.p2))
			    {
			      foundother = 1;
			      break;
			    }
			  
			  if (mesh.GetIdentifications().Get (seg2.p1, seg1.p2) &&
			      mesh.GetIdentifications().Get (seg2.p2, seg1.p1))
			    {
			      foundother = 1;
			      break;
			    }
			}
		    }
		  
		  if (!foundother)
		    {
		      idok = 0;
		      break;
		    }
		}
	    }
	  
	  
	  if (idok)
	    {
	      //(*testout) << "Identification " << nr << ", identify faces " << i << " and " << j << endl;
	      INDEX_2 fpair(i,j);
	      fpair.Sort();
	      identfaces.Set (fpair, 1);
	    }
	}
    }
}



void CloseSurfaceIdentification :: 
BuildSurfaceElements (ARRAY<Segment> & segs,
		      Mesh & mesh, const Surface * surf)
{
  bool found = 0;
  int cntquads = 0;

  ARRAY<int,PointIndex::BASE>  identmap;
  identmap = 0;

  mesh.GetIdentifications().GetMap (nr, identmap);
  
  for (int i = PointIndex::BASE; i < identmap.Size()+PointIndex::BASE; i++)
    if (identmap[i])  identmap[identmap[i]] = i;

    
  //(*testout) << "identification nr = " << nr << endl;
  //(*testout) << "surf = " << (*surf) << endl;
  //(*testout) << "domain = " << domain->GetSolid()->Name() << endl;
  //(*testout) << "segs = " << endl << segs << endl;
  //(*testout) << "identmap = " << endl << identmap << endl;
  
  //ARRAY<bool> foundseg(segs.Size());
  //foundseg = false;

  // insert quad layer:
  for (int i1 = 0; i1 < segs.Size(); i1++)
    {
      const Segment & s1 = segs[i1];
      if (identmap[s1.p1] && identmap[s1.p2])
	for (int i2 = 0; i2 < i1; i2++)
	  {
	    const Segment & s2 = segs[i2];
	    //(*testout) << "checking " << s1 << " and " << s2 << " for ident." << endl;

	    if(domain && !((s1.domin == dom_nr ||
			    s1.domout == dom_nr) &&
			   (s2.domin == dom_nr ||
			    s2.domout == dom_nr)))
	      continue;
	 
	    if ((mesh.GetIdentifications().Get (s1.p1, s2.p2, nr) && 
		 mesh.GetIdentifications().Get (s1.p2, s2.p1, nr))    || 
		(mesh.GetIdentifications().Get (s2.p1, s1.p2, nr) && 
		 mesh.GetIdentifications().Get (s2.p2, s1.p1, nr)))
	      {
		Element2d el(s1.p1, s1.p2, s2.p1, s2.p2);

		Vec<3> n = Cross (mesh[el[1]] - mesh[el[0]],
				  mesh[el[3]] - mesh[el[0]]);

		Vec<3> ns = surf->GetNormalVector (mesh[el[0]]);

		if (n * ns < 0)
		  {
		    Swap (el.PNum(1), el.PNum(2));
		    Swap (el.PNum(3), el.PNum(4));
		  }
			     
		mesh.AddSurfaceElement (el);
//  		(*testout) << "(id nr "<< nr <<") add rect element: "
//  			   << mesh.Point (el.PNum(1)) << " - "
//  			   << mesh.Point (el.PNum(2)) << " - "
//  			   << mesh.Point (el.PNum(3)) << " - "
//  			   << mesh.Point (el.PNum(4)) << endl;
		found = true;
		//foundseg[i1]=foundseg[i2] = true;
		cntquads++;
	      }
	  }
    }
  if (found)
    {
      PrintMessage(3, "insert quad layer of ", cntquads,
		   " elements at face ", segs.Get(1).si);
      //ARRAY<Segment> aux;
      //for(int i=0; i<segs.Size();i++)
      //	if(!foundseg[i])
      //	  aux.Append(segs[i]);
      segs.SetSize(0);
    }
  else
    {
      BuildSurfaceElements2 (segs, mesh, surf);
    }
}






void CloseSurfaceIdentification :: 
BuildSurfaceElements2 (ARRAY<Segment> & segs,
		       Mesh & mesh, const Surface * surf)
{
  // copy mesh


  //  (*testout) << "copy trig face, identnr = " << nr << endl;
  //  (*testout) << "identfaces = " << endl << identfaces << endl;

  if (!segs.Size()) return;

  bool found = 0;

  int fother;
  int facei = segs.Get(1).si;
  int surfnr = mesh.GetFaceDescriptor(facei).SurfNr();

  
  bool foundid = 0;
  for (INDEX_2_HASHTABLE<int>::Iterator it = identfaces.Begin();
       it != identfaces.End(); it++)
    {
      INDEX_2 i2;
      int data;
      identfaces.GetData (it, i2, data);
      if (i2.I1() == facei || i2.I2() == facei)
	foundid = 1;
    }

  /*
  for (int i = 1; i <= identfaces.GetNBags(); i++)
    for (int j = 1; j <= identfaces.GetBagSize(i); j++)
      {
	INDEX_2 i2;
	int data;
	identfaces.GetData (i, j, i2, data);
	if (i2.I1() == facei || i2.I2() == facei)
	  foundid = 1;

	(*testout) << "identface = " << i2 << endl;
	(*testout) << "face " << i2.I1() << " = " << mesh.GetFaceDescriptor(i2.I1()) << endl;
	(*testout) << "face " << i2.I2() << " = " << mesh.GetFaceDescriptor(i2.I2()) << endl;
      }
  */

  if (foundid)
    {
      //	  (*testout) << "surfaces found" << endl;
      // copy surface
      for (int i = 1; i <= mesh.GetNSE(); i++)
	{
	  const Element2d & sel = mesh.SurfaceElement(i);
	  INDEX_2 fpair (facei, sel.GetIndex());
	  fpair.Sort();
	  if (identfaces.Used (fpair))
	    {
	      found = 1;
	      fother = sel.GetIndex();
	      
	      // copy element
	      Element2d newel(sel.GetType());
	      newel.SetIndex (facei);
	      for (int k = 1; k <= sel.GetNP(); k++)
		{
		  newel.PNum(k) = 
		    GetIdentifiedPoint (mesh, sel.PNum(k));
		  //		      cout << "id-point = " << sel.PNum(k) << ", np = " << newel.PNum(k) << endl;
		}	  
	      
	      Vec<3> nt = Cross (Point<3> (mesh.Point (newel.PNum(2)))- 
				 Point<3> (mesh.Point (newel.PNum(1))),
				 Point<3> (mesh.Point (newel.PNum(3)))- 
				 Point<3> (mesh.Point (newel.PNum(1))));
	      Vec<3> nsurf;
	      nsurf = geom.GetSurface (surfnr)->GetNormalVector (mesh.Point(newel.PNum(1)));
	      if (nsurf * nt < 0)
		Swap (newel.PNum(2), newel.PNum(3));
	      
	      mesh.AddSurfaceElement (newel);
	    }
	}
    }
  
  if (found)
    {
      // (*mycout) << " copy face " << facei << " from face " << fother;
      PrintMessage (4, " copy face ", facei, " from face ", fother);
      segs.SetSize(0);
    }
}














void CloseSurfaceIdentification :: 
BuildVolumeElements (ARRAY<class Element2d> & surfels,
		     class Mesh & mesh)
{
  ;
}













/*   ***************** Close Edges Identification ********** */



CloseEdgesIdentification ::
CloseEdgesIdentification (int anr,
			  const CSGeometry & ageom,
			  const Surface * afacet,
			  const Surface * as1,
			  const Surface * as2)
  : Identification(anr, ageom)
{
  facet = afacet;
  s1 = as1;
  s2 = as2;
}

CloseEdgesIdentification :: ~CloseEdgesIdentification ()
{
  ;
}

void CloseEdgesIdentification :: Print (ostream & ost) const
{
  ost << "CloseEdges Identifiaction, facet = " 
      << facet->Name() << ", surfaces: " 
      << s1->Name() << " - " << s2->Name() << endl;
  facet->Print (ost);
  s1->Print (ost);
  s2->Print (ost);
  ost << endl;
}


void CloseEdgesIdentification :: GetData (ostream & ost) const
{
  ost << "closeedges " << facet->Name() << " " 
      << s1->Name() << " " << s2->Name();
}


/*
void CloseEdgesIdentification :: IdentifySpecialPoints 
(ARRAY<class SpecialPoint> & points)
{
  int i, j;
  int bestj;
  double bestval, val;

  for (i = 1; i <= points.Size(); i++)
    {
      Point<3> p1 = points.Get(i).p;
      Vec<3> n1;

      if (!s1->PointOnSurface (p1))
	continue;

	s1->GetNormalVector (p1, n1);
      n1 /= n1.Length();
      if ( fabs(n1 * points.Get(i).v) > 1e-3)
	continue;

      bestval = 1e8;
      bestj = 1;
      for (j = 1; j <= points.Size(); j++)
	{
	  Point<3> p2= points.Get(j).p;
	  if (!s2->PointOnSurface (p2))
	    continue;
	  
	  Vec<3> n2;
	  s2->GetNormalVector (p2, n2);
	  n2 /= n2.Length();
	  if ( fabs(n2 * points.Get(j).v) > 1e-3)
	    continue;


	  Vec<3> v(p1, p2);
	  double vl = v.Length();
	  double cl = fabs (v*n1);

	  val = 1 - cl*cl/(vl*vl);

	  val += (points.Get(i).v - points.Get(j).v).Length();

	  if (val < bestval)
	    {
	      bestj = j;
	      bestval = val;
	    }
	}

      (*testout) << "Identify close surfaces special points: pi = " 
		 << points.Get(i).p << ", vi = " << points.Get(i).v 
		 << " pj = " << points.Get(bestj).p 
		 << ", vj = " << points.Get(bestj).v 
		 << " bestval = " << bestval << endl;
    }
}
*/

int CloseEdgesIdentification :: 
Identifyable (const SpecialPoint & sp1, const SpecialPoint & sp2,
	      const TABLE<int> & specpoint2solid,
	      const TABLE<int> & specpoint2surface) const
{
  int i;
  double val;
  
  SpecialPoint hsp1 = sp1;
  SpecialPoint hsp2 = sp2;

  for (i = 1; i <= 1; i++)
    {
      if (!s1->PointOnSurface (hsp1.p))
	continue;

      Vec<3> n1;
      n1 = s1->GetNormalVector (hsp1.p);
      n1 /= n1.Length();
      if ( fabs(n1 * hsp1.v) > 1e-3)
	continue;


      if (!s2->PointOnSurface(hsp2.p))
	continue;

      Vec<3> n2;
      n2 = s2->GetNormalVector (hsp2.p);
      n2 /= n2.Length();
      if ( fabs(n2 * hsp2.v) > 1e-3)
	continue;


      Vec<3> v = hsp2.p - hsp1.p;
      double vl = v.Length();
      double cl = fabs (v*n1);
      

      val = 1 - cl*cl/(vl*vl);
      val += (hsp1.v - hsp2.v).Length();
    
      if (val < 1e-3)
	{
	  return 1;
	}
    }

  return 0;
}




void CloseEdgesIdentification :: IdentifyPoints (Mesh & mesh)
{
  int i, j;
  int i1, i2;

  int np = mesh.GetNP();
  for (i1 = 1; i1 <= np; i1++)
    for (i2 = 1; i2 <= np; i2++)
      {
	if (i2 == i1)
	  continue;
	
	const Point<3> p1 = mesh.Point(i1);
	const Point<3> p2 = mesh.Point(i2);
	Point<3> pp1 = p1;
	Point<3> pp2 = p2;
	
	s1->Project (pp1);
	facet->Project (pp1);
	s2->Project (pp2);
	facet->Project (pp2);

	if (Dist (p1, pp1) > 1e-6 || Dist (p2, pp2) > 1e-6)
	  continue;

	Vec<3> n1, nf, t;
	Vec<3> n = p2 - p1;
	n.Normalize();

	n1 = s1->GetNormalVector (p1);
	nf = facet->GetNormalVector (p1);
	t = Cross (n1, nf);
	t /= t.Length();

	if (fabs (n * t) < 0.5)
	  {
	    (*testout) << "close edges identify points " << p1 << " - " << p2 << endl;
	    mesh.GetIdentifications().Add (i1, i2, nr);
	    mesh.GetIdentifications().SetType(nr,Identifications::CLOSEEDGES);
	  }
      }
}

void CloseEdgesIdentification :: 
BuildSurfaceElements (ARRAY<Segment> & segs,
		      Mesh & mesh, const Surface * surf)
{
  int i1, i2;
  int found = 0;
  int i, j, k;

  if (surf != facet)
    return;

  for (i1 = 1; i1 <= segs.Size(); i1++)
    for (i2 = 1; i2 < i1; i2++)
      {
	const Segment & s1 = segs.Get(i1);
	const Segment & s2 = segs.Get(i2);
	if (mesh.GetIdentifications().Get (s1.p1, s2.p2) &&
	    mesh.GetIdentifications().Get (s1.p2, s2.p1))
	  {
	    Element2d el(QUAD);
	    el.PNum(1) = s1.p1;
	    el.PNum(2) = s1.p2;
	    el.PNum(3) = s2.p2;
	    el.PNum(4) = s2.p1;

	    Vec<3> n = Cross (Point<3> (mesh.Point(el.PNum(2)))-
			      Point<3> (mesh.Point(el.PNum(1))),
			      Point<3> (mesh.Point(el.PNum(3)))-
			      Point<3> (mesh.Point(el.PNum(1))));
	    Vec<3> ns;
	    ns = surf->GetNormalVector (mesh.Point(el.PNum(1)));
	    //(*testout) << "n = " << n << " ns = " << ns << endl;
	    if (n * ns < 0)
	      {
		//(*testout) << "Swap the quad" << endl;
		Swap (el.PNum(1), el.PNum(2));
		Swap (el.PNum(3), el.PNum(4));
	      }
			     
	    
	    Swap (el.PNum(3), el.PNum(4));
	    mesh.AddSurfaceElement (el);
//  	    (*testout) << "add rect element: "
//  		       << mesh.Point (el.PNum(1)) << " - "
//  		       << mesh.Point (el.PNum(2)) << " - "
//  		       << mesh.Point (el.PNum(3)) << " - "
//  		       << mesh.Point (el.PNum(4)) << endl;
	    found = 1;
	  }
      }

  if (found)
    segs.SetSize(0);
}

}

