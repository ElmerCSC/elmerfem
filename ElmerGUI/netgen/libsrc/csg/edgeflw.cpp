#include <mystdlib.h>
#include <meshing.hpp>
#include <csg.hpp>

#undef DEVELOP
// #define DEVELOP

namespace netgen
{

  EdgeCalculation :: 
  EdgeCalculation (const CSGeometry & ageometry,
		   ARRAY<SpecialPoint> & aspecpoints)
    : geometry(ageometry), specpoints(aspecpoints)
  {
    Box<3> bbox = geometry.BoundingBox();

    searchtree = new Point3dTree (bbox.PMin(), bbox.PMax());
    meshpoint_tree = new Point3dTree (bbox.PMin(), bbox.PMax());

    for (int i = 0; i < specpoints.Size(); i++)
      searchtree->Insert (specpoints[i].p, i);

    ideps = 1e-9;
  }

  EdgeCalculation :: ~EdgeCalculation()
  {
    delete searchtree;
    delete meshpoint_tree;
  }


  void EdgeCalculation :: Calc(double h, Mesh & mesh)
  {
    PrintMessage (1, "Find edges");
    PushStatus ("Find edges");


    for (int i = 1; i <= mesh.GetNP(); i++)
      meshpoint_tree->Insert (mesh.Point(i), i);


    // add all special points before edge points (important for periodic identification)
    // JS, Jan 2007
    const double di=1e-7*geometry.MaxSize();
    ARRAY<int> locsearch;

    for (int i = 0; i < specpoints.Size(); i++)
      if (specpoints[i].unconditional)
	{
	  Point<3> p = specpoints[i].p;
          meshpoint_tree -> GetIntersecting (p-Vec<3> (di,di,di),
                                             p+Vec<3> (di,di,di), locsearch);
          
	  if (locsearch.Size() == 0)
            {
              PointIndex pi = mesh.AddPoint (p, specpoints[i].GetLayer(), FIXEDPOINT);
              meshpoint_tree -> Insert (p, pi); 
            }
        }

    /*
      // slow version
    for (int i = 0; i < specpoints.Size(); i++)
      if (specpoints[i].unconditional)
	{
	  Point<3> p = specpoints[i].p;
	  bool found = false;
	  for (int j = 1; j <= mesh.GetNP(); j++)
	    if (Dist (p, mesh.Point(j)) < 1e-8)
	      found = true;
	  if (!found)
	    mesh.AddPoint (p, specpoints[i].GetLayer(), FIXEDPOINT);
	}
    */



    CalcEdges1 (h, mesh);
    SplitEqualOneSegEdges (mesh);
    FindClosedSurfaces (h, mesh);
    PrintMessage (3, cntedge, " edges found");

    PopStatus ();
  }





  void EdgeCalculation :: CalcEdges1 (double h, Mesh & mesh)
  {
    ARRAY<int> hsp(specpoints.Size());
    ARRAY<int> glob2hsp(specpoints.Size());
    ARRAY<int> startpoints, endpoints;


    int pos, ep;
    int layer;

    Point<3> p, np; 
    int pi1, s1, s2, s1_orig, s2_orig;

    ARRAY<Point<3> > edgepoints;
    ARRAY<double> curvelength;
    int copyedge = 0, copyfromedge = -1, copyedgeidentification = -1;

    ARRAY<int> locsurfind, locind;

    int checkedcopy = 0;

    // double size = geometry.MaxSize(); 
    // double epspointdist2 = sqr (size) * 1e-12;
    

    // copy special points to work with
    for (int i = 0; i < specpoints.Size(); i++)
      {
	hsp[i] = i;
	glob2hsp[i] = i;
      }

    //for(int i=0; i<hsp.Size(); i++)
    //  (*testout) << "hsp["<<i<<"] ... " << specpoints[hsp[i]].p << endl;
     

    cntedge = 0;
    INDEX_2_HASHTABLE<int> identification_used(100);  // identification i already used for startpoint j

    mesh.GetIdentifications().Delete();
    
    TABLE<int> specpoint2surface(specpoints.Size());
    if (geometry.identifications.Size())
      {
	for (int i = 0; i < specpoints.Size(); i++)
	  for (int j = 0; j < geometry.GetNSurf(); j++)
	    if (geometry.GetSurface(j)->PointOnSurface (specpoints[i].p))
	      specpoint2surface.Add (i, j);
      }

    TABLE<int> specpoint2tlo(specpoints.Size());
    if (geometry.identifications.Size())
      {
	for (int i = 0; i < specpoints.Size(); i++)
	  for (int j = 0; j < geometry.GetNTopLevelObjects(); j++)
	    {
	      const TopLevelObject * tlo = geometry.GetTopLevelObject (j);
	      if (tlo->GetSolid() && tlo->GetSolid()->VectorIn (specpoints[i].p,specpoints[i].v))
		//if (tlo->GetSolid() && tlo->GetSolid()->IsIn (specpoints[i].p))
		{
#ifdef DEVELOP
		  (*testout) << "point " << specpoints[i].p << " v " <<specpoints[i].v <<" is in " << tlo->GetSolid()->Name() << endl;
#endif
		  specpoint2tlo.Add (i, j);
		}
	    }
      }

    for (int i = 0; i < specpoints.Size(); i++)
      specpoints[i].nr = i;

    while (hsp.Size())
      {
	SetThreadPercent(100 - 100 * double (hsp.Size()) / specpoints.Size());

#ifdef DEVELOP
	(*testout) << "hsp.Size() " << hsp.Size() << " specpoints.Size() " << specpoints.Size() << endl;
	(*testout) << endl << "edge nr " << cntedge+1 << endl;
#endif

	edgepoints.SetSize (0);
	curvelength.SetSize (0);
      

	pi1 = 0;
	copyedge = 0;
	// identifyable point available ?

	
	for (int i = 0; i < geometry.identifications.Size() && !pi1; i++)
	  for (int j = checkedcopy; j < startpoints.Size() && !pi1; j++)
	    {
#ifdef DEVELOP
	      (*testout) << "checking point " << specpoints[startpoints[j]].p
			 << ", v = " << specpoints[startpoints[j]].v 
			 << " for copying (i,j = " << i << ", " << j << ")" << endl;	  
#endif
 	      if (geometry.identifications[i]->IdentifyableCandidate (specpoints[startpoints[j]]) &&
		  geometry.identifications[i]->IdentifyableCandidate (specpoints[endpoints[j]]))
		  
	      
		{
		  int pi1cand = 0;
		  double mindist = 1e10;
		
		  for (int k = 0; k < hsp.Size() && !pi1; k++)
		    {
		      //(*testout) << "   ? identifyable with " << specpoints[hsp[k]].p 
		      //<< ", v = " << specpoints[hsp[k]].v
		      //		 << endl;
		      if (identification_used.Used (INDEX_2(i, startpoints[j])) ||
			  identification_used.Used (INDEX_2(i, hsp[k])))
			{
			  //(*testout) << "failed at pos0" << endl;
			  continue;
			}
		      
		      if (geometry.identifications[i]
			  ->Identifyable(specpoints[startpoints[j]], specpoints[hsp[k]], specpoint2tlo, specpoint2surface) ||
			  geometry.identifications[i]
			  ->Identifyable(specpoints[hsp[k]], specpoints[startpoints[j]], specpoint2tlo, specpoint2surface))
			{
#ifdef DEVELOP
			  (*testout) << "identifyable: " << specpoints[hsp[k]].p << ", v = " << specpoints[hsp[k]].v
				     << " and " << specpoints[startpoints[j]].p << ", v = " << specpoints[startpoints[j]].v 
				     << " (identification " << i+1 << ")" << endl;
#endif

			  if (Dist (specpoints[startpoints[j]].p, specpoints[hsp[k]].p) < mindist)
			    {
			      mindist = Dist (specpoints[startpoints[j]].p, specpoints[hsp[k]].p);
			      pi1cand = k+1;
			    }
			}
		    }
	
	
		  if (pi1cand)
		    {
		      pi1 = pi1cand;
		      copyedge = 1;
		      copyfromedge = j+1;
		      copyedgeidentification = i+1;
		    
		      identification_used.Set (INDEX_2(i, startpoints[j]), 1);
		      identification_used.Set (INDEX_2(i, hsp.Get(pi1)), 1);
		    }
		}
	    }
	
      
	// cannot copy from other ege ?
	if (!pi1)
	  checkedcopy = startpoints.Size();
      
	// unconditional special point available ?
	if (!pi1)
	  for (int i = 1; i <= hsp.Size(); i++)
	    if (specpoints[hsp.Get(i)].unconditional == 1)
	      {
		pi1 = i;
		break;
	      }
 
     
	if (!pi1)
	  {
	    // no unconditional points available, choose first conitional
	    pi1 = 1;	     
	  }

	layer = specpoints[hsp.Get(pi1)].GetLayer();
      

	if (!specpoints[hsp.Get(pi1)].unconditional)
	  {
	    specpoints[hsp.Elem(pi1)].unconditional = 1;
	    for (int i = 1; i <= hsp.Size(); i++)
	      if (i != pi1 && 
		  Dist (specpoints[hsp.Get(pi1)].p, specpoints[hsp.Get(i)].p) < 1e-8*geometry.MaxSize() &&
		  (specpoints[hsp.Get(pi1)].v + specpoints[hsp.Get(i)].v).Length() < 1e-4)
		{
		  // opposite direction
		  specpoints[hsp.Elem(i)].unconditional = 1;
		}
	  }

	cntedge++;
	startpoints.Append (hsp.Get(pi1));

#ifdef DEVELOP
	(*testout) << "start followedge: p1 = " << specpoints[hsp.Get(pi1)].p 
		   << ", v = " << specpoints[hsp.Get(pi1)].v << endl;
#endif

	FollowEdge (pi1, ep, pos, hsp, h, mesh,
		    edgepoints, curvelength);


	if (multithread.terminate)
	  return;
      
	if (!ep)
	  {
	    // ignore starting point
	    hsp.DeleteElement (pi1);
	    cout << "yes, this happens" << endl;
	    continue;
	  }



	endpoints.Append (hsp.Get(ep));


	double elen = 0;
	for (int i = 1; i <= edgepoints.Size()-1; i++)
	  elen += Dist (edgepoints.Get(i), edgepoints.Get(i+1));


	int shortedge = 0;
	for (int i = 1; i <= geometry.identifications.Size(); i++)
	  if (geometry.identifications.Get(i)->ShortEdge(specpoints[hsp.Get(pi1)], specpoints[hsp.Get(ep)]))
	    shortedge = 1;
	// (*testout) << "shortedge = " << shortedge << endl;


	if (!shortedge)
	  {
	    mesh.RestrictLocalHLine (Point3d (specpoints[hsp.Get(pi1)].p), 
				     Point3d (specpoints[hsp.Get(ep)].p), 
				     elen / mparam.segmentsperedge);
	  }
      
	s1 = specpoints[hsp.Get(pi1)].s1;
	s2 = specpoints[hsp.Get(pi1)].s2;
	s1_orig = specpoints[hsp.Get(pi1)].s1_orig;
	s2_orig = specpoints[hsp.Get(pi1)].s2_orig;


	// delete initial, terminal and conditional points

#ifdef DEVELOP
	(*testout) << "terminal point: p = " << specpoints[hsp.Get(ep)].p 
		   << ", v = " << specpoints[hsp.Get(ep)].v << endl;      
#endif

	searchtree -> DeleteElement (hsp.Get(ep));
	searchtree -> DeleteElement (hsp.Get(pi1));

	if (ep > pi1)
	  {
	    glob2hsp[hsp[ep-1]] = -1;
	    glob2hsp[hsp.Last()] = ep-1;
	    hsp.DeleteElement (ep);

	    glob2hsp[hsp[pi1-1]] = -1;
	    glob2hsp[hsp.Last()] = pi1-1;
	    hsp.DeleteElement (pi1);
	  }
	else
	  {
	    glob2hsp[hsp[pi1-1]] = -1;
	    glob2hsp[hsp.Last()] = pi1-1;
	    hsp.DeleteElement (pi1);

	    glob2hsp[hsp[ep-1]] = -1;
	    glob2hsp[hsp.Last()] = ep-1;
	    hsp.DeleteElement (ep);
	  }


	for (int j = 1; j <= edgepoints.Size()-1; j++)
	  {
	    p = edgepoints.Get(j);
	    np = Center (p, edgepoints.Get(j+1));
	    double hd = Dist (p, np);
 

	    Box<3> boxp (np - (1.2 * hd) * Vec<3> (1, 1, 1),
			 np + (1.2 * hd) * Vec<3> (1, 1, 1));
	    searchtree -> GetIntersecting (boxp.PMin(), boxp.PMax(), locind);	    

	    for (int i = 0; i < locind.Size(); i++)
	      {
		if ( specpoints[locind[i]].HasSurfaces (s1, s2) &&
		     specpoints[locind[i]].unconditional == 0)
		  {
		    searchtree -> DeleteElement (locind[i]);

		    int li = glob2hsp[locind[i]];
		    glob2hsp[locind[i]] = -1;
		    glob2hsp[hsp.Last()] = li;
		    hsp.Delete (li);
		  }
	      }


	    /*
	    for (int i = 1; i <= hsp.Size(); i++)
	      if ( specpoints[hsp.Get(i)].HasSurfaces (s1, s2) &&
		   specpoints[hsp.Get(i)].unconditional == 0 &&
		   Dist2 (np, specpoints[hsp.Get(i)].p) < 1.2 * hd)
		{
		  searchtree -> DeleteElement (hsp.Get(i)+1);
		  hsp.DeleteElement (i);
		  i--;
		}
	    */
	  }

      
	ARRAY<Segment> refedges;
	ARRAY<bool> refedgesinv;
      

	AnalyzeEdge (s1_orig, s2_orig, s1, s2, pos, layer,
		     edgepoints,
		     refedges, refedgesinv);


	for (int i = 0; i < refedges.Size(); i++)
	  refedges[i].edgenr = cntedge;

	
#ifdef DEVELOP
	(*testout) << "edge " << cntedge << endl
		   << "startp: " << specpoints[startpoints.Last()].p 
		   << ", v = " << specpoints[startpoints.Last()].v << endl
		   << "copy = " << copyedge << endl
		   << refedges.Size() << " refedges: ";
	for (int i = 1; i <= refedges.Size(); i++)
	  (*testout) << " " << refedges.Get(i).si;
	(*testout) << endl;
	if (refedgesinv.Size())
	  (*testout) << "inv[1] = " << refedgesinv.Get(1) << endl;
#endif

	if (refedges.Size() == 0)
	  throw NgException ("Problem in edge detection");

      
	if (!copyedge)
	  {
	    // (*testout) << "store edge" << endl;
	    // int oldnseg = mesh.GetNSeg();

	    if (!shortedge)
	      StoreEdge (refedges, refedgesinv, 
			 edgepoints, curvelength, layer, mesh);
	    else
	      StoreShortEdge (refedges, refedgesinv, 
			      edgepoints, curvelength, layer, mesh);

	    for(int i = 0; i < refedges.Size(); i++)
	      {
		refedges[i].surfnr1 = geometry.GetSurfaceClassRepresentant(refedges[i].surfnr1);
		refedges[i].surfnr2 = geometry.GetSurfaceClassRepresentant(refedges[i].surfnr2);
	      }


	    /*
	      for (int i = oldnseg+1; i <= mesh.GetNSeg(); i++)
	      for (int j = 1; j <= oldnseg; j++)
	      {
	      const Point<3> & l1p1 = mesh.Point (mesh.LineSegment(i).p1);
	      const Point<3> & l1p2 = mesh.Point (mesh.LineSegment(i).p2);
	      const Point<3> & l2p1 = mesh.Point (mesh.LineSegment(j).p1);
	      const Point<3> & l2p2 = mesh.Point (mesh.LineSegment(j).p2);
	      Vec<3> vl1(l1p1, l1p2);
	      for (double lamk = 0; lamk <= 1; lamk += 0.1)
	      {
	      Point<3> l2p = l1p1 + lamk * vl1;
	      double dist = sqrt (MinDistLP2 (l2p1, l2p2, l2p));
	      if (dist > 1e-12)
	      mesh.RestrictLocalH (l2p, 3*dist);
	      }
	      }
	    */
	  }
	else
	  {
	    CopyEdge (refedges, refedgesinv,
		      copyfromedge, 
		      specpoints[startpoints.Get(copyfromedge)].p,
		      specpoints[endpoints.Get(copyfromedge)].p,
		      edgepoints.Get(1), edgepoints.Last(),
		      copyedgeidentification, 
		      layer,
		      mesh);
	  }
	
	
// 	for(int i=0; i<hsp.Size(); i++)
// 	  {
// 	    (*testout) << "pos2 hsp["<<i<<"] ... " << specpoints[hsp[i]].p << endl;
// 	  }
      }
  }
  





  /*
    If two or more edges share the same initial and end-points,
    then they need at least two segments 
  */
  void EdgeCalculation ::
  SplitEqualOneSegEdges (Mesh & mesh)
    {
    //    int i, j;
    SegmentIndex si;
    PointIndex pi;

    ARRAY<int> osedges(cntedge);
    INDEX_2_HASHTABLE<int> osedgesht (cntedge+1);

    osedges = 2;

    // count segments on edges
    for (si = 0; si < mesh.GetNSeg(); si++)
      {
	const Segment & seg = mesh[si];
	if (seg.seginfo && seg.edgenr >= 1 && seg.edgenr <= cntedge)
	  osedges.Elem(seg.edgenr)--;
      }

    // flag one segment edges
    for (int i = 0; i < cntedge; i++)
      osedges[i] = (osedges[i] > 0) ? 1 : 0;

    for (si = 0; si < mesh.GetNSeg(); si++)
      {
	const Segment & seg = mesh[si];
	if (seg.seginfo && seg.edgenr >= 1 && seg.edgenr <= cntedge)
	  {
	    if (osedges.Get(seg.edgenr))
	      {
		INDEX_2 i2(seg.p1, seg.p2);
		i2.Sort ();
		if (osedgesht.Used (i2))
		  osedgesht.Set (i2, 2);
		else
		  osedgesht.Set (i2, 1);
	      }
	  }
      }


    // one edge 1 segment, other 2 segments 
    // yes, it happens !
    point_on_edge_problem = 0;
    for (int i = 1; i <= osedgesht.GetNBags(); i++)
      for (int j = 1; j <= osedgesht.GetBagSize(i); j++)
	{
	  INDEX_2 i2; 
	  int val;
	  osedgesht.GetData (i, j, i2, val);

	  const Point<3> & p1 = mesh[PointIndex(i2.I1())];
	  const Point<3> & p2 = mesh[PointIndex(i2.I2())];
	  Vec<3> v = p2 - p1;
	  double vlen = v.Length();
	  v /= vlen;
	  for (pi = PointIndex::BASE; 
	       pi < mesh.GetNP()+PointIndex::BASE; pi++)

	    if (pi != i2.I1() && pi != i2.I2())
	      {
		const Point<3> & p = mesh[pi];
		Vec<3> v2 = p - p1;
		double lam = (v2 * v);
		if (lam > 0 && lam < vlen)
		  {
		    Point<3> hp = p1 + lam * v;
		    if (Dist (p, hp) < 1e-4 * vlen)
		      {
			PrintWarning ("Point on edge !!!");
			cout << "seg: " << i2 << ", p = " << pi << endl;
			osedgesht.Set (i2, 2);		      
			point_on_edge_problem = 1;

			(*testout) << "Point on edge" << endl
				   << "seg = " << i2 << ", p = " << pi << endl
				   << "pos = " << p << ", projected = " << hp << endl
				   << "seg is = " << mesh.Point(i2.I1()) << " - " << mesh.Point(i2.I2()) << endl;
		      }
		  }
	      }
	}


    // insert new points
    osedges = -1;

    int nseg = mesh.GetNSeg();
    for (si = 0; si < nseg; si++)
      {
	const Segment & seg = mesh[si];
	if (seg.seginfo && seg.edgenr >= 1 && seg.edgenr <= cntedge)
	  {
	    INDEX_2 i2(seg.p1, seg.p2);
	    i2.Sort ();
	    if (osedgesht.Used (i2) &&
		osedgesht.Get (i2) == 2 &&
		osedges.Elem(seg.edgenr) == -1)
	      {
		Point<3> newp = Center (mesh[PointIndex(seg.p1)],
					mesh[PointIndex(seg.p2)]);

		ProjectToEdge (geometry.GetSurface(seg.surfnr1), 
			       geometry.GetSurface(seg.surfnr2), 
			       newp);

		osedges.Elem(seg.edgenr) = 
		  mesh.AddPoint (newp, mesh[PointIndex(seg.p1)].GetLayer(), EDGEPOINT);
		meshpoint_tree -> Insert (newp, osedges.Elem(seg.edgenr));
	      }
	  }
      }


    for (int i = 1; i <= nseg; i++)
      {
	Segment & seg = mesh.LineSegment (i);
	if (seg.edgenr >= 1 && seg.edgenr <= cntedge)
	  {
	    if (osedges.Get(seg.edgenr) != -1)
	      {
		Segment newseg = seg;
		newseg.p1 = osedges.Get(seg.edgenr);
		seg.p2 = osedges.Get(seg.edgenr);
		mesh.AddSegment (newseg);
	      }
	  }
      }

  }



  void EdgeCalculation :: 
  FollowEdge (int pi1, int & ep, int & pos,
	      const ARRAY<int> & hsp,
	      double h, const Mesh & mesh,
	      ARRAY<Point<3> > & edgepoints,
	      ARRAY<double> & curvelength)
  {
    int s1, s2, s1_rep, s2_rep;
    double len, steplen, cursteplen, loch;
    Point<3> p, np, pnp;
    Vec<3> a1, a2, t;

    ARRAY<int> locind;

    double size = geometry.MaxSize();  
    double epspointdist2 = size * 1e-6;
    epspointdist2 = sqr (epspointdist2);
    int uselocalh = mparam.uselocalh;


    s1_rep = specpoints[hsp.Get(pi1)].s1;
    s2_rep = specpoints[hsp.Get(pi1)].s2;
    s1 = specpoints[hsp.Get(pi1)].s1_orig;
    s2 = specpoints[hsp.Get(pi1)].s2_orig;
  
    p = specpoints[hsp.Get(pi1)].p;
    //ProjectToEdge (geometry.GetSurface(s1), 
    //               geometry.GetSurface(s2), p);
    geometry.GetSurface(s1) -> CalcGradient (p, a1);
    geometry.GetSurface(s2) -> CalcGradient (p, a2);

    t = Cross (a1, a2);
    t.Normalize();

    pos = (specpoints[hsp.Get(pi1)].v * t) > 0;
    if (!pos) t *= -1;

  
    edgepoints.Append (p);
    curvelength.Append (0);
    len = 0;

    // (*testout) << "geometry.GetSurface(s1) -> LocH (p, 3, 1, h) " << geometry.GetSurface(s1) -> LocH (p, 3, 1, h)
    // << " geometry.GetSurface(s2) -> LocH (p, 3, 1, h) " << geometry.GetSurface(s2) -> LocH (p, 3, 1, h) << endl;

    loch = min2 (geometry.GetSurface(s1) -> LocH (p, 3, 1, h), 
		 geometry.GetSurface(s2) -> LocH (p, 3, 1, h));
  
  
  
    if (uselocalh)
      {
	double lh = mesh.GetH(p);
	// (*testout) << "lh " << lh << endl;
	if (lh < loch)
	  loch = lh;
      }

    steplen = 0.1 * loch;
  
    do
      {
	if (multithread.terminate)
	  return;
      
	if (fabs (p(0)) + fabs (p(1)) + fabs (p(2)) > 100000*size)
	  {
	    ep = 0;
	    PrintWarning ("Give up line");
	    break;
	  }

	if (steplen > 0.1 * loch) steplen = 0.1 * loch;
      
	steplen *= 2;
	do
	  {
	    steplen *= 0.5;
	    np = p + steplen * t;
	    pnp = np;
	    ProjectToEdge (geometry.GetSurface(s1), 
			   geometry.GetSurface(s2), pnp);
	  }
	while (Dist (np, pnp) > 0.1 * steplen);

      
	cursteplen = steplen;
	if (Dist (np, pnp) < 0.01 * steplen) steplen *= 2;
      
 
	np = pnp;
	ep = 0;
      
	double hvtmin = 1.5 * cursteplen;
      
	Box<3> boxp (p - (2 * cursteplen) * Vec<3> (1, 1, 1),
		     p + (2 * cursteplen) * Vec<3> (1, 1, 1));

	searchtree -> GetIntersecting (boxp.PMin(), boxp.PMax(), locind);
	
	for (int i = 0; i < locind.Size(); i++)
	  {
	    Vec<3> hv = specpoints[locind[i]].p - p;
	    if (hv.Length2() > 9 * cursteplen * cursteplen)
	      continue;

	    double hvt = hv * t;
	    hv -= hvt * t;

	    if (hv.Length() < 0.2 * cursteplen &&
		hvt > 0 && 
		//		  hvt < 1.5 * cursteplen &&
		hvt < hvtmin && 
		specpoints[locind[i]].unconditional == 1 &&
		(specpoints[locind[i]].v + t).Length() < 0.4  ) 
	      {
		Point<3> hep = specpoints[locind[i]].p;
		ProjectToEdge (geometry.GetSurface(s1), 
			       geometry.GetSurface(s2), hep);            
	      
	      
		if (Dist2 (hep, specpoints[locind[i]].p) < epspointdist2 )
		  {
		    geometry.GetSurface(s1) -> CalcGradient (hep, a1);
		    geometry.GetSurface(s2) -> CalcGradient (hep, a2);
		    Vec<3> ept = Cross (a1, a2);
		    ept /= ept.Length();
		    if (!pos) ept *= -1;
		  
		    if ( (specpoints[locind[i]].v + ept).Length() < 1e-4 )
		      {
			np = specpoints[locind[i]].p;

			for (int jj = 0; jj < hsp.Size(); jj++)
			  if (hsp[jj] == locind[i])
			    ep = jj+1;
			    
			if (!ep) 
			  cerr << "endpoint not found" << endl;
			  //			ep = i;
			hvtmin = hvt;
			//			  break;
		      }
		  }
	      }
	  }




	/*
	for (int i = 1; i <= hsp.Size(); i++)
	  {
	    if (!boxp.IsIn (specpoints[hsp.Get(i)].p))
	      continue;
	  
	    Vec<3> hv = specpoints[hsp.Get(i)].p - p;
	    if (hv.Length2() > 9 * cursteplen * cursteplen)
	      continue;

	    double hvt = hv * t;
	    hv -= hvt * t;
	  
	    if (hv.Length() < 0.2 * cursteplen &&
		hvt > 0 && 
		//		  hvt < 1.5 * cursteplen &&
		hvt < hvtmin && 
		specpoints[hsp.Get(i)].unconditional == 1 &&
		(specpoints[hsp.Get(i)].v + t).Length() < 0.4  ) 
	      {
		Point<3> hep = specpoints[hsp.Get(i)].p;
		ProjectToEdge (geometry.GetSurface(s1), 
			       geometry.GetSurface(s2), hep);            
	      
	      
		if (Dist2 (hep, specpoints[hsp.Get(i)].p) < epspointdist2 )
		  {
		    geometry.GetSurface(s1) -> CalcGradient (hep, a1);
		    geometry.GetSurface(s2) -> CalcGradient (hep, a2);
		    Vec<3> ept = Cross (a1, a2);
		    ept /= ept.Length();
		    if (!pos) ept *= -1;
		  
		    if ( (specpoints[hsp.Get(i)].v + ept).Length() < 1e-4 )
		      {
			np = specpoints[hsp.Get(i)].p;
			ep = i;
			hvtmin = hvt;
			//			  break;
		      }
		  }
	      }
	  }
	*/

	loch = min2 (geometry.GetSurface(s1_rep) -> LocH (np, 3, 1, h), 
		     geometry.GetSurface(s2_rep) -> LocH (np, 3, 1, h));
        loch = max2 (loch, mparam.minh);

	if (uselocalh)
	  {
	    double lh = mesh.GetH(np);
	    if (lh < loch)
	      loch = lh;
	  }
      
      
	len += Dist (p, np) / loch;
	edgepoints.Append (np);
	curvelength.Append (len);
      
	p = np;
      
	geometry.GetSurface(s1) -> CalcGradient (p, a1);
	geometry.GetSurface(s2) -> CalcGradient (p, a2);
	t = Cross (a1, a2);
	t.Normalize();
	if (!pos) t *= -1;
      }
    while (! ep);
  }







  void EdgeCalculation :: 
  AnalyzeEdge (int s1, int s2, int s1_rep, int s2_rep, int pos, int layer,
	       const ARRAY<Point<3> > & edgepoints,
	       ARRAY<Segment> & refedges,
	       ARRAY<bool> & refedgesinv)
  {
    int j, k, l;
    int hi;
    Point<3> hp;
    Vec<3> t, a1, a2, m, n;
    Segment seg;
    Solid * locsol;
    ARRAY<int> locsurfind, locsurfind2;

    ARRAY<int> edges_priority;

    double size = geometry.MaxSize();
    bool debug = 0;

#ifdef DEVELOP
    debug = 1;
#endif
    
    if (debug)
      {
	(*testout) << "debug edge !!!" << endl;
	(*testout) << "edgepoints = " << edgepoints << endl;
	(*testout) << "s1, s2 = " << s1_rep << " - " << s2_rep << endl;
      }

    refedges.SetSize(0);
    refedgesinv.SetSize(0);
    hp = Center (edgepoints[0], edgepoints[1]);
    ProjectToEdge (geometry.GetSurface(s1), geometry.GetSurface(s2), hp);
    
    if (debug)
      *testout << "hp = " << hp << endl;

    geometry.GetSurface(s1) -> CalcGradient (hp, a1);
    geometry.GetSurface(s2) -> CalcGradient (hp, a2);
    t = Cross (a1, a2);
    t.Normalize();
    if (!pos) t *= -1;    

  
    for (int i = 0; i < geometry.GetNTopLevelObjects(); i++)
      {
	if (geometry.GetTopLevelObject(i)->GetLayer() != layer) 
	  continue;
      
	const Solid * sol = geometry.GetTopLevelObject(i)->GetSolid();
	const Surface * surf = geometry.GetTopLevelObject(i)->GetSurface();

	sol -> TangentialSolid (hp, locsol, locsurfind, size*ideps);

	//*testout << "hp = " << hp << endl;
	//(*testout) << "locsol: " << endl;
	//if (locsol) locsol->Print(*testout);
	//(*testout) << endl;


	if (!locsol) continue;

	BoxSphere<3> boxp (hp, hp);
	boxp.Increase (1e-8*size);
	boxp.CalcDiamCenter();
      
	ReducePrimitiveIterator rpi(boxp);
	UnReducePrimitiveIterator urpi;
      
	((Solid*)locsol) -> IterateSolid (rpi);

	locsol -> CalcSurfaceInverse ();
      
	if (!surf)
	  {
	    locsol -> GetTangentialSurfaceIndices (hp,locsurfind,ideps*size);
	  }
	else
	  {
	    /*
	      if (fabs (surf->CalcFunctionValue (hp)) < 1e-6)
	      continue;
	    */
	    locsurfind.SetSize(1);
	    locsurfind[0] = -1;
	    for (j = 0; j < geometry.GetNSurf(); j++)
	      if (geometry.GetSurface(j) == surf)
		{
		  locsurfind[0] = j;
		  //		      geometry.GetSurfaceClassRepresentant(j);
		  break;
		}
	  }

	((Solid*)locsol) -> IterateSolid (urpi);

      
	if (debug)
	  (*testout) << "edge of tlo " << i << ", has " << locsurfind.Size() << " faces." << endl;
      

	for (j = locsurfind.Size()-1; j >= 0; j--)
	  if (fabs (geometry.GetSurface(locsurfind[j])
		    ->CalcFunctionValue (hp) ) > ideps*size)
	    locsurfind.Delete(j);
      
	if (debug)
	  (*testout) << locsurfind.Size() << " faces on hp" << endl;



	for (j = 0; j < locsurfind.Size(); j++)
	  {      
	    int lsi = locsurfind[j];
	    int rlsi = geometry.GetSurfaceClassRepresentant(lsi);
	  
	    Vec<3> rn;

	    // n is outer normal to solid
	    n = geometry.GetSurface(lsi) -> GetNormalVector (hp);
            if (debug)
              *testout << "n1 = " << n << endl;
	    if (geometry.GetSurface (lsi)->Inverse())
	      n *= -1;
	  
	    if (fabs (t * n) > 1e-4) continue;
	    if (debug)
	      {
		(*testout) << "face " << locsurfind[j] << ", rep = " << rlsi 
			   << " has (t*n) = " << (t*n) << endl;
		(*testout) << "n = " << n << endl;
	      }
	  
	    // rn is normal to class representant
	    rn = geometry.GetSurface(rlsi) -> GetNormalVector (hp);
	    if (debug)
	      {
		(*testout) << "rn = " << rn << endl;
	      }
	    
	    //if( n*rn < 0)
	    // rn *= -1;

	    bool sameasref = ((n * rn) > 0);

	    //m = Cross (t, rn);
	    m = Cross (t, n); if(!sameasref) m*=-1.;
	    
	    m.Normalize();

	    
	    if (debug)
	      (*testout) << "m = " << m << endl;


	    //bool founddirection = false;
	    //int k;
	    double eps = 1e-8*size;

	    ARRAY<bool> pre_ok(2);

 	    do
 	      {
 		eps *= 0.5;
 		pre_ok[0] = (locsol -> VectorIn2 (hp, m, n, eps) == IS_OUTSIDE &&
			     locsol -> VectorIn2 (hp, m, -1. * n, eps) == IS_INSIDE);
 		pre_ok[1] = (locsol -> VectorIn2 (hp, -1.*m, n, eps) == IS_OUTSIDE &&
			     locsol -> VectorIn2 (hp, -1.*m, -1. * n, eps) == IS_INSIDE);
 	      }
 	    while(pre_ok[0] && pre_ok[1] && eps > 1e-16*size);

            if (debug)
              {
                *testout << "eps = " << eps << ", size = " << size << endl;
                *testout << "pre_ok[0,1] = " << pre_ok[0] << "," << pre_ok[1] << endl;
              }

	    eps = 1e-8*size;
	    

	    for (k = 1; k <= 2; k ++)
	      {
		bool edgeinv = (k == 2);
	      
		if (debug)
		  {
		    (*testout) << "onface(" << hp << ", " << m << ")= " << flush;
		    (*testout) << locsol->OnFace (hp, m, eps) << flush;
		    (*testout) << " n " << n << flush;
		    (*testout) << " vec2in = "
			       << locsol -> VectorIn2 (hp, m, n, eps) << " and " 
			       << locsol -> VectorIn2 (hp, m, -1 * n, eps) << endl;
		  }

		//	      if (locsol -> OnFace (hp, m))
		

		// one side must be inside, the other must be outside
		bool ok = (pre_ok[k-1] || 
			   (locsol -> VectorIn2 (hp, m, n, eps) == IS_OUTSIDE &&
			    locsol -> VectorIn2 (hp, m, -1 * n, eps) == IS_INSIDE));

		if (debug)
		  (*testout) << "ok (before) " << ok <<  endl;

		// compute second order approximation
		// curve = hp + t m + t*t/2 m2
		Vec<3> grad, m2;
		Mat<3> hesse;
		geometry.GetSurface(lsi) -> CalcGradient (hp, grad);
		geometry.GetSurface(lsi) -> CalcHesse (hp, hesse);
		double fac = -(m * (hesse * m)) / (grad * grad);
		m2 = fac * grad;
		// (*testout) << "hp = " << hp << ", m = " << m << ", m2 = " << m2 << endl;

		Solid * locsol2;
		locsol -> TangentialSolid3 (hp, m, m2, locsol2, locsurfind2, ideps*size);
		if (!locsol2) ok = 0;
		delete locsol2;


		if (ok)
		  {
		    if (debug)
		      (*testout) << "is true" << endl;
		    hi = 0;
		    for (l = 1; !hi && l <= refedges.Size(); l++)
		      {
			   if (refedges.Get(l).si == rlsi &&     // JS sept 2006
			       // if (refedges.Get(l).si == lsi &&
			       refedgesinv.Get(l) == edgeinv)
			     {
			       hi = l;
			     }
		      }
		  
		    if (!hi)
		      {
			 seg.si = rlsi;  // JS Sept 2006
			 // seg.si = lsi;
			seg.domin = -1;
			seg.domout = -1;
			seg.tlosurf = -1;
			//seg.surfnr1 = s1_rep;
			//seg.surfnr2 = s2_rep;
			seg.surfnr1 = s1;
			seg.surfnr2 = s2;
			hi = refedges.Append (seg);
			refedgesinv.Append (edgeinv);
			edges_priority.Append((pre_ok[k-1]) ? 1 : 0);
		      }
		    else
		      {
			if(edges_priority[hi-1] / 10 == -i-1)
			  edges_priority[hi-1] = 10*(i+1);
			else
			  edges_priority[hi-1] = -10*(i+1);
		      }
		  
		    if (!surf)
		      {
			if (sameasref)
			  refedges.Elem(hi).domin = i;
			else 
			  refedges.Elem(hi).domout = i;
		      }
		    else
		      refedges.Elem(hi).tlosurf = i;

		    if(pre_ok[k-1])
		      edges_priority[hi-1] = 1;
		    

		    if (debug)
		      (*testout) << "add ref seg:" 
				 << "si = " << refedges.Get(hi).si
				 << ", domin = " << refedges.Get(hi).domin
				 << ", domout = " << refedges.Get(hi).domout
				 << ", surfnr1/2 = " << refedges.Get(hi).surfnr1
				 << ", " << refedges.Get(hi).surfnr2
				 << ", inv = " << refedgesinv.Get(hi) 
				 << ", refedgenr = " << hi
				 << ", priority = " << edges_priority[hi-1]
				 << ", hi = " << hi 
				 << endl;
		  }
		else
		  {
		    if (debug)
		      (*testout) << "is false" << endl;
		  }
		m *= -1;
	      } 
	  }
	delete locsol;          
      }

   
    if (debug)
      {
	*testout << "Refsegments, before delete: " << endl << refedges << endl;
	*testout << "inv: " << endl << refedgesinv << endl;
      }
    
    BitArray todelete(refedges.Size());
    todelete.Clear();


    for(int i=0; i<refedges.Size()-1; i++)
      {
	for(j=i+1; !todelete.Test(i) && j<refedges.Size(); j++)
	  {
	    if(todelete.Test(j))
	      continue;

	    if(refedges[i].si == refedges[j].si &&
	       refedges[i].domin == refedges[j].domin &&
	       refedges[i].domout == refedges[j].domout &&
	       geometry.GetSurfaceClassRepresentant(refedges[i].surfnr1) == geometry.GetSurfaceClassRepresentant(refedges[j].surfnr1) &&
	       geometry.GetSurfaceClassRepresentant(refedges[i].surfnr2) == geometry.GetSurfaceClassRepresentant(refedges[j].surfnr2)
	       // && refedgesinv[i] == refedgesinv[j] // JS, 20060802
	       )
	      {
		if(debug)
		  (*testout) << "equal segments: " << refedges[i] << " pri " << edges_priority[i] 
			     << " tlosurf " << refedges[i].tlosurf
			     << "\n and " << refedges[j] << " pri " << edges_priority[j]
			     << " tlosurf " << refedges[i].tlosurf << endl;
		
		if(edges_priority[i] < 10 && edges_priority[i] < edges_priority[j])
		  {
		    todelete.Set(i);
		  }
		else if (edges_priority[j] < 10 && edges_priority[i] > edges_priority[j])
		  {
		    todelete.Set(j);
		  }
	      }
	  }
	  
      }
    
    int num = refedges.Size();

    for(int i=refedges.Size()-1; num>2 && i>=0; i--)
      if(todelete.Test(i))
	{
	  refedges.Delete(i);
	  refedgesinv.Delete(i);
	  num--;
	}

    
    if (debug)
      {
	*testout << "Refsegments: " << endl << refedges << endl;
      }
  }



  void EdgeCalculation :: 
  StoreEdge (const ARRAY<Segment> & refedges,
	     const ARRAY<bool> & refedgesinv,
	     const ARRAY<Point<3> > & edgepoints,
	     const ARRAY<double> & curvelength,
	     int layer,
	     Mesh & mesh)
  {
  
    // Calculate optimal element-length
    int i, j, k;
    PointIndex pi;
    int ne;

    double len, corr, lam;
    PointIndex thispi, lastpi;
    Point<3> p, np;
    Segment seg;

    const Surface * surf1 = geometry.GetSurface (refedges.Get(1).surfnr1);
    const Surface * surf2 = geometry.GetSurface (refedges.Get(1).surfnr2);

    (*testout) << "s1 " << refedges.Get(1).surfnr1 << " s2 " << refedges.Get(1).surfnr2
	       << " rs1 " << geometry.GetSurfaceClassRepresentant(refedges.Get(1).surfnr1)
	       << " rs2 " << geometry.GetSurfaceClassRepresentant(refedges.Get(1).surfnr2) << endl;

    len = curvelength.Last();
    ne = int (len + 0.5);
    if (ne == 0) ne = 1;
    if (Dist (edgepoints.Get(1), edgepoints.Last()) < 1e-8*geometry.MaxSize() && 
	ne <= 6) 
      ne = 6;
    corr = len / ne;

    // generate initial point
    p = edgepoints.Get(1);
    lastpi = -1;

    /*
    for (pi = PointIndex::BASE; 
	 pi < mesh.GetNP()+PointIndex::BASE; pi++)
      if (Dist (mesh[pi], p) < 1e-6)
	{
	  lastpi = pi;
	  break;
	}
    */

    const double di=1e-7*geometry.MaxSize();

    ARRAY<int> locsearch;
    meshpoint_tree -> GetIntersecting (p-Vec<3> (di,di,di),
				       p+Vec<3> (di,di,di), locsearch);
    if (locsearch.Size())
      lastpi = locsearch[0];
				       


    if (lastpi == -1)
      {
	lastpi = mesh.AddPoint (p, layer, FIXEDPOINT);
	meshpoint_tree -> Insert (p, lastpi); 
	// (*testout) << "test1, store point " << lastpi << ", p = " << p << endl;
      }
  
    j = 1;
    for (i = 1; i <= ne; i++)
      {
	while (curvelength.Get(j) < i * corr && j < curvelength.Size()) j++;


	lam = (i * corr - curvelength.Get(j-1)) / 
	  (curvelength.Get(j) - curvelength.Get(j-1));

	np(0) = (1-lam) * edgepoints.Get(j-1)(0) + lam * edgepoints.Get(j)(0);
	np(1) = (1-lam) * edgepoints.Get(j-1)(1) + lam * edgepoints.Get(j)(1);
	np(2) = (1-lam) * edgepoints.Get(j-1)(2) + lam * edgepoints.Get(j)(2);
      
	thispi = -1;
	if (i == ne)
	  {
	    /*
	  for (pi = PointIndex::BASE; 
	       pi < mesh.GetNP()+PointIndex::BASE; pi++)
	    if (Dist(mesh[pi], np) < 1e-6)
	      thispi = pi;
	    */
	    
	    meshpoint_tree -> GetIntersecting (np-Vec<3> (di,di,di),
					       np+Vec<3> (di,di,di), locsearch);
	    if (locsearch.Size())
	      thispi = locsearch[0];
	  }

	if (thispi == -1)
	  {
	    ProjectToEdge (surf1, surf2, np);
	    thispi = mesh.AddPoint (np, layer, (i==ne) ? FIXEDPOINT : EDGEPOINT);
	   
	    meshpoint_tree -> Insert (np, thispi);
	    // (*testout) << "test2, store point " << thispi << ", p = " << np << endl;
	  }

	for (k = 1; k <= refedges.Size(); k++)
	  {
	    if (refedgesinv.Get(k))
	      {
		seg.p1 = lastpi;
		seg.p2 = thispi;
	      }
	    else
	      {
		seg.p1 = thispi;
		seg.p2 = lastpi;
	      }
	    seg.si = refedges.Get(k).si;
	    seg.domin = refedges.Get(k).domin;
	    seg.domout = refedges.Get(k).domout;
	    seg.tlosurf = refedges.Get(k).tlosurf;
	    seg.edgenr = refedges.Get(k).edgenr;
	    seg.surfnr1 = refedges.Get(k).surfnr1;
	    seg.surfnr2 = refedges.Get(k).surfnr2;
	    seg.seginfo = 0;
	    if (k == 1) seg.seginfo = (refedgesinv.Get(k)) ? 2 : 1;
	    mesh.AddSegment (seg);
	    //(*testout) << "add seg " << mesh[seg.p1] << "-" << mesh[seg.p2] << endl;
	    //(*testout) << "refedge " << k << " surf1 " << seg.surfnr1 << " surf2 " << seg.surfnr2 << " inv " << refedgesinv.Get(k) << endl;
	  
	    double maxh = min2 (geometry.GetSurface(seg.surfnr1)->GetMaxH(),
				geometry.GetSurface(seg.surfnr2)->GetMaxH());
			      
	    if (seg.domin != -1)
	      {
		const Solid * s1 = 
		  geometry.GetTopLevelObject(seg.domin) -> GetSolid();
		maxh = min2 (maxh, s1->GetMaxH());
		maxh = min2 (maxh, geometry.GetTopLevelObject(seg.domin)->GetMaxH());
		mesh.RestrictLocalH (p, maxh);
		mesh.RestrictLocalH (np, maxh);
	      }
	    if (seg.domout != -1)
	      {
		const Solid * s1 = 
		  geometry.GetTopLevelObject(seg.domout) -> GetSolid();
		maxh = min2 (maxh, s1->GetMaxH());
		maxh = min2 (maxh, geometry.GetTopLevelObject(seg.domout)->GetMaxH());
		mesh.RestrictLocalH (p, maxh);
		mesh.RestrictLocalH (np, maxh);
	      }
	    if (seg.tlosurf != -1)
	      {
		double hi = geometry.GetTopLevelObject(seg.tlosurf) -> GetMaxH();
		maxh = min2 (maxh, hi);
		mesh.RestrictLocalH (p, maxh);
		mesh.RestrictLocalH (np, maxh);
	      }	  
	  }
      
	p = np;
	lastpi = thispi;
      }

#ifdef DEVELOP
    (*testout) << " eplast = " << lastpi << " = " << p << endl;
#endif
  }
  





  void EdgeCalculation :: 
  StoreShortEdge (const ARRAY<Segment> & refedges,
		  const ARRAY<bool> & refedgesinv,
		  const ARRAY<Point<3> > & edgepoints,
		  const ARRAY<double> & curvelength,
		  int layer,
		  Mesh & mesh)
  {
  
    // Calculate optimal element-length
    PointIndex pi;
    // int ne;
    Segment seg;

    /*
      double len, corr, lam;
      int thispi, lastpi;
      Point<3> p, np;


      const Surface * surf1 = geometry.GetSurface (refedges.Get(1).surfnr1);
      const Surface * surf2 = geometry.GetSurface (refedges.Get(1).surfnr2);

      len = curvelength.Last();
      ne = int (len + 0.5);
      if (ne == 0) ne = 1;
      if (Dist2 (edgepoints[1], edgepoints.Last()) < 1e-8 && 
      ne <= 6) 
      ne = 6;
      corr = len / ne;
    */

    // generate initial point
    Point<3> p = edgepoints[0];
    PointIndex pi1 = -1;
    for (pi = PointIndex::BASE; 
	 pi < mesh.GetNP()+PointIndex::BASE; pi++)

      if (Dist (mesh[pi], p) < 1e-6*geometry.MaxSize())
	{
	  pi1 = pi;
	  break;
	}

    if (pi1 == -1) 
      {
	pi1 = mesh.AddPoint (p, layer, FIXEDPOINT);
	meshpoint_tree -> Insert (p, pi1);
	// (*testout) << "test3, store point " << pi1 << ", p = " << p << endl;
      }

    p = edgepoints.Last();
    PointIndex pi2 = -1;
    for (pi = PointIndex::BASE; 
	 pi < mesh.GetNP()+PointIndex::BASE; pi++)

      if (Dist (mesh[pi], p) < 1e-6*geometry.MaxSize())
	{
	  pi2 = pi;
	  break;
	}
    if (pi2==-1) 
      {
	pi2 = mesh.AddPoint (p, layer, FIXEDPOINT);
	meshpoint_tree -> Insert (p, pi2);
	// (*testout) << "test4, store point " << pi2 << ", p = " << p << endl;
      }

    /*
  
    j = 1;
    for (i = 1; i <= ne; i++)
    {
    while (curvelength[j] < i * corr && j < curvelength.Size()) j++;
      
    lam = (i * corr - curvelength[j-1]) / 
    (curvelength[j] - curvelength[j-1]);
      
    np(0) = (1-lam) * edgepoints[j-1](0) + lam * edgepoints[j](0);
    np(1) = (1-lam) * edgepoints[j-1](1) + lam * edgepoints[j](1);
    np(2) = (1-lam) * edgepoints[j-1](2) + lam * edgepoints[j](2);
      
      
    thispi = 0;
    if (i == ne)
    for (j = 1; j <= mesh.GetNP(); j++)
    if (Dist(mesh.Point(j), np) < 1e-6)
    thispi = j;
      
    if (!thispi)
    {
    ProjectToEdge (surf1, surf2, np);
    thispi = mesh.AddPoint (np);
    }
    */

    // (*testout) << "short edge " << pi1 << " - " << pi2 << endl;
  
    for (int k = 1; k <= refedges.Size(); k++)
      {
	if (refedgesinv.Get(k))
	  {
	    seg.p1 = pi1;
	    seg.p2 = pi2;
	  }
	else
	  {
	    seg.p1 = pi2;
	    seg.p2 = pi1;
	  }

	seg.si = refedges.Get(k).si;
	seg.domin = refedges.Get(k).domin;
	seg.domout = refedges.Get(k).domout;
	seg.tlosurf = refedges.Get(k).tlosurf;
	seg.edgenr = refedges.Get(k).edgenr;
	seg.surfnr1 = refedges.Get(k).surfnr1;
	seg.surfnr2 = refedges.Get(k).surfnr2;
	seg.seginfo = 0;
	if (k == 1) seg.seginfo = (refedgesinv.Get(k)) ? 2 : 1;
	mesh.AddSegment (seg);
	//	  (*testout) << "add seg " << seg.p1 << "-" << seg.p2 << endl;
      }
  }
  






  void EdgeCalculation :: 
  CopyEdge (const ARRAY<Segment> & refedges,
	    const ARRAY<bool> & refedgesinv,
	    int copyfromedge, 
	    const Point<3> & fromstart, const Point<3> & fromend,
	    const Point<3> & tostart, const Point<3> & toend,
	    int copyedgeidentification, 
	    int layer,
	    Mesh & mesh)
  {
    int k;
    PointIndex pi;

    double size = geometry.MaxSize();
    
    // copy start and end points
    for (int i = 1; i <= 2; i++)
      {
	Point<3> fromp =
	  (i == 1) ? fromstart : fromend;
	Point<3> top =
	  (i == 1) ? tostart : toend;
      
	PointIndex frompi = -1;
	PointIndex topi = -1;
	for (pi = PointIndex::BASE; 
	     pi < mesh.GetNP()+PointIndex::BASE; pi++)
	  {
	    if (Dist2 (mesh[pi], fromp) <= 1e-16*size)
	      frompi = pi;
	    if (Dist2 (mesh[pi], top) <= 1e-16*size)
	      topi = pi;
	  }

	
	if (topi == -1)
	  {
	    topi = mesh.AddPoint (top, layer, FIXEDPOINT);
	    meshpoint_tree -> Insert (top, topi);
	  }

	const Identification & csi = 
	  (*geometry.identifications.Get(copyedgeidentification));


	if (csi.Identifyable (mesh[frompi], mesh[topi]))
	  mesh.GetIdentifications().Add(frompi, topi, copyedgeidentification);
	else if (csi.Identifyable (mesh[topi], mesh[frompi]))
	  mesh.GetIdentifications().Add(topi, frompi, copyedgeidentification);
	else
	  {
	    cerr << "edgeflw.cpp: should identify, but cannot";
	    exit(1);
	  }
#ifdef DEVELOP
	(*testout) << "adding identification " << mesh[frompi] << "; " << mesh[topi]
		   << " (id " << copyedgeidentification <<")" << endl;
#endif


	/*
	  (*testout) << "Add Identification from CopyEdge, p1 = " 
	  << mesh[PointIndex(frompi)] << ", p2 = " 
	  << mesh[PointIndex(topi)] << endl;

	  mesh.GetIdentifications().Add(frompi, topi, copyedgeidentification);
	*/
      }

    int oldns = mesh.GetNSeg();
    for (int i = 1; i <= oldns; i++)
      {
	// real copy, since array might be reallocated !!
	const Segment oldseg = mesh.LineSegment(i);
	if (oldseg.edgenr != copyfromedge)
	  continue;
	if (oldseg.seginfo == 0)
	  continue;

	int pi1 = oldseg.p1;
	int pi2 = oldseg.p2;

	int npi1 = geometry.identifications.Get(copyedgeidentification)
	  -> GetIdentifiedPoint (mesh, pi1);
	int npi2 = geometry.identifications.Get(copyedgeidentification)
	  -> GetIdentifiedPoint (mesh, pi2);

	//(*testout) << "copy edge, pts = " << npi1 << " - " << npi2 << endl;

	Segment seg;

	for (k = 1; k <= refedges.Size(); k++)
	  {
	    bool inv = refedgesinv.Get(k);

	    // other edge is inverse
	    if (oldseg.seginfo == 1)
	      inv = !inv;

	    //	  (*testout) << "inv, now = " << inv << endl;

	    if (inv)
	      {
		seg.p1 = npi1;
		seg.p2 = npi2;
	      }
	    else
	      {
		seg.p1 = npi2;
		seg.p2 = npi1;
	      }
	    seg.si = refedges.Get(k).si;
	    seg.domin = refedges.Get(k).domin;
	    seg.domout = refedges.Get(k).domout;
	    seg.tlosurf = refedges.Get(k).tlosurf;
	    seg.edgenr = refedges.Get(k).edgenr;
	    seg.surfnr1 = refedges.Get(k).surfnr1;
	    seg.surfnr2 = refedges.Get(k).surfnr2;
	    seg.seginfo = 0;
	    if (k == 1) seg.seginfo = refedgesinv.Get(k) ? 2 : 1;
	    mesh.AddSegment (seg);
	    //	  (*testout) << "copy seg " << seg.p1 << "-" << seg.p2 << endl;
#ifdef DEVELOP

	    (*testout) << "copy seg, face = " << seg.si << ": " 
		       << " inv = " << inv << ", refinv = " << refedgesinv.Get(k)
		       << mesh.Point(seg.p1) << ", " << mesh.Point(seg.p2) << endl;
#endif

	  }
      
      }   
  }
  






  void EdgeCalculation :: 
  FindClosedSurfaces (double h, Mesh & mesh)
  {
    // if there is no special point at a sphere, one has to add a segment pair
  
    int i, j; 
    int nsol; 
    int nsurf = geometry.GetNSurf();
    int layer(0);

    BitArray pointatsurface (nsurf);
    Point<3> p1, p2;
    Vec<3> nv, tv;
    Solid * tansol;
    ARRAY<int> tansurfind;
    //  const Solid * sol;

    double size = geometry.MaxSize();
    nsol = geometry.GetNTopLevelObjects();
    

    pointatsurface.Clear();
  
    /*
      for (i = 1; i <= specpoints.Size(); i++)
      {
      int classrep;

      classrep = geometry.GetSurfaceClassRepresentant (specpoints[i].s1);
      pointatsurface.Set (classrep);
      classrep = geometry.GetSurfaceClassRepresentant (specpoints[i].s2);
      pointatsurface.Set (classrep);
      //      pointatsurface.Set (specpoints[i].s1);
      //      pointatsurface.Set (specpoints[i].s2);
      }
    */
    for (i = 1; i <= mesh.GetNSeg(); i++)
      {
	const Segment & seg = mesh.LineSegment(i);
	int classrep;

#ifdef DEVELOP      
	(*testout) << seg.surfnr1 << ", " << seg.surfnr2 << ", si = " << seg.si << endl;
#endif
	classrep = geometry.GetSurfaceClassRepresentant (seg.si);

	pointatsurface.Set (classrep);
      }

  
    for (i = 0; i < nsurf; i++)
      {
	int classrep = geometry.GetSurfaceClassRepresentant (i);

	if (!pointatsurface.Test(classrep))
	  {
	    const Surface * s = geometry.GetSurface(i);
	    p1 = s -> GetSurfacePoint();
	    nv = s -> GetNormalVector (p1);
		    
	    double hloc = 
	      min2 (s->LocH (p1, 3, 1, h), mesh.GetH(p1));

	    tv = nv.GetNormal ();
	    tv *=  (hloc / tv.Length());
	    p2 = p1 + tv;
	    s->Project (p2);
	  
		    
	    Segment seg1;
	    seg1.si = i;
	    seg1.domin = -1;
	    seg1.domout = -1;

	    Segment seg2;
	    seg2.si = i;
	    seg2.domin = -1;
	    seg2.domout = -1;

	    seg1.surfnr1 = i;
	    seg2.surfnr1 = i;
	    seg1.surfnr2 = i;
	    seg2.surfnr2 = i;

	    for (j = 0; j < nsol; j++)
	      {
		if (geometry.GetTopLevelObject(j)->GetSurface())
		  continue;
		  
		const Solid * sol = geometry.GetTopLevelObject(j)->GetSolid();
		sol -> TangentialSolid (p1, tansol, tansurfind, ideps*size);
		layer = geometry.GetTopLevelObject(j)->GetLayer();

		
		if (tansol)
		  {
		    tansol -> GetSurfaceIndices (tansurfind);
		
		    if (tansurfind.Size() == 1 && tansurfind.Get(1) == i)
		      {
			if (!tansol->VectorIn(p1, nv))
			  {
			    seg1.domin = j;
			    seg2.domin = j;
			    seg1.tlosurf = j;
			    seg2.tlosurf = j;
			  }
			else
			  {
			    seg1.domout = j;
			    seg2.domout = j;
			    seg1.tlosurf = j;
			    seg2.tlosurf = j;
			  }
			//        seg.s2 = i;
			//        seg.invs1 = surfaces[i] -> Inverse();
			//        seg.invs2 = ! (surfaces[i] -> Inverse());
		      }
		    delete tansol;
		  }
	      }


	    if (seg1.domin != -1 || seg1.domout != -1)
	      {
		mesh.AddPoint (p1, layer, EDGEPOINT);
		mesh.AddPoint (p2, layer, EDGEPOINT);
		seg1.p1 = mesh.GetNP()-1;
		seg1.p2 = mesh.GetNP();
		seg2.p2 = mesh.GetNP()-1;
		seg2.p1 = mesh.GetNP();
		seg1.geominfo[0].trignum = 1;
		seg1.geominfo[1].trignum = 1;
		seg2.geominfo[0].trignum = 1;
		seg2.geominfo[1].trignum = 1;
		mesh.AddSegment (seg1);
		mesh.AddSegment (seg2);

		PrintMessage (5, "Add line segment to smooth surface");

#ifdef DEVELOP
		(*testout) << "Add segment at smooth surface " << i;
		if (i != classrep) (*testout) << ", classrep = " << classrep;
		(*testout) << ": "
			   << mesh.Point (mesh.GetNP()-1) << " - "
			   << mesh.Point (mesh.GetNP()) << endl;
#endif
	      }
	  }
      }
  }

}
