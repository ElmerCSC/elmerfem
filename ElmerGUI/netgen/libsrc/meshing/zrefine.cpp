#include <mystdlib.h>
#include "meshing.hpp"

#include <csg.hpp>

namespace netgen
{

  // find singular edges
  void SelectSingularEdges (const Mesh & mesh, const CSGeometry & geom, 
			    INDEX_2_HASHTABLE<int> & singedges,
			    ZRefinementOptions & opt)
  {
    int i, j;

    // edges selected in csg input file
    for (i = 1; i <= geom.singedges.Size(); i++)
      {
	//if(geom.singedges.Get(i)->maxhinit > 0)
	//  continue; //!!!!

	const SingularEdge & se = *geom.singedges.Get(i);
	for (j = 1; j <= se.segms.Size(); j++)
	  {
	    INDEX_2 i2 = se.segms.Get(j);
	    singedges.Set (i2, 1);
	  }
      }

    // edges interactively selected
    for (i = 1; i <= mesh.GetNSeg(); i++)
      {
	const Segment & seg = mesh.LineSegment(i);
	if (seg.singedge_left || seg.singedge_right)
	  {
	    INDEX_2 i2(seg.p1, seg.p2);
	    i2.Sort();
	    singedges.Set (i2, 1);
	  }
      }
  }


  /**
     Convert elements (vol-tets, surf-trigs) into prisms/quads
  */
  void MakePrismsSingEdge (Mesh & mesh, INDEX_2_HASHTABLE<int> & singedges)
  {
    int i, j, k;

    // volume elements
    for (i = 1; i <= mesh.GetNE(); i++)
      {
	Element & el = mesh.VolumeElement(i);
	if (el.GetType() != TET) continue;

	for (j = 1; j <= 3; j++)
	  for (k = j+1; k <= 4; k++)
	    {
	      INDEX_2 edge(el.PNum(j), el.PNum(k));
	      edge.Sort();
	      if (singedges.Used (edge))
		{
		  int pi3 = 1, pi4 = 1;
		  while (pi3 == j || pi3 == k) pi3++;
		  pi4 = 10 - j - k - pi3;
		
		  int p3 = el.PNum(pi3);
		  int p4 = el.PNum(pi4);

		  el.SetType(PRISM);
		  el.PNum(1) = edge.I1();
		  el.PNum(2) = p3;
		  el.PNum(3) = p4;
		  el.PNum(4) = edge.I2();
		  el.PNum(5) = p3;
		  el.PNum(6) = p4;
		}
	    }
      }

    // surface elements
    for (i = 1; i <= mesh.GetNSE(); i++)
      {
	Element2d & el = mesh.SurfaceElement(i);
	if (el.GetType() != TRIG) continue;

	for (j = 1; j <= 3; j++)
	  {
	    k = (j % 3) + 1;
	    INDEX_2 edge(el.PNum(j), el.PNum(k));
	    edge.Sort();

	    if (singedges.Used (edge))
	      {
		int pi3 = 6-j-k;
		int p3 = el.PNum(pi3);
		int p1 = el.PNum(j);
		int p2 = el.PNum(k);

		el.SetType(QUAD);
		el.PNum(1) = p2;
		el.PNum(2) = p3;
		el.PNum(3) = p3;
		el.PNum(4) = p1;
	      }
	  }
      }
  }


  /*
    Convert tets and pyramids next to close (identified) points into prisms
  */
  void MakePrismsClosePoints (Mesh & mesh)
  {
    int i, j, k;
    for (i = 1; i <= mesh.GetNE(); i++)
      {
	Element & el = mesh.VolumeElement(i);
	if (el.GetType() == TET)
	  {
	    for (j = 1; j <= 3; j++)
	      for (k = j+1; k <= 4; k++)
		{
		  INDEX_2 edge(el.PNum(j), el.PNum(k));
		  edge.Sort();
		  if (mesh.GetIdentifications().GetSymmetric (el.PNum(j), el.PNum(k)))
		    {
		      int pi3 = 1, pi4 = 1;
		      while (pi3 == j || pi3 == k) pi3++;
		      pi4 = 10 - j - k - pi3;
		    
		      int p3 = el.PNum(pi3);
		      int p4 = el.PNum(pi4);
		    
		      el.SetType(PRISM);
		      el.PNum(1) = edge.I1();
		      el.PNum(2) = p3;
		      el.PNum(3) = p4;
		      el.PNum(4) = edge.I2();
		      el.PNum(5) = p3;
		      el.PNum(6) = p4;
		    }
		}
	  }

	if (el.GetType() == PYRAMID)
	  {
	    // pyramid, base face = 1,2,3,4
	  
	    for (j = 0; j <= 1; j++)
	      {
		int pi1 = el.PNum( (j+0) % 4 + 1);
		int pi2 = el.PNum( (j+1) % 4 + 1);
		int pi3 = el.PNum( (j+2) % 4 + 1);
		int pi4 = el.PNum( (j+3) % 4 + 1);
		int pi5 = el.PNum(5);

		INDEX_2 edge1(pi1, pi4);
		INDEX_2 edge2(pi2, pi3);
		edge1.Sort();
		edge2.Sort();
		if (mesh.GetIdentifications().GetSymmetric (pi1, pi4) &&
		    mesh.GetIdentifications().GetSymmetric (pi2, pi3))
		  {
		    //int p3 = el.PNum(pi3);
		    //int p4 = el.PNum(pi4);
		  
		    el.SetType(PRISM);
		    el.PNum(1) = pi1;
		    el.PNum(2) = pi2;
		    el.PNum(3) = pi5;
		    el.PNum(4) = pi4;
		    el.PNum(5) = pi3;
		    el.PNum(6) = pi5;
		  }
	      }
	  }
      }
  
    for (i = 1; i <= mesh.GetNSE(); i++)
      {
	Element2d & el = mesh.SurfaceElement(i);
	if (el.GetType() != TRIG) continue;

	for (j = 1; j <= 3; j++)
	  {
	    k = (j % 3) + 1;
	    INDEX_2 edge(el.PNum(j), el.PNum(k));
	    edge.Sort();

	    if (mesh.GetIdentifications().GetSymmetric (el.PNum(j), el.PNum(k)))
	      {
		int pi3 = 6-j-k;
		int p3 = el.PNum(pi3);
		int p1 = el.PNum(j);
		int p2 = el.PNum(k);

		el.SetType(QUAD);
		el.PNum(1) = p2;
		el.PNum(2) = p3;
		el.PNum(3) = p3;
		el.PNum(4) = p1;
	      }
	  }
      }
  }



#ifdef OLD
  void MakeCornerNodes (Mesh & mesh,
			INDEX_HASHTABLE<int> & cornernodes)
  {
    int i, j;
    int nseg = mesh.GetNSeg();
    ARRAY<int> edgesonpoint(mesh.GetNP());
    for (i = 1; i <= mesh.GetNP(); i++)
      edgesonpoint.Elem(i) = 0;

    for (i = 1; i <= nseg; i++)
      {
	for (j = 1; j <= 2; j++)
	  {
	    int pi = (j == 1) ? 
	      mesh.LineSegment(i).p1 :
	      mesh.LineSegment(i).p2;
	    edgesonpoint.Elem(pi)++;
	  }
      }

    /*
      cout << "cornernodes: ";
      for (i = 1; i <= edgesonpoint.Size(); i++)
      if (edgesonpoint.Get(i) >= 6)
      {
      cornernodes.Set (i, 1);
      cout << i << " ";
      }
      cout << endl;
    */
    //  cornernodes.Set (5, 1);
  }
#endif


  void RefinePrisms (Mesh & mesh, const CSGeometry * geom, 
		     ZRefinementOptions & opt)
  {
    int i, j;
    bool found, change;
    int cnt = 0;


    // markers for z-refinement:  p1, p2, levels  
    // p1-p2 is an edge to be refined
    ARRAY<INDEX_3> ref_uniform;
    ARRAY<INDEX_3> ref_singular;
    ARRAY<INDEX_4 > ref_slices;

    BitArray first_id(geom->identifications.Size());
    first_id.Set();

  
    INDEX_2_HASHTABLE<int> & identpts = 
      mesh.GetIdentifications().GetIdentifiedPoints ();

    if (&identpts)
      {
	for (i = 1; i <= identpts.GetNBags(); i++)
	  for (j = 1; j <= identpts.GetBagSize(i); j++)
	    {
	      INDEX_2 pair;
	      int idnr;
	      identpts.GetData(i, j, pair, idnr);
	      const CloseSurfaceIdentification * csid = 
		dynamic_cast<const CloseSurfaceIdentification*> 
		(geom->identifications.Get(idnr));
	      if (csid)
		{
		  if (!csid->GetSlices().Size())
		    {
		      if (first_id.Test (idnr))
			{
			  first_id.Clear(idnr);
			  ref_uniform.Append (INDEX_3 (pair.I1(), pair.I2(), csid->RefLevels()));
			  ref_singular.Append (INDEX_3 (pair.I1(), pair.I2(), csid->RefLevels1()));
			  ref_singular.Append (INDEX_3 (pair.I2(), pair.I1(), csid->RefLevels2()));
			}
		    }
		  else
		    {   
		      //const ARRAY<double> & slices = csid->GetSlices();
		      INDEX_4 i4;
		      i4[0] = pair.I1();
		      i4[1] = pair.I2();
		      i4[2] = idnr;
		      i4[3] = csid->GetSlices().Size();
		      ref_slices.Append (i4);
		    }
		}
	    }
      }

  
  
    ARRAY<EdgePointGeomInfo> epgi;

    while (1)
      {
	cnt++;
	PrintMessage (3, "Z-Refinement, level = ", cnt);
	INDEX_2_HASHTABLE<int> refedges(mesh.GetNSE()+1);


	found = 0;
	// mark prisms due to close surface flags:
	int oldsize = ref_uniform.Size();
	for (i = 1; i <= oldsize; i++)
	  {
	    int pi1 = ref_uniform.Get(i).I1();
	    int pi2 = ref_uniform.Get(i).I2();
	    int levels = ref_uniform.Get(i).I3();

	    if (levels > 0)
	      {
		const Point3d & p1 = mesh.Point(pi1);
		const Point3d & p2 = mesh.Point(pi2);
		int npi(0);
	      
		INDEX_2 edge(pi1, pi2);
		edge.Sort();
		if (!refedges.Used(edge))
		  {
		    Point3d np = Center (p1, p2);
		    npi = mesh.AddPoint (np);
		    refedges.Set (edge, npi);
		    found = 1;
		  }

		ref_uniform.Elem(i) = INDEX_3(pi1, npi, levels-1);
		ref_uniform.Append (INDEX_3(pi2, npi, levels-1));
	      }
	  }
	for (i = 1; i <= ref_singular.Size(); i++)
	  {
	    int pi1 = ref_singular.Get(i).I1();
	    int pi2 = ref_singular.Get(i).I2();
	    int levels = ref_singular.Get(i).I3();

	    if (levels > 0)
	      {
		const Point3d & p1 = mesh.Point(pi1);
		const Point3d & p2 = mesh.Point(pi2);
		int npi;
	      
		INDEX_2 edge(pi1, pi2);
		edge.Sort();
		if (!refedges.Used(edge))
		  {
		    Point3d np = Center (p1, p2);
		    npi = mesh.AddPoint (np);
		    refedges.Set (edge, npi);
		    found = 1;
		  }
		else
		  npi = refedges.Get (edge);

		ref_singular.Elem(i) = INDEX_3(pi1, npi, levels-1);
	      }
	  }

	for (i = 1; i <= ref_slices.Size(); i++)
	  {
	    int pi1 = ref_slices.Get(i)[0];
	    int pi2 = ref_slices.Get(i)[1];
	    int idnr = ref_slices.Get(i)[2];
	    int slicenr = ref_slices.Get(i)[3];

	    if (slicenr > 0)
	      {
		const Point3d & p1 = mesh.Point(pi1);
		const Point3d & p2 = mesh.Point(pi2);
		int npi;

		const CloseSurfaceIdentification * csid = 
		  dynamic_cast<const CloseSurfaceIdentification*> 
		  (geom->identifications.Get(idnr));

	      
		INDEX_2 edge(pi1, pi2);
		edge.Sort();
		if (!refedges.Used(edge))
		  {
		    const ARRAY<double> & slices = csid->GetSlices();
		    //(*testout) << "idnr " << idnr << " i " << i << endl;
		    //(*testout) << "slices " << slices << endl;
		    double slicefac = slices.Get(slicenr);
		    double slicefaclast = 
		      (slicenr == slices.Size()) ? 1 : slices.Get(slicenr+1);
		    
		    Point3d np = p1 + (slicefac / slicefaclast) * (p2-p1);
		    //(*testout) << "slicenr " << slicenr << " slicefac " << slicefac << " quot " << (slicefac / slicefaclast) << " np " << np << endl;
		    npi = mesh.AddPoint (np);
		    refedges.Set (edge, npi);
		    found = 1;
		  }
		else
		  npi = refedges.Get (edge);
		
		ref_slices.Elem(i)[1] = npi;
		ref_slices.Elem(i)[3] --;
	      }
	  }




	for (i = 1; i <= mesh.GetNE(); i++)
	  {
	    Element & el = mesh.VolumeElement (i);
	    if (el.GetType() != PRISM)
	      continue;

	    for (j = 1; j <= 3; j++)
	      {
		int pi1 = el.PNum(j);
		int pi2 = el.PNum(j+3);
		const Point3d & p1 = mesh.Point(pi1);
		const Point3d & p2 = mesh.Point(pi2);

		bool ref = 0;

		/*
		  if (Dist (p1, p2) > mesh.GetH (Center (p1, p2)))
		  ref = 1;
		*/

		/*
		  if (cnt <= opt.minref)
		  ref = 1;
		*/

		/*
		  if ((pi1 == 460 || pi2 == 460 ||
		  pi1 == 461 || pi2 == 461) && cnt <= 8) ref = 1;
		*/
		if (ref == 1)
		  {
		    INDEX_2 edge(pi1, pi2);
		    edge.Sort();
		    if (!refedges.Used(edge))
		      {
			Point3d np = Center (p1, p2);
			int npi = mesh.AddPoint (np);
			refedges.Set (edge, npi);
			found = 1;
		      }
		  }
	      }
	  }
      
	if (!found) break;

	// build closure:
	PrintMessage (5, "start closure");
	do
	  {
	    PrintMessage (5, "start loop");
	    change = 0;
	    for (i = 1; i <= mesh.GetNE(); i++)
	      {
		Element & el = mesh.VolumeElement (i);
		if (el.GetType() != PRISM)
		  continue;
	      
		bool hasref = 0, hasnonref = 0;
		for (j = 1; j <= 3; j++)
		  {
		    int pi1 = el.PNum(j);
		    int pi2 = el.PNum(j+3);
		    if (pi1 != pi2)
		      {
			INDEX_2 edge(pi1, pi2);
			edge.Sort();
			if (refedges.Used(edge))
			  hasref = 1;
			else 
			  hasnonref = 1;
		      }
		  }

		if (hasref && hasnonref)
		  {
		    //		  cout << "el " << i << " in closure" << endl;
		    change = 1;
		    for (j = 1; j <= 3; j++)
		      {
			int pi1 = el.PNum(j);
			int pi2 = el.PNum(j+3);
			const Point3d & p1 = mesh.Point(pi1);
			const Point3d & p2 = mesh.Point(pi2);
		      
			INDEX_2 edge(pi1, pi2);
			edge.Sort();
			if (!refedges.Used(edge))
			  {
			    Point3d np = Center (p1, p2);
			    int npi = mesh.AddPoint (np);
			    refedges.Set (edge, npi);
			  }
		      }
		  }
	      }
	  }
	while (change);

	PrintMessage (5, "Do segments");

	//      (*testout) << "closure formed, np = " << mesh.GetNP() << endl;

	int oldns = mesh.GetNSeg();

	for (i = 1; i <= oldns; i++)
	  {
	    const Segment & el = mesh.LineSegment(i);

	    INDEX_2 i2(el.p1, el.p2);
	    i2.Sort();
	  
	    int pnew;
	    EdgePointGeomInfo ngi;
      
	    if (refedges.Used(i2))
	      {
		pnew = refedges.Get(i2);
		//	      ngi = epgi.Get(pnew);
	      }
	    else
	      {
		continue;

		// 	      Point3d pb;

		// 	      /*
		// 	      geom->PointBetween (mesh.Point (el.p1),
		// 				  mesh.Point (el.p2),
		// 				  el.surfnr1, el.surfnr2,
		// 				  el.epgeominfo[0], el.epgeominfo[1],
		// 				  pb, ngi);
		// 	      */
		// 	      pb = Center (mesh.Point (el.p1), mesh.Point (el.p2));

		// 	      pnew = mesh.AddPoint (pb);
	      
		// 	      refedges.Set (i2, pnew);
	      
		// 	      if (pnew > epgi.Size())
		// 		epgi.SetSize (pnew);
		// 	      epgi.Elem(pnew) = ngi;
	      }
	  
	    Segment ns1 = el;
	    Segment ns2 = el;
	    ns1.p2 = pnew;
	    ns1.epgeominfo[1] = ngi;
	    ns2.p1 = pnew;
	    ns2.epgeominfo[0] = ngi;

	    mesh.LineSegment(i) = ns1;
	    mesh.AddSegment (ns2);
	  }
      
	PrintMessage (5, "Segments done, NSeg = ", mesh.GetNSeg());

	// do refinement
	int oldne = mesh.GetNE();
	for (i = 1; i <= oldne; i++)
	  {
	    Element & el = mesh.VolumeElement (i);
	    if (el.GetNP() != 6)
	      continue;

	    int npi[3];
	    for (j = 1; j <= 3; j++)
	      {
		int pi1 = el.PNum(j);
		int pi2 = el.PNum(j+3);

		if (pi1 == pi2)
		  npi[j-1] = pi1;
		else
		  {
		    INDEX_2 edge(pi1, pi2);
		    edge.Sort();
		    if (refedges.Used (edge))
		      npi[j-1] = refedges.Get(edge);
		    else
		      {
			/*
			  (*testout) << "ERROR: prism " << i << " has hanging node !!" 
			  << ", edge = " << edge << endl;
			  cerr << "ERROR: prism " << i << " has hanging node !!" << endl;
			*/
			npi[j-1] = 0;
		      }
		  }
	      }

	    if (npi[0])
	      {
		Element nel1(6), nel2(6);
		for (j = 1; j <= 3; j++)
		  {
		    nel1.PNum(j) = el.PNum(j);
		    nel1.PNum(j+3) = npi[j-1];
		    nel2.PNum(j) = npi[j-1];
		    nel2.PNum(j+3) = el.PNum(j+3);
		  }
		nel1.SetIndex (el.GetIndex());
		nel2.SetIndex (el.GetIndex());
		mesh.VolumeElement (i) = nel1;
		mesh.AddVolumeElement (nel2);
	      }
	  }

      
	PrintMessage (5, "Elements done, NE = ", mesh.GetNE());


	// do surface elements
	int oldnse = mesh.GetNSE();
	//      cout << "oldnse = " << oldnse << endl;
	for (i = 1; i <= oldnse; i++)
	  {
	    Element2d & el = mesh.SurfaceElement (i);
	    if (el.GetType() != QUAD)
	      continue;

	    int index = el.GetIndex();
	    int npi[2];
	    for (j = 1; j <= 2; j++)
	      {
		int pi1, pi2;

		if (j == 1)
		  {
		    pi1 = el.PNum(1);
		    pi2 = el.PNum(4);
		  }
		else
		  {
		    pi1 = el.PNum(2);
		    pi2 = el.PNum(3);
		  }

		if (pi1 == pi2)
		  npi[j-1] = pi1;
		else
		  {
		    INDEX_2 edge(pi1, pi2);
		    edge.Sort();
		    if (refedges.Used (edge))
		      npi[j-1] = refedges.Get(edge);
		    else
		      {
			npi[j-1] = 0;
		      }
		  }
	      }

	    if (npi[0])
	      {
		Element2d nel1(QUAD), nel2(QUAD);
		for (j = 1; j <= 4; j++)
		  {
		    nel1.PNum(j) = el.PNum(j);
		    nel2.PNum(j) = el.PNum(j);
		  }
		nel1.PNum(3) = npi[1];
		nel1.PNum(4) = npi[0];
		nel2.PNum(1) = npi[0];
		nel2.PNum(2) = npi[1];
		/*
		  for (j = 1; j <= 2; j++)
		  {
		  nel1.PNum(j) = el.PNum(j);
		  nel1.PNum(j+2) = npi[j-1];
		  nel2.PNum(j) = npi[j-1];
		  nel2.PNum(j+2) = el.PNum(j+2);
		  }
		*/
		nel1.SetIndex (el.GetIndex());
		nel2.SetIndex (el.GetIndex());

		mesh.SurfaceElement (i) = nel1;
		mesh.AddSurfaceElement (nel2);

		int si = mesh.GetFaceDescriptor (index).SurfNr();

		Point<3> hp = mesh.Point(npi[0]);
		geom->GetSurface(si)->Project (hp);
		mesh.Point (npi[0]).SetPoint (hp);

		hp = mesh.Point(npi[1]);
		geom->GetSurface(si)->Project (hp);
		mesh.Point (npi[1]).SetPoint (hp);

		//	      geom->GetSurface(si)->Project (mesh.Point(npi[0]));
		//	      geom->GetSurface(si)->Project (mesh.Point(npi[1]));
	      }
	  }

	PrintMessage (5, "Surface elements done, NSE = ", mesh.GetNSE());

      }
  }



  void ZRefinement (Mesh & mesh, const CSGeometry * geom,
		    ZRefinementOptions & opt)
  {
    INDEX_2_HASHTABLE<int> singedges(mesh.GetNSeg());

    SelectSingularEdges (mesh, *geom, singedges, opt);
    //MakePrismsSingEdge (mesh, singedges);
    MakePrismsClosePoints (mesh);

    RefinePrisms (mesh, geom, opt);
  }



  ZRefinementOptions :: ZRefinementOptions()
  {
    minref = 0;
  }

}
