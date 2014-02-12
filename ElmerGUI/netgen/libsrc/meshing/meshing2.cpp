#include <mystdlib.h>
#include "meshing.hpp"

namespace netgen
{
  static void glrender (int wait);


  // global variable for visualization

  static ARRAY<Point3d> locpoints;
  static ARRAY<int> legalpoints;
  static ARRAY<Point2d> plainpoints;
  static ARRAY<int> plainzones;
  static ARRAY<INDEX_2> loclines;
  static int geomtrig;
  //static const char * rname;
  static int cntelem, trials, nfaces;
  static int oldnl;
  static int qualclass;
 

  Meshing2 :: Meshing2 (const Box<3> & aboundingbox)
  {
    boundingbox = aboundingbox;

    LoadRules (NULL);
    // LoadRules ("rules/quad.rls");
    // LoadRules ("rules/triangle.rls");

    adfront = new AdFront2(boundingbox);
    starttime = GetTime();

    maxarea = -1;
  }


  Meshing2 :: ~Meshing2 ()
  {
    delete adfront;
    for (int i = 0; i < rules.Size(); i++)
      delete rules[i];
  }

  void Meshing2 :: AddPoint (const Point3d & p, PointIndex globind, 
			     MultiPointGeomInfo * mgi,
			     bool pointonsurface)
  {
    //(*testout) << "add point " << globind << endl;
    adfront ->AddPoint (p, globind, mgi, pointonsurface);
  }

  void Meshing2 :: AddBoundaryElement (int i1, int i2,
				       const PointGeomInfo & gi1, const PointGeomInfo & gi2)
  {
    //    (*testout) << "add line " << i1 << " - " << i2 << endl;
    if (!gi1.trignum || !gi2.trignum)
      {
	PrintSysError ("addboundaryelement: illegal geominfo");
      }
    adfront -> AddLine (i1-1, i2-1, gi1, gi2);
  }



  void Meshing2 :: StartMesh ()
  {
    foundmap.SetSize (rules.Size());
    canuse.SetSize (rules.Size());
    ruleused.SetSize (rules.Size());

    foundmap = 0;
    canuse = 0;
    ruleused = 0;

    cntelem = 0;
    trials = 0;
  }

  void Meshing2 :: EndMesh ()
  {
    for (int i = 0; i < ruleused.Size(); i++)
      (*testout) << setw(4) << ruleused[i]
		 << " times used rule " << rules[i] -> Name() << endl;
  }

  void Meshing2 :: SetStartTime (double astarttime)
  {
    starttime = astarttime;
  }

  
  void Meshing2 :: SetMaxArea (double amaxarea)
  {
    maxarea = amaxarea;
  }


  double Meshing2 :: CalcLocalH (const Point3d & /* p */, double gh) const
  {
    return gh;
  }

  // should be class variables !!(?)
  static Vec3d ex, ey;
  static Point3d globp1;

  void Meshing2 :: DefineTransformation (const Point3d & p1, const Point3d & p2,
					 const PointGeomInfo * geominfo1,
					 const PointGeomInfo * geominfo2)
  {
    globp1 = p1;
    ex = p2 - p1;
    ex /= ex.Length();
    ey.X() = -ex.Y();
    ey.Y() =  ex.X();
    ey.Z() = 0;
  }

  void Meshing2 :: TransformToPlain (const Point3d & locpoint, 
				     const MultiPointGeomInfo & geominf,
				     Point2d & plainpoint, double h, int & zone)
  {
    Vec3d p1p (globp1, locpoint);

    //    p1p = locpoint - globp1;
    p1p /= h;
    plainpoint.X() = p1p * ex;
    plainpoint.Y() = p1p * ey;
    zone = 0;
  }

  int Meshing2 :: TransformFromPlain (Point2d & plainpoint,
				      Point3d & locpoint, 
				      PointGeomInfo & gi, 
				      double h)
  {
    Vec3d p1p;
    gi.trignum = 1;

    p1p = plainpoint.X() * ex + plainpoint.Y() * ey;
    p1p *= h;
    locpoint = globp1 + p1p;
    return 0;
  }


  int Meshing2 :: BelongsToActiveChart (const Point3d & p, 
					const PointGeomInfo & gi)
  {
    return 1;
  }


  int Meshing2 :: ComputePointGeomInfo (const Point3d & p, PointGeomInfo & gi)
  {
    gi.trignum = 1;
    return 0;
  }


  int Meshing2 :: ChooseChartPointGeomInfo (const MultiPointGeomInfo & mpgi, 
					    PointGeomInfo & pgi)
  {
    pgi = mpgi.GetPGI(1);
    return 0;
  }



  int Meshing2 :: 
  IsLineVertexOnChart (const Point3d & p1, const Point3d & p2,
		       int endpoint, const PointGeomInfo & geominfo)
  {
    return 1;
  }

  void Meshing2 ::
  GetChartBoundary (ARRAY<Point2d> & points, 
		    ARRAY<Point3d> & points3d, 
		    ARRAY<INDEX_2> & lines, double h) const
  {
    points.SetSize (0);
    points3d.SetSize (0);
    lines.SetSize (0);
  }

  double Meshing2 :: Area () const
  {
    return -1;
  }





  MESHING2_RESULT Meshing2 :: GenerateMesh (Mesh & mesh, double gh, int facenr)
  {
    ARRAY<int> pindex, lindex;
    ARRAY<int> delpoints, dellines;

    ARRAY<PointGeomInfo> upgeominfo;  // unique info
    ARRAY<MultiPointGeomInfo> mpgeominfo;  // multiple info

    ARRAY<Element2d> locelements;

    int z1, z2, oldnp(-1);
    SurfaceElementIndex sei;
    bool found;
    int rulenr(-1);
    int globind;
    Point<3> p1, p2;

    const PointGeomInfo * blgeominfo1;
    const PointGeomInfo * blgeominfo2;

    bool morerisc;
    bool debugflag;

    double h, his, hshould;


    // test for 3d overlaps
    Box3dTree surfeltree (boundingbox.PMin(),
			  boundingbox.PMax());

    ARRAY<int> intersecttrias;
    ARRAY<Point3d> critpoints;

    // test for doubled edges
    //INDEX_2_HASHTABLE<int> doubleedge(300000);


    testmode = 0;

    StartMesh();

    ARRAY<Point2d> chartboundpoints;
    ARRAY<Point3d> chartboundpoints3d;
    ARRAY<INDEX_2> chartboundlines;

    // illegal points: points with more then 50 elements per node
    int maxlegalpoint(-1), maxlegalline(-1);
    ARRAY<int,PointIndex::BASE> trigsonnode;
    ARRAY<int,PointIndex::BASE> illegalpoint;

    trigsonnode.SetSize (mesh.GetNP());
    illegalpoint.SetSize (mesh.GetNP());

    trigsonnode = 0;
    illegalpoint = 0;
  

    double totalarea = Area ();
    double meshedarea = 0;

    // search tree for surface elements:
    for (sei = 0; sei < mesh.GetNSE(); sei++)
      {
	const Element2d & sel = mesh[sei];

	if (sel.IsDeleted()) continue;

	if (sel.GetIndex() == facenr)
	  {
	    Box<3> box;
	    box.Set ( mesh[sel[0]] );
	    box.Add ( mesh[sel[1]] );
	    box.Add ( mesh[sel[2]] );
	    surfeltree.Insert (box, sei);
	  }
      
	double trigarea = Cross ( mesh[sel[1]]-mesh[sel[0]],
				  mesh[sel[2]]-mesh[sel[0]] ).Length() / 2;


	if (sel.GetNP() == 4)
	  trigarea += Cross (Vec3d (mesh.Point (sel.PNum(1)),
				    mesh.Point (sel.PNum(3))),
			     Vec3d (mesh.Point (sel.PNum(1)),
				    mesh.Point (sel.PNum(4)))).Length() / 2;;
	meshedarea += trigarea;
      }


    const char * savetask = multithread.task;
    multithread.task = "Surface meshing";

    adfront ->SetStartFront ();


    int plotnexttrial = 999;

    double meshedarea_before = meshedarea;

    while (!adfront ->Empty() && !multithread.terminate)
      {
	if (multithread.terminate)
	  throw NgException ("Meshing stopped");

	// known for STL meshing
	if (totalarea > 0)
	  multithread.percent = 100 * meshedarea / totalarea;
	/*
	  else
	  multithread.percent = 0;
	*/

	locpoints.SetSize(0);
	loclines.SetSize(0);
	pindex.SetSize(0);
	lindex.SetSize(0);
	delpoints.SetSize(0);
	dellines.SetSize(0);
	locelements.SetSize(0);



	// plot statistics
	if (trials > plotnexttrial)
	  {
	    PrintMessage (5, 
			  "faces = ", nfaces,
			  " trials = ", trials,
			  " elements = ", mesh.GetNSE(),
			  " els/sec = ",
			  (mesh.GetNSE() / (GetTime() - starttime + 0.0001)));
	    plotnexttrial += 1000;
	  }


	// unique-pgi, multi-pgi
	upgeominfo.SetSize(0);
	mpgeominfo.SetSize(0);


	nfaces = adfront->GetNFL();
	trials ++;
    

	if (trials % 1000 == 0)
	  {
	    (*testout) << "\n";
	    for (int i = 1; i <= canuse.Size(); i++)
	      {
		(*testout) << foundmap.Get(i) << "/" 
			   << canuse.Get(i) << "/"
			   << ruleused.Get(i) << " map/can/use rule " << rules.Get(i)->Name() << "\n";
	      }
	    (*testout) << "\n";
	  }


	int baselineindex = adfront -> SelectBaseLine (p1, p2, blgeominfo1, blgeominfo2, qualclass);


	found = 1;

	his = Dist (p1, p2);

	Point3d pmid = Center (p1, p2);
	hshould = CalcLocalH (pmid, mesh.GetH (pmid));
	if (gh < hshould) hshould = gh;

	mesh.RestrictLocalH (pmid, hshould);

	h = hshould;

	double hinner = (3 + qualclass) * max2 (his, hshould);

	adfront ->GetLocals (baselineindex, locpoints, mpgeominfo, loclines, 
			     pindex, lindex, 2*hinner);
	//(*testout) << "h for locals: " << 2*hinner << endl;
	

	//(*testout) << "locpoints " << locpoints << endl;

	if (qualclass > mparam.giveuptol2d)
	  {
	    PrintMessage (3, "give up with qualclass ", qualclass);
	    PrintMessage (3, "number of frontlines = ", adfront->GetNFL());
	    // throw NgException ("Give up 2d meshing");
	    break;
	  }

	/*
	if (found && qualclass > 60)
	  {
	    found = 0;
	  }
	*/
	//      morerisc = ((qualclass > 20) && (qualclass % 2 == 1));
	//      morerisc = 1;
	morerisc = 0;


	PointIndex gpi1 = adfront -> GetGlobalIndex (pindex.Get(loclines[0].I1()));
	PointIndex gpi2 = adfront -> GetGlobalIndex (pindex.Get(loclines[0].I2()));


	debugflag = 
	  debugparam.haltsegment &&
	  ( (debugparam.haltsegmentp1 == gpi1) && 
	    (debugparam.haltsegmentp2 == gpi2) || 
	    (debugparam.haltsegmentp1 == gpi2) && 
	    (debugparam.haltsegmentp2 == gpi1)) ||
	  debugparam.haltnode &&
	  ( (debugparam.haltsegmentp1 == gpi1) ||
	    (debugparam.haltsegmentp2 == gpi1));

      
	if (debugparam.haltface && debugparam.haltfacenr == facenr)
	  {
	    debugflag = 1;
	    cout << "set debugflag" << endl;
	  }
	
	if (debugparam.haltlargequalclass && qualclass > 50)
	  debugflag = 1;


	// problem recognition !
	if (found && 
	    (gpi1 < illegalpoint.Size()+PointIndex::BASE) && 
	    (gpi2 < illegalpoint.Size()+PointIndex::BASE) )
	  {
	    if (illegalpoint[gpi1] || illegalpoint[gpi2])
	      found = 0;
	  }


	Point2d p12d, p22d;

	if (found)
	  {
	    oldnp = locpoints.Size();
	    oldnl = loclines.Size();
	  
	    if (debugflag)
	      (*testout) << "define new transformation" << endl;

	    DefineTransformation (p1, p2, blgeominfo1, blgeominfo2);
	  
	    plainpoints.SetSize (locpoints.Size());
	    plainzones.SetSize (locpoints.Size());

	    // (*testout) << endl;

	    //	    (*testout) << "3d->2d transformation" << endl;

	    for (int i = 1; i <= locpoints.Size(); i++)
	      {
		// (*testout) << "pindex(i) = " << pindex[i-1] << endl;
		TransformToPlain (locpoints.Get(i), 
				  mpgeominfo.Get(i),
				  plainpoints.Elem(i), h, plainzones.Elem(i));
		//		(*testout) << mpgeominfo.Get(i).GetPGI(1).u << " " << mpgeominfo.Get(i).GetPGI(1).v << " ";
		//		(*testout) << plainpoints.Get(i).X() << " " << plainpoints.Get(i).Y() << endl;
		//(*testout) << "transform " << locpoints.Get(i) << " to " << plainpoints.Get(i).X() << " " << plainpoints.Get(i).Y() << endl;
	      }
	    //	    (*testout) << endl << endl << endl;


	    p12d = plainpoints.Get(1);
	    p22d = plainpoints.Get(2);

	    /*
	    // last idea on friday
	    plainzones.Elem(1) = 0;
	    plainzones.Elem(2) = 0;
	    */


	    /*
	    // old netgen:
	    for (i = 2; i <= loclines.Size(); i++)  // don't remove first line
	    {
	    z1 = plainzones.Get(loclines.Get(i).I1());
	    z2 = plainzones.Get(loclines.Get(i).I2());
	      
	    if (z1 && z2 && (z1 != z2) || (z1 == -1) || (z2 == -1) )
	    {
	    loclines.DeleteElement(i);
	    lindex.DeleteElement(i);
	    oldnl--;
	    i--;
	    }
	    }

	    // 	  for (i = 1; i <= plainpoints.Size(); i++)
	    // 	    if (plainzones.Elem(i) == -1)
	    // 	      plainpoints.Elem(i) = Point2d (1e4, 1e4);
	    */
	  

	  
	    for (int i = 2; i <= loclines.Size(); i++)  // don't remove first line
	      {
		// (*testout) << "loclines(i) = " << loclines.Get(i).I1() << " - " << loclines.Get(i).I2() << endl;
		z1 = plainzones.Get(loclines.Get(i).I1());
		z2 = plainzones.Get(loclines.Get(i).I2());
	      
	      
		// one inner point, one outer
		if ( (z1 >= 0) != (z2 >= 0))
		  {
		    int innerp = (z1 >= 0) ? 1 : 2;
		    if (IsLineVertexOnChart (locpoints.Get(loclines.Get(i).I1()),
					     locpoints.Get(loclines.Get(i).I2()),
					     innerp,
					     adfront->GetLineGeomInfo (lindex.Get(i), innerp)))
		      // pgeominfo.Get(loclines.Get(i).I(innerp))))
		      {		

			if (!morerisc)
			  {
			    // use one end of line
			    int pini, pouti;
			    Vec2d v;
			  
			    pini = loclines.Get(i).I(innerp);
			    pouti = loclines.Get(i).I(3-innerp);
			  
			    Point2d pin (plainpoints.Get(pini));
			    Point2d pout (plainpoints.Get(pouti));
			    v = pout - pin;
			    double len = v.Length();
			    if (len <= 1e-6)
			      (*testout) << "WARNING(js): inner-outer: short vector" << endl;
			    else
			      v /= len;
			  
			    /*
			    // don't elongate line towards base-line !!
			    if (Vec2d (pin, p12d) * v > 0 && 
			    Vec2d (pin, p22d) * v > 0)
			    v *= -1;  
			    */

			    Point2d newpout = pin + 1000 * v;
			    newpout = pout;

			  
			    plainpoints.Append (newpout);
			    Point3d pout3d = locpoints.Get(pouti);
			    locpoints.Append (pout3d);

			    plainzones.Append (0);
			    pindex.Append (-1);
			    oldnp++;
			    loclines.Elem(i).I(3-innerp) = oldnp;
			  }
			else
			  plainzones.Elem(loclines.Get(i).I(3-innerp)) = 0;
			

			//		  (*testout) << "inner - outer correction" << endl;
		      }
		    else
		      {
			// remove line
			loclines.DeleteElement(i);
			lindex.DeleteElement(i);
			oldnl--;
			i--;
		      }			
		  }
	      
		else if (z1 > 0 && z2 > 0 && (z1 != z2) || (z1 < 0) && (z2 < 0) )
		  {
		    loclines.DeleteElement(i);
		    lindex.DeleteElement(i);
		    oldnl--;
		    i--;
		  }
	      }
	  




	    legalpoints.SetSize(plainpoints.Size());
	    for (int i = 1; i <= legalpoints.Size(); i++)
	      legalpoints.Elem(i) = 1;

	    double avy = 0;
	    for (int i = 1; i <= plainpoints.Size(); i++)
	      avy += plainpoints.Elem(i).Y();
	    avy *= 1./plainpoints.Size();
		

	    for (int i = 1; i <= plainpoints.Size(); i++)
	      {
		if (plainzones.Elem(i) < 0)
		  {
		    plainpoints.Elem(i) = Point2d (1e4, 1e4);
		    legalpoints.Elem(i) = 0;
		  }
		if (pindex.Elem(i) == -1)
		  {
		    legalpoints.Elem(i) = 0;
		  }
		    

		if (plainpoints.Elem(i).Y() < -1e-10*avy) // changed
		  {
		    legalpoints.Elem(i) = 0;
		  }
	      }
	    /*
	      for (i = 3; i <= plainpoints.Size(); i++)
	      if (sqr (plainpoints.Get(i).X()) + sqr (plainpoints.Get(i).Y())
	      > sqr (2 + 0.2 * qualclass))
	      legalpoints.Elem(i) = 0;
	    */  

	    /*	  
		 int clp = 0;
		 for (i = 1; i <= plainpoints.Size(); i++)
		 if (legalpoints.Get(i))
		 clp++;
		 (*testout) << "legalpts: " << clp << "/" << plainpoints.Size() << endl; 

		 // sort legal/illegal lines
		 int lastleg = 2;
		 int firstilleg = oldnl;

		 while (lastleg < firstilleg)
		 {
		 while (legalpoints.Get(loclines.Get(lastleg).I1()) &&
		 legalpoints.Get(loclines.Get(lastleg).I2()) &&
		 lastleg < firstilleg)
		 lastleg++;
		 while ( ( !legalpoints.Get(loclines.Get(firstilleg).I1()) ||
		 !legalpoints.Get(loclines.Get(firstilleg).I2())) &&
		 lastleg < firstilleg)
		 firstilleg--;
	      
		 if (lastleg < firstilleg)
		 {
		 swap (loclines.Elem(lastleg), loclines.Elem(firstilleg));
		 swap (lindex.Elem(lastleg), lindex.Elem(firstilleg));
		 }
		 }

		 (*testout) << "leglines " << lastleg << "/" << oldnl << endl;
	    */
	

	    GetChartBoundary (chartboundpoints, 
			      chartboundpoints3d,
			      chartboundlines, h);

	    oldnp = plainpoints.Size();

	    maxlegalpoint = locpoints.Size();
	    maxlegalline = loclines.Size();



	    if (mparam.checkchartboundary)
	      {
		for (int i = 1; i <= chartboundpoints.Size(); i++)
		  {
		    plainpoints.Append (chartboundpoints.Get(i));
		    locpoints.Append (chartboundpoints3d.Get(i));
		    legalpoints.Append (0);
		  }
	      

		for (int i = 1; i <= chartboundlines.Size(); i++)
		  {
		    INDEX_2 line (chartboundlines.Get(i).I1()+oldnp,
				  chartboundlines.Get(i).I2()+oldnp);
		    loclines.Append (line);
		    //	      (*testout) << "line: " << line.I1() << "-" << line.I2() << endl;
		  }
	      }

	    oldnl = loclines.Size();
	    oldnp = plainpoints.Size();
	  }


	/*
	  if (qualclass > 100)
	  {
	  multithread.drawing = 1;
	  glrender(1);
	  cout << "qualclass 100, nfl = " << adfront->GetNFL() << endl;
	  }
	*/

	if (found)
	  {
	    rulenr = ApplyRules (plainpoints, legalpoints, maxlegalpoint,
				 loclines, maxlegalline, locelements,
				 dellines, qualclass);
	    //	    (*testout) << "Rule Nr = " << rulenr << endl;
	    if (!rulenr)
	      {
		found = 0;
		if ( debugflag || debugparam.haltnosuccess )
		  PrintWarning ("no rule found");
	      }
	  }
      
	for (int i = 1; i <= locelements.Size() && found; i++)
	  {
	    const Element2d & el = locelements.Get(i);

	    for (int j = 1; j <= el.GetNP(); j++)
	      if (el.PNum(j) <= oldnp && pindex.Get(el.PNum(j)) == -1)
		{
		  found = 0;
		  PrintSysError ("meshing2, index missing");
		}
	  }


	if (found)
	  {
	    locpoints.SetSize (plainpoints.Size());
	    upgeominfo.SetSize(locpoints.Size());

	    for (int i = oldnp+1; i <= plainpoints.Size(); i++)
	      {
		int err =
		  TransformFromPlain (plainpoints.Elem(i), locpoints.Elem(i), 
				      upgeominfo.Elem(i), h);

		if (err)
		  {
		    found = 0;

		    if ( debugflag || debugparam.haltnosuccess )
		      PrintSysError ("meshing2, Backtransformation failed");

		    break;
		  }
	      }
	  }
	  

	//      for (i = 1; i <= oldnl; i++)
	//        adfront -> ResetClass (lindex[i]);


	/*
	  double violateminh;
	  if (qualclass <= 10)
	  violateminh = 3;
	  else
	  violateminh = 3 * qualclass;

	  if (uselocalh && found) //  && qualclass <= 10)
	  {
	  for (i = 1; i <= locelements.Size(); i++)
	  {
	  Point3d pmin = locpoints.Get(locelements.Get(i).PNum(1));
	  Point3d pmax = pmin;
	  for (j = 2; j <= 3; j++)
	  {
	  const Point3d & hp = 
	  locpoints.Get(locelements.Get(i).PNum(j));
	  pmin.SetToMin (hp);
	  pmax.SetToMax (hp);
	  }
	  double minh = mesh.GetMinH (pmin, pmax);
	  if (h > violateminh * minh)
	  {
	  found = 0;
	  loclines.SetSize (oldnl);
	  locpoints.SetSize (oldnp);
	  }
	  }
	  }
	*/


	if (found) 
	  {
	    double violateminh = 3 + 0.1 * sqr (qualclass);
	    double minh = 1e8;
	    double newedgemaxh = 0;
	    for (int i = oldnl+1; i <= loclines.Size(); i++)
	      {
		double eh = Dist (locpoints.Get(loclines.Get(i).I1()),
				  locpoints.Get(loclines.Get(i).I2()));

		// Markus (brute force method to avoid bad elements on geometries like \_/ )
		//if(eh > 4.*mesh.GetH(locpoints.Get(loclines.Get(i).I1()))) found = 0;
		//if(eh > 4.*mesh.GetH(locpoints.Get(loclines.Get(i).I2()))) found = 0;
		// Markus end

		if (eh > newedgemaxh)
		  newedgemaxh = eh;
	      }

	    for (int i = 1; i <= locelements.Size(); i++)
	      {
		Point3d pmin = locpoints.Get(locelements.Get(i).PNum(1));
		Point3d pmax = pmin;
		for (int j = 2; j <= locelements.Get(i).GetNP(); j++)
		  {
		    const Point3d & hp = 
		      locpoints.Get(locelements.Get(i).PNum(j));
		    pmin.SetToMin (hp);
		    pmax.SetToMax (hp);
		  }
		double eh = mesh.GetMinH (pmin, pmax);
		if (eh < minh)
		  minh = eh;
	      }

	    for (int i = 1; i <= locelements.Size(); i++)
	      for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		if (Dist2 (locpoints.Get(locelements.Get(i).PNum(j)), pmid) > hinner*hinner)
		  found = 0;

	    //	  cout << "violate = " << newedgemaxh / minh << endl;
	    static double maxviolate = 0;
	    if (newedgemaxh / minh > maxviolate)
	      {
		maxviolate = newedgemaxh / minh;
		//	      cout << "max minhviolate = " << maxviolate << endl;
	      }


	    if (newedgemaxh > violateminh * minh)
	      {
		found = 0;
		loclines.SetSize (oldnl);
		locpoints.SetSize (oldnp);

		if ( debugflag || debugparam.haltnosuccess )
		  PrintSysError ("meshing2, maxh too large");


	      }
	  }



	/*
	// test good ComputeLineGeoInfo
	if (found)
	{
	// is line on chart ?
	for (i = oldnl+1; i <= loclines.Size(); i++)
	{
	int gisize;
	void *geominfo;

	if (ComputeLineGeoInfo (locpoints.Get(loclines.Get(i).I1()),
	locpoints.Get(loclines.Get(i).I2()),
	gisize, geominfo))
	found = 0;
	}
	}
	*/


	// changed for OCC meshing
	if (found)
	  {
	    // take geominfo from dellines
	    // upgeominfo.SetSize(locpoints.Size());

	    /*
	      for (i = 1; i <= dellines.Size(); i++)
	      for (j = 1; j <= 2; j++)
	      {
	      upgeominfo.Elem(loclines.Get(dellines.Get(i)).I(j)) =
	      adfront -> GetLineGeomInfo (lindex.Get(dellines.Get(i)), j);
	      }
	    */


	    for (int i = 1; i <= locelements.Size(); i++)
	      for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		{
		  int pi = locelements.Get(i).PNum(j);
		  if (pi <= oldnp)
		    {
		    
		      if (ChooseChartPointGeomInfo (mpgeominfo.Get(pi), upgeominfo.Elem(pi)))
			{
			  // cannot select, compute new one
			  PrintWarning ("calc point geominfo instead of using");
			  if (ComputePointGeomInfo (locpoints.Get(pi), upgeominfo.Elem(pi)))
			    {
			      found = 0;
			      PrintSysError ("meshing2d, geominfo failed");
			    }
			}
		    }
		}

	    /*
	    // use upgeominfo from ProjectFromPlane
	    for (i = oldnp+1; i <= locpoints.Size(); i++)
	    {
	    if (ComputePointGeomInfo (locpoints.Get(i), upgeominfo.Elem(i)))
	    {
	    found = 0;
	    if ( debugflag || debugparam.haltnosuccess )
	    PrintSysError ("meshing2d, compute geominfo failed");
	    }
	    }
	    */
	  }


	if (found && mparam.checkoverlap)
	  {
	    // cout << "checkoverlap" << endl;
	    // test for overlaps
	  
	    Point3d hullmin(1e10, 1e10, 1e10);
	    Point3d hullmax(-1e10, -1e10, -1e10);
	  
	    for (int i = 1; i <= locelements.Size(); i++)
	      for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		{
		  const Point3d & p = locpoints.Get(locelements.Get(i).PNum(j));
		  hullmin.SetToMin (p);
		  hullmax.SetToMax (p);
		}
	    hullmin += Vec3d (-his, -his, -his);
	    hullmax += Vec3d ( his,  his,  his);

	    surfeltree.GetIntersecting (hullmin, hullmax, intersecttrias);

	    critpoints.SetSize (0);
	    for (int i = oldnp+1; i <= locpoints.Size(); i++)
	      critpoints.Append (locpoints.Get(i));

	    for (int i = 1; i <= locelements.Size(); i++)
	      {
		const Element2d & tri = locelements.Get(i);
		if (tri.GetNP() == 3)
		  {
		    const Point3d & tp1 = locpoints.Get(tri.PNum(1));
		    const Point3d & tp2 = locpoints.Get(tri.PNum(2));
		    const Point3d & tp3 = locpoints.Get(tri.PNum(3));
		  
		    Vec3d tv1 (tp1, tp2);
		    Vec3d tv2 (tp1, tp3);
		  
		    double lam1, lam2;
		    for (lam1 = 0.2; lam1 <= 0.8; lam1 += 0.2)
		      for (lam2 = 0.2; lam2 + lam1 <= 0.8; lam2 += 0.2)
			{
			  Point3d hp = tp1 + lam1 * tv1 + lam2 * tv2;
			  critpoints.Append (hp);
			}
		  }
		else if (tri.GetNP() == 4)
		  {
		    const Point3d & tp1 = locpoints.Get(tri.PNum(1));
		    const Point3d & tp2 = locpoints.Get(tri.PNum(2));
		    const Point3d & tp3 = locpoints.Get(tri.PNum(3));
		    const Point3d & tp4 = locpoints.Get(tri.PNum(4));
		  
		    double l1, l2;
		    for (l1 = 0.1; l1 <= 0.9; l1 += 0.1)
		      for (l2 = 0.1; l2 <= 0.9; l2 += 0.1)
			{
			  Point3d hp;
			  hp.X() = 
			    (1-l1)*(1-l2) * tp1.X() +
			    l1*(1-l2) * tp2.X() +
			    l1*l2 * tp3.X() +
			    (1-l1)*l2 * tp4.X();
			  hp.Y() = 
			    (1-l1)*(1-l2) * tp1.Y() +
			    l1*(1-l2) * tp2.Y() +
			    l1*l2 * tp3.Y() +
			    (1-l1)*l2 * tp4.Y();
			  hp.Z() = 
			    (1-l1)*(1-l2) * tp1.Z() +
			    l1*(1-l2) * tp2.Z() +
			    l1*l2 * tp3.Z() +
			    (1-l1)*l2 * tp4.Z();


			  critpoints.Append (hp);
			}
		  }
	      }
	    /*
	      for (i = oldnl+1; i <= loclines.Size(); i++)
	      {
	      Point3d hp = locpoints.Get(loclines.Get(i).I1());
	      Vec3d hv(hp, locpoints.Get(loclines.Get(i).I2()));
	      int ncp = 2;
	      for (j = 1; j <= ncp; j++)
	      critpoints.Append ( hp + (double(j)/(ncp+1)) * hv);
	      }
	    */


	    /*
	      for (i = oldnp+1; i <= locpoints.Size(); i++)
	      {
	      const Point3d & p = locpoints.Get(i);
	    */


	    for (int i = 1; i <= critpoints.Size(); i++)
	      {
		const Point3d & p = critpoints.Get(i);
		 

		/*
		  for (j = 1; j <= mesh.GetNSE(); j++)
		  {
		*/
		int jj;
		for (jj = 1; jj <= intersecttrias.Size(); jj++)
		  {
		    int j = intersecttrias.Get(jj);
		    const Element2d & el = mesh.SurfaceElement(j);
		  
		    int ntrig = (el.GetNP() == 3) ? 1 : 2;

		    int jl;
		    for (jl = 1; jl <= ntrig; jl++)
		      {
			Point3d tp1, tp2, tp3;

			if (jl == 1)
			  {
			    tp1 = mesh.Point(el.PNum(1));
			    tp2 = mesh.Point(el.PNum(2));
			    tp3 = mesh.Point(el.PNum(3));
			  }
			else
			  {
			    tp1 = mesh.Point(el.PNum(1));
			    tp2 = mesh.Point(el.PNum(3));
			    tp3 = mesh.Point(el.PNum(4));
			  }

			int onchart = 0;
			for (int k = 1; k <= el.GetNP(); k++)
			  if (BelongsToActiveChart (mesh.Point(el.PNum(k)),
						    el.GeomInfoPi(k)))
			    onchart = 1;
			if (!onchart)
			  continue;
		      
			Vec3d e1(tp1, tp2);
			Vec3d e2(tp1, tp3);
			Vec3d n = Cross (e1, e2);
			n /= n.Length();
			double lam1, lam2, lam3;
			lam3 = n * Vec3d (tp1, p);
			LocalCoordinates (e1, e2, Vec3d (tp1, p), lam1, lam2);
		      
			if (fabs (lam3) < 0.1 * hshould && 
			    lam1 > 0 && lam2 > 0 && (lam1 + lam2) < 1)
			  {
#ifdef DEVELOP
			    cout << "overlap" << endl;
			    (*testout) << "overlap:" << endl
				       << "tri = " << tp1 << "-" << tp2 << "-" << tp3 << endl
				       << "point = " << p << endl
				       << "lam1, 2 = " << lam1 << ", " << lam2 << endl
				       << "lam3 = " << lam3 << endl;
			  
			    //		      cout << "overlap !!!" << endl;
#endif
			    for (int k = 1; k <= 5; k++)
			      adfront -> IncrementClass (lindex.Get(1));

			    found = 0;
			  
			    if ( debugflag || debugparam.haltnosuccess )
			      PrintWarning ("overlapping");
			  
			  
			    if (debugparam.haltoverlap)
			      {
				debugflag = 1;
			      }
			  
			    /*
			      multithread.drawing = 1;
			      glrender(1);
			    */
			  }
		      }
		  }
	      }
	  }


	if (found)
	  {
	    // check, whether new front line already exists

	    for (int i = oldnl+1; i <= loclines.Size(); i++)
	      {
		int nlgpi1 = loclines.Get(i).I1();
		int nlgpi2 = loclines.Get(i).I2();
		if (nlgpi1 <= pindex.Size() && nlgpi2 <= pindex.Size())
		  {
		    nlgpi1 = adfront->GetGlobalIndex (pindex.Get(nlgpi1));
		    nlgpi2 = adfront->GetGlobalIndex (pindex.Get(nlgpi2));

		    int exval = adfront->ExistsLine (nlgpi1, nlgpi2);
		    if (exval)
		      {
			cout << "ERROR: new line exits, val = " << exval << endl;
			(*testout) << "ERROR: new line exits, val = " << exval << endl;
			found = 0;


			if (debugparam.haltexistingline)
			  debugflag = 1;

		      }
		  }
	      }
	  
	  }


	/*
	  if (found)
	  {
	  // check, whether new triangles insert edges twice
	  for (i = 1; i <= locelements.Size(); i++)
	  for (j = 1; j <= 3; j++)
	  {
	  int tpi1 = locelements.Get(i).PNumMod (j);
	  int tpi2 = locelements.Get(i).PNumMod (j+1);
	  if (tpi1 <= pindex.Size() && tpi2 <= pindex.Size())
	  {
	  tpi1 = adfront->GetGlobalIndex (pindex.Get(tpi1));
	  tpi2 = adfront->GetGlobalIndex (pindex.Get(tpi2));

	  if (doubleedge.Used (INDEX_2(tpi1, tpi2)))
	  {
	  if (debugparam.haltexistingline)
	  debugflag = 1;
	  cerr << "ERROR Insert edge "
	  << tpi1 << " - " << tpi2 << " twice !!!" << endl;
	  found = 0;
	  }
	  doubleedge.Set (INDEX_2(tpi1, tpi2), 1);
	  }
	  }
	  }
	*/


	if (found)
	  {
	    // everything is ok, perform mesh update

	    ruleused.Elem(rulenr)++;


	    pindex.SetSize(locpoints.Size());
	      
	    for (int i = oldnp+1; i <= locpoints.Size(); i++)
	      {
		globind = mesh.AddPoint (locpoints.Get(i));
		pindex.Elem(i) = adfront -> AddPoint (locpoints.Get(i), globind);
	      }
	      
	    for (int i = oldnl+1; i <= loclines.Size(); i++)
	      {
		/*
		  for (j = 1; j <= locpoints.Size(); j++)
		  {
		  (*testout) << j << ": " << locpoints.Get(j) << endl;
		  }
		*/
	      
		/*
		  ComputeLineGeoInfo (locpoints.Get(loclines.Get(i).I1()),
		  locpoints.Get(loclines.Get(i).I2()),
		  gisize, geominfo);
		*/		  

		if (pindex.Get(loclines.Get(i).I1()) == -1 || 
		    pindex.Get(loclines.Get(i).I2()) == -1)
		  {
		    (*testout) << "pindex is 0" << endl;
		  }

		if (!upgeominfo.Get(loclines.Get(i).I1()).trignum || 
		    !upgeominfo.Get(loclines.Get(i).I2()).trignum)
		  {
		    cout << "new el: illegal geominfo" << endl;
		  }

		adfront -> AddLine (pindex.Get(loclines.Get(i).I1()),
				    pindex.Get(loclines.Get(i).I2()),
				    upgeominfo.Get(loclines.Get(i).I1()),
				    upgeominfo.Get(loclines.Get(i).I2()));
	      }
	    for (int i = 1; i <= locelements.Size(); i++)
	      {
		Element2d mtri(locelements.Get(i).GetNP());
		mtri = locelements.Get(i);
		mtri.SetIndex (facenr);


		// compute triangle geominfo:
		//	      (*testout) << "triggeominfo: ";
		for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		  {
		    mtri.GeomInfoPi(j) = upgeominfo.Get(locelements.Get(i).PNum(j));
		    //		  (*testout) << mtri.GeomInfoPi(j).trignum << " ";
		  }
		//	      (*testout) << endl;

		for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		  {
		    mtri.PNum(j) = 
		      locelements.Elem(i).PNum(j) =
		      adfront -> GetGlobalIndex (pindex.Get(locelements.Get(i).PNum(j)));
		  }
	      
		
	      
	      
		mesh.AddSurfaceElement (mtri);
		cntelem++;
		//	      cout << "elements: " << cntelem << endl;


	      
		Box<3> box;
		box.Set (mesh[mtri[0]]);
		box.Add (mesh[mtri[1]]);
		box.Add (mesh[mtri[2]]);
		surfeltree.Insert (box, mesh.GetNSE());

		const Point3d & sep1 = mesh.Point (mtri.PNum(1));
		const Point3d & sep2 = mesh.Point (mtri.PNum(2));
		const Point3d & sep3 = mesh.Point (mtri.PNum(3));

		double trigarea = Cross (Vec3d (sep1, sep2), 
					 Vec3d (sep1, sep3)).Length() / 2;

		if (mtri.GetNP() == 4)
		  {
		    const Point3d & sep4 = mesh.Point (mtri.PNum(4));
		    trigarea += Cross (Vec3d (sep1, sep3), 
				       Vec3d (sep1, sep4)).Length() / 2;
		  }

		meshedarea += trigarea;

		if(maxarea > 0 && meshedarea-meshedarea_before > maxarea)
		  {
		    cerr << "meshed area = " << meshedarea-meshedarea_before << endl
			 << "maximal area = " << maxarea << endl
			 << "GIVING UP" << endl;
		    return MESHING2_GIVEUP;
		  }
	      


		for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		  {
		    int gpi = locelements.Get(i).PNum(j);

		    int oldts = trigsonnode.Size();
		    if (gpi >= oldts+PointIndex::BASE)
		      {
			trigsonnode.SetSize (gpi+1-PointIndex::BASE);
			illegalpoint.SetSize (gpi+1-PointIndex::BASE);
			for (int k = oldts+PointIndex::BASE; 
			     k <= gpi; k++)
			  {
			    trigsonnode[k] = 0;
			    illegalpoint[k] = 0;
			  }
		      }

		    trigsonnode[gpi]++;
		  
		    if (trigsonnode[gpi] > 20)
		      {
			illegalpoint[gpi] = 1;
			//		      cout << "illegal point: " << gpi << endl;
			(*testout) << "illegal point: " << gpi << endl;
		      }

		    static int mtonnode = 0;
		    if (trigsonnode[gpi] > mtonnode)
		      mtonnode = trigsonnode[gpi];
		  }
		//	      cout << "els = " << cntelem << " trials = " << trials << endl;
		//	      if (trials > 100)		return;
	      }
	      
	    for (int i = 1; i <= dellines.Size(); i++)
	      adfront -> DeleteLine (lindex.Get(dellines.Get(i)));
	      
	    //	  rname = rules.Get(rulenr)->Name();
#ifdef MYGRAPH
	    if (silentflag<3) 
	      {
		plotsurf.DrawPnL(locpoints, loclines);
		plotsurf.Plot(testmode, testmode);
	      }
#endif

	    if (morerisc)
	      {
		cout << "generated due to morerisc" << endl;
		//	      multithread.drawing = 1;
		//	      glrender(1);
	      }



	  
	    if ( debugparam.haltsuccess || debugflag )
	      {
		// adfront -> PrintOpenSegments (*testout);
		cout << "success of rule" << rules.Get(rulenr)->Name() << endl;
		multithread.drawing = 1;
		multithread.testmode = 1;
		multithread.pause = 1;


		/*
		  extern STLGeometry * stlgeometry;
		  stlgeometry->ClearMarkedSegs();
		  for (i = 1; i <= loclines.Size(); i++)
		  {
		  stlgeometry->AddMarkedSeg(locpoints.Get(loclines.Get(i).I1()),
		  locpoints.Get(loclines.Get(i).I2()));
		  }
		*/

		(*testout) << "success of rule" << rules.Get(rulenr)->Name() << endl;
		(*testout) << "trials = " << trials << endl;

		(*testout) << "locpoints " << endl;
		for (int i = 1; i <= pindex.Size(); i++)
		  (*testout) << adfront->GetGlobalIndex (pindex.Get(i)) << endl;

		(*testout) << "old number of lines = " << oldnl << endl;
		for (int i = 1; i <= loclines.Size(); i++)
		  {
		    (*testout) << "line ";
		    for (int j = 1; j <= 2; j++)
		      {
			int hi = 0;
			if (loclines.Get(i).I(j) >= 1 &&
			    loclines.Get(i).I(j) <= pindex.Size())
			  hi = adfront->GetGlobalIndex (pindex.Get(loclines.Get(i).I(j)));

			(*testout) << hi << " ";
		      }
		    (*testout) << " : " 
			       << plainpoints.Get(loclines.Get(i).I1()) << " - "
			       << plainpoints.Get(loclines.Get(i).I2()) << " 3d: "
			       << locpoints.Get(loclines.Get(i).I1()) << " - "
			       << locpoints.Get(loclines.Get(i).I2()) 
			       << endl;
		  }



		glrender(1);
	      }
	  }
	else
	  {
	    adfront -> IncrementClass (lindex.Get(1));

	    if ( debugparam.haltnosuccess || debugflag )
	      {
		cout << "Problem with seg " << gpi1 << " - " << gpi2
		     << ", class = " << qualclass << endl;

		(*testout) << "Problem with seg " << gpi1 << " - " << gpi2
			   << ", class = " << qualclass << endl;

		multithread.drawing = 1;
		multithread.testmode = 1;
		multithread.pause = 1;


		/*
		  extern STLGeometry * stlgeometry;
		  stlgeometry->ClearMarkedSegs();
		  for (i = 1; i <= loclines.Size(); i++)
		  {
		  stlgeometry->AddMarkedSeg(locpoints.Get(loclines.Get(i).I1()),
		  locpoints.Get(loclines.Get(i).I2()));
		  }
		*/

		for (int i = 1; i <= loclines.Size(); i++)
		  {
		    (*testout) << "line ";
		    for (int j = 1; j <= 2; j++)
		      {
			int hi = 0;
			if (loclines.Get(i).I(j) >= 1 &&
			    loclines.Get(i).I(j) <= pindex.Size())
			  hi = adfront->GetGlobalIndex (pindex.Get(loclines.Get(i).I(j)));

			(*testout) << hi << " ";
		      }
		    (*testout) << " : " 
			       << plainpoints.Get(loclines.Get(i).I1()) << " - "
			       << plainpoints.Get(loclines.Get(i).I2()) << " 3d: "
			       << locpoints.Get(loclines.Get(i).I1()) << " - "
			       << locpoints.Get(loclines.Get(i).I2()) 
			       << endl;
		  }


		/*
		  cout << "p1gi = " << blgeominfo[0].trignum 
		  << ", p2gi = " << blgeominfo[1].trignum << endl;
		*/

		glrender(1);
	      }

	  
#ifdef MYGRAPH      
	    if (silentflag<3)
	      {
		if (testmode || trials%2 == 0)
		  {
		    plotsurf.DrawPnL(locpoints, loclines);
		    plotsurf.Plot(testmode, testmode);
		  }
	      }
#endif
	  }

      }

    PrintMessage (3, "Surface meshing done");

    adfront->PrintOpenSegments (*testout);

    multithread.task = savetask;


    //  cout << "surfeltree.depth = " << surfeltree.Tree().Depth() << endl;
    EndMesh ();

    if (!adfront->Empty())
      return MESHING2_GIVEUP;
    
    return MESHING2_OK;
  }









}







#ifdef OPENGL

/* *********************** Draw Surface Meshing **************** */


#include <visual.hpp>
#include <stlgeom.hpp>

namespace netgen 
{

  extern STLGeometry * stlgeometry;
  extern Mesh * mesh;
  VisualSceneSurfaceMeshing vssurfacemeshing;



  void glrender (int wait)
  {
    //  cout << "plot adfront" << endl;

    if (multithread.drawing)
      {
	//      vssurfacemeshing.Render();
	Render ();
      
	if (wait || multithread.testmode)
	  {
	    multithread.pause = 1;
	  }
	while (multithread.pause);
      }
  }



  VisualSceneSurfaceMeshing :: VisualSceneSurfaceMeshing ()
    : VisualScene()
  {
    ;
  }

  VisualSceneSurfaceMeshing :: ~VisualSceneSurfaceMeshing ()
  {
    ;
  }

  void VisualSceneSurfaceMeshing :: DrawScene ()
  {
    int i, j, k;

    if (loclines.Size() != changeval)
      {
	center = Point<3>(0,0,-5);
	rad = 0.1;
  
	CalcTransformationMatrices();
	changeval = loclines.Size();
      }

  glClearColor(backcolor, backcolor, backcolor, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  SetLight();

  //  glEnable (GL_COLOR_MATERIAL);

  //  glDisable (GL_SHADING);
  //  glColor3f (0.0f, 1.0f, 1.0f);
  //  glLineWidth (1.0f);
  //  glShadeModel (GL_SMOOTH);

  //  glCallList (linelists.Get(1));

  //  SetLight();

  glPushMatrix();
  glMultMatrixf (transformationmat);

  glShadeModel (GL_SMOOTH);
  // glDisable (GL_COLOR_MATERIAL);
  glEnable (GL_COLOR_MATERIAL);
  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  //  glEnable (GL_LIGHTING);

  double shine = vispar.shininess;
  double transp = vispar.transp;

  glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, shine);
  glLogicOp (GL_COPY);



  /*

  float mat_col[] = { 0.2, 0.2, 0.8, 1 };
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col);

  glPolygonOffset (1, 1);
  glEnable (GL_POLYGON_OFFSET_FILL);

    float mat_colbl[] = { 0.8, 0.2, 0.2, 1 };
    float mat_cololdl[] = { 0.2, 0.8, 0.2, 1 };
    float mat_colnewl[] = { 0.8, 0.8, 0.2, 1 };


    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    glPolygonOffset (1, -1);
    glLineWidth (3);

    for (i = 1; i <= loclines.Size(); i++)
      {
	if (i == 1)
	  {
	    glEnable (GL_POLYGON_OFFSET_FILL);
	    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colbl);
	  }
	else if (i <= oldnl) 
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_cololdl);
	else 
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colnewl);

	int pi1 = loclines.Get(i).I1();
	int pi2 = loclines.Get(i).I2();

	if (pi1 >= 1 && pi2 >= 1)
	  {
	    Point3d p1 = locpoints.Get(pi1);
	    Point3d p2 = locpoints.Get(pi2);
	  
	    glBegin (GL_LINES);
	    glVertex3f (p1.X(), p1.Y(), p1.Z());
	    glVertex3f (p2.X(), p2.Y(), p2.Z());
	    glEnd();
	  }

	glDisable (GL_POLYGON_OFFSET_FILL);
      }
  

    glLineWidth (1);


    glPointSize (5);
    float mat_colp[] = { 1, 0, 0, 1 };
    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colp);
    glBegin (GL_POINTS);
    for (i = 1; i <= locpoints.Size(); i++)
      {
	Point3d p = locpoints.Get(i);
	glVertex3f (p.X(), p.Y(), p.Z());
      }
    glEnd();


    glPopMatrix();
  */

    float mat_colp[] = { 1, 0, 0, 1 };

    float mat_col2d1[] = { 1, 0.5, 0.5, 1 };
    float mat_col2d[] = { 1, 1, 1, 1 };
    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col2d);
  
    double scalex = 0.1, scaley = 0.1;

    glBegin (GL_LINES);
    for (i = 1; i <= loclines.Size(); i++)
      {
	glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col2d);
	if (i == 1)
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col2d1);

	int pi1 = loclines.Get(i).I1();
	int pi2 = loclines.Get(i).I2();

	if (pi1 >= 1 && pi2 >= 1)
	  {
	    Point2d p1 = plainpoints.Get(pi1);
	    Point2d p2 = plainpoints.Get(pi2);
	  
	    glBegin (GL_LINES);
	    glVertex3f (scalex * p1.X(), scaley * p1.Y(), -5);
	    glVertex3f (scalex * p2.X(), scaley * p2.Y(), -5);
	    glEnd();
	  }
      }
    glEnd ();


    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colp);
    glBegin (GL_POINTS);
    for (i = 1; i <= plainpoints.Size(); i++)
      {
	Point2d p = plainpoints.Get(i);
	glVertex3f (scalex * p.X(), scaley * p.Y(), -5);
      }
    glEnd();






  glDisable (GL_POLYGON_OFFSET_FILL);
 
  glPopMatrix();
  DrawCoordinateCross ();
  DrawNetgenLogo ();
  glFinish();  

  /*
    glDisable (GL_POLYGON_OFFSET_FILL);

    //  cout << "draw surfacemeshing" << endl;
    //
    //  if (changeval != stlgeometry->GetNT())
    //      BuildScene();
    //      changeval = stlgeometry->GetNT();
    

    glClearColor(backcolor, backcolor, backcolor, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    SetLight();

    glPushMatrix();
    glLoadMatrixf (transmat);
    glMultMatrixf (rotmat);

    glShadeModel (GL_SMOOTH);
    glDisable (GL_COLOR_MATERIAL);
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    float mat_spec_col[] = { 1, 1, 1, 1 };
    glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, mat_spec_col);

    double shine = vispar.shininess;
    double transp = vispar.transp;

    glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, shine);
    glLogicOp (GL_COPY);


    float mat_col[] = { 0.2, 0.2, 0.8, transp };
    float mat_colrt[] = { 0.2, 0.8, 0.8, transp };
    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col);

    glPolygonOffset (1, 1);
    glEnable (GL_POLYGON_OFFSET_FILL);

    glColor3f (1.0f, 1.0f, 1.0f);

    glEnable (GL_NORMALIZE);
    
    //  glBegin (GL_TRIANGLES);
    //      for (j = 1; j <= stlgeometry -> GetNT(); j++)
    //      {
    //      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col);
    //      if (j == geomtrig)
    //      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colrt);
	

    //      const STLReadTriangle & tria = stlgeometry -> GetReadTriangle(j);
    //      glNormal3f (tria.normal.X(),
    //      tria.normal.Y(),
    //      tria.normal.Z());
		  
    //      for (k = 0; k < 3; k++)
    //      {
    //      glVertex3f (tria.pts[k].X(),
    //      tria.pts[k].Y(),
    //      tria.pts[k].Z());
    //      }
    //      }    
    //      glEnd ();
    


    glDisable (GL_POLYGON_OFFSET_FILL);

    float mat_colbl[] = { 0.8, 0.2, 0.2, 1 };
    float mat_cololdl[] = { 0.2, 0.8, 0.2, 1 };
    float mat_colnewl[] = { 0.8, 0.8, 0.2, 1 };


    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    glPolygonOffset (1, -1);
    glLineWidth (3);

    for (i = 1; i <= loclines.Size(); i++)
      {
	if (i == 1)
	  {
	    glEnable (GL_POLYGON_OFFSET_FILL);
	    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colbl);
	  }
	else if (i <= oldnl) 
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_cololdl);
	else 
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colnewl);

	int pi1 = loclines.Get(i).I1();
	int pi2 = loclines.Get(i).I2();

	if (pi1 >= 1 && pi2 >= 1)
	  {
	    Point3d p1 = locpoints.Get(pi1);
	    Point3d p2 = locpoints.Get(pi2);
	  
	    glBegin (GL_LINES);
	    glVertex3f (p1.X(), p1.Y(), p1.Z());
	    glVertex3f (p2.X(), p2.Y(), p2.Z());
	    glEnd();
	  }

	glDisable (GL_POLYGON_OFFSET_FILL);
      }


    glLineWidth (1);


    glPointSize (5);
    float mat_colp[] = { 1, 0, 0, 1 };
    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colp);
    glBegin (GL_POINTS);
    for (i = 1; i <= locpoints.Size(); i++)
      {
	Point3d p = locpoints.Get(i);
	glVertex3f (p.X(), p.Y(), p.Z());
      }
    glEnd();


    glPopMatrix();


    float mat_col2d1[] = { 1, 0.5, 0.5, 1 };
    float mat_col2d[] = { 1, 1, 1, 1 };
    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col2d);
  
    double scalex = 0.1, scaley = 0.1;

    glBegin (GL_LINES);
    for (i = 1; i <= loclines.Size(); i++)
      {
	glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col2d);
	if (i == 1)
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col2d1);

	int pi1 = loclines.Get(i).I1();
	int pi2 = loclines.Get(i).I2();

	if (pi1 >= 1 && pi2 >= 1)
	  {
	    Point2d p1 = plainpoints.Get(pi1);
	    Point2d p2 = plainpoints.Get(pi2);
	  
	    glBegin (GL_LINES);
	    glVertex3f (scalex * p1.X(), scaley * p1.Y(), -5);
	    glVertex3f (scalex * p2.X(), scaley * p2.Y(), -5);
	    glEnd();
	  }
      }
    glEnd ();


    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colp);
    glBegin (GL_POINTS);
    for (i = 1; i <= plainpoints.Size(); i++)
      {
	Point2d p = plainpoints.Get(i);
	glVertex3f (scalex * p.X(), scaley * p.Y(), -5);
      }
    glEnd();

    glFinish();  
*/
  }


  void VisualSceneSurfaceMeshing :: BuildScene (int zoomall)
  {
    int i, j, k;
    /*
      center = stlgeometry -> GetBoundingBox().Center();
      rad = stlgeometry -> GetBoundingBox().Diam() / 2;

      CalcTransformationMatrices();
    */
  }

}


#else
namespace netgen
{
  void glrender (int wait)
  { ; }
}
#endif
