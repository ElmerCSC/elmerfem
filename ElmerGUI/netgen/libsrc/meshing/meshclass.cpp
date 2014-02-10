#include <mystdlib.h>
#include "meshing.hpp"

#ifdef PARALLEL
#include <parallel.hpp>
#endif

namespace netgen
{

  Mesh :: Mesh ()
  {
    volelements.SetName ("vol elements");
    surfelements.SetName ("surf elements");
    points.SetName ("meshpoints");

    boundaryedges = NULL;
    surfelementht = NULL; 
    segmentht = NULL;

    lochfunc = NULL;
    mglevels = 1;
    elementsearchtree = NULL;
    elementsearchtreets = NextTimeStamp();
    majortimestamp = timestamp = NextTimeStamp();
    hglob = 1e10;
    hmin = 0;
    numvertices = -1;
    dimension = 3;

    topology = new MeshTopology (*this);
    curvedelems = new CurvedElements (*this);
    clusters = new AnisotropicClusters (*this);
    ident = new Identifications (*this);

    hpelements = NULL;
    coarsemesh = NULL;

    ps_startelement = 0;

    geomtype = NO_GEOM;

    bcnames.SetSize(0);
#ifdef PARALLEL
    paralleltop = new ParallelMeshTopology (*this);
#endif
  }


  Mesh :: ~Mesh()
  {
    delete lochfunc;
    delete boundaryedges;
    delete surfelementht;
    delete segmentht;
    delete curvedelems;
    delete clusters;
    delete topology;
    delete ident;
    delete elementsearchtree;

    delete coarsemesh;
    delete hpelements;
    
    for (int i = 0; i < materials.Size(); i++)
      delete [] materials[i];

    for(int i = 0; i < userdata_int.Size(); i++)
      delete userdata_int[i];
    for(int i = 0; i < userdata_double.Size(); i++)
      delete userdata_double[i];

    for (int i = 0; i < bcnames.Size(); i++ )
      if ( bcnames[i] ) delete bcnames[i];

#ifdef PARALLEL
    delete paralleltop;
#endif
  }


  Mesh & Mesh :: operator= (const Mesh & mesh2)
  {
    points = mesh2.points;
    // eltyps = mesh2.eltyps;
    segments = mesh2.segments;
    surfelements = mesh2.surfelements;
    volelements = mesh2.volelements;
    lockedpoints = mesh2.lockedpoints;
    facedecoding = mesh2.facedecoding;
    dimension = mesh2.dimension;

    bcnames.SetSize( mesh2.bcnames.Size() );
    for ( int i = 0; i < mesh2.bcnames.Size(); i++ )
      if ( mesh2.bcnames[i] ) bcnames[i] = new string ( *mesh2.bcnames[i] );
      else bcnames[i] = 0;

    return *this;
  }


  void Mesh :: DeleteMesh()
  {
    points.SetSize(0);
    segments.SetSize(0);
    surfelements.SetSize(0);
    volelements.SetSize(0);
    lockedpoints.SetSize(0);
    surfacesonnode.SetSize(0);

    delete boundaryedges;
    boundaryedges = NULL;

    openelements.SetSize(0);
    facedecoding.SetSize(0);

    delete ident;
    ident = new Identifications (*this);
    delete topology;
    topology = new MeshTopology (*this);
    delete curvedelems;
    curvedelems = new CurvedElements (*this);
    delete clusters;
    clusters = new AnisotropicClusters (*this);

    for ( int i = 0; i < bcnames.Size(); i++ )
      if ( bcnames[i] ) delete bcnames[i];

#ifdef PARALLEL
    delete paralleltop;
    paralleltop = new ParallelMeshTopology (*this);
#endif


    timestamp = NextTimeStamp();
  }



  PointIndex Mesh :: AddPoint (const Point3d & p, int layer)
  { 
    NgLock lock(mutex);
    lock.Lock();

    timestamp = NextTimeStamp();

    PointIndex pi = points.Size() + PointIndex::BASE;
    points.Append ( MeshPoint (p, layer, INNERPOINT) ); 

#ifdef PARALLEL
    points.Last().SetGhost(0);
#endif

    lock.UnLock();

    return pi;
  }

  PointIndex Mesh :: AddPoint (const Point3d & p, int layer, POINTTYPE type)
  { 
    NgLock lock(mutex);
    lock.Lock();

    timestamp = NextTimeStamp();

    PointIndex pi = points.Size() + PointIndex::BASE;
    points.Append ( MeshPoint (p, layer, type) ); 

#ifdef PARALLEL
    points.Last().SetGhost(0);
#endif

    lock.UnLock();

    return pi;
  }


#ifdef PARALLEL
  PointIndex Mesh :: AddPoint (const Point3d & p, bool isghost,  int layer)
  { 
    NgLock lock(mutex);
    lock.Lock();

    timestamp = NextTimeStamp();

    PointIndex pi = points.Size() + PointIndex::BASE;
    points.Append ( MeshPoint (p, layer, INNERPOINT) ); 

    points.Last().SetGhost(isghost);

    lock.UnLock();

    return pi;
  }

  PointIndex Mesh :: AddPoint (const Point3d & p, bool isghost, int layer, POINTTYPE type)
  { 
    NgLock lock(mutex);
    lock.Lock();

    timestamp = NextTimeStamp();

    PointIndex pi = points.Size() + PointIndex::BASE;
    points.Append ( MeshPoint (p, layer, type) ); 

    points.Last().SetGhost(isghost);

    lock.UnLock();

    return pi;
  }

#endif



  SegmentIndex Mesh :: AddSegment (const Segment & s)
  { 
    NgLock lock(mutex);	
    lock.Lock();
    timestamp = NextTimeStamp();

    int maxn = max2 (s.p1, s.p2);
    maxn += 1-PointIndex::BASE;

    /*
      if (maxn > ptyps.Size())
      {
      int maxo = ptyps.Size();
      ptyps.SetSize (maxn);
      for (int i = maxo; i < maxn; i++)
      ptyps[i] = INNERPOINT;
      }

      if (ptyps[s.p1] > EDGEPOINT) ptyps[s.p1] = EDGEPOINT;
      if (ptyps[s.p2] > EDGEPOINT) ptyps[s.p2] = EDGEPOINT;
    */

    if (maxn <= points.Size())
      {
	if (points[s.p1].Type() > EDGEPOINT)
	  points[s.p1].SetType (EDGEPOINT);
	if (points[s.p2].Type() > EDGEPOINT)
	  points[s.p2].SetType (EDGEPOINT);
      }
    /*
      else
      {
      cerr << "edge points nrs > points.Size" << endl;
      }
    */

    SegmentIndex si = segments.Size();
    segments.Append (s); 
  
    lock.UnLock();
    return si;
  }

  SurfaceElementIndex Mesh :: AddSurfaceElement (const Element2d & el)
  {     
    NgLock lock(mutex);
    lock.Lock();
    timestamp = NextTimeStamp();

    int maxn = el[0];
    for (int i = 1; i < el.GetNP(); i++)
      if (el[i] > maxn) maxn = el[i];

    maxn += 1-PointIndex::BASE;

    /*
      if (maxn > ptyps.Size())
      {
      int maxo = ptyps.Size();
      ptyps.SetSize (maxn);
      for (i = maxo+PointIndex::BASE; 
      i < maxn+PointIndex::BASE; i++)
      ptyps[i] = INNERPOINT;
      
      }
    */
    if (maxn <= points.Size())
      {
	for (int i = 0; i < el.GetNP(); i++)
	  if (points[el[i]].Type() > SURFACEPOINT)
	    points[el[i]].SetType(SURFACEPOINT);
      }
    /*
      else
      {
      cerr << "surf points nrs > points.Size" << endl;      
      }
    */

    SurfaceElementIndex si = surfelements.Size();
    surfelements.Append (el); 

    if (el.index > facedecoding.Size())
      cerr << "has no facedecoding: fd.size = " << facedecoding.Size() << ", ind = " << el.index << endl;

    surfelements.Last().next = facedecoding[el.index-1].firstelement;
    facedecoding[el.index-1].firstelement = si;

#ifdef PARALLEL
    surfelements.Last().SetGhost ( el.IsGhost() );
#endif

    lock.UnLock();
    return si;
  }


  ElementIndex Mesh :: AddVolumeElement (const Element & el)
  { 
    NgLock lock(mutex);
    lock.Lock();

    int maxn = el[0];
    for (int i = 1; i < el.GetNP(); i++)
      if (el[i] > maxn) maxn = el[i];

    maxn += 1-PointIndex::BASE;

    /*
      if (maxn > ptyps.Size())
      {
      int maxo = ptyps.Size();
      ptyps.SetSize (maxn);
      for (i = maxo+PointIndex::BASE; 
      i < maxn+PointIndex::BASE; i++)
      ptyps[i] = INNERPOINT;
      }
    */
    /*
      if (maxn > points.Size())
      {
      cerr << "add vol element before point" << endl;
      }
    */

    int ve = volelements.Size();

    volelements.Append (el); 
    volelements.Last().flags.illegal_valid = 0;

#ifdef PARALLEL
    volelements.Last().SetGhost ( el.IsGhost() );
#endif

    // while (volelements.Size() > eltyps.Size())
    // eltyps.Append (FREEELEMENT);
  
    timestamp = NextTimeStamp();

    lock.UnLock();
    return ve;
  }






  void Mesh :: Save (const string & filename) const
  {
    
    ofstream outfile(filename.c_str());

    Save(outfile);
  }



  void Mesh :: Save (ostream & outfile) const
  {
    int i, j;

    double scale = 1;  // globflags.GetNumFlag ("scale", 1);
    int inverttets = 0;  // globflags.GetDefineFlag ("inverttets");
    int invertsurf = 0;  // globflags.GetDefineFlag ("invertsurfacemesh");



    outfile << "mesh3d" << "\n";

    outfile << "dimension\n" << GetDimension() << "\n";

    outfile << "geomtype\n" << int(geomtype) << "\n";


    outfile << "\n";
    outfile << "# surfnr    bcnr   domin  domout      np      p1      p2      p3"
	    << "\n";

    
    switch (geomtype)
      {
      case GEOM_STL:
        outfile << "surfaceelementsgi" << "\n";
        break;
      case GEOM_OCC: case GEOM_ACIS:
        outfile << "surfaceelementsuv" << "\n";
        break;
      default:
        outfile << "surfaceelements" << "\n";
      }

    outfile << GetNSE() << "\n";

    SurfaceElementIndex sei;
    for (sei = 0; sei < GetNSE(); sei++)
      {
	if ((*this)[sei].GetIndex())
	  {
	    outfile.width(8);
	    outfile << GetFaceDescriptor((*this)[sei].GetIndex ()).SurfNr()+1;
	    outfile.width(8);
	    outfile << GetFaceDescriptor((*this)[sei].GetIndex ()).BCProperty();
	    outfile.width(8);	  
	    outfile << GetFaceDescriptor((*this)[sei].GetIndex ()).DomainIn();
	    outfile.width(8);	  
	    outfile << GetFaceDescriptor((*this)[sei].GetIndex ()).DomainOut();
	  }
	else
	  outfile << "       0       0       0";

	Element2d sel = (*this)[sei];
	if (invertsurf)
	  sel.Invert();

	outfile.width(8);
	outfile << sel.GetNP();

	for (j = 0; j < sel.GetNP(); j++)
	  {
	    outfile.width(8);	  
	    outfile << sel[j];
	  }


        switch (geomtype)
          {
          case GEOM_STL:
            for (j = 1; j <= sel.GetNP(); j++)
              {
                outfile.width(7);	  
                outfile << " " << sel.GeomInfoPi(j).trignum;
              }
            break;
          case GEOM_OCC: case GEOM_ACIS:
            for (j = 1; j <= sel.GetNP(); j++)
              {
                outfile.width(7);	  
                outfile << " " << sel.GeomInfoPi(j).u;
                outfile << " " << sel.GeomInfoPi(j).v;
              }
            break;
          default:
            outfile << "\n";
          }


	outfile << endl;
      }

    outfile << "\n" << "\n";
    outfile << "#  matnr      np      p1      p2      p3      p4" << "\n";
    outfile << "volumeelements" << "\n";
    outfile << GetNE() << "\n";

    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      {
	outfile.width(8);
	outfile << (*this)[ei].GetIndex();
	outfile.width(8);
	outfile << (*this)[ei].GetNP();

	Element el = (*this)[ei];
	if (inverttets)
	  el.Invert();

	/*
	  for (j = 0; j < el.GetNP(); j++)
	  for (int k = 0; k < el.GetNP()-1; k++)
	  if (el[k] > el[k+1]) swap (el[k], el[k+1]);
	*/

	for (j = 0; j < el.GetNP(); j++)
	  {
	    outfile.width(8);
	    outfile << el[j];
	  }
	outfile << "\n";
      }


    outfile << "\n" << "\n";
//     outfile << "   surf1   surf2      p1      p2" << "\n";
    outfile << "# surfid  0   p1   p2   trignum1    trignum2   domin/surfnr1    domout/surfnr2   ednr1   dist1   ednr2   dist2 \n";
    outfile << "edgesegmentsgi2" << "\n";
    outfile << GetNSeg() << "\n";

    for (i = 1; i <= GetNSeg(); i++)
      {
	const Segment & seg = LineSegment (i);
	outfile.width(8);
	outfile << seg.si; // 2D: bc number, 3D: wievielte Kante
	outfile.width(8);
	outfile << 0;
	outfile.width(8);
	outfile << seg.p1;
	outfile.width(8);
	outfile << seg.p2;
	outfile << " ";
	outfile.width(8);
	outfile << seg.geominfo[0].trignum;  // stl dreiecke
	outfile << " ";
	outfile.width(8);
	outfile << seg.geominfo[1].trignum; // << endl;  // stl dreieck
	
	if (dimension == 3)
	  {
	    outfile << " ";
	    outfile.width(8);
	    outfile << seg.surfnr1+1;
	    outfile << " ";
	    outfile.width(8);
	    outfile << seg.surfnr2+1;
	  }
	else
	  {
	    outfile << " ";
	    outfile.width(8);
	    outfile << seg.domin;
	    outfile << " ";
	    outfile.width(8);
	    outfile << seg.domout;
	  }

	outfile << " ";
	outfile.width(8);
	outfile << seg.edgenr;
	outfile << " ";
	outfile.width(12);
	outfile.precision(16);
	outfile << seg.epgeominfo[0].dist;  // splineparameter (2D)
	outfile << " ";
	outfile.width(8);
	outfile.precision(16);
	outfile << seg.epgeominfo[1].edgenr;  // geometry dependent
	outfile << " ";
	outfile.width(12);
	outfile << seg.epgeominfo[1].dist;

	outfile << "\n";
      }


    outfile << "\n" << "\n";
    outfile << "#          X             Y             Z" << "\n";
    outfile << "points" << "\n";
    outfile << GetNP() << "\n";
    outfile.precision(16);
    outfile.setf (ios::fixed, ios::floatfield);
    outfile.setf (ios::showpoint);

    PointIndex pi;
    for (pi = PointIndex::BASE; 
	 pi < GetNP()+PointIndex::BASE; pi++)
      {
	outfile.width(22);
	outfile << (*this)[pi](0)/scale << "  ";
	outfile.width(22);
	outfile << (*this)[pi](1)/scale << "  ";
	outfile.width(22);
	outfile << (*this)[pi](2)/scale << "\n";
      }      

    if (ident -> GetMaxNr() > 0)
      {
	outfile << "identifications\n";
	ARRAY<INDEX_2> identpairs;
	int cnt = 0;
	for (i = 1; i <= ident -> GetMaxNr(); i++)
	  {
	    ident -> GetPairs (i, identpairs);
	    cnt += identpairs.Size();
	  }
	outfile << cnt << "\n";
	for (i = 1; i <= ident -> GetMaxNr(); i++)
	  {
	    ident -> GetPairs (i, identpairs);
	    for (j = 1; j <= identpairs.Size(); j++)
	      {
		outfile.width (8);
		outfile << identpairs.Get(j).I1();
		outfile.width (8);
		outfile << identpairs.Get(j).I2();
		outfile.width (8);
		outfile << i << "\n";
	      }
	  }

	outfile << "identificationtypes\n";
	outfile << ident -> GetMaxNr() << "\n";
	for (i = 1; i <= ident -> GetMaxNr(); i++)
	  {
	    int type = ident -> GetType(i);
	    outfile << " " << type;
	  }
	outfile << "\n";
      }

    int cntmat = 0;
    for (i = 1; i <= materials.Size(); i++)
      if (materials.Get(i) && strlen (materials.Get(i)))
	cntmat++;

    if (cntmat)
      {
	outfile << "materials" << endl;
	outfile << cntmat << endl;
	for (i = 1; i <= materials.Size(); i++)
	  if (materials.Get(i) && strlen (materials.Get(i)))
	    outfile << i << " " << materials.Get(i) << endl;
      }

    
    int cntbcnames = 0;
    for ( int ii = 0; ii < bcnames.Size(); ii++ )
      if ( bcnames[ii] ) cntbcnames++;

    if ( cntbcnames )
      {
	outfile << "\n\nbcnames" << endl << bcnames.Size() << endl;
	for ( i = 0; i < bcnames.Size(); i++ )
	  outfile << i+1 << "\t" << GetBCName(i) << endl;
	outfile << endl << endl;
      }

    /*
    if ( GetDimension() == 2 )
      {
	for (i = 1; i <= GetNSeg(); i++)
	  {
	    const Segment & seg = LineSegment (i);
	    if ( ! bcprops.Contains(seg.si) && seg.GetBCName() != "" )
	      {
		bcprops.Append(seg.si);
		cntbcnames++;
	      }
	  }
      }
    else
      {
	for (sei = 0; sei < GetNSE(); sei++)
	  {
	    if ((*this)[sei].GetIndex())
	      {
		int bcp = GetFaceDescriptor((*this)[sei].GetIndex ()).BCProperty();
		string name = GetFaceDescriptor((*this)[sei].GetIndex ()).BCName();
		if ( !bcprops.Contains(bcp) &&
		     name != "" )
		  {
		    bcprops.Append(bcp);
		    cntbcnames++;
		  }
	      }
	  }
      }

    bcprops.SetSize(0);
    if ( cntbcnames )
      {
	outfile << "\nbcnames" << endl << cntbcnames << endl;
	if ( GetDimension() == 2 )
	  {
	    for (i = 1; i <= GetNSeg(); i++)
	      {
		const Segment & seg = LineSegment (i);
		if ( ! bcprops.Contains(seg.si) && seg.GetBCName() != "" )
		  {
		    bcprops.Append(seg.si);
		    outfile << seg.si << "\t" << seg.GetBCName() << endl;
		  }
	      }
	  }
	else
	  {
	    for (sei = 0; sei < GetNSE(); sei++)
	      {
		if ((*this)[sei].GetIndex())
		  {
		    int bcp = GetFaceDescriptor((*this)[sei].GetIndex ()).BCProperty();
		    string name = GetFaceDescriptor((*this)[sei].GetIndex ()).BCName();
		    if ( !bcprops.Contains(bcp) &&
			 name != "" )
		      {
			bcprops.Append(bcp);
			outfile << bcp << "\t" << name << endl;
		      }
		  }
	      }
	  }
	outfile << endl << endl;
      }
    */

    int cnt_sing = 0;
    for (PointIndex pi = PointIndex::BASE; pi < GetNP()+PointIndex::BASE; pi++)
      if ((*this)[pi].Singularity()>=1.) cnt_sing++;
    
    if (cnt_sing)
      {
	outfile << "singular_points" << endl << cnt_sing << endl;
	for (PointIndex pi = PointIndex::BASE; pi < GetNP()+PointIndex::BASE; pi++)
	  if ((*this)[pi].Singularity()>=1.) 
	    outfile << int(pi) << "\t" << (*this)[pi].Singularity() << endl;
      }

    cnt_sing = 0;
    for (SegmentIndex si = 0; si < GetNSeg(); si++)
      if ( segments[si].singedge_left ) cnt_sing++;
    if (cnt_sing)
      {
	outfile << "singular_edge_left" << endl << cnt_sing << endl;
	for (SegmentIndex si = 0; si < GetNSeg(); si++)
	  if ( segments[si].singedge_left )
	    outfile << int(si) << "\t" << segments[si].singedge_left << endl;
      }

    cnt_sing = 0;
    for (SegmentIndex si = 0; si < GetNSeg(); si++)
      if ( segments[si].singedge_right ) cnt_sing++;
    if (cnt_sing)
      {
	outfile << "singular_edge_right" << endl << cnt_sing << endl;
	for (SegmentIndex si = 0; si < GetNSeg(); si++)
	  if ( segments[si].singedge_right  )
	    outfile << int(si) << "\t" << segments[si].singedge_right << endl;
      }


    cnt_sing = 0;
    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      if ( GetFaceDescriptor ((*this)[sei].GetIndex()).domin_singular) 
	cnt_sing++;

    if (cnt_sing)
      {
	outfile << "singular_face_inside" << endl << cnt_sing << endl;
	for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
	  if ( GetFaceDescriptor ((*this)[sei].GetIndex()).domin_singular) 
	    outfile << int(sei)  << "\t" << 
	      GetFaceDescriptor ((*this)[sei].GetIndex()).domin_singular  << endl;
      }

    cnt_sing = 0;
    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      if ( GetFaceDescriptor ((*this)[sei].GetIndex()).domout_singular) cnt_sing++;
    if (cnt_sing)
      {
	outfile << "singular_face_outside" << endl << cnt_sing << endl;
	for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
	  if ( GetFaceDescriptor ((*this)[sei].GetIndex()).domout_singular) 
	    outfile << int(sei) << "\t" 
		    << GetFaceDescriptor ((*this)[sei].GetIndex()).domout_singular << endl;
      }


  }


  
  void Mesh :: Load (const string & filename)
  {
    
    ifstream infile(filename.c_str());
    if (!infile.good())
      throw NgException ("mesh file not found");

    Load(infile);
  }



 
  void Mesh :: Load (istream & infile)
  {

    char str[100];
    int i, n;

    double scale = 1;  // globflags.GetNumFlag ("scale", 1);
    int inverttets = 0;  // globflags.GetDefineFlag ("inverttets");
    int invertsurf = 0;  // globflags.GetDefineFlag ("invertsurfacemesh");


    facedecoding.SetSize(0);

    bool endmesh = false;

    while (infile.good() && !endmesh)
      {
	infile >> str;

	if (strcmp (str, "dimension") == 0)
	  {
	    infile >> dimension;
	  }

	if (strcmp (str, "geomtype") == 0)
	  {
            int hi;
	    infile >> hi;
            geomtype = GEOM_TYPE(hi);
	  }


	if (strcmp (str, "surfaceelements") == 0 || strcmp (str, "surfaceelementsgi")==0 || strcmp (str, "surfaceelementsuv") == 0)
	  {
	    infile >> n;
	    PrintMessage (3, n, " surface elements");
	    for (i = 1; i <= n; i++)
	      {
		int j;
		int surfnr, bcp, domin, domout, nep, faceind = 0;

		infile >> surfnr >> bcp >> domin >> domout;
		surfnr--;

		for (j = 1; j <= facedecoding.Size(); j++)
		  if (GetFaceDescriptor(j).SurfNr() == surfnr &&
		      GetFaceDescriptor(j).BCProperty() == bcp &&
		      GetFaceDescriptor(j).DomainIn() == domin &&
		      GetFaceDescriptor(j).DomainOut() == domout)
		    faceind = j;

		if (!faceind)
		  {
		    faceind = AddFaceDescriptor (FaceDescriptor(surfnr, domin, domout, 0));
		    GetFaceDescriptor(faceind).SetBCProperty (bcp);
		  }

		infile >> nep;
		if (!nep) nep = 3;

		Element2d tri(nep);
		tri.SetIndex(faceind);

		for (j = 1; j <= nep; j++)
		  infile >> tri.PNum(j);

		if (strcmp (str, "surfaceelementsgi") == 0)
		  for (j = 1; j <= nep; j++)
		    infile >> tri.GeomInfoPi(j).trignum;

		if (strcmp (str, "surfaceelementsuv") == 0)
		  for (j = 1; j <= nep; j++)
		    infile >> tri.GeomInfoPi(j).u >> tri.GeomInfoPi(j).v;

		if (invertsurf)
		  tri.Invert();

		AddSurfaceElement (tri);
	      }
	  }


	if (strcmp (str, "volumeelements") == 0)
	  {
	    infile >> n;
	    PrintMessage (3, n, " volume elements");
	    for (i = 1; i <= n; i++)
	      {
		Element el;
		int hi, nep;
		infile >> hi;
		if (hi == 0) hi = 1;
		el.SetIndex(hi);
		infile >> nep;
		el.SetNP(nep);
	      
		for (int j = 0; j < nep; j++)
		  infile >> (int&)(el[j]);
	      
		if (inverttets)
		  el.Invert();

		AddVolumeElement (el);
	      }
	  }
    

	if (strcmp (str, "edgesegments") == 0)
	  {
	    infile >> n;
	    for (i = 1; i <= n; i++)
	      {
		Segment seg;
		int hi;
		infile >> seg.si >> hi >> seg.p1 >> seg.p2;
		AddSegment (seg);
	      }
	  }
      


	if (strcmp (str, "edgesegmentsgi") == 0)
	  {
	    infile >> n;
	    for (i = 1; i <= n; i++)
	      {
		Segment seg;
		int hi;
		infile >> seg.si >> hi >> seg.p1 >> seg.p2
		       >> seg.geominfo[0].trignum
		       >> seg.geominfo[1].trignum;
		AddSegment (seg);
	      }
	  }
	
	if (strcmp (str, "edgesegmentsgi2") == 0)
	  {
	    int a; 
	    infile >> a;
	    n=a; 

	    PrintMessage (3, n, " curve elements");

	    for (i = 1; i <= n; i++)
	      {
		Segment seg;
		int hi;
		infile >> seg.si >> hi >> seg.p1 >> seg.p2
		       >> seg.geominfo[0].trignum
		       >> seg.geominfo[1].trignum
		       >> seg.surfnr1 >> seg.surfnr2
		       >> seg.edgenr
		       >> seg.epgeominfo[0].dist
		       >> seg.epgeominfo[1].edgenr
		       >> seg.epgeominfo[1].dist;

		seg.epgeominfo[0].edgenr = seg.epgeominfo[1].edgenr;

		seg.domin = seg.surfnr1;
		seg.domout = seg.surfnr2;

		seg.surfnr1--;
		seg.surfnr2--;
	      
		AddSegment (seg);
	      }
	  }
      
	if (strcmp (str, "points") == 0)
	  {
	    infile >> n;
	    PrintMessage (3, n, " points");
	    for (i = 1; i <= n; i++)
	      {
		Point3d p;
		infile >> p.X() >> p.Y() >> p.Z();
		p.X() *= scale;
		p.Y() *= scale;
		p.Z() *= scale;
		AddPoint (p);
	      }
	  }

	if (strcmp (str, "identifications") == 0)
	  {
	    infile >> n;
	    for (i = 1; i <= n; i++)
	      {
		PointIndex pi1, pi2;
		int ind;
		infile >> pi1 >> pi2 >> ind;
		ident -> Add (pi1, pi2, ind);
	      }
	  }
       
	if (strcmp (str, "identificationtypes") == 0)
	  {
	    infile >> n;
	    for (i = 1; i <= n; i++)
	      {
		int type;
		infile >> type;
		ident -> SetType(i,Identifications::ID_TYPE(type));
	      }
	  }

	if (strcmp (str, "materials") == 0)
	  {
	    infile >> n;
	    for (i = 1; i <= n; i++)
	      {
		int nr;
		string mat;
		infile >> nr >> mat;
		SetMaterial (nr, mat.c_str());
	      }
	  }

	if ( strcmp (str, "bcnames" ) == 0 )
	  {
	    infile >> n;
	    ARRAY<int,0> bcnrs(n);

	    SetNBCNames(n);
	    for ( i = 1; i <= n; i++ )
	      {
		string nextbcname;
		infile >> bcnrs[i-1] >> nextbcname;
		bcnames[bcnrs[i-1]-1] = new string(nextbcname);
	      }
	    
	    if ( GetDimension() == 2 )
	      {
		for (i = 1; i <= GetNSeg(); i++)
		  {
		    Segment & seg = LineSegment (i);
		    if ( seg.si <= n )
		      seg.SetBCName (bcnames[seg.si-1]);
		    else
		      seg.SetBCName(0);
		  }
	      }
	    else
	      {
		for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
		  {
		    if ((*this)[sei].GetIndex())
		      {
			int bcp = GetFaceDescriptor((*this)[sei].GetIndex ()).BCProperty();
			if ( bcp <= n )
			  GetFaceDescriptor((*this)[sei].GetIndex ()).SetBCName(bcnames[bcp-1]);
			else
			  GetFaceDescriptor((*this)[sei].GetIndex ()).SetBCName(0);

		      }
		  }

	      }
	    

	  }
	
	if (strcmp (str, "singular_points") == 0)
	  {
	    infile >> n;
	    for (i = 1; i <= n; i++)
	      {
		PointIndex pi;
		double s; 
		infile >> pi;
		infile >> s; 
		(*this)[pi].Singularity (s);
	      }
	  }

	if (strcmp (str, "singular_edge_left") == 0)
	  {
	    infile >> n;
	    for (i = 1; i <= n; i++)
	      {
		SegmentIndex si;
		double s; 
		infile >> si;
		infile >> s; 
		(*this)[si].singedge_left = s;
	      }
	  }
	if (strcmp (str, "singular_edge_right") == 0)
	  {
	    infile >> n;
	    for (i = 1; i <= n; i++)
	      {
		SegmentIndex si;
		double s; 
		infile >> si;
		infile >> s; 
		(*this)[si].singedge_right = s;
	      }
	  }

	if (strcmp (str, "singular_face_inside") == 0)
	  {
	    infile >> n;
	    for (i = 1; i <= n; i++)
	      {
		SurfaceElementIndex sei;
		double s; 
		infile >> sei;
		infile >> s; 
		GetFaceDescriptor((*this)[sei].GetIndex()).domin_singular = s;
	      }
	  }

	if (strcmp (str, "singular_face_outside") == 0)
	  {
	    infile >> n;
	    for (i = 1; i <= n; i++)
	      {
		SurfaceElementIndex sei;
		double s; 
		infile >> sei;
		infile >> s; 
		GetFaceDescriptor((*this)[sei].GetIndex()).domout_singular = s;
	      }
	  }

	if (strcmp (str, "endmesh") == 0)
	  endmesh = true;
	  


	strcpy (str, "");
      }
  
    CalcSurfacesOfNode ();
    //  BuildConnectedNodes ();
    topology -> Update();
    clusters -> Update();
  
    SetNextMajorTimeStamp();
    //  PrintMemInfo (cout);


#ifdef PARALLEL
    if ( ntasks > 1 )
      {
	// for parallel processing
	Distribute ();
	return;
      }
#endif

  }
  

  


  void Mesh :: Merge (const string & filename, const int surfindex_offset)
  {
    ifstream infile(filename.c_str());
    if (!infile.good())
      throw NgException ("mesh file not found");

    Merge(infile,surfindex_offset);

  }



  void Mesh :: Merge (istream & infile, const int surfindex_offset)
  {
    char str[100];
    int i, n;

    
    int inverttets = 0;  // globflags.GetDefineFlag ("inverttets");

    int oldnp = GetNP();
    int oldne = GetNSeg();
    int oldnd = GetNDomains();
    
    for(SurfaceElementIndex si = 0; si < GetNSE(); si++)
      for(int j=1; j<=(*this)[si].GetNP(); j++) (*this)[si].GeomInfoPi(j).trignum = -1;
    
    int max_surfnr = 0;
    for (i = 1; i <= GetNFD(); i++)
      max_surfnr = max2 (max_surfnr, GetFaceDescriptor(i).SurfNr());
    max_surfnr++;

    if(max_surfnr < surfindex_offset) max_surfnr = surfindex_offset;


    bool endmesh = false;

    while (infile.good() && !endmesh)
      {
	infile >> str;

	if (strcmp (str, "surfaceelementsgi") == 0 || strcmp (str, "surfaceelements") == 0)
	  {
	    infile >> n;
	    PrintMessage (3, n, " surface elements");
	    for (i = 1; i <= n; i++)
	      {
		int j;
		int surfnr, bcp, domin, domout, nep, faceind = 0;
		infile >> surfnr >> bcp >> domin >> domout;

		surfnr--;

		if(domin > 0) domin += oldnd;
		if(domout > 0) domout += oldnd;
		surfnr += max_surfnr;
		
		
		for (j = 1; j <= facedecoding.Size(); j++)
		  if (GetFaceDescriptor(j).SurfNr() == surfnr &&
		      GetFaceDescriptor(j).BCProperty() == bcp &&
		      GetFaceDescriptor(j).DomainIn() == domin &&
		      GetFaceDescriptor(j).DomainOut() == domout)
		    faceind = j;

		if (!faceind)
		  {
		    faceind = AddFaceDescriptor (FaceDescriptor(surfnr, domin, domout, 0));
		    if(GetDimension() == 2) bcp++;
		    GetFaceDescriptor(faceind).SetBCProperty (bcp);
		  }

		infile >> nep;
		if (!nep) nep = 3;

		Element2d tri(nep);
		tri.SetIndex(faceind);

		for (j = 1; j <= nep; j++)
		  {
		    infile >> tri.PNum(j);
		    tri.PNum(j) = tri.PNum(j) + oldnp;
		  }

		
		if (strcmp (str, "surfaceelementsgi") == 0)
		  for (j = 1; j <= nep; j++)
		  {
		    infile >> tri.GeomInfoPi(j).trignum;
		    tri.GeomInfoPi(j).trignum = -1;
		  }

		AddSurfaceElement (tri);
	      }
	  }

	
	if (strcmp (str, "edgesegments") == 0)
	  {
	    infile >> n;
	    for (i = 1; i <= n; i++)
	      {
		Segment seg;
		int hi;
		infile >> seg.si >> hi >> seg.p1 >> seg.p2;
		seg.p1 = seg.p1 + oldnp;
		seg.p2 = seg.p2 + oldnp;
		AddSegment (seg);
	      }
	  }
      


	if (strcmp (str, "edgesegmentsgi") == 0)
	  {
	    infile >> n;
	    for (i = 1; i <= n; i++)
	      {
		Segment seg;
		int hi;
		infile >> seg.si >> hi >> seg.p1 >> seg.p2
		       >> seg.geominfo[0].trignum
		       >> seg.geominfo[1].trignum;
		seg.p1 = seg.p1 + oldnp;
		seg.p2 = seg.p2 + oldnp;
		AddSegment (seg);
	      }
	  }
	if (strcmp (str, "edgesegmentsgi2") == 0)
	  {
	    infile >> n;
	    PrintMessage (3, n, " curve elements");

	    for (i = 1; i <= n; i++)
	      {
		Segment seg;
		int hi;
		infile >> seg.si >> hi >> seg.p1 >> seg.p2
		       >> seg.geominfo[0].trignum
		       >> seg.geominfo[1].trignum
		       >> seg.surfnr1 >> seg.surfnr2
		       >> seg.edgenr
		       >> seg.epgeominfo[0].dist
		       >> seg.epgeominfo[1].edgenr
		       >> seg.epgeominfo[1].dist;
		seg.epgeominfo[0].edgenr = seg.epgeominfo[1].edgenr;

		seg.surfnr1--;
		seg.surfnr2--;

		if(seg.surfnr1 >= 0)  seg.surfnr1 = seg.surfnr1 + max_surfnr;
		if(seg.surfnr2 >= 0)  seg.surfnr2 = seg.surfnr2 + max_surfnr;
		seg.p1 = seg.p1 +oldnp;
		seg.p2 = seg.p2 +oldnp;
		seg.edgenr = seg.edgenr + oldne;
		seg.epgeominfo[1].edgenr = seg.epgeominfo[1].edgenr + oldne;

		AddSegment (seg);
	      }
	  }
 
	if (strcmp (str, "volumeelements") == 0)
	  {
	    infile >> n;
	    PrintMessage (3, n, " volume elements");
	    for (i = 1; i <= n; i++)
	      {
		Element el;
		int hi, nep;
		infile >> hi;
		if (hi == 0) hi = 1;
		el.SetIndex(hi+oldnd);
		infile >> nep;
		el.SetNP(nep);
	      
		for (int j = 0; j < nep; j++)
		  {
		    infile >> (int&)(el[j]);
		    el[j] = el[j]+oldnp;
		  }
	      
		if (inverttets)
		  el.Invert();

		AddVolumeElement (el);
	      }
	  }
         

	if (strcmp (str, "points") == 0)
	  {
	    infile >> n;
	    PrintMessage (3, n, " points");
	    for (i = 1; i <= n; i++)
	      {
		Point3d p;
		infile >> p.X() >> p.Y() >> p.Z();
		AddPoint (p);
	      }
	  }


	if (strcmp (str, "endmesh") == 0)
	  {
	    endmesh = true;
	  }


	if (strcmp (str, "materials") == 0)
	  {
	    infile >> n;
	    for (i = 1; i <= n; i++)
	      {
		int nr;
		string mat;
		infile >> nr >> mat;
		SetMaterial (nr+oldnd, mat.c_str());
	      }
	  }


	strcpy (str, "");
      }
  
    CalcSurfacesOfNode ();

    topology -> Update();
    clusters -> Update();
  
    SetNextMajorTimeStamp();
  }
  







   

  bool Mesh :: TestOk () const
  {
    for (ElementIndex ei = 0; ei < volelements.Size(); ei++)
      {
	for (int j = 0; j < 4; j++)
	  if ( (*this)[ei][j] <= PointIndex::BASE-1)
	    {
	      (*testout) << "El " << ei << " has 0 nodes: ";
	      for (int k = 0; k < 4; k++)
		(*testout) << (*this)[ei][k];
	      break;
	    }
      }
    CheckMesh3D (*this);
    return 1;
  }

  void Mesh :: SetAllocSize(int nnodes, int nsegs, int nsel, int nel)
  {
    points.SetAllocSize(nnodes);
    segments.SetAllocSize(nsegs);
    surfelements.SetAllocSize(nsel);
    volelements.SetAllocSize(nel);
  }
  

  void Mesh :: BuildBoundaryEdges(void)
  {
    delete boundaryedges;

    boundaryedges = new INDEX_2_CLOSED_HASHTABLE<int>
      (3 * (GetNSE() + GetNOpenElements()) + GetNSeg() + 1);


    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      {
	const Element2d & sel = surfelements[sei];
	if (sel.IsDeleted()) continue;

	int si = sel.GetIndex();
      
	for (int j = 0; j < sel.GetNP(); j++)
	  {
	    INDEX_2 i2;
	    i2.I1() = sel.PNumMod(j+1);
	    i2.I2() = sel.PNumMod(j+2);
	    i2.Sort();
	    if (sel.GetNP() <= 4)
	      boundaryedges->Set (i2, 1);
	  }
      }


    for (int i = 0; i < openelements.Size(); i++)
      {
	const Element2d & sel = openelements[i];
	for (int j = 0; j < sel.GetNP(); j++)
	  {
	    INDEX_2 i2;
	    i2.I1() = sel.PNumMod(j+1);
	    i2.I2() = sel.PNumMod(j+2);
	    i2.Sort();
	    boundaryedges->Set (i2, 1);

	    points[sel[j]].SetType(FIXEDPOINT);
	  }
      }

    for (int i = 0; i < GetNSeg(); i++)
      {
	const Segment & seg = segments[i];
	INDEX_2 i2(seg.p1, seg.p2);
	i2.Sort();

	boundaryedges -> Set (i2, 2);
	//segmentht -> Set (i2, i);
      }
    
    
  }

  void Mesh :: CalcSurfacesOfNode ()
  {
    int i, j, k;
    SurfaceElementIndex sei;

    surfacesonnode.SetSize (GetNP());

    delete boundaryedges;
    boundaryedges = NULL;

    delete surfelementht;
    delete segmentht;

    /*
      surfelementht = new INDEX_3_HASHTABLE<int> (GetNSE()/4 + 1);
      segmentht = new INDEX_2_HASHTABLE<int> (GetNSeg() + 1);
    */

    surfelementht = new INDEX_3_CLOSED_HASHTABLE<int> (3*GetNSE() + 1);
    segmentht = new INDEX_2_CLOSED_HASHTABLE<int> (3*GetNSeg() + 1);

    for (sei = 0; sei < GetNSE(); sei++)
      {
	const Element2d & sel = surfelements[sei];
	if (sel.IsDeleted()) continue;

	int si = sel.GetIndex();
      
	for (j = 0; j < sel.GetNP(); j++)
	  {
	    PointIndex pi = sel[j];
	    bool found = 0;
	    for (k = 0; k < surfacesonnode[pi].Size(); k++)
	      if (surfacesonnode[pi][k] == si)
		{
		  found = 1;
		  break;
		}
	  
	    if (!found)
	      surfacesonnode.Add (pi, si);

	  }
      }
    /*
      for (sei = 0; sei < GetNSE(); sei++)
      {
      const Element2d & sel = surfelements[sei];
      if (sel.IsDeleted()) continue;

      INDEX_3 i3;
      i3.I1() = sel.PNum(1);
      i3.I2() = sel.PNum(2);
      i3.I3() = sel.PNum(3);
      i3.Sort();
      surfelementht -> PrepareSet (i3);
      }

      surfelementht -> AllocateElements();
    */
    for (sei = 0; sei < GetNSE(); sei++)
      {
	const Element2d & sel = surfelements[sei];
	if (sel.IsDeleted()) continue;

	INDEX_3 i3;
	i3.I1() = sel.PNum(1);
	i3.I2() = sel.PNum(2);
	i3.I3() = sel.PNum(3);
	i3.Sort();
	surfelementht -> Set (i3, sei);   // war das wichtig ???    sel.GetIndex());
      }

    int np = GetNP();

    if (dimension == 3)
      {
        for (PointIndex pi = PointIndex::BASE; 
             pi < np+PointIndex::BASE; pi++)
          points[pi].SetType (INNERPOINT);
        
        if (GetNFD() == 0) 
          {
            for (sei = 0; sei < GetNSE(); sei++)
              {
                const Element2d & sel = surfelements[sei];
                if (sel.IsDeleted()) continue;
                for (j = 0;  j < sel.GetNP(); j++)
                  {
                    PointIndex pi = SurfaceElement(sei)[j];
                    points[pi].SetType(FIXEDPOINT);
                  }
              }
          }
        else
          {
            for (sei = 0; sei < GetNSE(); sei++)
              {
                const Element2d & sel = surfelements[sei];
                if (sel.IsDeleted()) continue;
                for (j = 0; j < sel.GetNP(); j++)
                  {
                    PointIndex pi = sel[j];
                    int ns = surfacesonnode[pi].Size();
                    if (ns == 1)
                      points[pi].SetType(SURFACEPOINT);
                    if (ns == 2)
                      points[pi].SetType(EDGEPOINT);
                    if (ns >= 3)
                      points[pi].SetType(FIXEDPOINT);
                  }      
              }
          }

    for (i = 0; i < segments.Size(); i++)
      {
	const Segment & seg = segments[i];
	for (j = 1; j <= 2; j++)
	  {
	    PointIndex hi = (j == 1) ? seg.p1 : seg.p2;
	  
	    if (points[hi].Type() == INNERPOINT ||
		points[hi].Type() == SURFACEPOINT)
	      points[hi].SetType(EDGEPOINT);
	  }
      }


    for (i = 0; i < lockedpoints.Size(); i++)
      points[lockedpoints[i]].SetType(FIXEDPOINT);
      }

  
    /*
    for (i = 0; i < openelements.Size(); i++)
      {
	const Element2d & sel = openelements[i];
	for (j = 0; j < sel.GetNP(); j++)
	  {
	    INDEX_2 i2;
	    i2.I1() = sel.PNumMod(j+1);
	    i2.I2() = sel.PNumMod(j+2);
	    i2.Sort();
	    boundaryedges->Set (i2, 1);

	    points[sel[j]].SetType(FIXEDPOINT);
	  }
      }
    */

    // eltyps.SetSize (GetNE());
    // eltyps = FREEELEMENT;

    for (i = 0; i < GetNSeg(); i++)
      {
	const Segment & seg = segments[i];
	INDEX_2 i2(seg.p1, seg.p2);
	i2.Sort();

	//boundaryedges -> Set (i2, 2);
	segmentht -> Set (i2, i);
      }
  }


  void Mesh :: FixPoints (const BitArray & fixpoints)
  {
    if (fixpoints.Size() != GetNP())
      {
	cerr << "Mesh::FixPoints: sizes don't fit" << endl;
	return;
      }
    int np = GetNP();
    for (int i = 1; i <= np; i++)
      if (fixpoints.Test(i))
	{
	  points.Elem(i).SetType (FIXEDPOINT);
	}
  }


  void Mesh :: FindOpenElements (int dom)
  {
    static int timer = NgProfiler::CreateTimer ("Mesh::FindOpenElements");
    static int timera = NgProfiler::CreateTimer ("Mesh::FindOpenElements A");
    static int timerb = NgProfiler::CreateTimer ("Mesh::FindOpenElements B");
    static int timerc = NgProfiler::CreateTimer ("Mesh::FindOpenElements C");
    static int timerd = NgProfiler::CreateTimer ("Mesh::FindOpenElements D");
    static int timere = NgProfiler::CreateTimer ("Mesh::FindOpenElements E");

    NgProfiler::RegionTimer reg (timer);

    int np = GetNP();
    int ne = GetNE();
    int nse = GetNSE();
    
    ARRAY<int,PointIndex::BASE> numonpoint(np);

    numonpoint = 0;

    NgProfiler::StartTimer (timera);
    for (ElementIndex ei = 0; ei < ne; ei++)
      {
	const Element & el = (*this)[ei];
	if (dom == 0 || dom == el.GetIndex())
	  {
	    if (el.GetNP() == 4)
	      {
		INDEX_4 i4(el[0], el[1], el[2], el[3]);
		i4.Sort();
		numonpoint[i4.I1()]++;
		numonpoint[i4.I2()]++;
	      }
	    else
	      for (int j = 0; j < el.GetNP(); j++)
		numonpoint[el[j]]++;
	  }
      }

    TABLE<ElementIndex,PointIndex::BASE> elsonpoint(numonpoint);
    for (ElementIndex ei = 0; ei < ne; ei++)
      {
	const Element & el = (*this)[ei];
	if (dom == 0 || dom == el.GetIndex())
	  {
	    if (el.GetNP() == 4)
	      {
		INDEX_4 i4(el[0], el[1], el[2], el[3]);
		i4.Sort();
		elsonpoint.Add (i4.I1(), ei);
		elsonpoint.Add (i4.I2(), ei);
	      }
	    else
	      for (int j = 0; j < el.GetNP(); j++)
		elsonpoint.Add (el[j], ei);
	  }
      }
    NgProfiler::StopTimer (timera);




    NgProfiler::StartTimer (timerb);	



    ARRAY<char, 1> hasface(GetNFD());

    int i;
    for (i = 1; i <= GetNFD(); i++)
      {
	int domin = GetFaceDescriptor(i).DomainIn();
	int domout = GetFaceDescriptor(i).DomainOut();
	hasface[i] = 
	  dom == 0 && (domin != 0 || domout != 0) ||
	  dom != 0 && (domin == dom || domout == dom);
      }

    numonpoint = 0;
    for (SurfaceElementIndex sii = 0; sii < nse; sii++)
      {
	int ind = surfelements[sii].GetIndex();
	/*
	  if (
	  GetFaceDescriptor(ind).DomainIn() && 
	  (dom == 0 || dom == GetFaceDescriptor(ind).DomainIn())
	  ||
	  GetFaceDescriptor(ind).DomainOut() && 
	  (dom == 0 || dom == GetFaceDescriptor(ind).DomainOut())
	  )
	*/
	if (hasface[ind])
	  {
	    /*
	      Element2d hel = surfelements[i];
	      hel.NormalizeNumbering();	  
	      numonpoint[hel[0]]++;
	    */
	    const Element2d & hel = surfelements[sii];
	    int mini = 0;
	    for (int j = 1; j < hel.GetNP(); j++)
	      if (hel[j] < hel[mini])
		mini = j;
	    numonpoint[hel[mini]]++;
	  }
      }

    TABLE<SurfaceElementIndex,PointIndex::BASE> selsonpoint(numonpoint);
    for (SurfaceElementIndex sii = 0; sii < nse; sii++)
      {
	int ind = surfelements[sii].GetIndex();

	/*
	  if (
	  GetFaceDescriptor(ind).DomainIn() && 
	  (dom == 0 || dom == GetFaceDescriptor(ind).DomainIn())
	  ||
	  GetFaceDescriptor(ind).DomainOut() && 
	  (dom == 0 || dom == GetFaceDescriptor(ind).DomainOut())
	  )
	*/
	if (hasface[ind])
	  {
	    /*
	      Element2d hel = surfelements[i];
	      hel.NormalizeNumbering();	  
	      selsonpoint.Add (hel[0], i);
	    */
	    const Element2d & hel = surfelements[sii];
	    int mini = 0;
	    for (int j = 1; j < hel.GetNP(); j++)
	      if (hel[j] < hel[mini])
		mini = j;
	    selsonpoint.Add (hel[mini], sii);
	  }
      }


    NgProfiler::StopTimer (timerb);

    int ii, j, k, l;
    PointIndex pi;
    SurfaceElementIndex sei;
    Element2d hel;

    NgProfiler::RegionTimer regc (timerc);


    INDEX_3_CLOSED_HASHTABLE<INDEX_2> faceht(100);   
    openelements.SetSize(0);
      
    for (pi = PointIndex::BASE; pi < np+PointIndex::BASE; pi++)
      if (selsonpoint[pi].Size()+elsonpoint[pi].Size())
	{
	  faceht.SetSize (2 * selsonpoint[pi].Size() + 4 * elsonpoint[pi].Size());

	  FlatArray<SurfaceElementIndex> row = selsonpoint[pi];
	  for (ii = 0; ii < row.Size(); ii++)
	    {
	      hel = SurfaceElement(row[ii]);
	      int ind = hel.GetIndex();	  
  
	      if (GetFaceDescriptor(ind).DomainIn() && 
		  (dom == 0 || dom == GetFaceDescriptor(ind).DomainIn()) )
		{
		  hel.NormalizeNumbering();
		  if (hel.PNum(1) == pi)
		    {
		      INDEX_3 i3(hel[0], hel[1], hel[2]);
		      INDEX_2 i2 (GetFaceDescriptor(ind).DomainIn(), 
				  (hel.GetNP() == 3) 
				  ? PointIndex (PointIndex::BASE-1)
				  : hel.PNum(4));
		      faceht.Set (i3, i2);
		    }
		}
	      if (GetFaceDescriptor(ind).DomainOut() &&
		  (dom == 0 || dom == GetFaceDescriptor(ind).DomainOut()) )
		{
		  hel.Invert();
		  hel.NormalizeNumbering();
		  if (hel.PNum(1) == pi)
		    {
		      INDEX_3 i3(hel[0], hel[1], hel[2]);
		      INDEX_2 i2 (GetFaceDescriptor(ind).DomainOut(), 
				  (hel.GetNP() == 3) 
				  ? PointIndex (PointIndex::BASE-1)
				  : hel.PNum(4));
		      faceht.Set (i3, i2);
		    }
		}
	    }

	  
	  FlatArray<ElementIndex> rowel = elsonpoint[pi];
	  for (ii = 0; ii < rowel.Size(); ii++)
	    {
	      const Element & el = VolumeElement(rowel[ii]);

	      if (dom == 0 || el.GetIndex() == dom)
		{
		  for (j = 1; j <= el.GetNFaces(); j++)
		    {
		      el.GetFace (j, hel);
		      hel.Invert();
		      hel.NormalizeNumbering();

		      if (hel[0] == pi)
			{
			  INDEX_3 i3(hel[0], hel[1], hel[2]);
			  
			  if (faceht.Used (i3))
			    {
			      INDEX_2 i2 = faceht.Get(i3);
			      if (i2.I1() == el.GetIndex())
				{
				  i2.I1() = PointIndex::BASE-1;
				  faceht.Set (i3, i2);
				}
			      else
				{
				  if (i2.I1() == 0)
				    {
				      PrintSysError ("more elements on face");
				      (*testout)  << "more elements on face!!!" << endl;
				      (*testout) << "el = " << el << endl;
				      (*testout) << "hel = " << hel << endl;
				      (*testout) << "face = " << i3 << endl;
				      (*testout) << "points = " << endl;
				      for (int jj = 1; jj <= 3; jj++)
					(*testout) << "p = " << Point(i3.I(jj)) << endl;
				    }
				}
			    }
			  else
			    {
			      hel.Invert();
			      hel.NormalizeNumbering();
			      INDEX_3 i3(hel[0], hel[1], hel[2]);
			      INDEX_2 i2(el.GetIndex(), 
					 (hel.GetNP() == 3) 
					 ? PointIndex (PointIndex::BASE-1)
					 : hel[3]);
			      faceht.Set (i3, i2);
			    }
			}
		    }
		}
	    }
	  for (i = 0; i < faceht.Size(); i++)
	    if (faceht.UsedPos (i))
	      {
		INDEX_3 i3;
		INDEX_2 i2;
		faceht.GetData (i, i3, i2);
		if (i2.I1() != PointIndex::BASE-1)
		  {
		    Element2d tri;
		    tri.SetType ( (i2.I2() == PointIndex::BASE-1) ? TRIG : QUAD);
		    for (l = 0; l < 3; l++)
		      tri[l] = i3.I(l+1);
		    tri.PNum(4) = i2.I2();
		    tri.SetIndex (i2.I1());

		    //	tri.Invert();
		    
		    openelements.Append (tri);
		  }
	      }
	}

    int cnt3 = 0;
    for (i = 0; i < openelements.Size(); i++)
      if (openelements[i].GetNP() == 3)
	cnt3++;

    int cnt4 = openelements.Size() - cnt3;


    MyStr treequad;
    if (cnt4)
      treequad = MyStr(" (") + MyStr(cnt3) + MyStr (" + ") + 
	MyStr(cnt4) + MyStr(")");

    PrintMessage (5, openelements.Size(), treequad, " open elements");

    BuildBoundaryEdges();


    NgProfiler::RegionTimer regd (timerd);

    for (i = 1; i <= openelements.Size(); i++)
      {
	const Element2d & sel = openelements.Get(i);

	if (boundaryedges)
	  for (j = 1; j <= sel.GetNP(); j++)
	    {
	      INDEX_2 i2;
	      i2.I1() = sel.PNumMod(j);
	      i2.I2() = sel.PNumMod(j+1);
	      i2.Sort();
	      boundaryedges->Set (i2, 1);
	    }
      
	for (j = 1; j <= 3; j++)
	  {
	    int pi = sel.PNum(j);
	    if (pi < points.Size()+PointIndex::BASE)
	      points[pi].SetType (FIXEDPOINT);
	  }
      }


    NgProfiler::RegionTimer rege (timere);

    /*
      for (i = 1; i <= GetNSeg(); i++)
      {
      const Segment & seg = LineSegment(i);
      INDEX_2 i2(seg.p1, seg.p2);
      i2.Sort();

      if (!boundaryedges->Used (i2))
      cerr << "WARNING: no boundedge, but seg edge: " << i2 << endl;

      boundaryedges -> Set (i2, 2);
      segmentht -> Set (i2, i-1);
      }
    */
  }

  bool Mesh :: HasOpenQuads () const
  {
    int no = GetNOpenElements();
    for (int i = 0; i < no; i++)
      if (openelements[i].GetNP() == 4)
	return true;
    return false;
  }





  void Mesh :: FindOpenSegments (int surfnr)
  {
    int i, j, k;

    // new version, general elemetns
    // hash index: pnum1-2
    // hash data : surfnr,  surfel-nr (pos) or segment nr(neg)
    INDEX_2_HASHTABLE<INDEX_2> faceht(4 * GetNSE()+GetNSeg()+1);   
  
    PrintMessage (5, "Test Opensegments");
    for (i = 1; i <= GetNSeg(); i++)
      {
	const Segment & seg = LineSegment (i);

	if (surfnr == 0 || seg.si == surfnr)
	  {
	    INDEX_2 key(seg.p1, seg.p2);
	    INDEX_2 data(seg.si, -i);

	    if (faceht.Used (key))
	      {
		cerr << "ERROR: Segment " << seg << " already used" << endl;
		(*testout) << "ERROR: Segment " << seg << " already used" << endl;
	      }

	    faceht.Set (key, data);
	  }
      }


    for (i = 1; i <= GetNSeg(); i++)
      {
	const Segment & seg = LineSegment (i);

	if (surfnr == 0 || seg.si == surfnr)
	  {
	    INDEX_2 key(seg.p2, seg.p1);
	    if (!faceht.Used(key))
	      {
		cerr << "ERROR: Segment " << seg << " brother not used" << endl;
		(*testout) << "ERROR: Segment " << seg << " brother not used" << endl;
	      }
	  }
      }


    for (i = 1; i <= GetNSE(); i++)
      {
	const Element2d & el = SurfaceElement(i);
	if (el.IsDeleted()) continue;

	if (surfnr == 0 || el.GetIndex() == surfnr)
	  {
	    for (j = 1; j <= el.GetNP(); j++)
	      {
		INDEX_2 seg (el.PNumMod(j), el.PNumMod(j+1));
		INDEX_2 data;

		if (seg.I1() <= 0 || seg.I2() <= 0)
		  cerr << "seg = " << seg << endl;

		if (faceht.Used(seg))
		  {
		    data = faceht.Get(seg);
		    if (data.I1() == el.GetIndex())
		      {
			data.I1() = 0;
			faceht.Set (seg, data);
		      }
		    else
		      {
			PrintSysError ("hash table si not fitting for segment: ",
				       seg.I1(), "-", seg.I2(), " other = ",
				       data.I2());
		      }
		  }
		else
		  {
		    Swap (seg.I1(), seg.I2());
		    data.I1() = el.GetIndex();
		    data.I2() = i;

		    faceht.Set (seg, data);
		  }
	      }
	  }
      }  

    (*testout) << "open segments: " << endl;
    opensegments.SetSize(0);
    for (i = 1; i <= faceht.GetNBags(); i++)
      for (j = 1; j <= faceht.GetBagSize(i); j++)
	{
	  INDEX_2 i2;
	  INDEX_2 data;
	  faceht.GetData (i, j, i2, data);
	  if (data.I1())  // surfnr
	    {
	      Segment seg;
	      seg.p1 = i2.I1();
	      seg.p2 = i2.I2();
	      seg.si = data.I1();

	      // find geomdata:
	      if (data.I2() > 0)
		{
		  // segment due to triangle
		  const Element2d & el = SurfaceElement (data.I2());
		  for (k = 1; k <= el.GetNP(); k++)
		    {
		      if (seg.p1 == el.PNum(k))
			seg.geominfo[0] = el.GeomInfoPi(k);
		      if (seg.p2 == el.PNum(k))
			seg.geominfo[1] = el.GeomInfoPi(k);
		    }

		  (*testout) << "trig seg: ";
		}
	      else
		{
		  // segment due to line
		  const Segment & lseg = LineSegment (-data.I2());
		  seg.geominfo[0] = lseg.geominfo[0];
		  seg.geominfo[1] = lseg.geominfo[1];

		  (*testout) << "line seg: ";
		}
	    
	      (*testout) << seg.p1 << " - " << seg.p2 
			 << " len = " << Dist (Point(seg.p1), Point(seg.p2))
			 << endl;

	      opensegments.Append (seg);
	      if (seg.geominfo[0].trignum <= 0 || seg.geominfo[1].trignum <= 0)
		{
		  (*testout) << "Problem with open segment: " << seg << endl;
		}

	    }
	}

    PrintMessage (3, opensegments.Size(), " open segments found");
    (*testout) << opensegments.Size() << " open segments found" << endl;
  
    /*
      ptyps.SetSize (GetNP());
      for (i = 1; i <= ptyps.Size(); i++)
      ptyps.Elem(i) = SURFACEPOINT;

      for (i = 1; i <= GetNSeg(); i++)
      {
      const Segment & seg = LineSegment (i);
      ptyps.Elem(seg.p1) = EDGEPOINT;
      ptyps.Elem(seg.p2) = EDGEPOINT;
      }
      for (i = 1; i <= GetNOpenSegments(); i++)
      {
      const Segment & seg = GetOpenSegment (i);
      ptyps.Elem(seg.p1) = EDGEPOINT;
      ptyps.Elem(seg.p2) = EDGEPOINT;
      }
    */
    for (i = 1; i <= points.Size(); i++)
      points.Elem(i).SetType(SURFACEPOINT);

    for (i = 1; i <= GetNSeg(); i++)
      {
	const Segment & seg = LineSegment (i);
	points[seg.p1].SetType(EDGEPOINT);
	points[seg.p2].SetType(EDGEPOINT);
      }
    for (i = 1; i <= GetNOpenSegments(); i++)
      {
	const Segment & seg = GetOpenSegment (i);
	points[seg.p1].SetType (EDGEPOINT);
	points[seg.p2].SetType (EDGEPOINT);
      }
  
  
  
    /*

    for (i = 1; i <= openelements.Size(); i++)
    {
    const Element2d & sel = openelements.Get(i);

    if (boundaryedges)
    for (j = 1; j <= sel.GetNP(); j++)
    {
    INDEX_2 i2;
    i2.I1() = sel.PNumMod(j);
    i2.I2() = sel.PNumMod(j+1);
    i2.Sort();
    boundaryedges->Set (i2, 1);
    }
      
    for (j = 1; j <= 3; j++)
    {
    int pi = sel.PNum(j);
    if (pi <= ptyps.Size())
    ptyps.Elem(pi) = FIXEDPOINT;
    }
    }
    */
  }


  void Mesh :: RemoveOneLayerSurfaceElements ()
  {
    int i, j;
    int np = GetNP();

    FindOpenSegments();
    BitArray frontpoints(np);

    frontpoints.Clear();
    for (i = 1; i <= GetNOpenSegments(); i++)
      {
	const Segment & seg = GetOpenSegment(i);
	frontpoints.Set (seg.p1);
	frontpoints.Set (seg.p2);
      }

    for (i = 1; i <= GetNSE(); i++)
      {
	Element2d & sel = surfelements.Elem(i);
	int remove = 0;
	for (j = 1; j <= sel.GetNP(); j++)
	  if (frontpoints.Test(sel.PNum(j)))
	    remove = 1;
	if (remove)
	  sel.PNum(1) = 0;
      }

    for (i = surfelements.Size(); i >= 1; i--)
      {
	if (surfelements.Elem(i).PNum(1) == 0)
	  {
	    surfelements.Elem(i) = surfelements.Last();
	    surfelements.DeleteLast();
	  }
      }

    for (int i = 0; i < facedecoding.Size(); i++)
      facedecoding[i].firstelement = -1;
    for (int i = surfelements.Size()-1; i >= 0; i--)
      {
        int ind = surfelements[i].GetIndex();
        surfelements[i].next = facedecoding[ind-1].firstelement;
        facedecoding[ind-1].firstelement = i;
      }


    timestamp = NextTimeStamp();
    //  Compress();
  }





  void Mesh :: FreeOpenElementsEnvironment (int layers)
  {
    int i, j, k;
    PointIndex pi;
    const int large = 9999;
    ARRAY<int,PointIndex::BASE> dist(GetNP());

    dist = large;

    for (int i = 1; i <= GetNOpenElements(); i++)
      {
	const Element2d & face = OpenElement(i);
	for (j = 0; j < face.GetNP(); j++)
	  dist[face[j]] = 1;
      }

    for (k = 1; k <= layers; k++)
      for (i = 1; i <= GetNE(); i++)
	{
	  const Element & el = VolumeElement(i);
	  if (el[0] == -1 || el.IsDeleted()) continue;

	  int elmin = large;
	  for (j = 0; j < el.GetNP(); j++)
	    if (dist[el[j]] < elmin)
	      elmin = dist[el[j]];

	  if (elmin < large)
	    {
	      for (j = 0; j < el.GetNP(); j++)
		if (dist[el[j]] > elmin+1)
		  dist[el[j]] = elmin+1;
	    }
	}

    int cntfree = 0;
    for (i = 1; i <= GetNE(); i++)
      {
	Element & el = VolumeElement(i);
	if (el[0] == -1 || el.IsDeleted()) continue;
	
	int elmin = large;
	for (j = 0; j < el.GetNP(); j++)
	  if (dist[el[j]] < elmin)
	    elmin = dist[el[j]];
      
        el.flags.fixed = elmin > layers;
	// eltyps.Elem(i) = (elmin <= layers) ? 
        // FREEELEMENT : FIXEDELEMENT;
	if (elmin <= layers)
	  cntfree++;
      }

    PrintMessage (5, "free: ", cntfree, ", fixed: ", GetNE()-cntfree);
    (*testout) << "free: " << cntfree << ", fixed: " << GetNE()-cntfree << endl;

    for (pi = PointIndex::BASE; 
	 pi < GetNP()+PointIndex::BASE; pi++)
      {
	if (dist[pi] > layers+1)
	  points[pi].SetType(FIXEDPOINT);
      }
  }



  void Mesh :: SetLocalH (const Point3d & pmin, const Point3d & pmax, double grading)
  {
    Point3d c = Center (pmin, pmax);
    double d = max3 (pmax.X()-pmin.X(),
		     pmax.Y()-pmin.Y(),
		     pmax.Z()-pmin.Z());
    d /= 2;
    Point3d pmin2 = c - Vec3d (d, d, d);
    Point3d pmax2 = c + Vec3d (d, d, d);
  

    delete lochfunc;
    lochfunc = new LocalH (pmin2, pmax2, grading);
  }

  void Mesh :: RestrictLocalH (const Point3d & p, double hloc)
  {
    if(hloc < hmin)
      hloc = hmin;

    //cout << "restrict h in " << p << " to " << hloc << endl;
    if (!lochfunc)
      {
	PrintWarning("RestrictLocalH called, creating mesh-size tree");

	Point3d boxmin, boxmax;
	GetBox (boxmin, boxmax);
	SetLocalH (boxmin, boxmax, 0.8);
      }

    lochfunc -> SetH (p, hloc);
  }

  void Mesh :: RestrictLocalHLine (const Point3d & p1, 
				   const Point3d & p2,
				   double hloc)
  {
    if(hloc < hmin)
      hloc = hmin;

    // cout << "restrict h along " << p1 << " - " << p2 << " to " << hloc << endl;
    int i;
    int steps = int (Dist (p1, p2) / hloc) + 2;
    Vec3d v(p1, p2);
  
    for (i = 0; i <= steps; i++)
      {
	Point3d p = p1 + (double(i)/double(steps) * v);
	RestrictLocalH (p, hloc);
      }
  }


  void Mesh :: SetMinimalH (double h)
  {
    hmin = h;
  }


  void Mesh :: SetGlobalH (double h)
  {
    hglob = h;
  }

  double Mesh :: MaxHDomain (int dom) const
  {
    if (maxhdomain.Size())
      return maxhdomain.Get(dom);
    else
      return 1e10;
  }

  void Mesh :: SetMaxHDomain (const ARRAY<double> & mhd)
  {
    maxhdomain.SetSize(mhd.Size());
    for (int i = 1; i <= mhd.Size(); i++)
      maxhdomain.Elem(i) = mhd.Get(i);
  }


  double Mesh :: GetH (const Point3d & p) const
  {
    double hmin = hglob;
    if (lochfunc)
      {
	double hl = lochfunc->GetH (p);
	if (hl < hglob)
	  hmin = hl;
      }
    return hmin;
  }

  double Mesh :: GetMinH (const Point3d & pmin, const Point3d & pmax)
  {
    double hmin = hglob;
    if (lochfunc)
      {
	double hl = lochfunc->GetMinH (pmin, pmax);
	if (hl < hmin)
	  hmin = hl;
      }
    return hmin;
  }





  double Mesh :: AverageH (int surfnr) const
  {
    int i, j, n;
    double hi, hsum;
    double maxh = 0, minh = 1e10;

    hsum = 0;
    n = 0;
    for (i = 1; i <= GetNSE(); i++)
      {
	const Element2d & el = SurfaceElement(i);
	if (surfnr == 0 || el.GetIndex() == surfnr)
	  {
	    for (j = 1; j <= 3; j++)
	      {
		hi = Dist (Point (el.PNumMod(j)), 
			   Point (el.PNumMod(j+1)));

		hsum += hi;

		if (hi > maxh) maxh = hi;
		if (hi < minh) minh = hi;
		n++;
	      }
	  }
      }

    PrintMessage (5, "minh = ", minh, " avh = ", (hsum/n), " maxh = ", maxh);
    return (hsum / n);
  }



  void Mesh :: CalcLocalH () 
  {
    if (!lochfunc)
      {
	Point3d pmin, pmax;
	GetBox (pmin, pmax);
	SetLocalH (pmin, pmax, mparam.grading);
      }

    PrintMessage (3,
		  "CalcLocalH: ", 
		  GetNP(), " Points ", 
		  GetNE(), " Elements ", 
		  GetNSE(), " Surface Elements");


    int i;
    for (i = 0; i < GetNSE(); i++)
      {
	const Element2d & el = surfelements[i];
	int j;

	if (el.GetNP() == 3)
	  {
	    double hel = -1;
	    for (j = 1; j <= 3; j++)
	      {
		const Point3d & p1 = points[el.PNumMod(j)];
		const Point3d & p2 = points[el.PNumMod(j+1)];
	      
		/*
		  INDEX_2 i21(el.PNumMod(j), el.PNumMod(j+1));
		  INDEX_2 i22(el.PNumMod(j+1), el.PNumMod(j));
		  if (! identifiedpoints->Used (i21) &&
		  ! identifiedpoints->Used (i22) )
		*/
		if (!ident -> UsedSymmetric (el.PNumMod(j),
					     el.PNumMod(j+1)))
		  {
		    double hedge = Dist (p1, p2);
		    if (hedge > hel)
		      hel = hedge;
		    //		  lochfunc->SetH (Center (p1, p2), 2 * Dist (p1, p2));
		    //		  (*testout) << "trigseth, p1,2 = " << el.PNumMod(j) << ", " << el.PNumMod(j+1) 
		    //			     << " h = " << (2 * Dist(p1, p2)) << endl;
		  }
	      }
	  
	    if (hel > 0)
	      {
		const Point3d & p1 = points[el.PNum(1)];
		const Point3d & p2 = points[el.PNum(2)];
		const Point3d & p3 = points[el.PNum(3)];
		lochfunc->SetH (Center (p1, p2, p3), hel);
	      }
	  }
	else
	  {
	    {
	      const Point3d & p1 = points[el.PNum(1)];
	      const Point3d & p2 = points[el.PNum(2)];
	      lochfunc->SetH (Center (p1, p2), 2 * Dist (p1, p2));
	    }
	    {
	      const Point3d & p1 = points[el.PNum(3)];
	      const Point3d & p2 = points[el.PNum(4)];
	      lochfunc->SetH (Center (p1, p2), 2 * Dist (p1, p2));
	    }
	  }
      }

    for (i = 0; i < GetNSeg(); i++)
      {
	const Segment & seg = segments[i];
	const Point3d & p1 = points[seg.p1];
	const Point3d & p2 = points[seg.p2];
	/*
	  INDEX_2 i21(seg.p1, seg.p2);
	  INDEX_2 i22(seg.p2, seg.p1);
	  if (identifiedpoints)
	  if (!identifiedpoints->Used (i21) && !identifiedpoints->Used (i22))
	*/
	if (!ident -> UsedSymmetric (seg.p1, seg.p2))
	  {
	    lochfunc->SetH (Center (p1, p2), Dist (p1, p2));
	  }
      }
    /*
      cerr << "do vol" << endl;
      for (i = 1; i <= GetNE(); i++)
      {
      const Element & el = VolumeElement(i);
      if (el.GetType() == TET)
      {
      int j, k;
      for (j = 2; j <= 4; j++)
      for (k = 1; k < j; k++)  
      {
      const Point3d & p1 = Point (el.PNum(j));
      const Point3d & p2 = Point (el.PNum(k));
      lochfunc->SetH (Center (p1, p2), 2 * Dist (p1, p2));
      (*testout) << "set vol h to " << (2 * Dist (p1, p2)) << endl;

      }
      }
      }
    */

    /*
      const char * meshsizefilename = 
      globflags.GetStringFlag ("meshsize", NULL);
      if (meshsizefilename)
      {
      ifstream msf(meshsizefilename);
      if (msf)
      {
      int nmsp;
      msf >> nmsp;
      for (i = 1; i <= nmsp; i++)
      {
      Point3d pi;
      double hi;
      msf >> pi.X() >> pi.Y() >> pi.Z();
      msf >> hi;
      lochfunc->SetH (pi, hi);
      }
      }
      }
    */
    //  lochfunc -> Convexify();
    //  lochfunc -> PrintMemInfo (cout);
  }


  void Mesh :: CalcLocalHFromPointDistances(void)
  {
    PrintMessage (3, "Calculating local h from point distances");
  
    if (!lochfunc)
      {
	Point3d pmin, pmax;
	GetBox (pmin, pmax);
      
	SetLocalH (pmin, pmax, mparam.grading);
      }

    PointIndex i,j;
    double hl;

  
    for (i = PointIndex::BASE; 
	 i < GetNP()+PointIndex::BASE; i++)
      {
	for(j=i+1; j<GetNP()+PointIndex::BASE; j++)
	  {
	    const Point3d & p1 = points[i];
	    const Point3d & p2 = points[j];
	    hl = Dist(p1,p2);
	    RestrictLocalH(p1,hl);
	    RestrictLocalH(p2,hl);
	    //cout << "restricted h at " << p1 << " and " << p2 << " to " << hl << endl;
	  }
      }
  

  }


  void Mesh :: CalcLocalHFromSurfaceCurvature (double elperr) 
  {
    PrintMessage (3, "Calculating local h from surface curvature");

    if (!lochfunc)
      {
	Point3d pmin, pmax;
	GetBox (pmin, pmax);
      
	SetLocalH (pmin, pmax, mparam.grading);
      }

  
    INDEX_2_HASHTABLE<int> edges(3 * GetNP() + 2);
    INDEX_2_HASHTABLE<int> bedges(GetNSeg() + 2);
    int i, j;

    for (i = 1; i <= GetNSeg(); i++)
      {
	const Segment & seg = LineSegment(i);
	INDEX_2 i2(seg.p1, seg.p2);
	i2.Sort();
	bedges.Set (i2, 1);
      }
    for (i = 1; i <= GetNSE(); i++)
      {
	const Element2d & sel = SurfaceElement(i);
	if (!sel.PNum(1))
	  continue;
	for (j = 1; j <= 3; j++)
	  {
	    INDEX_2 i2(sel.PNumMod(j), sel.PNumMod(j+1));
	    i2.Sort();
	    if (bedges.Used(i2)) continue;

	    if (edges.Used(i2))
	      {
		int other = edges.Get(i2);

		const Element2d & elother = SurfaceElement(other);

		int pi3 = 1;
		while ( (sel.PNum(pi3) == i2.I1()) || 
			(sel.PNum(pi3) == i2.I2()))
		  pi3++;
		pi3 = sel.PNum(pi3);

		int pi4 = 1;
		while ( (elother.PNum(pi4) == i2.I1()) || 
			(elother.PNum(pi4) == i2.I2()))
		  pi4++;
		pi4 = elother.PNum(pi4);

		double rad = ComputeCylinderRadius (Point (i2.I1()),
						    Point (i2.I2()),
						    Point (pi3), 
						    Point (pi4));
	      
		RestrictLocalHLine (Point(i2.I1()), Point(i2.I2()), rad/elperr);


		/*	      
		  (*testout) << "pi1,2, 3, 4 = " << i2.I1() << ", " << i2.I2() << ", " << pi3 << ", " << pi4
		  << " p1 = " << Point(i2.I1()) 
		  << ", p2 = " << Point(i2.I2()) 
		  //			 << ", p3 = " << Point(pi3) 
		  //			 << ", p4 = " << Point(pi4) 
		  << ", rad = " << rad << endl;
		*/
	      }
	    else
	      edges.Set (i2, i);
	  }
      }


    // Restrict h due to line segments

    for (i = 1; i <= GetNSeg(); i++)
      {
	const Segment & seg = LineSegment(i);
	const Point3d & p1 = Point(seg.p1);
	const Point3d & p2 = Point(seg.p2);
	RestrictLocalH (Center (p1, p2),  Dist (p1, p2));
      }



    /*


    int i, j;
    int np = GetNP();
    int nseg = GetNSeg();
    int nse = GetNSE();
  
    ARRAY<Vec3d> normals(np);
    BitArray linepoint(np);

    linepoint.Clear();
    for (i = 1; i <= nseg; i++)
    {
    linepoint.Set (LineSegment(i).p1);
    linepoint.Set (LineSegment(i).p2);
    }

    for (i = 1; i <= np; i++)
    normals.Elem(i) = Vec3d(0,0,0);

    for (i = 1; i <= nse; i++)
    {
    Element2d & el = SurfaceElement(i);
    Vec3d nf = Cross (Vec3d (Point (el.PNum(1)), Point(el.PNum(2))),
    Vec3d (Point (el.PNum(1)), Point(el.PNum(3))));
    for (j = 1; j <= 3; j++)
    normals.Elem(el.PNum(j)) += nf;
    }

    for (i = 1; i <= np; i++)
    normals.Elem(i) /= (1e-12 + normals.Elem(i).Length());

    for (i = 1; i <= nse; i++)
    {
    Element2d & el = SurfaceElement(i);
    Vec3d nf = Cross (Vec3d (Point (el.PNum(1)), Point(el.PNum(2))),
    Vec3d (Point (el.PNum(1)), Point(el.PNum(3))));
    nf /= nf.Length();
    Point3d c = Center (Point(el.PNum(1)),
    Point(el.PNum(2)),
    Point(el.PNum(3)));
			  
    for (j = 1; j <= 3; j++)
    {
    if (!linepoint.Test (el.PNum(j)))
    {
    double dist = Dist (c, Point(el.PNum(j)));
    double dn = (nf - normals.Get(el.PNum(j))).Length();
	  
    RestrictLocalH (Point(el.PNum(j)), dist / (dn+1e-12) /elperr);
    }
    }
    }
    */
  }


  void Mesh :: RestrictLocalH (resthtype rht, int nr, double loch)
  {
    int i;
    switch (rht)
      {
      case RESTRICTH_FACE:
	{
	  for (i = 1; i <= GetNSE(); i++)
	    {
	      const Element2d & sel = SurfaceElement(i);
	      if (sel.GetIndex() == nr)
		RestrictLocalH (RESTRICTH_SURFACEELEMENT, i, loch);
	    }
	  break;
	}
      case RESTRICTH_EDGE:
	{
	  for (i = 1; i <= GetNSeg(); i++)
	    {
	      const Segment & seg = LineSegment(i);
	      if (seg.edgenr == nr)
		RestrictLocalH (RESTRICTH_SEGMENT, i, loch);
	    }
	  break;
	}
      case RESTRICTH_POINT:
	{
	  RestrictLocalH (Point (nr), loch);
	  break;
	}

      case RESTRICTH_SURFACEELEMENT:
	{
	  const Element2d & sel = SurfaceElement(nr);
	  Point3d p = Center (Point(sel.PNum(1)),
			      Point(sel.PNum(2)),
			      Point(sel.PNum(3)));
	  RestrictLocalH (p, loch);
	  break;
	}
      case RESTRICTH_SEGMENT:
	{
	  const Segment & seg = LineSegment(nr);
	  RestrictLocalHLine (Point (seg.p1), Point(seg.p2), loch);
	  break;
	}
      }
  }


  void Mesh :: LoadLocalMeshSize (const char * meshsizefilename)
  {
    if (!meshsizefilename) return;

    ifstream msf(meshsizefilename);

    if (!msf) return;

    PrintMessage (3, "Load local mesh-size");
    int nmsp, nmsl;
    msf >> nmsp;
    for (int i = 0; i < nmsp; i++)
      {
	Point3d pi;
	double hi;
	msf >> pi.X() >> pi.Y() >> pi.Z();
	msf >> hi;
	if (!msf.good())
	  throw NgException ("problem in mesh-size file\n");
	RestrictLocalH (pi, hi);
      }
    msf >> nmsl;
    for (int i = 0; i < nmsl; i++)
      {
	Point3d p1, p2;
	double hi;
	msf >> p1.X() >> p1.Y() >> p1.Z();
	msf >> p2.X() >> p2.Y() >> p2.Z();
	msf >> hi;
	if (!msf.good())
	  throw NgException ("problem in mesh-size file\n");
	RestrictLocalHLine (p1, p2, hi);
      }  
  }



  void Mesh :: GetBox (Point3d & pmin, Point3d & pmax, int dom) const
  {
    if (points.Size() == 0)
      {
	pmin = pmax = Point3d(0,0,0);
	return;
      }

    if (dom <= 0)
      {
	pmin = Point3d (1e10, 1e10, 1e10);
	pmax = Point3d (-1e10, -1e10, -1e10); 

	for (PointIndex pi = PointIndex::BASE; 
	     pi < GetNP()+PointIndex::BASE; pi++)
	  {
	    pmin.SetToMin ( (*this) [pi] );
	    pmax.SetToMax ( (*this) [pi] );
	  }
      }
    else
      {
	int j, nse = GetNSE();
	SurfaceElementIndex sei;

	pmin = Point3d (1e10, 1e10, 1e10);
	pmax = Point3d (-1e10, -1e10, -1e10); 
	for (sei = 0; sei < nse; sei++)
	  {
	    const Element2d & el = (*this)[sei];
	    if (el.IsDeleted() ) continue;

	    if (dom == -1 || el.GetIndex() == dom)
	      {
		for (j = 0; j < 3; j++)
		  {
		    pmin.SetToMin ( (*this) [el[j]] );
		    pmax.SetToMax ( (*this) [el[j]] );
		  }
	      }
	  }
      }

    if (pmin.X() > 0.5e10)
      {
	pmin = pmax = Point3d(0,0,0);
      }
  }




  void Mesh :: GetBox (Point3d & pmin, Point3d & pmax, POINTTYPE ptyp) const
  {
    if (points.Size() == 0)
      {
	pmin = pmax = Point3d(0,0,0);
	return;
      }

    pmin = Point3d (1e10, 1e10, 1e10);
    pmax = Point3d (-1e10, -1e10, -1e10); 
  
    for (PointIndex pi = PointIndex::BASE; 
	 pi < GetNP()+PointIndex::BASE; pi++)
      if (points[pi].Type() <= ptyp)
	{
	  pmin.SetToMin ( (*this) [pi] );
	  pmax.SetToMax ( (*this) [pi] );
	}
  }




  double Mesh :: ElementError (int eli) const
  {
    const Element & el = volelements.Get(eli);
    return CalcTetBadness (points.Get(el[0]), points.Get(el[1]),
			   points.Get(el[2]), points.Get(el[3]), -1);
  }

  void Mesh :: AddLockedPoint (PointIndex pi)
  { 
    lockedpoints.Append (pi); 
  }

  void Mesh :: ClearLockedPoints ()
  { 
    lockedpoints.SetSize (0); 
  }



  void Mesh :: Compress ()
  {
    int i, j;
    ARRAY<int,PointIndex::BASE> op2np(GetNP());
    ARRAY<MeshPoint> hpoints;
    BitArrayChar<PointIndex::BASE> pused(GetNP());

    /*
      (*testout) << "volels: " << endl;
      for (i = 1; i <= volelements.Size(); i++)
      {
      for (j = 1; j <= volelements.Get(i).GetNP(); j++)
      (*testout) << volelements.Get(i).PNum(j) << " ";
      (*testout) << endl;
      }
      (*testout) << "np: " << GetNP() << endl;
    */

    for (i = 0; i < volelements.Size(); i++)
      if (volelements[i][0] <= PointIndex::BASE-1 ||
	  volelements[i].IsDeleted())
	{
	  volelements.Delete(i);
	  i--;
	}


    for (i = 0; i < surfelements.Size(); i++)
      if (surfelements[i].IsDeleted())
	{
	  surfelements.Delete(i);
	  i--;
	}

    for (i = 0; i < segments.Size(); i++)
      if (segments[i].p1 <= PointIndex::BASE-1)
	{
	  segments.Delete(i);
	  i--;
	}

    pused.Clear();
    for (i = 0; i < volelements.Size(); i++)
      {
	const Element & el = volelements[i];
	for (j = 0; j < el.GetNP(); j++)
	  pused.Set (el[j]);
      }

    for (i = 0; i < surfelements.Size(); i++)
      {
	const Element2d & el = surfelements[i];
	for (j = 0; j < el.GetNP(); j++)
	  pused.Set (el[j]);
      }

    for (i = 0; i < segments.Size(); i++)
      {
	const Segment & seg = segments[i];
	pused.Set (seg.p1);
	pused.Set (seg.p2);
      }

    for (i = 0; i < openelements.Size(); i++)
      {
	const Element2d & el = openelements[i];
	for (j = 0; j < el.GetNP(); j++)
	  pused.Set(el[j]);
      }

    for (i = 0; i < lockedpoints.Size(); i++)
      pused.Set (lockedpoints[i]);


    /*
    // compress points doesnt work for identified points !
    if (identifiedpoints)
    {
    for (i = 1; i <= identifiedpoints->GetNBags(); i++)
    if (identifiedpoints->GetBagSize(i))
    {
    pused.Set ();
    break;
    }
    }
    */
    //  pused.Set();


    int npi = PointIndex::BASE-1;

    for (i = PointIndex::BASE; 
	 i < points.Size()+PointIndex::BASE; i++)
      if (pused.Test(i))
	{
	  npi++;
	  op2np[i] = npi;
	  hpoints.Append (points[i]);
	}
      else
	op2np[i] = -1;



    points.SetSize(0);
    for (i = 0; i < hpoints.Size(); i++)
      points.Append (hpoints[i]);


    for (i = 1; i <= volelements.Size(); i++)
      {
	Element & el = VolumeElement(i);
	for (j = 0; j < el.GetNP(); j++)
	  el[j] = op2np[el[j]];
      }

    for (i = 1; i <= surfelements.Size(); i++)
      {
	Element2d & el = SurfaceElement(i);
	for (j = 0; j < el.GetNP(); j++)
	  el[j] = op2np[el[j]];
      }
  
    for (i = 0; i < segments.Size(); i++)
      {
	Segment & seg = segments[i];
	seg.p1 = op2np[seg.p1];
	seg.p2 = op2np[seg.p2];
      }

    for (i = 1; i <= openelements.Size(); i++)
      {
	Element2d & el = openelements.Elem(i);
	for (j = 0; j < el.GetNP(); j++)
	  el[j] = op2np[el[j]];
      }  


    for (i = 0; i < lockedpoints.Size(); i++)
      lockedpoints[i] = op2np[lockedpoints[i]];

    for (int i = 0; i < facedecoding.Size(); i++)
      facedecoding[i].firstelement = -1;
    for (int i = surfelements.Size()-1; i >= 0; i--)
      {
        int ind = surfelements[i].GetIndex();
        surfelements[i].next = facedecoding[ind-1].firstelement;
        facedecoding[ind-1].firstelement = i;
      }
    

    CalcSurfacesOfNode();


    //  FindOpenElements();
    timestamp = NextTimeStamp();

    /*
      (*testout) << "compress, done" << endl
      << "np = " << points.Size()
      << "ne = " << volelements.Size() << ", type.size = " << eltyps.Size()
      <<  "volelements = " << volelements << endl;
    */
  }


  int Mesh :: CheckConsistentBoundary () const
  {
    int nf = GetNOpenElements();
    INDEX_2_HASHTABLE<int> edges(nf+2);
    INDEX_2 i2, i2s, edge;
    int err = 0;

    for (int i = 1; i <= nf; i++)
      {
	const Element2d & sel = OpenElement(i);
        
	for (int j = 1; j <= sel.GetNP(); j++)
	  {
	    i2.I1() = sel.PNumMod(j);
	    i2.I2() = sel.PNumMod(j+1);

	    int sign = (i2.I2() > i2.I1()) ? 1 : -1;
	    i2.Sort();
	    if (!edges.Used (i2))
	      edges.Set (i2, 0);
	    edges.Set (i2, edges.Get(i2) + sign);
	  }
      }

    for (int i = 1; i <= edges.GetNBags(); i++)
      for (int j = 1; j <= edges.GetBagSize(i); j++)
	{
	  int cnt = 0;
	  edges.GetData (i, j, i2, cnt);
	  if (cnt)
	    {
	      PrintError ("Edge ", i2.I1() , " - ", i2.I2(), " multiple times in surface mesh");

	      (*testout) << "Edge " << i2 << " multiple times in surface mesh" << endl;
	      i2s = i2;
	      i2s.Sort();
	      for (int k = 1; k <= nf; k++)
		{
		  const Element2d & sel = OpenElement(k);
		  for (int l = 1; l <= sel.GetNP(); l++)
		    {
		      edge.I1() = sel.PNumMod(l);
		      edge.I2() = sel.PNumMod(l+1);
		      edge.Sort();

		      if (edge == i2s) 
			(*testout) << "edge of element " << sel << endl;
		    }
		}


	      err = 2;
	    }
	}

    return err;
  }



  int Mesh :: CheckOverlappingBoundary () 
  {
    int i, j, k;

    Point3d pmin, pmax;
    GetBox (pmin, pmax);
    Box3dTree setree(pmin, pmax);
    ARRAY<int> inters;
  
    bool overlap = 0;
    bool incons_layers = 0;


    for (i = 1; i <= GetNSE(); i++)
      SurfaceElement(i).badel = 0;


    for (i = 1; i <= GetNSE(); i++)
      {
	const Element2d & tri = SurfaceElement(i);
      
	Point3d tpmin (Point(tri[0]));
	Point3d tpmax (tpmin);

	for (k = 1; k < tri.GetNP(); k++)
	  {
	    tpmin.SetToMin (Point (tri[k]));
	    tpmax.SetToMax (Point (tri[k]));
	  }
	Vec3d diag(tpmin, tpmax);

	tpmax = tpmax + 0.1 * diag;
	tpmin = tpmin - 0.1 * diag;

	setree.Insert (tpmin, tpmax, i);
      }

    for (i = 1; i <= GetNSE(); i++)
      {
	const Element2d & tri = SurfaceElement(i);
      
	Point3d tpmin (Point(tri[0]));
	Point3d tpmax (tpmin);

	for (k = 1; k < tri.GetNP(); k++)
	  {
	    tpmin.SetToMin (Point (tri[k]));
	    tpmax.SetToMax (Point (tri[k]));
	  }

	setree.GetIntersecting (tpmin, tpmax, inters);

	for (j = 1; j <= inters.Size(); j++)
	  {
	    const Element2d & tri2 = SurfaceElement(inters.Get(j));	  

	    if ( (*this)[tri[0]].GetLayer() != (*this)[tri2[0]].GetLayer())
	      continue;

	    if ( (*this)[tri[0]].GetLayer() != (*this)[tri[1]].GetLayer() ||
		 (*this)[tri[0]].GetLayer() != (*this)[tri[2]].GetLayer())
	      {
		incons_layers = 1;
		cout << "inconsistent layers in triangle" << endl;
	      }


	    const netgen::Point<3> *trip1[3], *trip2[3];	  
	    for (k = 1; k <= 3; k++)
	      {
		trip1[k-1] = &Point (tri.PNum(k));
		trip2[k-1] = &Point (tri2.PNum(k));
	      }

	    if (IntersectTriangleTriangle (&trip1[0], &trip2[0]))
	      {
		overlap = 1;
		PrintWarning ("Intersecting elements " 
			      ,i, " and ", inters.Get(j));

		(*testout) << "Intersecting: " << endl;
		(*testout) << "openelement " << i << " with open element " << inters.Get(j) << endl;

		cout << "el1 = " << tri << endl;
		cout << "el2 = " << tri2 << endl;
		cout << "layer1 = " <<  (*this)[tri[0]].GetLayer() << endl;
		cout << "layer2 = " <<  (*this)[tri2[0]].GetLayer() << endl;


		for (k = 1; k <= 3; k++)
		  (*testout) << tri.PNum(k) << "  ";
		(*testout) << endl;
		for (k = 1; k <= 3; k++)
		  (*testout) << tri2.PNum(k) << "  ";
		(*testout) << endl;

		for (k = 0; k <= 2; k++)
		  (*testout) << *trip1[k] << "   ";
		(*testout) << endl;
		for (k = 0; k <= 2; k++)
		  (*testout) << *trip2[k] << "   ";
		(*testout) << endl;
		
		(*testout) << "Face1 = " << GetFaceDescriptor(tri.GetIndex()) << endl;
		(*testout) << "Face1 = " << GetFaceDescriptor(tri2.GetIndex()) << endl;

		/*
		  INDEX_3 i3(tri.PNum(1), tri.PNum(2), tri.PNum(3));
		  i3.Sort();
		  for (k = 1; k <= GetNSE(); k++)
		  {
		  const Element2d & el2 = SurfaceElement(k);
		  INDEX_3 i3b(el2.PNum(1), el2.PNum(2), el2.PNum(3));
		  i3b.Sort();
		  if (i3 == i3b)
		  {
		  SurfaceElement(k).badel = 1;
		  }
		  }
		*/
		SurfaceElement(i).badel = 1;
		SurfaceElement(inters.Get(j)).badel = 1;
	      }
	  }
      }

    // bug 'fix'
    if (incons_layers) overlap = 0;

    return overlap;
  }


  int Mesh :: CheckVolumeMesh () const
  {
    PrintMessage (3, "Checking volume mesh");
  
    int ne = GetNE();
    DenseMatrix dtrans(3,3);
    int i, j;

    PrintMessage (5, "elements: ", ne);
    for (i = 1; i <= ne; i++)
      {
	Element & el = (Element&) VolumeElement(i);
	el.flags.badel = 0;
	int nip = el.GetNIP();
	for (j = 1; j <= nip; j++)
	  {
	    el.GetTransformation (j, Points(), dtrans);
	    double det = dtrans.Det();
	    if (det > 0)
	      {
		PrintError ("Element ", i , " has wrong orientation");
		el.flags.badel = 1;
	      }
	  }
      }

    return 0;
  }


  bool Mesh :: LegalTrig (const Element2d & el) const
  {
    return 1;
    if ( /* hp */ 1)  // needed for old, simple hp-refinement
      { 
	// trigs with 2 or more segments are illegal
	int i;
	int nseg = 0;

	if (!segmentht)
	  {
	    cerr << "no segmentht allocated" << endl;
	    return 0;
	  }

	//      Point3d cp(0.5, 0.5, 0.5);
	for (i = 1; i <= 3; i++)
	  {
	    INDEX_2 i2(el.PNumMod (i), el.PNumMod (i+1));
	    i2.Sort();
	    if (segmentht -> Used (i2))
	      nseg++;
	  }
	if (nseg >= 2) 
	  return 0;
      }
    return 1;
  }




  ///
  bool Mesh :: LegalTet2 (Element & el) const
  {
    // static int timer1 = NgProfiler::CreateTimer ("Legaltet2");


    // Test, whether 4 points have a common surface plus
    // at least 4 edges at the boundary

    if(!boundaryedges)
      const_cast<Mesh *>(this)->BuildBoundaryEdges();

  
    // non-tets are always legal
    if (el.GetType() != TET)
      {
	el.SetLegal (1);
	return 1;
      }

    POINTTYPE pointtype[4];
    for(int i = 0; i < 4; i++)
      pointtype[i] = (*this)[el[i]].Type();

    

    // element has at least 2 inner points ---> legal
    int cnti = 0;
    for (int j = 0; j < 4; j++)
      if ( pointtype[j] == INNERPOINT)
	{
	  cnti++;
	  if (cnti >= 2)
	    {
	      el.SetLegal (1);
	      return 1;
	    }
	}



    // which faces are boundary faces ?
    int bface[4];
    for (int i = 0; i < 4; i++)
      {
	bface[i] = surfelementht->Used (INDEX_3::Sort(el[gftetfacesa[i][0]],
                                                      el[gftetfacesa[i][1]],
                                                      el[gftetfacesa[i][2]]));
      }

    int bedge[4][4];
    int segedge[4][4];
    static const int pi3map[4][4] = { { -1,  2,  1,  1 },
                                      {  2, -1,  0,  0 },
                                      {  1,  0, -1,  0 },
                                      {  1,  0,  0, -1 } };

    static const int pi4map[4][4] = { { -1,  3,  3,  2 },
                                      {  3, -1,  3,  2 },
                                      {  3,  3, -1,  1 },
                                      {  2,  2,  1, -1 } };

    
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < i; j++)
	{
	  bool sege = false, be = false;

          int pos = boundaryedges -> Position(INDEX_2::Sort(el[i], el[j]));
	  if (pos)
	    {
	      be = true;
	      if (boundaryedges -> GetData(pos) == 2)
		sege = true;
	    }

	  segedge[j][i] = segedge[i][j] = sege;
	  bedge[j][i] = bedge[i][j] = be;
	}

    // two boundary faces and no edge is illegal
    for (int i = 0; i < 3; i++)
      for (int j = i+1; j < 4; j++)
	{
	  if (bface[i] && bface[j])
            if (!segedge[pi3map[i][j]][pi4map[i][j]])
              {
                // 2 boundary faces withoud edge in between
                el.SetLegal (0);
                return 0;
              }
	}

    // three boundary edges meeting in a Surface point
    for (int i = 0; i < 4; i++)
      {
	bool alledges = 1;
	if ( pointtype[i] == SURFACEPOINT)
	  {
            bool alledges = 1;
	    for (int j = 0; j < 4; j++)
	      if (j != i && !bedge[i][j])
                {
                  alledges = 0;
                  break;
                }
	    if (alledges)
	      {
		// cout << "tet illegal due to unmarked node" << endl;
		el.SetLegal (0);
		return 0;
	      }
	  }
      }



    for (int fnr = 0; fnr < 4; fnr++)
      if (!bface[fnr])
        for (int i = 0; i < 4; i++)
          if (i != fnr)
            {
              int pi1 = pi3map[i][fnr];
              int pi2 = pi4map[i][fnr];

              if ( pointtype[i] == SURFACEPOINT)
                {
                  // two connected edges on surface, but no face
                  if (bedge[i][pi1] && bedge[i][pi2])
                    {
                      el.SetLegal (0);
                      return 0;
                    }
                }

              if ( pointtype[i] == EDGEPOINT)
                {
                  // connected surface edge and edge edge, but no face
                  if (bedge[i][pi1] && segedge[i][pi2] ||
                      bedge[i][pi2] && segedge[i][pi1])
                    {
                      el.SetLegal (0);
                      return 0;
                    }
                }

            }


    el.SetLegal (1);
    return 1;
  
  }



  int Mesh :: GetNDomains() const
  {
    int ndom = 0;

    for (int k = 0; k < facedecoding.Size(); k++)
      {
	if (facedecoding[k].DomainIn() > ndom)
	  ndom = facedecoding[k].DomainIn();
	if (facedecoding[k].DomainOut() > ndom)
	  ndom = facedecoding[k].DomainOut();
      }

    return ndom;
  }



  void Mesh :: SurfaceMeshOrientation ()
  {
    int i, j;
    int nse = GetNSE();
  
    BitArray used(nse);
    used.Clear();
    INDEX_2_HASHTABLE<int> edges(nse+1);

    bool haschanged = 0;


    const Element2d & tri = SurfaceElement(1);
    for (j = 1; j <= 3; j++)
      {
	INDEX_2 i2(tri.PNumMod(j), tri.PNumMod(j+1));
	edges.Set (i2, 1);
      }
    used.Set(1);

    bool unused;
    do
      {
	bool changed;
	do
	  {
	    changed = 0;
	    for (i = 1; i <= nse; i++)
	      if (!used.Test(i))
		{
		  Element2d & el = surfelements.Elem(i);
		  int found = 0, foundrev = 0;
		  for (j = 1; j <= 3; j++)
		    {
		      INDEX_2 i2(el.PNumMod(j), el.PNumMod(j+1));
		      if (edges.Used(i2))
			foundrev = 1;
		      swap (i2.I1(), i2.I2());
		      if (edges.Used(i2))
			found = 1;
		    }
		
		  if (found || foundrev)
		    {
		      if (foundrev)
			swap (el.PNum(2), el.PNum(3));
		    
		      changed = 1;
		      for (j = 1; j <= 3; j++)
			{
			  INDEX_2 i2(el.PNumMod(j), el.PNumMod(j+1));
			  edges.Set (i2, 1);
			}
		      used.Set (i);
		    }
		}
	    if (changed)
	      haschanged = 1;
	  }
	while (changed);
      

	unused = 0;
	for (i = 1; i <= nse; i++)
	  if (!used.Test(i))
	    {
	      unused = 1;
	      const Element2d & tri = SurfaceElement(i);
	      for (j = 1; j <= 3; j++)
		{
		  INDEX_2 i2(tri.PNumMod(j), tri.PNumMod(j+1));
		  edges.Set (i2, 1);
		}
	      used.Set(i);
	      break;
	    }
      }
    while (unused);

    if (haschanged)
      timestamp = NextTimeStamp();
  }
  

  void Mesh :: Split2Tets()
  {
    PrintMessage (1, "Split To Tets");
    bool has_prisms = 0;

    int oldne = GetNE(); 
    for (int i = 1; i <= oldne; i++)
      {
	Element el = VolumeElement(i);

	if (el.GetType() == PRISM)
	  {
	    // prism, to 3 tets

	    // make minimal node to node 1
	    int minpi=0;
	    PointIndex minpnum;
	    minpnum = GetNP() + 1;

	    for (int j = 1; j <= 6; j++)
	      {
		if (el.PNum(j) < minpnum)
		  {
		    minpnum = el.PNum(j);
		    minpi = j;
		  }
	      }

	    if (minpi >= 4)
	      {
		for (int j = 1; j <= 3; j++)
		  swap (el.PNum(j), el.PNum(j+3));
		minpi -= 3;
	      }

	    while (minpi > 1)
	      {
		int hi = 0;
		for (int j = 0; j <= 3; j+= 3)
		  {
		    hi = el.PNum(1+j);
		    el.PNum(1+j) = el.PNum(2+j);
		    el.PNum(2+j) = el.PNum(3+j);
		    el.PNum(3+j) = hi;
		  }
		minpi--;
	      }

	    /*
	      version 1: edge from pi2 to pi6,
	      version 2: edge from pi3 to pi5,
	    */

	    static const int ntets[2][12] =
	      { { 1, 4, 5, 6, 1, 2, 3, 6, 1, 2, 5, 6 },
		{ 1, 4, 5, 6, 1, 2, 3, 5, 3, 1, 5, 6 } };

	    const int * min2pi;

	    if (min2 (el.PNum(2), el.PNum(6)) <
		min2 (el.PNum(3), el.PNum(5)))
	      {
		min2pi = &ntets[0][0];
		// (*testout) << "version 1 ";
	      }
	    else
	      {
		min2pi = &ntets[1][0];
		// (*testout) << "version 2 ";
	      }

	  
	    int firsttet = 1;
	    for (int j = 1; j <= 3; j++)
	      {
		Element nel(TET);
		for (int k = 1; k <= 4; k++)
		  nel.PNum(k) = el.PNum(min2pi[4 * j + k - 5]);
		nel.SetIndex (el.GetIndex());

		int legal = 1;
		for (int k = 1; k <= 3; k++)
		  for (int l = k+1; l <= 4; l++)
		    if (nel.PNum(k) == nel.PNum(l))
		      legal = 0;

		// (*testout) << nel << " ";
		if (legal)
		  {
		    if (firsttet)
		      {
			VolumeElement(i) = nel;
			firsttet = 0;
		      }
		    else
		      {
			AddVolumeElement(nel);
		      }
		  }
	      }
	    if (firsttet) cout << "no legal";
	    (*testout) << endl;
	  }
      


	else if (el.GetType() == HEX)
	  {
	    // hex to A) 2 prisms or B) to 5 tets

	    // make minimal node to node 1
	    int minpi=0;
	    PointIndex minpnum;
	    minpnum = GetNP() + 1;

	    for (int j = 1; j <= 8; j++)
	      {
		if (el.PNum(j) < minpnum)
		  {
		    minpnum = el.PNum(j);
		    minpi = j;
		  }
	      }

	    if (minpi >= 5)
	      {
		for (int j = 1; j <= 4; j++)
		  swap (el.PNum(j), el.PNum(j+4));
		minpi -= 4;
	      }

	    while (minpi > 1)
	      {
		int hi = 0;
		for (int j = 0; j <= 4; j+= 4)
		  {
		    hi = el.PNum(1+j);
		    el.PNum(1+j) = el.PNum(2+j);
		    el.PNum(2+j) = el.PNum(3+j);
		    el.PNum(3+j) = el.PNum(4+j);
		    el.PNum(4+j) = hi;
		  }
		minpi--;
	      }



	    static const int to_prisms[3][12] =
	      { { 0, 1, 2, 4, 5, 6, 0, 2, 3, 4, 6, 7 },
		{ 0, 1, 5, 3, 2, 6, 0, 5, 4, 3, 6, 7 },
		{ 0, 7, 4, 1, 6, 5, 0, 3, 7, 1, 2, 6 },
	      };

	    const int * min2pi = 0;
	    if (min2 (el[4], el[6]) < min2 (el[5], el[7]))
	      min2pi = &to_prisms[0][0];
	    else if (min2 (el[3], el[6]) < min2 (el[2], el[7]))
	      min2pi = &to_prisms[1][0];
	    else if (min2 (el[1], el[6]) < min2 (el[2], el[5]))
	      min2pi = &to_prisms[2][0];

	    if (min2pi)
	      {
		has_prisms = 1;
		for (int j = 0; j < 2; j++)
		  {
		    Element nel(PRISM);
		    for (int k = 0; k < 6; k++)
		      nel[k] = el[min2pi[6*j + k]];
		    nel.SetIndex (el.GetIndex());
		  
		    if (j == 0)
		      VolumeElement(i) = nel;
		    else
		      AddVolumeElement(nel);
		  }
	      }
	    else
	      {
		// split to 5 tets
	      
		static const int to_tets[20] =
		  {
		    1, 2, 0, 5,
		    3, 0, 2, 7,
		    4, 5, 7, 0,
		    6, 7, 5, 2,
		    0, 2, 7, 5
		  };

		for (int j = 0; j < 5; j++)
		  {
		    Element nel(TET);
		    for (int k = 0; k < 4; k++)
		      nel[k] = el[to_tets[4*j + k]];
		    nel.SetIndex (el.GetIndex());
		  
		    if (j == 0)
		      VolumeElement(i) = nel;
		    else
		      AddVolumeElement(nel);
		  }
	      
	      }
	  }
      



      
	else if (el.GetType() == PYRAMID)
	  {
	    // pyramid, to 2 tets
	  
	    // cout << "pyramid: " << el << endl;
	    
	    static const int ntets[2][8] =
	      { { 1, 2, 3, 5, 1, 3, 4, 5 },
		{ 1, 2, 4, 5, 4, 2, 3, 5 }};

	    const int * min2pi;

	    if (min2 (el[0], el[2]) < min2 (el[1], el[3]))
	      min2pi = &ntets[0][0];
	    else
	      min2pi = &ntets[1][0];

	    bool firsttet = 1;
	    for (int j = 0; j < 2; j++)
	      {
		Element nel(TET);
		for (int k = 0; k < 4; k++)
		  nel[k] = el[min2pi[4*j + k]-1];
		nel.SetIndex (el.GetIndex());

		// cout << "pyramid-tet: " << nel << endl;

		bool legal = 1;
		for (int k = 0; k < 3; k++)
		  for (int l = k+1; l < 4; l++)
		    if (nel[k] == nel[l])
		      legal = 0;

		if (legal)
		  {
		    (*testout) << nel << " ";
		    if (firsttet)
		      VolumeElement(i) = nel;
		    else
		      AddVolumeElement(nel);

		    firsttet = 0;
		  }
	      }
	    if (firsttet) cout << "no legal";
	    (*testout) << endl;
	  }
      }

  
    int oldnse = GetNSE(); 
    for (int i = 1; i <= oldnse; i++)
      {
	Element2d el = SurfaceElement(i);
	if (el.GetNP() == 4)
	  {
	    (*testout) << "split el: " << el << " to ";
	  
	    static const int ntris[2][6] =
	      { { 1, 2, 3, 1, 3, 4 },
		{ 1, 2, 4, 4, 2, 3 }};

	    const int * min2pi;

	    if (min2 (el.PNum(1), el.PNum(3)) <
		min2 (el.PNum(2), el.PNum(4)))
	      min2pi = &ntris[0][0];
	    else
	      min2pi = &ntris[1][0];

	    for (int j = 0; j <6; j++)
	      (*testout) << min2pi[j] << " ";


	    int firsttri = 1;
	    for (int j = 1; j <= 2; j++)
	      {
		Element2d nel(3);
		for (int k = 1; k <= 3; k++)
		  nel.PNum(k) = el.PNum(min2pi[3 * j + k - 4]);
		nel.SetIndex (el.GetIndex());

		int legal = 1;
		for (int k = 1; k <= 2; k++)
		  for (int l = k+1; l <= 3; l++)
		    if (nel.PNum(k) == nel.PNum(l))
		      legal = 0;

		if (legal)
		  {
		    (*testout) << nel << " ";
		    if (firsttri)
		      {
			SurfaceElement(i) = nel;
			firsttri = 0;
		      }
		    else
		      {
			AddSurfaceElement(nel);
		      }
		  }
	      }
	    (*testout) << endl;

	  }
      }


    if (has_prisms)

      Split2Tets();
  
    else
      {
	for (int i = 1; i <= GetNE(); i++)
	  {
	    Element & el = VolumeElement(i);
	    const Point3d & p1 = Point (el.PNum(1));
	    const Point3d & p2 = Point (el.PNum(2));
	    const Point3d & p3 = Point (el.PNum(3));
	    const Point3d & p4 = Point (el.PNum(4));
	  
	    double vol = (Vec3d (p1, p2) * 
			  Cross (Vec3d (p1, p3), Vec3d(p1, p4)));
	    if (vol > 0)
	      swap (el.PNum(3), el.PNum(4));
	  }



	UpdateTopology();
	timestamp = NextTimeStamp();
      }
  }

  void Mesh :: BuildElementSearchTree ()
  {
    if (elementsearchtreets == GetTimeStamp())
      return;

    NgLock lock(mutex);
    lock.Lock();

    PrintMessage (4, "Rebuild element searchtree");

    if (elementsearchtree)
      delete elementsearchtree;
    elementsearchtree = NULL;

    Box3d box;
    int i, j;
    int ne = GetNE();
    if (!ne) 
      {
	lock.UnLock();
	return;
      }

    box.SetPoint (Point (VolumeElement(1).PNum(1)));
    for (i = 1; i <= ne; i++)
      {
	const Element & el = VolumeElement(i);
	for (j = 1; j <= el.GetNP(); j++)
	  box.AddPoint (Point (el.PNum(j)));
      }
  
    box.Increase (1.01 * box.CalcDiam());
    elementsearchtree = new Box3dTree (box.PMin(), box.PMax());
  


    for (i = 1; i <= ne; i++)
      {
	const Element & el = VolumeElement(i);
	box.SetPoint (Point (el.PNum(1)));
	for (j = 1; j <= el.GetNP(); j++)
	  box.AddPoint (Point (el.PNum(j)));

	elementsearchtree -> Insert (box.PMin(), box.PMax(), i);
      }

    elementsearchtreets = GetTimeStamp();

    lock.UnLock();
  }


  
  bool Mesh :: PointContainedIn2DElement(const Point3d & p,
					 double lami[3],
					 const int element,
					 bool consider3D) const
  {
    static Vec3d col1, col2, col3;
    static Vec3d rhs, sol;
    const double eps = 1e-6;

    static ARRAY<Element2d> loctrigs;


    //SZ 
    if(SurfaceElement(element).GetType()==QUAD)
      {
	const Element2d & el = SurfaceElement(element); 
	      
	const Point3d & p1 = Point(el.PNum(1)); 
	const Point3d & p2 = Point(el.PNum(2));
	const Point3d & p3 = Point(el.PNum(3));
	const Point3d & p4 = Point(el.PNum(4)); 

	// Coefficients of Bilinear Mapping from Ref-Elem to global Elem
	// X = a + b x + c y + d x y 
	Vec3d a = p1; 
	Vec3d b = p2 - a; 
	Vec3d c = p4 - a; 
	Vec3d d = p3 - a - b - c; 
	   
	double dxb = d.X()*b.Y()-d.Y()*b.X();
	double dxc = d.X()*c.Y()-d.Y()*c.X(); 
	double dxa = d.X()*a.Y()-d.Y()*a.X(); 
	double dxp = d.X()*p.Y()-d.Y()*p.X(); 
	      	      
	double c0,c1,c2,rt; 
	lami[2]=0.; 
	double eps = 1.E-12; 

	if(fabs(d.X()) <= eps && fabs(d.Y())<= eps)
	  {
	    //Solve Linear System
	    lami[0]=(c.Y()*(p.X()-a.X())-c.X()*(p.Y()-a.Y()))/
	      (b.X()*c.Y() -b.Y()*c.X()); 
	    lami[1]=(-b.Y()*(p.X()-a.X())+b.X()*(p.Y()-a.Y()))/
	      (b.X()*c.Y() -b.Y()*c.X()); 
	  } 
	else
	  if(fabs(dxb) <= eps) 
	    {
	      lami[1] = (dxp-dxa)/dxc;
	      if(fabs(b.X()-d.X()*lami[1])>=eps)
		lami[0] = (p.X()-a.X() - c.X()*lami[1])/(b.X()+d.X()*lami[1]); 
	      else
		lami[0] = (p.Y()-a.Y() - c.Y()*lami[1])/(b.Y()+d.Y()*lami[1]); 
	    }
	  else
	    if(fabs(dxc) <= eps)
	      {
		lami[0] = (dxp-dxa)/dxb;
		if(fabs(c.X()-d.X()*lami[0])>=eps)
		  lami[1] = (p.X()-a.X() - b.X()*lami[0])/(c.X()+d.X()*lami[0]); 
		else
		  lami[1] = (p.Y()-a.Y() - b.Y()*lami[0])/(c.Y()+d.Y()*lami[0]); 
	      }
	    else //Solve quadratic equation
	      {
		if(fabs(d.X()) >= eps)
		  {
		    c2 = d.X()*dxc;
		    c1 = d.X()*dxc - c.X()*dxb - d.X()*(dxp-dxa);
		    c0 = -b.X()*(dxp -dxa) - (a.X()-p.X())*dxb;
		  }
		else 
		  {
		    c2 = d.Y()*dxc;
		    c1 = d.Y()*dxc - c.Y()*dxb - d.Y()*(dxp-dxa);
		    c0 = -b.Y()*(dxp -dxa) - (a.Y()-p.Y())*dxb;
		  }

		double rt =  c1*c1 - 4*c2*c0;
		if (rt < 0.) return false; 
		lami[1] = (-c1 + sqrt(rt))/2/c2;
		if(lami[1]<=1. && lami[1]>=0.)
		  {
		    lami[0] = (dxp - dxa -dxc*lami[1])/dxb;
		    if(lami[0]<=1. && lami[0]>=0.)
		      return true;
		  }
		      
		lami[1] = (-c1 - sqrt(rt))/2/c2;
		lami[0] = (dxp - dxa -dxc*lami[1])/dxb;
	      }

	if( lami[0] <= 1.+eps  && lami[0] >= -eps && lami[1]<=1.+eps && lami[1]>=-eps)
	  {
	    if(consider3D)
	      {
		Vec3d n = Cross(b,c);
		lami[2] = 0;
		for(int i=1; i<=3; i++)
		  lami[2] +=(p.X(i)-a.X(i)-lami[0]*b.X(i)-lami[1]*c.X(i)) * n.X(i);
		if(lami[2] >= -eps && lami[2] <= eps)
		  return true;
	      }
	    else
	      return true;
	  }
		      
	return false;
	      
      }
    else
      {
	//	  SurfaceElement(element).GetTets (loctets);
	loctrigs.SetSize(1);
	loctrigs.Elem(1) = SurfaceElement(element);
	      
	      
	      
	for (int j = 1; j <= loctrigs.Size(); j++)
	  {
	    const Element2d & el = loctrigs.Get(j);
		  
		  
	    const Point3d & p1 = Point(el.PNum(1));
	    const Point3d & p2 = Point(el.PNum(2));
	    const Point3d & p3 = Point(el.PNum(3));
	    /*
	      Box3d box;
	      box.SetPoint (p1);
	      box.AddPoint (p2);
	      box.AddPoint (p3);
	      box.AddPoint (p4);
	      if (!box.IsIn (p))
	      continue;
	    */
	    col1 = p2-p1;
	    col2 = p3-p1;
	    col3 = Cross(col1,col2);
	    //col3 = Vec3d(0, 0, 1);
	    rhs = p - p1;
		  
	    int retval = SolveLinearSystem (col1, col2, col3, rhs, sol);

	    //(*testout) << "retval " << retval << endl;
		  
	    //(*testout) << "col1 " << col1 << " col2 " << col2 << " col3 " << col3 << " rhs " << rhs << endl;
	    //(*testout) << "sol " << sol << endl;

	    if (sol.X() >= -eps && sol.Y() >= -eps && 
		sol.X() + sol.Y() <= 1+eps)
	      {
		if(!consider3D || (sol.Z() >= -eps && sol.Z() <= eps))
		  {
		    lami[0] = sol.X();
		    lami[1] = sol.Y();
		    lami[2] = sol.Z();
		    
		    return true;
		  }
	      }
	  }
      }

    return false;

  }




  bool Mesh :: PointContainedIn3DElement(const Point3d & p,
					 double lami[3],
					 const int element) const
  {
    //bool oldresult = PointContainedIn3DElementOld(p,lami,element);
    //(*testout) << "old result: " << oldresult
    //       << " lam " << lami[0] << " " << lami[1] << " " << lami[2] << endl;

    //if(!curvedelems->IsElementCurved(element-1))
    //  return PointContainedIn3DElementOld(p,lami,element);


    const double eps = 1.e-4;
    const Element & el = VolumeElement(element);

    netgen::Point<3> lam;

    if (el.GetType() == TET)
      {
	lam = 0.25;
      }
    else if (el.GetType() == PRISM)
      {
	lam(0) = 0.33; lam(1) = 0.33; lam(2) = 0.5;
      }
    else if (el.GetType() == PYRAMID)
      {
	lam(0) = 0.4; lam(1) = 0.4; lam(2) = 0.2;
      }
    else if (el.GetType() == HEX)
      {
	lam = 0.5;
      }
    

    Vec<3> deltalam,rhs;
    netgen::Point<3> x;
    Mat<3,3> Jac,Jact;

    double delta=1;

    bool retval;
    
    int i = 0;

    const int maxits = 30;

    while(delta > 1e-16 && i<maxits)
      {
	curvedelems->CalcElementTransformation(lam,element-1,x,Jac);

	rhs = p-x;
	Jac.Solve(rhs,deltalam);

	lam += deltalam;

	delta = deltalam.Length2();

	i++;
	//(*testout) << "pcie i " << i << " delta " << delta << " p " << p << " x " << x << " lam " << lam << endl;
	//<< "Jac " << Jac << endl;
      }

    if(i==maxits)
      return false;


    for(i=0; i<3; i++)
      lami[i] = lam(i);



    if (el.GetType() == TET)
      {
	retval = (lam(0) > -eps && 
		  lam(1) > -eps && 
		  lam(2) > -eps && 
		  lam(0) + lam(1) + lam(2) < 1+eps);
      }
    else if (el.GetType() == PRISM)
      {
	retval = (lam(0) > -eps &&
		  lam(1) > -eps &&
		  lam(2) > -eps &&
		  lam(2) < 1+eps &&
		  lam(0) + lam(1) < 1+eps);
      }
    else if (el.GetType() == PYRAMID)
      {
	retval = (lam(0) > -eps &&
		  lam(1) > -eps &&
		  lam(2) > -eps &&
		  lam(0) + lam(2) < 1+eps &&
		  lam(1) + lam(2) < 1+eps);
      }
    else if (el.GetType() == HEX)
      {
	retval = (lam(0) > -eps && lam(0) < 1+eps &&
		  lam(1) > -eps && lam(1) < 1+eps &&
		  lam(2) > -eps && lam(2) < 1+eps);
      }
    else
      throw NgException("Da haun i wos vagessn");

    return retval;
  }
  

  
  bool Mesh :: PointContainedIn3DElementOld(const Point3d & p,
					    double lami[3],
					    const int element) const
  {

    static Vec3d col1, col2, col3;
    static Vec3d rhs, sol;
    const double eps = 1.e-4;
    
    static ARRAY<Element> loctets;

    VolumeElement(element).GetTets (loctets);
    
    for (int j = 1; j <= loctets.Size(); j++)
      {
	const Element & el = loctets.Get(j);
	
	const Point3d & p1 = Point(el.PNum(1));
	const Point3d & p2 = Point(el.PNum(2));
	const Point3d & p3 = Point(el.PNum(3));
	const Point3d & p4 = Point(el.PNum(4));
	
	Box3d box;
	box.SetPoint (p1);
	box.AddPoint (p2);
	box.AddPoint (p3);
	box.AddPoint (p4);
	if (!box.IsIn (p))
	  continue;
	
	col1 = p2-p1;
	col2 = p3-p1;
	col3 = p4-p1;
	rhs = p - p1;
	
	SolveLinearSystem (col1, col2, col3, rhs, sol);
	
	if (sol.X() >= -eps && sol.Y() >= -eps && sol.Z() >= -eps &&
	    sol.X() + sol.Y() + sol.Z() <= 1+eps)
	  {
	    ARRAY<Element> loctetsloc;
	    ARRAY<netgen::Point<3> > pointsloc;
	    
	    VolumeElement(element).GetTetsLocal (loctetsloc);
	    VolumeElement(element).GetNodesLocalNew (pointsloc);
	    
	    const Element & le = loctetsloc.Get(j);
	   

	    Point3d pp = 
	      pointsloc.Get(le.PNum(1)) 
	      + sol.X() * Vec3d (pointsloc.Get(le.PNum(1)), pointsloc.Get(le.PNum(2))) 
	      + sol.Y() * Vec3d (pointsloc.Get(le.PNum(1)), pointsloc.Get(le.PNum(3))) 
	      + sol.Z() * Vec3d (pointsloc.Get(le.PNum(1)), pointsloc.Get(le.PNum(4))) ;
	    
	    lami[0] = pp.X();
	    lami[1] = pp.Y();
	    lami[2] = pp.Z();
	    return true;
	  }
      }
    return false;
  }
 

  int Mesh :: GetElementOfPoint (const Point3d & p,
				 double lami[3],
				 bool build_searchtree,
				 const int index,
				 const bool allowindex) const
  {
    if(index != -1) 
      {
	ARRAY<int> dummy(1);
	dummy[0] = index;
	return GetElementOfPoint(p,lami,&dummy,build_searchtree,allowindex);
      }
    else
      return GetElementOfPoint(p,lami,NULL,build_searchtree,allowindex);
  }




  int Mesh :: GetElementOfPoint (const Point3d & p,
				 double lami[3],
				 const ARRAY<int> * const indices,
				 bool build_searchtree,
				 const bool allowindex) const
  {
    if (dimension == 2)
      {
	int i, j;
	int ne;
      
	
	if(ps_startelement != 0 && ps_startelement <= GetNSE() && PointContainedIn2DElement(p,lami,ps_startelement))
	  return ps_startelement;

	ARRAY<int> locels;
	if (0)
	  {
	    elementsearchtree->GetIntersecting (p, p, locels);
	    ne = locels.Size();
	  }
	else
	  ne = GetNSE();

	for (i = 1; i <= ne; i++)
	  {
	    int ii;

	    if (0)
	      ii = locels.Get(i);
	    else
	      ii = i;
	      
	    if(ii == ps_startelement) continue;

	    if(indices != NULL && indices->Size() > 0)
	      {
		bool contained = indices->Contains(SurfaceElement(ii).GetIndex());
		if((allowindex && !contained) || (!allowindex && contained)) continue;
	      }

	    if(PointContainedIn2DElement(p,lami,ii)) return ii;

	  }
	return 0;
      }
    else
      
      {
	int i, j;
	int ne;
      
	if(ps_startelement != 0 && PointContainedIn3DElement(p,lami,ps_startelement))
	  return ps_startelement;

	ARRAY<int> locels;
	if (elementsearchtree || build_searchtree)
	  {
	    // update if necessary:
	    const_cast<Mesh&>(*this).BuildElementSearchTree (); 
	    elementsearchtree->GetIntersecting (p, p, locels);
	    ne = locels.Size();
	  }
	else
	  ne = GetNE();
      
	for (i = 1; i <= ne; i++)
	  {
	    int ii;

	    if (elementsearchtree)
	      ii = locels.Get(i);
	    else
	      ii = i;
	      
	    if(ii == ps_startelement) continue;
	    
	    if(indices != NULL && indices->Size() > 0)
	      {
		bool contained = indices->Contains(VolumeElement(ii).GetIndex());
		if((allowindex && !contained) || (!allowindex && contained)) continue;
	      }
	    
	  
	    if(PointContainedIn3DElement(p,lami,ii)) 
	      {
		ps_startelement = ii;
		return ii;
	      }
	  }

	// Not found, try uncurved variant:
	for (i = 1; i <= ne; i++)
	  {
	    int ii;

	    if (elementsearchtree)
	      ii = locels.Get(i);
	    else
	      ii = i;
	    
	    if(indices != NULL && indices->Size() > 0)
	      {
		bool contained = indices->Contains(VolumeElement(ii).GetIndex());
		if((allowindex && !contained) || (!allowindex && contained)) continue;
	      }
	    
	  
	    if(PointContainedIn3DElementOld(p,lami,ii)) 
	      {
		ps_startelement = ii;
		(*testout) << "WARNING: found element of point " << p <<" only for uncurved mesh" << endl;
		return ii;
	      }
	  }

      
	return 0;
      }
  }



  int Mesh :: GetSurfaceElementOfPoint (const Point3d & p,
					double lami[3],
					bool build_searchtree,
					const int index,
					const bool allowindex) const
  {
    if(index != -1) 
      {
	ARRAY<int> dummy(1);
	dummy[0] = index;
	return GetSurfaceElementOfPoint(p,lami,&dummy,build_searchtree,allowindex);
      }
    else
      return GetSurfaceElementOfPoint(p,lami,NULL,build_searchtree,allowindex);
  }




  int Mesh :: GetSurfaceElementOfPoint (const Point3d & p,
					double lami[3],
					const ARRAY<int> * const indices,
					bool build_searchtree,
					const bool allowindex) const
  {
    if (dimension == 2)
      {
	throw NgException("GetSurfaceElementOfPoint not yet implemented for 2D meshes");
      }
    else
      {
	double vlam[3];
	int velement = GetElementOfPoint(p,vlam,NULL,build_searchtree,allowindex);

	//(*testout) << "p " << p << endl;
	//(*testout) << "velement " << velement << endl;

	ARRAY<int> faces;
	topology->GetElementFaces(velement,faces);

	//(*testout) << "faces " << faces << endl;

	for(int i=0; i<faces.Size(); i++)
	  faces[i] = topology->GetFace2SurfaceElement(faces[i]);

	//(*testout) << "surfel " << faces << endl;

	for(int i=0; i<faces.Size(); i++)
	  {
	    if(faces[i] == 0)
	      continue;

	    if(indices && indices->Size() != 0)
	      {
		if(indices->Contains(SurfaceElement(faces[i]).GetIndex()) &&
		   PointContainedIn2DElement(p,lami,faces[i],true))
		  return faces[i];
	      }
	    else
	      {
		if(PointContainedIn2DElement(p,lami,faces[i],true))
		  {
		    //(*testout) << "found point " << p << " in sel " << faces[i]
		    //	       << ", lam " << lami[0] << ", " << lami[1] << ", " << lami[2] << endl;
		    return faces[i];
		  }
	      }
	  }
	     
      }

    return 0;
  }


  void Mesh::GetIntersectingVolEls(const Point3d& p1, const Point3d& p2, 
				   ARRAY<int> & locels) const
  {
    elementsearchtree->GetIntersecting (p1, p2, locels);
  }

  void Mesh :: SplitIntoParts()
  {
    int i, j, dom;
    int ne = GetNE();
    int np = GetNP();
    int nse = GetNSE();

    BitArray surfused(nse);
    BitArray pused (np);

    surfused.Clear();

    dom = 0;
  
    while (1)
      {
	int cntd = 1;

	dom++;
      
	pused.Clear();

	int found = 0;
	for (i = 1; i <= nse; i++)
	  if (!surfused.Test(i))
	    {
	      SurfaceElement(i).SetIndex (dom);
	      for (j = 1; j <= 3; j++)
		pused.Set (SurfaceElement(i).PNum(j));
	      found = 1;
	      cntd = 1;
	      surfused.Set(i);
	      break;
	    }

	if (!found)
	  break;

	int change;
	do
	  {
	    change = 0;
	    for (i = 1; i <= nse; i++)
	      {
		int is = 0, isnot = 0;
		for (j = 1; j <= 3; j++)
		  if (pused.Test(SurfaceElement(i).PNum(j)))
		    is = 1;
		  else
		    isnot = 1;
	      
		if (is && isnot)
		  {
		    change = 1;
		    for (j = 1; j <= 3; j++)
		      pused.Set (SurfaceElement(i).PNum(j));
		  }

		if (is) 
		  {
		    if (!surfused.Test(i))
		      {
			surfused.Set(i);
			SurfaceElement(i).SetIndex (dom);
			cntd++;
		      }
		  }
	      }


	    for (i = 1; i <= ne; i++)
	      {
		int is = 0, isnot = 0;
		for (j = 1; j <= 4; j++)
		  if (pused.Test(VolumeElement(i).PNum(j)))
		    is = 1;
		  else
		    isnot = 1;
	      
		if (is && isnot)
		  {
		    change = 1;
		    for (j = 1; j <= 4; j++)
		      pused.Set (VolumeElement(i).PNum(j));
		  }

		if (is)
		  {
		    VolumeElement(i).SetIndex (dom);
		  }
	      }
	  }
	while (change);

	PrintMessage (3, "domain ", dom, " has ", cntd, " surfaceelements");
      }

    /*
      facedecoding.SetSize (dom);
      for (i = 1; i <= dom; i++)
      {
      facedecoding.Elem(i).surfnr = 0;
      facedecoding.Elem(i).domin = i;
      facedecoding.Elem(i).domout = 0;
      }
    */
    ClearFaceDescriptors();
    for (i = 1; i <= dom; i++)
      AddFaceDescriptor (FaceDescriptor (0, i, 0, 0));
    CalcSurfacesOfNode();
    timestamp = NextTimeStamp();
  }

  void Mesh :: SplitSeparatedFaces ()
  {
    PrintMessage (3, "SplitSeparateFaces");
    int fdi;
    int np = GetNP();

    BitArray usedp(np);
    ARRAY<SurfaceElementIndex> els_of_face;

    fdi = 1;
    while (fdi <= GetNFD())
      {
	GetSurfaceElementsOfFace (fdi, els_of_face);

	if (els_of_face.Size() == 0) continue;

	SurfaceElementIndex firstel = els_of_face[0];

	usedp.Clear();
	for (int j = 1; j <= SurfaceElement(firstel).GetNP(); j++)
	  usedp.Set (SurfaceElement(firstel).PNum(j));

	bool changed;
	do
	  {
	    changed = false;

	    for (int i = 0; i < els_of_face.Size(); i++)
	      {
		const Element2d & el = SurfaceElement(els_of_face[i]);

		bool has = 0;
		bool hasno = 0;
		for (int j = 0; j < el.GetNP(); j++)
		  {
		    if (usedp.Test(el[j]))
		      has = true;
		    else
		      hasno = true;
		  }

		if (has && hasno)
		  changed = true;

		if (has)
		  for (int j = 0; j < el.GetNP(); j++)
		    usedp.Set (el[j]);
	      }
	  }
	while (changed);

	int nface = 0;
	for (int i = 0; i < els_of_face.Size(); i++)
	  {
	    Element2d & el = SurfaceElement(els_of_face[i]);

	    int hasno = 0;
	    for (int j = 1; j <= el.GetNP(); j++)
	      if (!usedp.Test(el.PNum(j)))
		  hasno = 1;
	  
	    if (hasno)
	      {
		if (!nface)
		  {
		    FaceDescriptor nfd = GetFaceDescriptor(fdi);
		    nface = AddFaceDescriptor (nfd);
		  }

		el.SetIndex (nface);
	      }
	  }

        // reconnect list
        if (nface)
          {
            facedecoding[nface-1].firstelement = -1;
            facedecoding[fdi-1].firstelement = -1;

            for (int i = 0; i < els_of_face.Size(); i++)
              {
                int ind = SurfaceElement(els_of_face[i]).GetIndex();
                SurfaceElement(els_of_face[i]).next = facedecoding[ind-1].firstelement;
                facedecoding[ind-1].firstelement = els_of_face[i];
              }
          }

	fdi++;
      }


    /*
    fdi = 1;
    while (fdi <= GetNFD())
      {
	int firstel = 0;
	for (int i = 1; i <= GetNSE(); i++)
	  if (SurfaceElement(i).GetIndex() == fdi)
	    {
	      firstel = i;
	      break;
	    }
	if (!firstel) continue;

	usedp.Clear();
	for (int j = 1; j <= SurfaceElement(firstel).GetNP(); j++)
	  usedp.Set (SurfaceElement(firstel).PNum(j));

	int changed;
	do
	  {
	    changed = 0;
	    for (int i = 1; i <= GetNSE(); i++)
	      {
		const Element2d & el = SurfaceElement(i);
		if (el.GetIndex() != fdi)
		  continue;

		int has = 0;
		int hasno = 0;
		for (int j = 1; j <= el.GetNP(); j++)
		  {
		    if (usedp.Test(el.PNum(j)))
		      has = 1;
		    else
		      hasno = 1;
		  }
		if (has && hasno)
		  changed = 1;

		if (has)
		  for (int j = 1; j <= el.GetNP(); j++)
		    usedp.Set (el.PNum(j));
	      }
	  }
	while (changed);

	int nface = 0;
	for (int i = 1; i <= GetNSE(); i++)
	  {
	    Element2d & el = SurfaceElement(i);
	    if (el.GetIndex() != fdi)
	      continue;	  

	    int hasno = 0;
	    for (int j = 1; j <= el.GetNP(); j++)
	      {
		if (!usedp.Test(el.PNum(j)))
		  hasno = 1;
	      }
	  
	    if (hasno)
	      {
		if (!nface)
		  {
		    FaceDescriptor nfd = GetFaceDescriptor(fdi);
		    nface = AddFaceDescriptor (nfd);
		  }

		el.SetIndex (nface);
	      }
	  }
	fdi++;
      }
    */
  }


  void Mesh :: GetSurfaceElementsOfFace (int facenr, ARRAY<SurfaceElementIndex> & sei) const
  {
    static int timer = NgProfiler::CreateTimer ("GetSurfaceElementsOfFace");
    NgProfiler::RegionTimer reg (timer);

    /*
    sei.SetSize (0);
    for (SurfaceElementIndex i = 0; i < GetNSE(); i++)
      if ( (*this)[i].GetIndex () == facenr && (*this)[i][0] >= PointIndex::BASE &&
	   !(*this)[i].IsDeleted() )
	sei.Append (i);

    int size1 = sei.Size();
    */

    sei.SetSize(0);

    SurfaceElementIndex si = facedecoding[facenr-1].firstelement;
    while (si != -1)
      {
        if ( (*this)[si].GetIndex () == facenr && (*this)[si][0] >= PointIndex::BASE &&
             !(*this)[si].IsDeleted() )
          {
            sei.Append (si);
          }

        si = (*this)[si].next;
      }

    /*
     // *testout << "with list = " << endl << sei << endl;

    if (size1 != sei.Size()) 
      {
        cout << "size mismatch" << endl;
        exit(1);
      }
    */
  }




  void Mesh :: CalcMinMaxAngle (double badellimit, double * retvalues) 
  {
    int i, j;
    int lpi1, lpi2, lpi3, lpi4;
    double phimax = 0, phimin = 10;
    double facephimax = 0, facephimin = 10;
    int illegaltets = 0, negativetets = 0, badtets = 0;

    for (i = 1; i <= GetNE(); i++)
      {
	int badel = 0;

	Element & el = VolumeElement(i);

	if (el.GetType() != TET)
	  {
	    VolumeElement(i).flags.badel = 0;
	    continue;
	  }

	if (el.Volume(Points()) < 0)
	  {
	    badel = 1;
	    negativetets++;
	  }
      

	if (!LegalTet (el)) 
	  {
	    badel = 1;
	    illegaltets++;
	    (*testout) << "illegal tet: " << i << " ";
	    for (j = 1; j <= el.GetNP(); j++)
	      (*testout) << el.PNum(j) << " ";
	    (*testout) << endl;
	  }
      
	
	// angles between faces
	for (lpi1 = 1; lpi1 <= 3; lpi1++)
	  for (lpi2 = lpi1+1; lpi2 <= 4; lpi2++)
	    {
	      lpi3 = 1;
	      while (lpi3 == lpi1 || lpi3 == lpi2)
		lpi3++;
	      lpi4 = 10 - lpi1 - lpi2 - lpi3;

	      const Point3d & p1 = Point (el.PNum(lpi1));
	      const Point3d & p2 = Point (el.PNum(lpi2));
	      const Point3d & p3 = Point (el.PNum(lpi3));
	      const Point3d & p4 = Point (el.PNum(lpi4));

	      Vec3d n(p1, p2);
	      n /= n.Length();
	      Vec3d v1(p1, p3);
	      Vec3d v2(p1, p4);

	      v1 -= (n * v1) * n;
	      v2 -= (n * v2) * n;

	      double cosphi = (v1 * v2) / (v1.Length() * v2.Length());
	      double phi = acos (cosphi);
	      if (phi > phimax) phimax = phi;
	      if (phi < phimin) phimin = phi;

	      if ((180/M_PI) * phi > badellimit)
		badel = 1;
	    }


	// angles in faces
	for (j = 1; j <= 4; j++)
	  {
	    Element2d face;
	    el.GetFace (j, face);
	    for (lpi1 = 1; lpi1 <= 3; lpi1++)
	      {
		lpi2 = lpi1 % 3 + 1;
		lpi3 = lpi2 % 3 + 1;

		const Point3d & p1 = Point (el.PNum(lpi1));
		const Point3d & p2 = Point (el.PNum(lpi2));
		const Point3d & p3 = Point (el.PNum(lpi3));

		Vec3d v1(p1, p2);
		Vec3d v2(p1, p3);
		double cosphi = (v1 * v2) / (v1.Length() * v2.Length());
		double phi = acos (cosphi);
		if (phi > facephimax) facephimax = phi;
		if (phi < facephimin) facephimin = phi;

		if ((180/M_PI) * phi > badellimit)
		  badel = 1;

	      }
	  }

       
	VolumeElement(i).flags.badel = badel;
	if (badel) badtets++;
      }

    if (!GetNE())
      {
	phimin = phimax = facephimin = facephimax = 0;
      }

    if (!retvalues)
      {
	PrintMessage (1, "");
	PrintMessage (1, "between planes:  phimin = ", (180/M_PI) * phimin,
		      " phimax = ", (180/M_PI) *phimax);
	PrintMessage (1, "inside planes:   phimin = ", (180/M_PI) * facephimin,
		      " phimax = ", (180/M_PI) * facephimax);
	PrintMessage (1, "");      
      }
    else
      {
	retvalues[0] = (180/M_PI) * facephimin;
	retvalues[1] = (180/M_PI) * facephimax;
	retvalues[2] = (180/M_PI) * phimin;
	retvalues[3] = (180/M_PI) * phimax;
      }
    PrintMessage (3, "negative tets: ", negativetets);
    PrintMessage (3, "illegal tets:  ", illegaltets);
    PrintMessage (3, "bad tets:      ", badtets);
  }


  int Mesh :: MarkIllegalElements ()
  {
    int cnt = 0;
    int i;

    for (i = 1; i <= GetNE(); i++)
      {
	LegalTet (VolumeElement(i));

	/*
	  Element & el = VolumeElement(i);
	  int leg1 = LegalTet (el);
	  el.flags.illegal_valid = 0;
	  int leg2 = LegalTet (el);

	  if (leg1 != leg2) 
	  {
	  cerr << "legal differs!!" << endl;
	  (*testout) << "legal differs" << endl;
	  (*testout) << "elnr = " << i << ", el = " << el
	  << " leg1 = " << leg1 << ", leg2 = " << leg2 << endl;
	  }
      
	  //      el.flags.illegal = !LegalTet (el);
	  */
	cnt += VolumeElement(i).Illegal();
      }
    return cnt;
  }

// #ifdef NONE
//   void Mesh :: AddIdentification (int pi1, int pi2, int identnr)
//   {
//     INDEX_2 pair(pi1, pi2);
//     //  pair.Sort();
//     identifiedpoints->Set (pair, identnr);
//     if (identnr > maxidentnr)
//       maxidentnr = identnr;
//     timestamp = NextTimeStamp();
//   }

//   int Mesh :: GetIdentification (int pi1, int pi2) const
//   {
//     INDEX_2 pair(pi1, pi2);
//     if (identifiedpoints->Used (pair))
//       return identifiedpoints->Get(pair);
//     else
//       return 0;
//   }

//   int Mesh :: GetIdentificationSym (int pi1, int pi2) const
//   {
//     INDEX_2 pair(pi1, pi2);
//     if (identifiedpoints->Used (pair))
//       return identifiedpoints->Get(pair);

//     pair = INDEX_2 (pi2, pi1);
//     if (identifiedpoints->Used (pair))
//       return identifiedpoints->Get(pair);

//     return 0;
//   }


//   void Mesh :: GetIdentificationMap (int identnr, ARRAY<int> & identmap) const
//   {
//     int i, j;

//     identmap.SetSize (GetNP());
//     for (i = 1; i <= identmap.Size(); i++)
//       identmap.Elem(i) = 0;

//     for (i = 1; i <= identifiedpoints->GetNBags(); i++)
//       for (j = 1; j <= identifiedpoints->GetBagSize(i); j++)
// 	{
// 	  INDEX_2 i2;
// 	  int nr;
// 	  identifiedpoints->GetData (i, j, i2, nr);
	
// 	  if (nr == identnr)
// 	    {
// 	      identmap.Elem(i2.I1()) = i2.I2();
// 	    }
// 	}
//   }


//   void Mesh :: GetIdentificationPairs (int identnr, ARRAY<INDEX_2> & identpairs) const
//   {
//     int i, j;

//     identpairs.SetSize(0);

//     for (i = 1; i <= identifiedpoints->GetNBags(); i++)
//       for (j = 1; j <= identifiedpoints->GetBagSize(i); j++)
// 	{
// 	  INDEX_2 i2;
// 	  int nr;
// 	  identifiedpoints->GetData (i, j, i2, nr);
	
// 	  if (identnr == 0 || nr == identnr)
// 	    identpairs.Append (i2);
// 	}
//   }
// #endif

  

  void Mesh :: InitPointCurve(double red, double green, double blue) const
  {
    pointcurves_startpoint.Append(pointcurves.Size());
    pointcurves_red.Append(red);
    pointcurves_green.Append(green);
    pointcurves_blue.Append(blue);
  }
  void Mesh :: AddPointCurvePoint(const Point3d & pt) const
  {
    pointcurves.Append(pt);
  }
  int Mesh :: GetNumPointCurves(void) const
  {
    return pointcurves_startpoint.Size();
  }
  int Mesh :: GetNumPointsOfPointCurve(int curve) const
  {
    if(curve == pointcurves_startpoint.Size()-1)
      return (pointcurves.Size() - pointcurves_startpoint.Last());
    else
      return (pointcurves_startpoint[curve+1]-pointcurves_startpoint[curve]);
  }

  Point3d & Mesh :: GetPointCurvePoint(int curve, int n) const
  {
    return pointcurves[pointcurves_startpoint[curve]+n];
  }

  void Mesh :: GetPointCurveColor(int curve, double & red, double & green, double & blue) const
  {
    red = pointcurves_red[curve];
    green = pointcurves_green[curve];
    blue = pointcurves_blue[curve];
  }
  

  void Mesh :: ComputeNVertices ()
  {
    int i, j, nv;
    int ne = GetNE();
    int nse = GetNSE();

    numvertices = 0;
    for (i = 1; i <= ne; i++)
      {
	const Element & el = VolumeElement(i);
	nv = el.GetNV();
	for (j = 0; j < nv; j++)
	  if (el[j] > numvertices)
	    numvertices = el[j];
      }
    for (i = 1; i <= nse; i++)
      {
	const Element2d & el = SurfaceElement(i);
	nv = el.GetNV();
	for (j = 1; j <= nv; j++)
	  if (el.PNum(j) > numvertices)
	    numvertices = el.PNum(j);
      } 

    numvertices += 1- PointIndex::BASE;
  }

  int Mesh :: GetNV () const
  {
    if (numvertices < 0)
      return GetNP();
    else
      return numvertices;
  }

  void Mesh :: SetNP (int np)
  {
    points.SetSize(np);
    //  ptyps.SetSize(np);

    int mlold = mlbetweennodes.Size();
    mlbetweennodes.SetSize(np);
    if (np > mlold)
      for (int i = mlold+PointIndex::BASE; 
	   i < np+PointIndex::BASE; i++)
	{
	  mlbetweennodes[i].I1() = PointIndex::BASE-1;
	  mlbetweennodes[i].I2() = PointIndex::BASE-1;
	}

    GetIdentifications().SetMaxPointNr (np + PointIndex::BASE-1);
  }


  /*
    void Mesh :: BuildConnectedNodes ()
    {
    if (PureTetMesh())
    {
    connectedtonode.SetSize(0);
    return;
    }


    int i, j, k;
    int np = GetNP();
    int ne = GetNE();
    TABLE<int> conto(np);
    for (i = 1; i <= ne; i++)
    {
    const Element & el = VolumeElement(i);

    if (el.GetType() == PRISM)
    {
    for (j = 1; j <= 6; j++)
    {
    int n1 = el.PNum (j);
    int n2 = el.PNum ((j+2)%6+1);
    //	    if (n1 != n2)
    {
    int found = 0;
    for (k = 1; k <= conto.EntrySize(n1); k++)
    if (conto.Get(n1, k) == n2)
    {
    found = 1;
    break;
    }
    if (!found)
    conto.Add (n1, n2);
    }
    }
    }
    else if (el.GetType() == PYRAMID)
    {
    for (j = 1; j <= 4; j++)
    {
    int n1, n2;
    switch (j)
    {
    case 1: n1 = 1; n2 = 4; break;
    case 2: n1 = 4; n2 = 1; break;
    case 3: n1 = 2; n2 = 3; break;
    case 4: n1 = 3; n2 = 2; break;
    }

    int found = 0;
    for (k = 1; k <= conto.EntrySize(n1); k++)
    if (conto.Get(n1, k) == n2)
    {
    found = 1;
    break;
    }
    if (!found)
    conto.Add (n1, n2);
    }
    }
    }
  
    connectedtonode.SetSize(np);
    for (i = 1; i <= np; i++)
    connectedtonode.Elem(i) = 0;
  
    for (i = 1; i <= np; i++)
    if (connectedtonode.Elem(i) == 0)
    {
    connectedtonode.Elem(i) = i;
    ConnectToNodeRec (i, i, conto);
    }
  


    }

    void Mesh :: ConnectToNodeRec (int node, int tonode, 
    const TABLE<int> & conto)
    {
    int i, n2;
    //  (*testout) << "connect " << node << " to " << tonode << endl;
    for (i = 1; i <= conto.EntrySize(node); i++)
    {
    n2 = conto.Get(node, i);
    if (!connectedtonode.Get(n2))
    {
    connectedtonode.Elem(n2) = tonode;
    ConnectToNodeRec (n2, tonode, conto);
    }
    }
    }
  */


  bool Mesh :: PureTrigMesh (int faceindex) const
  {
    if (!faceindex)
      return !mparam.quad;

    int i;
    for (i = 1; i <= GetNSE(); i++)
      if (SurfaceElement(i).GetIndex() == faceindex &&
	  SurfaceElement(i).GetNP() != 3)
	return 0;
    return 1;
  }

  bool Mesh :: PureTetMesh () const
  {
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      if (VolumeElement(ei).GetNP() != 4)
	return 0;
    return 1;
  }

  void Mesh :: UpdateTopology()
  {
    topology->Update();
    clusters->Update();
  }


  void Mesh :: SetMaterial (int domnr, const char * mat)
  {
    if (domnr > materials.Size())
      {
	int olds = materials.Size();
	materials.SetSize (domnr);
	for (int i = olds; i < domnr; i++)
	  materials[i] = 0;
      }
    materials.Elem(domnr) = new char[strlen(mat)+1];
    strcpy (materials.Elem(domnr), mat);
  }

  const char * Mesh :: GetMaterial (int domnr) const
  {
    if (domnr <= materials.Size())
      return materials.Get(domnr);
    return 0;
  }

  void Mesh ::SetNBCNames ( int nbcn )
  {
    if ( bcnames.Size() )
      for ( int i = 0; i < bcnames.Size(); i++)
	if ( bcnames[i] ) delete bcnames[i];
    bcnames.SetSize(nbcn);
    bcnames = 0;
  }

  void Mesh ::SetBCName ( int bcnr, const string & abcname )
  {
    if ( bcnames[bcnr] ) delete bcnames[bcnr];
    if ( abcname != "default" )
      bcnames[bcnr] = new string ( abcname );
    else
      bcnames[bcnr] = 0;
  }

  string Mesh ::GetBCName ( int bcnr ) const
  {
    if ( !bcnames.Size() )
      return "default";
    if ( bcnames[bcnr] )
      return *bcnames[bcnr];
    else
      return "default";
  }

  void Mesh :: SetUserData(const char * id, ARRAY<int> & data)
  {
    if(userdata_int.Used(id))
      delete userdata_int.Get(id);

    ARRAY<int> * newdata = new ARRAY<int>(data);

    userdata_int.Set(id,newdata);      
  }
  bool Mesh :: GetUserData(const char * id, ARRAY<int> & data, int shift) const
  {
    if(userdata_int.Used(id))
      {
	if(data.Size() < (*userdata_int.Get(id)).Size()+shift)
	  data.SetSize((*userdata_int.Get(id)).Size()+shift);
	for(int i=0; i<(*userdata_int.Get(id)).Size(); i++)
	  data[i+shift] = (*userdata_int.Get(id))[i];
	return true;
      }
    else
      {
	data.SetSize(0);
	return false;
      }
  }
  void Mesh :: SetUserData(const char * id, ARRAY<double> & data)
  {
    if(userdata_double.Used(id))
      delete userdata_double.Get(id);

    ARRAY<double> * newdata = new ARRAY<double>(data);

    userdata_double.Set(id,newdata);      
  }
  bool Mesh :: GetUserData(const char * id, ARRAY<double> & data, int shift) const
  {
    if(userdata_double.Used(id))
      {
	if(data.Size() < (*userdata_double.Get(id)).Size()+shift)
	  data.SetSize((*userdata_double.Get(id)).Size()+shift);
	for(int i=0; i<(*userdata_double.Get(id)).Size(); i++)
	  data[i+shift] = (*userdata_double.Get(id))[i];
	return true;
      }
    else
      {
	data.SetSize(0);
	return false;
      }
  }
  


  void Mesh :: PrintMemInfo (ostream & ost) const
  {
    ost << "Mesh Mem:" << endl;

    ost << GetNP() << " Points, of size " 
	<< sizeof (Point3d) << " + " << sizeof(POINTTYPE) << " = "
	<< GetNP() * (sizeof (Point3d) + sizeof(POINTTYPE)) << endl;

    ost << GetNSE() << " Surface elements, of size " 
	<< sizeof (Element2d) << " = " 
	<< GetNSE() * sizeof(Element2d) << endl;

    ost << GetNE() << " Volume elements, of size " 
	<< sizeof (Element) << " = " 
	<< GetNE() * sizeof(Element) << endl;

    ost << "surfs on node:";
    surfacesonnode.PrintMemInfo (cout);
  
    ost << "boundaryedges: ";
    if (boundaryedges)
      boundaryedges->PrintMemInfo (cout);

    ost << "surfelementht: ";
    if (surfelementht)
      surfelementht->PrintMemInfo (cout);
  }
}
