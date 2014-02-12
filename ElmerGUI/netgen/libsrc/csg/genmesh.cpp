#include <mystdlib.h>


#include <myadt.hpp>

#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>


namespace netgen
{
  ARRAY<SpecialPoint> specpoints;
  static ARRAY<MeshPoint> spoints;

#define TCL_OK 0
#define TCL_ERROR 1



  static void FindPoints (CSGeometry & geom, Mesh & mesh)
  {
    PrintMessage (1, "Start Findpoints");

    const char * savetask = multithread.task;
    multithread.task = "Find points";

    for (int i = 0; i < geom.GetNUserPoints(); i++)
      {
	mesh.AddPoint(geom.GetUserPoint (i));
	mesh.Points().Last().Singularity (geom.GetUserPointRefFactor(i));
	mesh.AddLockedPoint (PointIndex (i+1));
      }

    SpecialPointCalculation spc;

    spc.SetIdEps(geom.GetIdEps());

    if (spoints.Size() == 0)
      spc.CalcSpecialPoints (geom, spoints);
    
    PrintMessage (2, "Analyze spec points");
    spc.AnalyzeSpecialPoints (geom, spoints, specpoints);
  
    PrintMessage (5, "done");

    (*testout) << specpoints.Size() << " special points:" << endl;
    for (int i = 0; i < specpoints.Size(); i++)
      specpoints[i].Print (*testout);

    /*
      for (int i = 1; i <= geom.identifications.Size(); i++)
      geom.identifications.Elem(i)->IdentifySpecialPoints (specpoints);
    */
    multithread.task = savetask;
  }






  static void FindEdges (CSGeometry & geom, Mesh & mesh, const bool setmeshsize = false)
  {
    EdgeCalculation ec (geom, specpoints);
    ec.SetIdEps(geom.GetIdEps());
    ec.Calc (mparam.maxh, mesh);

    for (int i = 0; i < geom.singedges.Size(); i++)
      {
	geom.singedges[i]->FindPointsOnEdge (mesh);
	if(setmeshsize)
	  geom.singedges[i]->SetMeshSize(mesh,10.*geom.BoundingBox().Diam());
      }
    for (int i = 0; i < geom.singpoints.Size(); i++)
      geom.singpoints[i]->FindPoints (mesh);

    for (int i = 1; i <= mesh.GetNSeg(); i++)
      {
	//(*testout) << "segment " << mesh.LineSegment(i) << endl;
	int ok = 0;
	for (int k = 1; k <= mesh.GetNFD(); k++)
	  if (mesh.GetFaceDescriptor(k).SegmentFits (mesh.LineSegment(i)))
	    {
	      ok = k;
	      //(*testout) << "fits to " << k << endl;
	    }

	if (!ok)
	  {
	    ok = mesh.AddFaceDescriptor (FaceDescriptor (mesh.LineSegment(i)));
	    //(*testout) << "did not find, now " << ok << endl;
	  }

	//(*testout) << "change from " << mesh.LineSegment(i).si;
	mesh.LineSegment(i).si = ok;
	//(*testout) << " to " << mesh.LineSegment(i).si << endl;
      }

    if (geom.identifications.Size())
      {
	PrintMessage (3, "Find Identifications");
	for (int i = 0; i < geom.identifications.Size(); i++)
	  {
	    geom.identifications[i]->IdentifyPoints (mesh);
	    //(*testout) << "identification " << i << " is " 
	    //	       << *geom.identifications[i] << endl;
	    
	  }
	for (int i = 0; i < geom.identifications.Size(); i++)
	  geom.identifications[i]->IdentifyFaces (mesh);
      }


    // find intersecting segments
    PrintMessage (3, "Check intersecting edges");
    
    Point3d pmin, pmax;
    mesh.GetBox (pmin, pmax);
    Box3dTree segtree (pmin, pmax);
    
    for (SegmentIndex si = 0; si < mesh.GetNSeg(); si++)
      {
	if (mesh[si].seginfo)
	  {
	    Box<3> hbox;
	    hbox.Set (mesh[mesh[si].p1]);
	    hbox.Add (mesh[mesh[si].p2]);
	    segtree.Insert (hbox.PMin(), hbox.PMax(), si);
	  }
      }

    ARRAY<int> loc;
    if (!ec.point_on_edge_problem)
      for (SegmentIndex si = 0; si < mesh.GetNSeg(); si++)
	{
	  if (!mesh[si].seginfo) continue;

	  Box<3> hbox;
	  hbox.Set (mesh[mesh[si].p1]);
	  hbox.Add (mesh[mesh[si].p2]);
	  hbox.Increase (1e-6);
	  segtree.GetIntersecting (hbox.PMin(), hbox.PMax(), loc);
	  	  
	  // for (SegmentIndex sj = 0; sj < si; sj++)
	  for (int j = 0; j < loc.Size(); j++)
	    {
	      SegmentIndex sj = loc[j];
	      if (sj >= si) continue;
	      if (!mesh[si].seginfo || !mesh[sj].seginfo) continue;
	      if (mesh[mesh[si].p1].GetLayer() != mesh[mesh[sj].p2].GetLayer()) continue;
	      
	      Point<3> pi1 = mesh[mesh[si].p1];
	      Point<3> pi2 = mesh[mesh[si].p2];
	      Point<3> pj1 = mesh[mesh[sj].p1];
	      Point<3> pj2 = mesh[mesh[sj].p2];
	      Vec<3> vi = pi2 - pi1;
	      Vec<3> vj = pj2 - pj1;
	      
	      if (sqr (vi * vj) > (1.-1e-6) * Abs2 (vi) * Abs2 (vj)) continue;
	      
	      // pi1 + vi t = pj1 + vj s
	      Mat<3,2> mat;
	      Vec<3> rhs;
	      Vec<2> sol;
	      
	      for (int jj = 0; jj < 3; jj++)
		{ 
		  mat(jj,0) = vi(jj); 
		  mat(jj,1) = -vj(jj); 
		  rhs(jj) = pj1(jj)-pi1(jj); 
		}
	      
	      mat.Solve (rhs, sol);

	      //(*testout) << "mat " << mat << endl << "rhs " << rhs << endl << "sol " << sol << endl;
	      
	      if (sol(0) > 1e-6 && sol(0) < 1-1e-6 &&
		  sol(1) > 1e-6 && sol(1) < 1-1e-6 &&
		  Abs (rhs - mat*sol) < 1e-6)
		{
		  Point<3> ip = pi1 + sol(0) * vi;
		  
		  //(*testout) << "ip " << ip << endl;

		  Point<3> pip = ip;
		  ProjectToEdge (geom.GetSurface (mesh[si].surfnr1),
				 geom.GetSurface (mesh[si].surfnr2), pip);
		  
		  //(*testout) << "Dist (ip, pip_si) " << Dist (ip, pip) << endl;
		  if (Dist (ip, pip) > 1e-6*geom.MaxSize()) continue;
		  pip = ip;
		  ProjectToEdge (geom.GetSurface (mesh[sj].surfnr1),
				 geom.GetSurface (mesh[sj].surfnr2), pip);

		  //(*testout) << "Dist (ip, pip_sj) " << Dist (ip, pip) << endl;
		  if (Dist (ip, pip) > 1e-6*geom.MaxSize()) continue;
		  
		  
		  
		  cout << "Intersection at " << ip << endl;
		  
		  geom.AddUserPoint (ip);
		  spoints.Append (MeshPoint (ip, mesh[mesh[si].p1].GetLayer()));
		  mesh.AddPoint (ip);
		  
		  (*testout) << "found intersection at " << ip << endl;
		  (*testout) << "sol = " << sol << endl;
		  (*testout) << "res = " << (rhs - mat*sol) << endl;
		  (*testout) << "segs = " << pi1 << " - " << pi2 << endl;
		  (*testout) << "and = " << pj1 << " - " << pj2 << endl << endl;
		}
	    }
	}  
  }






  static void MeshSurface (CSGeometry & geom, Mesh & mesh)
  {
    const char * savetask = multithread.task;
    multithread.task = "Surface meshing";
  
    ARRAY<Segment> segments;
    int noldp = mesh.GetNP();

    double starttime = GetTime();

    // find master faces from identified
    ARRAY<int> masterface(mesh.GetNFD());
    for (int i = 1; i <= mesh.GetNFD(); i++)
      masterface.Elem(i) = i;
  
    ARRAY<INDEX_2> fpairs;
    bool changed;
    do
      {
	changed = 0;
	for (int i = 0; i < geom.identifications.Size(); i++)
	  {
	    geom.identifications[i]->GetIdentifiedFaces (fpairs);

	    for (int j = 0; j < fpairs.Size(); j++)
	      {
		if (masterface.Get(fpairs[j].I1()) <
		    masterface.Get(fpairs[j].I2()))
		  {
		    changed = 1;
		    masterface.Elem(fpairs[j].I2()) =
		      masterface.Elem(fpairs[j].I1());
		  }
		if (masterface.Get(fpairs[j].I2()) <
		    masterface.Get(fpairs[j].I1()))
		  {
		    changed = 1;
		    masterface.Elem(fpairs[j].I1()) =
		      masterface.Elem(fpairs[j].I2());
		  }
	      }
	  }
      }
    while (changed);


    int bccnt=0;
    for (int k = 0; k < geom.GetNSurf(); k++)
      bccnt = max2 (bccnt, geom.GetSurface(k)->GetBCProperty());

    for (int k = 1; k <= mesh.GetNFD(); k++)
      {
	bool increased = false;

	FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
	const Surface * surf = geom.GetSurface(fd.SurfNr());

	if (fd.TLOSurface() && 
	    geom.GetTopLevelObject(fd.TLOSurface()-1) -> GetBCProp() > 0)
	  fd.SetBCProperty (geom.GetTopLevelObject(fd.TLOSurface()-1) -> GetBCProp());
	else if (surf -> GetBCProperty() != -1)
	  fd.SetBCProperty (surf->GetBCProperty());
	else
	  {
	    bccnt++;
	    fd.SetBCProperty (bccnt);
	    increased = true;
	  }      

	for (int l = 0; l < geom.bcmodifications.Size(); l++)
	  {
	    if (geom.GetSurfaceClassRepresentant (fd.SurfNr()) == 
		geom.GetSurfaceClassRepresentant (geom.bcmodifications[l].si) &&
		(fd.DomainIn() == geom.bcmodifications[l].tlonr+1 ||
		 fd.DomainOut() == geom.bcmodifications[l].tlonr+1))
	      {
		if(geom.bcmodifications[l].bcname == NULL)
		  fd.SetBCProperty (geom.bcmodifications[l].bcnr);
		else
		  {
		    if(!increased)
		      {
			bccnt++;
			fd.SetBCProperty (bccnt);
			increased = true;
		      }
		  }
	      }
	  }
      }

    mesh.SetNBCNames( bccnt );

    for (int k = 1; k <= mesh.GetNFD(); k++)
      {
	FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
	const Surface * surf = geom.GetSurface(fd.SurfNr());
	if (fd.TLOSurface() )
	  {
	    int bcp = fd.BCProperty();
	    string nextbcname = geom.GetTopLevelObject(fd.TLOSurface()-1) -> GetBCName();
	    if ( nextbcname != "default" )
	      mesh.SetBCName ( bcp - 1 , nextbcname );
	  }
	else // if (surf -> GetBCProperty() != -1)
	  {
	    int bcp = fd.BCProperty();
	    string nextbcname = surf->GetBCName();
	    if ( nextbcname != "default" )
	      mesh.SetBCName ( bcp - 1, nextbcname );
	  }
      }
    
    for (int k = 1; k <= mesh.GetNFD(); k++)
      {
	FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
	fd.SetBCName ( mesh.GetBCNamePtr ( fd.BCProperty() - 1 ) );
      }
    

    //!!
    
    for (int k = 1; k <= mesh.GetNFD(); k++)
      {
	FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
	//const Surface * surf = geom.GetSurface(fd.SurfNr());

	for (int l = 0; l < geom.bcmodifications.Size(); l++)
	  {
	    if (geom.GetSurfaceClassRepresentant (fd.SurfNr()) == 
		geom.GetSurfaceClassRepresentant (geom.bcmodifications[l].si) &&
		(fd.DomainIn() == geom.bcmodifications[l].tlonr+1 ||
		 fd.DomainOut() == geom.bcmodifications[l].tlonr+1) &&
		geom.bcmodifications[l].bcname != NULL
		)
	      {
		int bcp = fd.BCProperty();
		mesh.SetBCName ( bcp - 1, *(geom.bcmodifications[l].bcname) );
		fd.SetBCName ( mesh.GetBCNamePtr ( bcp - 1) );
	      }
	  }
      }

    for(int k = 0; k<geom.bcmodifications.Size(); k++)
      {
	delete geom.bcmodifications[k].bcname;
	geom.bcmodifications[k].bcname = NULL;
      }

    //!!


    for (int j = 0; j < geom.singfaces.Size(); j++)
      {
	ARRAY<int> surfs;
	geom.GetIndependentSurfaceIndices (geom.singfaces[j]->GetSolid(),
					   geom.BoundingBox(), surfs);
	for (int k = 1; k <= mesh.GetNFD(); k++)
	  {
	    FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
	    for (int l = 0; l < surfs.Size(); l++)
	      if (surfs[l] == fd.SurfNr())
		{
		  if (geom.singfaces[j]->GetDomainNr() == fd.DomainIn())
		    fd.domin_singular = 1;
		  if (geom.singfaces[j]->GetDomainNr() == fd.DomainOut())
		    fd.domout_singular = 1;
		}
	  }
      }
    

    // assemble edge hash-table
    mesh.CalcSurfacesOfNode();

    for (int k = 1; k <= mesh.GetNFD(); k++)
      {
	multithread.percent = 100.0 * k / (mesh.GetNFD()+1e-10);

	if (masterface.Get(k) != k)
	  continue;

	FaceDescriptor & fd = mesh.GetFaceDescriptor(k);

	(*testout) << "Surface " << k << endl;
	(*testout) << "Face Descriptor: " << fd << endl;
	PrintMessage (1, "Surface ", k, " / ", mesh.GetNFD());

	int oldnf = mesh.GetNSE();
      
	const Surface * surf =
	  geom.GetSurface((mesh.GetFaceDescriptor(k).SurfNr()));


	Meshing2Surfaces meshing(*surf, geom.BoundingBox());
	meshing.SetStartTime (starttime);


	for (PointIndex pi = PointIndex::BASE; pi < noldp+PointIndex::BASE; pi++)
	  { 
	    //if(surf->PointOnSurface(mesh[pi]))
	    meshing.AddPoint (mesh[pi], pi, NULL,
			      (surf->PointOnSurface(mesh[pi])!=0));
	  }

	segments.SetSize (0);

	for (SegmentIndex si = 0; si < mesh.GetNSeg(); si++)
	  if (mesh[si].si == k)
	    {
	      segments.Append (mesh[si]);
	      (*testout) << "appending segment " << mesh[si] << endl;
	      //<< " from " << mesh[mesh[si][0]]
	      //	 << " to " <<mesh[mesh[si][1]]<< endl;
	    }

	(*testout) << "num-segments " << segments.Size() << endl;

	for (int i = 1; i <= geom.identifications.Size(); i++)
	  {
	    geom.identifications.Get(i)->
	      BuildSurfaceElements(segments, mesh, surf);
	  }

	for (int si = 0; si < segments.Size(); si++)
	  {
	    PointGeomInfo gi;
	    gi.trignum = k;
	    meshing.AddBoundaryElement (segments[si].p1 + 1 - PointIndex::BASE, 
					segments[si].p2 + 1 - PointIndex::BASE, 
					gi, gi);
	  }

	double maxh = mparam.maxh;
	if (fd.DomainIn() != 0)
	  {
	    const Solid * s1 = 
	      geom.GetTopLevelObject(fd.DomainIn()-1) -> GetSolid();
	    if (s1->GetMaxH() < maxh)
	      maxh = s1->GetMaxH();
	    maxh = min2(maxh, geom.GetTopLevelObject(fd.DomainIn()-1)->GetMaxH());
	  }
	if (fd.DomainOut() != 0)
	  {
	    const Solid * s1 = 
	      geom.GetTopLevelObject(fd.DomainOut()-1) -> GetSolid();
	    if (s1->GetMaxH() < maxh)
	      maxh = s1->GetMaxH();
	    maxh = min2(maxh, geom.GetTopLevelObject(fd.DomainOut()-1)->GetMaxH());
	  }
	if (fd.TLOSurface() != 0)
	  {
	    double hi = geom.GetTopLevelObject(fd.TLOSurface()-1) -> GetMaxH();
	    if (hi < maxh) maxh = hi;
	  }

	(*testout) << "domin = " << fd.DomainIn() << ", domout = " << fd.DomainOut()
		   << ", tlo-surf = " << fd.TLOSurface()
		   << " mpram.maxh = " << mparam.maxh << ", maxh = " << maxh << endl;

	mparam.checkoverlap = 0;

	MESHING2_RESULT res =
	  meshing.GenerateMesh (mesh, maxh, k);

	if (res != MESHING2_OK)
	  {
	    PrintError ("Problem in Surface mesh generation");
	    throw NgException ("Problem in Surface mesh generation");
	  }

	if (multithread.terminate) return;
      
	for (int i = oldnf+1; i <= mesh.GetNSE(); i++)
	  mesh.SurfaceElement(i).SetIndex (k);


	//      mesh.CalcSurfacesOfNode();

	if (segments.Size())   
	  { 
	    // surface was meshed, not copied
	    PrintMessage (2, "Optimize Surface");
	    for (int i = 1; i <= mparam.optsteps2d; i++)
	      {
		if (multithread.terminate) return;

		{
		  MeshOptimize2dSurfaces meshopt(geom);
		  meshopt.SetFaceIndex (k);
		  meshopt.SetImproveEdges (0);
		  meshopt.SetMetricWeight (mparam.elsizeweight);
		  meshopt.SetWriteStatus (0);

		  meshopt.EdgeSwapping (mesh, (i > mparam.optsteps2d/2));
		}
		
		if (multithread.terminate) return;
		{
		  //		mesh.CalcSurfacesOfNode();
		
		  MeshOptimize2dSurfaces meshopt(geom);
		  meshopt.SetFaceIndex (k);
		  meshopt.SetImproveEdges (0);
		  meshopt.SetMetricWeight (mparam.elsizeweight);
		  meshopt.SetWriteStatus (0);

		  meshopt.ImproveMesh (mesh);
		}
		
		{
		  MeshOptimize2dSurfaces meshopt(geom);
		  meshopt.SetFaceIndex (k);
		  meshopt.SetImproveEdges (0);
		  meshopt.SetMetricWeight (mparam.elsizeweight);
		  meshopt.SetWriteStatus (0);

		  meshopt.CombineImprove (mesh);
		  //		mesh.CalcSurfacesOfNode();
		}
		
		if (multithread.terminate) return;
		{
		  MeshOptimize2dSurfaces meshopt(geom);
		  meshopt.SetFaceIndex (k);
		  meshopt.SetImproveEdges (0);
		  meshopt.SetMetricWeight (mparam.elsizeweight);
		  meshopt.SetWriteStatus (0);

		  meshopt.ImproveMesh (mesh);
		}
	      }
	  }


	PrintMessage (3, (mesh.GetNSE() - oldnf), " elements, ", mesh.GetNP(), " points");

#ifdef OPENGL
	extern void Render();
	Render();
#endif
      }
    
    mesh.Compress();

    do
      {
	changed = 0;
	for (int k = 1; k <= mesh.GetNFD(); k++)
	  {
	    multithread.percent = 100.0 * k / (mesh.GetNFD()+1e-10);
	  
	    if (masterface.Get(k) == k)
	      continue;

	    FaceDescriptor & fd = mesh.GetFaceDescriptor(k);

	    (*testout) << "Surface " << k << endl;
	    (*testout) << "Face Descriptor: " << fd << endl;
	    PrintMessage (2, "Surface ", k);

	    int oldnf = mesh.GetNSE();
      
	    const Surface * surf =
	      geom.GetSurface((mesh.GetFaceDescriptor(k).SurfNr()));

	    /*
	      if (surf -> GetBCProperty() != -1)
	      fd.SetBCProperty (surf->GetBCProperty());
	      else
	      {
	      bccnt++;
	      fd.SetBCProperty (bccnt);
	      }
	    */
  
	    segments.SetSize (0);
	    for (int i = 1; i <= mesh.GetNSeg(); i++)
	      {
		Segment * seg = &mesh.LineSegment(i);
		if (seg->si == k)
		  segments.Append (*seg);
	      }

	    for (int i = 1; i <= geom.identifications.Size(); i++)
	      {
		geom.identifications.Elem(i)->GetIdentifiedFaces (fpairs);
		int found = 0;
		for (int j = 1; j <= fpairs.Size(); j++)
		  if (fpairs.Get(j).I1() == k || fpairs.Get(j).I2() == k)
		    found = 1;

		if (!found)
		  continue;

		geom.identifications.Get(i)->
		  BuildSurfaceElements(segments, mesh, surf);
		if (!segments.Size())
		  break;
	      }

	  
	    if (multithread.terminate) return;

	    for (int i = oldnf+1; i <= mesh.GetNSE(); i++)
	      mesh.SurfaceElement(i).SetIndex (k);


	    if (!segments.Size())
	      {
		masterface.Elem(k) = k;
		changed = 1; 
	      }

	    PrintMessage (3, (mesh.GetNSE() - oldnf), " elements, ", mesh.GetNP(), " points");
	  }
      
#ifdef OPENGL
	extern void Render();
	Render();
#endif
      }
    while (changed);

    
    mesh.SplitSeparatedFaces();
    mesh.CalcSurfacesOfNode();

    multithread.task = savetask;
  }







  int GenerateMesh (CSGeometry & geom,
		    Mesh *& mesh,
		    int perfstepsstart, int perfstepsend,
		    const char * optstr)
  {

    if (mesh && mesh->GetNSE() &&
	!geom.GetNSolids())
      {
	if (perfstepsstart < MESHCONST_MESHVOLUME)
	  perfstepsstart = MESHCONST_MESHVOLUME;
      }



    if (perfstepsstart <= MESHCONST_ANALYSE)
      {
	delete mesh;
	mesh = new Mesh();

	mesh->SetGlobalH (mparam.maxh);
	mesh->SetMinimalH (mparam.minh);

	ARRAY<double> maxhdom(geom.GetNTopLevelObjects());
	for (int i = 0; i < maxhdom.Size(); i++)
	  maxhdom[i] = geom.GetTopLevelObject(i)->GetMaxH();

	mesh->SetMaxHDomain (maxhdom);

	if (mparam.uselocalh)
	  {
	    double maxsize = geom.MaxSize(); 
	    mesh->SetLocalH (Point<3>(-maxsize, -maxsize, -maxsize),
			     Point<3>(maxsize, maxsize, maxsize),
			     mparam.grading);

	    mesh -> LoadLocalMeshSize (mparam.meshsizefilename);
	  }

	spoints.SetSize(0);
	FindPoints (geom, *mesh);
      
	PrintMessage (5, "find points done");

#ifdef LOG_STREAM
	(*logout) << "Special points found" << endl
		  << "time = " << GetTime() << " sec" << endl
		  << "points: " << mesh->GetNP() << endl << endl;
#endif
      }


    if (multithread.terminate || perfstepsend <= MESHCONST_ANALYSE) 
      return TCL_OK;


    if (perfstepsstart <= MESHCONST_MESHEDGES)
      {
	FindEdges (geom, *mesh, true);
	if (multithread.terminate) return TCL_OK;
#ifdef LOG_STREAM      
	(*logout) << "Edges meshed" << endl
		  << "time = " << GetTime() << " sec" << endl
		  << "points: " << mesh->GetNP() << endl;
#endif
      
      
	if (multithread.terminate)
	  return TCL_OK;
  
	if (mparam.uselocalh)
	  {
	    mesh->CalcLocalH();
	    mesh->DeleteMesh();
	    
	    FindPoints (geom, *mesh);
	    if (multithread.terminate) return TCL_OK;
	    FindEdges (geom, *mesh, true);
	    if (multithread.terminate) return TCL_OK;
	    
	    mesh->DeleteMesh();
	  
	    FindPoints (geom, *mesh);
	    if (multithread.terminate) return TCL_OK;
	    FindEdges (geom, *mesh);
	    if (multithread.terminate) return TCL_OK;
	  }
      }
  
    if (multithread.terminate || perfstepsend <= MESHCONST_MESHEDGES)
      return TCL_OK;


    if (perfstepsstart <= MESHCONST_MESHSURFACE)
      {
	MeshSurface (geom, *mesh);  
	if (multithread.terminate) return TCL_OK;
      
#ifdef LOG_STREAM
	(*logout) << "Surfaces meshed" << endl
		  << "time = " << GetTime() << " sec" << endl
		  << "points: " << mesh->GetNP() << endl;
#endif      
      
	if (mparam.uselocalh && 0)
	  {
	    mesh->CalcLocalH();      
	    mesh->DeleteMesh();

	    FindPoints (geom, *mesh);
	    if (multithread.terminate) return TCL_OK;
	    FindEdges (geom, *mesh);
	    if (multithread.terminate) return TCL_OK;

	    MeshSurface (geom, *mesh);  
	    if (multithread.terminate) return TCL_OK;
	  }

#ifdef LOG_STREAM      
	(*logout) << "Surfaces remeshed" << endl
		  << "time = " << GetTime() << " sec" << endl
		  << "points: " << mesh->GetNP() << endl;
#endif      
      
#ifdef STAT_STREAM
	(*statout) << mesh->GetNSeg() << " & "
		   << mesh->GetNSE() << " & - &" 
		   << GetTime() << " & " << endl;
#endif  

	MeshQuality2d (*mesh);
	mesh->CalcSurfacesOfNode();
      }
  
    if (multithread.terminate || perfstepsend <= MESHCONST_OPTSURFACE)
      return TCL_OK;


    if (perfstepsstart <= MESHCONST_MESHVOLUME)
      {
	multithread.task = "Volume meshing";

	MESHING3_RESULT res =
	  MeshVolume (mparam, *mesh);

	if (res != MESHING3_OK) return TCL_ERROR;
      
	if (multithread.terminate) return TCL_OK;
      
	RemoveIllegalElements (*mesh);
	if (multithread.terminate) return TCL_OK;

	MeshQuality3d (*mesh);
      
	for (int i = 0; i < geom.GetNTopLevelObjects(); i++)
	  mesh->SetMaterial (i+1, geom.GetTopLevelObject(i)->GetMaterial().c_str());
      

#ifdef STAT_STREAM
	(*statout) << GetTime() << " & ";
#endif      
      
#ifdef LOG_STREAM
	(*logout) << "Volume meshed" << endl
		  << "time = " << GetTime() << " sec" << endl
		  << "points: " << mesh->GetNP() << endl;
#endif
      }

    if (multithread.terminate || perfstepsend <= MESHCONST_MESHVOLUME)
      return TCL_OK;


    if (perfstepsstart <= MESHCONST_OPTVOLUME)
      {
	multithread.task = "Volume optimization";
      
	OptimizeVolume (mparam, *mesh);
	if (multithread.terminate) return TCL_OK;
      
#ifdef STAT_STREAM
	(*statout) << GetTime() << " & "
		   << mesh->GetNE() << " & "
		   << mesh->GetNP() << " " << '\\' << '\\' << " \\" << "hline" << endl;
#endif      

#ifdef LOG_STREAM      
	(*logout) << "Volume optimized" << endl
		  << "time = " << GetTime() << " sec" << endl
		  << "points: " << mesh->GetNP() << endl;
#endif
      }

    return TCL_OK;
  }
}
