#include <mystdlib.h>

#include "meshing.hpp"

#ifdef PARALLEL
#include <parallel.hpp>
#endif

namespace netgen
{

MeshTopology ::  MeshTopology (const Mesh & amesh)
  : mesh(amesh)
{
  buildedges = 1;
  buildfaces = 1;
  vert2element = 0;
  vert2surfelement = 0;
  vert2segment = 0;
  timestamp = -1;

  edge2vert.SetName ("edge2vert");
  face2vert.SetName ("face2vert");
  edges.SetName ("el2edge");
  faces.SetName ("el2face");
  surfedges.SetName ("surfel2edge");
  segedges.SetName ("segment2edge");
  surffaces.SetName ("surfel2face");
  surf2volelement.SetName ("surfel2el");
  face2surfel.SetName ("face2surfel");
}

MeshTopology :: ~MeshTopology ()
{
  delete vert2element;
  delete vert2surfelement;
  delete vert2segment;
}

void MeshTopology :: Update()
{
  static int timer = NgProfiler::CreateTimer ("topology");
  NgProfiler::RegionTimer reg (timer);

#ifdef PARALLEL
  ParallelMeshTopology & paralleltop = mesh.GetParallelTopology();

  bool isparallel = 0;
#endif

  
  if (timestamp > mesh.GetTimeStamp()) return;
  
  int ne = mesh.GetNE();
  int nse = mesh.GetNSE();
  int nseg = mesh.GetNSeg();
  int np = mesh.GetNP();
  int nv = mesh.GetNV(); 
  int nfa = 0;
  int ned = edge2vert.Size();

  PrintMessage (3, "Update mesh topology");

   (*testout) << " UPDATE MESH TOPOLOGY " << endl; 
   (*testout) << "ne   = " << ne << endl;
   (*testout) << "nse  = " << nse << endl;
   (*testout) << "nseg = " << nseg << endl;
   (*testout) << "np   = " << np << endl;
   (*testout) << "nv   = " << nv << endl;
  
  delete vert2element;
  delete vert2surfelement;
  delete vert2segment;
  
  ARRAY<int,PointIndex::BASE> cnt(nv);
  ARRAY<int> vnums;

  /*
    generate:
    vertex to element 
    vertex to surface element
    vertex to segment 
   */


  cnt = 0;
  for (ElementIndex ei = 0; ei < ne; ei++)
    {
      const Element & el = mesh[ei];
      int nelv = el.GetNV();
      for (int j = 0; j < nelv; j++)
	cnt[el[j]]++;
    }

  vert2element = new TABLE<int,PointIndex::BASE> (cnt);
  for (ElementIndex ei = 0; ei < ne; ei++)
    {
      const Element & el = mesh[ei];
      int nelv = el.GetNV();
      for (int j = 0; j < nelv; j++)
	vert2element->AddSave (el[j], ei+1);
    }

  cnt = 0;
  for (SurfaceElementIndex sei = 0; sei < nse; sei++)
    {
      const Element2d & el = mesh[sei];
      int nelv = el.GetNV();
      for (int j = 0; j < nelv; j++)
	cnt[el[j]]++;
    }

  vert2surfelement = new TABLE<int,PointIndex::BASE> (cnt);
  for (SurfaceElementIndex sei = 0; sei < nse; sei++)
    {
      const Element2d & el = mesh[sei];
      int nelv = el.GetNV();
      for (int j = 0; j < nelv; j++)
	vert2surfelement->AddSave (el[j], sei+1);
    }

  cnt = 0;
  for (int i = 1; i <= nseg; i++)
    {
      const Segment & seg = mesh.LineSegment(i);
      cnt[seg.p1]++;
      cnt[seg.p2]++;
    }
 
  vert2segment = new TABLE<int,PointIndex::BASE> (cnt);
  for (int i = 1; i <= nseg; i++)
    {
      const Segment & seg = mesh.LineSegment(i);
      vert2segment->AddSave (seg.p1, i);
      vert2segment->AddSave (seg.p2, i);
    }

  if (buildedges)
    {
      static int timer1 = NgProfiler::CreateTimer ("topology::buildedges");
      NgProfiler::RegionTimer reg1 (timer1);

      PrintMessage (5, "Update edges ");
      
      edges.SetSize(ne);
      surfedges.SetSize(nse); 
      segedges.SetSize(nseg);

      for (int i = 0; i < ne; i++)
	for (int j = 0; j < 12; j++)
	  edges[i][j] = 0;
      for (int i = 0; i < nse; i++)
	for (int j = 0; j < 4; j++)
	  surfedges[i][j] = 0;

      // keep existing edges
      cnt = 0;
      for (int i = 0; i < edge2vert.Size(); i++)
	cnt[edge2vert[i][0]]++;
      TABLE<int,PointIndex::BASE> vert2edge (cnt);
      for (int i = 0; i < edge2vert.Size(); i++)
	vert2edge.AddSave (edge2vert[i][0], i+1);

      // ensure all coarse grid and intermediate level edges
      cnt = 0;
      for (int i = mesh.mlbetweennodes.Begin(); i < mesh.mlbetweennodes.End(); i++)
	{
	  int pa[2];
	  pa[0] = mesh.mlbetweennodes[i].I1();
	  pa[1] = mesh.mlbetweennodes[i].I2();
	  if (pa[0] > pa[1]) Swap (pa[0], pa[1]);
	  if (pa[0] > 0)
	    cnt.Elem(pa[0])++;
	}
      TABLE<int,PointIndex::BASE> vert2vertcoarse (cnt);
      for (int i = mesh.mlbetweennodes.Begin(); i < mesh.mlbetweennodes.End(); i++)
	{
	  int pa[2];
	  pa[0] = mesh.mlbetweennodes[i].I1();
	  pa[1] = mesh.mlbetweennodes[i].I2();
	  if (pa[0] > pa[1]) swap (pa[0], pa[1]);
	  if (pa[0] > 0)
	    vert2vertcoarse.AddSave1 (pa[0], pa[1]);
	}


      ARRAY<int,PointIndex::BASE> edgenr(nv);
      ARRAY<int,PointIndex::BASE> edgeflag(nv);
      edgeflag = 0;

      ned = edge2vert.Size();
      ARRAY<INDEX_3> missing;

      for (int i = 1; i <= nv; i++)
	{
	  for (int j = 1; j <= vert2edge.EntrySize(i); j++)
	    {
	      int ednr = vert2edge.Get(i,j);
	      int i2 = edge2vert.Get(ednr)[1];
	      edgeflag[i2] = i;
	      edgenr[i2] = ednr;
	    }
	  for (int j = 1; j <= vert2vertcoarse.EntrySize(i); j++)
	    {
	      int v2 = vert2vertcoarse.Get(i,j);
	      if (edgeflag[v2] < i)
		{
		  ned++;
		  edgenr[v2] = ned;
		  edgeflag[v2] = i;
		  missing.Append (INDEX_3(i,v2,ned));
		}
	    }

	  for (int j = 1; j <= vert2element->EntrySize(i); j++)
	    {
	      int elnr = vert2element->Get(i,j);
	      const Element & el = mesh.VolumeElement (elnr);

	      int neledges = GetNEdges (el.GetType());
	      const ELEMENT_EDGE * eledges = GetEdges (el.GetType());
	  
	      for (int k = 0; k < neledges; k++)
		{
		  INDEX_2 edge(el.PNum(eledges[k][0]),
			       el.PNum(eledges[k][1]));
	      
		  int edgedir = (edge.I1() > edge.I2());
		  if (edgedir) swap (edge.I1(), edge.I2());
	     
		  if (edge.I1() != i)
		    continue;
	     
		  if (edgeflag[edge.I2()] < i)
		    {
		      ned++;
		      edgenr[edge.I2()] = ned;
		      edgeflag[edge.I2()] = i;
		    }

		  int edgenum = edgenr[edge.I2()];
		  if (edgedir) edgenum *= -1;
		  edges.Elem(elnr)[k] = edgenum;
		}
	    }

	  for (int j = 1; j <= vert2surfelement->EntrySize(i); j++)
	    {
	      int elnr = vert2surfelement->Get(i,j);
	      const Element2d & el = mesh.SurfaceElement (elnr);

	      int neledges = GetNEdges (el.GetType());
	      const ELEMENT_EDGE * eledges = GetEdges (el.GetType());
	  
	      for (int k = 0; k < neledges; k++)
		{
		  INDEX_2 edge(el.PNum(eledges[k][0]),
			       el.PNum(eledges[k][1]));
	      
		  int edgedir = (edge.I1() > edge.I2());
		  if (edgedir) swap (edge.I1(), edge.I2());
	     
		  if (edge.I1() != i)
		    continue;
	     
		  if (edgeflag[edge.I2()] < i)
		    {
		      ned++;
		      edgenr[edge.I2()] = ned;
		      edgeflag[edge.I2()] = i;
		    }
	      
		  int edgenum = edgenr[edge.I2()];
		  if (edgedir) edgenum *= -1;
		  surfedges.Elem(elnr)[k] = edgenum;
		}
	    }

	  for (int j = 1; j <= vert2segment->EntrySize(i); j++)
	    {
	      int elnr = vert2segment->Get(i,j);
	      const Segment & el = mesh.LineSegment (elnr);

	      INDEX_2 edge(el.p1, el.p2);
	      
	      int edgedir = (edge.I1() > edge.I2());
	      if (edgedir) swap (edge.I1(), edge.I2());
	      
	      if (edge.I1() != i)
		continue;
	     
	      if (edgeflag[edge.I2()] < i)
		{
		  ned++;
		  edgenr[edge.I2()] = ned;
		  edgeflag[edge.I2()] = i;
		}   
 	      int edgenum = edgenr[edge.I2()];

	      if (edgedir) edgenum *= -1;
	      segedges.Elem(elnr) = edgenum;
	    }
	}


      edge2vert.SetSize (ned);
      for (int i = 1; i <= ne; i++)
	{
	  const Element & el = mesh.VolumeElement (i);
      
	  int neledges = GetNEdges (el.GetType());
	  const ELEMENT_EDGE * eledges = GetEdges (el.GetType());
	  
	  for (int k = 0; k < neledges; k++)
	    {
	      INDEX_2 edge(el.PNum(eledges[k][0]),
			   el.PNum(eledges[k][1]));
	  
	      int edgedir = (edge.I1() > edge.I2());
	      if (edgedir) swap (edge.I1(), edge.I2());

	      int edgenum = abs (edges.Elem(i)[k]);

	      edge2vert.Elem(edgenum)[0] = edge.I1();
	      edge2vert.Elem(edgenum)[1] = edge.I2();
	    }
	}
      for (int i = 1; i <= nse; i++)
	{
	  const Element2d & el = mesh.SurfaceElement (i);
      
	  int neledges = GetNEdges (el.GetType());
	  const ELEMENT_EDGE * eledges = GetEdges (el.GetType());
	  
	  for (int k = 0; k < neledges; k++)
	    {
	      INDEX_2 edge(el.PNum(eledges[k][0]),
			   el.PNum(eledges[k][1]));
	  
	      int edgedir = (edge.I1() > edge.I2());
	      if (edgedir) swap (edge.I1(), edge.I2());

	      int edgenum = abs (surfedges.Elem(i)[k]);

	      edge2vert.Elem(edgenum)[0] = edge.I1();
	      edge2vert.Elem(edgenum)[1] = edge.I2();
	    }
	}

      for (int i = 1; i <= nseg; i++)
	{
	  const Segment & el = mesh.LineSegment (i);
      
	  INDEX_2 edge(el.p1, el.p2);
	  int edgedir = (edge.I1() > edge.I2());
	  if (edgedir) swap (edge.I1(), edge.I2());
	  
	  int edgenum = abs (segedges.Elem(i));
	  
	  edge2vert.Elem(edgenum)[0] = edge.I1();
	  edge2vert.Elem(edgenum)[1] = edge.I2();
	}

      for (int i = 1; i <= missing.Size(); i++)
	{
	  INDEX_3 i3 = missing.Get(i);
	  edge2vert.Elem(i3.I3())[0] = i3.I1();
	  edge2vert.Elem(i3.I3())[1] = i3.I2();
	}
	
      
      /*
	(*testout) << "edge table:" << endl;
	(*testout) << "edge2vert:" << endl;
	for (int i = 1; i <= edge2vert.Size(); i++)
	(*testout) << "edge " << i << ", v1,2 = " << edge2vert.Elem(i)[0] << ", " << edge2vert.Elem(i)[1] << endl;
	(*testout) << "surfedges:" << endl;
	for (int i = 1; i <= surfedges.Size(); i++)
	(*testout) << "el " << i << ", edges = " 
	<< surfedges.Elem(i)[0] << ", "
	<< surfedges.Elem(i)[1] << ", "
	<< surfedges.Elem(i)[2] << endl;
      */ 
     
    
    }


  //  cout << "build edges done" << endl;

  // generate faces
  if (buildfaces) //  && mesh.GetDimension() == 3)
    {
      int i, j;

      static int timer2 = NgProfiler::CreateTimer ("topology::buildfaces");
      NgProfiler::RegionTimer reg2 (timer2);

      PrintMessage (5, "Update faces ");

      faces.SetSize(ne);
      surffaces.SetSize(nse);
      
      // face2vert.SetSize(0);  // keep old faces
      nfa = face2vert.Size();
      // INDEX_3_HASHTABLE<int> vert2face(ne+nse+1);
      INDEX_3_CLOSED_HASHTABLE<int> vert2face(8*ne+2*nse+nfa+2);

      for (i = 1; i <= face2vert.Size(); i++)
	{
	  INDEX_3 f;
	  f.I1() = face2vert.Get(i)[0];
	  f.I2() = face2vert.Get(i)[1];
	  f.I3() = face2vert.Get(i)[2];
	  vert2face.Set (f, i);
	}
     
      for (i = 1; i <= ne; i++)
	{
	  const Element & el = mesh.VolumeElement (i);
	  
	  int nelfaces = GetNFaces (el.GetType());
	  const ELEMENT_FACE * elfaces = GetFaces (el.GetType());
	  
	  for (j = 0; j < 6; j++)
	    faces.Elem(i)[j] = 0;
	  for (j = 0; j < nelfaces; j++)
	    if (elfaces[j][3] == 0)
	      
	      { // triangle
		
		int facenum;
		int facedir;
		
		INDEX_3 face(el.PNum(elfaces[j][0]),
			     el.PNum(elfaces[j][1]),
			     el.PNum(elfaces[j][2]));
		
		facedir = 0;
		if (face.I1() > face.I2())
		  {
		    swap (face.I1(), face.I2());
		    facedir += 1;
		  }
		if (face.I2() > face.I3())
		  {
		    swap (face.I2(), face.I3());
		    facedir += 2;
		  }
		if (face.I1() > face.I2())
		  {
		    swap (face.I1(), face.I2());
		    facedir += 4;
		  }
		
		if (vert2face.Used (face))
		  facenum = vert2face.Get(face);
		else
		  {
		    nfa++;
		    vert2face.Set (face, nfa);
		    facenum = nfa;
		    
		    INDEX_4 hface(face.I1(),face.I2(),face.I3(),0);
		    face2vert.Append (hface);
		    // face2vert.SetSize(face2vert.Size()+1);
		  }
		
		faces.Elem(i)[j] = 8*(facenum-1)+facedir+1;
	      }
	  
	    else
	      
	      {
		// quad
		int facenum;
		int facedir;
		INDEX_4Q face4(el.PNum(elfaces[j][0]),
			       el.PNum(elfaces[j][1]),
			       el.PNum(elfaces[j][2]),
			       el.PNum(elfaces[j][3]));
		
		facedir = 0;
		if (min2 (face4.I1(), face4.I2()) > 
		    min2 (face4.I4(), face4.I3())) 
		  {  // z - flip
		    facedir += 1; 
		    swap (face4.I1(), face4.I4());
		    swap (face4.I2(), face4.I3());
		  }
		if (min2 (face4.I1(), face4.I4()) >
		    min2 (face4.I2(), face4.I3())) 
		  {  // x - flip
		    facedir += 2; 
		    swap (face4.I1(), face4.I2());
		    swap (face4.I3(), face4.I4());
		  }
		if (face4.I2() > face4.I4())
		  {  // diagonal flip
		    facedir += 4; 
		    swap (face4.I2(), face4.I4());
		  }
		//		face4.Sort();
		
		INDEX_3 face(face4.I1(), face4.I2(), face4.I3());
		
		if (vert2face.Used (face))
		  {
		    facenum = vert2face.Get(face);
		  }
		else
		  {
		    nfa++;
		    vert2face.Set (face, nfa);
		    facenum = nfa;

		    // face2vert.SetSize(face2vert.Size()+1);

		    INDEX_4 hface(face4.I1(),face4.I2(),face4.I3(),face4.I4());
		    face2vert.Append (hface);
		  }
		
		faces.Elem(i)[j] = 8*(facenum-1)+facedir+1;
	      }
	}

      face2surfel.SetSize(nfa+nse);
      for (i = 1; i <= face2surfel.Size(); i++)
	face2surfel.Elem(i) = 0;

      for (i = 1; i <= nse; i++)
	{
	  const Element2d & el = mesh.SurfaceElement (i);
	  
	  const ELEMENT_FACE * elfaces = GetFaces (el.GetType());
	  
	  if (elfaces[0][3] == 0)
	    
	    { // triangle
	      
	      int facenum;
	      int facedir;
	      
	      INDEX_3 face(el.PNum(elfaces[0][0]),
			   el.PNum(elfaces[0][1]),
			   el.PNum(elfaces[0][2]));
	      
	      facedir = 0;
	      if (face.I1() > face.I2())
		{
		  swap (face.I1(), face.I2());
		  facedir += 1;
		}
	      if (face.I2() > face.I3())
		{
		  swap (face.I2(), face.I3());
		  facedir += 2;
		}
	      if (face.I1() > face.I2())
		{
		  swap (face.I1(), face.I2());
		  facedir += 4;
		}
	      
	      if (vert2face.Used (face))
		facenum = vert2face.Get(face);
	      else
		{
		  nfa++;
		  vert2face.Set (face, nfa);
		  facenum = nfa;
		  
		  // face2vert.SetSize(face2vert.Size()+1);
		  INDEX_4 hface(face.I1(),face.I2(),face.I3(),0);
		  face2vert.Append (hface);
		}
	      
	      surffaces.Elem(i) = 8*(facenum-1)+facedir+1;
	      face2surfel.Elem(facenum) = i;
	    }
	  
	  else
	    
	    {
	      // quad
	      int facenum;
	      int facedir;
	      
	      INDEX_4Q face4(el.PNum(elfaces[0][0]),
			     el.PNum(elfaces[0][1]),
			     el.PNum(elfaces[0][2]),
			     el.PNum(elfaces[0][3]));

	      facedir = 0;
	      if (min2 (face4.I1(), face4.I2()) > 
		  min2 (face4.I4(), face4.I3())) 
		{  // z - orientation
		  facedir += 1; 
		  swap (face4.I1(), face4.I4());
		  swap (face4.I2(), face4.I3());
		}
	      if (min2 (face4.I1(), face4.I4()) >
		  min2 (face4.I2(), face4.I3())) 
		{  // x - orientation
		  facedir += 2; 
		  swap (face4.I1(), face4.I2());
		  swap (face4.I3(), face4.I4());
		}
	      if (face4.I2() > face4.I4())
		{ 
		  facedir += 4; 
		  swap (face4.I2(), face4.I4());
		}
	      
	      INDEX_3 face(face4.I1(), face4.I2(), face4.I3());
	      
	      if (vert2face.Used (face))
		facenum = vert2face.Get(face);
	      else
		{
		  nfa++;
		  vert2face.Set (face, nfa);
		  facenum = nfa;
		  
		  // face2vert.SetSize(face2vert.Size()+1);
		  INDEX_4 hface(face4.I1(),face4.I2(),face4.I3(),face4.I3());
		  face2vert.Append (hface);
		  /*
		  face2vert.Last()[0] = face4.I1();
		  face2vert.Last()[1] = face4.I2();
		  face2vert.Last()[2] = face4.I3();
		  face2vert.Last()[3] = face4.I3();
		  */
		}
	      
	      surffaces.Elem(i) = 8*(facenum-1)+facedir+1;
	      face2surfel.Elem(facenum) = i;
	    }
	}


      surf2volelement.SetSize (nse);
      for (i = 1; i <= nse; i++)
	{
	  surf2volelement.Elem(i)[0] = 0;
	  surf2volelement.Elem(i)[1] = 0;
	}
      for (i = 1; i <= ne; i++)
	for (j = 0; j < 6; j++)
	  {
            int fnum = (faces.Get(i)[j]+7) / 8;
	    if (fnum > 0 && face2surfel.Elem(fnum))
	      {
		int sel = face2surfel.Elem(fnum);
		surf2volelement.Elem(sel)[1] = 
		  surf2volelement.Elem(sel)[0];
		surf2volelement.Elem(sel)[0] = i;
	      }
	  }

      face2vert.SetAllocSize (face2vert.Size());

      /*
	*testout << "face2vert: ";
	face2vert.PrintMemInfo(cout);
	*testout  << "faces: ";
	faces.PrintMemInfo(cout);
	*testout << "hashtable: ";
	vert2face.PrintMemInfo(cout);
      */

#ifdef PARALLEL
  (*testout) << " RESET Paralleltop" << endl;

      paralleltop.Reset ();
#endif

      ARRAY<char> face_els(nfa), face_surfels(nfa);
      face_els = 0;
      face_surfels = 0;
      ARRAY<int> hfaces;
      for (i = 1; i <= ne; i++)
	{
	  GetElementFaces (i, hfaces);
	  for (j = 0; j < hfaces.Size(); j++)
	    face_els[hfaces[j]-1]++;
	}
      for (i = 1; i <= nse; i++)
	face_surfels[GetSurfaceElementFace (i)-1]++;

      if (ne)
	{
	  int cnt_err = 0;
	  for (i = 0; i < nfa; i++)
	    {
	      /*
	      (*testout) << "face " << i << " has " << int(face_els[i]) << " els, " 
 			 << int(face_surfels[i]) << " surfels, tot = "
 			 << face_els[i] + face_surfels[i] << endl; 
	      */
	      if (face_els[i] + face_surfels[i] == 1)
		{
		  cnt_err++;
#ifdef PARALLEL
		  if ( ntasks > 1 )
		    {
		  if ( !paralleltop.DoCoarseUpdate() ) continue;
		  // "illegal" faces are exchange faces
		  /*
		  (*testout) << "exchange face : " << i << endl;
		  (*testout) << "points = " << face2vert[i] << endl;
		  */
// 		  (*testout) << "global points = ";
// 		  for ( int j = 0; j < 3; j++ )
// // 		    (*testout) << face2vert[i].I(j+1) << " -> " 
// // 			       << paralleltop.GetLoc2Glob_Vert( face2vert[i].I(j+1) ) << ",  ";
// 		  (*testout) << endl;
// 		  if ( !paralleltop.IsExchangeFace (i+1) )
// 		    paralleltop.SetRefinementFace (i+1);

		  paralleltop.SetExchangeFace (i+1);
		  
		  for (int j = 0; j < 4; j++)		    
		    {
		      if ( face2vert[i].I(j+1) > 0 )
			paralleltop.SetExchangeVert(face2vert[i].I(j+1));
		    }
		  
		  ARRAY<int> faceedges;
		  GetFaceEdges (i+1, faceedges);
		  for ( int j = 0; j < faceedges.Size(); j++)
		    {
		      paralleltop.SetExchangeEdge ( faceedges[j] );
		      int v1, v2;
		      GetEdgeVertices(faceedges[j], v1, v2 );
		    }
		  
		  /*
		  (*testout) << "pos = ";
		  for (int j = 0; j < 4; j++)
		    if (face2vert[i].I(j+1) >= 1)
		      (*testout) << mesh[(PointIndex)face2vert[i].I(j+1)] << " ";
		  (*testout) << endl;
		  */
		    }
		  else
		    {
#endif
		  (*testout) << "illegal face : " << i << endl;
		  (*testout) << "points = " << face2vert[i] << endl;
		  (*testout) << "pos = ";
		  for (j = 0; j < 4; j++)
		    if (face2vert[i].I(j+1) >= 1)
		      (*testout) << mesh[(PointIndex)face2vert[i].I(j+1)] << " ";
		  (*testout) << endl;

		  FlatArray<int> vertels = GetVertexElements (face2vert[i].I(1));
		  for (int k = 0; k < vertels.Size(); k++)
		    {
		      int elfaces[10], orient[10];
		      int nf = GetElementFaces (vertels[k], elfaces, orient);
		      for (int l = 0; l < nf; l++)
			if (elfaces[l] == i)
			  {
			    (*testout) << "is face of element " << vertels[k] << endl;
			    
			    if (mesh.coarsemesh && mesh.hpelements->Size() == mesh.GetNE() )
			      {
				const HPRefElement & hpref_el =
				  (*mesh.hpelements) [ mesh.VolumeElement (vertels[k]).hp_elnr];
				(*testout) << "coarse eleme = " << hpref_el.coarse_elnr << endl;
			      }

			  }
		    }
#ifdef PARALLEL
		    }
#endif
		}
	    }

#ifndef PARALLEL
	  if (cnt_err)
	    cout << cnt_err << " elements are not matching !!!" << endl;
#else
	  if (cnt_err && ntasks == 1)
	    cout << cnt_err << " elements are not matching !!!" << endl;
	  else if (cnt_err && ntasks > 1)
	    {
	      cout << "p" << id << ":  " << cnt_err << " elements are not local" << endl;
	      isparallel = 1;
	    }
	  else if ( ntasks > 1 )
	    cout << "p" << id << ":  " << "Partition " << id << " is totally local" << endl;
#endif

	}

#ifdef PARALLEL
     
      if ( isparallel )
	{
 	  paralleltop.Update();
	  if ( paralleltop.DoCoarseUpdate() )
	    {
	      paralleltop.UpdateCoarseGrid();
	    }
	  else
	    {
	      //  paralleltop.UpdateRefinement();
	    }
	  // paralleltop.Print();
	}

 
#endif




    }
 
 
  
  /* 
for (i = 1; i <= ne; i++)
    {
    (*testout) << "Element " << i << endl;
    (*testout) << "PNums " << endl; 
    for( int l=1;l<=8;l++) *testout << mesh.VolumeElement(i).PNum(l) << "\t"; 
    *testout << endl; 
    (*testout) << "edges: " << endl;
    for (j = 0; j < 9; j++)
    (*testout) << edges.Elem(i)[j] << " ";
    (*testout) << "faces: " << endl;
    for (j = 0; j < 6; j++)m
    (*testout) << faces.Elem(i)[j] << " ";
    }

    for (i = 1; i <= nse; i++)
    {
    (*testout) << "SElement " << i << endl;
    (*testout) << "PNums " << endl; 
    for( int l=1;l<=4;l++) *testout << mesh.SurfaceElement(i).PNum(l) << "\t"; 
    *testout << endl; 
    }
  */
  timestamp = NextTimeStamp();
}

  


int MeshTopology :: GetNVertices (ELEMENT_TYPE et)
{
  switch (et)
    {
    case SEGMENT:
    case SEGMENT3:
      return 2;

    case TRIG:
    case TRIG6:
      return 3;

    case QUAD:
    case QUAD6:
    case QUAD8:
      return 4;

    case TET:
    case TET10:
      return 4;

    case PYRAMID:
      return 5;

    case PRISM:
    case PRISM12:
      return 6;

    case HEX:
      return 8;

    default:
      cerr << "Ng_ME_GetNVertices, illegal element type " << et << endl;
    }
  return 0;
}

int MeshTopology :: GetNEdges (ELEMENT_TYPE et)
{
  switch (et)
    {
    case SEGMENT:
    case SEGMENT3:
      return 1;

    case TRIG:
    case TRIG6:
      return 3;

    case QUAD:
    case QUAD6:
    case QUAD8:
      return 4;

    case TET:
    case TET10:
      return 6;

    case PYRAMID:
      return 8;

    case PRISM:
    case PRISM12:
      return 9;

    case HEX:
      return 12;

    default:
      cerr << "Ng_ME_GetNEdges, illegal element type " << et << endl;
    }
  return 0;
}


int MeshTopology :: GetNFaces (ELEMENT_TYPE et)
{
  switch (et)
    {
    case SEGMENT:
    case SEGMENT3:
      return 0;

    case TRIG:
    case TRIG6:
      return 1;

    case QUAD:
    case QUAD6:
    case QUAD8:
      return 1;

    case TET:
    case TET10:
      return 4;

    case PYRAMID:
      return 5;

    case PRISM:
    case PRISM12:
      return 5;

    case HEX:
      return 6;

    default:
      cerr << "Ng_ME_GetNVertices, illegal element type " << et << endl;
    }
  return 0;
}




const Point3d * MeshTopology :: GetVertices (ELEMENT_TYPE et)
{
  static Point3d segm_points [] = 
    { Point3d (1, 0, 0),
      Point3d (0, 0, 0) };
  
  static Point3d trig_points [] = 
    { Point3d ( 1, 0, 0 ),
      Point3d ( 0, 1, 0 ),
      Point3d ( 0, 0, 0 ) };

  static Point3d quad_points [] = 
    { Point3d ( 0, 0, 0 ),
      Point3d ( 1, 0, 0 ),
      Point3d ( 1, 1, 0 ),
      Point3d ( 0, 1, 0 ) };

  static Point3d tet_points [] = 
    { Point3d ( 1, 0, 0 ),
      Point3d ( 0, 1, 0 ),
      Point3d ( 0, 0, 1 ),
      Point3d ( 0, 0, 0 ) };

  static Point3d pyramid_points [] =
    {
      Point3d ( 0, 0, 0 ),
      Point3d ( 1, 0, 0 ),
      Point3d ( 1, 1, 0 ),
      Point3d ( 0, 1, 0 ),
      Point3d ( 0, 0, 1-1e-7 ),
    };    
  
  static Point3d prism_points[] = 
    {
      Point3d ( 1, 0, 0 ),
      Point3d ( 0, 1, 0 ),
      Point3d ( 0, 0, 0 ),
      Point3d ( 1, 0, 1 ),
      Point3d ( 0, 1, 1 ),
      Point3d ( 0, 0, 1 )
    };


  static Point3d hex_points [] = 
    { Point3d ( 0, 0, 0 ),
      Point3d ( 1, 0, 0 ),
      Point3d ( 1, 1, 0 ),
      Point3d ( 0, 1, 0 ),
      Point3d ( 0, 0, 1 ),
      Point3d ( 1, 0, 1 ),
      Point3d ( 1, 1, 1 ),
      Point3d ( 0, 1, 1 ) };


  switch (et)
    {
    case SEGMENT:
    case SEGMENT3:
      return segm_points;

    case TRIG:
    case TRIG6:
      return trig_points;

    case QUAD:
    case QUAD6:
    case QUAD8:
      return quad_points;

    case TET:
    case TET10:
      return tet_points;

    case PYRAMID:
      return pyramid_points;

    case PRISM:
    case PRISM12:
      return prism_points;

    case HEX:
      return hex_points;
    default:
      cerr << "Ng_ME_GetVertices, illegal element type " << et << endl;
    }
  return 0;
}








void MeshTopology :: GetElementEdges (int elnr, ARRAY<int> & eledges) const
{
  int ned = GetNEdges (mesh.VolumeElement(elnr).GetType());
  eledges.SetSize (ned);
  for (int i = 0; i < ned; i++)
    eledges[i] = abs (edges.Get(elnr)[i]);
}
void MeshTopology :: GetElementFaces (int elnr, ARRAY<int> & elfaces, bool withorientation) const
{
  int i;
  int nfa = GetNFaces (mesh.VolumeElement(elnr).GetType());
  elfaces.SetSize (nfa);
  for (i = 1; i <= nfa; i++)
    {
      elfaces.Elem(i) = (faces.Get(elnr)[i-1]-1) / 8 + 1;
      if(withorientation)
	{
	  int orient = (faces.Get(elnr)[i-1]-1) % 8;
	  if(orient == 1 || orient == 2 || orient == 4 || orient == 7)
	    elfaces.Elem(i) *= -1;
	}
    }
}

void MeshTopology :: GetElementEdgeOrientations (int elnr, ARRAY<int> & eorient) const
{
  int i;
  int ned = GetNEdges (mesh.VolumeElement(elnr).GetType());
  eorient.SetSize (ned);
  for (i = 1; i <= ned; i++)
    eorient.Elem(i) = (edges.Get(elnr)[i-1] > 0) ? 1 : -1;
}

void MeshTopology :: GetElementFaceOrientations (int elnr, ARRAY<int> & forient) const
{
  int i;
  int nfa = GetNFaces (mesh.VolumeElement(elnr).GetType());
  forient.SetSize (nfa);
  for (i = 1; i <= nfa; i++)
    forient.Elem(i) = (faces.Get(elnr)[i-1]-1) % 8;
}



int MeshTopology :: GetElementEdges (int elnr, int * eledges, int * orient) const
{
  int i;
  //  int ned = GetNEdges (mesh.VolumeElement(elnr).GetType());

  if (mesh.GetDimension()==3 || 1)
    {
      if (orient)
	{
	  for (i = 0; i < 12; i++)
	    {
	      if (!edges.Get(elnr)[i]) return i;
	      eledges[i] = abs (edges.Get(elnr)[i]);
	      orient[i] = (edges.Get(elnr)[i] > 0 ) ? 1 : -1;
	    }
	}
      else
	{
	  for (i = 0; i < 12; i++)
	    {
	      if (!edges.Get(elnr)[i]) return i;
	      eledges[i] = abs (edges.Get(elnr)[i]);
	    }
	}
      return 12;
    }
  else
    {
		throw NgException("rethink implementation");
		/*
      if (orient)
	{
	  for (i = 0; i < 4; i++)
	    {
	      if (!surfedges.Get(elnr)[i]) return i;
	      eledges[i] = abs (surfedges.Get(elnr)[i]);
	      orient[i] = (surfedges.Get(elnr)[i] > 0 ) ? 1 : -1;
	    }
	}
      else
	{
	  if (!surfedges.Get(elnr)[i]) return i;
	  for (i = 0; i < 4; i++)
	    eledges[i] = abs (surfedges.Get(elnr)[i]);
	}
	*/
      return 4;
      //      return GetSurfaceElementEdges (elnr, eledges, orient);
    }
}

int MeshTopology :: GetElementFaces (int elnr, int * elfaces, int * orient) const
{
  int i;
  //  int nfa = GetNFaces (mesh.VolumeElement(elnr).GetType());
  if (orient)
    {
      for (i = 0; i < 6; i++)
	{
	  if (!faces.Get(elnr)[i]) return i;
	  elfaces[i] = (faces.Get(elnr)[i]-1) / 8 + 1;
	  orient[i] = (faces.Get(elnr)[i]-1) % 8;
	}
    }
  else
    {
      for (i = 0; i < 6; i++)
	{
	  if (!faces.Get(elnr)[i]) return i;
	  elfaces[i] = (faces.Get(elnr)[i]-1) / 8 + 1;
	}
    }
  return 6;
}

void MeshTopology :: GetSurfaceElementEdges (int elnr, ARRAY<int> & eledges) const
{
  int i;
  if (mesh.GetDimension()==3 || 1)
    {
      int ned = GetNEdges (mesh.SurfaceElement(elnr).GetType());
      eledges.SetSize (ned);
      for (i = 1; i <= ned; i++)
	eledges.Elem(i) = abs (surfedges.Get(elnr)[i-1]);
    }
  else
    {
      cout << "surfeledge(" << elnr << ") = " << flush;
      eledges.SetSize(1); 
      eledges.Elem(1) = abs (segedges.Get(elnr));
      cout << eledges.Elem(1) << endl;
    }
}

int MeshTopology :: GetSurfaceElementFace (int elnr) const
{
  return (surffaces.Get(elnr)-1) / 8 + 1;  
}

void MeshTopology :: 
GetSurfaceElementEdgeOrientations (int elnr, ARRAY<int> & eorient) const
{
  int ned = GetNEdges (mesh.SurfaceElement(elnr).GetType());
  eorient.SetSize (ned);
  for (int i = 1; i <= ned; i++)
    eorient.Elem(i) = (surfedges.Get(elnr)[i-1] > 0) ? 1 : -1;
}

int MeshTopology :: GetSurfaceElementFaceOrientation (int elnr) const
{
  return (surffaces.Get(elnr)-1) % 8;
}

int MeshTopology :: GetSurfaceElementEdges (int elnr, int * eledges, int * orient) const
{
  int i;
  if (mesh.GetDimension() == 3 || 1)
    {
      if (orient)
	{
	  for (i = 0; i < 4; i++)
	    {
	      if (!surfedges.Get(elnr)[i]) return i;
	      eledges[i] = abs (surfedges.Get(elnr)[i]);
	      orient[i] = (surfedges.Get(elnr)[i] > 0 ) ? 1 : -1;
	    }
	}
      else
	{
	  for (i = 0; i < 4; i++)
	    {
	      if (!surfedges.Get(elnr)[i]) return i;
	      eledges[i] = abs (surfedges.Get(elnr)[i]);
	    }
	}
      return 4;
    }
  else
    {
      eledges[0] = abs (segedges.Get(elnr));
      if (orient)
	orient[0] = segedges.Get(elnr) > 0 ? 1 : -1;
    }
  return 1;
}


void MeshTopology :: GetFaceVertices (int fnr, ARRAY<int> & vertices) const
{
  vertices.SetSize(4);
  int i;
  for (i = 1; i <= 4; i++)
    vertices.Elem(i) = face2vert.Get(fnr)[i-1];
  if (vertices.Elem(4) == 0)
    vertices.SetSize(3);
}

void MeshTopology :: GetFaceVertices (int fnr, int * vertices) const
{
  for (int i = 0; i <= 3; i++)
    vertices[i] = face2vert.Get(fnr)[i];
}


void MeshTopology :: GetEdgeVertices (int ednr, int & v1, int & v2) const
{
  v1 = edge2vert.Get(ednr)[0];
  v2 = edge2vert.Get(ednr)[1];
}


void MeshTopology :: GetFaceEdges (int fnr, ARRAY<int> & fedges, bool withorientation) const
{
  ArrayMem<int,4> pi(4);
  ArrayMem<int,12> eledges;
  
  fedges.SetSize (0);
  GetFaceVertices(fnr, pi);

  // Sort Edges according to global vertex numbers 
  // e1 = fmax, f2 
  // e2 = fmax, f1 
  // e3 = op e1(f2,f3) 
  // e4 = op e2(f1,f3) 

  /*  ArrayMem<int,4> fp; 
  fp[0] = pi[0]; 
  for(int k=1;k<pi.Size();k++) 
    if(fp[k]>fp[0]) swap(fp[k],fp[0]); 
  
    fp[1] = fp[0]+ */ 
  

  //  GetVertexElements (pi[0], els);
  FlatArray<int> els= GetVertexElements (pi[0]);

  // find one element having all vertices of the face
  for (int i = 0; i < els.Size(); i++)
    {
      const Element & el = mesh.VolumeElement(els[i]);
      int nref_faces = GetNFaces (el.GetType());
      const ELEMENT_FACE * ref_faces = GetFaces (el.GetType());
      int nfa_ref_edges = GetNEdges (GetFaceType(fnr));
      
      int cntv = 0,fa=-1; 
      for(int m=0;m<nref_faces;m++)
	{ 
	  cntv=0;
	  for(int j=0;j<nfa_ref_edges && ref_faces[m][j]>0;j++)
	    for(int k=0;k<pi.Size();k++)
	      {
		if(el[ref_faces[m][j]-1] == pi[k])
		  cntv++;
	      }
	  if (cntv == pi.Size())
	    {
	      fa=m;
	      break;
	    }
	}
     
      if(fa>=0)
	{
	  const ELEMENT_EDGE * fa_ref_edges = GetEdges(GetFaceType(fnr)); 
	  fedges.SetSize(nfa_ref_edges);
	  GetElementEdges (els[i], eledges);
	  
	  for (int j = 0; j < eledges.Size(); j++)
	    {
	      int vi1, vi2;
	      GetEdgeVertices (eledges[j], vi1, vi2);
	    
	      bool has1 = 0;
	      bool has2 = 0;
	      for (int k = 0; k < pi.Size(); k++)
		{
		  if (vi1 == pi[k]) has1 = 1;
		  if (vi2 == pi[k]) has2 = 1;
		  
		}
	      
	      if (has1 && has2) // eledges[j] is on face 
		{
		  // fedges.Append (eledges[j]);
		  for(int k=0;k<nfa_ref_edges;k++)
		    {
		      int w1 = el[ref_faces[fa][fa_ref_edges[k][0]-1]-1]; 
		      int w2 = el[ref_faces[fa][fa_ref_edges[k][1]-1]-1]; 

		      if(withorientation)
			{
			  if(w1==vi1 && w2==vi2)
			    fedges[k] = eledges[j];
			  if(w1==vi2 && w2==vi1)
			    fedges[k] = -eledges[j];
			}
		      else
			if((w1==vi1 && w2==vi2) || (w1==vi2 && w2==vi1))
			  fedges[k] = eledges[j];
		    }
		}
	    }
	  
	  // *testout << " Face " << fnr << endl; 
	  // *testout << " GetFaceEdges " << fedges << endl;
	  
	  return;
	}
    }   
}


ELEMENT_TYPE MeshTopology :: GetFaceType (int fnr) const
{
  if (face2vert.Get(fnr)[3] == 0) return TRIG; else return QUAD;
}


void MeshTopology :: GetVertexElements (int vnr, ARRAY<int> & elements) const
{
  if (vert2element)
    {
      int i; 
      int ne = vert2element->EntrySize(vnr);
      elements.SetSize(ne);
      for (i = 1; i <= ne; i++)
	elements.Elem(i) = vert2element->Get(vnr, i);
    }
}


FlatArray<int> MeshTopology :: GetVertexElements (int vnr) const
{
  if (vert2element)
    return (*vert2element)[vnr];
  return FlatArray<int> (0,0);
}

FlatArray<int> MeshTopology :: GetVertexSurfaceElements (int vnr) const
{
  if (vert2surfelement)
    return (*vert2surfelement)[vnr];
  return FlatArray<int> (0,0);
}


void MeshTopology :: GetVertexSurfaceElements( int vnr, 
					       ARRAY<int>& elements ) const
{
  if (vert2surfelement)
    {
      int i;
      int ne = vert2surfelement->EntrySize(vnr);
      elements.SetSize(ne);
      for (i = 1; i <= ne; i++)
	elements.Elem(i) = vert2surfelement->Get(vnr, i);
    }
}


int MeshTopology :: GetVerticesEdge ( int v1, int v2 ) const
{
  ARRAY<int> elements_v1, elementedges;
  GetVertexElements ( v1, elements_v1);
  int edv1, edv2;

  for ( int i = 0; i < elements_v1.Size(); i++ )
    {
      GetElementEdges( elements_v1[i], elementedges );
      for ( int ed = 0; ed < elementedges.Size(); ed ++)
	{
	  GetEdgeVertices( elementedges[ed], edv1, edv2 );
	  if ( ( edv1 == v1 && edv2 == v2 ) || ( edv1 == v2 && edv2 == v1 ) )
	    return elementedges[ed];
	}
    }

  return -1;
}



void MeshTopology :: GetSegmentVolumeElements ( int segnr, ARRAY<int> & volels ) const
{
  int v1, v2;
  GetEdgeVertices ( GetSegmentEdge (segnr), v1, v2 );
  ARRAY<int> volels1, volels2;
  GetVertexElements ( v1, volels1 );
  GetVertexElements ( v2, volels2 );
  volels.SetSize(0);

  for ( int eli1=1; eli1 <= volels1.Size(); eli1++)
    if ( volels2.Contains( volels1.Elem(eli1) ) )
      volels.Append ( volels1.Elem(eli1) );

}
}
