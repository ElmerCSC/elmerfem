#include <mystdlib.h>

#include "meshing.hpp"

#ifdef SOLIDGEOM
#include <csg.hpp>
#endif
#include <opti.hpp>

namespace netgen
{

/*
  Combine two points to one.
  Set new point into the center, if both are
  inner points.
  Connect inner point to boundary point, if one
  point is inner point.
*/

void MeshOptimize3d :: CombineImprove (Mesh & mesh,
				       OPTIMIZEGOAL goal)
{
  int np = mesh.GetNP();
  int ne = mesh.GetNE();

  TABLE<ElementIndex, PointIndex::BASE> elementsonnode(np); 
  ARRAY<ElementIndex> hasonepi, hasbothpi;

  ARRAY<double> oneperr;
  ARRAY<double> elerrs (ne);

  PrintMessage (3, "CombineImprove");
  (*testout)  << "Start CombineImprove" << "\n";

  //  mesh.CalcSurfacesOfNode ();
  const char * savetask = multithread.task;
  multithread.task = "Combine Improve";


  double totalbad = 0;
  for (ElementIndex ei = 0; ei < ne; ei++)
    {
      double elerr = CalcBad (mesh.Points(), mesh[ei], 0);
      totalbad += elerr;
      elerrs[ei] = elerr;
    }

  if (goal == OPT_QUALITY)
    {
      totalbad = CalcTotalBad (mesh.Points(), mesh.VolumeElements());
      (*testout) << "Total badness = " << totalbad << endl;
      PrintMessage (5, "Total badness = ", totalbad);
    }

  for (ElementIndex ei = 0; ei < ne; ei++)
    if (!mesh[ei].IsDeleted())
      for (int j = 0; j < mesh[ei].GetNP(); j++)
	elementsonnode.Add (mesh[ei][j], ei);
  
  INDEX_2_HASHTABLE<int> edgetested (np+1);

  int cnt = 0;

  for (ElementIndex ei = 0; ei < ne; ei++)
    {
      if (multithread.terminate)
	break;
      
      multithread.percent = 100.0 * (ei+1) / ne;

      if (mesh.ElementType(ei) == FIXEDELEMENT)
	continue;

      for (int j = 0; j < 6; j++)
	{
	  Element & elemi = mesh[ei];
	  if (elemi.IsDeleted()) continue;
	  
	  static const int tetedges[6][2] =
	    { { 0, 1 }, { 0, 2 }, { 0, 3 },
	      { 1, 2 }, { 1, 3 }, { 2, 3 } };

	  PointIndex pi1 = elemi[tetedges[j][0]];
	  PointIndex pi2 = elemi[tetedges[j][1]];

	  if (pi2 < pi1) Swap (pi1, pi2);
	  
	  INDEX_2 si2 (pi1, pi2);
	  si2.Sort();
	  
	  if (edgetested.Used (si2)) continue;
	  edgetested.Set (si2, 1);


	  // hasonepoint.SetSize(0);
	  //	  hasbothpoints.SetSize(0);
	  hasonepi.SetSize(0);
	  hasbothpi.SetSize(0);

	  FlatArray<ElementIndex> row1 = elementsonnode[pi1];
	  for (int k = 0; k < row1.Size(); k++)
	    {
	      Element & elem = mesh[row1[k]];
	      if (elem.IsDeleted()) continue;

	      if (elem[0] == pi2 || elem[1] == pi2 ||
		  elem[2] == pi2 || elem[3] == pi2)
		{
		  hasbothpi.Append (row1[k]);
		}
	      else
		{
		  hasonepi.Append (row1[k]);
		}
	    } 
	  
	  FlatArray<ElementIndex> row2 = elementsonnode[pi2];	  
	  for (int k = 0; k < row2.Size(); k++)
	    {
	      Element & elem = mesh[row2[k]];
	      if (elem.IsDeleted()) continue;

	      if (elem[0] == pi1 || elem[1] == pi1 ||
		  elem[2] == pi1 || elem[3] == pi1)
		;
	      else
		{
		  hasonepi.Append (row2[k]);
		}
	    } 
	  
	  double bad1 = 0;
	  for (int k = 0; k < hasonepi.Size(); k++)
	    bad1 += elerrs[hasonepi[k]];
	  for (int k = 0; k < hasbothpi.Size(); k++)
	    bad1 += elerrs[hasbothpi[k]];
	  
	  MeshPoint p1 = mesh[pi1];
	  MeshPoint p2 = mesh[pi2];
	  

	  // if (mesh.PointType(pi2) != INNERPOINT)
	  if (p2.Type() != INNERPOINT)
	    continue;
	  
	  MeshPoint pnew;
	  // if (mesh.PointType(pi1) != INNERPOINT)
	  if (p1.Type() != INNERPOINT)
	    pnew = p1;
	  else
	    pnew = Center (p1, p2);
	    
	  mesh[pi1] = pnew;
	  mesh[pi2] = pnew;

	  oneperr.SetSize (hasonepi.Size());

	  double bad2 = 0;
	  for (int k = 0; k < hasonepi.Size(); k++)
	    {
	      const Element & elem = mesh[hasonepi[k]];
	      double err = CalcTetBadness (mesh[elem[0]], mesh[elem[1]],  
					   mesh[elem[2]], mesh[elem[3]], 0);
	      bad2 += err;
	      oneperr[k] = err;
	    }
	  
	  mesh[pi1] = p1;
	  mesh[pi2] = p2;

	  // if (mesh.PointType(pi1) != INNERPOINT)
	  if (p1.Type() != INNERPOINT)
	    {
	      for (int k = 0; k < hasonepi.Size(); k++)
		{
		  Element & elem = mesh[hasonepi[k]];
		  int l;
		  for (l = 0; l < 4; l++)
		    if (elem[l] == pi2)
		      {
			elem[l] = pi1;
			break;
		      }
		      
		  elem.flags.illegal_valid = 0;
		  if (!mesh.LegalTet(elem))
		    bad2 += 1e4;
		  
		  if (l < 4)
		    {
		      elem.flags.illegal_valid = 0;
		      elem[l] = pi2;
		    }
		}
	    }

	  if (bad2 / hasonepi.Size()  < 
	      bad1 / (hasonepi.Size()+hasbothpi.Size()))
	    {
	      mesh[pi1] = pnew;
	      cnt++;

	      FlatArray<ElementIndex> row = elementsonnode[pi2];
	      for (int k = 0; k < row.Size(); k++)
		{
		  Element & elem = mesh[row[k]];
		  if (elem.IsDeleted()) continue;

		  elementsonnode.Add (pi1, row[k]);
		  for (int l = 0; l < elem.GetNP(); l++)
		    if (elem[l] == pi2)
		      elem[l] = pi1;
		  
		  elem.flags.illegal_valid = 0;
		  if (!mesh.LegalTet (elem))
		    (*testout) << "illegal tet " << elementsonnode[pi2][k] << endl;
		}

	      for (int k = 0; k < hasonepi.Size(); k++)
		elerrs[hasonepi[k]] = oneperr[k];
	      
	      for (int k = 0; k < hasbothpi.Size(); k++)
		{
		  mesh[hasbothpi[k]].flags.illegal_valid = 0;
		  mesh[hasbothpi[k]].Delete();
		}
	    }
	}
    }

  mesh.Compress();
  mesh.MarkIllegalElements();

  PrintMessage (5, cnt, " elements combined");
  (*testout) << "CombineImprove done" << "\n";

  totalbad = 0;
  for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
    totalbad += CalcBad (mesh.Points(), mesh[ei], 0);

  if (goal == OPT_QUALITY)
    {
      totalbad = CalcTotalBad (mesh.Points(), mesh.VolumeElements());
      (*testout) << "Total badness = " << totalbad << endl;

      int cntill = 0;
      for (ElementIndex ei = 0; ei < ne; ei++)
	if (!mesh.LegalTet (mesh[ei]))
	  cntill++;

      PrintMessage (5, cntill, " illegal tets");
    }
  multithread.task = savetask;
} 





/*
  Mesh improvement by edge splitting.
  If mesh quality is improved by inserting a node into an inner edge,
  the edge is split into two parts.
*/
void MeshOptimize3d :: SplitImprove (Mesh & mesh,
				     OPTIMIZEGOAL goal)
{
  int j, k, l;
  Point3d p1, p2, pnew;

  ElementIndex ei;
  SurfaceElementIndex sei;
  PointIndex pi1, pi2;

  double bad1, bad2, badmax, badlimit;


  int cnt = 0;
  int np = mesh.GetNP();
  int ne = mesh.GetNE();

  TABLE<ElementIndex,PointIndex::BASE> elementsonnode(np); 
  ARRAY<ElementIndex> hasbothpoints;

  BitArray origpoint(np), boundp(np);
  origpoint.Set();

  ARRAY<double> elerrs(ne);
  BitArray illegaltet(ne);
  illegaltet.Clear();

  const char * savetask = multithread.task;
  multithread.task = "Split Improve";


  PrintMessage (3, "SplitImprove");
  (*testout)  << "start SplitImprove" << "\n";

  ARRAY<INDEX_3> locfaces;

  INDEX_2_HASHTABLE<int> edgetested (np);

  bad1 = 0;
  badmax = 0;
  for (ei = 0; ei < ne; ei++)
    {
      elerrs[ei] = CalcBad (mesh.Points(), mesh[ei], 0);
      bad1 += elerrs[ei];
      if (elerrs[ei] > badmax) badmax = elerrs[ei];
    }

  PrintMessage (5, "badmax = ", badmax);
  badlimit = 0.5 * badmax;


  boundp.Clear();
  for (sei = 0; sei < mesh.GetNSE(); sei++)
    for (j = 0; j < 3; j++)
      boundp.Set (mesh[sei][j]);

  if (goal == OPT_QUALITY)
    {
      bad1 = CalcTotalBad (mesh.Points(), mesh.VolumeElements());
      (*testout) << "Total badness = " << bad1 << endl;
    }

  for (ei = 0; ei < ne; ei++)
    for (j = 0; j < mesh[ei].GetNP(); j++)
      elementsonnode.Add (mesh[ei][j], ei);


  mesh.MarkIllegalElements();
  if (goal == OPT_QUALITY || goal == OPT_LEGAL)
    {
      int cntill = 0;
      for (ei = 0; ei < ne; ei++)
	{
	  //	  if (!LegalTet (volelements.Get(i)))
	  if (mesh[ei].flags.illegal)
	    {
	      cntill++;
	      illegaltet.Set (ei+1);
	    }
	}
      //      (*mycout) << cntill << " illegal tets" << endl;
    }


  for (ei = 0; ei < ne; ei++)
    {
      if (multithread.terminate)
	break;

      multithread.percent = 100.0 * (ei+1) / ne;

      bool ltestmode = 0;


      if (elerrs[ei] < badlimit && !illegaltet.Test(ei+1)) continue;

      if ((goal == OPT_LEGAL) &&
	  !illegaltet.Test(ei+1) &&
	  CalcBad (mesh.Points(), mesh[ei], 0) < 1e3) 
	continue;

      
      Element & elem = mesh[ei];

      if (ltestmode)
	{
	  (*testout) << "test el " << ei << endl;
	  for (j = 0; j < 4; j++)
	    (*testout) << elem[j] << " ";
	  (*testout) << endl;
	}


      for (j = 0; j < 6; j++)
	{

	  static const int tetedges[6][2] =
	    { { 0, 1 }, { 0, 2 }, { 0, 3 },
	      { 1, 2 }, { 1, 3 }, { 2, 3 } };

	  pi1 = elem[tetedges[j][0]];
	  pi2 = elem[tetedges[j][1]];

	  if (pi2 < pi1) Swap (pi1, pi2);
	  if (pi2 > elementsonnode.Size()) continue;

	  if (!origpoint.Test(pi1) || !origpoint.Test(pi2))
	    continue;

  
	  INDEX_2 i2(pi1, pi2);
	  i2.Sort();

	  if (mesh.BoundaryEdge (pi1, pi2)) continue;

	  if (edgetested.Used (i2) && !illegaltet.Test(ei+1)) continue;
	  edgetested.Set (i2, 1);

	  hasbothpoints.SetSize (0);
	  for (k = 1; k <= elementsonnode.EntrySize(pi1); k++)
	    {
	      bool has1 = 0, has2 = 0;

	      ElementIndex elnr = elementsonnode.Get(pi1, k);
	      Element & el = mesh[elnr];

	      for (l = 0; l < el.GetNP(); l++)
		{
		  if (el[l] == pi1) has1 = 1;
		  if (el[l] == pi2) has2 = 1;
		}
	      if (has1 && has2) 
		{ // only once
		  for (l = 0; l < hasbothpoints.Size(); l++)
		    if (hasbothpoints[l] == elnr)
		      has1 = 0;
		  
		  if (has1)
		    hasbothpoints.Append (elnr);
		}
	    }
	  
	  bad1 = 0;
	  for (k = 0; k < hasbothpoints.Size(); k++)
	    bad1 += CalcBad (mesh.Points(), mesh[hasbothpoints[k]], 0);
	  

	  bool puretet = 1;
	  for (k = 0; k < hasbothpoints.Size(); k++)
	    if (mesh[hasbothpoints[k]].GetType() != TET)
	      puretet = 0;
	  if (!puretet) continue;

	  p1 = mesh[pi1];
	  p2 = mesh[pi2];

	  /*
	    pnew = Center (p1, p2);
	  
	    points.Elem(pi1) = pnew;
	    bad2 = 0;
	    for (k = 1; k <= hasbothpoints.Size(); k++)
	    bad2 += CalcBad (points, 
	    volelements.Get(hasbothpoints.Get(k)), 0);

	    points.Elem(pi1) = p1;
	    points.Elem(pi2) = pnew;
	      
	    for (k = 1; k <= hasbothpoints.Size(); k++)
	    bad2 += CalcBad (points, 
	    volelements.Get(hasbothpoints.Get(k)), 0);
	    points.Elem(pi2) = p2;
	  */


	  locfaces.SetSize (0);
	  for (k = 0; k < hasbothpoints.Size(); k++)
	    {
	      const Element & el = mesh[hasbothpoints[k]];

	      for (l = 0; l < 4; l++)
		if (el[l] == pi1 || el[l] == pi2)
		  {
		    INDEX_3 i3;
		    Element2d face;
		    el.GetFace (l+1, face);
		    for (int kk = 1; kk <= 3; kk++)
		      i3.I(kk) = face.PNum(kk);
		    locfaces.Append (i3);
		  }
	    }

	  PointFunction1 pf (mesh.Points(), locfaces, -1);
	  OptiParameters par;
	  par.maxit_linsearch = 50;
	  par.maxit_bfgs = 20;

	  pnew = Center (p1, p2);
	  Vector px(3);
	  px.Elem(1) = pnew.X();
	  px.Elem(2) = pnew.Y();
	  px.Elem(3) = pnew.Z();

	  if (elerrs[ei] > 0.1 * badmax)
	    BFGS (px, pf, par);

	  bad2 = pf.Func (px);

	  pnew.X() = px.Get(1);
	  pnew.Y() = px.Get(2);
	  pnew.Z() = px.Get(3);


	  int hpinew = mesh.AddPoint (pnew);
	  //	  ptyps.Append (INNERPOINT);

	  for (k = 0; k < hasbothpoints.Size(); k++)
	    {
	      Element & oldel = mesh[hasbothpoints[k]];
	      Element newel1 = oldel;
	      Element newel2 = oldel;

	      oldel.flags.illegal_valid = 0;
	      newel1.flags.illegal_valid = 0;
	      newel2.flags.illegal_valid = 0;
	      
	      for (l = 0; l < 4; l++)
		{
		  if (newel1[l] == pi2) newel1[l] = hpinew;
		  if (newel2[l] == pi1) newel2[l] = hpinew;
		}
	      
	      if (!mesh.LegalTet (oldel)) bad1 += 1e6;
	      if (!mesh.LegalTet (newel1)) bad2 += 1e6;
	      if (!mesh.LegalTet (newel2)) bad2 += 1e6;
	    }	  
	  
	  // mesh.PointTypes().DeleteLast();
	  mesh.Points().DeleteLast();

	  if (bad2 < bad1) 
	    /*	      (bad1 > 1e4 && boundp.Test(pi1) && boundp.Test(pi2)) */ 
	    {
	      cnt++;

	      PointIndex pinew = mesh.AddPoint (pnew);
	      
	      for (k = 0; k < hasbothpoints.Size(); k++)
		{
		  Element & oldel = mesh[hasbothpoints[k]];
		  Element newel = oldel;

		  newel.flags.illegal_valid = 0;
		  oldel.flags.illegal_valid = 0;

		  for (l = 0; l < 4; l++)
		    {
		      origpoint.Clear (oldel[l]);

		      if (oldel[l] == pi2) oldel[l] = pinew;
		      if (newel[l] == pi1) newel[l] = pinew;
		    }
		  mesh.AddVolumeElement (newel);
		}
	      
	      j = 10;
	    }
	}
    }


  mesh.Compress();
  PrintMessage (5, cnt, " splits performed");

  (*testout) << "Splitt - Improve done" << "\n";

  if (goal == OPT_QUALITY)
    {
      bad1 = CalcTotalBad (mesh.Points(), mesh.VolumeElements());
      (*testout) << "Total badness = " << bad1 << endl;

      int cntill = 0;
      ne = mesh.GetNE();
      for (ei = 0; ei < ne; ei++)
	{
	  if (!mesh.LegalTet (mesh[ei]))
	    cntill++;
	}
      //      cout << cntill << " illegal tets" << endl;
    }

  multithread.task = savetask;
}


      
  

void MeshOptimize3d :: SwapImprove (Mesh & mesh, OPTIMIZEGOAL goal,
				    const BitArray * working_elements)
{
  int j, k, l;

  ElementIndex ei;
  SurfaceElementIndex sei;

  PointIndex pi1(0), pi2(0), pi3(0), pi4(0), pi5(0), pi6(0);
  int cnt = 0;

  Element el21(TET), el22(TET), el31(TET), el32(TET), el33(TET);
  Element el1(TET), el2(TET), el3(TET), el4(TET);
  Element el1b(TET), el2b(TET), el3b(TET), el4b(TET);
  
  double bad1, bad2, bad3;

  int np = mesh.GetNP();
  int ne = mesh.GetNE();
  //int nse = mesh.GetNSE();
  
  // contains at least all elements at node
  TABLE<ElementIndex,PointIndex::BASE> elementsonnode(np);

  ARRAY<ElementIndex> hasbothpoints;

  PrintMessage (3, "SwapImprove ");
  (*testout) << "\n" << "Start SwapImprove" << endl;

  const char * savetask = multithread.task;
  multithread.task = "Swap Improve";
  
  //  mesh.CalcSurfacesOfNode ();
  /*
    for (i = 1; i <= GetNE(); i++)
    if (volelements.Get(i).PNum(1))
    if (!LegalTet (volelements.Get(i)))
    {
    cout << "detected illegal tet, 1" << endl;
    (*testout) << "detected illegal tet1: " << i << endl;
    }
  */	
  
    
  INDEX_3_HASHTABLE<int> faces(mesh.GetNOpenElements()/3 + 2);
  if (goal == OPT_CONFORM)
    {
      for (int i = 1; i <= mesh.GetNOpenElements(); i++)
	{
	  const Element2d & hel = mesh.OpenElement(i);
	  INDEX_3 face(hel[0], hel[1], hel[2]);
	  face.Sort();
	  faces.Set (face, 1);
	}
    }
  
  // Calculate total badness
  if (goal == OPT_QUALITY)
    {
      bad1 = CalcTotalBad (mesh.Points(), mesh.VolumeElements());
      (*testout) << "Total badness = " << bad1 << endl;
    }
  
  // find elements on node
  for (ei = 0; ei < ne; ei++)
    for (j = 0; j < mesh[ei].GetNP(); j++)
      elementsonnode.Add (mesh[ei][j], ei);

  /*
    BitArray illegaltet(GetNE());
    MarkIllegalElements();
    if (goal == OPT_QUALITY || goal == OPT_LEGAL)
    {
    int cntill = 0;
    for (i = 1; i <= GetNE(); i++)
    {
    //	  if (!LegalTet (volelements.Get(i)))
    if (VolumeElement(i).flags.illegal)
    {
    cntill++;
    illegaltet.Set (i);
    }
    }
    //      (*mycout) << cntill << " illegal tets" << endl;
    }
  */

  INDEX_2_HASHTABLE<int> edgeused(2 * ne + 5);

  for (ei = 0; ei < ne; ei++)
    {
      if (multithread.terminate)
	break;
      
      multithread.percent = 100.0 * (ei+1) / ne;

      if ((mesh.ElementType(ei)) == FIXEDELEMENT)
	continue;

      if(working_elements && 
	 ei < working_elements->Size() &&
	 !working_elements->Test(ei))
	continue;

      if (mesh[ei].IsDeleted())
	continue;

      if ((goal == OPT_LEGAL) && 
	  mesh.LegalTet (mesh[ei]) &&
	  CalcBad (mesh.Points(), mesh[ei], 0) < 1e3)
	continue;

      //      int onlybedges = 1;

      for (j = 0; j < 6; j++)
	{
	  // loop over edges

	  const Element & elemi = mesh[ei];
	  if (elemi.IsDeleted()) continue;


	  //	  (*testout) << "check element " << elemi << endl;

	  int mattyp = elemi.GetIndex();
	  
	  static const int tetedges[6][2] =
	    { { 0, 1 }, { 0, 2 }, { 0, 3 },
	      { 1, 2 }, { 1, 3 }, { 2, 3 } };

	  pi1 = elemi[tetedges[j][0]];
	  pi2 = elemi[tetedges[j][1]];
	  

	  if (pi2 < pi1) Swap (pi1, pi2);
	 	  
	  if (mesh.BoundaryEdge (pi1, pi2)) continue;


	  INDEX_2 i2 (pi1, pi2);
	  i2.Sort();
	  if (edgeused.Used(i2)) continue;
	  edgeused.Set (i2, 1);
	  
	  hasbothpoints.SetSize (0);
	  for (k = 0; k < elementsonnode[pi1].Size(); k++)
	    {
	      bool has1 = 0, has2 = 0;
	      ElementIndex elnr = elementsonnode[pi1][k];
	      const Element & elem = mesh[elnr];
	      
	      if (elem.IsDeleted()) continue;
	      
	      for (l = 0; l < elem.GetNP(); l++)
		{
		  if (elem[l] == pi1) has1 = 1;
		  if (elem[l] == pi2) has2 = 1;
		}

	      if (has1 && has2) 
		{ // only once
		  for (l = 0; l < hasbothpoints.Size(); l++)
		    if (hasbothpoints[l] == elnr)
		      has1 = 0;
		  
		  if (has1)
		    hasbothpoints.Append (elnr);
		}
	    }
	  
	  bool puretet = 1;
	  for (k = 0; k < hasbothpoints.Size(); k++)
	    if (mesh[hasbothpoints[k]].GetType () != TET)
	      puretet = 0;
	  if (!puretet)
	    continue;

	  int nsuround = hasbothpoints.Size();

	  if ( nsuround == 3 )
	    {
	      Element & elem = mesh[hasbothpoints[0]];
	      for (l = 0; l < 4; l++)
		if (elem[l] != pi1 && elem[l] != pi2)
		  {
		    pi4 = pi3;
		    pi3 = elem[l];
		  }
	      
	      el31[0] = pi1;
	      el31[1] = pi2;
	      el31[2] = pi3;
	      el31[3] = pi4;
	      el31.SetIndex (mattyp);
	      
	      if (WrongOrientation (mesh.Points(), el31))
		{
		  Swap (pi3, pi4);
		  el31[2] = pi3;
		  el31[3] = pi4;
		}
	      
	      pi5 = 0;
	      for (k = 1; k < 3; k++)
		{
		  const Element & elemk = mesh[hasbothpoints[k]];
		  bool has1 = 0;
		  for (l = 0; l < 4; l++)
		    if (elemk[l] == pi4)
		      has1 = 1;
		  if (has1)
		    {
		      for (l = 0; l < 4; l++)
			if (elemk[l] != pi1 && elemk[l] != pi2 && elemk[l] != pi4)
			  pi5 = elemk[l];
		    }
		}

	      if(pi5 == 0)
                throw NgException("Illegal state observed in SwapImprove");


	      
	      el32[0] = pi1;  
	      el32[1] = pi2;  
	      el32[2] = pi4;  
	      el32[3] = pi5;  
	      el32.SetIndex (mattyp);
	      
	      el33[0] = pi1;  
	      el33[1] = pi2;  
	      el33[2] = pi5;  
	      el33[3] = pi3;  
	      el33.SetIndex (mattyp);
	      
	      elementsonnode.Add (pi4, hasbothpoints[1]);
	      elementsonnode.Add (pi3, hasbothpoints[2]);
	      
	      bad1 = CalcBad (mesh.Points(), el31, 0) + 
		CalcBad (mesh.Points(), el32, 0) +
		CalcBad (mesh.Points(), el33, 0);

	      el31.flags.illegal_valid = 0;
	      el32.flags.illegal_valid = 0;
	      el33.flags.illegal_valid = 0;

	      if (!mesh.LegalTet(el31) ||
		  !mesh.LegalTet(el32) ||
		  !mesh.LegalTet(el33))
		bad1 += 1e4;
	      
	      el21[0] = pi3;
	      el21[1] = pi4;
	      el21[2] = pi5;
	      el21[3] = pi2;
	      el21.SetIndex (mattyp);
	      
	      el22[0] = pi5;
	      el22[1] = pi4;
	      el22[2] = pi3;
	      el22[3] = pi1;
	      el22.SetIndex (mattyp);	      

	      bad2 = CalcBad (mesh.Points(), el21, 0) + 
		CalcBad (mesh.Points(), el22, 0);

	      el21.flags.illegal_valid = 0;
	      el22.flags.illegal_valid = 0;

	      if (!mesh.LegalTet(el21) ||
		  !mesh.LegalTet(el22))
		bad2 += 1e4;


	      if (goal == OPT_CONFORM && bad2 < 1e4)
		{
		  INDEX_3 face(pi3, pi4, pi5);
		  face.Sort();
		  if (faces.Used(face))
		    {
		      // (*testout) << "3->2 swap, could improve conformity, bad1 = " << bad1
		      //				 << ", bad2 = " << bad2 << endl;
		      if (bad2 < 1e4)
			bad1 = 2 * bad2;
		    }
		  /*
		    else
		    {
		    INDEX_2 hi1(pi3, pi4);
		    hi1.Sort();
		    INDEX_2 hi2(pi3, pi5);
		    hi2.Sort();
		    INDEX_2 hi3(pi4, pi5);
		    hi3.Sort();

		    if (boundaryedges->Used (hi1) ||
		    boundaryedges->Used (hi2) ||
		    boundaryedges->Used (hi3) )
		    bad1 = 2 * bad2;
		    }
		  */
		}

	      if (bad2 < bad1)
		{
		  //		  (*mycout) << "3->2 " << flush;
		  //		  (*testout) << "3->2 conversion" << endl;
		  cnt++;
		  

		  /*
		  (*testout) << "3->2 swap, old els = " << endl
			     << mesh[hasbothpoints[0]] << endl
			     << mesh[hasbothpoints[1]] << endl
			     << mesh[hasbothpoints[2]] << endl
			     << "new els = " << endl
			     << el21 << endl
			     << el22 << endl;
		  */
                  
		  el21.flags.illegal_valid = 0;
		  el22.flags.illegal_valid = 0;
		  mesh[hasbothpoints[0]] = el21;
		  mesh[hasbothpoints[1]] = el22;
		  for (l = 0; l < 4; l++)
		    mesh[hasbothpoints[2]][l] = 0;
		  mesh[hasbothpoints[2]].Delete();

		  for (k = 0; k < 2; k++)
		    for (l = 0; l < 4; l++)
		      elementsonnode.Add (mesh[hasbothpoints[k]][l], hasbothpoints[k]);
		}
	    }
	  

	  if (nsuround == 4)
	    {
	      const Element & elem1 = mesh[hasbothpoints[0]];
	      for (l = 0; l < 4; l++)
		if (elem1[l] != pi1 && elem1[l] != pi2)
		  {
		    pi4 = pi3;
		    pi3 = elem1[l];
		  }
	      
	      el1[0] = pi1; el1[1] = pi2;
	      el1[2] = pi3; el1[3] = pi4;
	      el1.SetIndex (mattyp);

	      if (WrongOrientation (mesh.Points(), el1))
		{
		  Swap (pi3, pi4);
		  el1[2] = pi3;
		  el1[3] = pi4;
		}
	      
	      pi5 = 0;
	      for (k = 1; k < 4; k++)
		{
		  const Element & elem = mesh[hasbothpoints[k]];
		  bool has1 = 0;
		  for (l = 0; l < 4; l++)
		    if (elem[l] == pi4)
		      has1 = 1;
		  if (has1)
		    {
		      for (l = 0; l < 4; l++)
			if (elem[l] != pi1 && elem[l] != pi2 && elem[l] != pi4)
			  pi5 = elem[l];
		    }
		}
	      
	      pi6 = 0;
	      for (k = 1; k < 4; k++)
		{
		  const Element & elem = mesh[hasbothpoints[k]];
		  bool has1 = 0;
		  for (l = 0; l < 4; l++)
		    if (elem[l] == pi3)
		      has1 = 1;
		  if (has1)
		    {
		      for (l = 0; l < 4; l++)
			if (elem[l] != pi1 && elem[l] != pi2 && elem[l] != pi3)
			  pi6 = elem[l];
		    }
		}
	      
	      /*
	      INDEX_2 i22(pi3, pi5);
	      i22.Sort();
	      INDEX_2 i23(pi4, pi6);
	      i23.Sort();
	      */
	      
	      el1[0] = pi1; el1[1] = pi2;  
	      el1[2] = pi3; el1[3] = pi4; 
	      el1.SetIndex (mattyp);
	      
	      el2[0] = pi1; el2[1] = pi2;  
	      el2[2] = pi4; el2[3] = pi5;  
	      el2.SetIndex (mattyp);
	      
	      el3[0] = pi1; el3[1] = pi2;  
	      el3[2] = pi5; el3[3] = pi6;  
	      el3.SetIndex (mattyp);
	      
	      el4[0] = pi1; el4[1] = pi2;  
	      el4[2] = pi6; el4[3] = pi3;  
	      el4.SetIndex (mattyp);
	      
	      //        elementsonnode.Add (pi4, hasbothpoints.Elem(2));
	      //        elementsonnode.Add (pi3, hasbothpoints.Elem(3));
	      
	      bad1 = CalcBad (mesh.Points(), el1, 0) + 
		CalcBad (mesh.Points(), el2, 0) +
		CalcBad (mesh.Points(), el3, 0) + 
		CalcBad (mesh.Points(), el4, 0);
	      

	      el1.flags.illegal_valid = 0;
	      el2.flags.illegal_valid = 0;
	      el3.flags.illegal_valid = 0;
	      el4.flags.illegal_valid = 0;


	      if (goal != OPT_CONFORM)
		{
		  if (!mesh.LegalTet(el1) ||
		      !mesh.LegalTet(el2) ||
		      !mesh.LegalTet(el3) ||
		      !mesh.LegalTet(el4))
		    bad1 += 1e4;
		}
	      
	      el1[0] = pi3; el1[1] = pi5;  
	      el1[2] = pi2; el1[3] = pi4; 
	      el1.SetIndex (mattyp);
	 
	      el2[0] = pi3; el2[1] = pi5;  
	      el2[2] = pi4; el2[3] = pi1;  
	      el2.SetIndex (mattyp);
	      
	      el3[0] = pi3; el3[1] = pi5;  
	      el3[2] = pi1; el3[3] = pi6;  
	      el3.SetIndex (mattyp);
	      
	      el4[0] = pi3; el4[1] = pi5;  
	      el4[2] = pi6; el4[3] = pi2;  	
	      el4.SetIndex (mattyp);
      
	      bad2 = CalcBad (mesh.Points(), el1, 0) + 
		CalcBad (mesh.Points(), el2, 0) +
		CalcBad (mesh.Points(), el3, 0) + 
		CalcBad (mesh.Points(), el4, 0);

	      el1.flags.illegal_valid = 0;
	      el2.flags.illegal_valid = 0;
	      el3.flags.illegal_valid = 0;
	      el4.flags.illegal_valid = 0;

	      if (goal != OPT_CONFORM)
		{
		  if (!mesh.LegalTet(el1) ||
		      !mesh.LegalTet(el2) ||
		      !mesh.LegalTet(el3) ||
		      !mesh.LegalTet(el4))
		    bad2 += 1e4;
		}

	      
	      el1b[0] = pi4; el1b[1] = pi6;  
	      el1b[2] = pi3; el1b[3] = pi2; 
	      el1b.SetIndex (mattyp);
	      
	      el2b[0] = pi4; el2b[1] = pi6;  
	      el2b[2] = pi2; el2b[3] = pi5;  
	      el2b.SetIndex (mattyp);
	      
	      el3b[0] = pi4; el3b[1] = pi6;  
	      el3b[2] = pi5; el3b[3] = pi1;  
	      el3b.SetIndex (mattyp);
	      
	      el4b[0] = pi4; el4b[1] = pi6;  
	      el4b[2] = pi1; el4b[3] = pi3;  
	      el4b.SetIndex (mattyp);
	      
	      bad3 = CalcBad (mesh.Points(), el1b, 0) + 
		CalcBad (mesh.Points(), el2b, 0) +
		CalcBad (mesh.Points(), el3b, 0) +
		CalcBad (mesh.Points(), el4b, 0);
	      
	      el1b.flags.illegal_valid = 0;
	      el2b.flags.illegal_valid = 0;
	      el3b.flags.illegal_valid = 0;
	      el4b.flags.illegal_valid = 0;

	      if (goal != OPT_CONFORM)
		{
		  if (!mesh.LegalTet(el1b) ||
		      !mesh.LegalTet(el2b) ||
		      !mesh.LegalTet(el3b) ||
		      !mesh.LegalTet(el4b))
		    bad3 += 1e4;
		}


	      /*
	      int swap2 = (bad2 < bad1) && (bad2 < bad3);
	      int swap3 = !swap2 && (bad3 < bad1);
	      
	      if ( ((bad2 < 10 * bad1) || 
		    (bad2 < 1e6)) && mesh.BoundaryEdge (pi3, pi5))
		swap2 = 1;
	      else  if ( ((bad3 < 10 * bad1) || 
			  (bad3 < 1e6)) && mesh.BoundaryEdge (pi4, pi6))
		{
		  swap3 = 1;
		  swap2 = 0;
		}
	      */
	      bool swap2, swap3;

	      if (goal != OPT_CONFORM)
		{
		  swap2 = (bad2 < bad1) && (bad2 < bad3);
		  swap3 = !swap2 && (bad3 < bad1);
		}
	      else
		{
		  if (mesh.BoundaryEdge (pi3, pi5)) bad2 /= 1e6;
		  if (mesh.BoundaryEdge (pi4, pi6)) bad3 /= 1e6;

		  swap2 = (bad2 < bad1) && (bad2 < bad3);
		  swap3 = !swap2 && (bad3 < bad1);
		}
		

	      if (swap2 || swap3)
		{
		  // (*mycout) << "4->4 " << flush;
		  cnt++;
		  //		  (*testout) << "4->4 conversion" << "\n";
		  /*
		    (*testout) << "bad1 = " << bad1 
		    << " bad2 = " << bad2
		    << " bad3 = " << bad3 << "\n";
		  
		    (*testout) << "Points: " << pi1 << " " << pi2 << " " << pi3 
		    << " " << pi4 << " " << pi5 << " " << pi6 << "\n";
		    (*testout) << "Elements: " 
		    << hasbothpoints.Get(1) << "  "
		    << hasbothpoints.Get(2) << "  "
		    << hasbothpoints.Get(3) << "  "
		    << hasbothpoints.Get(4) << "  " << "\n";
		  */

		  /*
		    {
		    int i1, j1;
		    for (i1 = 1; i1 <= 4; i1++)
		    {
		    for (j1 = 1; j1 <= 4; j1++)
		    (*testout) << volelements.Get(hasbothpoints.Get(i1)).PNum(j1)
		    << "  ";
		    (*testout) << "\n";
		    }
		    }
		  */
		}
	      

	      if (swap2)
		{
		  //		  (*mycout) << "bad1 = " << bad1 << " bad2 = " << bad2 << "\n";


		  /*
		  (*testout) << "4->4 swap A, old els = " << endl
			     << mesh[hasbothpoints[0]] << endl
			     << mesh[hasbothpoints[1]] << endl
			     << mesh[hasbothpoints[2]] << endl
			     << mesh[hasbothpoints[3]] << endl
			     << "new els = " << endl
			     << el1 << endl
			     << el2 << endl
			     << el3 << endl
			     << el4 << endl;
		  */



		  el1.flags.illegal_valid = 0;
		  el2.flags.illegal_valid = 0;
		  el3.flags.illegal_valid = 0;
		  el4.flags.illegal_valid = 0;
		  
		  mesh[hasbothpoints[0]] = el1;
		  mesh[hasbothpoints[1]] = el2;
		  mesh[hasbothpoints[2]] = el3;
		  mesh[hasbothpoints[3]] = el4;
		  
		  for (k = 0; k < 4; k++)
		    for (l = 0; l < 4; l++)
		      elementsonnode.Add (mesh[hasbothpoints[k]][l], hasbothpoints[k]);
		}
	      else if (swap3)
		{
		  // (*mycout) << "bad1 = " << bad1 << " bad3 = " << bad3 << "\n";
		  el1b.flags.illegal_valid = 0;
		  el2b.flags.illegal_valid = 0;
		  el3b.flags.illegal_valid = 0;
		  el4b.flags.illegal_valid = 0;


		  /*
		  (*testout) << "4->4 swap A, old els = " << endl
			     << mesh[hasbothpoints[0]] << endl
			     << mesh[hasbothpoints[1]] << endl
			     << mesh[hasbothpoints[2]] << endl
			     << mesh[hasbothpoints[3]] << endl
			     << "new els = " << endl
			     << el1b << endl
			     << el2b << endl
			     << el3b << endl
			     << el4b << endl;
		  */


		  mesh[hasbothpoints[0]] = el1b;
		  mesh[hasbothpoints[1]] = el2b;
		  mesh[hasbothpoints[2]] = el3b;
		  mesh[hasbothpoints[3]] = el4b;


		  for (k = 0; k < 4; k++)
		    for (l = 0; l < 4; l++)
		      elementsonnode.Add (mesh[hasbothpoints[k]][l], hasbothpoints[k]);
		}
	    }

	  if (nsuround >= 5) 
	    {
	      Element hel(TET);

	      ArrayMem<PointIndex, 50> suroundpts(nsuround);
	      ArrayMem<char, 50> tetused(nsuround);

	      Element & elem = mesh[hasbothpoints[0]];

	      for (l = 0; l < 4; l++)
		if (elem[l] != pi1 && elem[l] != pi2)
		  {
		    pi4 = pi3;
		    pi3 = elem[l];
		  }

	      hel[0] = pi1;
	      hel[1] = pi2;
	      hel[2] = pi3;
	      hel[3] = pi4;
	      hel.SetIndex (mattyp);
	      
	      if (WrongOrientation (mesh.Points(), hel))
		{
		  Swap (pi3, pi4);
		  hel[2] = pi3;
		  hel[3] = pi4;
		}

	      
	      // suroundpts.SetSize (nsuround);
	      suroundpts[0] = pi3;
	      suroundpts[1] = pi4;

	      tetused = 0;
	      tetused[0] = 1;

	      for (l = 2; l < nsuround; l++)
		{
		  int oldpi = suroundpts[l-1];
		  int newpi = 0;

		  for (k = 0; k < nsuround && !newpi; k++)
		    if (!tetused[k])
		      {
			const Element & nel = mesh[hasbothpoints[k]];

			for (int k2 = 0; k2 < 4 && !newpi; k2++)
			  if (nel[k2] == oldpi)
			    {
			      newpi = 
				nel[0] + nel[1] + nel[2] + nel[3] 
				- pi1 - pi2 - oldpi;
			      
			      tetused[k] = 1; 
			      suroundpts[l] = newpi;
			    }			
		      }
		}

	      
	      bad1 = 0;
	      for (k = 0; k < nsuround; k++)
		{
		  hel[0] = pi1;
		  hel[1] = pi2;
		  hel[2] = suroundpts[k];
		  hel[3] = suroundpts[(k+1) % nsuround];
		  hel.SetIndex (mattyp);

		  bad1 += CalcBad (mesh.Points(), hel, 0);
		}

	      //  (*testout) << "nsuround = " << nsuround << " bad1 = " << bad1 << endl;


	      int bestl = -1;
	      int confface = -1;
	      int confedge = -1;
	      double badopt = bad1;

	      for (l = 0; l < nsuround; l++)
		{
		  bad2 = 0;

		  for (k = l+1; k <= nsuround + l - 2; k++)
		    {
		      hel[0] = suroundpts[l];
		      hel[1] = suroundpts[k % nsuround];
		      hel[2] = suroundpts[(k+1) % nsuround];
		      hel[3] = pi2;

		      bad2 += CalcBad (mesh.Points(), hel, 0);
		      hel.flags.illegal_valid = 0;
		      if (!mesh.LegalTet(hel)) bad2 += 1e4;

		      hel[2] = suroundpts[k % nsuround];
		      hel[1] = suroundpts[(k+1) % nsuround];
		      hel[3] = pi1;

		      bad2 += CalcBad (mesh.Points(), hel, 0);

		      hel.flags.illegal_valid = 0;
		      if (!mesh.LegalTet(hel)) bad2 += 1e4;
		    }
		  // (*testout) << "bad2," << l << " = " << bad2 << endl;
		  
		  if ( bad2 < badopt )
		    {
		      bestl = l;
		      badopt = bad2;
		    }
		  
		  
		  if (goal == OPT_CONFORM)
		       // (bad2 <= 100 * bad1 || bad2 <= 1e6))
		    {
		      bool nottoobad =
			(bad2 <= bad1) ||
			(bad2 <= 100 * bad1 && bad2 <= 1e18) ||
			(bad2 <= 1e8);
		      
		      for (k = l+1; k <= nsuround + l - 2; k++)
			{
			  INDEX_3 hi3(suroundpts[l],
				      suroundpts[k % nsuround],
				      suroundpts[(k+1) % nsuround]);
			  hi3.Sort();
			  if (faces.Used(hi3))
			    {
			      // (*testout) << "could improve face conformity, bad1 = " << bad1
			      // << ", bad 2 = " << bad2 << ", nottoobad = " << nottoobad << endl;
			      if (nottoobad)
				confface = l;
			    }
			}

		      for (k = l+2; k <= nsuround+l-2; k++)
			{
			  if (mesh.BoundaryEdge (suroundpts[l],
						 suroundpts[k % nsuround]))
			    {
			      /*
			      *testout << "could improve edge conformity, bad1 = " << bad1
				   << ", bad 2 = " << bad2 << ", nottoobad = " << nottoobad << endl;
			      */
			      if (nottoobad)
				confedge = l;
			    }
			}
		    }
		}
	      
	      if (confedge != -1)
		bestl = confedge;
	      if (confface != -1)
		bestl = confface;

	      if (bestl != -1)
		{
		  // (*mycout) << nsuround << "->" << 2 * (nsuround-2) << " " << flush;
		  cnt++;
		  
		  for (k = bestl+1; k <= nsuround + bestl - 2; k++)
		    {
		      int k1;

		      hel[0] = suroundpts[bestl];
		      hel[1] = suroundpts[k % nsuround];
		      hel[2] = suroundpts[(k+1) % nsuround];
		      hel[3] = pi2;
		      hel.flags.illegal_valid = 0;

		      /*
		      (*testout) << nsuround << "-swap, new el,top = "
				 << hel << endl;
		      */
		      mesh.AddVolumeElement (hel);

		      for (k1 = 0; k1 < 4; k1++)
			elementsonnode.Add (hel[k1], mesh.GetNE()-1);
		      

		      hel[2] = suroundpts[k % nsuround];
		      hel[1] = suroundpts[(k+1) % nsuround];
		      hel[3] = pi1;

		      /*
		      (*testout) << nsuround << "-swap, new el,bot = "
				 << hel << endl;
		      */

		      mesh.AddVolumeElement (hel);

		      for (k1 = 0; k1 < 4; k1++)
			elementsonnode.Add (hel[k1], mesh.GetNE()-1);
		    }		  
		  
		  for (k = 0; k < nsuround; k++)
		    {
		      Element & rel = mesh[hasbothpoints[k]];
		      /*
		      (*testout) << nsuround << "-swap, old el = "
				 << rel << endl;
		      */
		      rel.Delete();
		      for (int k1 = 0; k1 < 4; k1++)
			rel[k1] = 0;

		    }
		}
	    }
	}

      /*
	if (onlybedges)
	{
	(*testout) << "bad tet: " 
	<< volelements.Get(i)[0] 
	<< volelements.Get(i)[1] 
	<< volelements.Get(i)[2] 
	<< volelements.Get(i)[3] << "\n";

	if (!mesh.LegalTet (volelements.Get(i)))
	cerr << "Illegal tet" << "\n";
	}
      */
    }
  //  (*mycout) << endl;

  /*  
      cout << "edgeused: ";
      edgeused.PrintMemInfo(cout);
  */
  PrintMessage (5, cnt, " swaps performed");





  mesh.Compress ();

  /*
  if (goal == OPT_QUALITY)
    {
      bad1 = CalcTotalBad (mesh.Points(), mesh.VolumeElements());
      //      (*testout) << "Total badness = " << bad1 << endl;
    }
  */

  /*
    for (i = 1; i <= GetNE(); i++)
    if (volelements.Get(i)[0])
    if (!mesh.LegalTet (volelements.Get(i)))
    {
    cout << "detected illegal tet, 2" << endl;
    (*testout) << "detected illegal tet1: " << i << endl;
    }
  */

  multithread.task = savetask;
}
  





void MeshOptimize3d :: SwapImproveSurface (Mesh & mesh, OPTIMIZEGOAL goal,
					   const BitArray * working_elements,
					   const ARRAY< ARRAY<int,PointIndex::BASE>* > * idmaps)
{
  ARRAY< ARRAY<int,PointIndex::BASE>* > locidmaps;
  const ARRAY< ARRAY<int,PointIndex::BASE>* > * used_idmaps;

  if(idmaps)
    used_idmaps = idmaps;
  else
    {
      used_idmaps = &locidmaps;
      
      for(int i=1; i<=mesh.GetIdentifications().GetMaxNr(); i++)
	{
	  if(mesh.GetIdentifications().GetType(i) == Identifications::PERIODIC)
	    {
	      locidmaps.Append(new ARRAY<int,PointIndex::BASE>);
	      mesh.GetIdentifications().GetMap(i,*locidmaps.Last(),true);
	    }
	}
    }

  ElementIndex ei;
  SurfaceElementIndex sei;

  PointIndex pi1, pi2, pi3, pi4, pi5, pi6;
  PointIndex pi1other, pi2other;
  int cnt = 0;

  //double bad1, bad2, bad3, sbad;
  double bad1, sbad;
  double h;

  int np = mesh.GetNP();
  int ne = mesh.GetNE();
  int nse = mesh.GetNSE();

  int mattype, othermattype;

  
  // contains at least all elements at node
  TABLE<ElementIndex,PointIndex::BASE> elementsonnode(np);
  TABLE<SurfaceElementIndex,PointIndex::BASE> surfaceelementsonnode(np);
  TABLE<int,PointIndex::BASE> surfaceindicesonnode(np);

  ARRAY<ElementIndex> hasbothpoints;
  ARRAY<ElementIndex> hasbothpointsother;

  PrintMessage (3, "SwapImproveSurface ");
  (*testout) << "\n" << "Start SwapImproveSurface" << endl;

  const char * savetask = multithread.task;
  multithread.task = "Swap Improve Surface";
    
      
  
  // find elements on node
  for (ei = 0; ei < ne; ei++)
    for (int j = 0; j < mesh[ei].GetNP(); j++)
      elementsonnode.Add (mesh[ei][j], ei);

  for (sei = 0; sei < nse; sei++)
    for(int j=0; j<mesh[sei].GetNP(); j++)
      {
	surfaceelementsonnode.Add(mesh[sei][j], sei);
	if(!surfaceindicesonnode[mesh[sei][j]].Contains(mesh[sei].GetIndex()))
	  surfaceindicesonnode.Add(mesh[sei][j],mesh[sei].GetIndex());
      }

  bool periodic;
  int idnum(-1);

  INDEX_2_HASHTABLE<int> edgeused(2 * ne + 5);

  for (ei = 0; ei < ne; ei++)
    {
      if (multithread.terminate)
	break;
      
      multithread.percent = 100.0 * (ei+1) / ne;

      if (mesh.ElementType(ei) == FIXEDELEMENT)
	continue;
      
      if(working_elements && 
	 ei < working_elements->Size() &&
	 !working_elements->Test(ei))
	continue;

      if (mesh[ei].IsDeleted())
	continue;

      if ((goal == OPT_LEGAL) && 
	  mesh.LegalTet (mesh[ei]) &&
	  CalcBad (mesh.Points(), mesh[ei], 0) < 1e3)
	continue;

      const Element & elemi = mesh[ei];
      //Element elemi = mesh[ei];
      if (elemi.IsDeleted()) continue;


      mattype = elemi.GetIndex();

      bool swapped = false;

      for (int j = 0; !swapped && j < 6; j++)
	{
	  // loop over edges

	  
	  static const int tetedges[6][2] =
	    { { 0, 1 }, { 0, 2 }, { 0, 3 },
	      { 1, 2 }, { 1, 3 }, { 2, 3 } };

	  pi1 = elemi[tetedges[j][0]];
	  pi2 = elemi[tetedges[j][1]];

	  
	  if (pi2 < pi1)
	    Swap (pi1, pi2);
	    	  
	  
	  bool found = false;
	  for(int k=0; !found && k<used_idmaps->Size(); k++)
	    {
	      if(pi2 < (*used_idmaps)[k]->Size() + PointIndex::BASE)
		{
		  pi1other = (*(*used_idmaps)[k])[pi1];
		  pi2other = (*(*used_idmaps)[k])[pi2];
		  found = (pi1other != 0 && pi2other != 0 && pi1other != pi1 && pi2other != pi2);
		  if(found)
		    idnum = k;
		}
	    }
	  if(found)
	    periodic = true;
	  else
	    {
	      periodic = false;
	      pi1other = pi1; pi2other = pi2;
	    }


	 	  
	  if (!mesh.BoundaryEdge (pi1, pi2) ||
	      mesh.IsSegment(pi1, pi2)) continue;

	  othermattype = -1;

	  
	  INDEX_2 i2 (pi1, pi2);
	  i2.Sort();
	  if (edgeused.Used(i2)) continue;
	  edgeused.Set (i2, 1);
	  if(periodic)
	    {
	      i2.I1() = pi1other;
	      i2.I2() = pi2other;
	      i2.Sort();
	      edgeused.Set(i2,1);
	    }
	  
	  
	  hasbothpoints.SetSize (0);
	  hasbothpointsother.SetSize (0);
	  for (int k = 0; k < elementsonnode[pi1].Size(); k++)
	    {
	      bool has1 = false, has2 = false;
	      ElementIndex elnr = elementsonnode[pi1][k];
	      const Element & elem = mesh[elnr];
	      
	      if (elem.IsDeleted()) continue;
	      
	      for (int l = 0; l < elem.GetNP(); l++)
		{
		  if (elem[l] == pi1) has1 = true;
		  if (elem[l] == pi2) has2 = true;
		}

	      if (has1 && has2) 
		{ 
		  if(othermattype == -1 && elem.GetIndex() != mattype)
		    othermattype = elem.GetIndex();

		  if(elem.GetIndex() == mattype)
		    {
		      // only once
		      for (int l = 0; l < hasbothpoints.Size(); l++)
			if (hasbothpoints[l] == elnr)
			  has1 = 0;
		      
		      if (has1)
			hasbothpoints.Append (elnr);
		    }
		  else if(elem.GetIndex() == othermattype)
		    {
		      // only once
		      for (int l = 0; l < hasbothpointsother.Size(); l++)
			if (hasbothpointsother[l] == elnr)
			  has1 = 0;
		      
		      if (has1)
			hasbothpointsother.Append (elnr);
		    }
		  else
		    {
		      cout << "problem with domain indices" << endl;
		      (*testout) << "problem: mattype = " << mattype << ", othermattype = " << othermattype 
				 << " elem " << elem << " mt " << elem.GetIndex() << endl
				 << " pi1 " << pi1 << " pi2 " << pi2 << endl;
		      (*testout) << "hasbothpoints:" << endl;
		      for(int ii=0; ii < hasbothpoints.Size(); ii++)
			(*testout) << mesh[hasbothpoints[ii]] << endl;
		      (*testout) << "hasbothpointsother:" << endl;
		      for(int ii=0; ii < hasbothpointsother.Size(); ii++)
			(*testout) << mesh[hasbothpointsother[ii]] << endl;
		    }
		}
	    }

	  if(hasbothpointsother.Size() > 0 && periodic)
	    throw NgException("SwapImproveSurface: Assumption about interface/periodicity wrong!");

	  if(periodic)
	    {
	      for (int k = 0; k < elementsonnode[pi1other].Size(); k++)
		{
		  bool has1 = false, has2 = false;
		  ElementIndex elnr = elementsonnode[pi1other][k];
		  const Element & elem = mesh[elnr];
	      
		  if (elem.IsDeleted()) continue;
	      
		  for (int l = 0; l < elem.GetNP(); l++)
		    {
		      if (elem[l] == pi1other) has1 = true;
		      if (elem[l] == pi2other) has2 = true;
		    }
		  
		  if (has1 && has2) 
		    { 
		      if(othermattype == -1)
			othermattype = elem.GetIndex();

		      // only once
		      for (int l = 0; l < hasbothpointsother.Size(); l++)
			if (hasbothpointsother[l] == elnr)
			  has1 = 0;
		      
		      if (has1)
			hasbothpointsother.Append (elnr);
		    }
		}
	    }


	  //for(k=0; k<hasbothpoints.Size(); k++)
	  //  (*testout) << "hasbothpoints["<<k<<"]: " << mesh[hasbothpoints[k]] << endl;

	  
	  SurfaceElementIndex sel1=-1,sel2=-1;
	  SurfaceElementIndex sel1other=-1,sel2other=-1;
	  for(int k = 0; k < surfaceelementsonnode[pi1].Size(); k++)
	    {
	      bool has1 = false, has2 = false;
	      SurfaceElementIndex elnr = surfaceelementsonnode[pi1][k];
	      const Element2d & elem = mesh[elnr];

	      if (elem.IsDeleted()) continue;

	      for (int l = 0; l < elem.GetNP(); l++)
		{
		  if (elem[l] == pi1) has1 = true;
		  if (elem[l] == pi2) has2 = true;
		}

	      if(has1 && has2 && elnr != sel2)
		{
		  sel1 = sel2;
		  sel2 = elnr;
		}
	    }

	  if(periodic)
	    {
	      for(int k = 0; k < surfaceelementsonnode[pi1other].Size(); k++)
		{
		  bool has1 = false, has2 = false;
		  SurfaceElementIndex elnr = surfaceelementsonnode[pi1other][k];
		  const Element2d & elem = mesh[elnr];

		  if (elem.IsDeleted()) continue;

		  for (int l = 0; l < elem.GetNP(); l++)
		    {
		      if (elem[l] == pi1other) has1 = true;
		      if (elem[l] == pi2other) has2 = true;
		    }

		  if(has1 && has2 && elnr != sel2other)
		    {
		      sel1other = sel2other;
		      sel2other = elnr;
		    }
		}
	    }
	  else
	    {
	      sel1other = sel1; sel2other = sel2;
	    }

	  //(*testout) << "sel1 " << sel1 << " sel2 " << sel2 << " el " << mesh[sel1] << " resp. " << mesh[sel2] << endl;

	  PointIndex sp1(0), sp2(0);
	  PointIndex sp1other, sp2other;
	  for(int l=0; l<mesh[sel1].GetNP(); l++)
	    if(mesh[sel1][l] != pi1 && mesh[sel1][l] != pi2)
	      sp1 = mesh[sel1][l];
	  for(int l=0; l<mesh[sel2].GetNP(); l++)
	    if(mesh[sel2][l] != pi1 && mesh[sel2][l] != pi2)
	      sp2 = mesh[sel2][l];

	  if(periodic)
	    {
	      sp1other = (*(*used_idmaps)[idnum])[sp1];
	      sp2other = (*(*used_idmaps)[idnum])[sp2];

	      bool change = false;
	      for(int l=0; !change && l<mesh[sel1other].GetNP(); l++)
		change = (sp2other == mesh[sel1other][l]);
	      
	      if(change)
		{
		  SurfaceElementIndex aux = sel1other;
		  sel1other = sel2other;
		  sel2other = aux;
		}

	    }
	  else
	    {
	      sp1other = sp1; sp2other = sp2;
	    }
	  
	  Vec<3> v1 = mesh[sp1]-mesh[pi1],
	    v2 = mesh[sp2]-mesh[pi1],
	    v3 = mesh[sp1]-mesh[pi2],
	    v4 = mesh[sp2]-mesh[pi2];
	  double vol = 0.5*(Cross(v1,v2).Length() + Cross(v3,v4).Length());
	  h = sqrt(vol);
	  h = 0;

	  sbad = CalcTriangleBadness (mesh[pi1],mesh[pi2],mesh[sp1],0,0) + 
	    CalcTriangleBadness (mesh[pi2],mesh[pi1],mesh[sp2],0,0);
	  


	  bool puretet = true;
	  for (int k = 0; puretet && k < hasbothpoints.Size(); k++)
	    if (mesh[hasbothpoints[k]].GetType () != TET)
	      puretet = false;
	  for (int k = 0; puretet && k < hasbothpointsother.Size(); k++)
	    if (mesh[hasbothpointsother[k]].GetType () != TET)
	      puretet = false;
	  if (!puretet)
	    continue;

	  int nsuround = hasbothpoints.Size();
	  int nsuroundother = hasbothpointsother.Size();

	  ARRAY < int > outerpoints(nsuround+1);
	  outerpoints[0] = sp1;

	  for(int i=0; i<nsuround; i++)
	    {
	      bool done = false;
	      for(int jj=i; !done && jj<hasbothpoints.Size(); jj++)
		{
		  for(int k=0; !done && k<4; k++)
		    if(mesh[hasbothpoints[jj]][k] == outerpoints[i])
		      {
			done = true;
			for(int l=0; l<4; l++)
			  if(mesh[hasbothpoints[jj]][l] != pi1 &&
			     mesh[hasbothpoints[jj]][l] != pi2 &&
			     mesh[hasbothpoints[jj]][l] != outerpoints[i])
			    outerpoints[i+1] = mesh[hasbothpoints[jj]][l];
		      }
		  if(done)
		    {
		      ElementIndex aux = hasbothpoints[i];
		      hasbothpoints[i] = hasbothpoints[jj];
		      hasbothpoints[jj] = aux;
		    }
		}
	    }
	  if(outerpoints[nsuround] != sp2)
	    {
	      cerr << "OJE OJE OJE" << endl;
	      (*testout) << "OJE OJE OJE" << endl;
	      (*testout) << "hasbothpoints: " << endl;
	      for(int ii=0; ii < hasbothpoints.Size(); ii++)
		{
		  (*testout) << mesh[hasbothpoints[ii]] << endl;
		  for(int jj=0; jj<mesh[hasbothpoints[ii]].GetNP(); jj++)
		    if(mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][0] > 0)
		      (*testout) << mesh[hasbothpoints[ii]][jj] << " between "
				 << mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][0] << " and "
				 << mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][1] << endl;
		}
	      (*testout) << "outerpoints: " << outerpoints << endl;
	      (*testout) << "sel1 " << mesh[sel1] << endl
			 << "sel2 " << mesh[sel2] << endl;
	      for(int ii=0; ii<3; ii++)
		{
		  if(mesh.mlbetweennodes[mesh[sel1][ii]][0] > 0)
		    (*testout) << mesh[sel1][ii] << " between "
			       << mesh.mlbetweennodes[mesh[sel1][ii]][0] << " and "
			       << mesh.mlbetweennodes[mesh[sel1][ii]][1] << endl;
		  if(mesh.mlbetweennodes[mesh[sel2][ii]][0] > 0)
		    (*testout) << mesh[sel2][ii] << " between "
			       << mesh.mlbetweennodes[mesh[sel2][ii]][0] << " and "
			       << mesh.mlbetweennodes[mesh[sel2][ii]][1] << endl;
		}
	    }

	  
	  ARRAY < int > outerpointsother;

	  if(nsuroundother > 0)
	    {
	      outerpointsother.SetSize(nsuroundother+1);
	      outerpointsother[0] = sp2other;
	    }

	  for(int i=0; i<nsuroundother; i++)
	    {
	      bool done = false;
	      for(int jj=i; !done && jj<hasbothpointsother.Size(); jj++)
		{
		  for(int k=0; !done && k<4; k++)
		    if(mesh[hasbothpointsother[jj]][k] == outerpointsother[i])
		      {
			done = true;
			for(int l=0; l<4; l++)
			  if(mesh[hasbothpointsother[jj]][l] != pi1other &&
			     mesh[hasbothpointsother[jj]][l] != pi2other &&
			     mesh[hasbothpointsother[jj]][l] != outerpointsother[i])
			    outerpointsother[i+1] = mesh[hasbothpointsother[jj]][l];
		      }
		  if(done)
		    {
		      ElementIndex aux = hasbothpointsother[i];
		      hasbothpointsother[i] = hasbothpointsother[jj];
		      hasbothpointsother[jj] = aux;
		    }
		}
	    }
	  if(nsuroundother > 0 && outerpointsother[nsuroundother] != sp1other)
	    {
	      cerr << "OJE OJE OJE (other)" << endl;
	      (*testout) << "OJE OJE OJE (other)" << endl;
	      (*testout) << "pi1 " << pi1 << " pi2 " << pi2 << " sp1 " << sp1 << " sp2 " << sp2 << endl;
	      (*testout) << "hasbothpoints: " << endl;
	      for(int ii=0; ii < hasbothpoints.Size(); ii++)
		{
		  (*testout) << mesh[hasbothpoints[ii]] << endl;
		  for(int jj=0; jj<mesh[hasbothpoints[ii]].GetNP(); jj++)
		    if(mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][0] > 0)
		      (*testout) << mesh[hasbothpoints[ii]][jj] << " between "
				 << mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][0] << " and "
				 << mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][1] << endl;
		}
	      (*testout) << "outerpoints: " << outerpoints << endl;
	      (*testout) << "sel1 " << mesh[sel1] << endl
			 << "sel2 " << mesh[sel2] << endl;
	      for(int ii=0; ii<3; ii++)
		{
		  if(mesh.mlbetweennodes[mesh[sel1][ii]][0] > 0)
		    (*testout) << mesh[sel1][ii] << " between "
			       << mesh.mlbetweennodes[mesh[sel1][ii]][0] << " and "
			       << mesh.mlbetweennodes[mesh[sel1][ii]][1] << endl;
		  if(mesh.mlbetweennodes[mesh[sel2][ii]][0] > 0)
		    (*testout) << mesh[sel2][ii] << " between "
			       << mesh.mlbetweennodes[mesh[sel2][ii]][0] << " and "
			       << mesh.mlbetweennodes[mesh[sel2][ii]][1] << endl;
		}
		  
	      (*testout) << "pi1other " << pi1other << " pi2other " << pi2other << " sp1other " << sp1other << " sp2other " << sp2other << endl;
	      (*testout) << "hasbothpointsother: " << endl;
	      for(int ii=0; ii < hasbothpointsother.Size(); ii++)
		{
		  (*testout) << mesh[hasbothpointsother[ii]] << endl;
		  for(int jj=0; jj<mesh[hasbothpointsother[ii]].GetNP(); jj++)
		    if(mesh.mlbetweennodes[mesh[hasbothpointsother[ii]][jj]][0] > 0)
		      (*testout) << mesh[hasbothpointsother[ii]][jj] << " between "
				 << mesh.mlbetweennodes[mesh[hasbothpointsother[ii]][jj]][0] << " and "
				 << mesh.mlbetweennodes[mesh[hasbothpointsother[ii]][jj]][1] << endl;
		}
	      (*testout) << "outerpoints: " << outerpointsother << endl;
	      (*testout) << "sel1other " << mesh[sel1other] << endl
			 << "sel2other " << mesh[sel2other] << endl;
	      for(int ii=0; ii<3; ii++)
		{
		  if(mesh.mlbetweennodes[mesh[sel1other][ii]][0] > 0)
		    (*testout) << mesh[sel1other][ii] << " between "
			       << mesh.mlbetweennodes[mesh[sel1other][ii]][0] << " and "
			       << mesh.mlbetweennodes[mesh[sel1other][ii]][1] << endl;
		  if(mesh.mlbetweennodes[mesh[sel2other][ii]][0] > 0)
		    (*testout) << mesh[sel2other][ii] << " between "
			       << mesh.mlbetweennodes[mesh[sel2other][ii]][0] << " and "
			       << mesh.mlbetweennodes[mesh[sel2other][ii]][1] << endl;
		}
	    }

	  bad1=0;
	  for(int i=0; i<hasbothpoints.Size(); i++)
	    bad1 += CalcBad(mesh.Points(), mesh[hasbothpoints[i]],h);
	  for(int i=0; i<hasbothpointsother.Size(); i++)
	    bad1 += CalcBad(mesh.Points(), mesh[hasbothpointsother[i]],h);
	  bad1 /= double(hasbothpoints.Size() + hasbothpointsother.Size());

	  
	  int startpoints,startpointsother;


	  if(outerpoints.Size() == 3)
	    startpoints = 1;
	  else if(outerpoints.Size() == 4)
	    startpoints = 2;
	  else
	    startpoints = outerpoints.Size();
	  
	  if(outerpointsother.Size() == 3)
	    startpointsother = 1;
	  else if(outerpointsother.Size() == 4)
	    startpointsother = 2;
	  else
	    startpointsother = outerpointsother.Size();
	  

	  ARRAY < ARRAY < Element* > * > newelts(startpoints);
	  ARRAY < ARRAY < Element* > * > neweltsother(startpointsother);

	  double minbad = 1e50, minbadother = 1e50, currbad;
	  int minpos = -1, minposother = -1;

	  //(*testout) << "pi1 " << pi1 << " pi2 " << pi2 << " outerpoints " << outerpoints << endl;

	  for(int i=0; i<startpoints; i++)
	    {
	      newelts[i] = new ARRAY <Element*>(2*(nsuround-1));
	      
	      for(int jj=0; jj<nsuround-1; jj++)
		{
		  (*newelts[i])[2*jj] = new Element(TET);
		  (*newelts[i])[2*jj+1] = new Element(TET);
		  Element & newel1 = *((*newelts[i])[2*jj]);
		  Element & newel2 = *((*newelts[i])[2*jj+1]);

		  newel1[0] = pi1;
		  newel1[1] = outerpoints[i];
		  newel1[2] = outerpoints[(i+jj+1)%outerpoints.Size()];
		  newel1[3] = outerpoints[(i+jj+2)%outerpoints.Size()];

		  newel2[0] = pi2;
		  newel2[1] = outerpoints[i];
		  newel2[2] = outerpoints[(i+jj+2)%outerpoints.Size()];
		  newel2[3] = outerpoints[(i+jj+1)%outerpoints.Size()];
		  

		  //(*testout) << "j " << j << " newel1 " << newel1[0] << " "<< newel1[1] << " "<< newel1[2] << " "<< newel1[3] << endl
		  //     << " newel2 " << newel2[0] << " "<< newel2[1] << " "<< newel2[2] << " "<< newel2[3] << endl;
		  
		  newel1.SetIndex(mattype);
		  newel2.SetIndex(mattype);

		}

	      bool wrongorientation = true;
	      for(int jj = 0; wrongorientation && jj<newelts[i]->Size(); jj++)
		wrongorientation = wrongorientation && WrongOrientation(mesh.Points(), *(*newelts[i])[jj]);
	      
	      currbad = 0;

	      for(int jj=0; jj<newelts[i]->Size(); jj++)
		{
		  if(wrongorientation)
		    Swap((*(*newelts[i])[jj])[2],(*(*newelts[i])[jj])[3]);


		  // not two new faces on same surface
		  ARRAY<int> face_index;
		  for(int k = 0; k<surfaceindicesonnode[(*(*newelts[i])[jj])[0]].Size(); k++)
		    face_index.Append(surfaceindicesonnode[(*(*newelts[i])[jj])[0]][k]);

		  for(int k=1; k<4; k++)
		    {
		      for(int l=0; l<face_index.Size(); l++)
			{
			  if(face_index[l] != -1 && 
			     !(surfaceindicesonnode[(*(*newelts[i])[jj])[k]].Contains(face_index[l])))
			    face_index[l] = -1;
			}

		    }
		      
		  for(int k=0; k<face_index.Size(); k++)
		    if(face_index[k] != -1)
		      currbad += 1e12;


		  currbad += CalcBad(mesh.Points(),*(*newelts[i])[jj],h);


		}  

	      //currbad /= double(newelts[i]->Size());
		    


	      if(currbad < minbad)
		{
		  minbad = currbad;
		  minpos = i;
		}

	    }

	  if(startpointsother == 0)
	    minbadother = 0;

	  for(int i=0; i<startpointsother; i++)
	    {
	      neweltsother[i] = new ARRAY <Element*>(2*(nsuroundother));
	      
	      for(int jj=0; jj<nsuroundother; jj++)
		{
		  (*neweltsother[i])[2*jj] = new Element(TET);
		  (*neweltsother[i])[2*jj+1] = new Element(TET);
		  Element & newel1 = *((*neweltsother[i])[2*jj]);
		  Element & newel2 = *((*neweltsother[i])[2*jj+1]);

		  newel1[0] = pi1other;
		  newel1[1] = outerpointsother[i];
		  newel1[2] = outerpointsother[(i+jj+1)%outerpointsother.Size()];
		  newel1[3] = outerpointsother[(i+jj+2)%outerpointsother.Size()];

		  newel2[0] = pi2other;
		  newel2[1] = outerpointsother[i];
		  newel2[2] = outerpointsother[(i+jj+2)%outerpointsother.Size()];
		  newel2[3] = outerpointsother[(i+jj+1)%outerpointsother.Size()];
		  

		  //(*testout) << "j " << j << " newel1 " << newel1[0] << " "<< newel1[1] << " "<< newel1[2] << " "<< newel1[3] << endl
		  //	     << " newel2 " << newel2[0] << " "<< newel2[1] << " "<< newel2[2] << " "<< newel2[3] << endl;
		  
		  newel1.SetIndex(othermattype);
		  newel2.SetIndex(othermattype);

		}

	      bool wrongorientation = true;
	      for(int jj = 0; wrongorientation && jj<neweltsother[i]->Size(); jj++)
		wrongorientation = wrongorientation && WrongOrientation(mesh.Points(), *(*neweltsother[i])[jj]);
	      
	      currbad = 0;

	      for(int jj=0; jj<neweltsother[i]->Size(); jj++)
		{
		  if(wrongorientation)
		    Swap((*(*neweltsother[i])[jj])[2],(*(*neweltsother[i])[jj])[3]);

		  currbad += CalcBad(mesh.Points(),*(*neweltsother[i])[jj],h);
		}  

	      //currbad /= double(neweltsother[i]->Size());
		    


	      if(currbad < minbadother)
		{
		  minbadother = currbad;
		  minposother = i;
		}

	    }

	  //(*testout) << "minbad " << minbad << " bad1 " << bad1 << endl;

	  
	  double sbadnew = CalcTriangleBadness (mesh[pi1],mesh[sp2],mesh[sp1],0,0) + 
	    CalcTriangleBadness (mesh[pi2],mesh[sp1],mesh[sp2],0,0);
	  

	  int denom = newelts[minpos]->Size();
	  if(minposother >= 0)
	    denom += neweltsother[minposother]->Size();
	  

	  if((minbad+minbadother)/double(denom) < bad1 && 
	     sbadnew < sbad)
	    {
	      cnt++;

	      swapped = true;


	      int start1 = -1;
	      for(int l=0; l<3; l++)
		if(mesh[sel1][l] == pi1)
		  start1 = l;
	      if(mesh[sel1][(start1+1)%3] == pi2)
		{
		  mesh[sel1][0] = pi1;
		  mesh[sel1][1] = sp2;
		  mesh[sel1][2] = sp1;
		  mesh[sel2][0] = pi2;
		  mesh[sel2][1] = sp1;
		  mesh[sel2][2] = sp2;
		}
	      else
		{
		  mesh[sel1][0] = pi2;
		  mesh[sel1][1] = sp2;
		  mesh[sel1][2] = sp1;
		  mesh[sel2][0] = pi1;
		  mesh[sel2][1] = sp1;
		  mesh[sel2][2] = sp2;
		}
	      //(*testout) << "changed surface element " << sel1 << " to " << mesh[sel1] << ", " << sel2 << " to " << mesh[sel2] << endl;

	      for(int l=0; l<3; l++)
		{
		  surfaceelementsonnode.Add(mesh[sel1][l],sel1);
		  surfaceelementsonnode.Add(mesh[sel2][l],sel2);
		}
	      


	      if(periodic)
		{
		  start1 = -1;
		  for(int l=0; l<3; l++)
		    if(mesh[sel1other][l] == pi1other)
		      start1 = l;
		  


		  //(*testout) << "changed surface elements " << mesh[sel1other] << " and " << mesh[sel2other] << endl;
		  if(mesh[sel1other][(start1+1)%3] == pi2other)
		    {
		      mesh[sel1other][0] = pi1other;
		      mesh[sel1other][1] = sp2other;
		      mesh[sel1other][2] = sp1other;
		      mesh[sel2other][0] = pi2other;
		      mesh[sel2other][1] = sp1other;
		      mesh[sel2other][2] = sp2other;
		      //(*testout) << "       with rule 1" << endl;
		    }
		  else
		    {
		      mesh[sel1other][0] = pi2other;
		      mesh[sel1other][1] = sp2other;
		      mesh[sel1other][2] = sp1other;
		      mesh[sel2other][0] = pi1other;
		      mesh[sel2other][1] = sp1other;
		      mesh[sel2other][2] = sp2other;
		      //(*testout) << "       with rule 2" << endl;
		    }
		  //(*testout) << "         to " << mesh[sel1other] << " and " << mesh[sel2other] << endl;
		  
		  //(*testout) << "  and surface element " << sel1other << " to " << mesh[sel1other] << ", " << sel2other << " to " << mesh[sel2other] << endl;

		  for(int l=0; l<3; l++)
		    {
		      surfaceelementsonnode.Add(mesh[sel1other][l],sel1other);
		      surfaceelementsonnode.Add(mesh[sel2other][l],sel2other);
		    }
		}




	      for(int i=0; i<hasbothpoints.Size(); i++)
		{
		  mesh[hasbothpoints[i]] = *(*newelts[minpos])[i];

		  for(int l=0; l<4; l++)
		    elementsonnode.Add((*(*newelts[minpos])[i])[l],hasbothpoints[i]);
		}

	      for(int i=hasbothpoints.Size(); i<(*newelts[minpos]).Size(); i++)
		{
		  ElementIndex ni = mesh.AddVolumeElement(*(*newelts[minpos])[i]);
		  
		  for(int l=0; l<4; l++)
		    elementsonnode.Add((*(*newelts[minpos])[i])[l],ni);
		}

	      if(hasbothpointsother.Size() > 0)
		{
		  for(int i=0; i<hasbothpointsother.Size(); i++)
		    {
		      mesh[hasbothpointsother[i]] = *(*neweltsother[minposother])[i];
		      for(int l=0; l<4; l++)
			elementsonnode.Add((*(*neweltsother[minposother])[i])[l],hasbothpointsother[i]);
		    }
		  
		  for(int i=hasbothpointsother.Size(); i<(*neweltsother[minposother]).Size(); i++)
		    {
		      ElementIndex ni = mesh.AddVolumeElement(*(*neweltsother[minposother])[i]);
		      for(int l=0; l<4; l++)
			elementsonnode.Add((*(*neweltsother[minposother])[i])[l],ni);
		    }
		}

	      

	    }

	  for(int i=0; i<newelts.Size(); i++)
	    {
	      for(int jj=0; jj<newelts[i]->Size(); jj++)
		delete (*newelts[i])[jj];
	      delete newelts[i];
	    }

	  for(int i=0; i<neweltsother.Size(); i++)
	    {
	      for(int jj=0; jj<neweltsother[i]->Size(); jj++)
		delete (*neweltsother[i])[jj];
	      delete neweltsother[i];
	    }
	
	}
    }

  PrintMessage (5, cnt, " swaps performed");


  for(int i=0; i<locidmaps.Size(); i++)
    delete locidmaps[i];


  mesh.Compress ();

  multithread.task = savetask;
}
  







/*
  2 -> 3 conversion
*/



void MeshOptimize3d :: SwapImprove2 (Mesh & mesh, OPTIMIZEGOAL goal)
{
  int j, k, l;
  ElementIndex ei, eli1, eli2, elnr;
  SurfaceElementIndex sei;
  PointIndex pi1(0), pi2(0), pi3(0), pi4(0), pi5(0);

  int cnt = 0;

  Element el21(TET), el22(TET), el31(TET), el32(TET), el33(TET);

  double bad1, bad2;

  int np = mesh.GetNP();
  int ne = mesh.GetNE();
  int nse = mesh.GetNSE();

  if (goal == OPT_CONFORM) return;

  // contains at least all elements at node
  TABLE<ElementIndex, PointIndex::BASE> elementsonnode(np); 
  TABLE<SurfaceElementIndex, PointIndex::BASE> belementsonnode(np);

  PrintMessage (3, "SwapImprove2 ");
  (*testout) << "\n" << "Start SwapImprove2" << "\n";
  //  TestOk();



  /*
    CalcSurfacesOfNode ();
    for (i = 1; i <= GetNE(); i++)
    if (volelements.Get(i)[0])
    if (!mesh.LegalTet (volelements.Get(i)))
    {
    cout << "detected illegal tet, 1" << endl;
    (*testout) << "detected illegal tet1: " << i << endl;
    }
  */

  
  // Calculate total badness

  bad1 = CalcTotalBad (mesh.Points(), mesh.VolumeElements());
  (*testout) << "Total badness = " << bad1 << endl;
  //  cout << "tot bad = " << bad1 << endl;

  // find elements on node

  for (ei = 0; ei < ne; ei++)
    for (j = 0; j < mesh[ei].GetNP(); j++)
      elementsonnode.Add (mesh[ei][j], ei);

  for (sei = 0; sei < nse; sei++)
    for (j = 0; j < 3; j++)
      belementsonnode.Add (mesh[sei][j], sei);

  for (eli1 = 0; eli1 < ne; eli1++)
    {
      if (multithread.terminate)
	break;

      if (mesh.ElementType (eli1) == FIXEDELEMENT)
	continue;

      if (mesh[eli1].GetType() != TET)
	continue;

      if ((goal == OPT_LEGAL) && 
	  mesh.LegalTet (mesh[eli1]) &&
	  CalcBad (mesh.Points(), mesh[eli1], 0) < 1e3)
	continue;

      // cout << "eli = " << eli1 << endl;
      //      (*testout) << "swapimp2, eli = " << eli1 << "; el = " << mesh[eli1] << endl;

      for (j = 0; j < 4; j++)
	{
	  // loop over faces
      
	  Element & elem = mesh[eli1];
	  // if (elem[0] < PointIndex::BASE) continue;
	  if (elem.IsDeleted()) continue;

	  int mattyp = elem.GetIndex();
	  
	  switch (j)
	    {
	    case 0:
	      pi1 = elem.PNum(1); pi2 = elem.PNum(2); 
	      pi3 = elem.PNum(3); pi4 = elem.PNum(4);
	      break;
	    case 1:
	      pi1 = elem.PNum(1); pi2 = elem.PNum(4); 
	      pi3 = elem.PNum(2); pi4 = elem.PNum(3);
	      break;
	    case 2:
	      pi1 = elem.PNum(1); pi2 = elem.PNum(3); 
	      pi3 = elem.PNum(4); pi4 = elem.PNum(2);
	      break;
	    case 3:
	      pi1 = elem.PNum(2); pi2 = elem.PNum(4); 
	      pi3 = elem.PNum(3); pi4 = elem.PNum(1);
	      break;
	    }
	  

	  bool bface = 0;
	  for (k = 0; k < belementsonnode[pi1].Size(); k++)
	    {
	      const Element2d & bel = 
		mesh[belementsonnode[pi1][k]];

	      bool bface1 = 1;
	      for (l = 0; l < 3; l++)
		if (bel[l] != pi1 && bel[l] != pi2 && bel[l] != pi3)
		  {
		    bface1 = 0;
		    break;
		  }

	      if (bface1) 
		{
		  bface = 1;
		  break;
		}
	    }
	  
	  if (bface) continue;


	  FlatArray<ElementIndex> row = elementsonnode[pi1];
	  for (k = 0; k < row.Size(); k++)
	    {
	      eli2 = row[k];

	      // cout << "\rei1 = " << eli1 << ", pi1 = " << pi1 << ", k = " << k << ", ei2 = " << eli2 
	      // << ", getne = " << mesh.GetNE();

	      if ( eli1 != eli2 )
		{
		  Element & elem2 = mesh[eli2];
		  if (elem2.IsDeleted()) continue;
		  if (elem2.GetType() != TET)
		    continue;
		  
		  int comnodes=0;
		  for (l = 1; l <= 4; l++)
		    if (elem2.PNum(l) == pi1 || elem2.PNum(l) == pi2 ||
			elem2.PNum(l) == pi3)
		      {
			comnodes++;
		      }
		    else
		      {
			pi5 = elem2.PNum(l);
		      }
		  
		  if (comnodes == 3)
		    {
		      bad1 = CalcBad (mesh.Points(), elem, 0) + 
			CalcBad (mesh.Points(), elem2, 0); 
		      
		      if (!mesh.LegalTet(elem) || 
			  !mesh.LegalTet(elem2))
			bad1 += 1e4;

		      
		      el31.PNum(1) = pi1;
		      el31.PNum(2) = pi2;
		      el31.PNum(3) = pi5;
		      el31.PNum(4) = pi4;
		      el31.SetIndex (mattyp);
		      
		      el32.PNum(1) = pi2;
		      el32.PNum(2) = pi3;
		      el32.PNum(3) = pi5;
		      el32.PNum(4) = pi4;
		      el32.SetIndex (mattyp);
		      
		      el33.PNum(1) = pi3;
		      el33.PNum(2) = pi1;
		      el33.PNum(3) = pi5;
		      el33.PNum(4) = pi4;
		      el33.SetIndex (mattyp);
		      
		      bad2 = CalcBad (mesh.Points(), el31, 0) + 
			CalcBad (mesh.Points(), el32, 0) +
			CalcBad (mesh.Points(), el33, 0); 
		      

		      el31.flags.illegal_valid = 0;
		      el32.flags.illegal_valid = 0;
		      el33.flags.illegal_valid = 0;

		      if (!mesh.LegalTet(el31) || 
			  !mesh.LegalTet(el32) ||
			  !mesh.LegalTet(el33))
			bad2 += 1e4;


		      bool do_swap = (bad2 < bad1);

		      if ( ((bad2 < 1e6) || (bad2 < 10 * bad1)) &&
			   mesh.BoundaryEdge (pi4, pi5))
			do_swap = 1;
			   
		      if (do_swap)
			{
			  //			  cout << "do swap, eli1 = " << eli1 << "; eli2 = " << eli2 << endl;
			  //			  (*mycout) << "2->3 " << flush;
			  cnt++;

			  el31.flags.illegal_valid = 0;
			  el32.flags.illegal_valid = 0;
			  el33.flags.illegal_valid = 0;

			  mesh[eli1] = el31;
			  mesh[eli2] = el32;
			  
			  ElementIndex neli =
			    mesh.AddVolumeElement (el33);
			  
			  /*
			    if (!LegalTet(el31) || !LegalTet(el32) ||
			    !LegalTet(el33))
			    {
			    cout << "Swap to illegal tets !!!" << endl;
			    }
			  */
			  // cout << "neli = " << neli << endl;
			  for (l = 0; l < 4; l++)
			    {
			      elementsonnode.Add (el31[l], eli1);
			      elementsonnode.Add (el32[l], eli2);
			      elementsonnode.Add (el33[l], neli);
			    }

			  break;
			}
		    }
		}
	    }
	}
    }


  PrintMessage (5, cnt, " swaps performed");



  /*
    CalcSurfacesOfNode ();
    for (i = 1; i <= GetNE(); i++)
    if (volelements.Get(i).PNum(1))
    if (!LegalTet (volelements.Get(i)))
    {
    cout << "detected illegal tet, 2" << endl;
    (*testout) << "detected illegal tet2: " << i << endl;
    }
  */


  bad1 = CalcTotalBad (mesh.Points(), mesh.VolumeElements());
  (*testout) << "Total badness = " << bad1 << endl;
  (*testout) << "swapimprove2 done" << "\n";
  //  (*mycout) << "Vol = " << CalcVolume (points, volelements) << "\n";
}


/*
  void Mesh :: SwapImprove2 (OPTIMIZEGOAL goal)
  {
  int i, j;
  int eli1, eli2;
  int mattyp;

  Element el31(4), el32(4), el33(4);
  double bad1, bad2;


  INDEX_3_HASHTABLE<INDEX_2> elsonface (GetNE());

  (*mycout) << "SwapImprove2 " << endl;
  (*testout) << "\n" << "Start SwapImprove2" << "\n";

  // Calculate total badness

  if (goal == OPT_QUALITY)
  {
  double bad1 = CalcTotalBad (points, volelements);
  (*testout) << "Total badness = " << bad1 << endl;
  }

  // find elements on node


  Element2d face;
  for (i = 1; i <= GetNE(); i++)
  if ( (i > eltyps.Size()) || (eltyps.Get(i) != FIXEDELEMENT) )
  {
  const Element & el = VolumeElement(i);
  if (!el.PNum(1)) continue;

  for (j = 1; j <= 4; j++)
  {
  el.GetFace (j, face);
  INDEX_3 i3 (face.PNum(1), face.PNum(2), face.PNum(3));
  i3.Sort();


  int bnr, posnr;
  if (!elsonface.PositionCreate (i3, bnr, posnr))
  {
  INDEX_2 i2;
  elsonface.GetData (bnr, posnr, i3, i2);
  i2.I2() = i;
  elsonface.SetData (bnr, posnr, i3, i2);
  }
  else
  {
  INDEX_2 i2 (i, 0);
  elsonface.SetData (bnr, posnr, i3, i2);
  }

  //  	    if (elsonface.Used (i3))
  //  	      {
  //  		INDEX_2 i2 = elsonface.Get(i3);
  //  		i2.I2() = i;
  //  		elsonface.Set (i3, i2);
  //  	      }
  //  	    else
  //  	      {
  //  		INDEX_2 i2 (i, 0);
  //  		elsonface.Set (i3, i2);
  //  	      }

  }
  }

  BitArray original(GetNE());
  original.Set();

  for (i = 1; i <= GetNSE(); i++)
  {
  const Element2d & sface = SurfaceElement(i);
  INDEX_3 i3 (sface.PNum(1), sface.PNum(2), sface.PNum(3));
  i3.Sort();
  INDEX_2 i2(0,0);
  elsonface.Set (i3, i2);
  }


  for (i = 1; i <= elsonface.GetNBags(); i++)
  for (j = 1; j <= elsonface.GetBagSize(i); j++)
  {
  INDEX_3 i3;
  INDEX_2 i2;
  elsonface.GetData (i, j, i3, i2);


  int eli1 = i2.I1();
  int eli2 = i2.I2();

  if (eli1 && eli2 && original.Test(eli1) && original.Test(eli2) )
  {
  Element & elem = volelements.Elem(eli1);
  Element & elem2 = volelements.Elem(eli2);

  int pi1 = i3.I1();
  int pi2 = i3.I2();
  int pi3 = i3.I3();

  int pi4 = elem.PNum(1) + elem.PNum(2) + elem.PNum(3) + elem.PNum(4) - pi1 - pi2 - pi3;
  int pi5 = elem2.PNum(1) + elem2.PNum(2) + elem2.PNum(3) + elem2.PNum(4) - pi1 - pi2 - pi3;






  el31.PNum(1) = pi1;
  el31.PNum(2) = pi2;
  el31.PNum(3) = pi3;
  el31.PNum(4) = pi4;
  el31.SetIndex (mattyp);
	    
  if (WrongOrientation (points, el31))
  swap (pi1, pi2);


  bad1 = CalcBad (points, elem, 0) + 
  CalcBad (points, elem2, 0); 
	    
  //	    if (!LegalTet(elem) || !LegalTet(elem2))
  //	      bad1 += 1e4;

	    
  el31.PNum(1) = pi1;
  el31.PNum(2) = pi2;
  el31.PNum(3) = pi5;
  el31.PNum(4) = pi4;
  el31.SetIndex (mattyp);
	    
  el32.PNum(1) = pi2;
  el32.PNum(2) = pi3;
  el32.PNum(3) = pi5;
  el32.PNum(4) = pi4;
  el32.SetIndex (mattyp);
		      
  el33.PNum(1) = pi3;
  el33.PNum(2) = pi1;
  el33.PNum(3) = pi5;
  el33.PNum(4) = pi4;
  el33.SetIndex (mattyp);
	    
  bad2 = CalcBad (points, el31, 0) + 
  CalcBad (points, el32, 0) +
  CalcBad (points, el33, 0); 
	    
  //	    if (!LegalTet(el31) || !LegalTet(el32) ||
  //		!LegalTet(el33))
  //	      bad2 += 1e4;
	    
	    
  int swap = (bad2 < bad1);

  INDEX_2 hi2b(pi4, pi5);
  hi2b.Sort();
	    
  if ( ((bad2 < 1e6) || (bad2 < 10 * bad1)) &&
  boundaryedges->Used (hi2b) )
  swap = 1;
	    
  if (swap)
  {
  (*mycout) << "2->3 " << flush;
		
  volelements.Elem(eli1) = el31;
  volelements.Elem(eli2) = el32;
  volelements.Append (el33);
		
  original.Clear (eli1);
  original.Clear (eli2);
  }
  }
  }
  
  (*mycout) << endl;

  if (goal == OPT_QUALITY)
  {
  bad1 = CalcTotalBad (points, volelements);
  (*testout) << "Total badness = " << bad1 << endl;
  }

  //  FindOpenElements ();

  (*testout) << "swapimprove2 done" << "\n";
  }

*/
}
