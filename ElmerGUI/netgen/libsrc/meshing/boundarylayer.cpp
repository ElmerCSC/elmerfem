#include <mystdlib.h>
#include "meshing.hpp"

namespace netgen
{

void InsertVirtualBoundaryLayer (Mesh & mesh)
{
  cout << "Insert virt. b.l." << endl;
  
  int surfid;

  cout << "Boundary Nr:";
  cin >> surfid;

  int i, j;
  int np = mesh.GetNP();

  cout << "Old NP: " << mesh.GetNP() << endl;
  cout << "Trigs: " << mesh.GetNSE() << endl;

  BitArray bndnodes(np);
  ARRAY<int> mapto(np);

  bndnodes.Clear();
  for (i = 1; i <= mesh.GetNSeg(); i++)
    {
      int snr = mesh.LineSegment(i).edgenr;
      cout << "snr = " << snr << endl;
      if (snr == surfid)
	{
	  bndnodes.Set (mesh.LineSegment(i).p1);
	  bndnodes.Set (mesh.LineSegment(i).p2);
	}
    }
  for (i = 1; i <= mesh.GetNSeg(); i++)
    {
      int snr = mesh.LineSegment(i).edgenr;
      if (snr != surfid)
	{
	  bndnodes.Clear (mesh.LineSegment(i).p1);
	  bndnodes.Clear (mesh.LineSegment(i).p2);
	}
    }
  
  for (i = 1; i <= np; i++)
    {
      if (bndnodes.Test(i))
	mapto.Elem(i) = mesh.AddPoint (mesh.Point (i));
      else
	mapto.Elem(i) = 0;
    }

  for (i = 1; i <= mesh.GetNSE(); i++)
    {
      Element2d & el = mesh.SurfaceElement(i);
      for (j = 1; j <= el.GetNP(); j++)
	if (mapto.Get(el.PNum(j)))
	  el.PNum(j) = mapto.Get(el.PNum(j));
    }


  int nq = 0;
  for (i = 1; i <= mesh.GetNSeg(); i++)
    {
      int snr = mesh.LineSegment(i).edgenr;
      if (snr == surfid)
	{
	  int p1 = mesh.LineSegment(i).p1;
	  int p2 = mesh.LineSegment(i).p2;
	  int p3 = mapto.Get (p1);
	  if (!p3) p3 = p1;
	  int p4 = mapto.Get (p2);
	  if (!p4) p4 = p2;
	  
	  Element2d el(QUAD);
	  el.PNum(1) = p1;
	  el.PNum(2) = p2;
	  el.PNum(3) = p3;
	  el.PNum(4) = p4;
	  el.SetIndex (2);
	  mesh.AddSurfaceElement (el);
	  nq++;
	}
    }

  cout << "New NP: " << mesh.GetNP() << endl;
  cout << "Quads: " << nq << endl;
}

}

