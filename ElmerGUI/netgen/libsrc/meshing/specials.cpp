#include <mystdlib.h>
#include "meshing.hpp"


namespace netgen
{

// A special function for Hermann Landes, Erlangen


void CutOffAndCombine (Mesh & mesh, const Mesh & othermesh)
{
  int i, j;
  int nse = othermesh.GetNSE();
  int onp = othermesh.GetNP();

  int ne = mesh.GetNE();

  PrintMessage (1, "other mesh has ",
		othermesh.GetNP(), " points, ",
		othermesh.GetNSE(), " surface elements.");

  ARRAY<Box3d> otherbounds(nse);  
  Box3d otherbox;

  double maxh = 0;
  for (i = 1; i <= nse; i++)
    {
      const Element2d & sel = othermesh.SurfaceElement(i);
      sel.GetBox(othermesh.Points(), otherbounds.Elem(i));

      double loch = othermesh.GetH (othermesh.Point (sel.PNum(1)));
      otherbounds.Elem(i).Increase(loch);
      if (loch > maxh) maxh = loch;
    }

  otherbox.SetPoint (othermesh.Point(1));
  for (i = 1; i <= othermesh.GetNP(); i++)
    otherbox.AddPoint (othermesh.Point(i));
  otherbox.Increase (maxh);

  for (i = 1; i <= ne; i++)
    {
      Box3d box;
      int remove = 0;

      const Element & el = mesh.VolumeElement(i);
      el.GetBox(mesh.Points(), box);

      if (i % 10000 == 0)
	cout << "+" << flush;

      if (box.Intersect(otherbox))
	{
	  for (j = 1; j <= nse && !remove; j++)
	    if (box.Intersect(otherbounds.Get(j)))
	      remove = 1;
	}

      if (remove)
	mesh.VolumeElement(i).Delete();
    }
  cout << endl;

  BitArray connected(mesh.GetNP());
  connected.Clear();
  for (i = 1; i <= mesh.GetNSE(); i++)
    {
      const Element2d & el = mesh.SurfaceElement(i);
      for (j = 1; j <= 3; j++)
	connected.Set(el.PNum(j));
    }
  
  bool changed;
  do
    {
      changed = 0;
      for (i = 1; i <= mesh.GetNE(); i++)
	{
	  const Element & el = mesh.VolumeElement(i);
	  int has = 0, hasnot = 0;
	  if (el[0])
	    {
	      for (j = 0; j < 4; j++)
		{
		  if (connected.Test(el[j]))
		    has = 1;
		  else
		    hasnot = 1;
		}
	      if (has && hasnot)
		{
		  changed = 1;
		  for (j = 0; j < 4; j++)
		    connected.Set (el[j]);
		}
	    }
	}
      cout << "." << flush;
    }
  while (changed);
  cout << endl;

  for (i = 1; i <= mesh.GetNE(); i++)
    {
      const Element & el = mesh.VolumeElement(i);
      int hasnot = 0;
      if (el[0])
	{
	  for (j = 0; j < 4; j++)
	    {
	      if (!connected.Test(el[j]))
		hasnot = 1;
	    }
	  if (hasnot)
	    mesh.VolumeElement(i).Delete();
	}
    }

  mesh.Compress();
  
  mesh.FindOpenElements();
  BitArray locked(mesh.GetNP());
  locked.Set();
  for (i = 1; i <= mesh.GetNOpenElements(); i++)
    for (j = 1; j <= 3; j++)
      locked.Clear (mesh.OpenElement(i).PNum(j));

  for (i = 1; i <= locked.Size(); i++)
    if (locked.Test(i))
      {
	mesh.AddLockedPoint (i);
      }



  
  ARRAY<int> pmat(onp);

  for (i = 1; i <= onp; i++)
    pmat.Elem(i) = mesh.AddPoint (othermesh.Point(i));

  int fnum = 
    mesh.AddFaceDescriptor (FaceDescriptor(0,0,1,0));

  for (i = 1; i <= othermesh.GetNSE(); i++)
    {
      Element2d tri = othermesh.SurfaceElement(i);
      for (j = 1; j <= 3; j++)
	tri.PNum(j) = pmat.Get(tri.PNum(j));
      tri.SetIndex(fnum);
      mesh.AddSurfaceElement (tri);
    }

  for (i = 1; i <= onp; i++)
    mesh.AddLockedPoint (pmat.Elem(i));

  mesh.CalcSurfacesOfNode();
  mesh.CalcLocalH();
}




void HelmholtzMesh (Mesh & mesh)
{
  int i;
  double ri, ra, rinf;

  cout << "ri = ";
  cin >> ri;
  cout << "ra = ";
  cin >> ra;
  cout << "rinf = ";
  cin >> rinf;

  double det = ri * ra * rinf - ri * ri * rinf;
  double a = (ri - rinf) / det;
  double b = (ri*ri - ra * rinf) / det;
  for (i = 1; i <= mesh.GetNP(); i++)
    {
      Point<3> & p = mesh.Point(i);
      double rold = sqrt (sqr(p(0)) + sqr(p(1)) + sqr(p(2)));
      if (rold < ri) continue;

      double rnew = 1 / (a * rold - b);
      double fac = rnew / rold;
      p(0) *= fac;
      p(1) *= fac;
      p(2) *= fac;
    }
}
}
