#include <mystdlib.h>
#include "meshing.hpp"

namespace netgen
{

netrule :: netrule ()
{
  name = new char[1];
  name[0] = char(0);
  quality = 0;
}

netrule ::  ~netrule()
{
  // if(name != NULL) 
  delete [] name;
  for(int i=0; i<oldutofreearea_i.Size(); i++)
    delete oldutofreearea_i[i];
}


/*
void netrule :: GetFreeArea (ARRAY<Point2d> & afreearea)
  {
  int i;

  afreearea.SetSize (freearea.Size());
  for (i = 1; i <= freearea.Size(); i++)
    afreearea[i] = freearea[i];
  }
*/


void netrule :: SetFreeZoneTransformation (const Vector & devp, int tolclass)
{
  double lam1 = 1.0/tolclass;
  double lam2 = 1.-lam1;

  double mem1[100], mem2[100], mem3[100];

  int vs = oldutofreearea.Height();
  FlatVector devfree(vs, mem1);
  FlatVector devfree1(vs, mem2);
  FlatVector devfree2(vs, mem3);

  if (tolclass <= oldutofreearea_i.Size())
    {
      oldutofreearea_i[tolclass-1] -> Mult (devp, devfree);
    }
  else
    {
      oldutofreearea.Mult (devp, devfree1);
      oldutofreearealimit.Mult (devp, devfree2);
      devfree.Set2 (lam1, devfree1, lam2, devfree2);
    }

  
  int fzs = freezone.Size();
  transfreezone.SetSize (fzs);

  if (fzs > 0)
    {
      transfreezone[0].X() = lam1 * freezone[0].X() + lam2 * freezonelimit[0].X() + devfree[0];
      transfreezone[0].Y() = lam1 * freezone[0].Y() + lam2 * freezonelimit[0].Y() + devfree[1];
      fzmaxx = fzminx = transfreezone[0].X();
      fzmaxy = fzminy = transfreezone[0].Y();
    }

  for (int i = 1; i < fzs; i++)
    {
      transfreezone[i].X() = lam1 * freezone[i].X() + lam2 * freezonelimit[i].X() + devfree[2*i];
      transfreezone[i].Y() = lam1 * freezone[i].Y() + lam2 * freezonelimit[i].Y() + devfree[2*i+1];

      if (transfreezone[i].X() > fzmaxx) fzmaxx = transfreezone[i].X();
      if (transfreezone[i].X() < fzminx) fzminx = transfreezone[i].X();
      if (transfreezone[i].Y() > fzmaxy) fzmaxy = transfreezone[i].Y();
      if (transfreezone[i].Y() < fzminy) fzminy = transfreezone[i].Y();
    }

  for (int i = 0; i < fzs; i++)
    {
      Point2d p1 = transfreezone[i];
      Point2d p2 = transfreezone[(i+1) % fzs];

      Vec2d vn (p2.Y() - p1.Y(), p1.X() - p2.X());

      double len2 = vn.Length2();

      if (len2 < 1e-10)
	{
	  freesetinequ(i, 0) = 0;
	  freesetinequ(i, 1) = 0;
	  freesetinequ(i, 2) = -1;
	}
      else
	{
	  vn /= sqrt (len2);    // should not be necessary

	  freesetinequ(i,0) = vn.X(); 
	  freesetinequ(i,1) = vn.Y(); 
	  freesetinequ(i,2) = -(p1.X() * vn.X() + p1.Y() * vn.Y());
	}

      /*
      freesetinequ(i,0) = vn.X(); 
      freesetinequ(i,1) = vn.Y(); 
      freesetinequ(i,2) = -(p1.X() * vn.X() + p1.Y() * vn.Y());
      */
    }
}


/*
int netrule :: IsInFreeZone2 (const Point2d & p) const
{
  for (int i = 0; i < transfreezone.Size(); i++)
    {
      if (freesetinequ(i, 0) * p.X() + 
	  freesetinequ(i, 1) * p.Y() +
	  freesetinequ(i, 2) > 0) return 0;
    }
  return 1;
}
*/

int netrule :: IsLineInFreeZone2 (const Point2d & p1, const Point2d & p2) const
{
  int left, right, allleft, allright;

  if (p1.X() > fzmaxx && p2.X() > fzmaxx ||
      p1.X() < fzminx && p2.X() < fzminx ||
      p1.Y() > fzmaxy && p2.Y() > fzmaxy ||
      p1.Y() < fzminy && p2.Y() < fzminy) return 0;

  for (int i = 1; i <= transfreezone.Size(); i++)
    {
      if (freesetinequ.Get(i, 1) * p1.X() + freesetinequ.Get(i, 2) * p1.Y() +
	  freesetinequ.Get(i, 3) > -1e-8 &&    // -1e-6
	  freesetinequ.Get(i, 1) * p2.X() + freesetinequ.Get(i, 2) * p2.Y() +
	  freesetinequ.Get(i, 3) > -1e-8       // -1e-6
	  ) return 0;
    }

  double nx =  (p2.Y() - p1.Y());
  double ny = -(p2.X() - p1.X());
  double nl = sqrt (nx * nx + ny * ny);
  if (nl > 1e-8)
    {
      nx /= nl;
      ny /= nl;
      double c = - (p1.X() * nx + p1.Y() * ny);

      allleft = 1;
      allright = 1;

      for (int i = 1; i <= transfreezone.Size(); i++)
	{
	  left  = transfreezone.Get(i).X() * nx + transfreezone.Get(i).Y() + c <  1e-7;
	  right = transfreezone.Get(i).X() * nx + transfreezone.Get(i).Y() + c > -1e-7;

	  if (!left) allleft = 0;
	  if (!right) allright = 0;
	}
      if (allleft || allright) return 0;
    }

  return 1;
}

int netrule :: ConvexFreeZone () const
{
  int n = transfreezone.Size();
  for (int i = 1; i <= n; i++)
    {
      const bool counterclockwise = CCW (transfreezone.Get(i), 
					 transfreezone.Get(i % n + 1),
					 transfreezone.Get( (i+1) % n + 1 ),
					 1e-7);
      //(*testout) << "ccw " << counterclockwise << endl << " p1 " << transfreezone.Get(i) << " p2 " << transfreezone.Get(i % n + 1)
      //		 << " p3 " << transfreezone.Get( (i+1) % n + 1 ) << endl;
      if (!counterclockwise )
	return 0;
    }
  return 1;
}


/*
float netrule :: CalcPointDist (int pi, const Point2d & p) const
{
  float dx = p.X() - points.Get(pi).X();
  float dy = p.Y() - points.Get(pi).Y();
  const threefloat * tf = &tolerances.Get(pi);

  return tf->f1 * dx * dx + tf->f2 * dx * dy + tf->f3 * dy * dy;
}
*/

float netrule :: CalcLineError (int li, const Vec2d & v) const
{
  float dx = v.X() - linevecs.Get(li).X();
  float dy = v.Y() - linevecs.Get(li).Y();

  const threefloat * ltf = &linetolerances.Get(li);
  return ltf->f1 * dx * dx + ltf->f2 * dx * dy + ltf->f3 * dy * dy;
}




/*
int GetNRules ()
  {
  return rules.Size();
  }
*/











}
