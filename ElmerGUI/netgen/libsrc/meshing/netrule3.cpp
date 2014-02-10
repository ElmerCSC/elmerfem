#include <mystdlib.h>
#include "meshing.hpp"

// #define MARK
// #include <prof.h>


namespace netgen
{


vnetrule :: vnetrule ()
{
  name = new char[1];
  name[0] = char(0);
  quality = 0;
}

vnetrule :: ~vnetrule ()
{
  // if (strlen(name)) 
  delete [] name;
  for (int i = 1; i <= freefaces.Size(); i++)
    delete freefaces.Elem(i);
  for (int i = 1; i <= freesets.Size(); i++)
    delete freesets.Elem(i);
  for (int i = 1; i <= freeedges.Size(); i++)
    delete freeedges.Elem(i);
  for (int i = 1; i <= freefaceinequ.Size(); i++)
    delete freefaceinequ.Elem(i);
  delete oldutofreezone;
  delete oldutofreezonelimit;
}

int vnetrule :: TestFlag (char flag) const
{
  for (int i = 1; i <= flags.Size(); i++)
    if (flags.Get(i) == flag) return 1;
  return 0;
}


void vnetrule :: SetFreeZoneTransformation (const Vector & allp, int tolclass)
{
  int i, j;
  // double nx, ny, nz, v1x, v1y, v1z, v2x, v2y, v2z;
  double nl;
  const threeint * ti;
  int fs;

  double lam1 = 1.0/(2 * tolclass - 1);
  double lam2 = 1-lam1;

  transfreezone.SetSize (freezone.Size());
  
  int np = points.Size();
  int nfp = freezone.Size();
  Vector vp(np), vfp1(nfp), vfp2(nfp);


  for (i = 1; i <= 3; i++)
    {
      for (j = 1; j <= np; j++)
	vp.Elem(j) = allp.Get(i+3*j-3);

      oldutofreezone->Mult (vp, vfp1);
      oldutofreezonelimit->Mult (vp, vfp2);

      vfp1 *= lam1;
      vfp1.Add (lam2, vfp2);

      for (j = 1; j <= nfp; j++)
	transfreezone.Elem(j).X(i) = vfp1.Elem(j);
    }

  // MARK(setfz2);


  fzbox.SetPoint (transfreezone.Elem(1));
  for (i = 2; i <= freezone.Size(); i++)
    fzbox.AddPoint (transfreezone.Elem(i));
  
  
  // MARK(setfz3);


  for (fs = 1; fs <= freesets.Size(); fs++)
    {
      ARRAY<threeint> & freesetfaces = *freefaces.Get(fs);
      DenseMatrix & freesetinequ = *freefaceinequ.Get(fs);
      
      for (i = 1; i <= freesetfaces.Size(); i++)
	{
	  ti = &freesetfaces.Get(i);
	  const Point3d & p1 = transfreezone.Get(ti->i1);
	  const Point3d & p2 = transfreezone.Get(ti->i2);
	  const Point3d & p3 = transfreezone.Get(ti->i3);

	  Vec3d v1(p1, p2);   
	  Vec3d v2(p1, p3);   
	  Vec3d n;
	  Cross (v1, v2, n);

	  nl = n.Length();

	  if (nl < 1e-10)
	    {
	      freesetinequ.Set(1, 1, 0);
	      freesetinequ.Set(1, 2, 0);
	      freesetinequ.Set(1, 3, 0);
	      freesetinequ.Set(1, 4, -1);
	    }
	  else
	    {
	      //	      n /= nl;
	      
	      freesetinequ.Set(i, 1, n.X()/nl);
	      freesetinequ.Set(i, 2, n.Y()/nl);
	      freesetinequ.Set(i, 3, n.Z()/nl);
	      freesetinequ.Set(i, 4,
			       -(p1.X() * n.X() + p1.Y() * n.Y() + p1.Z() * n.Z()) / nl);
	    }
	}
    }

  /*
  (*testout) << "Transformed freezone: " << endl;
  for (i = 1; i <= transfreezone.Size(); i++)
    (*testout) << transfreezone.Get(i) << " ";
  (*testout) << endl;
  */
}

int vnetrule :: ConvexFreeZone () const
{
  int i, j, k, fs;

  // (*mycout) << "Convex free zone...\n";
  
  int ret1=1;
  // int ret2=1;

  for (fs = 1; fs <= freesets.Size(); fs++)
    {
      const DenseMatrix & freesetinequ = *freefaceinequ.Get(fs);

      // const ARRAY<int> & freeset = *freesets.Get(fs);
      const ARRAY<twoint> & freesetedges = *freeedges.Get(fs);
      // const ARRAY<threeint> & freesetfaces = *freefaces.Get(fs);
      
      for (i = 1; i <= freesetedges.Size(); i++)
	{
	  j = freesetedges.Get(i).i1;    //triangle j with opposite point k
	  k = freesetedges.Get(i).i2;
	  
	  if ( freesetinequ.Get(j, 1) * transfreezone.Get(k).X() +
	       freesetinequ.Get(j, 2) * transfreezone.Get(k).Y() +
	       freesetinequ.Get(j, 3) * transfreezone.Get(k).Z() +
	       freesetinequ.Get(j, 4) > 0 )
	    {
	      ret1=0;
	    }
	}
      
    }

  return ret1;
}


int vnetrule :: IsInFreeZone (const Point3d & p) const
{
  int i, fs;
  char inthis;
  
  
  for (fs = 1; fs <= freesets.Size(); fs++)
    {
      inthis = 1;
      ARRAY<threeint> & freesetfaces = *freefaces.Get(fs);
      DenseMatrix & freesetinequ = *freefaceinequ.Get(fs);
      
      for (i = 1; i <= freesetfaces.Size() && inthis; i++)
	{
	  if (freesetinequ.Get(i, 1) * p.X() + freesetinequ.Get(i, 2) * p.Y() +
	      freesetinequ.Get(i, 3) * p.Z() + freesetinequ.Get(i, 4) > 0)
	    inthis = 0;
	}
      
      if (inthis) return 1;
    }
  
  return 0;
}


int vnetrule :: IsTriangleInFreeZone (const Point3d & p1, 
				      const Point3d & p2,
				      const Point3d & p3, 
				      const ARRAY<int> & pi, int newone)
{
  int fs;
  int infreeset, cannot = 0;


  static ARRAY<int> pfi(3), pfi2(3);

  // convert from local index to freeset index
  int i, j;
  for (i = 1; i <= 3; i++)
    {
      pfi.Elem(i) = 0;
      if (pi.Get(i))
	{
	  for (j = 1; j <= freezonepi.Size(); j++)
	    if (freezonepi.Get(j) == pi.Get(i))
	      pfi.Elem(i) = j;
	}
    }

  for (fs = 1; fs <= freesets.Size(); fs++)
    {
      const ARRAY<int> & freeseti = *freesets.Get(fs);
      for (i = 1; i <= 3; i++)
	{
	  pfi2.Elem(i) = 0;
	  for (j = 1; j <= freeseti.Size(); j++)
	    if (pfi.Get(i) == freeseti.Get(j))
	      pfi2.Elem(i) = pfi.Get(i);
	}

      infreeset = IsTriangleInFreeSet(p1, p2, p3, fs, pfi2, newone);
      if (infreeset == 1) return 1;
      if (infreeset == -1) cannot = -1;
    }
  
  return cannot;
}



int vnetrule :: IsTriangleInFreeSet (const Point3d & p1, const Point3d & p2,
                                     const Point3d & p3, int fs,
				     const ARRAY<int> & pi, int newone)
{
  int i, ii;
  Vec3d n;
  int allleft, allright;
  int hos1, hos2, hos3, os1, os2, os3;
  double hf, lam1, lam2, f, c1, c2, alpha;
  double v1n, v2n, h11, h12, h22, dflam1, dflam2;
  double lam1old, lam2old, fold;
  double hpx, hpy, hpz, v1x, v1y, v1z, v2x, v2y, v2z;
  int act1, act2, act3, it;
  int cntout;
  static ARRAY<int> activefaces;
  int isin;
  

  // MARK(triinfz);
  
  ARRAY<threeint> & freesetfaces = *freefaces.Get(fs);
  DenseMatrix & freesetinequ = *freefaceinequ.Get(fs);
  

  int cnt = 0;
  for (i = 1; i <= 3; i++)
    if (pi.Get(i)) cnt++;

  /*
  (*testout) << "trig in free set : " << p1 << " - " << p2 << " - " << p3 << endl;
  (*testout) << "common points: " << cnt << endl;
  */
  if (!newone)
    cnt = 0;

  if (cnt == 1)
    {
      // MARK(triinfz1);

      int upi = 0, lpiu = 0;
      for (i = 1; i <= 3; i++)
	if (pi.Get(i))
	  {
	    upi = i;
	    lpiu = pi.Get(i);
	  }

      Vec3d v1, v2;
      switch (upi)
	{
	case 1:
	  {
	    v1 = p2 - p1;
	    v2 = p3 - p1;
	    break;
	  }
	case 2:
	  {
	    v1 = p3 - p2;
	    v2 = p1 - p2;
	    break;
	  }
	case 3:
	  {
	    v1 = p1 - p3;
	    v2 = p2 - p3;
	    break;
	  }
	}

      v1 /= v1.Length();
      v2 /= v2.Length();
      Cross (v1, v2, n);
      n /= n.Length();

      //      (*testout) << "Test new: " << endl;
      for (i = 1; i <= freesetfaces.Size(); i++)
	{
	  if ( (freesetfaces.Get(i).i1 == lpiu) || 
	       (freesetfaces.Get(i).i2 == lpiu) ||
	       (freesetfaces.Get(i).i3 == lpiu) )
	    {
	      // freeface has point


	      Vec3d a (freesetinequ.Get(i, 1),
		       freesetinequ.Get(i, 2),
		       freesetinequ.Get(i, 3));
	      
	      //	      if (1 - fabs (a * n) < 1e-8 ) 
	      //		continue;

	      Vec3d an;
	      Cross (a, n, an);
	      double lan = an.Length();
	      if (lan < 1e-10)
		continue;

	      an /= lan;
	      
	      int out1 = (a * v1) > 0;
	      int out2 = (a * v2) > 0;
	      //	      (*testout) << "out1, out2 = " << out1 << ", " << out2 << endl;
	      if (out1 && out2)
		return 0;

	      if (!out1 && !out2) 
		continue;


	      //	      if ( ( (an * v1) < 0) &&  ( (an * v2) < 0) )   // falsch !!!!
	      //		an *= -1;

	      // solve  an = lam1 v1 + lam2 v2
	      double vii11 = v1 * v1;
	      double vii12 = v1 * v2;
	      double vii22 = v2 * v2;
	      double det = vii11 * vii22 - vii12 * vii12;
	      if ( fabs (det) < 1e-10 )
		continue;
	      double rs1 = an * v1;
	      double rs2 = an * v2;
	      
	      double lambda1 = rs1 * vii22 - rs2 * vii12;
	      double lambda2 = rs2 * vii11 - rs1 * vii12;

	      if (fabs (lambda1) > fabs (lambda2))
		{
		  if (lambda1 < 0)
		    an *= -1;
		}
	      else
		{
		  if (lambda2 < 0)
		    an *= -1;
		}


	      if (lambda1 * lambda2 < 0 && 0)
		{
		  if (fabs (lambda1) > 1e-14 && fabs (lambda2) > 1e-14)
		    {
		      //		      (*mycout) << "lambda1 lambda2 < 0" << endl;
		      (*testout) << "lambdai different" << endl;
		      (*testout) << "v1 = " << v1 << endl;
		      (*testout) << "v2 = " << v2 << endl;
		      (*testout) << "n = " << n << endl;
		      (*testout) << "a = " << a << endl;
		      (*testout) << "an = " << an << endl;
		      (*testout) << "a * v1 = " << (a * v1) << endl;
		      (*testout) << "a * v2 = " << (a * v2) << endl;
		      (*testout) << "an * v1 = " << (an * v1) << endl;
		      (*testout) << "an * v2 = " << (an * v2) << endl;
		      
		      (*testout) << "vii = " << vii11 << ", " << vii12 << ", " << vii22 << endl;
		      (*testout) << "lambdai = " << lambda1 << ", " << lambda2 << endl;
		      (*testout) << "rs = " << rs1 << ", " << rs2 << endl;
		      continue;
		    }
		}

	      if (out1)
		v1 = an;
	      else
		v2 = an;
	    }
	}
      
      return 1;

      /*
      (*testout) << "overlap trig " << p1 << p2 << p3 << endl;
      (*testout) << "upi = " << upi << endl;
      (*testout) << "v1 = " << v1 << " v2 = " << v2 << endl;
      */

      switch (upi)
	{
	case 1:
	  {
	    v1 = p2 - p1;
	    v2 = p3 - p1;
	    break;
	  }
	case 2:
	  {
	    v1 = p3 - p2;
	    v2 = p1 - p2;
	    break;
	  }
	case 3:
	  {
	    v1 = p1 - p3;
	    v2 = p2 - p3;
	    break;
	  }
	}

      v1 /= v1.Length();
      v2 /= v2.Length();
      Cross (v1, v2, n);
      n /= n.Length();

      //      (*testout) << "orig v1, v2 = " << v1 << ", " << v2 << endl;

      
      for (i = 1; i <= freesetfaces.Size(); i++)
	{
	  if ( (freesetfaces.Get(i).i1 == lpiu) || 
	       (freesetfaces.Get(i).i2 == lpiu) ||
	       (freesetfaces.Get(i).i3 == lpiu) )
	    {
	      /*
	      (*testout) << "v1, v2, now = " << v1 << ", " << v2 << endl;

	      // freeface has point
	      (*testout) << "freesetface: "
			 << freesetfaces.Get(i).i1 << " "
			 << freesetfaces.Get(i).i2 << " "
			 << freesetfaces.Get(i).i3 << " ";
	      */

	      Vec3d a (freesetinequ.Get(i, 1),
		       freesetinequ.Get(i, 2),
		       freesetinequ.Get(i, 3));
	      //	      (*testout) << "a = " <<  a << endl;


	      Vec3d an;
	      Cross (a, n, an);
	      double lan = an.Length();
	      
	      //	      (*testout) << "an = " << an << endl;

	      if (lan < 1e-10)
		continue;

	      an /= lan;

	      //	      (*testout) << "a*v1 = " << (a*v1) << " a*v2 = " << (a*v2) << endl;
	      
	      int out1 = (a * v1) > 0;
	      // int out2 = (a * v2) > 0;


	      //	      (*testout) << "out1, 2 = " << out1 << ", " << out2 << endl;

	      
	      double vii11 = v1 * v1;
	      double vii12 = v1 * v2;
	      double vii22 = v2 * v2;
	      double det = vii11 * vii22 - vii12 * vii12;
	      if ( fabs (det) < 1e-10 )
		continue;
	      double rs1 = an * v1;
	      double rs2 = an * v2;
	      
	      double lambda1 = rs1 * vii22 - rs2 * vii12;
	      double lambda2 = rs2 * vii11 - rs1 * vii12;

	      //	      (*testout) << "lambda1, lambda2 = " << lambda1 << ", " << lambda2 << endl;


	      if (fabs (lambda1) > fabs (lambda2))
		{
		  if (lambda1 < 0)
		    an *= -1;
		}
	      else
		{
		  if (lambda2 < 0)
		    an *= -1;
		}


	      if (lambda1 * lambda2 < 0)
		{
		  if (fabs (lambda1) > 1e-14 && fabs (lambda2) > 1e-14)
		    {
		      //		      (*mycout) << "lambda1 lambda2 < 0" << endl;
		      (*testout) << "lambdai different" << endl;
		      (*testout) << "v1 = " << v1 << endl;
		      (*testout) << "v2 = " << v2 << endl;
		      (*testout) << "n = " << n << endl;
		      (*testout) << "a = " << a << endl;
		      (*testout) << "an = " << an << endl;
		      (*testout) << "a * v1 = " << (a * v1) << endl;
		      (*testout) << "a * v2 = " << (a * v2) << endl;
		      (*testout) << "an * v1 = " << (an * v1) << endl;
		      (*testout) << "an * v2 = " << (an * v2) << endl;
		      
		      (*testout) << "vii = " << vii11 << ", " << vii12 << ", " << vii22 << endl;
		      (*testout) << "lambdai = " << lambda1 << ", " << lambda2 << endl;
		      (*testout) << "rs = " << rs1 << ", " << rs2 << endl;
		      continue;
		    }
		}

	      if (out1)
		v1 = an;
	      else
		v2 = an;



	    }
	}

      return 1;
    }



  if (cnt == 2)
    {
      //      (*testout) << "tripoitns: " << p1 << " " << p2 << " " << p3 << endl;

      // MARK(triinfz2);

      int pi1 = 0, pi2 = 0, pi3 = 0;
      Vec3d a1, a2;  // outer normals
      Vec3d trivec;  // vector from common edge to third point of triangle
      for (i = 1; i <= 3; i++)
	if (pi.Get(i))
	  {
	    pi2 = pi1;
	    pi1 = pi.Get(i);
	  }
	else
	  pi3 = i;

      switch (pi3)
	{
	case 1: trivec = (p1 - p2); break;
	case 2: trivec = (p2 - p3); break;
	case 3: trivec = (p3 - p2); break;
	}

      ARRAY<int> lpi(freezonepi.Size());
      for (i = 1; i <= lpi.Size(); i++)
	lpi.Elem(i) = 0;
      lpi.Elem(pi1) = 1;
      lpi.Elem(pi2) = 1;
      
      int ff1 = 0, ff2 = 0;
      for (i = 1; i <= freesetfaces.Size(); i++)
	{
	  if (lpi.Get(freesetfaces.Get(i).i1) + 
	      lpi.Get(freesetfaces.Get(i).i2) + 
	      lpi.Get(freesetfaces.Get(i).i3) == 2)
	    {
	      ff2 = ff1;
	      ff1 = i;
	    }
	}

      if (ff2 == 0)
	return 1;

      a1 = Vec3d (freesetinequ.Get(ff1, 1),
		  freesetinequ.Get(ff1, 2),
		  freesetinequ.Get(ff1, 3));
      a2 = Vec3d (freesetinequ.Get(ff2, 1),
		  freesetinequ.Get(ff2, 2),
		  freesetinequ.Get(ff2, 3));

      if ( ( (a1 * trivec) > 0) || ( (a2 * trivec) > 0))
	return 0;

      return 1;
    }


  if (cnt == 3)
    {
      // MARK(triinfz3);  

      ARRAY<int> lpi(freezonepi.Size());
      for (i = 1; i <= lpi.Size(); i++)
	lpi.Elem(i) = 0;

      for (i = 1; i <= 3; i++)
	lpi.Elem(pi.Get(i)) = 1;
      
      for (i = 1; i <= freesetfaces.Size(); i++)
	{
	  if (lpi.Get(freesetfaces.Get(i).i1) + 
	      lpi.Get(freesetfaces.Get(i).i2) + 
	      lpi.Get(freesetfaces.Get(i).i3) == 3)
	    {
	      return 0;
	    }
	}
      return 1;
    }

  // MARK(triinfz0);  

  
  os1 = os2 = os3 = 0;
  activefaces.SetSize(0);

  // is point inside ?

  for (i = 1; i <= freesetfaces.Size(); i++)
    {
      hos1 = freesetinequ.Get(i, 1) * p1.X() +
	freesetinequ.Get(i, 2) * p1.Y() +
	freesetinequ.Get(i, 3) * p1.Z() +
	freesetinequ.Get(i, 4) > -1E-5;
      
      hos2 = freesetinequ.Get(i, 1) * p2.X() +
	freesetinequ.Get(i, 2) * p2.Y() +
	freesetinequ.Get(i, 3) * p2.Z() +
	freesetinequ.Get(i, 4) > -1E-5;
      
      hos3 = freesetinequ.Get(i, 1) * p3.X() +
	freesetinequ.Get(i, 2) * p3.Y() +
	freesetinequ.Get(i, 3) * p3.Z() +
	freesetinequ.Get(i, 4) > -1E-5;
      
      if (hos1 && hos2 && hos3) return 0;
      
      if (hos1) os1 = 1;
      if (hos2) os2 = 1;
      if (hos3) os3 = 1;
      
      if (hos1 || hos2 || hos3) activefaces.Append (i);
    }
  
  if (!os1 || !os2 || !os3) return 1;

  v1x = p2.X() - p1.X();
  v1y = p2.Y() - p1.Y();
  v1z = p2.Z() - p1.Z();

  v2x = p3.X() - p1.X();
  v2y = p3.Y() - p1.Y();
  v2z = p3.Z() - p1.Z();

  n.X() = v1y * v2z - v1z * v2y;
  n.Y() = v1z * v2x - v1x * v2z;
  n.Z() = v1x * v2y - v1y * v2x;
  n /= n.Length();

  allleft = allright = 1;
  for (i = 1; i <= transfreezone.Size() && (allleft || allright); i++)
    {
      const Point3d & p = transfreezone.Get(i);
      float scal = (p.X() - p1.X()) * n.X() +
	(p.Y() - p1.Y()) * n.Y() +
	(p.Z() - p1.Z()) * n.Z();

      if ( scal >  1E-8 ) allleft = 0;
      if ( scal < -1E-8 ) allright = 0;
    }

  if (allleft || allright) return 0;


  lam1old = lam2old = lam1 = lam2 = 1.0 / 3.0;


  //  testout << endl << endl << "Start minimizing" << endl;

  it = 0;
  int minit;
  minit = 1000;
  fold = 1E10;



  while (1)
    {
      it++;

      if (it > 1000) return -1;

      if (lam1 < 0) lam1 = 0;
      if (lam2 < 0) lam2 = 0;
      if (lam1 + lam2 > 1) lam1 = 1 - lam2;

      if (it > minit)
	{
	  (*testout) << "it = " << it << endl;
	  (*testout) << "lam1/2 = " << lam1 << "  " << lam2 << endl;
	}

      hpx = p1.X() + lam1 * v1x + lam2 * v2x;
      hpy = p1.Y() + lam1 * v1y + lam2 * v2y;
      hpz = p1.Z() + lam1 * v1z + lam2 * v2z;

      f = 0;

      h11 = h12 = h22 = dflam1 = dflam2 = 0;
      cntout = 0;

      isin = 1;

      for (i = 1; i <= activefaces.Size(); i++)
	{
	  ii = activefaces.Get(i);

	  hf = freesetinequ.Get(ii, 1) * hpx +
	    freesetinequ.Get(ii, 2) * hpy +
	    freesetinequ.Get(ii, 3) * hpz +
	    freesetinequ.Get(ii, 4);

	  if (hf > -1E-7) isin = 0;

	  hf += 1E-4;
	  if (hf > 0)
	    {
	      f += hf * hf;

	      v1n = freesetinequ.Get(ii, 1) * v1x +
		freesetinequ.Get(ii, 2) * v1y +
		freesetinequ.Get(ii, 3) * v1z;
	      v2n = freesetinequ.Get(ii, 1) * v2x +
		freesetinequ.Get(ii, 2) * v2y +
		freesetinequ.Get(ii, 3) * v2z;

	      h11 += 2 * v1n * v1n;
	      h12 += 2 * v1n * v2n;
	      h22 += 2 * v2n * v2n;
	      dflam1 += 2 * hf * v1n;
	      dflam2 += 2 * hf * v2n;
	      cntout++;
	    }
	}

      if (isin) return 1;

      if (it > minit)
	{
	  (*testout) << "f = " << f
		     << "  dfdlam = " << dflam1 << "  " << dflam2 << endl;
	  (*testout) << "h = " << h11 << "  " << h12 << "  " << h22 << endl;
	  (*testout) << "active: " << cntout << endl;
	  (*testout) << "lam1-lam1old = " << (lam1 - lam1old) << endl;
	  (*testout) << "lam2-lam2old = " << (lam2 - lam2old) << endl;
	}


      if (f >= fold)
	{
	  lam1 = 0.100000000000000 * lam1 + 0.9000000000000000 * lam1old;
	  lam2 = 0.100000000000000 * lam2 + 0.9000000000000000 * lam2old;
	}
      else
	{
	  lam1old = lam1;
	  lam2old = lam2;
	  fold = f;


	  if (f < 1E-9) return 1;

	  h11 += 1E-10;
	  h22 += 1E-10;
	  c1 = - ( h22 * dflam1 - h12 * dflam2) / (h11 * h22 - h12 * h12);
	  c2 = - (-h12 * dflam1 + h11 * dflam2) / (h11 * h22 - h12 * h12);
	  alpha = 1;


	  if (it > minit)
	    (*testout) << "c1/2 = " << c1 << "  " << c2 << endl;

	  act1 = lam1 <= 1E-6 && c1 <= 0;
	  act2 = lam2 <= 1E-6 && c2 <= 0;
	  act3 = lam1 + lam2 >= 1 - 1E-6 && c1 + c2 >= 0;

	  if (it > minit)
	    (*testout) << "act1,2,3 = " << act1 << act2 << act3 << endl;

	  if (act1 && act2 || act1 && act3 || act2 && act3) return 0;

	  if (act1)
	    {
	      c1 = 0;
	      c2 = - dflam2 / h22;
	    }

	  if (act2)
	    {
	      c1 = - dflam1 / h11;
	      c2 = 0;
	    }

	  if (act3)
	    {
	      c1 = - (dflam1 - dflam2) / (h11 + h22 - 2 * h12);
	      c2 = -c1;
	    }

	  if (it > minit)
	    (*testout) << "c1/2 now = " << c1 << "  " << c2 << endl;


	  if (f > 100 * sqrt (sqr (c1) + sqr (c2))) return 0;


	  if (lam1 + alpha * c1 < 0 && !act1)
	    alpha = -lam1 / c1;
	  if (lam2 + alpha * c2 < 0 && !act2)
	    alpha = -lam2 / c2;
	  if (lam1 + lam2 + alpha * (c1 + c2) > 1 && !act3)
	    alpha = (1 - lam1 - lam2) / (c1 + c2);

	  if (it > minit)
	    (*testout) << "alpha = " << alpha << endl;

	  lam1 += alpha * c1;
	  lam2 += alpha * c2;
	}
    }
}




int vnetrule :: IsQuadInFreeZone (const Point3d & p1, 
				  const Point3d & p2,
				  const Point3d & p3, 
				  const Point3d & p4, 
				  const ARRAY<int> & pi, int newone)
{
  int fs;
  int infreeset, cannot = 0;


  static ARRAY<int> pfi(4), pfi2(4);

  // convert from local index to freeset index
  int i, j;
  for (i = 1; i <= 4; i++)
    {
      pfi.Elem(i) = 0;
      if (pi.Get(i))
	{
	  for (j = 1; j <= freezonepi.Size(); j++)
	    if (freezonepi.Get(j) == pi.Get(i))
	      pfi.Elem(i) = j;
	}
    }

  for (fs = 1; fs <= freesets.Size(); fs++)
    {
      const ARRAY<int> & freeseti = *freesets.Get(fs);
      for (i = 1; i <= 4; i++)
	{
	  pfi2.Elem(i) = 0;
	  for (j = 1; j <= freeseti.Size(); j++)
	    if (pfi.Get(i) == freeseti.Get(j))
	      pfi2.Elem(i) = pfi.Get(i);
	}

      infreeset = IsQuadInFreeSet(p1, p2, p3, p4, fs, pfi2, newone);
      if (infreeset == 1) return 1;
      if (infreeset == -1) cannot = -1;
    }
  
  return cannot;
}


int vnetrule :: IsQuadInFreeSet (const Point3d & p1, const Point3d & p2,
				 const Point3d & p3, const Point3d & p4, 
				 int fs, const ARRAY<int> & pi, int newone)
{
  int i;
  
  int cnt = 0;
  for (i = 1; i <= 4; i++)
    if (pi.Get(i)) cnt++;
  
  /*
  (*testout) << "test quad in freeset: " << p1 << " - " << p2 << " - " << p3 << " - " << p4 << endl;
  (*testout) << "pi = ";
  for (i = 1; i <= pi.Size(); i++)
    (*testout) << pi.Get(i) << " ";
  (*testout) << endl;
  (*testout) << "cnt = " << cnt  << endl;
  */
  if (cnt == 4)
    {
      return 1;
    }

  if (cnt == 3)
    {
      return 1;
    }

  static ARRAY<int> pi3(3);
  int res;

  pi3.Elem(1) = pi.Get(1);
  pi3.Elem(2) = pi.Get(2);
  pi3.Elem(3) = pi.Get(3);
  res = IsTriangleInFreeSet (p1, p2, p3, fs, pi3, newone);
  if (res) return res;


  pi3.Elem(1) = pi.Get(2);
  pi3.Elem(2) = pi.Get(3);
  pi3.Elem(3) = pi.Get(4);
  res = IsTriangleInFreeSet (p2, p3, p4, fs, pi3, newone);
  if (res) return res;

  pi3.Elem(1) = pi.Get(3);
  pi3.Elem(2) = pi.Get(4);
  pi3.Elem(3) = pi.Get(1);
  res = IsTriangleInFreeSet (p3, p4, p1, fs, pi3, newone);
  if (res) return res;

  pi3.Elem(1) = pi.Get(4);
  pi3.Elem(2) = pi.Get(1);
  pi3.Elem(3) = pi.Get(2);
  res = IsTriangleInFreeSet (p4, p1, p2, fs, pi3, newone);
  return res;
}












float vnetrule :: CalcPointDist (int pi, const Point3d & p) const
{
  float dx = p.X() - points.Get(pi).X();
  float dy = p.Y() - points.Get(pi).Y();
  float dz = p.Z() - points.Get(pi).Z();
  
  //  const threefloat * tf = &tolerances.Get(pi);
  //  return tf->f1 * dx * dx + tf->f2 * dx * dy + tf->f3 * dy * dy;
  return tolerances.Get(pi) * (dx * dx + dy * dy + dz * dz);
}


int vnetrule :: TestOk () const
{
  ARRAY<int> cntpused(points.Size());
  ARRAY<int> edge1, edge2;
  ARRAY<int> delf(faces.Size());
  int i, j, k;
  int pi1, pi2;
  int found;

  for (i = 1; i <= cntpused.Size(); i++)
    cntpused.Elem(i) = 0;
  for (i = 1; i <= faces.Size(); i++)
    delf.Elem(i) = 0;
  for (i = 1; i <= delfaces.Size(); i++)
    delf.Elem(delfaces.Get(i)) = 1;


  for (i = 1; i <= faces.Size(); i++)
    if (delf.Get(i) || i > noldf)
      for (j = 1; j <= faces.Get(i).GetNP(); j++)
        cntpused.Elem(faces.Get(i).PNum(j))++;

  for (i = 1; i <= cntpused.Size(); i++)
    if (cntpused.Get(i) > 0 && cntpused.Get(i) < 2)
      {
	return 0;
      }


  //  (*testout) << endl;
  for (i = 1; i <= faces.Size(); i++)
    {
      //      (*testout) << "face " << i << endl;
      for (j = 1; j <= faces.Get(i).GetNP(); j++)
	{
	  pi1 = 0; pi2 = 0;
	  if (delf.Get(i))
	    {
	      pi1 = faces.Get(i).PNumMod(j);
	      pi2 = faces.Get(i).PNumMod(j+1);
	    }
	  if (i > noldf)
	    {
	      pi1 = faces.Get(i).PNumMod(j+1);
	      pi2 = faces.Get(i).PNumMod(j);
	    }

	  found = 0;
	  if (pi1)
	    {
	      for (k = 1; k <= edge1.Size(); k++)
		if (edge1.Get(k) == pi1 && edge2.Get(k) == pi2)
		  {
		    found = 1;
		    edge1.DeleteElement(k);
		    edge2.DeleteElement(k);
		    k--;
		    //		    (*testout) << "Del edge " << pi1 << "-" << pi2 << endl;
		  }
	      if (!found)
		{
		  edge1.Append (pi2);
		  edge2.Append (pi1);
		  //		  (*testout) << "Add edge " << pi1 << "-" << pi2 << endl;
		}
	    }
	}
    }


  if (edge1.Size() > 0)
    {
      return 0;
    }

  /*
    cntpused.SetSize(freezone.Size());
    for (i = 1; i <= cntpused.Size(); i++)
    cntpused[i] = 0;

    for (i = 1; i <= freefaces.Size(); i++)
    {
    cntpused[freefaces[i].i1]++;
    cntpused[freefaces[i].i2]++;
    cntpused[freefaces[i].i3]++;
    }

    for (i = 1; i <= cntpused.Size(); i++)
    if (cntpused[i] < 3)
    {
    (*mycout) << "Fall 3" << endl;
    return 0;
    }



    for (i = 1; i <= freefaces.Size(); i++)
    {
    for (j = 1; j <= 3; j++)
    {
    if (j == 1)
    {
    pi1 = freefaces[i].i1;
    pi2 = freefaces[i].i2;
    }
    if (j == 2)
    {
    pi1 = freefaces[i].i2;
    pi2 = freefaces[i].i3;
    }
    if (j == 3)
    {
    pi1 = freefaces[i].i3;
    pi2 = freefaces[i].i1;
    }

    found = 0;
    for (k = 1; k <= edge1.Size(); k++)
    if (edge1[k] == pi1 && edge2[k] == pi2)
    {
    found = 1;
    edge1.DeleteElement(k);
    edge2.DeleteElement(k);
    k--;
    }

    if (!found)
    {
    edge1.Append (pi2);
    edge2.Append (pi1);
    }
    }
    }

    if (edge1.Size() > 0)
    {
    (*mycout) << "Fall 4" << endl;
    return 0;
    }
    */
  return 1;
}


int vnetrule :: IsDelFace (int fn) const
{
  int i;
  for (i = 1; i <= GetNDelF(); i++)
    if (GetDelFace(i) == fn) return 1;
  return 0;
}

}
