// find inner point



inline void Minimize (const ARRAY<Vec3d> & a,
		      const ARRAY<double> & c,
		      int * act, 
		      Vec<3> & x, double & f,
		      int * sol)
{
  int act1[4];
  Mat<3> m, inv;
  Vec<3> rs, xmax, center;

  f = 1e99;

  for (int j = 0; j < 5; j++)
    {
      for (int hk = 0, k = 0; hk < 4; hk++)
	{
	  if (hk == j) k++;
	  act1[hk] = act[k];
	  k++;
	}

      for (int k = 0; k < 3; k++)
	{
	  m(k, 0) = a[act1[0]].X() - a[act1[k+1]].X();
	  m(k, 1) = a[act1[0]].Y() - a[act1[k+1]].Y();
	  m(k, 2) = a[act1[0]].Z() - a[act1[k+1]].Z();
	  rs(k) = c[act1[k+1]] - c[act1[0]];
	}

      /*
      (*testout) << "act1 = "
		 << act1[0] << " "
		 << act1[1] << " "
		 << act1[2] << " "
		 << act1[3] << endl;
      (*testout) << "Det = " << Det(m) << endl;
      */

      if (fabs (Det (m)) > 1e-10)
	{
	  CalcInverse (m, inv);
	  xmax = inv * rs;
	  
	  double fmax = -1e10;
	  for (int k = 0; k < 5; k++)
	    {
	      double hd = 
		xmax(0) * a[act[k]].X() + xmax(1) * a[act[k]].Y() + xmax(2) * a[act[k]].Z() + c[act[k]];
	      if (hd > fmax) fmax = hd;
	    }

	  if (fmax < f)
	    {
	      f = fmax;
	      x = xmax;
	      for (int k = 0; k < 4; k++)
		sol[k] = act1[k];
	    }
	}
    }
}




template <typename POINTARRAY, typename FACEARRAY>
inline int FindInnerPoint (POINTARRAY & points,
			   FACEARRAY & faces,
			   Point3d & p)
{
  static int timer = NgProfiler::CreateTimer ("FindInnerPoint");
  NgProfiler::RegionTimer reg (timer);

  ARRAY<Vec3d> a;
  ARRAY<double> c;
  Mat<3> m, inv;
  Vec<3> rs, x, center;
  double f;

  int nf = faces.Size();

  // minimize_x  max_i  a_i x + c_i

  a.SetSize (nf+4);
  c.SetSize (nf+4);

  for (int i = 0; i < nf; i++)
    {
      Point3d p1 = points.Get(faces[i][0]);
      a[i] = Cross (points.Get(faces[i][1]) - p1,
		    points.Get(faces[i][2]) - p1);
      a[i] /= a[i].Length();
      c[i] = - (a[i].X() * p1.X() + a[i].Y() * p1.Y() + a[i].Z() * p1.Z());
    }

  /*
  center = 0;
  for (int i = 0; i < points.Size(); i++)
    center += Vec<3> (points[i]);
  center /= points.Size();
  */

  center = 0;
  for (int i = 0; i < faces.Size(); i++)
    for (int j = 0; j < 3; j++)
      center += Vec<3> (points.Get(faces[i][j]));
  center /= (3*faces.Size());


  // (*testout) << "center = " << center << endl;

  double hmax = 0;
  for (int i = 0; i < nf; i++)
    {
      // const Element2d & el = faces[i];
      // (*testout) << "el[" << i << "] = " << el << endl;
      for (int j = 1; j <= 3; j++)
	{
	  double hi = Dist (points.Get(faces[i].PNumMod(j)),
			    points.Get(faces[i].PNumMod(j+1)));
	  if (hi > hmax) hmax = hi;
	}
    }
  
  // (*testout) << "hmax = " << hmax << endl;
  
  a[nf] = Vec<3> (1, 0, 0);
  c[nf] = -center(0) - hmax;
  a[nf+1] = Vec<3> (0, 1, 0);
  c[nf+1] = -center(1) - hmax;
  a[nf+2] = Vec<3> (0, 0, 1);
  c[nf+2] = -center(2) - hmax;
  a[nf+3] = Vec<3> (-1, -1, -1);
  c[nf+3] = center(0)+center(1)+center(2)-3*hmax;

  /*
  (*testout) << "findip, a now = " << endl << a << endl;
  (*testout) << "findip, c now = " << endl << c << endl;
  */

  int act[5] = { 0, nf, nf+1, nf+2, nf+3 };
  int sol[4];

  while (1)
    {
      /*
      (*testout) << "try ";
      for (int j = 0; j < 5; j++)
	(*testout)  << act[j] << " ";
      */

      Minimize (a, c, act, x, f, sol);

      /*
      (*testout) << endl << "sol = ";
      for (int j = 0; j < 4; j++)
	(*testout)  << sol[j] << " ";

      (*testout) << " fmin = " << f << endl;
      */
      for (int j = 0; j < 4; j++) act[j] = sol[j];
      
      bool found = 0;
      double maxval = f;
      for (int j = 0; j < nf; j++)
	{
	  double val = x(0) * a[j].X() + x(1) * a[j].Y() + x(2) * a[j].Z() + c[j];
	  if (val > maxval + hmax * 1e-6)
	    {
	      found = 1;
	      maxval = val;
	      act[4] = j;
	    }
	}
      
      // (*testout) << "maxval = " << maxval << endl;
      if (!found) break;
    }
  
  // cout << "converged, f = " << f << endl;
  
  p = Point3d (x(0), x(1), x(2));
  // (*testout) << "findip, f = " << f << ", hmax = " << hmax << endl;
  return (f < -1e-5 * hmax);
}



