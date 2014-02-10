// find inner point

template <typename POINTARRAY, typename FACEARRAY>
inline int FindInnerPoint2 (POINTARRAY & points,
			    FACEARRAY & faces,
			    Point3d & p)
{
  static int timer = NgProfiler::CreateTimer ("FindInnerPoint2");
  NgProfiler::RegionTimer reg (timer);

  ARRAY<Vec3d> a;
  ARRAY<double> c;
  Mat<3> m, inv;
  Vec<3> rs, x, pmin;

  int nf = faces.Size();

  a.SetSize (nf);
  c.SetSize (nf);

  for (int i = 0; i < nf; i++)
    {
      Point3d p1 = points.Get(faces[i][0]);
      a[i] = Cross (points.Get(faces[i][1]) - p1,
		    points.Get(faces[i][2]) - p1);
      a[i] /= a[i].Length();
      c[i] = - (a[i].X() * p1.X() + a[i].Y() * p1.Y() + a[i].Z() * p1.Z());
    }


  x = 0;
  
  
  double hmax = 0;
  for (int i = 0; i < nf; i++)
    {
      const Element2d & el = faces[i];
      for (int j = 1; j <= 3; j++)
	{
	  double hi = Dist (points.Get(el.PNumMod(j)),
			    points.Get(el.PNumMod(j+1)));
	  if (hi > hmax) hmax = hi;
	}
    }

  double fmin = 0;

  for (int i1 = 1; i1 <= nf; i1++)
    for (int i2 = i1+1; i2 <= nf; i2++)
      for (int i3 = i2+1; i3 <= nf; i3++)
        for (int i4 = i3+1; i4 <= nf; i4++)
          {
	    m(0, 0) = a.Get(i1).X() - a.Get(i2).X();
	    m(0, 1) = a.Get(i1).Y() - a.Get(i2).Y();
	    m(0, 2) = a.Get(i1).Z() - a.Get(i2).Z();
	    rs(0) = c.Get(i2) - c.Get(i1);

	    m(1, 0) = a.Get(i1).X() - a.Get(i3).X();
	    m(1, 1) = a.Get(i1).Y() - a.Get(i3).Y();
	    m(1, 2) = a.Get(i1).Z() - a.Get(i3).Z();
	    rs(1) = c.Get(i3) - c.Get(i1);

	    m(2, 0) = a.Get(i1).X() - a.Get(i4).X();
	    m(2, 1) = a.Get(i1).Y() - a.Get(i4).Y();
	    m(2, 2) = a.Get(i1).Z() - a.Get(i4).Z();
	    rs(2) = c.Get(i4) - c.Get(i1);


	    if (fabs (Det (m)) > 1e-10)
	      {
		CalcInverse (m, inv);
		x = inv * rs;

		double f = -1e10;
		for (int i = 0; i < nf; i++)
		  {
		    double hd = 
		      x(0) * a[i].X() + x(1) * a[i].Y() + x(2) * a[i].Z() + c[i];
		    if (hd > f) f = hd;
		    if (hd > fmin) break;
		  }

		if (f < fmin)
		  {
		    fmin = f;
		    pmin = x;
		  }
	      }
          }

  p = Point3d (pmin(0), pmin(1), pmin(2));
  (*testout) << "fmin = " << fmin << endl;
  return (fmin < -1e-3 * hmax);
}

