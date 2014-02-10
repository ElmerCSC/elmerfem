#include <algorithm>
#include <mystdlib.h>

#include <myadt.hpp>
#include <gprim.hpp>

namespace netgen
{
ostream & operator<<(ostream  & s, const Point3d & p)
  {
  return s << "(" << p.x[0] << ", " << p.x[1] << ", " << p.x[2] << ")";
  }

ostream & operator<<(ostream  & s, const Vec3d & v)
  {
  return s << "(" << v.x[0] << ", " << v.x[1] << ", " << v.x[2] << ")";
  }

double Angle (const Vec3d & v1, const Vec3d & v2)
{
  double co = (v1 * v2) / (v1.Length() * v2.Length());
  if (co > 1) co = 1;
  if (co < -1) co = -1;
  return acos ( co );
}


void Vec3d :: GetNormal (Vec3d & n) const
  {
  if (fabs (X()) > fabs (Z()))
    {
    n.X() = -Y();
    n.Y() = X();
    n.Z() = 0;
    }
  else
    {
    n.X() = 0;
    n.Y() = Z();
    n.Z() = -Y();
    }
  double len = n.Length();
  if (len == 0)
    {
    n.X() = 1;
    n.Y() = n.Z() = 0;
    }
  else
    n /= len;
  }

/*
ostream & operator<<(ostream  & s, const ROTDenseMatrix3D & r)
  {
  return s << "{ (" << r.txx << ", " << r.txy << ", " << r.txz << ") , ("
                    << r.tyx << ", " << r.tyy << ", " << r.tyz << ") , ("
                    << r.tzx << ", " << r.tzy << ", " << r.tzz << ") }";
  }
*/

/*
Vec3d operator- (const Point3d & p1, const Point3d & p2)
  {
  return Vec3d (p1.X() - p2.X(), p1.Y() - p2.Y(),p1.Z() - p2.Z());
  }

Point3d operator- (const Point3d & p1, const Vec3d & v)
  {
  return Point3d (p1.X() - v.X(), p1.Y() - v.Y(),p1.Z() - v.Z());
  }

Point3d operator+ (const Point3d & p1, const Vec3d & v)
  {
  return Point3d (p1.X() + v.X(), p1.Y() + v.Y(),p1.Z() + v.Z());
  }

Vec3d operator- (const Vec3d & v1, const Vec3d & v2)
  {
  return Vec3d (v1.X() - v2.X(), v1.Y() - v2.Y(),v1.Z() - v2.Z());
  }

Vec3d operator+ (const Vec3d & v1, const Vec3d & v2)
  {
  return Vec3d (v1.X() + v2.X(), v1.Y() + v2.Y(),v1.Z() + v2.Z());
  }

Vec3d operator* (double scal, const Vec3d & v)
  {
  return Vec3d (scal * v.X(), scal * v.Y(), scal * v.Z());
  }
*/
/*
double operator* (const Vec3d & v1, const Vec3d & v2)
  {
  return v1.X() * v2.X() + v1.Y() * v2.Y() + v1.Z() * v2.Z();
  }

double Cross (const Vec3d & v1, const Vec3d & v2)
  {
  return v1.X() * v2.Y() - v1.Y() * v2.X();
  }
*/

/*
void ROTDenseMatrix3D :: CalcRotMat(double ag, double bg, double lg, double size2, Vec3d r)
  {
  size = size2;
  txx=size * ( cos(bg) * cos(lg) );
  txy=size * ( cos(bg) * sin(lg) );
  txz=size * (-sin(bg)           );

  tyx=size * ( sin(ag) * sin(bg) * cos(lg) - cos(ag) * sin(lg) );
  tyy=size * ( sin(ag) * sin(bg) * sin(lg) + cos(ag) * cos(lg) );
  tyz=size * ( sin(ag) * cos(bg)                               );

  tzx=size * ( cos(ag) * sin(bg) * cos(lg) + sin(ag) * sin(lg) );
  tzy=size * ( cos(ag) * sin(bg) * sin(lg) - sin(ag) * cos(lg) );
  tzz=size * ( cos(ag) * cos(bg)                               );

  deltaR=r;
  }
ROTDenseMatrix3D :: ROTDenseMatrix3D(double ag, double bg, double lg, double size2, Vec3d r)
  {CalcRotMat(ag, bg, lg, size2, r); }

ROTDenseMatrix3D :: ROTDenseMatrix3D(Vec3d rot2)
  {
  Vec3d r2(0,0,0);
  CalcRotMat(rot2.X(), rot2.Y(), rot2.Z(), 1, r2);
  }

ROTDenseMatrix3D ROTDenseMatrix3D :: INV()
  {
  ROTDenseMatrix3D rinv(txx/sqr(size),tyx/sqr(size),tzx/sqr(size),
                   txy/sqr(size),tyy/sqr(size),tzy/sqr(size),
                   txz/sqr(size),tyz/sqr(size),tzz/sqr(size),
                   1/size,deltaR);
  return rinv;
  }

Vec3d operator* (const ROTDenseMatrix3D & r, const Vec3d & v)
  {
  return Vec3d (r.XX() * v.X() + r.XY() * v.Y() + r.XZ() * v.Z(),
                r.YX() * v.X() + r.YY() * v.Y() + r.YZ() * v.Z(),
                r.ZX() * v.X() + r.ZY() * v.Y() + r.ZZ() * v.Z() );
  }

Point3d operator* (const ROTDenseMatrix3D & r, const Point3d & p)
  {
  return Point3d (r.XX() * p.X() + r.XY() * p.Y() + r.XZ() * p.Z(),
                  r.YX() * p.X() + r.YY() * p.Y() + r.YZ() * p.Z(),
                  r.ZX() * p.X() + r.ZY() * p.Y() + r.ZZ() * p.Z() );
  }
*/







Box3d :: Box3d ( double aminx, double amaxx,
                 double aminy, double amaxy,
                 double aminz, double amaxz )
{
  minx[0] = aminx; maxx[0] = amaxx;
  minx[1] = aminy; maxx[1] = amaxy;
  minx[2] = aminz; maxx[2] = amaxz;
}

Box3d :: Box3d ( const Box3d & b2 )
{
  for (int i = 0; i < 3; i++)
    {
      minx[i] = b2.minx[i];
      maxx[i] = b2.maxx[i];
    }
}

Box3d :: Box3d ( const Box<3> & b2 )
{
  for (int i = 0; i < 3; i++)
    {
      minx[i] = b2.PMin()(i);
      maxx[i] = b2.PMax()(i);
    }
}


/*
int Box3d :: Intersect (const Box3d & box2) const
{
  int i;
  for (i = 0; i <= 2; i++)
    if (minx[i] > box2.maxx[i] || maxx[i] < box2.minx[i])
      return 0;
  return 1;
}
*/

/*
void Box3d :: SetPoint (const Point3d & p)
{
  minx[0] = maxx[0] = p.X();
  minx[1] = maxx[1] = p.Y();
  minx[2] = maxx[2] = p.Z();
}

void Box3d :: AddPoint (const Point3d & p)
{
  if (p.X() < minx[0]) minx[0] = p.X();
  if (p.X() > maxx[0]) maxx[0] = p.X();
  if (p.Y() < minx[1]) minx[1] = p.Y();
  if (p.Y() > maxx[1]) maxx[1] = p.Y();
  if (p.Z() < minx[2]) minx[2] = p.Z();
  if (p.Z() > maxx[2]) maxx[2] = p.Z();
}
*/

void Box3d :: GetPointNr (int i, Point3d & point) const
{
  i--;
  point.X() = (i & 1) ? maxx[0] : minx[0];
  point.Y() = (i & 2) ? maxx[1] : minx[1];
  point.Z() = (i & 4) ? maxx[2] : minx[2];
}


void Box3d :: Increase (double d)
{
  for (int i = 0; i <= 2; i++)
    {
      minx[i] -= d;
      maxx[i] += d;
    }
}

void Box3d :: IncreaseRel (double /* rel */)
{
  for (int i = 0; i <= 2; i++)
    {
      double d = 0.5 * (maxx[i] - minx[i]);
      minx[i] -= d;
      maxx[i] += d;
    }
}


Box3d :: Box3d (const Point3d& p1, const Point3d& p2)
{
  minx[0] = min2 (p1.X(), p2.X());
  minx[1] = min2 (p1.Y(), p2.Y());
  minx[2] = min2 (p1.Z(), p2.Z());
  maxx[0] = max2 (p1.X(), p2.X());
  maxx[1] = max2 (p1.Y(), p2.Y());
  maxx[2] = max2 (p1.Z(), p2.Z());
}

const Box3d& Box3d :: operator+=(const Box3d& b)
{
  minx[0] = min2 (minx[0], b.minx[0]);
  minx[1] = min2 (minx[1], b.minx[1]);
  minx[2] = min2 (minx[2], b.minx[2]);
  maxx[0] = max2 (maxx[0], b.maxx[0]);
  maxx[1] = max2 (maxx[1], b.maxx[1]);
  maxx[2] = max2 (maxx[2], b.maxx[2]);

  return *this;
}

Point3d Box3d :: MaxCoords() const
{
  return Point3d(maxx[0], maxx[1], maxx[2]);
}

Point3d Box3d :: MinCoords() const
{
  return Point3d(minx[0], minx[1], minx[2]);
}

/*
void Box3d :: CreateNegMinMaxBox()
{
  minx[0] = MAXDOUBLE;
  minx[1] = MAXDOUBLE;
  minx[2] = MAXDOUBLE;
  maxx[0] = MINDOUBLE;
  maxx[1] = MINDOUBLE;
  maxx[2] = MINDOUBLE;

}
*/

void Box3d :: WriteData(ofstream& fout) const
{
  for(int i = 0; i < 3; i++)
    {
      fout << minx[i] << " " << maxx[i] << " ";
    }
  fout << "\n";
}

void Box3d :: ReadData(ifstream& fin)
{
  for(int i = 0; i < 3; i++)
    {
      fin >> minx[i];
      fin >> maxx[i];
    }
}




Box3dSphere :: Box3dSphere ( double aminx, double amaxx,
			     double aminy, double amaxy,
			     double aminz, double amaxz )
  : Box3d (aminx, amaxx, aminy, amaxy, aminz, amaxz)
{
  CalcDiamCenter ();
}


void Box3dSphere :: CalcDiamCenter ()
{
  diam = sqrt( sqr (maxx[0] - minx[0]) +
	       sqr (maxx[1] - minx[1]) + 
	       sqr (maxx[2] - minx[2]));
  
  c.X() = 0.5 * (minx[0] + maxx[0]);
  c.Y() = 0.5 * (minx[1] + maxx[1]);
  c.Z() = 0.5 * (minx[2] + maxx[2]);
  
  inner = min2 ( min2 (maxx[0] - minx[0], maxx[1] - minx[1]), maxx[2] - minx[2]) / 2;
}


void Box3dSphere :: GetSubBox (int i, Box3dSphere & sbox) const
{
  i--;
  if (i & 1)
    {
      sbox.minx[0] = c.X();
      sbox.maxx[0] = maxx[0];
    }
  else
    {
      sbox.minx[0] = minx[0];
      sbox.maxx[0] = c.X();
    }
  if (i & 2)
    {
      sbox.minx[1] = c.Y();
      sbox.maxx[1] = maxx[1];
    }
  else
    {
      sbox.minx[1] = minx[1];
      sbox.maxx[1] = c.Y();
    }
  if (i & 4)
    {
      sbox.minx[2] = c.Z();
      sbox.maxx[2] = maxx[2];
    }
  else
    {
      sbox.minx[2] = minx[2];
      sbox.maxx[2] = c.Z();
    }
  
  //  sbox.CalcDiamCenter ();

  sbox.c.X() = 0.5 * (sbox.minx[0] + sbox.maxx[0]);
  sbox.c.Y() = 0.5 * (sbox.minx[1] + sbox.maxx[1]);
  sbox.c.Z() = 0.5 * (sbox.minx[2] + sbox.maxx[2]);
  sbox.diam = 0.5 * diam;
  sbox.inner = 0.5 * inner;
}




/*
double Determinant (const Vec3d & col1,
		    const Vec3d & col2,
		    const Vec3d & col3)
{
  return
    col1.x[0] * ( col2.x[1] * col3.x[2] - col2.x[2] * col3.x[1]) +
    col1.x[1] * ( col2.x[2] * col3.x[0] - col2.x[0] * col3.x[2]) +
    col1.x[2] * ( col2.x[0] * col3.x[1] - col2.x[1] * col3.x[0]);
}
*/

void Transpose (Vec3d & v1, Vec3d & v2, Vec3d & v3)
{
  Swap (v1.Y(), v2.X());
  Swap (v1.Z(), v3.X());
  Swap (v2.Z(), v3.Y());
}


int SolveLinearSystem (const Vec3d & col1, const Vec3d & col2,
		       const Vec3d & col3, const Vec3d & rhs,
		       Vec3d & sol)
{
  // changed by MW
  double matrix[3][3];
  double locrhs[3];
  int retval = 0;

  for(int i=0; i<3; i++)
    {
      matrix[i][0] = col1.X(i+1);
      matrix[i][1] = col2.X(i+1);
      matrix[i][2] = col3.X(i+1);
      locrhs[i] = rhs.X(i+1);
    }

  for(int i=0; i<2; i++)
    {
      int pivot = i;
      double maxv = fabs(matrix[i][i]);
      for(int j=i+1; j<3; j++)
	if(fabs(matrix[j][i]) > maxv)
	  {
	    maxv = fabs(matrix[j][i]);
	    pivot = j;
	  }

      if(fabs(maxv) > 1e-40)
	{
	  if(pivot != i)
	    {
	      swap(matrix[i][0],matrix[pivot][0]);
	      swap(matrix[i][1],matrix[pivot][1]);
	      swap(matrix[i][2],matrix[pivot][2]);
	      swap(locrhs[i],locrhs[pivot]);
	    }
	  for(int j=i+1; j<3; j++)
	    {
	      double fac = matrix[j][i] / matrix[i][i];
	      
	      for(int k=i+1; k<3; k++)
		matrix[j][k] -= fac*matrix[i][k];
	      locrhs[j] -= fac*locrhs[i];
	    }
	}
      else
	retval = 1;
    }

  if(fabs(matrix[2][2]) < 1e-40)
    retval = 1;

  if(retval != 0)
    return retval;
  

  for(int i=2; i>=0; i--)
    {
      double sum = locrhs[i];
      for(int j=2; j>i; j--)
	sum -= matrix[i][j]*sol.X(j+1);

      sol.X(i+1) = sum/matrix[i][i];
    }

  return 0;
  
  
  


  /*
  double det = Determinant (col1, col2, col3);
  if (fabs (det) < 1e-40)
    return 1;
  
  sol.X() = Determinant (rhs, col2, col3) / det;
  sol.Y() = Determinant (col1, rhs, col3) / det;
  sol.Z() = Determinant (col1, col2, rhs) / det;

  return 0;
  */
  /*
  Vec3d cr;
  Cross (col1, col2, cr);
  double det = cr * col3;

  if (fabs (det) < 1e-40)
    return 1;

  if (fabs(cr.Z()) > 1e-12)
    {
      // solve for 3. component
      sol.Z() = (cr * rhs) / det;
      
      // 2x2 system for 1. and 2. component
      double res1 = rhs.X() - sol.Z() * col3.X();
      double res2 = rhs.Y() - sol.Z() * col3.Y();
      
      sol.X() = (col2.Y() * res1 - col2.X() * res2) / cr.Z();
      sol.Y() = (col1.X() * res2 - col1.Y() * res1) / cr.Z();
  
    }
  else
    {
      det = Determinant (col1, col2, col3);
      if (fabs (det) < 1e-40)
	return 1;
      
      sol.X() = Determinant (rhs, col2, col3) / det;
      sol.Y() = Determinant (col1, rhs, col3) / det;
      sol.Z() = Determinant (col1, col2, rhs) / det;
    }

  return 0;
  */
}


int SolveLinearSystemLS (const Vec3d & col1,
			 const Vec3d & col2,
			 const Vec2d & rhs,
			 Vec3d & sol)
{
  double a11 = col1 * col1;
  double a12 = col1 * col2;
  double a22 = col2 * col2;
  
  double det = a11 * a22 - a12 * a12;

  if (det*det <= 1e-24 * a11 * a22)
    {
      sol = Vec3d (0, 0, 0);
      return 1;
    }
  
  Vec2d invrhs;
  invrhs.X() = ( a22 * rhs.X() - a12 * rhs.Y()) / det;
  invrhs.Y() = (-a12 * rhs.X() + a11 * rhs.Y()) / det;

  sol.X() = invrhs.X() * col1.X() + invrhs.Y() * col2.X();
  sol.Y() = invrhs.X() * col1.Y() + invrhs.Y() * col2.Y();
  sol.Z() = invrhs.X() * col1.Z() + invrhs.Y() * col2.Z();

  return 0;

  /*
  Vec3d inv1, inv2;
  int err = 
    PseudoInverse (col1, col2, inv1, inv2);

   sol = rhs.X() * inv1 + rhs.Y() * inv2;
   return err;
  */
}

int SolveLinearSystemLS2 (const Vec3d & col1,
			 const Vec3d & col2,
			 const Vec2d & rhs,
			 Vec3d & sol, double & x, double & y)
{
  double a11 = col1 * col1;
  double a12 = col1 * col2;
  double a22 = col2 * col2;
  
  double det = a11 * a22 - a12 * a12;

  if (fabs (det) <= 1e-12 * col1.Length() * col2.Length() || 
      col1.Length2() == 0 || col2.Length2() == 0)
    {
      sol = Vec3d (0, 0, 0);
      x = 0; y = 0;
      return 1;
    }
  
  Vec2d invrhs;
  invrhs.X() = ( a22 * rhs.X() - a12 * rhs.Y()) / det;
  invrhs.Y() = (-a12 * rhs.X() + a11 * rhs.Y()) / det;

  sol.X() = invrhs.X() * col1.X() + invrhs.Y() * col2.X();
  sol.Y() = invrhs.X() * col1.Y() + invrhs.Y() * col2.Y();
  sol.Z() = invrhs.X() * col1.Z() + invrhs.Y() * col2.Z();

  x = invrhs.X();
  y = invrhs.Y();

  return 0;

  /*
  Vec3d inv1, inv2;
  int err = 
    PseudoInverse (col1, col2, inv1, inv2);

   sol = rhs.X() * inv1 + rhs.Y() * inv2;
   return err;
  */
}

int PseudoInverse (const Vec3d & col1,
		   const Vec3d & col2,
		   Vec3d & inv1,
		   Vec3d & inv2)
{
  double a11 = col1 * col1;
  double a12 = col1 * col2;
  double a22 = col2 * col2;
  
  double det = a11 * a22 - a12 * a12;

  if (fabs (det) < 1e-12 * col1.Length() * col2.Length())
    {
      inv1 = Vec3d (0, 0, 0);
      inv2 = Vec3d (0, 0, 0);
      return 1;
    }

  double ia11 = a22 / det;
  double ia12 = -a12 / det;
  double ia22 = a11 / det;

  inv1 = ia11 * col1 + ia12 * col2;
  inv2 = ia12 * col1 + ia22 * col2;

  return 0;
}




QuadraticFunction3d :: 
QuadraticFunction3d (const Point3d & p, const Vec3d & v)
{
  Vec3d hv(v);
  hv /= (hv.Length() + 1e-12);
  Vec3d t1, t2;
  hv.GetNormal (t1);
  Cross (hv, t1, t2);
  
  double t1p = t1.X() * p.X() + t1.Y() * p.Y() + t1.Z() * p.Z();
  double t2p = t2.X() * p.X() + t2.Y() * p.Y() + t2.Z() * p.Z();
  c0 = sqr (t1p) + sqr (t2p);
  cx = -2 * (t1p * t1.X() + t2p * t2.X());
  cy = -2 * (t1p * t1.Y() + t2p * t2.Y());
  cz = -2 * (t1p * t1.Z() + t2p * t2.Z());

  cxx = t1.X() * t1.X() + t2.X() * t2.X();
  cyy = t1.Y() * t1.Y() + t2.Y() * t2.Y();
  czz = t1.Z() * t1.Z() + t2.Z() * t2.Z();

  cxy = 2 * t1.X() * t1.Y() + 2 * t2.X() * t2.Y();
  cxz = 2 * t1.X() * t1.Z() + 2 * t2.X() * t2.Z();
  cyz = 2 * t1.Y() * t1.Z() + 2 * t2.Y() * t2.Z();

  /*
  (*testout) << "c0 = " << c0
	     << " clin = " << cx << " " << cy << " " << cz 
	     << " cq = " << cxx << " " << cyy << " " << czz
	     << cxy << " " << cyz << " " << cyz << endl;
  */
}

// QuadraticFunction3d gqf (Point3d (0,0,0), Vec3d (1, 0, 0));





void referencetransform :: Set (const Point3d & p1, const Point3d & p2,
                                const Point3d & p3, double ah)
{
  ex = p2 - p1;
  ex /= ex.Length();
  ey = p3 - p1;
  ey -= (ex * ey) * ex;
  ey /= ey.Length();
  ez = Cross (ex, ey);
  rp = p1;
  h = ah;

  exh = ah * ex;
  eyh = ah * ey;
  ezh = ah * ez;
  ah = 1 / ah;
  ex_h = ah * ex;
  ey_h = ah * ey;
  ez_h = ah * ez;
}

void referencetransform :: ToPlain (const Point3d & p, Point3d & pp) const
{
  Vec3d v;
  v = p - rp;
  pp.X() = (ex_h * v);
  pp.Y() = (ey_h * v);
  pp.Z() = (ez_h * v);
}

void referencetransform :: ToPlain (const ARRAY<Point3d> & p,
                                    ARRAY<Point3d> & pp) const
{
  Vec3d v;
  int i;

  pp.SetSize (p.Size());
  for (i = 1; i <= p.Size(); i++)
    {
      v = p.Get(i) - rp;
      pp.Elem(i).X() = (ex_h * v);
      pp.Elem(i).Y() = (ey_h * v);
      pp.Elem(i).Z() = (ez_h * v);
    }
}

void referencetransform :: FromPlain (const Point3d & pp, Point3d & p) const
{
  Vec3d v;
  //  v = (h * pp.X()) * ex + (h * pp.Y()) * ey + (h * pp.Z()) * ez;
  //  p = rp + v;
  v.X() = pp.X() * exh.X() + pp.Y() * eyh.X() + pp.Z() * ezh.X();
  v.Y() = pp.X() * exh.Y() + pp.Y() * eyh.Y() + pp.Z() * ezh.Y();
  v.Z() = pp.X() * exh.Z() + pp.Y() * eyh.Z() + pp.Z() * ezh.Z();
  p.X() = rp.X() + v.X();
  p.Y() = rp.Y() + v.Y();
  p.Z() = rp.Z() + v.Z();
}


}
