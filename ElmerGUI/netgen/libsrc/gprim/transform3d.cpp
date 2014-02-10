#include <mystdlib.h>

#include <myadt.hpp>
#include <gprim.hpp>
#include <linalg.hpp>

namespace netgen
{

Transformation3d :: Transformation3d ()
{
  int i, j;
  for (i = 0; i < 3; i++)
    {
      offset[i] = 0;
      for (j = 0; j < 3; j++)
	lin[i][j] = 0;
    }
}

Transformation3d :: Transformation3d (const Vec3d & translate)
{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      lin[i][j] = 0;
  for (i = 0; i < 3; i++)
    {
      offset[i] = translate.X(i+1);
      lin[i][i] = 1;
    }
}


Transformation3d :: 
Transformation3d (const Point3d & c, double alpha, 
		  double beta, double gamma)
{
  // total = T_c x Rot_0 x T_c^{-1}
  // Use Euler angles, see many books from tech mech, e.g. 
  // Shabana "multibody systems"

  Transformation3d tc(c);
  Transformation3d tcinv;
  tc.CalcInverse (tcinv);

  Transformation3d r1, r2, r3, ht, ht2;
  r1.SetAxisRotation (3, alpha);
  r2.SetAxisRotation (1, beta);
  r3.SetAxisRotation (3, gamma);

  ht.Combine (tc, r3);
  ht2.Combine (ht, r2);
  ht.Combine (ht2, r1);
  Combine (ht, tcinv);

 cout << "Rotation - Transformation:" << (*this) << endl;
  //  (*testout) << "Rotation - Transformation:" << (*this) << endl;
}




Transformation3d :: Transformation3d (const Point3d ** pp)
{
  int i, j;
  for (i = 1; i <= 3; i++)
    {
      offset[i-1] = (*pp[0]).X(i);
      for (j = 1; j <= 3; j++)
	lin[i-1][j-1] = (*pp[j]).X(i) - (*pp[0]).X(i);
    }
}

Transformation3d :: Transformation3d (const Point3d pp[])
{
  int i, j;
  for (i = 1; i <= 3; i++)
    {
      offset[i-1] = pp[0].X(i);
      for (j = 1; j <= 3; j++)
	lin[i-1][j-1] = pp[j].X(i) - pp[0].X(i);
    }
}


void Transformation3d :: CalcInverse (Transformation3d & inv) const
{
  static DenseMatrix a(3), inva(3);
  static Vector b(3), sol(3);
  int i, j;
  
  for (i = 1; i <= 3; i++)
    {
      b.Elem(i) = offset[i-1];
      for (j = 1; j <= 3; j++)
	a.Elem(i, j) = lin[i-1][j-1];
    }

  ::netgen::CalcInverse (a, inva);
  inva.Mult (b, sol);

  for (i = 1; i <= 3; i++)
    {
      inv.offset[i-1] = -sol.Get(i);
      for (j = 1; j <= 3; j++)
	inv.lin[i-1][j-1] = inva.Elem(i, j);
    }
}


void  Transformation3d:: 
Combine (const Transformation3d & ta, const Transformation3d & tb)
{
  int i, j, k;

  // o = o_a+ m_a o_b
  // m = m_a m_b

  for (i = 0; i <= 2; i++)
    {
      offset[i] = ta.offset[i];
      for (j = 0; j <= 2; j++)
	offset[i] += ta.lin[i][j] * tb.offset[j];
    }
  
  for (i = 0; i <= 2; i++)
    for (j = 0; j <= 2; j++)
      {
	lin[i][j] = 0;
	for (k = 0; k <= 2; k++)
	  lin[i][j] += ta.lin[i][k] * tb.lin[k][j];
      }
}
void Transformation3d :: SetAxisRotation (int dir, double alpha)
{
  double co = cos(alpha);
  double si = sin(alpha);
  dir--;
  int pos1 = (dir+1) % 3;
  int pos2 = (dir+2) % 3;

  int i, j;
  for (i = 0; i <= 2; i++)
    {
      offset[i] = 0;
      for (j = 0; j <= 2; j++)
	lin[i][j] = 0;
    }

  lin[dir][dir] = 1;
  lin[pos1][pos1] = co;
  lin[pos2][pos2] = co;
  lin[pos1][pos2] = si;
  lin[pos2][pos1] = -si;
}

ostream & operator<< (ostream & ost, Transformation3d & trans)
{
  int i, j;
  ost << "offset = ";
  for (i = 0; i <= 2; i++)
    ost << trans.offset[i] << " ";
  ost << endl << "linear = " << endl;
  for (i = 0; i <= 2; i++)
    {
      for (j = 0; j <= 2; j++)
	ost << trans.lin[i][j] << " ";
      ost << endl;
    }
  return ost;
}
}
