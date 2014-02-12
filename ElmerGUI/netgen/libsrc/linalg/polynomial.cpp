#include <mystdlib.h>
#include <linalg.hpp>

namespace netgen
{

QuadraticPolynomial1V :: 
QuadraticPolynomial1V (double ac, double acx, double acxx)
{
  c = ac;
  cx = acx;
  cxx = acxx;
}

double QuadraticPolynomial1V :: 
Value (double x)
{
  return c + cx * x + cxx * x * x;
}

double QuadraticPolynomial1V ::  MaxUnitInterval ()
{
  // inner max
  if (cxx < 0 && cx > 0 && cx < -2 * cxx)
    {
      return c - 0.25 * cx * cx / cxx;
    }

  
  if (cx + cxx > 0)   // right edge
    return c + cx + cxx;

  // left end
  return c;
}




LinearPolynomial2V :: 
LinearPolynomial2V (double ac, double acx, double acy)
{
  c = ac;
  cx = acx;
  cy = acy;
};


QuadraticPolynomial2V ::   
QuadraticPolynomial2V ()
{
  ;
}


QuadraticPolynomial2V :: 
QuadraticPolynomial2V (double ac, double acx, double acy,
		       double acxx, double acxy, double acyy)
{
  c = ac;
  cx = acx;
  cy = acy;
  cxx = acxx;
  cxy = acxy;
  cyy = acyy;
}

void QuadraticPolynomial2V :: 
Square (const LinearPolynomial2V & lp)
{
  c = lp.c * lp.c;
  cx = 2 * lp.c * lp.cx;
  cy = 2 * lp.c * lp.cy;

  cxx = lp.cx * lp.cx;
  cxy = 2 * lp.cx * lp.cy;
  cyy = lp.cy * lp.cy;
}

void QuadraticPolynomial2V :: 
Add (double lam, const QuadraticPolynomial2V & qp2)
{
  c += lam * qp2.c;
  cx += lam * qp2.cx;
  cy += lam * qp2.cy;
  cxx += lam * qp2.cxx;
  cxy += lam * qp2.cxy;
  cyy += lam * qp2.cyy;
}

double QuadraticPolynomial2V :: 
Value (double x, double y)
{
  return c + cx * x + cy * y + cxx * x * x + cxy * x * y + cyy * y * y;
}

/*
double QuadraticPolynomial2V :: 
MinUnitSquare ()
{
  double x, y;
  double minv = 1e8;
  double val;
  for (x = 0; x <= 1; x += 0.1)
    for (y = 0; y <= 1; y += 0.1)
      {
	val = Value (x, y);
	if (val < minv)
	  minv = val;
      }
  return minv;
};
*/

double QuadraticPolynomial2V :: 
MaxUnitSquare ()
{
  // find critical point

  double maxv = c;
  double hv;

  double det, x0, y0;
  det = 4 * cxx * cyy - cxy * cxy;

  if (det > 0)
    {
      // definite surface
      
      x0 = (-2 * cyy * cx + cxy * cy) / det;
      y0 = (cxy * cx -2 * cxx * cy) / det;

      if (x0 >= 0 && x0 <= 1 && y0 >= 0 && y0 <= 1)
	{
	  hv = Value (x0, y0);
	  if (hv > maxv) maxv = hv;
	}
    }
  
  QuadraticPolynomial1V e1(c, cx, cxx);
  QuadraticPolynomial1V e2(c, cy, cyy);
  QuadraticPolynomial1V e3(c+cy+cyy, cx+cxy, cxx);
  QuadraticPolynomial1V e4(c+cx+cxx, cy+cxy, cyy);
  
  hv = e1.MaxUnitInterval();
  if (hv > maxv) maxv = hv;
  hv = e2.MaxUnitInterval();
  if (hv > maxv) maxv = hv;
  hv = e3.MaxUnitInterval();
  if (hv > maxv) maxv = hv;
  hv = e4.MaxUnitInterval();
  if (hv > maxv) maxv = hv;

  return maxv;

  //  (*testout) << "maxv = " << maxv << " =~= ";

  /*
  double x, y;
  maxv = -1e8;
  double val;
  for (x = 0; x <= 1.01; x += 0.1)
    for (y = 0; y <= 1.01; y += 0.1)
      {
	val = Value (x, y);
	if (val > maxv)
	  maxv = val;
      }

  //  (*testout) << maxv << endl;
  return maxv;
  */
};




double QuadraticPolynomial2V :: 
MaxUnitTriangle ()
{
  // find critical point
  
  double maxv = c;
  double hv;

  double det, x0, y0;
  det = 4 * cxx * cyy - cxy * cxy;

  if (cxx < 0 && det > 0)
    { 
      // definite surface
      
      x0 = (-2 * cyy * cx + cxy * cy) / det;
      y0 = (cxy * cx -2 * cxx * cy) / det;

      if (x0 >= 0 && y0 >= 0 && x0+y0 <= 1)
	{
	  return Value (x0, y0);
	}
    }
  
  
  QuadraticPolynomial1V e1(c, cx, cxx);
  QuadraticPolynomial1V e2(c, cy, cyy);
  QuadraticPolynomial1V e3(c+cy+cyy, cx-cy+cxy-2*cyy, cxx-cxy+cyy);
  
  hv = e1.MaxUnitInterval();
  if (hv > maxv) maxv = hv;
  hv = e2.MaxUnitInterval();
  if (hv > maxv) maxv = hv;
  hv = e3.MaxUnitInterval();
  if (hv > maxv) maxv = hv;

  return maxv;
}
}
