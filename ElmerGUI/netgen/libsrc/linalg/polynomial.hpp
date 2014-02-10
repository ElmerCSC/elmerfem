#ifndef FILE_POLYNOMIAL
#define FILE_POLYNOMIAL

/* *************************************************************************/
/* File:   polynomial.hh                                                   */
/* Author: Joachim Schoeberl                                               */
/* Date:   25. Nov. 99                                                     */
/* *************************************************************************/


class QuadraticPolynomial1V 
{
  double c, cx, cxx;
public:
  QuadraticPolynomial1V (double ac, double acx, double acxx);
  double Value (double x);
  double MaxUnitInterval ();
};

class LinearPolynomial2V
{
  double c, cx, cy;
public:
  LinearPolynomial2V (double ac, double acx, double acy);
  friend class QuadraticPolynomial2V;
};


class QuadraticPolynomial2V
{
  double c, cx, cy, cxx, cxy, cyy;
public:
  QuadraticPolynomial2V ();
  QuadraticPolynomial2V (double ac, double acx, double acy,
			 double acxx, double acxy, double acyy);
  void Square (const LinearPolynomial2V & lp);
  void Add (double lam, const QuadraticPolynomial2V & qp);

  double Value (double x, double y);
  //  double MinUnitSquare ();
  double MaxUnitSquare ();
  double MaxUnitTriangle ();
};

#endif
