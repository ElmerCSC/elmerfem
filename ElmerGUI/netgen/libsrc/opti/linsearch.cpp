/***************************************************************************/
/*                                                                         */
/* Problem:        Liniensuche                                             */
/*                                                                         */
/* Programmautor:  Joachim Schöberl                                        */
/* Matrikelnummer: 9155284                                                 */
/*                                                                         */
/* Algorithmus nach:                                                       */
/*                                                                         */
/*   Optimierung I, Gfrerer, WS94/95                                       */
/*   Algorithmus 2.1: Liniensuche Problem (ii)                             */
/*                                                                         */
/***************************************************************************/



#include <mystdlib.h>

#include <myadt.hpp>  // min, max, sqr

#include <linalg.hpp>
#include "opti.hpp"


namespace netgen
{
const double eps0 = 1E-15;

// Liniensuche


double MinFunction :: Func (const Vector & /* x */) const
{
  cerr << "Func of MinFunction called" << endl;
  return 0;
}

void MinFunction :: Grad (const Vector & /* x */, Vector & /* g */) const
{
  cerr << "Grad of MinFunction called" << endl;
}
  
double MinFunction :: FuncGrad (const Vector & x, Vector & g) const
{
  cerr << "Grad of MinFunction called" << endl;
  return 0;
  /*
  int n = x.Size();

  static Vector xr;
  static Vector xl;
  xr.SetSize(n);
  xl.SetSize(n);

  double eps = 1e-6;
  double fl, fr;
  
  for (int i = 1; i <= n; i++)
    {
      xr.Set (1, x);
      xl.Set (1, x);

      xr.Elem(i) += eps;
      fr = Func (xr);

      xl.Elem(i) -= eps;
      fl = Func (xl);

      g.Elem(i) = (fr - fl) / (2 * eps);
    }

  double f = Func(x);
  //  (*testout) << "f = " << f << " grad = " << g << endl;
  return f;
  */
}


double MinFunction :: FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const
{
  Vector g(x.Size());
  double f = FuncGrad (x, g);
  deriv = (g * dir);

  //  (*testout) << "g = " << g << ", dir = " << dir << ", deriv = " << deriv << endl;
  return f;
}

void MinFunction :: ApproximateHesse (const Vector & x,
				      DenseMatrix & hesse) const
{
  int n = x.Size();
  int i, j;

  static Vector hx;
  hx.SetSize(n);

  double eps = 1e-6;
  double f, f11, f12, f21, f22;
  
  for (i = 1; i <= n; i++)
    {
      for (j = 1; j < i; j++)
	{
	  hx = x;
	  hx.Elem(i) = x.Get(i) + eps;
	  hx.Elem(j) = x.Get(j) + eps;
	  f11 = Func(hx);
	  hx.Elem(i) = x.Get(i) + eps;
	  hx.Elem(j) = x.Get(j) - eps;
	  f12 = Func(hx);
	  hx.Elem(i) = x.Get(i) - eps;
	  hx.Elem(j) = x.Get(j) + eps;
	  f21 = Func(hx);
	  hx.Elem(i) = x.Get(i) - eps;
	  hx.Elem(j) = x.Get(j) - eps;
	  f22 = Func(hx);

	  hesse.Elem(i, j) = hesse.Elem(j, i) =
	    (f11 + f22 - f12 - f21) / (2 * eps * eps);
	}

      hx = x;
      f = Func(x);
      hx.Elem(i) = x.Get(i) + eps;
      f11 = Func(hx);
      hx.Elem(i) = x.Get(i) - eps;
      f22 = Func(hx);

      hesse.Elem(i, i) = (f11 + f22 - 2 * f) / (eps * eps);
    }
  //  (*testout) << "hesse = " << hesse << endl;
}







/// Line search, modified Mangasarien conditions
void lines (Vector & x,         // i: initial point of line-search
	    Vector & xneu,      // o: solution, if successful
	    Vector & p,         // i: search direction
	    double & f,         // i: function-value at x
	    // o: function-value at xneu, iff ifail = 0
	    Vector & g,         // i: gradient at x
	    // o: gradient at xneu, iff ifail = 0
	    const MinFunction & fun,  // function to minimize
	    const OptiParameters & par,
	    double & alphahat,  // i: initial value for alpha_hat
	    // o: solution alpha iff ifail = 0
	    double fmin,        // i: lower bound for f
	    double mu1,         // i: Parameter mu_1 of Alg.2.1
	    double sigma,       // i: Parameter sigma of Alg.2.1
	    double xi1,         // i: Parameter xi_1 of Alg.2.1
	    double xi2,         // i: Parameter xi_1 of Alg.2.1
	    double tau,         // i: Parameter tau of Alg.2.1
	    double tau1,        // i: Parameter tau_1 of Alg.2.1
	    double tau2,        // i: Parameter tau_2 of Alg.2.1
	    int & ifail)        // o: 0 on success
  //    -1 bei termination because lower limit fmin
  //     1 bei illegal termination due to different reasons

{
  double phi0, phi0prime, phi1, phi1prime, phihatprime;
  double alpha1, alpha2, alphaincr, c;
  char flag = 1;
  long it;

  alpha1 = 0;
  alpha2 = 1e50;
  phi0 = phi1 = f;

  phi0prime = g * p;

  if (phi0prime > 0)
    {
      ifail = 1;
      return;
    }

  ifail = 1;  // Markus

  phi1prime = phi0prime;

  //  (*testout) << "phi0prime = " << phi0prime << endl;

  //  it = 100000l;
  it = 0;

  while (it++ <= par.maxit_linsearch)
    {

      xneu.Set2 (1, x, alphahat, p);


      //    f = fun.FuncGrad (xneu, g);
      //      f = fun.Func (xneu);
      f = fun.FuncDeriv (xneu, p, phihatprime);

      // (*testout) << "lines, f = " << f << " phip = " << phihatprime << endl;

      if (f < fmin)
	{
	  ifail = -1;
	  break;
	}


      if (alpha2 - alpha1 < eps0 * alpha2)
	{
	  ifail = 0;
	  break;
	}

      // (*testout) << "i = " << it << " al = " << alphahat << " f = " << f << " fprime " << phihatprime << endl;;

      if (f - phi0 > mu1 * alphahat * phi1prime + eps0 * fabs (phi0))

	{

	  flag = 0;
	  alpha2 = alphahat;

	  c = 
	    (f - phi1 - phi1prime * (alphahat-alpha1)) / 
	    sqr (alphahat-alpha1);

	  alphahat = alpha1 - 0.5 * phi1prime / c;

	  if (alphahat > alpha2)
	    alphahat = alpha1 + 1/(4*c) *
	      ( (sigma+mu1) * phi0prime - 2*phi1prime
		+ sqrt (sqr(phi1prime - mu1 * phi0prime) -
			4 * (phi1 - phi0 - mu1 * alpha1 * phi0prime) * c));

	  alphahat = max2 (alphahat, alpha1 + tau * (alpha2 - alpha1));
	  alphahat = min2 (alphahat, alpha2 - tau * (alpha2 - alpha1));
	  
	  //	  (*testout) << " if-branch" << endl;

	}

      else

	{
	  /*
	  f = fun.FuncGrad (xneu, g);
	  phihatprime = g * p;
	  */
	  f = fun.FuncDeriv (xneu, p, phihatprime);

	  if (phihatprime < sigma * phi0prime * (1 + eps0))

	    {
	      if (phi1prime < phihatprime)   
		// Approximationsfunktion ist konvex

		alphaincr = (alphahat - alpha1) * phihatprime /
		  (phi1prime - phihatprime);

	      else
		alphaincr = 1e99; // MAXDOUBLE;

	      if (flag)
		{
		  alphaincr = max2 (alphaincr, xi1 * (alphahat-alpha1));
		  alphaincr = min2 (alphaincr, xi2 * (alphahat-alpha1));
		}
	      else
		{
		  alphaincr = max2 (alphaincr, tau1 * (alpha2 - alphahat));
		  alphaincr = min2 (alphaincr, tau2 * (alpha2 - alphahat));
		}

	      alpha1 = alphahat;
	      alphahat += alphaincr;
	      phi1 = f;
	      phi1prime = phihatprime;
	    }

	  else

	    {
	      ifail = 0;     // Erfolg !!
	      break;
	    }
	  
	  //	  (*testout) << " else, " << endl;

	}

    }

  //  (*testout) << "linsearch: it = " << it << " ifail = " << ifail << endl;

  fun.FuncGrad (xneu, g);


  if (it < 0)
    ifail = 1;

  //  (*testout) << "fail = " << ifail << endl;
}



















void SteepestDescent (Vector & x, const MinFunction & fun,
		      const OptiParameters & par)
{
  int it, n = x.Size();
  Vector xnew(n), p(n), g(n), g2(n);
  double val, alphahat;
  int fail;

  val = fun.FuncGrad(x, g);

  alphahat = 1;
  //  testout << "f = ";
  for (it = 0; it < 10; it++)
    {
      //    testout << val << " ";

      // p = -g;
      p.Set (-1, g);

      lines (x, xnew, p, val, g, fun, par, alphahat, -1e5,
	     0.1, 0.1, 1, 10, 0.1, 0.1, 0.6, fail);

      x = xnew;
    }
  //  testout << endl;
}
}
