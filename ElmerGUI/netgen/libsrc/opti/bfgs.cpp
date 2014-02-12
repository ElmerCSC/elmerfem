/***************************************************************************/
/*                                                                         */
/* Vorlesung Optimierung I, Gfrerer, WS94/95                               */
/* BFGS-Verfahren zur Lösung freier nichtlinearer Optimierungsprobleme     */
/*                                                                         */
/* Programmautor:  Joachim Schöberl                                        */
/* Matrikelnummer: 9155284                                                 */
/*                                                                         */
/***************************************************************************/

#include <mystdlib.h>
#include <myadt.hpp> 

#include <linalg.hpp>
#include "opti.hpp"


namespace netgen
{

void Cholesky (const DenseMatrix & a,
	       DenseMatrix & l, Vector & d)
{
  // Factors   A = L D L^T

  double x;

  int i, j, k;
  int n = a.Height();
  
  //  (*testout) << "a = " << a << endl;

  l = a;

  for (i = 1; i <= n; i++)
    {
      for (j = i; j <= n; j++)
	{
	  x = l.Get(i, j);

	  for (k = 1; k < i; k++)
	    x -= l.Get(i, k) * l.Get(j, k) * d.Get(k); 
	  
	  if (i == j)
	    {
	      d.Elem(i) = x;
	    }
	  else
	    {
	      l.Elem(j, i) = x / d.Get(k);
	    }
	}
    }

  for (i = 1; i <= n; i++)
    {
      l.Elem(i, i) = 1;
      for (j = i+1; j <= n; j++)
	l.Elem(i, j) = 0;
    }

  /*
  // Multiply:
  (*testout) << "multiplied factors: " << endl;
  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      {
	x = 0;
	for (k = 1; k <= n; k++)
	  x += l.Get(i, k) * l.Get(j, k) * d.Get(k);
	(*testout) << x << " ";
      }
  (*testout) << endl;
  */
}


void MultLDLt (const DenseMatrix & l, const Vector & d, const Vector & g, Vector & p)
{
  /*
  int i, j, n;
  double val;

  n = l.Height();
  p = g;
  for (i = 1; i <= n; i++)
    {
      val = 0;
      for (j = i; j <= n; j++)
	val += p.Get(j) * l.Get(j, i);
      p.Set(i, val);
    }
  for (i = 1; i <= n; i++)
    p.Elem(i) *= d.Get(i);

  for (i = n; i >= 1; i--)
    {
      val = 0;
      for (j = 1; j <= i; j++)
	val += p.Get(j) * l.Get(i, j);
      p.Set(i, val);
    }
  */



  double val;

  int n = l.Height();
  p = g;
  
  for (int i = 0; i < n; i++)
    {
      val = 0;
      for (int j = i; j < n; j++)
	val += p(j) * l(j, i);
      p(i) = val;
    }

  for (int i = 0; i < n; i++)
    p(i) *= d(i);

  for (int i = n-1; i >= 0; i--)
    {
      val = 0;
      for (int j = 0; j <= i; j++)
	val += p(j) * l(i, j);
      p(i) = val;
    }
}

void SolveLDLt (const DenseMatrix & l, const Vector & d, const Vector & g, Vector & p)
{
  double val;

  int n = l.Height();
  p = g;

  for (int i = 0; i < n; i++)
    {
      val = 0;
      for (int j = 0; j < i; j++)
	val += p(j) * l(i,j);
      p(i) -= val;
    }

  for (int i = 0; i < n; i++)
    p(i) /= d(i);
  
  for (int i = n-1; i >= 0; i--)
    {
      val = 0;
      for (int j = i+1; j < n; j++)
	val += p(j) * l(j, i);
      p(i) -= val;
    }
}

int LDLtUpdate (DenseMatrix & l, Vector & d, double a, const Vector & u)
{
  // Bemerkung: Es wird a aus R erlaubt
  // Rueckgabewert: 0 .. D bleibt positiv definit
  //                1 .. sonst

  int i, j, n;

  n = l.Height();

  Vector v(n);
  double t, told, xi;

  told = 1;
  v = u;

  for (j = 1; j <= n; j++)
    {
      t = told + a * sqr (v.Elem(j)) / d.Get(j);

      if (t <= 0) 
	{
	  (*testout) << "update err, t = " << t << endl;
	  return 1;
	}

      xi = a * v.Elem(j) / (d.Get(j) * t);

      d.Elem(j) *= t / told;

      for (i = j + 1; i <= n; i++)
	{
	  v.Elem(i) -= v.Elem(j) * l.Elem(i, j);
	  l.Elem(i, j) += xi * v.Elem(i);
	}

      told = t;
    }

  return 0;
}


double BFGS (
	     Vector & x,         // i: Startwert
	     // o: Loesung, falls IFAIL = 0
	     const MinFunction & fun,
	     const OptiParameters & par,
	     double eps
	     )


{
  int i, j, n = x.Size();
  long it;
  char a1crit, a3acrit;


  Vector d(n), g(n), p(n), temp(n), bs(n), xneu(n), y(n), s(n), x0(n);
  DenseMatrix l(n);
  DenseMatrix hesse(n);

  double /* normg, */ alphahat, hd, fold;
  double a1, a2;
  const double mu1 = 0.1, sigma = 0.1, xi1 = 1, xi2 = 10;
  const double tau = 0.1, tau1 = 0.1, tau2 = 0.6;

  Vector typx(x.Size());      // i: typische Groessenordnung der Komponenten
  double f, f0;
  double typf;               // i: typische Groessenordnung der Loesung
  double fmin = -1e5;           // i: untere Schranke fuer Funktionswert
  //  double eps = 1e-8;            // i: Abbruchschranke fuer relativen Gradienten
  double tauf = 0.1;            // i: Abbruchschranke fuer die relative Aenderung der
                                //    Funktionswerte
  int ifail;                    // o:  0 .. Erfolg
                                //    -1 .. Unterschreitung von fmin
                                //     1 .. kein Erfolg bei Liniensuche
                                //     2 .. Überschreitung von itmax

  typx = par.typx;
  typf = par.typf;


  l = 0;
  for (i = 1; i <= n; i++)
    l.Elem(i, i) = 1;

  f = fun.FuncGrad (x, g);
  f0 = f;
  x0 = x;

  it = 0;
  do
    {
      // Restart

      if (it % (5 * n) == 0)
	{

	  for (i = 1; i <= n; i++)
	    d.Elem(i) = typf/ sqr (typx.Get(i));   // 1;
	  for (i = 2; i <= n; i++)
	    for (j = 1; j < i; j++)
	      l.Elem(i, j) = 0;

	  /*
	  hesse = 0;
	  for (i = 1; i <= n; i++)
	    hesse.Elem(i, i) = typf / sqr (typx.Get(i));  

	  fun.ApproximateHesse (x, hesse);

	  Cholesky (hesse, l, d);
	  */
	}

      it++;
      if (it > par.maxit_bfgs)
	{
	  ifail = 2;
	  break;
	}


      // Solve with factorized B

      SolveLDLt (l, d, g, p);

 //      (*testout) << "l " << l << endl
// 		 << "d " << d << endl
// 		 << "g " << g << endl
// 		 << "p " << p << endl;


      p *= -1;
      y = g;

      fold = f;

      // line search

      alphahat = 1;
      lines (x, xneu, p, f, g, fun, par, alphahat, fmin,
	     mu1, sigma, xi1, xi2, tau, tau1, tau2, ifail);

      if(ifail == 1)
	(*testout) << "no success with linesearch" << endl;

       /*
      // if (it > par.maxit_bfgs/2)
	{
	  (*testout) << "x = " << x << endl;
	  (*testout) << "xneu = " << xneu << endl;
	  (*testout) << "f = " << f << endl;
	  (*testout) << "g = " << g << endl;
	}
      */

      //      (*testout) << "it = " << it << " f = " << f << endl;
      //      if (ifail != 0) break;

      s.Set2 (1, xneu, -1, x);
      y *= -1;
      y.Add (1,g); // y += g;

      x = xneu;

      // BFGS Update

      MultLDLt (l, d, s, bs);

      a1 = y * s;
      a2 = s * bs;

      if (a1 > 0 && a2 > 0)
	{
	  if (LDLtUpdate (l, d, 1 / a1, y) != 0)
	    {
              cerr << "BFGS update error1" << endl;
	      (*testout) << "BFGS update error1" << endl;
	      (*testout) << "l " << endl << l << endl
			 << "d " << d << endl;
	      ifail = 1;
	      break;
	    }

	  if (LDLtUpdate (l, d, -1 / a2, bs) != 0)
	    {
              cerr << "BFGS update error2" << endl;
	      (*testout) << "BFGS update error2" << endl;
	      (*testout) << "l " << endl << l << endl
			 << "d " << d << endl;
	      ifail = 1;
	      break;
	    }
	}

      // Calculate stop conditions

      hd = eps * max2 (typf, fabs (f));
      a1crit = 1;
      for (i = 1; i <= n; i++)
	if ( fabs (g.Elem(i)) * max2 (typx.Elem(i), fabs (x.Elem(i))) > hd)
	  a1crit = 0;


      a3acrit = (fold - f <= tauf * max2 (typf, fabs (f)));

      //    testout << "g = " << g << endl;
      //    testout << "a1crit, a3crit = " << int(a1crit) << ", " << int(a3acrit) << endl;

      /*
	// Output for tests

	normg = sqrt (g * g);

	testout << "it =" << setw (5) << it
	<< " f =" << setw (12) << setprecision (5) << f
	<< " |g| =" << setw (12) << setprecision (5) << normg;

	testout << " x = (" << setw (12) << setprecision (5) << x.Elem(1);
	for (i = 2; i <= n; i++)
	testout << "," << setw (12) << setprecision (5) << x.Elem(i);
	testout << ")" << endl;
	*/

      //(*testout) << "it = " << it << " f = " << f << " x = " << x << endl
      //	 << " g = " << g << " p = " << p << endl << endl;

      //      (*testout) << "|g| = " << g.L2Norm() << endl;

      if (g.L2Norm() < fun.GradStopping (x)) break;

    }
  while (!a1crit || !a3acrit);

  /*
  (*testout) << "it = " << it << " g = " << g << " f = " << f 
	     << " fail = " << ifail << endl;
  */
  if (f0 < f || (ifail == 1))
    {
      (*testout) << "fail, f = " << f << " f0 = " << f0 << endl;
      f = f0;
      x = x0;
    }

  //  (*testout) << "x = " << x << ", x0 = " << x0 << endl;
  return f;
}

}
