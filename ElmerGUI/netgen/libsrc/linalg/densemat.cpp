#include <mystdlib.h>

#include <linalg.hpp>


namespace netgen
{
  DenseMatrix :: DenseMatrix () 
  {
    data = NULL;
    height = 0;
    width = 0;
  }

  DenseMatrix :: DenseMatrix (int h, int w)
  {
    if (!w) w = h;
    width = w;
    height = h;
    if (h*w)
      data = new double[h*w];
    else 
      data = 0;

    for (int i = 0 ; i < (h * w); i++)
      data[i] = 0;
  }

  /*
  DenseMatrix :: DenseMatrix (int h, int w, const double * d) 
    : BaseMatrix (h, w)
  {
  int size = h * w;  
  int i;
  
  if (size)
    {
      data = new double[size]; 
      for (i = 0; i < size; i++)
	data[i] = d[i];
    }
  else
    data = NULL;
  }    
  */

  DenseMatrix :: DenseMatrix (const DenseMatrix & m2)
  {
    data = NULL; height = width = 0;
    SetSize (m2.Height(), m2.Width());
    memcpy (data, m2.data, sizeof(double) * Height() * Width());
  }

  DenseMatrix :: ~DenseMatrix ()
  {
    delete [] data;
  }
  
  
  void DenseMatrix :: SetSize (int h, int w)
  {
    if (!w) w = h;
    if (height == h && width == w)
      return;
          
    height = h;
    width = w;
    
    delete[] data;
    
    if (h*w)  
      data = new double[h*w];
    else
      data = NULL;
  }


  /*
DenseMatrix & DenseMatrix :: operator= (const BaseMatrix & m2)
  {
  int i, j;

  SetSize (m2.Height(), m2.Width());

  if (data)
    for (i = 1; i <= Height(); i++)
      for (j = 1; j <= Width(); j++)
        Set (i, j, m2(i, j));
  else
    (*myerr) << "DenseMatrix::Operator=: Matrix not allocated" << endl;

  return *this;
  }
  */


  DenseMatrix & DenseMatrix :: operator= (const DenseMatrix & m2)
  {
    SetSize (m2.Height(), m2.Width());
    
    if (data) memcpy (data, m2.data, sizeof(double) * m2.Height() * m2.Width());
    return *this;
  }


  DenseMatrix & DenseMatrix :: operator+= (const DenseMatrix & m2)
  {
    int i;
    double * p, * q;
    
    if (Height() != m2.Height() || Width() != m2.Width())
    {
      (*myerr) << "DenseMatrix::Operator+=: Sizes don't fit" << endl;
      return *this;
    }
    
    if (data)
      {
	p = data;
	q = m2.data;
	for (i = Width() * Height(); i > 0; i--)
      {
      *p += *q;
      p++;
      q++;
      }
    }
  else
    (*myerr) << "DenseMatrix::Operator+=: Matrix not allocated" << endl;

  return *this;
  }


DenseMatrix & DenseMatrix :: operator-= (const DenseMatrix & m2)
  {
  int i;
  double * p, * q;

  if (Height() != m2.Height() || Width() != m2.Width())
    {
    (*myerr) << "DenseMatrix::Operator-=: Sizes don't fit" << endl;
    return *this;
    }

  if (data)
    {
    p = data;
    q = m2.data;
    for (i = Width() * Height(); i > 0; i--)
      {
      *p -= *q;
      p++;
      q++;
      }
    }
  else
    (*myerr) << "DenseMatrix::Operator-=: Matrix not allocated" << endl;

  return *this;
  }




  /*
double & DenseMatrix :: operator() (int i, int j)
{
  if (i >= 1 && j >= 1 && i <= height && j <= width)
    return Elem(i,j);
  else (*myerr) << "DenseMatrix: index (" << i << "," << j << ") out of range (1.."
		<< height << ",1.." << width << ")\n";
  static double dummy = 0;
  return dummy;
}

  double DenseMatrix :: operator() (int i, int j) const
  {
    if (i >= 1 && j >= 1 && i <= height && j <= width)
      return Get(i,j);
    else (*myerr) << "DenseMatrix: index (" << i << "," << j << ") out of range (1.."
            << height << ",1.." << width << ")\n";

    static double dummy = 0;
    return dummy;
  }
  */

DenseMatrix & DenseMatrix :: operator= (double v)
  {
  int i;
  double * p = data;

  if (data)
    for (i = width*height; i > 0; i--, p++)
      *p = v;

  return *this;
  }



DenseMatrix & DenseMatrix :: operator*= (double v)
  {
  int i;
  double * p = data;

  if (data)
    for (i = width*height; i > 0; i--, p++)
      *p *= v;

  return *this;
  }


double DenseMatrix :: Det () const
  {
  if (width != height)
    {
    (*myerr) << "DenseMatrix :: Det: width != height" << endl;
    return 0;
    }

  switch (width)
    {
    case 1: return Get(1, 1);
    case 2: return Get(1) * Get(4) - Get(2) * Get(3);

    case 3: return Get(1) * Get(5) * Get(9)
                 + Get(2) * Get(6) * Get(7)
                 + Get(3) * Get(4) * Get(8)
                 - Get(1) * Get(6) * Get(8)
                 - Get(2) * Get(4) * Get(9)
                 - Get(3) * Get(5) * Get(7);
    default:
      {
      (*myerr) << "Matrix :: Det:  general size not implemented (size=" << width << ")" << endl;
      return 0;
      }
    }
  }


void CalcInverse (const DenseMatrix & m1, DenseMatrix & m2)
  {
    //  int i, j, k, n;
  double det;
  //  DenseMatrix m1 = hm1;

  if (m1.width != m1.height)
    {
    (*myerr) << "CalcInverse: matrix not symmetric" << endl;
    return;
    }
  if (m1.width != m2.width || m1.height != m2.height)
    {
    (*myerr) << "CalcInverse: dim(m2) != dim(m1)" << endl;
    return;
    }


  if (m1.Width() <= 3)
    {
    det = m1.Det();
    if (det == 0)
      {
      (*myerr) << "CalcInverse: Matrix singular" << endl;
      return;
      }

    det = 1e0 / det;
    switch (m1.width)
      {
      case 1:
        {
        m2.Set(1, 1, det);
        return;
        }
      case 2:
        {
        m2.Set(1, 1, det * m1.Get(4));
        m2.Set(2, 2, det * m1.Get(1));  
        m2.Set(1, 2, - det * m1.Get(2));
        m2.Set(2, 1, - det * m1.Get(3));
        return;
        }
      case 3:
        {
        m2.Set(1, 1,  det * (m1.Get(5) * m1.Get(9) - m1.Get(6) * m1.Get(8)));
        m2.Set(2, 1, -det * (m1.Get(4) * m1.Get(9) - m1.Get(6) * m1.Get(7)));
        m2.Set(3, 1,  det * (m1.Get(4) * m1.Get(8) - m1.Get(5) * m1.Get(7)));

        m2.Set(1, 2, -det * (m1.Get(2) * m1.Get(9) - m1.Get(3) * m1.Get(8)));
        m2.Set(2, 2,  det * (m1.Get(1) * m1.Get(9) - m1.Get(3) * m1.Get(7)));
        m2.Set(3, 2, -det * (m1.Get(1) * m1.Get(8) - m1.Get(2) * m1.Get(7)));

        m2.Set(1, 3,  det * (m1.Get(2) * m1.Get(6) - m1.Get(3) * m1.Get(5)));
        m2.Set(2, 3, -det * (m1.Get(1) * m1.Get(6) - m1.Get(3) * m1.Get(4)));
        m2.Set(3, 3,  det * (m1.Get(1) * m1.Get(5) - m1.Get(2) * m1.Get(4)));
        return;
        }
      }
    }
    
  else
    {
      int i, j, k, n;
      n = m1.Height();
      

#ifdef CHOL
      int dots = (n > 200);

      // Cholesky
      
      double x;
      Vector p(n);

      m2 = m1;
      /*
      m2.SetSymmetric();
      if (!m2.Symmetric())
	cerr << "m should be symmetric for Cholesky" << endl;
      */

      for (i = 1; i <= n; i++)
	for (j = 1; j < i; j++)
	  m2.Elem(j, i) = m2.Get(i, j);
      
      for (i = 1; i <= n; i++)
	{
	  if (dots && i % 10 == 0)
	    (*mycout) << "." << flush;

	  for (j = i; j <= n; j++)
	    {
	      x = m2.Get(i, j);

	      const double * pik = &m2.Get(i, 1);
	      const double * pjk = &m2.Get(j, 1);

	      for (k = i-2; k >= 0; --k, ++pik, ++pjk)
		x -= (*pik) * (*pjk);
		  
	      // for (k = i-1; k >= 1; --k)
	      //   x -= m2.Get(j, k) * m2.Get(i, k);

	      if (i == j)
		{
		  if (x <= 0)
		    {
		      cerr << "Matrix indefinite 1" << endl;
		      return;
		    }
		  
		  p.Elem(i) = 1 / sqrt(x);
		}
	      else
		{
		  m2.Elem(j, i) = x * p.Get(i);
		}
	    }
	}

      for (i = 1; i <= n; i++)
	m2.Elem(i, i) = 1 / p.Get(i);

      // check: A = L L^t

//       for (i = 1; i <= n; i++)
// 	for (j = 1; j <= n; j++)
// 	  {
// 	    x = 0;
// 	    for (k = 1; k <= i && k <= j; k++)
// 	      x += m2.Get(i, k) * m2.Get(j, k);
// 	    (*testout) << "err " << i << "," << j << " = " << (m1.Get(i, j) - x) << endl;
// 	  }


      
      // calc L^{-1}, store upper triangle
      
      //      DenseMatrix hm(n);
      //      hm = m2;

      for (i = 1; i <= n; i++)
	{
	  if (dots && i % 10 == 0)
	    (*mycout) << "+" << flush;

	  for (j = i; j <= n; j++)
	    {
	      x = 0;
	      if (j == i) x = 1;

	      const double * pjk = &m2.Get(j, i);
	      const double * pik = &m2.Get(i, i);
	      for (k = i; k < j; k++, ++pjk, ++pik)
		x -= *pik * *pjk;

	      //  for (k = i; k < j; k++)
	      //  x -= m2.Get(j, k) * m2.Get(i, k);

	      m2.Elem(i, j) = x / m2.Get(j, j);
	    }
	}
      
//      (*testout) << "check L^-1" << endl;
//      for (i = 1; i <= n; i++)
// 	for (j = 1; j <= n; j++)
// 	  {
// 	    x = 0;
// 	    for (k = j; k <= i; k++)
// 	      x += hm.Get(i, k) * m2.Get(j, k);
// 	    (*testout) << "i, j = " << i << "," << j << " x = " << x << endl;
// 	  }


      // calc A^-1 = L^-T * L^-1

      for (i = 1; i <= n; i++)
	{
	  if (dots && i % 10 == 0)
	    (*mycout) << "-" << flush;

	  for (j = 1; j <= i; j++)
	    {
	      x = 0;
	      k = i;
	      if (j > i) k = j;

	      const double * pik = &m2.Get(i, k);
	      const double * pjk = &m2.Get(j, k);

	      for ( ; k <= n; ++k, ++pik, ++pjk)
		x += *pik * *pjk;
	      // for (  ; k <= n; k++)
	      //   x += m2.Get(i, k) * m2.Get(j, k);
	      
	      m2.Elem(i, j) = x;
	    }
	}
	  
      for (i = 1; i <= n; i++)
	for (j = 1; j < i; j++)
	  m2.Elem(j, i) = m2.Get(i, j);
      
      if (dots) (*mycout) << endl;
#endif



      // Gauss - Jordan - algorithm
      
      int r, hi;
      double max, hr;
      

      ARRAY<int> p(n);   // pivot-permutation
      Vector hv(n);
    
      
      m2 = m1;

      /*      
      if (m2.Symmetric())
	for (i = 1; i <= n; i++)
	  for (j = 1; j < i; j++)
	    m2.Elem(j, i) = m2.Get(i, j);
      */
      
    // Algorithm of Stoer, Einf. i. d. Num. Math, S 145
      
      for (j = 1; j <= n; j++)
	p.Set(j, j);
      
      for (j = 1; j <= n; j++)
	{
	  // pivot search
	  
	  max = fabs(m2.Get(j, j));
	  r = j;
	  
	  for (i = j+1; i <= n ;i++)
	    if (fabs (m2.Get(i, j)) > max)
	      {
		r = i;
		max = fabs (m2.Get(i, j));
	      }
	  
	  if (max < 1e-20)
	    {
	      cerr << "Inverse matrix: matrix singular" << endl;
	      return;
	    }
	  
	  r = j;
	  
	  // exchange rows
	  if (r > j)
	    {
	      for (k = 1; k <= n; k++)
		{
		  hr = m2.Get(j, k);
		  m2.Elem(j, k) = m2.Get(r, k);
		  m2.Elem(r, k) = hr;
		}
	      hi = p.Get(j);
	      p.Elem(j) = p.Get(r);
	      p.Elem(r) = hi;
	    }
	  
	  
	  // transformation
	  
	  hr = 1 / m2.Get(j, j);
	  for (i = 1; i <= n; i++)
	    m2.Elem(i, j) *= hr;
	  m2.Elem(j, j) = hr;
	  
	  for (k = 1; k <= n; k++)
	    if (k != j)
	      {
		for (i = 1; i <= n; i++)
		  if (i != j)
		    m2.Elem(i, k) -= m2.Elem(i, j) * m2.Elem(j, k);
		m2.Elem(j, k) *= -hr;
	      }
	}
      
      // col exchange
      
      for (i = 1; i <= n; i++)
	{
	  for (k = 1; k <= n; k++)
	    hv.Elem(p.Get(k)) = m2.Get(i, k);
	  for (k = 1; k <= n; k++)
	    m2.Elem(i, k) = hv.Get(k);
	}



    /*
    if (m1.Symmetric())
      for (i = 1; i <= n; i++)
	for (j = 1; j < i; j++)
	  m1.Elem(j, i) = m1.Get(i, j);

    m2 = 0;
    
    for (i = 1; i <= n; i++)
      m2.Elem(i, i) = 1;
      
    for (i = 1; i <= n; i++)
      {
	//	(*mycout) << '.' << flush;
      q = m1.Get(i, i);
      for (k = 1; k <= n; k++)
        {
        m1.Elem(i, k) /= q;
        m2.Elem(i, k) /= q;
        }
        
      for (j = i+1; j <= n; j++)
        {
        q = m1.Elem(j, i);

	double * m1pi = &m1.Elem(i, i);
	double * m1pj = &m1.Elem(j, i);

	for (k = n; k >= i; --k, ++m1pi, ++m1pj)
	    *m1pj -= q * (*m1pi);

	double * m2pi = &m2.Elem(i, 1);
	double * m2pj = &m2.Elem(j, 1);

	for (k = i; k > 0; --k, ++m2pi, ++m2pj)
	    *m2pj -= q * (*m2pi);

	    //        for (k = 1; k <= n; k++)  
	    //          {
	    //          m1.Elem(j, k) -= q * m1.Elem(i, k);
	    //          m2.Elem(j, k) -= q * m2.Elem(i, k);
	    //          }
	  
        }
      }  
            
    for (i = n; i >= 1; i--)
      {
	//	(*mycout) << "+" << flush;
	for (j = 1; j < i; j++)
	  {
	    q = m1.Elem(j, i);

	    double * m2pi = &m2.Elem(i, 1);
	    double * m2pj = &m2.Elem(j, 1);

	    for (k = n; k > 0; --k, ++m2pi, ++m2pj)
	      *m2pj -= q * (*m2pi);	    

	    
	    //	    for (k = 1; k <= n; k++)
	    //	      {
	    //		m1.Elem(j, k) -= q * m1.Elem(i, k);
	    //		m2.Elem(j, k) -= q * m2.Elem(i, k);
	    //	      }    
	  }         
      }

    if (m2.Symmetric())
      {
	for (i = 1; i <= n; i++)
	  for (j = 1; j < i; j++)
	    m2.Elem(i, j) = m2.Elem(j, i);
      }
*/
    }
  }


void CalcAAt (const DenseMatrix & a, DenseMatrix & m2)
  {
  int n1 = a.Height();
  int n2 = a.Width();
  int i, j, k;
  double sum;
  const double *p, *q, *p0;

  if (m2.Height() != n1 || m2.Width() != n1)
    {
    (*myerr) << "CalcAAt: sizes don't fit" << endl;
    return;
    }

  for (i = 1; i <= n1; i++)
    {
    sum = 0;
    p = &a.ConstElem(i, 1);
    for (k = 1; k <= n2; k++)
      {
      sum += *p * *p;
      p++;
      }
    m2.Set(i, i, sum);

    p0 = &a.ConstElem(i, 1);
    q = a.data;
    for (j = 1; j < i; j++)
      {
      sum = 0;
      p = p0;

      for (k = 1; k <= n2; k++)
        {
        sum += *p * *q;
        p++;
        q++;
        }
      m2.Set(i, j, sum);
      m2.Set(j, i, sum);
      }
    }
  }



#ifdef ABC
BaseMatrix * DenseMatrix :: InverseMatrix (const BitArray * /* inner */) const
{
  if (Height() != Width())
    {
      (*myerr) << "BaseMatrix::InverseMatrix(): Matrix not symmetric" << endl;
      return new DenseMatrix(1);
    }
  else
    {
      if (Symmetric())
	{	
	  (*mycout) << "Invmat not available" << endl;
	  BaseMatrix * invmat = NULL;
	  return invmat;
	}

      DenseMatrix * invmat = new DenseMatrix (Height());

      CalcInverse (*this, *invmat);
      return invmat;
    }
}
#endif



void CalcAtA (const DenseMatrix & a, DenseMatrix & m2)
  {
  int n1 = a.Height();
  int n2 = a.Width();
  int i, j, k;
  double sum;

  if (m2.Height() != n2 || m2.Width() != n2)
    {
    (*myerr) << "CalcAtA: sizes don't fit" << endl;
    return;
    }

  for (i = 1; i <= n2; i++)
    for (j = 1; j <= n2; j++)
      {
      sum = 0;
      for (k = 1; k <= n1; k++)
        sum += a.Get(k, i) * a.Get(k, j);
      m2.Elem(i, j) = sum;
      }
  }






void CalcABt (const DenseMatrix & a, const DenseMatrix & b, DenseMatrix & m2)
  {
  int n1 = a.Height();
  int n2 = a.Width();
  int n3 = b.Height();
  int i, j, k;
  double sum;

  if (m2.Height() != n1 || m2.Width() != n3 || b.Width() != n2)
    {
    (*myerr) << "CalcABt: sizes don't fit" << endl;
    return;
    }

  double * pm2 = &m2.Elem(1, 1);
  const double * pa1 = &a.Get(1, 1);

  for (i = 1; i <= n1; i++)
    {
      const double * pb = &b.Get(1, 1);
      for (j = 1; j <= n3; j++)
	{
	  sum = 0;
	  const double * pa = pa1;
	  
	  for (k = 1; k <= n2; k++)
	    {
	      sum += *pa * *pb;
	      pa++; pb++;
	    }
	  
	  *pm2 = sum;
	  pm2++;
	}
      pa1 += n2;
    }
  }


void CalcAtB (const DenseMatrix & a, const DenseMatrix & b, DenseMatrix & m2)
  {
  int n1 = a.Height();
  int n2 = a.Width();
  int n3 = b.Width();
  int i, j, k;

  if (m2.Height() != n2 || m2.Width() != n3 || b.Height() != n1)
    {
    (*myerr) << "CalcAtB: sizes don't fit" << endl;
    return;
    }

  for (i = 1; i <= n2 * n3; i++)
    m2.data[i-1] = 0;

  for (i = 1; i <= n1; i++)
    for (j = 1; j <= n2; j++)
      {
	const double va = a.Get(i, j);
	double * pm2 = &m2.Elem(j, 1);
	const double * pb = &b.Get(i, 1);

	for (k = 1; k <= n3; ++k, ++pm2, ++pb)
	  *pm2 += va * *pb;
	//	for (k = 1; k <= n3; k++)
	//	  m2.Elem(j, k) += va * b.Get(i, k);
      }
  /*
  for (i = 1; i <= n2; i++)
    for (j = 1; j <= n3; j++)
      {
	sum = 0;
	for (k = 1; k <= n1; k++)
	  sum += a.Get(k, i) * b.Get(k, j);
	m2.Elem(i, j) = sum;
      }
      */
  }







DenseMatrix operator* (const DenseMatrix & m1, const DenseMatrix & m2)
  {
  DenseMatrix temp (m1.Height(), m2.Width());

  if (m1.Width() != m2.Height())
    {
    (*myerr) << "DenseMatrix :: operator*: Matrix Size does not fit" << endl;
    }
  else if (temp.Height() != m1.Height())
    {
    (*myerr) << "DenseMatrix :: operator*: temp not allocated" << endl;
    }
  else
    {
    Mult (m1, m2, temp);
    }
  return temp;
  }


void Mult (const DenseMatrix & m1, const DenseMatrix & m2, DenseMatrix & m3)
  {
  double sum;
  double *p1, *p1s, *p1sn, *p1snn, *p2, *p2s, *p2sn, *p3;

  if (m1.Width() != m2.Height() || m1.Height() != m3.Height() ||
       m2.Width() != m3.Width() )
    {
    (*myerr) << "DenseMatrix :: Mult: Matrix Size does not fit" << endl;
    (*myerr) << "m1: " << m1.Height() << " x " << m1.Width() << endl;
    (*myerr) << "m2: " << m2.Height() << " x " << m2.Width() << endl;
    (*myerr) << "m3: " << m3.Height() << " x " << m3.Width() << endl;
    return;
    }
  /*
  else if (m1.Symmetric() || m2.Symmetric() || m3.Symmetric())
    {
    (*myerr) << "DenseMatrix :: Mult: not implemented for symmetric matrices" << endl;
    return;
    }
  */
  else
    {
      //      int i, j, k;
      int n1 = m1.Height();
      int n2 = m2.Width();
      int n3 = m1.Width();

      /*
      for (i = n1 * n2-1; i >= 0; --i)
	m3.data[i] = 0;

      const double * pm1 = &m1.Get(1, 1);
      for (i = 1; i <= n1; i++)
	{
	  const double * pm2 = &m2.Get(1, 1);
	  double * pm3i = &m3.Elem(i, 1);

	  for (j = 1; j <= n3; j++)
	    {
	      const double vm1 = *pm1;
	      ++pm1;
	      //	      const double vm1 = m1.Get(i, j);
	      double * pm3 = pm3i;
	      //	      const double * pm2 = &m2.Get(j, 1);

	      for (k = 0; k < n2; k++)
		{
		  *pm3 += vm1 * *pm2;
		  ++pm2;
		  ++pm3;
		}

	    //	    for (k = 1; k <= n2; k++)
	    //	      m3.Elem(i, k) += m1.Get(i, j) * m2.Get(j, k);
	    }
	}
	*/

      /*
      for (i = 1; i <= n1; i++)
	for (j = 1; j <= n2; j++)
	  {
	    sum = 0;
	    for (k = 1; k <= n3; k++)
	      sum += m1.Get(i, k) * m2.Get(k, j);
	    m3.Set(i, j, sum);
	  }
	  */


      /*
      for (i = 1; i <= n1; i++)
	{
	  const double pm1i = &m1.Get(i, 1);
	  const double pm2j = &m2.Get(1, 1);

	  for (j = 1; j <= n2; j++)
	    {
	      double sum = 0;
	      const double * pm1 = pm1i;
	      const double * pm2 = pm2j;
	      pm2j++;

	      for (k = 1; k <= n3; k++)
		{
		  sum += *pm1 * *pm2;
		  ++pm1;
		  pm2 += n2;
		}
	      
	      m3.Set (i, j, sum);
	    }
	}
	*/


      p3 = m3.data;
      p1s = m1.data;
      p2sn = m2.data + n2;
      p1snn = p1s + n1 * n3;

      while (p1s != p1snn)
	{
	  p1sn = p1s + n3;
	  p2s = m2.data;
	  
	  while (p2s != p2sn)
	    {
	      sum = 0;
	      p1 = p1s;
	      p2 = p2s;
	      p2s++;

	      while (p1 != p1sn)
		{
		  sum += *p1 * *p2;
		  p1++;
		  p2 += n2;
		}
	      *p3++ = sum;
	    }
	  p1s = p1sn;
	}
    }
  }  



DenseMatrix operator+ (const DenseMatrix & m1, const DenseMatrix & m2)
  {
  DenseMatrix temp (m1.Height(), m1.Width());
  int i, j;

  if (m1.Width() != m2.Width() || m1.Height() != m2.Height())
    {
    (*myerr) << "BaseMatrix :: operator+: Matrix Size does not fit" << endl;
    }
  else if (temp.Height() != m1.Height())
    {
    (*myerr) << "BaseMatrix :: operator+: temp not allocated" << endl;
    }
  else
    {
    for (i = 1; i <= m1.Height(); i++)
      for (j = 1; j <= m1.Width(); j++)
        {
        temp.Set(i, j, m1.Get(i, j) + m2.Get(i, j));
        }
    }
  return temp;
  }




void Transpose (const DenseMatrix & m1, DenseMatrix & m2)
{
  int w = m1.Width();
  int h = m1.Height();
  int i, j;

  m2.SetSize (w, h);

  double * pm2 = &m2.Elem(1, 1);
  for (j = 1; j <= w; j++)
    {
      const double * pm1 = &m1.Get(1, j);
      for (i = 1; i <= h; i++)
	{
	  *pm2 = *pm1;
	  pm2 ++;
	  pm1 += w;
	}
    }
}


/*
void DenseMatrix :: Mult (const Vector & v, Vector & prod) const
  {
  double sum, val;
  const double * mp, * sp;
  double * dp;
  // const Vector & v = bv.CastToVector();
  // Vector & prod = bprod.CastToVector();
  

  int n = Height();
  int m = Width();

  if (prod.Size() != n)
    prod.SetSize (n);

#ifdef DEVELOP
  if (!n) 
    {
      cout << "DenseMatrix::Mult  mheight = 0" << endl;
    }
  if (!m) 
    {
      cout << "DenseMatrix::Mult mwidth = 0" << endl;
    }

  if (m != v.Size())
    {
    (*myerr) << "\nMatrix and Vector don't fit" << endl;
    }
  else if (Height() != prod.Size())
    {
    (*myerr) << "Base_Matrix::operator*(Vector): prod vector not ok" << endl;
    }
  else
#endif
    {
      if (Symmetric())
	{
	  int i, j;


	  for (i = 1; i <= n; i++)
	    {
	      sp = &v.Get(1);
	      dp = &prod.Elem(1);
	      mp = &Get(i, 1);

	      val = v.Get(i);
	      sum = Get(i, i) * val;

	      for (j = 1; j < i; ++j, ++mp, ++sp, ++dp)
		{
		  sum += *mp * *sp;
		  *dp += val * *mp;
		}

	      prod.Elem(i) = sum;
	    }
	}
      else
	{
	  mp = data;
	  dp = &prod.Elem(1);
	  for (int i = 1; i <= n; i++)
	    {
	      sum = 0;
	      sp = &v.Get(1);
	      
	      for (int j = 1; j <= m; j++)
		{
		  //        sum += Get(i,j) * v.Get(j);
		  sum += *mp * *sp;
		  mp++;
		  sp++;
		}
	      
	      //      prod.Set (i, sum);
	      *dp = sum;
	      dp++;
	    }
	}
    }
  }
*/

void DenseMatrix :: MultTrans (const Vector & v, Vector & prod) const
{
  // const Vector & v = (const Vector&)bv; // .CastToVector();
  // Vector & prod = (Vector & )bprod;     // .CastToVector();

  /*
  if (Height() != v.Size())
    {
    (*myerr) << "\nMatrix and Vector don't fit" << endl;
    }
  else if (Width() != prod.Size())
    {
    (*myerr) << "Base_Matrix::operator*(Vector): prod vector not ok" << endl;
    }
  else
  */
    {
      int i, j;
      int w = Width(), h = Height();
      if (prod.Size() != w)
	prod.SetSize (w);

      const double * pmat = &Get(1, 1);
      const double * pv = &v.Get(1);

      prod = 0;

      for (i = 1; i <= h; i++)
	{
	  double val = *pv;
	  ++pv;

	  double * pprod = &prod.Elem(1);

	  for (j = w-1; j >= 0; --j, ++pmat, ++pprod)
	    {
	      *pprod += val * *pmat;
	    }
	}
	
      /*
      double sum;

      for (i = 1; i <= Width(); i++)
	{
	  sum = 0;
	  
	  for (int j = 1; j <= Height(); j++)
	    sum += Get(j, i) * v.Get(j);
	  
	  prod.Set (i, sum);
	}
      */
    }
  }


void DenseMatrix :: Residuum (const Vector & x, const Vector & b,
      Vector & res) const
  {
  double sum;
  //   const Vector & x = bx.CastToVector();
  //  const Vector & b = bb.CastToVector();
  //  Vector & res = bres.CastToVector();

  res.SetSize (Height());

  if (Width() != x.Size() || Height() != b.Size())
    {
    (*myerr) << "\nMatrix and Vector don't fit" << endl;
    }
  else if (Height() != res.Size())
    {
    (*myerr) << "Base_Matrix::operator*(Vector): prod vector not ok" << endl;
    }
  else
    {
      int i, j;
      int h = Height(); 
      int w = Width();
      const double * mp = &Get(1, 1);

      for (i = 1; i <= h; i++)
	{
	  sum = b.Get(i);
	  const double * xp = &x.Get(1);

	  for (j = 1; j <= w; ++j, ++mp, ++xp)
	    sum -= *mp * *xp;
	  
	  res.Elem(i) = sum;
	}
    }
  }

#ifdef ABC
double DenseMatrix :: EvaluateBilinearform (const Vector & hx) const
  {
  double sum = 0, hsum;
  // const Vector & hx = x.CastToVector();
  int i, j;

  if (Width() != hx.Size() || Height() != hx.Size())
    {
    (*myerr) << "Matrix::EvaluateBilinearForm: sizes don't fit" << endl;
    }
  else
    {
    for (i = 1; i <= Height(); i++)
      {
      hsum = 0;
      for (j = 1; j <= Height(); j++)
        {
        hsum += Get(i, j) * hx.Get(j);
        }
      sum += hsum * hx.Get(i);
      }
    }

//  testout << "sum = " << sum << endl;
  return sum;
  }


void DenseMatrix :: MultElementMatrix (const ARRAY<int> & pnum, 
      const Vector & hx, Vector & hy)
  {
  int i, j;
  //  const Vector & hx = x.CastToVector();
  //  Vector & hy = y.CastToVector();

  if (Symmetric())
    {
    for (i = 1; i <= Height(); i++)
      {
      for (j = 1; j < i; j++)
        {
	hy.Elem(pnum.Get(i)) += Get(i, j) * hx.Get(pnum.Get(j));
	hy.Elem(pnum.Get(j)) += Get(i, j) * hx.Get(pnum.Get(i));
	}
      hy.Elem(pnum.Get(j)) += Get(i, i) * hx.Get(pnum.Get(i));	
      }
    }
  else
    for (i = 1; i <= Height(); i++)
      for (j = 1; j <= Width(); j++)
	hy.Elem(pnum.Get(i)) += Get(i, j) * hx.Get(pnum.Get(j));
    
  }
  
void DenseMatrix :: MultTransElementMatrix (const ARRAY<int> & pnum, 
      const Vector & hx, Vector & hy)
  {
  int i, j;
  //  const Vector & hx = x.CastToVector();
  //  Vector & hy = y.CastToVector();

  if (Symmetric())
    MultElementMatrix (pnum, hx, hy);
  else
    for (i = 1; i <= Height(); i++)
      for (j = 1; j <= Width(); j++)
	hy.Elem(pnum.Get(i)) += Get(j, i) * hx.Get(pnum.Get(j));
  }
#endif


void DenseMatrix :: Solve (const Vector & v, Vector & sol) const
{
  DenseMatrix temp (*this);
  temp.SolveDestroy (v, sol);
}


void DenseMatrix :: SolveDestroy (const Vector & v, Vector & sol)
  {
  double q;

  if (Width() != Height())
    {
    (*myerr) << "SolveDestroy: Matrix not square";
    return;
    }
  if (Width() != v.Size())
    {
    (*myerr) << "SolveDestroy: Matrix and Vector don't fit";
    return;
    }

  sol = v;
  if (Height() != sol.Size())
    {
    (*myerr) << "SolveDestroy: Solution Vector not ok";
    return;
    }


  if (0 /* Symmetric() */)
    {
      
      // Cholesky factorization

      int i, j, k, n;
      n = Height();
      
      // Cholesky
      
      double x;
      Vector p(n);

      for (i = 1; i <= n; i++)
	for (j = 1; j < i; j++)
	  Elem(j, i) = Get(i, j);
      
      for (i = 1; i <= n; i++)
	{
	  // (*mycout) << "." << flush;
	  for (j = i; j <= n; j++)
	    {
	      x = Get(i, j);

	      const double * pik = &Get(i, 1);
	      const double * pjk = &Get(j, 1);

	      for (k = i-2; k >= 0; --k, ++pik, ++pjk)
		x -= (*pik) * (*pjk);
		  
	      // for (k = i-1; k >= 1; --k)
	      //   x -= Get(j, k) * Get(i, k);

	      if (i == j)
		{
		  if (x <= 0)
		    {
		      cerr << "Matrix indefinite" << endl;
		      return;
		    }
		  
		  p.Elem(i) = 1 / sqrt(x);
		}
	      else
		{
		  Elem(j, i) = x * p.Get(i);
		}
	    }
	}

      for (i = 1; i <= n; i++)
        Elem(i, i) = 1 / p.Get(i);

      // A = L L^t 
      // L stored in left-lower triangle


      sol = v;

      // Solve L sol = sol

      for (i = 1; i <= n; i++)
	{
	  double val = sol.Get(i);

	  const double * pij = &Get(i, 1);
	  const double * solj = &sol.Get(1);

	  for (j = 1; j < i; j++, ++pij, ++solj)
	    val -= *pij * *solj;
	  //	  for (j = 1; j < i; j++)
	  //	    val -= Get(i, j) * sol.Get(j);

	  sol.Elem(i) = val / Get(i, i);
	}

      // Solve L^t sol = sol

      for (i = n; i >= 1; i--)
	{
	  double val = sol.Get(i) / Get(i, i);
	  sol.Elem(i) = val;

	  double * solj = &sol.Elem(1);
	  const double * pij = &Get(i, 1);

	  for (j = 1; j < i; ++j, ++pij, ++solj)
	    *solj -= val * *pij;
	  //	  for (j = 1; j < i; j++)
	  //	    sol.Elem(j) -= Get(i, j) * val;
	}


    }
  else
    {
      //      (*mycout) << "gauss" << endl;
      int i, j, k, n = Height();
      for (i = 1; i <= n; i++)
	{
	  for (j = i+1; j <= n; j++)
	    {
	      q = Get(j,i) / Get(i,i);
	      if (q)
		{
		  const double * pik = &Get(i, i+1);
		  double * pjk = &Elem(j, i+1);

		  for (k = i+1; k <= n; ++k, ++pik, ++pjk)
		    *pjk -= q * *pik;
		  
		  //  for (k = i+1; k <= Height(); k++)
		  //	Elem(j, k) -= q * Get(i,k);


		  sol.Elem(j) -= q * sol.Get(i);
		}
	    }
	}
      
      for (i = n; i >= 1; i--)
	{
	  q = sol.Get(i);
	  for (j = i+1; j <= n; j++)
	      q -= Get(i,j) * sol.Get(j);

	  sol.Elem(i) = q / Get(i,i);
	}
    }
  }


/*
BaseMatrix * DenseMatrix :: Copy () const
  {
  return new DenseMatrix (*this);
  }
*/




ostream & operator<< (ostream & ost, const DenseMatrix & m)
{
  for (int i = 0; i < m.Height(); i++)
    {
      for (int j = 0; j < m.Width(); j++)
	ost << m.Get(i+1,j+1) << " ";
      ost << endl;
    }
  return ost;
}



}
