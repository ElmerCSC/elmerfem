#ifndef FILE_GEOMFUNCS
#define FILE_GEOMFUNCS

/* *************************************************************************/
/* File:   geomfuncs.hpp                                                   */
/* Author: Joachim Schoeberl                                               */
/* Date:   20. Jul. 02                                                     */
/* *************************************************************************/


template <int D>
inline double Abs (const Vec<D> & v)
{
  double sum = 0;
  for (int i = 0; i < D; i++)
    sum += v(i) * v(i);
  return sqrt (sum);
}


template <int D>
inline double Abs2 (const Vec<D> & v)
{
  double sum = 0;
  for (int i = 0; i < D; i++)
    sum += v(i) * v(i);
  return sum;
}



template <int D>
inline double Dist (const Point<D> & a, const Point<D> & b)
{
  return Abs (a-b);
}

template <int D>
inline double Dist2 (const Point<D> & a, const Point<D> & b)
{
  return Abs2 (a-b);
}


template <int D>
inline Point<D> Center (const Point<D> & a, const Point<D> & b)
{
  Point<D> res;
  for (int i = 0; i < D; i++)
    res(i) = 0.5 * (a(i) + b(i));
  return res;
}

template <int D>
inline Point<D> Center (const Point<D> & a, const Point<D> & b, const Point<D> & c)
{
  Point<D> res;
  for (int i = 0; i < D; i++)
    res(i) = (1.0/3.0) * (a(i) + b(i) + c(i));
  return res;
}

template <int D>
inline Point<D> Center (const Point<D> & a, const Point<D> & b, const Point<D> & c, const Point<D> & d)
{
  Point<D> res;
  for (int i = 0; i < D; i++)
    res(i) = (1.0/4.0) * (a(i) + b(i) + c(i) + d(i));
  return res;
}



inline Vec<3> Cross (const Vec<3> & v1, const Vec<3> & v2)
{
  return Vec<3> 
    ( v1(1) * v2(2) - v1(2) * v2(1),
      v1(2) * v2(0) - v1(0) * v2(2),
      v1(0) * v2(1) - v1(1) * v2(0) );
}


inline double Determinant (const Vec<3> & col1,
			   const Vec<3> & col2,
			   const Vec<3> & col3)
{
  return
    col1(0) * ( col2(1) * col3(2) - col2(2) * col3(1)) +
    col1(1) * ( col2(2) * col3(0) - col2(0) * col3(2)) +
    col1(2) * ( col2(0) * col3(1) - col2(1) * col3(0));
}



template <>
inline  Vec<2> Vec<2> :: GetNormal () const
{
  return Vec<2> (-x[1], x[0]);
}

template <>
inline  Vec<3> Vec<3> :: GetNormal () const
{
  if (fabs (x[0]) > fabs (x[2]))
    return Vec<3> (-x[1], x[0], 0);
  else
    return Vec<3> (0, x[2], -x[1]);
}



// template <int H, int W>
inline void CalcInverse (const Mat<2,2> & m, Mat<2,2> & inv)
{
  double det = m(0,0) * m(1,1) - m(0,1) * m(1,0);
  if (det == 0) 
    {
      inv = 0;
      return;
    }

  double idet = 1.0 / det;
  inv(0,0) =  idet * m(1,1);
  inv(0,1) = -idet * m(0,1);
  inv(1,0) = -idet * m(1,0);
  inv(1,1) =  idet * m(0,0);
}

void CalcInverse (const Mat<3,3> & m, Mat<3,3> & inv);

inline void CalcInverse (const Mat<2,3> & m, Mat<3,2> & inv)
{
  Mat<2,2> a = m * Trans (m);
  Mat<2,2> ainv;
  CalcInverse (a, ainv);
  inv = Trans (m) * ainv;
}

void CalcInverse (const Mat<3,2> & m, Mat<2,3> & inv);

inline void CalcInverse (const Mat<3,2> & m, Mat<2,3> & inv)
{
  Mat<2,2> a = Trans (m) * m;
  Mat<2,2> ainv;
  CalcInverse (a, ainv);
  inv = ainv * Trans (m);
}


double Det (const Mat<2,2> & m);
double Det (const Mat<3,3> & m);

// eigenvalues of a symmetric matrix
void EigenValues (const Mat<3,3> & m, Vec<3> & ev);
void EigenValues (const Mat<2,2> & m, Vec<3> & ev);

#endif
