#ifndef FILE_AUTODIFF
#define FILE_AUTODIFF

/**************************************************************************/
/* File:   autodiff.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   24. Oct. 02                                                    */
/**************************************************************************/

// Automatic differentiation datatype


/**
   Datatype for automatic differentiation.
   Contains function value and D derivatives. Algebraic
   operations are overloaded by using product-rule etc. etc. 
**/
template <int D, typename SCAL = double>
class AutoDiff
{
  SCAL val;
  SCAL dval[D];
public:

  typedef AutoDiff<D,SCAL> TELEM;
  typedef SCAL TSCAL;


  /// elements are undefined
  AutoDiff  () throw() { }; 
  // { val = 0; for (int i = 0; i < D; i++) dval[i] = 0; }  // !

  /// copy constructor
  AutoDiff  (const AutoDiff & ad2) throw()
  {
    val = ad2.val;
    for (int i = 0; i < D; i++)
      dval[i] = ad2.dval[i];
  }

  /// initial object with constant value
  AutoDiff  (SCAL aval) throw()
  {
    val = aval;
    for (int i = 0; i < D; i++)
      dval[i] = 0;
  }

  /// init object with (val, e_diffindex)
  AutoDiff  (SCAL aval, int diffindex)  throw()
  {
    val = aval;
    for (int i = 0; i < D; i++)
      dval[i] = 0;
    dval[diffindex] = 1;
  }

  /// assign constant value
  AutoDiff & operator= (SCAL aval) throw()
  {
    val = aval;
    for (int i = 0; i < D; i++)
      dval[i] = 0;
    return *this;
  }

  /// returns value
  SCAL Value() const throw() { return val; }

  /// returns partial derivative
  SCAL DValue (int i) const throw() { return dval[i]; }

  /// access value
  SCAL & Value() throw() { return val; }

  /// accesses partial derivative 
  SCAL & DValue (int i) throw() { return dval[i]; }

  /// 
  AutoDiff<D,SCAL> & operator+= (const AutoDiff<D,SCAL> & y) throw()
  {
    val += y.val;
    for (int i = 0; i < D; i++)
      dval[i] += y.dval[i];
    return *this;
  }

  ///
  AutoDiff<D,SCAL> & operator-= (const AutoDiff<D,SCAL> & y) throw()
  {
    val -= y.val;
    for (int i = 0; i < D; i++)
      dval[i] -= y.dval[i];
    return *this;

  }

  ///
  AutoDiff<D,SCAL> & operator*= (const AutoDiff<D,SCAL> & y) throw()
  {
    for (int i = 0; i < D; i++)
      {
	// dval[i] *= y.val;
	// dval[i] += val * y.dval[i];
        dval[i] = dval[i] * y.val + val * y.dval[i];
      }
    val *= y.val;
    return *this;
  }

  ///
  AutoDiff<D,SCAL> & operator*= (const SCAL & y) throw()
  {
    val *= y;
    for (int i = 0; i < D; i++)
      dval[i] *= y;
    return *this;
  }

  ///
  AutoDiff<D,SCAL> & operator/= (const SCAL & y) throw()
  {
    SCAL iy = 1.0 / y;
    val *= iy;
    for (int i = 0; i < D; i++)
      dval[i] *= iy;
    return *this;
  }

  /// 
  bool operator== (SCAL val2) throw()
  {
    return val == val2;
  }

  ///
  bool operator!= (SCAL val2) throw()
  {
    return val != val2;
  }

  ///
  bool operator< (SCAL val2) throw()
  {
    return val < val2;
  }
  
  ///
  bool operator> (SCAL val2) throw()
  {
    return val > val2;
  }
};


//@{  AutoDiff helper functions.

/// prints AutoDiff
template<int D, typename SCAL>
inline ostream & operator<< (ostream & ost, const AutoDiff<D,SCAL> & x)
{
  ost << x.Value() << ", D = ";
  for (int i = 0; i < D; i++)
    ost << x.DValue(i) << " ";
  return ost;
}

/// AutoDiff plus AutoDiff
template<int D, typename SCAL>
inline AutoDiff<D,SCAL> operator+ (const AutoDiff<D,SCAL> & x, const AutoDiff<D,SCAL> & y) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value () = x.Value()+y.Value();
  // AutoDiff<D,SCAL> res(x.Value()+y.Value());
  for (int i = 0; i < D; i++)
    res.DValue(i) = x.DValue(i) + y.DValue(i);
  return res;
}


/// AutoDiff minus AutoDiff
template<int D, typename SCAL>
inline AutoDiff<D,SCAL> operator- (const AutoDiff<D,SCAL> & x, const AutoDiff<D,SCAL> & y) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x.Value()-y.Value();
  // AutoDiff<D,SCAL> res (x.Value()-y.Value());
  for (int i = 0; i < D; i++)
    res.DValue(i) = x.DValue(i) - y.DValue(i);
  return res;
}

/// double plus AutoDiff
template<int D, typename SCAL>
inline AutoDiff<D,SCAL> operator+ (double x, const AutoDiff<D,SCAL> & y) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x+y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = y.DValue(i);
  return res;
}

/// AutoDiff plus double
template<int D, typename SCAL>
inline AutoDiff<D,SCAL> operator+ (const AutoDiff<D,SCAL> & y, double x) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x+y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = y.DValue(i);
  return res;
}


/// minus AutoDiff
template<int D, typename SCAL>
inline AutoDiff<D,SCAL> operator- (const AutoDiff<D,SCAL> & x) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = -x.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = -x.DValue(i);
  return res;
}

/// AutoDiff minus double
template<int D, typename SCAL>
inline AutoDiff<D,SCAL> operator- (const AutoDiff<D,SCAL> & x, double y) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x.Value()-y;
  for (int i = 0; i < D; i++)
    res.DValue(i) = x.DValue(i);
  return res;
}

///
template<int D, typename SCAL>
inline AutoDiff<D,SCAL> operator- (double x, const AutoDiff<D,SCAL> & y) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x-y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = -y.DValue(i);
  return res;
}


/// double times AutoDiff
template<int D, typename SCAL>
inline AutoDiff<D,SCAL> operator* (double x, const AutoDiff<D,SCAL> & y) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x*y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = x*y.DValue(i);
  return res;
}

/// AutoDiff times double
template<int D, typename SCAL>
inline AutoDiff<D,SCAL> operator* (const AutoDiff<D,SCAL> & y, double x) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x*y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = x*y.DValue(i);
  return res;
}

/// AutoDiff times AutoDiff
template<int D, typename SCAL>
inline AutoDiff<D,SCAL> operator* (const AutoDiff<D,SCAL> & x, const AutoDiff<D,SCAL> & y) throw()
{
  AutoDiff<D,SCAL> res;
  SCAL hx = x.Value();
  SCAL hy = y.Value();

  res.Value() = hx*hy;
  for (int i = 0; i < D; i++)
    res.DValue(i) = hx*y.DValue(i) + hy*x.DValue(i);

  return res;
}

/// AutoDiff times AutoDiff
template<int D, typename SCAL>
inline AutoDiff<D,SCAL> sqr (const AutoDiff<D,SCAL> & x) throw()
{
  AutoDiff<D,SCAL> res;
  SCAL hx = x.Value();
  res.Value() = hx*hx;
  hx *= 2;
  for (int i = 0; i < D; i++)
    res.DValue(i) = hx*x.DValue(i);
  return res;
}

/// Inverse of AutoDiff
template<int D, typename SCAL>
inline AutoDiff<D,SCAL> Inv (const AutoDiff<D,SCAL> & x)
{
  AutoDiff<D,SCAL> res(1.0 / x.Value());
  for (int i = 0; i < D; i++)
    res.DValue(i) = -x.DValue(i) / (x.Value() * x.Value());
  return res;
}


/// AutoDiff div AutoDiff
template<int D, typename SCAL>
inline AutoDiff<D,SCAL> operator/ (const AutoDiff<D,SCAL> & x, const AutoDiff<D,SCAL> & y)
{
  return x * Inv (y);
}

/// AutoDiff div double
template<int D, typename SCAL>
inline AutoDiff<D,SCAL> operator/ (const AutoDiff<D,SCAL> & x, double y)
{
  return (1/y) * x;
}

/// double div AutoDiff
template<int D, typename SCAL>
inline AutoDiff<D,SCAL> operator/ (double x, const AutoDiff<D,SCAL> & y)
{
  return x * Inv(y);
}




template<int D, typename SCAL>
inline AutoDiff<D,SCAL> fabs (const AutoDiff<D,SCAL> & x)
{
  double abs = fabs (x.Value());
  AutoDiff<D,SCAL> res( abs );
  if (abs != 0.0)
    for (int i = 0; i < D; i++)
      res.DValue(i) = x.DValue(i) / abs;
  else
    for (int i = 0; i < D; i++)
      res.DValue(i) = 0.0;
  return res;
}

//@}

#endif
