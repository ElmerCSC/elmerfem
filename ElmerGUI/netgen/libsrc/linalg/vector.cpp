#ifdef abc
#include <mystdlib.h>
#include <linalg.hpp>
#include <algorithm>

namespace netgen
{

double BaseVector :: shit = 0;

// %FD Constructs a vector of length zero
BaseVector :: BaseVector ()
  {
  length = 0;
  }

// %FD Constructs a vector of given length
BaseVector :: BaseVector (
    INDEX alength  // length of the vector
    )
  {
  length = alength;
  }

// %FD Sets length of the vector, old vector will be destroyed
void
BaseVector :: SetLength (
    INDEX alength        // new length of the vector
    )
  {
  length = alength;
  }

// %FD Changes length of the vector, old values stay valid
void
BaseVector :: ChangeLength (
    INDEX alength        // new length of the vector
    )
  {
  length = alength;
  }



// %FD { Write a vector with the help of the '<<' operator onto a stream }
ostream &    // stream for further use
operator<< (
    ostream & s,            // stream to write vector onto
    const BaseVector & v   // vector to write
    )
  {
  return v.Print (s);
  }


// %FD{ Divides every component of the vector by the scalar c.
//      The function checks for division by zero }
BaseVector &      // result vector
BaseVector :: operator/= (
    double c       // scalar to divide by
    )
  {
  if (c)
    return (*this) *= (1/c);
  else
    {
    (*myerr) << "operator/=: division by zero" << endl;
    return *this;
    }
  }


// %FD Creates a copy of the object
BaseVector *      // pointer to the new vector
BaseVector :: Copy () const
  {
  cerr << "Base_vector::Copy called" << endl << flush;
  return NULL;
  }




void BaseVector :: GetElementVector (const ARRAY<INDEX> & pnum,
				 BaseVector & elvec) const
{
  int i;
  for (i = 1; i <= pnum.Size(); i++)
    elvec(i) = (*this)(pnum.Get(i));
}

void BaseVector :: SetElementVector (const ARRAY<INDEX> & pnum,
				 const BaseVector & elvec)
{
  int i;
  for (i = 1; i <= pnum.Size(); i++)
    (*this)(pnum.Get(i)) = elvec(i);
}


void BaseVector :: AddElementVector (const ARRAY<INDEX> & pnum,
				 const BaseVector & elvec)
{
  int i;
  for (i = 1; i <= pnum.Size(); i++)
    (*this)(pnum.Get(i)) += elvec(i);
}











TempVector :: ~TempVector ()
  {
  delete vec;
  }

TempVector BaseVector :: operator+ (const BaseVector & v2) const
  {
  return (*Copy()) += v2;
  }

TempVector BaseVector :: operator- (const BaseVector & v2) const
  {
  return (*Copy()) -= v2;
  }

TempVector BaseVector :: operator- () const
  {
  return (*Copy()) *= -1;
  }


TempVector operator* (const BaseVector & v1, double scal) 
  {
  return (*v1.Copy()) *= scal;
  }

TempVector operator/ (const BaseVector & v1, double scal) 
  {
  return (*v1.Copy()) /= scal;
  }


TempVector operator* (double scal, const BaseVector & v1)
  {
  return v1 * scal;
  }





BaseVector * TempVector :: Copy () const
  {
  return vec->Copy();
  }










double Vector :: shit = 0;

class clVecpool
{
public:
  ARRAY<double *> vecs;
  ARRAY<INDEX> veclens;

  ~clVecpool();
};

clVecpool :: ~clVecpool()
{
  int i;
  for (i = 1; i <= vecs.Size(); i++)
    delete vecs.Elem(i);
}

static clVecpool vecpool;



static double * NewDouble (INDEX len)
{
  if (len < 10)
    return new double[len];
  else
    {
      int i;
      for (i = 1; i <= vecpool.veclens.Size(); i++)
	if (vecpool.veclens.Get(i) == len)
	  {
	    double * hvec = vecpool.vecs.Get(i);
	    vecpool.vecs.DeleteElement(i);
	    vecpool.veclens.DeleteElement(i);
	    return hvec;
	  }

      return new double[len];
    }
}

static void DeleteDouble (INDEX len, double * dp)
{
  if (len < 10)
    delete [] dp;
  else
    {
      vecpool.vecs.Append (dp);
      vecpool.veclens.Append (len);
    }
}



Vector :: Vector () : BaseVector()
  {
  data = NULL;
  }

Vector :: Vector (INDEX alength) : BaseVector (alength)
  {
  if (length)
    {
      //    data = new double[length];
      data = NewDouble (length);

    if (!data)
      {
      length = 0;
      (*myerr) << "Vector not allocated" << endl;
      }
    }
  else
    data = NULL;
  }


Vector :: Vector (const Vector & v2)
  {
  length = v2.length;

  if (length)
    {
      //    data = new double[length];
      data = NewDouble (length);

    if (data)
      {
      memcpy (data, v2.data, length * sizeof (double));
      }
    else
      {
      length = 0;
      (*myerr) << "Vector::Vector : Vector not allocated" << endl;
      }
    }
  else
    data = NULL;
  }


Vector :: ~Vector ()
{
  //  veclenfile << "~Vector delete: " << length << endl;
  if (data) 
    {
      DeleteDouble (length, data);
      //      delete [] data;
    }

}

void Vector :: SetLength (INDEX alength)
  {
  if (length == alength) return;

  if (data) 
    {
      DeleteDouble (length, data);
      //      delete [] data;
    }
  data = NULL;
  length = alength;

  if (length == 0) return;
  //  data = new double[length];
  data = NewDouble (length);

  if (!data)
    {
    length = 0;
    (*myerr) << "Vector::SetLength: Vector not allocated" << endl;
    }
  }

void Vector :: ChangeLength (INDEX alength)
{
  (*mycout) << "Vector::ChangeLength called" << endl;
  if (length == alength) return;
  
  if (alength == 0)
    {
      //    delete [] data;
      DeleteDouble (length, data);
      length = 0;
      return;
    }
  
  double * olddata = data;

  data = NewDouble (alength);
  //  data = new double[alength];
  if (!data)
    {
    length = 0;
    (*myerr) << "Vector::SetLength: Vector not allocated" << endl;
    delete [] olddata;
    }

  memcpy (data, olddata, min2(alength, length));

  delete [] olddata;
  length = alength;
  }

/// NEW RM
void Vector::SetBlockLength (INDEX /* blength */)
{
  MyError("BaseVector::SetBlockLength was called for a Vector");
}


double & Vector :: operator() (INDEX i)
  {
  if (i >= 1 && i <= length) return Elem(i);
  else (*myerr) << "\nindex " << i << " out of range ("
                                << 1 << "," << Length() << ")\n";
  return shit;
  }

double Vector :: operator() (INDEX i) const
  {
  if (i >= 1 && i <= length) return Get(i);
  else (*myerr) << "\nindex " << i << " out of range ("
                                << 1 << "," << Length() << ")\n" << flush;
  return shit;
  }



double Vector :: SupNorm () const
  {
  double sup = 0;

  for (INDEX i = 1; i <= Length(); i++)
    if (fabs (Get(i)) > sup)
      sup = fabs(Get(i));

  return sup;
  }

double Vector :: L2Norm () const
  {
  double sum = 0;

  for (INDEX i = 1; i <= Length(); i++)
    sum += Get(i) * Get(i);

  return sqrt (sum);
  }

double Vector :: L1Norm () const
  {
  double sum = 0;

  for (INDEX i = 1; i <= Length(); i++)
    sum += fabs (Get(i));

  return sum;
  }

double Vector :: Max () const
  {
  if (!Length()) return 0;
  double m = Get(1);
  for (INDEX i = 2; i <= Length(); i++)
    if (Get(i) > m) m = Get(i);
  return m;
  }

double Vector :: Min () const
  {
  if (!Length()) return 0;
  double m = Get(1);
  for (INDEX i = 2; i <= Length(); i++)
    if (Get(i) < m) m = Get(i);
  return m;
  }


/*
ostream & operator<<(ostream & s, const Vector & v)
  {
  int w = s.width();
  if (v.Length())
    {
    s.width(0);
    s << '(';
    for (INDEX i = 1; i < v.Length(); i++)
      {
      s.width(w);
      s << v.Get(i) << ",";
      if (i % 8 == 0) s << endl << ' ';
      }
    s.width(w);
    s << v.Get(v.Length()) << ')';
    }
  else
    s << "(Vector not allocated)";

  return s;
  }
*/

ostream & Vector :: Print (ostream & s) const
  {
  int w = s.width();
  if (Length())
    {
    s.width(0);
    s << '(';
    for (INDEX i = 1; i < Length(); i++)
      {
      s.width(w);
      s << Get(i) << ",";
      if (i % 8 == 0) s << endl << ' ';
      }
    s.width(w);
    s << Get(Length()) << ')';
    }
  else
    s << "(Vector not allocated)";

  return s;
  }



BaseVector & Vector :: operator+= (const BaseVector & v2)
  {
  const Vector & hv2 = v2.CastToVector();

  if (Length() == hv2.Length())
    for (INDEX i = 1; i <= Length(); i++)
      Elem (i) += hv2.Get(i);
  else
    (*myerr) << "operator+= illegal dimension" << endl;
  return *this;
  }

BaseVector & Vector :: operator-= (const BaseVector & v2)
  {
  const Vector & hv2 = v2.CastToVector();

  if (Length() == hv2.Length())
    for (INDEX i = 1; i <= Length(); i++)
      Elem (i) -= hv2.Get(i);
  else
    (*myerr) << "operator-= illegal dimension" << endl;
  return *this;
  }

BaseVector & Vector :: operator*= (double c)
  {
  for (INDEX i = 1; i <= Length(); i++)
    Elem(i) *= c;
  return *this;
  }



BaseVector & Vector :: Add (double scal, const BaseVector & v2)
  {
  const Vector & hv2 = v2.CastToVector();

  if (Length() == hv2.Length())
    {
    double * p1 = data;
    double * p2 = hv2.data;

    for (INDEX i = Length(); i > 0; i--)
      {
      (*p1) += scal * (*p2);
      p1++; p2++;
      }
    }
  else
    (*myerr) << "Vector::Add: illegal dimension" << endl;
  return *this;
  }

BaseVector & Vector :: Add2 (double scal, const BaseVector & v2,
                             double scal3, const BaseVector & v3)
  {
  const Vector & hv2 = v2.CastToVector();
  const Vector & hv3 = v3.CastToVector();

  if (Length() == hv2.Length())
    {
    double * p1 = data;
    double * p2 = hv2.data;
    double * p3 = hv3.data;

    for (INDEX i = Length(); i > 0; i--)
      {
      (*p1) += (scal * (*p2) + scal3 * (*p3));
      p1++; p2++; p3++;
      }
    }
  else
    (*myerr) << "Vector::Add: illegal dimension" << endl;
  return *this;
  }

BaseVector & Vector :: Set (double scal, const BaseVector & v2)
  {
  const Vector & hv2 = v2.CastToVector();

  if (Length() == hv2.Length())
    {
    double * p1 = data;
    double * p2 = hv2.data;

    for (INDEX i = Length(); i > 0; i--)
      {
      (*p1) = scal * (*p2);
      p1++; p2++;
      }
    }
  else
    (*myerr) << "Vector::Set: illegal dimension" << endl;
  return *this;
  }


BaseVector & Vector :: Set2 (double scal , const BaseVector & v2,
      double scal3, const BaseVector & v3)
{
  const Vector & hv2 = v2.CastToVector();
  const Vector & hv3 = v3.CastToVector();

  if (Length() == hv2.Length() && Length() == hv3.Length())
    {
      double * p1 = data;
      double * p2 = hv2.data;
      double * p3 = hv3.data;
      
      for (INDEX i = Length(); i > 0; i--)
	{
	  (*p1) = scal * (*p2) + scal3 * (*p3);
	  p1++; p2++; p3++;
	}
    }
  else
    (*myerr) << "Vector::Set: illegal dimension" << endl;
  return *this;
}


void Vector :: GetPart (int startpos, BaseVector & v2) const
{
  Vector & hv2 = v2.CastToVector();

  if (Length() >= startpos + v2.Length() - 1)
    {
      const double * p1 = &Get(startpos);
      double * p2 = &hv2.Elem(1);
      
      memcpy (p2, p1, hv2.Length() * sizeof(double));
    }
  else
    MyError ("Vector::GetPart: Vector to short");
}


// NEW RM
void Vector :: SetPart (int startpos, const BaseVector & v2)
{
  const Vector & hv2 = v2.CastToVector();
  INDEX i;
  INDEX n = v2.Length();

  if (Length() >= startpos + n - 1)
    {
      double * p1 = &Elem(startpos);
      const double * p2 = &hv2.Get(1);

      for (i = 1; i <= n; i++)
	{
	  (*p1) = (*p2);
	  p1++;
	  p2++;
	}
    }
  else
    MyError ("Vector::SetPart: Vector to short");
}

void Vector :: AddPart (int startpos, double val, const BaseVector & v2)
{
  const Vector & hv2 = v2.CastToVector();
  INDEX i;
  INDEX n = v2.Length();

  if (Length() >= startpos + n - 1)
    {
      double * p1 = &Elem(startpos);
      const double * p2 = &hv2.Get(1);

      for (i = 1; i <= n; i++)
	{
	  (*p1) += val * (*p2);
	  p1++;
	  p2++;
	}
    }
  else
    MyError ("Vector::AddPart: Vector to short");
}




double Vector :: operator* (const BaseVector & v2) const
  {
  const Vector & hv2 = v2.CastToVector();

  double sum = 0;
  double * p1 = data;
  double * p2 = hv2.data;

  if (Length() == hv2.Length())
    {
    for (INDEX i = Length(); i > 0; i--)
      {
      sum += (*p1) * (*p2);
      p1++; p2++;
      }
    }
  else
    (*myerr) << "Scalarproduct illegal dimension" << endl;
  return sum;
  }

void Vector :: SetRandom ()
  {
  INDEX i;
  for (i = 1; i <= Length(); i++)
    Elem(i) = rand ();

  double l2 = L2Norm();
  if (l2 > 0)
    (*this) /= l2;
    //    Elem(i) = 1.0 / double(i);
    //    Elem(i) = drand48();
  }


/*
TempVector Vector :: operator- () const
  {
  Vector & sum = *(Vector*)Copy();

  if (sum.Length () == Length())
    {
    for (INDEX i = 1; i <= Length(); i++)
      sum.Set (i, Get(i));
    }
  else
    (*myerr) << "operator+ (Vector, Vector): sum.Length() not ok" << endl;
  return sum;
  }
*/

BaseVector & Vector::operator= (const Vector & v2)
  {
  SetLength (v2.Length());

  if (data == v2.data) return *this;
  
  if (v2.Length() == Length())
    memcpy (data, v2.data, sizeof (double) * Length());
  else
    (*myerr) << "Vector::operator=: not allocated" << endl;

  return *this;
  }

BaseVector & Vector::operator= (const BaseVector & v2)
  {
  const Vector & hv2 = v2.CastToVector();

  SetLength (hv2.Length());

  if (data == hv2.data) return *this;
  
  if (hv2.Length() == Length())
    memcpy (data, hv2.data, sizeof (double) * Length());
  else
    (*myerr) << "Vector::operator=: not allocated" << endl;

  return *this;
  }


BaseVector & Vector::operator= (double scal)
  {
  if (!Length()) (*myerr) << "Vector::operator= (double) : data not allocated"
                      << endl;

  for (INDEX i = 1; i <= Length(); i++)
    Set (i, scal);

  return *this;
  }


BaseVector * Vector :: Copy () const
  {
  return new Vector (*this);
  }


void Vector :: Swap (BaseVector & v2)
  {
  Vector & hv2 = v2.CastToVector();
  swap (length, hv2.length);
  swap (data, hv2.data);
  }




void Vector :: GetElementVector (const ARRAY<INDEX> & pnum,
				 BaseVector & elvec) const
{
  int i;
  Vector & helvec = elvec.CastToVector();
  for (i = 1; i <= pnum.Size(); i++)
    helvec.Elem(i) = Get(pnum.Get(i));
}

void Vector :: SetElementVector (const ARRAY<INDEX> & pnum,
				 const BaseVector & elvec)
{
  int i;
  const Vector & helvec = elvec.CastToVector();
  for (i = 1; i <= pnum.Size(); i++)
    Elem(pnum.Get(i)) = helvec.Get(i);
}


void Vector :: AddElementVector (const ARRAY<INDEX> & pnum,
				 const BaseVector & elvec)
{
  int i;
  const Vector & helvec = elvec.CastToVector();
  for (i = 1; i <= pnum.Size(); i++)
    Elem(pnum.Get(i)) += helvec.Get(i);
}
}
#endif
