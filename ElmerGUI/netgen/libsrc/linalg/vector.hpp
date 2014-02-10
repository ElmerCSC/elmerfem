#ifndef FILE_VECTOR
#define FILE_VECTOR

/* *************************************************************************/
/* File:   vector.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 94                                                    */
/* *************************************************************************/


class FlatVector
{
protected:
  int s;
  double *data;
public:
  FlatVector () { ; }
  FlatVector (int as, double * adata)
  { s = as; data = adata; }

  int Size () const
  { return s; }

  FlatVector & operator= (const FlatVector & v) 
  { memcpy (data, v.data, s*sizeof(double)); return *this; }

  FlatVector & operator= (double scal) 
  {
    for (int i = 0; i < s; i++) data[i] = scal; 
    return *this;
  }

  double & operator[] (int i) { return data[i]; }
  const double & operator[] (int i) const { return data[i]; }
  double & operator() (int i) { return data[i]; }
  const double & operator() (int i) const { return data[i]; }

  double & Elem (int i) { return data[i-1]; }
  const double & Get (int i) const { return data[i-1]; }
  void Set (int i, double val) { data[i-1] = val; }

  FlatVector & operator*= (double scal)
  {
    for (int i = 0; i < s; i++) data[i] *= scal;
    return *this;
  }

  FlatVector & Add (double scal, const FlatVector & v2)
  {
    for (int i = 0; i < s; i++) 
      data[i] += scal * v2[i];
    return *this;
  }

  FlatVector & Set (double scal, const FlatVector & v2)
  {
    for (int i = 0; i < s; i++) 
      data[i] = scal * v2[i];
    return *this;
  }

  FlatVector & Set2 (double scal1, const FlatVector & v1,
		 double scal2, const FlatVector & v2)
  {
    for (int i = 0; i < s; i++) 
      data[i] = scal1 * v1[i] + scal2 * v2[i];
    return *this;
  }
  
  double L2Norm() const
  {
    double sum = 0;
    for (int i = 0; i < s; i++)
      sum += data[i] * data[i];
    return sqrt (sum);
  }

  friend double operator* (const FlatVector & v1, const FlatVector & v2);
};



class Vector : public FlatVector
{

public:
  Vector () 
  { s = 0; data = 0; }
  Vector (int as)
  { s = as; data = new double[s]; }
  ~Vector ()
  { delete [] data; }

  Vector & operator= (const FlatVector & v) 
  { memcpy (data, &v.Get(1), s*sizeof(double)); return *this; }

  Vector & operator= (double scal) 
  {
    for (int i = 0; i < s; i++) data[i] = scal; 
    return *this;
  }

  void SetSize (int as)
  {
    if (s != as)
      {
	s = as;
	delete [] data;
	data = new double [s];
      }
  }

};


inline double operator* (const FlatVector & v1, const FlatVector & v2)
{
  double sum = 0;
  for (int i = 0; i < v1.s; i++)
    sum += v1.data[i] * v2.data[i];
  return sum;
}




inline ostream & operator<< (ostream & ost, const FlatVector & v)
{
  for (int i = 0; i < v.Size(); i++)
    ost << " " << setw(7) << v[i];
  return ost;
}



#endif


