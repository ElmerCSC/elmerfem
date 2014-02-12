#ifndef FILE_TEMPLATE
#define FILE_TEMPLATE

/**************************************************************************/
/* File:   template.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/*
   templates, global types, defines and variables
*/

///	The following value may be adapted to the hardware !
#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1000000
#endif


// #include <iostream>
/** output stream for testing.
  testout is opened by main */
extern ostream * testout;

/** use instead of cout */
extern ostream * mycout;

/** error output stream */
extern ostream * myerr;

/** Error messages display.
  Error messages are displayed by this function */
extern void MyError (const char * ch);


/** Rings the bell.
  Produces nr beeps. */
extern void MyBeep (int nr = 1);


template <class T>
inline void Swap (T & a, T & b)
{
  T temp = a;
  a = b;
  b = temp;
}

/*
template <class T>
inline void swap (T & a, T & b)
{
  T temp = a;
  a = b;
  b = temp;
}
*/



/**
  INDEX is a typedef for (at least) 4-byte integer
 */
typedef int INDEX;

/**
  BOOL is a typedef for boolean variables
  */
// typedef int BOOL;

typedef int ELIND;
typedef int PIND;


class twoint 
{ 
public: ///
  int i1, i2; ///
  twoint() {};
  ///
  twoint(int ii1, int ii2) {i1 = ii1; i2 = ii2;}
  friend int operator== (const twoint& t1, const twoint& t2);
  ///
  void Swap() {int x = i1; i1 = i2; i2 = x;}
  void Sort() {if (i1 > i2) {Swap();}}
};

inline int operator== (const twoint& t1, const twoint& t2) 
{
  return t1.i1 == t2.i1 && t1.i2 == t2.i2;
}

class threeint 
{ 
public: /// 
  int i1, i2, i3; ///
  threeint() {}; 
  ///
  threeint(int ii1, int ii2, int ii3) {i1 = ii1; i2 = ii2; i3 = ii3;}
};

///
class twodouble
{
public:
  ///
  double d1, d2;
  ///
  twodouble() {d1 = 0; d2 = 0;};
  ///
  twodouble(double id1, double id2) {d1 = id1; d2 = id2;}
  ///
  void Swap() {double x = d1; d1 = d2; d2 = x;}
};

class fourint { public: int i1, i2, i3, i4; fourint() {}; };


///
class INDEX_2;
ostream & operator<<(ostream  & s, const INDEX_2 & i2);


class INDEX_2
{
  ///
  INDEX i[2];

public:
  ///
  INDEX_2 () { }
  ///
  INDEX_2 (INDEX ai1, INDEX ai2)
    { i[0] = ai1; i[1] = ai2; }

  ///
  INDEX_2 (const INDEX_2 & in2)
    { i[0] = in2.i[0]; i[1] = in2.i[1]; }

  ///
  int operator== (const INDEX_2 & in2) const
    { return i[0] == in2.i[0] && i[1] == in2.i[1]; }

  ///


  INDEX_2 Sort ()
  {
    if (i[0] > i[1]) 
      {
	INDEX hi = i[0];
	i[0] = i[1];
	i[1] = hi;
      }
    return *this;
  }

  static INDEX_2 Sort (int i1, int i2)
  {
    if (i1 > i2)
      return INDEX_2 (i2,i1);
    else
      return INDEX_2 (i1,i2);
  }


  ///
  INDEX & I1 () { return i[0]; }
  ///
  INDEX & I2 () { return i[1]; }
  ///
  INDEX & I (int j) { return i[j-1]; }
  ///
  const INDEX & I1 () const { return i[0]; }
  ///
  const INDEX & I2 () const { return i[1]; }
  ///
  const INDEX & I (int j) const { return i[j-1]; }
  ///
  int & operator[] (int j) { return i[j]; }
  ///
  const int & operator[] (int j) const { return i[j]; }
  ///
  friend ostream & operator<<(ostream  & s, const INDEX_2 & i2);
};


///
class INDEX_3
{
  ///
  INDEX i[3];

public:
  ///
  INDEX_3 () { }
  ///
  INDEX_3 (INDEX ai1, INDEX ai2, INDEX ai3)
    { i[0] = ai1; i[1] = ai2; i[2] = ai3; }

  ///
  INDEX_3 (const INDEX_3 & in2)
    { i[0] = in2.i[0]; i[1] = in2.i[1]; i[2] = in2.i[2]; }


  static INDEX_3 Sort (INDEX_3 i3)
  {
    return i3.Sort();
  }

  static INDEX_3 Sort (int i1, int i2, int i3)
  {
    if (i1 > i2) Swap (i1, i2);
    if (i2 > i3) Swap (i2, i3);
    if (i1 > i2) Swap (i1, i2);
    return INDEX_3 (i1, i2, i3);
  }

  INDEX_3 Sort ()
  {
    if (i[0] > i[1]) Swap (i[0], i[1]);
    if (i[1] > i[2]) Swap (i[1], i[2]);
    if (i[0] > i[1]) Swap (i[0], i[1]);
    return *this;
  }

  int operator== (const INDEX_3 & in2) const
    { return i[0] == in2.i[0] && i[1] == in2.i[1] && i[2] == in2.i[2];}

  ///
  INDEX & I1 () { return i[0]; }
  ///
  INDEX & I2 () { return i[1]; }
  ///
  INDEX & I3 () { return i[2]; }
  ///
  INDEX & I (int j) { return i[j-1]; }
  ///
  const INDEX & I1 () const { return i[0]; }
  ///
  const INDEX & I2 () const { return i[1]; }
  ///
  const INDEX & I3 () const { return i[2]; }
  ///
  const INDEX & I (int j) const { return i[j-1]; }
  ///
  int & operator[] (int j) { return i[j]; }
  ///
  const int & operator[] (int j) const { return i[j]; }

  ///
  friend ostream & operator<<(ostream  & s, const INDEX_3 & i3);
};



///
class INDEX_4
{
  ///
  INDEX i[4];

public:
  ///
  INDEX_4 () { }
  ///
  INDEX_4 (INDEX ai1, INDEX ai2, INDEX ai3, INDEX ai4)
    { i[0] = ai1; i[1] = ai2; i[2] = ai3; i[3] = ai4; }

  ///
  INDEX_4 (const INDEX_4 & in2)
    { i[0] = in2.i[0]; i[1] = in2.i[1]; i[2] = in2.i[2]; i[3] = in2.i[3]; }

  ///
  void Sort ();

  ///
  int operator== (const INDEX_4 & in2) const
    { return 
	i[0] == in2.i[0] && i[1] == in2.i[1] && 
	i[2] == in2.i[2] && i[3] == in2.i[3]; }

  ///
  INDEX & I1 () { return i[0]; }
  ///
  INDEX & I2 () { return i[1]; }
  ///
  INDEX & I3 () { return i[2]; }
  ///
  INDEX & I4 () { return i[3]; }
  ///
  INDEX & I (int j) { return i[j-1]; }
  ///
  const INDEX & I1 () const { return i[0]; }
  ///
  const INDEX & I2 () const { return i[1]; }
  ///
  const INDEX & I3 () const { return i[2]; }
  ///
  const INDEX & I4 () const { return i[3]; }
  ///
  const INDEX & I (int j) const { return i[j-1]; }
  ///
  int & operator[] (int j) { return i[j]; }
  ///
  const int & operator[] (int j) const { return i[j]; }

  ///
  friend ostream & operator<<(ostream  & s, const INDEX_4 & i4);
};








/// The sort preserves quads !!!
class INDEX_4Q
{
  ///
  INDEX i[4];

public:
  ///
  INDEX_4Q () { }
  ///
  INDEX_4Q (INDEX ai1, INDEX ai2, INDEX ai3, INDEX ai4)
    { i[0] = ai1; i[1] = ai2; i[2] = ai3; i[3] = ai4; }

  ///
  INDEX_4Q (const INDEX_4Q & in2)
    { i[0] = in2.i[0]; i[1] = in2.i[1]; i[2] = in2.i[2]; i[3] = in2.i[3]; }

  ///
  void Sort ();

  ///
  int operator== (const INDEX_4Q & in2) const
    { return 
	i[0] == in2.i[0] && i[1] == in2.i[1] && 
	i[2] == in2.i[2] && i[3] == in2.i[3]; }

  ///
  INDEX & I1 () { return i[0]; }
  ///
  INDEX & I2 () { return i[1]; }
  ///
  INDEX & I3 () { return i[2]; }
  ///
  INDEX & I4 () { return i[3]; }
  ///
  INDEX & I (int j) { return i[j-1]; }
  ///
  const INDEX & I1 () const { return i[0]; }
  ///
  const INDEX & I2 () const { return i[1]; }
  ///
  const INDEX & I3 () const { return i[2]; }
  ///
  const INDEX & I4 () const { return i[3]; }
  ///
  const INDEX & I (int j) const { return i[j-1]; }
  ///
  friend ostream & operator<<(ostream  & s, const INDEX_4Q & i4);
};












///
template <class T>
inline T min2 (T a, T b)
{
  ///
  return (a < b) ? a : b;
}
///
template <class T>
inline T max2 (T a, T b)
{
  ///
  return (a > b) ? a : b;
}
///
template <class T>
inline T min3 (T a, T b, T c)
{
  ///
  return (a < b) ? (a < c) ? a : c
    : (b < c) ? b : c;
}
///
template <class T>
inline T max3 (T a, T b, T c)
{
  ///
  return (a > b) ? ((a > c) ? a : c)
    : ((b > c) ? b : c);
}

///

///
template <class T>
inline int sgn (T a)
{
  return (a > 0) ? 1 : (   ( a < 0) ? -1 : 0 );
}

///
template <class T>
inline T sqr (const T a)
{
  return a * a; 
}

///
template <class T>
inline T pow3 (const T a)
{
  return a * a * a; 
}



/*
template <class T>
void BubbleSort (int size, T * data);

template <class T>
void MergeSort (int size, T * data, T * help);
*/





#endif
