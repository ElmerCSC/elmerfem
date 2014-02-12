#ifndef FILE_BitArray
#define FILE_BitArray

/**************************************************************************/
/* File:   bitarray.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

#include <limits.h>

/**
   data type BitArray
   
   BitArray is a compressed array of Boolean information. By Set and Clear
   the whole array or one bit can be set or reset, respectively. 
   Test returns the state of the accoring bit.
   No range checking is done.

   index ranges from 0 to size-1
*/
class BitArray
{
  INDEX size;
  unsigned char * data;

public:
  BitArray ();
  ///
  BitArray (INDEX asize);
  ///
  ~BitArray ();

  /// 
  void SetSize (INDEX asize);
  ///
  INDEX Size () const
  {
    return size;
  }

  ///
  void Set ();
  ///
  void Set (INDEX i)
  {
    data[Addr(i)] |= Mask(i);
  }
  
  void Clear ();


  void Clear (INDEX i)
  {
    data[Addr(i)] &= ~Mask(i);
  }

  bool Test (INDEX i) const
  {
    return (data[i / CHAR_BIT] & (char(1) << (i % CHAR_BIT) ) ) ? true : false;
  }

  ///
  void Invert ();
  ///
  void And (const BitArray & ba2);
  ///
  void Or (const BitArray & ba2);
private:
  ///
  inline unsigned char Mask (INDEX i) const
  {
    return char(1) << (i % CHAR_BIT);
  }
  ///
  inline INDEX Addr (INDEX i) const
  {
  return (i / CHAR_BIT);
  }

  ///
  BitArray & operator= (BitArray &);
  ///
  BitArray (const BitArray &);
};



// print bitarray
inline ostream & operator<< (ostream & s, const BitArray & a)
{
  for (int i = 1; i <= a.Size(); i++)
    {
      s << int (a.Test(i));
      if (i % 40 == 0) s << "\n";
    }
  if (a.Size() % 40 != 0) s << "\n";
  return s;
}


/*
inline
INDEX BitArray :: Size () const
  {
  return size;
  }

inline
unsigned char BitArray :: Mask (INDEX i) const
  {
  return char(1) << (i % CHAR_BIT);
  }

inline
INDEX BitArray :: Addr (INDEX i) const
  {
  return (i / CHAR_BIT);
  }
inline
void BitArray :: Set (INDEX i)
  {
  data[Addr(i)] |= Mask(i);
  }

inline
void BitArray :: Clear (INDEX i)
  {
  data[Addr(i)] &= ~Mask(i);
  }


inline
int BitArray :: Test (INDEX i) const
  {
  return (data[i / CHAR_BIT] & (char(1) << (i % CHAR_BIT) ) ) ? 1 : 0;
  }

*/






/**
   data type BitArrayChar
   
   BitArray is an array of Boolean information. By Set and Clear
   the whole array or one bit can be set or reset, respectively. 
   Test returns the state of the accoring bit.
   No range checking is done.
*/
template <int BASE = 1>
class BitArrayChar
{
  ///
  ARRAY<char,BASE> data;

public:
  ///
  BitArrayChar ()
  { ; }
  ///
  BitArrayChar (int asize)
    : data(asize)
  { ; }
  ///
  ~BitArrayChar ()
  { ; }

  ///
  void SetSize (int asize)
  { data.SetSize(asize); }

  ///
  inline int Size () const
  { return data.Size(); }

  ///
  void Set ();
  ///
  inline void Set (int i)
  { data[i] = 1; }
  ///
  void Clear ();
  ///
  inline void Clear (int i)
  { data[i] = 0; }
  ///
  inline int Test (int i) const
  { return data[i]; }
  ///
  void Invert ();
  ///
  void And (const BitArrayChar & ba2);
  ///
  void Or (const BitArrayChar & ba2);
private:
  ///  copy bitarray is not supported
  BitArrayChar & operator= (BitArrayChar &) { return *this; }
  ///  copy bitarray is not supported
  BitArrayChar (const BitArrayChar &) { ; }
};




template <int BASE>
inline ostream & operator<< (ostream & s, const BitArrayChar<BASE> & a)
{
  for (int i = BASE; i < a.Size()+BASE; i++)
    {
      s << a.Test(i);
      if ( (i-BASE) % 40 == 39) s << "\n";
    }
  if (a.Size() % 40 != 0) s << "\n";
  return s;
}


#endif
