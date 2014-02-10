#ifndef FILE_SPBITA2D
#define FILE_SPBITA2D

/**************************************************************************/
/* File:   spbita2d.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/** 
   Implementation of sparse 2 dimensional bitarray
*/


class SPARSE_BIT_ARRAY_2D
  {
  class linestruct { public: INDEX size; INDEX maxsize; INDEX * col; };

  ///
  linestruct * lines;
  ///
  INDEX height, width;

  public:

  ///
  SPARSE_BIT_ARRAY_2D (INDEX ah = 0, INDEX aw = 0);
  ///
  ~SPARSE_BIT_ARRAY_2D ();

  ///
  void SetSize (INDEX ah, INDEX aw = 0);
  ///
  void DeleteElements ();

  ///
  int Get (INDEX i, INDEX j) const;

  ///
  INDEX Height () const { return height; }
  ///
  INDEX Width () const { return width; }

  ///
  void Set (INDEX i, INDEX j);
  ///
  int Test (INDEX i, INDEX j) const;

  ///
  INDEX BitsInLine (INDEX i) const { return lines[i-1].size; }
  ///
  INDEX GetIndex (INDEX i, INDEX nr) const { return lines[i-1].col[nr-1]; }
  };


#endif
