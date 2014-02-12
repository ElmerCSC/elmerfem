#ifndef FILE_AUTOPTR
#define FILE_AUTOPTR

/**************************************************************************/
/* File:   autoptr.hpp                                                    */
/* Author: STL, Joachim Schoeberl                                         */
/* Date:   29. Dec. 02                                                    */
/**************************************************************************/

template <typename T>
class AutoPtr
{
private:
  T * ptr;
public:
  typedef T* pT;
  explicit AutoPtr (T * p = 0)  { ptr = p; }
  ~AutoPtr () { delete ptr; }
  
  T & operator*() const { return *ptr; }
  T* operator->() const { return ptr; }
  T *& Ptr() { return ptr; }
  T * Ptr() const { return ptr; }
  void Reset(T * p = 0) { if (p != ptr) { delete ptr; ptr = p; } }
  operator bool () { return ptr != 0; }
private:
  AutoPtr (AutoPtr &) { ; }
  AutoPtr & operator= (AutoPtr &) { ; }
};

#endif
