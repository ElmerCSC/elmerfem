#ifndef FILE_DYNAMICMEM
#define FILE_DYNAMICMEM

/**************************************************************************/
/* File:   dynamicmem.hpp                                                 */
/* Author: Joachim Schoeberl                                              */
/* Date:   12. Feb. 2003                                                  */
/**************************************************************************/




class BaseDynamicMem
{
private:
  static BaseDynamicMem *first, *last;

  BaseDynamicMem *prev, *next;
  size_t size;
  char * ptr;
  char * name;

protected:
  BaseDynamicMem ();
  ~BaseDynamicMem ();
  void Alloc (size_t s);
  void ReAlloc (size_t s);
  void Free ();
  char * Ptr() { return ptr; }
  const char * Ptr() const { return ptr; }
  void Swap (BaseDynamicMem & m2);
public:
  void SetName (const char * aname);
  static void Print ();
  static void GetUsed (int nr, char * ch);
};


template <typename T>
class DynamicMem : public BaseDynamicMem
{
public:
  DynamicMem ()
    : BaseDynamicMem () 
  {
    ;
  }
  DynamicMem (size_t s)
    : BaseDynamicMem () 
  {
    Alloc (s);
  }
  void Alloc (size_t s)
  {
    BaseDynamicMem::Alloc (sizeof(T) * s);
  }
  void ReAlloc (size_t s)
  {
    BaseDynamicMem::ReAlloc (sizeof(T) * s);
  }
  void Free ()
  {
    BaseDynamicMem::Free ();
  }

  const T * Ptr() const
  {
    return reinterpret_cast<const T*> (BaseDynamicMem::Ptr());
  }

  T * Ptr()
  {
    return reinterpret_cast<T*> (BaseDynamicMem::Ptr());
  }

  operator const T* () const
  {
    return reinterpret_cast<const T*> (BaseDynamicMem::Ptr());
  }

  operator T* () 
  {
    return reinterpret_cast<T*> (BaseDynamicMem::Ptr());
  }

  void Swap (DynamicMem<T> & m2)
  {
    BaseDynamicMem::Swap (m2);
  }
protected:
  DynamicMem (const DynamicMem & m);
  DynamicMem & operator= (const DynamicMem & m);
};

#endif
