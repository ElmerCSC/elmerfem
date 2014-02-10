#ifndef FILE_MOVEABLEMEM
#define FILE_MOVEABLEMEM

/**************************************************************************/
/* File:   moveablemem.hpp                                                */
/* Author: Joachim Schoeberl                                              */
/* Date:   12. Feb. 2003                                                  */
/**************************************************************************/


extern NgMutex mem_mutex;

class BaseMoveableMem
{
public:
  static size_t totalsize;
  static size_t used;

private:
  static char * largeblock;
  static BaseMoveableMem *first, *last;

  BaseMoveableMem *prev, *next;
  size_t size, pos;
  char * ptr;
  char * name;

protected:
  BaseMoveableMem (size_t s = 0);
  ~BaseMoveableMem () throw();
  void Alloc (size_t s);
  void ReAlloc (size_t s);
  void MoveTo (size_t newpos);
  void Free () throw(); 
  char * Ptr() { return ptr; }
  char * Ptr() const { return ptr; }
  void Swap (BaseMoveableMem & m2) throw(); 
public:
  void SetName (const char * aname);
  static void Print ();

  friend class BaseDynamicMem;
};




template <typename T>
class MoveableMem : public BaseMoveableMem
{
public:
  MoveableMem (size_t s = 0)
    : BaseMoveableMem (sizeof(T) * s) 
  {
    ;
  }
  void Alloc (size_t s)
  {
    BaseMoveableMem::Alloc (sizeof(T) * s);
  }
  void ReAlloc (size_t s)
  {
    BaseMoveableMem::ReAlloc (sizeof(T) * s);
  }
  void Free ()
  {
    BaseMoveableMem::Free ();
  }

  const T * Ptr() const
  {
    return reinterpret_cast<const T*> (BaseMoveableMem::Ptr());
  }

  T * Ptr()
  {
    return reinterpret_cast<T*> (BaseMoveableMem::Ptr());
  }

  operator T* () const
  {
    return reinterpret_cast<T*> (BaseMoveableMem::Ptr());
  }

  operator T* () 
  {
    return reinterpret_cast<T*> (BaseMoveableMem::Ptr());
  }

  void Swap (MoveableMem<T> & m2)
  {
    BaseMoveableMem::Swap (m2);
  }
protected:
  MoveableMem (const MoveableMem & m) { ; }
  MoveableMem & operator= (const MoveableMem & m) { ; }
};

#endif
