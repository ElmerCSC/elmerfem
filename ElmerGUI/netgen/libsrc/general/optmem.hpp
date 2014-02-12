#ifndef FILE_OPTMEM
#define FILE_OPTMEM

/**************************************************************************/
/* File:   optmem.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   04. Apr. 97                                                    */
/**************************************************************************/

/** 
    Optimized Memory allocation classes
*/

class BlockAllocator
{
private:
  ///
  unsigned size, blocks;
  ///
  void * freelist;
  ///
  ARRAY<char*> bablocks;
public:
  ///
  BlockAllocator (unsigned asize, unsigned ablocks = 100);
  ///
  ~BlockAllocator ();
  ///

  void * Alloc ();
  /*
  {
    if (!freelist)
      Alloc2();

    void * p = freelist;
    // freelist = *(void**)freelist;
    freelist = *static_cast<void**> (freelist);

    return p;
  }
  */


  ///
  void Free (void * p)
  {
    *(void**)p = freelist;
    freelist = p;
  }
  

private:
  //  void Alloc2 ();
};



#endif
