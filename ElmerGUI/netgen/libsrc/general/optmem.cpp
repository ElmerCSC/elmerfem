/**************************************************************************/
/* File:   optmem.cc                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   04. Apr. 97                                                    */
/**************************************************************************/

/* 
   Abstract data type ARRAY
*/


#include <mystdlib.h>
#include <myadt.hpp>

namespace netgen
{
  //using namespace netgen;

  BlockAllocator :: BlockAllocator (unsigned asize, unsigned ablocks)
    : bablocks (0)
  {
    if (asize < sizeof(void*))
      asize = sizeof(void*);
    size = asize;
    blocks = ablocks;
    freelist = NULL;
  }

  BlockAllocator :: ~BlockAllocator ()
  {
    //for (unsigned i = 0; i < bablocks.Size(); i++)
    for (int i = 0; i < bablocks.Size(); i++)
      delete [] bablocks[i];
  }

  void * BlockAllocator :: Alloc ()
  {
    //  return new char[size];
    if (!freelist)
      {
	// cout << "freelist = " << freelist << endl;
	// cout << "BlockAlloc: " << size*blocks << endl;
	char * hcp = new char [size * blocks];
	bablocks.Append (hcp);
	bablocks.Last() = hcp;
	for (unsigned i = 0; i < blocks-1; i++)
	  *(void**)&(hcp[i * size]) = &(hcp[ (i+1) * size]);
	*(void**)&(hcp[(blocks-1)*size]) = NULL;
	freelist = hcp;
      }

    void * p = freelist;
    freelist = *(void**)freelist;
    return p;
  }

  /*
  void BlockAllocator :: Free (void * p)
  {
    *(void**)p = freelist;
    freelist = p;
  }
  */
}
