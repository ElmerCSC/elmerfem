/**************************************************************************/
/* File:   bitarray.cc                                                    */
/* Autho: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/* 
   data type BitArray
*/

#include <mystdlib.h>
#include <myadt.hpp>


namespace netgen
{
  //using namespace netgen;

  BitArray :: BitArray ()
  {
    size = 0;
    data = NULL;
  }

  BitArray :: BitArray (int asize)
  {
    size = 0;
    data = NULL;
    SetSize (asize);
  }

  BitArray :: ~BitArray ()
  {
    delete [] data;
  }

  void BitArray :: SetSize (int asize)
  {
    if (size == asize) return;
    delete [] data;

    size = asize;
    data = new unsigned char [Addr (size)+1];
  }

  void BitArray :: Set ()
  {
    if (!size) return;
    for (int i = 0; i <= Addr (size); i++)
      data[i] = UCHAR_MAX;
  }

  void BitArray :: Clear ()
  {
    if (!size) return;
    for (int i = 0; i <= Addr (size); i++)
      data[i] = 0;
  }



  void BitArray :: Invert ()
  {
    if (!size) return;
    for (int i = 0; i <= Addr (size); i++)
      data[i] ^= 255;
  }

  void BitArray :: And (const BitArray & ba2)
  {
    if (!size) return;
    for (int i = 0; i <= Addr (size); i++)
      data[i] &= ba2.data[i];
  }


  void BitArray :: Or (const BitArray & ba2)
  {
    if (!size) return;
    for (int i = 0; i <= Addr (size); i++)
      data[i] |= ba2.data[i];
  }











  template <int BASE>
  void BitArrayChar<BASE> :: Set ()
  {
    data = 1;
  }

  template <int BASE>
  void BitArrayChar<BASE> :: Clear ()
  {
    data = 0;
  }


  template <int BASE>
  void BitArrayChar<BASE> :: Invert ()
  {
    for (int i = BASE; i < data.Size()+BASE; i++)
      data[i] = 1 - data[i];
  }

  template <int BASE>
  void BitArrayChar<BASE> :: And (const BitArrayChar & ba2)
  {
    for (int i = BASE; i < data.Size()+BASE; i++)
      data[i] &= ba2.data[i];
  }
  

  template <int BASE>
  void BitArrayChar<BASE> :: Or (const BitArrayChar & ba2)
  {
    for (int i = BASE; i < data.Size()+BASE; i++)
      data[i] |= ba2.data[i];
  }
  

  template class BitArrayChar<0>;
  template class BitArrayChar<1>;
}
