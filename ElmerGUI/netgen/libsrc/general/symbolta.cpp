/**************************************************************************/
/* File:   symbolta.cc                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/* 
   Abstract data type Symbol Table
*/

#include <mystdlib.h>
#include <myadt.hpp>


#ifndef FILE_SYMBOLTABLECC
#define FILE_SYMBOLTABLECC
// necessary for SGI ????


namespace netgen
{
  //using namespace netgen;

  BASE_SYMBOLTABLE :: BASE_SYMBOLTABLE ()
  {
    ;
  }


  BASE_SYMBOLTABLE :: ~BASE_SYMBOLTABLE()
  {
    DelNames();
  }


  void BASE_SYMBOLTABLE :: DelNames()
  {
    for (int i = 0; i < names.Size(); i++)
      delete [] names[i];
    names.SetSize (0);
  }

  int BASE_SYMBOLTABLE :: Index (const char * name) const
  {
    if (!name) return 0;
    for (int i = 0; i < names.Size(); i++)
      if (strcmp (names[i], name) == 0) return i+1;
    return 0;
  }
}

#endif
