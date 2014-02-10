#ifndef FILE_LINALG
#define FILE_LINALG

/* *************************************************************************/
/* File:   linalg.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 94                                                    */
/* *************************************************************************/

/* 

   Data types for basic linear algebra
   more data types are found in linalgl.hpp
   
   The basic concepts include the data types 
   
    Vector
    SparseMatrix
    DenseMatrix

*/


#include "../include/myadt.hpp"
namespace netgen
{
#include "vector.hpp"
#include "densemat.hpp"
#include "polynomial.hpp"
}
#endif


