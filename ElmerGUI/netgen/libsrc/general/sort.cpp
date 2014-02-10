/**************************************************************************/
/* File:   sort.cc                                                        */
/* Author: Joachim Schoeberl                                              */
/* Date:   07. Jan. 00                                                    */
/**************************************************************************/

/* 
   Sorting
*/


#include <algorithm>
#include <mystdlib.h>
#include <myadt.hpp>

namespace netgen
{

void Sort (const ARRAY<double> & values,
	   ARRAY<int> & order)
{
  int n = values.Size();
  int i, j;

  order.SetSize (n);

  for (i = 1; i <= n; i++)
    order.Elem(i) = i;
  for (i = 1; i <= n-1; i++)
    for (j = 1; j <= n-1; j++)
      if (values.Get(order.Elem(j)) > values.Get(order.Elem(j+1)))
	{
	  Swap (order.Elem(j), order.Elem(j+1));
	}
}


void QickSortRec (const ARRAY<double> & values,
		  ARRAY<int> & order, 
		  int left, int right)
{
  int i, j;
  double midval;

  i = left;
  j = right;
  midval = values.Get(order.Get((i+j)/2));
  
  do
    {
      while (values.Get(order.Get(i)) < midval) i++;
      while (midval < values.Get(order.Get(j))) j--;
      
      if (i <= j)
	{
	  Swap (order.Elem(i), order.Elem(j));
	  i++; j--;
	}
    }
  while (i <= j);
  if (left < j) QickSortRec (values, order, left, j);
  if (i < right) QickSortRec (values, order, i, right);
}

void QickSort (const ARRAY<double> & values,
	       ARRAY<int> & order)
{
  int i, n = values.Size();
  order.SetSize (n);
  for (i = 1; i <= n; i++)
    order.Elem(i) = i;

  QickSortRec (values, order, 1, order.Size());
}
}
