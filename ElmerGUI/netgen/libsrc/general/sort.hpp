#ifndef FILE_SORT
#define FILE_SORT

/**************************************************************************/
/* File:   sort.hh                                                        */
/* Author: Joachim Schoeberl                                              */
/* Date:   07. Jan. 00                                                    */
/**************************************************************************/


// order(i) is sorted index of element i
extern void Sort (const ARRAY<double> & values,
		  ARRAY<int> & order);

extern void QickSort (const ARRAY<double> & values,
		      ARRAY<int> & order);




template <class T>
inline void BubbleSort (int size, T * data)
{
  T hv;
  for (int i = 0; i < size; i++)
    for (int j = i+1; j < size; j++)
      if (data[i] > data[j])
	{
	  hv = data[i];
	  data[i] = data[j];
	  data[j] = hv;
	}
}

template <class T>
inline void BubbleSort (ARRAY<T> & data)
{
  if(data.Size() > 0)
	  BubbleSort (data.Size(), &data[data.Begin()]);
}

#endif
