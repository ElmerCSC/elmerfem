/* femsolve.h */
/* This module includes the LU-decomposition algorithms for dense 
   linear matrices. There are also routines for symmetrisizing and 
   normalizing matrices. These are needed to correct the 
   discretization errors of the view factor calculations. */

void SortIndex( int N, Real *Key, int *Ord );
void Symmetrize(Real **vf,int sides);
void Normalize(Real **vf, const Real *b,int sides);
