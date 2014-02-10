/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  CISE Dept,   */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

GLOBAL Int UMF_kernel
(
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
#ifdef COMPLEX
    const double Az [ ],
#endif
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
) ;
