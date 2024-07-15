/* ========================================================================= */
/* === AMD_preprocess ====================================================== */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD Version 1.1 (Jan. 21, 2004), Copyright (c) 2004 by Timothy A. Davis,  */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README for License.         */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.           */
/* web: http://www.cise.ufl.edu/research/sparse/amd                          */
/* ------------------------------------------------------------------------- */

/* Sorts, removes duplicate entries, and transposes from the nonzero pattern of
 * a column-form matrix A, to obtain the matrix R.
 * See amd.h for a complete description of AMD_preprocess 
 */

#include "amd_internal.h"

GLOBAL Int AMD_preprocess   /* returns AMD_OK if input is OK, AMD_INVALID
			     * if the matrix is invalid, or AMD_OUT_OF_MEMORY
			     * if out of memory for the 2n workspace. */
(
    Int n,		/* input matrix: A is n-by-n */
    const Int Ap [ ],	/* size n+1 */
    const Int Ai [ ],	/* size nz = Ap [n] */

    /* output matrix R: */
    Int Rp [ ],		/* size n+1 */
    Int Ri [ ]		/* size nz (or less, if duplicates present) */
)
{
    /* --------------------------------------------------------------------- */
    /* local variables */
    /* --------------------------------------------------------------------- */

    Int *Flag, *W ;

    /* --------------------------------------------------------------------- */
    /* check inputs (note: fewer restrictions than AMD_order) */
    /* --------------------------------------------------------------------- */

    if (!AMD_preprocess_valid (n, Ap, Ai) || !Ri || !Rp)
    {
	return (AMD_INVALID) ;
    }

    /* --------------------------------------------------------------------- */
    /* allocate workspace */
    /* --------------------------------------------------------------------- */

    W = (Int *) ALLOCATE (MAX (n,1) * sizeof (Int)) ;
    if (!W)
    {
	return (AMD_OUT_OF_MEMORY) ;
    }
    Flag = (Int *) ALLOCATE (MAX (n,1) * sizeof (Int)) ;
    if (!Flag)
    {
	FREE (W) ;
	return (AMD_OUT_OF_MEMORY) ;
    }

    /* --------------------------------------------------------------------- */
    /* preprocess the matrix:  sort, remove duplicates, and transpose */
    /* --------------------------------------------------------------------- */

    AMD_wpreprocess (n, Ap, Ai, Rp, Ri, W, Flag) ;

    /* --------------------------------------------------------------------- */
    /* free the workspace */
    /* --------------------------------------------------------------------- */

    FREE (W) ;
    FREE (Flag) ;
    return (AMD_OK) ;
}


/* ========================================================================= */
/* === AMD_wpreprocess ===================================================== */
/* ========================================================================= */

/* The AMD_wpreprocess routine is not user-callable.  It does not check its
 * input for errors or allocate workspace (that is done by the user-callable
 * AMD_preprocess routine).  It does handle the n=0 case. */

GLOBAL void AMD_wpreprocess
(
    Int n,		/* input matrix: A is n-by-n */
    const Int Ap [ ],	/* size n+1 */
    const Int Ai [ ],	/* size nz = Ap [n] */

    /* output matrix R: */
    Int Rp [ ],		/* size n+1 */
    Int Ri [ ],		/* size nz (or less, if duplicates present) */

    Int W [ ],		/* workspace of size n */
    Int Flag [ ]	/* workspace of size n */
)
{

    /* --------------------------------------------------------------------- */
    /* local variables */
    /* --------------------------------------------------------------------- */

    Int i, j, p, p2 ;

    /* --------------------------------------------------------------------- */
    /* count the entries in each row of A (excluding duplicates) */
    /* --------------------------------------------------------------------- */

    for (i = 0 ; i < n ; i++)
    {
	W [i] = 0 ;		/* # of nonzeros in row i (excl duplicates) */
	Flag [i] = EMPTY ;	/* Flag [i] = j if i appears in column j */
    }
    for (j = 0 ; j < n ; j++)
    {
	p2 = Ap [j+1] ;
	for (p = Ap [j] ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    if (Flag [i] != j)
	    {
		/* row index i has not yet appeared in column j */
		W [i]++ ;	    /* one more entry in row i */
		Flag [i] = j ;	    /* flag row index i as appearing in col j*/
	    }
	}
    }

    /* --------------------------------------------------------------------- */
    /* compute the row pointers for R */
    /* --------------------------------------------------------------------- */

    Rp [0] = 0 ;
    for (i = 0 ; i < n ; i++)
    {
	Rp [i+1] = Rp [i] + W [i] ;
    }
    for (i = 0 ; i < n ; i++)
    {
	W [i] = Rp [i] ;
	Flag [i] = EMPTY ;
    }

    /* --------------------------------------------------------------------- */
    /* construct the row form matrix R */
    /* --------------------------------------------------------------------- */

    /* R = row form of pattern of A */
    for (j = 0 ; j < n ; j++)
    {
	p2 = Ap [j+1] ;
	for (p = Ap [j] ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    if (Flag [i] != j)
	    {
		/* row index i has not yet appeared in column j */
		Ri [W [i]++] = j ;  /* put col j in row i */
		Flag [i] = j ;	    /* flag row index i as appearing in col j*/
	    }
	}
    }

#ifndef NDEBUG
    for (j = 0 ; j < n ; j++)
    {
	ASSERT (W [j] == Rp [j+1]) ;
    }
    ASSERT (AMD_valid (n, n, Rp, Ri)) ;
#endif
}


/* ========================================================================= */
/* === AMD_preprocess_valid ================================================ */
/* ========================================================================= */

/* Not user-callable.  Checks a matrix and returns TRUE if it is valid as input
 * to AMD_wpreprocess, FALSE otherwise. */

GLOBAL Int AMD_preprocess_valid
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ]
)
{
    Int i, j, p, nz ;

    if (n < 0 || !Ai || !Ap)
    {
	return (FALSE) ;
    }
    nz = Ap [n] ;
    if (Ap [0] != 0 || nz < 0)
    {
	/* column pointers must start at Ap [0] = 0, and Ap [n] must be >= 0 */
	AMD_DEBUG0 (("column 0 pointer bad or nz < 0\n")) ;
	return (FALSE) ;
    }
    for (j = 0 ; j < n ; j++)
    {
	if (Ap [j] > Ap [j+1])
	{
	    /* column pointers must be ascending */
	    AMD_DEBUG0 (("column "ID" pointer bad\n", j)) ;
	    return (FALSE) ;
	}
    }
    for (p = 0 ; p < nz ; p++)
    {
	i = Ai [p] ;
	AMD_DEBUG3 (("row: "ID"\n", i)) ;
	if (i < 0 || i >= n)
	{
	    /* row index out of range */
	    AMD_DEBUG0 (("index out of range, col "ID" row "ID"\n", j, i)) ;
	    return (FALSE) ;
	}
    }
    return (TRUE) ;
}
