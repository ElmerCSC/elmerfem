/* ========================================================================== */
/* === umfpack_zi_demo ====================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  CISE Dept,   */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
  A demo of UMFPACK:   umfpack_zi_* version.

  First, factor and solve a 5-by-5 system, Ax=b, using default parameters.
  Then solve A'x=b using the factors of A.   Modify one entry (A (1,4) = 0,
  where the row and column indices range from 0 to 4.  The pattern of A
  has not changed (it has explicitly zero entry), so a reanalysis with
  umfpack_zi_symbolic does not need to be done.  Refactorize (with
  umfpack_zi_numeric), and solve Ax=b.  Note that the pivot ordering has
  changed.  Next, change all of the entries in A, but not the pattern.

  Finally, compute C = A', and do the symbolic and numeric factorization of C.
  Factorizing A' can sometimes be better than factorizing A itself (less work
  and memory usage).  Solve C'x=b twice; the solution is the same as the
  solution to Ax=b.

  A note about zero-sized arrays:  UMFPACK uses many user-provided arrays of
  size n (order of the matrix), and of size nz (the number of nonzeros in a
  matrix).  n cannot be zero; UMFPACK does not handle zero-dimensioned arrays.
  However, nz can be zero.  If you attempt to malloc an array of size nz = 0,
  however, malloc will return a null pointer which UMFPACK will report as a
  "missing argument."  Thus, nz1 in this code is set to MAX (nz,1), and
  similarly for lnz and unz.  Lnz can never be zero, however, since L is always
  unit diagonal.
*/

/* -------------------------------------------------------------------------- */
/* definitions */
/* -------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include "umfpack.h"

/* use a cheap approximate absolute value for complex numbers: */
#define ABS(x,z) ((x) >= 0 ? (x) : -(x)) + ((z) >= 0 ? (z) : -(z))

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

/* -------------------------------------------------------------------------- */
/* triplet form of the matrix.  The triplets can be in any order. */
/* -------------------------------------------------------------------------- */

static int    n = 5, nz = 12 ;
static int    Arow [ ] = { 0,  4,  1,  1,   2,   2,  0,  1,  2,  3,  4,  4} ;
static int    Acol [ ] = { 0,  4,  0,  2,   1,   2,  1,  4,  3,  2,  1,  2} ;
static double Aval [ ] = {2., 1., 3., 4., -1., -3., 3., 6., 2., 1., 4., 2.} ;
static double Avalz[ ] = {1., .4, .1, .2, -1., -.2, 0., 6., 3., 0., .3, .3} ;
static double b [ ] = {8., 45., -3., 3., 19.}, x [5], r [5] ;
static double bz[ ] = {1., -5., -2., 0., 2.2}, xz[5], rz[5] ;

/* Avalz, bz:  imaginary part of A and b */

/* -------------------------------------------------------------------------- */
/* error: print a message and exit */
/* -------------------------------------------------------------------------- */

static void error
(
    char *message
)
{
    printf ("\n\n====== error: %s =====\n\n", message) ;
    exit (1) ;
}


/* -------------------------------------------------------------------------- */
/* resid: compute the residual, r = Ax-b or r = A'x=b and return maxnorm (r) */
/* A' is the complex conjugate transpose, not the array transpose */
/* -------------------------------------------------------------------------- */

static double resid
(
    int transpose,
    int Ap [ ],
    int Ai [ ],
    double Ax [ ]
    , double Az [ ]
)
{
    int i, j, p ;
    double norm ;

    for (i = 0 ; i < n ; i++)
    {
	r [i] = -b [i] ;
	rz[i] = -bz[i] ;
    }
    if (transpose)
    {
	for (j = 0 ; j < n ; j++)
	{
	    for (p = Ap [j] ; p < Ap [j+1] ; p++)
	    {
		i = Ai [p] ;
		/* complex: r(j) += conj (Aij) * x (i) */
		r [j] += Ax [p] * x [i] ;
		r [j] += Az [p] * xz[i] ;
		rz[j] -= Az [p] * x [i] ;
		rz[j] += Ax [p] * xz[i] ;
	    }
	}
    }
    else
    {
	for (j = 0 ; j < n ; j++)
	{
	    for (p = Ap [j] ; p < Ap [j+1] ; p++)
	    {
		i = Ai [p] ;
		r [i] += Ax [p] * x [j] ;
		r [i] -= Az [p] * xz[j] ;
		rz[i] += Az [p] * x [j] ;
		rz[i] += Ax [p] * xz[j] ;
	    }
	}
    }
    norm = 0. ;
    for (i = 0 ; i < n ; i++)
    {
	norm = MAX (ABS (r [i], rz [i]), norm) ;
    }
    return (norm) ;
}


/* -------------------------------------------------------------------------- */
/* main program */
/* -------------------------------------------------------------------------- */

int main (int argc, char **argv)
{
    double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL], *Ax, *Cx, *Lx, *Ux,
	*W, t [2], *Dx, rnorm, *Rb, *y, *Rs ;
    double *Az, *Lz, *Uz, *Dz, *Cz, *Rbz, *yz ;
    int *Ap, *Ai, *Cp, *Ci, row, col, p, lnz, unz, nr, nc, *Lp, *Li, *Ui, *Up,
	*P, *Q, *Lj, i, j, k, anz, nfr, nchains, *Qinit, fnpiv, lnz1, unz1, nz1,
	status, *Front_npivcol, *Front_parent, *Chain_start, *Wi, *Pinit, n1,
	*Chain_maxrows, *Chain_maxcols, *Front_1strow, *Front_leftmostdesc,
	nzud, do_recip ;
    void *Symbolic, *Numeric ;

    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

    umfpack_tic (t) ;

    printf ("\n%s demo: _zi_ version\n", UMFPACK_VERSION) ;

    /* get the default control parameters */
    umfpack_zi_defaults (Control) ;

    /* change the default print level for this demo */
    /* (otherwise, nothing will print) */
    Control [UMFPACK_PRL] = 6 ;

    /* print the license agreement */
    umfpack_zi_report_status (Control, UMFPACK_OK) ;
    Control [UMFPACK_PRL] = 5 ;

    /* print the control parameters */
    umfpack_zi_report_control (Control) ;

    /* ---------------------------------------------------------------------- */
    /* print A and b, and convert A to column-form */
    /* ---------------------------------------------------------------------- */

    /* print the right-hand-side */
    printf ("\nb: ") ;
    (void) umfpack_zi_report_vector (n, b, bz, Control) ;

    /* print the triplet form of the matrix */
    printf ("\nA: ") ;
    (void) umfpack_zi_report_triplet (n, n, nz, Arow, Acol, Aval, Avalz,
	Control) ;

    /* convert to column form */
    nz1 = MAX (nz,1) ;	/* ensure arrays are not of size zero. */
    Ap = (int *) malloc ((n+1) * sizeof (int)) ;
    Ai = (int *) malloc (nz1 * sizeof (int)) ;
    Ax = (double *) malloc (nz1 * sizeof (double)) ;
    Az = (double *) malloc (nz1 * sizeof (double)) ;
    if (!Ap || !Ai || !Ax || !Az)
    {
	error ("out of memory") ;
    }

    status = umfpack_zi_triplet_to_col (n, n, nz, Arow, Acol, Aval, Avalz,
	Ap, Ai, Ax, Az, (int *) NULL) ;

    if (status < 0)
    {
	umfpack_zi_report_status (Control, status) ;
	error ("umfpack_zi_triplet_to_col failed") ;
    }

    /* print the column-form of A */
    printf ("\nA: ") ;
    (void) umfpack_zi_report_matrix (n, n, Ap, Ai, Ax, Az, 1, Control) ;

    /* ---------------------------------------------------------------------- */
    /* symbolic factorization */
    /* ---------------------------------------------------------------------- */

    status = umfpack_zi_symbolic (n, n, Ap, Ai, Ax, Az, &Symbolic,
	Control, Info) ;
    if (status < 0)
    {
	umfpack_zi_report_info (Control, Info) ;
	umfpack_zi_report_status (Control, status) ;
	error ("umfpack_zi_symbolic failed") ;
    }

    /* print the symbolic factorization */

    printf ("\nSymbolic factorization of A: ") ;
    (void) umfpack_zi_report_symbolic (Symbolic, Control) ;

    /* ---------------------------------------------------------------------- */
    /* numeric factorization */
    /* ---------------------------------------------------------------------- */

    status = umfpack_zi_numeric (Ap, Ai, Ax, Az, Symbolic, &Numeric,
	Control, Info) ;
    if (status < 0)
    {
	umfpack_zi_report_info (Control, Info) ;
	umfpack_zi_report_status (Control, status) ;
	error ("umfpack_zi_numeric failed") ;
    }

    /* print the numeric factorization */
    printf ("\nNumeric factorization of A: ") ;
    (void) umfpack_zi_report_numeric (Numeric, Control) ;

    /* ---------------------------------------------------------------------- */
    /* solve Ax=b */
    /* ---------------------------------------------------------------------- */

    status = umfpack_zi_solve (UMFPACK_A, Ap, Ai, Ax, Az, x, xz, b, bz,
	Numeric, Control, Info) ;
    umfpack_zi_report_info (Control, Info) ;
    umfpack_zi_report_status (Control, status) ;
    if (status < 0)
    {
	error ("umfpack_zi_solve failed") ;
    }
    printf ("\nx (solution of Ax=b): ") ;
    (void) umfpack_zi_report_vector (n, x, xz, Control) ;
    rnorm = resid (FALSE, Ap, Ai, Ax, Az) ;
    printf ("maxnorm of residual: %g\n\n", rnorm) ;

    /* ---------------------------------------------------------------------- */
    /* compute the determinant */
    /* ---------------------------------------------------------------------- */

    status = umfpack_zi_get_determinant (x, xz, r, Numeric, Info) ;
    umfpack_zi_report_status (Control, status) ;
    if (status < 0)
    {
	error ("umfpack_zi_get_determinant failed") ;
    }
    printf ("determinant: (%g", x [0]) ;
    printf ("+ (%g)i", xz [0]) ; /* complex */
    printf (") * 10^(%g)\n", r [0]) ;

    /* ---------------------------------------------------------------------- */
    /* solve Ax=b, broken down into steps */
    /* ---------------------------------------------------------------------- */

    /* Rb = R*b */
    Rb  = (double *) malloc (n * sizeof (double)) ;
    Rbz = (double *) malloc (n * sizeof (double)) ;
    y   = (double *) malloc (n * sizeof (double)) ;
    yz  = (double *) malloc (n * sizeof (double)) ;
    if (!Rb || !y) error ("out of memory") ;
    if (!Rbz || !yz) error ("out of memory") ;

    status = umfpack_zi_scale (Rb, Rbz, b, bz, Numeric) ;
    if (status < 0) error ("umfpack_zi_scale failed") ;
    /* solve Ly = P*(Rb) */
    status = umfpack_zi_solve (UMFPACK_Pt_L, Ap, Ai, Ax, Az, y, yz, Rb, Rbz,
	Numeric, Control, Info) ;
    if (status < 0) error ("umfpack_zi_solve failed") ;
    /* solve UQ'x=y */
    status = umfpack_zi_solve (UMFPACK_U_Qt, Ap, Ai, Ax, Az, x, xz, y, yz,
	Numeric, Control, Info) ;
    if (status < 0) error ("umfpack_zi_solve failed") ;
    printf ("\nx (solution of Ax=b, solve is split into 3 steps): ") ;
    (void) umfpack_zi_report_vector (n, x, xz, Control) ;
    rnorm = resid (FALSE, Ap, Ai, Ax, Az) ;
    printf ("maxnorm of residual: %g\n\n", rnorm) ;

    free (Rb) ;
    free (Rbz) ;
    free (y) ;
    free (yz) ;

    /* ---------------------------------------------------------------------- */
    /* solve A'x=b */
    /* ---------------------------------------------------------------------- */

    /* note that this is the complex conjugate transpose, A' */
    status = umfpack_zi_solve (UMFPACK_At, Ap, Ai, Ax, Az, x, xz, b, bz,
	Numeric, Control, Info) ;
    umfpack_zi_report_info (Control, Info) ;
    if (status < 0)
    {
	error ("umfpack_zi_solve failed") ;
    }
    printf ("\nx (solution of A'x=b): ") ;
    (void) umfpack_zi_report_vector (n, x, xz, Control) ;
    rnorm = resid (TRUE, Ap, Ai, Ax, Az) ;
    printf ("maxnorm of residual: %g\n\n", rnorm) ;

    /* ---------------------------------------------------------------------- */
    /* modify one numerical value in the column-form of A */
    /* ---------------------------------------------------------------------- */

    /* change A (1,4), look for row index 1 in column 4. */
    row = 1 ;
    col = 4 ;
    for (p = Ap [col] ; p < Ap [col+1] ; p++)
    {
	if (row == Ai [p])
	{
	    printf ("\nchanging A (%d,%d) to zero\n", row, col) ;
	    Ax [p] = 0.0 ;
	    Az [p] = 0.0 ;
	    break ;
	}
    }
    printf ("\nmodified A: ") ;
    (void) umfpack_zi_report_matrix (n, n, Ap, Ai, Ax, Az, 1, Control) ;

    /* ---------------------------------------------------------------------- */
    /* redo the numeric factorization */
    /* ---------------------------------------------------------------------- */

    /* The pattern (Ap and Ai) hasn't changed, so the symbolic factorization */
    /* doesn't have to be redone, no matter how much we change Ax. */

    /* We don't need the Numeric object any more, so free it. */
    umfpack_zi_free_numeric (&Numeric) ;

    /* Note that a memory leak would have occurred if the old Numeric */
    /* had not been free'd with umfpack_zi_free_numeric above. */
    status = umfpack_zi_numeric (Ap, Ai, Ax, Az, Symbolic, &Numeric,
	Control, Info) ;
    if (status < 0)
    {
	umfpack_zi_report_info (Control, Info) ;
	umfpack_zi_report_status (Control, status) ;
	error ("umfpack_zi_numeric failed") ;
    }
    printf ("\nNumeric factorization of modified A: ") ;
    (void) umfpack_zi_report_numeric (Numeric, Control) ;

    /* ---------------------------------------------------------------------- */
    /* solve Ax=b, with the modified A */
    /* ---------------------------------------------------------------------- */

    status = umfpack_zi_solve (UMFPACK_A, Ap, Ai, Ax, Az, x, xz, b, bz,
	Numeric, Control, Info) ;
    umfpack_zi_report_info (Control, Info) ;
    if (status < 0)
    {
	umfpack_zi_report_status (Control, status) ;
	error ("umfpack_zi_solve failed") ;
    }
    printf ("\nx (with modified A): ") ;
    (void) umfpack_zi_report_vector (n, x, xz, Control) ;
    rnorm = resid (FALSE, Ap, Ai, Ax, Az) ;
    printf ("maxnorm of residual: %g\n\n", rnorm) ;

    /* ---------------------------------------------------------------------- */
    /* modify all of the numerical values of A, but not the pattern */
    /* ---------------------------------------------------------------------- */

    for (col = 0 ; col < n ; col++)
    {
	for (p = Ap [col] ; p < Ap [col+1] ; p++)
	{
	    row = Ai [p] ;
	    printf ("changing ") ;
	    /* complex: */ printf ("real part of ") ;
	    printf ("A (%d,%d) from %g", row, col, Ax [p]) ;
	    Ax [p] = Ax [p] + col*10 - row ;
	    printf (" to %g\n", Ax [p]) ;
	}
    }
    printf ("\ncompletely modified A (same pattern): ") ;
    (void) umfpack_zi_report_matrix (n, n, Ap, Ai, Ax, Az, 1, Control) ;

    /* ---------------------------------------------------------------------- */
    /* save the Symbolic object to file, free it, and load it back in */
    /* ---------------------------------------------------------------------- */

    /* use the default filename, "symbolic.umf" */
    printf ("\nSaving symbolic object:\n") ;
    status = umfpack_zi_save_symbolic (Symbolic, (char *) NULL) ;
    if (status < 0)
    {
	umfpack_zi_report_status (Control, status) ;
	error ("umfpack_zi_save_symbolic failed") ;
    }
    printf ("\nFreeing symbolic object:\n") ;
    umfpack_zi_free_symbolic (&Symbolic) ;
    printf ("\nLoading symbolic object:\n") ;
    status = umfpack_zi_load_symbolic (&Symbolic, (char *) NULL) ;
    if (status < 0)
    {
	umfpack_zi_report_status (Control, status) ;
	error ("umfpack_zi_load_symbolic failed") ;
    }
    printf ("\nDone loading symbolic object\n") ;

    /* ---------------------------------------------------------------------- */
    /* redo the numeric factorization */
    /* ---------------------------------------------------------------------- */

    umfpack_zi_free_numeric (&Numeric) ;
    status = umfpack_zi_numeric (Ap, Ai, Ax, Az, Symbolic, &Numeric,
	Control, Info) ;
    if (status < 0)
    {
	umfpack_zi_report_info (Control, Info) ;
	umfpack_zi_report_status (Control, status) ;
	error ("umfpack_zi_numeric failed") ;
    }
    printf ("\nNumeric factorization of completely modified A: ") ;
    (void) umfpack_zi_report_numeric (Numeric, Control) ;

    /* ---------------------------------------------------------------------- */
    /* solve Ax=b, with the modified A */
    /* ---------------------------------------------------------------------- */

    status = umfpack_zi_solve (UMFPACK_A, Ap, Ai, Ax, Az, x, xz, b, bz,
	Numeric, Control, Info) ;
    umfpack_zi_report_info (Control, Info) ;
    if (status < 0)
    {
	umfpack_zi_report_status (Control, status) ;
	error ("umfpack_zi_solve failed") ;
    }
    printf ("\nx (with completely modified A): ") ;
    (void) umfpack_zi_report_vector (n, x, xz, Control) ;
    rnorm = resid (FALSE, Ap, Ai, Ax, Az) ;
    printf ("maxnorm of residual: %g\n\n", rnorm) ;

    /* ---------------------------------------------------------------------- */
    /* free the symbolic and numeric factorization */
    /* ---------------------------------------------------------------------- */

    umfpack_zi_free_symbolic (&Symbolic) ;
    umfpack_zi_free_numeric (&Numeric) ;

    /* ---------------------------------------------------------------------- */
    /* C = transpose of A */
    /* ---------------------------------------------------------------------- */

    Cp = (int *) malloc ((n+1) * sizeof (int)) ;
    Ci = (int *) malloc (nz1 * sizeof (int)) ;
    Cx = (double *) malloc (nz1 * sizeof (double)) ;
    Cz = (double *) malloc (nz1 * sizeof (double)) ;
    if (!Cp || !Ci || !Cx || !Cz)
    {
	error ("out of memory") ;
    }
    status = umfpack_zi_transpose (n, n, Ap, Ai, Ax, Az,
	(int *) NULL, (int *) NULL, Cp, Ci, Cx, Cz, TRUE) ;
    if (status < 0)
    {
	umfpack_zi_report_status (Control, status) ;
	error ("umfpack_zi_transpose failed: ") ;
    }
    printf ("\nC (transpose of A): ") ;
    (void) umfpack_zi_report_matrix (n, n, Cp, Ci, Cx, Cz, 1, Control) ;

    /* ---------------------------------------------------------------------- */
    /* symbolic factorization of C */
    /* ---------------------------------------------------------------------- */

    status = umfpack_zi_symbolic (n, n, Cp, Ci, Cx, Cz, &Symbolic,
	Control, Info) ;
    if (status < 0)
    {
	umfpack_zi_report_info (Control, Info) ;
	umfpack_zi_report_status (Control, status) ;
	error ("umfpack_zi_symbolic failed") ;
    }
    printf ("\nSymbolic factorization of C: ") ;
    (void) umfpack_zi_report_symbolic (Symbolic, Control) ;

    /* ---------------------------------------------------------------------- */
    /* copy the contents of Symbolic into user arrays print them */
    /* ---------------------------------------------------------------------- */

    printf ("\nGet the contents of the Symbolic object for C:\n") ;
    printf ("(compare with umfpack_zi_report_symbolic output, above)\n") ;
    Pinit = (int *) malloc ((n+1) * sizeof (int)) ;
    Qinit = (int *) malloc ((n+1) * sizeof (int)) ;
    Front_npivcol = (int *) malloc ((n+1) * sizeof (int)) ;
    Front_1strow = (int *) malloc ((n+1) * sizeof (int)) ;
    Front_leftmostdesc = (int *) malloc ((n+1) * sizeof (int)) ;
    Front_parent = (int *) malloc ((n+1) * sizeof (int)) ;
    Chain_start = (int *) malloc ((n+1) * sizeof (int)) ;
    Chain_maxrows = (int *) malloc ((n+1) * sizeof (int)) ;
    Chain_maxcols = (int *) malloc ((n+1) * sizeof (int)) ;
    if (!Pinit || !Qinit || !Front_npivcol || !Front_parent || !Chain_start ||
	!Chain_maxrows || !Chain_maxcols || !Front_1strow ||
	!Front_leftmostdesc)
    {
	error ("out of memory") ;
    }

    status = umfpack_zi_get_symbolic (&nr, &nc, &n1, &anz, &nfr, &nchains,
	Pinit, Qinit, Front_npivcol, Front_parent, Front_1strow,
	Front_leftmostdesc, Chain_start, Chain_maxrows, Chain_maxcols,
	Symbolic) ;

    if (status < 0)
    {
	error ("symbolic factorization invalid") ;
    }

    printf ("From the Symbolic object, C is of dimension %d-by-%d\n", nr, nc);
    printf ("   with nz = %d, number of fronts = %d,\n", nz, nfr) ;
    printf ("   number of frontal matrix chains = %d\n", nchains) ;

    printf ("\nPivot columns in each front, and parent of each front:\n") ;
    k = 0 ;
    for (i = 0 ; i < nfr ; i++)
    {
	fnpiv = Front_npivcol [i] ;
	printf ("    Front %d: parent front: %d number of pivot cols: %d\n",
		i, Front_parent [i], fnpiv) ;
	for (j = 0 ; j < fnpiv ; j++)
	{
	    col = Qinit [k] ;
	    printf (
	    "        %d-th pivot column is column %d in original matrix\n",
		k, col) ;
	    k++ ;
	}
    }

    printf ("\nNote that the column ordering, above, will be refined\n") ;
    printf ("in the numeric factorization below.  The assignment of pivot\n") ;
    printf ("columns to frontal matrices will always remain unchanged.\n") ;

    printf ("\nTotal number of pivot columns in frontal matrices: %d\n", k) ;

    printf ("\nFrontal matrix chains:\n") ;
    for (j = 0 ; j < nchains ; j++)
    {
	printf ("   Frontal matrices %d to %d are factorized in a single\n",
	    Chain_start [j], Chain_start [j+1] - 1) ;
	printf ("        working array of size %d-by-%d\n",
	    Chain_maxrows [j], Chain_maxcols [j]) ;
    }

    /* ---------------------------------------------------------------------- */
    /* numeric factorization of C */
    /* ---------------------------------------------------------------------- */

    status = umfpack_zi_numeric (Cp, Ci, Cx, Cz, Symbolic, &Numeric,
	Control, Info) ;
    if (status < 0)
    {
	error ("umfpack_zi_numeric failed") ;
    }
    printf ("\nNumeric factorization of C: ") ;
    (void) umfpack_zi_report_numeric (Numeric, Control) ;

    /* ---------------------------------------------------------------------- */
    /* extract the LU factors of C and print them */
    /* ---------------------------------------------------------------------- */

    if (umfpack_zi_get_lunz (&lnz, &unz, &nr, &nc, &nzud, Numeric) < 0)
    {
	error ("umfpack_zi_get_lunz failed") ;
    }
    /* ensure arrays are not of zero size */
    lnz1 = MAX (lnz,1) ;
    unz1 = MAX (unz,1) ;
    Lp = (int *) malloc ((n+1) * sizeof (int)) ;
    Lj = (int *) malloc (lnz1 * sizeof (int)) ;
    Lx = (double *) malloc (lnz1 * sizeof (double)) ;
    Lz = (double *) malloc (lnz1 * sizeof (double)) ;
    Up = (int *) malloc ((n+1) * sizeof (int)) ;
    Ui = (int *) malloc (unz1 * sizeof (int)) ;
    Ux = (double *) malloc (unz1 * sizeof (double)) ;
    Uz = (double *) malloc (unz1 * sizeof (double)) ;
    P = (int *) malloc (n * sizeof (int)) ;
    Q = (int *) malloc (n * sizeof (int)) ;
    Dx = (double *) NULL ;	/* D vector not requested */
    Dz = (double *) NULL ;
    Rs  = (double *) malloc (n * sizeof (double)) ;
    if (!Lp || !Lj || !Lx || !Lz || !Up || !Ui || !Ux || !Uz || !P || !Q || !Rs)
    {
	error ("out of memory") ;
    }
    status = umfpack_zi_get_numeric (Lp, Lj, Lx, Lz, Up, Ui, Ux, Uz,
	P, Q, Dx, Dz, &do_recip, Rs, Numeric) ;
    if (status < 0)
    {
	error ("umfpack_zi_get_numeric failed") ;
    }

    printf ("\nL (lower triangular factor of C): ") ;
    (void) umfpack_zi_report_matrix (n, n, Lp, Lj, Lx, Lz, 0, Control) ;
    printf ("\nU (upper triangular factor of C): ") ;
    (void) umfpack_zi_report_matrix (n, n, Up, Ui, Ux, Uz, 1, Control) ;
    printf ("\nP: ") ;
    (void) umfpack_zi_report_perm (n, P, Control) ;
    printf ("\nQ: ") ;
    (void) umfpack_zi_report_perm (n, Q, Control) ;
    printf ("\nScale factors: row i of A is to be ") ;
    if (do_recip)
    {
	printf ("multiplied by the ith scale factor\n") ;
    }
    else
    {
	printf ("divided by the ith scale factor\n") ;
    }
    for (i = 0 ; i < n ; i++) printf ("%d: %g\n", i, Rs [i]) ;

    /* ---------------------------------------------------------------------- */
    /* convert L to triplet form and print it */
    /* ---------------------------------------------------------------------- */

    /* Note that L is in row-form, so it is the row indices that are created */
    /* by umfpack_zi_col_to_triplet. */

    printf ("\nConverting L to triplet form, and printing it:\n") ;
    Li = (int *) malloc (lnz1 * sizeof (int)) ;
    if (!Li)
    {
	error ("out of memory") ;
    }
    if (umfpack_zi_col_to_triplet (n, Lp, Li) < 0)
    {
	error ("umfpack_zi_col_to_triplet failed") ;
    }
    printf ("\nL, in triplet form: ") ;
    (void) umfpack_zi_report_triplet (n, n, lnz, Li, Lj, Lx, Lz, Control) ;

    /* ---------------------------------------------------------------------- */
    /* save the Numeric object to file, free it, and load it back in */
    /* ---------------------------------------------------------------------- */

    /* use the default filename, "numeric.umf" */
    printf ("\nSaving numeric object:\n") ;
    status = umfpack_zi_save_numeric (Numeric, (char *) NULL) ;
    if (status < 0)
    {
	umfpack_zi_report_status (Control, status) ;
	error ("umfpack_zi_save_numeric failed") ;
    }
    printf ("\nFreeing numeric object:\n") ;
    umfpack_zi_free_numeric (&Numeric) ;
    printf ("\nLoading numeric object:\n") ;
    status = umfpack_zi_load_numeric (&Numeric, (char *) NULL) ;
    if (status < 0)
    {
	umfpack_zi_report_status (Control, status) ;
	error ("umfpack_zi_load_numeric failed") ;
    }
    printf ("\nDone loading numeric object\n") ;

    /* ---------------------------------------------------------------------- */
    /* solve C'x=b */
    /* ---------------------------------------------------------------------- */

    status = umfpack_zi_solve (UMFPACK_At, Cp, Ci, Cx, Cz, x, xz, b, bz,
	Numeric, Control, Info) ;
    umfpack_zi_report_info (Control, Info) ;
    if (status < 0)
    {
	umfpack_zi_report_status (Control, status) ;
	error ("umfpack_zi_solve failed") ;
    }
    printf ("\nx (solution of C'x=b): ") ;
    (void) umfpack_zi_report_vector (n, x, xz, Control) ;
    rnorm = resid (TRUE, Cp, Ci, Cx, Cz) ;
    printf ("maxnorm of residual: %g\n\n", rnorm) ;

    /* ---------------------------------------------------------------------- */
    /* solve C'x=b again, using umfpack_zi_wsolve instead */
    /* ---------------------------------------------------------------------- */

    printf ("\nSolving C'x=b again, using umfpack_zi_wsolve instead:\n") ;
    Wi = (int *) malloc (n * sizeof (int)) ;
    W = (double *) malloc (10*n * sizeof (double)) ;
    if (!Wi || !W)
    {
	error ("out of memory") ;
    }

    status = umfpack_zi_wsolve (UMFPACK_At, Cp, Ci, Cx, Cz, x, xz, b, bz,
	Numeric, Control, Info, Wi, W) ;
    umfpack_zi_report_info (Control, Info) ;
    if (status < 0)
    {
	umfpack_zi_report_status (Control, status) ;
	error ("umfpack_zi_wsolve failed") ;
    }
    printf ("\nx (solution of C'x=b): ") ;
    (void) umfpack_zi_report_vector (n, x, xz, Control) ;
    rnorm = resid (TRUE, Cp, Ci, Cx, Cz) ;
    printf ("maxnorm of residual: %g\n\n", rnorm) ;

    /* ---------------------------------------------------------------------- */
    /* free everything */
    /* ---------------------------------------------------------------------- */

    /* This is not strictly required since the process is exiting and the */
    /* system will reclaim the memory anyway.  It's useful, though, just as */
    /* a list of what is currently malloc'ed by this program.  Plus, it's */
    /* always a good habit to explicitly free whatever you malloc. */

    free (Ap) ;
    free (Ai) ;
    free (Ax) ;
    free (Az) ;

    free (Cp) ;
    free (Ci) ;
    free (Cx) ;
    free (Cz) ;

    free (Pinit) ;
    free (Qinit) ;
    free (Front_npivcol) ;
    free (Front_1strow) ;
    free (Front_leftmostdesc) ;
    free (Front_parent) ;
    free (Chain_start) ;
    free (Chain_maxrows) ;
    free (Chain_maxcols) ;

    free (Lp) ;
    free (Lj) ;
    free (Lx) ;
    free (Lz) ;

    free (Up) ;
    free (Ui) ;
    free (Ux) ;
    free (Uz) ;

    free (P) ;
    free (Q) ;

    free (Li) ;

    free (Wi) ;
    free (W) ;

    umfpack_zi_free_symbolic (&Symbolic) ;
    umfpack_zi_free_numeric (&Numeric) ;

    /* ---------------------------------------------------------------------- */
    /* print the total time spent in this demo */
    /* ---------------------------------------------------------------------- */

    umfpack_toc (t) ;
    printf ("\numfpack_zi_demo complete.\nTotal time: %5.2f seconds"
	" (CPU time), %5.2f seconds (wallclock time)\n", t [1], t [0]) ;
    return (0) ;
}
