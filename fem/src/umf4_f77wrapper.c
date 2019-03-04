
/* ----------------------------------------------  */
/* Modified file for use with Elmer.               */
/* Original copyright and usage information below: */


/* ========================================================================== */
/* === umf4_f77wrapper ====================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  CISE Dept,   */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* FORTRAN interface for the C-callable UMFPACK library (double / int version
 * only and double / long versions only).  This is HIGHLY non-portable.  You
 * will need to modify this depending on how your FORTRAN and C compilers
 * behave.  This has been tested in Linux, Sun Solaris, SGI IRIX, and IBM AIX,
 * with various compilers.  It has not been exhaustively tested on all possible
 * combinations of C and FORTRAN compilers.  The long version works on
 * Solaris, SGI IRIX, and IBM AIX when the UMFPACK library is compiled in
 * 64-bit mode.
 *
 * Only a subset of UMFPACK's capabilities are provided.  Refer to the UMFPACK
 * User Guide for details.
 *
 * For some C and FORTRAN compilers, the FORTRAN compiler appends a single
 * underscore ("_") after each routine name.  C doesn't do this, so the
 * translation is made here.  Other FORTRAN compilers treat underscores
 * differently.  For example, a FORTRAN call to a_b gets translated to a call
 * to a_b__ by g77, and to a_b_ by most other FORTRAN compilers.  Thus, the
 * FORTRAN names here do not use underscores.  The xlf compiler in IBM AIX
 * doesn't add an underscore.
 *
 * The matrix A is passed to UMFPACK in compressed column form, with 0-based
 * indices.  In FORTRAN, for an m-by-n matrix A with nz entries, the row
 * indices of the first column (column 1) are in Ai (Ap (1) + 1 ... Ap (2)),
 * with values in Ax (Ap (1) + 1 ... Ap (2)).  The last column (column n) is
 * in Ai (Ap (n) + 1 ... Ap (n+1)) and Ax (Ap (n) + 1 ... Ap (n+1)).  The row
 * indices in Ai are in the range 0 to m-1.  They must be sorted, with no
 * duplicate entries allowed.  Refer to umfpack_di_triplet_to_col for a more
 * flexible format for the input matrix.  The following definitions apply
 * for each of the routines in this file:
 *
 *	integer m, n, Ap (n+1), Ai (nz), symbolic, numeric, filenum, status
 *	double precision Ax (nz), control (20), info (90), x (n), b (n)
 *
 * UMFPACK's status is returned in either a status argument, or in info (1).
 * It is zero if everything is OK, 1 if the matrix is singular (this is a
 * warning, not an error), and negative if an error occurred.  See umfpack.h
 * for more details on the contents of the control and info arrays, and the
 * value of the sys argument.
 *
 * For the Numeric and Symbolic handles, it's probably safe to assume that a
 * FORTRAN integer is sufficient to store a C pointer.  If that doesn't work,
 * try defining numeric and symbolic as integer arrays of size 2, or as
 * integer*8, in the FORTRAN routine that calls these wrapper routines.
 * The latter is required on Solaris, SGI IRIX, and IBM AIX when UMFPACK is
 * compiled in 64-bit mode.
 *
 * If you want to use 64-bit integers, try compiling this file with the -DDLONG
 * compiler option (via "make fortran64").  First modify your Make/Make.include
 * and Make/Make.<arch> files to compile UMFPACK in LP64 mode (see the User
 * Guide for details).  Your FORTRAN code should use integer*8.  See umf4hb64.f
 * for an example.
 *
 * Tested with the following compilers:
 *	* Solaris with cc and f77 from Sun WorkShop 6 update 1
 *	    (32-bit and 64-bit modes)
 *	* SGI Irix with MIPSpro cc and f77 compilers version 7.4
 *	    (32-bit and 64-bit modes)
 *	* Linux with GNU gcc and Intel's icc, and GNU g77 and Intel's
 *	    ifc FORTRAN compiler.  See the comments above about g77 and
 *	    underscores.  Only supports 32-bit mode.
 *	* IBM AIX xlc and xlf compilers.
 *	    (32-bit and 64-bit modes)
 */

#include "../config.h"

#ifdef NULL
#undef NULL
#endif
#define NULL 0

/* Use no name mangling when ISO_C_BINDING is being used */
#ifdef USE_ISO_C_BINDINGS
#undef FC_FUNC
#define FC_FUNC(ARG1,ARG2) ARG1
#undef FC_FUNC_
#define FC_FUNC_(ARG1,ARG2) ARG1
#endif

/* -------------------------------------------------------------------------- */
/* integer type: int or long */
/* -------------------------------------------------------------------------- */

#define Int int
#define UMFPACK_defaults	 umfpack_di_defaults
#define UMFPACK_free_numeric	 umfpack_di_free_numeric
#define UMFPACK_free_symbolic	 umfpack_di_free_symbolic
#define UMFPACK_numeric		 umfpack_di_numeric
#define UMFPACK_report_control	 umfpack_di_report_control
#define UMFPACK_report_info	 umfpack_di_report_info
#define UMFPACK_save_numeric	 umfpack_di_save_numeric
#define UMFPACK_save_symbolic	 umfpack_di_save_symbolic
#define UMFPACK_load_numeric	 umfpack_di_load_numeric
#define UMFPACK_load_symbolic	 umfpack_di_load_symbolic
#define UMFPACK_scale		 umfpack_di_scale
#define UMFPACK_solve		 umfpack_di_solve
#define UMFPACK_symbolic	 umfpack_di_symbolic

#define UMFPACK_IRSTEP 7


/* -------------------------------------------------------------------------- */
/* umf4def: set default control parameters */
/* -------------------------------------------------------------------------- */

/* call umf4def (control) */

void STDCALLBULL FC_FUNC(umf4def,UMF4DEF) (double Control[])
{
#ifdef HAVE_UMFPACK
    UMFPACK_defaults (Control) ;
#endif
}


/* -------------------------------------------------------------------------- */
/* umf4sym: pre-ordering and symbolic factorization */
/* -------------------------------------------------------------------------- */

/* call umf4sym (m, n, Ap, Ai, Ax, symbolic, control, info) */

void STDCALLBULL FC_FUNC(umf4sym,UMF4SYM) (Int *m, Int *n, Int Ap [ ], Int Ai [ ],
    double Ax [ ], void **Symbolic,
    double Control[], double Info[])
{
#ifdef HAVE_UMFPACK
    (void) UMFPACK_symbolic (*m, *n, Ap, Ai, Ax, Symbolic, Control, Info) ;
#endif
}

/* -------------------------------------------------------------------------- */
/* umf4num: numeric factorization */
/* -------------------------------------------------------------------------- */

/* call umf4num (Ap, Ai, Ax, symbolic, numeric, control, info) */

void STDCALLBULL FC_FUNC(umf4num,UMF4NUM) (Int Ap [ ], Int Ai [ ], double Ax [ ],
    void **Symbolic, void **Numeric,
    double Control[], double Info[])
{
#ifdef HAVE_UMFPACK
    (void) UMFPACK_numeric (Ap, Ai, Ax, *Symbolic, Numeric, Control, Info);
#endif
}

/* -------------------------------------------------------------------------- */
/* umf4solr: solve a linear system with iterative refinement */
/* -------------------------------------------------------------------------- */

/* call umf4solr (sys, Ap, Ai, Ax, x, b, numeric, control, info) */

void STDCALLBULL FC_FUNC(umf4solr,UMF4SOLR) (Int *sys, Int Ap [ ], Int Ai [ ], double Ax [ ],
    double x [ ], double b [ ], void **Numeric,
    double Control[], double Info[])
{
#ifdef HAVE_UMFPACK
    (void) UMFPACK_solve (*sys, Ap, Ai, Ax, x, b, *Numeric, Control, Info) ;
#endif
}

/* -------------------------------------------------------------------------- */
/* umf4sol: solve a linear system without iterative refinement */
/* -------------------------------------------------------------------------- */

/* call umf4sol (sys, x, b, numeric, control, info) */

void STDCALLBULL FC_FUNC(umf4sol,UMF4SOL) (Int *sys, double x [ ], double b [ ], void **Numeric,
    double Control[], double Info[])
{
#ifdef HAVE_UMFPACK
    Control [UMFPACK_IRSTEP] = 0 ;
    (void) UMFPACK_solve (*sys, (Int *) NULL, (Int *) NULL, (double *) NULL,
	x, b, *Numeric, Control, Info) ;
#endif
}

/* -------------------------------------------------------------------------- */
/* umf4fnum: free the Numeric object */
/* -------------------------------------------------------------------------- */

/* call umf4fnum (numeric) */

void STDCALLBULL FC_FUNC(umf4fnum,UMF4FNUM) (void **Numeric)
{
#ifdef HAVE_UMFPACK
    UMFPACK_free_numeric (Numeric) ;
#endif
}

/* -------------------------------------------------------------------------- */
/* umf4fsym: free the Symbolic object */
/* -------------------------------------------------------------------------- */

/* call umf4fsym (symbolic) */

void STDCALLBULL FC_FUNC(umf4fsym,UMF4FSYM) (void **Symbolic)
{
#ifdef HAVE_UMFPACK
    UMFPACK_free_symbolic (Symbolic) ;
#endif
}


#undef Int
#undef UMFPACK_defaults
#undef UMFPACK_free_numeric
#undef UMFPACK_free_symbolic
#undef UMFPACK_numeric
#undef UMFPACK_report_control
#undef UMFPACK_report_info
#undef UMFPACK_save_numeric
#undef UMFPACK_save_symbolic
#undef UMFPACK_load_numeric
#undef UMFPACK_load_symbolic
#undef UMFPACK_scale
#undef UMFPACK_solve
#undef UMFPACK_symbolic

#define Int long
#define UMFPACK_defaults	 umfpack_dl_defaults
#define UMFPACK_free_numeric	 umfpack_dl_free_numeric
#define UMFPACK_free_symbolic	 umfpack_dl_free_symbolic
#define UMFPACK_numeric		 umfpack_dl_numeric
#define UMFPACK_report_control	 umfpack_dl_report_control
#define UMFPACK_report_info	 umfpack_dl_report_info
#define UMFPACK_save_numeric	 umfpack_dl_save_numeric
#define UMFPACK_save_symbolic	 umfpack_dl_save_symbolic
#define UMFPACK_load_numeric	 umfpack_dl_load_numeric
#define UMFPACK_load_symbolic	 umfpack_dl_load_symbolic
#define UMFPACK_scale		 umfpack_dl_scale
#define UMFPACK_solve		 umfpack_dl_solve
#define UMFPACK_symbolic	 umfpack_dl_symbolic


/* -------------------------------------------------------------------------- */
/* umf4def: set default control parameters */
/* -------------------------------------------------------------------------- */

/* call umf4def (control) */
void STDCALLBULL FC_FUNC_(umf4_l_def,UMF4_L_DEF) (double Control[])
{
#ifdef HAVE_UMFPACK
    UMFPACK_defaults (Control) ;
#endif
}

/* -------------------------------------------------------------------------- */
/* umf4sym: pre-ordering and symbolic factorization */
/* -------------------------------------------------------------------------- */

/* call umf4sym (m, n, Ap, Ai, Ax, symbolic, control, info) */

void STDCALLBULL FC_FUNC_(umf4_l_sym,UMF4_L_SYM) (Int *m, Int *n, Int Ap [ ], Int Ai [ ],
    double Ax [ ], void **Symbolic,
    double Control[], double Info[])
{
#ifdef HAVE_UMFPACK
    (void) UMFPACK_symbolic (*m, *n, Ap, Ai, Ax, Symbolic, Control, Info) ;
#endif
}

/* -------------------------------------------------------------------------- */
/* umf4num: numeric factorization */
/* -------------------------------------------------------------------------- */

/* call umf4num (Ap, Ai, Ax, symbolic, numeric, control, info) */

void STDCALLBULL FC_FUNC_(umf4_l_num,UMF4_L_NUM) (Int Ap [ ], Int Ai [ ], double Ax [ ],
    void **Symbolic, void **Numeric,
    double Control[], double Info[])
{
#ifdef HAVE_UMFPACK
    (void) UMFPACK_numeric (Ap, Ai, Ax, *Symbolic, Numeric, Control, Info);
#endif
}

/* -------------------------------------------------------------------------- */
/* umf4solr: solve a linear system with iterative refinement */
/* -------------------------------------------------------------------------- */

/* call umf4solr (sys, Ap, Ai, Ax, x, b, numeric, control, info) */

void STDCALLBULL FC_FUNC_(umf4_l_solr,UMF4_L_SOLR) (Int *sys, Int Ap [ ], Int Ai [ ], double Ax [ ],
    double x [ ], double b [ ], void **Numeric,
    double Control[], double Info[])
{
#ifdef HAVE_UMFPACK
    (void) UMFPACK_solve (*sys, Ap, Ai, Ax, x, b, *Numeric, Control, Info) ;
#endif
}

/* -------------------------------------------------------------------------- */
/* umf4sol: solve a linear system without iterative refinement */
/* -------------------------------------------------------------------------- */

/* call umf4sol (sys, x, b, numeric, control, info) */
void STDCALLBULL FC_FUNC_(umf4_l_sol,UMF4_L_SOL) (Int *sys, double x [ ], double b [ ], void **Numeric,
    double Control[], double Info[])
{
#ifdef HAVE_UMFPACK
    Control [UMFPACK_IRSTEP] = 0 ;
    (void) UMFPACK_solve (*sys, (Int *) NULL, (Int *) NULL, (double *) NULL,
	x, b, *Numeric, Control, Info) ;
#endif
}

/* -------------------------------------------------------------------------- */
/* umf4fnum: free the Numeric object */
/* -------------------------------------------------------------------------- */

/* call umf4fnum (numeric) */

void STDCALLBULL FC_FUNC_(umf4_l_fnum,UMF4FNUM) (void **Numeric)
{
#ifdef HAVE_UMFPACK
    UMFPACK_free_numeric (Numeric) ;
#endif
}

/* -------------------------------------------------------------------------- */
/* umf4fsym: free the Symbolic object */
/* -------------------------------------------------------------------------- */

/* call umf4fsym (symbolic) */

void STDCALLBULL FC_FUNC_(umf4_l_fsym,UMF4FSYM) (void **Symbolic)
{
#ifdef HAVE_UMFPACK
    UMFPACK_free_symbolic (Symbolic) ;
#endif
}
