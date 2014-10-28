
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */


#include "../config.h"
#ifdef HAVE_SUPERLU

// #include "slu_ddefs.h"
#include "pdsp_defs.h"
#define HANDLE_SIZE  8

/* kind of integer to hold a pointer.  Use 'long int'
   so it works on 64-bit systems. */
typedef long int fptr;  /* 64 bit */

typedef struct {
    SuperMatrix *L;
    SuperMatrix *U;
    int *perm_c;
    int *perm_r;
} factors_t;

void
FC_FUNC_(solve_superlu,SOLVE_SUPERLU)(int *iopt, int *nprocs, int *n, 
  int *nnz, int *nrhs, double *values, int *rowind, int *colptr,
      double *b, int *ldb, fptr *f_factors, int *info)

{
/* 
 * This routine can be called from Fortran.
 *
 * iopt (input) int
 *      Specifies the operation:
 *      = 1, performs LU decomposition for the first time
 *      = 2, performs triangular solve
 *      = 3, free all the storage in the end
 *
 * f_factors (input/output) fptr* 
 *      If iopt == 1, it is an output and contains the pointer pointing to
 *                    the structure of the factored matrices.
 *      Otherwise, it it an input.
 *
 */
 
    SuperMatrix A, AC, B;
    SuperMatrix *L, *U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    SCformat *Lstore;
    NCformat *Ustore;
    int i, panel_size, permc_spec, relax, usepr, lwork;
    fact_t fact;
    yes_no_t refact;
    trans_t  trans;
    double   drop_tol = 0.0, *work, diag_pivot_thresh;
    static superlumt_options_t options;
    Gstat_t stat;
    factors_t *LUfactors;

    trans = TRANS;
    fact = NO;
    refact = NO;
    usepr = NO;
    drop_tol=0;
    lwork = 0;
    work = NULL;
    diag_pivot_thresh=1;

    panel_size = sp_ienv(1);
    relax = sp_ienv(2);

    if ( *iopt == 1 ) { /* LU decomposition */

	/* Initialize the statistics variables. */
        StatAlloc(*n, *nprocs, panel_size, relax, &stat);
        StatInit(*n, *nprocs, &stat);
 
	/* Adjust to 0-based indexing */
	for (i = 0; i < *nnz; ++i) --rowind[i];
	for (i = 0; i <= *n; ++i)  --colptr[i];

	dCreate_CompCol_Matrix(&A, *n, *n, *nnz, values, rowind, colptr,
			       SLU_NC, SLU_D, SLU_GE);
	L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
    	perm_r = intMalloc(*n);
    	perm_c = intMalloc(*n);

	/*
	 * Get column permutation vector perm_c[], according to permc_spec:
	 *   permc_spec = 0: natural ordering 
	 *   permc_spec = 1: minimum degree on structure of A'*A
	 *   permc_spec = 2: minimum degree on structure of A'+A
	 *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	 */    	
        permc_spec = 1;
  	get_perm_c(permc_spec, &A, perm_c);
	
        pdgstrf_init(*nprocs, fact, trans, refact, panel_size, relax,
                 diag_pivot_thresh, usepr, drop_tol, perm_c, perm_r,
                 work, lwork, &A, &AC, &options, &stat);

        pdgstrf(&options, &AC, perm_r, L, U, &stat, info);

	if ( *info == 0 ) {
	    Lstore = (SCformat *) L->Store;
	    Ustore = (NCformat *) U->Store;
	    printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
	    printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
	    printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
	} else {
	    printf("dgstrf() error returns INFO= %d\n", *info);
	    if ( *info <= *n ) { /* factorization completes */
	    }
	}
	
	/* Restore to 1-based indexing */
	for (i = 0; i < *nnz; ++i) ++rowind[i];
	for (i = 0; i <= *n; ++i) ++colptr[i];

	/* Save the LU factors in the factors handle */
	LUfactors = (factors_t*) SUPERLU_MALLOC(sizeof(factors_t));
	LUfactors->L = L;
	LUfactors->U = U;
	LUfactors->perm_c = perm_c;
	LUfactors->perm_r = perm_r;
	*f_factors = (fptr) LUfactors;

	/* Free un-wanted storage */
	Destroy_SuperMatrix_Store(&A);
/*	Destroy_CompCol_Permuted(&AC); Instead, as Fabien suggested (untested): */
        pxgstrf_finalize(&options, &AC);
  	StatFree(&stat);

    } else if ( *iopt == 2 ) { /* Triangular solve */
	/* Initialize the statistics variables. */
        StatAlloc(*n, *nprocs, panel_size, relax, &stat);
        StatInit(*n, *nprocs, &stat);

	/* Extract the LU factors in the factors handle */
	LUfactors = (factors_t*) *f_factors;
	L = LUfactors->L;
	U = LUfactors->U;
	perm_c = LUfactors->perm_c;
	perm_r = LUfactors->perm_r;

	dCreate_Dense_Matrix(&B, *n, *nrhs, b, *ldb, SLU_DN, SLU_D, SLU_GE);

        /* Solve the system A*X=B, overwriting B with X. */
        dgstrs (trans, L, U, perm_r, perm_c, &B, &stat, info);

	Destroy_SuperMatrix_Store(&B);
        StatFree(&stat);

    } else if ( *iopt == 3 ) { /* Free storage */
	/* Free the LU factors in the factors handle */

	LUfactors = (factors_t*)*f_factors;
	SUPERLU_FREE (LUfactors->perm_r);
	SUPERLU_FREE (LUfactors->perm_c);
/*
	Destroy_CompCol_Matrix(LUfactors->U);
	Destroy_SuperNode_Matrix(LUfactors->L);
Instead,  as Fabien suggested (untested): */
        Destroy_SuperNode_SCP(LUfactors->L);
        Destroy_CompCol_NCP(LUfactors->U);
        SUPERLU_FREE (LUfactors->L);
        SUPERLU_FREE (LUfactors->U);
	SUPERLU_FREE (LUfactors);
    } else {
	fprintf(stderr,"Invalid iopt=%d passed to c_fortran_dgssv()\n",*iopt);
	exit(-1);
    }
}
#endif
