/* Simple test helper to print CRS row pointers from C
 * This function is called from Fortran for debugging.
 */
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <permonqps.h>
#include <petscmat.h>
#include <petscsys.h>

typedef struct {
    PetscBool initialized;
    MPI_Comm comm;
    int fcomm;
    PetscInt nrows;
    PetscInt nlocal;
    PetscInt *globaldofs_copy;

    Mat A;
    Vec b, c, x;
    Vec lb_fill, ub_fill;
    QP qp;
    QPS qps;

    IS is_from, is_to;
    Vec out;
    VecScatter scatter;

    PetscBool pattern_cached;
    PetscInt nnz_total;
    PetscInt maxnnz;
    PetscInt *cols0;
} PermonSolveCache;

static PermonSolveCache g_cache = {0};

static PetscErrorCode permon_cache_reset(PermonSolveCache *cache)
{
    PetscFunctionBegin;
    if (!cache) PetscFunctionReturn(0);

    PetscCall(MatDestroy(&cache->A));
    PetscCall(VecDestroy(&cache->b));
    PetscCall(VecDestroy(&cache->c));
    PetscCall(VecDestroy(&cache->x));
    PetscCall(VecDestroy(&cache->lb_fill));
    PetscCall(VecDestroy(&cache->ub_fill));
    PetscCall(QPSDestroy(&cache->qps));
    PetscCall(QPDestroy(&cache->qp));

    PetscCall(VecScatterDestroy(&cache->scatter));
    PetscCall(ISDestroy(&cache->is_from));
    PetscCall(ISDestroy(&cache->is_to));
    PetscCall(VecDestroy(&cache->out));

    PetscCall(PetscFree(cache->cols0));
    PetscCall(PetscFree(cache->globaldofs_copy));

    cache->initialized = PETSC_FALSE;
    cache->comm = MPI_COMM_NULL;
    cache->fcomm = 0;
    cache->nrows = 0;
    cache->nlocal = 0;
    cache->pattern_cached = PETSC_FALSE;
    cache->nnz_total = 0;
    cache->maxnnz = 0;
    PetscFunctionReturn(0);
}

static PetscBool permon_same_layout(const PermonSolveCache *cache, MPI_Comm comm, int fcomm, PetscInt nrows, PetscInt nlocal, const int *globaldofs)
{
    PetscInt i;

    if (!cache->initialized) return PETSC_FALSE;
    if (cache->fcomm != fcomm) return PETSC_FALSE;
    if (cache->comm != comm) return PETSC_FALSE;
    if (cache->nrows != nrows) return PETSC_FALSE;
    if (cache->nlocal != nlocal) return PETSC_FALSE;
    if (!cache->globaldofs_copy) return PETSC_FALSE;

    for (i = 0; i < nrows; ++i) {
        if (cache->globaldofs_copy[i] != (PetscInt)globaldofs[i]) return PETSC_FALSE;
    }
    return PETSC_TRUE;
}

static PetscErrorCode permon_setup_cached_objects(PermonSolveCache *cache, MPI_Comm comm, int fcomm, PetscInt nrows, PetscInt nlocal, const int *globaldofs)
{
    ISLocalToGlobalMapping ltog = NULL;
    PetscInt i;

    PetscFunctionBegin;
    PetscCall(permon_cache_reset(cache));

    cache->comm = comm;
    cache->fcomm = fcomm;
    cache->nrows = nrows;
    cache->nlocal = nlocal;

    PetscCall(PetscMalloc1(nrows, &cache->globaldofs_copy));
    for (i = 0; i < nrows; ++i) cache->globaldofs_copy[i] = (PetscInt)globaldofs[i];

    PetscCall(ISLocalToGlobalMappingCreate(comm, 1, nrows, cache->globaldofs_copy, PETSC_COPY_VALUES, &ltog));

    PetscCall(MatCreate(comm, &cache->A));
    PetscCall(MatSetType(cache->A, MATMPIAIJ));
    PetscCall(MatSetSizes(cache->A, nlocal, nlocal, PETSC_DECIDE, PETSC_DECIDE));
    PetscCall(MatSetFromOptions(cache->A));
    PetscCall(MatSetUp(cache->A));
    PetscCall(MatSetLocalToGlobalMapping(cache->A, ltog, ltog));
    PetscCall(ISLocalToGlobalMappingDestroy(&ltog));

    PetscCall(MatCreateVecs(cache->A, &cache->x, &cache->b));
    PetscCall(VecDuplicate(cache->x, &cache->c));

    PetscCall(VecDuplicate(cache->x, &cache->lb_fill));
    PetscCall(VecDuplicate(cache->x, &cache->ub_fill));
    PetscCall(VecSet(cache->lb_fill, PETSC_NINFINITY));
    PetscCall(VecSet(cache->ub_fill, PETSC_INFINITY));

    PetscCall(QPCreate(comm, &cache->qp));
    PetscCall(QPSetOperator(cache->qp, cache->A));
    PetscCall(QPSetRhs(cache->qp, cache->b));

    PetscCall(QPSCreate(PetscObjectComm((PetscObject)cache->qp), &cache->qps));
    PetscCall(QPSSetQP(cache->qps, cache->qp));

    PetscCall(ISCreateGeneral(comm, nrows, cache->globaldofs_copy, PETSC_COPY_VALUES, &cache->is_from));
    PetscCall(ISCreateStride(PETSC_COMM_SELF, nrows, 0, 1, &cache->is_to));
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, nrows, &cache->out));
    PetscCall(VecScatterCreate(cache->x, cache->is_from, cache->out, cache->is_to, &cache->scatter));

    cache->initialized = PETSC_TRUE;
    PetscFunctionReturn(0);
}

static PetscErrorCode permon_cache_pattern(PermonSolveCache *cache, PetscInt nrows, const int *rows_f, const int *cols_f)
{
    PetscInt i;

    PetscFunctionBegin;
    if (cache->pattern_cached) PetscFunctionReturn(0);

    cache->nnz_total = (PetscInt)rows_f[nrows] - 1;
    if (cache->nnz_total < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Invalid CRS row pointer end value");

    cache->maxnnz = 0;
    for (i = 0; i < nrows; ++i) {
        PetscInt rownnz = (PetscInt)rows_f[i + 1] - (PetscInt)rows_f[i];
        cache->maxnnz = PetscMax(cache->maxnnz, rownnz);
    }

    PetscCall(PetscMalloc1(cache->nnz_total, &cache->cols0));
    for (i = 0; i < cache->nnz_total; ++i) {
        /* Convert Fortran 1-based local columns to C 0-based once per cached layout. */
        cache->cols0[i] = (PetscInt)cols_f[i] - 1;
    }

    cache->pattern_cached = PETSC_TRUE;
    PetscFunctionReturn(0);
}


void mprgp_print_rows(void *cptr, intptr_t addr, int n)
{
    if (cptr == NULL) {
        printf("mprgp_print_rows: received NULL pointer (addr=%" PRIdPTR ")\n", (intptr_t)cptr);
        fflush(stdout);
        return;
    }

    int *rows = (int*) cptr;
    printf("mprgp_print_rows: C pointer address = %p (intptr=%" PRIdPTR ")\n", (void*)rows, (intptr_t)rows);
    int toprint = n;
    for (int i = 0; i < toprint; ++i) {
        printf(" C row %3d: %12d\n", i+1, rows[i]);
    }
    fflush(stdout);
}

void mprgp_print_vector(void *cptr, int n, char *name)
{
    if (cptr == NULL) {
        printf("mprgp_print_vector: received NULL pointer (addr=%" PRIdPTR ")\n", (intptr_t)cptr);
        fflush(stdout);
        return;
    }

    double *vec = (double*) cptr;
    printf("mprgp_print_vector: %s: C pointer address = %p (intptr=%" PRIdPTR ")\n", name, (void*)vec, (intptr_t)vec);
    int toprint = n;
    for (int i = 0; i < toprint; ++i) {
        printf(" %s row %3d: %12.5e\n", name, i+1, vec[i]);
    }
    fflush(stdout);
}

int permon_init(){
    // permonrc - default name for solver options file
    return PermonInitialize(NULL, NULL, NULL, NULL);
}

int permon_set_options_file(const char *options_file, int fcomm){
    MPI_Comm comm = MPI_Comm_f2c(fcomm);

    if (options_file == NULL || options_file[0] == '\0') {
        return 0;
    }

    PetscErrorCode ierr;
    ierr = PetscOptionsInsertFile(comm, NULL, options_file, PETSC_TRUE);

    /* If options were inserted after PetscInitialize, some options
       (e.g. log viewing) expect a log handler to be active. Start
       the default log handler.*/
    PetscCall(PetscLogDefaultBegin());

    return ierr;
}

int permon_finalize(){
    PetscCall(permon_cache_reset(&g_cache));
    return PermonFinalize();
}

// TODO check if the freeing of the arrays is correct
int permon_solve(void *rows_local, void *cols_local, void *vals_local, int nrows, void *b_ptr, void *c_ptr, void *x_ptr, int bound, int *globaldofs, int *owner, int fcomm) {
    Vec       b, c, x;
    Vec       lb = NULL, ub = NULL;
    Mat       A, AT = NULL;
    QP        qp;
    QPS       qps;
    PetscInt  i, rstart, rend, nnz;
    PetscBool converged, viewSol = PETSC_FALSE;
    PetscBool debugInit = PETSC_FALSE, pinInitToBound = PETSC_FALSE, pinInitToBoundAtFirst = PETSC_FALSE;
    PetscBool debugBounds = PETSC_FALSE;
    PetscBool checkSymmetry = PETSC_FALSE, symmetryStrict = PETSC_FALSE, isSymmetric = PETSC_FALSE;
    PetscBool symmetrizeOperator = PETSC_FALSE;
    PetscReal symmetryTol = 1e-12;
    PetscViewer viewer;
    static PetscInt solveCallCounter = 0;
    
    MPI_Comm comm=MPI_Comm_f2c(fcomm);
    
    int *rows_f = (int*)rows_local;
    int *cols_f = (int*)cols_local;

    double *vals = (double*)vals_local;

    // -----------------------------
    // 1. Compute local ownership range
    // -----------------------------
    PetscInt ilower = PETSC_MAX_INT, iupper = -1;
    PetscInt nlocal = 0;


    // Find number of local rows owned by this rank
    for (i = 0; i < nrows; i++) {
        if (owner[i]) {
            if (globaldofs[i] < ilower) ilower = globaldofs[i];
            if (globaldofs[i] > iupper) iupper = globaldofs[i];
            nlocal++;
        }
    }
    if (iupper == -1) { ilower = 0; iupper = -1; }
    
    if (!permon_same_layout(&g_cache, comm, fcomm, nrows, nlocal, globaldofs)) {
        PetscCall(permon_setup_cached_objects(&g_cache, comm, fcomm, nrows, nlocal, globaldofs));
    }

    A = g_cache.A;
    b = g_cache.b;
    c = g_cache.c;
    x = g_cache.x;
    qp = g_cache.qp;
    qps = g_cache.qps;

    PetscCall(permon_cache_pattern(&g_cache, nrows, rows_f, cols_f));

    PetscCall(MatZeroEntries(A));

    for (i = 0; i < nrows; i++) {
        nnz  = rows_f[i+1] - rows_f[i];
        PetscInt iloc = i; /* local row index */

        PetscCall(MatSetValuesLocal(A, 1, &iloc, nnz, &g_cache.cols0[rows_f[i] - 1],
                &vals[rows_f[i] - 1], ADD_VALUES));
    }


    PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

    /* Optional defect-correction style operator modification:
       A <- 0.5 * (A + A^T). This enforces symmetry for PERMON while
       keeping the RHS from the full assembled physics. */
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-permon_symmetrize_operator", &symmetrizeOperator, NULL));
    if (symmetrizeOperator) {
        PetscCall(MatTranspose(A, MAT_INITIAL_MATRIX, &AT));
        PetscCall(MatAXPY(A, 1.0, AT, DIFFERENT_NONZERO_PATTERN));
        PetscCall(MatScale(A, 0.5));
        PetscCall(MatDestroy(&AT));
        PetscCall(PetscPrintf(comm,
            "[permon_solve] Using symmetrized operator A_sym = 0.5*(A + A^T).\n"));
    }

    /* MPRGP/PERMON assumes a symmetric operator; this check verifies the
       assembled matrix and can optionally abort before solve. */
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-permon_check_symmetry", &checkSymmetry, NULL));
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-permon_check_symmetry_tol", &symmetryTol, NULL));
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-permon_check_symmetry_strict", &symmetryStrict, NULL));
    if (checkSymmetry) {
        PetscCall(MatIsSymmetric(A, symmetryTol, &isSymmetric));
        if (!isSymmetric) {
            PetscCall(PetscPrintf(comm,
                "[permon_solve] WARNING: matrix is NOT symmetric within tol=%g.\n",
                (double)symmetryTol));
            if (symmetryStrict) {
                PetscCall(PetscPrintf(comm,
                    "[permon_solve] Aborting because -permon_check_symmetry_strict is enabled.\n"));
                PetscCall(permon_cache_reset(&g_cache));
                return PETSC_ERR_ARG_WRONGSTATE;
            }
        }
    }

    PetscCall(VecSet(b, 0.0));
    PetscCall(VecSet(c, 0.0));
    PetscCall(VecSet(x, 0.0));

    if (b_ptr) {
        PetscCall(VecSetValues(b, nrows, globaldofs,
                       (PetscScalar*)b_ptr, ADD_VALUES));
        PetscCall(VecAssemblyBegin(b));
        PetscCall(VecAssemblyEnd(b));
    }

    if (x_ptr) {
        PetscCall(VecSetValues(x, nrows, globaldofs,
                       (PetscScalar*)x_ptr, INSERT_VALUES));
        PetscCall(VecAssemblyBegin(x));
        PetscCall(VecAssemblyEnd(x));
    }

    if (c_ptr) {
        PetscCall(VecSetValues(c, nrows, globaldofs,
                       (PetscScalar*)c_ptr, INSERT_VALUES));
        PetscCall(VecAssemblyBegin(c));
        PetscCall(VecAssemblyEnd(c));
    }



    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    * Setup QP: argmin 1/2 x'Ax -x'b s.t. c <= x
    *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */    

    solveCallCounter++;

    PetscCall(PetscOptionsGetBool(NULL, NULL, "-permon_debug_initial_guess", &debugInit, NULL));
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-permon_debug_bounds", &debugBounds, NULL));
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-permon_initial_guess_at_bound", &pinInitToBound, NULL));
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-permon_initial_guess_at_bound_first", &pinInitToBoundAtFirst, NULL));

    if (solveCallCounter == 1) {
        PetscCall(PetscPrintf(comm,"permon_solve: Saving A and b"));
        MatViewFromOptions(A,NULL,"-amat_view");
        VecViewFromOptions(b,NULL,"-bvec_view");
    }

    /* Re-attach updated matrix/rhs to reused QP objects. */
    PetscCall(QPSetOperator(qp, A));
    PetscCall(QPSetRhs(qp, b));
    /* Set box constraints.
    * For PERMON's QPCBox we must always provide a valid lower bound vector. */
    if (bound == 0) {
        /* Lower bound provided, fabricate an upper bound at +infinity. */
        lb = c;
        ub = g_cache.ub_fill;
    } else if (bound == 1) {
        /* Upper bound provided, fabricate a lower bound at -infinity. */
        ub = c;
        lb = g_cache.lb_fill;
    } else {
        return -1;
    }

    if (solveCallCounter == 1) {
        PetscCall(PetscPrintf(comm,"permon_solve: Saving c"));
        VecViewFromOptions(c,NULL,"-cvec_view");
    }

    if (debugInit) {
        PetscReal xNorm2 = 0.0, xMin = 0.0, xMax = 0.0;
        PetscReal cNorm2 = 0.0, cMin = 0.0, cMax = 0.0;
        PetscCall(VecNorm(x, NORM_2, &xNorm2));
        PetscCall(VecMin(x, NULL, &xMin));
        PetscCall(VecMax(x, NULL, &xMax));
        PetscCall(VecNorm(c, NORM_2, &cNorm2));
        PetscCall(VecMin(c, NULL, &cMin));
        PetscCall(VecMax(c, NULL, &cMax));
        PetscCall(PetscPrintf(comm,
            "[permon_solve #%" PetscInt_FMT "] initial guess from Elmer: ||x||_2=%.6e min=%.6e max=%.6e; bound=%s ||c||_2=%.6e min=%.6e max=%.6e\n",
            solveCallCounter,
            (double)xNorm2, (double)xMin, (double)xMax,
            (bound == 0) ? "lower" : "upper",
            (double)cNorm2, (double)cMin, (double)cMax));
    }

    if (pinInitToBound || (pinInitToBoundAtFirst && solveCallCounter == 1)) {
        /* Optional experiment: start exactly on active bound instead of Elmer's XVec. */
        PetscCall(PetscPrintf(comm,
            "[permon_solve #%" PetscInt_FMT "] pinning initial guess to bound\n",
            solveCallCounter));
        PetscCall(VecCopy(c, x));
        if (debugInit) {
            PetscReal xNorm2 = 0.0, xMin = 0.0, xMax = 0.0;
            PetscReal cNorm2 = 0.0, cMin = 0.0, cMax = 0.0;
            PetscCall(VecNorm(x, NORM_2, &xNorm2));
            PetscCall(VecMin(x, NULL, &xMin));
            PetscCall(VecMax(x, NULL, &xMax));
            PetscCall(VecNorm(c, NORM_2, &cNorm2));
            PetscCall(VecMin(c, NULL, &cMin));
            PetscCall(VecMax(c, NULL, &cMax));
            PetscCall(PetscPrintf(comm,
                "[permon_solve #%" PetscInt_FMT "] x after pinning: ||x||_2=%.6e min=%.6e max=%.6e; bound=%s ||c||_2=%.6e min=%.6e max=%.6e\n",
                solveCallCounter,
                (double)xNorm2, (double)xMin, (double)xMax,
                (bound == 0) ? "lower" : "upper",
                (double)cNorm2, (double)cMin, (double)cMax));
        }

    }

    /* Set initial guess.
    * THIS VECTOR WILL ALSO HOLD THE SOLUTION OF QP */
    PetscCall(QPSetInitialVector(qp, x));

    /* Optional diagnostic: print min/max of lb, ub and c to catch NaN/Inf or inverted bounds. */
    if (debugBounds) {
        PetscReal lb_min = 0.0, lb_max = 0.0, ub_min = 0.0, ub_max = 0.0, c_min = 0.0, c_max = 0.0;
        PetscCall(VecMin(lb, NULL, &lb_min));
        PetscCall(VecMax(lb, NULL, &lb_max));
        PetscCall(VecMin(ub, NULL, &ub_min));
        PetscCall(VecMax(ub, NULL, &ub_max));
        PetscCall(VecMin(c, NULL, &c_min));
        PetscCall(VecMax(c, NULL, &c_max));
        PetscCall(PetscPrintf(comm,
            "[permon_solve] Bounds check: lb[min,max]=%.12g,%.12g ub[min,max]=%.12g,%.12g c[min,max]=%.12g,%.12g\n",
            (double)lb_min, (double)lb_max, (double)ub_min, (double)ub_max, (double)c_min, (double)c_max));
    }

    PetscCall(QPSetBox(qp, NULL, lb, ub));
     /* Apply runtime options after constraints are attached so QPS type
         selection sees the final QP shape (boxed vs unconstrained). */
     PetscCall(QPSetFromOptions(qp));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    * Setup QPS, i.e. QP Solver
    *   Note the use of PetscObjectComm() to get the same comm as in qp object.
    *   We could specify the comm explicitly, in this case PETSC_COMM_WORLD.
    *   Also, all PERMON objects are PETSc objects as well :)
    *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /* Set QP to solve */
    PetscCall(QPSSetQP(qps, qp));
    PetscCall(QPSSetFromOptions(qps));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    * Solve QP
    *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(QPSSolve(qps));

    /* Check that QPS converged */
    PetscCall(QPIsSolved(qp, &converged));
    if (!converged) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "QPS did not converge!\n"));
    if(converged) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "QPS converged!\n"));

    /* Copy solution back into Fortran array x_ptr (if provided) so Elmer
       sees the computed solution for post-processing and norm checks. */
    if (x_ptr) {
        PetscCall(VecScatterBegin(g_cache.scatter, x, g_cache.out, INSERT_VALUES, SCATTER_FORWARD));
        PetscCall(VecScatterEnd(g_cache.scatter, x, g_cache.out, INSERT_VALUES, SCATTER_FORWARD));

        /* Copy out values into Fortran buffer */
        {
            const PetscScalar *valsbuf = NULL;
            PetscCall(VecGetArrayRead(g_cache.out, &valsbuf));
            PetscScalar *x_out = (PetscScalar*) x_ptr;
            for (i = 0; i < nrows; ++i) x_out[i] = valsbuf[i];
            PetscCall(VecRestoreArrayRead(g_cache.out, &valsbuf));
        }
    }

    return 0;
}
