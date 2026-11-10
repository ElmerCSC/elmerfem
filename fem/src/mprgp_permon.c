#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <permonqps.h>
#include <petscmat.h>
#include <petscsys.h>

// Cache is used as interface can be called mutliple times within Solver
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

// Destroy all cached PETSc/PERMON objects and mark cache as empty.
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

// Check if the cached objects remain the same and can be reused.
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
    // Rebuild cache from scratch when communicator/layout changed.
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
    PetscCall(MatSetOption(cache->A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
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


// Refresh cached CRS indices from Fortran row/col arrays.
static PetscErrorCode permon_cache_pattern(PermonSolveCache *cache, PetscInt nrows, const int *rows_f, const int *cols_f)
{
    PetscInt i, nnz_total_new;

    PetscFunctionBegin;
    nnz_total_new = (PetscInt)rows_f[nrows] - 1;
    if (nnz_total_new < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Invalid CRS row pointer end value");

    if (!cache->pattern_cached || nnz_total_new != cache->nnz_total) {
        PetscCall(PetscFree(cache->cols0));
        PetscCall(PetscMalloc1(nnz_total_new, &cache->cols0));
        cache->nnz_total = nnz_total_new;
    }

    cache->maxnnz = 0;
    for (i = 0; i < nrows; ++i) {
        PetscInt rownnz = (PetscInt)rows_f[i + 1] - (PetscInt)rows_f[i];
        cache->maxnnz = PetscMax(cache->maxnnz, rownnz);
    }

    for (i = 0; i < cache->nnz_total; ++i) {
        /* Convert Fortran 1-based local columns to C 0-based for this solve. */
        cache->cols0[i] = (PetscInt)cols_f[i] - 1;
    }

    cache->pattern_cached = PETSC_TRUE;
    PetscFunctionReturn(0);
}

int permon_init(){
    return PermonInitialize(NULL, NULL, NULL, NULL);
}

int permon_set_options_file(const char *options_file, int fcomm){
    MPI_Comm comm = MPI_Comm_f2c(fcomm);

    if (options_file == NULL || options_file[0] == '\0') {
        return 0;
    }

    PetscErrorCode ierr;
    ierr = PetscOptionsInsertFile(comm, NULL, options_file, PETSC_TRUE);
    PetscCall(PetscLogDefaultBegin());

    return ierr;
}

int permon_finalize(){
    PetscCall(permon_cache_reset(&g_cache));
    return PermonFinalize();
}

int permon_solve(void *rows_local, void *cols_local, void *vals_local, int nrows, void *b_ptr, void *c_ptr, void *x_ptr, int bound, int *globaldofs, int *owner, int fcomm) {
    Vec       b, c, x;
    Vec       lb = NULL, ub = NULL;
    Mat       A, AT = NULL;
    QP        qp;
    QPS       qps;
    PetscInt  i, rstart, rend, nnz;
    PetscBool converged, viewSol = PETSC_FALSE;
    PetscBool pinInitToBound = PETSC_FALSE, pinInitToBoundAtFirst = PETSC_FALSE;
    PetscViewer viewer;
    static PetscInt solveCallCounter = 0;
    
    MPI_Comm comm=MPI_Comm_f2c(fcomm);
    
    int *rows_f = (int*)rows_local;
    int *cols_f = (int*)cols_local;

    double *vals = (double*)vals_local;

    // -----------------------------
    // Compute local ownership range
    // -----------------------------
    PetscInt nlocal = 0;


    // Find number of local rows owned by this rank
    for (i = 0; i < nrows; i++) {
        if (owner[i]) {
            nlocal++;
        }
    }
    
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

    // Set Hessian values
    for (i = 0; i < nrows; i++) {
        nnz  = rows_f[i+1] - rows_f[i];
        PetscInt iloc = i; /* local row index */

        PetscCall(MatSetValuesLocal(A, 1, &iloc, nnz, &g_cache.cols0[rows_f[i] - 1],
                &vals[rows_f[i] - 1], ADD_VALUES));
    }


    PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));



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
    * Setup QP: argmin 1/2 x'Ax -x'b s.t. lb <= x <= ub 
    *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */    

    solveCallCounter++;

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

    if (pinInitToBound || (pinInitToBoundAtFirst && solveCallCounter == 1)) {
        /* Optional experiment: start exactly on active bound instead of Elmer's XVec. */
        PetscCall(PetscPrintf(comm,
            "[permon_solve #%" PetscInt_FMT "] pinning initial guess to bound\n",
            solveCallCounter));
        PetscCall(VecCopy(c, x));
    }

    /* Set initial guess.
    * THIS VECTOR WILL ALSO HOLD THE SOLUTION OF QP */
    PetscCall(QPSetInitialVector(qp, x));

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

    /* Copy solution back into Fortran array x_ptr */
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
