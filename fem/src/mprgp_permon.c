/* Simple test helper to print CRS row pointers from C
 * This function is called from Fortran for debugging.
 */
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <permonqps.h>
#include <petscmat.h>
#include <petscsys.h>


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
    return PermonFinalize();
}

// TODO check if the freeing of the arrays is correct
int permon_solve(void *rows_local, void *cols_local, void *vals_local, int nrows, void *b_ptr, void *c_ptr, void *x_ptr, int bound, int *globaldofs, int *owner, int fcomm) {
    Vec       b, c, x;
    Vec       lb_fill = NULL, ub_fill = NULL;
    Vec       lb = NULL, ub = NULL;
    Mat       A;
    QP        qp;
    QPS       qps;
    PetscInt  i, rstart, rend, nnz;
    PetscBool converged, viewSol = PETSC_FALSE;
    PetscBool debugInit = PETSC_FALSE, pinInitToBound = PETSC_FALSE, pinInitToBoundAfterFirst = PETSC_FALSE;
    PetscBool checkSymmetry = PETSC_TRUE, symmetryStrict = PETSC_FALSE, isSymmetric = PETSC_FALSE;
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
    
    // TODO maybe can be done simpler (look at previous code versions)
    ISLocalToGlobalMapping ltog;
    {
        PetscInt *gmap = NULL;
        PetscCall(PetscMalloc1(nrows, &gmap));
        for (i = 0; i < nrows; ++i) gmap[i] = (PetscInt) globaldofs[i];

        PetscCall(
            ISLocalToGlobalMappingCreate(
                comm,
                1,                 /* block size (1 = scalar FEM) */
                nrows,             /* number of local dofs on this rank */
                gmap,              /* local -> global map (PetscInt) */
                PETSC_COPY_VALUES,
                &ltog
            )
        );

        PetscCall(PetscFree(gmap));
    }

    // -----------------------------
    // 2. Create the MPI matrix
    // -----------------------------
    // TODO probably can be optimized by precomputing row nnz
    PetscCall(MatCreate(comm, &A));
    PetscCall(MatSetType(A, MATMPIAIJ));
    PetscCall(MatSetSizes(A, nlocal, nlocal, PETSC_DECIDE, PETSC_DECIDE));
    MatSetFromOptions(A);
    PetscCall(MatSetUp(A));
    /* Attach local->global mapping so we can insert with local indices */
    PetscCall(MatSetLocalToGlobalMapping(A, ltog, ltog));
    PetscCall(ISLocalToGlobalMappingDestroy(&ltog));

    PetscInt maxnnz = 0;
    for (i = 0; i < nrows; i++) {
        PetscInt nnz = rows_f[i+1] - rows_f[i];
        maxnnz = PetscMax(maxnnz, nnz);
    }
    PetscInt *rcols;
    PetscCall(PetscMalloc1(maxnnz, &rcols));

    for (i = 0; i < nrows; i++) {
        nnz  = rows_f[i+1] - rows_f[i];
        PetscInt iloc = i; /* local row index */

        for (PetscInt k = 0, j = rows_f[i]; j < rows_f[i+1]; j++, k++) {
            /* Local column index (Fortran->C): cols_f[j-1]-1 */
            rcols[k] = cols_f[j-1] - 1;
        }

        PetscCall(MatSetValuesLocal(A, 1, &iloc, nnz, rcols,
                &vals[rows_f[i] - 1], ADD_VALUES));
    }
    PetscCall(PetscFree(rcols));


    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

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
                PetscCall(MatDestroy(&A));
                return PETSC_ERR_ARG_WRONGSTATE;
            }
        }
    }

    // -----------------------------
    // 4. Create vectors (b, c, x)
    // -----------------------------
    PetscCall(MatCreateVecs(A, &x, &b));
    PetscCall(VecDuplicate(x, &c));

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
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-permon_initial_guess_at_bound", &pinInitToBound, NULL));
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-permon_initial_guess_at_bound_first", &pinInitToBoundAfFirst, NULL));

    if (solveCallCounter == 1) {
        PetscCall(PetscPrintf(comm,"permon_solve: Saving A and b"));
        MatViewFromOptions(A,NULL,"-amat_view");
        VecViewFromOptions(b,NULL,"-bvec_view");
    }

    PetscCall(QPCreate(comm, &qp));
    /* Set matrix representing QP operator */
    PetscCall(QPSetOperator(qp, A));
    /* Set right hand side */
    PetscCall(QPSetRhs(qp, b));
    /* Set box constraints.
    * For PERMON's QPCBox we must always provide a valid lower bound vector. */
    if (bound == 0) {
        /* Lower bound provided, fabricate an upper bound at +infinity. */
        lb = c;
        PetscCall(VecDuplicate(c, &ub_fill));
        PetscCall(VecSet(ub_fill, PETSC_INFINITY));
        ub = ub_fill;
    } else if (bound == 1) {
        /* Upper bound provided, fabricate a lower bound at -infinity. */
        ub = c;
        PetscCall(VecDuplicate(c, &lb_fill));
        PetscCall(VecSet(lb_fill, PETSC_NINFINITY));
        lb = lb_fill;
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

    if (pinInitToBound || (pinInitToBoundAfFirst && solveCallCounter == 1)) {
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

    PetscCall(QPSetBox(qp, NULL, lb, ub));
    /* Set runtime options, e.g
    *   -qp_chain_view_kkt */
    PetscCall(QPSetFromOptions(qp));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    * Setup QPS, i.e. QP Solver
    *   Note the use of PetscObjectComm() to get the same comm as in qp object.
    *   We could specify the comm explicitly, in this case PETSC_COMM_WORLD.
    *   Also, all PERMON objects are PETSc objects as well :)
    *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(QPSCreate(PetscObjectComm((PetscObject)qp), &qps));
    /* Set QP to solve */
    PetscCall(QPSSetQP(qps, qp));
    /* Set runtime options for solver, e.g,
    *   -qps_type <type> -qps_rtol <relative tolerance> -qps_view_convergence */
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
        /* Efficiently gather the requested global entries into a local sequential
           vector using VecScatter instead of calling VecGetValues for each index. */
        PetscInt *gids = NULL;
        IS is_from = NULL, is_to = NULL;
        Vec out = NULL;
        VecScatter scatter = NULL;

        PetscCall(PetscMalloc1(nrows, &gids));
        for (i = 0; i < nrows; ++i) gids[i] = globaldofs[i];

        PetscCall(ISCreateGeneral(comm, nrows, gids, PETSC_COPY_VALUES, &is_from));
        PetscCall(ISCreateStride(PETSC_COMM_SELF, nrows, 0, 1, &is_to));
        PetscCall(VecCreateSeq(PETSC_COMM_SELF, nrows, &out));

        PetscCall(VecScatterCreate(x, is_from, out, is_to, &scatter));
        PetscCall(VecScatterBegin(scatter, x, out, INSERT_VALUES, SCATTER_FORWARD));
        PetscCall(VecScatterEnd(scatter, x, out, INSERT_VALUES, SCATTER_FORWARD));

        /* Copy out values into Fortran buffer */
        {
            const PetscScalar *valsbuf = NULL;
            PetscCall(VecGetArrayRead(out, &valsbuf));
            PetscScalar *x_out = (PetscScalar*) x_ptr;
            for (i = 0; i < nrows; ++i) x_out[i] = valsbuf[i];
            PetscCall(VecRestoreArrayRead(out, &valsbuf));
        }

        PetscCall(VecScatterDestroy(&scatter));
        PetscCall(ISDestroy(&is_from));
        PetscCall(ISDestroy(&is_to));
        PetscCall(VecDestroy(&out));
        PetscCall(PetscFree(gids));
    }

    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&c));
    PetscCall(VecDestroy(&lb_fill));
    PetscCall(VecDestroy(&ub_fill));
    PetscCall(VecDestroy(&b));
    PetscCall(MatDestroy(&A));

    PetscCall(QPSDestroy(&qps));
    PetscCall(QPDestroy(&qp));

    return 0;
}
