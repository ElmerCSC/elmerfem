/* Simple test helper to print CRS row pointers from C
 * This function is called from Fortran for debugging.
 */
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <permonqps.h>


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
    return PermonInitialize(NULL, NULL, (char *)0, NULL);
}

int permon_finalize(){
    return PermonFinalize();
}

// TODO check if the freeing of the arrays is correct
int permon_solve(void *rows, void *cols, void *vals, int nrows, int ncols, void *b_ptr, void *c_ptr, void *x_ptr, int bound){
    Vec       b, c, x;
    Vec       lb_fill = NULL, ub_fill = NULL;
    Vec       lb = NULL, ub = NULL;
    Mat       A;
    QP        qp;
    QPS       qps;
    PetscInt  i;
    PetscBool converged, viewSol = PETSC_FALSE;
    PetscViewer viewer;
    
    /* Convert Fortran 1-based indices to C 0-based indices for PETSc */
    int *rows_f = (int*)rows;  /* Fortran 1-based row pointers */
    int *cols_f = (int*)cols;  /* Fortran 1-based column indices */
    int nrows_array = nrows + 1;  /* Rows array has nrows+1 elements */
    int nnz = rows_f[nrows] - 1;  /* Number of nonzeros (last element - 1, in 1-based) */

    mprgp_print_vector(c_ptr, 15, "c");
    
    /* Allocate temporary arrays for 0-based indices */
    PetscInt *rows_c, *cols_c;
    PetscMalloc1(nrows_array, &rows_c);
    PetscMalloc1(nnz, &cols_c);
    
    /* Convert row pointers: subtract 1 from each element */
    for (i = 0; i < nrows_array; i++) {
        rows_c[i] = (PetscInt)(rows_f[i] - 1);
    }
    
    /* Convert column indices: subtract 1 from each element */
    for (i = 0; i < nnz; i++) {
        cols_c[i] = (PetscInt)(cols_f[i] - 1);
    }


    /* Create matrix directly from arrays (MatCreateSeqAIJWithArrays creates a new matrix) */
    PetscCall(MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, nrows, ncols, rows_c, cols_c, vals, &A));
    
    /* Create vectors from Fortran arrays: b from RHS, c from LowerLimit */
    /* VecCreateSeqWithArray wraps existing array data (doesn't copy) */
    if (b_ptr != NULL) {
        PetscScalar *b_array = (PetscScalar*)b_ptr;
        PetscCall(VecCreateSeqWithArray(PETSC_COMM_SELF, 1, nrows, b_array, &b)); // does not copy the datam verify this
    } else {
        /* If b_ptr is NULL, create empty vector */
        PetscCall(MatCreateVecs(A, &b, NULL));
    }
    
    if (c_ptr != NULL) {
        PetscScalar *c_array = (PetscScalar*)c_ptr;
        PetscCall(VecCreateSeqWithArray(PETSC_COMM_SELF, 1, nrows, c_array, &c));
    } else {
        /* If c_ptr is NULL, create empty vector */
        PetscCall(MatCreateVecs(A, NULL, &c));
    }

    if (x_ptr != NULL) {
        PetscScalar *x_array = (PetscScalar*)x_ptr;
        PetscCall(VecCreateSeqWithArray(PETSC_COMM_SELF, 1, nrows, x_array, &x));
    } else {
        /* If x_ptr is NULL, create empty vector */
        PetscCall(MatCreateVecs(A, NULL, &x));
    }

    // MatView(A, PETSC_VIEWER_STDOUT_SELF);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    * Setup QP: argmin 1/2 x'Ax -x'b s.t. c <= x
    *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(QPCreate(PETSC_COMM_WORLD, &qp));
    /* Set matrix representing QP operator */
    PetscCall(QPSetOperator(qp, A));
    /* Set right hand side */
    PetscCall(QPSetRhs(qp, b));
    /* Set initial guess.
    * THIS VECTOR WILL ALSO HOLD THE SOLUTION OF QP */
    PetscCall(QPSetInitialVector(qp, x));
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

    /* Destroy matrix and vectors - MatCreateSeqAIJWithArrays takes ownership of the arrays
     * and will free them automatically when the matrix is destroyed, so we don't free them here */
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&c));
    PetscCall(VecDestroy(&lb_fill));
    PetscCall(VecDestroy(&ub_fill));
    PetscCall(VecDestroy(&b));
    PetscCall(MatDestroy(&A));
    /* Note: rows_c and cols_c are automatically freed by MatDestroy above */

    PetscCall(QPSDestroy(&qps));
    PetscCall(QPDestroy(&qp));

    return 0;
}
