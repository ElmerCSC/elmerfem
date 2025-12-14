/* Simple test helper to print CRS row pointers from C
 * This function is called from Fortran for debugging.
 */
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <permonqps.h>
#include <mpi.h>

typedef struct {
    PetscInt local;
    PetscInt global;
} OwnedRow;

static int compare_owned_rows(const void *a, const void *b)
{
    const OwnedRow *ra = (const OwnedRow*)a;
    const OwnedRow *rb = (const OwnedRow*)b;
    if (ra->global < rb->global) return -1;
    if (ra->global > rb->global) return 1;
    return 0;
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
    PetscErrorCode ierr;
    // permonrc - default name for solver options file
    ierr = PermonInitialize(NULL, NULL, "permonrc", NULL);
    printf("permon_init: PermonInitialize returned error code %d\n", ierr);
    if (ierr) {
        printf("permon_init: PermonInitialize failed with error code %d\n", ierr);
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, ierr);
    }
    return ierr;

}

int permon_finalize(){
    return PermonFinalize();
}

// TODO check if the freeing of the arrays is correct
int permon_solve(void *rows_local, void *cols_local, void *vals_local, int nrows, int ncols, void *b_ptr, void *c_ptr, void *x_ptr, int bound){
    Vec       b, c, x;
    Vec       lb_fill = NULL, ub_fill = NULL;
    Vec       lb = NULL, ub = NULL;
    Mat       A;
    QP        qp;
    QPS       qps;
    PetscBool converged, viewSol = PETSC_FALSE;
    PetscViewer viewer;
    
    /* Convert Fortran 1-based indices to C 0-based indices for PETSc */
    int *rows_f = (int*)rows_local;  /* Fortran 1-based row pointers */
    int *cols_f = (int*)cols_local;  /* Fortran 1-based column indices */
    int nrows_array = nrows + 1;  /* Rows array has nrows+1 elements */


    if (rows_f == NULL) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "rows array is NULL");
    }

    if (nrows < 1) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "nrows < 1");
    }

    int nnz = rows_f[nrows] - 1;


    mprgp_print_vector(c_ptr, 15, "c");
    mprgp_print_vector(b_ptr, 15, "b");
    
    /* Allocate temporary arrays for 0-based indices */
    PetscInt *rows_c, *cols_c;
    PetscMalloc1(nrows_array, &rows_c);
    PetscMalloc1(nnz, &cols_c);
    
    /* Convert row pointers: subtract 1 from each element */
    for (i = 0; i < nrows_array; i++) {
        rows_c[i] = (PetscInt)(rows_f[i] - 1);
    }

    if (nnz > 0) {
        PetscCall(PetscMalloc1(nnz, &cols_global));
        for (i = 0; i < nnz; ++i) {
            PetscInt loc = cols_local[i];
            if (loc >= 0 && loc < nrows) {
                cols_global[i] = gdofs[loc];
            } else {
                cols_global[i] = loc;
            }
        }
    }

    /* Validate computed global column indices to catch out-of-range errors early */
    if (cols_global) {
        for (i = 0; i < nnz; ++i) {
            if (cols_global[i] < 0 || cols_global[i] >= global_cols) {
                PetscPrintf(PETSC_COMM_WORLD, "permon_solve: INVALID cols_global[%d] = %d (nnz=%d, global_cols=%d)\n",
                            i, (int)cols_global[i], (int)nnz, (int)global_cols);
                PetscPrintf(PETSC_COMM_WORLD, "  sample rows_c[0..5]:");
                for (int j = 0; j < 6 && j < nrows+1; ++j) PetscPrintf(PETSC_COMM_WORLD, " %d", (int)rows_c[j]);
                PetscPrintf(PETSC_COMM_WORLD, "\n");
                PetscPrintf(PETSC_COMM_WORLD, "  sample gdofs[0..5]:");
                for (int j = 0; j < 6 && j < nrows; ++j) PetscPrintf(PETSC_COMM_WORLD, " %d", (int)gdofs[j]);
                PetscPrintf(PETSC_COMM_WORLD, "\n");
                PetscPrintf(PETSC_COMM_WORLD, "  owner_in[0..5]:");
                for (int j = 0; j < 6 && j < nrows; ++j) PetscPrintf(PETSC_COMM_WORLD, " %d", owner_in[j]);
                PetscPrintf(PETSC_COMM_WORLD, "\n");
                PetscPrintf(PETSC_COMM_WORLD, "  Aborting with PETSC error to avoid reallocation.\n");
                PetscMPIInt _ierr = 0;
                PetscMPIInt _rank;
                MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);
                MPI_Abort(PETSC_COMM_WORLD, 1);
            }
        }
    }

     /* Create matrix directly from arrays (MatCreateSeqAIJWithArrays creates a new matrix) */
     /* We're using PETSC_COMM_SELF here so the number of local rows is the
         full row count `nrows`. Define `local_nrows` accordingly. */
     PetscInt local_nrows = nrows;
     /* For a sequential communicator the local columns are also the full
         column count; use PETSC_DECIDE for MPI-aware behavior or explicitly
         pass `ncols` if desired. */
     PetscCall(MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, local_nrows, PETSC_DECIDE, nrows, ncols, rows_c, cols_c, vals_local, &A));
    // PetscCall(MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, nrows, ncols, rows_c, cols_c, vals_local, &A));
    
    /* Create vectors from Fortran arrays: b from RHS, c from LowerLimit */
    /* VecCreateSeqWithArray wraps existing array data (doesn't copy) */
    if (b_ptr != NULL) {
        PetscScalar *b_array = (PetscScalar*)b_ptr;
        // PetscCall(VecCreateSeqWithArray(PETSC_COMM_SELF, 1, nrows, b_array, &b)); // does not copy the datam verify this
        PetscCall(VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, nrows, PETSC_DETERMINE, b_array, &b));
    } else {
        PetscCall(MatCreateVecs(A, &b, NULL));
    }

    if (c_ptr != NULL) {
        PetscScalar *c_array = (PetscScalar*)c_ptr;
        // PetscCall(VecCreateSeqWithArray(PETSC_COMM_SELF, 1, nrows, c_array, &c));
        PetscCall(VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, nrows, PETSC_DETERMINE, c_array, &c));
    } else {
        PetscCall(MatCreateVecs(A, NULL, &c));
    }

    if (x_ptr != NULL) {
        PetscScalar *x_array = (PetscScalar*)x_ptr;
        // PetscCall(VecCreateSeqWithArray(PETSC_COMM_SELF, 1, nrows, x_array, &x));
        PetscCall(VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, nrows, PETSC_DETERMINE, x_array, &x));
    } else {
        PetscCall(MatCreateVecs(A, NULL, &x));
    }

    PetscCall(QPCreate(PETSC_COMM_WORLD, &qp));
    PetscCall(QPSetOperator(qp, A));
    PetscCall(QPSetRhs(qp, b));
    PetscCall(QPSetInitialVector(qp, x));
    if (bound == 0) {
        lb = c;
        PetscCall(VecDuplicate(c, &ub_fill));
        PetscCall(VecSet(ub_fill, PETSC_INFINITY));
        ub = ub_fill;
    } else if (bound == 1) {
        ub = c;
        PetscCall(VecDuplicate(c, &lb_fill));
        PetscCall(VecSet(lb_fill, PETSC_NINFINITY));
        lb = lb_fill;
    } else {
        return -1;
    }
    PetscCall(QPSetBox(qp, NULL, lb, ub));
    PetscCall(QPSetFromOptions(qp));

    PetscCall(QPSCreate(PetscObjectComm((PetscObject)qp), &qps));
    PetscCall(QPSSetQP(qps, qp));
    PetscCall(QPSSetFromOptions(qps));

    PetscCall(QPSSolve(qps));

    if (x_ptr != NULL && local_owned > 0 && x_owned != NULL) {
        for (i = 0; i < local_owned; ++i) {
            x_array[owned_rows[i].local] = x_owned[i];
        }
    }

    PetscCall(QPIsSolved(qp, &converged));
    if (!converged) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "QPS did not converge!\n"));
    if(converged) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "QPS converged!\n"));

    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&c));
    PetscCall(VecDestroy(&lb_fill));
    PetscCall(VecDestroy(&ub_fill));
    PetscCall(VecDestroy(&b));
    PetscCall(MatDestroy(&A));

    PetscCall(QPSDestroy(&qps));
    PetscCall(QPDestroy(&qp));

    PetscCall(PetscFree(rows_c));
    PetscCall(PetscFree(cols_local));
    PetscCall(PetscFree(cols_global));
    PetscCall(PetscFree(gdofs));
    PetscCall(PetscFree(owned_rows));
    PetscCall(PetscFree(d_nnz));
    PetscCall(PetscFree(o_nnz));
    PetscCall(PetscFree(b_owned));
    PetscCall(PetscFree(c_owned));
    PetscCall(PetscFree(x_owned));

    return 0;
}
