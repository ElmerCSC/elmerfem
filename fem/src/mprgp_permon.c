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
    // permonrc - default name for solver options file
    return PermonInitialize(NULL, NULL, "permonrc", NULL);
}

int permon_finalize(){
    return PermonFinalize();
}

// TODO check if the freeing of the arrays is correct
int permon_solve(void *rows, void *cols, void *vals, int nrows, int ncols, void *b_ptr, void *c_ptr, void *x_ptr, int bound, void *gdofs_ptr, void *owner_ptr){
    Vec       b = NULL, c = NULL, x = NULL;
    Vec       lb_fill = NULL, ub_fill = NULL;
    Vec       lb = NULL, ub = NULL;
    Mat       A;
    QP        qp;
    QPS       qps;
    PetscBool converged, viewSol = PETSC_FALSE;
    PetscViewer viewer;
    PetscMPIInt size, rank;
    PetscScalar *vals_array = (PetscScalar*)vals;
    PetscScalar *b_array = (PetscScalar*)b_ptr;
    PetscScalar *c_array = (PetscScalar*)c_ptr;
    PetscScalar *x_array = (PetscScalar*)x_ptr;
    PetscInt  nrows_array = nrows + 1;
    PetscInt  nnz = 0;
    PetscInt *rows_c = NULL, *cols_local = NULL, *cols_global = NULL;
    PetscInt *gdofs = NULL;
    int *rows_f = (int*)rows;
    int *cols_f = (int*)cols;
    int *owner_in = (int*)owner_ptr;
    int *gdofs_in = (int*)gdofs_ptr;
    PetscInt local_owned = 0;
    OwnedRow *owned_rows = NULL;
    PetscInt *d_nnz = NULL, *o_nnz = NULL;
    PetscScalar *b_owned = NULL, *c_owned = NULL, *x_owned = NULL;
    PetscInt global_rows = 0;
    PetscInt global_cols = ncols;
    PetscInt global_start = 0;
    PetscInt i;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

        printf("permon_solve: entry: nrows=%d, rows=%p, cols=%p, vals=%p, gdofs=%p, owner=%p\n",
               nrows, (void*)rows, (void*)cols, (void*)vals, (void*)gdofs_in, (void*)owner_in);
        fflush(stdout);

        /* Print only first few entries of the arrays to avoid huge output */
        int maxprint = 10;
        if (rows_f) {
            int rprint = (nrows + 1 < maxprint) ? (nrows + 1) : maxprint;
            printf(" rows_f[0..%d]:", rprint - 1);
            for (i = 0; i < rprint; ++i) printf(" %d", rows_f[i]);
            printf("\n");
        } else {
            printf(" rows_f: NULL\n");
        }

        if (cols_f) {
            int cprint = maxprint; /* we don't yet know nnz here; just peek */
            printf(" cols_f[0..%d]:", cprint - 1);
            for (i = 0; i < cprint; ++i) printf(" %d", cols_f[i]);
            printf("\n");
        } else {
            printf(" cols_f: NULL\n");
        }

        if (gdofs_in) {
            int gprint = (nrows < maxprint) ? nrows : maxprint;
            printf(" gdofs_in[0..%d]:", gprint - 1);
            for (i = 0; i < gprint; ++i) printf(" %d", gdofs_in[i]);
            printf("\n");
        } else {
            printf(" gdofs_in: NULL\n");
        }

        if (owner_in) {
            int oprint = (nrows < maxprint) ? nrows : maxprint;
            printf(" owner_in[0..%d]:", oprint - 1);
            for (i = 0; i < oprint; ++i) printf(" %d", owner_in[i]);
            printf("\n");
        } else {
            printf(" owner_in: NULL\n");
        }
        fflush(stdout);

    if (nrows > 0) {
        if (!rows || !cols || !gdofs_in || !owner_in) {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_NULL, "Missing matrix metadata");
        }
        nnz = rows_f[nrows] - 1;
        PetscCall(PetscMalloc1(nrows_array, &rows_c));
        PetscCall(PetscMalloc1(nnz, &cols_local));
        for (i = 0; i < nrows_array; ++i) {
            rows_c[i] = (PetscInt)(rows_f[i] - 1);
        }
        for (i = 0; i < nnz; ++i) {
            cols_local[i] = (PetscInt)(cols_f[i] - 1);
        }
        PetscCall(PetscMalloc1(nrows, &gdofs));
        for (i = 0; i < nrows; ++i) {
            gdofs[i] = (PetscInt)gdofs_in[i];
            if (owner_in[i] == 1) {
                local_owned++;
            }
        }
    } else {
        nnz = 0;
        local_owned = 0;
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

    if (local_owned > 0) {
        PetscCall(PetscMalloc1(local_owned, &owned_rows));
        PetscInt pos = 0;
        for (i = 0; i < nrows; ++i) {
            if (owner_in[i] == 1) {
                owned_rows[pos].local = i;
                owned_rows[pos].global = gdofs[i];
                pos++;
            }
        }
        qsort(owned_rows, (size_t)local_owned, sizeof(OwnedRow), compare_owned_rows);

        /* Determine PETSc ownership range for rows (diagonal column block) */
        PetscInt rstart = 0, rend = 0;
        Mat A_tmp;
        PetscCall(MatCreate(PETSC_COMM_WORLD, &A_tmp));
        PetscCall(MatSetSizes(A_tmp, local_owned, PETSC_DECIDE, global_rows, global_cols));
        PetscCall(MatSetType(A_tmp, MATMPIAIJ));
        PetscCall(MatSetUp(A_tmp));
        PetscCall(MatGetOwnershipRange(A_tmp, &rstart, &rend));
        PetscCall(MatDestroy(&A_tmp));

        PetscCall(PetscMalloc1(local_owned, &d_nnz));
        PetscCall(PetscMalloc1(local_owned, &o_nnz));
        for (i = 0; i < local_owned; ++i) {
            PetscInt row = owned_rows[i].local;
            PetscInt start = rows_c[row];
            PetscInt end   = rows_c[row+1];
            PetscInt d = 0, o = 0;
            for (PetscInt jj = start; jj < end; ++jj) {
                PetscInt gcol = cols_global[jj];
                if (gcol >= rstart && gcol < rend) d++;
                else o++;
            }
            d_nnz[i] = d;
            o_nnz[i] = o;
        }
    }

    PetscInt tmp = local_owned;
    PetscCallMPI(MPI_Allreduce(&tmp, &global_rows, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD));

    PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
    PetscCall(MatSetSizes(A, local_owned, PETSC_DECIDE, global_rows, global_cols));
    PetscCall(MatSetType(A, MATMPIAIJ));
     /* Allow PETSc to allocate new nonzeros if our preallocation is slightly off.
         This avoids hard aborts during debugging; for best performance compute
         accurate d_nnz/o_nnz using PETSc ownership ranges. */
     PetscCall(MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
    if (local_owned > 0) {
        PetscCall(MatMPIAIJSetPreallocation(A, 0, d_nnz, 0, o_nnz));
    } else {
        PetscCall(MatMPIAIJSetPreallocation(A, 0, NULL, 0, NULL));
    }
    PetscCall(MatSetUp(A));

    for (i = 0; i < local_owned; ++i) {
        PetscInt row = owned_rows[i].local;
        PetscInt grow = owned_rows[i].global;
        PetscInt start = rows_c[row];
        PetscInt end   = rows_c[row+1];
        PetscInt ncols_row = end - start;
        if (ncols_row > 0) {
            /* Sanity check: expected preallocation counts for this local row */
            PetscInt expected = 0;
            if (d_nnz && o_nnz) expected = d_nnz[i] + o_nnz[i];
            if (expected != ncols_row) {
                PetscPrintf(PETSC_COMM_WORLD, "permon_solve: prealloc mismatch for local row %d (grow=%d): ncols_row=%d, d_nnz=%d, o_nnz=%d, start=%d, end=%d\n",
                            row, (int)grow, (int)ncols_row, (d_nnz? (int)d_nnz[i] : -1), (o_nnz? (int)o_nnz[i] : -1), (int)start, (int)end);
                PetscPrintf(PETSC_COMM_WORLD, "  cols_global[start..end-1]:");
                for (PetscInt jj = start; jj < end; ++jj) PetscPrintf(PETSC_COMM_WORLD, " %d", (int)cols_global[jj]);
                PetscPrintf(PETSC_COMM_WORLD, "\n");
                PetscPrintf(PETSC_COMM_WORLD, "  rows_c sample [row-1..row+1]:");
                for (int jj = row-1; jj <= row+1; ++jj) if (jj >= 0 && jj < nrows+1) PetscPrintf(PETSC_COMM_WORLD, " %d", (int)rows_c[jj]);
                PetscPrintf(PETSC_COMM_WORLD, "\n");
                fflush(stdout);
            }
            PetscCall(MatSetValues(A, 1, &grow, ncols_row, &cols_global[start], &vals_array[start], INSERT_VALUES));
        }
    }

    PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

    if (b_ptr != NULL) {
        if (local_owned > 0) {
            PetscCall(PetscMalloc1(local_owned, &b_owned));
            for (i = 0; i < local_owned; ++i) {
                b_owned[i] = b_array[owned_rows[i].local];
            }
        }
        PetscScalar *arr = (local_owned > 0) ? b_owned : NULL;
        PetscCall(VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, local_owned, global_rows, arr, &b));
    } else {
        PetscCall(MatCreateVecs(A, &b, NULL));
    }

    if (c_ptr != NULL) {
        if (local_owned > 0) {
            PetscCall(PetscMalloc1(local_owned, &c_owned));
            for (i = 0; i < local_owned; ++i) {
                c_owned[i] = c_array[owned_rows[i].local];
            }
        }
        PetscScalar *arr = (local_owned > 0) ? c_owned : NULL;
        PetscCall(VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, local_owned, global_rows, arr, &c));
    } else {
        PetscCall(MatCreateVecs(A, NULL, &c));
    }

    if (x_ptr != NULL) {
        if (local_owned > 0) {
            PetscCall(PetscMalloc1(local_owned, &x_owned));
            for (i = 0; i < local_owned; ++i) {
                x_owned[i] = x_array[owned_rows[i].local];
            }
        }
        PetscScalar *arr = (local_owned > 0) ? x_owned : NULL;
        PetscCall(VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, local_owned, global_rows, arr, &x));
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
