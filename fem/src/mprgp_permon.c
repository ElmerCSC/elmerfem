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
    // permonrc - default name for solver options file
    return PermonInitialize(NULL, NULL, "permonrc", NULL);
}

int permon_finalize(){
    return PermonFinalize();
}

// TODO check if the freeing of the arrays is correct
int permon_solve(void *rows_local, void *cols_local, void *vals_local, int nrows, void *b_ptr, void *c_ptr, void *x_ptr, int bound, int *globaldofs, int *owner, int *fcomm) {
    Vec       b, c, x;
    Vec       lb_fill = NULL, ub_fill = NULL;
    Vec       lb = NULL, ub = NULL;
    Mat       A;
    QP        qp;
    QPS       qps;
    PetscInt  i, rstart, rend, nnz;
    PetscBool converged, viewSol = PETSC_FALSE;
    PetscViewer viewer;
    
    MPI_Comm comm=MPI_Comm_f2c(*fcomm);

    int *rows_f = (int*)rows_local;
    int *cols_f = (int*)cols_local;

    double *vals = (double*)vals_local;

    // -----------------------------
    // 1. Compute local ownership range
    // -----------------------------
    PetscInt ilower = PETSC_MAX_INT, iupper = -1;
    PetscInt nlocal = 0;

    /* Dump rows, cols and vals arrays to per-rank files for debugging */
    {
        int rank = 0;
        MPI_Comm_rank(comm, &rank);
        char fname[256];

        /* rows dump (existing) */
        snprintf(fname, sizeof(fname), "rows_dump_rank%d.txt", rank);
        FILE *fp = fopen(fname, "w");
        if (fp) {
            fprintf(fp, "rows_local pointer=%p\n", (void*)rows_local);
            fprintf(fp, "rows_f pointer=%p\n", (void*)rows_f);
            fprintf(fp, "rows_f (0..nloc):\n");
            for (i = 0; i <= nrows; ++i) {
                fprintf(fp, "%d\n", rows_f[i]);
            }
            fclose(fp);
        } else {
            fprintf(stderr, "Failed to open %s for writing\n", fname);
            fflush(stderr);
        }

        /* dump globaldofs */
        snprintf(fname, sizeof(fname), "globaldofs_dump_rank%d.txt", rank);
        fp = fopen(fname, "w");
        if (fp) {
            fprintf(fp, "globaldofs pointer=%p\n", (void*)globaldofs);
            fprintf(fp, "globaldofs (0..nrows-1):\n");
            for (i = 0; i < nrows; ++i) {
                fprintf(fp, "%d\n", globaldofs[i]);
            }
            fclose(fp);
        } else {
            fprintf(stderr, "Failed to open %s for writing\n", fname);
            fflush(stderr);
        }

        /* dump owner */
        snprintf(fname, sizeof(fname), "owner_dump_rank%d.txt", rank);
        fp = fopen(fname, "w");
        if (fp) {
            fprintf(fp, "owner pointer=%p\n", (void*)owner);
            fprintf(fp, "owner (0..nrows-1):\n");
            for (i = 0; i < nrows; ++i) {
                fprintf(fp, "%d\n", owner[i]);
            }
            fclose(fp);
        } else {
            fprintf(stderr, "Failed to open %s for writing\n", fname);
            fflush(stderr);
        }

        /* compute total number of nonzeros robustly */
        PetscInt total_nnz = 0;
        for (i = 0; i < nrows; ++i) {
            total_nnz += rows_f[i+1] - rows_f[i];
        }

        /* cols dump */
        snprintf(fname, sizeof(fname), "cols_dump_rank%d.txt", rank);
        fp = fopen(fname, "w");
        if (fp) {
            fprintf(fp, "cols_local pointer=%p\n", (void*)cols_local);
            fprintf(fp, "cols_f pointer=%p\n", (void*)cols_f);
            fprintf(fp, "total_nnz=%" PRIdPTR "\n", (intptr_t)total_nnz);
            for (PetscInt k = 0; k < total_nnz; ++k) {
                fprintf(fp, "%d\n", cols_f[k]);
            }
            fclose(fp);
        } else {
            fprintf(stderr, "Failed to open %s for writing\n", fname);
            fflush(stderr);
        }

        /* vals dump */
        snprintf(fname, sizeof(fname), "vals_dump_rank%d.txt", rank);
        fp = fopen(fname, "w");
        if (fp) {
            fprintf(fp, "vals_local pointer=%p\n", (void*)vals_local);
            fprintf(fp, "vals pointer=%p\n", (void*)vals);
            fprintf(fp, "total_nnz=%" PRIdPTR "\n", (intptr_t)total_nnz);
            for (PetscInt k = 0; k < total_nnz; ++k) {
                fprintf(fp, "%.18e\n", vals[k]);
            }
            fclose(fp);
        } else {
            fprintf(stderr, "Failed to open %s for writing\n", fname);
            fflush(stderr);
        }
    }
    for (i = 0; i < nrows; i++) {
        if (owner[i]) {
            if (globaldofs[i] < ilower) ilower = globaldofs[i];
            if (globaldofs[i] > iupper) iupper = globaldofs[i];
            nlocal++;
        }
    }
    if (iupper == -1) { ilower = 0; iupper = -1; }
    
    
    printf("permon_solve: rank info after ownership count: nrows=%d nlocal=%d ilower=%d iupper=%d\n", nrows, nlocal, ilower, iupper);
    // -----------------------------
    // 2. Create the MPI matrix
    // -----------------------------
    PetscCall(MatCreate(comm, &A));

    PetscCall(MatSetType(A, MATMPIAIJ));
    PetscCall(MatSetSizes(A, nlocal, nlocal, PETSC_DECIDE, PETSC_DECIDE));
    // MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetFromOptions(A);
    PetscCall(MatSetUp(A));

    printf("permon_solve: rank info: nrows=%d nlocal=%d ilower=%d iupper=%d\n", nrows, nlocal, ilower, iupper);


    PetscInt irow;
    PetscInt *rcols = NULL;

    for (i = 0; i < nrows; i++) {
        // if (!owner[i]) continue;

        nnz  = rows_f[i+1] - rows_f[i];
        irow = globaldofs[i];

        PetscMalloc1(nnz, &rcols);

        for (PetscInt k = 0, j = rows_f[i]; j < rows_f[i+1]; j++, k++) {
            rcols[k] = globaldofs[cols_f[j-1] - 1];
        }

        MatSetValues(A, 1, &irow, nnz, rcols,
                &vals[rows_f[i] - 1],
                ADD_VALUES);

        PetscFree(rcols);
    }


    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "matrix.txt", &viewer);
    MatView(A, viewer);
    PetscViewerDestroy(&viewer);

    // -----------------------------
    // 4. Create vectors (b, c, x)
    // -----------------------------
    if (b_ptr) {
        PetscScalar *b_array = (PetscScalar*)b_ptr;
        PetscCall(VecCreateMPIWithArray(comm, 1, nlocal, PETSC_DECIDE, b_array, &b));
    } else {
        PetscCall(MatCreateVecs(A, &b, NULL));
    }

    if (c_ptr) {
        PetscScalar *c_array = (PetscScalar*)c_ptr;
        PetscCall(VecCreateMPIWithArray(comm, 1, nlocal, PETSC_DECIDE, c_array, &c));
    } else {
        PetscCall(MatCreateVecs(A, NULL, &c));
    }

    if (x_ptr) {
        PetscScalar *x_array = (PetscScalar*)x_ptr;
        PetscCall(VecCreateMPIWithArray(comm, 1, nlocal, PETSC_DECIDE, x_array, &x));
    } else {
        PetscCall(MatCreateVecs(A, NULL, &x));
    }



    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    * Setup QP: argmin 1/2 x'Ax -x'b s.t. c <= x
    *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(QPCreate(comm, &qp));
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
