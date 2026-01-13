/* Simple test helper to print CRS row pointers from C
 * This function is called from Fortran for debugging.
 */
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <permonqps.h>
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
    return PermonInitialize(NULL, NULL, "permonrc", NULL);
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
    PetscViewer viewer;
    
    MPI_Comm comm=MPI_Comm_f2c(fcomm);
    
    /* Print communicator/rank info for debugging communicator mismatches */
    {
        int csize = -1, cmp = -1;
        int rank;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &csize);
        MPI_Comm_compare(comm, MPI_COMM_WORLD, &cmp);
        printf("permon_solve: C side: comm rank=%d size=%d compare_with_WORLD=%d\n", rank, csize, cmp);
        fflush(stdout);
    }
    int *rows_f = (int*)rows_local;
    int *cols_f = (int*)cols_local;

    double *vals = (double*)vals_local;

    /* Debug: print a small window of the RHS (b_ptr) as seen in C to
       verify that Fortran passed the expected values. Fortran arrays are
       1-based; here we print entries 60..70 (1-based) if available. */
    if (b_ptr) {
        int rank = 0;
        MPI_Comm_rank(comm, &rank);
        int start = 60, end = 70;
        if (rank == 0) {
            if (nrows < start) {
                printf("DBG: C side rank %d: nrows=%d < %d, skipping b[60..70]\n", rank, nrows, start);
            } else {
                if (nrows < end) end = nrows;
                printf("DBG: C side rank %d: b_ptr entries (1-based) %d..%d:\n", rank, start, end);
                for (int idx = start; idx <= end; ++idx) {
                    double val = ((double*)b_ptr)[idx-1];
                    printf("  b[%d] = %.18e\n", idx, val);
                }
            }
            fflush(stdout);
        }
    }

    /* Debug: print a small window of the constraint vector (c_ptr) as seen
       in C to verify Fortran->C passing. Use same indices as for `b_ptr`. */
    if (c_ptr) {
        int rank = 0;
        MPI_Comm_rank(comm, &rank);
        int start = 60, end = 70;
        if (rank == 0) {
            if (nrows < start) {
                printf("DBG: C side rank %d: nrows=%d < %d, skipping c[60..70]\n", rank, nrows, start);
            } else {
                if (nrows < end) end = nrows;
                printf("DBG: C side rank %d: c_ptr entries (1-based) %d..%d:\n", rank, start, end);
                for (int idx = start; idx <= end; ++idx) {
                    double val = ((double*)c_ptr)[idx-1];
                    printf("  c[%d] = %.18e\n", idx, val);
                }
            }
            fflush(stdout);
        }
    }

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

        /* c_ptr dump (constraint vector passed from Fortran) */
        snprintf(fname, sizeof(fname), "c_ptr_dump_rank%d.txt", rank);
        fp = fopen(fname, "w");
        if (fp) {
            fprintf(fp, "c_ptr pointer=%p\n", (void*)c_ptr);
            fprintf(fp, "c_ptr (1..nrows):\n");
            if (c_ptr) {
                for (i = 0; i < nrows; ++i) {
                    fprintf(fp, "%.18e\n", ((double*)c_ptr)[i]);
                }
            } else {
                fprintf(fp, "<NULL c_ptr>\n");
            }
            fclose(fp);
        } else {
            fprintf(stderr, "Failed to open %s for writing\n", fname);
            fflush(stderr);
        }
    }
    /* Debug checks: globaldofs min/max across ranks, sample entries, owner samples,
       and a couple of translated CRS rows for quick verification. */
    {
        int local_min = INT_MAX, local_max = -1;
        for (i = 0; i < nrows; ++i) {
            if (globaldofs[i] < local_min) local_min = globaldofs[i];
            if (globaldofs[i] > local_max) local_max = globaldofs[i];
        }
        int global_min = 0, global_max = 0;
        MPI_Allreduce(&local_min, &global_min, 1, MPI_INT, MPI_MIN, comm);
        MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, comm);

        int rank = -1;
        MPI_Comm_rank(comm, &rank);
        printf("DBG: rank %d globaldofs local_min/local_max=%d/%d global_min/global_max=%d/%d\n",
               rank, local_min, local_max, global_min, global_max);

        if (nrows > 0) {
            int mid = nrows / 2;
            printf("DBG: rank %d sample globaldofs[0]=%d globaldofs[mid]=%d\n",
                   rank, globaldofs[0], globaldofs[mid]);
        }

        /* print some owner samples (likely 0/1 flags) */
        int ns = nrows < 10 ? nrows : 10;
        printf("DBG: rank %d owner samples:", rank);
        for (int kk = 0; kk < ns; ++kk) printf(" %d", owner[kk]);
        printf("\n");

        /* Translate and print a couple of CRS rows (first and middle) */
        if (nrows > 0) {
            int rows_to_check[2];
            rows_to_check[0] = 0;
            rows_to_check[1] = nrows > 1 ? nrows / 2 : 0;
            for (int rr = 0; rr < 2; ++rr) {
                int idx = rows_to_check[rr];
                int start = rows_f[idx];
                int end = rows_f[idx+1];
                int local_nnz = end - start;
                printf("DBG: rank %d row idx=%d (rows_f start/end=%d/%d nnz=%d):\n", rank, idx, start, end, local_nnz);
                for (int jj = start; jj < end; ++jj) {
                    int col_local = cols_f[jj-1] - 1; /* Fortran->C */
                    int mapped = -999999;
                    if (col_local >= 0 && col_local < nrows) mapped = globaldofs[col_local];
                    printf("  col_local=%d mapped_global=%d raw_col_entry=%d\n", col_local, mapped, cols_f[jj-1]);
                }
            }
        }
    }

    // TODO probably not needed
    for (i = 0; i < nrows; i++) {
        if (owner[i]) {
            if (globaldofs[i] < ilower) ilower = globaldofs[i];
            if (globaldofs[i] > iupper) iupper = globaldofs[i];
            nlocal++;
        }
    }
    if (iupper == -1) { ilower = 0; iupper = -1; }
    
    
    printf("permon_solve: rank info after ownership count: nrows=%d nlocal=%d ilower=%d iupper=%d\n", nrows, nlocal, ilower, iupper);

    // TODO maybe can be done simpler (look at previous code versions)
    ISLocalToGlobalMapping ltog;
    /* PETSc's PetscInt may differ from C `int` (32 vs 64 bit). Create a
       temporary PetscInt copy of the Fortran-provided `globaldofs` so we
       pass correctly-typed data to PETSc and avoid undefined behaviour. */
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
    PetscCall(MatCreate(comm, &A));

    PetscCall(MatSetType(A, MATMPIAIJ));
    PetscCall(MatSetSizes(A, nlocal, nlocal, PETSC_DECIDE, PETSC_DECIDE));
    // MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetFromOptions(A);
    PetscCall(MatSetUp(A));
    /* Attach local->global mapping so we can insert with local indices */
    PetscCall(MatSetLocalToGlobalMapping(A, ltog, ltog));
    printf("permon_solve: rank info: nrows=%d nlocal=%d ilower=%d iupper=%d\n", nrows, nlocal, ilower, iupper);


    PetscInt irow;


    for (i = 0; i < nrows; i++) {
        nnz  = rows_f[i+1] - rows_f[i];
        PetscInt iloc = i; /* local row index */
        PetscInt *rcols;
        PetscMalloc1(nnz, &rcols);

        for (PetscInt k = 0, j = rows_f[i]; j < rows_f[i+1]; j++, k++) {
            /* Local column index (Fortran->C): cols_f[j-1]-1 */
            rcols[k] = cols_f[j-1] - 1;
        }

        PetscCall(MatSetValuesLocal(A, 1, &iloc, nnz, rcols,
                &vals[rows_f[i] - 1], ADD_VALUES));

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
    /* Create vectors with the same layout as A */
    // PetscCall(MatCreateVecs(A, &x, &b));
    // PetscCall(VecDuplicate(x, &c));

    // /* Initialize (important!) */
    // // PetscCall(VecSet(x, 0.0));
    // // PetscCall(VecSet(b, 0.0));
    // // PetscCall(VecSet(c, 0.0));

    // for (PetscInt i = 0; i < nrows; i++) {
    //     PetscScalar v = ((PetscScalar*)b_ptr)[i];
    //     PetscInt    iloc = i;
    //     PetscCall(VecSetValuesLocal(b, 1, &iloc, &v, INSERT_VALUES));
    // }

    // PetscCall(VecAssemblyBegin(b));
    // PetscCall(VecAssemblyEnd(b));

    // for (PetscInt i = 0; i < nrows; i++) {
    //     PetscScalar v = ((PetscScalar*)x_ptr)[i];
    //     PetscInt    iloc = i;
    //     PetscCall(VecSetValuesLocal(x, 1, &iloc, &v, INSERT_VALUES));
    // }

    // PetscCall(VecAssemblyBegin(x));
    // PetscCall(VecAssemblyEnd(x));

    // for (PetscInt i = 0; i < nrows; i++) {
    //     PetscScalar v = ((PetscScalar*)c_ptr)[i];
    //     PetscInt    iloc = i;
    //     PetscCall(VecSetValuesLocal(c, 1, &iloc, &v, INSERT_VALUES));
    // }

    // PetscCall(VecAssemblyBegin(c));
    // PetscCall(VecAssemblyEnd(c));




    /* Create vectors with the same layout as A and fill them from local arrays
       using VecSetValuesLocal so PETSc maps local ordering correctly. */
    PetscCall(MatCreateVecs(A, &x, &b));
    PetscCall(VecDuplicate(x, &c));

    if (b_ptr) {
        for (PetscInt ii = 0; ii < nrows; ++ii) {
            // if (!owner[ii]) continue; /* only owners set global b entries */
            PetscScalar v = ((PetscScalar*)b_ptr)[ii];
            PetscInt g = globaldofs[ii];
            PetscCall(VecSetValues(b, 1, &g, &v, ADD_VALUES));
        }
        PetscCall(VecAssemblyBegin(b));
        PetscCall(VecAssemblyEnd(b));
    }

    if (x_ptr) {
        for (PetscInt ii = 0; ii < nrows; ++ii) {
            // if (!owner[ii]) continue; /* only owners set initial x entries */
            PetscScalar v = ((PetscScalar*)x_ptr)[ii];
            PetscInt g = globaldofs[ii];
            PetscCall(VecSetValues(x, 1, &g, &v, INSERT_VALUES));
        }
        PetscCall(VecAssemblyBegin(x));
        PetscCall(VecAssemblyEnd(x));
    }

    if (c_ptr) {
        for (PetscInt ii = 0; ii < nrows; ++ii) {
            // if (!owner[ii]) continue; /* only owners set constraint vectors */
            PetscScalar v = ((PetscScalar*)c_ptr)[ii];
            PetscInt g = globaldofs[ii];
            PetscCall(VecSetValues(c, 1, &g, &v, INSERT_VALUES));
        }
        PetscCall(VecAssemblyBegin(c));
        PetscCall(VecAssemblyEnd(c));
    }



    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    * Setup QP: argmin 1/2 x'Ax -x'b s.t. c <= x
    *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /* Debug: report norms of A and b before creating QP */
    {
        PetscReal bnorm = 0.0, Anorm = 0.0;
        PetscCall(VecNorm(b, NORM_2, &bnorm));
        PetscCall(MatNorm(A, NORM_FROBENIUS, &Anorm));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "DBG: Pre-QP ||b||_2 = %22.16e  ||A||_F = %22.16e\n", (double)bnorm, (double)Anorm));
    }

    /* Dump assembled RHS vector b for comparison (per-rank ASCII) */
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "DBG: Writing assembled RHS to b_parallel.txt\n"));
        PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, "b_parallel.txt", &viewer));
        PetscCall(VecView(b, viewer));
        PetscCall(PetscViewerDestroy(&viewer));
    }

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

    /* Compute and print the 2-norm of the solution vector x */
    {
        PetscReal xnorm = 0.0;
        PetscCall(VecNorm(x, NORM_2, &xnorm));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "permon_solve: Solution ||x||_2 = %22.16e\n", (double)xnorm));
    }

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
    // PetscCall(ISLocalToGlobalMappingDestroy(&ltog));

    PetscCall(QPSDestroy(&qps));
    PetscCall(QPDestroy(&qp));

    return 0;
}
