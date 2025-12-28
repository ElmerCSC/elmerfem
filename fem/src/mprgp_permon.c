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
    // TODO check if nrows == local_size
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

    // CHANGE THIS, hardcoded plus nrows might not be correct, should be nlocal (not sure if theyre the same)
    // Is petsc_decide ok?
    PetscCall(MatSetType(A, MATMPIAIJ));
    PetscCall(MatSetSizes(A, nlocal, nlocal, PETSC_DECIDE, PETSC_DECIDE));
    // MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetFromOptions(A);
    PetscCall(MatSetUp(A));

    printf("permon_solve: rank info: nrows=%d nlocal=%d ilower=%d iupper=%d\n", nrows, nlocal, ilower, iupper);

    // new approach
    // PetscInt *d_nnz, *o_nnz;
    // PetscCall(PetscMalloc1(nlocal, &d_nnz));
    // PetscCall(PetscMalloc1(nlocal, &o_nnz));

    // for (PetscInt i = 0; i < nlocal; i++) {
    // d_nnz[i] = 0;
    // o_nnz[i] = 0;
    // }

    // for (PetscInt i = 0; i < nrows; i++) {
    // if (!owner[i]) continue;

    // PetscInt row = globaldofs[i];
    // PetscInt loc = row - ilower;
    // PetscInt nnz = rows_f[i+1] - rows_f[i];

    // for (PetscInt j = rows_f[i]; j < rows_f[i+1]; j++) {
    //     PetscInt col = globaldofs[cols_f[j-1]-1];
    //     if (col >= ilower && col <= iupper)
    //     d_nnz[loc]++;
    //     else
    //     o_nnz[loc]++;
    // }
    // }

    // MatMPIAIJSetPreallocation(A, 0, d_nnz, 0, o_nnz);
    // PetscCall(PetscFree(d_nnz));
    // PetscCall(PetscFree(o_nnz));

    // for (PetscInt i = 0; i < nrows; i++) {
    // PetscInt nnz = rows_f[i+1] - rows_f[i];
    // if (nnz == 0) continue;

    // PetscInt row = globaldofs[i];

    // PetscInt *cols_glob;
    // PetscScalar *vals_loc;

    // PetscCall(PetscMalloc1(nnz, &cols_glob));
    // PetscCall(PetscMalloc1(nnz, &vals_loc));

    // for (PetscInt k = 0, j = rows_f[i]; j < rows_f[i+1]; j++, k++) {
    //     cols_glob[k] = globaldofs[cols_f[j-1]-1];
    //     vals_loc[k] = vals[j-1];
    // }

    // PetscCall(MatSetValues(A, 1, &row, nnz, cols_glob, vals_loc, ADD_VALUES));

    // PetscCall(PetscFree(cols_glob));
    // PetscCall(PetscFree(vals_loc));
    // }

    // MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    // MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);



    // // with global mapping
    // PetscInt *lgmap;
    // PetscMalloc1(nlocal, &lgmap);

    // PetscInt cnt = 0;
    // for (i = 0; i < nrows; i++) {
    //     if (owner[i]) {
    //         lgmap[cnt++] = globaldofs[i];
    //     }
    // }

    // ISLocalToGlobalMapping l2g;
    // ISLocalToGlobalMappingCreate(comm, 1, nlocal, lgmap,
    //                             PETSC_OWN_POINTER, &l2g);

    // MatSetLocalToGlobalMapping(A, l2g, l2g);

    // cnt = 0;
    // int csize = 256;
    // int *rcols;
    // rcols = (int *)malloc( csize*sizeof(int) );

    // for (i = 0; i < nrows; i++) {
    //     if (!owner[i]) continue;

    //     PetscInt row = cnt++;  // LOCAL row index

    //     nnz = rows_f[i+1] - rows_f[i];
    //     for (PetscInt k = 0, j = rows_f[i]; j < rows_f[i+1]; j++, k++) {
    //         rcols[k] = globaldofs[cols_f[j-1] - 1];  // GLOBAL columns OK
    //     }

    //     PetscInt grow = lgmap[row];   // global row
    //     MatSetValues(A, 1, &grow, nnz, rcols,
    //                 &vals[rows_f[i]-1], ADD_VALUES);

    // }







    // int csize = 128;
    // {
    //   int nnz,irow,i,j,k,*rcols;

    //   rcols = (int *)malloc( csize*sizeof(int) );
    //   // todo check if nrows == local_size
    //   for (i = 0; i < nrows ; i++) {
    //     if (!owner[i]) continue;
        
    // 	nnz = rows_f[i+1]-rows_f[i];
    //     if ( nnz>csize ) {
    //         csize = nnz+csize;
    //         int *tmp = (int *)realloc(rcols, csize*sizeof(int));
    //         if (!tmp) { perror("realloc failed"); exit(1); }
    //         rcols = tmp;
    //     }
    //         irow=globaldofs[i];
    //         for( k=0,j=rows_f[i]; j<rows_f[i+1]; j++,k++) {
    //             rcols[k] = globaldofs[cols_f[j-1]-1];
    //         }

    //         /* Log the rcols being added to the matrix for indices 0..10 (if present)
    //          * Print the global row, number of nonzeros and up to the first 11 entries
    //          * with their corresponding values. This helps debug incorrect inserts. */
    //         {
    //             if (i < 10){
    //                 int tolog = (nnz < 11) ? nnz : 11;
    //                 printf("Mat insert: global row %d nnz=%d\n", (int)irow, (int)nnz);
    //                 for (int tt = 0; tt < tolog; ++tt) {
    //                     double v = vals[rows_f[i]-1 + tt];
    //                     printf("  entry %2d: col=%d val=%.18e\n", tt, rcols[tt], v);
    //                 }
    //                 if (nnz > 11) printf("  ... (%d more entries)\n", (int)(nnz-11));
    //                 fflush(stdout);
    //             }
 
    //         }
            
    //         MatSetValues(A, 1, &irow, nnz, rcols, &vals[rows_f[i]-1], ADD_VALUES);

    //     }
    //     free( rcols );
    // }

    // for (int i = 0; i < nrows; i++) {
    //     PetscInt irow = globaldofs[i];      /* global row index */

    //     /* if you want only owners to insert, uncomment the next line */
    //     // if (!owner[i]) continue; 

    //     /* convert Fortran 1-based row pointers into C 0-based indices */
    //     for (int kk = rows_f[i] - 1; kk < rows_f[i+1] - 1; kk++) {
    //         PetscInt local_col = cols_f[kk] - 1;   /* 1-based -> 0-based local col */
    //         PetscInt J = globaldofs[local_col];    /* global column index */
    //         PetscScalar v = vals[kk];              /* value at that entry */

    //         MatSetValue(A, irow, J, v, ADD_VALUES);
    //     }
    // }

    PetscInt irow;
    PetscInt *rcols = NULL;

    for (i = 0; i < nrows; i++) {
        if (!owner[i]) continue;

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


    // MatView(A, PETSC_VIEWER_STDOUT_SELF);

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


/* ---------- PETSc snippet (compile & run under MPI) ---------- */
/* Assumptions:
   - MPI initialized by Fortran and passed to C (comm available).
   - arrays available: local_size_all (number of local rows in arrays),
     owner[local_size_all] (0/1), globaldofs[local_size_all],
     rows[local_size_all+1], cols[nnz_total], vals[nnz_total].
*/


// void build_and_assemble(Mat *A_out, MPI_Comm comm,
//                         int local_size_all, int *owner, int *globaldofs,
//                         int *rows, int *cols, double *vals, int nvals_total,
//                         int global_nrows) {
//   PetscErrorCode ierr;
//   Mat A;
//   int i;

//   /* 1) Count how many actual owned rows this rank will insert */
//   int owned = 0;
//   for (i=0;i<local_size_all;i++) if (owner[i]) owned++;
//   /* Create PETSc matrix: tell PETSc how many local rows we own (owned) */
//   ierr = MatCreateMPIAIJ(comm, owned, owned, global_nrows, global_nrows, 
//                          0, NULL, 0, NULL, &A); CHKERRABORT(PETSC_COMM_SELF, ierr);

//   /* 2) Compute d_nnz and o_nnz arrays (per owned row) */
//   int *d_nnz = (int*)malloc(sizeof(int)*owned);
//   int *o_nnz = (int*)malloc(sizeof(int)*owned);

//   /* get PETSc ownership range assigned to this rank */
//   PetscInt rstart, rend;
//   ierr = MatGetOwnershipRange(A, &rstart, &rend); CHKERRABORT(PETSC_COMM_SELF, ierr);

//   int rowcount = 0;
//   for (i=0;i<local_size_all;i++) {
//     if (!owner[i]) continue;
//     int nnz = rows[i+1] - rows[i];
//     int dn=0, on=0, j;
//     for (j = rows[i]; j < rows[i+1]; ++j) {
//       int local_col_idx = cols[j-1] - 1;       /* if cols were 1-based like in your code */
//       int gcol = globaldofs[ local_col_idx ];
//       if (gcol >= rstart && gcol < rend) dn++; else on++;
//     }
//     d_nnz[rowcount] = dn;
//     o_nnz[rowcount] = on;
//     rowcount++;
//   }

//   /* 3) Preallocate using per-row arrays */
//   ierr = MatMPIAIJSetPreallocation(A, 0, NULL, 0, NULL); /* default zero — will be replaced */
//   /* PETSc expects global preallocation before any MatSetValues,
//      but convenience: call MatCreateMPIAIJ with zeros above then call MatMPIAIJSetPreallocationCSR
//      only if providing ai/aj. Simpler is MatSeqAIJSetPreallocation for seq, but here we use:
//   */
//   /* Alternatively, provide uniform estimate: */
//   /* ierr = MatMPIAIJSetPreallocation(A, max_diag_estimate, diag_counts, max_offdiag_estimate, offdiag_counts); */
//   /* For simplicity in this snippet we skip per-row preallocation call because
//      MatSetValues will still work; in production you should call MatMPIAIJSetPreallocation
//      with d_nnz/o_nnz arrays or construct CSR and use MatCreateMPIAIJWithArrays. */

//   /* 4) Insert values: iterate owned rows, map cols -> global columns and call MatSetValues */
//   rowcount = 0;
//   for (i=0;i<local_size_all;i++) {
//     if (!owner[i]) continue;
//     int grow = globaldofs[i];     /* global row number of this local row */
//     int nnz = rows[i+1] - rows[i];
//     int *gcols = (int*)malloc(sizeof(int)*nnz);
//     int j,k=0;
//     for (j=rows[i]; j<rows[i+1]; ++j,++k) {
//       int lcol = cols[j-1]-1;            /* convert 1-based->0-based if needed */
//       gcols[k] = globaldofs[lcol];
//     }
//     ierr = MatSetValues(A, 1, &grow, nnz, gcols, &vals[rows[i]-1], ADD_VALUES);
//     free(gcols);
//     rowcount++;
//   }

//   /* 5) Finalize */
//   ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//   ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);

//   free(d_nnz); free(o_nnz);
//   *A_out = A;
// }

/* ------------------------------------------------------------------ */
