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

    /* Dump rows array to a per-rank file for debugging */
    {
        int rank = 0;
        MPI_Comm_rank(comm, &rank);
        char fname[256];
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
    }
    /* Dump cols array to a per-rank file for debugging */
    // {
    //     int rank = 0;
    //     MPI_Comm_rank(comm, &rank);
    //     char fname[256];
    //     snprintf(fname, sizeof(fname), "cols_dump_rank%d.txt", rank);
    //     FILE *fp = fopen(fname, "w");
    //     if (fp) {
    //         fprintf(fp, "cols_local pointer=%p\n", (void*)cols_local);
    //         fprintf(fp, "cols_f pointer=%p\n", (void*)cols_f);
    //         fprintf(fp, "nloc=%d\n", nloc);
    //         fprintf(fp, "cols_f (0..nnz-1):\n");
    //         /* We don't have nnz here; attempt to print a portion around rows if possible */
    //         for (int rr = 0; rr < nloc; ++rr) {
    //             if (!owner[rr]) continue;
    //             int start = rows_f[rr];
    //             int end = rows_f[rr+1];
    //             fprintf(fp, "row %d nnz=%d:\n", rr, end - start);
    //             for (int j = start; j < end; ++j) {
    //                 fprintf(fp, "%d\n", cols_f[j-1]);
    //             }
    //         }
    //         fclose(fp);
    //     } else {
    //         fprintf(stderr, "Failed to open %s for writing\n", fname);
    //         fflush(stderr);
    //     }
    // }
    // TODO check if nrows == local_size
    for (i = 0; i < nrows; i++) {
        if (owner[i]) {
            if (globaldofs[i] < ilower) ilower = globaldofs[i];
            if (globaldofs[i] > iupper) iupper = globaldofs[i];
            nlocal++;
        } else {
            // printf("rank %p: skipping row %d (global %d) not owned\n", comm, i, globaldofs[i]);
        }
    }
    if (iupper == -1) { ilower = 0; iupper = -1; nlocal = 0; }  // rank owns nothing

    // -----------------------------
    // 2. Create the MPI matrix
    // -----------------------------
    PetscCall(MatCreate(comm, &A));

    // CHANGE THIS, hardcoded plus nrows might not be correct, should be nlocal (not sure if theyre the same)
    PetscCall(MatSetSizes(A, nlocal, nlocal, 4160, 4160));
    PetscCall(MatSetType(A, MATMPIAIJ));
    PetscCall(MatSetUp(A));

    printf("permon_solve: rank info: nrows=%d nlocal=%d ilower=%d iupper=%d\n", nrows, nlocal, ilower, iupper);
    int csize = 128;
    {
      int nnz,irow,i,j,k,*rcols;

      rcols = (int *)malloc( csize*sizeof(int) );
      // todo check if nrows == local_size
      for (i = 0; i < nrows ; i++) {
        
	    nnz = rows_f[i+1]-rows_f[i];
        if ( nnz>csize ) {
            csize = nnz+csize;
            int *tmp = (int *)realloc(rcols, csize*sizeof(int));
            if (!tmp) { perror("realloc failed"); exit(1); }
            rcols = tmp;
        }
        irow=globaldofs[i];
        for( k=0,j=rows_f[i]; j<rows_f[i+1]; j++,k++) {
	        rcols[k] = globaldofs[cols_f[j-1]-1];
        }
        MatSetValues(A, 1, &irow, nnz, rcols, &vals[rows_f[i]-1], ADD_VALUES);
        }
        free( rcols );
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "matrix.txt", &viewer);
    MatView(A, viewer);
    PetscViewerDestroy(&viewer);

    // PetscCall(MatView(A, PETSC_VIEWER_STDOUT_WORLD));
    // return 0;


    // PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    // PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

    // -----------------------------
    // 4. Create vectors (b, c, x)
    // -----------------------------
    if (b_ptr) {
        PetscScalar *b_array = (PetscScalar*)b_ptr;
        PetscCall(VecCreateMPIWithArray(comm, 1, nlocal, 4160, b_array, &b));
    } else {
        PetscCall(MatCreateVecs(A, &b, NULL));
    }

    if (c_ptr) {
        PetscScalar *c_array = (PetscScalar*)c_ptr;
        PetscCall(VecCreateMPIWithArray(comm, 1, nlocal, 4160, c_array, &c));
    } else {
        PetscCall(MatCreateVecs(A, NULL, &c));
    }

    if (x_ptr) {
        PetscScalar *x_array = (PetscScalar*)x_ptr;
        PetscCall(VecCreateMPIWithArray(comm, 1, nlocal, 4160, x_array, &x));
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
