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
int permon_solve(void *rows_local, void *cols_local, void *vals_local, int nrows, int ncols, void *b_ptr, void *c_ptr, void *x_ptr, int bound, int *globaldofs, int *owner, int *fcomm) {
    Vec       b, c, x;
    Vec       lb_fill = NULL, ub_fill = NULL;
    Vec       lb = NULL, ub = NULL;
    Mat       A;
    QP        qp;
    QPS       qps;
    PetscInt  i, rstart, rend;
    PetscBool converged, viewSol = PETSC_FALSE;
    PetscViewer viewer;
    
    MPI_Comm comm=MPI_Comm_f2c(*fcomm);

    printf("Communicator in permon_solve: %p\n", (void*)comm);
    fflush(stdout);


    /* Convert Fortran 1-based indices to C 0-based indices for PETSc */
    int *rows_f = (int*)rows_local;  /* Fortran 1-based row pointers */
    int *cols_f = (int*)cols_local;  /* Fortran 1-based column indices */
    double *vals = (double*)vals_local; /* matrix values */
    // int nrows_array = nrows + 1;  /* Rows array has nrows+1 elements */
    // int nnz = rows_f[nrows] - 1;  /* Number of nonzeros (last element - 1, in 1-based) */

    int ilower = 1000000000;
    int iupper = -1;

    for (i = 0; i < nrows; i++) {
        if (owner[i]) {
            if (globaldofs[i] < ilower) ilower = globaldofs[i];
            if (globaldofs[i] > iupper) iupper = globaldofs[i];
        }
    }
    /* handle rank owning no DOFs */
    if (iupper == -1) { ilower = 1; iupper = 0; }

    /* PETSc uses 0-based indexing normally */
    ilower--;  
    iupper--;


    // mprgp_print_vector(c_ptr, 15, "c");
    mprgp_print_vector(b_ptr, 15, "b");
    // mprgp_print_vector(x_ptr, 15, "x (initial)");
    
    /* Allocate temporary arrays for 0-based indices */
    // PetscInt *rows_c, *cols_c;
    // PetscMalloc1(nrows_array, &rows_c);
    // PetscMalloc1(nnz, &cols_c);
    
    // /* Convert row pointers: subtract 1 from each element */
    // for (i = 0; i < nrows_array; i++) {
    //     rows_c[i] = (PetscInt)(rows_f[i] - 1);
    // }
    
    // /* Convert column indices: subtract 1 from each element */
    // for (i = 0; i < nnz; i++) {
    //     cols_c[i] = (PetscInt)(cols_f[i] - 1);
    // }


    /* Create matrix directly from arrays (MatCreateSeqAIJWithArrays creates a new matrix) */
    // PetscCall(MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, nrows, ncols, rows_c, cols_c, vals_local, &A));
    MatCreate(comm, &A);
    MatSetType(A, MATAIJ);
    MatSetSizes(A,
            iupper-ilower+1,  // local rows
            iupper-ilower+1,  // local cols (square)
            4160,
            4160);  // global cols
    MatSetUp(A);

    for (i = 0; i < nrows; i++) {

        if (!owner[i]) continue;     // skip rows not owned

        int grow = globaldofs[i];   // PETSc row, Global dofs are already 0-based

        int start = rows_f[i]   - 1;  // convert Fortran 1-based → C 0-based
        int end   = rows_f[i+1] - 1;

        int nnz = end - start;

        /* build column array in global numbering */
        int colsPETSC[4096];   // (replace with malloc if needed)
        double valsPETSC[4096];
        
        for (int k = 0; k < nnz; k++) {
            colsPETSC[k] = globaldofs[ cols_f[start+k] - 1 ];  // Global dofs are already 0-based
            valsPETSC[k] = vals[start + k];
        }

        MatSetValues(A, 1, &grow, nnz, colsPETSC, valsPETSC, INSERT_VALUES);
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,   MAT_FINAL_ASSEMBLY);
    
{
  int myrank = -1;
  MPI_Comm_rank(comm, &myrank);

  /* make sure everyone reaches here (helps spot hangs) */
  MPI_Barrier(comm);

  /* Use PETSc synchronized printing so every rank's messages appear in order */
  PetscSynchronizedPrintf(comm, "=== DIAG: Hello from rank %d ===\n", myrank);

  /* Print PETSc ownership range for this matrix */
  PetscInt rstart=0, rend=0;
  MatGetOwnershipRange(A, &rstart, &rend);
  PetscSynchronizedPrintf(comm, "rank %d: PETSc ownership range = [%lld, %lld) owned = %lld  nslots(local)=%d\n",
                         myrank, (long long)rstart, (long long)rend, (long long)(rend-rstart), nrows);

  /* Print first few local slots (safe: guard by nrows) */
  for (int ii=0; ii < PetscMin(8, nrows); ++ii) {
    /* guard access if arrays might be NULL */
    int g = (globaldofs ? globaldofs[ii] : -999);
    int own = (owner ? owner[ii] : -1);
    double bval = (b_ptr ? ((double*)b_ptr)[ii] : 0.0);
    PetscSynchronizedPrintf(comm, "rank %d slot %d -> g=%d owner=%d bslot=%g\n",
                           myrank, ii, g, own, bval);
  }

  /* Check for NaN/Inf in matrix and print a norm */
  PetscReal anorm = 0.0;
  MatNorm(A, NORM_INFINITY, &anorm);
  PetscSynchronizedPrintf(comm, "rank %d: Mat infinity-norm = %g\n", myrank, (double)anorm);

  /* Check that globaldofs claimed by owner on this rank actually sit in PETSc range */
  for (int slot = 0; slot < 10; ++slot) {
    if (!owner[slot]) continue;
    int g = globaldofs[slot]; /* MUST be 0-based global index */
    if (g < (int)rstart || g >= (int)rend) {
      PetscSynchronizedPrintf(comm, "rank %d: WARNING slot %d -> g=%d not in PETSc range [%lld,%lld)\n",
                             myrank, slot, g, (long long)rstart, (long long)rend);
    }
  }

  /* Force flush in rank order so messages appear coherently */
  PetscSynchronizedFlush(comm, PETSC_STDOUT);
  MPI_Barrier(comm);  /* extra sync to ensure all ranks printed before proceeding */
}




    
    /* Create vectors from Fortran arrays: b from RHS, c from LowerLimit */
    /* VecCreateSeqWithArray wraps existing array data (doesn't copy) */
    if (b_ptr != NULL) {
        PetscScalar *b_array = (PetscScalar*)b_ptr;
        // PetscCall(VecCreateSeqWithArray(comm, 1, nrows, b_array, &b)); // does not copy the datam verify this
        PetscCall(VecCreateMPIWithArray(comm, 1, nrows, 4160, b_array, &b));
    } else {
        /* If b_ptr is NULL, create empty vector */
        PetscCall(MatCreateVecs(A, &b, NULL));
    }
    
    if (c_ptr != NULL) {
        PetscScalar *c_array = (PetscScalar*)c_ptr;
        // PetscCall(VecCreateSeqWithArray(comm, 1, nrows, c_array, &c));
        PetscCall(VecCreateMPIWithArray(comm, 1, nrows, 4160, c_array, &c));
        
    } else {
        /* If c_ptr is NULL, create empty vector */
        PetscCall(MatCreateVecs(A, NULL, &c));
    }

    if (x_ptr != NULL) {
        PetscScalar *x_array = (PetscScalar*)x_ptr;
        // PetscCall(VecCreateSeqWithArray(comm, 1, nrows, x_array, &x));
        PetscCall(VecCreateMPIWithArray(comm, 1, nrows, 4160, x_array, &x));
    } else {
        /* If x_ptr is NULL, create empty vector */
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
