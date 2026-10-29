 /* C wrapper around NVIDIA's cuDSS GPU sparse direct solver, callable
 * from Fortran (see CUDSS_SolveSystem in DirectSolve.F90). 
 * Three entry points (factorize, solve, free).
 *
 * Both real (CUDSS_R_64F) and complex (CUDSS_C_64F) double precision are
 * supported. The two are separate Fortran entry points -- cudss_ffactorize /
 * cudss_fsolve and cudss_zfactorize / cudss_zsolve. They share the implementation
 * below and differ only in the element size and the cuDSS value type. The caller
 * passes a matrix in *complex* CSR for the z entries (n complex rows, not the
 * 2n real rows Elmer stores internally).
 *
 * For Single-GPU and 1 MPI task. The CSR matrix is copied to device memory once, at
 * factorize time; only the right-hand side and solution are transferred on
 * every solve.
 *
 * For Multi-GPU Multi-node setup the additional entry points are 
 * cudss_dfactorize_par (double precision real), cudss_zdfactorize_par (double
 * precision complex).
 */

#include "../config.h"
#ifdef HAVE_CUDSS

#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cudss.h>
#ifdef HAVE_MPI
#include <mpi.h>
#endif

typedef struct {
  cudssHandle_t handle;
  cudssConfig_t config;
  cudssData_t   data;
  cudssMatrix_t Amat, bmat, xmat;

  int    n;
  size_t esz;  /* element size in bytes: 8 real, 16 complex */

  int   *d_rows;
  int   *d_cols;
  void  *d_vals;
  void  *d_b;
  void  *d_x;

#ifdef HAVE_MPI
  /* MGMN only: kept alive for the whole factorization because the cuDSS
   * communication layer holds a pointer to this MPI_Comm, not a copy. */
  MPI_Comm comm;
#endif
} ElmerCUDSS;

static int cuda_ok(cudaError_t st, const char *what)
{
  if (st != cudaSuccess) {
    fprintf(stderr, "CUDSS_SolveSystem: CUDA call '%s' failed: %s\n",
            what, cudaGetErrorString(st));
    return 0;
  }
  return 1;
}

static int cudss_ok(cudssStatus_t st, const char *what)
{
  if (st != CUDSS_STATUS_SUCCESS) {
    fprintf(stderr, "CUDSS_SolveSystem: cuDSS call '%s' failed with status %d\n",
            what, (int)st);
    return 0;
  }
  return 1;
}

static void cudss_teardown(ElmerCUDSS *ctx)
{
  if (!ctx) return;

  if (ctx->Amat) cudssMatrixDestroy(ctx->Amat);
  if (ctx->bmat) cudssMatrixDestroy(ctx->bmat);
  if (ctx->xmat) cudssMatrixDestroy(ctx->xmat);
  if (ctx->data) cudssDataDestroy(ctx->handle, ctx->data);
  if (ctx->config) cudssConfigDestroy(ctx->config);
  if (ctx->handle) cudssDestroy(ctx->handle);

  if (ctx->d_rows) cudaFree(ctx->d_rows);
  if (ctx->d_cols) cudaFree(ctx->d_cols);
  if (ctx->d_vals) cudaFree(ctx->d_vals);
  if (ctx->d_b) cudaFree(ctx->d_b);
  if (ctx->d_x) cudaFree(ctx->d_x);

  free(ctx);
}

/* mtype: 0 = general (nonsymmetric), 2 = symmetric positive definite. The
 * caller (CUDSS_SolveSystem) is responsible for passing only the upper
 * triangle (including the diagonal) of rows/cols/vals when mtype is 2.
 *
 * There is no value 1 (CUDSS_MTYPE_SYMMETRIC) implemented.
 */
static ElmerCUDSS *cudss_factorize_impl(int n, int nnz, const int *rows,
                                        const int *cols, const void *vals,
                                        int mtype, int is_complex)
{
  ElmerCUDSS *ctx;
  cudssMatrixType_t mt;
  cudssMatrixViewType_t mv;
  cudssDataType_t vt;

  ctx = (ElmerCUDSS *)calloc(1, sizeof(ElmerCUDSS));
  if (!ctx) {
    fprintf(stderr, "CUDSS_SolveSystem: out of host memory\n");
    return NULL;
  }
  ctx->n = n;
  ctx->esz = is_complex ? 2 * sizeof(double) : sizeof(double);
  vt = is_complex ? CUDSS_C_64F : CUDSS_R_64F;

  switch (mtype) {
    case 2:  
      mt = CUDSS_MTYPE_SPD;     
      mv = CUDSS_MVIEW_UPPER; 
      break;
    default: 
      mt = CUDSS_MTYPE_GENERAL; 
      mv = CUDSS_MVIEW_FULL;  
      break;
  }

  if (!cuda_ok(cudaMalloc((void **)&ctx->d_rows, (size_t)(n + 1) * sizeof(int)), "cudaMalloc(rows)") ||
      !cuda_ok(cudaMalloc((void **)&ctx->d_cols, (size_t)nnz * sizeof(int)), "cudaMalloc(cols)") ||
      !cuda_ok(cudaMalloc(&ctx->d_vals, (size_t)nnz * ctx->esz), "cudaMalloc(vals)") ||
      !cuda_ok(cudaMalloc(&ctx->d_b, (size_t)n * ctx->esz), "cudaMalloc(b)") ||
      !cuda_ok(cudaMalloc(&ctx->d_x, (size_t)n * ctx->esz), "cudaMalloc(x)")) {
    cudss_teardown(ctx);
    return NULL;
  }

  {
    /* convert 1-based (fortran) indices to 0-based (C), aka zero rows */
    int *z_rows = (int *)malloc((size_t)(n + 1) * sizeof(int));
    int *z_cols = (int *)malloc((size_t)nnz * sizeof(int));
    int ok, i;

    if (!z_rows || !z_cols) {
      fprintf(stderr, "CUDSS_SolveSystem: out of host memory for index conversion\n");
      free(z_rows);
      free(z_cols);
      cudss_teardown(ctx);
      return NULL;
    }

    for (i = 0; i <= n; i++) z_rows[i] = rows[i] - 1;
    for (i = 0; i < nnz; i++) z_cols[i] = cols[i] - 1;

    ok = cuda_ok(cudaMemcpy(ctx->d_rows, z_rows, (size_t)(n + 1) * sizeof(int), cudaMemcpyHostToDevice), "memcpy(rows)") &&
         cuda_ok(cudaMemcpy(ctx->d_cols, z_cols, (size_t)nnz * sizeof(int), cudaMemcpyHostToDevice), "memcpy(cols)") &&
         cuda_ok(cudaMemcpy(ctx->d_vals, vals, (size_t)nnz * ctx->esz, cudaMemcpyHostToDevice), "memcpy(vals)");

    free(z_rows);
    free(z_cols);

    if (!ok) {
      cudss_teardown(ctx);
      return NULL;
    }
  }

  if (!cudss_ok(cudssCreate(&ctx->handle), "cudssCreate") ||
      !cudss_ok(cudssConfigCreate(&ctx->config), "cudssConfigCreate") ||
      !cudss_ok(cudssDataCreate(ctx->handle, &ctx->data), "cudssDataCreate")) {
    cudss_teardown(ctx);
    return NULL;
  }

  /* Standard (gap-free) CSR: row i's entries run from d_rows[i] to
   * d_rows[i+1]-1, which is the only layout cuDSS supports -- rowEnd must be
   * NULL for it (a non-NULL rowEnd asks for the gapped CSR-3 variant). The
   * index arrays were rebased to zero above, so CUDSS_BASE_ZERO applies. */
  if (!cudss_ok(cudssMatrixCreateCsr(&ctx->Amat, n, n, nnz,
                    ctx->d_rows, NULL, ctx->d_cols, ctx->d_vals,
                    CUDSS_R_32I, CUDSS_R_32I, vt, mt, mv, CUDSS_BASE_ZERO),
                "cudssMatrixCreateCsr") ||
      !cudss_ok(cudssMatrixCreateDn(&ctx->bmat, n, 1, n, ctx->d_b, vt,
                    CUDSS_LAYOUT_COL_MAJOR), "cudssMatrixCreateDn(b)") ||
      !cudss_ok(cudssMatrixCreateDn(&ctx->xmat, n, 1, n, ctx->d_x, vt,
                    CUDSS_LAYOUT_COL_MAJOR), "cudssMatrixCreateDn(x)")) {
    cudss_teardown(ctx);
    return NULL;
  }

  if (!cudss_ok(cudssExecute(ctx->handle, CUDSS_PHASE_ANALYSIS, ctx->config, ctx->data,
                    ctx->Amat, ctx->xmat, ctx->bmat), "cudssExecute(ANALYSIS)") ||
      !cudss_ok(cudssExecute(ctx->handle, CUDSS_PHASE_FACTORIZATION, ctx->config, ctx->data,
                    ctx->Amat, ctx->xmat, ctx->bmat), "cudssExecute(FACTORIZATION)")) {
    cudss_teardown(ctx);
    return NULL;
  }

  return ctx;
}

static int cudss_solve_impl(ElmerCUDSS *ctx, int n, void *x, const void *b)
{
  if (!ctx) {
    fprintf(stderr, "CUDSS_SolveSystem: solve called without a factorization\n");
    return 0;
  }

  /* Guards against a real/complex mix-up: the complex path must hand us the
   * complex row count n, not Elmer's 2n real row count. */
  if (n != ctx->n) {
    fprintf(stderr, "CUDSS_SolveSystem: solve called with n=%d but the "
                    "factorization is of order %d\n", n, ctx->n);
    return 0;
  }

  if (!cuda_ok(cudaMemcpy(ctx->d_b, b, (size_t)n * ctx->esz, cudaMemcpyHostToDevice),
               "memcpy(b)"))
    return 0;

  if (!cudss_ok(cudssExecute(ctx->handle, CUDSS_PHASE_SOLVE, ctx->config, ctx->data,
                    ctx->Amat, ctx->xmat, ctx->bmat), "cudssExecute(SOLVE)"))
    return 0;

  return cuda_ok(cudaMemcpy(x, ctx->d_x, (size_t)n * ctx->esz, cudaMemcpyDeviceToHost),
                 "memcpy(x)");
}

ElmerCUDSS *FC_FUNC_(cudss_ffactorize,CUDSS_FFACTORIZE)
    (int *n, int *nnz, int *rows, int *cols, double *vals, int *mtype)
{
  return cudss_factorize_impl(*n, *nnz, rows, cols, vals, *mtype, 0);
}

/* vals holds 2*nnz doubles: interleaved (re,im) pairs, i.e. Fortran
 * COMPLEX(KIND=dp), which is layout-compatible with cuDoubleComplex. */
ElmerCUDSS *FC_FUNC_(cudss_zfactorize,CUDSS_ZFACTORIZE)
    (int *n, int *nnz, int *rows, int *cols, double *vals, int *mtype)
{
  return cudss_factorize_impl(*n, *nnz, rows, cols, vals, *mtype, 1);
}

void FC_FUNC_(cudss_fsolve,CUDSS_FSOLVE)(ElmerCUDSS **handle, int *n, double *x,
                                         double *b, int *status)
{
  *status = cudss_solve_impl(*handle, *n, x, b);
}

void FC_FUNC_(cudss_zsolve,CUDSS_ZSOLVE)(ElmerCUDSS **handle, int *n, double *x,
                                         double *b, int *status)
{
  *status = cudss_solve_impl(*handle, *n, x, b);
}

void FC_FUNC_(cudss_ffree,CUDSS_FFREE)(ElmerCUDSS **handle)
{
  cudss_teardown(*handle);
  *handle = NULL;
}

/* ------------------------------------------------------------------------- *
 * CuDSS multi-GPU multi-node mode
 *
 * Each rank is assumed hold complete fully assembled rows 
 * and there is no overlap and numbering must be contiguous.
 *
 * Reuses cudss_fsolve() and cudss_ffree() 
 * ------------------------------------------------------------------------- */
#ifdef HAVE_MPI

/* 
 * Distribute devices on node to node local ranks in round-robin manner.
 * */
static void cudss_bind_device(MPI_Comm comm)
{
  MPI_Comm shmcomm;
  int grank, lrank, ndev = 0;

  MPI_Comm_rank(comm, &grank);
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, grank, MPI_INFO_NULL, &shmcomm);
  MPI_Comm_rank(shmcomm, &lrank);
  MPI_Comm_free(&shmcomm);

  if (cudaGetDeviceCount(&ndev) == cudaSuccess && ndev > 0)
    cudaSetDevice(lrank % ndev);
}

/* nglob/nnz_glob: global order and global nonzero count (cuDSS descriptor
 * metadata, identical on every rank). nloc/nnzloc: the owned rows and their
 * nonzeros this rank actually provides. rows: local CSR offsets (nloc+1,
 * 1-based from Fortran). cols: GLOBAL column indices (1-based). first_row/
 * last_row: 0-based inclusive global range of the owned block (empty rank:
 * first_row > last_row). */
static ElmerCUDSS *cudss_dfactorize_par_impl(int nglob, int nnz_glob, int nloc,
    int nnzloc, int first_row, int last_row, const int *rows, const int *cols,
    const void *vals, int mtype, int is_complex, MPI_Comm comm)
{
  ElmerCUDSS *ctx;
  cudssMatrixType_t mt;
  cudssMatrixViewType_t mv;
  cudssDataType_t vt = is_complex ? CUDSS_C_64F : CUDSS_R_64F;

  ctx = (ElmerCUDSS *)calloc(1, sizeof(ElmerCUDSS));
  if (!ctx) { fprintf(stderr, "CUDSS_SolveSystem: out of host memory\n"); return NULL; }
  ctx->n    = nloc;              /* the solve works on the owned block */
  ctx->esz  = is_complex ? 2 * sizeof(double) : sizeof(double);
  ctx->comm = comm;

  switch (mtype) {
    case 2:  
      mt = CUDSS_MTYPE_SPD;
      mv = CUDSS_MVIEW_UPPER;
      break;
    default: 
      mt = CUDSS_MTYPE_GENERAL; 
      mv = CUDSS_MVIEW_FULL;
      break;
  }

  cudss_bind_device(comm);

  /* nloc/nnzloc can be 0 on an empty partition; keep device allocations >0. */
  if (!cuda_ok(cudaMalloc((void **)&ctx->d_rows, (size_t)(nloc + 1) * sizeof(int)), "cudaMalloc(rows)") ||
      !cuda_ok(cudaMalloc((void **)&ctx->d_cols, (size_t)(nnzloc > 0 ? nnzloc : 1) * sizeof(int)), "cudaMalloc(cols)") ||
      !cuda_ok(cudaMalloc(&ctx->d_vals, (size_t)(nnzloc > 0 ? nnzloc : 1) * ctx->esz), "cudaMalloc(vals)") ||
      !cuda_ok(cudaMalloc(&ctx->d_b, (size_t)(nloc > 0 ? nloc : 1) * ctx->esz), "cudaMalloc(b)") ||
      !cuda_ok(cudaMalloc(&ctx->d_x, (size_t)(nloc > 0 ? nloc : 1) * ctx->esz), "cudaMalloc(x)")) {
    cudss_teardown(ctx);
    return NULL;
  }

  /* Modify local rows and offsets to 0-based indexing and copy to device. 
   * TODO: Modify on device to get rid of extra memory copy. */
  {
    int *z_rows = (int *)malloc((size_t)(nloc + 1) * sizeof(int));
    int *z_cols = (int *)malloc((size_t)(nnzloc > 0 ? nnzloc : 1) * sizeof(int));
    int ok, i;

    if (!z_rows || !z_cols) {
      fprintf(stderr, "CUDSS_SolveSystem: out of host memory for index conversion\n");
      free(z_rows); free(z_cols); cudss_teardown(ctx);
      return NULL;
    }

    for (i = 0; i <= nloc; i++)  z_rows[i] = rows[i] - 1;
    for (i = 0; i < nnzloc; i++) z_cols[i] = cols[i] - 1;

    ok = cuda_ok(cudaMemcpy(ctx->d_rows, z_rows, (size_t)(nloc + 1) * sizeof(int), cudaMemcpyHostToDevice), "memcpy(rows)") &&
         (nnzloc == 0 || cuda_ok(cudaMemcpy(ctx->d_cols, z_cols, (size_t)nnzloc * sizeof(int), cudaMemcpyHostToDevice), "memcpy(cols)")) &&
         (nnzloc == 0 || cuda_ok(cudaMemcpy(ctx->d_vals, vals, (size_t)nnzloc * ctx->esz, cudaMemcpyHostToDevice), "memcpy(vals)"));

    free(z_rows); free(z_cols);
    if (!ok) { cudss_teardown(ctx); return NULL; }
  }

  if (!cudss_ok(cudssCreate(&ctx->handle), "cudssCreate") ||
      !cudss_ok(cudssConfigCreate(&ctx->config), "cudssConfigCreate") ||
      !cudss_ok(cudssDataCreate(ctx->handle, &ctx->data), "cudssDataCreate")) {
    cudss_teardown(ctx);
    return NULL;
  }

  /* The communication layer library name is determined compile time as macro
   * definition ELMER_CUDSS_COMMLAYER_NAME. If that is NULL, then fall back to
   * cuDSS default behavior: read the communication layer library name from
   * CUDSS_COMM_LIB environment variable.
   *
   * cuDSS splits the communicator into a device and a host one. With the MPI
   * (openmpi) comm layer both backends are MPI, so both get the same MPI
   * communicator.  For an NCCL device backend, CUDSS_DATA_COMM_DEVICE would
   * instead take an ncclComm_t and only CUDSS_DATA_COMM_HOST would stay MPI.
   *
   * If you build against an older cuDSS that still has the single
   * CUDSS_DATA_COMM, replace the two calls below with one on CUDSS_DATA_COMM.
   * */
  if (!cudss_ok(cudssSetCommLayer(ctx->handle, ELMER_CUDSS_COMMLAYER_NAME), "cudssSetCommLayer") ||
      !cudss_ok(cudssDataSet(ctx->handle, ctx->data, CUDSS_DATA_COMM_DEVICE,
                    &ctx->comm, sizeof(&ctx->comm)), "cudssDataSet(COMM_DEVICE)") ||
      !cudss_ok(cudssDataSet(ctx->handle, ctx->data, CUDSS_DATA_COMM_HOST,
                    &ctx->comm, sizeof(&ctx->comm)), "cudssDataSet(COMM_HOST)")) {
    cudss_teardown(ctx);
    return NULL;
  }

  /* Global descriptor (nglob x nglob, nnz_glob) with this rank's LOCAL arrays;
   * SetDistributionRow1d tells cuDSS which global rows the local arrays are. For
   * a complex system nglob/nloc/nnz are complex counts and vt is CUDSS_C_64F. */
  if (!cudss_ok(cudssMatrixCreateCsr(&ctx->Amat, nglob, nglob, nnz_glob,
                    ctx->d_rows, NULL, ctx->d_cols, ctx->d_vals,
                    CUDSS_R_32I, CUDSS_R_32I, vt, mt, mv, CUDSS_BASE_ZERO),
                "cudssMatrixCreateCsr") ||
      !cudss_ok(cudssMatrixSetDistributionRow1d(ctx->Amat, first_row, last_row),
                "cudssMatrixSetDistributionRow1d(A)")) {
    cudss_teardown(ctx);
    return NULL;
  }

  /* b and x share the matrix's owned-block distribution (no overlap). ld is the
   * local leading dimension; for a single RHS the stride is immaterial. */
  if (!cudss_ok(cudssMatrixCreateDn(&ctx->bmat, nglob, 1, nloc > 0 ? nloc : 1,
                    ctx->d_b, vt, CUDSS_LAYOUT_COL_MAJOR), "cudssMatrixCreateDn(b)") ||
      !cudss_ok(cudssMatrixSetDistributionRow1d(ctx->bmat, first_row, last_row),
                "cudssMatrixSetDistributionRow1d(b)") ||
      !cudss_ok(cudssMatrixCreateDn(&ctx->xmat, nglob, 1, nloc > 0 ? nloc : 1,
                    ctx->d_x, vt, CUDSS_LAYOUT_COL_MAJOR), "cudssMatrixCreateDn(x)") ||
      !cudss_ok(cudssMatrixSetDistributionRow1d(ctx->xmat, first_row, last_row),
                "cudssMatrixSetDistributionRow1d(x)")) {
    cudss_teardown(ctx);
    return NULL;
  }

  /* Analysis + factorization are collective over the communicator. */
  if (!cudss_ok(cudssExecute(ctx->handle, CUDSS_PHASE_ANALYSIS, ctx->config, ctx->data,
                    ctx->Amat, ctx->xmat, ctx->bmat), "cudssExecute(ANALYSIS)") ||
      !cudss_ok(cudssExecute(ctx->handle, CUDSS_PHASE_FACTORIZATION, ctx->config, ctx->data,
                    ctx->Amat, ctx->xmat, ctx->bmat), "cudssExecute(FACTORIZATION)")) {
    cudss_teardown(ctx);
    return NULL;
  }

  return ctx;
}
#endif /* HAVE_MPI */

ElmerCUDSS *FC_FUNC_(cudss_dfactorize_par,CUDSS_DFACTORIZE_PAR)
    (int *nglob, int *nnz_glob, int *nloc, int *nnzloc, int *first_row,
     int *last_row, int *rows, int *cols, double *vals, int *mtype, int *fcomm)
{
#ifdef HAVE_MPI
  MPI_Comm comm = MPI_Comm_f2c(*fcomm);
  return cudss_dfactorize_par_impl(*nglob, *nnz_glob, *nloc, *nnzloc, *first_row,
      *last_row, rows, cols, vals, *mtype, 0, comm);
#else
  (void)nglob; (void)nnz_glob; (void)nloc; (void)nnzloc; (void)first_row;
  (void)last_row; (void)rows; (void)cols; (void)vals; (void)mtype; (void)fcomm;
  fprintf(stderr, "CUDSS_SolveSystem: distributed (MGMN) cuDSS needs MPI\n");
  return NULL;
#endif
}

/* Complex counterpart: the matrix is a complex CSR of order n (not Elmer's
 * real 2n), vals holds 2*nnzloc interleaved (re,im) doubles -- i.e. Fortran
 * COMPLEX(KIND=dp), layout-compatible with cuDoubleComplex. */
ElmerCUDSS *FC_FUNC_(cudss_zdfactorize_par,CUDSS_ZDFACTORIZE_PAR)
    (int *nglob, int *nnz_glob, int *nloc, int *nnzloc, int *first_row,
     int *last_row, int *rows, int *cols, double *vals, int *mtype, int *fcomm)
{
#ifdef HAVE_MPI
  MPI_Comm comm = MPI_Comm_f2c(*fcomm);
  return cudss_dfactorize_par_impl(*nglob, *nnz_glob, *nloc, *nnzloc, *first_row,
      *last_row, rows, cols, vals, *mtype, 1, comm);
#else
  (void)nglob; (void)nnz_glob; (void)nloc; (void)nnzloc; (void)first_row;
  (void)last_row; (void)rows; (void)cols; (void)vals; (void)mtype; (void)fcomm;
  fprintf(stderr, "CUDSS_SolveSystem: distributed (MGMN) cuDSS needs MPI\n");
  return NULL;
#endif
}

#endif /* HAVE_CUDSS */
