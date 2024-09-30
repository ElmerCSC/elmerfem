/* ************************************************************************
 * Copyright (C) 2018-2021 Advanced Micro Devices, Inc. All rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * ************************************************************************ */
#include "../config.h"

#ifdef HAVE_ROCALUTION
#include <cstring>
#include <mpi.h>
#include <map>
#include <set>
#include <rocalution/rocalution.hpp>

using namespace rocalution;

typedef struct {
  int *Lcols, *Lrows, init;
  double *Lvals;
} ElmerRocalution;
ElmerRocalution S;


static void my_irecv(int* buf, int count, int source, int tag, MPI_Comm comm, MPI_Request* request)
{
    MPI_Irecv(buf, count, MPI_INT, source, tag, comm, request);
}

static void
    my_irecv(int64_t* buf, int count, int source, int tag, MPI_Comm comm, MPI_Request* request)
{
    MPI_Irecv(buf, count, MPI_INT64_T, source, tag, comm, request);
}

static void
    my_isend(const int* buf, int count, int dest, int tag, MPI_Comm comm, MPI_Request* request)
{
    MPI_Isend(buf, count, MPI_INT, dest, tag, comm, request);
}

static void
    my_isend(const int64_t* buf, int count, int dest, int tag, MPI_Comm comm, MPI_Request* request)
{
    MPI_Isend(buf, count, MPI_INT64_T, dest, tag, comm, request);
}

template <typename ValueType>
void elmer_distribute_matrix(const MPI_Comm*    comm,
                       GlobalMatrix<ValueType>* gmat,
                       int *rows, int *cols, double *vals,
		       int ln, int gn, int *index_offset, ParallelManager* pm)
{
    int rank;
    int num_procs;

    MPI_Comm_rank(*comm, &rank);
    MPI_Comm_size(*comm, &num_procs);

    int64_t global_nrow = ln;
    int64_t global_ncol = ln;
    int64_t global_nnz  = rows[ln+1];

    int*   global_row_offset = rows;
    int*       global_col    = cols;
    ValueType* global_val    = vals;

    // Compute local matrix sizes
    std::vector<int> local_size(num_procs);

    for(int i = 0; i < num_procs; ++i)
    {
        local_size[i] = index_offset[i+1] - index_offset[i];
    }

    // Read sub matrix - row_offset
    int                  local_nrow = global_nrow;
    std::vector<int> local_row_offset(local_nrow + 1);

    for(int k = 0; k < local_nrow + 1; ++k)
    {
        local_row_offset[k] = global_row_offset[k];
    }

    // Read sub matrix - col and val
    int                    local_nnz = local_row_offset[local_nrow] - local_row_offset[0];
    std::vector<int>       local_col(local_nnz);
    std::vector<ValueType> local_val(local_nnz);

    for(int k = 0; k < local_nnz; ++k)
    {
        local_col[k] = global_col[k];
        local_val[k] = global_val[k];
    }

    // Shift row_offset entries

    int interior_nnz = 0;
    int ghost_nnz    = 0;
    int boundary_nnz = 0;
    int neighbors    = 0;

    std::vector<std::vector<int>> boundary(num_procs, std::vector<int>());
    std::vector<bool>                 neighbor(num_procs, false);
    std::vector<std::map<int, bool>>  checked(num_procs, std::map<int, bool>());

    for(int i = 0; i < local_nrow; ++i)
    {
        for(int j = local_row_offset[i]; j < local_row_offset[i + 1]; ++j)
        {

            // Interior point
            if(local_col[j] >= index_offset[rank] && local_col[j] < index_offset[rank + 1])
            {
                ++interior_nnz;
            }
            else
            {
                // Boundary point above current process
                // Loop over ranks above current process
                for(int r = num_procs-1; r >= 0; --r)
                {
                    if ( r==rank ) continue;

                    // Check if boundary belongs to rank r
                    if(local_col[j] >= index_offset[r] && local_col[j] < index_offset[r + 1])
                    {
                        // Add boundary point to rank r if it has not been added yet
                        if(!checked[r][i + index_offset[rank]])
                        {
                            boundary[r].push_back(i + index_offset[rank]);
                            neighbor[r] = true;
                            ++boundary_nnz;
                            checked[r][i + index_offset[rank]] = true;
                        }
                        ++ghost_nnz;
                        // Rank for current boundary point local_col[j] has been found
                        // Continue with next boundary point
                        break;
                    }
                }
            }
        }
    }

    for(int i = 0; i < num_procs; ++i)
    {
       if(neighbor[i] == true) ++neighbors;
    }

    std::vector<MPI_Request> mpi_req(neighbors * 4);
    int                      n = 0;
    // Array to hold boundary size for each interface
    std::vector<int> boundary_size(neighbors);

    // MPI receive boundary sizes
    for(int i = 0; i < num_procs; ++i)
    {
        // If neighbor receive from rank i is expected...
        if(neighbor[i] == true)
        {
            // Receive size of boundary from rank i to current rank
            my_irecv(&(boundary_size[n]), 1, i, 0, *comm, &mpi_req[n]);
            ++n;
        }
    }


    // MPI send boundary sizes
    int size[num_procs];
    for(int i = 0; i < num_procs; ++i)
    {
        // Send required if boundary for rank i available
        if(boundary[i].size() > 0)
        {
            size[i] = boundary[i].size();
            // Send size of boundary from current rank to rank i

            my_isend(&size[i], 1, i, 0, *comm, &mpi_req[n]);
            ++n;
        }
    }
    // Wait to finish communication
    MPI_Waitall(n - 1, &(mpi_req[0]), MPI_STATUSES_IGNORE);

    n = 0;
    // Array to hold boundary offset for each interface
    int              k = 0;
    std::vector<int> recv_offset(neighbors + 1);
    std::vector<int> send_offset(neighbors + 1);
    recv_offset[0] = 0;
    send_offset[0] = 0;
    for(int i = 0; i < neighbors; ++i)
    {
        recv_offset[i + 1] = recv_offset[i] + boundary_size[i];
    }

    for(int i = 0; i < num_procs; ++i)
    {
        if(neighbor[i] == true)
        {
            send_offset[k + 1] = send_offset[k] + boundary[i].size();
            ++k;
        }
    }

    // Array to hold boundary for each interface
    std::vector<std::vector<int>> local_boundary(neighbors);
    for(int i = 0; i < neighbors; ++i)
    {
        local_boundary[i].resize(boundary_size[i]);
    }

    // MPI receive boundary
    for(int i = 0; i < num_procs; ++i)
    {
        // If neighbor receive from rank i is expected...
        if(neighbor[i] == true)
        {
            // Receive boundary from rank i to current rank
            my_irecv(local_boundary[n].data(), boundary_size[n], i, 0, *comm, &mpi_req[n]);
            ++n;
        }
    }

    // MPI send boundary
    for(int i = 0; i < num_procs; ++i)
    {
        // Send required if boundary for rank i is available
        if(boundary[i].size() > 0)
        {
            // Send boundary from current rank to rank i
            my_isend(&(boundary[i][0]), boundary[i].size(), i, 0, *comm, &mpi_req[n]);
            ++n;
        }
    }

    // Wait to finish communication
    MPI_Waitall(n - 1, &(mpi_req[0]), MPI_STATUSES_IGNORE);

    // Total boundary size
    int nnz_boundary = 0;
    for(int i = 0; i < neighbors; ++i)
    {
        nnz_boundary += boundary_size[i];
    }

    // Create local boundary index array
    k = 0;
    std::vector<int> bnd(boundary_nnz);

    for(int i = 0; i < num_procs; ++i)
    {
        for(unsigned int j = 0; j < boundary[i].size(); ++j)
        {
            bnd[k] = static_cast<int>(boundary[i][j] - index_offset[rank]);
            ++k;
        }
    }

    // Create boundary index array
    std::vector<int> boundary_index(nnz_boundary);

    k = 0;
    for(int i = 0; i < neighbors; ++i)
    {
        for(int j = 0; j < boundary_size[i]; ++j)
        {
            boundary_index[k] = local_boundary[i][j];
            ++k;
        }
    }


    // Create map with boundary index relations
    std::map<int, int> boundary_map;

    for(int i = 0; i < nnz_boundary; ++i)
    {
        boundary_map[boundary_index[i]] = i;
    }

    // Build up ghost and interior matrix
    int*       ghost_row = new int[ghost_nnz];
    int*       ghost_col = new int[ghost_nnz];
    ValueType* ghost_val = new ValueType[ghost_nnz];

    memset(ghost_row, 0, sizeof(int) * ghost_nnz);
    memset(ghost_col, 0, sizeof(int) * ghost_nnz);
    memset(ghost_val, 0, sizeof(ValueType) * ghost_nnz);

    int*   row_offset = new int[local_nrow + 1];
    int*       col        = new int[interior_nnz];
    ValueType* val        = new ValueType[interior_nnz];

    memset(row_offset, 0, sizeof(int) * (local_nrow + 1));
    memset(col, 0, sizeof(int) * interior_nnz);
    memset(val, 0, sizeof(ValueType) * interior_nnz);

    row_offset[0] = 0;
    k             = 0;
    int l         = 0;
    for(int i = 0; i < local_nrow; ++i)
    {
        for(int j = local_row_offset[i]; j < local_row_offset[i + 1]; ++j)
        {

            // Boundary point -- create ghost part
            if(local_col[j] < index_offset[rank] || local_col[j] >= index_offset[rank + 1])
            {
                ghost_row[k] = i;
                ghost_col[k] = boundary_map[local_col[j]];
                ghost_val[k] = local_val[j];
                ++k;
            }
            else
            {
                // Interior point -- create interior part
                int c = local_col[j] - index_offset[rank];

                col[l] = c;
                val[l] = local_val[j];
                ++l;
            }
        }
        row_offset[i + 1] = l;
    }

    std::vector<int> recv(neighbors);
    std::vector<int> sender(neighbors);

    int nbc = 0;
    for(int i = 0; i < num_procs; ++i)
    {
        if(neighbor[i] == true)
        {
            recv[nbc]   = i;
            sender[nbc] = i;
            ++nbc;
        }
    }

    pm->SetMPICommunicator(comm);
    pm->SetGlobalNrow(gn);
    pm->SetGlobalNcol(gn);
    pm->SetLocalNrow(global_nrow);
    pm->SetLocalNcol(global_nrow);
    pm->SetBoundaryIndex(boundary_nnz, bnd.data());
    pm->SetReceivers(neighbors, recv.data(), recv_offset.data());
    pm->SetSenders(neighbors, sender.data(), send_offset.data());

    gmat->SetParallelManager(*pm);
    gmat->SetLocalDataPtrCSR(&row_offset, &col, &val, "mat", interior_nnz);
    gmat->SetGhostDataPtrCOO(&ghost_row, &ghost_col, &ghost_val, "ghost", ghost_nnz);
    gmat->Sort();
}


extern "C" void ROCParallelSolve( int *gn, int *n, int *rows, int *cols, double *vals, double *b, double *x_out, 
   double *bnrm, int *gOffset,int *fcomm, int *imethod, int *prec, int *maxiter, double *TOL )
{
    int i, *Lrows, *Lcols, rank, nranks;
    double *Lvals;
    MPI_Comm comm;

    rank = 0; nranks = 1;
    if (*fcomm) {
      comm = MPI_Comm_f2c(*fcomm);
      MPI_Comm_size(comm, &nranks);
      MPI_Comm_rank(comm, &rank);
    }

        // Disable OpenMP thread affinity
//    set_omp_affinity_rocalution(false);

    // Initialize platform with rank and # of accelerator devices in the node
    init_rocalution(rank);

    // Disable OpenMP
    set_omp_threads_rocalution(1);

    // Start time measurement
    double tick, tack;
    tick = rocalution_time();

    ParallelManager manager;
    manager.SetMPICommunicator(&comm);

    GlobalMatrix<double> gmat;

    elmer_distribute_matrix(&comm, &gmat, rows, cols, vals, *n, *gn, gOffset, &manager);

    GlobalVector<double> x(manager);
    GlobalVector<double> rhs(manager);


    // Allocate vectors
    x.Allocate("x", gmat.GetN());
    rhs.Allocate("rhs", gmat.GetM());

    for(i=0; i<*n; i++ ) x[i]=x_out[i];
    for(i=0; i<*n; i++ ) rhs[i]=b[i];

    // Move objects to accelerator
    gmat.MoveToAccelerator();
    x.MoveToAccelerator();
    rhs.MoveToAccelerator();

    // Linear Solver
    IterativeLinearSolver<GlobalMatrix<double>, GlobalVector<double>, double> *ls;

    CG<GlobalMatrix<double>, GlobalVector<double>, double> ls_cg;
    BiCGStab<GlobalMatrix<double>, GlobalVector<double>, double> ls_bcg;
    BiCGStabl<GlobalMatrix<double>, GlobalVector<double>, double> ls_bcgl;
    GMRES<GlobalMatrix<double>, GlobalVector<double>, double> ls_gmres;
    FGMRES<GlobalMatrix<double>, GlobalVector<double>, double> ls_fgmres;

    switch(*imethod) {
      case(0): ls = &ls_cg;  break;
      case(1): ls = &ls_bcg; break;
      case(2): ls = &ls_bcgl; ls_bcgl.SetOrder(4); break;
      case(3): ls = &ls_gmres; ls_gmres.SetBasisSize(50); break;
      case(4): ls = &ls_fgmres; ls_fgmres.SetBasisSize(50); break;
      default: ls = &ls_bcg; break;
    }

    ls->Init(*TOL*(*bnrm),1e-20,1e20,*maxiter);

    // Preconditioner
    Jacobi<GlobalMatrix<double>, GlobalVector<double>, double> prec_j;
//  SGS<GlobalMatrix<double>, GlobalVector<double>, double> prec_g;

    switch(*prec) {
      case(0): ls->SetPreconditioner(prec_j); break;
//    case(1): ls->SetPreconditioner(prec_g); break;
    }

    ls->SetOperator(gmat);
    ls->Build();
    ls->Verbose(2);

    gmat.Info();

    ls->Solve(rhs, &x);

//  x.CopyToData(x_out);

    for(i=0; i<*n; i++ ) x_out[i]=x[i];

    ls->Clear();
}

extern "C" void ROCSerialSolve(int *n, int *rows, int *cols, double *vals, double *b, double *x_out,
      int *nonlin_update, int *imethod, int *prec, int *maxiter, double *TOL)
{

    int i, *Lrows, *Lcols, rank, nranks;
    double *Lvals;


    // Initialize rocALUTION
    init_rocalution();

    // Start time measurement
    double tick, tack;
    tick = rocalution_time();


    LocalMatrix<double> mat;
    LocalVector<double> x, rhs;

    Lrows = (int *)malloc((*n+1)*sizeof(int));
    Lcols = (int *)malloc((rows[*n])*sizeof(int));
    Lvals= (double *)malloc((rows[*n])*sizeof(double));

    for(i=0; i<=*n; i++ ) Lrows[i] = rows[i];
    for(i=0; i<rows[*n]; i++ ) Lcols[i] = cols[i];
    for(i=0; i<rows[*n]; i++ ) Lvals[i] = vals[i];

    mat.SetDataPtrCSR(&Lrows, &Lcols, &Lvals, "A", Lrows[*n], *n, *n);

    // Move objects to accelerator
    mat.MoveToAccelerator();
    x.MoveToAccelerator();
    rhs.MoveToAccelerator();

    // Allocate vectors
    x.Allocate("x", mat.GetN());
    rhs.Allocate("rhs", mat.GetM());

    rhs.CopyFromData(b);
    x.CopyFromData(x_out);

    // Linear Solver
    IterativeLinearSolver<LocalMatrix<double>, LocalVector<double>, double> *ls;

    CG<LocalMatrix<double>, LocalVector<double>, double> ls_cg;
    BiCGStab<LocalMatrix<double>, LocalVector<double>, double> ls_bcg;
    BiCGStabl<LocalMatrix<double>, LocalVector<double>, double> ls_bcgl;
    GMRES<LocalMatrix<double>, LocalVector<double>, double> ls_gmres;
    FGMRES<LocalMatrix<double>, LocalVector<double>, double> ls_fgmres;

    switch(*imethod) {
      case(0): ls = &ls_cg;  break;
      case(1): ls = &ls_bcg; break;
      case(2): ls = &ls_bcgl; ls_bcgl.SetOrder(4); break;
      case(3): ls = &ls_gmres; ls_gmres.SetBasisSize(50); break;
      case(4): ls = &ls_fgmres; ls_fgmres.SetBasisSize(50); break;
      default: ls = &ls_bcg; break;
    }

    {
      double bnrm;
      bnrm = 0.0;
      for(i=0; i<*n; i++ ) bnrm += b[i]*b[i];
      bnrm = sqrt(bnrm);
      if (bnrm<1e-16) bnrm = 1;
      ls->Init(*TOL*bnrm,1e-20,1e20,*maxiter);
    }

    // Preconditioner
    Jacobi<LocalMatrix<double>, LocalVector<double>, double> prec_j;
    SGS<LocalMatrix<double>, LocalVector<double>, double> prec_g;
    ILU<LocalMatrix<double>, LocalVector<double>, double> prec_i;

    switch(*prec) {
      case(0): ls->SetPreconditioner(prec_j); break;
      case(1): ls->SetPreconditioner(prec_g); break;
      case(2):
        ls->SetPreconditioner(prec_i);
        {
          int level=0;
          prec_i.Set(level);
        }
        break;
    }

    ls->SetOperator(mat);
    ls->Build();
    ls->Verbose(2);
    mat.Info();
    ls->Solve(rhs, &x);
    x.CopyToData(x_out);

    ls->Clear();

    // Stop time measurement
    tack = rocalution_time();
    std::cout << "Solver execution:" << (tack - tick) / 1e6 << " sec" << std::endl;

    // Stop rocALUTION platform
    stop_rocalution();

    free(Lrows); free(Lcols); free(Lvals);
} 
#else
#include <stdio.h>
#include <stdlib.h>
extern "C" void ROCSerialSolve(int *n, int *rows, int *cols, double *vals, double *b, double *x_out, int *nonlin_update)
{
  fprintf( stderr, "No serial ROCALUTION library included.\n" );
  exit(0);
}

extern "C" void ROCParallelSolve( int *gn, int *n, int *rows, int *cols, double *vals, double *b, double *x_out, 
		          int *gOffset,int *fcomm )
{
  fprintf( stderr, "No parallel ROCALUTION library included.\n" );
  exit(0);
}
#endif


