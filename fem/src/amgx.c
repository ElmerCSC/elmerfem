/* Copyright (c) 2011-2017, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "../config.h"

#ifdef HAVE_AMGX

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "amgx_c.h"

#include <mpi.h>

static int rank;

/* print callback (could be customized) */
void print_callback(const char *msg, int length)
  { if(rank==0) printf("%s", msg); }

typedef struct {
   AMGX_matrix_handle A;
   AMGX_config_handle cfg;
   AMGX_vector_handle b,x,u,v;
   AMGX_resources_handle rsrc;
   AMGX_solver_handle solver;
} ElmerAMGX;


#include "cuda_runtime.h"

/* CUDA error macro */
#define CUDA_SAFE_CALL(call) do {                                 \
  cudaError_t err = call;                                         \
  if(cudaSuccess != err) {                                        \
    fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", \
            __FILE__, __LINE__, cudaGetErrorString( err) );       \
    exit(EXIT_FAILURE);                                           \
  } } while (0)


void ElmerAMGXInitialize()
{
   static int first = 1;

   if ( first ) {
     first = 0;

     AMGX_SAFE_CALL(AMGX_initialize());
     AMGX_SAFE_CALL(AMGX_initialize_plugins());

     AMGX_SAFE_CALL(AMGX_register_print_callback(&print_callback));
     AMGX_SAFE_CALL(AMGX_install_signal_handler());
   }
}

void AMGXmv( int **a_in, int *n_in, int *rows, int *cols, double *vals,
           double *u_in, double *v_in,int *nonlin_update )
{
    int i,j,k,n = *n_in;

    AMGX_Mode mode;
    AMGX_SOLVE_STATUS status;

    ElmerAMGX *ptr;

    ElmerAMGXInitialize();

    ptr = (ElmerAMGX *)*a_in;
    if(!ptr)
    {
      ptr = (ElmerAMGX *)calloc(sizeof(ElmerAMGX),1);
      *a_in = (int *)ptr;

      mode = AMGX_mode_dDDI;

      AMGX_SAFE_CALL(AMGX_config_create(&ptr->cfg, "algorithm=AGGREGATION"));

      AMGX_resources_create_simple(&ptr->rsrc, ptr->cfg);
      AMGX_matrix_create(&ptr->A, ptr->rsrc, mode);
      AMGX_vector_create(&ptr->u, ptr->rsrc, mode);
      AMGX_vector_create(&ptr->v, ptr->rsrc, mode);

      AMGX_matrix_upload_all(ptr->A,n,rows[n],1,1,rows,cols,vals,NULL );
    }else  if(*nonlin_update) {
      AMGX_matrix_upload_all(ptr->A,n,rows[n],1,1,rows,cols,vals,NULL );
    }

    AMGX_vector_bind(ptr->u, ptr->A);
    AMGX_vector_upload( ptr->u, n, 1, u_in );
    AMGX_vector_bind(ptr->v, ptr->A);
    AMGX_vector_upload( ptr->v, n, 1, v_in );
    AMGX_matrix_vector_multiply(ptr->A, ptr->u, ptr->v);
    AMGX_vector_download(ptr->v, v_in);
}


/*
void AMGXSolve( int **a_in, int *n_in, int *rows, int *cols, double *vals,
  double *b_in, double *x_in,int *nonlin_update, char *config_name)
*/
void AMGXSolve( int **a_in, int *n_in, int *rows, int *cols, double *vals,
  double *b_in, double *x_in,int *nonlin_update, char *config_name, int *fcomm, 
     int *ng_in, int *part_vec, double *bnrm_in )
{
    int i,j,k,n = *n_in, ng=*ng_in, lrank, nranks, gpu_count;
    static MPI_Comm comm;
    double bnrm = *bnrm_in;

    AMGX_Mode mode;
    AMGX_SOLVE_STATUS status;

    ElmerAMGX *ptr;

    ElmerAMGXInitialize();

    ptr = (ElmerAMGX *)*a_in;
    if(!ptr)
    {
      ptr = (ElmerAMGX *)calloc(sizeof(ElmerAMGX),1);
      *a_in = (int *)ptr;

      mode = AMGX_mode_dDDI;

      AMGX_SAFE_CALL(AMGX_config_create_from_file(&ptr->cfg, config_name));

      if ( n==ng ) {
        nranks = 1; rank = 0;
        AMGX_resources_create_simple(&ptr->rsrc, ptr->cfg);
      } else {
        comm = MPI_Comm_f2c(*fcomm);
        MPI_Comm_size(comm, &nranks);
        MPI_Comm_rank(comm, &rank);

        CUDA_SAFE_CALL(cudaGetDeviceCount(&gpu_count));
        lrank = rank % gpu_count;
        CUDA_SAFE_CALL(cudaSetDevice(lrank));
        printf("Process %d selecting device %d\n", rank, lrank);

        AMGX_resources_create(&ptr->rsrc, ptr->cfg, &comm, 1, &lrank);
      }

      if ( rank==0 ) {
        fprintf( stderr, "----\n" );
        fprintf( stderr, "%s\n" , config_name );
        fprintf( stderr, "----\n" );
      }

      AMGX_matrix_create(&ptr->A, ptr->rsrc, mode);
      AMGX_vector_create(&ptr->x, ptr->rsrc, mode);
      AMGX_vector_create(&ptr->b, ptr->rsrc, mode);
      AMGX_solver_create(&ptr->solver, ptr->rsrc, mode, ptr->cfg);

      if ( n==ng ) {
        AMGX_matrix_upload_all(ptr->A,n,rows[n],1,1,rows,cols,vals,NULL );
      } else {
        AMGX_matrix_upload_all_global_32(ptr->A,ng,n,rows[n],1,1,rows,cols,vals,NULL,1,1,part_vec);
      }
    }else  if(*nonlin_update) {
      if ( n==ng ) {
        AMGX_matrix_upload_all(ptr->A,n,rows[n],1,1,rows,cols,vals,NULL );
      } else {
        AMGX_matrix_upload_all_global_32(ptr->A,ng,n,rows[n],1,1,rows,cols,vals,NULL,1,1,part_vec);
      }
    }

    if  ( n==ng )
    {
      double bnrm;

// scale by ||b|| to get comparable convergence criteria to other linear solvers.

      bnrm = 0.0;
      for(i<0; i<n; i++ ) bnrm += b_in[i]*b_in[i];
      bnrm = sqrt(bnrm);
      if ( bnrm < 1.e-16 ) bnrm = 1;


      for(i=0; i<n; i++ ) b_in[i] = b_in[i]/bnrm;
      AMGX_vector_upload( ptr->b, n, 1, b_in );
      for(i=0; i<n; i++ ) b_in[i] = b_in[i]*bnrm;

      for(i=0; i<n; i++ ) x_in[i] = x_in[i]/bnrm;
      AMGX_vector_upload( ptr->x, n, 1, x_in );

      // WRITE OUT MATRIX (REMOVE)
      //AMGX_write_system(ptr->A, ptr->b, ptr->x, "mymatrix.dat");
      //exit(0);
		  
      AMGX_solver_setup(ptr->solver, ptr->A);
      AMGX_solver_solve(ptr->solver, ptr->b, ptr->x);

      AMGX_vector_download(ptr->x, x_in);
      for(i=0; i<n; i++ ) x_in[i] = x_in[i]*bnrm;
    } else {
      double *b_t, *x_t;

      b_t = (double *)calloc(  ng,sizeof(double) );
      x_t = (double *)calloc(  ng,sizeof(double) );

      if ( bnrm < 1.e-16 ) bnrm = 1;

      for(i=0; i<n; i++ ) b_t[i] = b_in[i]/bnrm;
      AMGX_vector_bind(ptr->b, ptr->A);
      AMGX_vector_upload( ptr->b, n, 1, b_t );

      for(i=0; i<n; i++ ) x_t[i] = x_in[i]/bnrm;
      AMGX_vector_bind(ptr->x, ptr->A);
      AMGX_vector_upload( ptr->x, n, 1, x_t );

      MPI_Barrier(comm);
      AMGX_solver_setup(ptr->solver, ptr->A);
      MPI_Barrier(comm);
      AMGX_solver_solve(ptr->solver, ptr->b, ptr->x);

      AMGX_vector_download(ptr->x, x_t);
      for(i=0; i<n; i++ ) x_in[i] = x_t[i]*bnrm;

      free(b_t); free(x_t);
    }

//    AMGX_solver_get_status(ptr->solver, &status);

#if 0
    AMGX_solver_destroy(solver);
    AMGX_vector_destroy(x);
    AMGX_vector_destroy(b);

    AMGX_matrix_destroy(A);
    AMGX_resources_destroy(rsrc);

    /* destroy config (need to use AMGX_SAFE_CALL after this point) */
    AMGX_SAFE_CALL(AMGX_config_destroy(cfg));

    /* shutdown and exit */
    AMGX_SAFE_CALL(AMGX_finalize_plugins());
    AMGX_SAFE_CALL(AMGX_finalize());
#endif
}

#endif
