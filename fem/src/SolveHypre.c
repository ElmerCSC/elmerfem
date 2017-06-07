/*****************************************************************************
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! ******************************************************************************
! *
! *  Elmer interface for Hypre - High performance Preconditioners
! *
! *  For more information on Hypre see
! *  https://computation.llnl.gov/casc/linear_solvers/sls_hypre.html
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Thomas Zwinger, Jonas Thies, Peter Råback, Mika Malinen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 2000
! *
! *****************************************************************************/


#include "../config.h"

#ifdef HAVE_HYPRE
#include <math.h>
#include "_hypre_utilities.h"
#include "krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

#ifdef USE_ISO_C_BINDINGS
#define realtime_ FC_FUNC_(realtime,REALTIME)
#endif
typedef struct {

int ilower, iupper;

HYPRE_IJMatrix A;
HYPRE_IJMatrix Atilde;

int hypre_method;
HYPRE_Solver solver, precond;

} ElmerHypreContainer;

/* If the version of HYPRE is new enough, FlexGMRES and LGMRES solvers can be included
   by turning the following flag on */
#define HAVE_GMRES 1

/* there are two possible procedures of calling HYPRE here, 
  the standard one (does everything once), and a step-wise
  procedure of setup, solve and cleanup.
  The first one is obsolite. 
  TO DO: we should add the possibility to keep the precon-
  ditioner the same but update the system matrix (SolveHYPRE3), right now
  calling SolveHYPRE2 solves with the matrix passed into   
  SolveHYPRE1.

 standard call: - convert matrix
                - convert vector b
                - setup solver and preconditioner
                - solve system
                - convert vector x back
                - destroy all data structures
*/
void STDCALLBULL FC_FUNC(solvehypre,SOLVEHYPRE)
 (
  int *nrows,int *rows, int *cols, double *vals, int *perm,
  int *invperm, int *globaldofs, int *owner,  double *xvec,
  double *rhsvec, int *pe, int *ILUn, int *Rounds, double *TOL,
  int *hypre_method, int *hypre_intpara, double *hypre_dppara,int *fcomm
 )
{
   int i, j, k, *rcols;
   int myid, num_procs;
   int N, n;

   int ilower, iupper;
   int local_size, extra;

   int solver_id;
   int print_solution, print_system;


   double  *txvec, st, realtime_();

   HYPRE_IJMatrix A;
   HYPRE_ParCSRMatrix parcsr_A;
   HYPRE_IJVector b;
   HYPRE_ParVector par_b;
   HYPRE_IJVector x;
   HYPRE_ParVector par_x;

   HYPRE_Solver solver, precond;
   int verbosity = 10;

   MPI_Comm comm=MPI_Comm_f2c(*fcomm);
   
   st  = realtime_();
   /* How many rows do I have? */
   local_size = *nrows;

   ilower=1000000000;
   iupper=0;
   for( i=0; i<local_size; i++ )
   {
      if ( owner[i] ) {
        if ( iupper < globaldofs[i] ) iupper=globaldofs[i];
        if ( ilower > globaldofs[i] ) ilower=globaldofs[i];
      }
   }

  /* which process number am I? */
   MPI_Comm_rank(comm, &myid);
   /* Create the matrix.
      Note that this is a square matrix, so we indicate the row partition
      size twice (since number of rows = number of cols) */
   HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &A);

   /* Choose a parallel csr format storage (see the User's Manual) */
   HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

   /* Initialize before setting coefficients */
   HYPRE_IJMatrixInitialize(A);

   /* Now go through my local rows and set the matrix entries.
      Note that here we are setting one row at a time, though
      one could set all the rows together (see the User's Manual).
   */
   {
      int nnz,irow,i,j,k,*rcols,csize=32;

      rcols = (int *)malloc( csize*sizeof(int) );
      for (i = 0; i < local_size; i++)
      {
         nnz = rows[i+1]-rows[i];
         if ( nnz>csize ) {
           rcols = (int *)realloc( rcols, nnz*sizeof(int) );
           csize = nnz;
         }
         irow=globaldofs[i];
         for( k=0,j=rows[i]; j<rows[i+1]; j++,k++)
         {
           rcols[k] = globaldofs[cols[j-1]-1];
         }
         HYPRE_IJMatrixAddToValues(A, 1, &nnz, &irow, rcols, &vals[rows[i]-1]);
      }
        free( rcols );
   }

   /* Assemble after setting the coefficients */
   HYPRE_IJMatrixAssemble(A);

   /* Get the parcsr matrix object to use */
   HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);

   /* Create the rhs and solution */
   rcols = (int *)malloc( local_size*sizeof(int) );
   txvec = (double *)malloc( local_size*sizeof(double) );
   for( k=0,i=0; i<local_size; i++ ) rcols[k++] = globaldofs[i];

   HYPRE_IJVectorCreate(comm, ilower, iupper,&b);
   HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(b);
   for( i=0; i<local_size; i++ ) txvec[i] = rhsvec[i];
   HYPRE_IJVectorAddToValues(b, local_size, rcols, txvec );

   HYPRE_IJVectorCreate(comm, ilower, iupper,&x);
   HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(x);
   for( i=0; i<local_size; i++ ) txvec[i] = xvec[i];
   HYPRE_IJVectorSetValues(x, local_size, rcols, txvec );

   HYPRE_IJVectorAssemble(b);
   HYPRE_IJVectorGetObject(b, (void **) &par_b);

   HYPRE_IJVectorAssemble(x);
   HYPRE_IJVectorGetObject(x, (void **) &par_x);
   if( verbosity >= 12 ) fprintf( stderr, "ID no. %i: setup time: %g\n", myid, realtime_()-st );
   st = realtime_();


/*    fprintf(stderr,"HYRPE INT: %d %d  %d %d %d \n", hypre_intpara[0], hypre_intpara[1], hypre_intpara[2], hypre_intpara[3], hypre_intpara[4]);  */
/*    fprintf(stderr,"HYRPE DP: %d %d %d %d %d \n", hypre_dppara[0], hypre_dppara[1], hypre_dppara[2], hypre_dppara[3], hypre_dppara[4]);  */
   /* Choose a solver and solve the system */
   /* NB.: hypremethod = 0 ... BiCGStab + ILUn
                         1 ... BiCGStab + ParaSails
                         2 ... BiCGStab + BoomerAMG
                        10 ... BoomerAMG 
   */
   if ( *hypre_method < 10) { /* BiGSTAB methods */
     /* Create solver */
     HYPRE_ParCSRBiCGSTABCreate(comm, &solver);

     /* Set some parameters (See Reference Manual for more parameters) */
     HYPRE_ParCSRBiCGSTABSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_ParCSRBiCGSTABSetTol(solver, *TOL);       /* conv. tolerance */
     HYPRE_ParCSRBiCGSTABSetStopCrit(solver, 0);     /* use the two norm as the stopping criteria */
     HYPRE_ParCSRBiCGSTABSetPrintLevel(solver, 2);   /* print solve info */
     HYPRE_ParCSRBiCGSTABSetLogging(solver, 1);      /* needed to get run info later */
     if ( *hypre_method == 0 ) {
       HYPRE_EuclidCreate( comm, &precond );
       {
         static char *argv[5], str[3];
         argv[0] = "-level";
         sprintf( str, "%d", *ILUn );
	 if (myid == 0 & verbosity >= 5) fprintf( stderr,"SolveHypre: using BiCGStab + ILU%i\n",*ILUn); 
         argv[1] = str;
         HYPRE_EuclidSetParams( precond, 2, argv );
       }

       /* Set the PCG preconditioner */
       HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup, precond);
     } else if (*hypre_method == 1) { 
       if (myid == 0 & verbosity >= 5) fprintf( stderr,"SolveHypre: using BiCGStab + paraSails\n"); 

       /* Now set up the ParaSails preconditioner and specify any parameters */
       HYPRE_ParaSailsCreate(comm, &precond);
       {
	 /* Set some parameters (See Reference Manual for more parameters) */
         /* threshold = dppara[0]; maxlevels= intpara[1] */
	 HYPRE_ParaSailsSetParams(precond, hypre_dppara[0], hypre_intpara[1]);
	 /* filter = dppara[1] */
	 HYPRE_ParaSailsSetFilter(precond, hypre_dppara[1]);
         /* symmetry = intpara[0] */
	 HYPRE_ParaSailsSetSym(precond, hypre_intpara[0]);
	 HYPRE_ParaSailsSetLogging(precond, 3);
       }
       /* Set the PCG preconditioner */
       HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);

     } else if(*hypre_method == 2){
       if (myid == 0 ) {
	 if( verbosity >= 5 ) {
	   fprintf( stderr,"SolveHypre: using BiCGStab + boomerAMG\n");
	 }
	 if( verbosity >= 10 ) {
	   fprintf( stderr,"RelaxType = %d\n",hypre_intpara[0]); 
	   fprintf( stderr,"CoarsenType = %d\n",hypre_intpara[1]); 
	   fprintf( stderr,"NumSweeps = %d\n",hypre_intpara[2]); 
	   fprintf( stderr,"MaxLevels = %d\n",hypre_intpara[3]); 
	   fprintf( stderr,"Interpolation Type = %d\n",hypre_intpara[4]); 
	   fprintf( stderr,"Smooth Type = %d\n",hypre_intpara[5]);
	   fprintf( stderr,"Cycle Type = %d\n",hypre_intpara[6]);
	   fprintf( stderr,"DOFs = %d\n",hypre_intpara[7]);
	 }
       }
       HYPRE_BoomerAMGCreate(&precond); 

       /* Set some parameters (See Reference Manual for more parameters) */
       HYPRE_BoomerAMGSetNumFunctions(precond, hypre_intpara[7]); /* No. of PDE's */
       HYPRE_BoomerAMGSetPrintLevel(precond, 1); /* print amg solution info */
       HYPRE_BoomerAMGSetNumSweeps(precond, 1); /* fixed for preconditioner to 1 */
       HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
       HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */
       HYPRE_BoomerAMGSetRelaxType(precond, hypre_intpara[0]);   /* G-S/Jacobi hybrid relaxation */
       HYPRE_BoomerAMGSetCoarsenType(precond, hypre_intpara[1]);  /* coarsening type */
       
       HYPRE_BoomerAMGSetMaxLevels(precond, hypre_intpara[3]); /* levels of coarsening */
       HYPRE_BoomerAMGSetInterpType(precond, hypre_intpara[4]);  /* interpolation type */
       HYPRE_BoomerAMGSetSmoothType(precond, hypre_intpara[5]);  /* smoother type */

       /* Set the BiCGSTAB preconditioner */
       HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
     } else if (*hypre_method != 9) {
       fprintf( stderr,"Hypre preconditioning method not implemented\n");
       exit(EXIT_FAILURE);
     }
     

     /* Now setup and solve! */
     HYPRE_ParCSRBiCGSTABSetup(solver, parcsr_A, par_b, par_x);
     HYPRE_ParCSRBiCGSTABSolve(solver, parcsr_A, par_b, par_x);

     /* Destroy solver and preconditioner */
     HYPRE_ParCSRBiCGSTABDestroy(solver);
     if ( *hypre_method == 0 ) {
       HYPRE_EuclidDestroy(precond);
     } else if ( *hypre_method == 1 ) {
       HYPRE_ParaSailsDestroy(precond);
     } else {
       HYPRE_BoomerAMGDestroy(precond);
     }
   } else if ( *hypre_method == 10 ) { /* boomer AMG */
      int num_iterations;
      double final_res_norm;

      if (myid == 0) {
	if( verbosity >= 5 ) {
	  fprintf( stderr,"SolveHypre: using BoomerAMG\n"); 
	}
	if( verbosity >= 10 ) {
	  fprintf( stderr,"RelaxType = %d\n",hypre_intpara[0]); 
	  fprintf( stderr,"CoarsenType = %d\n",hypre_intpara[1]); 
	  fprintf( stderr,"NumSweeps = %d\n",hypre_intpara[2]); 
	  fprintf( stderr,"MaxLevels = %d\n",hypre_intpara[3]); 
	  fprintf( stderr,"Interpolation Type = %d\n",hypre_intpara[4]); 
	  fprintf( stderr,"Smooth Type = %d\n",hypre_intpara[5]);
	  fprintf( stderr,"Cycle Type = %d\n",hypre_intpara[6]);
	  fprintf( stderr,"DOFs = %d\n",hypre_intpara[7]);
	}
      }
      /* Create solver */
      HYPRE_BoomerAMGCreate(&solver);

      /* Set some parameters (See Reference Manual for more parameters) */
      HYPRE_BoomerAMGSetNumFunctions(solver, hypre_intpara[7]); /* No. of PDE's */
      HYPRE_BoomerAMGSetPrintLevel(solver, 3);  
      HYPRE_BoomerAMGSetRelaxType(solver, hypre_intpara[0]);   /* G-S/Jacobi hybrid relaxation */
      HYPRE_BoomerAMGSetCoarsenType(solver, hypre_intpara[1]);  /* coarsening type */
      HYPRE_BoomerAMGSetNumSweeps(solver, hypre_intpara[2]);   /* Sweeeps on each level */
      HYPRE_BoomerAMGSetMaxLevels(solver, hypre_intpara[3]); /* levels of coarsening */
      HYPRE_BoomerAMGSetInterpType(solver, hypre_intpara[4]);  /* interpolation type */
      HYPRE_BoomerAMGSetSmoothType(solver, hypre_intpara[5]);  /* smoother type */
      HYPRE_BoomerAMGSetCycleType(precond, hypre_intpara[6]);  /* coarsening type */
      HYPRE_BoomerAMGSetTol(solver, *TOL);      /* conv. tolerance */
      HYPRE_BoomerAMGSetMaxIter(solver, *Rounds); /* iteration rounds */

      /* Now setup and solve! */
      HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
      HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);

      /* Run info - needed logging turned on */
      HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
      HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
      if (myid == 0) {
	if( verbosity >= 5 ) {
	  fprintf(stderr,"BoomerAMG:\n");
	  fprintf(stderr,"Iterations = %d\n", num_iterations);
	  fprintf(stderr,"Final Relative Residual Norm = %e\n", final_res_norm);
	  fprintf(stderr,"\n");
	}
      }

      /* Destroy solver */
      HYPRE_BoomerAMGDestroy(solver);
   } else {
     fprintf( stderr,"Hypre solver not implemented\n");
     exit(EXIT_FAILURE);
   }

   for( k=0,i=0; i<local_size; i++ )
     if ( owner[i] ) rcols[k++] = globaldofs[i];
   
   HYPRE_IJVectorGetValues(x, k, rcols, txvec );
   
   for( i=0,k=0; i<local_size; i++ )
     if ( owner[i] ) xvec[i] = txvec[k++];
   
   if( myid == 0 && verbosity >= 5 ) {
     fprintf( stderr, "Hypre solve time: %g\n", realtime_()-st );
   }
   free( txvec );
   free( rcols );

   /* Clean up */
   HYPRE_IJMatrixDestroy(A);
   HYPRE_IJVectorDestroy(b);
   HYPRE_IJVectorDestroy(x);
}

/*///////////////////////////////////////////////////////////////////////////////////////////////*/

/* initialization for a new matrix.
      - convert matrix
      - setup solver and preconditioner
      - return a pointer 'Container' which the calling fortran
        program should not alter but pass back into subsequent
        SolveHYPRE2, ~3 and ~4 calls.

 This function has an additional feature compared to the SolveHYPRE call above,
 namely to use a block diagonal approximation of A for the preconditioner setup.
 This mimics the behavior of the BILUn preconditioners in Elmer, although any   
 preconditioner (like ParaSails or BoomerAMG) can still be used in combination  
 with block diagonal approximation. 
 BILU=0 or 1 - use A. 
 BILU=k - assume k equations and use block diagonal A with k blocks.
*/
void STDCALLBULL FC_FUNC(solvehypre1,SOLVEHYPRE1)
 (
  int *nrows,int *rows, int *cols, double *vals, int *precflag, double *precvals, 
    int *globaldofs, int *owner, int *ILUn, int *BILU, int *hypre_method,
      int *hypre_intpara, double *hypre_dppara,
         int *Rounds, double *TOL, int *verbosityPtr, int** ContainerPtr,
  int *fcomm
 )
{
   int i, j, k, *rcols;
   int myid, num_procs;
   int N, n, csize=128;

   int ilower, iupper;
   int local_size, extra;
   int hypre_sol, hypre_pre;
   MPI_Comm comm=MPI_Comm_f2c(*fcomm);
   ElmerHypreContainer* Container;

   HYPRE_IJMatrix A, Atilde;
   HYPRE_ParCSRMatrix parcsr_A;
   HYPRE_IJVector b;
   HYPRE_ParVector par_b;
   HYPRE_IJVector x;
   HYPRE_ParVector par_x;

   HYPRE_Solver solver, precond;

   double  *txvec, st, realtime_();
   
   int verbosity = *verbosityPtr;

   /* which process number am I? */
   MPI_Comm_rank(comm, &myid);
   
   if (myid==0 && verbosity >= 4) fprintf(stdout,"Performing HYPRE Setup\n");

   if (*ContainerPtr != 0) {
     fprintf( stderr, "ID no. %i: pointer passed into SolveHypre1 not NULL, possible memory leak.\n", myid);
   }
   
   Container = (ElmerHypreContainer*)malloc(sizeof(ElmerHypreContainer));
   
   *ContainerPtr=(int*)(Container);
   
   st  = realtime_();
   
   /* How many rows do I have? */
   local_size = *nrows;
   hypre_sol = *hypre_method / 10;
   hypre_pre = *hypre_method % 10;

   /* No preconditioner for BoomerAMG */
   if( hypre_sol == 1 ) hypre_pre = -1;

   ilower=1000000000;
   iupper=0;
   for( i=0; i<local_size; i++ ) {
     if ( owner[i] ) {
       if ( iupper < globaldofs[i] ) iupper=globaldofs[i];
       if ( ilower > globaldofs[i] ) ilower=globaldofs[i];
     }
   }
   
   /* Create the matrix.
      Note that this is a square matrix, so we indicate the row partition
      size twice (since number of rows = number of cols) */
   HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &A);

   /* Choose a parallel csr format storage (see the User's Manual) */
   HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

   /* Initialize before setting coefficients */
   HYPRE_IJMatrixInitialize(A);

   /* Now go through my local rows and set the matrix entries.
      Note that here we are setting one row at a time, though
      one could set all the rows together (see the User's Manual).
   */
   {
      int nnz,irow,i,j,k,*rcols;

      rcols = (int *)malloc( csize*sizeof(int) );
      for (i = 0; i < local_size; i++) {
	nnz = rows[i+1]-rows[i];
	if ( nnz>csize ) {
	  csize = nnz+csize;
	  rcols = (int *)realloc( rcols, csize*sizeof(int) );
	}
	irow=globaldofs[i];
	for( k=0,j=rows[i]; j<rows[i+1]; j++,k++) {
	  rcols[k] = globaldofs[cols[j-1]-1];
	}
	HYPRE_IJMatrixAddToValues(A, 1, &nnz, &irow, rcols, &vals[rows[i]-1]);
      }
      free( rcols );
   }

   /* Assemble after setting the coefficients */
   HYPRE_IJMatrixAssemble(A);

   if (!*precflag && *BILU <= 1) {
     Atilde = A;
   } else if ( *precflag ) {
     int nnz,irow,jcol,i,j,k,*rcols;
     double *dbuf;

     HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &Atilde);
     HYPRE_IJMatrixSetObjectType(Atilde, HYPRE_PARCSR);
     HYPRE_IJMatrixInitialize(Atilde);
     {
        int nnz,irow,i,j,k,*rcols;

        rcols = (int *)malloc( csize*sizeof(int) );
        for (i = 0; i < local_size; i++) {
          nnz = rows[i+1]-rows[i];
          if ( nnz>csize ) {
            csize = nnz+csize;
            rcols = (int *)realloc( rcols, csize*sizeof(int) );
          }
          irow=globaldofs[i];
          for( k=0,j=rows[i]; j<rows[i+1]; j++,k++) {
             rcols[k] = globaldofs[cols[j-1]-1];
          }
          HYPRE_IJMatrixAddToValues(Atilde, 1, &nnz, &irow, rcols, &precvals[rows[i]-1]);
        }
        free( rcols );
     }
     /* Assemble after setting the coefficients */
     HYPRE_IJMatrixAssemble(Atilde);     
   } else {
     int nnz,irow,jcol,i,j,k,*rcols;
     double *dbuf;
     if (myid==0 && verbosity >= 5) fprintf(stdout,"HYPRE: using BILU(%d) approximation for preconditioner\n",*BILU);
     
     HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &Atilde);
     HYPRE_IJMatrixSetObjectType(Atilde, HYPRE_PARCSR);
     HYPRE_IJMatrixInitialize(Atilde);
     
     rcols = (int *)malloc( csize*sizeof(int) );
     dbuf = (double *)malloc( csize*sizeof(double) );
     for (i = 0; i < local_size; i++) {
       irow=globaldofs[i];
       nnz = 0;
       for (j=rows[i];j<rows[i+1];j++) {
         jcol = globaldofs[cols[j-1]-1];
         /*TODO - is the block ordering preserved in the linear numbering?
	   Here we assume it is.
	 */
         if ((irow%*BILU)==(jcol%*BILU)) {
	   rcols[nnz] = jcol;
           dbuf[nnz] = vals[j-1];
           nnz++;
	 }
       }
       HYPRE_IJMatrixAddToValues(Atilde, 1, &nnz, &irow, rcols, dbuf);
     }
     free( rcols );
     free( dbuf );
     /* Assemble after setting the coefficients */
     HYPRE_IJMatrixAssemble(Atilde);     
   }

   /* Get the parcsr matrix object to use */
   /* note: this is only used for setup,  */
   /* so we put in the possibly approxima-*/
   /* ted matrix Atilde                   */
   HYPRE_IJMatrixGetObject(Atilde, (void**) &parcsr_A);
   
   HYPRE_IJVectorCreate(comm, ilower, iupper,&b);
   HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(b);
   HYPRE_IJVectorAssemble(b);
   HYPRE_IJVectorGetObject(b, (void **) &par_b);

   HYPRE_IJVectorCreate(comm, ilower, iupper,&x);
   HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(x);
   HYPRE_IJVectorAssemble(x);
   HYPRE_IJVectorGetObject(x, (void **) &par_x);

   if(myid==0 && verbosity >= 12) {
     fprintf(stderr,"HYRPE INT: %d %d  %d %d %d \n", 
	     hypre_intpara[0], hypre_intpara[1], hypre_intpara[2], hypre_intpara[3], hypre_intpara[4]);  
     fprintf(stderr,"HYRPE DP: %d %d %d %d %d \n", &
	     hypre_dppara[0], hypre_dppara[1], hypre_dppara[2], hypre_dppara[3], hypre_dppara[4]);  
   }
   /* Choose a solver and solve the system */
   /* NB.: hypremethod 
      0 ... BiCGStab + ILUn
      1 ... BiCGStab + ParaSails
      2 ... BiCGStab + BoomerAMG
      10 ... BoomerAMG 
      20 ... CG + ILUn 
      21 ... CG + ParaSails 
      22 ... CG + BoomerAMG
      30 ... GMRes + ILUn 
      31 ... GMRes + ParaSails 
      32 ... GMRes + BoomerAMG    (This choice disabled currently)

      with the flag HAVE_GMRES additionally
      40 ... FlexGMRes + ILUn 
      41 ... FlexGMRes + ParaSails 
      42 ... FlexGMRes + BoomerAMG
      50 ... LGMRes + ILUn 
      51 ... LGMRes + ParaSails 
      52 ... LGMRes + BoomerAMG


   */
   
   /* create preconditioner for Krylov methods */
   /* for Boomer as solver we create it as a   */
   /* preconditioner here and set the pointer  */
   if ( hypre_pre == 0) {
     HYPRE_EuclidCreate(comm, &precond );
     static char *argv[5], str[3];
     argv[0] = "-level";
     sprintf( str, "%d", *ILUn );
     if (myid == 0 && verbosity >= 4) fprintf( stderr,"SolveHypre: using ILU%i as preconditioner\n",*ILUn); 
     argv[1] = str;
     HYPRE_EuclidSetParams( precond, 2, argv );
   }

   else if ( hypre_pre == 1 ) {
     if (myid == 0 && verbosity >= 4) fprintf( stderr,"SolveHypre: using ParaSails as preconditioner\n"); 

     /* Now set up the ParaSails preconditioner and specify any parameters */
     HYPRE_ParaSailsCreate(comm, &precond);

     /* Set some parameters (See Reference Manual for more parameters) */
     /* threshold = dppara[0]; maxlevels= intpara[1] */
     HYPRE_ParaSailsSetParams(precond, hypre_dppara[0], hypre_intpara[1]);
     /* filter = dppara[1] */
     HYPRE_ParaSailsSetFilter(precond, hypre_dppara[1]);
     /* symmetry = intpara[0] */
     HYPRE_ParaSailsSetSym(precond, hypre_intpara[0]);

     i = 3*(verbosity >= 6 );
     HYPRE_ParaSailsSetLogging(precond, i);
   }

   else if ( hypre_pre == 2 || hypre_sol == 1 )  {
     if (myid == 0 ) {
       if( verbosity >= 5 ) {
	 fprintf( stderr,"SolveHypre: using BoomerAMG\n");
       }
       if( verbosity >= 10 ) {
	 fprintf( stderr,"RelaxType = %d\n",hypre_intpara[0]); 
	 fprintf( stderr,"CoarsenType = %d\n",hypre_intpara[1]); 
	 fprintf( stderr,"NumSweeps = %d\n",hypre_intpara[2]); 
	 fprintf( stderr,"MaxLevels = %d\n",hypre_intpara[3]); 
	 fprintf( stderr,"Interpolation Type = %d\n",hypre_intpara[4]); 
	 fprintf( stderr,"Smooth Type = %d\n",hypre_intpara[5]);
	 fprintf( stderr,"Cycle Type = %d\n",hypre_intpara[6]);
	 fprintf( stderr,"DOFs = %d\n",hypre_intpara[7]);
	 fprintf( stderr,"StrongThreshold = %g\n",hypre_dppara[0]);
       }
     }
     HYPRE_BoomerAMGCreate(&precond);
     /* Set some parameters (See Reference Manual for more parameters) */
     HYPRE_BoomerAMGSetNumFunctions(precond, hypre_intpara[7]); /* No. of PDE's */

     i = (verbosity >= 6);
     HYPRE_BoomerAMGSetPrintLevel(precond, i); /* print amg solution info */

     HYPRE_BoomerAMGSetNumSweeps(precond, hypre_intpara[2]); /* fixed for preconditioner to 1 */
     HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
     HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */
     HYPRE_BoomerAMGSetRelaxType(precond, hypre_intpara[0]);   /* G-S/Jacobi hybrid relaxation */
     HYPRE_BoomerAMGSetCoarsenType(precond, hypre_intpara[1]);  /* coarsening type */
     
     HYPRE_BoomerAMGSetMaxLevels(precond, hypre_intpara[3]); /* levels of coarsening */
     HYPRE_BoomerAMGSetInterpType(precond, hypre_intpara[4]);  /* interpolation type */
     HYPRE_BoomerAMGSetSmoothType(precond, hypre_intpara[5]);  /* smoother type */
     HYPRE_BoomerAMGSetCycleType(precond, hypre_intpara[6]);  /* coarsening type */
     /* threshold for strong coupling (default 0.25 recommended for 2D Laplace, 0.5-0.6 
	for 3D Laplace, 0.9 for elasticity) */
     HYPRE_BoomerAMGSetStrongThreshold(precond, hypre_dppara[0]);  	 

   } else if ( hypre_pre != 9 ) {
     fprintf( stderr,"Hypre preconditioning method not implemented\n");
     exit(EXIT_FAILURE);
   }

   /* create solver */
   if ( hypre_sol == 0) { /* BiGSTAB methods */
     /* Create solver */
     HYPRE_ParCSRBiCGSTABCreate(comm, &solver);

     /* Set some parameters (See Reference Manual for more parameters) */
     HYPRE_ParCSRBiCGSTABSetStopCrit(solver, 0);     /* use the two norm as the stopping criteria */
     HYPRE_ParCSRBiCGSTABSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_ParCSRBiCGSTABSetTol(solver, *TOL);       /* conv. tolerance */

     i = 2*(verbosity >= 6);
     HYPRE_ParCSRBiCGSTABSetPrintLevel(solver, i);   /* print solve info */
     i = (verbosity >= 6);
     HYPRE_ParCSRBiCGSTABSetLogging(solver, i);      /* needed to get run info later */


     /* Set the BiCGStabl preconditioner */
     if ( hypre_pre == 0) {
       HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup, precond);
     } else if (hypre_pre == 1) { 
       HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
     } else if(hypre_pre == 2) {
       HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
     }
     /* compute the preconditioner */
     if (myid==0 && verbosity >= 5 ) fprintf(stdout,"create preconditioner...");
     HYPRE_ParCSRBiCGSTABSetup(solver, parcsr_A, par_b, par_x);
     
   } else if ( hypre_sol == 1 ) { /* boomer AMG */
     int num_iterations;
     double final_res_norm;
     
     solver = precond;
     precond = NULL;
     
     /* Now setup - note that the input vectors are ignored so we can pass in NULLs */
     if (myid==0 && verbosity >= 5 ) {
       fprintf(stdout,"construct BoomerAMG solver");
     }
     HYPRE_BoomerAMGSetTol(solver, *TOL);      /* conv. tolerance */
     HYPRE_BoomerAMGSetMaxIter(solver, *Rounds); /* iteration rounds */
     HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
   }

   else if ( hypre_sol == 2) { /* CG */
     /* Create solver */
     HYPRE_ParCSRPCGCreate(comm, &solver);
     
     HYPRE_ParCSRPCGSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_ParCSRPCGSetTol(solver, *TOL);       /* conv. tolerance */

     /* Set some parameters (See Reference Manual for more parameters) */
     HYPRE_ParCSRPCGSetTwoNorm(solver, 1);     /* use the two norm as the stopping criteria */
     i = 2*(verbosity >= 6);
     HYPRE_ParCSRPCGSetPrintLevel(solver, i);   /* print solve info */

     i = (verbosity >= 6);
     HYPRE_ParCSRPCGSetLogging(solver, i);      /* needed to get run info later */
     
     /* Set the PCG preconditioner */
     if ( hypre_pre  == 0) {
       HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup, precond);
     } else if ( hypre_pre == 1) { 
       HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
     } else if( hypre_pre == 2){
       HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
     }
     /* compute the preconditioner */
     if (myid == 0 && verbosity >= 6 ) fprintf(stdout,"create preconditioner...");
     HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
   }

   else if ( hypre_sol == 3) { /* GMRES */
     /* Create solver */
     HYPRE_ParCSRGMRESCreate(comm, &solver);
     HYPRE_GMRESSetMaxIter(solver, *Rounds); /* max GMRES iterations */
     HYPRE_GMRESSetTol(solver, *TOL);        /* GMRES conv. tolerance */
    
     /* Set some parameters (See Reference Manual for more parameters) */
     i = 2*(verbosity >= 6);
     HYPRE_GMRESSetPrintLevel(solver, i);   /* print solve info */

     i = (verbosity >= 6);
     HYPRE_GMRESSetLogging(solver, i);      /* needed to get run info later */

     HYPRE_GMRESSetKDim(solver, hypre_intpara[8]);

     /* Set the GMRES preconditioner */
     /* NOTE: AMG preconditioning is not enabled currently as this choice */
     /* is more feasible in connection with FGMRES, which should be enabled */
     /* separately (TO DO) */
     if ( hypre_pre  == 0) {
       HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve,
			     (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup, precond);
     } else if (hypre_pre == 1) { 
       HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
			     (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
     }
     /* else if( hypre_pre == 2){
	HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
	(HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
	} */

     /* Pass the matrix and rhs into the solver */
     if (myid == 0 && verbosity >= 6 ) fprintf(stdout,"Passing matrix and rhs into GMRES solver...");
     HYPRE_ParCSRGMRESSetup(solver, parcsr_A, par_b, par_x);   
     
   }

#if HAVE_GMRES
   else if ( hypre_sol == 4) { /* FlexGMRes */
     /* Create solver */
     HYPRE_ParCSRFlexGMRESCreate(comm, &solver);
     HYPRE_ParCSRFlexGMRESSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_ParCSRFlexGMRESSetTol(solver, *TOL);       /* conv. tolerance */
     
     /* Set some parameters (See Reference Manual for more parameters) */
     i = 2*(verbosity >= 6);
     HYPRE_FlexGMRESSetPrintLevel(solver, i);   /* print solve info */

     i = (verbosity >= 6);
     HYPRE_FlexGMRESSetLogging(solver, i);      /* needed to get run info later */

     HYPRE_FlexGMRESSetKDim(solver,hypre_intpara[8]);

     /* Set the FlexGMRES preconditioner */
     if ( hypre_pre  == 0) {
       HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup, precond);
     } else if ( hypre_pre == 1) { 
       HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
     } else if( hypre_pre == 2){
       HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
     }
     /* Pass the matrix and rhs into the solver */
     if (myid == 0 && verbosity >= 6 ) fprintf(stdout,"Passing matrix and rhs into FlexGMRES solver...");
     HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_A, par_b, par_x);
   }

   else if ( hypre_sol == 5) { /* LGMRes */
     /* Create solver */
     HYPRE_ParCSRLGMRESCreate(comm, &solver);
     HYPRE_ParCSRLGMRESSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_ParCSRLGMRESSetTol(solver, *TOL);       /* conv. tolerance */
     
     /* Set some parameters (See Reference Manual for more parameters) */
     i = 2*(verbosity >= 6);
     HYPRE_ParCSRLGMRESSetPrintLevel(solver, i);   /* print solve info */

     i = (verbosity >= 6);
     HYPRE_ParCSRLGMRESSetLogging(solver, i);      /* needed to get run info later */

     HYPRE_ParCSRLGMRESSetKDim(solver,hypre_intpara[8]);

     /* Set the LGMRES preconditioner */
     if ( hypre_pre  == 0) {
       HYPRE_LGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup, precond);
     } else if ( hypre_pre == 1) { 
       HYPRE_LGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
     } else if( hypre_pre == 2){
       HYPRE_LGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
     }
     /* compute the preconditioner */
     if (myid == 0 && verbosity >= 6 ) fprintf(stdout,"Passing matrix and rhs into LGMRES solver...");
     HYPRE_ParCSRLGMRESSetup(solver, parcsr_A, par_b, par_x);
   }
#endif

   else {
     fprintf( stderr,"Hypre solver method not implemented\n");
     exit(EXIT_FAILURE);
   }

   Container->ilower = ilower;
   Container->iupper = iupper;     
   Container->hypre_method = *hypre_method;
   Container->A = A;
   Container->Atilde = Atilde;
   Container->solver = solver;
   Container->precond = precond;
   
   if( myid == 0 && verbosity >= 6 ) {
     fprintf( stdout, "Hypre setup time: %g\n", realtime_()-st ); 
   }
   
} /* SolveHypre1 - matrix conversion and solver setup */


/* Update the stopping tolerance of a previously constructed solver */
void STDCALLBULL FC_FUNC(updatehypre,UPDATEHYPRE)
     ( double *TOL,  int *hypre_method,  int** ContainerPtr )
{
  HYPRE_Solver solver;
  ElmerHypreContainer *Container;
  int hypre_sol;

  Container = (ElmerHypreContainer*)(*ContainerPtr);
  solver = Container->solver;

  hypre_sol = *hypre_method / 10;

  if ( hypre_sol == 0) { /* BiGSTAB method */
    HYPRE_ParCSRBiCGSTABSetTol(solver, *TOL);
  }
  else if ( hypre_sol == 1 ) { /* boomer AMG */
    HYPRE_BoomerAMGSetTol(solver, *TOL);
  }
  else if ( hypre_sol == 2) { /* CG */
    HYPRE_ParCSRPCGSetTol(solver, *TOL);
  }
  else if ( hypre_sol == 3) { /* GMRES */
    HYPRE_GMRESSetTol(solver, *TOL);
  }
#if HAVE_GMRES
  else if ( hypre_sol == 4) { /* FlexGMRes */
    HYPRE_ParCSRFlexGMRESSetTol(solver, *TOL);
  }
  else if ( hypre_sol == 5) { /* LGMRes */
    HYPRE_ParCSRLGMRESSetTol(solver, *TOL);
  }
#endif
  else {
    fprintf( stderr,"Hypre solver method not implemented\n");
    exit(EXIT_FAILURE);
  }
}


/*////////////////////////////////////////////////////////////////////////////////////////////////*/

/* solve a linear system with previously constructed solver and preconditioner */
void STDCALLBULL FC_FUNC(solvehypre2,SOLVEHYPRE2)
 (
  int *nrows, int *globaldofs, int *owner,  double *xvec,
  double *rhsvec, int *Rounds, double *TOL,
  int *verbosityPtr, int** ContainerPtr, int *fcomm
 )
{

   int i, j, k, *rcols;
   int myid, num_procs;
   int N, n;

   int ilower, iupper;
   int local_size, extra;

   int print_solution, print_system;

   double  *txvec, st, realtime_();
   
   HYPRE_Solver solver, precond;

   HYPRE_ParCSRMatrix parcsr_A;
   HYPRE_IJVector b;
   HYPRE_ParVector par_b;
   HYPRE_IJVector x;
   HYPRE_ParVector par_x;

   MPI_Comm comm=MPI_Comm_f2c(*fcomm);
   
   ElmerHypreContainer *Container;

   int verbosity = *verbosityPtr;
   int hypre_sol, hypre_pre;

   Container = (ElmerHypreContainer*)(*ContainerPtr);

   /* which process number am I? */
   MPI_Comm_rank(comm, &myid);

   if (myid==0 && verbosity >= 6) fprintf(stdout,"HYPRE Solve\n");

   if (Container==NULL) {
     fprintf( stderr, "ID no. %i: pointer passed into SolveHypre2 is NULL, not solving",myid);
     return;
   }
  
   st = realtime_();

   HYPRE_IJMatrixGetObject(Container->A, (void**) &parcsr_A);
   solver = Container->solver;
   precond = Container->precond;

   ilower = Container->ilower;
   iupper = Container->iupper;
   local_size = *nrows;

   hypre_sol = Container->hypre_method / 10;
   hypre_pre = Container->hypre_method % 10;
   /* No preconditioner for BoomerAMG */
   if( hypre_sol == 1 ) hypre_pre = -1;


   /* Create the rhs and solution */
   rcols = (int *)malloc( local_size*sizeof(int) );
   txvec = (double *)malloc( local_size*sizeof(double) );

   for( k=0,i=0; i<local_size; i++ ) rcols[k++] = globaldofs[i];

   for( i=0; i<local_size; i++ ) txvec[i] = rhsvec[i];
   
   HYPRE_IJVectorCreate(comm, ilower, iupper,&b);
   HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(b);

   for( i=0; i<local_size; i++ ) txvec[i] = rhsvec[i];
   HYPRE_IJVectorAddToValues(b, local_size, rcols, txvec );

   HYPRE_IJVectorCreate(comm, ilower, iupper,&x);
   HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(x);

   for( i=0; i<local_size; i++ ) txvec[i] = xvec[i];
   HYPRE_IJVectorSetValues(x, local_size, rcols, txvec );
   
   HYPRE_IJVectorAssemble(b);
   HYPRE_IJVectorGetObject(b, (void **) &par_b);
   
   HYPRE_IJVectorAssemble(x);
   HYPRE_IJVectorGetObject(x, (void **) &par_x);
   
   /* Now setup and solve! */
   if ( hypre_sol == 0) {
//     HYPRE_ParCSRBiCGSTABSetMaxIter(solver, *Rounds); /* max iterations */
//     HYPRE_ParCSRBiCGSTABSetTol(solver, *TOL);       /* conv. tolerance */
     HYPRE_ParCSRBiCGSTABSolve(Container->solver, parcsr_A, par_b, par_x);
   }

   else if ( hypre_sol == 1) {
     int num_iterations;
     double final_res_norm;
     
     HYPRE_BoomerAMGSolve(Container->solver, parcsr_A, par_b, par_x);
     
     /* Run info - needed logging turned on */
     HYPRE_BoomerAMGGetNumIterations(Container->solver, &num_iterations);
     HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
     if (myid == 0 && verbosity <= 5 ) {
       fprintf(stdout,"BoomerAMG:\n");
       fprintf(stdout,"Iterations = %d\n", num_iterations);
       fprintf(stdout,"Final Relative Residual Norm = %e\n", final_res_norm);
       fprintf(stdout,"\n");
     }
   }
 
   else if ( hypre_sol == 2) {
//     HYPRE_ParCSRPCGSetMaxIter(solver, *Rounds); /* max iterations */
//     HYPRE_ParCSRPCGSetTol(solver, *TOL);       /* conv. tolerance */
     HYPRE_ParCSRPCGSolve(Container->solver, parcsr_A, par_b, par_x);
   }

   else if ( hypre_sol == 3) {
//     HYPRE_GMRESSetMaxIter(solver, *Rounds); /* max GMRES iterations */
//     HYPRE_GMRESSetTol(solver, *TOL);        /* GMRES conv. tolerance */
     HYPRE_ParCSRGMRESSolve(Container->solver, parcsr_A, par_b, par_x);
   }

#if HAVE_GMRES
   else if ( hypre_sol == 4) {
//     HYPRE_ParCSRFlexGMRESSetMaxIter(solver, *Rounds); /* max iterations */
//     HYPRE_ParCSRFlexGMRESSetTol(solver, *TOL);       /* conv. tolerance */
     HYPRE_ParCSRFlexGMRESSolve(Container->solver, parcsr_A, par_b, par_x);
   }

   else if ( hypre_sol == 5) {
//     HYPRE_ParCSRLGMRESSetMaxIter(solver, *Rounds); /* max iterations */
//     HYPRE_ParCSRLGMRESSetTol(solver, *TOL);       /* conv. tolerance */
     HYPRE_ParCSRLGMRESSolve(Container->solver, parcsr_A, par_b, par_x);
   }
#endif


   for( k=0,i=0; i<local_size; i++ )
     if ( owner[i] ) rcols[k++] = globaldofs[i];
   
   HYPRE_IJVectorGetValues(x, k, rcols, txvec );
   
   for( i=0,k=0; i<local_size; i++ )
     if ( owner[i] ) xvec[i] = txvec[k++];
   
   if (myid==0 && verbosity >= 5) fprintf( stdout, "solve time: %g\n", realtime_()-st );
   free( txvec );
   free( rcols );
   
   HYPRE_IJVectorDestroy(x);
   HYPRE_IJVectorDestroy(b);
}


/*TODO - add function solvehypre3 that e..g updates the matrix in the
       Container and Krylov solver but leaves the preconditioner   
       unchanged.
*/

/* destroy HYPRE data structure stored in a fortran environment */
void STDCALLBULL FC_FUNC(solvehypre4,SOLVEHYPRE4)(int** ContainerPtr) {

   ElmerHypreContainer* Container = (ElmerHypreContainer*)(*ContainerPtr);

   int hypre_sol, hypre_pre;

   if (Container==0) return;

   hypre_sol = Container->hypre_method / 10;
   hypre_pre = Container->hypre_method % 10;
   if( hypre_sol == 1 ) hypre_pre = -1;


   if ( hypre_sol == 1) { /* boomer AMG */
     /* Destroy solver */
     HYPRE_BoomerAMGDestroy(Container->solver);
   }
   else {
     /* Destroy solver and preconditioner */
     if ( hypre_sol == 0) {
       HYPRE_ParCSRBiCGSTABDestroy(Container->solver);
     }
     else if ( hypre_sol == 2 ) {
       HYPRE_ParCSRPCGDestroy(Container->solver);
     }
     else if ( hypre_sol == 3 ) {
       HYPRE_ParCSRGMRESDestroy(Container->solver);
     }
#if HAVE_GMRES
     else if ( hypre_sol == 4 ) {
       HYPRE_ParCSRFlexGMRESDestroy(Container->solver);
     }
     else if ( hypre_sol == 5 ) {
       HYPRE_ParCSRLGMRESDestroy(Container->solver);
     }
#endif

     if ( hypre_pre == 0 ) {
       HYPRE_EuclidDestroy(Container->precond);
     } else if ( hypre_pre == 1 ) {
       HYPRE_ParaSailsDestroy(Container->precond);
     } else if ( hypre_pre == 2 ) {
       HYPRE_BoomerAMGDestroy(Container->precond);
     }
   }

   if (Container->Atilde != Container->A) {
     HYPRE_IJMatrixDestroy(Container->Atilde);
   }
   free(Container);
   *ContainerPtr = NULL;
}

void STDCALLBULL FC_FUNC(solvehypreams,SOLVEHYPREAMS)
 (
  int *nrows,int *rows, int *cols, double *vals, int *nnodes,
  int *grows, int *gcols, double *gvals, int *perm,
  int *invperm, int *globaldofs, int *owner,  int *globalnodes, 
  int *nodeowner, double *xvec,
  double *rhsvec, int *pe, int *ILUn, int *Rounds, double *TOL,
  double *xx_d, double *yy_d, double *zz_d, 
  int *hypre_method, int *hypre_intpara, double *hypre_dppara,
  int *verbosityPtr, int** ContainerPtr, int *fcomm
 )
{
   int i, j, k, *rcols;
   int myid, num_procs;
   int N, n;

   int ilower, iupper, nlower, nupper;
   int local_size, local_nodes, extra;

   int solver_id;
   int print_solution, print_system;

   double  *txvec, st, realtime_();

   HYPRE_ParCSRMatrix parcsr_A,parcsr_G;
   HYPRE_IJMatrix A,G;
   HYPRE_IJVector b;
   HYPRE_ParVector par_b;
   HYPRE_IJVector x;
   HYPRE_ParVector par_x;
   HYPRE_IJVector xx,yy,zz;
   HYPRE_ParVector par_xx,par_yy,par_zz;

   HYPRE_Solver solver, precond;
   int verbosity = 10;
   MPI_Comm comm=MPI_Comm_f2c(*fcomm);
   
   st  = realtime_();
   /* How many rows do I have? */
   local_size = *nrows;
   local_nodes = *nnodes;

   ilower=1000000000;
   iupper=0;
   for( i=0; i<local_size; i++ )
   {
      if ( owner[i] ) {
        if ( iupper < globaldofs[i] ) iupper=globaldofs[i];
        if ( ilower > globaldofs[i] ) ilower=globaldofs[i];
      }
   }


   nlower=1000000000;
   nupper=0;
   for( i=0; i<local_nodes; i++ )
   {
      if ( nodeowner[i] ) {
        if ( nupper < globalnodes[i] ) nupper=globalnodes[i];
        if ( nlower > globalnodes[i] ) nlower=globalnodes[i];
      }
   }

  /* which process number am I? */
   MPI_Comm_rank(comm, &myid);
   /* Create the matrix.
      Note that this is a square matrix, so we indicate the row partition
      size twice (since number of rows = number of cols) */
   HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &A);
   HYPRE_IJMatrixCreate(comm, ilower, iupper, nlower, nupper, &G);

   /* Choose a parallel csr format storage (see the User's Manual) */
   HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
   HYPRE_IJMatrixSetObjectType(G, HYPRE_PARCSR);

   /* Initialize before setting coefficients */
   HYPRE_IJMatrixInitialize(A);
   HYPRE_IJMatrixInitialize(G);

   /* Now go through my local rows and set the matrix entries.
      Note that here we are setting one row at a time, though
      one could set all the rows together (see the User's Manual).
   */
   {
      int nnz,irow,i,j,k,*rcols,csize=32;

      rcols = (int *)malloc( csize*sizeof(int) );
      for (i = 0; i < local_size; i++)
      {
         nnz = rows[i+1]-rows[i];
         if ( nnz>csize ) {
           rcols = (int *)realloc( rcols, nnz*sizeof(int) );
           csize = nnz;
         }
         irow=globaldofs[i];
         for( k=0,j=rows[i]; j<rows[i+1]; j++,k++)
         {
           rcols[k] = globaldofs[cols[j-1]-1];
         }
         HYPRE_IJMatrixAddToValues(A, 1, &nnz, &irow, rcols, &vals[rows[i]-1]);
      }
        free( rcols );
   }

   {
      int nnz,irow,i,j,k,*rcols,csize=32;

      rcols = (int *)malloc( csize*sizeof(int) );
      for (i = 0; i < local_size; i++)
      {
         nnz = grows[i+1]-grows[i];
         if ( nnz>csize ) {
           rcols = (int *)realloc( rcols, nnz*sizeof(int) );
           csize = nnz;
         }
         irow=globaldofs[i];
         for( k=0,j=grows[i]; j<grows[i+1]; j++,k++)
         {
           rcols[k] = globalnodes[gcols[j-1]-1];
         }
         HYPRE_IJMatrixAddToValues(G, 1, &nnz, &irow, rcols, &gvals[grows[i]-1]);
      }
      free( rcols );
   }

   /* Assemble after setting the coefficients */
   HYPRE_IJMatrixAssemble(A);
   HYPRE_IJMatrixAssemble(G);

   /* Get the parcsr matrix object to use */
   HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
   HYPRE_IJMatrixGetObject(G, (void**) &parcsr_G);

   /* Create the rhs and solution */
   rcols = (int *)malloc( local_size*sizeof(int) );
   txvec = (double *)malloc( local_size*sizeof(double) );
   for( k=0,i=0; i<local_size; i++ ) rcols[k++] = globaldofs[i];

   HYPRE_IJVectorCreate(comm, ilower, iupper,&b);
   HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(b);
   HYPRE_IJVectorAddToValues(b, local_size, rcols, rhsvec);
   HYPRE_IJVectorAssemble(b);
   HYPRE_IJVectorGetObject(b, (void **) &par_b);

   HYPRE_IJVectorCreate(comm, ilower, iupper,&x);
   HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(x);
   HYPRE_IJVectorSetValues(x, local_size, rcols, xvec);
   HYPRE_IJVectorAssemble(x);
   HYPRE_IJVectorGetObject(x, (void **) &par_x);

#if 0
   for( k=0,i=0; i<local_nodes; i++ ) rcols[k++] = globalnodes[i];

   HYPRE_IJVectorCreate(comm, nlower, nupper,&xx);
   HYPRE_IJVectorSetObjectType(xx, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(xx);
   HYPRE_IJVectorSetValues(xx, local_nodes, rcols,xx_d);
   HYPRE_IJVectorAssemble(xx);
   HYPRE_IJVectorGetObject(xx, (void **) &par_xx);

   HYPRE_IJVectorCreate(comm, nlower, nupper,&yy);
   HYPRE_IJVectorSetObjectType(yy, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(yy);
   HYPRE_IJVectorSetValues(yy, local_nodes, rcols, yy_d);
   HYPRE_IJVectorAssemble(yy);
   HYPRE_IJVectorGetObject(yy, (void **) &par_yy);

   HYPRE_IJVectorCreate(comm, nlower, nupper,&zz);
   HYPRE_IJVectorSetObjectType(zz, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(zz);
   HYPRE_IJVectorSetValues(zz, local_nodes, rcols, zz_d);
   HYPRE_IJVectorAssemble(zz);
   HYPRE_IJVectorGetObject(zz, (void **) &par_zz);
#else
   HYPRE_IJVectorCreate(comm, ilower, iupper,&xx);
   HYPRE_IJVectorSetObjectType(xx, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(xx);
   HYPRE_IJVectorSetValues(xx, local_size, rcols, xx_d);
   HYPRE_IJVectorAssemble(xx);
   HYPRE_IJVectorGetObject(xx, (void **) &par_xx);

   HYPRE_IJVectorCreate(comm, ilower, iupper,&yy);
   HYPRE_IJVectorSetObjectType(yy, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(yy);
   HYPRE_IJVectorSetValues(yy, local_size, rcols, yy_d);
   HYPRE_IJVectorAssemble(yy);
   HYPRE_IJVectorGetObject(yy, (void **) &par_yy);

   HYPRE_IJVectorCreate(comm, ilower, iupper,&zz);
   HYPRE_IJVectorSetObjectType(zz, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(zz);
   HYPRE_IJVectorSetValues(zz, local_size, rcols, zz_d);
   HYPRE_IJVectorAssemble(zz);
   HYPRE_IJVectorGetObject(zz, (void **) &par_zz);
#endif

   if( verbosity >= 12 ) fprintf( stderr, "ID no. %i: setup time: %g\n", myid, realtime_()-st );
   st = realtime_();

   if (myid == 0 ) 
      if ( verbosity >= 5 )
        fprintf( stderr,"SolveHypre: using BiCGStab + AMS\n");

   HYPRE_AMSCreate(&precond); 
   HYPRE_AMSSetMaxIter(precond,1);
   HYPRE_AMSSetDiscreteGradient(precond,parcsr_G);
// HYPRE_AMSSetCoordinateVectors(precond,par_xx,par_yy,par_zz);
   HYPRE_AMSSetEdgeConstantVectors(precond,par_xx,par_yy,par_zz);

   HYPRE_AMSSetCycleType(precond, hypre_intpara[6]);// 1-8
   HYPRE_AMSSetSmoothingOptions(precond, hypre_intpara[5], hypre_intpara[2], 1.0, 1.0);
// HYPRE_AMSSetAlphaAMGOptions(precond, 10, 1, 3, 0.25);

// HYPRE_AMSSetBetaAMGOptions(precond, 10, 1, 3, 0.25);
// HYPRE_AMSSetBetaPoissonMatrix(precond,NULL);

    /* Create solver */
     HYPRE_ParCSRBiCGSTABCreate(comm, &solver);

     /* Set some parameters (See Reference Manual for more parameters) */
     HYPRE_ParCSRBiCGSTABSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_ParCSRBiCGSTABSetTol(solver, *TOL);       /* conv. tolerance */
     HYPRE_ParCSRBiCGSTABSetStopCrit(solver, 0);     /* use the two norm as the stopping criteria */
     HYPRE_ParCSRBiCGSTABSetPrintLevel(solver, 2);   /* print solve info */
     HYPRE_ParCSRBiCGSTABSetLogging(solver, 1);      /* needed to get run info later */

   /* Set the BiCGSTAB preconditioner */
   HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_AMSSolve,
               (HYPRE_PtrToSolverFcn) HYPRE_AMSSetup, precond);

   /* Now setup and solve! */
   HYPRE_ParCSRBiCGSTABSetup(solver, parcsr_A, par_b, par_x);
   HYPRE_ParCSRBiCGSTABSolve(solver, parcsr_A, par_b, par_x);

   /* Destroy solver and preconditioner */
   HYPRE_AMSDestroy(precond);
   HYPRE_ParCSRBiCGSTABDestroy(solver);

   for( k=0,i=0; i<local_size; i++ )
     if ( owner[i] ) rcols[k++] = globaldofs[i];
   
   HYPRE_IJVectorGetValues(x, k, rcols, txvec );
   
   for( i=0,k=0; i<local_size; i++ )
     if ( owner[i] ) xvec[i] = txvec[k++];
   
   if( myid == 0 && verbosity >= 5 ) {
     fprintf( stderr, "Hypre solve time: %g\n", realtime_()-st );
   }
   free( txvec );
   free( rcols );

   /* Clean up */
   HYPRE_IJMatrixDestroy(A);
   HYPRE_IJMatrixDestroy(G);
   HYPRE_IJVectorDestroy(b);
   HYPRE_IJVectorDestroy(x);
}


#endif
