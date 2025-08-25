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

#define realtime_ FC_FUNC_(realtime,REALTIME)
typedef struct {

  int ilower, iupper;

  HYPRE_IJMatrix A;
  HYPRE_IJMatrix Atilde;

  HYPRE_IJVector x,b;
  
  int hypre_method;
  HYPRE_Solver solver, precond;

  /* AMS specific stuff */
  HYPRE_IJMatrix G, Pi;

} ElmerHypreContainer;


/* The interface calls Hypre step-wise separating phases for 
   procedure of setup, solve and cleanup.
   The original monolithic one that was not used for years has been removed. 

   TO DO: we should add the possibility to keep the precon-
   ditioner the same but update the system matrix (SolveHYPRE3), right now
   calling SolveHYPRE2 solves with the matrix passed into
   SolveHYPRE1. */
   
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

   int num_iterations;
   double final_res_norm;
  
   double  *txvec, st, realtime_();
   
   int verbosity = *verbosityPtr, myverb;

   /* which process number am I? */
   MPI_Comm_rank(comm, &myid);

   if( myid == 0 ) 
     myverb = verbosity;
   else
     myverb = 0;
   
   if (myverb > 8) fprintf(stdout,"SolveHypre: Performing HYPRE Setup\n");
   
   st  = realtime_();
   
   /* How many rows do I have? */
   local_size = *nrows;
   hypre_sol = *hypre_method / 100;
   hypre_pre = *hypre_method % 100;

   if(hypre_sol == 2 || hypre_pre == 2)  {
     if (*ContainerPtr == NULL) {
       fprintf( stdout, "Hypre pointer should not be zero at start for AMS\n");
     }   
     Container = (ElmerHypreContainer*)(*ContainerPtr);
   }
   else {   
     Container = (ElmerHypreContainer*)malloc(sizeof(ElmerHypreContainer));   
     *ContainerPtr=(int*)(Container);
   }
   
   ilower =  1000000000;
   iupper = -1;
   for( i=0; i<local_size; i++ ) {
     if ( owner[i] ) {
       if ( iupper < globaldofs[i] ) iupper = globaldofs[i];
       if ( ilower > globaldofs[i] ) ilower = globaldofs[i];
     }
   }


   /* if the partition doesn't own any of the dofs, apply null range (with valid indices) */
   if ( iupper == -1 ) { ilower = 1; iupper = 0; }

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
     /* Standard version - use A as preconditioner */
     Atilde = A;
   } else if ( *precflag ) {
     /* We have another matrix that is used as preconditioning matrix! */
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
     /* We only take the block diagonal values of the original matrix for our preconditioner */
     int nnz,irow,jcol,i,j,k,*rcols;
     double *dbuf;
     if (myverb > 6) fprintf(stdout,"HYPRE: using BILU(%d) approximation for preconditioner\n",*BILU);
     
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

   if(myverb > 12) {     
     fprintf(stdout,"SolveHypre: NonZero Hypre parameters\n");
     for(i=0;i<=12;i++) {
       if(hypre_intpara[i]) fprintf(stdout,"  intpara %d: %d\n",i+1,hypre_intpara[i]);
     }
     for(i=0;i<=5;i++) {
       if(fabs(hypre_dppara[i]) > 1.0e-20) fprintf(stdout,"   dppara %d: %lg\n",i+1,hypre_dppara[i]);
     }
   }
   
    /* This is copy-pasted from SParIterSolver for convenience.
    !---------------------------------------
    !              No  Prec
    ! none         0    x   -
    ! BoomerAMG    1    x   x 
    ! AMS          2    x   x
    ! ILU          3    x   x
    ! Parasails    4    x   -
    ! FSAI         5    x   -
    ! PCG          6    -   x
    ! BiCGStab     7    -   x
    ! GMRes        8    -   x
    ! FlexGMRes    9    -   x
    ! LGMRes       10   -   x
    ! COGMRes      11   -   x
    !---------------------------------------  */  

   
   /* Create preconditioner for Krylov methods. */
   /* Some methods may act as preconditioners and solvers and are also intialized here
      and later changes the pointer to solvers. */

   if ( hypre_pre == 1 || hypre_sol == 1 )  {
     if( myverb > 8 ) 
       fprintf( stdout,"SolveHypre: Creating Hypre BoomerAMG\n");
     if( myverb > 10 ) {
       fprintf( stdout,"RelaxType = %d\n",hypre_intpara[0]); 
       fprintf( stdout,"CoarsenType = %d\n",hypre_intpara[1]); 
       fprintf( stdout,"NumSweeps = %d\n",hypre_intpara[2]); 
       fprintf( stdout,"MaxLevels = %d\n",hypre_intpara[3]); 
       fprintf( stdout,"Interpolation Type = %d\n",hypre_intpara[4]); 
       fprintf( stdout,"Smooth Type = %d\n",hypre_intpara[5]);
       fprintf( stdout,"Cycle Type = %d\n",hypre_intpara[6]);
       fprintf( stdout,"DOFs = %d\n",hypre_intpara[7]);
       fprintf( stdout,"StrongThreshold = %g\n",hypre_dppara[0]);
     }
     HYPRE_BoomerAMGCreate(&precond);
     /* Set some parameters (See Reference Manual for more parameters) */
     HYPRE_BoomerAMGSetNumFunctions(precond, hypre_intpara[7]); /* No. of PDE's */

     i = (verbosity >= 6);
     if(verbosity >= 10) i=3;
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
     
     if( myverb > 10) fprintf(stdout,"SolveHypre: Created BoomerAMG preconditioner!\n");

   } else if ( hypre_pre == 2 || hypre_sol == 2 ) {
     precond = Container->precond;
     if(precond) {
       if(myverb > 10) fprintf( stdout,"SolveHypre: Using previously defined AMS!\n");
     } else {
       fprintf( stdout,"SolveHypre: Pointer to AMS preconditioner is NULL!\n");       
       exit(EXIT_FAILURE);
     }
   
   } else if ( hypre_pre == 3 || hypre_sol == 3 ) {
     int ilu_type, max_iter, reordering, print_level;
     int max_nnz_row, schur_max_iter,tri_solve, ljac_iters,ujac_iters; 
     double tol, threshold;
       
     if (myverb > 8) fprintf( stdout,"SolveHypre: using ILU%d as preconditioner\n",*ILUn); 
          
     /* (Required) Create ILU solver */
     HYPRE_ILUCreate(&precond);

     /* (Recommended) General solver options */
     HYPRE_ILUSetType(precond, ilu_type=*ILUn); /* 0, 1, 10, 11, 20, 21, 30, 31, 40, 41, 50 */
     HYPRE_ILUSetMaxIter(precond, max_iter=1);
#if 0 
     HYPRE_ILUSetTol(precond, tol);
     HYPRE_ILUSetLocalReordering(precond, reordering); /* 0: none, 1: RCM */
     HYPRE_ILUSetPrintLevel(precond, print_level);

     /* (Optional) Function calls for ILUK variants */
     HYPRE_ILUSetLevelOfFill(precond, fill);

     /* (Optional) Function calls for ILUT variants */
     HYPRE_ILUSetMaxNnzPerRow(precond, max_nnz_row);
     HYPRE_ILUSetDropThreshold(precond, threshold);

     /* (Optional) Function calls for GMRES-ILU or NSH-ILU */
     HYPRE_ILUSetNSHDropThreshold(precond, threshold);
     HYPRE_ILUSetSchurMaxIter(precond, schur_max_iter);

     /* (Optional) Function calls for iterative ILU variants */
     HYPRE_ILUSetTriSolve(precond, tri_solve);
     HYPRE_ILUSetLowerJacobiIters(precond, ljac_iters);
     HYPRE_ILUSetUpperJacobiIters(precond, ujac_iters);
     
     /* (Exclusively required) Function calls for using ILU as standalone solver */
     HYPRE_ILUSetup(precond, parcsr_M, b, x);
     HYPRE_ILUSolve(precond, parcsr_A, b, x);

     /* (Exclusively required) Function calls for using ILU as preconditioner to GMRES */
#endif

   }

   else if ( hypre_pre == 4 ) {
     if (myverb > 8) fprintf( stdout,"SolveHypre: using ParaSails as preconditioner\n"); 

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

   } else if ( hypre_pre == 5 ) {
     int max_steps=5,max_step_size=3;
     double kap_tolerance=1.0e-3;
#if 0
     if (myverb > 8) fprintf( stdout,"SolveHypre: using SPAI as preconditioner\n"); 
     HYPRE_FSAICreate(&precond);
     HYPRE_FSAISetMaxSteps(precond, max_steps);
     HYPRE_FSAISetMaxStepSize(precond, max_step_size);
     HYPRE_FSAISetKapTolerance(precond, kap_tolerance);     
     i = 3*(verbosity >= 6 );
     HYPRE_SPAISetLogging(precond, i);
#else
     fprintf( stdout,"Hypre preconditioning method FSAI not compiled with!\n");
     exit(EXIT_FAILURE);
#endif  
   } else if( hypre_pre != 0) {
     fprintf( stdout,"Hypre preconditioning method %d not implemented\n",hypre_pre);
     exit(EXIT_FAILURE);
   }

   /* Create Hypre solver */

   if ( hypre_sol == 1 ) { /* boomer AMG */     
     solver = precond;
     precond = NULL;
     /* Now setup - note that the input vectors are ignored so we can pass in NULLs */
     if (myverb > 10 ) fprintf(stdout,"SolveHypre: construct BoomerAMG solver");

     /* When BoomerAMG is used as a solver these override the definitions related to AMS as a preconditioner. */
     HYPRE_BoomerAMGSetTol(solver, *TOL);      
     HYPRE_BoomerAMGSetMaxIter(solver, *Rounds); 

     HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
     

   } else if ( hypre_sol == 2 ) { /* AMS */
     solver = precond;
     precond = NULL;
     
     /* Now setup - note that the input vectors are ignored so we can pass in NULLs */
     if (myverb > 10 ) fprintf(stdout,"SolveHypre: construct AMS solver");

     /* When AMS is used as a solver these override the definitions related to AMS as a preconditioner. */
     HYPRE_AMSSetMaxIter(solver,*Rounds);
     HYPRE_AMSSetTol(solver, *TOL);

     HYPRE_AMSSetup(solver, parcsr_A, par_b, par_x);


   } else if ( hypre_sol == 6) { /* PCG */
     /* Create solver */
     HYPRE_ParCSRPCGCreate(comm, &solver);
     
     HYPRE_PCGSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_PCGSetTol(solver, *TOL);        /* conv. tolerance */

     /* Set some parameters (See Reference Manual for more parameters) */
     HYPRE_PCGSetAbsoluteTol(solver, hypre_dppara[5]);       

     HYPRE_PCGSetTwoNorm(solver,hypre_intpara[10]);  /* use the two norm as the stopping criteria */
     //HYPRE_PCGSetFlex(solver,hypre_intpara[11]);     /* use flexible method for robustness */
     
     i = (verbosity >= 6);
     if(verbosity >= 10) i=3;
     HYPRE_PCGSetPrintLevel(solver, i);   /* print solve info */

     i = (verbosity >= 6);
     HYPRE_PCGSetLogging(solver, i);      /* needed to get run info later */
     
     /* Set the PCG preconditioner */
     if( hypre_pre == 1){
       HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
     } else if(hypre_pre == 2 ) {
       HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_AMSSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_AMSSetup, precond);
     } else if ( hypre_pre == 3) {
       HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ILUSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_ILUSetup, precond);
     } else if ( hypre_pre == 4) { 
       HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
     }
       
     if (myverb > 12 ) fprintf(stdout,"SolveHypre: Setting up PCG linear system");
     HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);


   } else if ( hypre_sol == 7) { /* BiGSTAB methods */
     /* Create solver */
     HYPRE_ParCSRBiCGSTABCreate(comm, &solver);

     HYPRE_BiCGSTABSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_BiCGSTABSetTol(solver, *TOL);       /* conv. tolerance */

     /* Set some parameters (See Reference Manual for more parameters) */
     HYPRE_BiCGSTABSetAbsoluteTol(solver, hypre_dppara[5]);       
     
     i = (verbosity >= 6);
     if(verbosity >= 10) i=3;
     HYPRE_BiCGSTABSetPrintLevel(solver, i);   /* print solve info */
     i = (verbosity >= 6);
     HYPRE_BiCGSTABSetLogging(solver, i);      /* needed to get run info later */
 
     /* Set the BiCGStabl preconditioner */
     if(hypre_pre == 1) {
       HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
     } else if(hypre_pre == 2 ) {
       HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_AMSSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_AMSSetup, precond);
     } else if ( hypre_pre == 3) {
       HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ILUSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_ILUSetup, precond);
     } else if (hypre_pre == 4) { 
       HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
     }
     if (myverb > 12 ) fprintf(stdout,"SolveHypre: Setting up BiCGStab linear system");
     HYPRE_ParCSRBiCGSTABSetup(solver, parcsr_A, par_b, par_x);     

     
   } else if ( hypre_sol == 8) { /* GMRES */
     /* Create solver */
     HYPRE_ParCSRGMRESCreate(comm, &solver);

     HYPRE_GMRESSetMaxIter(solver, *Rounds); /* max GMRES iterations */
     HYPRE_GMRESSetTol(solver, *TOL);        /* GMRES conv. tolerance */

     HYPRE_GMRESSetKDim(solver, hypre_intpara[9]);
     HYPRE_GMRESSetAbsoluteTol(solver, hypre_dppara[5]);       
     
     
     /* Set some parameters (See Reference Manual for more parameters) */
     i = (verbosity >= 6);
     if(verbosity >= 10) i=3;
     HYPRE_GMRESSetPrintLevel(solver, i);   /* print solve info */
     i = (verbosity > 6);
     HYPRE_GMRESSetLogging(solver, i);      /* needed to get run info later */

     if( hypre_pre == 1){
       HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
			     (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
     } else if(hypre_pre == 2 ) {
       HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_AMSSolve,
			     (HYPRE_PtrToSolverFcn) HYPRE_AMSSetup, precond);
     } else if ( hypre_pre  == 3) {
       HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ILUSolve,
			     (HYPRE_PtrToSolverFcn) HYPRE_ILUSetup, precond);
     } else if (hypre_pre == 4) { 
       HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
			     (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
     }
       
     if (myverb > 12 ) fprintf(stdout,"SolveHypre: Setting up GMRes linear system");
     HYPRE_ParCSRGMRESSetup(solver, parcsr_A, par_b, par_x);   

     
   } else if ( hypre_sol == 9) { /* FlexGMRes */
     /* Create solver */
     HYPRE_ParCSRFlexGMRESCreate(comm, &solver);

     HYPRE_FlexGMRESSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_FlexGMRESSetTol(solver, *TOL);       /* conv. tolerance */

     HYPRE_FlexGMRESSetKDim(solver,hypre_intpara[9]);
     HYPRE_FlexGMRESSetAbsoluteTol(solver, hypre_dppara[5]);       
     
     /* Set some parameters (See Reference Manual for more parameters) */
     i = (verbosity >= 6);
     if(verbosity >= 10) i=3;
     HYPRE_FlexGMRESSetPrintLevel(solver, i);   /* print solve info */
     i = (verbosity >= 6);
     HYPRE_FlexGMRESSetLogging(solver, i);      /* needed to get run info later */

     /* Set the FlexGMRES preconditioner */
     if( hypre_pre == 1) {
       HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
     } else if(hypre_pre == 2 ) {
       HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_AMSSolve,
				 (HYPRE_PtrToSolverFcn) HYPRE_AMSSetup, precond);
     } else if ( hypre_pre == 3) {
       HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ILUSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_ILUSetup, precond);
     } else if ( hypre_pre == 4) { 
       HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
     }
     if (myverb > 12 ) fprintf(stdout,"SolveHypre: Setting up FlexGMRes linear system");
     HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_A, par_b, par_x);

     
   } else if ( hypre_sol == 10) { /* LGMRes */
     /* Create solver */
     HYPRE_ParCSRLGMRESCreate(comm, &solver);

     HYPRE_LGMRESSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_LGMRESSetTol(solver, *TOL);       /* conv. tolerance */

     HYPRE_LGMRESSetKDim(solver,hypre_intpara[9]);     
     HYPRE_LGMRESSetAbsoluteTol(solver, hypre_dppara[5]);
     
     HYPRE_LGMRESSetAugDim(solver, hypre_intpara[10]);       
     
     /* Set some parameters (See Reference Manual for more parameters) */
     i = (verbosity >= 6);
     if(verbosity >= 10) i=3;
     HYPRE_LGMRESSetPrintLevel(solver, i);   /* print solve info */
     i = (verbosity >= 6);
     HYPRE_LGMRESSetLogging(solver, i);      /* needed to get run info later */

     /* Set the LGMRES preconditioner */
     if( hypre_pre == 1){
       HYPRE_LGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
     } else if(hypre_pre == 2 ) {
       HYPRE_LGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_AMSSolve,
				 (HYPRE_PtrToSolverFcn) HYPRE_AMSSetup, precond);
     }
     else if ( hypre_pre  == 3) {
       HYPRE_LGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ILUSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_ILUSetup, precond);
     } else if ( hypre_pre == 4) { 
       HYPRE_LGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
     }
       
     if (myverb > 12 ) fprintf(stdout,"SolveHypre: Setting up LGMRes linear system");
     HYPRE_ParCSRLGMRESSetup(solver, parcsr_A, par_b, par_x);

     
   } else if ( hypre_sol == 11) { /* COGMRES */
     /* Create solver */
     HYPRE_ParCSRCOGMRESCreate(comm, &solver);

     /* Set some parameters (See Reference Manual for more parameters) */
     //HYPRE_COGMRESSetStopCrit(solver, 0);     /* use the two norm as the stopping criteria */
     HYPRE_COGMRESSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_COGMRESSetTol(solver, *TOL);       /* conv. tolerance */

     HYPRE_COGMRESSetKDim(solver,hypre_intpara[9]);
     HYPRE_COGMRESSetAbsoluteTol(solver, hypre_dppara[5]);       

     HYPRE_COGMRESSetUnroll(solver, hypre_intpara[10]);       
     HYPRE_COGMRESSetCGS(solver, hypre_intpara[11]);       


     
     i = (verbosity >= 6);
     if(verbosity >= 10) i=3;
     HYPRE_COGMRESSetPrintLevel(solver, i);   /* print solve info */
     HYPRE_COGMRESSetLogging(solver,verbosity>=6);      /* needed to get run info later */

     /* Set the COGMRES preconditioner */
     if(hypre_pre == 1) {
       HYPRE_COGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
     } else if(hypre_pre == 2 ) {
       HYPRE_COGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_AMSSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_AMSSetup, precond);
     } else if ( hypre_pre == 3) {
       HYPRE_COGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ILUSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_ILUSetup, precond);
     } else if (hypre_pre == 4) { 
       HYPRE_COGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
				(HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
     }
       
     if (myverb > 12 ) fprintf(stdout,"SolveHypre: Setting up COGMRes linear system");
     HYPRE_ParCSRCOGMRESSetup(solver, parcsr_A, par_b, par_x);

     
   } else {
     fprintf( stdout,"SolveHypre: Hypre solver method %d not implemented!\n",hypre_sol);
     exit(EXIT_FAILURE);
   }

   Container->ilower = ilower;
   Container->iupper = iupper;     
   Container->hypre_method = *hypre_method;
   Container->A = A;
   Container->b = b;
   Container->x = x;
   Container->Atilde = Atilde;
   Container->solver = solver;
   Container->precond = precond;
   
   if (myverb > 5) fprintf( stdout, "SolveHypre: setup time (method %d): %g\n", 
			    Container->hypre_method, realtime_()-st );
   
} /* SolveHypre1 - matrix conversion and solver setup */


/* Update the stopping tolerance of a previously constructed solver */
void STDCALLBULL FC_FUNC(updatehypre,UPDATEHYPRE)
  ( double *TOL,  int *hypre_method,  int** ContainerPtr, int *verbosityPtr, int *fcomm )
{
  HYPRE_Solver solver;
  ElmerHypreContainer *Container;
  int hypre_sol;
  
  Container = (ElmerHypreContainer*)(*ContainerPtr);
  solver = Container->solver;

  hypre_sol = *hypre_method / 100;

  if ( hypre_sol == 1 ) { /* boomer AMG */
    HYPRE_BoomerAMGSetTol(solver, *TOL);
  }
  else if ( hypre_sol == 2 ) { /* AMS */
    HYPRE_AMSSetTol(solver, *TOL);
  }
  else if ( hypre_sol == 6) { /* CG */
    HYPRE_PCGSetTol(solver, *TOL);
  }
  else if ( hypre_sol == 7) { /* BiGSTAB method */
    HYPRE_BiCGSTABSetTol(solver, *TOL);
  }
  else if ( hypre_sol == 8) { /* GMRES */
    HYPRE_GMRESSetTol(solver, *TOL);
  }
  else if ( hypre_sol == 9) { /* FlexGMRes */
    HYPRE_FlexGMRESSetTol(solver, *TOL);
  }
  else if ( hypre_sol == 10) { /* LGMRes */
    HYPRE_LGMRESSetTol(solver, *TOL);
  }
  if ( hypre_sol == 11) { /* COGMRES method */
    HYPRE_COGMRESSetTol(solver, *TOL);
  }
  else {
    fprintf( stdout,"SolveHypre: Hypre solver method not implemented\n");
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

   int verbosity = *verbosityPtr, myverb;
   int hypre_sol, hypre_pre;

   int num_iterations;
   double final_res_norm;
   
   Container = (ElmerHypreContainer*)(*ContainerPtr);

   /* which process number am I? */
   MPI_Comm_rank(comm, &myid); 

   if(0) myid = 0;
   
   if( myid == 0 ) 
     myverb = verbosity;
   else
     myverb = 0;

   if (myverb > 10) fprintf(stdout,"SolveHypre: solving the linear system\n");
   if (Container==NULL) {
     if(myverb) fprintf( stdout, "ID no. %i: pointer passed into SolveHypre2 is NULL, not solving",myid);
     return;
   }
  
   st = realtime_();

   HYPRE_IJMatrixGetObject(Container->A, (void**) &parcsr_A);
   solver = Container->solver;
   precond = Container->precond;

   ilower = Container->ilower;
   iupper = Container->iupper;
   local_size = *nrows;

   hypre_sol = Container->hypre_method / 100;
   hypre_pre = Container->hypre_method % 100;


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

   HYPRE_IJVectorGetObject(b, (void **) &par_b);
   
   HYPRE_IJVectorAssemble(x);
   HYPRE_IJVectorGetObject(x, (void **) &par_x);
   
   /* Now setup and solve! */
   if( hypre_sol == 1) {
     if(myverb > 6) fprintf(stdout,"SolveHypre: Solving linear system with BoomerAMG (method %d)\n",Container->hypre_method);
     HYPRE_BoomerAMGSolve(Container->solver, parcsr_A, par_b, par_x);     
     if (myverb > 5 ) {
       HYPRE_BoomerAMGGetNumIterations(Container->solver, &num_iterations);
       HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
     }
   }
   
   else if( hypre_sol == 2 ) {
     if(myverb > 6) fprintf(stdout,"SolveHypre: Solving linear system with AMS (method %d)\n",Container->hypre_method);
     HYPRE_AMSSolve(Container->solver, parcsr_A, par_b, par_x);
     if (myverb > 5 ) {
       HYPRE_AMSGetNumIterations(Container->solver, &num_iterations);
       HYPRE_AMSGetFinalRelativeResidualNorm(solver, &final_res_norm);
     }
   }
   
   else if ( hypre_sol == 6) {
     if(myverb > 6) fprintf(stdout,"SolveHypre: Solving linear system with PCG (method %d)\n",Container->hypre_method);
//     HYPRE_ParCSRPCGSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_ParCSRPCGSolve(Container->solver, parcsr_A, par_b, par_x);

     if (myverb > 5 ) {
       HYPRE_PCGGetNumIterations(Container->solver, &num_iterations);
       HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
     }
   }

   if ( hypre_sol == 7) {
     if(myverb > 6) fprintf(stdout,"SolveHypre: Solving linear system with BiCGStab (method %d)\n",Container->hypre_method);
     //     HYPRE_ParCSRBiCGSTABSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_ParCSRBiCGSTABSolve(Container->solver, parcsr_A, par_b, par_x);

     if (myverb > 5 ) {
       HYPRE_BiCGSTABGetNumIterations(Container->solver, &num_iterations);
       HYPRE_BiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
     }
   }

   else if ( hypre_sol == 8) {
     if(myverb > 6) fprintf(stdout,"SolveHypre: Solving linear system with GMRes (method %d)\n",Container->hypre_method);
//     HYPRE_GMRESSetMaxIter(solver, *Rounds); /* max GMRES iterations */
     HYPRE_ParCSRGMRESSolve(Container->solver, parcsr_A, par_b, par_x);

     if (myverb > 5 ) {
       HYPRE_GMRESGetNumIterations(Container->solver, &num_iterations);
       HYPRE_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
     }
   }

   else if ( hypre_sol == 9) {
     if(myverb > 6) fprintf(stdout,"SolveHypre: Solving linear system with FlexGMRes (method %d)\n",Container->hypre_method);
//     HYPRE_ParCSRFlexGMRESSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_ParCSRFlexGMRESSolve(Container->solver, parcsr_A, par_b, par_x);
     
     if (myverb > 5 ) {
       HYPRE_FlexGMRESGetNumIterations(Container->solver, &num_iterations);
       HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
     }
   }

   else if ( hypre_sol == 10) {
     if(myverb > 6) fprintf(stdout,"SolveHypre: Solving linear system with LGMRes (method %d)\n",Container->hypre_method);     
//     HYPRE_ParCSRLGMRESSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_ParCSRLGMRESSolve(Container->solver, parcsr_A, par_b, par_x);

     if (myverb > 5 ) {
       HYPRE_LGMRESGetNumIterations(Container->solver, &num_iterations);
       HYPRE_LGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
     }
   }
   else if ( hypre_sol == 11) {
     if(myverb > 6) fprintf(stdout,"SolveHypre: Solving linear system with COGMRes (method %d)\n",Container->hypre_method);     
//     HYPRE_ParCSRLGMRESSetMaxIter(solver, *Rounds); /* max iterations */
     HYPRE_ParCSRCOGMRESSolve(Container->solver, parcsr_A, par_b, par_x);
          
     if (myverb > 5 ) {
       HYPRE_COGMRESGetNumIterations(Container->solver, &num_iterations);
       HYPRE_COGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
     }
   }
   
   for( k=0,i=0; i<local_size; i++ )
     if ( owner[i] ) rcols[k++] = globaldofs[i];
   
   HYPRE_IJVectorGetValues(x, k, rcols, txvec );
   
   for( i=0,k=0; i<local_size; i++ )
     if ( owner[i] ) xvec[i] = txvec[k++];

   if(myverb > 5) fprintf(stdout,"SolveHypre: Required iterations %d (method %d) to norm %lg\n",
			  num_iterations,Container->hypre_method, final_res_norm);   
   if (myverb > 4) fprintf( stdout, "SolveHypre: Solution time (method %d): %g\n", 
			    Container->hypre_method, realtime_()-st );
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
void STDCALLBULL FC_FUNC(solvehypre4,SOLVEHYPRE4)
  ( int** ContainerPtr, int *verbosityPtr ) {
  
   ElmerHypreContainer* Container = (ElmerHypreContainer*)(*ContainerPtr);
   int verbosity = *verbosityPtr, myverb;

   int hypre_sol, hypre_pre, myid;

   if (Container==0) return;

   myverb = verbosity;
   
   hypre_sol = Container->hypre_method / 100;
   hypre_pre = Container->hypre_method % 100;

   
   if(myverb > 10 ) fprintf(stdout,"SolveHypre: Detroying Hypre solver structures!\n");

   /* Destroy Hypre preconditioner */
   if ( hypre_pre == 1 ) {
     if(myverb > 10) fprintf(stdout,"SolveHypre: Detroying BoomerAMG preconditioner\n");
     HYPRE_BoomerAMGDestroy(Container->precond);
   } 
   else if ( hypre_pre == 2 ) {
     if(myverb > 10) fprintf(stdout,"SolveHypre: Detroying AMS preconditioner\n");
     HYPRE_AMSDestroy(Container->precond);
     if (Container->G) {
       // This leas to core dump...
       // HYPRE_IJMatrixDestroy(Container->G);
     }
   }
   else if ( hypre_pre == 3 ) {
     HYPRE_ILUDestroy(Container->precond);
   } else if ( hypre_pre == 4 ) {
     HYPRE_ParaSailsDestroy(Container->precond);
   }

   
   /* Destroy Hypre solver */
   if ( hypre_sol == 1) { /* boomer AMG */
     if(myverb > 10) fprintf(stdout,"SolveHypre: Detroying BoomerAMG solver\n");
     HYPRE_BoomerAMGDestroy(Container->solver);
   }
   else if ( hypre_sol == 2 ) {
     if(myverb > 10) fprintf(stdout,"SolveHypre: Detroying AMS solver\n");
     HYPRE_AMSDestroy(Container->solver);
     if (Container->G) {
       // This leas to core dump...
       // HYPRE_IJMatrixDestroy(Container->G);
     }
   }
   else if ( hypre_sol == 6 ) {
     HYPRE_ParCSRPCGDestroy(Container->solver);
   }
   else if ( hypre_sol == 7) {
     HYPRE_ParCSRBiCGSTABDestroy(Container->solver);
   }
   else if ( hypre_sol == 8 ) {
     HYPRE_ParCSRGMRESDestroy(Container->solver);
   }
   else if ( hypre_sol == 9 ) {
     HYPRE_ParCSRFlexGMRESDestroy(Container->solver);
   }
   else if ( hypre_sol == 10 ) {
     HYPRE_ParCSRLGMRESDestroy(Container->solver);
   }
   else if ( hypre_sol == 11 ) {
     HYPRE_ParCSRCOGMRESDestroy(Container->solver);
   }

   if (Container->Atilde != Container->A) {
     HYPRE_IJMatrixDestroy(Container->Atilde);
   }

   if (Container->A) {
     HYPRE_IJMatrixDestroy(Container->A);
   }
   if(Container->x) HYPRE_IJVectorDestroy(Container->x);
   if(Container->b) HYPRE_IJVectorDestroy(Container->b);

   free(Container);
   *ContainerPtr = NULL;
}



void STDCALLBULL FC_FUNC(createhypreams,CREATEHYPREAMS)
 (
  int *nrows,int *rows, int *cols, double *vals, int *nnodes,
  int *grows, int *gcols, double *gvals,
  int *pirows, int *picols, double *pivals, 
  int *perm, int *invperm, int *globaldofs, int *owner,  int *globalnodes, 
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
   ElmerHypreContainer *Container;
   HYPRE_ParCSRMatrix parcsr_A,parcsr_G, parcsr_Pi;
   HYPRE_IJMatrix A,G, Pi;
   HYPRE_IJVector b;
   HYPRE_ParVector par_b;
   HYPRE_IJVector x;
   HYPRE_ParVector par_x;
   HYPRE_IJVector xx,yy,zz;
   HYPRE_ParVector par_xx,par_yy,par_zz;

   HYPRE_Solver solver, precond;
   int verbosity = *verbosityPtr, myverb;
   MPI_Comm comm=MPI_Comm_f2c(*fcomm);

   
   Container = (ElmerHypreContainer*)malloc(sizeof(ElmerHypreContainer));   
   *ContainerPtr=(int*)(Container);
   
   st  = realtime_();
   MPI_Comm_rank(comm, &myid);

   if(myid==0)
     myverb = verbosity;
   else
     myverb = 0;
   
   if(myverb>6) fprintf(stdout,"SolveHypre: Setting up AMS preconditioner / solver\n");
   
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

#if 0
   nlower=1000000000;
   nupper=0;
   for( i=0; i<local_nodes; i++ )
     {
       if ( nodeowner[i] ) {
	 if ( nupper < globalnodes[i] ) nupper=globalnodes[i];
	 if ( nlower > globalnodes[i] ) nlower=globalnodes[i];
       }
     }
#else
   nlower=1000000000;
   nupper=0;
   for( i=0; i<local_nodes; i++ )
       if  (nodeowner[i])  {
         k = 3*globalnodes[i]+2;
         if ( nupper < k ) nupper = k;
         k = 3*globalnodes[i];
         if ( nlower > k ) nlower = k;
       }
#endif
//   fprintf( stderr, "%d %d %d %d\n", ilower, iupper, nlower, nupper );

   HYPRE_IJMatrixCreate(comm, ilower, iupper, nlower, nupper, &G);
   HYPRE_IJMatrixSetObjectType(G, HYPRE_PARCSR);
   HYPRE_IJMatrixInitialize(G);
   
   {
      int nnz,irow,i,j,k,l,p,q,*rcols,csize=32;

      rcols = (int *)malloc( csize*sizeof(int) );
      for (i = 0; i < local_size; i++)
      {
         if( !owner[i] ) continue;
         nnz = grows[i+1] - grows[i];
         if ( nnz>csize ) {
           rcols = (int *)realloc( rcols, nnz*sizeof(int) );
           csize = nnz;
         }
         irow = globaldofs[i];
         for( k=0,j=grows[i]; j<grows[i+1]; j++,k++)
         {
           l = gcols[j-1]-1;
           p = l % 3;
           q = l / 3;
           rcols[k] = 3*globalnodes[q]+p;
         }
         HYPRE_IJMatrixAddToValues(G, 1, &nnz, &irow, rcols, &gvals[grows[i]-1]);
      }
      free( rcols );
   }
   
   HYPRE_IJMatrixAssemble(G);
   HYPRE_IJMatrixGetObject(G, (void**) &parcsr_G);


   nlower=1000000000;
   nupper=0;
   for( i=0; i<local_nodes; i++ )
     {
       if  (nodeowner[i])  {
         k = 6*globalnodes[i]+5;
         if ( nupper < k ) nupper = k;
         k = 6*globalnodes[i];
         if ( nlower > k ) nlower = k;
       }
     }
//   fprintf( stderr, "%d %d %d %d\n", ilower, iupper, nlower, nupper );
   HYPRE_IJMatrixCreate(comm, ilower, iupper, nlower, nupper, &Pi);
   HYPRE_IJMatrixSetObjectType(Pi, HYPRE_PARCSR);
   HYPRE_IJMatrixInitialize(Pi);
   
   {
      int nnz,irow,i,j,k,l,p,q,*rcols,csize=32;

      rcols = (int *)malloc( csize*sizeof(int) );
      for (i = 0; i < local_size; i++)
      {
         if( !owner[i] ) continue;
         nnz =  pirows[i+1] - pirows[i];
         if ( nnz>csize ) {
           rcols = (int *)realloc( rcols, nnz*sizeof(int) );
           csize = nnz;
         }
         irow = globaldofs[i];
         for( k=0,j=pirows[i]; j<pirows[i+1]; j++,k++)
         {
           l = picols[j-1]-1;
           p = l % 6;
           q = l / 6;
           rcols[k] = 6*globalnodes[q]+p;
         }
         HYPRE_IJMatrixAddToValues(Pi, 1, &nnz, &irow, rcols, &pivals[pirows[i]-1]);
      }
      free( rcols );
   }
   
   HYPRE_IJMatrixAssemble(Pi);
   HYPRE_IJMatrixGetObject(Pi, (void**) &parcsr_Pi);

   
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

   HYPRE_AMSCreate(&precond); 
   HYPRE_AMSSetDiscreteGradient(precond,parcsr_G);
   HYPRE_AMSSetInterpolations(precond, parcsr_Pi, NULL, NULL, NULL);
//   HYPRE_AMSSetEdgeConstantVectors(precond,par_xx,par_yy,par_zz);
//   HYPRE_AMSSetCoordinateVectors(precond,par_xx,par_yy,par_zz);

   // AMS Parameters
   HYPRE_AMSSetMaxIter(precond,hypre_intpara[0]);
   HYPRE_AMSSetTol(precond,hypre_dppara[0]);
   
   HYPRE_AMSSetCycleType(precond, hypre_intpara[1]); // 1-14
   HYPRE_AMSSetSmoothingOptions(precond, hypre_intpara[2], hypre_intpara[3],
				hypre_dppara[1], hypre_dppara[2]);
   HYPRE_AMSSetAlphaAMGOptions(precond, 10, 1, 3, hypre_dppara[3], 0, 0);
   HYPRE_AMSSetBetaAMGOptions(precond, 10, 1, 3, hypre_dppara[4], 0, 0);

   if(hypre_intpara[4]) 
     HYPRE_AMSSetBetaPoissonMatrix(precond,NULL);
   
   i = (verbosity >= 6);
   if(verbosity >= 10) i=3;
   HYPRE_AMSSetPrintLevel(precond, i); 
   
   Container->precond = precond;
   Container->G = G; 
   Container->Pi = Pi; 
   
   if(myverb>5) fprintf( stdout, "SolveHypre: AMS preconditioner setup time: %g\n", realtime_()-st );
}

#endif
