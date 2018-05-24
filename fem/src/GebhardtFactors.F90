!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Jun 1997
! *
! *****************************************************************************/

!> \ingroup Programs
!> {
!> \defgroup GebhardtFactors Program GebhardtFactors
!> A stand-alone program that Solves for Gebhardt factors. 
!> This has been integrated into Elmer nowadays, maybe this one could go.

#define FULLMATRIX

   MODULE GebhardtFactorGlobals

#ifdef FULLMATRIX
      DOUBLE PRECISION, ALLOCATABLE :: GFactorMatrix(:,:)
#else
      USE CRSMatrix
      USE IterSolve

      TYPE(Matrix_t), POINTER :: GFactorMatrix
#endif
   END MODULE GebhardtFactorGlobals

   PROGRAM GebhardtFactors

     USE Types
     USE Lists

     USE GebhardtFactorGlobals

     USE CoordinateSystems
     USE SolverUtils
     USE MainUtils
     USE ModelDescription
     USE ElementDescription

     IMPLICIT NONE

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Model_t),POINTER :: Model
     TYPE(Mesh_t), POINTER :: Mesh

     INTEGER :: i,j,k,l,t,k1,k2,n,iter,Ndeg,Time,NSDOFs,istat

     DOUBLE PRECISION :: SimulationTime,dt

     TYPE(Element_t), POINTER :: CurrentElement
     TYPE(Factors_t), ALLOCATABLE :: ViewFactors(:)

     INTEGER :: BandSize,SubbandSize,RadiationSurfaces,Row,Col,MatrixElements

     DOUBLE PRECISION :: s,Norm,PrevNorm,Emissivity

     INTEGER, ALLOCATABLE, TARGET :: RowSpace(:),Reorder(:),ElementNumbers(:)

     DOUBLE PRECISION, POINTER :: Vals(:),Fac(:)
     INTEGER, POINTER :: Cols(:), NodeIndexes(:)

     CHARACTER(LEN=MAX_NAME_LEN) :: RadiationFlag,GebhardtFactorsFile, &
                            ViewFactorsFile,OutputName,ModelName,eq

     TYPE(Variable_t), POINTER :: var

     LOGICAL :: GotIt, Found, Radiation
     TYPE(Solver_t), POINTER :: Solver
     DOUBLE PRECISION, ALLOCATABLE :: SOL(:),RHS(:)

     EXTERNAL Matvec


     PRINT*, '----------------------------------------------------'
     PRINT*,' ELMER Gebhardt factor computation program, welcome'
     PRINT*, '----------------------------------------------------'
!------------------------------------------------------------------------------
!    Read element definition file, and initialize element types
!------------------------------------------------------------------------------
     CALL InitializeElementDescriptions
!------------------------------------------------------------------------------
!    Read Model from Elmer Data Base
!------------------------------------------------------------------------------
     PRINT*, ' '
     PRINT*, ' '
     PRINT*, '-----------------------'
     PRINT*,'Reading Model ...       '

!------------------------------------------------------------------------------
                                                                                                              
     OPEN( 1,file='ELMERSOLVER_STARTINFO', STATUS='OLD', ERR=10 )
     GOTO 20
                                                                                                              
10   CONTINUE
                                                                                                              
     CALL Fatal( 'GebhardtFactors', 'Unable to find ELMERSOLVER_STARTINFO, cant execute.' )
                                                                                                              
20   CONTINUE
       READ(1,'(a)') ModelName
     CLOSE(1)

     Model => LoadModel( ModelName,.FALSE. )
     CurrentModel => Model

     NULLIFY( Mesh )
     DO i=1,Model % NumberOfSolvers
       Solver => Model % Solvers(i)
       Radiation = ListGetLogical( Solver % Values, 'Radiation Solver', Found )
       eq = ListGetString( Solver % Values, 'Equation' )
       IF ( Radiation .OR. TRIM(eq) == 'heat equation' ) THEN
         Mesh => Solver % Mesh
         EXIT
       ENDIF
     END DO
     IF ( .NOT.ASSOCIATED( Mesh ) ) THEN
       PRINT*,'ERROR: GebhardtFactors: No heat equation definition. ' // &
                  'Cannot compute factors.'
       STOP
     END IF

     PRINT*,'... Done               '
     PRINT*, '-----------------------'

!------------------------------------------------------------------------------
!    Add coordinates to list of variables so that coordinate dependent
!    parameter computing routines can ask for them...
!------------------------------------------------------------------------------
    CALL VariableAdd(Mesh % Variables,Mesh,Solver, &
          'Coordinate 1',1,Mesh % Nodes % x )

    CALL VariableAdd(Mesh % Variables,Mesh,Solver, &
          'Coordinate 2',1,Mesh % Nodes % y )

    CALL VariableAdd(Mesh % Variables,Mesh,Solver, &
          'Coordinate 3',1,Mesh % Nodes % z )

    CALL SetCurrentMesh( Model, Mesh )

!------------------------------------------------------------------------------
!    Here we start...
!------------------------------------------------------------------------------
     RadiationSurfaces = 0
!------------------------------------------------------------------------------
!    loop to get the surfaces participating in radiation, discard the rest
!    of the elements...
!------------------------------------------------------------------------------
     ALLOCATE( ElementNumbers(Model % NumberOfBoundaryElements) )
     DO t=Model % NumberOfBulkElements+1, &
            Model % NumberOfBulkElements + Model % NumberOfBoundaryElements
       CurrentElement => Model % Elements(t)
       k = CurrentElement % BoundaryInfo % Constraint

       IF ( CurrentElement % TYPE % ElementCode /= 101 ) THEN
         DO i=1,Model % NumberOfBCs
           IF ( Model % BCs(i) % Tag == k ) THEN
             RadiationFlag = ListGetString( Model % BCs(i) % Values, &
                            'Radiation', GotIt )
             IF ( RadiationFlag == 'diffuse gray' ) THEN
               RadiationSurfaces = RadiationSurfaces + 1
               ElementNumbers(RadiationSurfaces) = t
               Model % Elements(RadiationSurfaces) = Model % Elements(t)
             END IF
           END IF
         END DO
       END IF
     END DO

     IF ( RadiationSurfaces == 0 ) THEN
       PRINT*,'WARNING: Gebhardtfactors: No surfaces participating in radiation?'
       PRINT*,'WARNING: Gebhardtfactors: Stopping cause nothing to be done...'
       STOP
     END IF

     Model % NumberOfBulkElements = 0
     Model % NumberOfBoundaryElements = RadiationSurfaces

     ALLOCATE( RowSpace( RadiationSurfaces ), Reorder( RadiationSurfaces ), &
                     STAT=istat )
     IF ( istat /= 0 ) THEN
       PRINT*,'ERROR: Gebhardtfactors: Memory allocation error. Aborting'
       STOP
     END IF

     ViewFactorsFile = ListGetString(Model % Simulation,'View Factors',GotIt)

     IF ( .NOT.GotIt ) THEN
       ViewFactorsFile = 'ViewFactors.dat'
     END IF

     IF ( LEN_TRIM(Model % Mesh % Name) > 0 ) THEN
       OutputName = TRIM(OutputPath) // '/' // TRIM(Model % Mesh % Name) // '/' // TRIM(ViewFactorsFile)
     ELSE
       OutputName = TRIM(ViewFactorsFile)
     END IF

     OPEN( 1,File=TRIM(OutputName) )

     ALLOCATE( ViewFactors(RadiationSurfaces) )
     MatrixElements = 0
     RowSpace       = 0
     DO i=1,RadiationSurfaces
       READ( 1,* ) RowSpace(i)

       ALLOCATE( ViewFactors(i) % Elements(RowSpace(i)) )
       ALLOCATE( ViewFactors(i) % Factors(RowSpace(i)) )

       ViewFactors(i) % NumberOfFactors = RowSpace(i)
       Vals => ViewFactors(i) % Factors
       Cols => ViewFactors(i) % Elements

       DO j=1,RowSpace(i)
         READ(1,*) t,Cols(j),Vals(j)
         MatrixElements = MatrixElements + 1
       END DO

#ifndef FULLMATRIX
       Reorder(i) = i
#endif
     END DO

     DO i=1,RadiationSurfaces

       Cols => ViewFactors(i) % Elements

       k = 0
       DO j=1,ViewFactors(i) % NumberOfFactors
         IF ( Cols(j) == i ) THEN
           k = 1
           EXIT
         END IF
       END DO
       IF ( k == 0 ) THEN
         RowSpace(i) = RowSpace(i) + 1
         MatrixElements = MatrixElements + 1
       END IF

     END DO

     CLOSE(1)

     ALLOCATE( RHS(RadiationSurfaces), SOL(RadiationSurfaces),Fac(RadiationSurfaces),STAT=istat )
     IF ( istat /= 0 ) THEN
       PRINT*,'ERROR: Gebhardtfactors: Memory allocation error. Aborting'
       STOP
     END IF
!
!    The coefficient matrix 
!
     PRINT*,' '
     PRINT*,'Computing factors...'

#ifdef FULLMATRIX
     ALLOCATE(GFactorMatrix(RadiationSurfaces,RadiationSurfaces),STAT=istat)

     IF ( istat /= 0 ) THEN
       PRINT*,'Gebhardt factors: Memory allocation error. Aborting'
       STOP
     END IF
     GFactorMatrix = 0.0D0
#else

     GFactorMatrix => CRS_CreateMatrix( RadiationSurfaces, &
          MatrixElements,RowSpace,1,Reorder,.TRUE. )

     DO t=1,RadiationSurfaces
       Cols => ViewFactors(t) % Elements

       DO j=1,ViewFactors(t) % NumberOfFactors
         CALL CRS_MakeMatrixIndex( GFactorMatrix,t,Cols(j) )
       END DO

       CALL CRS_MakeMatrixIndex( GFactorMatrix,t,t )
     END DO

     CALL CRS_SortMatrix( GFactorMatrix )
     CALL CRS_ZeroMatrix( GFactorMatrix )
#endif
!
!    Fill the matrix for gebhardt factors
!
     DO t=1,RadiationSurfaces
       CurrentElement => Model % Elements(t)
       n = CurrentElement % TYPE % NumberOfNodes
       k = CurrentElement % BoundaryInfo % Constraint
       NodeIndexes => CurrentElement % NodeIndexes

       DO i=1,Model % NumberOfBCs
         IF ( Model % BCs(i) % Tag  == k ) THEN
!          Emissivity = ListGetConstReal( Model % BCs(i) % Values, &
!                        'Emissivity' )
           Emissivity = SUM( ListGetReal( Model % BCs(i) % Values, &
                    'Emissivity', n, NodeIndexes ) ) / n

           Vals => ViewFactors(t) % Factors
           Cols => ViewFactors(t) % Elements

           DO j=1,ViewFactors(t) % NumberOfFactors
#ifdef FULLMATRIX
             GFactorMatrix(t,Cols(j)) = GFactorMatrix(t,Cols(j)) - &
                     (1-Emissivity)*Vals(j)
#else
             CALL CRS_AddToMatrixElement( GFactorMatrix,t,Cols(j), &
                    -(1-Emissivity)*Vals(j) )
#endif
           END DO
#ifdef FULLMATRIX
           GFactorMatrix(t,t) = GFactorMatrix(t,t) + 1.0D0
#else
           CALL CRS_AddToMatrixElement( GFactorMatrix,t,t,1.0D0 )
#endif
           EXIT
         END IF
       END DO
     END DO

     ALLOCATE( Solver )
     Solver % Values => ListAllocate()

     CALL ListAddString( Solver % Values, &
                  'Linear System Iterative Method', 'CGS' )

     CALL ListAddInteger( Solver % Values, &
                      'Linear System Max Iterations', 100 )

     CALL ListAddConstReal( Solver % Values, &
            'Linear System Convergence Tolerance', 1.0D-9 )

     CALL ListAddString( Solver % Values, &
                  'Linear System Preconditioning', 'None' )

     CALL ListAddInteger( Solver % Values, &
                       'Linear System Residual Output', 0 )

     RHS = 0.0D0
     SOL = 1.0D-4

     GebhardtFactorsFile = ListGetString(Model % Simulation, &
                   'Gebhardt Factors',GotIt )

     IF ( .NOT.GotIt ) THEN
       GebhardtFactorsFile = 'GebhardtFactors.dat'
     END IF


     IF ( LEN_TRIM(Model % Mesh % Name) > 0 ) THEN
       OutputName = TRIM(OutputPath) // '/' // TRIM(Model % Mesh % Name) // '/' // TRIM(GebhardtFactorsFile)
     ELSE
       OutputName = TRIM(GebhardtFactorsFile)
     END IF

     OPEN( 1,File=TRIM(OutputName) )

     WRITE( 1,* ) RadiationSurfaces

     DO t=1,RadiationSurfaces
       WRITE(1,*) t,ElementNumbers(t)
     END DO

     DO t=1,RadiationSurfaces

       RHS    = 0.0D0
       RHS(t) = 1.0D0

#ifdef FULLMATRIX
       CALL FIterSolver( RadiationSurfaces,SOL,RHS,Solver )
#else
       CALL IterSolver( GFactorMatrix,SOL,RHS,Solver )
#endif

       CurrentElement => Model % Elements(t)
       n = CurrentElement % TYPE % NumberOfNodes
       NodeIndexes => CurrentElement % NodeIndexes

       k = CurrentElement % BoundaryInfo % Constraint
       DO j=1,Model % NumberOfBCs
         IF ( Model % BCs(j) % Tag == k ) THEN
!          Emissivity = ListGetConstReal( Model % BCs(j) % Values, &
!                           'Emissivity' )
           Emissivity = SUM(ListGetReal( Model % BCs(j) % Values, &
               'Emissivity', n, NodeIndexes)) / n
           EXIT
         END IF
       END DO

       n = 0
       DO i=1,RadiationSurfaces
         Vals => ViewFactors(i) % Factors
         Cols => ViewFactors(i) % Elements

         s = 0.0d0
         DO k=1,ViewFactors(i) % NumberOfFactors
           s = s + Vals(k) * SOL(Cols(k))
         END DO
         Fac(i) = Emissivity * s
         IF ( Fac(i) > 1.0d-15 ) n = n + 1
       END DO
       
#if 0
       WRITE( 1,* ) n
       DO i=1,RadiationSurfaces
         IF ( Fac(i) > 1.0d-15 ) WRITE( 1,* ) t,i,Fac(i)
       END DO
#else
       WRITE( 1,* ) RadiationSurfaces
       DO i=1,RadiationSurfaces
         WRITE( 1,* ) t,i,Fac(i)
       END DO
#endif
     END DO

     ! Free allocated temporaries:
     ! ---------------------------
     DO i=1,RadiationSurfaces
        DEALLOCATE( ViewFactors(i) % Factors, ViewFactors(i) % Elements )
     END DO
     DEALLOCATE( ViewFactors )

     CLOSE(1)
     PRINT*,'...Done'
     PRINT*,' '
     PRINT*,'*** GebhardtFactors: ALL DONE ***'
!    CALL FLUSH(6)

#ifdef FULLMATRIX
   CONTAINS

#include "huti_fdefs.h"
      SUBROUTINE FIterSolver( N,x,b,SolverParam )

#ifdef USE_ISO_C_BINDINGS
        USE huti_sfe
#endif
        IMPLICIT NONE

        TYPE(Solver_t) :: SolverParam

        DOUBLE PRECISION, DIMENSION(:) CONTIG :: x,b

        DOUBLE PRECISION :: dpar(50)

        INTEGER :: ipar(50),wsize,N
        DOUBLE PRECISION, ALLOCATABLE :: work(:,:)
#ifndef USE_ISO_C_BINDINGS
        INTEGER  :: HUTI_D_CGS
        EXTERNAL :: HUTI_D_CGS
        INTEGER(KIND=AddrInt) :: AddrFunc
#else
        INTEGER(KIND=AddrInt) :: AddrFunc
        EXTERNAL :: AddrFunc
#endif
        INTEGER(KIND=addrInt) :: iterProc, mvProc, dProc=0

        ipar = 0
        dpar = 0.0D0

        HUTI_WRKDIM = HUTI_CGS_WORKSIZE
        wsize = HUTI_WRKDIM
          
        HUTI_NDIM     = N
        HUTI_DBUGLVL  = 0
        HUTI_MAXIT    = 100
 
        ALLOCATE( work(wsize,N) )

        IF ( ALL(x == 0.0) ) THEN
          HUTI_INITIALX = HUTI_RANDOMX
        ELSE
          HUTI_INITIALX = HUTI_USERSUPPLIEDX
        END IF

        HUTI_TOLERANCE = 1.0d-10

        HUTI_MAXTOLERANCE = 1.0d20
        
        mvProc   = AddrFunc( Matvec )
        iterProc = AddrFunc(HUTI_D_CGS)
        CALL iterCall( iterProc,x,b,ipar,dpar,work,mvProc, &
               dProc, dProc, dProc, dProc, dProc )
          
        DEALLOCATE( work )

      END SUBROUTINE FIterSolver 
#endif

  END PROGRAM GebhardtFactors


#ifdef FULLMATRIX
  SUBROUTINE Matvec( u,v,ipar )

     USE GebhardtFactorGlobals

     DOUBLE PRECISION :: u(*),v(*)
     INTEGER :: ipar(*)

     INTEGER :: i,j,n
     DOUBLE PRECISION :: s

     n = HUTI_NDIM

     DO i=1,n
       s = 0.0D0
       DO j=1,n
         s = s + GFactorMatrix(i,j)*u(j)
       END DO
       v(i) = s
     END DO

  END SUBROUTINE Matvec
#endif

!> \}




