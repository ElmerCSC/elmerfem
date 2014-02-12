!/*****************************************************************************
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
! * FreeSurfaceSolver: Solver for free surface evolution in 2d and 3d flows
! *                    with/without surface flux 
! *
! ******************************************************************************
! *
! *                    Author:  Juha Ruokolainen
! *
! *                    Address: CSC - IT Center for Science Ltd.
! *                                Keilaranta 14, P.O. BOX 405
! *                                  02101 Espoo, Finland
! *                                  Tel. +358 0 457 2723
! *                                Telefax: +358 0 457 2183
! *                              EMail: Thomas.Zwinger@csc.fi
! *
! *                       Date: 17  May 2002
! *
! *                Modified by: Peter Råback, Juha Ruokolainen, Mikko Lyly
! *
! *
!/******************************************************************************
! *
! *       Modified by: Thomas Zwinger
! *
! *       Date of modification: 4. May 2004
! *
! *****************************************************************************/
   SUBROUTINE FSurfSolverPenalty( Model,Solver,dt,TransientSimulation )
  !DEC$ATTRIBUTES DLLEXPORT :: FSurfSolverPenalty
     USE DefUtils
     USE Differentials
     IMPLICIT NONE


!------------------------------------------------------------------------------
!    external variables
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     TYPE(Solver_t):: Solver
     REAL(KIND=dp) :: dt
     LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     INTEGER :: & 
          i,j,t,N,NMAX,MMAX,KMAX,Nmatrix,bf_id,DIM,istat,&
          NSDOFs,NonlinearIter, iter, SubroutineVisited=0


     REAL(KIND=dp) :: &
          at,st,totat,totst,CPUTime,Norm,PrevNorm,LocalBottom, cv, &
          Relax, MaxDisp, maxdh, NonlinearTol, UzawaParameter, RelativeChange, &
          SteadyTol,SteadyDisp

     LOGICAL ::&
          firstTime=.TRUE., GotIt, AllocationsDone = .FALSE., stat, &
          NeedOldValues,  Bubbles = .FALSE.,&
          NormalFlux = .TRUE., SubstantialSurface = .TRUE.,&
          UseBodyForce = .TRUE., LimitedSolution(2) = .FALSE., &
          Converged = .FALSE., &
          LimitDisp,DispSteady



     CHARACTER(LEN=MAX_NAME_LEN)  ::  FlowSolverName


     TYPE(Nodes_t)   :: ElementNodes
     TYPE(Element_t),POINTER :: Element
     TYPE(Variable_t), POINTER :: FlowSol


     INTEGER, POINTER ::&
          FreeSurfPerm(:), FlowPerm(:), NodeIndexes(:)


     INTEGER, ALLOCATABLE :: BoundaryMapping(:), BoundaryPerm(:)


     REAL(KIND=dp), POINTER ::&
          FreeSurf(:), PreFreeSurf(:,:), TimeForce(:),&
          FlowSolution(:), PrevFlowSol(:,:), ElemFreeSurf(:), OldFreeSurf(:)
 
     REAL(KIND=dp), ALLOCATABLE :: &          
          STIFF(:,:),SourceFunc(:),FORCE(:), &
          MASS(:,:), Velo(:,:), Flux(:,:),  Limit(:,:), &
          auxReal(:), OldDy(:), LocalLimit(:,:)
     INTEGER, ALLOCATABLE :: AlreadyDone(:)
          
     TYPE(ValueList_t), POINTER :: BodyForce, SolverParams, Material
!-----------------------------------------------------------------------------
!      remember these variables
!----------------------------------------------------------------------------- 
     SAVE STIFF, MASS, SourceFunc, FORCE, &
          ElementNodes, AllocationsDone, Velo, OldFreeSurf, TimeForce, &
          SubroutineVisited, ElemFreeSurf, Flux, SubstantialSurface, NormalFlux,&
          UseBodyForce,  Limit, BoundaryPerm, BoundaryMapping, &
          auxReal, Norm, PrevNorm, OldDy, LocalLimit, AlreadyDone,KMAX


!------------------------------------------------------------------------------
!    Get constants
!------------------------------------------------------------------------------
     MMAX = Solver % Mesh % NumberOfNodes
     DIM = CoordinateSystemDimension()
     SolverParams => GetSolverParams()
     
     MaxDisp = GetConstReal( SolverParams, 'Maximum Displacement', LimitDisp)
     NeedOldValues=LimitDisp


     cv = GetConstReal( SolverParams, 'Velocity Implicity', GotIt)
     IF(.NOT. GotIt) cv = 1.0d0 
     WRITE(Message,'(a,F8.2)') 'Velocity implicity (1=fully implicit)=', cv
     CALL Info('FreeSurfaceSolver', Message, Level=4)


     Relax = GetConstReal( SolverParams, 'Relaxation Factor', GotIt)
     IF(.NOT. GotIt) Relax = 1.0d0
     NeedOldValues = (NeedOldValues.OR.(GotIt .AND. (Relax < 1.0d00)))


     NonlinearIter = GetInteger(   SolverParams, &
          'Nonlinear System Max Iterations', GotIt )
     IF ( .NOT.GotIt ) NonlinearIter = 1


     NonlinearTol  = GetConstReal( SolverParams, &
          'Nonlinear System Convergence Tolerance',    GotIt )
     IF ((NonlinearIter > 1) .AND. (.NOT.GotIt)) THEN
        CALL Fatal('FreeSurfaceSolver',&
             'No Value for Nonlinear System Convergence Tolerance found')
     END IF

! pour limiter en deplacement
     SteadyTol = GetConstReal(   SolverParams, &
                    'Steady State Convergence Tolerance', GotIt )
     SteadyDisp = GetConstReal(   SolverParams, &
                    'Steady State Displacement Tolerance', DispSteady )
                                                                                                                             
!!!!!!!!!!!!!!
 


    
     UzawaParameter = GetConstReal( SolverParams, &
               'Uzawa Parameter', GotIt)
     IF (.NOT.GotIt) THEN                
        LimitedSolution = .FALSE.
     END IF


     Bubbles = GetLogical( SolverParams,'Bubbles',GotIt )
     IF (Bubbles) THEN
        CALL Info('FreeSurfaceSolver', 'Using Bubble stabilization', Level=4)
     ELSE
        CALL Info('FreeSurfaceSolver', &
             'Using residual squared-stabilized formulation.', Level=4)
     END IF


     UseBodyForce =  GetLogical( SolverParams,'Use Accumulation',GotIt )
     IF (.NOT.Gotit) UseBodyForce = .TRUE.


     IF (UseBodyForce) THEN
        NormalFlux =  GetLogical( SolverParams,'Normal Flux',GotIt )
        IF (.NOT.Gotit) NormalFlux = .TRUE.


        IF (NormalFlux) THEN
           CALL Info('FreeSurfaceSolver', &
                'Using scalar value for accumulaton/ablation rate', Level=4)
        ELSE
           CALL Info('FreeSurfaceSolver', &
                'Computing accumulaton/ablation rate from input vector', Level=4)
        END IF
     ELSE
        NormalFlux = .TRUE.
        CALL Info('FreeSurfaceSolver', 'Zero accumulation/ablation', Level=4)
     END IF


!------------------------------------------------------------------------
!    Get Velocity from (Navier)-Stokes-Solver or AIFlow Solver
!------------------------------------------------------------------------
     FlowSolverName = GetString( SolverParams, 'Flow Solver Name', GotIt )    
     IF (.NOT.Gotit) FlowSolverName = 'Flow Solution'
     FlowSol => VariableGet( Solver % Mesh % Variables, FlowSolverName )
     IF ( ASSOCIATED( FlowSol ) ) THEN
       FlowPerm     => FlowSol % Perm
       NSDOFs       =  FlowSol % DOFs
       FlowSolution => FlowSol % Values
       PrevFlowSol => FlowSol % PrevValues
       WRITE(Message,'(a,i1,a,i1)') 'DIM=', DIM, ', NSDOFs=', NSDOFs
       CALL Info( 'FreeSurfaceSolver', Message, level=4) 
     ELSE
       CALL Info('FreeSurfaceSolver', 'No variable for velocity associated.', Level=4)
     END IF




!------------------------------------------------------------------------------
!    Get variables for the solution
!------------------------------------------------------------------------------
     FreeSurf     => Solver % Variable % Values     ! Nodal values for free surface displacement
     FreeSurfPerm => Solver % Variable % Perm       ! Permutations for free surface displacement
     PreFreeSurf  => Solver % Variable % PrevValues ! Nodal values for free surface displacement
                                                    !                     from previous timestep
!------------------------------------------------------------------------------
!    assign norm
!------------------------------------------------------------------------------
     IF (AllocationsDone) THEN 
        Norm = Solver % Variable % Norm
     ELSE
        Norm = 0.0d00
     END IF



!------------------------------------------------------------------------------
!    Allocate some permanent storage, 
!    this is done first time only or if there have been changes in the mesh
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
        NMAX = Model % MaxElementNodes
        IF (Bubbles) THEN
           Nmatrix = 2*NMAX
        ELSE
           Nmatrix = NMAX
        END IF


        ! get node mapping for boundary body
        KMAX = COUNT(FreeSurfPerm .NE. 0)
       
        IF (AllocationsDone ) THEN
           DEALLOCATE(  ElementNodes % x,    &
                ElementNodes % y,    &
                ElementNodes % z,    &
                TimeForce,           &
                FORCE,    &
                STIFF, &
                MASS,  &
                Velo, &
                OldDy, &
                AlreadyDone, &
                LocalLimit, &
                Flux, &
                ElemFreeSurf,&
                SourceFunc,&
                Limit, &
                auxReal, &
                BoundaryMapping, &
                BoundaryPerm)
        END IF
       
        ALLOCATE( ElementNodes % x( NMAX ),    &
             ElementNodes % y( NMAX ),    &
             ElementNodes % z( NMAX ),    &
             TimeForce( NMAX ),           &
             FORCE( Nmatrix ),    &
             STIFF( Nmatrix, Nmatrix ), &
             MASS( Nmatrix, Nmatrix ),  &
             Velo( 3, NMAX ), &
             OldDy( NMAX ),&
             AlreadyDone(KMAX),&
             LocalLimit(2, NMAX), &
             Flux( 3, NMAX), &
             ElemFreeSurf( NMAX ),&
             SourceFunc( NMAX ), &
             Limit( 2, kMAX ), &
             auxReal( NMAX ), &
             BoundaryMapping( KMAX ), &
             BoundaryPerm( MMAX ), &
             STAT=istat )
        IF ( istat /= 0 ) THEN
           CALL Fatal('FreeSurfaceSolver','Memory allocation error, Aborting.')
        END IF


        ! previous values
        !----------------
        IF(NeedOldValues) THEN


           IF (AllocationsDone) DEALLOCATE(OldFreeSurf)


           ALLOCATE(OldFreeSurf(KMAX), STAT=istat)
           
           IF ( istat /= 0 ) THEN
              CALL Fatal('FreeSurfaceSolver','Memory allocation error, Aborting.')
           END IF
           DO i=1,KMAX
              OldFreeSurf(i) = FreeSurf(i)
           END DO
        END IF


       
        IF (.NOT. AllocationsDone) THEN
           CALL Info('FreeSurfaceSolver','Memory allocations done', Level=3)
           AllocationsDone = .TRUE.
        ELSE
           CALL Info('FreeSurfaceSolver','Memory allocation update done', Level=3)
        END IF


        ! get boundary mapping and permutation
        !-------------------------------------
        j=0
        BoundaryPerm(1:MMAX) = 0
        DO i=1,MMAX
           IF (FreeSurfPerm(i) > 0) THEN
              j = j+1
              BoundaryMapping(j) = FreeSurfPerm(i) 
              BoundaryPerm(i) = j
           END IF
        END DO
        IF (j .NE. KMAX) THEN
           WRITE(Message,'(A,I5,A,I5)') 'Number of mapped nodes ', j,&
                ' does not match number of non-zero permutations ', KMAX
           CALL FATAL('FreeSurfaceSolver',Message)
        END IF
     END IF
     

      If (NeedOldValues) Then
           DO i=1,KMAX
              OldFreeSurf(i) = FreeSurf(i)
           END DO
       End if


!------------------------------------------------------------------------------
!    Nonlinear iteration loop
!------------------------------------------------------------------------------
     Converged = .FALSE.
     DO iter=1,NonlinearIter
        WRITE(Message,'(A,I5,A,I5)') 'Nonlinear iteration no.',iter,&
             ' of max. ', NonlinearIter
        CALL INFO('FreeSurfaceSolver',Message,level=4)


!------------------------------------------------------------------------------
!       Get limit for Upper and Lower surface position
!------------------------------------------------------------------------------
        LimitedSolution = .FALSE.
        ! Get lower and Upper limit
        ! (loop all elements)
        !--------------------------        
        DO t=1,Solver % NumberOfActiveElements
           Element => GetActiveElement(t)
           n = GetElementNOFNodes()
           CALL GetElementNodes( ElementNodes )
           Material => GetMaterial()
           IF (.NOT.ASSOCIATED(Material)) THEN
              CALL INFO('FreeSurfaceSolver','No material found', Level=3)
           ELSE
              auxReal(1:n) = ListGetReal(Material, 'Min Flow Depth', n, Element % NodeIndexes, GotIt)
              IF (GotIt) THEN
                 DO i= 1,n
                    Limit(1,BoundaryPerm(Element % NodeIndexes(i))) = auxReal(i)
                 END DO
                 LimitedSolution(1) = .TRUE.
              END IF
              auxReal(1:n) = ListGetReal(Material,'Max Flow Depth',n,Element % NodeIndexes, GotIt)
              IF (GotIt) THEN
                 DO i= 1,n          
                    Limit(2,BoundaryPerm(Element % NodeIndexes(i))) = auxReal(i)
                 END DO
                 LimitedSolution(2) = .TRUE.
              END IF
           END IF
        END DO
        
!------------------------------------------------------------------------------
!       Do some additional initialization, and go for it
!------------------------------------------------------------------------------
        AlreadyDone = 0
        totat = 0.0d0
        totst = 0.0d0
        at = CPUTime()
        CALL Info( 'FreeSurfaceSolver', 'start assembly', Level=4 )
        CALL DefaultInitialize()
!------------------------------------------------------------------------------
!       Do the assembly
!------------------------------------------------------------------------------
        DO t=1,Solver % NumberOfActiveElements
           Element => GetActiveElement(t)
           n = GetElementNOFNodes()
           NodeIndexes => Element % NodeIndexes


           ! set coords of highest occuring dimension to zero (to get correct path element)
           !-------------------------------------------------------------------------------
           ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
           IF (DIM == 2) THEN
              ElementNodes % y(1:n) = 0.0
              ElementNodes % z(1:n) = 0.0
           ELSE IF (DIM == 3) THEN
              ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
              ElementNodes % z(1:n) = 0.0d0
           ELSE
              WRITE(Message,'(a,i1,a)')&
                   'It is not possible to compute free-surface problems in DIM=',&
                   DIM, ' dimensions. Aborting'
              CALL Fatal( 'FreeSurfaceSolver', Message) 
              STOP   
           END IF


           ! get velocity profile
           IF ( ASSOCIATED( FlowSol ) ) THEN
              DO i=1,n
                 j = NSDOFs*FlowPerm(NodeIndexes(i))


                 IF(TransientSimulation .AND. ABS(cv-1.0) > 0.001) THEN
                    IF((DIM == 2) .AND. (NSDOFs == 3)) THEN
                       Velo(1,i) = cv * FlowSolution( j-2 ) + (1-cv) * PrevFlowSol(j-2,1)
                       Velo(2,i) = cv * FlowSolution( j-1 ) + (1-cv) * PrevFlowSol(j-1,1)
                       Velo(3,i) = 0.0d0
                    ELSE IF ((DIM == 3) .AND. (NSDOFs == 4)) THEN
                       Velo(1,i) = cv * FlowSolution( j-3 ) + (1-cv) * PrevFlowSol(j-3,1)
                       Velo(2,i) = cv * FlowSolution( j-2 ) + (1-cv) * PrevFlowSol(j-2,1)
                       Velo(3,i) = cv * FlowSolution( j-1 ) + (1-cv) * PrevFlowSol(j-1,1)
                    ELSE IF ((CurrentCoordinateSystem() == CylindricSymmetric) &
                         .AND. (DIM == 2) .AND. (NSDOFs == 4)) THEN  
                       Velo(1,i) = cv * FlowSolution( j-3 ) + (1-cv) * PrevFlowSol(j-3,1)
                       Velo(2,i) = cv * FlowSolution( j-2 ) + (1-cv) * PrevFlowSol(j-2,1)
                       Velo(3,i) = cv * FlowSolution( j-1 ) + (1-cv) * PrevFlowSol(j-1,1)
                    ELSE
                       WRITE(Message,'(a,i1,a,i1,a)')&
                            'DIM=', DIM, ' NSDOFs=', NSDOFs, ' does not combine. Aborting'
                       CALL Fatal( 'FreeSurfaceSolver', Message)               
                    END IF
                 ELSE
                    IF((DIM == 2) .AND. (NSDOFs == 3)) THEN
                       Velo(1,i) = FlowSolution( j-2 ) 
                       Velo(2,i) = FlowSolution( j-1 ) 
                       Velo(3,i) = 0.0d0
                    ELSE IF ((DIM == 3) .AND. (NSDOFs == 4)) THEN
                       Velo(1,i) = FlowSolution( j-3 ) 
                       Velo(2,i) = FlowSolution( j-2 ) 
                       Velo(3,i) = FlowSolution( j-1 ) 
                    ELSE IF ((CurrentCoordinateSystem() == CylindricSymmetric) &
                         .AND. (DIM == 2) .AND. (NSDOFs == 4)) THEN
                       Velo(1,i) = FlowSolution( j-3 ) 
                       Velo(2,i) = FlowSolution( j-2 ) 
                       Velo(3,i) = FlowSolution( j-1 ) 
                    ELSE
                       WRITE(Message,'(a,i1,a,i1,a)')&
                            'DIM=', DIM, ' NSDOFs=', NSDOFs, ' does not combine. Aborting'
                       CALL Fatal( 'FreeSurfaceSolver', Message)
                    END IF
                 END IF
              END DO
           ELSE
              Velo=0.0d0          
           END IF
!------------------------------------------------------------------------------
!          get the accumulation/ablation rate (i.e. normal surface flux)
!          from the body force section
!------------------------------------------------------------------------------
           SourceFunc = 0.0d0
           Flux  = 0.0d0
           SubstantialSurface = .TRUE.
           BodyForce => GetBodyForce()
           IF ( UseBodyForce.AND.ASSOCIATED( BodyForce ) ) THEN
              SubstantialSurface = .FALSE.
              ! Accumulation/ablation is given in normal direction of surface:
              !---------------------------------------------------------------
              IF (NormalFlux) THEN 
                 SourceFunc(1:n) = GetReal( BodyForce, &
                      'Accumulation Ablation' )
                 ! Accumulation/ablation has to be computed from given flux:
                 !----------------------------------------------------------
              ELSE 
                 Flux(1,1:n) = GetReal( BodyForce, 'Accumulation Flux 1')
                 IF (DIM == 2) THEN
                    Flux(2,1:n) = GetReal( BodyForce, 'Accumulation Flux 2' )
                 ELSE
                    Flux(2,1:n) = 0.0d0
                 END IF
                 IF (DIM == 3) THEN
                    Flux(3,1:n) = GetReal( BodyForce, 'Accumulation Flux 3' )
                 ELSE
                    Flux(3,1:n) = 0.0d0
                 END IF
                 SourceFunc = 0.0d0
              END IF
           END IF
           IF( TransientSimulation) THEN
              ElemFreeSurf(1:n) = PreFreeSurf(FreeSurfPerm(NodeIndexes),1)
           END IF

!---------------------------------------------------------------
!          Previous non linear iter solution for that element
!---------------------------------------------------------------
           OldDy(1:n) = FreeSurf( FreeSurfPerm( NodeIndexes(1:n)))

!---------------------------------------------------------------
!          LocalLimit from Limit
!---------------------------------------------------------------
           LocalLimit(1:2,1:n) = Limit(1:2,BoundaryPerm(NodeIndexes))           
!          
!------------------------------------------------------------------------------
!          Get element local matrix, and rhs vector
!------------------------------------------------------------------------------
           CALL LocalMatrix( STIFF, MASS, FORCE,&
                SourceFunc, ElemFreeSurf, Velo, Element,&
                n, ElementNodes, NodeIndexes, TransientSimulation,&
                Flux, NormalFlux, SubstantialSurface, OldDy,&
                LimitedSolution, LocalLimit, UzawaParameter)

!------------------------------------------------------------------------------
!        If time dependent simulation add mass matrix to stiff matrix
!------------------------------------------------------------------------------
           IF ( Bubbles ) TimeForce  = FORCE
           IF ( TransientSimulation ) THEN
!------------------------------------------------------------------------------
!          NOTE: This will replace STIFF and LocalForce with the
!                combined information...
!------------------------------------------------------------------------------
              IF ( Bubbles ) FORCE = 0.0d0
              CALL Default1stOrderTime( MASS, STIFF, FORCE )
           END IF
!------------------------------------------------------------------------------
!        Update global matrices from local matrices
!------------------------------------------------------------------------------
           IF (Bubbles) THEN
              CALL Condensate( N, STIFF, FORCE, TimeForce )
              IF ( TransientSimulation ) CALL DefaultUpdateForce( TimeForce )
           END IF
!------------------------------------------------------------------------------
!      Update global matrix and rhs vector from local matrix & vector
!------------------------------------------------------------------------------
           CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
        END DO


!------------------------------------------------------------------------------
!     Neumann & Newton boundary conditions
!------------------------------------------------------------------------------
!
! MIND: In weak formulation it is not possible to prescribe a contact angle on
!       a boundary in this solver. This has to be taken care of in the boundary
!       condition for the stress tensor in the Navier-Stokes Solver. Thus, in
!       generally it does not make sense to prescribe a van Neumann type of
!       condition here.


!------------------------------------------------------------------------------
!    FinishAssemebly must be called after all other assembly steps, but before
!    Dirichlet boundary settings. Actually no need to call it except for
!    transient simulations.
!------------------------------------------------------------------------------
        CALL DefaultFinishAssembly()
        CALL DefaultDirichletBCs()


!------------------------------------------------------------------------------
!    Solve System
!------------------------------------------------------------------------------
        at = CPUTime() - at
        st = CPUTime()


        PrevNorm = Norm
        Solver % Variable % Norm = Norm

        Norm = DefaultSolve()

        Solver % Variable % Norm = Norm
        


        IF(NeedOldValues) THEN
        IF(LimitDisp) THEN 
           maxdh = -HUGE(maxdh)         
           DO i=1, Model % NumberOfNodes
              j = FreeSurfPerm(i)
              IF(j > 0) THEN
                 maxdh = MAX(maxdh, ABS(FreeSurf(j)-OldFreeSurf(j)))
              END IF
           END DO
           IF(maxdh > MaxDisp) THEN
              Relax = Relax * MaxDisp/maxdh
           END IF
           WRITE(Message,'(a,E8.2)') 'Maximum displacement',maxdh
           CALL Info( 'FreeSurfaceSolver', Message, Level=4 )
        END IF
        WRITE(Message,'(a,F8.2)') 'pp Relaxation factor',Relax
        CALL Info( 'FreeSurfaceSolver', Message, Level=4 )
           DO i=1, Model % NumberOfNodes
              j = FreeSurfPerm(i)
              IF(j > 0) THEN
                 FreeSurf(j) = Relax * FreeSurf(j) + (1-Relax) * OldFreeSurf(j)
              END IF
           END DO
        END IF


        st = CPUTIme()-st
        totat = totat + at
        totst = totst + st


        WRITE(Message,'(a,F8.2,F8.2)') 'Assembly: (s)', at, totat
        CALL Info( 'FreeSurfaceSolver', Message, Level=4 )
        WRITE(Message,'(a,F8.2,F8.2)') ' Solve:    (s)', st, totst
        CALL Info( 'FreeSurfaceSolver', Message, Level=4 )
        SubroutineVisited = SubroutineVisited + 1
        
        ! is the solution converged?
        !---------------------------
        IF ( PrevNorm + Norm /= 0.0d0 ) THEN
           RelativeChange = 2.0d0 * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
        ELSE
           RelativeChange = 0.0d0
        END IF
        WRITE( Message, * ) 'Result Norm   : ',Norm
        CALL Info( 'Free Surface Solve', Message, Level=4 )
        WRITE( Message, * ) 'Relative Change : ',RelativeChange
        CALL Info( 'Free Surface Solve', Message, Level=4 )

!!!!!!!! !!  pour limiter le deplacement
     If (NonlinearIter.eq.1) then
     IF (DispSteady) Then
      IF(maxdh.GT.SteadyDisp) then
              SteadyTol= RelativeChange/2._dp
      else
              SteadyTol= RelativeChange*2._dp
      Endif
      CALL  ListAddConstReal( SolverParams,  &
                'Steady State Convergence Tolerance', SteadyTol )
     ENDIF
     Endif
     write(*,*) 'Steady State Convergence Tolerance', SteadyTol
!!!!!!!!!!!!!!!!


        IF ( RelativeChange < NonlinearTol ) EXIT
     END DO ! of the nonlinear iteration
!------------------------------------------------------------------------------
   CONTAINS


!------------------------------------------------------------------------------
!==============================================================================
     SUBROUTINE LocalMatrix( STIFF, MASS, FORCE,&
          SourceFunc, OldFreeSurf, Velo, &
          Element, nCoord, Nodes, NodeIndexes, TransientSimulation,&
          Flux, NormalFlux, SubstantialSurface, OldDy,&
          LimitedSolution, LocalLimit, UzawaParameter)
!------------------------------------------------------------------------------
!    INPUT:  SourceFunc(:)   nodal values of the accumulation/ablation function
!            
!            Element         current element
!            n               number of nodes
!            Nodes           current node points
!
!    OUTPUT: STIFF(:,:)
!            MASS(:,:)
!            FORCE(:)
!------------------------------------------------------------------------------
!      external variables:
!      ------------------------------------------------------------------------
       REAL(KIND=dp) ::&
             STIFF(:,:), MASS(:,:), FORCE(:), SourceFunc(:), &
             Velo(:,:), OldFreeSurf(:), Flux(:,:), OldDy(:), &
             LocalLimit(:,:), UzawaParameter
             


       INTEGER :: nCoord, NodeIndexes(:)
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
       LOGICAL :: TransientSimulation,NormalFlux,SubstantialSurface, &
                  LimitedSolution(:)
!------------------------------------------------------------------------------
!      internal variables:
!      ------------------------------------------------------------------------
       REAL(KIND=dp) ::&
          Basis(2*nCoord),dBasisdx(2*nCoord,3),ddBasisddx(2*nCoord,3,3),&
          Vgauss(3), VMeshGauss(3), Source, gradFreeSurf(3), normGradFreeSurf,&
          FluxGauss(3),X,Y,Z,U,V,W,S,SqrtElementMetric,&
          SU(2*nCoord),SW(2*nCoord),Tau,hK,UNorm, Oldh, LowerLimit, &
          UpperLimit
       LOGICAL :: Stat, Bed(2)
       INTEGER :: i,j,t,p,q, n
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------


       FORCE = 0.0d0
       STIFF = 0.0d0
       MASS  = 0.0d0


       IF (Bubbles) THEN
          n = nCoord * 2
       ELSE
          n = nCoord
       END IF


       hK = ElementDiameter( Element, Nodes )
!
!      Numerical integration:
!      ----------------------
       IF (Bubbles) THEN
          IntegStuff = GaussPoints( Element, Element % type % gausspoints2)
       ELSE
          IntegStuff = GaussPoints( Element )
       END IF


       SU = 0.0d0
       SW = 0.0d0


       DO t = 1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
         S = IntegStuff % s(t)
!
!        Basis function values & derivatives at the integration point:
!        -------------------------------------------------------------
         stat = ElementInfo( Element,Nodes,U,V,W,SqrtElementMetric, &
              Basis,dBasisdx,ddBasisddx,.FALSE., Bubbles )


!        Correction from metric
!        ----------------------
          S = S * SqrtElementMetric


         IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
            X = SUM( Nodes % x(1:nCoord) * Basis(1:nCoord) )
            Y = SUM( Nodes % y(1:nCoord) * Basis(1:nCoord) )
            Z = SUM( Nodes % z(1:nCoord) * Basis(1:nCoord) )
            S = S * X
         END IF
!
!        Velocities and (norm of) gradient of free surface and source function 
!        at Gauss point
!        ---------------------------------------------------------------------

         gradFreeSurf=0.0d0
         Vgauss=0.0d0
         VMeshGauss=0.0d0


         DO i=1,DIM-1
!GAG!      gradFreeSurf(i) = SUM(dBasisdx(1:nCoord,i)*OldFreeSurf(1:nCoord))
           gradFreeSurf(i) = SUM(dBasisdx(1:nCoord,i)*OldDy(1:nCoord))
         END DO

         gradFreeSurf(DIM) = 1.0d0

         DO i=1,DIM
           Vgauss(i) = SUM( Basis(1:nCoord)*Velo(i,1:nCoord) )
         END DO

         IF (DIM==3) THEN
           normGradFreeSurf = SQRT(1.0d0 + gradFreeSurf(1)**2 + &
                gradFreeSurf(2)**2)
         ELSE
           normGradFreeSurf = SQRT(1.0d0 + gradFreeSurf(1)**2)
         END IF


         UNorm = SQRT( SUM( Vgauss(1:dim-1)**2 ) )
         Tau = hK / ( 2*Unorm )


         ! Srabilized 
         IF ( .NOT. Bubbles ) THEN
            DO p=1,n
               SU(p) = 0.0d0
               DO i=1,dim-1
                  SU(p) = SU(p) + Vgauss(i) * dBasisdx(p,i)
               END DO


               SW(p) = 0.0d0
               DO i=1,dim-1
                  SW(p) = SW(p) + Vgauss(i) * dBasisdx(p,i)
               END DO
            END DO
         END IF


!        Stiffness matrix:
!        -----------------
         DO p=1,n
           DO q=1,n
             DO i=1,DIM-1
               STIFF(p,q) = STIFF(p,q) + &
                       s * Vgauss(i) * dBasisdx(q,i) * Basis(p)
             END DO
             STIFF(p,q) =  STIFF(p,q) + s * Tau * SU(q) * SW(p)
           END DO
           
         END DO



!        Mass Matrix:
!        ------------
         IF ( TransientSimulation ) THEN
            DO p=1,n
              DO q=1,n
                MASS(p,q) = MASS(p,q) +  &
                         S * Basis(q) * (Basis(p) + Tau*SW(p))
              END DO
            END DO
         END IF


!        Get accumulation/ablation function if flux input is given
!        (i.e., calculate vector product between flux and normal)
!        --------------------------------------------------------- 
         IF (.NOT.(SubstantialSurface)) THEN
            IF (NormalFlux) THEN 
               Source = normGradFreeSurf * SUM( SourceFunc(1:nCoord) &
                    * Basis(1:nCoord) )
            ELSE
               DO i=1,dim
                  FluxGauss(i) = SUM(Basis(1:nCoord)*Flux(i,1:nCoord))
               END DO
               Source = SUM(FluxGauss(1:DIM)*gradFreeSurf(1:DIM))
            END IF
         ELSE
            Source = 0.0d0
         END IF


!        Assemble force vector:
!        ---------------------
         FORCE(1:n) = FORCE(1:n) &
              + (Vgauss(dim)+Source) * (Basis(1:n) + Tau*SW(1:n)) * s
      END DO
!
! add the bed conditions using a penality if necessary
! 
      Do  p = 1, n  
         If (LimitedSolution(1)) THEN      
               IF  ((OldDy(p) - LocalLimit(1,p)).LE.0.0d0) Then
                       STIFF(p,p) = STIFF(p,p) +  &
                                      UzawaParameter
                    If (Abs(LocalLimit(1,p)).GT.1.0/UzawaParameter) Then               
                       FORCE(p) = FORCE(p) +  &
                        LocalLimit(1,p) * UzawaParameter 
                    End If    
               END IF
         End If

         If ((LimitedSolution(2)) .And.  &
                        ((OldDy(p) - LocalLimit(2,p)).GT. 0.0e0)) Then
                       STIFF(p,p) = STIFF(p,p) +  &
                                      UzawaParameter 
                       FORCE(p) = FORCE(p) +  &
                       UzawaParameter * LocalLimit(2,p) 
         End If
      End Do
!------------------------------------------------------------------------------
    END SUBROUTINE LocalMatrix


!==============================================================================
!------------------------------------------------------------------------------
  END SUBROUTINE FSurfSolverPenalty
!------------------------------------------------------------------------------ 
