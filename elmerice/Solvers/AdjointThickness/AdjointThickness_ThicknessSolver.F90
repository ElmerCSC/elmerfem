!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! ******************************************************************************
! *
! *  Authors: F. Gillet-Chaulet
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: Dec. 2020
! * 
! *****************************************************************************
!-----------------------------------------------------------------------------
SUBROUTINE AdjointThickness_ThicknessSolver( Model,Solver,dt,TransientSimulation )
  USE DefUtils
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
  TYPE(ValueList_t), POINTER :: BodyForce, SolverParams

  INTEGER :: NMAX,NSDOFs
  INTEGER :: istat
  INTEGER :: t,i,j,n
  LOGICAL :: ConvectionVar, Found

  REAL(KIND=dp) :: Norm

  CHARACTER(LEN=MAX_NAME_LEN)  :: FlowSolName
  TYPE(Variable_t), POINTER :: FlowSol
  REAL(KIND=dp), POINTER :: FlowSolution(:)
  INTEGER, POINTER :: FlowPerm(:)

  REAL(KIND=dp), ALLOCATABLE,SAVE :: STIFF(:,:),FORCE(:),LOAD(:),Velo(:,:)

  CHARACTER(LEN=MAX_NAME_LEN)  :: SolverName= 'AdjointThickness_ThicknessSolver'

  TYPE(Nodes_t),SAVE   :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement
  INTEGER, POINTER :: NodeIndexes(:)

  LOGICAL, SAVE :: AllocationsDone = .FALSE.

  !------------------------------------------------------------------------------
  !    Go
  !------------------------------------------------------------------------------
  SolverParams => GetSolverParams()

  !------------------------------------------------------------------------------
  !    check incompatibilities
  !------------------------------------------------------------------------------
  IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
          CALL FATAL(SolverName,'sorry only for cartesian coordinate system')
  END IF
  IF (TransientSimulation) THEN
          CALL FATAL(SolverName,'sorry works only in steady state')
  ENDIF

  !------------------------------------------------------------------------------
  !    Allocate some permanent storage, this is done first time only
  !------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
     NMAX = Model % MaxElementNodes

     IF ( AllocationsDone ) THEN
        DEALLOCATE( ElementNodes % x,    &
             ElementNodes % y,    &
             ElementNodes % z,    &
             FORCE,    &
             STIFF, &
             LOAD,&
             Velo)
     END IF

     ALLOCATE( ElementNodes % x( NMAX ),    &
          ElementNodes % y( NMAX ),    &
          ElementNodes % z( NMAX ),    &
          FORCE( NMAX ),    &
          STIFF( NMAX, NMAX ), &
          LOAD(NMAX) , &
          Velo( 2, NMAX ), &
          STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal(SolverName,'Memory allocation error 1, Aborting.')
     END IF

     CALL Info(SolverName,'Memory allocations done' )
     AllocationsDone = .TRUE.
  END IF

  !------------------------------------------------------------------------------
  !    Get Flow solution
  !------------------------------------------------------------------------------
   ConvectionVar=.True.
   FlowSolName =  GetString(SolverParams ,'Flow Solution Name', Found)
   IF(.NOT.Found) THEN        
        WRITE(Message,'(A)') &
           '<Flow Solution Name> Not Found; will look for <convection velocity> in body forces'
        CALL Info(SolverName,Message,level=10)
        ConvectionVar=.False.
        NSDOFS=GetInteger(SolverParams ,'Convection Dimension',Found)
        IF(.NOT.Found) &
            CALL Fatal(SolverName,'if <Flow Solution Name> not given prescribe <Convection Dimension>')
   ELSE
        FlowSol => VariableGet( Solver % Mesh % Variables, FlowSolName,UnFoundFatal=.TRUE.)
        FlowPerm     => FlowSol % Perm
        NSDOFs     =  FlowSol % DOFs
        FlowSolution => FlowSol % Values
   END IF


   CALL DefaultInitialize()
   !------------------------------------------------------------------------------
   !    Do the assembly
   !------------------------------------------------------------------------------
   DO t=1,Solver % NumberOfActiveElements

        CurrentElement => GetActiveElement(t)
        n = GetElementNOFNodes(CurrentElement)
        NodeIndexes => CurrentElement % NodeIndexes

        ! set coords of highest occurring dimension to zero (to get correct path element)
        !-------------------------------------------------------------------------------
        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        IF (NSDOFs == 1) THEN
           ElementNodes % y(1:n) = 0.0
           ElementNodes % z(1:n) = 0.0
        ELSE IF (NSDOFs == 2) THEN
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = 0.0_dp
        ELSE
           WRITE(Message,'(a,i0,a)')&
                'It is not possible to compute Thickness evolution if Flow Sol DOFs=',&
                NSDOFs, ' . Aborting'
           CALL Fatal( SolverName, Message) 
           STOP   
        END IF

        ! get pointers on BF
        !---------------------------------------------------------------
        BodyForce => GetBodyForce(CurrentElement)

        ! get flow soulution and velocity field from it
        !----------------------------------------------
        Velo = 0.0_dp
        !----------------------------------------------------
        ! get velocity profile
        IF (ConvectionVar) Then
          DO i=1,n
             j = NSDOFs*FlowPerm(NodeIndexes(i))
              !2D problem - 1D Thickness evolution
              IF(NSDOFs == 1) THEN 
                 Velo(1,i) = FlowSolution( j ) 
                 Velo(2,i) = 0.0_dp
              !2D problem 
              ELSE IF (NSDOFs == 2) THEN
                 Velo(1,i) = FlowSolution( j-1 ) 
                 Velo(2,i) = FlowSolution( j ) 
              ELSE
                 WRITE(Message,'(a)') 'Velocity should be size 1 or 2'
                 CALL Fatal( SolverName, Message)
              END IF
           END DO
       ELSE
        IF (ASSOCIATED( BodyForce ) ) THEN
          Velo(1,1:n) = ListGetReal( BodyForce, 'Convection Velocity 1',n, NodeIndexes,UnFoundFatal=.TRUE. )
          if (NSDOFs.eq.2) &
           Velo(2,1:n) = ListGetReal( BodyForce, 'Convection Velocity 2',n, NodeIndexes,UnFoundFatal=.TRUE. )
        END IF
       END IF
        !------------------------------------------------------------------------------
        !      get the accumulation/ablation rate (i.e. normal surface flux)
        !      from the body force section
        !------------------------------------------------------------------------------
        LOAD=0.0_dp
        IF (ASSOCIATED( BodyForce ) ) THEN
              LOAD(1:n) = LOAD(1:n) +   &
                      GetReal( BodyForce, 'Top Surface Accumulation', Found )
              LOAD(1:n) = LOAD(1:n) +   &
                      GetReal( BodyForce, 'Bottom Surface Accumulation', Found )
        END IF

        !------------------------------------------------------------------------------
        !      Get element local matrix, and rhs vector
        !------------------------------------------------------------------------------
        CALL LocalMatrix( STIFF, FORCE,LOAD,&
                Velo, NSDOFs, &
             CurrentElement, n, ElementNodes, NodeIndexes)

        !------------------------------------------------------------------------------
        CALL DefaultUpdateEquations( STIFF, FORCE )
        !------------------------------------------------------------------------------
        
     END DO ! End loop bulk elements
     CALL DefaultFinishBulkAssembly()

     !------------------------------------------------------------------------------
     !     Neumann & Newton boundary conditions
     !------------------------------------------------------------------------------
     !
     ! MIND: In weak formulation it is not possible to prescribe a contact angle on
     !       a boundary in this solver. This has to be taken care of in the boundary
     !       condition for the stress tensor in the Navier-Stokes Solver. Thus, in
     !       generally it does not make sense to prescribe a Neumann type of
     !       condition here.

     !------------------------------------------------------------------------------
     !    FinishAssemebly must be called after all other assembly steps, but before
     !    Dirichlet boundary settings. Actually no need to call it except for
     !    transient simulations.
     !------------------------------------------------------------------------------
     CALL DefaultFinishAssembly()
     CALL DefaultDirichletBCs()

     Norm = DefaultSolve()

     !------------------------------------------------------------------------------
   CONTAINS

     !------------------------------------------------------------------------------
     !==============================================================================
     SUBROUTINE LocalMatrix( STIFF, FORCE,LOAD,&
                     Velo, NSDOFs, &
                     Element, n, Nodes, NodeIndexes)
       !------------------------------------------------------------------------------
       !    INPUT:  LOAD(:)   nodal values of the accumulation/ablation function
       !            
       !            Element         current element
       !            n               number of nodes
       !            Nodes           current node points
       !
       !    OUTPUT: STIFF(:,:)
       !            FORCE(:)
       !------------------------------------------------------------------------------
       !      external variables:
       !      ------------------------------------------------------------------------
       REAL(KIND=dp) :: STIFF(:,:),FORCE(:), LOAD(:), Velo(:,:)

       INTEGER :: n, NodeIndexes(:), NSDOFs
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element

       !------------------------------------------------------------------------------
       !      internal variables:
       !      ------------------------------------------------------------------------
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
       REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
       REAL(KIND=dp) :: SU(n),SW(n)
       REAL(KIND=dp) :: Vgauss(2), UNorm,divu,Source
       REAL(KIND=dp) :: U,V,W,S,SqrtElementMetric
       REAL(KIND=dp) :: Tau,hk
       INTEGER :: i,j,t,p,q
       LOGICAL :: stat
       !------------------------------------------------------------------------------

       FORCE = 0.0_dp
       STIFF = 0.0_dp

       hK = ElementDiameter( Element, Nodes )

       !
       !      Numerical integration:
       !      ----------------------
       IntegStuff = GaussPoints( Element )

       SU = 0.0_dp
       SW = 0.0_dp

       DO t = 1,IntegStuff % n
          U = IntegStuff % u(t)
          V = IntegStuff % v(t)
          W = IntegStuff % w(t)
          S = IntegStuff % s(t)
          !
          !        Basis function values & derivatives at the integration point:
          !        -------------------------------------------------------------
          stat = ElementInfo( Element,Nodes,U,V,W,SqrtElementMetric, &
               Basis,dBasisdx, Bubbles=.False. )

          !        Correction from metric
          !        ----------------------
          S = S * SqrtElementMetric

          !
          !        Velocities and (norm of) gradient of free surface and source function 
          !        at Gauss point
          !        ---------------------------------------------------------------------

          Vgauss=0.0_dp
          DO i=1,NSDOFs
             Vgauss(i) = SUM( Basis(1:n)*(Velo(i,1:n)))
          END DO

          divu = 0.0_dp
          DO i=1,NSDOFs
             divu = divu +  SUM( dBasisdx(1:n,i)*(Velo(i,1:n)))
          END DO

          UNorm = SQRT( SUM( Vgauss(1:NSDOFs)**2 ) )
          Tau=0._dp
          IF (UNorm .GT. 0.0_dp) Tau = hK / ( 2*Unorm )

          DO p=1,n
             SU(p) = 0.0_dp
             SW(p) = 0.0_dp
             DO i=1,NSDOFs
                   SU(p) = SU(p) + Vgauss(i) * dBasisdx(p,i)
                   SW(p) = SW(p) + Vgauss(i) * dBasisdx(p,i)
             END DO
          END DO

          !        Stiffness matrix:
          !        -----------------
          DO p=1,n
             DO q=1,n
                DO i=1,NSDOFs
                   STIFF(p,q) = STIFF(p,q) + &
                        s * Vgauss(i) * dBasisdx(q,i) * Basis(p)
                END DO
                STIFF(p,q) =  STIFF(p,q) + s * Tau * SU(q) * SW(p)
                STIFF(p,q) =  STIFF(p,q) + s * divu * Basis(q) * (Basis(p) + Tau*SW(p))
             END DO
          END DO

          !        Get accumulation/ablation function 
          !        --------------------------------------------------------- 
          Source = 0.0_dp
          Source=SUM(Basis(1:n)*LOAD(1:n))

          !        Assemble force vector:
          !        ---------------------
          FORCE(1:n) = FORCE(1:n) &
               + Source * (Basis(1:n) + Tau*SW(1:n)) * s

       END DO

       !------------------------------------------------------------------------------
     END SUBROUTINE LocalMatrix

     !------------------------------------------------------------------------------
   END SUBROUTINE AdjointThickness_ThicknessSolver
!------------------------------------------------------------------------------
