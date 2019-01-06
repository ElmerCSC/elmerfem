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
! *****************************************************************************
!
!******************************************************************************
! *
! *  Authors: Mikko, Lyly, Juha Ruokolainen, Thomas Zwinger
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 04 Apr 2004
! *
! ****************************************************************************

!------------------------------------------------------------------------------
!>  Advection-reaction equation solver for scalar fields with discontinuous Galerkin method.
!> \ingroup Solvers
!------------------------------------------------------------------------------
  SUBROUTINE AdvectionReactionSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
     USE DefUtils

     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model            !< All model information (mesh, materials, BCs, etc...)
     TYPE(Solver_t), TARGET :: Solver  !< Linear & nonlinear equation solver options
     REAL(KIND=dp) :: dt               !< Timestep size for time dependent simulations
     LOGICAL :: Transient    !< Steady state or transient simulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: BC, BodyForce, Material,&
          Equation, Constants, SolverParams

     TYPE(Element_t), POINTER :: Element, Face, &
       ParentElement, P1, P2

     LOGICAL :: AllocationsDone = .FALSE., Found, Stat, &
          LimitSolution, FoundLowerLimit, FoundUpperLimit
     INTEGER :: Active, DIM,NonLinearIterMin,NonlinearIterMax,iter,&
          CorrectedLowerLimit,CorrectedUpperLimit
     INTEGER :: n1,nd1,nd2,n2, k, n, nd,t, istat, i, j, dummyInt, NumberOfFAces, Indexes(128)
#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: Norm,RelativeChange,at,at0,totat,st,totst,&
          OriginalValue 
#else
     REAL(KIND=dp) :: Norm,RelativeChange,at,at0,totat,st,totst,CPUTime,RealTime,&
          OriginalValue 
#endif
     REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), LOAD(:), &
              FORCE(:), Velo(:,:), MeshVelo(:,:), Gamma(:), Ref(:), &
              UpperLimit(:), LowerLimit(:)
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Variable_t), POINTER ::  Var
     CHARACTER(LEN=MAX_NAME_LEN) ::  VariableName, SolverName, ExpVariableName

     SAVE MASS, STIFF, LOAD, FORCE, Velo, MeshVelo, Gamma, &
          AllocationsDone, DIM, VariableName, SolverName, &
          UpperLimit, LowerLimit
!*******************************************************************************

     TYPE( Element_t ), POINTER :: Faces(:)
     TYPE(Nodes_t) :: ElementNodes

     INTEGER :: DOFs
!*******************************************************************************
     SolverName = 'AdvectionReaction ('// TRIM(Solver % Variable % Name) // ')'
     VariableName = TRIM(Solver % Variable % Name)
     WRITE(Message,'(A,A)')&
          'AdvectionReactionSolver for variable ', VariableName
     CALL INFO(SolverName,Message,Level=1)

     Mesh => GetMesh()

     IF ( CoordinateSystemDimension() == 2 ) THEN
        Faces => Mesh % Edges
        NumberOfFaces = Mesh % NumberOfEdges
     ELSE
        Faces => Mesh % Faces
        NumberOfFaces = Mesh % NumberOfFaces
     END IF

     ! Initialize & allocate some permanent storage, this is done first time only:
     !----------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN

        N = 2 * MAX(Mesh % MaxElementDOFs, Mesh % MaxElementNodes )
        ALLOCATE( FORCE(N), MASS(n,n), STIFF(N,N), LOAD(N),  &
                  Velo(3,N), MeshVelo( 3,N ), Gamma(n), &
                  UpperLimit(n), LowerLimit(n), STAT = istat )
        
       IF ( istat /= 0 ) THEN
          CALL FATAL(SolverName,'Memory allocation error.' )
       ELSE
          CALL INFO(SolverName,'Memory allocation done',Level=1 )
       END IF
       DIM = CoordinateSystemDimension()
       AllocationsDone = .TRUE.
     END IF

     !------------------------------------------------------------------------------
     !    Read physical and numerical constants and initialize 
     !------------------------------------------------------------------------------
     Constants => GetConstants()
     SolverParams => GetSolverParams()

     NonlinearIterMax = GetInteger(   SolverParams, &
                     'Nonlinear System Max Iterations', Found )
     IF ( .NOT.Found ) THEN
        CALL WARN(SolverName,'No > Nonlinear System Max Iterations < found. Setting 1')
        NonlinearIterMax = 1
     END IF

     NonlinearIterMin = GetInteger(   SolverParams, &
                     'Nonlinear System Min Iterations', Found )
     IF ( .NOT.Found ) THEN        
        CALL WARN(SolverName,'No >Nonlinear System Min Iterations< found. Setting 1')
        NonlinearIterMin = 1
     ELSE IF (NonlinearIterMin > NonlinearIterMax) THEN
        CALL WARN(SolverName,&
             '>Nonlinear System Min Iterations< is exceeding >Nonlinear System Max Iterations<.')
         CALL WARN(SolverName,&
             'First is being reset to the latter.')
        NonlinearIterMin = NonlinearIterMax
     END IF

     ExpVariableName = GetString(SolverParams , 'Exported Variable 1', Found )
     IF (.NOT.Found) &
        CALL FATAL(SolverName,'No value > Exported Variable 1 < found in Solver')

     !------------------------------------------------------------------------------
     !       non-linear system iteration loop
     !------------------------------------------------------------------------------
     totat = 0; totst = 0
     DO iter=1,NonlinearIterMax
        ! Assembly of the bulk elements:
        !-------------------------------
        !------------------------------------------------------------------------------
        ! print out some information
        !------------------------------------------------------------------------------
        at  = CPUTime()
        at0 = RealTime()

        CALL Info( SolverName, ' ', Level=4 )
        CALL Info( SolverName, ' ', Level=4 )
        CALL Info( SolverName, '-------------------------------------',Level=4 )
        WRITE( Message,'(A,I3,A,I3)') &
             'Nonlinear iteration no.', iter,' of max',NonlinearIterMax
        CALL Info( SolverName, Message, Level=4 )
        CALL Info( SolverName, '-------------------------------------',Level=4 )
        CALL Info( SolverName, ' ', Level=4 )
        CALL Info( SolverName, 'Starting Assembly...', Level=4 )


        CALL DefaultInitialize()
        Active = GetNOFActive()
        DO t = 1, Active  
           !------------------------------------------------------------------------------
           ! write some info on status of assembly
           !------------------------------------------------------------------------------
           IF ( RealTime() - at0 > 1.0 ) THEN
!             WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
!                  (Active-t) / &
!                  (1.0*Active)), ' % done'
!
!             CALL Info( SolverName, Message, Level=5 )

              at0 = RealTime()
           END IF
           !------------------------------------------------------------------------------
           ! assign pointers and get number of nodes in element
           !------------------------------------------------------------------------------  
           Element => GetActiveElement(t)
           n  = GetElementNOfNodes(Element)
           nd = GetElementNOfDOFs(Element)

           Material => GetMaterial(Element)
           Equation => GetEquation(Element)
           BodyForce => GetBodyForce(Element)

           !------------------------------------------------------------------------------         
           LOAD(1:n) = GetReal( BodyForce, TRIM(VariableName) // ' Source', Found )    

           !------------------------------------------------------------------------------
           ! Get convection and mesh velocity
           CALL GetLocalALEVelocity(Velo,MeshVelo,SolverName,Material,&
                Equation,Solver,Model,Element)

           !--------------------------
           ! get reaction coefficient
           !--------------------------
           Gamma(1:n)  = GetReal( Material, TRIM(VariableName)//' Gamma', Found )

           CALL LocalMatrix(MASS, STIFF, FORCE, LOAD, Velo, Gamma, Element, n,nd) 
           IF ( Transient ) CALL Default1stOrderTime(MASS, STIFF, FORCE)
           CALL DefaultUpdateEquations(STIFF, FORCE)
        END DO


        ! Assembly of the face terms:
        !----------------------------
        FORCE = 0.0_dp
        DO t=1,NumberOfFaces
           Face => Faces(t)
           IF ( .NOT. ActiveBoundaryElement(Face) ) CYCLE

           P1  => Face % BoundaryInfo % Left
           P2 => Face % BoundaryInfo % Right
           IF ( ASSOCIATED(P2) .AND. ASSOCIATED(P1) ) THEN
              n  = GetElementNOFNodes(Face)
              n1 = GetElementNOFNodes(P1)
              n2 = GetElementNOFNodes(P2)

              nd = GetElementNOFDOFs(Face)
              nd1 = GetElementNOFDOFs(P1)
              nd2 = GetElementNOFDOFs(P2)

              !------------------------------------------------------------------------------
              ! Get convection and mesh velocity
              !------------------------------------------------------------------------------ 
              Material => GetMaterial(P1)
              Equation => GetEquation(P1)          
              CALL GetLocalALEVelocity(Velo,MeshVelo,SolverName,Material,&
                   Equation,Solver,Model,Face)

              CALL LocalJumps( STIFF,Face,n,nd,P1,n1,nd1,P2,n2,nd2,Velo )
              CALL DefaultUpdateEquations( STIFF, FORCE, Face )
           END IF
        END DO

        ! Loop over the boundary elements:
        !---------------------------------
        DO t=1,Mesh % NumberOfBoundaryElements
           Element => GetBoundaryElement(t)
           IF( .NOT. ActiveBoundaryElement() )  CYCLE
           IF( GetElementFamily(Element) < dim ) CYCLE 

           ParentElement => Element % BoundaryInfo % Left
           IF ( .NOT. ASSOCIATED( ParentElement ) ) &
                ParentElement => Element % BoundaryInfo % Right

           n1  = GetElementNOFnodes(ParentElement)
           n   = GetElementNOFNodes(Element)

           nd  = GetElementNOFDOFs(Element)
           nd1 = GetElementNOFDOFs(ParentElement)
           !------------------------------------------------------------------------------
           ! Get convection and mesh velocity
           !------------------------------------------------------------------------------ 
           Material => GetMaterial(ParentElement)
           Equation => GetEquation(ParentElement)        
           CALL GetLocalALEVelocity(Velo,MeshVelo,SolverName,Material,&
                Equation,Solver,Model,Element)

           BC => GetBC()
           LOAD = 0.0_dp
           Found = .FALSE.
           IF ( ASSOCIATED(BC) ) THEN
             LOAD(1:n) = GetReal(BC, Solver % Variable % Name, Found)
           END IF

           MASS = 0.0_dp
           CALL LocalMatrixBoundary(  STIFF, FORCE, LOAD, &
                Element, n, nd, ParentElement, n1, nd1, Velo, Found )

           IF ( Transient ) CALL Default1stOrderTime(MASS, STIFF, FORCE)
           CALL DefaultUpdateEquations(STIFF, FORCE)
        END DO

        CALL DefaultFinishAssembly()
        CALL Info( SolverName, 'Assembly done', Level=4 )

        !------------------------------------------------------------------------------
        !     Solve the system and check for convergence
        !------------------------------------------------------------------------------
        at = CPUTime() - at
        st = CPUTime()

        ! Solve the system:
        !------------------
        Norm = DefaultSolve()

        st = CPUTIme()-st
        totat = totat + at
        totst = totst + st
        WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
        CALL Info( SolverName, Message, Level=4 )
        WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
        CALL Info( SolverName, Message, Level=4 )

        RelativeChange = Solver % Variable % NonlinChange 
        
        IF ( Solver % Variable % NonlinConverged == 1 )  THEN 
           WRITE(Message,'(A,I6,A,I6,A)') &
                'Nonlinear iteration converged after ', iter, &
                ' out of max ',NonlinearIterMax,' iterations'
           CALL INFO(SolverName,Message)
           EXIT
        ELSE IF ((iter .EQ. NonlinearIterMax) .AND. (NonlinearIterMax > 1)) THEN
           CALL WARN(SolverName,'Maximum nonlinear iterations reached, but system not converged')
        END IF

     !-----------------------------------------------
     END DO ! End of nonlinear iteration loop
     !----------------------------------------------


     ! Average the elemental results to nodal values:
     !-----------------------------------------------
     Var => VariableGet( Mesh % Variables,TRIM(ExpVariableName))
     IF ( ASSOCIATED( Var ) ) THEN
       n1 = Mesh % NumberOfNodes
       ALLOCATE( Ref(n1) )
       Ref = 0
       Var % Type = Variable_on_nodes

       IF ( ASSOCIATED( Var % Perm, Solver % Variable % Perm ) ) THEN
         ALLOCATE( Var % Perm(SIZE(Solver % Variable % Perm))  )
         Var % Perm = 0
         DO i = 1,n1
            Var % Perm(i) = i
         END DO
       END IF

       Var % Values = 0.0d0
       DO t=1,Active
          Element => GetActiveElement(t) 
          n = GetElementDOFs(Indexes)
          n = GetElementNOFNodes()
          DO i=1,n
             k = Element % NodeIndexes(i)
             Var % Values(k) = Var % Values(k) + &
               Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) )
             Ref(k) = Ref(k) + 1
          END DO
       END DO

       WHERE( Ref > 0 )
         Var % Values(1:n1) = Var % Values(1:n1) / Ref
       END WHERE
       DEALLOCATE( Ref )
    ELSE
       WRITE(Message,'(A,A,A)') 'Exported Variable >',TRIM(VariableName) //   'Nodal Result','< not found'
       CALL FATAL(SolverName,Message)
    END IF

   CONTAINS

!------------------------------------------------------------------------------      
     SUBROUTINE LocalMatrix(MASS, STIFF, FORCE, LOAD, Velo, Gamma, Element, n,nd)
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: MASS(:,:), STIFF(:,:), FORCE(:), &
                  LOAD(:), Velo(:,:), Gamma(:)
       INTEGER :: n,nd
       TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3)
       REAL(KIND=dp) :: detJ,U,V,W,S,A,L,cu(3),g
       LOGICAL :: Stat
       INTEGER :: i,p,q,t,dim
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
       TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
       dim = CoordinateSystemDimension()

       FORCE = 0.0_dp
       STIFF = 0.0_dp
       MASS  = 0.0_dp

       CALL GetElementNodes(Nodes, Element)
!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints(Element)
 
       DO t=1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
         S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
                 Basis, dBasisdx )

         S = S * detJ
         L = SUM( LOAD(1:n) *  Basis(1:n) )
         g = SUM( Basis(1:n) * Gamma(1:n) )

         cu = 0.0_dp
         DO i=1,dim
           cu(i) = SUM(Basis(1:n) * Velo(i,1:n))
         END DO
!------------------------------------------------------------------------------
!        The advection-reaction equation: dc/dt + grad(u . c) + gamma c = s
!------------------------------------------------------------------------------
         DO p=1,nd
            DO q=1,nd
              MASS(p,q)  = MASS(p,q)  + s * Basis(q) * Basis(p)
              STIFF(p,q) = STIFF(p,q) + s * g * Basis(q) * Basis(p)
              DO i=1,dim
                STIFF(p,q) = STIFF(p,q) - s * cu(i) * Basis(q) * dBasisdx(p,i)
              END DO
            END DO
         END DO
         FORCE(1:nd) = FORCE(1:nd) + s*L*Basis(1:nd)
!------------------------------------------------------------------------------
      END DO
!------------------------------------------------------------------------------
     END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE FindParentUVW(Face, nFace, Parent, nParent, U, V, W, Basis)
!------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(Element_t) :: Face, Parent
      INTEGER :: nFace, nParent
      REAL( KIND=dp ) :: U, V, W, Basis(:)
!------------------------------------------------------------------------------
      INTEGER :: i,j
      REAL(KIND=dp) :: ParentU(nFace), ParentV(nFace), ParentW(nFace)
!------------------------------------------------------------------------------
      DO i = 1,nFace
        DO j = 1,nParent
          IF (Face % NodeIndexes(i) == Parent % NodeIndexes(j)) THEN
            ParentU(i) = Parent % Type % NodeU(j)
            ParentV(i) = Parent % Type % NodeV(j)
            ParentW(i) = Parent % Type % NodeW(j)
            EXIT
          END IF
        END DO
      END DO
      U = SUM(Basis(1:nFace)*ParentU(1:nFace))
      V = SUM(Basis(1:nFace)*ParentV(1:nFace))
      W = SUM(Basis(1:nFace)*ParentW(1:nFace))
!------------------------------------------------------------------------------      
    END SUBROUTINE FindParentUVW
!------------------------------------------------------------------------------      


!------------------------------------------------------------------------------
    SUBROUTINE LocalJumps(STIFF,Face,n,nd,P1,n1,nd1,P2,n2,nd2,Velo)
!------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=dp) :: STIFF(:,:), Velo(:,:)
      INTEGER :: n,n1,n2,nd,nd1,nd2
      TYPE(Element_t), POINTER :: Face, P1, P2
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Basis(nd),P1Basis(nd1),P2Basis(nd2)
      REAL(KIND=dp) :: Jump(nd1+nd2), Average(nd1+nd2)
      REAL(KIND=dp) :: detJ, U, V, W, S, Udotn, xx, yy
      LOGICAL :: Stat
      INTEGER :: i, j, p, q, dim, t, nFace, nParent
      TYPE(GaussIntegrationPoints_t) :: IntegStuff
      REAL(KIND=dp) :: Normal(3), cu(3), LeftOut(3)

      TYPE(Nodes_t), SAVE :: FaceNodes, P1Nodes, P2Nodes
!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()
      STIFF = 0.0_dp

      CALL GetElementNodes(P1Nodes, P1)
      CALL GetElementNodes(P2Nodes, P2)
      CALL GetElementNodes(FaceNodes, Face)
!------------------------------------------------------------------------------
!     Numerical integration over the edge
!------------------------------------------------------------------------------
      IntegStuff = GaussPoints(Face)

      LeftOut(1) = SUM(P1Nodes % x(1:n1)) / n1
      LeftOut(2) = SUM(P1Nodes % y(1:n1)) / n1
      LeftOut(3) = SUM(P1Nodes % z(1:n1)) / n1
      LeftOut(1) = SUM(FaceNodes % x(1:n)) / n - LeftOut(1)
      LeftOut(2) = SUM(FaceNodes % y(1:n)) / n - LeftOut(2)
      LeftOut(3) = SUM(FaceNodes % z(1:n)) / n - LeftOut(3)

      DO t=1,IntegStuff % n
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)
        S = IntegStuff % s(t)

        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo(Face, FaceNodes, U, V, W, detJ,Basis)
        S = S * detJ

        Normal = NormalVector(Face, FaceNodes, U, V, .FALSE.)
        IF (SUM(LeftOut*Normal)<0) Normal = -Normal

        ! Find basis functions for the parent elements:
        ! ---------------------------------------------
        CALL FindParentUVW(Face,n,P1,n1,U,V,W,Basis)
        stat = ElementInfo(P1, P1Nodes, U, V, W, detJ, P1Basis)

        CALL FindParentUVW(Face,n,P2,n2,U,V,W,Basis)
        stat = ElementInfo(P2, P2Nodes, U, V, W, detJ, P2Basis)

        ! Integrate jump terms:
        ! ---------------------
        Jump(1:nd1) = P1Basis(1:nd1)
        Jump(nd1+1:nd1+nd2) = -P2Basis(1:nd2)

        Average(1:nd1) = P1Basis(1:nd1) / 2
        Average(nd1+1:nd1+nd2) = P2Basis(1:nd2) / 2

        cu = 0.0_dp
        DO i=1,dim
          cu(i) = SUM( Velo(i,1:n) * Basis(1:n) )
        END DO
        Udotn = SUM( Normal * cu )

        DO p=1,nd1+nd2
          DO q=1,nd1+nd2
            STIFF(p,q) = STIFF(p,q) + s * Udotn * Average(q) * Jump(p)
            STIFF(p,q) = STIFF(p,q) + s * ABS(Udotn)/2 * Jump(q) * Jump(p)
          END DO
        END DO
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE LocalJumps
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    SUBROUTINE LocalMatrixBoundary( STIFF, FORCE, LOAD, &
        Element, n, nd, ParentElement, np, npd,Velo, InFlowBC )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: STIFF(:,:),  FORCE(:), LOAD(:), Velo(:,:)
     INTEGER :: n, np,nd,npd
     LOGICAL :: InFlowBC
     TYPE(Element_t), POINTER :: Element, ParentElement
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(nd), ParentBasis(npd)
     INTEGER :: i,j,k,m,p,q,t,dim

     REAL(KIND=dp) :: Normal(3), g, L, Udotn, cu(3), detJ,U,V,W,S
     LOGICAL :: Stat, Inflow
     TYPE(GaussIntegrationPoints_t) :: IntegStuff

     TYPE(Nodes_t) :: Nodes, ParentNodes
!------------------------------------------------------------------------------
     dim = CoordinateSystemDimension()
     FORCE = 0.0_dp
     STIFF = 0.0_dp
 
     CALL GetElementNodes(Nodes, Element)
     CALL GetElementNodes(ParentNodes, ParentElement)

     Normal = NormalVector( Element, Nodes, 0.0_dp, 0.0_dp, .TRUE. ) 
     DO i=1,3
       cu(i) = SUM(Velo(i,1:n)) / n
     END DO
     Inflow = InFlowBC .AND. SUM( Normal * cu ) < 0.0_dp

     ! Numerical integration:
     !-----------------------
     IntegStuff = GaussPoints(Element)

     DO t=1,IntegStuff % n
       U = IntegStuff % u(t)
       V = IntegStuff % v(t)
       W = IntegStuff % w(t)
       S = IntegStuff % s(t)

       Normal = NormalVector(Element, Nodes, U, V, .TRUE.) 

       ! Basis function values & derivatives at the integration point:
       ! -------------------------------------------------------------
       stat = ElementInfo(Element, Nodes, U, V, W, detJ,Basis)
       S = S * detJ

       CALL FindParentUVW(Element, n, ParentElement, np, U, V, W, Basis)
       stat = ElementInfo( ParentElement, ParentNodes, U, V, W, &
            detJ, ParentBasis)

       L = SUM(LOAD(1:n) * Basis(1:n))
       cu = 0.0_dp
       DO i=1,dim
         cu(i)  = SUM(Velo(i,1:n) * Basis(1:n))
       END DO
       Udotn = SUM(Normal*cu)

       DO p = 1,npd
         IF ( Inflow ) THEN
            FORCE(p) = FORCE(p) - s*Udotn*L*ParentBasis(p)
         ELSE
           DO q=1,npd
             STIFF(p,q) = STIFF(p,q) + s*Udotn*ParentBasis(q)*ParentBasis(p)
           END DO
         END IF
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
   SUBROUTINE GetLocalALEVelocity(Velo,MeshVelo,SolverName,Material,&
        Equation,Solver,Model,Element)
     IMPLICIT NONE
     !------------------------------------------------------------------------------
     REAL (KIND=dp) :: Velo(:,:), MeshVelo(:,:)
     CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
     TYPE(ValueList_t), POINTER :: Material,Equation
     TYPE(Solver_t), TARGET :: Solver
     TYPE(Model_t), TARGET :: Model
     TYPE(Element_t),POINTER :: Element
     !------------------------------------------------------------------------------
     CHARACTER(LEN=MAX_NAME_LEN) :: ConvectionFlag, FlowSolName
     INTEGER :: i,j,k,N,FlowDOFs
     INTEGER, POINTER :: FlowPerm(:)
     REAL(KIND=dp), POINTER :: FlowSolution(:)
     TYPE(Variable_t), POINTER :: FlowSol
     LOGICAL :: Found
     !------------------------------------------------------------------------------
     n  = GetElementNOFNodes( Element )
     ConvectionFlag = GetString( Equation, 'Convection', Found )
     IF (.NOT. Found) &
          CALL FATAL(SolverName, 'No string for keyword >Convection< found in Equation')
     Velo = 0.0d00
     ! constant (i.e., in section Material given) velocity
     !---------------------------------------------------
     IF ( ConvectionFlag == 'constant' ) THEN
        Velo(1,1:N) = GetReal( Material, 'Convection Velocity 1', Found, Element )
        IF (.NOT.Found) Velo(1,1:N) = 0.0d0
        Velo(2,1:N) = GetReal( Material, 'Convection Velocity 2', Found, Element )
        IF (.NOT.Found) Velo(2,1:N) = 0.0d0
        Velo(3,1:N) = GetReal( Material, 'Convection Velocity 3', Found, Element )
        IF (.NOT.Found) Velo(3,1:N) = 0.0d0
        ! computed velocity
        !------------------Equation => GetEquation()
     ELSE IF (ConvectionFlag == 'computed' ) THEN
        FlowSolName =  GetString( Equation,'Flow Solution Name', Found)
        IF(.NOT.Found) THEN        
           CALL WARN(SolverName,'Keyword >Flow Solution Name< not found in section >Equation<')
           CALL WARN(SolverName,'Taking default value >Flow Solution<')
           WRITE(FlowSolName,'(A)') 'Flow Solution'
        END IF
        FlowSol => VariableGet( Solver % Mesh % Variables, FlowSolName )
        IF ( ASSOCIATED( FlowSol ) ) THEN
           FlowPerm     => FlowSol % Perm
           FlowDOFs     =  FlowSol % DOFs
           FlowSolution => FlowSol % Values
        ELSE
           WRITE(Message,'(A,A,A)') &
                'Convection flag set to >computed<, but no variable >',FlowSolName,'< found'
           CALL FATAL(SolverName,Message)              
        END IF


        DO i=1,n
           k = FlowPerm(Element % NodeIndexes(i))
           IF ( k > 0 ) THEN
              ! Pressure = FlowSolution(FlowDOFs*k)
              
              SELECT CASE( FlowDOFs )
              CASE(2)
                 Velo(1,i) = FlowSolution( FlowDOFs*k-1 )
                 Velo(2,i) = 0.0d0
                 Velo(3,i) = 0.0d0
              CASE(3)
                 Velo(1,i) = FlowSolution( FlowDOFs*k-2 )
                 Velo(2,i) = FlowSolution( FlowDOFs*k-1 )
                 Velo(3,i) = 0.0d0
              CASE(4)
                 Velo(1,i) = FlowSolution( FlowDOFs*k-3 )
                 Velo(2,i) = FlowSolution( FlowDOFs*k-2 )
                 Velo(3,i) = FlowSolution( FlowDOFs*k-1 )
              END SELECT
           END IF
        END DO
     ELSE IF (ConvectionFlag == 'none' ) THEN
        Velo = 0.0d0
     ELSE  
        WRITE(Message,'(A,A,A)') 'Convection flag >', ConvectionFlag ,'< not recognised'
        CALL FATAL(SolverName,Message) 
     END IF
     !-------------------------------------------------
     ! Add mesh velocity (if any) for ALE formulation
     !-------------------------------------------------
     MeshVelo = 0.0d0
     CALL GetVectorLocalSolution( MeshVelo, 'Mesh Velocity', Element)
     IF (ANY( MeshVelo /= 0.0d0 )) THEN
        DO i=1,FlowDOFs
           Velo(i,1:N) = Velo(i,1:N) - MeshVelo(i,1:N)
        END DO
     END IF
   END SUBROUTINE GetLocalALEVelocity

!------------------------------------------------------------------------------
 END SUBROUTINE AdvectionReactionSolver
!------------------------------------------------------------------------------
