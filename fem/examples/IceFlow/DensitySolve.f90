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
! * Module containing a solver for fabric parameter equations  
! *
! ******************************************************************************
! *
! *                     Author:       Juha Ruokolainen
! *
! *                    Address: CSC - IT Center for Science Ltd.
! *                            Keilaranta 14, P.O. Box 405
! *                              02101 Espoo, Finland
! *                              Tel. +358 0 457 2723
! *                            Telefax: +358 0 457 2302
! *                          EMail: Juha.Ruokolainen@csc.fi
! *
! *                       Date: 08 Jun 1997
! *
! *                Modified by: Olivier Gagliardini
! *
! *       Date of modification: avril 2004  
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE DensitySolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

      USE DefUtils

    IMPLICIT NONE
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve dating equation for one timestep (da/dt + a.c0 + da/dxi.ui = 1)
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations (NOTE: Not used
!            currently)
!
!******************************************************************************

     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver

     LOGICAL ::  TransientSimulation
     REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Solver_t), POINTER :: PSolver

     TYPE(Matrix_t),POINTER :: StiffMatrix

     INTEGER :: n1,n2,i,j,k,l,n,t,iter,NDeg,STDOFs,LocalNodes,istat

     TYPE(ValueList_t),POINTER :: Material, BC
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement, Element, &
              ParentElement, LeftParent, RightParent, Edge

     REAL(KIND=dp) :: RelativeChange,UNorm,PrevUNorm,Gravity(3), &
         Tdiff,Normal(3),NewtonTol,NonlinearTol,s

     INTEGER :: NewtonIter,NonlinearIter

     TYPE(Variable_t), POINTER :: DensSol,  DensVariable, FlowVariable, MeshVeloVariable

     REAL(KIND=dp), POINTER :: DensValues(:), FlowValues(:), Solution(:), Ref(:),  MeshVeloValues(:)
           
     CHARACTER(LEN=MAX_NAME_LEN)  ::  FlowSolverName

     INTEGER, POINTER :: DensPerm(:),NodeIndexes(:), FlowPerm(:),  MeshVeloPerm(:)
                        

     LOGICAL :: GotForceBC,GotIt,NewtonLinearization = .FALSE.

     INTEGER :: body_id,bf_id,eq_id, Indexes(128),dim
!
     INTEGER :: old_body = -1

     LOGICAL :: AllocationsDone = .FALSE., FreeSurface

     TYPE(Variable_t), POINTER :: TimeVar

     REAL(KIND=dp),ALLOCATABLE:: MASS(:,:), STIFF(:,:),&
       LOAD(:,:), Force(:), Alpha(:,:), Beta(:), &
          Velocity(:,:), MeshVelocity(:,:)

     REAL(KIND=dp) :: Bu, Bv, Bw, RM(3,3), SaveTime = -1
     REAL(KIND=dp), POINTER :: PrevDens(:), CurrDens(:)

     REAL(KIND=dp) :: at, at0, CPUTime, RealTime, MinD

     SAVE MASS, STIFF, LOAD, Force, ElementNodes, & 
         Alpha, Beta, AllocationsDone, Velocity, MeshVelocity, &
         old_body, PrevDens, CurrDens
       

!------------------------------------------------------------------------------
!  Read constants from constants section of SIF file
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
           
     Solution => Solver % Variable % Values
     STDOFs   =  Solver % Variable % DOFs

     DensSol   => VariableGet( Solver % Mesh % Variables, 'Relative Density' )
     DensPerm   => DensSol % Perm
     DensValues => DensSol % Values

!------------------------------------------------------------------------
!    Get Velocity from (Navier)-Stokes-Solver or AIFlow Solver
!------------------------------------------------------------------------
     FlowSolverName = GetString( Solver % Values, 'Flow Solver Name', GotIt )    
     IF (.NOT.Gotit) FlowSolverName = 'flow solution'
     FlowVariable => VariableGet( Solver % Mesh % Variables, FlowSolverName )

     IF ( ASSOCIATED( FlowVariable ) ) THEN
       FlowPerm    => FlowVariable % Perm
       FlowValues  => FlowVariable % Values
     ELSE
       CALL Info('DensitySolver', 'No variable for velocity associated.', Level=4)
     END IF
!!!!! Mesh Velo
     MeshVeloVariable => VariableGet( Solver % Mesh % Variables, &
            'Mesh Velocity' )

     IF ( ASSOCIATED( MeshVeloVariable ) ) THEN
       MeshVeloPerm    => MeshVeloVariable % Perm
       MeshVeloValues  => MeshVeloVariable % Values
     END IF
      
! Read Minimum Value for Relative Density
     MinD = ListGetConstReal( Solver % Values, 'Minimum Density', GotIt )    
     IF (.NOT.Gotit) MinD = 0.1_dp                     

     StiffMatrix => Solver % Matrix
     UNorm = Solver % Variable % Norm
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
        N = Model % MaxElementNodes
        
        IF ( AllocationsDone ) THEN
           DEALLOCATE (Force,     &
                Velocity, &
                MeshVelocity, &
                MASS,&
                STIFF,&
                LOAD,&
                Alpha,&
                Beta,&
                CurrDens)
        END IF

        ALLOCATE( Force( 2*STDOFs*N ), & 
             Velocity(4, N ),&
             MeshVelocity(3,N), &
             MASS( 2*STDOFs*N,2*STDOFs*N ),  &
             STIFF( 2*STDOFs*N,2*STDOFs*N ),  &
             LOAD( 4,N ),&
             Alpha( 3,N ),&
             Beta( N ),&
             CurrDens( SIZE(Solver % Variable % Values)),&
             STAT=istat )
        IF ( istat /= 0 ) THEN
           CALL Fatal( 'DensitySolve', 'Memory allocation error.' )
        END IF

        CurrDens = 0.0_dp

        IF ( TransientSimulation ) THEN
           IF ( AllocationsDone ) DEALLOCATE(PrevDens)             
           ALLOCATE( PrevDens( SIZE(Solver % Variable % Values)) )
           PrevDens = 0.0_dp
        END IF

        DO i=1,Solver % NumberOFActiveElements
           CurrentElement => GetActiveElement(i)   
           NodeIndexes => CurrentElement % NodeIndexes
           n = GetElementDOFs( Indexes )
           Indexes(1:n) = Solver % Variable % Perm( Indexes(1:n) )
           IF ( TransientSimulation ) THEN
              PrevDens(Indexes(1:n)) = &
                   DensValues(DensPerm(NodeIndexes))
           END IF
           CurrDens(Indexes(1:n)) = &
                DensValues(DensPerm(NodeIndexes))
        END DO

        AllocationsDone = .TRUE.
    
      END IF
!-----------------------------------------------------
        IF( TransientSimulation ) THEN
           TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )
           IF ( SaveTime /= TimeVar % Values(1) ) THEN
              SaveTime = TimeVar % Values(1)
              PrevDens = CurrDens
           END IF
        END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
     NonlinearTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Convergence Tolerance' )

     NewtonTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Newton After Tolerance' )

     NewtonIter = ListGetInteger( Solver % Values, &
        'Nonlinear System Newton After Iterations' )

     NonlinearIter = ListGetInteger( Solver % Values, &
         'Nonlinear System Max Iterations',GotIt )

     IF ( .NOT.GotIt ) NonlinearIter = 1
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

     old_body = -1
!------------------------------------------------------------------------------
     DO iter=1,NonlinearIter
!------------------------------------------------------------------------------
       at  = CPUTime()
       at0 = RealTime()

       CALL Info( 'DensitySolve', ' ', Level=4 )
       CALL Info( 'DensitySolve', ' ', Level=4 )
       CALL Info( 'DensitySolve', '-------------------------------------',Level=4 )
       WRITE( Message, * ) 'Density solver  iteration', iter
       CALL Info( 'DensitySolve', Message,Level=4 )
       CALL Info( 'DensitySolve', '-------------------------------------',Level=4 )
       CALL Info( 'DensitySolve', ' ', Level=4 )
       CALL Info( 'DensitySolve', 'Starting assembly...',Level=4 )

!------------------------------------------------------------------------------


       PrevUNorm = UNorm

       Solver % Variable % Values = CurrDens
       IF ( TransientSimulation ) THEN
          Solver % Variable % PrevValues(:,1) = PrevDens
       END IF

       CALL DefaultInitialize()
!------------------------------------------------------------------------------
       DO t=1,Solver % NumberOFActiveElements
!------------------------------------------------------------------------------

         IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
            (Solver % NumberOfActiveElements-t) / &
               (1.0*Solver % NumberOfActiveElements)), ' % done'
                       
           CALL Info( 'DensitySolve', Message, Level=5 )
           at0 = RealTime()
         END IF

         CurrentElement => GetActiveElement(t)
         CALL GetElementNodes( ElementNodes )
         n = GetElementDOFs( Indexes )
         n = GetElementNOFNodes()
         NodeIndexes => CurrentElement % NodeIndexes

!------------------------------------------------------------------------------
!        Read in material constants from Material section
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
!        Get element local stiffness & mass matrices
!------------------------------------------------------------------------------

         k = FlowVariable % DOFs
         Velocity = 0.0d0
         DO i=1,k-1
            Velocity(i,1:n) = FlowValues(k*(FlowPerm(NodeIndexes)-1)+i)
         END DO

!-------------------mesh velo
         MeshVelocity=0._dp
         IF (ASSOCIATED(MeshVeloVariable)) Then
           k = MeshVeloVariable % DOFs
           DO i=1,k
              MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(NodeIndexes)-1)+i)
           END DO
         EndIF
!--------------------------

         CALL LocalMatrix( MASS, STIFF, FORCE, LOAD, &
              Velocity, MeshVelocity, CurrentElement,  &
              n, ElementNodes)

!------------------------------------------------------------------------------
!        Update global matrices from local matrices 
!------------------------------------------------------------------------------
         IF ( TransientSimulation )  CALL Default1stOrderTime(MASS,STIFF,FORCE)
         CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
      END DO
      CALL Info( 'DensitySolve', 'Assembly done', Level=4 )
!------------------------------------------------------------------------------
   
!------------------------------------------------------------------------------
!     Assembly of the edge terms
!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()
!
! 3D Edge => Faces
      If (dim.eq.3) then
           DO t=1,Solver % Mesh % NumberOfFaces
              Edge => Solver % Mesh % Faces(t)
         IF ( .NOT. ActiveBoundaryElement(Edge) ) CYCLE
       
         LeftParent  => Edge % BoundaryInfo % Left
         RightParent => Edge % BoundaryInfo % Right
         IF ( ASSOCIATED(RightParent) .AND. ASSOCIATED(LeftParent) ) THEN
            n  = GetElementNOFNodes( Edge )
            n1 = GetElementNOFNodes( LeftParent )
            n2 = GetElementNOFNodes( RightParent )

            k = FlowVariable % DOFs
            Velocity = 0.0d0
            DO i=1,k-1
               Velocity(i,1:n) = FlowValues(k*(FlowPerm(Edge % NodeIndexes)-1)+i)
            END DO

   !-------------------mesh velo
         MeshVelocity=0._dp
      IF ( ASSOCIATED( MeshVeloVariable ) ) THEN
            k = MeshVeloVariable % DOFs
            DO i=1,k
               MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(Edge % NodeIndexes)-1)+i)
            END DO
      END IF
    !--------------------------

            FORCE = 0.0d0
            MASS  = 0.0d0
            CALL LocalJumps( STIFF,Edge,n,LeftParent,n1,RightParent,n2,Velocity,MeshVelocity )
            IF ( TransientSimulation )  CALL Default1stOrderTime(MASS, STIFF, FORCE)
            CALL DefaultUpdateEquations( STIFF, FORCE, Edge )
         END IF
      END DO
!
! Else 2D Edges
      ELSE
           DO t=1,Solver % Mesh % NumberOfEdges
              Edge => Solver % Mesh % Edges(t)
         IF ( .NOT. ActiveBoundaryElement(Edge) ) CYCLE
       
         LeftParent  => Edge % BoundaryInfo % Left
         RightParent => Edge % BoundaryInfo % Right
         IF ( ASSOCIATED(RightParent) .AND. ASSOCIATED(LeftParent) ) THEN
            n  = GetElementNOFNodes( Edge )
            n1 = GetElementNOFNodes( LeftParent )
            n2 = GetElementNOFNodes( RightParent )

            k = FlowVariable % DOFs
            Velocity = 0.0d0
            DO i=1,k-1
               Velocity(i,1:n) = FlowValues(k*(FlowPerm(Edge % NodeIndexes)-1)+i)
            END DO

   !-------------------mesh velo
         MeshVelocity=0._dp
       IF ( ASSOCIATED( MeshVeloVariable ) ) THEN
             k = MeshVeloVariable % DOFs
             DO i=1,k
                MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(Edge % NodeIndexes)-1)+i)
             END DO
       END IF
    !--------------------------

            FORCE = 0.0d0
            MASS  = 0.0d0
            CALL LocalJumps( STIFF,Edge,n,LeftParent,n1,RightParent,n2,Velocity,MeshVelocity )
            IF ( TransientSimulation )  CALL Default1stOrderTime(MASS, STIFF, FORCE)
            CALL DefaultUpdateEquations( STIFF, FORCE, Edge )
         END IF
      END DO

      End if


!------------------------------------------------------------------------------
!     Loop over the boundary elements
!------------------------------------------------------------------------------
      DO t = 1, Solver % Mesh % NumberOfBoundaryElements
!------------------------------------------------------------------------------

         Element => GetBoundaryElement(t)
         IF( .NOT. ActiveBoundaryElement() )  CYCLE
         IF( GetElementFamily(Element) == 1 ) CYCLE

         ParentElement => Element % BoundaryInfo % Left
         IF ( .NOT. ASSOCIATED( ParentElement ) ) &
            ParentElement => Element % BoundaryInfo % Right
          
         n  = GetElementNOFNodes( Element )
         n1 = GetElementNOFnodes( ParentElement )
       
         k = FlowVariable % DOFs
         Velocity = 0.0d0
         DO i=1,k-1
            Velocity(i,1:n) = FlowValues(k*(FlowPerm(Element % NodeIndexes)-1)+i)
         END DO
   !-------------------mesh velo
         MeshVelocity=0._dp
       IF ( ASSOCIATED( MeshVeloVariable ) ) THEN
         k = MeshVeloVariable % DOFs
         DO i=1,k
          MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(Element % NodeIndexes)-1)+i)
         End do
       END IF
    !--------------------------

         BC => GetBC()
         LOAD = 0.0d0
         GotIt = .FALSE.
         IF ( ASSOCIATED(BC) ) THEN
            LOAD(1,1:n) = GetReal( BC, 'Relative Density' , GotIt )
         END IF

         MASS = 0.0d0
         CALL LocalMatrixBoundary(  STIFF, FORCE, LOAD(1,1:n), &
                              Element, n, ParentElement, n1, Velocity,MeshVelocity, GotIt )

         IF ( TransientSimulation )  CALL Default1stOrderTime(MASS, STIFF, FORCE)
         CALL DefaultUpdateEquations( STIFF, FORCE )
      END DO

      CALL DefaultFinishAssembly()
!------------------------------------------------------------------------------
      CALL Info( 'DensitySolve', 'Set boundaries done', Level=4 )

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      Unorm = DefaultSolve()
      CurrDens = Solver % Variable % Values
      WRITE(Message,*) 'solve done', minval( solver % variable % values), maxval( Solver % variable % values)
      CALL Info( 'DensitySolve', Message, Level=4 )

      n1 = Solver % Mesh % NumberOfNodes
      ALLOCATE( Ref(n1) )
      Ref = 0.
      DensValues = 0.
      DO t=1,Solver % NumberOfActiveElements
         Element => GetActiveElement(t) 
         n = GetElementDOFs( Indexes )
         n = GetElementNOFNodes()
         DO i=1,n
            k = Element % NodeIndexes(i)
            DensValues( DensPerm(k) ) =    & 
               DensValues( DensPerm(k) ) + &
                  Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) )
           Ref(k) = Ref(k) + 1
        END DO
      END DO

      DO i=1,n1
         j = DensPerm(i)
         IF (j < 1) CYCLE
         IF ( Ref(i) > 0 ) THEN
            DensValues( DensPerm(i) ) = &
                DensValues( DensPerm(i) ) / Ref(i)
            DensValues(DensPerm(i)) = MAX (MIN(DensValues(DensPerm(i)), 1.0) , MinD)
         END IF
      END DO

      write(*,*)'Min, Max :',MINVAL(DensValues),MAXVAL(DensValues)
      
      DEALLOCATE( Ref )

      Unorm = SQRT( SUM( DensValues**2 ) / SIZE(DensValues) )

      IF ( PrevUNorm + UNorm /= 0.0d0 ) THEN
         RelativeChange = 2.0d0 * ABS( PrevUNorm - UNorm) / ( PrevUnorm + UNorm)
      ELSE
         RelativeChange = 0.0d0
      END IF

      WRITE( Message, * ) 'Result Norm   : ',UNorm
      CALL Info( 'DensitySolve', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'DensitySolve', Message, Level=4 )
!------------------------------------------------------------------------------
      IF ( RelativeChange < NewtonTol .OR. &
            iter > NewtonIter ) NewtonLinearization = .TRUE.

      IF ( RelativeChange < NonLinearTol ) EXIT
!------------------------------------------------------------------------------
    END DO ! of nonlinear iter
!------------------------------------------------------------------------------

CONTAINS


!------------------------------------------------------------------------------
   SUBROUTINE LocalMatrix(  MASS, STIFF, FORCE, LOAD, &
                            NodalVelo, NodMeshVelo, Element, &
                            n, Nodes)
          
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: STIFF(:,:),MASS(:,:)
     REAL(KIND=dp) :: LOAD(:,:), NodalVelo(:,:), NodMeshVelo(:,:)
     REAL(KIND=dp), DIMENSION(:) :: FORCE

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element

     INTEGER :: n
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(2*n),ddBasisddx(1,1,1)
     REAL(KIND=dp) :: dBasisdx(2*n,3),SqrtElementMetric

     REAL(KIND=dp) :: A,M, hK,tau,pe1,pe2,unorm,C0, SU(n), SW(n)
     REAL(KIND=dp) :: LoadAtIp,  Grad(3,3)

     INTEGER :: i,j,k,p,q,t,dim,NBasis,ind(3), DOFs = 1

     REAL(KIND=dp) :: s,u,v,w, Radius, B(6,3), G(3,6)
     REAL(KIND=dp) :: Velo(3), DivVelo
  
     INTEGER :: N_Integ
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat, CSymmetry,Fond
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()

     FORCE = 0.0D0
     MASS  = 0.0D0
     STIFF = 0.0D0
!    
!    Integration stuff:
!    ------------------
     NBasis = n
     IntegStuff = GaussPoints( Element  )

          
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

     hk = ElementDiameter( Element, Nodes )
!
!   Now we start integrating:
!   -------------------------
    DO t=1,N_Integ

      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)

!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
               Basis, dBasisdx, ddBasisddx, .FALSE. )

      s = SqrtElementMetric * S_Integ(t)
!------------------------------------------------------------------------------
!
      Velo = 0.0d0
      DO i=1,dim
         Velo(i) = SUM( Basis(1:n) * (NodalVelo(i,1:n)-NodMeshVelo(i,1:n)) )
      END DO
!
!     Reaction coefficient:
!     ---------------------
!        DivVelo = 0.0_dp 
!        DO i=1,dim
!           DivVelo = DivVelo + SUM(NodalVelo(i,1:n)*dBasisdx(1:n,i))
!        END DO
!        C0 = - DivVelo           


!     Loop over basis functions (of both unknowns and weights):
!     ---------------------------------------------------------
      DO p=1,NBasis
         DO q=1,NBasis
            A = 0.0d0
            M = Basis(p) * Basis(q) 
!
!           Reaction terms:
!
! No reaction term (it desapears from integrating by part int(ui dD/dxi Phi dV )
!
!           ---------------
!           A = A - C0 * Basis(q) * Basis(p)

            !
            ! Advection terms:
            ! ----------------
            DO j=1,dim
               A = A - Velo(j) * Basis(q) * dBasisdx(p,j) 
            END DO

!           Add nodal matrix to element matrix:
!           -----------------------------------
            MASS( p,q )  = MASS( p,q )  + s * M
            STIFF( p,q ) = STIFF( p,q ) + s * A
         END DO

!
!        The righthand side...:
!        ----------------------
         LoadAtIp = 0.0                        
         
         FORCE(p) = FORCE(p) + s*LoadAtIp

      END DO
   END DO 
!------------------------------------------------------------------------------
 END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    SUBROUTINE FindParentUVW( Edge, nEdge, Parent, nParent, U, V, W, Basis )
!------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(Element_t), POINTER :: Edge, Parent
      INTEGER :: nEdge, nParent
      REAL( KIND=dp ) :: U, V, W, Basis(:)
!------------------------------------------------------------------------------
      INTEGER :: i, j,l
      REAL(KIND=dp) :: NodalParentU(nEdge),NodalParentV(nEdge),NodalParentW(nEdge)
!------------------------------------------------------------------------------
      DO i = 1,nEdge
        DO j = 1,nParent
          IF ( Edge % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
            NodalParentU(i) = Parent % Type % NodeU(j)
            NodalParentV(i) = Parent % Type % NodeV(j)
            NodalParentW(i) = Parent % Type % NodeW(j)
            EXIT
          END IF
        END DO
      END DO
      U = SUM( Basis(1:nEdge) * NodalParentU(1:nEdge) )
      V = SUM( Basis(1:nEdge) * NodalParentV(1:nEdge) )
      W = SUM( Basis(1:nEdge) * NodalParentW(1:nEdge) )
!------------------------------------------------------------------------------      
    END SUBROUTINE FindParentUVW
!------------------------------------------------------------------------------      


!------------------------------------------------------------------------------
    SUBROUTINE LocalJumps( STIFF,Edge,n,LeftParent,n1,RightParent,n2,Velo,MeshVelo )
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: STIFF(:,:), Velo(:,:),MeshVelo(:,:)
      INTEGER :: n,n1,n2
      TYPE(Element_t), POINTER :: Edge, LeftParent, RightParent
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: EdgeBasis(n), EdgedBasisdx(n,3), EdgeddBasisddx(n,3,3)
      REAL(KIND=dp) :: LeftBasis(n1), LeftdBasisdx(n1,3), LeftddBasisddx(n1,3,3)
      REAL(KIND=dp) :: RightBasis(n2), RightdBasisdx(n2,3), RightddBasisddx(n2,3,3)
      REAL(KIND=dp) :: Jump(n1+n2), Average(n1+n2)
      REAL(KIND=dp) :: detJ, U, V, W, S, Udotn, xx, yy
      LOGICAL :: Stat
      INTEGER :: i, j, p, q, dim, t, nEdge, nParent
      TYPE(GaussIntegrationPoints_t) :: IntegStuff
      REAL(KIND=dp) :: hE, Normal(3), cu(3), LeftOut(3)

      TYPE(Nodes_t) :: EdgeNodes, LeftParentNodes, RightParentNodes

      SAVE EdgeNodes, LeftParentNodes, RightParentNodes

!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()
      STIFF = 0.0d0

      CALL GetElementNodes( EdgeNodes, Edge )
      CALL GetElementNodes( LeftParentNodes,  LeftParent )
      CALL GetElementNodes( RightParentNodes, RightParent )
!------------------------------------------------------------------------------
!     Numerical integration over the edge
!------------------------------------------------------------------------------
      IntegStuff = GaussPoints( Edge )

      LeftOut(1) = SUM( LeftParentNodes % x(1:n1) ) / n1
      LeftOut(2) = SUM( LeftParentNodes % y(1:n1) ) / n1
      LeftOut(3) = SUM( LeftParentNodes % z(1:n1) ) / n1
      LeftOut(1) = SUM( EdgeNodes % x(1:n) ) / n - LeftOut(1)
      LeftOut(2) = SUM( EdgeNodes % y(1:n) ) / n - LeftOut(2)
      LeftOut(3) = SUM( EdgeNodes % z(1:n) ) / n - LeftOut(3)

      DO t=1,IntegStuff % n
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)
        S = IntegStuff % s(t)

        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Edge, EdgeNodes, U, V, W, detJ, &
             EdgeBasis, EdgedBasisdx, EdgeddBasisddx, .FALSE. )

        S = S * detJ

        Normal = NormalVector( Edge, EdgeNodes, U, V, .FALSE. )
        IF ( SUM( LeftOut*Normal ) < 0 ) Normal = -Normal

        ! Find basis functions for the parent elements:
        ! ---------------------------------------------
        CALL FindParentUVW( Edge,n,LeftParent,n1,U,V,W,EdgeBasis )
        stat = ElementInfo( LeftParent, LeftParentNodes, U, V, W, detJ, &
                LeftBasis, LeftdBasisdx, LeftddBasisddx, .FALSE. )

        CALL FindParentUVW( Edge,n,RightParent,n2,U,V,W,EdgeBasis )
        stat = ElementInfo( RightParent, RightParentNodes, U, V, W, detJ, &
              RightBasis, RightdBasisdx, RightddBasisddx, .FALSE. )

        ! Integrate jump terms:
        ! ---------------------
        Jump(1:n1) = LeftBasis(1:n1)
        Jump(n1+1:n1+n2) = -RightBasis(1:n2)

        Average(1:n1) = LeftBasis(1:n1) / 2
        Average(n1+1:n1+n2) = RightBasis(1:n2) / 2

        cu = 0.0d0
        DO i=1,dim
          cu(i) = SUM( (Velo(i,1:n)-MeshVelo(i,1:n)) * EdgeBasis(1:n) )
        END DO
        Udotn = SUM( Normal * cu )

        DO p=1,n1+n2
          DO q=1,n1+n2
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
        Element, n, ParentElement, np, Velo,MeshVelo, InFlowBC )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: STIFF(:,:),  FORCE(:), LOAD(:), Velo(:,:),MeshVelo(:,:)
     INTEGER :: n, np
     LOGICAL :: InFlowBC
     TYPE(Element_t), POINTER :: Element, ParentElement
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3)
     REAL(KIND=dp) :: ParentBasis(np), ParentdBasisdx(np,3), ParentddBasisddx(np,3,3)
     INTEGER :: i,j,p,q,t,dim

     REAL(KIND=dp) :: Normal(3), g, L, Udotn,  cu(3), detJ,U,V,W,S
     LOGICAL :: Stat, Inflow
     TYPE(GaussIntegrationPoints_t) :: IntegStuff

     TYPE(Nodes_t) :: Nodes, ParentNodes
     SAVE Nodes, ParentNodes
!------------------------------------------------------------------------------
     dim = CoordinateSystemDimension()
     FORCE = 0.0d0
     STIFF = 0.0d0

     CALL GetElementNodes( Nodes, Element )
     CALL GetElementNodes( ParentNodes, ParentElement )

     Normal = NormalVector( Element, Nodes, 0.0d0, 0.0d0, .TRUE. )
     DO i=1,3
       cu(i) = SUM( Velo(i,1:n) ) / n
     END DO
     Inflow = InFlowBC .AND. SUM( Normal * cu ) < 0.0d0


     ! Numerical integration:
     !-----------------------
     IntegStuff = GaussPoints( Element )


     DO t=1,IntegStuff % n
       U = IntegStuff % u(t)
       V = IntegStuff % v(t)
       W = IntegStuff % w(t)
       S = IntegStuff % s(t)

       Normal = NormalVector( Element, Nodes, U, V, .TRUE. ) 

       ! Basis function values & derivatives at the integration point:
       ! -------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
               Basis, dBasisdx, ddBasisddx, .FALSE. )
       S = S * detJ

       CALL FindParentUVW( Element, n, ParentElement, np, U, V, W, Basis )
       stat = ElementInfo( ParentElement, ParentNodes, U, V, W, &
           detJ,  ParentBasis, ParentdBasisdx, ParentddBasisddx, .FALSE. )

       L = SUM( LOAD(1:n) * Basis(1:n) )
       cu = 0.0d0
       DO i=1,dim
          cu(i) = SUM( (Velo(i,1:n)-MeshVelo(i,1:n)) * Basis(1:n) )
       END DO
       Udotn = SUM( Normal * cu )


       DO p = 1,np
        IF ( InFlow ) THEN
            FORCE(p) = FORCE(p) - s * Udotn*L*ParentBasis(p)
         ELSE
           DO q=1,np
             STIFF(p,q) = STIFF(p,q) + s*Udotn*ParentBasis(q)*ParentBasis(p)
           END DO
         END IF
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 END SUBROUTINE DensitySolver
!------------------------------------------------------------------------------
