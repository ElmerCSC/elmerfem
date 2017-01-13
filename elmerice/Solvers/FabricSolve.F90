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
! *  Authors:  Juha Ruokolainen, Olivier Gagliardini, Fabien Gillet-Chaulet
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://elmerice.elmerfem.org
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *       Date of modification: 03/05
! *
! *****************************************************************************/
!>  Solver for fabric parameter equations 
!------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE FabricSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

      USE DefUtils

      IMPLICIT NONE
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve stress equations for one timestep
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

     INTEGER :: dim,n1,n2,i,j,k,l,n,t,iter,NDeg,STDOFs,LocalNodes,istat

     TYPE(ValueList_t),POINTER :: Material, BC
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement, Element, &
              ParentElement, LeftParent, RightParent, Edge

     REAL(KIND=dp) :: RelativeChange,UNorm,PrevUNorm,Gravity(3), &
         Tdiff,Normal(3),NewtonTol,NonlinearTol,s,Wn(7)


     INTEGER :: NewtonIter,NonlinearIter

     TYPE(Variable_t), POINTER :: FabricSol, TempSol, FabricVariable, FlowVariable, &
                                  EigenFabricVariable, MeshVeloVariable

     REAL(KIND=dp), POINTER :: Temperature(:),Fabric(:), &
           FabricValues(:), FlowValues(:), EigenFabricValues(:), &
           MeshVeloValues(:), Solution(:), Ref(:)

     INTEGER, POINTER :: TempPerm(:),FabricPerm(:),NodeIndexes(:), &
                        FlowPerm(:), MeshVeloPerm(:), EigenFabricPerm(:)

     REAL(KIND=dp) :: rho,lambda   !Interaction parameter,diffusion parameter
     REAL(KIND=dp) :: A1plusA2
     Real(KIND=dp), parameter :: Rad2deg=180._dp/Pi
     REAL(KIND=dp) :: a2(6)
     REAL(KIND=dp) :: ai(3), Angle(3)

     LOGICAL :: GotForceBC,GotIt,NewtonLinearization = .FALSE.,UnFoundFatal=.TRUE.

     INTEGER :: body_id,bf_id,eq_id, comp, Indexes(128)
!
     INTEGER :: old_body = -1

     REAL(KIND=dp) :: FabricGrid(4879)                   
                        
     LOGICAL :: AllocationsDone = .FALSE., FreeSurface

     TYPE(Variable_t), POINTER :: TimeVar

     REAL(KIND=dp), ALLOCATABLE:: MASS(:,:), STIFF(:,:),  &
       LocalFluidity(:), LOAD(:,:),Force(:), LocalTemperature(:), &
       Alpha(:,:),Beta(:), K1(:), K2(:), E1(:), E2(:), E3(:), &
       Velocity(:,:), MeshVelocity(:,:)

     SAVE MASS, STIFF, LOAD, Force,ElementNodes,Alpha,Beta, & 
          LocalTemperature, LocalFluidity,  AllocationsDone, K1, K2, &
          E1, E2, E3, Wn,  FabricGrid, rho, lambda, Velocity, &
          MeshVelocity, old_body, dim, comp
!------------------------------------------------------------------------------
     CHARACTER(LEN=MAX_NAME_LEN) :: viscosityFile

     REAL(KIND=dp) :: Bu,Bv,Bw,RM(3,3), SaveTime = -1
     REAL(KIND=dp), POINTER :: PrevFabric(:),CurrFabric(:),TempFabVal(:)

     SAVE  ViscosityFile, PrevFabric, CurrFabric,TempFabVal
#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at, at0
#else
     REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif
!------------------------------------------------------------------------------
      INTERFACE
        Subroutine R2Ro(a2,dim,ai,angle)
        USE Types
        REAL(KIND=dp),intent(in) :: a2(6)
        Integer :: dim
        REAL(KIND=dp),intent(out) :: ai(3), Angle(3)
       End Subroutine R2Ro
      End Interface                                                       
!------------------------------------------------------------------------------
!  Read constants from constants section of SIF file
!------------------------------------------------------------------------------
         
      Wn(7) = ListGetConstReal( Model % Constants, 'Gas Constant', GotIt,UnFoundFatal=UnFoundFatal )
      !Previous default value: Wn(7) = 8.314
      WRITE(Message,'(A,F10.4)')'Gas Constant =',Wn(7)
      CALL INFO('FabricSolve',Message,Level=4)
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
      IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

      Solution => Solver % Variable % Values
      STDOFs   =  Solver % Variable % DOFs

      FabricSol    => VariableGet( Solver % Mesh % Variables, 'Fabric' )
      FabricPerm   => FabricSol % Perm
      FabricValues => FabricSol % Values

      TempSol => VariableGet( Solver % Mesh % Variables, 'Temperature' )
      IF ( ASSOCIATED( TempSol) ) THEN
       TempPerm    => TempSol % Perm
       Temperature => TempSol % Values
      END IF

      FlowVariable => VariableGet( Solver % Mesh % Variables, 'AIFlow' )
      IF ( ASSOCIATED( FlowVariable ) ) THEN
       FlowPerm    => FlowVariable % Perm
       FlowValues  => FlowVariable % Values
      END IF
      
!!!!! Mesh Velo
     MeshVeloVariable => VariableGet( Solver % Mesh % Variables, &
            'Mesh Velocity' )

     IF ( ASSOCIATED( MeshVeloVariable ) ) THEN
       MeshVeloPerm    => MeshVeloVariable % Perm
       MeshVeloValues  => MeshVeloVariable % Values
     END IF
       
                                       
      StiffMatrix => Solver % Matrix
      ! UNorm = Solver % Variable % Norm
      Unorm = SQRT( SUM( FabricValues**2 ) / SIZE(FabricValues) )
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
      IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
        N = Model % MaxElementNodes

       dim = CoordinateSystemDimension()
       
       IF ( AllocationsDone ) THEN
         DEALLOCATE( LocalTemperature, &
                     K1,K2,E1,E2,E3, &
                     Force, LocalFluidity, &
                     Velocity,MeshVelocity, &
                     MASS,STIFF,      &
                     LOAD, Alpha, Beta, &
                     CurrFabric, TempFabVal )
       END IF

       ALLOCATE( LocalTemperature( N ), LocalFluidity( N ), &
                 K1( N ), K2( N ), E1( N ), E2( N ), E3( N ), &
                 Force( 2*STDOFs*N ), &
                 Velocity(4, N ),MeshVelocity(3,N), &
                 MASS( 2*STDOFs*N,2*STDOFs*N ),  &
                 STIFF( 2*STDOFs*N,2*STDOFs*N ),  &
                 LOAD( 4,N ), Alpha( 3,N ), Beta( N ), &
                 CurrFabric( 5*SIZE(Solver % Variable % Values)), &
                 TempFabVal( 5*SIZE(Solver % Variable % Values)), &
                 STAT=istat )


       IF ( istat /= 0 ) THEN
          CALL Fatal( 'FabricSolve', 'Memory allocation error.' )
       END IF

       CurrFabric = 0.
       TempFabVal = 0.
       IF ( TransientSimulation ) THEN
          IF (AllocationsDone ) DEALLOCATE (PrevFabric)
          ALLOCATE( PrevFabric( 5*SIZE(Solver % Variable % Values)) )
         PrevFabric = 0.
       END IF

       DO i=1,Solver % NumberOFActiveElements
          CurrentElement => GetActiveElement(i)   
          n = GetElementDOFs( Indexes )
          n = GetElementNOFNodes()
          NodeIndexes => CurrentElement % NodeIndexes
          Indexes(1:n) = Solver % Variable % Perm( Indexes(1:n) )
          DO COMP=1,5
            IF ( TransientSimulation ) THEN
               PrevFabric(5*(Indexes(1:n)-1)+COMP) = &
                   FabricValues(5*(FabricPerm(NodeIndexes(1:n))-1)+COMP)
            END IF
               CurrFabric(5*(Indexes(1:n)-1)+COMP) = &
                   FabricValues(5*(FabricPerm(NodeIndexes(1:n))-1)+COMP)
          END DO
       END DO

       AllocationsDone = .TRUE.
      END IF

      IF( TransientSimulation ) THEN
        TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )
        IF ( SaveTime /= TimeVar % Values(1) ) THEN
           SaveTime = TimeVar % Values(1)
           PrevFabric = CurrFabric
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

!------------------------------------------------------------------------------
      DO iter=1,NonlinearIter
!------------------------------------------------------------------------------
       at  = CPUTime()
       at0 = RealTime()


       CALL Info( 'FabricSolve', ' ', Level=4 )
       CALL Info( 'FabricSolve', ' ', Level=4 )
       CALL Info( 'FabricSolve', &
                    '-------------------------------------',Level=4 )
       WRITE( Message, * ) 'Fabric solver  iteration', iter
       CALL Info( 'FabricSolve', Message,Level=4 )
       CALL Info( 'FabricSolve', &
                     '-------------------------------------',Level=4 )
       CALL Info( 'FabricSolve', ' ', Level=4 )
       CALL Info( 'FabricSolve', 'Starting assembly...',Level=4 )

!------------------------------------------------------------------------------

       PrevUNorm = UNorm
       
       DO COMP=1,2*dim-1

       Solver % Variable % Values = CurrFabric( COMP::5 )
       IF ( TransientSimulation ) THEN
          Solver % Variable % PrevValues(:,1) = PrevFabric( COMP::5 )
       END IF

       CALL DefaultInitialize()
!------------------------------------------------------------------------------
       DO t=1,Solver % NumberOFActiveElements
!------------------------------------------------------------------------------

         IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
            (Solver % NumberOfActiveElements-t) / &
               (1.0*Solver % NumberOfActiveElements)), ' % done'
                       
           CALL Info( 'FabricSolve', Message, Level=5 )
           at0 = RealTime()
         END IF

         CurrentElement => GetActiveElement(t)
         CALL GetElementNodes( ElementNodes )
         n = GetElementDOFs( Indexes )
         n = GetElementNOFNodes()
         NodeIndexes => CurrentElement % NodeIndexes

         
         Material => GetMaterial()
         body_id = CurrentElement % BodyId
!------------------------------------------------------------------------------
!        Read in material constants from Material section
!------------------------------------------------------------------------------
         IF (body_id /= old_body) Then 
           old_body = body_id
           CALL GetMaterialDefs()
         END IF
      
         LocalFluidity(1:n) = ListGetReal( Material, &
                         'Fluidity Parameter', n, NodeIndexes, GotIt,&
                         UnFoundFatal=UnFoundFatal)
         !Previous default value: LocalFluidity(1:n) = 1.0
!------------------------------------------------------------------------------
!        Get element local stiffness & mass matrices
!------------------------------------------------------------------------------
         LocalTemperature = 0.0D0
         IF ( ASSOCIATED(TempSol) ) THEN
            DO i=1,n
               k = TempPerm(NodeIndexes(i))
               LocalTemperature(i) = Temperature(k)
            END DO
         ELSE
            LocalTemperature(1:n) = 0.0d0
         END IF

         K1(1:n) = CurrFabric( 5*(Solver % Variable % Perm(Indexes(1:n))-1)+1 )
         K2(1:n) = CurrFabric( 5*(Solver % Variable % Perm(Indexes(1:n))-1)+2 )
         E1(1:n) = CurrFabric( 5*(Solver % Variable % Perm(Indexes(1:n))-1)+3 )
         E2(1:n) = CurrFabric( 5*(Solver % Variable % Perm(Indexes(1:n))-1)+4 )
         E3(1:n) = CurrFabric( 5*(Solver % Variable % Perm(Indexes(1:n))-1)+5 )

         k = FlowVariable % DOFs
         Velocity = 0.0d0
         DO i=1,k
            Velocity(i,1:n) = FlowValues(k*(FlowPerm(NodeIndexes)-1)+i)
         END DO

!------------meshvelocity
         MeshVelocity=0._dp
         IF (ASSOCIATED(MeshVeloVariable)) Then
           k = MeshVeloVariable % DOFs
           DO i=1,k
              MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(NodeIndexes)-1)+i)
           END DO
         EndIF
!----------------------------------

         CALL LocalMatrix( COMP, MASS, STIFF, FORCE, LOAD, K1, K2, E1, &
           E2, E3, LocalTemperature, LocalFluidity,  Velocity, &
           MeshVelocity, CurrentElement, n, ElementNodes, Wn, rho, lambda )

!------------------------------------------------------------------------------
!        Update global matrices from local matrices 
!------------------------------------------------------------------------------
         IF ( TransientSimulation )  CALL Default1stOrderTime(MASS,STIFF,FORCE)
         CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
      END DO
      CALL Info( 'FabricSolve', 'Assembly done', Level=4 )
!------------------------------------------------------------------------------
   
!------------------------------------------------------------------------------
!     Assembly of the edge terms
!------------------------------------------------------------------------------
!
!3D => Edges => Faces
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
            DO i=1,k
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
!  2D
      Else

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
            DO i=1,k
               Velocity(i,1:n) = FlowValues(k*(FlowPerm(Edge % NodeIndexes(1:n))-1)+i)
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

      END IF

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
         DO i=1,k
            Velocity(i,1:n) = FlowValues(k*(FlowPerm(Element % NodeIndexes(1:n))-1)+i)
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
            LOAD(1,1:n) = GetReal( BC, ComponentName('Fabric', Comp) , GotIt )
         END IF

         MASS = 0.0d0
         CALL LocalMatrixBoundary(  STIFF, FORCE, LOAD(1,1:n), &
                              Element, n, ParentElement, n1, Velocity,MeshVelocity, GotIt )

         IF ( TransientSimulation )  CALL Default1stOrderTime(MASS, STIFF, FORCE)
         CALL DefaultUpdateEquations( STIFF, FORCE )
      END DO

      CALL DefaultFinishAssembly()
!------------------------------------------------------------------------------
      CALL Info( 'FabricSolve', 'Set boundaries done', Level=4 )

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      Unorm = DefaultSolve()
!      CurrFabric( COMP::5 ) = Solver % Variable % Values
      WRITE(Message,*) 'solve done', minval( solver % variable % values), maxval( Solver % variable % values)
      CALL Info( 'FabricSolve', Message, Level=4 )
      
      n1 = Solver % Mesh % NumberOfNodes
      ALLOCATE( Ref(n1) )
      Ref = 0
      !
      ! FabricValues( COMP::5 ) = 0 !fab sinon remet toutes les
      ! fabriques à 0 dans le cas ou on a 2 domaines dont l'un a la
      ! fabrique fixe
      TempFabVal(COMP::5 ) = 0. !fab
      
      DO t=1,Solver % NumberOfActiveElements
         Element => GetActiveElement(t) 
         n = GetElementDOFs( Indexes )
         n = GetElementNOFNodes()
         
         DO i=1,n
            k = Element % NodeIndexes(i)
            TempFabVal( 5*(FabricPerm(k)-1) + COMP ) =    & 
            TempFabVal( 5*(FabricPerm(k)-1) + COMP ) + &
            Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) )
            FabricValues( 5*(FabricPerm(k)-1) + COMP ) = &
                          TempFabVal(5*(FabricPerm(k)-1) + COMP ) 
            Ref(k) = Ref(k) + 1
         END DO
      END DO

      DO i=1,n1
         j=FabricPerm(i)
         IF (j < 1) CYCLE
         IF ( Ref(i) > 0 ) THEN
            FabricValues( 5*(j-1)+COMP ) = &
                   FabricValues( 5*(j-1)+COMP ) / Ref(i)
         END IF
      END DO

      DEALLOCATE( Ref )

      SELECT CASE( Comp ) 
      CASE(1)
      FabricValues( COMP:SIZE(FabricValues):5 ) = &
          MIN(MAX( FabricValues( COMP:SIZE(FabricValues):5 ) , 0._dp),1._dp)
          
      CASE(2)
       FabricValues( COMP:SIZE(FabricValues):5 ) = &
       MIN(MAX( FabricValues( COMP:SIZE(FabricValues):5 ) , 0._dp),1._dp)
       
     !  !DO i=1,SIZE(FabricValues),5 
     !  !  IF((FabricValues(i)+FabricValues(i+1)).GT.1._dp) THEN
     !  !      A1plusA2=FabricValues(i)+FabricValues(i+1)
     !  !      FabricValues(i)= &
     !  !        FabricValues(i)/A1plusA2
     !  !      FabricValues(i+1)= &
     !  !        FabricValues(i+1)/A1plusA2
     !  !   END IF
     !  ! END DO
     !   

!      CASE(3:5)
!       DO i=COMP,SIZE(FabricValues),5
!         IF(FabricValues(i).GT.0._dp) THEN
!               FabricValues(i) = MIN( FabricValues(i) , 0.5_dp)
!         ELSE
!               FabricValues(i) = MAX( FabricValues(i) , -0.5_dp)
!         END IF
!       END DO
      END SELECT

      END DO ! End DO Comp

       DO i=1,Solver % NumberOFActiveElements
          CurrentElement => GetActiveElement(i)   
          n = GetElementDOFs( Indexes )
          n = GetElementNOFNodes()
          NodeIndexes => CurrentElement % NodeIndexes
          Indexes(1:n) = Solver % Variable % Perm( Indexes(1:n) )
          DO COMP=1,5
            CurrFabric(5*(Indexes(1:n)-1)+COMP) = &
                        FabricValues(5*(FabricPerm(NodeIndexes(1:n))-1)+COMP)
          END DO
       END DO


      Unorm = SQRT( SUM( FabricValues**2 ) / SIZE(FabricValues) )
      Solver % Variable % Norm = Unorm  

      IF ( PrevUNorm + UNorm /= 0.0d0 ) THEN
         RelativeChange = 2.0d0 * ABS( PrevUNorm - UNorm) / ( PrevUnorm + UNorm)
      ELSE
         RelativeChange = 0.0d0
      END IF

      WRITE( Message, * ) 'Result Norm   : ',UNorm
      CALL Info( 'FabricSolve', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'FabricSolve', Message, Level=4 )

      
!------------------------------------------------------------------------------
      IF ( RelativeChange < NewtonTol .OR. &
            iter > NewtonIter ) NewtonLinearization = .TRUE.

      IF ( RelativeChange < NonLinearTol ) EXIT

!------------------------------------------------------------------------------
      EigenFabricVariable => &
       VariableGet( Solver % Mesh % Variables, 'EigenV' )
     IF ( ASSOCIATED( EigenFabricVariable ) ) THEN
         EigenFabricPerm    => EigenFabricVariable % Perm
         EigenFabricValues => EigenFabricVariable % Values
      
         DO t=1,Solver % NumberOFActiveElements

           CurrentElement => GetActiveElement(t)
           n = GetElementNOFNodes()
           NodeIndexes => CurrentElement % NodeIndexes

           K1(1:n) = FabricValues( 5*(FabricPerm(NodeIndexes(1:n))-1)+1 )
           K2(1:n) = FabricValues( 5*(FabricPerm(NodeIndexes(1:n))-1)+2 )
           E1(1:n) = FabricValues( 5*(FabricPerm(NodeIndexes(1:n))-1)+3 )
           E2(1:n) = FabricValues( 5*(FabricPerm(NodeIndexes(1:n))-1)+4 )
           E3(1:n) = FabricValues( 5*(FabricPerm(NodeIndexes(1:n))-1)+5 )
         
           Do i=1,n
            a2(1)=K1(i)
            a2(2)=K2(i)
            a2(3)=1._dp-a2(1)-a2(2)
            a2(4)=E1(i)
            a2(5)=E2(i)
            a2(6)=E3(i)

            call R2Ro(a2,dim,ai,angle)

            angle(:)=angle(:)*rad2deg
            If (angle(1).gt.90._dp) angle(1)=angle(1)-180._dp
            If (angle(1).lt.-90._dp) angle(1)=angle(1)+180._dp
          
           EigenFabricValues( 6 * (EigenFabricPerm(NodeIndexes(i))-1) + 1)=ai(1)
           EigenFabricValues( 6 * (EigenFabricPerm(NodeIndexes(i))-1) + 2)=ai(2)
           EigenFabricValues( 6 * (EigenFabricPerm(NodeIndexes(i))-1) + 3 )=ai(3)
           EigenFabricValues( 6 * (EigenFabricPerm(NodeIndexes(i))-1) + 4 )=angle(1)
           EigenFabricValues( 6 * (EigenFabricPerm(NodeIndexes(i))-1) + 5 )=angle(2)
           EigenFabricValues( 6 * (EigenFabricPerm(NodeIndexes(i))-1) + 6 )=angle(3)
           End do
        END DO

      END IF
!------------------------------------------------------------------------------
    END DO ! of nonlinear iter
!------------------------------------------------------------------------------
CONTAINS

      SUBROUTINE GetMaterialDefs()

      viscosityFile = ListGetString( Material ,'Viscosity File',GotIt, UnFoundFatal)
      OPEN( 1, File = viscosityFile)
      DO i=1,813
         READ( 1, '(6(e14.8))' ) FabricGrid( 6*(i-1)+1:6*(i-1)+6 )
      END DO
      READ(1 , '(e14.8)' ) FabricGrid(4879)
      CLOSE(1)

       rho = ListGetConstReal(Material, 'Interaction Parameter', GotIt )
       IF (.NOT.GotIt) THEN
           WRITE(Message,'(A)') 'Interaction  Parameter notfound. &
                         &Setting to the value in ViscosityFile'
           CALL INFO('AIFlowSolve', Message, Level = 20)
           rho = FabricGrid(4879)
       ELSE
           WRITE(Message,'(A,F10.4)') 'Interaction Parameter = ', rho
           CALL INFO('AIFlowSolve', Message, Level = 20)
       END IF

       lambda = ListGetConstReal( Material, 'Diffusion Parameter', GotIt,UnFoundFatal=UnFoundFatal)
           !Previous default value: lambda = 0.0_dp
      WRITE(Message,'(A,F10.4)') 'Diffusion Parameter = ', lambda
      CALL INFO('AIFlowSolve', Message, Level = 20)

      Wn(2) = ListGetConstReal( Material , 'Powerlaw Exponent', GotIt,UnFoundFatal=UnFoundFatal)
           !Previous default value: Wn(2) = 1.0
      WRITE(Message,'(A,F10.4)') 'Powerlaw Exponent = ',   Wn(2)
      CALL INFO('AIFlowSolve', Message, Level = 20)

      Wn(3) = ListGetConstReal( Material, 'Activation Energy 1', GotIt,UnFoundFatal=UnFoundFatal)
           !Previous default value: Wn(3) = 1.0
      WRITE(Message,'(A,F10.4)') 'Activation Energy 1 = ',   Wn(3)
      CALL INFO('AIFlowSolve', Message, Level = 20)

      Wn(4) = ListGetConstReal( Material, 'Activation Energy 2', GotIt,UnFoundFatal=UnFoundFatal)
           !Previous default value: Wn(4) = 1.0
      WRITE(Message,'(A,F10.4)') 'Activation Energy 2 = ',   Wn(4)
      CALL INFO('AIFlowSolve', Message, Level = 20)

      Wn(5) = ListGetConstReal(Material, 'Reference Temperature', GotIt,UnFoundFatal=UnFoundFatal)
           !Previous default value: Wn(5) = -10.0
      WRITE(Message,'(A,F10.4)') 'Reference Temperature = ',   Wn(5)
      CALL INFO('AIFlowSolve', Message, Level = 20)

      Wn(6) = ListGetConstReal( Material, 'Limit Temperature', GotIt,UnFoundFatal=UnFoundFatal)
           !Previous default value: Wn(6) = -10.0
      WRITE(Message,'(A,F10.4)') 'Limit Temperature = ',   Wn(6)
      CALL INFO('AIFlowSolve', Message, Level = 20)
!------------------------------------------------------------------------------
      END SUBROUTINE GetMaterialDefs
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE LocalMatrix( Comp, MASS, STIFF, FORCE, LOAD, &
          NodalK1, NodalK2, NodalEuler1, NodalEuler2, NodalEuler3, & 
          NodalTemperature, NodalFluidity, NodalVelo, NodMeshVel, &
          Element, n, Nodes, Wn, rho,lambda )
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: STIFF(:,:),MASS(:,:)
     REAL(KIND=dp) :: LOAD(:,:), NodalVelo(:,:),NodMeshVel(:,:)
     REAL(KIND=dp), DIMENSION(:) :: FORCE, NodalK1, NodalK2, NodalEuler1, &
               NodalEuler2, NodalEuler3, NodalTemperature, NodalFluidity

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
     INTEGER :: n, Comp
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(2*n),ddBasisddx(1,1,1)
     REAL(KIND=dp) :: dBasisdx(2*n,3),SqrtElementMetric

     REAL(KIND=dp) :: A1, A2, A3, E1, E2, E3, Theta

     REAL(KIND=dp) :: A,M, hK,tau,pe1,pe2,unorm,C0, SU(n), SW(n)
     REAL(KIND=dp) :: LoadAtIp, Temperature
     REAL(KIND=dp) :: rho,lambda,Deq, ai(6),a4(9),hmax

     INTEGER :: i,j,k,p,q,t,dim,NBasis,ind(3), DOFs = 1

     REAL(KIND=dp) :: s,u,v,w, Radius, B(6,3), G(3,6), C44,C55,C66
     REAL(KIND=dp) :: Wn(:),Velo(3),DStress(6),StrainR(6),Spin(3),SD(6)

     REAL(KIND=dp) :: LGrad(3,3),StrainRate(3,3),D(6),angle(3),epsi
     REAL(KIND=dp) :: ap(3),C(6,6),Spin1(3,3),Stress(3,3)
     Integer :: INDi(6),INDj(6)
     LOGICAL :: CSymmetry
     
     LOGICAL :: Fond
              
!              
     INTEGER :: N_Integ
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     INTERFACE
        FUNCTION BGlenT( Tc, W)
           USE Types
           REAL(KIND=dp) :: BGlenT,Tc,W(7)
        END FUNCTION
        
        SUBROUTINE IBOF(ai,a4)
           USE Types
           REAL(KIND=dp),intent(in) :: ai(6)
           REAL(KIND=dp),intent(out) :: a4(9)
        END SUBROUTINE

        Subroutine R2Ro(ai,dim,a2,angle)
         USE Types
         REAL(KIND=dp),intent(in) :: ai(6)
         Integer :: dim
         REAL(KIND=dp),intent(out) :: a2(3), Angle(3)
        End Subroutine R2Ro

        Subroutine OPILGGE_ai_nl(a2,Angle,etaI,eta36)
          USE Types
          REAL(kind=dp), INTENT(in),  DIMENSION(3)   :: a2
          REAL(kind=dp), INTENT(in),  DIMENSION(3)   :: Angle
          REAL(kind=dp), INTENT(in),  DIMENSION(:)   :: etaI
          REAL(kind=dp), INTENT(out), DIMENSION(6,6) :: eta36
        END SUBROUTINE OPILGGE_ai_nl
        
      END INTERFACE
!------------------------------------------------------------------------------
      Fond=.False.
      hmax = maxval (Nodes % y(1:n))
        
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
!     Orientation parameters at the integration point:
!     ------------------------------------------------
      A1 = SUM( NodalK1(1:n) * Basis(1:n) ) 
      A2 = SUM( NodalK2(1:n) * Basis(1:n) )
      A3 = 1._dp - A1 - A2

      E1 = SUM( NodalEuler1(1:n) * Basis(1:n) )
      E2 = SUM( NodalEuler2(1:n) * Basis(1:n) )
      E3 = SUM( NodalEuler3(1:n) * Basis(1:n) )
!
!      Fluidity  at the integration point:
!---------------------------------------------
!       Wn(1)=SUM( NodalFluidity(1:n)*Basis(1:n) )
!
!     Temperature at the integration point:
!     -------------------------------------
!      Temperature = SUM( NodalTemperature(1:n)*Basis(1:n) )
!
!     Get theta parameter: (the (fluidity in the basal plane)/2, 
!     function of the Temperature )
!     -----------------------------------------------------
      Theta = 1._dp / ( FabricGrid(5) + FabricGrid(6) )
      Theta = Theta 

!      Strain-Rate, Stresses and Spin

      CSymmetry = CurrentCoordinateSystem() == AxisSymmetric
      
      Stress = 0.0
      StrainRate = 0.0
      Spin1 = 0.0
!
!    Material parameters at that point
!    ---------------------------------
!
      ai(1)=A1
      ai(2)=A2
      ai(3)=A3
      ai(4)=E1
      ai(5)=E2
      ai(6)=E3

!    fourth order orientation tensor
      call IBOF(ai,a4)
      
!     A2 expressed in the orthotropic frame
!
      call R2Ro(ai,dim,ap,angle)

!     Get viscosity

      CALL OPILGGE_ai_nl(ap, Angle, FabricGrid, C)

!    Compute strainRate and Spin :
!    -----------------------------

      LGrad = MATMUL( NodalVelo(:,1:n), dBasisdx(1:n,:) )

      StrainRate = 0.5 * ( LGrad + TRANSPOSE(LGrad) )

      Spin1 = 0.5 * ( LGrad - TRANSPOSE(LGrad) )

      IF ( CSymmetry ) THEN
        StrainRate(1,3) = 0.0
        StrainRate(2,3) = 0.0
        StrainRate(3,1) = 0.0
        StrainRate(3,2) = 0.0
        StrainRate(3,3) = 0.0

        Radius = SUM( Nodes % x(1:n) * Basis(1:n) )

        IF ( Radius > 10*AEPS ) THEN
         StrainRate(3,3) = SUM( Nodalvelo(1,1:n) * Basis(1:n) ) / Radius
        END IF

        epsi = StrainRate(1,1)+StrainRate(2,2)+StrainRate(3,3)
        DO i=1,3   
          StrainRate(i,i) = StrainRate(i,i) - epsi/3.0
        END DO

      ELSE
        epsi = StrainRate(1,1)+StrainRate(2,2)+StrainRate(3,3)
        DO i=1,dim 
          StrainRate(i,i) = StrainRate(i,i) - epsi/dim
        END DO

      END IF
      
!
!    Compute deviatoric stresses: 
!    ----------------------------
      D(1) = StrainRate(1,1)
      D(2) = StrainRate(2,2)
      D(3) = StrainRate(3,3)
      D(4) = 2. * StrainRate(1,2)
      D(5) = 2. * StrainRate(2,3)
      D(6) = 2. * StrainRate(3,1)
      
      INDi(1:6) = (/ 1, 2, 3, 1, 2, 3 /)
      INDj(1:6) = (/ 1, 2, 3, 2, 3, 1 /)
      DO k = 1, 2*dim
       DO j = 1, 2*dim
        Stress( INDi(k),INDj(k) ) = &
        Stress( INDi(k),INDj(k) ) + C(k,j) * D(j)
       END DO
       IF (k > 3)  Stress( INDj(k),INDi(k) ) = Stress( INDi(k),INDj(k) )
      END DO


!     SD=(1-r)D + r psi/2 S :
!     -----------------------
      SD=0._dp
      DO i=1,2*dim
        SD(i)= (1._dp - rho)*StrainRate(INDi(i),INDj(i)) + rho *&
                                   Theta *  Stress(INDi(i),INDj(i))
      END DO
      Do i=1,2*dim-3
        Spin(i)=Spin1(INDi(i+3),INDj(i+3))
      End do

      Deq=sqrt(2._dp*(SD(1)*SD(1)+SD(2)*SD(2)+SD(3)*SD(3)+2._dp* &
                            (SD(4)*SD(4)+SD(5)*SD(5)+SD(6)*SD(6)))/3._dp)
!
!     Velocity :
!     ----------
      Velo = 0.0d0
      DO i=1,dim
         Velo(i) = SUM( Basis(1:n) * (NodalVelo(i,1:n) - NodMeshVel(i,1:n)) )
      END DO
      Unorm = SQRT( SUM( Velo**2._dp ) )

!
!     Reaction coefficient:
!     ---------------------
      SELECT CASE(comp)
      CASE(1)
        !C0 = -2._dp*(SD(1)-SD(3))
        C0=-2._dp*SD(1)-3._dp*lambda*Deq

      CASE(2)
        !C0 = -2._dp*(SD(2)-SD(3))
        C0 = -2._dp*SD(2)-3._dp*lambda*Deq
       
      CASE(3)
        !C0 = -(SD(1)+SD(2)-2._dp*SD(3))
        C0 = -(SD(1)+SD(2))-3._dp*lambda*Deq
        
      CASE(4)
        C0 = -(SD(2)-SD(3))-3._dp*lambda*Deq

      CASE(5)
        C0 = -(SD(1)-SD(3))-3._dp*lambda*Deq

        
      END SELECT

      If (Fond) C0=0._dp

!     Loop over basis functions (of both unknowns and weights):
!     ---------------------------------------------------------
      DO p=1,NBasis
         DO q=1,NBasis
            A = 0.0d0
            M = Basis(p) * Basis(q)
!
!           Reaction terms:
!           ---------------
            A = A - C0 * Basis(q) * Basis(p)

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
         SELECT CASE(comp)
         
         CASE(1)
         
          LoadAtIp = 2._dp*( Spin(1)*E1 + SD(4)*(2._dp*a4(7)-E1) + &
          SD(1)*a4(1) + SD(2)*a4(3) - SD(3)*(-A1+a4(1)+a4(3)) ) + & 
                      lambda*Deq

          IF(dim == 3) THEN
            LoadAtIp =  LoadAtIp + 2._dp*( -Spin(3)*E3 + &
            SD(6)*(2._dp*a4(6)-E3) + 2._dp*SD(5)*a4(4) )
          END IF
          
         
         CASE(2)
          LoadAtIp = 2._dp*( -Spin(1)*E1 + SD(4)*(2._dp*a4(9)-E1) + &
          SD(2)*a4(2) + SD(1)*a4(3) - SD(3)*(-A2+a4(2)+a4(3)) )  +  &
                        lambda*Deq
          
          IF(dim == 3) THEN
            LoadAtIp =  LoadAtIp + 2._dp*( Spin(2)*E2 + &
            SD(5)*(2._dp*a4(8)-E2) + 2._dp*SD(6)*a4(5) )
          END IF


         CASE(3)
          LoadAtIp = Spin(1)*(A2-A1)  +  SD(4)*(4._dp*a4(3)-A1-A2) + &
          2._dp* ( SD(1)*a4(7) + SD(2)*a4(9) - SD(3)*(-E1+a4(7)+a4(9)) )
          
          IF(dim == 3) THEN
           LoadAtIp =  LoadAtIp - Spin(3)*E2 + Spin(2)*E3  &
           + SD(6)*(4._dp*a4(4)-E2) + SD(5)*(4._dp*a4(5)-E3)  
          END IF
           
         CASE(4)
          LoadAtIp = Spin(2)*(A3-A2) +  Spin(3)*E1 - Spin(1)*E3 +&
          SD(4)*(4._dp*a4(5)-E3) + SD(5)*(3._dp*A2-A3-4._dp*(a4(2)+a4(3))) &
          + SD(6)*(3._dp*E1-4._dp*(a4(7)+a4(9))) + &
          2._dp*( SD(1)*a4(4) + SD(2)*a4(8) - SD(3)*(a4(4)+a4(8)) )

         CASE(5)
          LoadAtIp = Spin(3)*(A1-A3) +  Spin(1)*E2 - Spin(2)*E1 +&
           SD(4)*(4._dp*a4(4)-E2) + &
          SD(6)*(3._dp*A1-A3-4.*(a4(1)+a4(3))) + SD(5)*(3._dp*E1-4._dp*(a4(7)+a4(9))) + &
          2._dp*( SD(1)*a4(6) + SD(2)*a4(5) - SD(3)*(a4(6)+a4(5)) )



         END SELECT

        If (Fond) LoadAtIp=0._dp
        LoadAtIp= LoadAtIp * Basis(p)
        FORCE(p) = FORCE(p) + s*LoadAtIp
      END DO

              
      END DO 

 1000 FORMAT((a),x,i2,x,6(e13.5,x)) 
 1001 FORMAT(6(e13.5,x))
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

      Save EdgeNodes, LeftParentNodes, RightParentNodes
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

     REAL(KIND=dp) :: Normal(3), g, L, Udotn, UdotnA, cu(3), detJ,U,V,W,S
     LOGICAL :: Stat
     TYPE(GaussIntegrationPoints_t) :: IntegStuff

     TYPE(Nodes_t) :: Nodes, ParentNodes
     SAVE Nodes, ParentNodes
!------------------------------------------------------------------------------
     dim = CoordinateSystemDimension()
     FORCE = 0.0d0
     STIFF = 0.0d0

     CALL GetElementNodes( Nodes, Element )
     CALL GetElementNodes( ParentNodes, ParentElement )

     ! Numerical integration:
     !-----------------------
     IntegStuff = GaussPoints( Element )
!
! Compute the average velocity.dot.Normal        
!        
     UdotnA = 0.0   
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
       cu = 0.0d0
       DO i=1,dim
          cu(i) = SUM( (Velo(i,1:n)-MeshVelo(i,1:n)) * Basis(1:n) )
       END DO
       UdotnA = UdotnA + s*SUM( Normal * cu )

     END DO

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
        IF (InFlowBC .And. (UdotnA < 0.) ) THEN
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
 END SUBROUTINE FabricSolver
!------------------------------------------------------------------------------

