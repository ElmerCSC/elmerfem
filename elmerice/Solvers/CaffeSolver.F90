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
! *  Authors: Hakime Seddik
! *
! *  Address: Institute of Low Temperature Science
! *           Hokkaido University
! *           Kita 19, Nishi8, Kita-ku, Sapporo, 060-0819
! *           Japan    
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 06 October 2009
! * 
! *****************************************************************************
!> Module containing a solver for fabric equations using a diffusion       
!> term for numerical stabilization. Caffe based fabric evolution equation.
  RECURSIVE SUBROUTINE CaffeSolver( Model,Solver,dt,TransientSimulation )
!---------------------------------------------------------------------------------------------------------------
  USE DefUtils
  USE MaterialModels
  USE Integration
  USE Differentials

  IMPLICIT NONE

!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve fabric equations for one timestep
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

  INTEGER :: dim,n1,n2,i,j,k,l,n,t,iter,STDOFs, istat, COMP

  TYPE(ValueList_t),POINTER :: Material, BC, Equation, SolverParams
  TYPE(Nodes_t) :: ElementNodes, ParentNodes
  TYPE(Element_t),POINTER :: CurrentElement, Element, ParentElement

  REAL(KIND=dp) :: RelativeChange,UNorm,PrevUNorm,  &
                  NewtonTol,NonlinearTol,s, C0


  INTEGER :: NewtonIter,NonlinearIter

  TYPE(Variable_t), POINTER :: FabricSol, FabricVariable, FlowVariable, &
                                  MeshVeloVariable

  TYPE(Variable_t), POINTER :: a11, a22, a12, a23, a13 

  REAL(KIND=dp), POINTER :: Fabric(:), &
           FabricValues(:), FlowValues(:), &
           MeshVeloValues(:), Solution(:), Ref(:), Hwrk(:,:,:)

  INTEGER, POINTER :: FabricPerm(:),NodeIndexes(:), ParentNIndexes(:) , &
                         FlowPerm(:), MeshVeloPerm(:)

  LOGICAL :: GotIt, Found, NewtonLinearization = .FALSE., FlowSolutionFound, Bubbles = .FALSE., &
          UseBubbles, Stabilize = .TRUE., FluxBC, FreeBC, SpinShear
  INTEGER :: NSDOFs, body_id,bf_id,eq_id, material_id, other_body_id

  INTEGER :: old_body = -1
                        
  LOGICAL :: AllocationsDone = .FALSE., FreeSurface

  TYPE(Variable_t), POINTER :: TimeVar

  REAL(KIND=dp), ALLOCATABLE:: MASS(:,:), STIFF(:,:),  &
       LOAD(:),Force(:), &
       K1(:), K2(:), E1(:), E2(:), E3(:), &
       IceVeloU(:), IceVeloV(:), IceVeloW(:), Velocity(:,:), MeshVelocity(:,:), &
       iota(:), C1(:), CT(:), Zero(:), TimeForce(:), a2Diffuse(:,:,:)
  
  CHARACTER(LEN=MAX_NAME_LEN) :: ConvectionFlag, FlowSolName, SolverName

  SAVE MASS, STIFF, LOAD, Force, ElementNodes,& 
         AllocationsDone, K1, K2, &
         E1, E2, E3, iota, IceVeloU, IceVeloV, IceVeloW, C1, CT, a2Diffuse, Velocity, &
         MeshVelocity, Zero, Hwrk, TimeForce, old_body, dim, comp, SolverName, ConvectionFlag
!--------------------------------------------------------------------------------------------------------------------------------------------

  REAL(KIND=dp) :: SaveTime = -1
  REAL(KIND=dp), POINTER :: PrevFabric(:),CurrFabric(:),TempFabVal(:)

  SAVE  PrevFabric, CurrFabric,TempFabVal
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at, at0
#else
  REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif
                                                
!--------------------------------------------------------------------------------
 CALL Info( 'CaffeSolveFabric', 'version 4.0', Level=4 )
!--------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
 IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

  SolverName = 'CAFFE (Fabric solver)'

  Solution => Solver % Variable % Values
  STDOFs   =  Solver % Variable % DOFs

  FabricSol    => VariableGet( Solver % Mesh % Variables, 'Fabric' )
  IF ( ASSOCIATED( FabricSol) ) THEN
        FabricPerm   => FabricSol % Perm
        FabricValues => FabricSol % Values
  ELSE
        CALL Fatal( 'CAFFE (Fabric solver)', 'No Fabric associated!.' ) 
  END IF

  MeshVeloVariable => VariableGet( Solver % Mesh % Variables, &
            'Mesh Velocity' )
  IF ( ASSOCIATED( MeshVeloVariable ) ) THEN
       MeshVeloPerm    => MeshVeloVariable % Perm
       MeshVeloValues  => MeshVeloVariable % Values
  END IF
                 
  !--------------------------------------------------------------------------------------------------
  ! We are using exported variable so that the routine DefaultDirichletBCs can use those variables 
  ! name for setting dirichlet boundary conditions for each a22 components
  !--------------------------------------------------------------------------------------------------
  a11    => VariableGet( Solver % Mesh % Variables, 'a11' )
  IF ( .NOT. ASSOCIATED( a11) ) THEN
     CALL Fatal( 'CAFFE (Fabric solver)', 'No variable a11 associated!.' ) 
  END IF

  a22    => VariableGet( Solver % Mesh % Variables, 'a22' )
  IF ( .NOT. ASSOCIATED( a22) ) THEN
     CALL Fatal( 'CAFFE (Fabric solver)', 'No variable a22 associated!.' ) 
  END IF

  a12    => VariableGet( Solver % Mesh % Variables, 'a12' )
  IF ( .NOT. ASSOCIATED( a12) ) THEN
     CALL Fatal( 'CAFFE (Fabric solver)', 'No variable a12 associated!.' ) 
  END IF

  a23    => VariableGet( Solver % Mesh % Variables, 'a23' )
  IF ( .NOT. ASSOCIATED( a23) ) THEN
     CALL Fatal( 'CAFFE (Fabric solver)', 'No variable a23 associated!.' ) 
  END IF

  a13    => VariableGet( Solver % Mesh % Variables, 'a13' )
  IF ( .NOT. ASSOCIATED( a13) ) THEN
     CALL Fatal( 'CAFFE (Fabric solver)', 'No variable a13 associated!.' ) 
  END IF
                            
  ! UNorm = Solver % Variable % Norm
  Unorm = SQRT( SUM( FabricValues**2 ) / SIZE(FabricValues) )

!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
        N = Model % MaxElementNodes

       dim = CoordinateSystemDimension()
       
     IF ( AllocationsDone ) THEN
         DEALLOCATE(               &
               K1,               &
               K2,               &
               E1,               &
               E2,               &
               E3,               &      
               Force,            &
               Velocity,         &
               IceVeloU,         &
               IceVeloV,         &
               IceVeloW,         &
               MeshVelocity,     &
               MASS,             & 
               STIFF,            &
               LOAD,             &
               CurrFabric,       &
               TempFabVal,       &
               iota,             &
               C1,               &
               CT,               &
               a2Diffuse,        &
               Zero,             &
               TimeForce,        &
               ParentNodes % x,  &
               ParentNodes % y,  &
               ParentNodes % z  )
     END IF

     ALLOCATE(                   &   
           K1( N ),            &
           K2( N ),            &
           E1( N ),            &
           E2( N ),            &
           E3( N ),            &
           Force( 2*N ),       &
           Velocity(4, N),     &
           IceVeloU(N),        &
           IceVeloV(N),        &
           IceVeloW(N),        &
           MeshVelocity(3,N),  &
           MASS( 2*N,2*N ),    &
           STIFF( 2*N,2*N ),   &
           LOAD( N ),          &
           iota(N),            &
           C1(N),              &
           CT(N),              &
           Zero(N),            &
           TimeForce(2*N),     &
           a2Diffuse(3,3,N),   &
           CurrFabric( 5*SIZE(Solver % Variable % Values)), &
           TempFabVal( 5*SIZE(Solver % Variable % Values)), &
           ParentNodes % x (N), &
           ParentNodes % y (N), &
           ParentNodes % z (N), & 
           STAT=istat )

     IF ( istat /= 0 ) THEN
        CALL Fatal( 'CAFFE (Fabric solver)', 'Memory allocation error.' )
     END IF
     CALL Info('CAFFE (Fabric solver)','Memory allocations done', Level=3)  

     CurrFabric = 0.
     TempFabVal = 0.
     IF ( TransientSimulation ) THEN
        IF (AllocationsDone ) DEALLOCATE (PrevFabric)
        ALLOCATE( PrevFabric( 5*SIZE(Solver % Variable % Values)) )
        PrevFabric = 0.
     END IF

     DO i=1,Solver % NumberOFActiveElements
        CurrentElement => GetActiveElement(i)   
        n = GetElementNOFNodes()
        NodeIndexes => CurrentElement % NodeIndexes
        DO COMP=1,5
           IF ( TransientSimulation ) THEN
              PrevFabric(5*(Solver % Variable % Perm(NodeIndexes(1:n))-1)+COMP) = &
                 FabricValues(5*(FabricPerm(NodeIndexes(1:n))-1)+COMP)
           END IF
           CurrFabric(5*(Solver % Variable % Perm(NodeIndexes(1:n))-1)+COMP) = &
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
  NonlinearTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Convergence Tolerance' )

  NewtonTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Newton After Tolerance' )

  NewtonIter = ListGetInteger( Solver % Values, &
        'Nonlinear System Newton After Iterations' )

  NonlinearIter = ListGetInteger( Solver % Values, &
         'Nonlinear System Max Iterations',GotIt )
  IF ( .NOT.GotIt ) NonlinearIter = 1

  SolverParams => GetSolverParams()
  Stabilize = GetLogical( SolverParams,'Stabilize',Found )
  IF (.NOT. Found) Stabilize = .FALSE.
  UseBubbles = GetLogical( SolverParams,'Bubbles',Found )
  IF ( .NOT.Found .AND. (.NOT.Stabilize)) UseBubbles = .TRUE.

  SpinShear = GetLogical( SolverParams, 'Shearing in Fabric', Found )
  IF (.NOT. Found) SpinShear = .TRUE.

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  DO iter=1,NonlinearIter
!------------------------------------------------------------------------------
   at  = CPUTime()
   at0 = RealTime()

   CALL Info( 'CaffeSolveFabric', ' ', Level=4 )
   CALL Info( 'CaffeSolveFabric', ' ', Level=4 )
   CALL Info( 'CaffeSolveFabric', &
                '-------------------------------------',Level=4 )
   WRITE( Message, * ) 'Caffe Fabric solver  iteration', iter
   CALL Info( 'CaffeSolveFabric', Message,Level=4 )
   CALL Info( 'CaffeSolveFabric', &
                '-------------------------------------',Level=4 )

   PrevUNorm = UNorm
   
   !---------------------------------------------------   
   ! Loop over the number of fabric equations to solve
   !---------------------------------------------------
    DO COMP=1,2*dim-1
       
     Solver % Variable % Values = CurrFabric( COMP::5 )
     IF ( TransientSimulation ) THEN
          Solver % Variable % PrevValues(:,1) = PrevFabric( COMP::5 )
     END IF

     CALL DefaultInitialize()

     CALL Info( 'CaffeSolveFabric', ' ', Level=4 )
     CALL Info( 'CaffeSolveFabric', 'Starting assembly...',Level=4 )
     
     !------------------------------------------------------------------------------
     DO t=1,Solver % NumberOFActiveElements
     !------------------------------------------------------------------------------
      IF ( RealTime() - at0 > 1.0 ) THEN
          WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
           (Solver % NumberOfActiveElements-t) / &
               (1.0*Solver % NumberOfActiveElements)), ' % done'
          CALL Info( 'CaffeSolveFabric', Message, Level=5 )
          at0 = RealTime()
      END IF 

      CurrentElement => GetActiveElement(t)
      CALL GetElementNodes( ElementNodes, CurrentElement )
      n = GetElementNOFNodes(CurrentElement)
      NodeIndexes => CurrentElement % NodeIndexes

      body_id = CurrentElement % BodyId

      !------------------------------------------------------------------------------
      ! Read in material constants from Material section
      !------------------------------------------------------------------------------
      IF (body_id /= old_body) Then 
         old_body = body_id

         Material => GetMaterial()
         IF (.NOT. ASSOCIATED(Material)) THEN
            WRITE(Message,'(A,I5,A)') 'No Material for bulk element no. ',t,' found.'
            CALL FATAL('CAFFE (Fabric solver)',Message)
         END IF

         Equation => GetEquation()
         IF (.NOT.ASSOCIATED(Equation)) THEN
            WRITE (Message,'(A,I3)') 'No Equation for builk  element no. ', t, 'found.'
            CALL FATAL(SolverName,Message)
         END IF

         ConvectionFlag = GetString( Equation, 'Convection', Found )

         iota(1:n) = ListGetConstReal( Material,'Grain Rotation parameter', GotIt )
         IF (.NOT. GotIt) THEN
            WRITE(Message,'(A)') 'No grain rotation parameter (iota) found in Material '
            CALL Fatal('CAFFE (Fabric solver)', Message)
         ELSE 
            WRITE(Message,'(A,F10.4)') 'Iota parameter = ', iota(1)
            CALL INFO('CAFFE (Fabric solver)', Message, Level = 9)
         END IF

         CALL ListGetRealArray( Material, 'Fabric Diffusion',Hwrk,n, CurrentElement % NodeIndexes, GotIt)
         IF (.NOT. GotIt) THEN
             WRITE(Message,'(A)') 'No diffusion coefficient in Material '
            CALL Fatal('CAFFE (Fabric solver)', Message)
         END IF
         a2Diffuse = 0.0d0
         IF ( SIZE(Hwrk,1) == 1 ) THEN
            DO i=1,3
               a2Diffuse( i,i,1:N ) = Hwrk( 1,1,1:N)
            END DO
         ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN
            DO i=1,MIN(3,SIZE(Hwrk,1))
               a2Diffuse(i,i,1:N) = Hwrk(i,1,1:N)
            END DO
         ELSE
            DO i=1,MIN(3,SIZE(Hwrk,1))
              DO j=1,MIN(3,SIZE(Hwrk,2))
                 a2Diffuse(i,j,1:N) = Hwrk(i,j,1:N)
              END DO
            END DO
         END IF

         IF (.NOT. SpinShear) THEN 
            CALL INFO(SolverName, 'Shearing and rigid body rotation not included in fabric', Level=1)
         ELSE
           CALL INFO(SolverName, 'Shearing and rigid body rotation included in fabric', Level=1)
         END IF

         NULLIFY(Material)   

      END IF

      k = ListGetInteger( Model % Bodies(CurrentElement % BodyId) % Values, 'Equation', &
                   minv=1, maxv=Model % NumberOFEquations )

      SELECT CASE( ListGetString( Model % Equations(k) % Values, 'Convection', Found ) )

      CASE( 'computed' )
       FlowSolName =  GetString( Model % Equations(k) % Values,'Flow Solution Name', Found)
       IF(.NOT.Found) THEN        
         CALL WARN(SolverName,'Keyword >Flow Solution Name< not found in section >Equation<')
         CALL WARN(SolverName,'Taking default value >Flow Solution<')
         WRITE(FlowSolName,'(A)') 'Flow Solution'
       END IF

       FlowVariable => VariableGet( Solver % Mesh % Variables, FlowSolName )
       IF ( ASSOCIATED( FlowVariable ) ) THEN
            FlowPerm     => FlowVariable % Perm
            FlowValues => FlowVariable % Values
            NSDOFs       =  FlowVariable % DOFs
            FlowSolutionFound = .TRUE.
       ELSE
            CALL INFO(SolverName,'No Flow Solution associated',Level=1)
            FlowSolutionFound = .FALSE.
       END IF

      CASE( "none")
       FlowSolutionFound = .FALSE.

      END SELECT
      
      Material => GetMaterial()
      IF (.NOT. ASSOCIATED(Material)) THEN
         WRITE(Message,'(A,I5,A)') 'No Material for bulk element no. ',t,' found.'
         CALL FATAL('CAFFE (Fabric solver)',Message)
      ELSE 
         material_id = GetMaterialId( CurrentElement, Found)
         IF(.NOT.Found) THEN
           WRITE (Message,'(A,I3)') 'No Material ID found for boundary element no. ', t
           CALL FATAL(SolverName,Message)
         END IF
      END IF

      !------------------------------------------------------------------------------
      ! Get element local stiffness & mass matrices
      !------------------------------------------------------------------------------
      !-------------------------------
      ! Second order tensor components
      !-------------------------------
      K1(1:n) = CurrFabric( 5*(Solver % Variable % Perm( NodeIndexes(1:n) )-1)+1 )
      K2(1:n) = CurrFabric( 5*(Solver % Variable % Perm( NodeIndexes(1:n) )-1)+2 )
      E1(1:n) = CurrFabric( 5*(Solver % Variable % Perm( NodeIndexes(1:n) )-1)+3 )
      E2(1:n) = CurrFabric( 5*(Solver % Variable % Perm( NodeIndexes(1:n) )-1)+4 )
      E3(1:n) = CurrFabric( 5*(Solver % Variable % Perm( NodeIndexes(1:n) )-1)+5 )

      !-----------------------------------------------------
      ! We need that constant to be 1 so that the assembly 
      ! routine takes properly ino account the convection
      !-----------------------------------------------------
      C1 = 1.0D00

      !-------------------------------
      ! Mesh Velocoty
      !-------------------------------
      MeshVelocity=0._dp
      IF (ASSOCIATED(MeshVeloVariable)) Then
          k = MeshVeloVariable % DOFs
          DO i=1,k
             MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(NodeIndexes(1:n))-1)+i)
          END DO
      END IF
      
      !-------------------------------
      ! Velocity
      !-------------------------------
      IceVeloU = 0.0d00
      IceVeloV = 0.0d00
      IceVeloW = 0.0d00
      ! constant (i.e., in section Material given) velocity
      !---------------------------------------------------
      IF ( ConvectionFlag == 'constant' ) THEN
         IceVeloU(1:N)= GetReal( Material, 'Convection Velocity 1', Found )
         IceVeloV(1:N) = GetReal( Material, 'Convection Velocity 2', Found )
         IceVeloW(1:N) = GetReal( Material, 'Convection Velocity 3', Found )                 
      ! computed velocity
      !------------------
      ELSE IF (( ConvectionFlag == 'computed' ) .AND. FlowSolutionFound) THEN
         DO i=1,n
            k = FlowPerm(CurrentElement % NodeIndexes(i))
            IF ( k > 0 ) THEN
              SELECT CASE( NSDOFs )
               CASE(3)
                IceVeloU(i) = FlowValues( NSDOFs*k-2 )
                IceVeloV(i) = FlowValues( NSDOFs*k-1 )
                IceVeloW(i) = 0.0D0
               CASE(4)
                IceVeloU(i) = FlowValues( NSDOFs*k-3 )
                IceVeloV(i) = FlowValues( NSDOFs*k-2 )
                IceVeloW(i) = FlowValues( NSDOFs*k-1 )
              END SELECT
            END IF
         END DO
         WRITE(Message,'(a,i5, a, i5)') 'Convection in element ', t, &
                      ' material ',  material_id
      ELSE  ! Conduction and ALE contribution only
         IF (ANY( MeshVelocity /= 0.0d0 )) THEN
            WRITE(Message,'(a,i5, a, i5)') 'Only mesh deformation in element ', t,&
                         ' material ',  material_id
         ELSE ! neither convection nor ALE mesh deformation contribution -> all C1(1:N)=0
            C1 = 0.0D0 
            WRITE(Message,'(a,i5, a, i5)') 'No convection and mesh deformation in element ', t,&
                         ' material ',  material_id
         END IF                 
      END IF
      CALL INFO(SolverName,Message,Level=10)

      !----------------------------------------------------------------------
      ! For convenience, set a Velocity table that will be 
      ! used to compute the velocity gradient and the free boudary condition
      !----------------------------------------------------------------------
      k = FlowVariable % DOFs
      Velocity = 0.0d0
      DO i=1,k-1
         Velocity(i,1:n) = FlowValues(k*(FlowPerm(NodeIndexes(1:n))-1)+i)
      END DO

      !-----------------------------
      ! Dummy table for zero values.
      !-----------------------------
      Zero = 0.00D00
      
      !-----------------------------------------
      ! We will compute the reaction term later
      !-----------------------------------------
      C0 = 0.0D00
      
      !---------------------------------------------------------
      ! We don't have any coefficient with the time derivative
      !---------------------------------------------------------
      CT = 1.0D00

      !------------------------------------------------------------------------------
      ! Do we really need residual free Bubbles
      !------------------------------------------------------------------------------
      Bubbles = UseBubbles  .AND. &
                ( ConvectionFlag == 'computed' .OR. ConvectionFlag == 'constant' )     

      !------------------------------------------------------------------------------
      ! Get element local matrices, and RHS vectors
      !------------------------------------------------------------------------------
      MASS = 0.0d00
      STIFF = 0.0d00
      FORCE = 0.0D00
 
      CALL DiffuseConvectiveComposeCaffe( &
          COMP, MASS, STIFF, FORCE, LOAD, CT, C0, C1(1:N), a2Diffuse, &
          .FALSE., Zero, Zero, K1, K2, E1, E2, E3, Velocity, IceVeloU, IceVeloV, IceVeloW, &
          MeshVelocity(1,1:N), MeshVelocity(2,1:N),MeshVelocity(3,1:N),Zero, Zero, & 
          Zero, Zero, Zero, .FALSE., Stabilize, Bubbles, CurrentElement, n, ElementNodes, &
          iota, SpinShear)

      !------------------------------------------------------------------------------
      ! If time dependent simulation add mass matrix to stiff matrix
      !------------------------------------------------------------------------------
      TimeForce  = FORCE
      IF ( TransientSimulation ) THEN
         IF ( Bubbles ) FORCE = 0.0d0
            CALL Default1stOrderTime( MASS,STIFF,FORCE )
      END IF
      !------------------------------------------------------------------------------
      !  Update global matrices from local matrices
      !------------------------------------------------------------------------------
      IF (  Bubbles ) THEN
         CALL Condensate( N, STIFF, FORCE, TimeForce )
         IF (TransientSimulation) CALL DefaultUpdateForce( TimeForce )
      END IF

      CALL DefaultUpdateEquations( STIFF, FORCE )

     END DO     !  Bulk elements
     CALL Info( 'CaffeSolveFabric', 'Assembly done', Level=4 )

     !------------------------------------------------------------------------------
     !     Loop over the boundary elements
     !------------------------------------------------------------------------------
     DO t = 1, Solver % Mesh % NumberOfBoundaryElements
        
      Element => GetBoundaryElement(t)
      IF( .NOT. ActiveBoundaryElement() )  CYCLE
      IF( GetElementFamily(Element) == 1 ) CYCLE

      ParentElement => Element % BoundaryInfo % Left
      IF ( .NOT. ASSOCIATED( ParentElement ) ) &
         ParentElement => Element % BoundaryInfo % Right
      n  = GetElementNOFNodes( Element )
      n1 = GetElementNOFnodes( ParentElement )
      NodeIndexes => Element % NodeIndexes

      BC => GetBC()
      CALL GetElementNodes( ElementNodes, Element )
!      CALL GetElementNodes( ParentNodes, ParentElement )

      IF ( ASSOCIATED( BC ) ) THEN
         STIFF=0.0D00
         FORCE=0.0D00
         MASS=0.0D00
         LOAD=0.0D00     

         FluxBC = .FALSE.
         FreeBC = .FALSE.

         FluxBC = GetLogical(BC, 'Fabric Flux BC', Found)
         FreeBC = GetLogical(BC, 'Fabric Free BC', Found)
       
         IF (FluxBC .AND. .NOT. FreeBC) THEN
            !---------------
            !BC: -k@T/@n = q
            !---------------
            LOAD(1:N) = LOAD(1:N) + GetReal( BC, ComponentName('Fabric Flux', Comp) , GotIt )
            
            CALL DiffuseConvectiveBoundaryCaffe( STIFF,FORCE, &
                            LOAD, Zero, Element,n,ElementNodes )

         ELSE If (FreeBC .AND. .NOT. FluxBC) THEN
                   
            MeshVelocity=0._dp
            IF (ASSOCIATED(MeshVeloVariable)) Then
             k = MeshVeloVariable % DOFs
            DO i=1,k
              MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(NodeIndexes(1:n))-1)+i)
             END DO
            END IF   

            k = FlowVariable % DOFs
            Velocity = 0.0d0
            DO i=1,k-1
              Velocity(i,1:n) = ( FlowValues(k*(FlowPerm(NodeIndexes(1:n))-1)+i) )
            END DO

            CALL FreeBoundaryCaffe( STIFF,FORCE, &
                 LOAD, Element, n , ElementNodes, ParentElement, n1, ParentNodes, &
                 Velocity, MeshVelocity )

!            CALL FreeBoundaryCaffe2( Comp, MASS, STIFF, FORCE, LOAD, &
!              K1, K2, E1, E2, E3, Velocity, MeshVelocity, Element, n, ElementNodes, ParentElement, &
!              ParentNodes, n1, iota, SpinShear)
        
         ELSE IF (FluxBC .AND. FreeBC) THEN
            CALL Fatal(SolverName, 'Intend to use flux and free boundary condition at the same time')
         END IF
      END IF

      ! Update global matrices from local matrices
      !------------------------------------------------------------------------------
      IF ( TransientSimulation ) THEN
           MASS = 0.d0
           CALL Default1stOrderTime( MASS, STIFF, FORCE )
      END IF
      CALL DefaultUpdateEquations( STIFF, FORCE )
 
     END DO   ! Neumann & Newton BCs

     CALL DefaultFinishAssembly()
     !--------------------------------------------------------------------------------------------
     ! We need to pass independent variables for the Fabric so that this routine works as expected
     ! Here a11, a22, a12, a23, a13 are real exported variables not DOFs of the variable Fabric
     !--------------------------------------------------------------------------------------------
     SELECT CASE( Comp )
      CASE(1)
       CALL DefaultDirichletBCs(Solver, a11)
      CASE(2)
       CALL DefaultDirichletBCs(Solver, a22)
      CASE(3)
       CALL DefaultDirichletBCs(Solver, a12)
      CASE(4)
       CALL DefaultDirichletBCs(Solver, a23)
      CASE(5)
       CALL DefaultDirichletBCs(Solver, a13)
     END SELECT

     CALL Info( 'CaffeSolveFabric', 'Set boundaries done', Level=4 )

     !------------------------------------------------------------------------------
     ! Solve the system and check for convergence
     !------------------------------------------------------------------------------
     Unorm = DefaultSolve()
     WRITE(Message,*) 'solve done', minval( solver % variable % values), &
          maxval( Solver % variable % values)
     CALL Info( 'CaffeSolveFabric', Message, Level=4 )

     !-----------------------------------------------------------------------
     ! Copy solution of a2 components being solved to variable Fabric
     !----------------------------------------------------------------------
     DO t=1,Solver % NumberOfActiveElements
        Element => GetActiveElement(t) 
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
        DO i=1,n
          k = NodeIndexes(i)
          FabricValues( 5*(FabricPerm(k)-1) + COMP ) = & 
            Solver % Variable % Values( Solver % Variable % Perm(k) )
        END DO
     END DO

     SELECT CASE( Comp ) 
      CASE(1)
       FabricValues( COMP:SIZE(FabricValues):5 ) = &
          MIN(MAX( FabricValues( COMP:SIZE(FabricValues):5 ) , 0._dp),1._dp)
          
      CASE(2)
       FabricValues( COMP:SIZE(FabricValues):5 ) = &
       MIN(MAX( FabricValues( COMP:SIZE(FabricValues):5 ) , 0._dp),1._dp)

      CASE(3:5)
       DO i=COMP,SIZE(FabricValues),5
         IF(FabricValues(i).GT.0._dp) THEN
               FabricValues(i) = MIN( FabricValues(i) , 0.5_dp)
         ELSE
               FabricValues(i) = MAX( FabricValues(i) , -0.5_dp)
         END IF
       END DO
     END SELECT
!----------------------------------------------------------
! End of itreations for all a2 components
!----------------------------------------------------------
    END DO 
!----------------------------------------------------------

    DO i=1,Solver % NumberOFActiveElements
          CurrentElement => GetActiveElement(i)   
          n = GetElementNOFNodes()
          NodeIndexes => CurrentElement % NodeIndexes
          DO COMP=1,5
            CurrFabric(5*(Solver % Variable % Perm(NodeIndexes(1:n))-1)+COMP) = &
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
    CALL Info( 'CaffeSolveFabric', Message, Level=4 )
    WRITE( Message, * ) 'Relative Change : ',RelativeChange
    CALL Info( 'CaffeSolveFabric', Message, Level=4 )
  
!------------------------------------------------------------------------------
    IF ( RelativeChange < NewtonTol .OR. &
            iter > NewtonIter ) NewtonLinearization = .TRUE.

    IF ( RelativeChange < NonLinearTol ) EXIT

!------------------------------------------------------------
! End of non-linear iterations
!------------------------------------------------------------
  END DO
!------------------------------------------------------------

CONTAINS

!*******************************************************************************
! *
! * Diffuse-convective local matrix computing (cartesian coordinates)
! *
! *******************************************************************************
          
!------------------------------------------------------------------------------
   SUBROUTINE DiffuseConvectiveComposeCaffe( COMP, MassMatrix,StiffMatrix,ForceVector,  &
      LoadVector,NodalCT,NodalC0,NodalC1,NodalC2,PhaseChange,NodalTemperature, &
      Enthalpy, NodalK1, NodalK2, NodalEuler1, NodalEuler2, NodalEuler3, &
      NodalVelocity, Ux,Uy,Uz,MUx,MUy,MUz,Nodalmu,Nodalrho,NodalPressure, NodaldPressureDt, &
      NodalPressureCoeff, Compressible, Stabilize, UseBubbles, Element,n,Nodes, iota, SpinShear)
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Return element local matrices and RSH vector for diffusion-convection
!  equation: 
!
!  ARGUMENTS:
!
!  REAL(KIND=dp) :: MassMatrix(:,:)
!     OUTPUT: time derivative coefficient matrix
!
!  REAL(KIND=dp) :: StiffMatrix(:,:)
!     OUTPUT: rest of the equation coefficients
!
!  REAL(KIND=dp) :: ForceVector(:)
!     OUTPUT: RHS vector
!
!  REAL(KIND=dp) :: LoadVector(:)
!     INPUT:
!
!  REAL(KIND=dp) :: NodalCT,NodalC0,NodalC1
!     INPUT: Coefficient of the time derivative term, 0 degree term, and
!            the convection term respectively
!
!  REAL(KIND=dp) :: NodalC2(:,:,:)
!     INPUT: Nodal values of the diffusion term coefficient tensor
!
!  LOGICAL :: PhaseChange
!     INPUT: Do we model phase change here...
!
!  REAL(KIND=dp) :: NodalTemperature
!     INPUT: NodalTemperature from previous iteration
!
!  REAL(KIND=dp) :: Enthalpy
!     INPUT: Enthalpy from previous iteration, needed if we model
!            phase change
!
!  REAL(KIND=dp) :: UX(:),UY(:),UZ(:)
!     INPUT: Nodal values of velocity components from previous iteration
!           used only if coefficient of the convection term (C1) is nonzero
!
!  REAL(KIND=dp) :: Nodalmu(:)
!     INPUT: Nodal values of the viscosity
!
!  LOGICAL :: Stabilize
!     INPUT: Should stabilzation be used ? Used only if coefficient of the
!            convection term (C1) is nonzero
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!  INTEGER :: n
!       INPUT: Number of element nodes
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!
!******************************************************************************

     REAL(KIND=dp), DIMENSION(:)   :: ForceVector,UX,UY,UZ,MUX,MUY,MUZ,LoadVector
     REAL(KIND=dp), DIMENSION(:,:) :: MassMatrix,StiffMatrix
     REAL(KIND=dp) :: NodalTemperature(:),Enthalpy(:),Nodalmu(:), &
       NodaldPressureDt(:),NodalPressure(:),NodalPressureCoeff(:),Nodalrho(:)
     REAL(KIND=dp) :: NodalC0,NodalC1(:),NodalCT(:),NodalC2(:,:,:), NodalVelocity(:,:)
     REAL(KIND=dp), DIMENSION(:) :: NodalK1, NodalK2, NodalEuler1, &
               NodalEuler2, NodalEuler3, iota

     LOGICAL :: UseBubbles,PhaseChange,Compressible,Stabilize, SpinShear

     INTEGER :: n, COMP

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), POINTER :: Element

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
!
     CHARACTER(LEN=MAX_NAME_LEN) :: StabilizeFlag
     REAL(KIND=dp) :: dBasisdx(2*n,3),detJ
     REAL(KIND=dp) :: Basis(2*n)
     REAL(KIND=dp) :: ddBasisddx(n,3,3),dNodalBasisdx(n,n,3)

     REAL(KIND=dp) :: Velo(3),Grad(3,3),Force, Spin(3), Spin1(3,3), divu, iotaAtIp, LoadAtIp

     REAL(KIND=dp) :: A1, A2, A3, E1, E2, E3

     REAL(KIND=dp) :: LGrad(3,3),StrainRate(3,3), epsi, Radius

     REAL(KIND=dp) :: A,M
     REAL(KIND=dp) :: Load

     REAL(KIND=dp) :: VNorm,hK,mK
     REAL(KIND=dp) :: Lambda=1.0,Pe,Pe1,Pe2,Tau,x,y,z
     REAL(KIND=dp) :: ai(6),a4(9)

     REAL(KIND=dp) :: Tau_M, Tau_C, Gmat(3,3),Gvec(3), dt=0._dp, LC1, NodalVelo(4,n), &
       RM(n),LC(3,n), gradP(n), VRM(3), GradNodal(n,3,3), PVelo(3), NodalPVelo(3,n), &
       Work(3,n), Grav(3), expc, reft, temperature

     REAL(KIND=dp), POINTER :: gWrk(:,:)

     INTEGER :: i,j,k,c,p,q,t,dim,N_Integ,NBasis,Order, INDi(6),INDj(6)

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp) :: s,u,v,w,dEnth,dTemp,mu,DivVelo,Pressure,rho,&
                      Pcoeff, minl, maxl

     REAL(KIND=dp) :: C0,C00,C1,CT,C2(3,3),dC2dx(3,3,3),SU(n),SW(n)

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: Vms, Found, Transient, stat,Convection,ConvectAndStabilize,Bubbles, CSymmetry

!------------------------------------------------------------------------------

     StabilizeFlag = GetString( GetSolverParams(),'Stabilization Method',Found )
     Vms = StabilizeFlag == 'vms'

     Transient = GetString(GetSimulation(),'Simulation type',Found)=='transient'

     dim = CoordinateSystemDimension()
     c = dim + 1

     ForceVector = 0.0D0
     StiffMatrix = 0.0D0
     MassMatrix  = 0.0D0
     Load = 0.0D0

     Convection =  ANY( NodalC1 /= 0.0d0 )
     NBasis = n
     Bubbles = .FALSE.
     IF ( Convection .AND. .NOT. (Vms .OR. Stabilize) .AND. UseBubbles ) THEN
        NBasis = 2*n
        Bubbles = .TRUE.
     END IF

!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IF ( Bubbles ) THEN
        IntegStuff = GaussPoints( element, Element % Type % GaussPoints2 )
     ELSE
        IntegStuff = GaussPoints( element )
     END IF
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!    Stabilization parameters: hK, mK (take a look at Franca et.al.)
!    If there is no convection term we don't need stabilization.
!------------------------------------------------------------------------------
       hK = element % hK
       mK = element % StabilizationMK

     ConvectAndStabilize = .FALSE.
     IF  ( Vms ) THEN
       NodalVelo(1,1:n) = Ux(1:n)
       NodalVelo(2,1:n) = Uy(1:n)
       NodalVelo(3,1:n) = Uz(1:n)
       NodalVelo(dim+1,1:n) = NodalPressure(1:n)

       dNodalBasisdx = 0._dp
       GradNodal = 0._dp
       DO p=1,n
          u = Element % Type % NodeU(p)
          v = Element % Type % NodeV(p)
          w = Element % Type % NodeW(p)
          stat = ElementInfo( Element, Nodes, u,v,w, detJ, Basis, dBasisdx )
 
          dNodalBasisdx(1:n,p,:) = dBasisdx(1:n,:)
          GradNodal(p,1:dim,1:dim) = MATMUL( NodalVelo(1:dim,1:n), dBasisdx(1:n,1:dim) )
          GradNodal(p,dim+1,1:dim) = MATMUL( NodalVelo(dim+1,1:n), dBasisdx(1:n,1:dim) )
       END DO

       NodalPvelo = 0._dp
       IF ( Transient ) THEN
         dt = CurrentModel % Solver % dt
         Order = MIN(CurrentModel % Solver % DoneTime,CurrentModel % Solver % Order)
 
         CALL GetVectorLocalSolution( NodalPVelo, 'Flow Solution', tStep=-1 )
 
         IF ( Order<2 ) THEN
           NodalPVelo(1:dim,1:n)=(NodalVelo(1:dim,1:n)-NodalPVelo(1:dim,1:n))/dt
         ELSE
           CALL GetVectorLocalSolution( Work, 'Flow Solution', tStep=-2 )
           NodalPVelo(1:dim,1:n)=(1.5_dp*NodalVelo(1:dim,1:n) - &
                  2._dp*NodalPVelo(1:dim,1:n)+0.5_dp*Work(1:dim,1:n) )/dt
         END IF
       END IF

       expc = GetCReal(GetMaterial(),'Heat Expansion Coefficient',Found)
       reft = GetCReal(GetMaterial(),'Reference Temperature',Found)
       CALL GetConstRealArray( GetConstants(), gWrk, 'Grav',Found )
       IF ( Found ) THEN
         grav = gWrk(1:3,1)*gWrk(4,1)
       ELSE
         grav    =  0.00_dp
         grav(2) = -9.81_dp
       END IF

       LC1 = 2._dp/mK
       LC(1,1:n) = Element % Type % NodeU(1:n)
       LC(2,1:n) = Element % Type % NodeV(1:n)
       LC(3,1:n) = Element % Type % NodeW(1:n)

       DO i=1,Element % Type % Dimension
         minl=MINVAL(LC(i,1:n))
         maxl=MAXVAL(LC(i,1:n))
         LC(i,1:n) = 2*(LC(i,1:n)-minl)/(maxl-minl)-1
       END DO
     ELSE IF ( Stabilize .AND. Convection ) THEN
       ConvectAndStabilize = .TRUE.
       hK = element % hK
       mK = element % StabilizationMK
       dNodalBasisdx = 0._dp
       DO p=1,n
         u = Element % Type % NodeU(p)
         v = Element % Type % NodeV(p)
         w = Element % Type % NodeW(p)
         stat = ElementInfo( Element, Nodes, u,v,w, detJ, Basis, dBasisdx )
         dNodalBasisdx(1:n,p,:) = dBasisdx(1:n,:)
       END DO
     END IF

!------------------------------------------------------------------------------
!    Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ

       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
             Basis,dBasisdx, Bubbles=Bubbles )

       s = detJ * S_Integ(t)
!------------------------------------------------------------------------------
!      Coefficient of the convection and time derivative terms
!      at the integration point
!------------------------------------------------------------------------------
       C1 = SUM( NodalC1(1:n) * Basis(1:n) )
       CT = SUM( NodalCT(1:n) * Basis(1:n) )

!------------------------------------------------------------------------------
!      Orientation parameters at the integration point
!------------------------------------------------------------------------------
       A1 = SUM( NodalK1(1:n) * Basis(1:n) ) 
       A2 = SUM( NodalK2(1:n) * Basis(1:n) )
       A3 = 1._dp - A1 - A2

       E1 = SUM( NodalEuler1(1:n) * Basis(1:n) )
       E2 = SUM( NodalEuler2(1:n) * Basis(1:n) )
       E3 = SUM( NodalEuler3(1:n) * Basis(1:n) )

!------------------------------------------------------------------------------
!      iota parameter at the integration point
!------------------------------------------------------------------------------
       iotaAtIp = SUM( iota(1:n) * Basis(1:n) )

!------------------------------------------------------------------------------
!      Material parameters at that point
!------------------------------------------------------------------------------
       ai(1)=A1
       ai(2)=A2
       ai(3)=A3
       ai(4)=E1
       ai(5)=E2
       ai(6)=E3

!------------------------------------------------------------------------------
!      Fourth order orientation tensor
!------------------------------------------------------------------------------
       call IBOF(ai,a4)

!------------------------------------------------------------------------------
!      Compute strainRate and Spin :
!------------------------------------------------------------------------------
       CSymmetry = CurrentCoordinateSystem() == AxisSymmetric
      
       StrainRate = 0.0
       Spin1 = 0.0
    
       LGrad = MATMUL( NodalVelocity(:,1:n), dBasisdx(1:n,:) )

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

       INDi(1:6) = (/ 1, 2, 3, 1, 2, 3 /)
       INDj(1:6) = (/ 1, 2, 3, 2, 3, 1 /)

       Do i=1,2*dim-3
         Spin(i)=Spin1(INDi(i+3),INDj(i+3))
       End do

!------------------------------------------------------------------------------
!       If not used set shear and rigid rotation to 0 for the fabric:
!------------------------------------------------------------------------------
       IF (.NOT. SpinShear) THEN
          Spin(1) = 0.0
          Spin(2) = 0.0
          Spin(3) = 0.0
          StrainRate(1,2) = 0.0
          StrainRate(2,3) = 0.0
          StrainRate(1,3) = 0.0
          StrainRate(3,1) = 0.0
          StrainRate(2,1) = 0.0
          StrainRate(3,2) = 0.0
       END IF

!------------------------------------------------------------------------------
!       Reaction coefficient:
!------------------------------------------------------------------------------
       SELECT CASE(comp)
        CASE(1)
         C0 = 2._dp*iotaAtIp*(StrainRate(1,1)-StrainRate(3,3)) 

        CASE(2)
         C0 = 2._dp*iotaAtIp*(StrainRate(2,2)-StrainRate(3,3))
       
        CASE(3)
         C0 = iotaAtIp*(StrainRate(1,1)+StrainRate(2,2)-2._dp*StrainRate(3,3))
        
        CASE(4)
         C0 = iotaAtIp*(StrainRate(2,2)-StrainRate(3,3))

        CASE(5)
         C0 = iotaAtIp*(StrainRate(1,1)-StrainRate(3,3)) 
       END SELECT

!------------------------------------------------------------------------------
!       Compute effective heatcapacity, if modelling phase change,
!       at the integration point.
!       NOTE: This is for heat equation only, not generally for diff.conv. equ.
!------------------------------------------------------------------------------
        IF ( PhaseChange ) THEN
          dEnth = 0.0D0
          dTemp = 0.0D0
          DO i=1,3
            dEnth = dEnth + SUM( Enthalpy(1:n) * dBasisdx(1:n,i) )**2
            dTemp = dTemp + SUM( NodalTemperature(1:n) * dBasisdx(1:n,i) )**2
          END DO

          CT = SQRT( dEnth/dTemp )
        END IF
!------------------------------------------------------------------------------
!      Coefficient of the diffusion term & it s derivatives at the
!      integration point
!------------------------------------------------------------------------------
       rho = SUM( Nodalrho(1:n) * Basis(1:n) ) 

       DO i=1,dim
         DO j=1,dim
           C2(i,j) = SUM( NodalC2(i,j,1:n) * Basis(1:n) )
         END DO
       END DO

       DO i=1,dim
          C2(i,i) = EffectiveConductivity( C2(i,i), rho, Element, &
                 NodalTemperature, UX,UY,UZ, Nodes, n, n, u, v, w )
       END DO
!------------------------------------------------------------------------------
!      If there's no convection term we don't need the velocities, and
!      also no need for stabilization
!------------------------------------------------------------------------------
       Convection = .FALSE.
       IF ( C1 /= 0.0D0 ) THEN
          Convection = .TRUE.
          IF ( PhaseChange ) C1 = CT
!------------------------------------------------------------------------------
!         Velocity from previous iteration at the integration point
!------------------------------------------------------------------------------
          Velo = 0.0D0
          Velo(1) = SUM( (UX(1:n)-MUX(1:n))*Basis(1:n) )
          Velo(2) = SUM( (UY(1:n)-MUY(1:n))*Basis(1:n) )
          IF ( dim > 2 ) Velo(3) = SUM( (UZ(1:n)-MUZ(1:n))*Basis(1:n) )

          IF ( Compressible ) THEN
            Grad = 0.0D0
            DO i=1,3
              Grad(1,i) = SUM( UX(1:n)*dBasisdx(1:n,i) )
              Grad(2,i) = SUM( UY(1:n)*dBasisdx(1:n,i) )
              IF ( dim > 2 ) Grad(3,i) = SUM( UZ(1:n)*dBasisdx(1:n,i) )
            END DO

            Pressure = SUM( NodalPressure(1:n)*Basis(1:n) )
            DivVelo = 0.0D0
            DO i=1,dim
              DivVelo = DivVelo + Grad(i,i)
            END DO
          END IF


          IF ( Vms ) THEN
            mu = GetCReal( GetMaterial(), 'Viscosity', Found )
            mu = EffectiveViscosity( mu, rho, Ux, Uy, Uz, &
                   Element, Nodes, n, n, u,v,w )

            Grad = 0.0D0
            DO i=1,3
              Grad(1,i) = SUM( Ux(1:n)*dBasisdx(1:n,i) )
              Grad(2,i) = SUM( Uy(1:n)*dBasisdx(1:n,i) )
              Grad(3,i) = SUM( Uz(1:n)*dBasisdx(1:n,i) )
            END DO
            VNorm = SQRT( SUM(Velo(1:dim)**2) )

            Temperature = SUM(Basis(1:n)*NodalTemperature(1:n))

            DO i=1,dim
              GradP(i) = SUM( NodalPressure(1:n)*dBasisdx(1:n,i) )
            END DO


            Gmat = 0._dp
            Gvec = 0._dp
            DO i=1,dim
              DO j=1,dim
                Gvec(i) = Gvec(i) + SUM(LC(j,1:n)*dBasisdx(1:n,i))
                DO k=1,dim
                  Gmat(i,j) = Gmat(i,j) + SUM(LC(k,1:n)*dBasisdx(1:n,i)) * &
                                          SUM(LC(k,1:n)*dBasisdx(1:n,j))
                END DO
              END DO
            END DO

            IF ( Transient ) THEN
              Tau_M = 1._dp / SQRT( SUM(Velo*MATMUL(Gmat,Velo)) + &
                    LC1**2 * (mu/rho)**2*SUM(Gmat*Gmat)/dim + 4/dt**2 )
            ELSE
              Tau_M = 1._dp / SQRT( SUM(Velo*MATMUL(Gmat,Velo)) + &
                    LC1**2 * (mu/rho)**2 * SUM(Gmat*Gmat)/dim )
            END IF

            Pe  = MIN( 1.0_dp, mK*hK*C1*VNorm/(2*ABS(C2(1,1))) )
            IF ( VNorm /= 0.0 ) THEN
               Tau = hK * Pe / (2 * C1 * VNorm)
            END IF

            RM = 0._dp
            DO p=1,n
              RM(p) = C0 * Basis(p)
              DO i=1,dim
                RM(p) = RM(p) + C1 * Velo(i) * dBasisdx(p,i)
                DO j=1,dim
                  RM(p) = RM(p) - C2(i,j)*SUM(dNodalBasisdx(p,1:n,i)*dBasisdx(1:n,j))
                END DO
              END DO
            END DO

            VRM = 0._dp
            DO i=1,dim
              VRM(i) = SUM(NodalPVelo(i,1:n)*Basis(1:n))
              DO j=1,dim
                VRM(i) = VRM(i) + Velo(j) * Grad(i,j)
                VRM(i) = VRM(i) - (mu/rho)*SUM(GradNodal(1:n,i,j)*dBasisdx(1:n,j))
              END DO
              VRM(i) = VRM(i) + GradP(i)
              VRM(i) = VRM(i) + Grav(i)*ExpC*( Temperature - RefT )
            END DO
          ELSE IF ( Stabilize ) THEN
!------------------------------------------------------------------------------
!           Stabilization parameter Tau
!------------------------------------------------------------------------------
            VNorm = SQRT( SUM(Velo(1:dim)**2) )

            Pe  = MIN( 1.0D0, mK*hK*C1*VNorm/(2*ABS(C2(1,1))) )
            Tau = 0.0D0
            IF ( VNorm /= 0.0 ) THEN
               Tau = hK * Pe / (2 * C1 * VNorm)
            END IF
!------------------------------------------------------------------------------

            DO i=1,dim
              DO j=1,dim
                DO k=1,dim
                  dC2dx(i,j,k) = SUM( NodalC2(i,j,1:n)*dBasisdx(1:n,k) )
                END DO
              END DO
            END DO

!------------------------------------------------------------------------------
!           Compute residual & stablization vectors
!------------------------------------------------------------------------------
            DO p=1,N
              SU(p) = C0 * Basis(p)
              DO i = 1,dim
                SU(p) = SU(p) + C1 * dBasisdx(p,i) * Velo(i)
                DO j=1,dim
                  SU(p) = SU(p) - dC2dx(i,j,j) * dBasisdx(p,i)
                  SU(p) = SU(p) - C2(i,j) * SUM(dNodalBasisdx(p,1:n,i)*dBasisdx(1:n,j))
                END DO
              END DO

              SW(p) = C0 * Basis(p)
              DO i = 1,dim
                SW(p) = SW(p) + C1 * dBasisdx(p,i) * Velo(i)
                DO j=1,dim
                  SW(p) = SW(p) - dC2dx(i,j,j) * dBasisdx(p,i)
                  SW(p) = SW(p) - C2(i,j) * SUM(dNodalBasisdx(p,1:n,i)*dBasisdx(1:n,j))
                END DO
              END DO
            END DO
          END IF
        END IF

!------------------------------------------------------------------------------
!       Loop over basis functions of both unknowns and weights
!------------------------------------------------------------------------------
        DO p=1,NBasis
        DO q=1,NBasis
          A = 0.0D00
!------------------------------------------------------------------------------
!         The diffusive-convective equation without stabilization
!------------------------------------------------------------------------------
          M = CT * Basis(q) * Basis(p)
          A = A + C0 * Basis(q) * Basis(p)
!------------------------------------------------------------------------------
!         The diffusion term
!------------------------------------------------------------------------------
          DO i=1,dim
            DO j=1,dim
              A = A + C2(i,j) * dBasisdx(q,i) * dBasisdx(p,j)
            END DO
          END DO

          IF ( Convection ) THEN
!------------------------------------------------------------------------------
!           The convection term
!------------------------------------------------------------------------------
            DO i=1,dim
              A = A + C1 * Velo(i) * dBasisdx(q,i) * Basis(p)
            END DO
!------------------------------------------------------------------------------
!           Next we add the stabilization...
!------------------------------------------------------------------------------
            IF ( Vms ) THEN
              DO i=1,dim
                A = A - C1 * Tau_M * VRM(i) * dBasisdx(q,i) * Basis(p)

                A = A + C1 * Velo(i) * Tau * RM(q) * dBasisdx(p,i)
                M = M + C1 * Velo(i) * Tau * CT*Basis(q) * dBasisdx(p,i)

                A = A - C1 * Tau_M * VRM(i) * Tau * RM(q) * dBasisdx(p,i)
                M = M - C1 * Tau_M * VRM(i) * Tau * CT*Basis(q) * dBasisdx(p,i)
              END DO
            ELSE IF ( Stabilize ) THEN
              A = A + Tau * SU(q) * SW(p)
              M = M + Tau * CT * Basis(q) * SW(p)
            END IF
          END IF

          StiffMatrix(p,q) = StiffMatrix(p,q) + s * A
          MassMatrix(p,q)  = MassMatrix(p,q)  + s * M
        END DO
        END DO

!------------------------------------------------------------------------------
!       The righthand side...
!------------------------------------------------------------------------------
        SELECT CASE(comp)
         
         CASE(1)
          IF(dim == 3) THEN 
            LoadAtIp = 2._dp*(Spin(1)*E1-Spin(3)*E3)+ &
                 2._dp*iotaAtIp*(StrainRate(1,2)*(2._dp*a4(7)-E1)+ &
                 StrainRate(3,1)*(2._dp*a4(6)-E3)+a4(1)*StrainRate(1,1)+ &
                 a4(3)*StrainRate(2,2)-(a4(1)+a4(3))*StrainRate(3,3)+&
                 2._dp*a4(4)*StrainRate(2,3))
          ELSE
              LoadAtIp = 2._dp*Spin(1)*E1+ &
                 2._dp*iotaAtIp*(StrainRate(1,2)*(2._dp*a4(7)-E1)+ &
                 a4(1)*StrainRate(1,1)+a4(3)*StrainRate(2,2)- &
                 (a4(1)+a4(3))*StrainRate(3,3)) 
          END IF
         
         CASE(2)
          IF(dim == 3) THEN
            LoadAtIp = 2._dp*(-Spin(1)*E1+Spin(2)*E2)+ &
                 2._dp*iotaAtIp*(StrainRate(1,2)*(2._dp*a4(9)-E1)+ &
                 StrainRate(2,3)*(2._dp*a4(8)-E2)+a4(3)*StrainRate(1,1)+ &
                 a4(2)*StrainRate(2,2)-(a4(2)+a4(3))*StrainRate(3,3)+&
                 2._dp*a4(5)*StrainRate(3,1)) 
          ELSE
            LoadAtIp = -2._dp*Spin(1)*E1+ &
                 2._dp*iotaAtIp*(StrainRate(1,2)*(2._dp*a4(9)-E1)+ &
                 a4(3)*StrainRate(1,1)+a4(2)*StrainRate(2,2)- &
                 (a4(2)+a4(3))*StrainRate(3,3)) 
          END IF

         CASE(3)
          IF(dim == 3) THEN
            LoadAtIp = Spin(1)*(A2-A1)-Spin(3)*E2+Spin(2)*E3+ &
                 StrainRate(1,2)*iotaAtIp*(4._dp*a4(3)-(A1+A2))+ &
                 StrainRate(3,1)*iotaAtIp*(4._dp*a4(4)-E2)+ &
                 StrainRate(2,3)*iotaAtIp*(4._dp*a4(5)-E3)+ &
                 2._dp*iotaAtIp*(a4(7)*StrainRate(1,1)+a4(9)*StrainRate(2,2)- &
                 (a4(7)+a4(9))*StrainRate(3,3))  
          ELSE
              LoadAtIp = Spin(1)*(A2-A1)+ &
                 StrainRate(1,2)*iotaAtIp*(4._dp*a4(3)-(A1+A2))+ &
                 2._dp*iotaAtIp*(a4(7)*StrainRate(1,1)+a4(9)*StrainRate(2,2)- &
                 (a4(7)+a4(9))*StrainRate(3,3))
          END IF

         CASE(4)
          LoadAtIp = Spin(2)*(A3-A2)+Spin(3)*E1-Spin(1)*E3+ &
                 StrainRate(1,2)*iotaAtIP*(4._dp*a4(5)-E3)+ &
                 StrainRate(2,3)*iotaAtIP*(3._dp*A2-A3-4._dp*(a4(2)+a4(3)))+ &
                 StrainRate(3,1)*iotaAtIp*(3._dp*E1-4._dp*(a4(7)+a4(9)))+ &
                 2._dp*iotaAtIP*(a4(4)*StrainRate(1,1)+a4(8)*StrainRate(2,2)- &
                 (a4(4)+a4(8))*StrainRate(3,3))
         
         CASE(5)
          LoadAtIp = Spin(3)*(A1-A3)+Spin(1)*E2-Spin(2)*E1+ &
                 StrainRate(1,2)*iotaAtIP*(4._dp*a4(4)-E2) + &
                 StrainRate(3,1)*iotaAtIp*(3._dp*A1-A3-4._dp*(a4(1)+a4(3)))+ &
                 StrainRate(2,3)*iotaAtIp*(3._dp*E1-4._dp*(a4(7)+a4(9)))+ &
                 2._dp*iotaAtIp*(a4(6)*StrainRate(1,1)+a4(5)*StrainRate(2,2)- &
                 (a4(6)+a4(5))*StrainRate(3,3))

        END SELECT

!------------------------------------------------------------------------------
!       Force at the integration point
!------------------------------------------------------------------------------
        Force = LoadAtIp + &
          JouleHeat( Element, Nodes, u, v, w, n )

        IF ( Convection ) THEN
!         IF ( Compressible ) Force = Force - Pressure * DivVelo

          Pcoeff = SUM(NodalPressureCoeff(1:n)*Basis(1:n))
          IF ( Pcoeff /= 0.0_dp ) THEN
            Force = Force + Pcoeff * SUM(NodalDPressureDt(1:n)*Basis(1:n))
            DO i=1,dim
              Force = Force + Pcoeff*Velo(i)*SUM(NodalPressure(1:n)*dBasisdx(1:n,i))
            END DO
          END IF

          mu = SUM( Nodalmu(1:n) * Basis(1:n) )
          mu = EffectiveViscosity( mu, rho, Ux, Uy, Uz, &
                 Element, Nodes, n, n, u,v,w )
          IF ( mu > 0.0d0 ) THEN
            IF ( .NOT.Compressible ) THEN
              Grad = 0.0D0
              DO i=1,3
                Grad(1,i) = SUM( UX(1:n)*dBasisdx(1:n,i) )
                Grad(2,i) = SUM( UY(1:n)*dBasisdx(1:n,i) )
                IF ( dim > 2 ) Grad(3,i) = SUM( UZ(1:n)*dBasisdx(1:n,i) )
              END DO
            END IF
            Force = Force + 0.5d0 * mu*SecondInvariant(Velo,Grad)
          END IF
        END IF
!------------------------------------------------------------------------------
        DO p=1,NBasis
          Load = Force * Basis(p)
          IF ( Vms ) THEN
            DO i=1,dim
              Load = Load + C1 * Velo(i) * Tau  * Force * dBasisdx(p,i)
              Load = Load - C1 * Tau_M * VRM(i) * Tau * Force * dBasisdx(p,i)
            END DO
          ELSE IF ( ConvectAndStabilize ) THEN
            Load = Load + Tau * Force * SW(p)
          END IF
          ForceVector(p) = ForceVector(p) + s * Load
        END DO
      END DO
!------------------------------------------------------------------------------
   END SUBROUTINE DiffuseConvectiveComposeCaffe
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE DiffuseConvectiveBoundaryCaffe( BoundaryMatrix,BoundaryVector, &
               LoadVector,NodalAlpha,Element,n,Nodes )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Return element local matrices and RSH vector for boundary conditions
!  of diffusion convection equation: 
!
!  ARGUMENTS:
!
!  REAL(KIND=dp) :: BoundaryMatrix(:,:)
!     OUTPUT: coefficient matrix if equations
!
!  REAL(KIND=dp) :: BoundaryVector(:)
!     OUTPUT: RHS vector
!
!  REAL(KIND=dp) :: LoadVector(:)
!     INPUT: coefficient of the force term
!
!  REAL(KIND=dp) :: NodalAlpha
!     INPUT: coefficient for temperature dependent term
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!   INTEGER :: n
!       INPUT: Number  of element nodes
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!
!******************************************************************************

     REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:), &
                    LoadVector(:),NodalAlpha(:)

     TYPE(Nodes_t)   :: Nodes
     TYPE(Element_t) :: Element

     INTEGER :: n

     REAL(KIND=dp) :: ddBasisddx(n,3,3)
     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3),detJ

     REAL(KIND=dp) :: U,V,W,S
     REAL(KIND=dp) :: Force,Alpha
     REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

     INTEGER :: i,t,q,p,N_Integ

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat
!------------------------------------------------------------------------------

     BoundaryVector = 0.0D0
     BoundaryMatrix = 0.0D0
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( Element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ
       U = U_Integ(t)
       V = V_Integ(t)
       W = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
                  Basis,dBasisdx )

       S = detJ * S_Integ(t)
!------------------------------------------------------------------------------
       Force = SUM( LoadVector(1:n)*Basis )
       Alpha = SUM( NodalAlpha(1:n)*Basis )

       DO p=1,N
         DO q=1,N
           BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + &
              s * Alpha * Basis(q) * Basis(p)
         END DO
       END DO

       DO q=1,N
         BoundaryVector(q) = BoundaryVector(q) + s * Basis(q) * Force
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE DiffuseConvectiveBoundaryCaffe
!------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------------------------
 SUBROUTINE FreeBoundaryCaffe2( Comp, MASS, STIFF, FORCE, LOAD, &
          NodalK1, NodalK2, NodalEuler1, NodalEuler2, NodalEuler3, & 
          NodalVelo, NodMeshVel, Element, n, Nodes, ParentElement, ParentNodes, np, iota, SpinShear)
!-----------------------------------------------------------------------------------------------------------------------------------------------------

     REAL(KIND=dp) :: STIFF(:,:),MASS(:,:), FORCE(:)
     REAL(KIND=dp) :: LOAD(:), NodalVelo(:,:),NodMeshVel(:,:)
     REAL(KIND=dp), DIMENSION(:) :: NodalK1, NodalK2, NodalEuler1, &
               NodalEuler2, NodalEuler3, iota

     TYPE(Nodes_t) :: Nodes, ParentNodes
     TYPE(Element_t), POINTER :: Element, ParentElement
     INTEGER :: n, Comp, np
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(n),ddBasisddx(n,3,3)
     
     REAL(KIND=dp) :: dBasisdx(n,3), detJ
     REAL(KIND=dp) :: ParentBasis(np), ParentdBasisdx(np,3), ParentddBasisddx(np,3,3)

     REAL(KIND=dp) :: A1, A2, A3, E1, E2, E3

     REAL(KIND=dp) :: A, M, C0, epsi, Radius
     REAL(KIND=dp) :: ai(6),a4(9)

     INTEGER :: i,j,k,p,q,t,dim,ind(3)

     REAL(KIND=dp) :: s,u,v,w
     REAL(KIND=dp) :: Velo(3),Spin(3), divu, LoadAtIp, iotaAtIp

     REAL(KIND=dp) :: LGrad(3,3),StrainRate(3,3)
     REAL(KIND=dp) :: Spin1(3,3)
     Integer :: INDi(6),INDj(6)
     LOGICAL :: SpinShear, CSymmetry
     
     LOGICAL :: stat
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     INTERFACE
        
        SUBROUTINE IBOF(ai,a4)
           USE Types
           REAL(KIND=dp),intent(in) :: ai(6)
           REAL(KIND=dp),intent(out) :: a4(9)
        END SUBROUTINE
 
     END INTERFACE
!------------------------------------------------------------------------------
       
      dim = CoordinateSystemDimension()

      FORCE = 0.0D0
      MASS  = 0.0D0
      STIFF = 0.0D0

!----------------------
!    Integration stuff:
!----------------------
      IntegStuff = GaussPoints( Element  )

!----------------------------
!   Now we start integrating:
!----------------------------
      DO t=1,IntegStuff % n

      u = IntegStuff % u(t)
      v = IntegStuff % v(t)
      w = IntegStuff % w(t)
      S = IntegStuff % s(t) 

!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
               Basis, dBasisdx, ddBasisddx, .FALSE. )

      S = S * detJ

      CALL FindParentUVW( Element, n, ParentElement, np, U, V, W, Basis )
      stat = ElementInfo( ParentElement, ParentNodes, U, V, W, &
           detJ,  ParentBasis, ParentdBasisdx, ParentddBasisddx, .FALSE. )
      
!------------------------------------------------------
!     Orientation parameters at the integration point:
!------------------------------------------------------
      A1 = SUM( NodalK1(1:n) * Basis(1:n) ) 
      A2 = SUM( NodalK2(1:n) * Basis(1:n) )
      A3 = 1._dp - A1 - A2

      E1 = SUM( NodalEuler1(1:n) * Basis(1:n) )
      E2 = SUM( NodalEuler2(1:n) * Basis(1:n) )
      E3 = SUM( NodalEuler3(1:n) * Basis(1:n) )

!-----------------------------------------------------
!     iota parameter at the integration point:
!----------------------------------------------------
     iotaAtIp = SUM( iota(1:n) * Basis(1:n) )

!--------------------------------------
!      Strain-Rate, Stresses and Spin:
!-------------------------------------

      CSymmetry = CurrentCoordinateSystem() == AxisSymmetric
      
      StrainRate = 0.0
      Spin1 = 0.0

!-------------------------------------
!    Material parameters at that point
!------------------------------------
      ai(1)=A1
      ai(2)=A2
      ai(3)=A3
      ai(4)=E1
      ai(5)=E2
      ai(6)=E3

!------------------------------------
!    fourth order orientation tensor:
!------------------------------------
      call IBOF(ai,a4)

!----------------------------------
!    Compute strainRate and Spin :
!----------------------------------

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
      

      INDi(1:6) = (/ 1, 2, 3, 1, 2, 3 /)
      INDj(1:6) = (/ 1, 2, 3, 2, 3, 1 /)

   
      Do i=1,2*dim-3
        Spin(i)=Spin1(INDi(i+3),INDj(i+3))
      End do

!------------------------------------------------------------------------------
!    If not used set shear and rigid rotation to 0 for the fabric:
!------------------------------------------------------------------------------
      IF (.NOT. SpinShear) THEN
         Spin(1) = 0.0
         Spin(2) = 0.0
         Spin(3) = 0.0
         StrainRate(1,2) = 0.0
         StrainRate(2,3) = 0.0
         StrainRate(1,3) = 0.0
         StrainRate(3,1) = 0.0
         StrainRate(2,1) = 0.0
         StrainRate(3,2) = 0.0
      END IF

!---------------
!    Velocity :
! -------------
      Velo = 0.0d0
      divu = 0.0d0
      DO i=1,dim
         Velo(i) = SUM( Basis(1:n) * (NodalVelo(i,1:n) - NodMeshVel(i,1:n)) )
      END DO

!--------------------------
!     Reaction coefficient:
!--------------------------
      SELECT CASE(comp)
      CASE(1)
        C0 = 2._dp*iotaAtIp*(StrainRate(1,1)-StrainRate(3,3)) 

      CASE(2)
        C0 = 2._dp*iotaAtIp*(StrainRate(2,2)-StrainRate(3,3))
       
      CASE(3)
        C0 = iotaAtIp*(StrainRate(1,1)+StrainRate(2,2)-2._dp*StrainRate(3,3))
        
      CASE(4)
        C0 = iotaAtIp*(StrainRate(2,2)-StrainRate(3,3))

      CASE(5)
        C0 = iotaAtIp*(StrainRate(1,1)-StrainRate(3,3)) 

      END SELECT
      
!--------------------------------------------------------------
!     Loop over basis functions (of both unknowns and weights):
!--------------------------------------------------------------
      DO p=1,np
         DO q=1,np
            A = 0.0d0
            M = ParentBasis(p) * ParentBasis(q)

            !---------------
            !Reaction terms:
            !---------------
            A = A + C0 * ParentBasis(q) * ParentBasis(p)

            !-----------------
            ! Advection terms:
            !-----------------
            DO j=1,dim
               A = A - Velo(j) * ParentBasis(q) * ParentdBasisdx(p,j)
            END DO

            !-----------------------------------
            !Add nodal matrix to element matrix:
            !-----------------------------------
            MASS( p,q )  = MASS( p,q )  + S * M
            STIFF( p,q ) = STIFF( p,q ) + S * A
         END DO
      END DO

!----------------------------
!     The righthand side.:
!----------------------------
         SELECT CASE(comp)
         
         CASE(1)
          IF(dim == 3) THEN 
            LoadAtIp = 2._dp*(Spin(1)*E1-Spin(3)*E3)+ &
                 2._dp*iotaAtIp*(StrainRate(1,2)*(2._dp*a4(7)-E1)+ &
                 StrainRate(3,1)*(2._dp*a4(6)-E3)+a4(1)*StrainRate(1,1)+ &
                 a4(3)*StrainRate(2,2)-(a4(1)+a4(3))*StrainRate(3,3)+&
                 2._dp*a4(4)*StrainRate(2,3))
          ELSE
              LoadAtIp = 2._dp*Spin(1)*E1+ &
                 2._dp*iotaAtIp*(StrainRate(1,2)*(2._dp*a4(7)-E1)+ &
                 a4(1)*StrainRate(1,1)+a4(3)*StrainRate(2,2)- &
                 (a4(1)+a4(3))*StrainRate(3,3)) 
          END IF
         
         CASE(2)
          IF(dim == 3) THEN
            LoadAtIp = 2._dp*(-Spin(1)*E1+Spin(2)*E2)+ &
                 2._dp*iotaAtIp*(StrainRate(1,2)*(2._dp*a4(9)-E1)+ &
                 StrainRate(2,3)*(2._dp*a4(8)-E2)+a4(3)*StrainRate(1,1)+ &
                 a4(2)*StrainRate(2,2)-(a4(2)+a4(3))*StrainRate(3,3)+&
                 2._dp*a4(5)*StrainRate(3,1)) 
          ELSE
            LoadAtIp = -2._dp*Spin(1)*E1+ &
                 2._dp*iotaAtIp*(StrainRate(1,2)*(2._dp*a4(9)-E1)+ &
                 a4(3)*StrainRate(1,1)+a4(2)*StrainRate(2,2)- &
                 (a4(2)+a4(3))*StrainRate(3,3)) 
          END IF

         CASE(3)
          IF(dim == 3) THEN
            LoadAtIp = Spin(1)*(A2-A1)-Spin(3)*E2+Spin(2)*E3+ &
                 StrainRate(1,2)*iotaAtIp*(4._dp*a4(3)-(A1+A2))+ &
                 StrainRate(3,1)*iotaAtIp*(4._dp*a4(4)-E2)+ &
                 StrainRate(2,3)*iotaAtIp*(4._dp*a4(5)-E3)+ &
                 2._dp*iotaAtIp*(a4(7)*StrainRate(1,1)+a4(9)*StrainRate(2,2)- &
                 (a4(7)+a4(9))*StrainRate(3,3))  
          ELSE
              LoadAtIp = Spin(1)*(A2-A1)+ &
                 StrainRate(1,2)*iotaAtIp*(4._dp*a4(3)-(A1+A2))+ &
                 2._dp*iotaAtIp*(a4(7)*StrainRate(1,1)+a4(9)*StrainRate(2,2)- &
                 (a4(7)+a4(9))*StrainRate(3,3))
          END IF

         CASE(4)
          LoadAtIp = Spin(2)*(A3-A2)+Spin(3)*E1-Spin(1)*E3+ &
                 StrainRate(1,2)*iotaAtIP*(4._dp*a4(5)-E3)+ &
                 StrainRate(2,3)*iotaAtIP*(3._dp*A2-A3-4._dp*(a4(2)+a4(3)))+ &
                 StrainRate(3,1)*iotaAtIp*(3._dp*E1-4._dp*(a4(7)+a4(9)))+ &
                 2._dp*iotaAtIP*(a4(4)*StrainRate(1,1)+a4(8)*StrainRate(2,2)- &
                 (a4(4)+a4(8))*StrainRate(3,3))
         
         CASE(5)
          LoadAtIp = Spin(3)*(A1-A3)+Spin(1)*E2-Spin(2)*E1+ &
                 StrainRate(1,2)*iotaAtIP*(4._dp*a4(4)-E2) + &
                 StrainRate(3,1)*iotaAtIp*(3._dp*A1-A3-4._dp*(a4(1)+a4(3)))+ &
                 StrainRate(2,3)*iotaAtIp*(3._dp*E1-4._dp*(a4(7)+a4(9)))+ &
                 2._dp*iotaAtIp*(a4(6)*StrainRate(1,1)+a4(5)*StrainRate(2,2)- &
                 (a4(6)+a4(5))*StrainRate(3,3))

         END SELECT

        FORCE(1:np) = FORCE(1:np) + S*LoadAtIp*ParentBasis(1:np)
   END DO 

!------------------------------------------------------------------------------
    END SUBROUTINE FreeBoundaryCaffe2
!------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------------
   SUBROUTINE FreeBoundaryCaffe( BoundaryMatrix, BoundaryVector, &
         LoadVector, Element, n , Nodes, ParentElement, np, ParentNodes, Velo, MeshVelo )
!--------------------------------------------------------------------------------------------------
     REAL(KIND=dp) :: BoundaryMatrix(:,:),  BoundaryVector(:), LoadVector(:), Velo(:,:),MeshVelo(:,:)
     INTEGER :: n, np
     TYPE(Element_t), POINTER :: Element, ParentElement
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3)
     REAL(KIND=dp) :: ParentBasis(np), ParentdBasisdx(np,3), ParentddBasisddx(np,3,3)
     INTEGER :: i,j,p,q,t,dim

     REAL(KIND=dp) :: Normal(3), g, L, Udotn, cu(3), detJ,U,V,W,S
     LOGICAL :: Stat
     TYPE(GaussIntegrationPoints_t) :: IntegStuff

     TYPE(Nodes_t) :: Nodes, ParentNodes
!------------------------------------------------------------------------------
     
     dim = CoordinateSystemDimension()
     BoundaryMatrix = 0.0d0
     BoundaryVector = 0.0d0

!     CALL GetElementNodes( Nodes, Element )
!     CALL GetElementNodes( ParentNodes, ParentElement )

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

       L = SUM( LoadVector(1:n) * Basis(1:n) )
       cu = 0.0d0
       DO i=1,dim
          cu(i) = SUM( (Velo(i,1:n)-MeshVelo(i,1:n)) * Basis(1:n) )
       END DO
       Udotn = SUM( Normal * cu )

       DO p = 1,np
         DO q=1,np
            BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + s*Udotn*ParentBasis(q)*ParentBasis(p)
         END DO
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE FreeBoundaryCaffe
!------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------
    SUBROUTINE FindParentUVW( Edge, nEdge, Parent, nParent, U, V, W, Basis )
!---------------------------------------------------------------------------------------------------------
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

!---------------------------------------------------------------------------------------------
!       Compute fourth order tensor a4 from a2 with closure function IBOF (Chung, 2002)
!       a2 enters in the order : 11, 22, 33, 12, 23 ,13
!       Output for a4 is in the order : 1111, 2222, 1122, 1123, 2231, 1131, 1112, 2223, 2212
!       Code modified from Gillet-Chaullet source
!---------------------------------------------------------------------------------------------
   SUBROUTINE IBOF(a2,a4)
!---------------------------------------------------------------------------------------------

     USE Types
       
     implicit none

     Real(dp),dimension(6),intent(in):: a2  
     Real(dp),dimension(9),intent(out):: a4  
     Real(dp):: a_11,a_22,a_33,a_12,a_13,a_23
     Real(dp):: b_11,b_22,b_12,b_13,b_23
     Real(dp):: aPlusa

     Real(dp),dimension(21) :: vec
     Real(dp),dimension(3,21) :: Mat
     Real(dp),dimension(6) :: beta
     Real(dp) :: Inv2,Inv3
     integer :: i,j      
     
     a_11=a2(1)
     a_22=a2(2)
     a_33=a2(3)
     a_12=a2(4)
     a_23=a2(5)
     a_13=a2(6)

     !--------------------------------------------------
     !   Coefficiants 
     !--------------------------------------------------
     Mat(1,1)=0.217774509809788e+02_dp
     Mat(1,2)=-.297570854171128e+03_dp
     Mat(1,3)=0.188686077307885e+04_dp
     Mat(1,4)=-.272941724578513e+03_dp
     Mat(1,5)=0.417148493642195e+03_dp
     Mat(1,6)=0.152038182241196e+04_dp
     Mat(1,7)=-.137643852992708e+04_dp
     Mat(1,8)=-.628895857556395e+03_dp
     Mat(1,9)=-.526081007711996e+04_dp
     Mat(1,10)=-.266096234984017e+03_dp
     Mat(1,11)=-.196278098216953e+04_dp
     Mat(1,12)=-.505266963449819e+03_dp
     Mat(1,13)=-.110483041928547e+03_dp
     Mat(1,14)=0.430488193758786e+04_dp
     Mat(1,15)=-.139197970442470e+02_dp
     Mat(1,16)=-.144351781922013e+04_dp
     Mat(1,17)=-.265701301773249e+03_dp
     Mat(1,18)=-.428821699139210e+02_dp
     Mat(1,19)=-.443236656693991e+01_dp
     Mat(1,20)=0.309742340203200e+04_dp
     Mat(1,21)=0.386473912295113e+00_dp
     Mat(2,1)=-.514850598717222e+00_dp
     Mat(2,2)=0.213316362570669e+02_dp
     Mat(2,3)=-.302865564916568e+03_dp
     Mat(2,4)=-.198569416607029e+02_dp
     Mat(2,5)=-.460306750911640e+02_dp
     Mat(2,6)=0.270825710321281e+01_dp
     Mat(2,7)=0.184510695601404e+03_dp
     Mat(2,8)=0.156537424620061e+03_dp
     Mat(2,9)=0.190613131168980e+04_dp
     Mat(2,10)=0.277006550460850e+03_dp
     Mat(2,11)=-.568117055198608e+02_dp
     Mat(2,12)=0.428921546783467e+03_dp
     Mat(2,13)=0.142494945404341e+03_dp
     Mat(2,14)=-.541945228489881e+04_dp
     Mat(2,15)=0.233351898912768e+02_dp
     Mat(2,16)=0.104183218654671e+04_dp
     Mat(2,17)=0.331489412844667e+03_dp
     Mat(2,18)=0.660002154209991e+02_dp
     Mat(2,19)=0.997500770521877e+01_dp
     Mat(2,20)=0.560508628472486e+04_dp
     Mat(2,21)=0.209909225990756e+01_dp
     Mat(3,1)=0.203814051719994e+02_dp
     Mat(3,2)=-.283958093739548e+03_dp
     Mat(3,3)=0.173908241235198e+04_dp
     Mat(3,4)=-.195566197110461e+03_dp
     Mat(3,5)=-.138012943339611e+03_dp
     Mat(3,6)=0.523629892715050e+03_dp
     Mat(3,7)=0.859266451736379e+03_dp
     Mat(3,8)=-.805606471979730e+02_dp
     Mat(3,9)=-.468711180560599e+04_dp
     Mat(3,10)=0.889580760829066e+01_dp
     Mat(3,11)=-.782994158054881e+02_dp
     Mat(3,12)=-.437214580089117e+02_dp
     Mat(3,13)=0.112996386047623e+01_dp
     Mat(3,14)=0.401746416262936e+04_dp
     Mat(3,15)=0.104927789918320e+01_dp
     Mat(3,16)=-.139340154288711e+03_dp
     Mat(3,17)=-.170995948015951e+02_dp
     Mat(3,18)=0.545784716783902e+00_dp
     Mat(3,19)=0.971126767581517e+00_dp
     Mat(3,20)=0.141909512967882e+04_dp
     Mat(3,21)=0.994142892628410e+00_dp

     !--------------------------------------------------  
     !   Compute the invariants
     !--------------------------------------------------
     Inv2=0.5_dp*(1._dp-(a_11*a_11+a_22*a_22+a_33*a_33+ &
            2._dp*(a_12*a_12+a_13*a_13+a_23*a_23)))
            
     Inv3=a_11*(a_22*a_33-a_23*a_23)+a_12*(a_23*a_13-a_12*a_33)+ &
             a_13*(a_12*a_23-a_22*a_13)
     
     !----------------------------------------------------  
     !   Complete polynome of degree 5 for the 2 invariants.
     !---------------------------------------------------- 
     vec(1)=1._dp
     vec(2)=Inv2
     vec(3)=vec(2)*vec(2)
     vec(4)=Inv3
     vec(5)=vec(4)*vec(4)
     vec(6)=vec(2)*vec(4)
     vec(7)=vec(3)*vec(4)
     vec(8)=vec(2)*vec(5)
     vec(9)=vec(2)*vec(3)
     vec(10)=vec(5)*vec(4)
     vec(11)=vec(9)*vec(4)
     vec(12)=vec(3)*vec(5)
     vec(13)=vec(2)*vec(10)
     vec(14)=vec(3)*vec(3)
     vec(15)=vec(5)*vec(5)
     vec(16)=vec(14)*vec(4)
     vec(17)=vec(12)*vec(2)
     vec(18)=vec(12)*vec(4)
     vec(19)=vec(2)*vec(15)
     vec(20)=vec(14)*vec(2)
     vec(21)=vec(15)*vec(4)

      !-------------------------------------------------------------------------------
      !  Compites beta_bar (cf annexe C Chung)
      !  Warning: beta(1)=beta_bar_3 (Chung); beta(2)=beta_bar_4; beta(3)=beta_bar_6
      !           beta(4)=beta_bar_1        ; beta(5)=beta_bar_2; beta(6)=beta_bar_5
      !-------------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------------
      !  Calcul the three betas in terms of the polynomes
      !-------------------------------------------------------------------------------
      beta(:)=0._dp
      Do i=1,3
        Do j=1,21
           beta(i)=beta(i)+Mat(i,j)*vec(j)
        End do
      End do
       
       !-------------------------------------------------------------------------------   
       ! Calcul the other 3 to get the normalisation
       !-------------------------------------------------------------------------------
       beta(4)=3._dp*(-1._dp/7._dp+beta(1)* &
            (1._dp/7._dp+4._dp*Inv2/7._dp+8._dp*Inv3/3._dp)/5._dp- &
       
       beta(2)*(0.2_dp-8._dp*Inv2/15._dp-14._dp*Inv3/15._dp)- &
            beta(3)*(1._dp/35._dp-24._dp*Inv3/105._dp-4._dp*Inv2/35._dp+ &
            16._dp*Inv2*Inv3/15._dp+8._dp*Inv2*Inv2/35._dp))/5._dp

       beta(5)=6._dp*(1._dp-0.2_dp*beta(1)*(1._dp+4._dp*Inv2)+ &
            7._dp*beta(2)*(1._dp/6._dp-Inv2)/5._dp- &
            beta(3)*(-0.2_dp+2._dp*Inv3/3._dp+4._dp*Inv2/5._dp- &
            8._dp*Inv2*Inv2/5._dp))/7._dp

       beta(6)=-4._dp*beta(1)/5._dp-7._dp*beta(2)/5._dp- &
            6._dp*beta(3)*(1._dp-4._dp*Inv2/3._dp)/5._dp

       !----------------------------
       ! Beta_bar
       !----------------------------
       Do i=1,6
         beta(i)=beta(i)/3._dp
       End do
       beta(2)=beta(2)/2._dp
       beta(5)=beta(5)/2._dp
       beta(6)=beta(6)/2._dp

       !----------------------------
       ! Compute 5 b=a.a
       !----------------------------
       b_11=a_11*a_11+a_12*a_12+a_13*a_13
       b_22=a_22*a_22+a_12*a_12+a_23*a_23
       b_12=a_11*a_12+a_12*a_22+a_13*a_23
       b_13=a_11*a_13+a_12*a_23+a_13*a_33
       b_23=a_12*a_13+a_22*a_23+a_23*a_33

       !----------------------------
       ! Compute the  9 terms of a4
       !----------------------------
       a4(1)=3._dp*beta(4)+6._dp*beta(5)*a_11+3._dp*beta(1)*a_11*a_11+&
         6._dp*beta(2)*b_11+6._dp*beta(6)*a_11*b_11+3._dp*beta(3)*b_11*b_11
         a4(2)=3._dp*beta(4)+6._dp*beta(5)*a_22+3._dp*beta(1)*a_22*a_22+&
         6._dp*beta(2)*b_22+6._dp*beta(6)*a_22*b_22+3._dp*beta(3)*b_22*b_22

       a4(3)=beta(4)+beta(5)*(a_22+a_11)+beta(1)*(a_11*a_22+2._dp*a_12*a_12)+&
         beta(2)*(b_22+b_11)+beta(6)*(a_11*b_22+a_22*b_11+4._dp*a_12*b_12)+&
         beta(3)*(b_11*b_22+2._dp*b_12*b_12)


       a4(4)=beta(5)*a_23+beta(1)*(a_11*a_23+2._dp*a_12*a_13)+beta(2)*b_23+&
         beta(6)*(a_11*b_23+a_23*b_11+2._dp*(a_12*b_13+a_13*b_12))+beta(3)*&
         (b_11*b_23+2._dp*b_12*b_13)
         a4(5)=beta(5)*a_13+beta(1)*(a_22*a_13+2._dp*a_12*a_23)+beta(2)*b_13+&
         beta(6)*(a_22*b_13+a_13*b_22+2._dp*(a_12*b_23+a_23*b_12))+beta(3)*&
         (b_22*b_13+2._dp*b_12*b_23)


       a4(6)=3._dp*beta(5)*a_13+3._dp*beta(1)*a_11*a_13+3._dp*beta(2)*b_13+&
         3._dp*beta(6)*(a_11*b_13+a_13*b_11)+3._dp*beta(3)*b_11*b_13
         a4(7)=3._dp*beta(5)*a_12+3._dp*beta(1)*a_11*a_12+3._dp*beta(2)*b_12+&
         3._dp*beta(6)*(a_11*b_12+a_12*b_11)+3._dp*beta(3)*b_11*b_12
         a4(8)=3._dp*beta(5)*a_23+3._dp*beta(1)*a_22*a_23+3._dp*beta(2)*b_23+&
         3._dp*beta(6)*(a_22*b_23+a_23*b_22)+3._dp*beta(3)*b_22*b_23
         a4(9)=3._dp*beta(5)*a_12+3._dp*beta(1)*a_22*a_12+3._dp*beta(2)*b_12+&
         3._dp*beta(6)*(a_22*b_12+a_12*b_22)+3._dp*beta(3)*b_22*b_12

!------------------------------------------------------------------------------
    END SUBROUTINE IBOF
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 END  SUBROUTINE CaffeSolver
!------------------------------------------------------------------------------


FUNCTION a11Flux(Model, Node, depth) RESULT(flux)
  USE Types
  USE DefUtils
  USE CoordinateSystems
  USE SolverUtils
  USe ElementDescription
 !-------------------------------------------------------------------
 IMPLICIT NONE
 !-------------------------external variables-----------------------
 TYPE(Model_t) :: Model
 TYPE(Element_t), POINTER :: BoundaryElement, ParentElement
 TYPE(ValueList_t),POINTER :: Material
 INTEGER :: Node, N, dim, istat, i, j, NBoundary, BoundaryElementNode, other_body_id, & 
          body_id, material_id
 REAL(KIND=dp) :: depth, flux, gradient, diffus
 LOGICAL :: AllocationsDone, GotIt

 !----------------! get some information upon active boundary element and its parent  !----------------------------------------

 BoundaryElement => Model % CurrentElement

!-----------------------Do nothing if this is a halo-element of a parallel mesh?  !------------------------------------------

 IF (ParEnv % myPe .NE. BoundaryElement % partIndex) RETURN


 IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN     
 CALL FATAL('a11Flux','No boundary element found')  
 END IF

 NBoundary = BoundaryElement % Type % NumberOfNodes  
 DO BoundaryElementNode=1,NBoundary     
   IF (Node .EQ. BoundaryElement % NodeIndexes(BoundaryElementNode)) THEN        
      GotIt = .TRUE.        
      EXIT     
   END IF  
 END DO

 IF (.NOT.GotIt) THEN
   CALL WARN('a11Flux','Node not found in Current Element')     
   flux = 0.0D00
   RETURN
 END IF

 other_body_id = BoundaryElement % BoundaryInfo % outbody
  IF (other_body_id < 1) THEN ! only one body in calculation
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
  ELSE ! we are dealing with a body-body boundary and assume that the normal is pointing outwards
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
  END IF

! just to be on the save side, check again
!-----------------------------------------
  IF ( .NOT. ASSOCIATED(ParentElement) ) THEN
     WRITE(Message,'(A,I10,A)')&
          'Parent Element for Boundary element no. ',&
          BoundaryElement % ElementIndex, ' not found'
     CALL FATAL('a11Flux',Message)
  END IF

  body_id = ParentElement % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  Material => Model % Materials(material_id) % Values
  IF ((.NOT. ASSOCIATED(Material)) .OR. (.NOT. GotIt)) THEN
     WRITE(Message,'(A,I10,A,I10)')&
          'No material values found for body no ', body_id,&
          ' under material id ', material_id
     CALL FATAL('a11Flux',Message)
  END IF

 diffus = GetConstReal(Material, 'Fabric flux diffusion', GotIt)
IF (.NOT. GotIt) THEN
    CALL FATAL('a11Flux:','No flux diffusion found')
 END IF

 gradient = ( (1.0/3.0) - 0.05 ) / depth
 flux = -diffus * gradient
 
END FUNCTION a11Flux

FUNCTION a22Flux(Model, Node, depth) RESULT(flux)
  USE Types
  USE DefUtils
  USE CoordinateSystems
  USE SolverUtils
  USe ElementDescription
 !-------------------------------------------------------------------
 IMPLICIT NONE
 !-------------------------external variables-----------------------
 TYPE(Model_t) :: Model
 TYPE(Element_t), POINTER :: BoundaryElement, ParentElement
 TYPE(ValueList_t),POINTER :: Material
 INTEGER :: Node, N, dim, istat, i, j, NBoundary, BoundaryElementNode, other_body_id, & 
          body_id, material_id
 REAL(KIND=dp) :: depth, flux, gradient, diffus
 LOGICAL :: AllocationsDone, GotIt

 !----------------! get some information upon active boundary element and its parent  !----------------------------------------

 BoundaryElement => Model % CurrentElement

!-----------------------Do nothing if this is a halo-element of a parallel mesh?  !------------------------------------------

 IF (ParEnv % myPe .NE. BoundaryElement % partIndex) RETURN


 IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN     
 CALL FATAL('a22Flux','No boundary element found')  
 END IF

 NBoundary = BoundaryElement % Type % NumberOfNodes  
 DO BoundaryElementNode=1,NBoundary     
   IF (Node .EQ. BoundaryElement % NodeIndexes(BoundaryElementNode)) THEN        
      GotIt = .TRUE.        
      EXIT     
   END IF  
 END DO

 IF (.NOT.GotIt) THEN
   CALL WARN('a22Flux','Node not found in Current Element')     
   flux = 0.0D00
   RETURN
 END IF

 other_body_id = BoundaryElement % BoundaryInfo % outbody
  IF (other_body_id < 1) THEN ! only one body in calculation
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
  ELSE ! we are dealing with a body-body boundary and assume that the normal is pointing outwards
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
  END IF

! just to be on the save side, check again
!-----------------------------------------
  IF ( .NOT. ASSOCIATED(ParentElement) ) THEN
     WRITE(Message,'(A,I10,A)')&
          'Parent Element for Boundary element no. ',&
          BoundaryElement % ElementIndex, ' not found'
     CALL FATAL('a22Flux',Message)
  END IF

  body_id = ParentElement % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  Material => Model % Materials(material_id) % Values
  IF ((.NOT. ASSOCIATED(Material)) .OR. (.NOT. GotIt)) THEN
     WRITE(Message,'(A,I10,A,I10)')&
          'No material values found for body no ', body_id,&
          ' under material id ', material_id
     CALL FATAL('a22Flux',Message)
  END IF

 diffus = GetConstReal(Material, 'Fabric flux diffusion', GotIt)
IF (.NOT. GotIt) THEN
    CALL FATAL('a22Flux:','No flux diffusion found')
 END IF

 gradient = ( (1.0/3.0) - 0.05 ) / depth
 flux = -diffus * gradient
 
END FUNCTION a22Flux




























