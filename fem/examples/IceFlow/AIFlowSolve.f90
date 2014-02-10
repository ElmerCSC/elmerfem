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
! * Module containing a solver for (primarily thermal) anisotropic flow
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
! *                Modified by:  Fabien / OG 
! *
! *       Date of modification: 13/10/05 from version 1.5
! *                             of  AIFlowSolver_ai.
! *
! *****************************************************************************/
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE AIFlowSolver( Model,Solver,dt,TransientSimulation )
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

     INTEGER :: i, j, k, l, n, t, iter, NDeg, STDOFs, LocalNodes, istat
     INTEGER :: dim, comp 

     TYPE(ValueList_t),POINTER :: Material, BC, BodyForce
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement

     REAL(KIND=dp) :: RelativeChange, UNorm, PrevUNorm, Gravity(3), &
         Normal(3), NewtonTol, NonlinearTol, s, Wn(7), MinSRInvariant
         

     REAL(KIND=dp)  :: NodalStresses(3,3), &
       NodalStrainRate(3,3),  NodalSpin(3,3)   

     REAL(KIND=dp), ALLOCATABLE :: Basis(:),ddBasisddx(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: dBasisdx(:,:), SlipCoeff(:,:)
     REAL(KIND=dp) :: u,v,w,detJ
     
     LOGICAL :: stat, CSymmetry 
       
     INTEGER, PARAMETER :: INDi(1:6) = (/ 1, 2, 3, 1, 2, 3 /) ,&
           INDj(1:6)=(/ 1, 2, 3, 2, 3, 1 /)

     INTEGER :: NewtonIter, NonlinearIter

     TYPE(Variable_t), POINTER :: AIFlowSol, TempSol, Var, FabricVariable
     TYPE(Variable_t), POINTER :: SpinVar 
     REAL(KIND=dp), POINTER :: SpinValues(:)
     INTEGER, POINTER :: SpinPerm(:)

     TYPE(Variable_t), POINTER :: DevStressVar 
     REAL(KIND=dp), POINTER :: DSValues(:)
     INTEGER, POINTER :: DSPerm(:)

     TYPE(Variable_t), POINTER :: StrainRateVar 
     REAL(KIND=dp), POINTER :: SRValues(:)
     INTEGER, POINTER :: SRPerm(:)

     REAL(KIND=dp), POINTER :: Temperature(:),AIFlow(:),Work(:,:), &
           ForceVector(:), VonMises(:), NodalAIFlow(:), AIFlowComp(:), &
           FabricValues(:)

     CHARACTER(LEN=MAX_NAME_LEN) :: EquationName

     INTEGER, POINTER :: TempPerm(:), AIFlowPerm(:), NodeIndexes(:), &
                         FabricPerm(:)

     INTEGER :: AIFlowType
     LOGICAL :: GotForceBC, GotIt, NewtonLinearization = .FALSE., &
                NormalTangential=.FALSE.

     INTEGER :: body_id,bf_id
     INTEGER :: old_body = -1
     LOGICAL :: Isotropic, AllocationsDone = .FALSE., FreeSurface, &
                Requal0

     REAL(KIND=dp) :: FabricGrid(4878)

     INTEGER, DIMENSION(3:8) :: NnodeLinear=(/3, 4, 4, 5, 6, 8/)
           
     REAL(KIND=dp), ALLOCATABLE:: LocalMassMatrix(:,:), &
       LocalStiffMatrix(:,:), LoadVector(:,:), LocalForce(:), &
       LocalTemperature(:), Alpha(:,:), Beta(:), & 
       ReferenceTemperature(:), BoundaryDispl(:), K1(:), K2(:), E1(:), &
       E2(:), E3(:), TimeForce(:), RefS(:), RefD(:), RefSpin(:), &
       LocalVelo(:,:), LocalFluidity(:)
            
     INTEGER :: NumberOfBoundaryNodes
     INTEGER, POINTER :: BoundaryReorder(:)

     REAL(KIND=dp) :: Bu, Bv, Bw, RM(3,3)
     REAL(KIND=dp), POINTER :: BoundaryNormals(:,:), &
         BoundaryTangent1(:,:), BoundaryTangent2(:,:)
     CHARACTER(LEN=MAX_NAME_LEN) :: viscosityFile

     REAL(KIND=dp) :: at, at0, CPUTime, RealTime


!------------------------------------------------------------------------------
     SAVE NumberOfBoundaryNodes,BoundaryReorder,BoundaryNormals, &
              BoundaryTangent1, BoundaryTangent2, FabricGrid, viscosityFile

     SAVE TimeForce, Basis, dBasisdx, ddBasisddx
     SAVE LocalMassMatrix, LocalStiffMatrix, LoadVector, &
       LocalForce, ElementNodes, Alpha, Beta, LocalTemperature, &
       Isotropic,AllocationsDone,ReferenceTemperature,BoundaryDispl, &
       NodalAIFlow, K1, K2, E1, E2, E3, Wn, MinSRInvariant, old_body, &
       LocalFluidity

     SAVE RefD, RefS, RefSpin, LocalVelo, SlipCoeff 
              
!------------------------------------------------------------------------------
!  Read constants from constants section of SIF file
!------------------------------------------------------------------------------
      Wn(7) = GetConstReal( Model % Constants, 'Gas Constant', GotIt )
      IF (.NOT.GotIt) THEN
        WRITE(Message,'(A)') 'VariableGas Constant  not found. &
                     &Setting to 8.314'
        CALL INFO('AIFlowSolve', Message, level=20)
        Wn(7) = 8.314
      ELSE
        WRITE(Message,'(A,F10.4)') 'Gas Constant = ',   Wn(7)
        CALL INFO('AIFlowSolve', Message , level = 20)
      END IF
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
      IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

      AIFlowSol => Solver % Variable
      AIFlowPerm => AIFlowSol % Perm
      STDOFs =  AIFlowSol % DOFs
      AIFlow => AIFlowSol % Values

      LocalNodes = COUNT( AIFlowPerm > 0 )
      IF ( LocalNodes <= 0 ) RETURN

      TempSol => VariableGet( Solver % Mesh % Variables, 'Temperature' )
      IF ( ASSOCIATED( TempSol) ) THEN
        TempPerm    => TempSol % Perm
        Temperature => TempSol % Values
      END IF

      FabricVariable => VariableGet(Solver % Mesh % Variables, 'Fabric')
      IF ( ASSOCIATED( FabricVariable ) ) THEN
       FabricPerm    => FabricVariable % Perm
       FabricValues => FabricVariable % Values
      END IF

      SpinVar => VariableGet(Solver % Mesh %Variables,'Spin')
      IF ( ASSOCIATED( SpinVar ) ) THEN
      SpinPerm => SpinVar % Perm    
      SpinValues => SpinVar % Values  
      END IF
      
      StrainRateVar => VariableGet(Solver % Mesh %Variables,'StrainRate')
      IF ( ASSOCIATED( StrainRateVar ) ) THEN
      SRPerm => StrainRateVar % Perm    
      SRValues => StrainRateVar % Values  
      END IF

      DevStressVar => &
               VariableGet(Solver % Mesh % Variables,'DeviatoricStress')
      IF ( ASSOCIATED( DevStressVar ) ) THEN
      DSPerm => DevStressVar % Perm    
      DSValues => DevStressVar % Values  
      END IF

      StiffMatrix => Solver % Matrix
      ForceVector => StiffMatrix % RHS
      UNorm = Solver % Variable % Norm

!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
      IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
        N = Model % MaxElementNodes
        dim = CoordinateSystemDimension()

       IF ( AllocationsDone ) THEN
         DEALLOCATE( ElementNodes % x,     &
                     ElementNodes % y,     &
                     ElementNodes % z,     &
                     BoundaryDispl,        &
                     ReferenceTemperature, &
                     LocalTemperature,     &
                     LocalVelo,            &
                     LocalForce,           &
                     RefD, RefS, RefSpin,  &
                     LocalMassMatrix,      &
                     LocalStiffMatrix,     &
                     LoadVector, Alpha, Beta, &
                     SlipCoeff, LocalFluidity )
       END IF

       ALLOCATE( ElementNodes % x( N ), &
                 ElementNodes % y( N ), &
                 ElementNodes % z( N ), &
                 BoundaryDispl( N ), &
                 ReferenceTemperature( N ), &
                 LocalTemperature( N ), &
                 K1( N ), K2( N ), E1( N ), E2( N ), E3( N ), &
                 LocalForce( 2*STDOFs*N ),&
                 RefS(2*dim*LocalNodes ),&                              
                 RefD(2*dim*LocalNodes ),&                              
                 RefSpin((2*dim-3)*LocalNodes ),&                       
                 LocalVelo( 3,N ),&                                     
                 Basis( 2*N ),ddBasisddx(1,1,1), dBasisdx( 2*N,3 ), &
                 TimeForce( 2*STDOFs*N ), &
                 LocalMassMatrix( 2*STDOFs*N,2*STDOFs*N ),  &
                 LocalStiffMatrix( 2*STDOFs*N,2*STDOFs*N ),  &
                 LoadVector( 4,N ), Alpha( 3,N ), Beta( N ), &
                 SlipCoeff(3,N), LocalFluidity(N), STAT=istat )

       IF ( istat /= 0 ) THEN
          CALL Fatal( 'AIFlowSolve', 'Memory allocation error.' )
       END IF
!------------------------------------------------------------------------------
!    Check for normal/tangetial coordinate system defined velocities
!------------------------------------------------------------------------------
       CALL CheckNormalTangentialBoundary( Model, &
        'Normal-Tangential AIFlow',NumberOfBoundaryNodes, &
          BoundaryReorder, BoundaryNormals, BoundaryTangent1, &
             BoundaryTangent2, Model % Dimension )
!------------------------------------------------------------------------------

       AllocationsDone = .TRUE.
      END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      NonlinearTol = GetConstReal( Solver % Values, &
        'Nonlinear System Convergence Tolerance' )

      NewtonTol = GetConstReal( Solver % Values, &
        'Nonlinear System Newton After Tolerance' )

      NewtonIter = GetInteger( Solver % Values, &
        'Nonlinear System Newton After Iterations' )

      NonlinearIter = GetInteger( Solver % Values, &
         'Nonlinear System Max Iterations',GotIt )

      IF ( .NOT.GotIt ) NonlinearIter = 1
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

      EquationName = GetString( Solver % Values, 'Equation' )

      FreeSurface = .FALSE.

!------------------------------------------------------------------------------
      DO iter=1,NonlinearIter

       at  = CPUTime()
       at0 = RealTime()

       CALL Info( 'AIFlowSolve', ' ', Level=4 )
       CALL Info( 'AIFlowSolve', ' ', Level=4 )
       CALL Info( 'AIFlowSolve', &
                   '-------------------------------------',Level=4 )
       WRITE( Message, * ) 'ANISOTROPIC FLOW SOLVER ITERATION', iter
       CALL Info( 'AIFlowSolve', Message,Level=4 )
       CALL Info( 'AIFlowSolve', &
                   '-------------------------------------',Level=4 )
       CALL Info( 'AIFlowSolve', ' ', Level=4 )
       CALL Info( 'AIFlowSolve', 'Starting assembly...',Level=4 )
!------------------------------------------------------------------------------
!      Compute average normals for boundaries having the normal & tangential
!      field components specified on the boundaries
!------------------------------------------------------------------------------
       IF ( (iter == 1 .OR. FreeSurface) .AND. NumberOfBoundaryNodes > 0 ) THEN
          CALL AverageBoundaryNormals( Model, &
               'Normal-Tangential AIFlow', NumberOfBoundaryNodes, &
            BoundaryReorder, BoundaryNormals, BoundaryTangent1, &
               BoundaryTangent2, Model % Dimension )
       END IF
!------------------------------------------------------------------------------
       CALL DefaultInitialize()
!------------------------------------------------------------------------------
       DO t=1,Solver % NumberOFActiveElements

         IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ',  &
             INT(100.0 - 100.0 * (Solver % NumberOfActiveElements-t) / &
             (1.0*Solver % NumberOfActiveElements)), ' % done'
           CALL Info( 'AIFlowSolve', Message, Level=5 )
           at0 = RealTime()
         END IF

         CurrentElement => GetActiveElement(t)
         n = GetElementNOFNodes()
         NodeIndexes => CurrentElement % NodeIndexes

         ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes(1:n))
         ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes(1:n))
         ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes(1:n))

         Material => GetMaterial()

!------------------------------------------------------------------------------
!    Read in material constants from Material section
!------------------------------------------------------------------------------
         body_id = CurrentElement % BodyId
         IF ((body_id /= old_body) .OR. (t==1)) Then 
            old_body = body_id
            Call  GetMaterialDefs()
        END IF

        LocalFluidity(1:n) = ListGetReal( Material, &
                         'Fluidity Parameter', n, NodeIndexes, GotIt )
        IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Fluidity Parameter not found. &
                            &Setting to 1.0'
         CALL INFO('AIFlowSolve', Message, Level = 20)
         LocalFluidity(1:n) = 1.0
        END IF

!------------------------------------------------------------------------------
!        Set body forces
!------------------------------------------------------------------------------
         LoadVector = 0.0D0

         BodyForce => GetBodyForce()
         IF ( ASSOCIATED( BodyForce ) ) THEN
           LoadVector(1,1:n) = LoadVector(1,1:n) + ListGetReal( &
                   BodyForce, 'AIFlow Force 1', n, NodeIndexes, gotIt)
           LoadVector(2,1:n) = LoadVector(2,1:n) + ListGetReal( &
                   BodyForce, 'AIFlow Force 2', n, NodeIndexes, gotIt)
           LoadVector(3,1:n) = LoadVector(3,1:n) + ListGetReal( & 
                   BodyForce, 'AIFlow Force 3', n, NodeIndexes, gotIt)
         END IF
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

! fabric not needed if isotropic
         IF(.NOT.Isotropic) THEN
           K1(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 1 )
           K2(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 2 )
           E1(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 3 )
           E2(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 4 )
           E3(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 5 )
         ENDIF

         LocalVelo = 0.0d0
         DO i=1,STDOFs - 1
            LocalVelo(i,1:n) = AIFlow( STDOFs*(AIFlowPerm(NodeIndexes(1:n))-1) + i)
         END DO

         CALL LocalMatrix( LocalMassMatrix, LocalStiffMatrix, &
              LocalForce, LoadVector, K1, K2, E1, E2, E3, LocalVelo, &
              LocalTemperature, LocalFluidity, CurrentElement, n, &
              ElementNodes, Wn, MinSRInvariant, Isotropic)

        TimeForce = 0.0d0
        IF ( n<=NNodeLinear(GetElementFamily()) ) &
           CALL NSCondensate(N, N,STDOFs-1,LocalStiffMatrix,LocalForce,TimeForce )
!------------------------------------------------------------------------------
!        If boundary fields have been defined in normal/tangential
!        coordinate systems, we´ll have to rotate the matrix & force vector
!        to that coordinate system
!------------------------------------------------------------------------------
         IF ( NumberOfBoundaryNodes > 0 ) THEN
           CALL RotateMatrix( LocalStiffMatrix,LocalForce,n, &
           CoordinateSystemDimension(), STDOFs, &
            BoundaryReorder(NodeIndexes),BoundaryNormals,BoundaryTangent1, &
                              BoundaryTangent2 )
         END IF

!------------------------------------------------------------------------------
!        Update global matrices from local matrices 
!------------------------------------------------------------------------------
         CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )
      END DO

      CALL Info( 'AIFlowSolve', 'Assembly done', Level=4 )

!------------------------------------------------------------------------------
!     Neumann & Newton boundary conditions
!------------------------------------------------------------------------------
      DO t = 1, Model % NumberOFBoundaryElements

        CurrentElement => GetBoundaryElement(t)
        IF ( GetElementFamily() == 101 ) CYCLE
        IF ( .NOT. ActiveBoundaryElement() ) CYCLE

        n = GetElementNOFNodes()
        NodeIndexes => CurrentElement % NodeIndexes

        ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes(1:n))
        ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes(1:n))
        ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes(1:n))

        BC => GetBC()
        IF ( ASSOCIATED( BC ) ) THEN
            LoadVector = 0.0D0
            Alpha      = 0.0D0
            Beta       = 0.0D0
!------------------------------------------------------------------------------
!           Force in given direction BC: \tau\cdot n = F
!------------------------------------------------------------------------------
            GotForceBC = .FALSE.

            LoadVector(1,1:n) = &
                    ListGetReal( BC, 'Force 1', n, NodeIndexes,  GotIt )
            GotForceBC = GotForceBC .OR. gotIt

            LoadVector(2,1:n) = & 
                     ListGetReal( BC, 'Force 2', n, NodeIndexes, GotIt )
            GotForceBC = GotForceBC .OR. gotIt

            LoadVector(3,1:n) = &
                     ListGetReal( BC, 'Force 3', n, NodeIndexes, GotIt )
            GotForceBC = GotForceBC .OR. gotIt

            Beta(1:n) = &
                ListGetReal( BC, 'Normal Force', n, NodeIndexes, GotIt )
            GotForceBC = GotForceBC .OR. gotIt

!------------------------------------------------------------------------------
!             slip boundary condition BC: \tau\cdot n = R_k u_k
!------------------------------------------------------------------------------

              SlipCoeff = 0.0d0
              SlipCoeff(1,1:n) =  ListGetReal( BC, &
                    'AIFlow Slip Coeff 1', n, NodeIndexes, GotIt )
              GotForceBC = GotForceBC .OR. gotIt

              SlipCoeff(2,1:n) =  ListGetReal( BC, &
                    'AIFlow Slip Coeff 2', n, NodeIndexes, GotIt )
              GotForceBC = GotForceBC .OR. gotIt

              SlipCoeff(3,1:n) =  ListGetReal( BC, &
                    'AIFlow Slip Coeff 3', n, NodeIndexes, GotIt )
              GotForceBC = GotForceBC .OR. gotIt

              NormalTangential = ListGetLogical( BC, &
                     'Normal-Tangential AIFlow', GotIt )
               
            IF ( .NOT.GotForceBC ) CYCLE
!------------------------------------------------------------------------------
            CALL LocalMatrixBoundary( LocalStiffMatrix, LocalForce, &
                 LoadVector, Alpha, Beta, SlipCoeff, NormalTangential, &
                 CurrentElement, n, ElementNodes )
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!           If boundary fields have been defined in normal/tangetial coordinate
!           systems, we´ll have to rotate the matrix & force vector to that
!           coordinate system
!------------------------------------------------------------------------------
            IF ( NumberOfBoundaryNodes > 0 ) THEN
              CALL RotateMatrix( LocalStiffMatrix, LocalForce, n,&
                   CoordinateSystemDimension(), STDOFs, &
                   BoundaryReorder(NodeIndexes), BoundaryNormals, &
                   BoundaryTangent1, BoundaryTangent2 )
            END IF
!------------------------------------------------------------------------------
!           Update global matrices from local matrices (will also affect
!           LocalStiffMatrix and LocalForce if transientsimulation is on).
!------------------------------------------------------------------------------
            CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )
!------------------------------------------------------------------------------
         END IF
      END DO
!------------------------------------------------------------------------------

      CALL DefaultFinishAssembly()

!------------------------------------------------------------------------------
!     Dirichlet boundary conditions
!------------------------------------------------------------------------------
      CALL DefaultDirichletBCs()
!------------------------------------------------------------------------------

      CALL Info( 'AIFlowSolve', 'Set boundaries done', Level=4 )

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      PrevUNorm = UNorm

      UNorm = DefaultSolve()

      IF ( PrevUNorm + UNorm /= 0.0d0 ) THEN
         RelativeChange = 2.0d0 * ABS( PrevUNorm - UNorm) / ( PrevUnorm + UNorm)
      ELSE
         RelativeChange = 0.0d0
      END IF

      WRITE( Message, * ) 'Result Norm   : ',UNorm, PrevUNorm
      CALL Info( 'AIFlowSolve', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'AIFlowSolve', Message, Level=4 )
!------------------------------------------------------------------------------
!     If boundary fields have been defined in normal/tangential coordinate
!     systems, we´ll have to rotate the solution back to coordinate axis
!     directions
!------------------------------------------------------------------------------
      IF ( NumberOfBoundaryNodes > 0 ) THEN
        DO i=1,Model % NumberOfNodes
          k = BoundaryReorder(i)

          IF ( k > 0 ) THEN
            j = AIFlowPerm(i)

            IF ( j > 0 ) THEN
              IF ( CoordinateSystemDimension() < 3 ) THEN
                Bu = AIFlow( STDOFs*(j-1)+1 )
                Bv = AIFlow( STDOFs*(j-1)+2 )

                AIFlow( STDOFs*(j-1)+1) = BoundaryNormals(k,1) * Bu - &
                                BoundaryNormals(k,2) * Bv

                AIFlow( STDOFs*(j-1)+2) = BoundaryNormals(k,2) * Bu + &
                                BoundaryNormals(k,1) * Bv
              ELSE
                Bu = AIFlow( STDOFs*(j-1)+1 )
                Bv = AIFlow( STDOFs*(j-1)+2 )
                Bw = AIFlow( STDOFs*(j-1)+3 )

                RM(1,:) = BoundaryNormals(k,:)
                RM(2,:) = BoundaryTangent1(k,:)
                RM(3,:) = BoundaryTangent2(k,:)

                AIFlow(STDOFs*(j-1)+1) = RM(1,1)*Bu+RM(2,1)*Bv+RM(3,1)*Bw
                AIFlow(STDOFs*(j-1)+2) = RM(1,2)*Bu+RM(2,2)*Bv+RM(3,2)*Bw
                AIFlow(STDOFs*(j-1)+3) = RM(1,3)*Bu+RM(2,3)*Bv+RM(3,3)*Bw
              END IF
            END IF
          END IF
        END DO 
      END IF
!------------------------------------------------------------------------------
      IF ( RelativeChange < NewtonTol .OR. &
             iter > NewtonIter ) NewtonLinearization = .TRUE.

      IF ( RelativeChange < NonLinearTol ) EXIT

!------------------------------------------------------------------------------
    END DO ! of nonlinear iter
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!   Compute the StrainRate, Spin  and deviatoric Stress
!   Nodal values      
!------------------------------------------------------------------------------

     IF ((ASSOCIATED( StrainRateVar)).OR.(ASSOCIATED(DevStressVar))&
      .OR.(ASSOCIATED(SpinVar))) THEN
       RefD=0.
       RefS=0.
       RefSpin=0.
       IF (ASSOCIATED(StrainRateVar)) SRValues = 0.
       IF (ASSOCIATED(devStressVar)) DSValues = 0.
       IF (ASSOCIATED(SPinVar)) SpinValues = 0.

      DO t=1,Solver % NumberOFActiveElements

         CurrentElement => GetActiveElement(t)
         n = GetElementNOFNodes()
         NodeIndexes => CurrentElement % NodeIndexes

         body_id = CurrentElement % BodyId
         CSymmetry = CurrentCoordinateSystem() == AxisSymmetric
         dim = CoordinateSystemDimension()

!------------------------------------------------------------------------------
!    Read in material constants from Material section
!------------------------------------------------------------------------------
        IF (body_id /= old_body) Then 
              old_body = body_id
              Call  GetMaterialDefs()
        END IF

        LocalFluidity(1:n) = ListGetReal( Material, &
                      'Fluidity Parameter', n, NodeIndexes, GotIt )
        IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Fluidity Parameter not found. &
                            &Setting to 1.0'
         CALL INFO('AIFlowSolve', Message, Level = 20)
         LocalFluidity(1:n) = 1.0
        END IF


         ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes(1:n))
         ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes(1:n))
         ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes(1:n))

! n nodale values of the temperature
         
         LocalTemperature = 0.0D0
         IF ( ASSOCIATED(TempSol) ) THEN
           DO i=1,n
             k = TempPerm(NodeIndexes(i))
             LocalTemperature(i) = Temperature(k)
           END DO
         ELSE
           LocalTemperature(1:n) = 0.0d0
         END IF

! n nodales values of the 5 fabric parameters, not needed if isotropic
         IF(.NOT.Isotropic) Then
           K1(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 1 ) 
           K2(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 2 )
           E1(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 3 )
           E2(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 4 )
           E3(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 5 )
         END IF

! 2D U,V,p    STDOFs=3
! 3D U,V,W,p  STDOFs=4
         LocalVelo = 0.0d0
         DO i=1,STDOFs - 1
            LocalVelo(i,1:n) = AIFlow( STDOFs*(AIFlowPerm(NodeIndexes(1:n))-1) + i)
         END DO
! Go for all nodes of the element        
         Do i=1,n

! u, v, w local coord of node i
           u = CurrentElement % Type % NodeU(i)
           v = CurrentElement % Type % NodeV(i)
           w = CurrentElement % Type % NodeW(i)
       
            stat = ElementInfo(CurrentElement,ELementNodes,u,v,w,detJ, &
               Basis,dBasisdx,ddBasisddx,.FALSE.,.FALSE.)
! Axi symmetric case when R=0 strain, stress not calculated exactly in
! x=0 (I agree it is not very nice, better solution ???)
        Requal0 = .False.
        IF (( CSymmetry) .And. & 
          (SUM(ElementNodes % x(1:n) * Basis(1:n)) == 0.0)) THEN  
           Requal0 = .True.
            u= u + 0.0001  
            stat = ElementInfo(CurrentElement,ELementNodes,u,v,w,detJ, &
               Basis,dBasisdx,ddBasisddx,.FALSE.,.FALSE.)
        END IF

           CALL LocalSD(NodalStresses, NodalStrainRate, NodalSpin, & 
                 LocalVelo, LocalTemperature, LocalFluidity,  &
                K1, K2, E1, E2, E3, CSymmetry, Basis, dBasisdx, &
                CurrentElement, n, ElementNodes, dim, Wn, &
                MinSRInvariant, Isotropic)
                

        IF (Requal0) NodalSpin = 0. 

           IF (ASSOCIATED(StrainRateVar)) &
             RefD(2*dim*(SRPerm(NodeIndexes(i))-1)+1 : &
                                      2*dim*SRPerm(NodeIndexes(i))) &
             =RefD(2*dim*(SRPerm(NodeIndexes(i))-1)+1 : &
                                      2*dim*SRPerm(NodeIndexes(i))) + 1.

          IF (ASSOCIATED(DevStressVar)) &
            RefS(2*dim*(DSPerm(NodeIndexes(i))-1)+1 : &
                                      2*dim*DSPerm(NodeIndexes(i))) &
            =RefS(2*dim*(DSPerm(NodeIndexes(i))-1)+1 :  &
                                      2*dim*DSPerm(NodeIndexes(i))) + 1.

          IF (ASSOCIATED(SpinVar)) &
            RefSpin((2*dim-3)*(SpinPerm(NodeIndexes(i))-1)+1 :  &
                                (2*dim-3)*SpinPerm(NodeIndexes(i))) &
            =RefSpin((2*dim-3)*(SpinPerm(NodeIndexes(i))-1)+1 :  &
                                (2*dim-3)*SpinPerm(NodeIndexes(i))) + 1.


           IF (ASSOCIATED(StrainRateVar)) THEN
             comp=0
             DO j=1,2*dim
               comp=comp+1
               SRValues(2*dim*(SRPerm(NodeIndexes(i))-1)+comp)=&
               SRValues(2*dim*(SRPerm(NodeIndexes(i))-1)+comp) + &
                NodalStrainRate(INDi(j),INDj(j))
             END DO
           END IF

           IF (ASSOCIATED(DevStressVar)) THEN
             comp=0
             DO j=1,2*dim
               comp=comp+1
               DSValues(2*dim*(DSPerm(NodeIndexes(i))-1)+comp)=&
                DSValues(2*dim*(DSPerm(NodeIndexes(i))-1)+comp) + &
                NodalStresses(INDi(j),INDj(j))
             END DO
           END IF

           IF (ASSOCIATED(SpinVar)) THEN
             comp=0
             DO j=1,(2*dim-3)
             comp=comp+1
             SpinValues((2*dim-3)*(SpinPerm(NodeIndexes(i))-1)+comp)=&
             SPinValues((2*dim-3)*(SpinPerm(NodeIndexes(i))-1)+comp) + &
             NodalSpin(INDi(j+3),INDj(j+3))
             END DO
           END IF

          END DO
        END DO

        IF (ASSOCIATED(StrainRateVar)) THEN
              WHERE(RefD > 0.)
                SRVAlues = SRValues / RefD
              END WHERE
        END IF
        
        IF (ASSOCIATED(DevStressVar)) THEN
           WHERE(RefS > 0.)
                DSVAlues = DSValues / RefS
           END WHERE
        END IF

        IF (ASSOCIATED(SpinVar)) THEN
            WHERE(RefSpin > 0.)
                SpinVAlues = SpinValues / RefSpin
            END WHERE
        END IF
        
    END IF

!------------------------------------------------------------------------------
!  END  Compute the StrainRate and Deviatoric Stress
!------------------------------------------------------------------------------
      
CONTAINS

      SUBROUTINE GetMaterialDefs()
      ! check if we are isotropic or not
      Isotropic = ListGetLogical( Material , 'Isotropic',Gotit )
      IF (.NOT.Gotit) Then
          Isotropic = .False.
           WRITE(Message,'(A)') 'Isotropic set to False'
           CALL INFO('AIFlowSolve', Message, Level = 20)
      ELSE
           IF ( (ASSOCIATED( FabricVariable )).AND.Isotropic ) Then
              WRITE(Message,'(A)') 'Be carefull Isotropic is true &
                           & and Fabric is defined!'
              CALL INFO('AIFlowSolve', Message, Level = 1)
           END IF
      END IF

      IF (.NOT.Isotropic) Then
        ! Get the viscosity file and store the viscosities into FabricGrid
         viscosityFile = ListGetString( Material ,'Viscosity File',GotIt )
         IF (.NOT.GotIt) THEN
            WRITE(Message,'(3A)') &
                      'Viscosity File ', viscosityFile, ' not found'
           CALL FATAL('AIFlowSolve',Message)
         ELSE
             OPEN( 1, File = viscosityFile)
             DO i=1,813
                 READ( 1, '(6(e14.8))' ) FabricGrid( 6*(i-1)+1:6*(i-1)+6 )
             END DO
             CLOSE(1)
          END IF
      ENDIF


      Wn(2) = ListGetConstReal( Material , 'Powerlaw Exponent', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable  Powerlaw Exponent not found. &
                                    & Setting to 1.0'
         CALL INFO('AIFlowSolve', Message, Level = 20)
         Wn(2) = 1.0
      ELSE
       WRITE(Message,'(A,F10.4)') 'Powerlaw Exponent = ',   Wn(2)
       CALL INFO('AIFlowSolve', Message, Level = 20)
       END IF

      Wn(3) = ListGetConstReal( Material, 'Activation Energy 1', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Activation Energy 1 not found.&
                            & Setting to 1.0'
         CALL INFO('AIFlowSolve', Message, Level = 20)
         Wn(3) = 1.0
      ELSE
         WRITE(Message,'(A,F10.4)') 'Activation Energy 1 = ',   Wn(3)
         CALL INFO('AIFlowSolve', Message, Level = 20)
      END IF

      Wn(4) = ListGetConstReal( Material, 'Activation Energy 2', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Activation Energy 2 not found. &
                               &Setting to 1.0'
         CALL INFO('AIFlowSolve', Message, Level = 20)
         Wn(4) = 1.0
      ELSE
         WRITE(Message,'(A,F10.4)') 'Activation Energy 2 = ',   Wn(4)
         CALL INFO('AIFlowSolve', Message, Level = 20)
      END IF

      Wn(5) = ListGetConstReal(Material, 'Reference Temperature', GotIt)
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Reference Temperature not found. &
                               &Setting to -10.0 (Celsius)'
         CALL INFO('AIFlowSolve', Message, Level = 20)
         Wn(5) = -10.0
      ELSE
         WRITE(Message,'(A,F10.4)') 'Reference Temperature = ',   Wn(5)
         CALL INFO('AIFlowSolve', Message, Level = 20)
      END IF

      Wn(6) = ListGetConstReal( Material, 'Limit Temperature', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Limit Temperature not found. &
                               &Setting to -10.0 (Celsius)'
         CALL INFO('AIFlowSolve', Message, Level = 20)
         Wn(6) = -10.0
      ELSE
         WRITE(Message,'(A,F10.4)') 'Limit Temperature = ',   Wn(6)
         CALL INFO('AIFlowSolve', Message, Level = 20)
      END IF

! Get the Minimum value of the Effective Strain rate 
      MinSRInvariant = 1.0D-10
      IF ( Isotropic .AND. (Wn(2) > 1.0 ) ) THEN
        MinSRInvariant =  &
             ListGetConstReal( Material, 'Min Second Invariant', GotIt )
        IF (.NOT.GotIt) THEN
          WRITE(Message,'(A)') 'Variable Min Second Invariant not &
                    &found. Setting to 1.0D-10 )'
          CALL INFO('AIFlowSolve', Message, Level = 20)
        ELSE
          WRITE(Message,'(A,E14.8)') 'Min Second Invariant = ', MinSRInvariant
          CALL INFO('AIFlowSolve', Message, Level = 20)
        END IF
      END IF

!------------------------------------------------------------------------------
      END SUBROUTINE GetMaterialDefs
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      SUBROUTINE LocalMatrix( MassMatrix, StiffMatrix, ForceVector, &
              LoadVector, NodalK1, NodalK2, NodalEuler1, NodalEuler2, &
              NodalEuler3, NodalVelo, NodalTemperature, NodalFluidity, &
              Element, n, Nodes, Wn, MinSRInvariant, Isotropic )
              
              
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: StiffMatrix(:,:), MassMatrix(:,:)
     REAL(KIND=dp) :: LoadVector(:,:), NodalVelo(:,:)
     REAL(KIND=dp) :: Wn(7), MinSRInvariant
     REAL(KIND=dp), DIMENSION(:) :: ForceVector, NodalK1, NodalK2, &
             NodalEuler1, NodalEuler2, NodalEuler3, NodalTemperature, &
             NodalFluidity
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
     LOGICAL :: Isotropic
     INTEGER :: n
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(2*n),ddBasisddx(1,1,1)
     REAL(KIND=dp) :: dBasisdx(2*n,3),detJ, pBasis(n)

     REAL(KIND=dp) :: Force(3), K1, K2, Euler1, Euler2, Euler3
     Real(kind=dp) :: Bg, BGlenT

     REAL(KIND=dp), DIMENSION(4,4) :: A,M
     REAL(KIND=dp) :: Load(3), Temperature,  C(6,6)
     REAL(KIND=dp) :: nn, ss, LGrad(3,3), SR(3,3), epsi

     INTEGER :: i, j, k, p, q, t, dim, cc, NBasis, ind(3), LinearBasis

     REAL(KIND=dp) :: s,u,v,w, Radius, B(6,3), G(3,6), eta
  
     REAL(KIND=dp) :: dDispldx(3,3), ai(3), Angle(3), a2(6)
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     INTEGER, POINTER :: EdgeMap(:,:)
     INTEGER :: N_Integ

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat, CSymmetry, P2P1

     TYPE(ElementType_t), POINTER :: SaveType, LinearType

     INTERFACE
      Subroutine R2Ro(a2,dim,ai,angle)
         USE Types
         REAL(KIND=dp),intent(in) :: a2(6)
         Integer :: dim
         REAL(KIND=dp),intent(out) :: ai(3), Angle(3)
      End Subroutine R2Ro
                 
      Subroutine OPILGGE_ai(ai,Angle,Tc,W,etaI,eta36)
          USE Types
          REAL(KIND=dp) :: ai(3), Angle(3), Tc, W(7), EtaI(:), Eta36(6,6)
        END SUBROUTINE OPILGGE_ai

      END INTERFACE

!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()
      cc=dim+1

      ForceVector = 0.0D0
      StiffMatrix = 0.0D0
      MassMatrix  = 0.0D0

!    
!    Integration stuff
!    
      P2P1 = n>NNodeLinear(GetElementFamily())
      IF ( P2P1 ) THEN
        NBasis = n

        SELECT CASE(GetElementFamily())
        CASE(3)
          LinearBasis = 3
          LinearType => GetElementType(303)
        CASE(4)
          LinearBasis = 4
          LinearType => GetElementType(404)
        CASE(5)
          LinearBasis = 4
          LinearType => GetElementType(504)
        CASE(6)
          LinearBasis = 5
          LinearType => GetElementType(605)
        CASE(7)
          LinearBasis = 6
          LinearType => GetElementType(706)
        CASE(8)
          LinearBasis = 8
          LinearType => GetElementType(808)
        END SELECT
        SaveType => Element % Type
        IntegStuff = GaussPoints( Element )
      ELSE
        NBasis = 2*n
        IntegStuff = GaussPoints( Element, Element % Type % GaussPoints2 )
      END IF

      U_Integ => IntegStuff % u
      V_Integ => IntegStuff % v
      W_Integ => IntegStuff % w
      S_Integ => IntegStuff % s
      N_Integ =  IntegStuff % n
!
!   Now we start integrating
!
      DO t=1,N_Integ

      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)

!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      IF ( P2P1 ) THEN
         Element % Type => LinearType
         stat = ElementInfo( Element,Nodes,u,v,w,detJ,PBasis )
         Element % Type => SaveType
      END IF
      stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
         Basis,dBasisdx,ddBasisddx,.FALSE.,.TRUE. )

      s = detJ * S_Integ(t)
!------------------------------------------------------------------------------
!  
!     Force at integration point
!   
      Force = 0.0d0
      DO i=1,dim
         Force(i) = SUM( LoadVector(i,1:n)*Basis(1:n))
      END DO


      Radius = SUM( Nodes % x(1:n) * Basis(1:n) )
!
!     Temperature at the integration point
!
      Temperature = SUM( NodalTemperature(1:n)*Basis(1:n) )
      Wn(1) = SUM( NodalFluidity(1:n)*Basis(1:n) )
      Bg = BGlenT(Temperature,Wn)

! if not isotropic use GOLF
      C = 0.0_dp
      IF (.NOT.Isotropic) Then
         a2(1) = SUM( NodalK1(1:n) * Basis(1:n) ) 
         a2(2) = SUM( NodalK2(1:n) * Basis(1:n) ) 
         a2(3) = 1.d0 - a2(1) - a2(2)
         a2(4) = SUM( NodalEuler1(1:n) * Basis(1:n) )
         a2(5) = SUM( NodalEuler2(1:n) * Basis(1:n) )
         a2(6) = SUM( NodalEuler3(1:n) * Basis(1:n) )
      
         CALL R2Ro(a2,dim,ai,angle)
         CALL OPILGGE_ai(ai,Angle,Temperature,Wn,FabricGrid,C)

! else use isotropic law
      ELSE
          eta = 1.0_dp / Bg
      ENDIF
      
      CSymmetry = CurrentCoordinateSystem() == AxisSymmetric
      IF ( CSymmetry ) s = s * Radius

!
! Case non-linear and ISOTROPIC
! -----------------------------

      IF ( Isotropic .AND. (Wn(2) > 1.0) ) THEN
        LGrad = MATMUL( NodalVelo(:,1:n), dBasisdx(1:n,:) )
        SR = 0.5 * ( LGrad + TRANSPOSE(LGrad) )
      
        IF ( CSymmetry ) THEN
          SR(1,3) = 0.0
          SR(2,3) = 0.0
          SR(3,1) = 0.0
          SR(3,2) = 0.0
          SR(3,3) = 0.0
          IF ( Radius > 10*AEPS ) THEN
            SR(3,3) = SUM( Nodalvelo(1,1:n) * Basis(1:n) ) /Radius
                 
          END IF
          epsi = SR(1,1)+SR(2,2)+SR(3,3)
          DO i=1,3   
            SR(i,i) = SR(i,i) - epsi/3.0
          END DO
        ELSE
          epsi = SR(1,1)+SR(2,2)+SR(3,3)
          DO i=1,dim 
            SR(i,i) = SR(i,i) - epsi/dim
          END DO
        END IF

        ss = 0.0_dp
        DO i = 1, 3
          DO j = 1, 3
            ss = ss + SR(i,j)**2
          END DO
        END DO
        nn =(1.0 - Wn(2))/Wn(2)

        ss = SQRT(2.0_dp * ss)
        IF (ss < MinSRInvariant ) ss = MinSRInvariant

        eta  = eta * (ss/Bg)**nn 

      END IF

!
!    Loop over basis functions (of both unknowns and weights)
!

      M = 0.0d0
      DO p=1,NBasis

       
        IF (.Not.Isotropic) THEN
          G = 0.0d0
          IF ( CSymmetry ) THEN
            G(1,1) = dBasisdx(p,1)
            G(1,3) = Basis(p) / Radius
            G(1,4) = dBasisdx(p,2)
            G(2,2) = dBasisdx(p,2)
            G(2,4) = dBasisdx(p,1)
          ELSE
            G(1,1) = dBasisdx(p,1)
            G(2,2) = dBasisdx(p,2)
            G(3,3) = dBasisdx(p,3)
            G(1,4) = dBasisdx(p,2)
            G(2,4) = dBasisdx(p,1)
            G(2,5) = dBasisdx(p,3)
            G(3,5) = dBasisdx(p,2)
            G(1,6) = dBasisdx(p,3)
            G(3,6) = dBasisdx(p,1)
          END IF
          G = MATMUL( G, C )
        END IF 
 
       DO q=1,NBasis

         IF (.Not.Isotropic) THEN
           B = 0.0d0
           IF ( CSymmetry ) THEN
             B(1,1) = dBasisdx(q,1)
             B(2,2) = dBasisdx(q,2)
             B(3,1) = Basis(q) / Radius
             B(4,1) = dBasisdx(q,2)
             B(4,2) = dBasisdx(q,1)
           ELSE
             B(1,1) = dBasisdx(q,1)
             B(2,2) = dBasisdx(q,2)
             B(3,3) = dBasisdx(q,3)
             B(4,1) = dBasisdx(q,2)
             B(4,2) = dBasisdx(q,1)
             B(5,2) = dBasisdx(q,3)
             B(5,3) = dBasisdx(q,2)
             B(6,1) = dBasisdx(q,3)
             B(6,3) = dBasisdx(q,1)
           END IF
           A(1:3,1:3) = MATMUL( G, B )
         ELSE

         A = 0.0d0
           DO i=1,dim
              DO j=1,dim
                A(i,i)=A(i,i) + eta* dbasisdx(q,j)*dbasisdx(p,j)
                A(i,j)=A(i,j) + eta* dbasisdx(q,i)*dbasisdx(p,j)
              END DO
           END DO

         END IF

! Pressure gradient
         DO i=1,dim
            A(i,cc) = -dBasisdx(p,i) * Basis(q)
         END DO
         IF ( CSymmetry ) A(1,cc) =  A(1,cc) - Basis(p) * Basis(q) / Radius

! Continuity equation:
         DO i=1,dim
            IF ( P2P1 .AND. p <= LinearBasis ) THEN
              A(cc,i) = dBasisdx(q,i) * PBasis(p)
            ELSE IF (.NOT.P2P1) THEN
              A(cc,i) = dBasisdx(q,i) * Basis(p)
            END IF
         END DO
         IF ( P2P1 .AND. p <= LinearBasis ) THEN
            IF ( CSymmetry ) A(cc,1) =  A(cc,1) + PBasis(p) * Basis(q) / Radius
         ELSE IF (.NOT.P2P1) THEN
            IF ( CSymmetry ) A(cc,1) =  A(cc,1) + Basis(p) * Basis(q) / Radius
         END IF
         A(cc, cc) = 0.0d0

! Add nodal matrix to element matrix
         DO i=1,cc
            DO j=1,cc
               StiffMatrix( cc*(p-1)+i,cc*(q-1)+j ) =  &
                    StiffMatrix( cc*(p-1)+i,cc*(q-1)+j ) + s*A(i,j)
            END DO
         END DO

       END DO

! The righthand side...
        Load = 0.0d0
  
        DO i=1,dim
           Load(i) = Load(i) + Force(i) * Basis(p)
        END DO

        DO i=1,dim
           ForceVector(cc*(p-1)+i) = ForceVector(cc*(p-1)+i) + s*Load(i)
        END DO

      END DO

! Non linear element -> apply p2p1 element

      j = GetElementFamily()
      IF ( n > NNodeLinear(j) ) THEN
        EdgeMap => GetEdgeMap(j)
        DO i=j+1,n
          p = EdgeMap(i-j,1)
          q = EdgeMap(i-j,2)
          StiffMatrix( cc*i, : ) = 0.0d0
          MassMatrix(  cc*i, : ) = 0.0d0
          ForceVector( cc*i ) = 0.0d0
          StiffMatrix( cc*i, cc*i ) = 1.0d0
          StiffMatrix( cc*i, cc*p ) = -1.0d0/2.0d0
          StiffMatrix( cc*i, cc*q ) = -1.0d0/2.0d0
        END DO
      END IF 


      END DO 
!------------------------------------------------------------------------------
      END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      SUBROUTINE LocalMatrixBoundary( BoundaryMatrix, BoundaryVector, &
                 LoadVector, NodalAlpha, NodalBeta, NodalSlipCoeff, & 
                  NormalTangential, Element, n, Nodes )
                      
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:)
     REAL(KIND=dp) :: NodalAlpha(:,:),NodalBeta(:),LoadVector(:,:)
     REAL(KIND=dp) :: NodalSlipCoeff(:,:)
     TYPE(Element_t),POINTER  :: Element
     TYPE(Nodes_t)    :: Nodes
     LOGICAL :: NormalTangential
     INTEGER :: n
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(n),ddBasisddx(1,1,1)
     REAL(KIND=dp) :: dBasisdx(n,3),detJ

     REAL(KIND=dp) :: u,v,w,s
     REAL(KIND=dp) :: Force(3),Alpha(3),Beta,Normal(3)
     REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

     REAL(KIND=dp) :: Tangent(3),Tangent2(3),Vect(3), SlipCoeff
     INTEGER :: i,t,q,p,dim,N_Integ, c

     LOGICAL :: stat

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
!------------------------------------------------------------------------------

      dim = CoordinateSystemDimension()
      c=dim+1

      BoundaryVector = 0.0D0
      BoundaryMatrix = 0.0D0
!
!  Integration stuff
!
      IntegStuff = GaussPoints( element )
      U_Integ => IntegStuff % u
      V_Integ => IntegStuff % v
      W_Integ => IntegStuff % w
      S_Integ => IntegStuff % s
      N_Integ =  IntegStuff % n
!
!  Now we start integrating
!
      DO t=1,N_Integ

       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)

!------------------------------------------------------------------------------
!    Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
                Basis, dBasisdx, ddBasisddx, .FALSE. )

       s = detJ * S_Integ(t)
       IF ( CurrentCoordinateSystem() == AxisSymmetric ) &
        s = s * SUM( Nodes % x(1:n) * Basis(1:n) )
!------------------------------------------------------------------------------
       Force = 0.0D0
       DO i=1,dim
         Force(i) = SUM( LoadVector(i,1:n)*Basis(1:n) )
         Alpha(i) = SUM( NodalAlpha(i,1:n)*Basis(1:n) )
       END DO

       Normal = NormalVector( Element,Nodes,u,v,.TRUE. )
       Force = Force + SUM( NodalBeta(1:n)*Basis(1:n) ) * Normal

       SELECT CASE( Element % TYPE % DIMENSION )
       CASE(1)
        Tangent(1) =  Normal(2)
        Tangent(2) = -Normal(1)
        Tangent(3) =  0.0d0
       CASE(2)
        CALL TangentDirections( Normal, Tangent, Tangent2 ) 
       END SELECT
  
       IF ( ANY( NodalSlipCoeff(:,:) /= 0.0d0 ) ) THEN
         DO p=1,n
           DO q=1,n
             DO i=1,DIM
              SlipCoeff = SUM( NodalSlipCoeff(i,1:n) * Basis(1:n) )
  
              IF (NormalTangential ) THEN
                SELECT CASE(i)
                   CASE(1)
                     Vect = Normal
                   CASE(2)
                     Vect = Tangent
                   CASE(3)
                     Vect = Tangent2
                END SELECT
  
                DO j=1,DIM
                   DO k=1,DIM
                      BoundaryMatrix( (p-1)*c+j,(q-1)*c+k ) = &
                         BoundaryMatrix( (p-1)*c+j,(q-1)*c+k ) + &
                          s * SlipCoeff * Basis(q) * Basis(p) * Vect(j) * Vect(k)
                   END DO
                END DO
               ELSE
                 BoundaryMatrix( (p-1)*c+i,(q-1)*c+i ) = &
                     BoundaryMatrix( (p-1)*c+i,(q-1)*c+i ) + &
                          s * SlipCoeff * Basis(q) * Basis(p)
               END IF
             END DO
           END DO
         END DO
       END IF



      DO p=1,N
       DO q=1,N
         DO i=1,dim
           BoundaryMatrix((p-1)*(dim+1)+i,(q-1)*(dim+1)+i) =  &
             BoundaryMatrix((p-1)*(dim+1)+i,(q-1)*(dim+1)+i) + &
               s * Alpha(i) * Basis(q) * Basis(p)
         END DO
       END DO
      END DO

      DO q=1,N
       DO i=1,dim
         BoundaryVector((q-1)*(dim+1)+i) = BoundaryVector((q-1)*(dim+1)+i) + &
                   s * Basis(q) * Force(i)
       END DO
      END DO

      END DO
!------------------------------------------------------------------------------
      END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE LocalSD( Stress, StrainRate, Spin, &
        NodalVelo, NodalTemp, NodalFluidity, NodalK1,  &
        NodalK2, NodalE1, NodalE2, NodalE3, CSymmetry, &
        Basis, dBasisdx, Element, n,  Nodes, dim,  Wn, MinSRInvariant, &
        Isotropic)
!------------------------------------------------------------------------------
!    Subroutine to computre the nodal Strain-Rate, Stress, ...
!------------------------------------------------------------------------------
     LOGICAL ::  CSymmetry 
     INTEGER :: n, dim
     INTEGER :: INDi(6),INDj(6)
     REAL(KIND=dp) :: Stress(:,:), StrainRate(:,:), Spin(:,:)
     REAL(KIND=dp) :: NodalVelo(:,:), NodalTemp(:), NodalFluidity(:)
     REAL(KIND=dp) :: Basis(2*n), ddBasisddx(1,1,1)
     REAL(KIND=dp) :: dBasisdx(2*n,3)
     REAL(KIND=dp) :: detJ
     REAL(KIND=dp) :: NodalK1(:), NodalK2(:)
     REAL(KIND=dp) :: NodalE1(:), NodalE2(:), NodalE3(:)
     REAL(KIND=dp) :: u, v, w      
     REAL(KIND=dp) :: Wn(7),  D(6), MinSRInvariant
     LOGICAL :: Isotropic
      
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
     LOGICAL :: stat
     INTEGER :: i,j,k,p,q
     REAL(KIND=dp) :: LGrad(3,3),   Radius, Temp, ai(3), Angle(3),a2(6)
     REAL(KIND=dp) :: C(6,6), epsi
     Real(kind=dp) :: Bg, BGlenT, ss, nn
!------------------------------------------------------------------------------
     INTERFACE
      Subroutine R2Ro(a2,dim,ai,angle)
         USE Types
         REAL(KIND=dp),intent(in) :: a2(6)
         Integer :: dim
         REAL(KIND=dp),intent(out) :: ai(3), Angle(3)
      End Subroutine R2Ro
                 
      Subroutine OPILGGE_ai(ai,Angle,Tc,W,etaI,eta36)
          USE Types
          REAL(KIND=dp) :: ai(3), Angle(3), Tc, W(7), EtaI(:), Eta36(6,6)
        END SUBROUTINE OPILGGE_ai
      END INTERFACE
!------------------------------------------------------------------------------
     
      Stress = 0.0
      StrainRate = 0.0
      Spin = 0.0

!
!     Temperature at the integration point
      Temp = SUM( NodalTemp(1:n)*Basis(1:n) )
      Wn(1) = SUM( NodalFluidity(1:n)*Basis(1:n) )
      
      C = 0.0_dp
      IF (.Not.Isotropic) then
!    Material parameters at that point
        a2(1) = SUM( NodalK1(1:n) * Basis(1:n) ) 
        a2(2) = SUM( NodalK2(1:n) * Basis(1:n) ) 
        a2(3) = 1.d0 - a2(1) - a2(2)
        a2(4) = SUM( NodalE1(1:n) * Basis(1:n) )
        a2(5) = SUM( NodalE2(1:n) * Basis(1:n) )
        a2(6) = SUM( NodalE3(1:n) * Basis(1:n) )
      
        CALL R2Ro(a2,dim,ai,Angle)
        CALL OPILGGE_ai(ai,Angle,Temp,Wn,FabricGrid,C)

      ELSE
         Bg=BGlenT(Temp,Wn)
         Do i=1,3
            C(i,i)=2._dp/Bg
         End do
         Do i=4,6
            C(i,i)=1.0_dp/Bg
         End do
      END IF
!
!    Compute strainRate : 
!    -------------------

      LGrad = MATMUL( NodalVelo(:,1:n), dBasisdx(1:n,:) )
        
      StrainRate = 0.5 * ( LGrad + TRANSPOSE(LGrad) )
      
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
! Case non-linear and ISOTROPIC
! -----------------------------

      IF ( Isotropic .AND. (Wn(2) > 1.0) ) THEN
        ss = 0.0_dp
        DO i = 1, 3
          DO j = 1, 3
            ss = ss + StrainRate(i,j)**2
          END DO
        END DO
        nn = (1.0 - Wn(2))/Wn(2)
        ss = SQRT(2.0_dp * ss)
         
        IF (ss < MinSRInvariant ) ss = MinSRInvariant
        
        C = C * (ss/Bg)**nn 
      END IF


!
!    Compute Spin : 
!    --------------

        Spin = 0.5 * ( LGrad - TRANSPOSE(LGrad) )
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
!------------------------------------------------------------------------------
      END SUBROUTINE LocalSD      
!------------------------------------------------------------------------------
!        



!------------------------------------------------------------------------------
      END SUBROUTINE AIFlowSolver
!------------------------------------------------------------------------------
