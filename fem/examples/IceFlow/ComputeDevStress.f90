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
! * Module containing a solver for computing deviatoric or Cauchy 
! *          stress from flow solution 
! * 2D SDOFs = 4 (S11, S22, S33, S12)
! * 3D SDOFs = 6 (S11, S22, S33, S12, S23, S31)
! * Keywords : Cauchy (Logical), 
! *            Flow Solver Name (AIFlow, Flow Solution, Porous, ...)
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
! *
! *****************************************************************************/
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE ComputeDevStress( Model,Solver,dt,TransientSimulation )
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
     INTEGER :: dim

     TYPE(ValueList_t),POINTER :: Material, BC
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement

     REAL(KIND=dp) :: RelativeChange, UNorm, PrevUNorm, &
                      s, Wn(7), MinSRInvariant

     REAL(KIND=dp), ALLOCATABLE :: Basis(:),ddBasisddx(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: dBasisdx(:,:)
     REAL(KIND=dp) :: u, v, w, detJ
     
     LOGICAL :: stat, CSymmetry 

     INTEGER :: NewtonIter, NonlinearIter

     TYPE(Variable_t), POINTER :: StressSol, TempSol, FabricVariable, FlowVariable 

     REAL(KIND=dp), POINTER :: Temperature(:), Stress(:), &
           ForceVector(:),  NodalStress(:), StressComp(:), &
           FabricValues(:), FlowValues(:)


     INTEGER, POINTER :: TempPerm(:), StressPerm(:), NodeIndexes(:), &
                         FabricPerm(:), FlowPerm(:)

     INTEGER :: body_id
     INTEGER :: old_body = -1

     LOGICAL :: Isotropic, AllocationsDone = .FALSE.,  &
                Requal0
     LOGICAL :: GotIt,  Cauchy = .FALSE.

     REAL(KIND=dp) :: FabricGrid(4878)

     INTEGER, DIMENSION(3:8) :: NnodeLinear=(/3, 4, 4, 5, 6, 8/)
           
     REAL(KIND=dp), ALLOCATABLE:: LocalMassMatrix(:,:), &
       LocalStiffMatrix(:,:), LocalForce(:), &
       LocalTemperature(:), & 
       K1(:), K2(:), E1(:), &
       E2(:), E3(:), LocalP(:),  &
       LocalVelo(:,:), LocalFluidity(:)
            
     INTEGER :: NumberOfBoundaryNodes
     INTEGER, POINTER :: BoundaryReorder(:)

     REAL(KIND=dp), POINTER :: BoundaryNormals(:,:), &
         BoundaryTangent1(:,:), BoundaryTangent2(:,:)
     CHARACTER(LEN=MAX_NAME_LEN) :: viscosityFile, FlowSolverName

     REAL(KIND=dp) :: at, at0, CPUTime, RealTime


!------------------------------------------------------------------------------
     SAVE NumberOfBoundaryNodes, BoundaryReorder, BoundaryNormals, &
              BoundaryTangent1, BoundaryTangent2, FabricGrid, viscosityFile

     SAVE Basis, dBasisdx, ddBasisddx
     SAVE LocalMassMatrix, LocalStiffMatrix, LocalForce, &
          ElementNodes,  LocalTemperature, &
          Isotropic, AllocationsDone,  &
          K1, K2, E1, E2, E3, Wn, MinSRInvariant, old_body, &
          LocalFluidity, Cauchy
     SAVE LocalVelo, LocalP, dim

!------------------------------------------------------------------------------
!  Read the name of the Flow Solver (NS, AIFlow, Porous, ...)
!------------------------------------------------------------------------------
      
     FlowSolverName = GetString( Solver % Values, 'Flow Solver Name', GotIt )    
     IF (.NOT.Gotit) FlowSolverName = 'aiflow'
     FlowVariable => VariableGet( Solver % Mesh % Variables, FlowSolverName )
     IF ( ASSOCIATED( FlowVariable ) ) THEN
       FlowPerm    => FlowVariable % Perm
       FlowValues  => FlowVariable % Values
     ELSE
       CALL Info('ComputeDevStress', &
                      & 'No variable for velocity associated.', Level=4)
     END IF
!              
!------------------------------------------------------------------------------
!  Read constants from constants section of SIF file
!------------------------------------------------------------------------------
      Wn(7) = GetConstReal( Model % Constants, 'Gas Constant', GotIt )
      IF (.NOT.GotIt) THEN
        WRITE(Message,'(A)') 'VariableGas Constant  not found. &
                     &Setting to 8.314'
        CALL INFO('ComputeDevStress', Message, level=20)
        Wn(7) = 8.314
      ELSE
        WRITE(Message,'(A,F10.4)') 'Gas Constant = ',   Wn(7)
        CALL INFO('ComputeDevStress', Message , level = 20)
      END IF

      Cauchy = ListGetLogical( Solver % Values , 'Cauchy', Gotit )
      IF (.NOT.Gotit) THEN
          Cauchy = .FALSE.
           WRITE(Message,'(A)') 'Cauchy set to False'
           CALL INFO('ComputeDevStress', Message, Level = 20)
      END IF
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
      IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

      StressSol => Solver % Variable
      StressPerm => StressSol % Perm
      STDOFs =  StressSol % DOFs
      Stress => StressSol % Values

      LocalNodes = COUNT( StressPerm > 0 )
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
                     LocalTemperature,     &
                     K1, K2, E1, E2, E3,   &
                     LocalVelo, LocalP,    &                      
                     Basis, ddBasisddx,    &
                     dBasisdx,             &
                     LocalMassMatrix,      &
                     LocalStiffMatrix,     &
                     LocalForce,           &
                     LocalFluidity )
       END IF

       ALLOCATE( ElementNodes % x( N ), &
                 ElementNodes % y( N ), &
                 ElementNodes % z( N ), &
                 LocalTemperature( N ), &
                 K1( N ), K2( N ), E1( N ), E2( N ), E3( N ), &
                 LocalVelo( 3,N ), LocalP( N ), &                                     
                 Basis( 2*N ),ddBasisddx(1,1,1), dBasisdx( 2*N,3 ), &
                 LocalMassMatrix( 2*STDOFs*N,2*STDOFs*N ),  &
                 LocalStiffMatrix( 2*STDOFs*N,2*STDOFs*N ),  &
                 LocalForce( 2*STDOFs*N ),  &
                 LocalFluidity(N), STAT=istat )

       IF ( istat /= 0 ) THEN
          CALL Fatal( 'ComputeDevStress', 'Memory allocation error.' )
       END IF
!------------------------------------------------------------------------------

       AllocationsDone = .TRUE.
      END IF


!------------------------------------------------------------------------------
      NonlinearIter = 1
      DO iter=1,NonlinearIter

       at  = CPUTime()
       at0 = RealTime()

       CALL Info( 'ComputeDevStress', ' ', Level=4 )
       CALL Info( 'ComputeDevStress', ' ', Level=4 )
       CALL Info( 'ComputeDevStress', ' ', Level=4 )
       CALL Info( 'ComputeDevStress', ' ', Level=4 )
       CALL Info( 'ComputeDevStress', 'Starting assembly...',Level=4 )

!------------------------------------------------------------------------------
       CALL DefaultInitialize()
!------------------------------------------------------------------------------
       DO t=1,Solver % NumberOFActiveElements

         IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ',  &
             INT(100.0 - 100.0 * (Solver % NumberOfActiveElements-t) / &
             (1.0*Solver % NumberOfActiveElements)), ' % done'
           CALL Info( 'ComputeDevStress', Message, Level=5 )
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
         SELECT CASE( FlowSolverName)

         CASE('aiflow')

          body_id = CurrentElement % BodyId
          IF ((body_id /= old_body) .OR. (t==1)) THEN 
             old_body = body_id
             CALL  GetMaterialDefs()
          END IF

          LocalFluidity(1:n) = ListGetReal( Material, &
                         'Fluidity Parameter', n, NodeIndexes, GotIt )
          IF (.NOT.GotIt) THEN
            WRITE(Message,'(A)') 'Variable Fluidity Parameter not found. &
                            &Setting to 1.0'
            CALL INFO('ComputeStress', Message, Level = 20)
            LocalFluidity(1:n) = 1.0_dp
          END IF


!------------------------------------------------------------------------------
!        Get element local stiffness & mass matrices
!------------------------------------------------------------------------------
          LocalTemperature = 0.0_dp
          IF ( ASSOCIATED(TempSol) ) THEN
            DO i=1,n
              k = TempPerm(NodeIndexes(i))
              LocalTemperature(i) = Temperature(k)
            END DO
          ELSE
            LocalTemperature(1:n) = 0.0_dp
          END IF

! fabric not needed if isotropic
          IF(.NOT.Isotropic) THEN
            K1(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 1 )
            K2(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 2 )
            E1(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 3 )
            E2(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 4 )
            E3(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 5 )
          ENDIF

          LocalVelo = 0.0_dp
          DO i=1, dim
             LocalVelo(i,1:n) = FlowValues((dim+1)*(FlowPerm(NodeIndexes(1:n))-1) + i)
          END DO
          LocalP(1:n) = FlowValues((dim+1)*FlowPerm(NodeIndexes(1:n)))

          CALL LocalAIFlowMatrix( LocalMassMatrix, LocalStiffMatrix, &
              LocalForce,  K1, K2, E1, E2, E3, LocalVelo, LocalP, &
              LocalTemperature, LocalFluidity, CurrentElement, n, &
              ElementNodes, Wn, MinSRInvariant, Isotropic, Cauchy)

!!!! Restricted to the Power Law case
         CASE('flow solution')
          
          body_id = CurrentElement % BodyId
          IF ((body_id /= old_body) .OR. (t==1)) THEN 
             old_body = body_id
             CALL  GetMaterialDefsNS()
          END IF

          LocalFluidity(1:n) = ListGetReal( Material, &
                         'Viscosity', n, NodeIndexes, GotIt )
          IF (.NOT.GotIt) THEN
            WRITE(Message,'(A)') 'Variable Viscosity not found. &
                            &Setting to 1.0'
            CALL INFO('ComputeStress', Message, Level = 20)
            LocalFluidity(1:n) = 1.0_dp
          END IF
! B = 1/ eta^n
          LocalFluidity = 1.0_dp / (LocalFluidity**Wn(2)) 
          LocalVelo = 0.0_dp
          DO i=1, dim
             LocalVelo(i,1:n) = FlowValues((dim+1)*(FlowPerm(NodeIndexes(1:n))-1) + i)
          END DO
          LocalP(1:n) = FlowValues((dim+1)*FlowPerm(NodeIndexes(1:n)))

          CALL LocalNSMatrix( LocalMassMatrix, LocalStiffMatrix, &
              LocalForce,  LocalVelo, LocalP, &
              LocalFluidity, CurrentElement, n, &
              ElementNodes, Wn, MinSRInvariant, Cauchy)

! to be done
         CASE('porous')
         

         END SELECT 

!------------------------------------------------------------------------------
!        Update global matrices from local matrices 
!------------------------------------------------------------------------------
!        IF ( TransientSimulation ) THEN
!           CALL Default1stOrderTime( LocalMassMatrix, &
!                LocalStiffMatrix, LocalForce )
!        END IF
         CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )
      END DO

      CALL Info( 'ComputeDevStress', 'Assembly done', Level=4 )


      CALL DefaultFinishAssembly()

!------------------------------------------------------------------------------
!     Dirichlet boundary conditions
!------------------------------------------------------------------------------
      CALL DefaultDirichletBCs()

!------------------------------------------------------------------------------

      CALL Info( 'ComputeDevStress', 'Set boundaries done', Level=4 )

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      PrevUNorm = UNorm

      UNorm = DefaultSolve()

      unorm = 0
      n = Solver % Variable % DOFs
      DO i=1,n-1
        unorm = unorm + SUM( solver % variable % values(i::n)**2 )
      END DO
      unorm = SQRT( unorm )
      solver % variable % norm = unorm

      IF ( PrevUNorm + UNorm /= 0.0d0 ) THEN
         RelativeChange = 2.0d0 * ABS( PrevUNorm - UNorm) / ( PrevUnorm + UNorm)
      ELSE
         RelativeChange = 0.0d0
      END IF

      WRITE( Message, * ) 'Result Norm   : ',UNorm, PrevUNorm
      CALL Info( 'ComputeDevStress', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'ComputeDevStress', Message, Level=4 )


!------------------------------------------------------------------------------
    END DO ! of nonlinear iter
!------------------------------------------------------------------------------

      
CONTAINS

      SUBROUTINE GetMaterialDefs()
      ! check if we are isotropic or not
      Isotropic = ListGetLogical( Material , 'Isotropic',Gotit )
      IF (.NOT.Gotit) THEN
          Isotropic = .FALSE.
           WRITE(Message,'(A)') 'Isotropic set to False'
           CALL INFO('ComputeDevStress', Message, Level = 20)
      ELSE
           IF ( (ASSOCIATED( FabricVariable )).AND.Isotropic ) THEN
              WRITE(Message,'(A)') 'Be carefull Isotropic is true &
                           & and Fabric is defined!'
              CALL INFO('ComputeDevStress', Message, Level = 1)
           END IF
      END IF

      IF (.NOT.Isotropic) THEN
        ! Get the viscosity file and store the viscosities into FabricGrid
         viscosityFile = ListGetString( Material ,'Viscosity File',GotIt )
         IF (.NOT.GotIt) THEN
            WRITE(Message,'(3A)') &
                      'Viscosity File ', viscosityFile, ' not found'
           CALL FATAL('ComputeDevStress',Message)
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
         CALL INFO('ComputeDevStress', Message, Level = 20)
         Wn(2) = 1.0_dp
      ELSE
       WRITE(Message,'(A,F10.4)') 'Powerlaw Exponent = ',   Wn(2)
       CALL INFO('ComputeDevStress', Message, Level = 20)
       END IF

      Wn(3) = ListGetConstReal( Material, 'Activation Energy 1', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Activation Energy 1 not found.&
                            & Setting to 1.0'
         CALL INFO('ComputeDevStress', Message, Level = 20)
         Wn(3) = 1.0_dp
      ELSE
         WRITE(Message,'(A,F10.4)') 'Activation Energy 1 = ',   Wn(3)
         CALL INFO('ComputeDevStress', Message, Level = 20)
      END IF

      Wn(4) = ListGetConstReal( Material, 'Activation Energy 2', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Activation Energy 2 not found. &
                               &Setting to 1.0'
         CALL INFO('ComputeDevStress', Message, Level = 20)
         Wn(4) = 1.0_dp
      ELSE
         WRITE(Message,'(A,F10.4)') 'Activation Energy 2 = ',   Wn(4)
         CALL INFO('ComputeDevStress', Message, Level = 20)
      END IF

      Wn(5) = ListGetConstReal(Material, 'Reference Temperature', GotIt)
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Reference Temperature not found. &
                               &Setting to -10.0 (Celsius)'
         CALL INFO('ComputeDevStress', Message, Level = 20)
         Wn(5) = -10.0_dp
      ELSE
         WRITE(Message,'(A,F10.4)') 'Reference Temperature = ',   Wn(5)
         CALL INFO('ComputeDevStress', Message, Level = 20)
      END IF

      Wn(6) = ListGetConstReal( Material, 'Limit Temperature', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Limit Temperature not found. &
                               &Setting to -10.0 (Celsius)'
         CALL INFO('ComputeDevStress', Message, Level = 20)
         Wn(6) = -10.0_dp
      ELSE
         WRITE(Message,'(A,F10.4)') 'Limit Temperature = ',   Wn(6)
         CALL INFO('ComputeDevStress', Message, Level = 20)
      END IF

! Get the Minimum value of the Effective Strain rate 
      MinSRInvariant = 1.0e-10_dp
      IF ( Isotropic .AND. (Wn(2) > 1.0 ) ) THEN
        MinSRInvariant =  &
             ListGetConstReal( Material, 'Min Second Invariant', GotIt )
        IF (.NOT.GotIt) THEN
          WRITE(Message,'(A)') 'Variable Min Second Invariant not &
                    &found. Setting to 1.0D-10 )'
          CALL INFO('ComputeDevStress', Message, Level = 20)
        ELSE
          WRITE(Message,'(A,E14.8)') 'Min Second Invariant = ', MinSRInvariant
          CALL INFO('ComputeDevStress', Message, Level = 20)
        END IF
      END IF

!------------------------------------------------------------------------------
      END SUBROUTINE GetMaterialDefs
!------------------------------------------------------------------------------

      SUBROUTINE GetMaterialDefsNS()
      ! check if we are isotropic or not


      Wn(2) = ListGetConstReal( Material , 'Viscosity Exponent', GotIt )

      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable  Powerlaw Exponent not found. &
                                    & Setting to 1.0'
         CALL INFO('ComputeDevStress', Message, Level = 20)
         Wn(2) = 1.0_dp
      ELSE
         Wn(2) = 1.0_dp / Wn(2)
         WRITE(Message,'(A,F10.4)') 'Viscosity Exponent = ',   Wn(2)
         CALL INFO('ComputeDevStress', Message, Level = 20)
      END IF

! Get the Minimum value of the Effective Strain rate 
      MinSRInvariant = 1.0e-10_dp
      IF ( Wn(2) > 1.0  ) THEN
        MinSRInvariant =  &
             ListGetConstReal( Material, 'Critical Shear Rate', GotIt )
        IF (.NOT.GotIt) THEN
          WRITE(Message,'(A)') 'Variable Critical Shear Rate not &
                    &found. Setting to 1.0D-10 )'
          CALL INFO('ComputeDevStress', Message, Level = 20)
        ELSE
          WRITE(Message,'(A,E14.8)') 'Critical Shear Rate = ', MinSRInvariant
          CALL INFO('ComputeDevStress', Message, Level = 20)
        END IF
      END IF

!------------------------------------------------------------------------------
      END SUBROUTINE GetMaterialDefsNS
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      SUBROUTINE LocalAIFlowMatrix( MassMatrix, StiffMatrix, ForceVector, &
              NodalK1, NodalK2, NodalEuler1, NodalEuler2, &
              NodalEuler3, NodalVelo, NodalP, NodalTemperature, NodalFluidity, &
              Element, n, Nodes, Wn, MinSRInvariant, Isotropic, Cauchy )
              
              
              
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: StiffMatrix(:,:), MassMatrix(:,:)
     REAL(KIND=dp) ::  NodalVelo(:,:)
     REAL(KIND=dp) :: Wn(7), MinSRInvariant
     REAL(KIND=dp), DIMENSION(:) :: ForceVector, NodalK1, NodalK2, &
             NodalEuler1, NodalEuler2, NodalEuler3, NodalTemperature, &
             NodalFluidity, NodalP
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
     LOGICAL :: Isotropic, Cauchy
     INTEGER :: n
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(2*n),ddBasisddx(1,1,1)
     REAL(KIND=dp) :: dBasisdx(2*n,3),detJ, pBasis(n)

     REAL(KIND=dp) :: K1, K2, Euler1, Euler2, Euler3
     REAL(kind=dp) :: Bg, BGlenT, Stress(6)

     REAL(KIND=dp), DIMENSION(4,4) :: A,M
     REAL(KIND=dp) :: Temperature,  C(6,6), Pressure
     REAL(KIND=dp) :: nn, ss, LGrad(3,3), SR(3,3), D(6), epsi

     INTEGER :: i, j, k, p, q, t, dim, cc, NBasis, ind(3), LinearBasis

     REAL(KIND=dp) :: s, u, v, w, Radius, eta
  
     REAL(KIND=dp) :: dDispldx(3,3), ai(3), Angle(3), a2(6)
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     INTEGER, POINTER :: EdgeMap(:,:)
     INTEGER :: N_Integ, nd

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat, CSymmetry, P2P1, MiniElement

     TYPE(ElementType_t), POINTER :: SaveType, LinearType

     INTERFACE
      SUBROUTINE R2Ro(a2,dim,ai,angle)
         USE Types
         REAL(KIND=dp),INTENT(in) :: a2(6)
         INTEGER :: dim
         REAL(KIND=dp),INTENT(out) :: ai(3), Angle(3)
      END SUBROUTINE R2Ro
                 
      SUBROUTINE OPILGGE_ai(ai,Angle,Tc,W,etaI,eta36)
          USE Types
          REAL(KIND=dp) :: ai(3), Angle(3), Tc, W(7), EtaI(:), Eta36(6,6)
        END SUBROUTINE OPILGGE_ai

      END INTERFACE

!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()
      cc=2*dim


      ForceVector = 0.0_dp
      StiffMatrix = 0.0_dp
      MassMatrix  = 0.0_dp

      IntegStuff = GaussPoints( Element )

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
       stat = ElementInfo(Element,Nodes,u,v,w,detJ, &
                 Basis,dBasisdx,ddBasisddx,.FALSE.,.FALSE.)


       s = detJ * S_Integ(t)


       Radius = SUM( Nodes % x(1:n) * Basis(1:n) )
       CSymmetry = CurrentCoordinateSystem() == AxisSymmetric
       IF ( CSymmetry ) s = s * Radius
!
!     Temperature at the integration point
!
       Temperature = SUM( NodalTemperature(1:n)*Basis(1:n) )
       Wn(1) = SUM( NodalFluidity(1:n)*Basis(1:n) )
       Bg = BGlenT(Temperature,Wn)
       Pressure = SUM( NodalP(1:n)*Basis(1:n) )
!
! Strain-Rate
!
        LGrad = MATMUL( NodalVelo(:,1:n), dBasisdx(1:n,:) )
        SR = 0.5 * ( LGrad + TRANSPOSE(LGrad) )
        IF ( CSymmetry ) THEN
          SR(1,3) = 0.0_dp
          SR(2,3) = 0.0_dp
          SR(3,1) = 0.0_dp
          SR(3,2) = 0.0_dp
          SR(3,3) = 0.0_dp
          IF ( Radius > 10*AEPS ) THEN
            SR(3,3) = SUM( Nodalvelo(1,1:n) * Basis(1:n) ) /Radius
                 
          END IF
          epsi = SR(1,1)+SR(2,2)+SR(3,3)
          DO i=1,3   
            SR(i,i) = SR(i,i) - epsi/3.0_dp
          END DO
        ELSE
          epsi = SR(1,1)+SR(2,2)+SR(3,3)
          DO i=1,dim 
            SR(i,i) = SR(i,i) - epsi/dim
          END DO
        END IF

! if not isotropic use GOLF
      C = 0.0_dp
      IF (.NOT.Isotropic) THEN
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
         Do i = 1,3
           C(i,i) = 2.0_dp / Bg
           C(i+3,i+3) = 1.0_dp / Bg
         END DO
      ENDIF
      

!
! Case non-linear and ISOTROPIC
! -----------------------------

      IF ( Isotropic .AND. (Wn(2) > 1.0) ) THEN
      

        ss = 0.0_dp
        DO i = 1, 3
          DO j = 1, 3
            ss = ss + SR(i,j)**2
          END DO
        END DO
        nn =(1.0_dp - Wn(2))/Wn(2)

        ss = SQRT(2.0_dp * ss)
        IF (ss < MinSRInvariant ) ss = MinSRInvariant

        C  = C * (ss/Bg)**nn 

      END IF

!
!    Compute deviatoric stresses or Cauchy stresses: 
!    ----------------------------
      D(1) = SR(1,1)
      D(2) = SR(2,2)
      D(3) = SR(3,3)
      D(4) = 2. * SR(1,2)
      D(5) = 2. * SR(2,3)
      D(6) = 2. * SR(3,1)
      
      Stress = 0.0_dp
      DO k = 1, cc      
       DO j = 1, cc       
        Stress(k) = Stress(k) + C(k,j) * D(j)
       END DO
      END DO
      IF (Cauchy) THEN
        DO k=1,3
          Stress(k) = Stress(k) - Pressure
        END DO 
      END IF
      
      DO i = 1, cc
        DO p=1,n         
          DO q=1,n        
            StiffMatrix(cc*(p-1)+i,cc*(q-1)+i ) =  &
               StiffMatrix(cc*(p-1)+i,cc*(q-1)+i ) + s*Basis(q)*Basis(p)
          END DO

             ForceVector(cc*(p-1)+i) =  &
                     ForceVector(cc*(p-1)+i) + s*Stress(i)*Basis(p) 
        END DO
      END DO

      END DO 

!------------------------------------------------------------------------------
      END SUBROUTINE LocalAIFlowMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      SUBROUTINE LocalNSMatrix( MassMatrix, StiffMatrix, ForceVector, &
               NodalVelo, NodalP, NodalFluidity, &
              Element, n, Nodes, Wn, MinSRInvariant, Cauchy )
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: StiffMatrix(:,:), MassMatrix(:,:)
     REAL(KIND=dp) ::  NodalVelo(:,:)
     REAL(KIND=dp) :: Wn(7), MinSRInvariant
     REAL(KIND=dp), DIMENSION(:) :: ForceVector,  &
                              NodalFluidity, NodalP
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
     LOGICAL ::  Cauchy
     INTEGER :: n
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(2*n),ddBasisddx(1,1,1)
     REAL(KIND=dp) :: dBasisdx(2*n,3),detJ, pBasis(n)

     REAL(kind=dp) :: Bg,  Stress(6)

     REAL(KIND=dp), DIMENSION(4,4) :: A,M
     REAL(KIND=dp) :: Temperature,  C(6,6), Pressure
     REAL(KIND=dp) :: nn, ss, LGrad(3,3), SR(3,3), D(6), epsi

     INTEGER :: i, j, k, p, q, t, dim, cc, NBasis, ind(3), LinearBasis

     REAL(KIND=dp) :: s, u, v, w, Radius, eta
  
     REAL(KIND=dp) :: dDispldx(3,3), ai(3), Angle(3), a2(6)
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     INTEGER, POINTER :: EdgeMap(:,:)
     INTEGER :: N_Integ, nd

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat, CSymmetry, P2P1, MiniElement

     TYPE(ElementType_t), POINTER :: SaveType, LinearType


!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()
      cc=2*dim


      ForceVector = 0.0_dp
      StiffMatrix = 0.0_dp
      MassMatrix  = 0.0_dp

      IntegStuff = GaussPoints( Element )

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
       stat = ElementInfo(Element,Nodes,u,v,w,detJ, &
                 Basis,dBasisdx,ddBasisddx,.FALSE.,.FALSE.)


       s = detJ * S_Integ(t)


       Radius = SUM( Nodes % x(1:n) * Basis(1:n) )
       CSymmetry = CurrentCoordinateSystem() == AxisSymmetric
       IF ( CSymmetry ) s = s * Radius
!
!     Temperature at the integration point
!
       Wn(1) = SUM( NodalFluidity(1:n)*Basis(1:n) )
       Bg = Wn(1)                       
       Pressure = SUM( NodalP(1:n)*Basis(1:n) )
!
! Strain-Rate
!
        LGrad = MATMUL( NodalVelo(:,1:n), dBasisdx(1:n,:) )
        SR = 0.5 * ( LGrad + TRANSPOSE(LGrad) )
        IF ( CSymmetry ) THEN
          SR(1,3) = 0.0_dp
          SR(2,3) = 0.0_dp
          SR(3,1) = 0.0_dp
          SR(3,2) = 0.0_dp
          SR(3,3) = 0.0_dp
          IF ( Radius > 10*AEPS ) THEN
            SR(3,3) = SUM( Nodalvelo(1,1:n) * Basis(1:n) ) /Radius
                 
          END IF
          epsi = SR(1,1)+SR(2,2)+SR(3,3)
          DO i=1,3   
            SR(i,i) = SR(i,i) - epsi/3.0_dp
          END DO
        ELSE
          epsi = SR(1,1)+SR(2,2)+SR(3,3)
          DO i=1,dim 
            SR(i,i) = SR(i,i) - epsi/dim
          END DO
        END IF

        C = 0.0_dp
         Do i = 1,3
           C(i,i) = 2.0_dp / Bg
           C(i+3,i+3) = 1.0_dp / Bg
         END DO
      

!
! Case non-linear and ISOTROPIC
! -----------------------------

      IF ( Wn(2) > 1.0 ) THEN
        ss = 0.0_dp
        DO i = 1, 3
          DO j = 1, 3
            ss = ss + SR(i,j)**2
          END DO
        END DO
        nn =(1.0_dp - Wn(2))/Wn(2)

        ss = SQRT(2.0_dp * ss)
        IF (ss < MinSRInvariant ) ss = MinSRInvariant

        C  = C * (ss/Bg)**nn 
      END IF

!
!    Compute deviatoric stresses or Cauchy stresses: 
!    ----------------------------
      D(1) = SR(1,1)
      D(2) = SR(2,2)
      D(3) = SR(3,3)
      D(4) = 2. * SR(1,2)
      D(5) = 2. * SR(2,3)
      D(6) = 2. * SR(3,1)
      
      Stress = 0.0_dp
      DO k = 1, cc      
       DO j = 1, cc       
        Stress(k) = Stress(k) + C(k,j) * D(j)
       END DO
      END DO
      IF (Cauchy) THEN
        DO k=1,3
          Stress(k) = Stress(k) - Pressure
        END DO 
      END IF
      
      DO i = 1, cc
        DO p=1,n         
          DO q=1,n        
            StiffMatrix(cc*(p-1)+i,cc*(q-1)+i ) =  &
               StiffMatrix(cc*(p-1)+i,cc*(q-1)+i ) + s*Basis(q)*Basis(p)
          END DO

             ForceVector(cc*(p-1)+i) =  &
                     ForceVector(cc*(p-1)+i) + s*Stress(i)*Basis(p) 
        END DO
      END DO

      END DO 

!------------------------------------------------------------------------------
      END SUBROUTINE LocalNSMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      END SUBROUTINE ComputeDevStress
!------------------------------------------------------------------------------
