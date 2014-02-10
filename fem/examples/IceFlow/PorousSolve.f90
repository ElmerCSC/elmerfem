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
! * Module containing a solver for (primarily thermal) Porous material flow
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
! *                Modified by:  OG 
! *
! *       Date of modification: 
! *
! *****************************************************************************/
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE PorousSolver( Model,Solver,dt,TransientSimulation )
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
!
!    Solve:  -p,j + Sij,j + D.fi = 0 and ui,i + Kcp.p = 0
!    with 
!       Sij = 2 eta ( Eij - Ekk / 3 delta_ij) 
!       Kcp = b(D) B^(1/n) ED^{(n-1)/n}
!       ED^2 = 2 eij.eij/a(D) + Ekk^2/b(D)
!       eta = 2/a(D) B^{-1/n} ED^{(1-n)/n}
!       The relative density D is an input for this solver
!   From Gagliardini and Meyssonnier, 1997. Annals of Glaciology 24.
!
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
         Normal(3), NewtonTol, NonlinearTol, s, Wn(2), MinSRInvariant
         

     REAL(KIND=dp)  :: NodalStresses(3,3), &
       NodalStrainRate(3,3),  NodalSpin(3,3)   

     REAL(KIND=dp), ALLOCATABLE :: Basis(:),ddBasisddx(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: dBasisdx(:,:), SlipCoeff(:,:)
     REAL(KIND=dp) :: u,v,w,detJ
     
     LOGICAL :: stat, CSymmetry 
       
     INTEGER, PARAMETER :: INDi(1:6) = (/ 1, 2, 3, 1, 2, 3 /) ,&
           INDj(1:6)=(/ 1, 2, 3, 2, 3, 1 /)

     INTEGER :: NewtonIter, NonlinearIter

     TYPE(Variable_t), POINTER :: PorousSol, DensityVariable
     TYPE(Variable_t), POINTER :: SpinVar
     REAL(KIND=dp), POINTER :: SpinValues(:)
     INTEGER, POINTER :: SpinPerm(:)

     TYPE(Variable_t), POINTER :: DevStressVar 
     REAL(KIND=dp), POINTER :: DSValues(:)
     INTEGER, POINTER :: DSPerm(:)

     TYPE(Variable_t), POINTER :: StrainRateVar 
     REAL(KIND=dp), POINTER :: SRValues(:)
     INTEGER, POINTER :: SRPerm(:)

     REAL(KIND=dp), POINTER ::  Porous(:), Work(:,:), &
           ForceVector(:),  NodalPorous(:), PorousComp(:), &
           DensityValues(:)

     CHARACTER(LEN=MAX_NAME_LEN) :: EquationName

     INTEGER, POINTER :: PorousPerm(:), NodeIndexes(:), &
                         DensityPerm(:)

     INTEGER :: PorousType
     LOGICAL :: GotForceBC, GotIt, NewtonLinearization = .FALSE., &
                NormalTangential=.FALSE.

     INTEGER :: body_id,bf_id
     INTEGER :: old_body = -1
     LOGICAL :: AllocationsDone = .FALSE., FreeSurface, &
                Requal0
           
     REAL(KIND=dp), ALLOCATABLE:: LocalMassMatrix(:,:), &
       LocalStiffMatrix(:,:), LoadVector(:,:), LocalForce(:), &
       Alpha(:,:), Beta(:), & 
       BoundaryDispl(:), LocalDensity(:), &
       TimeForce(:), RefS(:), RefD(:), RefSpin(:), &
       LocalVelo(:,:), LocalFluidity(:), Localfa(:), Localfb(:)
            
     INTEGER :: NumberOfBoundaryNodes
     INTEGER, POINTER :: BoundaryReorder(:)

     REAL(KIND=dp) :: Bu, Bv, Bw, RM(3,3)
     REAL(KIND=dp), POINTER :: BoundaryNormals(:,:), &
         BoundaryTangent1(:,:), BoundaryTangent2(:,:)

     REAL(KIND=dp) :: at, at0, CPUTime, RealTime


!------------------------------------------------------------------------------
     SAVE NumberOfBoundaryNodes,BoundaryReorder,BoundaryNormals, &
              BoundaryTangent1, BoundaryTangent2

     SAVE TimeForce, Basis, dBasisdx, ddBasisddx
     SAVE LocalMassMatrix, LocalStiffMatrix, LoadVector, &
       LocalForce, ElementNodes, Alpha, Beta,  &
       AllocationsDone,BoundaryDispl, &
       NodalPorous, LocalDensity, Wn, MinSRInvariant, old_body, &
       LocalFluidity, Localfa, Localfb

     SAVE RefD, RefS, RefSpin, LocalVelo, SlipCoeff 
              
!------------------------------------------------------------------------------
!  Read constants from constants section of SIF file
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
      IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

      PorousSol => Solver % Variable
      PorousPerm => PorousSol % Perm
      STDOFs =  PorousSol % DOFs
      Porous => PorousSol % Values

      LocalNodes = COUNT( PorousPerm > 0 )
      IF ( LocalNodes <= 0 ) RETURN


      DensityVariable => &
              VariableGet(Solver % Mesh %Variables,'Relative Density')
      IF ( ASSOCIATED( DensityVariable ) ) THEN
       DensityPerm    => DensityVariable % Perm
       DensityValues => DensityVariable % Values
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
                     LocalVelo,            &
                     LocalDensity,         &
                     LocalForce,           &
                     RefD, RefS, RefSpin,  &
                     LocalMassMatrix,      &
                     LocalStiffMatrix,     &
                     LoadVector, Alpha, Beta, &
                     SlipCoeff, LocalFluidity , &
                     Localfa, Localfb )
       END IF

       ALLOCATE( ElementNodes % x( N ), &
                 ElementNodes % y( N ), &
                 ElementNodes % z( N ), &
                 BoundaryDispl( N ), &
                 LocalDensity( N ), &
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
                 SlipCoeff(3,N), LocalFluidity(N), &
                 Localfa( N ), Localfb( N ),  STAT=istat )

       IF ( istat /= 0 ) THEN
          CALL Fatal( 'PorousSolve', 'Memory allocation error.' )
       END IF
!------------------------------------------------------------------------------
!    Check for normal/tangetial coordinate system defined velocities
!------------------------------------------------------------------------------
       CALL CheckNormalTangentialBoundary( Model, &
        'Normal-Tangential Porous', NumberOfBoundaryNodes, &
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

       CALL Info( 'PorousSolve', ' ', Level=4 )
       CALL Info( 'PorousSolve', ' ', Level=4 )
       CALL Info( 'PorousSolve', &
                   '-------------------------------------',Level=4 )
       WRITE( Message, * ) 'POROUS MATERIAL FLOW SOLVER ITERATION', iter
       CALL Info( 'PorousSolve', Message,Level=4 )
       CALL Info( 'PorousSolve', &
                   '-------------------------------------',Level=4 )
       CALL Info( 'PorousSolve', ' ', Level=4 )
       CALL Info( 'PorousSolve', 'Starting assembly...',Level=4 )
!------------------------------------------------------------------------------
!      Compute average normals for boundaries having the normal & tangential
!      field components specified on the boundaries
!------------------------------------------------------------------------------
       IF ( (iter == 1 .OR. FreeSurface) .AND. NumberOfBoundaryNodes > 0 ) THEN
          CALL AverageBoundaryNormals( Model, &
               'Normal-Tangential Porous', NumberOfBoundaryNodes, &
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
           CALL Info( 'PorousSolve', Message, Level=5 )
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
         IF (body_id /= old_body) Then 
            old_body = body_id
            Call  GetMaterialDefs()
        END IF

        LocalFluidity(1:n) = ListGetReal( Material, &
                         'Fluidity Parameter', n, NodeIndexes, GotIt )
        IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Fluidity Parameter not found. &
                            &Setting to 1.0'
         CALL INFO('PorousSolve', Message, Level = 20)
         LocalFluidity(1:n) = 1.0
        END IF
!
! Function a(D) (=1 if incompressible)
! 
        Localfa(1:n) = ListGetReal( Material, &
                         'FunctionA', n, NodeIndexes, GotIt )
        IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'FunctionA not found. &
                            &Setting to 1.0'
         CALL INFO('PorousSolve', Message, Level = 20)
         Localfa(1:n) = 1.0
        END IF
!
! Function b(D) (=0 if incompressible)
! 
        Localfb(1:n) = ListGetReal( Material, &
                         'FunctionB', n, NodeIndexes, GotIt )
        IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'FunctionB not found. &
                            &Setting to 0.0'
         CALL INFO('PorousSolve', Message, Level = 20)
         Localfb(1:n) = 0.0
        END IF

!------------------------------------------------------------------------------
!        Set body forces
!------------------------------------------------------------------------------
         LoadVector = 0.0D0

         BodyForce => GetBodyForce()
         IF ( ASSOCIATED( BodyForce ) ) THEN
           LoadVector(1,1:n) = LoadVector(1,1:n) + ListGetReal( &
                   BodyForce, 'Porous Force 1', n, NodeIndexes, gotIt)
           LoadVector(2,1:n) = LoadVector(2,1:n) + ListGetReal( &
                   BodyForce, 'Porous Force 2', n, NodeIndexes, gotIt)
           LoadVector(3,1:n) = LoadVector(3,1:n) + ListGetReal( & 
                   BodyForce, 'Porous Force 3', n, NodeIndexes, gotIt)
         END IF
!------------------------------------------------------------------------------
!        Get element local stiffness & mass matrices
!------------------------------------------------------------------------------

         LocalDensity(1:n)= DensityValues(DensityPerm(NodeIndexes(1:n)))

         LocalVelo = 0.0d0
         DO i=1,STDOFs - 1
            LocalVelo(i,1:n) = Porous( STDOFs*(PorousPerm(NodeIndexes(1:n))-1) + i)
         END DO
         CALL LocalMatrix( LocalMassMatrix, LocalStiffMatrix, &
              LocalForce, LoadVector, LocalDensity, LocalVelo, &
              LocalFluidity, CurrentElement, n, &
              ElementNodes, Wn, MinSRInvariant, Localfa, Localfb )

!------------------------------------------------------------------------------
!        Update global matrices from local matrices 
!------------------------------------------------------------------------------
         IF ( TransientSimulation ) THEN 
           CALL Default1stOrderTime( LocalMassMatrix, &
                LocalStiffMatrix, LocalForce )
         END IF
         CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )
      END DO

      CALL Info( 'PorousSolve', 'Assembly done', Level=4 )

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
                    'Porous Slip Coeff 1', n, NodeIndexes, GotIt )
              GotForceBC = GotForceBC .OR. gotIt

              SlipCoeff(2,1:n) =  ListGetReal( BC, &
                    'Porous Slip Coeff 2', n, NodeIndexes, GotIt )
              GotForceBC = GotForceBC .OR. gotIt

              SlipCoeff(3,1:n) =  ListGetReal( BC, &
                    'Porous Slip Coeff 3', n, NodeIndexes, GotIt )
              GotForceBC = GotForceBC .OR. gotIt

              NormalTangential = ListGetLogical( BC, &
                     'Normal-Tangential Porous', GotIt )
               
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
            LocalMassMatrix = 0.0_dp
            IF ( TransientSimulation ) THEN
              CALL Default1stOrderTime( LocalMassMatrix, & 
                   LocalStiffMatrix, LocalForce )
            END IF
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

      CALL Info( 'PorousSolve', 'Set boundaries done', Level=4 )

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      PrevUNorm = UNorm

      UNorm = DefaultSolve()
      
      UNorm = 0.0_dp
      n = Solver % Variable % DOFs
      DO i=1,n-1
        UNorm = UNorm + SUM( Solver % Variable % Values(i::n)**2 )
      END DO
      UNorm = SQRT(UNorm)
      Solver % Variable % norm = UNorm

      IF ( PrevUNorm + UNorm /= 0.0d0 ) THEN
         RelativeChange = 2.0d0 * ABS( PrevUNorm - UNorm) / ( PrevUnorm + UNorm)
      ELSE
         RelativeChange = 0.0d0
      END IF

      WRITE( Message, * ) 'Result Norm   : ',UNorm, PrevUNorm
      CALL Info( 'PorousSolve', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'PorousSolve', Message, Level=4 )

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
         CALL INFO('PorousSolve', Message, Level = 20)
         LocalFluidity(1:n) = 1.0
        END IF


         ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes(1:n))
         ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes(1:n))
         ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes(1:n))

         
!
! Function a(D) (=1 if incompressible)
! 
        Localfa(1:n) = ListGetReal( Material, &
                         'FunctionA', n, NodeIndexes, GotIt )
        IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'FunctionA not found. &
                            &Setting to 1.0'
         CALL INFO('PorousSolve', Message, Level = 20)
         Localfa(1:n) = 1.0
        END IF
!
! Function b(D) (=0 if incompressible)
! 
        Localfb(1:n) = ListGetReal( Material, &
                         'FunctionB', n, NodeIndexes, GotIt )
        IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'FunctionB not found. &
                            &Setting to 0.0'
         CALL INFO('PorousSolve', Message, Level = 20)
         Localfb(1:n) = 0.0
        END IF

         LocalDensity(1:n) =DensityValues(DensityPerm(NodeIndexes(1:n)))

         LocalVelo = 0.0d0
         DO i=1,STDOFs - 1
            LocalVelo(i,1:n) = Porous( STDOFs*(PorousPerm(NodeIndexes(1:n))-1) + i)
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
                 LocalVelo, LocalFluidity,  &
                LocalDensity, CSymmetry, Basis, dBasisdx, &
                CurrentElement, n, ElementNodes, dim, Wn, &
                MinSRInvariant, Localfa, Localfb )
                

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

      Wn(2) = ListGetConstReal( Material , 'Powerlaw Exponent', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable  Powerlaw Exponent not found. &
                                    & Setting to 1.0'
         CALL INFO('PorousSolve', Message, Level = 20)
         Wn(2) = 1.0
      ELSE
       WRITE(Message,'(A,F10.4)') 'Powerlaw Exponent = ',   Wn(2)
       CALL INFO('PorousSolve', Message, Level = 20)
       END IF

! Get the Minimum value of the Effective Strain rate 
      MinSRInvariant = 10e-10_dp 
      IF ( Wn(2) > 1.0 ) THEN
        MinSRInvariant =  &
             ListGetConstReal( Material, 'Min Second Invariant', GotIt )
        IF (.NOT.GotIt) THEN
          WRITE(Message,'(A)') 'Variable Min Second Invariant not &
                    &found. Setting to 100.0*AEPS )'
          CALL INFO('PorousSolve', Message, Level = 20)
        ELSE
          WRITE(Message,'(A,E14.8)') 'Min Second Invariant = ', MinSRInvariant
          CALL INFO('PorousSolve', Message, Level = 20)
        END IF
      END IF

!------------------------------------------------------------------------------
      END SUBROUTINE GetMaterialDefs
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      SUBROUTINE LocalMatrix( MassMatrix, StiffMatrix, ForceVector, &
              LoadVector, NodalDensity, NodalVelo, &
              NodalFluidity, Element, n, Nodes, Wn, MinSRInvariant, &
              Nodalfa, Nodalfb )
              
              
              
              
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: StiffMatrix(:,:), MassMatrix(:,:)
     REAL(KIND=dp) :: LoadVector(:,:), NodalVelo(:,:)
     REAL(KIND=dp) :: Wn(2), MinSRInvariant
     REAL(KIND=dp), DIMENSION(:) :: ForceVector, NodalDensity,  &
                    NodalFluidity, Nodalfa, Nodalfb
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
     INTEGER :: n
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(2*n),ddBasisddx(1,1,1)
     REAL(KIND=dp) :: dBasisdx(2*n,3),SqrtElementMetric

     REAL(KIND=dp) :: Force(3), density, fa, fb 
     Real(kind=dp) :: Bg, BGlenT

     REAL(KIND=dp), DIMENSION(4,4) :: A,M
     REAL(KIND=dp) :: Load(3)
     REAL(KIND=dp) :: nn, ss, LGrad(3,3), SR(3,3)

     INTEGER :: i, j, k, p, q, t, dim, NBasis, ind(3)

     REAL(KIND=dp) :: s,u,v,w, Radius
     REAL(KIND=dp) :: ParameterA, ParameterB
  
     REAL(KIND=dp) :: Em, eta, Kcp
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     INTEGER :: N_Integ

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat, CSymmetry


!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()

      ForceVector = 0.0D0
      StiffMatrix = 0.0D0
      MassMatrix  = 0.0D0

!    
!    Integration stuff
!    
      NBasis = 2*n
      IntegStuff = GaussPoints( Element, Element % Type % GaussPoints2 )

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
      stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
            Basis,dBasisdx,ddBasisddx,.FALSE.,.TRUE. )

      s = SqrtElementMetric * S_Integ(t)
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
!    variables at the integration point
!
      Density = SUM( NodalDensity(1:n)*Basis(1:n) )
!
! Nodal value interpolated at the integration point
!
!     fa = SUM( Nodalfa(1:n)*Basis(1:n) )
!     fb = SUM( Nodalfb(1:n)*Basis(1:n) )
!
! Integration point value function of the integration pt density
      fa = ParameterA(Density)
      fb = ParameterB(Density)

      write(*,*)Density,fa,fb

      Wn(1) = SUM( NodalFluidity(1:n)*Basis(1:n) )

      Bg=Wn(1)                      

      CSymmetry = CurrentCoordinateSystem() == AxisSymmetric
      IF ( CSymmetry ) s = s * Radius

!
! Case non-linear calculate E_D^2 = gamma_e^2/fa + E_m^2/fb
! ----------------------------------------------------------
      ss = 1.0_dp
      IF ( Wn(2) > 1.0 ) THEN
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
        END IF

        Em = SR(1,1)+SR(2,2)+SR(3,3)
        DO i = 1, 3
          SR(i,i) = SR(i,i) - Em/3.0
        END DO

        ss = 0.0_dp
        DO i = 1, 3
          DO j = 1, 3
            ss = ss + SR(i,j)**2
          END DO
        END DO
        ss = 2.0*ss / fa
        IF ( fb /= 0.0_dp ) ss = ss + Em**2 / fb

        nn =(1.0 - Wn(2))/Wn(2)

        ss = SQRT(ss)
        IF (ss < MinSRInvariant ) ss = MinSRInvariant
        ss =  ss**nn 
      END IF
!
! Bulk effective viscosity and Compressibility parameter Kcp 
!     
      eta = ss / (fa * Bg**( 1.0 / Wn( 2 ) ) ) 
      Kcp = fb * Bg**( 1.0 / Wn( 2 ) ) / ss 

!
!    Loop over basis functions (of both unknowns and weights)
!

      DO p=1,NBasis
       DO q=1,NBasis

        A = 0.0d0
        M = 0.0d0

         DO i=1,dim
           DO j=1,dim

! terms 2 eta Eij dEij
             A(i,i)=A(i,i) + eta* dbasisdx(q,j)*dbasisdx(p,j)
             A(i,j)=A(i,j) + eta* dbasisdx(q,i)*dbasisdx(p,j)

! terms 2 / 3 eta Eii dEii
             A(i,j)=A(i,j) - 2.0/3.0 * eta* dbasisdx(q,j)*dbasisdx(p,i)

           END DO

! Pressure gradient
            A(i,dim+1) = -dBasisdx(p,i) * Basis(q)

! Continuity equation:
            A(dim+1,i) = A(dim+1,i) + dBasisdx(q,i) * Basis(p)
         END DO

         IF ( CSymmetry ) A(1,dim+1) =  A(1,dim+1) - Basis(p) * Basis(q) / Radius
         IF ( CSymmetry ) A(dim+1,1) =  A(dim+1,1) + Basis(p) * Basis(q) / Radius


         A(dim+1, dim+1) = A(dim+1,dim+1) + Kcp * basis(q) * basis(p)

! Add nodal matrix to element matrix
         DO i=1,dim+1
            DO j=1,dim+1
               StiffMatrix( (dim+1)*(p-1)+i,(dim+1)*(q-1)+j ) =  &
                    StiffMatrix( (dim+1)*(p-1)+i,(dim+1)*(q-1)+j ) + s*A(i,j)
            END DO
         END DO

       END DO

! The righthand side...
        Load = 0.0d0
  
        DO i=1,dim
           Load(i) = Load(i) + Force(i) * Basis(p)
        END DO

        DO i=1,dim
           ForceVector((dim+1)*(p-1)+i) = ForceVector((dim+1)*(p-1)+i) +  &
                     s * Load(i) * Density
        END DO
      END DO
      END DO 

      DO j=n+1,n+Element % bdofs
        i = (dim+1)*j
        StiffMatrix(i,:) = 0.0_dp
        StiffMatrix(:,i) = 0.0_dp
        StiffMatrix(i,i) = 1.0_dp
        MassMatrix(:,i) = 0.0_dp
        MassMatrix(i,:) = 0.0_dp
        ForceVector(i) = 0.0_dp
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
     REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric

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
       stat = ElementInfo( Element, Nodes, u, v, w, SqrtElementMetric, &
                Basis, dBasisdx, ddBasisddx, .FALSE. )

       s = SqrtElementMetric * S_Integ(t)
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
        NodalVelo,  NodalFluidity, Nodaldensity,  &
        CSymmetry, Basis, dBasisdx, Element, n,  Nodes, dim,  Wn, &
        MinSRInvariant, Nodalfa, Nodalfb )
       
!------------------------------------------------------------------------------
!    Subroutine to computre the nodal Strain-Rate, Stress, ...
!------------------------------------------------------------------------------
     LOGICAL ::  CSymmetry 
     INTEGER :: n, dim
     INTEGER :: INDi(6),INDj(6)
     REAL(KIND=dp) :: Stress(:,:), StrainRate(:,:), Spin(:,:)
     REAL(KIND=dp) :: NodalVelo(:,:), NodalFluidity(:)
     REAL(KIND=dp) :: Basis(2*n), ddBasisddx(1,1,1)
     REAL(KIND=dp) :: dBasisdx(2*n,3)
     REAL(KIND=dp) :: detJ
     REAL(KIND=dp) :: NodalDensity(:), Nodalfa(:), Nodalfb(:)
     REAL(KIND=dp) :: u, v, w      
     REAL(KIND=dp) :: Wn(2),  D(6), MinSRInvariant
     LOGICAL :: Isotropic
      
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
     LOGICAL :: stat
     INTEGER :: i,j,k,p,q
     REAL(KIND=dp) :: LGrad(3,3), Radius, Density, fa, fb 
     REAL(KIND=dp) :: DSR(3,3),  Em 
     Real(kind=dp) :: Bg, BGlenT, ss, nn
!------------------------------------------------------------------------------
     
      Stress = 0.0
      StrainRate = 0.0
      Spin = 0.0

!
      Density = SUM( NodalDensity(1:n)*Basis(1:n) )
      fa = SUM( Nodalfa(1:n)*Basis(1:n) )
      fb = SUM( Nodalfb(1:n)*Basis(1:n) )
      Wn(1) = SUM( NodalFluidity(1:n)*Basis(1:n) )
      
      Bg=Wn(1)                
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
      END IF
      Em = StrainRate(1,1)+StrainRate(2,2)+StrainRate(3,3)
      DSR = StrainRate
!
! Case non-linear calculate E_D^2 = gamma_e^2/fa + E_m^2/fb
! ----------------------------------------------------------
      ss = 1.0_dp
      IF ( Wn(2) > 1.0 ) THEN

        DO i = 1, 3
          DSR(i,i) = DSR(i,i) - Em/3.0
        END DO

        ss = 0.0_dp
        DO i = 1, 3
          DO j = 1, 3
            ss = ss + DSR(i,j)**2
          END DO
        END DO
        ss = 2.0*ss / fa
        IF ( fb /= 0.0_dp ) ss = ss + Em**2 / fb

        nn =(1.0 - Wn(2))/Wn(2)

        ss = SQRT(ss)
        IF (ss < MinSRInvariant ) ss = MinSRInvariant
        ss =  ss**nn 
      END IF 

!
!    Compute Spin : 
!    --------------

        Spin = 0.5 * ( LGrad - TRANSPOSE(LGrad) )
!
!    Compute deviatoric stresses: 
!    ----------------------------
        Stress = 2.0 * ss * DSR / (fa * Bg**(1.0 / Wn( 2 ) ) ) 
!------------------------------------------------------------------------------
      END SUBROUTINE LocalSD      
!------------------------------------------------------------------------------
!        
!------------------------------------------------------------------------------
      END SUBROUTINE PorousSolver
!------------------------------------------------------------------------------
