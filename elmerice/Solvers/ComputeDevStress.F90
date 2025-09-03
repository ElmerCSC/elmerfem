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
! *  Authors: Juha Ruokolainen, Fabien Gillet-Chaulet, Olivier Gagliardini
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 08 Jun 1997
! *  Date of modification: 16/01/2024
! *  Last modification by Julien Brondex (IGE) to make it compatible with 
! *  Porous Solver
! *****************************************************************************
!> Module containing a solver for computing deviatoric or Cauchy   
!>          stress from flow solution (NS or Porous solver)    
!> 2D SDOFs = 4 (S11, S22, S33, S12)                               
!> 3D SDOFs = 6 (S11, S22, S33, S12, S23, S31)                     
!> Keywords : Cauchy (Logical),                                    
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
  
  INTEGER :: i, j, k, l, n, t, iter, NDeg, iSolverName
  INTEGER :: dim, STDOFs, StressDOFs, LocalNodes, istat
  
  TYPE(ValueList_t),POINTER :: Material, BC, BodyForce, Constants
  TYPE(Nodes_t) :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement
  
  REAL(KIND=dp) :: RelativeChange, UNorm, PrevUNorm
  
  REAL(KIND=dp), ALLOCATABLE :: Basis(:),ddBasisddx(:,:,:)
  REAL(KIND=dp), ALLOCATABLE :: dBasisdx(:,:)
  REAL(KIND=dp) :: u, v, w, detJ
  
  LOGICAL :: stat, CSymmetry 
  
  INTEGER :: NewtonIter, NonlinearIter
  
  TYPE(Variable_t), POINTER :: StressSol, FlowVariable, DensityVariable, VMisesVariable
  
  REAL(KIND=dp), POINTER ::  Stress(:), Solution(:), &
       ForceVector(:), FlowValues(:), DensityValues(:), VMisesValues(:)
  
  INTEGER, POINTER :: StressPerm(:), NodeIndexes(:), &
       FlowPerm(:), DensityPerm(:), VMisesPerm(:)
  
  LOGICAL :: Isotropic, AllocationsDone = .FALSE.,  &
       Requal0, ComputeVMises
  LOGICAL :: GotIt, GotIt_CT, GotIt_TFV, OldKeyword,  Cauchy = .FALSE.,UnFoundFatal=.TRUE.,OutOfPlaneFlow
  
  REAL(KIND=dp), ALLOCATABLE:: LocalMassMatrix(:,:), &
       LocalStiffMatrix(:,:), LocalForce(:), &
       LocalP(:),  &
       LocalVelo(:,:), LocalViscosity(:), &
       LocalFluidity(:), LocalDensity(:), &
       RelativeT(:), ConstantT(:)
  
  INTEGER :: NumberOfBoundaryNodes, COMP
  INTEGER, POINTER :: BoundaryReorder(:)
  
  REAL(KIND=dp), POINTER :: BoundaryNormals(:,:), &
       BoundaryTangent1(:,:), BoundaryTangent2(:,:)
  CHARACTER(LEN=MAX_NAME_LEN) :: FlowSolverName, StressSolverName, TempName, DensityName, VMisesName
  
  REAL(KIND=dp) :: at, at0
!------------------------------------------------------------------------------
  SAVE NumberOfBoundaryNodes, BoundaryReorder, BoundaryNormals, &
       BoundaryTangent1, BoundaryTangent2
  
  SAVE Basis, dBasisdx, ddBasisddx
  SAVE LocalMassMatrix, LocalStiffMatrix, LocalForce, &
       ElementNodes,  &
       AllocationsDone,  &
       LocalViscosity, Cauchy, &
       LocalFluidity, LocalDensity, &
       OldKeyword, RelativeT, ConstantT
  
  SAVE LocalVelo, LocalP, dim
  
!  NULLIFY(StressSol, FlowVariable)

!------------------------------------------------------------------------------
!  Read the name of the Flow Solver (NS or Porous)
!------------------------------------------------------------------------------
  FlowSolverName = GetString( Solver % Values, 'Flow Solver Name', GotIt )    
  IF (.NOT.Gotit) THEN
     CALL FATAL('ComputeDevStress', '>Flow Solver Name< must be prescribed in Solver &
                 and set to "Flow Solution" or "Porous".') 
  ELSE
     FlowVariable => VariableGet( Solver % Mesh % Variables, FlowSolverName)
     IF ( ASSOCIATED( FlowVariable ) ) THEN
        FlowPerm    => FlowVariable % Perm
        FlowValues  => FlowVariable % Values
        OutOfPlaneFlow = GetLogical(Solver % Values , 'Out of Plane flow', GotIt)
        IF ( .NOT. GotIt ) OutOfPlaneFlow = .FALSE.
     ELSE
        CALL FATAL('ComputeDevStress', 'Flow Variable not associated. Check consistency between used flow solver &
                    and prescribed >Flow Solver Name< !')
     END IF
  END IF

  !!!Switch to integer indice to avoid case sensitivity issues
  SELECT CASE(FlowSolverName)
     CASE('Flow Solution','flow solution')
        iSolverName = 1
     CASE('Porous', 'porous')
        iSolverName = 2
     CASE DEFAULT
        CALL FATAL('ComputeDevStress', '>Flow Solver Name< must be either "Flow Solution" or "Porous"')
  END SELECT 
!------------------------------------------------------------------------------
!  Read the name of the Density Variable when using Porous solver
!------------------------------------------------------------------------------
  IF (iSolverName == 2) THEN
     Constants => GetConstants()
     DensityName = GetString(Constants,'Density Name', GotIt)
     IF (.NOT.GotIt) THEN
        CALL WARN('ComputeDevStress', 'No Keyword >Density Name< defined despite using Porous Solver.&
                   Using "Density" as default.')
        WRITE(DensityName,'(A)') 'Density'
     ELSE
        WRITE(Message,'(a,a)') 'Variable Name for density: ', DensityName
        CALL INFO('ComputeDevStress',Message,Level=12)
     END IF

     DensityVariable => VariableGet(Solver % Mesh %Variables, Densityname)
     IF ( ASSOCIATED( DensityVariable ) ) THEN
        DensityPerm    => DensityVariable % Perm
        DensityValues => DensityVariable % Values
     ELSE
        CALL FATAL('ComputeDevStress', 'Density not associated. Required as viscosity=f(D) for Porous !')
     END IF
  END IF

  ComputeVMises = GetLogical( Solver % Values, 'Compute von Mises Stress', GotIt )
  IF (.NOT.GotIt) THEN
    ComputeVMises = .FALSE.
  ELSE IF (ComputeVMises) THEN
    CALL INFO('ComputeDevStress', 'Computing von Mises stress', Level=5)
    VMisesName = GetString(Solver % Values,'Von Mises Stress Name', GotIt)
    IF (.NOT.GotIT) THEN
      WRITE(VMisesName,'(A)') 'Von Mises Stress'
    END IF
    VMisesVariable => VariableGet(Solver % Mesh % Variables, VMisesname)
    IF ( ASSOCIATED( VMisesVariable ) ) THEN
        VMisesPerm    => VMisesVariable % Perm
        VMisesValues => VMisesVariable % Values
        CALL INFO('ComputeDevStress', 'Output of von Mises stress to ' // TRIM(VMisesname), Level=3)
     ELSE
        CALL FATAL('ComputeDevStress', TRIM(VMisesName) // ' not associated, but output requested')
     END IF
  ELSE
    ComputeVMises = .FALSE.
  END IF
  
!------------------------------------------------------------------------------
!  Read constants from constants section of SIF file
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  
  Solution => Solver % Variable % Values
  STDOFs   =  Solver % Variable % DOFs
  
  IF ( STDOFs /=1 ) THEN
     CALL Fatal( 'ComputeDevStress', 'DOF must be equal to 1' )
  END IF
  
  StressSolverName = GetString( Solver % Values, 'Stress Variable Name', GotIt )    
  IF (.NOT.Gotit) CALL FATAL('ComputeDevStress', & 
       'Stress Variable Name not defined')
  
  StressSol => VariableGet( Solver % Mesh % Variables, TRIM(StressSolverName),&
       UnFoundFatal=UnFoundFatal )
  StressPerm => StressSol % Perm
  StressDOFs = StressSol % DOFs
  Stress => StressSol % Values
  
  dim = CoordinateSystemDimension()
  IF (StressDOfs /= 2*dim) THEN
     CALL Fatal( 'ComputeDesStress', 'Bad dimension of Stress Variable (4 in 2D, 6 in 3D)' )
  ENDIF
  
  StiffMatrix => Solver % Matrix
  ForceVector => StiffMatrix % RHS
  Unorm = SQRT( SUM( Stress**2 ) / SIZE(Stress) )
!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone .OR. Solver % MeshChanged) THEN
     N = Model % MaxElementNodes
     
     IF ( AllocationsDone ) THEN
        DEALLOCATE( ElementNodes % x,     &
             ElementNodes % y,     &
             ElementNodes % z,     &
             LocalVelo, LocalP,    &                      
             Basis, ddBasisddx,    &
             dBasisdx,             &
             LocalMassMatrix,      &
             LocalStiffMatrix,     &
             LocalForce,           &
             LocalViscosity,       &
             LocalFluidity,        &
             LocalDensity,         &
             RelativeT,            &
             ConstantT )
     END IF

     ALLOCATE( ElementNodes % x( N ), &
          ElementNodes % y( N ), &
          ElementNodes % z( N ), &
          LocalVelo( 3,N ), LocalP( N ), &                                     
          Basis( 2*N ),ddBasisddx(1,1,1), dBasisdx( 2*N,3 ), &
          LocalMassMatrix( 2*STDOFs*N,2*STDOFs*N ),  &
          LocalStiffMatrix( 2*STDOFs*N,2*STDOFs*N ),  &
          LocalForce( 2*STDOFs*N ),  &
          LocalViscosity(N), &
          LocalFluidity(N), LocalDensity(N), &
          RelativeT(N), ConstantT(N), STAT=istat )

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

! Loop over the Stress components [Sxx, Syy, Szz, Sxy, Syz, Szx] 

     PrevUNorm = UNorm

     DO COMP = 1, 2*dim

        WRITE(Message,'(a,i3)' ) ' Component : ', COMP  
        CALL Info( 'ComputeDevStress', Message, Level=5 )


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
           !!!! When the Flow Solver is Stokes:
           IF (iSolverName == 1) THEN       
              !! Check for the presence of keywords related to the new vectorized Stokes solver
              RelativeT(1:n) = ListGetReal( Material, 'Relative Temperature', n, NodeIndexes, GotIt )
              IF (GotIt) THEN
                 OldKeyword = ListGetLogical( Material, 'Glen Allow Old Keywords', GotIt)
                 IF (.NOT.(GotIt .AND. OldKeyword)) THEN
                    CALL FATAL('ComputeDevStress', 'When using ComputeDevStress with IncompressibleNSVec &
                         >Glen Allow Old Keywords< must be set to True in Material')
                 END IF
                 ConstantT(1:n) = ListGetReal( Material, 'Constant Temperature', n, NodeIndexes, GotIt_CT )
                 TempName = GetString( Material, 'Temperature Field Variable', GotIt_TFV)
                 IF (.NOT.(GotIt_CT .OR. GotIt_TFV)) THEN
                    CALL FATAL('ComputeDevStress', '>Constant Temperature< or >Temperature Field Variable< &
                            must be prescribed in Material')
                 END IF
                 !!! In the case of constant T check for consistency between prescribed relative and constant T
                 IF (GotIt_CT .AND. (.NOT. all(abs(RelativeT(1:n) - ConstantT(1:n))<AEPS))) THEN
                    CALL FATAL('ComputeDevStress', 'When considering constant temperature, >Constant Temperature< and >Relative &
                    Temperature< must be consistent')
                 END IF
              END IF
              !Get Viscosity at nodes
              LocalViscosity(1:n) = GetReal( Material, 'Viscosity', GotIt)
              IF(.NOT.GotIt) CALL FATAL('ComputeDevStress','Variable >Viscosity Parameter< not found.')
              !Give dummy default nodal values for rheology parameters associated to Porous: 
              LocalFluidity(1:n) = 1.0_dp
              LocalDensity(1:n) = 1.0_dp

             !!!! When the Flow Solver is Porous:
           ELSE IF (iSolverName ==2) THEN
              ! Get fluidity at nodes
              LocalFluidity(1:n) = GetReal( Material, 'Fluidity Parameter',  GotIt )
              IF(.NOT.GotIt) CALL FATAL('ComputeDevStress','Variable >Fluidity Parameter< not found.')
              ! Get Density at element nodes
              ! The Density can be a DG variable and it is then safe to call it 
              ! using the permutation vector
              IF (DensityVariable%TYPE == Variable_on_nodes_on_elements) THEN
                 LocalDensity(1:n) =  DensityValues(DensityPerm(CurrentElement % DGIndexes(1:n)))
              ELSE
                 LocalDensity(1:n) = DensityValues(DensityPerm(CurrentElement % NodeIndexes(1:n))) 
              END IF
              !Give dummy default nodal values for rheology parameters associated to Stokes: 
              LocalViscosity(1:n) = 1.0_dp
           END IF
!-------------------------------------------------------------------------------------
!Independently of Flow Solver, get nodal velo, nodal pressure and cauchy stress option
!-------------------------------------------------------------------------------------

! Do we want to return Cauchy or deviatoric stresses ? Deviatoric by default
	   Cauchy = ListGetLogical( Material , 'Cauchy', Gotit )
           IF (.NOT.Gotit) THEN
              Cauchy = .FALSE.
              WRITE(Message,'(A)') 'Cauchy set to False'
              CALL INFO('ComputeDevStress', Message, Level = 20)
           END IF

           LocalVelo = 0.0_dp
           DO i=1, dim
             LocalVelo(i,1:n) = FlowValues((dim+1)*(FlowPerm(NodeIndexes(1:n))-1) + i)
           END DO
           BodyForce => GetBodyForce()
           IF ( dim < 3 .AND. OutOfPlaneFlow ) THEN
             LocalVelo(DIM+1,1:n) = ListGetReal(BodyForce,'Out Of Plane Velocity',&
                  n, NodeIndexes(1:n),GotIt)
             IF (.NOT.GotIt) &
                  CALL WARN('ComputeDevStress',"Out of plane velocity not found")
           END IF
           
           ! Get Pressure at element nodes, i.e FlowValues((dim+1)*(FlowPerm(NodeIndexes(1:n))-1) + (dim +1))
           LocalP(1:n) = FlowValues((dim+1)*FlowPerm(NodeIndexes(1:n)))


           CALL LocalMatrix(COMP, LocalMassMatrix, LocalStiffMatrix, &
                LocalForce,  LocalVelo, LocalP, &
                LocalViscosity, LocalFluidity, LocalDensity, &
                CurrentElement, n, &
                ElementNodes, Cauchy, iSolverName)
            
!------------------------------------------------------------------------------
!        Update global matrices from local matrices 
!------------------------------------------------------------------------------
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

        DO t=1,Solver % NumberOfActiveElements
           CurrentElement => GetActiveElement(t) 
           n = GetElementNOFNodes()
           DO i=1,n
              k = CurrentElement % NodeIndexes(i)
              Stress( StressDOFs*(StressPerm(k)-1) + COMP ) =    & 
                   Solver % Variable % Values( Solver % Variable % Perm(k) )
           END DO
        END DO

     END DO ! End DO Comp


     IF (ComputeVMises) THEN
       VMisesValues = 0.0_dp
       DO k=1,Solver % Mesh % NumberOfNodes
         IF (VMisesPerm(k) > 0) THEN
           IF (DIM == 3) THEN
             DO COMP = 1, 3
               i = COMP + 1
               IF (COMP == 3) i = 1
               VMisesValues(VMisesPerm(k)) =  VMisesValues(VMisesPerm(k)) &
                    + 0.5_dp *( Stress( StressDOFs*(StressPerm(k)-1) + COMP ) &
                    - Stress( StressDOFs*(StressPerm(k)-1) + i ) )**2.0_dp &
                    + 3.0_dp * (Stress( StressDOFs*(StressPerm(k)-1) + COMP + 3))**2.0_dp
             END DO
           ELSE IF (DIM == 2) THEN
             DO COMP = 1, 2
               VMisesValues(VMisesPerm(k)) =  VMisesValues(VMisesPerm(k)) &
                    + (Stress( StressDOFs*(StressPerm(k)-1) + COMP))**2.0_dp
             END DO
             VMisesValues(VMisesPerm(k)) =  VMisesValues(VMisesPerm(k))  &
                  - (Stress( StressDOFs*(StressPerm(k)-1) + 1)) &
                  * (Stress( StressDOFs*(StressPerm(k)-1) + 2)) &
                  + 3.0_dp * (Stress( StressDOFs*(StressPerm(k)-1) + 4))**2.0_dp
           END IF                  
           VMisesValues(VMisesPerm(k)) = SQRT(VMisesValues(VMisesPerm(k)))
         END IF
       END DO
     END IF
     
     Unorm = SQRT( SUM( Stress**2 ) / SIZE(Stress) )
     Solver % Variable % Norm = Unorm  

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


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(COMP, MassMatrix, StiffMatrix, ForceVector, &
       NodalVelo, NodalP, NodalViscosity, NodalFluidity, NodalDensity, &
       Element, n, Nodes, Cauchy, iSolverName )
!------------------------------------------------------------------------------
    
    USE MaterialModels
    USE PorousMaterialModels
    
    REAL(KIND=dp) :: StiffMatrix(:,:), MassMatrix(:,:)
    REAL(KIND=dp) ::  NodalVelo(:,:)
    REAL(KIND=dp), DIMENSION(:) :: ForceVector,  &
         NodalViscosity, NodalP, NodalFluidity, NodalDensity
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    LOGICAL ::  Cauchy
    INTEGER :: n, COMP
    INTEGER :: iSolverName
!------------------------------------------------------------------------------
!
    REAL(KIND=dp) :: Basis(2*n),ddBasisddx(1,1,1)
    REAL(KIND=dp) :: dBasisdx(2*n,3),detJ, pBasis(n)
    
    REAL(KIND=dp) :: Stress, epsi
    
    REAL(KIND=dp) :: Pressure
    REAL(KIND=dp) :: LGrad(3,3), SR(3,3)
    
    INTEGER :: i, j, k, p, q, t, dim, cc, NBasis,  LinearBasis
    
    REAL(KIND=dp) :: s, u, v, w, Radius
    
    REAL(KIND=dp) :: Viscosity, Fluidity, Density, mu(2)
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    
    INTEGER :: N_Integ, nd
    INTEGER, DIMENSION(6), PARAMETER :: indx = (/1, 2, 3, 1, 2, 3/), &
                                         indy = (/1, 2, 3, 2, 3, 1/)
    
    REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
    
    LOGICAL :: stat, CSymmetry
    
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
! Deviatoric Strain-Rate at IP
!
       LGrad = 0.0_dp
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
             SR(i,i) = SR(i,i) - epsi/dim !Deviatoric SR
          END DO
       END IF
!
! Get Effective Viscosity at IP
!
       !!! For the Stokes
       IF (iSolverName == 1) THEN
          Viscosity = SUM( NodalViscosity(1:n)*Basis(1:n) )
          Viscosity = EffectiveViscosity( Viscosity, 1.0_dp, NodalVelo(1,1:n), NodalVelo(2,1:n), NodalVelo(3,1:n), &
                  Element, Nodes, n, n, u, v, w, LocalIP=t )
       !!! For the Porous
       ELSE IF (iSolverName == 2) THEN
          Fluidity = SUM( NodalFluidity(1:n)*Basis(1:n) )
          Density = SUM( NodalDensity(1:n)*Basis(1:n) )
          mu = PorousEffectiveViscosity( Fluidity, Density, SR, epsi, &
                  Element, Nodes, n, n, u, v, w, LocalIP=t )  !!mu=[eta, Kcp]
          Viscosity = mu(1)
       END IF
       
!
!    Compute deviatoric stresses or Cauchy stresses at IP for current COMP: 
!    ----------------------------
      
       Stress = 2.0 * Viscosity * SR(indx(COMP),indy(COMP))
       
       IF ((Cauchy).AND.(COMP.LE.3)) THEN
          Pressure = SUM( NodalP(1:n)*Basis(1:n) )
          Stress = Stress - Pressure
       END IF
       
       DO p=1,n         
          DO q=1,n        
             StiffMatrix(p,q) =  &
                  StiffMatrix(p,q) + s*Basis(q)*Basis(p)
          END DO
          ForceVector(p) =  &
               ForceVector(p) + s*Stress*Basis(p) 
       END DO
       
    END DO
     
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE ComputeDevStress
!------------------------------------------------------------------------------
