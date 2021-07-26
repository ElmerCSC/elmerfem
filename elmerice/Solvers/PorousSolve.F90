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
! * Authors: Juha Ruokolainen, Olivier Gagliardini       
! * Email:   Juha.Ruokolainen@csc.fi                     
! * Web:     http://www.csc.fi/elmer                     
! * Address: CSC - IT Center for Science Ltd.            
! *             Keilaranta 14                               
! *             02101 Espoo, Finland                        
! *                                                         
! *    Original Date: 08 Jun 1997                           
! *
! *****************************************************************************/
!>  Solver for (primarily thermal) Porous material flow
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

     INTEGER :: i, j, k, l, m, n, t, iter, NDeg, STDOFs, LocalNodes, istat
     INTEGER :: dim, comp, nd, nb 

     TYPE(ValueList_t),POINTER :: Material, BC, BodyForce, Constants
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement

     REAL(KIND=dp) :: RelativeChange, UNorm, PrevUNorm, Gravity(3), &
         Normal(3), NonlinearTol, s, Wn(2), MinSRInvariant
         

     REAL(KIND=dp)  :: NodalStresses(3,3), Stress(3,3), &
       NodalStrainRate(3,3),  NodalSpin(3,3) 

     REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
     REAL(KIND=dp), ALLOCATABLE :: SlipCoeff(:,:)
     REAL(KIND=dp) :: u,v,w,detJ
     
     LOGICAL :: stat, CSymmetry 
       
     INTEGER, PARAMETER :: INDi(1:6) = (/ 1, 2, 3, 1, 2, 3 /) ,&
           INDj(1:6)=(/ 1, 2, 3, 2, 3, 1 /)

     INTEGER :: NonlinearIter

     TYPE(Variable_t), POINTER :: PorousSol, DensityVariable
     TYPE(Variable_t), POINTER :: SpinVar
     REAL(KIND=dp), POINTER :: SpinValues(:)
     INTEGER, POINTER :: SpinPerm(:)

     TYPE(Variable_t), POINTER :: DevStressVar, DefHeatVar 
     REAL(KIND=dp), POINTER :: DSValues(:), DefHeatValues(:)
     INTEGER, POINTER :: DSPerm(:), DefHeatPerm(:)

     TYPE(Variable_t), POINTER :: StrainRateVar 
     REAL(KIND=dp), POINTER :: SRValues(:)
     INTEGER, POINTER :: SRPerm(:)

     REAL(KIND=dp), POINTER :: Porous(:), Work(:,:), &
           ForceVector(:),  NodalPorous(:),  &
           DensityValues(:)

     CHARACTER(LEN=MAX_NAME_LEN) :: EquationName, DensityName

     INTEGER, POINTER :: PorousPerm(:), NodeIndexes(:), &
                         DensityPerm(:)

     INTEGER :: PorousType
     LOGICAL :: GotForceBC, GotIt,  &
                NormalTangential=.FALSE.

     INTEGER :: body_id,bf_id
     INTEGER :: old_body = -1
     LOGICAL :: AllocationsDone = .FALSE., FreeSurface, &
                Requal0
           
     REAL(KIND=dp), ALLOCATABLE:: LocalMassMatrix(:,:), &
       LocalStiffMatrix(:,:), LoadVector(:,:), LocalForce(:), &
       Beta(:), & 
       LocalDensity(:), &
       TimeForce(:), RefS(:), RefD(:), RefSpin(:), &
       LocalVelo(:,:), LocalFluidity(:)

#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at, at0
#else
     REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif

!------------------------------------------------------------------------------

     SAVE TimeForce, Basis, dBasisdx

     SAVE LocalMassMatrix, LocalStiffMatrix, LoadVector, &
       LocalForce, ElementNodes, Beta, LocalFluidity, &
       NodalPorous, LocalDensity, Wn, MinSRInvariant, old_body, &
       AllocationsDone

     SAVE RefD, RefS, RefSpin, LocalVelo, SlipCoeff, dim 
              
!------------------------------------------------------------------------------
!  Read constants from constants section of SIF file
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
      IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

      Constants => GetConstants()

      PorousSol => Solver % Variable
      PorousPerm => PorousSol % Perm
      STDOFs =  PorousSol % DOFs
      Porous => PorousSol % Values

      LocalNodes = COUNT( PorousPerm > 0 )
      IF ( LocalNodes <= 0 ) RETURN


      DensityName = GetString(Constants,'Density Name', GotIt)
      IF (.NOT.GotIt) THEN
         CALL WARN('PorousSolve', 'No Keyword >Density Name< defined. Using >Density< as default.')
         WRITE(DensityName,'(A)') 'Density'
      ELSE
         WRITE(Message,'(a,a)') 'Variable Name for density: ', DensityName
         CALL INFO('PorousSolve',Message,Level=12)
      END IF


      DensityVariable => &
              VariableGet(Solver % Mesh %Variables,Densityname)
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

      DefHeatVar => &
               VariableGet(Solver % Mesh % Variables,'Porous Deformational Heat')
      IF ( ASSOCIATED( DefHeatVar ) ) THEN
         DefHeatPerm => DefHeatVar % Perm    
         DefHeatValues => DefHeatVar % Values  
      END IF

      StiffMatrix => Solver % Matrix
      ForceVector => StiffMatrix % RHS
      UNorm = Solver % Variable % Norm

!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
      IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
        N = Solver % Mesh % MaxElementDOFs
        dim = CoordinateSystemDimension()

       IF ( AllocationsDone ) THEN
         DEALLOCATE( ElementNodes % x,     &
                     ElementNodes % y,     &
                     ElementNodes % z,     &
                     LocalVelo,            &
                     Basis, dBasisdx,      &
                     LocalDensity,         &
                     LocalForce,           &
                     RefD, RefS, RefSpin,  &
                     TimeForce,            &
                     LocalMassMatrix,      &
                     LocalStiffMatrix,     &
                     LoadVector, Beta, &
                     SlipCoeff, LocalFluidity)
       END IF

       ALLOCATE( ElementNodes % x( N ), &
                 ElementNodes % y( N ), &
                 ElementNodes % z( N ), &
                 LocalVelo(4,n),&                                     
                 Basis(n), dBasisdx(n,3), &
                 LocalDensity( N ), &
                 LocalForce( 2*STDOFs*N ),&
                 RefD(2*dim*LocalNodes ),&                              
                 RefS(2*dim*LocalNodes ),&                              
                 RefSpin((2*dim-3)*LocalNodes ),&                       
                 TimeForce(2*STDOFs*n), &
                 LocalMassMatrix(2*STDOFs*n,2*STDOFs*n),  &
                 LocalStiffMatrix(2*STDOFs*n,2*STDOFs*n),  &
                 LoadVector(4,n), Beta(n), &
                 SlipCoeff(3,n), LocalFluidity(n), STAT=istat )

       IF ( istat /= 0 ) THEN
          CALL Fatal( 'PorousSolve', 'Memory allocation error.' )
       END IF
!------------------------------------------------------------------------------

       AllocationsDone = .TRUE.
      END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      NonlinearTol = GetConstReal( Solver % Values, &
        'Nonlinear System Convergence Tolerance' )


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
       CALL DefaultInitialize()
!------------------------------------------------------------------------------
       CALL StartAdvanceOutput('PorousSolve','Assembly:')
       DO t=1,Solver % NumberOFActiveElements

         CALL AdvanceOutput(t,GetNOFActive())

         CurrentElement => GetActiveElement(t)

         n = GetElementNOFNodes()
         nd = GetElementNOFDOFs()
         nb = GetElementNOFBDOFs()
         
         CALL GetElementNodes(ElementNodes,CurrentElement)

         Material => GetMaterial()
         IF (.NOT.ASSOCIATED(Material)) THEN
            CALL FATAL('PorousSolve','No Material found')
         END IF

!------------------------------------------------------------------------------
!    Read in material constants from Material section
!------------------------------------------------------------------------------
         body_id = CurrentElement % BodyId
         IF (body_id /= old_body) Then 
            old_body = body_id
            Call  GetMaterialDefs()
         END IF

         LocalFluidity(1:n) = GetReal( Material, &
                         'Fluidity Parameter',  GotIt )
         IF(.NOT.GotIt) THEN
            CALL FATAL('PorousSolve','Variable >Fluidity Parameter< not found.')
         END IF


!------------------------------------------------------------------------------
!        Set body forces
!------------------------------------------------------------------------------
         LoadVector = 0.0D0

         BodyForce => GetBodyForce()
         IF ( ASSOCIATED( BodyForce ) ) THEN
           LoadVector(1,1:n) = LoadVector(1,1:n) + GetReal( &
                   BodyForce, 'Porous Force 1',  gotIt)
           LoadVector(2,1:n) = LoadVector(2,1:n) + GetReal( &
                   BodyForce, 'Porous Force 2',  gotIt)
           LoadVector(3,1:n) = LoadVector(3,1:n) + GetReal( & 
                   BodyForce, 'Porous Force 3',  gotIt)
         END IF
!------------------------------------------------------------------------------
!        Get element local stiffness & mass matrices
!------------------------------------------------------------------------------
         ! The Density can be a DG variable and it is then safe to call it 
         ! using the permutation vector
         IF (DensityVariable%TYPE == Variable_on_nodes_on_elements) THEN
           LocalDensity(1:n) =  &
              DensityValues(DensityPerm(CurrentElement % DGIndexes(1:n)))
         ELSE
           LocalDensity(1:n) =  &
              DensityValues(DensityPerm(CurrentElement % NodeIndexes(1:n))) 
         ENDIF
         
         CALL GetVectorLocalSolution(LocalVelo)
          

         CALL LocalMatrix( LocalMassMatrix, LocalStiffMatrix, &
              LocalForce, LoadVector, LocalDensity, LocalVelo, &
              LocalFluidity, CurrentElement, n, nd, nd+nb, &
              ElementNodes, Wn, MinSRInvariant)

!------------------------------------------------------------------------------
!        Condensate in case of bubbles elements       
!------------------------------------------------------------------------------
         IF (nb > 0) THEN
             TimeForce = 0.0d0
             CALL NSCondensate( nd, nb, dim, LocalStiffMatrix, LocalForce, TimeForce )
         END IF


!------------------------------------------------------------------------------
!        Update global matrices from local matrices 
!------------------------------------------------------------------------------
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

        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()

        CALL GetElementNodes(ElementNodes,CurrentElement)

        BC => GetBC()
        IF ( ASSOCIATED( BC ) ) THEN
            LoadVector = 0.0D0
            Beta       = 0.0D0
!------------------------------------------------------------------------------
!           Force in given direction BC: \tau\cdot n = F
!------------------------------------------------------------------------------
            GotForceBC = .FALSE.

            LoadVector(1,1:n) = &
                    GetReal( BC, 'Force 1',   GotIt )
            GotForceBC = GotForceBC .OR. gotIt

            LoadVector(2,1:n) = & 
                     GetReal( BC, 'Force 2',  GotIt )
            GotForceBC = GotForceBC .OR. gotIt

            LoadVector(3,1:n) = &
                     GetReal( BC, 'Force 3',  GotIt )
            GotForceBC = GotForceBC .OR. gotIt

            Beta(1:n) = &
                GetReal( BC, 'Normal Force',  GotIt )
            GotForceBC = GotForceBC .OR. gotIt

!------------------------------------------------------------------------------
!             slip boundary condition BC: \tau\cdot n = R_k u_k
!------------------------------------------------------------------------------

              SlipCoeff = 0.0d0
              SlipCoeff(1,1:n) =  GetReal( BC, &
                    'Porous Slip Coeff 1',  GotIt )
              GotForceBC = GotForceBC .OR. gotIt

              SlipCoeff(2,1:n) =  GetReal( BC, &
                    'Porous Slip Coeff 2',  GotIt )
              GotForceBC = GotForceBC .OR. gotIt

              SlipCoeff(3,1:n) =  GetReal( BC, &
                    'Porous Slip Coeff 3',  GotIt )
              GotForceBC = GotForceBC .OR. gotIt

              IF ( .NOT.GotForceBC ) CYCLE

              NormalTangential = GetLogical( BC, &
                 'Normal-Tangential Porous', GotIt )

              IF (.NOT.GotIt) THEN
                NormalTangential = GetLogical( BC, &
                       'Normal-Tangential '//GetVarName(Solver % Variable), GotIt )
              END IF
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
            CALL LocalMatrixBoundary( LocalStiffMatrix, LocalForce, &
                 LoadVector, Beta, SlipCoeff, NormalTangential, &
                 CurrentElement, n, nd, ElementNodes )
!------------------------------------------------------------------------------

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

      CALL Info( 'PorousSolve', 'Set boundaries done', Level=4 )

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

      WRITE( Message, '(A,e12.4,e12.4)' ) 'Result Norm   : ',UNorm, PrevUNorm
      CALL Info( 'PorousSolve', Message, Level=4 )
      WRITE( Message, '(A,e12.4)' ) 'Relative Change : ',RelativeChange
      CALL Info( 'PorousSolve', Message, Level=4 )

!------------------------------------------------------------------------------

      IF ( RelativeChange < NonLinearTol ) EXIT

!------------------------------------------------------------------------------
    END DO ! of nonlinear iter
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!   Compute the StrainRate, Spin  and deviatoric Stress
!   Nodal values      
!------------------------------------------------------------------------------

     IF (ASSOCIATED( StrainRateVar).OR.ASSOCIATED(DevStressVar)&
      .OR.ASSOCIATED(SpinVar).OR.ASSOCIATED(DefHeatVar)) THEN
       RefD=0.0_dp
       RefS=0.0_dp
       RefSpin=0.0_dp
       IF (ASSOCIATED(StrainRateVar)) SRValues = 0.0_dp
       IF (ASSOCIATED(devStressVar)) DSValues = 0.0_dp
       IF (ASSOCIATED(SPinVar)) SpinValues = 0.0_dp

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

        LocalFluidity(1:n) = GetReal( Material, &
                      'Fluidity Parameter',  GotIt )
        IF (.NOT.GotIt) THEN
           CALL FATAL('PorousSolve', 'Variable Fluidity Parameter not found')
        END IF

         CALL GetElementNodes(ElementNodes,CurrentElement)

         CALL GetScalarLocalSolution(LocalDensity,DensityName)

         CALL GetVectorLocalSolution(LocalVelo)

! Go for all nodes of the element        
         Do i=1,n

! u, v, w local coord of node i
           u = CurrentElement % Type % NodeU(i)
           v = CurrentElement % Type % NodeV(i)
           w = CurrentElement % Type % NodeW(i)
       
           stat = ElementInfo( CurrentElement, ElementNodes, u, v, w, detJ, &
                Basis, dBasisdx )
! Axi symmetric case when R=0 strain, stress not calculated exactly in
! x=0 (I agree it is not very nice, better solution ???)
        Requal0 = .False.
        IF (( CSymmetry) .And. & 
          (SUM(ElementNodes % x(1:n) * Basis(1:n)) == 0.0)) THEN  
           Requal0 = .True.
            u= u + 0.0001  
            stat = ElementInfo( CurrentElement, ElementNodes, u, v, w, detJ, &
                Basis, dBasisdx )
        END IF

        CALL LocalSD(NodalStresses, NodalStrainRate, NodalSpin, & 
                 LocalVelo, LocalFluidity,  &
                LocalDensity, CSymmetry, Basis, dBasisdx, &
                CurrentElement, n, ElementNodes, dim, Wn, &
                MinSRInvariant)
                

        IF (Requal0) NodalSpin = 0.0_dp 


           IF (ASSOCIATED(StrainRateVar)) THEN
             RefD(2*dim*(SRPerm(NodeIndexes(i))-1)+1 : &
                                      2*dim*SRPerm(NodeIndexes(i))) &
                 =RefD(2*dim*(SRPerm(NodeIndexes(i))-1)+1 : &
                                      2*dim*SRPerm(NodeIndexes(i))) + 1.
             comp=0
             DO j=1,2*dim
               comp=comp+1
               SRValues(2*dim*(SRPerm(NodeIndexes(i))-1)+comp)=&
               SRValues(2*dim*(SRPerm(NodeIndexes(i))-1)+comp) + &
                NodalStrainRate(INDi(j),INDj(j))
             END DO
           END IF

           IF (ASSOCIATED(DevStressVar)) THEN
             RefS(2*dim*(DSPerm(NodeIndexes(i))-1)+1 : &
                                      2*dim*DSPerm(NodeIndexes(i))) &
                  =RefS(2*dim*(DSPerm(NodeIndexes(i))-1)+1 :  &
                                      2*dim*DSPerm(NodeIndexes(i))) + 1.
             comp=0
             DO j=1,2*dim
               comp=comp+1
               DSValues(2*dim*(DSPerm(NodeIndexes(i))-1)+comp)=&
                DSValues(2*dim*(DSPerm(NodeIndexes(i))-1)+comp) + &
                NodalStresses(INDi(j),INDj(j))
             END DO
           END IF

           IF (ASSOCIATED(SpinVar)) THEN
            RefSpin((2*dim-3)*(SpinPerm(NodeIndexes(i))-1)+1 :  &
                                (2*dim-3)*SpinPerm(NodeIndexes(i))) &
                =RefSpin((2*dim-3)*(SpinPerm(NodeIndexes(i))-1)+1 :  &
                                (2*dim-3)*SpinPerm(NodeIndexes(i))) + 1.
             comp=0
             DO j=1,(2*dim-3)
             comp=comp+1
             SpinValues((2*dim-3)*(SpinPerm(NodeIndexes(i))-1)+comp)=&
             SPinValues((2*dim-3)*(SpinPerm(NodeIndexes(i))-1)+comp) + &
             NodalSpin(INDi(j+3),INDj(j+3))
             END DO
           END IF

           IF (ASSOCIATED(defHeatVar)) THEN
              k = DefHeatPerm(NodeIndexes(i))
              defHeatValues(k) = 0.0_dp
              Stress = NodalStresses
              DO l=1,DIM
                 Stress(l,l) = Stress(l,l) - LocalVelo(DIM+1,i)  
                 DO m=1,DIM
                    defHeatValues(k) = defHeatValues(k) +  &
                            Stress(l,m)*NodalStrainRate(l,m)
                 END DO
              END DO
              defHeatValues(k) = Max(defHeatValues(k),0.0)
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

      Wn(2) = GetConstReal( Material , 'Powerlaw Exponent', GotIt )
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
      MinSRInvariant = 100.0*AEPS
      IF ( Wn(2) > 1.0 ) THEN
        MinSRInvariant =  &
             GetConstReal( Material, 'Min Second Invariant', GotIt )
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
              LoadVector, NodalDensity, NodalVelo,  &
              NodalFluidity, Element, n, nd, ntot, Nodes, Wn, MinSRInvariant)
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: StiffMatrix(:,:), MassMatrix(:,:)
     REAL(KIND=dp) :: LoadVector(:,:), NodalVelo(:,:)
     REAL(KIND=dp) :: Wn(2), MinSRInvariant
     REAL(KIND=dp), DIMENSION(:) :: ForceVector, NodalDensity,  &
                    NodalFluidity
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
     INTEGER :: n, nd, ntot
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(2*ntot)
     REAL(KIND=dp) :: dBasisdx(2*ntot,3),detJ

     REAL(KIND=dp) :: Force(3), density, fa, fb 

     REAL(KIND=dp), DIMENSION(4,4) :: A,M
     REAL(KIND=dp) :: Load(3)
     REAL(KIND=dp) :: nn, ss, LGrad(3,3), SR(3,3)

     INTEGER :: i, j, k, p, q, t, c, dim, NBasis, ind(3)

     REAL(KIND=dp) :: s,u,v,w, Radius
  
     REAL(KIND=dp) :: Em, eta, Kcp
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp) :: ParameterA, ParameterB
     INTEGER :: N_Integ

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat, CSymmetry


!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()

      ForceVector = 0.0_dp
      StiffMatrix = 0.0_dp
      MassMatrix  = 0.0_dp
!    
!    Integration stuff
!    
      
      NBasis = ntot
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
      stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
            Basis,dBasisdx )

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
!    variables at the integration point
!    Interpolate density first and the evaluate parameters
!    Density is the Relative Density
      Density = SUM( NodalDensity(1:n)*Basis(1:n) )
      fa = ParameterA(Density)
      fb = ParameterB(Density)


      Wn(1) = SUM( NodalFluidity(1:n)*Basis(1:n) )

      CSymmetry = CurrentCoordinateSystem() == AxisSymmetric
      IF ( CSymmetry ) s = s * Radius

!
! Case non-linear calculate E_D^2 = gamma_e^2/fa + E_m^2/fb
! ----------------------------------------------------------
      ss = 1.0_dp
      LGrad = 0.0_dp
      IF ( Wn(2) > 1.0 ) THEN
        LGrad(1:dim,1:dim) = MATMUL( NodalVelo(1:dim,1:nd), dBasisdx(1:nd,1:dim) )
        SR = 0.5 * ( LGrad + TRANSPOSE(LGrad) )
      
        IF ( CSymmetry ) THEN
          SR(1,3) = 0.0_dp
          SR(2,3) = 0.0_dp
          SR(3,1) = 0.0_dp
          SR(3,2) = 0.0_dp
          SR(3,3) = 0.0_dp
          IF ( Radius > 10*AEPS ) THEN
            SR(3,3) = SUM( Nodalvelo(1,1:nd) * Basis(1:nd) ) /Radius
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
        IF ( fb > 1.0e-8 ) ss = ss + Em**2 / fb

        nn =(1.0 - Wn(2))/Wn(2)

        ss = SQRT(ss)
        IF (ss < MinSRInvariant ) ss = MinSRInvariant
        ss =  ss**nn 
      END IF
!
! Bulk effective viscosity and Compressibility parameter Kcp 
!     
      eta =  ss / (fa * Wn( 1 )**( 1.0 / Wn( 2 ) ) ) 
      Kcp = fb * Wn( 1 )**( 1.0 / Wn( 2 ) ) / ss 

!
!    Loop over basis functions (of both unknowns and weights)
!

      DO p=1,NBasis
       DO q=1,NBasis

        A = 0.0_dp
        M = 0.0_dp

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

      c = dim + 1
      DO i=n+1,ntot
        StiffMatrix( c*i, : ) = 0._dp
        StiffMatrix( :, c*i ) = 0._dp
        ForceVector( c*i ) = 0._dp
        StiffMatrix( c*i, c*i ) = 1._dp
      END DO
      
!------------------------------------------------------------------------------
      END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE LocalMatrixBoundary( BoundaryMatrix, BoundaryVector, &
                 LoadVector, NodalBeta, NodalSlipCoeff, & 
                  NormalTangential, Element, n, nd, Nodes )
                      
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: BoundaryMatrix(:,:), BoundaryVector(:)
     REAL(KIND=dp) :: NodalBeta(:), LoadVector(:,:)
     REAL(KIND=dp) :: NodalSlipCoeff(:,:)
     TYPE(Element_t),POINTER  :: Element
     TYPE(Nodes_t)    :: Nodes
     LOGICAL :: NormalTangential
     INTEGER :: n,nd
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(nd)
     REAL(KIND=dp) :: dBasisdx(nd,3),detJ

     REAL(KIND=dp) :: u,v,w,s
     REAL(KIND=dp) :: Force(3), Beta, Normal(3)
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
                Basis, dBasisdx )

       s = detJ * S_Integ(t)
       IF ( CurrentCoordinateSystem() == AxisSymmetric ) &
        s = s * SUM( Nodes % x(1:n) * Basis(1:n) )
!------------------------------------------------------------------------------
       Force = 0.0D0
       DO i=1,dim
         Force(i) = SUM( LoadVector(i,1:n)*Basis(1:n) )
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
         DO p=1,nd
           DO q=1,nd
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

      DO q=1,nd
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
        NodalVelo, NodalFluidity, Nodaldensity,  &
        CSymmetry, Basis, dBasisdx, Element, n,  Nodes, dim,  Wn, &
        MinSRInvariant)
       
!------------------------------------------------------------------------------
!    Subroutine to compute the nodal Strain-Rate, Stress, ...
!------------------------------------------------------------------------------
     LOGICAL ::  CSymmetry 
     INTEGER :: n, dim
     REAL(KIND=dp) :: Stress(:,:), StrainRate(:,:), Spin(:,:)
     REAL(KIND=dp) :: NodalVelo(:,:), NodalFluidity(:)
     REAL(KIND=dp) :: Basis(:)
     REAL(KIND=dp) :: dBasisdx(:,:)
     REAL(KIND=dp) :: detJ
     REAL(KIND=dp) :: NodalDensity(:)
     REAL(KIND=dp) :: u, v, w      
     REAL(KIND=dp) :: Wn(2),  D(6), MinSRInvariant
     LOGICAL :: Isotropic
     REAL(KIND=dp) :: ParameterA, ParameterB
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
     LOGICAL :: stat
     INTEGER :: i,j,k,p,q
     REAL(KIND=dp) :: LGrad(3,3), Radius, Density, fa, fb 
     REAL(KIND=dp) :: DSR(3,3),  Em 
     Real(kind=dp) :: ss, nn
!------------------------------------------------------------------------------
     
      Stress = 0.0
      StrainRate = 0.0
      Spin = 0.0

!
      Density = SUM( NodalDensity(1:n)*Basis(1:n) )
      fa = ParameterA(Density)
      fb = ParameterB(Density)

      Wn(1) = SUM( NodalFluidity(1:n)*Basis(1:n) )
      
!
!    Compute strainRate : 
!    -------------------

      LGrad = 0.0_dp
      LGrad(1:dim,1:dim) = MATMUL( NodalVelo(1:dim,1:n), dBasisdx(1:n,1:dim) )
        
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
        IF ( fb > 1.0e-8 ) ss = ss + Em**2 / fb

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
        Stress = 2.0 * ss * DSR / (fa * Wn( 1 )**(1.0 / Wn( 2 ) ) ) 
!------------------------------------------------------------------------------
      END SUBROUTINE LocalSD      
!------------------------------------------------------------------------------
!        
!------------------------------------------------------------------------------
      END SUBROUTINE PorousSolver
!------------------------------------------------------------------------------
