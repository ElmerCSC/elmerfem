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
! *  include the mass conservation equation d D / dt + div (D.u) = 0
! *  variables u, v, (w), p, D
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
   RECURSIVE SUBROUTINE LumiSolver( Model,Solver,dt,TransientSimulation )
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
!        and d D/dt + div(D.u)  = 0
!    with 
!       Sij = 2 eta ( Eij - Ekk / 3 delta_ij) 
!       Kcp = b(D) B^(1/n) ED^{(n-1)/n}
!       ED^2 = 2 eij.eij/a(D) + Ekk^2/b(D)
!       eta = 2/a(D) B^{-1/n} ED^{(1-n)/n}
!       The relative density D is computed within this solver 
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
     INTEGER :: dim,  nComp, nd, nb 

     TYPE(ValueList_t),POINTER :: Material, BC, BodyForce
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement

     REAL(KIND=dp) :: RelativeChange, UNorm, PrevUNorm, &
         Normal(3), NewtonTol, NonlinearTol, s, Wn(2), MinSRInvariant, &
         IceDensity
         

     REAL(KIND=dp), ALLOCATABLE :: Basis(:),ddBasisddx(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: dBasisdx(:,:), SlipCoeff(:,:)
     REAL(KIND=dp) :: u,v,w,detJ
     
     LOGICAL :: stat, CSymmetry 
       
     INTEGER, PARAMETER :: INDi(1:6) = (/ 1, 2, 3, 1, 2, 3 /) ,&
           INDj(1:6)=(/ 1, 2, 3, 2, 3, 1 /)

     INTEGER :: NewtonIter, NonlinearIter

     TYPE(Variable_t), POINTER :: LumiSol, DensityVariable

     REAL(KIND=dp), POINTER ::  Lumi(:),  &
           ForceVector(:),  NodalLumi(:), LumiComp(:), &
           DensityValues(:)


     INTEGER, POINTER :: LumiPerm(:), NodeIndexes(:), &
                         DensityPerm(:)

     LOGICAL :: GotForceBC, GotIt, NewtonLinearization = .FALSE., &
                NormalTangential=.FALSE.

     INTEGER :: body_id,bf_id
     INTEGER :: old_body = -1
     LOGICAL :: AllocationsDone = .FALSE., FreeSurface, &
                Requal0
           
     REAL(KIND=dp), ALLOCATABLE:: LocalMassMatrix(:,:), &
       LocalStiffMatrix(:,:), LoadVector(:,:), LocalForce(:), &
       Alpha(:,:), Beta(:), SOL(:,:), & 
       BoundaryDispl(:), TimeForce(:), RefS(:), RefD(:), RefSpin(:), &
       LocalFluidity(:), Localfa(:), Localfb(:)
            
     INTEGER :: NumberOfBoundaryNodes
     INTEGER, POINTER :: BoundaryReorder(:)

     REAL(KIND=dp) :: Bu, Bv, Bw, RM(3,3)
     REAL(KIND=dp), POINTER :: BoundaryNormals(:,:), &
         BoundaryTangent1(:,:), BoundaryTangent2(:,:)

     REAL(KIND=dp) :: at, at0, CPUTime, RealTime

     LOGICAL, ALLOCATABLE :: ActiveList(:)
     REAL(KIND=dp), ALLOCATABLE ::  SaveValues(:), SaveRHS(:), Residual(:)

     TYPE(variable_t), POINTER :: var

!------------------------------------------------------------------------------
     SAVE NumberOfBoundaryNodes,BoundaryReorder,BoundaryNormals, &
              BoundaryTangent1, BoundaryTangent2

     SAVE TimeForce, Basis, dBasisdx, ddBasisddx
     SAVE LocalMassMatrix, LocalStiffMatrix, LoadVector, &
       LocalForce, ElementNodes, Alpha, Beta,  &
       AllocationsDone, BoundaryDispl, SOL, &
       NodalLumi, Wn, MinSRInvariant, old_body, &
       LocalFluidity, Localfa, Localfb, dim, IceDensity, nComp

     SAVE RefD, RefS, RefSpin, SlipCoeff, ActiveList, SaveValues, SaveRHS, Residual
              
!------------------------------------------------------------------------------
!  Read constants from constants section of SIF file
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
      IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

      LumiSol => Solver % Variable
      LumiPerm => LumiSol % Perm
      STDOFs =  LumiSol % DOFs
      Lumi => LumiSol % Values

      LocalNodes = COUNT( LumiPerm > 0 )
      IF ( LocalNodes <= 0 ) RETURN


      StiffMatrix => Solver % Matrix
      ForceVector => StiffMatrix % RHS
      UNorm = Solver % Variable % Norm


!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
      IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
        N = Solver % Mesh % MaxElementDOFs
        dim = CoordinateSystemDimension()
        nComp = STDOFs

       IF ( AllocationsDone ) THEN
         DEALLOCATE( ElementNodes % x,     &
                     ElementNodes % y,     &
                     ElementNodes % z,     &
                     BoundaryDispl,        &
                     LocalForce,           &
                     RefD, RefS, RefSpin,  &
                     LocalMassMatrix,      &
                     LocalStiffMatrix,     &
                     LoadVector, Alpha, Beta, &
                     SlipCoeff, LocalFluidity , &
                     Localfa, Localfb )
       END IF

       ALLOCATE( BoundaryDispl( N ), &
                 SOL( nComp, N ), &
                 LocalForce( STDOFs*N ),&
                 RefS(2*dim*LocalNodes ),&                              
                 RefD(2*dim*LocalNodes ),&                              
                 RefSpin((2*dim-3)*LocalNodes ),&                       
                 Basis( N ), dBasisdx( N,3 ), &
                 TimeForce( STDOFs*N ), &
                 LocalMassMatrix( STDOFs*N,STDOFs*N ),  &
                 LocalStiffMatrix( STDOFs*N,STDOFs*N ),  &
                 LoadVector( 4,N ), Alpha( 3,N ), Beta( N ), &
                 SlipCoeff(3,N), LocalFluidity(N), &
                 Localfa( N ), Localfb( N ),  STAT=istat )

       IF ( istat /= 0 ) THEN
          CALL Fatal( 'LumiSolve', 'Memory allocation error.' )
       END IF

       IF ( AllocationsDone ) THEN
         DEALLOCATE( ActiveList, SaveValues, SaveRHS, Residual )
       END IF
       ALLOCATE( ActiveList( Solver % Mesh % NumberOfNodes ),   &
                 SaveValues( SIZE( Solver % Matrix % Values) ), &
                 SaveRHS( Solver % Matrix % NumberOfRows ),     &
                 Residual( Solver % Matrix % NumberOfRows ), STAT=istat )

       IF ( istat /= 0 ) THEN
          CALL Fatal( 'LumiSolve', 'Memory allocation error.' )
       END IF

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

      FreeSurface = .FALSE.
      ActiveList = .FALSE.

      DO iter=1,NonlinearIter

       at  = CPUTime()
       at0 = RealTime()

       CALL Info( 'LumiSolve', ' ', Level=4 )
       CALL Info( 'LumiSolve', ' ', Level=4 )
       CALL Info( 'LumiSolve', &
                   '-------------------------------------',Level=4 )
       WRITE( Message, * ) 'POROUS MATERIAL FLOW SOLVER ITERATION', iter
       CALL Info( 'LumiSolve', Message,Level=4 )
       CALL Info( 'LumiSolve', &
                   '-------------------------------------',Level=4 )
       CALL Info( 'LumiSolve', ' ', Level=4 )
       CALL Info( 'LumiSolve', 'Starting assembly...',Level=4 )


!------------------------------------------------------------------------------
       CALL DefaultInitialize()
!------------------------------------------------------------------------------
       DO t=1,Solver % NumberOFActiveElements

         IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ',  &
             INT(100.0 - 100.0 * (Solver % NumberOfActiveElements-t) / &
             (1.0*Solver % NumberOfActiveElements)), ' % done'
           CALL Info( 'LumiSolve', Message, Level=5 )
           at0 = RealTime()
         END IF

         CurrentElement => GetActiveElement(t)
         n = GetElementNOFNodes()
         nd = GetElementNOFDOFs()
         nb = GetElementNOFBDOFs()

         CALL GetElementNodes( ElementNodes )
         Material => GetMaterial()

!------------------------------------------------------------------------------
!    Read in material constants from Material section
!------------------------------------------------------------------------------
         body_id = CurrentElement % BodyId
         IF (body_id /= old_body) Then 
            old_body = body_id
            Call  GetMaterialDefs()
        END IF

        LocalFluidity(1:n) = GetReal( Material, &
                         'Fluidity Parameter', GotIt )
        IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Fluidity Parameter not found. &
                            &Setting to 1.0'
         CALL INFO('LumiSolve', Message, Level = 20)
         LocalFluidity(1:n) = 1.0
        END IF

         Localfa(1:n) = 1.0
         Localfb(1:n) = 0.0
!
! Function a(D) (=1 if incompressible)
! 
!       Localfa(1:n) = GetReal( Material, &
!                        'FunctionA', n, NodeIndexes, GotIt )
!       IF (.NOT.GotIt) THEN
!        WRITE(Message,'(A)') 'FunctionA not found. &
!                           &Setting to 1.0'
!        CALL INFO('LumiSolve', Message, Level = 20)
!       END IF
!
! Function b(D) (=0 if incompressible)
! 
!       Localfb(1:n) = GetReal( Material, &
!                        'FunctionB', n, NodeIndexes, GotIt )
!       IF (.NOT.GotIt) THEN
!        WRITE(Message,'(A)') 'FunctionB not found. &
!                           &Setting to 0.0'
!        CALL INFO('LumiSolve', Message, Level = 20)
!       END IF

!------------------------------------------------------------------------------
!        Set body forces
!------------------------------------------------------------------------------
         LoadVector = 0.0D0

         BodyForce => GetBodyForce()
         IF ( ASSOCIATED( BodyForce ) ) THEN
           LoadVector(1,1:n) = LoadVector(1,1:n) + GetReal( &
                   BodyForce, 'Lumi Force 1', gotIt )
           LoadVector(2,1:n) = LoadVector(2,1:n) + GetReal( &
                   BodyForce, 'Lumi Force 2', gotIt )
           LoadVector(3,1:n) = LoadVector(3,1:n) + GetReal( & 
                   BodyForce, 'Lumi Force 3', gotIt )
         END IF
!------------------------------------------------------------------------------
!        Get element local stiffness & mass matrices
!------------------------------------------------------------------------------

         CALL GetVectorLocalSolution( SOL )

         CALL LocalMatrix( LocalMassMatrix, LocalStiffMatrix, &
              LocalForce, LoadVector, SOL, LocalFluidity, CurrentElement, n, &
              ElementNodes, Wn, MinSRInvariant, IceDensity, Localfa, Localfb, &
              nComp, dim )

          IF ( nb>0 ) CALL LCondensate( nd, nb, nComp, dim+1, &
                 LocalStiffMatrix, LocalForce )

!------------------------------------------------------------------------------
!        Update global matrices from local matrices 
!------------------------------------------------------------------------------
         IF ( TransientSimulation ) THEN 
           CALL Default1stOrderTime( LocalMassMatrix, &
                LocalStiffMatrix, LocalForce )
         END IF
         CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )
      END DO

      CALL Info( 'LumiSolve', 'Assembly done', Level=4 )

#if 1
!------------------------------------------------------------------------------
!     Neumann & Newton boundary conditions
!------------------------------------------------------------------------------
      DO t = 1, Model % NumberOFBoundaryElements

        CurrentElement => GetBoundaryElement(t)
        IF ( GetElementFamily() == 101 ) CYCLE
        IF ( .NOT. ActiveBoundaryElement() ) CYCLE

        n = GetElementNOFNodes()
        CALL GetElementNodes( ElementNodes )

        BC => GetBC()
        IF ( ASSOCIATED( BC ) ) THEN
            LoadVector = 0.0D0
            Alpha      = 0.0D0
            Beta       = 0.0D0
!------------------------------------------------------------------------------
!           Force in given direction BC: \tau\cdot n = F
!------------------------------------------------------------------------------
            GotForceBC = .FALSE.

            LoadVector(1,1:n) = GetReal( BC, 'Force 1', gotIt )
            GotForceBC = GotForceBC .OR. gotIt
            LoadVector(2,1:n) = GetReal( BC, 'Force 2', gotIt )
            GotForceBC = GotForceBC .OR. gotIt
            LoadVector(3,1:n) = GetReal( BC, 'Force 3', gotIt )
            GotForceBC = GotForceBC .OR. gotIt

            Beta(1:n) = GetReal( BC, 'Normal Force', gotIt )
            GotForceBC = GotForceBC .OR. gotIt

!------------------------------------------------------------------------------
!             slip boundary condition BC: \tau\cdot n = R_k u_k
!------------------------------------------------------------------------------

              SlipCoeff = 0.0d0
              SlipCoeff(1,1:n) =  GetReal( BC, 'Lumi Slip Coeff 1', gotit )
              GotForceBC = GotForceBC .OR. gotIt
              SlipCoeff(2,1:n) =  GetReal( BC, 'Lumi Slip Coeff 2', GotIt )
              GotForceBC = GotForceBC .OR. gotIt
              SlipCoeff(3,1:n) =  GetReal( BC, 'Lumi Slip Coeff 3', gotIt )
              GotForceBC = GotForceBC .OR. gotIt

              NormalTangential = GetLogical( BC, &
                     'Normal-Tangential Lumi', GotIt )
               
            IF ( .NOT.GotForceBC ) CYCLE
!------------------------------------------------------------------------------
            CALL LocalMatrixBoundary( LocalStiffMatrix, LocalForce, &
                 LoadVector, Alpha, Beta, SlipCoeff, NormalTangential, &
                 CurrentElement, n, ElementNodes, nComp )
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
#endif
!------------------------------------------------------------------------------

      CALL DefaultFinishAssembly()
      CALL DefaultDirichletBCs()

      n = Solver % Mesh % NumberOfNodes
      ActiveList = ActiveList .OR. Solver % Variable % Values( &
              nComp*Solver % Variable % Perm(1:n) )>1.0

      SaveValues = Solver % Matrix % Values
      SaveRHS    = Solver % Matrix % RHS
      DO i=1,Solver % Mesh % NumberOfNodes
        IF ( ActiveList(i) ) THEN
          j = nComp * Solver % Variable % Perm(i)
          CALL ZeroRow( Solver % Matrix, j )
          CALL SetMatrixElement( Solver % Matrix, j,j,1.0d0 )
          Solver % Matrix % RHS(j) = 1.0d0
        END IF
      END DO
      PRINT*, 'NOF Constrained values: ', COUNT(ActiveList)

      CALL Info( 'LumiSolve', 'Set boundaries done', Level=4 )

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      PrevUNorm = UNorm

      UNorm = DefaultSolve()

      Residual = 0
      DO i=1, Solver % Mesh % NumberOfNodes
         k = nComp*Solver % Variable % Perm(i)
         DO j=Solver % Matrix % Rows(k),Solver % Matrix % Rows(k+1)-1
           Residual(i) = Residual(i) + SaveValues(j) * &
             Solver % Variable % Values(Solver % Matrix % Cols(j))
         END DO
         Residual(i) = Residual(i) - SaveRHS(k)
         IF ( Residual(i) > 100*AEPS ) ActiveList(i) = .FALSE.
      END DO
      PRINT*,SUM(Residual**2)
      
      UNorm = 0.0_dp
      n = Solver % Variable % DOFs
      DO i=1,n-2
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
      CALL Info( 'LumiSolve', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'LumiSolve', Message, Level=4 )

!------------------------------------------------------------------------------
      IF ( RelativeChange < NewtonTol .OR. &
             iter > NewtonIter ) NewtonLinearization = .TRUE.

      IF ( RelativeChange < NonLinearTol ) EXIT

!------------------------------------------------------------------------------
    END DO ! of nonlinear iter
!------------------------------------------------------------------------------

var => variableGet( Solver % Mesh % Variables, 'Residual' )
var % values(var % Perm(1:Solver % mesh % numberofnodes)) = Residual

      
CONTAINS

      SUBROUTINE GetMaterialDefs()

      Wn(2) = GetConstReal( Material , 'Powerlaw Exponent', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable  Powerlaw Exponent not found. &
                                    & Setting to 1.0'
         CALL INFO('LumiSolve', Message, Level = 20)
         Wn(2) = 1.0
      ELSE
       WRITE(Message,'(A,F10.4)') 'Powerlaw Exponent = ',   Wn(2)
       CALL INFO('LumiSolve', Message, Level = 20)
       END IF

! Get the Minimum value of the Effective Strain rate 
      MinSRInvariant = 10e-20_dp 
      IF ( Wn(2) > 1.0 ) THEN
        MinSRInvariant =  &
             GetConstReal( Material, 'Min Second Invariant', GotIt )
        IF (.NOT.GotIt) THEN
          WRITE(Message,'(A)') 'Variable Min Second Invariant not &
                    &found. Setting to 10.0^-20 )'
          CALL INFO('LumiSolve', Message, Level = 20)
        ELSE
          WRITE(Message,'(A,E14.8)') 'Min Second Invariant = ', MinSRInvariant
          CALL INFO('LumiSolve', Message, Level = 20)
        END IF
      END IF

! Get the absolute ice density
      IceDensity = GetConstReal( Material , 'Ice Density', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable  Ice Density not found. &
                                    & Setting to 910.0'
         CALL INFO('LumiSolve', Message, Level = 20)
         IceDensity = 910.0_dp 
      ELSE
       WRITE(Message,'(A,F10.4)') ' Ice Density = ', IceDensity 
       CALL INFO('LumiSolve', Message, Level = 20)
       END IF

!------------------------------------------------------------------------------
      END SUBROUTINE GetMaterialDefs
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      SUBROUTINE LocalMatrix( MASS, STIFF, FORCE,       &
              LoadVector, SOL, NodalFluidity, Element, n, Nodes, Wn, MinSRInvariant, &
              IceDensity, Nodalfa, Nodalfb, nComp, dim )
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: STIFF(:,:), MASS(:,:), FORCE(:), SOL(:,:)
     REAL(KIND=dp) :: LoadVector(:,:)
     REAL(KIND=dp) :: Wn(2), MinSRInvariant, IceDensity
     REAL(KIND=dp), DIMENSION(:) :: NodalFluidity, Nodalfa, Nodalfb
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
     INTEGER :: n
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(nd)
     REAL(KIND=dp) :: dBasisdx(nd,3), DetJ

     REAL(KIND=dp) :: Load(3), density, fa, fb 
     REAL(kind=dp) :: Bg, BGlenT
     REAL(KIND=dp) :: A(nComp,nComp), M(nComp,nComp), F(nComp)

     REAL(KIND=dp) :: nn, ss, LGrad(3,3), SR(3,3), Velo(3)

     INTEGER :: i, j, k, p, q, t, dim, NBasis, ind(3)

     REAL(KIND=dp) :: s, u, v ,w, Radius
  
     REAL(KIND=dp) :: Em, eta, Kcp, ParameterA, ParameterB
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     INTEGER :: N_Integ, iP, iRho,nComp
     LOGICAL :: stat, CSymmetry
!------------------------------------------------------------------------------
      CSymmetry = CurrentCoordinateSystem() == AxisSymmetric
      dim = CoordinateSystemDimension()

      iP    = dim+1
      iRho  = dim+2

      FORCE = 0.0d0
      STIFF = 0.0d0
      MASS  = 0.0d0

!    
!    Integration stuff
!    
!     NBasis = 2*n
!     IntegStuff = GaussPoints( Element, Element % Type % GaussPoints2 )

      NBasis = nd+nb
      IntegStuff = GaussPoints( Element )

      !
      ! Now we start integrating:
      ! -------------------------
      DO t=1,IntegStuff % n

         u = IntegStuff % u(t)
         v = IntegStuff % v(t)
         w = IntegStuff % w(t)

         stat = ElementInfo( Element,Nodes,u,v,w,DetJ,Basis,dBasisdx )
         s = DetJ * IntegStuff % s(t)
         IF ( CSymmetry ) THEN
            Radius = SUM( Nodes % x(1:n) * Basis(1:n) )
            s = s * Radius
         END IF

         !  
         ! Force at integration point:
         ! ---------------------------
         Load = 0.0d0
         DO i=1,dim
           Load(i) = SUM( LoadVector(i,1:n)*Basis(1:n) )
         END DO

         !
         ! Variables at the integration point:
         ! -----------------------------------
         !     Density = SUM( NodalDensity(1:n)*Basis(1:n) )
         Density = SUM( SOL(iRho,1:n)*Basis(1:n) )
         Velo(1:dim) = MATMUL( SOL(1:dim,1:n), Basis(1:n) )

         !     fa = SUM( Nodalfa(1:n)*Basis(1:n) )
         !     fb = SUM( Nodalfb(1:n)*Basis(1:n) )
         fa = ParameterA(Density)
         fb = ParameterB(Density)

         Wn(1) = SUM( NodalFluidity(1:n)*Basis(1:n) )
         Bg=Wn(1)                      

         !
         ! Case non-linear calculate E_D^2 = gamma_e^2/fa + E_m^2/fb
         ! ----------------------------------------------------------
         ss = 1.0_dp
         IF ( Wn(2) > 1.0 ) THEN
           LGrad = 0
           LGrad(1:dim,1:dim) = MATMUL( SOL(1:dim,1:nd), dBasisdx(1:nd,1:dim) )
           SR = 0.5 * ( LGrad + TRANSPOSE(LGrad) )
      
              IF ( CSymmetry ) THEN
             SR(1,3) = 0.0
             SR(2,3) = 0.0
             SR(3,1) = 0.0
             SR(3,2) = 0.0
             SR(3,3) = 0.0
             IF ( Radius > 10*AEPS ) THEN
               SR(3,3) = SUM( SOL(1,1:nd) * Basis(1:nd) ) /Radius
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
         ! Bulk effective viscosity and Compressibility parameter Kcp :
         ! ------------------------------------------------------------     
         eta = ss / (fa * Bg**( 1.0 / Wn( 2 ) ) ) 
         Kcp = fb * Bg**( 1.0 / Wn( 2 ) ) / ss 

         !
         ! Loop over basis functions (of both unknowns and weights):
         ! ---------------------------------------------------------
         DO p=1,NBasis
          DO q=1,NBasis
            A = 0.0_dp
            M = 0.0_dp
            DO i=1,dim
              M(i,i) = M(i,i) + IceDensity*Density*Basis(q)*Basis(p) 
              DO j=1,dim

! terms 2 eta Eij dEij
                A(i,i) = A(i,i) + eta * dbasisdx(q,j)*dbasisdx(p,j)
                A(i,j) = A(i,j) + eta * dbasisdx(q,i)*dbasisdx(p,j)

! terms 2 / 3 eta Eii dEii
                A(i,j) = A(i,j) - 2.0d0/3.0d0 * eta*dbasisdx(q,j)*dbasisdx(p,i)
              END DO

! Pressure gradient, Continuity equation
              A(i,iP) = A(i,iP) - dBasisdx(p,i) * Basis(q)
              A(iP,i) = A(iP,i) + dBasisdx(q,i) * Basis(p)
            END DO
            A(iP,iP) = A(iP,iP) +  Kcp * Basis(q) * Basis(p)

            IF ( CSymmetry ) A(1,iP) =  A(1,iP) - Basis(p) * Basis(q) / Radius
            IF ( CSymmetry ) A(iP,1) =  A(iP,1) + Basis(p) * Basis(q) / Radius
!
! Mass conservation equation dD/dt + div(Du) = 0 (Integrated by parts)
! 0.5 + 0.5 = 1 (should try other decomposition ...
!
            M(iRho,iRho) = M(iRho,iRho) + Basis(q) * Basis(p) 
            DO i=1,dim
              A(iRho,i) = A(iRho,i) + Density * dBasisdx(q,i) * Basis(p)
              A(iRho,iRho) = A(iRho,iRho) + Velo(i) * dBasisdx(q,i) * Basis(p)
            END DO

! Add Nodal Matrix to Element Marix
            DO i=1, nComp
              DO j=1, nComp
                STIFF(nComp*(p-1)+i,nComp*(q-1)+j) = &
                     STIFF(nComp*(p-1)+i,nComp*(q-1)+j) + s * A(i,j)
                MASS(nComp*(p-1)+i,nComp*(q-1)+j) = &
                     MASS(nComp*(p-1)+i,nComp*(q-1)+j) + s * M(i,j)
              END DO
            END DO

          END DO !of q

! The righthand side...
           F = 0.0_dp
           DO i=1,dim
              F(i) = F(i) + Load(i) * Density * Basis(p)
           END DO
           DO i=1, nComp
             FORCE(nComp*(p-1)+i) = FORCE(nComp*(p-1)+i) + s * F(i)
           END DO
         END DO !of p
      END DO !of N_integ 

      DO j=n+1,nd
        i = nComp*(j-1) + iP
        STIFF(i,:) = 0.0_dp
        STIFF(:,i) = 0.0_dp
        STIFF(i,i) = 1.0_dp
        MASS(:,i)  = 0.0_dp
        MASS(i,:)  = 0.0_dp
        FORCE(i)   = 0.0_dp
!       i = nComp*(j-1) + iRho
!       STIFF(i,:) = 0.0_dp
!       STIFF(:,i) = 0.0_dp
!       STIFF(i,i) = 1.0_dp
!       MASS(:,i)  = 0.0_dp
!       MASS(i,:)  = 0.0_dp
!       FORCE(i)   = 0.0_dp
      END DO

!------------------------------------------------------------------------------
      END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      SUBROUTINE LocalMatrixBoundary( BoundaryMatrix, BoundaryVector, &
                 LoadVector, NodalAlpha, NodalBeta, NodalSlipCoeff, & 
                  NormalTangential, Element, n, Nodes, nComp )
                      
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:)
     REAL(KIND=dp) :: NodalAlpha(:,:),NodalBeta(:),LoadVector(:,:)
     REAL(KIND=dp) :: NodalSlipCoeff(:,:)
     TYPE(Element_t),POINTER  :: Element
     TYPE(Nodes_t)    :: Nodes
     LOGICAL :: NormalTangential
     INTEGER :: n, nComp
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(n),ddBasisddx(1,1,1)
     REAL(KIND=dp) :: dBasisdx(n,3),DetJ

     REAL(KIND=dp) :: u,v,w,s
     REAL(KIND=dp) :: Force(3), Alpha(3), Beta, Normal(3)
     REAL(KIND=dp), POINTER :: U_Integ(:), V_Integ(:), W_Integ(:), S_Integ(:)

     REAL(KIND=dp) :: Tangent(3),Tangent2(3),Vect(3), SlipCoeff
     INTEGER :: i,t,q,p,dim,N_Integ, c

     LOGICAL :: stat

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
!------------------------------------------------------------------------------

      dim = CoordinateSystemDimension()
      c=nComp

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
       stat = ElementInfo( Element, Nodes, u, v, w, DetJ, &
                Basis, dBasisdx, ddBasisddx, .FALSE. )

       s = DetJ * S_Integ(t)
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
           BoundaryMatrix((p-1)*c+i,(q-1)*c+i) =  &
             BoundaryMatrix((p-1)*c+i,(q-1)*c+i) + &
               s * Alpha(i) * Basis(q) * Basis(p)
         END DO
       END DO
      END DO

      DO q=1,N
       DO i=1,dim
         BoundaryVector((q-1)*c+i) = BoundaryVector((q-1)*c+i) + &
                   s * Basis(q) * Force(i)
       END DO
      END DO

      END DO
!------------------------------------------------------------------------------
      END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE LCondensate( N, nb, nComp, nRemove, K, F )
!------------------------------------------------------------------------------
      USE LinearAlgebra
      INTEGER :: N, nb, nComp, nRemove
      REAL(KIND=dp) :: K(:,:),F(:), Kbb(Nb*(nComp-1),Nb*(nComp-1)), &
       Kbl(nb*(nComp-1),n*nComp),Klb(n*nComp,nb*(nComp-1)),Fb(nb*(nComp-1))

      INTEGER :: m, i, j, l, p, Cdofs(nComp*n), Bdofs((nComp-1)*nb)

      m = 0
      DO p = 1,n
        DO i = 1,nComp
          m = m + 1
          Cdofs(m) = nComp*(p-1) + i
        END DO
      END DO
      
      m = 0
      DO p = 1,nb
        DO i = 1,nComp
          IF ( i==nRemove ) CYCLE
          m = m + 1
          Bdofs(m) = nComp*(p-1) + i + n*nComp
        END DO
      END DO

      Kbb = K(Bdofs,Bdofs)
      Kbl = K(Bdofs,Cdofs)
      Klb = K(Cdofs,Bdofs)
      Fb  = F(Bdofs)
      CALL InvertMatrix( Kbb,Nb*(nComp-1) )

      F(1:nComp*n) = F(1:nComp*n) - MATMUL( Klb, MATMUL( Kbb, Fb ) )
      K(1:nComp*n,1:nComp*n) = &
           K(1:nComp*n,1:nComp*n) - MATMUL( Klb, MATMUL( Kbb,Kbl ) )
!------------------------------------------------------------------------------
    END SUBROUTINE LCondensate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      END SUBROUTINE LumiSolver
!------------------------------------------------------------------------------
