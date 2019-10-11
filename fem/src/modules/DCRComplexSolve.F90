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
! *****************************************************************************/
!
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Mikko Lyly
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 04 Oct 2000
! *
! ****************************************************************************/

!------------------------------------------------------------------------------
!>  Solve the complex diffusion-convection-reaction equation.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE DCRComplexSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE Types
  USE Lists
  USE Adaptive
  USE Integration
  USE ElementDescription
  USE SolverUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(Nodes_t) :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement

  INTEGER, POINTER :: NodeIndexes(:)

  LOGICAL :: AllocationsDone = .FALSE., Bubbles, GotIt, notScalar = .TRUE., stat

  INTEGER, POINTER :: PressurePerm(:)
  REAL(KIND=dp), POINTER :: Pressure(:), ForceVector(:)

  INTEGER :: iter, i, j, k, n, t, istat, eq, LocalNodes
  REAL(KIND=dp) :: Norm, PrevNorm, RelativeChange

  TYPE(ValueList_t), POINTER :: Material

  INTEGER :: NonlinearIter
  REAL(KIND=dp) :: NonlinearTol,s

  REAL(KIND=dp), ALLOCATABLE :: LocalStiffMatrix(:,:), Load(:,:), Work(:), &
       LocalForce(:), &
       Amatrix(:,:,:), AvectorReal(:,:), AvectorImag(:,:), AscalarReal(:), &
       AscalarImag(:), Bvector(:,:), BscalarReal(:), BscalarImag(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: EquationName

  SAVE LocalStiffMatrix, Work, Load, LocalForce, ElementNodes, &
       AllocationsDone, &
       Amatrix, AvectorReal, AvectorImag, AscalarReal, AscalarImag, &
       Bvector, BscalarReal, BscalarImag
#ifdef USE_ISO_C_BINDINGS
   REAL(KIND=dp) :: at,at0,totat,st,totst,t1
#else
   REAL(KIND=dp) :: at,at0,totat,st,totst,t1,CPUTime,RealTime
#endif
!------------------------------------------------------------------------------
     INTERFACE
        FUNCTION DCRBoundaryResidual( Model,Edge,Mesh,Quant,Perm,Gnorm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Edge
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
          INTEGER :: Perm(:)
        END FUNCTION DCRBoundaryResidual

        FUNCTION DCREdgeResidual( Model,Edge,Mesh,Quant,Perm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Edge
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2)
          INTEGER :: Perm(:)
        END FUNCTION DCREdgeResidual

        FUNCTION DCRInsideResidual( Model,Element,Mesh,Quant,Perm, Fnorm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Element
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
          INTEGER :: Perm(:)
        END FUNCTION DCRInsideResidual
     END INTERFACE

!------------------------------------------------------------------------------
! Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT.ASSOCIATED( Solver % Matrix ) ) RETURN
  Solver % Matrix % Complex = .TRUE.

  Pressure     => Solver % Variable % Values
  PressurePerm => Solver % Variable % Perm

  LocalNodes = COUNT( PressurePerm > 0 )
  IF ( LocalNodes <= 0 ) RETURN

  StiffMatrix => Solver % Matrix
  ForceVector => StiffMatrix % RHS

  Norm = Solver % Variable % Norm

!------------------------------------------------------------------------------
! Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
     N = Solver % Mesh % MaxElementNodes

     IF ( AllocationsDone ) THEN
        DEALLOCATE(                 &
             ElementNodes % x,      &
             ElementNodes % y,      &
             ElementNodes % z,      &
             LocalForce,            &
             Work,                  &
             LocalStiffMatrix,      &
             Load, &
             Amatrix, AvectorReal, AvectorImag, AscalarReal, AscalarImag, &
             Bvector, BscalarReal, BscalarImag )
     END IF

     ALLOCATE( ElementNodes % x( N ),  &
          ElementNodes % y( N ),       &
          ElementNodes % z( N ),       &
          LocalForce( 2*N ),           &
          Work( N ),                   &
          LocalStiffMatrix( 2*N,2*N ), &
          Amatrix(3,3,N), AvectorReal(3,N), AvectorImag(3,N), &
          AscalarReal(N), AscalarImag(N), Bvector(3,N), &
          BscalarReal(N), BscalarImag(N), &
          Load( 2,N ), STAT=istat )

     IF ( istat /= 0 ) THEN
        CALL Fatal( 'DCRComplexSolve', 'Memory allocation error.' )
     END IF

     AllocationsDone = .TRUE.
  END IF
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Do some additional initialization, and go for it
!------------------------------------------------------------------------------
  NonlinearTol = ListGetConstReal( Solver % Values, &
       'Nonlinear System Convergence Tolerance',GotIt )

  NonlinearIter = ListGetInteger( Solver % Values, &
       'Nonlinear System Max Iterations', GotIt )

  IF ( .NOT.GotIt ) NonlinearIter = 1

  EquationName = ListGetString( Solver % Values, 'Equation' )

  Bubbles = ListGetLogical( Solver % Values, 'Bubbles', GotIt )

!------------------------------------------------------------------------------
! Iterate over any nonlinearity of material or source
!------------------------------------------------------------------------------
  totat = 0.0d0
  totst = 0.0d0

  DO iter=1,NonlinearIter
!------------------------------------------------------------------------------
     at  = CPUTime()
     at0 = RealTime()

     CALL Info( 'DCRComplexSolve', ' ', Level=4 )
     CALL Info( 'DCRComplexSolve', '-------------------------------------', Level=4 )
     WRITE( Message, * ) 'DCRComplex iteration', iter
     CALL Info( 'DCRComplexSolve', Message, Level=4 )
     CALL Info( 'DCRComplexSolve', '-------------------------------------', Level=4 )
     CALL Info( 'DCRComplexSolve', ' ', Level=4 )
     CALL Info( 'DCRComplexSolve', 'Starting Assmebly', Level=4 )

     CALL InitializeToZero( StiffMatrix, ForceVector )
!
!    Do the bulk assembly:
!    ---------------------

!------------------------------------------------------------------------------
     DO t=1,Solver % NumberOfActiveElements
!------------------------------------------------------------------------------
        IF ( RealTime() - at0 > 1.0 ) THEN
          WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
           (Solver % NumberOfActiveElements-t) / &
                       (1.0*Solver % NumberOfActiveElements)), ' % done'
          CALL Info( 'DCRComplexSolve', Message, Level=5 )
                      
          at0 = RealTime()
        END IF
!------------------------------------------------------------------------------
!       Check if this element belongs to a body where this equation
!       should be computed
!------------------------------------------------------------------------------
        CurrentElement => Solver % Mesh % Elements(Solver % ActiveElements(t))

!       IF ( .NOT. CheckElementEquation( Model, &
!            CurrentElement, EquationName ) ) CYCLE
!------------------------------------------------------------------------------
        Model % CurrentElement => CurrentElement

        n = CurrentElement % Type % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes

        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

!------------------------------------------------------------------------------
!       Get equation & material parameters
!------------------------------------------------------------------------------
        k = ListGetInteger( Model % Bodies( CurrentElement % &
         Bodyid ) % Values, 'Material',minv=1,maxv=Model % NumberOfMaterials )

        Material => Model % Materials(k) % Values

!------------------------------------------------------------------------------
!       Second and first order time derivative term coefficients on nodes
!------------------------------------------------------------------------------
        CALL InputTensor( Amatrix, notScalar, &
             'Amatrix', Material, n, NodeIndexes )

        CALL InputVector( AvectorReal, notScalar, &
             'Avector 1', Material, n, NodeIndexes )

        CALL InputVector( AvectorImag, notScalar, &
             'Avector 2', Material, n, NodeIndexes )

        AscalarReal(1:n) = ListGetReal( Material, &
             'Ascalar 1', n, NodeIndexes, GotIt)

        AscalarImag(1:n) = ListGetReal( Material, &
             'Ascalar 2', n, NodeIndexes, GotIt)

!------------------------------------------------------------------------------
!       The source term on nodes
!------------------------------------------------------------------------------
        k = ListGetInteger( Model % Bodies( CurrentElement % BodyId ) % &
              Values, 'Body Force', GotIt, 1, Model % NumberOFBodyForces )

        Load = 0.0d0
        IF ( k > 0 ) THEN
           Load(1,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
                'Pressure Source 1', n, NodeIndexes, GotIt )

           Load(2,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
                'Pressure Source 2', n, NodeIndexes, GotIt )
        END IF

!------------------------------------------------------------------------------
!       Get element local matrix and rhs vector
!------------------------------------------------------------------------------
        CALL LocalMatrix(  LocalStiffMatrix, LocalForce, & 
           Load, Bubbles, CurrentElement, n, ElementNodes, &
           Amatrix, AvectorReal, AvectorImag, AscalarReal, AscalarImag )

!------------------------------------------------------------------------------
!       Update global matrix and rhs vector from local matrix & vector
!------------------------------------------------------------------------------
        CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
             ForceVector, LocalForce, n, Solver % Variable % DOFs, &
                  PressurePerm(NodeIndexes) )
!------------------------------------------------------------------------------
     END DO
!------------------------------------------------------------------------------

!
!    Neumann & Newton BCs:
!    ---------------------

!------------------------------------------------------------------------------
     DO t = Solver % Mesh % NumberOfBulkElements + 1,  &
               Solver % Mesh % NumberOfBulkElements +  &
                  Solver % Mesh % NumberOfBoundaryElements
!------------------------------------------------------------------------------
        CurrentElement => Solver % Mesh % Elements(t)
        Model % CurrentElement => CurrentElement

!------------------------------------------------------------------------------
!       The element type 101 (point element) can only be used
!       to set Dirichlet BCs, so skip em at this stage.
!------------------------------------------------------------------------------
        IF ( CurrentElement % Type % ElementCode == 101 ) CYCLE

!------------------------------------------------------------------------------
        DO i=1,Model % NumberOfBCs
           IF ( CurrentElement % BoundaryInfo % Constraint == &
                Model % BCs(i) % Tag ) THEN
!------------------------------------------------------------------------------
              n = CurrentElement % Type % NumberOfNodes
              NodeIndexes => CurrentElement % NodeIndexes

              IF ( ANY( PressurePerm(NodeIndexes) == 0 ) ) CYCLE

              ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
              ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
              ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

              Load(1,1:n) = ListGetReal( Model % BCs(i) % Values, &
                   'Wave Flux 1', n, NodeIndexes, GotIt )

              Load(2,1:n) = ListGetReal( Model % BCs(i) % Values, &
                   'Wave Flux 2', n, NodeIndexes, GotIt )

              CALL InputVector( Bvector, notScalar, &
                   'Bvector', Model % BCs(i) % Values,  n, NodeIndexes )

              BscalarReal(1:n) = ListGetReal( Model % BCs(i) % Values, &
                   'Bscalar 1', n, NodeIndexes, GotIt)
              
              BscalarImag(1:n) = ListGetReal( Model % BCs(i) % Values, &
                   'Bscalar 2', n, NodeIndexes, GotIt)

!------------------------------------------------------------------------------
!             Get element local matrix and rhs vector
!------------------------------------------------------------------------------
              CALL LocalMatrixBoundary(  LocalStiffMatrix, LocalForce, &
                   Load, CurrentElement, n, ElementNodes, & 
                   Amatrix, AvectorReal, AvectorImag, AscalarReal, AscalarImag, &
                   Bvector, BscalarReal, BscalarImag )

!------------------------------------------------------------------------------
!             Update global matrix and rhs vector from local matrix & vector
!------------------------------------------------------------------------------
              CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
                  ForceVector, LocalForce, n, Solver % Variable % DOFs,  &
                      PressurePerm(NodeIndexes) )
!------------------------------------------------------------------------------
           END IF
        END DO
!------------------------------------------------------------------------------
     END DO
!------------------------------------------------------------------------------

     CALL FinishAssembly( Solver, ForceVector )
!
!    Dirichlet BCs:
!    --------------
     CALL SetDirichletBoundaries( Model, StiffMatrix, ForceVector, &
          ComponentName(Solver % Variable,1), 1, &
             Solver % Variable % DOFs, PressurePerm )

     CALL SetDirichletBoundaries( Model, StiffMatrix, ForceVector, &
          ComponentName(Solver % Variable,2), 2, &
             Solver % Variable % DOFs, PressurePerm )

     CALL Info( 'DCRComplexSolve', 'Assembly done', Level=4 )

!
!    Solve the system and we are done:
!    ---------------------------------
     PrevNorm = Norm
     at = CPUTime() - at
     st = CPUTime()

     CALL SolveSystem( StiffMatrix, ParMatrix, ForceVector, &
          Pressure, Norm, Solver % Variable % DOFs, Solver )

     st = CPUTIme()-st
     totat = totat + at
     totst = totst + st
     WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
     CALL Info( 'DCRComplexSolve', Message, Level=4 )
     WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
     CALL Info( 'DCRComplexSolve', Message, Level=4 )

!------------------------------------------------------------------------------
     IF ( PrevNorm + Norm /= 0.0d0 ) THEN
        RelativeChange = 2*ABS(PrevNorm - Norm) / (PrevNorm + Norm)
     ELSE
        RelativeChange = 0.0d0
     END IF

     CALL Info( 'DCRComplexSolve', ' ', Level=4 )
     WRITE( Message, * ) 'Result Norm    : ',Norm
     CALL Info( 'DCRComplexSolve', Message, Level=4 )
     WRITE( Message, * ) 'Relative Change: ',RelativeChange
     CALL Info( 'DCRComplexSolve', Message, Level=4 )

     IF ( RelativeChange < NonlinearTol ) EXIT
!------------------------------------------------------------------------------
  END DO ! of nonlinear iteration
!------------------------------------------------------------------------------
   IF ( ListGetLogical( Solver % Values, 'Adaptive Mesh Refinement', GotIt ) ) &
      CALL RefineMesh( Model,Solver,Pressure,PressurePerm, &
            DCRInsideResidual, DCREdgeResidual, DCRBoundaryResidual )

CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE InputTensor( Tensor, IsScalar, Name, Material, n, NodeIndexes )
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Tensor(:,:,:)
      INTEGER :: n, NodeIndexes(:)
      LOGICAL :: IsScalar
      CHARACTER(LEN=*) :: Name
      TYPE(ValueList_t), POINTER :: Material
!------------------------------------------------------------------------------
      LOGICAL :: FirstTime = .TRUE., stat
      REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

      SAVE FirstTime, Hwrk
!------------------------------------------------------------------------------
      IF ( FirstTime ) THEN
         NULLIFY( Hwrk )
         FirstTime = .FALSE.
      END IF

      Tensor = 0.0d0

      CALL ListGetRealArray( Material, Name, Hwrk, n, NodeIndexes, stat )
      IsScalar = SIZE(HWrk,1) == 1 .AND. SIZE(HWrk,2) == 1

      IF ( .NOT. stat ) RETURN

      IF ( SIZE(Hwrk,1) == 1 ) THEN

         DO i=1,MIN(3,SIZE(Hwrk,2))
            Tensor( i,i,1:n ) = Hwrk( 1,1,1:n )
         END DO

      ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN

         DO i=1,MIN(3,SIZE(Hwrk,1))
            Tensor(i,i,1:n) = Hwrk(i,1,1:n)
         END DO

      ELSE

        DO i=1,MIN(3,SIZE(Hwrk,1))
           DO j=1,MIN(3,SIZE(Hwrk,2))
              Tensor( i,j,1:n ) = Hwrk(i,j,1:n)
           END DO
        END DO

      END IF
!------------------------------------------------------------------------------
   END SUBROUTINE InputTensor
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE InputVector( Tensor, IsScalar, Name, Material, n, NodeIndexes )
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Tensor(:,:)
      INTEGER :: n, NodeIndexes(:)
      LOGICAL :: IsScalar
      CHARACTER(LEN=*) :: Name
      TYPE(ValueList_t), POINTER :: Material
!------------------------------------------------------------------------------
      LOGICAL :: FirstTime = .TRUE., stat
      REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

      SAVE FirstTime, Hwrk
!------------------------------------------------------------------------------
      IF ( FirstTime ) THEN
         NULLIFY( Hwrk )
         FirstTime = .FALSE.
      END IF

      Tensor = 0.0d0

      CALL ListGetRealArray( Material, Name, Hwrk, n, NodeIndexes, stat )
      IsScalar = SIZE(HWrk,1) == 1 .AND. SIZE(HWrk,2) == 1

      IF ( .NOT. stat ) RETURN

      IF ( SIZE(Hwrk,1) == 1 ) THEN

         DO i=1,MIN(3,SIZE(Hwrk,2))
            Tensor( i,1:n ) = Hwrk( 1,1,1:n )
         END DO

      ELSE

        DO i=1,MIN(3,SIZE(Hwrk,1))
           Tensor( i,1:n ) = Hwrk( i,1,1:n )
        END DO

      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE InputVector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------

  SUBROUTINE LocalMatrix(  StiffMatrix, Force, Load, Bubbles, Element, n, &
       Nodes, Amatrix, AvectorReal, AvectorImag, AscalarReal, AscalarImag )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: StiffMatrix(:,:), Force(:), Load(:,:), &
         Amatrix(:,:,:), AvectorReal(:,:), AvectorImag(:,:), &
         AscalarReal(:), AscalarImag(:)
    LOGICAL :: Bubbles
    INTEGER :: n
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(2*n),dBasisdx(2*n,3)
    REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,M,D,L1,L2
    REAL(KIND=dp) :: DiffCoef(3,3), Velo(3)
    COMPLEX(KIND=dp) :: LSTIFF(2*n,2*n), LFORCE(2*n), A
    LOGICAL :: Stat
    INTEGER :: i,p,q,t,dim, NBasis, CoordSys
    TYPE(GaussIntegrationPoints_t) :: IntegStuff

    REAL(KIND=dp) :: X,Y,Z,Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3)
    REAL(kind=dp) :: A2(3,3), A1r(3), A1i(3), A0r, A0i
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    CoordSys = CurrentCoordinateSystem()

    Metric = 0.0d0
    Metric(1,1) = 1.0d0
    Metric(2,2) = 1.0d0
    Metric(3,3) = 1.0d0

    LSTIFF = 0.0d0
    LFORCE = 0.0d0
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    IF ( Bubbles ) THEN
       IntegStuff = GaussPoints( Element, Element % Type % GaussPoints2 )
       NBasis = 2*n
    ELSE
       NBasis = n
       IntegStuff = GaussPoints( Element )
    END IF
!------------------------------------------------------------------------------
    DO t=1,IntegStuff % n
       U = IntegStuff % u(t)
       V = IntegStuff % v(t)
       W = IntegStuff % w(t)
       S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
            Basis, dBasisdx, Bubbles=Bubbles )

       s = s * SqrtElementMetric
       IF ( CoordSys /= Cartesian ) THEN
          X = SUM( Nodes % X(1:n) * Basis(1:n) )
          Y = SUM( Nodes % Y(1:n) * Basis(1:n) )
          Z = SUM( Nodes % Z(1:n) * Basis(1:n) )
          CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
          s = s * SqrtMetric
       END IF
!------------------------------------------------------------------------------
!      The source term and the coefficient of the time derivative and 
!      diffusion terms at the integration point
!------------------------------------------------------------------------------
!       D  =  WaveNumber * SUM( Damping(1:n) * Basis(1:n) )
!       M  = -WaveNumber**2

       L1 = SUM( Load(1,1:n) * Basis(1:n) )
       L2 = SUM( Load(2,1:n) * Basis(1:n) )

       A2 = 0.0d0
       A1r = 0.0d0
       A1i = 0.0d0
       A0r = 0.0d0
       A0i = 0.0d0

       A0r = SUM( AscalarReal(1:n) * Basis(1:n) )
       A0i = SUM( AscalarImag(1:n) * Basis(1:n) )
       do i = 1,dim
          A1r(i) = A1r(i) + SUM( AvectorReal(i,1:n) * Basis(1:n) )
          A1i(i) = A1i(i) + SUM( AvectorImag(i,1:n) * Basis(1:n) )
          do j = 1,dim
             A2(i,j) = A2(i,j) + SUM( Amatrix(i,j,1:n) * Basis(1:n) )
          end do
       end do


!      Stiffness matrix and load vector
!      --------------------------------
       DO p=1,NBasis
          DO q=1,NBasis
             A = CMPLX( A0r, A0i,KIND=dp ) * Basis(q) * Basis(p)
             DO i=1,dim
                A = A + CMPLX( A1r(i), A1i(i), KIND=dp ) * dBasisdx(q,i) * basis(p)
                DO j=1,dim
                   DO k = 1,dim
                      A = A + Metric(i,j) * A2(i,k) * dBasisdx(q,k) * dBasisdx(p,j)
                   END DO
                END DO
             END DO
             LSTIFF(p,q) = LSTIFF(p,q) + s*A
          END DO
          LFORCE(p) = LFORCE(p) + s * Basis(p) * CMPLX( L1,L2,KIND=dp )
       END DO
    END DO
!------------------------------------------------------------------------------

    IF ( Bubbles ) THEN
       CALL CondensateP( n, n, LSTIFF, LFORCE )
    END IF

    DO i=1,n
       Force( 2*(i-1)+1 ) = REAL( LFORCE(i) )
       Force( 2*(i-1)+2 ) = AIMAG( LFORCE(i) )

       DO j=1,n
         StiffMatrix( 2*(i-1)+1, 2*(j-1)+1 ) =  REAL( LSTIFF(i,j) )
         StiffMatrix( 2*(i-1)+1, 2*(j-1)+2 ) = -AIMAG( LSTIFF(i,j) )
         StiffMatrix( 2*(i-1)+2, 2*(j-1)+1 ) =  AIMAG( LSTIFF(i,j) )
         StiffMatrix( 2*(i-1)+2, 2*(j-1)+2 ) =  REAL( LSTIFF(i,j) )
       END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBoundary(  StiffMatrix, Force, & 
       Load, Element, n, Nodes, &
       Amatrix, AvectorReal, AvectorImag, AscalarReal, AscalarImag, &
       Bvector, BscalarReal, BscalarImag )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: StiffMatrix(:,:),Force(:),Load(:,:)
!    REAL(KIND=dp) :: ConvVelo(:,:)
    INTEGER :: n
    REAL(kind=dp) :: Amatrix(:,:,:), AvectorReal(:,:), AvectorImag(:,:), &
         AscalarReal(:), AscalarImag(:), Bvector(:,:), BscalarReal(:), &
         BscalarImag(:)
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,L1,L2
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),X,Y,Z
    REAL(KIND=dp) :: Normal(3), Velo(3), NormVelo, TangVelo(3)
    REAL(kind=dp) :: A2(3,3), A1r(3), A1i(3), A0r, A0i, B1(3), B0r, B0i, &
         C1(3), C0
    COMPLEX(KIND=dp) :: LSTIFF(n,n), LFORCE(n), A
    LOGICAL :: Stat
    INTEGER :: i,p,q,t,dim,CoordSys
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    CoordSys = CurrentCoordinateSystem()

    LSTIFF = 0.0d0
    LFORCE = 0.0d0
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( Element )
!------------------------------------------------------------------------------
    DO t=1,IntegStuff % n
       U = IntegStuff % u(t)
       V = IntegStuff % v(t)
       W = IntegStuff % w(t)
       S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
              Basis, dBasisdx )

       s = s * SqrtElementMetric

       Normal = Normalvector(Element, Nodes, U, V, .TRUE.)

       A2 = 0.0d0
       A1r = 0.0d0
       A1i = 0.0d0
       A0r = 0.0d0
       A0i = 0.0d0

       A0r = SUM( AscalarReal(1:n) * Basis(1:n) )
       A0i = SUM( AscalarImag(1:n) * Basis(1:n) )
       do i = 1,dim
          A1r(i) = A1r(i) + SUM( AvectorReal(i,1:n) * Basis(1:n) )
          A1i(i) = A1i(i) + SUM( AvectorImag(i,1:n) * Basis(1:n) )
          do j = 1,dim
             A2(i,j) = A2(i,j) + SUM( Amatrix(i,j,1:n) * Basis(1:n) )
          end do
       end do

       B1 = 0.0d0
       B0r = 0.0d0
       B0i = 0.0d0

       B0r = SUM( BscalarReal(1:n) * Basis(1:n) )
       B0i = SUM( BscalarImag(1:n) * Basis(1:n) )
       do i = 1,dim
          B1(i) = B1(i) + SUM( Bvector(i,1:n) * Basis(1:n) )
       end do
       B1(1:dim) = B1(1:dim) + Normal(1:dim)
       
       C1 = 0.0d0
       C0 = 0.0d0
       do i = 1,dim
          do j = 1,dim
             C0 = C0 + Normal(i)*A2(i,j)*Normal(j)
          end do
       end do
       C0 = C0 / SUM( B1(1:dim) * Normal(1:dim) )
       C1 = C0*B1
       do i = 1,dim
          do j = 1,dim
             C1(i) = C1(i) - A2(i,j)*Normal(j)
          end do
       end do


!------------------------------------------------------------------------------
       L1 = SUM( Load(1,1:n) * Basis )
       L2 = SUM( Load(2,1:n) * Basis )
!------------------------------------------------------------------------------
       DO p=1,n
          DO q=1,n
             A = CMPLX( B0r, B0i, KIND=dp ) * C0 * Basis(p)*Basis(q)
             A = A + SUM( C1(1:dim) * dBasisdx(q,1:dim) ) * Basis(p)
             LSTIFF(p,q) = LSTIFF(p,q) + s * A
          END DO
          LFORCE(p) = LFORCE(p) + s * Basis(p) * C0 * CMPLX( L1, L2, KIND=dp ) 
       END DO
!------------------------------------------------------------------------------
    END DO
!------------------------------------------------------------------------------
    DO i=1,n
       Force( 2*(i-1)+1 ) =  REAL( LFORCE(i) )
       Force( 2*(i-1)+2 ) = AIMAG( LFORCE(i) )

       DO j=1,n
         StiffMatrix( 2*(i-1)+1, 2*(j-1)+1 ) =  REAL( LSTIFF(i,j) )
         StiffMatrix( 2*(i-1)+1, 2*(j-1)+2 ) = -AIMAG( LSTIFF(i,j) )
         StiffMatrix( 2*(i-1)+2, 2*(j-1)+1 ) =  AIMAG( LSTIFF(i,j) )
         StiffMatrix( 2*(i-1)+2, 2*(j-1)+2 ) =  REAL( LSTIFF(i,j) )
       END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE DCRComplexSolver
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION DCRBoundaryResidual( Model, Edge, Mesh, Quant, Perm,Gnorm ) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
     TYPE( Mesh_t ), POINTER    :: Mesh
     TYPE( Element_t ), POINTER :: Edge
!------------------------------------------------------------------------------
     TYPE(Nodes_t) :: Nodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element, Bndry

     INTEGER :: i,j,k,n,l,t,DIM,Pn,En
     LOGICAL :: stat, GotIt, gotWIreal,gotWIimag 

     LOGICAL :: notScalar = .TRUE.

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
     REAL(KIND=dp) :: Grad(3,3), Normal(3), EdgeLength
     REAL(KIND=dp) :: Source, ResidualNorm, Area
     REAL(kind=dp) :: B1(3), B0r, B0i, Greal, Gimag
     REAL(KIND=dp) :: u, v, w, s, detJ, ResidualReal, ResidualImag, WaveFlux(2)

     REAL(KIND=dp), ALLOCATABLE :: EdgeBasis(:), Basis(:)
     REAL(KIND=dp), ALLOCATABLE :: dBasisdx(:,:)
     REAL(KIND=dp), ALLOCATABLE :: Flux(:)
     REAL(KIND=dp), ALLOCATABLE :: Work(:)
     REAL(KIND=dp), ALLOCATABLE :: x(:), y(:), z(:)
     REAL(KIND=dp), ALLOCATABLE :: Pressure(:,:), FluxReal(:), FluxImag(:)
     REAL(KIND=dp), ALLOCATABLE :: Bvector(:,:), BScalarReal(:), BscalarImag(:)

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: Dirichlet = .FALSE.

!    Initialize:
!    -----------
     Indicator = 0.0d0
     Gnorm     = 0.0d0

     Metric = 0.0d0
     DO i=1,3
        Metric(i,i) = 1.0d0
     END DO

     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           DIM = 3
        CASE DEFAULT
           DIM = CoordinateSystemDimension()
     END SELECT
!    
!    ---------------------------------------------
     Element => Edge % BoundaryInfo % Left
     IF ( .NOT. ASSOCIATED( Element ) ) THEN
        Element => Edge % BoundaryInfo % Right
     ELSE IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) THEN
        Element => Edge % BoundaryInfo % Right
     END IF

     IF ( .NOT. ASSOCIATED( Element ) ) RETURN
     IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN

     En = Edge % TYPE % NumberOfNodes
     Pn = Element % TYPE % NumberOfNodes

     ALLOCATE( EdgeNodes % x(En), EdgeNodes % y(En), EdgeNodes % z(En) )

     EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
     EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
     EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)

     ALLOCATE( Nodes % x(Pn), Nodes % y(Pn), Nodes % z(Pn) )

     Nodes % x = Mesh % Nodes % x(Element % NodeIndexes)
     Nodes % y = Mesh % Nodes % y(Element % NodeIndexes)
     Nodes % z = Mesh % Nodes % z(Element % NodeIndexes)

     ALLOCATE( x(En), y(En), z(En), EdgeBasis(En), Basis(Pn), dBasisdx(Pn,3),   &
           Flux(en), Pressure(2,Pn), FluxReal(En), FluxImag(En), BVector(3,Pn), &
           BScalarReal(Pn), BScalarImag(Pn) )

     DO l = 1,En
       DO k = 1,Pn
          IF ( Edge % NodeIndexes(l) == Element % NodeIndexes(k) ) THEN
             x(l) = Element % TYPE % NodeU(k)
             y(l) = Element % TYPE % NodeV(k)
             z(l) = Element % TYPE % NodeW(k)
             EXIT
          END IF
       END DO
     END DO

!
!    Integrate square of residual over boundary element:
!    ---------------------------------------------------
     Indicator    = 0.0d0
     EdgeLength   = 0.0d0
     ResidualNorm = 0.0d0
     Gnorm = 0.0d0

     DO j=1,Model % NumberOfBCs
        IF ( Edge % BoundaryInfo % Constraint /= Model % BCs(j) % Tag ) CYCLE
!
!       ...given parameters:
!       --------------------

        CALL InputVector( Bvector, notScalar, 'Bvector', &
             Model % BCs(j) % Values, Pn, Element % NodeIndexes )
        
        BscalarReal(1:Pn) = ListGetReal( Model % BCs(j) % Values, &
             'Bscalar 1', Pn, Element % NodeIndexes, GotIt)

        BscalarImag(1:Pn) = ListGetReal( Model % BCs(j) % Values, &
             'Bscalar 2', Pn, Element % NodeIndexes, GotIt)

        FluxReal(1:En) = ListGetReal( Model % BCs(j) % Values, &
          'Wave Flux 1', En, Edge % NodeIndexes, gotIt )
        if( .not.GotIt ) FluxReal(1:En) = 0.0d0

        FluxImag(1:En) = ListGetReal( Model % BCs(j) % Values, &
          'Wave Flux 2', En, Edge % NodeIndexes, gotIt )
        if( .not.GotIt ) FluxImag(1:En) = 0.0d0
!
!       get material parameters:
!       ------------------------
        k = ListGetInteger(Model % Bodies(Element % BodyId) % Values,'Material', &
                   minv=1, maxv=Model % NumberOfMaterials)
!
!       elementwise nodal solution:
!       ---------------------------
        Pressure(1,1:Pn) = Quant( 2*Perm(Element % NodeIndexes)-1 )
        Pressure(2,1:Pn) = Quant( 2*Perm(Element % NodeIndexes)-0 )

!       do the integration:
!       -------------------
        EdgeLength   = 0.0d0
        ResidualNorm = 0.0d0
        Gnorm = 0.0d0

        IntegStuff = GaussPoints( Edge )

        DO t=1,IntegStuff % n
           u = IntegStuff % u(t)
           v = IntegStuff % v(t)
           w = IntegStuff % w(t)

           stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
               EdgeBasis, dBasisdx )

           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              s = IntegStuff % s(t) * detJ

           ELSE
              u = SUM( EdgeBasis(1:En) * EdgeNodes % x(1:En) )
              v = SUM( EdgeBasis(1:En) * EdgeNodes % y(1:En) )
              w = SUM( EdgeBasis(1:En) * EdgeNodes % z(1:En) )
      
              CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                         Symb, dSymb, u, v, w )

              s = IntegStuff % s(t) * detJ * SqrtMetric

           END IF

!
!          Integration point in parent element local
!          coordinates:
!          -----------------------------------------
           u = SUM( EdgeBasis(1:En) * x(1:En) )
           v = SUM( EdgeBasis(1:En) * y(1:En) )
           w = SUM( EdgeBasis(1:En) * z(1:En) )

           stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
                 Basis, dBasisdx )
!
!          Flux and parameters at integration point:
!          ----------------------------------------

           Normal = NormalVector( Edge, EdgeNodes, u, v, .TRUE. )

           B1 = 0.0d0
           B0r = 0.0d0
           B0i = 0.0d0
           
           B0r = SUM( BscalarReal(1:Pn) * Basis(1:Pn) )
           B0i = SUM( BscalarImag(1:Pn) * Basis(1:Pn) )
           do i = 1,dim
              B1(i) = B1(i) + SUM( Bvector(i,1:Pn) * Basis(1:Pn) )
           end do
           B1(1:dim) = B1(1:dim) + Normal(1:dim)

           WaveFlux = 0.0d0
           WaveFlux(1) = SUM( FluxReal(1:En) * EdgeBasis(1:En) )
           WaveFlux(2) = SUM( FluxImag(1:En) * EdgeBasis(1:En) )
           
           ResidualReal = -WaveFlux(1)
           ResidualImag = -WaveFlux(2)

           Greal = 0.0d0
           Gimag = 0.0d0
!
!          flux given by the computed solution, and 
!          force norm for scaling the residual:
!          -----------------------------------------
           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              DO k=1,DIM

                 ResidualReal = ResidualReal + &
                      SUM( dBasisdx(1:Pn,k) * Pressure(1,1:Pn) ) * B1(k)

                 ResidualImag = ResidualImag + &
                      SUM( dBasisdx(1:Pn,k) * Pressure(2,1:Pn) ) * B1(k) 

                 GReal = GReal + &
                      SUM( dBasisdx(1:Pn,k) * Pressure(1,1:Pn) ) * B1(k)

                 Gimag = Gimag + &
                      SUM( dBasisdx(1:Pn,k) * Pressure(2,1:Pn) ) * B1(k) 

              END DO

              ResidualReal = ResidualReal + &
                    SUM( Basis(1:Pn) * Pressure(1,1:Pn) ) * B0r &
                   -SUM( Basis(1:Pn) * Pressure(2,1:Pn) ) * B0i
              
              ResidualImag = ResidualImag + &
                    SUM( Basis(1:Pn) * Pressure(1,1:Pn) ) * B0i &
                   +SUM( Basis(1:Pn) * Pressure(2,1:Pn) ) * B0r

              Greal = Greal + &
                    SUM( Basis(1:Pn) * Pressure(1,1:Pn) ) * B0r &
                   -SUM( Basis(1:Pn) * Pressure(2,1:Pn) ) * B0i
              
              Gimag = Gimag + &
                    SUM( Basis(1:Pn) * Pressure(1,1:Pn) ) * B0i &
                   +SUM( Basis(1:Pn) * Pressure(2,1:Pn) ) * B0r

           ELSE
!              DO k=1,DIM
!                 DO l=1,DIM
!                    Residual = Residual + Metric(k,l) * Conductivity  * &
!                       SUM( dBasisdx(1:Pn,k) * Temperature(1:Pn) ) * Normal(l)
!
!                    Gnorm = Gnorm + s * (Metric(k,l) * Conductivity * &
!                      SUM(dBasisdx(1:Pn,k) * Temperature(1:Pn) ) * Normal(l))**2
!                 END DO
!              END DO
           END IF

           EdgeLength   = EdgeLength + s
           Gnorm = Gnorm + s * ( Greal **2 + Gimag**2 )

           IF ( .NOT. Dirichlet ) THEN
              ResidualNorm = ResidualNorm + s * (ResidualReal**2 + ResidualImag**2)
           END IF

        END DO

        EXIT

     END DO

     IF ( CoordinateSystemDimension() == 3 ) THEN
        EdgeLength = SQRT(EdgeLength)
     END IF

     Gnorm = EdgeLength * Gnorm

     Indicator = EdgeLength * ResidualNorm

     DEALLOCATE( Nodes % x, Nodes % y, Nodes % z)
     DEALLOCATE( EdgeNodes % x, EdgeNodes % y, EdgeNodes % z)

     DEALLOCATE( x, y, z, EdgeBasis, Basis, dBasisdx,   &
           Flux, Pressure, FluxReal, FluxImag, BVector, &
           BScalarReal, BScalarImag )



!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
   SUBROUTINE InputVector( Tensor, IsScalar, Name, Material, n, NodeIndexes )
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Tensor(:,:)
      INTEGER :: n, NodeIndexes(:)
      LOGICAL :: IsScalar
      CHARACTER(LEN=*) :: Name
      TYPE(ValueList_t), POINTER :: Material
!------------------------------------------------------------------------------
      LOGICAL :: FirstTime = .TRUE., stat
      REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

      SAVE FirstTime, Hwrk
!------------------------------------------------------------------------------
      IF ( FirstTime ) THEN
         NULLIFY( Hwrk )
         FirstTime = .FALSE.
      END IF

      Tensor = 0.0d0

      CALL ListGetRealArray( Material, Name, Hwrk, n, NodeIndexes, stat )
      IsScalar = SIZE(HWrk,1) == 1 .AND. SIZE(HWrk,2) == 1

      IF ( .NOT. stat ) RETURN

      IF ( SIZE(Hwrk,1) == 1 ) THEN

         DO i=1,MIN(3,SIZE(Hwrk,2))
            Tensor( i,1:n ) = Hwrk( 1,1,1:n )
         END DO

      ELSE

        DO i=1,MIN(3,SIZE(Hwrk,1))
           Tensor( i,1:n ) = Hwrk( i,1,1:n )
        END DO

      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE InputVector
!------------------------------------------------------------------------------
   END FUNCTION DCRBoundaryResidual
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION DCREdgeResidual( Model, Edge, Mesh, Quant, Perm ) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE

     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2)
     TYPE( Mesh_t ), POINTER    :: Mesh
     TYPE( Element_t ), POINTER :: Edge
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element, Bndry

     INTEGER :: i,j,k,l,n,t,DIM,En,Pn
     LOGICAL :: stat, GotIt
!     REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

     REAL(KIND=dp) :: Grad(3,3), Normal(3), EdgeLength, Jump, JumpReal, JumpImag, &
                      GradReal(3,3),GradImag(3,3)
     REAL(KIND=dp) :: u, v, w, s, detJ
     REAL(KIND=dp) :: Residual, ResidualNorm, Area
     REAL(kind=dp) :: A2(3,3), A1r(3), A1i(3), A0r, A0i
     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

     REAL(KIND=dp), ALLOCATABLE :: AvectorImag(:,:), AscalarReal(:), AscalarImag(:)
     REAL(KIND=dp), ALLOCATABLE :: Flux(:), x(:), y(:), z(:)
     REAL(KIND=dp), ALLOCATABLE :: Amatrix(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: EdgeBasis(:)
     REAL(KIND=dp), ALLOCATABLE :: AvectorReal(:,:)
     REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
     REAL(KIND=dp), ALLOCATABLE :: Temperature(:), Pressure(:,:)

     LOGICAL :: notScalar = .TRUE.

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

!     LOGICAL :: First = .TRUE.
!     SAVE Hwrk, First
!------------------------------------------------------------------------------

!    Initialize:
!    -----------

!     IF ( First ) THEN
!        First = .FALSE.
!        NULLIFY( Hwrk )
!     END IF

     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           DIM = 3
        CASE DEFAULT
           DIM = CoordinateSystemDimension()
     END SELECT

     Metric = 0.0d0
     DO i = 1,3
        Metric(i,i) = 1.0d0
     END DO

     Grad = 0.0d0
     GradReal = 0.0d0
     GradImag = 0.0d0
!
!    ---------------------------------------------

     Element => Edge % BoundaryInfo % Left
     n = Element % TYPE % NumberOfNodes

     Element => Edge % BoundaryInfo % Right
     n = MAX( n, Element % TYPE % NumberOfNodes )

     ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )

     En = Edge % TYPE % NumberOfNodes
     ALLOCATE( EdgeNodes % x(En), EdgeNodes % y(En), EdgeNodes % z(En) )

     EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
     EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
     EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)

     ALLOCATE( AvectorImag(3,n), AscalarReal(n), AscalarImag(n), Flux(en),   &
       x(En), y(En), z(En), AMatrix(3,3,n), EdgeBasis(En), AvectorReal(3,n), &
       Basis(n), dBasisdx(n,3), Temperature(n), Pressure(2,n) )

!    Integrate square of jump over edge:
!    -----------------------------------
     ResidualNorm = 0.0d0
     EdgeLength   = 0.0d0
     Indicator    = 0.0d0

     IntegStuff = GaussPoints( Edge )

     DO t=1,IntegStuff % n

        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
             EdgeBasis, dBasisdx )

        Normal = NormalVector( Edge, EdgeNodes, u, v, .FALSE. )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM( EdgeBasis(1:En) * EdgeNodes % x(1:En) )
           v = SUM( EdgeBasis(1:En) * EdgeNodes % y(1:En) )
           w = SUM( EdgeBasis(1:En) * EdgeNodes % z(1:En) )

           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                      Symb, dSymb, u, v, w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF
! 
!       Compute flux over the edge as seen by elements
!       on both sides of the edge:
!       ----------------------------------------------
        DO i = 1,2
           SELECT CASE(i)
              CASE(1)
                 Element => Edge % BoundaryInfo % Left
              CASE(2)
                 Element => Edge % BoundaryInfo % Right
           END SELECT
!
!          Can this really happen (maybe it can...)  ?      
!          -------------------------------------------
           IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) CYCLE
!
!          Next, get the integration point in parent
!          local coordinates:
!          -----------------------------------------
           Pn = Element % TYPE % NumberOfNodes

           DO j = 1,En
              DO k = 1,Pn
                 IF ( Edge % NodeIndexes(j) == Element % NodeIndexes(k) ) THEN
                    x(j) = Element % TYPE % NodeU(k)
                    y(j) = Element % TYPE % NodeV(k)
                    z(j) = Element % TYPE % NodeW(k)
                    EXIT
                 END IF
              END DO
           END DO

           u = SUM( EdgeBasis(1:En) * x(1:En) )
           v = SUM( EdgeBasis(1:En) * y(1:En) )
           w = SUM( EdgeBasis(1:En) * z(1:En) )
!
!          Get parent element basis & derivatives at the
!          integration point:
!          ---------------------------------------------
           Nodes % x(1:Pn) = Mesh % Nodes % x(Element % NodeIndexes)
           Nodes % y(1:Pn) = Mesh % Nodes % y(Element % NodeIndexes)
           Nodes % z(1:Pn) = Mesh % Nodes % z(Element % NodeIndexes)

           stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
             Basis, dBasisdx )
!
!          Material parameters:
!          --------------------
           k = ListGetInteger( Model % Bodies( &
                Element % BodyId) % Values, 'Material', minv=1,maxv=Model % NumberOfMaterials )

           CALL InputTensor( Amatrix, notScalar, &
                'Amatrix', Model % Materials(k) % Values, Pn, Element % NodeIndexes )
           
           CALL InputVector( AvectorReal, notScalar, &
                'Avector 1', Model % Materials(k) % Values, Pn, Element % NodeIndexes )
           
           CALL InputVector( AvectorImag, notScalar, &
                'Avector 2', Model % Materials(k) % Values, Pn, Element % NodeIndexes )
           
           AscalarReal(1:Pn) = ListGetReal( Model % Materials(k) % Values, &
                'Ascalar 1', Pn, Element % NodeIndexes, GotIt)
           
           AscalarImag(1:Pn) = ListGetReal( Model % Materials(k) % Values, &
                'Ascalar 2', Pn, Element % NodeIndexes, GotIt)

           A2 = 0.0d0
           A1r = 0.0d0
           A1i = 0.0d0
           A0r = 0.0d0
           A0i = 0.0d0
           
           A0r = SUM( AscalarReal(1:Pn) * Basis(1:Pn) )
           A0i = SUM( AscalarImag(1:Pn) * Basis(1:Pn) )
           do j = 1,dim
              A1r(j) = A1r(j) + SUM( AvectorReal(j,1:Pn) * Basis(1:Pn) )
              A1i(j) = A1i(j) + SUM( AvectorImag(j,1:Pn) * Basis(1:Pn) )
              do k = 1,dim
                 A2(j,k) = A2(j,k) + SUM( Amatrix(j,k,1:Pn) * Basis(1:Pn) )
              end do
           end do
!
!          Pressure at element nodal points
!          ---------------------------------
           do k = 1,2
              Pressure(k,1:Pn) = Quant( 2*Perm(Element % NodeIndexes)-2+k )
           end do
!
!          Finally, the flux:
!          ------------------
           DO j=1,DIM
              do k = 1,dim
                 GradReal(j,i) = GradReal(j,i) &
                      + SUM( dBasisdx(1:Pn,k) * Pressure(1,1:Pn) )*A2(j,k)
                 GradImag(j,i) = GradImag(j,i) &
                      + SUM( dBasisdx(1:Pn,k) * Pressure(2,1:Pn) )*A2(j,k)
              end do
           END DO

        END DO

!       Compute squre of the flux jump:
!       -------------------------------   
        EdgeLength  = EdgeLength + s
        JumpReal = 0.0d0
        JumpImag = 0.0d0

        DO k=1,DIM
           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              JumpReal = JumpReal + (GradReal(k,1) - GradReal(k,2)) * Normal(k)
              JumpImag = JumpImag + (GradImag(k,1) - GradImag(k,2)) * Normal(k)
           ELSE
!              DO l=1,DIM
!                 Jump = Jump + &
!                       Metric(k,l) * (Grad(k,1) - Grad(k,2)) * Normal(l)
!              END DO
           END IF
        END DO
        ResidualNorm = ResidualNorm + s * ( JumpReal**2 + JumpImag**2 )
     END DO

     IF ( CoordinateSystemDimension() == 3 ) THEN
        EdgeLength = SQRT(EdgeLength)
     END IF
     Indicator = EdgeLength * ResidualNorm

     DEALLOCATE( Nodes % x, Nodes % y, Nodes % z)
     DEALLOCATE( EdgeNodes % x, EdgeNodes % y, EdgeNodes % z)

     DEALLOCATE( AvectorImag, AscalarReal, AscalarImag, Flux,    &
       x, y, z, AMatrix, EdgeBasis, AvectorReal, Basis, dBasisdx,&
       Temperature, Pressure )
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
   SUBROUTINE InputTensor( Tensor, IsScalar, Name, Material, n, NodeIndexes )
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Tensor(:,:,:)
      INTEGER :: n, NodeIndexes(:)
      LOGICAL :: IsScalar
      CHARACTER(LEN=*) :: Name
      TYPE(ValueList_t), POINTER :: Material
!------------------------------------------------------------------------------
      LOGICAL :: FirstTime = .TRUE., stat
      REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

      SAVE FirstTime, Hwrk
!------------------------------------------------------------------------------
      IF ( FirstTime ) THEN
         NULLIFY( Hwrk )
         FirstTime = .FALSE.
      END IF

      Tensor = 0.0d0

      CALL ListGetRealArray( Material, Name, Hwrk, n, NodeIndexes, stat )
      IsScalar = SIZE(HWrk,1) == 1 .AND. SIZE(HWrk,2) == 1

      IF ( .NOT. stat ) RETURN

      IF ( SIZE(Hwrk,1) == 1 ) THEN

         DO i=1,MIN(3,SIZE(Hwrk,2))
            Tensor( i,i,1:n ) = Hwrk( 1,1,1:n )
         END DO

      ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN

         DO i=1,MIN(3,SIZE(Hwrk,1))
            Tensor(i,i,1:n) = Hwrk(i,1,1:n)
         END DO

      ELSE

        DO i=1,MIN(3,SIZE(Hwrk,1))
           DO j=1,MIN(3,SIZE(Hwrk,2))
              Tensor( i,j,1:n ) = Hwrk(i,j,1:n)
           END DO
        END DO

      END IF
!------------------------------------------------------------------------------
   END SUBROUTINE InputTensor
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE InputVector( Tensor, IsScalar, Name, Material, n, NodeIndexes )
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Tensor(:,:)
      INTEGER :: n, NodeIndexes(:)
      LOGICAL :: IsScalar
      CHARACTER(LEN=*) :: Name
      TYPE(ValueList_t), POINTER :: Material
!------------------------------------------------------------------------------
      LOGICAL :: FirstTime = .TRUE., stat
      REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

      SAVE FirstTime, Hwrk
!------------------------------------------------------------------------------
      IF ( FirstTime ) THEN
         NULLIFY( Hwrk )
         FirstTime = .FALSE.
      END IF

      Tensor = 0.0d0

      CALL ListGetRealArray( Material, Name, Hwrk, n, NodeIndexes, stat )
      IsScalar = SIZE(HWrk,1) == 1 .AND. SIZE(HWrk,2) == 1

      IF ( .NOT. stat ) RETURN

      IF ( SIZE(Hwrk,1) == 1 ) THEN

         DO i=1,MIN(3,SIZE(Hwrk,2))
            Tensor( i,1:n ) = Hwrk( 1,1,1:n )
         END DO

      ELSE

        DO i=1,MIN(3,SIZE(Hwrk,1))
           Tensor( i,1:n ) = Hwrk( i,1,1:n )
        END DO

      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE InputVector
!------------------------------------------------------------------------------


   END FUNCTION DCREdgeResidual
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   FUNCTION DCRInsideResidual( Model, Element, Mesh, &
        Quant, Perm, Fnorm ) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
     TYPE( Mesh_t ), POINTER    :: Mesh
     TYPE( Element_t ), POINTER :: Element
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes

     INTEGER :: i,j,k,l,n,t,DIM
     LOGICAL :: stat, GotIt
     TYPE( Variable_t ), POINTER :: Var

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
     REAL(KIND=dp) :: Source, Residual, ResidualNorm, Area
     REAL(kind=dp) :: A2(3,3), A1r(3), A1i(3), A0r, A0i
     REAL(KIND=dp) :: u, v, w, s, detJ, ResidualReal, ResidualImag

     REAL(KIND=dp), ALLOCATABLE :: NodalSource(:,:), Basis(:)
     REAL(KIND=dp), ALLOCATABLE :: dBasisdx(:,:), ddBasisddx(:,:,:), Pressione(:,:)
     REAL(KIND=dp), ALLOCATABLE :: Amatrix(:,:,:), AvectorReal(:,:)
     REAL(KIND=dp), ALLOCATABLE :: AvectorImag(:,:), AscalarReal(:), AscalarImag(:)

     LOGICAL :: notScalar = .TRUE.
     TYPE( ValueList_t ), POINTER :: Material
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

!     LOGICAL :: First = .TRUE.
!     SAVE Hwrk, First
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     Indicator = 0.0d0
     Fnorm     = 0.0d0
!
!    Check if this eq. computed in this element:
!    -------------------------------------------
     IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN

     Metric = 0.0d0
     DO i=1,3
        Metric(i,i) = 1.0d0
     END DO

     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           DIM = 3
        CASE DEFAULT
           DIM = CoordinateSystemDimension()
     END SELECT
!
!    Element nodal points:
!    ---------------------
     n = Element % TYPE % NumberOfNodes

     ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )
     Nodes % x = Mesh % Nodes % x(Element % NodeIndexes)
     Nodes % y = Mesh % Nodes % y(Element % NodeIndexes)
     Nodes % z = Mesh % Nodes % z(Element % NodeIndexes)

     ALLOCATE( NodalSource(2,n), Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), &
         Pressione(2,n), AMatrix(3,3,n), AvectorReal(3,n), AvectorImag(3,n), &
         AscalarReal(n), AscalarImag(n) )
!
!    Elementwise nodal solution:
!    ---------------------------
     Pressione(1,1:n) = Quant( 2*Perm(Element % NodeIndexes)-1 )
     Pressione(2,1:n) = Quant( 2*Perm(Element % NodeIndexes)-0 )
!
!    Material parameters: sound speed:
!    ---------------------------------
     k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Material', &
                  minv=1, maxv=Model % NumberOfMaterials )

     Material => Model % Materials(k) % Values

     CALL InputTensor( Amatrix, notScalar, &
          'Amatrix', Material, n, Element % NodeIndexes )
     
     CALL InputVector( AvectorReal, notScalar, &
          'Avector 1', Material, n, Element % NodeIndexes )
     
     CALL InputVector( AvectorImag, notScalar, &
          'Avector 2', Material, n, Element % NodeIndexes )
     
     AscalarReal(1:n) = ListGetReal( Material, &
          'Ascalar 1', n, Element % NodeIndexes, GotIt)
     
     AscalarImag(1:n) = ListGetReal( Material, &
          'Ascalar 2', n, Element % NodeIndexes, GotIt)
!
!    Source term:
!    ------------
     k = ListGetInteger( &
         Model % Bodies(Element % BodyId) % Values,'Body Force',GotIt, &
                    1, Model % NumberOFBodyForces )

     NodalSource = 0.0d0
     IF ( GotIt .AND. k > 0  ) THEN
        NodalSource(1,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
             'Pressure Source 1', n, Element % NodeIndexes, GotIt )

        NodalSource(2,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
             'Pressure Source 2', n, Element % NodeIndexes, GotIt )
     END IF
!
!    Integrate square of residual over element:
!    ------------------------------------------

     ResidualNorm = 0.0d0
     Area = 0.0d0

     IntegStuff = GaussPoints( Element )

     DO t=1,IntegStuff % n
        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
            Basis, dBasisdx, ddBasisddx, .TRUE. )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ

        ELSE
           u = SUM( Basis(1:n) * Nodes % x(1:n) )
           v = SUM( Basis(1:n) * Nodes % y(1:n) )
           w = SUM( Basis(1:n) * Nodes % z(1:n) )

           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                       Symb, dSymb, u, v, w )
           s = IntegStuff % s(t) * detJ * SqrtMetric

        END IF

        ResidualReal = 0.0d0
        ResidualImag = 0.0d0

        A2 = 0.0d0
        A1r = 0.0d0
        A1i = 0.0d0
        A0r = 0.0d0
        A0i = 0.0d0
        
        A0r = SUM( AscalarReal(1:n) * Basis(1:n) )
        A0i = SUM( AscalarImag(1:n) * Basis(1:n) )
        do i = 1,dim
           A1r(i) = A1r(i) + SUM( AvectorReal(i,1:n) * Basis(1:n) )
           A1i(i) = A1i(i) + SUM( AvectorImag(i,1:n) * Basis(1:n) )
           do j = 1,dim
              A2(i,j) = A2(i,j) + SUM( Amatrix(i,j,1:n) * Basis(1:n) )
           end do
        end do

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN

! diffusion
!-----------
           do i = 1,dim
              do j = 1,dim
                 ResidualReal = ResidualReal &
                      - SUM( Pressione(1,1:n) * ddBasisddx(1:n,i,j) ) * A2(i,j)

                 ResidualImag = ResidualImag &
                      - SUM( Pressione(2,1:n) * ddBasisddx(1:n,i,j) ) * A2(i,j)
              end do
           end do

! convection
!------------
           do i = 1,dim
              ResidualReal = ResidualReal + &
                   SUM( dBasisdx(1:n,i) * Pressione(1,1:n) ) * A1r(i) &
                   -SUM( dBasisdx(1:n,i) * Pressione(2,1:n) ) * A1i(i)
              
              ResidualImag = ResidualImag + &
                   SUM( dBasisdx(1:n,i) * Pressione(1,1:n) ) * A1i(i) &
                   +SUM( dBasisdx(1:n,i) * Pressione(2,1:n) ) * A1r(i)
           end do

! reaction
!----------
           ResidualReal = ResidualReal + &
                SUM( Basis(1:n) * Pressione(1,1:n) ) * A0r &
                -SUM( Basis(1:n) * Pressione(2,1:n) ) * A0i
           
           ResidualImag = ResidualImag + &
                SUM( Basis(1:n) * Pressione(1,1:n) ) * A0i &
                +SUM( Basis(1:n) * Pressione(2,1:n) ) * A0r

        ELSE
           print *,'Only cartesian coordinates implemented at the moment!'
           stop

!           DO j=1,DIM
!              DO k=1,DIM
!
!                - g^{jk} C_{,k}T_{j}:
!                ---------------------
!
!                 Residual = Residual - Metric(j,k) * &
!                    SUM( Temperature(1:n) * dBasisdx(1:n,j) ) * &
!                    SUM( NodalConductivity(1:n) * dBasisdx(1:n,k) )

!
!                - g^{jk} C T_{,jk}:
!                -------------------
!
!                 Residual = Residual - Metric(j,k) * Conductivity * &
!                    SUM( Temperature(1:n) * ddBasisddx(1:n,j,k) )
!
!                + g^{jk} C {_jk^l} T_{,l}:
!                ---------------------------
!                 DO l=1,DIM
!                    Residual = Residual + Metric(j,k) * Conductivity * &
!                      Symb(j,k,l) * SUM( Temperature(1:n) * dBasisdx(1:n,l) )
!                 END DO
!              END DO
!           END DO
        END IF

!
!       Compute also force norm for scaling the residual:
!       -------------------------------------------------
        DO i=1,DIM
           Fnorm = Fnorm &
                + s * ( SUM( NodalSource(1,1:n) * Basis(1:n) ) ) ** 2 &
                + s * ( SUM( NodalSource(2,1:n) * Basis(1:n) ) ) ** 2
        END DO

        Area = Area + s
        ResidualNorm = ResidualNorm + s *  ( ResidualReal ** 2 + ResidualImag ** 2 )
     END DO

     Fnorm = Element % hk**2 * Fnorm
     Indicator = Element % hK**2 * ResidualNorm

     DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )
     DEALLOCATE( NodalSource, Basis, dBasisdx, ddBasisddx, Pressione, &
       AMatrix, AvectorReal, AvectorImag,AscalarReal, AscalarImag )

CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE InputTensor( Tensor, IsScalar, Name, Material, n, NodeIndexes )
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Tensor(:,:,:)
      INTEGER :: n, NodeIndexes(:)
      LOGICAL :: IsScalar
      CHARACTER(LEN=*) :: Name
      TYPE(ValueList_t), POINTER :: Material
!------------------------------------------------------------------------------
      LOGICAL :: FirstTime = .TRUE., stat
      REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

      SAVE FirstTime, Hwrk
!------------------------------------------------------------------------------
      IF ( FirstTime ) THEN
         NULLIFY( Hwrk )
         FirstTime = .FALSE.
      END IF

      Tensor = 0.0d0

      CALL ListGetRealArray( Material, Name, Hwrk, n, NodeIndexes, stat )
      IsScalar = SIZE(HWrk,1) == 1 .AND. SIZE(HWrk,2) == 1

      IF ( .NOT. stat ) RETURN

      IF ( SIZE(Hwrk,1) == 1 ) THEN

         DO i=1,MIN(3,SIZE(Hwrk,2))
            Tensor( i,i,1:n ) = Hwrk( 1,1,1:n )
         END DO

      ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN

         DO i=1,MIN(3,SIZE(Hwrk,1))
            Tensor(i,i,1:n) = Hwrk(i,1,1:n)
         END DO

      ELSE

        DO i=1,MIN(3,SIZE(Hwrk,1))
           DO j=1,MIN(3,SIZE(Hwrk,2))
              Tensor( i,j,1:n ) = Hwrk(i,j,1:n)
           END DO
        END DO

      END IF
!------------------------------------------------------------------------------
   END SUBROUTINE InputTensor
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE InputVector( Tensor, IsScalar, Name, Material, n, NodeIndexes )
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Tensor(:,:)
      INTEGER :: n, NodeIndexes(:)
      LOGICAL :: IsScalar
      CHARACTER(LEN=*) :: Name
      TYPE(ValueList_t), POINTER :: Material
!------------------------------------------------------------------------------
      LOGICAL :: FirstTime = .TRUE., stat
      REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

      SAVE FirstTime, Hwrk
!------------------------------------------------------------------------------
      IF ( FirstTime ) THEN
         NULLIFY( Hwrk )
         FirstTime = .FALSE.
      END IF

      Tensor = 0.0d0

      CALL ListGetRealArray( Material, Name, Hwrk, n, NodeIndexes, stat )
      IsScalar = SIZE(HWrk,1) == 1 .AND. SIZE(HWrk,2) == 1

      IF ( .NOT. stat ) RETURN

      IF ( SIZE(Hwrk,1) == 1 ) THEN

         DO i=1,MIN(3,SIZE(Hwrk,2))
            Tensor( i,1:n ) = Hwrk( 1,1,1:n )
         END DO

      ELSE

        DO i=1,MIN(3,SIZE(Hwrk,1))
           Tensor( i,1:n ) = Hwrk( i,1,1:n )
        END DO

      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE InputVector
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   END FUNCTION DCRInsideResidual
!------------------------------------------------------------------------------


