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
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Ville Savolainen
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
!>  A solver for the Fourier transformed time-dependent diffusion equation
!>  with a constant frequency forcing term,
!>  (nabla^2 - iw) B(x) = -f(x,w)
! 
!>  Here we are solving specifically the magnetic problem
!>  (nabla^2/(mu*sigma) - iw) B(x) = iw B_ac(x)
!  
!>  This version solves for a vector field, i.e., the AC field is
!>  B_ac(x,y,z,t) = B_ac(x,y,z) exp(iwt) is a vector
! 
!>  If requested, the time-averaged Lorentz force due to induced field B,
!>  < 1/mu (curl B) x B >, is calculated and written in the result files.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE FourierDiffusion3DSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE Types
  USE Lists
  USE Integration
  USE ElementDescription
  USE SolverUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(Variable_t), POINTER :: Var
  TYPE(Nodes_t) :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement

  INTEGER, POINTER :: NodeIndexes(:)

  LOGICAL :: AllocationsDone = .FALSE., Bubbles, GotIt

  INTEGER, POINTER :: SomeQuantityPerm(:)
  REAL(KIND=dp), POINTER :: SomeQuantity(:), ForceVector(:), &
       LrF(:), LrFr(:), LrFz(:), LrFp(:), &
       BRer(:), BImr(:), BRez(:), BImz(:), BRep(:), BImp(:)

  INTEGER :: iter, i, j, k, n, t, istat, eq, LocalNodes
  REAL(KIND=dp) :: Norm, RelativeChange, AngularFrequency

  TYPE(ValueList_t), POINTER :: Material

  INTEGER :: NonlinearIter
  REAL(KIND=dp) :: NonlinearTol,s

  REAL(KIND=dp), ALLOCATABLE :: LocalStiffMatrix(:,:), Load(:,:), Work(:), &
       LocalForce(:), Impedance(:,:), Conductivity(:), Permeability(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: EquationName

  SAVE LocalStiffMatrix, Work, Load, LocalForce, ElementNodes, &
       Impedance, AllocationsDone, Conductivity, Permeability, &
       LrF, LocalNodes

#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at,at0,totat,st,totst,t1
#else
  REAL(KIND=dp) :: at,at0,totat,st,totst,t1,CPUTime,RealTime
#endif

  TYPE(Solver_t), POINTER :: SolverPointer

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT.ASSOCIATED( Solver % Matrix ) ) RETURN

  Solver % Matrix % Complex = .TRUE.

  SomeQuantity     => Solver % Variable % Values
  SomeQuantityPerm => Solver % Variable % Perm

  LocalNodes = COUNT( SomeQuantityPerm > 0 )
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
        DEALLOCATE( ElementNodes % x, &
             ElementNodes % y,      &
             ElementNodes % z,      &
             LocalForce,            &
             Impedance,             &
             Work,                  &
             LocalStiffMatrix,      &
             Conductivity, Permeability, &
             Load )
     END IF

     ALLOCATE( ElementNodes % x( N ),  &
          ElementNodes % y( N ),       &
          ElementNodes % z( N ),       &
          LocalForce( 6*N ),           &
          Impedance( 2,N ),            &
          Work( N ),                   &
          LocalStiffMatrix( 6*N,6*N ), &
          Conductivity( N ),  &
          Permeability( N ),  &
          LrF(3 * Model % NumberOfNodes), & 
          Load( 6,N ), STAT=istat )

     IF ( istat /= 0 ) THEN
        CALL Fatal( 'FourierDiffusion3DSolve', 'Memory allocation error.' )
     END IF

! Add Lorentz force to *.result and *.ep files
!#if 1
     Var => VariableGet(Model % Variables, 'Lorentz Force')
     IF (.NOT. ASSOCIATED(Var) ) THEN
        SolverPointer => Solver
        CALL VariableAddVector(Solver % Mesh % Variables, Solver % Mesh, &
              SolverPointer, 'Lorentz Force', 3, LrF, SomeQuantityPerm)
      END IF
!#endif

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
  IF ( .NOT.(GotIt) ) Bubbles = .FALSE.

!------------------------------
! Figure out angular frequency:
!------------------------------
  AngularFrequency = GetAngularFrequency()

!------------------------------------------------------------------------------
! Iterate over any nonlinearity of material or source
!------------------------------------------------------------------------------
  totat = 0.0d0
  totst = 0.0d0

  DO iter=1,NonlinearIter
!------------------------------------------------------------------------------
      at  = CPUTime()
      at0 = RealTime()

      CALL Info( 'FourierDiffusion3DSolve', ' ', Level=4 )
      CALL Info( 'FourierDiffusion3DSolve', ' ', Level=4 )
      CALL Info( 'FourierDiffusion3DSolve', &
            '-------------------------------------', Level=4 )
      WRITE( Message, * ) 'Fourier Diffusion 3D iteration', iter
      CALL Info( 'FourierDiffusion3DSolve', Message, Level=4 )
      WRITE( Message, * ) 'Frequency (Hz): ', AngularFrequency/(2*PI)
      CALL Info( 'FourierDiffusion3DSolve', Message, Level=4 )
      CALL Info( 'FourierDiffusion3DSolve', &
            '-------------------------------------', Level=4 )
      CALL Info( 'FourierDiffusion3DSolve', ' ', Level=4 )
      CALL Info( 'FourierDiffusion3DSolve', 'Starting Assembly...', Level=4 )

     CALL DefaultInitialize( )
!
!    Do the bulk assembly:
!    ---------------------

!------------------------------------------------------------------------------
     DO t=1,Solver % Mesh % NumberOfBulkElements
!------------------------------------------------------------------------------
        IF ( RealTime() - at0 > 1.0 ) THEN
          WRITE(*,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
           (Solver % Mesh % NumberOfBulkElements-t) / &
              (1.0*Solver % Mesh % NumberOfBulkElements)), ' % done'
                      
          at0 = RealTime()
        END IF
!------------------------------------------------------------------------------
!       Check if this element belongs to a body where this equation
!       should be computed
!------------------------------------------------------------------------------
        CurrentElement => Solver % Mesh % Elements(t)

        IF ( .NOT. CheckElementEquation( Model, &
             CurrentElement, EquationName ) ) CYCLE
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
             Bodyid ) % Values, 'Equation', minv=1, maxv=Model % NumberOFEquations )

        Work(1:n) = ListGetReal(  Model % Equations(k) % Values, &
             'Angular Frequency', n, NodeIndexes, GotIt )

        IF ( GotIt ) THEN
           AngularFrequency = Work(1)
        ELSE
           Work(1:n) = ListGetReal(  Model % Equations(k) % Values, &
                 'Frequency', n, NodeIndexes, GotIt )
           IF ( GotIt ) AngularFrequency = 2*PI*Work(1)
        END IF

        k = ListGetInteger( Model % Bodies( CurrentElement % &
           Bodyid ) % Values, 'Material', minv=1, maxv=Model % NumberOfMaterials )

        Material => Model % Materials(k) % Values

        Conductivity(1:n) = ListGetReal(Material, &
             'Electrical Conductivity',n,NodeIndexes,GotIt)
        IF(GotIt) THEN
          CALL Warn('Fdiffusion3D','Use electric conductivity instead of electrical')
        END IF
 
        Conductivity(1:n) = ListGetReal(Material, &
             'Electric Conductivity',n,NodeIndexes)

        Permeability(1:n) = ListGetReal(Material, &
             'Magnetic Permeability',n,NodeIndexes)

!------------------------------------------------------------------------------
!       The source term on nodes
!------------------------------------------------------------------------------
        k = ListGetInteger( Model % Bodies( CurrentElement % BodyId ) % &
               Values, 'Body Force', GotIt, 1, Model % NumberOFBodyForces )

        Load = 0.0d0
        IF ( k > 0 ) THEN
           Load(1,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
                'AC Magnetic Field 1 Re', n, NodeIndexes, GotIt )

           Load(2,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
                'AC Magnetic Field 1 Im', n, NodeIndexes, GotIt )

           Load(3,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
                'AC Magnetic Field 2 Re', n, NodeIndexes, GotIt )

           Load(4,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
                'AC Magnetic Field 2 Im', n, NodeIndexes, GotIt )

           Load(5,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
                'AC Magnetic Field 3 Re', n, NodeIndexes, GotIt )

           Load(6,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
                'AC Magnetic Field 3 Im', n, NodeIndexes, GotIt )

        END IF

!------------------------------------------------------------------------------
!       Get element local matrix and rhs vector
!------------------------------------------------------------------------------
        CALL LocalMatrix(  LocalStiffMatrix, LocalForce, AngularFrequency, &
             Conductivity*Permeability, Load, Bubbles, CurrentElement, n, &
             ElementNodes )

!------------------------------------------------------------------------------
!       Update global matrix and rhs vector from local matrix & vector
!------------------------------------------------------------------------------
        CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
             ForceVector, LocalForce, n, Solver % Variable % DOFs, &
                  SomeQuantityPerm(NodeIndexes) )
!------------------------------------------------------------------------------
     END DO
!------------------------------------------------------------------------------

!
!    Neumann & Newton BCs:
!    ---------------------

! At the moment (for the magnetic problem),
! I can't think of any other BCs except Dirichlet and natural conditions.
! The scalar equation is handled in FDiffusion (LocalMatrixBoundary),
! and it (indices etc.) should be changed as LocalMatrix has been.

     CALL DefaultFinishAssembly( )
!
!    Dirichlet BCs:
!    --------------
     CALL DefaultDirichletBCs()

     CALL Info( 'FourierDiffusion3DSolve', 'Assembly done', Level=4 )

!
!    Solve the system and we are done:
!    ---------------------------------

     at = CPUTime() - at
     st = CPUTime()

     Norm = DefaultSolve()

     st = CPUTIme()-st
     totat = totat + at
     totst = totst + st
     WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
     CALL Info( 'FourierDiffusion3DSolve', Message, Level=4 )
     WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
     CALL Info( 'FourierDiffusion3DSolve', Message, Level=4 )

!------------------------------------------------------------------------------

     WRITE( Message, * ) 'Result Norm    : ',Norm
     CALL Info( 'FourierDiffusion3DSolve', Message, Level=4 )


     RelativeChange = Solver % Variable % NonlinChange
     WRITE( Message, * ) 'Relative Change: ',RelativeChange
     CALL Info( 'FourierDiffusion3DSolve', Message, Level=4 )

     IF( Solver % Variable % NonlinConverged == 1 ) EXIT
!------------------------------------------------------------------------------
  END DO ! of nonlinear iteration
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!   Compute the time-averaged Lorentz force
!------------------------------------------------------------------------------

!#if 1
    LrFr => LrF(1::3)
    LrFz => LrF(2::3)
    LrFp => LrF(3::3)

    BRer => SomeQuantity(1::6)
    BImr => SomeQuantity(2::6)
    BRez => SomeQuantity(3::6)
    BImz => SomeQuantity(4::6)
    BRep => SomeQuantity(5::6)
    BImp => SomeQuantity(6::6)

    CALL LorentzForceAve(LrFr,LrFz,LrFp,BRer,BImr,BRez,BImz,BRep,BImp,&
         SomeQuantityPerm)
!#endif

CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  StiffMatrix, Force, AngularFrequency, &
       NodalConductivity, Load, Bubbles, Element, n, Nodes )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: StiffMatrix(:,:), Force(:), AngularFrequency, &
         NodalConductivity(:), Load(:,:)
    LOGICAL :: Bubbles
    INTEGER :: n
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(2*n),dBasisdx(2*n,3)
    REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,Conductivity,M,D,LRe(3),LIm(3)
    COMPLEX(KIND=dp) :: LSTIFF(3*n,3*n), LFORCE(3*n), A(3,3)
    LOGICAL :: Stat
    INTEGER :: i,j,k,l,h,p,q,t,dim, NBasis, CoordSys
    TYPE(GaussIntegrationPoints_t) :: IntegStuff

    REAL(KIND=dp) :: X,Y,Z,Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3)
!------------------------------------------------------------------------------
    dim = 3
    CoordSys = CurrentCoordinateSystem()

    Metric = 0.0d0
    Metric(1,1) = 1.0d0
    Metric(2,2) = 1.0d0
    Metric(3,3) = 1.0d0

    LSTIFF = 0.0d0
    LFORCE = 0.0d0
    A = 0.0d0

!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    IF ( Bubbles ) THEN
       PRINT*,'FourierDiffusion3DSolver: LocalMatrix: Cannot handle bubbles.'
       STOP
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
!      Effective conductivity (times the magnetic permeability)
!------------------------------------------------------------------------------
       Conductivity = SUM( NodalConductivity(1:n)*Basis )

       D = -AngularFrequency
       M = 0._dp

       DO i=1,dim
          LRe(i) = AngularFrequency * SUM( Load(2*i-1,1:n) * Basis(1:n) )
          LIm(i) = AngularFrequency * SUM( Load(2*i,1:n) * Basis(1:n) )
       END DO
!------------------------------------------------------------------------------
!      The Helmholz equation
!------------------------------------------------------------------------------
       DO p=1,NBasis
          DO q=1,NBasis

             DO i=1,dim
                A(i,i) = CMPLX( M, D, KIND=dp ) * Basis(q) * Basis(p)
             END DO

             DO i=1,dim
                DO j=1,dim
                   DO k=1,dim
                      A(i,i) = A(i,i) + Metric(j,k) * dBasisdx(q,k) *  dBasisdx(p,j) / Conductivity
                      DO l=1,dim
                         A(l,i) = A(l,i) - Metric(j,k) * dBasisdx(q,k) * &
                              Symb(i,j,l) * Basis(p) / Conductivity
                         A(i,l) = A(i,l) + Metric(j,k) * Symb(l,k,i) * &
                              Basis(q) * dBasisdx(p,j) / Conductivity
                         DO h=1,dim
                            A(h,l) = A(h,l) - Metric(j,k) * Basis(q) * &
                                 Basis(p) * Symb(l,k,i) * Symb(i,j,h) / &
                                 Conductivity
                         END DO
                      END DO
                   END DO
                END DO
             END DO

             DO i=1,dim
                DO j=1,dim
                   LSTIFF(dim*(p-1)+i,dim*(q-1)+j) = LSTIFF(dim*(p-1)+i,dim*(q-1)+j) + s*A(i,j)
                END DO
             END DO

          END DO

          DO i=1,dim
             LFORCE(dim*(p-1)+i) = LFORCE(dim*(p-1)+i) + s * Basis(p) * CMPLX( LRe(i),LIm(i),KIND=dp )
          END DO

       END DO
    END DO
!------------------------------------------------------------------------------

    DO p=1,n
       DO i=1,dim
          Force( 2*dim*(p-1)+2*i-1 ) = REAL( LFORCE(dim*(p-1)+i) )
          Force( 2*dim*(p-1)+2*i ) = AIMAG( LFORCE(dim*(p-1)+i) )

          DO q=1,n
             DO j=1,dim

                StiffMatrix( 2*dim*(p-1)+2*i-1, 2*dim*(q-1)+2*j-1 ) = &
                     REAL( LSTIFF(dim*(p-1)+i,dim*(q-1)+j) )
                StiffMatrix( 2*dim*(p-1)+2*i-1, 2*dim*(q-1)+2*j ) = &
                     -AIMAG( LSTIFF(dim*(p-1)+i,dim*(q-1)+j) )
                StiffMatrix( 2*dim*(p-1)+2*i, 2*dim*(q-1)+2*j-1 ) = &
                     AIMAG( LSTIFF(dim*(p-1)+i,dim*(q-1)+j) )
                StiffMatrix( 2*dim*(p-1)+2*i, 2*dim*(q-1)+2*j ) = &
                     REAL( LSTIFF(dim*(p-1)+i,dim*(q-1)+j) )
             END DO
          END DO

       END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBoundary(  StiffMatrix, Force, AngularFrequency, &
              Impedance, Load, Element, n, Nodes )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: StiffMatrix(:,:),Force(:),Impedance(:,:),Load(:,:)
    REAL(KIND=dp) :: AngularFrequency
    INTEGER :: n
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,Impedance1,Impedance2,L1,L2
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),X,Y,Z
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
       IF ( CoordSys /= Cartesian ) THEN
          X = SUM( Nodes % x(1:n) * Basis(1:n) )
          Y = SUM( Nodes % y(1:n) * Basis(1:n) )
          Z = SUM( Nodes % z(1:n) * Basis(1:n) )
          s = s * CoordinateSqrtMetric( X,Y,Z )
       END IF
!------------------------------------------------------------------------------
       Impedance1 = SUM( Impedance(1,1:n) * Basis )
       IF ( ABS(Impedance1) > AEPS ) &
          Impedance1 = AngularFrequency / Impedance1

       Impedance2 = SUM( Impedance(2,1:n) * Basis )
       IF ( ABS(Impedance2) > AEPS ) &
          Impedance2 = AngularFrequency / Impedance2

       L1 = SUM( Load(1,1:n) * Basis )
       L2 = SUM( Load(2,1:n) * Basis )
!------------------------------------------------------------------------------
       DO p=1,n
          DO q=1,n
             A = CMPLX(Impedance1, Impedance2,KIND=dp) * Basis(q) * Basis(p)
             LSTIFF(p,q) = LSTIFF(p,q) + s * A
          END DO

          LFORCE(p) = LFORCE(p) + s * Basis(p) * CMPLX(L1,L2,KIND=dp)
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
   SUBROUTINE LorentzForceAve( LrFx,LrFy,LrFz,BRex,BImx,BRey,BImy,BRez,BImz,&
        Reorder )
!------------------------------------------------------------------------------

! Calculate the nodal, time-averaged, values of the Lorentz Force.
!
! Both induced field solved through the Fourier transform
! and the external AC field may have all components in any coordinate system

     USE Types
     IMPLICIT NONE
     REAL(KIND=dp) :: BRex(:),BImx(:),BRey(:),BImy(:),BRez(:),BImz(:)
     REAL(KIND=dp) :: LrFx(:),LrFy(:),LrFz(:), Lorentz(3)
     INTEGER :: Reorder(:)

     TYPE(Element_t), POINTER :: Element
     TYPE(Nodes_t) :: Nodes 

     LOGICAL :: Stat

     INTEGER, POINTER :: NodeIndexes(:),Visited(:)
     INTEGER :: p,q,i,t,n,k

     REAL(KIND=dp) :: u,v,w,x,y,z

     REAL(KIND=dp) :: B(3),dHdx(3,3),L(3)
     REAL(KIND=dp) :: SqrtElementMetric
     REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)

     REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
     REAL(KIND=dp), ALLOCATABLE :: Permeability(:), BacRe(:,:), BacIm(:,:)
!------------------------------------------

     ALLOCATE( Visited(CurrentModel % NumberOfNodes) )

     n = CurrentModel % Mesh % MaxElementDOFs
     ALLOCATE(Nodes % x(n),Nodes % y(n),Nodes % z(n))
     ALLOCATE( BacRe(3,n), BacIm(3,n) )
     ALLOCATE( Permeability(n), Basis(n), dBasisdx(n,3) )

     Visited = 0

     LrFx = 0._dp
     LrFy = 0._dp
     LrFz = 0._dp

     DO t=1,CurrentModel % NumberOfBulkElements

        Element => CurrentModel % Elements(t)
        n = Element % TYPE % NumberOfNodes
        NodeIndexes => Element % NodeIndexes

        Nodes % x(1:n) = CurrentModel % Nodes % x( NodeIndexes )
        Nodes % y(1:n) = CurrentModel % Nodes % y( NodeIndexes )
        Nodes % z(1:n) = CurrentModel % Nodes % z( NodeIndexes )

        Permeability(1:n) = ListGetReal(Material, &
             'Magnetic Permeability',n,NodeIndexes)

        k = ListGetInteger( Model % Bodies( Element % BodyId ) % &
             Values, 'Body Force', GotIt,1,Model % NumberOFBodyForces )

        BacRe = 0._dp
        BacIm = 0._dp

        IF ( k > 0 ) THEN

           BacRe(1,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
                'AC Magnetic Field 1 Im', n, NodeIndexes, GotIt )

           BacIm(1,1:n) = -ListGetReal( Model % BodyForces(k) % Values, &
                'AC Magnetic Field 1 Re', n, NodeIndexes, GotIt )

           BacRe(2,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
                'AC Magnetic Field 2 Im', n, NodeIndexes, GotIt )

           BacIm(2,1:n) = -ListGetReal( Model % BodyForces(k) % Values, &
                'AC Magnetic Field 2 Re', n, NodeIndexes, GotIt )

           BacRe(3,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
                'AC Magnetic Field 3 Im', n, NodeIndexes, GotIt )

           BacIm(3,1:n) = -ListGetReal( Model % BodyForces(k) % Values, &
                'AC Magnetic Field 3 Re', n, NodeIndexes, GotIt )

! The apparent inconsistency is due to the fact that the sif-file keyword
! actually refers to the time-derivated load with i (but not omega!).
! For calculating the average Lorentz force, we really need the phase of
! B^ac defined consistently with B^i.
        END IF


        IF ( MINVAL(Reorder(NodeIndexes)) > 0 ) THEN

           DO p=1,n

              B = 0._dp
              dHdx = 0._dp

              q = Reorder(NodeIndexes(p))
              u = Element % TYPE % NodeU(p)
              v = Element % TYPE % NodeV(p)

              IF ( Element % TYPE % DIMENSION == 3 ) THEN
                 w = Element % TYPE % NodeW(p)
              ELSE
                 w = 0.0D0
              END IF

              stat = ElementInfo( Element, Nodes, u, v, w, SqrtElementMetric, &
                   Basis, dBasisdx )
              x = SUM( Nodes % x(1:n) * Basis(1:n) )
              y = SUM( Nodes % y(1:n) * Basis(1:n) )
              z = SUM( Nodes % z(1:n) * Basis(1:n) )
              CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
              IF ( CurrentCoordinateSystem() /= Cartesian ) CALL InvertMatrix( Metric,3 )

! Call ComputeLorentz
! Contribution from real parts
              L = 0._dp

              B(1) = SUM( Basis(1:n)*BRex(Reorder(NodeIndexes)) ) + BacRe(1,p)
              B(2) = SUM( Basis(1:n)*BRey(Reorder(NodeIndexes)) ) + BacRe(2,p)
              B(3) = SUM( Basis(1:n)*BRez(Reorder(NodeIndexes)) ) + BacRe(3,p)

              DO i=1,3
                 dHdx(1,i) = SUM( dBasisdx(1:n,i)* &
                      BRex(Reorder(NodeIndexes)) ) / Permeability(p)
                 dHdx(2,i) = SUM( dBasisdx(1:n,i)* &
                      BRey(Reorder(NodeIndexes)) ) / Permeability(p)
                 dHdx(3,i) = SUM( dBasisdx(1:n,i)* &
                      BRez(Reorder(NodeIndexes)) ) / Permeability(p)
              END DO

              L = ComputeLorentz( B,dHdx,Permeability(p), &
                   SqrtMetric,Metric,Symb )
! Contribution from imaginary parts
              B(1) = SUM( Basis(1:n)*BImx(Reorder(NodeIndexes)) ) + BacIm(1,p)
              B(2) = SUM( Basis(1:n)*BImy(Reorder(NodeIndexes)) ) + BacIm(2,p)
              B(3) = SUM( Basis(1:n)*BImz(Reorder(NodeIndexes)) ) + BacIm(3,p)

              DO i=1,3
                 dHdx(1,i) = SUM( dBasisdx(1:n,i)* &
                      BImx(Reorder(NodeIndexes)) ) / Permeability(p)
                 dHdx(2,i) = SUM( dBasisdx(1:n,i)* &
                      BImy(Reorder(NodeIndexes)) ) / Permeability(p)
                 dHdx(3,i) = SUM( dBasisdx(1:n,i)* &
                      BImz(Reorder(NodeIndexes)) ) / Permeability(p)
              END DO

              L = L + ComputeLorentz( B,dHdx,Permeability(p), &
                   SqrtMetric,Metric,Symb )
! Add to the averaged Lorentz force
              LrFx(q) = LrFx(q) + 0.5_dp * L(1)
              LrFy(q) = LrFy(q) + 0.5_dp * L(2)
              IF ( CurrentCoordinateSystem() == CylindricSymmetric ) THEN
!#if 0
!! Convert the azimuthal force from SI Units (if it exists)
!                 IF (x > 1.0d-10) THEN
!                    L(3) = L(3) / x
!                 ELSE
!                    L(3) = 0._dp
!                 END IF
!#else
! If AC field doesn't have phi-component, then Lorentz force doesn't have it.
                 L(3) = 0._dp
!#endif
              END IF
              LrFz(q) = LrFz(q) + 0.5_dp * L(3)

              Visited(q) = Visited(q) + 1
           
           END DO
        END IF
      END DO

      DO i=1,CurrentModel % NumberOfNodes
         IF ( Visited(i) > 1 ) THEN
            LrFx(i) = LrFx(i) / Visited(i)
            LrFy(i) = LrFy(i) / Visited(i)
            LrFz(i) = LrFz(i) / Visited(i)
         END IF
      END DO

      DEALLOCATE( Visited )
      DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )
      DEALLOCATE ( Permeability )
      DEALLOCATE( BacRe, BacIm, Basis, dBasisdx )

    END SUBROUTINE LorentzForceAve

!------------------------------------------------------------------------------
    FUNCTION ComputeLorentz( B,dHdx,mu,SqrtMetric,Metric,Symb ) RESULT(LF)
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: B(:),dHdx(:,:),mu,LF(3),SqrtMetric,Metric(:,:), &
           Symb(:,:,:)
!------------------------------------------------------------------------------
      INTEGER :: i,j,k,l,m
      REAL(KIND=dp) :: Bc(3),Ji(3),Jc(3),s,Perm(3,3,3),r
!------------------------------------------------------------------------------

      IF ( CurrentCoordinateSystem() == Cartesian ) THEN
         Ji(1) = dHdx(3,2) - dHdx(2,3)
         Ji(2) = dHdx(1,3) - dHdx(3,1)
         Ji(3) = dHdx(2,1) - dHdx(1,2)
         LF(1) = Ji(2)*B(3) - Ji(3)*B(2)
         LF(2) = Ji(3)*B(1) - Ji(1)*B(3)
         LF(3) = Ji(1)*B(2) - Ji(2)*B(1)
         RETURN
      END IF

      r = SqrtMetric

      IF ( CurrentCoordinateSystem() == CylindricSymmetric ) THEN
         Ji(1) = -dHdx(3,2)
         Ji(2) =  dHdx(3,1)
         IF (r > 1.0d-10) THEN
            Ji(2) = Ji(2) + B(3)/(r*mu)
         ELSE
            Ji(2) = Ji(2) + Ji(2)
         END IF
         Ji(3) = dHdx(1,2) - dHdx(2,1)

         LF(1) = Ji(3)*B(2) - Ji(2)*B(3)
         LF(2) = Ji(1)*B(3) - Ji(3)*B(1)
! Use SI units for the azimuthal component or the Lorentz force,
! since it is computed at nodal points and possibly at symmetry axis,
! otherwise you divide by zero.
!#if 1
         LF(3) = Ji(2)*B(1) - Ji(1)*B(2)
!#else
!         IF (r > 1.0d-10) THEN
!            LF(3) = ( Ji(2)*B(1) - Ji(1)*B(2) ) / r
!         ELSE
!            LF(3) = 0.d0
!         END IF
!#endif
         RETURN
      END IF

      Perm = 0
      Perm(1,2,3) = -1.0d0 / SqrtMetric
      Perm(1,3,2) =  1.0d0 / SqrtMetric
      Perm(2,1,3) =  1.0d0 / SqrtMetric
      Perm(2,3,1) = -1.0d0 / SqrtMetric
      Perm(3,1,2) = -1.0d0 / SqrtMetric
      Perm(3,2,1) =  1.0d0 / SqrtMetric
!------------------------------------------------------------------------------

      Bc = 0.0d0
      DO i=1,3
         DO j=1,3
            Bc(i) = Bc(i) + Metric(i,j)*B(j)
         END DO
      END DO

!------------------------------------------------------------------------------

      Ji = 0.0d0
      DO i=1,3
         s = 0.0D0
         DO j=1,3
            DO k=1,3
               IF ( Perm(i,j,k) /= 0 ) THEN
                  DO l=1,3
                     s = s + Perm(i,j,k)*Metric(j,l)*dHdx(l,k)
                     DO m=1,3
                        s = s + Perm(i,j,k)*Metric(j,l)*Symb(k,m,l)*B(m)/mu
                     END DO
                  END DO
               END IF
            END DO
         END DO
         Ji(i) = s
      END DO
 
      Jc = 0.0d0
      DO i=1,3
         DO j=1,3
            Jc(i) = Jc(i) + Metric(i,j)*Ji(j)
         END DO
      END DO
!------------------------------------------------------------------------------

      LF = 0.0d0
      DO i=1,3
         s = 0.0D0
         DO j=1,3
            DO k=1,3
               IF ( Perm(i,j,k) /= 0 ) THEN
                  s = s + Perm(i,j,k)*Jc(k)*Bc(j)
               END IF
            END DO
         END DO
         LF(i) = s
      END DO
!------------------------------------------------------------------------------
    END FUNCTION ComputeLorentz
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE FourierDiffusion3DSolver
!------------------------------------------------------------------------------
