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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 2002
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
   MODULE GlobMatC
!------------------------------------------------------------------------------
      USE Types
      COMPLEX(KIND=dp), ALLOCATABLE :: Matrix(:,:)
!------------------------------------------------------------------------------
   END MODULE GlobMatC
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Solves the Helmholtz equation using BEM!
!> This solver can only deal with rather small problems as it does not use any
!> multilevel strategies. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE HelmholtzBEMSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
     USE GlobMatC
     USE DefUtils

     IMPLICIT NONE
!------------------------------------------------------------------------------
 
     TYPE(Model_t) :: Model
     TYPE(Solver_t):: Solver
 
     REAL(KIND=dp) :: dt
     LOGICAL :: TransientSimulation
 
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,n,t,istat,bf_id,BoundaryNodes
 
     TYPE(Matrix_t),POINTER  :: STIFF
     TYPE(Nodes_t)   :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement
 
     REAL(KIND=dp) :: Norm, PrevNorm
     INTEGER, POINTER :: NodeIndexes(:)

     LOGICAL :: AllocationsDone = .FALSE., GotIt
 
     COMPLEX(KIND=dp), POINTER CONTIG :: Potential(:),ForceVector(:), &
               Diagonal(:)
     INTEGER, POINTER :: PotentialPerm(:), BoundaryPerm(:)

     LOGICAL, ALLOCATABLE :: PotentialKnown(:)
 
     COMPLEX(KIND=dp), ALLOCATABLE ::  Flx(:), Pot(:), VolumeForce(:), &
                       Load(:)

     REAL(KIND=dp), ALLOCATABLE :: P1(:), P2(:)
 
#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at,st,s, AngularFrequency, Work(1)
#else
     REAL(KIND=dp) :: at,st,CPUTime,s, AngularFrequency, Work(1)
#endif
     TYPE(Variable_t), POINTER :: Var

     CHARACTER(LEN=MAX_NAME_LEN) :: EquationName

     SAVE Load, ElementNodes, AllocationsDone, &
        PotentialKnown, PotentialPerm, BoundaryPerm, BoundaryNodes, &
           Pot, Flx, P1, P2, Potential, ForceVector, Diagonal

     EXTERNAL Matvec, Precond

     IF( CurrentCoordinateSystem() /= Cartesian ) THEN
       CALL Fatal('HelmholtzBEMSolver','This solver is implemented only for cartesian coordinates!')
     END IF

!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN
!------------------------------------------------------------------------------
!       Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
!
!       Get permutation of mesh nodes, so that boundary nodes get
!       numbered from 1..nb:
!       ---------------------------------------------------------
        ALLOCATE( PotentialPerm( Solver % Mesh % NumberOfNodes ), &
                   BoundaryPerm( Solver % Mesh % NumberOfNodes ), STAT=istat)

        IF ( istat /= 0 ) THEN
           CALL Fatal( 'HelmholtzBEMSolver', 'Memory allocation error 1.' )
        END IF

        PotentialPerm = 0
        BoundaryPerm  = 0
        BoundaryNodes = 0

        DO t=1,Solver % NumberOfActiveElements
           CurrentElement => Solver % Mesh % Elements( Solver % ActiveElements(t) )

           IF ( .NOT. ASSOCIATED( CurrentElement % BoundaryInfo ) ) CYCLE
           IF (CurrentElement % Type % ElementCode == 101 ) CYCLE

           DO j=1, CurrentElement % Type % NumberOfNodes
              k = CurrentElement % NodeIndexes(j)
              IF ( PotentialPerm(k) == 0 ) THEN
                 BoundaryNodes = BoundaryNodes + 1
                 BoundaryPerm(BoundaryNodes) = k
                 PotentialPerm(k) = BoundaryNodes
              END IF
           END DO
        END DO

        N = Model % MaxElementNodes
 
        ALLOCATE( ElementNodes % x( N ),                  &
                  ElementNodes % y( N ),                  &
                  ElementNodes % z( N ),                  &
                  Flx( BoundaryNodes ),                   &
                  Pot( BoundaryNodes ),                   &
                  Load( BoundaryNodes ),                  &
                  Diagonal( BoundaryNodes ),              &
                  PotentialKnown( BoundaryNodes ),        &
                  Matrix( BoundaryNodes, BoundaryNodes ), STAT=istat )

        IF ( istat /= 0 ) THEN
           CALL Fatal( 'HelmholtzBEMSolver', 'Memory allocation error 2.' )
        END IF

        ALLOCATE( Potential( Solver % Mesh % NumberOfNodes ), P1(n), P2(n), &
             ForceVector( Solver % Mesh % NumberOfNodes ),STAT=istat ) 

        IF ( istat /= 0 ) THEN
           CALL Fatal( 'HelmholtzBEMSolver', 'Memory allocation error 3.' )
        END IF
 
        AllocationsDone = .TRUE.
     END IF

!
!------------------------------------------------------------------------------
! Figure out angular frequency:
!------------------------------------------------------------------------------
     AngularFrequency = GetAngularFrequency()


!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------
     at = CPUTime()
     EquationName = ListGetString( Model % Solver % Values, 'Equation' )

     Matrix      = 0.0d0
     Load        = 0.0d0
     Diagonal    = 0.0d0
     ForceVector = 0.0d0

     DO i=1,Solver % Mesh % NumberOfNodes
        Potential(i) = CMPLX( Solver % Variable % Values(2*(i-1)+1), &
                              Solver % Variable % Values(2*(i-1)+2),KIND=dp  )
     END DO
!------------------------------------------------------------------------------
!    Check the bndry conditions. For each node either flux or potential must be
!    given and the other is the unknown! After the loop a logical variable
!    PotentialKnown will be true if flux is the unknown and false if potential
!    is the unknown for each node. Also vector Load will contain the known value
!    for each node.
!------------------------------------------------------------------------------
     DO t=1,Solver % NumberOfActiveElements
        CurrentElement => Solver % Mesh % Elements( Solver % ActiveElements(t) )
        IF ( .NOT. ASSOCIATED( CurrentElement % BoundaryInfo ) ) CYCLE
        IF ( CurrentElement % Type % ElementCode == 101 ) CYCLE

        n = CurrentElement % Type % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes

        DO i=1,Model % NumberOfBCs
          IF ( CurrentElement % BoundaryInfo % Constraint /= Model % BCs(i) % Tag ) CYCLE

          
          P1(1:n) = ListGetReal( Model % BCs(i) % Values, &
               ComponentName(Solver % Variable,1), n, NodeIndexes, GotIt )

          IF ( .NOT. GotIt ) THEN
             P1(1:n) = ListGetReal( Model % BCs(i) % Values, 'Potential 1' , n, NodeIndexes, GotIt )
          END IF

          P2(1:n) = ListGetReal( Model % BCs(i) % Values, &
               ComponentName(Solver % Variable,2), n, NodeIndexes, GotIt )

          IF ( .NOT. GotIt ) THEN
             P2(1:n) = ListGetReal( Model % BCs(i) % Values, 'Potential 2' , n, NodeIndexes, GotIt )
          END IF

          Load( PotentialPerm(NodeIndexes) ) = CMPLX( P1(1:n),P2(1:n),KIND=dp )

          IF ( .NOT. GotIt ) THEN
             PotentialKnown( PotentialPerm(NodeIndexes) ) = .FALSE.

             P1(1:n) = ListGetReal( Model % BCs(i) % Values, 'Flux 1', n, NodeIndexes, GotIt )
             P2(1:n) = ListGetReal( Model % BCs(i) % Values, 'Flux 2', n, NodeIndexes, GotIt )
             Load( PotentialPerm(NodeIndexes) ) = CMPLX( P1(1:n),P2(1:n),KIND=dp )
          ELSE
             PotentialKnown( PotentialPerm(NodeIndexes) ) = .TRUE.
          END IF
          EXIT
        END DO
     END DO
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
!    Matrix assembly:
!    ----------------
     DO t=1,Solver % NumberOfActiveElements
        CurrentElement => Solver % Mesh % Elements( Solver % ActiveElements(t) )

        IF ( .NOT. ASSOCIATED( CurrentElement % BoundaryInfo ) ) CYCLE
        IF ( CurrentElement % Type % ElementCode == 101 ) CYCLE

        n = CurrentElement % Type % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes
 
        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x( NodeIndexes )
        ElementNodes % y(1:n) = Solver % Mesh % Nodes % y( NodeIndexes )
        ElementNodes % z(1:n) = Solver % Mesh % Nodes % z( NodeIndexes )

        CALL IntegrateMatrix( Matrix, Diagonal, ForceVector, Load, &
             PotentialKnown, CurrentElement, n, ElementNodes )
     END DO

!------------------------------------------------------------------------------
     DO i=1,BoundaryNodes
        IF ( PotentialKnown(i) ) THEN
           ForceVector(i) = ForceVector(i) - Load(i)*Diagonal(i)
        ELSE
           Matrix(i,i) = Diagonal(i)
        END IF
     END DO

!------------------------------------------------------------------------------

     at = CPUTime() - at
     PRINT*,'Assembly (s): ',at

!------------------------------------------------------------------------------
!    Solve the system and we are done.
!------------------------------------------------------------------------------
     st = CPUTime()
!
!    Solve system:
!    -------------
     CALL SolveFull( BoundaryNodes, Matrix, Potential, ForceVector, Solver )
!
!    Extract potential and fluxes for the boundary nodes:
!    ----------------------------------------------------
     DO i=1,BoundaryNodes
        IF ( PotentialKnown(i) ) THEN
           Flx(i) = Potential(i)
           Pot(i) = Load(i)
        ELSE
           Flx(i) = Load(i)
           Pot(i) = Potential(i)
        END IF
     END DO
     st = CPUTime() - st
     PRINT*,'Solve (s):    ',st
!
     st = CPUTime()
!    Now compute potential for all mesh points:
!    ------------------------------------------
     Potential = 0.0d0
     DO i=1,BoundaryNodes
        Potential(BoundaryPerm(i)) = Pot(i)
     END DO

     DO t=1,Solver % NumberOfActiveElements
        CurrentElement => Solver % Mesh % Elements( Solver % ActiveElements(t) )
        IF ( .NOT. ASSOCIATED( CurrentElement % BoundaryInfo ) ) CYCLE

        IF ( CurrentElement % Type % ElementCode == 101 ) CYCLE

        n = CurrentElement % Type % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes

        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x( NodeIndexes )
        ElementNodes % y(1:n) = Solver % Mesh % Nodes % y( NodeIndexes )
        ElementNodes % z(1:n) = Solver % Mesh % Nodes % z( NodeIndexes )

        CALL ComputePotential( Potential, Pot, Flx, CurrentElement, n, ElementNodes )
     END DO

     Solver % Variable % Values = 0.0d0
     DO i=1,Solver % Mesh % NumberOfNodes
        j = Solver % Variable % Perm(i)
        IF ( j > 0 ) THEN
           Solver % Variable % Values( 2*(j-1) + 1 ) =  REAL( Potential(i) )
           Solver % Variable % Values( 2*(j-1) + 2 ) = AIMAG( Potential(i) )
        END IF
     END DO

     Var => VariableGet( Solver % Mesh % Variables, 'Flux' )
     IF ( ASSOCIATED( Var ) ) THEN
        Var % Values = 0.0d0
        DO i=1,BoundaryNodes
           k = BoundaryPerm(i)
           Var % Values(2*(k-1)+1) =  REAL( Flx(i) )
           Var % Values(2*(k-1)+2) = AIMAG( Flx(i) )
        END DO
     END IF
!
!    All done, finalize:
!    -------------------
     Solver % Variable % Norm = SQRT( SUM( ABS(Potential)**2 ) ) / &
                Solver % Mesh % NumberOfNodes 

     CALL InvalidateVariable( Model % Meshes, &
                  Solver % Mesh, Solver % Variable % Name )
!------------------------------------------------------------------------------
     st = CPUTime() - st
     PRINT*,'Post Processing (s):    ',st
!------------------------------------------------------------------------------
 

   CONTAINS


!------------------------------------------------------------------------------
     SUBROUTINE IntegrateMatrix( STIFF, ADiagonal, Force, Source, &
                    PotentialKnown, Element, n, Nodes )
!------------------------------------------------------------------------------
       COMPLEX(KIND=dp) :: STIFF(:,:), ADiagonal(:)
       COMPLEX(KIND=dp) :: FORCE(:), Source(:)
       INTEGER :: n
       LOGICAL :: PotentialKnown(:)
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
       REAL(KIND=dp) :: LX,LY,LZ,x,y,z
       LOGICAL :: Stat, CheckNormals
       COMPLEX(KIND=dp) :: R, dGdN, G, GradG(3), i
       REAL(KIND=dp) :: detJ,U,V,W,S,A,L,Normal(3),rad

       INTEGER :: j,k,p,q,t,dim
 
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------
       dim = CoordinateSystemDimension()
!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------
       i = CMPLX( 0.0d0, 1.0d0, KIND=dp )

       SELECT CASE( Element % Type % ElementCode / 100 )
       CASE(2)
          IntegStuff = GaussPoints( Element,4 )
       CASE(3)
          IntegStuff = GaussPoints( Element,6 )
       CASE(4)
          IntegStuff = GaussPoints( Element,16 )
       END SELECT

       CheckNormals = ASSOCIATED( Element % BoundaryInfo )
       IF ( CheckNormals ) THEN
          CheckNormals = ASSOCIATED( Element % BoundaryInfo % Left  ) .OR. &
                         ASSOCIATED( Element % BoundaryInfo % Right )
       END IF

       DO t=1,IntegStuff % n
          U = IntegStuff % u(t)
          V = IntegStuff % v(t)
          W = IntegStuff % w(t)
          S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!         Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
          stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
                 Basis, dBasisdx )
 
          S = S * detJ

          Normal = NormalVector( Element, Nodes, u,v, CheckNormals )

          LX = SUM( Nodes % x(1:n) * Basis )
          LY = SUM( Nodes % y(1:n) * Basis )
          LZ = SUM( Nodes % z(1:n) * Basis )

          DO p=1,BoundaryNodes
             k = BoundaryPerm(p)

             x = LX - Solver % Mesh % Nodes % x(k)
             y = LY - Solver % Mesh % Nodes % y(k)
             z = LZ - Solver % Mesh % Nodes % z(k)

             CALL Green( dim,AngularFrequency,x,y,z,G,GradG )
             dGdN = SUM( GradG * Normal )

             DO j=1,N
                q = PotentialPerm( Element % NodeIndexes(j) )

                IF ( PotentialKnown(q) ) THEN
                   IF ( p /= q ) THEN
                      FORCE(p) = FORCE(p) - s * Source(q) * Basis(j) * dGdN
                   END IF
                   STIFF(p,q) = STIFF(p,q) - s * Basis(j) * G
                ELSE
                   FORCE(p) = FORCE(p) + s * Source(q) * Basis(j) * G
                   IF ( p /= q ) THEN
                      STIFF(p,q) = STIFF(p,q) + s * Basis(j) * dGdN
                   END IF
                END IF

                R = -AngularFrequency * SIN( AngularFrequency*x ) * G * Normal(1)
                IF ( p /= q ) THEN
                   R = R - COS( AngularFrequency*x ) * dGdN
                END IF
                ADiagonal(p) = ADiagonal(p) - s * Basis(j)  * R
             END DO
          END DO
!------------------------------------------------------------------------------
       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE IntegrateMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
     SUBROUTINE Green( dim,k,x,y,z,W,GradW )
!------------------------------------------------------------------------------
       IMPLICIT NONE
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: x,y,z
       INTEGER :: dim
       REAL(KIND=dp) :: k
       COMPLEX(KIND=dp) :: W
       COMPLEX(KIND=dp), OPTIONAL :: GradW(:)
!------------------------------------------------------------------------------
       COMPLEX(KIND=dp) :: i,dWdR
       REAL(KIND=dp) :: r,J0,Y0,dJ0,dY0
!------------------------------------------------------------------------------
       R = SQRT( X**2 + Y**2 + Z**2 )
       i = CMPLX( 0.0d0,1.0d0,KIND=dp )

       SELECT CASE(dim)
       CASE(2)
          CALL Bessel( k*R, j0, y0, dj0, dy0 )
          W = (J0 - i*Y0) / (i*4)
          dWdR = k * (dJ0 - i*dY0) / (i*4)

       CASE(3)
          W = EXP(-i*k*R) / (4*PI*R)
          dWdR = (-i*k - 1/R) * W
       END SELECT

       IF ( PRESENT(GradW) ) THEN
          GradW(1) = x * dWdR / R
          GradW(2) = y * dWdR / R
          GradW(3) = z * dWdR / R
       END IF
!------------------------------------------------------------------------------
     END SUBROUTINE Green
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
     SUBROUTINE Bessel( x, j0, y0, dj0, dy0 )
!------------------------------------------------------------------------------
       IMPLICIT NONE

       INTEGER, PARAMETER          :: maxrounds = 1000
       DOUBLE PRECISION, PARAMETER :: polylimit = 10.0d0
       DOUBLE PRECISION, PARAMETER :: accuracy  = 1.0d-8
       DOUBLE PRECISION, PARAMETER :: gamma = 0.577215664901532860606512090082D0

       INTEGER :: i, k
       REAL(KIND=dp) :: hk
       REAL(kind=dp) :: x, j0, y0, dj0, dy0, phi, res, p, f
       
       DOUBLE PRECISION :: A(7) = &
          (/ 0.79788456D0, -0.00000077D0, -0.00552740D0,  &
             0.00009512D0,  0.00137237D0, -0.00072805D0,  &
             0.00014476D0 /)

       DOUBLE PRECISION :: B(7) = &
          (/ -0.78539816D0, -0.04166397D0, -0.00003954D0, &
              0.00262573D0, -0.00054125D0, -0.00029333D0, &
              0.00013558D0 /)

       DOUBLE PRECISION :: C(7) = &
          (/ 0.79788456D0,  0.00000156D0, 0.01659667D0,   &
             0.00017105D0, -0.00249511D0, 0.00113653D0,   &
            -0.0020033D0 /)

       DOUBLE PRECISION :: D(7) = &
          (/ -2.35619449D0, 0.12499612D0, 0.00005650D0,   &
             -0.00637879D0, 0.00074348D0, 0.00079824D0,   &
             -0.00029166D0 /)

       IF ( ABS(x) > polylimit ) THEN
!         Use polynomials
!         ---------------
          P = x
          F = 0.0d0
          DO k=1,7 
             F = F + A(k)*(3/x)**(k-1.0d0)
             P = P + B(k)*(3/x)**(k-1.0d0)
          END DO
          
          j0 = F * COS(P) / SQRT(x)
          y0 = F * SIN(P) / SQRT(x)
          
          P = x
          F = 0.0d0
          DO k=1,7 
             F = F + C(k)*(3/x)**(k-1.0d0)
             P = P + D(k)*(3/x)**(k-1.0d0)
          END DO
          
          dj0 = -F * COS(P) / SQRT(x)
          dy0 = -F * SIN(P) / SQRT(x)
       ELSE
!         Use series
!         ----------
          j0 = 1.0d0
          y0 = 0.0d0
          
          dj0 = 0.0d0 ! = - j1
          dy0 = 0.0d0 ! = - y1
          
          hk = 0.0d0
          
          DO k = 1,maxrounds
             hk = hk + 1.0d0 / k
             
             res = 1.0d0
             DO i = 1,k
                res = res * ( x / (2.0d0 * i) )**2
             END DO
             
             j0 = j0 + (-1)**k * res
             y0 = y0 + (-1)**(k+1) * hk * res
             
             dj0 = dj0 + (-1)**k * k / (0.5d0 * x) * res
             dy0 = dy0 + (-1)**(k+1) * hk * k / (0.5d0 * x) * res
             
             IF ( ABS(k / (0.5d0 * x) * res) < accuracy ) EXIT
          END DO
          
          IF ( k >= maxrounds ) STOP 'Error in evaluating Bessel functions'

          y0 = y0 + ( LOG(0.5d0 * x) + gamma ) * j0
          y0 = y0 * 2.0d0 / PI
          
          dy0 = dy0 + (1.0d0 / x) * j0 + ( LOG(0.5d0 * x) + gamma ) * dj0
          dy0 = dy0 * 2.0d0 / PI
       END IF
!------------------------------------------------------------------------------       
     END SUBROUTINE Bessel
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
     SUBROUTINE ComputePotential( Potential, Pot, Flx, Element, n, Nodes )
!------------------------------------------------------------------------------
       COMPLEX(KIND=dp) :: Pot(:), Flx(:), Potential(:)
       INTEGER :: n
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
       REAL(KIND=dp) :: LX,LY,LZ,x,y,z
       LOGICAL :: Stat, CheckNormals
       REAL(KIND=dp) :: detJ,U,V,W,S,A,L,Normal(3)

       COMPLEX(KIND=dp) :: dGdN, G, GradG(3)

       INTEGER :: i,j,k,p,q,t,dim
 
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------
       dim = CoordinateSystemDimension()
!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------
       SELECT CASE( Element % Type % ElementCode / 100 )
       CASE(2)
          IntegStuff = GaussPoints( Element,4 )
       CASE(3)
          IntegStuff = GaussPoints( Element,6 )
       CASE(4)
          IntegStuff = GaussPoints( Element,16 )
       END SELECT

       CheckNormals = ASSOCIATED( Element % BoundaryInfo )
       IF ( CheckNormals ) THEN
          CheckNormals = ASSOCIATED( Element % BoundaryInfo % Left  ) .OR. &
                         ASSOCIATED( Element % BoundaryInfo % Right )
       END IF

       DO t=1,IntegStuff % n
          U = IntegStuff % u(t)
          V = IntegStuff % v(t)
          W = IntegStuff % w(t)
          S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!         Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
          stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
                 Basis, dBasisdx )

          S = S * detJ
          Normal = NormalVector( Element, Nodes, u,v, CheckNormals )

          LX = SUM( Nodes % x(1:n) * Basis )
          LY = SUM( Nodes % y(1:n) * Basis )
          LZ = SUM( Nodes % z(1:n) * Basis )

          DO i=1,Solver % Mesh % NumberOfNodes
             IF ( PotentialPerm(i) > 0 ) CYCLE
             IF ( Solver % Variable % Perm(i) <= 0 ) CYCLE

             x = LX - Solver % Mesh % Nodes % x(i)
             y = LY - Solver % Mesh % Nodes % y(i)
             z = LZ - Solver % Mesh % Nodes % z(i)

             CALL Green( dim,AngularFrequency,x,y,z,G,GradG )
             dGdN = SUM( GradG * Normal )

             DO j=1,n
                q = PotentialPerm( Element % NodeIndexes(j) )
                Potential(i) = Potential(i) - s * Basis(j) * &
                      ( Pot(q) * dGdN - Flx(q) * G )
             END DO
          END DO
!------------------------------------------------------------------------------
       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE ComputePotential
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
     SUBROUTINE SolveFull( N,A,x,b,Solver )
!------------------------------------------------------------------------------
       TYPE(Solver_t) :: Solver
!------------------------------------------------------------------------------
       INTERFACE SolveLapack_cmplx
          SUBROUTINE SolveLapack_cmplx( N,A,x )
             USE Types
             INTEGER N
             COMPLEX(KIND=dp) :: a(n,n), x(n)
           END SUBROUTINE SolveLapack_cmplx
        END INTERFACE
!------------------------------------------------------------------------------
       INTEGER ::  N
 
       COMPLEX(KIND=dp) CONTIG ::  A(:,:),x(:),b(:)
!------------------------------------------------------------------------------

       SELECT CASE( ListGetString( Solver % Values, 'Linear System Solver' ) )

       CASE( 'direct' )
          CALL SolveLapack_cmplx( N, A, b )
          x(1:n) = b(1:n)

       CASE( 'iterative' )
          CALL FullIterSolver( N, x, b, Solver )

       CASE DEFAULT
          CALL Fatal( 'SolveFull', 'Unknown solver type.' )

       END SELECT
!------------------------------------------------------------------------------
     END SUBROUTINE SolveFull
!------------------------------------------------------------------------------


#include "huti_fdefs.h"
!------------------------------------------------------------------------------
     SUBROUTINE FullIterSolver( N,x,b,SolverParam )
!------------------------------------------------------------------------------
#ifdef USE_ISO_C_BINDINGS
    USE huti_sfe
#endif
       IMPLICIT NONE
!------------------------------------------------------------------------------
       TYPE(Solver_t) :: SolverParam
       INTEGER :: N
       COMPLEX(KIND=dp), DIMENSION(:) CONTIG :: x,b
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: dpar(50)

       INTEGER :: ipar(50),wsize
       ! REAL(KIND=dp), ALLOCATABLE :: Work(:,:)
       COMPLEX(KIND=dp), ALLOCATABLE :: Work(:,:)

       COMPLEX :: s

       LOGICAL :: AbortNotConverged

#ifndef USE_ISO_C_BINDINGS
       INTEGER  :: HUTI_Z_BICGSTAB
       EXTERNAL :: HUTI_Z_BICGSTAB
       INTEGER(KIND=addrInt) :: AddrFunc
#else
       INTEGER(KIND=AddrInt) :: AddrFunc
       EXTERNAL :: AddrFunc
#endif 
       INTEGER(KIND=addrInt) :: iterProc, mvProc, pcondProc, dProc
!------------------------------------------------------------------------------
       ipar = 0; dpar = 0
       dProc = 0

       HUTI_WRKDIM = HUTI_BICGSTAB_WORKSIZE
       wsize = HUTI_WRKDIM
       HUTI_NDIM = N
       ! ALLOCATE( Work(wsize,2*N) )
       ALLOCATE( Work(wsize,N) )

       IF ( ALL(x == 0.0) ) THEN
          HUTI_INITIALX = HUTI_RANDOMX
       ELSE
          HUTI_INITIALX = HUTI_USERSUPPLIEDX
       END IF

       HUTI_TOLERANCE = ListGetConstReal( Solver % Values, &
            'Linear System Convergence Tolerance' )

       HUTI_MAXTOLERANCE = ListGetConstReal( Solver % Values, &
            'Linear System Divergence Limit', GotIt )
       IF(.NOT. GotIt) HUTI_MAXTOLERANCE = 1.0d20       
       
       HUTI_MAXIT = ListGetInteger( Solver % Values, &
            'Linear System Max Iterations' )

       HUTI_DBUGLVL  = ListGetInteger( SolverParam % Values, &
            'Linear System Residual Output', GotIt )

       IF ( .NOT.Gotit ) HUTI_DBUGLVL = 1

       AbortNotConverged = ListGetLogical( SolverParam % Values, &
            'Linear System Abort Not Converged', GotIt )
       IF ( .NOT. GotIt ) AbortNotConverged = .TRUE.

!------------------------------------------------------------------------------
       iterProc  = AddrFunc(HUTI_Z_BICGSTAB)
       mvProc    = AddrFunc(Matvec)
       pcondProc = AddrFunc(precond)
       CALL IterCall( iterProc,x,b,ipar,dpar,work,mvProc,pcondProc, &
                 dProc, dProc, dProc, dProc )
!------------------------------------------------------------------------------

       DEALLOCATE( Work )

       IF ( HUTI_INFO /= HUTI_CONVERGENCE ) THEN
          IF ( AbortNotConverged ) THEN
             CALL Fatal( 'IterSolve', 'Failed convergence tolerances.' )
          ELSE
             CALL Error( 'IterSolve', 'Failed convergence tolerances.' )
          END IF
       END IF
!------------------------------------------------------------------------------
     END SUBROUTINE FullIterSolver 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   END SUBROUTINE HelmholtzBEMSolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE Precond( u,v,ipar )
!------------------------------------------------------------------------------
     USE GlobMatC
!------------------------------------------------------------------------------
     COMPLEX(KIND=dp) :: u(*),v(*)
     INTEGER :: ipar(*)
!------------------------------------------------------------------------------
     DO i=1,HUTI_NDIM
        u(i) = v(i) / Matrix(i,i)
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Precond
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE Matvec( u,v,ipar )
!------------------------------------------------------------------------------
      USE GlobMatC
!------------------------------------------------------------------------------
      COMPLEX(KIND=dp) :: u(*),v(*)
      INTEGER :: ipar(*)
!------------------------------------------------------------------------------
      v(1:HUTI_NDIM) = MATMUL( Matrix, u(1:HUTI_NDIM) )

!     DO i=1,HUTI_NDIM
!        v(i) = ( 0.0d0, 0.0d0 )
!        DO j=1,HUTI_NDIM
!           v(i) = v(i) + Matrix(i,j) * u(j)
!        END DO
!     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Matvec
!------------------------------------------------------------------------------
