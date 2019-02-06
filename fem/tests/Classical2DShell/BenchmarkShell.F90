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
! *  A test case for solving the two-dimensional Reissner-Naghdi shell equations 
! *  over a straight cylindrical shell parametrized in (exact) lines of curvature 
! *  coordinates. 
! *
! *  Authors: Mika Malinen
! *  Email:   mika.malinen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: Feb 6, 2019
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
SUBROUTINE ShellSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE ElementDescription

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables:
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: SolverPars
  TYPE(Element_t), POINTER :: Element

  INTEGER :: ShellModelPar, csize
  INTEGER :: n, nb, nd, k, active 
  INTEGER :: ElementSought

  REAL(KIND=dp) :: t, RefWork, Work, Norm
  REAL(KIND=dp) :: Energy(5), ResVec(20)
  REAL(KIND=dp) :: LocalSol(6,9)

  REAL(KIND=dp) :: e1(3), e2(3), e3(3)
!------------------------------------------------------------------------------  
  SolverPars => GetSolverParams()
  Mesh => GetMesh()

  ShellModelPar = ListGetInteger(SolverPars, 'Variable DOFs', minv=5, maxv=6)
  IF (ShellModelPar == 6) THEN
    csize = 4
  ELSE
    csize = 3
  END IF

  ! ------------------------------------------------------------------------------
  ! Assembly loop for generating the shell stiffness matrix:
  ! ------------------------------------------------------------------------------
  CALL DefaultInitialize()
  Active = GetNOFActive()
  DO k=1,Active
    Element => GetActiveElement(k)

    n  = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()

    ! ------------------------------------------------------------------------------
    ! Generate the element stiffness matrix and assemble the local contribution:
    ! -----------------------------------------------------------------------------
    CALL ShellLocalMatrix(Element, n, nd+nb, ShellModelPar, csize)

  END DO

  CALL DefaultFinishBulkAssembly() 
  CALL DefaultFinishAssembly()

  CALL DefaultDirichletBCs()

  ! And finally, solve:
  !--------------------
  Norm = DefaultSolve()

  Work = 8.0d0*SUM(Solver % Variable % Values(:) * Solver % Matrix % RHS(:))
  t = 1.0d-1
!  RefWork = 12.0d0*(1.0d0-(1.0d0/3.0d0)**2)*(1.0d9)**2/7.0d10 * 2.688287959059254*t ! t=0.01
  RefWork = 12.0d0*(1.0d0-(1.0d0/3.0d0)**2)*(1.0d9)**2/7.0d10 * 1.828629366566552*t ! t=0.1
  PRINT *, 'Relative energy error = ', SQRT(ABS(RefWork-Work)/RefWork)

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE ShellLocalMatrix(Element, n, nd, m, csize)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd, m, csize 
    !------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BodyForce
    TYPE(Nodes_t) :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP

    LOGICAL :: Stat, Found
    INTEGER :: i, k, p, t

    REAL(KIND=dp) :: MASS(m*nd,m*nd), STIFF(m*nd,m*nd), FORCE(m*nd), LOAD(n)
    REAL(KIND=dp) :: BM(csize,m*nd), BS(2,m*nd), BB(3,m*nd)
    REAL(KIND=dp) :: PoissonRatio(n), YoungsMod(n), ShellThickness(n)
    REAL(KIND=dp) :: CMat(4,4), GMat(2,2)
    REAL(KIND=dp) :: nu, E, h, NormalTraction
    REAL(KIND=dp) :: R
    REAL(KIND=dp) :: Basis(nd), dBasisdx(nd,3), DetJ, LoadAtIP, Weight

    REAL(KIND=dp) :: y1, y2

    SAVE Nodes
    !------------------------------------------------------------------------------
    R = 1.0d0   ! R = principal radius of curvature, y1 = angular dir, y2 = axial dir  

    CALL GetElementNodes(Nodes)

    ! --------------------------------------------------------------------------
    ! Body forces, material parameters and the shell thickness:
    ! --------------------------------------------------------------------------
    PoissonRatio(1:n) = GetReal(GetMaterial(),'Poisson Ratio')
    YoungsMod(1:n) = GetReal(GetMaterial(),'Youngs Modulus')
    ShellThickness(1:n) = GetReal(GetMaterial(),'Shell Thickness')
    LOAD = 0.0_dp
    BodyForce => GetBodyForce()
    IF ( ASSOCIATED(BodyForce) ) &
        Load(1:n) = GetReal( BodyForce, 'Normal Pressure', Found )

    MASS = 0.0d0
    STIFF = 0.0d0
    FORCE = 0.0d0
    BM = 0.0d0
    BS = 0.0d0
    BB = 0.0d0

    IP = GaussPoints( Element )

    DO t=1,IP % n

      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )

      ! ------------------------------------------------
      ! Data interpolation:
      ! ------------------------------------------------
      nu = SUM( PoissonRatio(1:n) * Basis(1:n) )
      E = SUM( YoungsMod(1:n) * Basis(1:n) )
      h = SUM( ShellThickness(1:n) * Basis(1:n) )
      NormalTraction = SUM( Load(1:n) * Basis(1:n) )
      !
      ! Use a hard-coded load to avoid errors from representing the load: 
      !
      y1 = SUM( Nodes % x(1:n) * Basis(1:n) )
      NormalTraction = h * 1.0d9 * cos(2.0d0*y1)

      CALL ElasticityMatrix(CMat, GMat, 1.0D0, 1.0D0, E, nu)

      !---------------------------------------------------------------
      ! The part corresponding to the membrane strains:
      !---------------------------------------------------------------
      Weight = h * detJ * IP % s(t)
      DO p=1,nd
        BM(1,(p-1)*m+1) = dBasisdx(p,1) 
        BM(1,(p-1)*m+3) = Basis(p)/R 

        BM(2,(p-1)*m+2) = dBasisdx(p,2)

        BM(3,(p-1)*m+1) = dBasisdx(p,2)
        BM(3,(p-1)*m+2) = dBasisdx(p,1)

        IF (m==6) THEN
        !IF (.FALSE.) THEN
          !----------------------------------------------------------------
          ! Determine the thickness-stretch parameter via augmented energy:
          !----------------------------------------------------------------
          BM(4,(p-1)*m+6) = Basis(p)
          BM(4,(p-1)*m+1) = -nu/(1.0d0-nu) * dBasisdx(p,1)
          BM(4,(p-1)*m+2) = -nu/(1.0d0-nu) * dBasisdx(p,2)
          BM(4,(p-1)*m+3) = -nu/(1.0d0-nu) * Basis(p)/R
        END IF

      END DO
      CALL StrainEnergyDensity(Stiff, CMat(1:csize,1:csize), BM, csize, m*nd, Weight)

      !---------------------------------------------------------------
      ! The part corresponding to the shear strains:
      !---------------------------------------------------------------      
      DO p=1,nd       
        BS(1,(p-1)*m+1) = -Basis(p)/R        
        BS(1,(p-1)*m+3) = dBasisdx(p,1)
        BS(1,(p-1)*m+4) = -Basis(p) 

        BS(2,(p-1)*m+3) = dBasisdx(p,2)
        BS(2,(p-1)*m+5) = -Basis(p) 
      END DO
      CALL StrainEnergyDensity(Stiff, GMat, BS, 2, m*nd, Weight)

      !---------------------------------------------------------------
      ! The part corresponding to the bending strains (terms depending
      ! on the membrane strains dropped at the moment):
      !---------------------------------------------------------------      
      Weight = h**3/12.0d0 * detJ * IP % s(t)
      DO p=1,nd
        BB(1,(p-1)*m+4) = dBasisdx(p,1) 

        BB(2,(p-1)*m+5) = dBasisdx(p,2)

        BB(3,(p-1)*m+4) = dBasisdx(p,2)
        BB(3,(p-1)*m+5) = dBasisdx(p,1)
        BB(3,(p-1)*m+1) = -1.0d0/R * dBasisdx(p,2)
      END DO
      CALL StrainEnergyDensity(Stiff, CMat(1:3,1:3), BB, 3, m*nd, Weight)

      !IF (m==6) THEN
      IF (.FALSE.) THEN
        !----------------------------------------------------------------
        ! Determine the thickness-stretch parameter corresponding to the
        ! state of vanishing normal stress:
        !----------------------------------------------------------------
        Weight = detJ * IP % s(t)
        DO p=1,nd
          DO k=1,nd
            Stiff(p*6,k*6) = Basis(p) * Basis(k) * Weight
            Stiff(p*6,(k-1)*m+1) = -nu/(1.0d0-nu) * Basis(p) * dBasisdx(k,1) * Weight
            Stiff(p*6,(k-1)*m+3) = -nu/(1.0d0-nu) * Basis(p) * Basis(k)/R * Weight
            Stiff(p*6,(k-1)*m+2) = -nu/(1.0d0-nu) * Basis(p) * dBasisdx(k,2)  * Weight
          END DO
        END DO
      END IF

      !----------------------------------------------------------------
      ! RHS vector:
      !--------------------------------------------------------
      Weight = detJ * IP % s(t)
      DO p=1,nd
        i = m*(p-1)+3
        Force(i) = Force(i) + NormalTraction * Basis(p) * Weight
      END DO

    END DO

    CALL DefaultUpdateEquations(STIFF,FORCE)  

!------------------------------------------------------------------------------
  END SUBROUTINE ShellLocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE ElasticityMatrix(CMat, GMat, A1, A2, E, nu)
!------------------------------------------------------------------------------
! The matrix representation of the elasticity tensor with respect an orthogonal
! basis. The case A1 = A2 = 1 corresponds to an orthonormal basis.     
!------------------------------------------------------------------------------    
    IMPLICIT NONE

    REAL(KIND=dp) :: CMat(4,4), GMat(2,2), A1, A2, E, nu  
!------------------------------------------------------------------------------
    CMat = 0.0d0
    GMat = 0.0d0

    CMat = 0.0d0
    CMat(1,1) = 1.0d0    
    CMat(1,2) = nu
    CMat(2,1) = nu
    CMat(2,2) = 1.0d0
    CMat(3,3) = (1.0d0-nu)/2.0d0
    CMat = CMat * E / (1.0d0-nu**2)

    CMat(1,1) = CMat(1,1)/A1**4
    CMat(1,2) = CMat(1,2)/(A1**2 * A2**2)
    CMat(2,1) = CMat(2,1)/(A1**2 * A2**2)
    CMat(2,2) = CMat(2,2)/A2**4   
    CMat(3,3) = CMat(3,3)/(A1**2 * A2**2)

    ! The row corresponding to the normal stress:
    CMat(4,4) = (1.0d0-nu) * E /( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
    !CMat(4,1) = nu * E /(A1**2*(1.0d0+nu)*(1.0d0-2.0d0*nu))
    !CMat(4,2) = CMat(4,1) 
    !CMat(4,4) = (1.0d0-nu) * E /(A1**2*(1.0d0+nu)*(1.0d0-2.0d0*nu))

    GMat(1,1) = E/(2.0d0*(1.0d0 + nu)*A1**2)
    GMat(2,2) = E/(2.0d0*(1.0d0 + nu)*A2**2)
!------------------------------------------------------------------------------
  END SUBROUTINE ElasticityMatrix
!------------------------------------------------------------------------------    

!------------------------------------------------------------------------------
  SUBROUTINE StrainEnergyDensity(A, B, C, m, n, s)
!------------------------------------------------------------------------------
! Performs the operation
!
!    A = A + C' * B * C * s
!
! with
!
!    Size( A ) = n x n
!    Size( B ) = m x m
!    Size( C ) = m x n
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: A(n,n), B(m,m), C(m,n), s
    INTEGER :: m, n
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l
!------------------------------------------------------------------------------
    DO i=1,n
       DO j=1,n
          DO k=1,m
             DO l=1,m
                A(i,j) = A(i,j) + C(k,i)*B(k,l)*C(l,j) * s
             END DO
          END DO
       END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE StrainEnergyDensity
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE ShellSolver
!------------------------------------------------------------------------------
