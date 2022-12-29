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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 29.9.2022
! *
! *****************************************************************************/

!-----------------------------------------------------------------------------
!>  Solve the Cahn-Hillard interface equations as a strongly coupled system.
!-----------------------------------------------------------------------------
SUBROUTINE CahnHillardSolver( Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  REAL(KIND=dp) :: Norm
  INTEGER :: iter,  n, nb, nd, t, istat, active

  LOGICAL :: Found, NewtonLinearization=.FALSE.

  TYPE(ValueList_t), POINTER :: Params

  INTEGER :: NewtonIter,  NonlinIter
  REAL(KIND=dp):: NewtonTol, NonlinTol, RelativeChange
!------------------------------------------------------------------------------
  Params => GetSolverParams()

  NewtonIter = GetInteger(Params,'Nonlinear System Newton After Iterations',Found)
  NewtonTol  = GetConstReal(Params,'Nonlinear System Newton After Tolerance',Found)

  NonlinIter = GetInteger(Params,'Nonlinear System Max Iterations',Found)
  IF(.NOT.Found) NonlinIter=1

  NonlinTol = GetConstReal(Params,'Nonlinear System Convergence Tolerance',Found)

  DO iter=1,NonlinIter
    CALL DefaultInitialize()
    CALL BulkAssembly()
    CALL DefaultFinishBulkAssembly()
    CALL BoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    Norm = DefaultSolve()

    RelativeChange = Solver % Variable % NonlinChange

    WRITE( Message, * ) 'Result Norm   : ',Norm
    CALL Info( 'Cahn-Hillard interface Solver', Message, Level=4 )
    WRITE( Message, * ) 'Relative Change : ',RelativeChange
    CALL Info( 'Cahn-Hillard interface Solver', Message, Level=4 )

    IF ( RelativeChange < NewtonTol .OR. iter >= NewtonIter ) &
             NewtonLinearization = .TRUE.
    IF ( RelativeChange < NonlinTol ) EXIT
  END DO

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
!------------------------------------------------------------------------------
    INTEGER :: t,n,nd

!$omp parallel do private(Element,n,nd)
    DO t=1,GetNOFActive()
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      CALL LocalMatrix( Element, n, nd )
     END DO
!$omp end parallel do
!------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: MASS(2*nd,2*nd), STIFF(2*nd,2*nd), FORCE(2*nd), LOAD(2,nd)
    REAL(KIND=dp), POINTER :: A(:,:),M(:,:)

    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3), DetJ, SOL(2,nd), Phi, s
    LOGICAL :: Stat,isScalar,Found, Found_dens, Found_scoeff, Found_mobil
    INTEGER :: i,j,p,q,t,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: Velo(3), h, MixEnergyDensity, Mobility
    REAL(KIND=dp) :: Velo_n(3,n), h_n(n), MixEnergyDensity_n(n), Mobility_n(n), st_n(n)

    TYPE(ValueList_t), POINTER :: Material, BF

    TYPE(Nodes_t), SAVE :: Nodes
!$omp threadprivate(Nodes)
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes,Element )

    Material => GetMaterial(Element)

    dim = CoordinateSystemDimension()

    Velo_n = 0._dp
    DO i=1,dim
      Velo_n(i,:) = GetReal( Material, 'Convection Velocity '//i2s(i), Found )
    END DO

    h_n = GetReal( Material, 'Interface Thickness',Found)
    IF(.NOT.Found) THEN
      h_n = ElementDiameter(Element,Nodes)/2
    END IF

    Mobility_n = GetReal( Material, 'Mobility',Found_mobil,Element )

    mixEnergyDensity_n = GetReal( Material, 'Mixing Energy Density',Found_dens )
    IF(.NOT. Found_dens ) THEN
      st_n = Getreal( Material, 'Surface Tension Coefficient', Found_scoeff)
    END IF

    CALL GetVectorLocalSolution(SOL,UElement=Element)

    MASS  = 0.0_dp
    STIFF = 0.0_dp
    FORCE = 0.0_dp

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
       IP % W(t), detJ, Basis, dBasisdx )

      Velo = 0._dp
      DO i=1,dim
        Velo(i) = SUM(Velo_n(i,:)*Basis(1:n))
      END DO

      h = SUM(h_n*Basis(1:n))

      IF (Found_mobil) THEN
        Mobility = SUM(mobility_n*Basis(1:n))
      ELSE
        Mobility = h**2
      END IF

      MixEnergyDensity = 1._dp
      IF(Found_dens) THEN
        MixEnergyDensity = SUM(MixEnergyDensity_n*Basis(1:n))
      ELSE IF (Found_scoeff) THEN
        MixEnergyDensity = 3*h/2/SQRT(2._dp)*SUM(st_n*Basis(1:n))
      END IF

      Phi = SUM( SOL(1,1:nd)*Basis(1:nd) )

      s = IP % s(t) * DetJ

      DO p=1,nd
        i=2*(p-1)
        DO q=1,nd
          j=2*(q-1)
          M => MASS(i+1:i+2,j+1:j+2)
          A => STIFF(i+1:i+2,j+1:j+2)

          ! Phi-eq: DPhi/Dt - div(m*l/h**2*grad(theta)) = 0
 
          ! Time derivative
          ! ---------------
          M(1,1) = M(1,1) + s*Basis(q)*Basis(p)

          ! Convection
          ! ---------------
          A(1,1) = A(1,1) + s*SUM(Velo*dBasisdx(q,:))*Basis(p)


          ! diffusion
          ! -----------------
          A(1,2) = A(1,2) + s*Mobility*MixEnergyDensity/h**2* &
                  SUM(dBasisdx(q,:)*dBasisdx(p,:))

          ! Theta-eq: Theta + div(h**2*grad(Phi)) + Phi - Phi**3 = 0

          A(2,2) = A(2,2) + s*Basis(q)*Basis(p)
          A(2,1) = A(2,1) - s*h**2*SUM(dBasisdx(q,:)*dBasisdx(p,:))
          A(2,1) = A(2,1) + s*(1-3*Phi**2)*Basis(q)*Basis(p)
        END DO

        FORCE(i+2) = FORCE(i+2) - s * 2*Phi**3 * Basis(p)
      END DO
    END DO

    IF(Transient) THEN
      CALL Default1stOrderTime(MASS,STIFF,FORCE,UElement=Element)
    END IF
    CALL DefaultUpdateEquations( STIFF, FORCE,UElement=Element )
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE BoundaryAssembly()
!------------------------------------------------------------------------------
    INTEGER :: t,n,nd

!$omp parallel do private(Element,n,nd)
    DO t=1,GetNOFBoundaryElements()
      Element => GetBoundaryElement(t)
      IF(.NOT.ActiveBoundaryElement(Element)) CYCLE

      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      CALL BoundaryLocalMatrix( Element, n, nd )
     END DO
!$omp end parallel do
!------------------------------------------------------------------------------
  END SUBROUTINE BoundaryAssembly
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE BoundaryLocalMatrix( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: MASS(2*nd,2*nd), STIFF(2*nd,2*nd), FORCE(2*nd)
    REAL(KIND=dp), POINTER :: A(:,:),M(:,:)

    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3), s, DetJ
    LOGICAL :: Stat,Found
    INTEGER :: i,j,p,q,t
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: HF_n(n), ET_n(n), HT_n(n), EF_n(n), HF, ET, HT, EF

    TYPE(ValueList_t), POINTER :: BC

    TYPE(Nodes_t), SAVE :: Nodes
!$omp threadprivate(Nodes)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  END SUBROUTINE BoundaryLocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE CahnHillardSolver
!------------------------------------------------------------------------------
