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
! *  Original Date: 22.12.2011
! *
! *****************************************************************************/

!-----------------------------------------------------------------------------
!>  Solve the Thermoelectric equations as a strongly coupled system.
!-----------------------------------------------------------------------------
SUBROUTINE ThermoElectricSolver( Model,Solver,dt,Transient)
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
    CALL Info( 'ThermoElectricSolver', Message, Level=4 )
    WRITE( Message, * ) 'Relative Change : ',RelativeChange
    CALL Info( 'ThermoElectricSolver', Message, Level=4 )

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

    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3), DetJ, SOL(2,nd)
    LOGICAL :: Stat,isScalar,Found
    INTEGER :: i,j,p,q,t
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: EF(3),TG(3),CD(3),JH,HS,IC,SB(3,3),s, PeltSB(3,3), PeltSi(3,3), &
       sigma(3,3),alpha(3,3), pelt(3,3), hcond(3,3), epsil(3,3),rho,c_p,Temp

    REAL(KIND=dp) :: sigma_n(3,3,n),alpha_n(3,3,n), pelt_n(3,3,n), &
            hcond_n(3,3,n), epsil_n(3,3,n),rho_n(n),c_p_n(n),epsil0

    TYPE(ValueList_t), POINTER :: Material, BF

    TYPE(Nodes_t), SAVE :: Nodes
!$omp threadprivate(Nodes)
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes,Element )

    Material => GetMaterial(Element)

    CALL InputTensor(Material, hcond_n, isScalar,'Heat Conductivity', Element)
    CALL InputTensor(Material, alpha_n, isScalar,'Seebeck Coefficient', Element)
    CALL InputTensor(Material, sigma_n, isScalar,'Electric Conductivity', Element)

    IF( ListCheckPresent( Material,'Relative Permittivity') ) THEN
      epsil0 = GetConstReal( Model % Simulation,'Permittivity of Vacuum',Found)
      IF(.NOT. Found ) THEN
        CALL Fatal('ThermoElectricSolver',&
            '> Relative Permittivity < requires > Permittivity of Vacuum < !')
      END IF
      CALL InputTensor(Material, epsil_n, isScalar,'Electric Permittivity', Element)
    ELSE
      epsil0 = 1.0_dp
      CALL InputTensor(Material, epsil_n, isScalar,'Electric Permittivity', Element)
    END IF


    rho_n = GetReal( Material, 'Density',Found,Element )
    c_p_n = GetReal( Material, 'Heat Capacity',Found,Element )

    BF => GetBodyForce(Element)
    IF (ASSOCIATED(BF)) THEN      
      LOAD(1,1:n) = GetReal(BF,'Heat Source',Found,Element)
      IF( Found ) THEN
        Load(1,1:n) = rho_n * Load(1,1:n) 
      ELSE
        LOAD(1,1:n) = GetReal(BF,'Volumetric Heat Source',Found,Element)
      END IF
      LOAD(2,1:n) = GetReal(BF,'Current Source',Found,Element)
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

      DO i=1,3
        DO j=1,3
          hcond(i,j) = SUM(hcond_n(i,j,:)*Basis(1:n))
          alpha(i,j) = SUM(alpha_n(i,j,:)*Basis(1:n))
          sigma(i,j) = SUM(sigma_n(i,j,:)*Basis(1:n))
          epsil(i,j) = epsil0 * SUM(epsil_n(i,j,:)*Basis(1:n))
        END DO
      END DO
      rho = SUM(rho_n*Basis(1:n))
      c_p = SUM(c_p_n*Basis(1:n))

      IC = SUM( LOAD(2,1:n)*Basis(1:n) )    ! input current
      HS = SUM( LOAD(1,1:n)*Basis(1:n) )    ! input heat source

      EF = -MATMUL(SOL(2,1:nd), dBasisdx(1:nd,:)) ! electric field from prev.sol.
      TG =  MATMUL(SOL(1,1:nd), dBasisdx(1:nd,:)) ! temperature gradient  ""

      CD = MATMUL(sigma,EF-MATMUL(alpha,TG)) ! current density from prev.sol
      JH = SUM(CD*EF)                        ! joule heating       ""

      SB = MATMUL(sigma,alpha)

      Temp = SUM(Basis(1:nd)*SOL(1,1:nd)) ! temperature from previous iteration
      Pelt = Temp*alpha                   ! Peltier coefficient
      PeltSB = MATMUL(Pelt,SB)
      PeltSi = MATMUL(Pelt,Sigma)

      s = IP % s(t) * DetJ

      DO p=1,nd
        i=2*(p-1)
        DO q=1,nd
          j=2*(q-1)
          M => MASS(i+1:i+2,j+1:j+2)
          A => STIFF(i+1:i+2,j+1:j+2)

          ! thermal damping
          ! ---------------
          M(1,1) = M(1,1) + s*rho*c_p*Basis(q)*Basis(p)

          ! thermal diffusion
          ! -----------------
          A(1,1) = A(1,1) + s*SUM(MATMUL(hcond,dBasisdx(q,:))*dBasisdx(p,:))

          IF (NewtonLinearization) THEN

            ! Newton linarization of: div(pelt*J) (= div(T*alpha*J))
            ! ======================================================

            ! Peltier coefficient implicitly
            ! -------------------------------
            A(1,1) = A(1,1) - &
                     s*Basis(q)*SUM(MATMUL(alpha,CD)*dBasisdx(p,:))

            ! temperature gradient part of div(pelt_0*J)
            ! ------------------------------------------
            A(1,1) = A(1,1) - &
                  s*SUM(MATMUL(PeltSB,dBasisdx(q,:))*dBasisdx(p,:))

            ! electric field part of div(pelt_0*J)
            ! ------------------------------------
            A(1,2) = A(1,2) - &
                  s*SUM(MATMUL(PeltSi,dBasisdx(q,:))*dBasisdx(p,:))

            ! Newton linarization of Joule heating=-(J,E) ~
            ! -(J_0,E) + (J,E_0) - (J_0,E_0) (<-- to rhs)
            ! =============================================
            ! -(J_0,E)
            ! --------
            A(1,2) = A(1,2) + s*SUM(CD*dBasisdx(q,:))*Basis(p)

            ! temperature gradient part of (J,E_0)
            ! ------------------------------------
            A(1,1) = A(1,1) - &
                   s*SUM(MATMUL(SB,dBasisdx(q,:))*EF)*Basis(p)

            ! electric field part part of (J,E_0)
            ! ------------------------------------
            A(1,2) = A(1,2) - &
                s*SUM(MATMUL(sigma,dBasisdx(q,:))*EF)*Basis(p)
          END IF


          ! displacement current
          ! --------------------
          M(2,2) = M(2,2) + &
             s*SUM(MATMUL(epsil,dBasisdx(q,:))*dBasisdx(p,:))

          ! Seebeck
          ! -------
          A(2,1) = A(2,1) + &
                s*SUM(MATMUL(SB,dBasisdx(q,:))*dBasisdx(p,:))

          ! electric diffusion
          ! ------------------
          A(2,2) = A(2,2) + &
             s*SUM(MATMUL(sigma,dBasisdx(q,:))*dBasisdx(p,:))
        END DO

        ! heat load
        ! ---------
        FORCE(i+1) = FORCE(i+1) + s * (HS+JH)*Basis(p)

        ! div(pelt_0*J_0)
        ! ---------------
        FORCE(i+1) = FORCE(i+1) + &
                  s*SUM(MATMUL(Pelt,CD)*dBasisdx(p,:))


        ! electric load
        ! -------------
        FORCE(i+2) = FORCE(i+2) + s * IC * Basis(p)
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

    BC => GetBC(Element)
    IF(.NOT.ASSOCIATED(BC)) RETURN

    CALL GetElementNodes( Nodes,Element )

    HF_n = GetReal(BC, 'Heat Flux', Found, Element)
    ET_n = GetReal(BC, 'External Temperature', Found, Element)
    HT_n = GetReal(BC, 'Heat Transfer Coefficient', Found, Element)

    EF_n = GetReal(BC, 'Electric Flux', Found, Element)

    MASS  = 0.0_dp
    STIFF = 0.0_dp
    FORCE = 0.0_dp

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
       IP % W(t), detJ, Basis, dBasisdx )

      s = IP % s(t) * DetJ

      HF = SUM(HF_n*Basis(1:n)) ! heat flux
      ET = SUM(ET_n*Basis(1:n)) ! external temperature
      HT = SUM(HT_n*Basis(1:n)) ! heat transfer coeffiecient

      EF = SUM(EF_n*Basis(1:n)) ! electric flux

      DO p=1,nd
        i=2*(p-1)
        DO q=1,nd
          j=2*(q-1)
          A => STIFF(i+1:i+2,j+1:j+2)
          A(1,1) = A(1,1) + s*HT*Basis(q)*Basis(p)
        END DO
        FORCE(i+1) = FORCE(i+1) + s*(HF+HT*ET)*Basis(p)
        FORCE(i+2) = FORCE(i+2) + s*EF*Basis(p)
      END DO
    END DO

    IF(Transient) THEN
      CALL Default1stOrderTime(MASS,STIFF,FORCE,UElement=Element)
    END IF
    CALL DefaultUpdateEquations( STIFF, FORCE,UElement=Element )
!------------------------------------------------------------------------------
  END SUBROUTINE BoundaryLocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE InputTensor( Material, Tensor, IsScalar, Name, Element )
!------------------------------------------------------------------------------
      TYPE(ValueList_t), POINTER :: Material
      REAL(KIND=dp) :: Tensor(:,:,:)
      LOGICAL :: IsScalar
      CHARACTER(LEN=*) :: Name
      TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
      INTEGER :: i,j,n
      LOGICAL ::  stat

      REAL(KIND=dp), POINTER :: Hwrk(:,:,:) => NULL()
!$omp threadprivate(Hwrk)
!------------------------------------------------------------------------------
      Tensor = 0.0_dp
      IsScalar = .TRUE.

      n = Element % TYPE % NumberOfNodes
      CALL ListGetRealArray( Material, Name, Hwrk, n, &
              Element % NodeIndexes, stat )
      IF ( .NOT. stat ) RETURN

      IsScalar = SIZE(HWrk,1) == 1 .AND. SIZE(HWrk,2) == 1
      IF ( IsScalar ) THEN
        DO i=1,SIZE(Tensor,1)
          Tensor(i,i,1:n) = Hwrk(1,1,1:n)
        END DO
      ELSE
        IF ( SIZE(Hwrk,1) == 1 ) THEN
           DO i=1,MIN(3,SIZE(HWrk,2) )
              Tensor( i,i,1:n ) = Hwrk( 1,1,1:n )
           END DO
        ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN
           DO i=1,MIN(3,SIZE(Hwrk,1))
              Tensor( i,i,1:n ) = Hwrk( i,1,1:n )
           END DO
        ELSE
          DO i=1,MIN(3,SIZE(Hwrk,1))
             DO j=1,MIN(3,SIZE(Hwrk,2))
                Tensor( i,j,1:n ) = Hwrk( i,j,1:n )
             END DO
          END DO
        END IF
      END IF
!------------------------------------------------------------------------------
   END SUBROUTINE InputTensor
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE ThermoElectricSolver
!------------------------------------------------------------------------------
