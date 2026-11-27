!------------------------------------------------------------------------------
! A UMAT for the test UMAT_viscoelastic_adaptive: the standard linear solid
! (Zener model) for small strains, an elastic spring C_inf in parallel with a
! Maxwell branch of a spring C_1 and a dashpot of relaxation time tau. Both
! springs are isotropic with the same Poisson ratio.
!
! The state variables are the stress of the Maxwell branch at the previous
! time level. Its rate equation d(sigma_v)/dt = C_1 d(eps)/dt - sigma_v/tau
! is integrated by the implicit Euler method over the time increment DTIME.
!
! PROPS: density, E_inf, Poisson ratio, E_1, tau, and the largest DTIME/tau
! (0 for no limit). With the limit the UMAT asks for the timestep through
! PNEWDT: below one the step is repeated with a smaller one.
!------------------------------------------------------------------------------
  SUBROUTINE zener(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
       rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
       predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
       coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
       kstep, kinc)
!------------------------------------------------------------------------------
    USE Types
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDI, NSHR, NTENS, NSTATEV, NPROPS, NOEL, NPT, &
        layer, kspt, kstep, kinc
    REAL(KIND=dp), INTENT(INOUT) :: STRESS(NTENS), STATEV(NSTATEV), SSE, SPD, SCD, pnewdt
    REAL(KIND=dp), INTENT(OUT) :: DDSDDE(NTENS,NTENS), rpl, ddsddt(NTENS), &
        drplde(NTENS), drpldt
    REAL(KIND=dp), INTENT(IN) :: STRAN(NTENS), DSTRAN(NTENS), TIME(2), DTIME, &
        TEMP, dTemp, predef(1), dpred(1), PROPS(NPROPS), coords(3), drot(3,3), &
        celent, DFRGRD0(3,3), DFRGRD1(3,3)
    CHARACTER(len=80), INTENT(IN) :: CMNAME
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Cinf(NTENS,NTENS), C1(NTENS,NTENS), f, tau, MaxRatio

    CALL Elastic( PROPS(2), PROPS(3), Cinf )
    CALL Elastic( PROPS(4), PROPS(3), C1 )
    tau = PROPS(5)
    MaxRatio = PROPS(6)

    f = 1.0_dp / ( 1.0_dp + DTIME / tau )
    STATEV(1:NTENS) = f * ( STATEV(1:NTENS) + MATMUL( C1, DSTRAN ) )
    STRESS = MATMUL( Cinf, STRAN + DSTRAN ) + STATEV(1:NTENS)
    DDSDDE = Cinf + f * C1

    IF( MaxRatio > 0.0_dp ) pnewdt = MaxRatio * tau / DTIME

    rpl = 0.0_dp; ddsddt = 0.0_dp; drplde = 0.0_dp; drpldt = 0.0_dp

  CONTAINS

    SUBROUTINE Elastic( E, nu, C )
      REAL(KIND=dp) :: E, nu, C(:,:)
      REAL(KIND=dp) :: Lambda, Mu
      INTEGER :: i
      Lambda = E * nu / ( (1.0_dp+nu) * (1.0_dp-2.0_dp*nu) )
      Mu = E / (2.0_dp * (1.0_dp + nu))
      C = 0.0_dp
      C(1:NDI,1:NDI) = Lambda
      DO i=1,NTENS
        C(i,i) = C(i,i) + Mu
      END DO
      DO i=1,NDI
        C(i,i) = C(i,i) + Mu
      END DO
    END SUBROUTINE Elastic
!------------------------------------------------------------------------------
  END SUBROUTINE zener
!------------------------------------------------------------------------------
