!
! *  Authors: Mika Malinen
! *  Email:   mika.malinen@csc.fi
! *  Original Date: March 8, 2019
!  
! Adds plastic behviour models bi-linear, hsk=bi 
!bi_linear 
!
!small strain 6 material constants nu strain1, stress1, strain2, stress2, ultimate strain
!seven state vatiables max strain, mid strain, min strain,
! max principal, mid principal, min principle, mises
!
!bi-linear stress strain curve, ultimate stress at final strain value
!
!hsk_bi
!
!large strain 6 material constants nu strain1, stress1, strain2, stress2, ultimate strain
!seven state vatiables max strain, mid strain, min strain,
! max principal, mid principal, min principle, mises
!
!bi-linear stress strain curve, ultimate stress at final strain value
!
! plastic point wise stress strain curve  
! nu, strain, stress pairs up to nprops  
!  
!------------------------------------------------------------------------------
! The template for including a material model definition written in the form of
! an Abaqus user subroutine (UMAT). The arguments which can be supposed to be
! supported by Elmer are capitalized.
!------------------------------------------------------------------------------
  SUBROUTINE UMAT_template(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
       rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
       predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
       coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
       kstep, kinc)
!------------------------------------------------------------------------------
    USE Types
    USE Messages
    IMPLICIT NONE

    REAL(KIND=dp), INTENT(INOUT) :: STRESS(NTENS)
    ! Requirement for Elmer: At the time of calling the Cauchy stress T_n before
    ! the time/load increment is given
    ! Requirement for umat:  The stress T_{n+1}^{(k)} corresponding to the
    ! current approximation of the strain increment (DSTRAN) must be returned.
    ! If the strain increment is defined to be zero in the beginning of the
    ! nonlinear iteration, Elmer will generate a candidate for the strain increment
    ! by assuming purely elastic increment characterized by DDSDDE.

    REAL(KIND=dp), INTENT(INOUT) :: STATEV(NSTATEV)
    ! Requirement for Elmer: The state variables Q_n as specified at the
    ! previous time/load level for converged solution are given.
    ! Requirement for umat:  The state variables Q_{n+1}^{(k)} corresponding to
    ! the current approximation of the strain increment must be returned. If
    ! convergence is attained, these values will be saved and associated with the
    ! converged solution (cf. the input values)

    REAL(KIND=dp), INTENT(OUT) :: DDSDDE(NTENS,NTENS)
    ! The derivative of (Cauchy) stress response function with respect to the
    ! strain evaluated for the current approximation must be returned

    REAL(KIND=dp), INTENT(INOUT) :: SSE, SPD, SCD
    ! Requirement for Elmer: Provide specific strain energy (sse), plastic
    ! dissipation (spd) and creep dissipation (scd) at the previous time/load
    ! level (these are supposed to be declared to be state variables)
    ! Requirement for umat:  The values of the energy variables corresponding to
    ! the current approximation may be returned

    REAL(KIND=dp), INTENT(OUT) :: rpl
    ! The mechanical heating power (volumetric)

    REAL(KIND=dp), INTENT(OUT) :: ddsddt(NTENS), drplde(NTENS), drpldt

    REAL(KIND=dp), INTENT(IN) :: STRAN(NTENS)
    ! This gives the strains before the time/load increment.
    ! The strain can be computed from the deformation gradient, so this
    ! argument can be considered to be redundant. Elmer provides
    ! this information anyway. Abaqus assumes that the logarithmic strain
    ! is used, but Elmer may also use other strain measures.

    REAL(KIND=dp), INTENT(IN) :: DSTRAN(NTENS)
    ! The current candidate for the strain increment to obtain the current
    ! candidate for the stress. In principle this could be computed from the
    ! deformation gradient; cf. the variable stran.

    REAL(KIND=dp), INTENT(IN) :: TIME(2)
    ! Both entries give time before the time/load increment (the time for the last
    ! converged solution

    REAL(KIND=dp), INTENT(IN) :: DTIME
    ! The time increment

    REAL(KIND=dp), INTENT(IN) :: TEMP
    ! Temperature before the time/load increment

    REAL(KIND=dp), INTENT(IN) :: dtemp
    ! Temperature increment associated wíth the time/load increment. Currently
    ! Elmer assumes isothermal conditions during the load increment.

    REAL(KIND=dp), INTENT(IN) :: predef(1), dpred(1)
    ! These are just dummy variables for Elmer

    CHARACTER(len=80), INTENT(IN) :: CMNAME
    ! The material model name

    INTEGER, INTENT(IN) :: NDI
    ! The number of direct stress components

    INTEGER, INTENT(IN) :: NSHR
    ! The number of the engineering shear strain components

    INTEGER, INTENT(IN) :: NTENS
    ! The size of the array containing the stress or strain components

    INTEGER, INTENT(IN) :: NSTATEV
    ! The number of state variables associated with the material model

    REAL(KIND=dp), INTENT(IN) :: PROPS(NPROPS)
    ! An array of material constants

    INTEGER, INTENT(IN) :: NPROPS
    ! The number of the material constants

    REAL(KIND=dp), INTENT(IN) :: coords(3)
    ! The coordinates of the current point could be specified

    REAL(KIND=dp), INTENT(IN) :: drot(3,3)
    ! No support for keeping track of rigid body rotations
    ! (the variable is initialized to the identity)

    REAL(KIND=dp), INTENT(INOUT) :: pnewdt
    ! Currently, suggesting a new size of time increment does not make any impact

    REAL(KIND=dp), INTENT(IN) :: celent
    ! The element size is not yet provided by Elmer

    REAL(KIND=dp), INTENT(IN) :: DFRGRD0(3,3)
    ! The deformation gradient before the time/load increment (at the previous
    ! time/load level for converged solution)

    REAL(KIND=dp), INTENT(IN) :: DFRGRD1(3,3)
    ! The deformation gradient corresponding to the current approximation
    ! (cf. the return value of STRESS variable)

    INTEGER, INTENT(IN) :: NOEL
    ! The element number

    INTEGER, INTENT(IN) :: NPT
    ! The integration point number

    INTEGER, INTENT(IN) :: layer, kspt, kstep, kinc
    ! kstep and kinc could be provided to give information on the incrementation
    ! procedure
!------------------------------------------------------------------------------

    ! Local variables:

!------------------------------------------------------------------------------

    ! ADD THE MATERIAL MODEL DEFINITION HERE TO MAKE THIS FUNCTIONAL

!------------------------------------------------------------------------------
  END SUBROUTINE UMAT_template
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE linear_isotropic(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
       rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
       predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
       coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
       kstep, kinc)
!------------------------------------------------------------------------------
    USE Types
    IMPLICIT NONE

    REAL(KIND=dp), INTENT(INOUT) :: STRESS(NTENS)
    ! Requirement for Elmer: At the time of calling the Cauchy stress T_n before
    ! the time/load increment is given
    ! Requirement for umat:  The stress T_{n+1}^{(k)} corresponding to the
    ! current approximation of the strain increment (DSTRAN) must be returned.
    ! If the strain increment is defined to be zero in the beginning of the
    ! nonlinear iteration, Elmer will generate a candidate for the strain increment
    ! by assuming purely elastic increment characterized by DDSDDE.

    REAL(KIND=dp), INTENT(INOUT) :: STATEV(NSTATEV)
    ! Requirement for Elmer: The state variables Q_n as specified at the
    ! previous time/load level for converged solution are given.
    ! Requirement for umat:  The state variables Q_{n+1}^{(k)} corresponding to
    ! the current approximation of the strain increment must be returned. If
    ! convergence is attained, these values will be saved and associated with the
    ! converged solution (cf. the input values)

    REAL(KIND=dp), INTENT(OUT) :: DDSDDE(NTENS,NTENS)
    ! The derivative of (Cauchy) stress response function with respect to the
    ! strain evaluated for the current approximation must be returned

    REAL(KIND=dp), INTENT(INOUT) :: SSE, SPD, SCD
    ! Requirement for Elmer: Provide specific strain energy (sse), plastic
    ! dissipation (spd) and creep dissipation (scd) at the previous time/load
    ! level (these are supposed to be declared to be state variables)
    ! Requirement for umat:  The values of the energy variables corresponding to
    ! the current approximation may be returned

    REAL(KIND=dp), INTENT(OUT) :: rpl
    ! The mechanical heating power (volumetric)

    REAL(KIND=dp), INTENT(OUT) :: ddsddt(NTENS), drplde(NTENS), drpldt

    REAL(KIND=dp), INTENT(IN) :: STRAN(NTENS)
    ! This gives the strains before the time/load increment.
    ! The strain can be computed from the deformation gradient, so this
    ! argument can be considered to be redundant. Elmer provides
    ! this information anyway. Abaqus assumes that the logarithmic strain
    ! is used, but Elmer may also use other strain measures.

    REAL(KIND=dp), INTENT(IN) :: DSTRAN(NTENS)
    ! The current candidate for the strain increment to obtain the current
    ! candidate for the stress. In principle this could be computed from the
    ! deformation gradient; cf. the variable stran.

    REAL(KIND=dp), INTENT(IN) :: TIME(2)
    ! Both entries give time before the time/load increment (the time for the last
    ! converged solution

    REAL(KIND=dp), INTENT(IN) :: DTIME
    ! The time increment

    REAL(KIND=dp), INTENT(IN) :: TEMP
    ! Temperature before the time/load increment

    REAL(KIND=dp), INTENT(IN) :: dtemp
    ! Temperature increment associated wíth the time/load increment. Currently
    ! Elmer assumes isothermal conditions during the load increment.

    REAL(KIND=dp), INTENT(IN) :: predef(1), dpred(1)
    ! These are just dummy variables for Elmer

    CHARACTER(len=80), INTENT(IN) :: CMNAME
    ! The material model name

    INTEGER, INTENT(IN) :: NDI
    ! The number of direct stress components

    INTEGER, INTENT(IN) :: NSHR
    ! The number of the engineering shear strain components

    INTEGER, INTENT(IN) :: NTENS
    ! The size of the array containing the stress or strain components

    INTEGER, INTENT(IN) :: NSTATEV
    ! The number of state variables associated with the material model

    REAL(KIND=dp), INTENT(IN) :: PROPS(NPROPS)
    ! An array of material constants

    INTEGER, INTENT(IN) :: NPROPS
    ! The number of the material constants

    REAL(KIND=dp), INTENT(IN) :: coords(3)
    ! The coordinates of the current point could be specified

    REAL(KIND=dp), INTENT(IN) :: drot(3,3)
    ! No support for keeping track of rigid body rotations
    ! (the variable is initialized to the identity)

    REAL(KIND=dp), INTENT(INOUT) :: pnewdt
    ! Currently, suggesting a new size of time increment does not make any impact

    REAL(KIND=dp), INTENT(IN) :: celent
    ! The element size is not yet provided by Elmer

    REAL(KIND=dp), INTENT(IN) :: DFRGRD0(3,3)
    ! The deformation gradient before the time/load increment (at the previous
    ! time/load level for converged solution)

    REAL(KIND=dp), INTENT(IN) :: DFRGRD1(3,3)
    ! The deformation gradient corresponding to the current approximation
    ! (cf. the return value of STRESS variable)

    INTEGER, INTENT(IN) :: NOEL
    ! The element number

    INTEGER, INTENT(IN) :: NPT
    ! The integration point number

    INTEGER, INTENT(IN) :: layer, kspt, kstep, kinc
    ! kstep and kinc could be provided to give information on the incrementation
    ! procedure
!------------------------------------------------------------------------------
    ! Local variables:
    INTEGER :: i
    REAL(KIND=dp) :: nu, E, LambdaLame, MuLame
!------------------------------------------------------------------------------

    ! Get Young's modulus and the Poisson ratio:
    E = Props(2)
    nu = Props(3)

    LambdaLame = E * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
    MuLame = E / (2.0d0 * (1.0d0 + nu))

    ddsdde = 0.0d0
    ddsdde(1:ndi,1:ndi) = LambdaLame
    DO i=1,ntens
      ddsdde(i,i) = ddsdde(i,i) + MuLame
    END DO
    DO i=1,ndi
      ddsdde(i,i) = ddsdde(i,i) + MuLame
    END DO
    !
    ! We have a linear response function, so the following update is precise
    ! (no higher-order terms related to the notion of differentiability).
    ! Note that we could also define
    !
    !        stress = stress_response_function(stran + dstran)
    ! or
    !        stress = stress_response_function(dfrgrd1)
    !
    ! which may be the precise definition of the functionality required.
    stress = stress + MATMUL(ddsdde,dstran)
    ! So, for this model, the other way to return the stress:
    !stress = MATMUL(ddsdde,stran+dstran)

!------------------------------------------------------------------------------
  END SUBROUTINE linear_isotropic
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE stvenant_kirchhoff(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
       rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
       predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
       coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
       kstep, kinc)
!------------------------------------------------------------------------------
    USE Types
    IMPLICIT NONE

    REAL(KIND=dp), INTENT(INOUT) :: STRESS(NTENS)
    ! Requirement for Elmer: At the time of calling the Cauchy stress T_n before
    ! the time/load increment is given
    ! Requirement for umat:  The stress T_{n+1}^{(k)} corresponding to the
    ! current approximation of the strain increment (DSTRAN) must be returned.
    ! If the strain increment is defined to be zero in the beginning of the
    ! nonlinear iteration, Elmer will generate a candidate for the strain increment
    ! by assuming purely elastic increment characterized by DDSDDE.

    REAL(KIND=dp), INTENT(INOUT) :: STATEV(NSTATEV)
    ! Requirement for Elmer: The state variables Q_n as specified at the
    ! previous time/load level for converged solution are given.
    ! Requirement for umat:  The state variables Q_{n+1}^{(k)} corresponding to
    ! the current approximation of the strain increment must be returned. If
    ! convergence is attained, these values will be saved and associated with the
    ! converged solution (cf. the input values)

    REAL(KIND=dp), INTENT(OUT) :: DDSDDE(NTENS,NTENS)
    ! The derivative of (Cauchy) stress response function with respect to the
    ! strain evaluated for the current approximation must be returned

    REAL(KIND=dp), INTENT(INOUT) :: SSE, SPD, SCD
    ! Requirement for Elmer: Provide specific strain energy (sse), plastic
    ! dissipation (spd) and creep dissipation (scd) at the previous time/load
    ! level (these are supposed to be declared to be state variables)
    ! Requirement for umat:  The values of the energy variables corresponding to
    ! the current approximation may be returned

    REAL(KIND=dp), INTENT(OUT) :: rpl
    ! The mechanical heating power (volumetric)

    REAL(KIND=dp), INTENT(OUT) :: ddsddt(NTENS), drplde(NTENS), drpldt

    REAL(KIND=dp), INTENT(IN) :: STRAN(NTENS)
    ! This gives the strains before the time/load increment.
    ! The strain can be computed from the deformation gradient, so this
    ! argument can be considered to be redundant. Elmer provides
    ! this information anyway. Abaqus assumes that the logarithmic strain
    ! is used, but Elmer may also use other strain measures.

    REAL(KIND=dp), INTENT(IN) :: DSTRAN(NTENS)
    ! The current candidate for the strain increment to obtain the current
    ! candidate for the stress. In principle this could be computed from the
    ! deformation gradient; cf. the variable stran.

    REAL(KIND=dp), INTENT(IN) :: TIME(2)
    ! Both entries give time before the time/load increment (the time for the last
    ! converged solution

    REAL(KIND=dp), INTENT(IN) :: DTIME
    ! The time increment

    REAL(KIND=dp), INTENT(IN) :: TEMP
    ! Temperature before the time/load increment

    REAL(KIND=dp), INTENT(IN) :: dtemp
    ! Temperature increment associated wíth the time/load increment. Currently
    ! Elmer assumes isothermal conditions during the load increment.

    REAL(KIND=dp), INTENT(IN) :: predef(1), dpred(1)
    ! These are just dummy variables for Elmer

    CHARACTER(len=80), INTENT(IN) :: CMNAME
    ! The material model name

    INTEGER, INTENT(IN) :: NDI
    ! The number of direct stress components

    INTEGER, INTENT(IN) :: NSHR
    ! The number of the engineering shear strain components

    INTEGER, INTENT(IN) :: NTENS
    ! The size of the array containing the stress or strain components

    INTEGER, INTENT(IN) :: NSTATEV
    ! The number of state variables associated with the material model

    REAL(KIND=dp), INTENT(IN) :: PROPS(NPROPS)
    ! An array of material constants

    INTEGER, INTENT(IN) :: NPROPS
    ! The number of the material constants

    REAL(KIND=dp), INTENT(IN) :: coords(3)
    ! The coordinates of the current point could be specified

    REAL(KIND=dp), INTENT(IN) :: drot(3,3)
    ! No support for keeping track of rigid body rotations
    ! (the variable is initialized to the identity)

    REAL(KIND=dp), INTENT(INOUT) :: pnewdt
    ! Currently, suggesting a new size of time increment does not make any impact

    REAL(KIND=dp), INTENT(IN) :: celent
    ! The element size is not yet provided by Elmer

    REAL(KIND=dp), INTENT(IN) :: DFRGRD0(3,3)
    ! The deformation gradient before the time/load increment (at the previous
    ! time/load level for converged solution)

    REAL(KIND=dp), INTENT(IN) :: DFRGRD1(3,3)
    ! The deformation gradient corresponding to the current approximation
    ! (cf. the return value of STRESS variable)

    INTEGER, INTENT(IN) :: NOEL
    ! The element number

    INTEGER, INTENT(IN) :: NPT
    ! The integration point number

    INTEGER, INTENT(IN) :: layer, kspt, kstep, kinc
    ! kstep and kinc could be provided to give information on the incrementation
    ! procedure
!------------------------------------------------------------------------------

    INTEGER :: i, j, k

    REAL(KIND=dp) :: SymBasis(6,3,3)
    REAL(KIND=dp) :: Identity(3,3), B(3,3), C(3,3), Strain(3,3), S(3,3), Sigma(3,3)
    REAL(KIND=dp) :: WorkMat(3,3)
    REAL(KIND=dp) :: StrainVec(ntens), Stress2(ntens)
    REAL(KIND=dp) :: DetDefG
    REAL(KIND=dp) :: nu, E, LambdaLame, MuLame

!------------------------------------------------------------------------------

    SymBasis(1,1:3,1:3) = RESHAPE((/ 1,0,0,0,0,0,0,0,0 /),(/ 3,3 /))
    SymBasis(2,1:3,1:3) = RESHAPE((/ 0,0,0,0,1,0,0,0,0 /),(/ 3,3 /))
    SymBasis(3,1:3,1:3) = RESHAPE((/ 0,0,0,0,0,0,0,0,1 /),(/ 3,3 /))
    SymBasis(4,1:3,1:3) = RESHAPE((/ 0.0d0,0.5d0,0.0d0,0.5d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0 /),(/ 3,3 /))
    SymBasis(5,1:3,1:3) = RESHAPE((/ 0.0d0,0.0d0,0.5d0,0.0d0,0.0d0,0.0d0,0.5d0,0.0d0,0.0d0 /),(/ 3,3 /))
    SymBasis(6,1:3,1:3) = RESHAPE((/ 0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.5d0,0.0d0,0.5d0,0.0d0 /),(/ 3,3 /))
    Identity(1:3,1:3) = RESHAPE((/ 1,0,0,0,1,0,0,0,1 /),(/ 3,3 /))

    B = MATMUL(dfrgrd1, TRANSPOSE(dfrgrd1))
    C = MATMUL(TRANSPOSE(dfrgrd1), dfrgrd1)
    ! This example uses the Lagrangian (Green-St Venant) strain tensor:
    Strain = 0.5d0 * (C - Identity)

    DO i=1,ndi
      StrainVec(i) = Strain(i,i)
    END DO
    DO i=1,nshr
      SELECT CASE(i)
      CASE(1)
        StrainVec(ndi+i) = Strain(1,2)+Strain(2,1)
      CASE(2)
        StrainVec(ndi+i) = Strain(1,3)+Strain(3,1)
      CASE(3)
        StrainVec(ndi+i) = Strain(2,3)+Strain(3,2)
      END SELECT
    END DO

    DetDefG = Dfrgrd1(1,1) * ( Dfrgrd1(2,2)*Dfrgrd1(3,3) - Dfrgrd1(2,3)*Dfrgrd1(3,2) ) + &
        Dfrgrd1(1,2) * ( Dfrgrd1(2,3)*Dfrgrd1(3,1) - Dfrgrd1(2,1)*Dfrgrd1(3,3) ) + &
        Dfrgrd1(1,3) * ( Dfrgrd1(2,1)*Dfrgrd1(3,2) - Dfrgrd1(2,2)*Dfrgrd1(3,1) )

    ! Get Young's modulus and the Poisson ratio:
    E = Props(2)
    nu = Props(3)

    LambdaLame = E * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
    MuLame = E / (2.0d0 * (1.0d0 + nu))

    ! --------------------------------------------------------------------------------
    ! Here we compute the current stress directly by using the
    ! supplied deformation gradient, so that the strain increment is not used.
    ! In addition, since it seems that the exact differentiation of the response function
    ! for the Cauchy stress cannot be done in a straightforward manner, we now make only
    ! a partial approximation. The Cauchy stress is given by
    !
    !     sigma(F) = 1/det(F) F S(E(F)) F^T
    !
    ! We however consider only the depedence on the strain as
    !
    !     sigma(.) = 1/det(F) F S(.) F^T
    !
    ! This simplification makes the nonlinear iteration to be an inexact Newton
    ! method whose performance may deteriorate for large strains. If the convergence is
    ! attained, the solution nevertheless obeys the St. Venant-Kirchhoff law since
    ! there are no approximations in the computation of the residual.
    ! --------------------------------------------------------------------------------
    ! The constitutive matrix relating the second Piola-Kirchhoff stress and
    ! the strain tensor:
    ddsdde = 0.0d0
    ddsdde(1:ndi,1:ndi) = LambdaLame
    DO i=1,ntens
      ddsdde(i,i) = ddsdde(i,i) + MuLame
    END DO
    DO i=1,ndi
      ddsdde(i,i) = ddsdde(i,i) + MuLame
    END DO
    Stress2 = MATMUL(ddsdde,StrainVec)

    ! The second Piola-Kirchhoff stress in the tensor form:
    S = 0.0d0
    DO i=1,ndi
      S = S + Stress2(i)*SymBasis(i,:,:)
    END DO
    DO i=1,nshr
      S = S + 2.0d0 * Stress2(ndi+i) * SymBasis(ndi+i,:,:)
    END DO

    ! The Cauchy stress tensor:
    Sigma = 1.0d0/DetDefG * MATMUL(dfrgrd1, MATMUL(S,TRANSPOSE(dfrgrd1)))

    DO i=1,ndi
      Stress(i) = Sigma(i,i)
    END DO
    DO i=1,nshr
      SELECT CASE(i)
      CASE(1)
        Stress(ndi+i) = Sigma(1,2)
      CASE(2)
        Stress(ndi+i) = Sigma(1,3)
      CASE(3)
        Stress(ndi+i) = Sigma(2,3)
      END SELECT
    END DO

    ! The derivative: The part corresponding to lambda * tr(E) I
    ddsdde = 0.0d0
    WorkMat = LambdaLame * 1/DetDefG * B
    DO i=1,ndi
      DO j=1,ndi
        ddsdde(j,i) = ddsdde(j,i) + WorkMat(j,j)
      END DO
      DO j=1,nshr
        SELECT CASE(j)
        CASE(1)
          ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(1,2)
        CASE(2)
          ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(1,3)
        CASE(3)
          ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(2,3)
        END SELECT
      END DO
    END DO

    ! The rest corresponding to  2 * mu * E
    DO i=1,ndi
      WorkMat = 2.0d0 * MuLame * 1/DetDefG * MATMUL(dfrgrd1, MATMUL(SymBasis(i,:,:), TRANSPOSE(dfrgrd1)))
      DO j=1,ndi
        ddsdde(j,i) = ddsdde(j,i) + WorkMat(j,j)
      END DO
      DO j=1,nshr
        SELECT CASE(j)
        CASE(1)
          ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(1,2)
        CASE(2)
          ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(1,3)
        CASE(3)
          ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(2,3)
        END SELECT
      END DO
    END DO

    DO i=1,nshr
      WorkMat = 2.0d0 * MuLame * 1/DetDefG * MATMUL(dfrgrd1, MATMUL(SymBasis(ndi+i,:,:), TRANSPOSE(dfrgrd1)))
      DO j=1,ndi
        ddsdde(j,ndi+i) = ddsdde(j,ndi+i) + 1.0d0 * WorkMat(j,j)
      END DO
      DO j=1,nshr
        SELECT CASE(j)
        CASE(1)
          ddsdde(ndi+j,ndi+i) = ddsdde(ndi+j,ndi+i) + 1.0d0 * WorkMat(1,2)
        CASE(2)
          ddsdde(ndi+j,ndi+i) = ddsdde(ndi+j,ndi+i) + 1.0d0 * WorkMat(1,3)
        CASE(3)
          ddsdde(ndi+j,ndi+i) = ddsdde(ndi+j,ndi+i) + 1.0d0 * WorkMat(2,3)
        END SELECT
      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE stvenant_kirchhoff
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE hencky_stvenant_kirchhoff(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
       rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
       predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
       coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
       kstep, kinc)
!------------------------------------------------------------------------------
    USE Types
    USE Messages
    IMPLICIT NONE

    REAL(KIND=dp), INTENT(INOUT) :: STRESS(NTENS)
    ! Requirement for Elmer: At the time of calling the Cauchy stress T_n before
    ! the time/load increment is given
    ! Requirement for umat:  The stress T_{n+1}^{(k)} corresponding to the
    ! current approximation of the strain increment (DSTRAN) must be returned.
    ! If the strain increment is defined to be zero in the beginning of the
    ! nonlinear iteration, Elmer will generate a candidate for the strain increment
    ! by assuming purely elastic increment characterized by DDSDDE.

    REAL(KIND=dp), INTENT(INOUT) :: STATEV(NSTATEV)
    ! Requirement for Elmer: The state variables Q_n as specified at the
    ! previous time/load level for converged solution are given.
    ! Requirement for umat:  The state variables Q_{n+1}^{(k)} corresponding to
    ! the current approximation of the strain increment must be returned. If
    ! convergence is attained, these values will be saved and associated with the
    ! converged solution (cf. the input values)

    REAL(KIND=dp), INTENT(OUT) :: DDSDDE(NTENS,NTENS)
    ! The derivative of (Cauchy) stress response function with respect to the
    ! strain evaluated for the current approximation must be returned

    REAL(KIND=dp), INTENT(INOUT) :: SSE, SPD, SCD
    ! Requirement for Elmer: Provide specific strain energy (sse), plastic
    ! dissipation (spd) and creep dissipation (scd) at the previous time/load
    ! level (these are supposed to be declared to be state variables)
    ! Requirement for umat:  The values of the energy variables corresponding to
    ! the current approximation may be returned

    REAL(KIND=dp), INTENT(OUT) :: rpl
    ! The mechanical heating power (volumetric)

    REAL(KIND=dp), INTENT(OUT) :: ddsddt(NTENS), drplde(NTENS), drpldt

    REAL(KIND=dp), INTENT(IN) :: STRAN(NTENS)
    ! This gives the strains before the time/load increment.
    ! The strain can be computed from the deformation gradient, so this
    ! argument can be considered to be redundant. Elmer provides
    ! this information anyway. Abaqus assumes that the logarithmic strain
    ! is used, but Elmer may also use other strain measures.

    REAL(KIND=dp), INTENT(IN) :: DSTRAN(NTENS)
    ! The current candidate for the strain increment to obtain the current
    ! candidate for the stress. In principle this could be computed from the
    ! deformation gradient; cf. the variable stran.

    REAL(KIND=dp), INTENT(IN) :: TIME(2)
    ! Both entries give time before the time/load increment (the time for the last
    ! converged solution

    REAL(KIND=dp), INTENT(IN) :: DTIME
    ! The time increment

    REAL(KIND=dp), INTENT(IN) :: TEMP
    ! Temperature before the time/load increment

    REAL(KIND=dp), INTENT(IN) :: dtemp
    ! Temperature increment associated wíth the time/load increment. Currently
    ! Elmer assumes isothermal conditions during the load increment.

    REAL(KIND=dp), INTENT(IN) :: predef(1), dpred(1)
    ! These are just dummy variables for Elmer

    CHARACTER(len=80), INTENT(IN) :: CMNAME
    ! The material model name

    INTEGER, INTENT(IN) :: NDI
    ! The number of direct stress components

    INTEGER, INTENT(IN) :: NSHR
    ! The number of the engineering shear strain components

    INTEGER, INTENT(IN) :: NTENS
    ! The size of the array containing the stress or strain components

    INTEGER, INTENT(IN) :: NSTATEV
    ! The number of state variables associated with the material model

    REAL(KIND=dp), INTENT(IN) :: PROPS(NPROPS)
    ! An array of material constants

    INTEGER, INTENT(IN) :: NPROPS
    ! The number of the material constants

    REAL(KIND=dp), INTENT(IN) :: coords(3)
    ! The coordinates of the current point could be specified

    REAL(KIND=dp), INTENT(IN) :: drot(3,3)
    ! No support for keeping track of rigid body rotations
    ! (the variable is initialized to the identity)

    REAL(KIND=dp), INTENT(INOUT) :: pnewdt
    ! Currently, suggesting a new size of time increment does not make any impact

    REAL(KIND=dp), INTENT(IN) :: celent
    ! The element size is not yet provided by Elmer

    REAL(KIND=dp), INTENT(IN) :: DFRGRD0(3,3)
    ! The deformation gradient before the time/load increment (at the previous
    ! time/load level for converged solution)

    REAL(KIND=dp), INTENT(IN) :: DFRGRD1(3,3)
    ! The deformation gradient corresponding to the current approximation
    ! (cf. the return value of STRESS variable)

    INTEGER, INTENT(IN) :: NOEL
    ! The element number

    INTEGER, INTENT(IN) :: NPT
    ! The integration point number

    INTEGER, INTENT(IN) :: layer, kspt, kstep, kinc
    ! kstep and kinc could be provided to give information on the incrementation
    ! procedure
!------------------------------------------------------------------------------

    INTEGER :: i, j, k
    INTEGER :: PriLWork=102, PriInfo

    REAL(KIND=dp) :: SymBasis(6,3,3)
    REAL(KIND=dp) :: Identity(3,3), B(3,3), C(3,3), Strain(3,3), S(3,3), Sigma(3,3)
    REAL(KIND=dp) :: WorkMat(3,3)
    REAL(KIND=dp) :: StrainVec(ntens), Stress2(ntens)
    REAL(KIND=dp) :: DetDefG
    REAL(KIND=dp) :: nu, E, LambdaLame, MuLame
    REAL(KIND=dp) :: EigenVals(3), PriWork(102)
!------------------------------------------------------------------------------

    SymBasis(1,1:3,1:3) = RESHAPE((/ 1,0,0,0,0,0,0,0,0 /),(/ 3,3 /))
    SymBasis(2,1:3,1:3) = RESHAPE((/ 0,0,0,0,1,0,0,0,0 /),(/ 3,3 /))
    SymBasis(3,1:3,1:3) = RESHAPE((/ 0,0,0,0,0,0,0,0,1 /),(/ 3,3 /))
    SymBasis(4,1:3,1:3) = RESHAPE((/ 0.0d0,0.5d0,0.0d0,0.5d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0 /),(/ 3,3 /))
    SymBasis(5,1:3,1:3) = RESHAPE((/ 0.0d0,0.0d0,0.5d0,0.0d0,0.0d0,0.0d0,0.5d0,0.0d0,0.0d0 /),(/ 3,3 /))
    SymBasis(6,1:3,1:3) = RESHAPE((/ 0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.5d0,0.0d0,0.5d0,0.0d0 /),(/ 3,3 /))
    Identity(1:3,1:3) = RESHAPE((/ 1,0,0,0,1,0,0,0,1 /),(/ 3,3 /))

    B = MATMUL(dfrgrd1, TRANSPOSE(dfrgrd1))
    C = MATMUL(TRANSPOSE(dfrgrd1), dfrgrd1)

    ! -----------------------------------------------------------
    ! Compute the spectral decomposition of C
    ! -----------------------------------------------------------
    DO i=1,3
      k = i
      DO j=k,3
        WorkMat(i,j) = C(i,j)
      END DO
    END DO
    CALL DSYEV('V', 'U', 3, WorkMat, 3, EigenVals, PriWork, PriLWork, PriInfo)
    IF (PriInfo /= 0) THEN
      CALL Fatal( 'UMAT', 'DSYEV cannot generate eigen basis')
    END IF

    Strain = 0.0d0
    Strain(1,1) = LOG(SQRT(EigenVals(1)))
    Strain(2,2) = LOG(SQRT(EigenVals(2)))
    Strain(3,3) = LOG(SQRT(EigenVals(3)))
    ! Transform back to the original coordinates:
    Strain = MATMUL(WorkMat, MATMUL(Strain,TRANSPOSE(WorkMat)))

    DO i=1,ndi
      StrainVec(i) = Strain(i,i)
    END DO
    DO i=1,nshr
      SELECT CASE(i)
      CASE(1)
        StrainVec(ndi+i) = Strain(1,2)+Strain(2,1)
      CASE(2)
        StrainVec(ndi+i) = Strain(1,3)+Strain(3,1)
      CASE(3)
        StrainVec(ndi+i) = Strain(2,3)+Strain(3,2)
      END SELECT
    END DO

    DetDefG = Dfrgrd1(1,1) * ( Dfrgrd1(2,2)*Dfrgrd1(3,3) - Dfrgrd1(2,3)*Dfrgrd1(3,2) ) + &
        Dfrgrd1(1,2) * ( Dfrgrd1(2,3)*Dfrgrd1(3,1) - Dfrgrd1(2,1)*Dfrgrd1(3,3) ) + &
        Dfrgrd1(1,3) * ( Dfrgrd1(2,1)*Dfrgrd1(3,2) - Dfrgrd1(2,2)*Dfrgrd1(3,1) )

    ! Get Young's modulus and the Poisson ratio:
    E = Props(2)
    nu = Props(3)

    LambdaLame = E * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
    MuLame = E / (2.0d0 * (1.0d0 + nu))

    ! --------------------------------------------------------------------------------
    ! Here we compute the current stress directly by using the
    ! supplied deformation gradient, so that the strain increment is not used.
    ! In addition, since it seems that the exact differentiation of the response function
    ! for the Cauchy stress cannot be done in a straightforward manner, we now make only
    ! a partial approximation. The Cauchy stress is given by
    !
    !     sigma(F) = 1/det(F) F S(E(F)) F^T
    !
    ! We however consider only the depedence on the strain as
    !
    !     sigma(.) = 1/det(F) F S(.) F^T
    !
    ! This simplification makes the nonlinear iteration to be an inexact Newton
    ! method whose performance may deteriorate for large strains. If the convergence is
    ! attained, the solution nevertheless obeys the St. Venant-Kirchhoff law since
    ! there are no approximations in the computation of the residual.
    ! --------------------------------------------------------------------------------
    ! The constitutive matrix relating the second Piola-Kirchhoff stress and
    ! the strain tensor:
    ddsdde = 0.0d0
    ddsdde(1:ndi,1:ndi) = LambdaLame
    DO i=1,ntens
      ddsdde(i,i) = ddsdde(i,i) + MuLame
    END DO
    DO i=1,ndi
      ddsdde(i,i) = ddsdde(i,i) + MuLame
    END DO
    Stress2 = MATMUL(ddsdde,StrainVec)

    ! The second Piola-Kirchhoff stress in the tensor form:
    S = 0.0d0
    DO i=1,ndi
      S = S + Stress2(i)*SymBasis(i,:,:)
    END DO
    DO i=1,nshr
      S = S + 2.0d0 * Stress2(ndi+i) * SymBasis(ndi+i,:,:)
    END DO

    ! The Cauchy stress tensor:
    Sigma = 1.0d0/DetDefG * MATMUL(dfrgrd1, MATMUL(S,TRANSPOSE(dfrgrd1)))

    DO i=1,ndi
      Stress(i) = Sigma(i,i)
    END DO
    DO i=1,nshr
      SELECT CASE(i)
      CASE(1)
        Stress(ndi+i) = Sigma(1,2)
      CASE(2)
        Stress(ndi+i) = Sigma(1,3)
      CASE(3)
        Stress(ndi+i) = Sigma(2,3)
      END SELECT
    END DO

    ! The derivative: The part corresponding to lambda * tr(E) I
    ddsdde = 0.0d0
    WorkMat = LambdaLame * 1/DetDefG * B
    DO i=1,ndi
      DO j=1,ndi
        ddsdde(j,i) = ddsdde(j,i) + WorkMat(j,j)
      END DO
      DO j=1,nshr
        SELECT CASE(j)
        CASE(1)
          ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(1,2)
        CASE(2)
          ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(1,3)
        CASE(3)
          ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(2,3)
        END SELECT
      END DO
    END DO

    ! The rest corresponding to  2 * mu * E
    DO i=1,ndi
      WorkMat = 2.0d0 * MuLame * 1/DetDefG * MATMUL(dfrgrd1, MATMUL(SymBasis(i,:,:), TRANSPOSE(dfrgrd1)))
      DO j=1,ndi
        ddsdde(j,i) = ddsdde(j,i) + WorkMat(j,j)
      END DO
      DO j=1,nshr
        SELECT CASE(j)
        CASE(1)
          ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(1,2)
        CASE(2)
          ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(1,3)
        CASE(3)
          ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(2,3)
        END SELECT
      END DO
    END DO

    DO i=1,nshr
      WorkMat = 2.0d0 * MuLame * 1/DetDefG * MATMUL(dfrgrd1, MATMUL(SymBasis(ndi+i,:,:), TRANSPOSE(dfrgrd1)))
      DO j=1,ndi
        ddsdde(j,ndi+i) = ddsdde(j,ndi+i) + 1.0d0 * WorkMat(j,j)
      END DO
      DO j=1,nshr
        SELECT CASE(j)
        CASE(1)
          ddsdde(ndi+j,ndi+i) = ddsdde(ndi+j,ndi+i) + 1.0d0 * WorkMat(1,2)
        CASE(2)
          ddsdde(ndi+j,ndi+i) = ddsdde(ndi+j,ndi+i) + 1.0d0 * WorkMat(1,3)
        CASE(3)
          ddsdde(ndi+j,ndi+i) = ddsdde(ndi+j,ndi+i) + 1.0d0 * WorkMat(2,3)
        END SELECT
      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE hencky_stvenant_kirchhoff
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! UMATplastic_v2.F90  --  small-strain rate-independent J2 plasticity for Elmer
!
! Replaces the bi_linear and plastic routines of UMATplastic.F90 (hsk_bi is
! NOT included here, keep it in its own file).
!
! Model: additive strain split, von Mises yield, isotropic hardening, radial
! return (exact for a piecewise-linear hardening curve, no iteration needed),
! consistent algorithmic tangent returned in DDSDDE.
!
! The input curve is the uniaxial stress-TOTAL strain curve, as before; it is
! converted internally to yield stress vs equivalent plastic strain.
!   Young's modulus = stress1/strain1 (first point is the end of elasticity)
!   Beyond the last point the material is perfectly plastic.
!
! plastic   : PROPS = nu, strain1, stress1, strain2, stress2, ... (max 100 pts)
!
! STATEV (NSTATEV must be >= 8): Use OUTPUT STATE VARIABLES = LOGICAL TRUE
! for Paraview they import as cell data not point data
!   1-3 principal strains (max,mid,min)   4-6 principal stresses (max,mid,min)
!   7   von Mises stress                  8   equivalent plastic strain (history)
!
! Component order (Abaqus/Elmer): 11,22,33,12,[13,23]; shear strains are
! engineering strains. NDI must be 3 (3D, axisymmetric, plane strain).
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Eigenvalues of a symmetric 3x3 tensor, sorted lam(1)>=lam(2)>=lam(3).
! v = (a11,a22,a33,a12,a13,a23)
!------------------------------------------------------------------------------
SUBROUTINE pl_sym3_eig(v, lam)
  USE Types
  IMPLICIT NONE
  REAL(KIND=dp), INTENT(IN)  :: v(6)
  REAL(KIND=dp), INTENT(OUT) :: lam(3)
  REAL(KIND=dp) :: q, p1, p2, pp, r, phi, scale, detb
  REAL(KIND=dp) :: b11, b22, b33, b12, b13, b23
  REAL(KIND=dp), PARAMETER :: pl_pi = 3.14159265358979324_dp

  scale = MAXVAL(ABS(v))
  IF (scale == 0.0_dp) THEN
    lam = 0.0_dp
    RETURN
  END IF
  q  = (v(1)+v(2)+v(3))/3.0_dp
  p1 = v(4)**2 + v(5)**2 + v(6)**2
  p2 = (v(1)-q)**2 + (v(2)-q)**2 + (v(3)-q)**2 + 2.0_dp*p1
  IF (p2 <= 1.0e-24_dp*scale**2) THEN
    lam = q
    RETURN
  END IF
  pp  = SQRT(p2/6.0_dp)
  b11 = (v(1)-q)/pp; b22 = (v(2)-q)/pp; b33 = (v(3)-q)/pp
  b12 = v(4)/pp;     b13 = v(5)/pp;     b23 = v(6)/pp
  detb = b11*(b22*b33-b23*b23) - b12*(b12*b33-b23*b13) + b13*(b12*b23-b22*b13)
  r   = MIN(1.0_dp, MAX(-1.0_dp, 0.5_dp*detb))
  phi = ACOS(r)/3.0_dp
  lam(1) = q + 2.0_dp*pp*COS(phi)
  lam(3) = q + 2.0_dp*pp*COS(phi + 2.0_dp*pl_pi/3.0_dp)
  lam(2) = 3.0_dp*q - lam(1) - lam(3)
END SUBROUTINE pl_sym3_eig

!------------------------------------------------------------------------------
! Convert uniaxial (total strain, stress) points into yield stress Yk versus
! equivalent plastic strain pk, with hardening slope Hk on each segment.
! The last segment is perfectly plastic.  ierr /= 0 on invalid input.
!------------------------------------------------------------------------------
SUBROUTINE pl_table_setup(n, eps, sig, E, pk, Yk, Hk, ierr)
  USE Types
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL(KIND=dp), INTENT(IN)  :: eps(n), sig(n)
  REAL(KIND=dp), INTENT(OUT) :: E, pk(n), Yk(n), Hk(n)
  INTEGER, INTENT(OUT) :: ierr
  INTEGER :: i

  ierr = 0
  IF (n < 1 .OR. eps(1) <= 0.0_dp .OR. sig(1) <= 0.0_dp) THEN
    ierr = 1; RETURN
  END IF
  E = sig(1)/eps(1)
  DO i=1,n
    pk(i) = eps(i) - sig(i)/E
    Yk(i) = sig(i)
  END DO
  pk(1) = 0.0_dp
  DO i=1,n-1
    IF (pk(i+1) <= pk(i)) THEN
      ierr = 2; RETURN              ! slope not below E (or strain not increasing)
    END IF
    Hk(i) = (Yk(i+1)-Yk(i))/(pk(i+1)-pk(i))
    IF (Hk(i) < 0.0_dp) THEN
      ierr = 3; RETURN              ! softening not supported
    END IF
  END DO
  Hk(n) = 0.0_dp
END SUBROUTINE pl_table_setup

!------------------------------------------------------------------------------
! Stress update (radial return) and consistent tangent.
! stress: T_n on entry, T_{n+1} on exit.
!------------------------------------------------------------------------------
SUBROUTINE pl_j2_update(ntens, ndi, E, nu, npts, pk, Yk, Hk, stress, dstran, &
     pold, pnew, ddsdde, dissip)
  USE Types
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ntens, ndi, npts
  REAL(KIND=dp), INTENT(IN) :: E, nu, pk(npts), Yk(npts), Hk(npts)
  REAL(KIND=dp), INTENT(IN) :: dstran(ntens), pold
  REAL(KIND=dp), INTENT(INOUT) :: stress(ntens)
  REAL(KIND=dp), INTENT(OUT) :: pnew, ddsdde(ntens,ntens), dissip

  REAL(KIND=dp) :: G, Kb, lam, hyd, sn2, qtr, Yold, Ynew, Hnew, theta, thetab
  REAL(KIND=dp) :: sig(ntens), s(ntens), nvec(ntens), pcand
  INTEGER :: i, j, kseg

  G   = E/(2.0_dp*(1.0_dp+nu))
  Kb  = E/(3.0_dp*(1.0_dp-2.0_dp*nu))
  lam = Kb - 2.0_dp*G/3.0_dp

  ! Elastic tangent
  ddsdde = 0.0_dp
  ddsdde(1:ndi,1:ndi) = lam
  DO i=1,ndi
    ddsdde(i,i) = lam + 2.0_dp*G
  END DO
  DO i=ndi+1,ntens
    ddsdde(i,i) = G
  END DO

  ! Elastic predictor
  sig = stress + MATMUL(ddsdde, dstran)
  hyd = SUM(sig(1:ndi))/3.0_dp
  s = sig
  s(1:ndi) = s(1:ndi) - hyd
  sn2 = SUM(s(1:ndi)**2) + 2.0_dp*SUM(s(ndi+1:ntens)**2)   ! s:s (tensor shears)
  qtr = SQRT(1.5_dp*sn2)

  ! Current segment of the hardening curve and yield stress at p_n
  kseg = 1
  DO i=2,npts
    IF (pk(i) <= pold) kseg = i
  END DO
  Yold = Yk(kseg) + Hk(kseg)*(pold-pk(kseg))

  pnew = pold
  dissip = 0.0_dp
  stress = sig
  IF (qtr <= Yold*(1.0_dp+1.0e-10_dp)) RETURN        ! elastic step

  ! Plastic corrector: solve  qtr - 3G (p - pold) = Y(p)  segment by segment
  pcand = pold
  DO i=kseg,npts
    pcand = (qtr + 3.0_dp*G*pold - Yk(i) + Hk(i)*pk(i))/(3.0_dp*G + Hk(i))
    kseg = i
    IF (i == npts) EXIT
    IF (pcand <= pk(i+1)) EXIT
  END DO
  pnew = MAX(pcand, pold)
  Hnew = Hk(kseg)
  Ynew = Yk(kseg) + Hnew*(pnew-pk(kseg))

  theta  = Ynew/qtr
  thetab = 1.0_dp/(1.0_dp + Hnew/(3.0_dp*G))
  nvec   = s/SQRT(sn2)

  stress = theta*s
  stress(1:ndi) = stress(1:ndi) + hyd
  dissip = 0.5_dp*(Yold+Ynew)*(pnew-pold)

  ! Consistent tangent: K 1x1 + 2G theta I_dev - 2G (theta+thetab-1) n x n
  ddsdde = 0.0_dp
  ddsdde(1:ndi,1:ndi) = Kb - 2.0_dp*G*theta/3.0_dp
  DO i=1,ndi
    ddsdde(i,i) = ddsdde(i,i) + 2.0_dp*G*theta
  END DO
  DO i=ndi+1,ntens
    ddsdde(i,i) = G*theta
  END DO
  DO j=1,ntens
    DO i=1,ntens
      ddsdde(i,j) = ddsdde(i,j) - 2.0_dp*G*(theta+thetab-1.0_dp)*nvec(i)*nvec(j)
    END DO
  END DO
END SUBROUTINE pl_j2_update

!------------------------------------------------------------------------------
! Diagnostic state variables 1-7 (principal values, von Mises)
!------------------------------------------------------------------------------
SUBROUTINE pl_statev(ntens, ndi, stress, strain, statev)
  USE Types
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ntens, ndi
  REAL(KIND=dp), INTENT(IN) :: stress(ntens), strain(ntens)
  REAL(KIND=dp), INTENT(INOUT) :: statev(*)
  REAL(KIND=dp) :: v(6), lam(3), hyd, dev(3)

  v = 0.0_dp
  v(1:ndi) = stress(1:ndi)
  v(4:ntens) = stress(4:ntens)
  hyd = (v(1)+v(2)+v(3))/3.0_dp
  dev = v(1:3) - hyd
  statev(7) = SQRT(1.5_dp*(SUM(dev**2) + 2.0_dp*SUM(v(4:6)**2)))
  CALL pl_sym3_eig(v, lam)
  statev(4:6) = lam

  v = 0.0_dp
  v(1:ndi) = strain(1:ndi)
  v(4:ntens) = 0.5_dp*strain(4:ntens)        ! engineering -> tensor shear
  CALL pl_sym3_eig(v, lam)
  statev(1:3) = lam
END SUBROUTINE pl_statev

!------------------------------------------------------------------------------
! Common driver
!------------------------------------------------------------------------------
SUBROUTINE pl_drive(ntens, ndi, nstatev, nu, npts, eps, sig, stress, statev, &
     stran, dstran, ddsdde, spd)
  USE Types
  USE Messages
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ntens, ndi, nstatev, npts
  REAL(KIND=dp), INTENT(IN) :: nu, eps(npts), sig(npts), stran(ntens), dstran(ntens)
  REAL(KIND=dp), INTENT(INOUT) :: stress(ntens), statev(nstatev), spd
  REAL(KIND=dp), INTENT(OUT) :: ddsdde(ntens,ntens)

  REAL(KIND=dp) :: E, pk(npts), Yk(npts), Hk(npts), pnew, dissip
  INTEGER :: ierr

  IF (ndi /= 3 .OR. (ntens /= 4 .AND. ntens /= 6)) &
      CALL Fatal('UMATplastic', 'Only NDI=3 (NTENS=4 or 6) is supported')
  IF (nstatev < 8) &
      CALL Fatal('UMATplastic', 'At least 8 state variables are required')
  IF (nu <= -1.0_dp .OR. nu >= 0.5_dp) &
      CALL Fatal('UMATplastic', 'Poisson ratio must satisfy -1 < nu < 0.5')

  CALL pl_table_setup(npts, eps, sig, E, pk, Yk, Hk, ierr)
  IF (ierr == 1) CALL Fatal('UMATplastic', 'First strain/stress point must be positive')
  IF (ierr == 2) CALL Fatal('UMATplastic', 'Strains must increase and every tangent slope must be below E')
  IF (ierr == 3) CALL Fatal('UMATplastic', 'Strain softening is not supported')

  CALL pl_j2_update(ntens, ndi, E, nu, npts, pk, Yk, Hk, stress, dstran, &
       statev(8), pnew, ddsdde, dissip)
  statev(8) = pnew
  spd = spd + dissip
  CALL pl_statev(ntens, ndi, stress, stran+dstran, statev)
END SUBROUTINE pl_drive
!------------------------------------------------------------------------------
SUBROUTINE plastic(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
     rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
     predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
     coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
     kstep, kinc)
!------------------------------------------------------------------------------
  USE Types
  USE Messages
  IMPLICIT NONE
  REAL(KIND=dp), INTENT(INOUT) :: STRESS(NTENS), STATEV(NSTATEV), SSE, SPD, SCD
  REAL(KIND=dp), INTENT(OUT) :: DDSDDE(NTENS,NTENS), rpl, ddsddt(NTENS), drplde(NTENS), drpldt
  REAL(KIND=dp), INTENT(IN) :: STRAN(NTENS), DSTRAN(NTENS), TIME(2), DTIME, TEMP, dtemp
  REAL(KIND=dp), INTENT(IN) :: predef(1), dpred(1), PROPS(NPROPS), coords(3), drot(3,3)
  REAL(KIND=dp), INTENT(INOUT) :: pnewdt
  REAL(KIND=dp), INTENT(IN) :: celent, DFRGRD0(3,3), DFRGRD1(3,3)
  CHARACTER(len=80), INTENT(IN) :: CMNAME
  INTEGER, INTENT(IN) :: NDI, NSHR, NTENS, NSTATEV, NPROPS, NOEL, NPT, layer, kspt, kstep, kinc

  REAL(KIND=dp) :: eps(100), sig(100)
  INTEGER :: i, npts

  rpl = 0.0_dp; ddsddt = 0.0_dp; drplde = 0.0_dp; drpldt = 0.0_dp

  npts = (nprops-1)/2
  IF (npts < 1 .OR. npts > 100) &
      CALL Fatal('plastic', 'Need 1..100 strain/stress pairs after nu')
  DO i=1,npts
    eps(i) = PROPS(2*i)
    sig(i) = PROPS(2*i+1)
  END DO

  CALL pl_drive(ntens, ndi, nstatev, PROPS(1), npts, eps(1:npts), sig(1:npts), &
       stress, statev, stran, dstran, ddsdde, spd)
END SUBROUTINE plastic
!------------------------------------------------------------------------------
! UMATplastic_finite.F90  --  finite-strain J2 plasticity for Elmer
!
! Entry points:  fbi_linear   (bilinear table, same 6 constants as bi_linear)
!                fplastic     (table of up to 100 points, same layout as plastic)
!
! Model (Eterovic-Bathe / Simo): multiplicative split F = Fe Fp, Hencky (log)
! elastic law on the elastic left Cauchy-Green tensor be, von Mises yield,
! isotropic hardening, exact radial return in principal log-strain space.
! Everything is computed from DFRGRD0 and DFRGRD1 and the converged state
! variables, so it does not depend on which strain measure Elmer puts in
! STRAN/DSTRAN.  It is objective (rigid rotations rotate the stress).
!
! INPUT CURVE: TRUE stress against LOGARITHMIC TOTAL strain, uniaxial
!   (log strain = ln(1+e);  true stress = s (1+e)  from engineering e, s).
!   If you have true stress against plastic strain (as in CalculiX *PLASTIC),
!   add the elastic part:  total strain = plastic strain + stress/E.
!   E = stress1/strain1 (first point = end of elasticity).
!   Beyond the last point the material is perfectly plastic, so keep a
!   positive hardening slope up to the largest strain you expect.
!
! STATEV (NSTATEV must be >= 14): use OUTPUT STATE VARIABLES = LOGICAL TRUE
! They import to Paraview as cell data
!   1-3  principal logarithmic strains of the total deformation (max,mid,min)
!   4-6  principal Cauchy stresses (max,mid,min)
!   7    von Mises (Cauchy) stress
!   8    equivalent plastic strain
!   9-14 elastic left Cauchy-Green tensor be: 11,22,33,12,13,23 (history)
!        (all zero at start = no plastic history, be is then set from F_n)
!
! STRESS returned is the Cauchy stress.  DDSDDE is the spatial elastic-plastic
! tangent (1/J times the log-strain consistent tangent); it is an approximation
! of the exact finite-strain tangent, so Newton convergence is not quadratic at
! large strain, but the residual (the stress) is exact.
! Component order: 11,22,33,12,[13,23]; NDI must be 3.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Eigen-decomposition of a symmetric 3x3 matrix (cyclic Jacobi).
! Columns of V are the eigenvectors, w the eigenvalues (unsorted).
!------------------------------------------------------------------------------
SUBROUTINE pf_jacobi3(A, w, V)
  USE Types
  IMPLICIT NONE
  REAL(KIND=dp), INTENT(IN)  :: A(3,3)
  REAL(KIND=dp), INTENT(OUT) :: w(3), V(3,3)
  REAL(KIND=dp) :: a_(3,3), theta, t, c, s, app, aqq, apq, arp, arq
  REAL(KIND=dp) :: off, scale, vp(3), vq(3)
  INTEGER :: sweep, p, q, r, m
  INTEGER, PARAMETER :: ip(3) = (/1,1,2/), iq(3) = (/2,3,3/), ir(3) = (/3,2,1/)

  a_ = A
  V = 0.0_dp
  V(1,1) = 1.0_dp; V(2,2) = 1.0_dp; V(3,3) = 1.0_dp
  scale = SUM(a_*a_)
  DO sweep = 1, 60
    off = a_(1,2)**2 + a_(1,3)**2 + a_(2,3)**2
    IF (off <= 1.0e-30_dp*scale) EXIT
    DO m = 1, 3
      p = ip(m); q = iq(m); r = ir(m)
      apq = a_(p,q)
      IF (apq == 0.0_dp) CYCLE
      app = a_(p,p); aqq = a_(q,q)
      theta = (aqq-app)/(2.0_dp*apq)
      t = SIGN(1.0_dp,theta)/(ABS(theta)+SQRT(theta*theta+1.0_dp))
      c = 1.0_dp/SQRT(t*t+1.0_dp)
      s = t*c
      a_(p,p) = app - t*apq
      a_(q,q) = aqq + t*apq
      a_(p,q) = 0.0_dp; a_(q,p) = 0.0_dp
      arp = a_(r,p); arq = a_(r,q)
      a_(r,p) = c*arp - s*arq; a_(p,r) = a_(r,p)
      a_(r,q) = c*arq + s*arp; a_(q,r) = a_(r,q)
      vp = V(:,p); vq = V(:,q)
      V(:,p) = c*vp - s*vq
      V(:,q) = s*vp + c*vq
    END DO
  END DO
  w(1) = a_(1,1); w(2) = a_(2,2); w(3) = a_(3,3)
END SUBROUTINE pf_jacobi3

!------------------------------------------------------------------------------
! Inverse and determinant of a 3x3 matrix
!------------------------------------------------------------------------------
SUBROUTINE pf_inv3(A, Ai, det)
  USE Types
  IMPLICIT NONE
  REAL(KIND=dp), INTENT(IN)  :: A(3,3)
  REAL(KIND=dp), INTENT(OUT) :: Ai(3,3), det
  Ai(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
  Ai(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
  Ai(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)
  Ai(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
  Ai(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
  Ai(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)
  Ai(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
  Ai(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
  Ai(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)
  det = A(1,1)*Ai(1,1) + A(1,2)*Ai(2,1) + A(1,3)*Ai(3,1)
  IF (det /= 0.0_dp) Ai = Ai/det
END SUBROUTINE pf_inv3

!------------------------------------------------------------------------------
! Yield stress Yk versus equivalent plastic strain pk and slope Hk per segment
! from the uniaxial (log total strain, true stress) points.
!------------------------------------------------------------------------------
SUBROUTINE pf_table_setup(n, eps, sig, E, pk, Yk, Hk, ierr)
  USE Types
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL(KIND=dp), INTENT(IN)  :: eps(n), sig(n)
  REAL(KIND=dp), INTENT(OUT) :: E, pk(n), Yk(n), Hk(n)
  INTEGER, INTENT(OUT) :: ierr
  INTEGER :: i

  ierr = 0
  IF (n < 1 .OR. eps(1) <= 0.0_dp .OR. sig(1) <= 0.0_dp) THEN
    ierr = 1; RETURN
  END IF
  E = sig(1)/eps(1)
  DO i=1,n
    pk(i) = eps(i) - sig(i)/E
    Yk(i) = sig(i)
  END DO
  pk(1) = 0.0_dp
  DO i=1,n-1
    IF (pk(i+1) <= pk(i)) THEN
      ierr = 2; RETURN
    END IF
    Hk(i) = (Yk(i+1)-Yk(i))/(pk(i+1)-pk(i))
    IF (Hk(i) < 0.0_dp) THEN
      ierr = 3; RETURN
    END IF
  END DO
  Hk(n) = 0.0_dp
END SUBROUTINE pf_table_setup

!------------------------------------------------------------------------------
! Radial return for the piecewise-linear hardening curve (exact, no iteration)
!------------------------------------------------------------------------------
SUBROUTINE pf_radial(npts, pk, Yk, Hk, G, qtr, pold, Yold, pnew, Ynew, Hnew, yielded)
  USE Types
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npts
  REAL(KIND=dp), INTENT(IN) :: pk(npts), Yk(npts), Hk(npts), G, qtr, pold
  REAL(KIND=dp), INTENT(OUT) :: Yold, pnew, Ynew, Hnew
  LOGICAL, INTENT(OUT) :: yielded
  REAL(KIND=dp) :: pcand
  INTEGER :: i, kseg

  kseg = 1
  DO i=2,npts
    IF (pk(i) <= pold) kseg = i
  END DO
  Yold = Yk(kseg) + Hk(kseg)*(pold-pk(kseg))
  pnew = pold; Ynew = Yold; Hnew = Hk(kseg)
  yielded = .FALSE.
  IF (qtr <= Yold*(1.0_dp+1.0e-10_dp)) RETURN
  yielded = .TRUE.
  pcand = pold
  DO i=kseg,npts
    pcand = (qtr + 3.0_dp*G*pold - Yk(i) + Hk(i)*pk(i))/(3.0_dp*G + Hk(i))
    kseg = i
    IF (i == npts) EXIT
    IF (pcand <= pk(i+1)) EXIT
  END DO
  pnew = MAX(pcand, pold)
  Hnew = Hk(kseg)
  Ynew = Yk(kseg) + Hnew*(pnew-pk(kseg))
END SUBROUTINE pf_radial

!------------------------------------------------------------------------------
! Common driver
!------------------------------------------------------------------------------
SUBROUTINE pf_drive(ntens, ndi, nstatev, nu, npts, eps, sig, dfrgrd0, dfrgrd1, &
     stress, statev, ddsdde, spd)
  USE Types
  USE Messages
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ntens, ndi, nstatev, npts
  REAL(KIND=dp), INTENT(IN) :: nu, eps(npts), sig(npts), dfrgrd0(3,3), dfrgrd1(3,3)
  REAL(KIND=dp), INTENT(OUT) :: stress(ntens), ddsdde(ntens,ntens)
  REAL(KIND=dp), INTENT(INOUT) :: statev(nstatev), spd

  REAL(KIND=dp) :: E, pk(npts), Yk(npts), Hk(npts), G, Kb
  REAL(KIND=dp) :: F0i(3,3), f(3,3), be(3,3), bt(3,3), tmp(3,3), V(3,3), w(3)
  REAL(KIND=dp) :: detF0, detF1, ea(3), tr, edev(3), en, qtr, pold, pnew
  REAL(KIND=dp) :: Yold, Ynew, Hnew, theta, thetab, coef, ean(3), tau(3)
  REAL(KIND=dp) :: sg(3,3), ben(3,3), nt(3,3), nvec(ntens), vv(3,3), scl
  REAL(KIND=dp) :: sv(6), hyd, dev(3), lam(3), b1(3,3), w1(3), V1(3,3)
  LOGICAL :: yielded
  INTEGER :: ierr, i, j, a

  IF (ndi /= 3 .OR. (ntens /= 4 .AND. ntens /= 6)) &
      CALL Fatal('UMATplastic_finite', 'Only NDI=3 (NTENS=4 or 6) is supported')
  IF (nstatev < 14) &
      CALL Fatal('UMATplastic_finite', 'At least 14 state variables are required')
  IF (nu <= -1.0_dp .OR. nu >= 0.5_dp) &
      CALL Fatal('UMATplastic_finite', 'Poisson ratio must satisfy -1 < nu < 0.5')
  CALL pf_table_setup(npts, eps, sig, E, pk, Yk, Hk, ierr)
  IF (ierr == 1) CALL Fatal('UMATplastic_finite', 'First strain/stress point must be positive')
  IF (ierr == 2) CALL Fatal('UMATplastic_finite', 'Strains must increase and every slope must be below E')
  IF (ierr == 3) CALL Fatal('UMATplastic_finite', 'Strain softening is not supported')

  G  = E/(2.0_dp*(1.0_dp+nu))
  Kb = E/(3.0_dp*(1.0_dp-2.0_dp*nu))

  CALL pf_inv3(dfrgrd0, F0i, detF0)
  detF1 = dfrgrd1(1,1)*(dfrgrd1(2,2)*dfrgrd1(3,3)-dfrgrd1(2,3)*dfrgrd1(3,2)) &
        - dfrgrd1(1,2)*(dfrgrd1(2,1)*dfrgrd1(3,3)-dfrgrd1(2,3)*dfrgrd1(3,1)) &
        + dfrgrd1(1,3)*(dfrgrd1(2,1)*dfrgrd1(3,2)-dfrgrd1(2,2)*dfrgrd1(3,1))
  IF (detF0 <= 0.0_dp .OR. detF1 <= 0.0_dp) &
      CALL Fatal('UMATplastic_finite', 'Non-positive determinant of the deformation gradient')

  ! Elastic left Cauchy-Green tensor at the previous converged level
  IF (statev(9) == 0.0_dp .AND. statev(10) == 0.0_dp .AND. statev(11) == 0.0_dp) THEN
    be = MATMUL(dfrgrd0, TRANSPOSE(dfrgrd0))
  ELSE
    be(1,1) = statev(9);  be(2,2) = statev(10); be(3,3) = statev(11)
    be(1,2) = statev(12); be(2,1) = statev(12)
    be(1,3) = statev(13); be(3,1) = statev(13)
    be(2,3) = statev(14); be(3,2) = statev(14)
  END IF

  ! Trial elastic state: be_tr = f be f^T,  f = F_{n+1} F_n^{-1}
  f = MATMUL(dfrgrd1, F0i)
  tmp = MATMUL(f, be)
  bt = MATMUL(tmp, TRANSPOSE(f))
  bt = 0.5_dp*(bt + TRANSPOSE(bt))
  CALL pf_jacobi3(bt, w, V)
  IF (MINVAL(w) <= 0.0_dp) &
      CALL Fatal('UMATplastic_finite', 'Elastic left Cauchy-Green tensor is not positive definite')

  ea  = 0.5_dp*LOG(w)                 ! principal elastic log strains
  tr  = SUM(ea)
  edev = ea - tr/3.0_dp
  en  = SQRT(SUM(edev**2))
  qtr = 2.0_dp*G*SQRT(1.5_dp)*en      ! von Mises (Kirchhoff) trial stress

  pold = statev(8)
  CALL pf_radial(npts, pk, Yk, Hk, G, qtr, pold, Yold, pnew, Ynew, Hnew, yielded)

  IF (yielded) THEN
    theta  = Ynew/qtr
    thetab = 1.0_dp/(1.0_dp + Hnew/(3.0_dp*G))
    coef   = 2.0_dp*G*(theta + thetab - 1.0_dp)
  ELSE
    theta = 1.0_dp; thetab = 1.0_dp; coef = 0.0_dp
  END IF
  ean = tr/3.0_dp + theta*edev
  tau = Kb*tr + 2.0_dp*G*theta*edev   ! principal Kirchhoff stresses

  sg = 0.0_dp; ben = 0.0_dp; nt = 0.0_dp
  DO a=1,3
    DO j=1,3
      DO i=1,3
        vv(i,j) = V(i,a)*V(j,a)
      END DO
    END DO
    sg  = sg  + tau(a)*vv
    ben = ben + EXP(2.0_dp*ean(a))*vv
    IF (en > 1.0e-30_dp) nt = nt + (edev(a)/en)*vv
  END DO
  scl = 1.0_dp/detF1
  sg = scl*sg

  ! Cauchy stress in Voigt form
  DO i=1,3
    stress(i) = sg(i,i)
  END DO
  stress(4) = sg(1,2)
  IF (ntens == 6) THEN
    stress(5) = sg(1,3); stress(6) = sg(2,3)
  END IF

  ! Tangent (spatial, 1/J times log-strain consistent tangent)
  DO i=1,3
    nvec(i) = nt(i,i)
  END DO
  nvec(4) = nt(1,2)
  IF (ntens == 6) THEN
    nvec(5) = nt(1,3); nvec(6) = nt(2,3)
  END IF
  ddsdde = 0.0_dp
  ddsdde(1:3,1:3) = Kb - 2.0_dp*G*theta/3.0_dp
  DO i=1,3
    ddsdde(i,i) = ddsdde(i,i) + 2.0_dp*G*theta
  END DO
  DO i=4,ntens
    ddsdde(i,i) = G*theta
  END DO
  DO j=1,ntens
    DO i=1,ntens
      ddsdde(i,j) = ddsdde(i,j) - coef*nvec(i)*nvec(j)
    END DO
  END DO
  ddsdde = scl*ddsdde

  ! State variables
  spd = spd + scl*0.5_dp*(Yold+Ynew)*(pnew-pold)
  statev(8)  = pnew
  statev(9)  = ben(1,1); statev(10) = ben(2,2); statev(11) = ben(3,3)
  statev(12) = ben(1,2); statev(13) = ben(1,3); statev(14) = ben(2,3)

  ! Diagnostics: principal log strains of F, principal Cauchy stresses, Mises
  b1 = MATMUL(dfrgrd1, TRANSPOSE(dfrgrd1))
  b1 = 0.5_dp*(b1 + TRANSPOSE(b1))
  CALL pf_jacobi3(b1, w1, V1)
  lam = 0.5_dp*LOG(w1)
  statev(1) = MAXVAL(lam); statev(3) = MINVAL(lam)
  statev(2) = SUM(lam) - statev(1) - statev(3)
  lam = scl*tau
  statev(4) = MAXVAL(lam); statev(6) = MINVAL(lam)
  statev(5) = SUM(lam) - statev(4) - statev(6)
  sv = 0.0_dp
  sv(1:3) = stress(1:3); sv(4:ntens) = stress(4:ntens)
  hyd = (sv(1)+sv(2)+sv(3))/3.0_dp
  dev = sv(1:3) - hyd
  statev(7) = SQRT(1.5_dp*(SUM(dev**2) + 2.0_dp*SUM(sv(4:6)**2)))
END SUBROUTINE pf_drive
!------------------------------------------------------------------------------
SUBROUTINE fplastic(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
     rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
     predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
     coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
     kstep, kinc)
!------------------------------------------------------------------------------
  USE Types
  USE Messages
  IMPLICIT NONE
  REAL(KIND=dp), INTENT(INOUT) :: STRESS(NTENS), STATEV(NSTATEV), SSE, SPD, SCD
  REAL(KIND=dp), INTENT(OUT) :: DDSDDE(NTENS,NTENS), rpl, ddsddt(NTENS), drplde(NTENS), drpldt
  REAL(KIND=dp), INTENT(IN) :: STRAN(NTENS), DSTRAN(NTENS), TIME(2), DTIME, TEMP, dtemp
  REAL(KIND=dp), INTENT(IN) :: predef(1), dpred(1), PROPS(NPROPS), coords(3), drot(3,3)
  REAL(KIND=dp), INTENT(INOUT) :: pnewdt
  REAL(KIND=dp), INTENT(IN) :: celent, DFRGRD0(3,3), DFRGRD1(3,3)
  CHARACTER(len=80), INTENT(IN) :: CMNAME
  INTEGER, INTENT(IN) :: NDI, NSHR, NTENS, NSTATEV, NPROPS, NOEL, NPT, layer, kspt, kstep, kinc

  REAL(KIND=dp) :: eps(100), sig(100)
  INTEGER :: i, npts

  rpl = 0.0_dp; ddsddt = 0.0_dp; drplde = 0.0_dp; drpldt = 0.0_dp
  npts = (nprops-1)/2
  IF (npts < 1 .OR. npts > 100) &
      CALL Fatal('fplastic', 'Need 1..100 strain/stress pairs after nu')
  DO i=1,npts
    eps(i) = PROPS(2*i)
    sig(i) = PROPS(2*i+1)
  END DO
  CALL pf_drive(ntens, ndi, nstatev, PROPS(1), npts, eps(1:npts), sig(1:npts), &
       DFRGRD0, DFRGRD1, stress, statev, ddsdde, spd)
END SUBROUTINE fplastic
