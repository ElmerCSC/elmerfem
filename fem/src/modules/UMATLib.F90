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
  SUBROUTINE bi_linear(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
       rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
       predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
       coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
       kstep, kinc)
! the call from Elmer
!CALL UMATusersubrtn(UMATSubrtn, StressVec(1:ntens), StateV, StressDer(1:ntens,1:ntens), EnergyElast, &
!          EnergyPlast, EnergyVisc, rpl, ddsddt(1:ntens), drplde(1:ntens), drpldt, &
!          stran(1:ntens), dstran(1:ntens), TimeAtStep, dtime, Temp, dTemp, &
!          predef, dpred, cmname, ndi, nshr, ntens, NStateV, InProps, NrInProps, coords, &
!          drot, pnewdt, celent, DefG0, DefG, ElementIndex, t, layer, kspt, kstep, kinc)
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
!  Material Constants 5
!  1) Nu, 2) strain at yield, 3) yield stress, 4)strain at final, 5) stress at final
!  State Variables 7
! 3 prinipal strains, 3 prinipal stress, mises
!------------------------------------------------------------------------------
    ! Local variables:
    REAL(KIND=dp) :: nu, E1, E2, LambdaLame, MuLame, G, G1, G2, E
    LOGICAL :: exists,f1exist,f2exist
    REAl(KIND=dp) :: mises, term1, term2, term3, term4, eprinmax, eprinmin
    REAL(KIND=dp) :: eprinmid,sprinmax,sprinmid,sprinmin, mises1
    real(kind=dp) :: strain1,stressc1,strain2,stressc2,term
    REAL(KIND=dp) :: totstran(ntens), astress(6), prevstran(ntens)
    Real*8 AA(3,3), ZZ(3,3)
    integer :: i
!------------------------------------------------------------------------------
  ! Get Young's modulus and the Poisson ratio:
    E1 = Props(3)/Props(2)
    E2 = (Props(5)-Props(3))/(Props(4)-Props(2))
    
    nu = Props(1)

    G1 = E1/(2*(1+nu))
    G2 = E2/(2*(1+nu))

    strain1 = Props(2)
    stressc1 = Props(3)
    strain2 = Props(4)
    stressc2 = Props(5)
    
!
! check mises1
E=E1
totstran = stran + dstran
prevstran = stran - dstran
!
if (ntens.eq.4) then
 mises1 = sqrt(stran(1)**2-stran(1)*stran(2)+stran(2)**2+3.*stran(4)**2)
! mises = sqrt(totstran(1)**2-totstran(1)*totstran(2)+totstran(2)**2+3.*totstran(4)**2)
end if
if (ntens.eq.6) then
 term1=(stran(1)-stran(2))**2
 term2=(stran(2)-stran(3))**2
 term3=(stran(3)-stran(1))**2
 term4=3.0*(stran(4)**2+stran(5)**2+stran(6)**2)
 mises1=sqrt( ((term1+term2+term3)/2)+term4 )
endif
!
!
if (mises1.ge.props(2)) then
!
  E = E2
!else
!   if (mises.gt.props(2)) then
!    write(20,fmt='(f8.6,a7,f8.6,a8,f8.6)') props(2)," stran ",mises1," totstra ",mises
!    E=E1*0.85
!   endif       
end if 
if (mises1.ge.props(6)) then
  E = 0.0
end if
!
!
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
   stress = stress + MATMUL(ddsdde,dstran)
   totstran = stran + dstran

! Calculate invariants for statev
!
Select case(ntens)
     case(4) 
     mises = sqrt(stress(1)**2-stress(1)*stress(2)+stress(2)**2+3.*stress(4)**2)
     term1=stress(1)-stress(2) 
     term2 = (0.5*term1)**2
     term3 = (stress(4))**2
     sprinmax = term1+sqrt(term2+term3)
     sprinmin = term1-sqrt(term2+term3)
     term1=totstran(1)-totstran(2) 
     term2 = (0.5*term1)**2
     term3 = (totstran(4))**2
     eprinmax = term1+sqrt(term2+term3)
     eprinmin = term1-sqrt(term2+term3)
     statev(1) = eprinmax
     statev(2) = 0.0
     statev(3) = eprinmin
     statev(4) = sprinmax
     statev(5) = 0.0
     statev(6) = sprinmin
     statev(7) = mises
     case(6)
     term1=(stress(1)-stress(2))**2
    term2=(stress(2)-stress(3))**2
    term3=(stress(3)-stress(1))**2
    term4=3.0*(stress(4)**2+stress(5)**2+stress(6)**2)
    mises=sqrt( ((term1+term2+term3)/2)+term4 )
    AA(1,1)=stress(1)
    AA(2,2)=stress(2)
    AA(3,3)=stress(3)
    AA(1,2)=stress(4)
    AA(1,3)=stress(6)
    AA(2,3)=stress(5) 
    AA(2,1)=AA(1,2)
    AA(3,1)=AA(1,3)
    AA(3,2)=AA(2,3)
     call jacobi2(AA,3,ZZ,1.E-6,100)
     sprinmax=AA(1,1)
     sprinmid=AA(2,2)
     sprinmin=AA(3,3)
    AA(1,1)=totstran(1)
    AA(2,2)=totstran(2)
    AA(3,3)=totstran(3)
    AA(1,2)=totstran(4)
    AA(1,3)=totstran(6)
    AA(2,3)=totstran(5) 
    AA(2,1)=AA(1,2)
    AA(3,1)=AA(1,3)
    AA(3,2)=AA(2,3)
     call jacobi2(AA,3,ZZ,1.E-6,100)
     eprinmax=AA(1,1)
     eprinmid=AA(2,2)
     eprinmin=AA(3,3)
     statev(1) = eprinmax
     statev(2) = eprinmid
     statev(3) = eprinmin
     statev(4) = sprinmax
     statev(5) = sprinmid
     statev(6) = sprinmin
     statev(7) = mises
 end select
!------------------------------------------------------------------------------
  END SUBROUTINE bi_linear
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE hsk_bi(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
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
    INTEGER :: PriLWork=102, PriInfo

    REAL(KIND=dp) :: SymBasis(6,3,3)
    REAL(KIND=dp) :: Identity(3,3), B(3,3), C(3,3), Strain(3,3), S(3,3), Sigma(3,3)
    REAL(KIND=dp) :: WorkMat(3,3)
    REAL(KIND=dp) :: StrainVec(ntens), Stress2(ntens)
    REAL(KIND=dp) :: DetDefG
    REAL(KIND=dp) :: nu, E, LambdaLame, MuLame, E1, E2, G, G1, G2
    REAL(KIND=dp) :: EigenVals(3), PriWork(102)
    LOGICAL :: exists,f1exist,f2exist
    REAl(KIND=dp) :: mises, term1, term2, term3, term4, eprinmax, eprinmin
    REAL(KIND=dp) :: eprinmid,sprinmax,sprinmid,sprinmin
    real(kind=dp) :: strain1,stressc1,strain2,stressc2,term
    REAL(KIND=dp) :: totstran(ntens)
    Real*8 AA(3,3), ZZ(3,3)
!------------------------------------------------------------------------------
! Get Young's modulus and the Poisson ratio:
    E1 = Props(3)/Props(2)
    E2 = (Props(5)-Props(3))/(Props(4)-Props(2))
    
    nu = Props(1)

    G1 = E1/(2*(1+nu))
    G2 = E2/(2*(1+nu))

    strain1 = Props(2)
    stressc1 = Props(3)
    strain2 = Props(4)
    stressc2 = Props(5)
    
!
! check mises
E=E1
!
if (ntens.eq.4) then
 mises = sqrt(stran(1)**2-stran(1)*stran(2)+stran(2)**2+3.*stran(4)**2)
end if
if (ntens.eq.6) then
 term1=(stran(1)-stran(2))**2
 term2=(stran(2)-stran(3))**2
 term3=(stran(3)-stran(1))**2
 term4=3.0*(stran(4)**2+stran(5)**2+stran(6)**2)
 mises=sqrt( ((term1+term2+term3)/2)+term4 )
endif
!
!
if (mises.gt.props(2)) then
!
     E = E2   
     
else
     
    E = E1
     
end if 
if (mises.gt.props(6)) then
  E = 0.0
end if
!if (time(1).gt.0.0 .and. noel.eq.550) then
!
!write(20,fmt='(a3,Es13.6,a5,F4.2,1x,es13.6,1x,es13.6)') "E= ",E," nu= ",nu,props(3),mises
!
!end if
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
     ! CALL Fatal( 'UMAT', 'DSYEV cannot generate eigen basis')          
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

!    DO i=1,ndi
!      Stress(i) = Sigma(i,i)
!    END DO
!    DO i=1,nshr
!      SELECT CASE(i)
!      CASE(1)
!        Stress(ndi+i) = Sigma(1,2)
!      CASE(2)
!        Stress(ndi+i) = Sigma(1,3)
!      CASE(3)
!        Stress(ndi+i) = Sigma(2,3)
!      END SELECT
!    END DO

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
!
    stress = stress + MATMUL(ddsdde,dstran)
! Calculate invariants for statev
!
Select case(ntens)
     case(4) 
     mises = sqrt(stress(1)**2-stress(1)*stress(2)+stress(2)**2+3.*stress(4)**2)
     term1=stress(1)-stress(2) 
     term2 = (0.5*term1)**2
     term3 = (stress(4))**2
     sprinmax = term1+sqrt(term2+term3)
     sprinmin = term1-sqrt(term2+term3)
     term1=totstran(1)-totstran(2) 
     term2 = (0.5*term1)**2
     term3 = (totstran(4))**2
     eprinmax = term1+sqrt(term2+term3)
     eprinmin = term1-sqrt(term2+term3)
     statev(1) = eprinmax
     statev(2) = 0.0
     statev(3) = eprinmin
     statev(4) = sprinmax
     statev(5) = 0.0
     statev(6) = sprinmin
     statev(7) = mises
     case(6)
      term1=(stress(1)-stress(2))**2
      term2=(stress(2)-stress(3))**2
      term3=(stress(3)-stress(1))**2
      term4=3.0*(stress(4)**2+stress(5)**2+stress(6)**2)
       mises=sqrt( ((term1+term2+term3)/2)+term4 )
    AA(1,1)=stress(1)
    AA(2,2)=stress(2)
    AA(3,3)=stress(3)
    AA(1,2)=stress(4)
    AA(1,3)=stress(6)
    AA(2,3)=stress(5) 
    AA(2,1)=AA(1,2)
    AA(3,1)=AA(1,3)
    AA(3,2)=AA(2,3)
     call jacobi2(AA,3,ZZ,1.E-6,100)
     sprinmax=AA(1,1)
     sprinmid=AA(2,2)
     sprinmin=AA(3,3)
    AA(1,1)=totstran(1)
    AA(2,2)=totstran(2)
    AA(3,3)=totstran(3)
    AA(1,2)=totstran(4)
    AA(1,3)=totstran(6)
    AA(2,3)=totstran(5) 
    AA(2,1)=AA(1,2)
    AA(3,1)=AA(1,3)
    AA(3,2)=AA(2,3)
     call jacobi2(AA,3,ZZ,1.E-6,100)
     eprinmax=AA(1,1)
     eprinmid=AA(2,2)
     eprinmin=AA(3,3)
     statev(1) = eprinmax
     statev(2) = eprinmid
     statev(3) = eprinmin
     statev(4) = sprinmax
     statev(5) = sprinmid
     statev(6) = sprinmin
     statev(7) = mises
 end select
!------------------------------------------------------------------------------
  END SUBROUTINE hsk_bi
!------------------------------------------------------------------------------
! 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE plastic(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
       rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
       predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
       coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
       kstep, kinc)
! the call from Elmer
!CALL UMATusersubrtn(UMATSubrtn, StressVec(1:ntens), StateV, StressDer(1:ntens,1:ntens), EnergyElast, &
!          EnergyPlast, EnergyVisc, rpl, ddsddt(1:ntens), drplde(1:ntens), drpldt, &
!          stran(1:ntens), dstran(1:ntens), TimeAtStep, dtime, Temp, dTemp, &
!          predef, dpred, cmname, ndi, nshr, ntens, NStateV, InProps, NrInProps, coords, &
!          drot, pnewdt, celent, DefG0, DefG, ElementIndex, t, layer, kspt, kstep, kinc)
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
!  Material Constants variable 
!  1) Nu, strain, stress pairs up to nprops (material constants)
!  State Variables 7
! 3 prinipal strains, 3 prinipal stress, mises
!------------------------------------------------------------------------------
    ! Local variables:
    REAL(KIND=dp) :: nu, Ev(100), LambdaLame, MuLame, E
    LOGICAL :: exists,f1exist,f2exist
    REAl(KIND=dp) :: mises, term1, term2, term3, term4, eprinmax, eprinmin
    REAL(KIND=dp) :: eprinmid,sprinmax,sprinmid,sprinmin, mises1
    REAL(KIND=dp) :: totstran(ntens), astress(6), prevstran(ntens)
    Real*8 AA(3,3), ZZ(3,3)
    integer :: i,j,k,l,m,n
!------------------------------------------------------------------------------
  ! Get Young's modulus and the Poisson ratio:
    j=(nprops-1)/2
    Ev(1) = Props(3)/Props(2)
    do i=2,j
     k=i+(i-2)
     l=i*2
     n=i+(i-2)+1
     m=(i*2)+1
    Ev(i) = (Props(m)-Props(n))/(Props(l)-Props(k))
    end do
    nu = Props(1)
    
!
! check mises1
E=Ev(1)
totstran = stran + dstran
prevstran = stran - dstran
!
if (ntens.eq.4) then
 mises1 = sqrt(stran(1)**2-stran(1)*stran(2)+stran(2)**2+3.*stran(4)**2)
! mises = sqrt(totstran(1)**2-totstran(1)*totstran(2)+totstran(2)**2+3.*totstran(4)**2)
end if
if (ntens.eq.6) then
 term1=(stran(1)-stran(2))**2
 term2=(stran(2)-stran(3))**2
 term3=(stran(3)-stran(1))**2
 term4=3.0*(stran(4)**2+stran(5)**2+stran(6)**2)
 mises1=sqrt( ((term1+term2+term3)/2)+term4 )
endif
!
k=2
do i=1,j-1
 !
 if ((0.90*mises1).ge.Props(i*2)) then
 E=Ev(i+1)
 k=i*2
 end if
end do
!  if (noel.eq.125 .and. npt.eq.1) then
!   write(20,*) 'E=',E,' i=',i," Strain=",stran(1)
!   write(20,*) 'Mises1:',mises1," Props:",Props(k),' k:',k
!  endif
!
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
   stress = stress + MATMUL(ddsdde,dstran)
   totstran = stran + dstran

! Calculate invariants for statev
!
Select case(ntens)
     case(4) 
     mises = sqrt(stress(1)**2-stress(1)*stress(2)+stress(2)**2+3.*stress(4)**2)
     term1=stress(1)-stress(2) 
     term2 = (0.5*term1)**2
     term3 = (stress(4))**2
     sprinmax = term1+sqrt(term2+term3)
     sprinmin = term1-sqrt(term2+term3)
     term1=totstran(1)-totstran(2) 
     term2 = (0.5*term1)**2
     term3 = (totstran(4))**2
     eprinmax = term1+sqrt(term2+term3)
     eprinmin = term1-sqrt(term2+term3)
     statev(1) = eprinmax
     statev(2) = 0.0
     statev(3) = eprinmin
     statev(4) = sprinmax
     statev(5) = 0.0
     statev(6) = sprinmin
     statev(7) = mises
     case(6)
     term1=(stress(1)-stress(2))**2
    term2=(stress(2)-stress(3))**2
    term3=(stress(3)-stress(1))**2
    term4=3.0*(stress(4)**2+stress(5)**2+stress(6)**2)
    mises=sqrt( ((term1+term2+term3)/2)+term4 )
    AA(1,1)=stress(1)
    AA(2,2)=stress(2)
    AA(3,3)=stress(3)
    AA(1,2)=stress(4)
    AA(1,3)=stress(6)
    AA(2,3)=stress(5) 
    AA(2,1)=AA(1,2)
    AA(3,1)=AA(1,3)
    AA(3,2)=AA(2,3)
     call jacobi2(AA,3,ZZ,1.E-6,100)
     sprinmax=AA(1,1)
     sprinmid=AA(2,2)
     sprinmin=AA(3,3)
    AA(1,1)=totstran(1)
    AA(2,2)=totstran(2)
    AA(3,3)=totstran(3)
    AA(1,2)=totstran(4)
    AA(1,3)=totstran(6)
    AA(2,3)=totstran(5) 
    AA(2,1)=AA(1,2)
    AA(3,1)=AA(1,3)
    AA(3,2)=AA(2,3)
     call jacobi2(AA,3,ZZ,1.E-6,100)
     eprinmax=AA(1,1)
     eprinmid=AA(2,2)
     eprinmin=AA(3,3)
     statev(1) = eprinmax
     statev(2) = eprinmid
     statev(3) = eprinmin
     statev(4) = sprinmax
     statev(5) = sprinmid
     statev(6) = sprinmin
     statev(7) = mises
 end select
!------------------------------------------------------------------------------
  END SUBROUTINE plastic
!------------------------------------------------------------------------------
!
! Subroutine needed for calculating stress invariants to store in state variables    
!
      SUBROUTINE JACOBI2 (AO,N,A,EPS,ITMAX)
      Real*8 :: AO(3,3), A(3,3),Z,Y,dif,tang,cose,sine,zz,yy
      Integer iter,i,j,itmax,N,nm1,ip1,irow,icol
      Real eps
      ITER=0
      DO I=1,N
       DO J=1,N
        A(I,J)=0.0
        A(I,I)=1.0
       end do
      end do
    
    do while (iter.lt.itmax)
      Z=0.
      NM1=N-1 

      DO I=1,NM1
       IP1=I+1
       DO J=IP1,N
        IF (ABS(AO(I,J)) .gt. Z) then
         Z=ABS(AO(I,J))
         IROW=I
         ICOL=J
        end if
       end do
      end do
   
   IF (ITER .EQ. 0) Y=Z*EPS
 
     IF (Z .gt. Y) then
       DIF=AO(IROW,IROW)-AO(ICOL,ICOL)
       TANG= (-DIF+SQRT(DIF**2+4.0*Z**2))/(2.0*AO(IROW,ICOL))
       COSE=1.0/SQRT(1.0+TANG**2)
       SINE=COSE*TANG
       DO I=1,N
        ZZ=A(I,IROW)
        A(I,IROW)=COSE*ZZ+SINE*A(I,ICOL)
        A(I,ICOL)=COSE*A(I,ICOL)-SINE*ZZ
       end do
       I=1
       do while (I .ne. IROW)  
        YY=AO(I,IROW)
        AO(I,IROW)=COSE*YY+SINE*AO(I,ICOL)
        AO(I,ICOL)=COSE*AO(I,ICOL)-SINE*YY
        I=I+1
       end do
      I=IROW+1
      do while (I .ne. ICOL) 
       YY=AO(IROW,I)
       AO(IROW,I)=COSE*YY+SINE*AO(I,ICOL)
       AO(I,ICOL)=COSE*AO(I,ICOL)-SINE*YY
       I=I+1
      end do
      I=ICOL+1
      do while (I .le. N)  
      ZZ=AO(IROW,I)
      AO(IROW,I)=COSE*ZZ+SINE*AO(ICOL,I)
      AO(ICOL,I)=COSE*AO(ICOL,I)-SINE*ZZ
      I=I+1
      end do
      YY=AO(IROW,IROW)
      AO(IROW,IROW)=YY*COSE**2+AO(IROW,ICOL)*2.0*COSE*SINE+&
AO(ICOL,ICOL)*SINE**2
      AO(ICOL,ICOL)=AO(ICOL,ICOL)*COSE**2+YY*SINE**2-AO(IROW,ICOL)&
*2.0*COSE*SINE
      AO(IROW,ICOL)=0.0
      
   else
   iter=itmax 
   end if
 
 iter=iter+1 
 end do 
 END subroutine jacobi2
!     SUBROUTINE FOR SOLVING AX=B USING GAUSS-SEIDEL METHOD
      SUBROUTINE JACOBI3(A,B,X)
      Real*8 :: A(6,6),B(6),X(6),XO(6),DIFF,SUM,xnew,check,xcheck
      integer i,n,j,nn,max,iter,k
      N=6
      MAX=1000
      EPSI= 0.000001
      DO 50 I=1,N
	  XO(I)=0.0
      X(I)=0.0
50    continue
      NN=0
51    DIFF=0.0
      DO 56 I=1,N
      SUM=0.0
      DO 53 J=1,N
      k = j - i
      IF(k .lt. 0) SUM=SUM+A(I,J)*XO(J)
      IF(k .gt. 0) SUM=SUM+A(I,J)*XO(J)
53    CONTINUE
      XNEW=(-SUM+B(I))/A(I,I)
      CHECK=ABS((XO(I)-XNEW)/XNEW)
      xcheck = check - diff
      IF(xcheck .le. 0) goto 55
54    DIFF=CHECK
55    X(I)=XNEW
56    CONTINUE
      do i=1,n
	  XO(I)=X(I)
	  enddo
      NN=NN+1
      k = nn - max
      IF(k .ge. 0) goto 58
      xcheck = diff - epsi
57    IF(xcheck .le. 0) goto 59
58    ITER=1
      MAX=NN
      RETURN
59    ITER=0
      MAX=NN
      RETURN
      END
! 
