!
! *  Authors: Mika Malinen
! *  Email:   mika.malinen@csc.fi
! *  Original Date: March 8, 2019
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
