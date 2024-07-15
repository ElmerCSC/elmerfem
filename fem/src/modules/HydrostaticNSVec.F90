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
! *  Module for the solution of Stokes equation with hydrostatic 1st order assumption.
! *  Utilizes multithreading and vectorization features initially introduced by Mikko Byckling.
! *  Inherits most of the code from IncompresibleNSVec module.
! *
! *  For the equations look at (5.70) and (5.71) in 
! *  Ralf Greve & Heinz Blatter: Dynamics of Ice Sheets and Glaciers, Springer, 2009. 
! *
! *  Authors: Peter Råback
! *  Email:   elmeradm@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Created: 30.08.2023
! *
!/*****************************************************************************/

MODULE HydrostaticNSUtils

  USE DefUtils

  INTEGER, POINTER, SAVE :: UpPointer(:), DownPointer(:)

CONTAINS

  ! Compute effective viscosity. This is needed not only by the bulkassembly but also 
  ! by the postprocessing when estimating pressure.
  !----------------------------------------------------------------------------------
  FUNCTION EffectiveViscosityVec( ngp, ntot, BasisVec, dBasisdxVec, Element, NodalVelo, &
      InitHandles, ViscNewton, ViscDerVec, DetJVec ) RESULT ( EffViscVec ) 

    IMPLICIT NONE 
    
    INTEGER :: ngp,ntot
    REAL(KIND=dp) :: BasisVec(:,:), dBasisdxVec(:,:,:)
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: NodalVelo(:,:)
    LOGICAL :: InitHandles
    LOGICAL, OPTIONAL :: ViscNewton
    REAL(KIND=dp), ALLOCATABLE, OPTIONAL :: ViscDerVec(:)
    REAL(KIND=dp), POINTER  :: EffViscVec(:)
    REAL(KIND=dp), OPTIONAL, ALLOCATABLE :: DetJVec(:)

    INTEGER :: allocstat,i,j,k,dim,dofs,n
    LOGICAL :: Found     
    CHARACTER(LEN=MAX_NAME_LEN) :: ViscModel
    TYPE(ValueHandle_t), SAVE :: Visc_h, ViscModel_h, ViscExp_h, ViscCritical_h, &
        ViscNominal_h, ViscDiff_h, ViscTrans_h, ViscYasuda_h, ViscGlenExp_h, ViscGlenFactor_h, &
        ViscArrSet_h, ViscArr_h, ViscTLimit_h, ViscRate1_h, ViscRate2_h, ViscEne1_h, ViscEne2_h, &
        ViscTemp_h
    REAL(KIND=dp) :: R, vgrad(2,3), NewtonRelax
    REAL(KIND=dp) :: c1, c2, c3, c4, Tlimit, ArrheniusFactor, A1, A2, Q1, Q2, ViscCond
    LOGICAL, SAVE :: ConstantVisc = .FALSE., Visited = .FALSE., DoNewton, GotRelax = .FALSE.
    REAL(KIND=dp), ALLOCATABLE, SAVE :: ss(:), s(:), ArrheniusFactorVec(:)
    REAL(KIND=dp), POINTER, SAVE :: ViscVec0(:), ViscVec(:), TempVec(:), EhfVec(:) 
    TYPE(Variable_t), POINTER, SAVE :: ShearVar, ViscVar, WeightVar
    LOGICAL, SAVE :: SaveShear, SaveVisc, SaveWeight
    CHARACTER(*), PARAMETER :: Caller = 'EffectiveViscosityVec'

    SAVE NewtonRelax, R
    
    !$OMP THREADPRIVATE(ss,s,ViscVec0,ViscVec,ArrheniusFactorVec)

    dim = 3
    dofs = 2
    DoNewton = .FALSE.
    IF(PRESENT(ViscNewton)) DoNewton = ViscNewton
    IF(DoNewton .AND. .NOT. PRESENT(ViscDerVec)) THEN
      CALL Fatal(Caller,'Newton linearization requires "ViscDerVec"')
    END IF
    n = Element % Type % NumberOfNodes
    
    IF(InitHandles ) THEN
      CALL Info(Caller,'Initializing handles for viscosity models',Level=8)

      CALL ListInitElementKeyword( Visc_h,'Material','Viscosity')      
      CALL ListInitElementKeyword( ViscModel_h,'Material','Viscosity Model')      

      IF( ListGetElementSomewhere( ViscModel_h) ) THEN
        ViscCond = ListGetCReal( CurrentModel % Solver % Values,&
            'Newtonian Viscosity Condition',Found )      
        ConstantVisc = ( Found .AND. ViscCond > 0.0_dp ) 

        IF( ListGetLogical( CurrentModel % Solver % Values,&
            'Constant-Viscosity Start', Found) ) ConstantVisc = (.NOT. Visited ) 

        CALL ListInitElementKeyword( ViscExp_h,'Material','Viscosity Exponent')      
        CALL ListInitElementKeyword( ViscCritical_h,'Material','Critical Shear Rate')      
        CALL ListInitElementKeyword( ViscNominal_h,'Material','Nominal Shear Rate')      
        CALL ListInitElementKeyword( ViscDiff_h,'Material','Viscosity Difference')      
        CALL ListInitElementKeyword( ViscTrans_h,'Material','Viscosity Transition')      
        CALL ListInitElementKeyword( ViscYasuda_h,'Material','Yasuda Exponent')      

        ! Do these initializations for glen's model only
        IF ( ListCompareElementAnyString( ViscModel_h,'glen') ) THEN
          CALL ListInitElementKeyword( ViscGlenExp_h,'Material','Glen Exponent',DefRValue=3.0_dp)
          CALL ListInitElementKeyword( ViscGlenFactor_h,'Material','Glen Enhancement Factor',DefRValue=1.0_dp)           
          CALL ListInitElementKeyword( ViscArrSet_h,'Material','Set Arrhenius Factor',DefLValue=.FALSE.)
          CALL ListInitElementKeyword( ViscArr_h,'Material','Arrhenius Factor')            
          CALL ListInitElementKeyword( ViscTLimit_h,'Material','Limit Temperature',DefRValue=-10.0_dp)
          CALL ListInitElementKeyword( ViscRate1_h,'Material','Rate Factor 1',DefRValue=3.985d-13)
          CALL ListInitElementKeyword( ViscRate2_h,'Material','Rate Factor 2',DefRValue=1.916d3)
          CALL ListInitElementKeyword( ViscEne1_h,'Material','Activation Energy 1',DefRValue=60.0d03)
          CALL ListInitElementKeyword( ViscEne2_h,'Material','Activation Energy 2',DefRValue=139.0d03)       
          CALL ListInitElementKeyword( ViscTemp_h,'Material','Relative Temperature')            

          IF (.NOT.ListCheckPresentAnyMaterial( CurrentModel,'Glen Allow Old Keywords')) THEN
            IF( ListCheckPresentAnyMaterial( CurrentModel,'Constant Temperature') ) THEN
              CALL Warn(Caller,'Replace >Constant Temperature< with >Relative Temperature<')
            END IF            
            IF( ListCheckPresentAnyMaterial( CurrentModel,'Temperature Field Variable') ) THEN             
              CALL Warn(Caller,'Replace >Temperature Field Variable< with >Relative Temperature = Equals ...<')
            END IF
          END IF
          IF( ViscTemp_h % NotPresentAnywhere .AND. ViscArr_h % NotPresentAnywhere) THEN
            CALL Fatal(Caller,'Neither >Relative Temperature< nor >Arrhenius Factor< given for viscosity model "glen"')
          END IF

          IF( ListCheckPresentAnyMaterial( CurrentModel,'Glen Enhancement Factor Function')  ) THEN
            CALL Fatal(Caller,'No Glen function API yet!')
          END IF
          R = GetConstReal( CurrentModel % Constants,'Gas Constant',Found)
          IF (.NOT.Found) R = 8.314_dp

          NewtonRelax = ListGetCReal( CurrentModel % Solver % Values,&
              'Viscosity Newton Relaxation Factor',GotRelax )
          IF(.NOT. GotRelax) NewtonRelax = 1.0_dp
        END IF
      ELSE
        CALL Info(Caller,'Using constant viscosity!')
      END IF
      
      SaveShear = .FALSE.
      SaveVisc = .FALSE.
      IF( PRESENT( DetJVec ) ) THEN
        ShearVar => VariableGet( CurrentModel % Mesh % Variables,'Shearrate',ThisOnly=.TRUE.)
        SaveShear = ASSOCIATED(ShearVar)
        ViscVar => VariableGet( CurrentModel % Mesh % Variables,'Viscosity',ThisOnly=.TRUE.)
        SaveVisc = ASSOCIATED(ViscVar)        
      END IF
      Found = .FALSE.
      IF(SaveShear) THEN
        IF(ShearVar % TYPE == Variable_on_gauss_points ) THEN
          CALL Info(Caller,'Saving "Shearrate" on ip points!',Level=10)
        ELSE IF( ShearVar % TYPE == Variable_on_elements ) THEN
          CALL Info(Caller,'Saving "Shearrate" on elements!',Level=10)
        ELSE IF( ShearVar % TYPE == Variable_on_nodes ) THEN
          CALL Info(Caller,'Saving "Shearrate" on nodes!',Level=10)
          Found = .TRUE.
        ELSE
          CALL Fatal(Caller,'Invalid field type for "Shearrate"!')
        END IF
        ShearVar % Values = 0.0_dp
      END IF

      IF(SaveVisc) THEN
        IF(ViscVar % TYPE == Variable_on_gauss_points ) THEN
          CALL Info(Caller,'Saving "Viscosity" on ip points!',Level=10)
        ELSE IF( ViscVar % TYPE == Variable_on_elements ) THEN
          CALL Info(Caller,'Saving "Viscosity" on elements!',Level=10)
        ELSE IF( ViscVar % TYPE == Variable_on_nodes ) THEN
          CALL Info(Caller,'Saving "Viscosity" on nodes!',Level=10)
          Found = .TRUE.
        ELSE
          CALL Fatal(Caller,'Invalid field type for "Shearrate"!')
        END IF
        ViscVar % Values = 0.0_dp
      END IF

      SaveWeight = .FALSE.
      IF( Found ) THEN
        WeightVar => VariableGet( CurrentModel % Mesh % Variables,'Viscosity Weight',ThisOnly=.TRUE.)
        IF( ASSOCIATED( WeightVar ) ) THEN
          IF( WeightVar % TYPE /= Variable_on_nodes ) THEN
            CALL Fatal(Caller,'Invalid field type for "Viscosity Weight"!')
          END IF
          SaveWeight = .TRUE.
          WeightVar % Values = 0.0_dp
        END IF
      END IF

      Visited = .TRUE.
    END IF

    ViscVec0 => ListGetElementRealVec( Visc_h, ngp, BasisVec, Element )

    ViscModel = ListGetElementString( ViscModel_h, Element, Found ) 
    IF( .NOT. Found ) THEN
      ! Return the plain viscosity
      EffViscVec => ViscVec0
      RETURN
    END IF

    ! Initialize derivative of viscosity for when newtonian linearization is used
    IF( DoNewton ) ViscDerVec(1:ngp) = 0.0_dp

    ! This reverts the viscosity model to linear 
    IF( ConstantVisc ) THEN
      EffViscVec => ViscVec0        
      RETURN      
    END IF

    ! Deallocate too small storage if needed 
    IF (ALLOCATED(ss)) THEN
      IF (SIZE(ss) < ngp ) DEALLOCATE(ss, s, ViscVec, ArrheniusFactorVec )
    END IF

    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(ss)) THEN
      ALLOCATE(ss(ngp),s(ngp),ViscVec(ngp),ArrheniusFactorVec(ngp),STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed')
      END IF
    END IF

    ! For non-newtonian models compute the viscosity here
    EffViscVec => ViscVec

    ! Calculate the strain rate velocity at all integration points
    DO k=1,ngp
      DO i=1,dofs
        DO j=1,dim
          vgrad(i,j) = SUM( dBasisdxVec(k,1:ntot,j) * NodalVelo(i,1:ntot) )
        END DO
      END DO
      ! Factor 4 is to have this compatible with IncompressibleNSVec!
      ss(k) = 4.0 *  (vgrad(1,1)**2 + vgrad(2,2)**2 + vgrad(1,1)*vgrad(2,2) + &
          0.5_dp * vgrad(1,2)*vgrad(2,1) + &
          0.25_dp * (vgrad(1,2)**2 + vgrad(1,3)**2 + vgrad(2,1)**2 + vgrad(2,3)**2) )
    END DO
    !      ss(1:ngp) = 0.5_dp * ss(1:ngp)

    IF(SaveShear) THEN
      IF( ShearVar % TYPE == Variable_on_nodes ) THEN
        DO i=1,n
          ShearVar % Values(ShearVar % Perm(Element % NodeIndexes(i))) = &
              ShearVar % Values(ShearVar % Perm(Element % NodeIndexes(i))) + &
              SUM(BasisVec(1:ngp,i)*ss(1:ngp)*DetJVec(1:ngp))
        END DO
      ELSE
        i = Element % ElementIndex
        IF( ShearVar % TYPE == Variable_on_gauss_points ) THEN
          j = ShearVar % Perm(i+1) - ShearVar % Perm(i)
          IF(j /= ngp) THEN
            CALL Fatal(Caller,'Expected '//I2S(j)//' gauss point for "Shearrate" got '//I2S(ngp))
          END IF
          ShearVar % Values(ShearVar % Perm(i)+1:ShearVar % Perm(i+1)) = ss(1:ngp)
        ELSE IF( ShearVar % TYPE == Variable_on_elements ) THEN
          ShearVar % Values(ShearVar % Perm(i)) = SUM(ss(1:ngp)) / ngp
        END IF
      END IF
    END IF


    SELECT CASE( ViscModel )       

    CASE('glen')
      c2 = ListGetElementReal( ViscGlenExp_h,Element=Element,Found=Found)

      ! the second invariant is not taken from the strain rate tensor,
      ! but rather 2*strain rate tensor (that's why we divide by 4 = 2**2)        
      s(1:ngp) = ss(1:ngp)/4.0_dp

      c3 = ListGetElementReal( ViscCritical_h,Element=Element,Found=Found)
      IF( Found ) THEN
        c3 = c3**2
        WHERE( s(1:ngp) < c3 ) s(1:ngp) = c3
      END IF

      IF( ListGetElementLogical( ViscArrSet_h,Element,Found=Found) ) THEN
        ArrheniusFactor = ListGetElementReal( ViscArr_h,Element=Element)
        ViscVec(1:ngp) = 0.5_dp * (ArrheniusFactor)**(-1.0_dp/c2) * s(1:ngp)**(((1.0_dp/c2)-1.0_dp)/2.0_dp);                    

        IF( DoNewton ) THEN
          WHERE( s(1:ngp) > c3 ) ViscDerVec(1:ngp) = 0.5_dp * ArrheniusFactor**(-1.0_dp/c2) &
              * ((1.0_dp/c2)-1.0_dp)/2.0_dp * s(1:ngp)**(((1.0_dp/c2)-1.0_dp)/2.0_dp - 1.0_dp)/4.0_dp
        END IF
      ELSE         
        ! lets for the time being have this hardcoded
        Tlimit = ListGetElementReal( ViscTlimit_h,Element=Element)
        A1 = ListGetElementReal( ViscRate1_h,Element=Element)
        A2 = ListGetElementReal( ViscRate2_h,Element=Element)
        Q1 = ListGetElementReal( ViscEne1_h,Element=Element)
        Q2 = ListGetElementReal( ViscEne2_h,Element=Element)

        ! WHERE is faster than DO + IF
        TempVec => ListGetElementRealVec( ViscTemp_h, ngp, BasisVec, Element )

        WHERE( TempVec(1:ngp ) < Tlimit )
          ArrheniusFactorVec(1:ngp) = A1 * EXP( -Q1/(R * (273.15_dp + TempVec(1:ngp))))
        ELSE WHERE( TempVec(1:ngp) > 0.0_dp ) 
          ArrheniusFactorVec(1:ngp) = A2 * EXP( -Q2/(R * (273.15_dp)))
        ELSE WHERE
          ArrheniusFactorVec(1:ngp) = A2 * EXP( -Q2/(R * (273.15_dp + TempVec(1:ngp))))
        END WHERE

        EhfVec => ListGetElementRealVec( ViscGlenFactor_h, ngp, BasisVec,Element=Element )
        ViscVec(1:ngp) = 0.5_dp * (EhFVec(1:ngp) * ArrheniusFactorVec(1:ngp))**(-1.0_dp/c2) * &
            s(1:ngp)**(((1.0_dp/c2)-1.0_dp)/2.0_dp);

        IF( DoNewton ) THEN
          WHERE( s(1:ngp) > c3 ) 
            ViscDerVec(1:ngp) = 0.5_dp * (  EhFVec(1:ngp) * ArrheniusFactorVec(1:ngp))**(-1.0_dp/c2) &
                * ((1.0_dp/c2)-1.0_dp)/2.0_dp * s(1:ngp)**(((1.0_dp/c2)-1.0_dp)/2.0_dp - 1.0_dp)/4.0_dp
          END WHERE
        END IF
      END IF

    CASE('power law')
      c2 = ListGetElementReal( ViscExp_h,Element=Element)

      c3 = ListGetElementReal( ViscCritical_h,Element=Element,Found=Found)       
      IF( Found ) THEN
        c3 = c3**2
        WHERE( ss(1:ngp) < c3 ) ss(1:ngp) = c3
      END IF

      ViscVec(1:ngp) = ViscVec0(1:ngp) * ss(1:ngp)**((c2-1)/2)

      IF (DoNewton ) THEN
        WHERE(ss(1:ngp) /= 0) ViscDerVec(1:ngp) = &
            ViscVec0(1:ngp) * (c2-1)/2 * ss(1:ngp)**((c2-1)/2-1)
      END IF

      c4 = ListGetElementReal( ViscNominal_h,Element=Element,Found=Found)
      IF( Found ) THEN
        ViscVec(1:ngp) = ViscVec(1:ngp) / c4**(c2-1)
        IF (DoNewton ) THEN
          ViscDerVec(1:ngp) = ViscDerVec(1:ngp) / c4**(c2-1)
        END IF
      END IF

    CASE('power law too')
      c2 = ListGetElementReal( ViscExp_h,Element=Element)           
      ViscVec(1:ngp) = ViscVec0(1:ngp)**(-1/c2)* ss(1:ngp)**(-(c2-1)/(2*c2)) / 2

      IF (DoNewton ) THEN
        ViscDerVec(1:ngp) = ViscVec0(1:ngp)**(-1/c2)*(-(c2-1)/(2*c2))*ss(1:ngp)*(-(c2-1)/(2*c2)-1) / 2
      END IF

    CASE ('carreau')      
      c1 = ListGetElementReal( ViscDiff_h,Element=Element)
      c2 = ListGetElementReal( ViscExp_h,Element=Element)
      c3 = ListGetElementReal( ViscTrans_h,Element=Element)
      c4 = ListGetElementReal( ViscYasuda_h,Element=Element,Found=Found)
      IF( Found ) THEN
        ViscVec(1:ngp) = ViscVec0(1:ngp) + c1 * (1 + c3**c4*ss(1:ngp)**(c4/2))**((c2-1)/c4) 

        IF( DoNewton ) THEN
          ViscDerVec(1:ngp) = c1*(1+c3**c4*ss(1:ngp)**(c4/2))**((c2-1)/c4-1)*(c2-1)/2*c3**c4*&
              ss(1:ngp)**(c4/2-1)
        END IF
      ELSE
        ViscVec(1:ngp) = ViscVec0(1:ngp) + c1 * (1 + c3*c3*ss(1:ngp))**((c2-1)/2) 

        IF( DoNewton ) THEN
          ViscDerVec(1:ngp) = c1*(c2-1)/2*c3**2*(1+c3**2*ss(1:ngp))**((c2-1)/2-1)
        END IF
      END IF

    CASE ('cross')
      c1 = ListGetElementReal( ViscDiff_h,Element=Element)
      c2 = ListGetElementReal( ViscExp_h,Element=Element)
      c3 = ListGetElementReal( ViscTrans_h,Element=Element)

      ViscVec(1:ngp) = ViscVec0(1:ngp) + c1 / (1 + c3*ss(1:ngp)**(c2/2))

      IF( DoNewton ) THEN
        ViscDerVec(1:ngp) = -c1*c3*ss(1:ngp)**(c2/2)*c2 / (2*(1+c3*ss(1:ngp)**(c2/2))**2*ss(1:ngp))
      END IF

    CASE ('powell eyring')
      c1 = ListGetElementReal( ViscDiff_h,Element=Element)
      c2 = ListGetElementReal( ViscTrans_h,Element=Element)

      s(1:ngp) = SQRT(ss(1:ngp))

      IF( DoNewton ) THEN          
        WHERE( c2*s(1:ngp) < 1.0d-5 )
          ViscVec(1:ngp) = ViscVec0(1:ngp) + c1
          ViscDerVec(1:ngp) = 0.0_dp
        ELSE WHERE
          ViscVec(1:ngp) = ViscVec0(1:ngp) + c1 * LOG(c2*s(1:ngp)+SQRT(c2*c2*ss(1:ngp)+1))/(c2*ss(1:ngp))            
          ViscDerVec(1:ngp) = c1*(c2/(2*s(1:ngp))+c2**2/(2*SQRT(c2**2*ss(1:ngp)+1)))/ &
              ((c2*s(1:ngp)+SQRT(c2*ss(1:ngp)+1))*c2*s(1:ngp)) - &
              c1*LOG(c2*s(1:ngp)+SQRT(c2**2*ss(1:ngp)+1))/(c2*s(1:ngp)**3)/2
        END WHERE
      ELSE
        WHERE( c2*s(1:ngp) < 1.0d-5 )
          ViscVec(1:ngp) = ViscVec0(1:ngp) + c1
        ELSE WHERE
          ViscVec(1:ngp) = ViscVec0(1:ngp) + c1 * LOG(c2*s(1:ngp)+SQRT(c2*c2*ss(1:ngp)+1))/(c2*ss(1:ngp))            
        END WHERE
      END IF

    CASE DEFAULT 
      CALL Fatal(Caller,'Unknown material model')

    END SELECT

    IF( DoNewton ) THEN
      IF(GotRelax) ViscDerVec(1:ngp) = NewtonRelax * ViscDerVec(1:ngp)
    END IF

    ! If requested, save viscosity field (on nodes, ip points or elements). 
    IF(SaveVisc) THEN
      IF( ViscVar % TYPE == Variable_on_nodes ) THEN
        DO i=1,n
          ViscVar % Values(ViscVar % Perm(Element % NodeIndexes(i))) = &
              ViscVar % Values(ViscVar % Perm(Element % NodeIndexes(i))) + &
              SUM(BasisVec(1:ngp,i)*ViscVec(1:ngp)*detJVec(1:ngp))
        END DO
      ELSE
        i = Element % ElementIndex
        IF( ViscVar % TYPE == Variable_on_gauss_points ) THEN
          j = ViscVar % Perm(i+1) - ViscVar % Perm(i) 
          IF(j /= ngp) THEN
            CALL Fatal(Caller,'Expected '//I2S(j)//' gauss point for "Viscosity" got '//I2S(ngp))
          END IF
          ViscVar % Values(ViscVar % Perm(i)+1:ViscVar % Perm(i+1)) = ViscVec(1:ngp)
        ELSE
          ViscVar % Values(ViscVar % Perm(i)) = SUM(ViscVec(1:ngp)) / ngp
        END IF
      END IF
    END IF

    ! If requested, save normalization weight associated to viscosity (and shearrate).
    IF(SaveWeight) THEN
      DO i=1,n
        WeightVar % Values(WeightVar % Perm(Element % NodeIndexes(i))) = &
            WeightVar % Values(WeightVar % Perm(Element % NodeIndexes(i))) + &
            SUM(BasisVec(1:ngp,i)*detJVec(1:ngp))
      END DO
    END IF

  END FUNCTION EffectiveViscosityVec

  
!------------------------------------------------------------------------------
! Assemble local finite element matrix for a single bulk element and glue
! it to the global matrix. Routine is vectorized and multithreaded.
!------------------------------------------------------------------------------
  SUBROUTINE LocalBulkMatrix(Element, n, nd, ntot, &
       SpecificLoad, LinearAssembly, nb, Newton, InitHandles )
!------------------------------------------------------------------------------
    USE LinearForms
    IMPLICIT NONE

    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    INTEGER, INTENT(IN) :: n, nd, ntot, nb
    LOGICAL, INTENT(IN) :: SpecificLoad
    LOGICAL, INTENT(IN) :: LinearAssembly, Newton, InitHandles
!------------------------------------------------------------------------------
    INTEGER :: dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    LOGICAL :: Stat, Found, IsPelem

    REAL(KIND=dp), TARGET :: STIFF(ntot*2,ntot*2), FORCE(ntot*2)
    REAL(KIND=dp) :: NodalVelo(2,ntot),NodalHeight(ntot)
    REAL(KIND=dp) :: s, rho, crho
    REAL(KIND=dp), ALLOCATABLE, SAVE :: BasisVec(:,:), dBasisdxVec(:,:,:), DetJVec(:), &
        rhoVec(:), loadAtIpVec(:,:), ForcePart(:), GradVec(:,:,:), GradHeight(:,:), &
        weight_1(:), weight_2(:), weight_4(:), tauVec(:)        
    REAL(KIND=dp), POINTER :: muVec(:), LoadVec(:)
    REAL(KIND=dp), ALLOCATABLE :: muDerVec0(:),g(:,:,:),StrainRateVec(:,:,:)
    REAL(kind=dp) :: stifford(ntot,ntot,2,2), jacord(ntot,ntot,2,2), &
        JAC(ntot*2,ntot*2 ), SOL(ntot*2)

    INTEGER :: t, i, j, k, p, q, ngp, allocstat, dofs
    INTEGER, SAVE :: elemdim
    CHARACTER(LEN=MAX_NAME_LEN):: str

    TYPE(ValueHandle_t), SAVE :: Dens_h, Load_h(3)
    TYPE(Variable_t), POINTER, SAVE :: HeightVar 

!DIR$ ATTRIBUTES ALIGN:64 :: BasisVec, dBasisdxVec, DetJVec, rhoVec, loadAtIpVec
!DIR$ ATTRIBUTES ALIGN:64 :: STIFF, FORCE, weight_1, weight_2, weight_4
!$OMP THREADPRIVATE(BasisVec, dBasisdxVec, DetJVec, rhoVec, loadAtIpVec, ElemDim )
!$OMP THREADPRIVATE(ForcePart, weight_1, weight_2, weight_4)
!$OMP THREADPRIVATE(tauVec, GradVec, Nodes)

    SAVE Nodes
!------------------------------------------------------------------------------

    CALL GetElementNodesVec( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp
    JAC   = 0._dp
    JacOrd = 0._dp
    stifford = 0._dp

    dofs = 2
    dim = 3

    ! Numerical integration:
    !-----------------------
    IsPelem = isPElement(Element)

    IP = GaussPointsAdapt(Element, PReferenceElement = isPelem )
    ngp = IP % n

    ElemDim = Element % Type % Dimension
    
    ! Storage size depending ngp
    !-------------------------------------------------------------------------------
    
    ! Deallocate storage if needed 
    IF (ALLOCATED(BasisVec)) THEN
      IF (SIZE(BasisVec,1) < ngp .OR. SIZE(BasisVec,2) < ntot) &
          DEALLOCATE(BasisVec,dBasisdxVec, DetJVec, rhoVec, &
          LoadAtIpVec, weight_1, weight_2, weight_4, tauVec, &
          ForcePart, GradVec, GradHeight)
    END IF
    
    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(BasisVec)) THEN
      ALLOCATE(BasisVec(ngp,ntot), dBasisdxVec(ngp,ntot,3), DetJVec(ngp), &
          rhoVec(ngp), LoadAtIpVec(ngp,2), &
          weight_1(ngp), weight_2(ngp), weight_4(ngp), tauVec(ngp), &
          ForcePart(ntot), GradVec(ngp,3,3), GradHeight(ngp,dim), &
          STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal('HydrostaticNSSolver::LocalBulkMatrix','Local storage allocation failed')
      END IF
    END IF

    IF (Newton) THEN
      ALLOCATE(muDerVec0(ngp), g(ngp,ntot,dim), StrainRateVec(ngp,dim,dim))
      muDerVec0 = 0._dp
    END IF
    
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Dens_h,'Material','Density')      
      CALL ListInitElementKeyword( Load_h(1),'Body Force','Flow Bodyforce 1')
      CALL ListInitElementKeyword( Load_h(2),'Body Force','Flow Bodyforce 2')
      CALL ListInitElementKeyword( Load_h(3),'Body Force','Flow Bodyforce 3')

      str = ListGetString( CurrentModel % Solver % Values,'Height Variable Name',Found )
      IF(.NOT. Found) str = 'height'            
      HeightVar => VariableGet( CurrentModel % Mesh % Variables, str) 
    END IF
    
    ! Vectorized basis functions
    stat = ElementInfoVec(Element, nodes, ngp, IP % U, IP % V, &
        IP % W, detJvec, SIZE(basisVec, 2), BasisVec, dBasisdxVec)

    ! Weights at integration points
    DO t = 1, ngp
      DetJVec(t) = DetJVec(t) * IP % s(t)
    END DO

    ! We assume constant density so far:
    !-----------------------------------
    rho = ListGetElementReal( Dens_h, Element = Element )     
    
    ! Velocity from previous iteration at integration points
    IF(.NOT. LinearAssembly) THEN
      CALL GetLocalSolution( NodalVelo )
    END IF
    
    CALL GetLocalSolution( NodalHeight, UVariable = HeightVar ) 

    ! Return the effective viscosity. Currently only non-newtonian models supported.
    muvec => EffectiveViscosityVec( ngp, ntot, BasisVec, dBasisdxVec, Element, NodalVelo, &
        InitHandles, Newton, muDerVec0, DetJVec )
    
    ! Rho 
    rhovec(1:ngp) = rho
    IF( SpecificLoad ) THEN
      crho = 1.0_dp
    ELSE
      crho = rho
    END IF

    ! Flow bodyforce if present
    LoadAtIpVec = 0._dp
    DO i=1,dim
      LoadVec => ListGetElementRealVec( Load_h(i), ngp, BasisVec, Element, Found ) 
      IF( .NOT. Found ) CYCLE
      IF(i<3) THEN
        LoadAtIpVec(1:ngp,i) = crho * LoadVec(1:ngp)
      ELSE
        ! Project the 3rd force component to (x,y) plane.
        DO k=1,2
          DO j=1,ngp
            LoadAtIpVec(j,k) = LoadAtIpVec(j,k) + crho * SUM(dBasisdxVec(j,1:ntot,k)*NodalHeight(1:ntot)) * LoadVec(j)
          END DO
        END DO
      END IF
    END DO

    IF ( Newton ) THEN
      IF (ANY(muDerVec0(1:ngp)/=0)) THEN
        ! Note that it is quite delicate whether to use dim or dofs in the looping.
        ! We need to create du_x/dz and du_y/dz but don't have u_z at our disposal.
        DO i = 1,dim
          DO j = 1,dim
            IF(i<=dofs) THEN
              GradVec(1:ngp, i, j) = MATMUL(dBasisdxVec(1:ngp,1:ntot,j),NodalVelo(i,1:ntot))
            ELSE
              GradVec(1:ngp, i, j) = 0.0_dp
            END IF
          END DO
        END DO
        
        DO i = 1,dim
          DO j = 1,dim
            StrainRateVec(1:ngp,i,j) = ( GradVec(1:ngp,i,j) + GradVec(1:ngp,j,i) ) / 2
          END DO
        END DO
        
        DO i=1,dofs
          DO q = 1,ntot
            g(1:ngp,q,i) = SUM(StrainRateVec(1:ngp,i,1:dim)*dBasisdxvec(1:ngp,q,1:dim),2)
          END DO
        END DO
        
        muDerVec0(1:ngp) = muderVec0(1:ngp)*detJVec(1:ngp)*8
        DO i=1,dofs
          DO j=1,dofs
            CALL LinearForms_udotv(ngp,ntot,dim,g(:,:,j),g(:,:,i),mudervec0,jacord(:,:,j,i))
          END DO
        END DO
      END IF
    END IF
    
    weight_1(1:ngp) = muVec(1:ngp) * detJVec(1:ngp)
    weight_2(1:ngp) = 2*weight_1(1:ngp)
    weight_4(1:ngp) = 4*weight_1(1:ngp)

    i = 1; j = 1 
    CALL LinearForms_UdotV(ngp, ntot, elemdim, &
        dbasisdxvec(:, :, 1), dbasisdxvec(:,:,1), weight_4, stifford(:,:,i,j))
    CALL LinearForms_UdotV(ngp, ntot, elemdim, &
        dbasisdxvec(:, :, 2), dbasisdxvec(:,:,2), weight_1, stifford(:,:,i,j))
    CALL LinearForms_UdotV(ngp, ntot, elemdim, &
        dbasisdxvec(:, :, 3), dbasisdxvec(:,:,3), weight_1, stifford(:,:,i,j))
    
    i = 1; j = 2
    CALL LinearForms_UdotV(ngp, ntot, elemdim, &
        dbasisdxvec(:, :, i), dbasisdxvec(:,:,j), weight_2, stifford(:,:,i,j))
    CALL LinearForms_UdotV(ngp, ntot, elemdim, &
        dbasisdxvec(:, :, j), dbasisdxvec(:,:,i), weight_1, stifford(:,:,i,j))

    i = 2; j = 1
    CALL LinearForms_UdotV(ngp, ntot, elemdim, &
        dbasisdxvec(:, :, i), dbasisdxvec(:,:,j), weight_2, stifford(:,:,i,j))
    CALL LinearForms_UdotV(ngp, ntot, elemdim, &
        dbasisdxvec(:, :, j), dbasisdxvec(:,:,i), weight_1, stifford(:,:,i,j))
    
    i = 2; j = 2 
    CALL LinearForms_UdotV(ngp, ntot, elemdim, &
        dbasisdxvec(:, :, 1), dbasisdxvec(:,:,1), weight_1, stifford(:,:,i,j))
    CALL LinearForms_UdotV(ngp, ntot, elemdim, &
        dbasisdxvec(:, :, 2), dbasisdxvec(:,:,2), weight_4, stifford(:,:,i,j))
    CALL LinearForms_UdotV(ngp, ntot, elemdim, &
        dbasisdxvec(:, :, 3), dbasisdxvec(:,:, 3), weight_1, stifford(:,:,i,j))
    
    ! add loads
    DO i = 1,2
      ForcePart = 0._dp
      CALL LinearForms_UdotF(ngp, ntot, basisVec, detJVec, LoadAtIpVec(:,i), ForcePart)
      FORCE(i::dofs) = ForcePart(1:ntot)
    END DO

    DO i = 1, DOFS
      DO j = 1, DOFS
        Stiff(i::DOFS, j::DOFS) = StiffOrd(1:ntot, 1:ntot, i,j)
      END DO
    END DO

    IF ( Newton ) THEN
      SOL=0._dp
      DO i = 1, DOFS
        DO j = 1, DOFS
          JAC(i::DOFS, j::DOFS) = JacOrd(1:ntot, 1:ntot,i,j)
        END DO
        SOL(i::DOFs) = NodalVelo(i,1:ntot)
      END DO
      STIFF = STIFF + JAC
      FORCE = FORCE + MATMUL(JAC,SOL)
    END IF

    IF ( nb>0 ) THEN
      CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    END IF
    
    CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element, VecAssembly=.TRUE.)

  END SUBROUTINE LocalBulkMatrix
!------------------------------------------------------------------------------


  !------------------------------------------------------------------------------
  ! Assemble local finite element matrix for a single boundary element and glue
  ! it to the global matrix.
  !------------------------------------------------------------------------------
  SUBROUTINE LocalBoundaryMatrix( Element, n, nd, dim, InitHandles, FrictionNewton)
    !------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    INTEGER, INTENT(IN) :: n, nd, dim
    LOGICAL, INTENT(INOUT) :: InitHandles 
    LOGICAL :: FrictionNewton
    !------------------------------------------------------------------------------    
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp), TARGET :: STIFF(nd*2,nd*2), FORCE(nd*2),NodalHeight(nd)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    INTEGER :: c,i,j,k,l,p,q,t,ngp,dofs
    LOGICAL :: HaveSlip, HaveForce, HavePres, HaveFrictionW, HaveFrictionU, &
        HaveFriction, Found, Stat, GotRelax
    REAL(KIND=dp) :: ExtPressure, s, detJ, wut0, wexp, wcoeff, ut
    REAL(KIND=dp) :: SlipCoeff(3), SurfaceTraction(3), Normal(3), &
        Velo(3), TanFrictionCoeff, DummyVals(1)
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t), SAVE :: ExtPressure_h, IntPressure_h, SurfaceTraction_h, SlipCoeff_h, &
        WeertmanCoeff_h, WeertmanExp_h, FrictionUt0_h, FrictionCoeff_h
    TYPE(VariableHandle_t), SAVE :: Velo_v
    TYPE(Variable_t), POINTER, SAVE :: NrmSol, HeightVar
    TYPE(ValueList_t), POINTER :: BC    
    REAL(KIND=dp) :: TanFder,JAC(nd*2,nd*2),SOL(nd*2),NodalSol(2,nd),NewtonRelax
    TYPE(Variable_t), POINTER, SAVE :: SlipCoeffVar, SlipSpeedVar, SlipWeightVar
    LOGICAL, SAVE :: SaveSlipSpeed, SaveSlipCoeff, SaveSlipWeight

    
    SAVE Basis, dBasisdx, NewtonRelax, GotRelax
    
    !------------------------------------------------------------------------------

    IF( InitHandles ) THEN   
      CALL ListInitElementKeyword( ExtPressure_h,'Boundary Condition','Normal Surface Traction')
      IF( .NOT. ListGetElementSomewhere( ExtPressure_h) ) THEN
        CALL ListInitElementKeyword( ExtPressure_h,'Boundary Condition','External Pressure')      
      END IF
      CALL ListInitElementKeyword( IntPressure_h,'Boundary Condition','Internal Pressure')      
      CALL ListInitElementKeyword( SurfaceTraction_h,'Boundary Condition','Surface Traction',InitVec3D=.TRUE.)
      CALL ListInitElementKeyword( SlipCoeff_h,'Boundary Condition','Slip Coefficient',InitVec3D=.TRUE.)

      CALL ListInitElementKeyword( FrictionCoeff_h,'Boundary Condition','Friction Coefficient',&
          EvaluateAtIp=.TRUE., DummyCount=1)     
      CALL ListInitElementKeyword( FrictionUt0_h,'Boundary Condition','Friction Linear Velocity')
      IF(FrictionUt0_h % NotPresentAnywhere) THEN
        CALL ListInitElementKeyword( FrictionUt0_h,'Boundary Condition','Weertman Linear Velocity')
      END IF  
    
      CALL ListInitElementKeyword( WeertmanCoeff_h,'Boundary Condition','Weertman Friction Coefficient')
      CALL ListInitElementKeyword( WeertmanExp_h,'Boundary Condition','Weertman Exponent')

      NewtonRelax = ListGetCReal( CurrentModel % Solver % Values,&
          'Friction Newton Relaxation Factor',GotRelax )
      IF(.NOT. GotRelax) NewtonRelax = 1.0_dp    
      
      str = ListGetString( CurrentModel % Solver % Values,'Normal Vector Name',Found )
      IF(.NOT. Found) str = 'Normal Vector'
      NrmSol => VariableGet( CurrentModel % Solver % Mesh % Variables, str, ThisOnly = .TRUE.) 

      CALL ListInitElementVariable( Velo_v )

      str = ListGetString( CurrentModel % Solver % Values,'Height Variable Name',Found )
      IF(.NOT. Found) str = 'height'            
      HeightVar => VariableGet( CurrentModel % Mesh % Variables, str) 

      SlipCoeffVar => VariableGet( CurrentModel % Mesh % Variables,'Slip Coefficient',ThisOnly=.TRUE.)
      SaveSlipCoeff = ASSOCIATED(SlipCoeffVar)
      IF(SaveSlipCoeff) SlipCoeffVar % Values = 0.0_dp

      SaveSlipSpeed = .FALSE.
      IF( SaveSlipCoeff ) THEN
        SlipSpeedVar => VariableGet( CurrentModel % Mesh % Variables,'Slip Speed',ThisOnly=.TRUE.)
        SaveSlipSpeed = ASSOCIATED(SlipSpeedVar)
        IF(SaveSlipSpeed) SlipSpeedVar % Values = 0.0_dp
        
        SlipWeightVar => VariableGet( CurrentModel % Mesh % Variables,'Slip Weight',ThisOnly=.TRUE.)
        SaveSlipWeight = ASSOCIATED(SlipWeightVar)
        IF(SaveSlipWeight) SlipWeightVar % Values = 0.0_dp        
      END IF
      
      InitHandles = .FALSE.
    END IF

    IF( ALLOCATED( Basis ) ) THEN
      IF( SIZE( Basis ) < nd ) THEN
        DEALLOCATE( Basis, dBasisdx ) 
      END IF
    END IF

    IF( .NOT. ALLOCATED( Basis ) ) THEN
      ALLOCATE( Basis(nd), dBasisdx(nd,3) )
    END IF

    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    JAC = 0.0d0
    FORCE = 0.0d0
    dofs = 2
    c = dofs    
    
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    ngp = IP % n
    
    ! There is no elemental routine for this.
    ! So whereas this breaks the beauty it does not cost too much.
    BC => GetBC()     
    HaveFrictionW = ListCheckPresent( BC,'Weertman Friction Coefficient') 
    HaveFrictionU = ListCheckPresent( BC,'Friction Coefficient')
    HaveFriction = HaveFrictionU .OR. HaveFrictionW

    IF( HaveFriction ) THEN
      wut0 = ListGetElementReal( FrictionUt0_h, Element = Element )
    END IF


    CALL GetLocalSolution( NodalHeight, UVariable = HeightVar ) 

    
    DO t=1,ngp      
      !------------------------------------------------------------------------------
      !    Basis function values & derivatives at the integration point
      !------------------------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,IP % u(t),IP % v(t),IP % w(t), detJ, Basis )

      s = detJ * IP % s(t)

      ! Given force on a boundary componentwise
      !----------------------------------------
      SurfaceTraction = ListGetElementReal3D( SurfaceTraction_h, Basis, Element, HaveForce, GaussPoint = t )      

      ! Given force to the normal direction
      ! We have a special internal pressure bacause the hydrostatic version is different from the standard version. 
      !-----------------------------------------------------------------------------------------------------------
      ExtPressure = ListGetElementReal( ExtPressure_h, Basis, Element, HavePres, GaussPoint = t ) &
          - ListGetElementReal( IntPressure_h, Basis, Element, Found, GaussPoint = t )
      HavePres = HavePres .OR. Found
      
      ! Slip coefficient
      !----------------------------------
      SlipCoeff = ListGetElementReal3D( SlipCoeff_h, Basis, Element, HaveSlip, GaussPoint = t )      

      ! Nothing to do, exit the routine
      !---------------------------------
      IF(.NOT. (HaveForce .OR. HavePres .OR. HaveSlip .OR. HaveFriction )) RETURN

      ! Calculate normal vector only if needed.
      ! Slip and friction are assumed to be fully horizontal. 
      IF( HavePres ) THEN
        Normal = ConsistentNormalVector( CurrentModel % Solver, NrmSol, Element, Found, Basis = Basis )
        IF(.NOT. Found) Normal = NormalVector( Element, Nodes, IP % u(t), IP % v(t),.TRUE. )
      END IF

      !-----------------------------------------------------------------
      IF( HaveFriction .AND. .NOT. HaveSlip ) THEN
        IF( HaveSlip ) THEN
          CALL Fatal('IncompressibleNSVec','You cannot combine different friction models!')
        END IF

        ! Velocity at integration point for nonlinear friction laws
        Velo = ListGetElementVectorSolution( Velo_v, Basis, Element, dofs = dofs )

        ! For hyrdostatic model it is assumed that friction is always in (x,y) plane.
        ut = MAX(wut0, SQRT(SUM(Velo(1:dofs)**2)))

        IF( HaveFrictionW ) THEN
          ! Weertman friction law computed internally
          wcoeff = ListGetElementReal( WeertmanCoeff_h, Basis, Element, GaussPoint = t )
          wexp = ListGetElementReal( WeertmanExp_h, Basis, Element, GaussPoint = t )
          TanFrictionCoeff = MIN(wcoeff * ut**(wexp-1.0_dp),1.0e20)
          ! dTanFrictionCoeff/dut for Newton
          TanFder=0._dp
          IF ((ut > wut0).AND.(TanFrictionCoeff < 1.0e20)) &
              TanFder = (wexp-1.0_dp) * wcoeff * ut**(wexp-2.0_dp) 
        ELSE
          ! Else, user defined friction law
          DummyVals(1) = ut          
          TanFrictionCoeff = ListGetElementReal( FrictionCoeff_h, Basis, Element, &
              GaussPoint = t, DummyVals = DummyVals )             
        END IF
        SlipCoeff = TanFrictionCoeff
        HaveSlip = .TRUE.

        IF(SaveSlipSpeed) THEN
          DO i=1,n
            j = SlipSpeedVar % Perm(Element % NodeIndexes(i))
            IF(j>0) THEN
              SlipSpeedVar % Values(j) = SlipSpeedVar % Values(j) + &
                  Basis(i) * IP % s(t) * ut
            END IF
          END DO
        END IF
      END IF

      IF(SaveSlipCoeff) THEN
        DO i=1,n
          j = SlipCoeffVar % Perm(Element % NodeIndexes(i))
          IF(j>0) THEN
            SlipCoeffVar % Values(j) = SlipCoeffVar % Values(j) + &
                Basis(i) * IP % s(t) * MAXVAL(SlipCoeff(1:dim))
          END IF
        END DO
        IF(SaveSlipWeight) THEN
          DO i=1,n
            j = SlipWeightVar % Perm(Element % NodeIndexes(i))
            IF(j>0) SlipWeightVar % Values(j) = SlipWeightVar % Values(j) + Basis(i) * IP % s(t)
          END DO
        END IF
      END IF
                        
      ! Project external pressure to the normal direction
      ! Current 3rd component is neglected!
      IF( HavePres ) THEN
        SurfaceTraction = SurfaceTraction + ExtPressure * Normal
        HaveForce = .TRUE. 
      END IF

      ! Assemble the slip coefficients to the stiffness matrix
      ! This is always in (x,y) plane
      IF( HaveSlip ) THEN               
        DO p=1,nd
          DO q=1,nd
            DO i=1,dofs
              STIFF( (p-1)*c+i,(q-1)*c+i ) = &
                  STIFF( (p-1)*c+i,(q-1)*c+i ) + &
                  s * SlipCoeff(i) * Basis(q) * Basis(p)

              IF(HaveFrictionW .AND. FrictionNewton) THEN
                DO j=1,dofs
                  JAC((p-1)*c+i,(q-1)*c+j ) = JAC((p-1)*c+i,(q-1)*c+j ) + &
                      s * TanFder * Basis(q) * Basis(p) * velo(i) * velo(j) / ut
                END DO
              END IF

            END DO
          END DO
        END DO
      END IF

      ! Assemble given forces to r.h.s.
      IF( HaveForce .OR. HavePres ) THEN
        DO i=1,dim
          DO q=1,nd
            k = (q-1)*c + i
            IF(i <= 2 ) THEN
              ! Components in (x,y) plane
              FORCE(k) = FORCE(k) + s * Basis(q) * SurfaceTraction(i)
            ELSE
#if 0
              ! I thought I would follow the logic on how gravity is mapped but obviously
              ! this is a different story. 
              ! Project the z-force component also to (x,y) plane.
              DO j=1,2
                k = (q-1)*c + j
                FORCE(k) = FORCE(k) + s * Basis(q) * SUM(dBasisdx(1:n,j)*NodalHeight(1:n)) * SurfaceTraction(i)
              END DO
#endif
            END IF
          END DO
        END DO
      END IF
    END DO

    IF(HaveFrictionW .AND. FrictionNewton) THEN
      CALL GetLocalSolution( NodalSol )
      SOL=0._dp
      DO i = 1, c
        SOL(i::c) = NodalSol(i,1:nd)
      END DO

      IF(GotRelax) JAC = NewtonRelax * JAC

      STIFF=STIFF+JAC
      FORCE=FORCE + MATMUL(JAC,SOL)
    END IF

    CALL DefaultUpdateEquations( STIFF, FORCE )

  END SUBROUTINE LocalBoundaryMatrix


!------------------------------------------------------------------------------
! Compute vector for computing du_z/dz and corresponding weight.
!------------------------------------------------------------------------------
  SUBROUTINE LocalDuz(Element, n, ntot, duz, wuz, FirstElem, ub, dpr )
!------------------------------------------------------------------------------
    USE LinearForms
    IMPLICIT NONE

    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    INTEGER, INTENT(IN) :: n, ntot
    REAL(KIND=dp), POINTER :: duz(:), wuz(:), ub(:)
    LOGICAL :: FirstElem
    REAL(KIND=dp), POINTER, OPTIONAL :: dpr(:)
!------------------------------------------------------------------------------
    INTEGER :: dim, dofs
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp), POINTER :: muVec(:)
    TYPE(Nodes_t) :: Nodes
    LOGICAL :: Stat, Found, DoVisc
    INTEGER :: NodalPerm(ntot)
    REAL(KIND=dp) :: NodalVelo(2,ntot), duz_elem(ntot), wuz_elem(ntot), dp_elem(ntot), ub_elem(ntot), &
        vgrad(2), w, velo(2), zgrad(2)
    REAL(KIND=dp), ALLOCATABLE, SAVE :: BasisVec(:,:), dBasisdxVec(:,:,:), DetJVec(:)
    INTEGER :: t, i, j, k, ngp, allocstat


!DIR$ ATTRIBUTES ALIGN:64 :: BasisVec, dBasisdxVec, DetJVec
!$OMP THREADPRIVATE(BasisVec, dBasisdxVec, DetJVec, Nodes)

    SAVE Nodes
!------------------------------------------------------------------------------

    CALL GetElementNodesVec( Nodes )

    dofs = 2
    dim = 3
    
    ! Numerical integration:
    !-----------------------
    IP = GaussPointsAdapt(Element)
    ngp = IP % n

    DoVisc = PRESENT(dpr)
    
    ! Deallocate storage if needed 
    IF (ALLOCATED(BasisVec)) THEN
      IF (SIZE(BasisVec,1) < ngp .OR. SIZE(BasisVec,2) < ntot) &
          DEALLOCATE(BasisVec,dBasisdxVec, DetJVec )
    END IF
    
    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(BasisVec)) THEN
      ALLOCATE(BasisVec(ngp,ntot), dBasisdxVec(ngp,ntot,3), DetJVec(ngp), STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal('HydrostaticNSSolver::LocalBulkMatrix','Local storage allocation failed')
      END IF
    END IF
    
    ! Vectorized basis functions
    stat = ElementInfoVec(Element, nodes, ngp, IP % U, IP % V, &
        IP % W, detJvec, SIZE(basisVec, 2), BasisVec, dBasisdxVec)
    
    ! Weights at integration points
    DO t = 1, ngp
      DetJVec(t) = DetJVec(t) * IP % s(t)
    END DO
 
    ! Get the velocity at nodes
    CALL GetLocalSolution( NodalVelo )

    ! Return the effective viscosity. Currently only non-newtonian models supported.
    IF(DoVisc) THEN
      dp_elem = 0.0_dp
      muvec => EffectiveViscosityVec( ngp, ntot, BasisVec, dBasisdxVec, Element, NodalVelo, &
          InitHandles = FirstElem )
    END IF
      
    duz_elem = 0.0_dp
    wuz_elem = 0.0_dp
    ub_elem = 0.0_dp
    
    DO t = 1, ngp
      DO i = 1, 2
        vgrad(i) = SUM(dBasisdxVec(t,1:ntot,i)*NodalVelo(i,1:ntot))
      END DO
      DO i = 1, 2
        velo(i) = SUM(BasisVec(t,1:ntot)*NodalVelo(i,1:ntot))
        zgrad(i) = SUM(dBasisdxVec(t,1:ntot,i)*Nodes % z(1:ntot))
      END DO
      DO i=1,ntot
        ! It would seem that weighting with the DetJ would make sense, but maybe not...
        !w = DetJVec(t) * BasisVec(t,i)
        w = BasisVec(t,i) * IP % s(t)

        ! Negative sign comes from continuity equation: du_z/dz = -(du_x/dx + du_y/dy)
        duz_elem(i) = duz_elem(i) - w * SUM(vgrad(1:2))
        wuz_elem(i) = wuz_elem(i) + w

        ub_elem(i) = ub_elem(i) + w * SUM(velo*zgrad)
        
        ! This adds the correction to pressure from equation (5.63)
        IF(DoVisc) THEN
          dp_elem(i) = dp_elem(i) - w * 2 * muvec(t) * SUM(vgrad(1:2)) 
        END IF
      END DO
    END DO

    ! Global degrees of freedom
    NodalPerm(1:ntot) = CurrentModel % Solver % Variable % Perm(Element % NodeIndexes )

    ! Sum up the contribution
    duz(NodalPerm(1:ntot)) = duz(NodalPerm(1:ntot)) + duz_elem(1:ntot)
    wuz(NodalPerm(1:ntot)) = wuz(NodalPerm(1:ntot)) + wuz_elem(1:ntot)

    ub(NodalPerm(1:ntot)) = ub(NodalPerm(1:ntot)) + ub_elem(1:ntot)
    
    IF( DoVisc ) THEN
      dpr(NodalPerm(1:ntot)) = dpr(NodalPerm(1:ntot)) + dp_elem(1:ntot)
    END IF
    
    
  END SUBROUTINE LocalDuz

!------------------------------------------------------------------------------

  
  SUBROUTINE PopulateDerivedFields()

    IMPLICIT NONE

    TYPE(Variable_t), POINTER :: VarXY, VarFull, VarDuz, VarP, VarVx, VarVy
    CHARACTER(LEN=MAX_NAME_LEN):: str
    TYPE(ValueList_t), POINTER :: Params, Material
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: Found, PressureCorr, BaseVelo
    REAL(KIND=dp), POINTER :: duz(:), wuz(:), dpr(:), ub(:)
    INTEGER :: i,j,k,j1,j2,i1,i2,k1,k2,t,n,nd,nb,active, dofs, pdof, zdof
    TYPE(Element_t), POINTER :: Element
    TYPE(Solver_t), POINTER :: pSolver
    REAL(KIND=dp) :: dz, rho, g, Nrm(3)
    REAL(KIND=dp), POINTER :: gWork(:,:)
    
    
    SAVE :: duz, wuz, dpr, ub

    CALL Info('HydrostaticNSVec','Populating derived fields: vertial velocity and pressure')
    
    pSolver => CurrentModel % Solver    
    Params => pSolver % Values
    VarXY => pSolver % Variable
    Mesh => CurrentModel % Mesh
    dofs = 0
    
    str = ListGetString( Params,'Velocity Vector Name',Found )
    IF(.NOT. Found) str = ListGetString( Params,'Velocity Variable Name',Found )

    IF(.NOT. Found) str = 'Flow Solution'
    Found = .TRUE.
    
    IF(Found) THEN
      VarFull => VariableGet(Mesh % Variables, str, ThisOnly = .TRUE.)
      IF(ASSOCIATED(VarFull)) THEN
        dofs = VarFull % dofs
        IF(dofs < 3 .OR. dofs > 4 ) THEN
          CALL Fatal('HydrostaticNSVec','We need 3 or 4 components for velocity field, not '//I2S(dofs))
        END IF
        zdof = 3
      ELSE
        VarVx => VariableGet(Mesh % Variables, TRIM(str)//' 1', ThisOnly = .TRUE.)
        VarVy => VariableGet(Mesh % Variables, TRIM(str)//' 2', ThisOnly = .TRUE.)
        VarFull => VariableGet(Mesh % Variables, TRIM(str)//' 3', ThisOnly = .TRUE.)
        i = 0
        IF(ASSOCIATED(VarVx)) i=i+1
        IF(ASSOCIATED(VarVy)) i=i+1
        IF(ASSOCIATED(VarFull)) i=i+1
        IF(i<3) THEN                
          CALL Fatal('HydrostaticNSVec','Could not find components of velocity variable: '//TRIM(str))
        END IF
        CALL Info('HydrostaticNSVec','Setting velocity for each component separately',Level=10)
        dofs = 1          
        zdof = 1
      END IF
      VarFull % Values = 0.0_dp      
    END IF
      
    NULLIFY(VarP) 
    str = ListGetString( Params,'Pressure Variable Name',Found )
    IF( Found ) THEN
      VarP => VariableGet(Mesh % Variables, str, ThisOnly = .TRUE.)
      IF(.NOT. ASSOCIATED(VarP)) THEN
        CALL Fatal('HydrostaticNSVec','Could not find full pressure variable: '//TRIM(str))
      END IF
    END IF
    IF(.NOT. ASSOCIATED(VarP)) THEN
      IF(dofs == 4) THEN
        CALL Info('HydrostaticNSVec','Assuming Pressure to be last entry of: '//TRIM(VarFull % Name),Level=6)
        VarP => VarFull
      END IF
    END IF    
    
    ! Pressure can be 1st component of pressure or last compenent of 'flow solution'
    pdof = 0
    IF(ASSOCIATED(VarP)) pdof = VarP % Dofs
    IF(pdof == 0 .AND. dofs == 0 ) RETURN

    g = 1.0_dp
    PressureCorr = .FALSE.
    IF(pdof > 0) THEN
      gWork => ListGetConstRealArray( CurrentModel % Constants,'Gravity',Found)
      IF(Found) THEN
        g = ABS(gWork(SIZE(gWork,1),1))
      ELSE
        CALL Warn('HydrostaticNSVec','"Gravity" not found in simulation section, setting to 1')
      END IF      
      PressureCorr = ListGetLogical( Params,'Pressure Correction',Found ) 
    END IF

    ! Copy the velocity componts x & y
    IF( dofs > 0 ) THEN
      DO i=1,Mesh % NumberOfNodes
        j = VarFull % Perm(i)
        k = VarXY % Perm(i)
        IF(dofs == 1 ) THEN
          VarVx % Values(j) = VarXY % Values(2*k-1)
          VarVy % Values(j) = VarXY % Values(2*k)          
        ELSE
          VarFull % Values(dofs*(j-1)+1) = VarXY % Values(2*k-1)
          VarFull % Values(dofs*(j-1)+2) = VarXY % Values(2*k)
        END IF
      END DO
    END IF
      
    BaseVelo = .TRUE.
        
    n = SIZE(VarXY % Values) / 2
    ALLOCATE(duz(n),wuz(n),ub(n))      
    duz = 0.0_dp
    wuz = 0.0_dp
    ub = 0.0_dp
          
    IF(PressureCorr) THEN
      ALLOCATE(dpr(n))
      dpr = 0.0_dp
    END IF
      
    Active = GetNOFActive()
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes(Element)
      nb = GetElementNOFBDOFs(Element)
      nd = GetElementNOFDOFs(Element)
      
      ! Calculate elemental contribution to d(uz)/dz
      !---------------------------------------------
      IF( PressureCorr ) THEN 
        CALL LocalDuz(Element, n, n, duz, wuz, t==1, ub, dpr )      
      ELSE
        CALL LocalDuz(Element, n, n, duz, wuz, t==1, ub ) 
      END IF
        
      IF(t==1) THEN
        Material => GetMaterial(Element)
        rho = ListGetCReal(Material,'Density',UnfoundFatal=.TRUE.) 
      END IF
    END DO

    ! Communicate the values for the shared dofs prior to normalization
    IF( ParEnv % PEs > 1 ) THEN
      CALL ParallelSumNodalVector( Mesh, wuz, VarXY % Perm )
      CALL ParallelSumNodalVector( Mesh, duz, VarXY % Perm )
      CALL ParallelSumNodalVector( Mesh, ub, VarXY % Perm )
      IF( PressureCorr ) THEN
        CALL ParallelSumNodalVector( Mesh, dpr, VarXY % Perm )
      END IF
    END IF
            
    WHERE(wuz > EPSILON(dz) )
      duz = duz / wuz
    END WHERE
    WHERE(wuz > EPSILON(dz) )
      ub = ub / wuz
    END WHERE
    IF(PressureCorr) THEN
      WHERE( wuz > EPSILON(dz))
        dpr = dpr / wuz
      END WHERE
    END IF
    
    ! If we have exported variable "duz" we can save the value of the field
    VarDuz => VariableGet(Mesh % Variables,'duz', ThisOnly = .TRUE.)
    IF( ASSOCIATED( VarDuz ) ) VarDuz % Values = duz

    ! Integrate over structured mesh 
    DO i=1,Mesh % NumberOfNodes                   
      IF(DownPointer(i)==0) CYCLE
      i1 = i      
            
      ! Nodes at the bedrock
      IF(dofs > 0 .AND. DownPointer(i) == i) THEN
        ! Initailize values at the bedrock
        k1 = VarFull % Perm(i1)
        IF(k1>0) THEN
          VarFull % Values(dofs*(k1-1)+zdof) = 1.0_dp * ub(k1) 
        END IF
          
        DO WHILE(.TRUE.)
          i2 = UpPointer(i1)
          IF(i2==i1) EXIT

          j1 = VarXY % Perm(i1)
          j2 = VarXY % Perm(i2)
          
          k1 = VarFull % Perm(i1)
          k2 = VarFull % Perm(i2)

          IF(k1>0 .AND. k2>0) THEN
            ! Integrate for the z-velocity
            dz = Mesh % Nodes % z(i2) - Mesh % Nodes % z(i1)
            VarFull % Values(dofs*(k2-1)+zdof) = VarFull % Values(dofs*(k1-1)+zdof) + &
                dz * ( duz(j1) + duz(j2) ) / 2.0_dp
          END IF
          i1 = i2
        END DO
      END IF

      ! Nodes at the ice top
      IF(pdof > 0 .AND. UpPointer(i) == i) THEN
        ! Initailize values at ice top
        k1 = VarP % Perm(i1)
        IF(k1>0) VarP % Values(pdof*k1) = 0.0_dp
        
        DO WHILE(.TRUE.)
          i2 = DownPointer(i1)
          IF(i2==i1) EXIT

          k1 = VarP % Perm(i1)
          k2 = VarP % Perm(i2)

          IF(k1>0 .AND. k2>0) THEN
            ! Integrate for the hydrostatic pressure           
            dz = Mesh % Nodes % z(i1) - Mesh % Nodes % z(i2)            
            VarP % Values(pdof*k2) = VarP % Values(pdof*k1) + g * dz * rho
          END IF
          i1 = i2
        END DO
      END IF
    END DO

    DEALLOCATE(duz,wuz)
    
    IF( PressureCorr ) THEN
      DO i=1,Mesh % NumberOfNodes
        j = VarP % Perm(i)
        k = VarXY % Perm(i)
        VarP % Values(pdof*j) = VarP % Values(pdof*j) + dpr(k)
      END DO
      DEALLOCATE(dpr)
    END IF
     
  END SUBROUTINE PopulateDerivedFields

 
  ! Initialize height above z=0 such that we can access it from anywhere, not only at top.
  ! Structured mesh is assumed in z-direction.
  !---------------------------------------------------------------------------------------
  SUBROUTINE InitializeHeightField()

    IMPLICIT NONE

    TYPE(Variable_t), POINTER :: VarH
    CHARACTER(LEN=MAX_NAME_LEN):: str
    TYPE(ValueList_t), POINTER :: Params   
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: Found
    INTEGER :: i,j,k,i0,i1,i2,k0,k1,k2
    TYPE(Solver_t), POINTER :: pSolver
    REAL(KIND=dp) :: hloc
        

    pSolver => CurrentModel % Solver    
    Params => pSolver % Values
    Mesh => pSolver % Mesh

    str = ListGetString( Params,'Height Variable Name',Found )
    IF(.NOT. Found) str = 'height'
    VarH => VariableGet(Mesh % Variables, str, ThisOnly = .TRUE.)
    IF(.NOT. ASSOCIATED(VarH)) THEN
      CALL Fatal('HydrostaticNSVec','Could not find height variable: '//TRIM(str))
    END IF
    CALL Info('HydrostaticNSVec','Integrating for the height on the stuctured direction!')
              
    ! Integrate over structured mesh 
    DO i=1,Mesh % NumberOfNodes                   
      ! Found a node at the bedrock
      IF(DownPointer(i) == i) THEN
        hloc = 0.0_dp
        i1 = i
        ! Find the top and compute the height 
        DO WHILE(.TRUE.)
          i2 = UpPointer(i1)
          IF(i2==i1) THEN
            hloc = Mesh % Nodes % z(i1) ! - Mesh % Nodes % z(i)
            EXIT
          END IF
          i1 = i2
        END DO

        ! Map the height to all nodes on the stride
        i1 = i
        DO WHILE(.TRUE.)
          k1 = VarH % Perm(i1)
          IF(k1>0) VarH % Values(k1) = hloc
          i2 = UpPointer(i1)
          IF(i2==i1) EXIT
          i1 = i2
        END DO
      END IF
    END DO
    
  END SUBROUTINE InitializeHeightField

END MODULE HydrostaticNSUtils


!------------------------------------------------------------------------------
SUBROUTINE HydrostaticNSSolver_Init0(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------  
  LOGICAL :: Found, Serendipity

  Serendipity = GetLogical( GetSimulation(), 'Serendipity P Elements', Found)
  IF(.NOT.Found) Serendipity = .TRUE.

  CALL ListAddNewString(GetSolverParams(),'Element','p:1')
  
!------------------------------------------------------------------------------
END SUBROUTINE HydrostaticNSSolver_Init0
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
SUBROUTINE HydrostaticNSSolver_init(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------  
  TYPE(ValueList_t), POINTER :: Params 
  LOGICAL :: Found
  INTEGER :: dim
  CHARACTER(*), PARAMETER :: Caller = 'HydrostaticNSSolver_init'
!------------------------------------------------------------------------------ 
  Params => GetSolverParams() 
  
  IF( ListCheckPresentAnyBC( Model, 'Pressure 1' ) ) THEN
    CALL Fatal( Caller,'Use >Surface Traction 1< instead of >Pressure 1<')
  END IF
  IF( ListCheckPresentAnyBC( Model, 'Pressure 2' ) ) THEN
    CALL Fatal( Caller,'Use >Surface Traction 3< instead of >Pressure 2<')
  END IF
  IF( ListCheckPresentAnyBC( Model, 'Pressure 3' ) ) THEN
    CALL Fatal( Caller,'Use >Surface Traction 3< instead of >Pressure 3<')
  END IF

  IF( ListCheckPresentAnyBC( Model, 'Slip Coefficient 2' ) ) THEN
    IF( .NOT. ListCheckPresentAnyBC( Model, 'Slip Coefficient 1' ) ) THEN
      CALL Fatal(Caller,'This solver expects "Slip Coefficient i" for i={1,2}')   
    END IF
  END IF
  
  dim = CoordinateSystemDimension()
  IF( dim /= 3 ) THEN
    CALL Fatal(Caller,'This solver in only applicable in 3D!')
  END IF
    
  CALL ListAddNewString(Params, 'Variable','HorizontalVelocity[HorizontalVelocity:2]')

  ! The eliminated coordinate is always the 3rd one.
  ! This is needed by the extruded detection.
  CALL ListAddNewInteger(Params,'Active Coordinate',3)
  
  ! It makes sense to eliminate the bubbles to save memory and time
  CALL ListAddNewLogical(Params, 'Bubbles in Global System', .FALSE.)

  IF( GetLogical( Params,'Save Viscosity', Found ) ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),'Viscosity')
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),'Shearrate')
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),'-nooutput Viscosity Weight')
  END IF

  IF( GetLogical( Params,'Save Slip', Found ) ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),'Slip Coefficient')
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),'Slip Speed')
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),'-nooutput Slip Weight')
  END IF
  
!------------------------------------------------------------------------------ 
END SUBROUTINE HydrostaticNSSolver_Init
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE HydrostaticNSSolver(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  USE HydrostaticNSUtils
  USE MainUtils

  
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER :: Params 
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Solver_t), POINTER :: pSolver 
  TYPE(GaussIntegrationPoints_t) :: IP
  INTEGER :: Element_id
  INTEGER :: i, n, nb, nd, dim, Active, maxiter, iter
  REAL(KIND=dp) :: Norm
  LOGICAL :: Found, Converged
  LOGICAL :: SpecificLoad, InitBCHandles
  LOGICAL, SAVE :: InitDone = .FALSE.
  CHARACTER(*), PARAMETER :: Caller = 'HydrostaticNSSolver'

!------------------------------------------------------------------------------
! Local variables to be accessed by the contained subroutines:
!------------------------------------------------------------------------------
  LOGICAL :: LinearAssembly, Newton
!------------------------------------------------------------------------------ 

  dim = CoordinateSystemDimension()
  Mesh => GetMesh()
  Params => GetSolverParams() 
  pSolver => Solver
   
  IF(.NOT. InitDone ) THEN  
    ! We detect the extruded structure only for the first time.
    CALL DetectExtrudedStructure( Mesh, PSolver, &
        UpNodePointer = UpPointer, DownNodePointer = DownPointer)
    InitDone = .TRUE.
  END IF

  ! The height must be consistent with the current mesh. 
  IF( .NOT. ListGetLogical( Solver % Values,'Skip Height Initialization',Found) ) THEN
    CALL InitializeHeightField()
  END IF

  CALL DefaultStart()
  
  !-----------------------------------------------------------------------------
  ! Output the number of integration points as information.
  ! This in not fully informative if several element types are present.
  !-----------------------------------------------------------------------------
  Element => Mesh % Elements( Solver % ActiveElements(1) ) 
  IP = GaussPointsAdapt( Element )
  CALL Info(Caller,'Number of 1st integration points: '//I2S(IP % n), Level=5)
  
  !-----------------------------------------------------------------------------
  ! Set the flags/parameters which define how the system is assembled: 
  !-----------------------------------------------------------------------------
  LinearAssembly = GetLogical(Params, 'Linear Equation', Found )
  SpecificLoad = GetLogical(Params,'Specific Load',Found)
  
  Maxiter = GetInteger(Params, 'Nonlinear system max iterations', Found)
  IF (.NOT.Found) Maxiter = 1
  !-----------------------------------------------------------------------------

  
  DO iter=1,maxiter

    CALL Info(Caller,'--------------------------------------------------------', Level=4)
    CALL Info(Caller, 'Nonlinear iteration: '//I2S(iter), Level=4)
    CALL Info(Caller,'--------------------------------------------------------', Level=4)
    
100 CONTINUE
    
    Active = GetNOFActive()
    CALL DefaultInitialize()
    
    Newton = GetNewtonActive( Solver )

    DO Element_id=1,1
      Element => GetActiveElement(Element_id)
      n  = GetElementNOFNodes(Element)
      !
      ! When the number of bubbles is obtained with the Update=.TRUE. flag,
      ! we need to call GetElementNOFBDOFs before calling GetElementNOFDOFs.
      !
      nb = GetElementNOFBDOFs(Element)
      nd = GetElementNOFDOFs(Element)

      ! Get element local matrix and rhs vector:
      !-----------------------------------------
      CALL LocalBulkMatrix(Element, n, nd, nd+nb, &
          SpecificLoad, LinearAssembly, nb, Newton, .TRUE.)
    END DO

    !$OMP PARALLEL SHARED(Active, dim, SpecificLoad, &
    !$OMP                 dt, LinearAssembly, Newton ) &
    !$OMP PRIVATE(Element, Element_id, n, nd, nb)  DEFAULT(None)
    !$OMP DO    
    DO Element_id=2,Active
      Element => GetActiveElement(Element_id)
      n  = GetElementNOFNodes(Element)
      nb = GetElementNOFBDOFs(Element, Update=.TRUE.)
      nd = GetElementNOFDOFs(Element)
      
      ! Get element local matrix and rhs vector:
      !-----------------------------------------
      CALL LocalBulkMatrix(Element, n, nd, nd+nb, &
          SpecificLoad, LinearAssembly, nb, Newton, .FALSE.)
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    CALL DefaultFinishBulkAssembly()
    
    Active = GetNOFBoundaryElements()
    InitBCHandles = .TRUE.  
    DO Element_id=1,Active
      Element => GetBoundaryElement(Element_id) 
      IF (ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()

        IF ( GetElementDim() < 2 ) CYCLE

        ! Get element local matrix and rhs vector:
        !-----------------------------------------
        CALL LocalBoundaryMatrix(Element, n, nd, dim, InitBCHandles , Newton )
        InitBCHandles = .FALSE.
      END IF
    END DO
    
    CALL DefaultFinishBoundaryAssembly()

    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    ! Check stepsize for nonlinear iteration
    !------------------------------------------------------------------------------
    IF( DefaultLinesearch( Converged ) ) GOTO 100
    IF( Converged ) EXIT
    
    Norm = DefaultSolve()

    IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
  END DO

  CALL DefaultFinish()

  BLOCK
    TYPE(Variable_t), POINTER, SAVE :: pVar, wVar
    REAL(KIND=dp) :: minw
    minw = 1.0e-20
    DO i=1,4
      SELECT CASE(i)
      CASE( 1 ) 
        wVar => VariableGet( CurrentModel % Mesh % Variables,'Slip Weight',ThisOnly=.TRUE.)
        pVar => VariableGet( CurrentModel % Mesh % Variables,'Slip Coefficient',ThisOnly=.TRUE.)
      CASE( 2 ) 
        pVar => VariableGet( CurrentModel % Mesh % Variables,'Slip Speed',ThisOnly=.TRUE.)
      CASE( 3 )         
        wVar => VariableGet( CurrentModel % Mesh % Variables,'Viscosity Weight',ThisOnly=.TRUE.)
        pVar => VariableGet( CurrentModel % Mesh % Variables,'Viscosity',ThisOnly=.TRUE.)
      CASE( 4 ) 
        pVar => VariableGet( CurrentModel % Mesh % Variables,'Strainrate',ThisOnly=.TRUE.)
      END SELECT
      IF(ASSOCIATED(wVar) .AND. ASSOCIATED(pVar) ) THEN     
        CALL Info('IncompressibleNSSolver','Normalizing field number: '//I2S(i),Level=15)
        WHERE(wVar % Values > minw )
          pVar % Values = pVar % Values / wVar % Values
        END WHERE
      END IF
    END DO
  END BLOCK
  
  ! Compute all the other fields that can be computed from the solution.
  CALL PopulateDerivedFields()
  
  CALL Info( Caller,'All done',Level=10)
!------------------------------------------------------------------------------
END SUBROUTINE HydrostaticNSSolver
!------------------------------------------------------------------------------
