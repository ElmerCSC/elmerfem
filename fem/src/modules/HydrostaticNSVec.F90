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

CONTAINS

  ! Compute effective viscosity. This is needed not only by the bulkassembly but also 
  ! by the postprocessing when estimating pressure.
  !----------------------------------------------------------------------------------
  FUNCTION EffectiveViscosityVec( ngp, ntot, BasisVec, dBasisdxVec, Element, NodalVelo, &
      ViscDerVec, ViscNewton, InitHandles ) RESULT ( EffViscVec ) 

    IMPLICIT NONE 
    
    INTEGER :: ngp,ntot
    REAL(KIND=dp) :: BasisVec(:,:), dBasisdxVec(:,:,:)
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: NodalVelo(:,:)
    REAL(KIND=dp), ALLOCATABLE :: ViscDerVec(:)
    LOGICAL :: InitHandles , ViscNewton
    REAL(KIND=dp), POINTER  :: EffViscVec(:)

    INTEGER :: allocstat,i,j,k,dim,dofs
    LOGICAL :: Found     
    CHARACTER(LEN=MAX_NAME_LEN) :: ViscModel
    TYPE(ValueHandle_t), SAVE :: Visc_h, ViscModel_h, ViscExp_h, ViscCritical_h, &
        ViscNominal_h, ViscDiff_h, ViscTrans_h, ViscYasuda_h, ViscGlenExp_h, ViscGlenFactor_h, &
        ViscArrSet_h, ViscArr_h, ViscTLimit_h, ViscRate1_h, ViscRate2_h, ViscEne1_h, ViscEne2_h, &
        ViscTemp_h
    REAL(KIND=dp), SAVE :: R, vgrad(2,3)
    REAL(KIND=dp) :: c1, c2, c3, c4, Tlimit, ArrheniusFactor, A1, A2, Q1, Q2, ViscCond
    LOGICAL, SAVE :: ConstantVisc = .FALSE., Visited = .FALSE.
    REAL(KIND=dp), ALLOCATABLE, SAVE :: ss(:), s(:), ArrheniusFactorVec(:)
    REAL(KIND=dp), POINTER, SAVE :: ViscVec0(:), ViscVec(:), TempVec(:), EhfVec(:) 
    TYPE(Variable_t), POINTER, SAVE :: ShearVar, ViscVar
    LOGICAL, SAVE :: SaveShear, SaveVisc

    !$OMP THREADPRIVATE(ss,s,ViscVec0,ViscVec,ArrheniusFactorVec)

    dim = 3
    dofs = 2
    
    IF(InitHandles ) THEN
      CALL Info('EffectiveViscosityVec','Initializing handles for viscosity models',Level=8)

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
              CALL Fatal('EffectiveViscosityVec','Replace >Constant Temperature< with >Relative Temperature<')
            END IF
            IF( ListCheckPresentAnyMaterial( CurrentModel,'Temperature Field Variable') ) THEN
              CALL Fatal('EffectiveViscosityVec','Replace >Temperature Field Variable< with >Relative Temperature<')
            END IF
          END IF
          IF( ListCheckPresentAnyMaterial( CurrentModel,'Glen Enhancement Factor Function')  ) THEN
            CALL Fatal('EffectiveViscosityVec','No Glen function API yet!')
          END IF
          R = GetConstReal( CurrentModel % Constants,'Gas Constant',Found)
          IF (.NOT.Found) R = 8.314_dp
        END IF
      ELSE
        CALL Info('EffectiveViscosityVec','Using constant viscosity!')
      END IF

      ShearVar => VariableGet( CurrentModel % Mesh % Variables,'Shearrate',ThisOnly=.TRUE.)
      SaveShear = ASSOCIATED(ShearVar)
      IF(SaveShear) THEN
        IF(ShearVar % TYPE == Variable_on_gauss_points ) THEN
          CALL Info('EffectiveViscosityVec','Saving "Shearrate" on ip points!',Level=10)
        ELSE IF( ShearVar % TYPE == Variable_on_elements ) THEN
          CALL Info('EffectiveViscosityVec','Saving "Shearrate" on elements!',Level=10)
        ELSE
          CALL Fatal('EffectiveViscosityVec','"Shearrate" should be either ip or elemental field!')
        END IF
      END IF

      ViscVar => VariableGet( CurrentModel % Mesh % Variables,'Viscosity',ThisOnly=.TRUE.)
      SaveVisc = ASSOCIATED(ViscVar)        
      IF(SaveVisc) THEN
        IF(ViscVar % TYPE == Variable_on_gauss_points ) THEN
          CALL Info('EffectiveViscosityVec','Saving "Viscosity" on ip points!',Level=10)
        ELSE IF( ViscVar % TYPE == Variable_on_elements ) THEN
          CALL Info('EffectiveViscosityVec','Saving "Viscosity" on elements!',Level=10)
        ELSE
          CALL Fatal('EffectiveViscosityVec','"Viscosity" should be either ip or elemental field!')
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
    IF( ViscNewton ) THEN
      ViscDerVec(1:ngp) = 0.0_dp
    END IF

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
        CALL Fatal('EffectiveViscosityVec','Local storage allocation failed')
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
      i = Element % ElementIndex
      IF( ShearVar % TYPE == Variable_on_gauss_points ) THEN
        j = ShearVar % Perm(i+1) - ShearVar % Perm(i)
        IF(j /= ngp) THEN
          CALL Fatal('EffectiveViscosityVec','Expected '//I2S(j)//' gauss point for "Shearrate" got '//I2S(ngp))
        END IF
        ShearVar % Values(ShearVar % Perm(i)+1:ShearVar % Perm(i+1)) = ss(1:ngp)
      ELSE
        ShearVar % Values(ShearVar % Perm(i)) = SUM(ss(1:ngp)) / ngp
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

        IF( ViscNewton ) THEN
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

        IF( ViscNewton ) THEN
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

      IF (ViscNewton ) THEN
        WHERE(ss(1:ngp) /= 0) ViscDerVec(1:ngp) = &
            ViscVec0(1:ngp) * (c2-1)/2 * ss(1:ngp)**((c2-1)/2-1)
      END IF

      c4 = ListGetElementReal( ViscNominal_h,Element=Element,Found=Found)
      IF( Found ) THEN
        ViscVec(1:ngp) = ViscVec(1:ngp) / c4**(c2-1)
        IF (ViscNewton ) THEN
          ViscDerVec(1:ngp) = ViscDerVec(1:ngp) / c4**(c2-1)
        END IF
      END IF

    CASE('power law too')
      c2 = ListGetElementReal( ViscExp_h,Element=Element)           
      ViscVec(1:ngp) = ViscVec0(1:ngp)**(-1/c2)* ss(1:ngp)**(-(c2-1)/(2*c2)) / 2

      IF (ViscNewton ) THEN
        ViscDerVec(1:ngp) = ViscVec0(1:ngp)**(-1/c2)*(-(c2-1)/(2*c2))*ss(1:ngp)*(-(c2-1)/(2*c2)-1) / 2
      END IF

    CASE ('carreau')      
      c1 = ListGetElementReal( ViscDiff_h,Element=Element)
      c2 = ListGetElementReal( ViscExp_h,Element=Element)
      c3 = ListGetElementReal( ViscTrans_h,Element=Element)
      c4 = ListGetElementReal( ViscYasuda_h,Element=Element,Found=Found)
      IF( Found ) THEN
        ViscVec(1:ngp) = ViscVec0(1:ngp) + c1 * (1 + c3**c4*ss(1:ngp)**(c4/2))**((c2-1)/c4) 

        IF( ViscNewton ) THEN
          ViscDerVec(1:ngp) = c1*(1+c3**c4*ss(1:ngp)**(c4/2))**((c2-1)/c4-1)*(c2-1)/2*c3**c4*&
              ss(1:ngp)**(c4/2-1)
        END IF
      ELSE
        ViscVec(1:ngp) = ViscVec0(1:ngp) + c1 * (1 + c3*c3*ss(1:ngp))**((c2-1)/2) 

        IF( ViscNewton ) THEN
          ViscDerVec(1:ngp) = c1*(c2-1)/2*c3**2*(1+c3**2*ss(1:ngp))**((c2-1)/2-1)
        END IF
      END IF

    CASE ('cross')
      c1 = ListGetElementReal( ViscDiff_h,Element=Element)
      c2 = ListGetElementReal( ViscExp_h,Element=Element)
      c3 = ListGetElementReal( ViscTrans_h,Element=Element)

      ViscVec(1:ngp) = ViscVec0(1:ngp) + c1 / (1 + c3*ss(1:ngp)**(c2/2))

      IF( ViscNewton ) THEN
        ViscDerVec(1:ngp) = -c1*c3*ss(1:ngp)**(c2/2)*c2 / (2*(1+c3*ss(1:ngp)**(c2/2))**2*ss(1:ngp))
      END IF

    CASE ('powell eyring')
      c1 = ListGetElementReal( ViscDiff_h,Element=Element)
      c2 = ListGetElementReal( ViscTrans_h,Element=Element)

      s(1:ngp) = SQRT(ss(1:ngp))

      IF( ViscNewton ) THEN          
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
      CALL Fatal('EffectiveViscosityVec','Unknown material model')

    END SELECT

    IF(SaveVisc) THEN
      i = Element % ElementIndex
      IF( ViscVar % TYPE == Variable_on_gauss_points ) THEN
        j = ViscVar % Perm(i+1) - ViscVar % Perm(i) 
        IF(j /= ngp) THEN
          CALL Fatal('EffectiveViscosityVec','Expected '//I2S(j)//' gauss point for "Viscosity" got '//I2S(ngp))
        END IF
        ViscVar % Values(ViscVar % Perm(i)+1:ViscVar % Perm(i+1)) = ViscVec(1:ngp)
      ELSE
        ViscVar % Values(ViscVar % Perm(i)) = SUM(ViscVec(1:ngp)) / ngp
      END IF
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
    LOGICAL :: Stat, Found

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

    INTEGER :: t, i, j, k, p, q, ngp, allocstat
    INTEGER, SAVE :: elemdim
    INTEGER :: DOFs

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
    IP = GaussPointsAdapt(Element)
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
          ForcePart(ntot), GradVec(ngp,2,2), GradHeight(ngp,dim), &
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
      HeightVar => VariableGet( CurrentModel % Mesh % Variables,'height') 
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
        muDerVec0, Newton,  InitHandles )        
    
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
      DO i = 1,dofs
        DO j = 1,dofs
          DO k = 1, ngp
            GradVec(k, i, j) = SUM(dBasisdxVec(k,1:ntot,j)*NodalVelo(i,1:ntot))
          END DO
        END DO
      END DO
      
      IF (ANY(muDerVec0(1:ngp)/=0)) THEN
        DO i = 1,dofs
          DO j = 1,dofs
            StrainRateVec(1:ngp,i,j) = ( GradVec(1:ngp,i,j) + GradVec(1:ngp,j,i) ) / 2
          END DO
        END DO

        muDerVec0(1:ngp) = muderVec0(1:ngp)*detJVec(1:ngp)*8
        DO i=1,dofs
          DO q = 1,ntot
            g(1:ngp,q,i) = SUM(StrainRateVec(1:ngp,i,:)*dBasisdxvec(1:ngp,q,1:dim),2)
          END DO
        END DO

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
  SUBROUTINE LocalBoundaryMatrix( Element, n, nd, dim, InitHandles, Newton)
    !------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    INTEGER, INTENT(IN) :: n, nd, dim
    LOGICAL, INTENT(INOUT) :: InitHandles 
    LOGICAL :: Newton
    !------------------------------------------------------------------------------    
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp), TARGET :: STIFF(nd*2,nd*2), FORCE(nd*2),NodalHeight(nd)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    INTEGER :: c,i,j,k,l,p,q,t,ngp,dofs
    LOGICAL :: HaveSlip, HaveForce, HavePres, HaveFrictionW, HaveFrictionU, &
        HaveFriction, FrictionNewton, Found, Stat
    REAL(KIND=dp) :: ExtPressure, s, detJ, wut0, wexp, wcoeff, ut
    REAL(KIND=dp) :: SlipCoeff(3), SurfaceTraction(3), Normal(3), &
        Velo(3), ut_eps, TanFrictionCoeff, DummyVals(1)
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t), SAVE :: ExtPressure_h, SurfaceTraction_h, SlipCoeff_h, &
        WeertmanCoeff_h, WeertmanExp_h, &
        FrictionNewtonEps_h, FrictionUt0_h, FrictionNewton_h, FrictionCoeff_h
    TYPE(VariableHandle_t), SAVE :: Velo_v
    TYPE(Variable_t), POINTER, SAVE :: NrmSol, HeightVar
    TYPE(ValueList_t), POINTER :: BC    
    REAL(KIND=dp) :: TanFder,JAC(nd*2,nd*2),SOL(nd*2),NodalSol(2,nd)

    SAVE Basis, dBasisdx

    !------------------------------------------------------------------------------

    IF( InitHandles ) THEN   
      CALL ListInitElementKeyword( ExtPressure_h,'Boundary Condition','Normal Surface Traction')
      IF( .NOT. ListGetElementSomewhere( ExtPressure_h) ) THEN
        CALL ListInitElementKeyword( ExtPressure_h,'Boundary Condition','External Pressure')      
      END IF
      CALL ListInitElementKeyword( SurfaceTraction_h,'Boundary Condition','Surface Traction',InitVec3D=.TRUE.)
      CALL ListInitElementKeyword( SlipCoeff_h,'Boundary Condition','Slip Coefficient',InitVec3D=.TRUE.)

      CALL ListInitElementKeyword( FrictionCoeff_h,'Boundary Condition','Friction Coefficient',&
          EvaluateAtIp=.TRUE., DummyCount=1)     
      CALL ListInitElementKeyword( FrictionNewtonEps_h,'Boundary Condition','Friction Newton Epsilon')     
      CALL ListInitElementKeyword( FrictionNewton_h,'Boundary Condition','Friction Newton Linearization')
      CALL ListInitElementKeyword( FrictionUt0_h,'Boundary Condition','Friction Linear Velocity')

      CALL ListInitElementKeyword( WeertmanCoeff_h,'Boundary Condition','Weertman Friction Coefficient')
      CALL ListInitElementKeyword( WeertmanExp_h,'Boundary Condition','Weertman Exponent')

      str = ListGetString( CurrentModel % Solver % Values,'Normal Vector Name',Found )
      IF(.NOT. Found) str = 'Normal Vector'
      NrmSol => VariableGet( CurrentModel % Solver % Mesh % Variables, str, ThisOnly = .TRUE.) 

      CALL ListInitElementVariable( Velo_v )

      HeightVar => VariableGet( CurrentModel % Mesh % Variables,'height') 
      
      InitHandles = .FALSE.
    END IF

    IF( ALLOCATED( Basis ) ) THEN
      IF( SIZE( Basis ) < nd ) THEN
        DEALLOCATE( Basis ) 
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

    FrictionNewton = .FALSE.

    ! There is no elemental routine for this.
    ! So whereas this breaks the beauty it does not cost too much.
    BC => GetBC()     
    HaveFrictionW = ListCheckPresent( BC,'Weertman Friction Coefficient') 
    HaveFrictionU = ListCheckPresent( BC,'Friction Coefficient')
    HaveFriction = HaveFrictionU .OR. HaveFrictionW

    IF( HaveFriction ) THEN
      FrictionNewton = ListGetElementLogical(FrictionNewton_h, Element ) 
      IF( FrictionNewton ) THEN
        ut_eps = ListGetElementReal( FrictionNewtonEps_h, Element = Element, Found = Found )
        IF(.NOT. Found ) ut_eps = 1.0e-8
      END IF
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
      !------------------------------------
      ExtPressure = ListGetElementReal( ExtPressure_h, Basis, Element, HavePres, GaussPoint = t )      

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
      IF( HaveFriction ) THEN
        IF( HaveSlip ) THEN
          CALL Fatal('IncompressibleNSVec','You cannot combine different friction models!')
        END IF

        ! Velocity at integration point for nonlinear friction laws
        Velo = ListGetElementVectorSolution( Velo_v, Basis, Element, dofs = dofs )

        ! For hyrdostatic model it is assumed that friction is walsy in (x,y) plane.
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

              IF(HaveFrictionW .AND. Newton) THEN
                DO j=1,dofs
                  JAC((p-1)*c+i,(q-1)*c+j ) = &
                      JAC((p-1)*c+i,(q-1)*c+j ) + &
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

    IF(HaveFrictionW .AND. Newton) THEN
      CALL GetLocalSolution( NodalSol )
      SOL=0._dp
      DO i = 1, c
        SOL(i::c) = NodalSol(i,1:nd)
      END DO
      STIFF=STIFF+JAC
      FORCE=FORCE + MATMUL(JAC,SOL)
    END IF

    CALL DefaultUpdateEquations( STIFF, FORCE )

  END SUBROUTINE LocalBoundaryMatrix


!------------------------------------------------------------------------------
! Compute vector for computing du_z/dz and corresponding weight.
!------------------------------------------------------------------------------
  SUBROUTINE LocalDuz(Element, n, ntot, duz, wuz )
!------------------------------------------------------------------------------
    USE LinearForms
    IMPLICIT NONE

    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    INTEGER, INTENT(IN) :: n, ntot
    REAL(KIND=dp), POINTER :: duz(:), wuz(:)
!------------------------------------------------------------------------------
    INTEGER :: dim, dofs
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    LOGICAL :: Stat, Found
    INTEGER :: NodalPerm(ntot)
    REAL(KIND=dp) :: NodalVelo(2,ntot), duz_elem(ntot), wuz_elem(ntot), vgrad(2), w
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

    duz_elem = 0.0_dp
    wuz_elem = 0.0_dp
    
    DO t = 1, ngp
      DO i = 1, 2
        vgrad(i) = SUM(dBasisdxVec(t,1:ntot,i)*NodalVelo(i,1:ntot))
      END DO
      DO i=1,ntot
        ! It would seem that weighting with the DetJ would make sense, but maybe not...
        !w = DetJVec(t) * BasisVec(t,i)
        w = BasisVec(t,i) * IP % s(t)

        ! Negative sign comes from continuity equation: du_z/dz = -(du_x/dx + du_y/dy)
        duz_elem(i) = duz_elem(i) - w * SUM(vgrad(1:2))
        wuz_elem(i) = wuz_elem(i) + w
      END DO
    END DO

    ! Global degrees of freedom
    NodalPerm(1:ntot) = CurrentModel % Solver % Variable % Perm(Element % NodeIndexes )

    ! Sum up the contribution
    duz(NodalPerm(1:ntot)) = duz(NodalPerm(1:ntot)) + duz_elem(1:ntot)
    wuz(NodalPerm(1:ntot)) = wuz(NodalPerm(1:ntot)) + wuz_elem(1:ntot)
    
  END SUBROUTINE LocalDuz

!------------------------------------------------------------------------------

  
  SUBROUTINE PopulateDerivedFields()

    IMPLICIT NONE

    TYPE(Variable_t), POINTER :: VarXY, VarFull, VarDuz, VarP 
    CHARACTER(LEN=MAX_NAME_LEN):: str
    TYPE(ValueList_t), POINTER :: Params, Material
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: Found
    REAL(KIND=dp), POINTER :: duz(:), wuz(:)
    INTEGER, POINTER, SAVE :: UpPointer(:), DownPointer(:)
    INTEGER :: i,j,k,i1,i2,k1,k2,t,n,nd,nb,active, dofs, pdof
    TYPE(Element_t), POINTER :: Element
    TYPE(Solver_t), POINTER :: pSolver
    REAL(KIND=dp) :: dz, rho, g
    REAL(KIND=dp), POINTER :: gWork(:,:)
    
    
    SAVE :: duz, wuz

    pSolver => CurrentModel % Solver    
    Params => pSolver % Values
    VarXY => pSolver % Variable

    str = ListGetString( Params,'Velocity Vector Name',Found )
    IF(.NOT. Found) str = ListGetString( Params,'Velocity Variable Name',Found )
    IF(.NOT. Found) RETURN

    Mesh => CurrentModel % Mesh
    VarFull => VariableGet(Mesh % Variables, str, ThisOnly = .TRUE.)
    IF(.NOT. ASSOCIATED(VarFull)) THEN
      CALL Fatal('HydrostaticNSVec','Could not find full velocity variable: '//TRIM(str))
    END IF
    dofs = VarFull % dofs

    NULLIFY(VarP) 
    str = ListGetString( Params,'Pressure Variable Name',Found )
    IF( Found ) THEN
      VarP => VariableGet(Mesh % Variables, str, ThisOnly = .TRUE.)
      IF(.NOT. ASSOCIATED(VarP)) THEN
        CALL Fatal('HydrostaticNSVec','Could not find full pressure variable: '//TRIM(str))
      END IF
    END IF
    IF(.NOT. ASSOCIATED(VarP)) THEN
      IF(dofs == 4) VarP => VarFull
    END IF
        
    ! Pressure can be 1st component of pressure or last compenent of 'flow solution'
    pdof = 0
    IF(ASSOCIATED(VarP)) pdof = VarP % Dofs
    g = 1.0_dp
    IF(pdof > 0) THEN
      gWork => ListGetConstRealArray( CurrentModel % Constants,'Gravity',Found)
      IF(Found) THEN
        g = ABS(gWork(SIZE(gWork,1),1))
      ELSE
        CALL Warn('HydrostaticNSVec','"Gravity" not found in simulation section, setting to 1')
      END IF
    END IF       
    
    ! Copy the velocity componts x & y
    DO i=1,2
      VarFull % Values(i::dofs) = VarXY % Values(i::2)
    END DO

    n = SIZE(VarXY % Values) / 2
    ALLOCATE(duz(n),wuz(n))      
    duz = 0.0_dp; wuz = 0.0_dp
        
    Active = GetNOFActive()
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes(Element)
      nb = GetElementNOFBDOFs(Element)
      nd = GetElementNOFDOFs(Element)
      
      ! Calculate elemental contribution to d(uz)/dz
      !---------------------------------------------
      CALL LocalDuz(Element, n, n, duz, wuz )      

      IF(t==1) THEN
        Material => GetMaterial(Element)
        rho = ListGetCReal(Material,'Density',UnfoundFatal=.TRUE.) 
      END IF
    END DO

    ! Note: we are missing parallel communication here!
    WHERE(wuz > EPSILON(dz) )
      duz = duz / wuz
    END WHERE

    ! If we have exported variable "duz" we can save the value of the field
    VarDuz => VariableGet(Mesh % Variables,'duz', ThisOnly = .TRUE.)
    IF( ASSOCIATED( VarDuz ) ) VarDuz % Values = duz
    
    ! We need to save pointers 
    CALL DetectExtrudedStructure( Mesh, PSolver, &
        UpNodePointer = UpPointer, DownNodePointer = DownPointer)

    ! Integrate over structured mesh 
    DO i=1,Mesh % NumberOfNodes                   
      IF(DownPointer(i)==0) CYCLE
      i1 = i      
            
      ! Nodes at the bedrock
      IF(DownPointer(i) == i) THEN
        ! Initailize values at the bedrock
        k1 = VarFull % Perm(i1)
        IF(k1>0) VarFull % Values(dofs*(k1-1)+3) = 0.0_dp

        DO WHILE(.TRUE.)
          i2 = UpPointer(i1)
          IF(i2==i1) EXIT

          k1 = VarFull % Perm(i1)
          k2 = VarFull % Perm(i2)

          IF(k1>0 .AND. k2>0) THEN
            ! Integrate for the z-velocity
            dz = Mesh % Nodes % z(i2) - Mesh % Nodes % z(i1)
            VarFull % Values(dofs*(k2-1)+3) = VarFull % Values(dofs*(k1-1)+3) + &
                dz * ( duz(k1) + duz(k2) ) / 2.0_dp
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
    
  END SUBROUTINE PopulateDerivedFields
    
END MODULE HydrostaticNSUtils


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
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER :: Params 
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(GaussIntegrationPoints_t) :: IP
  INTEGER :: Element_id
  INTEGER :: i, n, nb, nd, dim, Active, maxiter, iter
  REAL(KIND=dp) :: Norm
  LOGICAL :: Found, Converged
  LOGICAL :: SpecificLoad, InitBCHandles
  CHARACTER(*), PARAMETER :: Caller = 'HydrostaticNSSolver'


!------------------------------------------------------------------------------
! Local variables to be accessed by the contained subroutines:
!------------------------------------------------------------------------------
  LOGICAL :: LinearAssembly, Newton
!------------------------------------------------------------------------------ 

  CALL DefaultStart()

  dim = CoordinateSystemDimension()
  Mesh => GetMesh()
  Params => GetSolverParams() 

  !-----------------------------------------------------------------------------
  ! Output the number of integration points as information.
  ! This in not fully informative if several element types are present.
  !-----------------------------------------------------------------------------
  Element => Mesh % Elements( Solver % ActiveElements(1) ) 
  IP = GaussPointsAdapt( Element, PReferenceElement = .TRUE. )
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

  ! Compute all the other fields that can be computed from the solution.
  CALL PopulateDerivedFields()
  
  CALL Info( Caller,'All done',Level=10)
!------------------------------------------------------------------------------
END SUBROUTINE HydrostaticNSSolver
!------------------------------------------------------------------------------
