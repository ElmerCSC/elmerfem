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
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/

!> \ingroup Solvers

!------------------------------------------------------------------------------
SUBROUTINE MagnetoDynamicsCalcFields_Init0(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE MagnetoDynamicsUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  CHARACTER(LEN=MAX_NAME_LEN) :: sname,pname
  LOGICAL :: Found, ElementalFields, RealField, FoundVar, Hcurl
  INTEGER, POINTER :: Active(:)
  INTEGER :: mysolver,i,j,k,l,n,m,vDOFs, soln,pIndex
  TYPE(ValueList_t), POINTER :: SolverParams, DGSolverParams
  TYPE(Solver_t), POINTER :: Solvers(:)

  ! This is really using DG so we don't need to make any dirty tricks to create DG fields
  ! as is done in this initialization. 
  SolverParams => GetSolverParams()
  
  ! The only purpose of this parsing of the variable name is to identify
  ! whether the field is real or complex. As the variable has not been
  ! created at this stage we have to do some dirty parsing. 
  pname = GetString(SolverParams, 'Potential variable', Found)
  vdofs = 0
  pIndex = 0 
  FoundVar = .FALSE.

  IF( Found ) THEN
    DO i=1,Model % NumberOfSolvers
      sname = GetString(Model % Solvers(i) % Values, 'Variable', Found)
      
      J=INDEX(sname,'[')-1
      IF ( j<=0 ) j=LEN_TRIM(sname)
      IF ( sname(1:j) == pname(1:LEN_TRIM(pname)) )THEN
        k = 0
        vDofs = 0
        j=INDEX(sname,':')
        DO WHILE(j>0)
          Vdofs=Vdofs+ICHAR(sname(j+k+1:j+k+1))-ICHAR('0')
          k = k+j
          IF(k<LEN(sname)) j=INDEX(sname(k+1:),':')
        END DO
        pIndex = i
        FoundVar = .TRUE.
        EXIT
      END IF
    END DO
    
    IF(.NOT. FoundVar ) THEN
      CALL Fatal('MagnetoDynamicsCalcFields_Init0','Could not find solver for variable: '//TRIM(sname))
    END IF
  END IF

  
  ! When we created the case for GUI where "av" is not given in sif then it is impossible to
  ! determine from the variable declaration what kind of solver we have. 
  IF( .NOT. FoundVar ) THEN
    Hcurl = .FALSE.
    DO i=Model % NumberOfSolvers,1,-1
      sname = GetString(Model % Solvers(i) % Values, 'Procedure', Found)
      
      j = INDEX( sname,'WhitneyAVHarmonicSolver')
      IF( j > 0 ) THEN
        Hcurl = .TRUE.
        vDofs = 2
        EXIT
      END IF

      j = INDEX( sname,'MagnetoDynamics2DHarmonic')
      IF( j > 0 ) THEN
        Vdofs = 2
        EXIT
      END IF

      j = INDEX( sname,'WhitneyAVSolver')
      IF( j > 0 ) THEN
        Hcurl = .TRUE.
        vDofs = 1
        EXIT
      END IF
      
      j = INDEX( sname,'MagnetoDynamics2D')
      IF( j > 0 ) THEN
        Vdofs = 1
        EXIT
      END IF
    END DO

    IF( Vdofs == 0 ) THEN
      CALL Fatal('MagnetoDynamicsCalcFields_Init0','Could not determine target variable type (real or complex)')
    END IF
    pIndex = i 
  END IF

  RealField = ( Vdofs /= 2 )
  IF( RealField ) THEN
    CALL Info('MagnetoDynamicsCalcFields_Init0','The target solver seems to be real valued',Level=12)
  ELSE
    CALL Info('MagnetoDynamicsCalcFields_Init0','The target solver seems to be complex valued',Level=12)
  END IF

  IF( GetLogical(Model % Solvers(pIndex) % Values, 'Eigen Analysis', Found) ) THEN
    CALL ListAddNewLogical( SolverParams,'Eigen Analysis',.TRUE.)
  END IF
  
  CALL ListAddNewLogical( SolverParams, 'Target Variable Real Field', RealField )   
  CALL Info('MagnetoDynamicsCalcFields_Init0','Target Variable Solver Index: '&
    //I2S(pIndex),Level=12)
  CALL ListAddNewInteger( SolverParams, 'Target Variable Solver Index', pIndex ) 
  
!------------------------------------------------------------------------------
END SUBROUTINE MagnetoDynamicsCalcFields_Init0
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE MagnetoDynamicsCalcFields_Init(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE MagnetoDynamicsUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------

  INTEGER  :: i, dim, fdim
  LOGICAL :: Found, FluxFound, NodalFields, ElementalFields, &
      RealField, ComplexField, LorentzConductivity, DoIt, VtuStyle
  TYPE(ValueList_t), POINTER :: EQ, SolverParams

  LorentzConductivity = ListCheckPrefixAnyBodyForce(Model, "Angular Velocity") .or. &
    ListCheckPrefixAnyBodyForce(Model, "Lorentz Velocity")

  IF(.NOT.ASSOCIATED(Solver % Values)) Solver % Values=>ListAllocate()
  SolverParams => GetSolverParams()

  dim = CoordinateSystemDimension()
  fdim = 3   
  IF( ListGetLogical( SolverParams,'2D result fields',Found ) ) THEN    
    IF( dim == 2 ) THEN
      CALL Info('MagnetoDynamicsCalcFields_init','Enforcing result fields to be 2D!',Level=12)
      fdim = dim
    ELSE
      CALL Info('MagnetoDynamicsCalcFields_init','Keyword "2D result fields" is not applicable in 3D!',Level=7)
    END IF
  END IF
      
  ! Inherit this from the _init0 solver. Hence we know it must exist!
  RealField = ListGetLogical( SolverParams,'Target Variable Real Field') 

  ! We have some challenges in plotting eigenvalues for vtu if the fields are of size
  ! 6 and have real and imaginary parts living separately. 
  VtuStyle = ListGetLogical( SolverParams,'Vtu Style', Found )
  IF(VtuStyle) RealField = .TRUE.

  ComplexField = .NOT. RealField  
   
  CALL ListAddString( SolverParams, 'Variable', '-nooutput hr_dummy' )

  !CALL ListAddNewLogical( SolverParams,'Skip Compute Nonlinear Change',.TRUE.)
  
  CALL ListAddLogical( SolverParams, 'Linear System refactorize', .FALSE.)

  ! add these in the beginning, so that SaveData sees these existing, even
  ! if executed before the actual computations...
  ! -----------------------------------------------------------------------
  CALL ListAddConstReal(Model % Simulation,'res: Eddy current power',0._dp)

  IF( ListGetLogical( SolverParams,'Separate Magnetic Energy',Found ) ) THEN
    CALL ListAddConstReal(Model % Simulation,'res: Electric Field Energy',0._dp)
    CALL ListAddConstReal(Model % Simulation,'res: Magnetic Field Energy',0._dp)
  ELSE
    CALL ListAddConstReal(Model % Simulation,'res: ElectroMagnetic Field Energy',0._dp)
  END IF
    
  IF (GetLogical(SolverParams,'Show Angular Frequency',Found)) &
    CALL ListAddConstReal(Model % Simulation,'res: Angular Frequency',0._dp)

  ! add these in the beginning only if the Magnetix Flux Average is computed
  ! -------------------------------------------------------------------------
  IF (ListGetLogicalAnyBC( Model,'Magnetic Flux Average')) THEN
    CALL ListAddConstReal( Model % Simulation,'res: Magnetic Flux Average', 0._dp)
    CALL ListAddConstReal(Model % Simulation, &
                           'res: Magnetic Flux Density Average',0._dp)

    IF (.NOT. RealField ) THEN 
      CALL ListAddConstReal(Model % Simulation,'res: Magnetic Flux im Average',0._dp)
      CALL ListAddConstReal(Model % Simulation, &
                   'res: Magnetic Flux Density im Average', 0._dp )
    END IF

    CALL ListAddConstReal(Model % Simulation,'res: Magnetic Flux Area',0._dp)
  END IF

  NodalFields = .NOT. GetLogical( SolverParams, 'Skip Nodal Fields', Found)
  IF(.NOT. Found ) NodalFields = GetLogical( SolverParams, 'Calculate Nodal Fields', Found)
  IF(.NOT. Found ) NodalFields = .TRUE.
  
  i=1
  DO WHILE(.TRUE.)
    IF ( .NOT. ListCheckPresent(SolverParams,"Exported Variable "//i2s(i)) ) EXIT
    i = i + 1
  END DO
  i = i - 1
  
  IF( NodalFields ) THEN
    DoIt = GetLogical(SolverParams,'Calculate Magnetic Flux Density',Found)
    IF(.NOT. Found) DoIt = .TRUE.

    IF( DoIt ) THEN
      i = i + 1

      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "Magnetic Flux Density[Magnetic Flux Density:"//I2S(fdim)//"]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "Magnetic Flux Density[Magnetic Flux Density re:"//I2S(fdim)//&
            " Magnetic Flux Density im:"//I2S(fdim)//"]" )
      END IF
    END IF
      
    IF (GetLogical(SolverParams,'Calculate Magnetic Vector Potential',Found)) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "Magnetic Vector Potential[Magnetic Vector Potential:3]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "Magnetic Vector Potential[Magnetic Vector Potential re:3 Magnetic Vector Potential im:3]")
      END IF
    END IF

    IF (GetLogical(SolverParams,'Calculate Magnetic Field Strength',Found)) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "Magnetic Field Strength[Magnetic Field Strength:"//I2S(fdim)//"]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "Magnetic Field Strength[Magnetic Field Strength re:"//I2S(fdim)//&
            " Magnetic Field Strength im:"//I2S(fdim)//"]")
      END IF
    END IF

    IF (GetLogical(SolverParams,'Calculate JxB',Found)) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "JxB[JxB:3]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "JxB[JxB re:3 JxB im:3]")
      END IF
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Maxwell Stress', Found ) ) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "Maxwell Stress[Maxwell Stress:6]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "Maxwell Stress[Maxwell Stress re:6 Maxwell Stress im:6]" )
      END IF
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Current Density', Found ) ) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "Current Density[Current Density:3]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "Current Density[Current Density re:3 Current Density im:3]" )
      END IF      
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Joule Heating', Found ) ) THEN
      i = i + 1
      CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
          "Joule Heating" )
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Harmonic Loss', Found ) ) THEN
      IF( RealField ) THEN
        CALL Warn('MagnetcDynamicsCalcFields',&
            'Harmonic loss computation only available for complex systems!')
      ELSE
        i = i + 1
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "Harmonic Loss Linear" )
        i = i + 1
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "Harmonic Loss Quadratic" )
      END IF
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Homogenization Loss', Found ) ) THEN
      IF( RealField ) THEN
        CALL Warn('MagnetoDynamicsCalcFields',&
            'Homogenization loss computation only available for complex systems!')
      ELSE
        i = i + 1
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "Proximity Loss" )
      END IF
    END IF

    IF ( Transient .OR. .NOT. RealField .OR. LorentzConductivity) THEN
      IF ( GetLogical( SolverParams, 'Calculate Electric Field', Found ) ) THEN
        i = i + 1
        IF ( RealField ) THEN
          CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
              "Electric Field[Electric Field:3]" )
        ELSE
          CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
              "Electric Field[Electric Field re:3 Electric Field im:3]" )
        END IF
      END IF

      IF ( GetLogical( SolverParams, 'Calculate Winding Voltage', Found ) ) THEN
        i = i + 1
        IF ( RealField ) THEN
          CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
              "Winding Voltage" )
        ELSE
          CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
              "Winding Voltage[Winding Voltage re:1 Winding Voltage im:1]" )
        END IF
      END IF
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Relative Permeability', Found ) ) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "Relative Permeability" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "Relative Permeability[Relative Permeability re:1 Relative Permeability im:1]" )
      END IF
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Electric Scalar Potential', Found ) ) THEN
      IF ( RealField ) THEN
        i = i + 1
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "Electric Scalar Potential" )
      END IF
    END IF
  END IF    

  IF ( GetLogical( SolverParams, 'Calculate Current Density', Found ) ) THEN
    IF( ListCheckPresentAnyBC( Model,'Layer Electric Conductivity') ) THEN 
      i = i + 1
      IF ( RealField ) THEN
        !CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
        !    "Surface Current[Surface Current:3]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "Surface Current[Surface Current re:3 Surface Current im:3]" )
      END IF
    END IF
  END IF

  IF ( GetLogical( SolverParams, 'Calculate Nodal Heating', Found ) ) THEN
    i = i + 1
    CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
        "Nodal Joule Heating" )
  END IF
    
  IF( ListGetLogicalAnyComponent(Model,'Calculate Magnetic Force') .OR. &
      ListGetLogicalAnyComponent(Model,'Calculate Magnetic Torque') .OR. &
      GetLogical(SolverParams, 'Calculate Nodal Forces', Found) ) THEN
    i = i + 1
    CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
        "Nodal Force[Nodal Force:"//I2S(fdim)//"]" )
  END IF
    
  ! If we have DG for the standard fields they are already elemental...
  IF (GetLogical(SolverParams,'Discontinuous Galerkin',Found)) RETURN

  ! Choose elemental if not otherwise specified. 
  ElementalFields = .NOT. GetLogical( SolverParams, 'Skip Elemental Fields', Found)
  IF(.NOT. Found ) ElementalFields = GetLogical( SolverParams, 'Calculate Elemental Fields', Found)
  IF(.NOT. Found ) ElementalFields = .TRUE.
  
  IF( ElementalFields ) THEN
    DoIt = GetLogical(SolverParams,'Calculate Magnetic Flux Density',Found)
    IF(.NOT. Found) DoIt = .TRUE.
    
    IF( DoIt ) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg Magnetic Flux Density E[Magnetic Flux Density E:"//I2S(fdim)//"]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg Magnetic Flux Density E[Magnetic Flux Density re E:"//I2S(fdim)//&
            " Magnetic Flux Density im E:"//I2S(fdim)//"]" )
      END IF
    END IF
      
    IF (GetLogical(SolverParams,'Calculate Magnetic Vector Potential',Found)) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg Magnetic Vector Potential E[Magnetic Vector Potential E:3]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg Magnetic Vector Potential E[Magnetic Vector Potential re E:3 Magnetic Vector Potential im E:3]" )
      END IF
    END IF

    IF (GetLogical(SolverParams,'Calculate Magnetic Field Strength',Found)) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg Magnetic Field Strength E[Magnetic Field Strength E:"//I2S(fdim)//"]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg Magnetic Field Strength E[Magnetic Field Strength re E:"//I2S(fdim)//&
            " Magnetic Field Strength im E:"//I2S(fdim)//"]" )
      END IF
    END IF

    IF (GetLogical(SolverParams,'Calculate JxB',Found)) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg JxB E[JxB E:3]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg JxB E[JxB re E:3 JxB im E:3]" )
      END IF
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Maxwell Stress', Found ) ) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg Maxwell Stress E[Maxwell Stress E:6]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg Maxwell Stress E[Maxwell Stress re E:6 Maxwell Stress im E:6]" )
      END IF
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Current Density', Found ) ) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg Current Density E[Current Density E:3]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg Current Density E[Current Density re E:3 Current Density im E:3]" )
      END IF
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Joule Heating', Found ) ) THEN
      i = i + 1
      CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
          "-dg Joule Heating E" )
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Harmonic Loss', Found ) ) THEN
      IF( RealField ) THEN
        CALL Warn('MagnetoDynamicsCalcFields',&
            'Harmonic loss computation only available for complex systems!')
      ELSE
        i = i + 1
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg Harmonic Loss Linear E" )
        i = i + 1
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg Harmonic Loss Quadratic E" )
      END IF
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Homogenization Loss', Found ) ) THEN
      IF( RealField ) THEN
        CALL Warn('MagnetoDynamicsCalcFields',&
            'Homogenization loss computation only available for complex systems!')
      ELSE
        i = i + 1
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg Proximity Loss E" )
      END IF
    END IF

    IF ( Transient .OR. ComplexField .OR. LorentzConductivity ) THEN
      IF ( GetLogical( SolverParams, 'Calculate Electric Field', Found ) ) THEN
        i = i + 1
        IF ( RealField ) THEN
          CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
              "-dg Electric Field E[Electric Field E:3]" )
        ELSE
          CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
              "-dg Electric Field E[Electric Field re E:3 Electric Field im E:3]" )
        END IF
      END IF

      IF ( GetLogical( SolverParams, 'Calculate Winding Voltage', Found ) ) THEN
        i = i + 1
        IF ( RealField ) THEN
          CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
              "-dg Winding Voltage E" )
        ELSE
          CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
              "-dg Winding Voltage E[Winding Voltage re E:1 Winding Voltage im E:1]" )
        END IF
      END IF

    END IF

    IF ( GetLogical( SolverParams, 'Calculate Relative Permeability', Found ) ) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg Relative Permeability E" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg Relative Permeability E[Relative Permeability re E:1 Relative Permeability im E:1]" )
      END IF
    END IF

    DoIt = ListGetLogicalAnyComponent(Model,'Calculate Magnetic Force') .OR. &
        ListGetLogicalAnyComponent(Model,'Calculate Magnetic Torque') .OR. &
        GetLogical( SolverParams,'Calculate Nodal Forces', Found)
    IF( DoIt ) THEN
      i = i + 1
      CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
          "-dg Nodal Force E[Nodal Force E:"//I2S(fdim)//"]" )
    END IF
  END IF
    
!------------------------------------------------------------------------------
END SUBROUTINE MagnetoDynamicsCalcFields_Init
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Calculate fields resulting from the edge element formulation of the magnetic 
!> field equations. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
 SUBROUTINE MagnetoDynamicsCalcFields(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
   USE MagnetoDynamicsUtils
   USE CircuitUtils
   USE Zirka
   use zirkautils
   
   IMPLICIT NONE
!------------------------------------------------------------------------------
   TYPE(Solver_t), TARGET :: Solver
   TYPE(Model_t) :: Model
   REAL(KIND=dp) :: dt
   LOGICAL :: Transient
!------------------------------------------------------------------------------
   REAL(KIND=dp), ALLOCATABLE :: WBasis(:,:), RotWBasis(:,:), Basis(:), lBasis(:), &
                dBasisdx(:,:)
   REAL(KIND=dp), ALLOCATABLE :: SOL(:,:), PSOL(:), ElPotSol(:,:), C(:)
   REAL(KIND=dp), ALLOCATABLE :: Wbase(:), alpha(:), NF_ip(:,:)
   REAL(KIND=dp), ALLOCATABLE :: omega_velo(:,:), lorentz_velo(:,:)
   COMPLEX(KIND=dp), ALLOCATABLE :: Magnetization(:,:), BodyForceCurrDens(:,:)
   COMPLEX(KIND=dp), ALLOCATABLE :: R_Z(:), PR(:)
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: s,Norm, Mult
   REAL(KIND=dp) :: B(2,3), E(2,3), JatIP(2,3), VP_ip(2,3), JXBatIP(2,3), CC_J(2,3), HdotB, LMSol(2)
   REAL(KIND=dp) :: ldetJ,detJ, C_ip, ST(3,3), Omega, ThinLinePower, Power, Energy(3), w_dens
   REAL(KIND=dp) :: localThickness
   REAL(KIND=dp) :: Freq, FreqPower(2), FieldPower(2), LossCoeff(2), ElemLoss(2), ValAtIP
   REAL(KIND=dp) :: ComponentLoss(2,2), rot_velo(3), angular_velo(3)
   REAL(KIND=dp) :: Coeff, TotalLoss(3), LumpedForce(3), localAlpha, localV(2), nofturns, coilthickness
   REAL(KIND=dp) :: Flux(2), AverageFluxDensity(2), Area, N_j, wvec(3)
   REAL(KIND=dp) :: R_ip, mu_r
   REAL(KIND=dp), SAVE :: mu0 = 1.2566370614359173e-6_dp

   COMPLEX(KIND=dp) :: MG_ip(3), BodyForceCurrDens_ip(3), PR_ip
   COMPLEX(KIND=dp) :: CST(3,3)
   COMPLEX(KIND=dp) :: CMat_ip(3,3)  
   COMPLEX(KIND=dp) :: imag_value, Zs
   COMPLEX(KIND=dp), ALLOCATABLE :: Tcoef(:,:,:), Nu_el(:,:,:)
   COMPLEX(KIND=dp), POINTER, SAVE :: Reluct_Z(:,:,:) => NULL()
   COMPLEX(KIND=dp) :: R_ip_Z, Nu(3,3)
   
   INTEGER, PARAMETER :: ind1(6) = [1,2,3,1,2,1]
   INTEGER, PARAMETER :: ind2(6) = [1,2,3,2,3,3]

   TYPE(Variable_t), POINTER :: Var, MFD, MFS, CD, SCD, EF, MST, ESP, &
                                JH, NJH, VP, FWP, MPerm, JXB, ML, ML2, &
                                LagrangeVar, NF, PL
                                
   TYPE(Variable_t), POINTER :: EL_MFD, EL_MFS, EL_CD, EL_EF, &
                                EL_MST, EL_JH, EL_VP, EL_FWP, EL_MPerm, EL_JXB, EL_ML, EL_ML2, &
                                EL_NF, EL_SL, EL_PL

   INTEGER :: Active,i,j,k,l,m,n,nd,np,p,q,DOFs,vDOFs,dim,BodyId,&
              VvarDofs,VvarId,IvarId,Reindex,Imindex,EdgeBasisDegree,eq_n, Indexes(100)

   TYPE(Solver_t), POINTER :: pSolver, ElPotSolver
   CHARACTER(LEN=MAX_NAME_LEN) :: Pname, CoilType, ElectricPotName, LossFile, CurrPathPotName, str

   TYPE(ValueList_t), POINTER :: Material, BC, BodyForce, BodyParams, SolverParams, PrevMaterial

   LOGICAL :: Found, FoundIm, FoundMagnetization, stat, LossEstimation, HomogenizationLoss, &
              CalcFluxLogical, CoilBody, PreComputedElectricPot, ImposeCircuitCurrent, &
              ItoJCoeffFound, ImposeBodyForceCurrent, HasVelocity, HasAngularVelocity, &
              HasLorenzVelocity, HaveAirGap, UseElementalNF, HasTensorReluctivity, &
              ImposeBodyForcePotential, JouleHeatingFromCurrent, HasZirka, DoAve, HomogenizationModel
   LOGICAL :: PiolaVersion, ElementalFields, NodalFields, RealField, pRef
   LOGICAL :: CSymmetry, HasHBCurve, LorentzConductivity, HasThinLines=.FALSE., NewMaterial
   
   TYPE(GaussIntegrationPoints_t) :: IP
   TYPE(Nodes_t), SAVE :: Nodes
   TYPE(Element_t), POINTER :: Element

   INTEGER, ALLOCATABLE :: Pivot(:), TorqueGroups(:)
   INTEGER, POINTER :: MasterBodies(:)

   REAL(KIND=dp), POINTER CONTIG :: Fsave(:)
   REAL(KIND=dp) :: Babs
   TYPE(Mesh_t), POINTER :: Mesh
   REAL(KIND=dp), ALLOCATABLE, TARGET :: Gforce(:,:), MASS(:,:), FORCE(:,:)
   REAL(KIND=dp), ALLOCATABLE :: BodyLoss(:,:), RotM(:,:,:), Torque(:)

   REAL(KIND=dp), ALLOCATABLE :: ThinLineCrossect(:),ThinLineCond(:),SheetThickness(:)

   REAL(KIND=DP), POINTER :: Cwrk(:,:,:)=>NULL(), Cwrk_im(:,:,:)=>NULL()

   REAL(KIND=dp) :: ItoJCoeff, CircuitCurrent, CircEqVoltageFactor
   TYPE(ValueList_t), POINTER :: CompParams
   REAL(KIND=dp) :: DetF, F(3,3), G(3,3), GT(3,3)
   REAL(KIND=dp), ALLOCATABLE :: EBasis(:,:), CurlEBasis(:,:) 

   REAL(KIND=dp) :: xcoord, grads_coeff, val
   REAL(KIND=dp) :: HarmPowerCoeff 
   REAL(KIND=dp) :: line_tangent(3)
   INTEGER :: IOUnit, pIndex
   REAL(KIND=dp) :: SaveNorm
   INTEGER :: NormIndex, fdim
   LOGICAL, SAVE :: ConstantMassMatrixInUse = .FALSE.
   LOGICAL :: Parallel, Erroneous
   LOGICAL :: CoilUseWvec, WvecInitHandle=.TRUE.
   CHARACTER(LEN=MAX_NAME_LEN) :: CoilWVecVarname
   TYPE(VariableHandle_t), SAVE :: Wvec_h
   INTEGER, POINTER, SAVE :: SetPerm(:) => NULL()
   LOGICAL :: LayerBC, CircuitDrivenBC
   REAL(KIND=dp) :: SurfPower
   INTEGER :: jh_k, ComponentId
   REAL, ALLOCATABLE :: SurfWeight(:)
   TYPE(ValueHandle_t), SAVE :: mu_h
   REAL(KIND=dp), POINTER :: muTensor(:,:)
   LOGICAL :: HasReluctivityFunction, HBIntegProblem, MaterialExponents
   REAL(KIND=dp) :: rdummy
   INTEGER :: mudim, ElementalMode, cdofs, LossN

   TYPE VariableArray_t
     TYPE(Variable_t), POINTER :: Field => Null()
   END TYPE VariableArray_t

   TYPE(VariableArray_t) :: NodalFieldPointers(32), ElementalFieldPointers(32)
   TYPE(Variable_t), POINTER :: FieldVariable
   LOGICAL :: EigenAnalysis, VtuStyle, OldLossKeywords
   INTEGER :: Field, FieldsToCompute, NOFEigen, MaxFields, NoSlices

!-------------------------------------------------------------------------------------------

   IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
   
   CALL Info('MagnetoDynamicsCalcFields','------------------------------',Level=6)
   CALL Info('MagnetoDynamicsCalcFields','Computing postprocessed fields',Level=5)
   
   SolverParams => GetSolverParams()

   Parallel = ( ParEnv % PEs > 1 )    
   NoSlices = MAX(1,ListGetInteger( Model % Simulation,'Number Of Slices', Found ) )   
   
   dim = CoordinateSystemDimension()
   fdim = 3   
   IF( ListGetLogical( SolverParams,'2D result fields',Found ) ) fdim = dim

   ElementalMode = ListGetInteger( SolverParams,'Calculate Elemental Mode',Found )
   
  ! This is a hack to be able to control the norm that is tested for
   NormIndex = GetInteger(SolverParams,'Show Norm Index',Found ) 
   SaveNorm = 0.0_dp
   
   IF (GetLogical(SolverParams, 'Calculate harmonic peak power', Found)) THEN
     HarmPowerCoeff = 1.0_dp
   ELSE
     HarmPowerCoeff = 0.5_dp
   END IF

   pIndex = ListGetInteger( SolverParams,'Target Variable Solver Index',UnfoundFatal=.TRUE.)
   pSolver => Model % Solvers(pIndex) 
   pname = getVarName(pSolver % Variable)

   CALL Info('MagnetoDynamicsCalcFields','Using potential variable: '//TRIM(Pname),Level=7)

   ! Inherit the solution basis from the primary solver
   vDOFs = pSolver % Variable % DOFs
   
   CALL EdgeElementStyle(pSolver % Values, PiolaVersion, BasisDegree = EdgeBasisDegree ) 
   IF (PiolaVersion) &
       CALL Info('MagnetoDynamicsCalcFields', &
       'Using Piola transformed finite elements',Level=7)
      
   ElectricPotName = GetString(SolverParams, 'Precomputed Electric Potential', PrecomputedElectricPot)
   IF (PrecomputedElectricPot) THEN
     DO i=1, Model % NumberOfSolvers
       ElPotSolver => Model % Solvers(i)
       IF (ElectricPotName==getVarName(ElPotSolver % Variable)) EXIT
     END DO
   END IF

   ! Do have impedance BCs? 
   LayerBC = (ListCheckPresentAnyBC(Model, 'Layer Electric Conductivity') .AND. vdofs==2) 
   jh_k = 0

   ! Do we have a real or complex valued primary field?
   RealField = ( vDofs == 1 ) 

   IF( ListCheckPresentAnyMaterial(Model,'Reluctivity Function') ) THEN
     IF(.NOT. RealField ) THEN
       CALL Fatal('MagnetoDynamicsCalcFields','Reluctivity Function not implemented in complex cases!')
     END IF
     CALL ListInitElementKeyword( mu_h,'Material','Reluctivity Function',&
         EvaluateAtIp=.TRUE.,DummyCount=3)
   END IF
   HbIntegProblem = .FALSE.

   
   LorentzConductivity = ListCheckPrefixAnyBodyForce(Model, "Angular Velocity") .or. &
       ListCheckPrefixAnyBodyForce(Model, "Lorentz Velocity")

   Mesh => GetMesh()

   LagrangeVar => NULL()
   str = LagrangeMultiplierName( pSolver )
   LagrangeVar => VariableGet( Mesh % Variables, str, ThisOnly = .TRUE.)

   MFD => VariableGet( Mesh % Variables, 'Magnetic Flux Density' )
   EL_MFD => VariableGet( Mesh % Variables, 'Magnetic Flux Density E' )

   MFS => VariableGet( Mesh % Variables, 'Magnetic Field Strength')
   EL_MFS => VariableGet( Mesh % Variables, 'Magnetic Field Strength E')

   VP => VariableGet( Mesh % Variables, 'Magnetic Vector Potential')
   EL_VP => VariableGet( Mesh % Variables, 'Magnetic Vector Potential E')

   ESP => VariableGet( Mesh % Variables, 'Electric Scalar Potential')
   
   IF( .NOT. PreComputedElectricPot ) THEN
     ImposeBodyForcePotential = GetLogical(SolverParams, 'Impose Body Force Potential', Found)
     IF (.NOT. Found) ImposeBodyForcePotential = &
         ListCheckPresentAnyBodyForce( Model,'Electric Potential')
   ELSE
     ImposeBodyForcePotential = .FALSE.
   END IF

   ImposeBodyForceCurrent = GetLogical(SolverParams, 'Impose Body Force Current', Found)
   IF (.NOT. Found) ImposeBodyForceCurrent = ListCheckPrefixAnyBodyForce( Model,'Current Density')

   ImposeCircuitCurrent = GetLogical(SolverParams, 'Impose Circuit Current', Found)
   CurrPathPotName = GetString(SolverParams, 'Circuit Current Path Potential Name', Found)
   IF (.NOT. Found) CurrPathPotName = 'W'

   EF  => NULL(); EL_EF => NULL(); 
   CD  => NULL(); EL_CD => NULL();
   JH  => NULL(); EL_JH => NULL();
   FWP => NULL(); EL_FWP => NULL();
   MPerm => NULL(); EL_MPerm => NULL();
   JXB => NULL(); EL_JXB => NULL();
   ML  => NULL(); EL_ML => NULL();
   ML2 => NULL(); EL_ML2 => NULL();
   PL => NULL(); EL_PL => NULL();
   NF => NULL(); EL_NF => NULL();
   NJH => NULL()
   
   IF ( Transient .OR. .NOT. RealField .OR. LorentzConductivity ) THEN
     EF => VariableGet( Mesh % Variables, 'Electric Field' )
     FWP => VariableGet( Mesh % Variables, 'Winding Voltage' )

     EL_EF => VariableGet( Mesh % Variables, 'Electric Field E' )
     EL_FWP => VariableGet( Mesh % Variables, 'Winding Voltage E' )
   END IF

   MPerm => VariableGet( Mesh % Variables, 'Relative Permeability' )
   EL_MPerm => VariableGet( Mesh % Variables, 'Relative Permeability E' )

   NF => VariableGet( Mesh % Variables, 'Nodal Force') 
   EL_NF => VariableGet( Mesh % Variables, 'Nodal Force E')

   CD => VariableGet( Mesh % Variables, 'Current Density' )
   EL_CD => VariableGet( Mesh % Variables, 'Current Density E' )

   JH => VariableGet( Mesh % Variables, 'Joule Heating' )
   EL_JH => VariableGet( Mesh % Variables, 'Joule Heating E' )

   NJH => VariableGet( Mesh % Variables, 'Nodal Joule Heating' )

   SCD => VariableGet( Mesh % Variables, 'Surface Current' )
   IF( ASSOCIATED( SCD ) ) THEN
     SCD % Values = 0.0_dp
     i = SIZE( SCD % Values ) / SCD % Dofs
     ALLOCATE( SurfWeight(i) )
     SurfWeight = 0.0_dp
   END IF
   
   IF(.NOT. RealField ) THEN
     ML => VariableGet( Mesh % Variables, 'Harmonic Loss Linear')
     EL_ML => VariableGet( Mesh % Variables, 'Harmonic Loss Linear E')
     ML2 => VariableGet( Mesh % Variables, 'Harmonic Loss Quadratic')
     EL_ML2 => VariableGet( Mesh % Variables, 'Harmonic Loss Quadratic E')
     PL => VariableGet( Mesh % Variables, 'Proximity Loss')
     EL_PL => VariableGet( Mesh % Variables, 'Proximity Loss E')
   END IF

   JXB => VariableGet( Mesh % Variables, 'JxB')
   EL_JXB => VariableGet( Mesh % Variables, 'JxB E')

   MST => variableGet( Mesh % Variables, 'Maxwell stress' )
   EL_MST => variableGet( Mesh % Variables, 'Maxwell stress E' )

   DOFs = 0
   IF ( ASSOCIATED(MFD) .OR. ASSOCIATED(EL_MFD) ) DOFs=DOFs+vDofs*fdim
   IF ( ASSOCIATED(MFS) .OR. ASSOCIATED(EL_MFS) ) DOFs=DOFs+vDofs*fdim
   IF ( ASSOCIATED(VP) .OR. ASSOCIATED(EL_VP) ) DOFs=DOFs+vDofs*3
   IF ( ASSOCIATED(EF) .OR. ASSOCIATED(EL_EF) ) DOFs=DOFs+vDofs*3
   IF ( ASSOCIATED(CD) .OR. ASSOCIATED(EL_CD) ) DOFs=DOFs+vDofs*3
   IF ( ASSOCIATED(JXB) .OR. ASSOCIATED(EL_JXB) ) DOFs=DOFs+vDofs*3
   IF ( ASSOCIATED(FWP) .OR. ASSOCIATED(EL_FWP) ) DOFs=DOFs+vDofs*1
   IF ( ASSOCIATED(MPerm) .OR. ASSOCIATED(EL_MPerm) ) DOFs=DOFs+vDofs*1
   IF ( ASSOCIATED(MST) .OR. ASSOCIATED(EL_MST) ) DOFs=DOFs+vDofs*6

   ! These have just one component even in complex equation
   IF ( ASSOCIATED(JH) .OR. ASSOCIATED(EL_JH) .OR. ASSOCIATED(NJH)) DOFs=DOFs+1   
   IF ( ASSOCIATED(ML) .OR. ASSOCIATED(EL_ML) ) DOFs=DOFs+1
   IF ( ASSOCIATED(ML2) .OR. ASSOCIATED(EL_ML2) ) DOFs=DOFs+1   
   IF ( ASSOCIATED(PL) .OR. ASSOCIATED(EL_PL) ) DOFs=DOFs+1   
   IF ( ASSOCIATED(NF) .OR. ASSOCIATED(EL_NF) ) DOFs=DOFs+fdim


   CALL Info('MagnetoDynamicsCalcFields',&
       'Number of components to compute: '//I2S(DOFs),Level=8)

   MaxFields = 15  !  

   NodalFieldPointers(1) % Field => MFD
   NodalFieldPointers(2) % Field => MFS
   NodalFieldPointers(3) % Field => VP
   NodalFieldPointers(4) % Field => CD
   NodalFieldPointers(5) % Field => FWP
   NodalFieldPointers(6) % Field => MPerm
   NodalFieldPointers(7) % Field => EF
   NodalFieldPointers(8) % Field => JXB
   NodalFieldPointers(9) % Field => MST
   NodalFieldPointers(10) % Field => NF
   NodalFieldPointers(11) % Field => NJH
   NodalFieldPointers(12) % Field => JH
   NodalFieldPointers(13) % Field => ML
   NodalFieldPointers(14) % Field => ML2
   NodalFieldPointers(15) % Field => PL
   
   ElementalFieldPointers(1) % Field => EL_MFD
   ElementalFieldPointers(2) % Field => EL_MFS
   ElementalFieldPointers(3) % Field => EL_VP
   ElementalFieldPointers(4) % Field => EL_CD
   ElementalFieldPointers(5) % Field => EL_FWP
   ElementalFieldPointers(6) % Field => EL_MPerm
   ElementalFieldPointers(7) % Field => EL_EF
   ElementalFieldPointers(8) % Field => EL_JXB
   ElementalFieldPointers(9) % Field => EL_MST
   ElementalFieldPointers(10) % Field => EL_NF
   ElementalFieldPointers(11) % Field => EL_JH
   ElementalFieldPointers(12) % Field => EL_ML
   ElementalFieldPointers(13) % Field => EL_ML2
   ElementalFieldPointers(14) % Field => EL_PL
   ElementalFieldPointers(15) % Field => Null()

   NodalFields = .FALSE.; ElementalFields = .FALSE.
   DO i=1,MaxFields
     NodalFields = NodalFields .OR. ASSOCIATED(NodalFieldPointers(i) % Field)
     ElementalFields = ElementalFields .OR. ASSOCIATED(ElementalFieldPointers(i) % Field)
   END DO
   
   IF(NodalFields ) THEN
     ALLOCATE(GForce(SIZE(Solver % Matrix % RHS),DOFs)); Gforce=0._dp
   END IF
   
   n = Mesh % MaxElementDOFs
   ALLOCATE( MASS(2*n,2*n), FORCE(n,DOFs), Tcoef(3,3,n), RotM(3,3,n), Pivot(n))

!------------------------------------------------------------------------------
   ALLOCATE( WBasis(n,3), RotWBasis(n,3), Basis(n), dBasisdx(n,3), lBasis(n) )
   ALLOCATE( SOL(2,n), PSOL(n), ElPotSol(1,n), C(n) )
   ALLOCATE( Wbase(n), alpha(n), NF_ip(n,3) )
   ALLOCATE( PR(n), omega_velo(3,n), lorentz_velo(3,n) )
   ALLOCATE( Magnetization(3,n), BodyForceCurrDens(3,n), R_Z(n) )
!------------------------------------------------------------------------------
   SOL = 0._dp; PSOL=0._dp

   IF ( ASSOCIATED(ESP) ) THEN
     ESP  % Values(ESP % Perm) = pSolver % Variable % Values( &
         pSolver % Variable % Perm(1:Mesh % NumberOfNodes))
   END IF

   LossEstimation = GetLogical(SolverParams,'Loss Estimation',Found) &
       .OR. ASSOCIATED( ML ) .OR. ASSOCIATED( EL_ML ) &
       .OR. ASSOCIATED( ML2 ) .OR. ASSOCIATED( EL_ML2 ) 

   IF (LossEstimation) THEN
     OldLossKeywords = ListCheckPrefixAnyMaterial( Model,'Harmonic Loss Linear Frequency Exponent') 
     MaterialExponents = ListCheckPrefixAnyMaterial( Model,'Harmonic Loss Frequency Exponent') 
     MaterialExponents = MaterialExponents .OR. OldLossKeywords

     ! Fixed for now. FourierLoss solver more generic here. 
     LossN = 2 
          
     IF(.NOT. MaterialExponents) THEN
       OldLossKeywords = .NOT. ListCheckPresent(SolverParams,'Harmonic Loss Frequency Exponent')
       CALL GetLossExponents(SolverParams,FreqPower,FieldPower,LossN,OldLossKeywords)
     END IF

     IF( OldLossKeywords ) THEN
       IF(.NOT. ListCheckPresentAnyMaterial( Model,'Harmonic Loss Linear Coefficient') ) THEN
         CALL Warn('MagnetoDynamicsCalcFields',&
             'Harmonic loss requires > Harmonic Loss Linear Coefficient < in material section!')
       END IF

       IF(.NOT. ListCheckPresentAnyMaterial( Model,'Harmonic Loss Quadratic Coefficient') ) THEN
         CALL Warn('MagnetoDynamicsCalcFields',&
             'Harmonic loss requires > Harmonic Loss Quadratic Coefficient < in material section!')
       END IF       

       CALL Info('MagnetoDynamicsCalcFields','Consider using more generic keywords for loss computation!')
     END IF
     
     ComponentLoss = 0.0_dp
     ALLOCATE( BodyLoss(3,Model % NumberOfBodies) )
     BodyLoss = 0.0_dp
     TotalLoss = 0._dp
   END IF

   HomogenizationLoss = ASSOCIATED(PL) .OR. ASSOCIATED(EL_PL)
   IF (HomogenizationLoss) ALLOCATE( Nu_el(3,3,n) )

   VtuStyle = .FALSE.
   cdofs = 1
   EigenAnalysis = GetLogical( pSolver % Values, 'Eigen Analysis', Found )
   IF(EigenAnalysis) THEN
     CALL Info('MagnetoDynamicsCalcFields','Ensure space for eigen analysis',Level=10)
     VtuStyle = ListGetLogical( SolverParams,'Vtu Style', Found )
     IF(VtuStyle) cdofs=2     
     NOFeigen = SIZE(pSolver % Variable % EigenValues)
     DO i=1,MaxFields

       FieldVariable => NodalFieldPointers(i) % Field
       IF(ASSOCIATED(FieldVariable)) THEN

         ALLOCATE( FieldVariable % EigenValues(NOFEigen) )
         FieldVariable % EigenValues = pSolver % Variable % EigenValues

         n = cdofs*SIZE(FieldVariable % Values)/2         
         ALLOCATE( FieldVariable % EigenVectors(NOFEigen,n) )
         FieldVariable % EigenVectors = 0
       END IF

       FieldVariable => ElementalFieldPointers(i) % Field
       IF(ASSOCIATED(FieldVariable)) THEN

         ALLOCATE( FieldVariable % EigenValues(NOFEigen) )
         FieldVariable % EigenValues = pSolver % Variable % EigenValues

         n = cdofs*SIZE(FieldVariable % Values)/2
         ALLOCATE( FieldVariable % EigenVectors(NOFEigen,n) )
         FieldVariable % EigenVectors = 0
       END IF
     END DO
     FieldsToCompute = NofEigen
   ELSE
     FieldsToCompute = 1
   END IF

   COMPUTE_FIELDS: DO Field=1,FieldsToCompute

   IF ( NodalFields ) GForce = 0._dp

   IF(EigenAnalysis) THEN
     DO i=1,pSolver % Matrix % NumberOfRows/2
       pSolver % Variable % Values(2*i-1) = REAL(pSolver % Variable % EigenVectors(Field,i))
       pSolver % Variable % Values(2*i) = AIMAG(pSolver % Variable % EigenVectors(Field,i))
     END DO
   END IF

   C = 0._dp; PR=0._dp
   Magnetization = 0._dp

   Power = 0._dp; Energy = 0._dp
   IF(.NOT. ConstantMassMatrixInUse ) THEN
     CALL DefaultInitialize()
   END IF
   
   PrevMaterial => NULL()
   
   DO i = 1, GetNOFActive()
     Element => GetActiveElement(i)

     n = GetElementNOFNodes()
     IF(dim==2) THEN
       eq_n = GetElementDOFs(Indexes)
     ELSE
       eq_n = n
       Indexes(1:n) = Element % NodeIndexes
     END IF

     np = n*pSolver % Def_Dofs(GetElementFamily(Element),Element % BodyId,1)
     nd = GetElementNOFDOFs(uSolver=pSolver)

     IF (SIZE(Tcoef,3) /= n) THEN
       DEALLOCATE(Tcoef)
       ALLOCATE(Tcoef(3,3,n))
     END IF
     
     CALL GetElementNodes( Nodes, USolver=pSolver )

     ! If potential is not available we have to use given current directly to estimate Joule losses
     JouleHeatingFromCurrent = ( np == 0 .AND. &
         .NOT. ( PreComputedElectricPot .OR. ImposeBodyForcePotential ) )
     
     BodyId = GetBody()
     Material => GetMaterial()
     BodyForce => GetBodyForce()
     
     NewMaterial = .NOT. ASSOCIATED(Material, PrevMaterial)
     IF (NewMaterial) THEN
       HasHBCurve = ListCheckPresent(Material, 'H-B Curve')
       HasReluctivityFunction = ListCheckPresent(Material,'Reluctivity Function')
       PrevMaterial => Material
     END IF
     
     ItoJCoeffFound = .FALSE.
     IF(ImposeCircuitCurrent) THEN
       ItoJCoeff = ListGetConstReal(GetBodyParams(Element), &
         'Current to density coefficient', ItoJCoeffFound)
       IF(ItoJCoeffFound) THEN 
         CALL GetLocalSolution(Wbase,CurrPathPotName)
         IvarId = GetInteger(GetBodyParams(Element), 'Circuit Current Variable Id', Found)
         IF (.NOT. Found) CALL Fatal ('MagnetoDynamicsCalcFields', 'Circuit Current Variable Id not found!')
         CircuitCurrent = LagrangeVar % Values(IvarId)
       END IF  
     END IF


     CALL GetVectorLocalSolution(SOL,Pname,uSolver=pSolver)
     IF (PrecomputedElectricPot) THEN
       CALL GetScalarLocalSolution(ElPotSol(1,:),ElectricPotName,uSolver=ElPotSolver)
     END IF
     
     IF( ImposeBodyForcePotential ) THEN
       ElPotSol(1,:) = GetReal(BodyForce,'Electric Potential',Found)
     END IF
       
     IF ( Transient ) THEN
       CALL GetLocalSolution(PSOL,Pname,uSolver=pSolver,Tstep=-1)
       IF (pSolver % TimeOrder==1) THEN
         PSOL(1:nd) = (SOL(1,1:nd)-PSOL(1:nd)) / dt
       END IF
     END IF

     Omega = GetAngularFrequency(pSOlver % Values,Found,Element)
     IF( .NOT. ( RealField .OR. Found ) ) THEN
!      CALL Fatal('MagnetoDynamicsCalcFields',&
!          '(Angular) Frequency must be given for complex fields!')
     END IF
     Freq = Omega / (2*PI)
     
     IF ( ASSOCIATED(MFS) .OR. ASSOCIATED(EL_MFS)) THEN
       FoundMagnetization = .FALSE.
       IF(ASSOCIATED(BodyForce)) THEN
         CALL GetComplexVector( BodyForce,Magnetization(1:3,1:n),'Magnetization',FoundMagnetization)
       END IF

       IF(.NOT.FoundMagnetization .AND. ASSOCIATED(Material)) THEN
         CALL GetComplexVector( Material,Magnetization(1:3,1:n),'Magnetization',FoundMagnetization)
       END IF
     END IF

     IF (ImposeBodyForceCurrent ) THEN
       BodyForceCurrDens = 0._dp
       IF ( ASSOCIATED(BodyForce) ) THEN
         SELECT CASE(DIM)
         CASE(3)
           CALL GetComplexVector( BodyForce, BodyForceCurrDens(1:3,1:n), 'Current Density', Found )
         CASE(2)
           IF (Vdofs == 2) THEN
             BodyForceCurrDens(3,1:n) = CMPLX( ListGetReal(BodyForce, 'Current Density', n, Element % NodeIndexes, Found), &
               ListGetReal(BodyForce, 'Current Density im', n, Element % NodeIndexes, Found), KIND=dp)
           ELSE
             BodyForceCurrDens(3,1:n) = CMPLX( ListGetReal(BodyForce, 'Current Density', & 
               n, Element % NodeIndexes, Found), 0, KIND=dp)
           END IF
         END SELECT
         
         Mult = ListGetCReal( BodyForce,'Current Density Multiplier', Found )
         IF(Found) BodyForceCurrDens(1:3,1:n) = Mult * BodyForceCurrDens(1:3,1:n)
       END IF
     END IF

     CALL GetPermittivity(Material,PR,n)

     CoilBody = .FALSE.
     CompParams => GetComponentParams( Element )
     CoilType = ''
     RotM = 0._dp
     IF (ASSOCIATED(CompParams)) THEN
       CoilType = GetString(CompParams, 'Coil Type', Found)
       IF (Found) CoilBody = .TRUE.
       CircEqVoltageFactor = GetConstReal(CompParams, 'Circuit Equation Voltage Factor', Found)
       IF (.NOT. Found) CircEqVoltageFactor = 1._dp
     END IF 
 
     !------------------------------------------------------------------------------
     !  Read conductivity values (might be a tensor)
     !------------------------------------------------------------------------------
     !C(1:n) = GetReal(Material,'Electric Conductivity',Found)
     Tcoef = GetCMPLXElectricConductivityTensor(Element, n, CoilBody, CoilType) 

     dim = CoordinateSystemDimension()
     CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
         CurrentCoordinateSystem() == CylindricSymmetric )
     
     IF (CoilBody) THEN

       CoilUseWvec = GetLogical(CompParams, 'Coil Use W Vector', Found)
       IF (.NOT. Found) CoilUseWvec = .FALSE.
     
       IF (CoilUseWvec) THEN
         IF( WvecInitHandle ) THEN
           CoilWVecVarname = GetString(CompParams, 'W Vector Variable Name', Found)
           IF ( .NOT. Found) CoilWVecVarname = 'W Vector E'
           CALL ListInitElementVariable( Wvec_h, CoilWVecVarname )
           WvecInitHandle = .FALSE.
         END IF
       ELSE
         Call GetWPotential(Wbase)
       END IF
  
       SELECT CASE (CoilType)
       CASE ('stranded')
         IvarId = GetInteger (CompParams, 'Circuit Current Variable Id', Found)
         IF (.NOT. Found) CALL Fatal ('MagnetoDynamicsCalcFields', 'Circuit Current Variable Id not found!')

         N_j = GetConstReal (CompParams, 'Stranded Coil N_j', Found)
         IF (.NOT. Found) CALL Fatal ('MagnetoDynamicsCalcFields', 'Stranded Coil N_j not found!')


         HomogenizationModel = GetLogical(CompParams, 'Homogenization Model', Found)

         IF (HomogenizationLoss .and. HomogenizationModel) THEN
           BLOCK
             REAL(KIND=dp) :: nu_11(n), nuim_11(n), &
                              nu_22(n), nuim_22(n), &
                              nu_33(n), nuim_33(n)
             REAL(KIND=dp) :: sigma_33(n), sigmaim_33(n)

             nu_11 = 0._dp
             nuim_11 = 0._dp
             nu_11 = GetReal(CompParams, 'nu 11', Found)
             nuim_11 = GetReal(CompParams, 'nu 11 im', FoundIm)
             IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal ('LocalMatrix', 'Homogenization Model nu 11 not found!')
             nu_22 = 0._dp
             nuim_22 = 0._dp
             nu_22 = GetReal(CompParams, 'nu 22', Found)
             nuim_22 = GetReal(CompParams, 'nu 22 im', FoundIm)
             IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal ('LocalMatrix', 'Homogenization Model nu 22 not found!')
             nu_33 = 0._dp
             nuim_33 = 0._dp
             nu_33 = GetReal(CompParams, 'nu 33', Found)
             nuim_33 = GetReal(CompParams, 'nu 33 im', FoundIm)
             IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal ('LocalMatrix', 'Homogenization Model nu 33 not found!')
             Nu_el = CMPLX(0.0d0, 0.0d0, kind=dp)
             Nu_el(1,1,1:n) = nu_11(1:n) + im * nuim_11(1:n)
             Nu_el(2,2,1:n) = nu_22(1:n) + im * nuim_22(1:n)
             Nu_el(3,3,1:n) = nu_33(1:n) + im * nuim_33(1:n)

             sigma_33 = GetReal(CompParams, 'sigma 33', Found)
             IF ( .NOT. Found ) sigma_33 = 0._dp
             sigmaim_33 = GetReal(CompParams, 'sigma 33 im', FoundIm)
             IF ( .NOT. FoundIm ) sigmaim_33 = 0._dp
             IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal ('LocalMatrix', 'Homogenization Model sigma 33 not found!')

             Tcoef = 0._dp
             Tcoef(1,1,1:n) = sigma_33(1:n) + im * sigmaim_33(1:n) ! stranded uses only sigma 33
             Tcoef(2,2,1:n) = Tcoef(1,1,1:n)
             Tcoef(3,3,1:n) = Tcoef(1,1,1:n)
             IF (dim == 3) THEN
               CALL GetElementRotM(Element, RotM, n)
               DO k = 1,n
                 Nu_el(1:3,1:3,k) = MATMUL(MATMUL(RotM(1:3,1:3,k), Nu_el(1:3,1:3,k)), TRANSPOSE(RotM(1:3,1:3,k)))
               END DO
             END IF
           END BLOCK
         END IF
         !nofturns = GetConstReal(CompParams, 'Number of Turns', Found)
         !IF (.NOT. Found) CALL Fatal('MagnetoDynamicsCalcFields','Stranded Coil: Number of Turns not found!')
       CASE ('massive')
         VvarId = GetInteger (CompParams, 'Circuit Voltage Variable Id', Found)
         IF (.NOT. Found) CALL Fatal ('MagnetoDynamicsCalcFields', 'Circuit Voltage Variable Id not found!')

       CASE ('foil winding')
         CALL GetLocalSolution(alpha,'Alpha')
         
         IF (dim == 3) CALL GetElementRotM(Element, RotM, n)

         VvarId = GetInteger (CompParams, 'Circuit Voltage Variable Id', Found)
         IF (.NOT. Found) CALL Fatal ('MagnetoDynamicsCalcFields', 'Circuit Voltage Variable Id not found!')

         coilthickness = GetConstReal(CompParams, 'Coil Thickness', Found)
         IF (.NOT. Found) CALL Fatal('MagnetoDynamicsCalcFields','Foil Winding: Coil Thickness not found!')

         !nofturns = GetConstReal(CompParams, 'Number of Turns', Found)
         !IF (.NOT. Found) CALL Fatal('MagnetoDynamicsCalcFields','Foil Winding: Number of Turns not found!')

         VvarDofs = GetInteger (CompParams, 'Circuit Voltage Variable dofs', Found)
         IF (.NOT. Found) CALL Fatal ('MagnetoDynamicsCalcFields', 'Circuit Voltage Variable dofs not found!')
         ! in case of a foil winding, transform the conductivity tensor:
         ! -------------------------------------------------------------
        
         IF (dim == 3) THEN
             DO k = 1,n
               Tcoef(1:3,1:3,k) = MATMUL(MATMUL(RotM(1:3,1:3,k), Tcoef(1:3,1:3,k)), TRANSPOSE(RotM(1:3,1:3,k)))
             END DO
         END IF
       CASE DEFAULT
         CALL Fatal ('MagnetoDynamicsCalcFields', 'Non existent Coil Type Chosen!')
       END SELECT
     END IF


     !---------------------------------------------------------------------------------------------
     R_Z = CMPLX(0.0_dp, 0.0_dp, kind=dp)
     HasTensorReluctivity = .FALSE.

     IF ( .NOT. ( HasHBCurve .OR. HasReluctivityFunction ) ) THEN
       ! Seek reluctivity as complex-valued: A given reluctivity can be a tensor 
       CALL GetReluctivity(Material,Reluct_Z,n,HasTensorReluctivity)
       IF (HasTensorReluctivity) THEN
         IF (SIZE(Reluct_Z,1)==1 .AND. SIZE(Reluct_Z,2)==1) THEN
           l = MIN(SIZE(R_Z), SIZE(Reluct_Z,3))
           R_Z(1:l) = Reluct_Z(1,1,1:l)
           HasTensorReluctivity = .FALSE.
         ELSE
           R_Z = CMPLX(0.0_dp, 0.0_dp, kind=dp)
         END IF
       ELSE
         ! Seek via a given permeability: In this case the reluctivity will be 
         ! a complex scalar:
         CALL GetReluctivity(Material,R_Z,n)
       END IF
     END IF

     HasVelocity = .FALSE.
     HasLorenzVelocity = .FALSE.
     HasAngularVelocity = .FALSE.
     IF(ASSOCIATED(BodyForce)) THEN
       CALL GetRealVector( BodyForce, omega_velo, 'Angular velocity', HasAngularVelocity)
       CALL GetRealVector( BodyForce, lorentz_velo, 'Lorentz velocity', HasLorenzVelocity)
       HasVelocity = HasAngularVelocity .OR. HasLorenzVelocity
     END IF
     
     ! Calculate nodal fields:
     ! -----------------------
     pRef = ( dim==3 .AND. PiolaVersion ) .OR. isPelement(element)
     IF( ElementalMode >= 3 ) THEN
       IF( ElementalMode == 3 ) THEN
         IP = CornerGaussPoints(Element, EdgeBasis=dim==3, PReferenceElement=pRef)
       ELSE
         IP = CenterGaussPoints(Element, EdgeBasis=dim==3, PReferenceElement=pRef)
       END IF
     ELSE 
       IP = GaussPoints(Element, EdgeBasis=dim==3, PReferenceElement=pRef, EdgeBasisDegree=EdgeBasisDegree)
     END IF

     MASS  = 0._dp
     FORCE = 0._dp
     E = 0._dp; B=0._dp

     haszirka = .FALSE.
     if(ASSOCIATED(MFS) .OR. ASSOCIATED(el_MFS)) THEN
       CALL GetHystereticMFS(Element, force(:,4:6), pSolver, HasZirka, CSymmetry=CSymmetry)
     end if

     DO j = 1,IP % n
       IF(dim == 2 ) THEN
         stat = ElementInfo(Element,Nodes,IP % u(j),IP % v(j),IP % w(j),&
             detJ,Basis,dBasisdx,USolver=pSolver)
       ELSE
         stat = ElementInfo( Element, Nodes, IP % U(j), IP % V(j), IP % W(j), &
             detJ, Basis, dBasisdx, &
             EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver ) 
       END IF         
       s = IP % s(j) * detJ

       grads_coeff = -1._dp/GetCircuitModelDepth()
       IF( CSymmetry ) THEN
         xcoord = SUM( Basis(1:n) * Nodes % x(1:n) )
         grads_coeff = grads_coeff/xcoord
         s = s * xcoord 
       END IF
                
       
       DO k=1,vDOFs
         SELECT CASE(dim)
         CASE(2)
            ! This has been done with the same sign convention as in MagnetoDynamics2D:
            ! -------------------------------------------------------------------------
            IF ( CSymmetry ) THEN
              B(k,1) = -SUM( SOL(k,1:nd) * dBasisdx(1:nd,2) )
              B(k,2) =  SUM( SOL(k,1:nd) * dBasisdx(1:nd,1) ) &
                         + SUM( SOL(k,1:nd) * Basis(1:nd) ) / xcoord
              B(k,3) = 0._dp
            ELSE
              B(k,1) =  SUM( SOL(k,1:nd) * dBasisdx(1:nd,2) )
              B(k,2) = -SUM( SOL(k,1:nd) * dBasisdx(1:nd,1) )
              B(k,3) = 0._dp
            END IF
         CASE(3)
            B(k,:) = MATMUL( SOL(k,np+1:nd), RotWBasis(1:nd-np,:) )
         END SELECT
       END DO

       IF(ImposeCircuitCurrent .and. ItoJCoeffFound) THEN
         IF (CoilUseWvec) THEN
           wvec = ListGetElementVectorSolution( Wvec_h, Basis, Element, dofs = dim )
         ELSE
           wvec = -MATMUL(Wbase(1:n), dBasisdx(1:n,:))
           IF(SUM(wvec**2._dp) > AEPS) THEN
             wvec = wvec/SQRT(SUM(wvec**2._dp))
           ELSE
             wvec = [0.0_dp, 0.0_dp, 1.0_dp]
           END IF
         END IF
       END IF

       ! Compute the velocity field in the form v + w x r:
       ! -------------------------------------------------
       IF(HasVelocity) THEN
         rot_velo = 0.0_dp
         angular_velo = 0.0_dp
         IF( HasAngularVelocity ) THEN
           angular_velo(1) = SUM(basis(1:n)*omega_velo(1,1:n))
           angular_velo(2) = SUM(basis(1:n)*omega_velo(2,1:n))
           angular_velo(3) = SUM(basis(1:n)*omega_velo(3,1:n))
           DO k=1,n
             rot_velo(1:3) = rot_velo(1:3) + CrossProduct(omega_velo(1:3,k), [ &
                 Basis(k) * Nodes % x(k), &
                 Basis(k) * Nodes % y(k), &
                 Basis(k) * Nodes % z(k)])
           END DO
         END IF
         IF( HasLorenzVelocity ) THEN
           rot_velo(1:3) = rot_velo(1:3) + [ &
               SUM(basis(1:n)*lorentz_velo(1,1:n)), &
               SUM(basis(1:n)*lorentz_velo(2,1:n)), &
               SUM(basis(1:n)*lorentz_velo(3,1:n))]
         END IF
       END IF
       !-------------------------------
       ! The conductivity as a tensor
       ! -------------------------------
       !C_ip  = SUM( Basis(1:n)*C(1:n) )
       DO k=1,3
         DO l=1,3
           CMat_ip(k,l) = SUM( Tcoef(k,l,1:n) * Basis(1:n) )
         END DO
       END DO
       BodyForceCurrDens_ip(1:3) = 0._dp
       IF(ImposeBodyForceCurrent) THEN
         DO l=1,3
           BodyForceCurrDens_ip(l) = SUM(BodyForceCurrDens(l,1:n)*Basis(1:n))
         END DO
       END IF
       
       IF (vDOFs > 1) THEN   ! Complex case (harmonic case)
         IF (CoilType /= 'stranded') THEN
           ! -j * Omega A
           SELECT CASE(dim)
           CASE(2)
             E(1,:) = 0._dp
             E(2,:) = 0._dp
             E(1,3) =  Omega*SUM(SOL(2,1:nd) * Basis(1:nd))
             E(2,3) = -Omega*SUM(SOL(1,1:nd) * Basis(1:nd))
           CASE(3)
             E(1,:) =  Omega*MATMUL(SOL(2,np+1:nd),WBasis(1:nd-np,:))
             E(2,:) = -Omega*MATMUL(SOL(1,np+1:nd),WBasis(1:nd-np,:))
           END SELECT
         ELSE
           E(1,:) = 0._dp
           E(2,:) = 0._dp
         END IF

         localV=0._dp
         SELECT CASE (CoilType)

         CASE ('stranded')
           SELECT CASE(dim)
           CASE(2)
             wvec = [0._dp, 0._dp, 1._dp]
           CASE(3)

             IF (CoilUseWvec) THEN
               wvec = ListGetElementVectorSolution( Wvec_h, Basis, Element, dofs = dim )
             ELSE
               wvec = -MATMUL(Wbase(1:np), dBasisdx(1:np,:))
             ! btw. This does not allow ununiform windings... I don't fix it now. -ettaka
               wvec = wvec/SQRT(SUM(wvec**2._dp)) !Why were we doing this? 04132021 -ettaka
             END IF
           END SELECT

           IF(CMat_ip(3,3) /= 0._dp ) THEN
             imag_value = LagrangeVar % Values(IvarId) + im * LagrangeVar % Values(IvarId+1)
             E(1,:) = E(1,:)+REAL(imag_value * N_j * wvec / CMat_ip(3,3))
             E(2,:) = E(2,:)+AIMAG(imag_value * N_j * wvec / CMat_ip(3,3))
           END IF

         CASE ('massive')
           localV(1) = localV(1) + LagrangeVar % Values(VvarId) * CircEqVoltageFactor
           localV(2) = localV(2) + LagrangeVar % Values(VvarId+1) * CircEqVoltageFactor
           SELECT CASE(dim)
           CASE(2)
             E(1,3) = E(1,3)-localV(1) * grads_coeff
             E(2,3) = E(2,3)-localV(2) * grads_coeff
           CASE(3)
             E(1,:) = E(1,:)-localV(1) * MATMUL(Wbase(1:np), dBasisdx(1:np,:))
             E(2,:) = E(2,:)-localV(2) * MATMUL(Wbase(1:np), dBasisdx(1:np,:))
           END SELECT
           
         CASE ('foil winding')
           localAlpha = coilthickness *SUM(alpha(1:np) * Basis(1:np)) 
           DO k = 1, VvarDofs-1
             Reindex = 2*k
             Imindex = Reindex+1
             localV(1) = localV(1) + LagrangeVar % Values(VvarId+Reindex) * localAlpha**(k-1) * CircEqVoltageFactor
             localV(2) = localV(2) + LagrangeVar % Values(VvarId+Imindex) * localAlpha**(k-1) * CircEqVoltageFactor
           END DO
           SELECT CASE(dim)
           CASE(2)
             E(1,3) = E(1,3)-localV(1) * grads_coeff
             E(2,3) = E(2,3)-localV(2) * grads_coeff
           CASE(3)
             IF (CoilUseWvec) THEN
               wvec = ListGetElementVectorSolution( Wvec_h, Basis, Element, dofs = dim )
             ELSE
               wvec = MATMUL(Wbase(1:np), dBasisdx(1:np,:))
             END IF

             E(1,:) = E(1,:)-localV(1) * wvec
             E(2,:) = E(2,:)-localV(2) * wvec
           END SELECT

         CASE DEFAULT
           SELECT CASE(dim)
           CASE(2)
             IF (HasLorenzVelocity) THEN
               ! Add v x curl A:
               IF (CSymmetry) THEN
                 E(1,3) = E(1,3) - rot_velo(1) * B(1,2) + rot_velo(2) * B(1,1)
                 E(2,3) = E(2,3) - rot_velo(1) * B(2,2) + rot_velo(2) * B(2,1)
               ELSE
                 E(1,3) = E(1,3) + rot_velo(1) * B(1,2) - rot_velo(2) * B(1,1)
                 E(2,3) = E(2,3) + rot_velo(1) * B(2,2) - rot_velo(2) * B(2,1)
               END IF
             END IF
             !
             ! To make this perfect, the electric field corresponding to the source
             ! should be returned on the source region
             !
           CASE(3)
             ! -Grad(V)
             E(1,:) = E(1,:) - MATMUL(SOL(1,1:np), dBasisdx(1:np,:))
             E(2,:) = E(2,:) - MATMUL(SOL(2,1:np), dBasisdx(1:np,:))

             IF (HasVelocity) THEN
               !
               ! Add v x curl A so as to handle the steady amplitude solution of
               ! the time harmonic equations. Multiplication with the electric 
               ! conductivity will give the current density with respect to 
               ! the fixed frame:
               !
               E(1,:) = E(1,:) + CrossProduct(rot_velo, &
                   MATMUL(SOL(1,np+1:nd), RotWBasis(1:nd-np,:)))
               E(2,:) = E(2,:) + CrossProduct(rot_velo, &
                   MATMUL(SOL(2,np+1:nd), RotWBasis(1:nd-np,:)))
             END IF

             IF( ImposeBodyForcePotential ) THEN
               E(1,:) = E(1,:) - MATMUL(ElPotSol(1,1:n), dBasisdx(1:n,:))
             END IF
           END SELECT

         END SELECT
         
       ELSE   ! Real case (transient case)
         E(1,:) = 0._dp
         IF (CoilType /= 'stranded') THEN 
           SELECT CASE(dim)
           CASE(2)
             E(1,3) = -SUM(PSOL(1:nd) * Basis(1:nd))
           CASE(3)
             E(1,:) = -MATMUL(PSOL(np+1:nd), Wbasis(1:nd-np,:))
           END SELECT
         END IF
         localV=0._dp

         SELECT CASE (CoilType)
         CASE ('stranded')
           SELECT CASE(dim)
           CASE(2)
             wvec = [0._dp, 0._dp, 1._dp]
             IF(CMat_ip(1,1) /= 0._dp ) &
               E(1,:) = E(1,:)+ LagrangeVar % Values(IvarId) * N_j * wvec / CMat_ip(1,1)
           CASE(3)
             IF (CoilUseWvec) THEN
               wvec = ListGetElementVectorSolution( Wvec_h, Basis, Element, dofs = dim )
             ELSE
               wvec = -MATMUL(Wbase(1:np), dBasisdx(1:np,:))
               wvec = wvec/SQRT(SUM(wvec**2._dp))
             END IF
             IF(CMat_ip(3,3) /= 0._dp ) THEN
               E(1,:) = E(1,:)+ LagrangeVar % Values(IvarId) * N_j * wvec / CMat_ip(3,3)
             ELSE IF (.NOT. ImposeCircuitCurrent) THEN
               CircuitCurrent = LagrangeVar % Values(IvarId)
               ItoJCoeff = N_j
               ItoJCoeffFound = .TRUE.
             END IF
           END SELECT

         CASE ('massive')
           localV(1) = localV(1) + LagrangeVar % Values(VvarId) * CircEqVoltageFactor
           SELECT CASE(dim)
           CASE(2)
             E(1,3) = E(1,3)-localV(1) * grads_coeff
           CASE(3)
             E(1,:) = E(1,:)-localV(1) * MATMUL(Wbase(1:np), dBasisdx(1:np,:))
           END SELECT

         CASE ('foil winding')
           localAlpha = coilthickness *SUM(alpha(1:np) * Basis(1:np)) 
           DO k = 1, VvarDofs-1
             localV(1) = localV(1) + LagrangeVar % Values(VvarId+k) * localAlpha**(k-1) * CircEqVoltageFactor
           END DO
           SELECT CASE(dim)
           CASE(2)
             E(1,3) = E(1,3)-localV(1) * grads_coeff
           CASE(3)
             E(1,:) = E(1,:)-localV(1) * MATMUL(Wbase(1:np), dBasisdx(1:np,:))
           END SELECT

         CASE DEFAULT
           SELECT CASE(dim)
           CASE(2)
             IF (HasLorenzVelocity) THEN
               ! Add v x curl A
               IF (CSymmetry) THEN
                 E(1,3) = E(1,3) - rot_velo(1) * B(1,2) + rot_velo(2) * B(1,1)
               ELSE
                 E(1,3) = E(1,3) + rot_velo(1) * B(1,2) - rot_velo(2) * B(1,1)
               END IF
             END IF
             !
             ! To make this perfect, the electric field corresponding to the source
             ! should be returned on the source region
             !
           CASE(3)
             IF (Transient) THEN
               E(1,:) = E(1,:) - MATMUL(SOL(1,1:np), dBasisdx(1:np,:))
             END IF

             IF (np > 0 .AND. .NOT. Transient) THEN
               E(1,:) = -MATMUL(SOL(1,1:np), dBasisdx(1:np,:))
               
               IF (HasVelocity) THEN
                 !
                 ! Add v x curl A so as to handle the steady state solution of
                 ! the evolutionary equations. Multiplication with the electric 
                 ! conductivity will give the current density with respect to 
                 ! the fixed frame:
                 !
                 E(1,:) = E(1,:) + CrossProduct(rot_velo, &
                     MATMUL(SOL(1,np+1:nd), RotWBasis(1:nd-np,:)))
               END IF
             ELSE IF ( PrecomputedElectricPot ) THEN
               E(1,:) = -MATMUL(ElPotSol(1,1:n), dBasisdx(1:n,:))
             END IF

             IF ( ImposeBodyForcePotential ) THEN
               E(1,:) = E(1,:) - MATMUL(ElPotSol(1,1:n), dBasisdx(1:n,:))
             END IF
           END SELECT


         END SELECT
       END IF
       
       Nu = CMPLX(0.0d0, 0.0d0, kind=dp)
       w_dens = 0._dp
       IF ( HasHBCurve ) THEN
         IF (RealField) THEN
           Babs=SQRT(SUM(B(1,:)**2))
         ELSE
           Babs = SQRT(SUM(B(1,:)**2 + B(2,:)**2))
         END IF
         Babs = MAX(Babs, 1.d-8)
         
         R_ip = ListGetFun(Material,'h-b curve', x=babs) / Babs

         BLOCK
           TYPE(ValueListEntry_t), POINTER :: ptr
           ptr => ListFind( Material,'h-b curve')
           w_dens = IntegrateCurve(ptr % TValues,ptr % FValues(1,1,:),ptr % CubicCoeff,0._dp,Babs,Found=Found)
           IF(.NOT. Found ) HbIntegProblem = .TRUE.
         END BLOCK

         DO k=1,3
           Nu(k,k) = CMPLX(R_ip, 0.0d0, kind=dp)
         END DO
       ELSE IF( HasReluctivityFunction ) THEN
         rdummy = ListGetElementReal( mu_h, Basis, Element, &
             GaussPoint = j, Rdim=mudim, Rtensor=MuTensor, DummyVals = B(1,:) )             
         Nu(1:3,1:3) = muTensor(1:3,1:3)                           
         w_dens = 0.5*SUM(B(1,:)*MATMUL(REAL(Nu), B(1,:)))
       ELSE IF (HomogenizationLoss .AND. CoilType == 'stranded' .and. HomogenizationModel) THEN
         DO k=1,3
           DO l=1,3
             Nu(k,l) = SUM( Nu_el(k,l,1:n) * Basis(1:n) )
           END DO
         END DO
       ELSE
         IF (HasTensorReluctivity) THEN
           IF (SIZE(Reluct_Z,2) == 1) THEN
             DO k = 1, MIN(3, SIZE(Reluct_Z,1))
               Nu(k,k) = SUM(Basis(1:n)*Reluct_Z(k,1,1:n))
             END DO
           ELSE
             DO k = 1, MIN(3, SIZE(Reluct_Z,1))
               DO l = 1, MIN(3, SIZE(Reluct_Z,2))
                 Nu(k,l) = sum(Basis(1:n)*Reluct_Z(k,l,1:n))
               END DO
             END DO
           END IF
           R_ip = 0.0d0
         ELSE
           R_ip_Z = SUM(Basis(1:n)*R_Z(1:n))
           DO k=1,3
             Nu(k,k) = R_ip_Z
           END DO
           ! 
           ! The calculation of the Maxwell stress tensor doesn't yet support
           ! a tensor-form reluctivity. Create the scalar reluctivity parameter
           ! so that the Maxwell stress tensor may be calculated. The complex 
           ! part will be ignored.
           !
           R_ip = REAL(R_ip_Z)
         END IF
         IF (RealField) THEN
           w_dens = 0.5*SUM(B(1,:)*MATMUL(REAL(Nu), B(1,:)))
         ELSE
           ! This yields twice the time average:
           w_dens = 0.5*( SUM(MATMUL(REAL(Nu), B(1,:)) * B(1,:)) + &
               SUM(MATMUL(REAL(Nu), B(2,:)) * B(2,:)) ) 
         END IF
       END IF
       PR_ip = SUM( Basis(1:n)*PR(1:n) )

       IF ( ASSOCIATED(MFS).OR.ASSOCIATED(EL_MFS) ) THEN
         DO l=1,3
           MG_ip(l) = SUM( Magnetization(l,1:n)*Basis(1:n) )
         END DO
       END IF

       IF( ASSOCIATED(VP).OR.ASSOCIATED(EL_VP) ) THEN
         DO l=1,vDOFs
           SELECT CASE(dim)
           CASE(2)
             VP_ip(l,1) = 0._dp
             VP_ip(l,2) = 0._dp
             VP_ip(l,3) = SUM(SOL(l,1:nd) * Basis(1:nd))
           CASE(3)
             VP_ip(l,:) = MATMUL(SOL(l,np+1:nd),WBasis(1:nd-np,:))
           END SELECT
         END DO
       END IF

       IF (RealField) THEN
         HdotB = SUM(MATMUL(REAL(Nu), B(1,:)) * B(1,:))
       ELSE
         HdotB = SUM(MATMUL(REAL(Nu), B(1,:)) * B(1,:)) + &
             SUM(MATMUL(REAL(Nu), B(2,:)) * B(2,:))
       END IF
       
       IF (ASSOCIATED(NF).OR.ASSOCIATED(EL_NF)) THEN
         NF_ip = 0._dp
         DO k=1,n
           DO l=1,3
             val = SUM(dBasisdx(k,1:3)*B(1,1:3))
             NF_ip(k,l) = NF_ip(k,l) - SUM(REAL(Nu(l,1:3)) * B(1,1:3)) * val + &
                 (HdotB-w_dens)*dBasisdx(k,l)
           END DO
         END DO

         IF (.NOT. RealField) THEN
           DO k=1,n
             DO l=1,3
               val = SUM(dBasisdx(k,1:3)*B(2,1:3))
               NF_ip(k,l) = NF_ip(k,l) - SUM(REAL(Nu(l,1:3)) * B(2,1:3)) * val
             END DO
           END DO
         END IF
       END IF
       
       Energy(1) = Energy(1) + s*0.5*ABS(PR_ip)*SUM(E**2)
       Energy(2) = Energy(2) + s*w_dens
       Energy(3) = Energy(3) + (HdotB - w_dens) * s

       IF (ElementalFields .OR. .NOT. ConstantMassMatrixInUse) THEN
         DO p=1,eq_n
           DO q=1,eq_n
             MASS(p,q)=MASS(p,q)+s*Basis(p)*Basis(q)
           END DO
         END DO
       END IF

       DO p=1,eq_n
         k = 0
         
         IF( ASSOCIATED(MFD) .OR. ASSOCIATED(EL_MFD) ) THEN
           DO l=1,vDOFs
             FORCE(p,k+1:fdim+k) = FORCE(p,k+1:fdim+k)+s*B(l,1:fdim)*Basis(p)
             k = k+fdim
           END DO
         END IF
           
         IF ( (ASSOCIATED(MFS) .OR. ASSOCIATED(EL_MFS)) ) THEN
           ! Don't really know how to compute MFS with Zirka
           ! Skipping it will certainly cause errors since k becomes invalid
           !IF(.NOT. HasZirka) then
           IF (RealField) THEN
             FORCE(p,k+1:k+fdim) = FORCE(p,k+1:k+fdim) + &
                 s * (MATMUL(REAL(Nu(1:fdim,1:fdim)), B(1,1:fdim)) - REAL(MG_ip(1:fdim))) * Basis(p) 
             k = k+fdim
           ELSE
             FORCE(p,k+1:k+fdim) = FORCE(p,k+1:k+fdim) + s * ( &
                 MATMUL(REAL(Nu(1:fdim,1:fdim)), B(1,1:fdim)) - &
                 MATMUL(AIMAG(Nu(1:fdim,1:fdim)), B(2,1:fdim)) - REAL(MG_ip(1:fdim))) * Basis(p)
             k = k+fdim
             FORCE(p,k+1:k+fdim) = FORCE(p,k+1:k+fdim) + s * ( &
                 MATMUL(AIMAG(Nu(1:fdim,1:fdim)), B(1,1:fdim)) + &
                 MATMUL(REAL(Nu(1:fdim,1:fdim)), B(2,1:fdim)) - AIMAG(MG_ip(1:fdim))) * Basis(p) 
             k = k+fdim
           END IF
           !ELSE
           ! Never here?
           !FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)-s*(REAL(MG_ip))*Basis(p)
           !END IF
         END IF

         IF ( ASSOCIATED(VP).OR.ASSOCIATED(EL_VP)) THEN
           DO l=1,vDOFs
             FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+s*VP_ip(l,:)*Basis(p)
             k = k+3
           END DO
         END IF

         IF ( ASSOCIATED(EF).OR.ASSOCIATED(EL_EF)) THEN
           DO l=1,vDOFs
             FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3) + s*E(l,:)*Basis(p)
             k = k+3
           END DO
         END IF

         IF ( ASSOCIATED(CD).OR.ASSOCIATED(EL_CD)) THEN
           IF (ItoJCoeffFound) THEN
             IF (Vdofs == 1) THEN
               DO l=1,3
                 CC_J(1,l) = ItoJCoeff*wvec(l)*CircuitCurrent
               END DO
             ELSE
               CALL Fatal('MagnetoDynamicsCalcFields','Complex circuit current imposing is not implemented')
             END IF
           ELSE 
             CC_J(1,:) = 0.0_dp
           END IF

           IF (Vdofs == 1) THEN
              DO l=1,3
                JatIP(1,l) =  SUM( REAL(CMat_ip(l,1:3)) * E(1,1:3) ) + CC_J(1,l) + REAL(BodyForceCurrDens_ip(l)) 
                !
                ! No need for a "HasVelocity" check: the effect of v x B is already inbuilt into 
                ! the definition of the E-field
                ! IF( HasVelocity ) THEN
                !   JatIP(1,l) = JatIP(1,l) + SUM( REAL(CMat_ip(l,1:3)) * CrossProduct(rot_velo, B(1,1:3)))
                ! END IF
                !
                FORCE(p,k+l) = FORCE(p,k+l)+s*JatIp(1,l)*Basis(p)
              END DO
              k = k+3
           ELSE
              DO l=1,3
                JatIp(1,l) = SUM( REAL(CMat_ip(l,1:3)) * E(1,1:3) ) - &
                             SUM( AIMAG(CMat_ip(l,1:3)) * E(2,1:3) ) + REAL(BodyForceCurrDens_ip(l))
                !
                ! No need for a "HasVelocity" check: the effect of v x B is already inbuilt into 
                ! the definition of the E-field
                ! IF( HasVelocity ) THEN
                !   JatIp(1,l) = JatIp(1,l) + SUM( REAL(CMat_ip(l,1:3)) * CrossProduct(rot_velo, B(1,1:3)) ) - &
                !                SUM( AIMAG(CMat_ip(l,1:3)) * CrossProduct(rot_velo, B(2,1:3)) )
                ! END IF
                !
                FORCE(p,k+l) = FORCE(p,k+l)+s*JatIp(1,l)*Basis(p)
              END DO
              k = k+3
              DO l=1,3
                JatIp(2,l) = SUM( AIMAG(CMat_ip(l,1:3)) * E(1,1:3) ) + &
                             SUM( REAL(CMat_ip(l,1:3)) * E(2,1:3) ) + AIMAG(BodyForceCurrDens_ip(l))
                !
                ! No need for a "HasVelocity" check: the effect of v x B is already inbuilt into 
                ! the definition of the E-field                
                ! IF( HasVelocity ) THEN
                !   JatIp(2,l) = JatIp(2,l) + SUM( AIMAG(CMat_ip(l,1:3)) * CrossProduct(rot_velo, B(1,1:3)) ) + &
                !                SUM( REAL(CMat_ip(l,1:3)) * CrossProduct(rot_velo, B(2,1:3)) )
                ! END IF
                !
                FORCE(p,k+l) = FORCE(p,k+l)+s*JatIp(2,l)*Basis(p)
              END DO
              k = k+3
           END IF
         END IF

         IF ( ASSOCIATED(JXB).OR.ASSOCIATED(EL_JXB)) THEN
           IF (.NOT. ASSOCIATED(CD) .AND. .NOT. ASSOCIATED(EL_CD)) THEN
             CALL Warn('MagnetoDynamicsCalcFields', 'Cannot Calculate JxB since Current Density is not calculated!')
           ELSE
             IF (Vdofs == 1) THEN
               JXBatIP(1,:) = crossproduct(JatIP(1,:),B(1,:))
               DO l=1,dim
                 FORCE(p,k+l) = FORCE(p,k+l)+s*JXBatIP(1,l)*Basis(p)
               END DO
               k = k+3
             ELSE
               IF( CSymmetry ) THEN
                 ! TODO: Have to figure out why cylindrical coords have opposite sign
                 JXBatIP(1,1) =   JatIP(2,3)*B(2,2) + JatIP(1,3)*B(1,2)
                 JXBatIP(1,2) = - JatIP(2,3)*B(2,1) - JatIP(1,3)*B(1,1)
                 JXBatIP(1,3) =   0.0_dp

                 JXBatIP(2,1) =   JatIP(2,3)*B(1,2) - JatIP(1,3)*B(2,2)
                 JXBatIP(2,2) = - JatIP(2,3)*B(1,1) + JatIP(1,3)*B(2,1)
                 JXBatIP(2,3) =   0.0_dp
               ELSE
                 JXBatIP(1,1) =   JatIP(2,2)*B(2,3) - JatIP(2,3)*B(2,2) + JatIP(1,2)*B(1,3) - JatIP(1,3)*B(1,2)
                 JXBatIP(1,2) =   JatIP(2,3)*B(2,1) - JatIP(2,1)*B(2,3) + JatIP(1,3)*B(1,1) - JatIP(1,1)*B(1,3)
                 JXBatIP(1,3) =   JatIP(2,1)*B(2,2) - JatIP(2,2)*B(2,1) + JatIP(1,1)*B(1,2) - JatIP(1,2)*B(1,1)

                 JXBatIP(2,1) =   JatIP(2,2)*B(1,3) - JatIP(2,3)*B(1,2) - JatIP(1,2)*B(2,3) + JatIP(1,3)*B(2,2)
                 JXBatIP(2,2) =   JatIP(2,3)*B(1,1) - JatIP(2,1)*B(1,3) - JatIP(1,3)*B(2,1) + JatIP(1,1)*B(2,3)
                 JXBatIP(2,3) =   JatIP(2,1)*B(1,2) - JatIP(2,2)*B(1,1) - JatIP(1,1)*B(2,2) + JatIP(1,2)*B(2,1)
               END IF

               JXBatIP = 0.5_dp*JXBatIP

               DO l=1,dim
                 FORCE(p,k+l) = FORCE(p,k+l)+s*JXBatIP(1,l)*Basis(p)
               END DO
               k = k+3
               DO l=1,dim
                 FORCE(p,k+l) = FORCE(p,k+l)+s*JXBatIP(2,l)*Basis(p)
               END DO
               k = k+3
             END IF
           END IF
         END IF

         IF ( ASSOCIATED(FWP).OR.ASSOCIATED(EL_FWP)) THEN
           IF (Vdofs == 1) THEN
              FORCE(p,k+1) = FORCE(p,k+1)+s*LocalV(1)*Basis(p)
              k = k+1
           ELSE
              FORCE(p,k+1) = FORCE(p,k+1)+s*LocalV(1) * Basis(p)
              k = k+1
              FORCE(p,k+1) = FORCE(p,k+1)+s*LocalV(2) * Basis(p)
              k = k+1
           END IF
         END IF

         IF ( ASSOCIATED(MPerm).OR.ASSOCIATED(EL_MPerm)) THEN
           IF (Vdofs == 1) THEN
              FORCE(p,k+1) = FORCE(p,k+1)+s*1_dp/Nu(1,1)/mu0*Basis(p)
              k = k+1
           ELSE
              FORCE(p,k+1) = FORCE(p,k+1)+s*REAL(1_dp/Nu(1,1)/mu0)*Basis(p)
              k = k+1
              FORCE(p,k+1) = FORCE(p,k+1)+s*AIMAG(1_dp/Nu(1,1)/mu0) * Basis(p)
              k = k+1
           END IF
         END IF

         IF (vDOFS == 1) THEN
           IF ( JouleHeatingFromCurrent .AND. (ASSOCIATED(CD).OR.ASSOCIATED(EL_CD)) ) THEN
             ! The Joule heating power per unit volume: J.E = J.J/sigma 
             Coeff = 0.0_dp
             DO l=1,3
               IF( REAL( CMat_ip(l,l) )  > EPSILON( Coeff ) ) THEN
                 Coeff = Coeff + JatIP(1,l) * JatIP(1,l) / REAL( CMat_ip(l,l) ) * &
                     Basis(p) * s 
               END IF
             END DO
           ELSE 
             ! The Joule heating power per unit volume: J.E = (sigma * E).E
             Coeff = SUM( MATMUL( REAL(CMat_ip(1:3,1:3)), TRANSPOSE(E(1:1,1:3)) ) * &
                 TRANSPOSE(E(1:1,1:3)) ) * Basis(p) * s
           END IF
           !
           !
           ! No need for a "HasVelocity" check: the effect of v x B is already inbuilt into 
           ! the definition of the J-field via J's dependence on the E-field
           !
           ! IF (HasVelocity) THEN
           !   Coeff = Coeff + SUM(MATMUL(REAL(CMat_ip), CrossProduct(rot_velo, B(1,:))) * &
           !       CrossProduct(rot_velo,B(1,:)))*Basis(p)*s
           ! END IF
           !
         ELSE
           ! Now Power = J.conjugate(E), with the possible imaginary component neglected.         
           Coeff = HarmPowerCoeff * (SUM( MATMUL( REAL(CMat_ip(1:3,1:3)), TRANSPOSE(E(1:1,1:3)) ) * &
               TRANSPOSE(E(1:1,1:3)) ) * Basis(p) * s - &
               SUM( MATMUL( AIMAG(CMat_ip(1:3,1:3)), TRANSPOSE(E(2:2,1:3)) ) * &
               TRANSPOSE(E(1:1,1:3)) ) * Basis(p) * s + &
               SUM( MATMUL( AIMAG(CMat_ip(1:3,1:3)), TRANSPOSE(E(1:1,1:3)) ) * &
               TRANSPOSE(E(2:2,1:3)) ) * Basis(p) * s + &               
               SUM( MATMUL( REAL(CMat_ip(1:3,1:3)), TRANSPOSE(E(2:2,1:3)) ) * &
               TRANSPOSE(E(2:2,1:3)) ) * Basis(p) * s)
           !
           ! Again, no need for a "HasVelocity" check
           ! IF (HasVelocity) THEN
           !   Coeff = Coeff + HarmPowerCoeff * (SUM( MATMUL( REAL(CMat_ip(1:3,1:3)), CrossProduct(rot_velo, B(1,:)) ) * &
           !     CrossProduct(rot_velo, B(1,:)) ) * Basis(p) * s - &
           !     SUM( MATMUL( AIMAG(CMat_ip(1:3,1:3)), CrossProduct(rot_velo, B(2,:)) ) * &
           !     CrossProduct(rot_velo, B(1,:)) ) * Basis(p) * s + &
           !     SUM( MATMUL( AIMAG(CMat_ip(1:3,1:3)), CrossProduct(rot_velo, B(1,:)) ) * &
           !     CrossProduct(rot_velo, B(2,:)) ) * Basis(p) * s + &
           !     SUM( MATMUL( REAL(CMat_ip(1:3,1:3)), CrossProduct(rot_velo, B(2,:)) ) * &
           !     CrossProduct(rot_velo, B(2,:)) ) * Basis(p) * s)
           ! END IF
           !
         END IF

         Power = Power + Coeff
         IF ( ASSOCIATED(JH) .OR. ASSOCIATED(EL_JH) .OR. ASSOCIATED(NJH) ) THEN           
           FORCE(p,k+1) = FORCE(p,k+1) + Coeff
           k = k+1
           jh_k = k
         END IF
         
         !-------------------------------------------------
         ! Compute a loss estimate for cos and sin modes:
         !-------------------------------------------------
         IF (LossEstimation) THEN
           BodyLoss(3,BodyId) = BodyLoss(3,BodyId) + Coeff

           IF( OldLossKeywords ) THEN
             LossCoeff(1) = ListGetFun( Material,'Harmonic Loss Linear Coefficient',Freq,Found ) 
             LossCoeff(2) = ListGetFun( Material,'Harmonic Loss Quadratic Coefficient',Freq,Found ) 
           ELSE
             DO l=1,LossN               
               LossCoeff(l) = ListGetFun( Material,&
                   'Harmonic Loss Coefficient '//I2S(l),Freq, Found)      
             END DO
           END IF

             
           IF(MaterialExponents) THEN
             CALL GetLossExponents(Material,FreqPower,FieldPower,LossN,OldLossKeywords)
           END IF
           
           ! No losses to add if loss coefficient is not given
           IF( Found .OR. MaterialExponents ) THEN
             ElemLoss = 0.0_dp
             DO l=1,2
               ValAtIP = SUM( B(l,1:3) ** 2 )
               ElemLoss(1) = ElemLoss(1) + s * Basis(p) * LossCoeff(1) * ( Freq ** FreqPower(1) ) * ( ValAtIp ** FieldPower(1) )
               ElemLoss(2) = ElemLoss(2) + s * Basis(p) * LossCoeff(2) * ( Freq ** FreqPower(2) ) * ( ValAtIp ** FieldPower(2) )
               ComponentLoss(:,l) = ComponentLoss(:,l) + ElemLoss
               BodyLoss(:,BodyId) = BodyLoss(:,BodyId) + ElemLoss
             END DO
           ELSE
             ElemLoss = 0.0_dp
           END IF

           IF ( ASSOCIATED(ML) .OR. ASSOCIATED(EL_ML) ) THEN
             FORCE(p,k+1) = FORCE(p,k+1) + ElemLoss(1)
             k = k + 1
           END IF
           IF ( ASSOCIATED(ML2) .OR. ASSOCIATED(EL_ML2) ) THEN
             FORCE(p,k+1) = FORCE(p,k+1) + ElemLoss(2)
             k = k + 1
           END IF
         END IF

         IF ( HomogenizationLoss .AND. CoilType == 'stranded' .and. HomogenizationModel) THEN
           ! homogenization loss should be real part of im omega b . conj(h)/2
           BLOCK
             COMPLEX(KIND=dp) :: Bloc(3)=0._dp
             COMPLEX(KIND=dp) :: Hloc(3)=0._dp

             Bloc = B(1, 1:3) + im * B(2, 1:3)
             Hloc = MATMUL(Nu(1:3, 1:3), Bloc)
             Coeff = s * Basis(p) * REAL(im * Omega * SUM(Bloc * CONJG(Hloc))/2._dp)
           END BLOCK

           IF ( ASSOCIATED(PL) .OR. ASSOCIATED(EL_PL) ) THEN
             FORCE(p,k+1) = FORCE(p,k+1) + Coeff
             k = k + 1
           END IF

         END IF

         IF ( ASSOCIATED(MST).OR.ASSOCIATED(EL_MST)) THEN
           IF ( Vdofs==1 ) THEN
             DO l=1,3
               DO m=l,3
                 ST(l,m)=REAL(PR_ip)*E(1,l)*E(1,m)+R_ip*B(1,l)*B(1,m)
               END DO
               ST(l,l)=ST(l,l)-(REAL(PR_ip)*SUM(E(1,:)**2)+R_ip*SUM(B(1,:)**2))/2
             END DO
             DO l=1,6
               FORCE(p,k+l)=FORCE(p,k+l) + s*ST(ind1(l),ind2(l))*Basis(p)
             END DO
             k = k + 6
           ELSE
             DO l=1,3
               DO m=l,3
                 CST(l,m) = PR_ip*CMPLX(E(1,l),E(2,l),KIND=dp) * &
                                  CMPLX(E(1,m),E(2,m),KIND=dp)
                 CST(l,m) = CST(l,m) + &
                             R_ip*CMPLX(B(1,l),B(2,l),KIND=dp) * &
                                  CMPLX(B(1,m),B(2,m),KIND=dp)
               END DO
               CST(l,l) = CST(l,l) - &
                      (PR_ip*SUM(ABS(CMPLX(E(1,:),E(2,:)))**2)+ &
                        R_ip*SUM(ABS(CMPLX(B(1,:),B(2,:)))**2))/2
             END DO
             DO l=1,6
               FORCE(p,k+l)=FORCE(p,k+l) + s*REAL(CST(ind1(l),ind2(l)))*Basis(p)
             END DO
             k = k + 6
             DO l=1,6
               FORCE(p,k+l)=FORCE(p,k+l) + s*AIMAG(CST(ind1(l),ind2(l)))*Basis(p)
             END DO
             k = k + 6
           END IF
         END IF
         IF (ASSOCIATED(NF).OR.ASSOCIATED(EL_NF)) THEN
           IF(RealField) THEN
             FORCE(p,k+1:k+fdim) = FORCE(p,k+1:k+fdim) + s*NF_ip(p,1:fdim)
           ELSE
             FORCE(p,k+1:k+fdim) = FORCE(p,k+1:k+fdim) + 0.5*s*NF_ip(p,1:fdim)
           END IF
           k = k + fdim
         END IF
       END DO ! p
     END DO ! j

     IF(NodalFields) THEN
       IF(.NOT. ConstantMassMatrixInUse ) THEN
         CALL DefaultUpdateEquations( MASS,Force(:,1))
       END IF
       Fsave => Solver % Matrix % RHS
       DO l=1,k
         Solver % Matrix % RHS => GForce(:,l)
         CALL DefaultUpdateForce(Force(:,l))
       END DO       
       Solver % Matrix % RHS => Fsave
     END IF

     IF(ElementalFields) THEN
       dofs = 0
       
       IF( ElementalMode == 1 ) THEN
         ! Perform classical mass lumping
         DO k=1,eq_n
           s = SUM(MASS(k,1:eq_n))
           MASS(k,1:eq_n) = 0.0_dp
           MASS(k,k) = s
         END DO
       END IF

       IF( ElementalMode /= 2 .AND. ElementalMode /= 4) THEN
         CALL LUdecomp(MASS,eq_n,pivot,Erroneous)
         IF (Erroneous) CALL Fatal('MagnetoDynamicsCalcFields', 'LU-decomposition fails')
       END IF

       CALL LocalSol(EL_MFD,  fdim*vdofs, n, eq_n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_MFS,  fdim*vdofs, n, eq_n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_VP,   3*vdofs, n, eq_n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_EF,   3*vdofs, n, eq_n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_CD,   3*vdofs, n, eq_n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_JXB,  3*vdofs, n, eq_n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_FWP,  1*vdofs, n, eq_n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_MPerm,  1*vdofs, n, eq_n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_JH,   1, n, eq_n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_ML,   1, n, eq_n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_ML2,  1, n, eq_n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_PL,  1, n, eq_n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_MST,  6*vdofs, n, eq_n, MASS, FORCE, pivot, Dofs)

       ! This is a nodal quantity
       CALL LocalCopy(EL_NF, fdim, eq_n, FORCE, Dofs)
     END IF
   END DO   

   !
   ! Some postprocessing of surface currents generated by skin BCs (in time-harmonic cases): 
   !----------------------------------------------------------------------------------------
   IF ( LayerBC ) THEN
     SurfPower = 0.0_dp
     DO i=1,GetNOFBoundaryElements()
       Element => GetBoundaryElement(i)
       BC => GetBC()
       IF (.NOT. ASSOCIATED(BC)) CYCLE

       IF( dim == 3 ) THEN
         SELECT CASE(GetElementFamily())
         CASE(1)
           CYCLE
         CASE(2)
           k = GetBoundaryEdgeIndex(Element,1)
           Element => Mesh % Edges(k)
         CASE(3,4)
           k = GetBoundaryFaceIndex(Element)
           Element => Mesh % Faces(k)
         END SELECT
       END IF
       IF (.NOT. ActiveBoundaryElement(Element)) CYCLE
       
       C_ip = ListGetCReal(BC, 'Layer Electric Conductivity', Found)
       IF(.NOT. Found) CYCLE

       mu_r = ListGetCReal(BC, 'Layer Relative Permeability', Found)
       IF (.NOT. Found) mu_r = 1.0_dp

       LMsol = 0.0_dp
       IF( ASSOCIATED( LagrangeVar ) ) THEN
         IF( ListGetLogical( BC,'Flux Integral BC', Found ) ) THEN
           k = GetBCId( Element ) 
           IF( ASSOCIATED( pSolver % MortarBCs ) ) THEN
             k = pSolver % MortarBCs(k) % rowoffset
             LMSol(1:2) = LagrangeVar % Values(k+1:k+2)
           END IF
         END IF
       END IF
       
       n = GetElementNOFNodes(Element)     
       nd = GetElementNOFDOFs(uElement=Element, uSolver=pSolver)
       np = n ! Note: the scalar potential should be present by default in the time-harmonic case
       
       CALL GetVectorLocalSolution(SOL, Pname, uElement=Element, uSolver=pSolver)
       CALL GetElementNodes(Nodes, Element)
       
       IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
           EdgeBasisDegree=EdgeBasisDegree)
       FORCE = 0.0_dp
       
       ComponentId=GetInteger( BC, 'Component', CircuitDrivenBC)
       IF (CircuitDrivenBC) THEN
         CompParams => GetComponentParams( Element )
         VvarId = GetInteger (CompParams, 'Circuit Voltage Variable Id', Found)
         IF (.NOT. Found) CALL Fatal ('MagnetoDynamicsCalcFields', 'Circuit Voltage Variable Id not found!')
       END IF

       DO j=1,IP % n
         IF( dim == 2 ) THEN
           stat = ElementInfo(Element, Nodes, IP % U(j), IP % V(j), IP % W(j), &
               detJ, Basis, dBasisdx)           
         ELSE 
           stat = ElementInfo( Element, Nodes, IP % U(j), IP % V(j), &
               IP % W(j), detJ, Basis, dBasisdx, &
               EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver ) 
         END IF
         
         val = SQRT(2.0_dp/(C_ip * Omega * 4.0d0 * PI * 1d-7 * mu_r)) ! The layer thickness
         Zs = CMPLX(1.0_dp, 1.0_dp, KIND=dp) / (C_ip*val)
         
         IF (.NOT. CircuitDrivenBC) THEN
           IF( dim == 2 ) THEN
             E(1,1:2) = 0.0_dp
             E(2,1:2) = 0.0_dp
             E(1,3) = Omega * SUM(SOL(2,1:nd) * Basis(1:nd))
             E(2,3) = -Omega * SUM(SOL(1,1:nd) * Basis(1:nd)) 
             
             ! Include the effect of constraint for surface current density.
             ! Currently only 2D model implemented where the effect goes into z-component only.
             E(1,3) = E(1,3) + Omega * LMSol(1)
             E(2,3) = E(2,3) + Omega * LMSol(2)
           ELSE
             E(1,:) = Omega * MATMUL(SOL(2,np+1:nd), WBasis(1:nd-np,:)) - MATMUL(SOL(1,1:np), dBasisdx(1:np,:))
             E(2,:) = -Omega * MATMUL(SOL(1,np+1:nd), WBasis(1:nd-np,:)) - MATMUL(SOL(2,1:np), dBasisdx(1:np,:))
           END IF
         ELSE
           ! we assume 3D massive coil here
           E(1,:) = Omega * MATMUL(SOL(2,np+1:nd), WBasis(1:nd-np,:))
           E(2,:) = -Omega * MATMUL(SOL(1,np+1:nd), WBasis(1:nd-np,:))

           localV(1) = localV(1) + LagrangeVar % Values(VvarId) * CircEqVoltageFactor
           localV(2) = localV(2) + LagrangeVar % Values(VvarId+1) * CircEqVoltageFactor
           E(1,:) = E(1,:)-localV(1) * MATMUL(Wbase(1:np), dBasisdx(1:np,:))
           E(2,:) = E(2,:)-localV(2) * MATMUL(Wbase(1:np), dBasisdx(1:np,:))
         END IF
         
         s = IP % s(j) * detJ
         IF( CSymmetry ) THEN
           xcoord = SUM( Basis(1:n) * Nodes % x(1:n) ) 
           s = s * xcoord 
         END IF
         
         ! Compute the (real) power to maintain the surface current j_S in terms of 
         ! the surface impedance from the power density P_S = 1/2 Real(1/Zs) E.conjugate(E)
         Coeff = HarmPowerCoeff * REAL(1.0_dp/Zs) * &
             (SUM(E(1,:)**2) + SUM(E(2,:)**2)) * s 
         
         SurfPower = SurfPower + Coeff

         ! Currently we only consider the effect of boundary to Joule heating
         IF ( ASSOCIATED(JH) .OR. ASSOCIATED(NJH) ) THEN           
           DO p=1,n
             FORCE(p,jh_k) = FORCE(p,jh_k) + Basis(p) * Coeff
           END DO
         END IF

         ! and surface current density (here SCD should give the tangential trace H x n)
         IF( ASSOCIATED( SCD ) ) THEN
           DO p=1,n
             k = SCD % Perm(Element % NodeIndexes(p))
             SCD % Values(6*k-5:6*k-3) = SCD % Values(6*k-5:6*k-3) + &
                 s * val * c_ip * Basis(p) * ( E(1,1:3) + E(2,1:3) ) / 2
             SCD % Values(6*k-2:6*k-0) = SCD % Values(6*k-2:6*k-0) + &
                 s * val * c_ip * Basis(p) * ( -E(1,1:3) + E(2,1:3) ) / 2
             SurfWeight(k) = SurfWeight(k) + s*Basis(p)
           END DO
         END IF

       END DO

       IF (NodalFields .AND. jh_k>0) THEN
         l = jh_k
         IF(dim==2) THEN
           CALL UpdateGlobalForce( GForce(:,l), &
               Force(1:eq_n,l), eq_n, 1, Solver % Variable % Perm(Element % NodeIndexes), UElement=Element)
         ELSE
           CALL UpdateGlobalForce( GForce(:,l), &
               Force(1:eq_n,l), eq_n, 1, Solver % Variable % Perm(Indexes(1:eq_n)), UElement=Element)
         END IF
       END IF
     END DO
   END IF
   
   DoAve = GetLogical(SolverParams,'Average Within Materials',Found)
      
   ! Assembly of the face terms:
   !----------------------------
   IF(.NOT. ConstantMassMatrixInUse ) THEN
     IF (GetLogical(SolverParams,'Discontinuous Galerkin',Found)) THEN
       IF (DoAve) THEN
         FORCE = 0.0_dp
         CALL AddLocalFaceTerms( MASS, FORCE(:,1) )
       END IF
     END IF
   END IF
   
   IF(NodalFields .OR. DoAve ) THEN
     FSave => NULL()
     IF( ASSOCIATED( Solver % Matrix ) ) THEN
       Fsave => Solver % Matrix % RHS
     END IF

     DOFs = 0
     CALL GlobalSol(MFD,  fdim*vdofs, Gforce, Dofs, EL_MFD)
     CALL GlobalSol(MFS,  fdim*vdofs, Gforce, Dofs, EL_MFS)
     CALL GlobalSol(VP ,  3*vdofs, Gforce, Dofs, EL_VP)
     CALL GlobalSol(EF,   3*vdofs, Gforce, Dofs, EL_EF)
     CALL GlobalSol(CD,   3*vdofs, Gforce, Dofs, EL_CD)
     CALL GlobalSol(JXB,  3*vdofs, Gforce, Dofs, EL_JXB)
     CALL GlobalSol(FWP,  1*vdofs, Gforce, Dofs, EL_FWP)
     CALL GlobalSol(MPerm,  1*vdofs, Gforce, Dofs, EL_MPerm)
     
     ! Nodal heating directly uses the loads 
     IF (ASSOCIATED(NJH)) THEN
       NJH % Values = Gforce(:,dofs+1)
       ! Update the dofs only if it not used as the r.h.s. for the following field
       IF(.NOT. (ASSOCIATED(JH) .OR. ASSOCIATED(EL_JH))) dofs = dofs + 1
     END IF
     CALL GlobalSol(JH ,  1      , Gforce, Dofs, EL_JH)

     CALL GlobalSol(ML ,  1      , Gforce, Dofs, EL_ML)
     CALL GlobalSol(ML2,  1      , Gforce, Dofs, EL_ML2)
     CALL GlobalSol(PL ,  1      , Gforce, Dofs, EL_PL)
     CALL GlobalSol(MST,  6*vdofs, Gforce, Dofs, EL_MST)

     IF (ASSOCIATED(NF)) THEN
       DO i=1,fdim
         dofs = dofs + 1
         NF % Values(i::fdim) = Gforce(:,dofs)
       END DO
     END IF
     IF(ASSOCIATED(FSave)) Solver % Matrix % RHS => Fsave
   END IF
   
   ! surface current density uses the loads
   IF (ASSOCIATED(SCD)) THEN
     DO i=1,3*vdofs
       WHERE( SurfWeight > 0 ) 
         SCD % Values(i::3*vdofs) = SCD % Values(i::3*vdofs) / SurfWeight(:)
       END WHERE
     END DO
   END IF
      

   IF( ListCheckPresentAnyComponent( Model, 'Flux linkage' ) ) THEN
     DO j=1,Model % NumberOfComponents
       CompParams => Model % Components(j) % Values       
       IF ( ListGetLogical( CompParams,'Flux linkage', Found ) ) THEN         
         s = ComponentStokesTheorem(Model, Mesh, CompParams, pSolver % Variable,.FALSE. )
         PRINT *,'Flux linkage:',j,s
         IF( ASSOCIATED(VP) ) THEN
           s = ComponentStokesTheorem(Model, Mesh, CompParams, VP,.FALSE. )
           PRINT *,'Flux linkage nodal:',j,s
         END IF
         s = ComponentStokesTheorem(Model, Mesh, CompParams, pSolver % Variable,.TRUE. )
         PRINT *,'Flux linkage averaged:',j,s
         IF( ASSOCIATED(VP) ) THEN
           s = ComponentStokesTheorem(Model, Mesh, CompParams, VP,.TRUE. )
           PRINT *,'Flux linkage nodal avereaged:',j,s
         END IF
       END IF       
     END DO
   END IF

   IF( ListCheckPresentAnyComponent( Model, 'Coil Energy' ) ) THEN          
     BLOCK 
       TYPE(Variable_t), POINTER :: CoilCurr
       REAL(KIND=dp), ALLOCATABLE :: CoilEnergy(:)
       INTEGER, POINTER :: MasterEntities(:)
            
       CoilCurr => VariableGet( Mesh % Variables,'CoilCurrent e',ThisOnly=.TRUE.)
       IF(.NOT. ASSOCIATED(CoilCurr)) THEN
         CoilCurr => VariableGet( Mesh % Variables,'CoilCurrent',ThisOnly=.TRUE.)
       END IF

       ALLOCATE(CoilEnergy(Model % NumberOfComponents))
       CoilEnergy = 0.0_dp
       
       DO j=1,Model % NumberOfComponents
         CompParams => Model % Components(j) % Values       
         IF ( ListGetLogical( CompParams,'Coil Energy', Found ) ) THEN         
           MasterEntities => ListGetIntegerArray( CompParams,'Master Bodies',Found )
           CoilEnergy(j) =  0.5_dp * ComponentCoilEnergy(Model, Mesh, MasterEntities, pSolver % Variable, CoilCurr )            
           PRINT *,'Coil Energy:',j,CoilEnergy(j)
           IF( ASSOCIATED(VP) ) THEN
             CoilEnergy(j) = 0.5_dp * ComponentCoilEnergy(Model, Mesh, MasterEntities, VP, CoilCurr) 
             PRINT *,'Coil Energy nodal A:',j,CoilEnergy(j)
           END IF
         END IF
       END DO
       PRINT *,'Total Coil Energy:',SUM(CoilEnergy)       
     END BLOCK
   END IF
   
   ! Lump componentwise forces and torques. 
   ! Prefer DG nodal force variable if air gap is present

   ! Warn if user has air gaps and no "nodal force e"
   HaveAirGap = ListCheckPresentAnyBC( Model, 'Air Gap Length' ) 
   UseElementalNF = ASSOCIATED( EL_NF ) .AND. ( .NOT. ASSOCIATED( NF ) .OR. HaveAirGap )
   
    
   IF( UseElementalNF ) THEN
      CALL Info('MagnetoDynamicsCalcFields','Doing elemental nodal force stuff!',Level=20)
            
     ! Collect nodal forces from airgaps
     CALL CalcBoundaryModels()

     ! Create a minimal discontinuous set such that discontinuity is only created
     ! when body has an air gap boundary condition. Only do the reduction for the 1st time.
     IF( .NOT. ASSOCIATED( SetPerm ) ) THEN
       CALL Info('MagnetoDynamicsCalcFields','Creating minimal elemental set',Level=10)
       SetPerm => MinimalElementalSet( Mesh,'db', Solver % Variable % Perm, &
         BcFlag = 'Air Gap Length', &
         NonGreedy = ListGetLogical( Solver % Values,'Nongreedy Jump',Found) ) 
     END IF

     ! Sum up (no averaging) the elemental fields such that each elemental nodes has also 
     ! the contributions of the related nodes in other elements
     CALL ReduceElementalVar( Mesh, EL_NF, SetPerm, TakeAverage = .FALSE.)

     DO j=1,Model % NumberOfComponents
       CompParams => Model % Components(j) % Values
       
       IF ( ListGetLogical( CompParams,'Calculate Magnetic Force', Found ) ) THEN

         CALL ComponentNodalForceReduction(Model, Mesh, CompParams, EL_NF, &
           Force = LumpedForce, SetPerm = SetPerm )

         WRITE( Message,'(A,3ES15.6)') 'Magnetic force reduced: > '&
           //TRIM(ListGetString(CompParams,'Name'))//' < :', LumpedForce
         CALL Info('MagnetoDynamicsCalcFields',Message,Level=6)           

         DO i=1,fdim
           CALL ListAddConstReal( CompParams,'res: magnetic force '//I2S(i), LumpedForce(i) )
         END DO

       END IF

       IF( ListGetLogical( CompParams,'Calculate Magnetic Torque', Found ) ) THEN         
         CALL ComponentNodalForceReduction(Model, Mesh, CompParams, EL_NF, &
           Torque = val, SetPerm = SetPerm )

         WRITE( Message,'(A,ES15.6)') 'Magnetic torque reduced: > '&
           //TRIM(ListGetString(CompParams,'Name'))//' < :', val
         CALL Info('MagnetoDynamicsCalcFields',Message,Level=6)           

         CALL ListAddConstReal( CompParams,'res: magnetic torque', val )
       END IF
     END DO
  ELSE
      CALL Info('MagnetoDynamicsCalcFields','Doing nodal nodal force stuff!',Level=20)
     
     DO j=1,Model % NumberOfComponents
       CompParams => Model % Components(j) % Values

       IF( ListGetLogical( CompParams,'Calculate Magnetic Force', Found ) ) THEN 
         CALL ComponentNodalForceReduction(Model, Mesh, CompParams, NF, &
           Force = LumpedForce )

         WRITE( Message,'(A,3ES15.6)') 'Magnetic force reduced: > '&
           //TRIM(ListGetString(CompParams,'Name'))//' < :', LumpedForce
         CALL Info('MagnetoDynamicsCalcFields',Message,Level=6)           

         DO i=1,fdim
           CALL ListAddConstReal( CompParams,'res: magnetic force '//I2S(i), LumpedForce(i) )
         END DO

       END IF

       IF( ListGetLogical( CompParams,'Calculate Magnetic Torque', Found ) ) THEN 

         ! Warn if user has air gaps and no "nodal force e" is available
         IF ( HaveAirGap ) THEN
           CALL Warn('MagnetoDynamicsCalcFields', 'Cannot calculate air gap &
             &forces correctly because elemental field "Nodal Force e" is not &
             &present.')
         END IF

         CALL ComponentNodalForceReduction(Model, Mesh, CompParams, NF, Torque = val )

         WRITE( Message,'(A,ES15.6)') 'Magnetic torque reduced: > '&
           //TRIM(ListGetString(CompParams,'Name'))//' < :', val
         CALL Info('MagnetoDynamicsCalcFields',Message,Level=6)           

         CALL ListAddConstReal( CompParams,'res: magnetic torque', val )
       END IF
     END DO
   END IF


   ! Perform parallel reductions 
   IF(Parallel) THEN
     Power = ParallelReduction(Power) / NoSlices
     IF( LayerBC ) SurfPower = ParallelReduction( SurfPower ) / NoSlices

     Energy(1) = ParallelReduction(Energy(1)) / NoSlices
     Energy(2) = ParallelReduction(Energy(2)) / NoSlices
     Energy(3) = ParallelReduction(Energy(3)) / NoSlices
     
     IF (LossEstimation) THEN
       DO j=1,2
         DO i=1,2
           ComponentLoss(j,i) = ParallelReduction(ComponentLoss(j,i)) / NoSlices
         END DO
       END DO

       TotalLoss = 0._dp
       DO j=1,3
         DO i=1,Model % NumberOfBodies
           BodyLoss(j,i) = ParallelReduction(BodyLoss(j,i)) / NoSlices
         END DO
         TotalLoss(j) = SUM( BodyLoss(j,:) )
       END DO
     END IF
   END IF

   IF( HbIntegProblem ) THEN
     CALL Warn('MagnetoDynamicsCalcFields','Could not integrate over H-B curve for magnetic energy!')
   END IF
   
   WRITE(Message,'(A,ES15.6)') 'Eddy current power: ', Power
   CALL Info( 'MagnetoDynamicsCalcFields', Message )
   CALL ListAddConstReal( Model % Simulation, 'res: Eddy current power', Power )

   IF( LayerBC ) THEN
     WRITE(Message,*) 'Surface current power (the Joule effect): ', SurfPower
     CALL Info( 'MagnetoDynamicsCalcFields', Message )
     CALL ListAddConstReal(Model % Simulation, 'res: Surface current power', SurfPower)
   END IF
   
   IF ( ListGetLogical( SolverParams,'Separate Magnetic Energy',Found ) ) THEN
     WRITE(Message,'(A,ES15.6)') 'Electric Field Energy: ', Energy(1)
     CALL Info( 'MagnetoDynamicsCalcFields', Message )
     WRITE(Message,'(A,ES15.6)') 'Magnetic Field Energy: ', Energy(2)
     CALL Info( 'MagnetoDynamicsCalcFields', Message )
     WRITE(Message,'(A,ES15.6)') 'Magnetic Coenergy: ', Energy(3)
     CALL Info( 'MagnetoDynamicsCalcFields', Message )
     CALL ListAddConstReal(Model % Simulation,'res: Electric Field Energy',Energy(1))
     CALL ListAddConstReal(Model % Simulation,'res: Magnetic Field Energy',Energy(2))
     CALL ListAddConstReal(Model % Simulation,'res: Magnetic Coenergy',Energy(3))
   ELSE
     WRITE(Message,'(A,ES15.6)') 'ElectroMagnetic Field Energy: ',SUM(Energy(1:2))
     CALL Info( 'MagnetoDynamicsCalcFields', Message )
     CALL ListAddConstReal(Model % Simulation,'res: ElectroMagnetic Field Energy',SUM(Energy(1:2)))
   END IF
   
   
   IF (LossEstimation) THEN
     CALL ListAddConstReal( Model % Simulation,'res: harmonic loss linear',TotalLoss(1) )
     CALL ListAddConstReal( Model % Simulation,'res: harmonic loss quadratic',TotalLoss(2) )
     CALL ListAddConstReal( Model % Simulation,'res: joule loss',TotalLoss(3) )

     DO k=1,2
       IF( k == 1 ) THEN
         CALL Info('MagnetoDynamicsCalcFields','Harmonic Loss Linear by components',Level=6)
       ELSE
         CALL Info('MagnetoDynamicsCalcFields','Harmonic Loss Quadratic by components',Level=6)
       END IF
       WRITE( Message,'(A,ES15.6)') 'Loss for cos mode: ', ComponentLoss(k,1)
       CALL Info('MagnetoDynamicsCalcFields', Message, Level=6 )
       WRITE( Message,'(A,ES15.6)') 'Loss for sin mode: ', ComponentLoss(k,2)
       CALL Info('MagnetoDynamicsCalcFields', Message, Level=6 )
       WRITE( Message,'(A,ES15.6)') 'Total loss: ',TotalLoss(k)
       CALL Info('MagnetoDynamicsCalcFields',Message, Level=5 )
     END DO

     DO k=1,3
       IF( TotalLoss(k) < TINY( TotalLoss(k) ) ) CYCLE
       IF( k == 1 ) THEN
         CALL Info('MagnetoDynamicsCalcFields','Harmonic Loss Linear by bodies',Level=6)
       ELSE IF( k == 2 ) THEN
         CALL Info('MagnetoDynamicsCalcFields','Harmonic Loss Quadratic by bodies',Level=6)
       ELSE
         CALL Info('MagnetoDynamicsCalcFields','Joule Loss by bodies',Level=6)
       END IF

       DO j=1,Model % NumberOfBodies
         IF( BodyLoss(k,j) < TINY( TotalLoss(k) ) ) CYCLE
         WRITE( Message,'(A,I0,A,ES15.6)') 'Body ',j,' : ',BodyLoss(k,j)
         CALL Info('MagnetoDynamicsCalcFields', Message, Level=6 )
       END DO

       ! Save losses to components if requested. 
       !---------------------------------------------------------------------------
       DO j=1,Model % NumberOfComponents
         CompParams => Model % Components(j) % Values
         IF( ListGetLogical( CompParams,'Calculate Magnetic Losses', Found ) ) THEN
           MasterBodies => ListGetIntegerArray( CompParams,'Master Bodies',Found ) 
           IF(.NOT. Found ) CYCLE
           val = SUM( BodyLoss(k,MasterBodies) )
           IF( k == 1 ) THEN
             CALL ListAddConstReal( CompParams,'res: harmonic loss linear',val )
           ELSE IF( k == 2 ) THEN
             CALL ListAddConstReal( CompParams,'res: harmonic loss quadratic',val )
           ELSE
             CALL ListAddConstReal( CompParams,'res: joule loss',val )
           END IF
         END IF
       END DO
     END DO

     IF( ParEnv % MyPe == 0 ) THEN
       LossFile = ListGetString(SolverParams,'Harmonic Loss Filename',Found )
       IF( Found ) THEN
         OPEN(NEWUNIT=IOUnit, FILE=LossFile)
         WRITE(IOUnit,'(A)')  '!body_id   harmonic(1)      harmonic(2)      joule'
         DO j=1,Model % NumberOfBodies
           IF( SUM(BodyLoss(1:3,j)) < TINY( TotalLoss(1) ) ) CYCLE
           WRITE(IOUnit,'(I0,T10,3ES17.9)') j, BodyLoss(1:3,j)
         END DO
         CALL Info('MagnetoDynamicsCalsFields', &
             'Harmonic loss for bodies was saved to file: '//TRIM(LossFile),Level=6 )
         CLOSE(IOUnit)
       END IF
     END IF

     DEALLOCATE( BodyLoss )      
   END IF

   IF (GetLogical(SolverParams,'Show Angular Frequency',Found)) THEN
     WRITE(Message,*) 'Angular Frequency: ', Omega
     CALL Info( 'MagnetoDynamicsCalcFields', Message )
     CALL ListAddConstReal(Model % Simulation,'res: Angular Frequency', Omega)
   END IF
   
   IF(ASSOCIATED(NF)) THEN
     CALL NodalTorque(Torque, TorqueGroups)
     DO i=1,SIZE(TorqueGroups)
       j = TorqueGroups(i)
       WRITE (Message,'(A)') 'res: Group '//I2S(j)//' torque'
       CALL ListAddConstReal(Model % Simulation, TRIM(Message), Torque(i))
       WRITE (Message,'(A,F0.8)') 'Torque Group '//I2S(j)//' torque:', Torque(i)
       CALL Info( 'MagnetoDynamicsCalcFields', Message)
     END DO

     IF( ListGetLogicalAnyBody( Model,'Calculate Torque over body') ) THEN
       CALL Fatal( 'MagnetoDynamicsCalcFields', &
           'Keyword "Calculate Torque over body" is deprecated, use Component with torque instead')
     END IF
   END IF

   
  ! Flux On Boundary:
  !------------------

  CalcFluxLogical = .FALSE.
  Flux = 0._dp
  Area = 0._dp
  AverageFluxDensity = 0._dp

  IF (ListGetLogicalAnyBC( Model,'Magnetic Flux Average')) THEN
    IF (PiolaVersion) THEN
      CALL Warn('MagnetoDynamicsCalcFields', &
          'Magnetic Flux Average: The feature is not yet available for Piola transformed basis functions')
    ELSE
      DO i=1,GetNOFBoundaryElements()
        Element => GetBoundaryElement(i)
        BC=>GetBC()
        IF (.NOT. ASSOCIATED(BC) ) CYCLE

        SELECT CASE(GetElementFamily())
        CASE(1)
          CYCLE
        CASE(2)
          k = GetBoundaryEdgeIndex(Element,1); Element => Mesh % Edges(k)
        CASE(3,4)
          k = GetBoundaryFaceIndex(Element)  ; Element => Mesh % Faces(k)
        END SELECT
        IF (.NOT. ActiveBoundaryElement(Element)) CYCLE

        IF (ASSOCIATED(Element % BoundaryInfo % Right)) THEN
          BodyId = Element % BoundaryInfo % Right % BodyID       
        ELSE IF (ASSOCIATED(Element % BoundaryInfo % Left)) THEN
          BodyId = Element % BoundaryInfo % Left % BodyID
        ELSE 
          CALL Fatal ('MagnetoDynamicsCalcFields', 'Magnetic Flux Average: Boundary Element has not got a parent element.')
        END IF

        n = GetElementNOFNodes()
        np = n*pSolver % Def_Dofs(GetElementFamily(Element),BodyId,1)
        nd = GetElementNOFDOFs(uElement=Element, uSolver=pSolver)
        CALL GetVectorLocalSolution(SOL,Pname,uElement=Element,uSolver=pSolver)

        CalcFluxLogical = GetLogical( BC, 'Magnetic Flux Average', Found)
        IF (Found .AND. CalcFluxLogical) CALL calcAverageFlux(Flux, Area, Element, n, nd, np, SOL, vDOFs)
      END DO

      IF( Parallel ) THEN
        Flux(1) = ParallelReduction(Flux(1)) / NoSlices
        Flux(2) = ParallelReduction(Flux(2)) / NoSlices
        Area = ParallelReduction(Area) / NoSlices
      END IF

      IF( Area < EPSILON( Area ) ) THEN
        CALL WARN('MagnetoDynamicsCalcFields', 'Magnetic Flux Average Computation: Area < Epsilon(Area)')
        RETURN
      END IF

      AverageFluxDensity = Flux / Area
 
    WRITE(Message,*) 'Magnetic Flux Average: ', Flux(1)
    CALL Info( 'MagnetoDynamicsCalcFields', Message )
    CALL ListAddConstReal( Model % Simulation, 'res: Magnetic Flux Average', Flux(1) )
 
    IF (vDOFs == 2) THEN 
      WRITE(Message,*) 'Magnetic Flux im Average: ', Flux(2)
      CALL Info( 'MagnetoDynamicsCalcFields', Message )
      CALL ListAddConstReal( Model % Simulation, 'res: Magnetic Flux im Average', Flux(2) )
    END IF

    WRITE(Message,*) 'Magnetic Flux Density Average: ', AverageFluxDensity(1)
    CALL Info( 'MagnetoDynamicsCalcFields', Message )
    CALL ListAddConstReal( Model % Simulation,'res: Magnetic Flux Density Average', &
                          AverageFluxDensity(1))

    IF (vDOFs == 2) THEN 
     WRITE(Message,*) 'Magnetic Flux Density im Average: ', AverageFluxDensity(2)
     CALL Info( 'MagnetoDynamicsCalcFields', Message )
     CALL ListAddConstReal( Model % Simulation,'res: Magnetic Flux Density im Average', &
                          AverageFluxDensity(2))
    END IF

    WRITE(Message,*) 'Magnetic Flux Area: ', Area
    CALL Info( 'MagnetoDynamicsCalcFields', Message )
    CALL ListAddConstReal( Model % Simulation,'res: Magnetic Flux Area', Area )
    END IF
  END IF

  
  HasThinLines = ListCheckPresentAnyBC(Model,'Thin Line Crossection Area')
  IF(HasThinLines) THEN
    ThinLinePower = 0._dp
    ALLOCATE(ThinLineCrossect(n), ThinLineCond(n))
    Active = GetNOFBoundaryElements()
    DO i=1,Active
       Element => GetBoundaryElement(i)
       BC=>GetBC()
       IF (.NOT. ASSOCIATED(BC) ) CYCLE
       
       ThinLineCrossect = GetReal( BC, 'Thin Line Crossection Area', Found)

       IF (Found) THEN
         CALL Info("CalcFields", "Found a Thin Line Element", level=10)
         ThinLineCond = GetReal(BC, 'Thin Line Conductivity', Found)
         IF (.NOT. Found) CALL Fatal('CalcFields','Thin Line Conductivity not found!')
       ELSE
         CYCLE
       END IF

       SELECT CASE(GetElementFamily())
       CASE(1)
         CYCLE
       CASE(2)
         k = GetBoundaryEdgeIndex(Element,1); Element => Mesh % Edges(k)
       CASE(3,4)
         CYCLE
       END SELECT
       IF (.NOT. ActiveBoundaryElement(Element)) CYCLE


       Model % CurrentElement => Element
       nd = GetElementNOFDOFs(Element)
       n  = GetElementNOFNodes(Element)
       CALL GetElementNodes(Nodes, Element)
  !     line_tangent = 0._dp
  !     line_tangent(1) = Nodes % x(2) - Nodes % x(1)
  !     line_tangent(2) = Nodes % y(2) - Nodes % y(1)
  !     line_tangent(3) = Nodes % z(2) - Nodes % z(1)
  !     line_tangent(:) = line_tangent(:) / SQRT(SUM(line_tangent(:)**2._dp))

       CALL GetVectorLocalSolution(SOL, Pname, uElement=Element, uSolver=pSolver)
       IF (Transient) THEN 
         CALL GetScalarLocalSolution(PSOL,Pname,uSolver=pSolver,Tstep=-1)
         PSOL(1:nd)=(SOL(1,1:nd)-PSOL(1:nd))/dt
       END IF

      ! Numerical integration:
      !-----------------------
      IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree)

      np = n*MAXVAL(Solver % Def_Dofs(GetElementFamily(Element),:,1))

      DO j=1,IP % n

         stat = EdgeElementInfo(Element, Nodes, IP % U(j), IP % V(j), &
             IP % W(j), DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, &
             dBasisdx = dBasisdx, BasisDegree = EdgeBasisDegree, &
             ApplyPiolaTransform = .TRUE.)
   
         C_ip  = SUM(Basis(1:n) * ThinLineCond(1:n))
         Area = SUM(Basis(1:n) * ThinLineCrossect(1:n))
         s = detJ*IP % s(j)

         IF (vDOFS == 1) THEN
           !da/dt part
           IF (Transient) THEN
             E(1,:) = -MATMUL(PSOL(np+1:nd), Wbasis(1:nd-np,:))
           END IF

           !grad V part
           E(1,:) = E(1,:)-MATMUL(SOL(1,1:np), dBasisdx(1:np,:))

           ! The Joule heating power per unit volume: J.E = (sigma * E).E
  !         Coeff = Area * C_ip * SUM(line_tangent(:) * E(1,:)) ** 2._dp * s
           Coeff = Area * C_ip * SUM(E(1,:) ** 2._dp) * s
         ELSE
           !da/dt part
           E(1,:) = Omega*MATMUL(SOL(2,np+1:nd),WBasis(1:nd-np,:))
           E(2,:) = -Omega*MATMUL(SOL(1,np+1:nd),WBasis(1:nd-np,:))

           !grad V part
           E(1,:) = E(1,:)-MATMUL(SOL(1,1:np), dBasisdx(1:np,:))
           E(2,:) = E(2,:)-MATMUL(SOL(2,1:np), dBasisdx(1:np,:))
           CALL Warn('CalcFields', 'Power loss not implemented for harmonic case')
           Coeff = 0._dp
           ! Now Power = J.conjugate(E), with the possible imaginary component neglected.         
           !Coeff = HarmPowerCoeff * (SUM( MATMUL( REAL(CMat_ip(1:3,1:3)), TRANSPOSE(E(1:1,1:3)) ) * &
           !    TRANSPOSE(E(1:1,1:3)) ) * Basis(p) * s - &
           !    SUM( MATMUL( AIMAG(CMat_ip(1:3,1:3)), TRANSPOSE(E(2:2,1:3)) ) * &
           !    TRANSPOSE(E(1:1,1:3)) ) * Basis(p) * s + &
           !    SUM( MATMUL( AIMAG(CMat_ip(1:3,1:3)), TRANSPOSE(E(1:1,1:3)) ) * &
           !    TRANSPOSE(E(2:2,1:3)) ) * Basis(p) * s + &               
           !    SUM( MATMUL( REAL(CMat_ip(1:3,1:3)), TRANSPOSE(E(2:2,1:3)) ) * &
           !    TRANSPOSE(E(2:2,1:3)) ) * Basis(p) * s)
         END IF

         ThinLinePower = ThinLinePower + Coeff

      END DO
    END DO

    
    IF( Parallel ) THEN
      ThinLinePower  = ParallelReduction(ThinLinePower) / NoSlices
    END IF
    WRITE(Message,*) 'Total thin line power (the Joule effect): ', ThinLinePower
    CALL Info( 'MagnetoDynamicsCalcFields', Message )
    CALL ListAddConstReal(Model % Simulation, 'res: thin line power', ThinLinePower)

    DEALLOCATE(ThinLineCrossect, ThinLineCond)
  END IF

  !
  ! postprocessing of surface currents generated by Thin sheet BCs (in time-harmonic cases): 
  !

  Found = ListCheckPresentAnyBC(Model,'Thin Sheet Thickness')
  IF(Found) THEN
    ALLOCATE(SheetThickness(n), ThinLineCond(n))
    Power = 0.0_dp
    DO i=1,GetNOFBoundaryElements()
      Element => GetBoundaryElement(i)
      BC => GetBC()
      IF (.NOT. ASSOCIATED(BC)) CYCLE

      SheetThickness = GetConstReal( BC, 'Thin Sheet Thickness', Found)
      IF (Found) THEN
        ! Technically, there is no skin but why create yet another conductivity variable?
        ThinLineCond = GetConstReal(BC, 'Thin Sheet Electric Conductivity', Found)
        IF (.NOT. Found) ThinLineCond = 1.0_dp ! if not found default to "air" property
      ELSE
        CYCLE
      END IF
      
      SELECT CASE(GetElementFamily())
      CASE(1)
        CYCLE
      CASE(2)
        k = GetBoundaryEdgeIndex(Element,1)
        Element => Mesh % Edges(k)
      CASE(3,4)
        k = GetBoundaryFaceIndex(Element)
        Element => Mesh % Faces(k)
      END SELECT
      IF (.NOT. ActiveBoundaryElement(Element)) CYCLE

      Model % CurrentElement => Element

      nd = GetElementNOFDOFs(Element)
      n  = GetElementNOFNodes(Element)
      CALL GetElementNodes(Nodes, Element)

      CALL GetLocalSolution(SOL, Pname, uElement=Element, uSolver=pSolver)
      IF (Transient) THEN 
        CALL GetLocalSolution(PSOL,Pname,uSolver=pSolver,Tstep=-1)
        PSOL(1:nd)=(SOL(1,1:nd)-PSOL(1:nd))/dt
      END IF

      ! Numerical integration:
      !-----------------------
      IP = GaussPoints(Element)

      np = n*MAXVAL(Solver % Def_Dofs(GetElementFamily(Element),:,1))

      DO j=1,IP % n
         stat = EdgeElementInfo( Element, Nodes, IP % U(j), IP % V(j), &
              IP % W(j), DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, &
              dBasisdx = dBasisdx, BasisDegree = 1, &
              ApplyPiolaTransform = .TRUE.)
   
         C_ip  = SUM(Basis(1:n) * ThinLineCond(1:n))
         localThickness = SUM(Basis(1:n) * SheetThickness(1:n))
         s = detJ*IP % s(j)

         IF (vDOFS == 1) THEN
           !da/dt part
           IF (Transient) THEN
             E(1,:) = -MATMUL(PSOL(np+1:nd), Wbasis(1:nd-np,:))
           END IF

           !grad V part
           E(1,:) = E(1,:)-MATMUL(SOL(1,1:np), dBasisdx(1:np,:))

           Coeff = localThickness * C_ip * SUM(E(1,:) ** 2._dp) * s
         ELSE
           !da/dt part
           E(1,:) = Omega*MATMUL(SOL(2,np+1:nd),WBasis(1:nd-np,:))
           E(2,:) = -Omega*MATMUL(SOL(1,np+1:nd),WBasis(1:nd-np,:))

           !grad V part
           E(1,:) = E(1,:)-MATMUL(SOL(1,1:np), dBasisdx(1:np,:))
           E(2,:) = E(2,:)-MATMUL(SOL(2,1:np), dBasisdx(1:np,:))
           ! Now Power = J.conjugate(E), with the possible imaginary component neglected.         
           Coeff = localThickness * HarmPowerCoeff * ( &
             SUM(C_ip * E(1,:) ** 2._dp ) +     &
             SUM(C_ip * E(2,:) ** 2._dp )       &
             ) * s
         END IF

         Power = Power + Coeff

      END DO
    END DO
    
    IF( Parallel ) THEN
      Power = ParallelReduction(Power) / NoSlices
    END IF
      
    WRITE(Message,*) 'Thin sheet current power (the Joule effect): ', Power
    CALL Info( 'MagnetoDynamicsCalcFields', Message )
    CALL ListAddConstReal(Model % Simulation, 'res: Thin sheet current power', Power)

  END IF

  IF (EigenAnalysis .AND. .NOT. VtuStyle ) THEN
    ! The primary variable with size 6 is problematic for Vtu output.
    ! Hence we operate directly on the eigenvectors if eigenmodes are active.
    DO i=1,MaxFields
      FieldVariable => NodalFieldPointers(i) % Field
      IF (ASSOCIATED(FieldVariable)) THEN
        n = SIZE(FieldVariable % Values)
        DO j=1,n/2
          FieldVariable % EigenVectors(field,j) = CMPLX( &
              FieldVariable % Values(2*j-1), FieldVariable % Values(2*j), KIND=dp )
        END DO
      END IF

      FieldVariable => ElementalFieldPointers(i) % Field
      IF (ASSOCIATED(FieldVariable)) THEN
        n = SIZE(FieldVariable % Values)
        DO j=1,n/2
          FieldVariable % EigenVectors(field,j) = CMPLX( &
              FieldVariable % Values(2*j-1), FieldVariable % Values(2*j), KIND=dp )
        END DO
      END IF
    END DO
  END IF

  END DO COMPUTE_FIELDS
  
  IF( NormIndex > 0 ) THEN
    WRITE(Message,*) 'Reverting norm to: ', SaveNorm
    CALL Info( 'MagnetoDynamicsCalcFields', Message )
    Solver % Variable % Norm = SaveNorm
  END IF

  IF(.NOT. ConstantMassMatrixInUse ) THEN
    ConstantMassMatrixInUse = ListGetLogical( SolverParams,'Constant Mass Matrix',Found )    
  END IF

CONTAINS

  
!-------------------------------------------------------------------
  SUBROUTINE SumElementalVariable(Var, Values, BodyId, uAdditive)
!-------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Variable_t), POINTER :: Var
    REAL(KIND=dp), OPTIONAL, TARGET :: Values(:)
    INTEGER, OPTIONAL :: BodyId
    LOGICAL, OPTIONAL :: uAdditive

    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), ALLOCATABLE :: NodeSum(:)
    INTEGER :: n, j, k, l, nodeind, dgind, bias
    LOGICAL, ALLOCATABLE :: AirGapNode(:)
    LOGICAL :: Additive
    REAL(KIND=dp), POINTER :: ValuesSource(:)


    Additive = .FALSE.
    IF(PRESENT(uAdditive)) Additive = uAdditive

    IF(PRESENT(Values)) THEN
      ValuesSource => Values
    ELSE 
      IF( .NOT. ASSOCIATED( Var ) ) RETURN
      ValuesSource => Var % Values
    END IF

    n = Mesh % NumberOFNodes
    ALLOCATE(NodeSum(n), AirGapNode(n))
    AirGapNode = .FALSE.

    ! Collect nodal sum of DG elements
    DO k=1,Var % Dofs
      NodeSum = 0.0_dp

      ! Collect DG data to nodal vector
      DO j=1, Mesh % NumberOfBulkElements
        Element => Mesh % Elements(j)
        IF(PRESENT(BodyID)) THEN
          IF(Element % BodyID /= BodyID) CYCLE
        END IF
        DO l = 1, Element % TYPE % NumberOfNodes
          nodeind = Element % NodeIndexes(l)
          dgind = Var % Perm(Element % DGIndexes(l))
          IF( dgind > 0 ) THEN
            NodeSum( nodeind ) = NodeSum( nodeind ) + &
              ValuesSource(Var % DOFs*( dgind-1)+k ) 
          END IF 
        END DO
      END DO

      ! Sum nodal data to elements
      DO j=1, Mesh % NumberOfBulkElements
        Element => Mesh % Elements(j)
        IF(PRESENT(BodyID)) THEN
          IF(Element % BodyID /= BodyID) CYCLE
        END IF
        DO l=1,Element%TYPE%NumberofNodes
          nodeind = Element % NodeIndexes(l)
          dgind = Var % Perm(Element % DGIndexes(l))
          IF( dgind > 0 ) THEN
            IF( Additive) THEN
              Var % Values( var % DOFs*(dgind-1)+k) = NodeSum(nodeind) + &
                Var % Values( var % DOFs*(dgind-1)+k)
            ELSE
              Var % Values( var % DOFs*(dgind-1)+k) = NodeSum(nodeind)
            END IF
          END IF
        END DO
      END DO

    END DO
!-------------------------------------------------------------------
  END SUBROUTINE SumElementalVariable
!-------------------------------------------------------------------


!-------------------------------------------------------------------
  SUBROUTINE CalcBoundaryModels( )
!-------------------------------------------------------------------
    IMPLICIT NONE
!-------------------------------------------------------------------
    REAL(KIND=dp) :: GapLength(27), AirGapMu(27)

!-------------------------------------------------------------------
    LOGICAL :: FirstTime = .TRUE.
    REAL(KIND=dp) :: B2, GapLength_ip, LeftCenter(3), &
      RightCenter(3), BndCenter(3), LeftNormal(3), RightNormal(3), &
      NF_ip_l(27,3), NF_ip_r(27,3), xcoord
    TYPE(Element_t), POINTER :: LeftParent, RightParent, BElement
    TYPE(Nodes_t), SAVE :: LPNodes, RPNodes
!    REAL(KIND=dp) :: F(3,3)
    INTEGER :: n_lp, n_rp, LeftBodyID, RightBodyID
    REAL(KIND=dp), ALLOCATABLE :: LeftFORCE(:,:), RightFORCE(:,:), &
      AirGapForce(:,:), ForceValues(:)
    INTEGER, ALLOCATABLE :: RightMap(:), LeftMap(:)
    REAL(KIND=dp) :: ParentNodalU(n), parentNodalV(n), ParentNodalW(n)
    REAL(KIND=dp) :: Normal(3)
    REAL(KIND=dp), SAVE :: mu0 = 1.2566370614359173e-6_dp
    LOGICAL, ALLOCATABLE :: BodyMask(:)
    LOGICAL :: HasLeft, HasRight

    n = Mesh % MaxElementDOFs

    ALLOCATE(LeftFORCE(n,3), RightForce(n,3), RightMap(n), LeftMap(n), &
      AirGapForce(3,Mesh % NumberOfNodes) )

    IF ( FirstTime ) THEN
      mu0 = GetConstReal(CurrentModel % Constants, &
        'Permeability of Vacuum', Found)
      IF(.NOT. Found) mu0 = 1.2566370614359173e-6
    END IF

    LeftBodyID = -1
    RightBodyID = -1

    DO i = 1, GetNOFBoundaryElements()
      BElement => GetBoundaryElement(i, uSolver=pSolver)
      BC => GetBC(BElement)
      IF (.NOT. ASSOCIATED(BC) ) CYCLE

      n = GetElementNOFNodes(BElement)
      GapLength(1:n) = GetReal(BC, 'Air Gap Length', Found)
      IF(.NOT. Found) CYCLE

      IF( .NOT. ASSOCIATED( BElement % BoundaryInfo ) ) CYCLE
      
      HasLeft = ASSOCIATED(BElement % BoundaryInfo % Left)
      HasRight = ASSOCIATED(BElement % BoundaryInfo % Right)
      IF( .NOT. (HasLeft .OR. HasRight)) THEN
        CALL Warn('MagnetoDynamicsCalcFields', 'Airgap Length given on orphan boundary')
        CYCLE
      END IF

      IF(.NOT. (HasLeft .AND. HasRight)) &
        CALL Warn('MagnetoDynamicsCalcFields', 'Onesided airgap force calculation is untested.')

      BElement => Mesh % Faces(GetBoundaryFaceIndex(BElement))
      IF(.NOT. ActiveBoundaryElement(BElement, uSolver=pSolver)) CYCLE

      LeftBodyID = BElement % BoundaryInfo % Left % BodyID
      RightBodyID = BElement % BoundaryInfo % Right % BodyID
      IF(LeftBodyID == RightBodyID) THEN
        CALL Warn('MagnetoDynamicsCalcFields', 'Airgap in the middle of single body Id')
        CYCLE
      END IF

      IF(HasLeft) LeftParent => BElement % BoundaryInfo % Left
      IF(HasRight) RightParent => BElement % BoundaryInfo % Right

      IF(HasLeft) n_lp = GetElementNOFNodes(LeftParent)
      if(HasRight) n_rp = GetElementNOFNodes(RightParent) 

      CALL GetElementNodes(Nodes, BElement)
      IF(HasLeft) CALL GetElementNodes(LPNodes, LeftParent)
      IF(HasRight) CALL GetElementNodes(RPNodes, RightParent)

      CALL GetVectorLocalSolution(SOL,Pname,uElement=BElement, uSolver=pSolver)

      IF(HasLeft) LeftCenter(1:3) = &
          [ SUM(LPNodes % x(1:n_lp)), SUM(LPNodes % y(1:n_lp)), SUM(LPNodes % z(1:n_lp)) ] / n_lp
      IF(HasRight) RightCenter(1:3) = &
          [ SUM(RPNodes % x(1:n_rp)), SUM(RPNodes % y(1:n_rp)), SUM(RPNodes % z(1:n_rp)) ] / n_rp
      BndCenter(1:3) = [ sum(Nodes % x(1:n)), sum(Nodes % y(1:n)), sum(Nodes % z(1:n)) ] / n

      np = n*MAXVAL(pSolver % Def_Dofs(GetElementFamily(BElement),:,1))
      nd = GetElementNOFDOFs(uElement=BElement, uSolver=pSolver)

      
      DO k = 1,n
        IF(HasLeft) THEN  
          DO l = 1,n_lp
            IF(LeftParent % NodeIndexes(l) == BElement % NodeIndexes(k)) LeftMap(k) = l
          END DO
        END IF
        IF(HasRight) THEN
          DO l = 1,n_rp
            IF(RightParent % NodeIndexes(l) == BElement % NodeIndexes(k)) RightMap(k) = l
          END DO
        END IF
      END DO

      AirGapMu(1:n) = GetReal(BC, 'Air Gap Relative Permeability', Found)
      IF(.NOT. Found) AirGapMu(1:n) = 1.0_dp

      LeftFORCE = 0.0_dp
      RightFORCE = 0.0_dp

      IP = GaussPoints(BElement, EdgeBasis=(dim==3), PReferenceElement=PiolaVersion, &
          EdgeBasisDegree=EdgeBasisDegree)
      
      DO j = 1,IP % n
        stat = ElementInfo( BElement, Nodes, IP % U(j), IP % V(j), &
            IP % W(j), detJ, Basis, dBasisdx, &
            EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver ) 

        R_ip = SUM( Basis(1:n)/(mu0*AirGapMu(1:n)) )
        GapLength_ip = SUM( Basis(1:n)*GapLength(1:n) )

        s = detJ * IP % s(j)        
        IF ( CSymmetry ) THEN
          xcoord = SUM( Basis(1:n) * Nodes % x(1:n) )
          s = s * xcoord 
        END IF
          
        Normal = NormalVector(BElement, Nodes, IP% U(j), IP % V(j))
        IF(HasLeft)  THEN
          IF( SUM(normal*(LeftCenter - bndcenter)) >= 0 ) THEN
            LeftNormal = -Normal
          ELSE
            LeftNormal = Normal
          END IF
        END IF

        IF(HasRight) THEN
          IF( SUM(normal*(RightCenter - bndcenter)) >= 0 ) THEN
            RightNormal = -Normal
          ELSE
            RightNormal = Normal
          END IF
        END IF


        DO k=1,vDOFs
          SELECT CASE(dim)
          CASE(2)
            ! This has been done with the same sign convention as in MagnetoDynamics2D:
            ! -------------------------------------------------------------------------
            IF ( CSymmetry ) THEN
              B(k,1) = -SUM( SOL(k,1:nd) * dBasisdx(1:nd,2) )
              B(k,2) = SUM( SOL(k,1:nd) * dBasisdx(1:nd,1) ) &
                + SUM( SOL(k,1:nd) * Basis(1:nd) ) / xcoord
              B(k,3) = 0._dp
            ELSE
              B(k,1) =  SUM( SOL(k,1:nd) * dBasisdx(1:nd,2) )
              B(k,2) = -SUM( SOL(k,1:nd) * dBasisdx(1:nd,1) )
              B(k,3) = 0._dp
            END IF
          CASE(3)
            IF (PiolaVersion) THEN
              B(k,:) = normal*SUM( SOL(k,np+1:nd) * RotWBasis(1:nd-np,3) )
            ELSE
              B(k,:) = normal*SUM(Normal(1:3) * MATMUL( SOL(k,np+1:nd), RotWBasis(1:nd-np,1:3) ))
            END IF
          END SELECT
        END DO



        B2 = SUM(B(1,:)*B(1,:) + B(2,:)*B(2,:))
        IF (ASSOCIATED(NF).OR.ASSOCIATED(EL_NF)) THEN
          NF_ip_r = 0._dp
          NF_ip_l = 0._dp
          DO k=1,n
            DO l=1,3
              DO m=1,3
                IF(HasLeft)  NF_ip_l(k,l) = NF_ip_l(k,l) + R_ip*B(1,l)*B(1,m)*(LeftNormal(m)*Basis(k))
                IF(HasRight) NF_ip_r(k,l) = NF_ip_r(k,l) + R_ip*B(1,l)*B(1,m)*(RightNormal(m)*Basis(k))
              END DO
              IF(HasLeft) NF_ip_l(k,l) = NF_ip_l(k,l) - 0.5*R_ip*B2*(LeftNormal(l)*Basis(k))
              IF(HasRight) NF_ip_r(k,l) = NF_ip_r(k,l) - 0.5*R_ip*B2*(RightNormal(l)*Basis(k))
            END DO
          END DO
        END IF

        Energy(2) = Energy(2) + GapLength_ip*s*0.5*R_ip*B2

        DO p=1,n
          IF(HasLeft) LeftFORCE(LeftMap(p), 1:3) = LeftFORCE(LeftMap(p), 1:3) + s*NF_ip_l(p,1:3)
          IF(HasRight) RightFORCE(RightMap(p), 1:3) = RightFORCE(RightMap(p), 1:3) + s*NF_ip_r(p,1:3)
        END DO
      END DO ! Integration points
      
      IF(ElementalFields) THEN
        IF(HasLeft) CALL LocalCopy(EL_NF, fdim, n_lp, LeftFORCE, 0, UElement=LeftParent, uAdditive=.TRUE.)
        IF(HasRight) CALL LocalCopy(EL_NF, fdim, n_rp, RightFORCE, 0, UElement=RightParent, uAdditive=.TRUE.)
      END IF
    END DO ! Boundary elements
    
    DEALLOCATE( LeftFORCE, RightFORCE, RightMap, LeftMap, AirGapForce )

    IF (HomogenizationLoss) DEALLOCATE (Nu_el)
!-------------------------------------------------------------------
  END SUBROUTINE CalcBoundaryModels
!-------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE NodalTorque(T, TorqueGroups)
!------------------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, ALLOCATABLE, INTENT(OUT) :: TorqueGroups(:)
   REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: T(:)
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------

   REAL(KIND=dp), POINTER :: origins(:,:), omegas(:,:)
   TYPE(Element_t), POINTER :: Element
   TYPE(Variable_t), POINTER :: CoordVar
   TYPE(ValueList_t), POINTER :: BodyParams, SolverParams
   TYPE(BodyArray_t), POINTER :: bodies(:)
   INTEGER, POINTER :: LocalGroups(:)
   REAL(KIND=dp), ALLOCATABLE :: axes(:,:)

   LOGICAL, ALLOCATABLE :: VisitedNode(:,:)
   INTEGER, ALLOCATABLE :: AllGroups(:)

   REAL(KIND=dp) :: origin(3), axisvector(3), P(3), F(3), v1(3), v2(3), nrm
   INTEGER :: pnodal, nnt, ElemNodeDofs(27), ndofs, globalnode, ng, ngroups, &
     n, maxngroups, m, k, num_origins, num_axes, pivot
   LOGICAL :: Found

   ! make union of body-wise declared torque groups. \TODO: make this abstract and move to generalutils
   bodies => Model % bodies
   maxngroups = 0
   DO n=1,size(bodies)
     LocalGroups => ListGetIntegerArray(bodies(n) % Values, "Torque Groups", Found)
     IF (Found) THEN
       maxngroups = maxngroups + size(LocalGroups)
     END IF
   END DO
   IF(maxngroups .eq. 0) THEN
     ALLOCATE(TorqueGroups(0), T(0))
     RETURN
   END IF

   ALLOCATE(AllGroups(maxngroups))
   AllGroups = -1
   ngroups = 0
   DO n=1,size(bodies)
     LocalGroups => ListGetIntegerArray(bodies(n) % Values, "Torque Groups", Found)
     IF (Found) THEN
       AllGroups((ngroups+1):(ngroups+size(LocalGroups))) = LocalGroups(1:size(LocalGroups))
       ngroups = ngroups + size(LocalGroups)
     END IF
   END DO
   call SORT(size(AllGroups), AllGroups)
   pivot = AllGroups(1)
   k = 1
   m = 1
   do while(pivot .ne. -1)
     AllGroups(k) = pivot
     DO n=m,size(AllGroups)
       IF (AllGroups(k) .ne. AllGroups(n)) then
         pivot = AllGroups(n)
         k = k + 1
         m = n
         exit
       end if
       pivot = -1
     END DO
   END DO
   ALLOCATE(TorqueGroups(k))
   IF(k == 0) RETURN

   TorqueGroups = AllGroups(1:k)
   ! done making union

   SolverParams => GetSolverParams()
   origins => ListGetConstRealArray(SolverParams, "Torque Group Origins", Found)
   IF (.NOT. Found) THEN
     num_origins = 0
   ELSE
     num_origins = SIZE(origins,1)
   END IF

   omegas => ListGetConstRealArray(SolverParams, "Torque Group Axes", Found)
   IF (.NOT. Found) THEN
     num_axes = 0
   ELSE
     num_axes = SIZE(omegas,1)
     ALLOCATE(axes(num_axes, size(omegas, 2)))
     axes = omegas
     DO k = 1, num_axes
       nrm = sqrt(sum(axes(k,:)*axes(k,:))) 
       IF (nrm == 0._dp) THEN
         CALL Warn('MagnetoDynamicsCalcFields',&
             'Axis for the torque group '//I2S(k)//' is a zero vector')
         CYCLE
       END IF
       axes(k,:) = axes(k,:) / nrm
     END DO
   END IF

   ng = size(TorqueGroups,1)
   ALLOCATE(T(ng*3))
   ALLOCATE(VisitedNode(Mesh % NumberOfNodes, ng))
   VisitedNode = .FALSE.
   T = 0._dp

   DO pnodal=1,GetNOFActive()
     Element => GetActiveElement(pnodal)
     BodyParams => GetBodyParams(Element)
     LocalGroups => ListGetIntegerArray(BodyParams, "Torque Groups", Found)
     IF(.NOT. Found) CYCLE
     ndofs = GetElementDOFs(ElemNodeDofs)
     DO nnt=1,ndofs
       globalnode = ElemNodeDofs(nnt)
       F(1) = NF % Values( 3*(NF % Perm((globalnode))-1) + 1)
       F(2) = NF % Values( 3*(NF % Perm((globalnode))-1) + 2)
       F(3) = NF % Values( 3*(NF % Perm((globalnode))-1) + 3)
       P(1) = Mesh % Nodes % x(globalnode)
       P(2) = Mesh % Nodes % y(globalnode)
       P(3) = Mesh % Nodes % z(globalnode)
       DO ng=1,size(LocalGroups)
         IF (.NOT. VisitedNode(globalnode, LocalGroups(ng))) THEN
           VisitedNode(globalnode, LocalGroups(ng)) = .TRUE.
           IF (LocalGroups(ng) > num_origins) THEN
             origin = 0._dp
           ELSE
             origin = origins(LocalGroups(ng),1:3)
           END IF
           IF (LocalGroups(ng) > num_axes) THEN
             axisvector = 0._dp
             axisvector(3) = 1._dp
           ELSE
             axisvector = axes(LocalGroups(ng), 1:3)
           END IF
           v1 = P - origin
           v1 = (1 - SUM(axisvector*v1))*v1
           v2 = CrossProduct(v1,F)
           T(LocalGroups(ng)) = T(LocalGroups(ng)) + sum(axisvector*v2)
         END IF
       END DO 
     END DO 

   END DO

   IF( Parallel ) THEN
     DO ng=1,SIZE(TorqueGroups)
       T(ng) = ParallelReduction(T(ng)) / NoSlices
     END DO
   END IF

!------------------------------------------------------------------------------
 END SUBROUTINE NodalTorque
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GlobalSol(Var, m, b, dofs,EL_Var )
!------------------------------------------------------------------------------
   IMPLICIT NONE
   REAL(KIND=dp), TARGET CONTIG :: b(:,:)
   INTEGER :: m, dofs
   TYPE(Variable_t), POINTER :: Var
   TYPE(Variable_t), POINTER, OPTIONAL :: EL_Var
!------------------------------------------------------------------------------
   INTEGER :: i
!------------------------------------------------------------------------------
   COMPLEX(KIND=dp), POINTER :: EigVec(:)
   INTEGER :: ic
!------------------------------------------------------------------------------

   IF(PRESENT(EL_Var)) THEN
     IF(ASSOCIATED(El_Var)) THEN
       El_Var % DgAveraged = .FALSE.
       IF( DoAve ) THEN
         CALL Info('MagnetoDynamicsCalcFields','Averaging for field: '//TRIM(El_Var % Name),Level=10)
         CALL CalculateBodyAverage(Mesh, El_Var, .FALSE.)              
       END IF
       IF(.NOT. (ASSOCIATED(var) .AND. NodalFields) ) THEN
         dofs = dofs+m
         RETURN
       END IF
     END IF
   END IF

   IF(.NOT. ASSOCIATED(Var) ) RETURN

   IF(VtuStyle) EigVec => Var % EigenVectors(Field,:)
        
   CALL Info('MagnetoDynamicsCalcFields','Solving for field: '//TRIM(Var % Name),Level=6)   
   DO i=1,m
     dofs = dofs+1

     Solver % Matrix % RHS => b(:,dofs)
     Solver % Variable % Values=0
     Norm = DefaultSolve()

     IF(VtuStyle) THEN
       ic=(i+1)/2
       IF(MODULO(i,2)==1) THEN
         EigVec(ic::m/2) = Solver % Variable % Values
       ELSE
         EigVec(ic::m/2) = CMPLX( REAL(EigVec(ic::m/2)), Solver % Variable % Values )
       END IF
     ELSE
       var % Values(i::m) = Solver % Variable % Values
     END IF
       
     IF( NormIndex == dofs ) SaveNorm = Norm
  END DO
!------------------------------------------------------------------------------
 END SUBROUTINE GlobalSol
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE LocalSol(Var, m, n, nd, A, b, pivot, dofs )
!------------------------------------------------------------------------------
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: Var
   REAL(KIND=dp) :: b(:,:), A(:,:)
   INTEGER :: pivot(:), m,n,nd,dofs
!------------------------------------------------------------------------------
   INTEGER :: ind(n), i
   REAL(KIND=dp) :: x(nd),s
   COMPLEX(KIND=dp), POINTER :: EigVec(:)
   INTEGER :: ic
!------------------------------------------------------------------------------
   IF(.NOT. ASSOCIATED(var)) RETURN

   ind(1:n) = Var % Perm(Element % DGIndexes(1:n))
   IF( ANY( ind(1:n) <= 0 ) ) RETURN

   ind(1:n) = Var % DOFs * (ind(1:n)-1)
   IF(VtuStyle) EigVec => Var % EigenVectors(Field,:)

   DO i=1,m
      dofs = dofs+1
      
      IF( ElementalMode == 2 .OR. ElementalMode == 4 ) THEN
        ! Perform total lumping 
        s = SUM(MASS(1:nd,1:nd))
        x(1:nd) = SUM(b(1:nd,dofs)) / s
      ELSE
        x(1:nd) = b(1:nd,dofs)
        CALL LUSolve(nd,MASS,x,pivot)
      END IF

      IF( VtuStyle ) THEN
        ic = (i+1)/2
          
        IF(MODULO(i,2)==1) THEN
          EigVec(ind(1:n)+ic) = x(1:n)
        ELSE
          EigVec(ind(1:n)+ic) = CMPLX( REAL(EigVec(ind(1:n)+ic)), x(1:n) )
        END IF
      ELSE      
        Var % Values(ind(1:n)+i) = x(1:n)
      END IF
   END DO
!------------------------------------------------------------------------------
 END SUBROUTINE LocalSol
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE LocalCopy(Var, m, n, b, bias, UElement, Values, uAdditive)
!------------------------------------------------------------------------------
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: Var
   INTEGER, INTENT(IN) :: m,n,bias
   INTEGER :: dofs
   REAL(KIND=dp) :: b(:,:)
   TYPE(Element_t), POINTER, OPTIONAL :: UElement
   REAL(KIND=dp), OPTIONAL, TARGET :: Values(:)
   LOGICAL, OPTIONAL :: uAdditive
!------------------------------------------------------------------------------
   INTEGER :: ind(n), i
   LOGICAL :: Additive
   REAL(KIND=dp), POINTER :: PValues(:)
!------------------------------------------------------------------------------
   IF(.NOT. ASSOCIATED(var)) RETURN

   IF( Var % TYPE == Variable_on_nodes_on_elements ) THEN
     IF(PRESENT(UElement)) THEN
       ind(1:n) = Var % Perm(UElement % DGIndexes(1:n))
     ELSE
       ind(1:n) = Var % Perm(Element % DGIndexes(1:n))
     END IF
   ELSE
     IF(PRESENT(UElement)) THEN
       ind(1:n) = Var % Perm(UElement % NodeIndexes(1:n))
     ELSE
       ind(1:n) = Var % Perm(Element % NodeIndexes(1:n))
     END IF
   END IF
     
   IF( ANY(ind(1:n) == 0 ) ) RETURN
   
   ind(1:n) = Var % Dofs * ( ind(1:n) - 1)
   
   IF(PRESENT(uAdditive)) THEN
     Additive = uAdditive
   ELSE
     Additive = .FALSE.
   END IF

   IF( PRESENT( Values ) ) THEN
     PValues => Values
   ELSE
     PValues => Var % Values
   END IF
   
   dofs = bias

   DO i=1,m
     dofs = dofs+1

     IF(Additive) THEN
       PValues(ind(1:n)+i) = PValues(ind(1:n)+i) + b(1:n,dofs)
     ELSE
       PValues(ind(1:n)+i) = b(1:n,dofs)
     END IF
   END DO
!------------------------------------------------------------------------------
 END SUBROUTINE LocalCopy
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE AddLocalFaceTerms(STIFF,FORCE)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     REAL(KIND=dp) :: STIFF(:,:), FORCE(:)

     TYPE(Element_t),POINTER :: P1,P2,Face,Faces(:)
     INTEGER ::t,n,n1,n2,NumberOfFaces,dim

     dim = CoordinateSystemDimension()

     IF (dim==2) THEN
       Faces => Solver % Mesh % Edges
       NumberOfFaces = Solver % Mesh % NumberOfEdges
     ELSE
       Faces => Solver % Mesh % Faces
       NumberOfFaces = Solver % Mesh % NumberOfFaces
     END IF

     DO t=1,NumberOfFaces
       Face => Faces(t)
       IF ( .NOT. ActiveBoundaryElement(Face) ) CYCLE

       P1 => Face % BoundaryInfo % Left
       P2 => Face % BoundaryInfo % Right
       IF ( ASSOCIATED(P2) .AND. ASSOCIATED(P1) ) THEN
          IF(.NOT.ASSOCIATED(GetMaterial(P1),GetMaterial(P2))) CYCLE

          n  = GetElementNOFNodes(Face)
          n1 = GetElementNOFNodes(P1)
          n2 = GetElementNOFNodes(P2)

          CALL LocalJumps( STIFF,Face,n,P1,n1,P2,n2)
          CALL DefaultUpdateEquations( STIFF, FORCE, Face )
       END IF
     END DO
!------------------------------------------------------------------------------
  END SUBROUTINE AddLocalFaceTerms
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE LocalJumps( STIFF,Face,n,P1,n1,P2,n2)
!------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=dp) :: STIFF(:,:)
      INTEGER :: n,n1,n2
      TYPE(Element_t), POINTER :: Face, P1, P2
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: FaceBasis(n), P1Basis(n1), P2Basis(n2)
      REAL(KIND=dp) :: Jump(n1+n2), detJ, U, V, W, S
      LOGICAL :: Stat
      INTEGER :: i, j, p, q, t, nFace, nParent
      TYPE(GaussIntegrationPoints_t) :: IntegStuff

      TYPE(Nodes_t) :: FaceNodes, P1Nodes, P2Nodes
      SAVE FaceNodes, P1Nodes, P2Nodes
!------------------------------------------------------------------------------
      STIFF = 0._dp

      CALL GetElementNodes(FaceNodes, Face)
      CALL GetElementNodes(P1Nodes, P1)
      CALL GetElementNodes(P2Nodes, P2)
!------------------------------------------------------------------------------
!     Numerical integration over the edge
!------------------------------------------------------------------------------
      IntegStuff = GaussPoints( Face )

      DO t=1,IntegStuff % n
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)
        S = IntegStuff % s(t)

        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo(Face, FaceNodes, U, V, W, detJ, FaceBasis)

        S = S * detJ

        ! Find basis functions for the parent elements:
        ! ---------------------------------------------
        CALL GetParentUVW(Face, n, P1, n1, U, V, W, FaceBasis)
        stat = ElementInfo(P1, P1Nodes, U, V, W, detJ, P1Basis)

        CALL GetParentUVW(Face, n, P2, n2, U, V, W, FaceBasis)
        stat = ElementInfo(P2, P2Nodes, U, V, W, detJ, P2Basis)

        ! Integrate jump terms:
        ! ---------------------
        Jump(1:n1) = P1Basis(1:n1)
        Jump(n1+1:n1+n2) = -P2Basis(1:n2)

        DO p=1,n1+n2
          DO q=1,n1+n2
            STIFF(p,q) = STIFF(p,q) + s * Jump(q)*Jump(p)
          END DO
        END DO
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE LocalJumps
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
    SUBROUTINE calcAverageFlux (Flux, Area, Element, n, nd, np, SOL, vDOFs)
!------------------------------------------------------------------------------
       IMPLICIT NONE
       INTEGER :: n, nd
       TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ,L(3),Normal(3)
       REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3), B(2, 3), &
                        SOL(2,32), Flux(2), Area
       LOGICAL :: Stat
       TYPE(GaussIntegrationPoints_t) :: IP
       INTEGER :: j, k, np, vDOFs

       TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
       CALL GetElementNodes( Nodes )
       IP = GaussPoints(Element)

       IF( dim == 2 ) THEN
         CALL Warn('CalcAverageFlux','Not implemented for 2D problems yet!')
       END IF

       B=0._dp

       DO j=1,IP % n
         stat = ElementInfo( Element, Nodes, IP % U(j), IP % V(j), &
             IP % W(j), detJ, Basis, dBasisdx, &
             EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver ) 
                  
         Normal = NormalVector( Element, Nodes, IP % U(j), IP % V(j), .TRUE. )         
         s = IP % s(j) * detJ

         DO k=1, vDOFs
           B(k,:) = MATMUL( SOL(k, np+1:nd), RotWBasis(1:nd-np,:) )
           Flux(k) = Flux(k) + s * SUM(Normal * B(k,:))
         END DO  

         Area = Area + s

      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE calcAverageFlux
!------------------------------------------------------------------------------

!------------------------------------------------------------------------
END SUBROUTINE MagnetoDynamicsCalcFields
!------------------------------------------------------------------------

