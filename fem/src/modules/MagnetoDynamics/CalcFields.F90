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
  LOGICAL :: Found, ElementalFields, RealField, LorentzConductivity, FoundVar
  INTEGER, POINTER :: Active(:)
  INTEGER :: mysolver,i,j,k,l,n,m,vDOFs, soln
  TYPE(ValueList_t), POINTER :: SolverParams, DGSolverParams
  TYPE(Solver_t), POINTER :: Solvers(:), PSolver

  LorentzConductivity = ListCheckPrefixAnyBodyForce(Model, "Angular Velocity") .OR. &
      ListCheckPrefixAnyBodyForce(Model, "Lorentz Velocity")

  ! This is really using DG so we don't need to make any dirty tricks to create DG fields
  ! as is done in this initialization. 
  SolverParams => GetSolverParams()

  ! The only purpose of this parsing of the variable name is to identify
  ! whether the field is real or complex. As the variable has not been
  ! created at this stage we have to do some dirty parsing. 
  pname = GetString(SolverParams, 'Potential variable', Found)
  vdofs = 0
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
        FoundVar = .TRUE.
        EXIT
      END IF
    END DO
  ELSE
    CALL ListAddString(SolverParams,'Potential Variable','av')
  END IF
    
  ! When we created the case for GUI where "av" is not given in sif then it is impossible to
  ! determine from the variable declaration what kind of solver we have. 
  IF( .NOT. FoundVar ) THEN
    DO i=1,Model % NumberOfSolvers
      sname = GetString(Model % Solvers(i) % Values, 'Procedure', Found)
      
      j = INDEX( sname,'WhitneyAVSolver')
      IF( j > 0 ) THEN
        CALL Info('MagnetoDynamicsCalcFields_Init0','The target solver seems to be real valued',Level=12)
        Vdofs = 1
        EXIT
      END IF

      j = INDEX( sname,'WhitneyAVHarmonicSolver')
      IF( j > 0 ) THEN
        CALL Info('MagnetoDynamicsCalcFields_Init0','The target solver seems to be complex valued',Level=12)
        Vdofs = 2
        EXIT
      END IF
    END DO

    IF( Vdofs == 0 ) THEN
      CALL Fatal('MagnetoDynamicsCalcFields_Init0','Could not determine target variable type (real or complex)')
    END IF
  END IF

  IF ( Vdofs==0 ) Vdofs=1

  RealField = ( Vdofs == 1 )
  CALL ListAddLogical( SolverParams, 'Target Variable Real Field', RealField ) 
  
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

  INTEGER  :: i
  LOGICAL :: Found, FluxFound, NodalFields, ElementalFields, &
      RealField, ComplexField, LorentzConductivity
  TYPE(ValueList_t), POINTER :: EQ, SolverParams

  LorentzConductivity = ListCheckPrefixAnyBodyForce(Model, "Angular Velocity") .or. &
    ListCheckPrefixAnyBodyForce(Model, "Lorentz Velocity")

  IF(.NOT.ASSOCIATED(Solver % Values)) Solver % Values=>ListAllocate()
  SolverParams => GetSolverParams()

  ! Inherit this from the _init0 solver. Hence we know it must exist!
  RealField = ListGetLogical( SolverParams,'Target Variable Real Field') 
  ComplexField = .NOT. RealField
  
  CALL ListAddString( SolverParams, 'Variable', '-nooutput hr_dummy' )
 
  CALL ListAddLogical( SolverParams, 'Linear System refactorize', .FALSE.)

  ! add these in the beginning, so that SaveData sees these existing, even
  ! if executed before the actual computations...
  ! -----------------------------------------------------------------------
  CALL ListAddConstReal(Model % Simulation,'res: Eddy current power',0._dp)
  CALL ListAddConstReal(Model % Simulation,'res: Magnetic Field Energy',0._dp)

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
    IF ( .NOT. ListCheckPresent(SolverParams,"Exported Variable "//TRIM(i2s(i))) ) EXIT
    i = i + 1
  END DO
  i = i - 1
  
  IF( NodalFields ) THEN
    i = i + 1
    IF ( RealField ) THEN
      CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
          "Magnetic Flux Density[Magnetic Flux Density:3]" )
    ELSE
      CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
          "Magnetic Flux Density[Magnetic Flux Density re:3 Magnetic Flux Density im:3]" )
    END IF

    IF (GetLogical(SolverParams,'Calculate Magnetic Vector Potential',Found)) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "Magnetic Vector Potential[Magnetic Vector Potential:3]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "Magnetic Vector Potential[Magnetic Vector Potential re:3 Magnetic Vector Potential im:3]")
      END IF
    END IF

    IF (GetLogical(SolverParams,'Calculate Magnetic Field Strength',Found)) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "Magnetic Field Strength[Magnetic Field Strength:3]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "Magnetic Field Strength[Magnetic Field Strength re:3 Magnetic Field Strength im:3]")
      END IF
    END IF

    IF (GetLogical(SolverParams,'Calculate JxB',Found)) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "JxB[JxB:3]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "JxB[JxB re:3 JxB im:3]")
      END IF
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Maxwell Stress', Found ) ) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "Maxwell Stress[Maxwell Stress:6]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "Maxwell Stress[Maxwell Stress re:6 Maxwell Stress im:6]" )
      END IF
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Current Density', Found ) ) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "Current Density[Current Density:3]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "Current Density[Current Density re:3 Current Density im:3]" )
      END IF
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Joule Heating', Found ) ) THEN
      i = i + 1
      CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
          "Joule Heating" )
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Harmonic Loss', Found ) ) THEN
      IF( RealField ) THEN
        CALL Warn('MagnetcDynamicsCalcFields',&
            'Harmonic loss computation only available for complex systems!')
      ELSE
        i = i + 1
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "Harmonic Loss Linear" )
        i = i + 1
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "Harmonic Loss Quadratic" )
      END IF
    END IF

    IF ( Transient .OR. .NOT. RealField .OR. LorentzConductivity) THEN
      IF ( GetLogical( SolverParams, 'Calculate Electric Field', Found ) ) THEN
        i = i + 1
        IF ( RealField ) THEN
          CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
              "Electric Field[Electric Field:3]" )
        ELSE
          CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
              "Electric Field[Electric Field re:3 Electric Field im:3]" )
        END IF
      END IF

      IF ( GetLogical( SolverParams, 'Calculate Winding Voltage', Found ) ) THEN
        i = i + 1
        IF ( RealField ) THEN
          CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
              "Winding Voltage" )
        ELSE
          CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
              "Winding Voltage[Winding Voltage re:1 Winding Voltage im:1]" )
        END IF
      END IF
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Nodal Heating', Found ) ) THEN
      i = i + 1
      CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
          "Nodal Joule Heating" )
    END IF

    IF (GetLogical(SolverParams, 'Calculate Nodal Forces', Found) ) THEN
      IF( RealField ) THEN
        i = i + 1
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "Nodal Force[Nodal Force:3]" )
      ELSE
        i = i + 1
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "Nodal Force[Nodal Force:3]" )
        CALL Warn('MagnetcDynamicsCalcFields',&
            'Calculating experimental average nodal forces. Use at own risk.')
      END IF
    END IF
  END IF
    
  ! If we have DG for the standard fields they are already elemental...
  IF (GetLogical(SolverParams,'Discontinuous Galerkin',Found)) RETURN

  ! Choose elemental if not otherwise specified. 
  ElementalFields = .NOT. GetLogical( SolverParams, 'Skip Elemental Fields', Found)
  IF(.NOT. Found ) ElementalFields = GetLogical( SolverParams, 'Calculate Elemental Fields', Found)
  IF(.NOT. Found ) ElementalFields = .TRUE.
  
  IF( ElementalFields ) THEN
    i = i + 1
    IF ( RealField ) THEN
      CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
          "-dg Magnetic Flux Density E[Magnetic Flux Density E:3]" )
    ELSE
      CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
          "-dg Magnetic Flux Density E[Magnetic Flux Density re E:3 Magnetic Flux Density im E:3]" )
    END IF

    IF (GetLogical(SolverParams,'Calculate Magnetic Vector Potential',Found)) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "-dg Magnetic Vector Potential E[Magnetic Vector Potential E:3]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "-dg Magnetic Vector Potential E[Magnetic Vector Potential re E:3 Magnetic Vector Potential im E:3]" )
      END IF
    END IF

    IF (GetLogical(SolverParams,'Calculate Magnetic Field Strength',Found)) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "-dg Magnetic Field Strength E[Magnetic Field Strength E:3]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "-dg Magnetic Field Strength E[Magnetic Field Strength re E:3 Magnetic Field Strength im E:3]" )
      END IF
    END IF

    IF (GetLogical(SolverParams,'Calculate JxB',Found)) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "-dg JxB E[JxB E:3]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "-dg JxB E[JxB re E:3 JxB im E:3]" )
      END IF
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Maxwell Stress', Found ) ) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "-dg Maxwell Stress E[Maxwell Stress E:6]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "-dg Maxwell Stress E[Maxwell Stress re E:6 Maxwell Stress im E:6]" )
      END IF
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Current Density', Found ) ) THEN
      i = i + 1
      IF ( RealField ) THEN
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "-dg Current Density E[Current Density E:3]" )
      ELSE
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "-dg Current Density E[Current Density re E:3 Current Density im E:3]" )
      END IF
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Joule Heating', Found ) ) THEN
      i = i + 1
      CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
          "-dg Joule Heating E" )
    END IF

    IF ( GetLogical( SolverParams, 'Calculate Harmonic Loss', Found ) ) THEN
      IF( RealField ) THEN
        CALL Warn('MagnetoDynamicsCalcFields',&
            'Harmonic loss computation only available for complex systems!')
      ELSE
        i = i + 1
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "-dg Harmonic Loss Linear E" )
        i = i + 1
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "-dg Harmonic Loss Quadratic E" )
      END IF
    END IF

    IF ( Transient .OR. ComplexField .OR. LorentzConductivity ) THEN
      IF ( GetLogical( SolverParams, 'Calculate Electric Field', Found ) ) THEN
        i = i + 1
        IF ( RealField ) THEN
          CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
              "-dg Electric Field E[Electric Field E:3]" )
        ELSE
          CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
              "-dg Electric Field E[Electric Field re E:3 Electric Field im E:3]" )
        END IF
      END IF

      IF ( GetLogical( SolverParams, 'Calculate Winding Voltage', Found ) ) THEN
        i = i + 1
        IF ( RealField ) THEN
          CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
              "-dg Winding Voltage E" )
        ELSE
          CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
              "-dg Winding Voltage E[Winding Voltage re E:1 Winding Voltage im E:1]" )
        END IF
      END IF
    END IF

    IF (GetLogical(SolverParams, 'Calculate Nodal Forces', Found) ) THEN
      IF( RealField ) THEN
        i = i + 1
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "-dg Nodal Force E[Nodal Force E:3]" )
      ELSE
        i = i + 1
        CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
            "-dg Nodal Force E[Nodal Force E:3]" )
        CALL Warn('MagnetcDynamicsCalcFields',&
            'Calculating experimental average nodal forces. Use at own risk.')
      END IF
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
   TYPE(Solver_t) :: Solver
   TYPE(Model_t) :: Model
   REAL(KIND=dp) :: dt
   LOGICAL :: Transient
!------------------------------------------------------------------------------
!  The following arrays have hard-coded sizes which may need to be altered if
!  new finite elements are added. Current defaults are 54 edge finite element 
!  DOFs and 27 nodal DOFs at maximum (obtained for the second-order brick over
!  a background element of type 827):
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: WBasis(54,3), RotWBasis(54,3), Basis(27), dBasisdx(27,3)
   REAL(KIND=dp) :: SOL(2,81), PSOL(81), ElPotSol(1,27), R(27), C(27)
   REAL(KIND=dp) :: Wbase(27), alpha(27), NF_ip(27,3)
   REAL(KIND=dp) :: PR(27), omega_velo(3,27), lorentz_velo(3,27)
   COMPLEX(KIND=dp) :: Magnetization(3,27), BodyForceCurrDens(3,27) 
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: s,u,v,w, Norm
   REAL(KIND=dp) :: B(2,3), E(2,3), JatIP(2,3), VP_ip(2,3), JXBatIP(2,3), CC_J(2,3), B2
   REAL(KIND=dp) :: detJ, C_ip, R_ip, PR_ip, ST(3,3), Omega, Power, Energy, w_dens, R_t_ip(3,3)
   REAL(KIND=dp) :: Freq, FreqPower, FieldPower, LossCoeff, ValAtIP
   REAL(KIND=dp) :: Freq2, FreqPower2, FieldPower2, LossCoeff2
   REAL(KIND=dp) :: ComponentLoss(2,2), rot_velo(3) 
   REAL(KIND=dp) :: Coeff, Coeff2, TotalLoss(3), LumpedForce(3), localAlpha, localV(2), nofturns, coilthickness
   REAL(KIND=dp) :: Flux(2), AverageFluxDensity(2), Area, N_j, wvec(3), PosCoord(3), TorqueDeprecated(3)

   COMPLEX(KIND=dp) :: MG_ip(3), BodyForceCurrDens_ip(3)
   COMPLEX(KIND=dp) :: CST(3,3)
   COMPLEX(KIND=dp) :: CMat_ip(3,3)  
   COMPLEX(KIND=dp) :: imag_value

   INTEGER, PARAMETER :: ind1(6) = [1,2,3,1,2,1]
   INTEGER, PARAMETER :: ind2(6) = [1,2,3,2,3,3]

   TYPE(Variable_t), POINTER :: Var, MFD, MFS, CD, EF, MST, &
                                JH, NJH, VP, FWP, JXB, ML, ML2, LagrangeVar, NF
   TYPE(Variable_t), POINTER :: EL_MFD, EL_MFS, EL_CD, EL_EF, &
                                EL_MST, EL_JH, EL_VP, EL_FWP, EL_JXB, EL_ML, EL_ML2, &
                                EL_NF

   INTEGER :: Active,i,j,k,l,m,n,nd,np,p,q,DOFs,vDOFs,dim,BodyId,&
              VvarDofs,VvarId,IvarId,Reindex,Imindex,EdgeBasisDegree

   TYPE(Solver_t), POINTER :: pSolver, ElPotSolver
   CHARACTER(LEN=MAX_NAME_LEN) :: Pname, CoilType, ElectricPotName, LossFile, CurrPathPotName

   TYPE(ValueList_t), POINTER :: Material, BC, BodyForce, BodyParams, SolverParams
   LOGICAL :: Found, FoundMagnetization, stat, Cubic, LossEstimation, &
              CalcFluxLogical, CoilBody, PreComputedElectricPot, ImposeCircuitCurrent, &
              ItoJCoeffFound, ImposeBodyForceCurrent, HasVelocity, HasAngularVelocity, &
              HasLorenzVelocity, HaveAirGap, UseElementalNF, HasTensorReluctivity, &
              ImposeBodyForcePotential, JouleHeatingFromCurrent, HasZirka
   
   TYPE(GaussIntegrationPoints_t) :: IP
   TYPE(Nodes_t), SAVE :: Nodes
   TYPE(Element_t), POINTER :: Element

   INTEGER, ALLOCATABLE :: Pivot(:), TorqueGroups(:)
   INTEGER, POINTER :: MasterBodies(:)

   REAL(KIND=dp), POINTER CONTIG :: Fsave(:)
   REAL(KIND=dp), POINTER :: HB(:,:)=>NULL(), CubicCoeff(:)=>NULL(), &
     HBBVal(:), HBCval(:), HBHval(:)
   REAL(KIND=dp) :: Babs
   TYPE(Mesh_t), POINTER :: Mesh
   REAL(KIND=dp), ALLOCATABLE, TARGET :: Gforce(:,:), MASS(:,:), FORCE(:,:)
   REAL(KIND=dp), ALLOCATABLE :: BodyLoss(:,:), RotM(:,:,:), Torque(:)

   REAL(KIND=DP), POINTER :: Cwrk(:,:,:)=>NULL(), Cwrk_im(:,:,:)=>NULL()
   COMPLEX(KIND=dp), ALLOCATABLE :: Tcoef(:,:,:)
   REAL(KIND=dp), POINTER :: R_t(:,:,:)

   LOGICAL :: PiolaVersion, ElementalFields, NodalFields, RealField, SecondOrder
   REAL(KIND=dp) :: ItoJCoeff, CircuitCurrent
   TYPE(ValueList_t), POINTER :: CompParams
   REAL(KIND=dp) :: DetF, F(3,3), G(3,3), GT(3,3)
   REAL(KIND=dp), ALLOCATABLE :: EBasis(:,:), CurlEBasis(:,:) 
   LOGICAL :: CSymmetry, HBCurve, LorentzConductivity
   REAL(KIND=dp) :: xcoord, grads_coeff, val
   TYPE(ValueListEntry_t), POINTER :: HBLst
   REAL(KIND=dp) :: HarmPowerCoeff = 0.5_dp
   
   INTEGER, POINTER, SAVE :: SetPerm(:) => NULL()
!-------------------------------------------------------------------------------------------
   IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
   
   CALL Info('MagnetoDynamicsCalcFields','------------------------------',Level=6)
   CALL Info('MagnetoDynamicsCalcFields','Computing postprocessed fields',Level=5)
   
   dim = CoordinateSystemDimension()
   SolverParams => GetSolverParams()

   IF (GetLogical(SolverParams, 'Calculate harmonic peak power', Found)) HarmPowerCoeff = 1.0_dp

   Pname = GetString(SolverParams, 'Potential Variable',Found)
   IF(.NOT. Found ) Pname = 'av'
   Found = .FALSE.
   DO i=1,Model % NumberOfSolvers
     pSolver => Model % Solvers(i)
     IF ( Pname == getVarName(pSolver % Variable)) THEN
       Found = .TRUE.
       EXIT
     END IF
   END DO

   IF(.NOT. Found ) THEN
     CALL Fatal('MagnetoDynamicsCalcFields','Solver associated to potential variable > '&
         //TRIM(Pname)//' < not found!')
   END IF

   ! Inherit the solution basis from the primary solver
   vDOFs = pSolver % Variable % DOFs
   SecondOrder = GetLogical( pSolver % Values, 'Quadratic Approximation', Found )  
   IF (SecondOrder) THEN
     EdgeBasisDegree = 2
   ELSE
     EdgeBasisDegree = 1
   END IF

   IF( SecondOrder ) THEN
     PiolaVersion = .TRUE.
   ELSE
     PiolaVersion = GetLogical( pSolver % Values,'Use Piola Transform', Found ) 
   END IF

   IF (PiolaVersion) &
       CALL Info('MagnetoDynamicsCalcFields', &
       'Using Piola transformed finite elements',Level=5)

   ElectricPotName = GetString(SolverParams, 'Precomputed Electric Potential', PrecomputedElectricPot)
   IF (PrecomputedElectricPot) THEN
     DO i=1, Model % NumberOfSolvers
       ElPotSolver => Model % Solvers(i)
       IF (ElectricPotName==getVarName(ElPotSolver % Variable)) EXIT
     END DO
   END IF
     
   ! Do we have a real or complex valued primary field?
   RealField = ( vDofs == 1 ) 

   LorentzConductivity = ListCheckPrefixAnyBodyForce(Model, "Angular Velocity") .or. &
       ListCheckPrefixAnyBodyForce(Model, "Lorentz Velocity")

   Mesh => GetMesh()
   LagrangeVar => VariableGet( Solver % Mesh % Variables,'LagrangeMultiplier', ThisOnly=.TRUE.)

   MFD => VariableGet( Mesh % Variables, 'Magnetic Flux Density' )
   EL_MFD => VariableGet( Mesh % Variables, 'Magnetic Flux Density E' )

   MFS => VariableGet( Mesh % Variables, 'Magnetic Field Strength')
   EL_MFS => VariableGet( Mesh % Variables, 'Magnetic Field Strength E')

   VP => VariableGet( Mesh % Variables, 'Magnetic Vector Potential')
   EL_VP => VariableGet( Mesh % Variables, 'Magnetic Vector Potential E')
   
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
   JXB => NULL(); EL_JXB => NULL();
   ML  => NULL(); EL_ML => NULL();
   ML2 => NULL(); EL_ML2 => NULL();
   NF => NULL(); EL_NF => NULL();
   NJH => NULL()
   
   IF ( Transient .OR. .NOT. RealField .OR. LorentzConductivity ) THEN
     EF => VariableGet( Mesh % Variables, 'Electric Field' )
     FWP => VariableGet( Mesh % Variables, 'Winding Voltage' )

     EL_EF => VariableGet( Mesh % Variables, 'Electric Field E' )
     EL_FWP => VariableGet( Mesh % Variables, 'Winding Voltage E' )
   END IF

   !IF( RealField ) THEN
     NF => VariableGet( Mesh % Variables, 'Nodal Force') 
     EL_NF => VariableGet( Mesh % Variables, 'Nodal Force E')
   !END IF

   CD => VariableGet( Mesh % Variables, 'Current Density' )
   EL_CD => VariableGet( Mesh % Variables, 'Current Density E' )

   JH => VariableGet( Mesh % Variables, 'Joule Heating' )
   EL_JH => VariableGet( Mesh % Variables, 'Joule Heating E' )

   NJH => VariableGet( Mesh % Variables, 'Nodal Joule Heating' )
   
   IF(.NOT. RealField ) THEN
     ML => VariableGet( Mesh % Variables, 'Harmonic Loss Linear')
     EL_ML => VariableGet( Mesh % Variables, 'Harmonic Loss Linear E')
     ML2 => VariableGet( Mesh % Variables, 'Harmonic Loss Quadratic')
     EL_ML2 => VariableGet( Mesh % Variables, 'Harmonic Loss Quadratic E')
   END IF

   JXB => VariableGet( Mesh % Variables, 'JxB')
   EL_JXB => VariableGet( Mesh % Variables, 'JxB E')

   MST => variableGet( Mesh % Variables, 'Maxwell stress' )
   EL_MST => variableGet( Mesh % Variables, 'Maxwell stress E' )

   DOFs = 0 
   IF ( ASSOCIATED(MFD) ) DOFs=DOFs+3
   IF ( ASSOCIATED(MFS) ) DOFs=DOFs+3
   IF ( ASSOCIATED(VP)  ) DOFs=DOFs+3
   IF ( ASSOCIATED(CD)  ) DOFs=DOFs+3
   IF ( ASSOCIATED(FWP) ) DOFs=DOFs+1
   IF ( ASSOCIATED(EF)  ) DOFs=DOFs+3
   IF ( ASSOCIATED(JXB) ) DOFs=DOFs+3
   IF ( ASSOCIATED(MST) ) DOFs=DOFs+6
   IF ( ASSOCIATED(NF)  ) DOFs=DOFs+3
   DOFs = DOFs*vDOFs
   IF ( ASSOCIATED(JH) .OR. ASSOCIATED(NJH)) DOFs=DOFs+1
   IF ( ASSOCIATED(ML) ) DOFs=DOFs+1
   IF ( ASSOCIATED(ML2) ) DOFs=DOFs+1
   NodalFields = DOFs > 0

   IF(NodalFields) THEN
     ALLOCATE(GForce(SIZE(Solver % Matrix % RHS),DOFs)); Gforce=0._dp
   ELSE
     DOFs = 0 
     IF ( ASSOCIATED(EL_MFD) ) DOFs=DOFs+3
     IF ( ASSOCIATED(EL_MFS) ) DOFs=DOFs+3
     IF ( ASSOCIATED(EL_VP)  ) DOFs=DOFs+3
     IF ( ASSOCIATED(EL_CD)  ) DOFs=DOFs+3
     IF ( ASSOCIATED(EL_FWP) ) DOFs=DOFs+1
     IF ( ASSOCIATED(EL_EF)  ) DOFs=DOFs+3
     IF ( ASSOCIATED(EL_JXB) ) DOFs=DOFs+3
     IF ( ASSOCIATED(EL_MST) ) DOFs=DOFs+6
     IF ( ASSOCIATED(EL_NF) ) DOFs=DOFs+3 
     DOFs = DOFs*vDOFs
     IF ( ASSOCIATED(EL_NF) ) DOFs=DOFs+3 
     IF ( ASSOCIATED(EL_JH) ) DOFs=DOFs+1
     IF ( ASSOCIATED(EL_ML) ) DOFs=DOFs+1
     IF ( ASSOCIATED(EL_ML2) ) DOFs=DOFs+1
   END IF

   ElementalFields = .FALSE.
   IF ( ASSOCIATED(EL_MFD) ) ElementalFields=.TRUE.
   IF ( ASSOCIATED(EL_MFS) ) ElementalFields=.TRUE.
   IF ( ASSOCIATED(EL_VP)  ) ElementalFields=.TRUE.
   IF ( ASSOCIATED(EL_CD)  ) ElementalFields=.TRUE.
   IF ( ASSOCIATED(EL_FWP) ) ElementalFields=.TRUE.
   IF ( ASSOCIATED(EL_EF)  ) ElementalFields=.TRUE.
   IF ( ASSOCIATED(EL_JXB) ) ElementalFields=.TRUE.
   IF ( ASSOCIATED(EL_MST) ) ElementalFields=.TRUE.
   IF ( ASSOCIATED(EL_NF)  ) ElementalFields=.TRUE.
   IF ( ASSOCIATED(EL_JH)  ) ElementalFields=.TRUE.
   IF ( ASSOCIATED(EL_ML)  ) ElementalFields=.TRUE.
   IF ( ASSOCIATED(EL_ML2)  ) ElementalFields=.TRUE.

   n = Mesh % MaxElementDOFs
   ALLOCATE( MASS(n,n), FORCE(n,DOFs), Tcoef(3,3,n), RotM(3,3,n), Pivot(n), R_t(3,3,n))

   SOL = 0._dp; PSOL=0._dp

   LossEstimation = GetLogical(SolverParams,'Loss Estimation',Found) &
       .OR. ASSOCIATED( ML ) .OR. ASSOCIATED( EL_ML ) &
       .OR. ASSOCIATED( ML2 ) .OR. ASSOCIATED( EL_ML2 ) 

   IF (LossEstimation) THEN
      FreqPower = GetCReal( SolverParams,'Harmonic Loss Linear Frequency Exponent',Found )
      IF( .NOT. Found ) FreqPower = 1.0_dp

      FreqPower2 = GetCReal( SolverParams,'Harmonic Loss Quadratic Frequency Exponent',Found )
      IF( .NOT. Found ) FreqPower2 = 2.0_dp

      FieldPower = GetCReal( SolverParams,'Harmonic Loss Linear Exponent',Found ) 
      IF( .NOT. Found ) FieldPower = 2.0_dp
      FieldPower = FieldPower / 2.0_dp

      FieldPower2 = GetCReal( SolverParams,'Harmonic Loss Quadratic Exponent',Found ) 
      IF( .NOT. Found ) FieldPower2 = 2.0_dp
      FieldPower2 = FieldPower2 / 2.0_dp

      IF(.NOT. ListCheckPresentAnyMaterial( Model,'Harmonic Loss Linear Coefficient') ) THEN
        CALL Warn('MagnetoDynamicsCalcFields',&
            'Harmonic loss requires > Harmonic Loss Linear Coefficient < in material section!')
      END IF

      IF(.NOT. ListCheckPresentAnyMaterial( Model,'Harmonic Loss Quadratic Coefficient') ) THEN
        CALL Warn('MagnetoDynamicsCalcFields',&
            'Harmonic loss requires > Harmonic Loss Quadratic Coefficient < in material section!')
      END IF

      ComponentLoss = 0.0_dp
      ALLOCATE( BodyLoss(3,Model % NumberOfBodies) )
      BodyLoss = 0.0_dp
   END IF


   C = 0._dp; R=0._dp; PR=0._dp
   Magnetization = 0._dp

   Power = 0._dp; Energy = 0._dp
   CALL DefaultInitialize()
   
   
   DO i = 1, GetNOFActive()
     Element => GetActiveElement(i)
     n = GetElementNOFNodes()
     np = n*pSolver % Def_Dofs(GetElementFamily(Element),Element % BodyId,1)
     nd = GetElementNOFDOFs(uSolver=pSolver)

     IF (SIZE(Tcoef,3) /= n) THEN
       DEALLOCATE(Tcoef)
       ALLOCATE(Tcoef(3,3,n))
     END IF
     
     CALL GetElementNodes( Nodes )

     ! If potential is not available we have to use given current directly to estimate Joule losses
     JouleHeatingFromCurrent = ( np == 0 .AND. &
         .NOT. ( PreComputedElectricPot .OR. ImposeBodyForcePotential ) )
     
     BodyId = GetBody()
     Material => GetMaterial()
     BodyForce => GetBodyForce()

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
       CALL GetScalarLocalSolution(PSOL,Pname,uSolver=pSolver,Tstep=-1)
       PSOL(1:nd)=(SOL(1,1:nd)-PSOL(1:nd))/dt
     END IF

     Omega = GetAngularFrequency(pSOlver % Values,Found,Element)
     IF( .NOT. ( RealField .OR. Found ) ) THEN
       CALL Fatal('MagnetoDynamicsCalcFields',&
           '(Angular) Frequency must be given for complex fields!')
     END IF
     Freq = Omega / (2*PI)
     
     IF ( ASSOCIATED(MFS) ) THEN
       FoundMagnetization = .FALSE.
       IF(ASSOCIATED(BodyForce)) THEN
         CALL GetComplexVector( BodyForce,Magnetization(1:3,1:n),'Magnetization',FoundMagnetization)
       END IF

       IF(.NOT.FoundMagnetization) THEN
         CALL GetComplexVector( BodyForce,Magnetization(1:3,1:n),'Magnetization',FoundMagnetization)
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
       
       !CALL GetLocalSolution(Wbase, 'w')
       Call GetWPotential(Wbase)
  
       SELECT CASE (CoilType)
       CASE ('stranded')
         IvarId = GetInteger (CompParams, 'Circuit Current Variable Id', Found)
         IF (.NOT. Found) CALL Fatal ('MagnetoDynamicsCalcFields', 'Circuit Current Variable Id not found!')

         N_j = GetConstReal (CompParams, 'Stranded Coil N_j', Found)
         IF (.NOT. Found) CALL Fatal ('MagnetoDynamicsCalcFields', 'Stranded Coil N_j not found!')

         nofturns = GetConstReal(CompParams, 'Number of Turns', Found)
         IF (.NOT. Found) CALL Fatal('MagnetoDynamicsCalcFields','Stranded Coil: Number of Turns not found!')
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

         nofturns = GetConstReal(CompParams, 'Number of Turns', Found)
         IF (.NOT. Found) CALL Fatal('MagnetoDynamicsCalcFields','Foil Winding: Number of Turns not found!')

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

     HasTensorReluctivity = .FALSE.
     CALL GetConstRealArray( Material, HB, 'H-B curve', Found )
     IF ( ASSOCIATED(HB) ) THEN
      Cubic = GetLogical( Material, 'Cubic spline for H-B curve', Found)
      l = SIZE(HB,1)
      HBBval => HB(:,1)
      HBHval => HB(:,2)
      IF(l>1) THEN
        IF (Cubic.AND..NOT.ASSOCIATED(CubicCoeff) ) THEN
          ALLOCATE(CubicCoeff(l))
          CALL CubicSpline(l,HB(:,1),HB(:,2),CubicCoeff)
        END IF
        HBCval => CubicCoeff
      END IF

      IF(l<=1) THEN
        HBLst => ListFind(Material,'H-B Curve',HBcurve)
        IF(HBcurve) THEN
          HBCval => HBLst % CubicCoeff
          HBBval => HBLst % TValues
          HBHval => HBLst % FValues(1,1,:)
        END IF
      END IF
     ELSE
       CALL GetReluctivity(Material,R_t,n,HasTensorReluctivity)
       IF(.NOT. HasTensorReluctivity) CALL GetReluctivity(Material,R,n)
     END IF

     HasVelocity = .FALSE.
     IF(ASSOCIATED(BodyForce)) THEN
       CALL GetRealVector( BodyForce, omega_velo, 'Angular velocity', HasAngularVelocity)
       CALL GetRealVector( BodyForce, lorentz_velo, 'Lorentz velocity', HasLorenzVelocity)
       HasVelocity = HasAngularVelocity .OR. HasLorenzVelocity
     END IF
     

     ! Calculate nodal fields:
     ! -----------------------
     IF (SecondOrder) THEN
        IP = GaussPoints(Element, EdgeBasis=dim==3, PReferenceElement=PiolaVersion, EdgeBasisDegree=EdgeBasisDegree)
     ELSE
        IP = GaussPoints(Element, EdgeBasis=dim==3, PReferenceElement=PiolaVersion)
     END IF

     MASS  = 0._dp
     FORCE = 0._dp
     E = 0._dp; B=0._dp

     haszirka = .false.
     if(ASSOCIATED(MFS) .OR. ASSOCIATED(el_MFS)) THEN
       CALL GetHystereticMFS(Element, force(:,4:6), pSolver, HasZirka, CSymmetry=CSymmetry)
     end if

     DO j = 1,IP % n
       u = IP % U(j)
       v = IP % V(j)
       w = IP % W(j)

       IF (PiolaVersion) THEN
          stat = EdgeElementInfo( Element, Nodes, u, v, w, DetF=DetJ, Basis=Basis, &
               EdgeBasis=WBasis, RotBasis=RotWBasis, dBasisdx=dBasisdx, &
               BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)
       ELSE
          stat=ElementInfo(Element,Nodes,u,v,w,detJ,Basis,dBasisdx)
          IF( dim == 3 ) THEN
            CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)
          END IF
       END IF


       grads_coeff = -1._dp/GetCircuitModelDepth()
       IF( CSymmetry ) THEN
         xcoord = SUM( Basis(1:n) * Nodes % x(1:n) )
         grads_coeff = grads_coeff/xcoord
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
            B(k,:) = MATMUL( SOL(k,np+1:nd), RotWBasis(1:nd-np,:) )
         END SELECT
       END DO
       IF(ImposeCircuitCurrent .and. ItoJCoeffFound) THEN
         wvec = -MATMUL(Wbase(1:n), dBasisdx(1:n,:))
         IF(SUM(wvec**2._dp) .GE. AEPS) THEN
           wvec = wvec/SQRT(SUM(wvec**2._dp))
         ELSE
           wvec = [0.0_dp, 0.0_dp, 1.0_dp]
         END IF
       END IF

       ! Compute convection type term coming from rotation
       ! -------------------------------------------------
       IF(HasVelocity) THEN
         rot_velo = 0.0_dp
         IF( HasAngularVelocity ) THEN
           DO k=1,n
             rot_velo(1:3) = rot_velo(1:3) + CrossProduct(omega_velo(1:3,k), [ &
                 basis(k) * Nodes % x(k), &
                 basis(k) * Nodes % y(k), &
                 basis(k) * Nodes % z(k)])
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
       
       IF (vDOFs > 1) THEN   ! Complex case
         IF (CoilType /= 'stranded') THEN
           ! -j * Omega A
           SELECT CASE(dim)
           CASE(2)
             E(1,:) = 0._dp
             E(2,:) = 0._dp
             E(1,3) =  Omega*SUM(SOL(2,1:nd) * Basis(1:nd))
             E(2,3) = -Omega*SUM(SOL(1,1:nd) * Basis(1:nd))
           CASE(3)
             E(1,:) = Omega*MATMUL(SOL(2,np+1:nd),WBasis(1:nd-np,:))
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
             wvec = -MATMUL(Wbase(1:np), dBasisdx(1:np,:))
             wvec = wvec/SQRT(SUM(wvec**2._dp))
           END SELECT
           IF(CMat_ip(3,3) /= 0._dp ) THEN
             imag_value = LagrangeVar % Values(IvarId) + im * LagrangeVar % Values(IvarId+1)
             E(1,:) = E(1,:)+REAL(imag_value * N_j * wvec / CMat_ip(3,3))
             E(2,:) = E(2,:)+AIMAG(imag_value * N_j * wvec / CMat_ip(3,3))
           END IF

         CASE ('massive')
           localV(1) = localV(1) + LagrangeVar % Values(VvarId)
           localV(2) = localV(2) + LagrangeVar % Values(VvarId+1)
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
             localV(1) = localV(1) + LagrangeVar % Values(VvarId+Reindex) * localAlpha**(k-1)
             localV(2) = localV(2) + LagrangeVar % Values(VvarId+Imindex) * localAlpha**(k-1)
           END DO
           SELECT CASE(dim)
           CASE(2)
             E(1,3) = E(1,3)-localV(1) * grads_coeff
             E(2,3) = E(2,3)-localV(2) * grads_coeff
           CASE(3)
             E(1,:) = E(1,:)-localV(1) * MATMUL(Wbase(1:np), dBasisdx(1:np,:))
             E(2,:) = E(2,:)-localV(2) * MATMUL(Wbase(1:np), dBasisdx(1:np,:))
           END SELECT

         CASE DEFAULT
           ! -Grad(V)
           IF(dim==3) THEN
             E(1,:) = E(1,:)-MATMUL(SOL(1,1:np), dBasisdx(1:np,:))
             E(2,:) = E(2,:)-MATMUL(SOL(2,1:np), dBasisdx(1:np,:))
           END IF

           IF( ImposeBodyForcePotential ) THEN
             E(1,:) = E(1,:) - MATMUL(ElPotSol(1,1:n), dBasisdx(1:n,:))
           END IF             
         END SELECT
         
       ELSE   ! Real case
         IF (CoilType /= 'stranded') THEN 
           SELECT CASE(dim)
           CASE(2)
             E(1,1) = 0._dp
             E(1,2) = 0._dp
             E(1,3) = -SUM(PSOL(1:nd) * Basis(1:nd))
           CASE(3)
             E(1,:) = -MATMUL(PSOL(np+1:nd), Wbasis(1:nd-np,:))
           END SELECT
         ELSE
           E(1,:) = 0._dp
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
             wvec = -MATMUL(Wbase(1:np), dBasisdx(1:np,:))
             wvec = wvec/SQRT(SUM(wvec**2._dp))
             IF(CMat_ip(3,3) /= 0._dp ) &
               E(1,:) = E(1,:)+ LagrangeVar % Values(IvarId) * N_j * wvec / CMat_ip(3,3)
           END SELECT

         CASE ('massive')
           localV(1) = localV(1) + LagrangeVar % Values(VvarId)
           SELECT CASE(dim)
           CASE(2)
             E(1,3) = E(1,3)-localV(1) * grads_coeff
           CASE(3)
             E(1,:) = E(1,:)-localV(1) * MATMUL(Wbase(1:np), dBasisdx(1:np,:))
           END SELECT

         CASE ('foil winding')
           localAlpha = coilthickness *SUM(alpha(1:np) * Basis(1:np)) 
           DO k = 1, VvarDofs-1
             localV(1) = localV(1) + LagrangeVar % Values(VvarId+k) * localAlpha**(k-1)
           END DO
           SELECT CASE(dim)
           CASE(2)
             E(1,3) = E(1,3)-localV(1) * grads_coeff
           CASE(3)
             E(1,:) = E(1,:)-localV(1) * MATMUL(Wbase(1:np), dBasisdx(1:np,:))
           END SELECT

         CASE DEFAULT
           IF(dim==3 .AND. Transient) THEN
             E(1,:) = E(1,:)-MATMUL(SOL(1,1:np), dBasisdx(1:np,:))
           END IF

           IF (np > 0 .AND. dim==3 .AND. .NOT. Transient) THEN
             E(1,:) = -MATMUL(SOL(1,1:np), dBasisdx(1:np,:))
           ELSE IF ( PrecomputedElectricPot ) THEN
             E(1,:) = -MATMUL(ElPotSol(1,1:n), dBasisdx(1:n,:))
           END IF

           IF( ImposeBodyForcePotential ) THEN
             E(1,:) = E(1,:) - MATMUL(ElPotSol(1,1:n), dBasisdx(1:n,:))
           END IF

         END SELECT
       END IF
       

       IF ( ASSOCIATED(HB) ) THEN
         Babs=SQRT(SUM(B(1,:)**2))
         R_ip = InterpolateCurve(HBBval,HBHval,Babs,HBCval)/Babs
         w_dens = IntegrateCurve(HBBval,HBHval,HBCval,0._dp,Babs)
       ELSE
         R_ip = SUM( Basis(1:n)*R(1:n) )
         IF(HasTensorReluctivity) THEN
           DO k = 1,3
             DO l = 1,3
               R_t_ip(k,l) = sum(Basis(1:n)*R_t(k,l,1:n))
             END DO
           END DO
           w_dens = 0.5*SUM(B(1,:)*MATMUL(R_t_ip,B(1,:)))
         END IF
         w_dens = 0.5*R_ip*SUM(B(1,:)**2)
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
             VP_ip(l,:)=MATMUL(SOL(l,np+1:nd),WBasis(1:nd-np,:))
           END SELECT
         END DO
       END IF
       
       IF (ASSOCIATED(NF).OR.ASSOCIATED(EL_NF)) THEN
         NF_ip = 0._dp
         B2 = sum(B(1,:)*B(1,:) + B(2,:)*B(2,:))
         DO k=1,n
           DO l=1,3
             DO m=1,3
               NF_ip(k,l) = NF_ip(k,l) - (R_ip*(B(1,l)*B(1,m)))*dBasisdx(k,m)
             END DO
             NF_ip(k,l) = NF_ip(k,l) + (R_ip*B2-w_dens)*dBasisdx(k,l)
           END DO
         END DO

         IF (.NOT. RealField) THEN
           DO k=1,n
             DO l=1,3
               DO m=1,3
                 NF_ip(k,l) = NF_ip(k,l) - (R_ip*(B(2,l)*B(2,m)))*dBasisdx(k,m)
               END DO
             END DO
           END DO
         END IF
       END IF

       s = IP % s(j) * detJ

       IF(ASSOCIATED(HB) .AND. RealField) THEN 
         Energy = Energy + s*(0.5*PR_ip*SUM(E**2) + w_dens)
       ELSE
         Energy = Energy + s*0.5*(PR_ip*SUM(E**2) + R_ip*SUM(B**2))
       END IF

       DO p=1,n
         DO q=1,n
           MASS(p,q)=MASS(p,q)+s*Basis(p)*Basis(q)
         END DO
         k = 0
         DO l=1,vDOFs
           FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+s*B(l,:)*Basis(p)
           k = k+3
         END DO

         IF ( (ASSOCIATED(MFS).OR.ASSOCIATED(EL_MFS)) .and. .not. HasZirka) THEN
           IF(.NOT. HasZirka) then
             FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+s*(R_ip*B(1,:)-REAL(MG_ip))*Basis(p)
             k = k+3
             IF ( Vdofs>1 ) THEN
               FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+s*(R_ip*B(2,:)-AIMAG(MG_ip))*Basis(p)
               k = k+3
             END IF
           ELSE
             FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)-s*(REAL(MG_ip))*Basis(p)
           END IF
         END IF
         IF ( ASSOCIATED(VP).OR.ASSOCIATED(EL_VP)) THEN
           DO l=1,vDOFs
             FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+s*VP_ip(l,:)*Basis(p)
             k = k+3
           END DO
         END IF
         IF ( ASSOCIATED(EF).OR.ASSOCIATED(EL_EF)) THEN
           DO l=1,vDOFs
             FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3)+s*E(l,:)*Basis(p)
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
                IF( HasVelocity ) THEN
                  JatIP(1,l) = JatIP(1,l) + SUM( REAL(CMat_ip(l,1:3)) * CrossProduct(rot_velo, B(1,1:3)))
                END IF
                FORCE(p,k+l) = FORCE(p,k+l)+s*JatIp(1,l)*Basis(p)
              END DO
              k = k+3
           ELSE
              DO l=1,3
                JatIp(1,l) = SUM( REAL(CMat_ip(l,1:3)) * E(1,1:3) ) - &
                             SUM( AIMAG(CMat_ip(l,1:3)) * E(2,1:3) ) + REAL(BodyForceCurrDens_ip(l))
                IF( HasVelocity ) THEN
                  JatIp(1,l) = JatIp(1,l) + SUM( REAL(CMat_ip(l,1:3)) * CrossProduct(rot_velo, B(1,1:3)) ) - &
                               SUM( AIMAG(CMat_ip(l,1:3)) * CrossProduct(rot_velo, B(2,1:3)) )
                END IF
                FORCE(p,k+l) = FORCE(p,k+l)+s*JatIp(1,l)*Basis(p)
              END DO
              k = k+3
              DO l=1,3
                JatIp(2,l) = SUM( AIMAG(CMat_ip(l,1:3)) * E(1,1:3) ) + &
                             SUM( REAL(CMat_ip(l,1:3)) * E(2,1:3) ) + AIMAG(BodyForceCurrDens_ip(l))
                IF( HasVelocity ) THEN
                  JatIp(2,l) = JatIp(2,l) + SUM( AIMAG(CMat_ip(l,1:3)) * CrossProduct(rot_velo, B(1,1:3)) ) + &
                               SUM( REAL(CMat_ip(l,1:3)) * CrossProduct(rot_velo, B(2,1:3)) )
                END IF
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

         IF (vDOFS == 1) THEN
           IF( JouleHeatingFromCurrent ) THEN
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
           IF (HasVelocity) THEN
             Coeff = Coeff + SUM(MATMUL(REAL(CMat_ip), CrossProduct(rot_velo, B(1,:))) * &
                 CrossProduct(rot_velo,B(1,:)))*Basis(p)*s
           END IF
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
           IF (HasVelocity) THEN
             Coeff = Coeff + HarmPowerCoeff * (SUM( MATMUL( REAL(CMat_ip(1:3,1:3)), CrossProduct(rot_velo, B(1,:)) ) * &
               CrossProduct(rot_velo, B(1,:)) ) * Basis(p) * s - &
               SUM( MATMUL( AIMAG(CMat_ip(1:3,1:3)), CrossProduct(rot_velo, B(2,:)) ) * &
               CrossProduct(rot_velo, B(1,:)) ) * Basis(p) * s + &
               SUM( MATMUL( AIMAG(CMat_ip(1:3,1:3)), CrossProduct(rot_velo, B(1,:)) ) * &
               CrossProduct(rot_velo, B(2,:)) ) * Basis(p) * s + &
               SUM( MATMUL( REAL(CMat_ip(1:3,1:3)), CrossProduct(rot_velo, B(2,:)) ) * &
               CrossProduct(rot_velo, B(2,:)) ) * Basis(p) * s)
           END IF
         END IF

         IF(ALLOCATED(BodyLoss)) BodyLoss(3,BodyId) = BodyLoss(3,BodyId) + Coeff
         Power = Power + Coeff
         IF ( ASSOCIATED(JH) .OR. ASSOCIATED(EL_JH) .OR. ASSOCIATED(NJH) ) THEN
           FORCE(p,k+1) = FORCE(p,k+1) + Coeff
           k = k+1
         END IF
         
         

         !-------------------------------------------------
         ! Compute a loss estimate for cos and sin modes:
         !-------------------------------------------------
         IF (LossEstimation) THEN
           LossCoeff = ListGetFun( Material,'Harmonic Loss Linear Coefficient',Freq,Found ) 
           LossCoeff2 = ListGetFun( Material,'Harmonic Loss Quadratic Coefficient',Freq,Found ) 
           ! No losses to add if loss coefficient is not given
           IF( Found ) THEN
             DO l=1,2
               ValAtIP = SUM( B(l,1:3) ** 2 )
               Coeff = s * Basis(p) * LossCoeff * ( Freq ** FreqPower ) * ( ValAtIp ** FieldPower )
               Coeff2 = s * Basis(p) * LossCoeff2 * ( Freq ** FreqPower2 ) * ( ValAtIp ** FieldPower2 )
               ComponentLoss(1,l) = ComponentLoss(1,l) + Coeff
               BodyLoss(1,BodyId) = BodyLoss(1,BodyId) + Coeff 
               ComponentLoss(2,l) = ComponentLoss(2,l) + Coeff2
               BodyLoss(2,BodyId) = BodyLoss(2,BodyId) + Coeff2
             END DO
           ELSE
             Coeff = 0.0_dp
             Coeff2 = 0.0_dp
           END IF

           IF ( ASSOCIATED(ML) .OR. ASSOCIATED(EL_ML) ) THEN
             FORCE(p,k+1) = FORCE(p,k+1) + Coeff
             k = k + 1
           END IF
           IF ( ASSOCIATED(ML2) .OR. ASSOCIATED(EL_ML2) ) THEN
             FORCE(p,k+1) = FORCE(p,k+1) + Coeff2
             k = k + 1
           END IF
         END IF

         IF ( ASSOCIATED(MST).OR.ASSOCIATED(EL_MST)) THEN
           IF ( Vdofs==1 ) THEN
             DO l=1,3
               DO m=l,3
                 ST(l,m)=PR_ip*E(1,l)*E(1,m)+R_ip*B(1,l)*B(1,m)
               END DO
               ST(l,l)=ST(l,l)-(PR_ip*SUM(E(1,:)**2)+R_ip*SUM(B(1,:)**2))/2
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
             FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3) + s*NF_ip(p,1:3)
           ELSE
             FORCE(p,k+1:k+3) = FORCE(p,k+1:k+3) + 0.5*s*NF_ip(p,1:3)
           END IF
           k = k + 3
         END IF
       END DO ! p
     END DO ! j


     IF(NodalFields) THEN
       CALL DefaultUpdateEquations( MASS,Force(:,1))
       Fsave => Solver % Matrix % RHS
       DO l=1,k
         Solver % Matrix % RHS => GForce(:,l)
         CALL DefaultUpdateForce(Force(:,l))
       END DO
       Solver % Matrix % RHS => Fsave
     END IF

     IF(ElementalFields) THEN
       dofs = 0
       CALL LUdecomp(MASS,n,pivot)
       CALL LocalSol(EL_MFD,  3*vdofs, n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_MFS,  3*vdofs, n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_VP,   3*vdofs, n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_EF,   3*vdofs, n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_CD,   3*vdofs, n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_JXB,  3*vdofs, n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_FWP,  1*vdofs, n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_JH,   1, n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_ML,   1, n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_ML2,  1, n, MASS, FORCE, pivot, Dofs)
       CALL LocalSol(EL_MST,  6*vdofs, n, MASS, FORCE, pivot, Dofs)

       ! This is a nodal quantity
       CALL LocalCopy(EL_NF, 3, n, FORCE, Dofs)
     END IF
   END DO
   
   
   ! Assembly of the face terms:
   !----------------------------
   IF (GetLogical(SolverParams,'Discontinuous Galerkin',Found)) THEN
     IF (GetLogical(SolverParams,'Average Within Materials',Found)) THEN
       FORCE = 0.0_dp
       CALL AddLocalFaceTerms( MASS, FORCE(:,1) )
     END IF
   END IF

   
   IF(NodalFields) THEN
     Fsave => Solver % Matrix % RHS
     DOFs = 0
     CALL GlobalSol(MFD,  3*vdofs, Gforce, Dofs)
     CALL GlobalSol(MFS,  3*vdofs, Gforce, Dofs)
     CALL GlobalSol(VP ,  3*vdofs, Gforce, Dofs)
     CALL GlobalSol(EF,   3*vdofs, Gforce, Dofs)
     CALL GlobalSol(CD,   3*vdofs, Gforce, Dofs)
     CALL GlobalSol(JXB,  3*vdofs, Gforce, Dofs)
     CALL GlobalSol(FWP,  1*vdofs, Gforce, Dofs)
     
     ! Nodal heating directly uses the loads 
     IF (ASSOCIATED(NJH)) THEN
       NJH % Values = Gforce(:,dofs+1)
       ! Update the dofs only if it not used as the r.h.s. for the following field
       IF(.NOT. ASSOCIATED(JH) ) dofs = dofs + 1
     END IF
     CALL GlobalSol(JH ,  1      , Gforce, Dofs)

     CALL GlobalSol(ML ,  1      , Gforce, Dofs)
     CALL GlobalSol(ML2,  1      , Gforce, Dofs)
     CALL GlobalSol(MST,  6*vdofs, Gforce, Dofs)
     IF (ASSOCIATED(NF)) THEN
       DO i=1,3
         dofs = dofs + 1
         NF % Values(i::3) = Gforce(:,dofs)
       END DO
     END IF
     Solver % Matrix % RHS => Fsave
   END IF


   ! Lump componentwise forces and torques. 
   ! Prefer DG nodal force variable if air gap is present

   ! Warn if user has air gaps and no "nodal force e"
   HaveAirGap = ListCheckPresentAnyBC( Model, 'Air Gap Length' ) 
   UseElementalNF = ASSOCIATED( EL_NF ) .AND. ( .NOT. ASSOCIATED( NF ) .OR. HaveAirGap )
   
    
   IF( UseElementalNF ) THEN

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

         DO i=1,3
           CALL ListAddConstReal( CompParams,'res: magnetic force '//TRIM(I2S(i)), LumpedForce(i) )
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
     DO j=1,Model % NumberOfComponents
       CompParams => Model % Components(j) % Values

       IF( ListGetLogical( CompParams,'Calculate Magnetic Force', Found ) ) THEN 

         ! fail if there is no nodal force variable available
         IF (.NOT. ASSOCIATED(NF)) THEN
           CALL Warn('MagnetoDynamicsCalcFields','Unable to calculated lumped &
             &forces because nodal forces are not present. Use keyword &
             &"Calculate Nodal Forces = true" in MagnetoDynamicsCalcFields solver.')
           EXIT
         END IF

         CALL ComponentNodalForceReduction(Model, Mesh, CompParams, NF, &
           Force = LumpedForce )

         WRITE( Message,'(A,3ES15.6)') 'Magnetic force reduced: > '&
           //TRIM(ListGetString(CompParams,'Name'))//' < :', LumpedForce
         CALL Info('MagnetoDynamicsCalcFields',Message,Level=6)           

         DO i=1,3
           CALL ListAddConstReal( CompParams,'res: magnetic force '//TRIM(I2S(i)), LumpedForce(i) )
         END DO

       END IF

       IF( ListGetLogical( CompParams,'Calculate Magnetic Torque', Found ) ) THEN 

         ! fail if there is no nodal force variable available
         IF (.NOT. ASSOCIATED(NF)) THEN
           CALL Warn('MagnetoDynamicsCalcFields','Unable to calculated lumped &
             &forces because nodal forces are not present. Use keyword &
             &"Calculate Nodal Forces = true" in MagnetoDynamicsCalcFields solver.')
           EXIT 
         END IF

         ! Warn if user has air gaps and no "nodal force e" is available
         IF ( HaveAirGap ) THEN
           CALL Warn('MagnetoDynamicsCalcFields', 'Cannot calculate air gap &
             &forces correctly because elemental field "Nodal Force e" is not &
             &present.')
         END IF

         CALL ComponentNodalForceReduction(Model, Mesh, CompParams, NF, &
           Torque = val )

         WRITE( Message,'(A,ES15.6)') 'Magnetic torque reduced: > '&
           //TRIM(ListGetString(CompParams,'Name'))//' < :', val
         CALL Info('MagnetoDynamicsCalcFields',Message,Level=6)           

         CALL ListAddConstReal( CompParams,'res: magnetic torque', val )
       END IF
     END DO
   END IF


   ! Perform parallel reductions 
   Power  = ParallelReduction(Power)
   Energy = ParallelReduction(Energy)
  
   IF (LossEstimation) THEN
     DO j=1,2
       DO i=1,2
         ComponentLoss(j,i) = ParallelReduction(ComponentLoss(j,i)) 
       END DO
     END DO

     DO j=1,3
       DO i=1,Model % NumberOfBodies
         BodyLoss(j,i) = ParallelReduction(BodyLoss(j,i))
       END DO
       TotalLoss(j) = SUM( BodyLoss(j,:) )
     END DO
   END IF

   
   WRITE(Message,*) 'Eddy current power: ', Power
   CALL Info( 'MagnetoDynamicsCalcFields', Message )
   CALL ListAddConstReal( Model % Simulation, 'res: Eddy current power', Power )

   WRITE(Message,*) '(Electro)Magnetic Field Energy: ', Energy
   CALL Info( 'MagnetoDynamicsCalcFields', Message )
   CALL ListAddConstReal(Model % Simulation,'res: Magnetic Field Energy',Energy)
   IF(ALLOCATED(Gforce)) DEALLOCATE(Gforce)
   DEALLOCATE( MASS,FORCE,Tcoef,RotM, R_t )

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
       WRITE( Message,'(A,ES12.3)') 'Loss for cos mode: ', ComponentLoss(k,1)
       CALL Info('MagnetoDynamicsCalcFields', Message, Level=6 )
       WRITE( Message,'(A,ES12.3)') 'Loss for sin mode: ', ComponentLoss(k,2)
       CALL Info('MagnetoDynamicsCalcFields', Message, Level=6 )
       WRITE( Message,'(A,ES12.3)') 'Total loss: ',TotalLoss(k)
       CALL Info('MagnetoDynamicsCalcFields',Message, Level=5 )
     END DO

     DO k=1,3
       IF( TotalLoss(k) > TINY( TotalLoss(k) ) ) CYCLE
       IF( k == 1 ) THEN
         CALL Info('MagnetoDynamicsCalcFields','Harmonic Loss Linear by bodies',Level=6)
       ELSE IF( k == 2 ) THEN
         CALL Info('MagnetoDynamicsCalcFields','Harmonic Loss Quadratic by bodies',Level=6)
       ELSE
         CALL Info('MagnetoDynamicsCalcFields','Joule Loss by bodies',Level=6)
       END IF

       DO j=1,Model % NumberOfBodies
         IF( BodyLoss(k,j) < TINY( TotalLoss(k) ) ) CYCLE
         WRITE( Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyLoss(k,j)
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
         OPEN (10, FILE=LossFile)
         WRITE( 10,'(A)')  '!body_id   harmonic(1)      harmonic(2)      joule'
         DO j=1,Model % NumberOfBodies
           IF( SUM(BodyLoss(1:3,j)) < TINY( TotalLoss(1) ) ) CYCLE
           WRITE( 10,'(I0,T10,3ES17.9)') j, BodyLoss(1:3,j)
         END DO
         CALL Info('MagnetoDynamicsCalsFields', &
             'Harmonic loss for bodies was saved to file: '//TRIM(LossFile),Level=6 )
         CLOSE(10)
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
       WRITE (Message,'(A)') 'res: Group '//TRIM(I2S(j))//' torque'
       CALL ListAddConstReal(Model % Simulation, TRIM(Message), Torque(i))
       WRITE (Message,'(A,F0.8)') 'Torque Group '//TRIM(I2S(j))//' torque:', Torque(i)
       CALL Info( 'MagnetoDynamicsCalcFields', Message)
     END DO

     CALL NodalTorqueDeprecated(TorqueDeprecated, Found)
     IF (Found) THEN
       WRITE(Message,*) 'Torque over defined bodies', TorqueDeprecated
       CALL Info( 'MagnetoDynamicsCalcFields', Message )
       CALL Warn( 'MagnetoDynamicsCalcFields', 'Keyword "Calculate Torque over body" is deprecated, use Torque Groups instead')
       CALL ListAddConstReal(Model % Simulation, 'res: x-axis torque over defined bodies', TorqueDeprecated(1))
       CALL ListAddConstReal(Model % Simulation, 'res: y-axis torque over defined bodies', TorqueDeprecated(2))
       CALL ListAddConstReal(Model % Simulation, 'res: z-axis torque over defined bodies', TorqueDeprecated(3))
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
    Flux(1) = ParallelReduction(Flux(1))
    Flux(2) = ParallelReduction(Flux(2))
    Area = ParallelReduction(Area)

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


CONTAINS

!-------------------------------------------------------------------
  SUBROUTINE SumElementalVariable(Var, Values, BodyId, Additive)
!-------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Variable_t), POINTER :: Var
    REAL(KIND=dp), OPTIONAL, TARGET :: Values(:)
    INTEGER, OPTIONAL :: BodyId
    LOGICAL, OPTIONAL :: Additive

    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), ALLOCATABLE :: NodeSum(:)
    INTEGER :: n, j, k, l, nodeind, dgind, bias
    LOGICAL, ALLOCATABLE :: AirGapNode(:)
    REAL(KIND=dp), POINTER :: ValuesSource(:)


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
        IF(PRESENT(BodyID) .AND. Element % BodyID /= BodyID) CYCLE
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
        IF(PRESENT(BodyID) .AND. Element % BodyID /= BodyID) CYCLE
        DO l=1,Element%TYPE%NumberofNodes
          nodeind = Element % NodeIndexes(l)
          dgind = Var % Perm(Element % DGIndexes(l))
          IF( dgind > 0 ) THEN
            IF (PRESENT(Additive) .AND. Additive) THEN
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
      NF_ip_l(27,3), NF_ip_r(27,3)
    TYPE(Element_t), POINTER :: LeftParent, RightParent, BElement
    TYPE(Nodes_t), SAVE :: LPNodes, RPNodes
    REAL(KIND=dp) :: F(3,3)
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

      IF (SecondOrder) THEN
        IP = GaussPoints(BElement, EdgeBasis=dim==3, PReferenceElement=PiolaVersion, EdgeBasisDegree=EdgeBasisDegree)
      ELSE
        IP = GaussPoints(BElement, EdgeBasis=dim==3, PReferenceElement=PiolaVersion)
      END IF

      
      DO j = 1,IP % n
        s = IP % s(j)

        IF ( PiolaVersion ) THEN
          stat = EdgeElementInfo( BElement, Nodes, IP % U(j), IP % V(j), IP % W(j), &
            F = F, DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, RotBasis = RotWBasis, &
            dBasisdx=dBasisdx, BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)
        ELSE
          stat = ElementInfo( BElement, Nodes, IP % U(j), IP % V(j), &
            IP % W(j), detJ, Basis, dBasisdx )

          CALL GetEdgeBasis(BElement, WBasis, RotWBasis, Basis, dBasisdx)
        END IF

        R_ip = SUM( Basis(1:n)/(mu0*AirGapMu(1:n)) )
        GapLength_ip = SUM( Basis(1:n)*GapLength(1:n) )

        s = s * detJ

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
            B(k,:) = normal*sum( SOL(k,np+1:nd)* RotWBasis(1:nd-np,3) )
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

        Energy = Energy + GapLength_ip*s*0.5*R_ip*B2

        DO p=1,n
          IF(HasLeft) LeftFORCE(LeftMap(p), 1:3) = LeftFORCE(LeftMap(p), 1:3) + s*NF_ip_l(p,1:3)
          IF(HasRight) RightFORCE(RightMap(p), 1:3) = RightFORCE(RightMap(p), 1:3) + s*NF_ip_r(p,1:3)
        END DO
      END DO ! Integration points
      
      IF(ElementalFields) THEN
        IF(HasLeft) CALL LocalCopy(EL_NF, 3, n_lp, LeftFORCE, 0, UElement=LeftParent, uAdditive=.TRUE.)
        IF(HasRight) CALL LocalCopy(EL_NF, 3, n_rp, RightFORCE, 0, UElement=RightParent, uAdditive=.TRUE.)
      END IF
    END DO ! Boundary elements
    
    DEALLOCATE( LeftFORCE, RightFORCE, RightMap, LeftMap, AirGapForce )
!-------------------------------------------------------------------
  END SUBROUTINE CalcBoundaryModels
!-------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE NodalTorqueDeprecated(T, FoundOne)
!------------------------------------------------------------------------------
   IMPLICIT NONE
   REAL(KIND=dp), INTENT(OUT) :: T(3)
   LOGICAL, INTENT(OUT) :: FoundOne
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: P(3), F(3)
   TYPE(Element_t), POINTER :: Element
   TYPE(Variable_t), POINTER :: CoordVar
   LOGICAL :: VisitedNode(Mesh % NumberOfNodes)
   INTEGER :: pnodal, nnt, ElemNodeDofs(27), ndofs, globalnode, m, n
   LOGICAL :: ONCE=.TRUE., Found
   
   VisitedNode = .FALSE.
   FoundOne = .FALSE.

   DO n=1,size(Model % bodies)
     IF(GetLogical(Model % bodies(n) % Values, 'Calculate Torque over body', FoundOne)) EXIT
   END DO
   IF(.not. FoundOne) RETURN
   T = 0._dp
   P = 0._dp

   DO pnodal=1,GetNOFActive()
     Element => GetActiveElement(pnodal)
     IF(GetLogical(GetBodyParams(Element), 'Calculate Torque over body', Found)) THEN
       ndofs = GetElementDOFs(ElemNodeDofs)
       DO nnt=1,ndofs
         globalnode = ElemNodeDofs(nnt)
         IF (.NOT. VisitedNode(globalnode)) THEN
           F(1) = NF % Values( 3*(NF % Perm((globalnode))-1) + 1)
           F(2) = NF % Values( 3*(NF % Perm((globalnode))-1) + 2)
           F(3) = NF % Values( 3*(NF % Perm((globalnode))-1) + 3)
           P(1) = Mesh % Nodes % x(globalnode)
           P(2) = Mesh % Nodes % y(globalnode)
           P(3) = Mesh % Nodes % z(globalnode)
           T(1) = T(1) + P(2)*F(3)-P(3)*F(2)
           T(2) = T(2) + P(3)*F(1)-P(1)*F(3)
           T(3) = T(3) + P(1)*F(2)-P(2)*F(1)
           VisitedNode(globalnode) = .TRUE.
         END IF
       END DO ! nnt
     END IF
   END DO ! pnodal
   T(1) = ParallelReduction(T(1))
   T(2) = ParallelReduction(T(2))
   T(3) = ParallelReduction(T(3))
!------------------------------------------------------------------------------
 END SUBROUTINE NodalTorqueDeprecated
!------------------------------------------------------------------------------

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
   IF(k .eq. 0) RETURN

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
       IF (nrm .EQ. 0._dp) THEN
         CALL Warn('MagnetoDynamicsCalcFields',&
             'Axis for the torque group '//TRIM(I2S(k))//' is a zero vector')
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
     IF(.not. Found) CYCLE
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
           IF (LocalGroups(ng) .gt. num_origins) THEN
             origin = 0._dp
           ELSE
             origin = origins(LocalGroups(ng),1:3)
           END IF
           IF (LocalGroups(ng) .gt. num_axes) THEN
             axisvector = 0._dp
             axisvector(3) = 1._dp
           ELSE
             axisvector = axes(LocalGroups(ng), 1:3)
           END IF
           v1 = P - origin
           v1 = (1 - sum(axisvector*v1))*v1
           v2 = CrossProduct(v1,F)
           T(LocalGroups(ng)) = T(LocalGroups(ng)) + sum(axisvector*v2)
         END IF
       END DO 
     END DO 

   END DO
   DO ng=1,size(TorqueGroups)
     T(ng) = ParallelReduction(T(ng))
   END DO

!------------------------------------------------------------------------------
  END SUBROUTINE NodalTorque
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GlobalSol(Var, m, b, dofs )
!------------------------------------------------------------------------------
   IMPLICIT NONE
   REAL(KIND=dp), TARGET CONTIG :: b(:,:)
   INTEGER :: m, dofs
   TYPE(Variable_t), POINTER :: Var
!------------------------------------------------------------------------------
   INTEGER :: i
!------------------------------------------------------------------------------
   IF(.NOT. ASSOCIATED(var)) RETURN

   CALL Info('MagnetoDynamicsCalcFields','Solving for field: '//TRIM(Var % Name),Level=6)
   
   DO i=1,m
     dofs = dofs+1
     Solver % Matrix % RHS => b(:,dofs)
     Solver % Variable % Values=0
     Norm = DefaultSolve()
     var % Values(i::m) = Solver % Variable % Values
   END DO
!------------------------------------------------------------------------------
 END SUBROUTINE GlobalSol
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE LocalSol(Var, m, n, A, b, pivot, dofs )
!------------------------------------------------------------------------------
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: Var
   INTEGER :: pivot(:), m,n,dofs
   REAL(KIND=dp) :: b(:,:), A(:,:)
!------------------------------------------------------------------------------
   INTEGER :: ind(n), i
   REAL(KIND=dp) :: x(n)
!------------------------------------------------------------------------------
   IF(.NOT. ASSOCIATED(var)) RETURN

   IF( ANY( Var % Perm( Element % DGIndexes(1:n) ) <= 0 ) ) THEN
     PRINT *,'size',SIZE( Var % Perm ), MAXVAL( Element % DGIndexes(1:n))
     PRINT *,'Perm zero:',m,n,dofs,Var % Perm( Element % DGIndexes(1:n) )
     PRINT *,'size values',SIZE(Var % Values)
     PRINT *,'Element index:',Element % ElementIndex
     PRINT *,'Element indexes:',Element % NodeIndexes
     STOP
   END IF

   ind(1:n) = Var % DOFs*(Var % Perm(Element % DGIndexes(1:n))-1)
 
   DO i=1,m
      dofs = dofs+1
      x = b(1:n,dofs)
      CALL LUSolve(n,MASS,x,pivot)
      Var % Values(ind(1:n)+i) = x(1:n)
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
   REAL(KIND=dp), OPTIONAL :: Values(:)
   LOGICAL, OPTIONAL :: uAdditive
!------------------------------------------------------------------------------
   INTEGER :: ind(n), i
   LOGICAL :: Additive
!------------------------------------------------------------------------------
   IF(.NOT. ASSOCIATED(var)) RETURN
   IF(PRESENT(UElement)) THEN
     ind(1:n) = Var % DOFs*(Var % Perm(UElement % DGIndexes(1:n))-1)
   ELSE
     ind(1:n) = Var % DOFs*(Var % Perm(Element % DGIndexes(1:n))-1)
   END IF
   
   IF(PRESENT(uAdditive)) THEN
     Additive = uAdditive
   ELSE
     Additive = .FALSE.
   END IF

   dofs = bias
   IF(PRESENT(Values)) THEN
     DO i=1,m
       dofs = dofs+1
       IF(Additive) THEN
         Values(ind(1:n)+i) = Values(ind(1:n)+i) + b(1:n,dofs)
       ELSE
         Values(ind(1:n)+i) = b(1:n,dofs)
       END IF
     END DO
   ELSE
     DO i=1,m
       dofs = dofs+1
       IF(Additive) THEN
         Var % Values(ind(1:n)+i) = Var % Values(ind(1:n)+i) + b(1:n,dofs)
       ELSE
         Var % Values(ind(1:n)+i) = b(1:n,dofs)
       END IF
     END DO
   END IF
!------------------------------------------------------------------------------
 END SUBROUTINE LocalCopy
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE GetElementRotM(Element,RotM,n)
!------------------------------------------------------------------------------
   IMPLICIT NONE
   TYPE(Element_t) :: Element
   INTEGER :: k, l, m, j, n
   REAL(KIND=dp) :: RotM(3,3,n)
   INTEGER, PARAMETER :: ind1(9) = [1,1,1,2,2,2,3,3,3]
   INTEGER, PARAMETER :: ind2(9) = [1,2,3,1,2,3,1,2,3]
   TYPE(Variable_t), POINTER, SAVE :: RotMvar
   LOGICAL, SAVE :: visited = .FALSE.
 

   IF(.NOT. visited) THEN
     visited = .TRUE.
     RotMvar => VariableGet( Mesh % Variables, 'RotM E')
     IF(.NOT. ASSOCIATED(RotMVar)) THEN
       CALL Fatal('GetElementRotM','RotM E variable not found')
     END IF
   END IF

   RotM = 0._dp
   DO j = 1, n
     DO k=1,RotMvar % DOFs
       RotM(ind1(k),ind2(k),j) = RotMvar % Values( &
             RotMvar % DOFs*(RotMvar % Perm(Element % DGIndexes(j))-1)+k)
     END DO
   END DO

!------------------------------------------------------------------------------
 END SUBROUTINE GetElementRotM
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
                  IP % W(j), detJ, Basis, dBasisdx )
         CALL GetEdgeBasis(Element, WBasis, RotWBasis, Basis, dBasisdx)
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

