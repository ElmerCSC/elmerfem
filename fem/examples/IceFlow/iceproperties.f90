!/*****************************************************************************
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
!*
! ******************************************************************************
! *
! *                    Author:  Juha Ruokolainen
! *
! *                    Address: CSC - IT Center for Science Ltd.
! *                                Keilaranta 14, P.O. BOX 405
! *                                  02101 Espoo, Finland
! *
! *
!/******************************************************************************
! *
! *  iceproperties.f90  material parameters and physical models for ice flow
! *
! *
! *       Module Author:           Thomas Zwinger
! *       Address:                 CSC - IT Center for Science Ltd.
! *                                Keilaranta 14, P.O. BOX 405
! *                                  02101 Espoo, Finland
! *                                  Tel. +358 0 457 2723
! *                                Telefax: +358 0 457 2183
! *                                EMail: Thomas.Zwinger@csc.fi
! *
! *       Modified by:             Thomas Zwinger
! *
! *       Date of modification: 11/11/2005
! *
! *****************************************************************************/
!
!
!

RECURSIVE SUBROUTINE getNetBoundaryHeatflux( Model,Solver,Timestep,TransientSimulation)
  USE DefUtils
  USE Materialmodels
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL :: TransientSimulation
  REAL(KIND=dp) :: Timestep
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: PointerToSolver
  TYPE(Nodes_t) :: Nodes
  TYPE(Element_t),POINTER :: Element
  TYPE(Variable_t), POINTER :: HFSol, TempSol
  TYPE(ValueList_t), POINTER :: Equation,Material,SolverParams,BodyForce,BC,Constants
  CHARACTER(LEN=MAX_NAME_LEN) :: TempName, SolverName
  INTEGER :: i, j, k, t, N, BN,  M, DIM, istat, material_id
  INTEGER, POINTER :: TempPerm(:), HFPerm(:),BoundaryReorder(:)
  INTEGER, ALLOCATABLE :: contributedElements(:)
  REAL(KIND=dp) ::  U, V, W, SqrtElementMetric
  REAL(KIND=dp), ALLOCATABLE ::  Basis(:),dBasisdx(:,:), ddBasisddx(:,:,:),&
       LatentHeat(:), HeatConductivity(:), Density(:), heatFlux(:,:)
  REAL(KIND=dp), POINTER :: HF(:), Temp(:), BoundaryNormals(:,:), BoundaryTangent1(:,:), BoundaryTangent2(:,:)
  LOGICAL :: GotIt, AllocationsDone=.FALSE., stat

  SAVE AllocationsDone, SolverName, DIM,  heatFlux, &
       Basis, dBasisdx, ddBasisddx, &
       LatentHeat, HeatConductivity, Density, contributedElements, &
       BoundaryNormals, BoundaryTangent1, BoundaryTangent2, BoundaryReorder, BN


  WRITE(SolverName, '(A)') 'iceproperties (getNetBoundaryHeatflux))'
  !-----------------------------------------------------------------------
  ! get solver variable
  !-----------------------------------------------------------------------
  HFSol => Solver % Variable
  IF (.NOT.ASSOCIATED(HFSol)) THEN
     CALL FATAL(SolverName,'No variable associated')
  END IF
  HFPerm  => HFSol % Perm
  HF => HFSol % Values

  
  HF = 0.0D00


  !-----------------------------------------------------------------------
  ! Allocations
  !-----------------------------------------------------------------------
  IF ( .NOT.AllocationsDone .OR. Solver % Mesh % Changed ) THEN
     DIM = CoordinateSystemDimension()
     N = Solver % Mesh % MaxElementNodes
     M = Model % Mesh % NumberOfNodes        
     IF ( AllocationsDone ) &
          DEALLOCATE( &
          Basis, &
          dBasisdx, &
          ddBasisddx, &
          LatentHeat, &
          HeatConductivity, &
          Density, &
          heatFlux,&
          Nodes % x, &
          Nodes % y, &
          Nodes % z, &
          contributedElements, &
          BoundaryReorder, &
          BoundaryNormals, &
          BoundaryTangent1, &
          BoundaryTangent2)
     CALL CheckNormalTangentialBoundary( Model, &
          'Basal Melting', BN, &
          BoundaryReorder, BoundaryNormals, BoundaryTangent1, &
          BoundaryTangent2, DIM )
      WRITE(Message,'(A,i6)') &
          'Number of boundary nodes on boundaries associated with melting:',&
           BN
      CALL INFO(SolverName,Message,Level=3)
      CALL AverageBoundaryNormals(Model, &
           'Basal Melting', BN, &
           BoundaryReorder, BoundaryNormals, BoundaryTangent1, &
           BoundaryTangent2, DIM )
     ALLOCATE( &
          Basis(N), &
          dBasisdx(N,3), &
          ddBasisddx(N,3,3), &
          LatentHeat(N), &
          HeatConductivity(N), &
          Density(N), &
          heatFlux(M,3),& 
          Nodes % x(N), &
          Nodes % y(N), &
          Nodes % z(N), &
          contributedElements(M), &
          STAT = istat)
     AllocationsDone = .TRUE.
  END IF
  !-------------------------------------------------------------------------
  ! Calculate element-wise heat flux contributions for each point of element
  !-------------------------------------------------------------------------
  contributedElements = 0
  heatFlux = 0.0D00  

  DO t=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     Model % CurrentElement => Element
     N = Element % Type % NumberOfNodes
     CALL GetElementNodes( Nodes )
     !-----------------------------------------------------------------------
     ! get material pointer
     !-----------------------------------------------------------------------
     Material => GetMaterial()
     IF (.NOT.ASSOCIATED(Material)) THEN
        WRITE (Message,'(A,I3)') 'No Material found for boundary element no. ', t
        CALL FATAL(SolverName,Message)
     ELSE
        material_id = GetMaterialId( Element, GotIt)
        IF(.NOT.GotIt) THEN
           WRITE (Message,'(A,I3)') 'No Material ID found for boundary element no. ', t
           CALL FATAL(SolverName,Message)
        END IF
     END IF
     !-----------------------------------------------------------------------
     ! get temperature variable
     !-----------------------------------------------------------------------
     TempName =  GetString(Material,'Temperature Name', GotIt)
     IF (.NOT.GotIt) THEN
        CALL FATAL(SolverName,'No Temperature Name found')
     ELSE
        WRITE(Message,'(a,a)') 'Variable Name for temperature: ', TempName
        CALL INFO(SolverName,Message,Level=12)
     END IF
     TempSol => VariableGet( Solver % Mesh % Variables, TRIM(TempName) )
     IF ( ASSOCIATED( TempSol ) ) THEN
        TempPerm => TempSol % Perm
        Temp => TempSol % Values
     ELSE
        WRITE(Message, '(A,A)') 'Could not find temperature field variable ',  TRIM(TempName)
        CALL FATAL(SolverName,Message)
     END IF 
     !-----------------------
     ! get material parameter
     !-----------------------
     HeatConductivity(1:N) =  ListGetReal( Material,  TRIM(TempName) // &
          ' Heat Conductivity', n, Element % NodeIndexes, GotIt )
     IF (.NOT.GotIt) THEN
        HeatConductivity = 0.0D00
        WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', TRIM(TempName) // &
          ' Heat Conductivity', '< not found for element ', t, ' material ', material_id
        CALL WARN(SolverName,Message)
     END IF
     !--------------------------
     ! loop all nodes in element
     !--------------------------
     DO i=1,N
        U = Element % Type % NodeU(i)
        V = Element % Type % NodeV(i)
        W = Element % Type % NodeW(i)
        stat = ElementInfo( Element,Nodes,U,V,W,SqrtElementMetric, &
             Basis,dBasisdx,ddBasisddx,.FALSE.,.FALSE. )
        !--------------------
        ! loop all dimensions
        !--------------------
        DO j=1,DIM
           DO k=1,N
              heatFlux(Element % NodeIndexes(i),j) = & 
                   heatFlux(Element % NodeIndexes(i),j) - HeatConductivity(k)*dBasisdx(k,j)*Temp(TempPerm(k))
           END DO
        END DO
        contributedElements(Element % NodeIndexes(i)) = contributedElements(Element % NodeIndexes(i)) + 1
     END DO
  END DO
  
  !----------------------------------------------------------
  ! averaged internal flux for each point of element
  !----------------------------------------------------------
  DO i=1,Solver % Mesh % NumberOfNodes
     IF (contributedElements(i) > 0) THEN
        heatFlux(i,1:DIM) = heatFlux(i,1:DIM)/contributedElements(i)
        j = BoundaryReorder(i)
        IF (j>0) THEN
           HF(HFPerm(i)) = SUM(heatFlux(i,1:DIM)*BoundaryNormals(j,1:DIM))
        ELSE
!           HF(HFPerm(i)) = SQRT(SUM(heatFlux(i,1:DIM)*heatFlux(i,1:DIM)))
           HF(HFPerm(i)) = 0.0D00
        END IF
     ELSE
        WRITE(Message,'(A,i6,A)') 'Node no. ',i,' appears not to be part of any element.'
        CALL WARN(SolverName,Message)
        HF(HFPerm(i)) = 0.0
     END IF     
  END DO

END SUBROUTINE getNetBoundaryHeatflux

!**************************************************************************
!*
!*  basal melting velocity  as a function of heat flux (provided externaly)
!*
!**************************************************************************
FUNCTION getBasalMeltingVelocity( Model, Node, HeatFlux ) RESULT(basalMeltingvelocity)

!-----------------------------------------------------------
  USE DefUtils
  USE SolverUtils
!-----------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------
  !external variables
  TYPE(Model_t), TARGET :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: HeatFlux, basalMeltingvelocity

  !internal variables
  TYPE(ValueList_t), POINTER :: ParentMaterial, BC
  TYPE(Element_t), POINTER :: BoundaryElement, ParentElement 
  TYPE(Variable_t), POINTER :: VarHomTemp
  INTEGER :: i, N, DIM, NBoundary, NParent, BoundaryElementNode, ParentElementNode, body_id, other_body_id, material_id, istat
  REAL(KIND=dp), ALLOCATABLE :: LatentHeat(:), Density(:), HomTemp(:),ExternalHF(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: TempName, SolverName
  LOGICAL ::  FirstTime = .TRUE., GotIt, stat

  SAVE FirstTime, LatentHeat, Density, HomTemp, DIM, N, ExternalHF

  WRITE(SolverName, '(A)') 'iceproperties (getBasalMeltingVelocity))'
  !--------------------------------
  ! Allocations
  !--------------------------------
  IF ( FirstTime .OR. Model % Mesh % Changed) THEN
      DIM = CoordinateSystemDimension()
      N = Model % MaxElementNodes 

      IF (.NOT.FirstTime) THEN
         DEALLOCATE( &
              LatentHeat,&        
              Density,&
              HomTemp,&
              ExternalHF)
      END IF
      ALLOCATE( &
           LatentHeat( N ),&
           Density( N ),&
           HomTemp( N ),&
           ExternalHF( N ),&
           STAT = istat)
      IF (istat /= 0) THEN
         CALL FATAL(SolverName,'Allocations failed')
      ELSE 
         CALL INFO(SolverName,'Allocations done')
      END IF

      FirstTime = .FALSE.
   END IF



  !-----------------------------------------------------------------
  ! get some information upon active boundary element and its parent
  !-----------------------------------------------------------------
  BoundaryElement => Model % CurrentElement
  IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN
     CALL FATAL(SolverName,'No boundary element found')
  END IF  
  NBoundary = BoundaryElement % Type % NumberOfNodes
  DO BoundaryElementNode=1,NBoundary
     IF (Node .EQ. BoundaryElement % NodeIndexes(BoundaryElementNode)) THEN
        GotIt = .TRUE.
        EXIT
     END IF
  END DO
  IF (.NOT.GotIt) THEN
     CALL WARN(SolverName,'Node not found in Current Element')
     basalMeltingVelocity = 0.0D00
     RETURN
  END IF
  other_body_id = BoundaryElement % BoundaryInfo % outbody
  IF (other_body_id < 1) THEN ! only one body in calculation
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
  ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
  END IF
  ! just to be on the save side, check again
  IF ( .NOT. ASSOCIATED(ParentElement) ) THEN
     WRITE(Message,'(A,I10,A)')&
          'Parent Element for Boundary element no. ',&
          BoundaryElement % ElementIndex, ' not found'
     CALL FATAL(SolverName,Message)
  END IF  

  body_id = ParentElement % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  ParentMaterial => Model % Materials(material_id) % Values
  IF ((.NOT. ASSOCIATED(ParentMaterial)) .OR. (.NOT. GotIt)) THEN
     WRITE(Message,'(A,I10,A,I10)')&
          'No material values found for body no ', body_id,&
          ' under material id ', material_id
     CALL FATAL(SolverName,Message)
  END IF 
  NParent = ParentElement % Type % NumberOfNodes
  DO ParentElementNode=1,NParent
     IF ( Node == ParentElement % NodeIndexes(ParentElementNode) ) EXIT
  END DO

  !-----------------------------------------------------------------------
  ! get temperature variable
  !-----------------------------------------------------------------------
  TempName =  GetString(ParentMaterial ,'Temperature Name', GotIt)
  IF (.NOT.GotIt) THEN
     CALL FATAL(SolverName,'No Temperature Name found')
  ELSE
     WRITE(Message,'(a,a)') 'Variable Name for temperature: ', TempName
     CALL INFO(SolverName,Message,Level=12)
  END IF

  !-------------------------
  ! Get material parameters
  !-------------------------
  Model % CurrentElement => ParentElement
  LatentHeat(1:NParent) = ListGetReal(ParentMaterial, 'Latent Heat', NParent, ParentElement % NodeIndexes, GotIt)
  IF (.NOT. GotIt) THEN
     CALL FATAL(SolverName,'No value for Latent Heat found')
  END IF
  Density(1:NParent) = ListGetReal( ParentMaterial, 'Density', NParent, ParentElement % NodeIndexes)
  IF (.NOT. GotIt) THEN
     CALL FATAL(SolverName,'No value for Density found')
  END IF
  !-------------------------
  ! Get boundary parameters
  !-------------------------
  Model % CurrentElement => BoundaryElement
  BC => GetBC()
  IF (.NOT.ASSOCIATED(BC)) THEN
     CALL FATAL(SolverName,'No Boundary Condition associated')
  ELSE
     ExternalHF(1:NBoundary) = GetReal(BC, TRIM(TempName) // ' Heat Flux', GotIt)
     IF (.NOT. GotIt) THEN
        WRITE(Message,'(a,a,a)') 'Keyword >', TRIM(TempName) // ' Heat Flux','< not found'
        CALL WARN(SolverName,Message)
        ExternalHF(1:NBoundary) = 0.0D00
     END IF
  END IF

  !---------------------------------
  ! Get homologous temperature field
  !---------------------------------
!  TempName =  GetString(Model % Constants ,'Temperature Name', GotIt)
!  IF (.NOT.GotIt) THEN
!     CALL FATAL(SolverName,'No Temperature Name found')
!  ELSE
!     WRITE(Message,'(a,a)') 'Variable Name for temperature: ', TempName
!     CALL INFO(SolverName,Message,Level=12)
!  END IF

  Model % CurrentElement => ParentElement
  VarHomTemp => VariableGet( Model % Variables, TRIM(TempName) // ' Homologous', .TRUE. )
  IF ( ASSOCIATED( VarHomTemp ) ) THEN
     DO i = 1, NParent
        HomTemp(i) = MIN(VarHomTemp % Values(VarHomTemp % Perm(ParentElement % NodeIndexes(i))),1.0D-06)
     END DO
  ELSE
    WRITE(Message,'(A,A,A)') 'Variable ', TRIM(TempName) // ' Homologous', ' not found'
     CALL INFO(SolverName,Message,Level=12)
     HomTemp(1:NParent) = -1.0D-06
  END IF
 
  !------------------------------------
  ! Get basal melting rate and velocity
  !------------------------------------
  IF (HomTemp(BoundaryElementNode) > -1.0D-03) THEN
     basalMeltingvelocity = MAX((HeatFlux + ExternalHF(BoundaryElementNode)) &
          /(LatentHeat(BoundaryElementNode) * Density(BoundaryElementNode)),1.0D-12)

  ELSE
     basalMeltingvelocity = 0.0D00
  END IF
  Model % CurrentElement =>  BoundaryElement
END FUNCTION getBasalMeltingVelocity

!*********************************************************************************************************************************
!*
!*  basal slip coefficient as a function of temperature
!*
!*********************************************************************************************************************************

FUNCTION basalSlip( Model, Node, Temperature ) RESULT(basalSlipCoefficient)
!-----------------------------------------------------------
  USE DefUtils
!-----------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------
  !external variables
  TYPE(Model_t), TARGET :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Temperature, basalSlipCoefficient
  !internal variables
  TYPE(Element_t), POINTER :: CurrentElementAtBeginning, BoundaryElement, ParentElement
  TYPE(ValueList_t), POINTER :: ParentMaterial, BC
  TYPE(Variable_t), POINTER :: varTemperature, varPressure
  INTEGER :: N, NBoundary, NParent, BoundaryElementNode, ParentElementNode, &
       i, DIM, other_body_id, body_id, material_id, istat, NSDOFs
  REAL(KIND=dp) :: TempHom, ThermalCoefficient, SlipCoefficientMax
  REAL(KIND=dp), ALLOCATABLE :: PressureMeltingPoint(:), TemperateSlipCoefficient(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: TempName
  LOGICAL ::  GotIt, stat, Jump=.FALSE.

  !---------------
  ! Initialization
  !---------------
  CurrentElementAtBeginning => Model % CurrentElement 
  N = Model % MaxElementNodes 
  ALLOCATE(TemperateSlipCoefficient(N),&  
       PressureMeltingPoint(N),&
       STAT = istat)
  IF (istat /= 0) THEN
     CALL FATAL('iceproperties (basalSlip)','Allocations failed')
  END IF

  TemperateSlipCoefficient = 1.0D30 ! high value - > no slip by default
  PressureMeltingPoint = 273.16D00

  !-----------------------------------------------------------------
  ! get some information upon active boundary element and its parent
  !-----------------------------------------------------------------
  BoundaryElement => Model % CurrentElement
  IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN
     CALL FATAL('iceproperties (basalMelting)','No boundary element found')
  END IF
  other_body_id = BoundaryElement % BoundaryInfo % outbody
  IF (other_body_id < 1) THEN ! only one body in calculation
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
  ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
  END IF
  ! just to be on the save side, check again
  IF ( .NOT. ASSOCIATED(ParentElement) ) THEN
     WRITE(Message,'(A,I10,A)')&
          'Parent Element for Boundary element no. ',&
          BoundaryElement % ElementIndex, ' not found'
     CALL FATAL('iceproperties (basalMelting)',Message)
  END IF  
  Model % CurrentElement => ParentElement
  body_id = ParentElement % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  ParentMaterial => Model % Materials(material_id) % Values
  IF ((.NOT. ASSOCIATED(ParentMaterial)) .OR. (.NOT. GotIt)) THEN
     WRITE(Message,'(A,I10,A,I10)')&
          'No material values found for body no ', body_id,&
          ' under material id ', material_id
     CALL FATAL('iceproperties (basalMelting)',Message)
  END IF
  ! number of nodes and node in elements
  NBoundary = BoundaryElement % Type % NumberOfNodes
  NParent = ParentElement % Type % NumberOfNodes
  DO BoundaryElementNode=1,Nboundary
     IF ( Node == BoundaryElement % NodeIndexes(BoundaryElementNode) ) EXIT
  END DO
  DO ParentElementNode=1,NParent
     IF ( Node == ParentElement % NodeIndexes(ParentElementNode) ) EXIT
  END DO
  !-------------------------
  ! Get Pressure Melting Point
  !-------------------------
  TempName =  GetString(ParentMaterial ,'Temperature Name', GotIt)
  PressureMeltingPoint(1:NParent) =&
       ListGetReal( ParentMaterial, TRIM(TempName) // ' Upper Limit',&
       NParent, ParentElement % NodeIndexes, GotIt)
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,a,a)') 'No entry for ', TRIM(TempName) // ' Upper Limit', ' found'
     CALL FATAL('iceproperties (basalMelting)',Message)
  END IF
  !---------------------------------------------------------------------
  ! get temperate slip coefficient and upper limit for slip coefficient
  !---------------------------------------------------------------------
  Model % CurrentElement => BoundaryElement
  BC => GetBC()

  IF (.NOT.ASSOCIATED(BC)) THEN
     CALL FATAL('iceproperties (basalSlip)','No Boundary Condition associated')
  ELSE
     TemperateSlipCoefficient(1:NBoundary) = GetReal(BC, 'Temperate Slip Coefficient', GotIt)
     IF (.NOT. GotIt) THEN
        CALL WARN('iceproperties (basalSlip)','Keyword >Temperate Slip Coefficient< not found')
        CALL WARN('iceproperties (basalSlip)','Asuming Default 5.63D08 [kg /(s m^2)]')
        TemperateSlipCoefficient(1:NBoundary) = 5.63D08
     END IF
     ThermalCoefficient = GetConstReal(BC, 'Thermal Coefficient', GotIt)
     IF (.NOT. GotIt) THEN
        CALL WARN('iceproperties (basalSlip)','Keyword >Thermal Coefficient< not found')
        CALL WARN('iceproperties (basalSlip)','Asuming Default 1 [1/K]')
        ThermalCoefficient =  1.0D00
     END IF
     SlipCoefficientMax = GetConstReal(BC, 'Max Slip Coefficient', GotIt)
     IF (.NOT. GotIt) THEN
        CALL WARN('iceproperties (basalSlip)','Keyword >Max Slip Coefficient< not found')
        CALL WARN('iceproperties (basalSlip)','Asuming Default 1.0E+16 [kg /(s m^2)]')
        SlipCoefficientMax = 1.0E+16
     END IF
  END IF
  !------------------------------
  ! check homologous temperature
  !------------------------------
  TempHom = MIN(Temperature - PressureMeltingPoint(ParentElementNode),0.0D00)
  !------------------------------
  ! get the result and check it
  !------------------------------
  basalSlipCoefficient = TemperateSlipCoefficient(BoundaryElementNode)*EXP(-1.0D00*TempHom*ThermalCoefficient)
  IF (basalSlipCoefficient < TemperateSlipCoefficient(BoundaryElementNode)) THEN
     WRITE(Message,'(A,e10.4)') 'Applied lower threshold of ', &
          TemperateSlipCoefficient(BoundaryElementNode)
     CALL INFO('iceproperties (basalSlip)',Message,Level=9)
     basalSlipCoefficient =  TemperateSlipCoefficient(BoundaryElementNode)
  ELSE IF(basalSlipCoefficient > SlipCoefficientMax) THEN
     WRITE(Message,'(A,e10.4)') 'Applied upper threshold of ', SlipCoefficientMax
     CALL INFO('iceproperties (basalSlip)',Message,Level=9)
     basalSlipCoefficient = SlipCoefficientMax
  END IF
  !------------------------------
  ! clean up
  !------------------------------
  DEALLOCATE(TemperateSlipCoefficient, PressureMeltingPoint)
END FUNCTION basalSlip

!*********************************************************************************************************************************
!*
!*  projecting vertical geothermal heat flux to boundary normal
!*
!*********************************************************************************************************************************
FUNCTION getNormalFlux( Model, Node, dummyArgument ) RESULT(NormalFlux)
!-----------------------------------------------------------
  USE DefUtils
!-----------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------
  !external variables
  TYPE(Model_t), TARGET :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: dummyArgument, NormalFlux
  !internal variables
  TYPE(Element_t), POINTER :: BoundaryElement
  TYPE(Nodes_t) :: Nodes
  TYPE(ValueList_t), POINTER :: BC
  INTEGER :: N, NBoundary, BoundaryElementNode, i, DIM, body_id, NumberOfBoundaryNodes, istat
  INTEGER, POINTER :: BoundaryReorder(:)
  REAL(KIND=dp) :: U, V, W, Normal(3), Gravity(3), direction(3),&
       HeatFlux, SqrtElementMetric
  REAL(KIND=dp), ALLOCATABLE ::  ExternalHeatFlux(:)
  REAL(KIND=dp), DIMENSION(:,:),  POINTER :: Work, BoundaryNormals,BoundaryTangent1, &
       BoundaryTangent2
  LOGICAL ::  FirstTime = .TRUE., GotIt, stat
!-----------------------------------------------------------
  SAVE FirstTime, NumberOfBoundaryNodes,BoundaryReorder,BoundaryNormals, &
       BoundaryTangent1, BoundaryTangent2, DIM

  IF (FirstTime) THEN 
     DIM = CoordinateSystemDimension()
     CALL CheckNormalTangentialBoundary( Model, &
          'Basal Melting', NumberOfBoundaryNodes, &
          BoundaryReorder, BoundaryNormals, BoundaryTangent1, &
          BoundaryTangent2, DIM )

     WRITE(Message,'(A,i6)') &
          'Number of boundary nodes on boundaries associated with melting:',&
          NumberOfBoundaryNodes
     CALL INFO('iceproperties (getNormalHeatFlux)','Message',Level=3)
     
     IF(NumberOfBoundaryNodes < 1) THEN
        CALL FATAL('iceproperties (getNormalHeatFlux)','No Points for basal melting found')
     END IF

     CALL AverageBoundaryNormals( Model, &
          'Basal Melting', NumberOfBoundaryNodes, &
          BoundaryReorder, BoundaryNormals, BoundaryTangent1, &
          BoundaryTangent2, DIM )
     FirstTime = .FALSE.
  END IF


  !--------------------------------
  ! Allocations
  !--------------------------------
  N = Model % MaxElementNodes 
  ALLOCATE(Nodes % x(N), Nodes % y(N), Nodes % z(N),&
       ExternalHeatFlux( N ),&
       STAT = istat)
  IF (istat /= 0) THEN
     CALL FATAL('iceproperties (etNormalHeatFlux)','Allocations failed')
  END IF

  !-----------------------------------------------------------------
  ! get some information upon active boundary element and its parent
  !-----------------------------------------------------------------
  BoundaryElement => Model % CurrentElement
  IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN
     CALL FATAL('iceproperties (normalFlux)','No boundary element found')
  END IF
  !-------------------------------------------
  ! Get normal of the boundary element at node
  !-------------------------------------------
  Nboundary = BoundaryElement % Type % NumberOfNodes
  DO BoundaryElementNode=1,Nboundary
     IF ( Node == BoundaryElement % NodeIndexes(BoundaryElementNode) ) EXIT
  END DO
  Normal(1:DIM) = BoundaryNormals(BoundaryReorder(Node),1:DIM)
  !-------------------------------
  ! get gravitational acceleration
  !-------------------------------
  Work => ListGetConstRealArray( Model % Constants,'Gravity',GotIt)
  IF ( GotIt ) THEN
     Gravity = Work(1:3,1)
  ELSE
     Gravity = 0.0D00
     CALL INFO('iceproperties (normalFlux)','No vector for Gravity (Constants) found', level=1)
     IF (DIM == 1) THEN
        Gravity(1) = -1.0D00
        CALL INFO('iceproperties (normalFlux)','setting direction to -1', level=1)
     ELSE IF (DIM == 2) THEN
        Gravity    =  0.00D0
        Gravity(2) = -1.0D00
        CALL INFO('iceproperties (normalFlux)','setting direction to (0,-1)', level=1)
     ELSE
        Gravity    =  0.00D00
        Gravity(3) = -1.0D00
        CALL INFO('iceproperties (normalFlux)','setting direction to (0,0,-1)', level=1)
     END IF
  END IF
  !------------------------
  ! get external heat flux
  !------------------------
  BC => GetBC()
  ExternalHeatFlux = 0.0D00  
  ExternalHeatFlux(1:NBoundary) = GetReal(BC, 'External Heat Flux', GotIt)
  IF (.NOT. GotIt) THEN
     CALL INFO('iceproperties (normalFlux)','No external heat flux given', Level=4)
  END IF
  !--------------------------------------------------------
  ! compute normal component of vertically aligned heatflux
  !--------------------------------------------------------
  NormalFlux = ExternalHeatFlux(BoundaryElementNode) * ABS(SUM(Gravity(1:DIM)*Normal(1:DIM)))
  !----------------------------------------------
  ! clean up before leaving
  !----------------------------------------------
  DEALLOCATE( Nodes % x,&
       Nodes % y,&
       Nodes % z,&
       ExternalHeatFlux)
END FUNCTION getNormalFlux

!*********************************************************************************************************************************
!*
!* heat conductivity of ice as a function of temperature (K):  k = c_1 * exp(c_2 * T[K]); c_2 < 0 
!*
!*********************************************************************************************************************************
FUNCTION getHeatConductivity( Model, N, temperature ) RESULT(conductivity)
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t) :: Model
  INTEGER :: N
  REAL(KIND=dp) :: temperature, conductivity
!------------ internal variables----------------------------
  TYPE(ValueList_t), POINTER :: Material
  INTEGER :: nMax,i,j,body_id,material_id,elementNodes,nodeInElement,istat
  REAL (KIND=dp), ALLOCATABLE :: conductivityExponentFactor(:), conductivityFactor(:)
  LOGICAL :: FirstTime = .TRUE., GotIt
!------------ remember this -------------------------------
  Save FirstTime, conductivityExponentFactor, conductivityFactor
  !-------------------------------------------
  ! Allocations 
  !------------------------------------------- 
  IF (FirstTime) THEN
     nMax = Model % MaxElementNodes
     ALLOCATE(conductivityExponentFactor(nMax),&
          conductivityFactor(nMax),&
          STAT=istat)
     IF ( istat /= 0 ) THEN
        CALL FATAL('iceproperties (getHeatConductivity)','Memory allocation error, Aborting.')
     END IF
     FirstTime = .FALSE.
     CALL INFO('iceproperties (getHeatConductivity)','Memory allocation done', level=3)
  END IF
  !-------------------------------------------
  ! get element properties
  !-------------------------------------------   
  IF ( .NOT. ASSOCIATED(Model % CurrentElement) ) THEN
     CALL FATAL('iceproperties (getHeatConductivity)', 'Model % CurrentElement not associated')
  END IF
  body_id = Model % CurrentElement % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  elementNodes = Model % CurrentElement % Type % NumberOfNodes
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2,a)') 'No material id for current element of node ',n,', body ',body_id,' found'
     CALL FATAL('iceproperties (getHeatConductivity)', Message)
  END IF
  DO nodeInElement=1,elementNodes
     IF ( N == Model % CurrentElement % NodeIndexes(nodeInElement) ) EXIT
  END DO
  Material => Model % Materials(material_id) % Values
  !-------------------------------------------
  ! get material properties
  !-------------------------------------------
  conductivityExponentFactor(1:elementNodes) = ListGetReal( Material,'Conductivity Exponent Factor', elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2,a)') 'No Conductivity Exponent Factor found in Material ', &
          material_id,' for node ', n, '.setting E=1'
     CALL FATAL('iceproperties (getHeatConductivity)', Message)
  END IF
  conductivityFactor(1:elementNodes) = ListGetReal( Material,'Conductivity Factor', elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2,a)') 'No Conductivity Factor found in Material ', material_id,' for node ', n, '.setting E=1'
     CALL FATAL('iceproperties (getHeatConductivity)', Message)
  END IF
  !-------------------------------------------
  ! compute heat conductivity
  !-------------------------------------------
  conductivity = conductivityFactor(nodeInElement)*EXP(conductivityExponentFactor(nodeInElement)*temperature)
END FUNCTION getHeatConductivity

!*********************************************************************************************************************************
!*
!* heat capacity of ice as a function of temperature (K):  k = c_1 + c_2 * T[C];
!*
!*********************************************************************************************************************************
FUNCTION getHeatCapacity( Model, N, temperature ) RESULT(capacity)
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t) :: Model
  INTEGER :: N
  REAL(KIND=dp) :: temperature, capacity
!------------ internal variables----------------------------
  REAL(KIND=dp) :: celsius
  INTEGER :: body_id, material_id
  TYPE(ValueList_t), POINTER :: Material
  LOGICAL :: UseCelsius = .FALSE., GotIt

  body_id = Model % CurrentElement % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  IF (.NOT.GotIt) CALL FATAL('iceproperties (getHeatCapacity)','No Material ID found')
  Material => Model % Materials(material_id) % Values
  UseCelsius = GetLogical( Material,'Use Celsius',GotIt )
  IF (.NOT.GotIt) UseCelsius = .FALSE.

  !-------------------------------------------
  ! compute celsius temperature and limit it 
  ! to 0 deg
  !------------------------------------------- 
  IF (UseCelsius) THEN
     celsius = MIN(temperature, 0.0d00)
  ELSE
     celsius = MIN(temperature - 2.7316D02,0.0d00)
  END IF
  !-------------------------------------------
  ! compute heat capacity
  !-------------------------------------------  
  capacity = 2.1275D03 + 7.253D00*celsius
END FUNCTION getHeatCapacity

!****************************************************************************************************************
!*
!* viscosity factor as a function of temperature
!*
!****************************************************************************************************************
FUNCTION getViscosityFactor( Model, n, temperature ) RESULT(visFact)
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: temperature, visFact
!------------ internal variables----------------------------
  TYPE(ValueList_t), POINTER :: Material
  INTEGER :: DIM,nMax,i,j,body_id,material_id,elementNodes,nodeInElement,istat
  REAL(KIND=dp) ::&
       rateFactor, aToMinusOneThird, gasconst, temphom
  REAL(KIND=dp), POINTER :: Hwrk(:,:,:)
  REAL (KIND=dp), ALLOCATABLE :: activationEnergy(:,:), arrheniusFactor(:,:),&
       enhancementFactor(:), viscosityExponent(:), PressureMeltingPoint(:),&
       Ux(:), Uy(:), Uz(:), LimitTemp(:)
  LOGICAL :: FirstTime = .TRUE., GotIt
  CHARACTER(LEN=MAX_NAME_LEN) :: TempName
!------------ remember this -------------------------------
  Save DIM, FirstTime, gasconst, activationEnergy, arrheniusFactor,&
       enhancementFactor, viscosityExponent, Hwrk, PressureMeltingPoint, &
       Ux, Uy, Uz, LimitTemp
!-----------------------------------------------------------
  !-----------------------------------------------------------
  ! Read in constants from SIF file and do some allocations
  !-----------------------------------------------------------
  IF (FirstTime) THEN
     ! inquire coordinate system dimensions  and degrees of freedom from NS-Solver
     ! ---------------------------------------------------------------------------
     DIM = CoordinateSystemDimension()
     ! inquire minimum temperature
     !------------------------- 
     gasconst = ListGetConstReal( Model % Constants,'Gas Constant',GotIt)
     IF (.NOT. GotIt) THEN
        gasconst = 8.314D00 ! m-k-s
        WRITE(Message,'(a,e10.4,a)') 'No entry for Gas Constant (Constants) in input file found. Setting to ',&
             gasconst,' (J/mol)'
        CALL INFO('iceproperties (getViscosityFactor)', Message, level=4)
     END IF
     nMax = Model % MaxElementNodes
     ALLOCATE(activationEnergy(2,nMax),&
          arrheniusFactor(2,nMax),&
          enhancementFactor(nMax),&
          LimitTemp( nMax),&
          PressureMeltingPoint( nMax ),&
          viscosityExponent(nMax),&
          Ux(nMax),&
          Uy(nMax),&
          Uz(nMax),&
          STAT=istat)
     IF ( istat /= 0 ) THEN
        CALL Fatal('iceproperties (getViscosityFactor)','Memory allocation error, Aborting.')
     END IF
     NULLIFY( Hwrk )
     FirstTime = .FALSE.
     CALL Info('iceproperties (getViscosityFactor)','Memory allocations done', Level=3)
  END IF
  !-------------------------------------------
  ! get element properties
  !-------------------------------------------   
  body_id = Model % CurrentElement % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  IF (.NOT.GotIt) CALL FATAL('iceproperties (getViscosityFactor)','No Material ID found')
  elementNodes = Model % CurrentElement % Type % NumberOfNodes
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2,a)') 'No material id for current element of node ',n,', body ',body_id,' found'
     CALL FATAL('iceproperties (getViscosityFactor)', Message)
  END IF
  DO nodeInElement=1,elementNodes
     IF ( N == Model % CurrentElement % NodeIndexes(nodeInElement) ) EXIT
  END DO
  Material => Model % Materials(material_id) % Values
  IF (.NOT.ASSOCIATED(Material)) THEN 
     WRITE(Message,'(a,I2,a,I2,a)') 'No Material for current element of node ',n,', body ',body_id,' found'
     CALL FATAL('iceproperties (getViscosityFactor)',Message)
  END IF
  !-------------------------------------------
  ! get material properties
  !-------------------------------------------
  ! activation energies
  !--------------------
  CALL ListGetRealArray( Material,'Activation Energies',Hwrk,elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2)') 'No Value for Activation Energy  found in Material ', material_id,' for node ', n
     CALL FATAL('iceproperties (getViscosityFactor)',Message)
  END IF
  IF ( SIZE(Hwrk,2) == 1 ) THEN
     DO i=1,MIN(3,SIZE(Hwrk,1))
        activationEnergy(i,1:elementNodes) = Hwrk(i,1,1:elementNodes)
     END DO
  ELSE
     WRITE(Message,'(a,I2,a,I2)') 'Incorrect array size for Activation Energy in Material ', material_id,' for node ', n
     CALL FATAL('iceproperties (getViscosityFactor)',Message)
  END IF
  ! Arrhenius Factors
  !------------------
  CALL ListGetRealArray( Material,'Arrhenius Factors',Hwrk,elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2)') 'No Value for Arrhenius Factors  found in Material ', material_id,' for node ', n
     CALL FATAL('iceproperties (getViscosityFactor)',Message)
  END IF
  IF ( SIZE(Hwrk,2) == 1 ) THEN
     DO i=1,MIN(3,SIZE(Hwrk,1))
        arrheniusFactor(i,1:elementNodes) = Hwrk(i,1,1:elementNodes)
     END DO
  ELSE
     WRITE(Message,'(a,I2,a,I2)') 'Incorrect array size for Arrhenius Factors in Material ', material_id,' for node ', n
     CALL FATAL('iceproperties (getViscosityFactor)',Message)
  END IF
  ! Enhancement Factor
  !-------------------
  enhancementFactor(1:elementNodes) = ListGetReal( Material,'Enhancement Factor', elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     enhancementFactor(1:elementNodes) = 1.0D00
     WRITE(Message,'(a,I2,a,I2,a)') 'No Enhancement Factor found in Material ', material_id,' for node ', n, '.setting E=1'
     CALL INFO('iceproperties (getViscosityFactor)', Message, level=9)
  END IF
  ! Threshold temperature for switching activation energies and Arrhenius factors
  !------------------------------------------------------------------------------
  LimitTemp(1:elementNodes) = ListGetReal( Material,'Limit Temperature', elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     LimitTemp(1:elementNodes) = -1.0D01
     WRITE(Message,'(a,I2,a,I2,a)') 'No keyword >Limit Temperature< found in Material ',&
          material_id,' for node ', n, '.setting to -10'
     CALL INFO('iceproperties (getFluidity)', Message, level=9)
  END IF
  ! Viscosity Exponent
  !-------------------
  viscosityExponent(1:elementNodes) = ListGetReal( Material,'Viscosity Exponent', elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     viscosityExponent(1:elementNodes) = 1.0D00/3.0D00
     WRITE(Message,'(a,I2,a,I2,a)') 'No Viscosity Exponent found in Material ', material_id,' for node ', n, '.setting k=1/3'
     CALL INFO('iceproperties (getViscosityFactor)', Message, level=9)
  END IF
  ! Pressure Melting Point and homologous temperature
  !--------------------------------------------------
  TempName =  GetString(Material  ,'Temperature Name', GotIt)
  IF (.NOT.GotIt) CALL FATAL('iceproperties (getViscosityFactor)','No Temperature Name found')
  PressureMeltingPoint(1:elementNodes) =&
       ListGetReal( Material, TRIM(TempName) // ' Upper Limit',&
       elementNodes, Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT.GotIt) THEN
     temphom = 0.0d00
     WRITE(Message,'(A,A,A,i3,A)') 'No entry for ',TRIM(TempName) // ' Upper Limit',&
          ' found in material no. ', material_id,'. Using 273.16 K.'
     CALL WARN('iceproperties (getViscosityFactor)',Message)
  ELSE
     temphom = MIN(temperature - PressureMeltingPoint(nodeInElement), 0.0d00)
  END IF
  !-------------------------------------------
  ! homologous Temperature is below 10 degrees
  !-------------------------------------------
  IF (temphom < LimitTemp(nodeInElement)) THEN
     i=1
     !-------------------------------------------
     ! homologous Temperature is above 10 degrees
     !-------------------------------------------
  ELSE
     i=2
  END IF
  rateFactor =&
       arrheniusFactor(i,nodeInElement)*exp(-1.0D00*activationEnergy(i,nodeInElement)/(gasconst*(2.7316D02 + temphom)))
  visFact = 0.5D00*(enhancementFactor(nodeInElement)*rateFactor)**(-1.0e00*viscosityExponent(nodeInElement))
END FUNCTION getViscosityFactor

!****************************************************************************************************************
!*
!* fluidity as a function of homologous temperature
!*
!****************************************************************************************************************
FUNCTION getFluidity( Model, n, temperature ) RESULT(fluidity)
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: temperature, fluidity
!------------ internal variables----------------------------
  TYPE(ValueList_t), POINTER :: Material
  INTEGER :: DIM,nMax,i,j,body_id,material_id,elementNodes,nodeInElement,istat
  REAL(KIND=dp) ::&
       rateFactor, aToMinusOneThird, gasconst, temphom
  REAL(KIND=dp), POINTER :: Hwrk(:,:,:)
  REAL (KIND=dp), ALLOCATABLE :: activationEnergy(:,:), arrheniusFactor(:,:),&
       enhancementFactor(:), viscosityExponent(:), PressureMeltingPoint(:),&
       LimitTemp(:)
  LOGICAL :: FirstTime = .TRUE., GotIt
  CHARACTER(LEN=MAX_NAME_LEN) :: TempName
!------------ remember this -------------------------------
  Save DIM, FirstTime, gasconst, activationEnergy, arrheniusFactor,&
       enhancementFactor,  Hwrk, PressureMeltingPoint, LimitTemp
!-----------------------------------------------------------
  !-----------------------------------------------------------
  ! Read in constants from SIF file and do some allocations
  !-----------------------------------------------------------
  IF (FirstTime) THEN
     ! inquire coordinate system dimensions  and degrees of freedom from NS-Solver
     ! ---------------------------------------------------------------------------
     DIM = CoordinateSystemDimension()
     ! inquire minimum temperature
     !------------------------- 
     gasconst = ListGetConstReal( Model % Constants,'Gas Constant',GotIt)
     IF (.NOT. GotIt) THEN
        gasconst = 8.314D00 ! m-k-s
        WRITE(Message,'(a,e10.4,a)') 'No entry for Gas Constant (Constants) in input file found. Setting to ',&
             gasconst,' (J/mol)'
        CALL INFO('iceproperties (getFluidity)', Message, level=4)
     END IF
     nMax = Model % MaxElementNodes
     ALLOCATE(activationEnergy(2,nMax),&
          arrheniusFactor(2,nMax),&
          enhancementFactor(nMax),&
          PressureMeltingPoint( nMax ),&
          LimitTemp( nMax ), &
          STAT=istat)
     IF ( istat /= 0 ) THEN
        CALL Fatal('iceproperties (getFluidity)','Memory allocation error, Aborting.')
     END IF
     NULLIFY( Hwrk )
     FirstTime = .FALSE.
     CALL Info('iceproperties (getFluidity)','Memory allocations done', Level=3)
  END IF
  !-------------------------------------------
  ! get element properties
  !-------------------------------------------   
  body_id = Model % CurrentElement % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  IF (.NOT.GotIt) CALL FATAL('iceproperties (getFluidity)','No Material ID found')
  elementNodes = Model % CurrentElement % Type % NumberOfNodes
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2,a)') 'No material id for current element of node ',n,', body ',body_id,' found'
     CALL FATAL('iceproperties (getFluidity)', Message)
  END IF
  DO nodeInElement=1,elementNodes
     IF ( N == Model % CurrentElement % NodeIndexes(nodeInElement) ) EXIT
  END DO
  Material => Model % Materials(material_id) % Values
  IF (.NOT.ASSOCIATED(Material)) THEN 
     WRITE(Message,'(a,I2,a,I2,a)') 'No Material for current element of node ',n,', body ',body_id,' found'
     CALL FATAL('iceproperties (getFluidity)',Message)
  END IF
  !-------------------------------------------
  ! get material properties
  !-------------------------------------------
  ! activation energies
  !--------------------
  CALL ListGetRealArray( Material,'Activation Energies',Hwrk,elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2)') 'No Value for Activation Energy  found in Material ', material_id,' for node ', n
     CALL FATAL('iceproperties (getFluidity)',Message)
  END IF
  IF ( SIZE(Hwrk,2) == 1 ) THEN
     DO i=1,MIN(3,SIZE(Hwrk,1))
        activationEnergy(i,1:elementNodes) = Hwrk(i,1,1:elementNodes)
     END DO
  ELSE
     WRITE(Message,'(a,I2,a,I2)') 'Incorrect array size for Activation Energy in Material ', material_id,' for node ', n
     CALL FATAL('iceproperties (getFluidity)',Message)
  END IF
  ! Arrhenius Factors
  !------------------
  CALL ListGetRealArray( Material,'Arrhenius Factors',Hwrk,elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2)') 'No Value for Arrhenius Factors  found in Material ', material_id,' for node ', n
     CALL FATAL('iceproperties (getFluidity)',Message)
  END IF
  IF ( SIZE(Hwrk,2) == 1 ) THEN
     DO i=1,MIN(3,SIZE(Hwrk,1))
        arrheniusFactor(i,1:elementNodes) = Hwrk(i,1,1:elementNodes)
     END DO
  ELSE
     WRITE(Message,'(a,I2,a,I2)') 'Incorrect array size for Arrhenius Factors in Material ', material_id,' for node ', n
     CALL FATAL('iceproperties (getFluidity)',Message)
  END IF
  ! Enhancement Factor
  !-------------------
  enhancementFactor(1:elementNodes) = ListGetReal( Material,'Enhancement Factor', elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     enhancementFactor(1:elementNodes) = 1.0D00
     WRITE(Message,'(a,I2,a,I2,a)') 'No Enhancement Factor found in Material ', material_id,' for node ', n, '.setting E=1'
     CALL INFO('iceproperties (getFluidity)', Message, level=9)
  END IF
  ! Threshold temperature for switching activation energies and Arrhenius factors
  !------------------------------------------------------------------------------
  LimitTemp(1:elementNodes) = ListGetReal( Material,'Limit Temperature', elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     LimitTemp(1:elementNodes) = -1.0D01
     WRITE(Message,'(a,I2,a,I2,a)') 'No keyword >Limit Temperature< found in Material ',&
          material_id,' for node ', n, '.setting to -10'
     CALL INFO('iceproperties (getFluidity)', Message, level=9)
  END IF
  ! Pressure Melting Point and homologous temperature
  !--------------------------------------------------
  TempName =  GetString(Material,'Temperature Name', GotIt)
  IF (.NOT.GotIt) CALL FATAL('iceproperties (getFluidity)','No Temperature Name found')
  PressureMeltingPoint(1:elementNodes) =&
       ListGetReal( Material, TRIM(TempName) // ' Upper Limit',&
       elementNodes, Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT.GotIt) THEN
     temphom = 0.0d00
     WRITE(Message,'(A,A,A,i3,A)') 'No entry for ',TRIM(TempName) // ' Upper Limit',&
          ' found in material no. ', material_id,'. Using 273.16 K.'
     CALL WARN('iceproperties (getFluidity)',Message)
  ELSE
     temphom = MIN(temperature - PressureMeltingPoint(nodeInElement), 0.0d00)
  END IF
  !-----------------------------------------------------
  ! homologous Temperature is below temperature treshold
  !----------------------------------------------------
  IF (temphom < LimitTemp(nodeInElement))THEN
     i=1
     !-----------------------------------------------------
     ! homologous Temperature is above temperature treshold
     !-----------------------------------------------------
  ELSE
     i=2
  END IF
  rateFactor =&
       arrheniusFactor(i,nodeInElement)*exp(-1.0D00*activationEnergy(i,nodeInElement)/(gasconst*(2.7316D02 + temphom)))
  fluidity = 2.0D00*enhancementFactor(nodeInElement)*rateFactor
END FUNCTION getFluidity


!*********************************************************************************************************************************
RECURSIVE SUBROUTINE getTotalViscosity(Model,Solver,Timestep,TransientSimulation)
  USE DefUtils
  USE Materialmodels
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL :: TransientSimulation
  REAL(KIND=dp) :: Timestep
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(GaussIntegrationPoints_t) :: IP
  TYPE(ValueList_t), POINTER :: Material, SolverParams
  TYPE(Variable_t), POINTER :: TotViscSol
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t) :: Nodes
  INTEGER :: DIM,nMax,i,p,q,t,N,body_id,material_id,nodeInElement,istat
  INTEGER, POINTER :: TotViscPerm(:)
  REAL(KIND=dp) ::  LocalDensity, LocalViscosity,effVisc,detJ,Norm
  REAL(KIND=dp), POINTER :: Hwrk(:,:,:)
  REAL(KIND=dp), POINTER :: TotVisc(:)
  REAL (KIND=dp), ALLOCATABLE ::  Ux(:), Uy(:), Uz(:), Density(:), Viscosity(:),&
       STIFF(:,:),FORCE(:),Basis(:),dBasisdx(:,:),ddBasisddx(:,:,:)
  LOGICAL :: FirstTime = .TRUE., GotIt, stat, limitation
!------------ remember this -------------------------------
  Save DIM, FirstTime, Nodes, Density, Ux, Uy, Uz, Viscosity, STIFF, FORCE, Basis,dBasisdx,ddBasisddx
  

  IF (FirstTime) THEN
     ! inquire coordinate system dimensions  and degrees of freedom from NS-Solver
     ! ---------------------------------------------------------------------------
     DIM = CoordinateSystemDimension()
     nMax = Model % MaxElementNodes 
     ALLOCATE(Ux(nMax),&
          Uy(nMax),&
          Uz(nMax),&
          Density(nMax),&
          Viscosity(nMax),&
          FORCE(nMax), &
          STIFF(nMax,nMax),&
          Nodes % x(nMax), &
          Nodes % y(nMax), &
          Nodes % z(nMax), &
          Basis(nMax),&
          dBasisdx(nMax,3),&
          ddBasisddx(nMax,3,3),&
          STAT=istat)
     IF ( istat /= 0 ) THEN
        CALL Fatal('iceproperties (getTotalViscosity)','Memory allocation error, Aborting.')
     END IF
     NULLIFY( Hwrk )
     FirstTime = .FALSE.
     CALL Info('iceproperties (getTotalViscosity)','Memory allocations done', Level=3)
  END IF

  CALL DefaultInitialize()

  DO t=1,Solver % NumberOfActiveElements 
     !-------------------------------------------
     ! get element properties
     !------------------------------------------- 
     Element => GetActiveElement(t)
     N = GetElementNOFNodes()
     body_id = Element % BodyId
     material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
     IF (.NOT. GotIt) THEN
        WRITE(Message,'(a,I2,a,I2,a)') &
             'No material id for current element of node ',n,', body ',body_id,' found'
        CALL FATAL('iceproperties (getTotalViscosity)', Message)
     
     END IF

     !-----------------------------------------
     ! get Material properties
     !-----------------------------------------
     Material => GetMaterial()
     IF (.NOT.ASSOCIATED(Material)) THEN 
        WRITE(Message,'(a,I2,a,I2,a)') &
             'No Material for current element of node ',n,', body ',body_id,' found'
        CALL FATAL('iceproperties (getTotalViscosit)',Message)
     END IF
     Density = ListGetReal( Material,'Density', N, Element % NodeIndexes, GotIt )
     IF (.NOT. GotIt) THEN
        WRITE(Message,'(a,I2,a,I2)') 'No Value for Density found in Material ', material_id,' for node ', n
        CALL FATAL('iceproperties (getTotalViscosity)',Message)
     END IF
     Viscosity = ListGetReal( Material,'Viscosity', N, Element % NodeIndexes, GotIt )
     IF (.NOT. GotIt) THEN
        WRITE(Message,'(a,I2,a,I2)') 'No Value for Viscosity found in Material ', material_id,' for node ', n
        CALL FATAL('iceproperties (getTotalViscosity)',Message)
     END IF

     !-----------------------------------------
     ! Velocity Field
     !----------------------------------------
     Ux = 0.0d00
     Uy = 0.0d00
     Uz = 0.0d00
  
     CALL GetScalarLocalSolution( Ux, 'Velocity 1')
     IF (DIM>1) THEN
        CALL GetScalarLocalSolution( Uy, 'Velocity 2')
        IF (DIM == 3) CALL GetScalarLocalSolution( Uz, 'Velocity 3')
     END IF
     

     STIFF = 0.0d00
     FORCE = 0.0d00

     IP = GaussPoints( Element )
     CALL GetElementNodes( Nodes )

     DO i=1,IP % n
        stat = ElementInfo( Element, Nodes, IP % U(i), IP % V(i), &
             IP % W(i),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )
        !
        ! get local Material parameters at Gauss-point
        !
        LocalDensity = SUM(Density(1:n)*Basis(1:n))
        LocalViscosity = SUM(Viscosity(1:n)*Basis(1:n))
        !
        ! get effective Viscosity at Integration point
        !
        effVisc = EffectiveViscosity( LocalViscosity, LocalDensity,&
             Ux, Uy, Uz, &
             Element, Nodes, N, N,&
             IP % U(i), IP % V(i), IP % W(i))
         IF (effVisc .le. 0.0E00) THEN
           WRITE(Message,'(A,i10,A,i10,A,e13.3)')&
                'effective viscosity for Gauss point no. ', i, ' in element no. ', t,' is negative:', effVisc
           CALL WARN('iceproperties (getTotalViscosity)',Message)
        END IF
        DO p=1,n
           FORCE(p) = FORCE(p) + IP % S(i) * DetJ * effVisc * Basis(p)
           DO q=1,n
              STIFF(p,q) = STIFF(p,q) + IP % S(i) * detJ * Basis(q)*Basis(p)
           END DO
        END DO
     END DO
     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
  CALL DefaultFinishAssembly()
!   CALL DefaultDirichletBCs()
  Norm = DefaultSolve()
  SolverParams => GetSolverParams()
  limitation = GetLogical( SolverParams,'Positive Values',GotIt )
  IF (.NOT. GotIt) limitation = .FALSE.
  IF (limitation) THEN
     CALL INFO('iceproperties (getTotalViscosity)','Results limited to positive values',Level=1)
     TotViscSol => Solver % Variable
     TotViscPerm  => TotViscSol % Perm
     TotVisc => TotViscSol % Values
     DO i= 1,Solver % Mesh % NumberOfNodes
        TotVisc(i) = MAX(TotVisc(i),0.0D00)
     END DO
  END IF
END SUBROUTINE getTotalViscosity



!*********************************************************************************************************************************
!*
!* viscosity as a function of the viscosity factor
!*
!*********************************************************************************************************************************
FUNCTION getCriticalShearRate( Model, n, temperature ) RESULT(critShear)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
!-----------------------------------------------------------
   IMPLICIT NONE
!------------ external variables ---------------------------
   TYPE(Model_t) :: Model
   INTEGER :: n
   REAL(KIND=dp) :: critShear, temperature
!------------ internal variables----------------------------
   TYPE(Element_t), POINTER :: Element
   TYPE(ValueList_t),POINTER :: Material
   INTEGER :: DIM, body_id, material_id
   REAL(KIND=dp) ::&
        visFact, cuttofViscosity, power, rateFactor, aToMinusOneThird,&
        activationEnergy, gasconst, enhancementFactor
   LOGICAL :: GotIt, FirstTime = .TRUE. 
!------------ remember this -------------------------------
   Save DIM, FirstTime
!-----------------------------------------------------------
! Read in constants from SIF file
!-----------------------------------------------------------
   IF (FirstTime) THEN
      ! inquire coordinate system dimensions  and degrees of freedom from NS-Solver
      ! ---------------------------------------------------------------------------
      DIM = CoordinateSystemDimension()
      gasconst = ListGetConstReal( Model % Constants,'Gas Constant',GotIt)
      IF (.NOT. GotIt) THEN
         gasconst = 8.314 ! m-k-s
         WRITE(Message,'(a,e10.4,a)') 'No entry for Gas Constant (Constants) in input file found. Setting to ',&
              gasconst,' (J/mol)'
         CALL INFO('iceproperties (getCriticalShearRate)', Message, level=9)
      END IF
      FirstTime = .FALSE.
   END IF
   ! inquire Model parameters
   !-------------------------
   Element => Model % CurrentElement
   body_id = Element % BodyId
   material_id = ListGetInteger( Model % Bodies(body_id) % Values,&
        'Material', Gotit, minv=1, maxv=Model % NumberOFMaterials)
   Material => Model % Materials(material_id) % Values
   power=  ListGetConstReal( Material, 'Viscosity Exponent', Gotit)
   IF (.NOT.Gotit) THEN
      CALL FATAL('iceproperties (getCriticalShearRate)', 'Viscosity Exponent not found')
   ELSE
      WRITE(Message,'(a,e10.4)') 'Viscosity Exponent =', power
      CALL INFO('iceproperties (getCriticalShearRate)', Message, level=9)
   END IF
   cuttofViscosity =  ListGetConstReal(Material, 'Cutoff Viscosity', Gotit)
   IF (.NOT.Gotit) THEN
      CALL FATAL('iceproperties (getCriticalShearRate)', 'Cutoff Viscosity not found')
   ELSE
      WRITE(Message,'(a,e10.4)') 'Cutoff Viscosity =', cuttofViscosity
      CALL INFO('iceproperties (getCriticalShearRate)', Message, level=9) 
   END IF

   ! get viscosity factor for local node
   ! -----------------------------------
!   visFact = ListGetReal(Material, 'Viscosity Factor',1, n, GotIt)
!   IF (.NOT.Gotit) THEN
!      WRITE(Message,'(a,i4,a)') 'Viscosity Factor for point no. ', n, ' not found' 
!      CALL FATAL('iceproperties (getCriticalShearRate)', 'Cutoff Viscosity not found'Message)
!   ELSE
!      WRITE(Message,'(a,i4,a,e10.4)') 'Viscosity Factor for point no. ',&
!           n, ' = ', visFact
!      CALL INFO('iceproperties (getCriticalShearRate)', Message, level=4) 
!   END IF

   IF (temperature < 263.15) THEN
      activationEnergy = 6.0e04 ! m-k-s
      aToMinusOneThird = 4.42577e16
      enhancementFactor = 1.0e00
   ELSE 
      activationEnergy = 1.39e05
      aToMinusOneThird = 2.62508e09
      enhancementFactor = 1.0e00
   END IF
   rateFactor = exp(-activationEnergy/(gasconst*temperature))
   visFact = 0.5*aToMinusOneThird*(enhancementFactor*rateFactor)**(-1.0e00/3.00e00)
   IF (visFact .NE. 0.0e0) THEN
      critShear = (cuttofViscosity/visFact)**(1.0e0/(1.0e0 - power))
   ELSE
      critShear = 0.0d0
   END IF
END FUNCTION getCriticalShearRate


!*********************************************************************************************************************************
!*
!* heat conductivity factor 
!*
!*********************************************************************************************************************************
FUNCTION getHeatConductivityFactorNew( Model, Node, relativeDensity ) RESULT(heatCondFact)
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: relativeDensity, heatCondFact
  REAL(KIND=dp) :: volumefraction
  heatCondFact = 9.828D00 * (7.1306D-02 - 4.7908D-01 * relativeDensity + 1.40778D00 * relativeDensity * relativeDensity)
END FUNCTION getHeatConductivityFactorNew


!*********************************************************************************************************************************
!*
!* pressure melting point
!*
!*********************************************************************************************************************************
FUNCTION getPressureMeltingPoint( Model, n, Pressure ) RESULT(PressureMeltingPoint)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
!-----------------------------------------------------------
   IMPLICIT NONE
!------------ external variables ---------------------------
   TYPE(Model_t) :: Model
   INTEGER :: n
   REAL(KIND=dp) :: PressureMeltingPoint, Pressure
!------------ internal variables----------------------------
   TYPE(Element_t), POINTER :: Element
   TYPE(ValueList_t),POINTER :: Material
   INTEGER :: body_id, material_id, nodeInElement, istat, elementNodes
   REAL(KIND=dp), ALLOCATABLE :: ClausiusClapeyron(:)
   LOGICAL :: GotIt, FirstTime = .TRUE., UseCelsius = .FALSE.
!------------ remember this -------------------------------
   Save FirstTime, ClausiusClapeyron
  !-------------------------------------------
  ! Allocations 
  !------------------------------------------- 
  IF (FirstTime) THEN
     ALLOCATE(ClausiusClapeyron(Model % MaxElementNodes),&
          STAT=istat)
     IF ( istat /= 0 ) THEN
        CALL FATAL('iceproperties (getPressureMeltingPoint)','Memory allocation error, Aborting.')
     END IF
     FirstTime = .FALSE.
     CALL INFO('iceproperties (getPressureMeltingPoint)','Memory allocation done', level=3)
  END IF
  !-------------------------------------------
  ! get element properties
  !-------------------------------------------   
  IF ( .NOT. ASSOCIATED(Model % CurrentElement) ) THEN
     CALL FATAL('iceproperties (getPressureMeltingPoint)', 'Model % CurrentElement not associated')
  ELSE
     Element => Model % CurrentElement
  END IF
  body_id = Element % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  elementNodes = Element % Type % NumberOfNodes
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2,a)') 'No material id for current element of node ',elementNodes,', body ',body_id,' found'
     CALL FATAL('iceproperties (getPressureMeltingPoint)', Message)
  END IF
  DO nodeInElement=1,elementNodes
     IF ( N == Model % CurrentElement % NodeIndexes(nodeInElement) ) EXIT
  END DO
  Material => Model % Materials(material_id) % Values
  !-------------------------
  ! get material parameters
  !-------------------------
  ClausiusClapeyron(1:elementNodes) = &
       ListGetReal( Material, 'Clausius Clapeyron', elementNodes, Element % NodeIndexes, GotIt)
  IF (.NOT. GotIt) THEN
     CALL FATAL('iceproperties (getPressureMeltingPoint)','No value for Clausius Clapeyron parameter found')
  END IF
  UseCelsius = GetLogical(Material,'Use Celsius',GotIt )
  IF (.NOT.GotIt) UseCelsius = .FALSE.
  !-------------------------------
  ! compute pressure melting point
  !-------------------------------
  PressureMeltingPoint = 2.7316D02 - ClausiusClapeyron(nodeInElement)*MAX(Pressure,0.0d00)
END FUNCTION getPressureMeltingPoint

!
!
!  Function a(D) and b(D) from Gagliardini and Meyssonier, 1997.
!  modified to fulfill the condition 3x a(D) >= 2x b(D) for D > 0.1 
!  for n=3 only 

FUNCTION OldParameterA ( Model, nodenumber, D ) RESULT(a)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber
   REAL(KIND=dp) :: a, D, DD     

    IF (D >= 1.0_dp) THEN
      a = 1.0_dp           

    ELSE IF (D <= 0.81) THEN
      DD = D
      IF (D < 0.1 ) DD = 0.1_dp
      a = EXP( 16.02784497_dp - 19.25_dp * DD )

    ELSE
      a =  1.0_dp  + 2.0/3.0 * (1.0_dp - D) 
      a = a / ( D**1.5 )

    End If
  END FUNCTION OldParameterA


FUNCTION OldParameterB ( Model, nodenumber, D ) RESULT(b)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber
   REAL(KIND=dp) :: b, D, DD 

    IF (D >= 1.0_dp) THEN
      b = 0.0_dp           

    ELSE IF (D <= 0.81) THEN
      DD = D
      IF (D < 0.1 ) DD = 0.1_dp
      b = EXP( 16.75024378_dp - 22.51_dp * DD )

    ELSE
      b = (1.0_dp - D)**(1.0/3.0) 
      b = 3.0/4.0 * ( b / (3.0 * (1 - b)) )**1.5

    End If
END FUNCTION OldParameterB

!
!  Function a(D) and b(D) from Gagliardini and Meyssonier, 1997.
!  modified to fulfill the condition 3x a(D) >= 2x b(D) for D > 0.4 
!  Such that a(0.4) = b(0.4) = 10^3
!  for n=3 only 

FUNCTION ParameterA ( Model, nodenumber, D ) RESULT(a)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber
   REAL(KIND=dp) :: a, D, DD     

    IF (D >= 1.0_dp) THEN
      a = 1.0_dp           

    ELSE IF (D <= 0.81) THEN
      DD = D
      IF (D < 0.4 ) DD = 0.4_dp
      a = EXP( 13.2224_dp - 15.78652_dp * DD )

    ELSE
      a =  1.0_dp  + 2.0/3.0 * (1.0_dp - D) 
      a = a / ( D**1.5 )

    End If
!    write(*,*)DD,a
END FUNCTION ParameterA 


FUNCTION ParameterB ( Model, nodenumber, D ) RESULT(b)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber
   REAL(KIND=dp) :: b, D, DD 

    IF (D >= 1.0_dp) THEN
      b = 0.0_dp           

    ELSE IF (D <= 0.81) THEN
      DD = D
      IF (D < 0.4 ) DD = 0.4_dp
      b = EXP( 15.09371_dp - 20.46489_dp * DD )

    ELSE
      b = (1.0_dp - D)**(1.0/3.0) 
      b = 3.0/4.0 * ( b / (3.0 * (1 - b)) )**1.5

    End If
!    write(*,*)DD,b
END FUNCTION ParameterB 


! *****************************************************************************
! calculates the SIA stress components on the free boundaries
!
! *****************************************************************************
RECURSIVE SUBROUTINE getSIAstress( Model,Solver,Timestep,TransientSimulation)
  USE DefUtils
  USE Materialmodels
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL :: TransientSimulation
  REAL(KIND=dp) :: Timestep
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: PointerToSolver
  TYPE(Nodes_t) :: Nodes
  TYPE(Element_t),POINTER :: Element
  TYPE(Variable_t), POINTER :: SIASol, DepthSol, SurfSol1, SurfSol2
  TYPE(ValueList_t), POINTER :: Equation,Material,SolverParams,BodyForce,BC,Constants
  CHARACTER(LEN=MAX_NAME_LEN) :: TempName, SolverName
  INTEGER :: i, j, k, l, t, N, BN,  M, DIM, SIADOFs, istat, material_id
  INTEGER, POINTER :: DepthSolPerm(:), SurfSol1Perm(:), SurfSol2Perm(:),&
       SIAPerm(:),BoundaryReorder(:)
  REAL(KIND=dp) ::  U, V, W, SqrtElementMetric, SIAstress(2), HydrostaticPressure
  REAL(KIND=dp), POINTER :: SIA(:), Depth(:), Surf1(:), Surf2(:), &
       BoundaryNormals(:,:), BoundaryTangent1(:,:), BoundaryTangent2(:,:)
  LOGICAL :: GotIt, FirstTimeAround=.TRUE., stat

  SAVE FirstTimeAround, SolverName, DIM, &
       BoundaryNormals, BoundaryTangent1, BoundaryTangent2, BoundaryReorder, BN

  ! assign solver name for communicative output
  !-----------------------------------------------------------------------
  WRITE(SolverName, '(A)') 'iceproperties (getSIAstress))'
  !-----------------------------------------------------------------------
  ! get solver variable
  !-----------------------------------------------------------------------
  SIASol => Solver % Variable
  IF (.NOT.ASSOCIATED(SIASol)) THEN
     CALL FATAL(SolverName,'No variable associated')
  END IF
  SIAPerm  => SIASol % Perm
  SIA => SIASol % Values
  SIADOFs =  SIASol % DOFs
  SIA = 0.0d00
  !-----------------------------------------------------------------------
  ! 
  !-----------------------------------------------------------------------
  IF ( FirstTimeAround .OR. Solver % Mesh % Changed ) THEN
     DIM = CoordinateSystemDimension()
     N = Solver % Mesh % MaxElementNodes
     M = Model % Mesh % NumberOfNodes        
     IF ( .NOT.FirstTimeAround ) &
          DEALLOCATE( &
          BoundaryReorder, &
          BoundaryNormals, &
          BoundaryTangent1, &
          BoundaryTangent2)
     ! check boundaries for calculation of SIA-stress and allocate necessary 
     ! space for averaged Normal and Tangentials
     !-----------------------------------------------------------------------
     CALL CheckNormalTangentialBoundary( Model, &
          'Calc SIA', BN, &
          BoundaryReorder, BoundaryNormals, BoundaryTangent1, &
          BoundaryTangent2, DIM )
     WRITE(Message,'(A,i6)') &
          'Number of boundary nodes on boundaries associated with SIA stresses:',&
          BN
     CALL INFO(SolverName,Message,Level=3)
     ! compute averaged normals and tangentials for  boundaries designated for
     ! SIA-stress boundary elements
     !-----------------------------------------------------------------------
     CALL AverageBoundaryNormals(Model, &
          'Calc SIA', BN, &
          BoundaryReorder, BoundaryNormals, BoundaryTangent1, &
          BoundaryTangent2, DIM )

     ! read in variables for flow depth and for free surface gradients
     !-----------------------------------------------------------------------
     DepthSol => VariableGet( Solver % Mesh % Variables, "Depth" )
     IF ( ASSOCIATED( DepthSol ) ) THEN
        DepthSolPerm => DepthSol % Perm
        Depth => DepthSol % Values
     ELSE
        WRITE(Message, '(A)') 'Could not find surface Gradient 1 field variable '
        CALL FATAL(SolverName,Message)
     END IF

     SurfSol1 => VariableGet( Solver % Mesh % Variables, "FreeSurfGrad1" )
     IF ( ASSOCIATED( SurfSol1 ) ) THEN
        SurfSol1Perm => SurfSol1 % Perm
        Surf1 => SurfSol1 % Values
     ELSE
        WRITE(Message, '(A)') 'Could not find Surface Gradient 1 field variable '
        CALL FATAL(SolverName,Message)
     END IF
     IF (DIM > 2) THEN
        SurfSol2 => VariableGet( Solver % Mesh % Variables, "FreeSurfGrad2" )
        IF ( ASSOCIATED( SurfSol2 ) ) THEN
           SurfSol2Perm => SurfSol2 % Perm
           Surf2 => SurfSol2 % Values
        ELSE
           WRITE(Message, '(A)') 'Could not find Surface Gradient 2 field variable '
           CALL FATAL(SolverName,Message)
        END IF
     END IF
     !-----------------------------------------------------------------------
     ! loop over all nodes in mesh
     !-----------------------------------------------------------------------
     DO i=1,Solver % Mesh % NumberOfNodes
        j = BoundaryReorder(i) ! projection from real space to SIA-stress 
        ! boundary space
        k = SIAPerm(i) ! projection from real space to solver matrix coordinate

        ! if boundary element node with SIA-stress condition enabled
        !-----------------------------------------------------------------------
        IF (j > 0) THEN
           HydrostaticPressure = -9.18D02 * 0.981E01 * Depth(DepthSolPerm(i))        
           SIAstress(1) =  HydrostaticPressure * Surf1(SurfSol1Perm(i))
           IF (DIM > 2) THEN
              SIAstress(2) = HydrostaticPressure * Surf2(SurfSol2Perm(i))
           ELSE
              SIAstress(2) = 0.0D00
           END IF

           ! vector product between Cauchy-stress tensor and surface normal
           !  = stress vector
           !-----------------------------------------------------------------------
           !Compute first element of the stress vector P*n(1)+txy*n(2) For  DIM=2 (two dimensions)
           !Compute first element of the stress vector P*n(1)+txz*n(3) for  Dim=3 (three-dimensions)
           SIA(SIADOFs*(k-1)+1)= HydrostaticPressure * BoundaryNormals(BoundaryReorder(i),1) &
                + SIAstress(1) * BoundaryNormals(BoundaryReorder(i),DIM)
           !Compute the second element of the stress vector in three-dimension P*n(2)+tyz*n(3)
           SIA(SIADOFs*(k-1)+2)= HydrostaticPressure * BoundaryNormals(BoundaryReorder(i),2) &
                + SIAstress(2) * BoundaryNormals(BoundaryReorder(i),3) ! this line doesn't contribute if 2d
           !If DIM=2, compute the second element of the stress vector in two-dimensions P*n(2)+txy*n(1)
           IF (DIM == 2) &
                SIA(SIADOFs*(k-1)+2)= SIA(SIADOFs*(k-1)+2) &
                + SIAstress(1) * BoundaryNormals(BoundaryReorder(i),1)
           !If DIM=3, compute the third element of the stress vector in three-dimensions P*n(3)+txz*n(1)+tyz*n(2)
           IF (DIM > 2) &
                SIA(SIADOFs*(k-1)+3)= HydrostaticPressure * BoundaryNormals(BoundaryReorder(i),3) &
                + SIAstress(1) * BoundaryNormals(BoundaryReorder(i),1) &
                + SIAstress(2) * BoundaryNormals(BoundaryReorder(i),2)   
        END IF
     END DO
     FirstTimeAround = .FALSE.
  END IF
END SUBROUTINE getSIAstress

!*********************************************************************************************************************************
!* Computation of component of the SIA stress vector. This code is  used for a rectangular type geometry where the automatic solution "getSIAstres" does not provide good results at the corner. Here we compute the components Pressure i in the sif file independantly for each side, North, Est, South and West. 
!*
!*********************************************************************************************************************************

!----------------------------------------------------------------------------------------------
!StressSIA: Compute the corresponding shallow-ice shear stress
!----------------------------------------------------------------------------------------------
FUNCTION StressSIA( Model, Node, surfaceGrad) RESULT(shearstress)
  USE Types
  USE CoordinateSystems
  USE SolverUtils
  USe ElementDescription
!-------------------------------------------------------------------
IMPLICIT NONE
!-------------------------external variables-----------------------
TYPE(Model_t) :: Model 
INTEGER :: Node
INTEGER, POINTER :: DepthPermutation(:)
REAL(KIND=dp) :: surfaceGrad, shearstress
TYPE(Variable_t), POINTER :: DepthSolution 

!Find variable for flow depth 
DepthSolution  => VariableGet(Model % Variables, 'Depth', .TRUE.)
DepthPermutation  => DepthSolution % Perm 
shearstress = -918.0D00*9.81D00*DepthSolution % Values(DepthPermutation(Node))*surfaceGrad
END FUNCTION StressSIA

!--------------------------------------------------------------------------------------------

FUNCTION invStressSIA( Model, Node, surfaceGrad) RESULT(shearstress)
  USE Types
  USE CoordinateSystems
  USE SolverUtils
  USe ElementDescription
!-------------------------------------------------------------------
IMPLICIT NONE
!-------------------------external variables-----------------------
TYPE(Model_t) :: Model 
INTEGER :: Node
INTEGER, POINTER :: DepthPermutation(:)
REAL(KIND=dp) :: surfaceGrad, shearstress
TYPE(Variable_t), POINTER :: DepthSolution 

!Find variable for flow depth 
DepthSolution  => VariableGet(Model % Variables, 'Depth', .TRUE.)
DepthPermutation  => DepthSolution % Perm 
shearstress = 918.0D00*9.81D00*DepthSolution % Values(DepthPermutation(Node))*surfaceGrad
END FUNCTION invStressSIA
