FUNCTION ISMIPHeinoAccumulation( Model, Node, DummyReal ) RESULT(accumulationRate)
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: DummyReal, accumulationRate
!------------ internal variables----------------------------
  REAL (KIND=dp) :: YCoord, XCoord,ZCoord, Radius

  XCoord = Model % Nodes % x(Node)
  YCoord = Model % Nodes % y(Node)

  Radius = SQRT(XCoord**2 + YCoord**2)

  accumulationRate = -1.0D00*(0.15D00 + 7.5D-08 * Radius)/3.1556926D07

END FUNCTION ISMIPHeinoAccumulation

FUNCTION ISMIPHeinoSurfaceTemp( Model, Node, DummyReal ) RESULT(surfaceTemp)
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: DummyReal, surfaceTemp
!------------ internal variables----------------------------
  REAL (KIND=dp) :: YCoord, XCoord, Radius

  XCoord = Model % Nodes % x(Node)
  YCoord = Model % Nodes % y(Node)

  Radius = SQRT(XCoord**2 + YCoord**2)

  surfaceTemp = 2.3315D02 + (2.5D-09 * (Radius/1.0D03)**3)

END FUNCTION ISMIPHeinoSurfaceTemp

FUNCTION getSlipCoefficient( Model, Node, Temperature ) RESULT(basalSlipCoefficient)
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
  REAL(KIND=dp), ALLOCATABLE :: PressureMeltingPoint(:), Slidingparameter(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: TempName, BedrockType
  LOGICAL ::  GotIt, stat, Jump=.FALSE.

  !---------------
  ! Initialization
  !---------------
  CurrentElementAtBeginning => Model % CurrentElement 
  N = Model % MaxElementNodes 
  ALLOCATE(SlidingParameter(N),&  
       PressureMeltingPoint(N),&
       STAT = istat)
  IF (istat /= 0) THEN
     CALL FATAL('iceproperties (basalSlip)','Allocations failed')
  END IF

  PressureMeltingPoint = 273.16D00 ! default

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
     Slidingparameter(1:NBoundary) = GetReal(BC, 'Sliding Parameter', GotIt)
     IF (.NOT. GotIt) THEN
        CALL WARN('iceproperties (basalSlip)','Keyword >Temperate Slip Coefficient< not found')
        CALL WARN('iceproperties (basalSlip)','Asuming Default 1E05 [1/a]')
        SlidingParameter(1:NBoundary) = 1.0D05/3.1556926D07
     END IF
     SlipCoefficientMax = GetConstReal(BC, 'Max Slip Coefficient', GotIt)
     IF (.NOT. GotIt) THEN
        CALL WARN('iceproperties (basalSlip)','Keyword >Max Slip Coefficient< not found')
        CALL WARN('iceproperties (basalSlip)','Asuming Default 1.0E+16 [kg /(s m^2)]')
        SlipCoefficientMax = 1.0E+16
     END IF
  END IF
  !---------------------------------------------------------------------
  ! get type of bedrock (to set slip law accordingly)
  !---------------------------------------------------------------------
  BedrockType =  GetString(BC ,'Bedrock Type', GotIt)
  !------------------------------
  ! check relative temperature
  !------------------------------
  TempHom = MIN(Temperature - PressureMeltingPoint(ParentElementNode),0.0D00)
  !------------------------------
  ! get the result and check it
  !------------------------------
  SELECT CASE(BedrockType)

  CASE("Hard Rock")

  CASE("Soft Sediment")
       

  CASE DEFAULT 
     STOP
  END SELECT

  IF(basalSlipCoefficient > SlipCoefficientMax) THEN
     WRITE(Message,'(A,e10.4)') 'Applied upper threshold of ', SlipCoefficientMax
     CALL INFO('iceproperties (basalSlip)',Message,Level=9)
     basalSlipCoefficient = SlipCoefficientMax
  END IF
  !------------------------------
  ! clean up
  !------------------------------
  DEALLOCATE(Slidingparameter, PressureMeltingPoint)
END FUNCTION getSlipCoefficient
