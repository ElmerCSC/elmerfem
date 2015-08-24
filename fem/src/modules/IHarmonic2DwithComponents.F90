MODULE CircuitsMod

CONTAINS 

!------------------------------------------------------------------------------
  FUNCTION GetComponentParams(Element) RESULT (ComponentParams)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    
    INTEGER :: i
    TYPE(Element_t), POINTER :: Element
    TYPE(Valuelist_t), POINTER :: BodyParams, ComponentParams
    LOGICAL :: Found
    
    BodyParams => GetBodyParams( Element )
    IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('GetCompParams', 'Component parameters not found')
    
    i = GetInteger(BodyParams, 'Component', Found)

    IF (.NOT. Found) CALL Fatal ('GetComponentParams', 'Body not associated to any Component!')

    ComponentParams => CurrentModel % Components(i) % Values
      
    IF (.NOT. ASSOCIATED(ComponentParams)) CALL Fatal ('CircuitsAndDynamics2DHarmonic', &
                                                         'Component parameters not found!')
    
!------------------------------------------------------------------------------
  END FUNCTION GetComponentParams
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AddComponentsToBodyLists()
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    
    LOGICAL :: Found
    INTEGER :: i, j, k
    
    ! Components and Bodies:
    ! ----------------------  
    INTEGER :: BodyId
    INTEGER, POINTER :: BodyAssociations(:) => Null()
    TYPE(Valuelist_t), POINTER :: BodyParams, ComponentParams
        
    DO i = 1, SIZE(CurrentModel % Components)
      ComponentParams => CurrentModel % Components(i) % Values
      
      IF (.NOT. ASSOCIATED(ComponentParams)) CALL Fatal ('CircuitsAndDynamics2DHarmonic', &
                                                         'Component parameters not found!')
      BodyAssociations => ListGetIntegerArray(ComponentParams, 'Body', Found)
      
      IF (.NOT. Found) CYCLE

      DO j = 1, SIZE(BodyAssociations)
        BodyId = BodyAssociations(j)
        BodyParams => CurrentModel % Bodies(BodyId) % Values
        IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('CircuitsAndDynamics2DHarmonic', &
                                                      'Body parameters not found!')
        k = GetInteger(BodyParams, 'Component', Found)
        IF (Found) CALL Fatal ('CircuitsAndDynamics2DHarmonic', &
                               'Body '//TRIM(i2s(BodyId))//' associated to two components!')
        CALL listAddInteger(BodyParams, 'Component', i)
        BodyParams => Null()
      END DO
    END DO

    DO i = 1, SIZE(CurrentModel % Bodies)
      BodyParams => CurrentModel % Bodies(i) % Values
      IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('CircuitsAndDynamics2DHarmonic', &
                                                    'Body parameters not found!')
      j = GetInteger(BodyParams, 'Component', Found)
      IF (.NOT. Found) CYCLE

      WRITE(Message,'(A,I2,A,I2)') 'Body',i,' associated to Component', j
      CALL Info('CircuitsAndDynamics2DHarmonic',Message,Level=3)
      BodyParams => Null()
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE AddComponentsToBodyLists
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION GetComponentBodyIds(Id) RESULT (BodyIds)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    
    LOGICAL :: Found
    INTEGER :: Id
    INTEGER, POINTER :: BodyIds(:)
    TYPE(Valuelist_t), POINTER :: ComponentParams
    
    ComponentParams => CurrentModel % Components(Id) % Values
    
    IF (.NOT. ASSOCIATED(ComponentParams)) CALL Fatal ('CircuitsAndDynamics2DHarmonic', &
                                                         'Component parameters not found!')
    BodyIds => ListGetIntegerArray(ComponentParams, 'Body', Found)
    IF (.NOT. Found) BodyIds => Null()
    
!------------------------------------------------------------------------------
  END FUNCTION GetComponentBodyIds
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION FindSolverWithKey(key, char_len) RESULT (Solver)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    
    LOGICAL :: Found
    INTEGER :: i, char_len
    TYPE(Solver_t), POINTER :: Solver
    CHARACTER(char_len) :: key
    
    ! Look for the solver we attach the circuit equations to:
    ! -------------------------------------------------------
    Found = .False.
    DO i=1, CurrentModel % NumberOfSolvers
      Solver => CurrentModel % Solvers(i)
      IF(ListCheckPresent(Solver % Values, key)) THEN 
        Found = .True. 
        EXIT
      END IF
    END DO
    
    IF (.NOT. Found) CALL Fatal('FindSolverWithKey', & 
       TRIM(Key)//' keyword not found in any of the solvers!')

!------------------------------------------------------------------------------
  END FUNCTION FindSolverWithKey
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE AllocateCircuitsList()
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    INTEGER :: slen,n_Circuits
    CHARACTER(LEN=MAX_NAME_LEN) :: cmd, name

    ! Read Circuit defintions from MATC:
    ! ----------------------------------
    cmd = "Circuits"
    slen = LEN_TRIM(cmd)
    CALL Matc( cmd, name, slen )
    READ(name(1:slen), *) n_Circuits
    
    CurrentModel%n_Circuits = n_Circuits
    
    ALLOCATE( CurrentModel%CMPLXCircuits(n_Circuits) )

!------------------------------------------------------------------------------
  END SUBROUTINE AllocateCircuitsList
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION CountNofCircVarsOfType(CId, Var_type) RESULT (nofc)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    INTEGER :: nofc, char_len, slen, CId, i
    CHARACTER(LEN=MAX_NAME_LEN) :: Var_type
    CHARACTER(LEN=MAX_NAME_LEN) :: name,cmd
    TYPE(CMPLXCircuit_t), POINTER :: Circuit

    
    Circuit => CurrentModel%CMPLXCircuits(CId)
    
    nofc = 0
    
    char_len = LEN_TRIM(Var_type)
    DO i=1,Circuit % n
      cmd = 'C.'//TRIM(i2s(CId))//'.name.'//TRIM(i2s(i))
      slen = LEN_TRIM(cmd)
      CALL Matc( cmd, name, slen )

      IF(name(1:char_len) == Var_type(1:char_len)) nofc = nofc + 1
    END DO

!------------------------------------------------------------------------------
  END FUNCTION CountNofCircVarsOfType
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION CountNofCircComponents(CId, nofvar) RESULT (nofc)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    INTEGER :: nofc, nofvar, slen, CId, i, j, CompId
    INTEGER :: ComponentIDs(nofvar)
    CHARACTER(LEN=MAX_NAME_LEN) :: name,cmd
    TYPE(CMPLXCircuit_t), POINTER :: Circuit

    nofc = 0
    ComponentIDs = -1
    
    Circuit => CurrentModel%CMPLXCircuits(CId)
    
   
    DO i=1,Circuit % n
      cmd = 'C.'//TRIM(i2s(CId))//'.name.'//TRIM(i2s(i))
      slen = LEN_TRIM(cmd)
      CALL Matc( cmd, name, slen )
      
      IF(name(1:12) == 'i_component(' .OR. name(1:12) == 'v_component(') THEN
        DO j=13,slen
          IF(name(j:j)==')') EXIT 
        END DO
        READ(name(13:j-1),*) CompId
        IF (.NOT. ANY(ComponentIDs == CompID)) nofc = nofc + 1
        ComponentIDs(i) = CompId
      END IF
    
    END DO

!------------------------------------------------------------------------------
  END FUNCTION CountNofCircComponents
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReadCircuitVariables(CId)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    INTEGER :: slen, ComponentId,i,j,CId, CompInd
    CHARACTER(LEN=MAX_NAME_LEN) :: cmd, name
    TYPE(CMPLXCircuit_t), POINTER :: Circuit
    TYPE(CMPLXCircuitVariable_t), POINTER :: CVar

    Circuit => CurrentModel%CMPLXCircuits(CId)
    
    CompInd = 0
    DO i=1,Circuit % n
      cmd = 'C.'//TRIM(i2s(CId))//'.name.'//TRIM(i2s(i))
      slen = LEN_TRIM(cmd)
      CALL Matc( cmd, name, slen )
      Circuit % names(i) = name(1:slen)
      
      CVar => Circuit % CircuitVariables(i)
      CVar % isIvar = .FALSE.
      CVar % isVvar = .FALSE.
      CVar % Component => Null()

      IF(name(1:12) == 'i_component(' .OR. name(1:12) == 'v_component(') THEN
        DO j=13,slen
          IF(name(j:j)==')') EXIT 
        END DO
        READ(name(13:j-1),*) ComponentId
        
        CVar % BodyId = ComponentId
        
        IF ( .NOT. ANY(Circuit % ComponentIds == ComponentId) ) THEN
          CompInd = CompInd + 1
          Circuit % ComponentIds(CompInd) = ComponentId
        END IF

        Cvar % Component => Circuit % Components(CompInd)
        Cvar % Component % ComponentId = ComponentId

        SELECT CASE (name(1:12))
        CASE('i_component(')
          CVar % isIvar = .TRUE.
          CVar % Component % ivar => CVar
        CASE('v_component(')
          CVar % isVvar = .TRUE.
          CVar % Component % vvar => CVar
        CASE DEFAULT
          CALL Fatal('Circuits_Init()', 'Circuit variable should be either i_component or v_component!')
        END SELECT
      ELSE
          CVar % isIvar = .FALSE.
          CVar % isVvar = .FALSE.
          CVar % dofs = 1
          CVar % pdofs = 0
          CVar % BodyId = 0
      END IF
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE ReadCircuitVariables
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION GetNofCircVariables(CId) RESULT(n)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    INTEGER :: CId, n, slen 
    CHARACTER(LEN=MAX_NAME_LEN) :: cmd, name
    TYPE(CMPLXCircuit_t), POINTER :: Circuit

    Circuit => CurrentModel%CMPLXCircuits(CId)
    
    cmd = 'C.'//TRIM(i2s(CId))//'.variables'
    slen = LEN_TRIM(cmd)
    CALL Matc( cmd, name, slen )
      
    READ(name(1:slen), *) Circuit % n

    n = Circuit % n

!------------------------------------------------------------------------------
  END FUNCTION GetNofCircVariables
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AllocateCircuit(CId)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    INTEGER :: slen,CId,n
    CHARACTER(LEN=MAX_NAME_LEN) :: cmd, name
    TYPE(CMPLXCircuit_t), POINTER :: Circuit

    Circuit => CurrentModel%CMPLXCircuits(CId)
    
    n = Circuit % n
    
    ALLOCATE( Circuit % ComponentIds(n), Circuit % names(n) )
    ALLOCATE( Circuit % sourceRe(n), Circuit % sourceIm(n), Circuit % sourcetype(n) )
    ALLOCATE( Circuit % CircuitVariables(n), Circuit % Perm(n) )
    ALLOCATE( Circuit % A(n,n), Circuit % B(n,n), &
              Circuit % Mre(n,n), Circuit % Mim(n,n)  )
    Circuit % ComponentIds = 0
    Circuit % names = ' '

!------------------------------------------------------------------------------
  END SUBROUTINE AllocateCircuit
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE ReadComponents(CId)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    INTEGER :: CId, CompInd
    TYPE(CMPLXCircuit_t), POINTER :: Circuit
    TYPE(CMPLXComponent_t), POINTER :: Comp
    TYPE(Valuelist_t), POINTER :: CompParams
    LOGICAL :: Found
    
    Circuit => CurrentModel%CMPLXCircuits(CId)

    Circuit % CvarDofs = 0
    DO CompInd=1,Circuit % n_comp
      Comp => Circuit % Components(CompInd)
      Comp % nofcnts = 0
!        Comp % ComponentId = Circuits(p) % body(CompInd)
      Comp % BodyIds => GetComponentBodyIds(Comp % ComponentId)

      IF (.NOT. ASSOCIATED(Comp % ivar) ) THEN
        CALL FATAL('Circuits_Init', 'Current Circuit Variable is not found for Component '//TRIM(i2s(Comp % ComponentId)))
      ELSE IF (.NOT. ASSOCIATED(Comp % vvar) ) THEN
        CALL FATAL('Circuits_Init', 'Voltage Circuit Variable is not found for Component '//TRIM(i2s(Comp % ComponentId)))
      END IF

      CompParams => CurrentModel % Components (Comp % ComponentId) % Values
      IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('Circuits_Init', 'Component parameters not found!')
      
      Comp % CoilType = GetString(CompParams, 'Coil Type', Found)
      IF (.NOT. Found) CALL Fatal ('Circuits_Init', 'Coil Type not found!')
      
      Comp % i_multiplier_re = GetConstReal(CompParams, 'Current Multiplier re', Found)
      IF (.NOT. Found) Comp % i_multiplier_re = 0._dp
      Comp % i_multiplier_im = GetConstReal(CompParams, 'Current Multiplier im', Found)
      IF (.NOT. Found) Comp % i_multiplier_im = 0._dp
      
      SELECT CASE (Comp % CoilType) 
      CASE ('stranded')
        Comp % nofturns = GetConstReal(CompParams, 'Number of Turns', Found)
        IF (.NOT. Found) CALL Fatal('Circuits_Init','Number of Turns not found!')

        Comp % ElBoundary = GetInteger(CompParams, 'Electrode Boundary 1', Found)
        IF (.NOT. Found) THEN 
          Comp % ElArea = GetConstReal(CompParams, 'Electrode Area', Found)
          IF (.NOT. Found) THEN
            CALL Fatal('Circuits_Init','Electrode Boundary 1 or Electrode Area not found!')
          END IF
        ELSE
          ! Compute Electrode Area Automatically:
          ! -------------------------------------
!              DO t=1,GetNOFBoundaryElements()
!              Element => GetBoundaryElement(t)
        END IF
        
        Comp % N_j = Comp % nofturns / Comp % ElArea

        ! Stranded coil has current and voltage 
        ! variables (which both have a dof):
        ! ------------------------------------
        Comp % ivar % dofs = 1
        Comp % vvar % dofs = 1
        Comp % ivar % pdofs = 0
        Comp % vvar % pdofs = 0

      CASE ('massive')
        ! Massive coil has current and voltage 
        ! variables (which both have a dof):
        ! ------------------------------------
        Comp % ivar % dofs = 1
        Comp % vvar % dofs = 1
        Comp % ivar % pdofs = 0
        Comp % vvar % pdofs = 0

      CASE ('foil winding')
        Comp % polord = GetInteger(CompParams, 'Foil Winding Voltage Polynomial Order', Found)
        IF (.NOT. Found) Comp % polord = 2

        ! Foil winding has current and voltage 
        ! variables. Current has one dof and 
        ! voltage has a polynom for describing the 
        ! global voltage. The polynom has 1+"polynom order"
        ! dofs. Thus voltage variable has 1+1+"polynom order"
        ! dofs (V=V0+V1*alpha+V2*alpha^2+..):
        ! dofs:
        ! V, V0, V1, V2, ...
        ! ------------------------------------
        Comp % ivar % dofs = 1
        Comp % ivar % pdofs = 0
        Comp % vvar % dofs = Comp % polord + 2
        ! polynom dofs:
        ! -------------
        Comp % vvar % pdofs = Comp % polord + 1

        Comp % coilthickness = GetConstReal(CompParams, 'Coil Thickness', Found)
        IF (.NOT. Found) CALL Fatal('Circuits_Init','Coil Thickness not found!')

        Comp % nofturns = GetConstReal(CompParams, 'Number of Turns', Found)
        IF (.NOT. Found) CALL Fatal('Circuits_Init','Number of Turns not found!')

      END SELECT
      CALL AddVariableToCircuit(Circuit, Comp % ivar, CId)
      CALL AddVariableToCircuit(Circuit, Comp % vvar, CId)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ReadComponents
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE AddVariableToCircuit(Circuit, Variable, k)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    TYPE(CMPLXCircuit_t) :: Circuit
    TYPE(CMPLXCircuitVariable_t) :: Variable
    INTEGER :: Owner=-1, k
    INTEGER, POINTER :: circuit_tot_n => Null()
    
    Circuit_tot_n => CurrentModel%Circuit_tot_n
    
    IF(k==1) THEN
      IF(Owner<=0) Owner = MAX(Parenv % PEs/2,1)
      Owner = Owner - 1
      Variable % Owner = Owner
variable % owner = 0
    ELSE
      IF(Owner<=ParEnv % PEs/2) Owner = ParEnv % PEs
      Owner = Owner - 1
      Variable % Owner = Owner
variable % owner = ParEnv % PEs-1
    END IF

    IF (Circuit % UsePerm) THEN
      Variable % valueId = Circuit % Perm(Circuit_tot_n + 1)
      Variable % ImValueId = Variable % valueId + 1
    ELSE
      Variable % valueId = Circuit_tot_n + 1
      Variable % ImValueId = Circuit_tot_n + 2
    END IF
    
    Circuit_tot_n = Circuit_tot_n + 2*Variable % dofs
!------------------------------------------------------------------------------
  END SUBROUTINE AddVariableToCircuit
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE AddComponentValuesToLists(CId)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    TYPE(CMPLXCircuit_t), POINTER :: Circuit
    TYPE(CMPLXComponent_t), POINTER :: Comp
    TYPE(Valuelist_t), POINTER :: CompParams
    INTEGER :: CId, CompInd
    
    Circuit => CurrentModel%CMPLXCircuits(CId)
    
    DO CompInd=1,Circuit % n_comp
 
      Comp => Circuit % Components(CompInd)   

      CompParams => CurrentModel % Components (Comp % ComponentId) % Values
      IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('Circuits_Init', 'Component Parameters not found!')

      CALL listAddInteger(CompParams, 'Circuit Voltage Variable Id', Comp % vvar % valueId)
      CALL listAddInteger(CompParams, 'Circuit Voltage Variable dofs', Comp % vvar % dofs)
      CALL listAddInteger(CompParams, 'Circuit Current Variable Id', Comp % ivar % valueId)
      CALL listAddInteger(CompParams, 'Circuit Current Variable dofs', Comp % ivar % dofs)
      CALL listAddConstReal(CompParams, 'Stranded Coil N_j', Comp % N_j)
      CurrentModel % Components (Comp % ComponentId) % Values => CompParams
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE AddComponentValuesToLists
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE AddBareCircuitVariables(CId)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    TYPE(CMPLXCircuit_t), POINTER :: Circuit
    TYPE(CMPLXCircuitVariable_t), POINTER :: CVar
    INTEGER :: CId, i
    
    Circuit => CurrentModel%CMPLXCircuits(CId)
    ! add variables that are not associated to components
    DO i=1,Circuit % n
      Cvar => Circuit % CircuitVariables(i)
      IF (Cvar % isIvar .OR. Cvar % isVvar) CYCLE
      CALL AddVariableToCircuit(Circuit, Circuit % CircuitVariables(i), CId)
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE AddBareCircuitVariables
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReadCoefficientMatrices(CId)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    INTEGER :: CId,n
    TYPE(CMPLXCircuit_t), POINTER :: Circuit

    Circuit => CurrentModel%CMPLXCircuits(CId)
    n = Circuit % n
    
    ! Read in the coefficient matrices for the circuit equations:
    ! Ax' + Bx = source:
    ! ------------------------------------------------------------

    CALL matc_get_array('C.'//TRIM(i2s(CId))//'.A'//CHAR(0),Circuit % A,n,n)
    CALL matc_get_array('C.'//TRIM(i2s(CId))//'.B'//CHAR(0),Circuit % B,n,n)
    
    ! Complex multiplier matrix is used for:
    ! B = times(M,B), where B times is the element-wise product
    ! ---------------------------------------------------------
    CALL matc_get_array('C.'//TRIM(i2s(CId))//'.Mre'//CHAR(0),Circuit % Mre,n,n)
    CALL matc_get_array('C.'//TRIM(i2s(CId))//'.Mim'//CHAR(0),Circuit % Mim,n,n)

!------------------------------------------------------------------------------
  END SUBROUTINE ReadCoefficientMatrices
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReadPermutationVector(CId)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    INTEGER :: CId,n,slen,i
    CHARACTER(LEN=MAX_NAME_LEN) :: cmd, name
    TYPE(CMPLXCircuit_t), POINTER :: Circuit

    Circuit => CurrentModel%CMPLXCircuits(CId)
    n = Circuit % n

    DO i=1,n
      cmd = 'C.'//TRIM(i2s(CId))//'.perm('//TRIM(i2s(i-1))//')'
      slen = LEN_TRIM(cmd)
      CALL Matc( cmd, name, slen )
      READ(name(1:slen),*) Circuit % Perm(i)
    END DO
    
    IF(ANY(Circuit % Perm /= 0)) THEN 
      Circuit % UsePerm = .TRUE.
      CALL Info( 'IHarmonic2D','Found Permutation vector for circuit '//i2s(CId), Level=4 )
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE ReadPermutationVector
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReadCircuitSources(CId)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    INTEGER :: CId,n,slen,i
    CHARACTER(LEN=MAX_NAME_LEN) :: cmd, name
    TYPE(CMPLXCircuit_t), POINTER :: Circuit

    Circuit => CurrentModel%CMPLXCircuits(CId)
    n = Circuit % n
    DO i=1,n
      ! Names of the source functions, these functions should be found
      ! in the "Body Force 1" block of the .sif file.
      ! (nc: is for 'no check' e.g. don't abort if the MATC variable is not found!)
      ! ---------------------------------------------------------------------------
      cmd = 'nc:C.'//TRIM(i2s(CId))//'.source.'//TRIM(i2s(i))//'.re'
      slen = LEN_TRIM(cmd)
      CALL Matc( cmd, name, slen )
      Circuit % SourceRe(i) = name(1:slen)
      
      cmd = 'nc:C.'//TRIM(i2s(CId))//'.source.'//TRIM(i2s(i))//'.im'
      slen = LEN_TRIM(cmd)
      CALL Matc( cmd, name, slen )
      Circuit % SourceIm(i) = name(1:slen)

      ! Types of the source functions: voltage or current
      ! (nc: is for 'no check' e.g. don't abort if the MATC variable is not found!)
      ! ---------------------------------------------------------------------------
      cmd = 'nc:C.'//TRIM(i2s(CId))//'.source.'//TRIM(i2s(i))//'.type'
      slen = LEN_TRIM(cmd)
      CALL Matc( cmd, name, slen )
      Circuit % sourcetype(i) = name(1:slen)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ReadCircuitSources
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE WriteCoeffVectorsForCircVariables(CId)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    INTEGER :: CId,n,slen,i,j,RowId
    TYPE(CMPLXCircuit_t), POINTER :: Circuit
    TYPE(CMPLXCircuitVariable_t), POINTER :: Cvar
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
  
    Circuit => CurrentModel%CMPLXCircuits(CId)
    n = Circuit % n

    DO i=1,n
      Cvar => Circuit % CircuitVariables(i)
      RowId = Cvar % ValueId

      ALLOCATE(Cvar % A(Circuit % n), &
               Cvar % B(Circuit % n), &
               Cvar % M(Circuit % n), &
               Cvar % Source(Circuit % n))
      Cvar % A = 0._dp
      Cvar % B = 0._dp
      Cvar % M = 0._dp
      Cvar % Source = 0._dp

      DO j=1,Circuit % n
        IF (Circuit % A(i,j)/=0) Cvar % A(j) = Circuit % A(i,j)
        IF (Circuit % B(i,j)/=0) Cvar % B(j) = Circuit % B(i,j)
        IF (Circuit % Mre(i,j)/=0 .OR. Circuit % Mim(i,j)/=0) &
          Cvar % M(j) = Circuit % Mre(i,j) + im * Circuit % Mim(i,j)
      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE WriteCoeffVectorsForCircVariables
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   FUNCTION IdInList(Id, List) RESULT (T)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     INTEGER :: List(:), Id
     LOGICAL :: T
     T = .FALSE.
     IF (ANY(List == Id)) T = .TRUE.
!------------------------------------------------------------------------------
   END FUNCTION IdInList
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   FUNCTION ElAssocToComp(Element, Component) RESULT (T)
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE
     TYPE(CMPLXComponent_t), POINTER :: Component
     TYPE(Element_t), POINTER :: Element
     LOGICAL :: T
     T = IdInList(Element % BodyId, Component % BodyIds)
!------------------------------------------------------------------------------
   END FUNCTION ElAssocToComp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   FUNCTION ElAssocToCvar(Element, Cvar) RESULT (T)
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE
     TYPE(CMPLXCircuitVariable_t), POINTER :: Cvar
     TYPE(Element_t), POINTER :: Element
     LOGICAL :: T
     T = .FALSE.
     IF (ASSOCIATED(Cvar % Component)) THEN
       IF (ASSOCIATED(Cvar % Component % BodyIds)) &
       T = IdInList(Element % BodyId, Cvar % Component % BodyIds)
     END IF
!------------------------------------------------------------------------------
   END FUNCTION ElAssocToCvar
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION ReIndex(Ind)
!------------------------------------------------------------------------------
    Integer :: Ind, ReIndex

    ReIndex = 2 * Ind - 1
!------------------------------------------------------------------------------
  END FUNCTION ReIndex
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION ImIndex(Ind)
!------------------------------------------------------------------------------
    Integer :: Ind, ImIndex

    ImIndex = 2 * Ind
!------------------------------------------------------------------------------
  END FUNCTION ImIndex
!------------------------------------------------------------------------------

END MODULE CircuitsMod

MODULE CircMatInitMod

CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE SetCircuitsParallelInfo()
!------------------------------------------------------------------------------
    USE DefUtils
    USE CircuitsMod
    IMPLICIT NONE
    TYPE(Matrix_t), POINTER :: CM
    TYPE(CMPLXCircuitVariable_t), POINTER :: Cvar
    TYPE(Solver_t), POINTER :: ASolver
    TYPE(CMPLXCircuit_t), POINTER :: Circuits(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: i, nm, Circuit_tot_n, p, j, &
               cnt(Parenv % PEs), r_cnt(ParEnv % PEs), &
               RowId, nn, l, k, n_Circuits
    
    CM => CurrentModel%CircuitMatrix
    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('SetCircuitsParallelInfo','ASolver not found!')
    nm = ASolver % Matrix % NumberOfRows
    Circuit_tot_n = CurrentModel%Circuit_tot_n
    Circuits => CurrentModel % CMPLXCircuits
    n_Circuits = CurrentModel % n_Circuits
    
    IF(.NOT.ASSOCIATED(CM % ParallelInfo)) THEN
      ALLOCATE(CM % ParallelInfo)
      ALLOCATE(CM % ParallelInfo % NeighbourList(nm+Circuit_tot_n))
      DO i=1,nm+Circuit_tot_n
        CM % ParallelInfo % NeighbourList(i) % Neighbours => Null()
      END DO
    END IF

    DO p = 1,n_Circuits
      DO i=1,Circuits(p) % n
        cnt  = 0
        Cvar => Circuits(p) % CircuitVariables(i)
        IF(ASSOCIATED(CVar%Component)) THEN
          DO j=1,GetNOFACtive()
             Element => GetActiveElement(j)
               IF(ElAssocToCvar(Element, Cvar)) THEN
                 cnt(ParEnv % mype+1)=cnt(ParEnv % mype+1)+1
               END IF
          END DO
        END IF
        CALL MPI_ALLREDUCE(cnt,r_cnt,ParEnv % PEs, MPI_INTEGER, &
                MPI_MAX,ASolver % Matrix % Comm,j)


        RowId = Cvar % ValueId + nm

        nn = COUNT(r_cnt>0)
        IF(r_cnt(Cvar % Owner+1)<=0) nn=nn+1

        DO j=1,Cvar % Dofs
          IF(.NOT.ASSOCIATED(CM % ParallelInfo % NeighbourList(RowId+2*(j-1))%Neighbours)) THEN
            ALLOCATE(CM % ParallelInfo % NeighbourList(RowId+2*(j-1)) % Neighbours(nn))
            ALLOCATE(CM % ParallelInfo % NeighbourList(RowId+2*(j-1)+1) % Neighbours(nn))
          END IF
          CM % ParallelInfo % NeighbourList(RowId+2*(j-1)) % Neighbours(1)   = CVar % Owner
          CM % ParallelInfo % NeighbourList(RowId+2*(j-1)+1) % Neighbours(1) = CVar % Owner
          l = 1
          DO k=0,ParEnv % PEs-1
            IF(k==CVar % Owner) CYCLE
            IF(r_cnt(k+1)>0) THEN
              l = l + 1
              CM % ParallelInfo % NeighbourList(RowId+2*(j-1)) % Neighbours(l) = k
              CM % ParallelInfo % NeighbourList(RowId+2*(j-1) + 1) % Neighbours(l) = k
            END IF
          END DO
          CM % RowOwner(RowId + 2*(j-1))     = Cvar % Owner
          CM % RowOwner(RowId + 2*(j-1) + 1) = Cvar % Owner
        END DO
      END DO
    END DO


    IF ( parenv % mype==0 ) THEN
      DO i=1,parenv % pes
        CALL INFO('SetCircuitsParallelInfo','owners: '//i2s(i)//' '//i2s(count(cm % rowowner==i-1)), Level=9)
      END DO
    END IF
!------------------------------------------------------------------------------
   END SUBROUTINE SetCircuitsParallelInfo
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE CountCmplxMatElement(Rows, Cnts, RowId, dofs)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: Rows(:), Cnts(:)
    INTEGER :: RowId, dofs

    ! Matrix element structure:
    !
    ! Re -Im
    ! Im Re
    !
    ! First do Re -Im:
    ! ----------------
    Cnts(RowId) = Cnts(RowId) + 2 * dofs

    ! Then do Im Re:
    ! --------------
    Cnts(RowId+1) = Cnts(RowId+1) + 2 * dofs

!------------------------------------------------------------------------------
   END SUBROUTINE CountCmplxMatElement
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CreateCmplxMatElement(Rows, Cols, Cnts, RowId, ColId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: Rows(:), Cols(:), Cnts(:)
    INTEGER :: RowId, ColId

    ! Matrix element structure:
    !
    ! Re -Im
    ! Im Re
    !
    ! First do Re (0,0):
    ! ------------------
    Cols(Rows(RowId) + Cnts(RowId)) = ColId
    Cnts(RowId) = Cnts(RowId) + 1

    ! Then do -Im (0,1):
    ! ------------------
    Cols(Rows(RowId) + Cnts(RowId)) = ColId + 1
    Cnts(RowId) = Cnts(RowId) + 1

    ! Then do Re (1,0):
    ! -----------------
    Cols(Rows(RowId+1) + Cnts(RowId+1)) = ColId
    Cnts(RowId+1) = Cnts(RowId+1) + 1

    ! Then do Im (1,1):
    ! -----------------
    Cols(Rows(RowId+1) + Cnts(RowId+1)) = ColId + 1
    Cnts(RowId+1) = Cnts(RowId+1) + 1
    
!------------------------------------------------------------------------------
   END SUBROUTINE CreateCmplxMatElement
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CountBasicCircuitEquations(Rows, Cnts)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    TYPE(CMPLXCircuit_t), POINTER :: Circuits(:)
    TYPE(CMPLXCircuitVariable_t), POINTER :: Cvar
    INTEGER :: i, j, p, nm, RowId, n_Circuits
    INTEGER, POINTER :: Rows(:), Cnts(:)
    
    Circuits => CurrentModel % CMPLXCircuits
    n_Circuits = CurrentModel % n_Circuits
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    
    ! Basic circuit equations...
    ! ---------------------------
    DO p = 1,n_Circuits
      DO i=1,Circuits(p) % n
        Cvar => Circuits(p) % CircuitVariables(i)
        IF(CVar % Owner /= ParEnv % myPE) CYCLE

        RowId = Cvar % ValueId + nm
        DO j=1,Circuits(p) % n
          IF(Cvar % A(j)/=0._dp.OR.Cvar % B(j)/=0._dp) &
             CALL CountCmplxMatElement(Rows, Cnts, RowId, 1)
        END DO
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountBasicCircuitEquations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CreateBasicCircuitEquations(Rows, Cols, Cnts)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    TYPE(CMPLXCircuit_t), POINTER :: Circuits(:)
    TYPE(CMPLXCircuitVariable_t), POINTER :: Cvar
    INTEGER :: i, j, p, nm, RowId, ColId, n_Circuits
    INTEGER, POINTER :: Rows(:), Cols(:), Cnts(:)
    
    Circuits => CurrentModel % CMPLXCircuits
    n_Circuits = CurrentModel % n_Circuits
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    
    ! Basic circuit equations...
    ! ---------------------------
    DO p = 1,n_Circuits
      DO i=1,Circuits(p) % n
        Cvar => Circuits(p) % CircuitVariables(i)
        IF(Cvar % Owner /= ParEnv % myPE) CYCLE

        RowId = Cvar % ValueId + nm
        DO j=1,Circuits(p) % n
          IF(Cvar % A(j)/=0._dp .OR. Cvar % B(j)/=0._dp) THEN
            ColId = Circuits(p) % CircuitVariables(j) % ValueId + nm
            CALL CreateCmplxMatElement(Rows, Cols, Cnts, RowId, ColId)
          END IF
        END DO
      END DO
    END DO

!------------------------------------------------------------------------------
   END SUBROUTINE CreateBasicCircuitEquations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CountComponentEquations(Rows, Cnts, Done, dofsdone)
!------------------------------------------------------------------------------
    USE DefUtils
    USE CircuitsMod
    IMPLICIT NONE
    TYPE(CMPLXCircuit_t), POINTER :: Circuits(:)
    TYPE(CMPLXCircuitVariable_t), POINTER :: Cvar
    TYPE(Solver_t), POINTER :: ASolver
    TYPE(Element_t), POINTER :: Element
    TYPE(CMPLXComponent_t), POINTER :: Comp
    INTEGER :: i, j, p, nm, nn, nd, &
               RowId, ColId, n_Circuits, &
               CompInd, q
    INTEGER, POINTER :: Rows(:), Cnts(:)
    LOGICAL :: dofsdone
    LOGICAL*1 :: Done(:)
    
    Circuits => CurrentModel % CMPLXCircuits
    n_Circuits = CurrentModel % n_Circuits
    Asolver => CurrentModel % Asolver
    nm = Asolver % Matrix % NumberOfRows

    DO p=1,n_Circuits
      DO CompInd=1,Circuits(p) % n_comp
        Done = .FALSE.
        Comp => Circuits(p) % Components(CompInd)
        Cvar => Comp % vvar
        RowId = Cvar % ValueId + nm
        ColId = Cvar % ValueId + nm
        SELECT CASE (Comp % CoilType)
        CASE('stranded')
          IF (Cvar % Owner == ParEnv % myPE) THEN
             CALL CountCmplxMatElement(Rows, Cnts, RowId, 1)
             CALL CountCmplxMatElement(Rows, Cnts, RowId, 1)
          END IF
        CASE('massive')
          IF (CVar % Owner == ParEnv % myPE) THEN
            CALL CountCmplxMatElement(Rows, Cnts, RowId, 1)
            CALL CountCmplxMatElement(Rows, Cnts, RowId, 1)
          END IF
        CASE('foil winding')
          IF (Cvar % Owner == ParEnv % myPE) THEN
            ! V = V0 + V1*alpha + V2*alpha^2 + ...
            CALL CountCmplxMatElement(Rows, Cnts, RowId, Cvar % dofs)

            ! Circuit eqns for the pdofs:
            ! I(Vj) - I = 0
              ! ------------------------------------
            DO j=1, Cvar % pdofs
              CALL CountCmplxMatElement(Rows, Cnts, RowId + 2*j, Cvar % dofs)
            END DO
          END IF
        END SELECT

!        temp = SUM(Cnts)
!print *, "Active elements", ParEnv % Mype, ":", GetNOFActive()
        DO q=GetNOFActive(),1,-1
          Element => GetActiveElement(q)
          IF (ElAssocToComp(Element, Comp)) THEN
            nn = GetElementNOFNodes(Element)
            nd = GetElementNOFDOFs(Element,ASolver)
            SELECT CASE (Comp % CoilType)
            CASE('stranded')              
              CALL CountAndCreateStranded(Element,nn,nd,RowId,Cnts,Done,Rows)
            CASE('massive')
              CALL CountAndCreateMassive(Element,nn,nd,RowId,Cnts,Done,Rows)
            CASE('foil winding')
              DO j = 1, Cvar % pdofs
                dofsdone = ( j==Cvar%pdofs )   
                CALL CountAndCreateFoilWinding(Element,nn,nd,2*j+RowId,Cnts,Done,dofsdone,Rows)
              END DO
            END SELECT
          END IF
        END DO
!        Comp % nofcnts = SUM(Cnts) - temp
!        print *, ParEnv % Mype, "CompInd:", CompInd, "Comp % nofcnts", Comp % nofcnts
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountComponentEquations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CreateComponentEquations(Rows, Cols, Cnts, Done, dofsdone)
!------------------------------------------------------------------------------
    USE DefUtils
    USE CircuitsMod
    IMPLICIT NONE
    TYPE(CMPLXCircuit_t), POINTER :: Circuits(:)
    TYPE(CMPLXCircuitVariable_t), POINTER :: Cvar
    TYPE(Solver_t), POINTER :: ASolver
    TYPE(Element_t), POINTER :: Element
    TYPE(CMPLXComponent_t), POINTER :: Comp
    INTEGER :: i, j, jj, p, nm, nn, nd, &
               RowId, ColId, n_Circuits, &
               CompInd, q
    INTEGER, POINTER :: Rows(:), Cols(:), Cnts(:)
    LOGICAL :: dofsdone
    LOGICAL*1 :: Done(:)
    
    Circuits => CurrentModel % CMPLXCircuits
    n_Circuits = CurrentModel % n_Circuits
    Asolver => CurrentModel % Asolver
    nm = Asolver % Matrix % NumberOfRows

    DO p=1,n_Circuits
      DO CompInd = 1, Circuits(p) % n_comp
        Done = .FALSE.
        Comp => Circuits(p) % Components(CompInd)
        Cvar => Comp % vvar
        RowId = Comp % vvar % ValueId + nm
        ColId = Comp % ivar % ValueId + nm

        SELECT CASE (Comp % CoilType)
        CASE('stranded')
          IF (Cvar % Owner == ParEnv % myPE) THEN
            CALL CreateCmplxMatElement(Rows, Cols, Cnts, RowId, ColId)
            CALL CreateCmplxMatElement(Rows, Cols, Cnts, RowId, RowId)
          END IF
        CASE('massive')
          IF (Cvar % Owner == ParEnv % myPE) THEN
            CALL CreateCmplxMatElement(Rows, Cols, Cnts, RowId, ColId)
            CALL CreateCmplxMatElement(Rows, Cols, Cnts, RowId, RowId)
          END IF
        CASE('foil winding')
          DO j=0, Cvar % pdofs
            IF (Cvar % Owner == ParEnv % mype) THEN
              ! V = V0 + V1*alpha + V2*alpha^2 + ...
              CALL CreateCmplxMatElement(Rows, Cols, Cnts, RowId, RowId + 2*j)
              IF (j/=0) THEN
                ! Circuit eqns for the pdofs:
                ! I(Vi) - I = 0
                ! ------------------------------------
                CALL CreateCmplxMatElement(Rows, Cols, Cnts, RowId + 2*j, ColId)
                DO jj = 1, Cvar % pdofs
                    CALL CreateCmplxMatElement(Rows, Cols, Cnts, RowId + 2*j, RowId + 2*jj)
                END DO
              END IF
            END IF
          END DO
        END SELECT

!        temp = SUM(Cnts)
!print *, "Active elements ", ParEnv % Mype, ":", GetNOFActive()
        DO q=GetNOFActive(),1,-1
          Element => GetActiveElement(q)
          IF (ElAssocToComp(Element, Comp)) THEN
            nn = GetElementNOFNodes(Element)
            nd = GetElementNOFDOFs(Element,ASolver)
            SELECT CASE (Comp % CoilType)
            CASE('stranded')              
              CALL CountAndCreateStranded(Element,nn,nd,RowId,Cnts,Done,Rows,Cols,ColId)
            CASE('massive')
              CALL CountAndCreateMassive(Element,nn,nd,RowId,Cnts,Done,Rows,Cols=Cols)
            CASE('foil winding')
              DO j = 1, Cvar % pdofs
                dofsdone = ( j==Cvar%pdofs )
                CALL CountAndCreateFoilWinding(Element,nn,nd,2*j+RowId,Cnts,Done,dofsdone,Rows,Cols=Cols)
              END DO
            END SELECT
          END IF
        END DO
!        Comp % nofcnts = SUM(Cnts) - temp
!        print *, ParEnv % Mype, "CompInd:", CompInd, "Coil Type:", Comp % CoilType, &
!                 "Comp % BodyId:", Comp % BodyId, "Comp % nofcnts", Comp % nofcnts

      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CreateComponentEquations
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE CountAndCreateStranded(Element,nn,nd,i,Cnts,Done,Rows,Cols,Jsind)
!------------------------------------------------------------------------------
    USE DefUtils
    USE CircuitsMod
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    INTEGER :: nn, nd
    OPTIONAL :: Cols
    INTEGER :: Rows(:), Cols(:), Cnts(:)
    INTEGER :: p,i,j,Indexes(nd)
    INTEGER, OPTIONAL :: Jsind
    INTEGER, POINTER :: PS(:)
    LOGICAL*1 :: Done(:)
  
    IF (.NOT. ASSOCIATED(CurrentModel % ASolver) ) CALL Fatal ('CountAndCreateStranded','ASolver not found!')
    PS => CurrentModel % Asolver % Variable % Perm
    nd = GetElementDOFs(Indexes,Element,CurrentModel % ASolver)
    DO p=1,nd
      j = Indexes(p)
      IF(.NOT.Done(j)) THEN
        Done(j) = .TRUE.
        j = ReIndex(PS(j))
        IF(PRESENT(Cols)) THEN
          CALL CreateCmplxMatElement(Rows, Cols, Cnts, i, j) 
          CALL CreateCmplxMatElement(Rows, Cols, Cnts, j, Jsind)
!          CALL CreateCmplxMatElement(Rows, Cols, Cnts, j, Jsind)
        ELSE
          CALL CountCmplxMatElement(Rows, Cnts, i, 1)
          CALL CountCmplxMatElement(Rows, Cnts, j, 1)
!          CALL CountCmplxMatElement(Rows, Cnts, j, 1)
        END IF
      END IF
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountAndCreateStranded
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CountAndCreateMassive(Element,nn,nd,i,Cnts,Done,Rows,Cols)
!------------------------------------------------------------------------------
    USE DefUtils
    USE CircuitsMod
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    INTEGER :: nn, nd
    OPTIONAL :: Cols
    INTEGER :: Rows(:), Cols(:), Cnts(:)
    INTEGER :: p,i,j,Indexes(nd)
    INTEGER, POINTER :: PS(:)
    LOGICAL*1 :: Done(:)
    
    IF (.NOT. ASSOCIATED(CurrentModel % ASolver) ) CALL Fatal ('CountAndCreateMassive','ASolver not found!')
    PS => CurrentModel % Asolver % Variable % Perm
    nd = GetElementDOFs(Indexes,Element,CurrentModel % ASolver)
    DO p=1,nd
      j = Indexes(p)
      IF(.NOT.Done(j)) THEN
        Done(j) = .TRUE.
        j = ReIndex(PS(j))
        IF(PRESENT(Cols)) THEN
          CALL CreateCmplxMatElement(Rows, Cols, Cnts, i, j)
          CALL CreateCmplxMatElement(Rows, Cols, Cnts, j, i)
        ELSE
          CALL CountCmplxMatElement(Rows, Cnts, i, 1)
          CALL CountCmplxMatElement(Rows, Cnts, j, 1)
        END IF
      END IF
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountAndCreateMassive
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CountAndCreateFoilWinding(Element,nn,nd,i,Cnts,Done,dofsdone,Rows,Cols)
!------------------------------------------------------------------------------
    USE DefUtils
    USE CircuitsMod
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    INTEGER :: nn, nd
    OPTIONAL :: Cols
    INTEGER :: Rows(:), Cols(:), Cnts(:)
    INTEGER :: p,i,j,Indexes(nd)
    LOGICAL :: dofsdone
    INTEGER, POINTER :: PS(:)
    LOGICAL*1 :: Done(:)
    
    IF (.NOT. ASSOCIATED(CurrentModel % ASolver) ) CALL Fatal ('CountAndCreateFoilWinding','ASolver not found!')
    PS => CurrentModel % Asolver % Variable % Perm
    nd = GetElementDOFs(Indexes,Element,CurrentModel % ASolver)
    DO p=1,nd
      j = Indexes(p)
      IF(.NOT.Done(j)) THEN
        Done(j) = dofsdone
        j = ReIndex(PS(j))
        IF(PRESENT(Cols)) THEN
          CALL CreateCmplxMatElement(Rows, Cols, Cnts, i, j)
          CALL CreateCmplxMatElement(Rows, Cols, Cnts, j, i)
        ELSE
          CALL CountCmplxMatElement(Rows, Cnts, i, 1)
          CALL CountCmplxMatElement(Rows, Cnts, j, 1)
        END IF
      END IF
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountAndCreateFoilWinding
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE Circuits_MatrixInit()
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    TYPE(Matrix_t), POINTER :: CM
    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:), Rows(:), Cols(:), Cnts(:)
    INTEGER :: nm, Circuit_tot_n, n, i
    LOGICAL :: dofsdone
    LOGICAL*1, ALLOCATABLE :: Done(:)
    REAL(KIND=dp), POINTER :: Values(:)

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Circuits_MatrixInit','ASolver not found!')
    Circuit_tot_n = CurrentModel%Circuit_tot_n
    
    ! Initialialize Circuit matrix:
    ! -----------------------------
    PS => Asolver % Variable % Perm
    nm =  Asolver % Matrix % NumberOfRows

    CM => AllocateMatrix()
    CurrentModel%CircuitMatrix=>CM
    
    CM % Format = MATRIX_CRS
    Asolver %  Matrix % AddMatrix => CM
    ALLOCATE(CM % RHS(nm + Circuit_tot_n)); CM % RHS=0._dp

    CM % NumberOfRows = nm + Circuit_tot_n
    n = CM % NumberOfRows
    ALLOCATE(Rows(n+1), Cnts(n)); Rows=0; Cnts=0
    ALLOCATE(Done(nm), CM % RowOwner(n)); Cm % RowOwner=-1

    CALL SetCircuitsParallelInfo()

    ! COUNT SIZES:
    ! ============
    dofsdone = .FALSE.
    
    CALL CountBasicCircuitEquations(Rows, Cnts)
    CALL CountComponentEquations(Rows, Cnts, Done, dofsdone)

    ! ALLOCATE CRS STRUCTURES (if need be):
    ! =====================================

    n = SUM(Cnts)

    IF (n<=0) THEN
      CM % NUmberOfRows = 0
      DEALLOCATE(Rows,Cnts,Done,CM); CM=>Null(); RETURN
    END IF

    ALLOCATE(Cols(n+1), Values(n+1))
    Cols = 0; Values = 0._dp

    ! CREATE ROW POINTERS:
    ! ====================

    CM % NumberOfRows = nm + Circuit_tot_n
    Rows(1) = 1
    DO i=2,CM % NumberOfRows+1
      Rows(i) = Rows(i-1) + Cnts(i-1)
    END DO

    Cnts = 0

    ! CREATE COLMUNS:
    ! ===============

    CALL CreateBasicCircuitEquations(Rows, Cols, Cnts)
    CALL CreateComponentEquations(Rows, Cols, Cnts, Done, dofsdone)
    

    IF (n /= SUM(Cnts)) THEN
      print *, "Counted Cnts:", n, "Applied Cnts:", SUM(Cnts)
      CALL Fatal('Circuits_MatrixInit', &
                 'There were different amount of matrix elements than was counted')
    END IF

    DEALLOCATE( Cnts, Done )
    CM % Rows => Rows
    CM % Cols => Cols
    CM % Values => Values
    CALL CRS_SortMatrix(CM)
    
    Asolver %  Matrix % AddMatrix => CM
!------------------------------------------------------------------------------
  END SUBROUTINE Circuits_MatrixInit
!------------------------------------------------------------------------------
END MODULE CircMatInitMod

MODULE CircEqWriteMod

CONTAINS
  
!------------------------------------------------------------------------------
   SUBROUTINE AddBasicCircuitEquations(p)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    TYPE(CMPLXCircuit_t), POINTER :: Circuit
    TYPE(CMPLXCircuitVariable_t), POINTER :: Cvar
    TYPE(ValueList_t), POINTER :: BF
    TYPE(Matrix_t), POINTER :: CM
    INTEGER :: p, i, nm, RowId, ColId, j
    REAL(KIND=dp) :: Omega, vphi
    COMPLEX(KIND=dp) :: cmplx_value
    LOGICAL :: Found
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    
    Circuit => CurrentModel % CMPLXCircuits(p)
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    Omega = GetAngularFrequency()
    BF => CurrentModel % BodyForces(1) % Values
    CM => CurrentModel%CircuitMatrix
    
    DO i=1,Circuit % n
      Cvar => Circuit % CircuitVariables(i)

      IF(Cvar % Owner /= ParEnv % myPE) CYCLE

      RowId = Cvar % ValueId + nm
      
      vphi=0._dp
      IF ( ASSOCIATED(BF) ) &
        vphi = GetCReal(BF, Circuit % SourceRe(i), Found)
      
      IF (Found) THEN 
        Cvar % Source(i) = vphi
      END IF
      
      IF (ASSOCIATED(BF) ) &
        vphi = GetCReal(BF, Circuit % SourceIm(i), Found)

      IF (Found) THEN 
        Cvar % Source(i) = Cvar % Source(i) + im * vphi
      END IF
      
      CM % RHS(RowId) = REAL(Cvar % Source(i))
      CM % RHS(RowId+1) = AIMAG(Cvar % Source(i))
        
      DO j=1,Circuit % n

        ColId = Circuit % CircuitVariables(j) % ValueId + nm

        ! - im * Omega * A x: (x could be voltage or current):
        !--------------------------------------------
        IF(Cvar % A(j) /= 0._dp) THEN
          CALL AddToCmplxMatrixElement(CM, RowId, ColId, 0._dp, -Omega * Cvar % A(j))
        END IF
        ! B x:
        ! ------
        IF(Cvar % B(j) /= 0._dp) THEN
          
          IF (Cvar % M(j) /= 0._dp) THEN
            cmplx_value = Cvar % M(j) * Cvar % B(j)
          ELSE
            cmplx_value = Cvar % B(j)
          END IF
          
          CALL AddToCmplxMatrixElement(CM, RowId, ColId, REAL(cmplx_value), AIMAG(cmplx_value))
        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE AddBasicCircuitEquations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE AddToCmplxMatrixElement(CM, RowId, ColId, Re, Im)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    TYPE(Matrix_t), POINTER :: CM
    INTEGER :: RowId, ColId
    REAL(KIND=dp) :: Re, Im

    CALL AddToMatrixElement(CM, RowId, ColId, Re)
    CALL AddToMatrixElement(CM, RowId, ColId+1, -Im)
    CALL AddToMatrixElement(CM, RowId+1, ColId, Im)
    CALL AddToMatrixElement(CM, RowId+1, ColId+1, Re)

!------------------------------------------------------------------------------
  END SUBROUTINE AddToCmplxMatrixElement
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE AddComponentEquationsAndCouplings(p, nn, CSymmetry)
!------------------------------------------------------------------------------
    USE DefUtils
    USE CircuitsMod
    IMPLICIT NONE
    INTEGER :: p, CompInd, nm, nn, nd
    TYPE(Solver_t), POINTER :: ASolver
    TYPE(CMPLXCircuit_t), POINTER :: Circuit
    TYPE(Matrix_t), POINTER :: CM
    TYPE(CMPLXComponent_t), POINTER :: Comp
    TYPE(CMPLXCircuitVariable_t), POINTER :: Cvar
    TYPE(Valuelist_t), POINTER :: CompParams
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: Omega
    INTEGER :: VvarId, IvarId, q, j
    COMPLEX(KIND=dp) :: i_multiplier, cmplx_value
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    COMPLEX(KIND=dp) :: Tcoef(3,3,nn), RotM(3,3,nn)
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType
    LOGICAL :: Found, CSymmetry

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('AddComponentEquationsAndCouplings','ASolver not found!')

    Circuit => CurrentModel % CMPLXCircuits(p)
    nm = Asolver % Matrix % NumberOfRows
    Omega = GetAngularFrequency()
    CM => CurrentModel%CircuitMatrix

    DO CompInd = 1, Circuit % n_comp
      Comp => Circuit % Components(CompInd)
      Cvar => Comp % vvar
      vvarId = Comp % vvar % ValueId + nm
      IvarId = Comp % ivar % ValueId + nm

      IF ( Cvar % Owner == ParEnv % myPE ) THEN
        SELECT CASE (Comp % CoilType)
        CASE('stranded')
          CALL AddToCmplxMatrixElement(CM, VvarId, VvarId, 1._dp, 0._dp)
        CASE('massive')
          i_multiplier = Comp % i_multiplier_re + im * Comp % i_multiplier_im
          IF (i_multiplier /= 0_dp) THEN
            CALL AddToCmplxMatrixElement(CM, VvarId, IvarId, -REAL(i_multiplier), -AIMAG(i_multiplier))
          ELSE
            CALL AddToCmplxMatrixElement(CM, VvarId, IvarId, -1._dp, 0._dp)
          END IF
        CASE('foil winding')
          ! Foil Winding voltage: 
          ! V + ...added next... = 0
          ! ----------------------
          i_multiplier = Comp % i_multiplier_re + im * Comp % i_multiplier_im
          CALL AddToCmplxMatrixElement(CM, VvarId, VvarId, 1._dp, 0._dp)
          
          IF (i_multiplier == 0_dp) i_multiplier = 1.0_dp
          
          DO j = 1, Cvar % pdofs 
            ! Foil Winding voltage: 
            !  ... - Nf/Lalpha * int_0^{Lalpha}(V_0+V_1*alpha+V_2*alpha**2+...) = 0
            !          => ... - Nf * (V_0*Lalpha^0 + V_1/2*Lalpha^1 + V_2/3*Lalpha^2 + ...) = 0
            ! where V_m is the mth dof of the polynomial
            ! --------------------------------------------------------------
            cmplx_value = -i_multiplier * REAL(Comp % nofturns) / REAL(j) * Comp % coilthickness**(j-1)
            CALL AddToCmplxMatrixElement(CM, VvarId, 2*j + VvarId, &
                REAL(cmplx_value), AIMAG(cmplx_value))

            ! Circuit eqns for the pdofs:
            ! - Nf/Lalpha * I * int_0^1(Vi'(alpha)) + ...added later... = 0
            ! ----------------------------------------------------------
            CALL AddToCmplxMatrixElement(CM, 2*j + VvarId, IvarId, &
               REAL(cmplx_value), AIMAG(cmplx_value))
          END DO
        END SELECT
      END IF

      DO q=GetNOFActive(),1,-1
        Element => GetActiveElement(q)
        IF (ElAssocToComp(Element, Comp)) THEN
          CompParams => GetComponentParams( Element )
          IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('Circuits_apply', 'Component parameters not found')

          CoilType = GetString(CompParams, 'Coil Type', Found)
          IF (.NOT. Found) CoilType = ''
          
          nn = GetElementNOFNodes(Element)
          nd = GetElementNOFDOFs(Element,ASolver)
          CALL GetConductivity(Element, Tcoef, nn)
          SELECT CASE(CoilType)
          CASE ('stranded')
            CALL Add_stranded(Element,Tcoef(1,1,1:nn),Comp,nn,nd,CSymmetry)
          CASE ('massive')
            CALL Add_massive(Element,Tcoef(1,1,1:nn),Comp,nn,nd,CSymmetry)
          CASE ('foil winding')
!              CALL GetElementRotM(Element, RotM, nn)
            CALL Add_foil_winding(Element,Tcoef(1,1,1:nn),Comp,nn,nd,CSymmetry)
          CASE DEFAULT
            CALL Fatal ('Circuits_apply', 'Non existent Coil Type Chosen!')
          END SELECT
        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE AddComponentEquationsAndCouplings
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE GetConductivity(Element, Tcoef, nn)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    TYPE(Valuelist_t), POINTER :: Material
    COMPLEX(KIND=dp) :: Tcoef(3,3,nn)
    REAL(KIND=dp), POINTER, SAVE :: Cwrk(:,:,:), Cwrk_im(:,:,:) 
    INTEGER :: nn, i, j
    LOGICAL, SAVE :: visited = .FALSE., Found

    IF (.NOT. visited) THEN
      NULLIFY( Cwrk, Cwrk_im )
    END IF

    Tcoef = cmplx(0.0d0,0.0d0)
    Material => GetMaterial( Element )
    IF (.NOT. ASSOCIATED(Material)) CALL Fatal('Circuits_apply','Material not found.')

    CALL ListGetRealArray( Material, &
           'Electric Conductivity', Cwrk, nn, Element % NodeIndexes, Found )

    IF (.NOT. Found) CALL Fatal('Circuits_apply', 'Electric Conductivity not found.')

    IF (Found) THEN
       IF ( SIZE(Cwrk,1) == 1 ) THEN
          DO i=1,3
             Tcoef( i,i,1:nn ) = Cwrk( 1,1,1:nn )
          END DO
       ELSE IF ( SIZE(Cwrk,2) == 1 ) THEN
          DO i=1,MIN(3,SIZE(Cwrk,1))
             Tcoef(i,i,1:nn) = Cwrk(i,1,1:nn)
          END DO
       ELSE
          DO i=1,MIN(3,SIZE(Cwrk,1))
             DO j=1,MIN(3,SIZE(Cwrk,2))
                Tcoef( i,j,1:nn ) = Cwrk(i,j,1:nn)
             END DO
          END DO
       END IF
    END IF

    CALL ListGetRealArray( Material, &
           'Electric Conductivity im', Cwrk_im, nn, Element % NodeIndexes, Found )

    IF (Found) THEN
       IF ( SIZE(Cwrk_im,1) == 1 ) THEN
          DO i=1,3
             Tcoef( i,i,1:nn ) = CMPLX( REAL(Tcoef( i,i,1:nn )), Cwrk_im( 1,1,1:nn ), KIND=dp)
          END DO
       ELSE IF ( SIZE(Cwrk_im,2) == 1 ) THEN
          DO i=1,MIN(3,SIZE(Cwrk_im,1))
             Tcoef(i,i,1:nn) = CMPLX( REAL(Tcoef( i,i,1:nn )), Cwrk_im( i,1,1:nn ), KIND=dp)
          END DO
       ELSE
          DO i=1,MIN(3,SIZE(Cwrk_im,1))
             DO j=1,MIN(3,SIZE(Cwrk_im,2))
                Tcoef( i,j,1:nn ) = CMPLX( REAL(Tcoef( i,j,1:nn )), Cwrk_im( i,j,1:nn ), KIND=dp)
             END DO
          END DO
       END IF
    END IF


!------------------------------------------------------------------------------
  END SUBROUTINE GetConductivity
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetElementRotM(Element,RotM,n)
!------------------------------------------------------------------------------
   USE DefUtils
   IMPLICIT NONE
   TYPE(Element_t) :: Element
   INTEGER :: k, l, m, j, n
   REAL(KIND=dp) :: RotM(3,3,n)
   TYPE(Variable_t), POINTER, SAVE :: RotMvar
   LOGICAL, SAVE :: visited = .FALSE.
   INTEGER, PARAMETER :: ind1(9) = [1,1,1,2,2,2,3,3,3]
   INTEGER, PARAMETER :: ind2(9) = [1,2,3,1,2,3,1,2,3]
  
   IF(.NOT. visited) THEN
     visited = .TRUE.
     RotMvar => VariableGet( CurrentModel % Mesh % Variables, 'RotM E')
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
   SUBROUTINE Add_stranded(Element,Tcoef,Comp,nn,nd,CSymmetry)
!------------------------------------------------------------------------------
    USE DefUtils
    USE CircuitsMod
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    COMPLEX(KIND=dp) :: Tcoef(nn)
    TYPE(CMPLXComponent_t) :: Comp
    INTEGER :: nn, nd, nm, Indexes(nd),VvarId,IvarId
    LOGICAL :: CSymmetry

    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Omega
    TYPE(Nodes_t), SAVE :: Nodes
    REAL(KIND=dp) :: Basis(nn), DetJ, x
    REAL(KIND=dp) :: dBasisdx(nn,3), wBase(nn), w(3)
    COMPLEX(KIND=dp) :: localC, i_multiplier, cmplx_value
    INTEGER :: p,t
    LOGICAL :: stat

    TYPE(GaussIntegrationPoints_t) :: IP
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
        
    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_stranded','ASolver not found!')
    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitMatrix
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    Omega = GetAngularFrequency()
    
    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)

    VvarId = Comp % vvar % ValueId + nm
    IvarId = Comp % ivar % ValueId + nm

    i_multiplier = Comp % i_multiplier_re + im * Comp % i_multiplier_im


    ! Numerical integration:
    ! ----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis,dBasisdx )

      IF( CSymmetry ) THEN
        x = SUM( Basis(1:nn) * Nodes % x(1:nn) )
        detJ = detJ * x
      END IF

      w = [0._dp, 0._dp, 1._dp]
      localC = SUM(Tcoef(1:nn) * Basis(1:nn))

      ! I * R, where 
      ! R = (1/sigma * js,js):
      ! ----------------------
      
      CALL AddToCmplxMatrixElement(CM, VvarId, IvarId, &
            REAL(Comp % N_j * IP % s(t)*detJ*SUM(w*w)/localC), &
           AIMAG(Comp % N_j * IP % s(t)*detJ*SUM(w*w)/localC))
            
      DO p=1,nd

        IF (Comp % N_j/=0._dp) THEN
          ! ( im * Omega a,w )
          CALL AddToCmplxMatrixElement(CM, VvarId, ReIndex(PS(Indexes(p))), &
                 REAL(-im * Omega * Comp % N_j * IP % s(t)*detJ*Basis(p)/localC), & 
                AIMAG(-im * Omega * Comp % N_j * IP % s(t)*detJ*Basis(p)/localC))

            IF (i_multiplier /= 0._dp) THEN
              cmplx_value = -i_multiplier*Comp % N_j*IP % s(t)*detJ*Basis(p)
            ELSE
              cmplx_value = -Comp % N_j*IP % s(t)*detJ*Basis(p)
            END IF
            
            CALL AddToCmplxMatrixElement(CM,ReIndex(PS(Indexes(p))), IvarId, &
               REAL(cmplx_value), AIMAG(cmplx_value))

        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Add_stranded
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE Add_massive(Element,Tcoef,Comp,nn,nd,CSymmetry)
!------------------------------------------------------------------------------
    USE DefUtils
    USE CircuitsMod
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    COMPLEX(KIND=dp) :: Tcoef(nn)
    TYPE(CMPLXComponent_t) :: Comp
    LOGICAL :: CSymmetry

    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Omega
    REAL(KIND=dp) :: Basis(nn), DetJ, x
    REAL(KIND=dp) :: dBasisdx(nn,3), wBase(nn)
    COMPLEX(KIND=dp) :: localC
    INTEGER :: nn, nd, j, t, nm, Indexes(nd),VvarId
    LOGICAL :: stat
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_massive','ASolver not found!')
    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitMatrix
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    Omega = GetAngularFrequency()

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)
    CALL GetLocalSolution(Wbase,'W')

    vvarId = Comp % vvar % ValueId + nm

    ! Numerical integration:
    ! ----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis,dBasisdx )

      IF( CSymmetry ) THEN
        x = SUM( Basis(1:nn) * Nodes % x(1:nn) )
        detJ = detJ * x
      END IF

      localC = SUM(Tcoef(1:nn) * Basis(1:nn))

      ! computing the source term Vi(sigma grad v0, grad si):
      ! ------------------------------------------------
      CALL AddToCmplxMatrixElement(CM, vvarId, vvarId, &
              REAL(IP % s(t)*detJ*localC), &
              AIMAG(IP % s(t)*detJ*localC))

      DO j=1,nd
        ! computing the mass term (sigma * im * Omega * a, grad si):
        ! ---------------------------------------------------------
        CALL AddToCmplxMatrixElement(CM, vvarId, ReIndex(PS(Indexes(j))), &
               REAL(im * Omega * IP % s(t)*detJ*localC*basis(j)), &
              AIMAG(im * Omega * IP % s(t)*detJ*localC*basis(j)))

          CALL AddToCmplxMatrixElement(CM, ReIndex(PS(indexes(j))), vvarId, &
                REAL(IP % s(t)*detJ*localC*basis(j)), &
               AIMAG(IP % s(t)*detJ*localC*basis(j)))

      END DO
    END DO

!------------------------------------------------------------------------------
   END SUBROUTINE Add_massive
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE Add_foil_winding(Element,Tcoef,Comp,nn,nd,CSymmetry)
!------------------------------------------------------------------------------
    USE DefUtils
    USE CircuitsMod
    IMPLICIT NONE
    INTEGER :: nn, nd
    TYPE(Element_t) :: Element
    COMPLEX(KIND=dp) :: Tcoef(:), C, value
    TYPE(CMPLXComponent_t) :: Comp
    LOGICAL :: CSymmetry

    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Basis(nn), DetJ, Omega, &
                     localAlpha, localV, localVtest, x, circ_eq_coeff, grads_coeff
    REAL(KIND=dp) :: dBasisdx(nn,3),alpha(nn)
    INTEGER :: nm,p,j,t,Indexes(nd),vvarId,vpolord_tot, &
               vpolord, vpolordtest, dofId, dofIdtest
    LOGICAL :: stat
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_foil_winding','ASolver not found!')
    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitMatrix
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    Omega = GetAngularFrequency()

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)
    CALL GetLocalSolution(alpha,'Alpha')

    vvarId = Comp % vvar % ValueId
    vpolord_tot = Comp % vvar % pdofs - 1

    ! Numerical integration:
    ! ----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis,dBasisdx )

      grads_coeff = 1._dp
      circ_eq_coeff = 1._dp
      IF( CSymmetry ) THEN
        x = SUM( Basis(1:nn) * Nodes % x(1:nn) )
        detJ = detJ * x
        grads_coeff = grads_coeff/(2._dp*pi*x)
        circ_eq_coeff = 2._dp * pi
      END IF

      localAlpha = SUM(alpha(1:nn) * Basis(1:nn))
      
      ! alpha is normalized to be in [0,1] thus, 
      ! it needs to be multiplied by the thickness of the coil 
      ! to get the real alpha:
      ! ------------------------------------------------------
      localAlpha = localAlpha * Comp % coilthickness


      ! Compute the conductivity tensor
      ! -------------------------------
!      DO i=1,3
!        DO j=1,3
          C = SUM( Tcoef(1:nn) * Basis(1:nn) )
!          RotMLoc(i,j) = SUM( RotM(i,j,1:nn) * Basis(1:nn) )
!        END DO
!      END DO

      ! Transform the conductivity tensor:
      ! ----------------------------------
!      C = MATMUL(MATMUL(RotMLoc, C),TRANSPOSE(RotMLoc))

      DO vpolordtest=0,vpolord_tot ! V'(alpha)
        localVtest = localAlpha**vpolordtest
        dofIdtest = 2*(vpolordtest + 1) + vvarId
        DO vpolord = 0, vpolord_tot ! V(alpha)

          localV = localAlpha**vpolord
          dofId = 2*(vpolord + 1) + vvarId
          
          ! Computing the stiff term (sigma V(alpha) grad v0, V'(alpha) grad si):
          ! ---------------------------------------------------------------------
          value = IP % s(t)*detJ*localV*localVtest*C*grads_coeff**2*circ_eq_coeff
          CALL AddToCmplxMatrixElement(CM, dofIdtest+nm, dofId+nm, REAL(value), AIMAG(value))
        END DO

        DO j=1,nd
          ! computing the mass term (sigma * im * Omega * a, V'(alpha) grad si):
          ! ---------------------------------------------------------
          value = im * Omega * IP % s(t)*detJ*localVtest*C*basis(j)*grads_coeff*circ_eq_coeff
          CALL AddToCmplxMatrixElement(CM, dofIdtest+nm, ReIndex(PS(Indexes(j))), REAL(value), AIMAG(value) )
        END DO

      END DO

      DO vpolord = 0, vpolord_tot ! V(alpha)
        localV = localAlpha**vpolord
        dofId = 2*(vpolord + 1) + vvarId

        DO j=1,nd
            value = IP % s(t)*detJ*localV*C*basis(j)*grads_coeff
            CALL AddToCmplxMatrixElement(CM, ReIndex(PS(indexes(j))), dofId+nm, REAL(value), AIMAG(value))
        END DO
      END DO

    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Add_foil_winding
!------------------------------------------------------------------------------



END MODULE CircEqWriteMod


!------------------------------------------------------------------------------
!> Initialization for the primary solver: CurrentSource
!------------------------------------------------------------------------------
SUBROUTINE CircuitsAndDynamics2DHarmonic_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params

  Params => Solver % Values

  ! When we introduce the variables in this way the variables are created
  ! so that they exist when the proper simulation cycle starts.
  ! This also keeps the command file cleaner.
  CALL ListAddString( Params,'Exported Variable 1',&
      '-global Rotor Angle')
  CALL ListAddString( Params,'Exported Variable 2',&
      '-global Rotor Velo')
  Solver % Values => Params

!------------------------------------------------------------------------------
END SUBROUTINE CircuitsAndDynamics2DHarmonic_init
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
SUBROUTINE CircuitsAndDynamics2DHarmonic( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE CircuitsMod
  USE CircMatInitMod
  USE CircEqWriteMod
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: Found, First=.TRUE.

  TYPE(Solver_t), POINTER :: Asolver => Null(), RotMSolver => Null()

  INTEGER :: p,q,bid,ind,nn,nd,nm,ind0
  TYPE(ValueList_t), POINTER :: Params,BC
  INTEGER, POINTER :: PS(:)
  REAL(KIND=dp)::vind,A,b,sigma,SumResistance
  INTEGER :: CompInd
  TYPE(Variable_t), POINTER :: LagrangeVar, AngVar, VeloVar
  REAL(KIND=dp), TARGET :: torq,imom=0,ang=0._dp,velo=0._dp,scale,sclA,sclB

  COMPLEX(KIND=dp), ALLOCATABLE :: Tcoef(:,:,:)
  REAL(KIND=dp), ALLOCATABLE :: RotM(:,:,:)
  TYPE(Mesh_t), POINTER :: Mesh  
  TYPE(Valuelist_t), POINTER :: Material, CompParams
  TYPE(Bodyarray_t), POINTER :: CircuitVariableBody

  integer :: slen,i,j,k,m,n,istat,BodyId, &
             ColId,RowId,ImRowId,jj, &
             ComponentId, circ_comp_count
  CHARACTER(LEN=MAX_NAME_LEN) :: name,cmd,CoilType,dofnumber
  REAL(KIND=dp) :: BodyY

  LOGICAL  :: owner, STAT
  REAL(KIND=dp) :: Omega
  COMPLEX(KIND=dp) :: cmplx_value
  COMPLEX(KIND=dp) :: i_multiplier

  TYPE(CMPLXComponent_t), POINTER :: Comp
  TYPE(CMPLXCircuitVariable_t), POINTER :: Cvar
  TYPE(Matrix_t), POINTER :: CM
  INTEGER, POINTER :: n_Circuits => Null(), circuit_tot_n => Null()
  TYPE(CMPLXCircuit_t), POINTER :: Circuits(:)
    
  LOGICAL, ALLOCATABLE :: Adirichlet(:)

  TYPE(Element_t), POINTER :: e, e_p

  TYPE(Variable_t), POINTER :: RotMvar

  COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)

  LOGICAL :: CSymmetry
  
  SAVE CSymmetry, RotMvar, Adirichlet, nm, Omega
!------------------------------------------------------------------------------

  IF (First) THEN
    First = .FALSE.
    
    CALL AddComponentsToBodyLists()
    
    CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )

    Mesh => Model % Mesh
    N = Mesh % MaxElementDOFs

    ALLOCATE( Model%Circuit_tot_n, Model%n_Circuits, STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'CircuitsAndDynamics2DHarmonic', 'Memory allocation error.' )
    END IF

    n_Circuits => Model%n_Circuits
    Circuit_tot_n => Model%Circuit_tot_n
    Circuit_tot_n = 0
  
    Omega = GetAngularFrequency()

    Model % ASolver => FindSolverWithKey('Export Lagrange Multiplier', 26)
    ASolver => Model % ASolver
    
    ! Initialize circuit matrices:
    ! ----------------------------
    CALL AllocateCircuitsList()
    Circuits => Model%CMPLXCircuits
    
    DO p=1,n_Circuits
      
      n = GetNofCircVariables(p)
      CALL AllocateCircuit(p)
      
      Circuits(p) % n_comp = CountNofCircComponents(p, n)
      ALLOCATE(Circuits(p) % Components(Circuits(p) % n_comp))
      
      CALL ReadCircuitVariables(p)
      CALL ReadComponents(p)
      CALL AddComponentValuesToLists(p)  ! Lists are used to communicate values to other solvers at the moment...
      CALL AddBareCircuitVariables(p)   ! these don't belong to any components
      CALL ReadCoefficientMatrices(p)
      CALL ReadPermutationVector(p)
      CALL ReadCircuitSources(p)
      CALL WriteCoeffVectorsForCircVariables(p)
    
    END DO

    ! Create CRS matrix strucures for the circuit equations:
    ! ------------------------------------------------------
    CALL Circuits_MatrixInit()
  END IF

  Circuits => Model%CMPLXCircuits
  n_Circuits => Model%n_Circuits
  Circuit_tot_n => Model%Circuit_tot_n
  CM=>Model%CircuitMatrix
  
  ! Initialialize Circuit matrix:
  ! -----------------------------
  PS => Asolver % Variable % Perm
  nm =  Asolver % Matrix % NumberOfRows
  IF(.NOT.ASSOCIATED(CM)) RETURN

  CM % RHS = 0._dp
  IF(ASSOCIATED(CM % Values)) CM % Values = 0._dp

  ! Write Circuit equations:
  ! ------------------------
  DO p = 1,n_Circuits
    CALL AddBasicCircuitEquations(p)
    CALL AddComponentEquationsAndCouplings(p, Model % Mesh % MaxElementDOFs, CSymmetry)
  END DO
  Asolver %  Matrix % AddMatrix => CM

  IF(ASSOCIATED(CM)) THEN
    IF(  CM % Format == MATRIX_LIST ) CALL List_toCRSMatrix(CM)
    IF(CM % NumberOfRows<=0)  THEN
      CALL FreeMatrix(CM)
      Asolver % Matrix % AddMatrix => Null()
    END IF
  ELSE
     ASolver % Matrix % AddMatrix => Null()
  END IF

!------------------------------------------------------------------------------
END SUBROUTINE CircuitsAndDynamics2DHarmonic
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE CircuitsOutput(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
   
   USE DefUtils
   USE CircuitsMod

   IMPLICIT NONE
   
   TYPE(Model_t) :: Model
   TYPE(Solver_t) :: Solver
   REAL(KIND=dp) :: dt
   LOGICAL :: Transient

   REAL(KIND=dp), POINTER :: LagrangeValues(:)
   TYPE(Variable_t), POINTER :: LagrangeVar
   REAL(KIND=dp), ALLOCATABLE  :: ip(:), ipt(:)
   INTEGER :: nm
   
   TYPE(Solver_t), POINTER :: ASolver

   CHARACTER(LEN=MAX_NAME_LEN) :: dofnumber
   INTEGER :: i,p,jj,j
   LOGICAL :: Found
   TYPE(CMPLXCircuitVariable_t), POINTER :: CVar

   TYPE(Matrix_t), POINTER :: CM    
   INTEGER, POINTER :: n_Circuits => Null(), circuit_tot_n => Null()
   TYPE(CMPLXCircuit_t), POINTER :: Circuits(:)

   Circuit_tot_n => Model%Circuit_tot_n
   n_Circuits => Model%n_Circuits
   CM => Model%CircuitMatrix
   Circuits => Model%CMPLXCircuits
   
    ! Look for the solver we attach the circuit equations to:
    ! -------------------------------------------------------
    Found = .False.
    DO i=1,Model % NumberOfSolvers
      Asolver => Model % Solvers(i)
      IF(ListCheckPresent(Asolver % Values,'Export Lagrange Multiplier')) THEN 
        Found = .True. 
        EXIT
      END IF
    END DO
    
    nm =  Asolver % Matrix % NumberOfRows

    ! Circuit variable values from previous timestep:
    ! -----------------------------------------------
    ALLOCATE(ip(circuit_tot_n), ipt(circuit_tot_n))
    ip = 0._dp
    ipt = 0._dp
    LagrangeVar => VariableGet( Solver % Mesh % Variables,'LagrangeMultiplier')
    IF(ASSOCIATED(LagrangeVar)) THEN
      IF(ParEnv % PEs>1) THEN
        DO i=1,SIZE(LagrangeVar % Values)
          IF( CM % RowOwner(nm+i)==Parenv%myPE) ipt(i) = LagrangeVar%Values(i)
        END DO
        CALL MPI_ALLREDUCE(ipt,ip,circuit_tot_n, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, ASolver % Matrix % Comm, j)
      ELSE
        ip(1:SIZE(LagRangeVar % Values)) = LagrangeVar % Values
      END IF
    END IF
     
    ! Export circuit & dynamic variables for "SaveScalars":
    ! -----------------------------------------------------

    CALL ListAddConstReal(GetSimulation(),'res: time', GetTime())
    
    DO p=1,n_Circuits
      DO i=1,Circuits(p) % n
        Cvar => Circuits(p) % CircuitVariables(i)
        
        CALL ListAddConstReal( GetSimulation(), 'res: '//TRIM(Circuits(p) % names(i))//' re', ip(Cvar % ValueId))
        CALL ListAddConstReal( GetSimulation(), 'res: '//TRIM(Circuits(p) % names(i))//' im', ip(Cvar % ImValueId))

        IF (Cvar % pdofs /= 0 ) THEN
          DO jj = 1, Cvar % pdofs
            write (dofnumber, "(I2)") jj
            CALL ListAddConstReal( GetSimulation(), 'res: '//TRIM(Circuits(p) % names(i))&
                                   //'re dof '//TRIM(dofnumber), ip(Cvar % ValueId + ReIndex(jj)))
            CALL ListAddConstReal( GetSimulation(), 'res: '//TRIM(Circuits(p) % names(i))&
                                   //'im dof '//TRIM(dofnumber), ip(Cvar % ValueId + ImIndex(jj)))
          END DO
        END IF

      END DO
    END DO
END SUBROUTINE CircuitsOutput
