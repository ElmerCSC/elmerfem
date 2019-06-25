
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
! *  Authors:   Eelis Takala(Trafotek Oy) and Juha Ruokolainen(CSC)
! *  Emails:    eelis.takala@trafotek.fi and Juha.Ruokolainen@csc.fi
! *  Web:       http://www.trafotek.fi and http://www.csc.fi/elmer
! *  Addresses: Trafotek Oy
! *             Kaarinantie 700
! *             Turku
! *
! *             and
! *
! *             CSC - IT Center for Science Ltd.
! *             Keilaranta 14
! *             02101 Espoo, Finland 
! *
! *  Original Date: October 2015
! *
! *****************************************************************************/
 
MODULE CircuitUtils

    USE DefUtils
    IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
  FUNCTION GetCircuitModelDepth() RESULT (Depth)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    TYPE(Valuelist_t), POINTER :: simulation
    REAL(KIND=dp) :: depth
    LOGICAL :: Found, CSymmetry

    CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )

    simulation => GetSimulation()
    IF (.NOT. ASSOCIATED(simulation)) CALL Fatal ('GetCircuitModelDepth', 'Simulation not found!')
   
    depth = GetConstReal(simulation, 'Circuit Model Depth', Found)
    
    IF (.NOT. Found) THEN
      depth = 1._dp
      IF (CSymmetry) depth = 2._dp * pi
    END IF
       
!------------------------------------------------------------------------------
  END FUNCTION GetCircuitModelDepth
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION GetComponentParams(Element) RESULT (ComponentParams)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: i
    TYPE(Element_t), POINTER :: Element
    TYPE(Valuelist_t), POINTER :: BodyParams, ComponentParams
    LOGICAL :: Found
    
    BodyParams => GetBodyParams( Element )
    IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('GetCompParams', 'Body Parameters not found')
   
    i = GetInteger(BodyParams, 'Component', Found)
    
    IF (.NOT. Found) THEN
      ComponentParams => Null()
    ELSE
      ComponentParams => CurrentModel % Components(i) % Values
    END IF
   
!------------------------------------------------------------------------------
  END FUNCTION GetComponentParams
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE GetWPotential(Wbase)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: Wbase(:)

    CALL GetLocalSolution(Wbase,'W Potential')
    IF(.NOT. ANY(Wbase/=0._dp)) CALL GetLocalSolution(Wbase,'W')
!------------------------------------------------------------------------------
  END SUBROUTINE GetWPotential
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE AddComponentsToBodyLists()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    LOGICAL :: Found
    INTEGER :: i, j, k
    LOGICAL, SAVE :: Visited=.FALSE.
    ! Components and Bodies:
    ! ----------------------  
    INTEGER :: BodyId
    INTEGER, POINTER :: BodyAssociations(:) => Null()
    TYPE(Valuelist_t), POINTER :: BodyParams, ComponentParams
     
    IF (Visited) RETURN

    Visited = .TRUE.
    DO i = 1, SIZE(CurrentModel % Components)
      ComponentParams => CurrentModel % Components(i) % Values

      IF( ListGetLogical( ComponentParams,'Passive Component', Found ) ) CYCLE 
      
      IF (.NOT. ASSOCIATED(ComponentParams)) CALL Fatal ('AddComponentsToBodyList', &
                                                         'Component parameters not found!')
      BodyAssociations => ListGetIntegerArray(ComponentParams, 'Body', Found)

      IF (.NOT. Found) BodyAssociations => ListGetIntegerArray(ComponentParams, 'Master Bodies', Found)
      
      IF (.NOT. Found) CYCLE

      DO j = 1, SIZE(BodyAssociations)
        BodyId = BodyAssociations(j)
        BodyParams => CurrentModel % Bodies(BodyId) % Values
        IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('AddComponentsToBodyList', &
                                                      'Body parameters not found!')
        k = GetInteger(BodyParams, 'Component', Found)
        IF (Found) CALL Fatal ('AddComponentsToBodyList', &
                              'Body '//TRIM(i2s(BodyId))//' associated to two components!')
        CALL listAddInteger(BodyParams, 'Component', i)
        BodyParams => Null()
      END DO
    END DO

    DO i = 1, SIZE(CurrentModel % Bodies)
      BodyParams => CurrentModel % Bodies(i) % Values
      IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('AddComponentsToBodyList', &
                                                   'Body parameters not found!')
      j = GetInteger(BodyParams, 'Component', Found)
      IF (.NOT. Found) CYCLE

      WRITE(Message,'(A,I2,A,I2)') 'Body',i,' associated to Component', j
      CALL Info('AddComponentsToBodyList',Message,Level=3)
      BodyParams => Null()
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE AddComponentsToBodyLists
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION GetComponentBodyIds(Id) RESULT (BodyIds)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    LOGICAL :: Found
    INTEGER :: Id
    INTEGER, POINTER :: BodyIds(:)
    TYPE(Valuelist_t), POINTER :: ComponentParams
    
    ComponentParams => CurrentModel % Components(Id) % Values
    
    IF (.NOT. ASSOCIATED(ComponentParams)) CALL Fatal ('GetComponentBodyIds', &
                                                         'Component parameters not found!')
    BodyIds => ListGetIntegerArray(ComponentParams, 'Body', Found)
    IF (.NOT. Found) BodyIds => ListGetIntegerArray(ComponentParams, 'Master Bodies', Found)
    IF (.NOT. Found) BodyIds => Null()
    
!------------------------------------------------------------------------------
  END FUNCTION GetComponentBodyIds
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION GetComponentHomogenizationBodyIds(Id) RESULT (BodyIds)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    LOGICAL :: Found
    INTEGER :: Id
    INTEGER, POINTER :: BodyIds(:)
    TYPE(Valuelist_t), POINTER :: ComponentParams
    
    ComponentParams => CurrentModel % Components(Id) % Values
    
    IF (.NOT. ASSOCIATED(ComponentParams)) CALL Fatal ('GetComponentHomogenizationBodyIds', &
                                                         'Component parameters not found!')
    BodyIds => ListGetIntegerArray(ComponentParams, 'Homogenization Parameters Body', Found)
    IF (.NOT. Found) BodyIds => GetComponentBodyIds(Id)

!------------------------------------------------------------------------------
  END FUNCTION GetComponentHomogenizationBodyIds
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION FindSolverWithKey(key, char_len) RESULT (Solver)
!------------------------------------------------------------------------------
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

END MODULE CircuitUtils


MODULE CircuitsMod

  USE DefUtils
  IMPLICIT NONE

CONTAINS 

!------------------------------------------------------------------------------
  SUBROUTINE AllocateCircuitsList()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: slen,n_Circuits
    CHARACTER(LEN=MAX_NAME_LEN) :: cmd, name

    ! Read Circuit definitions from MATC:
    ! ----------------------------------
    cmd = "Circuits"
    slen = LEN_TRIM(cmd)
    CALL Matc( cmd, name, slen )
    READ(name(1:slen), *) n_Circuits
    
    CurrentModel%n_Circuits = n_Circuits
    
    ALLOCATE( CurrentModel%Circuits(n_Circuits) )

!------------------------------------------------------------------------------
  END SUBROUTINE AllocateCircuitsList
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION CountNofCircVarsOfType(CId, Var_type) RESULT (nofc)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nofc, char_len, slen, CId, i
    CHARACTER(LEN=MAX_NAME_LEN) :: Var_type
    CHARACTER(LEN=MAX_NAME_LEN) :: name,cmd
    TYPE(Circuit_t), POINTER :: Circuit
    
    Circuit => CurrentModel%Circuits(CId)
    
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
    IMPLICIT NONE
    INTEGER :: nofc, nofvar, slen, CId, i, j, CompId
    INTEGER :: ComponentIDs(nofvar)
    CHARACTER(LEN=MAX_NAME_LEN) :: name,cmd
    TYPE(Circuit_t), POINTER :: Circuit

    nofc = 0
    ComponentIDs = -1
    
    Circuit => CurrentModel%Circuits(CId)
    
   
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
    IMPLICIT NONE
    INTEGER :: slen, ComponentId,i,j,CId, CompInd
    CHARACTER(LEN=MAX_NAME_LEN) :: cmd, name
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(CircuitVariable_t), POINTER :: CVar

    Circuit => CurrentModel%Circuits(CId)
    
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
    IMPLICIT NONE
    INTEGER :: CId, n, slen 
    CHARACTER(LEN=MAX_NAME_LEN) :: cmd, name
    TYPE(Circuit_t), POINTER :: Circuit

    Circuit => CurrentModel%Circuits(CId)
    
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
    IMPLICIT NONE
    INTEGER :: slen,CId,n
    CHARACTER(LEN=MAX_NAME_LEN) :: cmd, name
    TYPE(Circuit_t), POINTER :: Circuit

    Circuit => CurrentModel%Circuits(CId)
    
    n = Circuit % n
    
    ALLOCATE( Circuit % ComponentIds(n), Circuit % names(n) )
    ALLOCATE( Circuit % source(n) )
    ALLOCATE( Circuit % CircuitVariables(n), Circuit % Perm(n) )
    ALLOCATE( Circuit % A(n,n), Circuit % B(n,n), &
              Circuit % Mre(n,n), Circuit % Mim(n,n)  )
    Circuit % ComponentIds = 0
    Circuit % names = ' '
    Circuit % A = 0._dp
    Circuit % B = 0._dp
    Circuit % Mre = 0._dp
    Circuit % Mim = 0._dp

!------------------------------------------------------------------------------
  END SUBROUTINE AllocateCircuit
!------------------------------------------------------------------------------

!-------------------------------------------------------------------
 SUBROUTINE SetBoundaryAreasToValueLists()
!-------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Valuelist_t), POINTER :: BC
    REAL(KIND=dp) :: BoundaryAreas(CurrentModel % NumberOFBCs)
    INTEGER :: Active, t, i, BCid, n
    LOGICAL :: Found

    BoundaryAreas = 0._dp
    Mesh => CurrentModel % Mesh

    DO i=1, CurrentModel % NumberOfBcs
       BC => CurrentModel % BCs(i) % Values
       IF (.NOT. ASSOCIATED(BC) ) CALL Fatal('SetBoundaryAreasToValueLists', 'Boundary not found!')
       CALL ListAddInteger(BC, 'Boundary Id', i)
    END DO
    
    Active = GetNOFBoundaryElements()
    DO t=1,Active
       Element => GetBoundaryElement(t)
!       IF (.NOT. ActiveBoundaryElement()) CYCLE
       
       BC=>GetBC()
       IF (.NOT. ASSOCIATED(BC) ) CYCLE
     
       BCid = GetInteger(BC, 'Boundary Id', Found)
       n = GetElementNOFNodes() 
       BoundaryAreas(BCid) = BoundaryAreas(BCid) + ElementAreaNoAxisTreatment(Mesh, Element, n) 
    END DO

    DO i=1, CurrentModel % NumberOfBcs
       BC => CurrentModel % BCs(i) % Values
       IF (.NOT. ASSOCIATED(BC) ) CALL Fatal('ComputeCoilBoundaryAreas', 'Boundary not found!')
       BCid = GetInteger(BC, 'Boundary Id', Found)
       BoundaryAreas(BCid) = ParallelReduction(BoundaryAreas(BCid))
       CALL ListAddConstReal(BC, 'Area', BoundaryAreas(BCid))
    END DO
    
!-------------------------------------------------------------------
 END SUBROUTINE SetBoundaryAreasToValueLists
!-------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE ReadComponents(CId)
!------------------------------------------------------------------------------
    USE CircuitUtils
    IMPLICIT NONE
    INTEGER :: CId, CompInd
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(Component_t), POINTER :: Comp
    TYPE(Valuelist_t), POINTER :: CompParams
    LOGICAL :: Found

    Circuit => CurrentModel%Circuits(CId)
    
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

      Comp % ElBoundaries => ListGetIntegerArray(CompParams, 'Electrode Boundaries', Found)
      
      SELECT CASE (Comp % CoilType) 
      CASE ('stranded')
        Comp % nofturns = GetConstReal(CompParams, 'Number of Turns', Found)
        IF (.NOT. Found) CALL Fatal('Circuits_Init','Number of Turns not found!')

        Comp % ElArea = GetConstReal(CompParams, 'Electrode Area', Found)
        IF (.NOT. Found) CALL ComputeElectrodeArea(Comp, CompParams)

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

        Comp % ElArea = GetConstReal(CompParams, 'Electrode Area', Found)
        IF (.NOT. Found) CALL ComputeElectrodeArea(Comp, CompParams)

        Comp % N_j = Comp % nofturns / Comp % ElArea

      END SELECT
      CALL AddVariableToCircuit(Circuit, Comp % ivar, CId)
      CALL AddVariableToCircuit(Circuit, Comp % vvar, CId)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ReadComponents
!------------------------------------------------------------------------------

!-------------------------------------------------------------------
 SUBROUTINE ComputeElectrodeArea(Comp, CompParams)
!-------------------------------------------------------------------
  USE ElementUtils
  IMPLICIT NONE
  TYPE(Component_t), POINTER :: Comp
  TYPE(ValueList_t), POINTER :: CompParams, BC
  TYPE(Element_t), POINTER :: Element
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER :: t, n, BCid
  LOGICAL :: Found
  
  Mesh => CurrentModel % Mesh
  Comp % ElArea = 0._dp

  IF (CoordinateSystemDimension() == 2) THEN
    DO t=1,GetNOFActive()
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes() 
      IF (ElAssocToComp(Element, Comp)) THEN
        Comp % ElArea = Comp % ElArea + ElementAreaNoAxisTreatment(Mesh, Element, n) 
      END IF
    END DO
    Comp % ElArea = ParallelReduction(Comp % ElArea)
  ELSE
    IF (.NOT. ASSOCIATED(Comp % ElBoundaries)) &
      CALL Fatal('ComputeElectrodeArea','Electrode Boundaries not found')

    BCid = Comp % ElBoundaries(1)
    BC => CurrentModel % BCs(BCid) % Values
    IF (.NOT. ASSOCIATED(BC) ) CALL Fatal('ComputeElectrodeArea', 'Boundary not found!')

    Comp % ElArea = GetConstReal(BC, 'Area', Found)
    IF (.NOT. Found) CALL Fatal('ComputeElectrodeArea', 'Area not found!')
    
  END IF
!-------------------------------------------------------------------
 END SUBROUTINE ComputeElectrodeArea
!-------------------------------------------------------------------

! This function is originally from ElementUtils. However, there is 
! some kind of treatment regarding axisymmetric cases which fails 
! here since we don't want that.
!------------------------------------------------------------------------------
   FUNCTION ElementAreaNoAxisTreatment( Mesh,Element,N ) RESULT(A)
!------------------------------------------------------------------------------
     TYPE(Mesh_t), POINTER :: Mesh
     INTEGER :: N
     TYPE(Element_t) :: Element
!------------------------------------------------------------------------------

     REAL(KIND=dp), TARGET :: NX(N),NY(N),NZ(N)

     REAL(KIND=dp) :: A

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     INTEGER :: N_Integ,t

     REAL(KIND=dp) :: Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3), &
              SqrtMetric,SqrtElementMetric

     TYPE(Nodes_t) :: Nodes

     LOGICAL :: stat

     REAL(KIND=dp) :: Basis(n),u,v,w,x,y,z
     REAL(KIND=dp) :: dBasisdx(n,3)

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
!------------------------------------------------------------------------------
 
     Nodes % x => NX
     Nodes % y => NY
     Nodes % z => NZ

     Nodes % x = Mesh % Nodes % x(Element % NodeIndexes)
     Nodes % y = Mesh % Nodes % y(Element % NodeIndexes)
     Nodes % z = Mesh % Nodes % z(Element % NodeIndexes)

     IntegStuff = GaussPoints( element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ  = IntegStuff % n
!
!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
!
       A = 0.0
       DO t=1,N_Integ
!
!        Integration stuff
!
         u = U_Integ(t)
         v = V_Integ(t)
         w = W_Integ(t)
!
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                    Basis,dBasisdx )
!------------------------------------------------------------------------------
!        Coordinatesystem dependent info
!------------------------------------------------------------------------------
           A =  A + SqrtElementMetric * S_Integ(t)
       END DO
!------------------------------------------------------------------------------
   END FUNCTION ElementAreaNoAxisTreatment
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE AddVariableToCircuit(Circuit, Variable, k)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t) :: Circuit
    TYPE(CircuitVariable_t) :: Variable
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

    IF (Circuit % Harmonic) THEN
      IF (Circuit % UsePerm) THEN
        Variable % valueId = Circuit % Perm(Circuit_tot_n + 1)
        Variable % ImValueId = Variable % valueId + 1
      ELSE
        Variable % valueId = Circuit_tot_n + 1
        Variable % ImValueId = Circuit_tot_n + 2
      END IF
    
      Circuit_tot_n = Circuit_tot_n + 2*Variable % dofs
    ELSE
      IF (Circuit % UsePerm) THEN
        Variable % valueId = Circuit % Perm(Circuit_tot_n + 1)
      ELSE
        Variable % valueId = Circuit_tot_n + 1
      END IF
      
      Circuit_tot_n = Circuit_tot_n + Variable % dofs
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE AddVariableToCircuit
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE AddComponentValuesToLists(CId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(Component_t), POINTER :: Comp
    TYPE(Valuelist_t), POINTER :: CompParams
    INTEGER :: CId, CompInd
    
    Circuit => CurrentModel%Circuits(CId)
    
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
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(CircuitVariable_t), POINTER :: CVar
    INTEGER :: CId, i
    
    Circuit => CurrentModel%Circuits(CId)
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
    IMPLICIT NONE
    INTEGER :: CId,n
    TYPE(Circuit_t), POINTER :: Circuit

    Circuit => CurrentModel%Circuits(CId)
    n = Circuit % n

    ! Read in the coefficient matrices for the circuit equations:
    ! Ax' + Bx = source:
    ! ------------------------------------------------------------

    CALL matc_get_array('C.'//TRIM(i2s(CId))//'.A'//CHAR(0),Circuit % A,n,n)
    CALL matc_get_array('C.'//TRIM(i2s(CId))//'.B'//CHAR(0),Circuit % B,n,n)
    
    IF (Circuit % Harmonic) THEN
      ! Complex multiplier matrix is used for:
      ! B = times(M,B), where B times is the element-wise product
      ! ---------------------------------------------------------
      CALL matc_get_array('C.'//TRIM(i2s(CId))//'.Mre'//CHAR(0),Circuit % Mre,n,n)
      CALL matc_get_array('C.'//TRIM(i2s(CId))//'.Mim'//CHAR(0),Circuit % Mim,n,n)
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE ReadCoefficientMatrices
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReadPermutationVector(CId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: CId,n,slen,i
    CHARACTER(LEN=MAX_NAME_LEN) :: cmd, name
    TYPE(Circuit_t), POINTER :: Circuit

    Circuit => CurrentModel%Circuits(CId)
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
    IMPLICIT NONE
    INTEGER :: CId,n,slen,i
    CHARACTER(LEN=MAX_NAME_LEN) :: cmd, name
    TYPE(Circuit_t), POINTER :: Circuit

    Circuit => CurrentModel%Circuits(CId)
    n = Circuit % n
    DO i=1,n
      ! Names of the source functions, these functions should be found
      ! in the "Body Force 1" block of the .sif file.
      ! (nc: is for 'no check' e.g. don't abort if the MATC variable is not found!)
      ! ---------------------------------------------------------------------------
      cmd = 'nc:C.'//TRIM(i2s(CId))//'.source.'//TRIM(i2s(i))
      slen = LEN_TRIM(cmd)
      CALL Matc( cmd, name, slen )
      Circuit % Source(i) = name(1:slen)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ReadCircuitSources
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE WriteCoeffVectorsForCircVariables(CId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: CId,n,slen,i,j,RowId
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(CircuitVariable_t), POINTER :: Cvar
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
  
    Circuit => CurrentModel%Circuits(CId)
    n = Circuit % n

    DO i=1,n
      Cvar => Circuit % CircuitVariables(i)
      RowId = Cvar % ValueId

      ALLOCATE(Cvar % A(Circuit % n), &
               Cvar % B(Circuit % n), &
               Cvar % Mre(Circuit % n), &
               Cvar % Mim(Circuit % n), &
               Cvar % SourceRe(Circuit % n), &
               Cvar % SourceIm(Circuit % n))
      Cvar % A = 0._dp
      Cvar % B = 0._dp
      Cvar % Mre = 0._dp
      Cvar % Mim = 0._dp
      Cvar % SourceRe = 0._dp
      Cvar % SourceIm = 0._dp

      DO j=1,Circuit % n
        IF (Circuit % A(i,j)/=0) Cvar % A(j) = Circuit % A(i,j)
        IF (Circuit % B(i,j)/=0) Cvar % B(j) = Circuit % B(i,j)
        IF (Circuit % Mre(i,j)/=0 .OR. Circuit % Mim(i,j)/=0) THEN
          Cvar % Mre(j) = Circuit % Mre(i,j) 
          Cvar % Mim(j) = Circuit % Mim(i,j)
        END IF
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
     IMPLICIT NONE
     TYPE(Component_t), POINTER :: Component
     TYPE(Element_t), POINTER :: Element
     LOGICAL :: T
     T = IdInList(Element % BodyId, Component % BodyIds)
!------------------------------------------------------------------------------
   END FUNCTION ElAssocToComp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   FUNCTION ElAssocToCvar(Element, Cvar) RESULT (T)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     TYPE(CircuitVariable_t), POINTER :: Cvar
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
  FUNCTION AddIndex(Ind, Harmonic)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    Integer :: Ind, AddIndex
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm
    
    IF (.NOT. PRESENT(Harmonic)) THEN
      harm = CurrentModel % HarmonicCircuits
    ELSE
      harm = Harmonic
    END IF
 
    IF (harm) THEN
      AddIndex = 2 * Ind
    ELSE
      AddIndex = Ind
    END IF
!------------------------------------------------------------------------------
  END FUNCTION AddIndex 
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION AddImIndex(Ind)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: Ind
    Integer :: AddImIndex
    IF ( .NOT. CurrentModel % HarmonicCircuits ) CALL Fatal ('AddImIndex','Model is not of harmonic type!')
    
    AddImIndex = 2 * Ind + 1
!------------------------------------------------------------------------------
  END FUNCTION AddImIndex 
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION ReIndex(Ind, Harmonic)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: Ind, ReIndex
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm
    
    IF (.NOT. PRESENT(Harmonic)) THEN
      harm = CurrentModel % HarmonicCircuits
    ELSE
      harm = Harmonic
    END IF
 
    IF (harm) THEN
      ReIndex = 2 * Ind - 1
    ELSE
      ReIndex = Ind
    END IF
!------------------------------------------------------------------------------
  END FUNCTION ReIndex 
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION ImIndex(Ind)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    Integer :: Ind, ImIndex

    ImIndex = 2 * Ind
!------------------------------------------------------------------------------
  END FUNCTION ImIndex
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   FUNCTION HasSupport(Element, nn) RESULT(support)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nn, dim
    TYPE(Element_t) :: Element
    LOGICAL :: support, First=.TRUE.
    REAL(KIND=dp) :: wBase(nn)
    SAVE dim, First

    IF (First) THEN
      First = .FALSE.
      dim = CoordinateSystemDimension()
    END IF
    
    support = .TRUE. 
    IF (dim == 3) THEN
      support = .FALSE.
      CALL GetLocalSolution(Wbase,'W')
      IF ( ANY(Wbase .ne. 0d0) ) support = .TRUE.
    END IF
!------------------------------------------------------------------------------
   END FUNCTION HasSupport
!------------------------------------------------------------------------------

END MODULE CircuitsMod

MODULE CircMatInitMod

  USE CircuitsMod
  IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE SetCircuitsParallelInfo()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Matrix_t), POINTER :: CM
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(Solver_t), POINTER :: ASolver
    TYPE(Circuit_t), POINTER :: Circuits(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: i, nm, Circuit_tot_n, p, j, &
               cnt(Parenv % PEs), r_cnt(ParEnv % PEs), &
               RowId, nn, l, k, n_Circuits
    
    CM => CurrentModel%CircuitMatrix
    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('SetCircuitsParallelInfo','ASolver not found!')
    nm = ASolver % Matrix % NumberOfRows
    Circuit_tot_n = CurrentModel%Circuit_tot_n
    Circuits => CurrentModel % Circuits
    n_Circuits = CurrentModel % n_Circuits
    
#if 0
    IF(.NOT.ASSOCIATED(CM % ParallelInfo)) THEN
      ALLOCATE(CM % ParallelInfo)
      ALLOCATE(CM % ParallelInfo % NeighbourList(nm+Circuit_tot_n))
      DO i=1,nm+Circuit_tot_n
        CM % ParallelInfo % NeighbourList(i) % Neighbours => Null()
      END DO
    END IF
#endif

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
        IF(nn==0) THEN
          nn = Parenv % PEs
        ELSE
          IF( r_cnt(CVar % Owner+1)<=0 ) Nn=nn+1
        END IF
!         nn = Parenv % PEs

        IF (Circuits(p) % Harmonic) THEN
          DO j=1,Cvar % Dofs
#if 0
            IF(.NOT.ASSOCIATED(CM % ParallelInfo % NeighbourList(RowId+AddIndex(j-1))%Neighbours)) THEN
              ALLOCATE(CM % ParallelInfo % NeighbourList(RowId+AddIndex(j-1)) % Neighbours(nn))
              ALLOCATE(CM % ParallelInfo % NeighbourList(RowId+AddImIndex(j-1)) % Neighbours(nn))
            END IF
            CM % ParallelInfo % NeighbourList(RowId+AddIndex(j-1)) % Neighbours(1)   = CVar % Owner
            CM % ParallelInfo % NeighbourList(RowId+AddImIndex(j-1)) % Neighbours(1) = Cvar % Owner
            l = 1
            DO k=0,ParEnv % PEs-1
              IF(k==CVar % Owner) CYCLE
              IF(r_cnt(k+1)>0 .OR. nn==ParEnv % PEs) THEN
                l = l + 1
                CM % ParallelInfo % NeighbourList(RowId+AddIndex(j-1)) % Neighbours(l) = k
                CM % ParallelInfo % NeighbourList(RowId+AddImIndex(j-1)) % Neighbours(l) = k
              END IF
            END DO
#endif
            CM % RowOwner(RowId + AddIndex(j-1))   = Cvar % Owner
            CM % RowOwner(RowId + AddImIndex(j-1)) = Cvar % Owner
          END DO
        ELSE
          DO j=1,Cvar % Dofs
#if 0
            IF(.NOT.ASSOCIATED(CM % ParallelInfo % NeighbourList(RowId+j-1)%Neighbours)) THEN
              ALLOCATE(CM % ParallelInfo % NeighbourList(RowId+j-1) % Neighbours(nn))
            END IF
            CM % ParallelInfo % NeighbourList(RowId+j-1) % Neighbours(1) = CVar % Owner
            l = 1
            DO k=0,ParEnv % PEs-1
              IF(k==CVar % Owner) CYCLE
              IF(r_cnt(k+1)>0) THEN
                l = l + 1
                CM % ParallelInfo % NeighbourList(RowId+j-1) % Neighbours(l) = k
              END IF
            END DO
#endif
            CM % RowOwner(RowId+j-1) = Cvar % Owner
          END DO
        END IF
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
   SUBROUTINE CountMatElement(Rows, Cnts, RowId, dofs, Harmonic)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: Rows(:), Cnts(:)
    INTEGER :: RowId, dofs
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm
    
    IF (.NOT. PRESENT(Harmonic)) THEN
      harm = CurrentModel % HarmonicCircuits
    ELSE
      harm = Harmonic
    END IF
    
    IF (harm) THEN
      CALL CountCmplxMatElement(Rows, Cnts, RowId, dofs)
    ELSE
      Cnts(RowId) = Cnts(RowId) + dofs
    END IF

!------------------------------------------------------------------------------
   END SUBROUTINE CountMatElement
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
   SUBROUTINE CreateMatElement(Rows, Cols, Cnts, RowId, ColId, Harmonic)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: Rows(:), Cols(:), Cnts(:)
    INTEGER :: RowId, ColId
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm
    
    IF (.NOT. PRESENT(Harmonic)) THEN
      harm = CurrentModel % HarmonicCircuits
    ELSE
      harm = Harmonic
    END IF
    
    IF (harm) THEN
      CALL CreateCmplxMatElement(Rows, Cols, Cnts, RowId, ColId)
    ELSE
      Cols(Rows(RowId) + Cnts(RowId)) = ColId
      Cnts(RowId) = Cnts(RowId) + 1
    END IF

!------------------------------------------------------------------------------
   END SUBROUTINE CreateMatElement
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE CountBasicCircuitEquations(Rows, Cnts)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuits(:)
    TYPE(CircuitVariable_t), POINTER :: Cvar
    INTEGER :: i, j, p, nm, RowId, n_Circuits
    INTEGER, POINTER :: Rows(:), Cnts(:)
    
    Circuits => CurrentModel % Circuits
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
             CALL CountMatElement(Rows, Cnts, RowId, 1)
        END DO
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountBasicCircuitEquations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CreateBasicCircuitEquations(Rows, Cols, Cnts)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuits(:)
    TYPE(CircuitVariable_t), POINTER :: Cvar
    INTEGER :: i, j, p, nm, RowId, ColId, n_Circuits
    INTEGER, POINTER :: Rows(:), Cols(:), Cnts(:)
    
    Circuits => CurrentModel % Circuits
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
            CALL CreateMatElement(Rows, Cols, Cnts, RowId, ColId)
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
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuits(:)
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(Solver_t), POINTER :: ASolver
    TYPE(Element_t), POINTER :: Element
    TYPE(Component_t), POINTER :: Comp
    INTEGER :: i, j, p, nm, nn, nd, &
               RowId, ColId, n_Circuits, &
               CompInd, q
    INTEGER, POINTER :: Rows(:), Cnts(:)
    LOGICAL :: dofsdone
    LOGICAL*1 :: Done(:)
    
    Circuits => CurrentModel % Circuits
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
           CALL CountMatElement(Rows, Cnts, RowId, 1)
           CALL CountMatElement(Rows, Cnts, RowId, 1)
        CASE('massive')
           CALL CountMatElement(Rows, Cnts, RowId, 1)
           CALL CountMatElement(Rows, Cnts, RowId, 1)
        CASE('foil winding')
          ! V = V0 + V1*alpha + V2*alpha^2 + ...
          CALL CountMatElement(Rows, Cnts, RowId, Cvar % dofs)

          ! Circuit eqns for the pdofs:
          ! I(Vj) - I = 0
          ! ------------------------------------
          DO j=1, Cvar % pdofs
            CALL CountMatElement(Rows, Cnts, RowId + AddIndex(j), Cvar % dofs)
          END DO
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
              IF (.NOT. HasSupport(Element,nn)) CYCLE 
              CALL CountAndCreateMassive(Element,nn,nd,RowId,Cnts,Done,Rows)
           CASE('foil winding')
              IF (.NOT. HasSupport(Element,nn)) CYCLE 
              CALL CountAndCreateFoilWinding(Element,nn,nd,Comp,Cnts,Done,Rows)
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
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuits(:)
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(Solver_t), POINTER :: ASolver
    TYPE(Element_t), POINTER :: Element
    TYPE(Component_t), POINTER :: Comp
    INTEGER :: i, j, jj, p, nm, nn, nd, &
               VvarId, IvarId, n_Circuits, &
               CompInd, q
    INTEGER, POINTER :: Rows(:), Cols(:), Cnts(:)
    LOGICAL :: dofsdone
    LOGICAL*1 :: Done(:)
    
    Circuits => CurrentModel % Circuits
    n_Circuits = CurrentModel % n_Circuits
    Asolver => CurrentModel % Asolver
    nm = Asolver % Matrix % NumberOfRows

    DO p=1,n_Circuits
      DO CompInd = 1, Circuits(p) % n_comp
        Done = .FALSE.
        Comp => Circuits(p) % Components(CompInd)
        Cvar => Comp % vvar
        VvarId = Comp % vvar % ValueId + nm
        IvarId = Comp % ivar % ValueId + nm

        SELECT CASE (Comp % CoilType)
        CASE('stranded')
          CALL CreateMatElement(Rows, Cols, Cnts, VvarId, IvarId)
          CALL CreateMatElement(Rows, Cols, Cnts, VvarId, VvarId)
        CASE('massive')
          CALL CreateMatElement(Rows, Cols, Cnts, VvarId, IvarId)
          CALL CreateMatElement(Rows, Cols, Cnts, VvarId, VvarId)
        CASE('foil winding')
          DO j=0, Cvar % pdofs
            ! V = V0 + V1*alpha + V2*alpha^2 + ...
            CALL CreateMatElement(Rows, Cols, Cnts, VvarId, VvarId + AddIndex(j))
            IF (j/=0) THEN
              ! Circuit eqns for the pdofs:
              ! I(Vi) - I = 0
              ! ------------------------------------
              CALL CreateMatElement(Rows, Cols, Cnts, VvarId + AddIndex(j), IvarId)
              DO jj = 1, Cvar % pdofs
                  CALL CreateMatElement(Rows, Cols, Cnts, VvarId + AddIndex(j), VvarId + AddIndex(j))
              END DO
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
              CALL CountAndCreateStranded(Element,nn,nd,VvarId,Cnts,Done,Rows,Cols,IvarId)
            CASE('massive')
              IF (.NOT. HasSupport(Element,nn)) CYCLE 
              CALL CountAndCreateMassive(Element,nn,nd,VvarId,Cnts,Done,Rows,Cols=Cols)
           CASE('foil winding')
              IF (.NOT. HasSupport(Element,nn)) CYCLE   
              CALL CountAndCreateFoilWinding(Element,nn,nd,Comp,Cnts,Done,Rows,Cols=Cols)
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
   SUBROUTINE CountAndCreateStranded(Element,nn,nd,i,Cnts,Done,Rows,Cols,Jsind,Harmonic)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    INTEGER :: nn, nd, ncdofs1, ncdofs2, dim
    OPTIONAL :: Cols
    INTEGER :: Rows(:), Cols(:), Cnts(:)
    INTEGER :: p,i,j,Indexes(nd)
    INTEGER, OPTIONAL :: Jsind
    INTEGER, POINTER :: PS(:)
    LOGICAL*1 :: Done(:)
    LOGICAL :: First=.TRUE.
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm
    SAVE dim, First

    IF (First) THEN
      First = .FALSE.
      dim = CoordinateSystemDimension()
    END IF

    IF (.NOT. PRESENT(Harmonic)) THEN
      harm = CurrentModel % HarmonicCircuits
    ELSE
      harm = Harmonic
    END IF

    
    IF (.NOT. ASSOCIATED(CurrentModel % ASolver) ) CALL Fatal ('CountAndCreateStranded','ASolver not found!')
    PS => CurrentModel % Asolver % Variable % Perm
    nd = GetElementDOFs(Indexes,Element,CurrentModel % ASolver)
    IF(dim==2) THEN
      ncdofs1=1
      ncdofs2=nd
    ELSE IF(dim==3) THEN
      ncdofs1=nn
      ncdofs2=nd
    END IF
    
    DO p=ncdofs1,ncdofs2
      j = Indexes(p)
      IF(.NOT.Done(j)) THEN
        Done(j) = .TRUE.
        j = PS(j)
        IF (harm) j = ReIndex(j)
        IF(PRESENT(Cols)) THEN
          CALL CreateMatElement(Rows, Cols, Cnts, i, j, harm) 
          CALL CreateMatElement(Rows, Cols, Cnts, j, Jsind, harm)
!          CALL CreateMatElement(Rows, Cols, Cnts, j, Jsind)
        ELSE
          CALL CountMatElement(Rows, Cnts, i, 1, harm)
          CALL CountMatElement(Rows, Cnts, j, 1, harm)
!          CALL CountMatElement(Rows, Cnts, j, 1)
        END IF
      END IF
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountAndCreateStranded
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CountAndCreateMassive(Element,nn,nd,i,Cnts,Done,Rows,Cols,Harmonic)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    INTEGER :: nn, nd, ncdofs1, ncdofs2, dim
    OPTIONAL :: Cols
    INTEGER :: Rows(:), Cols(:), Cnts(:)
    INTEGER :: p,i,j,Indexes(nd)
    INTEGER, POINTER :: PS(:)
    LOGICAL*1 :: Done(:)
    LOGICAL :: First=.TRUE.
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm
    SAVE dim, First

    IF (First) THEN
      First = .FALSE.
      dim = CoordinateSystemDimension()
    END IF
    
    IF (.NOT. PRESENT(Harmonic)) THEN
      harm = CurrentModel % HarmonicCircuits
    ELSE
      harm = Harmonic
    END IF
    
    IF (.NOT. ASSOCIATED(CurrentModel % ASolver) ) CALL Fatal ('CountAndCreateMassive','ASolver not found!')
    PS => CurrentModel % Asolver % Variable % Perm
    nd = GetElementDOFs(Indexes,Element,CurrentModel % ASolver)
    IF(dim==2) THEN
      ncdofs1=1
      ncdofs2=nd
    ELSE IF(dim==3) THEN
      ncdofs1=nn
      ncdofs2=nd
    END IF
    DO p=ncdofs1,ncdofs2
      j = Indexes(p)
      IF(.NOT.Done(j)) THEN
        Done(j) = .TRUE.
        j = PS(j)
        IF (harm) j = ReIndex(j)
        IF(PRESENT(Cols)) THEN
          CALL CreateMatElement(Rows, Cols, Cnts, i, j, harm)
          CALL CreateMatElement(Rows, Cols, Cnts, j, i, harm)
        ELSE
          CALL CountMatElement(Rows, Cnts, i, 1, harm)
          CALL CountMatElement(Rows, Cnts, j, 1, harm)
        END IF
      END IF
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountAndCreateMassive
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
    SUBROUTINE CountAndCreateFoilWinding(Element,nn,nd,Comp,Cnts,Done,Rows,Cols,Harmonic)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    TYPE(Component_t), POINTER :: Comp
    INTEGER :: nn, nd, ncdofs, dim
    OPTIONAL :: Cols
    INTEGER :: Rows(:), Cols(:), Cnts(:)
    INTEGER :: Indexes(nd)
    INTEGER :: p,j,q,vpolord,vpolordtest,vpolord_tot,&
      dofId,dofIdtest,vvarId, nm
    LOGICAL :: dofsdone, First=.TRUE.
    INTEGER, POINTER :: PS(:)
    LOGICAL*1 :: Done(:)
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm
    SAVE dim, First

    IF (First) THEN
      First = .FALSE.
      dim = CoordinateSystemDimension()
    END IF
    
    IF (.NOT. PRESENT(Harmonic)) THEN
      harm = CurrentModel % HarmonicCircuits
    ELSE
      harm = Harmonic
    END IF
    
    IF (.NOT. ASSOCIATED(CurrentModel % ASolver) ) CALL Fatal ('CountAndCreateFoilWinding','ASolver not found!')
    PS => CurrentModel % Asolver % Variable % Perm
    nd = GetElementDOFs(Indexes,Element,CurrentModel % ASolver)
    nm = CurrentModel % ASolver % Matrix % NumberOfRows

    ncdofs=nd
    IF (dim == 3) ncdofs=nd-nn

    vvarId = Comp % vvar % ValueId
    vpolord_tot = Comp % vvar % pdofs - 1

    DO vpolordtest=0,vpolord_tot ! V'(alpha)
      dofIdtest = AddIndex(vpolordtest + 1) + vvarId
      DO vpolord = 0, vpolord_tot ! V(alpha)
        dofId = AddIndex(vpolord + 1) + vvarId
        IF (PRESENT(Cols)) THEN  
          CALL CreateMatElement(Rows, Cols, Cnts, dofIdtest+nm, dofId+nm, harm)
        ELSE
          CALL CountMatElement(Rows, Cnts, dofIdtest+nm, 1, harm)
        END IF
      END DO

      DO j=1,ncdofs
        q=j
        IF (dim == 3) q=q+nn
        IF (PRESENT(Cols)) THEN  
          q = PS(Indexes(q))
          IF (harm) q = ReIndex(q)
          CALL CreateMatElement(Rows, Cols, Cnts, dofIdtest+nm, q, harm)
        ELSE
          CALL CountMatElement(Rows, Cnts, dofIdtest+nm, 1, harm)
        END IF
      END DO
    END DO

    DO vpolord = 0, vpolord_tot ! V(alpha)
      dofId = AddIndex(vpolord + 1) + vvarId
      DO j=1,ncdofs
        q=j
        IF (dim == 3) q=q+nn
        q = PS(Indexes(q))
        IF (harm) q = ReIndex(q)
        IF (PRESENT(Cols)) THEN  
          CALL CreateMatElement(Rows, Cols, Cnts, q, dofId+nm, harm)
        ELSE
          CALL CountMatElement(Rows, Cnts, q, 1, harm)
        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountAndCreateFoilWinding
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE Circuits_MatrixInit()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Matrix_t), POINTER :: CM
    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:), Cnts(:)
    INTEGER, POINTER CONTIG :: Rows(:), Cols(:)
    INTEGER :: nm, Circuit_tot_n, n, i
    LOGICAL :: dofsdone
    LOGICAL*1, ALLOCATABLE :: Done(:)
    REAL(KIND=dp), POINTER CONTIG :: Values(:)

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Circuits_MatrixInit','ASolver not found!')
    Circuit_tot_n = CurrentModel%Circuit_tot_n
    
    ! Initialize Circuit matrix:
    ! -----------------------------
    PS => Asolver % Variable % Perm
    nm =  Asolver % Matrix % NumberOfRows

    CM => AllocateMatrix()
    CurrentModel%CircuitMatrix=>CM
    
    CM % Format = MATRIX_CRS
    Asolver % Matrix % AddMatrix => CM
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
      DEALLOCATE(Rows,Cnts,Done,CM); CM=>Null()
      Asolver %  Matrix % AddMatrix => CM
      CurrentModel%CircuitMatrix=>CM
      RETURN 
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


