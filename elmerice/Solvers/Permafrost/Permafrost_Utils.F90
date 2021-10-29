!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! ******************************************************************************
! *
! *  Authors: Thomas Zwinger, Denis Cohen, Juha Hartikainen
! *  Email:  thomas Zwinger [at] csc.fi 
! *  Web:     http://elmerice.elmerfem.org
! *  Address: CSC - Scientific Computing Ltd.  
! *               Keilaranta 14                    
! *               02101 Espoo, Finland             
! *                                                 
! *       Original Date:  January 2017  -               
! * 
! *****************************************************************************
!>  Auxiliary solvers for enhanced permafrost problem 
!---------------------------------------------------------------------------------------------


!==============================================================================
!>  initialization of IP variable to constant value
!> \ingroup Solvers
SUBROUTINE IPVariableInit(Model, Solver, Timestep, TransientSimulation )
  !==============================================================================

  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: Timestep
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: CurrentElement
  TYPE(Variable_t), POINTER :: IPVar
  REAL(KIND=dp), POINTER :: IPVarValue(:)
  TYPE(ValueHandle_t) :: InitialIPVar_h
  REAL(KIND=dp) :: InitValue, detJ
  INTEGER :: IPVarDOFs, I, t, ICId, N, istat
  INTEGER, POINTER :: IPVarPerm(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName="IPVariableInit"
  CHARACTER(LEN=MAX_NAME_LEN) :: IPVariableName
  TYPE(GaussIntegrationPoints_t), TARGET :: IP
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t) :: Nodes
  LOGICAL :: Visited = .FALSE., Found, ReadFromIC=.FALSE., stat
  REAL(KIND=dp), ALLOCATABLE :: Basis(:),dBasisdx(:,:)
  SAVE Visited

  IF (Visited) RETURN

  SolverParams => GetSolverParams()

  IPVariableName = ListGetString(SolverParams, &
       'IP Variable', Found )
  IF (.NOT.Found) THEN
    CALL FATAL(SolverName, ' "IP Variable" not found - you have to provide one')
  ELSE
    WRITE (Message,*) ' "IP Variable ": ', TRIM(IPVariableName),' found' 
    CALL INFO(SolverName, Message,Level=6)
  END IF
  IPVar => VariableGet( Solver % Mesh % Variables, IPVariableName,Found,UnfoundFatal=.TRUE. )
  
  IF ( ASSOCIATED( IPVar ) ) THEN
    IPVarPerm    => IPVar % Perm
    IPVarValue  => IPVar % Values
    IPVarDOFs = IPVar % DOFs
    InitValue = GetConstReal(SolverParams,TRIM(IPVariableName),Found)
    ReadFromIC = .NOT.(Found)
  ELSE
    CALL FATAL(SolverName, 'Could not find "IP Variable"')
  END IF
  
  IF (ReadFromIC) THEN
    CALL ListInitElementKeyword( InitialIPVar_h,'Initial Condition',TRIM(IPVariableName) )
    WRITE(Message,*) IPVariableName, ' from corresponding initial condition'
    N = 2 * MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
    ALLOCATE(Basis(N),dBasisdx(N,3),stat=istat)
  ELSE
    WRITE(Message,*) IPVariableName, ' to constant',  InitValue
    IPVarValue = InitValue
  END IF
  
  CALL INFO(SolverName, '-----------------------------------', Level=4)
  CALL INFO(SolverName, 'Initializing ip variable           ', Level=4)
  CALL INFO(SolverName, Message, Level=4)
  CALL INFO(SolverName, '-----------------------------------', Level=4)

  Visited = .TRUE.
  
  IF (ReadFromIC) THEN
    DO i = 1,  Solver % NumberOFActiveElements
      Element => GetActiveElement(i)

      IP = GaussPointsAdapt( Element )
      IF( Element % ElementIndex == 1 ) THEN
        CALL INFO(SolverName,'Number of Gauss points for 1st element:'&
            //TRIM(I2S(IP % n)),Level=7)
      END IF
      
      CALL GetElementNodes( Nodes )
      ICid = GetICId( Element, Found )
      IF (.NOT.Found) CALL FATAL(SolverName,'Corresponding "Initial Condition" not found')

      DO t=1,IP % n
        IF (ReadFromIC) THEN
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t), detJ, Basis, dBasisdx )
          InitValue = ListGetElementReal(InitialIPVar_h, Basis, Element, Found, GaussPoint=t)
          IF (.NOT.Found) CALL FATAL(SolverName,"Initial value not found in IC")
        END IF
        IPVarValue((IPVarPerm(i)*IPVarDOFs) + t*IPVarDOFs) = InitValue
!        IPVarValue((IPVarPerm(i)*IPVarDOFs) + t*IPVarDOFs) = 1.0
      END DO
    END DO
    DEALLOCATE(Basis, dBasisdx)
  END IF

  CALL INFO(SolverName,"Itialisation Done",Level=6)
  
END SUBROUTINE IPVariableInit

  
!==============================================================================
!>  initialization of arbitrary scalar given by nodal file
!> \ingroup Solvers
!==============================================================================
SUBROUTINE NodalVariableInit(Model, Solver, Timestep, TransientSimulation )
  !==============================================================================

  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: Timestep
  LOGICAL :: TransientSimulation

  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: NodalVariable
  TYPE(ValueList_t), POINTER :: SolverParams,Material
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER, POINTER :: NodalVariablePerm(:)
  INTEGER,PARAMETER :: io=26
  INTEGER, ALLOCATABLE :: GlobalToLocalPerm(:)
  REAL(KIND=dp), POINTER :: NodalVariableValues(:)
  REAL(KIND=dp) :: InputField, InitValue, ValueOffset
  INTEGER :: DIM, i, j, CurrentNode, NumberOfNodes, MaxNumberOfGNodes, MinNumberOfGNodes,&
       OK,  counter, localGlobalRange
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName="NodalVariableInit"
  CHARACTER(LEN=MAX_NAME_LEN) :: NodalVariableName,NodalVariableFileName
  LOGICAL :: Visited = .FALSE., Found, Parallel, GotIt, FromFile=.FALSE.

  !SAVE Visited
  !,DIM,NumberOfRockRecords

  !------------------------------------------------------------------------------

  ! Execute solver only once at beginning
  !if (Visited) RETURN

  CALL INFO(SolverName, '-----------------------------------', Level=4)
  CALL INFO(SolverName, 'Initializing variable to reference ', Level=4)
  CALL INFO(SolverName, 'levels (either file or IC)         ', Level=4)
  CALL INFO(SolverName, '-----------------------------------', Level=4)


  DIM = CoordinateSystemDimension()
  Parallel = (ParEnv % PEs > 1)
  Mesh => GetMesh()

  ! Get variable to fill in
  SolverParams => GetSolverParams()

  NodalVariableName = ListGetString(SolverParams, &
       'Nodal Variable', GotIt )
  IF (.NOT.GotIt) THEN
    CALL FATAL(SolverName, ' "Nodal Variable" not found')
  END IF
  NodalVariable => VariableGet( Mesh % Variables, NodalVariableName,GotIt )
  IF (.NOT.GotIt) CALL FATAL(SolverName,"Variable not found")

  IF ( ASSOCIATED( NodalVariable ) ) THEN
    NodalVariablePerm    => NodalVariable % Perm
    NodalVariableValues  => NodalVariable % Values
    WRITE (Message,*) 'Reading variable ',TRIM(NodalVariableName)
    CALL INFO(SolverName,Message,Level=5)
  ELSE
    WRITE (Message,*) 'Could not find ',TRIM(NodalVariableName)
    CALL FATAL(SolverName, Message)
  END IF
  NodalVariableValues = 0.0_dp

  NodalVariableFileName = ListGetString(SolverParams, &
       'Nodal Variable File', FromFile )

  ValueOffset = GetConstReal(SolverParams,'Variable Offset',GotIt)
  IF (GotIt) THEN
    WRITE (Message,*) ' "Variable Offset" found and set to: ', ValueOffset
    CALL INFO(SolverName,Message,Level=5)
  END IF

  IF (.NOT.FromFile) THEN
    InitValue = GetConstReal(SolverParams,TRIM(NodalVariableName),Found)
    IF (.NOT.Found) THEN
      WRITE(Message,*) 'No entry for ',TRIM(NodalVariableName),&
           ' found in Solver section (IC version not implemented)'
      CALL FATAL(SolverName,Message)
    END IF
    NodalVariableValues = InitValue + ValueOffset
  ELSE
    NumberOfNodes = Mesh % NumberOfNodes

    IF (Parallel) THEN
      MaxNumberOfGNodes = MAXVAL(Mesh % ParallelInfo % GlobalDOFs)
      MinNumberOfGNodes = MINVAL(Mesh % ParallelInfo % GlobalDOFs)
      !localGlobalRange = MaxNumberOfGNodes - MinNumberOfGNodes
      IF (MaxNumberOfGNodes <= MinNumberOfGNodes) CALL FATAL(SolverName,"No nodes in parallel domain")
      ALLOCATE(GlobalToLocalPerm(MinNumberOfGNodes:MaxNumberOfGNodes), STAT=OK)
      IF (OK /= 0) CALL FATAL(SolverName,"Allocation error of GlobalToLocalPerm")
      GlobalToLocalPerm = 0
      DO I=1,NumberOfNodes
        GlobalToLocalPerm(Mesh % ParallelInfo % GlobalDOFs(I)) = I
      END DO
      PRINT *, TRIM(SolverName),": ParENV:",ParEnv % MyPE,".  Global Nodal Numbers from",&
           MinNumberOfGNodes,"to",MaxNumberOfGNodes
    ELSE
      MinNumberOfGNodes = 1
      MaxNumberOfGNodes = NumberOfNodes   
    END IF

    OPEN(unit = io, file = TRIM(NodalVariableFileName), status = 'old',action='read',iostat = ok)
    IF (ok /= 0) THEN
      WRITE(Message,'(A,A)') 'Unable to open file ',TRIM(NodalVariableFileName)
      CALL FATAL(TRIM(SolverName),TRIM(message))
    ELSE
      !------------------------------------------------------------------------------
      ! Read in the number of records ordered in global node-numbering
      ! in file (first line integer)
      !------------------------------------------------------------------------------
      DO J=1,MaxNumberOfGNodes ! all or in parallel up to max global index
        READ (io, *, END=70, IOSTAT=OK, ERR=80) counter, InputField
        IF (counter .NE. J) CALL FATAL(SolverName,'No concecutive numbering in file')
        IF (J < MinNumberOfGNodes) CYCLE

        IF (Parallel) THEN
          I = GlobalToLocalPerm(J)
          IF (I == 0) CYCLE ! point in range, but not in partition        
        ELSE
          I=J
        END IF
        !IF ((NodalVariablePerm(I)<1) .OR. (NodalVariablePerm(I)>NumberOfNodes)) THEN
        !  PRINT *, "NodalVariableInit:", ParEnv % myPE, "NodalVariablePerm(",I,")=",&
        !       NodalVariablePerm(I),">",NumberOfNodes
        !  CALL FATAL(SolverName,'No corresponding entry of target variable')
        !END IF
        NodalVariableValues(NodalVariablePerm(I)) = InputField + ValueOffset
        ! PRINT *,i,counter
      END DO
      !PRINT *, "END", i,counter
70    IF (J-1 .NE. MaxNumberOfGNodes) THEN
        WRITE (Message,*) 'Number of records ',i,' in file ',&
             TRIM(NodalVariableFileName),' does not match number of nodes ',&
             NumberOfNodes, ' in mesh'
        CALL FATAL(SolverName,Message)
      END IF
      CLOSE (io)
      IF (Parallel) &
           DEALLOCATE(GlobalToLocalPerm)
      RETURN
80    CALL FATAL(SolverName,"I/O error")
    END IF
  END IF
END SUBROUTINE NodalVariableInit

!==============================================================================
!> output of material parameter at element
!==============================================================================
SUBROUTINE PermafrostElmntOutput_init( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  USE PermaFrostMaterials

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
    !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  ! 1 eta0 (=etat)
  ! 2 etak
  ! 3 alphaL
  ! 4 alphaT  
  ! 6 cs0
  !( 8 Es0
  ! 8 nus0
  ! 9 ks0
  ! 10 Kgwh0)

  LOGICAL :: WriteToFile(7)=.FALSE., Found, WriteAll
  TYPE(ValueList_t), POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostElmntOutput'
  
  CALL INFO( SolverName, '---------------------------------------',Level=4 )
  CALL INFO( SolverName, ' Assignment element material variables ',Level=4 )
  CALL INFO( SolverName, '---------------------------------------',Level=4 )
  SolverParams => GetSolverParams()
  WriteAll=ListGetLogical(SolverParams,"Export all",Found)
  IF (WriteAll) THEN
    WriteToFile(1:7)=.TRUE.
  ELSE
    WriteToFile(1)=ListGetLogical(SolverParams,"Export eta0",Found)
    WriteToFile(2)=ListGetLogical(SolverParams,"Export etak",Found)
    WriteToFile(3)=ListGetLogical(SolverParams,"Export alphaL",Found)
    WriteToFile(4)=ListGetLogical(SolverParams,"Export alphaT",Found)
    WriteToFile(5)=ListGetLogical(SolverParams,"Export cs0",Found)
    WriteToFile(6)=ListGetLogical(SolverParams,"Export qexp",Found)
    WriteToFile(7)=ListGetLogical(SolverParams,"Export Kgwh0",Found)
  END IF
  IF(WriteToFile(1)) THEN 
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         "-elem -dofs 1 eta0")
    CALL INFO(SolverName,'Added eta0 as variable',Level=5)
  END IF
  IF(WriteToFile(2)) THEN 
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         "-elem -dofs 1 etak")
    CALL INFO(SolverName,'Added etak as variable',Level=5)
  END IF
  IF(WriteToFile(3)) THEN 
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         "-elem -dofs 1 alphaL")
    CALL INFO(SolverName,'Added alphaL as variable',Level=5)
  END IF
  IF(WriteToFile(4)) THEN 
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         "-elem -dofs 1 alphaT")
    CALL INFO(SolverName,'Added alphaT as variable',Level=5)
  END IF
  IF(WriteToFile(5)) THEN 
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         "-elem -dofs 1 cs0")
    CALL INFO(SolverName,'Added cs0 as variable',Level=5)    
  END IF
  IF(WriteToFile(6)) THEN 
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         "-elem -dofs 1 qexp")
    CALL INFO(SolverName,'Added qexp as variable',Level=5)    
  END IF
  IF(WriteToFile(7)) THEN 
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         "-elem -dofs 1 Kgwh0_11")
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         "-elem -dofs 1 Kgwh0_22")
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         "-elem -dofs 1 Kgwh0_12")
    IF ( CoordinateSystemDimension() ==3) THEN
      CALL ListAddString( SolverParams,&
           NextFreeKeyword('Exported Variable',SolverParams),&
           "-elem -dofs 1 Kgwh0_33")
      CALL ListAddString( SolverParams,&
           NextFreeKeyword('Exported Variable',SolverParams),&
           "-elem -dofs 1 Kgwh0_13")
      CALL ListAddString( SolverParams,&
           NextFreeKeyword('Exported Variable',SolverParams),&
           "-elem -dofs 1 Kgwh0_23")
    END IF
    CALL INFO(SolverName,'Added Kgwh0 as variable',Level=5)    
  END IF
  CALL INFO( SolverName, 'assignment done',Level=6 )
  CALL INFO( SolverName, '---------------------------------------',Level=6 )
END SUBROUTINE PermafrostElmntOutput_init
!!!!!!!!!!!!!
!==============================================================================
!> output of material parameter at element
!> \ingroup Solvers
SUBROUTINE PermafrostElmntOutput( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  USE PermaFrostMaterials

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation  
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  LOGICAL :: WriteToFile(12)=.FALSE., FirstTime=.TRUE., WriteAll, Found,&
       ElementWiseRockMaterial
  INTEGER :: Active, t, J, RockMaterialID, CurrentValue, NumberOfRockRecords,&
       NumberOfExportedValues=0, DIM
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Params, Material,SolverParams
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostElmntOutput'
  CHARACTER(LEN=MAX_NAME_LEN) :: ElementRockMaterialName, ElmntVarName
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  TYPE(Variable_t), POINTER :: ElmntVar
  INTEGER, POINTER :: ElmntVarPerm(:)
  REAL(KIND=dp), POINTER :: ElmntVarVal(:)
 
  SAVE FirstTime, WriteToFile, NumberOfRockRecords,&
       NumberOfExportedValues, DIM
       
  
  CALL INFO( SolverName, '---------------------------------------',Level=4 )
  CALL INFO( SolverName, ' Assignment element material variables ',Level=4 )
  CALL INFO( SolverName, '---------------------------------------',Level=4 )

  SolverParams => GetSolverParams()
  
  IF (FirstTime) THEN
    DIM=CoordinateSystemDimension()
    WriteAll=ListGetLogical(SolverParams,"Export all",Found)
    IF (WriteAll) THEN
      WriteToFile(1:7)=.TRUE.
      NumberOfExportedValues=9
      WriteToFile(8) = .TRUE. 
      WriteToFile(9)= .TRUE.
      IF (DIM==3) THEN
        WriteToFile(10)= .TRUE.
        WriteToFile(11)= .TRUE.
        WriteToFile(12)= .TRUE.
        NumberOfExportedValues=12
      END IF
    ELSE
      WriteToFile(1)=ListGetLogical(SolverParams,"Export eta0",Found)      
      WriteToFile(2)=ListGetLogical(SolverParams,"Export etak",Found)
      WriteToFile(3)=ListGetLogical(SolverParams,"Export alphaL",Found)
      WriteToFile(4)=ListGetLogical(SolverParams,"Export alphaT",Found)
      WriteToFile(5)=ListGetLogical(SolverParams,"Export cs0",Found)
      WriteToFile(6)=ListGetLogical(SolverParams,"Export qexp",Found)
      WriteToFile(7)=ListGetLogical(SolverParams,"Export Kgwh0",Found)
      NumberOfExportedValues = 7
       IF (WriteToFile(7)) THEN
        WriteToFile(8) = .TRUE.
        WriteToFile(9) = .TRUE. 
        NumberOfExportedValues=9
        IF (DIM==3) THEN
          WriteToFile(10)= .TRUE.
          WriteToFile(11)= .TRUE.
          WriteToFile(12)= .TRUE.
          NumberOfExportedValues=12

        END IF
      END IF
    END IF
    CALL INFO(SolverName,'Exporting '//TRIM(I2S(NumberOfExportedValues))//' values',Level=4)
  END IF
  
  Active = GetNOFActive()
  DO CurrentValue=1,NumberOfExportedValues
    IF (.NOT.WriteToFile(CurrentValue)) CYCLE
    SELECT CASE(CurrentValue)
    CASE(1)
      WRITE (ElmntVarName,'(A)') "eta0"
    CASE(2)
      WRITE (ElmntVarName,'(A)') "etak"
    CASE(3)
      WRITE (ElmntVarName,'(A)') "alphaL"      
    CASE(4)
      WRITE (ElmntVarName,'(A)') "alphaT"
    CASE(5)
      WRITE (ElmntVarName,'(A)') "cs0"
    CASE(6)
      WRITE (ElmntVarName,'(A)') "qexp"
    CASE(7)
      WRITE (ElmntVarName,'(A)') "Kgwh0_11"
    CASE(8)
      WRITE (ElmntVarName,'(A)') "Kgwh0_22"
    CASE(9)
      WRITE (ElmntVarName,'(A)') "Kgwh0_12"
    CASE(10)
      WRITE (ElmntVarName,'(A)') "Kgwh0_33"
    CASE(11)
      WRITE (ElmntVarName,'(A)') "Kgwh0_13"
    CASE(12)
      WRITE (ElmntVarName,'(A)') "Kgwh0_23"
    END SELECT
    WRITE (Message,*) 'Writing ', TRIM(ElmntVarName), ' as variable'
    CALL INFO(SolverName,Message,Level=5)
    ElmntVar => VariableGet( Model % Mesh % Variables, ElmntVarName)
    IF (.NOT.ASSOCIATED(ElmntVar)) THEN
      WRITE(Message,*) 'Variable ',TRIM(ElmntVarName),' is not associated'
      CALL FATAL(SolverName,Message)
    END IF
    ElmntVarPerm => ElmntVar % Perm 
    ElmntVarVal  => ElmntVar % Values    
    DO t=1,Active
      Element => GetActiveElement(t)
      Material => GetMaterial()
      IF (FirstTime) THEN
        ! check, whether we have globally or element-wise defined values of rock-material parameters
        ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
        IF (ElementWiseRockMaterial) THEN
          WRITE (Message,*) 'Found "Element Rock Material File"'
          CALL INFO(SolverName,Message,Level=5)
          CALL INFO(SolverName,'Using element-wise rock material definition',Level=5)
        END IF
        IF (ElementWiseRockMaterial) THEN
          ! read element-wise material parameter (GlobalRockMaterial will have one entry each element)
          NumberOfRockRecords = &
               ReadPermafrostElementRockMaterial(ElementRockMaterialName,Solver,DIM)         
        ELSE
          NumberOfRockRecords =  ReadPermafrostRockMaterial( Material )
        END IF
        
        IF (NumberOfRockRecords < 1) THEN
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost Rock Material read',Level=6)
          FirstTime = .FALSE.
        END IF
      END IF
      IF (ElementWiseRockMaterial) THEN
        RockMaterialID = t  ! each element has it's own set of parameters
      ELSE
        RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
      END IF
      !IF (CurrentValue >= 6) THEN
      !  PRINT *,"eta0:", GlobalRockMaterial % eta0(RockMaterialID),&
      !       "aL:", GlobalRockMaterial % alphaL(RockMaterialID),"Kgwh0:", GlobalRockMaterial % Kgwh0(1:3,1:3,RockMaterialID)
      !END IF
      SELECT CASE(CurrentValue)
      CASE(1)
        ElmntVarVal(ElmntVarPerm(t)) = GlobalRockMaterial % eta0(RockMaterialID)
        !PRINT *,"eta0: ", ParEnv % MyPE, ":", RockMaterialID, ElmntVarPerm(t), GlobalRockMaterial % eta0(RockMaterialID), GlobalRockMaterial % Kgwh0(1,1,RockMaterialID)
      CASE(2)
        ElmntVarVal(ElmntVarPerm(t)) = GlobalRockMaterial % etak(RockMaterialID)
      CASE(3)
        ElmntVarVal(ElmntVarPerm(t)) = GlobalRockMaterial % alphaL(RockMaterialID)     
      CASE(4)
        ElmntVarVal(ElmntVarPerm(t)) = GlobalRockMaterial % alphaT(RockMaterialID)
      CASE(5)
        ElmntVarVal(ElmntVarPerm(t)) = GlobalRockMaterial % cs0(RockMaterialID)        
      CASE(7)
        ElmntVarVal(ElmntVarPerm(t)) = GlobalRockMaterial % Kgwh0(1,1,RockMaterialID)
!!$        PRINT *, "Kgwh0(", RockMaterialID, ")", GlobalRockMaterial % Kgwh0(1,1,RockMaterialID), &
!!$             GlobalRockMaterial % Kgwh0(2,2,RockMaterialID), GlobalRockMaterial % Kgwh0(3,3,RockMaterialID), &
!!$             GlobalRockMaterial % etak(RockMaterialID), GlobalRockMaterial % aas(0,RockMaterialID),&
!!$             GlobalRockMaterial % aas(1,RockMaterialID), GlobalRockMaterial % aas(2,RockMaterialID)
      CASE(8)
        ElmntVarVal(ElmntVarPerm(t)) = GlobalRockMaterial % Kgwh0(2,2,RockMaterialID)
      CASE(9)
        ElmntVarVal(ElmntVarPerm(t)) = GlobalRockMaterial % Kgwh0(1,2,RockMaterialID)
      CASE(10)
        ElmntVarVal(ElmntVarPerm(t)) = GlobalRockMaterial % Kgwh0(3,3,RockMaterialID)
      CASE(11)
        ElmntVarVal(ElmntVarPerm(t)) = GlobalRockMaterial % Kgwh0(1,3,RockMaterialID)
      CASE(12)
        ElmntVarVal(ElmntVarPerm(t)) = GlobalRockMaterial % Kgwh0(2,3,RockMaterialID)
    END SELECT
    END DO
  END DO
END SUBROUTINE PermafrostElmntOutput
!==============================================================================
!> output of material parameter at IP
!> \ingroup Solvers
!==============================================================================
SUBROUTINE PermafrostIPOutput_init( Model,Solver,dt,TransientSimulation )
  !==============================================================================
  USE DefUtils
  USE PermaFrostMaterials

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation  
  !------------------------------------------------------------------------------
  ! Local variables
  TYPE(ValueList_t), POINTER :: SolverParams
  LOGICAL :: WriteIPVar(6)=.FALSE., Found
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostIPOutput_init'
  INTEGER :: DIM, kGpeDOFs, XikG0hyDOFs,KGTTDOFs
  
  CALL INFO( SolverName, '-------------------------------------',Level=4 )
  CALL INFO( SolverName, ' Permafrost IP Output initialization ',Level=4 )
  CALL INFO( SolverName, '-------------------------------------',Level=4 )
  
  DIM = CoordinateSystemDimension()
  SolverParams => GetSolverParams()
  WriteIPVar(1)=ListGetLogical(SolverParams,"Export rhogw",Found)      
  WriteIPVar(2)=ListGetLogical(SolverParams,"Export mugw",Found)
  WriteIPVar(3)=ListGetLogical(SolverParams,"Export kGpe",Found)
  WriteIPVar(4)=ListGetLogical(SolverParams,"Export KGTT", Found)
  WriteIPVar(5)=ListGetLogical(SolverParams,"Export XikG0hy",Found)
  WriteIPVar(6)=ListGetLogical(SolverParams,"Export CgwTT",Found)
  !WriteIPVar(6)=ListGetLogical(SolverParams,"Export ",Found)

  IF(WriteIPVar(1)) THEN 
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         "-IP -dofs 1 rhogw")
    CALL INFO(SolverName,'Added rhogw as variable',Level=5)
  END IF
  IF(WriteIPVar(2)) THEN 
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         "-IP -dofs 1 mugw")
    CALL INFO(SolverName,'Added mugw as variable',Level=5)
  END IF
  IF(WriteIPVar(3)) THEN
    IF (DIM == 1) THEN
      CALL ListAddString( SolverParams,&
           NextFreeKeyword('Exported Variable',SolverParams),&
           "-IP -dofs 1 kGpe")
      kGpeDOFs = 1
    ELSE IF (DIM == 2) THEN
      CALL ListAddString( SolverParams,&
           NextFreeKeyword('Exported Variable',SolverParams),&
           "-IP -dofs 4 kGpe")
       kGpeDOFs = 4
    ELSE
      CALL ListAddString( SolverParams,&
           NextFreeKeyword('Exported Variable',SolverParams),&
           "-IP -dofs 9 kGpe")
       kGpeDOFs = 9
     END IF
     WRITE (Message,*) 'Added kGpe as variable with ',kGpeDOFs,' DOFs'
     CALL INFO(SolverName,Message,Level=5)
  END IF
  IF(WriteIPVar(4)) THEN
    IF (DIM == 1) THEN
      CALL ListAddString( SolverParams,&
           NextFreeKeyword('Exported Variable',SolverParams),&
           "-IP -dofs 1 KGTT")
      KGTTDOFs = 1
    ELSE IF (DIM == 2) THEN
      CALL ListAddString( SolverParams,&
           NextFreeKeyword('Exported Variable',SolverParams),&
           "-IP -dofs 4 KGTT")
      KGTTDOFs = 4
    ELSE
      CALL ListAddString( SolverParams,&
           NextFreeKeyword('Exported Variable',SolverParams),&
           "-IP -dofs 9 KGTT")
      KGTTDOFs = 9
    END IF
    WRITE (Message,*) 'Added KGTT as variable with ',KGTTDOFs,' DOFs'
    CALL INFO(SolverName,Message,Level=5)
  END IF
  IF(WriteIPVar(5)) THEN
    IF (DIM == 1) THEN
      CALL ListAddString( SolverParams,&
           NextFreeKeyword('Exported Variable',SolverParams),&
           "-IP -dofs 1 XikG0hy")
      XikG0hyDOFs = 1
    ELSE IF (DIM == 2) THEN
      CALL ListAddString( SolverParams,&
           NextFreeKeyword('Exported Variable',SolverParams),&
           "-IP -dofs 4 XikG0hy")
      XikG0hyDOFs = 4
    ELSE
      CALL ListAddString( SolverParams,&
           NextFreeKeyword('Exported Variable',SolverParams),&
           "-IP -dofs 9 XikG0hy")
      XikG0hyDOFs = 9
    END IF
    WRITE (Message,*) 'Added XikG0hy as variable with ',XikG0hyDOFs,' DOFs'
    CALL INFO(SolverName,Message,Level=5)
  END IF
  IF(WriteIPVar(6)) THEN 
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         "-IP -dofs 1 CgwTT")
    CALL INFO(SolverName,'Added CgwTT as variable',Level=5)
  END IF
END SUBROUTINE PermafrostIPOutput_init
SUBROUTINE PermafrostIPOutput( Model,Solver,dt,TransientSimulation )
  !==============================================================================
  USE DefUtils
  USE PermaFrostMaterials

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation  
  !------------------------------------------------------------------------------  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: SolverParams, Material
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter,&
       maxiter, istat, DepthDOFs, kGpeDOFs, XikG0hyDOFs, KGTTDOFs
  INTEGER,PARAMETER :: io=23
  REAL(KIND=dp) :: Norm
 
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE., FluxOutput = .FALSE.,&
       ComputeDt=.FALSE.,ElementWiseRockMaterial, DepthExists=.FALSE.,&
       InitializeSteadyState=.FALSE.,ActiveMassMatrix=.TRUE.,&
       WriteIPVar(6)=.FALSE.
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostIPOutput'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName, GWfluxName, PhaseChangeModel,&
       ElementRockMaterialName,VarName, DepthName, XiAtIPName
  TYPE(ValueHandle_t) :: Load_h, Temperature_h, Pressure_h, Salinity_h, Porosity_h,&
       PressureVelo_h, SalinityVelo_h, Depth_h,&
       Vstar1_h, Vstar2_h, Vstar3_h

  SAVE DIM,FirstTime,AllocationsDone,FluxOutput,DepthName,XiAtIPName,&
       CurrentSoluteMaterial,CurrentSolventMaterial,NumberOfRockRecords,&
       ElementWiseRockMaterial,ComputeDt,DepthExists,&
       Load_h, Temperature_h, Pressure_h, Salinity_h, Porosity_h,&
       PressureVelo_h, SalinityVelo_h, Depth_h, &
       Vstar1_h, Vstar2_h, Vstar3_h, kGpeDOFs, KGTTDOFs, XikG0hyDOFs
 

  CALL INFO( SolverName, '-------------------------------------',Level=4 )
  CALL INFO( SolverName, ' Permafrost IP Output                ',Level=4 )
  CALL INFO( SolverName, '-------------------------------------',Level=4 )


  Active = GetNOFActive()

  IF (FirstTime) THEN
    DIM = CoordinateSystemDimension()
    kGpeDOFs = DIM*DIM
    KGTTDOFs = DIM*DIM
    XikG0hyDOFs = DIM*DIM
    ! Handles to other system variables
    CALL ListInitElementKeyword( Temperature_h, 'Material', 'Temperature Variable' )
    CALL ListInitElementKeyword( Pressure_h, 'Material', 'Pressure Variable' )
    CALL ListInitElementKeyword( Salinity_h, 'Material', 'Salinity Variable' )
    CALL ListInitElementKeyword( Porosity_h, 'Material', 'Porosity Variable' )

    ! Handles to time derivatives of system variables
    CALL ListInitElementKeyword( PressureVelo_h, 'Material', 'Pressure Velocity Variable' )
    CALL ListInitElementKeyword( SalinityVelo_h, 'Material', 'Salinity Velocity Variable' )

    ! Handle to Heat Source (possible description of heat source at elements/IP's) 
    CALL ListInitElementKeyword( Load_h,'Body Force','Heat Source' )

    ! Handles to advection velocities
    CALL ListInitElementKeyword( Vstar1_h,'Material','Convection Velocity 1')
    CALL ListInitElementKeyword( Vstar2_h,'Material','Convection Velocity 2')
    IF (DIM > 2) &
         CALL ListInitElementKeyword( Vstar3_h,'Material','Convection Velocity 3')
    
    ! Handles to other variables
    CALL ListInitElementKeyword( Depth_h, 'Material', 'Depth Variable' )

  END IF
  
  !CALL DefaultStart()
  SolverParams => GetSolverParams()
  WriteIPVar(1)=ListGetLogical(SolverParams,"Export rhogw",Found)      
  WriteIPVar(2)=ListGetLogical(SolverParams,"Export mugw",Found)
  WriteIPVar(3)=ListGetLogical(SolverParams,"Export kGpe",Found)
  WriteIPVar(4)=ListGetLogical(SolverParams,"Export KGTT",Found)
  WriteIPVar(5)=ListGetLogical(SolverParams,"Export XikG0hy",Found)
  WriteIPVar(6)=ListGetLogical(SolverParams,"Export CgwTT",Found)
  
  Active = GetNOFActive()
    DO t=1,Active
      Element => GetActiveElement(t)
      Material => GetMaterial()


      IF (FirstTime) THEN
        ! check, whether we have globally or element-wise defined values of rock-material parameters
        ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
        IF (ElementWiseRockMaterial) THEN
          WRITE (Message,*) 'Found "Element Rock Material File"'
          CALL INFO(SolverName,Message,Level=6)
          CALL INFO(SolverName,'Using element-wise rock material definition',Level=6)
        END IF
        IF (ElementWiseRockMaterial) THEN
          ! read element-wise material parameter (GlobalRockMaterial will have one entry each element)
          NumberOfRockRecords = &
               ReadPermafrostElementRockMaterial(ElementRockMaterialName,Solver,DIM)
        ELSE
          NumberOfRockRecords =  ReadPermafrostRockMaterial( Material )
        END IF

        IF (NumberOfRockRecords < 1) THEN
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost Rock Material read',Level=6)
          FirstTime = .FALSE.
        END IF
        CALL ReadPermafrostSoluteMaterial( Material,Model % Constants,CurrentSoluteMaterial )
        CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
      END IF

      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()

      PhaseChangeModel = ListGetString(Material, &
           'Permafrost Phase Change Model', Found )
      IF (Found) THEN
        WRITE (Message,'(A,A)') '"Permafrost Phase Change Model" set to ', TRIM(PhaseChangeModel)
        CALL INFO(SolverName,Message,Level=9)
      END IF

      CALL SetIPValues(  Element, Element % ElementIndex, Active, n, nd+nb, WriteIPVar,&
           CurrentSoluteMaterial, CurrentSolventMaterial,&
           NumberOfRockRecords, PhaseChangeModel,ElementWiseRockMaterial)
    END DO
CONTAINS
  SUBROUTINE SetIPValues(Element, ElementID, NoElements, n, nd,&
       WriteIPVar, CurrentSoluteMaterial, CurrentSolventMaterial,&
       NumberOfRockRecords, PhaseChangeModel, ElementWiseRockMaterial)
     IMPLICIT NONE
    !------------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n, nd, ElementID, NoElements, NumberOfRockRecords
    TYPE(Element_t), POINTER :: Element
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    LOGICAL, INTENT(IN) :: ElementWiseRockMaterial, WriteIPVar(6)
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: DepthAtIP,RefDepth,CGTTAtIP, CGTpAtIP ,CGTycAtIP! CgwTTAtIP,KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,XiEtaAtIP,&
         ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,hiAtIP,hwAtIP  ! function values needed for C's and KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP, bijAtIP(2,2), bijYcAtIP(2,2),&
         gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) ::  gradTAtIP(3),gradPAtIP(3),JgwDAtIP(3),KgwAtIP(3,3),KgwpTAtIP(3,3),MinKgw,&
         KgwppAtIP(3,3),fwAtIP,DtdAtIP(3,3)!  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0, eps, Gravity(3) ! constants read only once
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,csAtIP,cwAtIP,ciAtIP,ccAtIP ! material properties at IP
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,&
         PressureVeloAtIP,SalinityVeloAtIP,&
         StiffPQ, meanfactor, vstarAtIP(3), auxtensor(3,3)
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp), POINTER :: gWork(:,:)
    INTEGER :: i,j,k,t,p,q,IPPerm,DIM, RockMaterialID, FluxDOFs
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.,ConstVal=.FALSE.,&
         CryogenicSuction=.FALSE.,HydroGeo=.FALSE.,ComputeFlux=.TRUE.,&
         NoSalinity=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN) :: MaterialFileName
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='PermafrostIPOutput(SetIPValues)'
    TYPE(Variable_t), POINTER :: XiAtIPVar, rhogwAtIPVar, mugwAtIPVar, kGpeAtIPVar, &
         XikG0hyAtIPVar,KGTTAtIPVar,CgwTTAtIPVar
    INTEGER, POINTER :: XiAtIPPerm(:),GWfluxPerm(:),rhogwAtIPPerm(:),&
         mugwAtIPPerm(:), kGpeAtIPPerm(:), XikG0hyAtIPPerm(:),KGTTAtIPPerm(:),CgwTTAtIPPerm(:)
    REAL(KIND=dp), POINTER :: XiAtIP(:), FluxAtElem(:), rhogwAtIP(:),&
         mugwATIP(:), kGpeAtIP(:), XikG0hyAtIP(:), KGTTAtIP(:), CgwTTAtIP(:)

    !------------------------------------------------------------------------------
    SAVE Nodes, ConstantsRead, ConstVal,DIM, GasConstant, N0,DeltaT, T0, p0, eps, Gravity
    !------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()
    gradTAtIP = 0.0_dp
    gradPAtIP = 0.0_dp
    IF(.NOT.ConstantsRead) THEN
      ConstantsRead = &
           ReadPermafrostConstants(Model, FunctionName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)
    END IF

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    XiAtIPVar => VariableGet( Solver % Mesh % Variables, 'Xi')
    IF (.NOT.ASSOCIATED(XiAtIPVar)) THEN
      WRITE(Message,*) 'Variable Xi is not associated'
      CALL FATAL(SolverName,Message)
    END IF
    XiAtIPPerm => XiAtIPVar % Perm
    XiAtIp => XiAtIPVar % Values

    IF (WriteIPVar(1)) THEN
      rhogwAtIPVar => VariableGet( Solver % Mesh % Variables, 'rhogw')
      IF (.NOT.ASSOCIATED(rhogwAtIPVar)) THEN
        WRITE(Message,*) 'Variable "rhogw" is not associated'
        CALL FATAL(SolverName,Message)
      END IF
      rhogwAtIPPerm => rhogwAtIPVar % Perm
      rhogwAtIP => rhogwAtIPVar % Values
    END IF

    IF (WriteIPVar(2)) THEN
      mugwAtIPVar => VariableGet( Solver % Mesh % Variables, 'mugw')
      IF (.NOT.ASSOCIATED(mugwAtIPVar)) THEN
        WRITE(Message,*) 'Variable "mugw" is not associated'
        CALL FATAL(SolverName,Message)
      END IF
      mugwAtIPPerm => mugwAtIPVar % Perm
      mugwAtIP => mugwAtIPVar % Values
    END IF

    IF (WriteIPVar(3)) THEN
      kGpeAtIPVar => VariableGet( Solver % Mesh % Variables, 'kGpe')
      IF (.NOT.ASSOCIATED(kGpeAtIPVar)) THEN
        WRITE(Message,*) 'Variable "kGpe" is not associated'
        CALL FATAL(SolverName,Message)
      END IF
      kGpeAtIPPerm => kGpeAtIPVar % Perm
      kGpeAtIP => kGpeAtIPVar % Values
    END IF
    IF (WriteIPVar(4)) THEN
      KGTTAtIPVar => VariableGet( Solver % Mesh % Variables, 'KGTT')
      IF (.NOT.ASSOCIATED(KGTTAtIPVar)) THEN
        WRITE(Message,*) 'Variable "KGTT" is not associated'
        CALL FATAL(SolverName,Message)
      END IF
      KGTTAtIPPerm => KGTTAtIPVar % Perm
      KGTTAtIP => KGTTAtIPVar % Values
    END IF
    IF (WriteIPVar(5)) THEN
      XikG0hyAtIPVar => VariableGet( Solver % Mesh % Variables, 'XikG0hy')
      IF (.NOT.ASSOCIATED(XikG0hyAtIPVar)) THEN
        WRITE(Message,*) 'Variable "XikG0hy" is not associated'
        CALL FATAL(SolverName,Message)
      END IF
      XikG0hyAtIPPerm => XikG0hyAtIPVar % Perm
      XikG0hyAtIP => XikG0hyAtIPVar % Values
    END IF
    IF (WriteIPVar(6)) THEN
      CgwTTAtIPVar => VariableGet( Solver % Mesh % Variables, 'CgwTT')
      IF (.NOT.ASSOCIATED(CgwTTAtIPVar)) THEN
        WRITE(Message,*) 'Variable "CgwTT" is not associated'
        CALL FATAL(SolverName,Message)
      END IF
      CgwTTAtIPPerm => CgwTTAtIPVar % Perm
      CgwTTAtIP => CgwTTAtIPVar % Values
    END IF
    ! Get stuff from SIF Material section
    Material => GetMaterial(Element)
    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = ElementID  ! each element has it's own set of parameters
    ELSE
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF

    HydroGeo = GetLogical(Material,'Hydrogeological Model',Found)

    ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
    IF (ConstVal) &
        CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
    NoSalinity = GetLogical(Material,'No Salinity',Found)
    meanfactor = GetConstReal(Material,"Conductivity Arithmetic Mean Weight",Found)
    IF (.NOT.Found) THEN
      CALL INFO(FunctionName,'"Conductivity Arithmetic Mean Weight" not found. Using default unity value.',Level=9)
      meanfactor = 1.0_dp
    END IF
    MinKgw = GetConstReal( Material, &
         'Hydraulic Conductivity Limit', Found)
    IF (.NOT.Found .OR. (MinKgw <= 0.0_dp))  &
         MinKgw = 1.0D-14

    deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)

    ! Numerical integration:
    !-----------------------
    IP = GaussPointsAdapt( Element )
    DO t=1,IP % n
      IPPerm = XiAtIPPerm(ElementID) + t
      
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )

      ! from coordinate system
      !Weight = IP % s(t) * DetJ

      ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
      TemperatureAtIP = ListGetElementReal( Temperature_h, Basis, Element, Found, GaussPoint=t)
      IF (.NOT.Found) CALL FATAL(SolverName,'Temperature not found')
      PorosityAtIP = ListGetElementReal( Porosity_h, Basis, Element, Found, GaussPoint=t)
      IF (.NOT.Found) CALL FATAL(SolverName,'Porosity not found')
      PressureAtIP = ListGetElementReal( Pressure_h, Basis, Element, Found, GaussPoint=t)
      IF (.NOT.Found) CALL FATAL(SolverName,'Pressure not found')
      IF (NoSalinity) THEN
        SalinityAtIP = 0.0_dp
      ELSE
        SalinityAtIP = ListGetElementReal( Salinity_h, Basis, Element, Found, GaussPoint=t)
        IF (.NOT.Found) THEN
          CALL INFO(SolverName,'Salinity not found - setting to zero',Level=7)
          NoSalinity=.TRUE.
        END IF
      END IF
      ! Materialproperties (basically densities and derivatives of it)
      ! needed at IP for Xi computation (anything ELSE thereafter)
      rhosAtIP = rhos(RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)      
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
      Xi0Tilde = GetXi0Tilde(RockMaterialID,PorosityAtIP)
      
      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('anderson') ! classic simpified Anderson model
        XiAtIP(IPPerm) = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,GlobalRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiTAtIP = &
             XiAndersonT(XiAtIP(IPPerm),0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,GlobalRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiPAtIP   = &
             XiAndersonP(XiAtIp(IPPerm),0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,GlobalRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)        
      CASE DEFAULT ! Hartikainen model
        CALL  GetXiHartikainen(RockMaterialID,&
             CurrentSoluteMaterial,CurrentSolventMaterial,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
             Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
             GasConstant,p0,T0,&
             XiAtIP(IPPerm),XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
             .FALSE.,.TRUE.,.FALSE.,.TRUE.,.FALSE.)
      END SELECT

      ! compute salt density
      IF (NoSalinity) THEN
        rhocAtIP    = 0.0_dp
      ELSE
        rhocAtIP    = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP(IPPerm),TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
      END IF

      ! compute values for output
      IF (WriteIPVar(1)) THEN
        rhowAtIP  = rhowupdate(CurrentSolventMaterial,rhowAtIP,XiAtIP(IPPerm),SalinityAtIP,ConstVal)
         
        rhogwAtIP(rhogwAtIPPerm(ElementID) +t) = rhogw(rhowAtIP,rhocAtIP,XiAtIP(IPPerm),SalinityAtIP)
      END IF
      IF (WriteIPVar(2)) THEN
        mugwAtIP(mugwATIPPerm(ElementID) +t) = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
             XiAtIP(IPPerm),T0,SalinityAtIP,TemperatureAtIP,ConstVal)
      END IF
      IF (WriteIPVar(3)) THEN
        auxtensor = &
                 GetKGpe(RockMaterialID,CurrentSolventMaterial,XiAtIp(IPPerm))
        K = 0
        DO I=1,DIM
          DO J=1,DIM
            K = K + 1
            kGpeAtIP( kGpeDOFs*( (kGpeATIPPerm(ElementID) + t) - 1) + K) = &
                 auxtensor(I,J)
          END DO
        END DO
      END IF
      IF (WriteIPVar(4)) THEN
        ksthAtIP = GetKalphath(GlobalRockMaterial % ks0th(RockMaterialID),&
             GlobalRockMaterial % bs(RockMaterialID),T0,TemperatureAtIP)
        !PRINT *,"ksthAtIP:",ksthAtIP
        kwthAtIP = GetKalphath(CurrentSolventMaterial % kw0th,CurrentSolventMaterial % bw,T0,TemperatureAtIP)
        !PRINT *,"kwthAtIP:",kwthAtIP
        kithAtIP = GetKalphath(CurrentSolventMaterial % ki0th,CurrentSolventMaterial % bi,T0,TemperatureAtIP)
        !PRINT *,"kithAtIP",kithAtIP
        kcthAtIP = GetKalphath(CurrentSoluteMaterial % kc0th,CurrentSoluteMaterial % bc,T0,TemperatureAtIP)
        !PRINT *,"kcthAtIP:",kcthAtIP
        auxtensor = &
             GetKGTT(ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,XiAtIP(IPPerm),&
             SalinityATIP,PorosityAtIP,meanfactor)
        !PRINT *, "KGTT:", auxtensor
        K = 0
        DO I=1,DIM
          DO J=1,DIM
            K = K + 1
            KGTTAtIP( KGTTDOFs*( (KGTTATIPPerm(ElementID) + t) - 1) + K) = &
                 auxtensor(I,J)
          END DO
        END DO
      END IF
      IF (WriteIPVar(5)) THEN
        auxtensor = &
             GetXikG0hy(RockMaterialID,XiAtIp(IPPerm))
        !PRINT *, "XikG0hy", auxtensor
        K = 0
        DO I=1,DIM
          DO J=1,DIM
            K = K + 1
            XikG0hyAtIP( XikG0hyDOFs*( (XikG0hyATIPPerm(ElementID) + t) - 1) + K) = &
                 auxtensor(I,J)
          END DO
        END DO
      END IF
      IF (WriteIPVar(6)) THEN
        rhowAtIP = rhowupdate(CurrentSolventMaterial,rhowAtIP,XiAtIP(IPPerm),SalinityAtIP,ConstVal)
        rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP(IPPerm),TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
        cwAtIP = cw(CurrentSolventMaterial,&
           T0,XiAtIP(IPPerm),TemperatureAtIP,SalinityAtIP,ConstVal)
        ccAtIP = cc(CurrentSoluteMaterial,&
           T0,TemperatureAtIP,SalinityAtIP,ConstVal)       
        CgwTTAtIP(CgwTTATIPPerm(ElementID) +t) = &
             GetCgwTT(rhowAtIP,rhocAtIP,cwAtIP,ccAtIP,XiAtIP(IPPerm),SalinityAtIP)
      END IF
    END DO ! end loop over IP points
  END SUBROUTINE SetIPValues
END SUBROUTINE PermafrostIPOutput
!------------------------------------------------------------------------------
!> Initialize unfrozen water content from state variables
!> (temperature, pressure, salinity and porosity)
!> \ingroup Solvers
SUBROUTINE InitiliazeXi( Model,Solver,dt,TransientSimulation )
  !------------------------------------------------------------------------------
  USE DefUtils
  USE PermaFrostMaterials

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Params, Material
  TYPE(Variable_t), POINTER :: TemperatureVar,PressureVar,PorosityVar,SalinityVar
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter, maxiter, istat,DepthDOFs
  INTEGER,PARAMETER :: io=23
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),&
       PorosityPerm(:),SalinityPerm(:),GWfluxPerm1(:),&
       TemperatureDtPerm(:), PressureDtPerm(:), SalinityDtPerm(:),&
       GWfluxPerm2(:),GWfluxPerm3(:), DepthPerm(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),&
       TemperatureDt(:), PressureDt(:), SalinityDt(:),&
       GWflux1(:),GWflux2(:),GWflux3(:), Depth(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE., FluxOutput = .FALSE.,&
       ComputeDt=.FALSE.,ElementWiseRockMaterial, DepthExists=.FALSE.,&
       InitializeSteadyState=.FALSE.,ActiveMassMatrix=.TRUE.
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='InitiliazeXi'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName, GWfluxName, PhaseChangeModel,&
       ElementRockMaterialName,VarName, DepthName, XiAtIPName
  TYPE(ValueHandle_t) :: Load_h, Temperature_h, Pressure_h, Salinity_h, Porosity_h

  SAVE DIM,FirstTime,AllocationsDone,FluxOutput,DepthName,XiAtIPName,&
       CurrentSoluteMaterial,CurrentSolventMaterial,NumberOfRockRecords,&
       ElementWiseRockMaterial,ComputeDt,DepthExists,&
       Load_h, Temperature_h, Pressure_h, Salinity_h, Porosity_h
       
  !------------------------------------------------------------------------------
  CALL INFO( SolverName, '-------------------------------------',Level=4 )
  CALL INFO( SolverName, 'Computing heat transfer              ',Level=4 )
  CALL INFO( SolverName, '-------------------------------------',Level=4 )

  IF (FirstTime) THEN
    DIM = CoordinateSystemDimension()
    ! Handles to other system variables
    CALL ListInitElementKeyword( Temperature_h, 'Material', 'Temperature Variable' )
    CALL ListInitElementKeyword( Pressure_h, 'Material', 'Pressure Variable' )
    CALL ListInitElementKeyword( Salinity_h, 'Material', 'Salinity Variable' )
    CALL ListInitElementKeyword( Porosity_h, 'Material', 'Porosity Variable' )

  END IF
  
  !CALL DefaultStart()

  VarName = Solver % Variable % Name
  Params => GetSolverParams()

  !CALL DefaultInitialize()
  Active = GetNOFActive()
  DO t=1,Active
    Element => GetActiveElement(t)
    Material => GetMaterial()


    IF (FirstTime) THEN
      ! check, whether we have globally or element-wise defined values of rock-material parameters
      ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
      IF (ElementWiseRockMaterial) THEN
        CALL INFO(SolverName,'Found "Element Rock Material File"',Level=5)
        CALL INFO(SolverName,'Using element-wise rock material definition',Level=5)
      END IF
      IF (ElementWiseRockMaterial) THEN
        ! read element-wise material parameter (GlobalRockMaterial will have one entry each element)
        NumberOfRockRecords = &
             ReadPermafrostElementRockMaterial(ElementRockMaterialName,Solver,DIM)
      ELSE
        NumberOfRockRecords =  ReadPermafrostRockMaterial( Material )
      END IF

      IF (NumberOfRockRecords < 1) THEN
        CALL FATAL(SolverName,'No Rock Material specified')
      ELSE
        CALL INFO(SolverName,'Permafrost Rock Material read',Level=6)
        FirstTime = .FALSE.
      END IF
      CALL ReadPermafrostSoluteMaterial( Material,Model % Constants,CurrentSoluteMaterial )
      CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
    END IF

    n  = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()

    PhaseChangeModel = ListGetString(Material, &
         'Permafrost Phase Change Model', Found )
    IF (Found) THEN
      WRITE (Message,'(A,A)') '"Permafrost Phase Change Model" set to ', TRIM(PhaseChangeModel)
      CALL INFO(SolverName,Message,Level=9)
    END IF

    CALL EvaluateXi(  Element, Element % ElementIndex, Active, n, nd+nb,&
         CurrentSoluteMaterial, CurrentSolventMaterial,&
         NumberOfRockRecords, PhaseChangeModel,ElementWiseRockMaterial,&
         ActiveMassMatrix,FluxOutput)
  END DO

    ! And finally, no need to solve:
    !--------------------
    !Norm = DefaultSolve()

    !IF( Solver % Variable % NonlinConverged > 0 ) EXIT


    !CALL DefaultFinish()

CONTAINS

  ! Assembly of the matrix entries arising from the bulk elements
  !------------------------------------------------------------------------------
  SUBROUTINE EvaluateXi(  Element, ElementID, NoElements, n, nd,&
       CurrentSoluteMaterial, CurrentSolventMaterial,&
       NumberOfRockRecords, PhaseChangeModel, ElementWiseRockMaterial,&
       ActiveMassMatrix,FluxOutput)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n, nd, ElementID, NoElements, NumberOfRockRecords
    TYPE(Element_t), POINTER :: Element
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    LOGICAL, INTENT(IN) :: ElementWiseRockMaterial,ActiveMassMatrix, FluxOutput
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: DepthAtIP,RefDepth,CGTTAtIP, CgwTTAtIP, CGTpAtIP, CGTycAtIP,KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,XiEtaAtIP,&
         ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,hiAtIP,hwAtIP  ! function values needed for C's and KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP, bijAtIP(2,2), bijYcAtIP(2,2),&
         gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) ::  gradTAtIP(3),gradPAtIP(3),JgwDAtIP(3),KgwAtIP(3,3),KgwpTAtIP(3,3),MinKgw,&
         KgwppAtIP(3,3),fwAtIP,mugwAtIP,DtdAtIP(3,3)!  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0, eps, Gravity(3) ! constants read only once
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,rhogwAtIP,csAtIP,cwAtIP,ciAtIP,ccAtIP ! material properties at IP
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,&
         PressureVeloAtIP,SalinityVeloAtIP,&
         StiffPQ, meanfactor, vstarAtIP(3)
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp), POINTER :: gWork(:,:)
    INTEGER :: i,t,p,q,IPPerm,DIM, RockMaterialID, FluxDOFs
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.,ConstVal=.FALSE.,&
         CryogenicSuction=.FALSE.,HydroGeo=.FALSE.,ComputeFlux=.TRUE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN) :: MaterialFileName
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='Permafrost(LocalMatrixHTEQ)'
    TYPE(Variable_t), POINTER :: XiAtIPVar, GWfluxVar1, GWfluxVar2, GWfluxVar3
    INTEGER, POINTER :: XiAtIPPerm(:),GWfluxPerm(:)
    REAL(KIND=dp), POINTER :: XiAtIP(:), FluxAtElem(:)

    !------------------------------------------------------------------------------
    SAVE Nodes, ConstantsRead, ConstVal,DIM, GasConstant, N0,DeltaT, T0, p0, eps, Gravity
    !------------------------------------------------------------------------------
    gradTAtIP = 0.0_dp
    gradPAtIP = 0.0_dp
    IF(.NOT.ConstantsRead) THEN
      ConstantsRead = &
           ReadPermafrostConstants(Model, FunctionName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)
    END IF

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    XiAtIPVar => VariableGet( Solver % Mesh % Variables, 'Xi')
    IF (.NOT.ASSOCIATED(XiAtIPVar)) THEN
      WRITE(Message,*) 'Variable Xi is not associated'
      CALL FATAL(SolverName,Message)
    END IF
    XiAtIPPerm => XiAtIPVar % Perm
    XiAtIp => XiAtIPVar % Values

    ! Get stuff from SIF Material section
    Material => GetMaterial(Element)
    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = ElementID  ! each element has it's own set of parameters
    ELSE
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF

    HydroGeo = GetLogical(Material,'Hydrogeological Model',Found)
    IF (.NOT.Found) HydroGeo = .FALSE.

    ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
    IF (ConstVal) &
        CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)

    meanfactor = GetConstReal(Material,"Conductivity Arithmetic Mean Weight",Found)
    IF (.NOT.Found) THEN
      CALL INFO(FunctionName,'"Conductivity Arithmetic Mean Weight" not found. Using default unity value.',Level=9)
      meanfactor = 1.0_dp
    END IF
    MinKgw = GetConstReal( Material, &
         'Hydraulic Conductivity Limit', Found)
    IF (.NOT.Found .OR. (MinKgw <= 0.0_dp))  &
         MinKgw = 1.0D-14

    deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)

    ! Numerical integration:
    !-----------------------
    IP = GaussPointsAdapt( Element )
    IF( Element % ElementIndex == 1 ) THEN
      CALL INFO(FunctionName,'Number of Gauss points for 1st element:'&
          //TRIM(I2S(IP % n)),Level=7)      
    END IF

    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )

      ! System variables (Temperature, Porosity, Pressure, Salinity) at IP
      PorosityAtIP = ListGetElementReal( Porosity_h, Basis, Element, Found, GaussPoint=t)
      IF (.NOT.Found) CALL FATAL(SolverName,'Porosity not found')
      PressureAtIP = ListGetElementReal( Pressure_h, Basis, Element, Found, GaussPoint=t)
      IF (.NOT.Found) CALL FATAL(SolverName,'Pressure not found')
      SalinityAtIP = ListGetElementReal( Salinity_h, Basis, Element, Found, GaussPoint=t)
      TemperatureAtIP = ListGetElementReal( Temperature_h, Basis, Element, Found, GaussPoint=t)
      IF (.NOT.Found) CALL FATAL(SolverName,'Temperature not found')
      !IF (.NOT.Found) CALL WARN(SolverName,'Salinity not found - setting to zero')

      !Materialproperties needed for computing Xi at IP
      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!      
      Xi0Tilde = GetXi0Tilde(RockMaterialID,PorosityAtIP)

      ! unfrozen pore-water content at IP
      IPPerm = XiAtIPPerm(ElementID) + t
      SELECT CASE(PhaseChangeModel)
      CASE('anderson')
        XiAtIP(IPPerm) = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,GlobalRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiTAtIP = &
             XiAndersonT(XiAtIP(IPPerm),0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,GlobalRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiPAtIP   = &
             XiAndersonP(XiAtIp(IPPerm),0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,GlobalRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)       
      CASE DEFAULT ! Hartikainen model
        CALL  GetXiHartikainen (RockMaterialID,&
             CurrentSoluteMaterial,CurrentSolventMaterial,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
             Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
             GasConstant,p0,T0,&
             XiAtIP(IPPerm),XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
             .TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.)
      END SELECT
    END DO
     !------------------------------------------------------------------------------
  END SUBROUTINE EvaluateXi
  !------------------------------------------------------------------------------
END SUBROUTINE InitiliazeXi
  











