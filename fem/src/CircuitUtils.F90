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

    ! Recorded when the circuit structures are built, compared on every later
    ! entry. See CircuitsCheckStale() for what these are good for.
    INTEGER, PRIVATE, SAVE :: BuiltNm = -1
    INTEGER, PRIVATE, SAVE :: BuiltTotN = -1
    TYPE(Mesh_t), POINTER, PRIVATE, SAVE :: BuiltMesh => NULL()

    !> Bumped by FreeCircuits(). Every routine in the package that caches
    !> something behind a "first time through" test compares its own saved copy
    !> against this instead of using a logical flag, so that one increment
    !> invalidates all of those caches at once without anything having to
    !> enumerate them. Values cached this way include dim, CSymmetry, the W
    !> potential variable and the element variable handles - all of which point
    !> at, or are derived from, structures that a rebuild replaces.
    INTEGER, SAVE :: CircuitsGeneration = 1

CONTAINS

!------------------------------------------------------------------------------
!> Release everything the circuit package allocated and invalidate its caches.
!>
!> Ownership matters here, since much of what a component points at is borrowed:
!>   owned    - Model % Circuits(:) and, per circuit, ComponentIds, Perm, names,
!>              source, A, B, Mre, Mim, CircuitVariables(:), Components(:); per
!>              circuit variable, A, B, Mre, Mim, SourceRe, SourceIm; per
!>              component, the deferred length CoilType and ComponentType; and
!>              Model % CircuitMatrix, n_Circuits and Circuit_tot_n.
!>   borrowed - Comp % BodyIds and Comp % ElBoundaries, which come straight from
!>              ListGetIntegerArray and belong to the value lists; Comp % ivar
!>              and vvar, which point into CircuitVariables; Cvar % Component,
!>              which points into Components; and every Solver pointer.
!> The borrowed ones are only nullified. Circuit % Area is in the type but is
!> never allocated, hence the guard.
!------------------------------------------------------------------------------
  SUBROUTINE FreeCircuits()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(Component_t), POINTER :: Comp
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(Solver_t), POINTER :: ASolver
    INTEGER :: p, i

    CALL Info('FreeCircuits','Releasing circuit structures',Level=7)

    ! The A solver matrix aliases the circuit matrix, so break that before the
    ! matrix goes away or it is left dangling.
    ASolver => CurrentModel % ASolver
    IF( ASSOCIATED(ASolver) ) THEN
      IF( ASSOCIATED(ASolver % Matrix) ) ASolver % Matrix % AddMatrix => NULL()
    END IF

    IF( ASSOCIATED(CurrentModel % Circuits) .AND. ASSOCIATED(CurrentModel % n_Circuits) ) THEN
      DO p=1,CurrentModel % n_Circuits
        Circuit => CurrentModel % Circuits(p)

        IF( ASSOCIATED(Circuit % CircuitVariables) ) THEN
          DO i=1,SIZE(Circuit % CircuitVariables)
            Cvar => Circuit % CircuitVariables(i)
            IF( ALLOCATED(Cvar % A) ) DEALLOCATE(Cvar % A)
            IF( ALLOCATED(Cvar % B) ) DEALLOCATE(Cvar % B)
            IF( ALLOCATED(Cvar % Mre) ) DEALLOCATE(Cvar % Mre)
            IF( ALLOCATED(Cvar % Mim) ) DEALLOCATE(Cvar % Mim)
            IF( ALLOCATED(Cvar % SourceRe) ) DEALLOCATE(Cvar % SourceRe)
            IF( ALLOCATED(Cvar % SourceIm) ) DEALLOCATE(Cvar % SourceIm)
            IF( ALLOCATED(Cvar % EqVarIds) ) DEALLOCATE(Cvar % EqVarIds)
            Cvar % Component => NULL()          ! borrowed
          END DO
        END IF

        IF( ASSOCIATED(Circuit % Components) ) THEN
          DO i=1,SIZE(Circuit % Components)
            Comp => Circuit % Components(i)
            IF( ALLOCATED(Comp % CoilType) ) DEALLOCATE(Comp % CoilType)
            IF( ALLOCATED(Comp % ComponentType) ) DEALLOCATE(Comp % ComponentType)
            Comp % BodyIds => NULL()            ! borrowed from a value list
            Comp % ElBoundaries => NULL()       ! borrowed from a value list
            Comp % ivar => NULL()               ! points into CircuitVariables
            Comp % vvar => NULL()               ! points into CircuitVariables
            IF( ALLOCATED(Comp % ElemIdx) ) DEALLOCATE(Comp % ElemIdx)
            IF( ALLOCATED(Comp % BCElemIdx) ) DEALLOCATE(Comp % BCElemIdx)
          END DO
          DEALLOCATE(Circuit % Components)
          Circuit % Components => NULL()
        END IF

        IF( ASSOCIATED(Circuit % CircuitVariables) ) THEN
          DEALLOCATE(Circuit % CircuitVariables)
          Circuit % CircuitVariables => NULL()
        END IF

        IF( ALLOCATED(Circuit % ComponentIds) ) DEALLOCATE(Circuit % ComponentIds)
        IF( ALLOCATED(Circuit % Perm) ) DEALLOCATE(Circuit % Perm)
        IF( ALLOCATED(Circuit % names) ) DEALLOCATE(Circuit % names)
        IF( ALLOCATED(Circuit % source) ) DEALLOCATE(Circuit % source)
        IF( ALLOCATED(Circuit % A) ) DEALLOCATE(Circuit % A)
        IF( ALLOCATED(Circuit % B) ) DEALLOCATE(Circuit % B)
        IF( ALLOCATED(Circuit % Mre) ) DEALLOCATE(Circuit % Mre)
        IF( ALLOCATED(Circuit % Mim) ) DEALLOCATE(Circuit % Mim)
        IF( ALLOCATED(Circuit % Area) ) DEALLOCATE(Circuit % Area)

        Circuit % ASolver => NULL()             ! borrowed
        Circuit % n = 0
        Circuit % n_comp = 0
        Circuit % UsePerm = .FALSE.
      END DO

      DEALLOCATE(CurrentModel % Circuits)
    END IF
    CurrentModel % Circuits => NULL()

    IF( ASSOCIATED(CurrentModel % CircuitMatrix) ) THEN
      CALL FreeMatrix(CurrentModel % CircuitMatrix)
      CurrentModel % CircuitMatrix => NULL()
    END IF

    IF( ASSOCIATED(CurrentModel % n_Circuits) ) DEALLOCATE(CurrentModel % n_Circuits)
    IF( ASSOCIATED(CurrentModel % Circuit_tot_n) ) DEALLOCATE(CurrentModel % Circuit_tot_n)
    CurrentModel % n_Circuits => NULL()
    CurrentModel % Circuit_tot_n => NULL()
    CurrentModel % ASolver => NULL()            ! borrowed

    ! Invalidate every cached value in the package, and forget what was built.
    CircuitsGeneration = CircuitsGeneration + 1
    BuiltNm = -1
    BuiltTotN = -1
    BuiltMesh => NULL()
!------------------------------------------------------------------------------
  END SUBROUTINE FreeCircuits
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Remember what the circuit structures were built against.
!------------------------------------------------------------------------------
  SUBROUTINE CircuitsRecordBuild()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Solver_t), POINTER :: ASolver

    ASolver => CurrentModel % ASolver
    IF(.NOT. ASSOCIATED(ASolver)) RETURN

    BuiltNm = ASolver % Matrix % NumberOfRows
    BuiltTotN = CurrentModel % Circuit_tot_n
    BuiltMesh => ASolver % Mesh

    CALL Info('CircuitsRecordBuild','Circuit structures built for '//I2S(BuiltNm)//&
        ' matrix rows and '//I2S(BuiltTotN)//' circuit dofs',Level=8)
!------------------------------------------------------------------------------
  END SUBROUTINE CircuitsRecordBuild
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Refuse to assemble into circuit structures that no longer fit the problem.
!>
!> Every circuit row is addressed as an offset from the number of rows of the A
!> solver matrix, and CircuitMatrix is sized nm + Circuit_tot_n once, in
!> Circuits_MatrixInit(). The whole package is then built exactly once, behind
!> saved flags, and never rebuilt - so if the A solver matrix or the mesh is
!> replaced underneath it (adaptivity, remeshing, a restart, a second mesh) the
!> row indices silently address the wrong rows.
!>
!> When that is detected the package is torn down through FreeCircuits(), which
!> also bumps CircuitsGeneration and so invalidates every cached value in it. The
!> caller then falls into its own build path again on this same entry.
!------------------------------------------------------------------------------
  SUBROUTINE CircuitsCheckStale()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Solver_t), POINTER :: ASolver
    INTEGER :: nm
    LOGICAL :: Stale
    CHARACTER(*), PARAMETER :: Caller = 'CircuitsCheckStale'

    IF( BuiltNm < 0 ) RETURN   ! nothing built yet

    ASolver => CurrentModel % ASolver
    IF(.NOT. ASSOCIATED(ASolver)) RETURN

    Stale = .FALSE.

    nm = ASolver % Matrix % NumberOfRows
    IF( nm /= BuiltNm ) THEN
      CALL Info(Caller,'The circuit equations were built for '//I2S(BuiltNm)//&
          ' matrix rows but the solver now has '//I2S(nm),Level=5)
      Stale = .TRUE.
    END IF

    IF( ASSOCIATED(BuiltMesh) .AND. ASSOCIATED(ASolver % Mesh) ) THEN
      IF( .NOT. ASSOCIATED(BuiltMesh, ASolver % Mesh) ) THEN
        CALL Info(Caller,'The A solver mesh is not the one the circuits were built on',Level=5)
        Stale = .TRUE.
      END IF
    END IF

    IF( .NOT. Stale ) RETURN

    CALL Info(Caller,'Circuit structures no longer fit the problem, rebuilding them',Level=4)
    CALL FreeCircuits()
!------------------------------------------------------------------------------
  END SUBROUTINE CircuitsCheckStale
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION GetCircuitModelDepth() RESULT (Depth)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    TYPE(Valuelist_t), POINTER :: simulation
    REAL(KIND=dp) :: depth
    LOGICAL :: Found, CSymmetry, Parallel
    INTEGER :: NoSlices
    
    CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )

    simulation => GetSimulation()
    depth = GetConstReal(simulation, 'Circuit Model Depth', Found)

    IF( Found ) THEN
      NoSlices = ListGetInteger(simulation,'Number of Slices',Found)
      IF(NoSlices > 1) THEN
        IF( CurrentModel % Solver % Parallel ) depth = depth / NoSlices
      END IF
    ELSE
      depth = 1._dp
      IF (CSymmetry) depth = 2._dp * pi
    END IF
        
!------------------------------------------------------------------------------
  END FUNCTION GetCircuitModelDepth
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Is the mesh really split over the partitions?
!>
!> Not the same question as Solver % Parallel, which says whether the linear
!> system is solved in parallel. With "Number of Slices" or parallel timestepping
!> every partition holds the whole mesh, and "Enforce Parallel" then makes the
!> linear system parallel while every element still exists on every partition.
!> Anything integrated over elements - an electrode area, a coil resistance, an
!> element count - is a partial sum only when this is true, and reducing it when
!> it is false multiplies the result by the number of partitions.
!------------------------------------------------------------------------------
  FUNCTION CircuitsPartitionedMesh() RESULT (Partitioned)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    LOGICAL :: Partitioned
    TYPE(Mesh_t), POINTER :: Mesh

    Partitioned = ( ParEnv % PEs > 1 )
    IF( .NOT. Partitioned ) RETURN

    Mesh => CurrentModel % Mesh
    IF( ASSOCIATED(CurrentModel % Solver) ) THEN
      IF( ASSOCIATED(CurrentModel % Solver % Mesh) ) Mesh => CurrentModel % Solver % Mesh
    END IF
    IF( ASSOCIATED(Mesh) ) Partitioned = .NOT. Mesh % SingleMesh
!------------------------------------------------------------------------------
  END FUNCTION CircuitsPartitionedMesh
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION GetComponentVoltageFactor(CompInd) RESULT (VoltageFactor)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: CompInd
    REAL(KIND=dp) :: VoltageFactor
    TYPE(Valuelist_t), POINTER :: CompParams
    LOGICAL :: Found
     
    CompParams => CurrentModel % Components(CompInd) % Values
    IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('GetComponentVoltageFactor',&
                                                        'Component parameters not found')
    VoltageFactor = GetConstReal(CompParams, 'Circuit Equation Voltage Factor', Found)
    IF (.NOT. Found) VoltageFactor = 1._dp
!------------------------------------------------------------------------------
  END FUNCTION GetComponentVoltageFactor
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION GetComponentParams(Element) RESULT (ComponentParams)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: i
    TYPE(Element_t) :: Element
    TYPE(Valuelist_t), POINTER :: ComponentParams, EntityParams
    LOGICAL :: Found
    
    EntityParams => GetBC(Element)
    IF (.NOT. ASSOCIATED(EntityParams)) THEN

      EntityParams => GetBodyParams( Element )
      IF (.NOT. ASSOCIATED(EntityParams)) CALL Fatal ('GetCompParams', 'Body Parameters not found')

    END IF
   
    i = GetInteger(EntityParams, 'Component', Found)
    
    IF (.NOT. Found) THEN
      ComponentParams => Null()
    ELSE
      ComponentParams => CurrentModel % Components(i) % Values
    END IF
   
!------------------------------------------------------------------------------
  END FUNCTION GetComponentParams
!------------------------------------------------------------------------------


! Get the current associated to the component from the solution of the
! constraint problem having lagrange multiplier vector. 
!------------------------------------------------------------------------------
  FUNCTION GetComponentCurrent(CompId,Found) RESULT ( Curr ) 
    INTEGER :: CompId
    LOGICAL :: Found
    COMPLEX(KIND=dp) :: Curr
    
    INTEGER :: i,j
    TYPE(CircuitVariable_t), POINTER :: iVar
    TYPE(Variable_t), POINTER :: LagrangeVar
    REAL(KIND=dp) :: CurrIm, CurrRe
    TYPE(Circuit_t), POINTER :: Circuit
    CHARACTER(LEN=MAX_NAME_LEN) :: str 
       
    Found = .FALSE.
    Curr = 0.0_dp

    IF(CurrentModel % n_Circuits == 0) RETURN
    
    CurrRe = 0.0_dp
    CurrIm = 0.0_dp
    
    DO i = 1, CurrentModel % n_Circuits     
      Circuit => CurrentModel % Circuits(i)
      
      str = LagrangeMultiplierName( CurrentModel % ASolver ) 
      LagrangeVar => VariableGet( CurrentModel % Mesh % Variables, str, ThisOnly = .TRUE.)
      IF(.NOT. ASSOCIATED(LagrangeVar) ) RETURN           
      
      DO j = 1, SIZE(Circuit % Components)
        ivar => Circuit % Components(j) % ivar
        IF(.NOT. ASSOCIATED(ivar)) CYCLE            
        IF(.NOT. iVar % isIvar ) CYCLE
        IF(iVar % BodyId /= CompId ) CYCLE
        IF(iVar % ValueId > 0 ) THEN
          Found = .TRUE.
          CurrRe = LagrangeVar % Values(iVar % ValueId)
        END IF
        ! Guarded on the circuit being harmonic, not just on the index looking
        ! plausible: ImValueId is only ever assigned for a harmonic circuit, so
        ! for a transient one this used to branch on uninitialized memory, and
        ! when that happened to be a valid index it read the real current as the
        ! imaginary part. A transient current is real by construction.
        IF( Circuit % Harmonic .AND. iVar % ImValueId > 0 ) THEN
          CurrIm = LagrangeVar % Values(iVar % ImValueId)
        END IF        
        IF(Found) EXIT
      END DO
      IF(Found) EXIT
    END DO

    IF(.NOT. Found) THEN
      CALL Fatal('GetComponentCurrent','Got circuits but no current for component: '//I2S(CompId))
    END IF
      
    !PRINT *,'Curr:',CompId,CurrRe,CurrIm
    Curr = CMPLX(CurrRe,CurrIm,KIND=dp)

  END FUNCTION GetComponentCurrent


  ! Get the current associated to the component from the solution of the
! constraint problem having lagrange multiplier vector. 
!------------------------------------------------------------------------------
  FUNCTION GetComponentArea(CompId,Found) RESULT ( Area ) 
    INTEGER :: CompId
    LOGICAL :: Found
    REAL(KIND=dp) :: Area
    
    INTEGER :: i,j
    TYPE(CircuitVariable_t), POINTER :: iVar
    TYPE(Variable_t), POINTER :: LagrangeVar
    REAL(KIND=dp) :: CurrIm, CurrRe
    TYPE(Circuit_t), POINTER :: Circuit
    CHARACTER(LEN=MAX_NAME_LEN) :: str 
    LOGICAL :: GotComp
       
    Found = .FALSE.
    GotComp = .FALSE.
    Area = 0.0_dp

    IF(CurrentModel % n_Circuits == 0) RETURN
    
    DO i = 1, CurrentModel % n_Circuits     
      Circuit => CurrentModel % Circuits(i)      
      DO j = 1, SIZE(Circuit % Components)
        ivar => Circuit % Components(j) % ivar
        IF(.NOT. ASSOCIATED(ivar)) CYCLE
        IF(iVar % BodyId /= CompId ) CYCLE

        GotComp = .TRUE.
        Area = Circuit % Components(j) % ElArea

        ! ElArea is only computed for the coil types that need it, so a zero here
        ! means "this component has no electrode area", not "the area is zero".
        ! Report that through Found rather than returning a number the caller
        ! cannot use - it is only ever wanted for optional postprocessing, so it
        ! is not worth aborting the run over.
        Found = ( Area > 0.0_dp )
        IF( .NOT. Found ) THEN
          CALL Warn('GetComponentArea','No electrode area for component '&
              //I2S(CompId)//' of coil type: '//TRIM(Circuit % Components(j) % CoilType))
        END IF
        EXIT
      END DO
      IF(GotComp) EXIT
    END DO

    IF(.NOT. GotComp) THEN
      CALL Fatal('GetComponentArea','Got circuits but no area for component: '//I2S(CompId))
    END IF
      
  END FUNCTION GetComponentArea

  
!------------------------------------------------------------------------------
  FUNCTION GetComponentId(Element) RESULT (ComponentId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: ComponentId
    TYPE(Element_t) :: Element
    TYPE(Valuelist_t), POINTER :: BodyParams
    LOGICAL :: Found
    
    BodyParams => GetBodyParams( Element )
    IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('GetCompParams', 'Body Parameters not found')
   
    ComponentId = GetInteger(BodyParams, 'Component', Found)
    
!------------------------------------------------------------------------------
  END FUNCTION GetComponentId
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
  SUBROUTINE GetWPotentialVar(pVar)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Variable_t), POINTER :: pVar

    pVar => VariableGet( CurrentModel % Mesh % Variables,'W Potential')
    IF(.NOT. ASSOCIATED(pVar) ) THEN
      pVar => VariableGet( CurrentModel % Mesh % Variables,'W')
    END IF
    IF(ASSOCIATED(pVar)) THEN
      CALL Info('GetWPotentialVar','Using gradient of field to define direction: '&
          //TRIM(pVar % Name),Level=7)
    ELSE
      CALL Warn('GetWPotentialVar','Could not obtain variable for potential "W"')
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GetWPotentialVar
!------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
  SUBROUTINE AddComponentsToBodyLists()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    LOGICAL :: Found
    INTEGER :: i, j, k
    INTEGER, SAVE :: MyGen = -1
    ! Components and Bodies:
    ! ----------------------  
    INTEGER :: BodyId, BoundaryId
    INTEGER, POINTER :: BodyAssociations(:) => Null()
    INTEGER, POINTER :: BCAssociations(:) => Null()
    TYPE(Valuelist_t), POINTER :: BodyParams, BCParams, ComponentParams
     
    IF (MyGen == CircuitsGeneration) RETURN

    MyGen = CircuitsGeneration
    DO i = 1, SIZE(CurrentModel % Components)
      ComponentParams => CurrentModel % Components(i) % Values

      IF( ListGetLogical( ComponentParams,'Passive Component', Found ) ) CYCLE 
      
      IF (.NOT. ASSOCIATED(ComponentParams)) CALL Fatal ('AddComponentsToBodyList', &
                                                         'Component parameters not found!')
      BodyAssociations => ListGetIntegerArray(ComponentParams, 'Body', Found)

      IF (.NOT. Found) BodyAssociations => ListGetIntegerArray(ComponentParams, 'Master Bodies', Found)

      IF (.NOT. Found) BCAssociations => ListGetIntegerArray(ComponentParams, 'Master BCs', Found)
      
      IF (.NOT. Found) CYCLE

      IF (ASSOCIATED(BodyAssociations)) THEN
        DO j = 1, SIZE(BodyAssociations)
          BodyId = BodyAssociations(j)
          BodyParams => CurrentModel % Bodies(BodyId) % Values
          IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('AddComponentsToBodyList', &
                                                        'Body parameters not found!')
          ! Idempotent on purpose. This writes into the body's value list, which
          ! outlives the circuit package, so on a rebuild the keyword is already
          ! there from the previous build. Only a body claimed by a *different*
          ! component is the error this check is for.
          k = GetInteger(BodyParams, 'Component', Found)
          IF (Found) THEN
            IF( k /= i ) CALL Fatal ('AddComponentsToBodyList', &
                'Body '//TRIM(i2s(BodyId))//' associated to two components!')
          ELSE
            CALL listAddInteger(BodyParams, 'Component', i)
          END IF
          BodyParams => Null()
        END DO
      END IF
      IF (ASSOCIATED(BCAssociations)) THEN
        DO j = 1, SIZE(BCAssociations)
          BoundaryId = BCAssociations(j)
          BCParams => CurrentModel % BCs(BoundaryId) % Values
          IF (.NOT. ASSOCIATED(BCParams)) CALL Fatal ('AddComponentsToBodyList', &
                         'Boundary Condition parameters not found!')

          ! Idempotent for the same reason as the body branch above.
          k = GetInteger(BCParams, 'Component', Found)

          IF (Found) THEN
            IF( k /= i ) CALL Fatal ('AddComponentsToBodyList', &
                'Boundary Condition '//TRIM(i2s(BoundaryId))//' associated to two components!')
          ELSE
            CALL ListAddInteger(BCParams, 'Component', i)
          END IF
          BCParams => Null()
        END DO
      END IF
      BodyAssociations => Null()
      BCAssociations => Null()
    END DO

    DO i = 1, CurrentModel % NumberOfBodies
      BodyParams => CurrentModel % Bodies(i) % Values
      IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('AddComponentsToBodyList', &
                          'Body parameters not found!')
      j = GetInteger(BodyParams, 'Component', Found)
      IF (.NOT. Found) CYCLE

      WRITE(Message,'(A)') '"Body '//TRIM(I2S(i))//'" associated to "Component '//TRIM(I2S(j))//'"' 
      CALL Info('AddComponentsToBodyList',Message,Level=5)
      BodyParams => Null()
    END DO

    DO i = 1, CurrentModel % NumberOfBCs
      BCParams => CurrentModel % BCs(i) % Values
      IF (.NOT. ASSOCIATED(BCParams)) CALL Fatal ('AddComponentsToBodyList', &
                  'Boundary Condition parameters not found!')
      j = GetInteger(BCParams, 'Component', Found)
      IF (.NOT. Found) CYCLE

      WRITE(Message,'(A)') '"Boundary Condition '//TRIM(I2S(i))// &
          '" associated to "Component '//TRIM(I2S(j))//'"' 
      CALL Info('AddComponentsToBodyList',Message,Level=5)
      BCParams => Null()
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE AddComponentsToBodyLists
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE CheckComponentVariables()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    LOGICAL :: Found
    INTEGER :: i, j, k
    INTEGER, SAVE :: MyGen = -1
    TYPE(Valuelist_t), POINTER :: ComponentParams
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType, VarName
   
    IF (MyGen == CircuitsGeneration) RETURN
    IF(CurrentModel % NumberOfComponents == 0) RETURN
    
    MyGen = CircuitsGeneration

    j = 0
    DO i = 1, CurrentModel % NumberOfComponents      
      ComponentParams => CurrentModel % Components(i) % Values                  
      IF( ListGetLogical( ComponentParams,'Passive Component', Found ) ) CYCLE 
      CoilType = GetString(ComponentParams, 'Coil Type', Found)
      IF(.NOT. Found) CYCLE

      SELECT CASE (CoilType)
      CASE ('stranded')
        VarName = 'Circuit Current Variable Id'
      CASE ('massive','foil winding')
        VarName = 'Circuit Voltage Variable Id'
      CASE DEFAULT
        CYCLE
      END SELECT

      ! The keyword is written by AddComponentValuesToLists(), which only walks
      ! the components that a circuit actually refers to. Its absence therefore
      ! means this component has a coil type but is not wired into any circuit.
      IF(.NOT. ListCheckPresent(ComponentParams, VarName) ) THEN
        CALL Warn('CheckComponentVariables','Component '//I2S(i)//' has a coil type but is '//&
            'not used by any circuit, so it has no: '//TRIM(VarName))
        j = j + 1
      END IF
    END DO

    IF(j > 0) THEN
      CALL Warn('CheckComponentVariables','Could not find variables for '//I2S(j)//' coils!')
      CALL Fatal('CheckComponentVariables','Check your circuit settings!')
    END IF

  END SUBROUTINE CheckComponentVariables
      

  
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
  FUNCTION FindSolverWithKey(key) RESULT (Solver)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    CHARACTER(*) :: key

    LOGICAL :: Found
    INTEGER :: i
    TYPE(Solver_t), POINTER :: Solver
    
    ! Look for the solver we attach the circuit equations to:
    ! -------------------------------------------------------
    Found = .FALSE.
    DO i=1, CurrentModel % NumberOfSolvers
      Solver => CurrentModel % Solvers(i)
      IF(ListCheckPresent(Solver % Values, key)) THEN 
        Found = .TRUE. 
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
  USE CircuitUtils
  IMPLICIT NONE

CONTAINS 

!------------------------------------------------------------------------------
  SUBROUTINE AllocateCircuitsList()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: slen,n_Circuits
    CHARACTER(:), ALLOCATABLE :: cmd
    CHARACTER(LEN=MAX_NAME_LEN) :: name

    ! Read Circuit definitions from MATC:
    ! ----------------------------------
    n_Circuits = NINT(GetMatcReal("Circuits"))
    CurrentModel % n_Circuits = n_Circuits

    IF( ASSOCIATED( CurrentModel % Circuits ) ) THEN
      IF( SIZE( CurrentModel % Circuits ) == n_Circuits ) THEN
        CALL Info('AllocateCircuitList','Circuit list already allocated!')
      ELSE
        ! This used to warn that it was deallocating the list and then keep it,
        ! leaving every later index out of range. There is no safe way to drop it
        ! either, since the circuits are referenced from Model % Circuits and from
        ! each Circuit % Asolver, so refuse instead.
        CALL Fatal('AllocateCircuitList','Circuit list already allocated with '//&
            I2S(SIZE(CurrentModel % Circuits))//' circuits, cannot resize to '//I2S(n_Circuits))
      END IF
    END IF

    IF(.NOT. ASSOCIATED(CurrentModel % Circuits ) ) THEN
      ALLOCATE( CurrentModel % Circuits(n_Circuits) )
    END IF
      
!------------------------------------------------------------------------------
  END SUBROUTINE AllocateCircuitsList
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION CountNofCircVarsOfType(CId, Var_type) RESULT (nofc)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nofc, char_len, slen, CId, i
    CHARACTER(LEN=*) :: Var_type
    CHARACTER(:), ALLOCATABLE :: cmd
    CHARACTER(LEN=MAX_NAME_LEN) :: name
    TYPE(Circuit_t), POINTER :: Circuit
    
    Circuit => CurrentModel % Circuits(CId)
    
    nofc = 0
    
    char_len = LEN_TRIM(Var_type)
    DO i=1,Circuit % n
      slen = Matc('C.'//i2s(CId)//'.name.'//i2s(i),name)
      IF(name(1:char_len) == Var_type(1:char_len)) nofc = nofc + 1
    END DO

!------------------------------------------------------------------------------
  END FUNCTION CountNofCircVarsOfType
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION CountNofCircComponents(CId, nofvar) RESULT (nofc)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nofc, nofvar, slen, CId, i, j, CompId, ibracket
    INTEGER :: ComponentIDs(nofvar)
    TYPE(Circuit_t), POINTER :: Circuit
    CHARACTER(:), ALLOCATABLE :: cmd
    CHARACTER(LEN=MAX_NAME_LEN) :: name

    nofc = 0
    ComponentIDs = -1
    
    Circuit => CurrentModel % Circuits(CId)
    
   
    DO i=1,Circuit % n
      slen = Matc('C.'//i2s(CId)//'.name.'//i2s(i),name)

      IF(isComponentName(name,slen)) THEN
        DO ibracket=1,slen
          IF(name(ibracket:ibracket)=='(') EXIT 
        END DO

        DO j=ibracket+1,slen
          IF(name(j:j)==')') EXIT 
        END DO

        READ(name(ibracket+1:j-1),*) CompId

        IF (.NOT. ANY(ComponentIDs == CompID)) nofc = nofc + 1
        ComponentIDs(i) = CompId
      END IF
    
    END DO

!------------------------------------------------------------------------------
  END FUNCTION CountNofCircComponents
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
FUNCTION isComponentName(name, len) RESULT(L)
!------------------------------------------------------------------------------
   CHARACTER(LEN=*) :: name
   INTEGER :: len
   LOGICAL :: L
   
   L = .FALSE.
   IF(len<12) RETURN
   IF(name(1:12)=='i_component(' .OR. &
      name(1:12)=='v_component(') THEN
     L = .TRUE.
   ELSE IF(len>=14) THEN
     ! Needs its own guard: 14 characters were read behind a check for 12.
     L = (name(1:14)=='phi_component(')
   END IF
!------------------------------------------------------------------------------
END FUNCTION isComponentName
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE ReadCircuitVariables(CId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: slen, ComponentId,i,j,CId, CompInd, nofc, ibracket
    LOGICAL :: Found
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(CircuitVariable_t), POINTER :: CVar
    ! Not initialized in the declaration: that would give it an implicit SAVE and
    ! carry a value over from the previous component, which combined with the
    ! aliasing below made the outcome depend on the order the variables are read.
    LOGICAL :: LondonEquations
    CHARACTER(:), ALLOCATABLE :: cmd
    CHARACTER(LEN=MAX_NAME_LEN) :: name

    LondonEquations = .FALSE.

    Circuit => CurrentModel % Circuits(CId)

    nofc = SIZE(Circuit % Components)
    DO i=1,nofc
      Circuit % Components(i) % ComponentId = 0
    END DO

    CompInd = 0
    DO i=1,Circuit % n
      slen = Matc('C.'//i2s(CId)//'.name.'//i2s(i),name)
      Circuit % names(i) = name(1:slen)

      CVar => Circuit % CircuitVariables(i)
      CVar % isIvar = .FALSE.
      CVar % isVvar = .FALSE.
      CVar % Component => Null()

      IF(isComponentName(name,slen)) THEN
        DO ibracket=1,slen
          IF(name(ibracket:ibracket)=='(') EXIT 
        END DO

        DO j=ibracket+1,slen
          IF(name(j:j)==')') EXIT 
        END DO
        READ(name(ibracket+1:j-1),*) ComponentId

        CVar % BodyId = ComponentId
        
        DO j=1,nofc
          Cvar % Component => Circuit % Components(j)
          IF(CVar % Component % ComponentId==ComponentId) EXIT
        END DO

        IF(CVar % Component % ComponentID /= ComponentId ) THEN
          CompInd = CompInd + 1
          CVar % Component => Circuit % Components(CompInd)
        END IF

        Cvar % Component % ComponentId = ComponentId

        SELECT CASE (name(1:ibracket))
        CASE('i_component(')
          CVar % isIvar = .TRUE.
          CVar % Component % ivar => CVar
        CASE('v_component(')
          ! Was passing LondonEquations as both the assignment target and the
          ! Found argument, i.e. aliasing an INTENT(OUT) dummy with the variable
          ! being assigned. It happened to work because the assignment lands
          ! after the call, but the value it read back was the Found flag rather
          ! than a default, so use a separate flag.
          LondonEquations = ListGetLogical(CurrentModel % Components (ComponentId) % Values, &
                                           'London Equations', Found)
          IF(.NOT. Found) LondonEquations = .FALSE.
          IF (.NOT. LondonEquations) THEN
            CVar % isVvar = .TRUE.
            CVar % Component % vvar => CVar
          ELSE
            Cvar % Component => Null()
            CVar % isIvar = .FALSE.
            CVar % isVvar = .FALSE.
            CVar % dofs = 1
            CVar % pdofs = 0
            CVar % BodyId = 0
          END IF
        CASE('phi_component(')
          ! London equations lead to driving the a-formulation 
          ! with the so called node flux. Thus we replace 'v_component'
          ! variable with phi_component:
          ! (beta a, phi') + phi_component(1) (beta grad phi_0, grad phi') = i_component(1)
          !--------------------------------------------------------------------------------
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
    TYPE(Circuit_t), POINTER :: Circuit

    Circuit => CurrentModel % Circuits(CId)
    Circuit % n = NINT(GetMatcReal('C.'//i2s(CId)//'.variables'))
    n = Circuit % n
!------------------------------------------------------------------------------
  END FUNCTION GetNofCircVariables
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AllocateCircuit(CId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: CId,n
    TYPE(Circuit_t), POINTER :: Circuit

    Circuit => CurrentModel % Circuits(CId)
    
    n = Circuit % n
    
    ALLOCATE(Circuit % ComponentIds(n))
    ALLOCATE(Circuit % CircuitVariables(n), Circuit % Perm(n))
    ALLOCATE(Circuit % names(n), Circuit % source(n))
    ALLOCATE(Circuit % A(n,n), Circuit % B(n,n), &
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
    REAL(KIND=dp), ALLOCATABLE :: BoundaryAreas(:)
    REAL(KIND=dp) :: area
    INTEGER :: Active, t0, t, i, BCid, n, nBC
    INTEGER, POINTER :: ChildBCs(:)
    LOGICAL :: Found
    LOGICAL :: Parallel 

    Mesh => CurrentModel % Mesh

    ! Boundary areas are element integrals, so this is the partitioned-mesh
    ! question and not Solver % Parallel, which it used to ask.
    Parallel = CircuitsPartitionedMesh()

    ! Not all boundary elements are associated to a BC.
    ! Some may also be created by extrusion. 
    t0 = Mesh % NumberOfBulkElements
    nBC = 0
    DO t=1,Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(t0+t)
      IF(ASSOCIATED(Element % BoundaryInfo) ) THEN
        nBC = MAX(nBC, Element % BoundaryInfo % Constraint )
      END IF
    END DO

    nBC = MAX(nBC,CurrentModel % NumberOfBCs)
    nBC = ParallelReduction(nBC,2)
    
    ALLOCATE( BoundaryAreas(nBC) )
    BoundaryAreas = 0.0_dp
        
    DO i=1, CurrentModel % NumberOfBcs
       BC => CurrentModel % BCs(i) % Values
       IF (.NOT. ASSOCIATED(BC) ) CALL Fatal('SetBoundaryAreasToValueLists', 'Boundary not found!')
       CALL ListAddInteger(BC, 'Boundary Id', i)
    END DO
    
    Active = GetNOFBoundaryElements()
    DO t=1,Active
       Element => GetBoundaryElement(t)
       
       BC=>GetBC()
       IF (ASSOCIATED(BC) ) THEN
         BCid = GetInteger(BC, 'Boundary Id', Found)
       ELSE
         BCid = Element % BoundaryInfo % Constraint
       END IF

       IF( BCid > 0 ) THEN
         n = GetElementNOFNodes() 
         BoundaryAreas(BCid) = BoundaryAreas(BCid) + ElementAreaNoAxisTreatment(Mesh, Element, n)
       END IF
     END DO
     
     IF( Parallel ) THEN
       DO i=1, nBC
         BoundaryAreas(i) = ParallelReduction(BoundaryAreas(i))
       END DO
     END IF

     IF( InfoActive(25) ) THEN
       DO i=1,nBC
         PRINT *,'A(i)',i,i<=CurrentModel % NumberOfBCs,BoundaryAreas(i)
       END DO
     END IF
     
     DO i=1, CurrentModel % NumberOfBcs
       BC => CurrentModel % BCs(i) % Values
       IF (.NOT. ASSOCIATED(BC) ) CALL Fatal('ComputeCoilBoundaryAreas', 'Boundary not found!')
       BCid = GetInteger(BC, 'Boundary Id', Found)
       CALL ListAddConstReal(BC, 'Area', BoundaryAreas(BCid))
     END DO
     
     DO i=1, CurrentModel % NumberOfBodies
       BC => CurrentModel % Bodies(i) % Values
       ChildBCs => ListGetIntegerArray( BC,'Extruded Child BCs',Found ) 
       IF(Found) THEN
         !PRINT *,'Child BCs area:',i,BoundaryAreas(ChildBCs)
         area = SUM(BoundaryAreas(ChildBCs)) / SIZE( ChildBCs )         
         CALL ListAddConstReal(BC,'Extruded Child Area', area )         
       END IF
     END DO
    
!-------------------------------------------------------------------
 END SUBROUTINE SetBoundaryAreasToValueLists
!-------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE ReadComponents(CId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: CId, CompInd
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(Component_t), POINTER :: Comp
    TYPE(Valuelist_t), POINTER :: CompParams
    LOGICAL :: Found
    INTEGER :: ExtMaster

    CALL Info('ReadComponents','Reading component: '//I2S(Cid),Level=20)

    Circuit => CurrentModel % Circuits(CId)
    
    Circuit % CvarDofs = 0
    DO CompInd=1,Circuit % n_comp
      Comp => Circuit % Components(CompInd)
      Comp % nofcnts = 0
!        Comp % ComponentId = Circuits(p) % body(CompInd)
      Comp % BodyIds => GetComponentBodyIds(Comp % ComponentId)

      IF (.NOT. ASSOCIATED(Comp % ivar) ) THEN
        CALL FATAL('Circuits_Init', 'Current Circuit Variable is not found for Component '//i2s(Comp % ComponentId))
      ELSE IF (.NOT. ASSOCIATED(Comp % vvar) ) THEN
        CALL FATAL('Circuits_Init', 'Voltage Circuit Variable is not found for Component '//i2s(Comp % ComponentId))
      END IF

      CompParams => CurrentModel % Components (Comp % ComponentId) % Values
      IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('Circuits_Init', 'Component parameters not found!')
      
      Comp % CoilType = GetString(CompParams, 'Coil Type', Found)
      IF (.NOT. Found) THEN
        CALL Info('Circuits_Init', 'Component '//i2s(Comp % ComponentId)//&
            ' is not a coil. Checking if it has a component type.', Level=7)
        Comp % ComponentType = GetString(CompParams, 'Component Type', Found)
        IF (.NOT. Found) CALL Fatal ('Circuits_Init', 'Component Type not found!')
      ELSE
        Comp % ComponentType = 'coil'
      END IF
      
      Comp % i_multiplier_re = GetConstReal(CompParams, 'Current Multiplier re', Found)
      IF (.NOT. Found) Comp % i_multiplier_re = 0._dp
      Comp % i_multiplier_im = GetConstReal(CompParams, 'Current Multiplier im', Found)
      IF (.NOT. Found) Comp % i_multiplier_im = 0._dp

      Comp % VoltageFactor = GetConstReal(CompParams, 'Circuit Equation Voltage Factor', Found)
      IF (.NOT. Found) Comp % VoltageFactor = 1._dp

      Comp % ElBoundaries => ListGetIntegerArray(CompParams, 'Electrode Boundaries', Found)
      
      ! This is a feature intended to make it easier to extruded meshes internally with
      ! ElmerSolver. The idea is that the code knows which are the BCs that were created
      ! from extruding this 2D body. 
      ExtMaster = 0
      IF(.NOT. Found ) THEN
        IF( ListGetLogical( CurrentModel % Solver % Values,'Extruded Child BC Electrode', Found ) ) THEN          

          IF( ListGetLogical( CurrentModel % Simulation,"Extruded BCs Collect",Found ) ) THEN
            CALL Fatal('Circuits_init',&
                'Conflicting keywords: "Extruded Child BC Electrode" vs. "Extruded BCs Collect"')
          END IF

          CALL Info('Circuits_init','Setting "Extruded Child BCs"',Level=10)
          BLOCK
            INTEGER :: body_id
            INTEGER, POINTER :: pIntArray(:) => NULL()
            pIntArray => ListGetIntegerArray(CompParams, 'Body', Found )
            IF(.NOT. Found) pIntArray => ListGetIntegerArray(CompParams, 'Master Bodies', Found )
            IF( Found ) THEN
              IF(SIZE(pIntArray)==1) THEN
                body_id = pIntArray(1)
                NULLIFY(pIntArray)
                pIntArray => ListGetIntegerArray(CurrentModel % Bodies(body_id) % Values,&
                    'Extruded Child BCs',Found )
                IF(Found) THEN
                  CALL Info('Circuits_init','Setting Component '//I2S(CompInd)//' "Electrode Boundaries" to '&
                      //I2S(pIntArray(1))//' '//I2S(pIntArray(2)),Level=10)
                  Comp % ElBoundaries => pIntArray
                  ExtMaster = body_id
                END IF
              END IF
            END IF
          END BLOCK
        END IF
      END IF

      IF (Comp % ComponentType == 'resistor') THEN
       Comp % ivar % dofs = 1
        Comp % vvar % dofs = 1
        Comp % ivar % pdofs = 0
        Comp % vvar % pdofs = 0
      ELSE
        SELECT CASE (Comp % CoilType) 
        CASE ('stranded')
          
          Comp % nofturns = GetConstReal(CompParams, 'Number of Turns', Found)
          IF (.NOT. Found) CALL Fatal('Circuits_Init','Number of Turns not found!')
          
          Comp % ElArea = GetConstReal(CompParams, 'Electrode Area', Found)
          IF (.NOT. Found) THEN
            CALL ComputeElectrodeArea(Comp, CompParams, ExtMaster )
            WRITE(Message,'(A,ES12.5)') 'Component '//I2S(CompInd)//' "Electrode Area" is ',Comp % ElArea
            CALL Info('Circuits_Init',Message,Level=10)
          END IF
            
          Comp % CoilThickness = GetConstReal(CompParams, 'Coil Thickness', Found)
          IF (.NOT. Found) Comp % CoilThickness = 1._dp

          Comp % SymmetryCoeff = GetConstReal(CompParams, 'Symmetry Coefficient', Found)
          IF (.NOT. Found) Comp % SymmetryCoeff = 1.0_dp
          
          Comp % N_j = Comp % CoilThickness * Comp % nofturns / Comp % ElArea
          
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
          IF (.NOT. Found) THEN
            CALL ComputeElectrodeArea(Comp, CompParams )
            WRITE(Message,'(A,ES12.5)') 'Component '//I2S(CompInd)//' "Electrode Area" is ',Comp % ElArea
            CALL Info('Circuits_Init',Message,Level=10)
          END IF
          
          Comp % N_j = Comp % nofturns / Comp % ElArea
        END SELECT
      END IF

      CALL AddVariableToCircuit(Circuit, Comp % ivar, CId)
      CALL AddVariableToCircuit(Circuit, Comp % vvar, CId)

    END DO
      
!------------------------------------------------------------------------------
  END SUBROUTINE ReadComponents
!------------------------------------------------------------------------------

!-------------------------------------------------------------------
 SUBROUTINE ComputeElectrodeArea(Comp, CompParams, ExtMaster )
!-------------------------------------------------------------------
  USE ElementUtils
  IMPLICIT NONE
  TYPE(Component_t), POINTER :: Comp
  TYPE(ValueList_t), POINTER :: CompParams
  INTEGER, OPTIONAL :: ExtMaster 

  TYPE(ValueList_t), POINTER :: BC
  TYPE(Element_t), POINTER :: Element
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER :: t, n, BCid, NoSlices
  LOGICAL :: Found
  LOGICAL :: Parallel 
  
  Mesh => CurrentModel % Mesh
  Comp % ElArea = 0._dp

  ! The area is summed over this partition's elements, so it only needs reducing
  ! when the mesh really is split. See CircuitsPartitionedMesh().
  Parallel = CircuitsPartitionedMesh()
    
  IF (CoordinateSystemDimension() == 2) THEN
    DO t=1,GetNOFActive()
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes() 
      IF (ElAssocToComp(Element, Comp)) THEN
        Comp % ElArea = Comp % ElArea + ElementAreaNoAxisTreatment(Mesh, Element, n) 
      END IF
    END DO
    
    IF( Parallel ) THEN
      Comp % ElArea = ParallelReduction(Comp % ElArea)
    END IF

    ! Add this to list since no need to compute this twice
    CALL ListAddConstReal(CompParams,'Electrode Area',Comp % ElArea )        
  ELSE
    BCid = 0

    ! This is a special case for extruded meshes where the area is computed and stored
    ! in the master body that the extruded boundary comes from. This is treated only
    ! when the extruded master is given as a parameter.
    IF( PRESENT( ExtMaster ) ) BCid = ExtMaster
    IF( BCid > 0 ) THEN
      BC => CurrentModel % Bodies(BCid) % Values
      IF (.NOT. ASSOCIATED(BC) ) CALL Fatal('ComputeElectrodeArea', 'Master body not found!')            
      Comp % ElArea = GetConstReal(BC, 'Extruded Child Area', Found ) 
      IF (.NOT. Found) CALL Fatal('ComputeElectrodeArea', '"Extruded Child Area" not found!')
    ELSE
      IF (.NOT. ASSOCIATED(Comp % ElBoundaries)) &
          CALL Fatal('ComputeElectrodeArea','Electrode Boundaries not found')      
      BCid = Comp % ElBoundaries(1)
      IF( BCid < 1 .OR. BCid > CurrentModel % NumberOfBCs ) &     
          CALL Fatal('ComputeElectrodeArea', 'BCid is beyond range: '//I2S(BCid))

      BC => CurrentModel % BCs(BCid) % Values
      IF (.NOT. ASSOCIATED(BC) ) CALL Fatal('ComputeElectrodeArea', 'Boundary not found!')
      Comp % ElArea = GetConstReal(BC, 'Area', Found)
      IF (.NOT. Found) CALL Fatal('ComputeElectrodeArea', '"Area" not found!')
    END IF
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
     TYPE(Mesh_t), TARGET :: Mesh
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
    INTEGER :: k                        ! index of the circuit, i.e. CId
    INTEGER, POINTER :: circuit_tot_n => Null()

    CALL Info('AddVariableToCircuit','Adding variables count to circuit!',Level=20)

    Circuit_tot_n => CurrentModel % Circuit_tot_n

    ! Pick the partition owning this circuit variable. The circuit dofs are global
    ! in the sense that every partition may assemble into them; the owner only
    ! tells which partition holds the row of the constraint matrix. It ends up in
    ! CM % RowOwner and as the first entry of the row's neighbour list, both set
    ! in SetCircuitsParallelInfo().
    !
    ! The choice is static on purpose: all variables of the first circuit go to
    ! partition 0 and those of any further circuit to the last partition. There
    ! used to be a round-robin scheme here that stepped a saved owner counter
    ! downwards from ParEnv % PEs/2 (circuit 1) or ParEnv % PEs (other circuits)
    ! over successive calls, but its result was immediately overwritten by the
    ! constant assignment that followed it, so it never had any effect. It is
    ! removed rather than revived: a distribution worth having should follow the
    ! partitions that actually carry the coil elements of the component - the
    ! element counts already gathered in SetCircuitsParallelInfo() - and not a
    ! counter over the order in which the variables happen to be read.
    ! --------------------------------------------------------------------------
    IF(k==1) THEN
      Variable % Owner = 0
    ELSE
      Variable % Owner = ParEnv % PEs-1
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
    
    Circuit => CurrentModel % Circuits(CId)

    CALL Info('AddComponentValuesToLists','Adding "Circuit Voltage Variable *" keywords for '&
        //I2S(Circuit % n_comp)//' components in Circuit '//I2S(CId),Level=20)
    
    DO CompInd=1,Circuit % n_comp
 
      Comp => Circuit % Components(CompInd)   

      CompParams => CurrentModel % Components (Comp % ComponentId) % Values
      IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('Circuits_Init', 'Component Parameters not found!')

      CALL listAddInteger(CompParams, 'Circuit Voltage Variable Id', Comp % vvar % valueId)
      CALL listAddInteger(CompParams, 'Circuit Voltage Variable dofs', Comp % vvar % dofs)
      CALL listAddInteger(CompParams, 'Circuit Current Variable Id', Comp % ivar % valueId)
      CALL listAddInteger(CompParams, 'Circuit Current Variable dofs', Comp % ivar % dofs)

      ! N_j is only computed for the coil types that have a current density
      ! shape function, i.e. 'stranded' and 'foil winding'. It used to be
      ! published for every component, which handed MagnetoDynamics2D and
      ! CalcFields - both of which read it and abort if it is missing - an
      ! uninitialized value for 'massive' coils and for resistors.
      SELECT CASE( Comp % CoilType )
      CASE( 'stranded', 'foil winding' )
        CALL listAddConstReal(CompParams, 'Stranded Coil N_j', Comp % N_j)
      END SELECT

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
    
    Circuit => CurrentModel % Circuits(CId)
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

    Circuit => CurrentModel % Circuits(CId)
    n = Circuit % n

    ! Read in the coefficient matrices for the circuit equations:
    ! Ax' + Bx = source:
    ! ------------------------------------------------------------

    CALL matc_get_array('C.'//i2s(CId)//'.A'//CHAR(0),Circuit % A,n,n)
    CALL matc_get_array('C.'//i2s(CId)//'.B'//CHAR(0),Circuit % B,n,n)

    IF (Circuit % Harmonic) THEN
      ! Complex multiplier matrix is used for:
      ! B = times(M,B), where B times is the element-wise product
      ! ---------------------------------------------------------
      CALL matc_get_array('C.'//i2s(CId)//'.Mre'//CHAR(0),Circuit % Mre,n,n)
      CALL matc_get_array('C.'//i2s(CId)//'.Mim'//CHAR(0),Circuit % Mim,n,n)
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE ReadCoefficientMatrices
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReadPermutationVector(CId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: CId,n,slen,i
    TYPE(Circuit_t), POINTER :: Circuit

    Circuit => CurrentModel % Circuits(CId)
    n = Circuit % n

    DO i=1,n
      Circuit % Perm(i) = NINT(GetMatcReal('C.'//i2s(CId)//'.perm('//i2s(i-1)//')'))
    END DO
    IF(ANY(Circuit % Perm /= 0)) THEN 
      Circuit % UsePerm = .TRUE.
      CALL Info( 'ReadPermutationVector','Found Permutation vector for circuit '//i2s(CId), Level=4 )
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE ReadPermutationVector
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReadCircuitSources(CId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: CId,n,slen,i
    TYPE(Circuit_t), POINTER :: Circuit
    CHARACTER(:), ALLOCATABLE :: cmd
    CHARACTER(LEN=MAX_NAME_LEN) :: name

    Circuit => CurrentModel % Circuits(CId)
    n = Circuit % n
    DO i=1,n
      ! Names of the source functions, these functions should be found
      ! in the "Body Force 1" block of the .sif file.
      ! (nc: is for 'no check' e.g. don't abort if the MATC variable is not found!)
      ! ---------------------------------------------------------------------------
      slen = Matc('nc:C.'//i2s(CId)//'.source.'//i2s(i),name)
      Circuit % Source(i) = name(1:slen)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ReadCircuitSources
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!> Give every circuit variable its own row of the circuit coefficient matrices.
!>
!> Note that this is duplicated state, deliberately left in place. Cvar % A, B,
!> Mre and Mim are copies of row i of Circuit % A, B, Mre and Mim, so the same
!> numbers live in two places - 4n^2 doubles of them - and anything that changes
!> one has to change the other. The assembly (AddBasicCircuitEquations and the
!> Count/CreateBasicCircuitEquations pair) could index Circuit % A(i,j) directly
!> and these per-variable vectors could go away entirely; that was measured
!> against the risk and judged to be churn, since even a 200 variable circuit is
!> only a couple of megabytes and the copies are written exactly once, here.
!>
!> SourceRe and SourceIm are looser still: they are allocated with n entries per
!> variable, but only entry i of variable i is ever read or written - see the
!> source handling in AddBasicCircuitEquations - so that is n^2 storage holding n
!> values. Worth collapsing to a scalar per variable if this is ever revisited.
!------------------------------------------------------------------------------
  SUBROUTINE WriteCoeffVectorsForCircVariables(CId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: CId,n,i
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(CircuitVariable_t), POINTER :: Cvar
  
    Circuit => CurrentModel % Circuits(CId)
    n = Circuit % n

    DO i=1,n
      Cvar => Circuit % CircuitVariables(i)

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

      ! Plain copies. The tests for a nonzero source value that used to guard
      ! these were redundant, the targets having been zeroed just above.
      Cvar % A(1:n) = Circuit % A(i,1:n)
      Cvar % B(1:n) = Circuit % B(i,1:n)
      Cvar % Mre(1:n) = Circuit % Mre(i,1:n)
      Cvar % Mim(1:n) = Circuit % Mim(i,1:n)
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
     INTEGER :: k
     LOGICAL :: T, Found

     k = GetInteger(GetBC(Element), 'Component', Found)
     
     IF (Found) THEN
       T = (k .eq. Component % ComponentId)
     ELSE IF (ASSOCIATED(Component % BodyIds)) THEN
       T = IdInList(Element % BodyId, Component % BodyIds)
     ELSE
       T = .False.
     END IF

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
    LOGICAL :: support, Gate=.FALSE.
    INTEGER :: MyGen = -1
    REAL(KIND=dp) :: wBase(nn)
    TYPE(Variable_t), POINTER, SAVE :: Wpot => NULL()
    SAVE dim, MyGen, Gate

    IF (MyGen /= CircuitsGeneration) THEN
      MyGen = CircuitsGeneration
      dim = CoordinateSystemDimension()

      ! Resolve the wire direction potential through GetWPotentialVar(), which is
      ! what the assembly routines use: "W Potential" first and "W" only as a
      ! fallback. This used to ask for "W" and nothing else, so a case naming its
      ! potential "W Potential" found no support in any element and had all of its
      ! massive and foil winding elements dropped - from the matrix structure as
      ! well as from the assembly, hence without so much as a warning.
      CALL GetWPotentialVar(Wpot)

      ! With no such field at all there is nothing to test against, so do not
      ! gate on it. The component may be driven by its own current density
      ! instead ("Coil Use W Vector", "Foil Winding Use J Vector"), and answering
      ! .FALSE. here would silently discard every one of its elements.
      Gate = ( dim == 3 .AND. ASSOCIATED(Wpot) )
      IF( dim == 3 .AND. .NOT. ASSOCIATED(Wpot) ) THEN
        CALL Info('HasSupport','No potential field to test element support with, '//&
            'accepting all elements of coil components',Level=7)
      END IF
    END IF

    support = .TRUE.
    IF (Gate) THEN
      Wbase = 0.0_dp
      CALL GetLocalSolution(Wbase,UElement=Element,UVariable=Wpot)
      support = ANY( Wbase(1:nn) /= 0.0_dp )
    END IF
!------------------------------------------------------------------------------
   END FUNCTION HasSupport
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Create a standard variable associated to the mesh that may be use for dependencies.
!------------------------------------------------------------------------------
  SUBROUTINE Circuits_ToMeshVariable(Solver,crt)
    
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: Crt(:)

    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(CircuitVariable_t), POINTER :: CVar
    TYPE(Variable_t), POINTER :: Var, VarIm
    INTEGER :: p,i,n,nv,ni,m,iv,nsize
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: Found 
    CHARACTER(:), ALLOCATABLE :: CrtName,VarName,VarnameIm
    
    IF( .NOT. ListGetLogical( Solver % Values,'Export Circuit Variables',Found ) ) RETURN

    CALL Info('Circuit_ToMeshVariable','Adding circuit variables to be mesh variables')
    
    Mesh => Solver % Mesh
            
    DO p=1,CurrentModel % n_Circuits
      CALL Info('Circuit_ToMeshVariable','Adding circuit: '//I2S(p),Level=12)

      Circuit => CurrentModel % Circuits(p)

      n = Circuit % n
      
      IF( CurrentModel % n_Circuits == 1) THEN
        crtName = 'crt'
      ELSE
        crtName = 'crt '//I2S(p)
      END IF
    
      ! Count the v and i variables of the circuit.
      nv = 0; ni = 0
      DO i=1,n
        Cvar => Circuit % CircuitVariables(i)
        IF(Cvar % isIvar) ni = ni + 1
        IF(Cvar % isVvar) nv = nv + 1
      END DO
      
      IF( nv + ni == 0 ) THEN
        CALL Warn('Circuits_ToMeshVariable','No voltage or current variables exists!')
        CYCLE
      END IF
      

      ! Go first through currents and then through voltages    
      DO iv=1,2      
        IF( Circuit % Harmonic ) THEN
          IF(iv==1) THEN
            varname =  crtname//' i re'
            varnameim =  crtname//' i im'
          ELSE
            varname = crtname//' v re'
            varnameim = crtname//' v im'
          END IF
        ELSE
          IF(iv==1) THEN
            varname =  crtname//' i'
          ELSE
            varname = crtname//' v'
          END IF
        END IF
        
        IF(iv==1) THEN
          nsize = ni
        ELSE
          nsize = nv
        END IF
                
        ! Get variable, if variable does not exist then we create here on-the-fly
        Var => VariableGet( Mesh % Variables,varname)
        IF(.NOT. ASSOCIATED( Var ) ) THEN
          CALL Info('Circuits_ToMeshVariable','Creating variable: '//TRIM(varname))
          CALL VariableAddVector( Mesh % Variables, Mesh, Solver,&
              varname,dofs=nsize,global=.TRUE.)
          Var => VariableGet( Mesh % Variables,varname)
        END IF
        CALL Info('Circuts_toMeshVariable','Filling variable: '//TRIM(VarName),Level=20)

        IF( Circuit % Harmonic ) THEN
          VarIm => VariableGet( Mesh % Variables,varnameim)
          IF(.NOT. ASSOCIATED( VarIm ) ) THEN
            CALL Info('Circuits_ToMeshVariable','Creating variable: '//TRIM(varnameim))
            CALL VariableAddVector( Mesh % Variables, Mesh, Solver,&
                varnameim,dofs=nsize,global=.TRUE.)
            VarIm => VariableGet( Mesh % Variables,varnameim)
          END IF
          CALL Info('Circuts_toMeshVariable','Filling variable: '//TRIM(VarNameim),Level=20)
       END IF
          
        
        ! Fill the currents or voltages
        m = 0
        DO i=1,n
          Cvar => Circuit % CircuitVariables(i)
          
          IF(iv==1 .AND. .NOT. CVar % isIvar ) CYCLE          
          IF(iv==2 .AND. .NOT. Cvar % isVvar) CYCLE
          
          CALL Info('Circuts_toMeshVariable','Inserting variable '//I2S(CVar % ValueId)//': '&
              //TRIM(Circuit % names(i)),Level=20)
                    
          m = m + 1
          Var % Values(m) = crt(Cvar % ValueId)          
          IF(Circuit % Harmonic) THEN
            VarIm % Values(m) = crt(Cvar % ImValueId)
          END IF
        END DO
      END DO
    END DO
          
  END SUBROUTINE Circuits_ToMeshVariable


!------------------------------------------------------------------------------
!> Work out once which elements belong to which component.
!>
!> The assembly and both matrix structure passes used to be shaped as a loop over
!> components containing a loop over every element, testing ElAssocToComp() on
!> each - so the mesh was walked once per component, with a GetBC() and a keyword
!> lookup at every step. That is O(n_components * n_elements) per assembly, per
!> timestep, per nonlinear iteration.
!>
!> ElAssocToComp() is still the only thing that decides membership here, so the
!> semantics are reproduced rather than restated: this simply asks it once. The
!> lists come out in ascending element order, and the users walk them backwards to
!> keep the accumulation order - and hence the floating point result - exactly what
!> it was.
!------------------------------------------------------------------------------
  SUBROUTINE BuildComponentElementLists()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(Component_t), POINTER :: Comp
    TYPE(Element_t), POINTER :: Element
    INTEGER :: p, CompInd, q, nA, nB, nActive, nBound
    INTEGER :: TotA, TotB
    CHARACTER(*), PARAMETER :: Caller = 'BuildComponentElementLists'

    nActive = GetNOFActive()
    nBound = GetNOFBoundaryElements()
    TotA = 0; TotB = 0

    DO p=1,CurrentModel % n_Circuits
      Circuit => CurrentModel % Circuits(p)

      DO CompInd=1,Circuit % n_comp
        Comp => Circuit % Components(CompInd)

        IF( ALLOCATED(Comp % ElemIdx) ) DEALLOCATE(Comp % ElemIdx)
        IF( ALLOCATED(Comp % BCElemIdx) ) DEALLOCATE(Comp % BCElemIdx)

        ! Count, then fill. Two passes rather than a growable list, since the
        ! whole point is to pay for this once.
        nA = 0
        DO q=1,nActive
          Element => GetActiveElement(q)
          IF( ElAssocToComp(Element, Comp) ) nA = nA + 1
        END DO
        nB = 0
        DO q=1,nBound
          Element => GetBoundaryElement(q)
          IF( ElAssocToComp(Element, Comp) ) nB = nB + 1
        END DO

        ALLOCATE( Comp % ElemIdx(nA), Comp % BCElemIdx(nB) )

        nA = 0
        DO q=1,nActive
          Element => GetActiveElement(q)
          IF( ElAssocToComp(Element, Comp) ) THEN
            nA = nA + 1
            Comp % ElemIdx(nA) = q
          END IF
        END DO
        nB = 0
        DO q=1,nBound
          Element => GetBoundaryElement(q)
          IF( ElAssocToComp(Element, Comp) ) THEN
            nB = nB + 1
            Comp % BCElemIdx(nB) = q
          END IF
        END DO

        TotA = TotA + nA
        TotB = TotB + nB
      END DO
    END DO

    CALL Info(Caller,'Element lists built for the components: '//I2S(TotA)//&
        ' bulk and '//I2S(TotB)//' boundary elements out of '//I2S(nActive)//&
        ' and '//I2S(nBound),Level=8)
!------------------------------------------------------------------------------
  END SUBROUTINE BuildComponentElementLists
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Check once that every named circuit source can actually be found.
!>
!> The assembly reads the sources with GetCReal() every time it runs and falls
!> back to zero without complaint when the name resolves nowhere. Since the names
!> come from MATC definitions ("C.1.source.1 = ...") and the values from a
!> different section of the sif, a typo in either is easy to make and produces a
!> circuit that is simply not excited - a plausible looking, entirely wrong,
!> solution. Look the names up once here instead and say so.
!------------------------------------------------------------------------------
  SUBROUTINE CheckCircuitSources()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(ValueList_t), POINTER :: Params, BF
    INTEGER :: p, i, nMissing
    LOGICAL :: FoundIt
    CHARACTER(LEN=MAX_NAME_LEN) :: sname
    CHARACTER(*), PARAMETER :: Caller = 'CheckCircuitSources'

    ! The same two lists the assembly consults, in the same order. Note that the
    ! body force is hardwired to the first one there, so a source value given in
    ! "Body Force 2" is not seen - which this check therefore also reports.
    Params => CurrentModel % Solver % Values
    BF => NULL()
    IF( CurrentModel % NumberOfBodyForces > 0 ) BF => CurrentModel % BodyForces(1) % Values

    nMissing = 0

    DO p=1,CurrentModel % n_Circuits
      Circuit => CurrentModel % Circuits(p)

      DO i=1,Circuit % n
        sname = Circuit % Source(i)
        IF( LEN_TRIM(sname) == 0 ) CYCLE

        IF( Circuit % Harmonic ) THEN
          ! Either part is enough: giving only the real part is normal.
          FoundIt = SourcePresent(TRIM(sname)//' re') .OR. SourcePresent(TRIM(sname)//' im')
        ELSE
          FoundIt = SourcePresent(TRIM(sname))
        END IF

        IF( FoundIt ) THEN
          CALL Info(Caller,'Circuit '//I2S(p)//' source "'//TRIM(sname)//'" found',Level=8)
        ELSE
          nMissing = nMissing + 1
          CALL Warn(Caller,'Circuit '//I2S(p)//' names a source "'//TRIM(sname)//&
              '" that is in neither the circuit solver section nor "Body Force 1"!')
        END IF
      END DO
    END DO

    IF( nMissing > 0 ) THEN
      CALL Warn(Caller,I2S(nMissing)//' circuit source(s) will be taken as zero.')
    END IF

  CONTAINS

    FUNCTION SourcePresent(key) RESULT(yes)
      CHARACTER(LEN=*) :: key
      LOGICAL :: yes
      yes = ListCheckPresent(Params,key)
      IF( .NOT. yes .AND. ASSOCIATED(BF) ) yes = ListCheckPresent(BF,key)
    END FUNCTION SourcePresent

!------------------------------------------------------------------------------
  END SUBROUTINE CheckCircuitSources
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Reject the component definitions that the transient assembly cannot honour.
!>
!> A component whose elements come from "Master BCs" is assembled by the harmonic
!> path only: its AddComponentEquationsAndCouplings() sweeps the boundary
!> elements as well as the bulk ones, and Add_massive() carries the impedance
!> boundary terms ("Layer Electric Conductivity") that go with them. The
!> transient path sweeps bulk elements only and has no counterpart for those
!> terms, so such a component would contribute nothing at all while the matrix
!> structure pass - which does sweep the boundary elements - still reserves rows
!> for it. The result is a silently uncoupled component, so say no instead.
!------------------------------------------------------------------------------
  SUBROUTINE CheckTransientComponents()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(Component_t), POINTER :: Comp
    TYPE(ValueList_t), POINTER :: CompParams
    INTEGER :: p, CompInd
    CHARACTER(*), PARAMETER :: Caller = 'CheckTransientComponents'

    DO p=1,CurrentModel % n_Circuits
      Circuit => CurrentModel % Circuits(p)
      DO CompInd=1,Circuit % n_comp
        Comp => Circuit % Components(CompInd)
        CompParams => CurrentModel % Components(Comp % ComponentId) % Values
        IF(.NOT. ASSOCIATED(CompParams)) CYCLE

        ! Only when the boundaries really are the element association. If bodies
        ! are given too then AddComponentsToBodyLists() uses those and ignores
        ! the boundaries, which is a different (and merely confusing) matter.
        IF( ListCheckPresent(CompParams,'Body') ) CYCLE
        IF( ListCheckPresent(CompParams,'Master Bodies') ) CYCLE
        IF( .NOT. ListCheckPresent(CompParams,'Master BCs') ) CYCLE

        CALL Warn(Caller,'Component '//I2S(Comp % ComponentId)//&
            ' is defined through "Master BCs", which only the harmonic circuit '//&
            'solver assembles.')
        CALL Fatal(Caller,'Transient circuits do not support "Master BCs" for '//&
            'component '//I2S(Comp % ComponentId)//'!')
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE CheckTransientComponents
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Report the circuit model as it was actually resolved, once, after the
!> circuits have been read. Most of the ways a circuit definition goes wrong are
!> silent - a component whose "Master Bodies" match no body still assembles,
!> it just contributes nothing - so print what was resolved and warn about the
!> combinations that cannot be intended.
!------------------------------------------------------------------------------
  SUBROUTINE CircuitsSummary()
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(Component_t), POINTER :: Comp
    TYPE(ValueList_t), POINTER :: CompParams
    TYPE(Element_t), POINTER :: Element
    TYPE(Variable_t), POINTER :: Wvar
    INTEGER :: p, CompInd, t, n, nElem, nBC, nSupp, dim
    LOGICAL :: Parallel, HaveW, CheckSupport
    CHARACTER(LEN=16) :: cType
    CHARACTER(LEN=8)  :: cSupp
    CHARACTER(LEN=12) :: cRes
    CHARACTER(LEN=11) :: cOwner
    ! One set of widths for the header and the rows so the columns line up.
    CHARACTER(*), PARAMETER :: HdrFmt = '(A6,2X,A16,A7,A7,A8,3A14,2X,A12,A11)'
    CHARACTER(*), PARAMETER :: RowFmt = '(I6,2X,A16,I7,I7,A8,3ES14.4,2X,A12,A11)'
    CHARACTER(*), PARAMETER :: Caller = 'CircuitsSummary'

    ! This early exit has to be taken by every partition or by none, because
    ! the element counts below are reduced with a collective. It is: the output
    ! level mask is set once from the sif in ModelDescription and is identical on
    ! all partitions - the non-root ones are silenced inside Info() through
    ! OutputPE, not through the mask. Keep it that way.
    IF( .NOT. InfoActive(5) ) RETURN

    dim = CoordinateSystemDimension()
    ! Element counts are partial sums only when the mesh is really split.
    Parallel = CircuitsPartitionedMesh()

    ! Probe the potential exactly as HasSupport() does, so the column reports what
    ! the assembly will really do. With no such field HasSupport() does not gate
    ! at all, so there is no support count to show.
    CALL GetWPotentialVar(Wvar)
    HaveW = ASSOCIATED(Wvar)

    CALL Info(Caller,'Resolved circuit model:',Level=5)

    DO p=1,CurrentModel % n_Circuits
      Circuit => CurrentModel % Circuits(p)

      WRITE(Message,'(A,I0,A,I0,A,I0,A)') 'Circuit ',p,': ',Circuit % n, &
          ' variables, ',Circuit % n_comp,' components'
      CALL Info(Caller,Message,Level=5)

      IF( Circuit % n_comp == 0 ) CYCLE

      WRITE(Message,HdrFmt) 'comp','type            ','elems','bcs','Wsupp', &
          'ElArea','N_j','turns','resistance  ','  owner i/v'
      CALL Info(Caller,Message,Level=5)

      DO CompInd=1,Circuit % n_comp
        Comp => Circuit % Components(CompInd)

        cType = Comp % CoilType
        IF( LEN_TRIM(cType) == 0 ) cType = Comp % ComponentType

        ! Only these two coil types are gated by HasSupport() in the assembly.
        CheckSupport = ( dim == 3 .AND. &
            ( Comp % CoilType == 'massive' .OR. Comp % CoilType == 'foil winding' ) )

        ! Straight off the lists built by BuildComponentElementLists(), so the
        ! numbers reported here are by construction the ones the assembly uses.
        nElem = SIZE(Comp % ElemIdx)
        nBC = SIZE(Comp % BCElemIdx)
        nSupp = 0
        IF( CheckSupport .AND. HaveW ) THEN
          DO t=1,nElem
            Element => GetActiveElement(Comp % ElemIdx(t))
            n = GetElementNOFNodes(Element)
            IF( HasSupport(Element,n) ) nSupp = nSupp + 1
          END DO
        END IF

        IF( Parallel ) THEN
          nElem = NINT( ParallelReduction( 1.0_dp*nElem ) )
          nBC   = NINT( ParallelReduction( 1.0_dp*nBC ) )
          nSupp = NINT( ParallelReduction( 1.0_dp*nSupp ) )
        END IF

        ! UseCoilResistance is only decided during assembly, so report the input
        ! that decides it rather than the flag itself.
        CompParams => CurrentModel % Components(Comp % ComponentId) % Values
        IF( ListCheckPresent(CompParams,'Resistance') ) THEN
          cRes = 'keyword'
        ELSE IF( Comp % CoilType == 'massive' ) THEN
          ! A massive coil integrates a conductance, which CircuitsOutput then
          ! inverts to get something to report as a resistance.
          cRes = 'conductance'
        ELSE
          cRes = 'integrated'
        END IF

        ! Support is a property of the bulk elements, so with none of those there
        ! is nothing to report.
        IF( .NOT. CheckSupport .OR. nElem == 0 ) THEN
          cSupp = '       -'
        ELSE IF( .NOT. HaveW ) THEN
          cSupp = '     n/a'
        ELSE
          WRITE(cSupp,'(I8)') nSupp
        END IF

        WRITE(cOwner,'(I5,A,I0)') Comp % ivar % Owner,'/',Comp % vvar % Owner

        WRITE(Message,RowFmt) Comp % ComponentId, cType, nElem, nBC, cSupp, &
            Comp % ElArea, Comp % N_j, Comp % nofturns, cRes, cOwner
        CALL Info(Caller,Message,Level=5)

        ! A component that names bodies or boundaries but matches no element is
        ! always a mistake: the equations are assembled and couple to nothing.
        IF( nElem + nBC == 0 ) THEN
          IF( ASSOCIATED(Comp % BodyIds) .OR. ASSOCIATED(Comp % ElBoundaries) ) THEN
            CALL Warn(Caller,'Component '//I2S(Comp % ComponentId)//&
                ' is associated to bodies or boundaries but matches no element!')
          END IF
        ELSE IF( CheckSupport .AND. HaveW .AND. nElem > 0 .AND. nSupp == 0 ) THEN
          ! Only meaningful when there are bulk elements to have support in the
          ! first place - a component defined through "Master BCs" has none, and
          ! that is not an error.
          CALL Warn(Caller,'Component '//I2S(Comp % ComponentId)//' has '//&
              I2S(nElem)//' bulk elements but none with "W" support, so they will be skipped!')
        END IF

        IF( Comp % ComponentType /= 'resistor' .AND. nBC > 0 .AND. &
            .NOT. CurrentModel % HarmonicCircuits ) THEN
          CALL Warn(Caller,'Component '//I2S(Comp % ComponentId)//&
              ' has boundary elements, but the transient assembly only visits bulk elements!')
        END IF
      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE CircuitsSummary
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
    INTEGER :: i, nm, Circuit_tot_n, p, j, &
               RowId, nn, l, k, n_Circuits
    
    CM => CurrentModel % CircuitMatrix
    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('SetCircuitsParallelInfo','ASolver not found!')
    nm = ASolver % Matrix % NumberOfRows
    Circuit_tot_n = CurrentModel % Circuit_tot_n
    Circuits => CurrentModel % Circuits
    n_Circuits = CurrentModel % n_Circuits
    
    IF(.NOT.ASSOCIATED(CM % ParallelInfo)) THEN
      ALLOCATE(CM % ParallelInfo)
      ALLOCATE(CM % ParallelInfo % NeighbourList(nm+Circuit_tot_n))
      DO i=1,nm+Circuit_tot_n
        CM % ParallelInfo % NeighbourList(i) % Neighbours => Null()
      END DO
    END IF

    ! Every circuit row is shared by all partitions: the row is owned by
    ! Cvar % Owner and every other partition is listed as its neighbour. That is
    ! a safe over-estimate - a partition that never assembles into the row only
    ! contributes zeros - but it does make the circuit block dense in the
    ! communication pattern, so the cost grows with the number of partitions.
    !
    ! There used to be code here that tried to do better: for each circuit
    ! variable it swept the active elements with ElAssocToCvar(), counted the
    ! ones belonging to the variable's component, and did an MPI_ALLREDUCE to
    ! learn which partitions actually touch it. The result was then discarded by
    ! an unconditional "r_cnt = 1; nn = ParEnv % PEs ! for now" before it was
    ! ever read, so the neighbour lists have always been the all-partitions ones
    ! built below. The dead sweep has been removed - it cost
    ! n_circuit_variables element loops plus one collective each, per startup,
    ! for nothing.
    !
    ! To restore the intended behaviour, gather the per-variable element counts
    ! into one array and reduce them in a single collective, then let nn be the
    ! number of partitions with a non-zero count (forcing in the owner and the
    ! local partition, as the old code did) and filter the neighbour loops below
    ! on that count. ElAssocToCvar() is kept for that purpose.
    ! --------------------------------------------------------------------------
    nn = ParEnv % PEs

    DO p = 1,n_Circuits
      DO i=1,Circuits(p) % n
        Cvar => Circuits(p) % CircuitVariables(i)

        RowId = Cvar % ValueId + nm

        IF (Circuits(p) % Harmonic) THEN
          DO j=1,Cvar % Dofs
            IF(.NOT.ASSOCIATED(CM % ParallelInfo % NeighbourList(RowId+AddIndex(j-1))%Neighbours)) THEN
              ALLOCATE(CM % ParallelInfo % NeighbourList(RowId+AddIndex(j-1)) % Neighbours(nn))
              ALLOCATE(CM % ParallelInfo % NeighbourList(RowId+AddImIndex(j-1)) % Neighbours(nn))
            END IF
            CM % ParallelInfo % NeighbourList(RowId+AddIndex(j-1)) % Neighbours(1)   = CVar % Owner
            CM % ParallelInfo % NeighbourList(RowId+AddImIndex(j-1)) % Neighbours(1) = Cvar % Owner
            l = 1
            DO k=0,ParEnv % PEs-1
              IF(k==CVar % Owner) CYCLE
              l = l + 1
              CM % ParallelInfo % NeighbourList(RowId+AddIndex(j-1)) % Neighbours(l) = k
              CM % ParallelInfo % NeighbourList(RowId+AddImIndex(j-1)) % Neighbours(l) = k
            END DO
            CM % RowOwner(RowId + AddIndex(j-1))   = Cvar % Owner
            CM % RowOwner(RowId + AddImIndex(j-1)) = Cvar % Owner
          END DO
        ELSE
          DO j=1,Cvar % Dofs
            IF(.NOT.ASSOCIATED(CM % ParallelInfo % NeighbourList(RowId+j-1)%Neighbours)) THEN
              ALLOCATE(CM % ParallelInfo % NeighbourList(RowId+j-1) % Neighbours(nn))
            END IF
            CM % ParallelInfo % NeighbourList(RowId+j-1) % Neighbours(1) = CVar % Owner
            l = 1
            DO k=0,ParEnv % PEs-1
              IF(k==CVar % Owner) CYCLE
              l = l + 1
              CM % ParallelInfo % NeighbourList(RowId+j-1) % Neighbours(l) = k
            END DO
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
    LOGICAL :: Parallel
    INTEGER, POINTER :: Rows(:), Cnts(:)
    
    Circuits => CurrentModel % Circuits
    n_Circuits = CurrentModel % n_Circuits
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    Parallel = CurrentModel % Solver % Parallel
    
    ! Basic circuit equations...
    ! ---------------------------
    DO p = 1,n_Circuits
      DO i=1,Circuits(p) % n
        Cvar => Circuits(p) % CircuitVariables(i)
        ! Guarded on Parallel exactly as AddBasicCircuitEquations() is. Without
        ! that guard a run with PEs > 1 but a non-parallel linear system - a
        ! replicated mesh, i.e. slices or parallel timestepping - created these
        ! rows only on the owning partition while every partition went on to
        ! assemble into them.
        IF( Parallel ) THEN
          ! Guarded on Parallel exactly as AddBasicCircuitEquations() is. Without
        ! that guard a run with PEs > 1 but a non-parallel linear system - a
        ! replicated mesh, i.e. slices or parallel timestepping - created these
        ! rows only on the owning partition while every partition went on to
        ! assemble into them.
        IF( Parallel ) THEN
          IF(Cvar % Owner /= ParEnv % myPE) CYCLE
        END IF
        END IF

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
    LOGICAL :: Parallel
    INTEGER, POINTER :: Rows(:), Cols(:), Cnts(:)
    
    Circuits => CurrentModel % Circuits
    n_Circuits = CurrentModel % n_Circuits
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    Parallel = CurrentModel % Solver % Parallel
    
    ! Basic circuit equations...
    ! ---------------------------
    DO p = 1,n_Circuits
      DO i=1,Circuits(p) % n
        Cvar => Circuits(p) % CircuitVariables(i)
        ! Guarded on Parallel exactly as AddBasicCircuitEquations() is. Without
        ! that guard a run with PEs > 1 but a non-parallel linear system - a
        ! replicated mesh, i.e. slices or parallel timestepping - created these
        ! rows only on the owning partition while every partition went on to
        ! assemble into them.
        IF( Parallel ) THEN
          IF(Cvar % Owner /= ParEnv % myPE) CYCLE
        END IF

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
               CompInd, q, qi
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
        IF (Comp % ComponentType == 'resistor') THEN
            CALL CountMatElement(Rows, Cnts, RowId, 1)
            CALL CountMatElement(Rows, Cnts, RowId, 1)
            CYCLE
        ELSE
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
        END IF

!        temp = SUM(Cnts)
!print *, "Active elements", ParEnv % Mype, ":", GetNOFActive()
        ! The component's own element lists, walked backwards, so the structure
        ! is created in exactly the order the assembly fills it. See
        ! BuildComponentElementLists().
        DO qi=SIZE(Comp % ElemIdx),1,-1
          Element => GetActiveElement(Comp % ElemIdx(qi))
          CALL CountComponentElements(Element, Comp, RowId, Rows, Cnts, Done, dofsdone)
        END DO

        DO qi=SIZE(Comp % BCElemIdx),1,-1
          Element => GetBoundaryElement(Comp % BCElemIdx(qi))
          CALL CountComponentElements(Element, Comp, RowId, Rows, Cnts, Done, dofsdone)
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
               CompInd, q, qi
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

        IF (Comp % ComponentType == 'resistor') THEN
            CALL CreateMatElement(Rows, Cols, Cnts, VvarId, IvarId)
            CALL CreateMatElement(Rows, Cols, Cnts, VvarId, VvarId)
            CYCLE
        ELSE
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
                    CALL CreateMatElement(Rows, Cols, Cnts, VvarId + AddIndex(j), VvarId + AddIndex(jj))
                END DO
              END IF
            END DO
          END SELECT
        END IF

!        temp = SUM(Cnts)
!print *, "Active elements ", ParEnv % Mype, ":", GetNOFActive()
        ! As in CountComponentEquations, and it has to match it exactly.
        DO qi=SIZE(Comp % ElemIdx),1,-1
          Element => GetActiveElement(Comp % ElemIdx(qi))
          CALL CreateComponentElements(Element, Comp, VvarId, IvarId, Rows, Cols, Cnts, Done, dofsdone)
        END DO

        DO qi=SIZE(Comp % BCElemIdx),1,-1
          Element => GetBoundaryElement(Comp % BCElemIdx(qi))
          CALL CreateComponentElements(Element, Comp, VvarId, IvarId, Rows, Cols, Cnts, Done, dofsdone)
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
   SUBROUTINE CountComponentElements(Element, Comp, RowId, Rows, Cnts, Done, dofsdone)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    TYPE(Component_t), POINTER :: Comp
    INTEGER :: nn, nd, RowId
    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: Rows(:), Cnts(:)
    LOGICAL*1 :: Done(:)
    LOGICAL :: dofsdone

    IF (ElAssocToComp(Element, Comp)) THEN
      Asolver => CurrentModel % Asolver
      nn = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element,ASolver)
      SELECT CASE (Comp % CoilType)
      CASE('stranded')           
        CALL CountAndCreateStranded(Element,nn,nd,RowId,Cnts,Done,Rows)
      CASE('massive')
        IF (HasSupport(Element,nn)) THEN
          CALL CountAndCreateMassive(Element,nn,nd,RowId,Cnts,Done,Rows)
        END IF
     CASE('foil winding')
        IF (HasSupport(Element,nn)) THEN
          CALL CountAndCreateFoilWinding(Element,nn,nd,Comp,Cnts,Done,Rows)
        END IF
      END SELECT
    END IF
!------------------------------------------------------------------------------
   END SUBROUTINE CountComponentElements
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CreateComponentElements(Element, Comp, VvarId, IvarId, Rows, Cols, Cnts, Done, dofsdone)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    TYPE(Component_t), POINTER :: Comp
    TYPE(Solver_t), POINTER :: ASolver
    INTEGER :: nn, nd, VvarId, IvarId
    INTEGER, POINTER :: Rows(:), Cols(:), Cnts(:)
    LOGICAL*1 :: Done(:)
    LOGICAL :: dofsdone
    
    IF (ElAssocToComp(Element, Comp)) THEN
      Asolver => CurrentModel % Asolver
      nn = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element,ASolver)
      SELECT CASE (Comp % CoilType)
      CASE('stranded')
        CALL CountAndCreateStranded(Element,nn,nd,VvarId,Cnts,Done,Rows,Cols,IvarId)
      CASE('massive')
        IF (HasSupport(Element,nn)) THEN
          CALL CountAndCreateMassive(Element,nn,nd,VvarId,Cnts,Done,Rows,Cols=Cols)
        END IF
     CASE('foil winding')
        IF (HasSupport(Element,nn)) THEN
          CALL CountAndCreateFoilWinding(Element,nn,nd,Comp,Cnts,Done,Rows,Cols=Cols)
        END IF
      END SELECT
    END IF
!------------------------------------------------------------------------------
   END SUBROUTINE CreateComponentElements
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CountAndCreateStranded(Element,nn,nd,i,Cnts,Done,Rows,Cols,Jsind,Harmonic)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    INTEGER :: nn, nd, ncdofs1, ncdofs2, dim
    OPTIONAL :: Cols
    INTEGER :: Rows(:), Cols(:), Cnts(:)
    INTEGER :: p,i,j,k,Indexes(nd)
    INTEGER, OPTIONAL :: Jsind
    INTEGER, POINTER :: PS(:)
    LOGICAL*1 :: Done(:)
    INTEGER :: MyGen = -1
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm
    SAVE dim, MyGen

    IF (MyGen /= CircuitsGeneration) THEN
      MyGen = CircuitsGeneration
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
      ncdofs1=nn+1
      ncdofs2=nd
    END IF

    DO p=ncdofs1,ncdofs2
      j = Indexes(p)

      IF( ASSOCIATED( CurrentModel % Mesh % PeriodicPerm ) ) THEN
        ! If we have periodicity eliminated only flag the master in Done
        k = CurrentModel % Mesh % PeriodicPerm(j)
        IF( k > 0 ) j = k
      END IF

      IF(.NOT.Done(j)) THEN
        Done(j) = .TRUE.
        j = PS(j)
        IF (harm) j = ReIndex(j)
        IF(PRESENT(Cols)) THEN
          CALL CreateMatElement(Rows, Cols, Cnts, i, j, harm) 
          CALL CreateMatElement(Rows, Cols, Cnts, j, Jsind, harm)
!         CALL CreateMatElement(Rows, Cols, Cnts, j, Jsind)
        ELSE
          CALL CountMatElement(Rows, Cnts, i, 1, harm)
          CALL CountMatElement(Rows, Cnts, j, 1, harm)
!         CALL CountMatElement(Rows, Cnts, j, 1)
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
    INTEGER :: p,i,j,k,Indexes(nd)
    INTEGER, POINTER :: PS(:)
    LOGICAL*1 :: Done(:)
    INTEGER :: MyGen = -1
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm
    SAVE dim, MyGen

    IF (MyGen /= CircuitsGeneration) THEN
      MyGen = CircuitsGeneration
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
      ncdofs1=nn+1
      ncdofs2=nd
    END IF
    DO p=ncdofs1,ncdofs2
      j = Indexes(p)

      IF( ASSOCIATED( CurrentModel % Mesh % PeriodicPerm ) ) THEN
        ! If we have periodicity eliminated only flag the master in Done
        k = CurrentModel % Mesh % PeriodicPerm(j)
        IF( k > 0 ) j = k
      END IF

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
    LOGICAL :: dofsdone
    INTEGER :: MyGen = -1
    INTEGER, POINTER :: PS(:)
    LOGICAL*1 :: Done(:)
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm
    SAVE dim, MyGen

    IF (MyGen /= CircuitsGeneration) THEN
      MyGen = CircuitsGeneration
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
    LOGICAL :: Parallel, Found
    
    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Circuits_MatrixInit','ASolver not found!')
    Circuit_tot_n = CurrentModel % Circuit_tot_n
    
    ! Initialize Circuit matrix:
    ! -----------------------------
    PS => Asolver % Variable % Perm
    nm =  Asolver % Matrix % NumberOfRows

    CM => AllocateMatrix()
    CurrentModel % CircuitMatrix=>CM
    
    CM % Format = MATRIX_CRS
    Asolver % Matrix % AddMatrix => CM
    ALLOCATE(CM % RHS(nm + Circuit_tot_n)); CM % RHS=0._dp

    CM % NumberOfRows = nm + Circuit_tot_n

    n = CM % NumberOfRows

    ALLOCATE(Rows(n+1), Cnts(n)); Rows=0; Cnts=0
    ALLOCATE(Done(SIZE(PS)), CM % RowOwner(n)); Cm % RowOwner=-1


    Parallel = CurrentModel % Solver % Parallel      
    IF( Parallel ) CALL SetCircuitsParallelInfo()

    ! COUNT SIZES:
    ! ============
    dofsdone = .FALSE.
    
    CALL CountBasicCircuitEquations(Rows, Cnts)
    CALL CountComponentEquations(Rows, Cnts, Done, dofsdone)

    ! ALLOCATE CRS STRUCTURES (if need be):
    ! =====================================

    n = SUM(Cnts)

    IF (n<=0) THEN
      CM % NumberOfRows = 0
      CALL CircuitsRecordBuild()
      DEALLOCATE(Rows,Cnts,Done)
      ! FreeMatrix rather than a bare DEALLOCATE of CM: RHS, RowOwner and, in
      ! parallel, ParallelInfo % NeighbourList were all allocated above and were
      ! being leaked with it.
      CALL FreeMatrix(CM); CM=>Null()
      Asolver %  Matrix % AddMatrix => CM
      CurrentModel % CircuitMatrix=>CM
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

    ! CREATE COLUMNS:
    ! ===============

    CALL CreateBasicCircuitEquations(Rows, Cols, Cnts)
    CALL CreateComponentEquations(Rows, Cols, Cnts, Done, dofsdone)
    
    IF (n /= SUM(Cnts)) THEN
      CALL Fatal('Circuits_MatrixInit', &
                 'Inconsistent number of matrix elements: '//I2S(n)//' vs. '//I2S(SUM(CNTs)))
    END IF

    DEALLOCATE( Cnts, Done )
    CM % Rows => Rows
    CM % Cols => Cols
    CM % Values => Values
    CALL CRS_SortMatrix(CM)
    
    Asolver %  Matrix % AddMatrix => CM

    CALL CircuitsRecordBuild()
!------------------------------------------------------------------------------
  END SUBROUTINE Circuits_MatrixInit
!------------------------------------------------------------------------------

      
END MODULE CircMatInitMod


