!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juhani Kataja
! *  Email:   Juhani.Kataja@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 12 Feb 2026
! *
! *****************************************************************************/


!-------------------------------------------------------------------------------
!>  Utilities for Adios2OutputSolver
!-------------------------------------------------------------------------------
MODULE AdiosOutputSolverUtils

USE ADIOS2Utils

IMPLICIT NONE

!> Per-field renumbering of the saved region. 
!> The permutations mirror those used by VtuOutputSolver and
!> are produced by GenerateSaveMask / GenerateSavePermutation.
TYPE FieldCache_t
  CHARACTER(len=:), ALLOCATABLE :: name    ! scalar field name / vector base name
  INTEGER, ALLOCATABLE :: NodePerm(:)      ! mesh node -> output point (0 if not saved)
  INTEGER, ALLOCATABLE :: InvNodePerm(:)   ! output point -> mesh node (coordinates)
  INTEGER, ALLOCATABLE :: DgPerm(:)        ! DG index -> output point (unused for nodal)
  INTEGER, ALLOCATABLE :: InvDgPerm(:)     ! output point -> DG index (unused for nodal)
  INTEGER, ALLOCATABLE :: GatherIdx(:)     ! output point -> index into var % values (0 => 0)
  LOGICAL, ALLOCATABLE :: ActiveElem(:)    ! per element in mesh (bulk+boundary)
  INTEGER :: NumberOfDofNodes = 0          ! number of output points
  INTEGER :: NumberOfElements = 0          ! number of active elements
  INTEGER :: ElemFirst = 0, ElemLast = 0
  LOGICAL :: NoPermutation = .FALSE.       ! .true. => output points are all mesh nodes
END TYPE

TYPE AdiosOutput_t
  TYPE(AdiosWriter_t) :: writer
  TYPE(FieldCache_t), ALLOCATABLE :: caches(:)
END TYPE

CONTAINS

!-------------------------------------------------------------------------------
!> Return the cache for a given field, building it once on first use.
!-------------------------------------------------------------------------------
FUNCTION GetFieldCache(holder, name, var, params) RESULT(c)
  USE DefUtils
  TYPE(AdiosOutput_t), POINTER :: holder
  CHARACTER(*), INTENT(IN) :: name
  TYPE(Variable_t), POINTER :: var
  TYPE(ValueList_t), POINTER :: params
  TYPE(FieldCache_t), POINTER :: c

  TYPE(FieldCache_t), ALLOCATABLE :: tmp(:)
  INTEGER :: i, nold

  c => NULL()
  IF (ALLOCATED(holder % caches)) THEN
    DO i = 1, SIZE(holder % caches)
      IF (holder % caches(i) % name == name) THEN
        c => holder % caches(i)
        RETURN
      END IF
    END DO
    nold = SIZE(holder % caches)
    CALL MOVE_ALLOC(holder % caches, tmp)
    ALLOCATE(holder % caches(nold+1))
    holder % caches(1:nold) = tmp(1:nold)
    DEALLOCATE(tmp)
  ELSE
    nold = 0
    ALLOCATE(holder % caches(1))
  END IF

  c => holder % caches(nold+1)
  c % name = name
  CALL BuildFieldCache(c, var, params)
END FUNCTION GetFieldCache

!-------------------------------------------------------------------------------
!> Populate a cache using the SaveUtils machinery. The saved region is
!> taken from an explicit mask on the solver (Mask Name / Mask Condition / Mask
!> Variable); if none is given it defaults to where this field is defined, i.e.
!> the variable's own perm. The resulting node/element permutations are identical
!> to those VtuOutputSolver builds.
!-------------------------------------------------------------------------------
SUBROUTINE BuildFieldCache(c, var, params)
  USE DefUtils
  USE SaveUtils
  TYPE(FieldCache_t) :: c
  TYPE(Variable_t), POINTER :: var
  TYPE(ValueList_t), POINTER :: params

  TYPE(Mesh_t), POINTER :: msh
  INTEGER :: NumberOfGeomNodes, ii, i, j
  LOGICAL :: parallel, injected, has_mask

  msh => var % primarymesh
  parallel = ( ParEnv % PEs > 1 )

  ! Default the saved region to where this field is defined (its perm) unless the
  ! user configured an explicit mask on the solver.
  has_mask = ListCheckPresent(params,'Mask Variable') .OR. &
             ListCheckPresent(params,'Mask Name') .OR. &
             ListCheckPresent(params,'2D Mask Name') .OR. &
             ListCheckPresent(params,'3D Mask Name') .OR. &
             ListCheckPresent(params,'Mask Condition')
  injected = .FALSE.
  IF (.NOT. has_mask .AND. ASSOCIATED(var % perm)) THEN
    CALL ListAddString(params,'Mask Variable', TRIM(var % Name))
    injected = .TRUE.
  END IF

  CALL GenerateSaveMask(msh, params, parallel, 0, .FALSE., &
      c % NodePerm, c % ActiveElem, NumberOfGeomNodes, c % NumberOfElements, &
      c % ElemFirst, c % ElemLast)

  ! Nodal only for now (DG=.false., DN=.false., LagN=0, SaveLinear=.false.); the
  ! DgPerm/InvDgPerm outputs are left unallocated until DG support is wired in.
  CALL GenerateSavePermutation(msh, .FALSE., .FALSE., 0, .FALSE., c % ActiveElem, &
      NumberOfGeomNodes, c % NoPermutation, c % NumberOfDofNodes, &
      c % DgPerm, c % InvDgPerm, c % NodePerm, c % InvNodePerm)

  IF (injected) CALL ListRemove(params,'Mask Variable')

  ! Precompute the per-point gather index into var % values (0 => write value 0,
  ! mirroring VtuOutputSolver for points inside the saved region but outside the
  ! field's own support).
  ALLOCATE(c % GatherIdx(c % NumberOfDofNodes))
  DO ii = 1, c % NumberOfDofNodes
    IF (c % NoPermutation) THEN
      i = ii
    ELSE
      i = c % InvNodePerm(ii)
    END IF
    IF (ASSOCIATED(var % perm)) THEN
      j = 0
      IF (i >= 1 .AND. i <= SIZE(var % perm)) j = var % perm(i)
    ELSE
      j = i
    END IF
    c % GatherIdx(ii) = j
  END DO
END SUBROUTINE BuildFieldCache

!-------------------------------------------------------------------------------
!> Gather a scalar (or vector component) over the cached output points.
!-------------------------------------------------------------------------------
FUNCTION GatherField(var, c) RESULT(vals)
  USE DefUtils
  TYPE(Variable_t), POINTER :: var
  TYPE(FieldCache_t) :: c
  REAL(KIND=dp) :: vals(c % NumberOfDofNodes)
  INTEGER :: ii, j

  DO ii = 1, c % NumberOfDofNodes
    j = c % GatherIdx(ii)
    IF (j > 0) THEN
      vals(ii) = var % values(j)
    ELSE
      vals(ii) = 0.0_dp
    END IF
  END DO
END FUNCTION GatherField

SUBROUTINE GetAdiosHolder(Solver, holder, found)
  USE DefUtils
  USE ADIOS2Utils
  USE iso_c_binding

  IMPLICIT NONE

  TYPE(Solver_t) :: Solver
  Type(AdiosOutput_t), POINTER, intent(out) :: holder
  LOGICAL :: Found
  INTEGER(Kind=AddrInt) :: WriterPtr
  TYPE(ValueList_t), POINTER :: params

  params => GetSolverParams(Solver)

  writerPtr = ListGetAddressInteger(Params, 'ADIOS2 Writer Ptr', Found)
  if (.not. found) return

  call C_F_Pointer(TRANSFER(WriterPtr, C_NULL_PTR), holder)

END SUBROUTINE

SUBROUTINE MakeFidesJson(Variables, params, holder, fname)
  USE DefUtils
  type(Variable_t), POINTER :: Variables
  TYPE(ValueList_T), POINTER :: params
  type(AdiosOutput_t) :: holder
  character(*), intent(in) :: fname
  character(len=:), allocatable :: pre, fields_pre, fields_post, post
  character :: NL
  character(len=:), allocatable :: fields

  NL = NEW_LINE('a')
  
  pre = ' {' // NL // &
  &'    "unstructured_grid": {' // NL // &
  &'        "data_sources": [' // NL // &
  &'            {' // NL // &
  &'                "name": "source",' // NL // &
  &'                "filename_mode": "relative",' // NL // &
  &'                "filename": "."' // NL // &
  &'            }' // NL // &
  &'        ],' // NL // &
  &'        "step_information": { "data_source": "source" }, ' // NL // &
  &'        "coordinate_system": {' // NL // &
  &'          "array": {' // NL // &
  &'            "array_type": "composite",' // NL // &
  &'              "x_array": {' // NL // &
  &'                "array_type": "basic", ' // NL // &
  &'                "data_source": "source", ' // NL // &
  &'                "variable": "points_x"' // NL // &
  &'              },' // NL // &
  &'              "y_array": {' // NL // &
  &'                "array_type": "basic", ' // NL // &
  &'                "data_source": "source", ' // NL // &
  &'                "variable": "points_y"' // NL // &
  &'            },' // NL // &
  &'              "z_array": {' // NL // &
  &'                "array_type": "basic", ' // NL // &
  &'                "data_source": "source", ' // NL // &
  &'                "variable": "points_z"' // NL // &
  &'              }' // NL // &
  &'          }' // NL // &
  &'        },' // NL // &
  &'        "cell_set": {' // NL // &
  &'            "cell_set_type": "explicit",' // NL // &
  &'            "connectivity": {' // NL // &
  &'                "array_type": "basic",' // NL // &
  &'                "data_source": "source",' // NL // &
  &'                "variable": "connectivity"' // NL // &
  &'            },' // NL // &
  &'            "cell_types": {' // NL // &
  &'                "array_type": "basic",' // NL // &
  &'                "data_source": "source",' // NL // &
  &'                "variable": "cell_types"' // NL // &
  &'            },' // NL // &
  &'            "number_of_vertices": {' // NL // &
  &'                "array_type": "basic",' // NL // &
  &'                "data_source": "source",' // NL // &
  &'                "variable": "num_verts"' // NL // &
  &'            }' // NL // &
  &'        },' // NL // &
  &'' // NL // &
  &'        "fields": ['

  post = NL // &
  &'        ]' // NL // &
  &'' // NL // &
  &'    } }' // NL 

  open(unit=10, file=fname, status='replace', action='write')
  write(10,*) pre

  block

    CHARACTER(len=:), ALLOCATABLE :: scalar_name, vector_name
    logical :: found_field_i, found_vfield_i
    integer ::  vector_var_ind, scalar_var_ind, round
    type(variable_t), pointer :: variable

    scalar_var_ind = 1
    vector_var_ind = 1
    round = 1
    var_ind_do: do while (.true.)
      scalar_name = trim(ListGetString(params, 'Scalar Field '//i2s(scalar_var_ind), found_field_i))
      if(found_field_i) then
        scalar_var_ind = scalar_var_ind + 1
        if(round>1) write(10,'(a)') ','
        round = round + 1
        write(10,'(a)', advance='no') trim(make_field(scalar_name, NL))
      end if

      vector_name = trim(ListGetString(params, 'Vector Field '//i2s(vector_var_ind), found_vfield_i))
      if(found_vfield_i) then
        vector_var_ind = vector_var_ind + 1
        if(round>1) write(10,'(a)') ','
        round = round + 1
        write(10,'(a)', advance='no') trim(make_field(vector_name, NL))
      end if
      if (.not. (found_vfield_i .or. found_field_i)) exit var_ind_do
    end do var_ind_do

  end block

  write(10,*) post
  close(10)

END SUBROUTINE

  function make_field(fieldname, NL) result(S)
    character(len=*), intent(in) :: fieldname
    character(len=8192) :: S
    character, intent(in) :: NL

    write(S, '(a,a,a,a,a)') '          {' // NL // &
&'            "name": "', fieldname, '",' // NL // &
&'            "association":"points",' // NL // &
&'            "array": {' // NL // &
&'              "array_type":"basic",' // NL // &
&'              "data_source": "source",' // NL // &
&'              "variable": "', fieldname, '"'//NL // &
&'            }' // NL // &
&'          }' 
  end function
END MODULE

!------------------------------------------------------------------------------
!> ADIOS2OutputSolver_Init initializes ADIOS2OutputSolver
!> Makes Fides json-file and sets up adios2 io-object AdiosWriter_t
!> 
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE ADIOS2OutputSolver_Init(Model, Solver, dt, TransientSimulation)

  USE DefUtils
  USE ADIOS2Utils
  USE iso_c_binding
  use AdiosOutputSolverUtils

  IMPLICIT NONE

  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

  Type(AdiosWriter_t), POINTER :: Writer
  TYPE(AdiosOutput_t), POINTER :: output_holder

  INTEGER(Kind=AddrInt) :: WriterPtr
  integer :: ierr
  LOGICAL :: Found
  TYPE(ValueList_t), POINTER :: params
  type(variable_t), pointer :: variable
  character(:), ALLOCATABLE :: output_fname

  output_holder => Null()

  params => ListGetSolverParams()
  
  CALL GetAdiosHolder(Solver, output_holder, found)

  if(found) return

  ALLOCATE(output_holder)
  writer => output_holder % writer 

  output_fname = trim(ListGetString(params, 'Output File name', UnfoundFatal=.true.))

  ierr = writer % init(output_fname, array_kind=ADIOS2_ARRAY_GLOBAL)
  call MakeFidesJson(model % variables, params, output_holder, output_fname // '/elmer_fields.json')
  writerPtr = TRANSFER(C_LOC(output_holder), WriterPtr)
  call ListAddAddressInteger(params, 'ADIOS2 Writer Ptr', WriterPtr)


END SUBROUTINE ADIOS2OutputSolver_Init

!------------------------------------------------------------------------------
!> ADIOS2OutputSolver main routine that checks if mesh needs to be updated each time step
!> Saves single mesh for all variables.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE ADIOS2OutputSolver(Model, Solver, dt, TransientSimulation)

  USE DefUtils
  USE ADIOS2Utils
  USE iso_c_binding
  use AdiosOutputSolverUtils
  USE SaveUtils

  IMPLICIT NONE

  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

  Type(AdiosWriter_t), POINTER :: Writer
  Type(AdiosOutput_t), POINTER :: output_holder
  INTEGER(Kind=AddrInt) :: WriterPtr
  integer :: ierr
  LOGICAL :: Found, found_field_i, save_mesh, found_save_mesh

  TYPE(ValueList_t), POINTER :: params
  type(variable_t), pointer :: variable

  integer :: nnodes, file_size, e_ind, n_ind, state, round, var_ind
  type(element_t), pointer :: elem
  integer(kind=4), allocatable :: elem_types(:), num_elem_nodes(:)
  integer(kind=4), allocatable :: connectivity(:)
  real(kind=dp), ALLOCATABLE :: debug_arr(:,:)
  integer :: lcon
  CHARACTER(:), ALLOCATABLE :: field_name

  TYPE(FieldCache_t), POINTER :: cache

  Writer => Null()

  params => ListGetSolverParams()
  
  CALL GetAdiosHolder(Solver, output_holder, found)

  IF(.not. Found) then
    call Fatal('ADIOS2OutputSolver', 'Writer ptr not found!')
  end if

  writer => output_holder % writer


  call writer % begin_step()

  ! Here loop over variables
  ! TODO This is fragile way to treat number of scalar fields since the field names must be
  !       continuously numbered starting from 1
  ! TODO: If mesh is saved twice on a timestep, then it will be overwritten

  save_mesh = .true.
  var_ind = 1
  var_ind_do: do while (var_ind > 0)
    field_name = trim(ListGetString(params, 'Scalar Field '//i2s(var_ind), found_field_i))
    if(.not. found_field_i) exit var_ind_do

    var_ind = var_ind + 1

    variable => VariableGet(model % variables, trim(field_name), DoInterp=.false., UnfoundFatal = .true.)

    if (variable % TYPE /= Variable_on_nodes) &
        call Fatal('ADIOS2OutputSolver', 'Only nodal fields are supported; scalar field "'// &
            trim(field_name)//'" has non-nodal TYPE '//i2s(variable % TYPE))

    cache => GetFieldCache(output_holder, field_name, variable, params)

    call writer % write_data(trim(field_name), GatherField(variable, cache))

    save_mesh = ListGetLogical(params, field_name // ' save mesh', found_save_mesh, defvalue = save_mesh)
    if(.not. found_save_mesh) call ListAddLogical(params, field_name // ' save mesh', .false.)

    if(save_mesh) then
      call LocalSaveMesh(cache)
      save_mesh = .false.
    end if
  end do var_ind_do

  block
    type(variable_t), pointer :: Vx, Vy, Vz
    real(kind=dp), allocatable :: V(:,:)
    logical :: dim
    ! treat vector variables
    var_ind = 1
    var_ind_do: do while(.true.)
      field_name = trim(ListGetString(params, 'Vector Field '//i2s(var_ind), found_field_i))
      if(.not. found_field_i) exit var_ind_do

      var_ind = var_ind + 1
      Vx => VariableGet(model % variables, trim(field_name) //" 1", DoInterp=.false., UnfoundFatal=.true.)
      Vy => VariableGet(model % variables, trim(field_name) //" 2", DoInterp=.false., UnfoundFatal=.true.)
      Vz => VariableGet(model % variables, trim(field_name) //" 3", DoInterp=.false., UnfoundFatal=.false.)

      if (Vx % TYPE /= Variable_on_nodes .or. Vy % TYPE /= Variable_on_nodes) &
          call Fatal('ADIOS2OutputSolver', 'Only nodal fields are supported; vector field "'// &
              trim(field_name)//'" is non-nodal')
      if (associated(Vz)) then
        if (Vz % TYPE /= Variable_on_nodes) &
            call Fatal('ADIOS2OutputSolver', 'Only nodal fields are supported; vector field "'// &
                trim(field_name)//'" is non-nodal')
      end if

      variable => Vx
      cache => GetFieldCache(output_holder, field_name, variable, params)

      if(allocated(V)) deallocate(V)
      if (associated(Vz)) then
        allocate(V(3,cache % NumberOfDofNodes))
        V(1,:) = GatherField(Vx, cache)
        V(2,:) = GatherField(Vy, cache)
        V(3,:) = GatherField(Vz, cache)
      else
        allocate(V(2,cache % NumberOfDofNodes))
        V(1,:) = GatherField(Vx, cache)
        V(2,:) = GatherField(Vy, cache)
      end if

      call writer % write_data(trim(field_name), V(:,:))

      save_mesh = ListGetLogical(params, field_name // ' save mesh', found_save_mesh, defvalue = save_mesh)
      if(.not. found_save_mesh) call ListAddLogical(params, field_name // ' save mesh', .false.)

      if(save_mesh) then
        call LocalSaveMesh(cache)
        save_mesh = .false.
      end if

    end do var_ind_do
  end block

  call writer % end_step()
contains

  subroutine LocalSaveMesh(cache)
    implicit none
    type(FieldCache_t) :: cache
    type(mesh_t), pointer :: msh
    type(element_t), pointer :: elem
    integer :: ii, i, jj, ncon, ea
    integer, allocatable :: TmpIndexes(:)
    real(kind=dp), allocatable :: px(:), py(:), pz(:)

    msh => variable % primarymesh

    ! Coordinates: one point per output dof node, mirroring VtuOutputSolver.
    allocate(px(cache % NumberOfDofNodes), py(cache % NumberOfDofNodes), pz(cache % NumberOfDofNodes))
    do ii = 1, cache % NumberOfDofNodes
      if (cache % NoPermutation) then
        i = ii
      else
        i = cache % InvNodePerm(ii)
      end if
      px(ii) = msh % nodes % x(i)
      py(ii) = msh % nodes % y(i)
      pz(ii) = msh % nodes % z(i)
    end do

    ! Connectivity: VTK-ordered node indexes remapped to output points.
    allocate(TmpIndexes(msh % MaxElementDOFs))
    ncon = 0
    do i = cache % ElemFirst, cache % ElemLast
      if (.not. cache % ActiveElem(i)) cycle
      ncon = ncon + msh % elements(i) % type % NumberOfNodes
    end do

    allocate(connectivity(ncon))
    allocate(elem_types(cache % NumberOfElements))
    allocate(num_elem_nodes(cache % NumberOfElements))

    lcon = 0
    ea = 0
    do i = cache % ElemFirst, cache % ElemLast
      if (.not. cache % ActiveElem(i)) cycle
      ea = ea + 1
      elem => msh % elements(i)
      call Elmer2VtkIndexes(elem, .false., .false., TmpIndexes)
      nnodes = elem % type % NumberOfNodes
      do n_ind = 1, nnodes
        if (cache % NoPermutation) then
          jj = TmpIndexes(n_ind)
        else
          jj = cache % NodePerm(TmpIndexes(n_ind))
        end if
        connectivity(lcon + n_ind) = jj - 1
      end do
      lcon = lcon + nnodes
      elem_types(ea) = Elmer2VtkElement(elem % type % elementcode, .false.)
      num_elem_nodes(ea) = nnodes
    end do

    call writer % write_data('num_verts', num_elem_nodes)
    call writer % write_data('cell_types', elem_types)
    call writer % write_data('connectivity', connectivity)
    call writer % write_data('points_x', px)
    call writer % write_data('points_y', py)
    call writer % write_data('points_z', pz)
  end subroutine
END SUBROUTINE

!------------------------------------------------------------------------------
!> ADIOS2OutputSolver_Finalize Finalizes ADIOS2OutputSolver
!> Finalizes the AdiosWriter_t object and deallocates structs 
!> that Elmer library is unaware
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE ADIOS2OutputSolver_Finalize(Model, Solver, dt, TransientSimulation)
  USE DefUtils
  USE ADIOS2Utils
  USE iso_c_binding
  use AdiosOutputSolverUtils

  IMPLICIT NONE

  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation, Found
  TYPE(AdiosWriter_t), POINTER :: Writer
  TYPE(AdiosOutput_t), POINTER :: output_holder
  integer :: ierr

  CALL GetAdiosHolder(Solver, output_holder, Found)

  if(Found) then
    ierr = output_holder % writer % finalize()
    deallocate(output_holder)
  end if

END SUBROUTINE
