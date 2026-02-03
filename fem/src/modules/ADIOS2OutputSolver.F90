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

TYPE AdiosOutput_t
  TYPE(AdiosWriter_t) :: writer
END TYPE

CONTAINS

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
  integer(kind=4), allocatable :: elem_types(:), offsets(:), num_elem_nodes(:)
  integer(kind=4), allocatable :: connectivity(:)
  real(kind=dp), ALLOCATABLE :: debug_arr(:,:)
  integer :: lcon 
  CHARACTER(:), ALLOCATABLE :: field_name

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

    variable => VariableGet(model % variables, trim(field_name), DoInterp=.false., UnfoundFatal = .false.)

    if (associated(variable%perm)) then
      call writer % write_data(trim(field_name), variable % values(variable%perm))
    else
      call writer % write_data(trim(field_name), variable % values(:))
    end if

    save_mesh = ListGetLogical(params, field_name // ' save mesh', found_save_mesh, defvalue = save_mesh)
    if(.not. found_save_mesh) call ListAddLogical(params, field_name // ' save mesh', .false.)
    
    if(save_mesh) then
      call LocalSaveMesh()
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
      if (associated(Vz)) then
        allocate(V(3,size(Vx % values)))
        call CopyPermMaybe(V(1,:), Vx%values, Vx%perm)
        call CopyPermMaybe(V(2,:), Vy%values, Vy%perm)
        call CopyPermMaybe(V(3,:), Vz%values, Vz%perm)
      else
        allocate(V(2,size(Vx%values)))
        call CopyPermMaybe(V(1,:), Vx%values, Vx%perm)
        call CopyPermMaybe(V(2,:), Vy%values, Vy%perm)
      end if

      call writer % write_data(trim(field_name), V(:,:))

      save_mesh = ListGetLogical(params, field_name // ' save mesh', found_save_mesh, defvalue = save_mesh)
      if(.not. found_save_mesh) call ListAddLogical(params, field_name // ' save mesh', .false.)

      if(save_mesh) then
        call LocalSaveMesh()
        save_mesh = .false.
      end if

    end do var_ind_do
  end block

  call writer % end_step()
contains

subroutine CopyPermMaybe(Y,X,perm)
  real(kind=dp), intent(out) :: Y(:)
  real(kind=dp), intent(in) :: X(:)
  integer, pointer :: perm(:)

  if(associated(perm)) then
    Y(:) = X(perm)
  else
    Y(:) = X(:)
  end if
end subroutine

  subroutine LocalSaveMesh()
    implicit none

        associate (meshelems => variable%primarymesh%elements)
          do state = 1,2
            lcon = 0
            do e_ind = 1,size(meshelems)
              nnodes = size(meshelems(e_ind)%nodeindexes)
              if (state == 2) then
                do n_ind = 1, nnodes
                  connectivity(lcon+n_ind) = meshelems(e_ind)%nodeindexes(n_ind)-1
                end do
                elem_types(e_ind) = Elmer2VTKElement(meshelems(e_ind)%type%elementcode, .false.)
                offsets(e_ind) = lcon
                num_elem_nodes(e_ind) = nnodes
              end if
              lcon = lcon + nnodes
            end do

            if (state == 2) offsets(size(meshelems)+1) = lcon

            if (state == 1) then
              allocate(connectivity(lcon))
              allocate(elem_types(size(meshelems)))
              allocate(num_elem_nodes(size(meshelems)))
              allocate(offsets(size(meshelems)+1))
            end if
          end do
        end associate
        call writer % write_data('num_verts', num_elem_nodes)
        call writer % write_data('cell_types', elem_types)
        call writer % write_data('connectivity', connectivity)
        call writer % write_data('points_x', variable % primarymesh % nodes % x)
        call writer % write_data('points_y', variable % primarymesh % nodes % y)
        call writer % write_data('points_z', variable % primarymesh % nodes % z)
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
