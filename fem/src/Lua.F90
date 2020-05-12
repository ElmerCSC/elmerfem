!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juhani Kataja
! *  Email:   juhani.kataja@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/
#define lua_upvalueindex(i)	(LUA_GLOBALSINDEX-(i))
#define LUA_GLOBALSINDEX	(-10002)
!-------------------------------------------------------------------------------
module Lua ! {{{
!-------------------------------------------------------------------------------
use ISO_C_BINDING
implicit none
private

!-Type declarations-------------------------------------------------------------
type, public :: LuaState_t
  private
  REAL(KIND=c_double), POINTER, PUBLIC :: tx(:) => NULL() ! This table will hold values for tx array
  type(c_ptr) :: L = c_null_ptr
  logical, public :: initialized=.false.
end type
!-------------------------------------------------------------------------------

type(LuaState_t), PUBLIC :: LuaState
!$OMP THREADPRIVATE(LuaState)

public :: lua_init, lua_close, lua_addfun, luaL_checkinteger, luaL_checknumber, &
    lua_pushnumber, luafun, lua_runfile, lua_dostring, &
    luaL_checkstring, lua_eval_f, lua_popnumber, lua_getnumber, lua_tolstring, &
    check_error, lua_getusertable, lua_poptensor, lua_popstring, lua_exec_fun, &
    lua_popvector

!-Interfaces-{{{----------------------------------------------------------------
interface ! 
  type(c_ptr) function lua_touserdata(L, n) bind(C)
    import
    type(c_ptr), value :: L
    integer(kind=c_int), value :: n
  end function

  function luaopen_array(L) result(n) bind(C, name="create_tx_table")
    import
    type(c_ptr), value :: L
    integer(kind=c_int) :: n
  end function

  function lua_tolstring_c(L, n, len) result(s) bind(C, name="lua_tolstring")
    import
    type(c_ptr), value :: L
    integer(kind=c_int), value :: n
    integer(kind=c_int) :: len
    type(c_ptr) :: s
  end function

  subroutine lua_pushnumber(L, x) bind(C)
    import
    type(c_ptr), value :: L
    real(kind=c_double), value :: x
  end subroutine

  function luaL_checknumber(L, n) result(r) bind(C, name="luaL_checknumber")
    import
    type(c_ptr), value :: L
    integer(kind=c_int), value :: n
    real(kind=c_double) :: r
  end function

  function lua_tonumber(L, n) result(r) bind(C)
    import
    type(c_ptr), value :: L
    integer(kind=c_int), value :: n
    real(kind=c_double) :: r
  end function

  function luaL_checkstring_c(L, n, len) result(s) bind(C, name="luaL_checklstring")
    import
    type(c_ptr), value :: L
    integer(kind=c_int), value :: n
    integer(kind=c_int) :: len
    type(c_ptr) :: s
  end function

  function luaL_checkinteger(L, n) result(r) bind(C, name="luaL_checkinteger")
    import
    type(c_ptr), value :: L
    integer(kind=c_int), value :: n
    integer(kind=c_int) :: r
  end function

  subroutine printfunloc(fn) bind(C)
    import
    type(c_funptr), value :: fn
  end subroutine

  function lua_init_c() result(L) bind(C, name="lua_init")
    import 
    type(c_ptr) :: L
  end function

  subroutine lua_close_c(L) bind(C, name = "lua_close")
    import
    type(c_ptr), value :: L
  end subroutine

  subroutine lua_pushcclosure(L, fun, n) bind(C, name  ="lua_pushcclosure")
    import
    type(c_ptr), value :: L
    type(c_funptr), value :: fun
    integer(kind=c_int), value :: n
  end subroutine

  function lua_cpcall(L, fun, ud) result(res) bind(C)
    import
    type(c_ptr), value :: L, fun, ud
    integer(kind=c_int) :: res
  end function

  subroutine lua_setfield(L,g_index, s) bind(C)
    import
    type(c_ptr), value :: L
    character(kind=c_char) :: s(*)
    integer(kind=c_int), value :: g_index
  end subroutine

  subroutine lua_getfield(L, g_index, s) bind(C)
    import
    type(c_ptr), value :: L
    character(kind=c_char) :: s(*)
    integer(kind=c_int), value :: g_index
  end subroutine

  subroutine lua_runfile_c(L, fname) bind(C, name="lua_runfile")
    import
    type(c_ptr), value :: L
    character(kind=c_char) :: fname(*)
  end subroutine

  integer(kind=c_int) function luaL_loadstring(L, s) bind(C, name="luaL_loadstring")
    import
    type(c_ptr), value :: L
    character(kind=c_char) :: s(*)
  end function

  integer(kind=c_int) function lua_pcall(L, a, b, c) bind(C)
    import
    type(c_ptr), value :: L
    integer(kind=c_int), value :: a, b, c
  end function

  subroutine lua_pop(L, n) bind(C, name="lua_pop_c")
    import
    type(c_ptr), value :: L
    integer(kind=c_int), value :: n
  end subroutine

  subroutine luaL_error(L, s) bind(C, name="luaL_error")
    import
    type(c_ptr), value :: L
    character(kind=c_char) :: s(*)
  end subroutine

  subroutine get_userdataptr(L, cp_raw, cp_data, len) bind(C)
    import
    type(c_ptr), value :: L, cp_raw
    type(c_ptr) :: cp_data
    integer(kind=c_int) :: len
  end subroutine

end interface 

abstract interface 
function luafun(L) result(n)
  import
  type(c_ptr), value :: L
  integer(kind=c_int) :: n
end function
end interface 
!-}}}---------------------------------------------------------------------------

CONTAINS

FUNCTION lua_getusertable(L, name) result(t)
  type(LuaState_t) :: L
  character(kind=c_char) :: name(*)
  real(kind=c_double), pointer :: t(:)

  type(c_ptr) :: cp_raw, cp_data
  integer(kind=c_int) :: len
  t => NULL()

  call lua_getfield(L%L, LUA_GLOBALSINDEX, name)
  cp_raw = lua_touserdata(L%L,-1)
  call get_userdataptr(L%L, cp_raw, cp_data, len) ! this should typecase cp_raw to "NumArray*" and extract size and len out of it
  call c_f_pointer(cp_data, t, shape=[len])
end function

function lua_tolstring(L, n, slen) result(sp)
  type(c_ptr) :: L
  integer(kind=c_int) :: n
  character(kind=c_char), pointer :: sp

  character(kind=c_char, len=:), allocatable :: s
  type(c_ptr) :: c_s
  integer(kind=c_int) :: slen

  c_s = lua_tolstring_c(L, n, slen)
  call c_f_pointer(c_s, sp)
end function

function luaL_checkstring(L, n, slen) result(sp)
  type(c_ptr) :: L
  integer(kind=c_int) :: n
  character(kind=c_char), pointer :: sp

  character(kind=c_char, len=:), allocatable :: s
  type(c_ptr) :: c_s
  integer(kind=c_int) :: slen

  c_s = luaL_checkstring_c(L, n, slen)
  call c_f_pointer(c_s, sp)
end function

subroutine lua_runfile(L, fname)
  type(LuaState_t) :: L
  character(kind=c_char) :: fname(*)
  call lua_runfile_c(L%L, fname)
end subroutine

function lua_dostring(L, s, m) result(n)
  type(LuaState_t) :: L
  character(kind=c_char) :: s(*)
  integer(kind=c_int) :: n, m_, load_error, pcall_error
  integer(kind=c_int), optional :: m
  if(.not. present(m)) then
    m_ = -1
  else
    m_ = m
  end if
  load_error = luaL_loadstring(L%L, s)
  call check_error(L, load_error)
  if(load_error == 0) then
    pcall_error = lua_pcall(L%L, 0, m_, 0)
    call check_error(L, pcall_error)
  end if
  n = IOR(load_error, pcall_error)
end function

real(kind=c_double) function lua_getnumber(L, s)
  type(LuaState_t) :: L
  character(kind=c_char) :: s(*)
  call lua_getfield(L%L, LUA_GLOBALSINDEX, s)
  lua_getnumber = lua_popnumber(L)
end function

subroutine lua_poptensor(L, t)
  type(LuaState_t) :: L
  real(kind=c_double), intent(out) :: t(:,:)

  integer :: n1, n2, i, j
  n1 = size(t, 1)
  n2 = size(t, 2)
  do i = n1,1, -1
    do j = n2,1,-1
      t(i,j) = lua_popnumber(L)
    end do
  end do
end subroutine

subroutine lua_popvector(L, t)
  type(LuaState_t) :: L
  real(kind=c_double), intent(out) :: t(:)
  integer :: n, i
  n = size(t, 1)
  do i = n, 1, -1
    t(i) = lua_popnumber(L)
  end do
end subroutine

subroutine lua_eval_f(L, fname, X, y)
  type(LuaState_t) :: L
  character(kind=c_char) :: fname(*)
  real(kind=c_double), intent(in) :: X(:)
  real(kind=c_double), intent(inout) :: Y(:)

  integer :: nx, ny, i, lstat
  nx = size(X,1)
  ny = size(Y,1)
  CALL lua_getfield(L%L, LUA_GLOBALSINDEX, fname)
  do i = 1,nx
    CALL lua_pushnumber(L%L, X(i))
  end do
  lstat = lua_pcall(L%L, nx, ny, 0)
  call check_error(L, lstat)
  if (lua_pcall(L%L, nx, ny, 0) /= 0) then
    CALL luaL_error(L%L, "error running '"//fname(1:len(fname))//"': ")
  end if
  do i = ny,1,-1
    Y(i) = lua_tonumber(L%L, -1) 
    CALL lua_pop(L%L,1)
  end do
end subroutine

!> Execute fname in lua state L, do not collect results from stack but expect user to collect them.
subroutine lua_exec_fun(L, fname, nin, nout)
  type(LuaState_t) :: L
  character(kind=c_char), intent(in) :: fname(*)
  integer, intent(in) :: nin, nout
  integer :: lstat

  CALL lua_getfield(L%L, LUA_GLOBALSINDEX, fname)
  lstat = lua_pcall(L%L, nin, nout, 0)
  call check_error(L, lstat)
end subroutine

real(kind=c_double) function lua_popnumber(L)
  type(LuaState_t) :: L
  lua_popnumber = lua_tonumber(L%L, -1)
  call lua_pop(L%L,1)
end function

function lua_init() result(L)
  type(LuaState_t) :: L
  type(c_ptr) :: ptr
  L%L = lua_init_c()
  if(c_associated(L%L)) L % initialized = .true.
end function

subroutine lua_addfun(L, fun, fname)
  type(LuaState_t) :: L
  procedure(luafun), pointer :: fun
  type(c_funptr) :: c_fun
  character(kind=c_char):: fname

  c_fun = c_funloc(fun)
  call lua_pushcclosure(L % L, c_fun, 0)
  call lua_setfield(L % L, LUA_GLOBALSINDEX ,fname)
end subroutine

subroutine lua_close(L)
  type(LuaState_t) :: L
  call lua_close_c(L%L)
  L % initialized = .false.
end subroutine

subroutine check_error(L, lstat)
  type(LuaState_t) :: L
  integer(kind=c_int), intent(in) :: lstat
  character(kind=c_char, len=:), pointer :: s
  integer(kind=c_int) :: slen
  if (lstat /= 0) then
    s => lua_tolstring(L%L, -1, slen)
    print *, 'Caught LUA error:', s(1:slen)
    call lua_pop(L%L,1);
  end if
end subroutine

function lua_popstring(L, slen) result(s)
  type(LuaState_t) :: L
  character(kind=c_char, len=:), pointer :: s
  integer :: slen
  s => lua_tolstring(L%L, -1, slen)
  call lua_pop(L%L, 1)
end function

!-------------------------------------------------------------------------------
end module ! Lua }}}
!-------------------------------------------------------------------------------
