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
! *  Authors: Mikko Byckling
! *  Email:   mikko.byckling@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *  Original Date: 17 Feb 2014
! *
! *****************************************************************************/

! Module for wrapping and replacing functionality in Load.c

MODULE LoadMod
    USE Types
    USE Messages
    USE, INTRINSIC :: ISO_C_BINDING
    USE huti_interfaces
    IMPLICIT NONE

#ifdef ARCH_32_BITS
#define CAddrInt c_int32_t
#else
#define CAddrInt c_int64_t
#endif

#include "huti_fdefs.h"
#include "../config.h"

    ! GENERAL ROUTINES

    ! TODO: getsolverhome and makedirectory can be implemented with Fortran
    INTERFACE
        SUBROUTINE getsolverhome(solverDir, len) BIND(C,name='getsolverhome')
            USE, INTRINSIC :: iso_c_binding
            CHARACTER(C_CHAR) :: solverDir(*)
            INTEGER(C_INT) :: len
        END SUBROUTINE getsolverhome
    END INTERFACE

    INTERFACE
        SUBROUTINE getexedir(exeDir, len) BIND(C,name='getexedir')
            USE, INTRINSIC :: iso_c_binding
            CHARACTER(C_CHAR) :: exeDir(*)
            INTEGER(C_INT) :: len
        END SUBROUTINE getexedir
    END INTERFACE

    INTERFACE
        SUBROUTINE makedirectory(name) BIND(C,name='makedirectory')
            USE, INTRINSIC :: iso_c_binding
            CHARACTER(C_CHAR) :: name(*)
        END SUBROUTINE makedirectory
    END INTERFACE

    INTERFACE
        SUBROUTINE ElmerSetNoMPI(set) BIND(C,name='elmersetnompi')
            USE, INTRINSIC :: iso_c_binding
            INTEGER(C_INT) :: set
        END SUBROUTINE ElmerSetNoMPI
    END INTERFACE

    ! MATC

    INTERFACE
        SUBROUTINE matc_get_array(name, values, nrows, ncols) &
                   BIND(C,name='matc_get_array')
            USE, INTRINSIC :: ISO_C_BINDING
            CHARACTER(C_CHAR) :: name(*)
            REAL(C_DOUBLE) :: values(*)
            INTEGER(C_INT) :: nrows, ncols
        END SUBROUTINE matc_get_array
    END INTERFACE

    INTERFACE
        SUBROUTINE matc_c(cmd,cmdlen,res,reslen) BIND(C,name='matc_c')
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cmdlen,reslen
            CHARACTER(C_CHAR) :: cmd(*), res(*)
        END SUBROUTINE matc_c
    END INTERFACE

    INTERFACE
        SUBROUTINE matc_c_cached(cmd,cmdlen,res,reslen) BIND(C,name='matc_c_cached')
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cmdlen,reslen
            CHARACTER(C_CHAR) :: cmd(*), res(*)
        END SUBROUTINE matc_c_cached
    END INTERFACE

    INTERFACE
        SUBROUTINE matc_c_set_params(name,namelen,values,n) BIND(C,name='matc_c_set_params')
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: namelen, n
            CHARACTER(C_CHAR) :: name(*)
            REAL(C_DOUBLE) :: values(*)
        END SUBROUTINE matc_c_set_params
    END INTERFACE

    ! CPUTime.c
    INTERFACE
        FUNCTION cputime() RESULT(dbl) BIND(C,name='cputime')
            USE, INTRINSIC :: ISO_C_BINDING
            REAL(C_DOUBLE) :: dbl
        END FUNCTION cputime
    END INTERFACE

    INTERFACE
        FUNCTION realtime() RESULT(dbl) BIND(C,name='realtime')
            USE, INTRINSIC :: ISO_C_BINDING
            REAL(C_DOUBLE) :: dbl
        END FUNCTION realtime
    END INTERFACE

    INTERFACE
        FUNCTION cpumemory() RESULT(dbl) BIND(C,name='cpumemory')
            USE, INTRINSIC :: ISO_C_BINDING
            REAL(C_DOUBLE) :: dbl
        END FUNCTION cpumemory
    END INTERFACE

    ! fft.c
    TYPE, BIND(C) :: C_FFTCMPLX
        REAL(C_DOUBLE) :: rval
        REAL(C_DOUBLE) :: ival
    END TYPE C_FFTCMPLX

    INTERFACE
        SUBROUTINE frfftb( N, F, T ) BIND(C,name='frfftb')
            USE, INTRINSIC :: ISO_C_BINDING
            IMPORT C_FFTCMPLX
            INTEGER(C_INT) :: N
            TYPE(C_FFTCMPLX) :: F(*)
            REAL(C_DOUBLE) :: T(*)
        END SUBROUTINE frfftb
    END INTERFACE

    INTERFACE
        SUBROUTINE fcfftb( N, F, T ) BIND(C,name='fcfftb')
            USE, INTRINSIC :: ISO_C_BINDING
            IMPORT C_FFTCMPLX
            INTEGER(C_INT) :: N
            TYPE(C_FFTCMPLX) :: F(*), T(*)
        END SUBROUTINE fcfftb
    END INTERFACE

#if 0
    ! Overloads for AddrFunc
    INTERFACE AddrFunc
        MODULE PROCEDURE AddrFuncSub, &
                         AddrFuncInt, AddrFuncLong, &
                         AddrFuncReal, AddrFuncDbl, &
                         AddrFuncCmp, AddrFuncDblCmp
    END INTERFACE AddrFunc
#endif

    ! Overloads for IterCall
    INTERFACE IterCall
        MODULE PROCEDURE IterCallFTNR, IterCallFTNC
    END INTERFACE IterCall

    CONTAINS

        FUNCTION matc(cmd, res, inlen) RESULT(reslen)
          INTEGER :: reslen
          CHARACTER(*) :: cmd, res
          INTEGER, OPTIONAL :: inlen

          INTEGER :: cmdlen

          reslen = LEN(res)
          IF(PRESENT(inlen)) THEN
            cmdlen = inlen
          ELSE
            cmdlen = LEN_TRIM(cmd)
          END IF
          CALL matc_c(cmd,cmdlen,res,reslen)
        END FUNCTION matc

        ! Set a MATC variable directly from a Fortran array, bypassing
        ! string formatting and parsing.  On first call allocates the variable;
        ! on subsequent calls with the same size it is just a memcpy in C.
        SUBROUTINE MatcSetParams(name, values, n)
          CHARACTER(*) :: name
          REAL(KIND=dp) :: values(*)
          INTEGER :: n
          INTEGER :: namelen
          namelen = LEN_TRIM(name)
          CALL matc_c_set_params(name, namelen, values, n)
        END SUBROUTINE MatcSetParams

        ! Like matc() but compiles the expression once and re-uses the
        ! parse tree on subsequent calls with the same expression string.
        ! Use for MATC material/BC expressions evaluated per element.
        FUNCTION MatcCached(cmd, res, inlen) RESULT(reslen)
          INTEGER :: reslen
          CHARACTER(*) :: cmd, res
          INTEGER, OPTIONAL :: inlen

          INTEGER :: cmdlen

          reslen = LEN(res)
          IF(PRESENT(inlen)) THEN
            cmdlen = inlen
          ELSE
            cmdlen = LEN_TRIM(cmd)
          END IF
          CALL matc_c_cached(cmd,cmdlen,res,reslen)
        END FUNCTION MatcCached

        SUBROUTINE systemc(cmd, Status)
            IMPLICIT NONE
            CHARACTER(LEN=*) :: cmd
            ! Nonzero if the command could not be run or exited nonzero.  The
            ! return used to be dropped here, so a child that died left no
            ! trace beyond a message and the caller carried on as if it had
            ! worked.
            INTEGER, OPTIONAL, INTENT(OUT) :: Status
            CHARACTER(LEN=40) :: buf
            INTEGER :: estat, cstat

            INTERFACE 
                ! int system(const char *command)
                FUNCTION system(command) RESULT(retval)
                    USE, INTRINSIC :: ISO_C_BINDING
                    CHARACTER(C_CHAR), INTENT(IN) :: command(*)
                    INTEGER(C_INT) :: retval
                END FUNCTION system
            END INTERFACE

#ifdef HAVE_EXECUTECOMMANDLINE
            estat = 0; cstat = 0
            CALL EXECUTE_COMMAND_LINE(cmd, .TRUE., EXITSTAT=estat, CMDSTAT=cstat)
#else
            ! Workaround for Fortran compilers which do not 
            ! support EXECUTE_COMMAND_LINE intrinsic function
            cstat = 0
            estat = system(TRIM(cmd) // C_NULL_CHAR)
            IF (estat == -1) THEN
                cstat = estat
                estat = 0
            END IF
#endif
            IF (estat /= 0) THEN
                WRITE (buf, '(A,I0)') 'Command exit status was ', estat
                CALL Error('systemc',TRIM(BUF))
            END IF
            IF (cstat /= 0) THEN
                CALL Error('systemc','Unable to execute system command')
            END IF

            IF (PRESENT(Status)) THEN
                IF (cstat /= 0) THEN
                    Status = cstat
                ELSE
                    Status = estat
                END IF
            END IF
        END SUBROUTINE systemc

        SUBROUTINE envir(name, value, len)
            IMPLICIT NONE
            CHARACTER(LEN=*) :: name
            CHARACTER(LEN=*) :: value
            INTEGER :: len
            INTEGER :: estat

            CALL GET_ENVIRONMENT_VARIABLE(name, value, len)
        END SUBROUTINE envir

#if 0
        ! FUNCTION ADDRESS (overloaded methods for different function return values)
        FUNCTION AddrFuncSub(fn) RESULT(addr)
            IMPLICIT NONE
            PROCEDURE() :: fn
            INTEGER(KIND=AddrInt) :: addr

            TYPE(C_FUNPTR) :: cptr
            cptr = C_FUNLOC(fn)
            addr = TRANSFER(cptr, addr)
        END FUNCTION AddrFuncSub

        FUNCTION AddrFuncInt(fn) RESULT(addr)
            IMPLICIT NONE
            PROCEDURE(INTEGER(KIND=selected_int_kind(9))) :: fn
            INTEGER(KIND=AddrInt) :: addr

            TYPE(C_FUNPTR) :: cptr
            cptr = C_FUNLOC(fn)
            addr = TRANSFER(cptr, addr)
        END FUNCTION AddrFuncInt

        FUNCTION AddrFuncLong(fn) RESULT(addr)
            IMPLICIT NONE
            PROCEDURE(INTEGER(KIND=selected_int_kind(18))) :: fn
            INTEGER(KIND=AddrInt) :: addr

            TYPE(C_FUNPTR) :: cptr
            cptr = C_FUNLOC(fn)
            addr = TRANSFER(cptr, addr)
        END FUNCTION AddrFuncLong

        FUNCTION AddrFuncReal(fn) RESULT(addr)
            IMPLICIT NONE
            PROCEDURE(REAL(KIND=SELECTED_REAL_KIND(6))) :: fn
            INTEGER(KIND=AddrInt) :: addr

            TYPE(C_FUNPTR) :: cptr
            cptr = C_FUNLOC(fn)
            addr = TRANSFER(cptr, addr)
        END FUNCTION AddrFuncReal

        FUNCTION AddrFuncDbl(fn) RESULT(addr)
            IMPLICIT NONE
            PROCEDURE(REAL(KIND=dp)) :: fn
            INTEGER(KIND=AddrInt) :: addr

            TYPE(C_FUNPTR) :: cptr
            cptr = C_FUNLOC(fn)
            addr = TRANSFER(cptr, addr)
        END FUNCTION AddrFuncDbl

        FUNCTION AddrFuncCmp(fn) RESULT(addr)
            IMPLICIT NONE
            PROCEDURE(COMPLEX(KIND=SELECTED_REAL_KIND(6))) :: fn
            INTEGER(KIND=AddrInt) :: addr

            TYPE(C_FUNPTR) :: cptr
            cptr = C_FUNLOC(fn)
            addr = TRANSFER(cptr, addr)
        END FUNCTION AddrFuncCmp

        FUNCTION AddrFuncDblCmp(fn) RESULT(addr)
            IMPLICIT NONE
            PROCEDURE(COMPLEX(KIND=dp)) :: fn
            INTEGER(KIND=AddrInt) :: addr

            TYPE(C_FUNPTR) :: cptr
            cptr = C_FUNLOC(fn)
            addr = TRANSFER(cptr, addr)
        END FUNCTION AddrFuncDblCmp
#endif

        ! DYNAMIC LOADING  (wrapper via module procedure for typecasting)
        FUNCTION loadfunction(quiet, abort_not_found, library, fname,mangle) RESULT(ptr)
            IMPLICIT NONE
            INTEGER :: quiet, abort_not_found, mangle
            CHARACTER :: library(*), fname(*)
            TYPE(C_FUNPTR) :: ptr

            INTERFACE
                FUNCTION loadfunction_c(quiet, abort_not_found, library, fname, mangle ) RESULT(cptr) &
                    BIND(C,name='loadfunction_c')
                    USE, INTRINSIC :: ISO_C_BINDING
                    INTEGER(C_INT) :: quiet
                    INTEGER(C_INT) :: abort_not_found, mangle
                    CHARACTER(C_CHAR) :: library(*), fname(*)
                    TYPE(C_FUNPTR) :: cptr
                END FUNCTION loadfunction_c
            END INTERFACE

            ptr = loadfunction_c(quiet, abort_not_found, library, fname, mangle)
        END FUNCTION loadfunction

        ! DYNAMIC FUNCTION CALLS (wrappers via module procedures)
        RECURSIVE FUNCTION execintfunction(fptr, model ) RESULT(intval)
            IMPLICIT NONE
            TYPE(C_FUNPTR) :: fptr
            TYPE(Model_t), POINTER :: model
            INTEGER :: intval

            INTERFACE
                FUNCTION ElmerIntFn(model) RESULT(intval)
                    IMPORT Model_t
                    TYPE(Model_t) :: model
                    INTEGER :: intval
                END FUNCTION ElmerIntFn
            END INTERFACE
            PROCEDURE(ElmerIntFn), POINTER :: pptr

            CALL C_F_PROCPOINTER(fptr, pptr)
            intval = pptr(model)
        END FUNCTION execintfunction

        RECURSIVE FUNCTION execconstrealfunction(fptr, model, x, y, z) RESULT(realval)
            IMPLICIT NONE
            TYPE(C_FUNPTR) :: fptr
            TYPE(Model_t), POINTER :: model
            REAL(KIND=dp) :: x, y, z
            REAL(KIND=dp) :: realval

            INTERFACE
                FUNCTION ElmerConstRealFn(model, x, y, z) RESULT(realval)
                    IMPORT Model_t, dp
                    TYPE(Model_t) :: model
                    REAL(KIND=dp) :: x, y, z
                    REAL(KIND=dp) :: realval
                END FUNCTION ElmerConstRealFn
            END INTERFACE
            PROCEDURE(ElmerConstRealFn), POINTER :: pptr

            CALL C_F_PROCPOINTER(fptr, pptr)
            realval = pptr(model, x, y, z)
        END FUNCTION execconstrealfunction

        RECURSIVE FUNCTION execrealfunction(fptr, model, node, val) RESULT(realval)
            IMPLICIT NONE
            TYPE(C_FUNPTR) :: fptr
            TYPE(Model_t), POINTER :: model
            INTEGER :: node
            REAL(KIND=dp) :: val(*)
            REAL(KIND=dp) :: realval

            INTERFACE
                FUNCTION ElmerRealFn(model, node, val ) RESULT(realval)
                    IMPORT Model_t, dp
                    TYPE(Model_t) :: model
                    INTEGER :: node
                    REAL(KIND=dp) :: val(*)
                    REAL(KIND=dp) :: realval
                END FUNCTION ElmerRealFn
            END INTERFACE
            PROCEDURE(ElmerRealFn), POINTER :: pptr

            CALL C_F_PROCPOINTER(fptr, pptr)
            realval = pptr(model, node, val)
        END FUNCTION execrealfunction

        RECURSIVE SUBROUTINE execrealarrayfunction(fptr, model, node, val, arr )
            IMPLICIT NONE
            TYPE(C_FUNPTR) :: fptr
            TYPE(Model_t), POINTER :: model
            INTEGER :: node
            REAL(KIND=dp) :: val(*)
            REAL(KIND=dp) :: arr(:,:)

            INTERFACE
                SUBROUTINE ElmerRealArrFn(model, node, val, arr)
                    IMPORT Model_t, dp
                    TYPE(Model_t) :: model
                    INTEGER :: node
                    REAL(KIND=dp) :: val(*)
                    REAL(KIND=dp) :: arr(:,:)
                END SUBROUTINE ElmerRealArrFn
            END INTERFACE
            PROCEDURE(ElmerRealArrFn), POINTER :: pptr

            CALL C_F_PROCPOINTER(fptr, pptr)
            CALL pptr(model, node, val, arr)
        END SUBROUTINE execrealarrayfunction

        RECURSIVE SUBROUTINE execrealvectorfunction(fptr, model, node, val, arr )
            IMPLICIT NONE
            TYPE(C_FUNPTR) :: fptr
            TYPE(Model_t), POINTER :: model
            INTEGER :: node
            REAL(KIND=dp) :: val(*), arr(:)

            INTERFACE
                SUBROUTINE ElmerRealArrFn(model, node, val, arr)
                    IMPORT Model_t, dp
                    TYPE(Model_t) :: model
                    INTEGER :: node
                    REAL(KIND=dp) :: val(*), arr(:)
                END SUBROUTINE ElmerRealArrFn
            END INTERFACE
            PROCEDURE(ElmerRealArrFn), POINTER :: pptr

            CALL C_F_PROCPOINTER(fptr, pptr)
            CALL pptr(model, node, val, arr)
        END SUBROUTINE execrealvectorfunction

        RECURSIVE SUBROUTINE execsolver(fptr, model, solver, dt, transient)
            IMPLICIT NONE
            TYPE(C_FUNPTR) :: fptr
            TYPE(Model_t) :: model
            TYPE(Solver_t) :: solver
            REAL(KIND=dp) :: dt
            LOGICAL :: transient

            INTERFACE
                SUBROUTINE ElmerSolverFn(model, solver, dt, transient)
                    IMPORT Solver_t, Model_t, dp
                    TYPE(Model_t) :: model
                    TYPE(Solver_t) :: solver
                    REAL(KIND=dp) :: dt
                    LOGICAL :: transient
                END SUBROUTINE ElmerSolverFn
            END INTERFACE
            PROCEDURE(ElmerSolverFn), POINTER :: pptr

            CALL C_F_PROCPOINTER(fptr, pptr)
            CALL pptr(model, solver, dt, transient)
        END SUBROUTINE execsolver


        SUBROUTINE execmortarprojector(fptr, mesh, slavemesh, mastermesh, bcind, projector )
            IMPLICIT NONE
            TYPE(C_FUNPTR) :: fptr
            TYPE(Mesh_t) :: mesh, slavemesh, mastermesh
            INTEGER :: bcind
            TYPE(Matrix_t) :: projector

            INTERFACE
                SUBROUTINE MortarProjectorFn(mesh, slavemesh, mastermesh, bcind, projector )
                    IMPORT Mesh_t, Matrix_t
                    TYPE(Mesh_t) :: mesh, slavemesh, mastermesh
                    INTEGER :: bcind
                    TYPE(Matrix_t) :: projector
                END SUBROUTINE MortarProjectorFn
            END INTERFACE
            PROCEDURE(MortarProjectorFn), POINTER :: pptr

            CALL C_F_PROCPOINTER(fptr, pptr)
            CALL pptr(mesh, slavemesh, mastermesh, bcind, projector )
          END SUBROUTINE execmortarprojector
          
          FUNCTION enhancementfactoruserfunction( fptr, model, element, nodes, n, nd, &
                                       Basis, dBasisdx, Viscosity,Velo, dVelodx,sinvsq,localip ) &
                                       RESULT(realval)
            IMPLICIT NONE
            TYPE(C_FUNPTR) :: fptr
            TYPE(Model_t) :: model
            TYPE(Element_t), POINTER :: element
            TYPE(Nodes_t) :: nodes
            INTEGER :: n,nd,localip
            REAL(KIND=dp) :: Basis(:),dBasisdx(:,:),Viscosity, &
                             Velo(:), dVelodx(:,:),sinvsq
            REAL(KIND=dp) :: realval

            INTERFACE
                FUNCTION ElmerEnhancemntFactorFn(model, element, nodes, n, nd, &
                               Basis, dBasisdx, Viscosity, Velo, dVelodx, sinvsq,localip) RESULT(realval)
                    IMPORT Model_t, Element_t, Nodes_t, dp
                    TYPE(Model_t) :: model
                    TYPE(Element_t), POINTER :: element
                    TYPE(Nodes_t) :: nodes
                    INTEGER :: n,nd,localip
                    REAL(KIND=dp) :: Basis(:),dBasisdx(:,:),Viscosity, &
                                     Velo(:), dVelodx(:,:), sinvsq
                    REAL(KIND=dp) :: realval
                  END FUNCTION ElmerEnhancemntFactorFn
            END INTERFACE
            PROCEDURE(ElmerEnhancemntFactorFn), POINTER :: pptr

            CALL C_F_PROCPOINTER(fptr, pptr)
            realval = pptr(model, element, nodes, n, nd, &
                           Basis, dBasisdx, Viscosity,Velo, dVelodx,sinvsq,localip)
        END FUNCTION enhancementfactoruserfunction  

        
        FUNCTION materialuserfunction( fptr, model, element, nodes, n, nd, &
                                       Basis, dBasisdx, Viscosity,Velo, dVelodx ) &
                                       RESULT(realval)
            IMPLICIT NONE
            TYPE(C_FUNPTR) :: fptr
            TYPE(Model_t) :: model
            TYPE(Element_t), POINTER :: element
            TYPE(Nodes_t) :: nodes
            INTEGER :: n,nd
            REAL(KIND=dp) :: Basis(:),dBasisdx(:,:),Viscosity, &
                             Velo(:), dVelodx(:,:)
            REAL(KIND=dp) :: realval

            INTERFACE
                FUNCTION ElmerMaterialFn(model, element, nodes, n, nd, &
                               Basis, dBasisdx, Viscosity, Velo, dVelodx) RESULT(realval)
                    IMPORT Model_t, Element_t, Nodes_t, dp
                    TYPE(Model_t) :: model
                    TYPE(Element_t), POINTER :: element
                    TYPE(Nodes_t) :: nodes
                    INTEGER :: n,nd
                    REAL(KIND=dp) :: Basis(:),dBasisdx(:,:),Viscosity, &
                                     Velo(:), dVelodx(:,:)
                    REAL(KIND=dp) :: realval
                END FUNCTION ElmerMaterialFn
            END INTERFACE
            PROCEDURE(ElmerMaterialFn), POINTER :: pptr

            CALL C_F_PROCPOINTER(fptr, pptr)
            realval = pptr(model, element, nodes, n, nd, &
                           Basis, dBasisdx, Viscosity,Velo, dVelodx)
        END FUNCTION materialuserfunction

        SUBROUTINE execsimulationproc(fptr, model)
            IMPLICIT NONE
            TYPE(C_FUNPTR) :: fptr
            TYPE(Model_t) :: model

            INTERFACE
                SUBROUTINE ElmerSimulationFn(model)
                    IMPORT Model_t
                    TYPE(Model_t)   :: model
                END SUBROUTINE ElmerSimulationFn
            END INTERFACE
            PROCEDURE(ElmerSimulationFn), POINTER :: pptr

            CALL C_F_PROCPOINTER(fptr, pptr)
            CALL pptr(model)
        END SUBROUTINE execsimulationproc

        RECURSIVE FUNCTION execlinsolveprocs(fptr, model, solver, mtr, b, x, n, DOFs, nrm) RESULT(intval)
            IMPLICIT NONE
            TYPE(C_FUNPTR) :: fptr
            TYPE(Model_t) :: model
            TYPE(Solver_t) :: solver
            TYPE(Matrix_t), POINTER :: mtr
            INTEGER :: n, DOFs
            REAL(KIND=dp) :: x(n), b(n), nrm
            INTEGER :: intval

            INTERFACE
                FUNCTION ElmerLinSolveFn(model, solver, mtr, b, x, n, DOFs, nrm) RESULT(intval)
                    IMPORT Solver_t, Model_t, Matrix_t, dp
                    TYPE(Model_t) :: model
                    TYPE(Solver_t) :: solver
                    TYPE(Matrix_t), POINTER :: mtr
                    INTEGER :: n, DOFs
                    REAL(KIND=dp) :: x(n),b(n), nrm
                    INTEGER :: intval
                END FUNCTION ElmerLinSolveFn
            END INTERFACE
            PROCEDURE(ElmerLinSolveFn), POINTER :: pptr

            CALL C_F_PROCPOINTER(fptr, pptr)
            intval =  pptr(model, solver, mtr, b, x, n, DOFs, nrm)
        END FUNCTION execlinsolveprocs

        SUBROUTINE execlocalproc(fptr, model, solver, G, F, element, n, nd)
            IMPLICIT NONE
            TYPE(C_FUNPTR) :: fptr
            TYPE(Model_t)   :: model
            TYPE(Solver_t)  :: solver
            REAL(KIND=dp) :: G(:,:), F(:)
            TYPE(Element_t) :: element
            INTEGER :: n, nd

            INTERFACE
                SUBROUTINE ElmerLocalFn(model, solver, G, F, element, n, nd)
                    IMPORT Model_t, Solver_t, Element_t, dp
                    TYPE(Model_t)   :: model
                    TYPE(Solver_t)  :: solver
                    REAL(KIND=dp) :: G(:,:), F(:)
                    TYPE(Element_t) :: element
                    INTEGER :: n, nd
                END SUBROUTINE ElmerLocalFn
            END INTERFACE
            PROCEDURE(ElmerLocalFn), POINTER :: pptr

            CALL C_F_PROCPOINTER(fptr, pptr)
            CALL pptr(model, solver, G, F, element, n, nd)
        END SUBROUTINE execlocalproc

        SUBROUTINE execlocalassembly(fptr, model, solver, dt, transient, &
                                     M, D, S, F, element, nrow, ncol)
            IMPLICIT NONE
            TYPE(C_FUNPTR) :: fptr
            TYPE(Model_t)   :: model
            TYPE(Solver_t)  :: solver
            REAL(KIND=dp)   :: dt
            LOGICAL :: transient
            REAL(KIND=dp) :: M(:,:), D(:,:), S(:,:), F(:)
            TYPE(Element_t) :: element
            INTEGER :: nrow, ncol

            INTERFACE
                SUBROUTINE ElmerLocalAssemblyFn(model, solver, dt, transient, &
                                                M, D, S, F, element, Nrow, Ncol )
                    IMPORT Model_t, Solver_t, Element_t, dp
                    TYPE(Model_t)   :: model
                    TYPE(Solver_t)  :: solver
                    REAL(KIND=dp)   :: dt
                    LOGICAL :: transient
                    REAL(KIND=dp) :: M(:,:), D(:,:), S(:,:), F(:)
                    TYPE(Element_t) :: element
                    INTEGER :: nrow, ncol
                END SUBROUTINE ElmerLocalAssemblyFn
            END INTERFACE
            PROCEDURE(ElmerLocalAssemblyFn), POINTER :: pptr

            CALL C_F_PROCPOINTER(fptr, pptr)
            CALL pptr(model, solver, dt, transient, &
                      M, D, S, F, element, nrow, ncol)
        END SUBROUTINE execlocalassembly

        SUBROUTINE matvecsubrext(fptr, spmv, n, rows, cols, vals, u, v, reinit)
            IMPLICIT NONE

            TYPE(C_FUNPTR) :: fptr
            INTEGER(KIND=AddrInt) :: spmv
            INTEGER :: n
            INTEGER, POINTER  CONTIG :: rows(:), cols(:)
            REAL(KIND=dp), POINTER  CONTIG :: vals(:)
            REAL(KIND=dp), DIMENSION(*) :: u
            REAL(KIND=dp), DIMENSION(*) :: v
            INTEGER :: reinit
            INTEGER(KIND=CAddrInt) :: ifptr

            INTERFACE
                SUBROUTINE matvecsubrext_c(fptr, spmv, n, rows, cols, vals, u, v, reinit) &
                    BIND(C,name='matvecsubrext_c')
                    USE, INTRINSIC :: ISO_C_BINDING
                    INTEGER(CAddrInt) :: fptr
                    INTEGER(CAddrInt) :: spmv
                    INTEGER(C_INT) :: n
                    INTEGER(C_INT) :: rows(*), cols(*)
                    REAL(C_DOUBLE) :: vals(*)
                    REAL(C_DOUBLE) :: u(*),v(*)
                    INTEGER(C_INT) :: reinit
                END SUBROUTINE
            END INTERFACE

            ifptr = TRANSFER(fptr, ifptr)
            CALL matvecsubrext_c(ifptr, spmv, n, rows, cols, vals, u, v, reinit)
        END SUBROUTINE matvecsubrext

        RECURSIVE SUBROUTINE itercallR(fptr, x, b, ipar, dpar, work, &
                            mvptr, pcondptr, pcondrptr, dotptr, normptr, stopcptr )
            IMPLICIT NONE

            INTEGER(KIND=AddrInt) :: fptr
            REAL(KIND=dp), DIMENSION(:) CONTIG :: x,b
            INTEGER :: ipar(50)
            REAL(KIND=dp) :: dpar(50)
            REAL(KIND=dp) :: work(:,:)
            INTEGER(KIND=Addrint) :: mvptr, pcondptr, pcondrptr, &
                                     dotptr, normptr, stopcptr
            INTERFACE
                SUBROUTINE itercall_c(fptr, x, b, ipar, dpar, work, &
                            mvptr, pcondptr, pcondrptr, dotptr, normptr, stopcptr ) &
                    BIND(C,name='itercall_c')
                    USE, INTRINSIC :: ISO_C_BINDING
                    INTEGER(CAddrInt) :: fptr
                    REAL(C_DOUBLE) :: x(*), b(*)
                    INTEGER(C_INT) :: ipar(50)
                    REAL(C_DOUBLE) :: dpar(50)
                    REAL(C_DOUBLE) :: work(*)
                    INTEGER(CAddrInt) :: mvptr, pcondptr, pcondrptr, dotptr, &
                                     normptr, stopcptr
                END SUBROUTINE itercall_c
            END INTERFACE

            CALL itercall_c(fptr, x, b, ipar, dpar, work, &
                            mvptr, pcondptr, pcondrptr, dotptr, normptr, stopcptr)
        END SUBROUTINE itercallR

        RECURSIVE SUBROUTINE itercallC(fptr, x, b, ipar, dpar, work, &
                            mvptr, pcondptr, pcondrptr, dotptr, normptr, stopcptr )
            IMPLICIT NONE

            INTEGER(KIND=AddrInt) :: fptr
            COMPLEX(KIND=dp), DIMENSION(:) CONTIG :: x,b
            INTEGER :: ipar(50)
            REAL(KIND=dp) :: dpar(50)
            REAL(KIND=dp) :: work(:,:)
            INTEGER(KIND=Addrint) :: mvptr, pcondptr, pcondrptr, &
                                     dotptr, normptr, stopcptr
            INTERFACE
                SUBROUTINE itercall_c(fptr, x, b, ipar, dpar, work, &
                            mvptr, pcondptr, pcondrptr, dotptr, normptr, stopcptr ) &
                    BIND(C,name='itercall_c')
                    USE, INTRINSIC :: ISO_C_BINDING
                    INTEGER(CAddrInt) :: fptr
                    COMPLEX(C_DOUBLE_COMPLEX) :: x(*), b(*)
                    INTEGER(C_INT) :: ipar(50)
                    REAL(C_DOUBLE) :: dpar(50)
                    REAL(C_DOUBLE) :: work(*)
                    INTEGER(CAddrInt) :: mvptr, pcondptr, pcondrptr, dotptr, &
                                     normptr, stopcptr
                END SUBROUTINE itercall_c
            END INTERFACE

            CALL itercall_c(fptr, x, b, ipar, dpar, work, &
                            mvptr, pcondptr, pcondrptr, dotptr, normptr, stopcptr)
        END SUBROUTINE itercallC

        RECURSIVE SUBROUTINE itercallFTNR(fptr, x, b, ipar, dpar, work, &
                            mvptr, pcondptr, pcondrptr, dotptr, normptr, stopcptr )
            IMPLICIT NONE

            INTEGER(KIND=AddrInt) :: fptr
            REAL(KIND=dp), DIMENSION(:) CONTIG :: x,b
            INTEGER :: ipar(HUTI_IPAR_DFLTSIZE)
            REAL(KIND=dp) :: dpar(HUTI_DPAR_DFLTSIZE)
            REAL(KIND=dp) :: work(:,:)
            INTEGER(KIND=Addrint) :: mvptr, pcondptr, pcondrptr, &
                                     dotptr, normptr, stopcptr


            TYPE(C_FUNPTR) :: cfptr
            PROCEDURE(mv_iface_d), POINTER :: mvfun
            PROCEDURE(pc_iface_d), POINTER :: pcondfun, pcondrfun
            PROCEDURE(dotp_iface_d), POINTER :: dotfun
            PROCEDURE(norm_iface_d), POINTER :: normfun
            PROCEDURE(stopc_iface_d), POINTER :: stopcfun
            PROCEDURE(huti_itercall_d), POINTER :: iterfun

            ! Initialize pointers
            mvfun => NULL()
            pcondfun => NULL()
            pcondrfun => NULL()
            dotfun => NULL()
            normfun => NULL()
            iterfun => NULL()

            ! Transfer address integers back to Fortran function pointers

            ! Matrix-vector operator
            cfptr = TRANSFER(mvptr, cfptr)
            CALL C_F_PROCPOINTER(cfptr, mvfun)
            ! Preconditioner operators
            cfptr = TRANSFER(pcondptr, cfptr)
            IF (C_ASSOCIATED(cfptr)) CALL C_F_PROCPOINTER(cfptr, pcondfun)
            cfptr = TRANSFER(pcondrptr, cfptr)
            IF (C_ASSOCIATED(cfptr)) CALL C_F_PROCPOINTER(cfptr, pcondrfun)
            ! Dot product operator
            cfptr = TRANSFER(dotptr, cfptr)
            IF (C_ASSOCIATED(cfptr)) CALL C_F_PROCPOINTER(cfptr, dotfun)
            ! Norm operator
            cfptr = TRANSFER(normptr, cfptr)
            IF (C_ASSOCIATED(cfptr)) CALL C_F_PROCPOINTER(cfptr, normfun)
            ! Stopping criterion operator
            cfptr = TRANSFER(stopcptr, cfptr)
            IF (C_ASSOCIATED(cfptr)) CALL C_F_PROCPOINTER(cfptr, stopcfun)

            ! Finally, do the itercall
            cfptr = TRANSFER(fptr, cfptr)
            CALL C_F_PROCPOINTER(cfptr, iterfun)
            CALL iterfun(x, b, ipar, dpar, work, &
                         mvfun, pcondfun, pcondrfun, dotfun, normfun, stopcfun)
        END SUBROUTINE itercallFTNR

        RECURSIVE SUBROUTINE itercallFTNC(fptr, x, b, ipar, dpar, work, &
                            mvptr, pcondptr, pcondrptr, dotptr, normptr, stopcptr )
            IMPLICIT NONE

            INTEGER(KIND=AddrInt) :: fptr
            COMPLEX(KIND=dp), DIMENSION(:) CONTIG :: x,b
            INTEGER :: ipar(HUTI_IPAR_DFLTSIZE)
            REAL(KIND=dp) :: dpar(HUTI_DPAR_DFLTSIZE)
            COMPLEX(KIND=dp) :: work(:,:)
            INTEGER(KIND=Addrint) :: mvptr, pcondptr, pcondrptr, &
                                     dotptr, normptr, stopcptr


            TYPE(C_FUNPTR) :: cfptr
            PROCEDURE(mv_iface_z), POINTER :: mvfun
            PROCEDURE(pc_iface_z), POINTER :: pcondfun, pcondrfun
            PROCEDURE(dotp_iface_z), POINTER :: dotfun
            PROCEDURE(norm_iface_z), POINTER :: normfun
            PROCEDURE(stopc_iface_z), POINTER :: stopcfun
            PROCEDURE(huti_itercall_z), POINTER :: iterfun

            ! Initialize pointers
            mvfun => NULL()
            pcondfun => NULL()
            pcondrfun => NULL()
            dotfun => NULL()
            normfun => NULL()
            iterfun => NULL()

            ! Transfer address integers back to Fortran function pointers

            ! Matrix-vector operator
            cfptr = TRANSFER(mvptr, cfptr)
            CALL C_F_PROCPOINTER(cfptr, mvfun)
            ! Preconditioner operators
            cfptr = TRANSFER(pcondptr, cfptr)
            IF (C_ASSOCIATED(cfptr)) CALL C_F_PROCPOINTER(cfptr, pcondfun)
            cfptr = TRANSFER(pcondrptr, cfptr)
            IF (C_ASSOCIATED(cfptr)) CALL C_F_PROCPOINTER(cfptr, pcondrfun)
            ! Dot product operator
            cfptr = TRANSFER(dotptr, cfptr)
            IF (C_ASSOCIATED(cfptr)) CALL C_F_PROCPOINTER(cfptr, dotfun)
            ! Norm operator
            cfptr = TRANSFER(normptr, cfptr)
            IF (C_ASSOCIATED(cfptr)) CALL C_F_PROCPOINTER(cfptr, normfun)
            ! Stopping criterion operator
            cfptr = TRANSFER(stopcptr, cfptr)
            IF (C_ASSOCIATED(cfptr)) CALL C_F_PROCPOINTER(cfptr, stopcfun)
            
            ! Finally, do the itercall
            cfptr = TRANSFER(fptr, cfptr)
            CALL C_F_PROCPOINTER(cfptr, iterfun)
            CALL iterfun(x, b, ipar, dpar, work, &
                         mvfun, pcondfun, pcondrfun, dotfun, normfun, stopcfun)
        END SUBROUTINE itercallFTNC


        SUBROUTINE UMATusersubrtn( fptr, &            
            STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
            rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
            predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
            coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
            kstep, kinc)
          
          IMPLICIT NONE
          
          TYPE(C_FUNPTR) :: fptr
          REAL(KIND=dp), INTENT(INOUT) :: STRESS(NTENS)
          REAL(KIND=dp), INTENT(INOUT) :: STATEV(NSTATEV)
          REAL(KIND=dp), INTENT(OUT) :: DDSDDE(NTENS,NTENS)
          REAL(KIND=dp), INTENT(INOUT) :: SSE, SPD, SCD
          REAL(KIND=dp), INTENT(OUT) :: rpl
          REAL(KIND=dp), INTENT(OUT) :: ddsddt(NTENS), drplde(NTENS), drpldt
          REAL(KIND=dp), INTENT(IN) :: STRAN(NTENS)
          REAL(KIND=dp), INTENT(IN) :: DSTRAN(NTENS)
          REAL(KIND=dp), INTENT(IN) :: TIME(2)
          REAL(KIND=dp), INTENT(IN) :: DTIME
          REAL(KIND=dp), INTENT(IN) :: TEMP
          REAL(KIND=dp), INTENT(IN) :: dtemp
          REAL(KIND=dp), INTENT(IN) :: predef(1), dpred(1)
          CHARACTER(len=80), INTENT(IN) :: CMNAME
          INTEGER, INTENT(IN) :: NDI
          INTEGER, INTENT(IN) :: NSHR
          INTEGER, INTENT(IN) :: NTENS 
          INTEGER, INTENT(IN) :: NSTATEV
          REAL(KIND=dp), INTENT(IN) :: PROPS(NPROPS)
          INTEGER, INTENT(IN) :: NPROPS
          REAL(KIND=dp), INTENT(IN) :: coords(3)
          REAL(KIND=dp), INTENT(IN) :: drot(3,3)
          REAL(KIND=dp), INTENT(INOUT) :: pnewdt
          REAL(KIND=dp), INTENT(IN) :: celent
          REAL(KIND=dp), INTENT(IN) :: DFRGRD0(3,3)
          REAL(KIND=dp), INTENT(IN) :: DFRGRD1(3,3)
          INTEGER, INTENT(IN) :: NOEL
          INTEGER, INTENT(IN) :: NPT
          INTEGER, INTENT(IN) :: layer, kspt, kstep, kinc

          INTERFACE
            SUBROUTINE UMATsubrtn( STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
                rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
                predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
                coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
                kstep, kinc)
              
              USE Types
              IMPLICIT NONE
              
              REAL(KIND=dp), INTENT(INOUT) :: STRESS(NTENS)
              REAL(KIND=dp), INTENT(INOUT) :: STATEV(NSTATEV)
              REAL(KIND=dp), INTENT(OUT) :: DDSDDE(NTENS,NTENS)
              REAL(KIND=dp), INTENT(INOUT) :: SSE, SPD, SCD
              REAL(KIND=dp), INTENT(OUT) :: rpl
              REAL(KIND=dp), INTENT(OUT) :: ddsddt(NTENS), drplde(NTENS), drpldt
              REAL(KIND=dp), INTENT(IN) :: STRAN(NTENS)
              REAL(KIND=dp), INTENT(IN) :: DSTRAN(NTENS)
              REAL(KIND=dp), INTENT(IN) :: TIME(2)
              REAL(KIND=dp), INTENT(IN) :: DTIME
              REAL(KIND=dp), INTENT(IN) :: TEMP
              REAL(KIND=dp), INTENT(IN) :: dtemp
              REAL(KIND=dp), INTENT(IN) :: predef(1), dpred(1)
              CHARACTER(len=80), INTENT(IN) :: CMNAME
              INTEGER, INTENT(IN) :: NDI
              INTEGER, INTENT(IN) :: NSHR
              INTEGER, INTENT(IN) :: NTENS 
              INTEGER, INTENT(IN) :: NSTATEV
              REAL(KIND=dp), INTENT(IN) :: PROPS(NPROPS)
              INTEGER, INTENT(IN) :: NPROPS
              REAL(KIND=dp), INTENT(IN) :: coords(3)
              REAL(KIND=dp), INTENT(IN) :: drot(3,3)
              REAL(KIND=dp), INTENT(INOUT) :: pnewdt
              REAL(KIND=dp), INTENT(IN) :: celent
              REAL(KIND=dp), INTENT(IN) :: DFRGRD0(3,3)
              REAL(KIND=dp), INTENT(IN) :: DFRGRD1(3,3)
              INTEGER, INTENT(IN) :: NOEL
              INTEGER, INTENT(IN) :: NPT
              INTEGER, INTENT(IN) :: layer, kspt, kstep, kinc
            END SUBROUTINE UMATsubrtn
          END INTERFACE

          
          PROCEDURE(UMATsubrtn), POINTER :: pptr

          CALL C_F_PROCPOINTER(fptr, pptr)
          
          CALL pptr( STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
              rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
              predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
              coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
              kstep, kinc )
        END SUBROUTINE UMATusersubrtn
    
        
END MODULE LoadMod
