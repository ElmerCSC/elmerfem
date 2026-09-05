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
! *  Authors: fabien Gillet-Chaulet
! *  Email:   fabien.gillet-chaulet@univ-grenoble-alpes.fr
! *  Web:     http://elmerice.elmerfem.org
! *    
! *  Original Date: 2 June 2026
! *    
! ******************************************************************************/
!--------------------------------------------------------------------------------
!>  Interface to the C library proj : 
!>   Using the PROJ6 API
!>   See src API description for reference
!>       https://proj.org/en/stable/development/reference/index.html
!>   Not everything has been interfaced - code more for other feratures
!--------------------------------------------------------------------------------

MODULE proj6_interface
USE,INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE

!> ---------------------------------------------------------------------------
!> DERIVED TYPES: see proj C API TYPES

!> pj type
TYPE,BIND(C) :: pj_t
  PRIVATE
  TYPE(c_ptr) :: ptr = C_NULL_PTR
END TYPE pj_t

!> pj_context type
TYPE,BIND(C) :: pj_context_t
  PRIVATE
  TYPE(c_ptr) :: ptr = C_NULL_PTR
END TYPE pj_context_t

!> PJ_COORD union data type, for generic 4D coordinate handling.
TYPE,BIND(C) :: pj_coord_t
  REAL(kind=c_double) :: x=HUGE(1.0_c_double)
  REAL(kind=c_double) :: y=HUGE(1.0_c_double)
  REAL(kind=c_double) :: z=0.0_c_double
  REAL(kind=c_double) :: t=HUGE(1.0_c_double)
END TYPE pj_coord_t

!> pj_area type
TYPE,BIND(C) :: pj_area_t
  PRIVATE
  TYPE(c_ptr) :: ptr = C_NULL_PTR
END TYPE pj_area_t

!> pj_info type
!> accessed with proj_info fucntion
TYPE,BIND(C) :: pj_info_t
 INTEGER(kind=c_int) :: major !< Major release number
 INTEGER(kind=c_int) :: minor !< Minor release number
 INTEGER(kind=c_int) :: patch !< Patch level
 TYPE(c_ptr) :: release !< Release info. Version + date
 TYPE(c_ptr) :: version !< Full version number
 TYPE(c_ptr) :: searchpath !< Paths where init and grid files are looked for
 TYPE(c_ptr) :: paths !< C char**, use c_ptr_ptr for decoding in fortran
 INTEGER(kind=c_size_t) :: path_count
END TYPE pj_info_t

!> pj_proj_info type
TYPE,BIND(C) :: pj_proj_info_t
  TYPE(c_ptr) :: id                  !< Name of the projection in question
  TYPE(c_ptr) :: description         !< Description of the projection
  TYPE(c_ptr) :: definition          !< Projection definition
  INTEGER(kind=c_int) :: has_inverse !< 1 if an inverse mapping exists, 0 otherwise
  REAL(kind=c_double) :: accuracy    !< Expected accuracy of the transformation. -1 if unknown
END TYPE pj_proj_info_t

!> ---------------------------------------------------------------------------
!> PARAMETERS: define null pointers
TYPE(pj_t),PARAMETER :: pj_null=pj_t(C_NULL_PTR)
TYPE(pj_context_t),PARAMETER :: pj_default_ctx=pj_context_t(C_NULL_PTR)
TYPE(pj_area_t),PARAMETER :: pj_area_null=pj_area_t(C_NULL_PTR)

!> ---------------------------------------------------------------------------
!> ENUMs

!> PJ_DIRECTION
!> /* Apply transformation to observation - in forward or inverse direction */
ENUM, BIND(c) 
  ENUMERATOR :: &
   PJ_FWD = 1, &
   PJ_IDENT= 0, &
   PJ_INV = -1
END ENUM

!> ---------------------------------------------------------------------------
!> INTERFACES: see proj C API functions

!> ############################
!> Transformation setup:
INTERFACE
  FUNCTION proj_create(ctx, definition) BIND(C,name='proj_create')
  IMPORT
  TYPE(pj_context_t),VALUE :: ctx
  CHARACTER(kind=c_char) :: definition(*) 
  TYPE(pj_t) :: proj_create
  END FUNCTION proj_create
END INTERFACE

INTERFACE
  FUNCTION proj_create_crs_to_crs(ctx, source_crs, target_crs, area) &
   BIND(C,name='proj_create_crs_to_crs')
  IMPORT
  TYPE(pj_context_t),VALUE :: ctx
  CHARACTER(kind=c_char) :: source_crs(*) 
  CHARACTER(kind=c_char) :: target_crs(*) 
  TYPE(pj_area_t),VALUE :: area
  TYPE(pj_t) :: proj_create_crs_to_crs
  END FUNCTION proj_create_crs_to_crs
END INTERFACE

INTERFACE
  FUNCTION proj_destroy(pj) BIND(C,name='proj_destroy')
  IMPORT
  TYPE(pj_t),VALUE :: pj
  TYPE(pj_t) :: proj_destroy
  END FUNCTION proj_destroy
END INTERFACE

!> ############################
!>  Coordinate transformation:
INTERFACE
  FUNCTION proj_trans(p, direction, coord) BIND(C,name='proj_trans')
  IMPORT
  TYPE(pj_t),VALUE :: p
  INTEGER(kind=kind(PJ_FWD)),VALUE :: direction 
  TYPE(pj_coord_t),VALUE :: coord
  TYPE(pj_coord_t) :: proj_trans
  END FUNCTION proj_trans
END INTERFACE

!> ############################
!>  Info functions:
INTERFACE
  FUNCTION proj_info() BIND(C,name='proj_info')
  IMPORT
  TYPE(pj_info_t) :: proj_info
  END FUNCTION proj_info
END INTERFACE

INTERFACE
  FUNCTION proj_pj_info(p) BIND(C,name='proj_pj_info')
  IMPORT
  TYPE(pj_t),VALUE :: p
  TYPE(pj_proj_info_t) :: proj_pj_info
  END FUNCTION proj_pj_info
END INTERFACE

!> ############################
!>  Error reporting:
INTERFACE
  FUNCTION proj_errno(p) BIND(C,name='proj_errno')
  IMPORT
  TYPE(pj_t),VALUE :: p
  INTEGER(kind=c_int) :: proj_errno
  END FUNCTION proj_errno
END INTERFACE

INTERFACE
  FUNCTION proj_context_errno(ctx) BIND(C,name='proj_context_errno')
  IMPORT
  TYPE(pj_context_t),VALUE :: ctx
  INTEGER(kind=c_int) :: proj_context_errno
  END FUNCTION proj_context_errno
END INTERFACE

INTERFACE
  FUNCTION proj_context_errno_string(ctx, err) BIND(C,name='proj_context_errno_string')
  IMPORT
  TYPE(pj_context_t),VALUE :: ctx
  INTEGER(kind=c_int),VALUE :: err
  TYPE(c_ptr) :: proj_context_errno_string
  END FUNCTION proj_context_errno_string
END INTERFACE

!> Interface to c strlen to get len of a string c_ptr
INTERFACE
  FUNCTION c_strlen(str_ptr) bind ( C, name = "strlen" ) result(len)
  IMPORT
  type(c_ptr), value              :: str_ptr
  integer(kind=c_size_t)          :: len
  END FUNCTION c_strlen
END INTERFACE

!>--------------------------------------------------------------
!> Generic interface to test if a pj pointer is associated
INTERFACE proj_associated
  MODULE PROCEDURE proj_associated_pj, proj_associated_context, &
   proj_associated_area
END INTERFACE proj_associated

!>--------------------------------------------------------------
!> FUNCTIONS AND MODULES        
!>--------------------------------------------------------------
CONTAINS

!> Convert a string c_ptr to a fortran character        
FUNCTION cstrtof(c_pointer) result(f_string)
        use :: Types, ONLY: MAX_STRING_LEN
        implicit none
        type(c_ptr), intent(in)                 :: c_pointer
        character(len=:), pointer               :: f_ptr
        character(len=MAX_STRING_LEN)           :: f_string
        integer(c_size_t)                       :: l_str
    
        f_string=''
        IF (C_ASSOCIATED(c_pointer)) THEN
           call c_f_pointer(c_pointer, f_ptr )
           l_str = c_strlen(c_pointer)

           l_str=MIN(l_str,MAX_STRING_LEN)
           f_string(1:l_str) = f_ptr(1:l_str)
        ENDIF
end function cstrtof

!> Generic interface to test if a pj pointer is associated
FUNCTION proj_associated_pj(object) RESULT(associated_)
   TYPE(pj_t),INTENT(in) :: object
   LOGICAL :: associated_
   associated_ = C_ASSOCIATED(object%ptr)
END FUNCTION proj_associated_pj

FUNCTION proj_associated_context(object) RESULT(associated_)
   TYPE(pj_context_t),INTENT(in) :: object
   LOGICAL :: associated_
   associated_ = C_ASSOCIATED(object%ptr)
END FUNCTION proj_associated_context

FUNCTION proj_associated_area(object) RESULT(associated_)
   TYPE(pj_area_t),INTENT(in) :: object
   LOGICAL :: associated_
   associated_ = C_ASSOCIATED(object%ptr)
END FUNCTION proj_associated_area

END MODULE proj6_interface
