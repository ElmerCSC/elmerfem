! *****************************************************************
! *                                                               *
! * iso_varying_string.f90                                        *
! *                                                               *
! * Copyright (C) 2003 Rich Townsend <rhdt@star.ucl.ac.uk>        *
! *                                                               *
! * This program is free software; you can redistribute it and/or *
! * modify it under the terms of the GNU Lesser General Public    *
! * License as published by the Free Software Foundation; either  *
! * version 2.1 of the License, or (at your option) any later     *
! * version.                                                      *
! *                                                               *
! * This program is distributed in the hope that it will be       *
! * useful, but WITHOUT ANY WARRANTY; without even the implied    *
! * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR       *
! * PURPOSE.  See the GNU Lesser General Public License for more  *
! * details.                                                      *
! *                                                               *
! * You should have received a copy of the GNU Lesser General     *
! * Public License along with this program; if not, write to the  *
! * Free Software Foundation, Inc., 59 Temple Place, Suite 330,   *
! * Boston, MA  02111-1307  USA                                   *
! *                                                               *
! *****************************************************************
!
! Developer : Rich Townsend <rhdt@star.ucl.ac.uk>
! Synopsis  : Definition of iso_varying_string module, conformant to 
!             the API specified in ISO/IEC 1539-2:2000 (varying-length 
!             strings for Fortran 95).
! Notes     : This implementation of iso_varying_string is designed to avoid
!             the possibility of memory leaks. To achieve this, it takes
!             advantage of language extensions specified in ISO/IEC
!             TR 15581 (enhanced data type facilities). Many vendors
!             support these extensions, and they will form a core part
!             of Fortran 2000.
! Version   : 1.2
! Thanks    : Lawrie Schonfelder's iso_varying_string module provided me
!             with much insight on how to go about writing this module,
!             for which I am very grateful. Furthermore, Lawrie helped
!             point out some subtle bugs in the module.

module iso_varying_string

! No implicit typing

  implicit none

! Parameter definitions

  integer, parameter :: GET_BUFFER_LEN = 256

! Type definitions

  type varying_string
     private
     character(LEN=1), dimension(:), allocatable :: chars
  end type varying_string

! Interface blocks

  interface assignment(=)
     module procedure op_assign_CH_VS
     module procedure op_assign_VS_CH
  end interface assignment(=)

  interface operator(//)
     module procedure op_concat_VS_VS
     module procedure op_concat_CH_VS
     module procedure op_concat_VS_CH
  end interface operator(//)

  interface operator(==)
     module procedure op_eq_VS_VS
     module procedure op_eq_CH_VS
     module procedure op_eq_VS_CH
  end interface operator(==)

  interface operator(/=)
     module procedure op_ne_VS_VS
     module procedure op_ne_CH_VS
     module procedure op_ne_VS_CH
  end interface operator (/=)
  
  interface operator(<)
     module procedure op_lt_VS_VS
     module procedure op_lt_CH_VS
     module procedure op_lt_VS_CH
  end interface operator (<)
  
  interface operator(<=)
     module procedure op_le_VS_VS
     module procedure op_le_CH_VS
     module procedure op_le_VS_CH
  end interface operator (<=)
  
  interface operator(>=)
     module procedure op_ge_VS_VS
     module procedure op_ge_CH_VS
     module procedure op_ge_VS_CH
  end interface operator (>=)

  interface operator(>)
     module procedure op_gt_VS_VS
     module procedure op_gt_CH_VS
     module procedure op_gt_VS_CH
  end interface operator (>)
  
  interface adjustl
     module procedure adjustl_
  end interface adjustl

  interface adjustr
     module procedure adjustr_
  end interface adjustr

  interface char
     module procedure char_auto
     module procedure char_fixed
  end interface char

  interface iachar
     module procedure iachar_
  end interface iachar

  interface ichar
     module procedure ichar_
  end interface ichar

  interface index
     module procedure index_VS_VS
     module procedure index_CH_VS
     module procedure index_VS_CH
  end interface index

  interface len
     module procedure len_
  end interface len

  interface len_trim
     module procedure len_trim_
  end interface len_trim

  interface lge
     module procedure lge_VS_VS
     module procedure lge_CH_VS
     module procedure lge_VS_CH
  end interface lge
  
  interface lgt
     module procedure lgt_VS_VS
     module procedure lgt_CH_VS
     module procedure lgt_VS_CH
  end interface lgt

  interface lle
     module procedure lle_VS_VS
     module procedure lle_CH_VS
     module procedure lle_VS_CH
  end interface lle

  interface llt
     module procedure llt_VS_VS
     module procedure llt_CH_VS
     module procedure llt_VS_CH
  end interface llt

  interface repeat
     module procedure repeat_
  end interface repeat

  interface scan
     module procedure scan_VS_VS
     module procedure scan_CH_VS
     module procedure scan_VS_CH
  end interface scan

  interface trim
     module procedure trim_
  end interface trim

  interface verify
     module procedure verify_VS_VS
     module procedure verify_CH_VS
     module procedure verify_VS_CH
  end interface verify

  interface var_str
     module procedure var_str_
  end interface var_str

  interface get
     module procedure get_
     module procedure get_unit
     module procedure get_set_VS
     module procedure get_set_CH
     module procedure get_unit_set_VS
     module procedure get_unit_set_CH
  end interface get

  interface put
     module procedure put_VS
     module procedure put_CH
     module procedure put_unit_VS
     module procedure put_unit_CH
  end interface put

  interface put_line
     module procedure put_line_VS
     module procedure put_line_CH
     module procedure put_line_unit_VS
     module procedure put_line_unit_CH
  end interface put_line

  interface extract
     module procedure extract_VS
     module procedure extract_CH
  end interface extract

  interface insert
     module procedure insert_VS_VS
     module procedure insert_CH_VS
     module procedure insert_VS_CH
     module procedure insert_CH_CH
  end interface insert

  interface remove
     module procedure remove_VS
     module procedure remove_CH
  end interface remove

  interface replace
     module procedure replace_VS_VS_auto
     module procedure replace_CH_VS_auto
     module procedure replace_VS_CH_auto
     module procedure replace_CH_CH_auto
     module procedure replace_VS_VS_fixed
     module procedure replace_CH_VS_fixed
     module procedure replace_VS_CH_fixed
     module procedure replace_CH_CH_fixed
     module procedure replace_VS_VS_VS_target
     module procedure replace_CH_VS_VS_target
     module procedure replace_VS_CH_VS_target
     module procedure replace_CH_CH_VS_target
     module procedure replace_VS_VS_CH_target
     module procedure replace_CH_VS_CH_target
     module procedure replace_VS_CH_CH_target
     module procedure replace_CH_CH_CH_target
  end interface

  interface split
     module procedure split_VS
     module procedure split_CH
  end interface split

! Access specifiers

  private

  public :: varying_string
  public :: assignment(=)
  public :: operator(//)
  public :: operator(==)
  public :: operator(/=)
  public :: operator(<)
  public :: operator(<=)
  public :: operator(>=)
  public :: operator(>)
  public :: adjustl
  public :: adjustr
  public :: char
  public :: iachar
  public :: ichar
  public :: index
  public :: len
  public :: len_trim
  public :: lge
  public :: lgt
  public :: lle
  public :: llt
  public :: repeat
  public :: scan
  public :: trim
  public :: verify
  public :: var_str
  public :: get
  public :: put
  public :: put_line
  public :: extract
  public :: insert
  public :: remove
  public :: replace
  public :: split

! Procedures

contains

!****

  elemental subroutine op_assign_CH_VS (var, exp)

    character(LEN=*), intent(out)    :: var
    type(varying_string), intent(in) :: exp

! Assign a varying string to a character string

    var = char(exp)

! Finish

    return

  end subroutine op_assign_CH_VS

!****

  elemental subroutine op_assign_VS_CH (var, exp)

    type(varying_string), intent(out) :: var
    character(LEN=*), intent(in)      :: exp

! Assign a character string to a varying string

    var = var_str(exp)

! Finish

    return

  end subroutine op_assign_VS_CH

!****

  elemental function op_concat_VS_VS (string_a, string_b) result (concat_string)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    type(varying_string)             :: concat_string

    integer                          :: len_string_a

! Concatenate two varying strings

    len_string_a = len(string_a)

    ALLOCATE(concat_string%chars(len_string_a+len(string_b)))

    concat_string%chars(:len_string_a) = string_a%chars
    concat_string%chars(len_string_a+1:) = string_b%chars

! Finish

    return

  end function op_concat_VS_VS

!****

  elemental function op_concat_CH_VS (string_a, string_b) result (concat_string)

    character(LEN=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    type(varying_string)             :: concat_string

! Concatenate a character string and a varying 
! string

    concat_string = op_concat_VS_VS(var_str(string_a), string_b)

! Finish

    return

  end function op_concat_CH_VS

!****

  elemental function op_concat_VS_CH (string_a, string_b) result (concat_string)

    type(varying_string), intent(in) :: string_a
    character(LEN=*), intent(in)     :: string_b
    type(varying_string)             :: concat_string

! Concatenate a varying string and a character
! string

    concat_string = op_concat_VS_VS(string_a, var_str(string_b))

! Finish

    return

  end function op_concat_VS_CH

!****

  elemental function op_eq_VS_VS (string_a, string_b) result (op_eq)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_eq

! Compare (==) two varying strings

    op_eq = char(string_a) == char(string_b)

! Finish

    return

  end function op_eq_VS_VS

!****

  elemental function op_eq_CH_VS (string_a, string_b) result (op_eq)

    character(LEN=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_eq

! Compare (==) a character string and a varying 
! string

    op_eq = string_a == char(string_b)

! Finish

    return

  end function op_eq_CH_VS

!****

  elemental function op_eq_VS_CH (string_a, string_b) result (op_eq)

    type(varying_string), intent(in) :: string_a
    character(LEN=*), intent(in)     :: string_b
    logical                          :: op_eq

! Compare (==) a varying string and a character
! string

    op_eq = char(string_a) == string_b

! Finish

    return

  end function op_eq_VS_CH

!****

  elemental function op_ne_VS_VS (string_a, string_b) result (op_ne)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_ne

! Compare (/=) two varying strings

    op_ne = char(string_a) /= char(string_b)

! Finish

    return

  end function op_ne_VS_VS

!****

  elemental function op_ne_CH_VS (string_a, string_b) result (op_ne)

    character(LEN=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_ne

! Compare (/=) a character string and a varying
! string

    op_ne = string_a /= char(string_b)

! Finish

    return

  end function op_ne_CH_VS

!****

  elemental function op_ne_VS_CH (string_a, string_b) result (op_ne)

    type(varying_string), intent(in) :: string_a
    character(LEN=*), intent(in)     :: string_b
    logical                          :: op_ne

! Compare (/=) a varying string and a character
! string

    op_ne = char(string_a) /= string_b

! Finish

    return

  end function op_ne_VS_CH

!****

  elemental function op_lt_VS_VS (string_a, string_b) result (op_lt)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_lt

! Compare (<) two varying strings

    op_lt = char(string_a) < char(string_b)

! Finish

    return

  end function op_lt_VS_VS

!****

  elemental function op_lt_CH_VS (string_a, string_b) result (op_lt)

    character(LEN=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_lt

! Compare (<) a character string and a varying
! string

    op_lt = string_a < char(string_b)

! Finish

    return

  end function op_lt_CH_VS

!****

  elemental function op_lt_VS_CH (string_a, string_b) result (op_lt)

    type(varying_string), intent(in) :: string_a
    character(LEN=*), intent(in)     :: string_b
    logical                          :: op_lt

! Compare (<) a varying string and a character 
! string

    op_lt = char(string_a) < string_b

! Finish

    return

  end function op_lt_VS_CH

!****

  elemental function op_le_VS_VS (string_a, string_b) result (op_le)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_le

! Compare (<=) two varying strings

    op_le = char(string_a) <= char(string_b)

! Finish

    return

  end function op_le_VS_VS

!****

  elemental function op_le_CH_VS (string_a, string_b) result (op_le)

    character(LEN=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_le

! Compare (<=) a character string and a varying 
! string

    op_le = string_a <= char(string_b)

! Finish

    return

  end function op_le_CH_VS

!****

  elemental function op_le_VS_CH (string_a, string_b) result (op_le)

    type(varying_string), intent(in) :: string_a
    character(LEN=*), intent(in)     :: string_b
    logical                          :: op_le

! Compare (<=) a varying string and a character 
! string

    op_le = char(string_a) <= string_b

! Finish

    return

  end function op_le_VS_CH

!****

  elemental function op_ge_VS_VS (string_a, string_b) result (op_ge)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_ge

! Compare (>=) two varying strings

    op_ge = char(string_a) >= char(string_b)

! Finish

    return

  end function op_ge_VS_VS

!****

  elemental function op_ge_CH_VS (string_a, string_b) result (op_ge)

    character(LEN=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_ge

! Compare (>=) a character string and a varying
! string

    op_ge = string_a >= char(string_b)

! Finish

    return

  end function op_ge_CH_VS

!****

  elemental function op_ge_VS_CH (string_a, string_b) result (op_ge)

    type(varying_string), intent(in) :: string_a
    character(LEN=*), intent(in)     :: string_b
    logical                          :: op_ge

! Compare (>=) a varying string and a character
! string

    op_ge = char(string_a) >= string_b

! Finish

    return

  end function op_ge_VS_CH

!****

  elemental function op_gt_VS_VS (string_a, string_b) result (op_gt)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_gt

! Compare (>) two varying strings

    op_gt = char(string_a) > char(string_b)

! Finish

    return

  end function op_gt_VS_VS

!****

  elemental function op_gt_CH_VS (string_a, string_b) result (op_gt)

    character(LEN=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_gt

! Compare (>) a character string and a varying
! string

    op_gt = string_a > char(string_b)

! Finish

    return

  end function op_gt_CH_VS

!****

  elemental function op_gt_VS_CH (string_a, string_b) result (op_gt)

    type(varying_string), intent(in) :: string_a
    character(LEN=*), intent(in)     :: string_b
    logical                          :: op_gt

! Compare (>) a varying string and a character
! string

    op_gt = char(string_a) > string_b

! Finish

    return

  end function op_gt_VS_CH

!****

  elemental function adjustl_ (string) result (adjustl_string)

    type(varying_string), intent(in) :: string
    type(varying_string)             :: adjustl_string

! Adjust the varying string to the left

    adjustl_string = ADJUSTL(CHAR(string))

! Finish

    return

  end function adjustl_

!****

  elemental function adjustr_ (string) result (adjustr_string)

    type(varying_string), intent(in) :: string
    type(varying_string)             :: adjustr_string

! Adjust the varying string to the right

    adjustr_string = ADJUSTR(CHAR(string))

! Finish

    return

  end function adjustr_

!****

  elemental function len_ (string) result (length)

    type(varying_string), intent(in) :: string
    integer                          :: length

! Get the length of a varying string

    if(ALLOCATED(string%chars)) then
       length = SIZE(string%chars)
    else
       length = 0
    endif

! Finish

    return

  end function len_

!****

  pure function char_auto (string) result (char_string)

    type(varying_string), intent(in) :: string
    character(LEN=len(string))       :: char_string

    integer                          :: i_char

! Convert a varying string into a character string
! (automatic length)
    
    forall(i_char = 1:len(string))
       char_string(i_char:i_char) = string%chars(i_char)
    end forall

! Finish

    return

  end function char_auto

!****

  pure function char_fixed (string, length) result (char_string)

    type(varying_string), intent(in) :: string
    integer, intent(in)              :: length
    character(LEN=length)            :: char_string

! Convert a varying string into a character string
! (fixed length)

    char_string = char(string)

! Finish

    return

  end function char_fixed

!****

  elemental function iachar_ (c) result (i)

    type(varying_string), intent(in) :: c
    integer                          :: i

! Get the position in the ISO 646 collating sequence
! of a varying string character

    i = IACHAR(char(c))

! Finish

    return

  end function iachar_

!****

  elemental function ichar_ (c) result (i)

    type(varying_string), intent(in) :: c
    integer                          :: i

! Get the position in the processor collating 
! sequence of a varying string character

    i = ICHAR(char(c))

! Finish

    return

  end function ichar_

!****

  elemental function index_VS_VS (string, substring, back) result (i_substring)

    type(varying_string), intent(in) :: string
    type(varying_string), intent(in) :: substring
    logical, intent(in), optional    :: back
    integer                          :: i_substring

! Get the index of a varying substring within a
! varying string

    i_substring = INDEX(char(string), char(substring), back)

! Finish

    return

  end function index_VS_VS

!****

  elemental function index_CH_VS (string, substring, back) result (i_substring)

    character(LEN=*), intent(in)     :: string
    type(varying_string), intent(in) :: substring
    logical, intent(in), optional    :: back
    integer                          :: i_substring

! Get the index of a varying substring within a
! character string

    i_substring = INDEX(string, char(substring), back)

! Finish

    return

  end function index_CH_VS

!****

  elemental function index_VS_CH (string, substring, back) result (i_substring)

    type(varying_string), intent(in) :: string
    character(LEN=*), intent(in)     :: substring
    logical, intent(in), optional    :: back
    integer                          :: i_substring

! Get the index of a character substring within a
! varying string

    i_substring = INDEX(char(string), substring, back)

! Finish

    return

  end function index_VS_CH

!****


  elemental function len_trim_ (string) result (length)

    type(varying_string), intent(in) :: string
    integer                          :: length

! Get the trimmed length of a varying string

    if(ALLOCATED(string%chars)) then
       length = LEN_TRIM(char(string))
    else
       length = 0
    endif

! Finish

    return

  end function len_trim_

!****

  elemental function lge_VS_VS (string_a, string_b) result (comp)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: comp

! Compare (LGE) two varying strings

    comp = LGE(char(string_a), char(string_b))

! Finish

    return

  end function lge_VS_VS

!****

  elemental function lge_CH_VS (string_a, string_b) result (comp)

    character(LEN=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: comp

! Compare (LGE) a character string and a varying
! string

    comp = LGE(string_a, char(string_b))

! Finish

    return

  end function lge_CH_VS

!****

  elemental function lge_VS_CH (string_a, string_b) result (comp)

    type(varying_string), intent(in) :: string_a
    character(LEN=*), intent(in)     :: string_b
    logical                          :: comp

! Compare (LGE) a varying string and a character
! string

    comp = LGE(char(string_a), string_b)

! Finish

    return

  end function lge_VS_CH

!****

  elemental function lgt_VS_VS (string_a, string_b) result (comp)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: comp

! Compare (LGT) two varying strings

    comp = LGT(char(string_a), char(string_b))

! Finish

    return

  end function lgt_VS_VS

!****

  elemental function lgt_CH_VS (string_a, string_b) result (comp)

    character(LEN=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: comp

! Compare (LGT) a character string and a varying
! string

    comp = LGT(string_a, char(string_b))

! Finish

    return

  end function lgt_CH_VS

!****

  elemental function lgt_VS_CH (string_a, string_b) result (comp)

    type(varying_string), intent(in) :: string_a
    character(LEN=*), intent(in)     :: string_b
    logical                          :: comp

! Compare (LGT) a varying string and a character
! string

    comp = LGT(char(string_a), string_b)

! Finish

    return

  end function lgt_VS_CH

!****

  elemental function lle_VS_VS (string_a, string_b) result (comp)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: comp

! Compare (LLE) two varying strings

    comp = LLE(char(string_a), char(string_b))

! Finish

    return

  end function lle_VS_VS

!****

  elemental function lle_CH_VS (string_a, string_b) result (comp)

    character(LEN=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: comp

! Compare (LLE) a character string and a varying
! string

    comp = LLE(string_a, char(string_b))

! Finish

    return

  end function lle_CH_VS

!****

  elemental function lle_VS_CH (string_a, string_b) result (comp)

    type(varying_string), intent(in) :: string_a
    character(LEN=*), intent(in)     :: string_b
    logical                          :: comp

! Compare (LLE) a varying string and a character
! string

    comp = LLE(char(string_a), string_b)

! Finish

    return

  end function lle_VS_CH

!****

  elemental function llt_VS_VS (string_a, string_b) result (comp)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: comp

! Compare (LLT) two varying strings

    comp = LLT(char(string_a), char(string_b))

! Finish

    return

  end function llt_VS_VS

!****

  elemental function llt_CH_VS (string_a, string_b) result (comp)

    character(LEN=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: comp

! Compare (LLT) a character string and a varying
! string

    comp = LLT(string_a, char(string_b))

! Finish

    return

  end function llt_CH_VS

!****

  elemental function llt_VS_CH (string_a, string_b) result (comp)

    type(varying_string), intent(in) :: string_a
    character(LEN=*), intent(in)     :: string_b
    logical                          :: comp

! Compare (LLT) a varying string and a character
! string

    comp = LLT(char(string_a), string_b)

! Finish

    return

  end function llt_VS_CH

!****

  elemental function repeat_ (string, ncopies) result (repeat_string)

    type(varying_string), intent(in) :: string
    integer, intent(in)              :: ncopies
    type(varying_string)             :: repeat_string

! Concatenate several copies of a varying string

    repeat_string = var_str(REPEAT(char(string), ncopies))

! Finish

    return

  end function repeat_

!****

  elemental function scan_VS_VS (string, set, back) result (i)

    type(varying_string), intent(in) :: string
    type(varying_string), intent(in) :: set
    logical, intent(in), optional    :: back
    integer                          :: i

! Scan a varying string for occurrences of 
! characters in a varying-string set

    i = SCAN(char(string), char(set), back)

! Finish

    return

  end function scan_VS_VS
    
!****

  elemental function scan_CH_VS (string, set, back) result (i)

    character(LEN=*), intent(in)     :: string
    type(varying_string), intent(in) :: set
    logical, intent(in), optional    :: back
    integer                          :: i

! Scan a character string for occurrences of 
! characters in a varying-string set

    i = SCAN(string, char(set), back)

! Finish

    return

  end function scan_CH_VS
    
!****

  elemental function scan_VS_CH (string, set, back) result (i)

    type(varying_string), intent(in) :: string
    character(LEN=*), intent(in)     :: set
    logical, intent(in), optional    :: back
    integer                          :: i

! Scan a varying string for occurrences of 
! characters in a character-string set

    i = SCAN(char(string), set, back)

! Finish

    return

  end function scan_VS_CH
    
!****

  elemental function trim_ (string) result (trim_string)

    type(varying_string), intent(in) :: string
    type(varying_string)             :: trim_string

! Remove trailing blanks from a varying string

    trim_string = TRIM(char(string))

! Finish

    return

  end function trim_

!****

  elemental function verify_VS_VS (string, set, back) result (i)

    type(varying_string), intent(in) :: string
    type(varying_string), intent(in) :: set
    logical, intent(in), optional    :: back
    integer                          :: i

! Verify a varying string for occurrences of
! characters in a varying-string set

    i = VERIFY(char(string), char(set), back)

! Finish

    return

  end function verify_VS_VS

!****

  elemental function verify_CH_VS (string, set, back) result (i)

    character(LEN=*), intent(in)     :: string
    type(varying_string), intent(in) :: set
    logical, intent(in), optional    :: back
    integer                          :: i

! Verify a character string for occurrences of 
! characters in a varying-string set

    i = VERIFY(string, char(set), back)

! Finish

    return

  end function verify_CH_VS

!****

  elemental function verify_VS_CH (string, set, back) result (i)

    type(varying_string), intent(in) :: string
    character(LEN=*), intent(in)     :: set
    logical, intent(in), optional    :: back
    integer                          :: i

! Verify a varying string for occurrences of 
! characters in a character-string set

    i = VERIFY(char(string), set, back)

! Finish

    return

  end function verify_VS_CH

!****

  elemental function var_str_ (char) result (string)

    character(LEN=*), intent(in) :: char
    type(varying_string)         :: string

    integer                      :: length
    integer                      :: i_char

! Convert a character string to a varying string

    length = LEN(char)

    ALLOCATE(string%chars(length))

    forall(i_char = 1:length)
       string%chars(i_char) = char(i_char:i_char)
    end forall

! Finish

    return

  end function var_str_

!****

  subroutine get_ (string, maxlen, iostat)

    type(varying_string), intent(out) :: string
    integer, intent(in), optional     :: maxlen
    integer, intent(out), optional    :: iostat

    integer                           :: n_chars_remain
    integer                           :: n_chars_read
    character(LEN=GET_BUFFER_LEN)     :: buffer

! Read from the default unit into a varying string

    string = ''

    if(PRESENT(maxlen)) then
       n_chars_remain = maxlen
    else
       n_chars_remain = HUGE(1)
    endif

    read_loop : do

       if(n_chars_remain <= 0) return

       n_chars_read = MIN(n_chars_remain, GET_BUFFER_LEN)

       if(PRESENT(iostat)) then
          read(*, FMT='(A)', ADVANCE='NO', IOSTAT=iostat, SIZE=n_chars_read) buffer(:n_chars_read)
          if(iostat < 0) exit read_loop
          if(iostat > 0) return
       else
          read(*, FMT='(A)', ADVANCE='NO', EOR=999, SIZE=n_chars_read) buffer(:n_chars_read)
       endif

       string = string//buffer(:n_chars_read)
       n_chars_remain = n_chars_remain - n_chars_read

    end do read_loop

999 continue

    string = string//buffer(:n_chars_read)

! Finish (end-of-record)

    return

  end subroutine get_

!****

  subroutine get_unit (unit, string, maxlen, iostat)

    integer, intent(in)               :: unit
    type(varying_string), intent(out) :: string
    integer, intent(in), optional     :: maxlen
    integer, intent(out), optional    :: iostat

    integer                           :: n_chars_remain
    integer                           :: n_chars_read
    character(LEN=GET_BUFFER_LEN)     :: buffer

! Read from the specified unit into a varying string

    string = ''

    if(PRESENT(maxlen)) then
       n_chars_remain = maxlen
    else
       n_chars_remain = HUGE(1)
    endif

    read_loop : do

       if(n_chars_remain <= 0) return

       n_chars_read = MIN(n_chars_remain, GET_BUFFER_LEN)

       if(PRESENT(iostat)) then
          read(unit, FMT='(A)', ADVANCE='NO', IOSTAT=iostat, SIZE=n_chars_read) buffer(:n_chars_read)
          if(iostat < 0) exit read_loop
          if(iostat > 0) return
       else
          read(unit, FMT='(A)', ADVANCE='NO', EOR=999, SIZE=n_chars_read) buffer(:n_chars_read)
       endif

       string = string//buffer(:n_chars_read)
       n_chars_remain = n_chars_remain - n_chars_read

    end do read_loop

999 continue

    string = string//buffer(:n_chars_read)

! Finish (end-of-record)

    return

  end subroutine get_unit

!****

  subroutine get_set_VS (string, set, separator, maxlen, iostat)

    type(varying_string), intent(out)           :: string
    type(varying_string), intent(in)            :: set
    type(varying_string), intent(out), optional :: separator
    integer, intent(in), optional               :: maxlen
    integer, intent(out), optional              :: iostat

! Read from the default unit into a varying string,
! with a custom varying-string separator

    call get(string, char(set), separator, maxlen, iostat)

! Finish

    return

  end subroutine get_set_VS

!****

  subroutine get_set_CH (string, set, separator, maxlen, iostat)

    type(varying_string), intent(out)           :: string
    character(LEN=*), intent(in)                :: set
    type(varying_string), intent(out), optional :: separator
    integer, intent(in), optional               :: maxlen
    integer, intent(out), optional              :: iostat

    integer                                     :: n_chars_remain
    character(LEN=1)                            :: buffer
    integer                                     :: i_set

! Read from the default unit into a varying string,
! with a custom character-string separator

    string = ''

    if(PRESENT(maxlen)) then
       n_chars_remain = maxlen
    else
       n_chars_remain = HUGE(1)
    endif

    if(PRESENT(separator)) separator = ''

    read_loop : do

       if(n_chars_remain <= 0) return

       if(PRESENT(iostat)) then
          read(*, FMT='(A1)', ADVANCE='NO', IOSTAT=iostat) buffer
          if(iostat /= 0) exit read_loop
       else
          read(*, FMT='(A1)', ADVANCE='NO', EOR=999) buffer
       endif

       i_set = SCAN(buffer, set)

       if(i_set == 1) then
          if(PRESENT(separator)) separator = buffer
          exit read_loop
       endif

       string = string//buffer
       n_chars_remain = n_chars_remain - 1

    end do read_loop

999 continue

! Finish

    return

  end subroutine get_set_CH

!****

  subroutine get_unit_set_VS (unit, string, set, separator, maxlen, iostat)

    integer, intent(in)                         :: unit
    type(varying_string), intent(out)           :: string
    type(varying_string), intent(in)            :: set
    type(varying_string), intent(out), optional :: separator
    integer, intent(in), optional               :: maxlen
    integer, intent(out), optional              :: iostat

! Read from the specified unit into a varying string,
! with a custom varying-string separator

    call get(unit, string, char(set), separator, maxlen, iostat)

! Finish

    return

  end subroutine get_unit_set_VS

!****

  subroutine get_unit_set_CH (unit, string, set, separator, maxlen, iostat)

    integer, intent(in)                         :: unit
    type(varying_string), intent(out)           :: string
    character(LEN=*), intent(in)                :: set
    type(varying_string), intent(out), optional :: separator
    integer, intent(in), optional               :: maxlen
    integer, intent(out), optional              :: iostat

    integer                                     :: n_chars_remain
    character(LEN=1)                            :: buffer
    integer                                     :: i_set

! Read from the default unit into a varying string,
! with a custom character-string separator

    string = ''

    if(PRESENT(maxlen)) then
       n_chars_remain = maxlen
    else
       n_chars_remain = HUGE(1)
    endif

    if(PRESENT(separator)) separator = ''

    read_loop : do

       if(n_chars_remain <= 0) return

       if(PRESENT(iostat)) then
          read(unit, FMT='(A1)', ADVANCE='NO', IOSTAT=iostat) buffer
          if(iostat /= 0) exit read_loop
       else
          read(unit, FMT='(A1)', ADVANCE='NO', EOR=999) buffer
       endif

       i_set = SCAN(buffer, set)

       if(i_set == 1) then
          if(PRESENT(separator)) separator = buffer
          exit read_loop
       endif

       string = string//buffer
       n_chars_remain = n_chars_remain - 1

    end do read_loop

999 continue

! Finish

    return

  end subroutine get_unit_set_CH

!****

  subroutine put_VS (string, iostat)

    type(varying_string), intent(in) :: string
    integer, intent(out), optional   :: iostat

! Append a varying string to the current record of 
! the default unit

    call put(char(string), iostat)

! Finish

  end subroutine put_VS
    
!****

  subroutine put_CH (string, iostat)

    character(LEN=*), intent(in)   :: string
    integer, intent(out), optional :: iostat

! Append a character string to the current record of 
! the default unit

    if(PRESENT(iostat)) then
       write(*, FMT='(A)', ADVANCE='NO', IOSTAT=iostat) string
    else
       write(*, FMT='(A)', ADVANCE='NO') string
    endif

! Finish

  end subroutine put_CH

!****

  subroutine put_unit_VS (unit, string, iostat)

    integer, intent(in)              :: unit
    type(varying_string), intent(in) :: string
    integer, intent(out), optional   :: iostat

! Append a varying string to the current record of 
! the specified unit

    call put(unit, char(string), iostat)

! Finish

    return

  end subroutine put_unit_VS

!****

  subroutine put_unit_CH (unit, string, iostat)

    integer, intent(in)            :: unit
    character(LEN=*), intent(in)   :: string
    integer, intent(out), optional :: iostat

! Append a character string to the current record of 
! the specified unit

    if(PRESENT(iostat)) then
       write(unit, FMT='(A)', ADVANCE='NO', IOSTAT=iostat) string
    else
       write(unit, FMT='(A)', ADVANCE='NO') string
    endif

! Finish

    return

  end subroutine put_unit_CH

!****

  subroutine put_line_VS (string, iostat)

    type(varying_string), intent(in) :: string
    integer, intent(out), optional   :: iostat

! Append a varying string to the current record of 
! the default unit, terminating the record

    call put_line(char(string), iostat)

! Finish

    return

  end subroutine put_line_VS

!****

  subroutine put_line_CH (string, iostat)

    character(LEN=*), intent(in)   :: string
    integer, intent(out), optional :: iostat

! Append a varying string to the current record of 
! the default unit, terminating the record

    if(PRESENT(iostat)) then
       write(*, FMT='(A,/)', ADVANCE='NO', IOSTAT=iostat) string
    else
       write(*, FMT='(A,/)', ADVANCE='NO') string
    endif

! Finish

    return

  end subroutine put_line_CH

!****

  subroutine put_line_unit_VS (unit, string, iostat)

    integer, intent(in)              :: unit
    type(varying_string), intent(in) :: string
    integer, intent(out), optional   :: iostat

! Append a varying string to the current record of 
! the specified unit, terminating the record

    call put_line(unit, char(string), iostat)

! Finish

    return

  end subroutine put_line_unit_VS

!****

  subroutine put_line_unit_CH (unit, string, iostat)

    integer, intent(in)            :: unit
    character(LEN=*), intent(in)   :: string
    integer, intent(out), optional :: iostat

! Append a varying string to the current record of 
! the specified unit, terminating the record

    if(PRESENT(iostat)) then
       write(unit, FMT='(A,/)', ADVANCE='NO', IOSTAT=iostat) string
    else
       write(unit, FMT='(A,/)', ADVANCE='NO') string
    endif

! Finish

    return

  end subroutine put_line_unit_CH

!****

  elemental function extract_VS (string, start, finish) result (ext_string)

    type(varying_string), intent(in) :: string
    integer, intent(in), optional    :: start
    integer, intent(in), optional    :: finish
    type(varying_string)             :: ext_string

! Extract a varying substring from a varying string

    ext_string = extract(char(string), start, finish)

! Finish

    return

  end function extract_VS

!****

  elemental function extract_CH (string, start, finish) result (ext_string)

    character(LEN=*), intent(in)  :: string
    integer, intent(in), optional :: start
    integer, intent(in), optional :: finish
    type(varying_string)          :: ext_string

    integer                       :: start_
    integer                       :: finish_

! Extract a varying substring from a character string

    if(PRESENT(start)) then
       start_ = MAX(1, start)
    else
       start_ = 1
    endif

    if(PRESENT(finish)) then
       finish_ = MIN(LEN(string), finish)
    else
       finish_ = LEN(string)
    endif

    ext_string = var_str(string(start_:finish_))

! Finish

    return

  end function extract_CH

!****

  elemental function insert_VS_VS (string, start, substring) result (ins_string)

    type(varying_string), intent(in) :: string
    integer, intent(in)              :: start
    type(varying_string), intent(in) :: substring
    type(varying_string)             :: ins_string

! Insert a varying substring into a varying string

    ins_string = insert(char(string), start, char(substring))

! Finish

    return

  end function insert_VS_VS

!****

  elemental function insert_CH_VS (string, start, substring) result (ins_string)

    character(LEN=*), intent(in)     :: string
    integer, intent(in)              :: start
    type(varying_string), intent(in) :: substring
    type(varying_string)             :: ins_string

! Insert a varying substring into a character string

    ins_string = insert(string, start, char(substring))

! Finish

    return

  end function insert_CH_VS

!****

  elemental function insert_VS_CH (string, start, substring) result (ins_string)

    type(varying_string), intent(in) :: string
    integer, intent(in)              :: start
    character(LEN=*), intent(in)     :: substring
    type(varying_string)             :: ins_string

! Insert a character substring into a varying string

    ins_string = insert(char(string), start, substring)

! Finish

    return

  end function insert_VS_CH

!****

  elemental function insert_CH_CH (string, start, substring) result (ins_string)

    character(LEN=*), intent(in) :: string
    integer, intent(in)          :: start
    character(LEN=*), intent(in) :: substring
    type(varying_string)         :: ins_string

    integer                      :: start_

! Insert a character substring into a character
! string

    start_ = MAX(1, MIN(start, LEN(string)+1))

    ins_string = var_str(string(:start_-1)//substring//string(start_:))

! Finish

    return

  end function insert_CH_CH

!****

  elemental function remove_VS (string, start, finish) result (rem_string)

    type(varying_string), intent(in) :: string
    integer, intent(in), optional    :: start
    integer, intent(in), optional    :: finish
    type(varying_string)             :: rem_string

! Remove a substring from a varying string

    rem_string = remove(char(string), start, finish)

! Finish

    return

  end function remove_VS

!****

  elemental function remove_CH (string, start, finish) result (rem_string)

    character(LEN=*), intent(in)  :: string
    integer, intent(in), optional :: start
    integer, intent(in), optional :: finish
    type(varying_string)          :: rem_string

    integer                       :: start_
    integer                       :: finish_

! Remove a substring from a character string

    if(PRESENT(start)) then
       start_ = MAX(1, start)
    else
       start_ = 1
    endif

    if(PRESENT(finish)) then
       finish_ = MIN(LEN(string), finish)
    else
       finish_ = LEN(string)
    endif

    if(finish_ >= start_) then
       rem_string = var_str(string(:start_-1)//string(finish_+1:))
    else
       rem_string = string
    endif

! Finish

    return

  end function remove_CH

!****

  elemental function replace_VS_VS_auto (string, start, substring) result (rep_string)

    type(varying_string), intent(in) :: string
    integer, intent(in)              :: start
    type(varying_string), intent(in) :: substring
    type(varying_string)             :: rep_string

! Replace part of a varying string with a varying
! substring

    rep_string = replace(char(string), start, MAX(start, 1)+len(substring)-1, char(substring))

! Finish

    return

  end function replace_VS_VS_auto

!****

  elemental function replace_CH_VS_auto (string, start, substring) result (rep_string)

    character(LEN=*), intent(in)     :: string
    integer, intent(in)              :: start
    type(varying_string), intent(in) :: substring
    type(varying_string)             :: rep_string

! Replace part of a character string with a varying
! substring

    rep_string = replace(string, start, MAX(start, 1)+len(substring)-1, char(substring))

! Finish

    return

  end function replace_CH_VS_auto

!****

  elemental function replace_VS_CH_auto (string, start, substring) result (rep_string)

    type(varying_string), intent(in) :: string
    integer, intent(in)              :: start
    character(LEN=*), intent(in)     :: substring
    type(varying_string)             :: rep_string

! Replace part of a varying string with a character
! substring

    rep_string = replace(char(string), start, MAX(start, 1)+LEN(substring)-1, substring)

! Finish

    return

  end function replace_VS_CH_auto

!****

  elemental function replace_CH_CH_auto (string, start, substring) result (rep_string)

    character(LEN=*), intent(in) :: string
    integer, intent(in)          :: start
    character(LEN=*), intent(in) :: substring
    type(varying_string)         :: rep_string

! Replace part of a character string with a character
! substring

    rep_string = replace(string, start, MAX(start, 1)+LEN(substring)-1, substring)

! Finish

    return

  end function replace_CH_CH_auto

!****

  elemental function replace_VS_VS_fixed (string, start, finish, substring) result (rep_string)

    type(varying_string), intent(in) :: string
    integer, intent(in)              :: start
    integer, intent(in)              :: finish
    type(varying_string), intent(in) :: substring
    type(varying_string)             :: rep_string

! Replace part of a varying string with a varying
! substring

    rep_string = replace(char(string), start, finish, char(substring))

! Finish

    return

  end function replace_VS_VS_fixed

!****

!****

  elemental function replace_CH_VS_fixed (string, start, finish, substring) result (rep_string)

    character(LEN=*), intent(in)     :: string
    integer, intent(in)              :: start
    integer, intent(in)              :: finish
    type(varying_string), intent(in) :: substring
    type(varying_string)             :: rep_string

! Replace part of a character string with a varying
! substring

    rep_string = replace(string, start, finish, char(substring))

! Finish

    return

  end function replace_CH_VS_fixed

!****

  elemental function replace_VS_CH_fixed (string, start, finish, substring) result (rep_string)

    type(varying_string), intent(in) :: string
    integer, intent(in)              :: start
    integer, intent(in)              :: finish
    character(LEN=*), intent(in)     :: substring
    type(varying_string)             :: rep_string

! Replace part of a varying string with a character
! substring

    rep_string = replace(char(string), start, finish, substring)

! Finish

    return

  end function replace_VS_CH_fixed

!****

  elemental function replace_CH_CH_fixed (string, start, finish, substring) result (rep_string)

    character(LEN=*), intent(in) :: string
    integer, intent(in)          :: start
    integer, intent(in)          :: finish
    character(LEN=*), intent(in) :: substring
    type(varying_string)         :: rep_string

    integer                      :: start_
    integer                      :: finish_

! Replace part of a character string with a character
! substring

    start_ = MAX(1, start)
    finish_ = MIN(LEN(string), finish)

    if(finish_ < start_) then
       rep_string = insert(string, start_, substring)
    else
       rep_string = var_str(string(:start_-1)//substring//string(finish_+1:))
    endif

! Finish

    return

  end function replace_CH_CH_fixed

!****

  elemental function replace_VS_VS_VS_target (string, target, substring, every, back) result (rep_string)

    type(varying_string), intent(in) :: string
    type(varying_string), intent(in) :: target
    type(varying_string), intent(in) :: substring
    logical, intent(in), optional    :: every
    logical, intent(in), optional    :: back
    type(varying_string)             :: rep_string

! Replace part of a varying string with a varying
! substring, at a location matching a varying-
! string target

    rep_string = replace(char(string), char(target), char(substring), every, back)

! Finish

    return

  end function replace_VS_VS_VS_target

!****

  elemental function replace_CH_VS_VS_target (string, target, substring, every, back) result (rep_string)

    character(LEN=*), intent(in)     :: string
    type(varying_string), intent(in) :: target
    type(varying_string), intent(in) :: substring
    logical, intent(in), optional    :: every
    logical, intent(in), optional    :: back
    type(varying_string)             :: rep_string

! Replace part of a character string with a varying
! substring, at a location matching a varying-
! string target

    rep_string = replace(string, char(target), char(substring), every, back)

! Finish

    return

  end function replace_CH_VS_VS_target

!****

  elemental function replace_VS_CH_VS_target (string, target, substring, every, back) result (rep_string)

    type(varying_string), intent(in) :: string
    character(LEN=*), intent(in)     :: target
    type(varying_string), intent(in) :: substring
    logical, intent(in), optional    :: every
    logical, intent(in), optional    :: back
    type(varying_string)             :: rep_string

! Replace part of a character string with a varying
! substring, at a location matching a character-
! string target

    rep_string = replace(char(string), target, char(substring), every, back)

! Finish

    return

  end function replace_VS_CH_VS_target

!****

  elemental function replace_CH_CH_VS_target (string, target, substring, every, back) result (rep_string)

    character(LEN=*), intent(in)     :: string
    character(LEN=*), intent(in)     :: target
    type(varying_string), intent(in) :: substring
    logical, intent(in), optional    :: every
    logical, intent(in), optional    :: back
    type(varying_string)             :: rep_string

! Replace part of a character string with a varying
! substring, at a location matching a character-
! string target

    rep_string = replace(string, target, char(substring), every, back)

! Finish

    return

  end function replace_CH_CH_VS_target

!****

  elemental function replace_VS_VS_CH_target (string, target, substring, every, back) result (rep_string)

    type(varying_string), intent(in) :: string
    type(varying_string), intent(in) :: target
    character(LEN=*), intent(in)     :: substring
    logical, intent(in), optional    :: every
    logical, intent(in), optional    :: back
    type(varying_string)             :: rep_string

! Replace part of a varying string with a character
! substring, at a location matching a varying-
! string target

    rep_string = replace(char(string), char(target), substring, every, back)

! Finish

    return

  end function replace_VS_VS_CH_target

!****

  elemental function replace_CH_VS_CH_target (string, target, substring, every, back) result (rep_string)

    character(LEN=*), intent(in)     :: string
    type(varying_string), intent(in) :: target
    character(LEN=*), intent(in)     :: substring
    logical, intent(in), optional    :: every
    logical, intent(in), optional    :: back
    type(varying_string)             :: rep_string

! Replace part of a character string with a character
! substring, at a location matching a varying-
! string target

    rep_string = replace(string, char(target), substring, every, back)

! Finish

    return

  end function replace_CH_VS_CH_target

!****

  elemental function replace_VS_CH_CH_target (string, target, substring, every, back) result (rep_string)

    type(varying_string), intent(in) :: string
    character(LEN=*), intent(in)     :: target
    character(LEN=*), intent(in)     :: substring
    logical, intent(in), optional    :: every
    logical, intent(in), optional    :: back
    type(varying_string)             :: rep_string

! Replace part of a varying string with a character
! substring, at a location matching a character-
! string target

    rep_string = replace(char(string), target, substring, every, back)

! Finish

    return

  end function replace_VS_CH_CH_target

!****

  elemental function replace_CH_CH_CH_target (string, target, substring, every, back) result (rep_string)

    character(LEN=*), intent(in)  :: string
    character(LEN=*), intent(in)  :: target
    character(LEN=*), intent(in)  :: substring
    logical, intent(in), optional :: every
    logical, intent(in), optional :: back
    type(varying_string)          :: rep_string

    logical                       :: every_
    logical                       :: back_
    type(varying_string)          :: work_string
    integer                       :: length_target
    integer                       :: i_target

! Handle special cases when LEN(target) == 0. Such
! instances are prohibited by the standard, but
! since this function is elemental, no error can be
! thrown. Therefore, it makes sense to handle them 
! in a sensible manner

    if(LEN(target) == 0) then
       if(LEN(string) /= 0) then
          rep_string = string
       else
          rep_string = substring
       endif
       return
    end if

! Replace part of a character string with a character
! substring, at a location matching a character-
! string target

    if(PRESENT(every)) then
       every_ = every
    else
       every_ = .false.
    endif

    if(PRESENT(back)) then
       back_ = back
    else
       back_ = .false.
    endif

    rep_string = ''

    work_string = string

    length_target = LEN(target)

    replace_loop : do

       i_target = index(work_string, target, back_)

       if(i_target == 0) exit replace_loop

       if(back_) then
          rep_string = substring//extract(work_string, start=i_target+length_target)//rep_string
          work_string = extract(work_string, finish=i_target-1)
       else
          rep_string = rep_string//extract(work_string, finish=i_target-1)//substring
          work_string = extract(work_string, start=i_target+length_target)
       endif

       if(.NOT. every_) exit replace_loop

    end do replace_loop

    if(back_) then
       rep_string = work_string//rep_string
    else
       rep_string = rep_string//work_string
    endif

! Finish

    return

  end function replace_CH_CH_CH_target

!****

  elemental subroutine split_VS (string, word, set, separator, back)

    type(varying_string), intent(inout)         :: string
    type(varying_string), intent(out)           :: word
    type(varying_string), intent(in)            :: set
    type(varying_string), intent(out), optional :: separator
    logical, intent(in), optional               :: back

! Split a varying string into two verying strings

    call split_CH(string, word, char(set), separator, back)

! Finish

    return

  end subroutine split_VS

!****

  elemental subroutine split_CH (string, word, set, separator, back)

    type(varying_string), intent(inout)         :: string
    type(varying_string), intent(out)           :: word
    character(LEN=*), intent(in)                :: set
    type(varying_string), intent(out), optional :: separator
    logical, intent(in), optional               :: back

    logical                                     :: back_
    integer                                     :: i_separator

! Split a varying string into two verying strings

    if(PRESENT(back)) then
       back_ = back
    else
       back_ = .false.
    endif

    i_separator = scan(string, set, back_)

    if(i_separator /= 0) then

       if(back_) then
          word = extract(string, start=i_separator+1)
          if(PRESENT(separator)) separator = extract(string, start=i_separator, finish=i_separator)
          string = extract(string, finish=i_separator-1)
       else
          word = extract(string, finish=i_separator-1)
          if(PRESENT(separator)) separator = extract(string, start=i_separator, finish=i_separator)
          string = extract(string, start=i_separator+1)
       endif

    else

       word = string
       if(PRESENT(separator)) separator = ''
       string = ''

    endif

! Finish

    return

  end subroutine split_CH

end module iso_varying_string
