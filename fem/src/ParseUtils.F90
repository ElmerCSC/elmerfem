!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! *
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! *
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU Lesser General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juhani Kataja
! *  Email:   juhani.kataja (at) csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Utilities for parsing strings
!> \ingroup Utils
!-------------------------------------------------------------------------------
MODULE ParseUtils

use Messages

IMPLICIT NONE

CONTAINS

!> Find substrings from str separated by sep
!> and produce substring array sub of size n_sub+1
!> so that i:th substring is str(sub(i)+1:sub(i+1)-1)
!> If sub is preallocated and its length is more 
!> than n_subs+1, then make the set sub(n_sub+2) = 0
pure subroutine splitsub(str, sep, subs, n_subs) ! {{{
  character(*), intent(in) :: str
  character, intent(in) :: sep
  integer, intent(inout) :: subs(:)
  integer, intent(out) :: n_subs
  integer, parameter :: state_counting = 0, state_filling = 1, state_finish = 2

  integer :: l, lstr, state, m, len_subs

  lstr = len(str)
  len_subs = size(subs,1)
  state = 0

  do while(state /= state_finish)
    m = 0
    l = 0
    while_l: do while(.true.)
      m = m + 1
      if (state == state_filling) subs(m) = l
      l = consumechar(str(1:lstr), l+1, sep, neg=.true.)
      if (l>lstr) exit while_l
    end do while_l

    if (state == state_filling) subs(m+1) = lstr+1

    ! state transitions: state_counting -> state filling -> state_finish
    if (state == state_counting) then
      if (len_subs < m+1) then
        subs(1) = -1
      elseif (size(subs,1) > m+1) then 
        subs(m+2) = -1
      end if

      state = state_filling
    elseif (state == state_filling) then
      state = state_finish
    end if

  end do
  n_subs = m
end subroutine ! }}}

!> Given string str and offset l0 on it, add 1 to offset l1
!> until str(l1) /= c or, if neg=.true. then add until str(l1) == c.
!> If a non-matching (or matching in case neg=.true.) character is
!> found, return len(str) + 1
pure function consumechar(str,  l0, c, neg) result(l1) ! {{{
  character(len=*), intent(in) :: str
  integer, intent(in) :: l0
  character, intent(in) :: c
  integer :: l1
  logical, optional, intent(in) :: neg
  logical :: neg_

  if (present(neg)) then
    neg_ = neg
  else 
    neg_ = .false.
  end if

  l1 = l0 
  do while(l0 <= l1 .and. l1 <= len(str) .and. ((str(l1:l1) == c) .neqv. neg_))
    l1 = l1 + 1
  end do
end function ! }}}

!> Parse dependency string part and return substring indices for 
!> * variable name(v_sub)
!> * method name (m_sub) if found and 
!> * method dimension (d_sub) if method is found
!> If method is not given then d_sub=m_sub=0
subroutine ParseDepString(str, m_sub, d_sub, v_sub) ! {{{
  integer :: ind1, ind2, ind3, l1, l0
  integer :: m_sub(2), d_sub(2), v_sub(2)
  integer :: sub(3)
  integer :: n_m_sub
  character, parameter :: openrank = '{', closerank = '}'

  character(*) :: str

  m_sub = 0
  d_sub = 0
  v_sub = 0
  l0 = consumechar(str, 1, ' ')
  
  if (l0 > len(str) ) return

  l1 = len(str)

  call splitsub(str(l0:l1), '%', sub, n_m_sub)
  if (n_m_sub > 2) then
    call Fatal('ParseDepString', 'Fatal: can''t chain % symbols.')
  end if

  ! ind1 = 0
  ! if (n_m_sub == 2) ind1 = sub(2)+1

  IF (n_m_sub == 1) THEN
    v_sub(1) = l0 
    v_sub(2) = l1

    ! remove preceding whitespace
    v_sub(1) = v_sub(1) + consumechar(str(v_sub(1):v_sub(2)), 1, ' ') - 1

    l0 = l1+2
  ELSE
    v_sub(1) = l0+sub(1)
    v_sub(2) = l0+sub(2)-2

    ! remove preceding whitespace
    v_sub(1) = v_sub(1) + consumechar(str(v_sub(1):v_sub(2)), 1, ' ') - 1

    if (v_sub(1) == v_sub(2)) then
      call Fatal('ParseDepString', 'Fatal: Empty variable name given')
    end if
    ind2 = consumechar(str(l0:l1), 1, openrank, .true.) ! this is ok

    call splitsub(str(v_sub(2)+2:l1), openrank, sub, n_m_sub)
    m_sub(1) = v_sub(2) + 2 + sub(1)
    m_sub(2) = v_sub(2) + 2 + sub(2)-2

    ! remove preceding whitespace
    m_sub(1) = m_sub(1) + consumechar(str(m_sub(1):m_sub(2)), 1, ' ') - 1

    if (ind2 == 0) then
      call Fatal('ParseDepString', 'Can''t find opening brace in method string > '// &
          TRIM(str(l0:l1))//' < for dependent variable:['//TRIM(str(v_sub(1):v_sub(2)))//']')
    end if
    if (m_sub(1) > m_sub(2)) then
      Call Fatal('ParseDepString','Method is missing for dependent variable:['//TRIM(str(v_sub(1):v_sub(2)))//']' )
    end if

    call splitsub(str(m_sub(2) + 2:l1), closerank, sub, n_m_sub)
    if (n_m_sub == 1) THEN
      Call Fatal('ParseDepString','Can''t find closing brace in method string > '// &
          TRIM(str(l0:l1))//' < for dependent variable:['//TRIM(str(v_sub(1):v_sub(2)))//']' )
    END IF
    d_sub(1) = m_sub(2) + 2 + sub(1) 
    d_sub(2) = m_sub(2) + 2 + sub(2) - 2

    ! remove preceding whitespace
    d_sub(1) = d_sub(1) + consumechar(str(d_sub(1):d_sub(2)), 1, ' ') - 1

    l0 = d_sub(2) + 3
  END IF
END subroutine ! }}}
END MODULE ParseUtils
