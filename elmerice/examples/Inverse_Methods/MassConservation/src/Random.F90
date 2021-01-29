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
! *  Authors: F. Gillet-Chaulet
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: Dec. 2020
! *
! * A user function to initialise variables with a random number
! * using an even or normal distribution
! *
! * INPUTS:
! *     - Argument: x(2) real values for a linear transformation of the
! random number: out = r*x(1)+x(2)
! *     - Constants: 
! *       - Random Seed = Integer [OPTIONAL] => give a seed for repetability
! *       - Random Function = String "even" for an even distribution between 0 1
! *                               or "normal" for a normal distribution of 0 mean and unit variance
! *****************************************************************************
      FUNCTION Random(Model,nodenumber,x) RESULT(r)
       USE DefUtils
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp),INTENT(IN) :: x(2)
       REAL(kind=dp) :: r

       INTEGER :: ssize
       INTEGER,ALLOCATABLE :: seed(:)
       LOGICAL, SAVE :: Firsttime=.TRUE.
       LOGICAL :: Found

       CHARACTER(LEN=MAX_NAME_LEN),SAVE :: Rtype

       IF (Firsttime) THEN
         CALL random_seed(size=ssize)
         allocate(seed(ssize))
         
         seed = ListGetInteger( Model % Constants, 'Random Seed',Found )
         IF (Found)  call random_seed( put=seed )
         CALL random_seed(get=seed)

         Rtype = ListGetString( Model % Constants, 'Random Function',UnFoundFatal=.TRUE.)
         deallocate(seed)
         Firsttime=.FALSE.
       ENDIF


       SELECT CASE(Rtype)
       CASE('even')
         r = EvenRandom()
       CASE('normal')
         r = NormalRandom()
       CASE DEFAULT
         CALL FATAL('Random Function','Random Function should be <even> or <normal>')
       END SELECT

       r = r*x(1)+x(2)

       End FUNCTION Random


