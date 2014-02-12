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

       SUBROUTINE SolveLapack( N,A,x )

       INTEGER  N,IPIV(N)
       DOUBLE PRECISION  A(n*n),x(n)

       IF ( N .LE. 0 ) RETURN
       CALL DGETRF( N,N,A,N,IPIV,INFO )
       IF ( info .NE. 0 ) PRINT*,'DGETRF: ', info

       CALL DGETRS( 'N',N,1,A,N,IPIV,X,N,INFO )
       IF ( info .NE. 0 ) PRINT*,'DGETRS: ', info

       END

! ******************************************************************************
