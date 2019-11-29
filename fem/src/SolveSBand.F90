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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!-------------------------------------------------------------------------------
!> Calls LAPACK symmetric band matrix solvers.
!-------------------------------------------------------------------------------
       SUBROUTINE SolveSBandLapack( N,M,A,X,Subband,Band )

       USE Types

       IMPLICIT NONE

       INTEGER :: N,M,Subband,Band
       REAL(KIND=dp) :: A(Band,N),X(N,M)

       INTEGER :: IPIV(N),iINFO

       IF ( N .LE. 0 ) RETURN

       iINFO = 0
       CALL DPBTRF( 'L',N,Subband,A,Band,iINFO )
        IF ( iinfo /= 0 ) THEN
         PRINT*,'ERROR: SolveSymmetricBand: singular matrix. LAPACK DPBTRF iinfo: ',iinfo
          STOP 1
        END IF

       iINFO = 0
       CALL DPBTRS( 'L',N,Subband,M,A,Band,X,N,iINFO )
        IF ( iinfo /= 0 ) THEN
         PRINT*,'ERROR: SolveSymmetricBand: singular matrix. LAPACK DPBTRS iinfo: ',iinfo
          STOP 1
        END IF

       END


!-------------------------------------------------------------------------------
!> Calls LAPACK complex symmetric band matrix solvers.
!-------------------------------------------------------------------------------
       SUBROUTINE SolveComplexSBandLapack( N,M,A,X,Subband,Band )

       USE Types

       IMPLICIT NONE

       INTEGER :: N,M,Subband,Band
       COMPLEX(KIND=dp) :: A(Band,N),X(N,M)

       INTEGER :: IPIV(N),iINFO

       IF ( N .LE. 0 ) RETURN

       iINFO = 0
       CALL ZPBTRF( 'L',N,Subband,A,Band,iINFO )
        IF ( iinfo /= 0 ) THEN
         PRINT*,'ERROR: SolveSymmetricBand: singular matrix. LAPACK ZPBTRF iinfo: ',iinfo
          STOP 1
        END IF

       iINFO = 0
       CALL ZPBTRS( 'L',N,Subband,M,A,Band,X,N,iINFO )
        IF ( iinfo /= 0 ) THEN
         PRINT*,'ERROR: SolveSymmetricBand: singular matrix. LAPACK ZPBTRS iinfo: ',iinfo
          STOP 1
        END IF

       END


!> \}
