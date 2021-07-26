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

!----------------------------------------------------------------------------
!>  Call LAPACK band matrix solvers.
!----------------------------------------------------------------------------


      SUBROUTINE SolveBandLapack( N,M,A,X,Subband,Band )

      USE Types
      IMPLICIT NONE
 

      INTEGER :: N,M,Subband,Band
      REAL(KIND=dp) :: A(Band,N),X(M,N)

      INTEGER :: IPIV(N),iINFO

      IF ( N .LE. 0 ) RETURN

      iINFO = 0
      CALL DGBTRF( N,N,Subband,Subband,A,Band,IPIV,iINFO )
      IF ( iinfo /= 0 ) THEN
         PRINT*,'ERROR: SolveBand: singular matrix. LAPACK DGBTRF info: ',iinfo
         STOP 1
      END IF

      iINFO = 0
      CALL DGBTRS( 'N',N,Subband,Subband,M,A,Band,IPIV,X,N,iINFO )
        IF ( iinfo /= 0 ) THEN
        PRINT*,'ERROR: SolveBand: singular matrix. LAPACK DGBTRS info: ',iinfo
          STOP 1
        END IF

      END


      SUBROUTINE SolveComplexBandLapack( N,M,A,X,Subband,Band )

      USE Types
      IMPLICIT NONE

      INTEGER :: N,M,Subband,Band
      COMPLEX(KIND=dp) :: A(Band,N),X(M,N)

      INTEGER :: IPIV(N),iINFO

      IF ( N .LE. 0 ) RETURN

      iINFO = 0
      CALL ZGBTRF( N,N,Subband,Subband,A,Band,IPIV,iINFO )
        IF ( iinfo /= 0 ) THEN
        PRINT*,'ERROR: SolveBand: singular matrix. LAPACK ZGBTRF info: ',iinfo
          STOP 1
        END IF

      iINFO = 0
      CALL ZGBTRS( 'N',N,Subband,Subband,M,A,Band,IPIV,X,N,iINFO )
        IF ( iinfo /= 0 ) THEN
        PRINT*,'ERROR: SolveBand: singular matrix. LAPACK ZGBTRS info: ',iinfo
          STOP 1
        END IF

      END

!> \}
	  
	  
