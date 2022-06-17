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
! *  Authors: Olivier Gagliardini
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!>  Create a Synthetic 2D or 3D mesh starting from a 
!>  1m thick mesh using the right contour of the final mesh
!>  and bedrock and surface elevations given as functions
!>  Work for partitioned mesh
!>  number of partitions in mesh_info.in
!> 
!>  Need to be compiled with 
!>     2D :  fbed(x) and fsurf(x) 
!>     3D :  fbed(x,y) and fsurf(x,y) 
!> 
PROGRAM  MshGlacierSynthetic 
USE types
!
IMPLICIT NONE 
!
INTEGER :: NtN
!
REAL(KIND=dp)  ::   znew, zb, zs
REAL(KIND=dp) :: fbed, fsurf
REAL(KIND=dp), ALLOCATABLE :: x(:), y(:), z(:)
INTEGER, ALLOCATABLE :: Node(:)
Character :: Rien*1
Character(len=10) :: iNp, iNN
CHARACTER :: NameMsh*30
Integer i, j, k, N, Np, dim
Logical :: Serial = .False.

!
! Read the name mesh and number of partitions
!
OPEN(10,file="mesh_info.in")
   READ(10,*)Rien
   READ(10,*)NameMsh 
   READ(10,*)Rien
   READ(10,*)Np 
CLOSE(10)
      
IF (Np==1) Serial=.True.

IF (.NOT.Serial) THEN
   IF (Np < 10) THEN
      WRITE(iNp,'(i1.1)')Np
   ELSE IF (Np < 100) THEN
      WRITE(iNp,'(i2.2)')Np
   ELSE IF (Np < 1000) THEN
      WRITE(iNp,'(i3.3)')Np
   ELSE IF (Np < 10000) THEN
      WRITE(iNp,'(i4.4)')Np
   ELSE
      WRITE(*,*)'Work for a number of partitions < 1000'
      STOP
   END IF
   iNp = ADJUSTL(iNp)
END IF

DO k = 1, Np
   WRITE(*,*)k,' of ', Np

   IF (.NOT.Serial) THEN
      IF (k < 10) THEN
         WRITE(iNN,'(i1.1)')k
      ELSE IF (k < 100) THEN
         WRITE(iNN,'(i2.2)')k 
      ELSE IF (k < 1000) THEN
         WRITE(iNN,'(i3.3)')k 
      ELSE IF (k < 10000) THEN
         WRITE(iNN,'(i4.4)')k 
      END IF
      iNN = ADJUSTL(iNN)
      WRITE(*,*)'iNN',iNN
 
      OPEN(11,file=TRIM(NameMsh)//"/partitioning."//TRIM(iNp)//"/part."//TRIM(iNN)//".header")
   ELSE
      OPEN(11,file=TRIM(NameMsh)//"/mesh.header")
   END IF 

   READ(11,*)NtN
   CLOSE(11)

   ALLOCATE (Node(NtN), x(NtN), y(NtN), z(NtN))
   WRITE(*,*)'Part ', k, ' NtN = ', NtN
        
   IF (.NOT.Serial) THEN
      OPEN(12,file=TRIM(NameMsh)//"/partitioning."//TRIM(iNp)//"/part."//TRIM(iNN)//".nodes")
   ELSE
      OPEN(12,file=TRIM(NameMsh)//"/mesh.nodes")
   END IF
   DO i=1,NtN
      READ(12,*)Node(i),j,x(i),y(i),z(i)
   END DO
    
   ! Check if we have a 2D or a 3D geometry
   IF (k==1) THEN
      IF (ANY(ABS(z)>AEPS)) THEN ! 3D Case
         IF (ANY(z > 1.0_dp).OR.ANY(z<0.0_dp)) THEN
            WRITE(*,*)'For 3D geometry, the initial mesh must fulfil 0 < z < 1'
            STOP
         END IF
         dim = 3
         WRITE(*,*)'Initial mesh verified, found to be a 3D geometry'
      ELSE ! 2D Case
         IF (ANY(y > 1.0_dp).OR.ANY(y<0.0_dp)) THEN
            WRITE(*,*)'For 2D geometry, the initial mesh must fulfil 0 < y < 1'
            STOP
         END IF
         dim = 2
         WRITE(*,*)'Initial mesh verified, found to be a 2D geometry'
      END IF
   END IF

   REWIND(12)

   DO i=1,NtN
      IF (dim==2) THEN
         zs = fsurf(x(i)) 
         zb = fbed(x(i))
      ELSE 
         zs = fsurf(x(i),y(i)) 
         zb = fbed(x(i),y(i))
      END IF 
      IF ((zs-zb).LE.0.0) THEN 
         WRITE(*,*)'NEGATIVE OR NULL THICKNESS!!!'
         STOP
      END IF 

      IF (dim==2) THEN
         znew = zb + y(i)*(zs - zb) 
         WRITE(12,1200)Node(i),j,x(i),znew,z(i)
      ELSE
         znew = zb + z(i)*(zs - zb) 
         WRITE(12,1200)Node(i),j,x(i),y(i),znew
      END IF
   END DO
   CLOSE(12)

   DEALLOCATE (Node, x, y, z)
END DO 
!
!
1000 Format(a1)
1001 Format(28x,3(i3,x))
1002 Format(x,e14.8)

1200 Format(i7,2x,i5,3(2x,e22.15)) 

END PROGRAM MshGlacierSynthetic
