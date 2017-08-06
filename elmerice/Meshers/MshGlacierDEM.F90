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
!> Create a 3D mesh given by its bed and surface topography 
!> given as structued Grid DEMs (x, y, z)
!> Add some testing to ensure the file structure does correspond to what is 
!> indicated in mesh_input.dat
! *****************************************************************************
PROGRAM  MshGlacierDEM 
USE types
USE DefUtils
 
IMPLICIT NONE
REAL(KIND=dp)  ::  x, y, z 
REAL(KIND=dp), ALLOCATABLE :: xb(:), yb(:), zb(:), xs(:), & 
                       ys(:), zs(:)
REAL(KIND=dp), ALLOCATABLE  :: xnode(:), ynode(:), znode(:)                       
INTEGER, ALLOCATABLE :: Node(:)
CHARACTER :: NameMsh*30, Rien*1, NameSurf*30, NameBed*30
Character(len=10) :: iNp, iNN
REAL(KIND=dp)  :: dsx, dsy, dbx, dby, x1, x2, y1, y2, zi(2,2)
REAL(KIND=dp)  :: xs0, ys0, xb0, yb0, zbed, zsurf, znew, Rmin, R
REAL(KIND=dp)  :: lsx, lsy, lbx, lby, hmin 
INTEGER :: NtN, i, j, Ns, Nsx, Nsy, Nb, Nbx, Nby, n, Npt, ix, iy, imin, is, ib
Integer k, Np
Logical :: Serial = .False.

!
!  Read data input from mesh_input.dat
!
OPEN(10,file="mesh_input.dat")
READ(10,*)Rien
READ(10,*)NameMsh
READ(10,*)Rien
READ(10,*)NameSurf
READ(10,*)Rien
READ(10,*)Nsx, Nsy
READ(10,*)Rien
READ(10,*)xs0, ys0
READ(10,*)Rien
READ(10,*)lsx, lsy
READ(10,*)Rien
READ(10,*)NameBed
READ(10,*)Rien
READ(10,*)Nbx, Nby
READ(10,*)Rien
READ(10,*)xb0, yb0
READ(10,*)Rien
READ(10,*)lbx, lby
READ(10,*)Rien
READ(10,*)hmin
READ(10,*)Rien
READ(10,*)Np   
CLOSE(10)

Ns = Nsx*Nsy
Nb = Nbx*Nby
ALLOCATE(xs(Ns), ys(Ns), zs(Ns))
ALLOCATE(xb(Nb), yb(Nb), zb(Nb))
      
dsx = lsx / (Nsx-1.0)
dsy = lsy / (Nsy-1.0)
dbx = lbx / (Nbx-1.0)
dby = lby / (Nby-1.0)

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
      
WRITE(*,*)'Ns, Nb',Ns,Nb
WRITE(*,*)'dx, dy',dsx,dsy,dbx,dby

!-------------------------------------------------------------
! Load Bedrock and Surface DEMs and make some verifications
!-------------------------------------------------------------
OPEN(10,file=TRIM(NameSurf))
READ(10,*)(xs(i), ys(i), zs(i), i=1,Ns)
CLOSE(10)

OPEN(10,file=TRIM(NameBed))
READ(10,*)(xb(i), yb(i), zb(i), i=1,Nb)
CLOSE(10)

k = 0 
DO j = 1, Nsy
   y = ys0 + dsy*(j-1)
   DO i = 1, Nsx 
      k = k + 1
      x = xs0 + dsx*(i-1)
      IF ((ABS(x-xs(k))>1.0e-6*dsx).OR.(ABS(y-ys(k))>1.0e-6*dsy)) THEN
         WRITE(*,*)'Structure of the DEM is not conforming to what is in mesh_input.dat for Surface DEM' 
         WRITE(*,*)'Found that point ',k,' coordinate is ',xs(k),ys(k),' whereas it should be ',x,y 
         STOP
      END IF
   END DO
END DO
k = 0
DO j = 1, Nby
   y = yb0 + dby*(j-1)
   DO i = 1, Nbx 
      k = k + 1
      x = xb0 + dbx*(i-1)
      IF ((ABS(x-xb(k))>1.0e-6*dbx).OR.(ABS(y-yb(k))>1.0e-6*dby)) THEN
         WRITE(*,*)'Structure of the DEM is not conforming to what is in mesh_input.dat for Surface DEM' 
         WRITE(*,*)'Found that point ',k,' coordinate is ',xb(k),yb(k),' whereas it should be ',x,y 
         STOP
      END IF
   END DO
END DO
!-------------------------------------------------------------
! Loop over the Np partitions
!-------------------------------------------------------------
DO k = 1, Np
   WRITE(*,*)k,' of ', Np, ' partitions'

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
 
      OPEN(11,file=TRIM(NameMsh)//"/partitioning."//TRIM(iNp)//"/part."//TRIM(iNN)//".header")
   ELSE
      OPEN(11,file=TRIM(NameMsh)//"/mesh.header")
   END IF 

   READ(11,*)NtN
   CLOSE(11)

   ALLOCATE (Node(NtN), xnode(NtN), ynode(NtN), znode(NtN))
   WRITE(*,*)'Part ', k, ' NtN = ', NtN
        
   IF (.NOT.Serial) THEN
      OPEN(12,file=TRIM(NameMsh)//"/partitioning."//TRIM(iNp)//"/part."//TRIM(iNN)//".nodes")
   ELSE
      OPEN(12,file=TRIM(NameMsh)//"/mesh.nodes")
   END IF
   READ(12,*)(Node(i), j, xnode(i), ynode(i), znode(i), i=1,NtN)
   REWIND(12)

! Make some verifications that all nodes are included in the DEMs
   IF (((MINVAL(xnode)<MINVAL(xs)).OR.(MAXVAL(xnode)>MAXVAL(xs))).OR. &
       ((MINVAL(ynode)<MINVAL(ys)).OR.(MAXVAL(ynode)>MAXVAL(ys)))) THEN  
      WRITE(*,*)'Some nodes are outside of the Surface DEM'
      WRITE(*,*)'DEM xmin, ymin, xmax, ymax: ', &
                       MINVAL(xs),MINVAL(ys),MAXVAL(xs),MAXVAL(ys)
      WRITE(*,*)'Mesh xmin, ymin, xmax, ymax: ', &
                       MINVAL(xnode),MINVAL(ynode),MAXVAL(xnode),MAXVAL(ynode)
      STOP
   END IF

   IF (((MINVAL(xnode)<MINVAL(xb)).OR.(MAXVAL(xnode)>MAXVAL(xb))).OR. &
       ((MINVAL(ynode)<MINVAL(yb)).OR.(MAXVAL(ynode)>MAXVAL(yb)))) THEN  
      WRITE(*,*)'Some nodes are outside of the bedrock DEM'
      WRITE(*,*)'DEM xmin, ymin, xmax, ymax: ', &
                       MINVAL(xb),MINVAL(yb),MAXVAL(xb),MAXVAL(yb)
      WRITE(*,*)'Mesh xmin, ymin, xmax, ymax: ', &
                       MINVAL(xnode),MINVAL(ynode),MAXVAL(xnode),MAXVAL(ynode)
      STOP
   END IF

   DO n=1, NtN
      x = xnode(n)
      y = ynode(n)
      z = znode(n)

      ! Find zbed for that point from the Bedrock MNT 
      ix = INT((x-xb0)/dbx)+1
      IF (x.EQ.(xb0+lbx)) ix=ix-1
      iy = INT((y-yb0)/dby)+1
      IF (y.EQ.(yb0+lby)) iy=iy-1
      ib = Nbx * (iy - 1) + ix
        
      x1 = xb(ib)
      x2 = xb(ib+1)
      y1 = yb(ib)
      y2 = yb(ib + Nbx)
        
      zi(1,1) = zb(ib)
      zi(2,1) = zb(ib+1)
      zi(2,2) = zb(ib + Nbx + 1)
      zi(1,2) = zb(ib + Nbx)

      IF ((zi(1,1)<-9990.0).OR.(zi(1,2)<-9990.0).OR.(zi(2,1)<-9990.0).OR.(zi(2,2)<-9990.0)) THEN
         IF ((zi(1,1)<-9990.0).AND.(zi(1,2)<-9990.0).AND.(zi(2,1)<-9990.0).AND.(zi(2,2)<-9990.0)) THEN
            ! Find the nearest point available
            Rmin = 9999.0
            DO i=1, Nb
               IF (zb(i)>-9990.0) THEN
                  R = SQRT((x-xb(i))**2.0+(y-yb(i))**2.0)
                  IF (R<Rmin) THEN
                     Rmin = R
                     imin = i
                  END IF
               END IF
            END DO
            zbed = zb(imin)
            WRITE(*,*)'No data in the DEM cell ',ix,iy,' . Using the closest point at distance ', Rmin, zbed 
            
         ELSE
            ! Mean value over the available data
            zbed = 0.0
            Npt = 0
            DO i=1, 2
               DO J=1, 2
                  IF (zi(i,j) > -9990.0) THEN 
                     zbed = zbed + zi(i,j)
                     Npt = Npt + 1
                  END IF   
               END DO
            END DO
            zbed = zbed / Npt
            WRITE(*,*)'Missing data in cell ',ix,iy,'. Value ',zbed,' computed using ', Npt, 'points.'
         END IF
      ELSE
         zbed = (zi(1,1)*(x2-x)*(y2-y)+zi(2,1)*(x-x1)*(y2-y)+zi(1,2)*(x2-x)*(y-y1)+zi(2,2)*(x-x1)*(y-y1))/(dbx*dby)
      END IF

      ! Find zsurf for that point from the Surface MNT
      ix = INT((x-xs0)/dsx)+1
      IF (x.EQ.(xs0+lsx)) ix=ix-1
      iy = INT((y-ys0)/dsy)+1
      IF (y.EQ.(ys0+lsy)) iy=iy-1
      is = Nsx * (iy - 1) + ix
        
      x1 = xs(is)
      x2 = xs(is+1)
      y1 = ys(is)
      y2 = ys(is + Nsx)
        
      zi(1,1) = zs(is)
      zi(2,1) = zs(is+1)
      zi(2,2) = zs(is + Nsx + 1)
      zi(1,2) = zs(is + Nsx)
        
      IF ((zi(1,1)<-9990.0).OR.(zi(1,2)<-9990.0).OR.(zi(2,1)<-9990.0).OR.(zi(2,2)<-9990.0)) THEN
         IF ((zi(1,1)<-9990.0).AND.(zi(1,2)<-9990.0).AND.(zi(2,1)<-9990.0).AND.(zi(2,2)<-9990.0)) THEN
            ! Find the nearest point available
            Rmin = 9999.0
            DO i=1, Ns
               IF (zs(i)>-9990.0) THEN
                  R = SQRT((x-xs(i))**2.0+(y-ys(i))**2.0)
                  IF (R<Rmin) THEN
                     Rmin = R
                     imin = i
                  END IF
               END IF
            END DO
            zsurf = zs(imin)
            WRITE(*,*)'No data in the DEM cell ',ix,iy,' . Using the closest point at distance ', Rmin, zsurf
            
         ELSE
            ! Mean value over the available data
            zsurf = 0.0
            Npt = 0
            DO i=1, 2
               DO J=1, 2
                  IF (zi(i,j) > -9990.0) THEN 
                     zsurf = zsurf + zi(i,j)
                     Npt = Npt + 1
                  END IF   
               END DO
            END DO
            zsurf = zsurf / Npt
            WRITE(*,*)'Missing data in cell ',ix,iy,'. Value ',zsurf,' computed using ', Npt, 'points.'
         END IF
      ELSE
         zsurf = (zi(1,1)*(x2-x)*(y2-y)+zi(2,1)*(x-x1)*(y2-y)+zi(1,2)*(x2-x)*(y-y1)+zi(2,2)*(x-x1)*(y-y1))/(dsx*dsy)
      END IF      
        
      znew = zbed + z * MAX((zsurf - zbed),hmin) 
        
      WRITE(12,1200)Node(n),j,x,y,znew
   END DO ! NtN
   CLOSE(12)
   DEALLOCATE (Node, xnode, ynode, znode)

END DO ! Np
   
WRITE(*,*)'END WITH NO TROUBLE ...'

!
1000 FORMAT(I6)
1200 FORMAT(i10,2x,i5,3(2x,e22.15)) 

END PROGRAM MshGlacierDEM

! --------------------------------------------------------

