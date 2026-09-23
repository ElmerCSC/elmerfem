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
!>  Create a 2D mesh given by its bed and surface topography 
!>  By deforming a square 1x1 mesh
!> 
!>  Input :
!>     NameMesh of the grd mesh (a square 1 x 1)
!>     x_start end x_end 
!>     yi=surf(xi) and yi=bed(xi) are in NameMesh_surf.dat and NameMesh_bed.dat
!> 
!>  Output :
!>     NameMesh/mesh.nodes
   PROGRAM  MshGlacier
       USE GeneralUtils
       USE types
       USE DefUtils
!
IMPLICIT NONE
!
REAL(KIND=dp) :: x0, x1, x, y, z, hmin, xnew, ynew, zb, zs
REAL(KIND=dp) :: Xconstraint, XCubic
REAL(KIND=dp), ALLOCATABLE :: xbed(:), ybed(:), xsurf(:), ysurf(:)
REAL(KIND=dp), TARGET, ALLOCATABLE ::  y2s(:), y2b(:)
REAL(KIND=dp), DIMENSION(:), POINTER :: y2s_p, y2b_p 
REAL(KIND=dp), ALLOCATABLE  :: xnode(:), ynode(:)                       
LOGICAL :: Constraint=.FALSE., Cubic=.FALSE.
CHARACTER :: NameMsh*20, NameSurf*20, NameBed*20, Rien*1
INTEGER :: NtN, i, j, NptS, NptB, n, NtNx, xi
INTEGER :: ioInput, ioHeader, ioSurf, ioBed, ioNodes
!
!
!  
!
!  Read data input from mesh_input.dat
!
!
      OPEN(NEWUNIT=ioInput,file="mesh_input.dat")
      READ(ioInput,*)Rien
      READ(ioInput,*)NameMsh
      READ(ioInput,*)Rien
      READ(ioInput,*)XConstraint
! Xconstraint = 0 -> xnode given by the mesh
! Xconstraint = 1 -> xnode given by the dataset
      IF (XConstraint > 0.5) Constraint = .TRUE.
! xCubic = 0 -> Linear interpolation
! xCubic = 1 -> Cubic Spline
      READ(ioInput,*)Rien
      READ(ioInput,*)xCubic
      IF (XCubic > 0.5) Cubic = .TRUE.
      READ(ioInput,*)Rien
      READ(ioInput,*)NameSurf
      READ(ioInput,*)Rien
      READ(ioInput,*)NptS
      READ(ioInput,*)Rien
      READ(ioInput,*)NameBed
      READ(ioInput,*)Rien
      READ(ioInput,*)NptB
      READ(ioInput,*)Rien
      READ(ioInput,*)hmin
      IF (.Not.Constraint) THEN
        READ(ioInput,*)Rien
        READ(ioInput,*)x0
        READ(ioInput,*)Rien
        READ(ioInput,*)x1
      END IF
      CLOSE(ioInput)

      ALLOCATE(xsurf(NptS), ysurf(NptS), y2s(NptS))
      ALLOCATE(xbed(NptB), ybed(NptB), y2b(NptB))

      OPEN(NEWUNIT=ioHeader,file=TRIM(NameMsh)//"/mesh.header")
        READ(ioHeader,1000)NtN
      CLOSE(ioHeader)
      ALLOCATE(xnode(NtN), ynode(NtN))

      OPEN(NEWUNIT=ioSurf,file=TRIM(NameSurf))
        READ(ioSurf,*)(xsurf(i), ysurf(i), i=1,NptS)
      CLOSE(ioSurf)

      OPEN(NEWUNIT=ioBed,file=TRIM(NameBed))
        READ(ioBed,*)(xbed(i), ybed(i), i=1,NptB)
      CLOSE(ioBed)
      
      IF (.Not.Constraint) THEN
        IF (((MINVAL(xbed)>x0) .OR. (MAXVAL(xbed)<x1)) ) THEN
           WRITE(*,*)'MUST BE : x0> MIN(xbed) AND x1 < MAX(xbed)',& 
                  &x0,MINVAL(xbed),x1,MAXVAL(xbed)
           STOP
        END IF
        IF (((MINVAL(xsurf)>x0) .OR. (MAXVAL(xsurf)<x1)) ) THEN
           WRITE(*,*)'MUST BE : x0 > MIN(xsurf) AND x1 < MAX(xsurf)'
           STOP
        END IF
      END IF
    
      IF (Cubic) THEN
        CALL CubicSpline(NptS,xsurf,ysurf,y2s)
        CALL CubicSpline(Nptb,xbed,ybed,y2b)
      END IF
      
      OPEN(NEWUNIT=ioNodes,file=TRIM(NameMsh)//"/mesh.nodes")
      READ(ioNodes,*)(N, j, xnode(i), ynode(i), z, i=1,NtN)
      REWIND(ioNodes)


      IF (Constraint) THEN
              WRITE(*,*)'Will have to check if number of nodes in x direction &
                    &equals number of data points'
! First test
          IF (NptS/=NptB) THEN
                  WRITE(*,*)'Surface and Bed data must have the same number of &
                  &points'
                  STOP
          END IF
          NtNx = 0
          DO n=1,NtN
             IF (ABS(ynode(n))<1.0e-4) NtNx = NtNx + 1 
          END DO
          IF (NptS/=NtNx) THEN
                  WRITE(*,*)'Mesh must have',NptS,' nodes in x direction' 
                  STOP
          END IF
          x0 = MINVAL(xbed)
          x1 = MAXVAL(xbed)
      END IF

      DO n=1, NtN

        x = xnode(n)
        y = ynode(n)


        IF (Constraint) THEN
           xi = ANINT(x / (1.0_dp / (NptS-1.0_dp))) + 1
           x = (xsurf(xi)-x0)/(x1-x0)
        END IF

        xnew = x0 + x * (x1 - x0) 

        IF (Cubic) THEN
          y2s_p => y2s
          y2b_p => y2b
          zs = InterpolateCurve(xsurf,ysurf,xnew,y2s_p)
          zb = InterpolateCurve(xbed,ybed,xnew,y2b_p)
        ELSE
          zs = InterpolateCurve(xsurf,ysurf,xnew)
          zb = InterpolateCurve(xbed,ybed,xnew)
        END IF

        ynew = zb + y * MAX((zs - zb),hmin) 
        
        WRITE(ioNodes,1200)N,j,xnew,ynew,z
      END DO
      WRITE(*,*)'END WITH NO TROUBLE ...'
!
1000 FORMAT(I6)
1200 FORMAT(i5,2x,i5,3(2x,e22.15)) 

End Program MshGlacier
