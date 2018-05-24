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
! *  Authors: 
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Interpolate data given on a regular 2D regular grid in an ASCII file (x y Value)
!    in the mesh nodes using bilinear interpolation
!    The data are ordered such that   
!    x1 y1 val11
!    x2 y1 val21
!    ...
!    xn y1 valn1
!    x1 y2 val12
!    ...
!    xn yn valnn 
!    
!    The grid is described by giving:
!    (x0, y0) the left-bottom corner coordinate
!    (lx, ly) the x and y lengths of the covered domain
!    (Nx, Ny) the number of cells in x and y directions 
!    No data are given by -9999 with a tolerance of 0.001
!    These can be over-ridden in the sif by 'no data' and 'no data tol'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Grid2DInterpolator( Model,Solver,dt,TransientSimulation )

   USE DefUtils

   IMPLICIT NONE
   TYPE(Solver_t), TARGET :: Solver
   TYPE(Model_t) :: Model
   REAL(KIND=dp) :: dt
   LOGICAL :: TransientSimulation

   TYPE(ValueList_t), POINTER :: Params
   TYPE(Variable_t), POINTER :: Var
   REAL(KIND=dp), POINTER :: Values(:)
   INTEGER, POINTER :: Perm(:)

   REAL(KIND=DP) :: Rmin, Rmax
   REAL(KIND=DP) :: x, y, z, x0, y0, lx, ly, dx, dy
   REAL(KIND=DP), ALLOCATABLE :: xb(:), yb(:), zb(:), xbaux(:), ybaux(:), zbaux(:)
   REAL(KIND=dp) :: noDataVal, noDataTol, posTol
   REAL(KIND=dp), PARAMETER :: noDataValDefault = -9999.0, noDataTolDefault = 0.001, posTolDefault= 1.0D-06

   INTEGER,parameter :: io=20
   INTEGER :: ok, Nx, Ny, Nb, Nbaux, OutNode
   INTEGER :: i, j, k, l, kmin, NoVar

   CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, DataF
   CHARACTER(LEN=MAX_NAME_LEN) :: Name, FName, ParaName
   CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='Grid2DInterpolator'

   LOGICAL :: GotVar, Found, InvertOrder, FillIn, UnFoundFatal=.TRUE.

   NULLIFY(Params,Var,Values,Perm)

   Params => GetSolverParams()

   ! Read variable to initialize and Data
   NoVar=0
   GotVar=.True.

   DO WHILE(GotVar)
      NoVar = NoVar + 1
      WRITE (Name,'(A,I0)') 'Variable ',NoVar

      VariableName = ListGetString( Params, TRIM(Name), GotVar )
      IF (.NOT.GotVar) EXIT

      Var => VariableGet(Model %  Mesh % Variables, VariableName,UnFoundFatal=UnFoundFatal )
      Values => Var % Values
      Perm => Var % Perm

      WRITE (FName,'(A,I0,A)') 'Variable ',NoVar,' Data File'
      DataF = ListGetString( Params, TRIM(FName), Found, UnFoundFatal )

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' x0'
      x0 = ListGetConstReal( Params, TRIM(ParaName), Found,UnFoundFatal=UnFoundFatal)

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' y0'
      y0 = ListGetConstReal( Params, TRIM(ParaName), Found,UnFoundFatal=UnFoundFatal )
            
      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' lx'
      lx = ListGetConstReal( Params, TRIM(ParaName), Found,UnFoundFatal=UnFoundFatal )
            
      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' ly'
      ly = ListGetConstReal( Params, TRIM(ParaName), Found,UnFoundFatal=UnFoundFatal )

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' Nx'
      Nx = ListGetInteger( Params, TRIM(ParaName), Found ,UnFoundFatal=UnFoundFatal)

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' Ny'
      Ny = ListGetInteger( Params, TRIM(ParaName), Found ,UnFoundFatal=UnFoundFatal)

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' Invert'
      InvertOrder = GetLogical( Params, TRIM(ParaName), Found )
      IF (.NOT.Found) THEN
         InvertOrder = .FALSE.
      END IF
      IF (InvertOrder) THEN
         WRITE(message,'(A,A,I0)')'Inverting order (row major) for variable ', 'Variable ',NoVar
         CALL INFO(Trim(SolverName),Trim(message),Level=1)
      END IF

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' Fill'
      FillIn = GetLogical( Params, TRIM(ParaName), Found )
      IF (.NOT.Found) THEN
         FillIn = .FALSE.
      END IF
      IF (FillIn) THEN
         WRITE(message,'(A,A,I0)')'Filling empty entries for ', 'Variable ',NoVar
         CALL INFO(Trim(SolverName),Trim(message),Level=1)
      END IF

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' no data'
      noDataVal = ListGetConstReal( Params, TRIM(ParaName), Found )
      IF (.NOT.Found) then
         noDataVal = noDataValDefault
         WRITE(message,'(A,A,A,e14.8)')'Keyword <',Trim(ParaName), & 
              '> not found, using default ',noDataValDefault
         CALL INFO(SolverName, Message, Level=3)
      END IF

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' no data tol'
      noDataTol = ListGetConstReal( Params, TRIM(ParaName), Found )
      IF (.NOT.Found) then
         noDataTol = noDataTolDefault
         WRITE(message,'(A,A,A,e14.8)')'Keyword <',Trim(ParaName), & 
              '> not found, using default ',noDataTolDefault
         CALL INFO(SolverName, Message, Level=3)
      END IF

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' position tol'
      posTol = ListGetConstReal( Params, TRIM(ParaName), Found )
      IF (.NOT.Found) then
         posTol = posTolDefault
         WRITE(message,'(A,A,A,e14.8)')'Keyword <',Trim(ParaName), & 
              '> not found, using default ',posTolDefault
         CALL INFO(SolverName, Message, Level=3)
      END IF

      OPEN(unit = io, file = TRIM(DataF), status = 'old',iostat = ok)

      IF (ok /= 0) THEN
         WRITE(message,'(A,A)') 'Unable to open file ',TRIM(DataF)
         CALL FATAL(Trim(SolverName),Trim(message))
      END IF
            
      Nb = Nx*Ny 
          
      ALLOCATE(xb(Nb), yb(Nb), zb(Nb), xbaux(Nb), ybaux(Nb), zbaux(Nb))

      ! read datas
      DO i = 1, Nb 
         READ(io,*,iostat = ok, end=100) xbaux(i), ybaux(i), zbaux(i)
      END DO
100   Nbaux = Nb - i
      IF (Nbaux > 0) THEN
         WRITE(message,'(I0,A,I0,A,A)') Nbaux,' out of ',Nb,' datasets in file ', TRIM(DataF)
         CALL INFO(Trim(SolverName),Trim(message))         
      END IF
      CLOSE(io)

      ! Make some verifications and - in case - manipulation 
      !on the DEM structure
      dx = lx / (Nx-1.0)
      dy = ly / (Ny-1.0)
      k = 0 
      l = 0
      IF (.NOT.InvertOrder) THEN
         DO j = 1, Ny
            y = y0 + dy*(j-1)
            DO i = 1, Nx 
               k = k + 1
               x = x0 + dx*(i-1)
               IF (.NOT.FillIn) THEN
                  xb(k) = xbaux(k)
                  yb(k) = ybaux(k)
                  zb(k) = zbaux(k)
                  IF ((ABS(x-xbaux(k))>posTol*dx).OR.(ABS(y-ybaux(k))>posTol*dy)) THEN
                     
                     WRITE(Message,'(A,A)')'Structure of the DEM is not conforming to what is given in the sif for ',TRIM(FName) 
                     CALL INFO(SolverName, Message, Level=1)
                     WRITE(Message,'(A,i4,A,i4,A,e14.8,2x,e14.8,A,e14.8,2x,e14.8,A)') &
                          'Variable', NoVar, ': Found that point ',k,&
                          ' coordinate is (',xbaux(k),ybaux(k),&
                          '), whereas it should be (',x,y,')' 
                     CALL FATAL(SolverName, Message)                      
                  END IF
               ELSE
                  IF ((ABS(x-xbaux(l+1))>posTol*dx).OR.(ABS(y-ybaux(l+1))>posTol*dy)) THEN
                     xb(k) = x
                     yb(k) = y
                     zb(k) = noDataVal ! setting to NaN                   
                  ELSE
                     l=l+1
                     xb(k) = xbaux(l)
                     yb(k) = ybaux(l)
                     zb(k) = zbaux(l)
                  END IF
               END IF
            END DO
         END DO
      ELSE ! inverse order
         DO i = 1, Nx 
            x = x0 + dx*(i-1) 
            DO j = 1, Ny
               k = k + 1
               y = y0 + dy*(j-1)

               IF (.NOT.FillIn) THEN
                  xb((j-1)*Nx + i) = xbaux(k)
                  yb((j-1)*Nx + i) = ybaux(k)
                  zb((j-1)*Nx + i) = zbaux(k)
                  IF ((ABS(x-xb((j-1)*Nx + i))>posTol*dx).OR.(ABS(y-yb((j-1)*Nx + i))>posTol*dy)) THEN
                     
                     WRITE(Message,'(A,A)')'Structure of the DEM is not conforming to what is given in the sif for ',TRIM(FName) 
                     CALL INFO(SolverName, Message, Level=1)
                     WRITE(Message,'(A,i4,A,i4,A,e14.8,2x,e14.8,A,e14.8,2x,e14.8,A)') &
                          'Variable', NoVar, ':Found that point ',k,&
                          ' coordinate is (',xb((j-1)*Nx),yb((j-1)*Nx + i),'),&
                          whereas 3 it should be (',x,y,')' 
                     CALL FATAL(SolverName, Message)                      
                  END IF
               ELSE
                  IF ((ABS(x-xbaux(l+1))>posTol*dx).OR.(ABS(y-ybaux(l+1))>posTol*dy)) THEN
                     xb((j-1)*Nx + i) = x
                     yb((j-1)*Nx + i) = y
                     zb((j-1)*Nx + i) = noDataVal ! setting to NaN
                  ELSE
                     l=l+1
                     xb((j-1)*Nx + i) = xbaux(l)
                     yb((j-1)*Nx + i) = ybaux(l)
                     zb((j-1)*Nx + i) = zbaux(l)
                  END IF
               END IF

            END DO
         END DO
      END IF

      OutNode = 0
      Rmax = 0.0
      DO i=1,Model % Mesh % NumberOfNodes
         x = Model % Mesh % Nodes % x(i)
         y = Model % Mesh % Nodes % y(i)
         Rmin = 0.0
         CALL InterpolateDEM(x,y,xb,yb,zb,Nx,Ny,x0,y0,lx,ly,Rmin,z,noDataVal,noDataTol)
         if ( perm(i) .eq. 0 ) CYCLE
         Values(Perm(i)) = z
         IF (Rmin > 0.0) THEN
            OutNode = OutNode + 1
            IF (Rmin > Rmax) Rmax = Rmin
         END IF
      END DO
          
      ! Give information on the number of Nodes which are outside of the
      ! DEM domain
      IF (OutNode > 0) THEN
         WRITE( Message, '(I0,A,A)' )OutNode,' nodes where found outside of &
                 the DEM domain in ',TRIM(DataF)
         CALL Info( TRIM(SolverName), Message, Level=3 )
         WRITE( Message, '(A,e14.8)' )'The farthest DEM point used to evaluate & 
                 the nodal value was: ', Rmax
         CALL Info( TRIM(SolverName), Message, Level=3 )
      END IF
            
      DEALLOCATE(xb, yb, zb, xbaux, ybaux, zbaux)
   END DO

   CALL INFO(Trim(SolverName), '----------ALL DONE----------',Level=5)

END SUBROUTINE Grid2DInterpolator


!!!!!!!!!!!!!!!!!!!
! Subroutine InterpolateDEM
!!------------------------------------------------------------------------------!!
SUBROUTINE InterpolateDEM (x, y, xb, yb, zb, Nbx, Nby, xb0, yb0, lbx, lby, Rmin, zbed, noDataVal, noDataTol)
  USE DefUtils
  IMPLICIT NONE
  REAL(KIND=dp),INTENT(IN) :: noDataVal, noDataTol
  INTEGER :: imin, Npt, t
  INTEGER :: NMAX, i, j, Nb, Nbx, Nby, ib, ix, iy
  REAL(KIND=dp) :: x, y, zbed, xb0, yb0, x1, x2, y1, y2, zi(2,2) 
  REAL(KIND=dp) :: R, Rmin, lbx, lby, dbx, dby
  REAL(KIND=dp) :: xb(Nbx*Nby), yb(Nbx*Nby), zb(Nbx*Nby)       

  ! Find zbed for that point from the Bedrock MNT 
  dbx = lbx / (Nbx-1.0)
  dby = lby / (Nby-1.0)
  Nb = Nbx*Nby

  ix = INT((x-xb0)/dbx)+1
  iy = INT((y-yb0)/dby)+1
  ib = Nbx * (iy - 1) + ix

  ! if we are already at the end of the domain then collapse the 2 by 2 interpolation 
  ! square to just 2 points at the end of the domain (else we get interpolation involving 
  ! points at the beginning of the domain).  This comment refers to the x direction.
  IF (MOD(ib,Nbx) .eq. 0.0) THEN
     zi(2,1) = noDataVal
     zi(2,2) = noDataVal
  ELSE
     IF ( (ib+1).gt.size(zb) ) THEN
        zi(2,1) = noDataVal
     ELSE
        zi(2,1) = zb(ib+1)
     END IF
     IF ( (ib+Nbx+1).gt.size(zb) ) THEN 
        zi(2,2) = noDataVal
     ELSE
        zi(2,2) = zb(ib + Nbx + 1)
     END IF
  END IF

  x1 = xb(ib)
  IF ( (ib+1).gt.size(xb) ) THEN 
     x2 = noDataVal
  ELSE
     x2 = xb(ib+1)
  END IF

  y1 = yb(ib)
  IF ( (ib+Nbx).gt.size(yb) ) THEN 
     y2 = noDataVal
  ELSE
     y2 = yb(ib + Nbx)
  END IF
  
  IF ( (ib).gt.size(zb) ) THEN
     zi(1,1) = noDataVal
  ELSE  
     zi(1,1) = zb(ib)
  END IF
  IF ( (ib+Nbx).gt.size(zb) ) THEN
     zi(1,2) = noDataVal
  ELSE
     zi(1,2) = zb(ib + Nbx)
  END IF

  IF ( (isNoData(zi(1,1))).OR. &
       (isNoData(zi(1,2))).OR. &
       (isNoData(zi(2,1))).OR. &
       (isNoData(zi(2,2))) ) THEN

     IF ( (isNoData(zi(1,1))).AND. &
          (isNoData(zi(1,2))).AND. &
          (isNoData(zi(2,1))).AND. &
          (isNoData(zi(2,2))) ) THEN

        ! Find the nearest point available if all neighbouring points have noData
        Rmin = 9999999.0
        DO i=1, Nb
           IF (.NOT.isNoData(zb(i))) THEN
              R = SQRT((x-xb(i))**2.0+(y-yb(i))**2.0)
              IF (R<Rmin) THEN
                 Rmin = R
                 imin = i
              END IF
           END IF
        END DO
        zbed = zb(imin)

     ELSE
        ! Mean value over the available data if only some points have noData
        zbed = 0.0
        Npt = 0
        DO i=1, 2
           DO J=1, 2
              IF (.NOT. isNoData(zi(i,j))) THEN 
                 zbed = zbed + zi(i,j)
                 Npt = Npt + 1
              END IF
           END DO
        END DO
        zbed = zbed / Npt
     END IF
  ELSE
     ! linear interpolation is only carried out if all 4 neighbouring points have data.
     zbed = (zi(1,1)*(x2-x)*(y2-y)+zi(2,1)*(x-x1)*(y2-y)+zi(1,2)*(x2-x)*(y-y1)+zi(2,2)*(x-x1)*(y-y1))/(dbx*dby)      
  END IF


CONTAINS

  LOGICAL FUNCTION isNoData(val)

    IMPLICIT NONE
    REAL(KIND=dp),INTENT(IN) :: val

    IF ((val .GT. noDataVal-noDataTol) .AND. (val .LT. noDataVal+noDataTol)) THEN
       isNoData = .TRUE.
    ELSE
       isNoData = .FALSE.
    END IF

    RETURN 

  END FUNCTION isNoData

END SUBROUTINE InterpolateDEM

