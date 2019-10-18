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
! *  Authors: F. Gillet-Chaulet (IGE-France)
! *  Web:     http://elmerice.elmerfem.org
! *  Original Date: 04/2019
! * 
! *****************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write the grid description file for CDO (https://code.mpimet.mpg.de/projects/cdo)
!   cf e.g. http://www.climate-cryosphere.org/wiki/index.php?title=Regridding_with_CDO
!   Limitations: restricted to meshes with 303 elements
!   Requires: fortrangis (http://fortrangis.sourceforge.net) with proj support
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE MakeCDOGridDesc( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      USE fortranc
      USE proj
      IMPLICIT NONE

      TYPE(Model_t) :: Model
      TYPE(Solver_t):: Solver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

      TYPE(ValueList_t), POINTER :: SolverParams
      INTEGER :: t,n,ne,i
      TYPE(Element_t),POINTER :: Element
      INTEGER, POINTER :: NodeIndexes(:)
      REAL(KIND=dp),ALLOCATABLE :: LatLon(:,:),Corner(:,:,:)
      REAL(KIND=dp) :: x(3),y(3),xg,yg,xt(3),yt(3)
      INTEGER, parameter :: io=12
      INTEGER :: nclock=0
      CHARACTER(len=100) :: Filename

      TYPE(pjuv_object) :: coordp,coordg
      TYPE(pj_object) :: pj
      CHARACTER(LEN=MAX_NAME_LEN) :: proj_string
    
! Get projection definition
      SolverParams => GetSolverParams()
      proj_string = ListGetString(SolverParams,'projection',UnFoundFatal=.TRUE.)
      pj = pj_init_plus(TRIM(proj_string)//CHAR(0))

      ne=GetNOFActive()
      allocate(LatLon(ne,2),Corner(ne,3,2))

      DO t = 1,ne
         Element => GetActiveElement(t)
         NodeIndexes => Element % NodeIndexes

         IF (Element % TYPE % ElementCode.NE.303) &
            CALL FATAL('MakeCDOGrid','Wrong element type')

         n = GetElementNOFNodes()
         ! element center
         xg=SUM(Solver%Mesh%Nodes%x(NodeIndexes(1:n)))/n
         yg=SUM(Solver%Mesh%Nodes%y(NodeIndexes(1:n)))/n

         ! convert to lon,lat
         coordp = pjuv_object(xg,yg)
         coordg = pj_inv(coordp, pj)
         LatLon(t,1)=coordg % v * pj_rad_to_deg
         LatLon(t,2)=coordg % u * pj_rad_to_deg

         ! element corners
         x(1:n)=Solver%Mesh%Nodes%x(NodeIndexes(1:n))
         y(1:n)=Solver%Mesh%Nodes%y(NodeIndexes(1:n))
         ! CDO requires anti-clokwise ordering
         IF(IsClockwise(x,y,n)) THEN
           nclock=nclock+1
           xt=x
           yt=y
           Do i=1,3
             x(4-i)=xt(i)
             y(4-i)=yt(i)
           End do
         END IF
         !convert corners to lon,lat
         Do i=1,n
          coordp = pjuv_object(x(i),y(i))
          coordg = pj_inv(coordp, pj)
          Corner(t,i,1)=coordg % v * pj_rad_to_deg
          Corner(t,i,2)=coordg % u * pj_rad_to_deg
         End do

      End do
      
      PRINT *,Model%Mesh%name 
      ! write the grid description file
      IF (ParEnv % Pes>1 ) THEN
         write(FileName,'(a,a,i0,a)') TRIM(Model%Mesh%name),'_CDOGrid_',ParEnv%MyPe,'.txt'
      ELSE
         write(FileName,'(a,a)') TRIM(Model%Mesh%name),'_CDOGrid.txt'
      END IF

      open(io,file=trim(Filename))
      write(io,*) 'gridtype = unstructured'
      write(io,*) 'gridsize = ',ne
      write(io,*) 'nvertex = ',3

      write(io,*) 'xvals = ',LatLon(1,2)
      Do t=2,ne
         write(io,*) LatLon(t,2)
      End do
      write(io,*) 'xbounds = ',Corner(1,1:3,2)
      Do t=2,ne
         write(io,*) Corner(t,1:3,2)
      End do

      write(io,*) 'yvals = ',LatLon(1,1)
      Do t=2,ne
         write(io,*) LatLon(t,1)
      End do
      write(io,*) 'ybounds = ',Corner(1,1:3,1)
      Do t=2,ne
         write(io,*) Corner(t,1:3,1)
      End do


      close(io)

      deallocate(LatLon,Corner)

      CONTAINS 
      ! Are nodes given clockwise?
      FUNCTION IsClockwise(x,y,n) RESULT(clokwise)
      LOGICAL :: clokwise
      REAL(KIND=dp) :: x(n),y(n)
      INTEGER :: n
      INTEGER :: i,ind2
      REAL(KIND=dp) :: sarea

      sarea=0._dp
      Do i=1,n-1
        ind2=mod(i+1,n)
        sarea=sarea+(x(ind2)-x(i))*(y(ind2)+y(i))
      End do
      clokwise=(sarea.GT.0)

      END FUNCTION IsClockwise 
      End
