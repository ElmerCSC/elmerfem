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
! *  Authors: Joe Todd
! *  Email:
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 18/10/19
! *
! *****************************************************************************
!>  Computes anisotropic target element size based on distance from calving front
SUBROUTINE GlacierMeshMetricAniso(Model, nodenumber, y, TargetLength)

  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(KIND=dp) :: y, TargetLength(:)
  !----------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Variable_t), POINTER :: TimeVar
  TYPE(Solver_t), POINTER :: Solver
  REAL(KIND=dp) :: xx,yy,zz,t,told, this_dist, MinDist2, Dist,&
       lc_mindist, lc_maxdist, lc_min, lc_max,s,dx,dz
  INTEGER, POINTER :: FrontPerm(:)
  INTEGER :: i,NNodes, NFront
  LOGICAL :: NewTime,FirstTime=.TRUE., Debug
  CHARACTER(LEN=MAX_NAME_LEN) :: FrontMaskName="Calving Front Mask"

  SAVE :: FirstTime, told, FrontPerm, Mesh, NNodes, nfront, lc_maxdist, lc_mindist,&
       lc_max, lc_min, dz, newtime

  Debug = .FALSE.
  Timevar => VariableGet( Model % Variables,'Time')
  t = TimeVar % Values(1)
  Solver => Model % Solver

  IF (FirstTime) THEN
    FirstTime = .FALSE.
    NewTime = .TRUE.
    told = t

    Mesh => Model % Mesh
    Material => GetMaterial(Mesh % Elements(1)) !TODO, this is not generalised

    lc_maxdist = ListGetConstReal(Material, "GlacierMeshMetric Max Distance",  Default=2000.0_dp)
    lc_mindist = ListGetConstReal(Material, "GlacierMeshMetric Min Distance",  Default=200.0_dp)
    lc_max = ListGetConstReal(Material, "GlacierMeshMetric Max LC",  Default=500.0_dp)
    lc_min = ListGetConstReal(Material, "GlacierMeshMetric Min LC",  Default=75.0_dp)
    dz = ListGetConstReal(Material, "GlacierMeshMetric Vertical LC",  Default=50.0_dp)
  END IF

  IF(t > told) NewTime = .TRUE. !TODO - replace this with Samuel's mesh counter logic
  IF(NewTime) THEN
    told = t
    NewTime = .FALSE.

    Mesh => Model % Mesh
    NNodes = Mesh % NumberOfNodes

    !mask is the most computationally expense element of this routine
    IF(ASSOCIATED(FrontPerm)) DEALLOCATE(FrontPerm)
    CALL MakePermUsingMask( Model, Solver, Mesh, FrontMaskName, &
         .FALSE., FrontPerm, nfront, ParallelComm = .FALSE.)

  END IF

  xx = Mesh % Nodes % x(nodenumber)
  yy = Mesh % Nodes % y(nodenumber)
  zz = Mesh % Nodes % z(nodenumber)

  MinDist2 = HUGE(MinDist2)
  DO i=1,NNodes
    IF(FrontPerm(i) == 0) CYCLE
    !Note - sqrt is an expensive FLOP so only do it once...
    this_dist = ((xx - Mesh % Nodes % x(i))**2.0) + &
         ((yy - Mesh % Nodes % y(i))**2.0) + &
         ((zz - Mesh % Nodes % z(i))**2.0)
    MinDist2 = MIN(this_dist, MinDist2)
  END DO
  Dist = SQRT(MinDist2)

  !Apply caps to this distance:
  Dist = MAX(Dist,lc_mindist)
  Dist = MIN(Dist,lc_maxdist)

  s = (Dist - lc_mindist) / (lc_maxdist - lc_mindist)
  dx = s*lc_max + (1-s)*lc_min

  TargetLength(1) = dx
  TargetLength(2) = dx
  TargetLength(3) = dz

  IF(Debug) PRINT *,'Node ',nodenumber,'TargetLength: ',TargetLength,' dist: ',Dist,' s: ',s

END SUBROUTINE GlacierMeshMetricAniso
