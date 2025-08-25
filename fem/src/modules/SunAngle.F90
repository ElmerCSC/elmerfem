!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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
!
! ******************************************************************************
! *
! *  Module for computing the critical sun angle. Not yet parallel!
! *  Also very stupid search algo but seems accurate.
! *
! *  Authors: Peter Råback
! *  Email:   peter.raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 18.07.2025
! *
! *****************************************************************************/


SUBROUTINE SunAngleSolver_init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
    
  CALL ListAddNewString( Solver % Values,'Variable','CriticalTan')
  CALL ListAddLogical( Solver % Values,'No Matrix',.TRUE.)
  
END SUBROUTINE SunAngleSolver_init


!------------------------------------------------------------------------------
!> Subroutine for extracting isosurfaces in 3d and isolines in 2d.
!> The result will be a new mesh which will be added to the list of meshes.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE SunAngleSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver  !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model    !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt       !< Timestep size for time dependent simulations
  LOGICAL :: Transient      !< Steady state or transient simulation
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER :: i,j,k,n,t,dofs,dim,idof
  INTEGER, POINTER :: SunPerm(:)
  INTEGER :: Inds(3)
  TYPE(Variable_t), POINTER :: SunVar
  LOGICAL :: Found, Singular
  INTEGER :: ElemFirst,ElemLast
  REAL(KIND=dp), POINTER :: SunAngle(:)
  REAL(KIND=dp) :: tanphi, ds, dz, Amat(2,2), b(2), c(2), x0(3), y0(3), z0(3), &
      Alpha, Eps, ElemHeight(4)
  TYPE(Element_t), POINTER  :: Element
  TYPE(ValueList_t), POINTER :: Params, Material
  REAL(KIND=dp) :: detA
  CHARACTER(*), PARAMETER :: Caller = 'SunAngleSolver'

  
  CALL Info( Caller,'-------------------------------------',Level=4 )
  CALL Info( Caller,'Determining the critical sun angle',Level=4 )
  CALL Info( Caller,'-------------------------------------',Level=4 )

  IF(ParEnv % PEs > 1) THEN
    CALL Fatal(Caller,'Routine has not been implemented in parallel!')
  END IF
  
  Mesh => GetMesh()
  Params => GetSolverParams()

  SunVar => Solver % Variable
  SunAngle => SunVar % Values
  SunPerm => SunVar % Perm
  dofs = SunVar % dofs
  
  ! Find mesh edges in order to define the intersection points
  !-----------------------------------------------------------
  IF (.NOT.ASSOCIATED(Mesh % Edges)) THEN
    CALL Info(Caller,'Creating mesh edges',Level=10)
    IF( dim == 2 ) THEN
      CALL FindMeshEdges2D(Mesh)
    ELSE
      CALL FindMeshEdges3D(Mesh)
    END IF
  END IF


  ! For testing purposes set the height here.
  !---------------------------------------------------------------
  IF( ListCheckPresentAnyMaterial(Model,'Test Height') ) THEN
    IF( Mesh % Elements(1) % TYPE % ElementCode > 500 ) THEN
      ElemFirst = Mesh % NumberOfBulkElements + 1
      ElemLast = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
    ELSE
      ElemFirst = 1
      ElemLast = Mesh % NumberOfBulkElements
    END IF
    
    DO t = ElemFirst,ElemLast 
      Element => Mesh % Elements( t ) 
      Model % CurrentElement => Element

      Material => GetMaterial(Element)
      n = Element % TYPE % NumberOfNodes
      ElemHeight(1:n) = GetReal(Material,'Test Height',Found) 
      
      IF(ANY(SunPerm(Element % NodeIndexes)==0)) CYCLE      
      Mesh % Nodes % z(Element % NodeIndexes) = ElemHeight(1:n)      
    END DO
  END IF
    

  eps = 1.0e-6
  SunAngle = -100.0_dp

  DO idof=1,dofs
    IF(dofs==1) THEN
      Alpha = ListGetCReal( Params,'Test Angle')
    ELSE
      Alpha = (idof-1)*2.0_dp*PI/dofs
    END IF

    Amat(1,1) = COS(Alpha)
    Amat(2,1) = SIN(Alpha)

    ! This is a stupid N^2 algo where we have nested loops over edges and nodes.
    DO t = 1, Mesh % NumberOfEdges
      Element => Mesh % Edges( t )
      Inds(1:2) = Element % NodeIndexes
      IF(ANY(SunPerm(Inds(1:2))==0)) CYCLE

      x0(1:2) = Mesh % Nodes % x(Inds(1:2))
      y0(1:2) = Mesh % Nodes % y(Inds(1:2))
      
      Amat(1,2) = x0(1)-x0(2)
      Amat(2,2) = y0(1)-y0(2)
      
      ! Note: for triangles and quads "n" is both the number of nodes and
      ! number of edges!      

      ! node which might be shaded by the edge.
      DO j=1, Mesh % NumberOfNodes
        k = SunPerm(j)
        IF(k==0) CYCLE

        ! We don't want to study nodes that are part of the edge.
        IF(ANY(Inds(1:2)==j)) CYCLE
        
        Inds(3) = j
        x0(3) = Mesh % Nodes % x(j)
        y0(3) = Mesh % Nodes % y(j)
        
        b(1) = x0(1)-x0(3)
        b(2) = y0(1)-y0(3)
        
        CALL Solve2x2(Amat,c,b,Singular)
        IF(Singular) CYCLE
        
        IF(c(1) < Eps ) CYCLE          
        IF(c(2) < -Eps .OR. c(2) > 1.0_dp+EPs ) CYCLE
        
        z0 = Mesh % Nodes % z(Inds)

        ds = c(1)
        dz = ( z0(1) + c(2) * (z0(2)-z0(1)) ) - z0(3)

        tanphi = dz/ds

        IF(dofs==1) THEN
          SunAngle(k) = MAX(SunAngle(k),tanphi)
        ELSE
          SunAngle(dofs*(k-1)+idof) = MAX(SunAngle(dofs*(k-1)+idof),tanphi)
        END IF
      END DO
    END DO
  END DO

  ! Set default value to unset nodes (=boundary at sunside)
  WHERE(SunAngle < -99.0 )
    SunAngle = 0.0_dp
  END WHERE


CONTAINS

!------------------------------------------------------------------------------
!> Solves a 2x2 linear system indicating if the system is singular.
!------------------------------------------------------------------------------
   SUBROUTINE Solve2x2( A, x, b, Singular )
!------------------------------------------------------------------------------
     REAL(KIND=dp), INTENT(out) :: x(:)
     REAL(KIND=dp), INTENT(in)  :: A(:,:),b(:)
     LOGICAL, INTENT(out) :: Singular
!------------------------------------------------------------------------------
!     REAL(KIND=dp) :: detA
!------------------------------------------------------------------------------
     detA = A(1,1) * A(2,2) - A(1,2) * A(2,1)

     Singular = ( ABS(detA) < EPSILON(detA) )
     IF(Singular) RETURN
     
     detA = 1.0d0 / detA
     x(1) = detA * (A(2,2) * b(1) - A(1,2) * b(2))
     x(2) = detA * (A(1,1) * b(2) - A(2,1) * b(1))
!------------------------------------------------------------------------------
   END SUBROUTINE Solve2x2
!------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
END SUBROUTINE SunAngleSolver
!------------------------------------------------------------------------------


