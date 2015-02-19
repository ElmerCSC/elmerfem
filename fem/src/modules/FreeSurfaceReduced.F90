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
!/******************************************************************************
! *
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 2001
! *
! ****************************************************************************/
   
!------------------------------------------------------------------------------
!> Solve the equation resulting from 1D drawing process. 
!> May be applied to 2D and axisymetric cases. Possible uses include drawing of 
!> viscous fibers or sheets.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE FreeSurfaceReduced( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
     USE DefUtils

     IMPLICIT NONE
!------------------------------------------------------------------------------
 
     TYPE(Model_t) :: Model
     TYPE(Solver_t), TARGET:: Solver
     REAL (KIND=DP) :: dt
     LOGICAL :: TransientSimulation
 
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Variable_t), POINTER :: Var, MoveCoord
     TYPE(Element_t),POINTER :: CurrentElement
     
     REAL (KIND=DP) :: Norm,Relax,MaxRad
     REAL (KIND=DP) :: xnew, xold, dx, fsum, fsum0, df 
     REAL (KIND=DP) :: a, b, v1, v2, x1, x2, y1, y2, v1y, v2y, v1x, v2x, &
         negerror, poserror, maxerror, x0
     REAL (KIND=DP) :: thickmin, thickmax, rnorm(3)
     REAL (KIND=DP), POINTER :: coords(:), coords2(:), OrigCoords(:)

     INTEGER, POINTER :: NodeIndexes(:), LeftNeighbours(:), RightNeighbours(:), &
         FreeSurfacePoints(:), UpwindPoints(:)
     INTEGER :: i,j,k,n,t,istart,iend,istat, DrawDirection, &
         DonePoints, CoordSystem, SubroutineVisited=0, NoFreeSurfaces, &
         Surface, tfirst

     LOGICAL :: AllocationsDone = .FALSE., axisymmetric, gotIt, IsFreeSurface, &
         PerformMapping

     SAVE FreeSurfacePoints, AllocationsDone, SubroutineVisited, &
         LeftNeighbours, RightNeighbours, DrawDirection, &
         coords, coords2, OrigCoords, NoFreeSurfaces, UpwindPoints
    
     CALL Info( 'FreeSurfaceReduced', '-------------------------------------',Level=4 )
     CALL Info( 'FreeSurfaceReduced', '2D Free Surface Solver:  ', Level=4 )
     CALL Info( 'FreeSurfaceReduced', '-------------------------------------',Level=4 )

     Mesh => Solver % Mesh
     MoveCoord => Solver % Variable 
     Var => VariableGet( Mesh % Variables, 'Flow Solution', .TRUE. )
     
     PerformMapping = ListGetLogical( Solver % Values, &
         'Perform Mapping',gotIt)
     IF(.NOT. gotIt) PerformMapping = .TRUE.
 
!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN
       N = Model % MaxElementNodes
       
       ALLOCATE( FreeSurfacePoints(Model%NumberOfNodes), &
           LeftNeighbours(Model%NumberOfNodes), &
           RightNeighbours(Model%NumberOfNodes), &
           OrigCoords(Model%NumberOfNodes), &
           STAT=istat )
       IF ( istat /= 0 ) CALL Warn('FreeSurfaceReduced','Memory allocation error')
       
       AllocationsDone = .TRUE.
     END IF

     CoordSystem = CurrentCoordinateSystem()
     IF(CoordSystem == CylindricSymmetric .OR. CoordSystem == AxisSymmetric) THEN
       axisymmetric = .TRUE.
     ELSE
       axisymmetric = .FALSE.
     ENDIF

     Relax = GetCReal( Solver % Values, &
         'Nonlinear System Relaxation Factor',gotIt)
     IF(.NOT. gotIt) Relax = 1.0

     ! Mark the nodes that are at the free surface
     IF(SubroutineVisited == 0) THEN
       
       FreeSurfacePoints = 0
       NoFreeSurfaces = 0
       
       DO t = Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + &
           Mesh % NumberOfBoundaryElements
         
         CurrentElement => Mesh % Elements(t)
         IF ( CurrentElement % TYPE % ElementCode == 101 ) CYCLE
         
         Model % CurrentElement => CurrentElement
         !------------------------------------------------------------------------------
         n = CurrentElement % TYPE % NumberOfNodes
         NodeIndexes => CurrentElement % NodeIndexes
         
         DO k=1, Model % NumberOfBCs
           IF ( Model % BCs(k) % Tag /= CurrentElement % BoundaryInfo % Constraint ) CYCLE
           
           IsFreeSurface = ListGetLogical(Model % BCs(k) % Values,'Free Surface Reduced',gotIt ) 
           
           IF(gotIt .AND. IsFreeSurface) THEN         
             Surface = ListGetInteger(Model % BCs(k) % Values,'Free Surface Number',GotIt)
             IF(.NOT. GotIt) Surface = 1
             FreeSurfacePoints(NodeIndexes(1:n)) = Surface
             NoFreeSurfaces = MAX(NoFreeSurfaces,Surface)
           END IF

           IsFreeSurface = ListGetLogical(Model % BCs(k) % Values,'Free Surface Bottom',gotIt ) 
           IF(gotIt .AND. IsFreeSurface) FreeSurfacePoints(NodeIndexes(1:n)) = -1

         END DO
         
       END DO

       WRITE(Message,'(A,I1)') 'Separate Free surfaces found ',NoFreeSurfaces
       CALL Info('FreeSurfaceReduced',Message,Level=5)

       ALLOCATE(UpwindPoints(NoFreeSurfaces))
       UpwindPoints = 0

       ! use simple decision rules to decide in which direction the free surface goes
       !-----------------------------------------------------------------------------
       GotIt = .FALSE.
       DO t=1,Model%NumberOfNodes
         IF(FreeSurfacePoints(t) <= 0) CYCLE
         IF(.NOT. GotIt) THEN
           GotIt = .TRUE.
           x1 = Mesh % Nodes % x(t)
           x2 = x1
           y1 = Mesh % Nodes % y(t)
           y2 = y1
         ELSE          
           x1 = MAX(x1, Mesh % Nodes % x(t))
           x2 = MIN(x2, Mesh % Nodes % x(t))
           y1 = MAX(y1, Mesh % Nodes % y(t))
           y2 = MIN(y2, Mesh % Nodes % y(t))         
         END IF
       END DO

       ! in x or y direction
       IF(x1-x2 > y1-y2) THEN
         DrawDirection = 1
       ELSE
         DrawDirection = 2
       END IF

       GotIt = .FALSE.
       DO t=1,Model%NumberOfNodes
         IF(FreeSurfacePoints(t) <= 0) CYCLE
         IF(.NOT. GotIt) THEN
           GotIt = .TRUE.
           v1 = Var % Values(Var % DOFs * (Var % Perm(t)-1)+DrawDirection)
           v2 = v1
         ELSE
           v1 = MAX(v1, Var % Values(Var % DOFs * (Var % Perm(t)-1)+DrawDirection))
           v2 = MIN(v2, Var % Values(Var % DOFs * (Var % Perm(t)-1)+DrawDirection))
         END IF
       END DO

       ! In positive or negative direction
       IF(-v2 > v1) THEN
         DrawDirection = -DrawDirection
       END IF


       WRITE(Message,'(A,I3)') 'DrawDirection was found to be',DrawDirection
       CALL Info('FreeSurfaceReduced',Message,Level=5)

       rnorm = 0.0d0
       IF(ABS(DrawDirection) == 2) THEN
         rnorm(1) = -1.0 
       ELSE
         rnorm(2) = -1.0
       END IF
 
       CALL FindNeighbourNodes( Mesh,rnorm,LeftNeighbours)
       CALL FindNeighbourNodes( Mesh,-1.0*rnorm,RightNeighbours)

       IF(ABS(DrawDirection) == 2) THEN
         coords => Mesh % Nodes % x
         coords2 => Mesh % Nodes % y
         OrigCoords = Mesh % Nodes % x
       ELSE 
         coords => Mesh % Nodes % y
         coords2 => Mesh % Nodes % x
         OrigCoords = Mesh % Nodes % y
       END IF
       
       ! Find the first nodes on the upwind side
       DO Surface = 1,NoFreeSurfaces

         gotIt = .FALSE.
      
         DO t=1,Model%NumberOfNodes

           IF(FreeSurfacePoints(t) == Surface) THEN
             IF(ABS(DrawDirection) == 2) THEN
               x2 = Mesh % Nodes % y(t)
             ELSE
               x2 = Mesh % Nodes % x(t)
             END IF
             IF(.NOT. GotIt) THEN
               x1 = x2
               GotIt = .TRUE.
             END IF

             IF(x2 * DrawDirection < x1 * DrawDirection) THEN
               x1 = x2
               UpwindPoints(Surface) = t
             END IF
           END IF
         END DO

         WRITE(Message,'(A,I6)') 'First upwinds node for surface is',UpwindPoints(Surface)
         CALL Info('FreeSurfaceReduced',Message,Level=5)

       END DO

       ! The initial displacement shoud be zero
       MoveCoord % Values = 0.0d0
     END IF


     DO Surface = 1,NoFreeSurfaces

       DonePoints = 0
       negerror = 0.0d0
       poserror = 0.0d0
       
       thickmin = HUGE(thickmin)
       thickmax = -HUGE(thickmax)
       tfirst = UpwindPoints(Surface)

       ! Go through all the free surface points
       !------------------------------------------
	   DO t=tfirst,Model%NumberOfNodes+tfirst-1
         
         iend = t
         IF(iend > Model%NumberOfNodes) iend = iend - Model%NumberOfNodes
         
         IF( FreeSurfacePoints(iend) /= Surface) CYCLE
       
         xold = coords(iend)
         x1 = coords(iend)

         ! Find the path from the free surface to the symmetry boundary
         j = iend
         DO 
           j = LeftNeighbours(j)
           
           IF(j==0) EXIT
           
           x2 = x1
           x1 = coords(j)
           
           IF(x1 > x2) THEN
             EXIT
           END IF
           
           IF(MoveCoord % Perm(j) <= 0) THEN
             CALL Warn('FreeSurfaceReduced','Problems as Perm(i) < 0')
             EXIT
           END IF
         
           istart = j

           ! Test if the next free, or the free surface bottom, has been reached.
           IF(FreeSurfacePoints(j) /= 0) THEN
             EXIT
           END IF

         END DO

         IF(iend < 1 .OR. iend > Model%NumberOfNodes .OR. &
             istart < 1 .OR. istart > Model%NumberOfNodes) THEN 
           CALL Warn('FreeSurfaceReduced','Invalid integration limits')
         END IF

         k = Var % PERM(istart)
         IF(k == 0) CALL Warn('FreeSurfaceReduced','no flow solution for given point')

         v2x = Var % Values(Var % DOFs * (k-1)+ABS(DrawDirection))
         v2y = Var % Values(Var % DOFs * (k-1)+3-ABS(DrawDirection))
         x2 = coords(istart)
         y2 = coords2(istart)
         fsum = 0.0d0

         ! Integrate from the fixed boundary until the desired flux is exceeded
         !----------------------------------------------------------------------
         i = istart
         DO 
           v1x = v2x
           v1y = v2y
           x1 = x2
           y1 = y2
           j = RightNeighbours(i)
           k = Var % Perm(j)
           IF(k==0) CALL Warn('FreeSurfaceReduced','no flow solution for given point')
           
           v2x = Var % Values(Var % DOFs * (k-1)+ABS(DrawDirection))
           v2y = Var % Values(Var % DOFs * (k-1)+3-ABS(DrawDirection))
           x2 = coords(j)
           y2 = coords2(j)
           
         
           ! Correct the speed in case the elements are not aligned so that dy is zero
           !---------------------------------------------------------------------------
           IF(ABS(x2-x1) < 1.0d-10) CALL Warn('FreeSurfaceReduced','Integration path has dx=0')
           
           ! Correct speeds by taking in account the angle of integration
           v1 = v1x - v1y * (y2-y1) / (x2-x1)
           v2 = v2x - v2y * (y2-y1) / (x2-x1)
           
           a = (v1*x2-v2*x1)/(x2-x1)
           b = (v2-v1)/(x2-x1)
           
           ! These fluxes are based on analytical integration over linear elements
           !----------------------------------------------------------------------
           IF(axisymmetric) THEN
             df = a*(x2**2-x1**2)/2.0d0 + b*(x2**3-x1**3)/3.0d0
           ELSE
             df = a*(x2-x1) + b*(x2**2-x1**2)/2.0d0
           END IF
           
           df = ABS(df)
           
           fsum = fsum + df
           
           IF(j == iend) EXIT
           
           IF(DonePoints > 0 .AND. fsum > fsum0) EXIT

           i = j
         END DO

         ! For the first node set the desired flux, 
         ! For other nodes move the whole chain to the corrected positions.
         !-----------------------------------------------------------------
		 IF(DonePoints == 0) THEN
           fsum0 = fsum
           
           xnew = x2 
           x0 = coords(istart)
           IF(thickmin > xnew - x0) thickmin = xnew - x0
           IF(thickmax < xnew - x0) thickmax = xnew - x0
           
         ELSE  
           
           IF(axisymmetric) THEN
             dx = (fsum0 - fsum) / ABS(v2*x2)
           ELSE 
             dx = (fsum0 - fsum) / ABS(v2)
           END IF
           
           xnew = x2 + dx
           x0 = coords(istart)
           
           IF(dx/(xnew-x0) > poserror) poserror = dx/(xnew-x0)
           IF(dx/(xnew-x0) < negerror) negerror = dx/(xnew-x0)
           IF(thickmin > xnew - x0) thickmin = xnew - x0
           IF(thickmax < xnew - x0) thickmax = xnew - x0
           
           IF(PerformMapping) THEN
             i=istart
             DO 
               i = RightNeighbours(i)
               x1 = x0 + (xnew-x0)*(coords(i)-x0)/(xold-x0)          
               
               j = MoveCoord % Perm(i)
               MoveCoord % Values(j) = (1.0-Relax) * MoveCoord % Values(j) + &
                   Relax * (x1 - OrigCoords(i))
               
               coords(i) = OrigCoords(i) + MoveCoord % Values(j)
               
               IF(i == iend) EXIT 
             END DO
           ELSE
             j = MoveCoord % Perm(iend)
             MoveCoord % Values(j) = (1.0-Relax) * MoveCoord % Values(j) + &
                 Relax * (xnew - OrigCoords(iend))  
           END IF

         END IF

         DonePoints = DonePoints + 1
       END DO

       maxerror = MAX(-negerror,poserror)

       IF(NoFreeSurfaces > 1)  THEN
         WRITE(Message,'(A,I1)') 'Free surface number ',Surface
         CALL Info('FreeSurfaceReduced',Message,Level=4)
       END IF

       WRITE(Message,'(A,T30,ES17.6E2)') 'Free surface flux:',fsum0
       CALL Info('FreeSurfaceReduced',Message,Level=4)
       WRITE(Message,'(A,T30,ES17.6E2)') 'Free surface error:',maxerror
       CALL Info('FreeSurfaceReduced',Message,Level=4)
       WRITE(Message,'(A,T30,ES17.6E2)') 'Minimum thickness:',thickmin
       CALL Info('FreeSurfaceReduced',Message,Level=4)
       WRITE(Message,'(A,T30,ES17.6E2)') 'Maximum thickness:',thickmax
       CALL Info('FreeSurfaceReduced',Message,Level=4)

       IF(Surface == 1) THEN
         Norm = 1.0d0 + maxerror
       ELSE 
         Norm = MAX(Norm,1.0d0+maxerror)
       END IF
     END DO

     Solver % Variable % Norm = Norm

     Solver % Variable % Norm = 1.0d0 + maxerror
     CALL ListAddConstReal( Model % Simulation, 'res: Free surface flux',fsum0)
     CALL ListAddConstReal( Model % Simulation, 'res: Free surface max. error',maxerror)
     CALL ListAddConstReal( Model % Simulation, 'res: Free thickness min',thickmin)
     CALL ListAddConstReal( Model % Simulation, 'res: Free thickness max',thickmax)
     
     SubroutineVisited = SubroutineVisited + 1

END SUBROUTINE FreeSurfaceReduced


