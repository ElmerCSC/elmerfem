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
!> Subroutine for computing critical angle of sun when it is not any more seen
!> in the horizon. The angle depends on the topograghy and on the direction of
!> the sun in the sky.
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
  TYPE(Variable_t), POINTER :: SunVar
  LOGICAL :: Found, Singular, Parallel, EnforceParallel, EnforceSerial, &
      GotMaxDist, Visited = .FALSE.
  INTEGER :: ElemFirst,ElemLast
  REAL(KIND=dp), POINTER :: SunAngle(:), x(:), y(:), z(:)
  REAL(KIND=dp) :: tanphi, ds, dz, Amat(2,2), b(2), c(2), &
      x1, y1, z1, Alpha, Eps, ElemHeight(4), MaxDist, tanphi2
  TYPE(Element_t), POINTER  :: Element
  TYPE(ValueList_t), POINTER :: Params, Material
  REAL(KIND=dp) :: detA

  TYPE RidgeData_t
    INTEGER :: n_recv = 0, n_send = 0
    REAL(KIND=dp), ALLOCATABLE :: Coords(:)
    REAL(KIND=dp) :: MeshBB(4)      
  END TYPE RidgeData_t
  TYPE(RidgeData_t),  ALLOCATABLE, TARGET, SAVE :: RidgeData(:)
  
  CHARACTER(*), PARAMETER :: Caller = 'SunAngleSolver'

  
  CALL Info( Caller,'-------------------------------------',Level=4 )
  CALL Info( Caller,'Determining the critical sun angle',Level=4 )
  CALL Info( Caller,'-------------------------------------',Level=4 )
  
  Mesh => GetMesh()
  x => Mesh % Nodes % x
  y => Mesh % Nodes % y
  z => Mesh % Nodes % z
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
    CALL SetTestHeight()
  END IF

  ! Enforce use of the same algo in serial as is used in parallel. 
  EnforceParallel = ListGetLogical( Params,'Enforce Parallel', Found )

  ! Even in parallel do not do any communication. 
  EnforceSerial = ListGetLogical( Params,'Enforce Serial', Found )
  
  Parallel = EnforceParallel .OR. (ParEnv % PEs > 1)
  IF(Parallel .AND. .NOT. Visited ) THEN
    MaxDist = ListGetCReal( Params,'SunAngle search distance', GotMaxDist )
    CALL Info(Caller,'Communicating ridges from different partitions',Level=6)
    CALL PopulateRidge()
  END IF

  CALL Info(Caller,'Computing the sun angles for real!',Level=7)
  
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

    ! We use different search algo locally unless not otherwise requested.
    ! This algo uses less memory. 
    IF(.NOT. EnforceParallel) CALL SetShadingAngle()

    IF( Parallel ) CALL SetShadingAngleParallel()          
  END DO
    
  ! Set default value to unset nodes (=boundary at sunside)
  WHERE(SunAngle < -99.0 )
    SunAngle = 0.0_dp
  END WHERE

  Visited = .TRUE.

  ! This is a heuristics since when we compute many sun angles at the same time
  ! we can generally continue and not return to this routine for some time.
  IF(dofs>1) CALL ReleaseRidge()
  
  CALL Info(Caller,'All done',Level=10)  

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
   !> This is just to make it easier to test the routine as we can here define
   !> the topography locally. Generally the topography would be externall provided
   !> with the mesh. 
   !------------------------------------------------------------------------------      
   SUBROUTINE SetTestHeight()

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
     
   END SUBROUTINE SetTestHeight


   !------------------------------------------------------------------------------
   !> Compute the shading angle in serial. We have the edges available here
   !> so no need to create any additional structures. The routine has been
   !> made to include so few operations as possible but uses stupid search algo.  
   !------------------------------------------------------------------------------   
   SUBROUTINE SetShadingAngle()

     REAL(KIND=dp) :: x1,y1,z1,x0(2),y0(2),z0(2),c(2),b(2),dz,dz0,maxangle,tanphi
     INTEGER :: t,j,k,Inds(2)
     LOGICAL :: Singular
     
     ! This is a stupid N^2 algo where we have nested loops over edges and nodes.
     DO t = 1, Mesh % NumberOfEdges
       Element => Mesh % Edges( t )
       Inds(1:2) = Element % NodeIndexes
       IF(ANY(SunPerm(Inds(1:2))==0)) CYCLE

       x0(1:2) = x(Inds(1:2))
       y0(1:2) = y(Inds(1:2))
       z0(1:2) = z(Inds(1:2))
       dz0 = z0(2)-z0(1)
       
       Amat(1,2) = x0(1)-x0(2)
       Amat(2,2) = y0(1)-y0(2)

       ! Note: for triangles and quads "n" is both the number of nodes and
       ! number of edges!      

       ! node which might be shaded by the edge.
       DO j=1, Mesh % NumberOfNodes
         k = SunPerm(j)
         IF(k==0) CYCLE

         ! We don't want to study nodes that are part of the edge.
         IF(ANY(Inds==j)) CYCLE

         b(1) = x0(1) - x(j)
         b(2) = y0(1) - y(j)

         CALL Solve2x2(Amat,c,b,Singular)
         IF(Singular .OR. c(1) < Eps .OR. c(2) < -Eps .OR. c(2) > 1.0_dp+EPs ) CYCLE

         dz = ( z0(1) + c(2) * dz0 ) - z(j)
         tanphi = dz/c(1)

         SunAngle(dofs*(k-1)+idof) = MAX(SunAngle(dofs*(k-1)+idof),tanphi)
       END DO
     END DO
   END SUBROUTINE SetShadingAngle


   !------------------------------------------------------------------------------
   !> Populate ridge date to compute shading angle in parallel.
   !------------------------------------------------------------------------------
   SUBROUTINE PopulateRidge()
     !------------------------------------------------------------------------------      
     REAL(KIND=dp) :: x0,y0,x1,y1,z0,z1
     INTEGER :: i,j,k,l,n,i0,i1,n_send,n_recv,n_max
     INTEGER :: MyPe, PEs, proc
     REAL(KIND=dp), ALLOCATABLE, TARGET :: SendCoords(:)
     REAL(KIND=dp), POINTER :: pCoords(:)
     INTEGER, ALLOCATABLE :: rPar(:)
     INTEGER :: comm, ierr, status(MPI_STATUS_SIZE)
     !------------------------------------------------------------------------------

     MyPe = ParEnv % MyPe + 1
     PEs = ParEnv % PEs

     ! Allocate data 
     ALLOCATE(RidgeData(PEs))
     n = Mesh % NumberOfEdges
     
     ! Bounding box for making quick and dirty search.
     IF(GotMaxDist) THEN
       RidgeData(MyPe) % MeshBB(1) = MINVAL(x)
       RidgeData(MyPe) % MeshBB(2) = MAXVAL(x)
       RidgeData(MyPe) % MeshBB(3) = MINVAL(y)
       RidgeData(MyPe) % MeshBB(4) = MAXVAL(y)
     END IF

     ! This is done to check that the serial and parallel search are identical.
     IF( EnforceParallel ) THEN
       RidgeData(MyPe) % n_recv = n
       ALLOCATE(RidgeData(MyPe) % Coords(6*n))
       RidgeData(MyPe) % Coords = 0.0_dp
       
       pCoords => RidgeData(MyPe) % Coords
       DO k=1,Mesh % NumberOfEdges 
         i0 = Mesh % Edges(k) % NodeIndexes(1)
         i1 = Mesh % Edges(k) % NodeIndexes(2)

         ! Just use same index as later.
         l = k

         ! Coordinates for the polyline. 
         pCoords(6*l-5) = x(i0)
         pCoords(6*l-4) = y(i0)
         pCoords(6*l-3) = z(i0)

         ! Compute the difference already here. 
         pCoords(6*l-2) = x(i0)-x(i1)
         pCoords(6*l-1) = y(i0)-y(i1)
         pCoords(6*l-0) = z(i0)-z(i1)

         !ss = (x(i0)-x(i1))**2 + (y(i0)-y(i1))**2
       END DO
     END IF

     
     IF(PEs > 1 ) THEN        
       ALLOCATE(rPar(4*PEs))
       comm = ELMER_COMM_WORLD

       ! Communiate the bounding box first. This is only needed in (x,y) plane.
       IF(GotMaxDist) THEN
         rPar = -HUGE(rPar)
         rPar(4*MyPe-3:4*MyPe) = RidgeData(MyPe) % MeshBB
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, rPar, 4*PEs, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
         DO i=1,PEs
           IF(i==MyPe) CYCLE
           RidgeData(i) % MeshBB = rPar(4*i-3:4*i)
         END DO
         CALL MPI_BARRIER( comm, ierr )         
         DEALLOCATE(rPar)
       END IF

       ! Calculate and sent the count of edges to send. 
       n_max = 0
       DO i=1,PEs
         IF(i==MyPe) CYCLE
         j = MyPe
         proc = i-1

         IF(GotMaxDist) THEN
           n_send = 0

           ! This communication check is symmetric:
           ! If "i" does not need to sent to "j" then also vice versa!
           IF(RidgeData(i) % MeshBB(1) > RidgeData(j) % MeshBB(2) + MaxDist ) CYCLE
           IF(RidgeData(j) % MeshBB(1) > RidgeData(i) % MeshBB(2) + MaxDist ) CYCLE
           IF(RidgeData(i) % MeshBB(3) > RidgeData(j) % MeshBB(4) + MaxDist ) CYCLE
           IF(RidgeData(j) % MeshBB(3) > RidgeData(i) % MeshBB(4) + MaxDist ) CYCLE

           ! Compute how many entries to sent from "j" to "i"
           ! "i" cannot compute this but it knows to expect something.
           DO k=1, Mesh % NumberOfEdges                           
             i0 = Mesh % Edges(k) % NodeIndexes(1)
             i1 = Mesh % Edges(k) % NodeIndexes(2)

             IF(RidgeData(i) % MeshBB(1) > x(i0) + MaxDist ) CYCLE
             IF(x(i0) > RidgeData(i) % MeshBB(2) + MaxDist  ) CYCLE
             IF(RidgeData(i) % MeshBB(3) > y(i0) + MaxDist ) CYCLE
             IF(y(i0) > RidgeData(i) % MeshBB(4) + MaxDist ) CYCLE
             n_send = n_send+1
           END DO
         ELSE
           n_send = Mesh % NumberOfEdges
         END IF

         RidgeData(i) % n_send = n_send
         n_max = MAX(n_max, n_send)

         CALL MPI_SEND( n_send, 1, MPI_INTEGER, proc, 999, comm, ierr )
       END DO

       ! Recieve the counts of edges to recieve.
       DO i=1,PEs
         IF(i==MyPe) CYCLE
         proc = i-1
         CALL MPI_RECV( n_recv, 1, MPI_INTEGER, proc, 999, comm, status, ierr )
         RidgeData(i) % n_recv = n_recv
         
         IF(n_recv > 0) THEN             
           ALLOCATE(RidgeData(i) % Coords(6*n_recv))
           RidgeData(i) % Coords = 0.0_dp
         END IF
       END DO

       CALL MPI_BARRIER( comm, ierr )          

       k = SUM( RidgeData(1:PEs) % n_send ) 
       l = SUM( RidgeData(1:PEs) % n_recv ) 

       IF(InfoActive(20)) THEN
         PRINT *,'Total number of edges to send and recieve: ',MyPe, k, l, n_max
       END IF

       ALLOCATE(SendCoords(6*n_max))
       SendCoords = 0.0_dp
       pCoords => SendCoords 

       ! Send the edge coordinates to other partitions. 
       DO i=1,PEs
         j = MyPe
         IF(i==MyPe) CYCLE
         proc = i-1

         n_send = RidgeData(i) % n_send
         IF(n_send > 0 ) THEN
           ! Tabulate the entries to be sent.
           l = 0
           DO k=1, Mesh % NumberOfEdges
             i0 = Mesh % Edges(k) % NodeIndexes(1)
             i1 = Mesh % Edges(k) % NodeIndexes(2)

             IF(GotMaxDist) THEN
               IF(RidgeData(i) % MeshBB(1) > x(i0) + MaxDist ) CYCLE
               IF(x(i0) > RidgeData(i) % MeshBB(2) + MaxDist  ) CYCLE
               IF(RidgeData(i) % MeshBB(3) > y(i0) + MaxDist ) CYCLE
               IF(y(i0) > RidgeData(i) % MeshBB(4) + MaxDist ) CYCLE
             END IF

             ! Coordinates for the polyline. 
             l = l+1

             pCoords(6*l-5) = x(i0)
             pCoords(6*l-4) = y(i0)
             pCoords(6*l-3) = z(i0)
             pCoords(6*l-2) = x(i0)-x(i1)
             pCoords(6*l-1) = y(i0)-y(i1)
             pCoords(6*l-0) = z(i0)-z(i1)
           END DO

           ! Sent data from partition MyPe to i           
           CALL MPI_BSEND( SendCoords, 6*n_send, MPI_DOUBLE_PRECISION,proc, 1001, comm, ierr )
         END IF
       END DO

       ! Recieve the edge coordinates from other partitions. 
       DO i=1,PEs
         IF(i==MyPe) CYCLE
         proc = i-1         
         n_recv = RidgeData(i) % n_recv
         IF(n_recv > 0) THEN
           ! Recieve data from partition i to MyPe
           CALL MPI_RECV( RidgeData(i) % Coords, 6*n_recv, MPI_DOUBLE_PRECISION, proc, 1001, comm, status, ierr )
           pCoords => RidgeData(i) % Coords
         END IF
       END DO
       CALL MPI_BARRIER( comm, ierr )          
     END IF
     
   END SUBROUTINE PopulateRidge


   !------------------------------------------------------------------------------
   !> Release the ridge since it can be quite a memory hog consisting of edges from
   !> many partitions.
   !------------------------------------------------------------------------------
   SUBROUTINE ReleaseRidge()
     !------------------------------------------------------------------------------      
     INTEGER :: i,n_recv
     INTEGER :: PEs
     !------------------------------------------------------------------------------
     IF(.NOT. ALLOCATED(RidgeData)) THEN
       CALL Fatal('ReleaseRidge','Cannot release an unallocated RidgeData!')
     END IF
     
     DO i=1,ParEnv % PEs 
       n_recv = RidgeData(i) % n_recv
       IF(n_recv > 0) THEN
         DEALLOCATE(RidgeData(i) % Coords)
       END IF
     END DO
     DEALLOCATE(RidgeData)
     
   END SUBROUTINE ReleaseRidge

   
   !------------------------------------------------------------------------------
   !> Compute the shading angle in parallel. This has a little bit different
   !> routine that the local mesh version that uses the existing edges.
   !------------------------------------------------------------------------------   
   SUBROUTINE SetShadingAngleParallel()

     REAL(KIND=dp) :: x1,y1,z1,maxangle
     REAL(KIND=dp) :: x0,y0,z0,dz
     REAL(KIND=dp) :: b(2),c(2)
     INTEGER :: i,j,k,t,pe,myPe
     LOGICAL :: Singular
     REAL(KIND=dp), POINTER :: pCoords(:)

     MyPe = ParEnv % MyPe + 1

     ! node which might be shaded by the edge.
     DO j=1, Mesh % NumberOfNodes
       k = SunPerm(j)
       IF(k==0) CYCLE

       x1 = x(j)
       y1 = y(j)
       z1 = z(j)

       ! This has a value that may have been set by the 
       maxangle = SunAngle(dofs*(k-1)+idof)       
       
       ! This is a stupid N^2 algo where we have nested loops over edges and nodes.
       DO i = 1, ParEnv % PEs
         IF(RidgeData(i) % n_recv == 0) CYCLE

         IF(EnforceSerial .AND. i /= MyPe ) CYCLE         
         pCoords => RidgeData(i) % Coords
         
         DO t = 1, RidgeData(i) % n_recv
           Amat(1,2) = pCoords(6*t-2)
           Amat(2,2) = pCoords(6*t-1)

           b(1) = pCoords(6*t-5)-x1
           b(2) = pCoords(6*t-4)-y1

           CALL Solve2x2(Amat,c,b,Singular)
           
           ! The code is faster when all IF's are on same line.
           IF(Singular .OR. c(1) < Eps .OR. c(2) < -Eps .OR. c(2) > 1.0_dp+Eps ) CYCLE

           dz = ( pCoords(6*t-3) - c(2) * pCoords(6*t) ) - z1                   
           maxangle = MAX(maxangle,dz/c(1))
         END DO
       END DO

       SunAngle(dofs*(k-1)+idof) = maxangle
     END DO     

   END SUBROUTINE SetShadingAngleParallel

   
!------------------------------------------------------------------------------
END SUBROUTINE SunAngleSolver
!------------------------------------------------------------------------------


