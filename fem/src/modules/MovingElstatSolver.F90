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
! *  Email:   Peter.Rabackcsc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 8.2.2006
! * 
! ******************************************************************************/
    
!------------------------------------------------------------------------------
!>  Solver for the computation the capacitance and sensitivity of capacitance for 
!>  a rigid system with rotations and translations. The solver includes features 
!> of electrostatics and mesh solver within one solver and its usability is limited
!> to cartesian coordinates. Probably of very limited use.
!> \ingroup Solvers
!------------------------------------------------------------------------------
    SUBROUTINE MovingElstatSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

     USE Types
     USE Lists
     USE Integration
     USE ElementDescription
     USE Differentials
     USE SolverUtils
     USE ElementUtils
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
     TYPE(Matrix_t), POINTER  :: StiffMatrix
     TYPE(Element_t), POINTER :: CurrentElement, Parent
     TYPE(Variable_t), POINTER :: Var
     TYPE(Nodes_t) :: ElementNodes, ParentNodes

     REAL (KIND=DP), POINTER :: Coordinate(:), Displacement(:,:), ForceVector(:), Coords(:), Component(:), &
         xorig(:), yorig(:), zorig(:), ElemPotential(:), Parray(:,:)
     REAL (KIND=DP), ALLOCATABLE ::  LocalStiffMatrix(:,:), LocalForce(:), &
         xnew(:), ynew(:), znew(:), MatrixValues(:), &
         Basis(:),dBasisdx(:,:),ddBasisddx(:,:,:), AplacResults(:,:), &
         ParentBasis(:),ParentdBasisdx(:,:),xelem(:),yelem(:),zelem(:)

#ifdef USE_ISO_C_BINDINGS
     REAL (KIND=DP) :: Norm, at, st, val, Normal(3), MinMaxEdge(3,2), Eps, aid, &
         da, ds, Dx(6), Dx0(6), Rotate(3,3), Point(3), Displace(3), Translate(3), Center(3), &
         TotCapac, ElemCapac, PermittivityOfVacuum, Base(6,6), Limits(6,2), Amplitude(6), &
         TotForce(3), TotMoment(3), TotArea, TotCharge, TotTime, MinMaxMoving(3,2), &
         SectionCapac(3,2), Dist, LengthScale
#else
     REAL (KIND=DP) :: Norm, at, st, CPUTime, val, Normal(3), MinMaxEdge(3,2), Eps, aid, &
         da, ds, Dx(6), Dx0(6), Rotate(3,3), Point(3), Displace(3), Translate(3), Center(3), &
         TotCapac, ElemCapac, PermittivityOfVacuum, Base(6,6), Limits(6,2), Amplitude(6), &
         TotForce(3), TotMoment(3), TotArea, TotCharge, TotTime, MinMaxMoving(3,2), &
         SectionCapac(3,2), Dist, LengthScale
#endif
     INTEGER, POINTER :: NodeIndexes(:), CoordinatePerm(:)
     INTEGER, ALLOCATABLE :: ElementNormals(:)
     INTEGER :: i, j, k, n, t, istat, bf_id, MeshNodes, DIM, coord, Visited = 0, &
         Intervals(6), Counter(6), i1, i2, i3, i4, i5, i6, TableDim, pn, Cases, &
         VerticalCoordinate, PeriodicCoordinate
     LOGICAL :: AllocationsDone = .FALSE., gotIt, stat, SaveCoords, MovingEdge(3,2), DoRotate(3), &
         CalculateMoment, CalculateForce, Debug, DoTranslate(3), MappedCoordinates, SetPoint
     LOGICAL, POINTER :: MovingNodes(:)
     CHARACTER(LEN=MAX_NAME_LEN) :: FileName

     CHARACTER(LEN=MAX_NAME_LEN) :: VariableName

     SAVE LocalStiffMatrix, LocalForce, ElementNodes, AllocationsDone, &
         xnew, ynew, znew, ElementNormals, ParentNodes, &
         Coords, SaveCoords, MinMaxEdge, MovingEdge, Eps, MeshNodes, &
         MatrixValues, Visited, Basis, dBasisdx, ddBasisddx, ElemPotential, &
         ParentBasis, ParentdBasisdx, MovingNodes, MinMaxMoving, &
         xelem, yelem, zelem, CalculateMoment, CalculateForce, VerticalCoordinate, &
         PeriodicCoordinate, LengthScale, Filename

!------------------------------------------------------------------------------

     CALL Info( 'MovingElstatSolver', '-------------------------------------',Level=4 )
     CALL Info( 'MovingElstatSolver', 'Rigid body coordinate update solver:  ', Level=4 )
     CALL Info( 'MovingElstatSolver', '-------------------------------------',Level=4 )


!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     Coordinate     => Solver % Variable % Values
     CoordinatePerm => Solver % Variable % Perm

     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS

     Norm = Solver % Variable % Norm
     DIM = CoordinateSystemDimension()

     xorig => Solver % Mesh % Nodes % x
     yorig => Solver % Mesh % Nodes % y
     zorig => Solver % Mesh % Nodes % z

!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN
       N = Model % MaxElementNodes
       MeshNodes = SIZE(Solver % Mesh % Nodes % x) 
 
       ALLOCATE( ElementNodes % x(N),   &
                 ElementNodes % y(N),   &
                 ElementNodes % z(N),   &
                 ParentNodes % x(N),   &
                 ParentNodes % y(N),   &
                 ParentNodes % z(N),   &
                 ElemPotential(N), &
                 LocalForce(N),         & 
                 LocalStiffMatrix(N,N), &
                 Basis(n), &
                 dBasisdx(n,3), &
                 ddBasisddx(n,3,3), &
                 ParentBasis(n), &
                 ParentdBasisdx(n,3), &
                 xelem(n), &
                 yelem(n), &
                 zelem(n), &
                 xnew(MeshNodes), &
                 ynew(MeshNodes), &
                 znew(MeshNodes), &
                 MovingNodes(MeshNodes), &
                 ElementNormals(Solver % Mesh % NumberOfBoundaryElements), &
                 MatrixValues(SIZE(StiffMatrix % Values)), &
                 STAT=istat )
 
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'MovingElstatSolver', 'Memory allocation error 1' )
       END IF

       !------------------------------------------------------------------------------
       !  Check whether the new coordinates should be saved 
       !------------------------------------------------------------------------------       
       SaveCoords = ListGetLogical(Solver % Values,'Save Displacements',GotIt)
       IF(SaveCoords) THEN
         ALLOCATE(Coords(3*MeshNodes))
         Coords = 0.0d0
         CALL VariableAddVector( Solver % Mesh % Variables, Solver % Mesh, Solver, 'dCoords',3,Coords)
       END IF

       AllocationsDone = .TRUE.      
     END IF


     TotTime  = CPUTime()

     IF(Visited == 0) THEN

       !------------------------------------------------------------------------------
       !    Read some control stuff
       !------------------------------------------------------------------------------

       CalculateMoment = ListGetLogical(Solver % Values,'Calculate Moment',GotIt)
       CalculateForce = CalculateMoment .OR. ListGetLogical(Solver % Values,'Calculate Force',GotIt)     

       VerticalCoordinate = ListGetInteger(Solver % Values,'Vertical Coordinate',GotIt)
       IF(.NOT. GotIt) VerticalCoordinate = 3
       PeriodicCoordinate = ListGetInteger(Solver % Values,'Periodic Coordinate',GotIt)
       IF(.NOT. GotIt) PeriodicCoordinate = 1
       
       LengthScale = ListGetConstReal(Solver % Values,'Length Scale',GotIt)
       IF(.NOT. GotIt) LengthScale = 1.0d0
       
       FileName = ListGetString(Solver % Values,'Filename',GotIt )
       IF(.NOT. GotIt) FileName = 'lumped.dat'

       !------------------------------------------------------------------------------
       !    Check directions of surface normals and memorize the moving nodes
       !------------------------------------------------------------------------------
       
       ElementNormals = 0
       MovingNodes = .FALSE.
       
       DO t = Solver % Mesh % NumberOfBulkElements + 1, &
           Solver % Mesh % NumberOfBulkElements + Solver % Mesh % NumberOfBoundaryElements
         
         CurrentElement => Solver % Mesh % Elements(t)
         NodeIndexes => CurrentElement % NodeIndexes
         n = CurrentElement % TYPE % NumberOfNodes

         ElementNodes % x(1:n) = xorig(NodeIndexes)
         ElementNodes % y(1:n) = yorig(NodeIndexes)
         ElementNodes % z(1:n) = zorig(NodeIndexes)

         Model % CurrentElement => CurrentElement
         j = t - Solver % Mesh % NumberOfBulkElements 

         Normal = NormalVector(CurrentElement,ElementNodes,0.0d0,0.0d0) 

         DO coord = 1, 3
           IF(ABS(Normal(coord)) > 1.0 - 1.0d-3) ElementNormals(j) = coord
         END DO
            
         DO j=1,Model % NumberOfBCs
           IF ( CurrentElement % BoundaryInfo % Constraint == &
               Model % BCs(j) % Tag ) THEN
             
             IF ( ListGetLogical(Model % BCs(j) % Values,'Moving Boundary',gotIt) ) THEN
               MovingNodes(NodeIndexes) = .TRUE.   
             END IF
           END IF
         END DO

       END DO

       !------------------------------------------------------------------------------
       !  Find the edges of the cartesian domain
       !------------------------------------------------------------------------------
       MinMaxEdge(:,1) = HUGE(MinMaxEdge(:,1))
       MinMaxEdge(:,2) = -HUGE(MinMaxEdge(:,2))
       MinMaxMoving = MinMaxEdge

       DO i = 1, MeshNodes
         IF(CoordinatePerm(i) == 0) CYCLE
         MinMaxEdge(1,1) = MIN( MinMaxEdge(1,1), xorig(i) )
         MinMaxEdge(1,2) = MAX( MinMaxEdge(1,2), xorig(i) )
         MinMaxEdge(2,1) = MIN( MinMaxEdge(2,1), yorig(i) )
         MinMaxEdge(2,2) = MAX( MinMaxEdge(2,2), yorig(i) )
         MinMaxEdge(3,1) = MIN( MinMaxEdge(3,1), zorig(i) )
         MinMaxEdge(3,2) = MAX( MinMaxEdge(3,2), zorig(i) )
         
         IF(.NOT. MovingNodes(i) ) CYCLE
         MinMaxMoving(1,1) = MIN( MinMaxMoving(1,1), xorig(i) )
         MinMaxMoving(1,2) = MAX( MinMaxMoving(1,2), xorig(i) )
         MinMaxMoving(2,1) = MIN( MinMaxMoving(2,1), yorig(i) )
         MinMaxMoving(2,2) = MAX( MinMaxMoving(2,2), yorig(i) )
         MinMaxMoving(3,1) = MIN( MinMaxMoving(3,1), zorig(i) )
         MinMaxMoving(3,2) = MAX( MinMaxMoving(3,2), zorig(i) )        
       END DO

       Eps = 1.0d-8 * SQRT( (MinMaxEdge(1,1)-MinMaxEdge(1,2))**2 + &
           (MinMaxEdge(2,1)-MinMaxEdge(2,2))**2 + &
           (MinMaxEdge(3,1)-MinMaxEdge(3,2))**2 )

       !------------------------------------------------------------------------------
       !  If edge is not moving with Displacement it is assumed to be fixed
       !  Also periodic coordinate direction is assumed to be fixed.
       !------------------------------------------------------------------------------       
       MovingEdge = .FALSE.
       DO i=1,3
         DO j=1,2
           IF( ABS( MinMaxEdge(i,j) - MinMaxMoving(i,j) ) < Eps) MovingEdge(i,j) = .TRUE.
         END DO
       END DO

       IF(PeriodicCoordinate > 0) THEN
         MovingEdge(PeriodicCoordinate,1:2) = .FALSE.
       END IF

       PRINT *,'Epsilon',Eps
       PRINT *,'MinMaxEdge',MinMaxEdge
       PRINT *,'MinMaxMoving',MinMaxMoving
       PRINT *,'Moving',MovingEdge

       !------------------------------------------------------------------------------
       !    Assembly the Laplace operator and save it for repetitious use
       !------------------------------------------------------------------------------       
       at  = CPUTime()
       LocalForce = 0.0d0
       CALL InitializeToZero( StiffMatrix, ForceVector )

       DO t = 1, Solver % NumberOfActiveElements
         
         CurrentElement => Solver % Mesh % Elements(Solver % ActiveElements(t))
         NodeIndexes => CurrentElement % NodeIndexes
         n = CurrentElement % TYPE % NumberOfNodes
         
         ElementNodes % x(1:n) = xorig(NodeIndexes)
         ElementNodes % y(1:n) = yorig(NodeIndexes)
         ElementNodes % z(1:n) = zorig(NodeIndexes)

         Model % CurrentElement => CurrentElement

         CALL LaplaceCompose( LocalStiffMatrix,CurrentElement,n,ElementNodes)

         CALL UpdateGlobalEquations( StiffMatrix,LocalStiffMatrix, &
             ForceVector, LocalForce,n,1,CoordinatePerm(NodeIndexes) )
       END DO

       CALL FinishAssembly( Solver,ForceVector )
       
       WRITE( Message, * ) 'Bulk Assembly (s)     :',CPUTime() - at
       CALL Info( 'MovingElstatSolver', Message, Level=4 )

       MatrixValues = StiffMatrix % Values
     END IF
     
     !------------------------------------------------------------------------------
     !  Get the rules for rigid displacement
     !------------------------------------------------------------------------------      
     Center(1) = ListGetConstReal( Model % Simulation,'Moment About 1', GotIt )
     IF(.NOT. GotIt) Center(1) = ListGetConstReal( Solver % Values,'Moment About 1', GotIt )
     IF(.NOT. GotIt) Center(1) = ListGetConstReal( Model % Simulation,'res: Center of Mass 1')
     Center(2) = ListGetConstReal( Model % Simulation,'Moment About 2', GotIt )
     IF(.NOT. GotIt) Center(2) = ListGetConstReal( Solver % Values,'Moment About 2', GotIt )
     IF(.NOT. GotIt) Center(2) = ListGetConstReal( Model % Simulation,'res: Center of Mass 2')
     Center(3) = ListGetConstReal( Model % Simulation,'Moment About 3', GotIt )
     IF(.NOT. GotIt) Center(3) = ListGetConstReal( Solver % Values,'Moment About 3', GotIt )
     IF(.NOT. GotIt) Center(3) = ListGetConstReal( Model % Simulation,'res: Center of Mass 3')
     
     Base = 0.0d0
     Limits = 0.0d0
     Intervals = 1
     TableDim = 0
     Cases = 1

     DO i=1,6
       PArray => ListGetConstRealArray( Solver % Values,'Lumping Basis '//CHAR(i+ICHAR('0')),GotIt)
       IF(GotIt) THEN
         Base(i,1:6) = Parray(1:6,1)
       ELSE
         Base(i,i) = 1.0d0
       END IF

       PArray => ListGetConstRealArray( Solver % Values,'Lumping Limits '//CHAR(i+ICHAR('0')),GotIt)
       IF(GotIt) Limits(i,1:2) = Parray(1:2,1)

       Intervals(i) = ListGetInteger( Solver % Values,'Lumping Points '//CHAR(i+ICHAR('0')),GotIt)
       Intervals(i) = MAX(1,Intervals(i))
       Cases = Cases * Intervals(i)
       IF(Intervals(i) > 1) TableDim = TableDim + 1
     END DO

     IF(TableDim == 1) THEN
       ALLOCATE( AplacResults(Cases, 2) )
     END IF

     !------------------------------------------------------------------------------
     !  Write a file with metainformation of the lumping procedure
     !------------------------------------------------------------------------------           
     OPEN (10, FILE=TRIM(FileName)//'.info')
     WRITE(10,*) 'Information for lumped electrostatics'
     WRITE(10,*) '*************************************'
     WRITE(10,*)
     WRITE(10,*) 'Lumping base'
     DO i=1,6
       WRITE(10,*) Base(i,1:6)
     END DO

     WRITE(10,*)
     WRITE(10,*) 'Lumping Interval and Limits'
     DO i=1,6
       WRITE(10,*) Intervals(i), Limits(i,1:2)
     END DO
     
     WRITE(10,*)
     WRITE(10,*) 'Center of coordinate'
     WRITE(10,*) Center

     j = 0
     WRITE(10,*)
     WRITE(10,*) 'Columns in file: ',TRIM(FileName)     
     DO i=1,6
       IF(Intervals(i) > 1) THEN
         j = j + 1
         WRITE(10,*) j,': Amplitude of base',i
       END IF
     END DO
     WRITE(10,*) j+1,': C'
     WRITE(10,*) j+2,': Cup'
     WRITE(10,*) j+3,': Cdown'
     IF(CalculateForce) THEN
       WRITE(10,*) j+4,': C (surface)'
       WRITE(10,*) j+5,': dC/dx'
       WRITE(10,*) j+6,': dC/dy'
       WRITE(10,*) j+7,': dC/dz'
     END IF
     IF(CalculateMoment) THEN
       WRITE(10,*) j+8,': Mx'
       WRITE(10,*) j+9,': My'
       WRITE(10,*) j+10,': Mz'
     END IF
      
     WRITE(10,*)
     WRITE(10,*) 'Degrees of Freedom',SIZE(Coordinate)

     CLOSE(10)


     !------------------------------------------------------------------------------
     !  Loop over all possible combinations of displacement
     !------------------------------------------------------------------------------           
     Cases = 0
     DO i1 = 1, Intervals(1)
     DO i2 = 1, Intervals(2)
     DO i3 = 1, Intervals(3)
     DO i4 = 1, Intervals(4)
     DO i5 = 1, Intervals(5)
     DO i6 = 1, Intervals(6)
       
       Counter(1) = i1
       Counter(2) = i2
       Counter(3) = i3
       Counter(4) = i4
       Counter(5) = i5
       Counter(6) = i6

       Cases = Cases + 1
         
       Dx = 0.0d0
       DO i=1,6
         Amplitude(i) = Limits(i,1) 
         IF(Counter(i) > 1) THEN
           Amplitude(i) = Amplitude(i) + (Counter(i)-1) * (Limits(i,2)-Limits(i,1)) / (Intervals(i)-1) 
         END IF         
         Dx = Dx + Base(i,:) * Amplitude(i)
       END DO

       Translate(1:3) = Dx(1:3)
       Rotate = 0.0d0
       Rotate(1,2) = Dx(6)
       Rotate(1,3) = -Dx(5)
       Rotate(2,3) = Dx(4)
       
       DO i=1,3
         DO j=i+1,3
           Rotate(j,i) = -Rotate(i,j)
         END DO
       END DO
       
       DoTranslate = (ABS(Dx(1:3)) > 1.0d-20) 
       DO i=1,3
         DoRotate(i) = ANY( ABS(Rotate(i,1:3)) > 1.0d-20 )
       END DO

!------------------------------------------------------------------------------
!    Go through the coordinates separately and set the following BCs
!    1) Set nodes attached to rigid bodies
!    2) Set the normal component of all other boundaries to zero
!    3) Set the normal component of the whole frame to zero
!------------------------------------------------------------------------------
         
       MappedCoordinates = .FALSE.

       DO coord = 1, DIM
         
         IF(.NOT. ( DoTranslate(coord) .OR. DoRotate(coord) ) ) THEN
           IF(coord == 1) xnew = xorig
           IF(coord == 2) ynew = yorig
           IF(coord == 3) znew = zorig
           WRITE(Message,'(a,i1)' ) 'No need to map coordinate ',coord        
           CALL Info( 'MovingElstatSolver', Message, Level=4 )         
           CYCLE
         END IF
         
         MappedCoordinates = .TRUE.
         WRITE(Message,'(a,i1)' ) 'Mapping coordinate ',coord
         CALL Info( 'MovingElstatSolver', Message, Level=4 )
         at  = CPUTime()
         
         StiffMatrix % Values = MatrixValues
         ForceVector = 0.0d0         

         !------------------------------------------------------------------------------
         ! Fix displacements to normal direction 
         ! and for the whole moving part in case of rotation
         !------------------------------------------------------------------------------           

         DO t = Solver % Mesh % NumberOfBulkElements + 1, &
             Solver % Mesh % NumberOfBulkElements + Solver % Mesh % NumberOfBoundaryElements
           
           CurrentElement => Solver % Mesh % Elements(t)
           NodeIndexes => CurrentElement % NodeIndexes
           n = CurrentElement % TYPE % NumberOfNodes
           k = t - Solver % Mesh % NumberOfBulkElements 
           
           DO j=1,Model % NumberOfBCs
             IF ( CurrentElement % BoundaryInfo % Constraint == &
                 Model % BCs(j) % Tag ) THEN

               IF ( ListGetLogical(Model % BCs(j) % Values,'Moving Boundary',gotIt) ) THEN

                 IF(ElementNormals(k) == 0 .OR. ElementNormals(k) == coord) THEN
                   ! Fix displacement in normal direction and to unclear direction
                   SetPoint = .TRUE.
                 ELSE IF( DoRotate(ElementNormals(k)) ) THEN
                   ! If normal direction is not constant fix the others too
                   SetPoint = .TRUE.
                 ELSE
                   SetPoint = .FALSE.
                 END IF

                 IF( SetPoint ) THEN
                   DO i = 1, n
                     Point(1) = xorig(NodeIndexes(i))
                     Point(2) = yorig(NodeIndexes(i))
                     Point(3) = zorig(NodeIndexes(i))
                     
                     Point = Point - Center
                     Displace = Translate + MATMUL(Rotate,Point)           
                     
                     val = Displace(coord)                    

                     CALL SetDirichletPoint( StiffMatrix, ForceVector,1,1, &
                         CoordinatePerm, NodeIndexes(i), val)           
                   END DO
                 END IF
               ELSE IF( ListGetLogical(Model % BCs(j) % Values,'Fixed Boundary',gotIt) ) THEN
                 IF( ElementNormals(k) == 0 .OR. ElementNormals(k) == coord ) THEN
                   val = 0.0d0
                   DO i = 1, n
                     CALL SetDirichletPoint( StiffMatrix, ForceVector,1,1, &
                         CoordinatePerm, NodeIndexes(i), val)                
                   END DO
                 END IF
               END IF
               
             END IF
           END DO
         END DO
         
         !------------------------------------------------------------------------------
         !  Fix the movement of the whole frame in the normal direction
         !------------------------------------------------------------------------------        

         DO i = 1, MeshNodes
           IF( CoordinatePerm(i) == 0 ) CYCLE

           Point(1) = xorig(i)
           Point(2) = yorig(i)
           Point(3) = zorig(i)
                     
           DO j= 1,2
             IF( ABS( MinMaxEdge(coord,j) - Point(coord) ) < Eps) THEN

               ! In case of displaced frame only account for the translation
               ! This will ensure that the box remains a box 
               IF( MovingEdge(coord,j)) THEN
                 val = Translate(coord)               
               ELSE
                 val = 0.0d0
               END IF
               
               CALL SetDirichletPoint( StiffMatrix, ForceVector,1,1, &
                   CoordinatePerm, i, val)                          
             END IF
           END DO
         END DO
         
         
         at = CPUTime() - at
         WRITE( Message, * ) 'Boundary Assembly (s) :',at
         CALL Info( 'MovingElstatSolver', Message, Level=4 )
         
         !------------------------------------------------------------------------------
         !    Solve the system 
         !------------------------------------------------------------------------------
         st = CPUTime()
         CALL SolveSystem( StiffMatrix, ParMatrix, ForceVector, &
             Coordinate, Norm, 1, Solver )
         st = CPUTime() - st
         WRITE( Message, * ) 'Solve (s)             :',st
         CALL Info( 'MovingElstatSolver', Message, Level=4 )

         !------------------------------------------------------------------------------
         ! Add the displacement to the coordinate values
         !------------------------------------------------------------------------------       
         DO i=1,MeshNodes
           j = CoordinatePerm(i)
           IF(j > 0) THEN
             IF(coord == 1) xnew(i) = xorig(i) + Coordinate(j)
             IF(coord == 2) ynew(i) = yorig(i) + Coordinate(j)
             IF(coord == 3) znew(i) = zorig(i) + Coordinate(j)
             IF(SaveCoords) Coords(DIM*(i-1)+coord) = Coordinate(j)
           END IF
         END DO
         
       END DO
       
       CALL Info('MovingElstatSolver','Coordinates updated')
             
       Visited = Visited + 1
       
       !------------------------------------------------------------------------------
       ! Assembly the electrostatic Equation
       !------------------------------------------------------------------------------       
       at = CPUTime()

       IF(.NOT. MappedCoordinates) THEN
         CALL Info('MovingElstatSolver','Using original matrix for potential equation')
         StiffMatrix % Values = MatrixValues
       ELSE
         CALL Info('MovingElstatSolver','Assembling the potential equation')
         LocalForce = 0.0d0
         CALL InitializeToZero( StiffMatrix, ForceVector )
         
         DO t = 1, Solver % NumberOfActiveElements
           
           CurrentElement => Solver % Mesh % Elements(Solver % ActiveElements(t))
           NodeIndexes => CurrentElement % NodeIndexes
           n = CurrentElement % TYPE % NumberOfNodes
           
           ElementNodes % x(1:n) = xnew(NodeIndexes)
           ElementNodes % y(1:n) = ynew(NodeIndexes)
           ElementNodes % z(1:n) = znew(NodeIndexes)
           
           Model % CurrentElement => CurrentElement
           
           CALL LaplaceCompose( LocalStiffMatrix,CurrentElement,n,ElementNodes)
           
           CALL UpdateGlobalEquations( StiffMatrix,LocalStiffMatrix, &
               ForceVector, LocalForce,n,1,CoordinatePerm(NodeIndexes) )
         END DO
       END IF

       
       !------------------------------------------------------------------------------
       !  Set the Boundary conditions to be either zero or one
       !------------------------------------------------------------------------------
       
       DO t = Solver % Mesh % NumberOfBulkElements + 1, &
           Solver % Mesh % NumberOfBulkElements + Solver % Mesh % NumberOfBoundaryElements
         
         CurrentElement => Solver % Mesh % Elements(t)
         NodeIndexes => CurrentElement % NodeIndexes
         n = CurrentElement % TYPE % NumberOfNodes        
        
         DO j=1,Model % NumberOfBCs
           IF ( CurrentElement % BoundaryInfo % Constraint == &
               Model % BCs(j) % Tag ) THEN
             
             IF( ListGetLogical(Model % BCs(j) % Values, &
                 'Moving Boundary',gotIt) ) THEN
               val = 1.0d0
               DO i = 1, n                 
                CALL SetDirichletPoint( StiffMatrix, ForceVector,1,1, &
                     CoordinatePerm, NodeIndexes(i), val)                
               END DO
             ELSE IF( ListGetLogical(Model % BCs(j) % Values, &
                 'Fixed Boundary',gotIt) ) THEN
               val = 0.0d0
               DO i = 1, n                 
                 CALL SetDirichletPoint( StiffMatrix, ForceVector,1,1, &
                     CoordinatePerm, NodeIndexes(i), val)                
               END DO
             END IF
             
           END IF
         END DO
       END DO
       
       
       CALL FinishAssembly( Solver,ForceVector )

       !------------------------------------------------------------------------------
       !    Dirichlet boundary conditions (symmetry)
       !------------------------------------------------------------------------------
       CALL SetDirichletBoundaries( Model,StiffMatrix,ForceVector, &
           ComponentName(Solver % Variable),1,1,CoordinatePerm )

       
       WRITE( Message, * ) 'Bulk Assembly (s)     :',CPUTime() - at
       CALL Info( 'MovingElstatSolver', Message, Level=4 )
       
       !------------------------------------------------------------------------------
       !    Solve the system 
       !------------------------------------------------------------------------------
       st = CPUTime()
       CALL SolveSystem( StiffMatrix, ParMatrix, ForceVector, Coordinate, Norm, 1, Solver )
       st = CPUTime() - st
       WRITE( Message, * ) 'Solve (s)             :',st
       CALL Info( 'MovingElstatSolver', Message, Level=4 )
       
       !------------------------------------------------------------------------------
       ! Calculate Electric Capacitance From Volume Integral
       !------------------------------------------------------------------------------
       
       TotCapac = 0.0d0       
       SectionCapac = 0.0d0

       DO t = 1, Solver % NumberOfActiveElements         
         CurrentElement => Solver % Mesh % Elements(Solver % ActiveElements(t))
         NodeIndexes => CurrentElement % NodeIndexes
         n = CurrentElement % TYPE % NumberOfNodes
         
         ElementNodes % x(1:n) = xnew(NodeIndexes)
         ElementNodes % y(1:n) = ynew(NodeIndexes)
         ElementNodes % z(1:n) = znew(NodeIndexes)
         
         Model % CurrentElement => CurrentElement         
         ElemPotential(1:n) = Coordinate(CoordinatePerm(NodeIndexes))
         
         ElemCapac = ElementEnergy(CurrentElement, n, ElementNodes, ElemPotential)
         TotCapac = TotCapac + ElemCapac
       END DO

       PermittivityOfVacuum = ListGetConstReal( Model % Constants, &
           'Permittivity Of Vacuum',gotIt )
       IF ( .NOT.gotIt ) PermittivityOfVacuum = 1.0d0

       TotCapac = PermittivityOfVacuum * TotCapac
       PRINT *,'Dx',Dx
       PRINT *,'Capacitance',TotCapac
       PRINT *,'SectionCapac',SectionCapac
       
       !------------------------------------------------------------------------------
       ! Calculate Electric Force and Charge from Surface Integral
       !------------------------------------------------------------------------------

       IF(CalculateForce .OR. CalculateMoment) THEN
         TotForce = 0.0d0
         TotMoment = 0.0d0 
         TotArea = 0.0d0
         TotCharge = 0.0d0
         
         DO t = Solver % Mesh % NumberOfBulkElements + 1, &
             Solver % Mesh % NumberOfBulkElements + Solver % Mesh % NumberOfBoundaryElements
           
           CurrentElement => Solver % Mesh % Elements(t)
           
           DO j=1,Model % NumberOfBCs
             IF ( CurrentElement % BoundaryInfo % Constraint == Model % BCs(j) % Tag ) THEN
               
               IF( .NOT. ListGetLogical(Model % BCs(j) % Values,'Moving Boundary',gotIt)) CYCLE
               
               NodeIndexes => CurrentElement % NodeIndexes
               IF( ANY (CoordinatePerm(NodeIndexes) == 0 ) ) CYCLE
               
               n = CurrentElement % TYPE % NumberOfNodes        
               Model % CurrentElement => CurrentElement
               
               ElementNodes % x(1:n) = xnew(NodeIndexes)
               ElementNodes % y(1:n) = ynew(NodeIndexes)
               ElementNodes % z(1:n) = znew(NodeIndexes)
               
               !------------------------------------------------------------------------------
               !     Need parent element to determine the derivative at the surface
               !------------------------------------------------------------------------------
               Parent => CurrentElement % BoundaryInfo % Left    
               stat = ASSOCIATED( Parent )
               IF ( stat ) stat = ALL (CoordinatePerm(Parent % NodeIndexes) > 0 ) 
               
               IF ( .NOT. stat ) THEN
                 Parent => CurrentElement % BoundaryInfo % Right
                 stat = ASSOCIATED( Parent) 
                 IF(stat) stat = ALL (CoordinatePerm(Parent % NodeIndexes) > 0 ) 
               END IF
               
               IF(.NOT. stat) CYCLE
               
               pn = Parent % TYPE % NumberOfNodes      
               ParentNodes % x(1:pn) = xnew(Parent % NodeIndexes)
               ParentNodes % y(1:pn) = ynew(Parent % NodeIndexes)
               ParentNodes % z(1:pn) = znew(Parent % NodeIndexes)
               
               ElemPotential(1:pn) = Coordinate( CoordinatePerm( Parent % NodeIndexes ) )
               CALL ElectricForceIntegrate( TotForce, TotMoment, TotArea, TotCharge )

             END IF
           END DO
         END DO
         
         TotCharge = -PermittivityOfVacuum * TotCharge
         TotForce = 2.0d0 * PermittivityOfVacuum * TotForce
         TotMoment = PermittivityOfVacuum * TotMoment
         TotArea = TotArea

         PRINT *,'Charge',TotCharge
         PRINT *,'dC/dx',TotForce
         IF(CalculateMoment) PRINT *,'Moment',TotMoment
         PRINT *,'Area',TotArea
       END IF


       !------------------------------------------------------------------------------
       ! Save results in matrix format
       !------------------------------------------------------------------------------       
       IF(Cases == 1) THEN
         OPEN (10, FILE=TRIM(FileName))
       ELSE
         OPEN (10, FILE=TRIM(FileName), POSITION='APPEND')
       END IF

       DO i=1,6
         IF(Intervals(i) > 1) WRITE (10,'(ES22.12E3)',advance='no') Amplitude(i)
       END DO

       WRITE (10,'(ES20.10E3)',advance='no') TotCapac
       WRITE (10,'(2ES20.10E3)',advance='no') SectionCapac(VerticalCoordinate, 1:2)       
       IF(CalculateForce) THEN
         WRITE (10,'(ES20.10E3)',advance='no') TotCharge
         WRITE (10,'(3ES20.10E3)',advance='no') TotForce         
       END IF
       IF(CalculateMoment) THEN
         WRITE (10,'(3ES20.10E3)',advance='no') TotMoment
       END IF
       WRITE (10,*) ' '
       CLOSE(10)

       !------------------------------------------------------------------------------
       ! Save results in Aplac format
       !------------------------------------------------------------------------------       
       IF(TableDim == 1) THEN
         DO i=1,6
           IF(Intervals(i) > 1) AplacResults(Cases,1) = Amplitude(i)           
         END DO
         AplacResults(Cases,2) = TotCapac
       END IF      
       
     END DO
     END DO
     END DO
     END DO
     END DO
     END DO

     IF(TableDim == 1) THEN
       DO i=1,6
         IF(Intervals(i) > 1) EXIT
       END DO

       Dist = ListGetConstReal(Solver % Values,'Aplac Distance',GotIt)

       OPEN (10, FILE=TRIM(FileName)//'.aplac')
       
       WRITE(10,'(A)',advance='no') 'vector vh'
       DO i=1,Cases
         WRITE (10,'(ES16.8E2,A)',advance='no') &
             1.0d6 * LengthScale * AplacResults(i,1),'u'
       END DO
       WRITE(10,'(A)') ' '

       WRITE(10,'(A)',advance='no') 'vector vc'
       DO i=1,Cases
         WRITE (10,'(ES16.8E2,A)',advance='no') &
             1.0d12 * LengthScale**2 * AplacResults(i,2),'p'
       END DO
       WRITE(10,'(A)') ' '

       WRITE(10,'(A,ES16.8E2,A,I3)') 'GeneralTransducer "Elmer" nv 0 nz 0 nU 0 D=',&
           1.0d6*LengthScale*Dist,'u H=vh C=vc N=',Cases

     END IF
 
     
     TotTime  = CPUTime() - TotTime
     OPEN (10, FILE=TRIM(FileName)//'.info',POSITION='append')
     WRITE(10,*) 'Number of Cases',Cases
     WRITE(10,*) 'Total CPU Time (s)',TotTime


!------------------------------------------------------------------------------
 
   CONTAINS

 
!------------------------------------------------------------------------------
     SUBROUTINE LaplaceCompose( StiffMatrix,Element,n,Nodes )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: StiffMatrix(:,:)
       INTEGER :: n
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
 
       REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,A,L,x,y,z,aid
       LOGICAL :: Stat
       INTEGER :: i,j,p,q,t 
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
 
!------------------------------------------------------------------------------

       StiffMatrix = 0.0d0
 
!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints( Element )
 
       DO t=1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Element,Nodes,U,V,W,SqrtElementMetric, &
             Basis,dBasisdx,ddBasisddx,.FALSE. )
         S = IntegStuff % s(t) * SqrtElementMetric

!------------------------------------------------------------------------------
!        The Poisson equation
!------------------------------------------------------------------------------

         DO p=1,N
           DO q=1,N
             DO i=1,DIM
               StiffMatrix(p,q) = StiffMatrix(p,q) + s * dBasisdx(p,i) * dBasisdx(q,i)
             END DO
           END DO
         END DO

       END DO
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
     END SUBROUTINE LaplaceCompose
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
     FUNCTION ElementEnergy( Element,n,Nodes,ElemPotential) RESULT (Wetot)
!------------------------------------------------------------------------------
       INTEGER :: n
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
       REAL(KIND=dp), POINTER :: ElemPotential(:)
       REAL(KIND=dp) :: Wetot
!------------------------------------------------------------------------------
 
       REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,A,L,x,y,z,Grad(3),We
       LOGICAL :: Stat
       INTEGER :: i,j,p,q,t 
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
 
!------------------------------------------------------------------------------

       Wetot = 0.0d0
 
!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints( Element )
 
       DO t=1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
         S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Element,Nodes,U,V,W,SqrtElementMetric, &
                    Basis,dBasisdx,ddBasisddx,.FALSE. )
         S = S * SqrtElementMetric

         x = SUM( Basis(1:n) * Nodes % x(1:n) )
         y = SUM( Basis(1:n) * Nodes % y(1:n) )
         z = SUM( Basis(1:n) * Nodes % z(1:n) )         

!------------------------------------------------------------------------------
!        The Poisson equation
!------------------------------------------------------------------------------

         DO j = 1, DIM
           Grad(j) = SUM( dBasisdx(1:n,j) * ElemPotential(1:n) )
         END DO
         
         We = s * SUM( Grad(1:DIM) * Grad(1:DIM) )
         
         Wetot = Wetot + We

         IF(x < MinMaxMoving(1,1) ) SectionCapac(1,1) = SectionCapac(1,1) + We
         IF(x > MinMaxMoving(1,2) ) SectionCapac(1,2) = SectionCapac(1,2) + We
         IF(y < MinMaxMoving(2,1) ) SectionCapac(2,1) = SectionCapac(2,1) + We
         IF(y > MinMaxMoving(2,2) ) SectionCapac(2,2) = SectionCapac(2,2) + We
         IF(z < MinMaxMoving(3,1) ) SectionCapac(3,1) = SectionCapac(3,1) + We
         IF(z > MinMaxMoving(3,2) ) SectionCapac(3,2) = SectionCapac(3,2) + We

       END DO

!------------------------------------------------------------------------------
     END FUNCTION ElementEnergy
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
     SUBROUTINE ElectricForceIntegrate( Force, Moment, Area, Charge )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Force(3), Moment(3), Area, Charge
!------------------------------------------------------------------------------

       TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
       REAL(KIND=dp), POINTER :: U_Integ(:), V_Integ(:), W_Integ(:), S_Integ(:)
       REAL(KIND=dp) :: u,v,w,s, detJ, Tensor(3,3), Normal(3), EField(3), DFlux(3), Lforce(3), &
           LMoment(3), Radius(3), ElemForce(3), ElemMoment(3), ElemArea, ElemCharge, &
           ElementArea, ElementForce(3), xpos, ypos, zpos, LCharge
       INTEGER :: N_Integ,i,j,l
       LOGICAL :: stat
       
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints( CurrentElement )
       
       U_Integ => IntegStuff % u
       V_Integ => IntegStuff % v
       W_Integ => IntegStuff % w
       S_Integ => IntegStuff % s
       N_Integ =  IntegStuff % n

       ElemForce = 0.0d0
       ElemMoment = 0.0d0
       ElemArea = 0.0d0
       ElemCharge = 0.0d0

!------------------------------------------------------------------------------
!     Over integration points
!------------------------------------------------------------------------------
       DO l=1,N_Integ
 !------------------------------------------------------------------------------
         u = U_Integ(l)
         v = V_Integ(l)
         w = W_Integ(l)
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( CurrentElement, ElementNodes, u, v, w, &
             detJ, Basis, dBasisdx, ddBasisddx, .FALSE., .FALSE. )         
         s = detJ * S_Integ(l)

         Model % CurrentElement => CurrentElement         

!------------------------------------------------------------------------------
! The checking of normal direction is problematic since it requires that 
! densities are given and new coordinates are updated. For this reason
! the normal direction is set by the PlusMinus sign later on.
!------------------------------------------------------------------------------


         Normal = NormalVector( CurrentElement, ElementNodes, u, v, .FALSE. )
!------------------------------------------------------------------------------
! Need parent element basis etc., for computing normal derivatives on boundary.
!------------------------------------------------------------------------------
         DO i = 1,n
           DO j = 1,pn
             IF ( CurrentElement % NodeIndexes(i) == &
                 Parent % NodeIndexes(j) ) THEN
               xelem(i) = Parent % TYPE % NodeU(j)
               yelem(i) = Parent % TYPE % NodeV(j)
               zelem(i) = Parent % TYPE % NodeW(j)
               EXIT
             END IF
           END DO
         END DO
         
         u = SUM( Basis(1:n) * xelem(1:n) )
         v = SUM( Basis(1:n) * yelem(1:n) )
         w = SUM( Basis(1:n) * zelem(1:n) )
         
         stat = ElementInfo( Parent, ParentNodes, u, v, w, detJ, ParentBasis, &
             ParentdBasisdx, ddBasisddx, .FALSE., .FALSE. )

!------------------------------------------------------------------------------
         DO i=1,DIM
           EField(i) = SUM( ParentdBasisdx(1:pn,i) * ElemPotential(1:pn) )
         END DO
         DFlux = EField
         
         DO i=1, DIM
           DO j=1, DIM
             Tensor(i,j) = - DFlux(i) * EField(j)
           END DO
         END DO
         DO i=1, DIM
           Tensor(i,i) = Tensor(i,i) + SUM( DFlux * EField ) / 2.0d0
         END DO

         ElemArea = ElemArea + s 

         LCharge = SUM( DFlux(1:DIM) * Normal(1:DIM) )
         ElemCharge = ElemCharge + s * LCharge

         LForce = - MATMUL( Tensor, Normal )
         ElemForce  = ElemForce + s * LForce
        
         IF(CalculateMoment) THEN
           Radius(1) = SUM( ElementNodes % x(1:n) * Basis(1:n) ) - Center(1)
           Radius(2) = SUM( ElementNodes % y(1:n) * Basis(1:n) ) - Center(2)
           Radius(3) = SUM( ElementNodes % z(1:n) * Basis(1:n) ) - Center(3)
           
           LMoment(1) = Radius(2) * LForce(3) - Radius(3) * LForce(2)
           LMoment(2) = Radius(3) * LForce(1) - Radius(1) * LForce(3)
           LMoment(3) = Radius(1) * LForce(2) - Radius(2) * LForce(1)

           ElemMoment = ElemMoment + s * Lmoment
         END IF
         
!------------------------------------------------------------------------------
       END DO

       IF(ElemCharge > 0) THEN
         ElemCharge = -ElemCharge
         ElemForce = -ElemForce
         ElemMoment = -ElemMoment
       END IF
                
       Force = Force + ElemForce
       Area = Area + ElemArea
       Charge = Charge + ElemCharge
       Moment = Moment + ElemMoment

!------------------------------------------------------------------------------
     END SUBROUTINE ElectricForceIntegrate
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   END SUBROUTINE MovingElstatSolver
!------------------------------------------------------------------------------

