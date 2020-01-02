!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Jun 1997
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!----------------------------------------------------------------------
!> Module for solving Gebhardt factors inside the Elmer code
!----------------------------------------------------------------------


   MODULE RadiationFactorGlobals

      USE CRSMatrix
      USE IterSolve
 
      DOUBLE PRECISION, ALLOCATABLE :: GFactorFull(:,:)
      TYPE(Matrix_t), POINTER :: GFactorSP

   END MODULE RadiationFactorGlobals


   SUBROUTINE RadiationFactors( TSolver, FirstTime )

     USE DefUtils
     USE RadiationFactorGlobals

     IMPLICIT NONE

     LOGICAL :: FirstTime
     TYPE(Solver_t) :: TSolver

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Model_t), POINTER :: Model
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Solver_t), POINTER :: Solver
     TYPE(Element_t), POINTER :: Element
     TYPE(Matrix_t), POINTER :: Amatrix
     TYPE(Factors_t), POINTER :: GebhardtFactors, ViewFactors(:)

     INTEGER :: i,j,k,l,t,n,n2,istat,Colj,j0,couplings, sym 

     REAL (KIND=dp) :: SimulationTime,dt,MinFactor,MaxOmittedFactor, ConsideredSum
     REAL (KIND=dp) :: s,r,Norm,PrevNorm,Transmittivity, &
         maxds,refds,ds,x,y,z,&
         dx,dy,dz,x0(1),y0(1),MeshU(TSolver % Mesh % MaxElementNodes)
     REAL (KIND=dp), POINTER :: Vals(:), Wrk(:,:)
     REAL (KIND=dp), ALLOCATABLE :: SOL(:), RHS(:), Fac(:), RowSums(:), &
             Reflectivity(:), Emissivity(:), Areas(:), RelAreas(:), Diag(:)
     REAL (KIND=dp) :: MinSum, MaxSum, SolSum, PrevSelf, FactorSum, &
         ImplicitSum, ImplicitLimit, NeglectLimit, SteadyChange, FactorsFixedTol, &
         GeometryFixedTol, BackScale(3), Coord(3)
     REAL (KIND=dp) :: at, at0, st
     INTEGER :: BandSize,SubbandSize,RadiationSurfaces, GeometryFixedAfter, &
         Row,Col,MatrixElements, MatrixFormat, maxind, TimesVisited=0, &
         RadiationBody, MaxRadiationBody, MatrixEntries, ImplicitEntries, &
         FactorsFixedAfter
     INTEGER, POINTER :: Cols(:), NodeIndexes(:), TempPerm(:), NewPerm(:)
     INTEGER, ALLOCATABLE :: FacPerm(:)
     INTEGER, ALLOCATABLE, TARGET :: RowSpace(:),Reorder(:),ElementNumbers(:), &
         InvElementNumbers(:)

     CHARACTER(LEN=MAX_NAME_LEN) :: RadiationFlag, GebhardtFactorsFile, &
         ViewFactorsFile,OutputName, OutputName2, SolverType 
     CHARACTER(LEN=100) :: cmd

     LOGICAL :: GotIt, SaveFactors, UpdateViewFactors, UpdateGebhardtFactors, &
         ComputeViewFactors, OptimizeBW, TopologyTest, TopologyFixed, &
         FilesExist, FullMatrix, ImplicitLimitIs, IterSolveGebhardt, &
         ConstantEmissivity, Found, Debug, gSymm, gTriv, BinaryMode
     LOGICAL, POINTER :: ActiveNodes(:)
     LOGICAL :: DoScale
     INTEGER, PARAMETER :: VFUnit = 10

     TYPE(ValueList_t), POINTER :: Params, BC
     
     SAVE TimesVisited 

     EXTERNAL RMatvec

     Model => CurrentModel

     Found = .FALSE.
     DO i=1,Model % NumberOfBCs
       RadiationFlag = GetString( Model % BCs(i) % Values, 'Radiation', GotIt )
        IF (GotIt) THEN
          IF ( RadiationFlag == 'diffuse gray' ) Found = .TRUE.
        END IF
     END DO
     IF(.NOT. Found) RETURN

     Mesh => TSolver % Mesh
     CALL SetCurrentMesh( Model, Mesh )

     IF (.NOT. ASSOCIATED(Model)) THEN
       CALL Fatal('RadiationFactors','No pointer to model')
     END IF
     IF (.NOT. ASSOCIATED(Mesh) ) THEN
       CALL Fatal('RadiationFactors','No pointer to mesh')
     END IF

     IF(.NOT. FirstTime) TimesVisited = TimesVisited + 1

     Params => TSolver % Values

     UpdateViewFactors = GetLogical( Params, 'Update View Factors', GotIt )

     GeometryFixedAfter = GetInteger( Params, &
         'View Factors Fixed After Iterations',GotIt)
     IF(.NOT. GotIt) GeometryFixedAfter = HUGE(GeometryFixedAfter)

     IF( UpdateViewFactors ) THEN
       IF(GeometryFixedAfter < TimesVisited) UpdateViewFactors = .FALSE.
       IF(TimesVisited > 1 ) THEN
         SteadyChange = TSolver % Variable % SteadyChange
         GeometryFixedTol = GetConstReal( Params, 'View Factors Fixed Tolerance',GotIt)
         IF(GotIt .AND. SteadyChange < GeometryFixedTol) UpdateViewFactors = .FALSE.
       END IF
     END IF 

     UpdateGebhardtFactors = GetLogical( Params, 'Update Gebhardt Factors',GotIt )

     FactorsFixedAfter = GetInteger( Params, &
         'Gebhardt Factors Fixed After Iterations',GotIt)
     IF(.NOT. GotIt) FactorsFixedAfter = HUGE(FactorsFixedAfter)

     IF( UpdateGebhardtFactors ) THEN
       IF(FactorsFixedAfter < TimesVisited) UpdateGebhardtFactors = .FALSE.
       IF(TimesVisited > 1 ) THEN
         SteadyChange = TSolver % Variable % SteadyChange
         FactorsFixedTol = ListGetConstReal(TSolver % Values, &
           'Gebhardt Factors Fixed Tolerance',GotIt)
         IF( GotIt .AND. SteadyChange < FactorsFixedTol) UpdateGebhardtFactors = .FALSE.
       END IF
     END IF

     IF(.NOT. (FirstTime .OR. UpdateViewFactors .OR. UpdateGebhardtFactors)) THEN
       RETURN
     END IF

!------------------------------------------------------------------------------
!    Go for it
!------------------------------------------------------------------------------

     at0 = CPUTime()
     CALL Info('RadiationFactors','----------------------------------------------------',Level=5)
     CALL Info('RadiationFactors','Computing radiation factors for heat transfer',       Level=5)
     CALL Info('RadiationFactors','----------------------------------------------------',Level=5)

     ConstantEmissivity = FirstTime .AND. ( UpdateGebhardtFactors .OR. UpdateViewFactors )

     IF(GeometryFixedAfter == TimesVisited) THEN
       x0(1) = 1.0
       y0(1) = 1.0
     END IF

     FullMatrix = GetLogical( Params, 'Gebhardt Factors Solver Full',GotIt) 
     IF( FullMatrix ) THEN
       CALL Info('RadiationFactors','Using full matrix for Gebhardt factors',Level=6)
     ELSE
       CALL Info('RadiationFactors','Using sparse matrix for Gebhardt factors',Level=6)
     END IF

     IterSolveGebhardt =  GetLogical( Params, 'Gebhardt Factors Solver Iterative',GotIt) 
     IF(.NOT. GotIt) THEN
       SolverType = GetString( Params, 'radiation: Linear System Solver', GotIt )
       IF( GotIt ) THEN
         IF( SolverType == 'iterative' ) IterSolveGebhardt = .TRUE. 
       END IF
     END IF
     IF( IterSolveGebhardt ) THEN
       CALL Info('RadiationFactors','Using iterative solver for Gebhardt factors',Level=6)
     ELSE
       CALL Info('RadiationFactors','Using direct solver for Gebhardt factors',Level=6)
     END IF
       
     ComputeViewFactors = GetLogical( Params, 'Compute View Factors',GotIt )

!------------------------------------------------------------------------------
!    Compute the number of elements at the surface and check if the 
!    geometry has really changed.
!------------------------------------------------------------------------------

     RadiationSurfaces = 0
     MaxRadiationBody = 1
     ALLOCATE( ActiveNodes(Model % NumberOfNodes ), STAT=istat )
     IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 1.')
     ActiveNodes = .FALSE.

     DO t=Model % NumberOfBulkElements+1, &
         Model % NumberOfBulkElements + Model % NumberOfBoundaryElements
       
       Element => Model % Elements(t)
       k = Element % BoundaryInfo % Constraint
       
       IF ( Element % TYPE % ElementCode /= 101 ) THEN
         DO i=1,Model % NumberOfBCs
           IF ( Model % BCs(i) % Tag == k ) THEN
             RadiationFlag = GetString( Model % BCs(i) % Values, 'Radiation', GotIt )
             IF ( RadiationFlag == 'diffuse gray' ) THEN
               l = MAX(1, GetInteger( Model % BCs(i) % Values,'Radiation Boundary',GotIt) )
               MaxRadiationBody = MAX(l, MaxRadiationBody)

               NodeIndexes =>  Element % NodeIndexes
               n = Element % TYPE % NumberOfNodes

               ActiveNodes(NodeIndexes) = .TRUE.
               RadiationSurfaces = RadiationSurfaces + 1
                 
               IF(GeometryFixedAfter == TimesVisited) THEN
                 MeshU(1:n) = GetReal(Model % BCs(i) % Values, 'Mesh Update 1',GotIt, Element)
                 IF(.NOT. GotIt) THEN
                   WRITE (Message,'(A,I3)') 'Freezing Mesh Update 1 for bc',i
                   CALL Info('RadiationFactors',Message)
                   CALL ListAddDepReal( Model % BCs(i) % Values,'Mesh Update 1', &
                       'Mesh Update 1',1, x0, y0 )
                 END IF
                 MeshU(1:n) = GetReal(Model % BCs(i) % Values,'Mesh Update 2',GotIt, Element)
                 IF(.NOT. GotIt) THEN
                   WRITE (Message,'(A,I3)') 'Freezing Mesh Update 2 for bc',i
                   CALL Info('RadiationFactors',Message)
                   CALL ListAddDepReal( Model % BCs(i) % Values,'Mesh Update 2', &
                       'Mesh Update 2',1, x0, y0 )
                 END IF
               END IF
             END IF
           END IF
         END DO
       END IF
     END DO

     IF ( RadiationSurfaces == 0 ) THEN
       CALL Info('RadiationFactors','No surfaces participating in radiation',Level=5)
       DEALLOCATE(ActiveNodes)
       RETURN
     ELSE
       CALL Info('RadiationFactors','Total number of Radiation Surfaces '//TRIM(I2S(RadiationSurfaces))// &
           ' out of '//TRIM(I2S(Model % NumberOfBoundaryElements)),Level=5)
     END IF

     ! Check that the geometry has really changed before computing the viewfactors 
     IF(.NOT. FirstTime .AND. UpdateViewFactors) THEN

       ! This is a dirty thrick where the input file is stampered
       CALL Info('RadiationFactors','Checking changes in mesh.nodes file!',Level=5)
       
       OutputName = TRIM(OutputPath) // '/' // TRIM(Mesh % Name) // '/mesh.nodes.new'
       
       INQUIRE(FILE=TRIM(OutputName),EXIST=GotIt)
       IF(.NOT. GotIt) THEN
         OutputName = TRIM(OutputPath) // '/' // TRIM(Mesh % Name) // '/mesh.nodes'
       END IF

       OPEN( VFUnit,File=TRIM(OutputName) )
       dx = MAXVAL(Mesh % Nodes % x) - MINVAL(Mesh % Nodes % x)
       dy = MAXVAL(Mesh % Nodes % y) - MINVAL(Mesh % Nodes % y)
       dz = MAXVAL(Mesh % Nodes % z) - MINVAL(Mesh % Nodes % z)
       refds = SQRT(dx*dx+dy*dy+dz*dz)
       
       maxds = 0.0       
       maxind = 0
       GotIt = .FALSE.

       DO i=1,Mesh % NumberOfNodes
         READ(VFUnit,*,ERR=10,END=10) j,k,x,y,z
         IF(i == Mesh % NumberOfNodes) GotIt = .TRUE.
         IF(ActiveNodes(i)) THEN
           dx = Mesh % Nodes % x(i) - x
           dy = Mesh % Nodes % y(i) - y
           dz = Mesh % Nodes % z(i) - z
           ds = SQRT(dx*dx+dy*dy+dz*dz)
           IF(ds > maxds) THEN
             maxds = ds
             maxind = i
           END IF
         END IF
       END DO

10     CONTINUE

       CLOSE(VFUnit)

       IF(.NOT. GotIt) THEN
         UpdateViewFactors = .TRUE.        
         CALL Info('RadiationFactors','Mismatch in coordinates compared to file: '//TRIM(OutputName))
       ELSE
         WRITE(Message,'(A,E15.5)') 'Maximum geometry alteration on radiation BCs:',maxds
         CALL Info('RadiationFactors',Message)
         
         x = ListGetConstReal(TSolver % Values,'View Factors Geometry Tolerance',GotIt)
         IF(.NOT. GotIt) x = 1.0d-8
         
         IF(maxds <= refds * x) THEN
           CALL Info('RadiationFactors','Geometry change is neglected and old view factors are used')
           UpdateViewFactors = .FALSE.
         ELSE
           CALL Info('RadiationFactors','Geometry change requires recomputation of view factors')
         END IF
       END IF
     END IF
     
     DEALLOCATE(ActiveNodes)

     ! If the geometry has not changed and Gebhardt factors are fine return
     IF(.NOT. (FirstTime .OR. UpdateViewFactors .OR. UpdateGebhardtFactors)) THEN
       RETURN
     END IF     

!------------------------------------------------------------------------------
!    Check that the needed files exist if os assumed, if not recompute view factors
!------------------------------------------------------------------------------

     FilesExist = .TRUE.
     DO RadiationBody = 1, MaxRadiationBody 

       ViewFactorsFile = GetString(Model % Simulation,'View Factors',GotIt)
       
       IF ( .NOT.GotIt ) THEN
         ViewFactorsFile = 'ViewFactors.dat'
       END IF
       
       IF ( LEN_TRIM(Model % Mesh % Name) > 0 ) THEN
         OutputName = TRIM(OutputPath) // '/' // TRIM(Model % Mesh % Name) &
             // '/' // TRIM(ViewFactorsFile)
       ELSE
         OutputName = TRIM(ViewFactorsFile)
       END IF
       IF(RadiationBody > 1) THEN
         OutputName2 = OutputName
         WRITE(OutputName,'(A,I1)') TRIM(OutputName2),RadiationBody
       END IF
       INQUIRE(FILE=TRIM(OutputName),EXIST=GotIt)
       IF(.NOT. GotIt) FilesExist = .FALSE.
     END DO
     IF(.NOT. FilesExist) ComputeViewFactors = .TRUE.

!------------------------------------------------------------------------------
!    Rewrite the nodes for view factor computations if they have changed
!    and compute the ViewFactors with an external function call.
!------------------------------------------------------------------------------

     IF(ComputeViewFactors .OR. (.NOT. FirstTime .AND. UpdateViewFactors)) THEN
       
       ! This is a dirty thrick where the input mesh is scaled after loading.
       ! We need to perform scaling and backscaling then here too. 
       IF(.NOT. FirstTime) THEN
         CALL Info('RadiationFactors','Temporarely updating the mesh.nodes file!',Level=5)

         OutputName = TRIM(OutputPath) // '/' // TRIM(Mesh % Name) // '/mesh.nodes'         
         OutputName2 = TRIM(OutputPath) // '/' // TRIM(Mesh % Name) // '/mesh.nodes.orig'         
         CALL Rename(OutputName, OutputName2)         

         DoScale = ListCheckPresent( Model % Simulation,'Coordinate Scaling')

         IF( DoScale ) THEN
           Wrk => ListGetConstRealArray( Model % Simulation,'Coordinate Scaling',GotIt )    
           BackScale = 1.0_dp
           DO i=1,Mesh % MeshDim 
             j = MIN( i, SIZE(Wrk,1) )
             BackScale(i) = 1.0_dp / Wrk(j,1)
           END DO
         END IF

         OPEN( VFUnit,FILE=TRIM(OutputName), STATUS='unknown' )                 
         DO i=1,Mesh % NumberOfNodes
           Coord(1) = Mesh % Nodes % x(i)
           Coord(2) = Mesh % Nodes % y(i)
           Coord(3) = Mesh % Nodes % z(i)

           IF( DoScale ) Coord = BackScale * Coord
           
           WRITE( VFUnit,'(i7,i3,3f20.12)' ) i,-1, Coord
         END DO
         CLOSE(VFUnit)
       END IF
       
       ! Compute the factors using an external program call
       CALL SystemCommand( 'ViewFactors' )
       
       ! Set back the original node coordinates to prevent unwanted user errors
       IF(.NOT. FirstTime) THEN
         OutputName = TRIM(OutputPath) // '/' // TRIM(Mesh % Name) // '/mesh.nodes'         
         OutputName2 = TRIM(OutputPath) // '/' // TRIM(Mesh % Name) // '/mesh.nodes.new'         
         CALL Rename(OutputName, OutputName2)
         
         OutputName = TRIM(OutputPath) // '/' // TRIM(Mesh % Name) // '/mesh.nodes.orig'         
         OutputName2 = TRIM(OutputPath) // '/' // TRIM(Mesh % Name) // '/mesh.nodes'         
         CALL Rename(OutputName, OutputName2)
       END IF
     END IF


!------------------------------------------------------------------------------
!    Load the different view factors
!------------------------------------------------------------------------------

     TopologyFixed = GetLogical( Params, 'Matrix Topology Fixed',GotIt)
     MinFactor = GetConstReal( Params, 'Minimum Gebhardt Factor',GotIt )
     IF(.NOT. GotIt) MinFactor = 1.0d-20
     
     SaveFactors = ListGetLogical( Params, 'Save Gebhardt Factors',GotIt )
   
     MaxOmittedFactor = 0.0d0   
     TopologyTest = .TRUE.
     IF(FirstTime) TopologyTest = .FALSE. 

     RadiationBody = 0
     ALLOCATE( ElementNumbers(Model % NumberOfBoundaryElements), &
         RelAreas(Model % NumberOfBoundaryElements), &
         Areas(Model % NumberOfBoundaryElements), STAT=istat )
     IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 2.')

20   RadiationBody = RadiationBody + 1
     RadiationSurfaces = 0
     ElementNumbers = 0


     DO t=Model % NumberOfBulkElements+1, &
            Model % NumberOfBulkElements + Model % NumberOfBoundaryElements

       Element => Model % Elements(t)
       Model % CurrentElement => Element
       IF ( Element % TYPE % ElementCode == 101 ) CYCLE

       BC => GetBC(Element)
       IF(.NOT.ASSOCIATED(BC)) CYCLE

       RadiationFlag = GetString( BC, 'Radiation', GotIt )
       IF ( RadiationFlag == 'diffuse gray' ) THEN
         l = MAX(1, GetInteger( BC,'Radiation Boundary',GotIt) )
         n = GetElementNOFNodes(Element)
         IF(l == RadiationBody) THEN
           RadiationSurfaces = RadiationSurfaces + 1
           ElementNumbers(RadiationSurfaces) = t
           Areas(RadiationSurfaces) = ElementArea(Mesh,Element,n)
         END IF
       END IF
     END DO
     n = RadiationSurfaces
     RelAreas(1:n) = Areas(1:n) / MAXVAL(Areas(1:n))

     IF(MaxRadiationBody > 1) THEN
       CALL Info('RadiationFactors','Number of Radiation Surfaces '//TRIM(I2S(RadiationSurfaces))// &
           ' for boundary '//TRIM(I2S(RadiationBody)),Level=5)
       IF(RadiationSurfaces == 0) GOTO 30
     END IF

     ALLOCATE( RowSpace( RadiationSurfaces ), Reorder( RadiationSurfaces ), &
         InvElementNumbers(Model % NumberOfBoundaryElements), STAT=istat )
     IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 3.')
     RowSpace = 0
     ReOrder = 0

     ! Make the inverse of the list of element numbers of boundaries
     InvElementNumbers = 0
     j0 = Model % NumberOfBulkElements
     DO i=1,RadiationSurfaces
       InvElementNumbers(ElementNumbers(i)-j0) = i
     END DO

     ViewFactors => TSolver % Mesh % ViewFactors
     ! Open the file for ViewFactors

     
     IF( FirstTime .OR. UpdateViewFactors ) THEN

       IF ( .NOT.ASSOCIATED(ViewFactors) ) THEN
         ALLOCATE( ViewFactors(RadiationSurfaces), STAT=istat )
         IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 4.')
       ELSE
         IF (RadiationSurfaces /= SIZE(ViewFactors) ) THEN
            DO i=1,SIZE(ViewFactors)
               IF (ASSOCIATED(ViewFactors(i) % Factors )) DEALLOCATE(ViewFactors(i) % Factors )
               IF (ASSOCIATED(ViewFactors(i) % Elements)) DEALLOCATE(ViewFactors(i) % Elements)
            END DO
            DEALLOCATE( ViewFactors )

            ALLOCATE( ViewFactors(RadiationSurfaces), STAT=istat )
            IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 5.')
         END IF
       END IF
       TSOlver % Mesh % ViewFactors => ViewFactors

       ViewFactorsFile = GetString(Model % Simulation,'View Factors',GotIt)
       
       IF ( .NOT.GotIt ) THEN
         ViewFactorsFile = 'ViewFactors.dat'
       END IF
       
       IF ( LEN_TRIM(Model % Mesh % Name) > 0 ) THEN
         OutputName = TRIM(OutputPath) // '/' // TRIM(Model % Mesh % Name) // &
             '/' // TRIM(ViewFactorsFile)
       ELSE
         OutputName = TRIM(ViewFactorsFile)
       END IF
       IF(RadiationBody > 1) THEN
         OutputName2 = OutputName
         WRITE(OutputName,'(A,I1)') TRIM(OutputName2),RadiationBody
       END IF

       INQUIRE(FILE=TRIM(OutputName),EXIST=GotIt)
       IF(.NOT. GotIt) THEN
         WRITE(Message,*) 'View Factors File does NOT exist:',TRIM(OutputName)
         CALL Info('RadiationFactors','Message')
         GOTO 30
       END IF


       BinaryMode = ListGetLogical( Params,'Viewfactor Binary Output',Found ) 
         
       IF( BinaryMode ) THEN
         CALL Info('RadiationFactors','Loading view factors in binary mode',Level=5)

         OPEN( UNIT=VFUnit, FILE=TRIM(OutputName), FORM = 'unformatted', &
             ACCESS = 'stream', STATUS='old', ACTION='read' )         
         READ( VFUnit ) n
         IF( n /= RadiationSurfaces ) THEN
           CALL Fatal('RadiationFactors','Mismatch in viewfactor file size: '&
               //TRIM(I2S(n))//' vs. '//TRIM(I2S(RadiationSurfaces)))
         END IF
       ELSE
         OPEN( VFUnit,File=TRIM(OutputName) )
       END IF
                
       ! Read in the ViewFactors
       DO i=1,RadiationSurfaces
         IF( BinaryMode ) THEN
           READ( VFUnit ) n
         ELSE
           READ( VFUnit,* ) n
         END IF
           
         IF(FirstTime) THEN
           ViewFactors(i) % NumberOfFactors = n
           ALLOCATE( ViewFactors(i) % Elements(n) )
           ALLOCATE( ViewFactors(i) % Factors(n) )
         ELSE 
	   n2 = SIZE( ViewFactors(i) % Factors) 
	   IF(n /=  n2 ) THEN
             IF(ASSOCIATED(ViewFactors(i) % Elements)) DEALLOCATE( ViewFactors(i) % Elements )
	     IF(ASSOCIATED(ViewFactors(i) % Factors))  DEALLOCATE( ViewFactors(i) % Factors )
             ALLOCATE( ViewFactors(i) % Elements(n), ViewFactors(i) % Factors(n) )
           END IF
           ViewFactors(i) % NumberOfFactors = n
         END IF

         Vals => ViewFactors(i) % Factors
         Cols => ViewFactors(i) % Elements

	 Vals = 0; Cols = 0
  
         DO j=1,n
           IF( BinaryMode ) THEN
             READ(VFUnit) Cols(j),Vals(j)         
           ELSE
             READ(VFUnit,*) t,Cols(j),Vals(j)         
           END IF
           Vals(j) = RelAreas(i) * Vals(j)  ! Scale by area to make symmetric
           Cols(j) = ElementNumbers(Cols(j))
         END DO
       END DO
       CLOSE(VFUnit)
     END IF

     ! Check whether the element already sees itself, it will when 
     ! gebhardt factors are computed. Also compute matrix size.
     DO i=1,RadiationSurfaces
       Cols => ViewFactors(i) % Elements
       n = ViewFactors(i) % NumberOfFactors
       Rowspace(i) = n
       
       k = 0
       DO j=1,n
         IF ( Cols(j) == ElementNumbers(i) ) THEN
           k = 1
           EXIT
         END IF
       END DO
       IF ( k == 0 ) RowSpace(i) = RowSpace(i) + 1
     END DO
     MatrixElements = SUM(RowSpace(1:RadiationSurfaces))

     ALLOCATE( RHS(RadiationSurfaces), SOL(RadiationSurfaces),&
         Fac(RadiationSurfaces), FacPerm(RadiationSurfaces), STAT=istat )
     IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 6.')

     ! The coefficient matrix 
     CALL Info('RadiationFactors','Computing factors...',Level=5)

     IF(FullMatrix) THEN
       ALLOCATE(GFactorFull(RadiationSurfaces,RadiationSurfaces),STAT=istat)
       IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 7.')
       GFactorFull = 0.0_dp
     ELSE 
       ! Assembly the matrix form
       DO i=1,RadiationSurfaces
         Reorder(i) = i
       END DO

       GFactorSP => CRS_CreateMatrix( RadiationSurfaces, &
           MatrixElements,RowSpace,1,Reorder,.TRUE. )

       DO t=1,RadiationSurfaces
         Cols => ViewFactors(t) % Elements
         
         DO j=1,ViewFactors(t) % NumberOfFactors
           Colj = InvElementNumbers(Cols(j)-j0)
           CALL CRS_MakeMatrixIndex( GFactorSP,t,Colj )
         END DO

         CALL CRS_MakeMatrixIndex( GFactorSP,t,t )
       END DO

       CALL CRS_SortMatrix( GFactorSP )
       CALL CRS_ZeroMatrix( GFactorSP )
       MatrixEntries = SIZE(GFactorSP % Cols)
       WRITE(Message,'(A,T35,ES15.4)') 'View factors filling (%)',(100.0 * MatrixEntries) / &
           (RadiationSurfaces**2)
       CALL Info('RadiationFactors',Message,Level=5)
     END IF

     MatrixEntries = 0
     ImplicitEntries = 0

     ALLOCATE(Emissivity(RadiationSurfaces), Reflectivity(RadiationSurfaces))
     Emissivity = GetConstReal( Params,'Constant Emissivity',GotIt )
     IF(.NOT. GotIt) Emissivity = 0.5_dp
     Reflectivity = 1-Emissivity

     IF(.NOT. ConstantEmissivity) THEN
       DO i=1,RadiationSurfaces
         Element => Model % Elements(ElementNumbers(i))
         n = GetElementNOFNodes(Element)
         BC => GetBC(Element)
         Emissivity(i) = SUM(GetReal(BC, 'Emissivity', Found, Element))/n
         IF( Found ) THEN
           Transmittivity = SUM(GetReal(BC,'Transmissivity', Found, Element))/n
           Reflectivity(i) = 1 - Emissivity(i) - Transmittivity
         ELSE
            Emissivity(i) = SUM(GetParentMatProp( 'Emissivity', Element, Found))/n
            IF(.NOT. Found) THEN
              WRITE( Message,'(A,I3,I3)') 'Emissivity not found in BC: ',i,k
              CALL Fatal('RadiationFactors',Message)
            END IF
            Transmittivity = SUM(GetParentMatProp('Transmissivity',Element, Found))/n
            Reflectivity(i) = 1 - Emissivity(i) - Transmittivity
         END IF
       END DO
     END IF


     ALLOCATE(Diag(RadiationSurfaces)); Diag=0

     ! Fill the matrix for gebhardt factors
     ! Scale by (1-Emissivity) to get a symmetric system, if all emissivities are not equal to unity anywhere

     gTriv = ALL(ABS(Reflectivity)<=AEPS)
     gSymm = ALL(ABS(Reflectivity)>AEPS) .OR. gTriv

     IF(.NOT. gTriv) THEN
       r=1._dp
       DO i=1,RadiationSurfaces

         Element => Model % Elements(ElementNumbers(i))

         Vals => ViewFactors(i) % Factors
         Cols => ViewFactors(i) % Elements

         if ( gSymm ) r = Reflectivity(i)
         DO j=1,ViewFactors(i) % NumberOfFactors
           Colj = InvElementNumbers(Cols(j)-j0)

           IF (gSymm) THEN
             s = Reflectivity(i)*Reflectivity(colj)
           ELSE
             s = Reflectivity(i)
           END IF
           IF(FullMatrix) THEN
             GFactorFull(i,colj) = GFactorFull(i,colj)-s*Vals(j)
           ELSE
             CALL CRS_AddToMatrixElement(GFactorSP,i,colj,-s*Vals(j))
           END IF
         END DO

         Diag(i) = Diag(i) + r*RelAreas(i)
         IF(FullMatrix) THEN
           GFactorFull(i,i) = GFactorFull(i,i) + r*RelAreas(i)
         ELSE
           CALL CRS_AddToMatrixElement( GFactorSP,i,i,r*RelAreas(i) )
         END IF
       END DO

       ! Scale matrix to unit diagonals
       Diag = SQRT(1._dp/Diag)
       DO i=1,RadiationSurfaces
         IF(FullMatrix) THEN
           DO j=1,RadiationSurfaces
             GFactorFull(i,j) = GFactorFull(i,j)*Diag(i)*Diag(j)
           END DO
         ELSE
           DO j=GFactorSP % Rows(i),GFactorSP % Rows(i+1)-1
             GFactorSP % Values(j) = GFactorSP % Values(j)*Diag(i)*Diag(GFactorSP % Cols(j))
           END DO
         END IF
       END DO
     END IF

     ALLOCATE( Solver )
     CALL InitFactorSolver(TSolver, Solver)
     Params => TSolver % Values
     
     RHS = 0.0D0
     SOL = 1.0D-4

     ALLOCATE(RowSums(RadiationSurfaces))
     RowSums=0
          
     st = RealTime()

     n = 0
     DO t=1,RadiationSurfaces


       IF ( gTriv ) THEN

         Vals => ViewFactors(t) % Factors
         Cols => ViewFactors(t) % Elements
         Fac = 0._dp
         DO k=1,ViewFactors(t) % NumberOfFactors
           Colj = InvElementNumbers(Cols(k)-j0)
           Fac(Colj) = Vals(k) / RelAreas(t)
         END DO

         DO i=1,RadiationSurfaces
           s = Fac(i)*RelAreas(t)/(RelAreas(i)*Emissivity(i))
           IF (s>MinFactor ) n=n+1
           RowSums(i) = RowSums(i) + s
         END DO

       ELSE
         RHS(t) = Diag(t)

         ! It may be a good initial start that the Gii is 
         ! the same as previously
         IF (t>1) THEN
           PrevSelf = SOL(t)
           SOL(t) = SOL(t-1)
           SOL(t-1) = PrevSelf
           RHS(t-1) = 0.0_dp
         END IF

         SOL = SOL/Diag
         IF(FullMatrix) THEN
           CALL FIterSolver( RadiationSurfaces, SOL, RHS, Solver )
         ELSE
           IF(IterSolveGebhardt) THEN
             Solver % Matrix => GFactorSP
             CALL IterSolver( GFactorSP, SOL, RHS, Solver )
!------------------------------------------------------------------------------
           ELSE           
             IF (t==1) THEN
               CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', .TRUE. )
               CALL ListAddLogical( Solver % Values, 'Linear System Free Factorization', .FALSE. )
             ELSE IF(t==2) THEN
               CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', .FALSE. )
             END IF
  
             IF(.NOT.gSymm) THEN
               IF(GetString(Solver % Values,'Linear System Direct Method', Found)=='cholmod')THEN
                 CALL Warn('RadiationFactors', 'Can not use Cholesky solver if any emissivity==1')
                 CALL ListAddString( Solver % Values, 'Linear System Direct Method', 'UMFpack' )
               END IF
             END IF
             
             CALL DirectSolver( GFactorSP, SOL, RHS, Solver )
           END IF
         END IF
         SOL = SOL*Diag
         CALL ListRemove(Solver % Values,'Linear System Free Factorization')

         n = 0
         DO i=1,RadiationSurfaces
           Vals => ViewFactors(i) % Factors
           Cols => ViewFactors(i) % Elements

           s = 0.0_dp
           DO k=1,ViewFactors(i) % NumberOfFactors
             Colj = InvElementNumbers(Cols(k)-j0)
             IF(gSymm) THEN
               s = s + Reflectivity(colj)*Vals(k)*SOL(colj)
             ELSE
               s = s + Vals(k)*SOL(colj)
             END IF 
           END DO
           Fac(i) = s*Emissivity(t)*Emissivity(i)

           ! rowsums should add up to 1
           s = Fac(i)*RelAreas(t)/(RelAreas(i)*Emissivity(i))
           IF (s>MinFactor ) n=n+1
           RowSums(i) = RowSums(i) + s
         END  DO
       END IF

       FactorSum = SUM(Fac)
       ConsideredSum = 0.0_dp

       ImplicitLimit = GetConstReal( Params, 'Implicit Gebhardt Factor Fraction', ImplicitLimitIs) 
       NeglectLimit  = GetConstReal( Params, 'Neglected Gebhardt Factor Fraction', GotIt) 
       IF(.NOT. GotIt) NeglectLimit = 1.0d-6
       
       IF(ImplicitLimitIs) THEN
         DO i=1,RadiationSurfaces
           FacPerm(i) = i
         END DO
         ! Ensure that the self vision is always implicit to avoid trouble in the future!
         Fac(t) = Fac(t) + FactorSum
         CALL SortR( RadiationSurfaces, FacPerm, Fac) 
         Fac(1) = Fac(1) - FactorSum        

         ConsideredSum = 0.0_dp
         n = 0
         DO i=1,RadiationSurfaces
           IF(ConsideredSum < (1.0_dp-NeglectLimit) * FactorSum) THEN
             ConsideredSum = ConsideredSum + Fac(i)
             n = i
           END IF           
         END DO
       ELSE
         n = 0
         DO i=1,RadiationSurfaces
           IF ( Fac(i) > MinFactor ) n = n + 1
         END DO
       END IF

       MatrixEntries = MatrixEntries + n
       Element => Model % Elements(ElementNumbers(t))
       GebhardtFactors => Element % BoundaryInfo % GebhardtFactors       
       IF ( .NOT. ASSOCIATED( GebhardtFactors ) ) THEN
         ALLOCATE( GebhardtFactors )
         Element % BoundaryInfo % GebhardtFactors => GebhardtFactors
       END IF

       IF(FirstTime) THEN
         GebhardtFactors % NumberOfFactors = n
         GebhardtFactors % NumberOfImplicitFactors = n
         ALLOCATE( GebhardtFactors % Elements(n), GebhardtFactors % Factors(n) )
       ELSE IF(ImplicitLimitIs) THEN 
         IF( TopologyFixed ) THEN
           CALL Warn('RadiationFactors','Matrix topology cannot be fixed with implicit Gebhardt factors')
         END IF
         TopologyFixed = .FALSE.
         TopologyTest = .FALSE.
         DEALLOCATE( GebhardtFactors % Elements, GebhardtFactors % Factors )
         GebhardtFactors % NumberOfFactors = n
         ALLOCATE( GebhardtFactors % Elements(n), GebhardtFactors % Factors(n) )
         GebhardtFactors % NumberOfImplicitFactors = 0
       ELSE IF(GebhardtFactors % NumberOfFactors /= n .AND. .NOT. TopologyFixed) THEN         
         TopologyTest = .FALSE.
         DEALLOCATE( GebhardtFactors % Elements, GebhardtFactors % Factors )
         GebhardtFactors % NumberOfFactors = n
         ALLOCATE( GebhardtFactors % Elements(n), GebhardtFactors % Factors(n) )
         GebhardtFactors % NumberOfImplicitFactors = n
       END IF
       
       Vals => GebhardtFactors % Factors
       Cols => GebhardtFactors % Elements


       IF( ImplicitLimitIs ) THEN

         ImplicitSum = 0.0d0
         DO i=1,n
           Cols(i) = ElementNumbers(FacPerm(i)) 
           Vals(i) = Fac(i)

           IF(ImplicitSum < ImplicitLimit * FactorSum) THEN
             ImplicitSum = ImplicitSum + Fac(i)
             GebhardtFactors % NumberOfImplicitFactors = i
           END IF
         END DO
         Vals(2:n) = Vals(2:n) * (FactorSum - Vals(1)) / (ConsideredSum - Vals(1))

         IF(ImplicitLimit < TINY(ImplicitLimit)) GebhardtFactors % NumberOfImplicitFactors = 0       
         ImplicitEntries = ImplicitEntries + GebhardtFactors % NumberOfImplicitFactors

       ELSE IF(FirstTime .OR. .NOT. TopologyFixed) THEN
         n = 0
         DO i=1,RadiationSurfaces
           IF ( Fac(i) > MinFactor ) THEN
             n = n + 1
             IF(TopologyTest .AND. Cols(n) /= ElementNumbers(i)) TopologyTest = .FALSE.
             Cols(n) = ElementNumbers(i) 
             Vals(n) = Fac(i)
             ConsideredSum = ConsideredSum + Fac(i)
           END IF
         END DO
       ELSE
         ! If the topology is fixed the values are put only according to the existing structure
         ! and others are neglected
         n = GebhardtFactors % NumberOfFactors         
         Vals => GebhardtFactors % Factors
         Cols => GebhardtFactors % Elements
         
         DO i=1,n
           j = InvElementNumbers(Cols(i)-j0)
           Vals(i) = Fac(j)
           ConsideredSum = ConsideredSum + Fac(j)
         END DO
       END IF

       MaxOmittedFactor = MAX(MaxOmittedFactor,(FactorSum-ConsideredSum)/FactorSum) 
       
       IF ( RealTime() - st > 10.0 ) THEN
         WRITE(Message,'(a,i3,a)' ) '   Solution: ', &
         INT((100.0*t)/RadiationSurfaces),' % done'
         CALL Info( 'RadiationFactors', Message, Level=5 )
         st = RealTime()
       END IF
     END DO

     MinSum = MINVAL(RowSums)
     MaxSum = MAXVAL(RowSums)

     WRITE(Message,'(A,T35,2ES15.4)') 'Minimum Gebhardt factors sum',MINVAL(RowSums)
     CALL Info('RadiationFactors',Message,Level=5)
     WRITE(Message,'(A,T35,2ES15.4)') 'Maximum Gebhardt factors sum',MAXVAL(RowSums)
     CALL Info('RadiationFactors',Message,Level=5)
     WRITE(Message,'(A,T35,ES15.4)') 'Maximum share of omitted factors',MaxOmittedFactor
     CALL Info('RadiationFactors',Message,Level=5)
     WRITE(Message,'(A,T35,ES15.4)') 'Gebhardt factors filling (%)',(100.0 * MatrixEntries) / &
         (RadiationSurfaces**2)
     CALL Info('RadiationFactors',Message,Level=5)
     IF(ImplicitEntries > 0) THEN
       WRITE(Message,'(A,T35,ES15.4)') 'Implicit factors filling (%)',(100.0 * ImplicitEntries) / &
           (RadiationSurfaces**2)
       CALL Info('RadiationFactors',Message,Level=5)
     END IF


     ! Save factors is mainly for debugging purposes
     IF(SaveFactors) THEN
       GebhardtFactorsFile = GetString(Model % Simulation, 'Gebhardt Factors',GotIt )
       
       IF ( .NOT.GotIt ) THEN
         GebhardtFactorsFile = 'GebhardtFactors.dat'
       END IF
       
       IF ( LEN_TRIM(Model % Mesh % Name) > 0 ) THEN
         OutputName = TRIM(OutputPath) // '/' // TRIM(Model % Mesh % Name) // &
             '/' // TRIM(GebhardtFactorsFile)
       ELSE
         OutputName = TRIM(GebhardtFactorsFile) 
       END IF

       IF(RadiationBody > 1) THEN
         OutputName2 = OutputName
         WRITE(OutputName,'(A,I1)') TRIM(OutputName2), RadiationBody
       END IF

       OPEN( VFUnit,File=TRIM(OutputName) )
       
       WRITE (Message,'(A,A)') 'Writing Gephardt Factors to file: ',TRIM(OutputName)
       CALL Info('RadiationFactors',Message,Level=5)
       
       WRITE( VFUnit,* ) RadiationSurfaces
       
       DO t=1,RadiationSurfaces
         WRITE(VFUnit,*) t,ElementNumbers(t)
       END DO
       
       DO t=1,RadiationSurfaces
         Element => Model % Elements(ElementNumbers(t))
         GebhardtFactors => Element % BoundaryInfo % GebhardtFactors
         
         n = GebhardtFactors % NumberOfFactors 
         Vals => GebhardtFactors % Factors
         Cols => GebhardtFactors % Elements
         
         WRITE( VFUnit,* ) n
         DO i=1,n
           WRITE(VFUnit,*) t,InvElementNumbers(Cols(i)-j0),Vals(i)
         END DO
       END DO
       
       CLOSE(VFUnit)
     END IF

     DEALLOCATE(Solver, InvElementNumbers,RowSpace,Reorder,RHS,SOL,Fac,FacPerm,RowSums, &
                   Emissivity, Reflectivity,Diag )
     IF(FullMatrix) THEN
       DEALLOCATE(GFactorFull)
     ELSE
       CALL FreeMatrix( GFactorSP )     
     END IF

     WRITE (Message,'(A,T35,ES15.4)') 'Gebhardt factors determined (s)',CPUTime()-at0
     CALL Info('RadiationFactors',Message)


30   IF(RadiationBody < MaxRadiationBody) GOTO 20
     
     DEALLOCATE(ElementNumbers) 


     IF(.NOT. (FirstTime .OR. TopologyTest .OR. TopologyFixed) ) THEN

       CALL Info('RadiationFactors','Reorganizing the matrix structure',Level=5)

       MatrixFormat =  Tsolver % Matrix % FORMAT
       OptimizeBW = ListGetLogical(TSolver % Values,'Optimize Bandwidth',GotIt) 
       IF(.NOT. GotIt) OptimizeBW = .TRUE.

       CALL FreeMatrix( TSolver % Matrix)

       CALL Info('RadiationFactors','Creating new matrix topology')

       IF ( OptimizeBW ) THEN
         ALLOCATE( NewPerm( SIZE(Tsolver % Variable % Perm)) )
         TempPerm => Tsolver % Variable % Perm         
       ELSE
         NewPerm => Tsolver % Variable % Perm
       END IF

       AMatrix => CreateMatrix( CurrentModel,TSolver,TSolver % Mesh, &
           NewPerm, 1, MatrixFormat, OptimizeBW,  &
            ListGetString( TSolver % Values, 'Equation', Found ) )       
      
       ! Reorder the primary variable for bandwidth optimization:
       ! --------------------------------------------------------
       IF ( OptimizeBW ) THEN
         WHERE( NewPerm > 0 )
           TSolver % Variable % Values( NewPerm ) = &
               TSolver % Variable % Values( TempPerm )
         END WHERE
         
         IF ( ASSOCIATED( TSolver % Variable % PrevValues ) ) THEN
           DO j=1,SIZE( TSolver % Variable % PrevValues,2 )
             WHERE( NewPerm > 0 )
               TSolver % Variable % PrevValues( NewPerm,j) = &
                   TSolver % Variable % PrevValues(TempPerm,j)
             END WHERE
           END DO
         END IF

         Tsolver % Variable % Perm = NewPerm
         DEALLOCATE( NewPerm )
       END IF

       ! TODO: CreateMatrix should do these:
       ! -----------------------------------
       AMatrix % Lumped = GetLogical( Params, 'Lumped Mass Matrix', GotIt )
       AMatrix % Symmetric = ListGetLogical( Params, 'Linear System Symmetric', GotIt )       
       
       n = AMatrix % NumberOFRows
       ALLOCATE( AMatrix % RHS(n) )

       ! Transient case additional allocations:
       ! --------------------------------------
       IF ( ListGetString( CurrentModel % Simulation,'Simulation Type' ) == 'transient' ) THEN
         ALLOCATE( Amatrix % Force(n, TSolver % TimeOrder+1) )
         Amatrix % Force = 0.0d0
       END IF

       TSolver % Matrix => Amatrix
       CALL ParallelInitMatrix( TSolver, AMatrix )
     END IF
 
     WRITE (Message,'(A,T35,ES15.4)') 'All done time (s)',CPUTime()-at0
     CALL Info('RadiationFactors',Message)
     CALL Info('RadiationFactors','----------------------------------------------------',Level=5)
     
     
   CONTAINS

#include "huti_fdefs.h"
     SUBROUTINE FIterSolver( N,x,b,SolverParam )
       USE huti_sfe
       IMPLICIT NONE
       
       TYPE(Solver_t) :: SolverParam
       REAL (KIND=dp), DIMENSION(:) CONTIG :: x,b
       INTEGER :: N

       REAL (KIND=dp) :: dpar(50)       
       INTEGER :: ipar(50),wsize
       REAL (KIND=dp), ALLOCATABLE :: work(:,:)
!------------------------------------------------------------------------------
       INTEGER(KIND=AddrInt) :: AddrFunc
       EXTERNAL :: AddrFunc
       LOGICAL :: Found
       INTEGER(KIND=addrInt) :: iterProc, mvProc, dProc=0
!------------------------------------------------------------------------------
       
       ipar = 0
       dpar = 0.0_dp
       
       HUTI_WRKDIM = HUTI_CGS_WORKSIZE
       mvProc    = AddrFunc(rMatvec)

       iterProc  = AddrFunc(HUTI_D_CGS)
       wsize = HUTI_WRKDIM
       
       HUTI_NDIM     = N
       HUTI_DBUGLVL  = GetInteger( SOlverParam % Values, &
                'Linear System Residual Output', Found)
       HUTI_MAXIT    = GetInteger( SolverParam % Values, &
                'Linear System Max Iterations', Found )
       IF(.NOT.Found) HUTI_MAXIT = 100
       
       ALLOCATE( work(wsize,n) )
       
       IF ( ALL(x == 0.0) ) THEN
         HUTI_INITIALX = HUTI_RANDOMX
       ELSE
         HUTI_INITIALX = HUTI_USERSUPPLIEDX
       END IF
       
       HUTI_TOLERANCE = GetCReal( SolverParam % Values, &
               'Linear System Convergence Tolerance', Found)
       IF(.NOT.Found) HUTI_TOLERANCE = 1.d-10

       HUTI_MAXTOLERANCE = 1d20
       
       CALL IterCall( iterProc,x,b,ipar,dpar,work,mvProc, &
             dProc, dProc, dProc, dProc, dProc )

       DEALLOCATE( work )
       
     END SUBROUTINE FIterSolver
     

     SUBROUTINE InitFactorSolver(TSolver, Solver)

       TYPE(Solver_t) :: TSolver
       TYPE(Solver_t) :: Solver

       Solver % Values => ListAllocate()

       CALL ListCopyPrefixedKeywords( TSolver % Values, Solver % Values, 'radiation:' )

       
       CALL ListAddNewString( Solver % Values, &
           'Linear System Iterative Method', 'CGS' )

       CALL ListAddNewString( Solver % Values, &
           'Linear System Direct Method', 'Umfpack' )

       CALL ListAddNewInteger( Solver % Values, &
           'Linear System Max Iterations', 500 )
       
       CALL ListAddNewConstReal( Solver % Values, &
           'Linear System Convergence Tolerance', 1.0D-9 )
       
       CALL ListAddNewString( Solver % Values, &
           'Linear System Preconditioning', 'None' )
       
       CALL ListAddNewInteger( Solver % Values, &
           'Linear System Residual Output', 10 )
       
     END SUBROUTINE InitFactorSolver

   END SUBROUTINE RadiationFactors


   SUBROUTINE RMatvec( u,v,ipar )
     
     USE RadiationFactorGlobals
     
     REAL (KIND=dp) :: u(*),v(*)
     INTEGER :: ipar(*)
     
     INTEGER :: i,j,n
     REAL (KIND=dp) :: s

     
     n = HUTI_NDIM
     CALL DGEMV('N',n,n,1._dp,GFactorFull,n,u,1,0._dp,v,1)
     
!    DO i=1,n
!      s = 0.0D0
!      DO j=1,n
!        s = s + GFactorFull(i,j) * u(j)
!      END DO
!      v(i) = s
!    END DO
     
   END SUBROUTINE RMatvec

!> \}   
   
