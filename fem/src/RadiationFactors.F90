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
!> Module for solving Gebhart factors inside the Elmer code
!----------------------------------------------------------------------


   MODULE RadiationFactorGlobals
      USE Types
      TYPE(Matrix_t), POINTER :: G
      REAL(KIND=dp), ALLOCATABLE :: Gfull(:,:)
   END MODULE RadiationFactorGlobals


   SUBROUTINE RadiationFactors( TSolver, FirstTime, Newton )

     USE DefUtils
     USE RadiationFactorGlobals

     IMPLICIT NONE

     LOGICAL :: FirstTime
     LOGICAL :: Newton
     TYPE(Solver_t) :: TSolver

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Model_t), POINTER :: Model
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Element_t), POINTER :: Element
     TYPE(Solver_t), POINTER, SAVE :: Solver => Null()
     TYPE(Factors_t), POINTER :: GebhartFactors, ViewFactors(:), RadiosityFactors
          
     INTEGER :: i,j,k,l,t,n,istat,j0
     REAL (KIND=dp) :: Transmittivity

     REAL (KIND=dp), ALLOCATABLE :: SOL(:), SOL_d(:), RHS(:), RHS_d(:), Fac(:), &
        Reflectivity(:), Emissivity(:), Areas(:), RelAreas(:), Diag(:), SurfaceTemperature(:)
 
     REAL(KIND=dp), POINTER :: Temperature(:)

     REAL (KIND=dp) :: SteadyChange, Tol, Sigma

     REAL (KIND=dp) :: at, at0, st
     INTEGER :: RadiationSurfaces, GeometryFixedAfter, &
         TimesVisited=0, RadiationBody, MaxRadiationBody, FactorsFixedAfter

     INTEGER, ALLOCATABLE :: FacPerm(:)
     INTEGER, POINTER :: TempPerm(:)
     INTEGER, ALLOCATABLE, TARGET :: ElementNumbers(:), InvElementNumbers(:)

     CHARACTER(:), ALLOCATABLE :: RadiationFlag, SolverType

     LOGICAL :: GotIt, SaveFactors, UpdateViewFactors, UpdateGebhartFactors,  &
         ComputeViewFactors, TopologyTest, TopologyFixed, FullMatrix,         &
         IterSolveFactors, ConstantEmissivity, Found,  UpdateRadiatorFactors, &
         ComputeRadiatorFactors, RadiatorsFound, DiffuseGrayRadiationFound,   &
         UpdateGeometry, Radiosity, Spectral, Failed

     INTEGER, PARAMETER :: VFUnit = 10
     LOGICAL, ALLOCATABLE :: ActiveNodes(:)

     TYPE(Variable_t), POINTER :: Var
     TYPE(ValueList_t), POINTER :: Params, BC
     
     SAVE TimesVisited 

     Model  => CurrentModel
     Params => TSolver % Values

     Found = .FALSE.
     RadiatorsFound = .FALSE.
     DiffuseGrayRadiationFound = .FALSE.
     
     DO i=1,Model % NumberOfBCs
       BC => Model % BCs(i) % Values
       RadiatorsFound = RadiatorsFound .OR. GetLogical( BC, 'Radiator BC', Gotit )
       RadiationFlag = GetString( BC, 'Radiation', GotIt )
       IF ( RadiationFlag == 'diffuse gray' ) DiffuseGrayRadiationFound = .TRUE.
     END DO
     IF(.NOT. DiffuseGrayRadiationFound .AND. .NOT. RadiatorsFound) RETURN

     Radiosity = GetLogical( Params, 'Radiosity Model', Found )
     Spectral = GetLogical( Params, 'Spectral Model' ,Found )
     IF( Spectral ) Radiosity = .TRUE.

     CALL InitFactorSolver(TSolver,Solver)
     Mesh => TSolver % Mesh
     CALL SetCurrentMesh( Model, Mesh )

     IF (.NOT. ASSOCIATED(Model)) THEN
       CALL Fatal('RadiationFactors','No pointer to model')
     END IF
     IF (.NOT. ASSOCIATED(Mesh) ) THEN
       CALL Fatal('RadiationFactors','No pointer to mesh')
     END IF

     IF(.NOT. FirstTime) TimesVisited = TimesVisited + 1

     UpdateViewFactors = GetLogical( Params, 'Update View Factors', GotIt )

     GeometryFixedAfter = GetInteger( Params, &
         'View Factors Fixed After Iterations',GotIt)
     IF(.NOT. GotIt) GeometryFixedAfter = HUGE(GeometryFixedAfter)

     IF( UpdateViewFactors ) THEN
       IF(GeometryFixedAfter < TimesVisited) UpdateViewFactors = .FALSE.
       IF(TimesVisited > 1 ) THEN
         SteadyChange = TSolver % Variable % SteadyChange
         Tol = GetConstReal( Params, 'View Factors Fixed Tolerance',GotIt)
         IF(GotIt .AND. SteadyChange < Tol) UpdateViewFactors = .FALSE.
       END IF
     END IF 

     UpdateGebhartFactors = GetLogical( Params, 'Update Gebhart Factors',GotIt )
     IF(.NOT.GotIt ) &
       UpdateGebhartFactors = GetLogical( Params, 'Update Gebhardt Factors',GotIt )

     IF( UpdateGebhartFactors ) THEN       
       FactorsFixedAfter = GetInteger( Params, &
         'Gebhart Factors Fixed After Iterations',GotIt)

       IF(.NOT.GotIt) &
         FactorsFixedAfter = GetInteger( Params, &
           'Gebhardt Factors Fixed After Iterations',GotIt)

       IF( GotIt ) THEN       
         IF(FactorsFixedAfter < TimesVisited) UpdateGebhartFactors = .FALSE.
       END IF
       
       FactorsFixedAfter = GetInteger( Params, &
         'Gebhart Factors Fixed After Nonlinear Iterations',GotIt)
       IF (.NOT. GotIt) &
          FactorsFixedAfter = GetInteger( Params, &
             'Gebhardt Factors Fixed After Nonlinear Iterations',GotIt)

       IF( GotIt ) THEN                
         Var => VariableGet( Mesh % Variables, 'nonlin iter' )
         IF( ASSOCIATED( Var ) ) THEN
           k = NINT( Var % Values(1) ) 
           IF(FactorsFixedAfter < k ) UpdateGebhartFactors = .FALSE.
         END IF
       END IF
      
       Tol = ListGetConstReal(TSolver % Values, &
           'Gebhart Factors Fixed Tolerance',GotIt)
       IF (.NOT. GotIt) &
         Tol = ListGetConstReal(TSolver % Values, &
             'Gebhardt Factors Fixed Tolerance',GotIt)

       IF( GotIt ) THEN
         IF(TimesVisited > 1 ) THEN
           SteadyChange = TSolver % Variable % SteadyChange
           IF( SteadyChange < Tol) UpdateGebhartFactors = .FALSE.
         END IF
       END IF
     END IF

     UpdateRadiatorFactors = GetLogical(Params,'Update Radiator Factors',GotIt)

     IF(.NOT. (FirstTime .OR. UpdateViewFactors .OR. UpdateGebhartFactors .OR. &
             UpdateRadiatorFactors .OR. Radiosity)) RETURN

!------------------------------------------------------------------------------
!    Go for it
!------------------------------------------------------------------------------

     at0 = CPUTime()
     CALL Info('RadiationFactors','----------------------------------------------------',Level=5)
     CALL Info('RadiationFactors','Computing radiation factors for heat transfer',       Level=5)
     CALL Info('RadiationFactors','----------------------------------------------------',Level=5)

     ConstantEmissivity = FirstTime .AND. &
            ( UpdateGebhartFactors .OR. UpdateViewFactors .OR. UpdateRadiatorFactors )
     
     FullMatrix = GetLogical( Params, 'Radiation Factors Solver Full',GotIt) 
     IF(.NOT.GotIt) &
       FullMatrix = GetLogical( Params, 'Gebhart Factors Solver Full',GotIt) 
     IF(.NOT.GotIt) &
       FullMatrix = GetLogical( Params, 'Gebhardt Factors Solver Full',GotIt) 
     IF( FullMatrix ) THEN
       CALL Fatal('RadiationFactors', &
             'Using full matrix format for radiation problems not available anymore.')
     ELSE
       CALL Info('RadiationFactors','Using sparse matrix format for factor computations.',Level=6)
     END IF

     IterSolveFactors = GetLogical( Params, 'Radiation Factors Solver Iterative',GotIt) 
     IF(.NOT.GotIt) &
       IterSolveFactors  =  GetLogical( Params, 'Gebhart Factors Solver Iterative',GotIt) 
     IF(.NOT.GotIt) &
        IterSolveFactors =  GetLogical( Params, 'Gebhardt Factors Solver Iterative',GotIt) 
     IF(.NOT. GotIt) THEN
       SolverType = GetString( Params, 'radiation: Linear System Solver', GotIt )
       IF( GotIt ) THEN
         IF( SolverType == 'iterative' ) IterSolveFactors = .TRUE. 
       END IF
     END IF
     IF( IterSolveFactors ) THEN
       CALL Info('RadiationFactors','Using iterative solver for radiation factors',Level=6)
     ELSE
       CALL Info('RadiationFactors','Using direct solver for radiation factors',Level=6)
     END IF
       
     ComputeViewFactors = GetLogical( Params, 'Compute View Factors',GotIt )
     ComputeRadiatorFactors = GetLogical( Params, 'Compute Radiator Factors',GotIt )
!------------------------------------------------------------------------------
!    Compute the number of elements at the surface and check if the 
!    geometry has really changed.
!------------------------------------------------------------------------------

     RadiationSurfaces = 0
     MaxRadiationBody  = 1

     ALLOCATE(ActiveNodes(Mesh % NumberOfNodes), STAT=istat )
     IF (istat/= 0) CALL Fatal('RadiationFactors','Memory allocation error 1.')
     ActiveNodes = .FALSE.

     DO t=1,GetNOFBoundaryElements()
       Element => GetBoundaryElement(t)
       IF ( GetElementFamily(Element)<=1 ) CYCLE

       BC => GetBC(Element)
       i = GetBCId(Element)
       IF(i<=0) CYCLE

       RadiationFlag = GetString( BC, 'Radiation', GotIt )
       IF ( RadiationFlag == 'diffuse gray' .OR. GetLogical(BC, 'Radiator BC', Found) ) THEN
         l = MAX(1, GetInteger( BC,'Radiation Boundary',GotIt) )
         MaxRadiationBody = MAX(l, MaxRadiationBody)

         n = GetElementNOFNodes(Element)
         ActiveNodes(Element % NodeIndexes) = .TRUE.

         RadiationSurfaces = RadiationSurfaces + 1
         IF(GeometryFixedAfter == TimesVisited) CALL FixGeometryAfter(n,Element,BC,i)
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
     IF(.NOT. FirstTime .AND. (UpdateViewFactors .OR. UpdateRadiatorFactors)) THEN
       IF( .NOT. CheckMeshHasChanged() ) THEN
         UpdateViewFactors = .FALSE.
         UpdateRadiatorFactors = .FALSE.
       END IF         
     END IF
     
     DEALLOCATE(ActiveNodes)

     ! If the geometry has not changed and Gebhart factors are fine return
     IF(.NOT. (FirstTime .OR. UpdateViewFactors .OR. UpdateGebhartFactors .OR. &
              UpdateRadiatorFactors .OR. Radiosity)) RETURN

!-----------------------------------------------------------------------------------
!    Check that the needed files exist if os assumed, if not recompute view factors
!-----------------------------------------------------------------------------------
     CALL CheckFactorsFilesExist()
     
!------------------------------------------------------------------------------
!    Rewrite the nodes for view factor computations if they have changed
!    and compute the ViewFactors with an external function call.
!------------------------------------------------------------------------------

     UpdateGeometry = ListGetLogical(Params,'Update Factors Geometry',Found )  
     IF(.NOT. Found ) THEN 
       UpdateGeometry = ComputeViewFactors .OR. (ComputeRadiatorFactors.AND.RadiatorsFound) .OR. &
           (.NOT. FirstTime .AND. (UpdateViewFactors .OR. UpdateRadiatorFactors))
     END IF
     CALL ComputeViewFactorsAndRadiators()

!------------------------------------------------------------------------------
!    Load the different view factors
!------------------------------------------------------------------------------

     TopologyFixed = GetLogical( Params, 'Matrix Topology Fixed',GotIt)
     SaveFactors = ListGetLogical( Params, 'Save Gebhart Factors',GotIt )
     IF(.NOT. GotIt) &
       SaveFactors = ListGetLogical( Params, 'Save Gebhardt Factors',GotIt )
   
     TopologyTest = .TRUE.
     IF(FirstTime) TopologyTest = .FALSE. 

     RadiationBody = 0
     n = Mesh % NumberOfBoundaryElements
     ALLOCATE(ElementNumbers(n),InvElementNumbers(n),RelAreas(n),Areas(n),STAT=istat)
     IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 2.')

20   RadiationBody = RadiationBody + 1
     ElementNumbers = 0
     RadiationSurfaces = 0

     j0 = Mesh % NumberOfBulkElements
     DO t=1,GetNOFBoundaryElements()
       Element => GetBoundaryElement(t)
       IF ( GetElementFamily(Element)<=1 ) CYCLE

       BC => GetBC(Element)
       IF(.NOT.ASSOCIATED(BC)) CYCLE

       RadiationFlag = GetString( BC, 'Radiation', GotIt )
       IF ( RadiationFlag == 'diffuse gray' .OR. GetLogical( BC, 'Radiator BC', Found) ) THEN
         l = MAX(1, GetInteger( BC,'Radiation Boundary',GotIt) )
         n = GetElementNOFNodes(Element)
         IF(l == RadiationBody) THEN
           RadiationSurfaces = RadiationSurfaces + 1
           ElementNumbers(RadiationSurfaces) = t+j0
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


     ! Make the inverse of the list of element numbers of boundaries
     InvElementNumbers = 0
     DO i=1,RadiationSurfaces
       InvElementNumbers(ElementNumbers(i)-j0) = i
     END DO

!-----------------------------------------------

     IF( RadiatorsFound .AND. (FirstTime .OR. UpdateRadiatorFactors) ) THEN
       Failed = ReadRadiatorFactorsFromFile()
       IF(Failed) GOTO 30 
     END IF

!--------------------------------------------------

     IF( .NOT. DiffuseGrayRadiationFound ) RETURN

     ViewFactors => TSolver % Mesh % ViewFactors
     ! Open the file for ViewFactors

     IF( FirstTime .OR. UpdateViewFactors ) THEN
       Failed = ReadViewFactorsFromFile()
       IF(Failed) GOTO 30
     END IF

     ! This is just to enable quicker testing. It applies only when emissivity is one everywhere...
     IF( ListGetLogical( Params,'Use ViewFactors As Gebhart Factors',Found ) ) THEN
       CALL Warn('RadiationFactors','Used ViewFactors for RadiationFactors (assumes eps=1)')
       CALL UseViewFactorsAsGebhartFactors()
       IF(SaveFactors) CALL SaveGebhartFactors()       
       RETURN       
     END IF

     
     ALLOCATE( RHS(RadiationSurfaces), SOL(RadiationSurfaces), &
          Fac(RadiationSurfaces), FacPerm(RadiationSurfaces), STAT=istat )
     IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 8.')

     ! The coefficient matrix 
     CALL Info('RadiationFactors','Computing factors...',Level=5)

     CALL CreateRadiationMatrix(RadiationSurfaces)

     ALLOCATE(Emissivity(RadiationSurfaces), Reflectivity(RadiationSurfaces), STAT=istat)
     IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 10.')
     CALL TabulateEmissivity()

     ALLOCATE(Diag(RadiationSurfaces), STAT=istat)
     IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 11.')

     Params => TSolver % Values

     RHS  = 0.0_dp
     Diag = 0.0_dp
          
     IF( Radiosity ) THEN
       IF(FirstTime) RETURN
       CALL CalculateRadiosity()
     ELSE
       ! Fill the matrix for gebhardt factors
       CALL CalculateGebhartFactors()
       CALL FreeMatrix(G);G => Null()
     END IF
            
     DEALLOCATE(Diag,RHS,SOL,Fac,FacPerm,Emissivity,Reflectivity)

     IF( Radiosity ) THEN
       WRITE (Message,'(A,T35,ES15.4)') 'Radiosity vector determined (s)',CPUTime()-at0
     ELSE
       WRITE (Message,'(A,T35,ES15.4)') 'Gebhart factors determined (s)',CPUTime()-at0
     END IF
     CALL Info('RadiationFactors',Message)

30   IF(RadiationBody < MaxRadiationBody) GOTO 20
     
     IF(.NOT. (FirstTime .OR. TopologyTest .OR. TopologyFixed) ) THEN
       CALL UpdateMatrixTopologyWithFactors()
     END IF
 
     WRITE (Message,'(A,T35,ES15.4)') 'All done time (s)',CPUTime()-at0
     CALL Info('RadiationFactors',Message)
     CALL Info('RadiationFactors','----------------------------------------------------',Level=5)


   CONTAINS


     SUBROUTINE FixGeometryAfter(n,Element, BC, BCind)
       INTEGER :: n, BCind
       TYPE(ValueList_t), POINTER :: BC
       TYPE(Element_t), POINTER :: Element

       LOGICAL :: GotIt
       REAL(KIND=dp) :: x0(1), y0(1), MeshU(n)
          
       x0(1) = 1.0; y0(1) = 1.0
       MeshU(1:n) = GetReal(BC, 'Mesh Update 1',GotIt, Element)
       IF(.NOT. GotIt) THEN
         WRITE (Message,'(A,I3)') 'Freezing Mesh Update 1 for bc',BCind
         CALL Info('RadiationFactors',Message)
         CALL ListAddDepReal(BC,'Mesh Update 1', 'Mesh Update 1',1,x0,y0 )
       END IF
       MeshU(1:n) = GetReal(BC,'Mesh Update 2',GotIt, Element)
       IF(.NOT. GotIt) THEN
         WRITE (Message,'(A,I3)') 'Freezing Mesh Update 2 for bc',BCind
         CALL Info('RadiationFactors',Message)
         CALL ListAddDepReal( BC,'Mesh Update 2', 'Mesh Update 2',1,x0,y0 )
       END IF
       MeshU(1:n) = GetReal(BC,'Mesh Update 3',GotIt, Element)
       IF(.NOT. GotIt) THEN
         WRITE (Message,'(A,I3)') 'Freezing Mesh Update 3 for bc',BCind
         CALL Info('RadiationFactors',Message)
         CALL ListAddDepReal( BC, 'Mesh Update 3', 'Mesh Update 3',1,x0,y0 )
       END IF
     END SUBROUTINE FixGeometryAfter


     FUNCTION CheckMeshHasChanged() RESULT ( HasChanged )

       CHARACTER(:), ALLOCATABLE :: OutputName
       LOGICAL :: HasChanged,GotIt
       INTEGER :: i
       REAL(KIND=dp) :: ds,dx,dy,dz,maxds,refds,maxind,x,y,z

       HasChanged = .FALSE.
       
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
         HasChanged = .TRUE.
         CALL Info('RadiationFactors','Mismatch in coordinates compared to file: '//TRIM(OutputName))
       ELSE
         WRITE(Message,'(A,E15.5)') 'Maximum geometry alteration on radiation BCs:',maxds
         CALL Info('RadiationFactors',Message)

         x = ListGetConstReal(TSolver % Values,'View Factors Geometry Tolerance',GotIt)
         IF(.NOT. GotIt) x = 1.0d-8

         IF(maxds <= refds * x) THEN
           CALL Info('RadiationFactors','Geometry change is neglected and old view factors are used')
         ELSE
           HasChanged = .TRUE.
           CALL Info('RadiationFactors','Geometry change requires recomputation of view factors')
         END IF
       END IF
     END FUNCTION CheckMeshHasChanged


     SUBROUTINE ComputeViewFactorsAndRadiators()

       CHARACTER(:), ALLOCATABLE :: cmd, OutputName, OutputName2
       LOGICAL :: DoScale
       INTEGER :: i,j
       REAL(KIND=dp), POINTER :: Wrk(:,:)
       REAL(KIND=dp) :: BackScale(3), Coord(3)

       ! This is a dirty thrick where the input mesh is scaled after loading.
       ! We need to perform scaling and backscaling then here too. 
       IF( UpdateGeometry ) THEN
         CALL Info('RadiationFactors','Temporarily updating the mesh.nodes file!',Level=5)
         
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
       IF (ComputeViewFactors .OR.  .NOT.FirstTime .AND. UpdateViewFactors ) THEN
         cmd = 'ViewFactors '//TRIM(GetSifName())
         CALL SystemCommand( cmd )
       END IF

       IF( RadiatorsFound ) THEN
         IF (ComputeRadiatorFactors .OR. .NOT.FirstTime .AND. UpdateRadiatorFactors ) THEN
           cmd = 'Radiators '//TRIM(GetSifName())
           CALL SystemCommand( cmd )
         END IF
       END IF
     
       ! Set back the original node coordinates to prevent unwanted user errors
       IF( UpdateGeometry ) THEN
         OutputName = TRIM(OutputPath) // '/' // TRIM(Mesh % Name) // '/mesh.nodes'         
         OutputName2 = TRIM(OutputPath) // '/' // TRIM(Mesh % Name) // '/mesh.nodes.new'         
         CALL Rename(OutputName, OutputName2)

         OutputName = TRIM(OutputPath) // '/' // TRIM(Mesh % Name) // '/mesh.nodes.orig'         
         OutputName2 = TRIM(OutputPath) // '/' // TRIM(Mesh % Name) // '/mesh.nodes'         
         CALL Rename(OutputName, OutputName2)
       END IF
     END SUBROUTINE ComputeViewfactorsAndRadiators
     

     SUBROUTINE CheckFactorsFilesExist()
       LOGICAL :: FilesExist, GotIt
       INTEGER :: RadiationBody
       CHARACTER(:), ALLOCATABLE :: ViewFactorsFile, RadiatorFactorsFile, &
            OutputName

       IF( DiffuseGrayRadiationFound ) THEN 
         FilesExist = .TRUE.
         DO RadiationBody = 1, MaxRadiationBody 
           ViewFactorsFile = GetString(Model % Simulation,'View Factors',GotIt)
           IF ( .NOT.GotIt ) ViewFactorsFile = 'ViewFactors.dat'
       
           IF ( LEN_TRIM(Model % Mesh % Name) > 0 ) THEN
             OutputName = TRIM(OutputPath) // '/' // TRIM(Model % Mesh % Name) &
                     // '/' // TRIM(ViewFactorsFile)
           ELSE
             OutputName = TRIM(ViewFactorsFile)
           END IF
           IF(RadiationBody > 1) OutputName = OutputName//TRIM(I2S(RadiationBody))
           INQUIRE(FILE=TRIM(OutputName),EXIST=GotIt)
           IF(.NOT. GotIt) FilesExist = .FALSE.
         END DO
         IF(.NOT. FilesExist) ComputeViewFactors = .TRUE.
       END IF

       IF( RadiatorsFound ) THEN
         FilesExist = .TRUE.
         DO RadiationBody = 1, MaxRadiationBody 
           RadiatorFactorsFile = GetString(Model % Simulation,'Radiator Factors',GotIt)
           IF ( .NOT.GotIt ) RadiatorFactorsFile = 'RadiatorFactors.dat'
         
           IF ( LEN_TRIM(Model % Mesh % Name) > 0 ) THEN
             OutputName = TRIM(OutputPath) // '/' // TRIM(Model % Mesh % Name) &
                     // '/' // TRIM(RadiatorFactorsFile)
           ELSE
             OutputName = TRIM(RadiatorFactorsFile)
           END IF
           IF(RadiationBody > 1) OutputName = OutputName//TRIM(I2S(RadiationBody))
           INQUIRE(FILE=TRIM(OutputName),EXIST=GotIt)
           IF(.NOT. GotIt) FilesExist = .FALSE.
         END DO
         IF(.NOT. FilesExist) ComputeRadiatorFactors = .TRUE.
       END IF
     END SUBROUTINE CheckFactorsFilesExist


     ! This is an add-on for including point like radiators into the system.
     ! They act as heat sources with known total power that is distributed among the
     ! surface elements that the point sees.
     !----------------------------------------------------------------------------------     
     FUNCTION ReadRadiatorFactorsFromFile() RESULT ( Failed ) 

       LOGICAL :: Failed 
       
       TYPE(BoundaryInfo_t), POINTER :: BoundaryInfo
       REAL(KIND=dp), ALLOCATABLE :: Vals(:)
       INTEGER, ALLOCATABLE ::  Cols(:)
       INTEGER :: NofRadiators, i,n
       LOGICAL :: BinaryMode, Found
       REAL(KIND=dp), POINTER :: Radiators(:,:)
       TYPE(ValueList_t), POINTER :: RadList
       CHARACTER(:), ALLOCATABLE :: RadiatorFactorsFile, OutputName

       CALL Info('RadiationFactors','Loading radiator factors!',Level=7)
       Failed = .FALSE.
       
       RadiatorFactorsFile = GetString(Model % Simulation,'Radiator Factors',GotIt)       
       IF ( .NOT.GotIt ) RadiatorFactorsFile = 'RadiatorFactors.dat'
       
       IF ( LEN_TRIM(Model % Mesh % Name) > 0 ) THEN
         OutputName = TRIM(OutputPath) // '/' // TRIM(Model % Mesh % Name) // &
                 '/' // TRIM(RadiatorFactorsFile)
       ELSE
            OutputName = TRIM(RadiatorFactorsFile)
       END IF

       IF(RadiationBody > 1) OutputName = OutputName//TRIM(I2S(RadiationBody))

       INQUIRE(FILE=TRIM(OutputName),EXIST=GotIt)
       IF(.NOT. GotIt) THEN
         WRITE(Message,*) 'Radiator Factors File does NOT exist:',TRIM(OutputName)
         CALL Info('RadiationFactors','Message')
         Failed = .TRUE.
         RETURN
       END IF

       BinaryMode = ListGetLogical( Params,'Radiatorfactor Binary Output',Found ) 
       IF(.NOT. Found) BinaryMode = ListGetLogical( Params,'Viewfactor Binary Output',Found ) 
         
       IF( BinaryMode ) THEN
         CALL Info('RadiationFactors','Loading radiator factors from binary file: '//TRIM(OutputName),Level=5)
         OPEN( UNIT=VFUnit, FILE=TRIM(OutputName), FORM = 'unformatted', &
             ACCESS = 'stream', STATUS='old', ACTION='read' )         
         READ( VFUnit ) n
         IF( n /= RadiationSurfaces ) THEN
           CALL Fatal('RadiationFactors','Mismatch in radiation factor file size: '&
               //TRIM(I2S(n))//' vs. '//TRIM(I2S(RadiationSurfaces)))
         END IF
       ELSE
         CALL Info('RadiationFactors','Loading radiator factors from ascii file: '//TRIM(OutputName),Level=5)
         OPEN( VFUnit,File=TRIM(OutputName) )
       END IF

       
       IF( .NOT. ListCheckPresentAnyBodyForce( Model,'Radiator Coordinates',RadList ) ) &
           RadList => Params
       
       CALL GetConstRealArray( RadList, Radiators, 'Radiator Coordinates', Found )
       IF(.NOT. Found ) CALL Fatal( 'RadiationFactors', 'No radiators present, quitting' )

       NofRadiators = SIZE(Radiators,1)

       ! Read in the RadiatorFactors
       DO i=1,NofRadiators
         IF( BinaryMode ) THEN
           READ( VFUnit ) n
         ELSE
           READ( VFUnit,* ) n
         END IF

         ALLOCATE( Vals(n), Cols(n) )
         Vals = 0; Cols = 0

         DO j=1,n
           IF( BinaryMode ) THEN
             READ(VFUnit) Cols(j),Vals(j)         
           ELSE
             READ(VFUnit,*) t,Cols(j),Vals(j)         
           END IF
           Vals(j) = Vals(j) / Areas(Cols(j))
           Cols(j) = ElementNumbers(Cols(j))
         END DO

         DO j=1,n
           BoundaryInfo => Mesh % Elements(Cols(j)) % BoundaryInfo
           IF ( .NOT.ALLOCATED( BoundaryInfo % Radiators ) ) THEN
             ALLOCATE( BoundaryInfo % Radiators(NofRadiators) )
             BoundaryInfo % Radiators = 0
           END IF

           BoundaryInfo % Radiators(i) = Vals(j)
         END DO
         DEALLOCATE( Cols, Vals )
       END DO
       CLOSE(VFUnit)
              
     END FUNCTION ReadRadiatorFactorsFromFile


     FUNCTION ReadViewFactorsFromFile() RESULT ( Failed )

       LOGICAL :: Failed

       LOGICAL :: Found, BinaryMode
       INTEGER :: i,j,n,n2
       INTEGER, POINTER :: Cols(:)
       REAL(KIND=dp), POINTER :: Vals(:)
       CHARACTER(:), ALLOCATABLE :: ViewFactorsFile, OutputName

       CALL Info('RadiationFactors','Loading view factors!',Level=7)

       Failed = .FALSE.
       
       IF ( .NOT.ASSOCIATED(ViewFactors) ) THEN
         ALLOCATE( ViewFactors(RadiationSurfaces), STAT=istat )
         IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 4.')
       ELSE
         IF (RadiationSurfaces /= SIZE(ViewFactors)) THEN
           DO i=1,SIZE(ViewFactors)
             IF (ASSOCIATED(ViewFactors(i) % Factors)) DEALLOCATE(ViewFactors(i) % Factors)
             IF (ASSOCIATED(ViewFactors(i) % Elements)) DEALLOCATE(ViewFactors(i) % Elements)
           END DO
           DEALLOCATE( ViewFactors )
           ALLOCATE( ViewFactors(RadiationSurfaces), STAT=istat )
           IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 5.')
         END IF
       END IF

       TSolver % Mesh % ViewFactors => ViewFactors

       ViewFactorsFile = GetString(Model % Simulation,'View Factors',GotIt)
       IF ( .NOT.GotIt ) ViewFactorsFile = 'ViewFactors.dat'

       IF ( LEN_TRIM(Model % Mesh % Name) > 0 ) THEN
         OutputName = TRIM(OutputPath) // '/' // TRIM(Model % Mesh % Name) // &
             '/' // TRIM(ViewFactorsFile)
       ELSE
         OutputName = TRIM(ViewFactorsFile)
       END IF

       IF(RadiationBody > 1) OutputName = OutputName//TRIM(I2S(RadiationBody))

       INQUIRE(FILE=TRIM(OutputName),EXIST=GotIt)
       IF(.NOT. GotIt) THEN
         WRITE(Message,*) 'View Factors File does NOT exist:',TRIM(OutputName)
         CALL Info('RadiationFactors','Message')
         Failed = .TRUE.
         RETURN
       END IF

       BinaryMode = ListGetLogical( Params,'Viewfactor Binary Output',Found ) 

       IF( BinaryMode ) THEN
         CALL Info('RadiationFactors','Loading view factors from binary file: '//TRIM(OutputName),Level=5)

         OPEN( UNIT=VFUnit, FILE=TRIM(OutputName), FORM = 'unformatted', &
             ACCESS = 'stream', STATUS='old', ACTION='read' )         
         READ( VFUnit ) n
         IF( n /= RadiationSurfaces ) THEN
           CALL Fatal('RadiationFactors','Mismatch in viewfactor file size: '&
               //TRIM(I2S(n))//' vs. '//TRIM(I2S(RadiationSurfaces)))
           Failed = .TRUE.
           RETURN
         END IF
       ELSE
         CALL Info('RadiationFactors','Loading view factors from ascii file: '//TRIM(OutputName),Level=5)
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
           ALLOCATE( ViewFactors(i) % Elements(n), ViewFactors(i) % Factors(n), STAT=istat )
           IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 6.')
         ELSE 
	   n2 = SIZE( ViewFactors(i) % Factors) 
	   IF(n /=  n2 ) THEN
             IF(ASSOCIATED(ViewFactors(i) % Elements)) DEALLOCATE( ViewFactors(i) % Elements )
	     IF(ASSOCIATED(ViewFactors(i) % Factors))  DEALLOCATE( ViewFactors(i) % Factors )
             ALLOCATE( ViewFactors(i) % Elements(n), ViewFactors(i) % Factors(n), STAT=istat )
             IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 7.')
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
!          Cols(j) = ElementNumbers(Cols(j))
         END DO
       END DO
       CLOSE(VFUnit)
       
     END FUNCTION ReadViewFactorsFromFile



     ! This uses the view factors mainly to provide a simple test case where
     ! we assumes that emissivity is everywhere 1.
     !------------------------------------------------------------------------
     SUBROUTINE UseViewFactorsAsGebhartFactors()
       INTEGER :: i,j,n
       REAL(KIND=dp) :: s

       DO i=1,RadiationSurfaces
         j = ElementNumbers(i)
         Element => Model % Mesh % Elements(j)

         GebhartFactors => Element % BoundaryInfo % RadiationFactors
         IF ( .NOT. ASSOCIATED( GebhartFactors ) ) THEN
           ALLOCATE( GebhartFactors )
           Element % BoundaryInfo % RadiationFactors => GebhartFactors
         END IF

         n = ViewFactors(i) % NumberOfFactors
         GebhartFactors % NumberOfFactors = n
         GebhartFactors % NumberOfImplicitFactors = n
         ALLOCATE( GebhartFactors % Elements(n), GebhartFactors % Factors(n), STAT=istat)
         IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 20.')

         s = SUM(Viewfactors(i) % Factors)
         GebhartFactors % Factors = ViewFactors(i) % Factors/s
         GebhartFactors % Elements = ElementNumbers(Viewfactors(i) % Elements)
!        GebhartFactors % Elements = Viewfactors(i) % Elements
       END DO       
     END SUBROUTINE UseViewFactorsAsGebhartFactors


     SUBROUTINE CreateRadiationMatrix(n)
       INTEGER :: n

       INTEGER, ALLOCATABLE :: RowSpace(:),Reorder(:)
       INTEGER, POINTER :: Cols(:)
       INTEGER :: i,j,k,nf,t,colj,previ,MatrixEntries

       ALLOCATE(RowSpace(n), Reorder(n))

       ! Check whether the element already sees itself, it will when 
       ! gebhardt factors are computed. Also compute matrix size.
       DO i=1,n
         RowSpace(i) = ViewFactors(i) % NumberOfFactors
         Cols => ViewFactors(i) % Elements
         IF (ALL(Cols/=i)) RowSpace(i) = RowSpace(i)+1
       END DO
       MatrixEntries = SUM(RowSpace(1:n))

       IF(.NOT.ASSOCIATED(G) .OR. UpdateViewFactors) THEN
         IF(ASSOCIATED(G)) CALL FreeMatrix(G)

         ! Assembly the matrix form
         Reorder = [(i, i=1,n)]
         G => CRS_CreateMatrix(n,MatrixEntries,RowSpace,1,Reorder,.TRUE. )

         DO t=1,n
           Cols => ViewFactors(t) % Elements         
           previ = G % Rows(t)-1
           DO j=1,ViewFactors(t) % NumberOfFactors
             CALL CRS_MakeMatrixIndex(G,t,Cols(j),previ)
           END DO
           CALL CRS_MakeMatrixIndex(G,t,t)
         END DO

         CALL CRS_SortMatrix(G)
         CALL CRS_ZeroMatrix(G)
         MatrixEntries = SIZE(G % Cols)
         WRITE(Message,'(A,T35,ES15.4)') 'View factors filling (%)',(100.0*MatrixEntries)/(n**2)
         CALL Info('RadiationFactors',Message,Level=5)
       ELSE
         CALL CRS_ZeroMatrix(G)
       END IF
     END SUBROUTINE CreateRadiationMatrix
     

     SUBROUTINE TabulateSurfaceTemperatures()
       TYPE(Element_t), POINTER :: Element
       LOGICAL :: DG, Found
       INTEGER :: i,n
       INTEGER, POINTER :: Inds(:)
       INTEGER, TARGET :: DGInds(27)

       DG = ListGetLogical(Params, 'Discontinuous Galerkin',Found ) .OR. & 
            ListGetLogical(Params, 'DG Reduced Basis',Found ) 
 
       DO i=1,RadiationSurfaces
         Element => Mesh % Elements(ElementNumbers(i))
         n = GetElementNOFNodes(Element)
         IF( DG ) THEN
           CALL DgRadiationIndexes(Element,n,DGInds,.TRUE.)
           Inds => DGInds(1:n)
         ELSE
           Inds => Element % NodeIndexes(1:n)
         END IF
         SurfaceTemperature(i) = SUM(Temperature(TempPerm(Inds)))/n
       END DO
     END SUBROUTINE TabulateSurfaceTemperatures


     SUBROUTINE TabulateEmissivity()
       REAL(KIND=dp) :: Transmittivity
       LOGICAL ::GotIt, Found
       INTEGER :: i,n
       TYPE(ValueList_t), POINTER :: BC
       TYPE(Element_t), POINTER :: Element

       Emissivity = GetConstReal( Params,'Constant Emissivity',GotIt )
       IF(.NOT. GotIt) Emissivity = 0.5_dp

       Reflectivity = 1-Emissivity

       IF(.NOT. (ConstantEmissivity .OR. Spectral .OR. Radiosity .AND. FirstTime ) ) THEN
         DO i=1,RadiationSurfaces
           Element => Mesh % Elements(ElementNumbers(i))
           n = GetElementNOFNodes(Element)
           BC => GetBC(Element)
           Emissivity(i) = SUM(GetReal(BC, 'Emissivity', Found, Element))/n

           IF( Found ) THEN
             Transmittivity = SUM(GetReal(BC,'Transmissivity', Found, Element))/n
             Reflectivity(i) = 1 - Emissivity(i) - Transmittivity
           ELSE
              Emissivity(i) = SUM(GetParentMatProp( 'Emissivity', Element, Found))/n
              IF(.NOT. Found) THEN
                WRITE( Message,'(A,I3,I3)') 'Emissivity not found in BC: ',i,GetBCId(Element)
                CALL Fatal('RadiationFactors',Message)
              END IF
              Transmittivity = SUM(GetParentMatProp('Transmissivity',Element, Found))/n
              Reflectivity(i) = 1 - Emissivity(i) - Transmittivity
           END IF
         END DO
       END IF
     END SUBROUTINE TabulateEmissivity

     !------------------------------------------------------------------------------
     ! To save some time tabulate the spectral emissivity data for each temperature.
     !------------------------------------------------------------------------------
     SUBROUTINE TabulateSpectralEmissivity(Emissivity,Trad,Emissivity_d)
       REAL(KIND=dp) :: Trad
       REAL(KIND=dp) :: Emissivity(:)
       REAL(KIND=dp), OPTIONAL :: Emissivity_d(:)
       
       TYPE(ValueList_t), POINTER :: Vlist
       TYPE(Element_t), POINTER :: Element, Parent
       INTEGER :: i,k,ent_id
       
       DO i=1,RadiationSurfaces         
         Element => Mesh % Elements(ElementNumbers(i))
         
         DO ent_id=1,CurrentModel % NumberOfBCs
           IF ( Element % BoundaryInfo % Constraint == CurrentModel % BCs(ent_id) % Tag ) EXIT
         END DO
         IF ( ent_id > CurrentModel % NumberOfBCs ) CALL Fatal('TabulateSpectralEmissivity','Could not find BC!')
         Vlist => CurrentModel % BCs(ent_id) % Values
         
         IF( .NOT. ListCheckPresent(Vlist,'Emissivity') ) THEN
           DO k=1,2
             IF(k==1) THEN
               Parent => Element % BoundaryInfo % Left
             ELSE
               Parent => Element % BoundaryInfo % Right
             END IF
             IF(ASSOCIATED(Parent) ) THEN
               IF( Parent % BodyId > 0 .AND. Parent % BodyId <= CurrentModel % NumberOfBodies ) THEN
                 ent_id = ListGetInteger( CurrentModel % Bodies(Parent % BodyId) % Values,'Material',Found)
                 IF(Found) THEN
                   Vlist => CurrentModel % Materials(ent_id) % Values
                   IF(ListCheckPresent(Vlist,'Emissivity') ) EXIT
                 END IF
               END IF
             END IF
           END DO
         END IF

         IF(PRESENT(Emissivity_d))  THEN
           Emissivity(i) = ListGetFun( Vlist,'Emissivity',Trad,minv=0.0_dp,maxv=1.0_dp, &
                         dfdx=Emissivity_d(i), eps=0.1_dp)
         ELSE
           Emissivity(i) = ListGetFun( Vlist,'Emissivity',Trad,minv=0.0_dp,maxv=1.0_dp)
         END IF
       END DO
       
     END SUBROUTINE TabulateSpectralEmissivity
            

     ! Calculate Gebhart factors for radiation. These result to good convergence with the cost of
     ! computing many smaller linear systems.
     !--------------------------------------------------------------------------------------------
     SUBROUTINE CalculateGebhartFactors()

       REAL(KIND=dp) :: MinFactor, MaxOmittedFactor, ConsideredSum
       INTEGER :: Colj,i,j,t, ImplicitEntries, MatrixEntries
       LOGICAL :: gTriv, gSymm, ImplicitLimitIs
       REAL(KIND=dp) :: r,s,st,PrevSelf, MinSum,MaxSum,SolSum,FactorSum,ImplicitSum,&
           ImplicitLimit, NeglectLimit

       REAL(KIND=dp), ALLOCATABLE :: RowSums(:)
       INTEGER, POINTER :: Cols(:)
       REAL(KIND=dp), POINTER :: Vals(:)

       MaxOmittedFactor = 0._dp
       MatrixEntries = 0
       ImplicitEntries = 0
       
       MinFactor = GetConstReal( Params, 'Minimum Gebhart Factor',GotIt )
       IF (.NOT. GotIt) &
           MinFactor = GetConstReal( Params, 'Minimum Gebhardt Factor',GotIt )
       IF(.NOT. GotIt) MinFactor = 1.0d-20

       ! Scale by (1-Emissivity) to get a symmetric system, if all emissivities are not equal to unity anywhere       
       gTriv = ALL(ABS(Reflectivity)<=AEPS)
       gSymm = ALL(ABS(Reflectivity)>AEPS) .OR. gTriv

              
       IF(.NOT. gTriv) THEN
         r=1._dp
         DO i=1,RadiationSurfaces
           Vals => ViewFactors(i) % Factors
           Cols => ViewFactors(i) % Elements

           if ( gSymm ) r = Reflectivity(i)
           DO j=1,ViewFactors(i) % NumberOfFactors
             IF (gSymm) THEN
               s = Reflectivity(i)*Reflectivity(Cols(j))
             ELSE
               s = Reflectivity(i)
             END IF
             CALL CRS_AddToMatrixElement(G,i,Cols(j),-s*Vals(j))
           END DO
           Diag(i) = Diag(i) + r*RelAreas(i)
           CALL CRS_AddToMatrixElement( G,i,i,r*RelAreas(i) )
         END DO
       END IF
       
       ! Scale matrix to unit diagonals
       Diag = SQRT(1._dp/ABS(Diag))
       DO i=1,RadiationSurfaces
         DO j=G % Rows(i),G % Rows(i+1)-1
           G % Values(j) = G % Values(j)*Diag(i)*Diag(G % Cols(j))
         END DO
       END DO
       
       SOL = 1.0d-4
       st = RealTime()

       n = 0
       ALLOCATE(RowSums(RadiationSurfaces), STAT=istat)
       IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 12.')       
       RowSums=0

       DO t=1,RadiationSurfaces
         IF ( gTriv ) THEN

           Vals => ViewFactors(t) % Factors
           Cols => ViewFactors(t) % Elements
           Fac = 0._dp
           DO k=1,ViewFactors(t) % NumberOfFactors
             Fac(Cols(k)) = Vals(k) / RelAreas(t)
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
           IF(IterSolveFactors) THEN
             Solver % Matrix => G
             CALL IterSolver( G, SOL, RHS, Solver )
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

             CALL DirectSolver( G, SOL, RHS, Solver )
           END IF

           SOL = SOL*Diag
           CALL ListRemove(Solver % Values,'Linear System Free Factorization')

           n = 0
           DO i=1,RadiationSurfaces
             Vals => ViewFactors(i) % Factors
             Cols => ViewFactors(i) % Elements

             s = 0.0_dp
             DO k=1,ViewFactors(i) % NumberOfFactors
               IF(gSymm) THEN
                 s = s + Reflectivity(Cols(k))*Vals(k)*SOL(Cols(k))
               ELSE
                 s = s + Vals(k)*SOL(Cols(k))
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

         ImplicitLimit = GetConstReal( Params, 'Implicit Gebhart Factor Fraction', ImplicitLimitIs) 
         IF  (.NOT. ImplicitLimitIs) &
            ImplicitLimit = GetConstReal( Params, 'Implicit Gebhardt Factor Fraction', ImplicitLimitIs) 

         NeglectLimit  = GetConstReal( Params, 'Neglected Gebhart Factor Fraction', GotIt) 
         IF(.NOT.Gotit) &
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
         Element => Mesh % Elements(ElementNumbers(t))
         GebhartFactors => Element % BoundaryInfo % RadiationFactors
         IF ( .NOT. ASSOCIATED( GebhartFactors ) ) THEN
           ALLOCATE( GebhartFactors )
           Element % BoundaryInfo % RadiationFactors => GebhartFactors
         END IF

         IF(FirstTime) THEN
           GebhartFactors % NumberOfFactors = n
           GebhartFactors % NumberOfImplicitFactors = n
           ALLOCATE( GebhartFactors % Elements(n), GebhartFactors % Factors(n), STAT=istat)
           IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 13.')
         ELSE IF(ImplicitLimitIs) THEN 
           IF( TopologyFixed ) THEN
             CALL Warn('RadiationFactors','Matrix topology cannot be fixed with implicit Gebhart factors')
           END IF
           TopologyFixed = .FALSE.
           TopologyTest = .FALSE.
           DEALLOCATE( GebhartFactors % Elements, GebhartFactors % Factors )
           GebhartFactors % NumberOfFactors = n
           ALLOCATE( GebhartFactors % Elements(n), GebhartFactors % Factors(n), STAT=istat )
           IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 14.')
           GebhartFactors % NumberOfImplicitFactors = 0
         ELSE IF(GebhartFactors % NumberOfFactors /= n .AND. .NOT. TopologyFixed) THEN         
           TopologyTest = .FALSE.
           DEALLOCATE( GebhartFactors % Elements, GebhartFactors % Factors )
           GebhartFactors % NumberOfFactors = n         
           ALLOCATE( GebhartFactors % Elements(n), GebhartFactors % Factors(n), STAT=istat )
           IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 15.')

           GebhartFactors % NumberOfImplicitFactors = n
         END IF

         Vals => GebhartFactors % Factors
         Cols => GebhartFactors % Elements

         IF( ImplicitLimitIs ) THEN

           ImplicitSum = 0.0d0
           DO i=1,n
             Cols(i) = ElementNumbers(FacPerm(i)) 
             Vals(i) = Fac(i)

             IF(ImplicitSum < ImplicitLimit * FactorSum) THEN
               ImplicitSum = ImplicitSum + Fac(i)
               GebhartFactors % NumberOfImplicitFactors = i
             END IF
           END DO
           Vals(2:n) = Vals(2:n) * (FactorSum - Vals(1)) / (ConsideredSum - Vals(1))

           IF(ImplicitLimit < TINY(ImplicitLimit)) GebhartFactors % NumberOfImplicitFactors = 0       
           ImplicitEntries = ImplicitEntries + GebhartFactors % NumberOfImplicitFactors

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
           n = GebhartFactors % NumberOfFactors         
           Vals => GebhartFactors % Factors
           Cols => GebhartFactors % Elements

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
       END  DO

       MinSum = MINVAL(RowSums)
       MaxSum = MAXVAL(RowSums)

       WRITE(Message,'(A,T35,2ES15.4)') 'Minimum Gebhart factors sum',MINVAL(RowSums)
       CALL Info('RadiationFactors',Message,Level=5)
       WRITE(Message,'(A,T35,2ES15.4)') 'Maximum Gebhart factors sum',MAXVAL(RowSums)
       CALL Info('RadiationFactors',Message,Level=5)
       WRITE(Message,'(A,T35,ES15.4)') 'Maximum share of omitted factors',MaxOmittedFactor
       CALL Info('RadiationFactors',Message,Level=5)
       WRITE(Message,'(A,T35,ES15.4)') 'Gebhart factors filling (%)',(100.0 * MatrixEntries) / &
           (RadiationSurfaces**2)
       CALL Info('RadiationFactors',Message,Level=5)
       IF(ImplicitEntries > 0) THEN
         WRITE(Message,'(A,T35,ES15.4)') 'Implicit factors filling (%)',(100.0 * ImplicitEntries) / &
             (RadiationSurfaces**2)
         CALL Info('RadiationFactors',Message,Level=5)
       END IF

       IF(SaveFactors) THEN
         CALL SaveGebhartFactors()       
       END IF
       DEALLOCATE(RowSums)
     END SUBROUTINE CalculateGebhartFactors

     
     ! When Gebhart factors may have changed also modify the matrix topology so that
     ! when we assemble the matrices we are not hitting non-existing entries.
     !--------------------------------------------------------------------------------
     SUBROUTINE UpdateMatrixTopologyWithFactors()

       TYPE(Matrix_t), POINTER :: AMatrix
       LOGICAL :: OptimizeBW, GotIt
       INTEGER :: n,MatrixFormat
       INTEGER, POINTER :: NewPerm(:)
     
       CALL Info('RadiationFactors','Reorganizing the matrix structure',Level=5)

       MatrixFormat =  Tsolver % Matrix % FORMAT
       OptimizeBW = ListGetLogical(TSolver % Values,'Optimize Bandwidth',GotIt) 
       IF(.NOT. GotIt) OptimizeBW = .TRUE.

       CALL FreeMatrix( TSolver % Matrix)

       CALL Info('RadiationFactors','Creating new matrix topology')

       IF ( OptimizeBW ) THEN
         ALLOCATE( NewPerm( SIZE(Tsolver % Variable % Perm)), STAT=istat)
         IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 15.')

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
       ALLOCATE( AMatrix % RHS(n), STAT=istat)
       IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 16.')


       ! Transient case additional allocations:
       ! --------------------------------------
       IF ( ListGetString( CurrentModel % Simulation,'Simulation Type' ) == 'transient' ) THEN
         ALLOCATE( Amatrix % Force(n, TSolver % TimeOrder+1), STAT=istat )
         IF ( istat /= 0 ) CALL Fatal('RadiationFactors','Memory allocation error 17.')         
         Amatrix % Force = 0.0d0
       END IF

       TSolver % Matrix => Amatrix
       CALL ParallelInitMatrix( TSolver, AMatrix )

     END SUBROUTINE UpdateMatrixTopologyWithFactors


     SUBROUTINE CalculateRadiosity()
       IF (Newton) THEN
         ALLOCATE( RHS_d(RadiationSurfaces), SOL_d(RadiationSurfaces) )
         RHS_d = 0.0_dp
       END IF

       Temperature => Null()
       IF(ASSOCIATED(TSolver % Variable))  THEN
         TempPerm => TSolver % Variable % Perm
         Temperature => TSolver % Variable % Values
       END IF

       IF(.NOT.ASSOCIATED(Temperature)) &
         CALL Fatal('RadiationFactors', &
              "Radiosity solution can't be completed without the temperature field.")
       
       ALLOCATE(SurfaceTemperature(RadiationSurfaces))
       CALL TabulateSurfaceTemperatures()

       Sigma = ListGetConstReal( Model % Constants,&
         'Stefan Boltzmann',UnfoundFatal=.TRUE. )

       IF( Spectral ) THEN
         CALL SpectralRadiosity()
       ELSE
         CALL ConstantRadiosity()
       END IF

       ! This for debugging now!
       !PRINT *,'Range of SOL:',MINVAL(SOL),MAXVAL(SOL),SUM(SOL)
       !IF(Newton) THEN
       !  PRINT *,'Range of SOL_d:',MINVAL(SOL_d),MAXVAL(SOL_d), SUM(SOL_d)        
       !  DEALLOCATE(RHS_d, SOL_d)
       !END IF

       DEALLOCATE(SurfaceTemperature)
       IF(Newton) DEALLOCATE(RHS_d, SOL_d)
     END SUBROUTINE CalculateRadiosity
       

     ! Compute radiosity vector in the case that the emissivity is constant with temperature:
     !
     ! I.e. Solve (r*F-1)J = -e*Sigma*T^4 - r*R*P for J
     !  J: Radiosity
     !  e: Emissivity
     !  r: Reflectivity (=1-Emissivity)
     !  F: The "viewfactor" matrix
     !  Sigma: Stefan-Boltzmann constant
     !  T: Temperature
     !  P: "Radiator" power
     !  R: Elementwise visibility of the "radiators"
     !----------------------------------------------------------------------------------------
     SUBROUTINE ConstantRadiosity()
 
       LOGICAL :: RBC
       REAL(KIND=dp) :: r, e, a, Temp, Black
       REAL(KIND=dp), ALLOCATABLE :: RadiatorPowers(:)


       ! Assemble the equations, first coefficient matrix:
       ! -------------------------------------------------
       CALL RadiosityAssembly(RadiationSurfaces,G,Diag)

       ! ... and then the RHS:
       ! ---------------------

       DO i=1,RadiationSurfaces
         e = Emissivity(i)
         r = 1-e
         a = RelAreas(i) * (r/e)
         Temp = SurfaceTemperature(i)
         Black = Sigma*Temp**4
         RHS(i) = -a*e*Black
         IF(Newton) RHS_d(i) = -a*e*Black*4/Temp
       END DO

       ! Check for radiation sources:
       RBC = CheckForRadiators(RadiatorPowers)
       IF( RBC) THEN
         DO i=1,RadiationSurfaces
           Element => Mesh % Elements(ElementNumbers(i))
           IF ( ALLOCATED(Element % BoundaryInfo % Radiators)) THEN
             e = Emissivity(i)
             r = 1-e
             a = RelAreas(i) * (r/e)
             RHS(i) = RHS(i) - a*r* & 
                      SUM(Element % BoundaryInfo % Radiators*RadiatorPowers)
           END IF
         END DO
       END IF

       ! Solve for the radiosity and its derivative with respect to temperature
       !-----------------------------------------------------------------------
       CALL RadiationLinearSolver(RadiationSurfaces,G,SOL,RHS,Diag,Solver)
       IF( Newton ) THEN
         CALL RadiationLinearSolver(RadiationSurfaces,G,SOL_d,RHS_d, &
                      Diag, Solver, Scaling=.FALSE.)
       END IF

       ! Store the results for access by e.g. heat equation solvers:
       !------------------------------------------------------------
       CALL UpdateRadiosityFactors(SOL,SOL_d)
     END SUBROUTINE ConstantRadiosity
       
     
     ! Divide temperature into intervals.
     ! This is needed in order to compute problems where emissivity depends on temperature.
     !-------------------------------------------------------------------------------------
     SUBROUTINE SpectralRadiosity()
       REAL(KIND=dp) :: Tmin, Tmax, dT, Trad
       INTEGER :: k,kmin,kmax
       REAL(KIND=dp) :: q, qsum, totsum, a, r, e, s, Temp, Black
       REAL(KIND=dp), ALLOCATABLE :: tmpSOL(:), tmpSOL_d(:), EffAbs(:), EffTemp(:)

       LOGICAL :: RBC, ApproxNewton, AccurateNewton, UsedEdT
       INTEGER, ALLOCATABLE :: RadiatorSet(:)
       REAL(KIND=dp), ALLOCATABLE :: RadiatorPowers(:), RadiatorTemps(:)

       ApproxNewton = .FALSE.
       AccurateNewton = .FALSE.      
       IF( Newton ) THEN
         AccurateNewton = ListGetLogical( TSolver % Values,'Accurate Spectral Newton',Found ) 
         ApproxNewton = .NOT. AccurateNewton
       END IF
         
       IF(.NOT. ASSOCIATED(Temperature)) RETURN

       Tmin = MINVAL(SurfaceTemperature)
       Tmax = MAXVAL(SurfaceTemperature)
       
       WRITE(Message,'(A,ES12.3)') 'Minimum boundary temperature: ',Tmin
       CALL Info('SpectralRadiosity',Message,Level=10)
       WRITE(Message,'(A,ES12.3)') 'Maximum boundary temperature: ',Tmax
       CALL Info('SpectralRadiosity',Message,Level=10)
       
       IF( Tmin < 0.0_dp ) THEN
         CALL Fatal('SpectralRadiosity','Negative temperature not a good starting point!')
       END IF
       
       ! We have a fixed dT instead of having variable one related to Tmin and Tmax since
       ! adaptive intervals could generate funny attractors.
       dT = ListGetCReal( TSolver % Values,'Spectral dT',UnfoundFatal=.TRUE.) 
       
       kmin = FLOOR( Tmin / dT )
       kmax = CEILING( Tmax / dT )

       CALL Info('SpectralRadiosity','Going through discrete intervals: '&
           //TRIM(I2S(kmin))//'-'//TRIM(I2S(kmax)))

       SOL = 0.0_dp
       IF(Newton) SOL_d = 0.0_dp
       
       ALLOCATE( tmpSOL(RadiationSurfaces) )
       IF(Newton) ALLOCATE(tmpSOL_d(RadiationSurfaces))

       ALLOCATE(EffAbs(RadiationSurfaces),EffTemp(RadiationSurfaces))
       EffAbs = 0.0_dp
       EffTemp = 0.0_dp
       
       totsum = 0.0_dp

       DO k = kmin, kmax

         qsum = 0.0_dp         
         DO i=1,RadiationSurfaces           
           q = ( SurfaceTemperature(i) / dT - k )
           IF( ABS(q) < 1 ) THEN
             q = 1-ABS(q)
             qsum = qsum + q
           END IF
         END DO

         ! There is nothing to compute here
         ! So no need to resolve equations for this interval.        
         IF(qsum < 1.0d-6 ) THEN
           CALL Info('SpectralRadiosity','Skipping interval '//TRIM(I2S(k)),Level=12)
           CYCLE
         END IF

         WRITE(Message,'(A,G12.5)') 'Spectral radiosity sources '//TRIM(I2S(k))//': ',qsum
         CALL Info('SpectralRadiosity',Message,Level=10) 
                    
         ! Initialize matrix equation
         Diag = 0.0_dp
         RHS  = 0.0_dp         
         IF(Newton) RHS_d = 0.0_dp

         G % Values = 0.0_dp
         
         ! This is the temperature under study for which we will get the emissivities for. 
         Trad = k*dT         
         CALL TabulateSpectralEmissivity(Emissivity,Trad)
         CALL RadiosityAssembly(RadiationSurfaces,G,Diag)
         DO i=1,RadiationSurfaces
           e = Emissivity(i)
           r = 1-e
           a = RelAreas(i) * (r/e)
           Temp = SurfaceTemperature(i)
           Black = Sigma*Temp**4
           ! The portion of the emissivity to consider for this radiating element
           q = ( SurfaceTemperature(i) / dT - k )
           IF( ABS(q) < 1 ) THEN
             ! As a weight we use linear interpolation.
             ! Perfect hit get weight 1 that goes to zero when hitting next temperature interval. 
             q = 1-ABS(q)
             RHS(i) = -q*a*e*Black
             IF (AccurateNewton) RHS_d(i) = 4*RHS(i)/Temp

             S = -q*e**2/r*Black
             SOL(i) = SOL(i) + S
             IF(Newton) SOL_d(i) = SOL_d(i) + 4*S/Temp
           END IF
         END DO

         ! This is a checksum since integration over all temperature intervals should go through all the
         ! participating surface elements. 
         totsum = totsum + qsum
         !PRINT *,'Trad:',k,Trad,qsum,totsum
         CALL RadiationLinearSolver(RadiationSurfaces,G,tmpSOL,RHS,Diag,Solver)

         ! Newton linearization including only "self"
         IF( ApproxNewton ) THEN
           tmpSOL_d = (4.0_dp/Trad) * tmpSOL
         ELSE IF( AccurateNewton ) THEN
           CALL RadiationLinearSolver(RadiationSurfaces,G,tmpSOL_d, &
                        RHS_d,Diag,Solver,Scaling=.FALSE.)
         END IF
         
         ! Cumulative radiosity         
         SOL = SOL + tmpSOL
         IF( Newton ) SOL_d = SOL_d + tmpSOL_d

         EffTemp = EffTemp + Trad * tmpSOL
         EffAbs = EffAbs + Emissivity * tmpSOL 
       END DO

       ! This should be exactly one!
       WRITE(Message,'(A,G12.5)') 'Checksum for radiosity sources: ',totsum / RadiationSurfaces 
       CALL Info('SpectralRadiosity',Message,Level=5) 

       ! Check for radiation sources:
       RBC = CheckForRadiators(RadiatorPowers,RadiatorTemps)
       IF(RBC) THEN
         Tmin = MINVAL(RadiatorTemps)
         Tmax = MAXVAL(RadiatorTemps)

         IF(ABS(Tmin-Tmax) < 1.0e-6 ) THEN           
           WRITE(Message,'(A,ES12.3)') 'Only radiator temperature: ',Tmin
           CALL Info('SpectralRadiosity',Message,Level=10)
         ELSE           
           WRITE(Message,'(A,ES12.3)') 'Minimum radiator temperature: ',Tmin
           CALL Info('SpectralRadiosity',Message,Level=10)
           WRITE(Message,'(A,ES12.3)') 'Maximum radiator temperature: ',Tmax
           CALL Info('SpectralRadiosity',Message,Level=10)
         END IF

         ALLOCATE( RadiatorSet(SIZE(RadiatorTemps)) )
         RadiatorSet = 0
         k = 0
         DO i=1,SIZE(RadiatorTemps)
           IF(RadiatorSet(i) > 0) CYCLE
           k=k+1
           RadiatorSet(i) = k
           DO j=i+1,SIZE(RadiatorTemps)
             IF(ABS(RadiatorTemps(i)-RadiatorTemps(j)) < 1.0e-6) RadiatorSet(j) = RadiatorSet(i)
           END DO
         END DO
         kmax = k
         CALL Info('SpectralRadiosity','Going through radiators in '//TRIM(I2S(kmax))//' sets')
                           
         DO k = 1, kmax
           DO j=1,SIZE(RadiatorSet)
             IF(RadiatorSet(j) == k) Trad = RadiatorTemps(j)
           END DO
           
           WRITE(Message,'(A,G12.5)') 'Spectral radiosity radiators '//TRIM(I2S(k))//' at: ',Trad
           CALL Info('SpectralRadiosity',Message,Level=10) 
           
           ! Initialize matrix equation
           Diag = 0.0_dp
           RHS = 0.0_dp         
           G % Values = 0.0_dp
           
           CALL TabulateSpectralEmissivity(Emissivity,Trad)
           CALL RadiosityAssembly(RadiationSurfaces,G,Diag)
           DO i=1,RadiationSurfaces
             Element => Mesh % Elements(ElementNumbers(i))
             IF(ALLOCATED(Element % BoundaryInfo % Radiators)) THEN
               e = Emissivity(i)
               r = 1-e
               a = RelAreas(i) * (r/e)
               DO j=1,SIZE(RadiatorSet)
                 IF(RadiatorSet(j) == k) THEN
                   RHS(i) = RHS(i) - a * r * &
                       Element % BoundaryInfo % Radiators(j) * RadiatorPowers(j)
                 END IF
               END DO
             END IF
           END DO

           CALL RadiationLinearSolver(RadiationSurfaces,G,tmpSOL,RHS,Diag,Solver)          
           
           ! Cumulative radiosity
           SOL = SOL + tmpSOL
           EffTemp = EffTemp + Trad * tmpSOL
           EffAbs = EffAbs + Emissivity * tmpSOL 
         END DO
       END IF
       
       ! Normalize with weight i.e. incoming heat flux
       EffAbs = EffAbs / SOL
       EffTemp = EffTemp / SOL

       ! Store the results for access by e.g. heat equation solvers:
       !------------------------------------------------------------
       CALL UpdateRadiosityFactors(SOL,SOL_d,EffAbs,EffTemp)
     END SUBROUTINE SpectralRadiosity
     

     FUNCTION CheckForRadiators(RadiatorPowers,RadiatorTemps) RESULT(RBC)
       LOGICAL :: RBC
       REAL(KIND=dp), ALLOCATABLE :: RadiatorPowers(:)
       REAL(KIND=dp), ALLOCATABLE, OPTIONAL :: RadiatorTemps(:)

       TYPE(ValueList_t), POINTER :: RadList
       INTEGER :: n
       LOGICAL :: Found
       REAL(KIND=dp), POINTER :: RadiatorCoords(:,:), rWrk(:,:)

       ! If radiator is in body force section then use it:
       ! This will make it easier to make GUIs etc.
       IF( .NOT. ListCheckPresentAnyBodyForce( Model,'Radiator Coordinates',RadList ) ) &
             RadList => TSolver % Values
       CALL GetConstRealArray( RadList, RadiatorCoords, 'Radiator Coordinates',RBC)

       IF(RBC) THEN
         n = SIZE(RadiatorCoords,1)
         ALLOCATE( RadiatorPowers(n))
         CALL GetConstRealArray( RadList, rWrk, 'Radiator Power', Found ) 
         IF( Found ) THEN
           IF(SIZE(rWrk,1)==1) THEN
             RadiatorPowers(1:n) = rWrk(1,1)
           ELSE IF(SIZE(rWrk,1)==n) THEN
             RadiatorPowers(1:n) = rWrk(1:n,1)
           ELSE
             CALL Fatal('ConstantRadiosity','Mismatch between size of "Radiator Coordinates" and "Radiator Power"')
           END IF
         ELSE
           DO t=1,n
             RadiatorPowers(t) = ListGetCReal(RadList, 'Radiator Power '//TRIM(I2S(t)),UnfoundFatal=.TRUE.)
           END DO
         END IF

         IF(PRESENT(RadiatorTemps)) THEN
           ALLOCATE( RadiatorTemps(n))
           CALL GetConstRealArray( RadList, rWrk, 'Radiator Temperature', Found )
           IF( Found ) THEN
             IF(SIZE(rWrk,1)==1) THEN
               RadiatorTemps(1:n) = rWrk(1,1)
             ELSE IF(SIZE(rWrk,1)==n) THEN
                 RadiatorTemps(1:n) = rWrk(1:n,1)
             ELSE
               CALL Fatal('SpectralRadiosity','Mismatch between size of "Radiator Coordinates" and "Radiator Power"')
             END IF
           ELSE
             DO t=1,n
               RadiatorTemps(t) = ListGetCReal(RadList, 'Radiator Temperature '//TRIM(I2S(t)),UnfoundFatal=.TRUE.)
             END DO
           END IF
         END IF
       END IF
     END FUNCTION CheckForRadiators


     SUBROUTINE RadiosityAssembly(n,G,Diag)
       TYPE(Matrix_t) :: G
       INTEGER :: n
       REAL(KIND=dp) :: Diag(:)

       REAL(KIND=dp), POINTER :: Vals(:) 
       INTEGER, POINTER :: Cols(:) 
       REAL(KIND=dp) :: s,r,e,rj,ej,a
       INTEGER :: i, j, nf, colj, previ

       DO i=1,n
         nf = ViewFactors(i) % NumberOfFactors
         Vals => ViewFactors(i) % Factors
         Cols => ViewFactors(i) % Elements
         e = Emissivity(i)
         r = 1-e
         a = RelAreas(i) * (r/e)**2
         previ = G % Rows(i)-1
         DO j=1,nf
           ej = Emissivity(Cols(j))
           rj = 1-ej
           s = r*Vals(j) * (r/e*rj/ej)
           CALL CRS_AddToMatrixElement(G,i,Cols(j),s,previ)
         END DO
         CALL CRS_AddToMatrixElement(G,i,i,-a)
       END DO
       Diag = G % Values(G % Diag)
     END SUBROUTINE RadiosityAssembly


     
     ! Scale linear system & solve
     !-----------------------------
     SUBROUTINE RadiationLinearSolver(n, A, x, b, Diag,  Solver, Scaling)

       INTEGER :: n
       REAL(KIND=dp) :: x(:), b(:), Diag(:)
       TYPE(Matrix_t), POINTER :: A
       LOGICAL, OPTIONAL :: Scaling
       TYPE(Solver_t), POINTER :: Solver

       LOGICAL :: Scal,Found
       REAL(KIND=dp) :: bscal

       scal = .TRUE.
       IF(PRESENT(Scaling)) scal = Scaling

       ! Scale matrix to unit diagonals (if not done already)
       IF(scal) THEN
         Diag = SQRT(1._dp/ABS(Diag))
         DO i=1,n
           DO j=A % Rows(i),A % Rows(i+1)-1
             A % Values(j) = A % Values(j)*Diag(i)*Diag(A % Cols(j))
           END DO
         END DO
       END IF
       b = b * Diag

       ! Scale rhs to one!
       bscal = SQRT(SUM(b**2))
       b = b / bscal

       x = 0.0_dp
       IF(IterSolveFactors) THEN
         Solver % Matrix => A
         CALL IterSolver( A, x, b, Solver )
       ELSE           
         CALL DirectSolver( A, x, b, Solver )
       END IF
       x = x * bscal * Diag
     END SUBROUTINE RadiationLinearSolver

#if TESTCG
     ! Tailored local CG algo for speed testing (somewhat faster than any of the 
     ! library routines but not so much...)
     !-------------------------------------------------------------------------
     SUBROUTINE  RadiationCG( n, A, x, b, eps)
       REAL(KIND=dp) :: x(:),b(:), eps
       INTEGER :: n
       TYPE(Matrix_t), POINTER :: A

       REAL(KIND=dp):: alpha, beta, rho, oldrho
       REAL(KIND=dp) :: r(n), p(n), q(n), z(n)
       INTEGER :: iter, i
       REAL(KIND=dp) :: residual, eps2

       eps2 = eps*eps

       CALL CRS_MatrixVectorMultiply(A,x,r) 
       r = b - r
       residual = SUM(r*r)
       IF(residual<eps2) RETURN

       DO iter=1,100
         rho = SUM(r*r)
         IF(rho==0.0_dp) STOP 'CG, rho=0'
  
         IF ( iter==1 ) THEN
           p = r
         ELSE
           beta = rho / oldrho
           p = r + beta * p
         END IF

         CALL CRS_MatrixVectorMultiply(A,p,q) 
         alpha = rho/SUM(p*q)

         x = x + alpha * p
         r = r - alpha * q
         residual = SUM(r*r)
         IF ( residual < eps2) EXIT

         oldrho = rho
       END DO

       CALL CRS_MatrixVectorMultiply(A,x,r) 
       r = b - r
       residual = SQRT(SUM(r*r))
       WRITE (*, '(I8, E11.4)') iter, residual
     END SUBROUTINE RadiationCG
#endif


     SUBROUTINE UpdateRadiosityFactors(SOL,SOL_d,EffAbs,EffTemp)
       REAL(KIND=dp) :: SOL(:), SOL_d(:)
       REAL(KIND=dp), OPTIONAL :: EffAbs(:), EffTemp(:)

       TYPE(Element_t), POINTER :: Element
       INTEGER :: n,i
       TYPE(Factors_t), POINTER :: RadiationFactors

       n=2
       IF(PRESENT(EffAbs))  n=n+1
       IF(PRESENT(EffTemp)) n=n+1

       DO i=1,RadiationSurfaces
         Element => Mesh % Elements(ElementNumbers(i))

         RadiosityFactors => Element % BoundaryInfo % RadiationFactors       
         IF ( .NOT. ASSOCIATED( RadiosityFactors ) ) THEN
           ALLOCATE(RadiosityFactors)
           Element % BoundaryInfo % RadiationFactors => RadiosityFactors
           Element % BoundaryInfo % RadiationFactors % Elements => Null()
         END IF

         IF (.NOT.ASSOCIATED(RadiosityFactors % Elements)) THEN
           ALLOCATE( RadiosityFactors % Elements(1) )
           ALLOCATE( RadiosityFactors % Factors(4) )
           RadiosityFactors % Factors = 0.0_dp
           RadiosityFactors % NumberOfFactors = 1
           RadiosityFactors % Elements(1) = ElementNumbers(i)
         END IF

         RadiosityFactors % Factors(1) = SOL(i)
         IF(Newton) RadiosityFactors % Factors(2) = SOL_d(i)
         IF(PRESENT(EffAbs))  RadiosityFactors % Factors(3) = EffAbs(i)
         IF(PRESENT(EffTemp)) RadiosityFactors % Factors(4) = EffTemp(i)
       END DO
     END SUBROUTINE UpdateRadiosityFactors



     ! Save factors is mainly for debugging purposes
     !-------------------------------------------------------------------
     SUBROUTINE SaveGebhartFactors()

       CHARACTER(:), ALLOCATABLE :: OutputName, GebhartFactorsFile
       INTEGER ::  i,t,n
       INTEGER, POINTER :: Cols(:)
       REAL(KIND=dp), POINTER :: Vals(:)

       GebhartFactorsFile = GetString(Model % Simulation, 'Gebhart Factors',GotIt )
       IF (.NOT.GotIt) &
         GebhartFactorsFile = GetString(Model % Simulation, 'Gebhardt Factors',GotIt )

       IF ( .NOT.GotIt ) THEN
         GebhartFactorsFile = 'GebhartFactors.dat'
       END IF

       IF ( LEN_TRIM(Model % Mesh % Name) > 0 ) THEN
         OutputName = TRIM(OutputPath) // '/' // TRIM(Model % Mesh % Name) // &
             '/' // TRIM(GebhartFactorsFile)
       ELSE
         OutputName = TRIM(GebhartFactorsFile) 
       END IF

       IF(RadiationBody > 1) OutputName = OutputName//TRIM(I2S(RadiationBody))

       OPEN( VFUnit,File=TRIM(OutputName) )

       WRITE (Message,'(A,A)') 'Writing Gephardt Factors to file: ',TRIM(OutputName)
       CALL Info('RadiationFactors',Message,Level=5)

       WRITE( VFUnit,* ) RadiationSurfaces

       DO t=1,RadiationSurfaces
         WRITE(VFUnit,*) t,ElementNumbers(t)
       END DO

       DO t=1,RadiationSurfaces
         Element => Mesh % Elements(ElementNumbers(t))
         GebhartFactors => Element % BoundaryInfo % RadiationFactors

         n = GebhartFactors % NumberOfFactors 
         Vals => GebhartFactors % Factors
         Cols => GebhartFactors % Elements

         WRITE( VFUnit,* ) n
         DO i=1,n
           WRITE(VFUnit,*) t,InvElementNumbers(Cols(i)-j0),Vals(i)
         END DO
       END DO

       CLOSE(VFUnit)
     END SUBROUTINE SaveGebhartFactors


     SUBROUTINE InitFactorSolver(TSolver, Solver)
       TYPE(Solver_t) :: TSolver
       TYPE(Solver_t), POINTER :: Solver

       IF(.NOT. ASSOCIATED(Solver)) ALLOCATE(Solver)
       IF(.NOT. ASSOCIATED(Solver % Values)) Solver % Values => ListAllocate()

       CALL ListCopyPrefixedKeywords( TSolver % Values, Solver % Values, 'radiation:' )

       CALL ListAddNewString(Solver % Values,'Linear System Iterative Method', 'CGS')
       CALL ListAddNewString(Solver % Values,'Linear System Direct Method','Umfpack')
       CALL ListAddNewInteger(Solver % Values,'Linear System Max Iterations',500)
       CALL ListAddNewInteger(Solver % Values,'Linear System Residual Output',10 )
       CALL ListAddNewString(Solver % Values,'Linear System Preconditioning','None' )
       CALL ListAddNewConstReal(Solver % Values,'Linear System Convergence Tolerance',1.0d-9)
     END SUBROUTINE InitFactorSolver

   END SUBROUTINE RadiationFactors
