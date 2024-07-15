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
! * ELMER/FEM Viewfactor computation
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Jun 1997
! *
! *****************************************************************************/

!> \ingroup Programs
!> \{

!> \defgroup ViewFactors Program RadiatorFactors
!> \{

!------------------------------------------------------------------------------
!> A separate program that computes the radiatiave source coefficients to an
!> external file. 
!>
!> This file is later used within the ElmerSolver. If the radiative files
!> do not exist, a system call for this program is performed. 
!------------------------------------------------------------------------------
      
   PROGRAM RadiatorFactors
   
     USE DefUtils

     IMPLICIT NONE

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Model_t), POINTER :: Model
     TYPE(Mesh_t), POINTER  :: Mesh
     TYPE(Solver_t), POINTER  :: Solver

     INTEGER :: i,j,k,l,t,k1,k2,n,iter,Ndeg,Time,NSDOFs,MatId,istat

     REAL(KIND=dp) :: SimulationTime,dt,s,a1,a2,FMin,FMax

     INTEGER, ALLOCATABLE ::  Surfaces(:), TYPE(:)
     REAL(KIND=dp), ALLOCATABLE :: Coords(:),Normals(:),Factors(:)

     TYPE(Element_t),POINTER :: Element, Parent

     INTEGER :: BandSize,SubbandSize,RadiationSurfaces,Row,Col
     INTEGER, DIMENSION(:), POINTER :: Perm

     REAL(KIND=dp) :: Norm,PrevNorm,MinFactor, Normal_in, Plane
 
     TYPE(Nodes_t) :: ElementNodes
     TYPE(ValueList_t), POINTER :: BC, Material, Params, RadList

     INTEGER :: LeftNode,RightNode,LeftBody,RightBody,RadBody
     REAL(KIND=dp) :: NX,NY,NZ,NRM(3),DensL,DensR

     INTEGER :: divide, nprob, NofRadiators
     REAL(KIND=dp) :: AreaEPS, RayEPS, FactEPS
     REAL(KIND=dp) :: at0, rt0
     CHARACTER(*), PARAMETER :: Caller = 'Radiators'

     INTERFACE
        ! void radiatiors3d
        ! ( int *EL_N,  int *EL_Topo, int *EL_Type, double *EL_Coord, double *EL_Normals,
        ! int *RT_N0, int *RT_Topo0, int *RT_Type, double *RT_Coord, double *RT_Normals,
        ! double *Factors, double *Feps, double *Aeps, double *Reps, int *Nr, 
        ! int *NInteg,int *NInteg3, int  *Combine )
        SUBROUTINE radiatorfactors3d(EL_N,  EL_Topo, EL_Type, EL_Coord, EL_Normals, &
                                 RT_N0, RT_Topo0, RT_Type, RT_Coord, RT_Normals, &
                                 NofRadiators, RadiatorCoords, &
                                 Factors, Feps, Aeps, Reps, Nr, NInteg, NInteg3, Combine) BIND(C)

            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT) :: EL_N, EL_Topo(*), EL_Type(*)
            REAL(KIND=C_DOUBLE) :: EL_Coord(*), EL_Normals(*)
            INTEGER(C_INT) :: RT_N0, RT_Topo0(*), RT_Type(*), NofRadiators
            REAL(KIND=C_DOUBLE) :: RT_Coord(*), RT_Normals(*), Factors(*), RadiatorCoords(*)
            REAL(KIND=C_DOUBLE) :: Feps, Aeps, Reps
            INTEGER(C_INT) :: Nr, NInteg, NInteg3, Combine
        END SUBROUTINE radiatorfactors3d
        
      END INTERFACE

     INTEGER, POINTER :: Timesteps(:)
     INTEGER :: TimeIntervals,interval,timestep,combineInt
     
     LOGICAL :: CylindricSymmetry,GotIt, Found, Radiation, LeftEmis, RightEmis

     CHARACTER(LEN=MAX_NAME_LEN) :: eq, &
           ViewFactorsFile,OutputName,ModelName,LMessage, TempString

     TYPE(Element_t), POINTER :: RadElements(:)
     INTEGER :: RadiationBody, MaxRadiationBody, Nrays
     LOGICAL :: RadiationOpen, Combine, RadiationFlag
     INTEGER, PARAMETER :: VFUnit = 10
     INTEGER :: iostat, NoArgs
     REAL(KIND=dp), POINTER :: Radiators(:,:)
     

     CALL Info( Caller, ' ', Level=3 )
     CALL Info( Caller, '==============================================', Level=3 )
     CALL Info( Caller, ' E L M E R  R A D I A T O R S,  W E L C O M E',  Level=3  )
     CALL Info( Caller, '==============================================', Level=3 )

!------------------------------------------------------------------------------
!    Read element definition file, and initialize element types
!------------------------------------------------------------------------------
     CALL InitializeElementDescriptions
!------------------------------------------------------------------------------
!    Read Model from Elmer Data Base
!------------------------------------------------------------------------------
     CALL Info( Caller, ' ', Level=3 )
     CALL Info( Caller, ' ', Level=3 )
     CALL Info( Caller, 'Reading Model... ', Level=3 )

!------------------------------------------------------------------------------
     NoArgs = COMMAND_ARGUMENT_COUNT()
     IF ( NoArgs > 0 ) THEN
       CALL GET_COMMAND_ARGUMENT(1, ModelName)
     ELSE
       OPEN( 1,file='ELMERSOLVER_STARTINFO', STATUS='OLD', IOSTAT=iostat )
       IF( iostat /= 0 ) THEN
         CALL Fatal( Caller, 'Unable to find ELMERSOLVER_STARTINFO, cannot execute.' )
       END IF
       READ(1,'(a)') ModelName
       CLOSE(1)
     END IF

     Model => LoadModel( ModelName,.FALSE.,1,0 )

     CurrentModel => Model
          
     NULLIFY( Mesh )
     DO i=1,Model % NumberOfSolvers
       Solver => Model % Solvers(i)
       Radiation = ListGetLogical( Solver % Values, 'Radiation Solver', Found )
       IF ( Radiation ) THEN
         Mesh => Solver % Mesh
         Model % Solver => Solver
         EXIT
       ENDIF
     END DO
     
     IF ( .NOT. ASSOCIATED(Mesh) ) THEN
       CALL Fatal(Caller,'No heat equation definition. Cannot compute factors.')
     END IF

     Params => GetSolverParams()

     IF( .NOT. ListCheckPresentAnyBodyForce( Model,'Radiator Coordinates',RadList ) ) &
         RadList => Params
     
     CALL GetConstRealArray( RadList, Radiators, 'Radiator Coordinates', Found )
     IF(.NOT. Found ) CALL Fatal( 'RadiatorFactors', 'No radiators present, quitting' )

     NofRadiators = SIZE(Radiators,1)
     CALL Info( 'RadiatorFactors', 'Computing flux coefficients for ' // &
             I2S(NofRadiators) // ' radiative sources', LEVEL=5 )

     
#define SYMMETRY_NOW .TRUE.
     IF(SYMMETRY_NOW) THEN
       DO i=1,6
         SELECT CASE(i)
         CASE(1)
           Plane = GetCReal( Params, 'Viewfactor Symmetry x min', Found );
           IF(.NOT. Found ) Found = ListGetLogical( Params, 'Viewfactor Symmetry x', GotIt );
         CASE(2)
           Plane = GetCReal( Params, 'Viewfactor Symmetry x max', Found );
         CASE(3)
           Plane = GetCReal( Params, 'Viewfactor Symmetry y min', Found );
           IF(.NOT. Found ) Found = ListGetLogical( Params, 'Viewfactor Symmetry y', GotIt );
         CASE(4)
           Plane = GetCReal( Params, 'Viewfactor Symmetry y max', Found );
         CASE(5)
           Plane = GetCReal( Params, 'Viewfactor Symmetry z min', Found );
           IF(.NOT. Found ) Found = ListGetLogical( Params, 'Viewfactor Symmetry z', GotIt );
         CASE(6)
           Plane = GetCReal( Params, 'Viewfactor Symmetry z max', Found );
         END SELECT

         IF(.NOT. Found ) CYCLE

         CALL Info(Caller,'Duplicating mesh in coordinate direction: '//I2S((i+1)/2))

         CALL MirrorMesh(Mesh, i, Plane)
       END DO
     END  IF

     CALL SetCurrentMesh( Model,Mesh )

     CALL Info( Caller,'Number of nodes in mesh: '//I2S(Mesh % NumberOfNodes),Level=7)

!------------------------------------------------------------------------------
!    Figure out requested coordinate system
!------------------------------------------------------------------------------
     eq = GetString( GetSimulation(), 'Coordinate System' )
     SELECT CASE(eq)
     CASE('cartesian','cartesian 2d','cartesian 3d')
       Coordinates = Cartesian

     CASE DEFAULT
       CALL Error( Caller, &
         'Unknown Global Coordinate System for Viewfactor computation ')
       CALL Error( Caller, TRIM(eq) )
       CALL Fatal( Caller, &
         'Only Cartesian 2D/3D coordinates allowed. Aborting' )
     END SELECT

     CylindricSymmetry = (Coordinates == AxisSymmetric) .OR. (Coordinates==CylindricSymmetric)
!------------------------------------------------------------------------------

     ALLOCATE( ElementNodes % x(Model % MaxElementNodes), &
         ElementNodes % y(Model % MaxElementNodes), &
         ElementNodes % z(Model % MaxElementNodes),STAT=istat )
     
     IF ( CylindricSymmetry ) THEN
       ALLOCATE( Coords(2 * Mesh % NumberOfNodes), STAT=istat )
       DO i=1,Mesh % NumberOfNodes
         Coords(2*(i-1)+1) = Mesh % Nodes % x(i)
         Coords(2*(i-1)+2) = Mesh % Nodes % y(i)
       END DO
     ELSE
       ALLOCATE( Coords(3 * Mesh % NumberOfNodes), STAT=istat )
       DO i=1,Mesh % NumberOfNodes
         Coords(3*(i-1)+1) = Mesh % Nodes % x(i)
         Coords(3*(i-1)+2) = Mesh % Nodes % y(i)
         Coords(3*(i-1)+3) = Mesh % Nodes % z(i)
       END DO
     END IF

     IF ( istat /= 0 ) THEN
       CALL Fatal('Radiatorfactors', 'Memory allocation error. Aborting' )
     END IF

     ! The routine originally plays with the element list and therefore
     ! when several radiation boundaries are needed both the original and
     ! the new elementlist needs to be in the memory. Thus the hassle.

     MinFactor = ListGetConstReal(Params,'Minimum Radiator Factor',GotIt)
     IF(.NOT. GotIt) MinFactor = 1.0d-20

     CALL AllocateVector( RadElements, Mesh % NumberOfBoundaryElements, Caller )

     IF( Mesh % NumberOfBoundaryElements == 0) THEN
       CALL Warn(Caller,'There are no boundary elements at all!')
       STOP
     END IF

! Check the maximum radiation body
     MaxRadiationBody = 0

#if 0
     DO t= 1, Mesh % NumberOfBoundaryElements

       Element => GetBoundaryElement(t)
       IF ( GetElementFamily() == 1 ) CYCLE
       BC => GetBC()
       IF ( .NOT. ASSOCIATED( BC ) ) CYCLE
 
       RadiationFlag = GetLogical( BC, 'Radiator BC',GotIt )

       IF ( RadiationFlag ) THEN
         i = MAX(1, GetInteger( BC, 'Radiation Boundary', GotIt ) )
         MaxRadiationBody = MAX(i, MaxRadiationBody)
       END IF
     END DO
#else
     MaxRadiationBody = 1
#endif

     DO RadiationBody = 1, MaxRadiationBody
       WRITE( LMessage,'(A,I2)') 'Computing view factors for radiation body',RadiationBody
       CALL Info(Caller,LMessage,Level=3)
    
!------------------------------------------------------------------------------
!    Here we start...
!------------------------------------------------------------------------------
       RadiationSurfaces = 0
       RadiationOpen = .FALSE.

!------------------------------------------------------------------------------
!    loop to get the surfaces participating in radiation, discard the rest
!    of the elements...
!------------------------------------------------------------------------------

       DO t=1,Mesh % NumberOfBoundaryElements
         Element => GetBoundaryElement(t)
         IF ( GetElementFamily() == 1 ) CYCLE

         BC => GetBC()
         IF ( .NOT. ASSOCIATED( BC ) ) CYCLE
         
         RadiationFlag = GetLogical( BC, 'Radiator BC', GotIt )
         IF ( RadiationFlag ) THEN
#if 0
           i = MAX(1, GetInteger( BC, 'Radiation Boundary', GotIt ))
           IF(i == RadiationBody) THEN
             RadiationOpen = RadiationOpen .OR. GetLogical( BC, 'Radiation Boundary Open', GotIt )
             RadiationSurfaces = RadiationSurfaces + 1
             j = t + Mesh % NumberOFBulkElements
             RadElements(RadiationSurfaces) = Mesh % Elements(j)
           END IF
#else
           RadiationOpen = RadiationOpen .OR. GetLogical(BC,'Radiation Boundary Open', GotIt)
           RadiationSurfaces = RadiationSurfaces + 1
           j = t + Mesh % NumberOFBulkElements
           RadElements(RadiationSurfaces) = Mesh % Elements(j)
#endif
         END IF
       END DO
       
       
       N = RadiationSurfaces
       
       IF ( N == 0 ) THEN
         CALL Warn( 'RadiatorFactors', 'No surfaces participating in radiation?' )
         IF(RadiationBody < MaxRadiationBody) THEN
           CYCLE
         ELSE
           CALL Warn( 'RadiatorFctors', 'Stopping cause nothing to be done...' )
           STOP
         END IF
       END IF
       
       CALL Info(Caller,'Number of surfaces participating in radiation: '//I2S(N))
       
       IF ( CylindricSymmetry ) THEN
         ALLOCATE( Surfaces(2*N), Factors(N*N), STAT=istat )
       ELSE
         ALLOCATE( Normals(3*n), Factors(NofRadiators*n),Surfaces(4*n), TYPE(n), STAT=istat )
       END IF
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'RadiatonFactors', 'Memory allocation error. Aborting' )
       END IF
       
       DO t=1,N

         Element => RadElements(t)
         Model % CurrentElement => Element
         k = GetElementNOFNodes()
         CALL GetElementNodes(ElementNodes)
         
         IF ( GetElementFamily() == 3 ) THEN
           nrm = NormalVector( Element,ElementNodes,1/3._dp, 1/3._dp )
         ELSE
           nrm = NormalVector( Element, ElementNodes, 0.0_dp, 0.0_dp )
         END IF
         
         LeftBody = 0
         LeftNode = -1
	 LeftEmis = .FALSE.
         Parent => Element % BoundaryInfo % Left
         IF ( ASSOCIATED(Parent) ) THEN
           LeftBody  = Parent % BodyId
           DO i=1,Parent % TYPE % NumberOfNodes
             gotIt =.TRUE.
             DO j=1,Element % TYPE % NumberOfNodes
               IF ( Element % NodeIndexes(j)==Parent % NodeIndexes(i)) THEN
                 gotIt=.FALSE.
                 EXIT
               END IF
             END DO
             IF (gotIt) THEN
               LeftNode = Parent % NodeIndexes(i)
               EXIT
             END IF
           END DO
	   IF( LeftBody > 0 ) THEN 
             MatId = GetInteger( Model % Bodies(LeftBody) % Values,'Material', GotIt)
	     IF( MatId == 0 ) THEN
               CALL Warn(Caller,'Invalid material index in body, perhaps none')
             END IF 
           ELSE
             CALL Warn(Caller,'LeftBody not associated')
           END IF
           LeftEmis = ListCheckPresent(Model % Materials(MatId) % Values,'Emissivity') 
         END IF
         
         RightBody = 0
         RightNode=-1
         RightEmis = .FALSE.
         Parent => Element % BoundaryInfo % Right
         IF ( ASSOCIATED(Parent) ) THEN
           RightBody = Parent % BodyId
           DO i=1,Parent % TYPE % NumberOfNodes
             gotIt =.TRUE.
             DO j=1,Element % TYPE % NumberOfNodes
               IF (Element % NodeIndexes(j)==Parent % NodeIndexes(i)) THEN
                 gotIt=.FALSE.
                 EXIT
               END IF
             END DO
             IF (gotIt) THEN
               RightNode = Parent % NodeIndexes(i)
               EXIT
             END IF
           END DO
	   
           MatId = GetInteger( Model % Bodies(RightBody) % Values,'Material', GotIt)
           RightEmis=ListCheckPresent(Model % Materials(MatId) % Values,'Emissivity') 
         END IF


         BC => GetBC()
         RadBody = GetInteger( BC, 'Radiation Target Body',GotIt )
         IF ( .NOT. GotIt ) RadBody = GetInteger( BC, 'Normal Target Body',GotIt )

         IF ( .NOT. Gotit ) THEN
	   ! If emissivity given in either side then make the radiation target body to be the other.
           !----------------------------------------------------------------------------------------
           IF( RightEmis .AND. LeftEmis ) THEN
             CALL Warn(Caller,'Emissivity defined on both sides!')
             WRITE(Message,'(A,I3,I3,A,I3,A,I5)') 'Bodies:',RightBody,LeftBody,' BC:', &
                            GetBCId( Element ),' Ind:',t
             CALL Info(Caller,Message)
             IF( ASSOCIATED(Element % BoundaryInfo % Left, Element % BoundaryInfo % Right)) THEN
              CALL Warn(Caller,'Parents of the boundary element are the same')
              RadBody = LeftBody
             END IF
           ELSE IF( RightEmis ) THEN
             RadBody = LeftBody
           ELSE IF( LeftEmis ) THEN
             RadBody = RightBody
           END IF
         END IF

         IF ( RadBody < 0 ) RadBody = 0
         
         IF ( RadBody>0 .AND. (RadBody /= RightBody .AND. RadBody /= LeftBody) ) THEN
           CALL Error( Caller, 'Inconsistent direction information (Radiation Target Body)' )
           WRITE( LMessage, * ) 'Radiation Target: ', RadBody, ' Left, Right: ', LeftBody, RightBody
           CALL Fatal( Caller, LMessage )
         END IF
         
         IF ( LeftNode <= 0 .OR. (RadBody>0 .AND. RadBody==RightBody) ) &
           LeftNode = RightNode
         Normal_in = -1.0
         IF ( RadBody <= 0 ) Normal_in = 1.0
         
         nx = SUM(ElementNodes % x)/k - Mesh % Nodes % x(LeftNode)
         ny = SUM(ElementNodes % y)/k - Mesh % Nodes % y(LeftNode)
         nz = SUM(ElementNodes % z)/k - Mesh % Nodes % z(LeftNode)
         
         IF ( CylindricSymmetry ) THEN
           IF ( Normal_in * (Nrm(1)*nx + Nrm(2)*Ny + Nrm(3)*nz) > 0 ) THEN
             Surfaces(2*(t-1)+1) = Element % NodeIndexes(2)-1
             Surfaces(2*(t-1)+2) = Element % NodeIndexes(1)-1
           ELSE
             Surfaces(2*(t-1)+1) = Element % NodeIndexes(1)-1
             Surfaces(2*(t-1)+2) = Element % NodeIndexes(2)-1
           END IF
         ELSE
           j = Element % TYPE % ElementCode/100
           SELECT CASE(j)
             CASE(2)
                TYPE(t) = 202
             CASE(3)
                TYPE(t) = 303
             CASE(4)
                TYPE(t) = 404
           END SELECT
           DO i=1,j
             Surfaces(j*(t-1)+i) = Element % NodeIndexes(i)-1
           END DO
           
           IF (Normal_in*(Nrm(1)*Nx + Nrm(2)*Ny + Nrm(3)*nz)>0) THEN
             Normals(3*(t-1)+1:3*(t-1)+3) =  Nrm
           ELSE
             Normals(3*(t-1)+1:3*(t-1)+3) = -Nrm
           END IF
         END IF
       END DO

       CALL Info( Caller, 'Computing radiator view coefficients...', Level=4 )

       at0 = CPUTime(); rt0 = RealTime()

       Combine = GetLogical( Params, 'Viewfactor combine elements',GotIt)
       IF ( .NOT. GotIt ) Combine = .TRUE.
       IF( Combine ) THEN
         CombineInt = 1
       ELSE
         CombineInt = 0
       END IF
       
       AreaEPS = GetConstReal( Params, 'Viewfactor Area Tolerance',  GotIt )
       IF ( .NOT. GotIt ) AreaEPS = 1.0d-1

       FactEPS = GetConstReal( Params, 'Viewfactor Factor Tolerance ', GotIt )
       IF ( .NOT. GotIt ) FactEPS = 1.0d-2

       RayEPS = GetConstReal( Params, 'Viewfactor Raytrace Tolerace',  GotIt )
       IF ( .NOT. GotIt ) RayEPS = 1.0d-5

       Nrays = GetInteger( Params, 'Viewfactor Number of Rays ',  GotIt )
       IF ( .NOT. GotIt ) Nrays = 1

       CALL RadiatorFactors3d( N, Surfaces, Type, Coords, Normals, 0, Surfaces, Type, Coords, Normals, &
           NofRadiators, Radiators, Factors, AreaEPS, FactEPS, RayEPS, Nrays, 4, 3, CombineInt )
       
       WRITE (Message,'(A,F8.2,F8.2)') 'Radiator factors computed in time (s):',CPUTime()-at0, RealTime()-rt0
       CALL Info( Caller,Message, Level=3 )

       IF(SYMMETRY_NOW) THEN
         DO l=6,1,-1
           SELECT CASE(l)
           CASE(1)
             Plane = GetCReal( Params, 'Viewfactor Symmetry x min', Found );
             IF(.NOT. Found ) Found = ListGetLogical( Params, 'Viewfactor Symmetry x', GotIt );
             ! Note that Plane is zero if 1st keyword not found!
           CASE(2)
             Plane = GetCReal( Params, 'Viewfactor Symmetry x max', Found );
           CASE(3)
             Plane = GetCReal( Params, 'Viewfactor Symmetry y min', Found );
             IF(.NOT. Found ) Found = ListGetLogical( Params, 'Viewfactor Symmetry y', GotIt );
           CASE(4)
             Plane = GetCReal( Params, 'Viewfactor Symmetry y max', Found );
           CASE(5)
             Plane = GetCReal( Params, 'Viewfactor Symmetry z min', Found );
             IF(.NOT. Found ) Found = ListGetLogical( Params, 'Viewfactor Symmetry z', GotIt );
           CASE(6)
             Plane = GetCReal( Params, 'Viewfactor Symmetry z max', Found );
           END SELECT

           IF(.NOT.Found) CYCLE

           CALL Info(Caller,'Symmetry reduction in coordinate direction: '//I2S((l+1)/2))
           
           k = 0
           DO i=1,NofRadiators
             DO j=1,n/2
               k = k + 1
               Factors(k) = Factors((i-1)*n+j) + Factors((i-1)*n+j+n/2)
             END DO
           END DO
           n = n/2
         END DO
       END IF
       

       nprob = 0
       DO i=1,NofRadiators
         s = 0.0_dp
         DO j=1,N
           IF(Factors((i-1)*N+j) < MinFactor) Factors((i-1)*N+j) = 0.0d0         
           s = s + Factors((i-1)*N+j)
         END DO
         
	 IF( .NOT. RadiationOpen .AND. s < 0.5 ) nprob = nprob + 1

         IF(i == 1) THEN
           Fmin = s 
           Fmax = s
         ELSE         
           FMin = MIN( FMin,s )
           FMax = MAX( FMax,s )
         END IF
       END DO       
       
       CALL Info( Caller, ' ', Level=3 )
       CALL info( Caller, 'Radiator factors before manipulation: ', Level=3 )
       WRITE( Message,'(A,ES12.3)') 'Minimum row sum: ',FMin
       CALL Info( Caller, Message )
       WRITE( Message,'(A,ES12.3)') 'Maximum row sum: ',Fmax
       CALL Info( Caller, Message )
       IF(nprob>0) CALL info( Caller, 'Number of rowsums below 0.5 is: '&
           //I2S(nprob)//' (out of '//I2S(n)//')')
       
       IF( .NOT. RadiationOpen ) THEN 
         at0 = CPUTime()

         CALL Info( Caller,'Normalizing Factors...')

         DO i=1,NofRadiators
           s = 0.0_dp
           DO j=1,N
             s = s + Factors((i-1)*N+j)
           END DO

           DO j=1,N
             Factors((i-1)*N+j) = Factors((i-1)*N+j) / s
           END DO
         END DO

         WRITE (Message,'(A,F8.2)') 'Radiator factors manipulated in time (s):',CPUTime()-at0
         CALL Info( Caller,Message, Level=3 )

         CALL Info( Caller, ' ', Level=3 )
         CALL info( Caller, 'Radiator factors after manipulation: ')
         WRITE( Message,'(A,ES12.3)') 'Minimum row sum: ',FMin
         CALL Info( Caller, Message )
         WRITE( Message,'(A,ES12.3)') 'Maximum row sum: ',Fmax
         CALL Info( Caller, Message )
         IF( FMax > 1.001 ) THEN
           CALL Warn(Caller,'Rowsum of view factors should not be larger than one!')
         END IF
         IF( FMin < 0.999 ) THEN
           CALL Warn(Caller,'Rowsum of view factors should not be smaller than one!')
         END IF
       END IF

       ViewFactorsFile = GetString( GetSimulation(),'Raditor Factors',GotIt)
       IF ( .NOT.GotIt ) ViewFactorsFile = 'RadiatorFactors.dat'

       IF(RadiationBody > 1) THEN
         TempString = ViewFactorsFile
         WRITE(ViewFactorsFile, '(A,I1)') TRIM(TempString),RadiationBody
       END IF
       
       IF ( LEN_TRIM(Model % Mesh % Name) > 0 ) THEN
         OutputName = TRIM(OutputPath) // '/' // &
             TRIM(Model % Mesh % Name) // '/' // TRIM(ViewFactorsFile)
       ELSE
         OutputName = TRIM(ViewFactorsFile)
       END IF
       
       BLOCK
         LOGICAL :: BinaryMode
         LOGICAL, ALLOCATABLE :: SaveMask(:)
         ALLOCATE( SaveMask(SIZE(Factors)), STAT=istat)
         IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error for SaveMask.')
         
         ! Use loser constraint for MinFactor as the errors can't be renormalized any more 
         MinFactor = MinFactor / 10.0
         
         BinaryMode = ListGetLogical( Params,'Radiatorfactor Binary Output',Found ) 
         IF(.NOT. Found) BinaryMode = ListGetLogical( Params,'Viewfactor Binary Output',Found ) 
         
         SaveMask = ( Factors(1:n*NofRadiators) > MinFactor )

         IF( BinaryMode ) THEN
           CALL Info(Caller,'Saving radiator factors in binary mode',Level=5)

           OPEN( UNIT=VFUnit, FILE=TRIM(OutputName), FORM = 'unformatted', &
               ACCESS = 'stream', STATUS='replace', ACTION='write' )
           
           WRITE( VFUnit ) N

           DO i=1,NofRadiators
             k = COUNT( SaveMask((i-1)*n+1:i*n) )
             WRITE( VFUnit ) k 
             DO j=1,n
               IF( SaveMask((i-1)*N+j ) ) THEN
                 WRITE( VFUnit ) j,Factors((i-1)*n+j)
               END IF
             END DO
           END DO
         ELSE
           CALL Info(Caller,'Saving radiator factors in ascii mode',Level=5)

           OPEN( UNIT=VFUnit, FILE=TRIM(OutputName), STATUS='unknown' )
           
           DO i=1,NofRadiators
             k = COUNT( SaveMask((i-1)*n+1:i*n) )
             WRITE( VFUnit,* ) k
             DO j=1,n
               IF ( SaveMask((i-1)*n+j) ) THEN
                 WRITE( VFUnit,* ) i,j,Factors((i-1)*n+j)
               END IF
             END DO
           END DO
         END IF
           
         CLOSE(VFUnit)

         DEALLOCATE( SaveMask ) 
         
       END BLOCK
         
       IF ( CylindricSymmetry ) THEN
         DEALLOCATE( Surfaces, Factors)
       ELSE
         DEALLOCATE( Normals, Factors, Surfaces, TYPE)
       END IF
       
     END DO  ! Of radiation RadiationBody

     CALL Info( Caller, '*** ALL DONE ***' )
     CALL FLUSH(6)

   CONTAINS
   


    SUBROUTINE MirrorMesh(Mesh,c,Plane )
       IMPLICIT NONE

       TYPE(Mesh_t) :: Mesh
       INTEGER :: c
       REAL(KIND=dp) :: Plane


       TYPE(Element_t), POINTER :: el(:)
       INTEGER :: i,j,ne, nd, nn, nb, nv
       REAL(KIND=dp), POINTER :: ox(:), oy(:), oz(:)

       nv  = Mesh % NumberOfBulkElements
       nb  = Mesh % NumberOfBoundaryElements
       ne  = nb+nv

       nd = Mesh % NumberOfNodes
       ox => Mesh % Nodes % x
       oy => Mesh % Nodes % y
       oz => Mesh % Nodes % z
       ALLOCATE(Mesh % Nodes % x(2*nd), Mesh % Nodes % y(2*nd), Mesh % Nodes % z(2*nd), STAT=istat )
       IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation for MirroMesh nodes.')

       
       DO i=1,nd
         Mesh % Nodes % x(i) = ox(i)
         Mesh % Nodes % y(i) = oy(i)
         Mesh % Nodes % z(i) = oz(i)
         SELECT CASE(c)
         CASE(1,2)
           Mesh % Nodes % x(i+nd) = 2*Plane - ox(i)
           Mesh % Nodes % y(i+nd) = oy(i)
           Mesh % Nodes % z(i+nd) = oz(i)
         CASE(3,4)
           Mesh % Nodes % x(i+nd) = ox(i)
           Mesh % Nodes % y(i+nd) = 2*Plane - oy(i)
           Mesh % Nodes % z(i+nd) = oz(i)
         CASE(5,6)
           Mesh % Nodes % x(i+nd) = ox(i)
           Mesh % Nodes % y(i+nd) = oy(i)
           Mesh % Nodes % z(i+nd) = 2*Plane - oz(i)
         END SELECT
       END DO

       el => Mesh % Elements
       ALLOCATE(Mesh % Elements(2*ne), STAT=istat)
       IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation for MirroMesh elements.')
       
       DO i=1,nv
         Mesh % Elements(i)    = el(i)
         Mesh % Elements(i+nv) = el(i)
         nn = el(i) % Type % NumberOfNodes
         ALLOCATE(Mesh % Elements(i+nv) % NodeIndexes(nn),STAT=istat)
         IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation for MirroMesh node indexes.')
         
         Mesh % Elements(i+nv) % NodeIndexes = el(i) % NodeIndexes+nd
         Mesh % Elements(i+nv) % ElementIndex = i+nv
       END DO

       DO i=nv+1,ne
         j = i+nv
         Mesh % Elements(j)    = el(i)
         Mesh % Elements(j+nb) = el(i)
         nn = el(i) % TYPE % NumberOfNodes
         
         ALLOCATE(Mesh % Elements(j+nb) % NodeIndexes(nn),STAT=istat)
         IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation for MirroMesh NodeIndexes.')         
         Mesh % Elements(j+nb) % NodeIndexes = el(i) % NodeIndexes+nd

         ALLOCATE(Mesh % Elements(j) % BoundaryInfo)
         Mesh % Elements(j) % BoundaryInfo    = el(i) % BoundaryInfo

         ALLOCATE(Mesh % Elements(j+nb) % BoundaryInfo)
         Mesh % Elements(j+nb) % BoundaryInfo = el(i) % BoundaryInfo

         IF(ASSOCIATED(Mesh % Elements(j) % BoundaryInfo % Left)) THEN
           Mesh % Elements(j+nb) % BoundaryInfo % Left => &
               Mesh % Elements(el(i) % BoundaryInfo % Left % ElementIndex+nv)
         END IF

         IF(ASSOCIATED(Mesh % Elements(j) % BoundaryInfo % Right)) THEN
           Mesh % Elements(j+nb) % BoundaryInfo % Right => &
               Mesh % Elements(el(i) % BoundaryInfo % Right % ElementIndex+nv)
         END IF
       END DO

       DEALLOCATE(ox,oy,oz)

       Mesh % NumberOfNodes = 2*nd
       Mesh % NumberOfBulkElements = 2*Mesh % NumberOfBulkElements
       Mesh % NumberOfBoundaryElements = 2*Mesh % NumberOfBoundaryElements

       Model % NumberOfNodes = Mesh % NUmberOfNodes
       Model % NumberOfBulkElements = Mesh % NumberOfBulkElements
       Model % NumberOfBoundaryElements = Mesh % NumberOfBoundaryElements
     END SUBROUTINE MirrorMesh

  END PROGRAM RadiatorFactors

  
!> \}
!> \}  
