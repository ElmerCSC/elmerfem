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

!> \defgroup ViewFactors Program ViewFactors
!> \{

   MODULE ViewFactorGlobals
     USE Types
     REAL(KIND=dp), ALLOCATABLE :: Jdiag(:), Jacobian(:,:)
   END MODULE ViewFactorGlobals

!------------------------------------------------------------------------------
!> A separate program that computes the view factors to an external file. 
!> This file is later used within the ElmerSolver. If the view factor files
!> do not exist, a system call for this program is performed. 
!------------------------------------------------------------------------------
      
   PROGRAM ViewFactors
   
     USE DefUtils
     USE ViewFactorGlobals

     IMPLICIT NONE

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Model_t), POINTER :: Model
     TYPE(Mesh_t), POINTER  :: Mesh
     TYPE(Solver_t), POINTER  :: Solver

     INTEGER :: i,j,k,l,t,k1,k2,n,iter,Ndeg,Time,NSDOFs,MatId,istat

     REAL(KIND=dp) :: SimulationTime,dt,s,a1,a2,FMin,FMax,Fave

     INTEGER, ALLOCATABLE ::  Surfaces(:), TYPE(:)
     REAL(KIND=dp), ALLOCATABLE :: Coords(:),Normals(:),Factors(:)

     TYPE(Element_t),POINTER :: Element, LParent, RParent

     INTEGER :: BandSize,SubbandSize,RadiationSurfaces,Row,Col
     INTEGER, DIMENSION(:), POINTER :: Perm

     REAL(KIND=dp) :: Norm,PrevNorm,MinFactor, Normal_in, Plane
 
     TYPE(Nodes_t) :: ElementNodes
     TYPE(ValueList_t), POINTER :: BC, Material, Params

     INTEGER :: LeftNode,RightNode,LeftBody,RightBody,RadBody
     REAL(KIND=dp) :: NX,NY,NZ,NRM(3),NrmB(3),DensL,DensR
     INTEGER, ALLOCATABLE :: VF_cohorts(:)
     
     INTEGER :: divide, nprob
     REAL(KIND=dp) :: AreaEPS, RayEPS, FactEPS
     REAL(KIND=dp) :: at0, rt0
     CHARACTER(*), PARAMETER :: Caller = 'ViewFactors'

     INTERFACE
        ! void viewfactors3d
        ! ( int *EL_N,  int *EL_Topo, int *EL_Type, double *EL_Coord, double *EL_Normals,
        ! int *RT_N0, int *RT_Topo0, int *RT_Type, double *RT_Coord, double *RT_Normals,
        ! double *Factors, double *Feps, double *Aeps, double *Reps, int *Nr, 
        ! int *NInteg,int *NInteg3, int  *Combine )
        SUBROUTINE viewfactors3d(EL_N,  EL_Topo, EL_Type, EL_Coord, EL_Normals, &
                                 RT_N0, RT_Topo0, RT_Type, RT_Coord, RT_Normals, &
                                 Factors, Feps, Aeps, Reps, Nr, NInteg, NInteg3, Combine) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT) :: EL_N, EL_Topo(*), EL_Type(*)
            REAL(KIND=C_DOUBLE) :: EL_Coord(*), EL_Normals(*)
            INTEGER(C_INT) :: RT_N0, RT_Topo0(*), RT_Type(*)
            REAL(KIND=C_DOUBLE) :: RT_Coord(*), RT_Normals(*), Factors(*)
            REAL(KIND=C_DOUBLE) :: Feps, Aeps, Reps
            INTEGER(C_INT) :: Nr, NInteg, NInteg3, Combine
        END SUBROUTINE viewfactors3d
        
        ! extern "C" void STDCALLBULL viewfactorsaxis
        ! (int *n,int *surf, Real *crd, Real *vf, int *idiv, int *fast)
        ! CALL ViewFactorsAxis( N, Surfaces, Coords, Factors, divide, CombineInt )
        SUBROUTINE viewfactorsaxis(n, surf, crd, vf, idiv, fast) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE 
            INTEGER(C_INT) :: n, surf(*)
            REAL(KIND=C_DOUBLE) :: crd(*), vf(*)
            INTEGER :: idiv, fast 
        END SUBROUTINE viewfactorsaxis

      END INTERFACE

     INTEGER, POINTER :: Timesteps(:)
     INTEGER :: TimeIntervals,interval,timestep,combineInt
     
     LOGICAL :: CylindricSymmetry,GotIt, Found, Radiation, LeftEmis, RightEmis

     CHARACTER(:), ALLOCATABLE :: eq, RadiationFlag, ViewFactorsFile, OutputName, &
                  LMessage, TempString
     CHARACTER(LEN=MAX_NAME_LEN) :: ModelName

     TYPE(Element_t), POINTER :: RadElements(:)
     INTEGER :: RadiationBody, MaxRadiationBody, Nrays
     LOGICAL :: RadiationOpen, Combine
     INTEGER, PARAMETER :: VFUnit = 10
     INTEGER :: iostat, NoArgs
     
     EXTERNAL MatvecViewFact,DiagPrecViewFact



     CALL Info( Caller, ' ', Level=3 )
     CALL Info( Caller, '==================================================', Level=3 )
     CALL Info( Caller, ' E L M E R  V I E W F A C T O R S,  W E L C O M E',  Level=3  )
     CALL Info( Caller, '==================================================', Level=3 )

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
     CALL Info(Caller,'Computing view factors as defined in file: '//TRIM(ModelName),Level=5)
     
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
     CASE('axi symmetric')
       Coordinates = AxisSymmetric
     CASE('cylindric symmetric')
       Coordinates = CylindricSymmetric
     CASE DEFAULT
       CALL Error( Caller, &
         'Unknown Global Coordinate System for Viewfactor computation ')
       CALL Error( Caller, TRIM(eq) )
       CALL Fatal( Caller, &
         'Only Cartesian 3D or Axi/Cylindrical Symmetric coordinates allowed. Aborting' )
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
       CALL Fatal('Viewfactors', 'Memory allocation error. Aborting' )
     END IF

     ! The routine originally plays with the element list and therefore
     ! when several radiation boundaries are needed both the original and
     ! the new elementlist needs to be in the memory. Thus the hassle.

     MinFactor = ListGetConstReal(Params,'Minimum View Factor',GotIt)
     IF(.NOT. GotIt) MinFactor = 1.0d-20

     CALL AllocateVector( RadElements, Mesh % NumberOfBoundaryElements, Caller )

     IF( Mesh % NumberOfBoundaryElements == 0) THEN
       CALL Warn(Caller,'There are no boundary elements at all!')
       STOP
     END IF

! Check the maximum radiation body
     MaxRadiationBody = 0

     DO t= 1, Mesh % NumberOfBoundaryElements

       Element => GetBoundaryElement(t)
       IF ( GetElementFamily() == 1 ) CYCLE
       BC => GetBC()
       IF ( .NOT. ASSOCIATED( BC ) ) CYCLE
 
       RadiationFlag = GetString( BC, 'Radiation',GotIt )

       IF ( GotIt .AND. RadiationFlag == 'diffuse gray' ) THEN
         i = MAX(1, GetInteger( BC, 'Radiation Boundary', GotIt ) )
         MaxRadiationBody = MAX(i, MaxRadiationBody)
       END IF
     END DO

     IF( Mesh % NumberOfBoundaryElements == 0) THEN
       CALL Warn(Caller,'There are no radiation boundary elements!')
       STOP
     END IF

     RadiationBody = 0
     DO RadiationBody = 1, MaxRadiationBody
       LMessage = 'Computing view factors for radiation body'//I2S(RadiationBody)
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
         
         RadiationFlag = GetString( BC, 'Radiation', GotIt )
         IF ( GotIt .AND. RadiationFlag == 'diffuse gray' ) THEN
           i = MAX(1, GetInteger( BC, 'Radiation Boundary', GotIt ))
           IF(i == RadiationBody) THEN
             RadiationOpen = RadiationOpen .OR. GetLogical( BC, 'Radiation Boundary Open', GotIt )
             RadiationSurfaces = RadiationSurfaces + 1
             j = t + Mesh % NumberOFBulkElements
             RadElements(RadiationSurfaces) = Mesh % Elements(j)
           END IF
         END IF
       END DO
       
       
       N = RadiationSurfaces
       
       IF ( N == 0 ) THEN
         CALL Warn( 'Viewfactors', 'No surfaces participating in radiation?' )
         IF(RadiationBody < MaxRadiationBody) THEN
           CYCLE
         ELSE
           CALL Warn( 'Viewfactors', 'Stopping cause nothing to be done...' )
           STOP
         END IF
       END IF
       
       CALL Info(Caller,'Number of surfaces participating in radiation: '//I2S(N))
       
       IF ( CylindricSymmetry ) THEN
         ALLOCATE( Surfaces(2*N), Factors(N*N), STAT=istat )
       ELSE
         ALLOCATE( Normals(3*N), Factors(N*N),Surfaces(4*N), TYPE(N), STAT=istat )
       END IF
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'Viewfactors', 'Memory allocation error. Aborting' )
       END IF
       
       DO t=1,N

         Element => RadElements(t)
         Model % CurrentElement => Element
         k = GetElementNOFNodes()
         CALL GetElementNodes(ElementNodes)
                  
         nrm = NormalVector( Element, ElementNodes )
         
         LeftBody = 0
         LeftNode = -1
	 LeftEmis = .FALSE.
         LParent => Element % BoundaryInfo % Left
         IF ( ASSOCIATED(LParent) ) THEN
           LeftBody  = LParent % BodyId
           k1 = 0
           DO i=1,LParent % TYPE % NumberOfNodes
             gotIt =.TRUE.
             ! Find a node in parent element that is not a node in boundary element.
             DO j=1,Element % TYPE % NumberOfNodes
               IF ( Element % NodeIndexes(j) == LParent % NodeIndexes(i)) THEN
                 k1 = k1 + 1
                 gotIt=.FALSE.
                 EXIT
               END IF
             END DO
             IF (gotIt) THEN
               LeftNode = LParent % NodeIndexes(i)
             END IF
           END DO
           IF(k1 /= Element % TYPE % NumberOfNodes ) THEN
             CALL Fatal(Caller,'Boundary element '//I2S(Element % ElementIndex - Mesh % NumberOfBulkElements)//&
                 ' not included in left parent '//I2S(LParent % ElementIndex)//'!')               
           END IF
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
         RParent => Element % BoundaryInfo % Right
         IF ( ASSOCIATED(RParent) ) THEN
           RightBody = RParent % BodyId
           k1 = 0
           DO i=1,RParent % TYPE % NumberOfNodes
             gotIt =.TRUE.
             DO j=1,Element % TYPE % NumberOfNodes
               IF (Element % NodeIndexes(j) == RParent % NodeIndexes(i)) THEN
                 k1 = k1 + 1
                 gotIt=.FALSE.
                 EXIT
               END IF
             END DO
             IF (gotIt) THEN
               RightNode = RParent % NodeIndexes(i)
             END IF
           END DO
           IF(k1 /= Element % TYPE % NumberOfNodes ) THEN
             CALL Fatal(Caller,'Boundary element '//I2S(Element % ElementIndex - Mesh % NumberOfBulkElements)//&
                 ' not included in right parent '//I2S(RParent % ElementIndex)//'!')               
           END IF
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
             CALL Fatal(Caller,'Cannot continue!')
           ELSE IF( RightEmis ) THEN
             RadBody = LeftBody
           ELSE IF( LeftEmis ) THEN
             RadBody = RightBody
           END IF
         END IF

         IF ( RadBody < 0 ) RadBody = 0
         
         IF ( RadBody>0 .AND. (RadBody /= RightBody .AND. RadBody /= LeftBody) ) THEN
           CALL Error( Caller, 'Inconsistent direction information (Radiation Target Body)' )
           LMessage = 'Radiation Target: '//I2S(RadBody)//' Left, Right: '//&
                        I2S(LeftBody)//I2S(RightBody)
           CALL Fatal( Caller, LMessage )
         END IF
         
         IF ( LeftNode <= 0 .OR. (RadBody>0 .AND. RadBody==RightBody) ) THEN
           LeftNode = RightNode
         END IF
         Normal_in = -1.0
         IF ( RadBody <= 0 ) Normal_in = 1.0
                  
         BLOCK
           REAL(KIND=dp) :: r1(3), r2(3)           
           r1(1) = SUM(ElementNodes % x)/k
           r1(2) = SUM(ElementNodes % y)/k
           r1(3) = SUM(ElementNodes % z)/k

           r2(1) = Mesh % Nodes % x(LeftNode)
           r2(2) = Mesh % Nodes % y(LeftNode)
           r2(3) = Mesh % Nodes % z(LeftNode)

           ! Direction to test the normal
           NrmB = r1 - r2
         END BLOCK
         
         IF ( CylindricSymmetry ) THEN
           IF ( Normal_in * SUM(Nrm*NrmB) > 0 ) THEN
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
           
           IF (Normal_in * SUM(Nrm*NrmB)>0) THEN
             Normals(3*(t-1)+1:3*(t-1)+3) = Nrm
           ELSE
             Normals(3*(t-1)+1:3*(t-1)+3) = -Nrm
           END IF
         END IF
       END DO

       CALL Info( Caller, 'Computing viewfactors...', Level=4 )

       at0 = CPUTime(); rt0 = RealTime()
       
       Combine = GetLogical( Params, 'Viewfactor combine elements',GotIt)
       IF ( .NOT. GotIt ) Combine = .TRUE.
       IF( Combine ) THEN
         CombineInt = 1
       ELSE
         CombineInt = 0
       END IF

       IF ( CylindricSymmetry ) THEN
         divide = GetInteger( Params, 'Viewfactor divide',GotIt)
         IF ( .NOT. GotIt ) Divide = 1
         CALL ViewFactorsAxis( N, Surfaces, Coords, Factors, divide, CombineInt )
       ELSE
         AreaEPS = GetConstReal( Params, 'Viewfactor Area Tolerance',  GotIt )
         IF ( .NOT. GotIt ) AreaEPS = 1.0d-1
         FactEPS = GetConstReal( Params, 'Viewfactor Factor Tolerance ', GotIt )
         IF ( .NOT. GotIt ) FactEPS = 1.0d-2
         RayEPS = GetConstReal( Params, 'Viewfactor Raytrace Tolerace',  GotIt )
         IF ( .NOT. GotIt ) RayEPS = 1.0d-5
         Nrays = GetInteger( Params, 'Viewfactor Number of Rays ',  GotIt )
         IF ( .NOT. GotIt ) Nrays = 1

         CALL ViewFactors3D( &
             N, Surfaces, Type, Coords, Normals, &
             0, Surfaces, Type, Coords, Normals, &
             Factors, AreaEPS, FactEPS, RayEPS, Nrays, 4, 3, CombineInt )
       END IF
       
       WRITE (Message,'(A,F8.2,F8.2)') 'View factors computed in time (s):',CPUTime()-at0, RealTime()-rt0
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
           DO i=1,n/2
           DO j=1,n/2
             k = k + 1
             Factors(k) = Factors((i-1)*n+j) + Factors((i-1)*n+j+n/2)
           END DO
           END DO
           n = n/2
         END DO
       END IF
       
#if 0
       write(1,*) 2*n,n, 1,4,'vector: nrm scalar:view'
       do i=1,n
         Element => RadElements(i)
         write(1,*) model % nodes % x(element % nodeindexes(1)), mesh % nodes % y(element % nodeindexes(1)),0
         write(1,*) model % nodes % x(element % nodeindexes(2)), mesh % nodes % y(element % nodeindexes(2)),0
       end do
       do i=1,n
         Element => RadElements(i)
         write(1,*) element % boundaryinfo % constraint, ' 202 ', 2*(i-1),2*(i-1)+1
       end do
       do i=1,n
         write(1,*) normals(3*(i-1)+1:3*(i-1)+3), factors((3-1)*n+i)
         write(1,*) normals(3*(i-1)+1:3*(i-1)+3), factors((3-1)*n+i)
       end do
#endif

       nprob = 0
       Fave = 0.0_dp
       IF(.NOT. ALLOCATED(VF_cohorts)) THEN
         ALLOCATE(VF_cohorts(100))
       END IF
       VF_cohorts = 0

       k = 0
       DO i=1,N
         s = 0.0_dp
         DO j=1,N
           IF(Factors((i-1)*N+j) < MinFactor) Factors((i-1)*N+j) = 0.0d0         
           s = s + Factors((i-1)*N+j)
         END DO
         
	 IF( .NOT. RadiationOpen .AND. s < 0.5 ) nprob = nprob + 1

         IF(i == 1) THEN
           Fmin = s 
           Fmax = s
           k = 1
         ELSE
           IF(s < Fmin) THEN
             Fmin = s
             k = i
           END IF
           FMax = MAX( FMax,s )
         END IF

         j = CEILING(100*s)
         j = MIN(100,MAX(1,j))
         VF_cohorts(j) = VF_cohorts(j) + 1
         
         Fave = Fave + s         
       END DO
       Fave = Fave / N
       
       CALL Info( Caller, ' ', Level=3 )
       CALL info( Caller, 'Viewfactors before manipulation: ', Level=3 )
       WRITE( Message,'(A,ES14.6)') 'Minimum row sum: ',FMin
       CALL Info( Caller, Message )
       WRITE( Message,'(A,ES14.6)') 'Maximum row sum: ',Fmax
       CALL Info( Caller, Message )
       WRITE( Message,'(A,ES14.6)') 'Average row sum: ',Fave
       CALL Info( Caller, Message )
       IF(nprob>0) CALL info( Caller, 'Number of rowsums below 0.5 is: '&
           //I2S(nprob)//' (out of '//I2S(n)//')')
       
       IF( InfoActive(10) ) THEN
         ! Report on the most problematic element which has too small viewfactors.
         IF( Fmin < 0.5 .AND. .NOT. RadiationOpen ) THEN
           Element => RadElements(k)
           j = Element % ElementIndex - Model % Mesh % NumberOfBulkElements
           CALL Info( Caller,'Location of minimum rowsum '//I2S(k)//' element '//I2S(j))
           PRINT *,'Indexes:',element % nodeindexes
           PRINT *,'X coord:',model % nodes % x(element % nodeindexes)
           PRINT *,'Y coord:',model % nodes % y(element % nodeindexes)
           PRINT *,'Z coord:',model % nodes % z(element % nodeindexes)
           DO i=1,3
             j = element % nodeindexes(i)
             PRINT *,'r:',sqrt(model % nodes % x(j)**2 + model % nodes % y(j)**2 + model % nodes % z(j)**2)
           END DO

         END IF
         
         DO i=1,100
           j = VF_cohorts(i)
           IF(j==0) CYCLE
           CALL Info(Caller,'VF percentile '//I2S(i-1)//'-'//I2S(i)//' count: '//I2S(j))
         END DO
       END IF

       
       at0 = CPUTime()

       IF( RadiationOpen ) THEN
         CALL Info( Caller,'Symmetrizing Factors... ')
       ELSE
         CALL Info( Caller,'Normalizaing Factors...')
         IF( Fmin < EPSILON( Fmin ) ) THEN
           CALL Fatal(Caller,'Invalid view factors for normalization, check your geometry!')
         END IF       
       END IF

       CALL NormalizeFactors( Model )

       DO i=1,N
         s = 0.0D0
         DO j=1,N
           s = s + Factors((i-1)*N+j)
         END DO
         IF(i == 1) THEN
           Fmin = s; Fmax = s
         ELSE
           FMin = MIN(FMin,s)
           FMax = MAX(FMax,s)
         END IF
       END DO
       
       WRITE (Message,'(A,F8.2)') 'View factors manipulated in time (s):',CPUTime()-at0
       CALL Info( Caller,Message, Level=3 )
       
       CALL Info( Caller, ' ', Level=3 )
       CALL info( Caller, 'Viewfactors after manipulation: ')
       WRITE( Message,'(A,ES12.3)') 'Minimum row sum: ',FMin
       CALL Info( Caller, Message )
       WRITE( Message,'(A,ES12.3)') 'Maximum row sum: ',Fmax
       CALL Info( Caller, Message )
       IF( FMax > 1.001 ) THEN
         CALL Warn(Caller,'Rowsum of view factors should not be larger than one!')
       END IF
       IF( FMin < 0.999 ) THEN
         ! For open BCs the view factor sum may be much less than one, otherwise not.
         IF(.NOT. ListCheckPresentAnyBC( Model,'Radiation Boundary Open') ) THEN          
           CALL Warn(Caller,'Rowsum of view factors should not be smaller than one!')
         END IF
       END IF

       IF(InfoActive(7)) THEN
         CALL ViewFactorsLumping()
       END IF
         
       ViewFactorsFile = GetString( GetSimulation(),'View Factors',GotIt)
       IF ( .NOT.GotIt ) ViewFactorsFile = 'ViewFactors.dat'
       IF(RadiationBody > 1) THEN
         TempString = ViewFactorsFile
         ViewFactorsFile = TRIM(TempString)//I2S(RadiationBody)
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
         
         BinaryMode = ListGetLogical( Params,'Viewfactor Binary Output',Found ) 
         
         SaveMask = ( Factors > MinFactor )

         IF( BinaryMode ) THEN
           CALL Info(Caller,'Saving view factors in binary mode',Level=5)

           OPEN( UNIT=VFUnit, FILE=TRIM(OutputName), FORM = 'unformatted', &
               ACCESS = 'stream', STATUS='replace', ACTION='write' )
           
           WRITE( VFUnit ) N

           DO i=1,N
             k = COUNT( SaveMask((i-1)*N+1:i*N) )
             WRITE( VFUnit ) k 
             DO j=1,N
               IF( SaveMask((i-1)*N+j ) ) THEN
                 WRITE( VFUnit ) j,Factors((i-1)*N+j)
               END IF
             END DO
           END DO           
         ELSE
           CALL Info(Caller,'Saving view factors in ascii mode',Level=5)

           OPEN( UNIT=VFUnit, FILE=TRIM(OutputName), STATUS='unknown' )
           
           DO i=1,N
             k = COUNT( SaveMask((i-1)*N+1:i*N) )
             WRITE( VFUnit,* ) k
             DO j=1,N
               IF ( SaveMask((i-1)*N+j) ) THEN
                 WRITE( VFUnit,* ) i,j,Factors((i-1)*N+j)
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


     ! Provide useful information on the boundary-to-boundary view factors that is
     ! obtained as the area-weigted average of elemental view factor sums.
     !-----------------------------------------------------------------------------
     SUBROUTINE ViewFactorsLumping()

       INTEGER :: i,j,k,m,MaxRadBC,bc_id
       INTEGER, ALLOCATABLE :: BCNumbering(:), VFPerm(:)
       REAL(KIND=dp), ALLOCATABLE :: LumpedVF(:,:), Areas(:), LumpedAreas(:)

       CALL Info(Caller,'Printing some lumped information of view factors')
       
       ALLOCATE(BCNumbering(Model % NumberOfBCs),VFPerm(N),Areas(N))
       BCNumbering = 0
       VFPerm = 0
       Areas = 0.0_dp
              
       DO i=1,N
         Element => RadElements(i)
         Areas(i) = ElementArea( Mesh, Element, Element % TYPE % NumberOfNodes)
         bc_id = GetBCId( Element ) 
         IF(bc_id < 0 .OR. bc_id > Model % NumberOfBCs) THEN
           CALL Warn(Caller,'BC index out of bounds: '//I2S(bc_id))
           CYCLE
         END IF
         BCNumbering(bc_id) = BCNumbering(bc_id) + 1
         VFPerm(i) = bc_id
       END DO

       j = 0
       DO i=1,Model % NumberOfBCs
         m = BCNumbering(i)
         IF(m>0) THEN
           j = j+1
           BCNumbering(i) = j
           CALL Info(Caller,'BC '//I2S(i)//' with '//I2S(m)//' elems perm: '//I2S(j))
         END IF         
       END DO
       MaxRadBC = j

       DO i=1,N
         VFPerm(i) = BCNumbering(VFPerm(i))
       END DO
       
       ALLOCATE(LumpedVF(MaxRadBC,MaxRadBC),LumpedAreas(MaxRadBC))
       LumpedVF = 0.0_dp
       LumpedAreas = 0.0_dp
       
       DO i=1,N
         IF(VFPerm(i) > 0) THEN
           LumpedAreas(VFPerm(i)) = LumpedAreas(VFPerm(i)) + Areas(i)
         END IF
       END DO
         
       CALL Info(Caller,'Lumped areas:')
       WRITE(Message,*) LumpedAreas(:)
       CALL Info(Caller, Message ) 
       
       DO i=1,N
         DO j=1,N
           k = (i-1)*N+j
           LumpedVF(VFPerm(i),VFPerm(j)) = LumpedVF(VFPerm(i),VFPerm(j)) + Areas(i) * Factors(k)
         END DO
       END DO
       DO i=1,MaxRadBC
         DO j=1,MaxRadBC
           LumpedVF(i,j) = LumpedVF(i,j) / LumpedAreas(i)
         END DO
       END DO

       CALL Info(Caller,'Lumped View Factor Matrix:')
       DO i=1,MaxRadBC 
         WRITE(Message,*) LumpedVF(i,:)
         CALL Info(Caller, Message ) 
       END DO
       
     END SUBROUTINE ViewFactorsLumping



     
!> View factors are normalized in order to improve the numberical accuracy. With 
!> normalization it is ensured that all boundary elements see exactly half 
!> space. 
!------------------------------------------------------------------------------

      SUBROUTINE NormalizeFactors( Model )
        IMPLICIT NONE
        TYPE(Model_t), POINTER :: Model
!------------------------------------------------------------------------------
        INTEGER :: itmax,it,i,j,k
        LOGICAL :: li,lj
        REAL(KIND=dp), ALLOCATABLE :: RHS(:),SOL(:),Areas(:),PSOL(:)
        REAL(KIND=dp) :: cum,s,si,sj
        REAL(KIND=dp), PARAMETER :: eps=1.0D-20
        REAL(KIND=dp) :: at1

        itmax = 20
        it = 0
        cum = 0.0_dp
        
        ALLOCATE( Areas(n),STAT=istat )
        IF ( istat /= 0 ) THEN
          CALL Fatal(Caller,'Memory allocation error in NormalizeFactors for Areas.' )
        END IF

!------------------------------------------------------------------------------
!       First force the matrix (before dividing by area) to be symmetric
!------------------------------------------------------------------------------
        DO i=1,n
          Element => RadElements(i)
          Areas(i) = ElementArea( Mesh, Element, Element % Type % NumberOfNodes)
        END DO
        
        DO i=1,n
          DO j=i,n
            si = Areas(i) * Factors((i-1)*n+j)
            sj = Areas(j) * Factors((j-1)*n+i)

            li = (ABS(si) < HUGE(si)) 
            lj = (ABS(sj) < HUGE(sj)) 

            IF(li .AND. lj) THEN 
              s = (si+sj)/2.0
            ELSE IF(li) THEN
              s = si
            ELSE IF(lj) THEN
              s = sj
            ELSE 
              s = 0.0
            END IF

            Factors((i-1)*n+j) = s
            Factors((j-1)*n+i) = s
          END DO
        END DO

!------------------------------------------------------------------------------
!       Next we solve the equation DFD = A by Newton iteration (this is a very
!       well behaved equation (symmetric, diagonal dominant), no need for any
!       tricks...)
!------------------------------------------------------------------------------
        IF(.NOT. RadiationOpen ) THEN
          
          ALLOCATE( RHS(n),SOL(n),PSOL(n),Jdiag(n),Jacobian(n,n),STAT=istat )
          IF ( istat /= 0 ) THEN
            CALL Fatal( Caller,'Memory allocation error in NormalizeFactors for RHS etc.' )
          END IF

          SOL = 1.0_dp
          cum = 1.0_dp
          
          DO it=1,itmax
            DO i=1,n
              cum = 0.0_dp
              DO j=1,n
                cum = cum + Factors((i-1)*n+j) * SOL(j)
              END DO
              cum = cum * SOL(i)
              RHS(i) = Areas(i) - cum
            END DO
            
            cum = SUM( RHS*RHS/Areas ) / n
            
            WRITE (Message,'(A,ES12.3)') &
                'Normalization iteration '//I2S(it)//': ',cum
            CALL Info( Caller,Message, Level=3 )
            
            IF ( cum <= eps ) EXIT
            
            DO i=1,n
              DO j=1,n
                Jacobian(i,j) = Factors((i-1)*n+j) * SOL(i)
              END DO
              DO j=1,n
                Jacobian(i,i) = Jacobian(i,i) + Factors((i-1)*n+j) * SOL(j)
              END DO
              Jdiag(i) = 1._dp / Jacobian(i,i)
            END DO

            PSOL = SOL
            CALL IterSolv( n,SOL,RHS )
            SOL = PSOL + SOL
          END DO
          
!------------------------------------------------------------------------------
!       Normalize the factors and (re)divide by areas
!------------------------------------------------------------------------------
          DO i=1,N
            DO j=1,N
              Factors((i-1)*N+j) = Factors((i-1)*N+j)*SOL(i)*SOL(j)/Areas(i)
            END DO
          END DO
          DEALLOCATE( SOL,RHS,PSOL,Jdiag,Jacobian )

        ELSE
         DO i=1,N
           DO j=1,N
             Factors((i-1)*N+j) = Factors((i-1)*N+j)/Areas(i)
           END DO
         END DO
       END IF

       DEALLOCATE( Areas )

    END SUBROUTINE NormalizeFactors


#include "huti_fdefs.h"

!> Local handle to the iterative methods for linear systems. 
!------------------------------------------------------------------------------
    SUBROUTINE IterSolv( N,x,b )
      IMPLICIT NONE

      INTEGER :: N
      REAL(KIND=dp), DIMENSION(n) :: x,b

      REAL(KIND=dp) :: dpar(50)
      INTEGER :: ipar(50),wsize
      REAL(KIND=dp), ALLOCATABLE :: work(:,:)
      INTEGER(KIND=addrInt) :: iterProc, mvProc, pcondProc, dProc
      INTEGER(KIND=AddrInt) :: AddrFunc
      EXTERNAL :: AddrFunc
!------------------------------------------------------------------------------
      HUTI_NDIM = N
      dProc = 0
    
      ipar = 0
      dpar = 0.0_dp

      HUTI_WRKDIM = HUTI_CG_WORKSIZE
      wsize = HUTI_WRKDIM
          
      HUTI_NDIM     = N
      HUTI_DBUGLVL  = 0
      HUTI_MAXIT    = 100
 
      ALLOCATE( work(N, wsize), STAT=istat )
      IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error for IterSolv work.')

      
      work = 0D0
      HUTI_TOLERANCE = 1.0D-12
      HUTI_MAXTOLERANCE = 1.0d20
      HUTI_INITIALX = HUTI_USERSUPPLIEDX
      HUTI_STOPC = HUTI_TRESID_SCALED_BYB
      
      iterProc  = AddrFunc(HUTI_D_CG)
      mvProc    = AddrFunc(MatvecViewFact)
      pcondProc = AddrFunc(DiagPrecViewFact)
      CALL IterCall( iterProc,x,b,ipar,dpar,work,mvProc,pcondProc, &
                dProc, dProc, dProc, dProc )
          
      DEALLOCATE( work )
    END SUBROUTINE IterSolv 


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

  END PROGRAM ViewFactors


  SUBROUTINE DiagPrecViewFact( u,v,ipar )
    USE ViewFactorGlobals
    IMPLICIT NONE

    REAL(KIND=dp) :: u(*),v(*)
    INTEGER :: ipar(*)

    INTEGER :: n

    n = HUTI_NDIM
    u(1:n) = v(1:n)*Jdiag(1:n)
  END SUBROUTINE DiagPrecViewFact


  SUBROUTINE MatvecViewFact( u,v,ipar )
    USE ViewFactorGlobals
    IMPLICIT NONE

    INTEGER :: ipar(*)
    REAL(KIND=dp) :: u(*),v(*)

    INTEGER :: n

    n = HUTI_NDIM
    CALL DGEMV('N',n,n,1.0_dp,Jacobian,n,u,1,0.0_dp,v,1)
  END SUBROUTINE MatvecViewFact

  
!> \}
!> \}  
