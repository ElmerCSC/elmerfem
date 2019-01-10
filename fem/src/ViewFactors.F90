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
! *  V0.0a ELMER/FEM Viewfactor computation
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

     REAL(KIND=dp) :: SimulationTime,dt,s,a1,a2,FMin,FMax

     INTEGER, ALLOCATABLE ::  Surfaces(:), TYPE(:)
     REAL(KIND=dp), ALLOCATABLE :: Coords(:),Normals(:),Factors(:)

     TYPE(Element_t),POINTER :: Element, Parent

     INTEGER :: BandSize,SubbandSize,RadiationSurfaces,Row,Col
     INTEGER, DIMENSION(:), POINTER :: Perm

     REAL(KIND=dp) :: Norm,PrevNorm,MinFactor, Normal_in
 
     TYPE(Nodes_t) :: ElementNodes
     TYPE(ValueList_t), POINTER :: BC, Material

     INTEGER :: LeftNode,RightNode,LeftBody,RightBody,RadBody
     REAL(KIND=dp) :: NX,NY,NZ,NRM(3),DensL,DensR

     INTEGER :: divide
     REAL(KIND=dp) :: AreaEPS, RayEPS, FactEPS
#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at0

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
#else
     REAL(KIND=dp) :: at0, CPUTime
#endif

     INTEGER, POINTER :: Timesteps(:)
     INTEGER :: TimeIntervals,interval,timestep,combineInt
     
     LOGICAL :: CylindricSymmetry,GotIt, Found, Radiation, LeftEmis, RightEmis

     CHARACTER(LEN=MAX_NAME_LEN) :: eq,RadiationFlag, &
           ViewFactorsFile,OutputName,ModelName,LMessage, TempString

     TYPE(Element_t), POINTER :: RadElements(:)
     INTEGER :: RadiationBody, MaxRadiationBody, Nrays
     LOGICAL :: RadiationOpen, Combine

     EXTERNAL MatvecViewFact,DiagPrecViewFact


     CALL Info( 'ViewFactors', ' ', Level=3 )
     CALL Info( 'ViewFactors', '==================================================', Level=3 )
     CALL Info( 'ViewFactors', ' E L M E R  V I E W F A C T O R S,  W E L C O M E',  Level=3  )
     CALL Info( 'ViewFactors', '==================================================', Level=3 )

!------------------------------------------------------------------------------
!    Read element definition file, and initialize element types
!------------------------------------------------------------------------------
     CALL InitializeElementDescriptions
!------------------------------------------------------------------------------
!    Read Model from Elmer Data Base
!------------------------------------------------------------------------------
     CALL Info( 'ViewFactors', ' ', Level=3 )
     CALL Info( 'ViewFactors', ' ', Level=3 )
     CALL Info( 'ViewFactors', 'Reading Model... ', Level=3 )

!------------------------------------------------------------------------------
     OPEN( 1,file='ELMERSOLVER_STARTINFO', STATUS='OLD', ERR=10 )
     GOTO 20


10   CONTINUE
     CALL Fatal( 'ElmerSolver', 'Unable to find ELMERSOLVER_STARTINFO, cannot execute.' )

20   CONTINUE
       READ(1,'(a)') ModelName
     CLOSE(1)

     Model => LoadModel( ModelName,.FALSE.,1,0 )

     CurrentModel => Model

     NULLIFY( Mesh )
     DO i=1,Model % NumberOfSolvers
       Solver => Model % Solvers(i)
       Radiation = ListGetLogical( Solver % Values, 'Radiation Solver', Found )
       eq = ListGetString( Solver % Values, 'Equation' )
       IF ( Radiation .OR. TRIM(eq) == 'heat equation' ) THEN
         Mesh => Solver % Mesh
         Model % Solver => Solver
         EXIT
       ENDIF
     END DO
  
     IF ( .NOT. ASSOCIATED(Mesh) ) THEN
       CALL Fatal( 'ViewFactors', 'No heat equation definition. ' // &
                  'Cannot compute factors.' )
     END IF
     CALL SetCurrentMesh( Model,Mesh )

     CALL Info( 'ViewFactors', '... Done',Level=3 )
     CALL Info( 'ViewFactors', ' ',Level=3 )

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
       CALL Error( 'ViewFactors', &
         'Unknown Global Coordinate System for Viewfactor computation ')
       CALL Error( 'ViewFactors', TRIM(eq) )
       CALL Fatal( 'ViewFactors', &
         'Only Cartesian 3D or Axi/Cylindrical Symmetric coordinates allowed. Aborting' )
     END SELECT

     CylindricSymmetry = (Coordinates == AxisSymmetric) .OR. (Coordinates==CylindricSymmetric)
!------------------------------------------------------------------------------

     ALLOCATE( ElementNodes % x(Model % MaxElementNodes), &
         ElementNodes % y(Model % MaxElementNodes), &
         ElementNodes % z(Model % MaxElementNodes),STAT=istat )
     
     IF ( CylindricSymmetry ) THEN
       ALLOCATE( Coords(2 * Model % NumberOfNodes), STAT=istat )
       DO i=1,Model % NumberOfNodes
         Coords(2*(i-1)+1) = Model % Nodes % x(i)
         Coords(2*(i-1)+2) = Model % Nodes % y(i)
       END DO
     ELSE
       ALLOCATE( Coords(3 * Model % NumberOfNodes), STAT=istat )
       DO i=1,Model % NumberOfNodes
         Coords(3*(i-1)+1) = Model % Nodes % x(i)
         Coords(3*(i-1)+2) = Model % Nodes % y(i)
         Coords(3*(i-1)+3) = Model % Nodes % z(i)
       END DO
     END IF

     IF ( istat /= 0 ) THEN
       CALL Fatal('Viewfactors', 'Memory allocation error. Aborting' )
     END IF

     ! The routine originally plays with the element list and therefore
     ! when several radiation boundaries are needed both the original and
     ! the new elementlist needs to be in the memory. Thus the hassle.

     MinFactor = ListGetConstReal(Solver % Values,'Minimum View Factor',GotIt)
     IF(.NOT. GotIt) MinFactor = 1.0d-20

     CALL AllocateVector( RadElements, Model % NumberOfBoundaryElements, 'ViewFactors' )


     IF( Model % NumberOfBoundaryElements == 0) THEN
       CALL Warn('ViewFactors','There are no boundary elements at all!')
       STOP
     END IF

! Check the maximum radiation body
     MaxRadiationBody = 0

     DO t= 1, Model % NumberOfBoundaryElements

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

     IF( Model % NumberOfBoundaryElements == 0) THEN
       CALL Warn('ViewFactors','There are no radiation boundary elements!')
       STOP
     END IF

     RadiationBody = 0
     DO RadiationBody = 1, MaxRadiationBody
       WRITE( LMessage,'(A,I2)') 'Computing view factors for radiation body',RadiationBody
       CALL Info('ViewFactors',LMessage,Level=3)
    
!------------------------------------------------------------------------------
!    Here we start...
!------------------------------------------------------------------------------
       RadiationSurfaces = 0
       RadiationOpen = .FALSE.

!------------------------------------------------------------------------------
!    loop to get the surfaces participating in radiation, discard the rest
!    of the elements...
!------------------------------------------------------------------------------

       DO t=1,Model % NumberOfBoundaryElements
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
             j = t + Model % NumberOFBulkElements
             RadElements(RadiationSurfaces) = Model % Elements(j)
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
       
       WRITE( LMessage,'(A,I9)' ) 'Number of surfaces participating in radiation',N
       CALL Info('ViewFactors',LMessage)
       
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
         
         IF ( GetElementFamily() == 3 ) THEN
           nrm = NormalVector( Element,ElementNodes, &
               1.0d0 / 3.0d0, 1.0d0 / 3.0d0 )
         ELSE
           nrm = NormalVector( Element, ElementNodes, 0.0d0, 0.0d0 )
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
               CALL Warn('ViewFactors','Invalid material index in body, perhaps none')
             END IF 
           ELSE
             CALL Warn('ViewFactors','LeftBody not associated')
           END IF
           LeftEmis=ListCheckPresent(Model % Materials(MatId) % Values,'Emissivity') 
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
             CALL Warn('ViewFactors','Emissivity defined on both sides!')
             WRITE(Message,'(A,I3,I3,A,I3,A,I5)') 'Bodies:',RightBody,LeftBody,' BC:', &
                            GetBCId( Element ),' Ind:',t
             CALL Info('ViewFactors',Message)
             IF( ASSOCIATED(Element % BoundaryInfo % Left, Element % BoundaryInfo % Right)) THEN
              CALL Warn('ViewFactors','Parents of the boundary element are the same')
              RadBody = LeftBody
             END IF
           ELSE IF( RightEmis ) THEN
             RadBody = LeftBody
           ELSE IF( LeftEmis ) THEN
             RadBody = RightBody
           END IF
         END IF

         IF ( RadBody < 0 ) RadBody = 0
         
         IF ( RadBody > 0 .AND. (RadBody /= RightBody .AND. RadBody /= LeftBody) ) THEN
           CALL Error( 'ViewFactors', 'Inconsistent direction information (Radiation Target Body)' )
           WRITE( LMessage, * ) 'Radiation Target: ', RadBody, ' Left, Right: ', LeftBody, RightBody
           CALL Fatal( 'ViewFactors', LMessage )
         END IF
         
         IF ( LeftNode <= 0 .OR. (RadBody>0 .AND. RadBody==RightBody) ) &
           LeftNode = RightNode
         Normal_in = -1.0
         IF ( RadBody <= 0 ) Normal_in = 1.0
         
         nx = SUM(ElementNodes % x)/k - Model % Nodes % x(LeftNode)
         ny = SUM(ElementNodes % y)/k - Model % Nodes % y(LeftNode)
         nz = SUM(ElementNodes % z)/k - Model % Nodes % z(LeftNode)
         
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
             Surfaces(j*(t-1)+i) = Element % NodeIndexes(i) - 1
           END DO
           
           IF (Normal_in*(Nrm(1)*Nx + Nrm(2)*Ny + Nrm(3)*nz)>0) THEN
             Normals(3*(t-1)+1:3*(t-1)+3) =  Nrm
           ELSE
             Normals(3*(t-1)+1:3*(t-1)+3) = -Nrm
           END IF
         END IF
         
       END DO
       
       CALL Info( 'ViewFactors', 'Computing viewfactors...', Level=4 )

       at0 = CPUTime()
       
       Combine = GetLogical( GetSolverParams(), 'Viewfactor combine elements',GotIt)
       IF ( .NOT. GotIt ) Combine = .TRUE.
       IF( Combine ) THEN
         CombineInt = 1
       ELSE
         CombineInt = 0
       END IF

       IF ( CylindricSymmetry ) THEN
         divide = GetInteger( GetSolverParams(), 'Viewfactor divide',GotIt)
         IF ( .NOT. GotIt ) Divide = 1
         CALL ViewFactorsAxis( N, Surfaces, Coords, Factors, divide, CombineInt )
       ELSE
         AreaEPS = GetConstReal( GetSolverParams(), 'Viewfactor Area Tolerance',  GotIt )
         IF ( .NOT. GotIt ) AreaEPS = 1.0d-1
         FactEPS = GetConstReal( GetSolverParams(), 'Viewfactor Factor Tolerance ', GotIt )
         IF ( .NOT. GotIt ) FactEPS = 1.0d-2
         RayEPS = GetConstReal( GetSolverParams(), 'Viewfactor Raytrace Tolerace',  GotIt )
         IF ( .NOT. GotIt ) RayEPS = 1.0d-5
         Nrays = GetInteger( GetSolverParams(), 'Viewfactor Number of Rays ',  GotIt )
         IF ( .NOT. GotIt ) Nrays = 1

         CALL ViewFactors3D( &
             N, Surfaces, TYPE, Coords, Normals, &
             0, Surfaces, TYPE, Coords, Normals, &
             Factors, AreaEPS, FactEPS, RayEPS, Nrays, 4, 3, CombineInt )
       END IF
       
       WRITE (Message,'(A,F8.2)') 'View factors computed in time (s):',CPUTime()-at0
       CALL Info( 'ViewFactors',Message, Level=3 )
       
       
       DO i=1,N
         s = 0.0D0
         DO j=1,N
           IF(Factors((i-1)*N+j) < MinFactor) Factors((i-1)*N+j) = 0.0d0         
           s = s + Factors((i-1)*N+j)
         END DO
         
	 IF( .NOT. RadiationOpen .AND. s < 0.1 ) THEN
	  PRINT *,'Problematic row sum',i,s,n
	  j = Surfaces(2*i-1)
	  PRINT *,'coord 1:',j,Coords(2*j-1),Coords(2*j)
	  j = Surfaces(2*i)
	  PRINT *,'coord 2:',j,Coords(2*j-1),Coords(2*j)
         END IF

         IF(i == 1) THEN
           Fmin = s 
           Fmax = s
         ELSE         
           FMin = MIN( FMin,s )
           FMax = MAX( FMax,s )
         END IF
       END DO
       
       
       CALL Info( 'ViewFactors', ' ', Level=3 )
       CALL info( 'ViewFactors', 'Viewfactors before manipulation: ', Level=3 )
       CALL Info( 'ViewFactors', ' ', Level=3 )
       WRITE( LMessage, * ) '        Minimum row sum: ',FMin
       CALL Info( 'ViewFactors', LMessage, Level=3 )
       WRITE( LMessage, * ) '        Maximum row sum: ',FMax
       CALL Info( 'ViewFactors', LMessage, Level=3 )
       CALL Info( 'ViewFactors', ' ', Level=3 )
       
       
       at0 = CPUTime()

       IF( RadiationOpen ) THEN
         CALL Info( 'ViewFactors','Symmetrizing Factors... ', Level=3 )
       ELSE
         CALL Info( 'ViewFactors','Normalizaing Factors...',Level=3)
       END IF

       CALL NormalizeFactors( Model )

       DO i=1,N
         s = 0.0D0
         DO j=1,N
           s = s + Factors((i-1)*N+j)
         END DO
         IF(i == 1) THEN
           Fmin = s
           Fmax = s
         ELSE
           FMin = MIN( FMin,s )
           FMax = MAX( FMax,s )
         END IF
       END DO
       
       WRITE (Message,'(A,F8.2)') 'View factors manipulated in time (s):',CPUTime()-at0
       CALL Info( 'ViewFactors',Message, Level=3 )
       
       CALL Info( 'ViewFactors', ' ', Level=3 )
       CALL info( 'ViewFactors', 'Viewfactors after manipulation: ', Level=3 )
       CALL Info( 'ViewFactors', ' ', Level=3 )
       WRITE( LMessage, * ) '        Minimum row sum: ',FMin
       CALL Info( 'ViewFactors', LMessage, Level=3 )
       WRITE( LMessage, * ) '        Maximum row sum: ',FMax
       CALL Info( 'ViewFactors', LMessage, Level=3 )
       IF( FMax > 1.0_dp ) THEN
         CALL Warn('ViewFactors','Rowsum of view factors should not be larger than one!')
       END IF
       CALL Info( 'ViewFactors', ' ', Level=3 )


       ViewFactorsFile = GetString( GetSimulation(),'View Factors',GotIt)
       IF ( .NOT.GotIt ) ViewFactorsFile = 'ViewFactors.dat'
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
       
       OPEN( 1,File=TRIM(OutputName),STATUS='UNKNOWN' )
       
       ! Use loser constraint for MinFactor as the errors can't be renormalized any more 
       MinFactor = MinFactor / 10.0
       
       DO i=1,N
         k = 0
         DO j=1,N
           IF ( Factors((i-1)*N+j) > MinFactor ) k = k + 1
         END DO
         WRITE( 1,* ) k
         DO j=1,N
           IF ( Factors((i-1)*N+j) > MinFactor ) THEN
             WRITE( 1,* ) i,j,Factors((i-1)*N+j)
           END IF
         END DO
       END DO
       
       CLOSE(1)

       IF ( CylindricSymmetry ) THEN
         DEALLOCATE( Surfaces, Factors)
       ELSE
         DEALLOCATE( Normals, Factors, Surfaces, TYPE)
       END IF
       
     END DO  ! Of radiation RadiationBody

     CALL Info( 'ViewFactors', '*** ALL DONE ***' )
     CALL FLUSH(6)

   CONTAINS
   
!> View factors are normalized in order to improve the numberical accuracy. With 
!> normalization it is ensured that all boundary elements see exactly half 
!> space. 
!------------------------------------------------------------------------------

      SUBROUTINE NormalizeFactors( Model )
        IMPLICIT NONE
        TYPE(Model_t), POINTER :: Model

        INTEGER :: itmax,it,i,j,k

        LOGICAL :: li,lj

        REAL(KIND=dp), ALLOCATABLE :: RHS(:),SOL(:),Areas(:),PSOL(:)

        REAL(KIND=dp) :: cum,s,si,sj
        REAL(KIND=dp), PARAMETER :: eps=1.0D-20
#ifdef USE_ISO_C_BINDINGS
        REAL(KIND=dp) :: at1
#else
        REAL(KIND=dp) :: at1, CPUTime
#endif

        itmax = 20
        it = 0
        cum = 0.0_dp
        
        ALLOCATE( Areas(n),STAT=istat )
        IF ( istat /= 0 ) THEN
          CALL Fatal( 'Viewfactors', &
              'Memory allocation error 1 in NormalizeFactors.Aborting.' )
        END IF

!------------------------------------------------------------------------------
!       First force the matrix (before dividing by area) to be symmetric
!------------------------------------------------------------------------------
        DO i=1,n
          Areas(i) = ElementArea( Model % Mesh, RadElements(i), &
               RadElements(i) % TYPE % NumberOfNodes )
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
            CALL Fatal( 'Viewfactors', &
                'Memory allocation error 2 in NormalizeFactors.Aborting.' )
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
            
            WRITE (Message,'(A,I2,A,ES12.3)') 'Normalization iteration',it,': ',cum
            CALL Info( 'ViewFactors',Message, Level=3 )
            
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

      REAL(KIND=dp), DIMENSION(n) :: x,b

      REAL(KIND=dp) :: dpar(50)

      INTEGER :: N,ipar(50),wsize
      REAL(KIND=dp), ALLOCATABLE :: work(:,:)

      INTEGER(KIND=addrInt) :: iterProc, mvProc, pcondProc, dProc
#ifndef USE_ISO_C_BINDINGS
      INTEGER  :: HUTI_D_CG
      EXTERNAL :: HUTI_D_CG
      INTEGER(KIND=AddrInt) :: AddrFunc
#else
      INTEGER(KIND=AddrInt) :: AddrFunc
      EXTERNAL :: AddrFunc
#endif
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
 
      ALLOCATE( work(N, wsize) )

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

  END PROGRAM ViewFactors


  SUBROUTINE DiagPrecViewFact( u,v,ipar )
    USE ViewFactorGlobals

    REAL(KIND=dp) :: u(*),v(*)
    INTEGER :: ipar(*)

    INTEGER :: n

    n = HUTI_NDIM
    u(1:n) = v(1:n)*Jdiag(1:n)
  END SUBROUTINE DiagPrecViewFact


  SUBROUTINE MatvecViewFact( u,v,ipar )
    USE ViewFactorGlobals

    INTEGER :: ipar(*)
    REAL(KIND=dp) :: u(*),v(*)

    INTEGER :: n

    n = HUTI_NDIM
    CALL DGEMV('N',n,n,1.0_dp,Jacobian,n,u,1,0.0_dp,v,1)
  END SUBROUTINE MatvecViewFact

  
!> \}
!> \}  
