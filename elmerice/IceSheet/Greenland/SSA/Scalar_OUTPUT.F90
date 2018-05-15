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
! *  Author: F. Gillet-Chaulet (IGE)
! *  Email:  fabien.gillet-chaulet@univ-grenoble-alpes.fr
! *  Web:    http://elmerice.elmerfem.org
! *
! *  Original Date: 03-2018, 
! *****************************************************************************
!!! Compute standard 1D variables:
!   1: Time
!   2: Volume
!   3: Volume Above Floatation
!   4: Volume rate of change
!   5: SMB Flux
!   6: BMB Flux
!   7: Residual Flux
!   8: Ice Discharge
!   9: Ice flux at Grounding Line
!  10: Grounded ice area
!  11: Floating ice area
!  12: Ice Free area
! *****************************************************************************     
      SUBROUTINE Scalar_OUTPUT( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      IMPLICIT NONE

      TYPE(Model_t) :: Model
      TYPE(Solver_t):: Solver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

      TYPE(Variable_t),POINTER :: GMVAR,FlowVar,HVar,HRVar,BedVar,DHDTVar
      INTEGER,POINTER :: Permutation(:)

      REAL(KIND=dp) :: Volume,VAF
      REAL(KIND=dp) :: DHDTFlux,SMBFlux,BMBFlux,HMinFlux
      REAL(KIND=dp) :: GroundedArea,FloatingArea,FreeArea
      REAL(KIND=dp) :: CalvingFlux
      REAL(KIND=dp) :: GLFlux
      REAL(KIND=dp),SAVE :: zsea,rhow

      REAL (KIND=dp), ALLOCATABLE, DIMENSION(:),SAVE :: NodeArea
      REAL (KIND=dp), ALLOCATABLE, DIMENSION(:),SAVE :: LocalArea
      REAL (KIND=dp), ALLOCATABLE, DIMENSION(:),SAVE :: NodalH,NodalDHDT,MinH,NodalGM
      REAL (KIND=dp), ALLOCATABLE, DIMENSION(:),SAVE :: NodalSMB,NodalBMB,NodalMB
      REAL (KIND=dp), ALLOCATABLE, DIMENSION(:),SAVE :: NodalHf
      REAL (KIND=dp), ALLOCATABLE, DIMENSION(:),SAVE :: rhoi
      REAL (KIND=dp), ALLOCATABLE, DIMENSION(:),SAVE :: Val,ParVal
      REAL(KIND=dp),ALLOCATABLE,SAVE :: Basis(:), dBasisdx(:,:)

      INTEGER :: i
      INTEGER :: ierr
      INTEGER,PARAMETER :: NVal=11
      INTEGER, PARAMETER :: io=12
      INTEGER,PARAMETER :: DIM=2 !dimension of the pb restricted to 2 currently

      INTEGER :: FlowDofs

      LOGICAL,SAVE :: Firsttime=.TRUE.

      CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='INITMIP_Scalar_OUTPUT'
      CHARACTER(LEN=MAX_NAME_LEN),SAVE :: OUTPUT_FName
      CHARACTER(LEN=MAX_NAME_LEN),ALLOCATABLE,SAVE :: ValueNames(:)

      CALL GET_VARIABLES()

      IF (Firsttime.OR.Solver%Mesh%Changed) THEN

        IF (.NOT.ASSOCIATED(Solver%Variable)) & 
         CALL FATAL(SolverName,'Solver%Variable Not associated')
        IF (.NOT.ASSOCIATED(Solver%Matrix)) &
         CALL FATAL(SolverName,'Solver%Matrix Not associated')

        IF ( CurrentCoordinateSystem() /= Cartesian )  &
          CALL FATAL(SolverName,'Only For cartesian system')

        IF ( Model % Mesh % MeshDim /= DIM ) &
          CALL FATAL(SolverName,'Only For 2D plan view')

       !## DO SOME ALLOCATION
        CALL DO_ALLOCATION(Firsttime)

       !## Name of Saved variables
        ValueNames(1)='Volume'
        ValueNames(2)='Volume Above Floatation'
        ValueNames(3)='Volume rate of change'
        ValueNames(4)='SMB Flux'
        ValueNames(5)='BMB Flux'
        ValueNames(6)='Residual Flux'
        ValueNames(7)='Ice Discharge'
        ValueNames(8)='Ice flux at Grounding Line'
        ValueNames(9)='Grounded ice area'
        ValueNames(10)='Floating ice area'
        ValueNames(11)='Ice Free area'

        IF (Firsttime) CALL GET_CONSTANTS(zsea,rhow)
        IF (Firsttime) CALL INIT_OUTPUT_FILE(OUTPUT_FName)
        IF (Firsttime) CALL COMPUTE_NodeArea(NodeArea)
        Firsttime=.FALSE.          
      END IF


      CALL BODY_INTEGRATION(Volume,VAF,DHDTFlux,SMBFlux,BMBFlux,HMinFlux,&
               GroundedArea,FloatingArea,FreeArea)

      CALL BC_INTEGRATION(CalvingFlux)

      CALL COMPUTE_GL_FLUX(GLFlux)


      Val(1)=Volume
      Val(2)=VAF

      Val(3)=DHDTFlux
      Val(4)=SMBFlux
      Val(5)=BMBFlux
      Val(6)=HMinFlux

      Val(7)=CalvingFlux
      Val(8)= GLFlux

      Val(9)=GroundedArea
      Val(10)=FloatingArea
      Val(11)=FreeArea

      IF (ParEnv % PEs > 1 ) THEN
        DO i=1,NVal
         CALL MPI_ALLREDUCE(Val(i),ParVal(i),1,&
                       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        END DO
        Val(1:NVal)=ParVal(1:NVal)
      END IF

      IF ((ParEnv % PEs > 1 ).AND.(ParEnv % MyPe.NE.0)) RETURN

      IF( Solver % TimesVisited > 0 ) THEN
        OPEN(io,file=TRIM(OUTPUT_FName),position='append')
      ELSE
        OPEN(io,file=TRIM(OUTPUT_FName))
      END IF

      write(io,'(ES22.12E3)',advance='no') GetTime()
      Do i=1,NVal-1
        write(io,'(ES22.12E3)',advance='no') Val(i)
      End do
        write(io,'(ES22.12E3)') Val(NVal)
      CLOSE(io)

      CONTAINS

      SUBROUTINE DO_ALLOCATION(Firsttime)
      LOGICAL,INTENT(IN) :: Firsttime
      INTEGER :: M
      INTEGER :: N
         IF (.NOT.Firsttime) &
            DEALLOCATE(NodalH,NodalDHDT,NodalHf,MinH,NodalGM,NodalSMB,NodalBMB,NodalMB,rhoi,&
                      LocalArea, &
                      Basis,dBasisdx,&
                      Val,ParVal,ValueNames,&
                      NodeArea)
          N=Model % Mesh % NumberOfNodes
          M=Model % MaxElementNodes
          ALLOCATE(Basis(M),&
                   dBasisdx(M,3),&
                   NodalH(M),&
                   NodalDHDT(M),&
                   NodalHf(M),&
                   MinH(M),&
                   NodalGM(M),&
                   NodalSMB(M),&
                   NodalBMB(M),&
                   NodalMB(M),&
                   rhoi(M),&
                   LocalArea(M),&
                   Val(NVal),ParVal(NVal),ValueNames(NVal),&
                   NodeArea(N))
      END SUBROUTINE DO_ALLOCATION

      SUBROUTINE GET_CONSTANTS(zsea,rhow)
      IMPLICIT NONE
      REAL(KIND=dp),INTENT(OUT) :: zsea,rhow
      LOGICAL :: Found

        zsea = GetCReal( Model % Constants, 'Sea Level', Found )
        IF (.NOT.Found) CALL FATAL(SolverName,'<Sea Level> not found')
        rhow = GetCReal( Model % Constants, 'water density', Found )
        IF (.NOT.Found) CALL FATAL(SolverName,'<water density not found')
      END SUBROUTINE GET_CONSTANTS

      SUBROUTINE INIT_OUTPUT_FILE(OUTPUT_FName)
      USE GeneralUtils
      IMPLICIT NONE
      CHARACTER(LEN=MAX_NAME_LEN),INTENT(OUT) :: OUTPUT_FName

      CHARACTER(LEN=MAX_NAME_LEN) ::NamesFile,&
                   OUTPUT_FName_D='INITMIP_Scalar_OUTPUT.dat'

      CHARACTER(LEN=MAX_NAME_LEN) :: DateStr 
      TYPE(ValueList_t), POINTER :: SolverParams
      LOGICAL :: Found
      INTEGER :: i

       SolverParams=>GetSolverParams(Solver)
       OUTPUT_FName = ListGetString(SolverParams,'File Name',Found)
       IF (.NOT.Found) OUTPUT_FName=OUTPUT_FName_D

       NamesFile = TRIM(OUTPUT_FName) // '.' // TRIM("names")
         
       IF ((ParEnv % PEs >1).AND.(ParEnv%MyPe.NE.0)) RETURN

       DateStr = FormatDate()

       OPEN(io,file=TRIM(NamesFile))
       WRITE(io,'(A)') 'File started at: '//TRIM(DateStr)
       WRITE(io,'(A)') ' '
       WRITE(io,'(A)') 'Elmer version: '//TRIM(GetVersion())
       WRITE(io,'(A)') 'Elmer revision: '//TRIM(GetRevision())
       WRITE(io,'(A)') 'Elmer Compilation Date: '//TRIM(GetCompilationDate())
       WRITE(io,'(A)') ' '
       WRITE(io,'(A)') 'Variables in columns of matrix:'//TRIM(OUTPUT_FName)
       WRITE(io,'(I4,": ",A)') 1,'Time'
       DO i=1,NVal
          WRITE(io,'(I4,": ",A)') i+1,TRIM(ValueNames(i))
       END DO
       CLOSE(io)
      END SUBROUTINE INIT_OUTPUT_FILE

      SUBROUTINE GET_VARIABLES()
       HVar    => VariableGet(Solver%Mesh%Variables,'h',UnfoundFatal=.TRUE.)

       DHDTVar => VariableGet(Solver%Mesh%Variables,'dhdt',UnfoundFatal=.TRUE.)

       HRVar   => VariableGet(Solver%Mesh%Variables,'h loads')
       IF (.NOT.ASSOCIATED(HRVar)) &
         HRVar   => VariableGet(Solver%Mesh%Variables,'h residual',UnfoundFatal=.TRUE.)

       BedVar  => VariableGet(Solver%Mesh%Variables,'bedrock',UnfoundFatal=.TRUE.)

       GMVar   => VariableGet(Solver%Mesh%Variables,'GroundedMask',UnfoundFatal=.TRUE.)

       FlowVar => VariableGet(Solver%Mesh%Variables,'SSAVelocity',UnfoundFatal=.TRUE.)
       FlowDofs = FlowVar % DOFs

       Permutation => Solver%Variable%Perm
      END SUBROUTINE GET_VARIABLES

      SUBROUTINE COMPUTE_NodeArea(NodeArea)
      IMPLICIT NONE
      REAL(KIND=dp),INTENT(OUT) :: NodeArea(:)

      TYPE(Element_t), POINTER :: Element
      TYPE(Nodes_t),SAVE :: ElementNodes
      TYPE(GaussIntegrationPoints_t) :: IntegStuff
      REAL(KIND=dp) :: U,V,W,SqrtElementMetric
      INTEGER,POINTER :: Indexes(:)
      INTEGER :: n
      INTEGER :: t,i
      LOGICAL :: stat


      NodeArea=0._dp

      Do t=1,GetNOFActive()

         Element => GetActiveElement(t)
         n = GetElementNOFNodes(Element)
         Indexes => Element % NodeIndexes
         CALL GetElementNodes( ElementNodes, Element )
         IntegStuff = GaussPoints( Element )

         Do i=1,IntegStuff % n
            U = IntegStuff % u(i)
            V = IntegStuff % v(i)
            W = IntegStuff % w(i)
            stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric, &
                        Basis,dBasisdx )

            NodeArea(Permutation(Indexes(1:n)))=NodeArea(Permutation(Indexes(1:n)))+&
                SqrtElementMetric*IntegStuff % s(i) * Basis(1:n)

         End do
      End do
      IF (ParEnv % PEs > 1 ) CALL ParallelSumVector( Solver % Matrix, NodeArea, 0 )
  
      END SUBROUTINE COMPUTE_NodeArea

      SUBROUTINE BODY_INTEGRATION(Volume,VAF,DHDTFlux,SMBFlux,BMBFlux,HMinFlux,&
                     GroundedArea,FloatingArea,FreeArea)
      IMPLICIT NONE
      REAL(KIND=dp),INTENT(OUT) :: Volume,VAF,&
                       DHDTFlux,SMBFlux,BMBFlux,HMinFlux,&
                       GroundedArea,FloatingArea,FreeArea
      REAL(KIND=dp),parameter :: Fsmall=100.0*EPSILON(1.0)

      TYPE(Element_t),POINTER :: Element
      TYPE(ValueList_t), POINTER :: BodyForce,Material
      TYPE(GaussIntegrationPoints_t) :: IntegStuff
      TYPE(Nodes_t),SAVE :: ElementNodes

      REAL(KIND=dp) :: U,V,W,SqrtElementMetric
      REAL(KIND=dp) :: Normal(3),Flow(3)
      REAL(KIND=dp) :: cellarea
      REAL(KIND=dp) :: HAtIP,SMBAtIP,BMBAtIP

      LOGICAL :: CalvingFront
      LOGICAL :: IsFloating,IceFree
      LOGICAL :: stat

      INTEGER,POINTER :: NodeIndexes(:)
      INTEGER :: t
      INTEGER :: i
      INTEGER :: n
      INTEGER :: ne

      ne=GetNOFActive()

      Volume=0._dp
      VAF=0._dp

      DHDTFlux=0._dp
      SMBFlux=0._dp
      BMBFlux=0._dp
      HMinFlux=0._dp

      GroundedArea=0._dp
      FloatingArea=0._dp
      FreeArea=0._dp

      DO t = 1,ne

         Element => GetActiveElement(t)

         IF ( CheckPassiveElement(Element) )  CYCLE

         n = GetElementNOFNodes(Element)
         NodeIndexes => Element % NodeIndexes
         CALL GetElementNodes( ElementNodes )

         BodyForce => GetBodyForce(Element)
         IF (.NOT.ASSOCIATED(BodyForce)) &
            CALL FATAL(SolverName,'No BodyForce Found')
         Material => GetMaterial(Element)
         IF (.NOT.ASSOCIATED(Material)) &
            CALL FATAL(SolverName,'No Material Found')

         NodalH(1:n) = HVar%Values(HVar%Perm(NodeIndexes(1:n)))
         NodalDHDT(1:n) = DHDTVar%Values(DHDTVar%Perm(NodeIndexes(1:n)))
         
         rhoi(1:n) = ListGetReal(Material,'SSA Mean Density',n,NodeIndexes,UnfoundFatal=.TRUE. )

         Do i=1,n
           NodalHf(i)=Max(0._dp,NodalH(i)-&
            Max(0._dp,(zsea-BedVar%Values(BedVar%Perm(NodeIndexes(i))))*rhow/rhoi(i)))
         End do

         MinH=0._dp
         MinH(1:n) = ListGetReal(Material,'Min H',n,NodeIndexes,UnfoundFatal=.TRUE. )

        
        ! FLOATING OR GROUNDED CELL
         NodalGM(1:n) = GMVar%Values(GMVar%Perm(NodeIndexes(1:n)))
         IsFloating=ANY(NodalGM(1:n).LT.0._dp)

        ! TOP ACCUMULATION
         NodalSMB=0._dp
         NodalSMB(1:n) = &
            ListGetReal(BodyForce,'Top Surface Accumulation', n,NodeIndexes,UnfoundFatal=.TRUE. )

        ! Bottom ACCUMULATION
         NodalBMB=0._dp
         NodalBMB(1:n) = &
            ListGetReal(BodyForce,'Bottom Surface Accumulation', n,NodeIndexes,UnfoundFatal=.TRUE. )

        ! Total ACCUMULATION
         NodalMB(1:n) = NodalSMB(1:n) + NodalBMB(1:n)

        ! CELL IS NOT ACTIVE ALL H VALUES BELOW MinH Value
         IceFree=.FALSE.
         IF (ALL((NodalH(1:n)-MinH(1:n)).LE.Fsmall).AND.ALL(NodalMB(1:n).LT.0._dp)) IceFree=.TRUE.

        ! GO TO INTEGRATION
         cellarea=0._dp
         LocalArea=0._dp

         IntegStuff = GaussPoints( Element )
         DO i=1,IntegStuff % n
            U = IntegStuff % u(i)
            V = IntegStuff % v(i)
            W = IntegStuff % w(i)

            stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric, &
                        Basis,dBasisdx )
           ! cell area
           cellarea=cellarea+SqrtElementMetric*IntegStuff % s(i)
           ! the area seen by each node
           LocalArea(1:n)=LocalArea(1:n)+SqrtElementMetric*IntegStuff % s(i) * Basis(1:n)

           IF (IceFree) CYCLE

           ! Integrate H
           HAtIP=SUM(NodalH(1:n)*Basis(1:n))
           Volume=Volume+HAtIP*SqrtElementMetric*IntegStuff % s(i)

           VAF=VAF+SUM(NodalHf(1:n)*Basis(1:n))*SqrtElementMetric*IntegStuff % s(i)

           SMBAtIP=SUM(NodalSMB(1:n)*Basis(1:n))
           BMBAtIP=SUM(NodalBMB(1:n)*Basis(1:n))

           DHDTFlux=DHDTFlux+SUM(NodalDHDT(1:n)*Basis(1:n))*SqrtElementMetric*IntegStuff % s(i)
           SMBFlux=SMBFlux+SMBAtIP*SqrtElementMetric*IntegStuff % s(i)
           BMBFlux=BMBFlux+BMBAtIP*SqrtElementMetric*IntegStuff % s(i)

         End DO

         IF (.NOT.IceFree) THEN
           Do i=1,n
             HMinFlux = HMinFlux + &
              HRVar%Values(HRVar%Perm(NodeIndexes(i)))*LocalArea(i)/NodeArea(Permutation(NodeIndexes(i)))
           End do
           IF (IsFloating) THEN
              FloatingArea=FloatingArea+cellarea
           ELSE
              GroundedArea=GroundedArea+cellarea
           END IF
         ELSE
           FreeArea=FreeArea+cellarea
         END IF
      End do
      END SUBROUTINE BODY_INTEGRATION

      SUBROUTINE BC_INTEGRATION(CalvingFlux)
      IMPLICIT NONE
      REAL(KIND=dp),INTENT(OUT) :: CalvingFlux

      TYPE(Element_t),POINTER :: Element
      TYPE(ValueList_t), POINTER :: BC
      TYPE(GaussIntegrationPoints_t) :: IntegStuff
      TYPE(Nodes_t),SAVE :: ElementNodes

      REAL(KIND=dp) :: U,V,W,SqrtElementMetric
      REAL(KIND=dp) :: Normal(3),Flow(3)

      LOGICAL :: CalvingFront
      LOGICAL :: Found
      LOGICAL :: stat

      INTEGER,POINTER :: NodeIndexes(:)
      INTEGER :: t
      INTEGER :: i,j
      INTEGER :: n

      CalvingFlux=0._dp
      DO t = 1,GetNOFBoundaryElements()
         Element => GetBoundaryElement(t)

         IF ( .NOT. ActiveBoundaryElement() ) CYCLE
         IF ( GetElementFamily() == 1 ) CYCLE

         BC => GetBC()
         IF ( .NOT. ASSOCIATED(BC) ) CYCLE
         CalvingFront=.FALSE.
         CalvingFront=ListGetLogical(BC,'Calving Front', Found)
         IF (.NOT.CalvingFront) CYCLE

         n = GetElementNOFNodes()
         NodeIndexes => Element % NodeIndexes
         CALL GetElementNodes( ElementNodes )

         NodalH(1:n) = HVar%Values(HVar%Perm(NodeIndexes(1:n)))

         IntegStuff = GaussPoints( Element )
         DO i=1,IntegStuff % n
           U = IntegStuff % u(i)
           V = IntegStuff % v(i)
           W = IntegStuff % w(i)
           stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric, &
                                   Basis,dBasisdx )
           Normal=0._dp
           Normal = NormalVector( Element,ElementNodes,u,v,.TRUE. )

           Flow=0._dp
           DO j=1,FlowDofs
              Flow(j) = SUM( FlowVar % Values(FlowDofs*(FlowVar % Perm(NodeIndexes(1:n)) -1)+j) * Basis(1:n) )
           END DO
           CalvingFlux=CalvingFlux+&
                       SUM(NodalH(1:n)*Basis(1:n))*SUM(Normal * Flow)*SqrtElementMetric*IntegStuff % s(i)
         END DO
      END DO
      END SUBROUTINE BC_INTEGRATION

    
      SUBROUTINE COMPUTE_GL_FLUX( GLFlux  )
      IMPLICIT NONE
      REAL(KIND=DP),INTENT(OUT) :: GLFlux
      
      TYPE(Mesh_t),POINTER :: Mesh
      TYPE(Element_t),DIMENSION(:),POINTER :: Edges
      TYPE(Element_t),POINTER :: Edge,Parent
    
      REAL(KIND=DP),ALLOCATABLE,SAVE :: LocalGM(:),LocalFlow(:,:),LocalH(:)
    
      INTEGER, POINTER :: NodeIndexes(:)
      INTEGER ::  Ne 
      INTEGER :: n
      INTEGER :: i,j
      INTEGER :: M
    
      LOGICAL, SAVE :: Firsttime=.TRUE.
    
    
      IF (Firsttime) THEN
        M=Model % MaxElementNodes
        ALLOCATE(LocalGM(M),LocalFlow(3,M),LocalH(M))
        Firsttime=.FALSE.
      END IF

      Mesh => GetMesh()


      CALL FindMeshEdges(Mesh,.FALSE.) 
  
      Edges => Mesh % Edges
      Ne = Mesh % NumberOfEdges
 
      GLFlux=0._dp
    
      DO i=1,Ne
        Edge => Edges(i)
        n=Edge % TYPE % NumberOfNodes
        NodeIndexes(1:n) => Edge % NodeIndexes(1:n)
      
        LocalGM(1:n)=GMVar % Values ( GMVar % Perm ( NodeIndexes(1:n) ) )
        ! Edge is GL if all GM=0
        IF ( ANY( abs(LocalGM(1:n)) .GT. AEPS ) ) CYCLE
      
        DO j=1,DIM
           LocalFlow(j,1:n) = FlowVar % Values ( DIM*(FlowVar % Perm ( NodeIndexes(1:n) ) - 1) + j )
        END DO
        LocalH(1:n) = HVar % Values ( HVar % Perm ( NodeIndexes(1:n) ) )
      
        Parent => Edge % BoundaryInfo % Right
        CALL AddLocalFlux(GLFlux,LocalFlow,LocalH,Edge,Parent,Mesh)
        Parent => Edge % BoundaryInfo % Left 
        CALL AddLocalFlux(GLFlux,LocalFlow,LocalH,Edge,Parent,Mesh)
      END DO
    
      END SUBROUTINE COMPUTE_GL_FLUX
     
      SUBROUTINE AddLocalFlux(Flux,LocalFlow,LocalH,Edge,Parent,Mesh) 
      IMPLICIT NONE
      TYPE(Element_t),POINTER :: Edge,Parent
      TYPE(Mesh_t), POINTER :: Mesh
      REAL(KIND=dp),INTENT(INOUT) :: Flux
      REAL(KIND=dp),INTENT(IN) :: LocalFlow(:,:),LocalH(:)
      
      TYPE(Nodes_t),SAVE :: EdgeNodes
      TYPE(GaussIntegrationPoints_t) :: IntegStuff
      REAL(KIND=dp) :: U,V,W,SqrtElementMetric
      REAL(KIND=dp),ALLOCATABLE,SAVE :: Basis(:), dBasisdx(:,:)
      
      REAL(KIND=dp),DIMENSION(3) :: dx, Normal,Flow
      REAL(KIND=dp) :: h
      
      INTEGER :: i,j
      INTEGER :: n,np 
      INTEGER :: M
      
      LOGICAL :: stat
      LOGICAL,SAVE :: FirstVisit=.TRUE.
      
      IF (FirstVisit) THEN
        M=Model % MaxElementNodes
        ALLOCATE(Basis(M),dBasisdx(M,3))
        FirstVisit=.FALSE.
      END IF
     
       ! Parent not associated
       IF (.NOT.ASSOCIATED(Parent)) RETURN
       ! Parent is halo element
       IF( Parent % PartIndex /= ParEnv % MyPe ) RETURN
       !
       n = Edge % TYPE % NumberOfNodes
       np = Parent % TYPE % NumberOfNodes

       ! GL could be a boundary; compute flux from gounded parent
       ! Parent is Grounded if ALL GM>=0
       IF (.NOT.( ALL( GMVar % Values( GMVar % Perm( Parent % NodeIndexes(1:np) ) ) .GT. -AEPS ) ) ) RETURN
       ! a vector from the center of the edge to the center of the parent to check that normal points outside parent 
       dx = ElementCenter(Mesh,Parent) - ElementCenter(Mesh,Edge)
     
       CALL GetElementNodes(EdgeNodes, Edge)
        
       IntegStuff = GaussPoints( Edge ) 
       DO i=1,IntegStuff % n
          U = IntegStuff % u(i)
          V = IntegStuff % v(i)
          W = IntegStuff % w(i)
          stat = ElementInfo( Edge,EdgeNodes,U,V,W,SqrtElementMetric, &
              Basis,dBasisdx )
              
         DO j=1,DIM
                Flow(j) = SUM(LocalFlow(j,1:n)*Basis(1:n))
         END DO  
         h = SUM(LocalH(1:n)*Basis(1:n))

         Normal = NormalVector( Edge,EdgeNodes,U,V,.FALSE. )
         
         IF ( SUM(dx(1:DIM)*Normal(1:DIM)).GT.0._dp ) Normal=-Normal
         
         Flux = Flux + SqrtElementMetric * IntegStuff % s(i) * h * SUM(Normal(1:DIM)*Flow(1:DIM))
       END DO
     
      END SUBROUTINE AddLocalFlux
!
!      
      FUNCTION ElementCenter(Mesh,Element) RESULT(xc)
      IMPLICIT NONE
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Element_t),POINTER :: Element
      REAL(KIND=dp),DIMENSION(3) :: xc
      INTEGER :: n
        
       n = Element % TYPE % NumberOfNodes
       xc=0._dp
       SELECT CASE( Element % TYPE % ElementCode / 100 )
         CASE(2,4,8)
          xc(1) = InterpolateInElement( Element, Mesh % Nodes % x(Element % NodeIndexes(1:n)), 0.0d0, 0.0d0, 0.0d0 )
          xc(2) = InterpolateInElement( Element, Mesh % Nodes % y(Element % NodeIndexes(1:n)), 0.0d0, 0.0d0, 0.0d0 )
        CASE(3)
          xc(1) = InterpolateInElement( Element, Mesh % Nodes % x(Element % NodeIndexes(1:n)), 1.0d0/3, 1.0d0/3, 0.0d0 )
          xc(2) = InterpolateInElement( Element, Mesh % Nodes % y(Element % NodeIndexes(1:n)), 1.0d0/3, 1.0d0/3, 0.0d0 )
         CASE DEFAULT
            CALL FATAL(SolverName,'Element type not supported')
       END SELECT
     
      END FUNCTION ElementCenter
     
      


      END
