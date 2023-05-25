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
! *  11-2020
! *  O. Gagliardini (IGE) modified the initial version from Fabien to adapt it for Glacier
! *****************************************************************************
!!! Compute standard 1D variables:
!   1: Time
!   2: Volume
!   3: Area    
!   4: Ablation Area    
!   5: Accumulation Area    
!   6: SMB Total
!   7: SMB Ablation 
!   8: SMB Accumulation 
!   9: Front elevation
!  10: X Front
!  11: Y Front
!  12: Front distance
! *****************************************************************************     
      SUBROUTINE Scalar_OUTPUT( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      IMPLICIT NONE

      TYPE(Model_t) :: Model
      TYPE(Solver_t):: Solver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

      TYPE(Variable_t),POINTER :: SMBVar,HVar,BedVar,ZsVar,ZsLoads,DzDtVar,EVar
      INTEGER,POINTER :: Permutation(:)

      REAL(KIND=dp) :: Volume,Vchange
      REAL(KIND=dp) :: TotalArea,AblaArea,AccuArea
      REAL(KIND=dp) :: TotalSMB,AblaSMB,AccuSMB
      REAL(KIND=dp) :: Hfront,XFront,YFront,Dist
      REAL(KIND=dp) :: ZsLoadsFlux

      REAL (KIND=dp), ALLOCATABLE, DIMENSION(:),SAVE :: NodeArea
      REAL (KIND=dp), ALLOCATABLE, DIMENSION(:),SAVE :: LocalArea
      REAL (KIND=dp), ALLOCATABLE, DIMENSION(:),SAVE :: NodalH,NodalBed,MinH
      REAL (KIND=dp), ALLOCATABLE, DIMENSION(:),SAVE :: NodalSMB,NodalZs,NodalDzDt
      REAL (KIND=dp), ALLOCATABLE, DIMENSION(:),SAVE :: Val,ParVal
      REAL(KIND=dp),ALLOCATABLE,SAVE :: Basis(:), dBasisdx(:,:)

      INTEGER :: i
      INTEGER :: ierr
      INTEGER,PARAMETER :: NVal=13
      INTEGER, PARAMETER :: io=12
      INTEGER,PARAMETER :: DIM=3 !dimension of the pb restricted to 3 currently

      INTEGER :: FlowDofs

      LOGICAL,SAVE :: Firsttime=.TRUE.

      CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='Scalar_OUTPUT_Glacier'
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
        ValueNames(2)='Volume Change'

        ValueNames(3)='Area'
        ValueNames(4)='Ablation Area'
        ValueNames(5)='Accumulation Area'

        ValueNames(6)='SMB Total'
        ValueNames(7)='SMB Ablation'
        ValueNames(8)='SMB Accumulation'
        ValueNames(9)='Residual Flux'


        ValueNames(10)='Front Elevation'
        ValueNames(11)='X Front'
        ValueNames(12)='Y Front'
        ValueNames(13)='Front Distance'

        IF (Firsttime) CALL INIT_OUTPUT_FILE(OUTPUT_FName)
        IF (Firsttime) CALL COMPUTE_NodeArea(NodeArea)
        Firsttime=.FALSE.          
      END IF

      CALL BODY_INTEGRATION(Volume,TotalArea,AblaARea,AccuArea,TotalSMB, &
               AblaSMB,AccuSMB,Hfront,XFront,YFront,Dist,ZsLoadsFlux,Vchange) 

      Val(1)=Volume
      Val(2)=Vchange

      Val(3)=TotalArea
      Val(4)=AblaArea
      Val(5)=AccuArea

      Val(6)=TotalSMB
      Val(7)=AblaSMB
      Val(8)=AccuSMB
      Val(9)= ZsLoadsFlux

      Val(10)= Hfront
      Val(11)= Xfront
      Val(12)= Yfront
      Val(13)= Dist

      IF (ParEnv % PEs > 1 ) THEN
         
        DO i=1,9
           CALL MPI_ALLREDUCE(Val(i),ParVal(i),1,&
                       MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
        END DO
        CALL MPI_ALLREDUCE(Val(10),ParVal(10),1,&
                       MPI_DOUBLE_PRECISION,MPI_MIN,ELMER_COMM_WORLD,ierr)

        IF (Val(10).NE.ParVal(10)) THEN
           Val(11) = Huge(1.0)
           Val(12) = Huge(1.0)
        END IF
        DO i=11,12
           CALL MPI_ALLREDUCE(Val(i),ParVal(i),1,&
                       MPI_DOUBLE_PRECISION,MPI_MIN,ELMER_COMM_WORLD,ierr)
        END DO

        CALL MPI_ALLREDUCE(Val(13),ParVal(13),1,&
                               MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)

        Val(1:NVal)=ParVal(1:NVal)
      END IF

      IF ((ParEnv % PEs > 1 ).AND.(ParEnv % MyPe.NE.0)) RETURN

      IF( Solver % TimesVisited > 0 ) THEN
        OPEN(io,file=TRIM(OUTPUT_FName),position='append')
      ELSE
        OPEN(io,file=TRIM(OUTPUT_FName))
      END IF

      WRITE(io,'(ES22.12E3)',advance='no') GetTime()
      Do i=1,NVal-1
        WRITE(io,'(ES22.12E3)',advance='no') Val(i)
      End do
        WRITE(io,'(ES22.12E3)') Val(NVal)
      CLOSE(io)

      CONTAINS

      SUBROUTINE DO_ALLOCATION(Firsttime)
      LOGICAL,INTENT(IN) :: Firsttime
      INTEGER :: M
      INTEGER :: N
         IF (.NOT.Firsttime) &
            DEALLOCATE(NodalH,NodalBed,NodalZs,MinH,NodalSMB,NodalDzDt,&
                      LocalArea, &
                      Basis,dBasisdx,&
                      Val,ParVal,ValueNames,&
                      NodeArea)
          N=Model % Mesh % NumberOfNodes
          M=Model % MaxElementNodes
          ALLOCATE(Basis(M),&
                   MinH(M),&
                   dBasisdx(M,3),&
                   NodalH(M),NodalZs(M),NodalBed(M),&
                   NodalSMB(M),NodalDzDt(M),&
                   LocalArea(M),&
                   Val(NVal),ParVal(NVal),ValueNames(NVal),&
                   NodeArea(N))
      END SUBROUTINE DO_ALLOCATION

      SUBROUTINE INIT_OUTPUT_FILE(OUTPUT_FName)
        USE GeneralUtils
        IMPLICIT NONE
        CHARACTER(LEN=MAX_NAME_LEN),INTENT(OUT) :: OUTPUT_FName

        CHARACTER(LEN=MAX_NAME_LEN) ::NamesFile,&
                   OUTPUT_FName_D='glacier_Scalar_OUTPUT.dat'

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
       LOGICAL :: Found
       CHARACTER(LEN=MAX_NAME_LEN) :: PVarName
       TYPE(ValueList_t), POINTER :: SolverParams

       SolverParams=>GetSolverParams(Solver)

       HVar    => VariableGet(Solver%Mesh%Variables,'Thickness',UnfoundFatal=.TRUE.)

       BedVar => VariableGet(Solver%Mesh%Variables,'BedDEM',UnfoundFatal=.TRUE.)

       ZsVar => VariableGet(Solver%Mesh%Variables,'Zs',UnfoundFatal=.TRUE.)

       ZsLoads => VariableGet(Solver%Mesh%Variables,'Zs Loads',UnfoundFatal=.TRUE.)

       DzDtVar => VariableGet(Solver%Mesh%Variables,'Zs Velocity',UnfoundFatal=.TRUE.)

       Permutation => Solver%Variable%Perm

       PVarName = ListGetString(SolverParams,'Passive Variable Name',Found)
       IF (Found) THEN
        EVar => VariableGet(Solver%Mesh%Variables,TRIM(PVarName),UnFoundFatal=.TRUE.)
        IF (EVar % Type .NE. Variable_on_elements) &
          CALL FATAL(TRIM(SolverName),'Passive Variable should be on elements')
        EVar % Values = +1._dp
       ELSE
         EVar => NULL()
       END IF

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
         ElementNodes % z =0._dp


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

      SUBROUTINE BODY_INTEGRATION(Volume,TotalArea,AblaARea,AccuArea,TotalSMB, &
               AblaSMB,AccuSMB,Hfront,Xfront,Yfront,Dist,ZsLoadsFlux,Vchange)
      IMPLICIT NONE
      REAL(KIND=dp),INTENT(OUT) :: Volume,TotalArea,&
                       AblaArea,AccuArea,TotalSMB,AblaSMB,AccuSMB,&
                       Hfront,Xfront,Yfront,Dist,ZsLoadsFlux,Vchange

      REAL(KIND=dp),parameter :: Fsmall=100.0*EPSILON(1.0)

      TYPE(Mesh_t),POINTER :: Mesh
      TYPE(Element_t),POINTER :: Element
      TYPE(ValueList_t), POINTER :: Material,BodyForce
      TYPE(GaussIntegrationPoints_t) :: IntegStuff
      TYPE(Nodes_t),SAVE :: ElementNodes
      TYPE(ValueList_t), POINTER :: SolverParams

      REAL(KIND=dp) :: U,V,W,SqrtElementMetric
      REAL(KIND=dp) :: cellarea,XRefFront,YRefFront
      REAL(KIND=dp) :: HAtIP,SMBAtIP
      REAL(KIND=dp) :: s

      LOGICAL :: IceFree,FoundXRef,FoundYRef
      LOGICAL :: stat
      LOGICAL :: Found

      INTEGER,POINTER :: NodeIndexes(:)
      INTEGER :: t
      INTEGER :: i
      INTEGER :: n
      INTEGER :: ne, nf

      ne=GetNOFActive()
      SolverParams=>GetSolverParams(Solver)
      XRefFront = GetConstReal(SolverParams,'XRefFrontPosition',FoundXRef)
      YRefFront = GetConstReal(SolverParams,'YRefFrontPosition',FoundYRef)
 
      Volume=0._dp
      
      TotalArea = 0._dp
      AblaArea = 0._dp
      AccuArea = 0._dp
      
      TotalSMB = 0._dp
      AblaSMB = 0._dp
      AccuSMB = 0._dp
      
      Hfront = Huge(1.0)
      Dist = Huge(1.0)

      ZsLoadsFlux = 0._dp
      Vchange = 0._dp

      DO t = 1,ne

         Element => GetActiveElement(t)


         n = GetElementNOFNodes(Element)
         NodeIndexes => Element % NodeIndexes
         CALL GetElementNodes( ElementNodes )
         ElementNodes % z = 0._dp

         IF ( CheckPassiveElement(Element) ) THEN
                 IF (ASSOCIATED(EVar)) &
                   EVar % Values (EVar % Perm (element % elementIndex)) = -1.0_dp
                 CYCLE
         ENDIF

         Material => GetMaterial(Element,Found)
         IF (.NOT.Found) &
            CALL FATAL(SolverName,'No Material Found')
         BodyForce => GetBodyForce( Element, Found )
         IF (.NOT.Found) &
            CALL FATAL(SolverName,'No Body Force Found')

         NodalH(1:n) = HVar%Values(HVar%Perm(NodeIndexes(1:n)))
         NodalBed(1:n) = BedVar%Values(BedVar%Perm(NodeIndexes(1:n)))
         NodalZs(1:n) = ZsVar%Values(ZsVar%Perm(NodeIndexes(1:n)))

         NodalDzDt(1:n) = DzDtVar%Values(DzDtVar%Perm(NodeIndexes(1:n)))
         
         !NodalSMB(1:n) = SMBVar%Values(SMBVar%Perm(NodeIndexes(1:n)))
         NodalSMB(1:n) = ListGetReal(BodyForce, 'Zs Accumulation Flux 3',n,NodeIndexes,UnfoundFatal=.TRUE.)

         MinH=0._dp
         MinH(1:n) = ListGetReal(Material,'Min Zs',n,NodeIndexes,UnfoundFatal=.TRUE.)
         MinH(1:n) = MinH(1:n) - NodalBed(1:n)
         
        
        ! CELL IS NOT ACTIVE ALL H VALUES BELOW MinH Value
         IceFree=.FALSE.
         ! IF (ALL((NodalH(1:n)-MinH(1:n)).LE.Fsmall).AND.ALL(NodalSMB(1:n).LT.0._dp)) IceFree=.TRUE.
         IF (ALL((NodalH(1:n)-MinH(1:n)).LE.Fsmall)) IceFree=.TRUE.

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
            s=SqrtElementMetric*IntegStuff % s(i)
           ! cell area
           cellarea=cellarea+s
           ! the area seen by each node
           LocalArea(1:n)=LocalArea(1:n)+ s * Basis(1:n)

           IF (IceFree) CYCLE

           ! Integrate H
           HAtIP=SUM(NodalH(1:n)*Basis(1:n))
           Volume=Volume+HAtIP*s
           Vchange=Vchange+SUM(NodalDzDt(1:n)*Basis(1:n))*s

           SMBAtIP=SUM(NodalSMB(1:n)*Basis(1:n))

           IF (SMBAtIP>0.0_dp) THEN
             AccuSMB = AccuSMB + SMBAtIP*s
             AccuArea = AccuArea + s
           ELSE
             AblaSMB = AblaSMB + SMBAtIP*s
             AblaArea = AblaArea + s
           END IF
         END DO
         TotalArea = AblaArea + AccuArea
         TotalSMB = AblaSMB + AccuSMB

         ! find the front elevation (min z for Not IceFree)
         IF (.NOT.IceFree) THEN
           Do i=1,n
             ZsLoadsFlux = ZsLoadsFlux + &
               ZsLoads%Values(ZsLoads%Perm(NodeIndexes(i)))*LocalArea(i)/NodeArea(Permutation(NodeIndexes(i)))
           End do

           IF (ANY(NodalZs(1:n)<Hfront)) THEN
             Hfront = MINVAL(NodalZs(1:n))
             nf = MINLOC(NodalZs(1:n), DIM=1)
             Xfront = Solver % Mesh % Nodes % x(NodeIndexes(nf))
             Yfront = Solver % Mesh % Nodes % y(NodeIndexes(nf))
             IF (FoundXRef.and.FoundYRef) THEN
                Dist = sqrt((Solver % Mesh % Nodes % x(NodeIndexes(nf))-XRefFront)**2&
                  +(Solver % Mesh % Nodes % y(NodeIndexes(nf))-YRefFront)**2)
             ENDIF
           END IF
         END IF
      END DO 
      END SUBROUTINE BODY_INTEGRATION

     

      END
