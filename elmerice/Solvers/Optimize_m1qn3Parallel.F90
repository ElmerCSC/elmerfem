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
! *  Authors: 
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!  Optimize a cost function 
!    using quasi-Newton M1QN3 Routine in Reverse Communication
!    Using Euclidian inner product
!
!  Parallel Only   2D/3D
!
!  Need:
!  - Value of the Cost function
!  - Value of the variable to optimize
!  - Value of gradient of cost function with respect to Variable
!      (sum the contribution of each partition shared node)
!  - Optimisation Mask Variable (Optional): name of a mask variable. If
!  mask.lt.0 the variable is considered fixed
!
! => Update the new value of variable to optimize
!
! *****************************************************************************
SUBROUTINE Optimize_m1qn3Parallel( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!  
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
  TYPE(Element_t),POINTER ::  Element
! Variables Beta,DJDbeta and Cost  
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: BetaVar,CostVar,GradVar,MaskVar,TimeVar
  REAL(KIND=dp), POINTER :: BetaValues(:),CostValues(:),GradValues(:),MaskValues(:)
  INTEGER, POINTER :: BetaPerm(:),GradPerm(:),NodeIndexes(:),MaskPerm(:)

  REAL(KIND=dp),allocatable :: x(:),g(:),xx(:),gg(:),xtot(:),gtot(:)
  REAL(KIND=dp) :: f,Normg
  real :: dumy

  integer :: i,j,t,n,NMAX,NActiveNodes,NPoints,ni,ind
  INTEGER :: status(MPI_STATUS_SIZE)
  integer,allocatable :: ActiveNodes(:),NodePerPe(:)
  integer,allocatable :: NewNode(:)
  integer, allocatable :: LocalToGlobalPerm(:),nodePerm(:),TestPerm(:)

  Logical :: FirstVisit=.true.,Firsttime=.true.,Found,UseMask,ComputeNormG=.False.,UnFoundFatal=.TRUE.
  logical,allocatable :: VisitedNode(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName,VarSolName,GradSolName,NormM1QN3,MaskVarName,NormFile
  CHARACTER*10 :: date,temps


!Variables for m1qn3
  external simul_rc,euclid,ctonbe,ctcabe
  character*3 normtype
  REAL(KIND=dp) :: dxmin,df1,epsrel,dzs(1)
  real(kind=dp), allocatable :: dz(:)
  REAL :: rzs(1)
  integer :: imp,io=20,imode(3),omode=-1,niter,nsim,iz(5),ndz,reverse,indic,izs(1)
  integer :: ierr,Npes,ntot
  CHARACTER(LEN=MAX_NAME_LEN) :: IOM1QN3
  logical :: DISbool
!
  save NActiveNodes,Npes,NPoints,ntot
   
  save x,g,xx,gg,xtot,gtot
  save ActiveNodes,NodePerPe

  save TestPerm

  save normtype,dxmin,df1,epsrel,dz,dzs,rzs,imp,io,imode,omode,niter,nsim,iz,ndz,reverse,indic,izs
  save SolverName
  save FirstVisit,Firsttime
  save ComputeNormG,NormFile
  save CostSolName,VarSolName,GradSolName,IOM1QN3


!  Read Constant from sif solver section
      IF(FirstVisit) Then
            FirstVisit=.FALSE.
            CALL DATE_AND_TIME(date,temps)
            WRITE(SolverName, '(A)') 'Optimize_m1qn3Parallel'

       ! Check we have a parallel run
          IF(.NOT.ASSOCIATED(Solver %  Matrix % ParMatrix)) Then
             CALL FATAL(SolverName,'ParMatrix not associated! This solver for parallel only!!')
          End if

            SolverParams => GetSolverParams()

            CostSolName =  GetString( SolverParams,'Cost Variable Name', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Cost Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >CostValue<')
                    WRITE(CostSolName,'(A)') 'CostValue'
                END IF
            VarSolName =  GetString( SolverParams,'Optimized Variable Name', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Optimized Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >Beta<')
                    WRITE(VarSolName,'(A)') 'Beta'
                END IF
            GradSolName =  GetString( SolverParams,'Gradient Variable Name', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Gradient Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >DJDB<')
                    WRITE(GradSolName,'(A)') 'DJDB'
                END IF

          MaskVarName = GetString( SolverParams,'Optimisation Mask Variable',UseMask)
            IF (UseMask) Then
                MaskVar => VariableGet( Solver % Mesh % Variables, MaskVarName,UnFoundFatal=UnFoundFatal) 
                MaskValues => MaskVar % Values 
                MaskPerm => MaskVar % Perm 
            ENDIF

           If (ParEnv % MyPe.EQ.0) then
              NormFile=GetString( SolverParams,'gradient Norm File',Found)
              IF(Found)  Then
                    ComputeNormG=.True.
                    open(io,file=trim(NormFile))
                    write(io,'(a1,a2,a1,a2,a1,a4,5x,a2,a1,a2,a1,a2)')'#',date(5:6),'/',date(7:8),'/',date(1:4), &
                                 temps(1:2),':',temps(3:4),':',temps(5:6)
                    close(io)
              END IF

!!  initialization of m1qn3 variables
            dxmin=GetConstReal( SolverParams,'M1QN3 dxmin', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 dxmin< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >1.e-10<')
                    dxmin=1.e-10
                END IF
            epsrel=GetConstReal( SolverParams,'M1QN3 epsg', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 epsg< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >1.e-06<')
                    epsrel=1.e-6
                END IF
            niter=GetInteger(SolverParams,'M1QN3 niter', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 niter< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >200<')
                    niter=200
                END IF
            nsim=GetInteger(SolverParams,'M1QN3 nsim', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 nsim< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >200<')
                    nsim=200
                END IF
            imp=GetInteger(SolverParams,'M1QN3 impres', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 impres< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >5<')
                    imp=5
                END IF
            DISbool=GetLogical( SolverParams, 'M1QN3 DIS Mode', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 DIS Mode< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >FALSE<')
                    DISbool=.False.
                END IF
                if(DISbool) then
                    imode(1)=0 !DIS Mode
                else
                    imode(1)=1 !SIS Mode
                End if
            df1=GetConstReal( SolverParams,'M1QN3 df1', Found)
                IF(.NOT.Found) THEN
                   CALL WARN(SolverName,'Keyword >M1QN3 df1< not found  in section >Solver<')
                   CALL WARN(SolverName,'Taking default value >0.2<')
                   df1=0.2
                End if
                CostVar => VariableGet( Solver % Mesh % Variables, CostSolName,UnFoundFatal=UnFoundFatal)
                CostValues => CostVar % Values
                df1=CostValues(1)*df1
             NormM1QN3 = GetString( SolverParams,'M1QN3 normtype', Found)
                 IF((.NOT.Found).AND.((NormM1QN3(1:3).ne.'dfn').OR.(NormM1QN3(1:3).ne.'sup') &
                     .OR.(NormM1QN3(1:3).ne.'two'))) THEN
                       CALL WARN(SolverName,'Keyword >M1QN3 normtype< not good in section >Solver<')
                       CALL WARN(SolverName,'Taking default value >dfn<')
                       PRINT *,'M1QN3 normtype  ',NormM1QN3(1:3)
                       normtype = 'dfn'
                  ELSE
                       PRINT *,'M1QN3 normtype  ',NormM1QN3(1:3)
                       normtype = NormM1QN3(1:3)
                  END IF

              IOM1QN3 = GetString( SolverParams,'M1QN3 OutputFile', Found)
                 IF(.NOT.Found) THEN
                       CALL WARN(SolverName,'Keyword >M1QN3 OutputFile< not found  in section >Solver<')
                       CALL WARN(SolverName,'Taking default value >M1QN3.out<')
                       WRITE(IOM1QN3,'(A)') 'M1QN3.out'
                 END IF
                 open(io,file=trim(IOM1QN3))
                    write(io,*) '******** M1QN3 Output file ************'
                    write(io,'(a2,a1,a2,a1,a4,5x,a2,a1,a2,a1,a2)') date(5:6),'/',date(7:8),'/',date(1:4), &
                                 temps(1:2),':',temps(3:4),':',temps(5:6)
                    write(io,*) '*****************************************'
                 close(io)
              ndz=GetInteger( SolverParams,'M1QN3 ndz', Found)
                  IF(.NOT.Found) THEN
                       CALL WARN(SolverName,'Keyword >M1QN3 ndz< not found  in section >Solver<')
                       CALL WARN(SolverName,'Taking default value >5< update')
                       ndz=5
                   END IF

                    imode(2)=0 
                    imode(3)=0 
                    reverse=1 
                    omode=-1 
                    dzs=0.0 
                    rzs=0.0
                    izs=0

            End if !ParEnv % MyPe=0

        End if ! FirtsVisit


! Omode from previous iter; if > 0 m1qn3 has terminated => return 
     IF (omode.gt.0) then 
             WRITE(Message,'(a,I1)') 'm1qn3 finished; omode=',omode 
             CALL Info(SolverName, Message, Level=1) 
             return  
     End if

!  Get Variables CostValue, Beta and DJDBeta
    CostVar => VariableGet( Solver % Mesh % Variables, CostSolName,UnFoundFatal=UnFoundFatal)
    CostValues => CostVar % Values 
    f=CostValues(1)

     BetaVar => VariableGet( Solver % Mesh % Variables, VarSolName,UnFoundFatal=UnFoundFatal) 
     BetaValues => BetaVar % Values 
     BetaPerm => BetaVar % Perm 

     GradVar => VariableGet( Solver % Mesh % Variables, GradSolName,UnFoundFatal=UnFoundFatal) 
     GradValues   => GradVar % Values 
     GradPerm => GradVar % Perm 

! Do some allocation etc if first iteration
  If (Firsttime) then 
          
     Firsttime = .False.

     NMAX=Solver % Mesh % NumberOfNodes
     allocate(VisitedNode(NMAX),NewNode(NMAX))

!!!!!!!!!!!!find active nodes 
      VisitedNode=.false.  
      NewNode=-1

      NActiveNodes=0 
      DO t=1,Solver % NumberOfActiveElements
         Element => GetActiveElement(t)
         n = GetElementNOFNodes()
         NodeIndexes => Element % NodeIndexes
         Do i=1,n
           if (VisitedNode(NodeIndexes(i))) then
                   cycle
           else
              VisitedNode(NodeIndexes(i))=.true.
              IF (UseMask) Then
                       IF (MaskValues(MaskPerm(NodeIndexes(i))).lt.0) cycle
              END IF
              NActiveNodes=NActiveNodes+1
              NewNode(NActiveNodes)=NodeIndexes(i)
           endif
        End do
      End do

     if (NActiveNodes.eq.0) THEN
            WRITE(Message,'(A)') 'NActiveNodes = 0 !!!'
            CALL FATAL(SolverName,Message)
     End if
  
    allocate(ActiveNodes(NActiveNodes),LocalToGlobalPerm(NActiveNodes),x(NActiveNodes),g(NActiveNodes))
    ActiveNodes(1:NActiveNodes)=NewNode(1:NActiveNodes)

    deallocate(VisitedNode,NewNode)

!! Gather number of active nodes in each partition and compute total number of
!active nodes in partion 0
    Npes=Solver %  Matrix % ParMatrix % ParEnv % PEs
    allocate(NodePerPe(Npes))

    call MPI_Gather(NActiveNodes,1,MPI_Integer,NodePerPe,1,MPI_Integer,0,ELMER_COMM_WORLD,ierr)

    if (Solver %  Matrix % ParMatrix % ParEnv % MyPE.eq.0) then
          ntot=0
          Do i=1,Npes
               ntot=ntot+NodePerPe(i)
          End do
          allocate(xtot(ntot),gtot(ntot))
          allocate(NodePerm(ntot))
    End if

! Send global node  numbering to partition 0

    LocalToGlobalPerm(1:NActiveNodes)=Model % Mesh % ParallelInfo % GlobalDOFs(ActiveNodes(1:NActiveNodes))

   if (Solver %  Matrix % ParMatrix % ParEnv % MyPE .ne.0) then
             call MPI_BSEND(LocalToGlobalPerm(1),NActiveNodes,MPI_INTEGER,0,8001,ELMER_COMM_WORLD,ierr)
   else
           NodePerm(1:NActiveNodes)=LocalToGlobalPerm(1:NActiveNodes)
           ni=1+NActiveNodes
           Do i=2,Npes
             call   MPI_RECV(NodePerm(ni),NodePerPe(i),MPI_INTEGER,i-1,8001,ELMER_COMM_WORLD, status, ierr )
             ni=ni+NodePerPe(i)
           End do

! Create a permutation table from NodePerm
           allocate(TestPerm(ntot))
           ind=1
           TestPerm(1)=ind
           Do i=2,ntot
              Do j=1,i-1
                 if (NodePerm(j).eq.NodePerm(i)) exit
               End do
               if (j.eq.i) then
                       ind=ind+1
                       TestPerm(i)=ind
               else
                       TestPerm(i)=TestPerm(j)
               end if
           End do

           NPoints=ind
           allocate(xx(Npoints),gg(Npoints))
           deallocate(NodePerm,LocalToGlobalPerm)

 ! M1QN3 allocation of dz function of Npoints nd requested number of updates
           IF (DISbool) then
                   ndz=4*NPoints+ndz*(2*NPoints+1)+10
           else
                   ndz=3*NPoints+ndz*(2*NPoints+1)+10
           end if
           allocate(dz(ndz))

    End if
   
  END IF


     x(1:NActiveNodes)=BetaValues(BetaPerm(ActiveNodes(1:NActiveNodes)))
     g(1:NActiveNodes)=GradValues(GradPerm(ActiveNodes(1:NActiveNodes)))

    ! Send variables to partition 0
    ! and receive results from partion 0

     if (Solver %  Matrix % ParMatrix % ParEnv % MyPE .ne.0) then

                     call MPI_SEND(x(1),NActiveNodes,MPI_DOUBLE_PRECISION,0,8003,ELMER_COMM_WORLD,ierr)
                     call MPI_SEND(g(1),NActiveNodes,MPI_DOUBLE_PRECISION,0,8004,ELMER_COMM_WORLD,ierr)
                     call MPI_RECV(x(1),NActiveNodes,MPI_DOUBLE_PRECISION,0,8005,ELMER_COMM_WORLD,status, ierr )
                     call MPI_RECV(omode,1,MPI_Integer,0,8006,ELMER_COMM_WORLD,status,ierr )

                     ! Update Beta Values 
                     BetaValues(BetaPerm(ActiveNodes(1:NActiveNodes)))=x(1:NActiveNodes)
     else
                     xtot(1:NActiveNodes)=x(1:NActiveNodes)
                     gtot(1:NActiveNodes)=g(1:NActiveNodes)
                     ni=1+NActiveNodes
                     Do i=2,Npes
                       call MPI_RECV(xtot(ni),NodePerPe(i),MPI_DOUBLE_PRECISION,i-1,8003,ELMER_COMM_WORLD, status, ierr )
                       call MPI_RECV(gtot(ni),NodePerPe(i),MPI_DOUBLE_PRECISION,i-1,8004,ELMER_COMM_WORLD, status, ierr )
                       ni=ni+NodePerPe(i)
                     End do
                     
                     xx=0.0
                     gg=0.0
                     Do i=1,ntot
                         xx(TestPerm(i))=xtot(i)  ! same Beta Value for same node
                         gg(TestPerm(i))=gg(TestPerm(i))+gtot(i)  ! gather the contribution to DJDB 
                     End do 

                     If (ComputeNormG) then
                             TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )
                             Normg=0.0_dp

                             Do i=1,NPoints
                                Normg=Normg+gg(i)*gg(i)
                             End do
                             open(io,file=trim(NormFile),position='append')
                              write(io,'(e13.5,2x,e15.8)') TimeVar % Values(1),sqrt(Normg)
                             close(io)
                    End if

            ! go to minimization
            open(io,file=trim(IOM1QN3),position='append')
            call m1qn3 (simul_rc,Euclid,ctonbe,ctcabe,NPoints,xx,f,gg,dxmin,df1, &
                        epsrel,normtype,imp,io,imode,omode,niter,nsim,iz, &
                        dz,ndz,reverse,indic,izs,rzs,dzs)

            close(io)
            WRITE(Message,'(a,e15.8,x,I2)') 'm1qn3: Cost,omode= ',f,omode
            CALL Info(SolverName, Message, Level=3)

            ! Put new Beta Value in xtot and send to each partition
            xtot=0.0
            Do i=1,ntot
               xtot(i)=xx(TestPerm(i))
            End do 

            ! Update Beta Values 
            BetaValues(BetaPerm(ActiveNodes(1:NActiveNodes)))=xtot(1:NActiveNodes)
                      
            ni=1+NActiveNodes
            Do i=2,Npes
                  call MPI_SEND(xtot(ni),NodePerPe(i),MPI_DOUBLE_PRECISION,i-1,8005,ELMER_COMM_WORLD,ierr)
                  call MPI_SEND(omode,1,MPI_Integer,i-1,8006,ELMER_COMM_WORLD,ierr)
                  ni=ni+NodePerPe(i)
           End do
  endif

   Return
!------------------------------------------------------------------------------
END SUBROUTINE Optimize_m1qn3Parallel
!------------------------------------------------------------------------------


