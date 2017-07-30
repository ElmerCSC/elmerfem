!   Compare the total derivative of the cost function computed as:
!     (1) dJ=P.G  where P is a perturbation vector of the variable of interest
!                     G is the gradient of the cost function computed by an inverse method
!     (2) [J(V+hP)-J(V)]/h  : forward finite diffrence computation of the derivative
!                             V is the variable of interest
!                             h is the step size 
!
!
!  Compute (1) from at the first iteration and update V=Vini+hP, h=1
!  Compute (2) for all the other iteration with h^i+1=h^i/2
!
!  Serial/parallel   2D/3D
!
!  Keyword in Solver section of the .sif:
!           Cost Variable Name
!           Optimized Variable Name
!           Perturbed Variable Name
!           Gradient Variable Name
!           Result File
!
!  Output: in result File: h , abs((1)-(2))/(1) , (1), (2)
!
!
! *****************************************************************************
SUBROUTINE GradientValidation ( Model,Solver,dt,TransientSimulation )
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
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: Var,PVar,CostVar,GradVar
  REAL(KIND=dp), POINTER :: Values(:),PValues(:),CostValues(:),GradValues(:)
  INTEGER, POINTER :: Perm(:),PPerm(:),GradPerm(:)
  INTEGER, POINTER :: NodeIndexes(:)

  REAL(KIND=dp),allocatable :: x(:),xp(:),g(:)
  REAL(KIND=dp) :: J,J0,h,dJ,dJd

  integer :: i,t,n,NMAX,NActiveNodes
  integer :: ierr
  integer,allocatable :: ActiveNodes(:)
  integer,allocatable :: NewNode(:)
  integer,parameter :: io=20
  integer :: MyPe=-1

  Logical :: FirstVisit=.true.,Found,UnFoundFatal
  Logical :: Parallel
  logical,allocatable :: VisitedNode(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName,VarSolName,GradSolName,PSolName,ResultFile
  CHARACTER*10 :: date,temps

!
  save FirstVisit
  save SolverName
  save CostSolName,VarSolName,GradSolName,ResultFile
  save ActiveNodes,NActiveNodes
  save x,xp
  save J0,h,dJ
  save MyPe
  save Parallel


!  Read Constant from sif solver section
      IF(FirstVisit) Then

            WRITE(SolverName, '(A)') 'GradientValidation'

           ! Check we have a parallel run
           Parallel = .FALSE.
           IF(ASSOCIATED(Solver %  Matrix % ParMatrix)) Then
             IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 ) Parallel =.True.
             MyPe=ParEnv % MyPe
           End if


          SolverParams => GetSolverParams()

!!!!!!!!!!!!find active nodes 
           NMAX=Solver % Mesh % NumberOfNodes
           allocate(VisitedNode(NMAX),NewNode(NMAX))
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
                     NActiveNodes=NActiveNodes+1
                     NewNode(NActiveNodes)=NodeIndexes(i)
                 endif
             End do
           End do

           if (NActiveNodes.eq.0) THEN
              WRITE(Message,'(A)') 'NActiveNodes = 0 !!!'
              CALL FATAL(SolverName,Message)
           End if

           allocate(ActiveNodes(NActiveNodes),x(NActiveNodes),xp(NActiveNodes),g(NActiveNodes))
           ActiveNodes(1:NActiveNodes)=NewNode(1:NActiveNodes)

           deallocate(VisitedNode,NewNode)

!!!!!!!  Solver Params

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
            PSolName =  GetString( SolverParams,'Perturbed Variable Name', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Perturbed Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >BetaP<')
                    WRITE(VarSolName,'(A)') 'BetaP'
                END IF
            GradSolName =  GetString( SolverParams,'Gradient Variable Name', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Gradient Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >DJDB<')
                    WRITE(GradSolName,'(A)') 'DJDB'
                END IF
            ResultFile=GetString( SolverParams,'Result File',Found)
            IF(Found)  Then
               IF ((Parallel.AND.(MyPe.EQ.0)).OR.(.NOT.Parallel)) Then
                     open(io,file=trim(ResultFile))
                     CALL DATE_AND_TIME(date,temps)
                     write(io,'(a1,a2,a1,a2,a1,a4,5x,a2,a1,a2,a1,a2)')'#',date(5:6),'/',date(7:8),'/',date(1:4), &
                                 temps(1:2),':',temps(3:4),':',temps(5:6)
                     write(io,'(A)') '# step size, relative error, Adjoint total der., FD total der.'
                     close(io)
               ENDIF
            ELSE
                    CALL FATAL(SolverName,'Keyword <Result File> Not Found')
            ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             CostVar => VariableGet( Solver % Mesh % Variables, CostSolName,UnFoundFatal=UnFoundFatal)
             CostValues => CostVar % Values 

             Var => VariableGet( Solver % Mesh % Variables, VarSolName,UnFoundFatal=UnFoundFatal) 
             Values => Var % Values 
             Perm => Var % Perm 
             
             PVar => VariableGet( Solver % Mesh % Variables, PSolName,UnFoundFatal=UnFoundFatal) 
             PValues => PVar % Values 
             PPerm => PVar % Perm 

             GradVar => VariableGet( Solver % Mesh % Variables, GradSolName,UnFoundFatal=UnFoundFatal)
             GradValues   => GradVar % Values 
             GradPerm => GradVar % Perm 


             x(1:NActiveNodes)=Values(Perm(ActiveNodes(1:NActiveNodes)))
             g(1:NActiveNodes)=GradValues(GradPerm(ActiveNodes(1:NActiveNodes)))
             xp(1:NActiveNodes)=PValues(PPerm(ActiveNodes(1:NActiveNodes)))

             !!!!!! On calcul la dervivee totale à partir de la valeur du
             !gradient
             dJ=0._dp
             Do i=1,NActiveNodes
                dJ=dJ+xp(i)*g(i)
             End do
             deallocate(g)

             IF (Parallel) THEN
                CALL MPI_ALLREDUCE(dJ,dJ,1,&
                     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
             ENDIF


             !!!!! On sauve la valeur de J0 pour le calcul diffrence fini.
             J0=CostValues(1)
               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             h=1.0_dp
             Do i=1,NActiveNodes
                Values(Perm(ActiveNodes(i)))=x(i)+h*xp(i)
             End do


            FirstVisit=.FALSE.
            Return
        End if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CostVar => VariableGet( Solver % Mesh % Variables, CostSolName,UnFoundFatal=UnFoundFatal)
        CostValues => CostVar % Values 

        Var => VariableGet( Solver % Mesh % Variables, VarSolName,UnFoundFatal=UnFoundFatal) 
        Values => Var % Values 
        Perm => Var % Perm 

       
        J=CostValues(1)

        dJd=(J-J0)/h

        IF (Parallel) MyPe=ParEnv % MyPe
        IF ((Parallel.AND.(MyPe.EQ.0)).OR.(.NOT.Parallel)) Then
           open(io,file=trim(ResultFile),position='append')
                write(io,'(4(e15.8,2x))') h,abs(dJ-dJd)/abs(dJ),dJ,dJd
           close(io)
        ENDIF

        h=h/2.0_dp
        Do i=1,NActiveNodes
           Values(Perm(ActiveNodes(i)))=x(i)+h*xp(i)
        End Do

   Return
!------------------------------------------------------------------------------
END SUBROUTINE GradientValidation
!------------------------------------------------------------------------------


