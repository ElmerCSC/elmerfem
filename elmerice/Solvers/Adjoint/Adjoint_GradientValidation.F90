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
! *  Authors: f. Gillet-Chaulet (IGE, Grenoble,France)
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: April 2020
! * 
! *****************************************************************************
SUBROUTINE Adjoint_GradientValidation_init0(Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  CHARACTER(LEN=MAX_NAME_LEN) :: Name

  Name = ListGetString( Solver % Values, 'Equation',UnFoundFatal=.TRUE.)
  CALL ListAddNewString( Solver % Values,'Variable',&
          '-nooutput '//TRIM(Name)//'_var')
  CALL ListAddLogical(Solver % Values, 'Optimize Bandwidth',.FALSE.)
  CALL ListAddInteger(Solver % Values, 'Nonlinear System Norm Degree',0)
END SUBROUTINE Adjoint_GradientValidation_init0
! *****************************************************************************
SUBROUTINE Adjoint_GradientValidation ( Model,Solver,dt,TransientSimulation )
! *****************************************************************************
!   Compare the total derivative of the cost function computed as:
!     (1) dJ=P.G  where P is a perturbation vector of the variable of interest
!                     G is the gradient of the cost function computed by an inverse method
!     (2) [J(V+hP)-J(V)]/h  : forward finite difference computation of the derivative
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
!           Perturbed Variable Name !optional: take -g if not given
!           Gradient Variable Name
!           Result File
!
!  Output: in result File: h , abs((1)-(2))/(1) , (1), (2)
!
!
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

  integer :: i,c,t,n,NMAX,NActiveNodes
  integer :: ierr
  integer,allocatable :: ActiveNodes(:)
  integer,allocatable :: NewNode(:)
  integer,parameter :: io=20
  integer :: MyPe=-1
  integer,SAVE :: DOFs

  Logical :: FirstVisit=.true.,Found,UnFoundFatal=.TRUE.
  Logical :: Parallel
  Logical :: haveP
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

          CostSolName =  ListGetString( SolverParams,'Cost Variable Name', UnFoundFatal=.TRUE.)
          CostVar => VariableGet( Solver % Mesh % Variables, CostSolName,UnFoundFatal=UnFoundFatal)
          CostValues => CostVar % Values
 
          VarSolName =  ListGetString( SolverParams,'Optimized Variable Name', UnFoundFatal=.TRUE.)
          Var => VariableGet( Solver % Mesh % Variables, VarSolName,UnFoundFatal=UnFoundFatal) 
          Values => Var % Values 
          Perm => Var % Perm
          DOFs = Var % DOFs

          GradSolName =  ListGetString( SolverParams,'Gradient Variable Name', UnFoundFatal=.TRUE.)
          GradVar => VariableGet( Solver % Mesh % Variables, GradSolName,UnFoundFatal=UnFoundFatal) 
          GradValues   => GradVar % Values 
          GradPerm => GradVar % Perm 
          IF (GradVar%DOFs.NE.DOFs) &
                  CALL FATAL(SolverName,'DOFs not corresponding for Gradient Variable')

          PSolName =  ListGetString( SolverParams,'Perturbation Variable Name', Found=HaveP)
          IF (HAVEP) THEN
               PVar => VariableGet( Solver % Mesh % Variables, PSolName,UnFoundFatal=UnFoundFatal) 
               PValues => PVar % Values 
               PPerm => PVar % Perm
               IF (PVar%DOFs.NE.DOFs) &
                  CALL FATAL(SolverName,'DOFs not corresponding for Perturbation Variable')
          ENDIF

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

           allocate(ActiveNodes(NActiveNodes),x(DOFS*NActiveNodes),xp(DOFs*NActiveNodes),g(DOFs*NActiveNodes))
           ActiveNodes(1:NActiveNodes)=NewNode(1:NActiveNodes)

           deallocate(VisitedNode,NewNode)

!!!!!!!  Solver Params


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
            DO i=1,NActiveNodes
             DO c=1,DOFs
              x(DOFs*(i-1)+c)=Values(DOFs*(Perm(ActiveNodes(i))-1)+c)
              g(DOFs*(i-1)+c)=GradValues(DOFs *(GradPerm(ActiveNodes(i))-1)+c)
              IF (HaveP) THEN
               xp(DOFs*(i-1)+c)=PValues(DOFs*(PPerm(ActiveNodes(i))-1)+c)
              ELSE
               xp(DOFs*(i-1)+c)=-g(DOFs*(i-1)+c)
              ENDIF
             END DO
            END DO

             !!!!!!  total derivative from gradient
             dJ=0._dp
             dJ=SUM(xp(:)*g(:)) 
             deallocate(g)

             IF (Parallel) THEN
                CALL MPI_ALLREDUCE(dJ,dJ,1,&
                     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
             ENDIF

             !!!!! Store cost value at first iteration
             J0=CostValues(1)
               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !!! new values for x=x+h*xp
             h=1.0_dp
             Do i=1,NActiveNodes
              Do c=1,DOFs
              Values(DOFs*(Perm(ActiveNodes(i))-1)+c)=x(DOFs*(i-1)+c)+h*xp(DOFs*(i-1)+c)
              End do
             END DO


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
          Do c=1,DOFs
              Values(DOFs*(Perm(ActiveNodes(i))-1)+c)=x(DOFs*(i-1)+c)+h*xp(DOFs*(i-1)+c)
          End do
        END DO

   Return
!------------------------------------------------------------------------------
END SUBROUTINE Adjoint_GradientValidation
!------------------------------------------------------------------------------


