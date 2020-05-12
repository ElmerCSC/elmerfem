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
! *  Original Date: 
! * 
! *****************************************************************************
SUBROUTINE Adjoint_CostContSolver_init0(Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  CHARACTER(LEN=MAX_NAME_LEN) :: Name
  
  Name = ListGetString( Solver % Values, 'Equation',UnFoundFatal=.TRUE.)
  CALL ListAddNewString( Solver % Values,'Variable',&
          '-nooutput '//TRIM(Name)//'_var')
  CALL ListAddLogical(Solver % Values, 'Optimize Bandwidth',.FALSE.)
  CALL ListAddInteger(Solver % Values, 'Nonlinear System Norm Degree',0)
END SUBROUTINE Adjoint_CostContSolver_init0
!----------------------------------------------------------------------
! *****************************************************************************
SUBROUTINE Adjoint_CostContSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!Compute a Cost function  as: Integral Node_Cost ds
!
!  OUTPUT are : J and xb (sensitivity of J w.r.t. u)

!    !! Be careful this solver will reset the cost and xb to 0;
!     so it has to be used as the first cost solver in an inverse problem sequence
!
!  see documentation under : elmerice/Solvers/Documentation/Adjoint_CostContSolver.md
!******************************************************************************
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!  
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="Adjoint_CostContSolver"
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: DefaultCostFile = 'CostOfT.dat'
  CHARACTER(LEN=MAX_NAME_LEN),SAVE :: CostFile
  CHARACTER(LEN=MAX_NAME_LEN),SAVE :: CostSolName
  CHARACTER(LEN=MAX_NAME_LEN),SAVE :: VbName
  CHARACTER(LEN=MAX_NAME_LEN) :: DerName

  TYPE(Variable_t), POINTER :: TimeVar,CostVar
  TYPE(Variable_t), POINTER :: VbSol
  REAL(KIND=dp), POINTER :: Vb(:)
  INTEGER, POINTER :: VbPerm(:)
  INTEGER,SAVE :: VDOFs

  TYPE(ValueList_t), POINTER :: SolverParams,BodyForce

  TYPE(Element_t),POINTER ::  Element
  TYPE(Nodes_t),SAVE :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  INTEGER, POINTER :: NodeIndexes(:)
  REAL(KIND=dp),ALLOCATABLE,SAVE :: Basis(:),dBasisdx(:,:)
  REAL(KIND=dp) :: u,v,w,SqrtElementMetric,x,s
  INTEGER :: n

  LOGICAL,SAVE :: Firsttime=.TRUE.
  LOGICAL,SAVE :: Parallel
  LOGICAL :: BoundarySolver

  LOGICAL :: Found,stat
  integer :: i,j,k,t
  INTEGER :: ierr
  real(kind=dp) :: Cost,Cost_S
  real(kind=dp) :: Area,Area_S
  real(kind=dp) :: coeff

  REAL(KIND=dp),ALLOCATABLE,SAVE :: NodeCost(:)
  REAL(KIND=dp),ALLOCATABLE,SAVE :: NodeCostb(:),NodeCost_der(:,:)

  INTEGER, SAVE :: DIM

  CHARACTER*10 :: date,temps


  SolverParams => GetSolverParams()

  If (Firsttime) then
    N = model % MaxElementNodes
    allocate(ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N))
    allocate(Basis(N),dBasisdx(N,3))

!!!!!!! Check for parallel run 
    Parallel = .FALSE.
    IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
      IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
          Parallel = .TRUE.
      END IF
    END IF

    !! check if we are on a boundary or in the bulk
    BoundarySolver = ( Solver % ActiveElements(1) > Model % Mesh % NumberOfBulkElements )
    IF (BoundarySolver) THEN 
      DIM = CoordinateSystemDimension() - 1
    ELSE
      DIM = CoordinateSystemDimension()
    ENDIF

!!!!!!!!!!! get Solver Variables
  CostFile = ListGetString(Solver % Values,'Cost Filename',Found )
    IF (.NOT. Found) CostFile = DefaultCostFile
    CALL DATE_AND_TIME(date,temps)
    If (Parallel) then
      if (ParEnv % MyPe.EQ.0) then
        OPEN (12, FILE=CostFile)
             write(12,1000) date(5:6),date(7:8),date(1:4),temps(1:2),temps(3:4),temps(5:6)
             write(12,'(A)') '#, 1.0'
             write(12,'(A)') '# iter, J0, sqrt(2J0/Area)'
        CLOSE(12)
      End if
    Else
        OPEN (12, FILE=CostFile)
             write(12,1000) date(5:6),date(7:8),date(1:4),temps(1:2),temps(3:4),temps(5:6)
             write(12,'(A)') '#, 1.0'
             write(12,'(A)') '# iter, J0, sqrt(2J0/Area)'
        CLOSE(12)
    End if

   CostSolName =  GetString( SolverParams,'Cost Variable Name', Found)
   IF(.NOT.Found) THEN
         CALL WARN(SolverName,'Keyword >Cost Variable Name< not found  in section >Solver<')
         CALL WARN(SolverName,'Taking default value >CostValue<')
         WRITE(CostSolName,'(A)') 'CostValue'
   END IF

   VbName = ListGetString(SolverParams,'Sensitivity Variable Name', UnFoundFatal=.TRUE.)
   VbSol => VariableGet( Solver % Mesh % Variables,TRIM(VbName), UnFoundFatal=.TRUE.  )
   VDOFs = VbSol % DOFs
   allocate(NodeCost(N),NodeCostb(N),NodeCost_der(VDOFs,N))
  
  !!! End of First visit
    Firsttime=.false.
  Endif

    VbSol => VariableGet( Solver % Mesh % Variables,TRIM(VbName), UnFoundFatal=.TRUE.  )
    Vb => VbSol % Values
    VbPerm => VbSol % Perm
    Vb=0.0_dp

    Cost=0._dp
    Area=0._dp

    DO t=1,Solver % NumberOfActiveElements
       Element => GetActiveElement(t)
       IF (CheckPassiveElement(Element)) CYCLE
       n = GetElementNOFNodes()

       NodeIndexes => Element % NodeIndexes

 ! set coords of highest occurring dimension to zero (to get correct path element)
        !-------------------------------------------------------------------------------
        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        IF (DIM == 1) THEN !1D 
           ElementNodes % y(1:n) = 0.0_dp
           ElementNodes % z(1:n) = 0.0_dp
        ELSE IF (DIM == 2) THEN !2D 
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = 0.0_dp
        ELSE IF (DIM == 3) THEN 
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
        END IF

        BodyForce => GetBodyForce(Element)

        NodeCost=0.0_dp
        NodeCost(1:n) = ListGetReal( BodyForce, 'Adjoint Cost', n, NodeIndexes, UnFoundFatal=.TRUE.)
        NodeCost_der=0.0_dp
        
        IF (VDOFs.EQ.1) THEN
          DerName='Adjoint Cost der'
          IF(.NOT. ListCheckPresent( BodyForce,TRIM(DerName))) &
                  DerName='Adjoint Cost der ' // I2S(1)
           NodeCost_der(1,1:n)=ListGetReal( BodyForce,DerName, n, NodeIndexes, UnFoundFatal=.TRUE.)
        ELSE
          DO k=1,VDOFs
            DerName=ComponentName("Adjoint Cost der",k)
            NodeCost_der(k,1:n)=ListGetReal( BodyForce,TRIM(DerName), n, NodeIndexes, Found)
            IF (.NOT.Found) NodeCost_der(k,1:n)=0._dp
          END DO
        END IF

!------------------------------------------------------------------------------
!    Numerical integration
!------------------------------------------------------------------------------
        IntegStuff = GaussPoints( Element )

        NodeCostb=0.0_dp

        DO i=1,IntegStuff % n
          U = IntegStuff % u(i)
          V = IntegStuff % v(i)
          W = IntegStuff % w(i)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
          stat = ElementInfo( Element,ElementNodes,U,V,W,SqrtElementMetric, &
              Basis,dBasisdx )

          x = SUM( ElementNodes % x(1:n) * Basis(1:n) )
          s = 1.0d0

          IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
             s = 2.0d0 * PI * x 
          END IF
          s = s * SqrtElementMetric * IntegStuff % s(i)
          
          coeff = SUM(NodeCost(1:n)  * Basis(1:n))
          Cost=Cost+coeff*s
          Area=Area+s

          NodeCostb(1:n)=NodeCostb(1:n) + s*Basis(1:n)

        End do

         Do j=1,n
          Do i=1,VDOFs
            k=(VbPerm(NodeIndexes(j))-1)*VDOFs+i
            Vb(k)=Vb(k)+NodeCostb(j)*NodeCost_der(i,j)
          End Do
         End Do
    End do

    TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )

    IF (Parallel) THEN
     CALL MPI_ALLREDUCE(Cost,Cost_S,1,&
            MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)

     CALL MPI_ALLREDUCE(Area,Area_S,1,&
            MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)

     CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
     IF (ASSOCIATED(CostVar)) THEN
         CostVar % Values(1)=Cost_S
     END IF
     IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) then
        OPEN (12, FILE=CostFile,POSITION='APPEND')
           write(12,'(3(ES20.11E3))') TimeVar % Values(1),Cost_S,sqrt(2*Cost_S/Area_S)
        CLOSE(12)
     End if
   ELSE
     CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
     IF (ASSOCIATED(CostVar)) THEN
        CostVar % Values(1)=Cost
     END IF
     OPEN (12, FILE=CostFile,POSITION='APPEND')
        write(12,'(3(ES20.11E3))') TimeVar % Values(1),Cost,sqrt(2*Cost/Area)
     close(12)
     Cost_S=Cost
   END IF

   Solver % Variable % Values(1)=Cost_S
   
   Return

 1000  format('#date,time,',a2,'/',a2,'/',a4,',',a2,':',a2,':',a2)
!------------------------------------------------------------------------------
END SUBROUTINE Adjoint_CostContSolver
!------------------------------------------------------------------------------
! *****************************************************************************
