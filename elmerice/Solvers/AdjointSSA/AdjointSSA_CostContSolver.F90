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
! *  Authors: f. Gillet-Chaulet (LGGE, Grenoble,France)
! *  Email:   gillet-chaulet@lgge.obs.ujf-grenoble.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!
! *****************************************************************************
SUBROUTINE AdjointSSA_CostContSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!Compute the Cost function of the SSA Adjoint inverse Problem  as: Integral Node_Cost ds
!
!   Serial/Parallel    2D/3D
!
!     OUTPUT are : J and DJDu (==Velocityb variable used as forcing of the SSA adjoint problem)
!
!    !! Be careful this solver will reset the cost and DJDu to 0; so it has to
!    be used as the first cost solver if regularistaion of flux divergence cost
!    solvers are used in the simulation!!
!
!
!     INPUT PARAMETERS are:
!
!      In solver section:
!               Problem Dimension = Integer (default = Coordinate system dimension)
!               Cost Filename = File (default = 'CostOfT.dat')
!               Cost Variable Name = String (default= 'CostValue')
!
!      Variables
!                Velocityb (forcing for the adjoint pb;DOFs== Pb dimension)
!
!     In body Forces:
!                'Adjoint Cost = Real ...' : The nodal value of the cost
!                'Adjoint Cost der 1 = Real ...' : The derivative of 'Adjoint Cost' w.r.t. u-velocity
!                'Adjoint Cost der 2 = Real ...' : The derivative of 'Adjoint Cost' w.r.t. v-velocity
!
!******************************************************************************
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!  
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: DefaultCostFile = 'CostOfT.dat'
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,CostFile
  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: TimeVar,CostVar
  TYPE(Variable_t), POINTER :: VelocitybSol
  TYPE(ValueList_t), POINTER :: SolverParams,BodyForce
  TYPE(Nodes_t) :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  REAL(KIND=dp), POINTER :: Vb(:)
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER, POINTER :: VbPerm(:)
  Logical :: Firsttime=.true.,Found,Parallel,stat,Gotit
  integer :: i,j,k,l,t,n,NMAX,DIM,ierr,c
  real(kind=dp) :: Cost,Cost_surf,Cost_S
  real(kind=dp) :: u,v,w,s,coeff,SqrtElementMetric,x
  REAL(KIND=dp) :: NodeCost(Model % MaxElementNodes)
  REAL(KIND=dp) :: Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)
  REAL(KIND=dp) :: NodeCostb(Model % MaxElementNodes),NodeCost_der(3,Model %MaxElementNodes)
  CHARACTER*10 :: date,temps

  save Firsttime,Parallel 
  save SolverName,CostSolName,CostFile
  save ElementNodes


  SolverParams => GetSolverParams()
  DIM=GetInteger(SolverParams ,'Problem Dimension',Found)
  If (.NOT.Found) then
     CALL WARN(SolverName,'Keyword >Problem Dimension< not found, assume DIM = CoordinateSystemDimension()')
     DIM = CoordinateSystemDimension()
  Endif

  If (Firsttime) then
    N = model % MaxElementNodes
    allocate(ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N))

!!!!!!! Check for parallel run 
    Parallel = .FALSE.
    IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
            IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
                    Parallel = .TRUE.
            END IF
    END IF

     WRITE(SolverName, '(A)') 'CostSolver_Adjoint'

!!!!!!!!!!! get Solver Variables

  CostFile = ListGetString(Solver % Values,'Cost Filename',Found )
    IF (.NOT. Found) CostFile = DefaultCostFile
    CALL DATE_AND_TIME(date,temps)
    If (Parallel) then
        if (ParEnv % MyPe.EQ.0) then
           OPEN (12, FILE=CostFile)
                   write(12,1000) date(5:6),date(7:8),date(1:4),temps(1:2),temps(3:4),temps(5:6)
                   write(12,'(A)') '#, 1.0'
                   write(12,'(A)') '# iter, J0'
           CLOSE(12)
         End if
    Else
           OPEN (12, FILE=CostFile)
                   write(12,1000) date(5:6),date(7:8),date(1:4),temps(1:2),temps(3:4),temps(5:6)
                   write(12,'(A)') '#, 1.0'
                   write(12,'(A)') '# iter, J0'
           CLOSE(12)
    End if

   CostSolName =  GetString( SolverParams,'Cost Variable Name', Found)
          IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Cost Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >CostValue<')
                    WRITE(CostSolName,'(A)') 'CostValue'
          END IF
  
  !!! End of First visit
    Firsttime=.false.
  Endif


    VelocitybSol => VariableGet( Solver % Mesh % Variables, 'Velocityb'  )
    IF ( ASSOCIATED( VelocitybSol ) ) THEN
            Vb => VelocitybSol % Values
            VbPerm => VelocitybSol % Perm
    ELSE
            WRITE(Message,'(A)') &
                               'No variable > Velocityb < found'
            CALL FATAL(SolverName,Message)
    END IF  
    c=DIM  ! size of the velocity variable
    IF (VelocitybSol % DOFs.NE.c) then
           WRITE(Message,'(A,I1,A,I1)') &
            'Variable Velocityb has ',VelocitybSol % DOFs,' DOFs, should be',c
            CALL FATAL(SolverName,Message)
    End If
    Vb=0.0_dp


    Cost=0._dp
    Cost_surf=0._dp

    DO t=1,Solver % NumberOfActiveElements
       Element => GetActiveElement(t)
       IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
       n = GetElementNOFNodes()

     NodeIndexes => Element % NodeIndexes

 ! set coords of highest occurring dimension to zero (to get correct path element)
        !-------------------------------------------------------------------------------
        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        IF (DIM == 1) THEN !1D SSA
           ElementNodes % y(1:n) = 0.0_dp
           ElementNodes % z(1:n) = 0.0_dp
        ELSE IF (DIM == 2) THEN !2D SSA
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = 0.0_dp
        ELSE
           WRITE(Message,'(a,i1,a)')&
                'It is not possible to compute SSA problems with DOFs=',&
                DIM, ' . Aborting'
           CALL Fatal( SolverName, Message)
           STOP
        END IF

      BodyForce => GetBodyForce()

      NodeCost=0.0_dp
      NodeCost(1:n) = ListGetReal( BodyForce, 'Adjoint Cost', n, NodeIndexes, GotIt)
          IF (.NOT.GotIt) Then
                  WRITE(Message,'(A)') &
                     'No variable >Adjoint Cost< found in "Body Forces"'
                  CALL FATAL(SolverName,Message)
          END IF 
      NodeCost_der=0.0_dp
          
      NodeCost_der(1,1:n)=ListGetReal( BodyForce, 'Adjoint Cost der 1', n, NodeIndexes, GotIt)
      NodeCost_der(2,1:n)=ListGetReal( BodyForce, 'Adjoint Cost der 2', n, NodeIndexes, GotIt)

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
          Cost_surf=Cost_surf+coeff*s
          NodeCostb(1:n)=NodeCostb(1:n) + s*Basis(1:n)

        End do

        c=DIM ! size of the velocity variable
         Do j=1,n
          Do i=1,DIM
            k=(VbPerm(NodeIndexes(j))-1)*c+i
            Vb(k)=Vb(k)+NodeCostb(j)*NodeCost_der(i,j)
          End Do
         End Do
    End do

   Cost=Cost_surf

    TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )

    IF (Parallel) THEN
           CALL MPI_ALLREDUCE(Cost,Cost_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
          CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
          IF (ASSOCIATED(CostVar)) THEN
                 CostVar % Values(1)=Cost_S
          END IF
         IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) then
                 OPEN (12, FILE=CostFile,POSITION='APPEND')
                 write(12,'(e13.5,2x,e15.8)') TimeVar % Values(1),Cost_S
                 CLOSE(12)
         End if
   ELSE
            CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
            IF (ASSOCIATED(CostVar)) THEN
                    CostVar % Values(1)=Cost
            END IF
                    OPEN (12, FILE=CostFile,POSITION='APPEND')
                       write(12,'(e13.5,2x,e15.8)') TimeVar % Values(1),Cost
                    close(12)
   END IF
   
   Return

 1000  format('#date,time,',a1,'/',a1,'/',a4,',',a2,':',a2,':',a2)
!------------------------------------------------------------------------------
END SUBROUTINE AdjointSSA_CostContSolver
!------------------------------------------------------------------------------
! *****************************************************************************
