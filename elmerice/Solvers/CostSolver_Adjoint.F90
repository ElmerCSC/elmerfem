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
!Compute the Cost function of the Adjoint inverse Problem
!      as Integral_Surface Node_Cost ds
!      with a regularization as Sum_bedrock 0.5 Lambda (dBeta/dx)^2
!
!   Serial/Parallel    2D/3D
!
! Need :
!   - Name of the Cost Variable
!   - Lambda and Beta for regularization
!   - define in the sif Name='surface' and Name='bed' in appropriate BC.
!   - define in the sif in the surface BC:
!                'Adjoint Cost = Real ...' : The nodal value of the cost
!                'Adjoint Cost der 1 = Real ...' : The derivative of 'Adjoint Cost' w.r.t. u-velocity
!                'Adjoint Cost der 2 = Real ...' : The derivative of 'Adjoint Cost' w.r.t. v-velocity
!                'Adjoint Cost der 3 = Real ...' : The derivative of 'Adjoint Cost' w.r.t. w-velocity
!
! *****************************************************************************
SUBROUTINE CostSolver_Adjoint_init( Model,Solver,dt,TransientSimulation )
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
  INTEGER :: NormInd, LineInd, i
  LOGICAL :: GotIt, MarkFailed, AvoidFailed
  CHARACTER(LEN=MAX_NAME_LEN) :: Name

  Name = ListGetString( Solver % Values, 'Equation',GotIt)
  IF( .NOT. ListCheckPresent( Solver % Values,'Variable') ) THEN
      CALL ListAddString( Solver % Values,'Variable',&
          '-nooutput -global '//TRIM(Name)//'_var')
 END IF
END
!******************************************************************************
SUBROUTINE CostSolver_Adjoint( Model,Solver,dt,TransientSimulation )
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
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: DefaultCostFile = 'CostOfT.dat'
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,CostFile
  CHARACTER(LEN=MAX_NAME_LEN) :: BCName,CostSolName,VarSolName
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: TimeVar,CostVar,BetaSol
  TYPE(Variable_t), POINTER :: VelocitybSol
  TYPE(ValueList_t), POINTER :: BC,SolverParams
  TYPE(Nodes_t) :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  REAL(KIND=dp), POINTER :: Beta(:)
  REAL(KIND=dp), POINTER :: Vb(:)
  INTEGER, POINTER :: NodeIndexes(:), BetaPerm(:)
  INTEGER, POINTER :: VbPerm(:)
  Logical :: Firsttime=.true.,Found,Parallel,stat,Gotit,UnFoundFatal=.TRUE.
  integer :: i,j,k,l,t,n,NMAX,DIM,ierr,c
  real(kind=dp) :: Cost,Cost_bed,Cost_surf,Cost_S,Cost_bed_S,Cost_surf_S,Lambda,Change
  real(kind=dp),Save :: Oldf=0._dp
  real(kind=dp) :: Bu,Bv,u,v,w,s,coeff,SqrtElementMetric,x
  REAL(KIND=dp) :: NodeCost(Model % MaxElementNodes),Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)
  REAL(KIND=dp) :: NodeCostb(Model % MaxElementNodes),NodeCost_der(3,Model %MaxElementNodes)
  CHARACTER*10 :: date,temps

  save Firsttime,Parallel 
  save SolverName,CostSolName,VarSolName,Lambda,CostFile
  save ElementNodes

  DIM = CoordinateSystemDimension()

  If (Firsttime) then

    WRITE(SolverName, '(A)') 'CostSolver_Adjoint'

!!!!!!! Check for parallel run 
    Parallel = (ParEnv % PEs > 1)

!!!!!!!!!!! get Solver Variables
  SolverParams => GetSolverParams()

  CostFile = ListGetString(Solver % Values,'Cost Filename',Found )
    IF (.NOT. Found) CostFile = DefaultCostFile
    CALL DATE_AND_TIME(date,temps)
    If (Parallel) then
        if (ParEnv % MyPe.EQ.0) then
           OPEN (12, FILE=CostFile)
                    write(12,'(a1,a2,a1,a2,a1,a4,5x,a2,a1,a2,a1,a2)') '#',date(5:6),'/',date(7:8),'/',date(1:4), &
                                 temps(1:2),':',temps(3:4),':',temps(5:6)
           CLOSE(12)
         End if
    Else
           OPEN (12, FILE=CostFile)
                    write(12,'(a1,a2,a1,a2,a1,a4,5x,a2,a1,a2,a1,a2)') '#',date(5:6),'/',date(7:8),'/',date(1:4), &
                                 temps(1:2),':',temps(3:4),':',temps(5:6)
           CLOSE(12)
    End if
    

  VarSolName =  GetString( SolverParams,'Optimized Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Optimized Variable Name< not found in section >Solver<')
              CALL WARN(SolverName,'Taking default value >Beta<')
              WRITE(VarSolName,'(A)') 'Beta'
      END IF

   CostSolName =  GetString( SolverParams,'Cost Variable Name', Found)
          IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Cost Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >CostValue<')
                    WRITE(CostSolName,'(A)') 'CostValue'
          END IF

   Lambda =  GetConstReal( SolverParams,'Lambda', Found)
   IF(.NOT.Found) THEN
           CALL WARN(SolverName,'Keyword >Lambda< not found  in section >Solver<')
           CALL WARN(SolverName,'Taking default value Lambda=0.0')
           Lambda = 0.0
   End if
  
  !!! End of First visit
    Firsttime=.false.
  Endif

    BetaSol => VariableGet( Model % Mesh % Variables, VarSolName,UnFoundFatal=UnFoundFatal)
    Beta => BetaSol % Values
    BetaPerm => BetaSol % Perm

    VelocitybSol => VariableGet( Model % Mesh % Variables, 'Velocityb',UnFoundFatal=UnFoundFatal)
    Vb => VelocitybSol % Values
    VbPerm => VelocitybSol % Perm
    c=DIM + 1 ! size of the velocity variable
    IF (VelocitybSol % DOFs.NE.c) then
           WRITE(Message,'(A,I1,A,I1)') &
            'Variable Velocityb has ',VelocitybSol % DOFs,' DOFs, should be',c
            CALL FATAL(SolverName,Message)
    End If
    Vb=0.0_dp


    Cost=0._dp
    Cost_surf=0._dp
    Cost_bed=0._dp
    DO t=1,Model % Mesh % NumberOfBoundaryElements

      Element => GetBoundaryElement(t)

      BC => GetBC()
      IF ( .NOT. ASSOCIATED(BC) ) CYCLE

      BCName =  ListGetString( BC,'Name', Found)
      IF((BCName /= 'surface').AND.(BCName /= 'bed')) CYCLE

      CALL GetElementNodes( ElementNodes )
      n = GetElementNOFNodes()
      NodeIndexes => Element % NodeIndexes

      NodeCost=0.0_dp
      IF (BCName == 'surface') THEN
          NodeCost(1:n) = ListGetReal(BC, 'Adjoint Cost', n, NodeIndexes, GotIt,&
               UnFoundFatal=UnFoundFatal)
          NodeCost_der=0.0_dp
          
          NodeCost_der(1,1:n)=ListGetReal(BC, 'Adjoint Cost der 1', n, NodeIndexes, GotIt)
          NodeCost_der(2,1:n)=ListGetReal(BC, 'Adjoint Cost der 2', n, NodeIndexes, GotIt)
          NodeCost_der(3,1:n)=ListGetReal(BC, 'Adjoint Cost der 3', n, NodeIndexes, GotIt)
          
      Else IF (BCName == 'bed') Then
          NodeCost(1:n)=Beta(BetaPerm(NodeIndexes(1:n)))
      End if

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
          
           IF (BCName == 'surface') Then   
            coeff = SUM(NodeCost(1:n)  * Basis(1:n))
            Cost_surf=Cost_surf+coeff*s
            NodeCostb(1:n)=NodeCostb(1:n) + s*Basis(1:n)
           else IF (BCName == 'bed') Then
            coeff = SUM(NodeCost(1:n) * dBasisdx(1:n,1))
            coeff =  coeff * coeff
            IF (DIM.eq.3) then
                    coeff=coeff+ & 
                    SUM(NodeCost(1:n)*dBasisdx(1:n,2))*SUM(NodeCost(1:n) * dBasisdx(1:n,2))
            END IF
            Cost_bed=Cost_bed+coeff*s
           else 
            coeff = 0.0
           End if

        End do
        IF (BCName == 'surface') Then
        c=DIM + 1 ! size of the velocity variable
         Do j=1,n
          Do i=1,DIM
            k=(VbPerm(NodeIndexes(j))-1)*c+i
            Vb(k)=Vb(k)+NodeCostb(j)*NodeCost_der(i,j)
          End Do
         End Do
        END if
    End do

   Cost=Cost_surf+0.5*Lambda*Cost_bed

    TimeVar => VariableGet( Model % Mesh % Variables, 'Time' )

    IF (Parallel) THEN
           CALL MPI_ALLREDUCE(Cost,Cost_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
           CALL MPI_ALLREDUCE(Cost_surf,Cost_surf_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
           CALL MPI_ALLREDUCE(Cost_bed,Cost_bed_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
          CostVar => VariableGet( Model % Mesh % Variables, CostSolName )
          IF (ASSOCIATED(CostVar)) THEN
                 CostVar % Values(1)=Cost_S
          END IF
         IF (ParEnv % MyPE == 0) then
                 OPEN (12, FILE=CostFile,POSITION='APPEND')
                 write(12,'(e13.5,2x,e15.8,2x,e15.8,2x,e15.8)') TimeVar % Values(1),Cost_S,Cost_surf_S,Cost_bed_S
                 CLOSE(12)
         End if
   ELSE
            CostVar => VariableGet( Model % Mesh % Variables, CostSolName )
            IF (ASSOCIATED(CostVar)) THEN
                    CostVar % Values(1)=Cost
            END IF
            OPEN (12, FILE=CostFile,POSITION='APPEND')
              write(12,'(e13.5,2x,e15.8,2x,e15.8,2x,e15.8)') TimeVar % Values(1),Cost,Cost_surf,Cost_bed
            close(12)
            Cost_S=Cost
   END IF
   
   Solver % Variable % Values(1)=Cost_S
   Solver % Variable % Norm = Cost_S
   IF (SIZE(Solver%Variable % Values) == Solver%Variable % DOFs) THEN
      !! MIMIC COMPUTE CHANGE STYLE
      Change=2.*(Cost_S-Oldf)/(Cost_S+Oldf)
      Change=abs(Change)
      WRITE( Message, '(a,g15.8,g15.8,a)') &
              'SS (ITER=1) (NRM,RELC): (',Cost_S, Change,&
              ' ) :: Cost'
      CALL Info( 'ComputeChange', Message, Level=3 )
      Oldf=Cost_S
   ENDIF

   Return
!------------------------------------------------------------------------------
END SUBROUTINE CostSolver_Adjoint
!------------------------------------------------------------------------------
! *****************************************************************************
