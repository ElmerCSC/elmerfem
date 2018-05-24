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
! Compute the integrated gradient of the Cost function for the
! Arthern/Gudmindsson Inverse Problem
!
! Serial/Parallel   and 2D/3D
!
! Need:
! - Dirichlet and Neumann Problems Solutions
! - Name of Beta variable and Grad Variable
! - Power formulation if 'Slip Coef'=10^Beta and optimization on Beta
! - Beta2 formulation if 'Slip Coef'=Beta^2  and optimization on Beta
! - Lambda : the regularization coefficient (Jr=0.5 Lambda (dBeta/dx)^2)
!
! *****************************************************************************
SUBROUTINE DJDBeta_Robin( Model,Solver,dt,TransientSimulation )
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
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,NeumannSolName,DirichletSolName
  CHARACTER(LEN=MAX_NAME_LEN) :: VarSolName,GradSolName
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: PointerToVariable, BetaVariable, VeloSolN,VeloSolD,TimeVar
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Nodes_t) :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  REAL(KIND=dp), POINTER ::  VariableValues(:),VelocityN(:),VelocityD(:),BetaValues(:)
  INTEGER, POINTER :: Permutation(:), VeloNPerm(:),VeloDPerm(:),BetaPerm(:),NodeIndexes(:)
  Logical ::  Firsttime=.true.,Found,stat,UnFoundFatal=.TRUE.
  integer :: i,j,t,n,NMAX,NActiveNodes,DIM
  real(kind=dp),allocatable :: VisitedNode(:),db(:),Basis(:),dBasisdx(:,:)
  real(kind=dp),allocatable :: NodeDJ(:)
  real(kind=dp) :: vn(3),vd(3)
  real(kind=dp) :: u,v,w,SqrtElementMetric
  real(kind=dp) :: Lambda,IPGrad
  logical :: PowerFormulation,Beta2Formulation


  save Firsttime,DIM
  save ElementNodes
  save SolverName
  save NeumannSolName,DirichletSolName,VarSolName,GradSolName
  save Lambda
  save PowerFormulation,Beta2Formulation
  
  save VisitedNode,db,NodeDJ,Basis,dBasisdx

  If (Firsttime) then

      DIM = CoordinateSystemDimension()
      WRITE(SolverName, '(A)') 'DJDBeta_Robin'

      NMAX=Solver % Mesh % NumberOfNodes
      allocate(VisitedNode(NMAX),db(NMAX), NodeDJ(Model %  MaxElementNodes), &
               Basis(Model % MaxElementNodes),  &
               dBasisdx(Model % MaxElementNodes,3))

!!!!!!!!!!! get Solver Variables
      SolverParams => GetSolverParams()

      NeumannSolName =  GetString( SolverParams,'Neumann Solution Name', Found)
          IF(.NOT.Found) THEN        
               CALL WARN(SolverName,'Keyword >Neumann Solution Name< not found in section >Equation<')
               CALL WARN(SolverName,'Taking default value >Flow Solution<')
               WRITE(NeumannSolName,'(A)') 'Flow Solution'
          END IF
      DirichletSolName =  GetString( SolverParams,'Dirichlet Solution Name', Found)
          IF(.NOT.Found) THEN        
               CALL WARN(SolverName,'Keyword >Dirichlet Solution Name< not found in section >Equation<')
               CALL WARN(SolverName,'Taking default value >VeloD<')
               WRITE(DirichletSolName,'(A)') 'VeloD'
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

       Lambda =  GetConstReal( SolverParams,'Lambda', Found)
           IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Lambda< not found  in section >Equation<')
              CALL WARN(SolverName,'Taking default value Lambda=0.0')
              Lambda = 0.0
           End if

       PowerFormulation=GetLogical( SolverParams, 'PowerFormulation', Found)
           IF(.NOT.Found) THEN
                   CALL WARN(SolverName,'Keyword >PowerFormulation< not found  in section >Equation<')
                   CALL WARN(SolverName,'Taking default value >FALSE<')
                   PowerFormulation=.FALSE.
           END IF
       Beta2Formulation=GetLogical( SolverParams, 'Beta2Formulation', Found)
           IF(.NOT.Found) THEN
                   CALL WARN(SolverName,'Keyword >Beta2Formulation< not found  in section >Equation<')
                   CALL WARN(SolverName,'Taking default value >FALSE<')
                   Beta2Formulation=.FALSE.
           END IF
       IF (PowerFormulation.and.Beta2Formulation) then
               WRITE(Message,'(A)') 'Can t be PowerFormulation and Beta2Formulation in the same time'
               CALL FATAL(SolverName,Message)
       End if
  
  !!! End of First visit
    Firsttime=.false.
  Endif

 ! Get variables needed by the Solver

        PointerToVariable => VariableGet( Solver % Mesh % Variables, GradSolName,UnFoundFatal=UnFoundFatal)
        VariableValues => PointerToVariable % Values
        Permutation => PointerToVariable % Perm
        VariableValues=0._dp

        BetaVariable => VariableGet( Solver % Mesh % Variables, VarSolName,UnFoundFatal=UnFoundFatal)
        BetaValues => BetaVariable % Values
        BetaPerm => BetaVariable % Perm

        VeloSolN => VariableGet( Solver % Mesh % Variables, NeumannSolName,UnFoundFatal=UnFoundFatal)
        VelocityN => VeloSolN % Values
        VeloNPerm => VeloSolN % Perm

        VeloSolD => VariableGet( Solver % Mesh % Variables, DirichletSolName,UnFoundFatal=UnFoundFatal)
        VelocityD => VeloSolD % Values
        VeloDPerm => VeloSolD % Perm


    VisitedNode=0.0_dp
    db=0.0_dp

    DO t=1,Solver % NumberOfActiveElements

          Element => GetActiveElement(t)
          CALL GetElementNodes( ElementNodes )
          n = GetElementNOFNodes()
          NodeIndexes => Element % NodeIndexes

          ! Compute Nodal Value of DJDBeta
          Do i=1,n
             VisitedNode(NodeIndexes(i))=VisitedNode(NodeIndexes(i))+1.0_dp

             vn=0.0
             vd=0.0
             vn(1) = VelocityN((DIM+1)*(VeloNPerm(NodeIndexes(i))-1)+1)
             vn(2) = VelocityN((DIM+1)*(VeloNPerm(NodeIndexes(i))-1)+2)
             vd(1) = VelocityD((DIM+1)*(VeloDPerm(NodeIndexes(i))-1)+1)
             vd(2) = VelocityD((DIM+1)*(VeloDPerm(NodeIndexes(i))-1)+2)
             if (DIM.eq.3) then
                 vn(3)=VelocityN((DIM+1)*(VeloNPerm(NodeIndexes(i))-1)+3)
                 vd(3)=VelocityD((DIM+1)*(VeloDPerm(NodeIndexes(i))-1)+3)
             endif

              NodeDJ(i)=(calcNorm(vD)-calcNorm(vN))
              IF (PowerFormulation) then
                 NodeDJ(i)=NodeDJ(i)*(10**(BetaValues(BetaPerm(NodeIndexes(i)))))*log(10.0)
              END IF
              IF (Beta2Formulation) then
                      NodeDJ(i)=NodeDJ(i)*2.0_dp*BetaValues(BetaPerm(NodeIndexes(i)))
              ENDIF
           End do

           ! Compute Integrated Nodal Value of DJDBeta
           IntegStuff = GaussPoints( Element )
           DO j=1,IntegStuff % n
              U = IntegStuff % u(j)
              V = IntegStuff % v(j)
              W = IntegStuff % w(j)
              stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric, &
                             Basis,dBasisdx )
              Do i=1,n
                    IPGrad=NodeDJ(i)*Basis(i)
                    If (Lambda /= 0.0) then
                        IPGrad=IPGrad+Lambda*(SUM(dBasisdx(1:n,1)*BetaValues(BetaPerm(NodeIndexes(1:n))))*dBasisdx(i,1))
                        IF (DIM.eq.3) then
                                IPGrad=IPGrad+Lambda*(SUM(dBasisdx(1:n,2)*BetaValues(BetaPerm(NodeIndexes(1:n))))*dBasisdx(i,2))
                        End if     
                    End if
                    db(NodeIndexes(i)) = db(NodeIndexes(i)) + &
                                   SqrtElementMetric*IntegStuff % s(j)*IPGrad
              End do
            End Do
    End do

   Do t=1,Solver % Mesh % NumberOfNodes
     if (VisitedNode(t).lt.1.0_dp) cycle
     VariableValues(Permutation(t))=db(t) 
   End do

   Return

   CONTAINS

           function calcNorm(v) result(v2)
             implicit none
             real(kind=dp) :: v(3),v2

             v2=v(1)*v(1)+v(2)*v(2)+v(3)*v(3)
           end function calcNorm
!------------------------------------------------------------------------------
END SUBROUTINE DJDBeta_Robin
!------------------------------------------------------------------------------


