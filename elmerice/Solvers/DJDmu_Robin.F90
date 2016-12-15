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
! Compute the nodal gradient of the Cost function with respect to the ice viscosity
!  for the Robin inverse method  (Arthern & Gudmundsson, J. Glaciol., 2010)
!
! Serial/Parallel   and 2D/3D
!
!
! .sif parameters:
!    In the Solver section:
!     - Neumann Solution Name = String ; name of the variable for the neumann problem ("Flow Solution" default)
!     - Dirichlet Solution Name = String ; name of the variable for the Dirichlet problem ("VeloD"  default)
!     - Optimized Variable Name = String ;  ("Mu" default)
!     - Gradient Variable Name = String ; t ("DJDMu" default)
!     - SquareFormulation = Logical ; True if the viscosity ios defined as alpha^2
!                                               and optimisation on alpha to insure Mu> 0
!
!
! In the  Material Section:
!     if SquareFormulation = False:
!          Viscosity = Equals mu
!     If SquareFormulation = True:
!         Viscosity = Variable mu
!           Real MATC "tx*tx"
!
!
!  Execute this solver in the main ice body in conjunction with :
!       -CostSolver_Robin.f90: To compute the cost function;
!          !!ATTENTION!! : No regularistaion yet put Lamda=Real 0.0 for the computation of the cost function
!       -Optimise_m1qn3[Serial/Parallel].f90: for the optimization
!
! *****************************************************************************
SUBROUTINE DJDMu_Robin( Model,Solver,dt,TransientSimulation )
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
!!!!  Variables utiles pour les elements et les fonctions de base
  TYPE(Element_t),POINTER ::  Element
  TYPE(Nodes_t) :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  TYPE(ValueList_t), POINTER :: SolverParams, Material
  real(kind=dp),allocatable :: Basis(:),dBasisdx(:,:)
  real(kind=dp) :: u,v,w,SqrtElementMetric
  INTEGER, POINTER :: NodeIndexes(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

!!!!! variables Elmer
   TYPE(Variable_t), POINTER :: Variable, GradVariable,VeloSolN,VeloSolD
   REAL(KIND=dp), POINTER :: Values(:),GradValues(:),VelocityN(:),VelocityD(:)
   INTEGER, POINTER :: Perm(:),GradPerm(:),VeloNPerm(:),VeloDPerm(:)
   CHARACTER(LEN=MAX_NAME_LEN) :: VarSolName,GradSolName,NeumannSolName,DirichletSolName

!! autres variables
  real(kind=dp),allocatable :: VisitedNode(:),db(:)
  real(kind=dp),allocatable :: NodeDJ(:),NodalVeloN(:,:),NodalVeloD(:,:)
  real(kind=dp),allocatable :: m(:),cs(:)
  real(kind=dp) :: vn(3),vd(3),LGradN(3,3),LGradD(3,3),SRD(3,3),SRN(3,3)
  real(kind=dp) :: IPGrad,SecInv

  integer :: i,j,t,n,NMAX,NActiveNodes,DIM

  Logical :: Firsttime=.true.,Found,stat,UnFoundFatal=.TRUE.
  logical :: SquareFormulation


  save Firsttime,DIM
  save ElementNodes
  save SolverName
  save NeumannSolName,DirichletSolName,VarSolName,GradSolName
  save SquareFormulation
  save VisitedNode,db,NodeDJ,Basis,dBasisdx,NodalVeloN,NodalVeloD
  save  m,cs

  If (Firsttime) then

      DIM = CoordinateSystemDimension()
      WRITE(SolverName, '(A)') 'DJDMu_Robin'

      NMAX=Solver % Mesh % NumberOfNodes
      allocate(VisitedNode(NMAX),db(NMAX), NodeDJ(Model %  MaxElementNodes), &
               Basis(Model % MaxElementNodes),  &
               dBasisdx(Model % MaxElementNodes,3), &
               NodalVeloN(3,Model % MaxElementNodes),NodalVeloD(3,Model % MaxElementNodes), &
               m(Model % MaxElementNodes),cs(Model % MaxElementNodes))

!!!!!!!!!!! get Solver Variables
      SolverParams => GetSolverParams()

      NeumannSolName =  GetString( SolverParams,'Neumann Solution Name', Found)
          IF(.NOT.Found) THEN        
               CALL WARN(SolverName,'Keyword >Neumann Solution Name< not found in section >Solver<')
               CALL WARN(SolverName,'Taking default value >Flow Solution<')
               WRITE(NeumannSolName,'(A)') 'Flow Solution'
          END IF
      DirichletSolName =  GetString( SolverParams,'Dirichlet Solution Name', Found)
          IF(.NOT.Found) THEN        
               CALL WARN(SolverName,'Keyword >Dirichlet Solution Name< not found in section >Solver<')
               CALL WARN(SolverName,'Taking default value >VeloD<')
               WRITE(DirichletSolName,'(A)') 'VeloD'
          END IF

      VarSolName =  GetString( SolverParams,'Optimized Variable Name', Found)
             IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Optimized Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >Mu<')
                    WRITE(VarSolName,'(A)') 'Mu'
              END IF
      GradSolName =  GetString( SolverParams,'Gradient Variable Name', Found)
             IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Gradient Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >DJDMu<')
                    WRITE(GradSolName,'(A)') 'DJDmu'
             END IF
       SquareFormulation=GetLogical( SolverParams, 'SquareFormulation', Found)
           IF(.NOT.Found) THEN
                   CALL WARN(SolverName,'Keyword >SquareFormulation< not found  in section >Solver<')
                   CALL WARN(SolverName,'Taking default value >FALSE<')
                   SquareFormulation=.FALSE.
           END IF
  
  !!! End of First visit
    Firsttime=.false.
  Endif

 ! Get variables needed by the Solver

        GradVariable => VariableGet( Solver % Mesh % Variables, GradSolName,UnFoundFatal=UnFoundFatal)
        GradValues => GradVariable % Values
        GradPerm => GradVariable % Perm

        Variable => VariableGet( Solver % Mesh % Variables, VarSolName,UnFoundFatal=UnFoundFatal)
        Values => Variable % Values
        Perm => Variable % Perm

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
          Material => GetMaterial()
          CALL GetElementNodes( ElementNodes )
          n = GetElementNOFNodes()
          NodeIndexes => Element % NodeIndexes

          NodalVeloN = 0.0d0
          NodalVeloD = 0.0d0
          DO i=1, dim
             NodalVeloN(i,1:n) = VelocityN((DIM+1)*(VeloNPerm(NodeIndexes(1:n))-1)+i)
             NodalVeloD(i,1:n) = VelocityD((DIM+1)*(VeloDPerm(NodeIndexes(1:n))-1)+i)
          END DO

          !! exposant nodal
          m=ListGetReal(Material, 'Viscosity Exponent', n, NodeIndexes)
          cs=ListGetReal(Material, 'Critical Shear Rate', n, NodeIndexes)

          ! Compute Nodal Value of DJDmu
          Do i=1,n
             VisitedNode(NodeIndexes(i))=VisitedNode(NodeIndexes(i))+1.0_dp

             u=Element % Type % NodeU(i)
             v=Element % Type % NodeV(i)
             w=Element % Type % NodeW(i)

             stat=ElementInfo(Element,ElementNodes,u,v,w, &
                                SqrtElementMetric,Basis,dBasisdx)

             LGradN=0.0_dp
             LGradD=0.0_dp
             LGradN = MATMUL( NodalVeloN(:,1:n), dBasisdx(1:n,:) )
             SRN = 0.5 * ( LGradN + TRANSPOSE(LGradN) )
             LGradD = MATMUL( NodalVeloD(:,1:n), dBasisdx(1:n,:) )
             SRD = 0.5 * ( LGradD + TRANSPOSE(LGradD) )

              NodeDJ(i)=2._dp*(calcNorm2(SRD)-calcNorm2(SRN))
              if (m(i).ne.1.0_dp) then
                      SecInv=2.0_dp*calcNorm2(SRN) ! le carre du Second invariant de D
                      If (SecInv.lt.cs(i)*cs(i)) SecInv=cs(i)*cs(i)
                      NodeDJ(i)=NodeDJ(i)*SecInv**(0.5_dp*(m(i)-1._dp))
              end if

              if (SquareFormulation) then
                      NodeDJ(i)=NodeDJ(i)*2.0_dp*Values(Perm(NodeIndexes(i)))
              End if
           End do

           ! Compute Integrated Nodal Value of DJDmu
           IntegStuff = GaussPoints( Element )
           DO j=1,IntegStuff % n
              U = IntegStuff % u(j)
              V = IntegStuff % v(j)
              W = IntegStuff % w(j)
              stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric, &
                             Basis,dBasisdx )
              Do i=1,n
                    IPGrad=NodeDJ(i)*Basis(i)

                    db(NodeIndexes(i)) = db(NodeIndexes(i)) + &
                                   SqrtElementMetric*IntegStuff % s(j)*IPGrad
              End do
            End Do
    End do

   Do t=1,Solver % Mesh % NumberOfNodes
     if (VisitedNode(t).lt.1.0_dp) cycle
     GradValues(GradPerm(t))=db(t) 
   End do

   Return

   CONTAINS

           function calcNorm(v) result(v2)
             implicit none
             real(kind=dp) :: v(3),v2

             v2=v(1)*v(1)+v(2)*v(2)+v(3)*v(3)
           end function calcNorm

           function calcNorm2(v) result(v2)
             implicit none
             real(kind=dp) :: v(3,3),v2
             integer :: i,j

             v2=0._dp
             Do i=1,3
               Do j=1,3
                 v2=v2+v(i,j)*v(j,i)
               End do
             End do

           end function calcNorm2
!------------------------------------------------------------------------------
END SUBROUTINE DJDMu_Robin
!------------------------------------------------------------------------------


