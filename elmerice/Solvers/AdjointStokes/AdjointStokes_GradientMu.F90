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
!  Compute the nodal gradient of the Cost function with respect to the ice
!  viscosity for the control inverse method
!    (cf Morlighem, M., Ice sheet properties inferred by combining numerical
!    modeling and remote sensing data)
!
! Serial/Parallel   and 2D/3D
!
!
! .sif parameters:
!    In the Solver section:
!     - Flow Solution Name = String ; name of the variable for the direct problem ("Flow Solution" default)
!     - Adjoint Solution Name = String ; name of the variable for the Adjoint system ("Adjoint"  default)
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
!       -CostSolver_Adjoint.f90: To compute the cost function;
!          !!ATTENTION!! : No regularistaion yet put Lamda=Real 0.0 for the computation of the cost function
!       -Optimise_m1qn3[Serial/Parallel].f90: for the optimization
!
! *****************************************************************************
SUBROUTINE DJDMu_Adjoint( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
  USE DefUtils
  USE MaterialModels
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
  TYPE(ValueList_t), POINTER :: SolverParams,Material
  real(kind=dp),allocatable :: Basis(:),dBasisdx(:,:)
  real(kind=dp) :: u,v,w,SqrtElementMetric
  INTEGER, POINTER :: NodeIndexes(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

!!!!! variables Elmer
  TYPE(Variable_t), POINTER :: GradVariable, Variable, VeloSolN,VeloSolD
  REAL(KIND=dp), POINTER ::  GradValues(:),VelocityN(:),VelocityD(:),Values(:)
  INTEGER, POINTER :: GradPerm(:), VeloNPerm(:),VeloDPerm(:),Perm(:)
  CHARACTER(LEN=MAX_NAME_LEN) ::GradSolName,NeumannSolName,DirichletSolName,VarSolName

!! autres variables
  real(kind=dp),allocatable :: VisitedNode(:),db(:)
  real(kind=dp),allocatable,dimension(:) :: Ux,Uy,Uz
  real(kind=dp) :: Velo(3),dVelodx(3,3)
  real(kind=dp) :: s,ss,c2,c3
  real(kind=dp) :: mub,Viscosityb
  real(kind=dp),allocatable,dimension(:) :: c2n,c3n
  real(kind=dp),allocatable,dimension(:) :: NodalViscosityb


  integer :: i,j,t,n,NMAX,NpN,NActiveNodes,DIM,e,p,q

  CHARACTER(LEN=MAX_NAME_LEN) :: ViscosityFlag

  logical :: SquareFormulation
  Logical ::  Firsttime=.true.,Found,stat,gotit,UnFoundFatal=.TRUE.


  save Firsttime,DIM
  save ElementNodes
  save SolverName
  save NeumannSolName,DirichletSolName,VarSolName,GradSolName
  save SquareFormulation
  save VisitedNode,db,Basis,dBasisdx
  save Ux,Uy,Uz
  save c2n,c3n
  save NodalViscosityb

  !!!! Firsttime Do some allocation and initialisation
  If (Firsttime) then

      DIM = CoordinateSystemDimension()
      WRITE(SolverName, '(A)') 'DJDMu_Adjoint'

      NMAX=Solver % Mesh % NumberOfNodes
      NpN=Model % MaxElementNodes

      allocate(VisitedNode(NMAX),db(NMAX), &
               Basis(NpN),  &
               dBasisdx(NpN,3), &
               Ux(NpN),Uy(NpN),Uz(NpN),&
               c2n(NpN),c3n(NpN),&
               NodalViscosityb(NpN))

!!!!!!!!!!! get Solver Variables
      SolverParams => GetSolverParams()

      NeumannSolName =  GetString( SolverParams,'Flow Solution Name', Found)
          IF(.NOT.Found) THEN        
               CALL WARN(SolverName,'Keyword >Flow Solution Name< not found in section >Solver<')
               CALL WARN(SolverName,'Taking default value >Flow Solution<')
               WRITE(NeumannSolName,'(A)') 'Flow Solution'
          END IF
      DirichletSolName =  GetString( SolverParams,'Adjoint Solution Name', Found)
          IF(.NOT.Found) THEN        
               CALL WARN(SolverName,'Keyword >Adjoint Solution Name< not found in section >Solver<')
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
                   CALL WARN(SolverName,'Logical Keyword >SquareFormulation< not found  in section >Solver<')
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
        GradValues=0._dp

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

    DO e=1,Solver % NumberOfActiveElements

          Element => GetActiveElement(e)
          Material => GetMaterial()
          CALL GetElementNodes( ElementNodes )
          n = GetElementNOFNodes()
          NodeIndexes => Element % NodeIndexes

          VisitedNode(NodeIndexes(1:n))=VisitedNode(NodeIndexes(1:n))+1.0_dp

          Ux=0.0_dp
          Uy=0.0_dp
          Uz=0.0_dp
          Ux(1:n)=VelocityN((DIM+1)*(VeloNPerm(NodeIndexes(1:n))-1)+1)
          Uy(1:n)=VelocityN((DIM+1)*(VeloNPerm(NodeIndexes(1:n))-1)+2)
          If (DIM.eq.3) Uz(1:n)=VelocityN((DIM+1)*(VeloNPerm(NodeIndexes(1:n))-1)+3)

          !!!!
          nodalViscosityb=0.0_dp

          IntegStuff = GaussPoints( Element )

          DO t=1,IntegStuff%n

             u = IntegStuff % u(t)
             v = IntegStuff % v(t)
             w =IntegStuff % w(t)

             stat = ElementInfo( Element, ElementNodes, u, v, w,SqrtElementMetric, &
                           Basis, dBasisdx) !removed bubbles 

             s = SqrtElementMetric * IntegStuff % s(t)

             mub=0.0_dp
             Do p=1,n
                Do q=1,n
                  Do i=1,DIM
                    Do j=1,DIM
                       mub=mub+ s * dBasisdx(q,j) * dBasisdx(p,j) * &
                         (- VelocityN((DIM+1)*(VeloNPerm(NodeIndexes(q))-1)+i) * &
                          VelocityD((DIM+1)*(VeloDPerm(NodeIndexes(p))-1)+i))

                       mub=mub+ s * dBasisdx(q,i) * dBasisdx(p,j) * &
                         (- VelocityN((DIM+1)*(VeloNPerm(NodeIndexes(q))-1)+j) * &
                           & VelocityD((DIM+1)*(VeloDPerm(NodeIndexes(p))-1)+i))
                    End Do !j
                  End Do !i
                 End Do !q
              End Do !p

              ViscosityFlag = ListGetString( Material,'Viscosity Model', GotIt,UnFoundFatal)

              SELECT CASE( ViscosityFlag )
                CASE('power law')
                DO j=1,3
                   dVelodx(1,j) = SUM( Ux(1:n)*dBasisdx(1:n,j) )
                   dVelodx(2,j) = SUM( Uy(1:n)*dBasisdx(1:n,j) )
                   dVelodx(3,j) = SUM( Uz(1:n)*dBasisdx(1:n,j) )
                END DO

                Velo(1) = SUM( Basis(1:n) * Ux(1:n) )
                Velo(2) = SUM( Basis(1:n) * Uy(1:n) )
                Velo(3) = SUM( Basis(1:n) * Uz(1:n) )

                ss = SecondInvariant(Velo,dVelodx)/2
        
                c2n = ListGetReal( Material, 'Viscosity Exponent', n, NodeIndexes )
                c2 = SUM( Basis(1:n) * c2n(1:n) )

                s = ss

                c3n = ListGetReal( Material, 'Critical Shear Rate',n, NodeIndexes ,gotIt )
                IF (GotIt) THEN
                  c3 = SUM( Basis(1:n) * c3n(1:n) )
                  IF(s < c3**2) THEN
                     s = c3**2
                  END IF
                END IF

                Viscosityb=mub*s**((c2-1)/2)

                CASE default
                    CALL FATAL(SolverName,'Viscosity Model has to be power Law')
              END SELECT 

              nodalViscosityb(1:n)=nodalViscosityb(1:n)+Viscosityb*Basis(1:n)
          End Do !on IPs

          IF (SquareFormulation) then
               nodalViscosityb(1:n)=nodalViscosityb(1:n)*2.0_dp*Values(Perm(NodeIndexes(1:n)))
          END IF

          db(NodeIndexes(1:n)) = db(NodeIndexes(1:n)) + nodalViscosityb(1:n)
       End Do ! on elements

   Do t=1,Solver % Mesh % NumberOfNodes
     if (VisitedNode(t).lt.1.0_dp) cycle
     GradValues(GradPerm(t))=db(t) 
   End do

   Return

!------------------------------------------------------------------------------
END SUBROUTINE DJDMu_Adjoint
!------------------------------------------------------------------------------


