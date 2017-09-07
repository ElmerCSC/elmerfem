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
!Compute the Cost function of the Arhtern/Gudmundsson inverse Problem
!      as Sum_Surface (vn-vd).(sigma_n-sigma_d).n
!      with a regularization as Sum_bedrock 0.5 Lambda (dBeta/dx)^2
!
!   Serial/Parallel    2D/3D
!
! Need : 
!   - Name of the Cost Variable
!   - Solutions of Neumann and Dirchlet problem
!       (Velocities Only, Stresses are computed here)
!   - Lambda and Beta for regularization
!   - define in the sif Name='surface' and Name='bed' in appropriate BC.
!
! *****************************************************************************
SUBROUTINE CostSolver_Robin( Model,Solver,dt,TransientSimulation )
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
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: DefaultCostFile = 'CostOfT.dat'
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,NeumannSolName,DirichletSolName,CostFile
  CHARACTER(LEN=MAX_NAME_LEN) :: BCName
  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName, VarSolName 
  TYPE(Solver_t), POINTER :: ParSolver
  TYPE(Element_t),POINTER ::  Element,Parent
  TYPE(Variable_t), POINTER :: TimeVar,CostVar
  TYPE(ValueList_t), POINTER :: BC,SolverParams
  TYPE(Nodes_t) :: ElementNodes,ParentNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  INTEGER, POINTER :: NodeIndexes(:)

  Logical :: Firsttime=.true.,Found,Parallel,stat
  integer :: i,j,k,l,t,n,np,NMAX,DIM,ierr

  real(kind=dp) :: Cost,Cost_surf,Cost_bed,Cost_S,Cost_surf_S,Cost_bed_S
  real(kind=dp) :: coeff,sTimesN
  real(kind=dp) :: Viscosity,Viscosityn,Viscosityd
  real(kind=dp) :: pressuren,pressured
  real(kind=dp),dimension(3) :: Normal,vn,vd
  real(kind=dp),dimension(3,3) :: LGradn,LGradd,SRn,SRd,Sn,Sd
  real(kind=dp) :: Lambda
  real(kind=dp) :: u,v,w,s,SqrtElementMetric,PSqrtElementMetric,x,y,z
  REAL(KIND=dp),allocatable :: NodalBeta(:),NodalViscosity(:)
  REAL(KIND=dp),allocatable :: Nodalvn(:,:),Nodalvd(:,:),Nodalvelon(:,:),Nodalvelod(:,:)
  REAL(KIND=dp),allocatable :: Basis(:), PBasis(:),dBasisdx(:,:),PdBasisdx(:,:)

  CHARACTER*10 :: date,temps

  save Firsttime,Parallel,CostFile,DIM,ElementNodes,ParentNodes
  save SolverName,NeumannSolName,DirichletSolName,VarSolname,CostSolName
  save Lambda
  save NodalBeta,NodalViscosity,Nodalvn,Nodalvd,Nodalvelon,Nodalvelod
  save Basis,PBasis,dBasisdx,PdBasisdx


  If (Firsttime) then
     DIM = CoordinateSystemDimension()
     WRITE(SolverName, '(A)') 'CostSolver_Robin'

!!!!!!! Check for parallel run 
    Parallel = .FALSE.
      IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
        IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
          Parallel = .TRUE.
        END IF
      END IF

     NMAX= Model % MaxElementNodes
     allocate(NodalBeta(NMAX),NodalViscosity(NMAX), &
              Nodalvn(DIM+1,NMAX),Nodalvd(DIM+1,NMAX), Nodalvelon(3,NMAX),Nodalvelod(3,NMAX), &
              Basis(NMAX),PBasis(NMAX),&
              dBasisdx(NMAX,3),PdBasisdx(NMAX,3))


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


   NeumannSolName =  GetString( SolverParams,'Neumann Solution Name', Found)
         IF(.NOT.Found) THEN
                 CALL WARN(SolverName,'Keyword >Neumann Solution Name< not found  in section >Equation<')
                 CALL WARN(SolverName,'Taking default value >Flow Solution<')
                 WRITE(NeumannSolName,'(A)') 'Flow Solution'
         END IF

   DirichletSolName =  GetString( SolverParams,'Dirichlet Solution Name', Found)
       IF(.NOT.Found) THEN
           CALL WARN(SolverName,'Keyword >Dirichlet Solution Name< not found  in section >Equation<')
           CALL WARN(SolverName,'Taking default value >Flow Solution<')
           WRITE(NeumannSolName,'(A)') 'Flow Solution'
       End if

   Lambda =  GetConstReal( SolverParams,'Lambda', Found)
       IF(.NOT.Found) THEN
           CALL WARN(SolverName,'Keyword >Lambda< not found  in section >Equation<')
           CALL WARN(SolverName,'Taking default value Lambda=0.0')
           Lambda = 0.0
       End if
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


  
  !!! End of First visit
    Firsttime=.false.
  Endif


    Cost=0._dp
    Cost_surf=0.0_dp
    Cost_bed=0.0_dp

    DO t=1,Solver % Mesh % NumberOfBoundaryElements
      Element => GetBoundaryElement(t)

      BC => GetBC()
      IF ( .NOT. ASSOCIATED(BC) ) CYCLE

      BCName =  ListGetString( BC,'Name', Found)
      IF((BCName /= 'surface').AND.(BCName /= 'bed')) CYCLE


      CALL GetElementNodes( ElementNodes )
      n = GetElementNOFNodes()
      NodeIndexes => Element % NodeIndexes


      IF (BCName == 'surface') THEN
         Parent => Element % BoundaryInfo % Left
         IF ( .NOT. ASSOCIATED(Parent) ) &
         Parent => Element % BoundaryInfo % Right

         np = GetElementNOFDOFs(Parent)
         CALL GetElementNodes( ParentNodes, Parent )

         NodalViscosity(1:np) = GetReal( GetMaterial(Parent),'Viscosity',UElement=Parent )
         CALL GetVectorLocalSolution(Nodalvn,NeumannSolName,UElement=Parent )
         CALL GetVectorLocalSolution(Nodalvd,DirichletSolName,UElement=Parent )
         NodalVelon=0._dp
         NodalVelod=0._dp
         Do k=1,DIM
            NodalVelon(k,1:np)=Nodalvn(k,1:np)
            NodalVelod(k,1:np)=Nodalvd(k,1:np)
         End do

       ELSE IF (BCName == 'bed') Then
         CALL GetScalarLocalSolution(NodalBeta,VarSolName,Element)
       END IF

!------------------------------------------------------------------------------
!    Numerical integration
!------------------------------------------------------------------------------
        IntegStuff = GaussPoints( Element )

        DO i=1,IntegStuff % n
          U = IntegStuff % u(i)
          V = IntegStuff % v(i)
          W = IntegStuff % w(i)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
          stat = ElementInfo( Element,ElementNodes,U,V,W,SqrtElementMetric, &
              Basis,dBasisdx )

          s =  SqrtElementMetric * IntegStuff % s(i)
          IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
                  x = SUM( ElementNodes % x(1:n)*Basis(1:n) )
                  y = SUM( ElementNodes % y(1:n)*Basis(1:n) )
                  z = SUM( ElementNodes % z(1:n)*Basis(1:n) )
                  s = s *  CoordinateSqrtMetric(x,y,z)
          END IF

          IF (BCName == 'surface') THEN
            Normal = NormalVector( Element, ElementNodes, U, V, .TRUE. )

            CALL GetParentUVW( Element,n,Parent,np,U,V,W,Basis )
            stat = ElementInfo( Parent,ParentNodes,U,V,W,PSqrtElementMetric, &
                         PBasis,PdBasisdx )

            Viscosity = SUM( NodalViscosity(1:np)*PBasis(1:np) )

            Do k=1,3
               vn(k)=SUM( Nodalvelon(k,1:np)*PBasis(1:np) )
               vd(k)=SUM( Nodalvelod(k,1:np)*PBasis(1:np) )
            End do

            Viscosityn = EffectiveViscosity( Viscosity, 1.0_dp, NodalVelon(1,1:np), Nodalvelon(2,1:np), Nodalvelon(3,1:np), &
                         Parent, ParentNodes, np, np, u, v, w )
            LGradn = MATMUL( NodalVelon(:,1:np), PdBasisdx(1:np,:) )
            SRn = 0.5 * ( LGradn + TRANSPOSE(LGradn) )
            Pressuren = SUM( Nodalvn(DIM+1,1:np)*PBasis(1:np) )
            sn=2.0*Viscosityn*SRn
            Do k=1,3
               sn(k,k)=sn(k,k)-Pressuren
            End do


            Viscosityd = EffectiveViscosity( Viscosity, 1.0_dp, Nodalvelod(1,1:np), Nodalvelod(2,1:np), Nodalvelod(3,1:np), &
                         Parent, ParentNodes, np, np, u, v, w )
            LGradd = MATMUL( NodalVelod(:,1:np), PdBasisdx(1:np,:) )
            SRd = 0.5 * ( LGradd + TRANSPOSE(LGradd) )
            Pressured = SUM( Nodalvd(DIM+1,1:np)*PBasis(1:np) )
            sd=2.0*Viscosityd*SRd
            Do k=1,3
               sd(k,k)=sd(k,k)-Pressured
            End do

            coeff=0._dp
            Do k=1,3
               sTimesN=0.0_dp
               Do l=1,3
                 sTimesN=sTimesN+(sn(k,l)-sd(k,l))*Normal(l)
               End do
               coeff=coeff+(Vn(k)-Vd(k))*sTimesN
            End do

            Cost_surf=Cost_surf+coeff*s

           ELSE IF (BCName == 'bed') THEN
             coeff=SUM(NodalBeta(1:n) * dBasisdx(1:n,1))
             coeff=coeff*coeff
             IF (DIM.eq.3) then
                    coeff=coeff+ & 
                    SUM(NodalBeta(1:n)*dBasisdx(1:n,2))*SUM(NodalBeta(1:n) * dBasisdx(1:n,2))
            END IF
            Cost_bed=Cost_bed+coeff*s

           END IF

        End do
    End do

   Cost=Cost_surf+0.5*Lambda*Cost_bed

   TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )


   IF (Parallel) THEN
           CALL MPI_ALLREDUCE(Cost,Cost_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
           CALL MPI_ALLREDUCE(Cost_surf,Cost_surf_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
           CALL MPI_ALLREDUCE(Cost_bed,Cost_bed_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
          CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
          IF (ASSOCIATED(CostVar)) THEN
                 CostVar % Values(1)=Cost_S
          END IF
         IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) then
                 OPEN (12, FILE=CostFile,POSITION='APPEND')
                 write(12,'(e13.5,2x,e15.8,2x,e15.8,2x,e15.8)') TimeVar % Values(1),Cost_S,Cost_surf_S,Cost_bed_S
                 CLOSE(12)
         End if
   ELSE
              CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
              IF (ASSOCIATED(CostVar)) THEN
                   CostVar % Values(1)=Cost
              END IF
               OPEN (10, FILE=CostFile,POSITION='APPEND')
               write(10,'(e13.5,2x,e15.8,2x,e15.8,2x,e15.8)') TimeVar % Values(1),Cost,Cost_surf,Cost_bed
               close(10)
   END IF

   Return


!------------------------------------------------------------------------------
END SUBROUTINE CostSolver_Robin
!------------------------------------------------------------------------------
! *****************************************************************************
