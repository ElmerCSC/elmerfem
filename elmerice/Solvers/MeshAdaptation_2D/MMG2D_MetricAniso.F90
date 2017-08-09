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
! *  Author: F. Gillet-Chaulet (IGE)
! *  Email:  fabien.gillet-chaulet@univ-grenoble-alpes.fr
! *  Web:    http://elmerice.elmerfem.org
! *
! *  Original Date: 13-07-2017, 
! *
! *****************************************************************************
!> Module to compute the 2D Metric (or only 2D Hessian) for
!>  anisotropic remeshing with MMG2D
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE MMG2D_MetricAniso( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver
     LOGICAL ::  TransientSimulation
     REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='MMG_MetricAniso'
    TYPE(Matrix_t),POINTER :: StiffMatrix
    REAL(KIND=dp), POINTER :: ForceVector(:)

    TYPE(ValueList_t),POINTER :: SolverParams, BodyForce

    TYPE(Variable_t), POINTER :: Variable
    TYPE(Variable_t), POINTER :: GradVariable,TensorVariable,HVariable
    REAL(KIND=dp),DIMENSION(:),POINTER :: Values
    REAL(KIND=dp),DIMENSION(:),POINTER :: GradValues,TensorValues,HValues
    INTEGER,DIMENSION(:),POINTER :: Perm
    INTEGER,DIMENSION(:),POINTER :: GradPerm,TensorPerm,HPerm
    CHARACTER(LEN=MAX_NAME_LEN) :: GradName,TensorName,HName

    TYPE(Element_t),POINTER :: CurrentElement
    INTEGER, POINTER ::  NodeIndexes(:)
    TYPE(Nodes_t),SAVE :: ElementNodes

    REAL(KIND=dp), ALLOCATABLE, SAVE:: LocalMassMatrix(:,:), &
       LocalStiffMatrix(:,:), LocalForce(:), &
       LocalGrad(:,:)

    REAL(KIND=dp) :: Diffusivity,kappa
    REAL(KIND=dp) :: VNorm
    REAL(KIND=dp),Dimension(3) :: Param,Hessian,Metric

    INTEGER :: COMP,t,i
    INTEGER :: n
    INTEGER :: k,node,bf_id
    INTEGER :: istat

    LOGICAL :: HessianOnly=.FALSE.,HessianInitialized=.FALSE.

    LOGICAL :: UnFoundFatal=.TRUE.,GotIt
    LOGICAL,SAVE :: AllocationsDone = .FALSE.
    LOGICAL :: stat

!------------------------------------------------------------------------------
     SolverParams => GetSolverParams(Solver)

     HName = ListGetString( SolverParams, 'Hessian Variable Name',  UnFoundFatal=UnFoundFatal )    
     HVariable => VariableGet( Solver % Mesh % Variables, HName, UnFoundFatal=UnFoundFatal ) 
     HPerm => HVariable % Perm
     HValues => HVariable % Values
     IF (HVariable % DOFs /= 3) &
          CALL Fatal( SolverName, 'Bad dimension of Hessian Variable; should be 3' )

     HessianInitialized=ListGetLogical(SolverParams,'Hessian Initialized',Gotit)
    !--------------------------------------------------------------------
    ! Hessian not provided solve the diffusive equation for each component
    !--------------------------------------------------------------------
     IF (.NOT.HessianInitialized) THEN

       IF (.NOT.ASSOCIATED(Solver % Matrix )) &
          CALL FATAL(SolverName,'Matrix Not Associated')
       StiffMatrix => Solver % Matrix
       ForceVector => StiffMatrix % RHS

       Variable => Solver % Variable
       Perm => Variable % Perm
       Values => Variable % Values
       IF ( Variable % DOFs /= 1 ) &
         CALL Fatal( SolverName, 'Variable DOFs must be equal to 1' )

       GradName = ListGetString( SolverParams, 'Gradient Name', UnFoundFatal=UnFoundFatal )    
       GradVariable => VariableGet( Solver % Mesh % Variables, GradName, UnFoundFatal=UnFoundFatal )
       GradPerm    => GradVariable % Perm
       GradValues  => GradVariable % Values
       IF ( GradVariable % DOFs /=2 ) &
         CALL Fatal( SolverName, 'Gradient Variable DOFs must be equal to 2' )

       Diffusivity = ListGetCReal(SolverParams,'Diffusivity',UnFoundFatal=UnFoundFatal )
!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
       IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN

         IF(CoordinateSystemDimension() /= 2) &
            CALL Fatal( SolverName,'Sorry only 2D Solver')
         IF ( CurrentCoordinateSystem() /= Cartesian ) &
            CALL Fatal( SolverName,'Sorry only For Cartesian coordinate system')

         N = Model % MaxElementNodes
         IF ( AllocationsDone ) THEN
           DEALLOCATE( ElementNodes % x,     &
                       ElementNodes % y,     &
                       ElementNodes % z,     &
                       LocalGrad,            &                      
                       LocalMassMatrix,      &
                       LocalStiffMatrix,     &
                       LocalForce )
         END IF
         ALLOCATE( ElementNodes % x( N ), &
                   ElementNodes % y( N ), &
                   ElementNodes % z( N ), &
                   LocalGrad( 2,N ), &                                    
                   LocalMassMatrix( 2*N,2*N ),  &
                   LocalStiffMatrix( 2*N,2*N ),  &
                   LocalForce( 2*N ),  STAT=istat )

         IF ( istat /= 0 ) &
            CALL Fatal( SolverName , 'Memory allocation error.' )

         AllocationsDone = .TRUE.
        END IF
!------------------------------------------------------------------------------

! Loop over the Tensor components [Exx,Eyy,Exy] 
        DO COMP = 1, 3

          CALL DefaultInitialize()

          DO t=1,Solver % NumberOFActiveElements

            CurrentElement => GetActiveElement(t)
            n = GetElementNOFNodes(CurrentElement)
            NodeIndexes => CurrentElement % NodeIndexes

            ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes(1:n))
            ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes(1:n))
            ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes(1:n))

            LocalGrad = 0.0d0
            DO i=1, 2
               LocalGrad(i,1:n) = GradValues(2*(GradPerm(NodeIndexes(1:n))-1) + i)
            END DO

           kappa=Diffusivity*ElementArea(Solver%Mesh,CurrentElement,n)

           CALL LocalMatrix(COMP, LocalMassMatrix, LocalStiffMatrix, &
                LocalForce, LocalGrad,kappa, CurrentElement, n, ElementNodes )
        !------------------------------------------------------------------------------
        !        Update global matrices from local matrices 
        !------------------------------------------------------------------------------
           CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )
         END DO

         CALL DefaultFinishAssembly()

      !------------------------------------------------------------------------------
      !     Dirichlet boundary conditions
      !------------------------------------------------------------------------------
         CALL DefaultDirichletBCs()

      !------------------------------------------------------------------------------
      !     Solve the system 
      !------------------------------------------------------------------------------
         VNorm = DefaultSolve()

        DO node=1,Solver % Mesh % NumberOfNodes
          k=Perm( node )
          IF( k == 0 ) CYCLE
          HValues(3*(HPerm(node)-1) + COMP) = Values(k)  
        END DO

      END DO ! End DO Comp

      HessianOnly = ListGetLogical(SolverParams,'Compute Hessian Only', Gotit)
      IF (HessianOnly) RETURN
  !------------------------------------------------------------------------------
  !  Hessian is initialized
  !------------------------------------------------------------------------------
      END IF

  !------------------------------------------------------------------------------
  ! Now Compute Metric
  !------------------------------------------------------------------------------
    TensorName = ListGetString( SolverParams, 'Metric Variable Name',  UnFoundFatal=UnFoundFatal )    
    TensorVariable => VariableGet( Solver % Mesh % Variables, TensorName, UnFoundFatal=UnFoundFatal ) 
    TensorPerm => TensorVariable % Perm
    TensorValues => TensorVariable % Values
    IF (TensorVariable % DOFs /= 3) &
          CALL Fatal( SolverName, 'Bad dimension of Tensor Variable; should be 3' )

    DO bf_id=1,CurrentModel % NumberOfBodyForces
       BodyForce => CurrentModel % BodyForces(bf_id) % Values
       IF( ListCheckPresent(BodyForce,TRIM(TensorName) // ' err')&
      .AND.ListCheckPresent(BodyForce,TRIM(TensorName) // ' hmin')&
      .AND.ListCheckPresent(BodyForce,TRIM(TensorName) // ' hmax')) EXIT
    END DO
    if (bf_id.GT.CurrentModel % NumberOfBodyForces) &
       CALL FATAL(SolverName,'Metric keywords not found in any Body Force')

    DO node=1,Solver % Mesh % NumberOfNodes
       Param(1)=ListGetRealAtNode(BodyForce,TRIM(TensorName) // ' hmin', node)
       Param(2)=ListGetRealAtNode(BodyForce,TRIM(TensorName) // ' hmax', node)
       Param(3)= ListGetRealAtNode(BodyForce,TRIM(TensorName) // ' err', node)
       Hessian(1:3)=HValues(3*(HPerm(node)-1) + 1:3)
       CALL ComputeMetric(Hessian,Param,Metric)
       TensorValues(3*(TensorPerm(node)-1) + 1:3)=Metric
    END DO
      
CONTAINS
!------------------------------------------------------------------------------
      SUBROUTINE LocalMatrix(COMP, MassMatrix, StiffMatrix, ForceVector, &
              NodalGrad,kappa, Element, n, Nodes )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: StiffMatrix(:,:), MassMatrix(:,:)
     REAL(KIND=dp) :: NodalGrad(:,:), ForceVector(:)
     REAL(KIND=dp) :: kappa
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
     INTEGER :: n, COMP
!------------------------------------------------------------------------------
!
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ, V_Integ, W_Integ, S_Integ
     REAL(KIND=dp) :: s, u, v, w
     REAL(KIND=dp) :: Basis(2*n), ddBasisddx(1,1,1)
     REAL(KIND=dp) :: dBasisdx(2*n,3), detJ
     REAL(KIND=dp) :: LGrad(2,2), SR(2,2),  Eij(3)
     INTEGER :: N_Integ
     INTEGER :: DIM
     INTEGER :: t,p, q
     LOGICAL :: stat

!------------------------------------------------------------------------------
      DIM = CoordinateSystemDimension()

      ForceVector = 0.0D0
      StiffMatrix = 0.0D0
      MassMatrix  = 0.0D0

      IntegStuff = GaussPoints( Element )

      U_Integ => IntegStuff % u
      V_Integ => IntegStuff % v
      W_Integ => IntegStuff % w
      S_Integ => IntegStuff % s
      N_Integ =  IntegStuff % n
!
!   Now we start integrating
!
      DO t=1,N_Integ

        u = U_Integ(t)
        v = V_Integ(t)
        w = W_Integ(t)

!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
        stat = ElementInfo(Element,Nodes,u,v,w,detJ, &
                 Basis,dBasisdx,ddBasisddx,.FALSE.,.FALSE.)
        s = detJ * S_Integ(t)
!
! Hij= 1/2 (Grad(GradVar) + T^Grad(GradVar))
        SR = 0.0
        Eij = 0.0
        LGrad = MATMUL( NodalGrad(1:DIM,1:n), dBasisdx(1:n,1:DIM) )
        SR = 0.5 * ( LGrad + TRANSPOSE(LGrad) )

        Eij(1) = SR(1,1)
        Eij(2) = SR(2,2)        
        Eij(3) = SR(1,2)        

        DO p=1,n         
          DO q=1,n        
            StiffMatrix(p,q) =  &
               StiffMatrix(p,q) + s*Basis(q)*Basis(p) + &
                s*kappa*SUM(dBasisdx(p,1:DIM)*dBasisdx(q,1:DIM))
          END DO
          ForceVector(p) =  &
                     ForceVector(p) + s*Eij(COMP)*Basis(p) 
        END DO
      END DO

!------------------------------------------------------------------------------
      END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------
      SUBROUTINE ComputeMetric(Hessian,Param,Metric)
         implicit none
         REAL(KIND=dp),dimension(3),intent(in) :: Hessian !(1,1);(2,2);(1,2)
         REAL(KIND=dp),dimension(3),intent(in) :: Param !hmin,hmax,err
         REAL(KIND=dp),dimension(3),intent(out) :: Metric !(1,1);(2,2);(1,2)

         REAL(KIND=dp),parameter :: c=2._dp/9._dp
         REAL(KIND=dp) :: e1,e2,en
         REAL(KIND=dp) :: Delta,l1,l2
         REAL(KIND=dp) :: hmin,hmax,err
         REAL(KIND=dp) :: Lambd1,Lambd2

         hmin=Param(1)
         hmax=Param(2)
         IF (hmax.LT.hmin) CALL FATAL('SolverName','Error Metric Hmax<Hmin')
         err=Param(3)

! compute eigenvalues of Hessian
         IF (abs(Hessian(3)).lt.TINY(1.0_dp)) then
                 l1=Hessian(1)
                 l2=Hessian(2)
                 e1=1.0
                 e2=0.0
         ELSE
           Delta=(Hessian(1)-Hessian(2))*(Hessian(1)-Hessian(2))
           Delta=Delta+4.0*Hessian(3)*Hessian(3)
           if (Delta.lt.0) then
                 CALL FATAL(SolverName,'Metric Delta<0')
           endif
           l1=(Hessian(1)+Hessian(2))+sqrt(Delta)
           l1=l1/2.0
           l2=(Hessian(1)+Hessian(2))-sqrt(Delta)
           l2=l2/2.0

! compute normilized eignevector 1 of Hessian
           e1=Hessian(3)
           e2=(l1-Hessian(1))
           en=sqrt(e1*e1+e2*e2)
           e1=e1/en
           e2=e2/en

         END IF

! compute Metric in eigenframe        
         Lambd1=min(max(c*abs(l1)/err,1.0/(hmax*hmax)),1.0/(hmin*hmin))
         Lambd2=min(max(c*abs(l2)/err,1.0/(hmax*hmax)),1.0/(hmin*hmin))

! Metric in reference frame
         IF (abs(Lambd1-Lambd2).LT.AEPS) THEN
            Metric(1)=(Lambd1+Lambd2)/2
            Metric(2)=(Lambd1+Lambd2)/2
            Metric(3)=0._dp
         ELSE
            Metric(1)=Lambd1*e1*e1+Lambd2*e2*e2  !Metric(1,1)
            Metric(2)=Lambd1*e2*e2+Lambd2*e1*e1  !Metric(2,2)
            Metric(3)=Lambd1*e1*e2-Lambd2*e1*e2  !Metric(1,2)=Metric(2,1)
         ENDIF

         Return
!------------------------------------------------------------------------------
        End !Metric
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      END SUBROUTINE MMG2D_MetricAniso
!------------------------------------------------------------------------------
