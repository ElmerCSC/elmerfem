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
! Adjoint Inverse Problem
!
! Serial/Parallel(not halo?)   and 2D/3D
!
! Need:
! - Navier-stokes and Adjoint Problems Solutions
! - Name of Beta variable and Grad Variable
! - Power formulation if 'Slip Coef'=10^Beta and optimization on Beta
! - Beta2 formulation if 'Slip Coef'=Beta^2  and optimization on Beta
! - Lambda : the regularization coefficient (Jr=0.5 Lambda (dBeta/dx)^2)
!
! If DJDBeta should be set to zero for floating ice shelves the following is 
! to be used (default is not to do this):
! - FreeSlipShelves (logical, default false)
! - mask name (string, default GroundedMask)
! Note that if FreeSlipShelves is true not only is DJDBeta set to zero for 
! floating ice, but also the regularisation term NodalRegb.
!
! *****************************************************************************
SUBROUTINE DJDBeta_Adjoint( Model,Solver,dt,TransientSimulation )
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
  TYPE(ValueList_t), POINTER :: BC

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,NeumannSolName,AdjointSolName
  CHARACTER(LEN=MAX_NAME_LEN) :: VarSolName,GradSolName
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: PointerToVariable, BetaVariable, VeloSolN,VeloSolD
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Nodes_t) :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff

  REAL(KIND=dp), POINTER ::  VariableValues(:),VelocityN(:),VelocityD(:),BetaValues(:)
  INTEGER, POINTER :: Permutation(:), VeloNPerm(:),VeloDPerm(:),BetaPerm(:),NodeIndexes(:)

  real(kind=dp),allocatable :: VisitedNode(:),db(:),Basis(:),dBasisdx(:,:)
  real(kind=dp),allocatable :: nodalbetab(:),NodalRegb(:)
  real(kind=dp) :: betab
  real(kind=dp) :: u,v,w,SqrtElementMetric,s
  real(kind=dp) :: Lambda
  REAL(KIND=dp) :: Normal(3),Tangent(3),Tangent2(3),Vect(3)

  integer :: i,j,k,e,t,n,NMAX,NActiveNodes,DIM
  integer :: p,q

  logical :: PowerFormulation,Beta2Formulation
  Logical ::  Firsttime=.true.,Found,stat,UnFoundFatal=.TRUE.
  Logical :: NormalTangential1,NormalTangential2

  ! Variables for setting DJDBeta to zero for ice shelves.
  TYPE(Variable_t), POINTER   :: PointerToMask => NULL()
  REAL(KIND=dp), POINTER      :: MaskValues(:) => NULL()
  INTEGER, POINTER            :: MaskPerm(:) => NULL()
  LOGICAL                     :: FreeSlipShelves
  CHARACTER(LEN=MAX_NAME_LEN) :: MaskName

  save FreeSlipShelves,MaskName
  save SolverName,NeumannSolName,AdjointSolName,VarSolName,GradSolName
  save VisitedNode,db,Basis,dBasisdx,nodalbetab,NodalRegb
  save Firsttime,DIM,Lambda
  save ElementNodes
  save PowerFormulation,Beta2Formulation

  If (Firsttime) then

     DIM = CoordinateSystemDimension()
     WRITE(SolverName, '(A)') 'DJDBeta_Adjoint'
     
     NMAX=Solver % Mesh % NumberOfNodes
     allocate(VisitedNode(NMAX),db(NMAX),  &
          nodalbetab(Model %  MaxElementNodes),&
          NodalRegb(Model %  MaxElementNodes),&
          Basis(Model % MaxElementNodes),  &
          dBasisdx(Model % MaxElementNodes,3))
     
!!!!!!!!!!! get Solver Variables
     SolverParams => GetSolverParams()
     
     NeumannSolName =  GetString( SolverParams,'Flow Solution Name', Found)
     IF(.NOT.Found) THEN        
        CALL WARN(SolverName,'Keyword >Neumann Solution Name< not found in section >Solver<')
        CALL WARN(SolverName,'Taking default value >Flow Solution<')
        WRITE(NeumannSolName,'(A)') 'Flow Solution'
     END IF
     AdjointSolName =  GetString( SolverParams,'Adjoint Solution Name', Found)
     IF(.NOT.Found) THEN        
        CALL WARN(SolverName,'Keyword >Adjoint Solution Name< not found in section >Solver<')
        CALL WARN(SolverName,'Taking default value >Adjoint<')
        WRITE(AdjointSolName,'(A)') 'Adjoint'
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
     Lambda =  GetConstReal( SolverParams,'Lambda', Found)
     IF(.NOT.Found) THEN
        CALL WARN(SolverName,'Keyword >Lambda< not found  in section  >Solver<')
        CALL WARN(SolverName,'Taking default value Lambda=0.0')
        Lambda = 0.0
     End if
     
     FreeSlipShelves=GetLogical( SolverParams, 'FreeSlipShelves', Found)
     IF(.NOT.Found) THEN
        CALL WARN(SolverName,'Keyword >FreeSlipShelves< not found in solver params')
        CALL WARN(SolverName,'Taking default value >FALSE<')
        FreeSlipShelves=.FALSE.
     END IF
     IF (FreeSlipShelves) THEN
        MaskName =  GetString( SolverParams,'mask name', Found)
        IF(.NOT.Found) THEN
           CALL WARN(SolverName,'Keyword >mask name< not found in solver section')
           CALL WARN(SolverName,'Taking default value >GroundedMask<')
           WRITE(MaskName,'(A)') 'GroundedMask'
        END IF
     END IF

!!! End of First visit
     Firsttime=.false.
  Endif
  
  PointerToVariable => VariableGet( Model % Mesh % Variables, GradSolName,UnFoundFatal=UnFoundFatal)
  VariableValues => PointerToVariable % Values
  Permutation => PointerToVariable % Perm
  VariableValues=0._dp
  
  BetaVariable => VariableGet( Model % Mesh % Variables, VarSolName,UnFoundFatal=UnFoundFatal)
  BetaValues => BetaVariable % Values
  BetaPerm => BetaVariable % Perm
  
  VeloSolN => VariableGet( Model % Mesh % Variables, NeumannSolName,UnFoundFatal=UnFoundFatal)
  VelocityN => VeloSolN % Values
  VeloNPerm => VeloSolN % Perm
  
  VeloSolD => VariableGet( Model % Mesh % Variables, AdjointSolName,UnFoundFatal=UnFoundFatal)
  VelocityD => VeloSolD % Values
  VeloDPerm => VeloSolD % Perm
  
  IF (FreeSlipShelves) THEN
     PointerToMask => VariableGet( Model % Variables, MaskName, UnFoundFatal=.TRUE.)
     MaskValues => PointerToMask % Values
     MaskPerm => PointerToMask % Perm
  END IF
  
  VisitedNode=0.0_dp
  db=0.0_dp
  
  DO e=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(e)
     CALL GetElementNodes( ElementNodes )
     n = GetElementNOFNodes()
     NodeIndexes => Element % NodeIndexes
     
     ! Compute Nodal Value of DJDBeta
     BC => GetBC()
     if (.NOT.ASSOCIATED(BC)) CYCLE
     
     NormalTangential1 = GetLogical( BC, &
          'Normal-Tangential Velocity', Found )
     IF (.NOT.Found) then
        NormalTangential1 = GetLogical( BC, &
             'Normal-Tangential '//trim(NeumannSolName), Found)
     END IF
     NormalTangential2 = GetLogical( BC, &
          'Normal-Tangential '//trim(AdjointSolName), Found)
     IF (NormalTangential1.NEQV.NormalTangential2) then
        WRITE(Message,'(A,I1,A,I1)') &
             'NormalTangential Velocity is : ',NormalTangential1, &
             'But NormalTangential Adjoint is : ',NormalTangential2
        CALL FATAL(SolverName,Message)
     ENDIF
     IF (.NOT.NormalTangential1) then 
        WRITE(Message,'(A)') &
             'ALWAYS USE Normal-Tangential COORDINATES with SlipCoef 2=SlipCoef 3'
        CALL FATAL(SolverName,Message)
     ENDIF
     
     VisitedNode(NodeIndexes(1:n))=VisitedNode(NodeIndexes(1:n))+1.0_dp
     
     ! Compute Integrated Nodal Value of DJDBeta
     nodalbetab=0.0_dp
     NodalRegb=0.0_dp
     
     
     IntegStuff = GaussPoints( Element )
     DO t=1,IntegStuff % n
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)
        stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric, &
             Basis,dBasisdx )
        
        s = SqrtElementMetric * IntegStuff % s(t) 
        
        ! compute gradient from Stokes and adjoint computation
        ! follow the compuation of the stiffMatrix as done in the NS solver
        Normal = NormalVector( Element, ElementNodes, u,v,.TRUE. )
        SELECT CASE( Element % TYPE % DIMENSION )
        CASE(1)
           Tangent(1) =  Normal(2)
           Tangent(2) = -Normal(1)
           Tangent(3) =  0.0_dp
           Tangent2   =  0.0_dp
        CASE(2)
           CALL TangentDirections( Normal, Tangent, Tangent2 )
        END SELECT
        
        
        betab=0.0_dp
        Do p=1,n
           Do q=1,n
              
              Do i=2,dim
                 SELECT CASE(i)
                 CASE(2)
                    Vect = Tangent
                 CASE(3)
                    Vect = Tangent2
                 END SELECT
                 
                 Do j=1,DIM
                    Do k=1,DIM
                       betab = betab + s *  Basis(q) * Basis(p) * Vect(j) * Vect(k) * &
                            (- VelocityN((DIM+1)*(VeloNPerm(NodeIndexes(q))-1)+k) * &
                            VelocityD((DIM+1)*(VeloDPerm(NodeIndexes(p))-1)+j))
                    End Do !on k
                 End Do !on j
              End Do !on i
           End Do !on q
        End Do !on p
        
        nodalbetab(1:n)=nodalbetab(1:n)+betab*Basis(1:n)
        
        If (Lambda /= 0.0) then
           NodalRegb(1:n)=NodalRegb(1:n)+&
                s*Lambda*(SUM(dBasisdx(1:n,1)*BetaValues(BetaPerm(NodeIndexes(1:n))))*dBasisdx(1:n,1))
           IF (DIM.eq.3) then
              NodalRegb(1:n)=NodalRegb(1:n)+&
                   s*Lambda*(SUM(dBasisdx(1:n,2)*BetaValues(BetaPerm(NodeIndexes(1:n))))*dBasisdx(1:n,2))
           End if
        End if
     End DO ! on IPs
     
     IF (PowerFormulation) then
        nodalbetab(1:n)=nodalbetab(1:n)*(10**(BetaValues(BetaPerm(NodeIndexes(1:n)))))*log(10.0)
     ENDIF
     IF (Beta2Formulation) then
        nodalbetab(1:n)=nodalbetab(1:n)*2.0_dp*BetaValues(BetaPerm(NodeIndexes(1:n)))
     END IF

     ! Set regularisation to zero for floating points
     IF (FreeSlipShelves) THEN
        DO t=1,n
           IF ( MaskValues(MaskPerm(NodeIndexes(t))).LT.0.0_dp ) THEN
              NodalRegb(t) = 0.0_dp
           END IF
        END DO
     END IF

     db(NodeIndexes(1:n)) = db(NodeIndexes(1:n)) + nodalbetab(1:n) + NodalRegb(1:n)
  End do ! on elements
  
  Do t=1,Solver % Mesh % NumberOfNodes
     if (VisitedNode(t).lt.1.0_dp) cycle
     VariableValues(Permutation(t)) = db(t) 
     IF (FreeSlipShelves) THEN
        IF ( MaskValues(MaskPerm(t)).LT.0.0_dp ) THEN
           VariableValues(Permutation(t)) = 0.0_dp
        END IF
     END IF
  End do
  
  IF (FreeSlipShelves) THEN
     NULLIFY(PointerToMask)
     NULLIFY(MaskValues)
     NULLIFY(MaskPerm)
  END IF
  
  Return
  
CONTAINS
  
  function calcNorm(v) result(v2)
    implicit none
    real(kind=dp) :: v(3),v2
    
    v2=v(1)*v(1)+v(2)*v(2)+v(3)*v(3)
  end function calcNorm

  function scalar(v1,v2) result(vr)
    implicit none
    real(kind=dp) :: v2(3),v1(3),vr
    
    vr=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
  end function scalar
  
!------------------------------------------------------------------------------
END SUBROUTINE DJDBeta_Adjoint
!------------------------------------------------------------------------------
         

