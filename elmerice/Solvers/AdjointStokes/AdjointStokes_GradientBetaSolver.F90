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
! Adjoint of the Stokes problem for the slip coefficient
! Provide the derivative with respect to the slip coef.
!
! Serial/Parallel(not halo?)   and 2D/3D
!
! Need:
! - Navier-stokes and Adjoint Problems Solutions
! - Grad Variable
!
! *****************************************************************************
SUBROUTINE AdjointStokes_GradientBetaSolver_init0(Model,Solver,dt,TransientSimulation )
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

END SUBROUTINE AdjointStokes_GradientBetaSolver_init0
! *****************************************************************************
SUBROUTINE AdjointStokes_GradientBetaSolver( Model,Solver,dt,TransientSimulation )
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
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(ValueList_t), POINTER :: BC

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="AdjointStokes_GradientBeta"

  INTEGER,SAVE :: DIM

  CHARACTER(LEN=MAX_NAME_LEN),SAVE :: NeumannSolName,AdjointSolName
  TYPE(Variable_t), POINTER :: VeloSolN,VeloSolD
  REAL(KIND=dp), POINTER ::  VelocityN(:),VelocityD(:)
  INTEGER, POINTER :: VeloNPerm(:),VeloDPerm(:)

  CHARACTER(LEN=MAX_NAME_LEN),SAVE :: GradSolName
  TYPE(Variable_t), POINTER :: PointerToVariable
  REAL(KIND=dp), POINTER ::  VariableValues(:)
  INTEGER, POINTER :: Permutation(:)

  TYPE(Element_t),POINTER ::  Element
  TYPE(Nodes_t),SAVE :: ElementNodes
  INTEGER, POINTER :: NodeIndexes(:)
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  REAL(kind=dp) :: u,v,w,SqrtElementMetric,s
  REAL(kind=dp),allocatable,SAVE :: Basis(:),dBasisdx(:,:)
  INTEGER :: n

  LOGICAL :: NormalTangential1,NormalTangential2
  REAL(KIND=dp) :: Normal(3),Tangent(3),Tangent2(3),Vect(3)
  REAL(kind=dp) :: betab
  REAL(kind=dp),allocatable,SAVE :: NodalDer(:),NodalGrad(:)
  LOGICAL :: HaveDer

  INTEGER :: NMAX

  integer :: i,j,k,p,q,e,t

  LOGICAL ::  Found,stat
  LOGICAL :: Reset
  LOGICAL, SAVE :: Firsttime=.true.
  LOGICAL,PARAMETER :: UnFoundFatal=.TRUE.



  SolverParams => GetSolverParams()

  If (Firsttime) then

     DIM = CoordinateSystemDimension()
     
     NMAX=Model % MaxElementNodes
     allocate(Basis(NMAX),dBasisdx(NMAX,3))
     allocate(NodalDer(NMAX),NodalGrad(NMAX))
!!!!!!!!!!! get Solver Variables
     
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

     GradSolName =  GetString( SolverParams,'Gradient Variable Name', Found)
     IF(.NOT.Found) THEN
        CALL WARN(SolverName,'Keyword >Gradient Variable Name< not found  in section >Solver<')
        CALL WARN(SolverName,'Taking default value >DJDB<')
        WRITE(GradSolName,'(A)') 'DJDB'
     END IF

!!! End of First visit
     Firsttime=.false.
  Endif
  
  PointerToVariable => VariableGet( Model % Mesh % Variables, GradSolName,UnFoundFatal=UnFoundFatal)
  VariableValues => PointerToVariable % Values
  Permutation => PointerToVariable % Perm
  Reset =  ListGetLogical( SolverParams,'Reset Gradient Variable', Found)
  if (Reset.OR.(.NOT.Found)) VariableValues=0._dp
  
  VeloSolN => VariableGet( Model % Mesh % Variables, NeumannSolName,UnFoundFatal=UnFoundFatal)
  VelocityN => VeloSolN % Values
  VeloNPerm => VeloSolN % Perm
  
  VeloSolD => VariableGet( Model % Mesh % Variables, AdjointSolName,UnFoundFatal=UnFoundFatal)
  VelocityD => VeloSolD % Values
  VeloDPerm => VeloSolD % Perm
  
  
  
  DO e=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(e)
     CALL GetElementNodes( ElementNodes )
     n = GetElementNOFNodes(Element)
     NodeIndexes => Element % NodeIndexes
     
     ! Compute Nodal Value of DJDBeta
     BC => GetBC(Element)
     if (.NOT.ASSOCIATED(BC)) &
        CALL FATAL(SolverName,'This solver is intended to be executed on a BC')
     
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

     NodalDer(1:n) = ListGetReal(BC,'Slip Coefficient derivative',n,NodeIndexes,Found=HaveDer)
     
     ! Compute Integrated Nodal Value of DJDBeta
     IntegStuff = GaussPoints( Element )
     DO t=1,IntegStuff % n
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)
        stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric, &
             Basis,dBasisdx )
        
        s = SqrtElementMetric * IntegStuff % s(t) 
        
        ! compute gradient from Stokes and adjoint computation
        ! follow the computation of the stiffMatrix as done in the NS solver
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
        
        IF (HaveDer) THEN
          NodalGrad(1:n)=Basis(1:n)*NodalDer(1:n)
        ELSE
          NodalGrad(1:n)=Basis(1:n)
        ENDIF
        VariableValues(Permutation(NodeIndexes(1:n))) = VariableValues(Permutation(NodeIndexes(1:n))) &
                 + betab*NodalGrad(1:n)
        
     End DO ! on IPs
     
  End do ! on elements
  
  Return
  
!------------------------------------------------------------------------------
END SUBROUTINE AdjointStokes_GradientBetaSolver
!------------------------------------------------------------------------------
         

