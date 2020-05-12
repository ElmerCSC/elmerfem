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
! *  Original Date: April 2020; Adapted from AdjointSSA_CostRegSolver
! * 
! *****************************************************************************
SUBROUTINE Adjoint_CostRegSolver_init0(Model,Solver,dt,TransientSimulation )
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
END SUBROUTINE Adjoint_CostRegSolver_init0
! *****************************************************************************
SUBROUTINE Adjoint_CostRegSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Compute a regularisation term and update the
!  gradient of the cost function with respect to the regularized variable.
!
!   Regularisation by default is: int_{Pb dimension} 0.5 * (d(var)/dx)**2 
!   A priori regularisation can also be used ( A priori Regularisation=True) :
!                                 int_{Pb dimension} 0.5 *(1/sigma**2)*(var-var{a_priori})**2
!
!     OUTPUT are : J and DJDvar
!                      
!
!    !!!!! BE careful it will reset Cost and DJ to 0 by default !!!!
!      !!! If other cost and gradient are computed before (i.e. from the adjoint pb, 
!       use "<Reset Cost Value> = False" to add cost and gradient to previously computed values !!!
!
!
!       Required Sif parameters are:
!
!          In the solver section:
!               Cost Filename=File (default: CostOfT.dat),
!               Optimized Variable Name= String (default='Beta'),
!               Gradient Variable Name= String (default = 'DJDBeta'),
!               Cost Variable Name= String (default='Cost Value'),
!               Lambda= Real (default 1.0),
!               Reset Cost Value= Logical (default = True),
!               A priori Regularisation= Logical (default = .False.),
!
!          In Body Force section:
!               <VarSolName> a priori value = Real (default =0.0),
!               <VarSolName> variance = real (default=1.0)
!
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
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='Adjoint_CostRegSolver'
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: DefaultCostFile = 'CostOfT.dat'
  CHARACTER(LEN=MAX_NAME_LEN) :: CostFile
  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName,VarSolName,GradSolName,varname

  TYPE(Variable_t), POINTER :: TimeVar,CostVar

  TYPE(Variable_t), POINTER :: Variable,DJDVariable
  REAL(KIND=dp), POINTER :: Values(:),DJDValues(:)
  INTEGER, POINTER :: Perm(:),DJDPerm(:)

  TYPE(ValueList_t), POINTER :: SolverParams,BodyForce

  TYPE(Element_t),POINTER ::  Element
  TYPE(Nodes_t),SAVE :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  REAL(KIND=dp),ALLOCATABLE,SAVE :: Basis(:), dBasisdx(:,:)
  INTEGER, POINTER :: NodeIndexes(:)

  LOGICAL, SAVE :: Firsttime=.true.,Parallel
  LOGICAL :: Found,stat,Gotit
  LOGICAL :: BoundarySolver
  integer :: i,j,k,l,t,n,NMAX,DIM,ierr,c

  real(kind=dp) :: Cost,Cost_S,Lambda
  real(kind=dp) :: Area,Area_S
  real(kind=dp) :: u,v,w,s,coeff_reg,SqrtElementMetric,x

  REAL(KIND=dp),dimension(:),allocatable,SAVE :: NodeAp,NodeRMS,NodalRegb
  REAL(KIND=dp),dimension(:),allocatable,SAVE :: NodeValues,NodalDer,NodalGrad
  REAL(KIND=dp) :: IPerr,IPvar

  LOGICAL :: Apriori,Reset
  LOGICAL :: HaveNodalVariable
  LOGICAL :: HaveDer

  CHARACTER*10 :: date,temps


  SolverParams => GetSolverParams()

  !! check if we are on a boundary or in the bulk
  BoundarySolver = ( Solver % ActiveElements(1) > Model % Mesh % NumberOfBulkElements )
  IF (BoundarySolver) THEN
     DIM = CoordinateSystemDimension() - 1
  ELSE
     DIM = CoordinateSystemDimension()
  ENDIF

! get some needed solver parameters
!! Cost File for Output
   CostFile = ListGetString(Solver % Values,'Cost Filename',Found )
   IF (.NOT. Found) CostFile = DefaultCostFile

!! Name of the variable to regularise
  VarSolName =  ListGetString( SolverParams,'Optimized Variable Name', Found=HaveNodalVariable)

!! Name of the variable to regularise
   GradSolName =  ListGetString( SolverParams,'Gradient Variable Name', UnFoundFatal=.TRUE.)
!! Name of the variable with the cost function
   CostSolName =  ListGetString( SolverParams,'Cost Variable Name', Found )
   IF(.NOT.Found) THEN
       CALL WARN(SolverName,'Keyword >CostSolName< not found  in section >Solver<')
       CALL WARN(SolverName,'Taking default value CostValue')
       CostSolName='CostValue'
   End if

!! Optional weighting term
   Lambda =  GetConstReal( SolverParams,'Lambda', Found)
   IF(.NOT.Found) THEN
       CALL WARN(SolverName,'Keyword >Lambda< not found  in section >Solver<')
       CALL WARN(SolverName,'Taking default value Lambda=1.0')
       Lambda = 1.0
   End if

!! Do we need to reset cost and DJDVar to 0? Default YES
   Reset =  GetLogical( SolverParams,'Reset Cost Value', Found)
   IF(.NOT.Found) Reset=.True.

!! What type of regularistaion ? Default penalise 1st derivative
   Apriori =  GetLogical( SolverParams,'A priori Regularisation', Found)
            IF(.NOT.Found) Apriori=.False.

!!! SOME INITIALISATION AT FIRST TIME
  If (Firsttime) then
    N = model % MaxElementNodes
    allocate(ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N))
    allocate(Basis(N),dBasisdx(N,3))
    allocate(NodeAp(N),NodeRMS(N),NodeValues(N),NodalRegb(N),NodalDer(N),NodalGrad(N))

!!!!!!! Check for parallel run 
    Parallel = .FALSE.
    IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
            IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
                    Parallel = .TRUE.
            END IF
    END IF

!!!!!!!!!!!  initiate Cost File

    CALL DATE_AND_TIME(date,temps)
    If (Parallel) then
        if (ParEnv % MyPe.EQ.0) then
           OPEN (12, FILE=CostFile)
                   write(12,1000) date(5:6),date(7:8),date(1:4),temps(1:2),temps(3:4),temps(5:6)
                   write(12,1001) Lambda
                   write(12,'(A)') '# iter, Jreg'
           CLOSE(12)
         End if
    Else
           OPEN (12, FILE=CostFile)
                   write(12,1000) date(5:6),date(7:8),date(1:4),temps(1:2),temps(3:4),temps(5:6)
                   write(12,1001) Lambda
                   write(12,'(A)') '# iter, Jreg'
           CLOSE(12)
    End if
  
  !!! End of First visit
    Firsttime=.false.
  Endif
!!!! INITIALISATION DONE
    IF (HaveNodalVariable) THEN
      Variable => VariableGet( Solver % Mesh % Variables, VarSolName  , UnFoundFatal=.TRUE.)
      Values => Variable % Values
      Perm => Variable % Perm
    END IF

    DJDVariable => VariableGet( Solver % Mesh % Variables, GradSolName , UnFoundFatal=.TRUE. )
    DJDValues => DJDVariable % Values
    DJDPerm => DJDVariable % Perm
    IF (Reset) DJDValues=0.0_dp

    Cost=0._dp
    Area=0._dp

    DO t=1,Solver % NumberOfActiveElements
       Element => GetActiveElement(t)
       IF (CheckPassiveElement(Element)) THEN
          CYCLE
       END IF
       BodyForce => GetBodyForce(Element)
       IF (.NOT.ASSOCIATED(BodyForce)) THEN
          IF (Apriori.OR.(.NOT.HaveNodalVariable)) &
             CALL FATAL(SolverName,'Body force should be associated for this regularisation')
       ENDIF

       n = GetElementNOFNodes()

       NodeIndexes => Element % NodeIndexes

 ! set coords of highest occurring dimension to zero (to get correct path element)
        !-------------------------------------------------------------------------------
        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        SELECT CASE (DIM)
         CASE (1) 
           ElementNodes % y(1:n) = 0.0_dp
           ElementNodes % z(1:n) = 0.0_dp
         CASE (2)
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = 0.0_dp
         CASE (3)
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
         END SELECT

 ! Compute integrated cost
      IF (Apriori) then
         NodeAp(1:n) = ListGetReal( BodyForce, 'CostReg Nodal Prior', n, NodeIndexes, GotIt)
          IF (.NOT.GotIt) Then
            CALL WARN(SolverName,'No value for the prior found in <Body Force>, default is 0')
            NodeAp(1:n) = 0._dp
          END IF 
          NodeRMS(1:n)=ListGetReal( BodyForce,'CostReg Nodal std', n, NodeIndexes, GotIt)
          IF (.NOT.GotIt) Then
            CALL WARN(SolverName,'No value for the standard deviation found in <Body Force>, default is 1')
            NodeRMS=1.0_dp
          END IF 
     END IF

 ! Nodal values of the variable        
     IF (HaveNodalVariable) THEN
       NodeValues(1:n)=Values(Perm(NodeIndexes(1:n)))
       HaveDer=.FALSE.
     ELSE
       NodeValues(1:n)=ListGetReal( BodyForce,'CostReg Nodal Variable',n, NodeIndexes, UnFoundFatal=.TRUE.)
       NodalDer(1:n) = ListGetReal( BodyForce,'CostReg Nodal Variable derivative',n,NodeIndexes,Found=HaveDer)
     END IF

!------------------------------------------------------------------------------
!    Numerical integration
!------------------------------------------------------------------------------
        NodalRegb = 0.0_dp

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

          x = SUM( ElementNodes % x(1:n) * Basis(1:n) )
          s = 1.0d0

          IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
             s = 2.0d0 * PI * x 
          END IF
          s = s * SqrtElementMetric * IntegStuff % s(i)
          


          IF (Apriori) then
             IPerr = SUM((NodeValues(1:n)-NodeAp(1:n))*Basis(1:n))
             IPvar = SUM(NodeRMS(1:n)*Basis(1:n))
             coeff_reg=IPerr/IPvar
             coeff_reg =  0.5*coeff_reg*coeff_reg 

             IF (HaveDer) THEN
              NodalGrad(1:n)=Basis(1:n)*NodalDer(1:n)
             ELSE
              NodalGrad(1:n)=Basis(1:n)
             ENDIF

             !Now compute the derivative
               NodalRegb(1:n)=NodalRegb(1:n)+&
                    s*Lambda*IPerr*NodalGrad(1:n)/(IPVar**2.0)
          Else
             coeff_reg=0._dp
             DO k=1,DIM
               coeff_reg = coeff_reg + 0.5*SUM(NodeValues(1:n) * dBasisdx(1:n,k))**2
             END DO
       
             !Now compute the derivative
             DO k=1,DIM

               IF (HaveDer) THEN
                NodalGrad(1:n)=dBasisdx(1:n,k)*NodalDer(1:n)
               ELSE
                NodalGrad(1:n)=dBasisdx(1:n,k)
               ENDIF
               
               NodalRegb(1:n)=NodalRegb(1:n)+&
                    s*Lambda*(SUM(dBasisdx(1:n,k)*NodeValues(1:n))*NodalGrad(1:n))
             END DO
          Endif

          Cost=Cost+coeff_reg*s
          Area=Area+s

        End do !IP

        DJDValues(DJDPerm(NodeIndexes(1:n)))=DJDValues(DJDPerm(NodeIndexes(1:n))) + NodalRegb(1:n)

    End do !Elements

    TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )

    IF (Parallel) THEN
           CALL MPI_ALLREDUCE(Cost,Cost_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
           CALL MPI_ALLREDUCE(Area,Area_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
          CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
          IF (ASSOCIATED(CostVar)) THEN
               IF (Reset) then
                 CostVar % Values(1)=Lambda*Cost_S
               Else
                 CostVar % Values(1)=CostVar % Values(1)+Lambda*Cost_S
               Endif
          END IF
         IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) then
                 OPEN (12, FILE=CostFile,POSITION='APPEND')
                 write(12,'(3(ES20.11E3))') TimeVar % Values(1),Cost_S,sqrt(2*Cost_S/Area_S)
                 CLOSE(12)
         End if
   ELSE
         CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
         IF (ASSOCIATED(CostVar)) THEN
            IF (Reset) then
                CostVar % Values(1)=Lambda*Cost
            Else
                CostVar % Values(1)=CostVar % Values(1)+Lambda*Cost
            Endif
         END IF
         OPEN (12, FILE=CostFile,POSITION='APPEND')
           write(12,'(3(ES20.11E3))') TimeVar % Values(1),Cost,sqrt(2*Cost/Area)
         close(12)
         Cost_S=Cost
   END IF
 
   Solver % Variable % Values(1)=Cost_S

   Return

 1000  format('#date,time,',a2,'/',a2,'/',a4,',',a2,':',a2,':',a2)
 1001  format('#lambda,',e15.8)
!------------------------------------------------------------------------------
END SUBROUTINE Adjoint_CostRegSolver
!------------------------------------------------------------------------------
! *****************************************************************************
