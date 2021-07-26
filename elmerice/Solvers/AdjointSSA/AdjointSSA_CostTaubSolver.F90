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
SUBROUTINE AdjointSSA_CostTaubSolver_init0(Model,Solver,dt,TransientSimulation )
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

END SUBROUTINE AdjointSSA_CostTaubSolver_init0
! ******************************************************************************
! *****************************************************************************
SUBROUTINE AdjointSSA_CostTaubSolver( Model,Solver,dt,TransientSimulation )
! *****************************************************************************
!------------------------------------------------------------------------------
!
!  Compute a cost function that penalises 1rst derivative of the basal shear stress
!     J=int_{Pb dimension} 0.5 * ( dTau_b/dx )**2
!
!     The basal friction is computed at the nodes following:
!       Tau_b=beta Velocity_nom**fm
!          with fm the friction exponent
!
!   and provides the derivatives with respect to beta en velocity
!
!   Be Careful, by default, this solver 
!      - reset derivatives wrt beta to 0 (Set Reset DJDBeta = Logical False,
!                   if values have been computed in a previous cost function)
!      - do not reset Cost and derivatives wrt velocity to 0 (Set Reset Cost value = Logical True,
!                    if values have been computed in a previous cost function)
!
!     OUTPUT are : J ; DJDBeta ; Velocityb
!
!     INPUT PARAMETERS are:
!
!      In solver section:
!               Reset Cost Value = Logical (default = .FALSE.)
!               Reset DJDBeta = Logical (default = .True.)
!               Cost Filename = File (default = 'CostTaub.dat')
!               Cost Variable Name = String (default= 'CostValue')
!               DJDBeta Name = String (default= 'DJDBeta')
!               Lambda = Real (default = 1.0)
!
!
!      Variables
!                SSAVelocity (solution of the SSA pb)
!                Velocityb (forcing for the adjoint pb)
!
!     In Material:
!       Keywords related to SSA Friction law (only linear and Weertman)
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
  TYPE(Variable_t), POINTER :: CostVar
  TYPE(Variable_t), POINTER :: TimeVar
  TYPE(Variable_t), POINTER :: DJDVariable !derivée/beta
  TYPE(Variable_t), POINTER :: VVar ! ssavelocity
  TYPE(Variable_t), POINTER :: VbVar ! velocityb
  REAL(KIND=dp),POINTER :: DJDValues(:),VValues(:),VbValues(:)
  INTEGER,POINTER :: DJDPerm(:),VPerm(:),VbPerm(:)
  INTEGER :: DOFs 
  TYPE(Nodes_t),SAVE :: ElementNodes
  INTEGER :: N
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Element_t), POINTER :: Element
  INTEGER, POINTER :: NodeIndexes(:)
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  REAL(KIND=dp) :: U,V,W,SqrtElementMetric
  REAL(KIND=dp),ALLOCATABLE,SAVE :: Basis(:), dBasisdx(:,:)
  REAL(KIND=dp),ALLOCATABLE,SAVE :: Beta(:),Vnode(:,:),Vn(:),Taub(:),NodalBetaDer(:)
  REAL(KIND=dp),ALLOCATABLE,SAVE :: Betab(:),Vnodeb(:,:),Vnb(:),Taubb(:)
  REAL(KIND=dp) :: Cost,Cost_S
  REAL(KIND=dp) :: coeff,coeffb,s,Lambda
  REAL(KIND=dp) :: fm

  CHARACTER(LEN=MAX_NAME_LEN),SAVE :: CostSolName
  CHARACTER(LEN=MAX_NAME_LEN),SAVE :: CostFile
  CHARACTER(LEN=MAX_NAME_LEN) :: DefaultCostFile="CostTaub.dat"
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="CostTaub"
  CHARACTER(LEN=MAX_NAME_LEN) :: SName
  CHARACTER(LEN=MAX_NAME_LEN) :: Friction
  CHARACTER*10 :: date,temps

  LOGICAL :: Reset 
  LOGICAL :: ResetCost
  LOGICAL :: Found
  LOGICAL :: stat
  LOGICAL :: HaveBetaDer
  LOGICAL, SAVE :: Parallel
  LOGICAL, SAVE :: Firsttime=.TRUE.

  INTEGER :: i,t,p,j
  INTEGER :: ierr


  SolverParams => GetSolverParams()

  Lambda =  GetConstReal( SolverParams,'Lambda', Found)
  IF(.NOT.Found) THEN
     CALL WARN(SolverName,'Keyword >Lambda< not found  in section >Solver<')
     CALL WARN(SolverName,'Taking default value Lambda=1.0')
     Lambda = 1.0
  End if

  ResetCost =  GetLogical( SolverParams,'Reset Cost Value', Found)
  IF(.NOT.Found) ResetCost=.FALSE.

!!! SOME INITIALISATION AT FIRST TIME
  If (Firsttime) then
    N = model % MaxElementNodes
    allocate(ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N))
    allocate( Basis(N),dBasisdx(N,3))
    allocate( Beta(N),Vnode(N,3),Vn(N),Taub(N),NodalBetaDer(N))
    allocate( Betab(N),Vnodeb(N,3),Vnb(N),Taubb(N))

!!!!!!! Check for parallel run 
    Parallel = .FALSE.
    IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
            IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
                    Parallel = .TRUE.
            END IF
    END IF

    !! check some kws
    SName =  GetString( SolverParams,'DJDBeta Name', Found)
    IF(.NOT.Found) THEN
          CALL WARN(SolverName,'Keyword >DJDBeta Name< not found in section >Solver<')
          CALL WARN(SolverName,'Taking default value >DJDBeta<')
          WRITE(SName,'(A)') 'DJDBeta'
          CALL ListAddString(  SolverParams, 'DJDBeta Name', TRIM(SName))
    END IF
    !!
    CostSolName =  GetString( SolverParams,'Cost Variable Name', Found)
    IF(.NOT.Found) THEN
       CALL WARN(SolverName,'Keyword >Cost Variable Name< not found  in section >Solver<')
       CALL WARN(SolverName,'Taking default value >CostValue<')
       WRITE(CostSolName,'(A)') 'CostValue'
    END IF


    CostFile = ListGetString(Solver % Values,'Cost Filename',Found )
    IF (.NOT. Found) CostFile = DefaultCostFile
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

    SName =  ListGetString( SolverParams,'DJDBeta Name', UnFoundFatal=.TRUE.)
    DJDVariable => VariableGet( Solver % Mesh % Variables,TRIM(SName),UnFoundFatal=.TRUE.)
    DJDValues => DJDVariable % Values
    DJDPerm => DJDVariable % Perm
    Reset =  ListGetLogical( SolverParams,'Reset DJDBeta', Found)
    if (Reset.OR.(.NOT.Found)) DJDValues=0.0_dp ! reset to zero herz!
    

    VVar => VariableGet( Solver % Mesh % Variables, 'ssavelocity',UnFoundFatal=.TRUE.)
    VValues => VVar % Values
    VPerm => VVar % Perm
    DOFs = VVar % DOFs

    VbVar => VariableGet( Solver % Mesh % Variables, 'velocityb',UnFoundFatal=.TRUE.)
    VbValues => VbVar % Values
    VbPerm => VbVar % Perm
    IF (ResetCost) VbValues = 0._dp
    IF (VbVar%DOFs.NE.DOFs) CALL FATAL(SolverName,'Dimension error')

    Cost=0._dp

    DO t=1,Solver % NumberOfActiveElements
       Element => GetActiveElement(t)
       IF (CheckPassiveElement(Element)) THEN
          CYCLE
       END IF
       n = GetElementNOFNodes(Element)

       NodeIndexes => Element % NodeIndexes

! set coords of highest occurring dimension to zero (to get correct path element)
       !-------------------------------------------------------------------------------
       ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes(1:n))
       IF (DOFs.EQ.2 ) THEN
         ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes(1:n))
       ELSE
         ElementNodes % y(1:n) = 0.0_dp
       ENDIF
       ElementNodes % z(1:n) = 0.0_dp

      ! Get Friction coefficient
       Material => GetMaterial(Element)
       Friction = ListGetString(Material, 'SSA Friction Law', Found,UnFoundFatal=.TRUE.)

       SELECT CASE(Friction)
         CASE('linear')
          fm = 1.0_dp
         CASE('weertman')
          fm = ListGetConstReal( Material, 'SSA Friction Exponent', UnFoundFatal=.TRUE.)
         CASE DEFAULT
           CALL FATAL(SolverName,'Friction should be linear or Weertman')
      END SELECT
      
       Beta(1:n) = ListGetReal( Material, 'SSA Friction Parameter', n, NodeIndexes,UnFoundFatal=.TRUE.)
       NodalBetaDer(1:n) = ListGetReal( Material, 'SSA Friction Parameter Derivative',n, NodeIndexes,Found=HaveBetaDer)
       DO i=1,n
         Vn(i)=0._dp
         Do j=1,DOFs
           Vnode(i,j)=VValues(VVar%DOFs*(VPerm(NodeIndexes(i))-1)+j)
           Vn(i)=Vn(i)+Vnode(i,j)*Vnode(i,j)
         END DO
         IF (Vn(i).GT.AEPS) THEN
           Taub(i)=Beta(i)*Vn(i)**(fm/2)
         ELSE
           Taub(i)=0._dp
         END IF
       END DO

!------------------------------------------------------------------------------
!    Numerical integration
!------------------------------------------------------------------------------
        IntegStuff = GaussPoints( Element )

        Taubb=0._dp
        DO p=1,IntegStuff % n
          U = IntegStuff % u(p)
          V = IntegStuff % v(p)
          W = IntegStuff % w(p)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
          stat = ElementInfo( Element,ElementNodes,U,V,W,SqrtElementMetric, &
              Basis,dBasisdx )

          s =  SqrtElementMetric * IntegStuff % s(p)

          DO j=1,DOFs
           coeff=0._dp
           DO i=1,n
            coeff=coeff+Taub(i)*dbasisdx(i,j)
           END DO
           Cost=Cost+0.5*s*coeff*coeff

           coeffb=s*Lambda*coeff
           DO i=1,n
            Taubb(i)=Taubb(i)+coeffb*dbasisdx(i,j)
           END DO

          END DO

        End do !IP

        DO i=1,n
          IF (Vn(i).GT.AEPS) THEN
            Vnb(i)=Taubb(i)*Beta(i)*0.5*fm*Vn(i)**(fm/2-1._dp)
            Betab(i)=Taubb(i)*Vn(i)**(fm/2)
            IF (HaveBetaDer) Betab(i)=Betab(i)*NodalBetaDer(i)
          ELSE
            Vnb(i)=0._dp
            Betab(i)=0._dp
          END IF

          DO j=1,DOFs
            Vnodeb(i,j)=Vnb(i)*2*Vnode(i,j)
            VbValues(VbVar%DOFs*(VbPerm(NodeIndexes(i))-1)+j)=&
              VbValues(VbVar%DOFs*(VbPerm(NodeIndexes(i))-1)+j)+Vnodeb(i,j)
          END DO

          DJDValues(DJDPerm(NodeIndexes(i)))=&
              DJDValues(DJDPerm(NodeIndexes(i)))+Betab(i)
        END DO

    End do !Elements

    TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )
    IF (Parallel) THEN
         CALL MPI_ALLREDUCE(Cost,Cost_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
         CostVar => VariableGet( Solver % Mesh % Variables, TRIM(CostSolName) ,UnFoundFatal=.TRUE.)
         IF (ResetCost) then
             CostVar % Values(1)=Lambda*Cost_S
         Else
             CostVar % Values(1)=CostVar % Values(1)+Lambda*Cost_S
         Endif
         IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) then
            OPEN (12, FILE=TRIM(CostFile),POSITION='APPEND')
              WRITE(12,'(e13.5,2x,e15.8)') TimeVar % Values(1),Cost_S
            CLOSE(12)
         End if
   ELSE
         CostVar => VariableGet( Solver % Mesh % Variables, TRIM(CostSolName),UnFoundFatal=.TRUE. )
         IF (ResetCost) then
              CostVar % Values(1)=Lambda*Cost
         Else
              CostVar % Values(1)=CostVar % Values(1)+Lambda*Cost
         Endif
         OPEN (12, FILE=TRIM(CostFile),POSITION='APPEND')
           WRITE(12,'(e13.5,2x,e15.8)') TimeVar % Values(1),Cost
         CLOSE(12)
   END IF
   
   RETURN

 1000  format('#date,time,',a2,'/',a2,'/',a4,',',a2,':',a2,':',a2)
 1001  format('#lambda,',e15.8)
!------------------------------------------------------------------------------
END SUBROUTINE AdjointSSA_CostTaubSolver
!------------------------------------------------------------------------------
! *****************************************************************************
