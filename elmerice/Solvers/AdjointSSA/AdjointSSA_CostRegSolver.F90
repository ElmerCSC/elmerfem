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
! *****************************************************************************
SUBROUTINE AdjointSSA_CostRegSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Compute a regularisation term for SSA inverse problems and update the
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
!               Problem Dimension=Integer (default:coordinate system dimension),
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
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: DefaultCostFile = 'CostOfT.dat'
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,CostFile
  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName,VarSolName,GradSolName,varname
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: TimeVar,CostVar

  TYPE(Variable_t), POINTER :: Variable,DJDVariable
  REAL(KIND=dp), POINTER :: Values(:),DJDValues(:)
  INTEGER, POINTER :: Perm(:),DJDPerm(:)

  TYPE(ValueList_t), POINTER :: SolverParams,BodyForce
  TYPE(Nodes_t) :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  INTEGER, POINTER :: NodeIndexes(:)

  Logical :: Firsttime=.true.,Found,Parallel,stat,Gotit
  integer :: i,j,k,l,t,n,NMAX,DIM,ierr,c

  real(kind=dp) :: Cost,Cost_S,Lambda
  real(kind=dp) :: u,v,w,s,coeff_reg,SqrtElementMetric,x

  REAL(KIND=dp),dimension(:),allocatable,SAVE :: NodeAp,NodeRMS,NodeValues,NodalRegb
  REAL(KIND=dp) :: Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)
  REAL(KIND=dp) :: IPerr,IPvar

  LOGICAL :: Apriori,Reset

  CHARACTER*10 :: date,temps

  save Firsttime,Parallel 
  save SolverName,CostSolName,VarSolName,Lambda,CostFile
  save ElementNodes

   WRITE(SolverName, '(A)') 'CostSolver_Regular'

  SolverParams => GetSolverParams()

!! Dimension of the pb; ie with SSA we can be 1D or 2D on a 2D mesh, or 2D on a 3D mesh
  DIM=GetInteger(SolverParams ,'Problem Dimension',Found)
  If (.NOT.Found) then
     CALL WARN(SolverName,'Keyword >Problem Dimension< not found, assume DIM = CoordinateSystemDimension()')
     DIM = CoordinateSystemDimension()
  Endif

! get some needed solver parameters
!! Cost File for Output
  CostFile = ListGetString(Solver % Values,'Cost Filename',Found )
    IF (.NOT. Found) CostFile = DefaultCostFile

!! Name of the variable to regularise
  VarSolName =  GetString( SolverParams,'Optimized Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Optimized Variable Name< not found in section >Solver<')
              CALL WARN(SolverName,'Taking default value >Beta<')
              WRITE(VarSolName,'(A)') 'Beta'
      END IF

!! Name of the variable to regularise
  GradSolName =  GetString( SolverParams,'Gradient Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Optimized Variable Name< not found in section >Solver<')
              CALL WARN(SolverName,'Taking default value >DJDBeta<')
              WRITE(GradSolName,'(A)') 'DJDBeta'
      END IF


!! Name of the variable with the cost function
   CostSolName =  GetString( SolverParams,'Cost Variable Name', Found)
          IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Cost Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >CostValue<')
                    WRITE(CostSolName,'(A)') 'CostValue'
          END IF

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
    allocate(NodeAp(N),NodeRMS(N),NodeValues(N),NodalRegb(N))

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

    Variable => VariableGet( Solver % Mesh % Variables, VarSolName  )
    IF ( ASSOCIATED( Variable ) ) THEN
            Values => Variable % Values
            Perm => Variable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable >',VarSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF  
    DJDVariable => VariableGet( Solver % Mesh % Variables, GradSolName  )
    IF ( ASSOCIATED( DJDVariable ) ) THEN
            DJDValues => DJDVariable % Values
            DJDPerm => DJDVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable >',VarSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF  
    IF (Reset) DJDValues=0.0_dp


    Cost=0._dp

    DO t=1,Solver % NumberOfActiveElements
       Element => GetActiveElement(t)
       IF (CheckPassiveElement(Element)) THEN
          !PRINT *,ParEnv%myPe,'REG: PASSIVE ELEEMNT'
          CYCLE
       END IF
       IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
       n = GetElementNOFNodes()

       NodeIndexes => Element % NodeIndexes

 ! set coords of highest occurring dimension to zero (to get correct path element)
        !-------------------------------------------------------------------------------
        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        IF (DIM == 1) THEN !1D SSA
           ElementNodes % y(1:n) = 0.0_dp
           ElementNodes % z(1:n) = 0.0_dp
        ELSE IF (DIM == 2) THEN !2D SSA
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = 0.0_dp
        ELSE
           WRITE(Message,'(a,i1,a)')&
                'It is not possible to compute SSA problems with DOFs=',&
                DIM, ' . Aborting'
           CALL Fatal( SolverName, Message)
           STOP
        END IF

 ! Compute inetgrated cost


      IF (Apriori) then
         BodyForce => GetBodyForce()
         write(varname,'(A,A)') trim(VarSolName),' a priori value'
         NodeAp(1:n) = 0._dp
         NodeAp(1:n) = ListGetReal( BodyForce, trim(varname), n, NodeIndexes, GotIt)
          IF (.NOT.GotIt) Then
                  WRITE(Message,'(A,A,A)') &
                     'No variable >',trim(varname),'< found in "Body Forces" default is 0'
                  CALL Info(SolverName,Message,level=6)
          END IF 
          write(varname,'(A,A)') trim(VarSolName),' variance'
          NodeRMS(1:n)=ListGetReal( BodyForce, trim(varname), n, NodeIndexes, GotIt)
          IF (.NOT.GotIt) Then
                  WRITE(Message,'(A,A,A)') &
                     'No variable >',trim(varname),'< found in "Body Forces" default is 1'
                  CALL Info(SolverName,Message,level=6)
                  NodeRMS=1.0_dp
          END IF 
     END IF

 ! Nodal values of the variable        
      NodeValues(1:n)=Values(Perm(NodeIndexes(1:n)))

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
             coeff_reg =  coeff_reg*coeff_reg 

             !Now compute the derivative
               NodalRegb(1:n)=NodalRegb(1:n)+&
                    s*Lambda*IPerr*Basis(1:n)/(IPVar**2.0)
          Else
             coeff_reg = SUM(NodeValues(1:n) * dBasisdx(1:n,1))
             coeff_reg =  coeff_reg*coeff_reg 
             IF (DIM.eq.2) then
                  coeff_reg=coeff_reg+ & 
                  SUM(NodeValues(1:n)*dBasisdx(1:n,2))*SUM(NodeValues(1:n) * dBasisdx(1:n,2))
             END IF
       
             !Now compute the derivative
               NodalRegb(1:n)=NodalRegb(1:n)+&
                    s*Lambda*(SUM(dBasisdx(1:n,1)*NodeValues(1:n))*dBasisdx(1:n,1))
               IF (DIM.eq.2) then
                  NodalRegb(1:n)=NodalRegb(1:n)+&
                          s*Lambda*(SUM(dBasisdx(1:n,2)*NodeValues(1:n))*dBasisdx(1:n,2))
               End if
          Endif


          Cost=Cost+0.5*coeff_reg*s

        End do !IP

        DJDValues(DJDPerm(NodeIndexes(1:n)))=DJDValues(DJDPerm(NodeIndexes(1:n))) + NodalRegb(1:n)

    End do !Elements

    TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )

    IF (Parallel) THEN
           CALL MPI_ALLREDUCE(Cost,Cost_S,1,&
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
                 write(12,'(e13.5,2x,e15.8)') TimeVar % Values(1),Cost_S
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
                       write(12,'(e13.5,2x,e15.8)') TimeVar % Values(1),Cost
                    close(12)
   END IF
   
   Return

 1000  format('#date,time,',a1,'/',a1,'/',a4,',',a2,':',a2,':',a2)
 1001  format('#lambda,',e15.8)
!------------------------------------------------------------------------------
END SUBROUTINE AdjointSSA_CostRegSolver
!------------------------------------------------------------------------------
! *****************************************************************************
