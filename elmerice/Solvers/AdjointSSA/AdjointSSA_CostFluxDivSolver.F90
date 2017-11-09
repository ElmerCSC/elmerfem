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
!
! *****************************************************************************
SUBROUTINE AdjointSSA_CostFluxDivSolver( Model,Solver,dt,TransientSimulation )
! *****************************************************************************
!------------------------------------------------------------------------------
!  
!  Compute a cost function that measure the ice flux divergence anomaly
!                 and the required forcing for the SSA adjoint problem
!     J=int_{Pb dimension} 0.5 * (dhdt{obs} + u.grad(H) + H* div(u) - MB )**2
!
!     OUTPUT are : J ; DJDu (==Velocityb variable used as forcing of the SSA adjoint problem)
!                      DJDZb (optional); DJDZs(optional)
!
!     TODO : add a variance term to regularise the cost
!
!    !!!!! BE careful it will reset Cost , Velocityb, and DJZb; DJDZs to 0 by default !!!!
!      !!! If other cost and gradient are computed before , 
!       use "<Reset Cost Value> = False" to add cost and gradient to previously computed values !!!
!
!     INPUT PARAMETERS are:
!
!      In solver section:
!               Problem Dimension = Integer (default = Coordinate system dimension)
!               Compute DJDZb = Logical (default= .True.)
!               Compute DJDZs = Logical (default= .True.)
!               Reset Cost Value = Logical (default = .True.)
!               Cost Filename = File (default = 'CostOfT.dat')
!               Cost Variable Name = String (default= 'CostValue')
!
!      Variables
!                SSAVelocity (solution of the SSA pb;DOFs== Pb Dimension)
!                Velocityb (forcing for the adjoint pb;DOFs== Pb dimension)
!                Zb (bottom elevation)
!                DJDZb (gradient of J wr. to Zb, as J is fn of H=Zs-Zb)
!                Zs (surface elevation)
!                DJDZs (gradient of J wr. to Zs, as J is fn of H=Zs-Zb)
!
!     In body Forces:
!               Top Surface Accumulation = Real (default= 0)
!               Bottom Surface Accumulation = Real (default = 0)
!               Observed dhdt = Real (default = 0)
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
  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName,VarSolName
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: TimeVar,CostVar

  TYPE(Variable_t), POINTER :: VelocitySol,VelocitybSol
  REAL(KIND=dp), POINTER :: Velocity(:),Vb(:)
  INTEGER, POINTER :: VeloPerm(:),VbPerm(:)

  TYPE(Variable_t), POINTER :: ZbSol,ZsSol
  TYPE(Variable_t), POINTER :: DJDZbSol,DJDZsSol
  REAL(KIND=dp), POINTER :: Zb(:),Zs(:)
  REAL(KIND=dp), POINTER :: DJDZb(:),DJDZs(:)
  INTEGER, POINTER :: ZbPerm(:),ZsPerm(:)
  INTEGER, POINTER :: DJDZbPerm(:),DJDZsPerm(:)

  TYPE(ValueList_t), POINTER :: SolverParams,BodyForce

  TYPE(Nodes_t) :: ElementNodes

  TYPE(GaussIntegrationPoints_t) :: IntegStuff

  INTEGER, POINTER :: NodeIndexes(:)
  integer :: i,j,k,l,t,n,NMAX,DIM,ierr,c

  real(kind=dp) :: Cost,Cost_S,Costb,Lambda
  real(kind=dp) :: u,v,w,s,coeff,SqrtElementMetric,x
  REAL(KIND=dp) :: Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)
  REAL(KIND=dp),dimension(:),allocatable,SAVE ::  NodeSMB,NodeDHDT,NodeH
  REAL(KIND=dp),dimension(:,:),allocatable,SAVE :: Velo
  REAL(KIND=dp) :: smb,dhdt,h,ugrdh,divu,gradh(2),Vgauss(2)

  LOGICAL :: ComputeDJDZb,ComputeDJDZs,ResetCost
  Logical :: Firsttime=.true.,Found,Parallel,stat,Gotit
  CHARACTER*10 :: date,temps

  save Firsttime,Parallel 
  save SolverName,CostSolName,CostFile
  save ElementNodes


  SolverParams => GetSolverParams()
  DIM=GetInteger(SolverParams ,'Problem Dimension',Found)
  If (.NOT.Found) then
     CALL WARN(SolverName,'Keyword >Problem Dimension< not found, assume DIM = CoordinateSystemDimension()')
     DIM = CoordinateSystemDimension()
  Endif

   Lambda =  GetConstReal( SolverParams,'Lambda', Found)
   IF(.NOT.Found) THEN
           CALL WARN(SolverName,'Keyword >Lambda< not found  in section >Solver<')
           CALL WARN(SolverName,'Taking default value Lambda=1.0')
           Lambda = 1.0
   End if

  ComputeDJDZb =  GetLogical( SolverParams,'Compute DJDZb', Found)
                    IF(.NOT.Found) ComputeDJDZb=.True.
  ComputeDJDZs =  GetLogical( SolverParams,'Compute DJDZs', Found)
                    IF(.NOT.Found) ComputeDJDZb=.True.
  ResetCost =  GetLogical( SolverParams,'Reset Cost Value', Found)
                    IF(.NOT.Found) ResetCost=.True.
  

  If (Firsttime) then
    N = model % MaxElementNodes
    allocate(ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N))
    allocate(NodeSMB(N),NodeDHDT(N),Velo(2,N),NodeH(N))

!!!!!!! Check for parallel run 
    Parallel = .FALSE.
    IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
            IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
                    Parallel = .TRUE.
            END IF
    END IF

     WRITE(SolverName, '(A)') 'CostSolver_Adjoint'

!!!!!!!!!!! get Solver Variables

  CostFile = ListGetString(Solver % Values,'Cost Filename',Found )
    IF (.NOT. Found) CostFile = DefaultCostFile
    CALL DATE_AND_TIME(date,temps)
    If (Parallel) then
        if (ParEnv % MyPe.EQ.0) then
           OPEN (12, FILE=CostFile)
                   write(12,1000) date(5:6),date(7:8),date(1:4),temps(1:2),temps(3:4),temps(5:6)
                   write(12,1001) Lambda
                   write(12,'(A)') '# iter, Jdiv'
           CLOSE(12)
         End if
    Else
           OPEN (12, FILE=CostFile)
                   write(12,1000) date(5:6),date(7:8),date(1:4),temps(1:2),temps(3:4),temps(5:6)
                   write(12,1001) Lambda
                   write(12,'(A)') '# iter, Jdiv'
           CLOSE(12)
    End if

   CostSolName =  GetString( SolverParams,'Cost Variable Name', Found)
          IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Cost Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >CostValue<')
                    WRITE(CostSolName,'(A)') 'CostValue'
          END IF
  
  !!! End of First visit
    Firsttime=.false.
  Endif


    VelocitySol => VariableGet( Solver % Mesh % Variables, 'SSAVelocity'  )
    IF ( ASSOCIATED( VelocitySol ) ) THEN
            Velocity => VelocitySol % Values
            VeloPerm => VelocitySol % Perm
    ELSE
            WRITE(Message,'(A)') &
                               'No variable > SSAVelocity < found'
            CALL FATAL(SolverName,Message)
    END IF  
    c=DIM  ! size of the velocity variable
    IF (VelocitySol % DOFs.NE.c) then
           WRITE(Message,'(A,I1,A,I1)') &
            'Variable SSAVelocity has ',VelocitySol % DOFs,' DOFs, should be',c
            CALL FATAL(SolverName,Message)
    End If

    VelocitybSol => VariableGet( Solver % Mesh % Variables, 'Velocityb'  )
    IF ( ASSOCIATED( VelocitybSol ) ) THEN
            Vb => VelocitybSol % Values
            VbPerm => VelocitybSol % Perm
    ELSE
            WRITE(Message,'(A)') &
                               'No variable > Velocityb < found'
            CALL FATAL(SolverName,Message)
    END IF  
    IF (VelocitybSol % DOFs.NE.c) then
           WRITE(Message,'(A,I1,A,I1)') &
            'Variable Velocityb has ',VelocitybSol % DOFs,' DOFs, should be',c
            CALL FATAL(SolverName,Message)
    End If
    if (ResetCost) Vb=0.0_dp

    ZbSol => VariableGet( Solver % Mesh % Variables, 'Zb' )
    IF (ASSOCIATED(ZbSol)) THEN
           Zb => ZbSol % Values
           ZbPerm => ZbSol % Perm
    ELSE
           CALL FATAL(SolverName,'Could not find variable >Zb<')
    END IF
    IF (ComputeDJDZb) Then
       DJDZbSol => VariableGet( Solver % Mesh % Variables, 'DJDZb' )
       IF (ASSOCIATED(DJDZbSol)) THEN
           DJDZb => DJDZbSol % Values
           DJDZbPerm => DJDZbSol % Perm
       ELSE
           CALL FATAL(SolverName,'Could not find variable >DJDZb<')
       END IF
       !!!!! Reset DJDZ to 0 HERE
       DJDZb=0._dp
    ENDIF

    ZsSol => VariableGet( Solver % Mesh % Variables, 'Zs' )
    IF (ASSOCIATED(ZsSol)) THEN
         Zs => ZsSol % Values
         ZsPerm => ZsSol % Perm
    ELSE
        CALL FATAL(SolverName,'Could not find variable >Zs<')
    END IF
    IF (ComputeDJDZs) Then
       DJDZsSol => VariableGet( Solver % Mesh % Variables, 'DJDZs' )
       IF (ASSOCIATED(DJDZsSol)) THEN
           DJDZs => DJDZsSol % Values
           DJDZsPerm => DJDZsSol % Perm
       ELSE
           CALL FATAL(SolverName,'Could not find variable >DJDZs<')
       END IF
       !!!!! Reset DJDZ to 0 HERE
       DJDZs=0._dp
    ENDIF


    Cost=0._dp

    DO t=1,Solver % NumberOfActiveElements
       Element => GetActiveElement(t)
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

      BodyForce => GetBodyForce()

     ! nodal values of SMB
      NodeSMB=0.0_dp
      NodeSMB(1:n) = ListGetReal( BodyForce, 'Top Surface Accumulation', n, NodeIndexes, GotIt)
          IF (.NOT.GotIt) Then
                  WRITE(Message,'(A)') &
                     'No variable >Top Surface Accumulation< found in "Body Forces" default is 0'
                  CALL Info(SolverName,Message,level=6)
          END IF 
      NodeSMB(1:n) = NodeSMB(1:n) + ListGetReal( BodyForce, 'Bottom Surface Accumulation', n, NodeIndexes, GotIt)
          IF (.NOT.GotIt) Then
                  WRITE(Message,'(A)') &
                     'No variable >Top Surface Accumulation< found in "Body Forces" default is 0'
                  CALL Info(SolverName,Message,level=6)
          END IF 

    ! nodal values of observed dhdt
      NodeDHDT=0._dp
      NodeDHDT(1:n) = ListGetReal( BodyForce, 'Observed dhdt', n, NodeIndexes, GotIt)
          IF (.NOT.GotIt) Then
                  WRITE(Message,'(A)') & 
                     'No variable >Observed dhdt< found in "Body Forces"; default is 0'
                  CALL Info(SolverName,Message,level=6)
          END IF 

     ! nodal values of velocities
      Velo=0._dp
      Do i=1,DIM
         Velo(i,1:n)=Velocity(DIM*(VeloPerm(NodeIndexes(1:n))-1)+i)
      End DO
     ! get Nodal values of H
      NodeH(1:n)=Zs(ZsPerm(NodeIndexes(1:n)))-Zb(ZbPerm(NodeIndexes(1:n)))

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

          x = SUM( ElementNodes % x(1:n) * Basis(1:n) )
          s = 1.0d0

          IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
             s = 2.0d0 * PI * x 
          END IF
          s = s * SqrtElementMetric * IntegStuff % s(i)
          
          !get smb at IP
          smb = SUM(NodeSMB(1:n)  * Basis(1:n))
          !observed dhdt at IP
          dhdt = SUM(NodeDHDT(1:n) * Basis(1:n))
          !h at IP
          h = SUM(NodeH(1:n)  * Basis(1:n))

          !u.grad H
          ugrdh=0._dp
          Do j=1,DIM
             gradh(j) = SUM( dBasisdx(1:n,j)*NodeH(1:n) )
             Vgauss(j) = SUM( Basis(1:n)*(Velo(j,1:n)) )
             ugrdh=ugrdh+gradh(j)*Vgauss(j)
          End do

          divu=0._dp
          Do j=1,DIM
             divu = divu + SUM( dBasisdx(1:n,j)*(Velo(j,1:n)))
          End Do

          coeff=dhdt+ugrdh+h*divu-smb

          Cost=Cost+0.5*coeff*coeff*s
         
         ! compute the derivatives (NOTE Cost=Lambda*Cost at the end for output
         ! reasons)
          Costb=0.5*s*Lambda*2.0*coeff
          
          Do l=1,n
            IF (ComputeDJDZb) &
             DJDZb(DJDZbPerm(NodeIndexes(l)))=DJDZb(DJDZbPerm(NodeIndexes(l))) - & ! - car h=Zs-Zb
                                             Costb*divu*Basis(l)
            IF (ComputeDJDZs) &
             DJDZs(DJDZsPerm(NodeIndexes(l)))=DJDZs(DJDZsPerm(NodeIndexes(l))) + & ! + car h=Zs-Zb
                                             Costb*divu*Basis(l)
             Do j=1,DIM
                k=(VbPerm(NodeIndexes(l))-1)*c+j
                Vb(k)=Vb(k) + Costb * Basis(l) * gradh(j) +&
                              Costb * h * dBasisdx(l,j)

               IF (ComputeDJDZb) &
                DJDZb(DJDZbPerm(NodeIndexes(l)))=DJDZb(DJDZbPerm(NodeIndexes(l))) - &  ! - car h=Zs-Zb
                   Costb*Vgauss(j)*dBasisdx(l,j)
               IF (ComputeDJDZs) &
                DJDZs(DJDZsPerm(NodeIndexes(l)))=DJDZs(DJDZsPerm(NodeIndexes(l))) + &  ! + car h=Zs-Zb
                   Costb*Vgauss(j)*dBasisdx(l,j)
             End do !j
          End do !l

        End do !IPs

    End do !elements

   

    TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )

    IF (Parallel) THEN
           CALL MPI_ALLREDUCE(Cost,Cost_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)

          IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) then
                 OPEN (12, FILE=CostFile,POSITION='APPEND')
                 write(12,'(e13.5,2x,e15.8)') TimeVar % Values(1),Cost_S
                 CLOSE(12)
          End if

          Cost_S = Lambda * Cost_S

          CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
          IF (ASSOCIATED(CostVar)) THEN
               IF (ResetCost) then
                 CostVar % Values(1)=Cost_S
               Else
                 CostVar % Values(1)=CostVar % Values(1)+Cost_S
               Endif
          END IF
   ELSE
            OPEN (12, FILE=CostFile,POSITION='APPEND')
                  write(12,'(e13.5,2x,e15.8)') TimeVar % Values(1),Cost
            close(12)

            Cost = Lambda * Cost

            CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
            IF (ASSOCIATED(CostVar)) THEN
                IF (ResetCost) then
                    CostVar % Values(1)=Cost
                Else
                     CostVar % Values(1)=CostVar % Values(1)+Cost
                Endif
            END IF
   END IF
   
   Return

 1000  format('#date,time,',a1,'/',a1,'/',a4,',',a2,':',a2,':',a2)
 1001  format('#lambda,',e15.8)
!------------------------------------------------------------------------------
END SUBROUTINE AdjointSSA_CostFluxDivSolver
!------------------------------------------------------------------------------
! *****************************************************************************
