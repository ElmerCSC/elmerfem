!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> The current density given as a source must be divergence free to allow a 
!> hope for a solution. This solver may be used to enforce a given current 
!> density to be divergence free by solving for an equivalent potential
!> so that the vector field is a gradient of the potential.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE JfixPotentialSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL ::  Transient
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: SolverParams
  INTEGER :: i,j,k,n,m,dim,dofs
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Matrix_t), POINTER :: A,B
  REAL(KIND=dp) :: Norm
  LOGICAL:: AutomatedBCs, Found
  INTEGER, ALLOCATABLE :: Def_Dofs(:,:,:)
  CHARACTER(LEN=MAX_NAME_LEN):: Equation
  INTEGER, POINTER :: Perm(:)
  REAL(KIND=dp), POINTER :: fixpot(:)
  TYPE(Variable_t), POINTER :: fixJpot, svar, IterV 
  
  dim = CoordinateSystemDimension()
  Mesh => GetMesh()
  B => GetMatrix()
  SolverParams => GetSolverParams()

  fixJpot => VariableGet( Mesh % Variables, 'Jfix')

  IF( ASSOCIATED(fixJPot)) THEN
    Perm => fixJPot % Perm
  ELSE
    ALLOCATE(Perm(SIZE(Solver % Variable % Perm)))
  END IF
  Perm = 0
  dofs = Solver % Variable % DOFs

  Equation=GetString(SolverParams,'Equation',Found)

  n=SIZE(Solver % Def_Dofs,1)
  m=SIZE(Solver % Def_Dofs,2)
  k=SIZE(Solver % Def_Dofs,3)
  ALLOCATE(Def_Dofs(n,m,k))
  Def_Dofs=Solver % Def_Dofs
  Solver % Def_Dofs=0
  Solver % Def_Dofs(:,:,1)=1

  A => CreateMatrix( CurrentModel, Solver, Solver % Mesh, &
     Perm, dofs, MATRIX_CRS, .TRUE., Equation, .FALSE., .FALSE.,&
     NodalDofsOnly = .TRUE.)
  n = A % NumberOfRows
  IF (dofs>1) A % COMPLEX = .TRUE.

  IF (.NOT.ASSOCIATED(fixJPot)) THEN
    ALLOCATE(fixpot(n)); fixpot=0._dp

    CALL VariableAddVector( Mesh % Variables, Mesh, &
          Solver,'Jfix',dofs,fixpot,Perm)

    fixJpot => VariableGet(Mesh % Variables, 'Jfix')
  END IF

  AutomatedBCs = GetLogical( SolverParams, &
       'Automated Source Projection BCs', Found )
  IF (.NOT. Found) AutomatedBCs = .TRUE.

  svar => Solver % Variable
  Solver % Variable => fixJpot

  Solver % Matrix => A
  ALLOCATE(A % RHS(n))
  IF(ParEnv % PEs>1) CALL ParallelInitMatrix(Solver,A)

  CALL ListSetNameSpace('jfix:')
  CALL ListAddString(SolverParams,'Jfix: Linear System Solver', 'Iterative')
  CALL ListAddString(SolverParams,'Jfix: Linear System Iterative Method', 'BiCGStab')

  CALL ListAddLogical(SolverParams,'Jfix: Linear System Use HYPRE', .FALSE.)
  CALL ListAddLogical(SolverParams,'Jfix: Use Global Mass Matrix',.FALSE.)

  IF (.NOT.ListCheckPresent(SolverParams,'Jfix: Linear System Preconditioning')) &
    CALL ListAddString(SolverParams,'Jfix: Linear System Preconditioning', 'Ilu')

  IF (.NOT.ListCheckPresent(SolverParams,'Jfix: Linear System Convergence Tolerance')) &
    CALL ListAddConstReal(SolverParams,'Jfix: Linear System Convergence Tolerance', &
      0.001_dp*GetCReal(SolverParams,'Linear System Convergence Tolerance', Found))


  CALL ListAddLogical(SolverParams,'Jfix: Skip Compute Nonlinear Change',.TRUE.)

  CALL DefaultInitialize()
  CALL BulkAssembly()
  IF ( (.NOT. A % COMPLEX) .AND. (.NOT. AutomatedBCs) ) CALL BCVectorAssembly()
  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  Norm = DefaultSolve()

  Solver % Matrix => B
  Solver % Variable => svar
  CALL ListSetNameSpace('')

  CALL FreeMatrix(A)
  Solver % Def_Dofs=Def_Dofs
  
  IterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
  IterV % Values(1) = 1
CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
!------------------------------------------------------------------------------
    IMPLICIT NONE       
    COMPLEX(KIND=dp), ALLOCATABLE :: STIFF_C(:,:), FORCE_C(:)
    REAL(KIND=dp), ALLOCATABLE :: STIFF_R(:,:), FORCE_R(:)
    INTEGER :: elem,t,i,j,k,p,q,n,nd, Rank
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    LOGICAL :: Found
    TYPE(Element_t), POINTER :: Element, P1,P2
    TYPE(ValueList_t), POINTER ::  BodyForce

    REAL(KIND=dp) :: weight,detJ,Source,L_R(3),Nrm(3), Jfluxeps,rabs
    REAL(KIND=dp), ALLOCATABLE :: Load_R(:,:)

    COMPLEX(KIND=dp) :: L_C(3)
    COMPLEX(KIND=dp), ALLOCATABLE :: Load_C(:,:)

    TYPE(Matrix_t), POINTER :: gm
    INTEGER, ALLOCATABLE :: narr(:)
    INTEGER :: l,m,ncols,ierr,status(MPI_STATUS_SIZE)

    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    
    SAVE Nodes
    
    n = MAX(Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes)
    IF (A % COMPLEX) THEN
      ALLOCATE( Load_c(3,n), STIFF_C(n,n), FORCE_C(n) )
    ELSE
      ALLOCATE( Load_r(3,n), STIFF_R(n,n), FORCE_R(n) )
    END IF
    ALLOCATE( Basis(n), dBasisdx(n,3) )

    DO elem = 1,GetNOFActive()
      ! Element information
      ! ---------------------
      Element => GetActiveElement(elem)
      CALL GetElementNodes( Nodes )
      n  = GetElementNOFNodes()
      nd = n

      IF (A % COMPLEX ) THEN
        Load_C =0._dp
        STIFF_C=0._dp
        FORCE_C=0._dp
      ELSE
        Load_R =0._dp
        STIFF_R=0._dp
        FORCE_R=0._dp
      END IF

      BodyForce => GetBodyForce()
      IF (ASSOCIATED(BodyForce)) THEN
        IF (A % COMPLEX ) THEN
          CALL GetComplexVector(BodyForce,Load_C(1:3,1:n),'Current Density',Found)
        ELSE
          CALL GetRealVector(BodyForce,Load_R(1:3,1:n),'Current Density',Found)
        END IF
      END IF

      IntegStuff = GaussPoints( Element )
      DO t=1,IntegStuff % n
        Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
                IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )
        
        Weight = IntegStuff % s(t) * detJ
        IF (A % COMPLEX) THEN
          DO p=1,n
            DO q=1,n
              STIFF_C(p,q) = STIFF_C(p,q) + Weight * SUM(dBasisdx(q,:)*dBasisdx(p,:))
            END DO
          END DO
          L_C = MATMUL( Load_C(:,1:n), Basis(1:n) )
          FORCE_C(1:n) = FORCE_C(1:n) + MATMUL(dBasisdx(1:n,:),L_C)*Weight
        ELSE
          DO p=1,n
            DO q=1,n
              STIFF_R(p,q) = STIFF_R(p,q) + Weight * SUM(dBasisdx(q,:)*dBasisdx(p,:))
            END DO
          END DO
          L_R = MATMUL( Load_R(:,1:n), Basis(1:n) )
          FORCE_R(1:n) = FORCE_R(1:n) + MATMUL(dBasisdx(1:n,:),L_R)*Weight
        END IF
      END DO

      IF (A % COMPLEX) THEN
        DO i=n+1,nd
          STIFF_C(i,i)=1._dp
        END DO
        CALL DefaultUpdateEquations(STIFF_C, FORCE_C)
      ELSE
        DO i=n+1,nd
          STIFF_R(i,i)=1._dp
        END DO
        CALL DefaultUpdateEquations(STIFF_R, FORCE_R)
      END IF
    END DO


    Jfluxeps = GetCReal(SolverParams, 'J normal eps', Found)
    IF (.NOT. Found) Jfluxeps = 1.0d-2

    IF (ParEnv%Pes>1) THEN
      gm => AllocateMatrix()
      gm % FORMAT = MATRIX_LIST
    END IF

    ! At outer boundaries, if J has component normal to the
    ! surface, set Fix_pot=0:
    ! THIS MAY NOT DO WHAT IS PHYSICALLY SOUND: It may be better
    ! to disable this by the command 
    ! Automated Projection BCs = Logical False
    ! and set BCs explicitly in the sif file    
    ! ------------------------------------------------------
    IF (AutomatedBCs) THEN
       DO i=1,GetNOFBoundaryElements()
          Element => GetBoundaryElement(i)
          n = GetElementNOFNodes()
          IF (.NOT.ActiveBoundaryElement()) CYCLE


          P1 => Element % BoundaryInfo % Left
          P2 => Element % BoundaryInfo % Right
          IF(ASSOCIATED(P1).AND.ASSOCIATED(P2)) THEN
             IF (ALL(Perm(P1 % NodeIndexes)>0).AND.&
                  ALL(Perm(P2 % NodeIndexes)>0)) CYCLE
          ELSE IF (.NOT.ASSOCIATED(P1).AND..NOT.ASSOCIATED(P2)) THEN
             CYCLE
          END IF

          IF (ASSOCIATED(P1)) THEN
             BodyForce => GetBodyForce(P1)
             IF (Element % TYPE % DIMENSION<P1 % TYPE % DIMENSION-1) CYCLE
          ELSE
             BodyForce => GetBodyForce(P2)
             IF (Element % TYPE % DIMENSION<P2 % TYPE % DIMENSION-1) CYCLE
          END IF
          IF (.NOT. ASSOCIATED(BodyForce)) CYCLE

          CALL GetElementNodes(Nodes)
          Nrm = NormalVector(Element,Nodes,0._dp,0._dp)

          IF (A % COMPLEX ) THEN
             CALL GetComplexVector( BodyForce, Load_C(1:3,1:n), 'Current Density', Found )
             L_C(1) = SUM(Load_C(1,1:n))/n
             L_C(2) = SUM(Load_C(2,1:n))/n
             L_C(3) = SUM(Load_C(3,1:n))/n

             L_R = REAL(L_C)
             IF (ANY(L_R /= 0.0d0)) L_R = L_R / SQRT(SUM(L_R**2))
             IF (ABS(SUM(L_R*Nrm))<Jfluxeps) THEN
                L_R = AIMAG(L_C)
                IF (ANY(L_R /= 0.0d0)) L_R = L_R / SQRT(SUM(L_R**2))
                IF (ABS(SUM(L_R*Nrm))<Jfluxeps) CYCLE
             END IF
          ELSE
             CALL GetRealVector( BodyForce, Load_R(1:3,1:n), 'Current Density', Found )

             L_R(1) = SUM(Load_R(1,1:n))/n
             L_R(2) = SUM(Load_R(2,1:n))/n
             L_R(3) = SUM(Load_R(3,1:n))/n
             rabs = SQRT(SUM(L_R**2))
             IF(rabs/=0) L_R = L_R / rabs
             IF (ABS(SUM(L_R*Nrm))<Jfluxeps) CYCLE
          END IF

          DO j=1,dofs
             DO k=1,n
                p = Element % NodeIndexes(k) 
                p = dofs*(Perm(p)-1)+j

#ifdef OLD_DIRICHLET
                ! In parallel collect shared constrained dofs and send to
                ! other owners...
                ! -------------------------------------------------------
                IF (Parenv % PEs>1) THEN
                   IF(A % ParallelInfo % INTERFACE(p)) THEN
                      DO l=1,SIZE(A % ParallelInfo % Neighbourlist(p) % Neighbours)
                         m=A % ParallelInfo % NeighbourList(p) % Neighbours(l)
                         IF(m/=ParEnv % MyPE) THEN
                            CALL List_AddToMatrixElement( gm % ListMatrix, m+1, &
                                 A % ParallelInfo % GlobalDOFs(p), 1._dp)
                         ELSE
                            CALL CRS_SetSymmDirichlet(A,A % RHS,p,1._dp)
                         END IF
                      END DO
                   ELSE
                      CALL CRS_SetSymmDirichlet(A,A % RHS,p,1._dp)
                   END IF
                ELSE
                   CALL CRS_SetSymmDirichlet(A,A % RHS,p,1._dp)
                END IF
#else
                CALL UpdateDirichletDof(A, p, 0._dp)
#endif
             END DO
          END DO
       END DO
    END IF

    ! In parallel apply constraints spotted by other partitions
    ! ---------------------------------------------------------
#ifdef OLD_DIRICHLET
    IF(ParEnv % PEs>1) THEN
      IF(ASSOCIATED(gm % ListMatrix)) CALL List_toCRSMatrix(gm)

      DO i=1,ParEnv % PEs
        IF(i-1==Parenv % myPE.OR..NOT.ParEnv % Active(i)) CYCLE

        IF(i<=gm % NumberOfRows) THEN
          ncols=gm % Rows(i+1)-gm % Rows(i)
          CALL MPI_BSEND(ncols,1,MPI_INTEGER,i-1,2120,ELMER_COMM_WORLD,ierr)
          IF(ncols>0) THEN
            CALL MPI_BSEND(gm%cols(gm%rows(i):gm%rows(i+1)-1),ncols, &
                 MPI_INTEGER,i-1,2121,ELMER_COMM_WORLD,ierr)
          END IF
        ELSE
           ncols=0
          CALL MPI_BSEND(ncols,1,MPI_INTEGER,i-1,2120,ELMER_COMM_WORLD,ierr)
        END IF
      END DO

      CALL FreeMatrix(gm)

      DO i=1,ParEnv % PEs
        IF(i-1==Parenv % myPE.OR..NOT.ParEnv % Active(i)) CYCLE

        CALL MPI_RECV(ncols,1,MPI_INTEGER,i-1,2120,ELMER_COMM_WORLD,status,ierr)
        IF(ncols>0) THEN
          ALLOCATE(narr(ncols))
          CALL MPI_RECV(narr,ncols,MPI_INTEGER,i-1,2121,ELMER_COMM_WORLD,status,ierr)
          DO j=1,ncols
            k = SearchIAItem(A % NumberOfRows,A % ParallelInfo % GlobalDOFs, narr(j))
            IF(k>0) THEN
              CALL CRS_SetSymmDirichlet(A,A % RHS,k,1._dp)
            END IF
          END DO
          DEALLOCATE(narr)
        END IF
      END DO
    END IF
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE BCVectorAssembly()
!------------------------------------------------------------------------------
!   This is intended to alter the natural boundary conditions of the potential associated
!   with the Helmholtz projection if the source electric current density is
!   generated by applying the StatCurrentSolve module.
!-------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), ALLOCATABLE :: FORCE_R(:), LOAD(:), Basis(:), dBasisdx(:,:)
    TYPE(GaussIntegrationPoints_t), TARGET :: IP  
    TYPE(Element_t), POINTER :: Element
    TYPE(ValueList_t), POINTER :: BC
    INTEGER :: Active, k, p, t, nd, n
    REAL(KIND=dp) :: detJ, Jn
    LOGICAL :: Found, Stat

    TYPE(Nodes_t), SAVE :: Nodes
!-------------------------------------------------------------------------------

    n = MAX(Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes)  
    ALLOCATE( FORCE_R(n), LOAD(n), Basis(n), dBasisdx(n,3) )    

    Active = GetNOFBoundaryElements()
    DO k=1,Active
       Element => GetBoundaryElement(k)
       IF (.NOT. ActiveBoundaryElement()) CYCLE
       BC=>GetBC()
       IF (.NOT. ASSOCIATED(BC) ) CYCLE
       IF (GetElementFamily()==1) CYCLE

       nd = GetElementNOFDOFs(Element)
       n  = GetElementNOFNodes(Element)

       Load = 0.0d0
       FORCE_R =0.0d0

       Load(1:n) = GetReal( BC, 'Current Density', Found )
       Jn = SUM( Load(1:n)*Basis(1:n) ) 

       CALL GetElementNodes( Nodes )
       IP = GaussPoints( Element )
       DO t=1,IP % n
          stat = ElementInfo( Element, Nodes, IP % u(t), &
               IP % v(t), IP % w(t), detJ, Basis, dBasisdx )
          DO p=1,n
             FORCE_R(p) = FORCE_R(p) + Jn*Basis(p)*detJ*IP % s(t)
          END DO
       END DO

       CALL DefaultUpdateForce(FORCE_R,Element)
    END DO

    DEALLOCATE(FORCE_R,LOAD,Basis,dBasisdx)
!------------------------------------------------------------------------------
END SUBROUTINE BCVectorAssembly
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
END SUBROUTINE JfixPotentialSolver
!------------------------------------------------------------------------------
