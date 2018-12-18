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
! *  Original Date: 10 May 2000
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Initialization for the primary solver i.e. MeshSolver.
!------------------------------------------------------------------------------
 SUBROUTINE MeshSolver_Init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  LOGICAL ::  TransientSimulation
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  INTEGER :: dim
  LOGICAL :: Found, Calculate

  Params => GetSolverParams()
  dim = CoordinateSystemDimension()

  Calculate = ListGetLogical( Params,'Compute Mesh Velocity',Found ) 
  IF(.NOT. Found ) Calculate = .TRUE.

  IF( Calculate ) THEN
    IF( TransientSimulation ) THEN
      IF( dim == 2 ) THEN
        CALL ListAddString( Params,&
            NextFreeKeyword('Exported Variable',Params),&
            '-dofs 2 Mesh Velocity')        
      ELSE
        CALL ListAddString( Params,&
            NextFreeKeyword('Exported Variable',Params),&
            '-dofs 3 Mesh Velocity')                  
      END IF
    END IF    
  END IF

END SUBROUTINE MeshSolver_Init


!------------------------------------------------------------------------------
!> Subroutine for extending displacement in mesh smoothly over 
!> the domain. The intended use of the solver is in fluid-structure interaction, 
!> for example. In transient cases the solver also computes the mesh velocity. 
!> This is a dynamically loaded solver with a standard interface.
!> May be also loaded internally to mimic the old static implementation. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
 SUBROUTINE MeshSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  LOGICAL ::  TransientSimulation
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  INTEGER :: i,j,k,l,m,n,nd,nb,t,DOFs,STDOFs,LocalNodes,istat

  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t),POINTER :: Material, BC

  REAL(KIND=dp) :: RelativeChange, UNorm, PrevUNorm,  maxu

  TYPE(Variable_t), POINTER :: StressSol, MeshSol

  REAL(KIND=dp), POINTER :: MeshUpdate(:),Displacement(:), &
       MeshVelocity(:)

  INTEGER, POINTER :: TPerm(:), MeshPerm(:), StressPerm(:), MeshVeloPerm(:)

  LOGICAL :: AllocationsDone = .FALSE., Isotropic = .TRUE., &
            GotForceBC, Found, ComputeMeshVelocity, DisplaceFirst, &
            SkipFirstMeshVelocity = .FALSE., FirstTime = .TRUE., &
            SkipDisplace 
  REAL(KIND=dp),ALLOCATABLE:: STIFF(:,:),&
       LOAD(:,:),FORCE(:), ElasticModulus(:,:,:),PoissonRatio(:), &
       Alpha(:,:), Beta(:)

  INTEGER :: dim
  
  SAVE STIFF, LOAD, FORCE, MeshVelocity, MeshVeloPerm, AllocationsDone, &
       ElasticModulus, PoissonRatio, TPerm, Alpha, Beta, &
       SkipFirstMeshVelocity, FirstTime

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at,at0
#else
  REAL(KIND=dp) :: at,at0,CPUTime,RealTime
#endif
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

  NULLIFY( MeshVelocity ) 
  IF ( TransientSimulation ) THEN
    MeshSol      => VariableGet( Solver % Mesh % Variables, 'Mesh Velocity' )
    IF( ASSOCIATED( MeshSol ) ) THEN
      MeshVelocity => MeshSol % Values
      MeshVeloPerm => MeshSol % Perm
    END IF
  END IF

  MeshSol => Solver % Variable
  MeshPerm      => MeshSol % Perm
  STDOFs        =  MeshSol % DOFs
  MeshUpdate    => MeshSol % Values

  LocalNodes = COUNT( MeshPerm > 0 )
  IF ( LocalNodes <= 0 ) RETURN

  dim = CoordinateSystemDimension()
  IF( STDOFs < dim ) THEN
    CALL Warn('MeshSolve','Displacement dofs smaller than dim, assuming reduced dimension!')

    ! We can solve the equation also on a DIM-1 dimensional boundary of DIM dimensional object
    dim = STDOFs
  END IF
    
  
!------------------------------------------------------------------------------

  StressSol => VariableGet( Solver % Mesh % Variables, 'Displacement' )
  IF( ListGetLogical( Solver % Values,'Ignore Displacement',Found) ) THEN
    NULLIFY( StressSol )
  END IF

  DisplaceFirst = ListGetLogical( Solver % Values,'First Time Non-Zero', Found)
  SkipDisplace = ListGetLogical( Solver % Values,'Skip Displace Mesh',Found )
  
  IF( SkipDisplace ) THEN
    CALL Info('MeshSolver','Skipping the displacement of mesh!',Level=5)
  ELSE IF ( ASSOCIATED( StressSol ) )  THEN
     StressPerm   => StressSol % Perm
     STDOFs       =  StressSol % DOFs
     Displacement => StressSol % Values

     IF( .NOT.AllocationsDone .OR. Solver % Mesh % Changed ) THEN
        IF ( AllocationsDone ) DEALLOCATE( TPerm )

        ALLOCATE( TPerm( SIZE(MeshPerm) ), STAT=istat )
        IF ( istat /= 0 ) THEN
           CALL Fatal( 'MeshSolve', 'Memory allocation error.' )
        END IF
     END IF

     TPerm = MeshPerm
     DO i=1,SIZE( MeshPerm )
        IF ( StressPerm(i) /= 0 .AND. MeshPerm(i) /= 0 ) TPerm(i) = 0
     END DO

     IF ( AllocationsDone .OR. DisplaceFirst ) THEN
        CALL DisplaceMesh( Solver % Mesh, MeshUpdate, -1, TPerm,   STDOFs )
     END IF
     CALL DisplaceMesh( Solver % Mesh, Displacement,  -1, StressPerm, STDOFs )
  ELSE
     IF ( AllocationsDone .OR. DisplaceFirst ) THEN
        CALL DisplaceMesh( Solver % Mesh, MeshUpdate, -1, MeshPerm, STDOFs )
     END IF
  END IF

!------------------------------------------------------------------------------

  UNorm = Solver % Variable % Norm
!------------------------------------------------------------------------------
! Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
     N = Solver % Mesh % MaxElementDOFs

     IF ( AllocationsDone ) THEN
        DEALLOCATE(  ElasticModulus, PoissonRatio, &
             FORCE, Alpha, Beta, STIFF, LOAD, STAT=istat )
     END IF

     ALLOCATE( &
          Alpha(3,N), Beta(N), &
          ElasticModulus( 6,6,N ), PoissonRatio( N ), &
          FORCE( STDOFs*N ), STIFF( STDOFs*N,STDOFs*N ),  &
          LOAD( 4,N ),STAT=istat )

     IF ( istat /= 0 ) THEN
        CALL Fatal( 'MeshSolve', 'Memory allocation error.' )
     END IF

!------------------------------------------------------------------------------
     AllocationsDone = .TRUE.
!------------------------------------------------------------------------------
  END IF
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Do some additional initialization, and go for it
!------------------------------------------------------------------------------
  at  = CPUTime()
  at0 = RealTime()

  CALL Info( 'MeshSolve', ' ', Level=4 )
  CALL Info( 'MeshSolve', '-------------------------------------', Level=4 )
  CALL Info( 'MeshSolve', 'MESH UPDATE SOLVER:', Level=4 )
  CALL Info( 'MeshSolve', '-------------------------------------', Level=4 )
  CALL Info( 'MeshSolve', ' ', Level=4 )
  CALL Info( 'MeshSolve', 'Starting assembly...', Level=4 )
!------------------------------------------------------------------------------
  CALL DefaultInitialize()
!------------------------------------------------------------------------------

  
  DO t=1,Solver % NumberOfActiveElements

     IF ( RealTime() - at0 > 1.0 ) THEN
        WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
             (Solver % NumberOfActiveElements-t) / &
             (1.0*Solver % NumberOfActiveElements)), ' % done'

        CALL Info( 'MeshSolve', Message, Level=5 )
                     
        at0 = RealTime()
     END IF

     Element => GetActiveElement(t)
     nd = GetElementNOFDOFs()
     nb = GetElementNOFBDOFs()
     n  = GetElementNOFNodes()

     Material => GetMaterial()

     ElasticModulus(1,1,1:n) = GetReal( Material, &
          'Mesh Elastic Modulus', Found )
     IF ( .NOT. Found ) THEN
        ElasticModulus(1,1,1:n) = GetReal( Material, &
             'Youngs Modulus', Found )
     END IF
     IF ( .NOT. Found ) ElasticModulus(1,1,1:n) = 1.0d0

     PoissonRatio(1:n) = GetReal( Material, &
          'Mesh Poisson Ratio', Found )
     IF ( .NOT. Found ) THEN
        PoissonRatio(1:n) = GetReal( Material, &
          'Poisson Ratio', Found )
     END IF
     IF ( .NOT. Found ) PoissonRatio(1:n) = 0.25d0

!------------------------------------------------------------------------------
!    Get element local stiffness & mass matrices
!------------------------------------------------------------------------------
     CALL LocalMatrix( STIFF, FORCE, ElasticModulus, &
       PoissonRatio, .FALSE., Isotropic, Element, n, nd, nb )

!------------------------------------------------------------------------------
!    Update global matrices from local matrices 
!------------------------------------------------------------------------------
     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
  CALL DefaultFinishBulkAssembly()


!------------------------------------------------------------------------------
!     Neumann & Newton boundary conditions
!------------------------------------------------------------------------------
  DO t = 1, Solver % Mesh % NumberOfBoundaryElements

    Element => GetBoundaryElement(t)
    IF ( .NOT.ActiveBoundaryElement() ) CYCLE

    BC => GetBC()
    IF ( .NOT. ASSOCIATED(BC) ) CYCLE

!------------------------------------------------------------------------------
!        Force in given direction BC: \tau\cdot n = F
!------------------------------------------------------------------------------
     nd = GetElementNOFDOFs()
     n  = GetElementNOFNodes()
     nb = GetElementNOFBDOFs()

     LOAD = 0.0D0
     Alpha      = 0.0D0
     Beta       = 0.0D0

     GotForceBC = .FALSE.
     LOAD(1,1:n) =  GetReal( BC, 'Mesh Force 1', Found )
     GotForceBC = GotForceBC.OR.Found
     LOAD(2,1:n) =  GetReal( BC, 'Mesh Force 2', Found )
     GotForceBC = GotForceBC.OR.Found
     LOAD(3,1:n) =  GetReal( BC, 'Mesh Force 3', Found )
     GotForceBC = GotForceBC.OR.Found

     Beta(1:n) = GetReal( BC, 'Mesh Normal Force',Found )
     GotForceBC = GotForceBC.OR.Found

     IF ( .NOT.GotForceBC ) CYCLE

     CALL MeshBoundary( STIFF,FORCE, LOAD,Alpha,Beta,Element,n,nd,nb )

!------------------------------------------------------------------------------

     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
!------------------------------------------------------------------------------

  CALL DefaultFinishAssembly()
  CALL Info( 'MeshSolve', 'Assembly done', Level=4 )

!------------------------------------------------------------------------------
! Dirichlet boundary conditions
!------------------------------------------------------------------------------
  CALL DefaultDirichletBCs()

!------------------------------------------------------------------------------
  CALL Info( 'MeshSolve', 'Set boundaries done', Level=4 )
!------------------------------------------------------------------------------
! Solve the system and check for convergence
!------------------------------------------------------------------------------
  PrevUNorm = UNorm

  UNorm = DefaultSolve()

  IF ( UNorm + PrevUNorm  /= 0.0d0 ) THEN
     RelativeChange = 2*ABS( PrevUNorm - UNorm ) / (PrevUNorm + UNorm)
  ELSE
     RelativeChange = 0.0d0
  END IF

  WRITE( Message, * ) 'Result Norm   : ',UNorm
  CALL Info( 'MeshSolve', Message, Level=4 )
  WRITE( Message, * ) 'Relative Change : ',RelativeChange
  CALL Info( 'MeshSolve', Message, Level=4 )

  
  IF ( TransientSimulation ) THEN
    ComputeMeshVelocity = ListGetLogical( Solver % Values, 'Compute Mesh Velocity', Found )
    IF ( .NOT. Found ) ComputeMeshVelocity = .TRUE.
    SkipFirstMeshVelocity = .FALSE.
    IF (ComputeMeshVelocity .AND. FirstTime) THEN
       SkipFirstMeshVelocity = ListGetLogical( Solver % Values, 'Skip First Mesh Velocity', Found )
       IF (.NOT. Found ) THEN 
          SkipFirstMeshVelocity = .FALSE.
       ELSE
          CALL INFO('MeshSolve', 'Skipping computation of initial Mesh Velocity', Level=3)
       END IF
       FirstTime = .FALSE.
    END IF
    
    IF ( ComputeMeshVelocity .AND. (.NOT.(SkipFirstMeshVelocity)) ) THEN
      k = MIN( SIZE(Solver % Variable % PrevValues,2), Solver % DoneTime )
      
      j = ListGetInteger( Solver % Values,'Compute Mesh Velocity Order', Found)
      IF( Found ) THEN
        k = MIN( k, j )        
      ELSE
        k = 1
      END IF
      
      DOFs = MeshSol % DOFs

      DO i=1,Solver % Mesh % NumberOfNodes
         IF(MeshVeloPerm(i) <= 0) CYCLE

         j = MeshPerm(i)*DOFs
         l = MeshVeloPerm(i)*DOFs

         SELECT CASE(k)
         CASE(0)
            MeshVelocity = 0._dp
         CASE(1)
            DO m=0,DOFs-1
               MeshVelocity(l-m) = ( MeshUpdate(j-m) - &
                    Solver % Variable % PrevValues(j-m,1) ) / dt
            END DO
         CASE(2)
            DO m=0,DOFs-1
               MeshVelocity(l-m) = ( &
                    MeshUpdate(j-m) - (4.0d0/3.0d0)*Solver % Variable % PrevValues(j-m,1) &
                    + (1.0d0/3.0d0)*Solver % Variable % PrevValues(j-m,2) ) / dt
            END DO
         CASE DEFAULT
            DO m=0,DOFs-1
               MeshVelocity(l-m) = ( &
                    MeshUpdate(j-m) - (18.0d0/11.0d0)*Solver % Variable % PrevValues(j-m,1) &
                    + ( 9.0d0/11.0d0)*Solver % Variable % PrevValues(j-m,2) &
                    - ( 2.0d0/11.0d0)*Solver % Variable % PrevValues(j-m,3) ) / dt
            END DO
         END SELECT
      END DO
    ELSE IF( ASSOCIATED( MeshVelocity ) .AND. (.NOT.(SkipFirstMeshVelocity)) ) THEN 
      MeshVelocity = 0.0d0
    END IF
  END IF


  IF( SkipDisplace ) THEN
    CONTINUE
  ELSE IF ( ASSOCIATED( StressSol ) ) THEN
    CALL DisplaceMesh( Solver % Mesh, MeshUpdate,   1, TPerm,      STDOFs )
    CALL DisplaceMesh( Solver % Mesh, Displacement, 1, StressPerm, STDOFs, .FALSE.)
  ELSE
    CALL DisplaceMesh( Solver % Mesh, MeshUpdate,   1, MeshPerm,   STDOFs )
  END IF

  CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE LocalMatrix( STIFF,FORCE,NodalYoung, NodalPoisson, &
              PlaneStress, Isotropic, Element,n, nd, nb )
!------------------------------------------------------------------------------
     IMPLICIT NONE

     REAL(KIND=dp) :: NodalPoisson(:), NodalYoung(:,:,:)
     REAL(KIND=dp), TARGET :: STIFF(:,:), FORCE(:)

     INTEGER :: n,nd,nb

     TYPE(Element_t) :: Element
     LOGICAL :: PlaneStress, Isotropic
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(nd)
     REAL(KIND=dp) :: dBasisdx(nd,3),detJ

     REAL(KIND=dp) :: NodalLame1(n),NodalLame2(n),Lame1,Lame2, &
                      Poisson, Young

     REAL(KIND=dp), POINTER :: A(:,:)
     REAL(KIND=dp) :: s,u,v,w
     INTEGER :: i,j,k,p,q,t
  
     LOGICAL :: stat
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     TYPE(Nodes_t) :: Nodes
     SAVE  Nodes
!------------------------------------------------------------------------------

     CALL GetElementNodes( Nodes )
     
     IF ( PlaneStress ) THEN
        NodalLame1(1:n) = NodalYoung(1,1,1:n) * NodalPoisson(1:n) / &
               ((1.0d0 - NodalPoisson(1:n)**2))
     ELSE
        NodalLame1(1:n) = NodalYoung(1,1,1:n) * NodalPoisson(1:n) /  &
           ((1.0d0 + NodalPoisson(1:n)) * (1.0d0 - 2.0d0*NodalPoisson(1:n)))
     END IF

     NodalLame2(1:n) = NodalYoung(1,1,1:n) / (2* (1.0d0 + NodalPoisson(1:n)))

     STIFF = 0.0d0
     FORCE = 0.0d0

     ! Integration stuff:
     ! ------------------  
     IntegStuff = GaussPoints( Element )

     ! Now we start integrating:
     ! -------------------------
     DO t=1,IntegStuff % n
       u = IntegStuff % u(t)
       v = IntegStuff % v(t)
       w = IntegStuff % w(t)

       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element,Nodes, u, v, w, detJ, &
             Basis, dBasisdx )

       s = detJ * IntegStuff % s(t)

       ! Lame parameters at the integration point:
       ! -----------------------------------------
       Lame1 = SUM( NodalLame1(1:n)*Basis(1:n) )
       Lame2 = SUM( NodalLame2(1:n)*Basis(1:n) )


       ! Loop over basis functions (of both unknowns and weights):
       ! ---------------------------------------------------------
       DO p=1,nd
       DO q=p,nd
          A => STIFF( dim*(p-1)+1:dim*p,dim*(q-1)+1:dim*q )
          DO i=1,dim
             DO j = 1,dim
                A(i,j) = A(i,j) + s * Lame1 * dBasisdx(q,j) * dBasisdx(p,i)
                A(i,i) = A(i,i) + s * Lame2 * dBasisdx(q,j) * dBasisdx(p,j)
                A(i,j) = A(i,j) + s * Lame2 * dBasisdx(q,i) * dBasisdx(p,j)
             END DO
          END DO
       END DO
       END DO
     END DO 

     ! Assign the symmetric block:
     ! ---------------------------
     DO p=1,dim*nd
       DO q=1,p-1
         STIFF(p,q)=STIFF(q,p)
       END DO
     END DO

     IF ( nb == 0 ) THEN
       DO p=MAX(n+1,nd-Element % BDOFs+1),nd
         DO i=1,dim
            j = (p-1)*dim + i
            STIFF( j,: ) = 0.0d0
            STIFF( :,j ) = 0.0d0
            STIFF( j,j ) = 1.0d0
            FORCE( j )   = 0.0d0
         END DO
       END DO
     END IF
!------------------------------------------------------------------------------
 END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE MeshBoundary( STIFF,FORCE,LOAD,NodalAlpha,NodalBeta,Element,n,nd,nb )
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: STIFF(:,:),FORCE(:)
   REAL(KIND=dp) :: NodalAlpha(:,:),NodalBeta(:),LOAD(:,:)
   TYPE(Element_t),POINTER  :: Element
   INTEGER :: n,nd,nb
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: Basis(nd)
   REAL(KIND=dp) :: dBasisdx(nd,3),detJ

   REAL(KIND=dp) :: u,v,w,s
   REAL(KIND=dp) :: Alpha(3),Beta,Normal(3),LoadAtIP(3)

   INTEGER :: i,t,q,p

   LOGICAL :: stat

   TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

   TYPE(Nodes_t)    :: Nodes
   SAVE Nodes
!------------------------------------------------------------------------------

   CALL GetElementNodes( Nodes )

   FORCE = 0.0D0
   STIFF = 0.0D0
!
!  Integration stuff
!
   IntegStuff = GaussPoints( Element )
!
!  Now we start integrating
!
   DO t=1,IntegStuff % n

     u = IntegStuff % u(t)
     v = IntegStuff % v(t)
     w = IntegStuff % w(t)

!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
                 Basis, dBasisdx )

      s = detJ * IntegStuff % s(t)
!------------------------------------------------------------------------------
     LoadAtIP = 0.0D0
     DO i=1,dim
       LoadAtIP(i) = SUM( LOAD(i,1:n)*Basis )
       Alpha(i) = SUM( NodalAlpha(i,1:n)*Basis )
     END DO

     Normal = NormalVector( Element,Nodes,u,v,.TRUE. )
     LoadAtIP = LoadAtIP + SUM( NodalBeta(1:n)*Basis ) * Normal

     DO p=1,nd
       DO q=1,nd
         DO i=1,dim
           STIFF((p-1)*dim+i,(q-1)*dim+i) =  &
             STIFF((p-1)*dim+i,(q-1)*dim+i) + &
               s * Alpha(i) * Basis(q) * Basis(p)
         END DO
       END DO
     END DO

     DO q=1,nd
       DO i=1,dim
         FORCE((q-1)*dim+i) = FORCE((q-1)*dim+i) + &
                   s * Basis(q) * LoadAtIP(i)
       END DO
     END DO

   END DO
 
   IF ( nb == 0 ) THEN
     DO p=MAX(n+1,nd-Element % BDOFs+1),nd
       DO i=1,dim
          j = (p-1)*dim + i
          STIFF( j,: ) = 0.0d0
          STIFF( :,j ) = 0.0d0
          STIFF( j,j ) = 1.0d0
          FORCE( j )   = 0.0d0
       END DO
     END DO
   END IF

!------------------------------------------------------------------------------
 END SUBROUTINE MeshBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE MeshSolver
!------------------------------------------------------------------------------
