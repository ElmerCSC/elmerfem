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
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 10 May 2000
! *  Edited:         9.5.2012
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!>  An alternative fork for mesh update solver. The intended use for this
!>  is, for example, in geometry optimization when the normal MeshUpdate is
!>  not always suitable as it may be needed to extent real physical displacements.
!> This is a dynamically loaded solver with a standard interface.
!> \ingroup Solvers
!------------------------------------------------------------------------------
 SUBROUTINE MeshSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  TransientSimulation
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  INTEGER :: i,j,k,n,nd,nb,t,ind,STDOFs,LocalNodes,istat,dim,NoActive
  INTEGER :: VisitedTimes = 0

  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t),POINTER :: Params, Material, BC
  REAL(KIND=dp) :: UNorm, maxu, TargetCoeff, val, Relax
  TYPE(Variable_t), POINTER :: MeshSol, TargetSol, GradSol
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), POINTER :: MeshUpdate(:), PrevMeshUpdate(:)
  INTEGER, POINTER :: MeshPerm(:), NodeIndexes(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: TargetFieldName

  LOGICAL :: AllocationsDone = .FALSE., Isotropic = .TRUE., &
            GotForceBC, Found, MovingMesh, Cumulative,GotTargetField, &
            GotTargetSurface, GotGradSol
			
  REAL(KIND=dp),ALLOCATABLE:: STIFF(:,:),&
       LOAD(:,:),FORCE(:), ElasticModulus(:),PoissonRatio(:), &
		Alpha(:,:), Beta(:), Gamma(:), RefSurface(:)
  REAL(KIND=dp), POINTER CONTIG :: OrigX(:), OrigY(:), OrigZ(:), &
      TrueX(:), TrueY(:), TrueZ(:)

  SAVE STIFF, LOAD, FORCE, AllocationsDone, &
       ElasticModulus, PoissonRatio, &
       OrigX, OrigY, OrigZ, TrueX, TrueY, TrueZ, PrevMeshUpdate, &
       VisitedTimes,dim,Alpha,Beta,Gamma,RefSurface

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at,at0
#else
  REAL(KIND=dp) :: at,at0,CPUTime,RealTime
#endif
!------------------------------------------------------------------------------

 CALL Info( 'MeshSolve', '-------------------------------------', Level=4 )
 CALL Info( 'MeshSolve', 'Nonphysical Mesh Solver:', Level=4 )
 CALL Info( 'MeshSolve', '-------------------------------------', Level=4 )
 
 
!------------------------------------------------------------------------------
! Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

  VisitedTimes = VisitedTimes + 1

  dim = CoordinateSystemDimension()
  MeshSol => Solver % Variable
  MeshPerm => MeshSol % Perm
  STDOFs =  MeshSol % DOFs
  MeshUpdate => MeshSol % Values
  Params => GetSolverParams()
  Mesh => Solver % Mesh
  
  LocalNodes = COUNT( MeshPerm > 0 )
  
  IF ( LocalNodes <= 0 ) RETURN
  Cumulative = ListGetLogical( Params,'Cumulative Displacement',Found)
  
!------------------------------------------------------------------------------
! Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone  ) THEN
    N = Mesh % MaxElementDOFs
    
    IF ( AllocationsDone ) THEN
      DEALLOCATE(  ElasticModulus, PoissonRatio, &
          FORCE, STIFF, Alpha, Beta, Gamma, RefSurface, LOAD, STAT=istat )
    END IF
    
    ALLOCATE( &
        Alpha(3,N), Beta(N), Gamma(N), RefSurface(N), &
        ElasticModulus( N ), PoissonRatio( N ), &
        FORCE( STDOFs*N ), STIFF( STDOFs*N,STDOFs*N ),  &
        LOAD( 4,N ),STAT=istat )
    
    IF(.NOT. Cumulative) THEN
      n = SIZE( Mesh % Nodes % x )
      ALLOCATE( OrigX(n), OrigY(n), OrigZ(n) )
      OrigX = Mesh % Nodes % x
      OrigY = Mesh % Nodes % y
      OrigZ = Mesh % Nodes % z
    END IF
    
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'MeshSolve', 'Memory allocation error.' )
    END IF

    AllocationsDone = .TRUE.
!------------------------------------------------------------------------------
  END IF
!------------------------------------------------------------------------------

  TrueX => Mesh % Nodes % x
  TrueY => Mesh % Nodes % y
  TrueZ => Mesh % Nodes % z
  
  IF( .NOT. Cumulative ) THEN
    Mesh % Nodes % x => OrigX
    Mesh % Nodes % y => OrigY
    Mesh % Nodes % z => OrigZ
  END IF
  
! This refers to an another mesh deformation routine moving the mesh 
! such that these two must coexist. For generality this is set true
! by default.
!--------------------------------------------------------------------
  MovingMesh = ListGetLogical( Params,'Moving Mesh', Found ) 
  IF(.NOT. Found) MovingMesh = .TRUE.
  MovingMesh = MovingMesh .AND. ( VisitedTimes > 1 )

! This variable is needed to enable that the mesh is changed externally independent of this solver.
! Then the variations of this solver must always be relative to the previous update.
!---------------------------------------------------------------------------------
  IF( MovingMesh ) THEN
    IF( VisitedTimes == 2 ) THEN
      n = SIZE( MeshUpdate ) 
      ALLOCATE( PrevMeshUpdate(n) )
      PrevMeshUpdate = 0.0_dp
    END IF
    PrevMeshUpdate = MeshUpdate
  END IF
  
!------------------------------------------------------------------------------
! Do some additional initialization, and go for it
!------------------------------------------------------------------------------
  at  = CPUTime()
  at0 = RealTime()

  CALL DefaultInitialize()
  CALL StartAdvanceOutput( 'NonphysicalMeshSolve', 'Assembly: ' )
!------------------------------------------------------------------------------

  NoActive = GetNOFActive()
  
  
  DO t=1,NoActive
    
    CALL AdvanceOutput( t, NoActive )
    
    Element => GetActiveElement(t)
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()
    n  = GetElementNOFNodes()
    
    Material => GetMaterial()
    
    ElasticModulus(1:n) = GetReal( Material,'Mesh Elastic Modulus', Found )
    IF ( .NOT. Found ) THEN
      ElasticModulus(1:n) = GetReal( Material,'Youngs Modulus', Found )
    END IF
    IF ( .NOT. Found ) ElasticModulus(1:n) = 1.0_dp
    
    PoissonRatio(1:n) = GetReal( Material,'Mesh Poisson Ratio', Found )
    IF ( .NOT. Found ) THEN
      PoissonRatio(1:n) = GetReal( Material,'Poisson Ratio', Found )
    END IF
    IF ( .NOT. Found ) PoissonRatio(1:n) = 0.25_dp
     
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
! Nodal displacements may be set by penalty also. This is the target field
! for the displacements.
!------------------------------------------------------------------------------

  TargetFieldName = GetString(Params,'Target Field',GotTargetField )
  IF( GotTargetField ) THEN
    TargetSol => VariableGet( Model % Variables, TargetFieldName )     
    IF( .NOT. ASSOCIATED(TargetSol)) THEN
      CALL Fatal('NonphysicalMeshSolve',&
          'Given > Target Field < does not exist: '//TRIM( TargetFieldName ) )
    END IF
  END IF
  
  TargetFieldName = GetString(Params,'Target Surface',GotTargetSurface )
  GotGradSol = .FALSE.
  IF( GotTargetSurface ) THEN
    IF( GotTargetField ) THEN
      CALL Fatal('NonphysicalMeshSolve','Cannot have > Target Field < and > Target Surface < at same time!')
    END IF
    TargetSol => VariableGet( Model % Variables, TargetFieldName )     
    IF( .NOT. ASSOCIATED(TargetSol)) THEN
      CALL Fatal('NonphysicalMeshSolve',&
          'Given > Target Surface < does not exist: '//TRIM( TargetFieldName ) )
    END IF
    
    TargetFieldName = GetString(Params,'Grad Surface',GotGradSol ) 
    IF( GotGradSol ) THEN
      GradSol => VariableGet( Model % Variables, TargetFieldName )     
      IF( .NOT. ASSOCIATED(TargetSol)) THEN
        CALL Fatal('NonphysicalMeshSolve',&
            'Given > Grad Surface < does not exist: '//TRIM( TargetFieldName ) )
      END IF
    END IF
  ELSE
    RefSurface = 0.0_dp
  END IF

  
!------------------------------------------------------------------------------
!     Neumann & Newton boundary conditions
!------------------------------------------------------------------------------
  DO t = 1, Mesh % NumberOfBoundaryElements
    
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
    
    NodeIndexes => Element % NodeIndexes
    
    LOAD = 0.0_dp
    Alpha = 0.0_dp
    Beta = 0.0_dp
    Gamma = 0.0_dp    
    RefSurface = 0.0_dp

    Alpha(1,1:n) =  GetReal( BC, 'Mesh Coefficient 1', Found )
    GotForceBC = Found
    Alpha(2,1:n) =  GetReal( BC, 'Mesh Coefficient 2', Found )
    GotForceBC = GotForceBC .OR. Found
    Alpha(3,1:n) =  GetReal( BC, 'Mesh Coefficient 3', Found )
    GotForceBC = GotForceBC .OR. Found
    
    LOAD(1,1:n) =  GetReal( BC, 'Mesh Force 1', Found )
    GotForceBC = GotForceBC .OR. Found
    LOAD(2,1:n) =  GetReal( BC, 'Mesh Force 2', Found )
    GotForceBC = GotForceBC .OR. Found
    LOAD(3,1:n) =  GetReal( BC, 'Mesh Force 3', Found )
    GotForceBC = GotForceBC .OR. Found
    
    Beta(1:n) = GetReal( BC, 'Mesh Normal Force',Found )
    GotForceBC = GotForceBC .OR. Found
    
    Gamma(1:n) =  GetReal( BC, 'Mesh Penalty Factor', Found )
    GotForceBC = GotForceBC .OR. Found

    IF( GotTargetSurface ) THEN
      RefSurface(1:n) = GetReal( BC,'Reference Surface',Found)
    END IF

    IF ( .NOT. GotForceBC ) CYCLE
     
    CALL MeshBoundary( STIFF,FORCE,LOAD,Alpha,Beta,Gamma,RefSurface,Element,n,nd,nb )
    
!------------------------------------------------------------------------------

    CALL DefaultUpdateEquations( STIFF, FORCE )
   END DO
!------------------------------------------------------------------------------

   CALL DefaultFinishAssembly()
   CALL Info( 'MeshSolve', 'Assembly done', Level=4 )
   
!------------------------------------------------------------------------------
! Set the nodal displacement for the nodes by using weights, if requested
!------------------------------------------------------------------------------
   CALL NodalDisplacementPenalty()
      
!------------------------------------------------------------------------------
! Dirichlet boundary conditions
!------------------------------------------------------------------------------
   CALL DefaultDirichletBCs()
   
!------------------------------------------------------------------------------
   CALL Info( 'MeshSolve', 'Set boundaries done', Level=4 )
!------------------------------------------------------------------------------
! Solve the system and check for convergence
!------------------------------------------------------------------------------

   UNorm = DefaultSolve()
      
   Relax = ListGetCReal( Params,'Nonlinear System Relaxation Factor',Found)
   IF( Found ) MeshUpdate = Relax * MeshUpdate

   n = SIZE( MeshPerm )
   Relax = ListGetCReal( Params,'Mesh Relaxation Factor',Found)
   IF( .NOT. Found ) Relax = 1.0_dp
   
   IF( MovingMesh ) THEN
     DO i=1,n
       j = MeshPerm(i)
       IF( j == 0 ) CYCLE
       IF( dim == 2 ) THEN
         Truex(i) = Truex(i) + Relax * ( MeshUpdate(2*j-1) - PrevMeshUpdate(2*j-1) )
         Truey(i) = Truey(i) + Relax * ( MeshUpdate(2*j) - PrevMeshUpdate(2*j) )
       ELSE
         Truex(i) = Truex(i) + Relax * ( MeshUpdate(3*j-2) - PrevMeshUpdate(3*j-2) )
         Truey(i) = Truey(i) + Relax * ( MeshUpdate(3*j-1) - PrevMeshUpdate(3*j-1) )
         Truez(i) = Truez(i) + Relax * ( MeshUpdate(3*j) - PrevMeshUpdate(3*j) )
       END IF
     END DO
   ELSE
     DO i=1,n
       j = MeshPerm(i)
       IF( j == 0 ) CYCLE
       IF( dim == 2 ) THEN
         Truex(i) = Truex(i) + Relax * MeshUpdate(2*j-1)
         Truey(i) = Truey(i) + Relax * MeshUpdate(2*j) 
       ELSE
         Truex(i) = Truex(i) + Relax * MeshUpdate(3*j-2)
         Truey(i) = Truey(i) + Relax * MeshUpdate(3*j-1)
         Truez(i) = Truez(i) + Relax * MeshUpdate(3*j)
       END IF
     END DO
   END IF
   
   IF(.NOT. Cumulative) THEN
     Mesh % Nodes % x => TrueX
     Mesh % Nodes % y => TrueY
     Mesh % Nodes % z => TrueZ
   END IF
   
   

  CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE LocalMatrix( STIFF,FORCE,NodalYoung, NodalPoisson, &
              PlaneStress, Isotropic, Element,n, nd, nb )
!------------------------------------------------------------------------------
     IMPLICIT NONE

     REAL(KIND=dp) :: NodalPoisson(:), NodalYoung(:)
     REAL(KIND=dp), TARGET :: STIFF(:,:), FORCE(:)
     INTEGER :: n,nd,nb
     TYPE(Element_t) :: Element
     LOGICAL :: PlaneStress, Isotropic
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(nd)
     REAL(KIND=dp) :: dBasisdx(nd,3),detJ

     REAL(KIND=dp) :: NodalLame1(n),NodalLame2(n),Lame1,Lame2, &
                      Poisson, Young, Coeff

     REAL(KIND=dp), POINTER :: A(:,:)
     REAL(KIND=dp) :: s,u,v,w
     INTEGER :: i,j,k,p,q,t,dim  
     LOGICAL :: stat
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     TYPE(Nodes_t) :: Nodes
     SAVE  Nodes
!------------------------------------------------------------------------------

     CALL GetElementNodes( Nodes )
     dim = CoordinateSystemDimension()

     Coeff = ListGetConstReal( Params,'Mass Coefficient',Stat)
     
     IF ( PlaneStress ) THEN
        NodalLame1(1:n) = NodalYoung(1:n) * NodalPoisson(1:n) / &
               ((1.0d0 - NodalPoisson(1:n)**2))
     ELSE
        NodalLame1(1:n) = NodalYoung(1:n) * NodalPoisson(1:n) /  &
           ((1.0d0 + NodalPoisson(1:n)) * (1.0d0 - 2.0d0*NodalPoisson(1:n)))
     END IF

     NodalLame2(1:n) = NodalYoung(1:n) / (2* (1.0d0 + NodalPoisson(1:n)))

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
             A(i,i) = A(i,i) + s * Coeff * dBasisdx(q,i) * dBasisdx(p,i)
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
       DO p=nd-Element % BDOFs+1,nd
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
 SUBROUTINE MeshBoundary( STIFF,FORCE,LOAD,NodalAlpha,NodalBeta,NodalGamma,&
	NodalRefSurface,Element,n,nd,nb )
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: STIFF(:,:),FORCE(:)
   REAL(KIND=dp) :: NodalAlpha(:,:),NodalBeta(:),NodalGamma(:),NodalRefSurface(:),LOAD(:,:)
   TYPE(Element_t),POINTER  :: Element
   INTEGER :: n,nd,nb
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: Basis(nd), dBasisdx(nd,3),detJ
   REAL(KIND=dp) :: u,v,w,s
   REAL(KIND=dp) :: Alpha(3),Beta,Gamma,Normal(3),LoadAtIP(3),ExtDisp(3),RefSurface, &
       GradSurface(3), AbsGradSurface2, Surface
   REAL(KIND=dp), ALLOCATABLE :: ParentBasis(:), dParentBasisdx(:,:)
   REAL(KIND=dp) :: ParentdetJ

   INTEGER :: i,t,q,p,dim,np
   LOGICAL :: stat, AllocationsDone
   TYPE(Element_t),POINTER :: Parent

   TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

   TYPE(Nodes_t)    :: Nodes, ParentNodes
   SAVE Nodes, ParentNodes, AllocationsDone, ParentBasis, dParentBasisdx
!------------------------------------------------------------------------------

   dim = Element % TYPE % DIMENSION + 1
   CALL GetElementNodes( Nodes )

   FORCE = 0.0D0
   STIFF = 0.0D0

   IntegStuff = GaussPoints( Element )

   IF( GotTargetSurface .AND. .NOT. GotGradSol ) THEN
     IF( .NOT. AllocationsDone ) THEN
       np = Mesh % MaxElementNodes
       ALLOCATE( ParentBasis(np), dParentBasisdx(np,3))
       AllocationsDone = .TRUE.
     END IF

     Parent => Element % BoundaryInfo % Left
     IF ( .NOT. ASSOCIATED(Parent) ) &
         Parent => Element % BoundaryInfo % Right
     IF(.NOT.ASSOCIATED(Parent)) RETURN

     np = GetElementNOFDOFs(Parent)
     CALL GetElementNodes( ParentNodes, Parent )
   END IF


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
     LoadAtIP = LoadAtIP + SUM( NodalBeta(1:n) * Basis ) * Normal

     IF( GotTargetSurface .OR. GotTargetField ) THEN
       IF( GotTargetField ) THEN
         DO i=1,dim
           ExtDisp(i) = SUM( Basis(1:n) * TargetSol % Values( &
               dim * ( TargetSol % Perm( Element % NodeIndexes(1:n) ) - 1) +  i ) )
         END DO
       ELSE
         RefSurface = SUM( Basis(1:n) * NodalRefSurface(1:n))
         Surface = SUM( Basis(1:n) * TargetSol % Values( Element % NodeIndexes(1:n) ) )

         IF( GotGradSol ) THEN
           DO i=1,dim
             GradSurface(i) = &
                 SUM( Basis(1:n) * GradSol % Values( dim*(Element % NodeIndexes(1:n)-1)+i ) )
           END DO
         ELSE
           CALL GetParentUVW( Element,n,Parent,np,U,V,W,Basis )
           stat = ElementInfo( Parent,ParentNodes,U,V,W,ParentdetJ, &
               ParentBasis,dParentBasisdx )
           DO i=1,dim
             GradSurface(i) = &
                 SUM( dParentBasisdx(1:np,i) * TargetSol % Values( Parent % NodeIndexes(1:np) ) )
           END DO
         END IF
         ! A suggested displacement to the direction of the normal. 
         ! The absolute gradient is used both in the normalization and 
         ! when making the unit vector and hence the SQRT operator is not 
         ! performed on purpose.
         AbsGradSurface2 = SUM( GradSurface(1:dim)**2) 
         DO i=1,dim
           ExtDisp(i) = (RefSurface-Surface) * GradSurface(i) / AbsGradSurface2
         END DO
       END IF
       Gamma = SUM( NodalGamma(1:n) * Basis )
       LoadAtIP = LoadAtIP + ExtDisp * Gamma
     ELSE        
       Gamma = 0.0_dp
     END IF


     DO p=1,nd
       DO q=1,nd
         DO i=1,dim
           STIFF((p-1)*dim+i,(q-1)*dim+i) =  &
             STIFF((p-1)*dim+i,(q-1)*dim+i) + &
               s * ( Alpha(i) + Gamma ) * Basis(q) * Basis(p)
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
     DO p=nd-Element % BDOFs+1,nd
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

!-----------------------------------------------------------   
! Set the nodal coordinates by penalty
!------------------------------------------------------------
 SUBROUTINE NodalDisplacementPenalty()

   INTEGER :: i,j,k,j2
   REAL(KIND=dp) :: ExtDisp(3),GradSurface(3),AbsGradSurface2,Surface,&
       RefSurface,TargetCoeff,Coeff
   LOGICAL :: LocalPenalty, GotWeightSol,Visited=.FALSE.
   INTEGER, POINTER :: PenaltyPerm(:),WeightPerm(:)
   TYPE(Variable_t), POINTER :: WeightSol

   SAVE Visited

   IF(.NOT. (GotTargetSurface .OR. GotTargetField) ) RETURN

   TargetCoeff = GetCReal(Params,'Nodal Penalty Factor',Found)     
   IF( ABS( TargetCoeff ) < TINY(TargetCoeff) ) RETURN

   GotWeightSol = .FALSE.
   IF( GetLogical( Params,'Use Boundary Weights',Found) ) THEN
     IF( .NOT. Visited ) THEN
       CALL CalculateNodalWeights( Solver,.TRUE.,VarName = 'Boundary Weights') 
       WeightSol => VariableGet( Model % Variables,'Boundary Weights')
       WeightSol % Output = .FALSE.
     ELSE
       WeightSol => VariableGet( Model % Variables,'Boundary Weights')
     END IF       
     GotWeightSol = ASSOCIATED(WeightSol)
     IF( GotWeightSol ) THEN
       CALL Info('NonphysicalMeshSolver','Using Boundary Weights for setting nodal penalty',Level=8)
     ELSE
       CALL Fatal('NonphysicalMeshSolver','Boundary Weights not present!')
     END IF
   END IF

   LocalPenalty = ListCheckPresentAnyBC( Model,'Apply Nodal Penalty') 
   IF(.NOT. LocalPenalty ) THEN
     LocalPenalty = ListCheckPresentAnyBodyForce( Model,'Apply Nodal Penalty')     
   END IF

   IF( LocalPenalty ) THEN
     CALL Info('NonphysicalMeshSolve','Setting permutation for local penalties')
     ALLOCATE( PenaltyPerm( Mesh % NumberOfNodes ) )
     CALL MakePermUsingMask( Model, Solver, Mesh,'Apply Nodal Penalty', .FALSE., PenaltyPerm, i )
     ! PRINT *,'Number of penalty nodes',i
   END IF
   
   IF( .NOT. GotTargetField ) THEN
     RefSurface = ListGetCReal( Params,'Reference Surface')
   END IF

   DO i=1,Mesh % NumberOfNodes
         
     j = MeshPerm(i)
     IF(j==0) CYCLE
     
     j2 = TargetSol % Perm(i)
     IF(j2==0) CYCLE

     IF(LocalPenalty) THEN
       IF( PenaltyPerm(i) == 0 ) CYCLE
     END IF

     IF( GotTargetField ) THEN
       DO k=1,dim
         ExtDisp(k) = TargetSol % Values( dim * (j2-1) + k )
       END DO
     ELSE
       Surface = TargetSol % Values( j2 )
       IF( GotGradSol ) THEN
         j2 = GradSol % Perm(i)
         DO k=1,dim
           GradSurface(k) = GradSol % Values( dim * (j2-1) + k )  
         END DO
       ELSE
         CALL Fatal('NonphysicalMeshSolve','You should provide GradSol!')
       END IF
       AbsGradSurface2 = SUM( GradSurface(1:dim)**2) 
       DO k=1,dim
         ExtDisp(k) = (RefSurface-Surface) * GradSurface(k) / AbsGradSurface2
       END DO
     END IF

     DO k=1,dim
       ind = 3 * (j-1) + k
       Coeff = TargetCoeff
       IF( GotWeightSol ) Coeff = Coeff * WeightSol % Values( WeightSol % Perm(i) )
       Solver % Matrix % rhs(ind) = Solver % Matrix % rhs(ind) + TargetCoeff * ExtDisp(k)
       CALL CRS_AddToMatrixElement( Solver % Matrix,ind,ind,TargetCoeff )
     END DO
   END DO

   IF( LocalPenalty ) DEALLOCATE( PenaltyPerm ) 

   Visited = .TRUE.

 END SUBROUTINE NodalDisplacementPenalty



!------------------------------------------------------------------------------
END SUBROUTINE MeshSolver
!------------------------------------------------------------------------------
