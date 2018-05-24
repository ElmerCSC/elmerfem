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
! *  Authors: Peter R�back
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 24.2.2009
! *
! *****************************************************************************/


!------------------------------------------------------------------------------
!>  Subroutine for mapping the mesh using analytical commands of scaling, 
!> rotation, translation and smoothing. Additionally may include a grading field in [0,1]
!> that may be solved from a Laplace equation. Provides often the most economical way
!> of distorting the mesh.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE RigidMeshMapper( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE CoordinateSystems
  USE ElementUtils
  USE ElementDescription
  USE ParallelUtils
  USE Types
  USE Lists
  USE DefUtils

  
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: SolverParams, BC
  INTEGER :: i,j,k,n,m,t,dim,elem,bf_id,istat,node,nodei,RotateOrder(3)
  INTEGER, TARGET :: CurrentNode(1)
  INTEGER :: NonlinIter, MaxNonlinIter, NoNodes
  INTEGER, POINTER :: RelaxPerm(:), VeloPerm(:), NodeIndex(:), IntArray(:) => NULL()
  REAL(KIND=dp) :: x0(4), x1(4), RotMatrix(4,4),TrsMatrix(4,4),SclMatrix(4,4), &
      TrfMatrix(4,4),Identity(4,4), Origin(4),Angles(3),Scaling(3),alpha, dCoord(3), Norm, dx(3)
  REAL(KIND=dp) :: at0,at1,at2,Coeff,Source,relax(1),MaxDeform
#ifndef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: CPUTime,RealTime
#endif
  REAL(KIND=dp), POINTER :: Xorig(:),Yorig(:),Zorig(:),Xnew(:),Ynew(:),Znew(:),&
      RelaxField(:),VeloVal(:), PArray(:,:) => NULL()
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)
  TYPE(Variable_t), POINTER :: VeloVar, RelaxVar
  TYPE(Mesh_t), POINTER :: Mesh
  LOGICAL :: Found,GotMatrix,GotRotate,GotTranslate,GotScale,Visited=.FALSE.,&
      UseOriginalMesh, Cumulative, GotRelaxField=.FALSE., &
      CalculateVelocity,TranslateBeforeRotate
  LOGICAL :: AnyMeshMatrix,AnyMeshRotate,AnyMeshTranslate,AnyMeshScale,&
      AnyMeshOrigin
  LOGICAL, POINTER :: NodeDone(:)
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t), SAVE :: Nodes
  TYPE(ValueList_t),POINTER :: BodyForce

  
  SAVE Xorig,Yorig,Zorig,Parray,Visited
   
  CALL Info( 'RigidMeshMapper','---------------------------------------',Level=4 )
  CALL Info( 'RigidMeshMapper','Performing analytic mesh mapping ',Level=4 )
  CALL Info( 'RigidMeshMapper','---------------------------------------',Level=4 )

  SolverParams => GetSolverParams()
  Mesh => Solver % Mesh

  Cumulative = GetLogical( SolverParams,'Cumulative Displacements',Found)
  UseOriginalMesh = .NOT. Cumulative

  dim = CoordinateSystemDimension()
  
  Xnew => Mesh % Nodes % x
  Ynew => Mesh % Nodes % y
  Znew => Mesh % Nodes % z

  NoNodes = Mesh % NumberOfNodes
  ALLOCATE( NodeDone(NoNodes) )

  NodeDone = .FALSE.
  NodeIndex => CurrentNode


  CalculateVelocity = GetLogical( SolverParams,'Calculate Mesh Velocity',Found)
  IF( CalculateVelocity ) THEN
    VeloVar => VariableGet( Mesh % Variables,'Mesh Velocity')
    IF(ASSOCIATED(VeloVar)) THEN
      VeloVal => VeloVar % Values
      VeloPerm => VeloVar % Perm
    ELSE
      CALL Info('RigidMeshMapper','Creating variable for mesh velocity')
      n = SIZE( Xnew )
      ALLOCATE( VeloVal(dim*NoNodes), VeloPerm(NoNodes) )
      VeloVal = 0.0_dp
      DO i=1,NoNodes
        VeloPerm(i) = i
      END DO
      IF( dim == 2 ) THEN
        CALL VariableAddVector( Mesh % variables, Mesh, Solver,&
            'Mesh Velocity[Mesh Velocity:2]',&
            dim, VeloVal, VeloPerm ) 
      ELSE
        CALL VariableAddVector( Mesh % variables, Mesh, Solver,&
            'Mesh Velocity[Mesh Velocity:3]',&
            dim, VeloVal, VeloPerm ) 
      END IF
    END IF

    CALL InvalidateVariable( CurrentModel % Meshes, Solver % Mesh,&
        'Mesh Velocity' )	
  END IF

  TranslateBeforeRotate = GetLogical( SolverParams,&
      'Translate Before Rotate',Found )

  ! Permit the user to specify the order of rotation
  ! The order is reversed because the order of matrix 
  ! multiplication is such that current default order is
  ! 3, 2, 1
  !----------------------------------------------------------
  RotateOrder = (/1, 2, 3/)
  IntArray => ListGetIntegerArray( SolverParams,'Mesh Rotation Axis Order', Found)
  IF(Found) THEN
     j = SIZE(IntArray)
     IF(j /= dim) THEN
        CALL Fatal("RigidMeshMapper","Size of Mesh Rotation Axis Order must match dimension.")
     END IF
     DO i=1,j
        RotateOrder(i) = IntArray(j+1-i) !reverse the order
     END DO
  END IF

  ! If using original mesh as a reference mesh it must be saved,
  ! otherwise the analytic mapping does not require two meshes
  !------------------------------------------------------------
  IF( UseOriginalMesh ) THEN
    IF( .NOT. Visited ) THEN
      WRITE(Message,* ) 'Allocating new nodes of size: ',NoNodes
      CALL Info('RigidMeshMapper',Message,Level=6)
      ALLOCATE(Xorig(NoNodes),Yorig(NoNodes),Zorig(NoNodes))      
      Xorig = Xnew(1:NoNodes)
      Yorig = Ynew(1:NoNodes)
      Zorig = Znew(1:NoNodes)
    END IF
  ELSE
    Xorig => Xnew
    Yorig => Ynew
    Zorig => Znew
  END IF


  IF( .NOT. Visited .AND. ASSOCIATED(Solver % Matrix) ) THEN
    N = Solver % Mesh % MaxElementNodes 
    ALLOCATE( FORCE(N), STIFF(N,N), STAT=istat )

    ! Implement moving and fixed BCs
    ! ------------------------------
    DO i=1,Model % NumberOFBCs
      BC => Model % BCs(i) % Values
      IF ( GetLogical(  BC, 'Moving Boundary', Found ) ) THEN
        CALL ListAddConstReal( BC,Solver % Variable % Name, 1.0_dp )
      ELSE IF ( GetLogical(  BC, 'Fixed Boundary', Found ) ) THEN
        CALL ListAddConstReal( BC,Solver % Variable % Name, 0.0_dp )
      END IF
    END DO

    CALL Info('RigidMeshMapper','Solving mesh relaxation field using Laplace',Level=6)
    
    MaxNonlinIter = GetInteger( SolverParams,&
       'Nonlinear System Max Iterations',Found)
    IF(.NOT. Found) MaxNonlinIter = 1
    
    Coeff = GetCReal( SolverParams,'Nonlinear Conductivity Coeffient',Found)
    Source = GetCReal( SolverParams,'Mesh Relax Source',Found)


    
    DO NonlinIter = 1, MaxNonlinIter
      CALL DefaultInitialize()
      
      DO t=1, GetNOFActive()
        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        CALL LocalMatrix(  STIFF, FORCE, Element, n )
        CALL DefaultUpdateEquations( STIFF, FORCE )
      END DO
      CALL DefaultFinishBulkAssembly()

      ! No Flux BCs
      CALL DefaultFinishAssembly()
      CALL DefaultDirichletBCs()
      Norm = DefaultSolve()      
      
      IF( Solver % Variable % NonlinConverged == 1 ) EXIT
    END DO

    IF( ListGetLogical(SolverParams,'Mesh Relax Normalize') ) THEN
      MaxDeform = MAXVAL( ABS( Solver % Variable % Values ) )
      MaxDeform = ParallelReduction( MaxDeform, 2 )      
      WRITE(Message,'(A,ES12.3)') 'Normalizing deformation by:',MaxDeform
      CALL Info('RigidMeshMapper',Message,Level=6)
      Solver % Variable % Values = Solver % Variable % Values / MaxDeform
    END IF
    
    DEALLOCATE( FORCE, STIFF )
  END IF

  RelaxVar => Solver % Variable
  IF( ASSOCIATED( RelaxVar ) ) THEN
    IF( ASSOCIATED( RelaxVar % Values ) ) THEN
      IF( SIZE( RelaxVar % Values ) > 0 ) THEN
        GotRelaxField = .TRUE.
        RelaxField => Solver % Variable % Values
        RelaxPerm => Solver % Variable % Perm
      END IF
    END IF
  END IF

    
  ! Initialize the mapping matrices
  Identity = 0.0d0
  DO i=1,4
    Identity(i,i) = 1.0d0
  END DO
  
  at0 = CPUTime()

   AnyMeshTranslate = ListCheckPresentAnyBodyForce( Model,'Mesh Displacement') .OR. &
      ListCheckPresentAnyBodyForce( Model,'Mesh Displacement 1') .OR. &
      ListCheckPresentAnyBodyForce( Model,'Mesh Displacement 2') .OR. &
      ListCheckPresentAnyBodyForce( Model,'Mesh Displacement 3')
   IF( AnyMeshTranslate ) THEN
     CALL Info('RigidMeshMapper','> Mesh Displacement < is an obsolite keyword',Level=3)
     CALL Warn('RigidMeshMapper','Replace with > Mesh Translate < ')
     AnyMeshTranslate = .FALSE.
   END IF


  AnyMeshMatrix = ListCheckPresentAnyBodyForce( Model,'Mesh Matrix')
  AnyMeshRotate = ListCheckPresentAnyBodyForce( Model,'Mesh Rotate') .OR. &
      ListCheckPresentAnyBodyForce( Model,'Mesh Rotate 1') .OR. &
      ListCheckPresentAnyBodyForce( Model,'Mesh Rotate 2') .OR. &
      ListCheckPresentAnyBodyForce( Model,'Mesh Rotate 3')
  AnyMeshTranslate = ListCheckPresentAnyBodyForce( Model,'Mesh Translate') .OR. &
      ListCheckPresentAnyBodyForce( Model,'Mesh Translate 1') .OR. &
      ListCheckPresentAnyBodyForce( Model,'Mesh Translate 2') .OR. &
      ListCheckPresentAnyBodyForce( Model,'Mesh Translate 3')
  AnyMeshScale = ListCheckPresentAnyBodyForce( Model,'Mesh Scale') .OR. &
      ListCheckPresentAnyBodyForce( Model,'Mesh Scale 1') .OR. &
      ListCheckPresentAnyBodyForce( Model,'Mesh Scale 2') .OR. &
      ListCheckPresentAnyBodyForce( Model,'Mesh Scale 3')
  AnyMeshOrigin = ListCheckPresentAnyBodyForce( Model,'Mesh Origin') .OR. &
      ListCheckPresentAnyBodyForce( Model,'Mesh Origin 1') .OR. &
      ListCheckPresentAnyBodyForce( Model,'Mesh Origin 2') .OR. &
      ListCheckPresentAnyBodyForce( Model,'Mesh Origin 3')

    
  GotRotate = .FALSE.
  GotTranslate = .FALSE.
  GotScale = .FALSE.
  GotMatrix = .FALSE.
  
  DO elem = 1,Solver % Mesh % NumberOfBulkElements      

    Element => Solver % Mesh % Elements(elem)
    Model % CurrentElement => Element
    n = GetElementNOFNodes()
    
    bf_id = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Body Force',Found )
    IF(.NOT. Found) CYCLE
    
    BodyForce => Model % BodyForces(bf_id) % Values
    
    DO Node=1,n
      NodeIndex(1) = Element % NodeIndexes(Node)
      NodeI = NodeIndex(1)
     
      IF(NodeDone(NodeI)) CYCLE
      NodeDone(NodeI) = .TRUE.

      ! Generic transformation matrix
      !--------------------------------
      GotMatrix = .FALSE.
      IF( AnyMeshMatrix ) THEN
        Parray => ListGetConstRealArray( BodyForce,'Mesh Matrix', GotMatrix )
        IF ( GotMatrix ) THEN
          DO i=1,SIZE(Parray,1)
            DO j=1,SIZE(Parray,2)
              TrfMatrix(i,j) = Parray(j,i)
            END DO
          END DO
        END IF
      END IF

      IF(.NOT. GotMatrix ) THEN

        ! Rotations around main axis:
        !----------------------------        
        GotRotate = .FALSE.
        IF( AnyMeshRotate ) THEN
          Parray => ListGetConstRealArray( BodyForce,'Mesh Rotate', GotRotate )                
          IF ( GotRotate ) THEN
            Angles = 0.0_dp
            DO i=1,SIZE(Parray,1)
              Angles(i) = Parray(i,1) 
            END DO
          ELSE 
            Angles = 0.0_dp
            Angles(1:1) = ListGetReal( BodyForce,'Mesh Rotate 1', 1, NodeIndex, Found )
            IF( Found ) GotRotate = .TRUE.
            Angles(2:2) = ListGetReal( BodyForce,'Mesh Rotate 2', 1, NodeIndex, Found )
            IF( Found ) GotRotate = .TRUE.
            Angles(3:3) = ListGetReal( BodyForce,'Mesh Rotate 3', 1, NodeIndex, Found )
            IF( Found ) GotRotate = .TRUE.
          END IF
        END IF

        ! Scaling:
        !---------
        GotScale = .FALSE.
        IF( AnyMeshScale ) THEN
          Parray => ListGetConstRealArray( BodyForce,'Mesh Scale', GotScale )
          IF ( GotScale ) THEN
            Scaling = 0.0_dp
            DO i=1,SIZE(Parray,1)
              Scaling(i) = Parray(i,1)
            END DO
          ELSE 
            Scaling(1:1) = ListGetReal( BodyForce,'Mesh Scale 1',1,NodeIndex,GotScale) 
	    IF(.NOT. GotScale ) Scaling(1) = 1.0_dp
            Scaling(2:2) = ListGetReal( BodyForce,'Mesh Scale 2',1,NodeIndex,Found) 
	    IF(.NOT. Found ) Scaling(2) = 1.0_dp
            GotScale = GotScale .OR. Found
            Scaling(3:3) = ListGetReal( BodyForce,'Mesh Scale 3',1,NodeIndex,Found) 
	    IF(.NOT. Found ) Scaling(3) = 1.0_dp
            GotScale = GotScale .OR. Found
          END IF
        END IF

        ! Translations:
        !---------------
        GotTranslate = .FALSE.
        IF( AnyMeshTranslate ) THEN
          Parray => ListGetConstRealArray( BodyForce,'Mesh Translate', GotTranslate )
          IF ( GotTranslate ) THEN
            dCoord = 0.0_dp
            DO i=1,SIZE(Parray,1)
              dCoord(i) = Parray(i,1)
            END DO
          ELSE 
            dCoord(1:1) = ListGetReal( BodyForce,'Mesh Translate 1', 1, NodeIndex, GotTranslate) 
            dCoord(2:2) = ListGetReal( BodyForce,'Mesh Translate 2', 1, NodeIndex, Found) 
            GotTranslate = GotTranslate .OR. Found
            dCoord(3:3) = ListGetReal( BodyForce,'Mesh Translate 3', 1, NodeIndex, Found) 
            GotTranslate = GotTranslate .OR. Found
          END IF
        END IF

        GotMatrix = GotRotate .OR. GotScale

        IF(GotMatrix) THEN
          TrsMatrix = Identity
          SclMatrix = Identity
        
          ! Origin:
          !---------
          IF( GotRotate ) THEN
            RotMatrix = Identity
            
            DO i=1,3
              j = RotateOrder(i)
              Alpha = Angles(j) * PI / 180.0_dp
              
              IF( ABS(Alpha) < TINY(Alpha) ) CYCLE
              TrfMatrix = Identity
              
              SELECT CASE(j)
              CASE(1)
                TrfMatrix(2,2) =  COS(Alpha)
                TrfMatrix(2,3) = -SIN(Alpha)
                TrfMatrix(3,2) =  SIN(Alpha)
                TrfMatrix(3,3) =  COS(Alpha)
              CASE(2)
                TrfMatrix(1,1) =  COS(Alpha)
                TrfMatrix(1,3) = -SIN(Alpha)
                TrfMatrix(3,1) =  SIN(Alpha)
                TrfMatrix(3,3) =  COS(Alpha)
              CASE(3)
                TrfMatrix(1,1) =  COS(Alpha)
                TrfMatrix(1,2) = -SIN(Alpha)
                TrfMatrix(2,1) =  SIN(Alpha)
                TrfMatrix(2,2) =  COS(Alpha)
              END SELECT              
              RotMatrix = MATMUL( RotMatrix, TrfMatrix )
            END DO
          END IF

          IF( GotTranslate ) THEN
            DO i=1,3
              TrsMatrix(i,4) = dCoord(i)
            END DO
          END IF
                    
          ! It may be easier to first translate the matrix to origin 
          ! and only the do the rotation than vice versa. 
          IF( TranslateBeforeRotate ) THEN
            TrfMatrix = MATMUL( RotMatrix, TrsMatrix )
          ELSE
            TrfMatrix = MATMUL( TrsMatrix, RotMatrix )
          END IF

          IF( GotScale ) THEN
            DO i=1,3
              SclMatrix(i,i) = Scaling(i)
            END DO
            TrsMatrix = TrfMatrix
            TrfMatrix = MATMUL( SClMatrix, TrsMatrix )
          END IF
        END IF
      END IF

      IF( GotMatrix ) THEN        

        x0(1) = Xorig(NodeI)
        x0(2) = Yorig(NodeI)
        x0(3) = Zorig(NodeI)
        x0(4) = 1.0_dp

        ! Get mesh origin
        !----------------------------------------------------
        IF( AnyMeshOrigin ) THEN
          Origin = 0._dp
          Parray => ListGetConstRealArray( BodyForce,'Mesh Origin', Found )
          IF ( Found ) THEN
            DO i=1,SIZE(Parray,1)
              Origin(i) = Parray(i,1)
            END DO
          ELSE
            Origin(1:1) = ListGetReal( BodyForce,'Mesh Origin 1', 1, NodeIndex, Found) 
            Origin(2:2) = ListGetReal( BodyForce,'Mesh Origin 2', 1, NodeIndex, Found) 
            Origin(3:3) = ListGetReal( BodyForce,'Mesh Origin 3', 1, NodeIndex, Found) 
          END IF
          x1 = MATMUL( TrfMatrix, x0 - Origin ) + Origin
        ELSE
          x1 = MATMUL( TrfMatrix, x0 ) 
        END IF
          
        dx(1:3) = x1(1:3) / x1(4) - x0(1:3)

      ELSE IF( GotTranslate ) THEN
        dx(1:3) = dCoord(1:3)
      ELSE
        CYCLE
      END IF

      ! Find the relaxation parameters that may interpolate the displacement between 
      ! moving and fixed walls.
      !------------------------------------------------------------------------------
      IF( GotRelaxField ) THEN
        IF( RelaxPerm(NodeI) == 0 ) THEN
          Found = .FALSE.
        ELSE
          Relax(1) = RelaxField( RelaxPerm( NodeI ) )
          Found = .TRUE.
        END IF
        dx = Relax(1) * dx 
      ELSE
        Relax(1:1) = ListGetReal( BodyForce,'Mesh Relax',1,NodeIndex,Found)
        IF( Found ) dx = Relax(1) * dx
      END IF

      ! Compute mesh velocity if requested
      !---------------------------------------------------------------
      IF( CalculateVelocity ) THEN
        k = NodeI
        IF( ASSOCIATED( VeloPerm) ) k = VeloPerm(NodeI)
        IF( k > 0 ) THEN  
	  IF( dim == 2 ) THEN
	    VeloVal(2*k-1) = ( Xorig(NodeI) + dx(1) - Xnew(NodeI) ) / dt
	    VeloVal(2*k) = ( Yorig(NodeI) + dx(2) - Ynew(NodeI) ) / dt	
          ELSE
	    VeloVal(3*k-2) = ( Xorig(NodeI) + dx(1) - Xnew(NodeI) ) / dt
	    VeloVal(3*k-1) = ( Yorig(NodeI) + dx(2) - Ynew(NodeI) ) / dt	
	    VeloVal(3*k) = ( Zorig(NodeI) + dx(3) - Znew(NodeI) ) / dt		
          END IF
        END IF
      END IF

      Xnew(NodeI) = dx(1) + Xorig(NodeI)
      Ynew(NodeI) = dx(2) + Yorig(NodeI)
      Znew(NodeI) = dx(3) + Zorig(NodeI)

      NodeDone(NodeI) = .TRUE.
    END DO

  END DO

  IF(.NOT. Visited ) THEN
    WRITE(Message,* ) 'Number of nodes mapped: ',COUNT( NodeDone )
    CALL Info('RigidMeshMapper',Message)
    at1 = CPUTime()
    IF( at1-at0 > 0.1_dp ) THEN
      WRITE(Message,* ) 'Coordinate mapping time: ',at1-at0
      CALL Info('RigidMeshMapper',Message)
    END IF  
    CALL Info('RigidMeshMapper','All done' ) 
  END IF

  DEALLOCATE( NodeDone )
  
  Visited = .TRUE.


CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, n )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ,Grad(3),Cond,LocalRelax(n)
    LOGICAL :: Stat
    INTEGER :: i,j,t
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0
    
    CALL GetScalarLocalSolution( LocalRelax )

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t),  detJ, Basis, dBasisdx )
      DO i=1,3
        Grad(i) = SUM( dBasisdx(:,i) * LocalRelax(1:n) )
      END DO

      Cond = 1.0_dp + Coeff * SQRT( SUM( Grad * Grad ) )
      
      ! Laplace operator
      !------------------
      STIFF(1:n,1:n) = STIFF(1:n,1:n) + Cond * IP % s(t) * DetJ * &
          MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )

      FORCE(1:n) = FORCE(1:n) + Source * IP % s(t) * DetJ * Basis(1:n)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE RigidMeshMapper
!------------------------------------------------------------------------------
