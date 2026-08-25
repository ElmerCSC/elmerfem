!/******************************************************************************
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! *
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! *
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write
! *  to the Free Software Foundation, Inc., 51 Franklin Street,
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Peter Raback
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!>  Boundary condition enforcement: Dirichlet, periodic, constraint modes,
!>  normal-tangential, and passive boundary utilities.
!>  Extracted from SolverUtils.
!------------------------------------------------------------------------------

MODULE BoundaryConditionUtils

    USE SolverBasics
    IMPLICIT NONE

CONTAINS

!> Search faces between passive / non-passive domains; add to boundary
!> elements with given bc-id.
!------------------------------------------------------------------------------
  SUBROUTINE GetPassiveBoundary(Model,Mesh,BcId)
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    INTEGER :: BcId
    TYPE(Mesh_t) :: Mesh 

    INTEGER, ALLOCATABLE :: arr(:)
    INTEGER :: i,j,n,cnt,ind, sz
    LOGICAL :: L1,L2
    TYPE(Element_t), POINTER :: Faces(:), Telems(:), Face, P1, P2

    CALL FindMeshEdges(Mesh,.FALSE.)
    SELECT CASE(Mesh % MeshDim)
    CASE(2)
      Faces => Mesh % Edges
      n = Mesh % NumberOfEdges
    CASE(3)
      Faces => Mesh % Faces
      n = Mesh % NumberOfFaces
    END SELECT

    ALLOCATE(arr(n)); cnt=0
    DO i=1,n
      P1 => Faces(i) % BoundaryInfo % Right
      P2 => Faces(i) % BoundaryInfo % Left
      IF ( .NOT. ASSOCIATED(P1) .OR. .NOT. ASSOCIATED(P2) ) CYCLE

      L1 = CheckPassiveElement(P1)
      L2 = CheckPassiveElement(P2)

      IF ( L1.NEQV.L2) THEN
        cnt = cnt+1
        arr(cnt) = i
      END IF
    END DO

    sz = Mesh % NumberOfBulkElements + Mesh % NumberOFBoundaryElements - &
             Mesh % PassBCcnt
    IF ( sz+cnt>SIZE(Mesh % Elements) ) THEN
      Telems => Mesh % Elements
      ALLOCATE(Mesh % Elements(sz+cnt))
      IF ( ASSOCIATED(Model % Elements,Telems) ) &
        Model % Elements => Mesh % Elements

      Mesh % Elements(1:sz) = Telems

      ! fix boundary element parent pointers to use new array ...
      ! --------------------------------------------------------
      DO i=1,Mesh % NumberOfBoundaryElements-Mesh % PassBCcnt
        ind = i+Mesh % NumberOfBulkElements
        Face => Mesh % Elements(ind)
        IF ( ASSOCIATED(Face % BoundaryInfo % Left) ) &
          Face % BoundaryInfo % Left  => &
             Mesh % Elements(Face % BoundaryInfo % Left % ElementIndex)
        IF ( ASSOCIATED(Face % BoundaryInfo % Right ) ) &
          Face % BoundaryInfo % Right => &
             Mesh % Elements(Face % BoundaryInfo % Right % ElementIndex)
      END DO

      ! ...likewise for  faces (edges).
      ! -------------------------------
      DO i=1,n
        Face => Faces(i)
        IF ( ASSOCIATED(Face % BoundaryInfo % Left) ) &
          Face % BoundaryInfo % Left  => &
             Mesh % Elements(Face % BoundaryInfo % Left % ElementIndex)
        IF ( ASSOCIATED(Face % BoundaryInfo % Right ) ) &
          Face % BoundaryInfo % Right => &
             Mesh % Elements(Face % BoundaryInfo % Right % ElementIndex)
      END DO

      DEALLOCATE(Telems)
    END IF

    DO i=1,cnt
      sz = sz+1
      Mesh % Elements(sz) = Faces(arr(i))
      Mesh % Elements(sz) % Copy = .TRUE.
      Mesh % Elements(sz) % ElementIndex = sz
      Mesh % Elements(sz) % BoundaryInfo % Constraint = BcId
    END DO
    Mesh % NumberOfBoundaryElements = Mesh % NumberOfBoundaryElements - &
                Mesh % PassBCcnt + cnt
    Mesh % PassBCcnt = cnt
    IF ( ASSOCIATED(Model % Elements,Mesh % Elements) ) &
      Model % NumberOfBoundaryElements = Mesh % NumberOfBoundaryElements
!------------------------------------------------------------------------------
  END SUBROUTINE GetPassiveBoundary

!> Sets one Dirichlet condition to the desired value
!------------------------------------------------------------------------------
   SUBROUTINE UpdateDirichletDof( A, dof, dval )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    INTEGER :: dof
    REAL(KIND=dp) :: dval

    IF(.NOT. ALLOCATED(A % ConstrainedDOF)) THEN
      ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
      A % ConstrainedDOF = .FALSE.
    END IF
    
    IF(.NOT. ALLOCATED(A % Dvalues)) THEN
      ALLOCATE(A % Dvalues(A % NumberOfRows))
      A % Dvalues = 0._dp
    END IF
    
    A % Dvalues( dof ) = dval
    A % ConstrainedDOF( dof ) = .TRUE.
    
  END SUBROUTINE UpdateDirichletDof
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE UpdateDirichletDofC( A, dof, cval )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    INTEGER :: dof
    COMPLEX(KIND=dp) :: cval

    IF(.NOT. ALLOCATED(A % ConstrainedDOF)) THEN
      ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
      A % ConstrainedDOF = .FALSE.
    END IF
    
    IF(.NOT. ALLOCATED(A % Dvalues)) THEN
      ALLOCATE(A % Dvalues(A % NumberOfRows))
      A % Dvalues = 0._dp
    END IF
    
    A % Dvalues( 2*dof-1 ) = REAL( cval )
    A % ConstrainedDOF( 2*dof-1 ) = .TRUE.

    A % Dvalues( 2*dof ) = AIMAG( cval )
    A % ConstrainedDOF( 2*dof ) = .TRUE.
    
  END SUBROUTINE UpdateDirichletDofC
!------------------------------------------------------------------------------



  
!> Releases one Dirichlet condition 
!------------------------------------------------------------------------------
   SUBROUTINE ReleaseDirichletDof( A, dof )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    INTEGER :: dof
    REAL(KIND=dp) :: dval
      
    IF(.NOT. ALLOCATED(A % ConstrainedDOF)) THEN
      ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
      A % ConstrainedDOF = .FALSE.
    END IF
    
    IF(.NOT.ALLOCATED(A % Dvalues)) THEN
      ALLOCATE(A % Dvalues(A % NumberOfRows))
      A % Dvalues = 0._dp
    END IF
    
    A % ConstrainedDOF( dof ) = .FALSE.
    
  END SUBROUTINE ReleaseDirichletDof
!------------------------------------------------------------------------------


  
!> Release the range or min/max values of Dirichlet values.
!------------------------------------------------------------------------------
  FUNCTION DirichletDofsRange( Solver, Oper ) RESULT ( val ) 
!------------------------------------------------------------------------------
    TYPE(Solver_t), OPTIONAL :: Solver
    CHARACTER(LEN=*), OPTIONAL :: Oper 
    REAL(KIND=dp) :: val
    
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) :: minv,maxv
    LOGICAL :: FindMin, FindMax
    INTEGER :: i,OperNo
    
    IF( PRESENT( Solver ) ) THEN
      A => Solver % Matrix
    ELSE
      A => CurrentModel % Solver % Matrix
    END IF
    
    val = 0.0_dp
    
    ! Defaulting to range
    OperNo = 0

    IF( PRESENT( Oper ) ) THEN
      IF( Oper == 'range' ) THEN
        OperNo = 0
      ELSE IF( Oper == 'min' ) THEN
        OperNo = 1 
      ELSE IF( Oper == 'max' ) THEN
        OperNo = 2
      ELSE
        CALL Fatal('DirichletDofsRange','Unknown operator: '//TRIM(Oper))
      END IF
    END IF
          
    IF(.NOT. ALLOCATED(A % ConstrainedDOF)) THEN
      RETURN
    END IF
  
    IF( OperNo == 0 .OR. OperNo == 1 ) THEN
      minv = HUGE( minv ) 
      DO i=1,SIZE( A % ConstrainedDOF )
        IF( A % ConstrainedDOF(i) ) minv = MIN( A % DValues(i), minv ) 
      END DO
      minv = ParallelReduction( minv, 1 ) 
    END IF

    IF( OperNo == 0 .OR. OperNo == 2 ) THEN
      maxv = -HUGE( maxv ) 
      DO i=1,SIZE( A % ConstrainedDOF )
        IF( A % ConstrainedDOF(i) ) maxv = MAX( A % DValues(i), maxv ) 
      END DO
      maxv = ParallelReduction( maxv, 2 ) 
    END IF
    
    IF( OperNo == 0 ) THEN    
      val = maxv - minv
    ELSE IF( OperNo == 1 ) THEN
      val = minv
    ELSE
      val = maxv
    END IF
      
  END FUNCTION DirichletDofsRange

!------------------------------------------------------------------------------
!> Set Dirichlet boundary condition for given dof. The conditions are
!> set based on the given name and applied directly to the matrix structure
!> so that a row is zeroed except for the diagonal which is set to one. 
!> Then the r.h.s. value determines the value of the field variable 
!> in the solution of the linear system.
!------------------------------------------------------------------------------
   SUBROUTINE SetDirichletBoundaries( Model, A, b, Name, DOF, NDOFs, Perm, &
       PermOffSet, OffDiagonalMatrix )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model          !< The current model structure
    TYPE(Matrix_t), TARGET :: A     !< The global matrix
    REAL(KIND=dp) :: b(:)           !< The global RHS vector
    CHARACTER(LEN=*) :: Name        !< Name of the dof to be set
    INTEGER :: DOF                  !< The order number of the dof
    INTEGER :: NDOFs                !< The total number of DOFs for this equation
    INTEGER :: Perm(:)              !< The node reordering info, this has been generated at the beginning of the 
                                    !< simulation for bandwidth optimization
    INTEGER, OPTIONAL :: PermOffSet  !< If the matrix and permutation vectors are not in sync the offset may used as a remedy. 
                                     !< Needed in fully coupled systems.
    LOGICAL, OPTIONAL :: OffDiagonalMatrix  !< For block systems the only the diagonal matrix should be given non-zero 
                                            !< entries for matrix and r.h.s., for off-diagonal matrices just set the row to zero.
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:), IndNodes(:), BCOrder(:)
    INTEGER, POINTER :: Indexes(:), PassPerm(:)
    INTEGER :: BC,i,j,j2,k,l,m,n,nd,p,t,k1,k2,OffSet
    LOGICAL :: GotIt, periodic, OrderByBCNumbering, ReorderBCs
    REAL(KIND=dp), POINTER :: WorkA(:,:,:) => NULL()
    REAL(KIND=dp) ::  s

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), POINTER :: Solver

    LOGICAL :: Conditional
    LOGICAL, ALLOCATABLE :: DonePeriodic(:)
    CHARACTER(:), ALLOCATABLE :: CondName, DirName, PassName, PassCondName, EqName

    INTEGER :: NoNodes,NoDims,bf_id,nlen, NOFNodesFound, dim, &
        bndry_start, bndry_end, Upper
    REAL(KIND=dp), POINTER :: CoordNodes(:,:), Condition(:), Work(:)
    REAL(KIND=dp) :: MinDist,Dist, Eps
    LOGICAL, ALLOCATABLE :: ActivePart(:), ActiveCond(:), ActivePartAll(:)
    TYPE(ValueList_t), POINTER :: ValueList, Params
    LOGICAL :: NodesFound, Passive, PassiveCond, OffDiagonal, ApplyLimiter
    LOGICAL, POINTER :: LimitActive(:)
    TYPE(Variable_t), POINTER :: Var

    TYPE(Element_t), POINTER :: Parent

    INTEGER :: ind, ElemFirst, ElemLast, bf, BCstr, BCend, BCinc, dgind
    REAL(KIND=dp) :: SingleVal
    LOGICAL :: AnySingleBC, AnySingleBF
    LOGICAL, ALLOCATABLE :: LumpedNodeSet(:)
    LOGICAL :: NeedListMatrix
    INTEGER :: DirCount
    CHARACTER(*), PARAMETER :: Caller = 'SetDirichletBoundaries'
    LOGICAL, ALLOCATABLE :: CandNodes(:)
    INTEGER, ALLOCATABLE :: LumpedIndx(:)
    INTEGER, POINTER :: PlaneInds(:)
    LOGICAL :: Parallel

    TYPE(NormalTangential_t), POINTER :: NT
    LOGICAL, ALLOCATABLE, SAVE :: NTzeroing_done(:,:)
    INTEGER, ALLOCATABLE, SAVE :: NTelement(:,:)

    REAL(KIND=dp) :: Mult(Model % MaxElementNodes), MaxMult, ParMaxMult, MoveCoeff
    LOGICAL :: GotMult
    INTEGER :: maxind
    
!------------------------------------------------------------------------------
! These logical vectors are used to minimize extra effort in setting up different BCs
!------------------------------------------------------------------------------

    nlen = LEN_TRIM(Name)
    n = MAX( Model % NumberOfBodyForces,Model % NumberOfBCs)
    IF( n == 0 ) THEN
      CALL Info(Caller,'No BCs or Body Forces present, exiting early...',Level=12)
    END IF

    ALLOCATE( ActivePart(n), ActivePartAll(n), ActiveCond(n))
    CondName = Name(1:nlen) // ' Condition'
    PassName = Name(1:nlen) // ' Passive'
    PassCondName = Name(1:nlen) // ' Condition' // ' Passive'

    OffSet = 0
    OffDiagonal = .FALSE.
    IF( PRESENT( PermOffSet) ) OffSet = PermOffSet
    IF( PRESENT( OffDiagonalMatrix ) ) OffDiagonal = OffDiagonalMatrix

    Mesh => Model % Mesh
    ALLOCATE( Indexes(Mesh % MaxElementDOFs) )

    NT => Model % Solver % NormalTangential
    n = NT % NormalTangentialNOFNodes 
    IF( n > 0 ) THEN
      ! We need to have these available for different components of the same vector.
      ! Hence a dirty compromise between localility and saving values. 
      m = 0
      IF( ALLOCATED( NTElement ) ) THEN
        m = SIZE( NTElement, 1 )
      END IF
      IF( m /= n .AND. m > 0 ) THEN
        DEALLOCATE( NTzeroing_done, NTelement )
      END IF
      IF( m /= n ) THEN      
        ALLOCATE( NTzeroing_done(n,3), NTelement(n,3) )
        NTZeroing_done = .FALSE.
        NTelement = 0
      END IF
    END IF
    
    Parallel = ( ParEnv % PEs > 1 ) .AND. ( .NOT. Mesh % SingleMesh ) 
    
!------------------------------------------------------------------------------
! Go through the periodic BCs and set the linear dependence
!------------------------------------------------------------------------------

   ActivePart = .FALSE.
   DO BC=1,Model % NumberOfBCs
     IF ( ListGetLogical( Model % BCs(BC) % Values, &
         'Periodic BC ' // Name(1:nlen), GotIt ) ) ActivePart(BC) = .TRUE.
     IF ( ListGetLogical( Model % BCs(BC) % Values, &
         'Anti Periodic BC ' // Name(1:nlen), GotIt ) ) ActivePart(BC) = .TRUE.
     IF ( ListCheckPresent( Model % BCs(BC) % Values, &
         'Periodic BC Scale ' // Name(1:nlen) ) ) ActivePart(BC) = .TRUE.
   END DO
   
   IF( ANY(ActivePart) ) THEN    
     IF( Offset > 0 ) THEN
       CALL Fatal(Caller,'Periodicity not considered with offset')
     END IF

     ALLOCATE( DonePeriodic( Mesh % NumberOFNodes ) )
     DonePeriodic = .FALSE.
     DO BC=1,Model % NumberOfBCs
       IF( ActivePart(BC) ) THEN
         CALL SetPeriodicBoundariesPass1( Model, A, b, Name, DOF, &
             NDOFs, Perm, BC, DonePeriodic )
       END IF
     END DO
     
     DonePeriodic = .FALSE.
     DO BC=1,Model % NumberOfBCs
       IF(ActivePart(BC)) THEN       
         CALL SetPeriodicBoundariesPass2( Model, A, b, Name, DOF, &
             NDOFs, Perm, BC, DonePeriodic )
       END IF
     END DO

     IF( InfoActive(12) ) THEN
       CALL Info(Caller,'Number of periodic points set: '&
           //I2S(COUNT(DonePeriodic)),Level=12)
     END IF

     DEALLOCATE( DonePeriodic ) 

   END IF
   

! Add the possible friction coefficient
!----------------------------------------------------------
   IF ( ListCheckPresentAnyBC( Model,'Friction BC ' // Name(1:nlen) ) ) THEN
     CALL SetFrictionBoundaries( Model, A, b, Name, NDOFs, Perm )
   END IF


! Add the possible nodal jump in case of mortar projectors
!---------------------------------------------------------------
   IF( ListGetLogical( Model % Solver % Values,'Apply Mortar BCs',GotIt ) ) THEN
     CALL SetWeightedProjectorJump( Model, A, b, &
                      Name, DOF, NDOFs, Perm )
   END IF


!------------------------------------------------------------------------------
! Go through the normal Dirichlet BCs applied on the boundaries
!------------------------------------------------------------------------------

    ActivePart = .FALSE.
    ActiveCond = .FALSE.
    ActivePartAll = .FALSE.
    DO BC=1,Model % NumberOfBCs
      ActivePartAll(BC) = ListCheckPresent( &
            Model % BCs(bc) % Values, Name(1:nlen) // ' DOFs' )
      ActivePart(BC) = ListCheckPresent( Model % BCs(bc) % Values, Name ) 
      ActiveCond(BC) = ListCheckPresent( Model % BCs(bc) % Values, CondName )      

      IF(ActivePart(BC)) THEN
        CALL ListPrepareRealDependence( Model % BCs(bc) % Values, Name ) 
      END IF
      IF(ActiveCond(BC)) THEN
        CALL ListPrepareRealDependence( Model % BCs(bc) % Values, CondName ) 
      END IF
    END DO

    
    OrderByBCNumbering = ListGetLogical( Model % Simulation, &
       'Set Dirichlet BCs by BC Numbering', gotIt)

    BCOrder => ListGetIntegerArray( Model % Solver % Values, &
         'Dirichlet BC Order', ReorderBCs)
    IF(ReorderBCs) THEN
       IF(.NOT. OrderByBCNumbering) THEN
          CALL Warn(Caller,"Requested 'Dirichlet BC Order' but &
               &not 'Set Dirichlet BCs by BC Numbering', ignoring...")
       ELSE IF(SIZE(BCOrder) /= Model % NumberOfBCs) THEN
          CALL Fatal(Caller,"'Dirichlet BC Order' is the wrong length!")
       END IF
    END IF

    bndry_start = Mesh % NumberOfBulkElements+1
    bndry_end   = bndry_start + Mesh % NumberOfBoundaryElements-1
    DirCount = 0

    ! check and set some flags for nodes belonging to n-t boundaries
    ! getting set by other bcs:
    ! --------------------------------------------------------------
    IF ( NT % NormalTangentialNOFNodes>0 ) THEN
      IF ( OrderByBCNumbering ) THEN
        DO i=1,Model % NumberOfBCs
          BC = i
          IF(ReorderBCs) BC = BCOrder(BC)
          IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE
          Conditional = ActiveCond(BC)

          DO t = bndry_start, bndry_end
            Element => Mesh % Elements(t)
            IF ( Element % BoundaryInfo % Constraint /= &
                   Model % BCs(BC) % Tag ) CYCLE

            ValueList => Model % BCs(BC) % Values
            Model % CurrentElement => Element

            IF ( ActivePart(BC) ) THEN
              n = Element % TYPE % NumberOfNodes
              IF ( Model % Solver % DG ) THEN
                 Parent => Element % BoundaryInfo % Left
                 DO p=1,Parent % Type % NumberOfNodes
                   DO j=1,n
                      IF (Parent % NodeIndexes(p) == Element % NodeIndexes(j) ) THEN
                        Indexes(j) = Parent % DGIndexes(p); EXIT
                      END IF
                   END DO
                 END DO
              ELSE
                Indexes(1:n) = Element % NodeIndexes
              END IF
            ELSE
              n = mGetElementDOFs( Indexes )
            END IF
            CALL CheckNTelement(n,t)
          END DO
        END DO
      ELSE
        DO t = bndry_start, bndry_end
          DO BC=1,Model % NumberOfBCs
            IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE
            Conditional = ActiveCond(BC)
          
            Element => Mesh % Elements(t)
            IF ( Element % BoundaryInfo % Constraint /= &
                 Model % BCs(BC) % Tag ) CYCLE
          
            ValueList => Model % BCs(BC) % Values
            Model % CurrentElement => Element
            IF ( ActivePart(BC) ) THEN
              n = Element % TYPE % NumberOfNodes
              IF ( Model % Solver % DG ) THEN
                 Parent => Element % BoundaryInfo % Left
                 DO p=1,Parent % Type % NumberOfNodes
                   DO j=1,n
                      IF (Parent % NodeIndexes(p) == Element % NodeIndexes(j) ) THEN
                        Indexes(j) = Parent % DGIndexes(p); EXIT
                      END IF
                   END DO
                 END DO
              ELSE
                Indexes(1:n) = Element % NodeIndexes
              END IF
            ELSE
              n = mGetElementDOFs( Indexes )
            END IF
            CALL CheckNTelement(n,t)
          END DO
        END DO
      END IF

      IF ( DOF<= 0 ) THEN
        DO t=bndry_start,bndry_end
          Element => Mesh % Elements(t)
          n = Element % TYPE % NumberOfNodes
          DO j=1,n
            k = NT % BoundaryReorder(Element % NodeIndexes(j))
            IF (k>0) THEN
              NTelement(k,:)=0
              NTzeroing_done(k,:) = .FALSE.
            END IF
          END DO
        END DO
      END IF
    END IF

    
    ! Set the Dirichlet BCs from active boundary elements, if any...:
    !----------------------------------------------------------------
    IF( ANY(ActivePart) .OR. ANY(ActivePartAll) ) THEN    
      IF ( OrderByBCNumbering ) THEN
        DO i=1,Model % NumberOfBCs
          BC = i
          IF(ReorderBCs) BC = BCOrder(BC)
          IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE
          Conditional = ActiveCond(BC)

          DO t = bndry_start, bndry_end
            Element => Mesh % Elements(t)
            IF ( Element % BoundaryInfo % Constraint /= &
                Model % BCs(BC) % Tag ) CYCLE
            Model % CurrentElement => Element
            IF ( ActivePart(BC) ) THEN
              n = Element % TYPE % NumberOfNodes
              IF ( Model % Solver % DG ) THEN
                 Parent => Element % BoundaryInfo % Left
                 DO p=1,Parent % Type % NumberOfNodes
                   DO j=1,n
                      IF (Parent % NodeIndexes(p) == Element % NodeIndexes(j) ) THEN
                        Indexes(j) = Parent % DGIndexes(p); EXIT
                      END IF
                   END DO
                 END DO
              ELSE
                Indexes(1:n) = Element % NodeIndexes
              END IF
            ELSE
              n = mGetElementDOFs( Indexes )
            END IF
            ValueList => Model % BCs(BC) % Values
            CALL SetElementValues(n)
          END DO
        END DO
      ELSE
        DO t = bndry_start, bndry_end
          DO BC=1,Model % NumberOfBCs
            IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE
            Conditional = ActiveCond(BC)
            
            Element => Mesh % Elements(t)
            IF ( Element % BoundaryInfo % Constraint /= Model % BCs(BC) % Tag ) CYCLE
            
            Model % CurrentElement => Element
            IF ( ActivePart(BC) ) THEN
              n = Element % TYPE % NumberOfNodes
              IF ( Model % Solver % DG ) THEN
                 Parent => Element % BoundaryInfo % Left
                 DO p=1,Parent % Type % NumberOfNodes
                   DO j=1,n
                      IF (Parent % NodeIndexes(p)  == Element % NodeIndexes(j) ) THEN
                        Indexes(j) = Parent % DGIndexes(p); EXIT
                      END IF
                   END DO
                 END DO
              ELSE
                Indexes(1:n) = Element % NodeIndexes
              END IF
            ELSE
              n = mGetElementDOFs( Indexes )
            END IF
            ValueList => Model % BCs(BC) % Values
            CALL SetElementValues(n)
          END DO
        END DO
      END IF
    END IF

    BLOCK
      INTEGER,ALLOCATABLE :: ChildBCs(:)
      INTEGER,POINTER :: BCInds(:)
      REAL(KIND=dp) :: zero = 0.0_dp
      IF( ListGetLogical( Model % Solver % Values,'Extruded Child BC Zero',GotIt ) ) THEN
        CALL Info(Caller,'Setting extruded BCs (start & end) to zero!',Level=10)

        ALLOCATE(ChildBCs(2*Model % NumberOfBodies))
        ChildBCs=0

        ! Collect a list of child BCs that were generated when extruding
        m = 0
        DO i=1,Model % NumberOfBodies
          BCInds => ListGetIntegerArray(Model % Bodies(i) % Values,'Extruded Child BCs',GotIt)
          IF(GotIt) THEN
            ChildBCs(m+1:m+2) = BCInds(1:2)
            m=m+2
          END IF          
        END DO        
        IF(m==0) CALL Fatal(Caller,'No "Extruded Child BCs" to set')

        IF( InfoActive(20) ) THEN
          PRINT *,'ChildBCs:',m,ChildBCs(1:m)
        END IF
          
        ! Set the extruded BCs to zero. Note that only this value is available currently.
        DO t = bndry_start, bndry_end
          Element => Mesh % Elements(t)
          IF(ANY(ChildBCs(1:m) == Element % BoundaryInfo % Constraint ) ) THEN          
            Model % CurrentElement => Element
            n = mGetElementDOFs( Indexes, Element, Model % Solver )            
            DO i=1,n
              CALL SetSinglePoint(Indexes(i),DOF,zero,.TRUE.)
            END DO
          END IF
        END DO
      END IF
    END BLOCK
      
!------------------------------------------------------------------------------
! Go through the Dirichlet conditions in the body force lists
!------------------------------------------------------------------------------
    
    ActivePart = .FALSE.
    ActiveCond = .FALSE.
    ActivePartAll = .FALSE.
    Passive = .FALSE.
    PassiveCond = .FALSE.
    DO bf_id=1,Model % NumberOFBodyForces
      ValueList => Model % BodyForces(bf_id) % Values

      ActivePartAll(bf_id) = ListCheckPresent(ValueList, Name(1:nlen) // ' DOFs' ) 
      ActiveCond(bf_id) = ListCheckPresent( ValueList,CondName )      
      ActivePart(bf_id) = ListCheckPresent(ValueList, Name(1:nlen) ) 

      Passive = Passive .OR. ListCheckPresent(ValueList, PassName)
      PassiveCond = PassiveCond .OR. ListCheckPresent(ValueList, PassCondName)
    END DO
       
    IF ( ANY(ActivePart) .OR. ANY(ActivePartAll) ) THEN
      Solver => Model % Solver
      Mesh   => Solver % Mesh
      EqName = ListGetString( Solver % Values, 'Equation', GotIt )

      IF( PassiveCond ) THEN
        ALLOCATE(PassPerm(Mesh % NumberOfNodes),NodeIndexes(1))
        PassPerm = 0
        DO i=0,Mesh % PassBCCnt-1
          j=Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements-i
          PassPerm(Mesh % Elements(j) % NodeIndexes)=1
        END DO
      END IF
        
      DO t=1,Solver % Mesh % NumberOfBulkElements
        Element => Mesh % Elements(t)
        IF( Element % BodyId <= 0 .OR. Element % BodyId > Model % NumberOfBodies ) THEN
          CALL Warn(Caller,'Element body id beyond body table!')
          CYCLE
        END IF
                    
        bf_id = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Body Force', GotIt)
        
        IF(.NOT. GotIt) CYCLE
        IF(.NOT. ActivePart(bf_id) .AND. .NOT. ActivePartAll(bf_id)) CYCLE
        Conditional = ActiveCond(bf_id)

        Model % CurrentElement => Element

        GotIt = CheckElementEquation( Model, Element, EqName ) 
        
        n = Element % TYPE % NumberOfNodes
        Indexes(1:n) = Element % NodeIndexes

        ValueList => Model % BodyForces(bf_id) % Values
        IF(.NOT. ASSOCIATED( ValueList ) ) CYCLE
        
        IF (ListGetLogical(ValueList,PassCondName,GotIt)) THEN
          IF (.NOT.CheckPassiveElement(Element)) CYCLE

          IF (ParEnv % PEs > 1) THEN
             IF( Element % PartIndex /= ParEnv % MyPe ) CYCLE
          END IF

          DO j=1,n
            k=Indexes(j)
            IF (k<=0) CYCLE

            k=Perm(k)
            IF (k<=0) CYCLE

            s=0._dp
            DO l=1,NDOFs
              m=NDOFs*(k-1)+l
              s=s+ABS(A % Values(A % Diag(m)))
            END DO
            IF (s>EPSILON(s)) CYCLE

            NodeIndexes(1) = Indexes(j)
            IF(PassPerm(NodeIndexes(1))==0) CALL SetPointValues(1)
          END DO
        ELSE
          CALL SetElementValues(n)
        END IF
        
        ! Set the higher p-dofs related to Dirichlet BC's to zero. 
        IF ( isActivePElement(Element, Solver) ) THEN
          nd = mGetElementDOFs( Indexes, Uelement = Element, USolver = Model % Solver )
          DO i=n+1,nd            
            CALL SetSinglePoint(Indexes(i),DOF,0.0_dp,.TRUE.)            
          END DO
        END IF
        
      END DO
      
      IF(PassiveCond) DEALLOCATE(PassPerm,NodeIndexes)
    END IF
    
    DEALLOCATE(ActivePart, ActiveCond)

    
!------------------------------------------------------------------------------
! Go through the pointwise Dirichlet BCs that are created on-the-fly
! Note that it is best that the coordinates are transformed to nodes using 
! the right variable. Otherwise it could point to nodes that are not active.
!------------------------------------------------------------------------------
     
    DO BC=1,Model % NumberOfBCs
      
      ValueList => Model % BCs(BC) % Values
      IF( .NOT. ListCheckPresent( ValueList,Name )) CYCLE
      NodesFound = ListCheckPresent( ValueList,'Target Nodes' )

      ! The coordinates are only requested for a body that has no list of nodes.
      ! At the first calling the list of coordinates is transformed to list of nodes.
      IF(.NOT. NodesFound) THEN
        IF( ListCheckPresent( ValueList,'Target Coordinates' ) ) THEN
          CALL TargetCoordinatesToTargetNodes( Mesh, ValueList, NodesFound )
        END IF
      END IF
                  
      ! If the target coordinates has already been assigned to an empty list 
      ! cycle over it by testing the 1st node. 
      IF( NodesFound ) THEN
        NodeIndexes => ListGetIntegerArray( ValueList,'Target Nodes')
        IF( NodeIndexes(1) == 0 ) NodesFound = .FALSE. 
      END IF

      IF(NodesFound) THEN           
        Conditional = ListCheckPresent( ValueList, CondName )      
        n = SIZE(NodeIndexes)
        CALL SetPointValues(n)
      END IF
    END DO

    ! Check the boundaries and body forces for possible single nodes BCs that are used to fixed 
    ! the domain for undetermined equations. The loop is slower than optimal in the case that there is 
    ! a large amount of different boundaries that have a node to set. 
    !--------------------------------------------------------------------------------------------
    DirName = TRIM(Name)//' Single Node'
    AnySingleBC = ListCheckPresentAnyBC( Model, DirName )
    AnySingleBF = ListCheckPresentAnyBodyForce( Model, DirName )

    IF( AnySingleBC .OR. AnySingleBF ) THEN
      Solver => Model % Solver
      Mesh   => Solver % Mesh

      DO bc = 1,Model % NumberOfBCs  + Model % NumberOfBodyForces    

        ! Make a distinction between BCs and BFs. 
        ! These are treated in the same loop because most of the logic is still the same. 
        IF( bc <= Model % NumberOfBCs ) THEN
          IF(.NOT. AnySingleBC ) CYCLE
          ValueList => Model % BCs(BC) % Values
          ElemFirst =  Mesh % NumberOfBulkElements + 1 
          ElemLast =  Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
        ELSE
          IF( .NOT. AnySingleBF ) CYCLE
          ValueList => Model % BodyForces(bc - Model % NumberOfBCs) % Values
          ElemFirst =  1
          ElemLast =  Mesh % NumberOfBulkElements
        END IF

        SingleVal = ListGetCReal( ValueList,DirName, GotIt) 
        IF( .NOT. GotIt ) CYCLE
        ind = ListGetInteger( ValueList,TRIM(Name)//' Single Node Index',GotIt )     
        
        ! On the first time find a one single uniquely defined node for setting 
        ! the value. In parallel it will be an unshared node with the highest possible 
        ! node number 
        IF(.NOT. GotIt ) THEN                  
          ind = 0
          DO t = ElemFirst, ElemLast
            Element => Mesh % Elements(t)
            n = Element % TYPE % NumberOfNodes
            NodeIndexes => Element % NodeIndexes

            IF( bc <= Model % NumberOfBCs ) THEN
              IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
            ELSE
              j = Element % BodyId
              IF( j < 0 .OR. j > Model % NumberOfBodies ) CYCLE
              bf = ListGetInteger( Model % Bodies(j) % Values,'Body Force',GotIt)
              IF(.NOT. GotIt) CYCLE
              IF( bc - Model % NumberOfBCs /= bf ) CYCLE
            END IF
            
            DO i=1,n
              j = NodeIndexes(i)
              IF( Perm(j) == 0) CYCLE
              IF( Parallel ) THEN
                IF( SIZE( Mesh % ParallelInfo % NeighbourList(j) % Neighbours) > 1 ) CYCLE               
                IF( Mesh % ParallelInfo % NeighbourList(j) % Neighbours(1) /= ParEnv % MyPe ) CYCLE               
              END IF
              ind = j 
              EXIT
            END DO
            IF( ind > 0 ) EXIT
          END DO

          k = ind
          IF( Parallel ) THEN
            k = ParallelReduction( ind, 2 ) 
            
            ! Find the maximum partition that owns a suitable node. 
            ! It could be minimum also, just some convection is needed. 
            k = -1
            IF( ind > 0 ) k = ParEnv % MyPe          
            k = ParallelReduction( k, 2 ) 
            IF( k == -1 ) THEN
              CALL Warn(Caller,'Could not find node to set: '//TRIM(DirName))
            ELSE
              IF( k /= ParEnv % MyPe ) ind = 0                         
              IF( InfoActive(8) ) THEN
                ind = ParallelReduction( ind, 2 )                
                CALL Info(Caller,'Fixing single node '&
                    //I2S(ind)//' at partition '//I2S(k),Level=8)
                IF( k /= ParEnv % MyPe ) ind = 0
              END IF
            END IF
          ELSE
            IF( ind == 0 ) THEN
              CALL Warn(Caller,'Could not find node to set: '//TRIM(DirName))
            ELSE              
              CALL Info(Caller,'Fixing single node '//I2S(ind),Level=8)
            END IF
          END IF
            
          CALL ListAddInteger( ValueList,TRIM(Name)//' Single Node Index', ind )          
        END IF

        ! Ok, if this is the partition where the single node to eliminate the floating should 
        ! be eliminated then set it here. Index equal to zero tells that we are in a wrong partition.        
        IF( ind > 0 ) THEN
          CALL SetSinglePoint(ind,DOF,SingleVal,.TRUE.)
        END IF
      END DO
    END IF

!------------------------------------------------------------------------------
!   Take care of the matrix entries of passive elements
!------------------------------------------------------------------------------

    IF ( Passive ) THEN
      Solver => Model % Solver
      Mesh => Solver % Mesh

      ALLOCATE(PassPerm(Mesh % NumberOfNodes),NodeIndexes(1))
      PassPerm=0

      ! Mark the interface, don't know what the idea is but it seems to set the
      ! flag to "1" so that we can avoid it when setting Dirichlet conditions. 
      DO i=0,Mesh % PassBCCnt-1
        j=Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements-i
        PassPerm(Mesh % Elements(j) % NodeIndexes) = 1
      END DO

      ! Here set the flag to "2" for all nodes that are active anywhere.
      ! This is basically redundant with the above.
      DO i=1,Solver % NumberOfActiveElements
        Element => Mesh % Elements(Solver % ActiveElements(i))
        IF (.NOT. CheckPassiveElement(Element)) THEN
          PassPerm(Element % NodeIndexes) = 2
        END IF
      END DO

      ! Is is essential to communicate the parallel tag to avoid problems when
      ! passive interface and partition interface match. 
      BLOCK
        TYPE(ParallelInfo_t), POINTER :: ParallelInfo=>NULL()
        ParallelInfo => Mesh % ParallelInfo
        CALL CommunicateParallelSystemTag(ParallelInfo,Itag=PassPerm,ParOper=2)
      END BLOCK
      
      DO i=1,Solver % NumberOfActiveElements
        Element => Mesh % Elements(Solver % ActiveElements(i))
        IF (CheckPassiveElement(Element)) THEN
          nd = mGetElementDOFs(Indexes,UElement=Element)
          n  = Element % Type % NumberOfNodes
          DO j=1,nd
            k = Indexes(j)
            IF (k<=0) CYCLE

            IF(k<=SIZE(PassPerm)) THEN
              IF(PassPerm(k) > 0) CYCLE
            END IF
            
            k=Perm(k)
            IF (k<=0) CYCLE

            DO l=1,NDOFs
              m=NDOFs*(k-1)+l
#if 0 
              ! I don't trust this piece of code for parallel interfaces
              s=ABS(A % Values(A % Diag(m)))
              IF (s>EPSILON(s)) CYCLE
#endif
              
              m = NDOFs*(k-1)+l
              IF(A % ConstrainedDOF(m)) CYCLE
              CALL SetSinglePoint(k,l,Solver % Variable % Values(m),.FALSE.)
            END DO
          END DO
        END IF
      END DO
      DEALLOCATE(PassPerm,NodeIndexes)
    END IF

    ! Check the boundaries and body forces for possible single nodes BCs that must have a constant
    ! value on that boundary / body force.
    !--------------------------------------------------------------------------------------------
    NeedListMatrix = .FALSE.
    
    DO l = 0, 1
      IF( l == 0 ) THEN
        DirName = TRIM(Name)//' Constant'
      ELSE
        DirName = TRIM(Name)//' Profile'
      END IF
      
      AnySingleBC = ListCheckPresentAnyBC( Model, DirName )
      AnySingleBF = ListCheckPresentAnyBodyForce( Model, DirName )

      IF( AnySingleBC .OR. AnySingleBF ) THEN
        GotMult = (l == 1 )
        ALLOCATE( LumpedNodeSet( SIZE( Perm ) ) )

        IF( AnySingleBC ) CALL Info(Caller,'Found BC constraint: '//TRIM(DirName),Level=6)
        IF( AnySingleBF ) CALL Info(Caller,'Found BodyForce constraint: '//TRIM(DirName),Level=6)

        ! Improve the logic in future
        ! Now we assume that if the "supernode" has been found then also the matrix has the correct topology. 
        IF( AnySingleBC ) THEN
          NeedListMatrix = .NOT. ListCheckPresentAnyBC( Model, TRIM(Name)//' Constant Node Index')
        ELSE 
          NeedListMatrix = .NOT. ListCheckPresentAnyBodyForce( Model, TRIM(Name)//' Constant Node Index')
        END IF

        ! Move the list matrix because of its flexibility. Register the initial topology.
        IF( NeedListMatrix ) THEN
          CALL CRS_ChangeTopology(A, Init=.TRUE.)
          CALL List_toListMatrix(A)
        END IF
      
        DO bc = 1,Model % NumberOfBCs + Model % NumberOfBodyForces
          ! Make a distinction between BCs and BFs. 
          ! These are treated in the same loop because most of the logic is still the same. 
          IF( bc <= Model % NumberOfBCs ) THEN
            IF(.NOT. AnySingleBC ) CYCLE
            ValueList => Model % BCs(BC) % Values
            ElemFirst =  Mesh % NumberOfBulkElements + 1 
            ElemLast =  Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
          ELSE
            IF(.NOT. AnySingleBF ) CYCLE
            ValueList => Model % BodyForces(bc - Model % NumberOfBCs) % Values
            ElemFirst =  1
            ElemLast =  Model % NumberOfBulkElements
          END IF

          IF(GotMult) THEN
            IF(.NOT. ListCheckPresent( ValueList,DirName) ) CYCLE
          ELSE            
            IF(.NOT. ListGetLogical( ValueList,DirName, GotIt) ) CYCLE
          END IF
                    
          MoveCoeff = ListGetCReal( ValueList,TRIM(DirName)//' Resistance',GotIt)
          IF(.NOT. GotIt) MoveCoeff = 1.0_dp

          ! This tells us that this has been visited before            
          ind = ListGetInteger( ValueList,TRIM(DirName)//' Node Index',GotIt )               
         
          ! On the first time find a one single uniquely defined node for setting 
          ! the value. In parallel it will be an unshared node with the highest possible 
          ! node number 
          IF( GotIt ) THEN        
            dgind = ListGetInteger( ValueList,TRIM(DirName)//' DG Node Index',GotIt )               
            IF( GotMult ) THEN
              MaxMult = ListGetConstReal( ValueList,TRIM(DirName)//' Max Value',UnfoundFatal=.TRUE.) 
            END IF
          ELSE
            MaxMult = 0.0_dp          
            ind = 0
            dgind = 0
            maxind = 0
            
            DO t = ElemFirst, ElemLast
              Element => Mesh % Elements(t)

              IF( bc <= Model % NumberOfBCs ) THEN
                IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
              ELSE
                bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values,'Body Force',GotIt)
                IF( bc - Model % NumberOfBCs /= bf ) CYCLE
              END IF

              n = Element % TYPE % NumberOfNodes
              NodeIndexes => Element % NodeIndexes
                
              IF(GotMult) THEN
                Mult(1:n) = ListGetReal( ValueList,TRIM(DirName),n,NodeIndexes,UnfoundFatal=.TRUE.)
              END IF

              DO i=1,n
                j = NodeIndexes(i)

                IF(Model % Solver % DG) THEN                                   
                  CALL PickDgIndexes(Element,Indexes)                  
                  IF( Perm(Indexes(i)) == 0) CYCLE
                ELSE
                  IF( Perm(j) == 0) CYCLE
                END IF
                  
                IF( Parallel ) THEN
                  IF( SIZE( Mesh % ParallelInfo % NeighbourList(j) % Neighbours) > 1 ) CYCLE               
                  IF( Mesh % ParallelInfo % NeighbourList(j) % Neighbours(1) /= ParEnv % MyPe ) CYCLE               
                END IF
                IF( GotMult) THEN
                  ! Find the point with maximum value of the multiplier
                  IF( ABS(Mult(i)) > ABS(MaxMult) ) THEN
                    maxind = j
                    MaxMult = Mult(i)
                  END IF
                ELSE
                  ind = j 
                  IF( Model % Solver % DG ) dgind = Indexes(i)
                  EXIT
                END IF
              END DO
              IF( ind > 0 ) EXIT
            END DO

            IF(GotMult) ind = maxind

            ! Find the maximum partition that owns the node. 
            ! It could be minimum also, just some convention is needed. 
            IF( Parallel ) THEN
              IF( GotMult ) THEN
                ParMaxMult = ABS(MaxMult)
                ParMaxMult = ParallelReduction( ParMaxMult, 2 ) 
                IF(ABS(ABS(MaxMult)-ParMaxMult) > 1.0e-3*ParMaxMult) ind = 0 
              END IF
              k = -1
              IF( ind > 0 ) k = ParEnv % MyPe          
              k = ParallelReduction( k, 2 ) 
              IF( k == -1 ) THEN
                CALL Warn(Caller,'Could not find node to set: '//TRIM(DirName))
              END IF
              j = 0
              IF( k /= ParEnv % MyPe .AND. ind > 0) THEN
                ind = 0
                dgind = 0
                j = 1 
              END IF
              ! Just a counter for partitions that have hits but do not own the index.
              j = ParallelReduction(j)
              IF(j>0) THEN
                CALL Warn(Caller,'Problems: boundary extends over '//I2S(j+1)//' partitions')
              END IF
            ELSE
              IF( ind == 0 ) THEN
                CALL Warn(Caller,'Could not find node to set: '//TRIM(DirName))
              END IF
            END IF

            CALL ListAddInteger( ValueList,TRIM(DirName)//' Node Index', ind )
            IF( Model % Solver % DG ) THEN
              CALL ListAddInteger( ValueList,TRIM(DirName)//' DG Node Index', dgind )
            END IF

            IF(GotMult) THEN
              CALL ListAddConstReal( ValueList,TRIM(DirName)//' Max Value', MaxMult )
            END IF
            WRITE(Message,'(A,ES12.3)') TRIM(DirName)//' Max Value', MaxMult
            CALL Info(Caller,Message,Level=15)

            NeedListMatrix = .TRUE.
          END IF

          ! This is probably in parallel a passive partition. 
          IF( ind == 0 ) CYCLE

          ! Ok, now sum up the rows to the corresponding nodal index
          ! We go through elements and need to mark the nodes in order not to visit them twice. 
          LumpedNodeSet = .FALSE.

          ! Actually from this on we do not need the node index if we have a DG field.
          ! Replace it with the DG index.
          IF( Model % Solver % DG ) THEN
            IF( dgind == 0 ) THEN
              CALL Fatal(Caller,'Could not determine DG index associated to node index!')
            END IF
            ind = dgind
          END IF

          ! If we need to do scaling then play with the supernode too!
          IF( Ndofs == 1 ) THEN
            BLOCK
              INTEGER :: k0
              k0 = Offset + NDOFs * (Perm(ind)-1) + DOF
              IF(ABS(MoveCoeff-1.0) > EPSILON(MoveCoeff)) THEN
                CALL MoveRow( A, k0, k0, MoveCoeff )
                b(k0) = MoveCoeff * b(k0)
              END IF
            END BLOCK
          END IF
          ! Supernode has been set, if needed. 
          LumpedNodeSet(ind) = .TRUE.
                    
          IF(.NOT. ALLOCATED(LumpedIndx)) THEN
            ALLOCATE(LumpedIndx(Model % NumberOfBCs + Model % NumberOfBodyForces))
            LumpedIndx = 0
          END IF
          LumpedIndx(bc) = ind
            
          DO t = ElemFirst, ElemLast
            Element => Mesh % Elements(t)

            IF( bc <= Model % NumberOfBCs ) THEN
              IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
            ELSE
              bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values,'Body Force',GotIt)
              IF( bc - Model % NumberOfBCs /= bf ) CYCLE
            END IF

            n = Element % TYPE % NumberOfNodes
            IF ( Model % Solver % DG ) THEN              
              CALL PickDgIndexes(Element,Indexes)                  
            ELSE
              Indexes(1:n) = Element % NodeIndexes
            END IF
                
            IF(GotMult) Mult(1:n) = ListGetReal( ValueList,TRIM(Name)//' Profile',n,&
                Element % NodeIndexes,UnfoundFatal=.TRUE.)
                       
            CALL SetLumpedRows(ind,n)
          END DO

          IF( NDOFs == 1 ) THEN
            SingleVal = ListGetCReal( ValueList,TRIM(DirName)//' Flux',GotIt)
            IF( GotIt ) THEN
              t = Offset + Perm(ind)
              b(t) = b(t) + SingleVal
            END IF

            SingleVal = ListGetCReal( ValueList,TRIM(DirName)//' Coefficient',GotIt)
            IF( GotIt ) THEN                        
              t = Offset + Perm(ind)
              CALL AddToMatrixElement(A,t,t,SingleVal) 
            END IF                        
          END IF

          n = COUNT( LumpedNodeSet ) 
          CALL Info(Caller,'Number of lumped nodes set: '//I2S(n),Level=10)
        END DO
      END IF
    END DO
      
    IF( NeedListMatrix ) THEN
      DEALLOCATE( LumpedNodeSet )
      
      ! Revert back to CRS matrix and change to the new topology. 
      CALL List_ToCRSMatrix(A)
      CALL CRS_ChangeTopology(A, Init=.FALSE.)

      CALL Info(Caller,'Modified matrix non-zeros: '&
          //I2S(SIZE( A % Cols )),Level=8)
    END IF
      
    ! We are back to CRS matrix.
    ! If we have a fixed point iteration we may add the flux multiplied by resistance to rhs as well. 
    BLOCK
      INTEGER :: k0
      REAL(KIND=dp) :: prevFlux
      TYPE(ValueList_t), POINTER :: vList
      IF(ALLOCATED(LumpedIndx) ) THEN
        DO bc = 1,SIZE(LumpedIndx) 
          ind = LumpedIndx(bc)
          IF(ind>0) THEN
            k0 = Offset + NDOFs * (Perm(ind)-1) + DOF
            prevFlux = CRS_MatrixRowVectorMultiply(A,Model % Solver % Variable % Values,k0)            

            ! In case this was added, remove it from the flux.
            SingleVal = ListGetCReal( ValueList,TRIM(Name)//' Constant Coefficient',GotIt)
            prevFlux = prevFlux - SingleVal * Model % Solver % Variable % Values(k0)
            
            IF( bc <= Model % NumberOfBCs ) THEN
              Vlist => Model % BCs(bc) % Values
              WRITE(Message,'(A,ES12.3)') 'Previous bc '//I2S(bc)//' lumped flux: ',prevFlux 
            ELSE
              i = bc-Model % NumberOfBCs
              Vlist => Model % BodyForces(i) % Values
              WRITE(Message,'(A,ES12.3)') 'Previous bf '//I2S(i)//' lumped flux: ',prevFlux 
            END IF            
            CALL Info(Caller,Message)
            DirName = TRIM(Name)//' Constant Prev Flux'
            CALL ListAddConstReal(Vlist, DirName, prevFlux )
          END IF
        END DO
      END IF
    END BLOCK

    ! Check the boundaries for a curve bc.
    !--------------------------------------------------------------------------------------------
    BLOCK
      INTEGER :: l,dofs,c1,c2,c3,ivec
      REAL(KIND=dp) :: Coeff,d(3),Dfdx(3),f,x(3),a11,a22
      LOGICAL, ALLOCATABLE :: NodeDone(:)
      TYPE(Variable_t), POINTER :: Dvar
      INTEGER, POINTER :: GivenComps(:)
      INTEGER :: Comps(3)
      LOGICAL :: AnyHingeBC, HingeBC
      REAL(KIND=dp) :: Normal(3),Tan1(3),Tan2(3),xt(2),cfit(7)
      
      
      DirName = TRIM(Name)//' Curve'
      AnySingleBC = ListCheckPresentAnyBC( Model, DirName )
      AnyHingeBC = ListGetLogicalAnyBC( Model, TRIM(Name)//' Hinge' )
      
      IF( AnySingleBC .OR. AnyHingeBC ) THEN
        IF( AnySingleBC ) THEN
          CALL Info(Caller,'Found BC constraint for curve: '//TRIM(DirName),Level=6)
        END IF          
        IF( AnySingleBC ) THEN
          CALL Info(Caller,'Found BC constraint for hinge: '//TRIM(Name)//' Hinge',Level=6)
        END IF
          
        Solver => Model % Solver
                
        Dvar => Solver % Variable
        IF(.NOT. ASSOCIATED( DVar ) ) THEN
          CALL Fatal(Caller,'Solver variable not associated for Curve constraint!')
        END IF
        
        dim = CoordinateSystemDimension()        
        dofs = MIN( dim, DVar % Dofs )
        IF( dofs < 2 ) THEN
          CALL Fatal(Caller,'Curve constraint only makes sense for vector fields!')
        END IF

        ALLOCATE(NodeDone(Mesh % NumberOfNodes))
        NodeDone = .FALSE.
        
        DO bc = 1,Model % NumberOfBCs
          
          ValueList => Model % BCs(BC) % Values

          HingeBC = .FALSE.
          IF( AnyHingeBC ) THEN
            HingeBC = ListGetLogical( ValueList,TRIM(Name)//' Hinge',GotIt)
          END IF
          
          IF(HingeBC ) THEN            
            IF( dim == 2 ) THEN
              CALL CylinderFit(Mesh, ValueList, bc, dim, cfit )
            ELSE IF( dim == 3 ) THEN
              CALL CylinderFit(Mesh, ValueList, bc, dim, cfit )
              Normal = cfit(4:6)
              CALL TangentDirections(Normal,Tan1,Tan2)
            END IF
            !PRINT *,'Hinge Params',cfit
          ELSE           
            IF( .NOT. ListCheckPresent( ValueList,DirName) ) CYCLE
          END IF

          Comps = 0
          GivenComps => ListGetIntegerArray( ValueList, TRIM(Dirname)//' Components',GotIt)
          IF(GotIt) THEN
            dofs = MIN(dofs, SIZE(GivenComps) )
            Comps(1:dofs) = GivenComps(1:dofs)
          ELSE
            DO i=1,dofs
              Comps(i) = i
            END DO
          END IF
                    
          ElemFirst = Mesh % NumberOfBulkElements + 1 
          ElemLast = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
                         
          DO t = ElemFirst, ElemLast
            Element => Mesh % Elements(t)

            IF(.NOT. ASSOCIATED(Element % BoundaryInfo) ) CYCLE
            IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
            n = Element % TYPE % NumberOfNodes
            Indexes(1:n) = Element % NodeIndexes
            
            DO i=1,n
              j = Indexes(i)
              IF(NodeDone(j)) CYCLE
              NodeDone(j) = .TRUE.
              k = DVar % Perm(j)

              IF(k==0) CYCLE

              l = DVar % dofs*(k-1)
              d(1:dofs) = DVar % Values(l+1:l+dofs)
              IF(dofs==2) d(3) = 0.0_dp
              
              x(1) = Mesh % Nodes % x(j) 
              x(2) = Mesh % Nodes % y(j)
              x(3) = Mesh % Nodes % z(j)

              x = x + d
              
              IF( HingeBC ) THEN
                ! We can analytically derive the case of 2d hinge
                IF( dim == 2 ) THEN
                  x(1:2) = x(1:2) - cfit(1:2)
                  f = x(1)**2 + x(2)**2 - cfit(3)**2
                  DfDx(1) = 2*x(1)
                  DfDx(2) = 2*x(2)
                ELSE
                  x = x - cfit(1:3)
                  xt(1) = SUM(x * Tan1)
                  xt(2) = SUM(x * Tan2)
                                    
                  f = xt(1)**2 + xt(2)**2 - cfit(7)**2
                  
                  DfDx(1) = 2*xt(1)*Tan1(1) + 2*xt(2)*Tan2(1)
                  DfDx(2) = 2*xt(1)*Tan1(2) + 2*xt(2)*Tan2(2)
                  DfDx(3) = 2*xt(1)*Tan1(3) + 2*xt(2)*Tan2(3)                  
                END IF
              ELSE                
                f = ListGetFunVec( ValueList, DirName, x(1:dofs), dofs, DfDx=dfdx(1:dofs) )
              END IF

              ! Check whether this is an normal-tangential node. If it is then m>0. 
              ! We already know that this has active perm for DVar. 
              m = 0
              IF ( NT % NormalTangentialNOFNodes > 0 ) THEN
                m = NT % BoundaryReorder(j) 
              END IF
                                         
              IF( m == 0 ) THEN                
                ! Let us take the most sensitive component to be the one for
                ! which the curve constraint is applied ensuring maximum diagonal entry.              ´              
                ! Then choose 2nd (and 3rd) components in order. 
                c1 = comps(1)
                DO ivec=2,dofs
                  IF( ABS(DfDx(comps(ivec))) > ABS(DfDx(c1))) c1 = comps(ivec)
                END DO
                DO ivec=1,dofs
                  IF(comps(ivec) /= c1) THEN
                    c2 = comps(ivec)
                    EXIT
                  END IF
                END DO

                IF( dofs == 3 ) THEN
                  c3 = 6 - c1 - c2
                END IF
               
                ! It may happen that the rows are linearly dependent on each other.
                ! Then a solution of type (x+y) for both is not good. Rather use then
                ! (x-y) for the other to have a unique solution. 
                a11 = GetMatrixElement(A,l+c1,l+c1)
                a22 = GetMatrixElement(A,l+c2,l+c2)
                
                ! This is a simple sign rule that avoids two equations being redundant.
                ! Don't multiply too numbers that could be almost zero!
                Coeff = -SIGN(1.0_dp,a11) * SIGN(1.0_dp,a22) * &
                    SIGN(1.0_dp, DfDx(c1)) * SIGN(1.0_dp, DfDx(c2))
              ELSE
                ! For normal-tangential system the normal component (1st one) should
                ! by construction be most sensitive to deviations from the curve. 
                c1 = 1
                c2 = 2
                c3 = 3
                Coeff = 1.0_dp
              END IF
                                                      
              ! Move all the entries from "c1" to "c2" and nullify the row.
              CALL MoveRow( A, l+c1, l+c2, Coeff )
              b(l+c2) = b(l+c2) + Coeff * b(l+c1)
              b(l+c1) = 0.0_dp

              ! In parallel only the owner makes the curve constraint
              IF( ParEnv % PEs > 1 ) THEN
                IF( A % ParallelInfo % NeighbourList(l+1) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
              END IF

              ! Residual mode is never active at this stage
              b(l+c1) = -f              
              b(l+c1) = b(l+c1) + DfDx(1)*d(1) + DfDx(2)*d(2)
              IF(dofs==3) b(l+c1) = b(l+c1) + DfDx(3)*d(3)
              
              IF(m>0) THEN                           
                CALL RotateNTSystem( DfDx, j )
              END IF

              CALL AddToMatrixElement( A, l+c1, l+1, DfDx(1) )
              CALL AddToMatrixElement( A, l+c1, l+2, DfDx(2) )
              IF(dofs==3) CALL AddToMatrixElement( A, l+c1, l+3, DfDx(3) )
            END DO
          END DO
        END DO
        
        n = COUNT( NodeDone ) 
        CALL Info(Caller,'Number of curved nodes set: '//I2S(n),Level=10)        
      END IF
      
    END BLOCK

    ! Check the boundaries and body forces for possible single nodes BCs that must have a constant
    ! value on that boundary / body force.
    !--------------------------------------------------------------------------------------------
    DirName = TRIM(Name)//' Plane'
    AnySingleBC = ListCheckPresentAnyBC( Model, DirName )
    
    IF( AnySingleBC ) THEN
      dim = CoordinateSystemDimension()
      
      ALLOCATE( LumpedNodeSet( SIZE( Perm ) ) )

      CALL Info(Caller,'Found BC constraint: '//TRIM(DirName),Level=6)

      ! Improve the logic in future
      ! Now we assume that if the "supernode" has been found then also the matrix has the correct topology. 
      NeedListMatrix = .NOT. ListCheckPresentAnyBC( Model, TRIM(Name)//' Plane Node Indices')
      
      ! Move the list matrix because of its flexibility
      IF( NeedListMatrix ) THEN
        CALL CRS_ChangeTopology(A,Init=.TRUE.)
        CALL List_toListMatrix(A)
      END IF

      ElemFirst = Mesh % NumberOfBulkElements + 1 
      ElemLast = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      DO bc = 1,Model % NumberOfBCs 

        ValueList => Model % BCs(BC) % Values
        IF( .NOT. ListGetLogical( ValueList,DirName, GotIt) ) CYCLE

        PlaneInds => ListGetIntegerArray( ValueList,TRIM(Name)//' Plane Node Indices',GotIt )     

        IF(.NOT. GotIt ) THEN
          IF(.NOT. ALLOCATED(CandNodes) ) THEN
            ALLOCATE( CandNodes( Mesh % NumberOfNodes ) )        
          END IF
          CandNodes = .FALSE.

          ! Add nodes to the set that are associated with this BC only.
          DO t = ElemFirst, ElemLast
            Element => Mesh % Elements(t)            
            IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
              NodeIndexes => Element % NodeIndexes
              CandNodes(NodeIndexes) = .TRUE.
            END IF
          END DO

          ! Remove nodes from the set that may be set by other BCs also. 
          DO t = ElemFirst, ElemLast
            Element => Mesh % Elements(t)            
            IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) THEN
              NodeIndexes => Element % NodeIndexes
              CandNodes(NodeIndexes) = .FALSE.
            END IF
          END DO

          ALLOCATE(PlaneInds(3))
          CALL FindExtremumNodes(Mesh,CandNodes,dim,PlaneInds) 
          
          CALL ListAddIntegerArray( ValueList,TRIM(Name)//' Plane Node Indices',dim, PlaneInds )
          NeedListMatrix = .TRUE.
        END IF

        IF( Parallel ) CALL Warn(Caller,'Node index perhaps not set properly in parallel')
        ! IF( ind == 0 ) CYCLE

        ! Ok, now sum up the rows to the corresponding nodal index
        LumpedNodeSet = .FALSE.

        ! Don't lump the "supernodes" and therefore mark it set already
        LumpedNodeSet(PlaneInds) = .TRUE.

        DO t = ElemFirst, ElemLast
          Element => Mesh % Elements(t)

          IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
            n = Element % TYPE % NumberOfNodes
            Indexes(1:n) = Element % NodeIndexes
            CALL SetRigidRows(PlaneInds,bc,n)
          END IF
        END DO

        n = COUNT( LumpedNodeSet ) 
        CALL Info(Caller,'Number of lumped nodes set: '//I2S(n),Level=10)
      END DO

      IF( NeedListMatrix ) THEN
        DEALLOCATE( LumpedNodeSet )
        
        ! Revert back to CRS matrix
        CALL List_ToCRSMatrix(A)
        CALL CRS_ChangeTopology(A,Init=.FALSE.)
      END IF
    END IF

    IF( InfoActive(12) )  THEN
      IF( Parallel ) THEN
        DirCount = ParallelReduction( DirCount ) 
      END IF
      CALL Info(Caller,'Number of dofs set for '//TRIM(Name)//': '&
          //I2S(DirCount),Level=12)
    END IF
!------------------------------------------------------------------------------

  CONTAINS

    ! For boundary element in a DG mesh get the DG indexes from itself, or
    ! from the parent elements.
    !-------------------------------------------------------------------------
    SUBROUTINE PickDgIndexes(Element,DgIndexes)
      TYPE(Element_t) :: Element
      INTEGER :: DGIndexes(:)
      
      TYPE(Element_t), POINTER :: Parent
      INTEGER :: i,j,lr,n,m

      n = Element % Type % NumberOfNodes
       
      IF( ASSOCIATED( Element % DgIndexes ) ) THEN
        DGIndexes(1:n) = Element % DGIndexes(1:n)
        IF(ANY(DgIndexes(1:n) == 0) ) THEN
          ! If this fails, maybe used the 2nd strategy always...
          CALL Fatal('PickDgIndexes','Some of BC element DG indexes is zero!')
        END IF
      ELSE IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
        m = 0
        DO lr=1,2
          IF(lr==1) THEN
            Parent => Element % BoundaryInfo % Left
          ELSE
            Parent => Element % BoundaryInfo % Right
          END IF
          IF(.NOT. ASSOCIATED( Parent ) ) CYCLE
          IF(.NOT. ASSOCIATED( Parent % DGIndexes ) ) CYCLE
          IF(ANY(Perm(Parent % DGIndexes) == 0)) CYCLE
          
          DO i=1,Element % TYPE % NumberOfNodes
            DO j=1,Parent % TYPE % NumberOfNodes
              IF( Element % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
                DGIndexes(i) = Parent % DGIndexes(j)
                m = m+1
                EXIT
              END IF
            END DO
          END DO
          EXIT
        END DO
        IF(m /= n) THEN
          CALL Fatal('PickDgIndexes','Could not find all DG Indexes for BC element')
        END IF
      END IF
      
    END SUBROUTINE PickDgIndexes

    
     ! Check n-t node setting element
     !-------------------------------
    SUBROUTINE CheckNTElement(n,elno)
      INTEGER :: n,elno
      INTEGER :: i,j,k,l,m,dim,kmax
      LOGICAL :: found
      REAL(KIND=dp) :: Condition(n), RotVec(3)
      
      dim = CoordinateSystemDimension()

      IF ( DOF <= 0 ) RETURN
      
      IF ( NT % NormalTangentialNOFNodes == 0 ) RETURN
      IF ( ALL(NT % BoundaryReorder(Indexes(1:n))<1) ) RETURN
      IF ( .NOT. ListCheckPresent(ValueList, Name) ) RETURN
      IF ( ListGetLogical(ValueList,NT % NormalTangentialName,Found) ) RETURN

      IF ( Conditional ) THEN
        Condition(1:n) = ListGetReal( ValueList, CondName, n, Indexes, gotIt )
        Conditional = Conditional .AND. GotIt
      END IF

      !
      ! Check for nodes belonging to n-t boundary getting set by other bcs.
      ! -------------------------------------------------------------------
      DO j=1,n
        IF ( Conditional .AND. Condition(j)<0.0_dp ) CYCLE
        k = Perm(Indexes(j))
        IF ( k > 0 ) THEN          
          k = k + OffSet
          m = NT % BoundaryReorder(Indexes(j))
          IF ( m>0 ) THEN
            RotVec = 0._dp
            RotVec(DOF) = 1._dp
            CALL RotateNTSystem( RotVec, Indexes(j) )
            kmax = 1
            DO k=1,dim
              IF ( ABS(RotVec(k)) > ABS(RotVec(kmax)) ) kmax = k
            END DO
            NTelement(m,kmax)=elno
          END IF
        END IF
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE CheckNTElement
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Set values related to a specific boundary or bulk element.
!------------------------------------------------------------------------------
    SUBROUTINE SetElementValues(n)
!------------------------------------------------------------------------------
      INTEGER :: n
      INTEGER :: i,j,k,l,m,dim,kmax,lmax
      LOGICAL :: CheckNT,found
      REAL(KIND=dp) :: Condition(n), Work(n), RotVec(3)
      
      dim = CoordinateSystemDimension()

      IF ( DOF > 0 ) THEN
        IF (Model % Solver % DG) THEN
          Work(1:n)  = ListGetReal( ValueList, Name, n, Element % NodeIndexes, gotIt )
        ELSE
          Work(1:n)  = ListGetReal( ValueList, Name, n, Indexes, gotIt )
        END IF
        IF ( .NOT. GotIt ) THEN
          Work(1:n)  = ListGetReal( ValueList, Name(1:nlen) // ' DOFs', n, Indexes, gotIt )
        END IF
      ELSE
        CALL ListGetRealArray( ValueList, Name, WorkA, n, Indexes, gotIt )
      END IF
      
      IF ( gotIt ) THEN
        IF ( Conditional ) THEN
          IF (Model % Solver % DG) THEN
            Condition(1:n) = ListGetReal( ValueList, CondName, n, Element % NodeIndexes, gotIt )
          ELSE
            Condition(1:n) = ListGetReal( ValueList, CondName, n, Indexes, gotIt )
          END IF
          Conditional = Conditional .AND. GotIt
        END IF

        !
        ! Check for nodes belonging to n-t boundary getting set by other bcs.
        ! Here we don't need to track p-bubbles as they are not shared by
        ! many BCs. 
        ! -------------------------------------------------------------------
        CheckNT = .FALSE.
        IF ( NT % NormalTangentialNOFNodes>0 .AND. DOF>0 ) THEN
          CheckNT = .TRUE.
          IF ( ALL(NT % BoundaryReorder(Indexes(1:n))<1) ) CheckNT = .FALSE.
          IF ( ListGetLogical(ValueList,NT % NormalTangentialName,Found)) CheckNT=.FALSE.
        END IF
        
        DO j=1,n
          IF ( Conditional ) THEN
            IF( Condition(j) < 0.0_dp ) CYCLE
          END IF

          k = Perm(Indexes(j))
          IF ( k > 0 ) THEN
            
            IF ( DOF>0 ) THEN
              m = 0
              IF ( NT % NormalTangentialNOFNodes>0 ) m = NT % BoundaryReorder(Indexes(j))
              IF ( m>0 .AND. CheckNT ) THEN
                RotVec = 0._dp
                RotVec(DOF) = 1._dp
                CALL RotateNTSystem( RotVec, Indexes(j) )

                ! When cartesian component "DOF" is defined set the N-T component
                ! closest to its direction. 
                kmax = 1 
                DO k=2,dim
                  IF ( ABS(RotVec(k)) > ABS(RotVec(kmax)) ) THEN
                    kmax = k
                  END IF
                END DO

                lmax = NDOFs * (Perm(Indexes(j))-1) + kmax
                IF ( .NOT. NTZeroing_done(m,kmax) ) THEN
                  b(lmax) = 0._dp

                  IF( .NOT. OffDiagonal ) THEN
                    b(lmax) = b(lmax) + Work(j) 
                  END IF

                  ! Consider all components of the cartesian vector mapped to the 
                  ! N-T coordinate system. Should this perhaps have scaling included?
                  DirCount = DirCount + 1
                  CALL ZeroRow( A,lmax )
                  IF( .NOT. OffDiagonal) THEN
                    DO k=1,dim
                      l = NDOFs * (Perm(Indexes(j))-1) + k
                      CALL SetMatrixElement( A,lmax,l,RotVec(k) )
                    END DO
                  END IF

                  NTZeroing_done(m,kmax)   = .TRUE.
                  A % ConstrainedDOF(lmax) = .FALSE.
                END IF
              ELSE
                CALL SetSinglePoint(k,DOF,Work(j),.FALSE.)
              END IF
            ELSE
              DO l=1,MIN( NDOFs, SIZE(WorkA,1) )
                CALL SetSinglePoint(k,l,WorkA(l,1,j),.FALSE.)
              END DO
            END IF
          END IF
        END DO
      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE SetElementValues
!------------------------------------------------------------------------------
  


!------------------------------------------------------------------------------
!> Set values related to a specific boundary or bulk element.
!------------------------------------------------------------------------------
    SUBROUTINE SetLumpedRows(ind0,n)
!------------------------------------------------------------------------------
      INTEGER :: ind0,n
      INTEGER :: ind,i,j,k,k0,l
      REAL(KIND=dp) :: Coeff
      ! -------------------------------------------------------------------        
      
      DO j=1,n
        ind = Indexes(j)

        IF( LumpedNodeSet(ind) ) CYCLE
        LumpedNodeSet(ind) = .TRUE.

        IF ( DOF > 0 ) THEN
          k0 = Offset + NDOFs * (Perm(ind0)-1) + DOF
          k = OffSet + NDOFs * (Perm(ind)-1) + DOF

          IF( GotMult ) THEN
            Coeff = Mult(j) / MaxMult
          ELSE
            Coeff = 1.0_dp
          END IF
          
          CALL MoveRow( A, k, k0, MoveCoeff )
          b(k0) = b(k0) + MoveCoeff * b(k)

          CALL AddToMatrixElement( A, k, k, 1.0_dp )
          CALL AddToMatrixElement( A, k, k0, -Coeff )
          b(k) = 0.0_dp
        ELSE
          DO l = 1, NDOFs
            k0 = Offset + NDOFs + (Perm(ind0)-1) * DOF
            k = OffSet + NDOFs * (Perm(ind)-1) + l

            IF( GotMult ) THEN
              Coeff = Mult(j) / MaxMult
            ELSE
              Coeff = 1.0_dp
            END IF
              
            CALL MoveRow( A, k, k0, 1.0_dp )
            b(k0) = b(k0) + Coeff * b(k)
          
            CALL AddToMatrixElement( A, k, k, 1.0_dp )
            CALL AddToMatrixElement( A, k, k0, -Coeff )
          END DO
        END IF
      END DO

!------------------------------------------------------------------------------
    END SUBROUTINE SetLumpedRows
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Set values related to a rigid plane boundary such that any node on the boundary
!> is expressed as linear combination of the selected three (or two if on line)
!> nodes.
!------------------------------------------------------------------------------
    SUBROUTINE SetRigidRows(inds0,bcind,n)
!------------------------------------------------------------------------------
      INTEGER :: inds0(:)
      INTEGER :: bcind
      INTEGER :: n

      INTEGER :: bcind0 =  0
      INTEGER :: ind,i,j,k,k0
      REAL(KIND=dp) :: Coeff, Weights(3)
      REAL(KIND=dp) :: BaseCoord(3,3),r1(3),r2(3),Coord(3),dCoord(3),Amat(2,2),A0mat(2,2),bvec(2)
      
      SAVE bcind0, BaseCoord, A0mat, r1, r2
!-------------------------------------------------------------------        

      IF(bcind /= bcind0 ) THEN
        BaseCoord = 0.0_dp
        DO i=1,dim
          j = inds0(i)
          BaseCoord(i,1) = Mesh % Nodes % x(j)
          BaseCoord(i,2) = Mesh % Nodes % y(j)
          BaseCoord(i,3) = Mesh % Nodes % z(j)
        END DO
        bcind0 = bcind

        r1 = BaseCoord(2,:) - BaseCoord(1,:)
        Amat(1,1) = SUM(r1*r1)
        
        IF( dim == 3 ) THEN
          r2 = BaseCoord(3,:) - BaseCoord(1,:)
          Amat(1,2) = SUM(r1*r2)
          Amat(2,1) = Amat(1,2)
          Amat(2,2) = SUM(r2*r2)
        END IF

        A0mat = Amat
        bcind0 = bcind
      END IF
                   
      DO j=1,n
        ind = Indexes(j)

        IF( LumpedNodeSet(ind) ) CYCLE
        LumpedNodeSet(ind) = .TRUE.
        
        Coord(1) = Mesh % Nodes % x(ind)
        Coord(2) = Mesh % Nodes % y(ind)
        Coord(3) = Mesh % Nodes % z(ind)

        dCoord = Coord - BaseCoord(1,:)
        
        bvec(1) = SUM( dCoord * r1 )
        IF( dim == 3 ) THEN
          bvec(2) = SUM( dCoord * r2 )
        END IF

        IF( dim == 2 ) THEN
          bvec(1) = bvec(1) / A0mat(1,1)
          Weights(2) = bvec(1)
          Weights(1) = 1.0_dp - Weights(2)
        ELSE
          Amat = A0mat          
          CALL LUSolve(2,Amat,bvec)          
          Weights(2:3) = bvec(1:2)
          Weights(1) = 1.0_dp - SUM(bvec(1:2))
        END IF

        DO l = 1, dim
          k = OffSet + NDOFs * (Perm(ind)-1) + l    

          ! Distribute row in accordance with the weights
          DO m = 1, dim
            k0 = Offset + NDOFs * (Perm(inds0(m))-1) + l          
            Coeff = Weights(m)

            b(k0) = b(k0) + Coeff * b(k)
            IF( m < dim ) THEN
              ! This does not nullify the row
              CALL MoveRow( A, k, k0, Coeff, 1.0_dp )
            ELSE
              ! Now also nullify the row
              CALL MoveRow( A, k, k0, Coeff )
              b(k) = 0.0_dp
            END IF
          END DO

          ! Express the node as linear combination of the base nodes
          DO m = 1,dim
            k0 = Offset + NDOFs * (Perm(inds0(m))-1) + l          
            Coeff = Weights(m)            
            CALL AddToMatrixElement( A, k, k0, -Coeff )
          END DO
          CALL AddToMatrixElement( A, k, k, 1.0_dp )
        END DO

      END DO

!------------------------------------------------------------------------------
    END SUBROUTINE SetRigidRows
!------------------------------------------------------------------------------

    
!------------------------------------------------------------------------------
!> Set values related to individual points.
!------------------------------------------------------------------------------
    SUBROUTINE SetPointValues(n)
!------------------------------------------------------------------------------
      INTEGER :: n
      REAL(KIND=dp) :: Work(n), Condition(n)        

      INTEGER :: i,j,k,k1,l

      IF ( DOF > 0 ) THEN
        Work(1:n) = ListGetReal( ValueList, Name, n, NodeIndexes, gotIt )
      ELSE
        CALL ListGetRealArray( ValueList, Name, WorkA, n, NodeIndexes, gotIt )
      END IF

      IF ( gotIt ) THEN

        Condition(1:n) = 1.0d0
        IF ( Conditional ) THEN
          Condition(1:n) = ListGetReal( ValueList, CondName, n, NodeIndexes, gotIt )
          Conditional = Conditional .AND. GotIt
        END IF

        DO j=1,n
          ! This can (for example) happen with "Target Node" - keyword used in parallel...
          IF(NodeIndexes(j)>SIZE(Perm) .OR. NodeIndexes(j)<=0) CYCLE

          k = Perm(NodeIndexes(j))
          IF( k == 0 ) CYCLE

          IF ( Conditional .AND. Condition(j) < 0.0d0 ) CYCLE

          IF ( DOF>0 ) THEN
            CALL SetSinglePoint(k,DOF,Work(j),.FALSE.)
          ELSE
            DO l=1,MIN( NDOFs, SIZE(Worka,1) )
              CALL SetSinglePoint(k,l,WorkA(l,1,j),.FALSE.)
            END DO
          END IF

        END DO
      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE SetPointValues
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Set values related to one single point
!------------------------------------------------------------------------------
    SUBROUTINE SetSinglePoint(ind,DOF,val,ApplyPerm)
!------------------------------------------------------------------------------
      LOGICAL :: ApplyPerm
      INTEGER :: ind, DOF
      REAL(KIND=dp) :: val

      REAL(KIND=dp) :: s
      INTEGER :: i,j,k,k1,l


      IF(.NOT. ALLOCATED(A % ConstrainedDOF)) THEN
        ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
        A % ConstrainedDOF = .FALSE.
      END IF
      
      IF(.NOT. ALLOCATED(A % Dvalues)) THEN
        ALLOCATE(A % Dvalues(A % NumberOfRows))
        A % Dvalues = 0._dp
      END IF
      
      k = ind
      IF (ApplyPerm) k = Perm(ind)
      IF( k == 0 ) RETURN
      
      k = OffSet + NDOFs * (k-1) + DOF

      ! Do not add non-zero entries to pure halo nodes which are not associated with the partition.
      ! These are nodes could be created by the -halobc flag in ElmerGrid.
      IF( Parallel ) THEN
        IF( ASSOCIATED( A % ParallelInfo ) ) THEN
          IF( .NOT. ANY( A % ParallelInfo % NeighbourList(k) % Neighbours == ParEnv % MyPe ) ) THEN
            RETURN
          END IF
        END IF
      END IF

      DirCount = DirCount + 1
      
      IF( .NOT. OffDiagonal ) THEN
        A % Dvalues(k) =  val
      END IF
      A % ConstrainedDOF(k) = .TRUE.

!------------------------------------------------------------------------------
    END SUBROUTINE SetSinglePoint
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Set values related to upper and lower limiters.
!------------------------------------------------------------------------------
    SUBROUTINE SetLimiterValues(n)
!------------------------------------------------------------------------------
      INTEGER :: n
      REAL(KIND=dp) :: Work(n)

      Work(1:n)  = ListGetReal( ValueList, CondName, n, NodeIndexes, gotIt )

      IF ( gotIt ) THEN
        DO j=1,n
          k = Perm(NodeIndexes(j))
          IF( k == 0 ) CYCLE

          IF( .NOT. LimitActive(nDofs*(k-1)+dof)) CYCLE
          CALL SetSinglePoint(k,DOF,Work(j),.FALSE.)
        END DO
      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE SetLimiterValues
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> At first pass sum together the rows related to the periodic dofs.
!------------------------------------------------------------------------------
   SUBROUTINE SetPeriodicBoundariesPass1( Model, A, b, &
                      Name, DOF, NDOFs, Perm, This, Done )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model        !< The current model structure
    TYPE(Matrix_t), TARGET :: A   !< The global matrix
    REAL(KIND=dp) :: b(:)         !< The global RHS vector
    CHARACTER(LEN=*) :: Name      !< name of the dof to be set
    LOGICAL :: Done(:)            !< Has the node already been done. 
    INTEGER :: This               !< Number of the current boundary.
    INTEGER :: DOF                !< The order number of the dof
    INTEGER :: NDOFs              !< the total number of DOFs for this equation
    INTEGER, TARGET :: Perm(:)    !< The node reordering info, this has been generated at the
                                  !< beginning of the simulation for bandwidth optimization
!------------------------------------------------------------------------------
    INTEGER :: p,q,i,j,k,l,m,n,nn,ii,nlen,jmp,size0
    INTEGER, POINTER :: PPerm(:)
    LOGICAL :: GotIt, Found, Jump
    REAL(KIND=dp) :: Scale, weight, coeff
    TYPE(Matrix_t), POINTER :: F, G, Projector, Projector1
    TYPE(Variable_t), POINTER :: Var, WeightVar
    TYPE(ValueList_t), POINTER :: BC
    TYPE(MortarBC_t), POINTER :: MortarBC
    LOGICAL :: Parallel
!------------------------------------------------------------------------------

    nlen = LEN_TRIM(Name)
    BC => Model % BCs(This) % Values

    IF ( ListGetLogical( BC,& 
        'Periodic BC ' // Name(1:nlen), GotIt ) ) THEN
      IF( ListGetLogical( BC,'Antisymmetric BC',GotIt ) ) THEN
        Scale = 1.0_dp
      ELSE
        Scale = -1.0_dp
      END IF
    ELSE IF ( ListGetLogical( BC, &
        'Anti Periodic BC ' // Name(1:nlen), GotIt ) ) THEN
      Scale = 1.0d0
    ELSE 
      Scale = ListGetConstReal( BC, &
          'Periodic BC Scale ' // Name(1:nlen), GotIt) 
      IF(.NOT. GotIt ) RETURN      
    END IF
    
    Projector => Model % BCs(This) % PMatrix
    IF ( .NOT. ASSOCIATED(Projector) ) RETURN
    
!   For explicit conditions just create the dependency almost like a normal Dirichlet BC, 
!   For implicit one (otherwise) do the assembly of the projector:
!   ---------------------------------
    IF ( ListGetLogical( BC, &
        'Periodic BC Explicit', Found ) ) THEN
      
      Var => VariableGet( Model % Variables,Name(1:nlen) ) 
      
      DO i=1,Projector % NumberOfRows
        ii = Projector % InvPerm(i)
        IF( ii == 0 ) CYCLE
        k = Perm(ii)
        IF ( .NOT. Done(ii) .AND. k>0 ) THEN
          k = NDOFs * (k-1) + DOF
          A % Dvalues(k) = 0._dp
          A % ConstrainedDOF(k) = .TRUE.
          
          DO l = Projector % Rows(i), Projector % Rows(i+1)-1
            IF ( Projector % Cols(l) <= 0 ) CYCLE
            m = Perm( Projector % Cols(l) )
            IF ( m > 0 ) THEN
              m = NDOFs * (m-1) + DOF
              A % Dvalues(k) = A % Dvalues(k) - Scale * Projector % Values(l) * &
                  Var % Values(m)
            END IF
          END DO
        END IF
      END DO
      
    ELSE IF ( ListGetLogical( BC, &
        'Periodic BC Use Lagrange Coefficient', Found ) ) THEN

      Jump = ListCheckPresent( BC, &
          'Periodic BC Coefficient '//Name(1:nlen))
      
      IF( .NOT. ASSOCIATED( Model % Solver % MortarBCs ) ) THEN
        CALL Info('SetPeriodicBoundariesPass1',&
            'Allocating mortar BCs for solver',Level=8)
        ALLOCATE( Model % Solver % MortarBCs( Model % NumberOfBCs + Model % NumberOfBodyForces ) )
        DO i=1, Model % NumberOfBCs + Model % NumberOfBodyForces
          Model % Solver % MortarBCs(i) % Projector => NULL()
        END DO
      END IF
      
      IF( ASSOCIATED( Projector, &
          Model % Solver % MortarBCs(This) % Projector) ) THEN
        CALL Info('SetPeriodicBoundariesPass1','Using existing projector: '&
            //I2S(This),Level=8)
        RETURN
      END IF
      
      Model % Solver % MortarBCs(This) % Projector => Projector
      CALL Info('SetPeridociBoundariesPass1','Using projector as mortar constraint: '&
          //I2S(This),Level=8)

      MortarBC => Model % Solver % MortarBCs(This)      
      IF( Jump ) THEN
        IF( ASSOCIATED( MortarBC % Diag ) ) THEN
          IF( SIZE( MortarBC % Diag ) < NDofs * Projector % NumberOfRows ) THEN
            DEALLOCATE( MortarBC % Diag ) 
          END IF
        END IF
        IF( .NOT. ASSOCIATED( MortarBC % Diag ) ) THEN
          CALL Info('SetPeriodicBoundariesPass1','Allocating projector mortar diag',Level=10)
          ALLOCATE( MortarBC % Diag( NDofs * Projector % NumberOfRows ) )
          MortarBC % Diag = 0.0_dp
        ELSE
          MortarBC % Diag( DOF::NDOFs ) = 0.0_dp
        END IF

        IF( ASSOCIATED( MortarBC % Rhs ) ) THEN
          IF( SIZE( MortarBC % Rhs ) < NDofs * Projector % NumberOfRows ) THEN
            DEALLOCATE( MortarBC % Rhs ) 
          END IF
        END IF
        IF( .NOT. ASSOCIATED( MortarBC % Rhs ) ) THEN
          CALL Info('SetPeriodicBoundariesPass1','Allocating projector mortar rhs',Level=10)
          ALLOCATE( MortarBC % Rhs( NDofs * Projector % NumberOfRows ) )
          MortarBC % Rhs = 0.0_dp
        ELSE
          MortarBC % Rhs( DOF::NDOFs ) = 0.0_dp
        END IF
      END IF

      ! Create the permutation that is later need in putting the diag and rhs to correct position
      IF( ASSOCIATED( MortarBC % Perm ) ) THEN
        IF( SIZE( MortarBC % Perm ) < SIZE( Perm ) ) THEN
          DEALLOCATE( MortarBC % Perm ) 
        END IF
      END IF
      IF( .NOT. ASSOCIATED( MortarBC % Perm ) ) THEN
        CALL Info('SetPeriodicBoundariesPass1','Allocating projector mortar perm',Level=10)
        ALLOCATE( MortarBC % Perm( SIZE( Perm ) ) )
      END IF
      
      MortarBC % Perm = 0
      DO i=1,SIZE( Projector % InvPerm )
        j = Projector % InvPerm(i) 
        IF( j > 0 .AND. j <= SIZE( Perm ) ) THEN
          MortarBC % Perm( j ) = i
        END IF
      END DO
      
      ! We can use directly the nodal projector
      MortarBC % Projector => Projector
      MortarBC % SlaveScale = -Scale
      MortarBC % MasterScale = -1.0_dp
 
      IF( Jump ) THEN
        PPerm => Perm
        CALL CalculateNodalWeights(Model % Solver,.TRUE.,&
            PPerm,'Periodic BC Coefficient '//Name(1:nlen),WeightVar )
        IF(.NOT. ASSOCIATED( WeightVar ) ) THEN
          CALL Fatal('SetPeriodicBoundariesPass1','Nodal weights needed for setting jumps!')
        END IF
        
        DO i=1,Projector % NumberOfRows
          k = Projector % InvPerm(i)
          IF ( k<=0 ) CYCLE
          
          ! Add the diagonal unity projector (scaled)
          weight = WeightVar % Values( PPerm( k ) )
          coeff = ListGetRealAtNode( BC,'Periodic BC Coefficient '&
              //Name(1:nlen), k, Found )

          ! For Nodal projector the entry is 1/(weight*coeff)
          ! For Galerkin projector the is weight/coeff 
          IF( Found ) THEN
            MortarBC % Diag( NDOFS* (i-1) + DOF ) = 1.0_dp / ( weight * coeff ) 
          END IF
        END DO
      END IF

      Model % Solver % MortarBCsChanged = .TRUE.
      
    ELSE

      ALLOCATE(F)
      F = A
      F % RHS => F % BulkRHS
      F % Values => F % BulkValues

      DO i=1,Projector % NumberOfRows
         ii = Projector % InvPerm(i)
         IF( ii == 0 ) CYCLE
         k = Perm(ii)
         IF ( .NOT. Done(ii) .AND. k>0 ) THEN
            k = NDOFs * (k-1) + DOF
            DO l=Projector % Rows(i),Projector % Rows(i+1)-1
              IF ( Projector % Cols(l) <= 0 .OR. Projector % Values(l)==0.0d0 ) CYCLE

              m = Perm(Projector % Cols(l))
              IF ( m > 0 ) THEN
                m = NDOFs*(m-1) + DOF
                DO nn=A % Rows(k),A % Rows(k+1)-1
                   CALL AddToMatrixElement( A, m, A % Cols(nn), &
                          -scale*Projector % Values(l) * A % Values(nn) ) 
                   IF (ASSOCIATED(F % Values)) THEN
                     CALL AddToMatrixElement( F, m, F % Cols(nn), &
                          -scale*Projector % Values(l) * F % Values(nn) )
                   END IF
                END DO
                b(m)=b(m) - scale*Projector % Values(l)*b(k) 
                IF (ASSOCIATED(F % RHS)) THEN
                  F % RHS(m) = F % RHS(m) - scale*Projector % Values(l)*F % RHS(k)
                END IF
              END IF
            END DO
         END IF
         Done(ii) = .TRUE.
      END DO
      DEALLOCATE(F)
    END IF

!------------------------------------------------------------------------------
   END SUBROUTINE SetPeriodicBoundariesPass1
!------------------------------------------------------------------------------


!> At second pass add the [...1 .. -1 ...] row structure that results to the 
!> equality of the periodic dofs. 
!------------------------------------------------------------------------------
   SUBROUTINE SetPeriodicBoundariesPass2( Model, A, b, &
                      Name, DOF, NDOFs, Perm, This, Done )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model        !< The current model structure
    TYPE(Matrix_t), TARGET :: A   !< The global matrix
    REAL(KIND=dp) :: b(:)         !< The global RHS vector
    CHARACTER(LEN=*) :: Name      !< name of the dof to be set
    LOGICAL :: Done(:)            !< Has the node already been done. 
    INTEGER :: This               !< Number of the current boundary.
    INTEGER :: DOF                !< The order number of the dof
    INTEGER :: NDOFs              !< the total number of DOFs for this equation
    INTEGER, TARGET :: Perm(:)    !< The node reordering info, this has been generated at the
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,m,n,nn,ii,nlen
    LOGICAL :: GotIt, Jump, Found
    REAL(KIND=dp) :: Scale,s,ValueOffset,val,coeff,weight
    TYPE(Matrix_t), POINTER :: Projector
    INTEGER, POINTER :: PPerm(:)
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Variable_t), POINTER :: WeightVar
!------------------------------------------------------------------------------


    BC =>  Model % BCs(This) % Values
    IF ( ListGetLogical( BC, &
         'Periodic BC Use Lagrange Coefficient', GotIt ) ) RETURN

    IF ( ListGetLogical( BC, &
         'Periodic BC Explicit', GotIt ) ) RETURN

    nlen = LEN_TRIM(Name)

    IF ( ListGetLogical( BC, &
       'Periodic BC ' // Name(1:nlen), GotIt ) ) THEN
      IF( ListGetLogical( BC,'Antisymmetric BC',GotIt ) ) THEN
        Scale = 1.0_dp
      ELSE
        Scale = -1.0_dp
      END IF
    ELSE IF ( ListGetLogical( BC, &
        'Anti Periodic BC ' // Name(1:nlen), GotIt ) ) THEN
      Scale = 1.0d0
    ELSE 
      Scale = ListGetCReal( BC, &
          'Periodic BC Scale ' // Name(1:nlen), GotIt) 
      IF(.NOT. GotIt ) RETURN      
    END IF

    ValueOffset = ListGetCReal( BC, &
          'Periodic BC Offset ' // Name(1:nlen), GotIt) 

    Jump = ListCheckPresent( BC, &
        'Periodic BC Coefficient '//Name(1:nlen))
    IF( Jump ) THEN
      PPerm => Perm
      CALL CalculateNodalWeights(Model % Solver,.TRUE.,&
          PPerm,'Periodic BC Coefficient '//Name(1:nlen),WeightVar )
      IF(.NOT. ASSOCIATED( WeightVar ) ) THEN
        CALL Fatal('SetPeriodicBoundariesPass1','Nodal weights needed for setting jumps!')
      END IF
    END IF

    Projector => Model % BCs(This) % PMatrix
    IF ( .NOT. ASSOCIATED(Projector) ) RETURN


!   Do the assembly of the projector:
!   ---------------------------------
    DO i=1,Projector % NumberOfRows
       ii = Projector % InvPerm(i)
       IF( ii == 0 ) CYCLE

       k = Perm( ii )
       IF ( .NOT. Done(ii) .AND. k > 0 ) THEN

         IF( Jump ) THEN
           weight = WeightVar % Values( k )
           coeff = ListGetRealAtNode( BC,'Periodic BC Coefficient '&
               //Name(1:nlen),ii, Found )
           val = -weight * coeff 
           scale = -1.0
         ELSE         
           val = 1.0_dp
         END IF

          k = NDOFs * (k-1) + DOF
          IF(.NOT. Jump) THEN
            CALL ZeroRow( A,k )
            b(k) = 0.0_dp
          END IF

          DO l = Projector % Rows(i), Projector % Rows(i+1)-1
             IF ( Projector % Cols(l) <= 0 ) CYCLE
             m = Perm( Projector % Cols(l) )
             IF ( m > 0 ) THEN
               m = NDOFs * (m-1) + DOF
               CALL AddToMatrixElement( A,k,m,val * Projector % Values(l) )
             END IF
          END DO

          b(k) = b(k) - ValueOffset 
          CALL AddToMatrixElement( A,k,k,scale*val )
          
        END IF
       Done(ii) = .TRUE.
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE SetPeriodicBoundariesPass2
!------------------------------------------------------------------------------




!> At first pass sum together the rows related to the periodic dofs.
!------------------------------------------------------------------------------
   SUBROUTINE SetFrictionBoundaries( Model, A, b, &
                      Name, NDOFs, Perm )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model        !< The current model structure
    TYPE(Matrix_t), TARGET :: A   !< The global matrix
    REAL(KIND=dp) :: b(:)         !< The global RHS vector
    CHARACTER(LEN=*) :: Name      !< name of the dof to be set
    INTEGER :: NDOFs              !< the total number of DOFs for this equation
    INTEGER, TARGET :: Perm(:)    !< The node reordering info, this has been generated at the
                                  !< beginning of the simulation for bandwidth optimization
!------------------------------------------------------------------------------
    INTEGER :: t,u,j,k,k2,l,l2,n,bc_id,nlen,NormalInd
    LOGICAL :: Found
    REAL(KIND=dp),ALLOCATABLE :: Coeff(:)
    LOGICAL, ALLOCATABLE :: NodeDone(:)
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:)
!------------------------------------------------------------------------------


    CALL Info('SetFrictionBoundaries','Setting friction boundaries for variable: '//TRIM(Name),&
        Level=8)

    IF( NDOFs /= 2 ) THEN
      CALL Warn('SetFrictionBoundaries','Assumes friction only in 2D system')
    END IF

    nlen = LEN_TRIM(Name)
    Mesh => Model % Mesh

    ALLOCATE( NodeDone( SIZE( Perm ) ) )
    ALLOCATE( Coeff( Mesh % MaxElementNodes ) )

    NodeDone = .FALSE.
    Coeff = 0.0_dp
    
    DO t = Mesh % NumberOfBulkElements+1, &
        Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(t)
      
      Model % CurrentElement => Element
            
      DO bc_id = 1,Model % NumberOfBCs
        IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc_id) % Tag ) EXIT
      END DO
      IF ( bc_id > CurrentModel % NumberOfBCs ) CYCLE     
      BC => Model % BCs(bc_id) % Values

      IF( .NOT. ListGetLogical( BC,& 
          'Friction BC ' // Name(1:nlen), Found ) ) CYCLE

      NodeIndexes => Element % NodeIndexes
      n = Element % TYPE % NumberOfNodes
      
      Coeff(1:n) = ListGetReal( BC,& 
          'Friction Coefficient ' // Name(1:nlen), n, NodeIndexes )
      IF( ListGetLogical( BC,& 
          'Normal-Tangential ' // Name(1:nlen), Found ) ) THEN
        NormalInd = 1 
      ELSE
        NormalInd = ListGetInteger( BC,& 
            'Friction Normal Component ' // Name(1:nlen) )
      END IF

      DO i = 1, n
        j = Perm( Nodeindexes(i) )
        IF( NodeDone( j ) ) CYCLE

        k = NDOFs * (j-1) + NormalInd
        k2 = NDOFs * (j-1) + ( 3 - NormalInd ) 

        DO l = A % Rows(k),A % Rows(k+1)-1
          DO l2 = A % Rows(k2), A % Rows(k2+1)-1
            IF( A % Cols(l2) == A % Cols(l) ) EXIT
          END DO
          A % Values(l2) = A % Values(l2) - Coeff(i) * A % Values(l)
        END DO
        A % Rhs(k2) = A % Rhs(k2) - Coeff(i) * A % Rhs(k)
        NodeDone( j ) = .TRUE.
      END DO
    END DO

    n = COUNT( NodeDone ) 
    CALL Info('SetFrictionBoundaries','Number of friction nodes: '//I2S(n),Level=10)

    DEALLOCATE( NodeDone, Coeff )

!------------------------------------------------------------------------------
  END SUBROUTINE SetFrictionBoundaries
!------------------------------------------------------------------------------


!> Set the diagonal entry related to mortar BCs.
!> This implements the implicit jump condition. 
!------------------------------------------------------------------------------
   SUBROUTINE SetWeightedProjectorJump( Model, A, b, &
       Name, DOF, NDOFs, Perm )
     !------------------------------------------------------------------------------
     TYPE(Model_t) :: Model        !< The current model structure
     TYPE(Matrix_t), TARGET :: A   !< The global matrix
     REAL(KIND=dp) :: b(:)         !< The global RHS vector
     CHARACTER(LEN=*) :: Name      !< name of the dof to be set
     INTEGER :: DOF                !< The order number of the dof
     INTEGER :: NDOFs              !< the total number of DOFs for this equation
     INTEGER, TARGET :: Perm(:)    !< The node reordering info
     !------------------------------------------------------------------------------
     INTEGER :: i,j,k,i2,j2,k2,n,u,v,node,totsize,nodesize,bc_ind,&
         nnodes,nlen,TargetBC
     INTEGER, POINTER :: PPerm(:)
     INTEGER, ALLOCATABLE :: InvInvPerm(:)
     LOGICAL :: Found, AddRhs, AddCoeff, AddRes
     LOGICAL, ALLOCATABLE :: NodeDone(:)
     REAL(KIND=dp) :: coeff, weight, voff, res
     TYPE(Matrix_t), POINTER :: Projector
     TYPE(ValueList_t), POINTER :: BC
     TYPE(Element_t), POINTER :: Element, Left, Right
     LOGICAL :: SomethingDone
     TYPE(MortarBC_t), POINTER :: MortarBC
     !------------------------------------------------------------------------------

     ! If there is no mortar projector then nothing to do
     SomethingDone = .FALSE.

     ! Go through the projectors and check for jumps
     ! If there is a jump add an entry to the diagonal-to-be
     DO bc_ind=1,Model % NumberOFBCs

       MortarBC => Model % Solver % MortarBCs(bc_ind) 

       Projector => MortarBC % Projector
       IF( .NOT. ASSOCIATED( Projector ) ) CYCLE

       ! For this boundary there should also be a coefficient 
       ! otherwise nothing needs to be done. 
       nlen = LEN_TRIM(Name)
       BC => Model % BCs(bc_ind) % Values

       AddCoeff = ListCheckPresent( BC,'Mortar BC Coefficient '//Name(1:nlen))
       AddRes = ListCheckPresent( BC,'Mortar BC Resistivity '//Name(1:nlen))
       AddRhs = ListCheckPresent( BC,'Mortar BC Offset '//Name(1:nlen))

       IF( .NOT. (AddCoeff .OR. AddRes .OR. AddRhs) ) CYCLE

       Model % Solver % MortarBCsChanged = .TRUE.
       
       IF( .NOT. ASSOCIATED( Projector % InvPerm ) ) THEN
         CALL Fatal('SetWeightedProjectorJump','The > Projector % InvPerm < is really needed here!')
       END IF

       totsize = MAXVAL( Perm )
       nodesize = MAXVAL( Perm(1:Model % Mesh % NumberOfNodes ) )

       IF( AddCoeff .OR. AddRes ) THEN
         IF( ASSOCIATED( MortarBC % Diag ) ) THEN
           IF( SIZE( MortarBC % Diag ) < NDofs * Projector % NumberOfRows ) THEN
             DEALLOCATE( MortarBC % Diag ) 
           END IF
         END IF
         IF( .NOT. ASSOCIATED( MortarBC % Diag ) ) THEN
           CALL Info('SetWeightedProjectorJump','Allocating projector mortar diag',Level=10)
           ALLOCATE( MortarBC % Diag( NDofs * Projector % NumberOfRows ) )
           MortarBC % Diag = 0.0_dp
         ELSE
           MortarBC % Diag(DOF::NDOFs) = 0.0_dp
         END IF
       END IF

       IF( AddRhs ) THEN
         IF( ASSOCIATED( MortarBC % Rhs ) ) THEN
           IF( SIZE( MortarBC % Rhs ) < NDofs * Projector % NumberOfRows ) THEN
             DEALLOCATE( MortarBC % Rhs ) 
           END IF
         END IF
         IF( .NOT. ASSOCIATED( MortarBC % Rhs ) ) THEN
           CALL Info('SetWeightedProjectorJump','Allocating projector mortar rhs',Level=10)
           ALLOCATE( MortarBC % Rhs( NDofs * Projector % NumberOfRows ) )
           MortarBC % Rhs = 0.0_dp
         ELSE
           MortarBC % Rhs(DOF::NDOFs) = 0.0_dp
         END IF
       END IF

       ! Create the permutation that is later need in putting the diag and rhs to correct position
       IF( ASSOCIATED( MortarBC % Perm ) ) THEN
         IF( SIZE( MortarBC % Perm ) < SIZE( Perm ) ) THEN
           DEALLOCATE( MortarBC % Perm ) 
         END IF
       END IF
       IF( .NOT. ASSOCIATED( MortarBC % Perm ) ) THEN
         CALL Info('SetWeightedProjectorJump','Allocating projector mortar perm',Level=10)
         ALLOCATE( MortarBC % Perm( SIZE( Perm ) ) )
       END IF

       MortarBC % Perm = 0
       DO i=1,SIZE( Projector % InvPerm )
         j = Projector % InvPerm(i) 
         IF( j > 0 .AND. j <= nodesize ) THEN
           MortarBC % Perm( j ) = i
         END IF
       END DO


       TargetBC = ListGetInteger( BC,'Mortar BC',Found ) 

       CALL Info('SetWeightedProjectorJump','Setting jump to mortar projector in BC '&
           //I2S(bc_ind),Level=7)
    
       ! Create a table that shows how the additional degrees of freedom map
       ! to their corresponding regular dof. This is needed when creating the jump.
       ALLOCATE( NodeDone( Projector % NumberOfRows ) )
       NodeDone = .FALSE.
       
       ! Looping through elements rather than looping through projector rows directly
       ! is done in order to be able to refer to boundary properties associated 
       ! with the element. 
       DO t=1,Model % Mesh % NumberOfBoundaryElements
         Element => Model % Mesh % Elements( t + Model % Mesh % NumberOfBulkElements )

         IF( Element % BoundaryInfo % Constraint /= Model % BCs(bc_ind) % Tag ) CYCLE

         ! Outside code this tells the active element
         Model % CurrentElement => Element

         Left => Element % BoundaryInfo % Left
         Right => Element % BoundaryInfo % Right 
        
         IF( TargetBC > 0 ) THEN
           IF( ASSOCIATED( Left ) ) THEN
             IF( Left % PartIndex /= ParEnv % myPE ) CYCLE
           ELSE IF ( ASSOCIATED( Right ) ) THEN
             IF( Left % PartIndex /= ParEnv % myPE ) CYCLE
           ELSE
             CYCLE
           END IF
         ELSE
           ! This case is for the case when TargetBC = 0 i.e. for Discontinuous BC
           ! These are conditions that resulted to creation of zero 
           ! constraint matrix entries in this partition so no need to do them.
           IF(.NOT. ASSOCIATED( Left ) .OR. .NOT. ASSOCIATED( Right ) ) THEN
             CYCLE
           END IF
           
           ! For this we have a zero mass matrix entry so don't bother to add zero
!          IF( Left % PartIndex /= ParEnv % myPE .AND. &
!              Right % PartIndex /= ParEnv % myPe ) THEN
!            CYCLE
!          END IF
         END IF

         nnodes = Element % TYPE % NumberOfNodes
         DO u=1, nnodes
           node = Element % NodeIndexes(u)

           IF( Perm( node ) == 0 ) CYCLE

           i = MortarBC % Perm( node ) 
           IF( i == 0 ) CYCLE

           IF( NodeDone( i ) ) CYCLE
           NodeDone( i ) = .TRUE. 

           Found = .FALSE.

           IF( AddCoeff ) THEN
             coeff = ListGetRealAtNode( BC,'Mortar BC Coefficient '&
                 //Name(1:nlen),node, Found )        
             res = 1.0_dp / coeff
           END IF

           IF( AddRes ) THEN
             res = ListGetRealAtNode( BC,'Mortar BC Resistivity '&
                 //Name(1:nlen),node, Found )
           END IF

           ! For Nodal projector the entry is 1/(weight*coeff)
           ! For Galerkin projector the is weight/coeff 
           IF( Found ) THEN 
             IF( AddCoeff .OR. Addres ) THEN
               MortarBC % Diag(NDOFs*(i-1)+DOF) = res
             END IF
           END IF

           IF( AddRhs ) THEN
             voff = ListGetRealAtNode( BC,'Mortar BC Offset '&
                 //Name(1:nlen),node, Found )        
             IF( Found ) THEN
               MortarBC % Rhs(NDofs*(i-1)+DOF) = voff
             END IF
           END IF

         END DO
       END DO
       
       SomethingDone = .TRUE.

       DEALLOCATE( NodeDone )
     END DO

     IF( SomethingDone ) THEN
       CALL Info('setWeightedProjectorJump','Created a jump for weighted projector',Level=7)
     END IF
 
!------------------------------------------------------------------------------
   END SUBROUTINE SetWeightedProjectorJump
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  END SUBROUTINE SetDirichletBoundaries
!------------------------------------------------------------------------------



  
!------------------------------------------------------------------------------
!> Prepare to set Dirichlet conditions for attachment DOFs in the case of
!> component mode synthesis
!------------------------------------------------------------------------------
  SUBROUTINE SetConstraintModesBoundaries( Model, Solver, A, b, &
      Name, NDOFs, Perm )
    !------------------------------------------------------------------------------
    TYPE(Model_t) :: Model              !< current model structure
    TYPE(Solver_t), TARGET :: Solver    !< current solver structure 
    TYPE(Matrix_t), TARGET :: A         !< global matrix
    REAL(KIND=dp) :: b(:)               !< global RHS vector
    CHARACTER(LEN=*) :: Name            !< name of the dof to be set
    INTEGER :: NDOFs                    !< the total number of DOFs for this equation
    INTEGER, TARGET :: Perm(:)          !< The node reordering info, this has been generated at the
                                        !< beginning of the simulation for bandwidth optimization
!------------------------------------------------------------------------------
    INTEGER :: i,t,t1,t2,u,j,k,k2,l,l2,n,ent_id,nlen,NormalInd,Ncomplex
    LOGICAL :: Found, IgnoreP, HaveP, Passive, Parallel, DoIt
    TYPE(Solver_t), POINTER :: pSolver
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Element
    TYPE(Variable_t), POINTER :: Var
    INTEGER :: Indexes(50),nb,poffset
    INTEGER, POINTER :: NodeIndexes(:), MasterBodies(:)
    INTEGER, ALLOCATABLE :: BCPerm(:)
    INTEGER :: NoModes, NoEntities
    LOGICAL :: ExternalLoop, BFMode, BCMode, LumpedMode, CompMode, &
        EmWaveMode, RhsMode, CoilMode, ComplexMode, SingleMode, Stat, GotDir
    REAL(KIND=dp) :: ModeDir(NDofs)
    REAL(KIND=dp), POINTER :: HelperArray(:,:)
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    CHARACTER(*), PARAMETER :: Caller = 'SetConstraintModesBoundaries'

!------------------------------------------------------------------------------

    CALL Info(Caller,'Checking constraint modes boundaries!',Level=20)
    
    nlen = LEN_TRIM(Name)
    Mesh => Solver % Mesh
    Var => Solver % Variable     

    ! Initially this is -1 and and hence the 2nd call is fast if no modes are present
    IF( Solver % NumberOfConstraintModes == 0 ) RETURN

    ! Check in which list we have the constraint modes set
    BCMode = ListCheckPrefixAnyBC( Model,'Constraint Mode')
    BFMode = ListCheckPrefixAnyBodyForce( Model,'Constraint Mode')
    CompMode = ListCheckPrefixAnyComponent( Model,'Constraint Mode')

    IF( .NOT. (BCMode .OR. BFMode .OR. CompMode) ) THEN
      Solver % NumberOfConstraintModes = 0
      RETURN
    END IF
    
    ExternalLoop = ListGetLogical( Solver % Values,'Nonlinear System Constraint Modes', Found ) .OR. &
        ListGetLogical( Solver % Values,'Steady State Constraint Modes', Found ) .OR. &
        ListGetLogical( Solver % Values,'Run Control Constraint Modes', Found ) .OR. &
        ListGetLogical( Model % Control,'Constraint Modes Analysis', Found )
       
    EmWaveMode = ListGetLogical( Solver % Values,'Constraint Modes EM Wave',Found )
    CoilMode = ListGetLogical( Solver % Values,'Constraint Modes Coils',Found ) 
    RhsMode = ListGetLogical(Solver % Values,'Constraint Modes Rhs',Found ) 
    LumpedMode = ListGetLogical(Solver % Values,'Constraint Modes Lumped',Found )
    
    ! These work on the rhs vector, not Dirichlet values.
    RhsMode = RhsMode .OR. EmWaveMode .OR. CoilMode
    LumpedMode = LumpedMode .OR. EmWaveMode .OR. CoilMode
    
    Element => Mesh % Elements(1)
    pSolver => Solver
    HaveP = isActivePElement(Element,pSolver)
    
    IgnoreP = .FALSE.
    IF( HaveP ) THEN
      IgnoreP = ListGetLogical( Solver % Values,'Ignore Constraint Modes p',Found )  
    END IF      

    CALL Info(Caller,'Setting constraint modes boundaries for variable: '&
        //TRIM(Name),Level=7)

    Parallel = ( ParEnv % PEs > 1 ) .AND. ( .NOT. Mesh % SingleMesh ) 
    
    ! Allocate the indeces for the constraint modes.
    ! We may be revisiting the routine, and the mesh may have changed...
    IF( ASSOCIATED( Var % ConstraintModesIndeces ) ) THEN
      IF( SIZE( Var % ConstraintModesIndeces ) /= A % NumberOfRows ) THEN
        DEALLOCATE( Var % ConstraintModesIndeces )
        IF( RhsMode ) THEN
          IF(ASSOCIATED( Var % ConstraintModesWeights ) ) THEN
            DEALLOCATE( Var % ConstraintModesWeights )
          END IF
        END IF
      END IF
      DoIt = ListGetLogical( Solver % Values,'Constraint Modes Fluxes Norm', Found )
      IF(.NOT. Found) DoIt = ListGetLogical( Solver % Values,'Constraint Modes Matrix Norm', Found )
      IF(DoIt) CALL ListAddConstReal( CurrentModel % Simulation,&
          'res: '//TRIM(Solver % Variable % Name)//' lumped matrix norm',0.0_dp)
    END IF
    IF(.NOT. ASSOCIATED( Var % ConstraintModesIndeces ) ) THEN
      ALLOCATE( Var % ConstraintModesIndeces( A % NumberOfRows ) )
      IF( RhsMode ) THEN
        ALLOCATE( Var % ConstraintModesWeights( A % NumberOfRows ) )
      END IF
    END IF
    
    Var % ConstraintModesIndeces = 0
    IF( RhsMode ) Var % ConstraintModesWeights = 0.0_dp
    
    IF( BCMode ) THEN
      NoEntities = Model % NumberOfBCs
    ELSE IF( BFMode ) THEN
      NoEntities = Model % NumberOfBodyForces
    ELSE IF(CompMode) THEN
      NoEntities = Model % NumberOfComponents
    ELSE
      CALL Fatal(Caller,'Uknown list for constraint modes!')
    END IF
    CALL Info(Caller,'Number of list entities to check: '//I2S(NoEntities),Level=20)
    
    ALLOCATE( BCPerm( NoEntities ) )
    BCPerm = 0    
    j = 0

    DO ent_id = 1,NoEntities
      IF( BCMode ) THEN
        BC => Model % BCs(ent_id) % Values
      ELSE IF(BFMode ) THEN
        BC => Model % BodyForces(ent_id) % Values
      ELSE IF(CompMode) THEN
        BC => Model % Components(ent_id) % Values
      END IF
        
      k = ListGetInteger( BC,'Constraint Mode', Found )
      IF(.NOT. Found ) THEN
        k = ListGetInteger( BC,&
            'Constraint Mode '// Name(1:nlen), Found )
      END IF
      IF( Found ) THEN
        IF( k == 0 ) k = -1  ! Ground gets negative value
        BCPerm(ent_id) = k        
      ELSE
        DoIt = ListGetLogical( BC,'Constraint Modes', Found )
        IF(.NOT. Found ) THEN
          DoIt = ListGetLogical( BC,'Constraint Modes ' // Name(1:nlen), Found ) 
        END IF
        IF(DoIt) THEN
          j = j + 1
          BCPerm(ent_id) = j
        END IF
      END IF
    END DO
    
    j = MAXVAL( BCPerm )
    CALL Info(Caller,'Number of active constraint modes boundaries: '&
        //I2S(j),Level=7)
    IF( j == 0 ) THEN
      CALL Fatal(Caller,&
          'Constraint Modes Analysis requested but no constrained BCs given!')
    END IF

    
    ComplexMode = LIstGetLogical( Solver % Values,'Linear System Complex',Found)
    IF( ComplexMode ) THEN
      Ncomplex = 2
      CALL Info(Caller,'Assuming complex valued system for constraint modes',Level=12)
    ELSE
      Ncomplex = 1
      CALL Info(Caller,'Assuming real valued system for constraint modes',Level=12)
    END IF

    SingleMode = ListGetLogical( Solver % Values,'Constraint Modes Single',Found ) 
    IF(.NOT. Found ) THEN
      SingleMode = ListCheckPresentAnyBC( Model,'Constraint Mode Direction' )            
    END IF
    IF(SingleMode) THEN
      CALL Info(Caller,'Setting constraint modes for all components at once!',Level=12)
    END IF
    
    IF( SingleMode ) THEN
      NoModes = j
    ELSE      
      NoModes = NDOFS * j  / Ncomplex
    END IF
      
    IF( BcMode ) THEN
      t1 = Mesh % NumberOfBulkElements+1
      t2 = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
    ELSE
      t1 = 1
      t2 = Mesh % NumberOfBulkElements
    END IF
    
    DO t = t1, t2
      Element => Mesh % Elements(t)      
      
      ModeDir = 1.0_dp
      IF( BCMode ) THEN
        DO ent_id = 1,Model % NumberOfBCs
          IF ( Element % BoundaryInfo % Constraint == Model % BCs(ent_id) % Tag ) EXIT
        END DO
        IF( ent_id > Model % NumberOfBCs ) CYCLE        
        HelperArray => ListGetConstRealArray( Model % BCs(ent_id) % Values,'Constraint Mode Direction',Found)
        IF(Found) ModeDir(1:ndofs) = HelperArray(1:ndofs,1)
      ELSE IF( BFMode ) THEN
        ent_id = ListGetInteger( Model % Bodies( Element % BodyId ) % Values,'Body Force',Found)
        IF( ent_id == 0) CYCLE
        HelperArray => ListGetConstRealArray( Model % BodyForces(ent_id) % Values,'Constraint Mode Direction',Found)
        IF(Found) ModeDir(1:ndofs) = HelperArray(1:ndofs,1)
      ELSE        
        DO ent_id=1,Model % NumberOfComponents
          IF(BCPerm(ent_id) == 0) CYCLE
          MasterBodies => ListGetIntegerArray( Model % Components(ent_id) % Values,'Master Bodies',Found)
          IF(.NOT. Found) CYCLE
          IF( ANY( MasterBodies == Element % BodyId ) ) EXIT          
        END DO
        IF(ent_id > Model % NumberOfComponents ) CYCLE
      END IF

      ! Look-up table for quick determiation
      IF( BCPerm(ent_id) == 0 ) CYCLE

      n = Element % Type % NumberOfNodes

      ! This is used in standard setting of Dirichlet BCs
      nb = mGetElementDOFs( Indexes, Element, Solver )
            
      IF( RhsMode ) THEN
        BLOCK 
          REAL(KIND=dp) :: Basis(n), Weight(n), detJ, s          
          IP = GaussPoints( Element )
          CALL CopyElementNodesFromMesh( Nodes, Solver % Mesh, n, Element % NodeIndexes)
                    
          DO i=1,IP % n 
            stat = ElementInfo( Element, Nodes, IP % U(i), IP % V(i), &
              IP % W(i), detJ, Basis )          
            s = IP % s(i) * DetJ
            DO k=1,NDOFS                            
              Var % ConstraintModesWeights(NDOfs*(Perm(Indexes(1:n))-1)+k) = &
                  Var % ConstraintModesWeights(NDOfs*(Perm(Indexes(1:n))-1)+k) + ModeDir(k) * s * Basis(1:n)
            END DO
          END DO
        END BLOCK
      END IF
            
      ! If for some reason we do not want to set the P dofs to zero
      IF(IgnoreP) nb = n

      !PRINT *,'BCPerm:',MINVAL(BCPerm), MAXVAL(BCPerm), COUNT(BCPerm>0)
      
      ! For vector valued problems treat each component as separate dof
      DO k=1,NDOFs       
        IF( SingleMode ) THEN
          j = BCPerm(ent_id)
        ELSE
          j = NDOFS*(BCPerm(ent_id)-1)+k
        END IF
        DO l=1,nb
          ! The index to constrain
          IF( Perm(Indexes(l)) == 0 ) CYCLE
          l2 = NDOFS*(Perm(Indexes(l))-1)+k         
          Var % ConstraintModesIndeces(l2) = j
        END DO
      END DO
      !PRINT *,'Single range:',SingleMode,MAXVAL(Var % ConstraintModesIndeces), COUNT(Var % ConstraintModesIndeces > 0)
    END DO
    DEALLOCATE(BCPerm)

    ! Some single node or edge could stretch to the surface even though it is not
    ! part of any boundary element in parallel. Hence we need to communicate the tag. 
    IF( Parallel ) THEN
      CALL Info(Caller,'Communicating tags for constraint modes',Level=20)
      CALL CommunicateParallelSystemTag(A % ParallelInfo,Itag = Var % ConstraintModesIndeces,ParOper=2)
    END IF
      
    ! Set the p dofs to negative since we don't ever want to set them to one!
    ! Note that there are some dofs related to ground that are already negative.
    ! Hence ground and p-pubbles are treated alike. 
    IF(HaveP .AND. .NOT. IgnoreP) THEN
      poffset = 2*(NoModes + 1)
      DO l=Mesh % NumberOfNodes+1, SIZE(Perm)
        j = Perm(l)
        IF(j==0) CYCLE
        DO k=1,NDOFs       
          l2 = NDOFS*(j-1)+k                 
          
          ! Subtract a big enough number of the constraint modes so that they are always negative. 
          IF( Var % ConstraintModesIndeces(l2) /= 0 ) THEN
            Var % ConstraintModesIndeces(l2) = Var % ConstraintModesIndeces(l2) - poffset 
          END IF
        END DO
      END DO
    END IF
    
    ! This is just for information.
    ! We show how Dirichlet conditions are split among nodal and p dofs, and ground. 
    IF( InfoActive(12) ) THEN
      DO i=0,NoModes
        j = i
        IF(i==0) j=-1

        k = COUNT( Var % ConstraintModesIndeces == j )
        IF(k>0) THEN
          CALL Info('SetConstaintModesBoundaries',&
              'Mode '//I2S(i)//' has '//I2S(k)//' dofs')
          IF(.NOT. IgnoreP .AND. HaveP ) THEN
            j = j - poffset
            k = COUNT( Var % ConstraintModesIndeces == j )
            IF(k>0) THEN
              CALL Info('SetConstaintModesBoundaries',&
                  'Mode '//I2S(i)//' has additional '//I2S(k)//' p-dofs')
            END IF
          END IF
        END IF
      END DO
    END IF

    
    ! The constraint modes can be either lumped or not.
    ! If they are not lumped then mark each individually
    IF( .NOT. LumpedMode ) THEN
      j = 0
      DO i=1,A % NumberOfRows
        IF( Var % ConstraintModesIndeces(i) > 0 ) THEN
          j = j + 1
          Var % ConstraintModesIndeces(i) = j
        END IF
      END DO
      CALL Info(Caller,'Number of active constraint modes: '&
          //I2S(j),Level=7)
      NoModes = j 
    END IF
        
    ! Manipulate the boundaries such that we need to modify only the r.h.s. in the actual linear solver
    ! Do not manipulate if we are setting fluxes!
    IF( .NOT. (ComplexMode .OR. RhsMode ) ) THEN
      WHERE( Var % ConstraintModesIndeces /= 0 ) 
        A % ConstrainedDOF = .TRUE.
        A % DValues = 0.0_dp
      END WHERE
    END IF
    
    Solver % NumberOfConstraintModes = NoModes

    ! We may want to save the results for postprocessing even when we do this in
    ! one sweep. If we iterate over nonlinear, steady state of run control then automatically
    ! we can have access to all components in saving. 
    IF(.NOT. ExternalLoop ) THEN
      Var % NumberOfConstraintModes = NoModes
      ! This routine is visited on every solver call, so the previous table has
      ! to go first -- as ConstraintModesIndeces and ConstraintModesWeights above
      ! already do. Without this the whole NoModes x NumberOfRows table was
      ! dropped unfreed once per timestep.
      IF( ASSOCIATED( Var % ConstraintModes ) ) DEALLOCATE( Var % ConstraintModes )
      ALLOCATE( Var % ConstraintModes( Var % NumberOfConstraintModes, A % NumberOfRows ) )
      Var % ConstraintModes = 0.0_dp
    END IF
          
    CALL Info(Caller,'All done',Level=10)

!------------------------------------------------------------------------------
  END SUBROUTINE SetConstraintModesBoundaries
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Sets just one Dirichlet point in contrast to setting the whole field.
!> This is a lower order routine that the previous one. 
!------------------------------------------------------------------------------
  SUBROUTINE SetDirichletPoint( A, b,DOF, NDOFs, Perm, NodeIndex, NodeValue) 
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Matrix_t), TARGET :: A
    REAL(KIND=dp) :: b(:)
    REAL(KIND=dp) :: NodeValue
    INTEGER :: DOF, NDOFs, Perm(:), NodeIndex
!------------------------------------------------------------------------------

    REAL(KIND=dp) :: s
    INTEGER :: PermIndex

!------------------------------------------------------------------------------
    PermIndex = Perm(NodeIndex)
    IF ( PermIndex > 0 ) THEN
      PermIndex = NDOFs * (PermIndex-1) + DOF

      ! Allocate on demand, as UpdateDirichletDof above does. These arrays are
      ! otherwise created by the higher level Dirichlet routines, so a caller
      ! reaching this one before DefaultDirichletBCs has run -- which is exactly
      ! when setting a single point is useful -- would write through unallocated
      ! arrays and segfault.
      IF(.NOT. ALLOCATED(A % ConstrainedDOF)) THEN
        ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
        A % ConstrainedDOF = .FALSE.
      END IF
      IF(.NOT. ALLOCATED(A % Dvalues)) THEN
        ALLOCATE(A % Dvalues(A % NumberOfRows))
        A % Dvalues = 0._dp
      END IF

      A % ConstrainedDOF(PermIndex) = .TRUE.
      A % DValues(PermIndex) = NodeValue
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE SetDirichletPoint
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------------
  SUBROUTINE CommunicateDirichletBCs(A)
  !-------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A

     REAL(KIND=dp), ALLOCATABLE :: d_e(:,:), g_e(:)
     INTEGER, ALLOCATABLE :: s_e(:,:), r_e(:), fneigh(:), ineigh(:)
     INTEGER :: i,j,k,l,n,nn,ii(ParEnv % PEs), ierr, status(MPI_STATUS_SIZE)

     IF( ParEnv % PEs<=1 ) RETURN
     IF( A % ParallelInfo % NothingShared ) RETURN
     
     ALLOCATE( fneigh(ParEnv % PEs), ineigh(ParEnv % PEs) )

     nn = 0
     ineigh = 0
     DO i=0, ParEnv % PEs-1
       k = i+1
       IF(.NOT.ParEnv % Active(k) ) CYCLE
       IF(i==ParEnv % myPE) CYCLE
       IF(.NOT.ParEnv % IsNeighbour(k) ) CYCLE
       nn = nn + 1
       fneigh(nn) = k
       ineigh(k) = nn
     END DO

     n = COUNT(A % ConstrainedDOF .AND. A % ParallelInfo % GInterface)
     ALLOCATE( s_e(n, nn ), r_e(n) )
     ALLOCATE( d_e(n, nn ), g_e(n) )

     CALL CheckBuffer( nn*3*n )

     ii = 0
     DO i=1, A % NumberOfRows
       IF(A % ConstrainedDOF(i) .AND. A % ParallelInfo % GInterface(i) ) THEN
          DO j=1,SIZE(A % ParallelInfo % Neighbourlist(i) % Neighbours)
            k = A % ParallelInfo % Neighbourlist(i) % Neighbours(j)
            IF ( k == ParEnv % MyPE ) CYCLE
            k = k + 1
            k = ineigh(k)
            IF ( k> 0) THEN
              ii(k) = ii(k) + 1
              d_e(ii(k),k) = A % DValues(i)
              s_e(ii(k),k) = A % ParallelInfo % GlobalDOFs(i)
            END IF
          END DO
       END IF
     END DO

     DO i=1, nn
       j = fneigh(i) 

       CALL MPI_BSEND( ii(i),1,MPI_INTEGER,j-1,110,ELMER_COMM_WORLD,ierr )
       IF( ii(i) > 0 ) THEN
         CALL MPI_BSEND( s_e(1:ii(i),i),ii(i),MPI_INTEGER,j-1,111,ELMER_COMM_WORLD,ierr )
         CALL MPI_BSEND( d_e(1:ii(i),i),ii(i),MPI_DOUBLE_PRECISION,j-1,112,ELMER_COMM_WORLD,ierr )
       END IF
     END DO

     DO i=1, nn
       j = fneigh(i)
       CALL MPI_RECV( n,1,MPI_INTEGER,j-1,110,ELMER_COMM_WORLD, status,ierr )
       IF ( n>0 ) THEN
         IF( n>SIZE(r_e)) THEN
           DEALLOCATE(r_e,g_e)
           ALLOCATE(r_e(n),g_e(n))
         END IF

         CALL MPI_RECV( r_e,n,MPI_INTEGER,j-1,111,ELMER_COMM_WORLD,status,ierr )
         CALL MPI_RECV( g_e,n,MPI_DOUBLE_PRECISION,j-1,112,ELMER_COMM_WORLD, status,ierr )
         DO j=1,n
           k = SearchNode( A % ParallelInfo, r_e(j), Order=A % ParallelInfo % Gorder )
           IF ( k>0 ) THEN
             IF(.NOT. A % ConstrainedDOF(k)) THEN
               A % Dvalues(k) = g_e(j)
               A % ConstrainedDOF(k) = .TRUE.
             END IF
           END IF
         END DO
       END IF
     END DO
     DEALLOCATE(s_e, r_e, d_e, g_e)
  !-------------------------------------------------------------------------------
  END SUBROUTINE CommunicateDirichletBCs
  !-------------------------------------------------------------------------------

  
!-------------------------------------------------------------------------------
  SUBROUTINE EnforceDirichletConditions( Solver, A, b, OffDiagonal ) 
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t), TARGET :: A
    REAL(KIND=dp) :: b(:)
    LOGICAL, OPTIONAL :: OffDiagonal
    
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: ScaleSystem, DirichletComm, Found, NoDiag
    REAL(KIND=dp) :: dval, s
    INTEGER :: i,j,k,n,n2,n3
    CHARACTER(*), PARAMETER :: Caller = 'EnforceDirichletConditions'
    LOGICAL :: Parallel, DoDiagScale
    
    
    Params => Solver % Values

    IF(.NOT. ALLOCATED( A % ConstrainedDOF ) ) THEN
      CALL Info(Caller,'ConstrainedDOF not associated, returning...',Level=8)
      RETURN
    END IF

    Parallel = ( ParEnv % PEs > 1 ) .AND. ( .NOT. Solver % Mesh % SingleMesh ) 

    n = COUNT( A % ConstrainedDOF )
    IF( Parallel ) n = ParallelReduction( n )
      
    IF( n == 0 ) THEN
      CALL Info(Caller,'No Dirichlet conditions to enforce, exiting!',Level=10)
      RETURN
    ELSE
      CALL Info(Caller,'Enforcing total of '//I2S(n)//' Dirichlet conditions.',Level=10)
    END IF    
        
    ! Communicate the Dirichlet conditions for parallel cases since there may be orphans      
    IF ( Parallel ) THEN
      DirichletComm = ListGetLogical( CurrentModel % Simulation, 'Dirichlet Comm', Found)
      IF(.NOT. Found) DirichletComm = .TRUE.
      IF( DirichletComm) THEN
        CALL Info('EnforceDirichletConditions','Communicating Dirichlet conditions',Level=20)
        CALL CommunicateDirichletBCs(A)
      END IF
    END IF
    
    IF( PRESENT( OffDiagonal ) ) THEN
      NoDiag = OffDiagonal
    ELSE
      NoDiag = .FALSE.
    END IF

    IF( NoDiag ) THEN
      ScaleSystem = .FALSE.
    ELSE
      ScaleSystem = ListGetLogical(Params,'Linear System Dirichlet Scaling',Found)
      IF(.NOT.Found) THEN
        ScaleSystem = ListGetLogical(Params,'Linear System Scaling',Found)
        IF(.NOT.Found) ScaleSystem=.TRUE.
      END IF
    END IF

    DoDiagScale = .FALSE.
    IF( ScaleSystem ) THEN
      CALL Info(Caller,'Applying Dirichlet conditions using scaled diagonal',Level=8)
      CALL ScaleLinearSystem(Solver,A,b,ApplyScaling=.FALSE.) !,scalingStr='diagonal')
      DoDiagScale = ASSOCIATED(A % DiagScaling)
    END IF
    
    n2 = 0
    IF( Parallel ) THEN
      n = 0
      DO i=1,A % NumberOfRows
        IF( A % ConstrainedDOF(i) ) THEN
          IF( A % ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % MyPE ) THEN
            n = n + 1
          ELSE
            n2 = n2 + 1
          END IF
        END IF
      END DO
      n = ParallelReduction( n )
      n2 = ParallelReduction( n2 ) 
    END IF
        
    ! Eliminate all entries in matrix that may be eliminated in one sweep
    ! If this is an offdiagonal entry this cannot be done.
    ! Also, if we want to do swap the Dirichlet conditions without
    ! rebuilding the matrix (as in capacitance matrix computation) this
    ! cannot be done. 
    IF ( A % Symmetric .AND. .NOT. NoDiag ) CALL CRS_ElimSymmDirichlet(A,b)
    
    DO k=1,A % NumberOfRows

      IF ( A % ConstrainedDOF(k) ) THEN
        
        dval = A % Dvalues(k) 

        s = 1.0_dp
        IF( ScaleSystem ) THEN
          SELECT CASE(A % ScalingMethod ) 
          CASE( 1 )
            s = A % DiagScaling(k)            
            IF( ABS(s) <= TINY(s) ) s = 1.0_dp
            s = 1._dp / s**2
          CASE( 2 )
            s = 1.0_dp / A % DiagScaling(k)
          CASE( 3 ) 
            s = A % AveScaling
          END SELECT
        END IF
          
        CALL ZeroRow(A, k)

        ! Off-diagonal entries for a block matrix are neglected since the code will
        ! also go through the diagonal entries where the r.h.s. target value will be set.
        IF(.NOT. NoDiag ) THEN
#if 0
          ! It is not clear whether we want to set the Dirichlet conditions only
          ! for the owner partition, or to all partitions. The latter may perhaps
          ! help when using partition-specific preconditioners...
          IF( Parallel ) THEN
            IF( SIZE( A % ParallelInfo % NeighbourList(k) % Neighbours) > 1 ) THEN
              IF( A % ParallelInfo % NeighbourList(k) % Neighbours(1) /= ParEnv % MyPe ) THEN
                b(k) = 0.0_dp
                CYCLE
              END IF
            END IF
          END IF
#endif          
          CALL SetMatrixElement(A,k,k,s)
          b(k) = s * dval
        END IF

      END IF
    END DO

    ! Deallocate scaling since otherwise it could be misused out of context
    IF (DoDiagScale) DEALLOCATE( A % DiagScaling ) 
        
    CALL Info(Caller,'Dirichlet conditions enforced for dofs: '//I2S(n), Level=6)
    IF(n2>0) CALL Info(Caller,'Dirichlet conditions shared count: '//I2S(n2), Level=12)
    
  END SUBROUTINE EnforceDirichletConditions
!-------------------------------------------------------------------------------
   

!------------------------------------------------------------------------------
!> Check if Normal / Tangential vector boundary conditions present and
!> allocate space for normals, and if in 3D for two tangent direction
!> vectors.
!------------------------------------------------------------------------------
   SUBROUTINE CheckNormalTangentialBoundary( Model, VariableName, &
       NumberOfBoundaryNodes, BoundaryReorder, BoundaryNormals,     &
       BoundaryTangent1, BoundaryTangent2, dim )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    CHARACTER(LEN=*) :: VariableName
    INTEGER, POINTER :: BoundaryReorder(:)
    INTEGER :: NumberOfBoundaryNodes,dim
    REAL(KIND=dp), POINTER :: BoundaryNormals(:,:),BoundaryTangent1(:,:), &
        BoundaryTangent2(:,:)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER :: i,j,k,n,np,t,ierr,iter, proc
    LOGICAL :: GotIt, Found, Conditional, Rotational
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: Indexes(:)
    REAL(KIND=dp), ALLOCATABLE :: Condition(:)
    TYPE buff_t
      INTEGER, ALLOCATABLE :: buff(:)
    END TYPE buff_t
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    INTEGER, TARGET :: pIndexes(12)
    INTEGER, POINTER :: nlist(:)
    TYPE(Buff_t), ALLOCATABLE, TARGET :: n_index(:)
    INTEGER, ALLOCATABLE :: n_count(:), gbuff(:)
    TYPE(Variable_t), POINTER :: DispVar
    LOGICAL :: pDisp
    REAL(KIND=dp) :: rad
    CHARACTER(*), PARAMETER :: Caller = 'CheckNormalTangentialBoundary'
!------------------------------------------------------------------------------

    ! need an early initialization to average normals across partitions:
    !-------------------------------------------------------------------
    IF ( Parenv  % PEs >1 ) THEN
      IF (.NOT. ASSOCIATED(Model % Solver % Matrix % ParMatrix) ) &
         CALL ParallelInitMatrix( Model % Solver, Model % Solver % Matrix )
    END IF

    NumberOfBoundaryNodes = 0

    Found = .FALSE.
    DO i=1,Model % NumberOfBCs
      IF ( ListGetLogical(Model % BCs(i) % Values, VariableName, Gotit) ) THEN
        Found = ListGetLogical( Model % BCs(i) % Values, &
           TRIM(VariableName) // ' Rotate',Gotit )
        IF (.NOT. Gotit ) Found = .TRUE.
        IF ( Found ) EXIT
      END IF
    END DO
    IF ( .NOT. Found ) RETURN

    Mesh => Model % Mesh
    n = Mesh % NumberOFNodes

    pDisp = .FALSE.
    NULLIFY( DispVar )
    DO i=1,Model % NumberOfSolvers
      IF( ListGetLogical( Model % Solvers(i) % Values,'Use p Normals',GotIt ) ) THEN
        DispVar => Model % Solvers(i) % Variable
        pDisp = .TRUE.
        IF( SIZE( DispVar % Perm ) > n ) THEN
          n = SIZE( DispVar % Perm )        
        END IF
        EXIT
      END IF
    END DO
            
    IF ( ASSOCIATED( BoundaryReorder ) ) THEN
      IF ( SIZE(BoundaryReorder) < n ) DEALLOCATE( BoundaryReorder )
    END IF
    
    IF ( .NOT. ASSOCIATED( BoundaryReorder ) ) THEN
      CALL Info( Caller,'Allocating BoundaryOrder of size: '//I2S(n),Level=12)
      IF( pDisp ) THEN
        CALL Info(Caller,'Creating normals for p-dofs as well!',Level=12)
      END IF
      ALLOCATE( BoundaryReorder(n) )
    END IF
    
    BoundaryReorder = 0
    
!------------------------------------------------------------------------------
    DO t=Mesh % NumberOfBulkElements + 1, Mesh % NumberOfBulkElements + &
                  Mesh % NumberOfBoundaryElements

      Element => Model % Elements(t)
      IF ( Element % TYPE % ElementCode == 101 )  CYCLE

      Indexes => Element % NodeIndexes
      n = Element % TYPE % NumberOfNodes      
      ALLOCATE( Condition(n)  )
      
      DO i=1,Model % NumberOfBCs
        IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN
          IF ( ListGetLogical( Model % BCs(i) % Values,VariableName, gotIt) ) THEN
            Found = ListGetLogical( Model % BCs(i) % Values, &
                TRIM(VariableName) // ' Rotate',gotIt)
            IF ( Found .OR. .NOT. GotIt ) THEN
              Condition(1:n) = ListGetReal( Model % BCs(i) % Values, &
                  TRIM(VariableName) // ' Condition', n, Indexes, Conditional )
              Rotational = ListGetLogical( Model % BCs(i) % Values,'Rotational Normals',GotIt)
                            
              DO j=1,n
                IF ( Conditional .AND. Condition(j)<0._dp ) CYCLE                
                k = Indexes(j)

                ! If we are using rotational normals then avoid origin.
                IF(Rotational) THEN
                  rad = SQRT(Mesh % Nodes % x(k)**2 + Mesh % Nodes % y(k)**2)
                  IF(rad < EPSILON(rad)) CYCLE
                END IF
                  
                IF ( BoundaryReorder(k)==0 ) THEN
                  NumberOfBoundaryNodes = NumberOfBoundaryNodes + 1
                  BoundaryReorder(k) = NumberOfBoundaryNodes
                END IF
              END DO
            END IF
          END IF
        END IF
      END DO
      DEALLOCATE( Condition )
    END DO
        
    IF (ParEnv % PEs>1 )  THEN
!------------------------------------------------------------------------------
!   If parallel execution, check for parallel matrix initializations
!------------------------------------------------------------------------------
      ALLOCATE( n_count(ParEnv% PEs),n_index(ParEnv % PEs) )
      n_count = 0
      IF ( NumberOfBoundaryNodes>0 ) THEN
        DO i=1,Mesh % NumberOfNodes
          IF (BoundaryReorder(i)<=0 ) CYCLE
          IF (.NOT.Mesh % ParallelInfo % GInterface(i) ) CYCLE

          nlist => Mesh % ParallelInfo % NeighbourList(i) % Neighbours
          DO j=1,SIZE(nlist)
            k = nlist(j)+1
            IF ( k-1 == ParEnv % myPE ) CYCLE
            n_count(k) = n_count(k)+1
          END DO
        END DO
        DO i=1,ParEnv % PEs
          IF ( n_count(i)>0 ) &
            ALLOCATE( n_index(i) % buff(n_count(i)) )
        END DO
        n_count = 0
        DO i=1,Mesh % NumberOfNodes
          IF (BoundaryReorder(i)<=0 ) CYCLE
          IF (.NOT.Mesh % ParallelInfo % GInterface(i) ) CYCLE

          nlist => Mesh % ParallelInfo % NeighbourList(i) % Neighbours
          DO j=1,SIZE(nlist)
            k = nlist(j)+1
            IF ( k == ParEnv % myPE+1 ) CYCLE
            n_count(k) = n_count(k)+1
            n_index(k) % buff(n_count(k)) = Mesh % Parallelinfo % &
                 GlobalDOFs(i)
          END DO
        END DO
      END IF

      DO i=1,ParEnv % PEs
        IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
           CALL MPI_BSEND( n_count(i), 1, MPI_INTEGER, i-1, &
                800, ELMER_COMM_WORLD, ierr )
           IF ( n_count(i)>0 ) &
             CALL MPI_BSEND( n_index(i) % buff, n_count(i), MPI_INTEGER, i-1, &
                 801, ELMER_COMM_WORLD, ierr )
        END IF
      END DO

      DO i=1,ParEnv % PEs
        IF ( n_count(i)>0 ) DEALLOCATE( n_index(i) % Buff)

        IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
           CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                800, ELMER_COMM_WORLD, status, ierr )
           IF ( n>0 ) THEN
             ALLOCATE( gbuff(n) )
             proc = status(MPI_SOURCE)
             CALL MPI_RECV( gbuff, n, MPI_INTEGER, proc, &
                 801, ELMER_COMM_WORLD, status, ierr )

             DO j=1,n
               k = SearchNodeL( Mesh % ParallelInfo, gbuff(j), Mesh % NumberOfNodes )
               IF ( k>0 ) THEN
                 IF ( BoundaryReorder(k)<= 0 ) THEN
                   NumberOfBoundaryNodes = NumberOfBoundaryNodes + 1
                   BoundaryReorder(k) = NumberOfBoundaryNodes
                 END IF
               END IF
             END DO
             DEALLOCATE(gbuff)
           END IF
        END IF
      END DO
      DEALLOCATE( n_index, n_count )
    END IF


    ! We add the normals for p-elements at 2nd stage to get higher indexes for them.
    ! This way the lower indexes should be the same as for non p-elements.
    ! Also these dofs do not to be communicated.
    !-------------------------------------------------------------------------------
    IF( pDisp ) THEN
      IF( ListCheckPresentAnyBC( Model,TRIM(VariableName) // ' Condition') ) THEN
        CALL Fatal(Caller,'Cannot deal with conditional n-t condition and p-elements')
      END IF      
      
      DO t=Mesh % NumberOfBulkElements + 1, Mesh % NumberOfBulkElements + &
          Mesh % NumberOfBoundaryElements
        
        Element => Model % Elements(t)
        IF ( Element % TYPE % ElementCode < 200 )  CYCLE
        
        n = Element % TYPE % NumberOfNodes
        np = mGetElementDOFs(pIndexes,Element,USolver=DispVar % Solver)
        Indexes => pIndexes

        DO i=1,Model % NumberOfBCs
          IF ( Element % BoundaryInfo % Constraint == &
              Model % BCs(i) % Tag ) THEN
            IF ( ListGetLogical( Model % BCs(i) % Values,VariableName, gotIt) ) THEN
              Found = ListGetLogical( Model % BCs(i) % Values, &
                  TRIM(VariableName) // ' Rotate',gotIt)
              IF ( Found .OR. .NOT. GotIt ) THEN
                DO j=n+1,np
                  k = Indexes(j)
                  IF ( BoundaryReorder(k)==0 ) THEN
                    NumberOfBoundaryNodes = NumberOfBoundaryNodes + 1
                    BoundaryReorder(k) = NumberOfBoundaryNodes
                  END IF
                END DO
              END IF
            END IF
          END IF
        END DO
      END DO
    END IF

    CALL Info(Caller,'Number of normal-tangential dofs: '&
        //I2S(NumberOfBoundaryNodes),Level=10)
    
!------------------------------------------------------------------------------

    n = 0
    IF( ASSOCIATED( BoundaryNormals ) ) THEN
      n = SIZE( BoundaryNormals, 1)
    END IF

    IF( n > 0 .AND. NumberOfBoundaryNodes /= n ) THEN
      DEALLOCATE( BoundaryNormals, BoundaryTangent1, BoundaryTangent2 )
    END IF

    IF ( NumberOfBoundaryNodes == 0 ) THEN        
      DEALLOCATE( BoundaryReorder ) 
    ELSE
      IF( n /= NumberOfBoundaryNodes ) THEN
        ALLOCATE( BoundaryNormals(NumberOfBoundaryNodes,3), & 
            BoundaryTangent1(NumberOfBoundaryNodes,3), &
            BoundaryTangent2(NumberOfBoundaryNodes,3) )
      END IF

      BoundaryNormals  = 0.0_dp
      BoundaryTangent1 = 0.0_dp
      BoundaryTangent2 = 0.0_dp
    END IF
    
!------------------------------------------------------------------------------
  END SUBROUTINE CheckNormalTangentialBoundary
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Average boundary normals for nodes. The average boundary normals
!> may be beneficial as they provide more continuous definition of normal
!> over curved boundaries. 
!------------------------------------------------------------------------------
   SUBROUTINE AverageBoundaryNormals( Model, VariableName,    &
       NumberOfBoundaryNodes, BoundaryReorder, BoundaryNormals, &
       BoundaryTangent1, BoundaryTangent2, dim )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    INTEGER, POINTER :: BoundaryReorder(:)
    INTEGER :: NumberOfBoundaryNodes, dim
    REAL(KIND=dp), POINTER :: BoundaryNormals(:,:),BoundaryTangent1(:,:), &
        BoundaryTangent2(:,:)
    CHARACTER(LEN=*) :: VariableName
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: ElementNodes
    INTEGER :: i,j,k,l,m,n,np,t,iBC,ierr,proc,i1,i2,k1,k2
    LOGICAL :: GotIt, Found, PeriodicNormals, Conditional
    REAL(KIND=dp) :: s,Bu,Bv,Nrm(3),Basis(32),DetJ
    INTEGER, POINTER :: Indexes(:)
    TYPE(Matrix_t), POINTER :: Projector
    REAL(KIND=dp), ALLOCATABLE :: Condition(:)
    TYPE(Variable_t), POINTER :: NrmVar, Tan1Var, Tan2Var
    LOGICAL, ALLOCATABLE :: Done(:), NtMasterBC(:), NtSlaveBC(:)
    REAL(KIND=dp), POINTER :: SetNormal(:,:), Rot(:,:)
    REAL(KIND=dp), TARGET :: x(Model % MaxElementNodes)
    REAL(KIND=dp), TARGET :: y(Model % MaxElementNodes)
    REAL(KIND=dp), TARGET :: z(Model % MaxElementNodes)

    TYPE buff_t
      INTEGER, ALLOCATABLE :: buff(:)
      REAL(KIND=dp), ALLOCATABLE :: normals(:)
    END TYPE buff_t
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    INTEGER, POINTER :: nlist(:)
    TYPE(Buff_t), ALLOCATABLE :: n_index(:)
    REAL(KIND=dp), ALLOCATABLE :: nbuff(:)
    INTEGER, ALLOCATABLE :: n_count(:), gbuff(:), n_comp(:)
    LOGICAL :: MassConsistent, LhsSystem, RotSystem, RotationalNormals
    LOGICAL, ALLOCATABLE :: LhsTangent(:),RhsTangent(:)
    INTEGER :: LhsConflicts, NormalConflicts, ConflictCount
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: Origin(3),Axis(3)
    INTEGER, TARGET :: pIndexes(12)
    REAL(KIND=dp), POINTER :: Pwrk(:,:)
    LOGICAL :: OneSidedNormals,pDisp
    LOGICAL :: NtBoss, AnyNtBoss, ThisBoss, NeedToAverage
    INTEGER :: NtBossCount
    LOGICAL, ALLOCATABLE :: NtBossTag(:)
    TYPE(Variable_t), POINTER :: DispVar
    CHARACTER(*), PARAMETER :: Caller = 'AverageBoundaryNormals'
    !------------------------------------------------------------------------------

    CALL Info(Caller,'Setting boundary normals for n-t conditions',Level=8)
    
    pDisp = .FALSE.
    NULLIFY( DispVar )
    DO i=1,Model % NumberOfSolvers
      IF( ListGetLogical( Model % Solvers(i) % Values,'Use p Normals',GotIt ) ) THEN
        DispVar => Model % Solvers(i) % Variable
        pDisp = .TRUE.
        EXIT
      END IF
    END DO

    Mesh => Model % Mesh
    ElementNodes % x => x
    ElementNodes % y => y
    ElementNodes % z => z

    NeedToAverage = .FALSE.
    
    
    ! Tag all nodes that have priority over conflicting normal-tangential BCs.
    AnyNtBoss = ListGetLogicalAnyBC( Model,'Normal-Tangential Priority')
    IF(AnyNtBoss) THEN
      ALLOCATE(NtBossTag(Mesh % NumberOfNodes) ) 
      NtBossTag = .FALSE.
      DO t=Mesh % NumberOfBulkElements + 1, Mesh % NumberOfBulkElements + &
          Mesh % NumberOfBoundaryElements
        Element => Model % Elements(t)
        DO i=1,Model % NumberOfBCs
          IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN
            ThisBoss = ListGetLogical( Model % BCs(i) % Values,'Normal-Tangential Priority', gotIt)
            IF( ThisBoss ) NtBossTag(Element % NodeIndexes) = .TRUE.
            EXIT
          END IF
        END DO
      END DO
      n = COUNT( NtBossTag )
      CALL Info(Caller,'Number of nodes on with normal-tangential priority is: '//I2S(n),Level=10)
      NtBossCount = 0
    END IF

    
    Mesh => Model % Mesh
    NrmVar => VariableGet( Mesh % Variables, 'Normals' )
    
    IF ( ASSOCIATED(NrmVar) ) THEN
      
      IF ( NumberOfBoundaryNodes >0 ) THEN
        BoundaryNormals = 0._dp
        DO i=1,Model % NumberOfNodes
           k = BoundaryReorder(i)
           IF (k>0 ) THEN
             DO l=1,NrmVar % DOFs
                BoundaryNormals(k,l) = NrmVar % Values( NrmVar % DOFs* &
                             (NrmVar % Perm(i)-1)+l)
             END DO
           END IF
         END DO
      END IF

    ELSE

!------------------------------------------------------------------------------
!   Compute sum of elementwise normals for nodes on boundaries
!------------------------------------------------------------------------------
      
      IF ( NumberOfBoundaryNodes>0 ) THEN
        ALLOCATE( n_comp(SIZE(BoundaryReorder)) )
        n_comp = 0

        BoundaryNormals = 0._dp
        ConflictCount = 0

        DO t=Model % NumberOfBulkElements + 1, Model % NumberOfBulkElements + &
                      Model % NumberOfBoundaryElements
          Element => Model % Elements(t)
          IF ( Element % TYPE  % ElementCode < 200 ) CYCLE

          n = Element % TYPE % NumberOfNodes
          IF(pDisp) THEN
            np = mGetElementDOFs(pIndexes,Element,USolver=DispVar % Solver)
            Indexes => pIndexes
          ELSE
            np = n
            Indexes => Element % NodeIndexes
          END IF

          ElementNodes % x(1:n) = Model % Nodes % x(Indexes(1:n))
          ElementNodes % y(1:n) = Model % Nodes % y(Indexes(1:n))
          ElementNodes % z(1:n) = Model % Nodes % z(Indexes(1:n))

          ALLOCATE(Condition(n))

          DO i=1,Model % NumberOfBCs
            IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN
              BC => Model % BCs(i) % Values

              IF ( ListGetLogical( BC, VariableName, gotIt) ) THEN
                ThisBoss = .FALSE.
                IF( AnyNtBoss ) THEN
                  ThisBoss = ListGetLogical(BC,'Normal-Tangential Priority', gotIt)
                END IF
                                
                RotationalNormals = ListGetLogical(BC,'Rotational Normals',gotIt)

                Found = ListGetLogical( BC, TRIM(VariableName) // ' Rotate',gotIt)
                IF ( Found .OR. .NOT. Gotit ) THEN
                  MassConsistent = ListGetLogical( BC,'Mass Consistent Normals',gotIt)

                  IF( RotationalNormals ) THEN
                    Pwrk => ListGetConstRealArray(BC,'Normals Origin',GotIt )
                    IF( GotIt ) THEN
                      IF( SIZE(Pwrk,1) /= 3 .OR. SIZE(Pwrk,2) /= 1 ) THEN
                        CALL Fatal(Caller,'Size of > Normals Origin < should be 3!')
                      END IF
                      Origin = Pwrk(1:3,1)
                    ELSE
                      ! Default origin is the origin.
                      Origin = 0.0_dp
                    END IF
                    Pwrk => ListGetConstRealArray(BC,'Normals Axis',GotIt )
                    IF( GotIt ) THEN
                      IF( SIZE(Pwrk,1) /= 3 .OR. SIZE(Pwrk,2) /= 1 ) THEN
                        CALL Fatal(Caller,'Size of > Normals Axis < should be 3!')
                      END IF
                      Axis = Pwrk(1:3,1)
                      ! Normalize axis is it should just be used for the direction
                      Axis = Axis / SQRT( SUM( Axis*Axis ) )
                    ELSE
                      ! Default axis is z-axis.
                      Axis = 0.0_dp
                      Axis(3) = 1.0_dp
                    END IF
                  END IF
                  
                  Condition(1:n) = ListGetReal( BC,&
                       TRIM(VariableName) // ' Condition', n, Indexes, Conditional )

                  DO j=1,n
                    IF ( Conditional ) THEN
                      IF( Condition(j) < 0._dp ) CYCLE
                    END IF

                    IF( AnyNtBoss ) THEN
                      IF(NtBossTag(Indexes(j)) .AND. .NOT. ThisBoss) THEN
                        NtBossCount = NtBossCount + 1
                        CYCLE
                      END IF
                    END IF                    
                    
                    k = BoundaryReorder( Indexes(j) )
                    IF (k>0) THEN
                      nrm = 0._dp
                      IF( RotationalNormals ) THEN
                        nrm(1) = ElementNodes % x(j)
                        nrm(2) = ElementNodes % y(j)
                        nrm(3) = ElementNodes % z(j)

                        !PRINT *,'nrm:',j,nrm
                        
                        nrm = nrm - Origin
                        nrm = nrm - SUM( nrm * Axis ) * Axis
                        nrm = nrm / SQRT( SUM( nrm * nrm ) )
                      ELSE 
                        IF (MassConsistent) THEN
                          IF(j>n) CYCLE
                          CALL IntegMassConsistent(j,n,nrm)
                        ELSE
                          Bu = Element % TYPE % NodeU(j)
                          Bv = Element % TYPE % NodeV(j)
                          nrm = NormalVector(Element,ElementNodes,Bu,Bv,.TRUE.)
                        END IF
                        NeedToAverage = .TRUE.
                        
                        l = n_comp(Indexes(j))
                        n_comp(Indexes(j)) = l + 1
                        IF( l > 0 ) THEN
                          IF( SUM( BoundaryNormals(k,:) * nrm ) < -EPSILON(s) ) THEN
                            ConflictCount = ConflictCount + 1
                            !CALL Warn(Caller,'Node '//I2S(Indexes(j))//' has conflicting normal directions!')
                          END IF
                        END IF
                      END IF
                      
                      BoundaryNormals(k,:) = BoundaryNormals(k,:) + nrm
                    END IF
                  END DO
                END IF
              END IF
            END IF
          END DO
          DEALLOCATE(Condition)
        END DO

        IF(AnyNtBoss ) THEN
          CALL Info(Caller,'Number of priority nodes for normal-tangential dofs: '&
              //I2S(NtBossCount),Level=10)
        END IF
                
        IF( ConflictCount > 0 ) THEN
          CALL Info(Caller,'There are '//I2S(ConflictCount)//' conflicting normal directions!',Level=8)
        END IF
        
        ! Here we go through the periodic projectors and average the normals
        ! such that the normals are the same where the nodes are the same.
        !--------------------------------------------------------------------
        DO iBC=1,Model % NumberOfBCs
          Projector => Model % BCs(iBC) % PMatrix
          IF ( .NOT. ASSOCIATED( Projector ) ) CYCLE

          ! This is the legacy periodic projector.
          ! The mortars etc. should be treated differently. 
          IF( Projector % ProjectorType /= PROJECTOR_TYPE_NODAL ) CYCLE
          BC => Model % BCs(iBC) % Values

          ! This is already exact.
          IF( ListGetLogical(BC,'Rotational Normals',gotIt) ) CYCLE
          
          ! TODO: consistent normals, if rotations given:
          ! ---------------------------------------------
          Rot => ListGetConstRealArray(BC,'Periodic BC Rotate', Found )
          IF ( Found .AND. ASSOCIATED(Rot) ) THEN
            IF ( ANY(Rot/=0) ) THEN
              ALLOCATE( Done(SIZE(BoundaryNormals,1)) )
              Done=.FALSE.
              DO i=1,Projector % NumberOfRows
                 k = BoundaryReorder(Projector % InvPerm(i))
                 IF ( k <= 0 ) CYCLE
                 DO l=Projector % Rows(i),Projector % Rows(i+1)-1
                   IF ( Projector % Cols(l) <= 0 ) CYCLE
                   m = BoundaryReorder(Projector % Cols(l))
                   IF ( m>0 ) THEN
                     IF ( .NOT.Done(m) ) THEN
                       Done(m) = .TRUE.
                       BoundaryNormals(m,:) = -BoundaryNormals(m,:)
                     END IF
                   END IF
                 END DO
              END DO
              DEALLOCATE(Done)
              CYCLE
            END IF
          END IF
          
          ! Here we are projecting with transpose of the projector which is not
          ! really exact generally, but is usually better than not considering the values
          ! at all!
          !-------------------------------------------------------------------------------
          OneSidedNormals = ListGetLogical(BC,'One Sided Normals',Found ) 
          IF(.NOT. OneSidedNormals ) THEN
            NeedToAverage = .TRUE.
            DO i=1,Projector % NumberOfRows
              k = BoundaryReorder(Projector % InvPerm(i))
              IF ( k <= 0 ) CYCLE
              DO l=Projector % Rows(i),Projector % Rows(i+1)-1
                IF ( Projector % Cols(l) <= 0 ) CYCLE
                m = BoundaryReorder(Projector % Cols(l))
                IF ( m>0 ) THEN
                  BoundaryNormals(k,:) = BoundaryNormals(k,:) + &
                      Projector % Values(l) * BoundaryNormals(m,:)
                END IF
              END DO
            END DO
          END IF
            
          ! Ok, now we need to nullify the values so that we can apply the projector
          ! in the next sequence. This used to be done before without the upper part
          !--------------------------------------------------------------------------
          DO i=1,Projector % NumberOfRows
            k = BoundaryReorder(Projector % InvPerm(i))
            IF ( k <= 0 ) CYCLE
            DO l=Projector % Rows(i),Projector % Rows(i+1)-1
              IF ( Projector % Cols(l) <= 0 ) CYCLE
              m = BoundaryReorder(Projector % Cols(l))
              IF ( m>0 ) BoundaryNormals(m,:) = 0.0_dp
            END DO
          END DO
        END DO

        ! Here we use the values of the master side to have same values
        ! on the slave side as well. So this is hierarchical.
        ! For rotational normals everything is already ok, so don't do those!
        !----------------------------------------------------------------
        DO iBC=1,Model % NumberOfBCs
          IF( ListGetLogical(Model % BCs(iBc) % Values,'Rotational Normals',gotIt) ) CYCLE
          
          Projector => Model % BCs(iBC) % PMatrix
           IF ( .NOT. ASSOCIATED( Projector ) ) CYCLE
           IF( Projector % ProjectorType /= PROJECTOR_TYPE_NODAL ) CYCLE
          
           ! TODO: consistent normals, if rotations given:
           ! ---------------------------------------------
           BC => Model % BCs(iBC) % Values
           Rot => ListGetConstRealArray(BC,'Periodic BC Rotate', Found )
           IF ( Found .AND. ASSOCIATED(Rot) ) THEN
             IF ( ANY(Rot/=0) ) CYCLE
           END IF

           NeedToAverage = .TRUE.
           
           DO i=1,Projector % NumberOfRows
              k = BoundaryReorder(Projector % InvPerm(i))
              IF ( k <= 0 ) CYCLE
              DO l=Projector % Rows(i),Projector % Rows(i+1)-1
                IF ( Projector % Cols(l) <= 0 ) CYCLE
                m = BoundaryReorder(Projector % Cols(l))
                IF ( m > 0 ) THEN
                  BoundaryNormals(m,:) = BoundaryNormals(m,:) + &
                      Projector % Values(l) * BoundaryNormals(k,:)
                END IF
              END DO
            END DO
        END DO
      END IF

      ! Communicate normals in parallel case so that they are consistent
      ! over the interfaces
      !-----------------------------------------------------------------
      i = 0
      IF( ParEnv % PEs > 1 ) THEN
        IF(NeedToAverage) i = 1
        i = ParallelReduction(i)
      END IF

      IF (i > 0 ) THEN
        ALLOCATE( n_count(ParEnv% PEs),n_index(ParEnv % PEs) )
        n_count = 0

        IF ( NumberOfBoundaryNodes>0 ) THEN
          DO i=1,Mesh % NumberOfNodes
            IF (BoundaryReorder(i)<=0 .OR. n_comp(i)<=0 ) CYCLE
            IF (.NOT.Mesh % ParallelInfo % GInterface(i) ) CYCLE
  
            nlist => Mesh % ParallelInfo % NeighbourList(i) % Neighbours
            DO j=1,SIZE(nlist)
              k = nlist(j)+1
              IF ( k-1 == ParEnv % myPE ) CYCLE
              n_count(k) = n_count(k)+1
            END DO
          END DO
          DO i=1,ParEnv % PEs
            IF ( n_count(i)>0 ) &
                ALLOCATE( n_index(i) % buff(n_count(i)), &
                        n_index(i) % normals(3*n_count(i)) )
          END DO

          n_count = 0
          DO i=1,Model % NumberOfNodes
            IF (BoundaryReorder(i)<=0 .OR. n_comp(i)<=0 ) CYCLE
            IF (.NOT.Mesh % ParallelInfo % GInterface(i) ) CYCLE

            nlist => Mesh % ParallelInfo % NeighbourList(i) % Neighbours
            DO j=1,SIZE(nlist)
              k = nlist(j)+1
              IF ( k-1 == ParEnv % myPE ) CYCLE
              n_count(k) = n_count(k)+1
              n_index(k) % buff(n_count(k)) = Mesh % Parallelinfo % &
                 GlobalDOFs(i)
              l = BoundaryReorder(i)
              n_index(k) % normals(3*n_count(k)-2)=BoundaryNormals(l,1)
              n_index(k) % normals(3*n_count(k)-1)=BoundaryNormals(l,2)
              n_index(k) % normals(3*n_count(k)-0)=BoundaryNormals(l,3)
            END DO
          END DO
        END IF

        DO i=1,ParEnv % PEs
          IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
            CALL MPI_BSEND( n_count(i), 1, MPI_INTEGER, i-1, &
                900, ELMER_COMM_WORLD, ierr )
            IF ( n_count(i)>0 ) THEN
              CALL MPI_BSEND( n_index(i) % buff, n_count(i), MPI_INTEGER, i-1, &
                  901, ELMER_COMM_WORLD, ierr )
              CALL MPI_BSEND( n_index(i) % normals, 3*n_count(i), MPI_DOUBLE_PRECISION, &
                    i-1,  902, ELMER_COMM_WORLD, ierr )
            END IF
          END IF
        END DO
        DO i=1,ParEnv % PEs
          IF ( n_count(i)>0 ) DEALLOCATE( n_index(i) % Buff, n_index(i) % Normals)

          IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
             CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                    900, ELMER_COMM_WORLD, status, ierr )
             IF ( n>0 ) THEN
               proc = status(MPI_SOURCE)
               ALLOCATE( gbuff(n), nbuff(3*n) )
               CALL MPI_RECV( gbuff, n, MPI_INTEGER, proc, &
                   901, ELMER_COMM_WORLD, status, ierr )

               CALL MPI_RECV( nbuff, 3*n, MPI_DOUBLE_PRECISION, proc, &
                    902, ELMER_COMM_WORLD, status, ierr )

               DO j=1,n
                 k = SearchNodeL( Mesh % ParallelInfo, gbuff(j), Mesh % NumberOfNodes )
                 IF ( k>0 ) THEN
                   n_comp(k) = n_comp(k)+1
                   l = BoundaryReorder(k)
                   IF ( l>0 ) THEN
                     BoundaryNormals(l,1)=BoundaryNormals(l,1)+nbuff(3*j-2)
                     BoundaryNormals(l,2)=BoundaryNormals(l,2)+nbuff(3*j-1)
                     BoundaryNormals(l,3)=BoundaryNormals(l,3)+nbuff(3*j-0)
                   END IF
                 END IF
               END DO
               DEALLOCATE(gbuff, nbuff)
             END IF
          END IF
        END DO
        DEALLOCATE( n_index, n_count )
      END IF
    END IF

!------------------------------------------------------------------------------
!   normalize 
!------------------------------------------------------------------------------
    IF ( NumberOfBoundaryNodes>0 ) THEN

      RotSystem = ListGetLogical(Model % Simulation,'Use Cylinder System',Found) 
      LhsSystem = ListGetLogical(Model % Simulation,'Use Lhs System',Found) 
      IF(.NOT. Found ) LhsSystem = ( dim == 3 .AND. .NOT. RotSystem )

      IF( LhsSystem ) THEN
        ALLOCATE( NtMasterBC( Model % NumberOfBCs ), NtSlaveBC( Model % NumberOfBCs ) )
        NtMasterBC = .FALSE.; NtSlaveBC = .FALSE.

        DO i = 1, Model % NumberOfBcs
          IF( .NOT. ListCheckPrefix( Model % BCs(i) % Values,'Normal-Tangential') ) CYCLE
          
          j = ListGetInteger( Model % BCs(i) % Values,'Mortar BC',Found )
          IF( .NOT. Found ) THEN
            j = ListGetInteger( Model % BCs(i) % Values,'Contact BC',Found )
          END IF
          IF( j == 0 .OR. j > Model % NumberOfBCs ) CYCLE

          NtSlaveBC( i ) = .TRUE.
          NtMasterBC( j ) = .TRUE.
        END DO
        LhsSystem = ANY( NtMasterBC )
      END IF

      IF( LhsSystem ) THEN
        DO i = 1, Model % NumberOfBcs
          IF( NtSlaveBC( i ) .AND. NtMasterBC( i ) ) THEN
            CALL Warn(Caller,'BC '//I2S(i)//' is both N-T master and slave!')
          END IF
        END DO

        ALLOCATE( LhsTangent( Model % NumberOfNodes ) )
        LhsTangent = .FALSE.

        ALLOCATE( RhsTangent( Model % NumberOfNodes ) )
        RhsTangent = .FALSE. 

        DO t=Model % NumberOfBulkElements + 1, Model % NumberOfBulkElements + &
            Model % NumberOfBoundaryElements
          Element => Model % Elements(t)
          IF ( Element % TYPE  % ElementCode < 200 ) CYCLE
          
          n = Element % TYPE % NumberOfNodes
          Indexes => Element % NodeIndexes
          
          DO i=1,Model % NumberOfBCs
            IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN
              IF( NtMasterBC(i) ) LhsTangent( Indexes ) = .TRUE.
              IF( NtSlaveBC(i) ) RhsTangent( Indexes ) = .TRUE.
              EXIT
            END IF
          END DO
        END DO

        LhsConflicts = COUNT( LhsTangent .AND. RhsTangent )
        IF( LhsConflicts > 0 ) THEN
          CALL Warn(Caller,&
              'There are '//I2S(LhsConflicts)//' nodes that could be both rhs and lhs!')
        END IF
      END IF
            
      
      ! Normalize the normals and compute the tangent directions.
      !----------------------------------------------------------
      DO i=1,Model % NumberOfNodes
        k = BoundaryReorder(i) 
        IF ( k > 0 ) THEN
          s = SQRT( SUM( BoundaryNormals(k,:)**2 ) )
          IF ( s > TINY(s) ) THEN
            BoundaryNormals(k,:) = BoundaryNormals(k,:) / s
            IF ( dim > 2 ) THEN
              CALL TangentDirections( BoundaryNormals(k,:),  &
                  BoundaryTangent1(k,:), BoundaryTangent2(k,:) )
              IF( RotSystem ) THEN
                ! We want to always have n_z as positive in cylinder system.
                IF( BoundaryTangent2(k,3) < 0.0 ) THEN
                  BoundaryTangent2(k,:) = -BoundaryTangent2(k,:)
                END IF
              ELSE IF( LhsSystem ) THEN
                IF( LhsTangent(i) ) THEN
                  BoundaryTangent2(k,:) = -BoundaryTangent2(k,:)
                END IF
              END IF
            END IF
          ELSE
            CALL Warn(Caller,'Suspiciously small normal for node: '//I2S(i))
          END IF
        END IF        
      END DO
      

      ! Inherit the normal direction for 2nd order p-elements from the nodes.
      !----------------------------------------------------------------------
      IF( pDisp ) THEN
        DO t=Mesh % NumberOfBulkElements + 1, Mesh % NumberOfBulkElements + &
            Mesh % NumberOfBoundaryElements
          
          Element => Model % Elements(t)
          IF ( Element % TYPE % ElementCode < 200 )  CYCLE
          
          n = Element % TYPE % NumberOfNodes
          np = mGetElementDOFs(pIndexes,Element,USolver=DispVar % Solver)
          Indexes => pIndexes

          IF(ALL(BoundaryReorder(Indexes(1:n)) == 0 ) ) CYCLE
          
          DO i=n+1,np
            i1 = i-n
            i2 = i+1-n
            IF(i2>n) i2=1

            k1 = BoundaryReorder(Indexes(i1)) 
            IF(k1==0) CYCLE

            k2 = BoundaryReorder(Indexes(i2))
            IF(k2==0) CYCLE
            
            k = BoundaryReorder(Indexes(i)) 

            BoundaryNormals(k,:) = ( BoundaryNormals(k1,:) + BoundaryNormals(k2,:) ) / 2

            ! Even though the two normals have unit length their mean may not have unit length. 
            s = SQRT( SUM( BoundaryNormals(k,:)**2 ) )
            IF ( s > TINY(s) ) THEN
              BoundaryNormals(k,:) = BoundaryNormals(k,:) / s
            ELSE
              CALL Warn(Caller,'Starnegly small normal for dofs: '//I2S(i))
            END IF
              
            IF( dim > 2 ) THEN
              BoundaryTangent1(k,:) = ( BoundaryTangent1(k1,:) + BoundaryTangent1(k2,:) ) / 2
              BoundaryTangent2(k,:) = ( BoundaryTangent2(k1,:) + BoundaryTangent2(k2,:) ) / 2

              s = SQRT( SUM( BoundaryTangent1(k,:)**2 ) )
              IF ( s > TINY(s) ) BoundaryTangent1(k,:) = BoundaryTangent1(k,:) / s

              s = SQRT( SUM( BoundaryTangent2(k,:)**2 ) )
              IF ( s > TINY(s) ) BoundaryTangent2(k,:) = BoundaryTangent2(k,:) / s
            END IF
          END DO
        END DO
      END IF

      ! Save the normals and tangents as fields if requested. 
      !----------------------------------------------------------
      IF( ListGetLogical( Model % Simulation,'Save Averaged Normals',Found ) ) THEN
        CALL Info(Caller,'Saving averaged boundary normals to variable: Averaged Normals')
        NrmVar => VariableGet( Mesh % Variables, 'Averaged Normals' )
        
        IF(.NOT. ASSOCIATED( NrmVar ) ) THEN
          CALL VariableAddVector( Mesh % Variables, Mesh, Model % Solver,'Averaged Normals',3,&
              Perm = BoundaryReorder )
          NrmVar => VariableGet( Mesh % Variables, 'Averaged Normals' )
        END IF
            
        DO i=1,Model % NumberOfNodes
          k = BoundaryReorder(i)
          IF (k>0 ) THEN
            DO l=1,NrmVar % DOFs
              NrmVar % Values( NrmVar % DOFs* &
                  (NrmVar % Perm(i)-1)+l)  = BoundaryNormals(k,l)
            END DO
          END IF
        END DO

        IF( dim > 2 .AND. ListGetLogical( Model % Simulation,'Save Averaged Tangents',Found ) ) THEN
          Tan1Var => VariableGet( Mesh % Variables, 'Averaged First Tangent' )
          Tan2Var => VariableGet( Mesh % Variables, 'Averaged Second Tangent' )

          IF(.NOT. ASSOCIATED( Tan1Var ) ) THEN
            CALL VariableAddVector( Mesh % Variables, Mesh, Model % Solver,&
                'Averaged First Tangent',3, Perm = BoundaryReorder )
            Tan1Var => VariableGet( Mesh % Variables, 'Averaged First Tangent' )
            CALL VariableAddVector( Mesh % Variables, Mesh, Model % Solver,&
                'Averaged Second Tangent',3, Perm = BoundaryReorder )
            Tan2Var => VariableGet( Mesh % Variables, 'Averaged Second Tangent' )
          END IF
          
          DO i=1,Model % NumberOfNodes
            k = BoundaryReorder(i)
            IF (k>0 ) THEN
              DO l=1,Tan1Var % DOFs
                Tan1Var % Values( Tan1Var % DOFs* &
                    (Tan1Var % Perm(i)-1)+l)  = BoundaryTangent1(k,l)
                Tan2Var % Values( Tan2Var % DOFs* &
                    (Tan2Var % Perm(i)-1)+l)  = BoundaryTangent2(k,l)
              END DO
            END IF
          END DO
        END IF
      END IF
    END IF

    IF( InfoActive(25) ) THEN
      DO i=1,3        
        CALL VectorValuesRange(BoundaryNormals(:,i),NumberOfBoundaryNodes,'Normal '//I2S(i))
      END DO
      IF( dim > 2 ) THEN
        DO i=1,3        
          CALL VectorValuesRange(BoundaryTangent1(:,i),NumberOfBoundaryNodes,'Tangent1 '//I2S(i))
          CALL VectorValuesRange(BoundaryTangent2(:,i),NumberOfBoundaryNodes,'Tangent2 '//I2S(i))
        END DO
      END IF
    END IF

 CONTAINS

    SUBROUTINE IntegMassConsistent(j,n,nrm)
      INTEGER :: t,j,n
      LOGICAL :: stat
      REAL(KIND=dp) :: detJ,Basis(n),nrm(:),lnrm(3)

      TYPE(GaussIntegrationPoints_t) :: IP

      !----------------------
      IP = GaussPoints(Element)
      DO t=1,IP % n
        stat = ElementInfo(Element, ElementNodes, IP % U(t), &
               IP % v(t), IP % W(t), detJ, Basis)

        lnrm = NormalVector(Element,ElementNodes, &
              IP % U(t),IP % v(t),.TRUE.)

        nrm = nrm + IP % s(t) * lnrm * detJ * Basis(j)
      END DO
    END SUBROUTINE IntegMassConsistent

!------------------------------------------------------------------------------
  END SUBROUTINE AverageBoundaryNormals
!------------------------------------------------------------------------------


END MODULE BoundaryConditionUtils
