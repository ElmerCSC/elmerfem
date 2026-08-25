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
! *  Authors: Juha Ruokolainen
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
!>  Projector and constraint matrix generation utilities.
!>  Extracted from SolverUtils.
!------------------------------------------------------------------------------

MODULE ProjectorUtils

    USE SolverBasics
    IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------



  ! Create Linear constraints from mortar BCs:
  ! -------------------------------------------   
  SUBROUTINE GenerateProjectors(Model,Solver,Nonlinear,SteadyState) 
    
     TYPE(Model_t) :: Model
     TYPE(Solver_t), TARGET :: Solver
     LOGICAL, OPTIONAL :: Nonlinear, SteadyState

     LOGICAL :: IsNonlinear,IsSteadyState,Timing, RequireNonlinear, ContactBC, &
         MortarBC, IntegralBC
     LOGICAL :: ApplyMortar, ApplyContact, ApplyIntegral, StoreCyclic, Found, &
         StaticProj, IsBodyForce
     INTEGER :: i,j,k,l,n,dsize,size0,col,row,dim
     TYPE(ValueList_t), POINTER :: BC
     TYPE(Matrix_t), POINTER :: CM, CMP, CM0, CM1
     TYPE(Variable_t), POINTER :: DispVar
     REAL(KIND=dp) :: t0,rt0,rst,st,ct
     CHARACTER(*), PARAMETER :: Caller = 'GenerateProjectors'
     TYPE(Solver_t), POINTER :: PSolver
     TYPE(Matrix_t), POINTER :: Proj
     
     ApplyIntegral = ListGetLogical( Solver % Values,'Apply Integral BCs',Found)
     ApplyMortar = ListGetLogical(Solver % Values,'Apply Mortar BCs',Found) 
     ApplyContact = ListGetLogical(Solver % Values,'Apply Contact BCs',Found)     
     
     IF( .NOT. ( ApplyMortar .OR. ApplyContact .OR. ApplyIntegral ) ) RETURN
     
     ! Here we give the option to block out cyclic projector if not wanted. 
     StoreCyclic = ListGetLogical( Solver % Values,'Store Cyclic Projector', Found)
     IF(.NOT. Found ) StoreCyclic = ListGetLogical( Solver % Values,'Store Cyclic System', Found)
     PSolver => Solver
     
     i = ListGetInteger( Solver % Values,'Mortar BC Master Solver',Found ) 
     IF( Found ) THEN
       Solver % MortarBCs => CurrentModel % Solvers(i) % MortarBCs
       IF( .NOT. ASSOCIATED( Solver % MortarBCs ) ) THEN
         CALL Fatal(Caller,'Could not reuse projectors from solver: '//I2S(i))
       END IF
       CALL Info(Caller,'Reusing projectors from solver: '//I2S(i),Level=8)
       RETURN
     END IF

     CALL Info(Caller,'Generating various boundary projectors',Level=8)

     Timing = ListCheckPrefix(Solver % Values,'Projector Timing')
     IF( Timing ) THEN
       t0 = CPUTime(); rt0 = RealTime()      
     END IF

     IsNonlinear = .FALSE.
     IF( PRESENT( Nonlinear ) ) IsNonlinear = Nonlinear
     IsSteadyState = .NOT. IsNonlinear

     IF( .NOT. ASSOCIATED( Solver % MortarBCs ) ) THEN
       ALLOCATE( Solver % MortarBCs( Model % NumberOfBCs + Model % NumberOfBodyForces ) )
       DO i=1, Model % NumberOfBCs + Model % NumberOfBodyForces
         Solver % MortarBCs(i) % Projector => NULL()
       END DO
     END IF
     
     dim = CoordinateSystemDimension()

     DO i=1,Model % NumberOFBCs + Model % NumberOfBodyForces
       IF(i > Model % NumberOfBCs ) THEN
         BC => Model % BodyForces(i-Model % NumberOfBCs) % Values
         IsBodyForce = .TRUE.
       ELSE         
         BC => Model % BCs(i) % Values
         IsBodyForce = .FALSE.
       END IF
         
       k = 0
       j = ListGetInteger( BC,'Mortar BC',MortarBC)       
       j = j + ListGetInteger( BC,'Contact BC',ContactBC)       
       IntegralBC = ListGetLogical( BC,'Integral BC',Found ) 

       IF( MortarBC ) k = k+1
       IF( ContactBC ) k = k+1
       IF( IntegralBC ) k = k+1
       IF(k==0) CYCLE
         
       IF( k > 1 ) THEN
         CALL Fatal(Caller,'Boundary '//I2S(i)//' can only be one of mortar, contact and integral!')
       END IF     

       IF(IsBodyForce .AND. .NOT. IntegralBC ) THEN
         CALL Fatal(Caller,'Body Force '//I2S(i)//' can only have integral bc!')
       END IF     

       IF( InfoActive(10) ) THEN
         IF( MortarBC ) CALL Info(Caller,'Generating mortar conditions for BC: '//I2S(i))
         IF( ContactBC ) CALL Info(Caller,'Generating contact conditions for BC: '//I2S(i))
         IF( IntegralBC ) CALL Info(Caller,'Generating integral conditions for BC: '//I2S(i))
       END IF
       
       RequireNonlinear = ListGetLogical( BC,'Mortar BC Nonlinear',Found)
       IF( .NOT. Found ) THEN
         RequireNonlinear = ContactBC .AND. .NOT. ListGetLogical( BC,'Tie Contact',Found )
       END IF

       IF( IntegralBC ) RequireNonlinear = .FALSE.
       
       IF( IsNonlinear ) THEN
         IF( .NOT. RequireNonlinear ) CYCLE
       ELSE
         IF( RequireNonlinear ) CYCLE
       END IF             

       StaticProj = ListGetLogical( BC,'Mortar BC Static',Found)
       
       Proj => Solver % MortarBCs(i) % Projector
       IF( ASSOCIATED( Proj ) ) THEN
         IF( StaticProj ) CYCLE         

         IF( StoreCyclic ) THEN
           ! Don't release projectors in case they are cyclic 
           ! Instead reassign the pointer.
           CALL StoreCyclicProjector(PSolver,Proj,Found)
           IF(Found) THEN
             Solver % MortarBCs(i) % Projector => Proj
             Solver % MortarBCsChanged = .TRUE.
             CYCLE
           END IF
         ELSE  
           IF( ASSOCIATED( Proj % Ematrix ) ) THEN
             CALL FreeMatrix( Proj % Ematrix )
           END IF
           CALL FreeMatrix( Proj ) 
         END IF
       END IF

       ! Compute new projector
       IF( IntegralBC ) THEN         
         Proj => IntegralProjector(Model,Solver % Mesh, i, IsBodyForce ) 
       ELSE
         ! This is the same for mortar and contact!
         Proj => PeriodicProjector(Model,Solver % Mesh,i,j,dim,.TRUE.)
       END IF
         
       Solver % MortarBCs(i) % Projector => Proj       
       IF( ASSOCIATED( Proj ) ) THEN
         Solver % MortarBCsChanged = .TRUE.
       END IF

       ! Store new projector to the cyclic set
       IF( StoreCyclic ) THEN
         IF(.NOT. StaticProj ) CALL StoreCyclicProjector(PSolver,Proj)
       END IF
       
     END DO


     IF( Timing ) THEN
       st  = CPUTime() - t0;
       rst = RealTime() - rt0
       
       WRITE(Message,'(a,f8.2,f8.2,a)') 'Projector creation time (CPU,REAL) for '&
           //GetVarName(Solver % Variable)//': ',st,rst,' (s)'
       CALL Info(Caller,Message,Level=6)    
       
       IF( ListGetLogical(Solver % Values,'Projector Timing',Found)) THEN
         CALL ListAddConstReal(CurrentModel % Simulation,'res: projector cpu time '&
             //GetVarName(Solver % Variable),st)
         CALL ListAddConstReal(CurrentModel % Simulation,'res: projector real time '&
             //GetVarName(Solver % Variable),rst)
       END IF

       IF( ListGetLogical(Solver % Values,'Projector Timing Cumulative',Found)) THEN
         ct = ListGetConstReal(CurrentModel % Simulation,'res: cum projector cpu time '&
             //GetVarName(Solver % Variable),Found)
         st = st + ct
         ct = ListGetConstReal(CurrentModel % Simulation,'res: cum projector real time '&
             //GetVarName(Solver % Variable),Found)
         rst = rst + ct
         CALL ListAddConstReal(CurrentModel % Simulation,'res: cum projector cpu time '&
             //GetVarName(Solver % Variable),st)
         CALL ListAddConstReal(CurrentModel % Simulation,'res: cum projector real time '&
             //GetVarName(Solver % Variable),rst)
       END IF
     END IF
     
   END SUBROUTINE GenerateProjectors



  ! Create special type of projector after the linear system has been created. 
  ! -------------------------------------------------------------------------
   SUBROUTINE GenerateRobinProjectors(Model,Solver)

     TYPE(Model_t) :: Model
     TYPE(Solver_t), TARGET :: Solver

     TYPE(Solver_t), POINTER :: PSolver
     INTEGER :: i
     LOGICAL :: Found, GotSome, DoIt, IsBodyForce
     CHARACTER(*), PARAMETER :: Caller="GenerateRobinProjector"

     
     IF(.NOT. ListGetLogical( Solver % Values,'Apply Integral BCs',Found) ) RETURN

     PSolver => Solver
     GotSome = .FALSE.
     DO i=1,Model % NumberOFBCs + Model % NumberOfBodyForces
       IF(i<=Model % NumberOfBCs) THEN
         IF(.NOT. ListGetLogical( Model % BCs(i) % Values,'Flux Integral BC',Found ) ) CYCLE
         IsBodyForce = .FALSE.
       ELSE
         IF(.NOT. ListGetLogical( Model % BodyForces(i-Model % NumberOfBCs) % Values,'Flux Integral BC',Found ) ) CYCLE
         IsBodyForce = .TRUE.
       END IF
       
       CALL Info(Caller,'Generating flux integral conditions for BC: '//I2S(i))
       CALL RobinProjector(Model,PSolver, i, IsBodyForce) 
       GotSome = .TRUE.
     END DO

     IF( GotSome ) THEN
       Solver % MortarBCsChanged = .TRUE.
       ! We may want to export the lagrange multiplier as it has a physical meaning. 
       CALL ListAddNewLogical(Solver % Values,'Export Lagrange Multiplier',.TRUE.) 
     END IF
            
   END SUBROUTINE GenerateRobinProjectors



   ! This creates a projector that integrates over the BCs on the boundary such that
   ! an integral constraint for Robin type of BCs may be applied.
   !--------------------------------------------------------------------------------------
   SUBROUTINE RobinProjector(Model, Solver, BCInd, IsBodyForce )

     TYPE(Model_t) :: Model
     TYPE(Solver_t) :: Solver    
     INTEGER :: BCInd
     LOGICAL :: IsBodyForce

     TYPE(Matrix_t), POINTER :: Proj        
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(ValueList_t), POINTER :: BC
     LOGICAL :: Found
     INTEGER :: dofs
     TYPE(MortarBC_t), POINTER :: MortarBC
     CHARACTER(*), PARAMETER :: Caller="RobinProjector"

     IF(IsBodyForce) THEN
       BC => Model % BodyForces(BCind-Model % NumberOfBCs) % Values
     ELSE
       BC => Model % BCs(BCInd) % Values
     END IF
     IF( .NOT. ListGetLogical( BC,'Flux Integral BC', Found ) ) RETURN

     Proj => Solver % MortarBCs(BCInd) % Projector
     IF(.NOT. ASSOCIATED(Proj)) THEN
       Proj => AllocateMatrix()
       Proj % FORMAT = MATRIX_LIST
       Proj % ProjectorType = PROJECTOR_TYPE_ROBIN
       CALL Info(Caller,'Creating Robin flux constraint matrix for boundary: '//I2S(BCind),Level=6)
     ELSE
       CALL Info(Caller,'Updating Robin flux constraint matrix for boundary: '//I2S(BCind),Level=6)
       Proj % Values = 0.0_dp
     END IF

     Mesh => Solver % Mesh

     CALL CreateRobinProjector()

     IF(Proj % FORMAT == MATRIX_LIST ) THEN
       CALL List_toCRSMatrix(Proj)
     END IF

     IF( InfoActive(20) ) THEN
       WRITE(Message,'(A,ES12.3)') 'Sum of constraint matrix entries: ',SUM(Proj%Values)
       CALL Info(Caller,Message)
       CALL Info(Caller,'Constraint matrix cols min: '//I2S(MINVAL(Proj%Cols)))
       CALL Info(Caller,'Constraint matrix cols max: '//I2S(MAXVAL(Proj%Cols)))
       CALL Info(Caller,'Constraint matrix rows min: '//I2S(MINVAL(Proj%Rows)))
       CALL Info(Caller,'Constraint matrix rows max: '//I2S(MAXVAL(Proj%Rows)))
     END IF

     MortarBC => Solver % MortarBCs(BCind) 
     MortarBC % Projector => Proj       

     IF(.NOT. ASSOCIATED(MortarBC % Diag ) ) THEN      
       dofs = Solver % Variable % Dofs
       ALLOCATE(MortarBC % Diag(dofs))
       MortarBC % Diag = 1.0_dp
     END IF

     
   CONTAINS

     SUBROUTINE CreateRobinProjector()

       INTEGER :: i,j,k,n,t,dofs,idof,t1,t2
       REAL(KIND=dp) :: dval
       TYPE(Matrix_t), POINTER :: A
       TYPE(Variable_t), POINTER :: Var
       TYPE(Element_t), POINTER :: Element
       LOGICAL, ALLOCATABLE :: ActiveDof(:)
       INTEGER, POINTER :: Indexes(:)
       
       A => Solver % Matrix
       Var => Solver % Variable       
       n = A % NumberOfRows
       dofs = Var % dofs

       IF(.NOT. ASSOCIATED(A % BulkValues)) THEN
         CALL Fatal(Caller,'BulkValues are needed to aveluate Robin terms!')
       END IF

       IF(IsBodyForce) THEN
         t1 = 1
         t2 = Mesh % NumberOfBulkElements 
       ELSE
         t1 = Mesh % NumberOfBulkElements + 1
         t2 = (t1-1) + Mesh % NumberOfBoundaryElements
       END IF
       
       ! Mark the dofs of the matrix that are on the boundary.
       ! ActiveDof table will directly refer to the indexes of the matrix.
       ALLOCATE(ActiveDof(n))
       ActiveDof = .FALSE.       
       DO t = t1, t2
         Element => Mesh % Elements(t)
         IF(IsBodyForce) THEN
           i = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Body Force',Found)
           IF(i /= BCind - Model % NumberOfBCs) CYCLE
         ELSE
           IF ( Element % BoundaryInfo % Constraint /= Model % BCs(BCInd) % Tag ) CYCLE
         END IF
         Indexes => Element % NodeIndexes      

         IF(ANY(Var % Perm(Indexes) == 0 ) ) CYCLE
         DO i=1,dofs
           ActiveDof(dofs*(Var % Perm(Indexes)-1)+i) = .TRUE.
         END DO
       END DO
              
       DO i=1, n / dofs 
         DO idof = 1, dofs
           j = dofs*(i-1)+idof
           IF(.NOT. ActiveDof(j)) CYCLE
           
           DO k=A % Rows(j),A % Rows(j+1)-1
             IF(.NOT. ActiveDof(A % Cols(k))) CYCLE
             dval = A % Values(k) - A % BulkValues(k)
             ! We use the negative sign so that we get desired sign for the lagrange multiplier!
             CALL AddToMatrixElement(Proj, idof, A % Cols(k), -dval )
           END DO
         END DO
       END DO       
       
     END SUBROUTINE CreateRobinProjector

   END SUBROUTINE RobinProjector
   

   ! Generate constraint matrix from mortar projectors. 
   ! This routine takes each boundary projector and applies it 
   ! to the current field variable (scalar or vector) merging 
   ! all into one single projector. 
   !---------------------------------------------------------
   SUBROUTINE GenerateConstraintMatrix( Model, Solver )

     TYPE(Model_t) :: Model
     TYPE(Solver_t) :: Solver

     INTEGER, POINTER :: Perm(:)
     INTEGER :: i,j,j2,k,k2,l,l2,dofs,maxperm,permsize,bc_ind,constraint_ind,row,col,col2,&
         mcount,bcount,kk,cdofs,dim
     TYPE(Matrix_t), POINTER :: Atmp,Btmp, Ctmp
     LOGICAL :: AllocationsDone, CreateSelf, ComplexMatrix, TransposePresent, Found, &
         SetDof, SomeSet, SomeSkip, SumProjectors, NewRow, SumThis
     INTEGER, ALLOCATABLE :: SumPerm(:),SumCount(:)
     LOGICAL, ALLOCATABLE :: ActiveComponents(:), SetDefined(:)
     TYPE(ValueList_t), POINTER :: BC
     TYPE(MortarBC_t), POINTER :: MortarBC
     REAL(KIND=dp) :: wsum, Scale
     INTEGER :: rowoffset, arows, sumrow, EliminatedRows, NeglectedRows, sumrow0, k20
     LOGICAL :: ThisIsMortar, ThisIsRobin, Reorder
     LOGICAL :: AnyPriority
     INTEGER :: Priority, PrevPriority
     INTEGER, ALLOCATABLE :: BCOrdering(:), BCPriority(:)
     LOGICAL :: NeedToGenerate, ComplexSumRow 

     LOGICAL :: HaveMortarDiag, LumpedDiag, PerFlipActive, SkipConstrained
     LOGICAL, POINTER :: ConstrainedDof(:)
             
     REAL(KIND=dp) :: MortarDiag, val, valsum, EpsVal
     LOGICAL, POINTER :: PerFlip(:)
     CHARACTER(*), PARAMETER :: Caller = 'GenerateConstraintMatrix'

     LOGICAL :: IntegralBC
     REAL(KIND=dp) :: SetVal(6)
     REAL(KIND=dp), ALLOCATABLE :: PrevValues(:)
     INTEGER, ALLOCATABLE :: PrevInvPerm(:)
     TYPE(Variable_t), POINTER :: Var
     CHARACTER(:), ALLOCATABLE :: Str,MultName
     REAL(KIND=dp), ALLOCATABLE :: rsum(:)
     LOGICAL :: IsDg, IsBodyForce
     INTEGER, ALLOCATABLE :: DgSome(:)
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Element_t), POINTER :: Element
     
     ! Should we genarete the matrix
     NeedToGenerate = Solver % MortarBCsChanged
     
     Mesh => Solver % Mesh 
     IsDg = Solver % DG     

     PerFlipActive = Solver % PeriodicFlipActive
     IF( PerFlipActive ) THEN
       CALL Info(Caller,'Periodic flip is active',Level=8)
       PerFlip => Mesh % PeriodicFlip           
     END IF
     
     ! Set pointers to save the initial constraint matrix
     ! We assume that the given ConstraintMatrix is constant but we have consider it the 1st time
     IF(.NOT. Solver % ConstraintMatrixVisited ) THEN       
       IF( ASSOCIATED( Solver % Matrix % ConstraintMatrix ) ) THEN
         CALL Info(Caller,'Saving initial constraint matrix to Solver',Level=12)
         Solver % ConstraintMatrix => Solver % Matrix % ConstraintMatrix
         Solver % Matrix % ConstraintMatrix => NULL()
         NeedToGenerate = .TRUE. 
       END IF
       Solver % ConstraintMatrixVisited = .TRUE.
     END IF
     
     IF( NeedToGenerate ) THEN
       CALL Info(Caller,'Building constraint matrix',Level=12)
     ELSE     
       CALL Info(Caller,'Nothing to do for now',Level=12)
       RETURN
     END IF

     SkipConstrained = ListGetLogical( Solver % Values, 'Skip Already Constrained Dofs', Found)
     ConstrainedDof => Solver % Matrix % ConstrainedDof
     
     ! Compute the number and size of initial constraint matrices
     !-----------------------------------------------------------
     row    = 0
     mcount = 0
     bcount = 0
     Ctmp => Solver % ConstraintMatrix
     IF( ASSOCIATED( Ctmp ) ) THEN
       DO WHILE(ASSOCIATED(Ctmp))
         mcount = mcount + 1
         row = row + Ctmp % NumberOfRows
         Ctmp => Ctmp % ConstraintMatrix

         IF( InfoActive(32) .AND. ASSOCIATED(Ctmp) ) THEN
           CALL VectorValuesRange(Ctmp % Values,SIZE(Ctmp % Values),'Ctmp'//I2S(mcount))
           IF( ASSOCIATED( Ctmp % InvPerm ) ) THEN
             PRINT *,'InvPerm range:',MINVAL(Ctmp % InvPerm), MAXVAL(Ctmp % InvPerm), SUM(Ctmp % InvPerm)
           END IF
         END IF
    
       END DO
       CALL Info(Caller,'Number of initial constraint matrices: '//I2S(mcount),Level=12)       
       CALL Info(Caller,'Number of rows in constraint matrices: '//I2S(row),Level=20)       
     END IF
       
     
     ! Compute the number and size of mortar matrices
     !-----------------------------------------------
     IF( ASSOCIATED( Solver % MortarBCs ) ) THEN
       DO bc_ind=1,Model % NumberOFBCs + Model % NumberOfBodyForces
         Atmp => Solver % MortarBCs(bc_ind) % Projector         
         IF( .NOT. ASSOCIATED( Atmp ) ) CYCLE
         IF( Atmp % ProjectorType == PROJECTOR_TYPE_NITSCHE ) CYCLE
         bcount = bcount + 1
         row = row + Atmp % NumberOfRows
       END DO
       CALL Info(Caller,'Number of mortar matrices: '//I2S(bcount),Level=12)       
     END IF
     
     IF( row==0 ) THEN
       CALL Info(Caller,'Nothing to do since there are no constrained dofs!',Level=12)       
       RETURN
     END IF

     MortarDiag = ListGetCReal( Solver % Values,'Mortar Diag',HaveMortarDiag )
     LumpedDiag = ListGetLogical( Solver % Values,'Lumped Diag',Found )

     IF( HaveMortarDiag ) THEN
       CALL Info(Caller,'Adding diagonal entry to mortar constraint!',Level=12)              
     END IF
     
     IF( mcount == 1 .AND. bcount == 0 .AND. .NOT. HaveMortarDiag ) THEN
       CALL Info(Caller,'Using initial constraint matrix',Level=12)       
       Solver % Matrix % ConstraintMatrix => Solver % ConstraintMatrix
       RETURN
     END IF

     ! Now we are generating something more complex and different than last time
     IF( ASSOCIATED( Solver % Matrix % ConstraintMatrix ) ) THEN
       IF ( ListGetLogical( Solver % Values,'Apply Contact BCs', Found ) ) THEN
         CALL Info(Caller,'Remember the previous InvPerm for contact mechanics',Level=20)
         ALLOCATE( PrevInvPerm( SIZE( Solver % Matrix % ConstraintMatrix % InvPerm ) ) )
         PrevInvPerm = Solver % Matrix % ConstraintMatrix % InvPerm       
       END IF
       
       CALL Info(Caller,'Releasing previous constraint matrix',Level=12)     
       CALL FreeMatrix(Solver % Matrix % ConstraintMatrix)
       Solver % Matrix % ConstraintMatrix => NULL()
     END IF
       
     EpsVal = ListGetConstReal( Solver % Values,&
         'Minimum Projector Value', Found )
     IF(.NOT. Found ) EpsVal = 1.0d-8
     
     
     SumProjectors = ListGetLogical( Solver % Values,&
         'Mortar BCs Additive', Found )
     IF( .NOT. Found ) THEN
       IF( bcount > 1 .AND. ListGetLogical( Solver % Values, &
           'Eliminate Linear Constraints',Found ) ) THEN
         CALL Info(Caller,'Enforcing > Mortar BCs Additive < to True to enable elimination',Level=8)
         SumProjectors = .TRUE.
       END IF       
       IF( .NOT. SumProjectors .AND. ListGetLogical( Solver % Values, &
           'Apply Conforming BCs',Found ) ) THEN
         CALL Info(Caller,'Enforcing > Mortar BCs Additive < to True because of conforming BCs',Level=8)
         SumProjectors = .TRUE.
       END IF
     END IF
     EliminatedRows = 0

     CALL Info(Caller,'There are '&
         //I2S(row)//' initial rows in constraint matrices',Level=10)

     dim = Mesh % MeshDim              
     dofs = Solver % Variable % DOFs
     
     Perm => Solver % Variable % Perm
     permsize = SIZE( Perm )
     maxperm  = MAXVAL( Perm )
     AllocationsDone = .FALSE.
     arows = Solver % Matrix % NumberOfRows

     ! Create a table that shows one way how continuous nodal dofs maps to
     ! DG nodal dofs. Only one is needed since we assume reduced basis!
     IF( IsDG ) THEN
       ALLOCATE(DgSome(permsize))
       DgSome = 0

       DO i=1,Mesh % NumberOfBulkElements
         Element => Mesh % Elements(i)
         DO j=1,Element % TYPE % NumberOfNodes
           k = Element % NodeIndexes(j)
           k2 = Element % DGIndexes(j)
           DgSome(k) = k2
         END DO
       END DO       
     END IF

       
     ComplexMatrix = Solver % Matrix % Complex
     ComplexSumRow = .FALSE.
     
     IF( ComplexMatrix ) THEN
       IF( MODULO( Dofs,2 ) /= 0 ) CALL Fatal(Caller,&
           'Complex matrix should have even number of components!')
     ELSE
       ! Currently complex matrix is enforced if there is an even number of 
       ! entries since it seems that we cannot rely on the flag to be set.
       ComplexMatrix = ListGetLogical( Solver % Values,'Linear System Complex',Found )
       IF( .NOT. Found ) ComplexMatrix = ( Dofs == 2*dim) 
     END IF

     IF( ComplexMatrix ) THEN
       IF(dofs==dim .OR. dofs == 2) THEN
         cdofs = dofs
       ELSE
         CALL Fatal(Caller,'Invalid number of dofs for field: '//I2S(dofs))
       END IF
     ELSE     
       IF(dofs==dim .OR. dofs == 1) THEN
         cdofs = dofs
       ELSE IF(dofs==dim+1) THEN
         ! For contact mechanics we want to ignore the pressure. 
         IF( ListGetLogical( Solver % Values,'Apply Contact BCs',Found ) ) THEN
           cdofs = dim
         ELSE
           cdofs = dofs
         END IF
       ELSE
         CALL Fatal(Caller,'Invalid number of dofs for field: '//I2S(dofs))
       END IF
     END IF
     
     ALLOCATE( ActiveComponents(dofs), SetDefined(dofs), rsum(dofs) ) 
     
     IF( SumProjectors ) THEN
       ALLOCATE( SumPerm( dofs * permsize ) )
       SumPerm = 0
       ALLOCATE( SumCount( arows ) )
       SumCount = 0
     END IF
          
     AnyPriority = ListCheckPresentAnyBC( Model,'Projector Priority') 
     IF( AnyPriority ) THEN
       IF(.NOT. SumProjectors ) THEN
         CALL Warn(Caller,'Priority has effect only in additive mode!')
         AnyPriority = .FALSE.
       ELSE
         CALL Info(Caller,'Using priority for projector entries',Level=7)
         ALLOCATE( BCPriority(Model % NumberOfBCs), BCOrdering( Model % NumberOfBCs) )
         BCPriority = 0; BCOrdering = 0
         DO bc_ind=1, Model % NumberOFBCs
           Priority = ListGetInteger( Model % BCs(bc_ind) % Values,'Projector Priority',Found)
           BCPriority(bc_ind) = -bc_ind + Priority * Model % NumberOfBCs 
           BCOrdering(bc_ind) = bc_ind
         END DO
         CALL SortI( Model % NumberOfBCs, BCPriority, BCOrdering )
       END IF
     END IF
     NeglectedRows = 0


100  sumrow = 0
     k2 = 0
     rowoffset = 0
     Priority = -1
     PrevPriority = -1
     sumrow0 = 0
     k20 = 0
     
     TransposePresent = .FALSE.
     Ctmp => Solver % ConstraintMatrix

     DO constraint_ind = Model % NumberOFBCs + Model % NumberOfBodyForces + mcount,1,-1
       
       ! This is the default i.e. all components are applied mortar BCs
       ActiveComponents = .TRUE.
       ThisIsRobin = .FALSE.
       
       IF(constraint_ind > Model % NumberOfBCs + Model % NumberOfBodyForces ) THEN
         ThisIsMortar = .FALSE.
         Reorder = .FALSE.                
         SumThis = .FALSE.
         Atmp => Ctmp
         IF( .NOT. ASSOCIATED( Atmp ) ) CYCLE
         Ctmp => Ctmp % ConstraintMatrix
         IF( .NOT. ASSOCIATED( Atmp % InvPerm ) ) THEN
           IF(.NOT. AllocationsDone ) THEN
             CALL Warn(Caller,'InvPerm is expected, using identity!')
           END IF
         END IF
         CALL Info(Caller,'Adding initial constraint matrix: '&
             //I2S(constraint_ind - Model % NumberOfBCs),Level=8)         
       ELSE
         ThisIsMortar = .TRUE.
         
         ! Assume the mortar matrices refer to unordered mesh dofs
         ! and existing ConstraintMatrix to already ordered entities. 
         Reorder = ThisIsMortar

         SumThis = SumProjectors
         IsBodyForce = .FALSE.
         IF( constraint_ind > Model % NumberOfBCs ) THEN
           IsBodyForce = .TRUE.
           bc_ind = constraint_ind
         ELSE IF( AnyPriority ) THEN           
           bc_ind = BCOrdering(constraint_ind)
         ELSE
           bc_ind = constraint_ind 
         END IF
         
         MortarBC => Solver % MortarBCs(bc_ind) 
         Atmp => MortarBC % Projector
         
         ! Add the number of rows already populated in the constraint matrix so we can associate
         ! the single constraint to the correct entry in the constraint matrix.
         MortarBC % rowoffset = sumrow

         IF( .NOT. ASSOCIATED( Atmp ) ) CYCLE
         IF( Atmp % ProjectorType == PROJECTOR_TYPE_NITSCHE ) CYCLE

         IF(IsBodyForce ) THEN
           BC => Model % BodyForces(bc_ind - Model % NumberOfBCs) % Values
           Priority = 0
         ELSE
           BC => Model % BCs(bc_ind) % Values         
           IF( AnyPriority ) THEN
             Priority = ListGetInteger( BC,'Projector Priority',Found)
           END IF
         END IF
           
         IntegralBC = ListGetLogical( BC,'Integral BC',Found ) 

         IF( Atmp % ProjectorType == PROJECTOR_TYPE_ROBIN ) THEN
           ThisIsRobin = .TRUE.
           Reorder = .FALSE.         
         ELSE
           IF( .NOT. ASSOCIATED( Atmp % InvPerm ) ) THEN
             CALL Fatal(Caller,'InvPerm is required for geometric projector!')
           END IF
         END IF

         IF(.NOT. AllocationsDone ) THEN
           IF( ThisIsRobin ) THEN
             CALL Info(Caller,'Adding flux constraint for BC: '//I2S(bc_ind),Level=8)
           ELSE IF( IntegralBC ) THEN
             CALL Info(Caller,'Adding integral constraint for BC: '//I2S(bc_ind),Level=8)
           ELSE             
             CALL Info(Caller,'Adding mortar projector for BC: '//I2S(bc_ind),Level=8)
           END IF
           CALL Info(Caller,'Adding projector rows: '//I2S(Atmp % NumberOfRows),Level=12)
         END IF           
         
         ! Enable that the user can for vector valued cases either set some 
         ! or skip some field components. 
         SomeSet = .FALSE.
         SomeSkip = .FALSE.
         DO i=1,cDofs
           IF( cDofs > 1 ) THEN
             str = ComponentName( Solver % Variable, i )
           ELSE
             str = TRIM(Solver % Variable % Name)
           END IF

           IF( IntegralBC ) THEN
             SetVal(i) = ListGetCReal( BC,'Integral BC '//TRIM(str),SetDof )
             IF(ThisIsRobin) SetDof = .TRUE.
           ELSE
             SetDof = ListGetLogical( BC,'Mortar BC '//TRIM(str),Found )
           END IF

           SetDefined(i) = Found
           IF(Found) THEN
             ActiveComponents(i) = SetDof
             IF( SetDof ) THEN
               SomeSet = .TRUE.
             ELSE
               SomeSkip = .TRUE.
             END IF
           END IF
         END DO
         
         ! By default all components are applied mortar BC and some are turned off.
         ! If the user does the opposite then the default for other components is True.
         IF( SomeSet .AND. .NOT. ALL(SetDefined(1:cdofs)) ) THEN
           IF( SomeSkip ) THEN
             CALL Fatal(Caller,'Do not know what to do with all '//I2S(cdofs)//' components')
           ELSE
             CALL Info(Caller,'Unspecified components will not be set for BC '//I2S(bc_ind),Level=10)
             DO i=1,cDofs
               IF( .NOT. SetDefined(i) ) ActiveComponents(i) = .FALSE.
             END DO
           END IF
         END IF
       END IF

       TransposePresent = TransposePresent .OR. ASSOCIATED(Atmp % Child)
       IF( TransposePresent ) THEN
         IF(ASSOCIATED(Atmp % Child) ) THEN
           IF(.NOT. ASSOCIATED(Atmp % Child % Values)) THEN
             CALL Fatal(Caller,'Atmp has Child which has no values!')
           END IF
         END IF
         CALL Info(Caller,'Transpose matrix is present',Level=8)
       END IF

       ! If the projector is of type x_s=P*x_m then generate a constraint matrix
       ! of type [D-P]x=0 where D is diagonal unit matrix. 
       CreateSelf = ( Atmp % ProjectorType == PROJECTOR_TYPE_NODAL ) 
       
       IF( SumThis .AND. CreateSelf ) THEN
         CALL Fatal(Caller,'It is impossible to sum up nodal projectors!')
       END IF

       ComplexSumRow = ListGetLogical( Solver % Values,'Complex Sum Row ', Found )
       IF(.NOT. Found ) THEN       
         ComplexSumRow = ( dofs == 2 .AND. ComplexMatrix .AND. .NOT. CreateSelf .AND. &
             SumThis .AND. .NOT. (ASSOCIATED( MortarBC % Diag ) .OR. HaveMortarDiag ) )
       END IF
       
       
       ! We deal with the Robin Flux cBC's here even though they would be associated 
       ! to vector or complex valued field. 
       IF( Dofs == 1 .OR. ThisIsRobin ) THEN         
         IF( .NOT. ActiveComponents(1) ) CYCLE
         CALL AddScalarConstraint()                  
       ELSE IF( ComplexSumRow ) THEN
         CALL AddComplexConstraint()        
       ELSE
         CALL AddVectorConstraint()
       END IF 
       
       IF( .NOT. SumThis ) THEN
         rowoffset = rowoffset + Arows
         IF( SumProjectors ) THEN
           CALL Info(Caller,'Not summed up size is ' &
           //I2S(sumrow)//' rows and '//I2S(k2)//' nonzeros',Level=8)
           sumrow0 = sumrow
           k20 = k2
         END IF
       END IF
         
       PrevPriority = Priority 
     END DO ! constrain_ind

     IF( k2 == 0 ) THEN
       CALL Info(Caller,'No entries in constraint matrix!',Level=8)
!      Solver % Matrix % ConstraintMatrix => NULL()
       RETURN
     END IF

     ! Allocate the united matrix of all the boundary matrices
     !-------------------------------------------------------
     IF( .NOT. AllocationsDone ) THEN
       CALL Info(Caller,'Allocating '//&
           I2S(sumrow)//' rows and '//I2S(k2)//' nonzeros',Level=8)

       IF( ComplexSumRow ) THEN
         sumrow = 2 * sumrow
       END IF
       
       Btmp => AllocateMatrix()
       ALLOCATE( Btmp % RHS(sumrow), Btmp % Rows(sumrow+1), &
           Btmp % Cols(k2), Btmp % Values(k2), &
           Btmp % InvPerm(sumrow) )

       Btmp % Rhs = 0.0_dp
       Btmp % Rows = 0
       Btmp % Cols = 0
       Btmp % Values = 0.0_dp
       Btmp % NumberOFRows = sumrow 
       Btmp % InvPerm = 0
       Btmp % Rows(1) = 1

       IF(TransposePresent) THEN
         ALLOCATE(Btmp % TValues(k2))
         Btmp % Tvalues = 0._dp
       END IF

       IF( SumProjectors ) THEN
         Btmp % Rows(sumrow0+1) = k20+1 
         DO i=sumrow0+2,sumrow+1
           Btmp % Rows(i) = Btmp % Rows(i-1) + SumCount(i-1)
         END DO
         SumPerm = 0
         DEALLOCATE( SumCount ) 
       END IF

       AllocationsDone = .TRUE.

       GOTO 100
     END IF
     
     CALL Info(Caller,'Used '//I2S(sumrow)//&
         ' rows and '//I2S(k2)//' nonzeros',Level=7)

     ! Eliminate entries
     IF( SumProjectors ) THEN
       CALL Info(Caller,'Number of eliminated rows: '//I2S(EliminatedRows),Level=6)
       IF( EliminatedRows > 0 ) CALL CRS_PackMatrix( Btmp ) 
     END IF

     IF( NeglectedRows > 0 ) THEN
       CALL Info(Caller,'Number of neglected rows: '//I2S(NeglectedRows),Level=6)
     END IF

     i = COUNT(Btmp % Cols == 0 )
     IF(i>0) CALL Fatal(Caller,'Number of zero Cols in constraint matrix: '//I2S(i))

     
     IF( InfoActive(30) ) THEN
       BLOCK
         REAL(KIND=dp), POINTER :: px(:)
         px => Btmp % Values
         CALL VectorValuesRange(px,SIZE(px),'ConstraintMatrix')
       END BLOCK
     END IF

     Solver % Matrix % ConstraintMatrix => Btmp     
     Solver % MortarBCsChanged = .FALSE.

     IF( InfoActive(20) ) THEN
       WRITE(Message,'(A,ES15.6)') 'Sum of constraint matrix entries: ',SUM(Btmp % Values)
       CALL Info(Caller,Message)
       WRITE(Message,'(A,ES15.6)') 'Abs sum of constraint matrix entries: ',SUM(ABS(Btmp % Values))
       CALL Info(Caller,Message)
       WRITE(Message,'(A,ES15.6)') 'Sum of constraint matrix rhs: ',SUM(Btmp % Rhs)
       CALL Info(Caller,Message)
       CALL Info(Caller,'Constraint matrix cols min:'//I2S(MINVAL(Btmp%Cols)))
       CALL Info(Caller,'Constraint matrix cols max:'//I2S(MAXVAL(Btmp%Cols)))
       CALL Info(Caller,'Constraint matrix rows min:'//I2S(MINVAL(Btmp%Rows)))
       CALL Info(Caller,'Constraint matrix rows max:'//I2S(MAXVAL(Btmp%Rows)))
     END IF

     ! For contact mechanics the number of lagrange multipliers may change.
     ! Hence redistribute the old values to the new initial guess using the InvPerm
     ! to identify the correct location.
     !---------------------------------------------------------------------------------     
     IF ( ListGetLogical( Solver % Values,'Apply Contact BCs', Found ) .AND. &
         ALLOCATED( PrevInvPerm ) ) THEN
       MultName = LagrangeMultiplierName( Solver, SetUnfound = .TRUE. ) 
       Var => VariableGet(Solver % Mesh % Variables, MultName)
       IF( ASSOCIATED( Var ) ) THEN
         ALLOCATE( PrevValues( SIZE( Var % Values ) ) )
         PrevValues = Var % Values
         
         k = 0
         l = SIZE(Btmp % InvPerm) 
         Var % Values = 0.0_dp
         
         DO i=1,l
           DO j=1,SIZE(PrevInvPerm)
             IF( Btmp % InvPerm(i) == PrevInvPerm(j) ) THEN
               k = k + 1
               Var % Values(i) = PrevValues(j)
               EXIT
             END IF
           END DO
         END DO

         CALL Info(Caller,'Previous Lagrange multipliers utilized: '&
             //I2S(k)//' out of '//I2S(l),Level=8)
       END IF
     END IF

     CALL Info(Caller,'Finished creating constraint matrix',Level=12)

   CONTAINS

     SUBROUTINE AddScalarConstraint()
       IF( .NOT. ActiveComponents(1) ) THEN
         CALL Info(Caller,'Skipping component: '//I2S(1),Level=12)
         RETURN
         !CYCLE
       END IF

       ! Number the rows. 
       IF( SumThis ) THEN
         DO i=1,Atmp % NumberOfRows                               
           ! Skip empty row
           IF( Atmp % Rows(i) >= Atmp % Rows(i+1) ) CYCLE 

           ! If the mortar boundary is not active at this round don't apply it
           IF( ThisIsMortar ) THEN
             IF( ASSOCIATED( MortarBC % Active ) ) THEN
               IF( .NOT. MortarBC % Active(i) ) CYCLE
             END IF
           END IF

           ! Use InvPerm if it is present
           IF( ASSOCIATED( Atmp % InvPerm ) ) THEN
             k = Atmp % InvPerm(i)
             ! Node does not have an active dof to be constrained
             IF( k == 0 ) CYCLE
           ELSE
             k = i
           END IF

           kk = k             
           IF( Reorder ) THEN
             IF(IsDG) THEN
               kk = Perm(DgSome(k))
             ELSE
               kk = Perm(k)
             END IF
             IF( kk == 0 ) CYCLE
           END IF

           NewRow = ( SumPerm(kk) == 0 )
           IF( NewRow ) THEN
             sumrow = sumrow + 1                
             SumPerm(kk) = sumrow 
           ELSE IF(.NOT. AllocationsDone ) THEN
             IF( Priority /= PrevPriority .AND. SumPerm(kk) < 0 ) THEN
               NeglectedRows = NeglectedRows + 1
             ELSE
               EliminatedRows = EliminatedRows + 1
             END IF
           END IF
         END DO
         CALL Info(Caller,'Number of rows: '//I2S(sumrow),Level=20)         
       END IF

       IF( ASSOCIATED( MortarBC % Diag ) .OR. HaveMortarDiag) THEN
         CALL Info(Caller,'MotarBC diag exists!',Level=30)
         IF( ASSOCIATED(Atmp % InvPerm) ) THEN
           CALL Info(Caller,'MotarBC InvPerm exists!',Level=30)
           IF( .NOT. ASSOCIATED( MortarBC % Perm ) ) THEN                   
             k = MAXVAL( Atmp % Cols )
             ALLOCATE( MortarBC % Perm(k) )
             MortarBC % Perm = 0
             DO k=1,SIZE(Atmp % InvPerm )
               j = Atmp % InvPerm(k)
               MortarBC % Perm( j ) = k
             END DO
           END IF
         END IF
       END IF


       DO i=1,Atmp % NumberOfRows                     

         IF( Atmp % Rows(i) >= Atmp % Rows(i+1) ) CYCLE ! skip empty rows

         ! If the mortar boundary is not active at this round don't apply it
         IF( ThisIsMortar ) THEN
           IF( ASSOCIATED( MortarBC % Active ) ) THEN
             IF( .NOT. MortarBC % Active(i) ) CYCLE
           END IF
         END IF

         IF( ASSOCIATED( Atmp % InvPerm ) ) THEN
           k = Atmp % InvPerm(i)
           IF( k == 0 ) CYCLE
         ELSE
           k = i
         END IF

         kk = k
         IF( Reorder ) THEN
           IF(IsDg) THEN
             kk = Perm(DgSome(k)) 
           ELSE
             kk = Perm(k)
           END IF
           IF( kk == 0 ) CYCLE
           IF ( SkipConstrained .AND. ConstrainedDof(kk) ) CYCLE
         END IF

         IF( SumThis ) THEN             
           row = SumPerm(kk)

           ! Mark this for future contributions so we know this is already set
           ! and can skip this above.
           IF( AnyPriority ) THEN
             IF( row < 0 ) CYCLE
             IF( Priority /= PrevPriority ) SumPerm(kk) = -SumPerm(kk)
           END IF

           IF( row <= 0 ) THEN
             CALL Fatal(Caller,'Invalid row index: '//I2S(row))
           END IF
         ELSE
           sumrow = sumrow + 1
           row = sumrow
         END IF

         ! This fixes a problem for Robin type of constraints.
         ! Rather than understanding this just omits the problem...
         IF( AllocationsDone .AND. .NOT. ThisIsRobin ) THEN
           Btmp % InvPerm(row) = rowoffset + kk
         END IF

         wsum = 0.0_dp
         rsum = 0.0_dp

         valsum = 0.0_dp
         DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1             
           valsum = valsum + ABS( Atmp % Values(l) ) 
         END DO

         DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1

           col = Atmp % Cols(l) 
           val = Atmp % Values(l)

           IF( ABS( val ) < EpsVal * valsum ) CYCLE

           IF( Reorder ) THEN               
             IF( col <= permsize ) THEN
               IF(IsDg) THEN
                 col2 = Perm(DgSome(col))
               ELSE
                 col2 = Perm(col)
               END IF
               IF( col2 == 0 ) CYCLE
             ELSE
               CALL Fatal(Caller,'col index too large: '//I2S(col)//' vs '//I2S(permsize))
             END IF
           ELSE
             col2 = col
           END IF

           IF( AllocationsDone ) THEN
             ! By Default there is no scaling
             Scale = 1.0_dp
             IF( ThisIsMortar ) THEN
               IF( CreateSelf ) THEN
                 ! We want to create [D-P] hence the negative sign
                 Scale = MortarBC % MasterScale
                 wsum = wsum + val
               ELSE IF( ASSOCIATED( MortarBC % Perm ) ) THEN
                 ! Look if the component refers to the slave
                 IF( MortarBC % Perm( col ) > 0 ) THEN
                   Scale = MortarBC % SlaveScale 
                   wsum = wsum + val
                 ELSE
                   Scale = MortarBC % MasterScale
                 END IF
               ELSE IF( ThisIsRobin ) THEN
                 j = MODULO(col-1,dofs)+1
                 rsum(j) = rsum(j) + val
               ELSE
                 wsum = wsum + val
               END IF

               ! If we sum up to anti-periodic dof then use different sign
               ! - except if the target is also antiperiodic.
               IF( PerFlipActive ) THEN
                 IF(  PerFlip(col) .NEQV. PerFlip(k) ) Scale = -Scale
               END IF
             END IF

             ! Add a new column index to the summed up row               
             ! At the first sweep we need to find the first unset position
             IF( SumThis ) THEN
               k2 = Btmp % Rows(row)
               DO WHILE( Btmp % Cols(k2) > 0 )
                 k2 = k2 + 1
               END DO
             ELSE
               k2 = k2 + 1
             END IF

             Btmp % Cols(k2) = col2
             Btmp % Values(k2) = Scale * val
             IF(ASSOCIATED(Btmp % TValues)) THEN
               IF(ASSOCIATED(Atmp % Child)) THEN
                 Btmp % TValues(k2) = Scale * Atmp % Child % Values(l)
               ELSE
                 Btmp % TValues(k2) = Scale * val
               END IF
             END IF
           ELSE
             k2 = k2 + 1
             IF( SumThis ) THEN
               SumCount(row) = SumCount(row) + 1
             END IF
           END IF
         END DO

         ! Add the self entry as in 'D'
         IF( CreateSelf ) THEN
           k2 = k2 + 1
           IF( AllocationsDone ) THEN
             IF(IsDG) THEN
               Btmp % Cols(k2) = Perm( DGSome( Atmp % InvPerm(i) ) )
             ELSE
               Btmp % Cols(k2) = Perm( Atmp % InvPerm(i) )
             END IF
             Btmp % Values(k2) = MortarBC % SlaveScale * wsum
           ELSE               
             IF( SumThis) SumCount(row) = SumCount(row) + 1
           END IF
         END IF

         ! Add a diagonal entry if requested. When this is done at the final stage
         ! all the hassle with the right column index is easier.
         IF( ThisIsMortar .OR. ThisIsRobin ) THEN
           diag: IF( ASSOCIATED( MortarBC % Diag ) .OR. HaveMortarDiag ) THEN
             IF( .NOT. HaveMortarDiag ) THEN
               MortarDiag = MortarBC % Diag(i)
               LumpedDiag = MortarBC % LumpedDiag
             END IF

             IF( ThisIsRobin ) THEN
               DO j=1,dofs
                 k2 = k2 + 1
                 IF( AllocationsDone ) THEN
                   Btmp % Cols(k2) = j + arows                      
                   Btmp % Values(k2) = Btmp % Values(k2) - MortarDiag * rsum(j)
                 END IF
               END DO

             ELSE IF( LumpedDiag ) THEN
               k2 = k2 + 1
               IF( AllocationsDone ) THEN
                 Btmp % Cols(k2) = row + arows 
                 Btmp % Values(k2) = Btmp % Values(k2) - MortarDiag * wsum
               ELSE
                 IF( SumThis) SumCount(row) = SumCount(row) + 1
               END IF
             ELSE
               IF(ThisIsMortar .AND. .NOT. ASSOCIATED( MortarBC % Perm ) ) THEN                   
                 CALL Fatal(Caller,'MortarBC % Perm required, try lumped')
               END IF

               DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1                 
                 col = Atmp % Cols(l) 

                 IF( Reorder ) THEN
                   IF( col > permsize ) THEN
                     PRINT *,'col too large',col,permsize
                     CYCLE
                   END IF
                   IF(IsDg) THEN
                     col2 = Perm(DgSome(col))
                   ELSE
                     col2 = Perm(col)
                   END IF
                 ELSE
                   col2 = col
                 END IF
                 IF( col2 == 0 ) CYCLE

                 IF( CreateSelf ) THEN
                   Scale = -MortarBC % MasterScale
                 ELSE
                   IF( MortarBC % Perm( col ) > 0 ) THEN
                     Scale = MortarBC % SlaveScale 
                   ELSE
                     CYCLE                     
                   END IF
                 END IF

                 k2 = k2 + 1
                 IF( AllocationsDone ) THEN                                        
                   IF( SumThis ) THEN
                     l2 = ABS( SumPerm( col2) )
                   ELSE
                     l2 = MortarBC % Perm(col)
                   END IF

                   Btmp % Cols(k2) = l2 + arows + rowoffset
                   Btmp % Values(k2) = Btmp % Values(k2) - val * MortarDiag
                 ELSE
                   IF( SumThis) SumCount(row) = SumCount(row) + 1
                 END IF
               END DO
             END IF
           END IF diag
         END IF

         IF( AllocationsDone ) THEN
           IF( IntegralBC ) THEN
             Btmp % Rhs(row) = SetVal(1)
           ELSE IF( ThisIsMortar ) THEN
             IF( ASSOCIATED( MortarBC % Rhs ) ) THEN
               Btmp % Rhs(row) = Btmp % Rhs(row) + wsum * MortarBC % rhs(i)
             END IF
           END IF

           ! If every component is uniquely summed we can compute the row indexes simply
           IF( .NOT. SumThis ) THEN
             Btmp % Rows(row+1) = k2 + 1
           END IF
         END IF
       END DO

     END SUBROUTINE AddScalarConstraint

     
     SUBROUTINE AddComplexConstraint()
       IF(IsDG) CALL Fatal(Caller,'DG not implemented for complex systems!')

       CALL Info(Caller,'Using simplified complex summing!',Level=8)
       ComplexSumRow = .TRUE.

       ! In case of a vector valued problem create a projector that acts on all 
       ! components of the vector. Otherwise follow the same logic.
       IF( SumThis ) THEN
         DO i=1,Atmp % NumberOfRows                        

           IF( ASSOCIATED( Atmp % InvPerm ) ) THEN
             k = Atmp % InvPerm(i)
             IF( k == 0 ) CYCLE
           ELSE
             k = i
           END IF

           kk = k
           IF( Reorder ) THEN
             kk = Perm(k)
             IF( kk == 0 ) CYCLE
           END IF

           NewRow = ( SumPerm(kk) == 0 )
           IF( NewRow ) THEN
             sumrow = sumrow + 1                
             SumPerm(kk) = sumrow 
           ELSE IF(.NOT. AllocationsDone ) THEN
             EliminatedRows = EliminatedRows + 1
           END IF
         END DO
       END IF


       DO i=1,Atmp % NumberOfRows           

         IF( ASSOCIATED( Atmp % InvPerm ) ) THEN
           k = Atmp % InvPerm(i)
           IF( k == 0 ) CYCLE
         ELSE
           k = i
         END IF

         kk = k
         IF( Reorder ) THEN
           kk = Perm(k) 
           IF( kk == 0 ) CYCLE
         END IF

         IF( SumThis ) THEN             
           row = SumPerm(kk)
         ELSE
           sumrow = sumrow + 1
           row = sumrow
         END IF

         ! For complex matrices 
         IF( AllocationsDone ) THEN
           Btmp % InvPerm(2*row-1) = rowoffset + 2*(kk-1)+1
           Btmp % InvPerm(2*row) = rowoffset + 2*kk
         END IF

         wsum = 0.0_dp


         DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1

           col = Atmp % Cols(l) 
           val = Atmp % Values(l)

           IF( Reorder ) THEN
             col2 = Perm(col)
             IF( col2 == 0 ) CYCLE
           ELSE
             col2 = col
           END IF

           IF( AllocationsDone ) THEN
             ! By Default there is no scaling
             Scale = 1.0_dp
             IF( ThisIsMortar ) THEN
               IF( ASSOCIATED( MortarBC % Perm ) ) THEN
                 ! Look if the component refers to the slave
                 IF( MortarBC % Perm( col ) > 0 ) THEN
                   Scale = MortarBC % SlaveScale 
                   wsum = wsum + val
                 ELSE
                   Scale = MortarBC % MasterScale
                 END IF
               ELSE
                 wsum = wsum + val
               END IF

               ! If we sum up to anti-periodic dof then use different sign
               ! - except if the target is also antiperiodic.
               IF( PerFlipActive ) THEN
                 IF(  PerFlip(col) .NEQV. PerFlip(k) ) Scale = -Scale
               END IF

             END IF

             ! Add a new column index to the summed up row               
             ! At the first sweep we need to find the first unset position
             ! Real part
             IF( SumThis ) THEN
               k2 = Btmp % Rows(2*row-1)
               DO WHILE( Btmp % Cols(k2) > 0 )
                 k2 = k2 + 1
               END DO
             ELSE
               k2 = k2 + 1
             END IF

             Btmp % Cols(k2) = 2 * col2 - 1
             Btmp % Values(k2) = Scale * val

             k2 = k2 + 1
             Btmp % Cols(k2) = 2 * col2
             Btmp % Values(k2) = 0.0

             ! Complex part
             IF( SumThis ) THEN
               k2 = Btmp % Rows(2*row)
               DO WHILE( Btmp % Cols(k2) > 0 )
                 k2 = k2 + 1
               END DO
             ELSE
               k2 = k2 + 1
             END IF

             Btmp % Cols(k2) = 2 * col2 - 1 
             Btmp % Values(k2) = 0.0

             k2 = k2 + 1
             Btmp % Cols(k2) = 2 * col2 
             Btmp % Values(k2) = Scale * val
           ELSE
             k2 = k2 + 4
             IF( SumThis ) THEN
               SumCount(2*row-1) = SumCount(2*row-1) + 2
               SumCount(2*row) = SumCount(2*row) + 2
             END IF
           END IF
         END DO

         IF( AllocationsDone ) THEN
           IF( ThisIsMortar ) THEN
             IF( ASSOCIATED( MortarBC % Rhs ) ) THEN
               Btmp % Rhs(2*row-1) = Btmp % Rhs(2*row-1) + wsum * MortarBC % rhs(i)
             END IF
           END IF
         END IF
       END DO

     END SUBROUTINE AddComplexConstraint


     SUBROUTINE AddVectorConstraint()
       
       IF(IsDG) CALL Fatal(Caller,'DG not implemented for vector systems!')

       ! dofs > 1
       ! In case of a vector valued problem create a projector that acts on all 
       ! components of the vector. Otherwise follow the same logic.
       DO i=1,Atmp % NumberOfRows           

         IF( Atmp % Rows(i) >= Atmp % Rows(i+1) ) CYCLE 

         DO j=1,cDofs

           IF( .NOT. ActiveComponents(j) ) THEN
             CALL Info(Caller,'Skipping component: '//I2S(j),Level=12)
             CYCLE
           END IF

           ! For complex matrices both entries mist be created
           ! since preconditioning benefits from 
           IF( ComplexMatrix ) THEN
             IF( MODULO( j, 2 ) == 0 ) THEN
               j2 = j-1
             ELSE 
               j2 = j+1
             END IF
           ELSE
             j2 = 0
           END IF

           IF( ThisIsMortar ) THEN
             IF( ASSOCIATED( MortarBC % Active ) ) THEN
               IF( .NOT. MortarBC % Active(cDofs*(i-1)+j) ) CYCLE
             END IF
           END IF

           IF( ASSOCIATED( Atmp % InvPerm ) ) THEN
             k = Atmp % InvPerm(i)
             IF( k == 0 ) CYCLE
           ELSE
             k = i
           END IF

           kk = k
           IF( Reorder ) THEN
             kk = Perm(k)
             IF( kk == 0 ) CYCLE
           END IF

           IF( SumThis ) THEN
             IF( cDofs*(k-1)+j > SIZE(SumPerm) ) THEN
               CALL Fatal(Caller,'Index out of range')
             END IF
             NewRow = ( SumPerm(cDofs*(kk-1)+j) == 0 )
             IF( NewRow ) THEN
               sumrow = sumrow + 1                
               IF( Priority /= 0 ) THEN
                 ! Use negative sign to show that this has already been set by priority
                 SumPerm(cDofs*(kk-1)+j) = -sumrow 
               ELSE
                 SumPerm(cDofs*(kk-1)+j) = sumrow 
               END IF
             ELSE IF( Priority /= PrevPriority .AND. SumPerm(cDofs*(kk-1)+j) < 0 ) THEN
               IF(.NOT. AllocationsDone ) THEN
                 NeglectedRows = NeglectedRows + 1
               END IF
               CYCLE
             ELSE
               IF(.NOT. AllocationsDone ) THEN
                 EliminatedRows = EliminatedRows + 1
               END IF
             END IF
             row = ABS( SumPerm(cDofs*(kk-1)+j) )
           ELSE
             sumrow = sumrow + 1
             row = sumrow
           END IF

           IF( AllocationsDone ) THEN
             Btmp % InvPerm(row) = rowoffset + Dofs * ( kk - 1 ) + j
           END IF
             
           wsum = 0.0_dp

           valsum = 0.0_dp
           DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1             
             valsum = valsum + ABS( Atmp % Values(l) ) 
           END DO

           
           DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1             

             col = Atmp % Cols(l)                
             val = Atmp % Values(l)                

             IF( ABS( val ) < EpsVal * valsum ) CYCLE
             
             IF( Reorder ) THEN                 
               IF( col <= permsize ) THEN
                 col2 = Perm(col)
                 IF( col2 == 0 ) CYCLE
               ELSE 
                 PRINT *,'col too large',col,permsize
                 CYCLE
               END IF
             ELSE
               col2 = col
             END IF

             k2 = k2 + 1

             IF( AllocationsDone ) THEN
               Scale = 1.0_dp
               IF( ThisIsMortar ) THEN
                 IF( CreateSelf ) THEN
                   Scale = MortarBC % MasterScale
                   wsum = wsum + val
                 ELSE IF( ASSOCIATED( MortarBC % Perm ) ) THEN
                   IF( MortarBC % Perm(col) > 0 ) THEN
                     Scale = MortarBC % SlaveScale 
                     wsum = wsum + val
                   ELSE
                     Scale = MortarBC % MasterScale
                   END IF
                 END IF

                 ! If we sum up to anti-periodic dof then use different sign
                 ! - except if the target is also antiperiodic.
                 IF( PerFlipActive ) THEN
                   IF( PerFlip(col) .NEQV. PerFlip(k) ) Scale = -Scale
                 END IF
               END IF

               IF(Btmp % Cols(k2) /= 0) CALL Fatal('','b1')
               
               Btmp % Cols(k2) = Dofs * ( col2 - 1) + j
               Btmp % Values(k2) = Scale * val

               IF(Btmp % Cols(k2) == 0) CALL Fatal('','zero k2')


               IF(ASSOCIATED(Btmp % Tvalues)) THEN
                 IF(ASSOCIATED(Atmp % Child)) THEN
                   Btmp % TValues(k2) = Scale * Atmp % Child % Values(l)
                 ELSE
                   Btmp % TValues(k2) = Scale * val
                 END IF
               END IF
             ELSE
               IF( SumThis ) THEN
                 SumCount(row) = SumCount(row) + 1
               END IF
             END IF
           END DO

           ! Add the self entry as in 'D'
           IF( CreateSelf ) THEN
             k2 = k2 + 1
             IF( AllocationsDone ) THEN

               IF(Btmp % Cols(k2) /= 0) CALL Fatal('','b2')

               Btmp % Cols(k2) = Dofs * ( Perm( Atmp % InvPerm(i) ) -1 ) + j
               Btmp % Values(k2) = MortarBC % SlaveScale * wsum

               IF(Btmp % Cols(k2) == 0) CALL Fatal('','zero k22')

             END IF
           END IF

           ! Create the imaginary part (real part) corresponding to the 
           ! real part (imaginary part) of the projector. 
           IF( j2 /= 0 ) THEN
             DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1             

               col = Atmp % Cols(l)                

               IF( Reorder ) THEN
                 IF( col <= permsize ) THEN
                   col2 = Perm(col)
                   IF( col2 == 0 ) CYCLE
                 END IF
               ELSE
                 col2 = col
               END IF

               k2 = k2 + 1
               IF( AllocationsDone ) THEN
                 IF(Btmp % Cols(k2) /= 0) CALL Fatal('','b3')


                 Btmp % Cols(k2) = Dofs * ( col2 - 1) + j2
                 IF(Btmp % Cols(k2) == 0) CALL Fatal('','zero k223')

               ELSE
                 IF( SumThis ) THEN
                   SumCount(row) = SumCount(row) + 1
                 END IF
               END IF
             END DO

             IF( CreateSelf ) THEN
               k2 = k2 + 1
               IF( AllocationsDone ) THEN
                 IF(Btmp % Cols(k2) /= 0) CALL Fatal('','b4')


                 Btmp % Cols(k2) = Dofs * ( Perm( Atmp % InvPerm(i) ) -1 ) + j2

                 IF(Btmp % Cols(k2) == 0) CALL Fatal('','zero k24')

               END IF
             END IF
           END IF  ! ComplexMatrix


           IF( ThisIsMortar ) THEN
             IF( ASSOCIATED( MortarBC % Diag ) .OR. HaveMortarDiag ) THEN
               IF( .NOT. HaveMortarDiag ) THEN
                 MortarDiag = MortarBC % Diag(cDofs*(i-1)+j)
                 LumpedDiag = MortarBC % LumpedDiag
               END IF

               IF( LumpedDiag ) THEN
                 k2 = k2 + 1
                 IF( AllocationsDone ) THEN
                   IF(Btmp % Cols(k2) /= 0) CALL Fatal('','b5')


                   Btmp % Cols(k2) = row + arows
                   Btmp % Values(k2) = -wsum * MortarDiag
                   IF(Btmp % Cols(k2) == 0) CALL Fatal('','zero k25')
                 END IF

                 
               ELSE
                 DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1                 
                   col = Atmp % Cols(l) 

                   IF( col > permsize ) CYCLE
                   col2 = Perm(col)

                   IF( CreateSelf ) THEN
                     Scale = -MortarBC % MasterScale
                   ELSE 
                     IF( MortarBC % Perm( col ) > 0 ) THEN
                       Scale = MortarBC % SlaveScale 
                     ELSE
                       CYCLE                     
                     END IF
                   END IF

                   k2 = k2 + 1
                   IF( AllocationsDone ) THEN                   
                     IF(Btmp % Cols(k2) /= 0) CALL Fatal('','b6')


                     Btmp % Cols(k2) = Dofs*(MortarBC % Perm( col )-1) + j + arows + rowoffset
                     IF(Btmp % Cols(k2) == 0) CALL Fatal('','zero k26')

                     Btmp % Values(k2) = -Atmp % Values(l) * MortarDiag
                   END IF
                 END DO
               END IF
             END IF
           END IF


           IF( AllocationsDone ) THEN
             IF( IntegralBC ) THEN
               Btmp % Rhs(row) = SetVal(j)
             ELSE IF( ThisIsMortar ) THEN
               IF( ASSOCIATED( MortarBC % Rhs ) ) THEN
                 Btmp % Rhs(row) = wsum * MortarBC % rhs(cDofs*(i-1)+j)
               END IF
             END IF
             IF(.NOT. SumThis ) THEN
               Btmp % Rows(row+1) = k2 + 1
             END IF
           END IF

         END DO
       END DO
         
     END SUBROUTINE AddVectorConstraint     
     
   END SUBROUTINE GenerateConstraintMatrix
     

   SUBROUTINE ReleaseConstraintMatrix(Solver) 
     TYPE(Solver_t) :: Solver

     CALL FreeMatrix(Solver % Matrix % ConstraintMatrix)
     Solver % Matrix % ConstraintMatrix => NULL()

   END SUBROUTINE ReleaseConstraintMatrix



   ! Generate add matrix from mortar projectors of Nistche type. 
   ! This routine takes each boundary projector and applies it 
   ! to the current field variable (scalar or vector) merging 
   ! all into one single matrix. 
   !---------------------------------------------------------
   SUBROUTINE GenerateAddMatrix( Model, Solver )

     TYPE(Model_t) :: Model
     TYPE(Solver_t), TARGET :: Solver

     INTEGER, POINTER :: Perm(:)
     INTEGER :: i,j,k,k2,l,dofs,permsize,bc_ind,row,col,col2,&
         bcount,kk,cdofs,dim
     TYPE(Matrix_t), POINTER :: Atmp,Btmp 
     LOGICAL :: Found, ComplexMatrix, SomeSet, SomeSkip, SetDof
     LOGICAL, ALLOCATABLE :: ActiveComponents(:), SetDefined(:)
     TYPE(ValueList_t), POINTER :: BC
     TYPE(MortarBC_t), POINTER :: MortarBC
     REAL(KIND=dp) :: Scale
     INTEGER :: arows
     LOGICAL :: Reorder
     LOGICAL :: PerFlipActive, SkipConstrained
     LOGICAL, POINTER :: ConstrainedDof(:)             
     REAL(KIND=dp) :: val, valsum, EpsVal
     LOGICAL, POINTER :: PerFlip(:)
     CHARACTER(:), ALLOCATABLE :: Str 
     LOGICAL :: IsDg, IsBodyForce
     INTEGER, ALLOCATABLE :: DgSome(:)
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Element_t), POINTER :: Element
     
     CHARACTER(*), PARAMETER :: Caller = 'GenerateAddMatrix'

     
     ! Should we genarete the matrix
     IF(.NOT. Solver % MortarBCsChanged ) THEN
       CALL Info(Caller,'Nothing to do for now',Level=20)
       RETURN
     END IF

     IF(.NOT. ASSOCIATED( Solver % MortarBCs ) ) THEN
       CALL Fatal(Caller,'Mortar BCs changed but structure not present!?')
     END IF

     ! Compute the number and size of mortar matrices of type Nitshce
     !---------------------------------------------------------------
     row = 0
     bcount = 0
     DO bc_ind=1,Model % NumberOFBCs + Model % NumberOfBodyForces
       Atmp => Solver % MortarBCs(bc_ind) % Projector         
       IF( .NOT. ASSOCIATED( Atmp ) ) CYCLE
       IF( Atmp % ProjectorType == PROJECTOR_TYPE_NITSCHE ) THEN
         bcount = bcount + 1
         row = row + Atmp % NumberOfRows
       END IF
     END DO
     IF( row==0 ) RETURN
     CALL Info(Caller,'Number of Nitsche matrices: '//I2S(bcount),Level=12)       
     CALL Info(Caller,'Number of primary rows: '//I2S(row),Level=12)       
     
     ! Set pointers to save the initial constraint matrix
     ! We assume that the given ConstraintMatrix is constant but we have consider it the 1st time
     IF( ASSOCIATED( Solver % Matrix % AddMatrix ) ) THEN
       CALL Info(Caller,'Releasing previous AddMatrix!',Level=12)
       CALL FreeMatrix( Solver % Matrix % AddMatrix ) 
     END IF
     
     Mesh => Solver % Mesh 
     IsDg = Solver % DG     

     PerFlipActive = Solver % PeriodicFlipActive
     IF( PerFlipActive ) THEN
       CALL Info(Caller,'Periodic flip is active',Level=8)
       PerFlip => Mesh % PeriodicFlip           
     END IF
          
     SkipConstrained = ListGetLogical( Solver % Values, 'Skip Already Constrained Dofs', Found)
     ConstrainedDof => Solver % Matrix % ConstrainedDof
     
     EpsVal = ListGetConstReal( Solver % Values,&
         'Minimum Projector Value', Found )
     IF(.NOT. Found ) EpsVal = 1.0d-8
          
     dim = Mesh % MeshDim              
     dofs = Solver % Variable % DOFs
     
     Perm => Solver % Variable % Perm
     permsize = SIZE( Perm )
     arows = Solver % Matrix % NumberOfRows

     ! Use list matrix type since it saves us from many headaches. 
     Btmp => AllocateMatrix()
     Btmp % Format = MATRIX_LIST
     CALL AddToMatrixElement( Btmp, arows, arows, 0.0_dp )          

     
     ! Create a table that shows one way how continuous nodal dofs maps to
     ! DG nodal dofs. Only one is needed since we assume reduced basis!
     IF( IsDG ) THEN
       ALLOCATE(DgSome(permsize))
       DgSome = 0

       DO i=1,Mesh % NumberOfBulkElements
         Element => Mesh % Elements(i)
         DO j=1,Element % TYPE % NumberOfNodes
           k = Element % NodeIndexes(j)
           k2 = Element % DGIndexes(j)
           DgSome(k) = k2
         END DO
       END DO       
     END IF
       
     ComplexMatrix = Solver % Matrix % Complex
     IF( ComplexMatrix ) THEN
       IF( MODULO( Dofs,2 ) /= 0 ) CALL Fatal(Caller,&
           'Complex matrix should have even number of components!')
     ELSE
       ! Currently complex matrix is enforced if there is an even number of 
       ! entries since it seems that we cannot rely on the flag to be set.
       ComplexMatrix = ListGetLogical( Solver % Values,'Linear System Complex',Found )
       IF( .NOT. Found ) ComplexMatrix = ( Dofs == 2*dim) 
     END IF

     IF( ComplexMatrix ) THEN
       cdofs = MIN(dofs,2*dim)
     ELSE
       cdofs = MIN(dofs,dim)
     END IF
       
     ALLOCATE( ActiveComponents(cdofs), SetDefined(cdofs) )
     

     
     DO bc_ind = Model % NumberOFBCs + Model % NumberOfBodyForces,1,-1
       
       ! This is the default i.e. all components are applied mortar BCs
       ActiveComponents = .TRUE.
       Reorder = .TRUE.

       IsBodyForce = ( bc_ind > Model % NumberOfBCs )        
       
       MortarBC => Solver % MortarBCs(bc_ind) 
       Atmp => MortarBC % Projector
       IF( .NOT. ASSOCIATED( Atmp ) ) CYCLE

       IF( Atmp % ProjectorType /= PROJECTOR_TYPE_NITSCHE ) CYCLE         

       IF(IsBodyForce ) THEN
         BC => Model % BodyForces(bc_ind - Model % NumberOfBCs) % Values
       ELSE
         BC => Model % BCs(bc_ind) % Values         
       END IF
         
       CALL Info(Caller,'Adding mortar projector of type Nitsche for BC: '//I2S(bc_ind),Level=8)
       
       ! Enable that the user can for vector valued cases either set some 
       ! or skip some field components. 
       SomeSet = .FALSE.
       SomeSkip = .FALSE.
       DO i=1,cDofs
         IF( cDofs > 1 ) THEN
           str = ComponentName( Solver % Variable, i )
         ELSE
           str = TRIM(Solver % Variable % Name)
         END IF
         SetDof = ListGetLogical( BC,'Mortar BC '//TRIM(str),Found )

         SetDefined(i) = Found
         IF(Found) THEN
           ActiveComponents(i) = SetDof
           IF( SetDof ) THEN
             SomeSet = .TRUE.
           ELSE
             SomeSkip = .TRUE.
           END IF
         END IF
       END DO

       ! By default all components are applied mortar BC and some are turned off.
       ! If the user does the opposite then the default for other components is True.
       IF( SomeSet .AND. .NOT. ALL(SetDefined(1:cdofs)) ) THEN
         IF( SomeSkip ) THEN
           CALL Fatal(Caller,'Do not know what to do with all '//I2S(cdofs)//' components')
         ELSE
           CALL Info(Caller,'Unspecified components will not be set for BC '//I2S(bc_ind),Level=10)
           DO i=1,cDofs
             IF( .NOT. SetDefined(i) ) ActiveComponents(i) = .FALSE.
           END DO
         END IF
       END IF

       ! Add one single matrix of type Nitshce to AddMatrix.
       CALL AddNitscheMatrix(Atmp)
     END DO

     
     CALL List_ToCRSMatrix( Btmp )
     k2 = SIZE(Btmp % Values)         
     
     CALL Info(Caller,'Used '//I2S(Btmp % NumberOfRows)//&
         ' rows and '//I2S(k2)//' nonzeros',Level=7)
     
     Solver % Matrix % AddMatrix => Btmp     
     Solver % MortarBCsChanged = .FALSE.

     IF( InfoActive(20) ) THEN
       WRITE(Message,'(A,ES15.6)') 'Sum of add matrix entries: ',SUM(Btmp % Values)
       CALL Info(Caller,Message)
       WRITE(Message,'(A,ES15.6)') 'Abs sum of add matrix entries: ',SUM(ABS(Btmp % Values))
       CALL Info(Caller,Message)
       CALL Info(Caller,'Add matrix cols range: '&
           //I2S(MINVAL(Btmp%Cols))//':'//I2S(MAXVAL(Btmp%Cols)))
       CALL Info(Caller,'Add matrix rows range: '&
           //I2S(MINVAL(Btmp%Rows))//':'//I2S(MAXVAL(Btmp%Rows)))
     END IF

     Solver % Matrix % AddMatrix => Btmp
     
     CALL Info(Caller,'Finished creating add matrix',Level=12)

   CONTAINS

     SUBROUTINE AddNitscheMatrix(Atmp)

       TYPE(Matrix_t) :: Atmp
       
       REAL(KIND=dp), ALLOCATABLE :: Vals(:), Zeros(:)
       INTEGER, ALLOCATABLE :: Cols(:)
       INTEGER :: dofi,jc,n
       INTEGER :: NormalInd(1), NodeCols(3)
       REAL(KIND=dp) :: NodeMatrix(3,3), NodeForce(3)
       TYPE(NormalTangential_t), POINTER :: NT
       LOGICAL :: RotateNT
       
       IF( .NOT. ANY(ActiveComponents(1:dofs)) ) RETURN

       RotateNT = .FALSE.
       NT => Solver % NormalTangential
       IF(ASSOCIATED(NT)) THEN
         RotateNT =  ( NT % NormalTangentialNOFNodes > 0 .AND. dofs > 1)
       END IF
                     
       j = 0
       DO i=1,Atmp % NumberOfRows                     
         j = MAX(j, Atmp % Rows(i+1) - Atmp % Rows(i))
       END DO
       n = j * dofs
       ALLOCATE(Vals(n),Cols(n),Zeros(n))
       Zeros = 0.0_dp


       DO i=1,Atmp % NumberOfRows                     

         IF( Atmp % Rows(i) >= Atmp % Rows(i+1) ) CYCLE ! skip empty rows directly

         ! If the mortar boundary is not active at this round don't apply it.
         ! Every projector here comes from Solver % MortarBCs, so unlike in
         ! GenerateConstraintMatrix there is nothing to tell apart here: the
         ! test is just whether an active set has been given at all.
         IF( ASSOCIATED( MortarBC % Active ) ) THEN
           IF( .NOT. MortarBC % Active(i) ) CYCLE
         END IF

         ! Relate the constraint to geometric entity
         ! We have created Nitsche projector directly using the size of full matrix, not reduced system. 
         
         kk = i
         IF( Reorder ) THEN
           IF(IsDg) THEN
             kk = Perm(DgSome(i)) 
           ELSE
             kk = Perm(i)
           END IF
           IF( kk == 0 ) CYCLE
           IF ( SkipConstrained .AND. ConstrainedDof(kk) ) CYCLE
         END IF

         valsum = 0.0_dp
         DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1             
           valsum = valsum + ABS( Atmp % Values(l) ) 
         END DO

         k2 = 0
         DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1

           col = Atmp % Cols(l) 
           val = Atmp % Values(l)

           IF( ABS( val ) < EpsVal * valsum ) CYCLE

           IF( Reorder ) THEN               
             IF( col <= permsize ) THEN
               IF(IsDg) THEN
                 col2 = Perm(DgSome(col))
               ELSE
                 col2 = Perm(col)
               END IF
               IF( col2 == 0 ) CYCLE
             ELSE
               CALL Fatal(Caller,'col index too large: '//I2S(col)//' vs '//I2S(permsize))
             END IF
           ELSE
             col2 = col
           END IF

           ! By Default there is no scaling
           Scale = 1.0_dp

           ! If we sum up to anti-periodic dof then use different sign
           ! - except if the target is also antiperiodic.
           ! Both indexes are unpermuted mesh dofs, and the row here is i:
           ! GenerateConstraintMatrix has a k of its own for this, this routine
           ! does not.
           IF( PerFlipActive ) THEN
             IF(  PerFlip(col) .NEQV. PerFlip(i) ) Scale = -Scale
           END IF

           k2 = k2+1
           Cols(k2) = col2
           Vals(k2) = Scale * val           
         END DO

         IF(RotateNT) THEN
           DO j=1,k2
             NodeMatrix = 0.0_dp
             NodeForce = 0.0_dp
             DO k=1,dim
               ! We copy the Nitsche entry of the matrix in a vector valued entry. 
               NodeMatrix(k,k) = Vals(k2)
               NodeCols(k) = dofs*(Cols(k2)-1)+j
             END DO

             ! This is not correct since this is not a local square matrix where we would
             ! operate the matrix from left and right WITH same rotation matrix.
             NormalInd(1) = NT % BoundaryReorder(kk)
             IF(NormalInd(1) > 0) THEN             
               CALL RotateMatrix( NodeMatrix, NodeForce, 1, dim, dofs, &
                   NormalInd, NT % BoundaryNormals, NT % BoundaryTangent1, NT % BoundaryTangent2 )
             END IF
             DO k=1,dim
               CALL List_AddMatrixRow(Btmp % ListMatrix,dofs*(kk-1)+i,dim,NodeCols,NodeMatrix(k,1:dim),SortedInput=.TRUE.)
             END DO
           END DO         
         ELSE IF(dofs > 1) THEN           
           Cols(1:k2) = dofs*(Cols(1:k2)-1)           
           DO dofi=1,dofs
             Cols(1:k2) = Cols(1:k2)+1
             IF( .NOT. ActiveComponents(dofi) ) CYCLE
             CALL List_AddMatrixRow(Btmp % ListMatrix,dofs*(kk-1)+dofi,k2,Cols,Vals,SortedInput=.TRUE.)
             ! We should ensure by construction that the complex matrix includes all the entries to allow
             ! all linear solvers. 
             IF(ComplexMatrix) THEN
               IF(MODULO(dofi,2)==0) THEN
                 jc = -1
               ELSE
                 jc = +1
               END IF
               CALL List_AddMatrixRow(Btmp % ListMatrix,dofs*(kk-1)+dofi+jc,k2,Cols,Zeros,SortedInput=.TRUE.)
             END IF
           END DO
         ELSE
           CALL List_AddMatrixRow(Btmp % ListMatrix,kk,k2,Cols,Vals,SortedInput=.TRUE.)
         END IF
       END DO
       
     END SUBROUTINE AddNitscheMatrix
          
   END SUBROUTINE GenerateAddMatrix



   

   SUBROUTINE ReleaseProjectors(Model, Solver) 

     TYPE(Model_t) :: Model
     TYPE(Solver_t) :: Solver

     TYPE(ValueList_t), POINTER :: BC
     TYPE(Matrix_t), POINTER :: Projector
     INTEGER :: i
     

     IF( .NOT. ASSOCIATED( Solver % MortarBCs ) ) RETURN

     DO i=1,Model % NumberOFBCs
       BC => Model % BCs(i) % Values
       Projector => Solver % MortarBCs(i) % Projector 
       IF( ASSOCIATED( Projector ) ) THEN
         IF( ASSOCIATED( Projector % EMatrix ) ) THEN
           CALL FreeMatrix( Projector % Ematrix ) 
         END IF
         CALL FreeMatrix( Projector )
         Solver % MortarBCs(i) % Projector => NULL()
       END IF
     END DO

   END SUBROUTINE ReleaseProjectors

   !> This subroutine saves a projector assuming time-periodic system.
   !> There are two operation modes.
   !> a) Fetching a precomputed projector when GotProj argument is provided.
   !> b) Storing a projector when no GotProj argument is provided. 
   !----------------------------------------------------------------------
   SUBROUTINE StoreCyclicProjector(Solver,Proj,GotProj)
     TYPE ProjTable_t
       TYPE(Matrix_t), POINTER :: Proj
     END TYPE ProjTable_t
     TYPE(Solver_t) :: Solver
     TYPE(Matrix_t), POINTER :: Proj
     LOGICAL, OPTIONAL :: GotProj
     
     TYPE(Variable_t), POINTER :: v
     TYPE(Matrix_t), POINTER :: A     
     TYPE(Model_t), POINTER :: Model
     LOGICAL :: Found
     TYPE(ProjTable_t), POINTER :: ProjTable(:)
     INTEGER :: n, i, Ncycle, Ntime, Nstore, Ntimes
     LOGICAL :: SetProj
     
     SAVE ProjTable
     
     Model => CurrentModel 
     Ncycle = ListGetInteger( Model % Simulation,'Periodic Timesteps')
     Ntimes = ListGetInteger( Model % Simulation,'Number Of Times',Found )
     IF(Found ) Ncycle = Ncycle / Ntimes     

     v => VariableGet( Solver % Mesh % Variables, 'timestep' )
     Ntime = NINT(v % Values(1))

     A => Solver % Matrix
     n = A % NumberOfRows

     ! allocate space for projectors 
     IF(.NOT. ASSOCIATED( ProjTable ) ) THEN
       ALLOCATE( ProjTable(Ncycle) )
       DO i=1,Ncycle
         ProjTable(i) % Proj => NULL()
       END DO
     END IF

     ! Nstrore in [1,Ncycle]
     Nstore = MODULO( Ntime-1,Ncycle)+1

     IF( PRESENT( GotProj ) ) THEN
       ! getting projector
       IF( Ntime <= Ncycle ) THEN
         Proj => NULL()
       ELSE
         Proj => ProjTable(Nstore) % Proj
       END IF
       GotProj = ASSOCIATED( Proj ) 
       IF( InfoActive(20) ) THEN
         PRINT *,'Getting cyclic projector:',GotProj,Ntime,Nstore,Ncycle,ASSOCIATED(Proj)
       END IF
     ELSE
       ! storing projector
       SetProj = .NOT. ASSOCIATED( ProjTable(Nstore) % Proj )       
       IF( SetProj ) ProjTable(Nstore) % Proj => Proj
       IF( InfoActive(20) ) THEN
         PRINT *,'Setting cyclic projector:',SetProj,Ntime,Nstore,Ncycle,ASSOCIATED(Proj)
       END IF
     END IF
         
   END SUBROUTINE StoreCyclicProjector

END MODULE ProjectorUtils
