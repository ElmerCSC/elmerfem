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
!>  Contact mechanics: contact set determination and slave solver merging.
!>  Extracted from SolverUtils.
!------------------------------------------------------------------------------

MODULE ContactUtils

    USE SolverBasics
    USE SolveCore, ONLY : CalculateLoads
    IMPLICIT NONE

CONTAINS

!> Subroutine for determine the contact set and create the necessary data
!> for setting up the contact conditions. As input the mortar projectors,
!> the current solution, and the stiffness matrix are used.  
!------------------------------------------------------------------------------
   SUBROUTINE DetermineContact( Solver )
!------------------------------------------------------------------------------
     TYPE(Solver_t), POINTER :: Solver
!-----------------------------------------------------------------------------
     TYPE(Model_t), POINTER :: Model
     TYPE(variable_t), POINTER :: Var, LoadVar, IterVar
     TYPE(Variable_t), POINTER :: DistVar, NormalLoadVar, SlipLoadVar, VeloVar, &
         WeightVar, NormalActiveVar, StickActiveVar, GapVar, ContactLagrangeVar
     TYPE(Element_t), POINTER :: Element
     TYPE(Mesh_t), POINTER :: Mesh
     INTEGER :: i,j,k,l,n,m,t,ind,dofs,cdofs,dim, bf, Upper, &
         ElemFirst, ElemLast, totsize, i2, j2, ind2, bc_ind, master_ind, &
         DistSign, LimitSign, DofN, DofT1, DofT2, Limited, LimitedMin, TimeStep
     REAL(KIND=dp), POINTER :: FieldValues(:), LoadValues(:), ElemLimit(:),pNormal(:,:),&
         RotatedField(:)
     REAL(KIND=dp) :: ValEps, LoadEps, val, ContactNormal(3), &
         ContactT1(3), ContactT2(3), LocalT1(3), LocalT2(3), &
         LocalNormal(3), NodalForce(3), wsum, coeff, &
         Dist, DistN, DistT1, DistT2, NTT(3,3), RotVec(3), dt
     INTEGER, POINTER :: FieldPerm(:), NodeIndexes(:)
     LOGICAL :: Found,AnyLimitBC, AnyLimitBF
     LOGICAL, ALLOCATABLE :: LimitDone(:),InterfaceDof(:)
     LOGICAL, POINTER :: LimitActive(:)
     TYPE(ValueList_t), POINTER :: Params
     CHARACTER(:), ALLOCATABLE ::  Str,VarName, ContactType
     INTEGER :: ConservativeAfterIters, ActiveDirection, NonlinIter, CoupledIter, EasyIters
     LOGICAL :: ConservativeAdd, ConservativeRemove, &
         DoAdd, DoRemove, DirectionActive, Rotated, FlatProjector, PlaneProjector, &
         RotationalProjector, NormalProjector, FirstTime = .TRUE., &
         AnyRotatedContact, ThisRotatedContact, StickContact, TieContact, FrictionContact, SlipContact, &
         CalculateVelocity, NodalNormal, ResidualMode, AddDiag, SkipFriction, DoIt
     TYPE(MortarBC_t), POINTER :: MortarBC
     TYPE(Matrix_t), POINTER :: Projector, DualProjector
     TYPE(ValueList_t), POINTER :: BC, MasterBC
     REAL(KIND=dp), POINTER :: nWrk(:,:)
     LOGICAL :: CreateDual, pContact
     INTEGER, TARGET :: pIndexes(12)
     TYPE(Variable_t), POINTER :: UseLoadVar
     LOGICAL :: UseLagrange
     TYPE(NormalTangential_t), POINTER :: NT
     CHARACTER(*), PARAMETER :: Caller = 'DetermineContact'

     
     SAVE FirstTime

     CALL Info(Caller,'Setting up contact conditions',Level=8)
     
     Model => CurrentModel
     Var => Solver % Variable
     VarName = GetVarName( Var ) 
     Mesh => Solver % Mesh
     NT => Model % Solver % NormalTangential
     
     ! Is any boundary rotated or not
     AnyRotatedContact = ( NT % NormalTangentialNOFNodes > 0 ) 

     ! The variable to be constrained by the contact algorithm
     ! Here it is assumed to be some "displacement" i.e. a vector quantity
     FieldValues => Var % Values
     FieldPerm => Var % Perm
     totsize = SIZE( FieldValues )
     dofs = Var % Dofs
     dim = Mesh % MeshDim 
     Params => Solver % Values

     IF(dofs == dim) THEN
       cdofs = dim
     ELSE IF(dofs == dim+1) THEN
       CALL Info(Caller,'We seem to have mixed formulation, ignoring pressure!',Level=7)
       cdofs = dim
     ELSE 
       CALL Fatal(Caller,'Invalid number of dofs for contact problem: '//I2S(dofs))
     END IF     

     pContact = IsActivePelement(Mesh % Elements(1), Solver)
     IF( pContact ) THEN
       ! We only have to deal with the middle dofs if they are not condensated away!
       pContact = Solver % GlobalBubbles 
     END IF
     IF( ListGetLogical( Params,'Contact Linear Basis',Found ) ) THEN
       pContact = .FALSE.
     END IF
    
     IF( pContact ) THEN
       CALL Info(Caller,'Using p-elements for contact, if available in projector!',Level=8)
     END IF
     
     IterVar => VariableGet( Model % Variables,'coupled iter',UnfoundFatal=.TRUE.)
     CoupledIter = NINT( IterVar % Values(1) )

     IterVar => VariableGet( Model % Variables,'nonlin iter',UnfoundFatal=.TRUE.)
     NonlinIter = NINT( IterVar % Values(1) )
     
     IterVar => VariableGet( Model % Variables,'timestep',UnfoundFatal=.TRUE.)
     Timestep = NINT( IterVar % Values(1) )

     IterVar => VariableGet( Mesh % Variables,'timestep size')
     IF( ASSOCIATED( IterVar ) ) THEN
       dt = IterVar % Values(1)
     ELSE
       dt = 1.0_dp
     END IF

     !FirstTime = ( NonlinIter == 1 .AND. CoupledIter == 1 ) 

     ConservativeAfterIters = ListGetInteger(Params,&
         'Apply Limiter Conservative Add After Iterations',ConservativeAdd ) 
     IF( ConservativeAdd ) THEN
       IF( CoupledIter == 1 ) ConservativeAdd = ( ConservativeAfterIters < NonlinIter )
       IF( ConservativeAdd ) THEN
         CALL Info(Caller,'Adding dofs in conservative fashion',Level=8)
       END IF
     END IF

     ConservativeAfterIters = ListGetInteger(Params,&
         'Apply Limiter Conservative Remove After Iterations',ConservativeRemove ) 
     IF( ConservativeRemove ) THEN
       IF( CoupledIter == 1 ) ConservativeRemove = ( ConservativeAfterIters < NonlinIter )
       IF( ConservativeRemove ) THEN
         CALL Info(Caller,'Removing dofs in conservative fashion',Level=8)
       END IF
     END IF
         
     EasyIters = ListGetInteger(Params,'Nonlinear System Initial Friction Iterations',Found )
     IF(.NOT. Found ) EasyIters = 1          

     ResidualMode = ListGetLogical(Params,&
         'Linear System Residual Mode',Found )

     CalculateVelocity = ListGetLogical(Params,&
         'Apply Contact Velocity',Found )
     IF(.NOT. Found ) THEN
       Str = ListGetString( CurrentModel % Simulation, 'Simulation Type' )
       CalculateVelocity =  ( Str == 'transient' ) 
     END IF

     NodalNormal = ListGetLogical(Params,&
         'Use Nodal Normal',Found )

     LoadEps = ListGetConstReal(Params,'Limiter Load Tolerance',Found ) 
     IF(.NOT. Found ) LoadEps = EPSILON( LoadEps )
         
     ValEps = ListGetConstReal(Params,'Limiter Value Tolerance',Found ) 
     IF(.NOT. Found ) ValEps = EPSILON( ValEps )

     IF( .NOT. ASSOCIATED( Model % Solver % MortarBCs ) ) THEN
       CALL Fatal(Caller,'Cannot apply contact without projectors!')
     END IF

     ! a) Create rotateted contact if needed
     CALL RotatedDisplacementField() 

     CALL PickLagrangeMultiplier()

     ! b) Create and/or obtain pointers to boundary variables 
     CALL GetContactFields( FirstTime )

     ! c) Calculate the contact loads to the normal direction
     LoadVar => CalculateContactLoad() 
     LoadValues => LoadVar % Values

     UseLagrange = ListGetLogical( Params,'Use Lagrange Multiplier for Contact',Found )
     IF( UseLagrange ) THEN
       CALL Info(Caller,'Using Lagrange multiplier to determine contact condition!')
       UseLoadVar => ContactLagrangeVar 
     ELSE
       UseLoadVar => NormalLoadVar 
     END IF     
     
     ! Loop over each contact pair
     !--------------------------------------------------------------
     DO bc_ind = 1, Model % NumberOfBCs
       
       MortarBC => Model % Solver % MortarBCs(bc_ind)  
       IF( .NOT. ASSOCIATED( MortarBC ) ) CYCLE

       Projector => MortarBC % Projector
       IF(.NOT. ASSOCIATED(Projector) ) CYCLE

       BC => Model % BCs(bc_ind) % Values

       CALL Info(Caller,'Set contact for boundary: '&
           //I2S(bc_ind),Level=8)
       Model % Solver % MortarBCsChanged = .TRUE.
       
       FlatProjector = ListGetLogical( BC, 'Flat Projector',Found ) 
       PlaneProjector = ListGetLogical( BC, 'Plane Projector',Found )
       RotationalProjector = ListGetLogical( BC, 'Rotational Projector',Found ) .OR. &
           ListGetLogical( BC, 'Cylindrical Projector',Found )
       NormalProjector = ListGetLogical( BC, 'Normal Projector',Found )
       
       ! Is the current boundary rotated or not
       ThisRotatedContact = ListGetLogical( BC,'Normal-Tangential '//TRIM(VarName),Found)

       IF( FlatProjector ) THEN
         ActiveDirection = ListGetInteger( BC, 'Flat Projector Coordinate',Found )
         IF( .NOT. Found ) ActiveDirection = cdofs       
       ELSE IF( PlaneProjector ) THEN
         pNormal => ListGetConstRealArray( BC,'Plane Projector Normal',Found)
         IF( ThisRotatedContact ) THEN
           ActiveDirection = 1
         ELSE
           ActiveDirection = 1
           DO i=2,3
             IF( ABS( pnormal(i,1) ) > ABS( pnormal(ActiveDirection,1) ) ) THEN
               ActiveDirection = i
             END IF
           END DO
           CALL Info(Caller,'Active direction set to: '//I2S(ActiveDirection),Level=6)
         END IF
       ELSE IF( RotationalProjector .OR. NormalProjector ) THEN
         ActiveDirection = 1
         IF( .NOT. ThisRotatedContact ) THEN
           CALL Warn(Caller,'Rotational and normal projectors should only work with N-T coordinates!')
         END IF
       ELSE
         CALL Fatal(Caller,'Projector must be current either flat, plane, cylindrical or rotational!')
       END IF
      

       ! Get the pointer to the other side i.e. master boundary  
       master_ind = ListGetInteger( BC,'Mortar BC',Found )
       IF( .NOT. Found ) master_ind = ListGetInteger( BC,'Contact BC',Found )
       MasterBC => Model % BCs(master_ind) % Values
       
       ! If we have dual projector we may use it to map certain quantities directly to master nodes
       DualProjector => Projector % Ematrix
       CreateDual = ASSOCIATED( DualProjector )
       IF( CreateDual ) THEN
         CALL Info(Caller,'Using also the dual projector',Level=8)
       END IF
      
       ! If we have N-T system then the mortar condition for the master side
       ! should have reverse sign as both normal displacement diminish the gap.
       IF( ThisRotatedContact ) THEN
         IF( master_ind > 0 ) THEN
           IF( .NOT. ListGetLogical( MasterBC, &
               'Normal-Tangential '//TRIM(VarName),Found) ) THEN
             CALL Fatal(Caller,'Master boundary '//I2S(master_ind)//&
                 ' should also have N-T coordinates!')
           END IF
         END IF

         CALL Info(Caller,'We have a normal-tangential system',Level=6)
         MortarBC % MasterScale = -1.0_dp
         DofN = 1
       ELSE                 
         DofN = ActiveDirection 
       END IF

       ! Get the degrees of freedom related to the normal and tangential directions
       DofT1 = 0; DofT2 = 0
       DO i=1,cdofs
         IF( i == DofN ) CYCLE
         IF( DofT1 == 0 ) THEN
           DofT1 = i 
           CYCLE
         END IF
         IF( DofT2 == 0 ) THEN
           DofT2 = i
           CYCLE
         END IF
       END DO

       ! This is the normal that is used to detect the signed distance
       ! and tangent vectors used to detect surface velocity
       IF( PlaneProjector ) THEN
         ContactNormal = pNormal(1:3,1)
       ELSE
         ContactNormal = 0.0_dp
         ContactNormal(ActiveDirection) = 1.0_dp
       END IF
       ContactT1 = 0.0_dp
       ContactT1(DofT1) = 1.0_dp
       ContactT2 = 0.0_dp
       IF(DofT2>0) ContactT2(DofT2) = 1.0_dp

       ! Get the contact type. There are four possibilities currently. 
       ! Only one is active at a time while others are false. 
       StickContact = .FALSE.; TieContact = .FALSE.
       FrictionContact = .FALSE.; SlipContact = .FALSE.

       ContactType = ListGetString( BC,'Contact Type',Found ) 
       IF( Found ) THEN
         SELECT CASE ( ContactType )
         CASE('stick')
           StickContact = .TRUE.
         CASE('tie')
           TieContact = .TRUE.
         CASE('friction')
           FrictionContact = .TRUE.
         CASE('slip','slide')
           ! Both spellings: the keyword is "Slip Contact" but the string form
           ! has always been "slide", and the default warning says "slip".
           SlipContact = .TRUE.
         CASE Default
           CALL Fatal(Caller,'Unknown contact type: '//TRIM(ContactType))
         END SELECT
       ELSE
         StickContact = ListGetLogical( BC,'Stick Contact',Found )
         IF(.NOT. Found ) TieContact = ListGetLogical( BC,'Tie Contact',Found )
         IF(.NOT. Found ) FrictionContact = ListGetLogical( BC,'Friction Contact',Found )
         IF(.NOT. Found ) SlipContact = ListGetLogical( BC,'Slip Contact',Found )
         IF(.NOT. Found ) SlipContact = ListGetLogical( BC,'Slide Contact',Found )
         IF(.NOT. Found ) THEN 
           CALL Warn(Caller,'No contact type given, assuming > Slip Contact <')
           SlipContact = .TRUE.
         END IF
       END IF

       IF( StickContact ) CALL Info(Caller,'Using "Stick Contact" for displacement',Level=10)
       IF( TieContact ) CALL Info(Caller,'Using "Tie Contact" for displacement',Level=10)
       IF( FrictionContact ) CALL Info(Caller,'Using "Friction Contact" for displacement',Level=10)
       IF( SlipContact ) CALL Info(Caller,'Using "Slip Contact" for displacement',Level=10)
       

       ! At the start of nonlinear or transient iteration it may be beneficial to assume
       ! less challenging initial contat model
       ! friction models -> stick contact /slip contact -> tie contact
       BLOCK
         LOGICAL :: InitialTie, InitialStick, InitialSlip 

         InitialTie = .FALSE.
         InitialStick = .FALSE.
         InitialSlip = .FALSE.
         
         DO i=1,3
           IF(i==1) THEN
             Str = 'Nonlinear System'
             IF( NonlinIter > EasyIters ) CYCLE
           ELSE IF(i==2) THEN
             Str = 'Transient System'
             IF( TimeStep > 1 ) CYCLE 
           ELSE
             Str = ''
             IF( NonlinIter > 1 .OR. Timestep > 1 ) CYCLE
           END IF

           IF( ListGetLogical( BC,TRIM(Str)//' Initial Tie Contact', Found ) ) InitialTie = .TRUE.
           IF( ListGetLogical( BC,TRIM(Str)//' Initial Stick Contact', Found ) ) InitialStick = .TRUE.
           IF( ListGetLogical( BC,TRIM(Str)//' Initial Slip Contact', Found ) ) InitialSlip = .TRUE.
         END DO
         
         SkipFriction = .FALSE.
         IF( InitialTie .OR. InitialStick .OR. InitialSlip ) THEN
           FrictionContact = .FALSE.
           StickContact = .FALSE.
           SlipContact = .FALSE.
           TieContact = .FALSE.
           SkipFriction = .TRUE.
         END IF

         IF( InitialTie ) THEN
           TieContact = .TRUE.
           CALL Info(Caller,'Reverting to "Tie Contact" for the start',Level=10)                      
         ELSE IF(InitialStick) THEN
           StickContact = .TRUE.
           CALL Info(Caller,'Reverting to "Stick Contact" for the start',Level=10)                      
         ELSE IF(InitialSlip ) THEN
           SlipContact = .TRUE.
           CALL Info(Caller,'Reverting to "Slip Contact" for the start',Level=10)                      
         END IF
       END BLOCK

       ! If we have stick contact then create a diagonal entry to the projection matrix.
       IF( StickContact .OR. FrictionContact ) THEN
         AddDiag = ListCheckPresent( BC,'Stick Contact Coefficient')      
       ELSE
         AddDiag = .FALSE.
       END IF

       IF(InfoActive(30)) THEN
         PRINT *,'Contact Flags:',TieContact, FrictionContact, StickContact, SlipContact, SkipFriction, AddDiag 
       END IF
              
       ! d) allocate and initialize all necessary vectors for the contact 
       !------------------------------------------------------------------
       CALL InitializeMortarVectors()
     
       ! e) If the contact set is set up in a conservative fashion we need to mark interface nodes
       !------------------------------------------------------------------
       IF( ConservativeAdd .OR. ConservativeRemove ) THEN
         CALL MarkInterfaceDofs()
       END IF

       ! f) Compute the normal load used to determine whether contact should be released.
       !    Also check the direction to which the signed distance should be computed
       !------------------------------------------------------------------
       CALL CalculateContactPressure()

       
       ! g) Calculate the distance used to determine whether contact should be added
       !------------------------------------------------------------------
       CALL CalculateMortarDistance()
        
       ! h) Determine the contact set in normal direction
       !------------------------------------------------------------------
       CALL NormalContactSet()

       ! i) If requested ensure a minimum number of contact nodes
       !-------------------------------------------------------------------
       LimitedMin = ListGetInteger( BC,'Contact Active Set Minimum',Found)
       IF( Found ) CALL IncreaseContactSet( LimitedMin )

       ! j) Determine the stick set in tangent direction
       !------------------------------------------------------------------
       CALL TangentContactSet()
       
       ! k) Add the stick coefficient if present
       !------------------------------------------------------------------
       IF( AddDiag ) THEN
         CALL StickCoefficientSet()
       END IF

       ! l) We can map information from slave to master either by creating a dual projector
       !    or using the transpose of the original projector to map field from slave to master.
       !-----------------------------------------------------------------------------------
       IF(.NOT. CreateDual ) THEN
         CALL ProjectFromSlaveToMaster()
       END IF

       ! m) If we have dynamic friction then add it 
       IF( .NOT. SkipFriction .AND. ( SlipContact .OR. FrictionContact ) ) THEN
         CALL SetSlideFriction()
       END IF

       IF( ConservativeAdd .OR. ConservativeRemove ) THEN
         DEALLOCATE( InterfaceDof )
       END IF
     END DO
     
     ! Use N-T coordinate system for the initial guess
     ! This is mandatory if using the residual mode linear solvers 
     IF( AnyRotatedContact ) THEN
       DEALLOCATE( RotatedField ) 
     END IF
     

     FirstTime = .FALSE.
     CALL Info(Caller,'All done',Level=10)

   CONTAINS


     ! Given the cartesian solution compute the rotated solution.
     !-------------------------------------------------------------------------
     SUBROUTINE RotatedDisplacementField( ) 

       REAL(KIND=dp) :: RotVec(3)
       INTEGER :: i,j,k,n,m

       IF( .NOT. AnyRotatedContact ) RETURN

       CALL Info(Caller,'Rotating displacement field',Level=8)
       ALLOCATE( RotatedField(Solver % Matrix % NumberOfRows ) )
       RotatedField = Var % Values

       n = SIZE( FieldPerm ) 
       m = SIZE( NT % BoundaryReorder )
       IF( n > m ) THEN
         i = COUNT(FieldPerm(m+1:n) > 0 )
         IF( i > 0 ) THEN
           CALL Fatal(Caller,'Number of potential untreated rotations: '//I2S(i))
         END IF
       END IF
       
       DO i=1,SIZE(FieldPerm)
         j = FieldPerm(i)
         IF( j == 0 ) CYCLE
         m = NT % BoundaryReorder(i)
         IF( m == 0 ) CYCLE
         
         RotVec = 0._dp
         DO k=1,cdofs
           RotVec(k) = RotatedField(dofs*(j-1)+k)
         END DO
         CALL RotateNTSystem( RotVec, i )
         DO k=1,cdofs
           RotatedField(dofs*(j-1)+k) = RotVec( k )
         END DO
       END DO

       IF( InfoActive(30) ) THEN
         CALL VectorValuesRange(RotatedField, SIZE(RotatedField),'RotatedField')
       END IF
       
     END SUBROUTINE RotatedDisplacementField


     ! Given the previous solution and the current stiffness matrix 
     ! computes the load normal to the surface i.e. the contact load.
     ! If we have normal-tangential coordinate system then also the load is in 
     ! the same coordinate system. 
     !-------------------------------------------------------------------------
     FUNCTION CalculateContactLoad( ) RESULT ( LoadVar )

       TYPE(Variable_t), POINTER :: LoadVar
       REAL(KIND=dp), POINTER :: TempX(:)
       REAL(KIND=dp) :: RotVec(3)
       INTEGER :: i,j,k,m


       CALL Info(Caller,'Determining reaction forces for contact problem',Level=10)

       LoadVar => VariableGet( Model % Variables, &
           TRIM(VarName) // ' Contact Load',ThisOnly = .TRUE. )
       IF( .NOT. ASSOCIATED( LoadVar ) ) THEN
         CALL Fatal(Caller, &
             'No Loads associated with variable: '//GetVarName(Var) )
       END IF

       IF( AnyRotatedContact ) THEN
         TempX => RotatedField 
       ELSE
         TempX => FieldValues
       END IF

       CALL CalculateLoads( Solver, Solver % Matrix, TempX, dofs, .FALSE., LoadVar ) 

       IF( InfoActive(30) ) THEN
         CALL VectorValuesRange(LoadVar % Values, SIZE(LoadVar % Values),'ContactLoad')
       END IF
       
     END FUNCTION CalculateContactLoad

     
     ! Given the previous solution and the related Lagrange multiplier pick the
     ! new multiplier such that it may be visualized as a field.
     !-------------------------------------------------------------------------
     SUBROUTINE PickLagrangeMultiplier( ) 

       TYPE(Variable_t), POINTER :: LinSysVar, ContactSysVar, ActiveVar
       INTEGER :: i,j,k,l,n
       INTEGER, POINTER :: InvPerm(:)
       LOGICAL :: Stat
       
       CALL Info(Caller,'Pick lagrange coefficient from the active set to whole set',Level=10)

       LinSysVar => VariableGet( Model % Variables, &
           LagrangeMultiplierName(Solver),ThisOnly = .TRUE. )
       IF( .NOT. ASSOCIATED( LinSysVar ) ) THEN
         ! Having nothing to pick from is the normal case, not a fault: the
         ! linear system multiplier is created in SolveWithLinearRestriction
         ! only when "Export Lagrange Multiplier" is given, and even then not
         ! before the first solve -- whereas this runs from
         ! DefaultFinishBoundaryAssembly, ahead of it. The exception is when
         ! the contact condition is itself asked to be based on the
         ! multiplier, as it would then be determined from a field that is
         ! never filled, i.e. from zeros.
         ! Note the test is on the keywords, not on this call: with both given
         ! the field simply does not exist yet on the first round, and picking
         ! from it starts working on the next one.
         IF( ListGetLogical( Params,'Use Lagrange Multiplier for Contact', Stat ) .AND. &
             .NOT. ListGetLogical( Params,'Export Lagrange Multiplier', Stat ) ) THEN
           CALL Fatal(Caller,'"Use Lagrange Multiplier for Contact" requires also '// &
               '"Export Lagrange Multiplier": '//GetVarName(Var) )
         ELSE
           CALL Info(Caller,'No Lagrange multiplier field to pick from, skipping: '// &
               GetVarName(Var), Level=8 )
         END IF
         RETURN
       END IF
       
       ContactSysVar => VariableGet( Model % Variables, &
           TRIM(VarName) // ' Lagrange Multiplier',ThisOnly = .TRUE. )
       IF( .NOT. ASSOCIATED( ContactSysVar ) ) THEN
         CALL Fatal(Caller, &
             'No Lagrange multiplier field associated with: '//GetVarName(Var) )
       END IF
       ContactSysVar % Values = 0.0_dp

       IF(.NOT. ASSOCIATED( Solver % Matrix % ConstraintMatrix ) ) THEN
          CALL Fatal(Caller, &
             'No constraint matrix associated with: '//GetVarName(Var) )          
       END IF            
       
       InvPerm => Solver % Matrix % ConstraintMatrix % InvPerm
       n = Solver % Matrix % ConstraintMatrix % NumberOfRows
       
       DO i=1,SIZE(InvPerm)
         ! This is related to the full matrix equation
         j = InvPerm(i)

         IF( MODULO(j,dofs) /= 1 ) CYCLE
         l = (j-1)/dofs+1
         
         IF( l > 0 .AND. l <= SIZE( ContactSysVar % Perm ) ) THEN
           k = ContactSysVar % Perm(l)
           IF( k > 0 ) THEN
             ContactSysVar % Values(k) = LinSysVar % Values(i)
           END IF
         END IF
       END DO

       IF(InfoActive(30)) THEN
         CALL VectorValuesRange( LinsysVar % Values,SIZE(LinsysVar % Values),'LinsysValues')
         CALL VectorValuesRange( ContactSysVar % Values,SIZE(ContactSysVar % Values),'ContactSysValues')
       END IF
                     
     END SUBROUTINE PickLagrangeMultiplier

     
     ! Create fields where the contact information will be saved.
     ! Create the fields both for slave and master nodes at each 
     ! contact pair.
     !--------------------------------------------------------------
     SUBROUTINE GetContactFields( DoAllocate )

       LOGICAL :: DoAllocate
       INTEGER, POINTER :: BoundaryPerm(:), Indexes(:)
       INTEGER :: i,j,k,t,n
       TYPE(Element_t), POINTER :: Element
       LOGICAL, ALLOCATABLE :: ActiveBCs(:)
       

       IF( DoAllocate ) THEN
         CALL Info(Caller,'Creating contact fields',Level=8)

         n = SIZE( FieldPerm ) 
         ALLOCATE( BoundaryPerm(n) )
         BoundaryPerm = 0
         
         ALLOCATE( ActiveBCs(Model % NumberOfBcs ) )
         ActiveBCs = .FALSE.

         DO i=1,Model % NumberOfBCs 
           j = ListGetInteger( Model % BCs(i) % Values,'Mortar BC',Found ) 
           IF(.NOT. Found ) THEN
             j = ListGetInteger( Model % BCs(i) % Values,'Contact BC',Found ) 
           END IF
           IF( j > 0 ) THEN
             ActiveBCs(i) = .TRUE.
             ActiveBCs(j) = .TRUE. 
           END IF
         END DO

         DO t=Mesh % NumberOfBulkElements + 1, &
             Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

           Element => Mesh % Elements( t )                 
           DO i = 1, Model % NumberOfBCs
             IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN
               IF( ActiveBCs(i) ) THEN
                 IF( pContact ) THEN
                   n = mGetElementDOFs(pIndexes,Element)                   
                   BoundaryPerm(pIndexes(1:n)) = 1
                 ELSE                 
                   BoundaryPerm( Element % NodeIndexes ) = 1
                 END IF
               END IF
             END IF
           END DO
         END DO

         DEALLOCATE( ActiveBCs )

         j = 0
         DO i=1,SIZE(BoundaryPerm)
           IF( BoundaryPerm(i) > 0 ) THEN
             j = j + 1
             BoundaryPerm(i) = j
           END IF
         END DO

         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Distance',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Gap',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Normalload',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Slipload',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Weight',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Active',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Stick',1,Perm = BoundaryPerm )
         IF( CalculateVelocity ) THEN
           CALL VariableAddVector( Model % Variables,Mesh,Solver,&
               TRIM(VarName)//' Contact Velocity',cDofs,Perm = BoundaryPerm )
         END IF
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Lagrange Multiplier',1,Perm = BoundaryPerm )
       END IF

       DistVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Distance')
       GapVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Gap')
       NormalLoadVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Normalload')
       SlipLoadVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Slipload')
       WeightVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Weight')
       NormalActiveVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Active')
       StickActiveVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Stick') 
       IF( CalculateVelocity ) THEN
         VeloVar => VariableGet( Model % Variables,&
             TRIM(VarName)//' Contact Velocity')
       END IF
       
       ContactLagrangeVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Lagrange Multiplier')       
       
       NormalActiveVar % Values = -1.0_dp
       StickActiveVar % Values = -1.0_dp

     END SUBROUTINE GetContactFields



     ! Allocates the vectors related to the mortar contact surface, if needed.
     ! Initialize the mortar vectors and mortar permutation future use. 
     ! As the geometry changes the size of the projectors may also change. 
     !----------------------------------------------------------------------------
     SUBROUTINE InitializeMortarVectors()

       INTEGER :: onesize, totsize
       INTEGER, POINTER :: Perm(:)
       LOGICAL, POINTER :: Active(:)
       REAL(KIND=dp), POINTER :: Diag(:)
       LOGICAL :: SamePerm, SameSize
       
       onesize = Projector % NumberOfRows
       totsize = cDofs * onesize

       IF( .NOT. AddDiag .AND. ASSOCIATED(MortarBC % Diag) ) THEN
         DEALLOCATE( MortarBC % Diag ) 
       END IF

       ! Create the permutation that is later needed in putting the diag and rhs to correct position.
       ALLOCATE( Perm( SIZE( FieldPerm ) ) )
       Perm = 0
       DO i=1,SIZE( Projector % InvPerm )
         j = Projector % InvPerm(i) 
         IF( j == 0 ) CYCLE
         Perm( j ) = i
       END DO

       ! First time nothing is allocated
       IF( .NOT. ASSOCIATED( MortarBC % Perm ) ) THEN
         CALL Info(Caller,'Allocating projector mortar vectors of size: '//I2S(totsize),Level=10)
         ALLOCATE( MortarBC % Active( totsize ), MortarBC % Rhs( totsize) )
         MortarBC % Active = .FALSE.
         MortarBC % Rhs = 0.0_dp
         MortarBC % Perm => Perm 

         IF( AddDiag ) THEN
           ALLOCATE( MortarBC % Diag( totsize ) )
           MortarBC % Diag = 0.0_dp
         END IF

         RETURN
       END IF

       
       ! If permutation has changed we need to change the vectors also
       SamePerm = ALL( Perm == MortarBC % Perm )

       ! Permutation unchanged, just return
       IF( SamePerm ) THEN
         DEALLOCATE( Perm ) 
         RETURN
       END IF
       
       ! Permutation changes, and also sizes changed?
       SameSize = ( SIZE(MortarBC % Rhs) == totsize )
       IF(.NOT. SameSize ) THEN
         DEALLOCATE( MortarBC % Rhs )
         ALLOCATE( MortarBC % Rhs( totsize ) )
         MortarBC % Rhs = 0.0_dp
       END IF

       ! .NOT. SamePerm
       ALLOCATE(Active(totsize))
       Active = .FALSE.

       IF( AddDiag ) THEN
         ALLOCATE( Diag(totsize) )
         Diag = 0.0_dp
       END IF


       DO i=1,SIZE( Perm ) 
         j = Perm(i)
         IF( j == 0 ) CYCLE

         k = MortarBC % Perm(i)
         IF( k == 0 ) CYCLE

         DO l=1,cDofs
           Active(cDofs*(j-1)+l) = MortarBC % Active(cDofs*(k-1)+l)
         END DO
       END DO

       IF( ASSOCIATED( MortarBC % Active ) ) DEALLOCATE( MortarBC % Active )
       IF( ASSOCIATED( MortarBC % Perm ) ) DEALLOCATE( MortarBC % Perm )
       MortarBC % Active => Active 
       MortarBC % Perm => Perm 

       IF( AddDiag ) THEN
         IF( ASSOCIATED( MortarBC % Diag ) ) DEALLOCATE( MortarBC % Diag ) 
         MortarBC % Diag => Diag 
       END IF

       CALL Info(Caller,'Copied > Active < flag to changed projector',Level=8)

     END SUBROUTINE InitializeMortarVectors

     

     ! Make a list of interface dofs to allow conservative algorithms. 
     ! There only nodes that are at the interface are added or removed from the set.
     !------------------------------------------------------------------------------
     SUBROUTINE MarkInterfaceDofs()
       
       INTEGER :: i,j,i2,j2,k,k2,l,n,ind,ind2,elem
       INTEGER, POINTER :: Indexes(:)
       TYPE(Element_t), POINTER :: Element
       
       CALL Info(Caller,'Marking interface dofs for conservative adding/removal',Level=8)

       IF(.NOT. ALLOCATED( InterfaceDof ) ) THEN
         ALLOCATE( InterfaceDof( SIZE(MortarBC % Active) ) )
       END IF
       InterfaceDof = .FALSE. 


       DO elem=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
         
         Element => Mesh % Elements( elem )         
         IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc_ind) % Tag ) CYCLE

         IF( pContact ) THEN
           n = mGetElementDOFs(pIndexes,Element)                   
           Indexes => pIndexes
         ELSE         
           n = Element % TYPE % NumberOfNodes         
           Indexes => Element % NodeIndexes
         END IF
           
         DO i=1,n
           j = FieldPerm( Indexes(i) )
           IF( j == 0 ) CYCLE
           k = MortarBC % Perm( Indexes(i) )
           
           DO i2 = i+1,n
             j2 = FieldPerm( Indexes(i2) )
             IF( j2 == 0 ) CYCLE
             k2 = MortarBC % perm( Indexes(i2) )
             
             DO l=1,cDofs             
               ind = cDofs * ( k - 1 ) + l
               ind2 = cDofs * ( k2 - 1) + l
               
               IF( MortarBC % Active(ind) .NEQV. MortarBC % Active(ind2) ) THEN
                 InterfaceDof(ind) = .TRUE.
                 InterfaceDof(ind2) = .TRUE.
               END IF
             END DO
           END DO
         END DO
       END DO

       n = COUNT(InterfaceDof)
       CALL Info(Caller,'Number of interface dofs: '//I2S(n),Level=8)
       
     END SUBROUTINE MarkInterfaceDofs
     

     ! Calculates the signed distance that is used to define whether we have contact or not.
     ! If distance is negative then we can later add the corresponding node to the contact set
     ! Also computes the right-hand-side of the mortar equality constrained which is the 
     ! desired distance in the active direction. Works also for residual mode which greatly 
     ! improves the convergence for large displacements.  
     !----------------------------------------------------------------------------------------
     SUBROUTINE CalculateMortarDistance()

       REAL(KIND=dp) :: Disp(3), Coord(3), PrevDisp(3), Velo(3), ContactVec(3), ContactVelo(3), &
           LocalNormal0(3), SlipCoord(3), CartVec(3), ContactDist
       REAL(KIND=dp), POINTER :: DispVals(:), PrevDispVals(:) 
       REAL(KIND=dp) :: MinDist, MaxDist, wsum, wsumM, mult
       TYPE(Matrix_t), POINTER :: ActiveProjector
       LOGICAL :: IsSlave, IsMaster, DistanceSet
       LOGICAL, ALLOCATABLE :: SlaveNode(:), MasterNode(:), NodeDone(:)
       INTEGER, POINTER :: Indexes(:)
       INTEGER :: elemcode, CoeffSign
       REAL(KIND=dp), ALLOCATABLE :: CoeffTable(:)
       INTEGER :: l2,elem,i1,i2,j1,j2,n
       LOGICAL :: LinearContactGap, DebugNormals
       
       
       CALL Info('CalculateMortarDistance','Computing distance between mortar boundaries',Level=14)

       DebugNormals = .TRUE.
       
       DispVals => Solver % Variable % Values
       IF( .NOT. ASSOCIATED( DispVals ) ) THEN
         CALL Fatal('CalculateMortarDistance','Displacement variable not associated!')
       END IF

       IF( CalculateVelocity ) THEN
         IF( .NOT. ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
           CALL Fatal('CalculateMortarDistance','Displacement PrevValues not associated!')         
         END IF
         IF( Solver % TimeOrder == 1 ) THEN
           PrevDispVals => Solver % Variable % PrevValues(:,1)
         ELSE
           PrevDispVals => Solver % Variable % PrevValues(:,3)
         END IF
         IF(.NOT. ASSOCIATED( PrevDispVals ) ) CALL Fatal('CalculateMortarDistance',&
             'Previous displacement field required!')
       END IF

       LinearContactGap = ListGetLogical( Model % Simulation,&
           'Contact BCs linear gap', Found )      

       ALLOCATE( SlaveNode( SIZE( FieldPerm ) ) )
       SlaveNode = .FALSE.

       IF( CreateDual ) THEN
         ALLOCATE( MasterNode( SIZE( FieldPerm ) ) )
         MasterNode = .FALSE.
       END IF

       DO i=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

         Element => Mesh % Elements( i )                  
         IF( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag ) THEN
           IF( pContact ) THEN
             n = mGetElementDOFs(pIndexes,Element)                   
             SlaveNode(pIndexes(1:n)) = .TRUE.
           ELSE
             SlaveNode( Element % NodeIndexes ) = .TRUE.
           END IF
         END IF
         IF( CreateDual ) THEN
           IF ( Element % BoundaryInfo % Constraint == Model % BCs(master_ind) % Tag ) THEN
             IF( pContact ) THEN
               n = mGetElementDOFs(pIndexes,Element)                   
               MasterNode(pIndexes(1:n)) = .TRUE.
             ELSE               
               MasterNode( Element % NodeIndexes ) = .TRUE.
             END IF
           END IF
         END IF
       END DO

       ! First create the master, then the slave if needed
       IsSlave = .TRUE.
       IsMaster = .NOT. IsSlave
       ActiveProjector => Projector
       MinDist = HUGE(MinDist)
       MaxDist = -HUGE(MaxDist)

       IF( .NOT. ASSOCIATED( ActiveProjector ) ) THEN
         CALL Fatal('CalculateMortarDistance','Projector not associated!')
       END IF

       DebugNormals = ListGetLogical( Params,'Debug Normals',Found ) 

       IF( DebugNormals ) THEN
         PRINT *,'Flags:',TieContact,ResidualMode,ThisRotatedContact,NodalNormal,StickContact,RotationalProjector
       END IF

       
100    CONTINUE

       DO i = 1,ActiveProjector % NumberOfRows

         j = ActiveProjector % InvPerm(i)

         IF( j == 0 ) CYCLE
         
         wsum = 0.0_dp
         wsumM = 0.0_dp
         Dist = 0.0_dp
         DistN = 0.0_dp
         DistT1 = 0.0_dp
         DistT2 = 0.0_dp
         ContactVelo = 0.0_dp
         ContactVec = 0.0_dp
         DistanceSet = .FALSE.
         ContactDist = 0.0_dp
         CartVec = 0.0_dp
         
         ! This is the most simple contact condition. We just want no slip on the contact.
         IF( TieContact .AND. .NOT. ResidualMode ) GOTO 200

         ! Get the normal of the slave surface.
         IF( ThisRotatedContact ) THEN
           Rotated = GetSolutionRotation(NTT, j )
           LocalNormal = NTT(:,1)
           LocalNormal0 = LocalNormal
           LocalT1 = NTT(:,2)
           IF( cDofs == 3 ) LocalT2 = NTT(:,3)
         ELSE
           LocalNormal = ContactNormal
           LocalT1 = ContactT1
           IF( cDofs == 3 ) LocalT2 = ContactT2 
         END IF

         ! Compute normal of the master surface from the average sum of normals
         IF( NodalNormal ) THEN
           LocalNormal = 0.0_dp
           LocalT1 = 0.0_dp
           LocalT2 = 0.0_dp

           DO j = ActiveProjector % Rows(i),ActiveProjector % Rows(i+1)-1
             k = ActiveProjector % Cols(j)
             IF(DebugNormals) THEN
               IF(k<1 .OR. k>SIZE(FieldPerm)) THEN
                 PRINT *,'k is beyond FieldPerm:',i,j,k,SIZE(FieldPerm)
                 CYCLE
               END IF
             END IF
               
             l = FieldPerm( k ) 
             IF( l == 0 ) CYCLE

             coeff = ActiveProjector % Values(j)             
             Rotated = GetSolutionRotation(NTT, k )

             ! Weighted direction for the unit vectors
             LocalNormal = LocalNormal + coeff * NTT(:,1)
             LocalT1 = LocalT1 + coeff * NTT(:,2)
             IF( cDofs == 3 ) LocalT2 = LocalT2 + coeff * NTT(:,3)
           END DO

           ! Normalize the unit vector length to one
           LocalNormal = LocalNormal / SQRT( SUM( LocalNormal**2 ) )
           LocalT1 = LocalT1 / SQRT( SUM( LocalT1**2 ) )
           IF( cDofs == 3 ) LocalT2 = LocalT2 / SQRT( SUM( LocalT1**2 ) )

           !PRINT *,'NodalNormal:',i,j,LocalNormal0,LocalNormal
         END IF

         ! For debugging reason, check that normals are roughly opposite
         IF( DebugNormals ) THEN
           DO j = ActiveProjector % Rows(i),ActiveProjector % Rows(i+1)-1
             k = ActiveProjector % Cols(j)

             IF(DebugNormals) THEN
               IF(k<1 .OR. k>SIZE(FieldPerm)) THEN
                 PRINT *,'k2 is beyond FieldPerm:',i,j,k,SIZE(FieldPerm)
                 CYCLE
               END IF
             END IF
             
             l = FieldPerm( k ) 
             IF( l == 0 ) CYCLE

             Rotated = GetSolutionRotation(NTT, k )
             coeff = SUM( LocalNormal * NTT(:,1) ) + SUM( LocalT1*NTT(:,2)) + SUM(LocalT2*NTT(:,3))
             IF( SlaveNode(k) .AND. coeff < 2.5_dp ) THEN
               Found = .TRUE.
               !PRINT *,'Slave Normal:',i,j,k,Rotated,coeff
             ELSE IF( .NOT. SlaveNode(k) .AND. coeff > -2.5_dp ) THEN
               Found = .TRUE.
               !PRINT *,'Master Normal:',i,j,k,Rotated,coeff
             ELSE
               Found = .FALSE.
             END IF
             IF( Found ) THEN
               !PRINT *,'Prod:',SUM( LocalNormal * NTT(:,1) ), SUM( LocalT1*NTT(:,2)), SUM(LocalT2*NTT(:,3))
               !PRINT *,'N:',LocalNormal,NTT(:,1)
               !PRINT *,'T1:',LocalT1,NTT(:,2)
               !PRINT *,'T2:',LocalT2,NTT(:,3)
             END IF
           END DO
         END IF

         ! Compute the weigted distance in the normal direction.
         DO j = ActiveProjector % Rows(i),ActiveProjector % Rows(i+1)-1
           k = ActiveProjector % Cols(j)

           IF(DebugNormals) THEN
             IF(k<1 .OR. k>SIZE(FieldPerm)) THEN
               PRINT *,'k3 is beyond FieldPerm:',i,j,k,SIZE(FieldPerm)
               CYCLE
             END IF
           END IF
           
           l = FieldPerm( k ) 
           IF( l == 0 ) CYCLE

           coeff = ActiveProjector % Values(j)
                  
           ! Only compute the sum related to the active projector
           IF( SlaveNode(k) ) THEN
             wsum = wsum + coeff
           ELSE
             wsumM = wsumM + coeff
           END IF           
         END DO

         IF( ABS( wsum ) <= TINY( wsum ) ) THEN
           CALL Fatal('CalculateMortarDistance','wsum seems to be almost zero!')
         END IF
         IF( ABS( wsumM ) <= TINY( wsumM ) ) THEN
           CALL Fatal('CalculateMortarDistance','wsumM seems to be almost zero!')
         END IF

         ! Slave and master multipliers should sum up to same value
         mult = ABS( wsum / wsumM ) 
         
         ! Compute the weigted distance in the normal direction.
         DO j = ActiveProjector % Rows(i),ActiveProjector % Rows(i+1)-1
           k = ActiveProjector % Cols(j)

           IF(DebugNormals) THEN
             IF(k<1 .OR. k>SIZE(FieldPerm)) THEN
               PRINT *,'k4 is beyond FieldPerm:',i,j,k,SIZE(FieldPerm)
               CYCLE
             END IF
           END IF
           
           l = FieldPerm( k ) 
           IF( l == 0 ) CYCLE

           ! This includes only the coordinate since the displacement
           ! is added to the coordinate!
           coeff = ActiveProjector % Values(j)

           CoeffSign = 1
           
           ! Only compute the sum related to the active projector
           IF( .NOT. SlaveNode(k) ) THEN
             coeff = mult * coeff
             IF( ThisRotatedContact ) CoeffSign = -1
           END IF
           
           disp(1) = DispVals( dofs * (l-1) + 1)
           disp(2) = DispVals( dofs * (l-1) + 2 )
           IF( cdofs == 2 ) THEN
             disp(3) = 0.0_dp
           ELSE
             disp(3) = DispVals( dofs * (l-1) + 3 )
           END IF

           ! If nonlinear analysis is used we may need to cancel the introduced gap due to numerical errors 
           IF( TieContact .AND. ResidualMode ) THEN !.AND. k <= dofs * Mesh % NumberOfNodes ) THEN
             IF( ThisRotatedContact ) THEN
               ContactVec(1) = ContactVec(1) + coeff * SUM( LocalNormal * Disp )
               ContactVec(2) = ContactVec(2) + coeff * SUM( LocalT1 * Disp )
               IF( cDofs == 3) ContactVec(3) = ContactVec(3) + coeff * SUM( LocalT2 * Disp )
             ELSE
               ContactVec(1) = ContactVec(1) + coeff * SUM( ContactNormal * Disp )
               ContactVec(2) = ContactVec(2) + coeff * SUM( ContactT1 * Disp )
               IF( cDofs == 3 ) ContactVec(3) = ContactVec(3) + coeff * SUM( ContactT2 * Disp ) 
             END IF
             CYCLE
           END IF

           coord(1) = Mesh % Nodes % x( k ) 
           coord(2) = Mesh % Nodes % y( k ) 
           coord(3) = Mesh % Nodes % z( k ) 

           PrevDisp = 0._dp
           IF( CalculateVelocity ) THEN
             PrevDisp(1) = PrevDispVals( dofs * (l-1) + 1)
             PrevDisp(2) = PrevDispVals( dofs * (l-1) + 2 )
             IF( cdofs == 2 ) THEN
               PrevDisp(3) = 0.0_dp
             ELSE
               PrevDisp(3) = PrevDispVals( dofs * (l-1) + 3 )
             END IF
           END IF

           ! If the linear system is in residual mode also set the current coordinate in residual mode too!
           ! Note that displacement field is given always in cartesian coordinates!
           IF( ResidualMode ) THEN
             Coord = Coord + Disp
           END IF

           ! DistN is used to give the distance that we need to move the original coordinates
           ! in the wanted direction in order to have contact.
           IF( ThisRotatedContact ) THEN
             ContactVec(1) = ContactVec(1) + coeff * SUM( LocalNormal * Coord )
           ELSE
             ContactVec(1) = ContactVec(1) + coeff * SUM( ContactNormal * Coord )
           END IF

           ! Tangential distances needed to move the original coordinates to the contact position
           ! If stick is required then we want to keep the tangential slip zero. 
           IF( StickContact ) THEN             
             SlipCoord = -PrevDisp 
             IF( ResidualMode ) SlipCoord = SlipCoord + Disp 

             IF( ThisRotatedContact ) THEN
               ContactVec(2) = ContactVec(2) + coeff * SUM( LocalT1 * SlipCoord )
               IF( cDofs == 3) ContactVec(3) = ContactVec(3) + coeff * SUM( LocalT2 * SlipCoord )
             ELSE
               ContactVec(2) = ContactVec(2) + coeff * SUM( ContactT1 * SlipCoord )
               IF( cDofs == 3 ) ContactVec(3) = ContactVec(3) + coeff * SUM( ContactT2 * SlipCoord )
             END IF
           END IF

           ! If not in the residual mode still take into account the displacement for the condition
           IF( .NOT. ResidualMode ) Coord = Coord + Disp

           ! Dist is used to compute the current signed distance that is used to determine
           ! whether we have contact or not. 
           IF( RotationalProjector ) THEN
             Dist = Dist + coeff * SQRT( SUM( Coord**2 ) )
           ELSE IF( NormalProjector ) THEN
             Dist = Dist + coeff * SUM( LocalNormal * Coord )
           ELSE             
             Dist = Dist + coeff * SUM( ContactNormal * Coord )
           END IF

           CartVec = CartVec + coeff * Coord
           
           IF( CalculateVelocity ) THEN
             Velo = ( Disp - PrevDisp ) !/ dt
             ContactVelo(1) = ContactVelo(1) + coeff * SUM( Velo * LocalNormal ) 
             ContactVelo(2) = ContactVelo(2) + coeff * SUM( Velo * LocalT1 )
             ContactVelo(3) = ContactVelo(3) + coeff * SUM( Velo * LocalT2 ) 
           END IF
           DistanceSet = .TRUE.
         END DO

         ! Divide by weight to get back to real distance in the direction of the normal
         ContactVec = ContactVec / wsum 
         Dist = DistSign * Dist / wsum
         IF( CalculateVelocity ) THEN
           ContactVelo = ContactVelo / wsum
         END IF
         CartVec = CartVec / wsum
         
200      IF( IsSlave ) THEN

           MortarBC % Rhs(cDofs*(i-1)+DofN) = -ContactVec(1)
           IF( StickContact .OR. TieContact ) THEN
             MortarBC % Rhs(cDofs*(i-1)+DofT1) = -ContactVec(2) 
             IF( cDofs == 3 ) THEN
               MortarBC % Rhs(cDofs*(i-1)+DofT2) = -ContactVec(3)
             END IF
           END IF
           
           MinDist = MIN( Dist, MinDist ) 
           MaxDist = MAX( Dist, MaxDist )
         END IF

         IF( IsMaster ) THEN
           Dist = -Dist
           ContactVelo = -ContactVelo
         END IF

         ! We use the same permutation for all boundary variables
         IF(ActiveProjector % InvPerm(i) <= 0 ) CYCLE
         j = DistVar % Perm( ActiveProjector % InvPerm(i) )

         DistVar % Values( j ) = Dist

         GapVar % Values( j ) = ContactVec(1)

         IF( CalculateVelocity ) THEN
           DO k=1,cDofs             
             VeloVar % Values( cDofs*(j-1)+k ) = ContactVelo(k) 
           END DO
         END IF
       END DO
       
       IF( IsSlave ) THEN
         IF( CreateDual ) THEN
           IsSlave = .FALSE.
           IsMaster = .NOT. IsSlave
           ActiveProjector => DualProjector
           GOTO 100
         END IF
       END IF

       
       IF( LinearContactGap .OR. pContact ) THEN       
         DO elem=Mesh % NumberOfBulkElements + 1, &
             Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
           
           Element => Mesh % Elements( elem )         
           
           IsSlave = ( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag ) 
           IsMaster = ( Element % BoundaryInfo % Constraint == Model % BCs(master_ind) % Tag ) 
           IF( .NOT. ( IsSlave .OR. ( CreateDual .AND. IsMaster ) ) ) CYCLE
           
           ElemCode = Element % TYPE % ElementCode           
           IF( pContact ) THEN
             n = mGetElementDOFs(pIndexes,Element)                   
             Indexes => pIndexes
           ELSE            
             n = Element % TYPE % NumberOfNodes
             Indexes => Element % NodeIndexes         
           END IF
           
           SELECT CASE ( ElemCode )
           CASE( 202, 203 )
             i=3
             
             j = DistVar % Perm(Indexes(i))
             IF( j > 0 ) THEN
               IF( pContact ) THEN
                 DistVar % Values(j) = 0.0_dp
                 GapVar % Values(j) = 0.0_dp
                 IF( CalculateVelocity ) THEN
                   DO k=1,cDofs
                     VeloVar % Values( cDofs*(j-1)+k ) = 0.0_dp
                   END DO
                 END IF
               ELSE
                 i1=1
                 i2=2
                 j1 = DistVar % Perm(Indexes(i1))
                 j2 = DistVar % Perm(Indexes(i2))
                 
                 DistVar % Values(j) = 0.5_dp * &
                     ( DistVar % Values(j1) + DistVar % Values(j2))
                 GapVar % Values(j) = 0.5_dp * &
                     ( GapVar % Values(j1) + GapVar % Values(j2))
                 
                 IF( CalculateVelocity ) THEN
                   DO k=1,cDofs
                     VeloVar % Values( cDofs*(j-1)+k ) = 0.5_dp * &
                         ( VeloVar % Values(cDofs*(j1-1)+k) + VeloVar % Values(cDofs*(j2-1)+k))  
                   END DO
                 END IF
               END IF
             END IF

           CASE( 404, 408 )
             DO i=5,8
               
               j = DistVar % Perm(Indexes(i))
               IF( j == 0 ) CYCLE
               
               IF( pContact ) THEN
                 DistVar % Values(j) = 0.0_dp
                 GapVar % Values(j) = 0.0_dp
                 IF( CalculateVelocity ) THEN
                   DO k=1,cDofs
                     VeloVar % Values( cDofs*(j-1)+k ) = 0.0_dp
                   END DO
                 END IF
               ELSE
                 i1=i-4
                 i2=i1+1
                 IF(i2==5) i2=1
                 j1 = DistVar % Perm(Indexes(i1))
                 j2 = DistVar % Perm(Indexes(i2))
                 
                 DistVar % Values(j) = 0.5_dp * &
                     ( DistVar % Values(j1) + DistVar % Values(j2))
                 GapVar % Values(j) = 0.5_dp * &
                     ( GapVar % Values(j1) + GapVar % Values(j2))
                 
                 IF( CalculateVelocity ) THEN
                   DO k=1,cDofs
                     VeloVar % Values( cDofs*(j-1)+k ) = 0.5_dp * &
                         ( VeloVar % Values(cDofs*(j1-1)+k) + VeloVar % Values(cDofs*(j2-1)+k))  
                   END DO
                 END IF
               END IF
             END DO
               
           CASE DEFAULT
             CALL Fatal('CalculateMortarDistance','Implement linear gaps for: '//I2S(ElemCode))
           END SELECT

         END DO
       END IF
       
       DEALLOCATE( SlaveNode )
       IF( CreateDual ) DEALLOCATE( MasterNode )
       
       IF( InfoActive(30) ) THEN
         ! We don't know if other partitions are here, so let us not make parallel reductions!
         CALL VectorValuesRange(DistVar % Values,SIZE(DistVar % Values),'Dist',.TRUE.)
         CALL VectorValuesRange(GapVar % Values,SIZE(GapVar % Values),'Gap',.TRUE.)
         CALL VectorValuesRange(MortarBC % rhs,SIZE(MortarBC % rhs),'Mortar Rhs',.TRUE.)       
       END IF

       CALL Info(Caller,'Finished computing mortar distance',Level=25)
       
     END SUBROUTINE CalculateMortarDistance



     ! Calculates the contact pressure in the normal direction from the nodal loads.
     ! The nodal loads may be given either in cartesian or n-t coordinate system. 
     !-------------------------------------------------------------------------------
     SUBROUTINE CalculateContactPressure()
       
       INTEGER :: elem
       INTEGER, POINTER :: Indexes(:)
       TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
       TYPE(Nodes_t) :: Nodes
       INTEGER :: i,j,k,t,CoordSys, NormalSign0, NormalSign, NormalCount
       REAL(KIND=dp) :: s, x, DetJ, u, v, w, Normal(3),NodalForce(3),DotProd, &
           NormalForce, SlipForce
       REAL(KIND=dp), ALLOCATABLE :: Basis(:)
       LOGICAL :: Stat, IsSlave, IsMaster
       TYPE(Matrix_t), POINTER :: ActiveProjector
       LOGICAL, ALLOCATABLE :: NodeDone(:)
       LOGICAL :: LinearContactLoads
       INTEGER :: i1,i2,j1,j2,ElemCode,m
       
       n = Mesh % MaxElementNodes
       ALLOCATE(Basis(2*n), Nodes % x(2*n), Nodes % y(2*n), Nodes % z(2*n) )
       Nodes % x = 0.0_dp; Nodes % y = 0.0_dp; Nodes % z = 0.0_dp

       CALL Info(Caller,'Computing pressure for contact problems',Level=20)
       
       CoordSys = CurrentCoordinateSystem()
       NodalForce = 0.0_dp

       NormalSign0 = 0
       NormalCount = 0
       
       ALLOCATE( NodeDone( SIZE( FieldPerm ) ) )
       NodeDone = .FALSE.

       LinearContactLoads = ListGetLogical( Model % Simulation,&
           'Contact BCs linear loads', Found )

       
100    DO elem=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
         
         Element => Mesh % Elements( elem )         

         IsSlave = ( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag ) 
         IsMaster = ( Element % BoundaryInfo % Constraint == Model % BCs(master_ind) % Tag ) 

         IF( .NOT. ( IsSlave .OR. ( CreateDual .AND. IsMaster ) ) ) CYCLE
                  
         IF( pContact ) THEN
           n = mGetElementDOFs(pIndexes,Element)                   
           Indexes => pIndexes
         ELSE         
           n = Element % TYPE % NumberOfNodes         
           Indexes => Element % NodeIndexes
         END IF
         
         Nodes % x(1:n) = Mesh % Nodes % x(Indexes(1:n))
         Nodes % y(1:n) = Mesh % Nodes % y(Indexes(1:n))
         Nodes % z(1:n) = Mesh % Nodes % z(Indexes(1:n))
         
         IntegStuff = GaussPoints( Element )

         DO t=1,IntegStuff % n        
           U = IntegStuff % u(t)
           V = IntegStuff % v(t)
           W = IntegStuff % w(t)
           
           stat = ElementInfo( Element, Nodes, U, V, W, detJ, Basis )
           S = DetJ * IntegStuff % s(t)
           
           IF ( CoordSys /= Cartesian ) THEN
             X = SUM( Nodes % X(1:n) * Basis(1:n) )
             s = s * x
           END IF
           
           Normal = NormalVector( Element,Nodes,u,v,.TRUE. )

           ! Check the consistency of sign in the projector
           IF( IsSlave .AND. ( FlatProjector .OR. PlaneProjector .OR. NormalProjector ) ) THEN
             DotProd = SUM( Normal * ContactNormal ) 
             IF( DotProd < 0.0 ) THEN
               NormalSign = 1
             ELSE
               NormalSign = -1 
             END IF
             IF( NormalSign0 == 0 ) THEN
               NormalSign0 = NormalSign
             ELSE
               IF( NormalSign0 /= NormalSign ) NormalCount = NormalCount + 1
             END IF
           END IF

           DO i=1,n
             j = NormalLoadVar % Perm( Indexes(i) )
             IF( j == 0 ) CYCLE             
             
             IF( .NOT. NodeDone( Indexes(i) ) ) THEN             
               NodeDone( Indexes(i) ) = .TRUE.
               WeightVar % Values(j) = 0.0_dp
               NormalLoadVar % Values(j) = 0.0_dp
               SlipLoadVar % Values(j) = 0.0_dp
             END IF

             k = FieldPerm( Indexes(i) )
             IF( k == 0 ) CYCLE
             
             DO l=1,cdofs
               NodalForce(l) = LoadValues(dofs*(k-1)+l)
             END DO

             IF( ThisRotatedContact ) THEN
               NormalForce = NodalForce(1)
             ELSE
               NormalForce = SUM( NodalForce * Normal ) 
             END IF
             ! By construction the expression should be positive but due to rounding errors we have to ensure it!
             SlipForce = SQRT( MAX(0.0_dp, SUM( NodalForce**2 ) - NormalForce**2 ) )
             
             NormalLoadVar % Values(j) = NormalLoadVar % Values(j) - &
                 s * Basis(i) * NormalForce
             SlipLoadVar % Values(j) = SlipLoadVar % Values(j) + &
                 s * Basis(i) * SlipForce
             
             WeightVar % Values(j) = WeightVar % Values(j) + s * Basis(i)
           END DO
           
         END DO
       END DO
       
       ! Normalize the computed normal loads such that the unit will be that of pressure
       DO i=1,SIZE(FieldPerm)
         IF( NodeDone( i ) ) THEN             
           j = WeightVar % Perm(i)
           IF(j==0) CYCLE
           s = WeightVar % Values(j)
           IF( s /= s ) CYCLE
           IF( ABS(s) > EPSILON(s) ) THEN
             SlipLoadVar % Values(j) = SlipLoadVar % Values(j) / s**2
             NormalLoadVar % Values(j) = NormalLoadVar % Values(j) / s**2
           END IF
         END IF
       END DO

       IF( LinearContactLoads ) THEN       
         DO elem=Mesh % NumberOfBulkElements + 1, &
             Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
           
           Element => Mesh % Elements( elem )         

           IsSlave = ( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag ) 
           IsMaster = ( Element % BoundaryInfo % Constraint == Model % BCs(master_ind) % Tag ) 
           IF( .NOT. ( IsSlave .OR. ( CreateDual .AND. IsMaster ) ) ) CYCLE
           
           Indexes => Element % NodeIndexes         
           n = Element % TYPE % NumberOfNodes
           ElemCode = Element % TYPE % ElementCode           

           SELECT CASE ( ElemCode )
             
           CASE( 203, 306, 408 )
             n = MODULO( ElemCode, 100 )
             m = ElemCode / 100             

             DO i=m+1,n
               i1=i-m
               IF(i1==m) THEN                 
                 i2=m+1
               ELSE
                 i2=1
               END IF
               j = SlipLoadVar % Perm(Indexes(i))
               j1 = SlipLoadVar % Perm(Indexes(i1))
               j2 = SlipLoadVar % Perm(Indexes(i2))
               SlipLoadVar % Values(j) = 0.5_dp * &
                   ( SlipLoadVar % Values(j1) + SlipLoadVar % Values(j2))
               NormalLoadVar % Values(j) = 0.5_dp * &
                   ( NormalLoadVar % Values(j1) + NormalLoadVar % Values(j2))
             END DO
               
           CASE DEFAULT
             CALL Fatal(Caller,'Implement linear loads for: '//I2S(ElemCode))
           END SELECT
         END DO
       END IF

       IF( InfoActive(20) ) THEN
         CALL VariableValuesRange(SlipLoadVar,'SlipLoad',AlwaysSerial=.TRUE.)       
         CALL VariableValuesRange(NormalLoadVar,'NormalLoad',AlwaysSerial=.TRUE.)       
       END IF
       
       IF( FlatProjector .OR. PlaneProjector .OR. NormalProjector ) THEN
         IF( NormalCount == 0 ) THEN
           CALL Info(Caller,'All normals are consistently signed',Level=10)
         ELSE
           CALL Warn(Caller,'There are normals with conflicting signs: '&
               //I2S(NormalCount)  )
           NormalSign = 1
         END IF
         CALL Info(Caller,'Normal direction for distance measure: '&
             //I2S(NormalSign),Level=8)
         DistSign = NormalSign 
       END IF

       ! Check whether the normal sign has been enforced
       IF( ListGetLogical( BC,'Normal Sign Negative',Found ) ) DistSign = -1
       IF( ListGetLogical( BC,'Normal Sign Positive',Found ) ) DistSign = 1

       DEALLOCATE( Basis, Nodes % x, Nodes % y, Nodes % z, NodeDone )

       CALL Info(Caller,'Finished computing contact pressure',Level=25)
       
     END SUBROUTINE CalculateContactPressure

  

     ! Sets the contact in the normal direction by looking at the signed distance and 
     ! contact force. The initial contact set may be enlarged to eliminate null-space
     ! related to rigid-body motion.
     !----------------------------------------------------------------------------------
     SUBROUTINE NormalContactSet() 
       
       INTEGER :: LimitSign, Removed, Added
       REAL(KIND=dp) :: DistOffSet, MinLoad, MaxLoad, NodeLoad, MinDist, MaxDist, NodeDist
       INTEGER :: i,j,k,ind
       LOGICAL :: Found

       CALL Info('NormalContactSet','Defining normal contact set',Level=20)
       
       ! This is related to the formulation of the PDE and is probably fixed for all elasticity solvers      
       LimitSign = -1

       Removed = 0
       Added = 0        
       MinLoad = HUGE(MinLoad)
       MaxLoad = -HUGE(MaxLoad)
       MinDist = HUGE(MinDist)
       MaxDist = -HUGE(MaxDist)

       Found = .FALSE.
       IF( FirstTime ) THEN
         DistOffset = ListGetCReal( BC,&
             'Mortar BC Initial Contact Depth',Found)
         IF(.NOT. Found ) DistOffset = ListGetCReal( BC,&
             'Contact Depth Offset Initial',Found)
       END IF
       IF( .NOT. Found ) DistOffset = ListGetCReal( BC,&
           'Contact Depth Offset',Found)

       IF( InfoActive(30) ) THEN      
         PRINT *,'InitialActiveSet:',COUNT( MortarBC % Active ), SIZE( MortarBC % Active )
       END IF
       
       ! Determine now whether we have contact or not
       DO i = 1,Projector % NumberOfRows
         j = Projector % InvPerm( i ) 
         IF( j == 0 ) CYCLE
         k = FieldPerm( j ) 
         IF( k == 0 ) CYCLE
         k = UseLoadVar % Perm(j)

         ind = cDofs * (i-1) + DofN

         ! Tie contact should always be in contact - if we have found a counterpart
         IF( TieContact ) THEN
           MortarBC % Active(ind) = .TRUE.
           CYCLE
         END IF

         ! Enforce contact 
         !------------------------------------------------------
         coeff = ListGetRealAtNode( BC,'Contact Active Condition', j, Found )
         IF( Found .AND. coeff > 0.0_dp ) THEN
           MortarBC % Active(ind) = .TRUE.
           CYCLE
         END IF

         ! Enforce no contact
         !------------------------------------------------------
         coeff = ListGetRealAtNode( BC,'Contact Passive Condition', j, Found )
         IF( Found .AND. coeff > 0.0_dp ) THEN
           MortarBC % Active(ind) = .FALSE.
           CYCLE
         END IF

         ! Free nodes with wrong sign in contact force
         !--------------------------------------------------------------------------       
         IF( MortarBC % Active( ind ) ) THEN
           NodeLoad = UseLoadVar % Values(k)
           MaxLoad = MAX( MaxLoad, NodeLoad )
           MinLoad = MIN( MinLoad, NodeLoad )
           DoRemove = ( LimitSign * NodeLoad > LimitSign * LoadEps ) 
           IF( DoRemove .AND. ConservativeRemove ) THEN
             DoRemove = InterfaceDof(ind) 
           END IF
           IF( DoRemove ) THEN
             removed = removed + 1
             MortarBC % Active(ind) = .FALSE.
           END IF
         ELSE 
           NodeDist = DistVar % Values(k)
           MaxDist = MAX( MaxDist, NodeDist ) 
           MinDist = MIN( MinDist, NodeDist )

           DoAdd = ( NodeDist < -ValEps + DistOffset )            
           IF( DoAdd .AND. ConservativeAdd ) THEN
             DoAdd = InterfaceDof(ind)
           END IF
           IF( DoAdd ) THEN
             added = added + 1
             MortarBC % Active(ind) = .TRUE.
           END IF
         END IF
       END DO

       IF( InfoActive(20) ) THEN
         IF ( -HUGE(MaxDist) /= MaxDist ) THEN
           IF( MaxDist - MinDist >= 0.0_dp ) THEN
             PRINT *,'NormalContactSet Dist:',MinDist,MaxDist
           END IF
         END IF
         IF ( -HUGE(MaxLoad) /= MaxLoad) THEN
           IF( MaxLoad - MinLoad >= 0.0_dp ) THEN
             PRINT *,'NormalContactSet Load:',MinLoad,MaxLoad
           END IF
         END IF
         PRINT *,'NormalContactSet active count:',COUNT(MortarBC % Active)
         PRINT *,'NormalContactSet passive count:',COUNT(.NOT. MortarBC % Active)
       END IF

       IF(added > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Added ',added,' nodes to the set'
         CALL Info(Caller,Message,Level=6)
       END IF
       
       IF(removed > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Removed ',removed,' nodes from the set'
         CALL Info(Caller,Message,Level=6)
       END IF

       IF( InfoActive(30) ) THEN      
         PRINT *,'ModifiedActiveSet:',COUNT( MortarBC % Active ), SIZE( MortarBC % Active )
       END IF

       CALL Info(Caller,'Finished definition of Normal contact set',Level=25)
       
     END SUBROUTINE NormalContactSet



     ! If requested add new nodes to the contact set
     ! This would be typically done in order to make the elastic problem well defined
     ! Without any contact the bodies may float around.
     !---------------------------------------------------------------------------------
     SUBROUTINE IncreaseContactSet( LimitedMin )
       INTEGER :: LimitedMin

       REAL(KIND=dp), ALLOCATABLE :: DistArray(:)
       INTEGER, ALLOCATABLE :: IndArray(:)
       REAL(KIND=dp) :: Dist
       INTEGER :: i,j,ind,LimitedNow,NewNodes

       ! Nothing to do 
       IF( LimitedMin <= 0 ) RETURN

       LimitedNow = COUNT( MortarBC % active(DofN::cDofs) )      
       NewNodes = LimitedMin - LimitedNow
       IF( NewNodes <= 0 ) RETURN

       WRITE(Message,'(A,I0)') 'Initial number of contact nodes for '&
           //TRIM(VarName)//': ',LimitedNow 
       CALL Info(Caller,Message,Level=5)

       CALL Info(Caller,&
           'Setting '//I2S(NewNodes)//' additional contact nodes',Level=5)

       ALLOCATE( DistArray( NewNodes ), IndArray( NewNodes ) ) 
       DistArray = HUGE( DistArray ) 
       IndArray = 0

       ! Find additional contact nodes from the closest non-contact nodes
       DO i = 1,Projector % NumberOfRows
         ind = cDofs * (i-1) + DofN
         IF( MortarBC % Active(ind)  ) CYCLE

         IF( Projector % InvPerm(i) == 0 ) CYCLE
         j = DistVar % Perm(Projector % InvPerm(i))
         Dist = DistVar % Values(j)
 
         IF( Dist < DistArray(NewNodes) ) THEN
           DistArray(NewNodes) = Dist
           IndArray(NewNodes) = i

           ! Order the new nodes such that the last node always has the largest distance
           ! This way we only need to compare to the one distance when adding new nodes.
           DO j=1,NewNodes-1
             IF( DistArray(j) > DistArray(NewNodes) ) THEN
               Dist = DistArray(NewNodes)
               DistArray(NewNodes) = DistArray(j)
               DistArray(j) = Dist                
               ind = IndArray(NewNodes)
               IndArray(NewNodes) = IndArray(j)
               IndArray(j) = ind                
             END IF
           END DO
         END IF
       END DO

       IF( ANY( IndArray == 0 ) ) THEN
         CALL Fatal(Caller,'Could not define sufficient number of new nodes!')
       END IF

       WRITE(Message,'(A,ES12.4)') 'Maximum distance needed for new nodes:',DistArray(NewNodes)
       CALL Info(Caller,Message,Level=8)

       MortarBC % Active( cDofs*(IndArray-1)+DofN ) = .TRUE.

       DEALLOCATE( DistArray, IndArray ) 

     END SUBROUTINE IncreaseContactSet



     ! Sets the contact in the tangent direction(s) i.e. the stick condition.
     ! Stick condition in 1st and 2nd tangent condition are always the same. 
     !----------------------------------------------------------------------------------
     SUBROUTINE TangentContactSet() 
       
       INTEGER :: Removed0, Removed, Added
       REAL(KIND=dp) :: NodeLoad, TangentLoad, mustatic, mudynamic, stickcoeff, &
           Fstatic, Fdynamic, Ftangent, du(3), Slip
       INTEGER :: i,j,k,l,ind,IndN, IndT1, IndT2
       LOGICAL :: Found

       
       CALL Info(Caller,'Setting Tangent contact set',Level=20)
       
       IF( FrictionContact .AND. &
           ListGetLogical( BC,'Stick Contact Global',Found ) ) THEN
        
         ! Sum up global normal and slide forces
         DO i = 1,Projector % NumberOfRows
           j = Projector % InvPerm( i ) 
           IF( j == 0 ) CYCLE
           k = FieldPerm( j ) 
           IF( k == 0 ) CYCLE
           k = UseLoadVar % Perm(j)
                      
           ! If there is no contact there can be no stick either
           indN = cDofs * (i-1) + DofN
           IF( .NOT. MortarBC % Active(indN) ) CYCLE

           NodeLoad = UseLoadVar % Values(k)
           TangentLoad = SlipLoadVar % Values(k)
         
           mustatic = ListGetRealAtNode( BC,'Static Friction Coefficient', j )
           mudynamic = ListGetRealAtNode( BC,'Dynamic Friction Coefficient', j )
           IF( ( mustatic - mudynamic ) < -EPSILON(mustatic) ) THEN
             CALL Warn('TangentContactSet','Static friction coefficient should be larger than dynamic!')
           END IF
           
           Fstatic = Fstatic + mustatic * ABS( NodeLoad ) 
           Fdynamic = Fdynamic + mudynamic * ABS( NodeLoad )
           Ftangent = Ftangent + ABS( TangentLoad ) 
           IF( Ftangent > Fstatic ) THEN
             SlipContact = .TRUE.
             FrictionContact = .FALSE.
           ELSE 
             GOTO 100
           END IF
         END DO
       END IF

       
       ! For stick and tie contact inherit the active flag from the normal component
       IF( SlipContact ) THEN
         MortarBC % Active( DofT1 :: cDofs ) = .FALSE.
         IF( cDofs == 3 ) THEN
            MortarBC % Active( DofT2 :: cDofs ) = .FALSE.
          END IF
          GOTO 100 
       ELSE IF( StickContact .OR. TieContact ) THEN
         MortarBC % Active( DofT1 :: cDofs ) = MortarBC % Active( DofN :: cDofs )
         IF( cDofs == 3 ) THEN
           MortarBC % Active( DofT2 :: cDofs ) = MortarBC % Active( DofN :: cDofs ) 
         END IF
         GOTO 100
       END IF

       CALL Info('TangentContactSet','Setting the stick set tangent components',Level=10)

       Removed0 = 0
       Removed = 0
       Added = 0        

       ! Determine now whether we have contact or not
       DO i = 1,Projector % NumberOfRows
         j = Projector % InvPerm( i ) 
         IF( j == 0 ) CYCLE
         k = FieldPerm( j ) 
         IF( k == 0 ) CYCLE
         k = UseLoadVar % Perm(j)

         indN = cDofs * (i-1) + DofN
         indT1 = ind - DofN + DofT1
         IF(cDofs == 3 ) indT2 = ind - DofN + DofT2

         ! If there is no contact there can be no stick either
         IF( .NOT. MortarBC % Active(indN) ) THEN
           IF( MortarBC % Active(indT1) ) THEN
             removed0 = removed0 + 1
             MortarBC % Active(indT1) = .FALSE.
             IF( cDofs == 3 ) MortarBC % Active(indT2) = .FALSE.
           END IF
           CYCLE
         END IF

         ! Ok, we have normal contact what about stick
         ! Enforce stick condition
         !------------------------------------------------------
         coeff = ListGetRealAtNode( BC,'Stick Active Condition', j, Found )
         IF( Found .AND. coeff > 0.0_dp ) THEN
           IF( .NOT. MortarBC % Active(indT1) ) added = added + 1
           MortarBC % Active(indT1) = .TRUE.
           IF( cDofs == 3 ) MortarBC % Active(indT2) = .TRUE.
           CYCLE
         END IF

         ! Enforce no-stick condition (=slip)
         !------------------------------------------------------
         coeff = ListGetRealAtNode( BC,'Stick Passive Condition', j, Found )
         IF( Found .AND. coeff > 0.0_dp ) THEN
           IF( MortarBC % Active(IndT1) ) removed = removed + 1
           MortarBC % Active(indT1) = .FALSE.
           IF( cDofs == 3 ) MortarBC % Active(indT2) = .FALSE.
           CYCLE
         END IF

         ! Remove nodes with too large tangent force
         !--------------------------------------------------------------------------       

         NodeLoad = UseLoadVar % Values(k)
         TangentLoad = SlipLoadVar % Values(k)
         
         mustatic = ListGetRealAtNode( BC,'Static Friction Coefficient', j )
         mudynamic = ListGetRealAtNode( BC,'Dynamic Friction Coefficient', j )

         IF( ( mustatic - mudynamic ) < -EPSILON(mustatic) ) THEN
           CALL Warn('TangentContactSet','Static friction coefficient should be larger than dynamic!')
         END IF

         IF( MortarBC % Active(IndT1) ) THEN
           IF( TangentLoad > mustatic * ABS( NodeLoad ) ) THEN
             removed = removed + 1
             MortarBC % Active(indT1) = .FALSE.
             IF( cDofs == 3 ) MortarBC % Active(indT2) = .FALSE.
           END IF
         ELSE              
           stickcoeff = ListGetRealAtNode( BC,'Stick Contact Coefficient', j, Found )
           IF( Found ) THEN
             DO l=1,cDofs             
               du(l) = VeloVar % Values( cDofs*(k-1)+l ) 
             END DO
             IF( cDofs == 3 ) THEN
               Slip = SQRT(du(dofT1)**2 + du(DofT2)**2)
             ELSE
               Slip = ABS( du(dofT1) )
             END IF
             IF( stickcoeff * slip  < mudynamic * ABS( NodeLoad ) ) THEN
               added = added + 1
               MortarBC % Active(indT1) = .TRUE.
               IF( cDofs == 3 ) MortarBC % Active(indT2) = .TRUE.
             END IF
           END IF
         END IF
       END DO

       IF(added > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Added ',added,' nodes to the stick set'
         CALL Info(Caller,Message,Level=6)
       END IF
       
       IF(removed0 > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Removed ',removed0,' non-contact nodes from the stick set'
         CALL Info(Caller,Message,Level=6)
       END IF

       IF(removed > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Removed ',removed,' sliding nodes from the stick set'
         CALL Info(Caller,Message,Level=6)
       END IF


100    CALL Info(Caller,'Creating fields out of normal and stick contact sets',Level=10)

       DO i = 1, Projector % NumberOfRows
         j = Projector % InvPerm(i)
         IF( j == 0 ) CYCLE
         k = NormalActiveVar % Perm(j)

         IF( MortarBC % Active(cDofs*(i-1)+DofN) ) THEN
           NormalActiveVar % Values(k) = 1.0_dp
         ELSE
           NormalActiveVar % Values(k) = -1.0_dp
         END IF

         IF( MortarBC % Active(cDofs*(i-1)+DofT1) ) THEN
           StickActiveVar % Values(k) = 1.0_dp
         ELSE
           StickActiveVar % Values(k) = -1.0_dp
         END IF
       END DO

       IF( InfoActive(30) ) THEN
         PRINT *,'Active Tangent set:',COUNT( MortarBC % Active ) 
         CALL VariableValuesRange(NormalActiveVar,'NormalActive',AlwaysSerial=.TRUE.)
         CALL VariableValuesRange(StickActiveVar,'StickActive',AlwaysSerial=.TRUE.)
       END IF
         
     END SUBROUTINE TangentContactSet



     ! Sets the diagonal entry for slip in the tangent direction(s).
     ! This coefficient may be used to relax the stick condition, and also to
     ! revert back nodes from slip to stick set. 
     !----------------------------------------------------------------------------------
     SUBROUTINE StickCoefficientSet() 
       
       REAL(KIND=dp) :: NodeLoad, TangentLoad
       INTEGER :: i,j,k,ind,IndN, IndT1, IndT2
       LOGICAL :: Found

       CALL Info('StickCoefficientSet','Setting the stick coefficient entry for tangent components at stick',Level=10)

       ! Determine now whether we have contact or not
       DO i = 1,Projector % NumberOfRows
         j = Projector % InvPerm( i ) 
         IF( j == 0 ) CYCLE
         k = FieldPerm( j ) 
         IF( k == 0 ) CYCLE
         k = UseLoadVar % Perm(j)

         indN = cDofs * (i-1) + DofN
         indT1 = cDofs * (i-1) + DofT1
         IF(cDofs == 3 ) indT2 = cDofs * (i-1) + DofT2

         IF( .NOT. MortarBC % Active(indN) ) THEN
           ! If there is no contact there can be no stick either
           coeff = 0.0_dp            
         ELSE IF( .NOT. MortarBC % Active(indT1) ) THEN
           ! If there is no stick there can be no stick coefficient either
           coeff = 0.0_dp
         ELSE
           ! Get the stick coefficient
           coeff = ListGetRealAtNode( BC,'Stick Contact Coefficient', j )
         END IF

         MortarBC % Diag(indT1) = coeff
         IF( cDofs == 3 ) MortarBC % Diag(indT2) = coeff
       END DO

       IF(InfoActive(30)) THEN
         CALL VectorValuesRange( MortarBC % Diag, SIZE( MortarBC % Diag),'MortarBC Diag')
       END IF
       
     END SUBROUTINE StickCoefficientSet



     ! Here we eliminate the middle nodes from the higher order elements if they 
     ! are different than both nodes of which they are associated with.
     ! There is no way geometric information could be accurate enough to allow
     ! such contacts to exist.
     !---------------------------------------------------------------------------
     SUBROUTINE QuadraticContactSet()

       LOGICAL :: ElemActive(8)
       INTEGER :: i,j,k,n,added, removed, elem, elemcode, ElemInds(8)
       INTEGER, POINTER :: Indexes(:)

       added = 0
       removed = 0

       DO elem=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

         Element => Mesh % Elements( elem )         

         IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc_ind) % Tag ) CYCLE

         Indexes => Element % NodeIndexes         
         n = Element % TYPE % NumberOfNodes
         elemcode = Element % Type % ElementCode

         DO i=1,n
           ElemActive(i) = MortarBC % Active( ElemInds(i) ) 
           IF(j>0) THEN
             ElemInds(i) = cDofs * ( j - 1) + DofN
             ElemActive(i) = MortarBC % Active( ElemInds(i) ) 
           ELSE
             ElemActive(i) = .FALSE.
           END IF
         END DO

         SELECT CASE ( elemcode ) 

         CASE( 202, 303, 404 ) 
           CONTINUE

         CASE( 203 )
           IF( ( ElemActive(1) .EQV. ElemActive(2) ) &
               .AND. ( ElemActive(1) .NEQV. ElemActive(3) ) ) THEN
             MortarBC % Active( ElemInds(3) ) = ElemActive(1) 
             IF( ElemActive(1) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

         CASE( 306 ) 
           IF( ( ElemActive(1) .EQV. ElemActive(2) ) &
               .AND. ( ElemActive(1) .NEQV. ElemActive(4) ) ) THEN
             MortarBC % Active( ElemInds(4) ) = ElemActive(1) 
             IF( ElemActive(1) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

           IF( ( ElemActive(2) .EQV. ElemActive(3) ) &
               .AND. ( ElemActive(2) .NEQV. ElemActive(5) ) ) THEN
             MortarBC % Active( ElemInds(5) ) = ElemActive(2) 
             IF( ElemActive(2) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

           IF( ( ElemActive(3) .EQV. ElemActive(1) ) &
               .AND. ( ElemActive(3) .NEQV. ElemActive(6) ) ) THEN
             MortarBC % Active( ElemInds(6) ) = ElemActive(3) 
             IF( ElemActive(3) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

         CASE( 408 ) 
           IF( ( ElemActive(1) .EQV. ElemActive(2) ) &
               .AND. ( ElemActive(1) .NEQV. ElemActive(5) ) ) THEN
             MortarBC % Active( ElemInds(5) ) = ElemActive(1) 
             IF( ElemActive(1) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

           IF( ( ElemActive(2) .EQV. ElemActive(3) ) &
               .AND. ( ElemActive(2) .NEQV. ElemActive(6) ) ) THEN
             MortarBC % Active( ElemInds(6) ) = ElemActive(2) 
             IF( ElemActive(2) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

           IF( ( ElemActive(3) .EQV. ElemActive(4) ) &
               .AND. ( ElemActive(3) .NEQV. ElemActive(7) ) ) THEN
             MortarBC % Active( ElemInds(7) ) = ElemActive(3) 
             IF( ElemActive(3) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

           IF( ( ElemActive(4) .EQV. ElemActive(1) ) &
               .AND. ( ElemActive(4) .NEQV. ElemActive(8) ) ) THEN
             MortarBC % Active( ElemInds(8) ) = ElemActive(4) 
             IF( ElemActive(4) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

         CASE DEFAULT
           CALL Fatal(Caller,'Cannot deal with element: '//I2S(elemcode))

         END SELECT
       END DO

       IF(added > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Added ',added,' quadratic nodes to contact set'
         CALL Info(Caller,Message,Level=6)
       END IF

       IF(removed > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Removed ',removed,' quadratic nodes from contact set'
         CALL Info(Caller,Message,Level=6)
       END IF
         
     END SUBROUTINE QuadraticContactSet


     ! Project contact fields from slave to master
     !----------------------------------------------------------------------------------------
     SUBROUTINE ProjectFromSlaveToMaster()

       REAL(KIND=dp) :: Disp(3), Coord(3), PrevDisp(3), Velo(3), ContactVelo(3), &
           LocalNormal0(3), SlipCoord(3)
       REAL(KIND=dp), POINTER :: DispVals(:), PrevDispVals(:) 
       REAL(KIND=dp) :: MinDist, MaxDist, CoeffEps
       LOGICAL, ALLOCATABLE :: SlaveNode(:), NodeDone(:)
       REAL(KIND=dp), ALLOCATABLE :: CoeffTable(:), RealActive(:)
       INTEGER :: i,j,k,l,l2,Indexes(12)

       CALL Info(Caller,'Mapping entities from slave to master',Level=10)

       n = SIZE( FieldPerm )
       ALLOCATE( SlaveNode( n ) )
       SlaveNode = .FALSE.
       
       DO i=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

         Element => Mesh % Elements( i )                  
         IF( Element % BoundaryInfo % Constraint /= Model % BCs(bc_ind) % Tag ) CYCLE

         CurrentModel % CurrentElement => Element
         
         IF( pContact ) THEN
           n = mGetElementDOFs(pIndexes,Element)
           SlaveNode(pIndexes(1:n)) = .TRUE.
         ELSE
           SlaveNode( Element % NodeIndexes ) = .TRUE.
         END IF
       END DO

       IF( InfoActive(20) ) THEN
         n = COUNT( SlaveNode )
         CALL Info(Caller,'Number of dofs on slave side: '//I2S(n))
       END IF
         
       n = SIZE( DistVar % Values )
       ALLOCATE( CoeffTable( n ), NodeDone( n ) )
           
       CoeffTable = 0.0_dp
       NodeDone = .FALSE.
       

       DO i = 1,Projector % NumberOfRows             
         
         IF( Projector % InvPerm(i) == 0 ) CYCLE          
         l = DistVar % Perm( Projector % InvPerm(i) )
         
         IF(.NOT. pContact ) THEN
           IF( l > Mesh % NumberOfNodes ) CYCLE
         END IF

         DO j = Projector % Rows(i),Projector % Rows(i+1)-1
           k = Projector % Cols(j)

           IF(.NOT. pContact ) THEN
             IF( k > Mesh % NumberOfNodes ) CYCLE
           END IF                  
           
           IF( FieldPerm( k ) == 0 ) CYCLE               
           IF( SlaveNode( k ) ) CYCLE
           
           coeff = Projector % Values(j)
           
           l2 = DistVar % Perm( k )                
           
           IF(.NOT. NodeDone( l2 ) ) THEN
             DistVar % Values( l2 ) = 0.0_dp
             GapVar % Values( l2 ) = 0.0_dp
             NormalActiveVar % Values( l2 ) = 0.0_dp
             StickActiveVar % Values( l2 ) = 0.0_dp
             NormalLoadVar % Values( l2 ) = 0.0_dp
             SlipLoadVar % Values( l2 ) = 0.0_dp
             IF( CalculateVelocity ) THEN
               DO k=1,cDofs             
                 VeloVar % Values( cDofs*(l2-1)+k ) = 0.0_dp
               END DO
             END IF
             NodeDone( l2 ) = .TRUE.
           END IF

           CoeffTable( l2 ) = CoeffTable( l2 ) + coeff           
           DistVar % Values( l2 ) = DistVar % Values( l2 ) + coeff * DistVar % Values( l ) 
           GapVar % Values( l2 ) = GapVar % Values( l2 ) + coeff * GapVar % Values( l )
           NormalActiveVar % Values( l2 ) = NormalActiveVar % Values( l2 ) + coeff * NormalActiveVar % Values( l ) 
           StickActiveVar % Values( l2 ) = StickActiveVar % Values( l2 ) + coeff * StickActiveVar % Values( l ) 
           NormalLoadVar % Values( l2 ) = NormalLoadVar % Values( l2 ) + coeff * NormalLoadVar % Values( l ) 
           SlipLoadVar % Values( l2 ) = SlipLoadVar % Values( l2 ) + coeff * SlipLoadVar % Values( l ) 
           IF( CalculateVelocity ) THEN
             DO k=1,cDofs             
               VeloVar % Values( cDofs*(l2-1)+k ) = VeloVar % Values( cDofs*(l2-1)+k ) + &
                   coeff * VeloVar % Values( cDofs*(l-1)+k)
             END DO
           END IF
         END DO
       END DO
       
       CoeffEps = 1.0d-8 * MAXVAL( ABS( CoeffTable ) )
       DO i=1,SIZE( CoeffTable )            
         IF( NodeDone( i ) .AND. ( ABS( CoeffTable(i) ) > CoeffEps ) ) THEN
           DistVar % Values( i ) = DistVar % Values( i ) / CoeffTable( i ) 
           GapVar % Values( i ) = GapVar % Values( i ) / CoeffTable( i ) 
           NormalActiveVar % Values( i ) = NormalActiveVar % Values( i ) / CoeffTable( i ) 
           StickActiveVar % Values( i ) = StickActiveVar % Values( i ) / CoeffTable( i ) 

           IF( NormalActiveVar % Values( i ) >= 0.0_dp ) THEN
             NormalLoadVar % Values( i ) = NormalLoadVar % Values( i ) / CoeffTable( i ) 
             SlipLoadVar % Values( i ) = SlipLoadVar % Values( i ) / CoeffTable( i ) 
             IF( CalculateVelocity ) THEN
               DO k=1,cDofs
                 VeloVar % Values( cDofs*(i-1)+k ) = VeloVar % Values( cDofs*(i-1)+k ) / CoeffTable( i ) 
               END DO
             END IF
           ELSE
             NormalLoadVar % Values( i ) = 0.0_dp
             SlipLoadVar % Values( i ) = 0.0_dp
             IF( CalculateVelocity ) THEN
               DO k=1,cDofs
                 VeloVar % Values( cDofs*(i-1)+k ) = 0.0_dp
               END DO
             END IF             
           END IF

         END IF
       END DO

       DO i = 1, Projector % NumberOfRows
         j = Projector % InvPerm(i)
         IF( j == 0 ) CYCLE
         
         IF( .NOT. pContact ) THEN
           IF( j > Mesh % NumberOfNodes ) CYCLE
         END IF
         
         k = NormalActiveVar % Perm(j)
         
         IF( NormalActiveVar % Values( k ) < 0.0_dp ) THEN
           IF( CalculateVelocity ) THEN
             DO l=1,cDofs
               VeloVar % Values( cDofs*(k-1)+l ) = 0.0_dp
             END DO
           END IF
         END IF
       END DO
       
       IF( InfoActive(30) ) THEN
         CALL Info('PojectFromSlaveToMaster','Projecting fields')
         CALL VariableValuesRange(NormalLoadVar,'NormalLoadVar',AlwaysSerial=.TRUE.)
         CALL VariableValuesRange(SlipLoadVar,'SlipLoadVar',AlwaysSerial=.TRUE.)
         IF( CalculateVelocity ) THEN
           CALL VariableValuesRange(VeloVar,'VeloVar',AlwaysSerial=.TRUE.)
         END IF
       END IF

     END SUBROUTINE ProjectFromSlaveToMaster
   


     ! Set the friction in an implicit manner by copying matrix rows of the normal component
     ! to matrix rows of the tangential component multiplied by friction coefficient and 
     ! direction vector. 
     !---------------------------------------------------------------------------------------
     SUBROUTINE SetSlideFriction()

       REAL(KIND=dp), POINTER :: Values(:)
       LOGICAL, ALLOCATABLE :: NodeDone(:)
       REAL(KIND=dp) :: Coeff, ActiveLimit
       TYPE(Element_t), POINTER :: Element
       INTEGER, POINTER :: NodeIndexes(:)
       INTEGER :: i,j,k,k2,k3,l,l2,l3,n,t
       TYPE(Matrix_t), POINTER :: A
       LOGICAL :: Slave, Master, GivenDirection
       REAL(KIND=dp), POINTER :: VeloDir(:,:)
       REAL(KIND=dp) :: VeloCoeff(3),AbsVeloCoeff
       INTEGER :: VeloSign = 1
       CHARACTER(LEN=MAX_NAME_LEN) :: ContactVeloName


       IF(.NOT. ListCheckPresent( BC, 'Dynamic Friction Coefficient') ) RETURN
      
       CALL Info(Caller,'Setting contact friction for boundary',Level=10)

       ContactVeloName = 'Contact Velocity'
       GivenDirection = ListCheckPresent( BC, ContactVeloName )        
       IF( .NOT. GivenDirection .AND. TimeStep == 1 ) THEN
         ContactVeloName = 'Initial Contact Velocity'
         GivenDirection = ListCheckPresent( BC, ContactVeloName ) 
       END IF
       
       IF(.NOT. GivenDirection ) THEN
         IF(.NOT. ASSOCIATED( VeloVar ) ) THEN
           CALL Info(Caller,'Contact velocity not defined: Give "Contact Velocity" or',Level=3)
           CALL Info(Caller,'"Initial Contact Velocity" or use e.g. "Initial Stick Model = True"')
           CALL Fatal(Caller,'Cannot continue with friction without direction!')
         END IF
       END IF

       ActiveLimit = 0.0_dp
      
       Values => Solver % Matrix % values              
       ALLOCATE( NodeDone( SIZE( FieldPerm ) ) )
       A => Solver % Matrix        

       NodeDone = .FALSE.
       Coeff = 0.0_dp

       DO t = Mesh % NumberOfBulkElements+1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
         Element => Mesh % Elements(t)
         
         Model % CurrentElement => Element
         
         Slave = ( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag )
         Master = ( Element % BoundaryInfo % Constraint == Model % BCs(master_ind) % Tag )
         
         IF( .NOT. ( Slave .OR. Master ) ) CYCLE

         NodeIndexes => Element % NodeIndexes
         n = Element % TYPE % NumberOfNodes
         
         DO i = 1, n
           j = Nodeindexes(i) 

           IF( NodeDone( j ) ) CYCLE
           IF( FieldPerm( j ) == 0 ) CYCLE

           ! Skipping the nodes not in the boundary
           k = NormalActiveVar % Perm( j )
           IF( k == 0 ) CYCLE
           
           ! Skipping the nodes not in contact. 
           IF( NormalActiveVar % Values( k ) <= -ActiveLimit ) CYCLE
           
           ! skipping the nodes in tangent stick
           IF( StickActiveVar % Values( k ) >= ActiveLimit ) CYCLE

           NodeDone( j ) = .TRUE.

           IF( Slave ) THEN
             Coeff = ListGetRealAtNode( BC,& 
                 'Dynamic Friction Coefficient', j, Found )
           ELSE
             Coeff = ListGetRealAtNode( MasterBC,& 
                 'Dynamic Friction Coefficient', j, Found )
             ! If friction not found in master then use the friction coefficient of the slave 
             ! Ideally they should be the same. 
             IF(.NOT. Found ) THEN
               Coeff = ListGetRealAtNode( BC,& 
                   'Dynamic Friction Coefficient', j, Found )               
             END IF
           END IF

           ! There is no point of setting too small friction coefficient
           IF(ABS(Coeff) < 1.0d-10) CYCLE

           IF( ThisRotatedContact ) THEN
             Rotated = GetSolutionRotation(NTT, j )
             LocalNormal = NTT(:,1)
             LocalT1 = NTT(:,2)
             IF( cDofs == 3 ) LocalT2 = NTT(:,3)
           ELSE
             Rotated = .FALSE.
             LocalNormal = ContactNormal
             LocalT1 = ContactT1
             IF( cDofs == 3 ) LocalT2 = ContactT2 
           END IF
           
           VeloCoeff = 0.0_dp           
           VeloSign = 1

           IF( GivenDirection ) THEN
             IF( Slave ) THEN
               VeloDir => ListGetConstRealArray( BC, ContactVeloName )
             ELSE
               VeloDir => ListGetConstRealArray( MasterBC, ContactVeloName, Found)
               IF(.NOT. Found ) THEN
                 ! If velocity direction not found in master then use the opposite velocity of the slave
                 VeloDir => ListGetConstRealArray( BC, ContactVeloName )
                 VeloSign = -1
               END IF
             END IF
             VeloCoeff(DofT1) = SUM( VeloDir(1:3,1) * LocalT1 )
             IF( cDofs == 3 ) THEN
               VeloCoeff(DofT2) = SUM( VeloDir(1:3,1) * LocalT2 )
             END IF
           ELSE
             VeloCoeff(DofT1) = VeloVar % Values(cDofs*(k-1)+DofT1) 
             IF(cDofs==3) VeloCoeff(DofT2) = VeloVar % Values(cDofs*(k-1)+DofT2) 
             IF( .NOT. Slave .AND. .NOT. Rotated ) THEN
               VeloSign = -1
             END IF
           END IF

           ! Normalize coefficient to unity so that it only represents the direction of the force
           AbsVeloCoeff = SQRT( SUM( VeloCoeff**2 ) )
           IF( AbsVeloCoeff > TINY(AbsVeloCoeff) ) THEN
             VeloCoeff = VeloSign * VeloCoeff / AbsVeloCoeff
           ELSE
             CYCLE
           END IF

           ! Add the friction coefficient 
           VeloCoeff = Coeff * VeloCoeff 

           j = FieldPerm( j ) 
           k = cDOFs * (j-1) + DofN 

           k2 = cDOFs * (j-1) + DofT1 
           A % Rhs(k2) = A % Rhs(k2) - VeloCoeff(DofT1) * A % Rhs(k)

           IF( cDofs == 3 ) THEN
             k3 = cDOFs * (j-1) + DofT2
             A % Rhs(k3) = A % Rhs(k3) - VeloCoeff(DofT2) * A % Rhs(k)             
           END IF

           DO l = A % Rows(k),A % Rows(k+1)-1
             DO l2 = A % Rows(k2), A % Rows(k2+1)-1
               IF( A % Cols(l2) == A % Cols(l) ) EXIT
             END DO

             A % Values(l2) = A % Values(l2) - VeloCoeff(DofT1) * A % Values(l)
             
             IF( cDofs == 3 ) THEN
               DO l3 = A % Rows(k3), A % Rows(k3+1)-1
                 IF( A % Cols(l3) == A % Cols(l) ) EXIT
               END DO
               A % Values(l3) = A % Values(l3) - VeloCoeff(DofT2) * A % Values(l)
             END IF
           END DO
         END DO
       END DO
       
       n = COUNT( NodeDone ) 
       CALL Info(Caller,'Number of friction nodes: '//I2S(n),Level=10)
       
       DEALLOCATE( NodeDone )       

       IF( InfoActive(30) ) THEN
         CALL VectorValuesRange(A % Values,SIZE(A % Values),'A-friction')
       END IF
       
     END SUBROUTINE SetSlideFriction
     
   END SUBROUTINE DetermineContact


!------------------------------------------------------------------------------
!> Just toggles the initial system to harmonic one and back
!------------------------------------------------------------------------------
SUBROUTINE MergeSlaveSolvers( Solver, PreSolve )
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  LOGICAL :: PreSolve
  !------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Solver_t), POINTER :: Solver2
  TYPE(Matrix_t), POINTER :: A1, A2, A
  TYPE(Variable_t), POINTER :: betaVar
  REAL(KIND=dp), POINTER :: vals1(:), vals2(:), vals(:)
  REAL(KIND=dp) :: alpha, invAlpha
  INTEGER :: i,j,k,n,dofs1,dofs2
  INTEGER, POINTER :: SlaveSolverIndexes(:)
  INTEGER, POINTER :: perm1(:), perm2(:), perm(:)
  CHARACTER(:), ALLOCATABLE :: str
  LOGICAL :: Found

  SAVE :: A1, perm1, dofs1, vals1, &
      A2, perm2, dofs2, vals2, A, perm, vals, &
      alpha, invAlpha, betaVar

  CALL Info('MergeSlaveSolvers','Monolithic treatment of solvers')
  
  IF( .NOT. ASSOCIATED( Solver % Variable ) ) THEN
    CALL Fatal('MergeSlaveSolvers','Not applicable without a variable')
    RETURN    
  END IF
  IF( .NOT. ASSOCIATED( Solver % Matrix ) ) THEN
    CALL Fatal('MergeSlaveSolvers','Not applicable without a matrix')
    RETURN    
  END IF
  
  Params => Solver % Values  
  SlaveSolverIndexes => ListGetIntegerArray( Params,'Slave Solvers',Found )
  IF(.NOT. Found ) RETURN
  
  IF(SIZE(SlaveSolverIndexes) > 1 ) THEN
    CALL Warn('MergeSlaveSolvers','Cannot current only deal with one slave solver!')   
  END IF
  i = SlaveSolverIndexes(1)

  IF(PreSolve) THEN
    A1 => Solver % Matrix
    perm1 => Solver % Variable % Perm
    dofs1 = Solver % Variable % Dofs
    vals1 => Solver % Variable % Values

    Solver2 => CurrentModel % Solvers(i)
    A2 => Solver2 % Matrix
    perm2 => Solver2 % Variable % Perm
    dofs2 = Solver2 % Variable % Dofs
    vals2 => Solver2 % Variable % Values

    IF(dofs1 /= 1 .OR. dofs2 /= 1 ) THEN
      CALL Fatal('MergeSlaveSolver','Enabled for 1 dof only: '//I2S(dofs1)//','//I2S(dofs2))
    END IF

    alpha = ListGetCReal( Params,'Slave Field Scaling',Found)
    IF(.NOT. Found) alpha = 1.0_dp
    invAlpha = 1.0_dp / alpha

    betaVar => NULL()
    str = ListGetString( Params,'Slave Field Offset Variable',Found )
    IF( Found ) THEN
      betaVar => VariableGet(Solver % Mesh % Variables,str)
      IF(.NOT. ASSOCIATED(betaVar)) THEN
        CALL Fatal('MergeSlaveSolver','Could not find offset variable: '//TRIM(str))
      END IF
      IF(ANY(betaVar % Perm /= perm2 ) ) THEN
        CALL Fatal('MergeSlaveSolver','The offset field should have same permutation!')
      END IF      
    END IF
      
    A => AllocateMatrix()
    CALL CRS_MergeMatrix(A1, A2, C = A, PermA = perm1, PermB = perm2, PermC = perm) 

    ALLOCATE( A % Diag(A % NumberOfRows) )
    A % diag = 0
    CALL CRS_SortMatrix( A, .TRUE. )
    
    Solver % Matrix => A
    Solver % Variable % Perm => Perm

    ALLOCATE( vals( A % NumberOfRows) )
    Solver % Variable % Values => vals

    CALL MergeRhsAndSolutions()
  ELSE
    CALL UnmergeSolutions()

    
    Solver % Matrix => A1
    Solver % Variable % Values => vals1
    Solver % Variable % Perm => perm1

    DEALLOCATE(vals)
    DEALLOCATE(perm)
    CALL FreeMatrix(A)
  END IF

CONTAINS

  SUBROUTINE MergeRhsAndSolutions()

    INTEGER :: i,j,j1,j2
    REAL(KIND=dp), ALLOCATABLE :: rhsadd(:)
    REAL(KIND=dp) :: c1,x2
       
    CALL Info('MergeSlaveSolvers','Merging rhs and initial guess for monolithic solution!',Level=10)
    
    IF(.NOT. ASSOCIATED(A % rhs)) THEN
      ALLOCATE(A % rhs(A % NumberOfRows))
      A % RHS = 0._dp
    END IF

    IF(ASSOCIATED(betaVar)) THEN
      ALLOCATE(rhsadd(SIZE(betaVar % Values)))
      rhsadd = 0.0_dp
      CALL MatrixVectorMultiply( A2, betaVar % Values, rhsadd )
      rhsAdd = -InvAlpha * rhsAdd 
    END IF
      
    DO i=1,SIZE(perm)
      j = perm(i)
      IF(j==0) CYCLE

      j1 = perm1(i)
      j2 = perm2(i)

      IF(j1>0) A % rhs(j) = A1 % rhs(j1) 
      IF(j2>0) THEN
        A % rhs(j) = A % rhs(j) + InvAlpha * A2 % rhs(j2) 
        IF(ASSOCIATED(betaVar)) A % rhs(j) = A % rhs(j) + rhsAdd(j2)
      END IF

      vals(j) = 0.0_dp
      c1 = 0.0_dp      
      IF(j1>0) THEN
        c1 = 1.0_dp
        vals(j) = vals1(j1) 
      END IF      
      IF(j2>0) THEN
        x2 = vals2(j2)
        IF(ASSOCIATED(betaVar)) x2 = x2 - betaVar % Values(j2)
        vals(j) = (vals(j) + InvAlpha * x2)/(c1+1.0_dp) 
      END IF
    END DO    
    
  END SUBROUTINE MergeRhsAndSolutions

  
  SUBROUTINE UnmergeSolutions()

    INTEGER :: i,j,j1,j2

    CALL Info('MergeSlaveSolvers','Unmerging the solution back to composite solvers!',Level=10)
    
    DO i=1,SIZE(perm)
      j = perm(i)
      IF(j==0) CYCLE
      j1 = perm1(i)
      j2 = perm2(i)
      IF(j1>0) vals1(j1) = vals(j)
      IF(j2>0) THEN
        vals2(j2) = Alpha * vals(j)
      END IF
    END DO

    ! The permutation has been checked for these
    IF(ASSOCIATED(BetaVar)) vals2 = vals2 + betaVar % Values

    
  END SUBROUTINE UnmergeSolutions
    
END SUBROUTINE MergeSlaveSolvers



END MODULE ContactUtils
