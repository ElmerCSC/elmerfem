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
! *  Authors: f. Gillet-Chaulet (IGE-Grenoble/France)
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!------------------------------------------------------------------------------
   SUBROUTINE Adjoint_LinearSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!> Compute the adjoint of the linear Solver and dirichlet conditions.
!
!     OUTPUT is : Solver % Variable the adjoint variable
!
!     INPUT PARAMETERS are:
!
!      In solver section:
!               Direct Equation Name = String 
!
!      Variables:
!                VarName_b
!
!******************************************************************************
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
   USE DefUtils
   IMPLICIT NONE
!------------------------------------------------------------------------------
   TYPE(Solver_t) :: Solver
   TYPE(Model_t) :: Model
   REAL(KIND=dp) :: dt
   LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
   TYPE(Solver_t),Pointer :: DSolver  ! Pointer to the Direct Solver
   TYPE(Variable_t), POINTER :: Sol   ! Solution Variable
   INTEGER :: DOFs

   TYPE(Matrix_t),POINTER :: InitMat,TransMat,StiffMatrix
   REAL(KIND=dp),POINTER :: ForceVector(:)
   INTEGER, POINTER :: Perm(:)

   TYPE(ValueList_t),POINTER ::  BC,BF,SolverParams
   TYPE(Element_t),POINTER :: Element
   INTEGER, POINTER :: NodeIndexes(:)

   ! Variables related to the frocing
   TYPE(Variable_t), POINTER :: VbSol
   REAL(KIND=dp), POINTER :: Vb(:)
   INTEGER, POINTER :: VbPerm(:)

   integer :: i,j,k,t,n
   INTEGER :: kmax
   Logical :: Gotit
   integer :: p,q,dim,c,m
   Real(KIND=dp) :: Unorm
   Real(KIND=dp) :: RotVec(3)

   ! Var related to Normal-Tangential stuff and Dirchlet BCs
   TYPE(ValueListEntry_t),POINTER :: NormalTangential,NormalTangentialC
   REAL(KIND=dp), POINTER :: BoundaryNormals(:,:)=> NULL(), &
                       BoundaryTangent1(:,:)=> NULL(), &
                       BoundaryTangent2(:,:)=> NULL()
   REAL(KIND=dp),ALLOCATABLE,SAVE :: Condition(:)
   INTEGER, POINTER,SAVE :: BoundaryReorder(:)=> NULL()
   INTEGER,SAVE :: NormalTangentialNOFNodes
   CHARACTER(LEN=MAX_NAME_LEN),SAVE :: NormalTangentialName,NormalTangentialNameb
   LOGICAL :: AnyNT
   LOGICAL :: Conditional,CheckNT

   CHARACTER(LEN=MAX_NAME_LEN),SAVE :: SolverName="Adjoint Solver" ! SolverName for messages
   CHARACTER(LEN=MAX_NAME_LEN) :: DEqName  ! Equation Name of Direct Solver
   CHARACTER(LEN=MAX_NAME_LEN) :: DVarName ! Var Name of The Direct Solver
   CHARACTER(LEN=MAX_NAME_LEN) :: FVarName ! Sensitivity Var Name (velocityb or DVarNameb)
   INTEGER, SAVE :: SolverInd ! Indice of the direct solver in the solver list
   LOGICAL, SAVE :: Firsttime=.TRUE.


   StiffMatrix => Solver % Matrix
   ForceVector => StiffMatrix % RHS

   Sol => Solver % Variable
   DOFs   =  Sol % DOFs
   Perm => Sol % Perm

   DIM = CoordinateSystemDimension()
  ! IF DIM = 3 and DOFs=2; e.g. solving SSA on a 3D mesh
  !Normal-Tangential can not be used => trick temporary set Model Dimension to 2
   IF (DIM.eq.(DOFs+1)) CurrentModel % Dimension = DOFs

   IF (Firsttime) then
      Firsttime=.False.

      N = Model % MaxElementNodes
      ALLOCATE(Condition(N))

      ! Get solver associated to the direct problem
      SolverParams => GetSolverParams(Solver)
      DEqName = ListGetString( SolverParams,'Direct Solver Equation Name',UnFoundFatal=.TRUE.)
      DO i=1,Model % NumberOfSolvers
          if (TRIM(DEqName) == ListGetString(Model % Solvers(i) % Values, 'Equation')) exit
      End do
      if (i.eq.(Model % NumberOfSolvers+1)) &
          CALL FATAL(SolverName,'Could not find Equation Name ' // TRIM(DEqName))
      SolverInd=i
      DSolver => Model % Solvers(SolverInd)

      !# Check For NT boundaries
      IF (DOFs>1) THEN
        NormalTangentialName = 'normal-tangential'
        IF ( SEQL(DSolver % Variable % Name, 'flow solution') ) THEN
         NormalTangentialName = TRIM(NormalTangentialName) // ' velocity'
        ELSE
         NormalTangentialName = TRIM(NormalTangentialName) // ' ' // &
                   GetVarName(DSolver % Variable)
        END IF

        AnyNT = ListGetLogicalAnyBC( Model, TRIM(NormalTangentialName) )
        IF (AnyNT) THEN

         !# We will stay in NT to impose dirichlet conditions
         CALL ListAddLogical(SolverParams,'Back Rotate N-T Solution',.FALSE.)

         !# !!! Has to be in lower case to copy item
         NormalTangentialNameb='normal-tangential' // ' ' // GetVarName(Solver % Variable)

         DO i=1,Model % NumberOfBCs
           BC => Model % BCs(i) % Values
           NormalTangential => ListFind( BC ,  NormalTangentialName , GotIt )
           IF (.NOT.Gotit) CYCLE

           WRITE(Message,'(A,I0)') 'Copy Normal-Tangential keyword in BC ',i
           CALL Info(SolverName,Message,level=4)
           CALL ListCopyItem( NormalTangential, BC, NormalTangentialNameb)
        
           NormalTangentialC => ListFind( BC ,  TRIM(NormalTangentialName) // ' condition' , GotIt )
           IF (.NOT.Gotit) CYCLE

           WRITE(Message,'(A,I0)') 'Copy Normal-Tangential Condition keyword in BC ',i
           CALL Info(SolverName,Message,level=4)
           CALL ListCopyItem( NormalTangentialC, BC, TRIM(NormalTangentialNameb) // ' condition' ) 
         END DO

         CALL CheckNormalTangentialBoundary( Model, TRIM(NormalTangentialNameb), &
            NormalTangentialNOFNodes, BoundaryReorder, &
            BoundaryNormals, BoundaryTangent1, BoundaryTangent2, dim )

        END IF
     END IF
   END IF

   ! Get Matrix from the Direct Solver
   DSolver => Model % Solvers(SolverInd)

   IF ( SEQL(DSolver % Variable % Name, 'flow solution') ) THEN
        DVarName = 'Velocity'
   ELSE
        DVarName = GetVarName(DSolver % Variable)
   END IF

   CALL DefaultInitialize()

   InitMat => AllocateMatrix()
   InitMat % NumberOfRows =   DSolver % Matrix % NumberOfRows
   InitMat % Values => DSolver % Matrix % Values
   InitMat % Rows => DSolver % Matrix % Rows 
   InitMat % Cols => DSolver % Matrix % Cols
   InitMat % Diag => DSolver % Matrix % Diag

   IF ( SEQL(DVarName,'velocity').OR.SEQL(DVarName,'ssavelocity')) THEN
    FVarName = 'Velocityb'
   ELSE
    FVarName = TRIM(DVarName) // 'b'
   ENDIF
   VbSol => VariableGet( Solver % Mesh % Variables, TRIM(FVarName),UnFoundFatal=.TRUE. )
   Vb => VbSol % Values
   VbPerm => VbSol % Perm
   IF (VbSol % DOFs.NE.DOFs) then
       WRITE(Message,'(A,I1,A,I1)') &
            'Variable Vb has ',VbSol % DOFs,' DOFs, should be',DOFs
       CALL FATAL(SolverName,Message)
   END IF

   TransMat => NULL()
   TransMat => CRS_Transpose(InitMat)

   NULLIFY( InitMat % Rows, InitMat % Cols, InitMat % Diag, InitMat % Values )
   CALL FreeMatrix( InitMat )

   CALL CRS_SortMatrix( TransMat , .true. )

   StiffMatrix % Values = TransMat % Values
   StiffMatrix % Rows = TransMat % Rows
   StiffMatrix % Cols = TransMat % Cols
   IF(ASSOCIATED(TransMat % Diag)) StiffMatrix % Diag = TransMat % Diag
   ForceVector = 0.0
   Perm = DSolver % Variable % Perm

   deallocate( TransMat % Rows, TransMat % Cols , TransMat % Values)
   IF(ASSOCIATED(TransMat % Diag)) DEALLOCATE(TransMat % Diag)
   nullify(TransMat)
      
   !forcing of the adjoint system comes from the Vb variable computed
   !with the cost function

   !! Vb is expressed in the model coodinate system => Rotate to NT
   CALL RotateNTSystemAll( Vb, VbPerm, DOFs )

   c = DOFs
   Do t=1,Solver%Mesh%NumberOfNodes
      Do i=1,c
         p=(Perm(t)-1)*c+i
         q=(VbPerm(t)-1)*c+i
         ForceVector(p)=Vb(q)
      End Do
   EndDo

   CALL FinishAssembly( Solver, ForceVector )

   Unorm = DefaultSolve()

  ! Go through BC to check if dirichlet was applied to direct solver
  ! Go to the boundary elements
   Do t=1,GetNOFBoundaryElements()
      Element => GetBoundaryElement(t)
      BC => GetBC(Element)
      IF (.NOT.ASSOCIATED( BC ) ) CYCLE

      NodeIndexes => Element % NodeIndexes
      n = Element % TYPE % NumberOfNodes
 
      ! Cf correspond lines in SolverUtils
      ! Check for nodes belonging to n-t boundary getting set by other bcs.
      ! -------------------------------------------------------------------
      CheckNT = .FALSE.
      IF ( NormalTangentialNOFNodes>0 .AND. DOFs>0 ) THEN
          CheckNT = .TRUE.
          IF ( ALL(BoundaryReorder(NodeIndexes(1:n))<1) ) CheckNT = .FALSE.
          IF ( ListGetLogical(BC,'normal-tangential' // ' ' // DVarName,Gotit)) CheckNT=.FALSE.
      END IF

      ! set BC for DOFs=1 or applied to the whole vector
      IF( ListCheckPresent( BC, DVarName ) ) THEN
         Condition(1:n) = ListGetReal( BC, &
               TRIM(DVarName) // ' Condition', n, NodeIndexes, Conditional )
         DO j=1,n
            IF ( Conditional .AND. Condition(j)<0._dp ) CYCLE
            DO i=1,DOFS
              Sol%Values(DOFs*(Perm(NodeIndexes(j))-1)+i)=0._dp
            END DO
         END DO
      ENDIF

      ! Do the same for each component
      IF (DOFs>1) THEN 
        DO i=1,DOFS
           IF( ListCheckPresent( BC,  ComponentName(DVarName,i)) ) THEN
             Condition(1:n) = ListGetReal( BC, &
               TRIM(ComponentName(DVarName,i)) // ' Condition', n, NodeIndexes, Conditional )

             DO j=1,n
               IF ( Conditional .AND. Condition(j)<0._dp ) CYCLE

               k = Perm(NodeIndexes(j))
               IF ( k > 0 ) THEN

                 m = 0
                 IF ( NormalTangentialNOFNodes>0 ) m=BoundaryReorder(NodeIndexes(j))
                 !! set in the NT system
                 IF ( m>0 .AND. CheckNT ) THEN
                 RotVec = 0._dp
                 RotVec(i) = 1._dp
                 CALL RotateNTSystem( RotVec, NodeIndexes(j) )

                 ! When cartesian component "DOF" is defined set the N-T component
                 ! closest to its direction. 
                 kmax = 1
                 DO k=2,dim
                   IF ( ABS(RotVec(k)) > ABS(RotVec(kmax)) ) THEN
                     kmax = k
                   END IF
                 END DO
                 Sol%Values(DOFs*(Perm(NodeIndexes(j))-1) + kmax)=0._dp
                 !else set the given DOF
                ELSE
                  Sol%Values(DOFs*(Perm(NodeIndexes(j))-1)+i)=0._dp
                END IF

               ENDIF
             END DO

           ENDIF

        END DO
      ENDIF
   END DO

   ! Go through the Body forces
   Do t=1,Solver % NumberOfActiveElements
      Element => GetActiveElement(t)

      BF => GetBodyForce(Element)
      IF (.NOT.ASSOCIATED( BF ) ) CYCLE

      NodeIndexes => Element % NodeIndexes
      n = Element % TYPE % NumberOfNodes

      ! set BC for DOFs=1 or applied to the whole vector
      IF( ListCheckPresent( BF, DVarName ) ) THEN
         Condition(1:n) = ListGetReal( BF, &
               TRIM(DVarName) // ' Condition', n, NodeIndexes, Conditional )
         DO j=1,n
            IF ( Conditional .AND. Condition(j)<0._dp ) CYCLE
            DO i=1,DOFS
              Sol%Values(DOFs*(Perm(NodeIndexes(j))-1)+i)=0._dp
            END DO
         END DO
      ENDIF

      ! Do the same for each component
      IF (DOFs>1) THEN
        DO i=1,DOFS
           IF( ListCheckPresent( BF,  ComponentName(DVarName,i)) ) THEN
             Condition(1:n) = ListGetReal( BF, &
               TRIM(ComponentName(DVarName,i)) // ' Condition', n, NodeIndexes, Conditional )

             DO j=1,n
               IF ( Conditional .AND. Condition(j)<0._dp ) CYCLE
               Sol%Values(DOFs*(Perm(NodeIndexes(j))-1)+i)=0._dp
             END DO
           END IF
        END DO
      END IF

   END DO

   ! Back Rotate to model coordinate system
   CALL BackRotateNTSystem(Sol%Values,Perm,DOFs)

   ! reset Dimension to DIM
   IF (DIM.eq.(DOFs+1)) CurrentModel % Dimension = DIM

   RETURN
!------------------------------------------------------------------------------
   END SUBROUTINE Adjoint_LinearSolver
!------------------------------------------------------------------------------

