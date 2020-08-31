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
! *  Authors: 
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!------------------------------------------------------------------------------
   SUBROUTINE AdjointSSA_AdjointSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!> Compute the adjoint state of the SSA equations.
!
!     OUTPUT is : Solver % Variable the adjoint state of the SSA problem
!
!     INPUT PARAMETERS are:
!
!      In solver section:
!               Flow Solution Equation Name = String (default 'SSA')
!
!      Variables
!                Velocityb (forcing for the adjoint pb)
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
   TYPE(Solver_t),Pointer :: NSSolver
   TYPE(Matrix_t),POINTER :: InitMat,TransMat,StiffMatrix
   TYPE(ValueList_t),POINTER ::  BC,BF,SolverParams
   TYPE(ValueListEntry_t),POINTER :: NormalTangential,NormalTangentialC
   character(LEN=MAX_NAME_LEN),SAVE :: NormalTangentialName,NormalTangentialNameb
   LOGICAL :: AnyNT
   TYPE(Element_t),POINTER :: Element
   TYPE(Variable_t), POINTER :: Sol
   TYPE(Variable_t), POINTER :: VelocitybSol
   REAL(KIND=dp), POINTER :: Vb(:)
   INTEGER, POINTER :: VbPerm(:)
   REAL(KIND=dp),POINTER :: ForceVector(:)
   integer :: i,j,k,t,n,NSDOFs
   INTEGER :: kmax
   Logical :: Gotit
   INTEGER, POINTER :: NodeIndexes(:),Perm(:)
   integer :: p,q,dim,c,m
   Real(KIND=dp) :: Unorm
   Real(KIND=dp) :: RotVec(3)
   character(LEN=MAX_NAME_LEN) :: SolName,NSVarName
   LOGICAL :: Conditional,CheckNT
   REAL(KIND=dp),ALLOCATABLE,SAVE :: Condition(:)

   INTEGER, POINTER,SAVE :: BoundaryReorder(:)=> NULL()
   INTEGER,SAVE :: NormalTangentialNOFNodes
   REAL(KIND=dp), POINTER :: BoundaryNormals(:,:)=> NULL(), &
                       BoundaryTangent1(:,:)=> NULL(), &
                       BoundaryTangent2(:,:)=> NULL()

   CHARACTER(LEN=MAX_NAME_LEN),SAVE :: SolverName="Adjoint Solver"
   INTEGER, SAVE :: SolverInd
   LOGICAL, SAVE :: Firsttime=.TRUE.


   DIM = CoordinateSystemDimension()

   StiffMatrix => Solver % Matrix
   ForceVector => StiffMatrix % RHS

   Sol => Solver % Variable
   NSDOFs   =  Sol % DOFs
   Perm => Sol % Perm

  ! IF DIM = 3 and NSDOFs=2; Normal-Tangential can not be used => trick temporary set
  ! Model Dimension to 2
   IF (DIM.eq.(NSDOFs+1)) CurrentModel % Dimension = NSDOFs


   IF (Firsttime) then
      Firsttime=.False.

      N = Model % MaxElementNodes
      ALLOCATE(Condition(N))

      SolverParams => GetSolverParams(Solver)
      SolName = ListGetString( SolverParams,'Flow Solution Equation Name',GotIt)
      IF (.NOT.Gotit) Then
        CALL WARN(SolverName,'Keyword >Flow Solution Equation Name< not found in SolverParams')
        CALL WARN(SolverName,'Taking default value >SSA<')
        WRITE(SolName,'(A)') 'SSA'
      Endif
      DO i=1,Model % NumberOfSolvers
          if (TRIM(SolName) == ListGetString(Model % Solvers(i) % Values, 'Equation')) exit
      End do
      if (i.eq.(Model % NumberOfSolvers+1)) &
          CALL FATAL(SolverName,'Could not find Equation Name ' // TRIM(SolName))

      SolverInd=i
      NSSolver => Model % Solvers(SolverInd)
 
      !# Check For NT boundaries
      NormalTangentialName = 'normal-tangential'
      IF ( SEQL(NSSolver % Variable % Name, 'flow solution') ) THEN
       NormalTangentialName = TRIM(NormalTangentialName) // ' velocity'
      ELSE
       NormalTangentialName = TRIM(NormalTangentialName) // ' ' // &
                   GetVarName(NSSolver % Variable)
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


   NSSolver => Model % Solvers(SolverInd)

   IF ( SEQL(NSSolver % Variable % Name, 'flow solution') ) THEN
        NSVarName = 'Velocity'
   ELSE
        NSVarName = GetVarName(NSSolver % Variable)
   END IF

   CALL DefaultInitialize()

   InitMat => AllocateMatrix()
   InitMat % NumberOfRows =   NSSolver % Matrix % NumberOfRows
   InitMat % Values => NSSolver % Matrix % Values
   InitMat % Rows => NSSolver % Matrix % Rows 
   InitMat % Cols => NSSolver % Matrix % Cols
   InitMat % Diag => NSSolver % Matrix % Diag


   VelocitybSol => VariableGet( Solver % Mesh % Variables, 'Velocityb'  )
   IF ( ASSOCIATED( VelocitybSol ) ) THEN
     Vb => VelocitybSol % Values
     VbPerm => VelocitybSol % Perm
   ELSE
     WRITE(Message,'(A)') 'No variable > Velocityb < found'
     CALL FATAL(SolverName,Message)
   END IF  
   IF (VelocitybSol % DOFs.NE.NSDOFs) then
       WRITE(Message,'(A,I1,A,I1)') &
            'Variable Velocityb has ',VelocitybSol % DOFs,' DOFs, should be',NSDOFs
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
   Perm = NSSolver % Variable % Perm

   deallocate( TransMat % Rows, TransMat % Cols , TransMat % Values)
   IF(ASSOCIATED(TransMat % Diag)) DEALLOCATE(TransMat % Diag)
   nullify(TransMat)
      
   !forcing of the adjoint system comes from the Velocityb variable computed
   !with the cost function

   !! Vb is expressed in the model coodinate system => Rotate to NT
   CALL RotateNTSystemAll( Vb, VbPerm, NSDOFs )

   c = NSDOFs
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
      IF ( NormalTangentialNOFNodes>0 .AND. NSDOFs>0 ) THEN
          CheckNT = .TRUE.
          IF ( ALL(BoundaryReorder(NodeIndexes(1:n))<1) ) CheckNT = .FALSE.
          IF ( ListGetLogical(BC,'normal-tangential' // ' ' // NSVarName,Gotit)) CheckNT=.FALSE.
      END IF

      ! set BC for NSDOFs=1 or applied to the whole vector
      IF( ListCheckPresent( BC, NSVarName ) ) THEN
         Condition(1:n) = ListGetReal( BC, &
               TRIM(NSVarName) // ' Condition', n, NodeIndexes, Conditional )
         DO j=1,n
            IF ( Conditional .AND. Condition(j)<0._dp ) CYCLE
            DO i=1,NSDOFS
              Sol%Values(NSDOFs*(Perm(NodeIndexes(j))-1)+i)=0._dp
            END DO
         END DO
      ENDIF

      ! Do the same for each component
      IF (NSDOFs>1) THEN 
        DO i=1,NSDOFS
           IF( ListCheckPresent( BC,  ComponentName(NSVarName,i)) ) THEN
             Condition(1:n) = ListGetReal( BC, &
               TRIM(ComponentName(NSVarName,i)) // ' Condition', n, NodeIndexes, Conditional )

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
                 Sol%Values(NSDOFs*(Perm(NodeIndexes(j))-1) + kmax)=0._dp
                 !else set the given DOF
                ELSE
                  Sol%Values(NSDOFs*(Perm(NodeIndexes(j))-1)+i)=0._dp
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

      ! set BC for NSDOFs=1 or applied to the whole vector
      IF( ListCheckPresent( BF, NSVarName ) ) THEN
         Condition(1:n) = ListGetReal( BF, &
               TRIM(NSVarName) // ' Condition', n, NodeIndexes, Conditional )
         DO j=1,n
            IF ( Conditional .AND. Condition(j)<0._dp ) CYCLE
            DO i=1,NSDOFS
              Sol%Values(NSDOFs*(Perm(NodeIndexes(j))-1)+i)=0._dp
            END DO
         END DO
      ENDIF

      ! Do the same for each component
      IF (NSDOFs>1) THEN
        DO i=1,NSDOFS
           IF( ListCheckPresent( BF,  ComponentName(NSVarName,i)) ) THEN
             Condition(1:n) = ListGetReal( BF, &
               TRIM(ComponentName(NSVarName,i)) // ' Condition', n, NodeIndexes, Conditional )

             DO j=1,n
               IF ( Conditional .AND. Condition(j)<0._dp ) CYCLE
               Sol%Values(NSDOFs*(Perm(NodeIndexes(j))-1)+i)=0._dp
             END DO
           END IF
        END DO
      END IF

   END DO

   ! Back Rotate to model coordinate system
   CALL BackRotateNTSystem(Sol%Values,Perm,NSDOFs)

   ! reset Dimension to DIM
   IF (DIM.eq.(NSDOFs+1)) CurrentModel % Dimension = DIM

   RETURN
!------------------------------------------------------------------------------
   END SUBROUTINE AdjointSSA_AdjointSolver
!------------------------------------------------------------------------------

