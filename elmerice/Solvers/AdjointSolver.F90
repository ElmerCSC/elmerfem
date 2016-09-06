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
!> Compute the adjoint state of the Stokes equations.
!> 
SUBROUTINE AdjointSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!
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
  USE NavierStokes

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
  TYPE(ValueList_t),POINTER ::  BC,SolverParams
  TYPE(Nodes_t) :: ElementNodes
  TYPE(Element_t),POINTER :: Element
  TYPE(Variable_t), POINTER :: Sol
  TYPE(Variable_t), POINTER :: VelocitybSol
  REAL(KIND=dp), POINTER :: Vb(:)
  INTEGER, POINTER :: VbPerm(:)
  REAL(KIND=dp),POINTER :: ForceVector(:)
  integer :: t,n,NSDOFs,NVals,SolverInd
  REAL(KIND=dp),ALLOCATABLE :: STIFF(:,:),FORCE(:),ExtPressure(:),LoadVector(:,:),Alpha(:),Beta(:),SlipCoeff(:,:),w(:)
  Logical :: Gotit,GotForceBC,NormalTangential,Firsttime=.true.,UnFoundFatal=.TRUE.
  INTEGER, POINTER :: NodeIndexes(:),Perm(:)
  integer :: p,q,dim,c

  integer :: i,iii,jjj,k
  Real(KIND=dp) :: Unorm
  REAL(KIND=dp), allocatable :: FullMat(:,:)
  character(len=50) :: fo1
  character(LEN=MAX_NAME_LEN) :: SolName,SolverName


  save SolverName,Firsttime,SolverInd,STIFF,FORCE,ExtPressure,LoadVector,Alpha,Beta,SlipCoeff,w

  DIM = CoordinateSystemDimension()

  StiffMatrix => Solver % Matrix
  ForceVector => StiffMatrix % RHS

  Sol => Solver % Variable
  NSDOFs   =  Sol % DOFs
  Perm => Sol % Perm

  CALL InitializeToZero( StiffMatrix, ForceVector )

  if (Firsttime) then
          Firsttime=.False.
          SolverName = "Adjoint Solver"
          N = Solver % Mesh % MaxElementDOFs
          allocate(FORCE(  2*NSDOFs*N ),STIFF( 2*NSDOFs*N,2*NSDOFs*N ),ExtPressure(N), & 
                    SlipCoeff(3,N),LoadVector(4,N),Alpha(N),Beta(N),w(N))

          SolverParams => GetSolverParams()
          SolName = GetString( SolverParams,'Flow Solution Equation Name',Gotit)
          IF (.NOT.Gotit) Then
             CALL WARN(SolverName,'Keyword >Flow Solution Equation Name< not found in SolverParams')
             CALL WARN(SolverName,'Taking default value >Flow Solution<')
             WRITE(SolName,'(A)') 'Flow Solution'
          Endif

          DO i=1,Model % NumberOfSolvers
             if (TRIM(SolName) == ListGetString(Model % Solvers(i) % Values, 'Equation')) exit
          End do
          if (i.eq.(Model % NumberOfSolvers+1)) CALL FATAL(SolverName,'Could not find Flow Solver Equation Name')
          SolverInd=i
  end if



        NSSolver => Model % Solvers(SolverInd)
        IF(.NOT.ASSOCIATED(NSSolver % Matrix % BulkValues)) CALL FATAL(SolverName,'Flow Solver BulkValues not associated. & 
            Add >Calculate Loads = Logical true< In the flow solver section')

         InitMat => AllocateMatrix()
         InitMat % NumberOfRows =   NSSolver % Matrix % NumberOfRows
         InitMat % Values => NSSolver % Matrix % BulkValues
         InitMat % Rows => NSSolver % Matrix % Rows 
         InitMat % Cols => NSSolver % Matrix % Cols


         VelocitybSol => VariableGet( Solver % Mesh % Variables, 'Velocityb',UnFoundFatal=UnFoundFatal )
         Vb => VelocitybSol % Values
         VbPerm => VelocitybSol % Perm

         IF (VelocitybSol % DOFs.NE.(dim+1)) then
           WRITE(Message,'(A,I1,A,I1)') &
            'Variable Velocityb has ',VelocitybSol % DOFs,' DOFs, should be',dim+1
            CALL FATAL(SolverName,Message)
         End If

        TransMat => NULL()
        TransMat => CRS_Transpose(InitMat)

        NULLIFY( InitMat % Rows, InitMat % Cols, InitMat % Diag, InitMat % Values )
        CALL FreeMatrix( InitMat )

        CALL CRS_SortMatrix( TransMat , .true. )

        IF ( SIZE(StiffMatrix % Values) .NE. SIZE(TransMat % Values) ) THEN
              CALL WARN(SolverName,'StiffMatrix different size to TransMat. Is this the correct the body?')
        END IF

        StiffMatrix % Values = TransMat % Values
        StiffMatrix % Rows = TransMat % Rows
        StiffMatrix % Cols = TransMat % Cols
        StiffMatrix % Diag = TransMat % Diag
        ForceVector = 0.0
        Perm = NSSolver % Variable % Perm

        deallocate( TransMat % Rows, TransMat % Cols , TransMat % Values, TransMat %  Diag)
        nullify(TransMat)

      DO t = 1,Solver % Mesh % NumberOfBoundaryElements

        Element => GetBoundaryElement(t)
        IF ( .NOT. ActiveBoundaryElement() ) CYCLE

        n = GetElementNOFNodes()
!
!       The element type 101 (point element) can only be used
!       to set Dirichlet BCs, so skip ´em at this stage.
!
        IF ( GetElementFamily() == 1 ) CYCLE

        CALL GetElementNodes( ElementNodes )
        NodeIndexes => Element % NodeIndexes

        BC => GetBC()
        IF ( .NOT. ASSOCIATED(BC) ) CYCLE

        GotForceBC = GetLogical( BC, 'Adjoint Force BC',gotIt )

        IF (GotForceBC) Then
        LoadVector=0.0d0
        Alpha=0.0d0
        Beta=0.0d0
        STIFF=0.0
        FORCE=0.0
        ExtPressure=0.0
        SlipCoeff = 0.0d0

        ! I only see 1 case where we have to impose Neumann condition to the
        ! Adjoint system; the slip BC
        NormalTangential = GetLogical( BC, &
                         'Normal-Tangential Adjoint', GotIt )

        SlipCoeff(1,1:n) =  GetReal( BC, 'Slip Coefficient 1',GotIt )
        SlipCoeff(2,1:n) =  GetReal( BC, 'Slip Coefficient 2',GotIt )
        SlipCoeff(3,1:n) =  GetReal( BC, 'Slip Coefficient 3',GotIt )

        CALL NavierStokesBoundary(  STIFF, FORCE, &
             LoadVector, Alpha, Beta, ExtPressure, SlipCoeff, NormalTangential,   &
                Element, n, ElementNodes )

        CALL DefaultUpdateEquations( STIFF, FORCE )
       end if

      END DO
      
      !forcing of the adjoint system comes from the Velocityb variable computed
      !with the cost function
      c = dim + 1
      Do t=1,Solver%Mesh%NumberOfNodes
         Do i=1,c
           p=(Perm(t)-1)*c+i
           q=(VbPerm(t)-1)*c+i
           ForceVector(p)=Vb(q)
        End Do
      EndDo

      CALL FinishAssembly( Solver, ForceVector )
       
      CALL DefaultDirichletBCs()

      Unorm = DefaultSolve()

      Return
!------------------------------------------------------------------------------
END SUBROUTINE AdjointSolver
!------------------------------------------------------------------------------

