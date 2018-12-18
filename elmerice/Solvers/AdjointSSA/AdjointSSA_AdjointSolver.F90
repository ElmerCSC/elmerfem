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
!> 
SUBROUTINE AdjointSSA_AdjointSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!> Compute the adjoint state of the SSA equations.
!
!     OUTPUT is : Solver % Variable the adjoint sate of the SSA problem
!
!     INPUT PARAMETERS are:
!
!      In solver section:
!               Flow Solution Equation Name = String (defualt 'SSA')
!
!      Variables
!                Velocityb (forcing for the adjoint pb)
!                Bulk Values of the SSA pb (need 'compute Loads = True' in the SSA solver
!                
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
  Logical :: Gotit,GotForceBC,NormalTangential,Firsttime=.true.
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

  ! IF DIM = 3 and NSDOFs=2; Normal-Tangential can not be used => trick temporary set
  ! Model Dimension to 2
  IF (DIM.eq.(NSDOFs+1)) CurrentModel % Dimension = NSDOFs

  !CALL InitializeToZero( StiffMatrix, ForceVector )
  CALL DefaultInitialize()

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
             CALL WARN(SolverName,'Taking default value >SSA<')
             WRITE(SolName,'(A)') 'SSA'
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
         InitMat % Diag => NSSolver % Matrix % Diag


         VelocitybSol => VariableGet( Solver % Mesh % Variables, 'Velocityb'  )
         IF ( ASSOCIATED( VelocitybSol ) ) THEN
            Vb => VelocitybSol % Values
            VbPerm => VelocitybSol % Perm
         ELSE
            WRITE(Message,'(A)') &
                               'No variable > Velocityb < found'
            CALL FATAL(SolverName,Message)
         END IF  
         IF (VelocitybSol % DOFs.NE.NSDOFs) then
           WRITE(Message,'(A,I1,A,I1)') &
            'Variable Velocityb has ',VelocitybSol % DOFs,' DOFs, should be',NSDOFs
            CALL FATAL(SolverName,Message)
         End If

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
      c = NSDOFs
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

     ! reset Dimension to DIM
      IF (DIM.eq.(NSDOFs+1)) CurrentModel % Dimension = DIM

      Return
!------------------------------------------------------------------------------
END SUBROUTINE AdjointSSA_AdjointSolver
!------------------------------------------------------------------------------

