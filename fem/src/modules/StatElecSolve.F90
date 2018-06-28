!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under The terms of the GNU General Public License
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
! *  Authors: Juha Ruokolainen, Leila Puska, Antti Pursula, Peter R�back
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 20 Jun 2002
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Initialization of the primary solver, i.e. StatElecSolver.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE StatElecSolver_Init( Model,Solver,dt,TransientSimulation)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Model_t)  :: Model
    TYPE(Solver_t), TARGET :: Solver
    LOGICAL ::  TransientSimulation
    REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
    LOGICAL :: Found, Calculate, CalculateCapMatrix
    TYPE(ValueList_t), POINTER :: Params
    CHARACTER(LEN=MAX_NAME_LEN) :: VariableName
    INTEGER :: i,dim

    Params => GetSolverParams()
    dim = CoordinateSystemDimension()

    IF (ListGetLogical(Params,'Calculate Electric Energy',Found)) &
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        'Electric Energy Density' )
    
    Calculate = ListGetLogical(Params,'Calculate Electric Field',Found)
    IF( Calculate ) THEN
      IF( Dim == 2 ) THEN
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
            'elfield[Electric Field:2]' )
      ELSE
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
            'elfield[Electric Field:3]' )
      END IF
    END IF
    
    Calculate = ListGetLogical(Params,'Calculate Electric Flux',Found)
    IF( Calculate ) THEN
      IF( Dim == 2 ) THEN
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
            'elflux[Electric Flux:2]' )
      ELSE
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
            'elflux[Electric Flux:3]' )
      END IF
    END IF


    ! If computation of capacitance matrix is requested then compute 
    ! set the flag for load computation also.
    !------------------------------------------------------------------
    CalculateCapMatrix = ListGetLogical( Params, &
        'Calculate Capacitance Matrix', Found )
    IF(.NOT. Found ) THEN
      DO i = 1, Model % NumberOfEquations 
        CalculateCapMatrix = ListGetLogical( Model % Equations(i) % Values, &
            'Calculate Capacitance Matrix', Found )
        IF ( CalculateCapMatrix ) THEN
          CALL ListAddLogical( Params,'Calculate Capacitance Matrix',.TRUE.)
          EXIT
        END IF
      END DO
    END IF
    IF( CalculateCapMatrix ) THEN
      CALL ListAddLogical( Params,'Calculate Loads',.TRUE.)
    END IF

    CALL ListAddInteger( Params,'Time Derivative Order', 0 )

!------------------------------------------------------------------------------
END SUBROUTINE StatElecSolver_Init
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>  Solves the Poisson equation for the electric potential and compute the 
!>  electric field, flux, energy and capacitance as requested using nodal averaging.
!
!>  Note that the permittivity of vacuum is divided into the right hand side of the 
!>  equation. This has to be accounted for in setting the body forces and
!>  assigning flux boundary conditions
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE StatElecSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE Types
  USE Lists 
  USE Integration
  USE ElementDescription
  USE Differentials
  USE SolverUtils
  USE ElementUtils
  USE Adaptive
  USE DefUtils
!$ USE omp_lib ! Include module conditionally
  
  IMPLICIT NONE
!------------------------------------------------------------------------------
 
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET:: Solver
  REAL (KIND=DP) :: dt
  LOGICAL :: TransientSimulation
      
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(Matrix_t), POINTER  :: StiffMatrix
  TYPE(Element_t), POINTER :: CurrentElement
  TYPE(Variable_t), POINTER :: TimeVar, Var
  TYPE(Nodes_t) :: ElementNodes
  
  REAL (KIND=DP), POINTER :: ForceVector(:), Potential(:), Displacement(:,:)
  REAL (KIND=DP), POINTER :: Field(:),Flux(:),Energy(:),PermIso(:)
  REAL (KIND=dp), POINTER CONTIG :: PValues(:)
  REAL (KIND=dp), POINTER :: Charges(:)
  REAL (KIND=DP), POINTER :: Pwrk(:,:,:), Pz_w(:,:,:)
  REAL (KIND=DP), ALLOCATABLE :: CapMatrix(:,:),CapMatrixPara(:,:)
  REAL (KIND=DP), ALLOCATABLE ::  Permittivity(:,:,:), PiezoCoeff(:,:,:), &
      LocalStiffMatrix(:,:), Load(:), LocalForce(:), PotDiff(:), &
      Alpha(:), Beta(:),LayerH(:),LayerV(:), Basis(:), dBasisdx(:,:)
  
  REAL(KIND=dp) :: RelPerm1, RelPerm2
  REAL(KIND=dp) :: PermittivityOfVacuum, Norm, RelativeChange
  
#ifdef USE_ISO_C_BINDINGS
  REAL (KIND=DP) :: Wetot, at0, ss
  REAL (KIND=DP) :: at, st, PotentialDifference, Capacitance
#else
  REAL (KIND=DP) :: Wetot, at0, RealTime, ss
  REAL (KIND=DP) :: at, st, CPUTime, PotentialDifference, Capacitance
#endif
  REAL (KIND=DP) :: MinPotential, MaxPotential
  
  INTEGER, POINTER :: NodeIndexes(:), CapBodyIndex(:), Ivals(:)
  INTEGER, POINTER :: PotentialPerm(:), EnergyPerm(:), SurfPerm(:)
  INTEGER, POINTER :: FieldPerm(:), FluxPerm(:)
  INTEGER :: CapBodies, CapBody, Permi, Permj, iter, MaxIterations
  INTEGER :: i, j, k, l, m, istat, bf_id, LocalNodes, DIM, NonlinearIter, &
      RelIntegOrder, nsize, N, ntot, t, TID
  
  LOGICAL :: AllocationsDone = .FALSE., gotIt, FluxBC, OpenBc, LayerBC
  LOGICAL :: CalculateField, CalculateFlux, CalculateEnergy
  LOGICAL :: CalculateCapMatrix, ConstantWeights
  LOGICAL :: PiezoMaterial
  LOGICAL :: ConstantBulk = .FALSE., AssemblyDone = .FALSE.

  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Matrix_t), POINTER :: CM
  LOGICAL, ALLOCATABLE :: Done(:)
  LOGICAL :: DoneL
  
  CHARACTER(LEN=MAX_NAME_LEN) :: CapMatrixFile, Name, VarName
  TYPE(ValueList_t), POINTER :: Params, BC 
  
  SAVE LocalStiffMatrix, Load, LocalForce, PotDiff, &
      ElementNodes, CalculateFlux, CalculateEnergy, &
      AllocationsDone, Permittivity, &
      CapBodies, CalculateCapMatrix, CapBodyIndex, &
      CapMatrix, CalculateField, CapMatrixFile, ConstantWeights, &
      PiezoCoeff, PiezoMaterial, Displacement, Pwrk, Pz_w, &
      ConstantBulk, AssemblyDone, Alpha, Beta, LayerH, LayerV, &
      PermIso, Charges, Basis, dBasisdx
  
  ! Variables private to the thread (local storage)
  !$omp threadprivate(ElementNodes, Permittivity, LocalForce, Alpha, &
  !$omp               Beta, LayerV, LayerH, PermIso, LocalStiffMatrix, &
  !$omp               Load, PotDiff, Displacement, PiezoCoeff, &
  !$omp               CapBodyIndex, Charges, CapMatrix, Pwrk, Pz_w, &
  !$omp               Basis, dBasisdx, PiezoMaterial)
  
  INTERFACE
    FUNCTION ElectricBoundaryResidual( Model,Edge,Mesh,Quant,Perm,Gnorm ) RESULT(Indicator)
      USE Types
      TYPE(Element_t), POINTER :: Edge
      TYPE(Model_t) :: Model
      TYPE(Mesh_t), POINTER :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
      INTEGER :: Perm(:)
    END FUNCTION ElectricBoundaryResidual
    
    FUNCTION ElectricEdgeResidual( Model,Edge,Mesh,Quant,Perm ) RESULT(Indicator)
      USE Types
      TYPE(Element_t), POINTER :: Edge
      TYPE(Model_t) :: Model
      TYPE(Mesh_t), POINTER :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2)
      INTEGER :: Perm(:)
    END FUNCTION ElectricEdgeResidual
    
    FUNCTION ElectricInsideResidual( Model,Element,Mesh,Quant,Perm, Fnorm ) RESULT(Indicator)
      USE Types
      TYPE(Element_t), POINTER :: Element
      TYPE(Model_t) :: Model
      TYPE(Mesh_t), POINTER :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
      INTEGER :: Perm(:)
    END FUNCTION ElectricInsideResidual
  END INTERFACE
  
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------

  Params => GetSolverParams()

  PotentialPerm => Solver % Variable % Perm
  IF ( COUNT(PotentialPerm > 0) == 0 ) RETURN

  Potential => Solver % Variable % Values
  StiffMatrix => Solver % Matrix
  ForceVector => StiffMatrix % RHS
  VarName = GetVarName( Solver % Variable )
  Mesh => Solver % Mesh

  Norm = Solver % Variable % Norm
  DIM = CoordinateSystemDimension()

!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
    N = Mesh % MaxElementDOFs

    !$omp parallel shared(dim, n) private(istat) default(none)
    ALLOCATE( ElementNodes % x(N),   &
        ElementNodes % y(N),   &
        ElementNodes % z(N),   &
        Permittivity(3,3,N),   &
        LocalForce(N),         &
        Alpha(N),              &
        Beta(N),               &
        LayerV(N),             &
        LayerH(N),             &
        PermIso(N),            &
        LocalStiffMatrix(N,N), &
        Load(N),               &
        PotDiff(N),            &
        Displacement(N,Dim),   &
        Basis(N),              &
        dBasisdx(N,DIM),       &
        STAT=istat )
    
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'StatElecSolve', 'Memory allocation error 1' )
    END IF
    
    ALLOCATE( PiezoCoeff(DIM,2*DIM,N), STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'StatElecSolve', 'Memory allocation error 2' )
    END IF

    NULLIFY( Pwrk )
    NULLIFY( Pz_w )
    !$omp end parallel

    CalculateField = ListGetLogical( Params,'Calculate Electric Field', GotIt )    
    CalculateFlux = ListGetLogical( Params,'Calculate Electric Flux', GotIt )    
    CalculateEnergy = ListGetLogical( Params,'Calculate Electric Energy', GotIt )     

    ConstantBulk = ListGetLogical( Params,'Constant Bulk System',GotIt)
    ConstantBulk = ConstantBulk .OR. ListGetLogical( Params,'Save Bulk System',GotIt)

    ConstantWeights = ListGetLogical( Params,'Constant Weights', GotIt )

    CalculateCapMatrix = ListGetLogical( Params, &
        'Calculate Capacitance Matrix', GotIt )
        
    IF(CalculateCapMatrix) THEN
      ConstantBulk = .TRUE.
      CapBodies = ListGetInteger( Params, 'Capacitance Bodies',GotIt)
      IF(.NOT. GotIt) THEN
        DO i = 1, Model % NumberOfBCs
          j = ListGetInteger( Model % BCs(i) % Values, &
            'Capacitance Body', GotIt )
          IF( j > CapBodies ) CapBodies = j
        END DO
      END IF	
      
      IF( CapBodies == 0 ) THEN
        CALL Fatal('StatElecSolve',&
            'Capacitance calculation requested without any > Capacitance Body <')
      END IF

      nsize = SIZE( Potential )
      !$omp parallel shared(nsize, CapBodies) private(istat) default(none)
      ALLOCATE(CapBodyIndex(nsize), &
          Charges(nsize), &
          CapMatrix(CapBodies,CapBodies), &
          STAT = istat)
      IF ( istat /= 0 ) THEN
        CALL Fatal( 'StatElecSolve', 'Memory allocation error 3' )
      END IF
      
      CapMatrix = 0.0_dp
      CapBodyIndex = 0
      !$omp end parallel
    END IF

    IF ( .NOT.ASSOCIATED( StiffMatrix % MassValues ) ) THEN
      ALLOCATE( StiffMatrix % Massvalues( Model % NumberOfNodes ) )
      StiffMatrix % MassValues = 0.0_dp
    END IF

    AllocationsDone = .TRUE.
  END IF

!------------------------------------------------------------------------------
!  Get the result fields
!  These should have been automatically created by the 
!  'Exported Variables' defined in the _init section
!------------------------------------------------------------------------------

  IF ( CalculateField ) THEN
    Var => VariableGet( Solver % Mesh % Variables,'elfield')
    IF( ASSOCIATED( Var) ) THEN
      Field => Var % Values
    ELSE
      CALL Fatal('StatElecSolver','Electric Field does not exist')
    END IF
  END IF
  
  IF ( CalculateFlux ) THEN
    Var => VariableGet( Solver % Mesh % Variables,'elflux')
    IF( ASSOCIATED( Var ) ) THEN
      Flux => Var % Values
    ELSE
      CALL Fatal('StatElecSolver','Electric Flux does not exist')
    END IF
  END IF

  IF ( CalculateEnergy ) THEN
    Var => VariableGet( Solver % Mesh % Variables,'Electric Energy Density')
    IF( ASSOCIATED( Var ) ) THEN    
      Energy => Var % Values
    ELSE
      CALL Fatal('StatElecSolver','Electric Energy Density does not exist')
    END IF
  END IF
   
  RelIntegOrder = ListGetInteger( Params,'Relative Integration Order',GotIt)

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

  PermittivityOfVacuum = ListGetConstReal( Model % Constants, &
      'Permittivity Of Vacuum',gotIt )
  IF ( .NOT.gotIt ) PermittivityOfVacuum = 1.0_dp
  
  NonlinearIter = ListGetInteger( Params, &
      'Nonlinear System Max Iterations', GotIt )
  IF ( .NOT. GotIt ) NonlinearIter = 1
 
  IF(CalculateCapMatrix) NonlinearIter = CapBodies
  
  CALL DefaultInitialize()
!------------------------------------------------------------------------------
  CALL Info( 'StatElecSolve', '-------------------------------------',Level=4 )
  CALL Info( 'StatElecSolve', 'STATELEC SOLVER:  ', Level=4 )
  CALL Info( 'StatElecSolve', '-------------------------------------',Level=4 )

  CALL DefaultStart()

  
  DO iter = 1, NonlinearIter
     at  = CPUTime()
     at0 = RealTime()

     IF ( NonlinearIter > 1 ) THEN
        WRITE( Message, '(a,I0)' ) 'Electrostatic iteration: ', iter
        CALL Info( 'StatElecSolve', ' ', LEVEL=4 )
        CALL Info( 'StatElecSolve', Message, LEVEL=4 )
     END IF
     CALL Info( 'StatElecSolve', 'Starting Assembly...', Level=4 )

!------------------------------------------------------------------------------
!    Do the assembly
!------------------------------------------------------------------------------
    IF ( ConstantBulk .AND. AssemblyDone ) THEN
      Solver % Matrix % RHS = Solver % Matrix % BulkRHS
      Solver % Matrix % Values = Solver % Matrix % BulkValues
    ELSE
      CALL BulkAssembly()
      CALL DefaultFinishBulkAssembly()
      AssemblyDone = .TRUE.
    END IF

    CALL BoundaryAssembly()

!------------------------------------------------------------------------------
!    Solve the system and we are done.
!------------------------------------------------------------------------------
    st = CPUTime()

    Norm = DefaultSolve()
    RelativeChange = Solver % Variable % NonlinChange
    st = CPUTime() - st
    WRITE( Message, * ) 'Solve (s)             :',st
    CALL Info( 'StatElecSolve', Message, Level=4 )

!------------------------------------------------------------------------------
!    Compute the electric field from the potential: E = -grad Phi
!    Compute the electric flux: D = epsilon (-grad Phi)
!    Compute the total electric energy: W_e,tot = Integral (E . D)dV
!------------------------------------------------------------------------------

    IF ( CalculateField .OR. CalculateFlux .OR. CalculateEnergy ) THEN
      CALL GeneralElectricFlux( Mesh, Potential )
    END IF

    IF ( CalculateEnergy ) THEN
      Wetot = ParallelReduction(Wetot)      
      WRITE( Message, * ) 'Tot. Electric Energy  :', Wetot
      CALL Info( 'StatElecSolve', Message, Level=4 )
      CALL ListAddConstReal( Model % Simulation, &
          'RES: Electric Energy', Wetot )
    END IF

!------------------------------------------------------------------------------
!    Try to find a potential difference for scalar capacitance calculation
!------------------------------------------------------------------------------
    IF ( .NOT. CalculateCapMatrix ) THEN      
      PotentialDifference = ListGetConstReal( Params, &
           'Potential Difference',gotIt )
      IF ( .NOT.gotIt )  PotentialDifference = ListGetConstReal( &
          Model % Simulation, 'Potential Difference',gotIt )
      IF ( .NOT. gotIt ) THEN
        DO i = 1, Model % NumberOfMaterials
          PotentialDifference = &
              ListGetConstReal( Model % Materials(i) % Values, &
              'Potential Difference', GotIt )
          IF ( GotIt )  EXIT
        END DO
      END IF
      
      IF(.NOT. GotIt) THEN
        ! parallel reduction needed
        MinPotential = ParallelReduction(MinPotential,1)
        MaxPotential = ParallelReduction(MaxPotential,2)
        PotentialDifference = MaxPotential - MinPotential
      END IF
      
      IF(PotentialDifference > TINY(PotentialDifference)) THEN
        CALL ListAddConstReal( Model % Simulation, &
            'RES: Potential Difference', PotentialDifference )
        
        IF(CalculateEnergy) THEN
          Capacitance = 2*Wetot / (PotentialDifference*PotentialDifference)
          WRITE( Message,* ) 'Potential difference: ',PotentialDifference
          CALL Info( 'StatElecSolve', Message, Level=8 )
          
          WRITE( Message, * ) 'Capacitance           :', Capacitance
          CALL Info( 'StatElecSolve', Message, Level=4 )
          
          CALL ListAddConstReal( Model % Simulation, &
              'res: Capacitance', Capacitance )
        END IF
      END IF
    END IF
    
!------------------------------------------------------------------------------

    IF(CalculateCapMatrix) THEN

      PValues => Solver % Matrix % Values
      Solver % Matrix % Values => Solver % Matrix % BulkValues
      CALL MatrixVectorMultiply( Solver % Matrix, Potential, Charges)
      Solver % Matrix % Values => PValues
      Charges = Charges * PermittivityOfVacuum
      
      Permi = iter
      DO i=1,Mesh % NumberOfNodes
        Permj = CapBodyIndex(i)
        IF(Permj > 0) THEN
          j = PotentialPerm(i)
          CapMatrix(Permi, Permj) = CapMatrix(Permi, Permj) + Charges(j)
        END IF
      END DO
      DO Permj = 1, CapBodies
        IF(Permi == Permj) CYCLE
        CapMatrix(Permi, Permi) = CapMatrix(Permi, Permi) + CapMatrix(Permi, Permj)
        CapMatrix(Permi, Permj) = -CapMatrix(Permi, Permj)
      END DO

    ELSE

      WRITE( Message, * ) 'Result Norm   : ',Norm
      CALL Info( 'StatElecSolve', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'StatElecSolve', Message, Level=4 )       
      CALL Info( 'StatElecSolve', ' ', Level=4 )       
      
      IF( Solver % Variable % NonlinConverged == 1 ) EXIT
    END IF
    
  END DO
   
   IF(CalculateCapMatrix) THEN
     ! Symmetrisize
     
     IF( ParEnv % PEs > 1 ) THEN
       ALLOCATE( CapMatrixPara( CapBodies, CapBodies ) )
       CapMatrixPara = CapMatrix
       CALL MPI_ALLREDUCE(CapMatrixPara, CapMatrix, CapBodies**2, MPI_DOUBLE_PRECISION, MPI_SUM, &
           ELMER_COMM_WORLD, i)
       DEALLOCATE( CapMatrixPara )
     END IF

     IF( ParEnv % MyPE == 0 ) THEN
       CapMatrix = 0.5_dp * (CapMatrix + TRANSPOSE(CapMatrix))
       
       CALL Info('StatElecSolve','Capacitance matrix computation performed (i,j,C_ij)',Level=4)
       DO i=1, CapBodies 
         DO j = i, CapBodies
           WRITE( Message, '(I3,I3,ES15.5)' ) i,j,CapMatrix(i,j)
           CALL Info( 'StatElecSolve', Message, Level=4 )
         END DO
       END DO
      
       CapMatrixFile = ListGetString(Params,'Capacitance Matrix Filename',GotIt )
       IF( GotIt ) THEN
         OPEN (10, FILE=CapMatrixFile)
         DO i=1,CapBodies
           DO j=1,CapBodies
             WRITE (10,'(ES17.9)',advance='no') CapMatrix(i,j)
           END DO
           WRITE(10,'(A)') ' '
         END DO
         CLOSE(10)     
         WRITE(Message,'(A,A)') 'Capacitance matrix was saved to file ',CapMatrixFile
         CALL Info('StatElecSolve',Message)
       END IF
     END IF


   END IF
   

   IF ( ListGetLogical( Params, 'Adaptive Mesh Refinement', GotIt ) ) &
       CALL RefineMesh( Model, Solver, Potential, PotentialPerm, &
       ElectricInsideResidual, ElectricEdgeResidual, ElectricBoundaryResidual )
   
   CALL InvalidateVariable( Model % Meshes, Solver % Mesh, 'Potential')

   IF ( CalculateField ) THEN
     CALL InvalidateVariable( Model % Meshes, Solver % Mesh, 'Electric Field')
   END IF
   
   IF ( CalculateFlux ) THEN
     CALL InvalidateVariable( Model % Meshes, Solver % Mesh, 'Electric Flux')
   END IF
   
   IF ( CalculateEnergy ) THEN
     CALL InvalidateVariable( Model % Meshes, Solver % Mesh, 'Electric Energy Density')
   END IF

   CALL DefaultFinish()

   
!------------------------------------------------------------------------------
 
   CONTAINS

!------------------------------------------------------------------------------
     SUBROUTINE BulkAssembly()
       !------------------------------------------------------------------------------
       !$omp parallel shared(Solver, Model, dim, at0, &
       !$omp                 PermittivityOfVacuum, Message) &
       !$omp          private(t, CurrentElement, i, j, k, n, ntot, NodeIndexes, &
       !$omp                  bf_id, gotIt, Var, TID) default(none)

       TID = 1
       !$ TID = omp_get_thread_num()+1

       !$omp do 
       DO t = 1,GetNOFActive()

         !------------------------------------------------------------------------------
         !        Check if this element belongs to a body where potential
         !        should be calculated
         !------------------------------------------------------------------------------           
         CurrentElement => GetActiveElement(t)
         n = GetElementNOFNOdes(CurrentElement)
         ntot = GetElementNOFDOFs(CurrentElement)

         NodeIndexes => CurrentElement % NodeIndexes
         ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
         ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
         ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
         !------------------------------------------------------------------------------

         bf_id = ListGetInteger( Model % Bodies(CurrentElement % BodyId) % &
               Values, 'Body Force',gotIt, minv=1, maxv=Model % NumberOfBodyForces )
         Load  = 0.0_dp
         PiezoMaterial = .FALSE.
         IF ( gotIt ) THEN
           Load(1:n) = ListGetReal( Model % BodyForces(bf_id) % Values, &
                 'Charge Density', n, NodeIndexes, GotIt )        
           Load(1:n) = Load(1:n) / PermittivityOfVacuum        
           PiezoMaterial = GetLogical( Model % BodyForces(bf_id) % Values, &
                 'Piezo Material', GotIt ) 
         END IF

         k = ListGetInteger( Model % Bodies(CurrentElement % BodyId) % &
               Values, 'Material', minv=1, maxv=Model % NumberOfMaterials )

         !------------------------------------------------------------------------------
         !      Read permittivity values (might be a tensor)
         !------------------------------------------------------------------------------
         CALL ListGetRealArray( Model % Materials(k) % Values, &
               'Relative Permittivity', Pwrk,n,NodeIndexes, gotIt )
         IF ( .NOT. gotIt ) &
               CALL ListGetRealArray( Model % Materials(k) % Values, &
               'Permittivity', Pwrk, n, NodeIndexes, gotIt )

         IF ( .NOT. gotIt ) CALL Fatal( 'StatElecSolve', &
               'No > Relative permittivity < found!' )

         Permittivity = 0.0_dp
         IF ( SIZE(Pwrk,1) == 1 ) THEN
           DO i=1,3
             Permittivity( i,i,1:n ) = Pwrk( 1,1,1:n )
           END DO
         ELSE IF ( SIZE(Pwrk,2) == 1 ) THEN
           DO i=1,MIN(3,SIZE(Pwrk,1))
             Permittivity(i,i,1:n) = Pwrk(i,1,1:n)
           END DO
         ELSE
           DO i=1,MIN(3,SIZE(Pwrk,1))
             DO j=1,MIN(3,SIZE(Pwrk,2))
               Permittivity( i,j,1:n ) = Pwrk(i,j,1:n)
             END DO
           END DO
         END IF

         !------------------------------------------------------------------------------
         !      Read piezo material coefficients if applicable
         !------------------------------------------------------------------------------
         IF ( PiezoMaterial ) THEN
           PiezoCoeff = 0.0_dp
           CALL GetRealArray( Model % Materials(k) % Values, Pz_w, &
                 'Piezo Material Coefficients', gotIt, CurrentElement )
           IF ( .NOT. GotIt )  CALL Fatal( 'StatElecSolve', &
                 'No > Piezo Material Coefficients < defined!' )        
           DO i=1, Dim
             DO j=1, 2*Dim
               PiezoCoeff( i,j,1:n ) = Pz_w(i,j,1:n)
             END DO
           END DO

           !------------------------------------------------------------------------------
           !      Read also the local displacement
           !------------------------------------------------------------------------------         
           Displacement = 0.0_dp
           NULLIFY (Var)
           Var => VariableGet( Model % Variables, 'Displacement' )
           IF ( .NOT. ASSOCIATED( Var ) )  THEN
             CALL Fatal('StatElecSolve', 'No displacements' )
           END IF
           DO i = 1, Var % DOFs
             Displacement(1:n,i) = &
                   Var % Values( Var % DOFs * ( Var % Perm( NodeIndexes ) - 1 ) + i )
           END DO
         END IF

         !------------------------------------------------------------------------------
         !      Get element local matrix, and rhs vector
         !------------------------------------------------------------------------------
         CALL StatElecCompose( LocalStiffMatrix,LocalForce, PiezoMaterial, &
               PiezoCoeff, Permittivity,Load,CurrentElement,n,ntot,ElementNodes, &
               Displacement )

         !------------------------------------------------------------------------------
         !      Update global matrix and rhs vector from local matrix & vector
         !------------------------------------------------------------------------------
         CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce, UElement=CurrentElement)

         !------------------------------------------------------------------------------
         !     Print the state of the assembly to stdout
         !------------------------------------------------------------------------------
         IF ( TID == 1 .AND. RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
                 (Solver % Mesh % NumberOfBulkElements-t) / &
                 (1.0*Solver % Mesh % NumberOfBulkElements)), ' % done'          
           CALL Info( 'StatElecSolve', Message, Level=5 )          
           at0 = RealTime()
         END IF
       END DO
       !$omp end do
       !$omp end parallel

       !------------------------------------------------------------------------------
     END SUBROUTINE BulkAssembly
     !------------------------------------------------------------------------------

     !------------------------------------------------------------------------------
     SUBROUTINE BoundaryAssembly()
       !------------------------------------------------------------------------------
       !------------------------------------------------------------------------------
       !     Neumann boundary conditions
       !------------------------------------------------------------------------------
       CM => Solver % Matrix % ConstraintMatrix
       IF (ASSOCIATED(CM) ) THEN
         CM % Values = 0._dp
         IF ( .NOT. CM % Ordered ) CALL CRS_SortMatrix(CM)

         ALLOCATE( Done(Solver % Mesh % NumberOfNodes) )
         Done = .FALSE.

         Ivals => ListGetIntegerArray( Params, &
               'Constraint DOF 1 Body', Gotit )
         ! IF ( .NOT. ASSOCIATED(Ivals) ) CONTINUE ! This did not seem right..
         IF (ASSOCIATED(Ivals)) THEN

           Solver % Matrix % ConstraintMatrix % RHS(1) = &
                 GetCReal( Params, 'Constraint DOF 1 Value', GotIt )

           !$omp parallel shared(Ivals, PotentialPerm, Done, CM, Solver) &
           !$omp private(i, j, k, l, m, DoneL, CurrentElement, NodeIndexes) default(none)

           !$omp do
           DO i=1,GetNOFBoundaryElements()
             CurrentElement => GetBoundaryElement(i)
             IF ( .NOT. ActiveBoundaryElement(CurrentElement) ) CYCLE

             j = -1
             IF ( ASSOCIATED(CurrentElement % BoundaryInfo % Left) ) &
                   j = CurrentElement % BoundaryInfo % Left % BodyId

             k = -1
             IF ( ASSOCIATED(CurrentElement % BoundaryInfo % Right) ) &
                   k = CurrentElement % BoundaryInfo % Right % BodyId

             IF ( ANY(Ivals==j.OR.Ivals==k) ) THEN
               NodeIndexes => CurrentElement % NodeIndexes
               DO j=1,GetElementNOFNodes()
                 l = PotentialPerm(NodeIndexes(j))
                 
                 ! IF ( Done(l) ) CYCLE
                 ! Done(l) = .TRUE.
                                
                 ! NOTE: critical could be replace by atomic capture
                 !$omp critical(StatElecSolvePermDone)
                 DoneL = Done(l)
                 Done(l) = .TRUE.
                 !$omp end critical(StatElecSolvePermDone)
                 IF (DoneL) CYCLE

                 !$omp atomic
                 CM % RHS(1) = CM % RHS(1)+Solver % Matrix % RHS(l)
                 m = CM % Rows(1)
                 DO k=Solver % Matrix % Rows(l),Solver % Matrix % Rows(l+1)-1
                   DO WHILE(m<CM % Rows(2))
                     IF ( CM % Cols(m)>=Solver % Matrix % Cols(k) ) EXIT
                     m = m+1  
                   END DO
                   IF ( m>=CM % Rows(2) ) EXIT
                   IF ( CM % Cols(m) == Solver % Matrix % Cols(k) ) THEN
                     !$omp atomic
                     CM % Values(m) = CM % Values(m)+Solver % Matrix % Values(k)
                   END IF
                 END DO
               END DO
             END IF
           END DO
           !$omp end do
           !$omp end parallel
         END IF

         DEALLOCATE(Done)
       END IF

       MinPotential = HUGE(MinPotential)
       MaxPotential = -HUGE(MaxPotential)

       !$omp parallel shared(Solver, Mesh, PermittivityOfVacuum, VarName) &
       !$omp          private(BC, t, n, ntot, gotit, FluxBC, NodeIndexes, &
       !$omp                  LayerBC, OpenBC, CurrentElement) &
       !$omp          reduction(min:MinPotential) reduction(max:MaxPotential) default(none)

       !$omp do
       DO t= 1, Mesh % NumberOfBoundaryElements

         CurrentElement => GetBoundaryElement(t)
         IF ( .NOT. ActiveBoundaryElement(CurrentElement) ) CYCLE 
         BC => GetBC(CurrentElement)
         IF ( .NOT.ASSOCIATED( BC ) ) CYCLE 

         n = GetElementNOFNodes(CurrentElement)
         ntot = GetElementNOFDOFs(CurrentElement)
         NodeIndexes => CurrentElement % NodeIndexes

         !------------------------------------------------------------------------------
         ! Memorize the min and max potential as given by Dirichtlet BCs
         ! These are not needed here so for lower dimensional elements we may 
         ! cycle thereafter. 
         !------------------------------------------------------------------------------          

         Load(1:n) = ListGetReal( BC, &
               ComponentName(Solver % Variable), n, NodeIndexes, gotIt)
         IF(GotIt) THEN
           MinPotential = MIN(MinPotential, MINVAL(Load(1:n)))
           MaxPotential = MAX(MaxPotential, MAXVAL(Load(1:n)))             
         END IF

         ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
         ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
         ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

         !------------------------------------------------------------------------------
         !             BC: epsilon@Phi/@n = g
         !------------------------------------------------------------------------------
         Load = 0.0_dp
         Load(1:n) = ListGetReal( BC,'Electric Flux', &
               n,NodeIndexes,FluxBC )
         IF ( .NOT. FluxBC )  Load(1:n) = ListGetReal( BC, &
               'Surface Charge Density', n,NodeIndexes, FluxBC )
         IF(FluxBC) THEN
           Load(1:n) = Load(1:n) / PermittivityOfVacuum
         END IF

         !------------------------------------------------------------------------------
         !             BC: -epsilon@Phi/@n = -alpha Phi + beta
         !------------------------------------------------------------------------------          
         Alpha(1:n) = ListGetReal( BC, &
               'Layer Relative Permittivity',n, NodeIndexes,LayerBC )
         Beta(1:n) = 0.0_dp

         OpenBC = ListGetLogical( BC,'Electric Infinity BC',GotIt)
         IF(.NOT. GotIt) OpenBC = ListGetLogical( BC,'Infinity BC '//TRIM(VarName),GotIt)

         IF(.NOT. ( LayerBC .OR. FluxBC .OR. OpenBC) ) CYCLE

         IF ( LayerBC ) THEN
           LayerH(1:n) = ListGetReal( BC, &
                 'Layer Thickness', n, NodeIndexes, gotit )
           IF ( .NOT. gotit ) THEN
             CALL Fatal( 'StatElecSolve','Charge > Layer thickness < not given!' )
           END IF
           Alpha(1:n) = Alpha(1:n) / LayerH(1:n)

           LayerV(1:n) = ListGetReal( BC, &
                 'Electrode Potential', n, NodeIndexes, gotit )
           Beta(1:n) = ListGetReal( BC, &
                 'Layer Charge Density', n, NodeIndexes, gotit )
           Beta(1:n) = Alpha(1:n)*LayerV(1:n) + 0.5_dp*Beta(1:n)*LayerH(1:n) / PermittivityOfVacuum            
         END IF

         !------------------------------------------------------------------------------
         !             BC: -epsilon@Phi/@n = epsilon*Phi*(r \cdot n)/r^2
         !------------------------------------------------------------------------------         
         IF( OpenBC ) THEN
           PermIso(1:n) = GetParentMatProp('Relative Permittivity',&
                 CurrentElement,GotIt)
           IF(.NOT. GotIt) THEN
             PermIso(1:n) = GetParentMatProp('Relative Permittivity',&
                   CurrentElement,GotIt)
           END IF
           IF(.NOT. GotIt) THEN
             CALL Fatal( 'StatElecSolve','Could not find > Relative Permittivity < for parent!' )           
           END IF
         END IF

         !------------------------------------------------------------------------------
         !             Get element matrix and rhs due to boundary conditions ...
         !------------------------------------------------------------------------------
         CALL StatElecBoundary( LocalStiffMatrix, LocalForce,  &
               Load, Alpha, Beta, OpenBC, PermIso, CurrentElement, &
               n, ElementNodes )

         !------------------------------------------------------------------------------
         !             Update global matrices from local matrices
         !------------------------------------------------------------------------------
         CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce, UElement=CurrentElement )

         !------------------------------------------------------------------------------   
       END DO   ! Neumann BCs
       !------------------------------------------------------------------------------
       !$omp end do
       !$omp end parallel

       !------------------------------------------------------------------------------
       !    FinishAssembly must be called after all other assembly steps, but before
       !    Dirichlet boundary settings. Actually no need to call it except for
       !    transient simulations.
       !------------------------------------------------------------------------------
       CALL DefaultFinishAssembly()

       !------------------------------------------------------------------------------
       !   This sets the BC flags so that the potential form a permulation
       !------------------------------------------------------------------------------
       IF(CalculateCapMatrix) THEN
         CALL SetPermutationBoundaries( Model, StiffMatrix, ForceVector, &
               'Capacitance Body',PotentialPerm, iter)
       END IF

       !------------------------------------------------------------------------------
       !    Dirichlet boundary conditions
       !------------------------------------------------------------------------------
       CALL DefaultDirichletBCs()

       at = CPUTime() - at
       WRITE( Message, * ) 'Assembly (s)          :',at
       CALL Info( 'StatElecSolve', Message, Level=4 )
       !------------------------------------------------------------------------------
     END SUBROUTINE BoundaryAssembly
     !------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE TotalChargeBC(F,Element,n,Nodes)
!------------------------------------------------------------------------------
     TYPE(Nodes_t) :: Nodes
     INTEGER :: n
     REAL(KIND=dp) :: F(:)
     TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
     TYPE(GaussIntegrationPoints_t) :: IntegStuff
     INTEGER :: i, j, pn
     LOGICAL :: stat
     TYPE(Nodes_t), SAVE :: Pnodes
     TYPE(Element_t), POINTER ::  Parent
     REAL(KIND=dp) :: s,u,v,w,detJ,pdetJ,pBasis(10),Basis(10),dBasisdx(10,3),Normal(3)
     !$omp threadprivate(Pnodes)

     Parent => Element % BoundaryInfo % Left
     CALL GetElementNodes( PNodes, Parent )
     pn = Parent % TYPE % NumberOfNodes

     IntegStuff = GaussPoints(Element,RelOrder = RelIntegOrder)

     F = 0._dp
     DO i=1,IntegStuff % N
       u = IntegStuff % u(i)
       v = IntegStuff % v(i)
       w = IntegStuff % w(i)
       Normal = NormalVector(Element,Nodes,u,v,.TRUE.)
       stat = ElementInfo(Element,Nodes,u,v,w,detJ,Basis )
       CALL GetParentUVW( Element,n,Parent,pn,u,v,w,Basis )
       stat = ElementInfo(Parent,PNodes,u,v,w,pdetJ,Basis,dBasisdx )
       DO j=1,3
         F(1:pn) = F(1:pn) - IntegStuff % s(i) * detJ * &
               dBasisdx(1:pn,j)  * Normal(j)
       END DO
     END DO
     
!------------------------------------------------------------------------------
   END SUBROUTINE TotalChargeBC
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Compute the Electric Flux, Electric Field and Electric Energy at model nodes.
!------------------------------------------------------------------------------
   SUBROUTINE GeneralElectricFlux( Mesh, Potential )
!------------------------------------------------------------------------------
     TYPE(Mesh_t) :: Mesh
     REAL(KIND=dp) :: Potential(:)
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element, Parent
     TYPE(Nodes_t) :: Nodes, BoundaryNodes
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     REAL(KIND=dp), POINTER :: U_Integ(:), V_Integ(:), W_Integ(:), S_Integ(:)
     REAL(KIND=dp), ALLOCATABLE :: SumOfWeights(:), SurfWeights(:), x(:), y(:), z(:)
     REAL(KIND=dp) :: PermittivityOfVacuum
     ! REAL(KIND=dp) :: Permittivity(3,3,Mesh % MaxElementNodes)
     ! REAL(KIND=dp) :: Basis(Mesh % MaxElementDofs)
     ! REAL(KIND=dp) :: dBasisdx(Mesh % MaxElementDofs,3)
     REAL(KIND=DP) :: SqrtElementMetric, detJ
     REAL(KIND=dp), ALLOCATABLE :: ElementPot(:)
     REAL(KIND=dp) :: EnergyDensity, Sigma, Normal(3)
     REAL(KIND=dp) :: NodalFlux(3), NodalField(3), ElemVol
     REAL(KIND=dp) :: s, ug, vg, wg, Grad(3), EpsGrad(3)
     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
     REAL(KIND=dp) :: xpos, ypos, zpos
     INTEGER, ALLOCATABLE :: Indexes(:), PotIndexes(:)
     INTEGER :: n, N_Integ, t, tg, i, j, k, DIM, p, nPar, matId, nd, istat
     LOGICAL :: Stat

!------------------------------------------------------------------------------
     IF(CalculateEnergy) Energy = 0.0_dp
     IF(CalculateFlux) Flux = 0.0_dp
     IF(CalculateField) Field = 0.0_dp
     
     DIM = CoordinateSystemDimension()
     Wetot = 0.0_dp

     ALLOCATE(SumOfWeights(SIZE(Potential)), STAT=istat)
     IF (istat /= 0) CALL Fatal('GeneralElectricFlux',&
                                'Memory allocation failed 1')
     SumOfWeights = 0.0_dp
     
     PermittivityOfVacuum = ListGetConstReal( Model % Constants, &
           'Permittivity Of Vacuum',gotIt )
     IF ( .NOT.gotIt ) PermittivityOfVacuum = 1

     !$omp parallel shared(Mesh, Model, Solver, DIM, Energy, Flux, &
     !$omp                 Field, PermittivityOfVacuum, PotentialPerm, &
     !$omp                 Potential, ConstantWeights, CalculateEnergy, &
     !$omp                 CalculateFlux, CalculateField, SumOfWeights) &
     !$omp          private(tg, t, i, j, k, n, s, nd, ug, vg, wg, Element, &
     !$omp                  Indexes, NodeIndexes, xpos, ypos, zpos, &
     !$omp                  PotIndexes, ElementPot, Nodes, &
     !$omp                  IntegStuff, U_Integ, V_Integ, Symb, dSymb, &
     !$omp                  SqrtMetric, Metric, &
     !$omp                  W_Integ, S_Integ, N_Integ, EnergyDensity, &
     !$omp                  NodalFlux, NodalField, ElemVol, &
     !$omp                  SqrtElementMetric, &
     !$omp                  Grad, EpsGrad, istat, gotIt, stat) &
     !$omp                  reduction(+:Wetot) default(none)
     
     n = Mesh % MaxElementNodes
     Permittivity = 0_dp
     ! Allocate thread local workspace
     ALLOCATE(Nodes % x(n), Nodes % y(n), Nodes % z(n), &
              Indexes(Mesh % MaxElementDofs), &
              PotIndexes(Mesh % MaxElementDofs), &
              ElementPot(Mesh % MaxElementDofs), STAT=istat)
     IF (istat /= 0) CALL Fatal('GeneralElectricFlux',&
                                'Memory allocation failed 2')

!------------------------------------------------------------------------------
!   Go through model elements, we will compute on average of elementwise
!   fluxes to nodes of the model
!------------------------------------------------------------------------------

     !$omp do
     DO t=1,Solver % NumberOfActiveElements

!------------------------------------------------------------------------------
!        Check if this element belongs to a body where electrostatics
!        should be calculated
!------------------------------------------------------------------------------

       Element => GetActiveElement(t)
       NodeIndexes => Element % NodeIndexes
       n = GetElementNOFNOdes(Element)
       nd = GetElementDOFs(Indexes, Element)

       PotIndexes(1:nd) = PotentialPerm( Indexes(1:nd) )
       ElementPot(1:nd) = Potential( PotIndexes(1:nd) )

       Nodes % x(1:n) = Mesh % Nodes % x( NodeIndexes )
       Nodes % y(1:n) = Mesh % Nodes % y( NodeIndexes )
       Nodes % z(1:n) = Mesh % Nodes % z( NodeIndexes )

!------------------------------------------------------------------------------
!    Gauss integration stuff
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints( Element )
       U_Integ => IntegStuff % u
       V_Integ => IntegStuff % v
       W_Integ => IntegStuff % w
       S_Integ => IntegStuff % s
       N_Integ =  IntegStuff % n
!------------------------------------------------------------------------------

       k = ListGetInteger( Model % Bodies( Element % BodyId ) % &
             Values, 'Material', minv=1, maxv=Model % NumberOfMaterials )

       CALL ListGetRealArray( Model % Materials(k) % Values, &
             'Relative Permittivity', Pwrk, n, NodeIndexes, gotIt )
       IF ( .NOT. gotIt ) &
             CALL ListGetRealArray( Model % Materials(k) % Values, &
             'Permittivity', Pwrk, n, NodeIndexes, gotIt )

       Permittivity = 0.0_dp
       IF ( SIZE(Pwrk,1) == 1 ) THEN
         DO i=1,3
           Permittivity( i,i,1:n ) = Pwrk( 1,1,1:n )
         END DO
       ELSE IF ( SIZE(Pwrk,2) == 1 ) THEN
         DO i=1,MIN(3,SIZE(Pwrk,1))
           Permittivity(i,i,1:n) = Pwrk(i,1,1:n)
         END DO
       ELSE
         DO i=1,MIN(3,SIZE(Pwrk,1))
           DO j=1,MIN(3,SIZE(Pwrk,2))
             Permittivity( i,j,1:n ) = Pwrk(i,j,1:n)
           END DO
         END DO
       END IF

       EnergyDensity = 0.0_dp
       NodalFlux = 0.0_dp
       NodalField = 0.0_dp
       ElemVol = 0.0_dp

!------------------------------------------------------------------------------
! Loop over Gauss integration points
!------------------------------------------------------------------------------
       DO tg=1,N_Integ

         ug = U_Integ(tg)
         vg = V_Integ(tg)
         wg = W_Integ(tg)

!------------------------------------------------------------------------------
! Need SqrtElementMetric and Basis at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Element, Nodes,ug,vg,wg, &
               SqrtElementMetric,Basis,dBasisdx )

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
         IF (CurrentCoordinateSystem() /= Cartesian ) THEN
           xpos = SUM( Nodes % x(1:n) * Basis(1:n) )
           ypos = SUM( Nodes % y(1:n) * Basis(1:n) )
           zpos = SUM( Nodes % z(1:n) * Basis(1:n) )
         END IF

         CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,xpos,ypos,zpos )

         s = SqrtMetric * SqrtElementMetric * S_Integ(tg)

         !------------------------------------------------------------------------------

         EpsGrad = 0.0_dp
         DO j=1, DIM
           Grad(j) = SUM( dBasisdx(1:nd,j) * ElementPot(1:nd) )
           DO i = 1, DIM
             EpsGrad(j) = EpsGrad(j) + SUM( Permittivity(j,i,1:n) * &
                   Basis(1:n) ) * SUM( dBasisdx(1:nd,i) * ElementPot(1:nd) )
           END DO
         END DO

         Wetot = Wetot + s * SUM( Grad(1:DIM) * EpsGrad(1:DIM) )

         EnergyDensity = EnergyDensity + &
               s * SUM(Grad(1:DIM) * EpsGrad(1:DIM))

         DO j = 1,DIM
           NodalFlux(j) = NodalFlux(j) - EpsGrad(j) * s
           NodalField(j) = NodalField(j) - Grad(j) * s
         END DO

         ElemVol = ElemVol + s
       END DO ! Gauss point integration

!------------------------------------------------------------------------------
!   Weight with element area if required
!------------------------------------------------------------------------------
       IF ( ConstantWeights ) THEN
         EnergyDensity = EnergyDensity / ElemVol
         NodalFlux(1:DIM) = NodalFlux(1:DIM) / ElemVol
         NodalField(1:DIM) = NodalField(1:DIM) / ElemVol
         DO j=1,nd
           !$omp atomic
           SumOfWeights( PotIndexes(j) ) = &
                 SumOfWeights( PotIndexes(j) ) + 1
         END DO
       ELSE
         DO j=1,nd
           !$omp atomic
           SumOfWeights( PotIndexes(j) )  = &
                 SumOfWeights( PotIndexes(j) ) + ElemVol
         END DO
       END IF

!------------------------------------------------------------------------------

       IF(CalculateEnergy) THEN
         DO j=1,nd
           !$omp atomic
           Energy( PotIndexes(j) ) = Energy( PotIndexes(j) ) + EnergyDensity
         END DO
       END IF

       IF(CalculateFlux) THEN
         NodalFlux = NodalFlux * PermittivityOfVacuum
         DO i=1,DIM
           DO j=1,nd
             !$omp atomic
             Flux( DIM * ( PotIndexes(j)-1) + i ) = &
                   Flux( DIM * ( PotIndexes(j)-1) + i ) + NodalFlux(i)
           END DO
         END DO
       END IF

       IF(CalculateField) THEN
         DO i=1,DIM
           DO j=1,nd
             !$omp atomic
             Field( DIM * ( PotIndexes(j)-1) + i ) = &
                   Field( DIM * ( PotIndexes(j)-1) + i ) + NodalField(i)
           END DO
         END DO
       END IF

     END DO ! element loop
     !$omp end do
     

!------------------------------------------------------------------------------
!   Finally, compute average of the fluxes at nodes
!------------------------------------------------------------------------------
     !$omp do 
     DO j = 1, SIZE( Potential )
       IF ( ABS( SumOfWeights(j) ) > 0.0_dp ) THEN
         IF ( CalculateEnergy )  Energy(j) = Energy(j) / SumOfWeights( j )

         IF ( CalculateField ) THEN
           DO k=1,DIM
             Field( DIM*(j-1)+k ) = Field( DIM*(j-1)+k) / SumOfWeights( j )
           END DO
         END IF

         IF ( CalculateFlux ) THEN
           DO k=1,DIM
             Flux( DIM*(j-1)+k ) = Flux( DIM*(j-1)+k) / SumOfWeights( j )
           END DO
         END IF
       END IF
     END DO
     !$omp end do

     IF(CalculateEnergy) THEN
       !$omp do
       DO j=1,SIZE(Potential)
         Energy(j) = PermittivityOfVacuum * Energy(j) / 2.0_dp
       END DO
       !$omp end do
     END IF

     DEALLOCATE( Nodes % x, Nodes % y, Nodes % z, &
                 Indexes, PotIndexes, ElementPot )
     !$omp end parallel

     DEALLOCATE(SumOfWeights)

     ! Compute total energy
     Wetot = PermittivityOfVacuum * Wetot / 2.0_dp
!------------------------------------------------------------------------------
   END SUBROUTINE GeneralElectricFlux
!------------------------------------------------------------------------------

 
!------------------------------------------------------------------------------
     SUBROUTINE StatElecCompose( StiffMatrix,Force,PiezoMaterial, PiezoCoeff, &
                            Permittivity,Load,Element,n,ntot,Nodes, Displacement )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: StiffMatrix(:,:),Force(:),Load(:), Permittivity(:,:,:)
       REAL(KIND=dp) :: PiezoCoeff(:,:,:), Displacement(:,:)
       INTEGER :: n
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
       LOGICAL :: PiezoMaterial
!------------------------------------------------------------------------------
 
       REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)
       ! REAL(KIND=dp) :: Basis(ntot),dBasisdx(ntot,3)
       REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,A,L,C(3,3),x,y,z
       REAL(KIND=dp) :: PiezoForce(ntot), LocalStrain(6), PiezoLoad(3)

       LOGICAL :: Stat

       INTEGER :: i,j,p,q,t,DIM,ntot,Nbasis
 
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------

       DIM = CoordinateSystemDimension()
       NBasis = ntot

       PiezoForce = 0.0_dp
       Force = 0.0_dp
       StiffMatrix = 0.0_dp
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints( Element, RelOrder = RelIntegOrder )

       DO t=1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
         S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Element,Nodes,U,V,W,SqrtElementMetric, &
                    Basis,dBasisdx )
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
         IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
           x = SUM( ElementNodes % x(1:n) * Basis(1:n) )
           y = SUM( ElementNodes % y(1:n) * Basis(1:n) )
           z = SUM( ElementNodes % z(1:n) * Basis(1:n) )
         END IF

         CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
 
         S = S * SqrtElementMetric * SqrtMetric
!------------------------------------------------------------------------------
!        The piezo force term
!------------------------------------------------------------------------------

         IF ( PiezoMaterial ) THEN
           ! So far only plane strain in 2D  (LocalStrain(3) = 0)           
           LocalStrain = 0.0_dp
           DO i = 1, Dim
             LocalStrain(i) = SUM( dBasisdx(1:n,i) * Displacement(1:n,i) )
           END DO
           LocalStrain(4) = 0.5_dp * ( SUM( dBasisdx(1:n,1) * Displacement(1:n,2) ) &
               + SUM( dBasisdx(1:n,2) * Displacement(1:n,1) ) )
           IF ( Dim == 3 ) THEN
             LocalStrain(5) = 0.5_dp * ( SUM( dBasisdx(1:n,2) * Displacement(1:n,3) ) &
                 + SUM( dBasisdx(1:n,3) * Displacement(1:n,2) ) )
             LocalStrain(6) = 0.5_dp * ( SUM( dBasisdx(1:n,1) * Displacement(1:n,3) ) &
                 + SUM( dBasisdx(1:n,3) * Displacement(1:n,1) ) )
           END IF
           
           PiezoLoad = 0.0_dp
           DO i = 1, Dim
             DO j = 1, 2*Dim
               PiezoLoad(i) = PiezoLoad(i) + SUM( Basis(1:n) * PiezoCoeff(i,j,1:n) ) * &
                   LocalStrain(j)
             END DO
           END DO
         END IF
         
!------------------------------------------------------------------------------
         L = SUM( Load(1:n) * Basis )
         DO i=1,DIM
           DO j=1,DIM
             C(i,j) = SUM( Permittivity(i,j,1:n) * Basis(1:n) )
           END DO
         END DO
!------------------------------------------------------------------------------
!        The Poisson equation
!------------------------------------------------------------------------------
         DO p=1,Nbasis
           DO q=1,Nbasis
             A = 0._dp
             DO i=1,DIM
               DO J=1,DIM
                 A = A + C(i,j) * dBasisdx(p,i) * dBasisdx(q,j)
               END DO
             END DO
             StiffMatrix(p,q) = StiffMatrix(p,q) + S*A
           END DO
           Force(p) = Force(p) + S*L*Basis(p)

           IF ( PiezoMaterial ) THEN
             PiezoForce(p) = PiezoForce(p) + S * SUM( dBasisdx(p,1:Dim) * PiezoLoad(1:Dim) )
           END IF

        END DO
!------------------------------------------------------------------------------
       END DO

       IF ( PiezoMaterial )  Force = Force + PiezoForce
!       IF ( PiezoMaterial )  Force = Force + PiezoForce / PermittivityOfVacuum ?

!------------------------------------------------------------------------------
     END SUBROUTINE StatElecCompose
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> To compute the capacitance matrix for m bodies, also m permutations are 
!> required. This subroutine sets at permutation i, the body i to 1, 
!> and all else to 0.
!------------------------------------------------------------------------------
  SUBROUTINE SetPermutationBoundaries( Model, StiffMatrix, ForceVector, &
      Name, Perm, Permutation)
!------------------------------------------------------------------------------

    TYPE(Model_t) :: Model
    TYPE(Matrix_t), POINTER :: StiffMatrix
    REAL(KIND=dp) :: ForceVector(:)
    CHARACTER(LEN=*) :: Name 
    INTEGER :: Perm(:), Permutation
!------------------------------------------------------------------------------

    TYPE(Element_t), POINTER :: CurrentElement
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: i,j,k,n,t,k1,k2, Body, MaxBody
    LOGICAL :: GotIt
    REAL(KIND=dp) :: val,s

!------------------------------------------------------------------------------
!  Manipulates the list structure of the BCs so that the permutations BCs
!  are translated into corresponding Dirichlet conditions.
!------------------------------------------------------------------------------
    MaxBody = 0
    DO i=1,Model % NumberOfBCs
      BC => Model % BCs(i) % Values
      Body = ListGetInteger( BC, Name, gotIt )
      
      IF ( gotIt ) THEN           
        IF(Body == Permutation) THEN
          val = 1.0_dp
        ELSE
          val = 0.0_dp
        END IF
        
        CALL ListAddConstReal( BC,TRIM(VarName),val)
        MaxBody = MAX( MaxBody, Body ) 
      END IF
    END DO

    IF( Permutation == 1 ) THEN
      DO t = Mesh % NumberOfBulkElements + 1, &
          Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
        
        CurrentElement => Mesh % Elements(t)
!------------------------------------------------------------------------------
!      Set the current element pointer in the model structure to
!      reflect the element being processed
!------------------------------------------------------------------------------
        Model % CurrentElement => Mesh % Elements(t)
!------------------------------------------------------------------------------
        n = CurrentElement % TYPE % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes
        
        DO i=1,Model % NumberOfBCs
          IF ( CurrentElement % BoundaryInfo % Constraint == &
              Model % BCs(i) % Tag ) THEN
            
            BC => Model % BCs(i) % Values
            Body = ListGetInteger( BC, Name, gotIt )
            IF( GotIt ) CapBodyIndex(NodeIndexes) = Body
          END IF
        END DO
      END DO

      PRINT *,'Number of permutation BCs'
      DO i=1,MaxBody
        PRINT *,'Capacitance body:',i,'no',COUNT(CapBodyIndex == i )
      END DO
    END IF

!------------------------------------------------------------------------------
   END SUBROUTINE SetPermutationBoundaries
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE StatElecBoundary( BoundaryMatrix, BoundaryVector, &
        LoadVector, Alpha, Beta, OpenBC, Permittivity, Element, n, Nodes )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: BoundaryMatrix(:,:) 
     REAL(KIND=dp) :: BoundaryVector(:)  
     REAL(KIND=dp) :: LoadVector(:) 
     REAL(KIND=dp) :: Alpha(:)  !< Coefficient of the Robin BC: g = alpha * u + beta
     REAL(KIND=dp) :: Beta(:)   !< Coefficient of the Robin BC: g = alpha * u + beta
     LOGICAL :: OpenBC
     REAL(KIND=dp) :: Permittivity(:)
     TYPE(Element_t), POINTER :: Element
     TYPE(Element_t), POINTER :: PElement
     INTEGER :: n
     TYPE(Nodes_t)   :: Nodes
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric
     REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3),&
         Coord(3), Normal(3)
     REAL(KIND=dp) :: u,v,w,s,x,y,z
     REAL(KIND=dp) :: Force, AlphaAtIP, BetaAtIP, PermAtIP
     REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)
     INTEGER :: t,p,q,N_Integ
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     LOGICAL :: stat
!------------------------------------------------------------------------------

     BoundaryVector = 0.0_dp
     BoundaryMatrix = 0.0_dp
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( Element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivates at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                 Basis,dBasisdx )

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
      IF ( CurrentCoordinateSystem() /= Cartesian .OR. OpenBC ) THEN
        x = SUM( ElementNodes % x(1:n)*Basis(1:n) )
        y = SUM( ElementNodes % y(1:n)*Basis(1:n) )
        z = SUM( ElementNodes % z(1:n)*Basis(1:n) )
      END IF
      
      CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
      
      s = S_Integ(t) * SqrtElementMetric * SqrtMetric
      
      AlphaAtIP = SUM( Basis(1:n) * Alpha(1:n))
      BetaAtIp = SUM( Basis(1:n) * Beta(1:n))

      IF( OpenBC ) THEN
        PermAtIP = SUM( Basis(1:n) * Permittivity(1:n) ) 
        Normal = NormalVector( Element, Nodes, u, v, .TRUE. )
        Coord(1) = x
        Coord(2) = y
        Coord(3) = z
        AlphaAtIP = PermAtIP * SUM( Coord * Normal ) / SUM( Coord * Coord ) 
      END IF

!------------------------------------------------------------------------------
      Force = SUM( LoadVector(1:n) * Basis )
      Force = Force + BetaAtIp
      DO p=1,N
        DO q=1,N
          BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + &
              s * AlphaAtIp * Basis(q) * Basis(p)
        END DO
      END DO
      
      DO q=1,N
        BoundaryVector(q) = BoundaryVector(q) + s * Basis(q) * Force
      END DO
     END DO
   END SUBROUTINE StatElecBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 END SUBROUTINE StatElecSolver
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION ElectricBoundaryResidual( Model, Edge, Mesh, Quant, Perm, Gnorm ) &
       RESULT( Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
     TYPE( Mesh_t ), POINTER    :: Mesh
     TYPE( Element_t ), POINTER :: Edge
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element


     INTEGER :: i,j,k,n,l,t,DIM,Pn,En
     LOGICAL :: stat, Found

     REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
     REAL(KIND=dp) :: Grad(3,3), Normal(3), EdgeLength
     REAL(KIND=dp) :: u, v, w, s, detJ
     REAL(KIND=dp) :: Residual, ResidualNorm, Permittivity

     REAL(KIND=dp), ALLOCATABLE :: EdgeBasis(:), Basis(:)
     REAL(KIND=dp), ALLOCATABLE :: Flux(:)
     REAL(KIND=dp), ALLOCATABLE :: x(:), y(:), z(:), dBasisdx(:,:)
     REAL(KIND=dp), ALLOCATABLE :: NodalPermittivity(:), Potential(:)

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: First = .TRUE., Dirichlet

     SAVE Hwrk, First
     !$omp threadprivate(First, Hwrk)
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     IF ( First ) THEN
        First = .FALSE.
        NULLIFY( Hwrk )
     END IF

     Gnorm     = 0.0_dp

     Metric = 0.0_dp
     DO i=1,3
        Metric(i,i) = 1.0_dp
     END DO

     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           DIM = 3
        CASE DEFAULT
           DIM = CoordinateSystemDimension()
     END SELECT
!    
!    ---------------------------------------------

     Element => Edge % BoundaryInfo % Left

     IF ( .NOT. ASSOCIATED( Element ) ) THEN

        Element => Edge % BoundaryInfo % Right

     ELSE IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) THEN

        Element => Edge % BoundaryInfo % Right

     END IF

     IF ( .NOT. ASSOCIATED( Element ) ) RETURN
     IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN

     En = Edge % TYPE % NumberOfNodes
     Pn = Element % TYPE % NumberOfNodes

     ALLOCATE( EdgeNodes % x(En), EdgeNodes % y(En), EdgeNodes % z(En) )

     EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
     EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
     EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)

     ALLOCATE( Nodes % x(Pn), Nodes % y(Pn), Nodes % z(Pn) )

     Nodes % x = Mesh % Nodes % x(Element % NodeIndexes)
     Nodes % y = Mesh % Nodes % y(Element % NodeIndexes)
     Nodes % z = Mesh % Nodes % z(Element % NodeIndexes)

     ALLOCATE( EdgeBasis(En), Basis(Pn), dBasisdx(Pn,3), Flux(En), &
      x(En), y(En), z(En), NodalPermittivity(En), Potential(Pn) )

     DO l = 1,En
       DO k = 1,Pn
          IF ( Edge % NodeIndexes(l) == Element % NodeIndexes(k) ) THEN
             x(l) = Element % TYPE % NodeU(k)
             y(l) = Element % TYPE % NodeV(k)
             z(l) = Element % TYPE % NodeW(k)
             EXIT
          END IF
       END DO
     END DO
!
!    Integrate square of residual over boundary element:
!    ---------------------------------------------------

     Indicator    = 0.0_dp
     EdgeLength   = 0.0_dp
     ResidualNorm = 0.0_dp

     DO j=1,Model % NumberOfBCs
        IF ( Edge % BoundaryInfo % Constraint /= Model % BCs(j) % Tag ) CYCLE

!       IF ( .NOT. ListGetLogical( Model % BCs(j) % Values, &
!                 'Heat Flux BC', Found ) ) CYCLE

!
!       Check if dirichlet BC given:
!       ----------------------------
        s = ListGetConstReal( Model % BCs(j) % Values, &
              ComponentName(Model % Solver % Variable), Dirichlet )

!       Get various flux bc options:
!       ----------------------------

!       ...given flux:
!       --------------
        Flux(1:En) = ListGetReal( Model % BCs(j) % Values, &
          'Electric Flux', En, Edge % NodeIndexes, Found )


!       get material parameters:
!       ------------------------
        k = ListGetInteger(Model % Bodies(Element % BodyId) % Values,'Material', &
                    minv=1, maxv=Model % NumberOFMaterials)

        CALL ListGetRealArray( Model % Materials(k) % Values, &
               'Relative Permittivity', Hwrk, En, Edge % NodeIndexes,stat )
        IF ( .NOT. stat )  &
             CALL ListGetRealArray( Model % Materials(k) % Values, &
             'Permittivity', Hwrk, En, Edge % NodeIndexes )

        NodalPermittivity( 1:En ) = Hwrk( 1,1,1:En )

!       elementwise nodal solution:
!       ---------------------------
        Potential(1:Pn) = Quant( Perm(Element % NodeIndexes) )

!       do the integration:
!       -------------------
        EdgeLength   = 0.0_dp
        ResidualNorm = 0.0_dp

        IntegStuff = GaussPoints( Edge )

        DO t=1,IntegStuff % n
           u = IntegStuff % u(t)
           v = IntegStuff % v(t)
           w = IntegStuff % w(t)

           stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
               EdgeBasis, dBasisdx )

           Normal = NormalVector( Edge, EdgeNodes, u, v, .TRUE. )

           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              s = IntegStuff % s(t) * detJ
           ELSE
              u = SUM( EdgeBasis(1:En) * EdgeNodes % x(1:En) )
              v = SUM( EdgeBasis(1:En) * EdgeNodes % y(1:En) )
              w = SUM( EdgeBasis(1:En) * EdgeNodes % z(1:En) )
      
              CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                         Symb, dSymb, u, v, w )

              s = IntegStuff % s(t) * detJ * SqrtMetric
           END IF

!
!          Integration point in parent element local
!          coordinates:
!          -----------------------------------------
           u = SUM( EdgeBasis(1:En) * x(1:En) )
           v = SUM( EdgeBasis(1:En) * y(1:En) )
           w = SUM( EdgeBasis(1:En) * z(1:En) )

           stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
                 Basis, dBasisdx )
!
!          Heat conductivity at the integration point:
!          --------------------------------------------
           Permittivity = SUM( NodalPermittivity(1:En) * EdgeBasis(1:En) )
!
!          given flux at integration point:
!          --------------------------------
           Residual = -SUM( Flux(1:En) * EdgeBasis(1:En) )


!          flux given by the computed solution, and 
!          force norm for scaling the residual:
!          -----------------------------------------
           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              DO k=1,DIM
                 Residual = Residual + Permittivity  * &
                    SUM( dBasisdx(1:Pn,k) * Potential(1:Pn) ) * Normal(k)

                 Gnorm = Gnorm + s * (Permittivity * &
                       SUM(dBasisdx(1:Pn,k) * Potential(1:Pn)) * Normal(k))**2
              END DO
           ELSE
              DO k=1,DIM
                 DO l=1,DIM
                    Residual = Residual + Metric(k,l) * Permittivity  * &
                       SUM( dBasisdx(1:Pn,k) * Potential(1:Pn) ) * Normal(l)

                    Gnorm = Gnorm + s * (Metric(k,l) * Permittivity * &
                      SUM(dBasisdx(1:Pn,k) * Potential(1:Pn) ) * Normal(l))**2
                 END DO
              END DO
           END IF

           EdgeLength   = EdgeLength + s
           IF ( .NOT. Dirichlet ) THEN
              ResidualNorm = ResidualNorm + s * Residual ** 2
           END IF
        END DO
        EXIT
     END DO

     IF ( CoordinateSystemDimension() == 3 ) THEN
        EdgeLength = SQRT(EdgeLength)
     END IF

!    Gnorm = EdgeLength * Gnorm
     Indicator = EdgeLength * ResidualNorm

     DEALLOCATE( Nodes % x, Nodes % y, Nodes % z)
     DEALLOCATE( EdgeNodes % x, EdgeNodes % y, EdgeNodes % z)

     DEALLOCATE( EdgeBasis, Basis, dBasisdx, Flux, x, y, z, &
             NodalPermittivity, Potential )
!------------------------------------------------------------------------------
  END FUNCTION ElectricBoundaryResidual
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  FUNCTION ElectricEdgeResidual( Model, Edge, Mesh, Quant, Perm ) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE

     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2)
     TYPE( Mesh_t ), POINTER    :: Mesh
     TYPE( Element_t ), POINTER :: Edge
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element

     INTEGER :: i,j,k,l,n,t,DIM,En,Pn
     LOGICAL :: stat, Found

     REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
     REAL(KIND=dp) :: Permittivity
     REAL(KIND=dp) :: u, v, w, s, detJ
     REAL(KIND=dp) :: Grad(3,3), Normal(3), EdgeLength, Jump

     REAL(KIND=dp), ALLOCATABLE :: NodalPermittivity(:)
     REAL(KIND=dp), ALLOCATABLE :: x(:),y(:),z(:)
     REAL(KIND=dp), ALLOCATABLE :: EdgeBasis(:)
     REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), Potential(:)

     REAL(KIND=dp) :: Residual, ResidualNorm

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: First = .TRUE.

     SAVE Hwrk, First
     !$omp threadprivate(First, Hwrk)
!------------------------------------------------------------------------------

!    Initialize:
!    -----------

     IF ( First ) THEN
        First = .FALSE.
        NULLIFY( Hwrk )
     END IF

     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           DIM = 3
        CASE DEFAULT
           DIM = CoordinateSystemDimension()
     END SELECT

     Metric = 0.0_dp
     DO i = 1,3
        Metric(i,i) = 1.0_dp
     END DO

     Grad = 0.0_dp
!
!    ---------------------------------------------

     Element => Edge % BoundaryInfo % Left
     n = Element % TYPE % NumberOfNodes

     Element => Edge % BoundaryInfo % Right
     n = MAX( n, Element % TYPE % NumberOfNodes )

     ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )

     En = Edge % TYPE % NumberOfNodes
     ALLOCATE( EdgeNodes % x(En), EdgeNodes % y(En), EdgeNodes % z(En) )

     EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
     EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
     EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)

     ALLOCATE( x(En), y(En), z(En), NodalPermittivity(En), EdgeBasis(En), &
             Basis(n), dBasisdx(n,3), Potential(n) )

!    Integrate square of jump over edge:
!    -----------------------------------
     ResidualNorm = 0.0_dp
     EdgeLength   = 0.0_dp
     Indicator    = 0.0_dp

     IntegStuff = GaussPoints( Edge )

     DO t=1,IntegStuff % n

        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
             EdgeBasis, dBasisdx )

        Normal = NormalVector( Edge, EdgeNodes, u, v, .FALSE. )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM( EdgeBasis(1:En) * EdgeNodes % x(1:En) )
           v = SUM( EdgeBasis(1:En) * EdgeNodes % y(1:En) )
           w = SUM( EdgeBasis(1:En) * EdgeNodes % z(1:En) )

           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                      Symb, dSymb, u, v, w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        ! 
        ! Compute flux over the edge as seen by elements
        ! on both sides of the edge:
        ! ----------------------------------------------
        DO i = 1,2
           SELECT CASE(i)
              CASE(1)
                 Element => Edge % BoundaryInfo % Left
              CASE(2)
                 Element => Edge % BoundaryInfo % Right
           END SELECT
!
!          Can this really happen (maybe it can...)  ?      
!          -------------------------------------------
           IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) CYCLE
!
!          Next, get the integration point in parent
!          local coordinates:
!          -----------------------------------------
           Pn = Element % TYPE % NumberOfNodes

           DO j = 1,En
              DO k = 1,Pn
                 IF ( Edge % NodeIndexes(j) == Element % NodeIndexes(k) ) THEN
                    x(j) = Element % TYPE % NodeU(k)
                    y(j) = Element % TYPE % NodeV(k)
                    z(j) = Element % TYPE % NodeW(k)
                    EXIT
                 END IF
              END DO
           END DO

           u = SUM( EdgeBasis(1:En) * x(1:En) )
           v = SUM( EdgeBasis(1:En) * y(1:En) )
           w = SUM( EdgeBasis(1:En) * z(1:En) )
!
!          Get parent element basis & derivatives at the integration point:
!          -----------------------------------------------------------------
           Nodes % x(1:Pn) = Mesh % Nodes % x(Element % NodeIndexes)
           Nodes % y(1:Pn) = Mesh % Nodes % y(Element % NodeIndexes)
           Nodes % z(1:Pn) = Mesh % Nodes % z(Element % NodeIndexes)

           stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
             Basis, dBasisdx )
!
!          Material parameters:
!          --------------------
           k = ListGetInteger( Model % Bodies( &
                    Element % BodyId) % Values, 'Material', &
                     minv=1, maxv=Model % NumberOFMaterials )

           CALL ListGetRealArray( Model % Materials(k) % Values, &
                'Relative Permittivity', Hwrk, En, Edge % NodeIndexes,stat )
           IF ( .NOT. stat )  &
                CALL ListGetRealArray( Model % Materials(k) % Values, &
                'Permittivity', Hwrk, En, Edge % NodeIndexes )

           NodalPermittivity( 1:En ) = Hwrk( 1,1,1:En )
           Permittivity = SUM( NodalPermittivity(1:En) * EdgeBasis(1:En) )
!
!          Potential at element nodal points:
!          ------------------------------------
           Potential(1:Pn) = Quant( Perm(Element % NodeIndexes) )
!
!          Finally, the flux:
!          ------------------
           DO j=1,DIM
              Grad(j,i) = Permittivity * SUM( dBasisdx(1:Pn,j) * Potential(1:Pn) )
           END DO
        END DO

!       Compute square of the flux jump:
!       -------------------------------   
        EdgeLength  = EdgeLength + s
        Jump = 0.0_dp
        DO k=1,DIM
           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              Jump = Jump + (Grad(k,1) - Grad(k,2)) * Normal(k)
           ELSE
              DO l=1,DIM
                 Jump = Jump + &
                       Metric(k,l) * (Grad(k,1) - Grad(k,2)) * Normal(l)
              END DO
           END IF
        END DO
        ResidualNorm = ResidualNorm + s * Jump ** 2
     END DO

     IF ( CoordinateSystemDimension() == 3 ) THEN
        EdgeLength = SQRT(EdgeLength)
     END IF
     Indicator = EdgeLength * ResidualNorm

     DEALLOCATE( Nodes % x, Nodes % y, Nodes % z)
     DEALLOCATE( EdgeNodes % x, EdgeNodes % y, EdgeNodes % z)

     DEALLOCATE( x, y, z, NodalPermittivity, EdgeBasis, Basis, &
                dBasisdx, Potential )
!------------------------------------------------------------------------------
  END FUNCTION ElectricEdgeResidual
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   FUNCTION ElectricInsideResidual( Model, Element, Mesh, &
        Quant, Perm, Fnorm ) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
     TYPE( Mesh_t ), POINTER    :: Mesh
     TYPE( Element_t ), POINTER :: Element
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes

     LOGICAL :: stat, Found
     INTEGER :: i,j,k,l,n,t,DIM

     REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
     REAL(KIND=dp) :: Permittivity
     REAL(KIND=dp) :: u, v, w, s, detJ
     REAL(KIND=dp) :: Source, Residual, ResidualNorm, Area

     REAL(KIND=dp), ALLOCATABLE :: NodalPermittivity(:)
     REAL(KIND=dp), ALLOCATABLE :: PrevPot(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalSource(:), Potential(:)
     REAL(KIND=dp), ALLOCATABLE :: Basis(:)
     REAL(KIND=dp), ALLOCATABLE :: dBasisdx(:,:), ddBasisddx(:,:,:)

     TYPE( ValueList_t ), POINTER :: Material

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: First = .TRUE.

     SAVE Hwrk, First
     !$omp threadprivate(First, Hwrk)
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     Indicator = 0.0_dp
     Fnorm     = 0.0_dp
!
!    Check if this eq. computed in this element:
!    -------------------------------------------
     IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN

     IF ( First ) THEN
        First = .FALSE.
        NULLIFY( Hwrk )
     END IF

     Metric = 0.0_dp
     DO i=1,3
        Metric(i,i) = 1.0_dp
     END DO

     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           DIM = 3
        CASE DEFAULT
           DIM = CoordinateSystemDimension()
     END SELECT
!
!    Element nodal points:
!    ---------------------
     n = Element % TYPE % NumberOfNodes

     ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n), &
       NodalPermittivity(n), Basis(n), dBasisdx(n,3),    &
       ddBasisddx(n,3,3), PrevPot(n), NodalSource(n), Potential(n) )

     Nodes % x = Mesh % Nodes % x(Element % NodeIndexes)
     Nodes % y = Mesh % Nodes % y(Element % NodeIndexes)
     Nodes % z = Mesh % Nodes % z(Element % NodeIndexes)
!
!    Elementwise nodal solution:
!    ---------------------------
     Potential(1:n) = Quant( Perm(Element % NodeIndexes) )
!
!    Material parameters: relative permittivity
!    ------------------------------------------
     k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Material', &
                     minv=1, maxv=Model % NumberOfMaterials )

     Material => Model % Materials(k) % Values

     CALL ListGetRealArray( Model % Materials(k) % Values, &
          'Relative Permittivity', Hwrk, n, Element % NodeIndexes,stat )
     IF ( .NOT. stat )  &
          CALL ListGetRealArray( Model % Materials(k) % Values, &
          'Permittivity', Hwrk, n, Element % NodeIndexes )

     NodalPermittivity( 1:n ) = Hwrk( 1,1,1:n )

!
!    Charge density (source):
!    ------------------------
!
     k = ListGetInteger( &
         Model % Bodies(Element % BodyId) % Values,'Body Force',Found, &
                 1, Model % NumberOFBodyForces)

     NodalSource = 0.0_dp
     IF ( Found .AND. k > 0  ) THEN
        NodalSource(1:n) = ListGetReal( Model % BodyForces(k) % Values, &
             'Charge Density', n, Element % NodeIndexes, stat )
        IF ( .NOT. stat )  &
             NodalSource(1:n) = ListGetReal( Model % BodyForces(k) % Values, &
             'Source', n, Element % NodeIndexes )
     END IF
!
!    Integrate square of residual over element:
!    ------------------------------------------

     ResidualNorm = 0.0_dp
     Area = 0.0_dp

     IntegStuff = GaussPoints( Element )

     DO t=1,IntegStuff % n
        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
            Basis, dBasisdx, ddBasisddx, .TRUE., .FALSE. )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM( Basis(1:n) * Nodes % x(1:n) )
           v = SUM( Basis(1:n) * Nodes % y(1:n) )
           w = SUM( Basis(1:n) * Nodes % z(1:n) )

           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                       Symb, dSymb, u, v, w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        Permittivity = SUM( NodalPermittivity(1:n) * Basis(1:n) )
!
!       Residual of the electrostatic equation:
!
!        R = -div(e grad(u)) - s
!       ---------------------------------------------------
!
!       or more generally:
!
!        R = -g^{jk} (C T_{,j}}_{,k} - s
!       ---------------------------------------------------
!
        Residual = -SUM( NodalSource(1:n) * Basis(1:n) )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           DO j=1,DIM
!
!             - grad(e).grad(T):
!             --------------------
!
              Residual = Residual - &
                 SUM( Potential(1:n) * dBasisdx(1:n,j) ) * &
                 SUM( NodalPermittivity(1:n) * dBasisdx(1:n,j) )

!
!             - e div(grad(u)):
!             -------------------
!
              Residual = Residual - Permittivity * &
                 SUM( Potential(1:n) * ddBasisddx(1:n,j,j) )
           END DO
        ELSE
           DO j=1,DIM
              DO k=1,DIM
!
!                - g^{jk} C_{,k}T_{j}:
!                ---------------------
!
                 Residual = Residual - Metric(j,k) * &
                    SUM( Potential(1:n) * dBasisdx(1:n,j) ) * &
                    SUM( NodalPermittivity(1:n) * dBasisdx(1:n,k) )

!
!                - g^{jk} C T_{,jk}:
!                -------------------
!
                 Residual = Residual - Metric(j,k) * Permittivity * &
                    SUM( Potential(1:n) * ddBasisddx(1:n,j,k) )
!
!                + g^{jk} C {_jk^l} T_{,l}:
!                ---------------------------
                 DO l=1,DIM
                    Residual = Residual + Metric(j,k) * Permittivity * &
                      Symb(j,k,l) * SUM( Potential(1:n) * dBasisdx(1:n,l) )
                 END DO
              END DO
           END DO
        END IF

!
!       Compute also force norm for scaling the residual:
!       -------------------------------------------------
        DO i=1,DIM
           Fnorm = Fnorm + s * ( SUM( NodalSource(1:n) * Basis(1:n) ) ) ** 2
        END DO

        Area = Area + s
        ResidualNorm = ResidualNorm + s *  Residual ** 2
     END DO

!    Fnorm = Element % hk**2 * Fnorm
     Indicator = Element % hK**2 * ResidualNorm

     DEALLOCATE( Nodes % x, Nodes % y, Nodes % z, NodalPermittivity, &
        Basis, dBasisdx, ddBasisddx, PrevPot, NodalSource, Potential )
!------------------------------------------------------------------------------
  END FUNCTION ElectricInsideResidual
!------------------------------------------------------------------------------

