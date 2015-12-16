!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Ville Savolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \}


!------------------------------------------------------------------------------
!>  Utility routines for the elmer main program.
!------------------------------------------------------------------------------
MODULE MainUtils

!------------------------------------------------------------------------------
  Use BlockSolve
#ifdef USE_ISO_C_BINDINGS
  USE LoadMod
#endif

!------------------------------------------------------------------------------
  IMPLICIT NONE
!------------------------------------------------------------------------------
    
    LOGICAL, PRIVATE :: isParallel=.FALSE.

CONTAINS

!------------------------------------------------------------------------------
!> Routine checks the feasibility of solver options.
!------------------------------------------------------------------------------
  SUBROUTINE CheckLinearSolverOptions( Solver ) 
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: Found
    CHARACTER(LEN=MAX_NAME_LEN) :: str
!------------------------------------------------------------------------------

    Params => ListGetSolverParams(Solver)
    str = ListGetString( Params,'Linear System Solver', Found )

    IF ( str == 'direct' ) THEN
      str = ListGetString( Params,'Linear System Direct Method', Found )
      
      IF( Found ) THEN        
        IF ( ParEnv % PEs > 1 ) THEN
          IF ( str /= 'mumps' .AND. str /= 'cpardiso' ) THEN
            CALL Warn( 'CheckLinearSolverOptions', 'Only MUMPS and CPardiso direct solver' // &
                ' interface implemented in parallel, trying MUMPS!')
            str = 'mumps' 
            CALL ListAddString( Params,'Linear System Direct Method', str)
          END IF
        ELSE
!         IF ( str == 'mumps' ) THEN
!           CALL Warn( 'CheckSolverOptions', 'Currently no serial interface' // &
!               ' to the MUMPS solver implemented, trying UMFPACK!')
!           str = 'umfpack'    
!           CALL ListAddString( Params,'Linear System Direct Method', str)
!         END IF
        END IF
        
        SELECT CASE( str )
        CASE('banded' )
          
        CASE( 'umfpack', 'big umfpack' )
#include "../config.h"
#ifndef HAVE_UMFPACK
          CALL Fatal( 'GetMatrixFormat', 'UMFPACK solver has not been installed.' )
#endif
        CASE( 'mumps', 'mumpslocal' )
#ifndef HAVE_MUMPS
          CALL Fatal( 'CheckLinearSolverOptions', 'MUMPS solver has not been installed.' )
#endif
        CASE( 'superlu' )
#ifndef HAVE_SUPERLU
          CALL Fatal( 'CheckLinearSolverOptions', 'SuperLU solver has not been installed.' )
#endif
        CASE( 'pardiso' )
#if !defined(HAVE_PARDISO) && !defined(HAVE_MKL)
          CALL Fatal( 'CheckLinearSolverOptions', 'Pardiso solver has not been installed.' )
#endif
        CASE( 'cpardiso')
#if !defined(HAVE_CPARDISO) || !defined(HAVE_MKL)
        CALL Fatal( 'CheckSolverOptions', ' Cluster Pardiso solver has not been installed.' )
#endif
        CASE( 'cholmod','spqr' )
#ifndef HAVE_CHOLMOD
          CALL Fatal( 'CheckLinearSolverOptions', 'Cholmod solver has not been installed.' )
#endif
        CASE DEFAULT
          CALL Fatal( 'CheckLinearSolverOptions', 'Unknown direct solver method: ' // TRIM(str) )
        END SELECT
        
      ELSE
        IF ( ParEnv % PEs <= 1 ) THEN
#ifdef HAVE_UMFPACK
          str = 'umfpack'
#else
          str = 'banded'
#endif
        ELSE
#ifdef HAVE_MUMPS
          str = 'mumps'
#else
          CALL Fatal( 'CheckLinearSolverOptions', 'There is no direct parallel solver available (MUMPS)')
#endif
        END IF
        CALL Info('CheckLinearSolverOptions','Setting > Linear System Direct Method < to:'//TRIM(str) )
        CALL ListAddString( Params,'Linear System Direct Method', str )
      END IF

    ELSE IF ( str == 'feti' ) THEN
      IF( ParEnv % PEs <= 1 ) THEN
        CALL Fatal('CheckLinearSolverOptions','Feti not usable in serial!')
      END IF

    ELSE
      IF (ListGetLogical( Params,  &
          'Linear System Use Hypre', Found )) THEN
        IF( ParEnv % PEs <= 1 ) THEN
          CALL Fatal('CheckLinearSolverOptions','Hypre not usable in serial!')
        END IF
#ifndef HAVE_HYPRE
        CALL Fatal('CheckLinearSolverOptions','Hypre requested but not compiled with!')
#endif
      END IF

      IF (ListGetLogical( Params,  &
          'Linear System Use Trilinos', Found )) THEN        
        IF( ParEnv % PEs <= 1 ) THEN
          CALL Fatal('CheckLinearSolverOptions','Trilinos not usable in serial!')
        END IF
#ifndef HAVE_TRILINOS
        CALL Fatal('CheckLinearSolverOptions','Trilinos requested but not compiled with!')
#endif
      END IF
    END IF
  
!------------------------------------------------------------------------------
  END SUBROUTINE CheckLinearSolverOptions
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! This subroutine enables the scaling of keywords by geometric entities. 
! If any Real type keyword has an associated keyword with suffix 'Normalize by Area'
! or 'Normalize by Volume' the keyword will be divided with the the area/volume.
! This way any quantity may be made distributed on-the-fly. The current implementation
! assumes fixed geometry. The inverse of area/volume will be located at ValueList % Coeff.
! There is some small additional cost for fetching the value. 
!------------------------------------------------------------------------------
   SUBROUTINE SetNormalizedKeywords(Model, Mesh)
     TYPE(Model_t) :: Model
     TYPE(Mesh_t), POINTER :: Mesh

     LOGICAL :: AnyBC, AnyBodyForce, AnyMat, AnyBody
     CHARACTER(LEN=MAX_NAME_LEN) :: Suffix
     INTEGER :: list_ind
     REAL(KIND=dp) :: EntityWeight
     LOGICAL :: Found
     TYPE(ValueList_t), POINTER :: List

     ! If the weights have been already computed for linear cases no need to redo
     IF( Mesh % EntityWeightsComputed ) THEN
       IF( .NOT. ListGetLogical( Model % Simulation,&
           'Update Keyword Normalization',Found ) ) RETURN
     END IF


     Suffix = 'Normalize By Area'
     AnyBC = ListCheckSuffixAnyBC( Model, Suffix )

     Suffix = 'Normalize By Volume'
     AnyBodyForce = ListCheckSuffixAnyBodyForce( Model, Suffix )
     AnyMat = ListCheckSuffixAnyMaterial( Model, Suffix )
     AnyBody = ListCheckSuffixAnyBody( Model, Suffix )

     IF( .NOT. ( AnyBC .OR. AnyBodyForce .OR. AnyMat .OR. AnyBody ) ) RETURN

     ! This computes all the weights as the cost is not that large and
     ! places the 'Entity Weight' to each list. 
     CALL CalculateEntityWeights( Model, Mesh ) 

     Suffix = 'Normalize By Area'
     IF( AnyBC ) THEN
       Found = .FALSE.
       DO list_ind = 1,Model % NumberOfBCs
         List => Model % BCs(list_ind) % Values
         EntityWeight = Mesh % BCWeight(list_ind)
         CALL ListSetCoefficients( List, Suffix,1.0_dp / EntityWeight )
       END DO
     END IF

     Suffix = 'Normalize By Volume'
     IF( AnyBodyForce ) THEN
       Found = .FALSE.
       DO list_ind = 1,Model % NumberOfBodyForces
         List => Model % BodyForces(list_ind) % Values
         EntityWeight = Mesh % BodyForceWeight(list_ind)
         CALL ListSetCoefficients( List, Suffix,1.0_dp / EntityWeight )
       END DO
     END IF

     IF( AnyMat ) THEN
       Found = .FALSE.
       DO list_ind = 1,Model % NumberOfMaterials
         List => Model % Materials(list_ind) % Values
         EntityWeight = Mesh % MaterialWeight(list_ind)
         CALL ListSetCoefficients( List, Suffix,1.0_dp / EntityWeight )
       END DO
     END IF

     IF( AnyBody ) THEN
       Found = .FALSE.
       DO list_ind = 1,Model % NumberOfBodies
         List => Model % Bodies(list_ind) % Values
         EntityWeight = Mesh % BodyWeight(list_ind)
         CALL ListSetCoefficients( List, Suffix,1.0_dp / EntityWeight )
       END DO
     END IF

   END SUBROUTINE SetNormalizedKeywords
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! This subroutine enables the rotation of vector valued keywords.
!------------------------------------------------------------------------------
   SUBROUTINE SetRotatedProperties(Model, Mesh)
     TYPE(Model_t) :: Model
     TYPE(Mesh_t), POINTER :: Mesh

     LOGICAL :: AnyBC, AnyBodyForce, AnyMat, AnyBody
     CHARACTER(LEN=MAX_NAME_LEN) :: Suffix
     INTEGER :: list_ind
     REAL(KIND=dp) :: RotMatrix(3,3), Angles(3,1)
     REAL(KIND=dp), POINTER :: Pmatrix(:,:)
     INTEGER, TARGET :: NodeIndex(1)
     INTEGER, POINTER :: NodeIndexes(:)
     REAL(KIND=dp) :: EntityWeight
     LOGICAL :: Found
     TYPE(ValueList_t), POINTER :: List
     LOGICAL :: Visited = .FALSE.

     SAVE Visited, Suffix, AnyBC, AnyBodyForce, AnyMat, AnyBody

     IF(.NOT. Visited ) THEN
       Suffix = 'Property Rotate'
       AnyBC = ListCheckSuffixAnyBC( Model, Suffix )
       AnyBodyForce = ListCheckSuffixAnyBodyForce( Model, Suffix )
       AnyMat = ListCheckSuffixAnyMaterial( Model, Suffix )
       AnyBody = ListCheckSuffixAnyBody( Model, Suffix )
       Visited = .TRUE.
     END IF
       
     IF( .NOT. ( AnyBC .OR. AnyBodyForce .OR. AnyMat .OR. AnyBody ) ) RETURN

     NodeIndex(1) = 1
     NodeIndexes => NodeIndex

     IF( AnyBody ) THEN
       DO list_ind = 1,Model % NumberOfBodies
         List => Model % Bodies(list_ind) % Values                  
         IF(.NOT. ListCheckSuffix( List, Suffix ) ) CYCLE

         CALL ListGetRealVector( List,'Property Rotate',Angles,1,NodeIndexes,Found )
         IF(.NOT. Found ) THEN
           CALL Fatal('SetRotatedProperties','> Property Rotate < suffix requires angles!')
         END IF
         
         RotMatrix = AnglesToRotationMatrix( Angles(1:3,1) )

         PMatrix => ListGetConstRealArray( List,'Property Rotation Matrix',Found )
         IF( Found ) THEN
           PMatrix = RotMatrix(3,3) 
         ELSE
           CALL ListAddConstRealArray( List,'Property Rotation Matrix', 3, 3, RotMatrix )
         END IF
       END DO
     END IF


     IF( AnyBodyForce ) THEN
       DO list_ind = 1,Model % NumberOfBodyForces
         List => Model % BodyForces(list_ind) % Values                  
         IF(.NOT. ListCheckSuffix( List, Suffix ) ) CYCLE

         CALL ListGetRealVector( List,'Property Rotate',Angles,1,NodeIndexes,Found )
         IF(.NOT. Found ) THEN
           CALL Fatal('SetRotatedProperties','> Property Rotate < suffix requires angles!')
         END IF
         
         RotMatrix = AnglesToRotationMatrix( Angles(1:3,1) )

         PMatrix => ListGetConstRealArray( List,'Property Rotation Matrix',Found )
         IF( Found ) THEN
           PMatrix = RotMatrix(3,3) 
         ELSE
           CALL ListAddConstRealArray( List,'Property Rotation Matrix', 3, 3, RotMatrix )
         END IF
       END DO
     END IF


     IF( AnyMat ) THEN
       DO list_ind = 1,Model % NumberOfMaterials
         List => Model % Materials(list_ind) % Values                  
         IF(.NOT. ListCheckSuffix( List, Suffix ) ) CYCLE

         CALL ListGetRealVector( List,'Property Rotate',Angles,1,NodeIndexes,Found )
         IF(.NOT. Found ) THEN
           CALL Fatal('SetRotatedProperties','> Property Rotate < suffix requires angles!')
         END IF
         
         RotMatrix = AnglesToRotationMatrix( Angles(1:3,1) )

         PMatrix => ListGetConstRealArray( List,'Property Rotation Matrix',Found )
         IF( Found ) THEN
           PMatrix = RotMatrix(3,3) 
         ELSE
           CALL ListAddConstRealArray( List,'Property Rotation Matrix', 3, 3, RotMatrix )
         END IF
       END DO
     END IF


     IF( AnyBC ) THEN
       DO list_ind = 1,Model % NumberOfBCs
         List => Model % BCs(list_ind) % Values                  
         IF(.NOT. ListCheckSuffix( List, Suffix ) ) CYCLE

         CALL ListGetRealVector( List,'Property Rotate',Angles,1,NodeIndexes,Found )
         IF(.NOT. Found ) THEN
           CALL Fatal('SetRotatedProperties','> Property Rotate < suffix requires angles!')
         END IF
         
         RotMatrix = AnglesToRotationMatrix( Angles(1:3,1) )

         PMatrix => ListGetConstRealArray( List,'Property Rotation Matrix',Found )
         IF( Found ) THEN
           PMatrix = RotMatrix(3,3) 
         ELSE
           CALL ListAddConstRealArray( List,'Property Rotation Matrix', 3, 3, RotMatrix )
         END IF
       END DO
     END IF

   CONTAINS 

     FUNCTION AnglesToRotationMatrix( Angles, RotateOrder ) RESULT ( RotMatrix )

       REAL(KIND=dp) :: Angles(3), RotMatrix(3,3)
       INTEGER, OPTIONAL :: RotateOrder(3)

       INTEGER :: i,j
       REAL(KIND=dp) :: TrfMatrix(3,3), IdentityMatrix(3,3), Alpha

       IdentityMatrix = 0.0_dp
       DO i=1,3
         IdentityMatrix(i,i) = 1.0_dp
       END DO

       RotMatrix = IdentityMatrix

       DO i=1,3
         j = i
         IF( PRESENT( RotateOrder ) ) j = RotateOrder(i)
         Alpha = Angles(j) * PI / 180.0_dp

         IF( ABS(Alpha) < TINY(Alpha) ) CYCLE
         TrfMatrix = IdentityMatrix

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

     END FUNCTION AnglesToRotationMatrix


   END SUBROUTINE SetRotatedProperties
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!> Get calling address of the procedure and add it to the Solver structure.
!------------------------------------------------------------------------------
  SUBROUTINE AddSolverProcedure( Solver,PROCEDURE  )
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    EXTERNAL :: PROCEDURE
    INTEGER  :: PROCEDURE
!------------------------------------------------------------------------------
#ifndef USE_ISO_C_BINDINGS
   INTEGER(KIND=AddrInt) :: AddrFunc
#else
    INTEGER(KIND=AddrInt) :: AddrFunc
    EXTERNAL :: AddrFunc
#endif
!------------------------------------------------------------------------------
    Solver % PROCEDURE = AddrFunc( PROCEDURE )
!------------------------------------------------------------------------------
  END SUBROUTINE AddSolverProcedure
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
   SUBROUTINE SwapMesh(Model,Mesh,Name)
!------------------------------------------------------------------------------
     CHARACTER(LEN=*) :: Name
     TYPE(Model_t) :: Model
     TYPE(Mesh_t), POINTER :: Mesh, Newmesh, tmesh
!------------------------------------------------------------------------------
     INTEGER :: Def_Dofs(10,6), i,j,k
     LOGICAL :: Found, Transient
     TYPE(Solver_t), POINTER :: Solver
!------------------------------------------------------------------------------

     Def_Dofs = -1;
     DO i=1,Model % NumberOfSolvers
       DO j=1,10
         DO k=1,6
           Def_Dofs(j,k) = MAX(Def_Dofs(j,k), MAXVAL(Model % Solvers(i) % Def_Dofs(j,:,k)))
         END DO
       END DO
     END DO

     Newmesh => LoadMesh2( Model, Name, Name, &
       .FALSE., Parenv % PEs, ParEnv % myPE, Def_Dofs )
     IF(.NOT.ASSOCIATED(NewMesh)) RETURN

     NewMesh % Next => Mesh % Next
     IF(ASSOCIATED(Mesh,  Model % Meshes)) THEN
       Model % Meshes => Newmesh
     ELSE
       Tmesh => Model % Meshes
       DO WHILE(ASSOCIATED(Tmesh % next))
         IF(ASSOCIATED(Mesh, Tmesh % next)) THEN
           Tmesh % Next => Newmesh
           EXIT
         END IF
       END DO
     END IF

     NewMesh % Name = Name
     CALL AddCoordAndTime(Mesh,NewMesh)

     DO i=1,Model % NumberOfSolvers
       Solver => Model % Solvers(i)
       IF(ASSOCIATED(Solver % Mesh, Mesh)) Solver % Mesh => Newmesh
     END DO

     IF(Mesh % DiscontMesh) CALL CreateDiscontMesh(Model,Newmesh,.TRUE.)

     IF(ASSOCIATED(Model % Mesh, Mesh)) Model % Mesh => NewMesh
     IF(ASSOCIATED(Model % Variables, Mesh % Variables)) Model % Variables => NewMesh % Variables

     Mesh % Next => Null()
     CALL ReleaseMesh( Mesh )

     Transient = ListGetString( Model % Simulation, 'Simulation Type' ) == 'transient'

     DO i=1,Model % NumberOfSolvers
       Solver => Model % Solvers(i)
       IF(ASSOCIATED(Solver % Mesh, NewMesh)) THEN
         CALL FreeMatrix(Solver % Matrix)
         Model % Solver => Solver

         CALL AddEquationBasics( Solver, ListGetString(Solver % Values, &
                  'Variable', Found), Transient )
         CALL AddEquationSolution( Solver, Transient )
         IF ( Transient .AND. Solver % PROCEDURE /= 0 ) CALL InitializeTimestep(Solver)
       END IF
     END DO

     CALL MeshStabParams( Newmesh )
     NewMesh % Changed = .TRUE.

CONTAINS


   SUBROUTINE AddCoordAndTime(M1,M2)
    TYPE(Solver_t), POINTER :: Solver => Null()
    TYPE(Mesh_t) :: M1,M2
    TYPE(Variable_t), POINTER :: DtVar, V

     CALL VariableAdd( M2 % Variables, M2,Solver, &
           'Coordinate 1',1, M2 % Nodes % x )

     CALL VariableAdd(M2 % Variables,M2,Solver, &
           'Coordinate 2',1, M2 % Nodes % y )

     CALL VariableAdd(M2 % Variables,M2,Solver, &
          'Coordinate 3',1,M2 % Nodes % z )

     V => VariableGet( M1 % Variables, 'Time' )
     CALL VariableAdd( M2 % Variables, M2, Solver, 'Time', 1, V % Values )

     V => VariableGet( M1 % Variables, 'Periodic Time' )
     CALL VariableAdd( M2 % Variables, M2, Solver, 'Periodic Time', 1, V % Values)

     V => VariableGet( M1 % Variables, 'Timestep' )
     CALL VariableAdd( M2 % Variables, M2, Solver, 'Timestep', 1, V % Values )

     V => VariableGet( M1 % Variables, 'Timestep size' )
     CALL VariableAdd( M2 % Variables, M2, Solver, 'Timestep size', 1, V % Values )

     V => VariableGet( M1 % Variables, 'Timestep interval' )
     CALL VariableAdd( M2 % Variables, M2, Solver, 'Timestep interval', 1, V % Values )

     ! Save some previous timesteps for variable timestep multistep methods
     V => VariableGet( M1 % Variables, 'Timestep size' )
     DtVar => VariableGet( M2 % Variables, 'Timestep size' )
     DtVar % PrevValues => V % PrevValues

     V => VariableGet( M1 % Variables, 'nonlin iter' )
     CALL VariableAdd( M2 % Variables, M2, Solver, &
             'nonlin iter', 1, V % Values )

     V => VariableGet( M1 % Variables, 'coupled iter' )
     CALL VariableAdd( M2 % Variables, M2, Solver, &
             'coupled iter', 1, V % Values )

     V => VariableGet( M1 % Variables, 'partition' )
     CALL VariableAdd( M2 % Variables, M2, Solver, 'Partition', 1, V % Values )
!------------------------------------------------------------------------------
    END SUBROUTINE AddCoordAndTime
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   END SUBROUTINE SwapMesh
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!> Add the generic stuff related to each Solver. 
!> A few solvers are for historical reasons given a special treatment. 
!------------------------------------------------------------------------------
  SUBROUTINE AddEquationBasics( Solver, Name,Transient )
!------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: Solver
    LOGICAL :: Transient
    CHARACTER(LEN=*) :: Name
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: Solution(:)

    INTEGER, POINTER :: Perm(:)

    INTEGER(KIND=AddrInt) :: InitProc

    INTEGER :: MaxDGDOFs, MaxNDOFs, MaxEDOFs, MaxFDOFs, MaxBDOFs
    INTEGER :: i,j,k,l,NDeg,Nrows,nSize,n,m,DOFs,dim,MatrixFormat,istat,Maxdim

    LOGICAL :: Found, Stat, BandwidthOptimize, EigAnal, ComplexFlag, &
    MultigridActive, VariableOutput, GlobalBubbles, HarmonicAnal, MGAlgebraic, &
    VariableGlobal, NoMatrix, IsAssemblySolver, IsCoupledSolver, IsBlockSolver, &
    IsProcedure, LegacySolver

    CHARACTER(LEN=MAX_NAME_LEN) :: str,eq,var_name,proc_name,tmpname

    TYPE(ValueList_t), POINTER :: SolverParams
    TYPE(Mesh_t),   POINTER :: NewMesh,OldMesh
    TYPE(Element_t), POINTER :: CurrentElement
    TYPE(Matrix_t), POINTER :: NewMatrix, Oldmatrix

    TYPE(Variable_t), POINTER :: Var
    TYPE(Variable_t), POINTER :: NewVariable

    REAL(KIND=dp) :: tt, InitValue
    REAL(KIND=dp), POINTER :: Component(:)


    ! Set pointer to the list of solver parameters
    !------------------------------------------------------------------------------
    SolverParams => ListGetSolverParams(Solver)

    !------------------------------------------------------------------------------
    ! Check the historical solvers that may be built-in on some .sif files
    ! Therefore some special strategies are used for them.
    !------------------------------------------------------------------------------
    IsProcedure = ListCheckPresent( SolverParams, 'Procedure' )
    Dim = CoordinateSystemDimension()        
    Dofs = 1
    InitValue = 0.0_dp
    LegacySolver = .TRUE.

    SELECT CASE( Name )       
      !------------------------------------------------------------------------------
      
      !------------------------------------------------------------------------------
    CASE('navier-stokes')
      !------------------------------------------------------------------------------
      dofs = dim 
      IF ( CurrentCoordinateSystem() == CylindricSymmetric ) DOFs = DOFs + 1
      IF( dofs == 3 ) THEN
        var_name = 'Flow Solution[Velocity:3 Pressure:1]'
      ELSE
        var_name = 'Flow Solution[Velocity:2 Pressure:1]'
      END IF
      proc_name = 'FlowSolve FlowSolver'
      InitValue = 1.0d-6
      ! We don't want to use automated setting of dofs later so set this back to one!
      dofs = 1

      !------------------------------------------------------------------------------
    CASE('magnetic induction')
      !------------------------------------------------------------------------------
      var_name = 'Magnetic Field'
      proc_name = 'MagneticSolve MagneticSolver'
      dofs = 3
      CALL ListAddString( SolverParams,&
          NextFreeKeyword('Exported Variable',SolverParams),&
          'Electric Current[Electric Current:3]')                  
      
      !------------------------------------------------------------------------------
    CASE('stress analysis')
      !------------------------------------------------------------------------------
      dofs = dim
      var_name = 'Displacement'
      proc_name = 'StressSolve StressSolver'      
            
      !------------------------------------------------------------------------------
    CASE('mesh update')
      !------------------------------------------------------------------------------
      dofs = dim
      var_name = 'Mesh Update'
      proc_name = 'MeshSolve MeshSolver'

      IF( Transient ) THEN
        IF( Dofs == 2 ) THEN
          CALL ListAddString( SolverParams,&
              NextFreeKeyword('Exported Variable',SolverParams),&
              '-dofs 2 Mesh Velocity')        
        ELSE
          CALL ListAddString( SolverParams,&
              NextFreeKeyword('Exported Variable',SolverParams),&
              '-dofs 3 Mesh Velocity')                  
        END IF
      END IF

      !------------------------------------------------------------------------------
    CASE('heat equation')
      !------------------------------------------------------------------------------
      var_name = 'Temperature'
      proc_name = 'HeatSolve HeatSolver'
      
      IF( .NOT. ListCheckPresent( SolverParams,'Radiation Solver') ) THEN
        CALL ListAddLogical( SolverParams,'Radiation Solver',.TRUE.)
      END IF
      !------------------------------------------------------------------------------

    CASE DEFAULT
      LegacySolver = .FALSE.
      
    END SELECT

    IF( LegacySolver ) THEN
      CALL Info('AddEquationBasics','Setting up keywords internally for legacy solver: '&
          //TRIM(Name),Level=10)
      IF( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
        CALL ListAddString( SolverParams,'Variable',var_name )
        IF( Dofs > 1 ) CALL ListAddInteger( SolverParams,'Variable Dofs',dofs )
      END IF
      IF( .NOT. IsProcedure ) THEN      
        CALL ListAddString(SolverParams, 'Procedure', proc_name,.FALSE.)
      END IF
    END IF

    ! We should have the procedure 
    proc_name = ListGetString( SolverParams, 'Procedure',IsProcedure)
    IF( IsProcedure ) THEN
      CALL Info('AddEquationBasics','Using procedure: '//TRIM(proc_name),Level=10)
    END IF


    ! If there is a matrix level Flux Corrected Transport and/or nonlinear timestepping
    ! then you must use global matrices for time integration.
    !----------------------------------------------------------------------------------
    IF( ListGetLogical( SolverParams,'Linear System FCT',Found ) ) THEN
      IF( ParEnv % PEs > 1 ) THEN
        CALL Fatal('AddEquationBasics','FCT scheme not implemented in parallel yet!')
      END IF
      CALL ListAddLogical( SolverParams,'Use Global Mass Matrix',.TRUE.)
    END IF
    IF( ListGetLogical( SolverParams,'Nonlinear Timestepping',Found ) ) THEN
      CALL ListAddLogical( SolverParams,'Use Global Mass Matrix',.TRUE.)
    END IF

    ! Compute the mesh dimension for this solver
    !----------------------------------------------------------------------------
    eq = ListGetString( SolverParams, 'Equation', Found )
    IF( Found ) THEN
      CALL Info('AddEquationBasics','Setting up solver: '//TRIM(eq),Level=8)
    END IF

    IF ( Found ) THEN
      MAXdim = 0
      DO i=1,Solver % Mesh % NumberOfBulkElements+Solver % Mesh % NumberOFBoundaryElements
        CurrentElement => Solver % Mesh % Elements(i)
        IF ( CheckElementEquation( CurrentModel, CurrentElement, eq ) ) THEN
          Maxdim = MAX( CurrentElement % TYPE % DIMENSION, Maxdim )
        END IF
      END DO
      CALL ListAddInteger( SolverParams, 'Active Mesh Dimension', Maxdim )
    END IF

    ! Check the solver for initialization
    ! This is utilized only by some solvers.
    !-----------------------------------------------------------------
    IF( IsProcedure ) THEN
      InitProc = GetProcAddr( TRIM(proc_name)//'_Init', abort=.FALSE. )
      IF ( InitProc /= 0 ) THEN
        CALL ExecSolver( InitProc, CurrentModel, Solver, &
            Solver % dt, Transient )
      END IF
    END IF

    Solver % SolverMode = SOLVER_MODE_DEFAULT
    IF( ListGetLogical( SolverParams, 'Auxiliary Solver', Found ) ) &
        Solver % SolverMode = SOLVER_MODE_AUXILIARY
    IF( ListGetLogical( SolverParams, 'Coupled Solver', Found ) ) &
        Solver % SolverMode = SOLVER_MODE_COUPLED
    IF( ListGetLogical( SolverParams, 'Block Solver', Found ) ) &
        Solver % SolverMode = SOLVER_MODE_BLOCK
    IF( ListGetLogical( SolverParams, 'Assembly Solver', Found ) ) &
        Solver % SolverMode = SOLVER_MODE_ASSEMBLY

    IF( Solver % SolverMode == SOLVER_MODE_DEFAULT ) THEN
      eq = ListGetString( Solver  % Values, 'Equation', Found )
      IF(.NOT. Found ) Solver % SolverMode = SOLVER_MODE_AUXILIARY 
    END IF

    IsCoupledSolver = ( Solver % SolverMode == SOLVER_MODE_COUPLED ) 
    IsBlockSolver = ( Solver % SolverMode == SOLVER_MODE_BLOCK ) 
    IsAssemblySolver = ( Solver % SolverMode == SOLVER_MODE_ASSEMBLY ) 
    IsAssemblySolver = IsAssemblySolver .OR. &
	( IsCoupledSolver .AND. .NOT. IsProcedure ) .OR. &
        ( IsBlockSolver .AND. .NOT. IsProcedure ) 

    ! Default order of equation
    !--------------------------
    Solver % Order = 1
    Solver % TimeOrder = 1
    
    ! Set up time-stepping strategies for transient problems
    !------------------------------------------------------------------------------
    IF ( Transient ) THEN
      str = ListGetString( SolverParams, 'Timestepping Method',Found )
      IF ( .NOT. Found ) THEN
        str = ListGetString( CurrentModel % Simulation, 'Timestepping Method',Found )
        IF ( Found ) THEN
          CALL ListAddString( SolverParams, 'Timestepping Method', str )
        END IF
      END IF
      
      IF ( Found ) THEN
        IF (str=='bdf') THEN
          Solver % Order = ListGetInteger( SolverParams, &
              'BDF Order', Found, minv=1, maxv=5 )
          IF ( .NOT. Found ) THEN
            Solver % Order = ListGetInteger( CurrentModel % &
                Simulation, 'BDF Order', Found, minv=1, maxv=5 )
          END IF
          IF ( .NOT.Found ) THEN
            Solver % Order = 2
            CALL Warn( 'AddEquation', 'BDF order defaulted to 2.' )
          END IF
        ELSE IF ( str=='runge-kutta') THEN
          Solver % Order = ListGetInteger( CurrentModel % &
              Simulation, 'Runge-Kutta Order', Found, minv=2, maxv=4 )
          IF ( .NOT.Found ) Solver % Order = 2
        END IF
      ELSE
        CALL Warn( 'AddEquation', '> Timestepping method < defaulted to > Implicit Euler <' )
        CALL ListAddString( SolverParams, 'Timestepping Method', 'Implicit Euler' )
      END IF
    END IF

    ! Get the procudure that really runs the solver
    !------------------------------------------------------------------------------
    IF( IsProcedure ) THEN
      Solver % PROCEDURE = GetProcAddr(proc_name)
    END IF
    
    ! Initialize and get the variable 
    !------------------------------------------------------------------------    
    Solver % TimeOrder = 0
    NULLIFY( Solver % Matrix )    
    var_name = ListGetString( SolverParams, 'Variable', Found )

    IF(.NOT. Found ) THEN
      ! Variable does not exist
      !------------------------------------------------------

      ALLOCATE( Solver % Variable )
      Solver % Variable % Name = ''
      Solver % Variable % NameLen = 0
      Solver % Variable % Norm = 0.0d0
      NULLIFY( Solver % Variable % Perm )
      NULLIFY( Solver % Variable % Values )
      
      
    ELSE IF( IsCoupledSolver .AND. .NOT. IsProcedure ) THEN
      ! Coupled solver may inherit the matrix only if procedure is given
      !-----------------------------------------------------------------

    ELSE IF( IsBlockSolver .AND. .NOT. IsProcedure ) THEN
      ! Block solver may inherit the matrix only if procedure is given
      !-----------------------------------------------------------------
      
    ELSE
      ! It may be a normal field variable or a global (0D) variable
      !------------------------------------------------------------------------
      VariableGlobal = ListGetLogical( SolverParams, 'Variable Global', Found )
      
      VariableOutput = ListGetLogical( SolverParams, 'Variable Output', Found )
      IF ( .NOT. Found ) VariableOutput = .TRUE.
      
      DOFs = ListGetInteger( SolverParams, 'Variable DOFs', Found, minv=1 )
      IF ( .NOT. Found ) THEN
        j = 0
        DOFs = 0
        DO WHILE( .TRUE. )
          i = INDEX( var_name(j+1:), ':' ) + j
          IF ( i<=j ) EXIT
          READ( var_name(i+1:),'(i1)' ) k
          DOFs = DOFs + k
          j = i + 1
        END DO
      END IF
      
      DO WHILE( var_name(1:1) == '-' )
        IF ( SEQL(var_name, '-nooutput ') ) THEN
          VariableOutput = .FALSE.
          var_name(1:LEN(var_name)-10) = var_name(11:)
        END IF
        
        IF ( SEQL(var_name, '-global ') ) THEN
          VariableGlobal = .TRUE.
          var_name(1:LEN(var_name)-8) = var_name(9:)
        END IF
        
        IF ( SEQL(var_name, '-dofs ') ) THEN
          READ( var_name(7:), * ) DOFs
          i = 7
          j = LEN_TRIM( var_name )
          DO WHILE( var_name(i:i) /= ' '  )
            i = i + 1
            IF ( i > j ) EXIT
          END DO
          var_name(1:LEN(var_name)-i) = var_name(i+1:)
        END IF
      END DO
      IF ( DOFs == 0 ) DOFs = 1
      
      n = LEN_TRIM(var_name)
      
      ! If the variable is 'global' it has nothing to do with the mesh and it may be simply 
      ! allocated. 
      !------------------------------------------------------------------------------------
      IF( VariableGlobal ) THEN
        CALL Info('AddEquationBasics','Creating global variable: '//var_name(1:n),Level=8)

        Solver % SolverMode = SOLVER_MODE_GLOBAL
        ALLOCATE( Solution( DOFs ) )
        Solution = 0.0_dp
        
        CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver, &
            var_name(1:n), DOFs, Solution )
        Solver % Variable => VariableGet( Solver % Mesh % Variables, var_name(1:n) )
        IF( DOFs > 1 ) THEN
          DO i=1,DOFs
            tmpname = ComponentName( var_name(1:n), i )
            Component => Solution( i:i )
            CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver,&
                tmpname, 1, Component )
          END DO
        END IF

        IF( ListGetLogical( SolverParams,'Ode Matrix',Found ) ) THEN
          CALL Info('AddEquationBasics','Creating dense matrix for ODE: '&
              //TRIM(I2S(Dofs)),Level=8)
          ALLOCATE( Solver % Variable % Perm(1) )
          Solver % Variable % Perm(1) = 1
          Solver % Matrix => CreateOdeMatrix( CurrentModel, Solver, Dofs )
        END IF

      ELSE        

        ! If the variable is a field variable create a permutation and matrix related to it
        !----------------------------------------------------------------------------------
        eq = ListGetString( SolverParams, 'Equation', Found )
        IF(.NOT. Found) THEN
          CALL Fatal('AddEquationBasics','Variable exists but > Equation < is not defined')
        END IF
        Found = .FALSE.
        DO i=1, CurrentModel % NumberOfEquations
          IF( ListGetLogical( CurrentModel % Equations(i) % Values, TRIM(eq), Stat)) THEN
            Found = .TRUE.
            EXIT
          END IF
        END DO
        IF(.NOT. Found ) THEN
          CALL Fatal('AddEquationBasics','Variable > '//var_name(1:n)//&
                     ' < exists but it is not associated to any equation')
        END IF
        
        ! Computate the size of the permutation vector
        !-----------------------------------------------------------------------------------------
        Ndeg = 0
!        IF( Solver % SolverMode == SOLVER_MODE_DEFAULT .OR. &
!            Solver % SolverMode == SOLVER_MODE_ASSEMBLY ) THEN
        IF(.TRUE.) THEN
          eq = ListGetString( Solver  % Values, 'Equation', Found )
          MaxNDOFs  = 0
          MaxDGDOFs = 0
          DO i=1,Solver % Mesh % NumberOFBulkElements
            CurrentElement => Solver % Mesh % Elements(i)
            MaxDGDOFs = MAX( MaxDGDOFs, CurrentElement % DGDOFs )
            MaxNDOFs  = MAX( MaxNDOFs,  CurrentElement % NDOFs )
          END DO
          
          MaxEDOFs = 0
          DO i=1,Solver % Mesh % NumberOFEdges 
            CurrentElement => Solver % Mesh % Edges(i)
            MaxEDOFs  = MAX( MaxEDOFs,  CurrentElement % BDOFs )
          END DO
          
          MaxFDOFs = 0
          DO i=1,Solver % Mesh % NumberOFFaces 
            CurrentElement => Solver % Mesh % Faces(i)
            MaxFDOFs  = MAX( MaxFDOFs,  CurrentElement % BDOFs )
          END DO
          
          MaxBDOFs = 0
          DO i=1,Solver % Mesh % NumberOFBulkElements
            CurrentElement => Solver % Mesh % Elements(i)
            MaxBDOFs  = MAX( MaxBDOFs,  CurrentElement % BDOFs )
          END DO
          
          GlobalBubbles = ListGetLogical( SolverParams, 'Bubbles in Global System', Found )
          IF (.NOT.Found) GlobalBubbles = .TRUE.

          Ndeg = Ndeg + Solver % Mesh % NumberOfNodes 
          IF ( MaxEDOFs > 0 ) Ndeg = Ndeg + MaxEDOFs * Solver % Mesh % NumberOFEdges
          IF ( MaxFDOFs > 0 ) Ndeg = Ndeg + MaxFDOFs * Solver % Mesh % NumberOFFaces
          IF ( GlobalBubbles ) &
              Ndeg = Ndeg + MaxBDOFs * Solver % Mesh % NumberOfBulkElements
          IF ( ListGetLogical( SolverParams, 'Discontinuous Galerkin', Found ) ) &
              Ndeg = MAX( NDeg, MaxDGDOFs * (Solver % Mesh % NumberOfBulkElements+ &
              Solver % Mesh % NumberOfBoundaryElements) )
        END IF
        
        IF( ListGetLogical( SolverParams,'Radiation Solver',Found ) ) THEN
          CALL RadiationFactors( Solver, .TRUE.)
        END IF
        
        BandwidthOptimize = ListGetLogical( SolverParams, &
            'Optimize Bandwidth', Found )
        IF ( .NOT. Found ) BandwidthOptimize = .TRUE.
        CALL CheckLinearSolverOptions( Solver )

        ALLOCATE( Perm(Ndeg) )

        Perm = 0
        MatrixFormat = MATRIX_CRS
        Solver % Matrix => CreateMatrix( CurrentModel, Solver, Solver % Mesh, &
            Perm, DOFs, MatrixFormat, BandwidthOptimize, eq(1:LEN_TRIM(eq)), &
            ListGetLogical( SolverParams,'Discontinuous Galerkin', Found ), &
            GlobalBubbles=GlobalBubbles )          
        Nrows = DOFs * Ndeg
        IF (ASSOCIATED(Solver % Matrix)) Nrows = Solver % Matrix % NumberOfRows
       
        ! Basically the solver could be matrix free but still the matrix
        ! is used here temperarily since it is needed when making the 
        ! permutation vector
        !-----------------------------------------------------------------
        IF( ListGetLogical( SolverParams, 'No Matrix', Found ) ) THEN
          Solver % SolverMode = SOLVER_MODE_MATRIXFREE
          CALL FreeMatrix( Solver % Matrix )
        END IF
       
        IF (Nrows>0) THEN
          ALLOCATE(Solution(Nrows))
          Solution = InitValue
          
          CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver, &
              var_name(1:n), DOFs, Solution, Perm, Output=VariableOutput )          
          Solver % Variable => VariableGet( Solver % Mesh % Variables, var_name(1:n) )
          
          IF ( DOFs > 1 ) THEN
            DO i=1,DOFs
              tmpname = ComponentName( var_name(1:n), i )
              Component => Solution( i:Nrows-DOFs+i:DOFs )
              CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver,&
                  tmpname, 1, Component, Perm, Output=VariableOutput )
            END DO
          END IF
        END IF

        IF (ASSOCIATED(Solver % Matrix)) Solver % Matrix % Comm = MPI_COMM_WORLD
        IF ( ListGetLogical( SolverParams, 'Discontinuous Galerkin', Found) ) &
          Solver % Variable % TYPE = Variable_on_nodes_on_elements
      END IF
      !------------------------------------------------------------------------------
    END IF

    !------------------------------------------------------------------------------
    ! Add the exported variables which are typically auxiliary variables derived
    ! from the solution without their own matrix equation.  
    !------------------------------------------------------------------------------
    l = 0
    DO WHILE( .TRUE. )
      l = l + 1
      str = ComponentName( 'exported variable', l )
      var_name = ListGetString( SolverParams, str, Found )
      
      IF(.NOT. Found) EXIT
      
      str = TRIM( ComponentName( 'exported variable', l ) ) // ' Output'
      VariableOutput = ListGetLogical( SolverParams, str, Found )
      IF ( .NOT. Found ) VariableOutput = .TRUE.
      
      str = TRIM( ComponentName( 'exported variable', l ) ) // ' DOFs'
      DOFs = ListGetInteger( SolverParams, str, Found )
      IF ( .NOT. Found ) THEN
        j = 0
        DOFs = 0
        DO WHILE( .TRUE. )
          i = INDEX( var_name(j+1:), ':' ) + j
          IF ( i<=j ) EXIT
          READ( var_name(i+1:),'(i1)' ) k
          DOFs = DOFs + k
          j = i + 1
        END DO
      END IF
      
      VariableOutput = .TRUE.
      VariableGlobal = .FALSE.
      
      DO WHILE( var_name(1:1) == '-' )
        IF ( SEQL(var_name, '-nooutput ') ) THEN
          VariableOutput = .FALSE.
          var_name(1:LEN(var_name)-10) = var_name(11:)
        END IF
        
        IF ( SEQL(var_name, '-global ') ) THEN
          VariableGlobal = .TRUE.
          var_name(1:LEN(var_name)-8) = var_name(9:)
        END IF
        
        IF ( SEQL(var_name, '-dofs ') ) THEN
          READ( var_name(7:), * ) DOFs 
          j = LEN_TRIM( var_name )
          k = 7
          DO WHILE( var_name(k:k) /= ' '  )
            k = k + 1
            IF ( k > j ) EXIT
          END DO
          var_name(1:LEN(var_name)-(k+2)) = var_name(k+1:)
        END IF
      END DO
      IF ( DOFs == 0 ) DOFs = 1
      
      NewVariable => VariableGet( Solver % Mesh % Variables, Var_name )
      
      IF ( .NOT. ASSOCIATED(NewVariable) ) THEN
        IF( VariableGlobal ) THEN
          nSize = DOFs
          NULLIFY( Perm )
        ELSE
          nSize = DOFs * SIZE(Solver % Variable % Values) / Solver % Variable % DOFs
          Perm => Solver % Variable % Perm
        END IF
        
        ALLOCATE( Solution(nSize) )
        Solution = 0.0d0
        IF( ASSOCIATED(Perm) ) THEN
          CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver,&
              var_name, DOFs, Solution, Solver % Variable % Perm, &
              Output=VariableOutput, TYPE=Solver % Variable % TYPE )
        ELSE          
          CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver,&
              var_name, DOFs, Solution, TYPE=Solver % Variable % TYPE )
        END IF
        
        IF ( DOFs > 1 ) THEN
!.AND. .NOT. VariableGlobal ) THEN
          n = LEN_TRIM( var_name )
          DO j=1,DOFs
            tmpname = ComponentName( var_name(1:n), j )
            Component => Solution( j:nSize-DOFs+j:DOFs )
            IF( ASSOCIATED(Perm) ) THEN
              CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver,&
                  tmpname, 1, Component, Perm,  &
                  Output=VariableOutput, TYPE=Solver % Variable % TYPE )
            ELSE
              CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver,&
                  tmpname, 1, Component, TYPE=Solver % Variable % TYPE )
            END IF
          END DO
        END IF
      END IF
    END DO


    !------------------------------------------------------------------
    ! Check for special solvers, to be executed only 
    ! at a certain instances during the simulation:
    !------------------------------------------------------------------

    Solver % SolverExecWhen = SOLVER_EXEC_ALWAYS

    SELECT CASE( ListGetString( SolverParams, 'Exec Solver', Found )  )
      CASE( 'never' )
      Solver % SolverExecWhen = SOLVER_EXEC_NEVER
      CASE( 'always' )
      Solver % SolverExecWhen = SOLVER_EXEC_ALWAYS
      CASE( 'after simulation', 'after all' )
      Solver % SolverExecWhen = SOLVER_EXEC_AFTER_ALL
      CASE( 'before simulation', 'before all' )
         Solver % SolverExecWhen = SOLVER_EXEC_AHEAD_ALL
      CASE( 'before timestep' )
         Solver % SolverExecWhen = SOLVER_EXEC_AHEAD_TIME
      CASE( 'after timestep' )
         Solver % SolverExecWhen = SOLVER_EXEC_AFTER_TIME
      CASE( 'before saving' )
         Solver % SolverExecWhen = SOLVER_EXEC_AHEAD_SAVE
      CASE( 'after saving' )
         Solver % SolverExecWhen = SOLVER_EXEC_AFTER_SAVE
      CASE DEFAULT
         Solver % SolverExecWhen = SOLVER_EXEC_ALWAYS
    END SELECT

    IF ( ListGetLogical( SolverParams, 'Before All', Found ) ) THEN
       Solver % SolverExecWhen = SOLVER_EXEC_AHEAD_ALL
    ELSE IF ( ListGetLogical( SolverParams, 'Before Simulation', Found ) ) THEN
       Solver % SolverExecWhen = SOLVER_EXEC_AHEAD_ALL
    ELSE IF ( ListGetLogical( SolverParams, 'After All', Found ) ) THEN
       Solver % SolverExecWhen = SOLVER_EXEC_AFTER_ALL
    ELSE IF ( ListGetLogical( SolverParams, 'After Simulation', Found ) ) THEN
       Solver % SolverExecWhen = SOLVER_EXEC_AFTER_ALL
    ELSE IF ( ListGetLogical( SolverParams, 'Before Timestep', Found ) ) THEN
       Solver % SolverExecWhen = SOLVER_EXEC_AHEAD_TIME
    ELSE IF ( ListGetLogical( SolverParams, 'After Timestep', Found ) ) THEN
       Solver % SolverExecWhen = SOLVER_EXEC_AFTER_TIME
    ELSE IF ( ListGetLogical( SolverParams, 'Before Saving', Found ) ) THEN
       Solver % SolverExecWhen = SOLVER_EXEC_AHEAD_SAVE
    ELSE IF ( ListGetLogical( SolverParams, 'After Saving', Found ) ) THEN
       Solver % SolverExecWhen = SOLVER_EXEC_AFTER_SAVE
    END IF

    Solver % LinBeforeProc = 0
    str = ListGetString( Solver % Values, 'Before Linsolve', Found )
    IF ( Found ) Solver % LinBeforeProc = GetProcAddr( str )

    Solver % LinAfterProc = 0
    str = ListGetString( Solver % Values, 'After Linsolve', Found )
    IF ( Found ) Solver % LinAfterProc = GetProcAddr( str )

    IF( ASSOCIATED( Solver % Matrix ) ) THEN
      Solver % Matrix % MatVecSubr = 0
      str = ListGetString( Solver % Values, 'Matrix Vector Proc', Found )
      IF ( Found ) Solver % Matrix % MatVecSubr = GetProcAddr( str )
    END IF

    Solver % MortarProc = 0
    str = ListGetString( Solver % Values, 'External Projector Procedure', Found )
    IF ( Found ) Solver % MortarProc = GetProcAddr( str )

  END SUBROUTINE AddEquationBasics


!------------------------------------------------------------------------------  
!> Add information that is typically only needed if there's a matrix equation
!> to work with. This should be called only after both the solution vector and
!> matrix have been created.
!------------------------------------------------------------------------------
  SUBROUTINE AddEquationSolution(Solver, Transient )
!------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: Solver
    LOGICAL :: Transient
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var
    CHARACTER(LEN=MAX_NAME_LEN) :: str, var_name, tmpname
    INTEGER :: i,j,k,n,nrows,DOFs
    REAL(KIND=dp), POINTER :: Solution(:)
    INTEGER, POINTER :: Perm(:)
    LOGICAL :: Found, Stat, ComplexFlag, HarmonicAnal, EigAnal, &
         VariableOutput, MGAlgebraic,MultigridActive
    INTEGER :: MgLevels
    REAL(KIND=dp), POINTER :: Component(:)
    REAL(KIND=dp), POINTER :: freqv(:,:)
    REAL(KIND=dp) :: Freq
    TYPE(Mesh_t),   POINTER :: NewMesh,OldMesh
    TYPE(Element_t), POINTER :: CurrentElement
    TYPE(Matrix_t), POINTER :: NewMatrix, Oldmatrix, SaveMatrix

    !------------------------------------------------------------------------------
    Solver % DoneTime = 0
    IF ( .NOT. ASSOCIATED( Solver % Variable ) ) RETURN
    IF ( .NOT. ASSOCIATED( Solver % Variable % Values ) ) RETURN
    !------------------------------------------------------------------------------
	
    !------------------------------------------------------------------------------
    ! If soft limiters are applied then also loads must be computed
    !------------------------------------------------------------------------------
     IF( ListGetLogical( Solver % Values,'Calculate Boundary Fluxes',Found) ) THEN
       CALL ListAddLogical( Solver % Values,'Calculate Loads',.TRUE.)
     END IF
	    
    !------------------------------------------------------------------------------
    ! Create the variable needed for the computation of nodal loads: r=b-Ax
    !------------------------------------------------------------------------------
    IF ( ListGetLogical( Solver % Values,'Calculate Loads', Found ) ) THEN
      Var_name = GetVarName(Solver % Variable) // ' Loads'
      Var => VariableGet( Solver % Mesh % Variables, var_name )
      IF ( .NOT. ASSOCIATED(Var) ) THEN
        ALLOCATE( Solution(SIZE(Solver % Variable % Values)) )
        DOFs = Solver % Variable % DOFs
        Solution = 0.0d0
        nrows = SIZE( Solution ) 
        Perm => Solver % Variable % Perm
        VariableOutput = Solver % Variable % Output

        CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver,&
            var_name, Solver % Variable % DOFs, Solution, &
            Solver % Variable % Perm, Output=VariableOutput )

        IF ( DOFs > 1 ) THEN
          n = LEN_TRIM( Var_name )
          DO j=1,DOFs
            tmpname = ComponentName( var_name(1:n), j )
            Component => Solution( j:nRows-DOFs+j:DOFs )
            CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver,&
                tmpname, 1, Component, Perm, Output=VariableOutput )
          END DO
        END IF
        NULLIFY( Solution )
      END IF
    END IF

    ! Create a additional variable for the limiters. For elasticity, for example
    ! the variable will be the contact load. 
    IF ( ListGetLogical( Solver % Values,'Apply Limiter', Found ) .OR. &
        ListGetLogical( Solver % Values,'Apply Contact BCs',Found ) ) THEN
      Var_name = GetVarName(Solver % Variable) // ' Contact Load'
      Var => VariableGet( Solver % Mesh % Variables, var_name )
      IF ( .NOT. ASSOCIATED(Var) ) THEN
        ALLOCATE( Solution(SIZE(Solver % Variable % Values)) )
        DOFs = Solver % Variable % DOFs
        Solution = 0.0d0
        nrows = SIZE( Solution ) 
        Perm => Solver % Variable % Perm
        VariableOutput = Solver % Variable % Output

        CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver,&
            var_name, Solver % Variable % DOFs, Solution, &
            Solver % Variable % Perm, Output=VariableOutput )

        IF ( DOFs > 1 ) THEN
          n = LEN_TRIM( Var_name )
          DO j=1,DOFs
            tmpname = ComponentName( var_name(1:n), j )
            Component => Solution( j:nRows-DOFs+j:DOFs )
            CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver,&
                tmpname, 1, Component, Perm, Output=VariableOutput )
          END DO
        END IF
        NULLIFY( Solution )
      END IF
    END IF

    Solver % NOFEigenValues = 0
    Solver % MultiGridLevel = 1
    Solver % MultiGridTotal = 0
    Solver % MultiGridSolver = .FALSE.
    Solver % MultiGridEqualSplit = .FALSE.


    HarmonicAnal = ListGetLogical( Solver % Values, 'Harmonic Analysis', Found )
    
    IF ( ASSOCIATED( Solver % Matrix ) ) THEN
      IF(.NOT. ASSOCIATED(Solver % Matrix % RHS)) THEN
        ALLOCATE( Solver % Matrix % RHS(Solver % Matrix % NumberOFRows) )
        Solver % Matrix % RHS = 0.0d0
        
        Solver % Matrix % RHS_im => NULL()
        IF ( HarmonicAnal ) THEN
          ALLOCATE( Solver % Matrix % RHS_im(Solver % Matrix % NumberOFRows) )
          Solver % Matrix % RHS_im = 0.0d0
        END IF
      END IF
    END IF
!------------------------------------------------------------------------------

    EigAnal = ListGetLogical( Solver % Values, 'Eigen Analysis', Found )

    IF ( Transient .AND. .NOT. EigAnal .AND. .NOT. HarmonicAnal ) THEN
      k = ListGetInteger( Solver % Values, 'Time Derivative Order', Found, &
          minv=0, maxv=2 )
      Solver % TimeOrder = 1
      IF ( Found ) Solver % TimeOrder = MIN(MAX(1,k),2)
      
      IF ( ASSOCIATED( Solver % Matrix ) ) THEN
        ALLOCATE( Solver % Matrix % Force(Solver % Matrix % NumberOFRows, &
            Solver % TimeOrder+1) )
        Solver % Matrix % Force = 0.0d0
      END IF
      
      IF ( .NOT. ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
        IF ( Solver % TimeOrder == 2 ) THEN
          ALLOCATE( Solver % Variable % PrevValues( &
              SIZE(Solver % Variable % Values),5) )
        ELSE IF ( Solver % Order > Solver % TimeOrder ) THEN
          ALLOCATE(Solver % Variable % PrevValues( &
              SIZE(Solver % Variable % Values), Solver % Order))
        ELSE
          ALLOCATE(Solver % Variable % PrevValues( &
              SIZE(Solver % Variable % Values), Solver % TimeOrder))
        END IF
        Solver % Variable % PrevValues = 0.0d0
        
        IF ( Solver % Variable % DOFs > 1 ) THEN
          IF ( GetVarName( Solver % Variable ) == 'flow solution' ) THEN
            DO k=1,Solver % Variable % DOFs-1
              str = 'Velocity ' // CHAR(k+ICHAR('0'))
              Var => VariableGet( Solver % Mesh % Variables, str, .TRUE. )
              IF(.NOT. ASSOCIATED(Var)) THEN
                 WRITE (Message, '(a,a,a)')  "Failed to get variable:&
                      & ",TRIM(str), ". Try specifying Variable =&
                      & Flow Solution[velocity:DIM pressure:1] in .SIF"
                 CALL Fatal("AddEquationSolution",Message)
              END IF
              Var % PrevValues =>  &
                  Solver % Variable % PrevValues(k::Solver % Variable % DOFs,:)
            END DO
            Var => VariableGet( Solver % Mesh % Variables, 'Pressure', .TRUE. )
            IF(.NOT. ASSOCIATED(Var)) THEN
               CALL Fatal("AddEquationSolution","Failed to get variable: &
                    &pressure. Try specifying Variable = Flow Solution &
                    &[velocity:DIM pressure:1] in .SIF")
            END IF
            Var % PrevValues =>  &
                Solver % Variable % PrevValues(k::Solver % Variable % DOFs,:)
          ELSE
            DO k=1,Solver % Variable % DOFs
              str = ComponentName( Solver % Variable % Name, k ) 
              Var => VariableGet( Solver % Mesh % Variables, str, .TRUE. )
              Var % PrevValues =>  &
                  Solver % Variable % PrevValues(k::Solver % Variable % DOFs,:)
            END DO
          END IF
        END IF
      END IF


      IF( ListGetLogical( Solver % Values,'Calculate Velocity',Found) .OR. &
        ListGetLogical( Solver % Values,'Nonlinear Calculate Velocity',Found) ) THEN
	IF( Solver % TimeOrder < 1 ) THEN
          CALL Warn('AddEquationSolution',&
              'Velocity computation implemented only for time-dependent equations')
        ELSE IF ( Solver % TimeOrder == 1 ) THEN
          str = TRIM( Solver % Variable % Name ) // ' Velocity'
          CALL VariableAddVector( Solver % Mesh % Variables, Solver % Mesh, Solver, &
              str, Solver % Variable % Dofs, Perm = Solver % Variable % Perm )
        ELSE IF ( Solver % TimeOrder >= 2 ) THEN
          Component => Solver % Variable % PrevValues(:,1)
          str = TRIM( Solver % Variable % Name ) // ' Velocity'
          CALL VariableAddVector( Solver % Mesh % Variables, Solver % Mesh, Solver, &
              str, Solver % Variable % Dofs, Component, Solver % Variable % Perm, &
              Secondary = .TRUE. )
        END IF
      END IF

      IF( ListGetLogical( Solver % Values,'Calculate Acceleration',Found) ) THEN
        IF ( Solver % TimeOrder == 1 ) THEN
          CALL Warn('AddEquationSolution',&
              'Acceleration computation implemented only for 2nd order time equations')
        ELSE IF ( Solver % TimeOrder >= 2 ) THEN
          Component => Solver % Variable % PrevValues(:,2)
          str = TRIM( Solver % Variable % Name ) // ' Acceleration'
          CALL VariableAddVector( Solver % Mesh % Variables, Solver % Mesh, Solver, &
              str, Solver % Variable % Dofs, Component, Solver % Variable % Perm, &
              Secondary = .TRUE. )
        END IF
      END IF

    ELSE
      Solver % TimeOrder = 0

      IF( ListGetLogical( Solver % Values,'Calculate Derivative',Found) ) THEN
        str = TRIM( Solver % Variable % Name ) // ' Derivative'
        CALL VariableAddVector( Solver % Mesh % Variables, Solver % Mesh, Solver, &
            str, Solver % Variable % Dofs, Perm = Solver % Variable % Perm )
      END IF

      IF ( EigAnal ) THEN
        !ComplexFlag = ListGetLogical( Solver % Values,  'Eigen System Complex', Found )
        !IF ( .NOT. Found ) ComplexFlag = .FALSE.
        
        n = ListGetInteger( Solver % Values,  'Eigen System Values', Found )
        IF ( Found .AND. n > 0 ) THEN
          Solver % NOFEigenValues = n
          IF ( .NOT. ASSOCIATED( Solver % Variable % EigenValues ) ) THEN
            ALLOCATE( Solver % Variable % EigenValues(n) )
            ALLOCATE( Solver % Variable % EigenVectors(n, &
                SIZE( Solver % Variable % Values ) ) )
            
            Solver % Variable % EigenValues  = 0.0d0
            Solver % Variable % EigenVectors = 0.0d0
            
            DO k=1,Solver % Variable % DOFs
              str = ComponentName( Solver % Variable % Name, k )
              Var => VariableGet( Solver % Mesh % Variables, str, .TRUE. )
              IF ( ASSOCIATED( Var ) ) THEN
                Var % EigenValues => Solver % Variable % EigenValues
                !IF ( ComplexFlag ) THEN
                !  Var % EigenVectors => Solver % Variable % EigenVectors
                !ELSE
                Var % EigenVectors =>  & 
                     Solver % Variable % EigenVectors(:,k::Solver % Variable % DOFs )
                !END IF
              END IF
            END DO
          END IF
          ALLOCATE( Solver % Matrix % MassValues(SIZE(Solver % Matrix % Values)) )
          Solver % Matrix % MassValues = 0.0d0
        END IF
      ELSE IF ( HarmonicAnal ) THEN

        n = ListGetInteger( Solver % Values,'Harmonic System Values',Found )
        IF( n > 1 ) THEN
          freqv => ListGetConstRealArray( Solver % Values, 'Frequency', Found )
          IF(Found ) THEN
            IF( SIZE( Freqv,1) < n ) THEN
              CALL Fatal( 'AddEquation', 'Frequency must be at least same size as > Harmonic System Values <')
            END IF
          ELSE
            CALL Fatal( 'AddEquation', '> Frequency < must be given for harmonic analysis.' )
          END IF
        ELSE
          n = 1
        END IF

        Solver % NOFEigenValues = n
        IF ( .NOT. ASSOCIATED( Solver % Variable % EigenValues ) ) THEN
          ALLOCATE( Solver % Variable % EigenValues(n) )
          ALLOCATE( Solver % Variable % EigenVectors(n, &
              SIZE( Solver % Variable % Values ) ) )
          
          Solver % Variable % EigenValues  = 0.0d0
          Solver % Variable % EigenVectors = 0.0d0
          
          DO k=1,Solver % Variable % DOFs
            str = ComponentName( Solver % Variable % Name, k )
            Var => VariableGet( Solver % Mesh % Variables, str, .TRUE. )
            IF ( ASSOCIATED( Var ) ) THEN
              Var % EigenValues => Solver % Variable % EigenValues
              Var % EigenVectors => &
                  Solver % Variable % EigenVectors(:,k::Solver % Variable % DOFs)
            END IF
          END DO
        END IF
        ALLOCATE( Solver % Matrix % MassValues(SIZE(Solver % Matrix % Values)) )
        Solver % Matrix % MassValues = 0.0d0
      END IF
    END IF


!------------------------------------------------------------------------------

    IF ( ASSOCIATED( Solver % Matrix ) ) THEN
      Solver % Matrix % Symmetric = ListGetLogical( Solver % Values, &
          'Linear System Symmetric', Found )
      
      Solver % Matrix % Lumped = ListGetLogical( Solver % Values, &
          'Lumped Mass Matrix', Found )
      
      MultigridActive = &
          ListGetString( Solver % Values, 'Linear System Solver', Found ) == 'multigrid' .OR. &
          ListGetString( Solver % Values, 'Linear System Preconditioning', Found ) == 'multigrid'


!      Check for multigrid solver:
!      ---------------------------
       IF ( MultigridActive ) THEN

         ! Multigrid may be either solver or preconditioner, it it solver?
         Solver % MultiGridSolver = ListGetString(Solver % Values, &
             'Linear System Solver', Found ) == 'multigrid'
 
         ! There are four different methods: algrabraic, cluster, p and geometric
         str = ListGetString( Solver % Values,'MG Method',Found) 
         IF( Found ) THEN
           MGAlgebraic = ( str == 'algebraic' ) .OR. ( str == 'cluster') .OR. ( str == 'p')
         ELSE    
           MGAlgebraic = ListGetLogical( Solver % Values, 'MG Algebraic', Found ) &
               .OR. ListGetLogical( Solver % Values, 'MG Cluster', Found ) &
               .OR. ListGetLogical( Solver % Values, 'MG Pelement', Found ) 
         END IF
         
         
         MgLevels = ListGetInteger( Solver % Values, &
             'MG Levels', Found, minv=1 )         
         IF ( .NOT. Found ) THEN
           MgLevels = ListGetInteger( Solver % Values, &
               'Multigrid Levels', Found, minv=1 )
         END IF
         IF ( .NOT. Found ) THEN
           IF( MGAlgebraic ) THEN 
             MgLevels = 10
           ELSE
             CALL Fatal('AddEquationSolution','> MG Levels < must be defined for geometric multigrid!')
           END IF
         END IF
         Solver % MultiGridTotal = MgLevels
 

!         In case of geometric multigrid make the hierarchical meshes
!         -----------------------------------------------------------
         IF(.NOT. MGAlgebraic ) THEN 

!         Check if h/2 splitting of mesh requested:
!         ------------------------------------------
           Solver % MultiGridEqualSplit = ListGetLogical( &
               Solver % Values, 'MG Equal Split', Found )
         
           IF ( Solver % MultiGridEqualSplit ) THEN
             CALL ParallelInitMatrix( Solver, Solver % Matrix )
             Solver % MultiGridLevel = 1
             
             DO WHILE( Solver % MultiGridLevel < Solver % MultiGridTotal )
               IF ( ASSOCIATED( Solver % Mesh % Child ) ) THEN
                 NewMesh => Solver % Mesh % Child
                 
                 OldMesh   => Solver % Mesh
                 OldMatrix => Solver % Matrix
                 
                 CALL UpdateSolverMesh( Solver, NewMesh )
                 Solver % Mesh % Changed = .FALSE.
               ELSE
                 NewMesh => SplitMeshEqual( Solver % Mesh )
                 NewMesh % Next => CurrentModel % Meshes
                 CurrentModel % Meshes => NewMesh
                 
                 OldMesh   => Solver % Mesh
                 OldMatrix => Solver % Matrix
                 
                 CALL UpdateSolverMesh( Solver, NewMesh )
                 Solver % Mesh % Changed = .FALSE.
                 
                 NewMesh % Parent => OldMesh
                 OldMesh % Child  => NewMesh
                 NewMesh % Name = OldMesh % Name
               END IF
               Newmesh % OutputActive = .TRUE.
               OldMesh % OutputActive = .FALSE.
               
               NewMatrix => Solver % Matrix
               NewMatrix % Parent => OldMatrix
               OldMatrix % Child  => NewMatrix
               CALL ParallelInitMatrix( Solver, Solver % Matrix )
               Solver % MultiGridLevel = Solver % MultiGridLevel + 1
             END DO
           ELSE
             CALL ParallelInitMatrix( Solver, Solver % Matrix )
             OldMesh   => Solver % Mesh
             Var => Solver % Variable
             SaveMatrix => Solver % Matrix
             
             Solver % MultiGridLevel = 1
             DO WHILE( Solver % MultiGridLevel < Solver % MultigridTotal )
               IF ( ASSOCIATED(Solver % Mesh % Parent) ) THEN
                 NewMesh => Solver % Mesh % Parent
                 OldMatrix => Solver % Matrix
                 CALL UpdateSolverMesh( Solver, NewMesh )
                 Solver % Mesh % Changed = .FALSE.
                 NewMatrix => Solver % Matrix
                 NewMatrix % Child => OldMatrix
                 OldMatrix % Parent  => NewMatrix
                 CALL ParallelInitMatrix(Solver, Solver % Matrix )
               END IF
               Solver % MultiGridLevel = Solver % MultiGridLevel+1
             END DO
             
             Solver % Mesh => OldMesh
             Solver % Variable => Var
             Solver % Matrix => SaveMatrix
             CALL SetCurrentMesh(CurrentModel,Solver % Mesh)
           END IF

           CALL MeshStabParams( Solver % Mesh )
         END IF

         Solver % MultiGridLevel = Solver % MultigridTotal

        END IF

       ! Set the default verbosity of the iterative solvers accordingly with the 
       ! global verbosity.
       !-----------------------------------------------------------------------------
       IF( .NOT. ListCheckPresent( Solver % Values,'Linear System Residual Output') ) THEN
         k = 1
         IF( .NOT. OutputLevelMask(4) ) THEN
           k = 0
         ELSE IF( .NOT. OutputLevelMask(5) ) THEN
           k = 10
         END IF
         IF( k /= 1 ) THEN
           CALL ListAddInteger( Solver % Values,'Linear System Residual Output',k)
         END IF
       END IF
     END IF
     
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   END SUBROUTINE AddEquationSolution
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solve the equations one-by-one. 
!------------------------------------------------------------------------------
  SUBROUTINE SolveEquations( Model, dt, TransientSimulation, &
      CoupledMinIter, CoupledMaxIter, SteadyStateReached, RealTimestep )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: dt
    INTEGER, OPTIONAL :: RealTimestep
    INTEGER :: CoupledMinIter, CoupledMaxIter
    LOGICAL :: TransientSimulation, SteadyStateReached, TransientSolver
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: RelativeChange, Tolerance, PrevDT = 0.0d0, Relaxation, SSCond
    INTEGER :: i,j,k,l,n,ierr,istat,Visited=0, RKorder=0, nSolvers
    LOGICAL :: Found, Stat, AbsNorm, Scanning, Convergence, RungeKutta, MeActive, &
       NeedSol, CalculateDerivative, TestConvergence
    LOGICAL, ALLOCATABLE :: DoneThis(:), AfterConverged(:)
    TYPE(Solver_t), POINTER :: Solver
    TYPE(Mesh_t),   POINTER :: Mesh
    CHARACTER(LEN=max_name_len) :: When
    TYPE(Variable_t), POINTER :: IterV, TimestepV
    REAL(KIND=dp), POINTER :: steadyIt,nonlnIt
    REAL(KIND=dp), POINTER :: k1(:), k2(:), k3(:), k4(:), sTime

    TYPE(Variable_t), POINTER :: TimeVar
    REAL(KIND=dp) :: RK2_err
    REAL(KIND=dp), ALLOCATABLE :: RK2_ErrorEstimate(:)
    TYPE RungeKutta_t
      REAL(KIND=dp), ALLOCATABLE :: k1(:),k2(:),k3(:),k4(:)
    END TYPE RungeKutta_t
    TYPE(RungeKutta_t), ALLOCATABLE, TARGET :: RKCoeff(:)
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!   Intialize equation solvers for new timestep
!------------------------------------------------------------------------------
    nSolvers = Model % NumberOfSolvers

    Scanning = &
      ListGetString( CurrentModel % Simulation, 'Simulation Type', Found ) == 'scanning'

    IF ( TransientSimulation ) THEN
       DO k=1,nSolvers
          Solver => Model % Solvers(k)
          IF ( Solver % PROCEDURE /= 0 ) CALL InitializeTimestep(Solver)
       END DO
    END IF

    IF( TransientSimulation .OR. Scanning ) THEN
       IterV => VariableGet(Model % Solvers(1) % Mesh % Variables, 'coupled iter')
       steadyIt => IterV % Values(1)
    END IF
!------------------------------------------------------------------------------

    DO k=1,nSolvers
      Solver => Model % Solvers(k)
      IF ( Solver % PROCEDURE==0 ) CYCLE

      when  = ListGetString( Solver % Values, 'Exec Solver', Found )
      IF ( Found ) THEN
        IF ( When == 'before timestep' ) THEN
          CALL SolverActivate( Model,Solver,dt,TransientSimulation )
          CALL ParallelBarrier
        END IF
      ELSE
        IF ( Solver % SolverExecWhen == SOLVER_EXEC_AHEAD_TIME ) THEN
          CALL SolverActivate( Model,Solver,dt,TransientSimulation )
          CALL ParallelBarrier
        END IF
      END IF
    END DO

!------------------------------------------------------------------------------

    ALLOCATE( DoneThis(nSolvers), AfterConverged(nSolvers) )

    DO i=1,nSolvers
      Solver => Model % Solvers(i)
      AfterConverged(i) = ListGetLogical( Solver % Values, &
             'Coupled System After Others Converged', Found )
    END DO

!------------------------------------------------------------------------------
    RungeKutta = ListGetString( Model % Simulation, &
           'Timestepping Method', Found ) == 'runge-kutta'

    TimeVar => VariableGet( Model % Variables, 'Time')
    sTime => TimeVar % Values(1)
    IF(RungeKutta) sTime = sTime - dt
    CALL SolveCoupled()
    IF(RungeKutta) sTime = sTime + dt

    DO i=1,nSolvers
      Solver => Model % Solvers(i)
      IF ( Solver % PROCEDURE==0 ) CYCLE

      RungeKutta = .FALSE.
      IF ( TransientSimulation .AND. Solver % TimeOrder == 1 ) THEN
          RungeKutta = ListGetString( Solver % Values, &
              'Timestepping Method', Found ) == 'runge-kutta'
      END IF
      IF ( .NOT. RungeKutta ) CYCLE

      IF ( .NOT. ALLOCATED(RKCoeff) ) THEN
        ALLOCATE(RKCoeff(nSolvers), RK2_ErrorEstimate(nSolvers))
      END IF

      n = SIZE(Solver % Variable % Values)
      ALLOCATE(RKCoeff(i) % k1(n), RKCoeff(i) % k2(n), RKCoeff(i) % k3(n), &
                RKCoeff(i) % k4(n) )

      RKorder = Solver % Order
      k1 => RKCoeff(i) % k1
      k1 = Solver % Variable % Values-Solver % Variable % PrevValues(:,1)
      Solver % Variable % Values = Solver % Variable % PrevValues(:,1) + 2*k1/RKorder
    END DO

    IF ( ALLOCATED(RKCoeff) ) THEN
      IF(RKorder==4) THEN
        dt = dt / 2
        sTime = sTime - dt
      END IF
      CALL SolveCoupled()

      RK2_ErrorEstimate = 0._dp
      DO i=1,nSolvers
        Solver => Model % Solvers(i)
        IF ( .NOT. ALLOCATED(RKCoeff(i) % k1)) CYCLE

        k1 => RKCoeff(i) % k1
        k2 => RKCoeff(i) % k2
        k2 = Solver % Variable % Values - Solver % Variable % PrevValues(:,1)
        SELECT CASE(RKorder)
        CASE(2)
          Solver % Variable % Values = Solver % Variable % PrevValues(:,1) + (k1+k2)/2
          RK2_Errorestimate(i) = SUM(((k2-k1)/2)**2)
        CASE(4)
          Solver % Variable % Values = Solver % Variable % PrevValues(:,1) + k2
          k2 = 2*k2
        END SELECT
      END DO

      IF ( RKorder==2 ) THEN
        ! Provide error measure for adaptive timestepping = ||RK2 (Heun) - Explicit Euler|| !
        RK2_err = 0._dp
        DO i=1,nSolvers
          RK2_ErrorEstimate(i) = SQRT( ParallelReduction(RK2_ErrorEstimate(i)) )
          RK2_err = MAX(RK2_err, RK2_ErrorEstimate(i) )
        END DO 
        CALL ListAddConstReal( Model % Simulation, 'Adaptive Error Measure', &
                   RK2_err / Solver % Variable % Norm )
      ELSE ! RK4
        CALL SolveCoupled()

        DO i=1,nSolvers
          Solver => Model % Solvers(i)
          IF ( .NOT. ALLOCATED(RKCoeff(i) % k1)) CYCLE

          k3 => RKCoeff(i) % k3
          k3 = 2*(Solver % Variable % Values - Solver % Variable % PrevValues(:,1))
          Solver % Variable % Values = Solver % Variable % PrevValues(:,1) + k3
        END DO

        sTime = sTime + dt
        dt = 2 * dt
        CALL SolveCoupled()

        DO i=1,nSolvers
          Solver => Model % Solvers(i)
          IF ( .NOT. ALLOCATED(RKCoeff(i) % k1)) CYCLE

          k1 => RKCoeff(i) % k1
          k2 => RKCoeff(i) % k2
          k3 => RKCoeff(i) % k3
          k4 => RKCoeff(i) % k4
          k4 = Solver % Variable % Values - Solver % Variable % PrevValues(:,1)

          Solver % Variable % Values = Solver % Variable % PrevValues(:,1) + &
                             ( k1 + 2*k2 + 2*k3 + k4 ) / 6
        END DO
      END IF
        
      DO i=1,nSolvers
        Solver => Model % Solvers(i)
        IF ( ALLOCATED(RKCoeff(i) % k1) ) THEN
          DEALLOCATE(RKCoeff(i) % k1, RKCoeff(i) % k2, RKCoeff(i) % k3, RKCoeff(i) % k4)
        END IF
        Solver % Variable % Norm = ComputeNorm(Solver, n, Solver % Variable % Values)
      END DO
      DEALLOCATE(RKCoeff)
    END IF
!------------------------------------------------------------------------------

    DO k=1,nSolvers
      Solver => Model % Solvers(k)
      IF ( Solver % PROCEDURE==0 ) CYCLE

      When = ListGetString( Solver % Values, 'Exec Solver', Found )
      IF (  Found ) THEN
         IF ( When == 'after timestep' ) THEN
           CALL SolverActivate( Model,Solver,dt,TransientSimulation )
           CALL ParallelBarrier
         END IF
      ELSE
         IF ( Solver % SolverExecWhen == SOLVER_EXEC_AFTER_TIME ) THEN
           CALL SolverActivate( Model,Solver,dt,TransientSimulation )
           CALL ParallelBarrier
         END IF
      END IF
    END DO

!------------------------------------------------------------------------------
    IF ( .NOT.TransientSimulation ) SteadyStateReached = ALL(DoneThis)
!------------------------------------------------------------------------------
    DEALLOCATE( DoneThis, AfterConverged )
!------------------------------------------------------------------------------

CONTAINS

    SUBROUTINE SolveCoupled()
!------------------------------------------------------------------------------
    DO i=1,CoupledMaxIter
       IF ( TransientSimulation .OR. Scanning ) THEN
         IF( CoupledMaxIter > 1 ) THEN
           CALL Info( 'SolveEquations', '-------------------------------------', Level=3 )
           WRITE( Message, * ) 'Coupled system iteration: ', i
           CALL Info( 'SolveEquations', Message, Level=3 )
           CALL Info( 'SolveEquations', '-------------------------------------', Level=3 )         
         END IF
         steadyIt = i
       END IF

       DoneThis = .FALSE.

!      Initialize the mesh output flag to FALSE here, reactivated
!      later for meshes connected to active solvers.
!      ----------------------------------------------------------
       Mesh => Model % Meshes
       DO WHILE( ASSOCIATED( Mesh ) )
          Mesh % OutputActive = .FALSE.
          Mesh => Mesh % Next
       END DO

!------------------------------------------------------------------------------
!      Go trough number of solvers (heat,laminar or turbulent flow, etc...)
!------------------------------------------------------------------------------
       DO k=1,Model % NumberOfSolvers
!------------------------------------------------------------------------------
          Solver => Model % Solvers(k)

          IF ( Solver % PROCEDURE == 0 ) THEN
            IF( .NOT. ( Solver % SolverMode == SOLVER_MODE_COUPLED .OR. &
              Solver % SolverMode == SOLVER_MODE_ASSEMBLY .OR. &
              Solver % SolverMode == SOLVER_MODE_BLOCK ) ) THEN

              CALL Warn('SolveEquations','No routine related to solver!')
              DoneThis(k) = .TRUE.
              CYCLE
            END IF
          END IF

          When = ListGetString( Solver % Values, 'Exec Solver', Found )
          IF ( Found ) THEN
            IF ( When /= 'always' ) THEN
              DoneThis(k) = .TRUE.
              CYCLE
            END IF
          ELSE
            IF ( Solver % SolverExecWhen /= SOLVER_EXEC_ALWAYS ) THEN
              DoneThis(k) = .TRUE.
              CYCLE
            END IF
          END IF

          IF ( AfterConverged(k) .AND. .NOT. ALL(AfterConverged .OR. DoneThis) ) CYCLE
!------------------------------------------------------------------------------

          n = 0
          IF ( ASSOCIATED(Solver % Variable) ) THEN
            IF ( ASSOCIATED(Solver % Variable % Values) ) &
                n = SIZE(Solver % Variable % Values)
            Solver % Variable % PrevNorm = Solver % Variable % Norm
          END IF

          ! There are some operations that require that the previous steady state values
          ! are present. Check for these operations.
          !------------------------------------------------------------------------------
          CalculateDerivative = ListGetLogical( Solver % Values, &
                'Calculate Derivative', Stat )
          IF( CalculateDerivative .AND. .NOT. Scanning ) THEN
            CALL Fatal('SolveEquations','> Calculate Derivative < should be active only for scanning!')
          END IF

          NeedSol = CalculateDerivative
          IF(.NOT. NeedSol ) THEN
            NeedSol = ( ListGetString( Solver % Values, &
                'Steady State Convergence Measure', Stat ) /= 'norm')  
            NeedSol = NeedSol .AND. Stat
          END IF
          IF(.NOT. NeedSol ) THEN
            NeedSol = ListCheckPresent( Solver % Values, &
                'Steady State Relaxation Factor')
          END IF

          IF ( NeedSol .AND. n > 0 ) THEN
            Stat = ASSOCIATED(Solver % Variable % SteadyValues)
            IF(Stat .AND. SIZE(Solver % Variable % SteadyValues) /= n) THEN
              DEALLOCATE(Solver % Variable % SteadyValues)
              Stat = .FALSE.
            END IF
            IF(.NOT. Stat) THEN
              ALLOCATE( Solver % Variable % SteadyValues(n), STAT=istat ) 
              IF ( istat /= 0 ) CALL Fatal( 'SolveEquations', 'Memory allocation error.' )
            END IF
            Solver % Variable % SteadyValues(1:n) = Solver % Variable % Values(1:n)
          END IF

          ! The solver may conditionally be turned into steady state 
          ! if the > Steady State Condition < is set positive.
          !------------------------------------------------------------------------
          TransientSolver = TransientSimulation
          IF( TransientSolver ) THEN
            SSCond = ListGetCReal( Solver % Values,'Steady State Condition',Found )
            IF( Found .AND. SSCond > 0.0_dp ) THEN
              TransientSolver = .FALSE.
              CALL Info('SolveEquations','Running solver in steady state',Level=6)
            END IF
          END IF

          CALL SolverActivate(Model,Solver,dt,TransientSolver)

          IF ( ASSOCIATED(Solver % Variable) ) THEN
            IF ( ASSOCIATED(Solver % Variable % Values) ) &
                n = SIZE(Solver % Variable % Values)
          END IF
!------------------------------------------------------------------------------
!         check for coupled system convergence
!------------------------------------------------------------------------------

           IF ( Scanning .OR. TransientSimulation ) THEN             
             TestConvergence = ( i >= CoupledMinIter .AND. i /= CoupledMaxIter )
           ELSE    ! Steady-state
             TestConvergence = .TRUE.
           END IF
           
           IF( TestConvergence .OR. CalculateDerivative ) THEN
             IF ( ParEnv % PEs > 1 ) THEN
               IF ( ParEnv % Active(ParEnv % MyPE+1) ) THEN
                 CALL ComputeChange(Solver,.TRUE., n)
               ELSE
                 Solver % Variable % SteadyConverged = 1
               END IF
             ELSE
               CALL ComputeChange(Solver,.TRUE.)
             END IF
           END IF
           
           ! The ComputeChange subroutine sets a flag to zero if not yet
           ! converged (otherwise -1/1)
           !------------------------------------------------------------
           IF( TestConvergence ) THEN
             DoneThis(k) = ( Solver % Variable % SteadyConverged /= 0 ) 
           END IF
          
           CALL ParallelAllReduceAnd( DoneThis(k) )
           IF( ALL(DoneThis) ) EXIT
!------------------------------------------------------------------------------
         END DO
!------------------------------------------------------------------------------
       Model % Mesh % Changed = .FALSE.
       IF ( ALL(DoneThis) ) EXIT
    END DO

    IF ( TransientSimulation .AND. .NOT. ALL(DoneThis) ) THEN
       IF ( ListGetLogical( Model % Simulation,  &
               'Coupled System Abort Not Converged', Found ) ) THEN
          CALL Error( 'SolveEquations', ' ' )
          WRITE( Message, * ) 'Coupled system iteration: ', i
          CALL Error( 'SolveEquations', Message )
          CALL Fatal( 'SolveEquations', ' ' )
       ELSE
!         CALL Error( 'SolveEquations', ' ' )
!         WRITE( Message, * ) 'Coupled system iteration: ', i
!         CALL Error( 'SolveEquations', Message )
!         CALL Error( 'SolveEquations', ' ' )
       END IF
    END IF

    END SUBROUTINE SolveCoupled
!------------------------------------------------------------------------------
  END SUBROUTINE SolveEquations
!------------------------------------------------------------------------------
  

!------------------------------------------------------------------------------
!> This is a line of monolithic solvers where different physics and constraints 
!> are assemblied to the same matrix.
!> Provide assembly loop and solution of linear and nonlinear systems
!> This routine uses minimalistic assmebly routines to create the 
!> matrices. Often the results to less labour in coding which may be 
!> comprimized by less flexibility.
! 
! There are currently two modes
! IsCoupledSolver:  the variable consists of several pieces and the matrix equation
!                   is created within this solver.
! IsAssemblySolver: only one assembly routine is used and the matrix equation may
!                   be formed externally using standard procedures.
! TODO:
! non-nodal elements (internal allocation & right elementtype)
! solution (and perm) vectors could be reused by setting pointers to correct positions
! check the time-dependency
! check for eigenmode analysis
! different lumped functionalities (energy,...)
! bandwidth optimization for coupled system
!------------------------------------------------------------------------------
  SUBROUTINE CoupledSolver( Model, Solver, dt, Transient )
!------------------------------------------------------------------------------    
    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Solver_t), TARGET :: Solver
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: PSolver, PSolver2
    TYPE(Variable_t), POINTER :: Var
    TYPE(Element_t),POINTER :: Element
    INTEGER :: i,j,k,l,t,n,nd,NoIterations,iter,istat
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), DAMP(:,:), MASS(:,:), FORCE(:)
    LOGICAL :: GotIt, GotIt2, GotProc, BulkMode, ThisConstraint, &
        IsCoupledSolver, IsAssemblySolver, IsProcedure, IsListMatrix
    INTEGER :: Row, Col, ColDofs, RowDofs, ColInd0, RowInd0, MaxDofs, &
         ColVar, RowVar, Nrow, Ncol, NoVar, NoCons, TotSize, ConDofs, OffSet(20), &
         VarSizes(20),VarDofs(20)
    INTEGER :: ElementsFirst, ElementsLast, bf_id, bc_id, body_id, eq_id
    INTEGER, POINTER :: VarPerm(:), ColPerm(:), RowPerm(:), ColInds(:), RowInds(:), DirPerm(:)
    REAL(KIND=dp), POINTER :: ConsValues(:)
    REAL(KIND=dp) :: NonlinearTol, Norm, PrevNorm, ConsValue, ConsCoeff, ConsVolume
    CHARACTER(LEN=max_name_len) :: str, VarName, ColName, RowName, ConsType
    LOGICAL :: Coupling, Equality, AssemblySymmetric, AssemblyAntisymmetric, &
        AllDirFlag, Robust, ReduceStep
    INTEGER, POINTER :: Rows(:),Cols(:),Indexes(:),AllPerm(:)
    TYPE(ListMatrix_t), POINTER :: Alist(:) => NULL()
    REAL(KIND=dp), POINTER :: ForceVector(:),AllValues(:)
    LOGICAL, POINTER :: AllDir(:)
    TYPE (Matrix_t), POINTER :: Amat
    REAL(KIND=dp), POINTER :: Component(:)
    INTEGER, POINTER :: VarInds(:)
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: SolverParams

#ifndef USE_ISO_C_BINDINGS
    INTERFACE 
      SUBROUTINE ExecLocalAssembly( Proc, Model, Solver, dt, Transient, &
          M, D, S, F, Element, Nrow, Ncol )
        USE Types
        INTEGER(KIND=AddrInt) :: Proc
        TYPE(Model_t)   :: Model
        TYPE(Solver_t)  :: Solver
        REAL(KIND=dp)   :: dt
        LOGICAL :: Transient
        REAL(KIND=dp) :: S(:,:), D(:,:), M(:,:), F(:)
        TYPE(Element_t) :: Element
        INTEGER :: Nrow, Ncol
      END SUBROUTINE ExecLocalAssembly
    END INTERFACE
#endif

    SolverParams => ListGetSolverParams(Solver)

    IsCoupledSolver = .FALSE.
    IsAssemblySolver = .FALSE.

    IsProcedure = ListCheckPresent( SolverParams,'Procedure')

    CALL Info('CoupledSolver','---------------------------------------',Level=5)
    IF( IsProcedure ) THEN
      CALL Info('CoupledSolver','Inheriting matrix from procedure',Level=5)     
    ELSE IF( Solver % SolverMode == SOLVER_MODE_COUPLED ) THEN
      CALL Info('CoupledSolver','Solving a system of equations',Level=5)
      IsCoupledSolver = .TRUE.
    ELSE IF( Solver % SolverMode == SOLVER_MODE_ASSEMBLY ) THEN
      CALL Info('CoupledSolver','Solving one equation',Level=5)      
      IsAssemblySolver = .TRUE.
    ELSE
      CALL Fatal('CoupledSolver','You should maybe not be here?')
    END IF
    CALL Info('CoupledSolver','---------------------------------------',Level=5)


    Mesh => GetMesh()
    PSolver => Solver

    !------------------------------------------------------------------------------
    ! Check out which variables the coupled model includes
    ! and compure size information related to the new coupled dof.
    !------------------------------------------------------------------------------
    Offset = 0
    VarSizes = 0    
    NoVar = 0
    NoCons = 0
    VarDofs = 0

    AllDirFlag = .FALSE.
    IF( IsProcedure .OR. IsAssemblySolver ) THEN
      CALL Info('CoupledSolver','Using existing variable and matrix')
      IF(.NOT. ASSOCIATED(Solver % Matrix)) THEN
        CALL Fatal('CoupledSolver','In legacy solver mode the Matrix should exist!')
      END IF
      IF(.NOT. ASSOCIATED(Solver % Variable)) THEN
        CALL Fatal('CoupledSolver','In legacy solver mode the Variable should exist!')
      END IF
      NoVar = 1
      Var => Solver % Variable
      VarDofs(1) = Var % Dofs
      VarSizes(1) = SIZE( Var % Values )
      TotSize = VarSizes(1)
      MaxDofs = VarDofs(1)
    ELSE

      DO i = 1,100
        WRITE (str,'(A,I0)') 'Variable ',i
        VarName = ListGetString( SolverParams, TRIM(str), GotIt )
        IF(.NOT. GotIt) EXIT
        Var => VariableGet( Mesh % Variables, TRIM(VarName) )
        
        IF(.NOT. ASSOCIATED( Var )) THEN
          CALL Info('CoupledSolver','Variable '//TRIM(VarName)//' does not exist, creating')
          PSolver => Solver
          Var => CreateBlockVariable(PSolver, i, VarName )
        END IF

        NoVar = NoVar + 1
        VarDofs(NoVar) = Var % Dofs
        VarSizes(NoVar) = SIZE( Var % Values )
        Offset(NoVar+1) = Offset(NoVar) + VarSizes(NoVar)
      END DO

      !------------------------------------------------------------------------------------
      ! Here is a hack for taking constraints into account where dofs are created on-the-fly
      ! might be better to create special solvers somewhere else.
      ! The size of constraint is deduced directly from its type.
      !-------------------------------------------------------------------------------------
      DO i = 1,100
        WRITE (str,'(A,I0)') 'Constraint ',i
        VarName = ListGetString( SolverParams, TRIM(str), GotIt )
        IF(.NOT. GotIt) EXIT
        
        NoCons = NoCons + 1      
        
        WRITE (str,'(A,I0,A)') 'Constraint ',i,' Type'
        ConsType = ListGetString( SolverParams,TRIM(str))
        
        SELECT CASE( ConsType )
          
        CASE('integral')
          ConDofs = 1
          
        CASE('floating')
          ConDofs = 1
          AllDirFlag = .TRUE.
          
        CASE('equality')
          ConDofs = 0
          AllDirFlag = .TRUE.
          
        CASE DEFAULT          
          CALL Warn('CoupledSolver','Coupled constraint does not really work yet')
          WRITE (str,'(A,I0,A)') 'Constraint ',i,' DOFs'
          ConDofs = ListGetInteger( SolverParams, str)

        END SELECT

        ! Create the constrained variables for possible other use
        !-----------------------------------------------------------
        IF( ConDofs > 0 ) THEN
          Var => VariableGet( Mesh % Variables, VarName )
          IF( .NOT. ASSOCIATED(Var) ) THEN
            CALL Info('CoupledSolver','Constraint '//TRIM(VarName)//' does not exist, creating')
            ALLOCATE(ConsValues(ConDofs))
            CALL VariableAdd( Mesh % Variables, Mesh, Solver, &
                VarName, ConDofs, ConsValues, Output = .FALSE. )          
            Var => VariableGet( Mesh % Variables, VarName )         
          END IF
        END IF
        
        j = NoVar + NoCons
        VarDofs(j) = ConDofs
        VarSizes(j) = 1
        Offset(j+1) = Offset(j) + VarSizes(j)
      END DO

      DO j=1,NoVar+NoCons
        WRITE(Message,'(A,I0,A,T35,I0)') 'Permutation offset ',j,': ',OffSet(j)
        CALL Info('CoupledSolver',Message)
      END DO
      
      TotSize = SUM( VarSizes )
      MaxDofs = MAXVAL( VarDofs )
      
      WRITE(Message,'(A,T35,I0)') 'Number of coupled variables: ',NoVar
      CALL Info('CoupledSolver',Message)
      
      WRITE(Message,'(A,T35,I0)') 'Number of constraints: ',NoCons
      CALL Info('CoupledSolver',Message)
      
      WRITE(Message,'(A,T35,I0)') 'Size of coupled system: ',TotSize
      CALL Info('CoupledSolver',Message)
      
      ! For the 1st time the matrix should be a list, later CRS
      !--------------------------------------------------------
      IF(.NOT. ASSOCIATED(Solver % Matrix)) THEN
        Amat => AllocateMatrix()
        Amat % ListMatrix => Alist
        Amat % FORMAT = MATRIX_LIST      
        Solver % Matrix => Amat
      END IF

      ! This actual variable is not saved, so a dummy name suffices
      !------------------------------------------------------------
      VarName = ListGetString( SolverParams,'Variable', GotIt )
      IF(.NOT. GotIt) THEN
        CALL Info('CoupledSolver','New coupled variable added: Alldofs', Level=5)
        VarName = 'alldofs'
      END IF

      ! If the variable hasn't been created do it now
      !----------------------------------------------
      Var => VariableGet(Mesh % Variables, VarName )
      IF(.NOT. ASSOCIATED(Var) ) THEN

        ALLOCATE(AllPerm(TotSize),AllValues(TotSize),ForceVector(TotSize))
        DO i=1,TotSize
          AllPerm(i) = i
          AllValues(i) = 0.0_dp
        END DO
        CALL VariableAdd( Mesh % Variables, Mesh, Solver, &
            VarName,1,AllValues,AllPerm,Output=.FALSE.)        
        Amat % Rhs => ForceVector
        Solver % Variable => VariableGet(Mesh % Variables, VarName )
        Var => Solver % Variable
       
        ! Map the original vectors to the monolithic vector if requested
        !----------------------------------------------------------------------------------
        IF( ListGetLogical( SolverParams,'Coupled Initial Guess',GotIt)) THEN
          CALL SingleToCoupledVector()
        END IF
      END IF

    END IF  ! IsProcedure
      

!------------------------------------------------------------------------------
! Do some initial stuff common for both modes
!------------------------------------------------------------------------------
  
    N = Mesh % MaxElementDOFs
    
    ALLOCATE( FORCE( MaxDofs*N ),      &
        STIFF( MaxDofs*N, MaxDofs*N ), &
        DAMP( MaxDofs*N, MaxDofs*N ),  &
        MASS(  MaxDofs*N, MaxDofs*N ), &
        ColInds( N ), RowInds( N ),    &
        Indexes( n ), STAT=istat )
    IF ( istat /= 0 ) CALL FATAL('CoupledSolver','Memory allocation error')
    
    NoIterations = GetInteger( SolverParams,'Nonlinear System Max Iterations',GotIt)
    IF(.NOT. GotIt) NoIterations = 1
    NonlinearTol = GetCReal( SolverParams,'Nonlinear System Convergence Tolerance',gotIt)
    Robust = ListGetLogical(SolverParams,'Nonlinear System Linesearch',GotIt)
    Amat => GetMatrix()
    ForceVector => Amat % Rhs
    IF( AllDirFlag ) ALLOCATE( AllDir(TotSize) ) 
   
!------------------------------------------------------------------------------
! Iterate over any nonlinearity of material or source
!------------------------------------------------------------------------------
    CALL Info('CoupledSolver','-------------------------------------------------',Level=5)

    IsListMatrix = ( Amat % FORMAT == MATRIX_LIST )
    
    DO iter = 1,NoIterations

      WRITE(Message,'(A,T35,I5)') 'Coupled iteration:',iter
      CALL Info('CoupledSolver',Message,Level=5)
   
100   IF( IsProcedure ) THEN

        ! Perform the normal solver procedure with just one iteration 
        ! linear system solver being disabled i.e. just do the assembly.
        ! Also the possible internal linesearch is disabled.
        !---------------------------------------------------------------------
        CALL ListAddInteger( SolverParams,'Nonlinear System Max Iterations',1)
        CALL ListAddLogical( SolverParams,'Linear System Solver Disabled',.TRUE.)
        CALL ListAddLogical( SolverParams,'Skip Compute Nonlinear Change',.TRUE.)
        IF( Robust ) THEN
          CALL ListRemove( SolverParams,'Nonlinear System Linesearch')
        END IF

        CALL SingleSolver( Model, PSolver, dt, Transient )

        CALL ListAddInteger( SolverParams,'Nonlinear System Max Iterations',NoIterations)
        CALL ListRemove( SolverParams,'Linear System Solver Disabled')
        IF(Robust ) THEN
          CALL ListAddLogical(SolverParams,'Nonlinear System Linesearch',.TRUE.)
        ELSE
          CALL ListAddLogical(SolverParams,'Skip Compute Nonlinear Change',.FALSE.)
        END IF          
     
      ELSE
        IF(.NOT. IsListMatrix ) CALL DefaultInitialize()
        IF( AllDirFlag ) AllDir = .FALSE.

!----------------------------------------------------------------
! Here we use the same assembly for coupled system and block system
! approach.
!----------------------------------------------------------------
        DO RowVar = 1,NoVar 
          DO ColVar = 1,NoVar          
            CALL BlockSystemAssembly(PSolver,dt,Transient,RowVar,ColVar,&
                    Offset(RowVar),Offset(ColVar))
          END DO
        END DO
    
!-------------------------------------------------------------------------
! Tailored assembly routines for setting the constraints to the same matrix
!-------------------------------------------------------------------------

        IF( IsCoupledSolver .AND. NoCons > 0 ) THEN
          CALL CoupledConstraintAssembly()
        END IF

!----------------------------------------------------------------------
! The CRS matrix may be created only when the matrix structure is known
! and some initialization may be done only when the matrix exists
!----------------------------------------------------------------------
        IF( IsListMatrix ) THEN
          CALL List_ToCRSMatrix(Amat)
          IsListMatrix = .FALSE.
          CALL AddEquationSolution(PSolver, Transient )
          GOTO 100  
        END IF
        CALL DefaultFinishAssembly()
      
!------------------------------------------------------------------------------
!    Do the Dirichlet conditions using offset
!------------------------------------------------------------------------------          
        IF( IsCoupledSolver ) THEN
          CALL CoupledSystemDirichlet()
        ELSE      
          CALL DefaultDirichletBCs()
        END IF
      END IF

!------------------------------------------------------------------------------
!    Check the stepsize of nonlinear iteration using the Armijo-GoldStein 
!    criterion for the stepsize.
!------------------------------------------------------------------------------          
      IF( Robust  ) THEN   
        IF( CheckStepSize( Solver, iter == 1) ) THEN
          IF( IsCoupledSolver ) CALL CoupledToSingleVector()
          GOTO 100
        ELSE
          IF(Solver % Variable % NonlinConverged == 1) EXIT
        END IF
      END IF

!------------------------------------------------------------------------------
!   Finally solve the system
!------------------------------------------------------------------------------          
      Norm = DefaultSolve()     
      CALL Info('CoupledSolver','-------------------------------------------------',Level=5)

      IF( .NOT. ( IsAssemblySolver .OR. IsProcedure ) ) THEN
        CALL CoupledToSingleVector()
      END IF

      ! The non-robust solver will use one assembly routine less 
      IF(.NOT. Robust) THEN
        IF(Solver % Variable % NonlinChange < NonlinearTol) EXIT
      END IF
    END DO

    
    DEALLOCATE( FORCE, STIFF, DAMP, MASS, ColInds, RowInds, Indexes )
    IF( AllDirFlag ) DEALLOCATE( AllDir ) 
    
    CALL Info('CoupledSolver','All done')
    CALL Info('CoupledSolver','-------------------------------------------------',Level=5)


  CONTAINS 

!------------------------------------------------------------------------------
!> Integration routine for integral type of constraints i.e. body or boundary
!> integral of some dof is known a priori.
!------------------------------------------------------------------------------
   SUBROUTINE IntegralConstraint( Mass, Damp, Stiff, Force, Element, n )
!------------------------------------------------------------------------------
     
      IMPLICIT NONE
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Stiff(:,:), Damp(:,:), Mass(:,:), Force(:)
      TYPE(Element_t), POINTER :: Element
      INTEGER :: n
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      TYPE(GaussIntegrationPoints_t) :: IntegStuff
      REAL(KIND=dp) :: Basis(n)
      REAL(KIND=dp) :: Weight, DetJ
      INTEGER :: i,j,t,p,q
      LOGICAL :: Found
      
      TYPE(Nodes_t),SAVE :: Nodes
      
      CALL GetElementNodes( Nodes ) 
      IntegStuff = GaussPoints( Element )
      
      DO t=1,IntegStuff % n
        Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
            IntegStuff % v(t), IntegStuff % w(t), detJ, Basis )     
        Weight = IntegStuff % s(t) * detJ
        DO p=1,n
          STIFF(1,p) = STIFF(1,p) + Weight * Basis(p)
          FORCE(p) = FORCE(p) + ConsValue * Weight * Basis(p)
        END DO
        ConsVolume = ConsVolume + Weight 
      END DO

!------------------------------------------------------------------------------
    END SUBROUTINE IntegralConstraint
!------------------------------------------------------------------------------
 
   
    !------------------------------------------------------------------------------          
    !> Subroutine sets the Dirichlet conditions for the coupled system 
    !> using the single system vector names and sizes.
    !----------------------------------------------------------------------------------
    SUBROUTINE CoupledSystemDirichlet()
      CALL Info( 'CoupledSolver', 'Setting coupled system Dirichlet conditions', Level=4 )
      
      DO i = 1,NoVar
        WRITE (str,'(A,I0)') 'Variable ',i
        VarName = ListGetString( SolverParams, TRIM(str), GotIt )
        IF(.NOT. GotIt) EXIT
        
        Var => VariableGet( Mesh % Variables, TRIM(VarName) )
        IF (.NOT. ASSOCIATED(Var)) EXIT 
        
        CALL DefaultDirichletBCs( Ux=Var, UOffset=Offset(i) )
      END DO
 
    END SUBROUTINE CoupledSystemDirichlet


    !------------------------------------------------------------------------------          
    ! Subroutine copies results from the single system vectors to a coupled system
    ! initial guess.
    !----------------------------------------------------------------------------------
    SUBROUTINE SingleToCoupledVector()
      CALL Info('CoupledSolver','Copying an initial guess',Level=5)
      
      DO i = 1,NoVar + NoCons
        IF( VarDofs(i) == 0 ) CYCLE
        
        IF( i <= NoVar ) THEN
          WRITE (str,'(A,I0)') 'Variable ',i
        ELSE
          WRITE (str,'(A,I0)') 'Constraint ',i-NoVar        
        END IF
        VarName = ListGetString( SolverParams, TRIM(str), GotIt )
        Var => VariableGet( Mesh % Variables, TRIM(VarName) )
        
        ! CALL Info('CoupledSolver','Copying from variable: '//TRIM(VarName),Level=5)
        DO j=1,SIZE(Var % Values)
          Solver % Variable % Values(Offset(i)+j) = Var % Values(j)
        END DO
      END DO

    END SUBROUTINE SingleToCoupledVector


    !------------------------------------------------------------------------------          
    !> Subroutine copies results from the coupled system vector back to the 
    !> original vectors.
    !----------------------------------------------------------------------------------
    SUBROUTINE CoupledToSingleVector()
      CALL Info('CoupledSolver','Copying results into original variables',Level=5)
      DO i = 1,NoVar + NoCons
        IF( VarDofs(i) == 0 ) CYCLE
        
        IF( i <= NoVar ) THEN
          WRITE (str,'(A,I0)') 'Variable ',i
        ELSE
          WRITE (str,'(A,I0)') 'Constraint ',i-NoVar        
        END IF
        VarName = ListGetString( SolverParams, TRIM(str), GotIt )
        Var => VariableGet( Mesh % Variables, TRIM(VarName) )
        
        ! CALL Info('CoupledSolver','Copying to variable: '//TRIM(VarName),Level=5)
        DO j=1,SIZE(Var % Values)
          Var % Values(j) = Solver % Variable % Values(Offset(i)+j)
        END DO

        CALL InvalidateVariable( Model % Meshes, Solver % Mesh, VarName )
      END DO

    END SUBROUTINE CoupledToSingleVector


    !-----------------------------------------------------------------------------------
    !> Perform constraints assembly
    !> Some constraints are assembled by integration while others are done by elimination.
    !> It is important that constraints are applied after the normal assembly process. 
    !-----------------------------------------------------------------------------------
    SUBROUTINE CoupledConstraintAssembly()
      !-----------------------------------------------------------------------------------
      INTEGER :: TargetDof, TmpInds(4) 

      CALL Info('CoupledSolver','Starting constraint assembly',Level=5)
      
      BulkMode = .TRUE.
200   IF(BulkMode) THEN
        ! CALL Info('CoupledSolver','Starting constraing bulk assembly',Level=5)
        ElementsFirst = 1
        ElementsLast = Mesh % NumberOfBulkElements 
      ELSE
        ! CALL Info('CoupledSolver','Starting boundary assembly',Level=5)
        ElementsFirst = Mesh % NumberOfBulkElements + 1
        ElementsLast =  Mesh % NumberOfBulkElements + &
              Mesh % NumberOfBoundaryElements
      END IF
      
      ! Variables over rows
      !-------------------------------------------
      DO RowVar = 1,NoCons        

        RowDofs = VarDofs(NoVar + RowVar) 
        RowInd0 = Offset(NoVar + RowVar)

        AssemblySymmetric = .FALSE.
        AssemblyAntiSymmetric = .FALSE.

        WRITE (str,'(A,I0)') 'Constraint ',RowVar
        RowName = ListGetString( Solver % Values, TRIM(str), GotIt )

        WRITE (str,'(A,I0,A)') 'Constraint ',RowVar,' Variables'
        VarInds => ListGetIntegerArray( Solver % Values, TRIM(str) )

        IF( ASSOCIATED(VarInds)) THEN
          ColVar =  VarInds(1) 
        ELSE        
          CALL Fatal('CoupledSolver','Cannot continue without pointer to variables')
        END IF
        
        WRITE (str,'(A,I0,A)') 'Constraint ',RowVar,' Components'
        TargetDof = ListGetInteger( Solver % Values, TRIM(str), GotIt )
      
        WRITE (str,'(A,I0,A)') 'Constraint ',RowVar,' Value'
        ConsValue = GetCReal(Solver % Values,TRIM(str), GotIt)

        WRITE (str,'(A,I0,A)') 'Constraint ',RowVar,' Coeff'
        ConsCoeff = GetCReal(Solver % Values,TRIM(str), GotIt)
        IF(.NOT. GotIt) ConsCoeff = 1.0_dp


        IF ( ConsType == 'equality' ) THEN
          WRITE (str,'(A,I0)') 'Variable ',VarInds(2)
          VarName = ListGetString( SolverParams, TRIM(str) )
          Var => VariableGet( Mesh % Variables, TRIM(VarName) ) 
          RowPerm => Var % Perm
          RowDofs = VarDofs(VarInds(2))
          RowInd0 = Offset(VarInds(2)) 
        ELSE
          Nrow = 1
          RowInds(1:Nrow) = 1         
          Var => VariableGet( Mesh % Variables, TRIM(RowName) )      
          
          ! Add the diagonal entry expected by some subroutines
          !-------------------------------------------------------          
          DO i=1,Nrow
            DO j=1,RowDofs
              Row = RowInd0 + RowDofs * (RowInds(i)-1) + j
              CALL AddToMatrixElement( Amat, Row, Row, 0.0_dp )
            END DO
          END DO
        END IF
         
        IF( ConsType == 'integral') THEN
          AssemblySymmetric = .TRUE.
          ConsVolume = 0.0_dp
        END IF

        WRITE (str,'(A,I0)') 'Variable ',ColVar
        ColName = ListGetString( Solver % Values, TRIM(str), GotIt )
        Var => VariableGet( Mesh % Variables, TRIM(ColName) )
          
        ! Constrain only the target variable
        !-----------------------------------------------
        
        ColPerm => Var % Perm
        ColDofs = Var % Dofs
        ColInd0 = Offset(ColVar)
          
        ! The assembly loop for a submatrix starts here
        !------------------------------------------------------------------------------             
        DO t=ElementsFirst,ElementsLast

          Element => Mesh % Elements(t)
          Model % CurrentElement => Element
          
          ! How to treat non-nodal elements must be rethought 
          ! nd = GetElementNOFDOFs( Element, Solver )                  
          !-----------------------------------------------------------------
          n  = GetElementNOFNodes()
          nd = GetElementDOFs(Indexes)
                        
          ! Set the permutations for equality constraint
          !----------------------------------------------------
          IF( ConsType == 'equality') THEN
            Nrow = nd
            RowInds(1:Nrow) = RowPerm(Indexes)
            IF(.NOT. ALL(RowInds(1:Nrow) > 0)) CYCLE
          END IF
            
          Ncol = nd
          ColInds(1:n) = ColPerm(Indexes)
          IF(.NOT. ALL(ColInds(1:Ncol) > 0)) CYCLE                 

            
          ! Check where constraints are active, both bodies and BCs
          !-------------------------------------------------------------------                         
          Coupling = .FALSE.              
          IF( BulkMode ) THEN
            ! Check coupling to bodies using Body/Equation section
            Coupling = GetLogical( GetBodyForce(), RowName, gotIt)
          ELSE 
            Coupling = GetLogical( GetBC(), RowName, gotIt)
          END IF
          IF( .NOT. Coupling ) CYCLE
              
          ! These two constraints are based on moving already assembled information
          ! rather than assembling new information. If the row is already treated cycle.
          ! The constrainst have been tested only with one-component cases.
          !-----------------------------------------------------------------------------
          IF( ConsType == 'floating' .OR. ConsType == 'equality') THEN

            DO i=1,Ncol                
              DO k=0,ColDofs-1
                
                ! Note that for this type the column is rather also a row
                !--------------------------------------------------------
                Col  = ColInd0 + ColDofs * ColInds(i) - k                  
                IF( AllDir(Col) ) CYCLE
                AllDir(Col) = .TRUE.
                
                IF( ConsType == 'floating') THEN
                  Row = RowInd0 + 1
                ELSE IF( ConsType == 'equality') THEN
                  Row = RowInd0 + RowDofs * RowInds(i) - k
                END IF

                IF( IsListMatrix ) THEN
                  CALL MoveRow( Amat, Col, Row )
                  CALL SetMatrixElement( Amat,Col,Col,0.0_dp )
                  CALL SetMatrixElement( Amat,Col,Row,0.0_dp )
                  CYCLE
                END IF

                CALL MoveRow( Amat, Col, Row, ConsCoeff ) 
                ForceVector(Row) = ForceVector(Row) + ForceVector(Col)
                ForceVector(Col) = 0.0_dp
                
                ForceVector(Row) = ForceVector(Row) + ConsValue
                CALL SetMatrixElement( Amat,Col,Col,1.0_dp )
                CALL SetMatrixElement( Amat,Col,Row,-ConsCoeff)
              END DO
            END DO
            CYCLE
          END IF


          ! Do the assembly, now really (active only for some constraints)
          !----------------------------------------------------------------
            
          STIFF = 0.0_dp
          DAMP = 0.0_dp
          MASS = 0.0_dp
          FORCE = 0.0_dp
          
          CALL IntegralConstraint( MASS, DAMP, STIFF, FORCE, Element, Ncol )
          
          IF ( Transient ) THEN
            IF( Solver % TimeOrder == 1 ) THEN
              CALL Default1stOrderTime( MASS,STIFF,FORCE)
            ELSE IF( Solver % TimeOrder == 2) THEN
              CALL Default2ndOrderTime( MASS,DAMP,STIFF,FORCE )
            END IF
          END IF

           
          ! Assemble the matrix with offset
          ! Because we want to have constraints component-wise
          ! There is somewhat dirty hack for calling the gluematrix
          !---------------------------------------------------------
          DO k=1,ColDofs
            IF( TargetDof /= 0 .AND. TargetDof /= k) CYCLE
            TmpInds(1:Ncol) = ColDofs * (ColInds(1:Ncol)-1) + k

            IF( IsListMatrix ) THEN            
              CALL GlueLocalSubMatrix( Amat, &
                  RowInd0,ColInd0,Nrow,Ncol,RowInds,TmpInds,&
                  RowDofs,1,STIFF )
              IF( AssemblySymmetric .OR. AssemblyAntisymmetric ) THEN
                CALL GlueLocalSubMatrix( Amat, &
                    ColInd0,RowInd0,Ncol,Nrow,TmpInds,RowInds,&
                    1,RowDofs,STIFF )               
              END IF
              CYCLE
            END IF

            CALL GlueLocalSubMatrix( Amat, &
                RowInd0,ColInd0,Nrow,Ncol,RowInds,TmpInds,&
                RowDofs,1,STIFF )

            ! For some constraints assemble also the transpose
            !--------------------------------------------------
            IF( AssemblySymmetric ) THEN
              CALL GlueLocalSubMatrix( Amat, &
                  ColInd0,RowInd0,Ncol,Nrow,TmpInds,RowInds,&
                  1,RowDofs,TRANSPOSE(STIFF) )               
            ELSE IF( AssemblyAntisymmetric ) THEN
              CALL GlueLocalSubMatrix( Amat, &
                  ColInd0,RowInd0,Ncol,Nrow,TmpInds,RowInds,&
                  1,RowDofs,-TRANSPOSE(STIFF) )               
            END IF
            
            ! Assemble the r.h.s with offset
            !-----------------------------------------------
            DO i=1,Nrow
              DO j=1,RowDofs
                Row = RowInd0 + RowDofs * (RowInds(i)-1) + j
                ForceVector(Row) = ForceVector(Row) + &
                    FORCE(RowDofs*(i-1)+j)
              END DO
            END DO
          
          END DO          
        END DO
      
        ! For constraints do some special setting
        !---------------------------------------------------
        IF( ConsType == 'integral' ) THEN
          ! PRINT *,'Integral constraint sum of weights: ',ConsVolume
          DO i=1,Nrow
            DO j=1,RowDofs
              Row = RowInd0 + RowDofs * (RowInds(i)-1) + j
              ForceVector(Row) = ConsValue
            END DO
          END DO
        END IF
      END DO
      
      IF(BulkMode) THEN
        ! CALL Info( 'CoupledSolver', 'Bulk assembly done for constraints', Level=4 )
        BulkMode = .FALSE.
        GOTO 200
      ELSE 
        ! CALL Info( 'CoupledSolver', 'Boundary assembly done for constraints', Level=4 )
      END IF

    END SUBROUTINE CoupledConstraintAssembly

!------------------------------------------------------------------------------
  END SUBROUTINE CoupledSolver
!------------------------------------------------------------------------------
 



!------------------------------------------------------------------------------
!> This is a line of solvers where a matrix is matrices and a vector of vectors 
!> are created to allow different kinds of block strategies for the solvers.
!> This strategy has optimal memory consumption even if block strategies are 
!> employed on the linear system level. 
!------------------------------------------------------------------------------
  SUBROUTINE BlockSolver( Model, Solver, dt, Transient )
!------------------------------------------------------------------------------
    IMPLICIT NONE
 !------------------------------------------------------------------------------
    TYPE(Solver_t), TARGET :: Solver
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------

    TYPE(Solver_t), POINTER :: PSolver
    TYPE(Variable_t), POINTER :: Var, SolverVar
    INTEGER :: i,j,k,l,n,nd,NonLinIter,tests,NoTests,iter
    LOGICAL :: GotIt, GotIt2, BlockPrec, BlockGS
    INTEGER :: ColVar, RowVar, NoVar, BlockDofs

    REAL(KIND=dp) :: NonlinearTol, Norm, PrevNorm, Residual, PrevResidual, &
        TotNorm, MaxChange, alpha, beta, omega, rho, oldrho, s, r, PrevTotNorm, &
        Coeff
    CHARACTER(LEN=max_name_len) :: str, VarName, ColName, RowName
    LOGICAL :: Robust, LinearSearch, AcceptStep, IsProcedure, ScaleSystem,&
        ReuseMatrix, LS, InitDone
    INTEGER, POINTER :: VarPerm(:)
    
    TYPE (Matrix_t), POINTER :: Amat, SolverMatrix
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: SolverParams

    CALL Info('BlockSolver','---------------------------------------',Level=5)
    IF( Solver % SolverMode /= SOLVER_MODE_BLOCK ) THEN
      CALL Fatal('BlockSolver','You should maybe not be here?')
    ELSE
      CALL Info('BlockSolver','Solving system of equations utilizing block strategies')
    END IF
    CALL Info('BlockSolver','---------------------------------------',Level=5)

    SolverParams => ListGetSolverParams(Solver)
    Mesh => Solver % Mesh
    PSolver => Solver

    isParallel = ParEnv % PEs > 1
    
         
    ! Determine some parameters related to the block strategy
    !------------------------------------------------------------------------------
    BlockPrec = GetLogical( SolverParams,'Block Preconditioner',GotIt)
    IF(.NOT. GotIt ) BlockPrec = .TRUE.
    BlockGS = GetLogical( SolverParams,'Block Gauss-Seidel',GotIt)    
    IsProcedure = ListCheckPresent( SolverParams,'Procedure')

    ! Initialize the matrix
    ! 1) Matrix is inherited from a monolithic matrix equation, or
    ! 2) Matrix is assemblied using the compact assembly process
    !------------------------------------------------------------------------------
    IF( IsProcedure ) THEN 
      BlockDofs = Solver % Variable % Dofs      

      CALL BlockInitMatrix( Solver, TotMatrix, BlockDofs )

      NoVar = TotMatrix % NoVar

      TotMatrix % Solver => Solver
      SolverMatrix => Solver % Matrix
      SolverVar => Solver % Variable

    ELSE
      TotMatrix => Solver % BlockMatrix
      IF( .NOT.ASSOCIATED(TotMatrix)) THEN
        DO i = 1,9
          WRITE (str,'(A,I0)') 'Variable ',i
          IF(.NOT. ListCheckPresent( SolverParams, TRIM(str)) ) EXIT
        END DO
        NoVar = i-1
        
        CALL BlockInitMatrix( Solver, TotMatrix, NoVar )
        TotMatrix % Solver => Solver
        
        WRITE(Message,'(A,I0)') 'Number of coupled variables: ',NoVar
        CALL Info('BlockSolver',Message)
        IF( NoVar == 0 ) CALL Fatal('BlockSolver','No variables, nothing to do!')

        ! Create the matrix structures using the list structure as 
        !------------------------------------------------------------------------------    
        CALL Info('BlockSolver','Creating matrix structured using list matrices')
        DO RowVar=1,NoVar
          Solver % Variable => TotMatrix % SubVector(RowVar) % Var
          
          DO ColVar=1,NoVar            
            Solver % Variable => TotMatrix % SubVector(ColVar) % Var        
            Solver % Matrix => TotMatrix % Submatrix(RowVar,ColVar) % Mat
            
            Amat => Solver % Matrix        
            Amat % ListMatrix => NULL()
            Amat % FORMAT = MATRIX_LIST      
            Amat % NumberOfRows = 0
            
            CALL BlockSystemAssembly(PSolver,dt,Transient,RowVar,ColVar)
            
            IF( .NOT. ASSOCIATED( Amat % ListMatrix ) ) THEN
              Amat % FORMAT = MATRIX_CRS
            ELSE
              CALL List_ToCRSMatrix(Amat)
              CALL AddEquationSolution(PSolver, Transient )            
            END IF

          END DO
        END DO
      END IF

      ! The variable and solver pointer should direct somewhere as they are used to 
      ! monitor the steady-state solution, for example. Therefore the user may
      ! choose the component if not preferring one. 
      !---------------------------------------------------------------------------
      RowVar = ListGetInteger( SolverParams,'Primary Variable',GotIt)
      IF(.NOT. GotIt) RowVar = 1
      SolverVar => TotMatrix % SubVector(RowVar) % Var
      SolverMatrix => TotMatrix % Submatrix(RowVar,RowVar) % Mat
    END IF

    !------------------------------------------------------------------------------
    ! This is the nonlinear loop that may be used either for true 
    ! nonlinearities.
    !------------------------------------------------------------------------------    
    NonLinIter = GetInteger( SolverParams,'Nonlinear System Max Iterations',GotIt)
    IF(.NOT. GotIt) NonLinIter = 1
    NonlinearTol = GetCReal( SolverParams,'Nonlinear System Convergence Tolerance',gotIt)

    Robust = ListGetLogical(SolverParams,'Nonlinear System Linesearch',GotIt)
    IF( Robust ) THEN
      LinearSearch = ListGetLogical( SolverParams,'Nonlinear System Linesearch Linear')
      NoTests = GetInteger( SolverParams,'Nonlinear System Linesearch Iterations',GotIt)
      IF(.NOT. GotIt) NoTests = NonLinIter
      CALL ListAddString(SolverParams,'Nonlinear System Convergence Measure','residual')
      CALL ListAddLogical(SolverParams,'Skip Compute Nonlinear Change',.TRUE.)
    END IF       

    CALL Info('BlockSolver','-------------------------------------------------',Level=5)
    Residual = -1.0_dp
    PrevResidual = -1.0_dp
    
    DO iter = 1,NonLinIter
      
      WRITE(Message,'(A,T35,I0)') 'Coupled iteration: ',iter
      CALL Info('BlockSolver',Message,Level=5)

      tests = 0

      ! Assembly of matrices either using legacy or compact strategy
      !------------------------------------------------------------------
100   IF( IsProcedure ) THEN

        ! Perform the normal solver procedure with just one iteration 
        ! linear system solver being disabled i.e. just do the assembly.
        !---------------------------------------------------------------------
        CALL ListAddInteger( SolverParams,'Nonlinear System Max Iterations',1)
        CALL ListAddLogical( SolverParams,'Linear System Solver Disabled',.TRUE.)

        CALL SingleSolver( Model, PSolver, dt, Transient )

        ! Scaling of linear system. This could also be fused with the submatrix
        ! picking procedure.
        !---------------------------------------------------------------------
        ScaleSystem = ListGetLogical( SolverParams,'Linear System Scaling',GotIt) 
        IF(.NOT. GotIt) ScaleSystem = .TRUE.
        
        IF( ScaleSystem ) THEN
	  CALL Info('BlockSolver','Applying scaling')
          CALL ScaleLinearSystem(Solver, SolverMatrix, &
              SolverMatrix % Rhs, Solver % Variable % Values )
          CALL ListAddLogical( Solver % Values,'Linear System Scaling',.FALSE.)
        END IF

        CALL BlockPickMatrix( Solver, NoVar )

        CALL ListAddInteger( SolverParams,'Nonlinear System Max Iterations',NonLinIter)
        CALL ListRemove( SolverParams,'Linear System Solver Disabled')
      ELSE
                
        DO RowVar=1,NoVar          
          Solver % Variable => TotMatrix % SubVector(RowVar) % Var          
          DO ColVar=1,NoVar            
            
            Solver % Matrix => TotMatrix % Submatrix(RowVar,ColVar) % Mat
            IF( Solver % Matrix % NumberOfRows == 0 ) CYCLE
            Solver % Variable => TotMatrix % SubVector(ColVar) % Var
            CALL InitializeToZero(Solver % Matrix, Solver % Matrix % rhs)
            
            CALL ListPushNameSpace('block '//TRIM(i2s(RowVar))//TRIM(i2s(ColVar))//':')
            CALL BlockSystemAssembly(PSolver,dt,Transient,RowVar,ColVar)
            
            ! Mainly sets the r.h.s. in transient case correctly
            CALL DefaultFinishAssembly()                    
            
            CALL BlockSystemDirichlet(TotMatrix,RowVar,ColVar)
            CALL ListPopNameSpace()
          END DO
        END DO
      END IF

      ! The user may give a user defined preconditioner matrix
      !-----------------------------------------------------------
      CALL BlockPrecMatrix( Solver, NoVar ) 
      
      IF (isParallel) THEN
        DO RowVar=1,NoVar
          DO ColVar=1,NoVar
            CALL ParallelActive( .TRUE.)
            Amat => TotMatrix % SubMatrix(RowVar,ColVar) % Mat
            Amat % Comm = MPI_COMM_WORLD
            Parenv % ActiveComm = Amat % Comm
            Solver % Variable => TotMatrix % SubVector(ColVar) % Var
            CALL ParallelInitMatrix(Solver,Amat)

            IF(ASSOCIATED(Amat % ParMatrix )) THEN
              Amat % ParMatrix % ParEnv % ActiveComm = &
                       Amat % Comm
              ParEnv = Amat % ParMatrix % ParEnv
            END IF
          END DO
        END DO
     END IF


!------------------------------------------------------------------------------
!    Check the stepsize of nonlinear iteration using the Armijo-GoldStein 
!    criterion for the stepsize, if requested
!------------------------------------------------------------------------------          
      IF( Robust  ) THEN   

200     AcceptStep = CheckStepSizeBlock(TotMatrix,tests==0,PrevResidual,Residual)
        tests = tests + 1
        
        IF( tests >  NoTests ) THEN
          CALL Fatal('BlockSolver','Maximum number of linesearch steps exceeded')
        END IF
        
        ! Update the reference residual only when new step is accepted
        IF( iter == 1 ) THEN
          PrevResidual = Residual
        ELSE
          IF( AcceptStep ) THEN 
            PrevResidual = Residual
            IF(Solver % Variable % NonlinChange < NonlinearTol) EXIT
          ELSE
            IF( LinearSearch ) THEN
              GOTO 100
            ELSE
              GOTO 200
            END IF
          END IF
        END IF
      END IF
      

      !------------------------------------------------------------------------------
      ! Finally solve the system using 'outer: ' as the optional namespace
      ! for the linear system setting.
      !------------------------------------------------------------------------------          
      
      TotNorm = 0.0_dp
      MaxChange = 0.0_dp

      CALL ListPushNameSpace('outer:')

      ! The case with one block is mainly for testing and developing features
      ! related to nonlinearity and assembly.
      !----------------------------------------------------------------------
      IF( NoVar == 1 ) THEN
	CALL Info('BlockSolver','Solving in standard manner',Level=6)

        Solver % Variable => TotMatrix % SubVector(1) % Var
        Solver % Matrix => TotMatrix % Submatrix(1,1) % Mat

        TotNorm = DefaultSolve()
        MaxChange = Solver % Variable % NonlinChange 

      ELSE IF( BlockPrec ) THEN
	CALL Info('BlockSolver','Using block precontioning strategy',Level=6)        
        CALL BlockKrylovIter( Solver, MaxChange )

      ELSE
        Solver % Variable => TotMatrix % SubVector(1) % Var
        Solver % Matrix => TotMatrix % Submatrix(1,1) % Mat

	CALL Info('BlockSolver','Using block solution strategy',Level=6)
        CALL BlockStandardIter( Solver, MaxChange )
      END IF      
      CALL ListPopNameSpace()

      ! For legacy matrices do the backmapping 
      !------------------------------------------
      IF( IsProcedure ) THEN
        Solver % Matrix => SolverMatrix
        Solver % variable => SolverVar
        IF( ScaleSystem ) THEN
          CALL BackScaleLinearSystem(Solver, SolverMatrix, &
              SolverMatrix % Rhs, Solver % Variable % Values )
          CALL ListAddLogical( Solver % Values,'Linear System Scaling',.TRUE.)
        END IF
      END IF
      
      IF(.NOT. Robust ) THEN
        IF( MaxChange < NonlinearTol) EXIT
      END IF
    END DO
       
    ! Before returning set the pointers appropriately as saved before
    !---------------------------------------------------------------------------
    Solver % Variable => SolverVar
    Solver % Matrix => SolverMatrix

    CALL Info('BlockSolver','All done')
    CALL Info('BlockSolver','-------------------------------------------------',Level=5)


  CONTAINS 


    !------------------------------------------------------------------------------          
    !> Subroutine sets the Dirichlet conditions for the block system 
    !> using the single system vector names and sizes. For that aim there is an 
    !> additional flag that is used to detect that the matrix cannot have the digonal 
    !> entry that is by default set to unity and the r.h.s. to the target value.
    !> Now for off-diagonal matrices they will be both omitted. 
    !> The routine assumes that the r.h.s. of the off diagonal matrices is zero.
    !----------------------------------------------------------------------------------
    SUBROUTINE BlockSystemDirichlet(BlockMatrix,NoRow,NoCol)
      
      TYPE(BlockMatrix_t) :: BlockMatrix
      INTEGER :: NoRow, NoCol

      LOGICAL :: OffDiagonal

      CALL Info( 'BlockSolver', 'Setting block system Dirichlet conditions', Level=4 )
      
      Solver % Matrix => BlockMatrix % SubMatrix( NoRow, NoCol ) % Mat
      IF( Solver % Matrix % NumberOfRows == 0 ) RETURN
       
      WRITE (str,'(A,I0)') 'Variable ',NoRow
      VarName = ListGetString( SolverParams, TRIM(str), GotIt )
      IF(.NOT. GotIt) RETURN

      Var => VariableGet( Mesh % Variables, TRIM(VarName) )
      IF (.NOT. ASSOCIATED(Var)) RETURN
 
      Solver % Variable => Var
      OffDiagonal = ( NoRow /= NoCol )
      CALL DefaultDirichletBCs( Ux=Var, OffDiagonalMatrix = OffDiagonal )
 
    END SUBROUTINE BlockSystemDirichlet


    !------------------------------------------------------------------------------          
    !> Compute the block norm, optionally. 
    !----------------------------------------------------------------------------------
    FUNCTION ComputeBlockNorm(BlockMatrix, DiagonalOnly, MatrixOnly ) RESULT ( Residual )
      
      TYPE(BlockMatrix_t), TARGET :: BlockMatrix
      LOGICAL, OPTIONAL :: DiagonalOnly, MatrixOnly
      REAL(KIND=dp) :: Residual

      REAL(KIND=dp) :: rnorm, xnorm, bnorm, TotXnorm, TotRnorm, TotBnorm, TotRelNorm
      REAL(KIND=dp), POINTER :: x(:),res(:),b(:),rtmp(:)
      INTEGER :: n, NoRow,NoCol, NoVar
      TYPE(Matrix_t), POINTER :: A


      CALL Info('BlockSolver','Computing block matrix norm',Level=5)
      
      NoVar = BlockMatrix % NoVar
      ALLOCATE( rtmp( BlockMatrix % MaxSize ), res( BlockMatrix % MaxSize) )

      DO NoRow = 1,NoVar 
        Var => BlockMatrix % SubVector(NoRow) % Var
        n = SIZE( Var % Values )

        x => Var % Values
        xNorm = ComputeNorm(Solver, n, x)

        A => BlockMatrix % SubMatrix(NoRow,NoRow) % Mat
        Solver % Matrix => A
        b => A % rhs
        
        xNorm = ComputeNorm(Solver, n, x)
        bNorm = ComputeNorm(Solver, n, b)

	res = 0.0_dp        
        rtmp = 0.0_dp  
        
        DO NoCol = 1,NoVar           
          IF( PRESENT( DiagonalOnly ) ) THEN
            IF( DiagonalOnly .AND. NoCol /= NoRow ) CYCLE
          END IF
          
          Var => BlockMatrix % SubVector(NoCol) % Var
          x => Var % Values

          A => BlockMatrix % SubMatrix( NoRow, NoCol )  % Mat
          IF( A % NumberOfRows == 0 ) CYCLE
          b => A % rhs
          
          CALL MatrixVectorMultiply( A, x, rtmp)      
          res = res + rtmp
        END DO
        
        IF( PRESENT( MatrixOnly ) ) THEN
          IF( .NOT. MatrixOnly ) THEN              
            res = res - b
          END IF
        ELSE
          res = res - b
        END IF
	
        rNorm = ComputeNorm(Solver, n, res)

        ! PRINT *,'comp norms:',NoRow,Rnorm,Xnorm

        BlockMatrix % SubVector(NoRow) % rnorm = rnorm
        BlockMatrix % SubVector(NoRow) % xnorm = xnorm
        BlockMatrix % SubVector(NoRow) % bnorm = bnorm
      END DO

      TotRnorm = SUM( BlockMatrix % SubVector(1:NoVar) % rnorm ** 2) 
      TotXnorm = SUM( BlockMatrix % SubVector(1:NoVar) % xnorm ** 2) 
      TotBnorm = SUM( BlockMatrix % SubVector(1:NoVar) % bnorm ** 2) 

      TotRnorm = SQRT( TotRnorm / NoVar )
      TotXnorm = SQRT( TotXnorm / NoVar )
      TotBnorm = SQRT( TotBnorm / NoVar )

      BlockMatrix % Rnorm = TotRNorm
      BlockMatrix % Xnorm = TotXnorm
      BlockMatrix % Bnorm = TotBnorm

      ! PRINT *,'tot norms:',TotRnorm,TotXnorm
      DEALLOCATE( rtmp, res )

      Residual = BlockMatrix % Rnorm

    END FUNCTION ComputeBlockNorm
    



!------------------------------------------------------------------------------
    FUNCTION CheckStepSizeBlock(BlockMatrix,FirstTrial,PrevResidual,Residual) &
        RESULT (Success) 
!------------------------------------------------------------------------------
      TYPE(BlockMatrix_t) :: BlockMatrix
      REAL(KIND=dp) :: PrevResidual, Residual
      LOGICAL :: FirstTrial,Success
!------------------------------------------------------------------------------
      INTEGER :: i,n,niter
      TYPE(Matrix_t), POINTER :: A
      REAL(KIND=dp), POINTER :: b(:), x(:), x0(:), r(:)
      REAL(KIND=dp) :: Norm, PrevNorm, rNorm, bNorm, Relaxation, Alpha, Myy
      TYPE(Variable_t), POINTER :: iterV
      LOGICAL :: Stat
      CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
      
      
      SAVE Alpha, Relaxation, Myy
      
      !--------------------------------------------------------------------------
      ! This is the real residual r=b-Ax
      !--------------------------------------------------------------------------
      Residual = ComputeBlockNorm( BlockMatrix ) 
      
      
      ! Negative (impossible) value may be used as indicator that's its the 1st step
      !-----------------------------------------------------------------------------
      IF( PrevResidual < 0.0 ) THEN
        Success = .FALSE.
        RETURN
      END IF
      
      ! At the first step set the relaxation
      !-----------------------------------------------------------------------------
      IF( FirstTrial ) THEN
        Alpha = 1.0_dp
        Relaxation = ListGetConstReal( Solver % Values, &
            'Nonlinear System Linesearch Factor', Stat )
        IF(.NOT. Stat) Relaxation = 0.5_dp
        Myy = ListGetConstReal( Solver % Values, &
            'Nonlinear System Linesearch Limit', Stat )
        IF(.NOT. Stat) Myy = 0.5_dp
      END IF
      
      ! Armijo GoldStein Criterion for accepting stepsize
      !-----------------------------------------------------------------
      Success = ( PrevResidual - Residual > Myy * Alpha * PrevResidual)
      
      IF( Success ) THEN      
        iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
        niter = NINT(iterV % Values(1))
        
        DO i=1,BlockMatrix % NoVar
          Var => BlockMatrix % SubVector(i) % Var
          b => BlockMatrix % SubMatrix(i,i) % Mat % rhs
          x => Var % Values
          n = SIZE( b )
          
          Solver % Variable => Var
          bNorm = ComputeNorm(Solver, n, b)
          
          Var % NonlinChange = Residual / bNorm
          
          Norm = ComputeNorm(Solver, n, x)
          Solver % Variable % Norm = Norm
          
          SolverName = ListGetString( Solver % Values, 'Equation',Stat)
          IF(.NOT. Stat) SolverName = Solver % Variable % Name
          
          WRITE( Message, '(a,g15.8,g15.8,a,I2)') &
              'NS (ITER='//TRIM(i2s(niter))//') (NRM,RELC): (',Norm, Residual / bNorm,&
              ' ) :: '// TRIM(SolverName),i
          CALL Info( 'CheckStepSize', Message, Level=3 )       
        END DO
        iterV % Values(1) = niter + 1 
      ELSE
        
        DO i=1,BlockMatrix % NoVar
          Var => BlockMatrix % SubVector(i) % Var
          IF(.NOT. ASSOCIATED(Var % NonlinValues)) &
              CALL Fatal('CheckStepSize','Previous nonlinear solution is needed')       
          x0 => Var % NonlinValues

          x => Var % Values

          ! PRINT *,'Before Range:',i,MINVAL(x),MAXVAL(x),MINVAL(x0),MAXVAL(x0)

          x = (1-Relaxation) * x0 + Relaxation * x

          ! PRINT *,'After Range:',i,MINVAL(x),MAXVAL(x),MINVAL(x0),MAXVAL(x0)

        END DO
        
        Alpha = Alpha * Relaxation
        CALL Info( 'CheckStepSize','Step rejected, increasing relaxation', Level=3 )      
        ! PRINT *,'Residual',Residual,PrevResidual,Alpha

      END IF

      RowVar = ListGetInteger( SolverParams,'Primary Variable',GotIt)
      IF(.NOT. GotIt) RowVar = 1
      Solver % Matrix => TotMatrix % Submatrix(RowVar,RowVar) % Mat
      Solver % Variable => TotMatrix % SubVector(RowVar) % Var

      
 !------------------------------------------------------------------------------
    END FUNCTION CheckStepSizeBlock
 !------------------------------------------------------------------------------
       
!------------------------------------------------------------------------------
  END SUBROUTINE BlockSolver
!------------------------------------------------------------------------------
 
!---------------------------------------------------
!> Perform assembly for the block system linear
!> system of equations.
!---------------------------------------------------
  SUBROUTINE BlockSystemAssembly(Solver,dt,Transient,RowVar,ColVar,&
      RowIndOffset,ColIndOffset)
!---------------------------------------------------
    TYPE(Solver_t), POINTER :: Solver
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient
    INTEGER :: RowVar, ColVar
    INTEGER, OPTIONAL :: RowIndOffset,ColIndOffset

    INTEGER :: RowInd0, ColInd0
    REAL(KIND=dp) :: SymmCoeff
    TYPE(Element_t), POINTER :: Element
    TYPE(Variable_t), POINTER :: Var
    TYPE(Matrix_t), POINTER :: Amat
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Valuelist_t), POINTER :: SolverParams
    INTEGER :: i,j,t,n,nd,istat
    INTEGER :: ElementsFirst, ElementsLast, Row, Col, RowDofs, ColDofs, Ncol, Nrow
    INTEGER :: AllocCols, AllocRows
    INTEGER, POINTER :: ColPerm(:), RowPerm(:), Indexes(:),ColInds(:),RowInds(:) 
    LOGICAL :: GotIt
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), DAMP(:,:), MASS(:,:), FORCE(:)
    REAL(KIND=dp), POINTER :: ForceVector(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: ProcName, RowName, ColName, str
    INTEGER(KIND=AddrInt) :: ProcPntr    
    LOGICAL :: BulkMode, AssemblySymmetric, AssemblyAntiSymmetric, IsListMatrix
    LOGICAL :: AllocationsDone = .FALSE., Diagonal

#ifndef USE_ISO_C_BINDINGS
    INTERFACE 
      SUBROUTINE ExecLocalAssembly( Proc, Model, Solver, dt, Transient, &
          M, D, S, F, Element, Nrow, Ncol )
        USE Types
        INTEGER(KIND=AddrInt) :: Proc
        TYPE(Model_t)   :: Model
        TYPE(Solver_t)  :: Solver
        REAL(KIND=dp)   :: dt
        LOGICAL :: Transient
        REAL(KIND=dp) :: S(:,:), D(:,:), M(:,:), F(:)
        TYPE(Element_t) :: Element
        INTEGER :: Nrow, Ncol
      END SUBROUTINE ExecLocalAssembly
    END INTERFACE
#endif
    SAVE :: AllocationsDone, AllocCols, AllocRows, &
        FORCE, STIFF, DAMP, MASS, ColInds, RowInds, indexes

    Mesh => Solver % Mesh
    SolverParams => Solver % Values
    Amat => Solver % Matrix
    ForceVector => Amat % rhs

    RowInd0 = 0
    ColInd0 = 0
    IF( PRESENT( RowIndOffset ) ) RowInd0 = RowIndOffset
    IF( PRESENT( ColIndOffset ) ) ColInd0 = ColIndOffset

    ! Row Variable
    !-------------------------------------------
    WRITE (str,'(A,I0)') 'Variable ',RowVar
    RowName = ListGetString( SolverParams, TRIM(str), GotIt )
    Var => VariableGet( Mesh % Variables, TRIM(RowName) ) 
    IF(.NOT. ASSOCIATED( Var ) .AND. RowVar == 1 .AND. ColVar == 1 ) THEN
      Var => Solver % Variable
    END IF
    RowDofs = Var % Dofs
    RowPerm => Var % Perm
    
    ! Column variable
    !------------------------------------------
    IF( ColVar /= RowVar ) THEN
      WRITE (str,'(A,I0)') 'Variable ',ColVar
      ColName = ListGetString( SolverParams, TRIM(str), GotIt )
      Var => VariableGet( Mesh % Variables, TRIM(ColName) )
    END IF          
    ColDofs = Var % Dofs
    ColPerm => Var % Perm

    ! These could be user provided for each block
    !-----------------------------------------
    AssemblySymmetric = ListGetLogical( SolverParams,'Symmetric Assembly',GotIt)
    AssemblyAntiSymmetric = ListGetLogical( SolverParams,'AntiSymmetric Assembly',GotIt)

    IsListMatrix = ( Solver % Matrix % FORMAT == MATRIX_LIST ) 
    N = Mesh % MaxElementDOFs    

    IF( AllocationsDone ) THEN
      IF( RowDofs * n > AllocRows .OR. ColDofs * n > AllocCols )  THEN
        DEALLOCATE( FORCE, STIFF, DAMP, MASS, ColInds, RowInds )        
        AllocationsDone = .FALSE.
      END IF
    END IF
    
    IF( .NOT. AllocationsDone ) THEN
      AllocRows = RowDofs * n
      AllocCols = ColDofs * n
      ALLOCATE(Indexes(AllocRows))

      ALLOCATE( FORCE( AllocRows ),      &
          STIFF( AllocRows, AllocCols ), &
          DAMP( AllocRows, AllocCols ),  &
          MASS( AllocRows, AllocCols ), &
          ColInds( N ), RowInds( N ),    &
          STAT=istat )
      IF ( istat /= 0 ) CALL FATAL('BlockSystemAssembly','Memory allocation error')
      AllocationsDone = .TRUE.
    END IF
      
    CALL Info('BlockSystemAssembly','Starting block system assembly',Level=5)
    
    BulkMode = .TRUE.
    
100 IF(BulkMode) THEN
      ! CALL Info('BlockSystemAssembly','Starting bulk assembly',Level=5)
      ElementsFirst = 1
      ElementsLast = Mesh % NumberOFBulkElements
    ELSE
      ! CALL Info('BlockSystemAssembly','Starting boundary assembly',Level=5)
      ElementsFirst = Mesh % NumberOFBulkElements+1
      ElementsLast = Mesh % NumberOfBulkElements + &
          Mesh % NumberOFBoundaryElements
    END IF
        
      
    ! Load the assembly procudure
    !-----------------------------------------
    IF( BulkMode ) THEN
      WRITE (str,'(A,I0,I0)') 'Bulk Assembly Procedure ',RowVar,ColVar
    ELSE
      WRITE (str,'(A,I0,I0)') 'Boundary Assembly Procedure ',RowVar,ColVar
    END IF    
    ProcName = ListGetString( SolverParams, TRIM(str), GotIt )

    ! Test if the 11 block is not given with 'ij' indexes
    !--------------------------------------------------------
    IF(.NOT. GotIt .AND. RowVar == 1 .AND. ColVar == 1) THEN
      IF( BulkMode ) THEN
        WRITE (str,'(A)') 'Bulk Assembly Procedure'
      ELSE
        WRITE (str,'(A)') 'Boundary Assembly Procedure'
      END IF
      ProcName = ListGetString( SolverParams, TRIM(str), GotIt )
      IF(.NOT. GotIt) THEN
         CALL Fatal('BlockSystemAssembly','Bulk Assembly Precedure not given!')
      END IF
    END IF

    Diagonal = (RowVar == ColVar) 
    
    IF( .NOT. GotIt ) THEN
      IF( BulkMode .AND. Diagonal ) THEN
        CALL Warn('BlockSystemAssembly','Diagonal bulk entries should be assembled!')
      END IF
      RETURN
    END IF

    ProcPntr = GetProcAddr( TRIM(ProcName), abort=.FALSE.)
    IF ( ProcPntr == 0 ) THEN
      CALL Fatal('BlockSystemAssembly','Assembly routine not found: '//TRIM(ProcName))
    ELSE
      CALL Info('BlockSystemAssembly','Using assembly routine: '//TRIM(ProcName))
    END IF
    
    ! These may be fetched within the assembly routine, if needed
    CALL ListAddInteger( SolverParams,'Block Matrix Row',RowVar )
    CALL ListAddInteger( SolverParams,'Block Matrix Column',ColVar )


    ! The assembly loop for a submatrix starts here
    !------------------------------------------------------------------------------             
    DO t=ElementsFirst,ElementsLast
      
      Element => Mesh % Elements(t)
      CurrentModel % CurrentElement => Element
      
      !-----------------------------------------------------------------
      n  = GetElementNOFNodes()
      nd = GetElementDOFs(Indexes)
!     Indexes => Element % NodeIndexes
!     n = Element % TYPE % NumberOfnodes
!     nd = n
      Nrow = nd

      RowInds(1:Nrow) = RowPerm(Indexes(1:nd))
      IF(.NOT. ALL(RowInds(1:Nrow) > 0)) CYCLE

      Ncol = nd
      ColInds(1:Ncol) = ColPerm(Indexes(1:nd))
      IF(.NOT. ALL(ColInds(1:Ncol) > 0)) CYCLE                 

      ! Here just the matrix structure is set, no values
      !---------------------------------------------------------------
      IF( IsListMatrix ) THEN
        CALL GlueLocalSubMatrix( Amat, &
            RowInd0,ColInd0,Nrow,Ncol,RowInds,ColInds,&
            RowDofs,ColDofs,STIFF )
        IF( AssemblySymmetric .OR. AssemblyAntiSymmetric ) THEN
          CALL GlueLocalSubMatrix( Amat, &
              ColInd0,RowInd0,Ncol,Nrow,ColInds,RowInds,&
              ColDofs,RowDofs,STIFF )          
        END IF
        CYCLE
      END IF
      
      ! Do the assembly, now really
      !---------------------------------------------------------------
      STIFF = 0.0_dp
      DAMP = 0.0_dp
      MASS = 0.0_dp
      FORCE = 0.0_dp
      
      CALL ExecLocalAssembly( ProcPntr, CurrentModel, Solver, &
          dt, Transient, MASS, DAMP, STIFF, FORCE, Element, &
          Nrow, Ncol )

      IF ( Transient ) THEN
        IF( Solver % TimeOrder == 1 ) THEN
          CALL Default1stOrderTime( MASS,STIFF,FORCE)
        ELSE IF( Solver % TimeOrder == 2) THEN
          CALL Default2ndOrderTime( MASS,DAMP,STIFF,FORCE )
        END IF
      END IF
      
      IF ( Solver % NOFEigenValues > 0 ) THEN
        IF( Solver % TimeOrder == 1 ) THEN
          CALL DefaultUpdateMass(MASS)
        ELSE IF( Solver % TimeOrder == 2) THEN
          CALL DefaultUpdateDamp(DAMP)
          CALL DefaultUpdateMass(MASS)
        END IF
      END IF

      
      ! Assemble the matrix with offset
      !-----------------------------------------------
      IF( .TRUE. .OR. ColInd0 > 0 .OR. RowInd0 > 0 ) THEN
        CALL GlueLocalSubMatrix( Amat, &
            RowInd0,ColInd0,Nrow,Ncol,RowInds,ColInds,&
            RowDofs,ColDofs,STIFF )
        
        ! Assemble the r.h.s with offset
        !-----------------------------------------------
        DO i=1,Nrow
          DO j=1,RowDofs
            Row = RowInd0 + RowDofs * (RowInds(i)-1) + j
            ForceVector(Row) = ForceVector(Row) + &
                FORCE(RowDofs*(i-1)+j)
          END DO
        END DO
        
        ! For some constraints assemble also the transpose
        !--------------------------------------------------
        IF( AssemblySymmetric ) THEN
          CALL GlueLocalSubMatrix( Amat, &
              ColInd0,RowInd0,Ncol,Nrow,ColInds,RowInds,&
              ColDofs,RowDofs,TRANSPOSE(STIFF) )               
        ELSE IF( AssemblyAntisymmetric ) THEN
          CALL GlueLocalSubMatrix( Amat, &
              ColInd0,RowInd0,Ncol,Nrow,ColInds,RowInds,&
              ColDofs,RowDofs,-TRANSPOSE(STIFF) )               
        END IF
      ELSE              
        CALL DefaultUpdateEquations( STIFF, FORCE )
      END IF
      
    END DO
    
    IF(BulkMode) THEN
      ! CALL Info( 'BlockSystemAssembly', 'Bulk assembly done for blocks', Level=4 )
      BulkMode = .FALSE.
!      GOTO 100
    ELSE 
      ! CALL Info( 'BlockSystemAssembly', 'Boundary assembly done for blocks', Level=4 )
    END IF

  END SUBROUTINE BlockSystemAssembly
!-----------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> This executes the original line of solvers (legacy solvers) where each solver 
!> includes looping over elements and the convergence control. From generality
!> point of view this misses some opportunities to have control of the nonlinear
!> system. 
!------------------------------------------------------------------------------
  SUBROUTINE SingleSolver( Model, Solver, dt, TransientSimulation )
!------------------------------------------------------------------------------
     TYPE(Model_t)  :: Model
     TYPE(Solver_t),POINTER :: Solver
     LOGICAL :: TransientSimulation
     REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
     LOGICAL :: stat, Found, GB, MeActive, GotIt
     INTEGER :: i, j, k, l, col, row, n, SolverAddr, BDOFs, maxdim, dsize, size0
     TYPE(Element_t), POINTER :: CurrentElement
     TYPE(ValueList_t), POINTER :: SolverParams
     CHARACTER(LEN=MAX_NAME_LEN) :: EquationName, str, Message, Equation

     INTEGER, ALLOCATABLE :: memb(:)
     TYPE(Matrix_t), POINTER :: M
     INTEGER :: comm_active, group_active, group_world, ierr

     LOGICAL :: ApplyMortar, FoundMortar, SlaveNotParallel
     TYPE(Matrix_t), POINTER :: CM, CM0, CM1, CMP

!------------------------------------------------------------------------------
     IF ( Solver % Mesh % Changed .OR. Solver % NumberOfActiveElements <= 0 ) THEN
       Solver % NumberOFActiveElements = 0
       EquationName = ListGetString( Solver % Values, 'Equation', Found)

       IF ( Found ) THEN
          IF ( ASSOCIATED(Solver % ActiveElements)) DEALLOCATE( Solver % ActiveElements )
          ALLOCATE( Solver % ActiveElements( Solver % Mesh % NumberOfBulkElements + &
                       Solver % Mesh % NumberOFBoundaryElements ) )

          Maxdim = 0
          DO i=1,Solver % Mesh % NumberOfBulkElements+Solver % Mesh % NumberOFBoundaryElements
             CurrentElement => Solver % Mesh % Elements(i)
             IF( CurrentElement % PartIndex /= ParEnv % myPE ) CYCLE
             IF ( CheckElementEquation( Model, CurrentElement, EquationName ) ) THEN
                Solver % NumberOfActiveElements = Solver % NumberOFActiveElements + 1
                Solver % ActiveElements( Solver % NumberOFActiveElements ) = i
                Maxdim = MAX( CurrentElement % TYPE % DIMENSION, Maxdim )
             END IF
          END DO
          CALL ListAddInteger( Solver % Values, 'Active Mesh Dimension', Maxdim )

          IF( ListGetLogical( Solver % Values,'Calculate Weights',Found )) THEN
            CALL CalculateNodalWeights(Solver,.FALSE.)
          END IF
          IF( ListGetLogical( Solver % Values,'Calculate Boundary Weights',Found )) THEN
            CALL CalculateNodalWeights(Solver,.TRUE.) 
          END IF
       END IF
     END IF
!------------------------------------------------------------------------------

     SlaveNotParallel = ListGetLogical( Solver % Values, 'Slave not parallel',Found )

     MeActive = ASSOCIATED(Solver % Matrix)
     IF ( MeActive ) &
        MeActive = MeActive .AND. (Solver % Matrix % NumberOfRows > 0)
     IF(.NOT.SlaveNotParallel) CALL ParallelActive( MeActive )

     IF ( ParEnv % PEs > 1 .AND. .NOT.SlaveNotParallel ) THEN
       ! Check that the solver is active in some of the active output solvers
       IF( ANY( ParEnv % Active(MinOutputPE+1:MaxOutputPE+1) ) ) THEN
         IF( ParEnv % MyPe >= MinOutputPE .AND. &
             ParEnv % MyPe <= MaxOutputPE ) THEN 
           OutputPE = ParEnv % MyPE
         ELSE
           OutputPE = -1
         END IF
       ELSE         
        ! Otherwise get the first active partition for this solver
        DO i=1,ParEnv % PEs
           IF ( ParEnv % Active(i) ) THEN
             EXIT
           END IF
         END DO

         OutputPE = -1
         IF ( i-1 == ParEnv % MyPE ) THEN
           OutputPE = i-1 
         ELSE IF( i > ParEnv % PEs .AND. ParEnv % myPE == 0 ) THEN
           OutputPE = 0
         END IF
       END IF


       n = COUNT(ParEnv % Active)
       IF ( n>0 .AND. n<ParEnv % PEs ) THEN
         IF ( ASSOCIATED(Solver % Matrix) ) THEN
           IF ( Solver % Matrix % Comm /= MPI_COMM_WORLD ) &
              CALL MPI_Comm_Free( Solver % Matrix % Comm, ierr )
         END IF

         CALL MPI_Comm_group( MPI_COMM_WORLD, group_world, ierr )
         ALLOCATE(memb(n))
         n = 0
         DO i=1,ParEnv % PEs
           IF ( ParEnv % Active(i) ) THEN
             n=n+1
             memb(n)=i-1
           END IF
         END DO
         CALL MPI_Group_incl( group_world, n, memb, group_active, ierr)
         DEALLOCATE(memb)
         CALL MPI_Comm_create( MPI_COMM_WORLD, group_active, &
                 comm_active, ierr)

         M => Solver % Matrix
         DO WHILE(ASSOCIATED(M))
           M % Comm = comm_active
           M => M % Parent
         END DO
       ELSE
         M => Solver % Matrix
         DO WHILE( ASSOCIATED(M) )
           M % Comm = MPI_COMM_WORLD
           M => M % Parent
         END DO
       END IF
     END IF

     IF ( ASSOCIATED(Solver % Matrix) ) THEN
       ParEnv % ActiveComm = Solver % Matrix % Comm
       IF ( ParEnv % PEs>1 .AND. MeActive ) THEN
         IF ( ASSOCIATED(Solver % Mesh % ParallelInfo % INTERFACE) ) THEN
           IF (.NOT. ASSOCIATED(Solver % Matrix % ParMatrix) ) &
             CALL ParallelInitMatrix(Solver, Solver % Matrix )

           Solver % Matrix % ParMatrix % ParEnv % ActiveComm = &
                    Solver % Matrix % Comm
           ParEnv = Solver % Matrix % ParMatrix % ParEnv
         END IF
       END IF
     ELSE IF (.NOT.SlaveNotParallel) THEN
       Parenv % ActiveComm = MPI_COMM_WORLD
     END IF


     ! Linear constraints from mortar BCs:
     ! -----------------------------------
!     FoundMortar = .FALSE.
!     IF(.NOT.GetLogical(ListGetSolverParams(),'Mortar Projector Nonlinear',Found)) &
!     FoundMortar = 
     CALL GenerateProjectors(Model,Solver,Nonlinear = .FALSE. )

     CALL INFO("SingleSolver", "Attempting to call solver", level=5)
     SolverParams => ListGetSolverParams(Solver)
     Equation = GetString(SolverParams, 'Equation', GotIt)
     IF (GotIt) THEN
        WRITE(Message,'(A,A)') 'Solver Equation string is: ', TRIM(Equation)
        CALL INFO("SingleSolver", Message, level=5)
     END IF
#ifdef SGIn32
     SolverAddr = Solver % PROCEDURE
     CALL ExecSolver( SolverAddr, Model, Solver, dt, TransientSimulation)
#else
     CALL ExecSolver( Solver % PROCEDURE, &
         Model, Solver, dt, TransientSimulation)
#endif

!    IF(FoundMortar) CALL ReleaseProjectors(Solver)
!------------------------------------------------------------------------------
   END SUBROUTINE SingleSolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Runs a solver as defined in the command file. There are several
!> ways how to skip the execution. The are also three different ways
!> how the matrices may be assembled and solved: standard (single), coupled and
!> block. 
!------------------------------------------------------------------------------
  SUBROUTINE SolverActivate( Model, Solver, dt, TransientSimulation )
!------------------------------------------------------------------------------
     TYPE(Model_t)  :: Model
     TYPE(Solver_t),POINTER :: Solver
     LOGICAL :: TransientSimulation
     REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: OrigDT, DTScal
     LOGICAL :: stat, Found, TimeDerivativeActive, Timing, IsPassiveBC, &
         UpdateExported, GotCoordTransform, NamespaceFound
     INTEGER :: i, j, n, BDOFs, timestep, timei, timei0, PassiveBcId
     INTEGER, POINTER :: ExecIntervals(:),ExecIntervalsOffset(:)
     REAL(KIND=dp) :: tcond, t0, rt0, st, rst, ct
#ifndef USE_ISO_C_BINDINGS
     CPUTime,RealTime
#endif
     TYPE(Variable_t), POINTER :: TimeVar, IterV
     CHARACTER(LEN=MAX_NAME_LEN) :: str, CoordTransform
     TYPE(ValueList_t), POINTER :: Params

     SAVE TimeVar
!------------------------------------------------------------------------------
     CALL SetCurrentMesh( Model, Solver % Mesh )
     Model % Solver => Solver
     Params => ListGetSolverParams(Solver)

     CoordTransform = ListGetString(Params,'Coordinate Transformation',&
         GotCoordTransform )
     IF( GotCoordTransform ) THEN
       CALL CoordinateTransformation( Solver % Mesh, CoordTransform, &
           Params, .FALSE. )
     END IF

     ! Some properties might be rotated by "Property Rotate"
     !--------------------------------------------------------------
     CALL SetRotatedProperties(Model, Solver % Mesh)

     ! Some keywords may be normalized by area or volume
     !--------------------------------------------------------------
     CALL SetNormalizedKeywords(Model, Solver % Mesh)

!------------------------------------------------------------------------------
! The solver may be skipped if the exec condition is negative. 
! This allows for complicated conditions, and the earlier simple ones 
! have therefore been replaced by this keyword.
!------------------------------------------------------------------------------
     IF( ListCheckPresent( Params,'Start Time') ) &
       CALL Fatal('SolverActivate','Use > Exec Condition = Real < instead of > Start Time <')

     IF( ListCheckPresent( Params,'Stop Time') ) & 
       CALL Fatal('SolverActivate','Use > Exec Condition = Real < instead of > Stop Time <')

     st = ListGetCReal( Params,'Exec Condition',Found) 
     IF( Found .AND. st < 0.0 ) RETURN

!---------------------------------------------------------------------------------   
! There may also be predefined discrete intervals for the execution of the solver.
!---------------------------------------------------------------------------------   
     ExecIntervals => ListGetIntegerArray( Params,'Exec Intervals', Found )
     IF( .NOT. Found ) THEN
       ExecIntervals =>  ListGetIntegerArray( Params,'Exec Interval', Found )
     END IF
     IF ( Found ) THEN
       TimeVar => VariableGet( Model % Variables, 'Timestep Interval' )
       timei = NINT(Timevar % Values(1))

       IF( ExecIntervals(timei) == 0 ) RETURN

       TimeVar => VariableGet( Model % Variables, 'Timestep' )
       timestep = NINT(TimeVar % Values(1))

       ExecIntervalsOffset =>  ListGetIntegerArray( Params,&
           'Exec Intervals Offset', Found )
       IF( Found ) THEN
         timei0 = ExecIntervalsOffset(timei)
       ELSE
         timei0 = 0
       END IF

       IF( MOD( timestep-1-timei0, ExecIntervals(timei)) /= 0 ) RETURN
     END IF

!------------------------------------------------------------------------------
! If solver timing is requested start the watches
!------------------------------------------------------------------------------
     Timing = ListGetLogical(Params,'Solver Timing',Found)
     IF( Timing ) THEN
       t0 = CPUTime()
       rt0 = RealTime()
     END IF

     Solver % Mesh % OutputActive = .TRUE.
     TimeDerivativeActive = TransientSimulation

!----------------------------------------------------------------------
! This is to avoid resetting of certain info that could be interesting
! when saving data i.e. using an auxiliary solver.
!----------------------------------------------------------------------
     IF(.NOT. ListGetLogical( Params,'Auxiliary Solver',Found)) THEN
       DTScal = ListGetConstReal( Params, 'Timestep Scale', Found )
       IF ( .NOT. Found ) DTScal = 1.0_dp
       Solver % dt = DtScal * dt 

       IF ( TransientSimulation ) THEN
         TimeDerivativeActive = &
           ListGetLogical( Params, 'Time Derivative Active', Found )

         IF ( .NOT. Found ) THEN
           TimeDerivativeActive = .TRUE.
           tcond = ListGetCReal(Params,'Time Derivative Condition',Found)
           IF ( Found ) TimeDerivativeActive = TimeDerivativeActive .AND. tcond>0
         END IF
       END IF

       iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
       iterV % Values(1) = 1

       str = ListGetString( Params, 'Namespace', NamespaceFound )
       IF (NamespaceFound) CALL ListPushNamespace(TRIM(str))
     END IF

 !------------------------------------------------------------------------------
 ! Check for passive-active boundaries
 !------------------------------------------------------------------------------
     PassiveBcId = 0
     IsPassiveBC = .FALSE.
     DO j=1,Model % NumberOfBCs
       IsPassiveBC = ListGetLogical( Model % BCs(j) % Values, &
           'Passive Target',Found)
       IF (IsPassiveBC) THEN
         PassiveBcId = j
         EXIT
       END IF
     END DO
     IF ( IsPassiveBC ) THEN
       CALL GetPassiveBoundary( Model, Model % Mesh, PassiveBcId )
       WRITE(Message, '(A,I0,A,I0,A)' ) &
           'Passive element BC no. ',j, ' assigned to BC-ID no. ', &
           PassiveBcId
       CALL INFO('MainUtils',Message,Level=5)
     END IF

 !---------------------------------------------------------------------
 ! Call the correct type of solver: standard (single), coupled or block
 ! This is where everything really happens!
 !---------------------------------------------------------------------
     IF( Solver % SolverMode == SOLVER_MODE_COUPLED .OR. &
         Solver % SolverMode == SOLVER_MODE_ASSEMBLY ) THEN
       CALL CoupledSolver( Model, Solver, DTScal * dt, TimeDerivativeActive )
     ELSE IF( Solver % SolverMode == SOLVER_MODE_BLOCK ) THEN
       CALL BlockSolver( Model, Solver, DTScal * dt, TimeDerivativeActive )
     ELSE 
       CALL SingleSolver( Model, Solver, DTScal * dt, TimeDerivativeActive )
     END IF
     IF(.NOT. ListGetLogical( Params,'Auxiliary Solver',Found)) THEN
       IF(NamespaceFound) CALL ListPopNamespace()
     END IF
     Solver % dt = dt
     Solver % TimesVisited = Solver % TimesVisited + 1

     IF( GotCoordTransform ) THEN
       CALL BackCoordinateTransformation( Solver % Mesh )
     END IF

!------------------------------------------------------------------------------
! After solution register the timing, if requested
!------------------------------------------------------------------------------
     IF( Timing ) THEN
       st  = CPUTime() - t0
       rst = RealTime() - rt0

       str = ListGetString( Params,'Equation',Found)
       CALL ListAddConstReal(CurrentModel % Simulation,'res: solver cpu time '&
              //TRIM(str),st)
       CALL ListAddConstReal(CurrentModel % Simulation,'res: solver real time '&
              //TRIM(str),rst)
       WRITE(Message,'(a,f8.2,f8.2,a)') 'Solver time (CPU,REAL) for '&
        //TRIM(str)//': ',st,rst,' (s)'
       CALL Info('SolverActivate',Message)    

       IF( ListGetLogical(Params,'Solver Timing Cumulative',Found)) THEN
          ct = ListGetConstReal(CurrentModel % Simulation,'res: cum solver cpu time '&
                //TRIM(str),Found)
          st = st + ct
          ct = ListGetConstReal(CurrentModel % Simulation,'res: cum solver real time '&
                //TRIM(str),Found)
          rst = rst + ct
          CALL ListAddConstReal(CurrentModel % Simulation,'res: cum solver cpu time '&
              //TRIM(str),st)
          CALL ListAddConstReal(CurrentModel % Simulation,'res: cum solver real time '&
              //TRIM(str),rst)
        END IF 
      END IF

!------------------------------------------------------------------------------
   END SUBROUTINE SolverActivate
!------------------------------------------------------------------------------

END MODULE MainUtils

!> \} ElmerLib
