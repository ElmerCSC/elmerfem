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
! *  Module for computing the wire direction for coils
! *
! *  Author: Eelis Takala
! *  Email:   eelis.takala@trafotek.fi
! *  Web:     http://www.trafotek.fi
! *  Address: Trafotek
! *           Kaarinantie 700
! *           20540 Turku
! *
! *  Original Date: December 2015
! *
! *****************************************************************************/
 
!> \ingroup Solvers
!> \{
!------------------------------------------------------------------------------
SUBROUTINE Wsolve_Init0(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  CHARACTER(LEN=MAX_NAME_LEN) :: varname, jvarname
  LOGICAL :: Found
  INTEGER, POINTER :: Active(:)
  INTEGER :: mysolver,i,j,k,l,n,m,vDOFs, soln
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Solver_t), POINTER :: Solvers(:), PSolver
  INTEGER, SAVE :: visited=0

  visited = visited + 1

  ! This is really using DG so we don't need to make any dirty tricks to create DG fields
  ! as is done in this initialization. 
  IF (GetLogical(GetSolverParams(),'Discontinuous Galerkin',Found)) RETURN

  PSolver => Solver
  DO mysolver=1,Model % NumberOfSolvers
    IF ( ASSOCIATED(PSolver,Model % Solvers(mysolver)) ) EXIT
  END DO

  varname = GetString(GetSolverParams(), 'Variable', Found)
  IF (.NOT. Found) varname = 'W'
  

  n = Model % NumberOfSolvers
  DO i=1,Model % NumberOFEquations
    Active => ListGetIntegerArray(Model % Equations(i) % Values, &
                'Active Solvers', Found)
    m = SIZE(Active)
    IF ( ANY(Active==mysolver) ) &
      CALL ListAddIntegerArray( Model % Equations(i) % Values,  &
           'Active Solvers', m+1, [Active, n+1] )
  END DO

  ! Create DG solver structures on-the-fly without actually solving the matrix
  ! equations. It is assumed that the DG field within each element is independent
  ! and hence no coupling between elemental fields is needed. 
  ALLOCATE(Solvers(n+1))
  Solvers(1:n) = Model % Solvers
  Solvers(n+1) % Values => ListAllocate()
  SolverParams => Solvers(n+1) % Values
  CALL ListAddLogical( SolverParams, 'Discontinuous Galerkin', .TRUE. )
  Solvers(n+1) % DG = .TRUE.
  Solvers(n+1) % PROCEDURE = 0
  Solvers(n+1) % ActiveElements => NULL()
  CALL ListAddString( SolverParams, 'Exec Solver', 'never' )
  CALL ListAddLogical( SolverParams, 'No Matrix',.TRUE.)
  CALL ListAddLogical( SolverParams, 'Optimize Bandwidth',.FALSE.)
  CALL ListAddString( SolverParams, 'Equation', &
  'elementaladd'//i2s(visited) )
  CALL ListAddString( SolverParams, 'Procedure', &
              'WPotentialSolver Wsolve_Dummy',.FALSE. )
  CALL ListAddString( SolverParams, 'Variable', '-nooutput '//TRIM(varname)//'_dummy' )


!  pname = ListGetString( Model % Solvers(soln) % Values, 'Mesh', Found )
!  IF(Found) THEN
!    CALL ListAddString( SolverParams, 'Mesh', pname )
!  END IF

  i = 1
  DO WHILE(.TRUE.)
    IF(ListCheckPresent(SolverParams, "Exported Variable "//i2s(i))) THEN
      i=i+1
    ELSE
      EXIT
    END IF
  END DO

  CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
        TRIM(varname)//" Potential" )

  i=i+1

  CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
          "-dg "//TRIM(varname)//" Vector E["//TRIM(varname)//" Vector E:3]" )

  IF (GetLogical(GetSolverParams(), 'Compute J Vector', Found)) THEN
    jvarname = GetString(GetSolverParams(), 'J Variable', Found)
    IF (.NOT. Found) jvarname = 'J'
    i=i+1

    CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
            "-dg "//TRIM(jvarname)//" Vector E["//TRIM(jvarname)//" Vector E:3]" )
  END IF

  DEALLOCATE(Model % Solvers)
  Model % Solvers => Solvers
  Model % NumberOfSolvers = n+1
!------------------------------------------------------------------------------
END SUBROUTINE Wsolve_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE Wsolve_Dummy(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
END SUBROUTINE Wsolve_Dummy
!------------------------------------------------------------------------------

SUBROUTINE Wsolve( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the Poisson equation!
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
  USE LinearAlgebra
  USE CircuitUtils
  USE MGDynMaterialUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: AllocationsDone = .FALSE., Found, PosEl, NegEl, NoRotM=.False.
  TYPE(Element_t),POINTER :: Element

  REAL(KIND=dp) :: Norm, Wnorms(Model%numberofbodies)
  INTEGER :: n, nb, nd, t, istat, active, i, j, MaxN
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: BodyForce, BC, BodyParams, Material
  TYPE(ValueList_t), POINTER :: CompParams
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
                                RotM(:,:,:), Tcoef(:,:,:)
  CHARACTER(LEN=MAX_NAME_LEN):: CoilType
  LOGICAL :: CoilBody


  SAVE STIFF, LOAD, FORCE, Tcoef, RotM, AllocationsDone
!------------------------------------------------------------------------------

  IF (.NOT.ASSOCIATED(Solver % Matrix)) RETURN
  CALL AddComponentsToBodyLists()
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()

  IF ( .NOT. AllocationsDone ) THEN
     N = Solver % Mesh % MaxElementDOFs  ! just big enough for elemental arrays
     ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N), Tcoef(3,3,N), RotM(3,3,N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'Wsolve', 'Memory allocation error.' )
     END IF

     MaxN = N
     AllocationsDone = .TRUE.
  END IF

   !System assembly:
   !----------------
   Active = GetNOFActive()
   CALL DefaultInitialize()
   DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()

      IF (SIZE(Tcoef,3) /= n) THEN
        DEALLOCATE(Tcoef,RotM)
        ALLOCATE(Tcoef(3,3,N), RotM(3,3,N), STAT=istat)
        IF ( istat /= 0 ) THEN
          CALL Fatal( 'Wsolve', 'Memory allocation error.' )
        END IF
      END IF
      
      LOAD = 0.0d0
      BodyForce => GetBodyForce()
      IF ( ASSOCIATED(BodyForce) ) &
         Load(1:n) = GetReal( BodyForce, 'Source', Found )

     CoilBody = .FALSE.
     CompParams => GetComponentParams( Element )
     CoilType = ''
     RotM = 0._dp
     IF (ASSOCIATED(CompParams)) THEN
       CoilType = GetString(CompParams, 'Coil Type', Found)
       IF (Found) CoilBody = .TRUE.
     END IF 

      IF (CoilBody) THEN
        SELECT CASE (CoilType)
        CASE ('stranded')
          CoilBody = .True.
          NoRotM = GetLogical(CompParams, 'Without RotM', Found)
          IF (.NOT. NoRotM) CALL GetElementRotM(Element, RotM, n)
        CASE ('massive')
          CoilBody = .True.
        CASE ('foil winding')
          CoilBody = .True.
          CALL GetElementRotM(Element, RotM, n)
        CASE DEFAULT
          CALL Fatal ('Wsolve', 'Non existent Coil Type Chosen!')
        END SELECT
      END IF

      Tcoef = 0.0d0
      Tcoef = GetElectricConductivityTensor(Element,n,'re',CoilBody,CoilType)

      !Get element local matrix and rhs vector:
      !----------------------------------------
      CALL LocalMatrix(  STIFF, FORCE, LOAD, Element, CoilBody, CoilType, Tcoef, RotM, n, nd+nb, NoRotM)
      CALL CondensateP( nd, nb, STIFF, FORCE )
      CALL DefaultUpdateEquations( STIFF, FORCE )
   END DO

   DO t=1, Solver % Mesh % NumberOfBoundaryElements
     CYCLE
     ! get element and BC info
     ! -----------------------
     Element => GetBoundaryElement(t)
     IF ( .NOT.ActiveBoundaryElement() ) CYCLE
     n = GetElementNOFNodes()
     ! no evaluation of Neumann BC’s on points
     IF ( GetElementFamily() == 1 ) CYCLE
     BC => GetBC()
     FORCE = 0.0d00
     STIFF = 0.0d00

     PosEl = .False.
     NegEl = .False.
     PosEl = GetLogical(BC,'Positive Electrode', Found)
     NegEl = GetLogical(BC,'Negative Electrode', Found)

     LOAD = 0._dp
     IF (PosEl) LOAD = 1._dp
     IF (NegEl) LOAD = -1._dp
     
!     IF (Solver % Variable % name == 'w') CALL BoundaryCondition(LOAD, FORCE, Element, n)

     CALL DefaultUpdateEquations( STIFF, FORCE )
   END DO

   CALL DefaultFinishAssembly()
   CALL DefaultDirichletBCs()

   ! And finally, solve:
   !--------------------
   Norm = DefaultSolve()

  ! CALL SaveWPotSolution(Model % numberofbodies, Tcoef)

  Wnorms = GetWNormsForBodies(Model % numberofbodies)
  DO t=1,Active
     Element => GetActiveElement(t)
     n = GetElementNOFNodes()

     CoilBody = .FALSE.
     CompParams => GetComponentParams( Element )
     CoilType = ''
     RotM = 0._dp
     IF (ASSOCIATED(CompParams)) THEN
       CoilType = GetString(CompParams, 'Coil Type', Found)
       IF (Found) CoilBody = .TRUE.
     END IF 

      IF (CoilBody) THEN
        SELECT CASE (CoilType)
        CASE ('stranded')
          CoilBody = .True.
          NoRotM = GetLogical(CompParams, 'Without RotM', Found)
          IF (.NOT. NoRotM) CALL GetElementRotM(Element, RotM, n)
        CASE ('massive')
          CoilBody = .True.
        CASE ('foil winding')
          CoilBody = .True.
          CALL GetElementRotM(Element, RotM, n)
        CASE DEFAULT
          CALL Fatal ('Wsolve', 'Non existent Coil Type Chosen!')
        END SELECT
      END IF

      Tcoef = 0.0d0
      Tcoef = GetElectricConductivityTensor(Element,n,'re',CoilBody,CoilType)
      CALL SaveElementWSolution(Element, n, Wnorms(Element%BodyId), RotM, Tcoef, NoRotM)

  END DO
CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, LOAD, Element, CoilBody, CoilType, Tcoef, RotM, n, nd, NoRotM )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:), Tcoef(:,:,:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP, C(3,3), &
                     RotMLoc(3,3), RotM(3,3,n)
    CHARACTER(LEN=MAX_NAME_LEN):: CoilType
    LOGICAL :: Stat, CoilBody, NoRotM
    INTEGER :: i,j,t
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    C(1,1:3) = [0,0,0]
    C(2,1:3) = [0,0,0]
    C(3,1:3) = [0,0,1]

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
       IP % W(t), detJ, Basis, dBasisdx )

      ! Compute the conductivity tensor
      ! -------------------------------
      DO i=1,3
        DO j=1,3
          C(i,j) = SUM( Tcoef(i,j,1:n) * Basis(1:n) )
          IF (CoilBody .AND. CoilType /= 'massive' .AND. .NOT. NoRotM ) THEN
            RotMLoc(i,j) = SUM( RotM(i,j,1:n) * Basis(1:n) )
          END IF
        END DO
      END DO
      
      ! Transform the conductivity tensor (in case of a foil winding):
      ! --------------------------------------------------------------
      IF (CoilBody .AND. CoilType /= 'massive' .AND. .NOT. NoRotM) THEN
        C = MATMUL(MATMUL(RotMLoc, C),TRANSPOSE(RotMLoc))
      END IF

      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      DO i=1,nd
        DO j=1,nd
          ! Finally, the elemental matrix & vector:
          !----------------------------------------
          STIFF(i,j) = STIFF(i,j) + SUM(MATMUL(C, dBasisdx(i,:)) * dBasisdx(j,:))*detJ*IP % s(t)
        END DO
      END DO

      FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * LoadAtIP * Basis(1:nd)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LCondensate( N, Nb, K, F )
!------------------------------------------------------------------------------
    INTEGER :: N, Nb
    REAL(KIND=dp) :: K(:,:),F(:),Kbb(Nb,Nb), &
         Kbl(Nb,N), Klb(N,Nb), Fb(Nb)

    INTEGER :: m, i, j, l, p, Ldofs(N), Bdofs(Nb)

    IF ( Nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:n,1:n) = &
         K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
  END SUBROUTINE LCondensate
!------------------------------------------------------------------------------

!----------------------------------------------------------------
  SUBROUTINE BoundaryCondition(LOAD, FORCE, Element, n)
!----------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), DIMENSION(:) :: FORCE, LOAD
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!----------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
    REAL(KIND=dp) :: detJ, LoadAtIP,&
    LocalHeatCapacity, LocalDensity
    LOGICAL :: stat, getSecondDerivatives
    INTEGER :: t,j
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!----------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    FORCE = 0.0d0
    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    !-----------------------------------------------------------------
    ! Loop over Gauss-points (boundary element Integration)
    !-----------------------------------------------------------------
    DO t=1,IP % n
      !Basis function values & derivatives at the integration point:
      !-------------------------------------------------------------
      getSecondDerivatives = .FALSE.
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
      IP % W(t), detJ, Basis, dBasisdx)

      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      DO j=1,n
        FORCE(j) = FORCE(j) + IP % s(t)*DetJ*LoadAtIP*Basis(j)
      END DO
    END DO
  END SUBROUTINE BoundaryCondition

!!------------------------------------------------------------------------------
!  SUBROUTINE SaveWPotSolution(nofbodies, Tcoef)
!!------------------------------------------------------------------------------
!  IMPLICIT NONE
!  INTEGER :: Active, n, t, nn, ns_iter, nofbodies
!  REAL(KIND=dp) :: Tcoef(:,:,:)
!  REAL(KIND=dp) :: Wnorms(nofbodies)
!  TYPE(Element_t), POINTER :: Element
!!------------------------------------------------------------------------------
!!------------------------------------------------------------------------------
!  END SUBROUTINE SaveWPotSolution
!!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE SaveElementWSolution(Element, n, Wnorm, RotM, Tcoef, NoRotM)
!------------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: n, j, k, t
  TYPE(Element_t), POINTER :: Element
  TYPE(Valuelist_t), POINTER :: Solverparams
  TYPE(Variable_t), POINTER, SAVE :: wpotvar
  TYPE(Variable_t), POINTER, SAVE :: wvecvar
  TYPE(Variable_t), POINTER, SAVE :: jvecvar
  LOGICAL, SAVE :: First=.TRUE.
  LOGICAL :: stat, ComputeJvec, NoRotM
  REAL(KIND=dp) :: Basis(n), DetJ
  REAL(KIND=dp) :: dBasisdx(n,3)
  CHARACTER(LEN=MAX_NAME_LEN), SAVE :: varname, jvarname
  REAL(KIND=dp) :: wpot(n), wvec(3), jvec(3)
  REAL(KIND=dp) :: Wnorm
  TYPE(GaussIntegrationPoints_t) :: IP
  TYPE(Nodes_t) :: Nodes
  REAL(KIND=dp) :: Tcoef(:,:,:)
  REAL(KIND=dp) :: C(3,3), RotMLoc(3,3), RotM(3,3,n)
  SAVE Nodes

  CALL GetElementNodes(Nodes)
  varname = GetString(GetSolverParams(), 'Variable', Found)
  IF (.NOT. Found) varname = 'W'
  ComputeJVec = GetLogical(GetSolverParams(), 'Compute J Vector', Found)
  IF (.NOT. Found) ComputeJVec = .False.

  IF (ComputeJVec) THEN
    jvarname = GetString(GetSolverParams(), 'J Variable', Found)
    IF (.NOT. Found) jvarname = 'J'
  END IF

  IF (First) THEN
    First = .FALSE.
    wpotvar => VariableGet( Mesh % Variables, TRIM(varname)//' Potential')
    IF(.NOT. ASSOCIATED(wpotvar)) THEN
      CALL Fatal('SaveElementWSolution()','W Potential variable not found')
    END IF
    wvecvar => VariableGet( Mesh % Variables, TRIM(varname)//' Vector E')
    IF(.NOT. ASSOCIATED(wvecvar)) THEN
      CALL Fatal('SaveElementWSolution()','W Vector variable not found')
    END IF

    IF (ComputeJVec) THEN
      jvecvar => VariableGet( Mesh % Variables, TRIM(jvarname)//' Vector E')
      IF(.NOT. ASSOCIATED(jvecvar)) THEN
        CALL Fatal('SaveElementWSolution()','J Vector variable not found')
      END IF
    END IF
  END IF
   
  CALL GetLocalSolution(wpot, varname)
  IP = GaussPoints(Element)
  DO t=1,n
    IF (ASSOCIATED(wpotvar)) THEN
      DO k=1,wpotvar % DOFs
        IF( CoilType/='stranded' ) Wnorm = 1._dp
        IF (Wnorm > EPSILON(Wnorm)) THEN
!          print *, ParEnv % MyPe, "Wnorm:", Wnorm
          wpotvar % Values( wpotvar % DOFs*(wpotvar % Perm( &
              Element % DGIndexes(t))-1)+k) = wpot(t)/Wnorm
        END IF
      END DO
    END IF

    IF (ASSOCIATED(wvecvar)) THEN
      ! ----------------------
      ! Basis function values & derivatives at the integration point
      ! That hopefully are at the vardof locations:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis,dBasisdx )

      wvec = MATMUL(wpot(1:n), dBasisdx(1:n,:))

      DO k=1,wvecvar % DOFs
        IF (Wnorm > EPSILON(Wnorm)) THEN
          IF( CoilType/='stranded' ) Wnorm = 1._dp
!          print *, ParEnv % MyPe, "Wnorm:", Wnorm
          wvecvar % Values( wvecvar % DOFs*(wvecvar % Perm( &
              Element % DGIndexes(t))-1)+k) = wvec(k)/Wnorm
        END IF
      END DO

      IF (ComputeJVec) THEN
        IF (ASSOCIATED(jvecvar)) THEN

          DO i=1,3
            DO j=1,3
              C(i,j) = SUM( Tcoef(i,j,1:n) * Basis(1:n) )
              IF (CoilBody .AND. CoilType /= 'massive' .AND. .NOT. NoRotM ) THEN
                RotMLoc(i,j) = SUM( RotM(i,j,1:n) * Basis(1:n) )
              END IF
            END DO
          END DO
          
          ! Transform the conductivity tensor (in case of a foil winding):
          ! --------------------------------------------------------------
          IF (CoilBody .AND. CoilType /= 'massive' .AND. .NOT. NoRotM ) THEN
            C = MATMUL(MATMUL(RotMLoc, C),TRANSPOSE(RotMLoc))
          END IF
          jvec = MATMUL(C, wvec)

          DO k=1,wvecvar % DOFs
            IF (Wnorm > EPSILON(Wnorm)) THEN
    !          print *, ParEnv % MyPe, "Wnorm:", Wnorm
              jvecvar % Values( jvecvar % DOFs*(jvecvar % Perm( &
                  Element % DGIndexes(t))-1)+k) = jvec(k)
            END IF
          END DO
        END IF
      END IF
    END IF
  END DO 
!------------------------------------------------------------------------------
  END SUBROUTINE SaveElementWSolution
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION GetWNormsForBodies(nofbodies) RESULT (Wnorms)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nofbodies
    REAL(KIND=dp) :: Volumes(nofbodies)
    REAL(KIND=dp) :: Wnorms(nofbodies)
    REAL(KIND=dp) :: WnormCoeffs(nofbodies)
    INTEGER :: Active, t, n, i
    TYPE(Element_t),POINTER :: Element
    LOGICAL :: CoilBody, Found, stat
    CHARACTER(LEN=MAX_NAME_LEN):: CoilType
    TYPE(ValueList_t), POINTER :: CompParams

    Wnorms = 0._dp
    Volumes = 0._dp
    WnormCoeffs = 0._dp
    Active = GetNOFActive()
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()      
      nd = GetElementNOFDOFs()
      
      CoilBody = .FALSE.
      CompParams => GetComponentParams( Element )
      CoilType = ''
      IF (ASSOCIATED(CompParams)) THEN
        CoilType = GetString(CompParams, 'Coil Type', Found)
        IF (Found) CoilBody = .TRUE.
      END IF 

      IF (CoilBody) THEN
        CALL AddElementWNormAndVolume(Element, n, nd, &
            WnormCoeffs(Element % BodyId), &
            Volumes(Element % BodyId))
      END IF 
    END DO
   
    DO i=1, nofbodies
       WnormCoeffs(i) = ParallelReduction(WnormCoeffs(i)) 
       Volumes(i) = ParallelReduction(Volumes(i))
       IF (Volumes(i) /= 0.0_dp) THEN 
         Wnorms(i) = WnormCoeffs(i) &
               /   Volumes(i)
       ELSE
         Wnorms(i) = 0.0_dp
       END IF
    END DO

    CALL Info('GetWnormsForBodies', &
         'W norm computed in components.', level=9)

!------------------------------------------------------------------------------
  END FUNCTION GetWNormsForBodies
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE AddElementWNormAndVolume(Element, n, nd, WnormCoeff, Volume)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: WnormCoeff, Volume
    
    INTEGER :: Active, n, nd, j
    TYPE(Element_t),POINTER :: Element
    LOGICAL :: CoilBody, Found, stat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ
    REAL(KIND=dp) :: Wbase(n), s, w(3)

    CALL GetElementNodes( Nodes )
    CALL GetLocalSolution(Wbase,'W')
  
    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )

    DO j=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(j), IP % V(j), &
          IP % W(j), detJ, Basis, dBasisdx )
      
      ! Compute the Element Volume
      ! -----------------------------
      s = IP % s(j) * detJ
      Volume = Volume + s

      ! Compute the norm of W
      ! ---------------------
      w = -MATMUL(WBase(1:n), dBasisdx(1:n,:))
      WnormCoeff = WnormCoeff + SQRT(SUM(w**2))*s
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE AddElementWNormAndVolume
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
END SUBROUTINE Wsolve
!------------------------------------------------------------------------------
