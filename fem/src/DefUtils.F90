!/*****************************************************************************/
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
! *  Original Date: 12 Jun 2003
! *
! ******************************************************************************/

!> \defgroup ElmerLib Elmer library  
!> \{
!> \defgroup DefUtils Default API 
!> \{

!--------------------------------------------------------------------------------
!>  Module containing utility subroutines with default values for various
!>  system subroutine arguments. For user defined solvers and subroutines
!> this module should provide all the needed functionality for typical finite
!> element procedures.  
!--------------------------------------------------------------------------------
MODULE DefUtils

#include "../config.h"

   USE Adaptive
   USE MeshGenerate 
   USE ElementUtils
   USE SolverUtils

   IMPLICIT NONE

   INTERFACE DefaultUpdateEquations
     MODULE PROCEDURE DefaultUpdateEquationsR, DefaultUpdateEquationsC, &
         DefaultUpdateEquationsDiagC
   END INTERFACE

   INTERFACE DefaultUpdatePrec
     MODULE PROCEDURE DefaultUpdatePrecR, DefaultUpdatePrecC
   END INTERFACE

   INTERFACE DefaultUpdateMass
     MODULE PROCEDURE DefaultUpdateMassR, DefaultUpdateMassC
   END INTERFACE

   INTERFACE DefaultUpdateBulk
     MODULE PROCEDURE DefaultUpdateBulkR, DefaultUpdateBulkC
   END INTERFACE

   INTERFACE DefaultUpdateDamp
     MODULE PROCEDURE DefaultUpdateDampR, DefaultUpdateDampC
   END INTERFACE

   INTERFACE DefaultUpdateForce
     MODULE PROCEDURE DefaultUpdateForceR, DefaultUpdateForceC
   END INTERFACE

   INTERFACE DefaultUpdateTimeForce
     MODULE PROCEDURE DefaultUpdateTimeForceR, DefaultUpdateTimeForceC
   END INTERFACE

   INTERFACE Default1stOrderTime
     MODULE PROCEDURE Default1stOrderTimeR, Default1stOrderTimeC
   END INTERFACE

   INTERFACE Default2ndOrderTime
     MODULE PROCEDURE Default2ndOrderTimeR, Default2ndOrderTimeC
   END INTERFACE

   INTERFACE GetLocalSolution
     MODULE PROCEDURE GetScalarLocalSolution, GetVectorLocalSolution
   END INTERFACE

   INTERFACE GetLocalEigenmode
     MODULE PROCEDURE GetScalarLocalEigenmode, GetVectorLocalEigenmode
   END INTERFACE

   INTEGER, ALLOCATABLE, TARGET, PRIVATE :: IndexStore(:), VecIndexStore(:)
   REAL(KIND=dp), ALLOCATABLE, TARGET, PRIVATE  :: ValueStore(:)
   !$OMP THREADPRIVATE(IndexStore, VecIndexStore, ValueStore)

   TYPE(Element_t), POINTER :: CurrentElementThread => NULL()
   !$OMP THREADPRIVATE(CurrentElementThread)
   
   ! TODO: Get actual values for these from mesh
   INTEGER, PARAMETER, PRIVATE :: ISTORE_MAX_SIZE = 1024
   INTEGER, PARAMETER, PRIVATE :: VSTORE_MAX_SIZE = 1024
   PRIVATE :: GetIndexStore, GetPermIndexStore, GetValueStore
CONTAINS


   FUNCTION GetVersion() RESULT(ch)
     CHARACTER(LEN=:), ALLOCATABLE :: ch
     ch = VERSION
   END FUNCTION GetVersion

   FUNCTION GetSifName(Found) RESULT(ch)
     CHARACTER(LEN=:), ALLOCATABLE :: ch
     LOGICAL, OPTIONAL :: Found     
     ch = GetString(CurrentModel % Simulation,'Solver Input File',Found)
   END FUNCTION GetSifName
    
   FUNCTION GetRevision(Found) RESULT(ch)
     CHARACTER(LEN=:), ALLOCATABLE :: ch
     LOGICAL, OPTIONAL :: Found
#ifdef REVISION
     ch = REVISION
     IF(PRESENT(Found)) Found = .TRUE.
#else
     ch = "unknown"
     IF(PRESENT(Found)) Found = .FALSE.
#endif
   END FUNCTION GetRevision

   FUNCTION GetCompilationDate(Found) RESULT(ch)
     CHARACTER(LEN=:), ALLOCATABLE :: ch
     LOGICAL, OPTIONAL :: Found
#ifdef COMPILATIONDATE
     ch = COMPILATIONDATE
     IF(PRESENT(Found)) Found = .TRUE.
#else
     ch = "unknown"
     IF(PRESENT(Found)) Found = .FALSE.
#endif
   END FUNCTION GetCompilationDate
  
  FUNCTION GetIndexStore() RESULT(ind)
    IMPLICIT NONE
    INTEGER, POINTER CONTIG :: ind(:)
    INTEGER :: istat

    IF ( .NOT. ALLOCATED(IndexStore) ) THEN
        ALLOCATE( IndexStore(ISTORE_MAX_SIZE), STAT=istat )
        IndexStore = 0
        IF ( Istat /= 0 ) CALL Fatal( 'GetIndexStore', &
                'Memory allocation error.' )
    END IF
    ind => IndexStore
  END FUNCTION GetIndexStore

  FUNCTION GetPermIndexStore() RESULT(ind)
    IMPLICIT NONE
    INTEGER, POINTER CONTIG :: ind(:)
    INTEGER :: istat
     
    IF ( .NOT. ALLOCATED(VecIndexStore) ) THEN
      ALLOCATE( VecIndexStore(ISTORE_MAX_SIZE), STAT=istat )
      VecIndexStore = 0
      IF ( istat /= 0 ) CALL Fatal( 'GetPermIndexStore', &
              'Memory allocation error.' )
    END IF
    ind => VecIndexStore
  END FUNCTION GetPermIndexStore

  FUNCTION GetValueStore(n) RESULT(val)
    IMPLICIT NONE
    REAL(KIND=dp), POINTER CONTIG :: val(:)
    INTEGER :: n, istat

    IF ( .NOT.ALLOCATED(ValueStore) ) THEN
      ALLOCATE( ValueStore(VSTORE_MAX_SIZE), STAT=istat )
      ValueStore = REAL(0, dp)
      IF ( Istat /= 0 ) CALL Fatal( 'GetValueStore', &
              'Memory allocation error.' )
    END IF
    IF (n > VSTORE_MAX_SIZE) THEN
      CALL Fatal( 'GetValueStore', 'Not enough memory allocated for store.' )
    END IF
    val => ValueStore(1:n)
  END FUNCTION GetValueStore

!> Returns handle to the active solver
  FUNCTION GetSolver() RESULT( Solver )
     TYPE(Solver_t), POINTER :: Solver
     Solver => CurrentModel % Solver
  END FUNCTION GetSolver

!> Returns handle to the active matrix 
  FUNCTION GetMatrix( USolver ) RESULT( Matrix )
     TYPE(Matrix_t), POINTER :: Matrix 
     TYPE(Solver_t), OPTIONAL, TARGET :: USolver

     IF ( PRESENT( USolver ) ) THEN
        Matrix => USolver % Matrix
     ELSE
        Matrix => CurrentModel % Solver % Matrix
     END IF
  END FUNCTION GetMatrix

!> Returns handle to the active mesh
  FUNCTION GetMesh( USolver ) RESULT( Mesh )
     TYPE(Mesh_t), POINTER :: Mesh 
     TYPE(Solver_t), OPTIONAL, TARGET :: USolver

     IF ( PRESENT( USolver ) ) THEN
        Mesh => USolver % MEsh
     ELSE
        Mesh => CurrentModel % Solver % Mesh
     END IF
  END FUNCTION GetMesh

!> Returns handle to the active element
  FUNCTION GetCurrentElement(Element) RESULT(Ret_Element)
    IMPLICIT NONE
    TYPE(Element_t), OPTIONAL, TARGET :: Element
    TYPE(Element_t), POINTER :: Ret_Element

    IF (PRESENT(Element)) THEN
      Ret_Element=>Element
    ELSE
#ifdef _OPENMP
      IF (omp_in_parallel()) THEN
        Ret_Element=>CurrentElementThread
      ELSE
        Ret_Element=>CurrentModel % CurrentElement
      END IF
#else
      Ret_Element => CurrentModel % CurrentElement
#endif
    END IF
  END FUNCTION GetCurrentElement

!> Sets handle to the active element of the current thread. 
!> Old handle is given as a return value as what would be returned
!> by a call to GetCurrentElement
  FUNCTION SetCurrentElement(Element) RESULT(OldElement)
    IMPLICIT NONE
    TYPE(Element_t), TARGET :: Element
    TYPE(Element_t), POINTER :: OldElement

#ifdef _OPENMP
    IF (omp_in_parallel()) THEN
      OldElement => CurrentElementThread
      CurrentElementThread => Element
    ELSE
      OldElement => CurrentModel % CurrentElement
      CurrentModel % CurrentElement => Element
    END IF
#else
    OldElement => CurrentModel % CurrentElement
    CurrentModel % CurrentElement => Element
#endif
  END FUNCTION SetCurrentElement

!> Returns handle to the index of the current element
  FUNCTION GetElementIndex(Element) RESULT(Indx)
    TYPE(Element_t), OPTIONAL :: Element
    INTEGER :: Indx
    TYPE(Element_t), POINTER :: CurrElement
    
    CurrElement => GetCurrentElement(Element)
    Indx = CurrElement % ElementIndex
  END FUNCTION GetElementIndex

  SUBROUTINE GetElementNodeIndex(i, Element, n, FOUND)
    IMPLICIT None
 
    ! variables in function header
    INTEGER :: i, n
    TYPE(Element_t), POINTER :: Element
    Logical :: FOUND
    
    DO i=1, SIZE(Element%NodeIndexes)
      IF (n == Element%NodeIndexes(i)) THEN
        FOUND=.TRUE.
        EXIT
      END IF
    END DO
  END SUBROUTINE GetElementNodeIndex

  FUNCTION GetIPIndex( LocalIp, USolver, Element, IpVar ) RESULT ( GlobalIp ) 
    INTEGER :: LocalIp, GlobalIp

    TYPE(Solver_t), OPTIONAL, TARGET :: USolver
    TYPE(Element_t), OPTIONAL :: Element
    TYPE(Variable_t), POINTER, OPTIONAL :: IpVar

    TYPE(Solver_t), POINTER :: Solver
    TYPE(Element_t), POINTER :: CurrElement
    INTEGER :: n, m
    INTEGER, POINTER :: IpPerm(:)
    
    IF ( PRESENT( USolver ) ) THEN
      Solver => USolver
    ELSE
      Solver => CurrentModel % Solver 
    END IF
    
    CurrElement => GetCurrentElement(Element)
    n = CurrElement % ElementIndex
    GlobalIp = 0
    
    IF( PRESENT( IpVar ) ) THEN
      IF( IpVar % TYPE /= Variable_on_gauss_points ) THEN
        CALL Fatal('GetIpIndex','Variable is not of type gauss points!')
      END IF
      
      IpPerm => IpVar % Perm
      m = IpPerm(n+1) - IpPerm(n)

      ! This is a sign that the variable is not active at the element
      IF( m == 0 ) RETURN
    ELSE
      IF( .NOT. ASSOCIATED( Solver % IpTable ) ) THEN
        CALL Fatal('GetIpIndex','Cannot access index of gaussian point!')
      END IF

      IpPerm => Solver % IpTable % IpOffset 
      m = IpPerm(n+1) - IpPerm(n)
    END IF

    ! There are not sufficient number of gauss points in the permutation table to have a
    ! local index this big. 
    IF( m < LocalIp ) THEN      
      CALL Warn('GetIpIndex','Inconsistent number of IP points!')
      RETURN
    END IF
    
    GlobalIp = IpPerm(n) + LocalIp

  END FUNCTION GetIPIndex

  
  
  FUNCTION GetIPCount( USolver, IpVar ) RESULT ( IpCount ) 
    INTEGER :: IpCount
    TYPE(Solver_t), OPTIONAL, TARGET :: USolver
    TYPE(Variable_t), OPTIONAL, POINTER :: IpVar
    
    TYPE(Solver_t), POINTER :: Solver
    INTEGER, POINTER :: IpPerm 

    IF ( PRESENT( USolver ) ) THEN
      Solver => USolver
    ELSE
      Solver => CurrentModel % Solver 
    END IF

    IF( PRESENT( IpVar ) ) THEN
      IF( IpVar % TYPE /= Variable_on_gauss_points ) THEN
        CALL Fatal('GetIpCount','Variable is not of type gauss points!')
      END IF
      IpCount = SIZE( IpVar % Values ) / IpVar % Dofs
    ELSE    
      IF( .NOT. ASSOCIATED( Solver % IpTable ) ) THEN
        CALL Fatal('GetIpCount','Gauss point table not initialized')
      END IF    
      IpCount = Solver % IpTable % IpCount 
    END IF
      
  END FUNCTION GetIPCount

  
!> Returns the number of active elements for the current solver
  FUNCTION GetNOFActive( USolver ) RESULT(n)
     INTEGER :: n
     TYPE(Solver_t), OPTIONAL, TARGET :: USolver
     TYPE(Solver_t), POINTER :: Solver

     IF ( PRESENT( USolver ) ) THEN
       Solver => USolver
     ELSE
       Solver => CurrentModel % Solver 
     END IF

     IF( ASSOCIATED( Solver % ColourIndexList ) ) THEN
       Solver % CurrentColour = Solver % CurrentColour + 1
       n = Solver % ColourIndexList % ptr(Solver % CurrentColour+1) &
           - Solver % ColourIndexList % ptr(Solver % CurrentColour)
       CALL Info('GetNOFActive','Number of active elements: '&
           //I2S(n)//' in colour '//I2S(Solver % CurrentColour),Level=20)
     ELSE
       n = Solver % NumberOfActiveElements
       CALL Info('GetNOFActive','Number of active elements: '&
           //I2S(n),Level=20)
     END IF

  END FUNCTION GetNOFActive

  !> Return number of boundary elements of the current boundary colour
  !> and increments the colour counter
  FUNCTION GetNOFBoundaryActive( USolver ) RESULT(n)
     INTEGER :: n
     TYPE(Solver_t), OPTIONAL, TARGET :: USolver
     TYPE(Solver_t), POINTER :: Solver

     IF ( PRESENT( USolver ) ) THEN
       Solver => USolver
     ELSE
       Solver => CurrentModel % Solver 
     END IF

     IF( ASSOCIATED( Solver % BoundaryColourIndexList ) ) THEN
       Solver % CurrentBoundaryColour = Solver % CurrentBoundaryColour + 1
       n = Solver % BoundaryColourIndexList % ptr(Solver % CurrentBoundaryColour+1) &
           - Solver % BoundaryColourIndexList % ptr(Solver % CurrentBoundaryColour)
       CALL Info('GetNOFBoundaryActive','Number of boundary elements: '&
           //I2S(n)//' in colour '//I2S(Solver % CurrentBoundaryColour),Level=20)
     ELSE
       n = Solver % Mesh % NumberOfBoundaryElements
       CALL Info('GetNOFBoundaryActive','Number of active elements: '&
           //I2S(n),Level=20)
     END IF

  END FUNCTION GetNOFBoundaryActive

!> Returns the current time
  FUNCTION GetTime() RESULT(st)
     REAL(KIND=dp) :: st
     TYPE(Variable_t), POINTER :: v

     v => CurrentModel % Solver % Mesh % Variables
     v => VariableGet( v, 'time' )
     st = v % Values(1)
  END FUNCTION GetTime

!> Returns the current periodic time
  FUNCTION GetPeriodicTime() RESULT(st)
     REAL(KIND=dp) :: st
     TYPE(Variable_t), POINTER :: v

     v => CurrentModel % Solver % Mesh % Variables
     v => VariableGet( v, 'periodic time' )
     st = v % Values(1)
  END FUNCTION GetPeriodicTime

!> Returns the current timestep  
  FUNCTION GetTimeStep() RESULT(st)
     INTEGER :: st
     TYPE(Variable_t), POINTER :: v

     v => CurrentModel % Solver % Mesh % Variables
     v => VariableGet( v, 'timestep' )
     st = NINT(v % Values(1))
  END FUNCTION GetTimestep

!> Returns the current timestep interval
  FUNCTION GetTimeStepInterval() RESULT(st)
     INTEGER :: st
     TYPE(Variable_t), POINTER :: v

     v => CurrentModel % Solver % Mesh % Variables
     v => VariableGet( v, 'timestep interval')
     st = NINT(v % Values(1))
  END FUNCTION GetTimestepInterval

!> Returns the current timestep size  
  FUNCTION GetTimestepSize() RESULT(st)
     REAL(KIND=dp) :: st
     TYPE(Variable_t), POINTER :: v

     v => CurrentModel % Solver % Mesh % Variables
     v => VariableGet( v, 'timestep size')
     st = v % Values(1)
  END FUNCTION GetTimestepSize
 
!> Returns the angular frequency  
  FUNCTION GetAngularFrequency(ValueList,Found, UElement) RESULT(w)
    REAL(KIND=dp) :: w
    TYPE(ValueList_t), POINTER, OPTIONAL :: ValueList
    LOGICAL, OPTIONAL :: Found
    TYPE(Element_t), OPTIONAL :: UElement

    w = ListGetAngularFrequency( ValueList, Found, UElement )
  END FUNCTION GetAngularFrequency

!> Returns the current coupled system iteration loop count
  FUNCTION GetCoupledIter() RESULT(st)
     INTEGER :: st
     TYPE(Variable_t), POINTER :: v

     v => CurrentModel % Solver % Mesh % Variables
     v => VariableGet( v, 'coupled iter')
     st = NINT(v % Values(1))
  END FUNCTION GetCoupledIter

!> Returns the current nonlinear system iteration loop count
  FUNCTION GetNonlinIter() RESULT(st)
     INTEGER :: st
     TYPE(Variable_t), POINTER :: v

     v => CurrentModel % Solver % Mesh % Variables
     v => VariableGet( v, 'nonlin iter')
     st = NINT(v % Values(1))
  END FUNCTION GetNonlinIter

!> Returns the number of boundary elements  
  FUNCTION GetNOFBoundaryElements( UMesh ) RESULT(n)
    INTEGER :: n
    TYPE(Mesh_t), OPTIONAL :: UMesh
  
    IF ( PRESENT( UMesh ) ) THEN
       n = UMesh % NumberOfBoundaryElements
    ELSE
       n = CurrentModel % Mesh % NumberOfBoundaryElements
    END IF
  END FUNCTION GetNOFBoundaryElements

!> Returns a scalar field in the nodes of the element
  SUBROUTINE GetScalarLocalSolution( x,name,UElement,USolver,tStep, UVariable, Found)
     REAL(KIND=dp) :: x(:)
     CHARACTER(LEN=*), OPTIONAL :: name
     TYPE(Solver_t)  , OPTIONAL, TARGET :: USolver
     TYPE(Element_t),  OPTIONAL, TARGET :: UElement
     TYPE(Variable_t), OPTIONAL, TARGET :: UVariable
     INTEGER, OPTIONAL :: tStep
     LOGICAL, OPTIONAL :: Found

     REAL(KIND=dp), POINTER :: Values(:)
     TYPE(Variable_t), POINTER :: Variable
     TYPE(Solver_t)  , POINTER :: Solver
     TYPE(Element_t),  POINTER :: Element, Parent

     INTEGER :: i, j, k, n, lr
     INTEGER, POINTER :: Indexes(:)
     LOGICAL :: Found0
     
     Solver => CurrentModel % Solver
     IF ( PRESENT(USolver) ) Solver => USolver

     x = 0.0d0
     IF(PRESENT(Found)) Found = .FALSE.

     IF(.NOT. PRESENT(UVariable)) THEN
       Variable => Solver % Variable
     ELSE
       Variable => UVariable
     END IF
     
     IF ( PRESENT(name) ) THEN
        Variable => VariableGet( Solver % Mesh % Variables, name )
     END IF
     IF ( .NOT. ASSOCIATED( Variable ) ) RETURN

     Values => Variable % Values
     IF ( PRESENT(tStep) ) THEN
       IF ( tStep<0 ) THEN
         IF ( ASSOCIATED(Variable % PrevValues) ) THEN
           IF( -tStep<=SIZE(Variable % PrevValues,2)) &
              Values => Variable % PrevValues(:,-tStep)
         END IF
       END IF
     END IF

     Element => GetCurrentElement(UElement)
     Found0 = .FALSE.
     
     ! Some variables do not really follow the numbering
     ! nodes + faces + edges etc. of the standard solver.
     ! For example, if we want to request DG values from a variable
     ! that is not called by a DG solver we have to treat the DG variable
     ! separately. As is the case for Gauss variables.
     ! If variable is defined on gauss points return that instead
     !-------------------------------------------------------------
     IF( Variable % TYPE == Variable_on_gauss_points ) THEN
       j = Element % ElementIndex
       n = Variable % Perm(j+1) - Variable % Perm(j)
       DO i=1,n
         x(i) = Values(Variable % Perm(j) + i)
       END DO
       IF(PRESENT(Found)) Found = (n>1)
       RETURN
     ELSE IF( Variable % TYPE == Variable_on_nodes_on_elements ) THEN
       n = Element % TYPE % NumberOfNodes
       Indexes => Element % DGIndexes       
       IF(ASSOCIATED( Indexes ) ) THEN
         DO i=1,n
           j = Variable % Perm(Indexes(i))
           IF(j>0) THEN
             Found0 = .TRUE.
             x(i) = Values(j)
           END IF
         END DO
       ELSE IF ( ASSOCIATED( Element % BoundaryInfo ) ) THEN
         DO lr=1,2
           IF(lr==1) THEN
             Parent => Element % BoundaryInfo % Left
           ELSE
             Parent => Element % BoundaryInfo % Right
           END IF
           IF(.NOT. ASSOCIATED( Parent ) ) CYCLE
           IF( ANY( Variable % Perm( Parent % DGIndexes ) == 0) ) CYCLE                                     
           DO i=1,n
             DO j=1,Parent % TYPE % NumberOfNodes
               IF( Element % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
                 k = Variable % Perm( Parent % DGIndexes(j) )
                 IF(k>0) THEN
                   Found0 = .TRUE.
                   x(i) = Values(k)
                 END IF
                 EXIT
               END IF
             END DO
           END DO
           EXIT
         END DO
       END IF
       IF(PRESENT(Found)) Found = Found0
       RETURN       
     END IF

     Indexes => GetIndexStore()
     IF ( ASSOCIATED(Variable % Solver) ) THEN
       n = GetElementDOFs( Indexes, Element, Variable % Solver )
     ELSE
       n = GetElementDOFs( Indexes, Element, Solver )
     END IF
     n = MIN( n, SIZE(x) )


     IF ( ASSOCIATED( Variable % Perm ) ) THEN
       IF( Variable % PeriodicFlipActive ) THEN
         DO i=1,n
           j = Indexes(i)
           IF ( j>0 .AND. j<=SIZE(Variable % Perm) ) THEN
             k = Variable % Perm(j)
             IF ( k>0 ) THEN
               Found0 = .TRUE.
               x(i) = Values(k)
               IF( CurrentModel % Mesh % PeriodicFlip(j) ) x(i) = -x(i)
             END IF
           END IF
         END DO
       ELSE
         DO i=1,n
           j = Indexes(i)
           IF ( j>0 .AND. j<=SIZE(Variable % Perm) ) THEN
             j = Variable % Perm(j)
             IF ( j>0 ) THEN
               Found0 = .TRUE.
               x(i) = Values(j)
             END IF
           END IF
         END DO
       END IF
     ELSE
        DO i=1,n
          j = Indexes(i)
          IF ( j>0 .AND. j<=SIZE(Variable % Values) ) THEN
            Found0 = .TRUE.
            x(i) = Values(Indexes(i))
          END IF
        END DO
      END IF

      IF(PRESENT(Found)) Found = Found0 
      
  END SUBROUTINE GetScalarLocalSolution



!> Returns a vector field in the nodes of the element
  SUBROUTINE GetVectorLocalSolution( x,name,UElement,USolver,tStep, UVariable, Found)
     REAL(KIND=dp) :: x(:,:)
     CHARACTER(LEN=*), OPTIONAL :: name
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     TYPE(Variable_t), OPTIONAL, TARGET :: UVariable
     INTEGER, OPTIONAL :: tStep
     LOGICAL, OPTIONAL :: Found
     
     TYPE(Variable_t), POINTER :: Variable
     TYPE(Solver_t)  , POINTER :: Solver
     TYPE(Element_t),  POINTER :: Element, Parent

     INTEGER :: i, j, k, l, lr, n
     INTEGER, POINTER :: Indexes(:)
     REAL(KIND=dp), POINTER ::  Values(:)
     LOGICAL :: Found0
     
     Solver => CurrentModel % Solver
     IF ( PRESENT(USolver) ) Solver => USolver

     x = 0.0d0
     IF(PRESENT(Found)) Found = .FALSE.
     
     IF(.NOT. PRESENT(UVariable)) THEN
       Variable => Solver % Variable
     ELSE
       Variable => UVariable
     END IF

     IF ( PRESENT(name) ) THEN
        Variable => VariableGet( Solver % Mesh % Variables, name )
     END IF
     IF ( .NOT. ASSOCIATED( Variable ) ) RETURN

     Values => Variable % Values
     IF ( PRESENT(tStep) ) THEN
       IF ( tStep<0 ) THEN
         IF ( ASSOCIATED(Variable % PrevValues) ) THEN
           IF ( -tStep<=SIZE(Variable % PrevValues,2)) &
             Values => Variable % PrevValues(:,-tStep)
         END IF
       END IF
     END IF
     
     Element => GetCurrentElement(UElement)
     Found0 = .FALSE.

       
     ! If variable is defined on gauss points return that instead
     IF( Variable % TYPE == Variable_on_gauss_points ) THEN
       ASSOCIATE(dofs => variable % dofs)
         j = Element % ElementIndex
         n = Variable % Perm(j+1) - Variable % Perm(j)
         IF (SIZE(x,1) < dofs .OR. SIZE(x,2) < n) THEN
           WRITE (message,*) 'Attempting to get IP solution to a too small array of size', &
               SHAPE(x), '. Required size:', dofs, n
           CALL Fatal('GetVectorLocalSolution', message)
         END IF
         DO i=1,n
           ASSOCIATE(p => variable % perm(j) + i)
             DO k=1,dofs
               x(k, i) = Values((p-1)*dofs + k)
             END DO
           END ASSOCIATE
         END DO
         IF(PRESENT(Found)) Found = (n>1)
         RETURN
       END ASSOCIATE
     ELSE IF(  Variable % TYPE == Variable_on_nodes_on_elements ) THEN
       n = Element % TYPE % NumberOfNodes
       Indexes => Element % DGIndexes
       IF(ASSOCIATED( Indexes ) ) THEN
         ASSOCIATE(dofs => variable % dofs)
           DO i=1,n
             j = variable % perm(indexes(i))
             IF( j==0 ) CYCLE
             Found0 = .TRUE.
             DO k=1,dofs
               x(k,i) = Values((j-1)*dofs + k)
             END DO
           END DO
         END ASSOCIATE
         IF(PRESENT(Found)) Found = Found0
         RETURN         
       ELSE IF ( ASSOCIATED( Element % BoundaryInfo ) ) THEN
         DO lr=1,2
           IF(lr==1) THEN
             Parent => Element % BoundaryInfo % Left
           ELSE
             Parent => Element % BoundaryInfo % Right
           END IF
           IF(.NOT. ASSOCIATED( Parent ) ) CYCLE
           IF( ANY( Variable % Perm( Parent % DGIndexes ) == 0) ) CYCLE                          

           ASSOCIATE(dofs => variable % dofs)
             DO i=1,n
               DO j=1,Parent % TYPE % NumberOfNodes
                 IF( Element % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
                   l = Variable % Perm( Parent % DGIndexes(j) )
                   IF(l>0) THEN
                     Found0 = .TRUE.
                     DO k=1,dofs
                       x(k,i) = Values((l-1)*dofs + k )
                     END DO
                   END IF
                   EXIT
                 END IF
               END DO
             END DO
           END ASSOCIATE
           EXIT
         END DO
       END IF
       IF(PRESENT(Found)) Found = Found0
       RETURN       
     END IF

     
     Indexes => GetIndexStore()
     IF ( ASSOCIATED(Variable % Solver ) ) THEN
       n = GetElementDOFs( Indexes, Element, Variable % Solver )
     ELSE
       n = GetElementDOFs( Indexes, Element, Solver )
     END IF
     n = MIN( n, SIZE(x,2) )
     
     DO i=1,Variable % DOFs
       IF ( ASSOCIATED( Variable % Perm ) ) THEN
         IF( Variable % PeriodicFlipActive ) THEN
           DO j=1,n
             k = Indexes(j)
             IF ( k>0 .AND. k<=SIZE(Variable % Perm) ) THEN
               l = Variable % Perm(k)
               IF( l>0 ) THEN
                 Found0 = .TRUE.
                 x(i,j) = Values(Variable % DOFs*(l-1)+i)
                 IF( CurrentModel % Mesh % PeriodicFlip(k) ) x(i,j) = -x(i,j)
               END IF
             END IF
           END DO
         ELSE           
           DO j=1,n
             k = Indexes(j)
             IF ( k>0 .AND. k<=SIZE(Variable % Perm) ) THEN
               l = Variable % Perm(k)
               IF (l>0) THEN
                 Found0 = .TRUE.
                 x(i,j) = Values(Variable % DOFs*(l-1)+i)
               END IF
             END IF
           END DO
         END IF
       ELSE
         DO j=1,n
           IF ( Variable % DOFs*(Indexes(j)-1)+i <= &
               SIZE( Variable % Values ) ) THEN
             Found0 = .TRUE.
             x(i,j) = Values(Variable % DOFs*(Indexes(j)-1)+i)
           END IF
         END DO
       END IF
     END DO
     IF( PRESENT(Found)) Found = Found0
     
  END SUBROUTINE GetVectorLocalSolution


! Eigenmodes may be used as a basis of sensitivity analysis, model reduction etc.
! then these subroutines may be used to obtain the local eigenmodes
!-------------------------------------------------------------------------------

!> Returns the number of eigenmodes
  FUNCTION GetNofEigenModes( name,USolver) RESULT (NofEigenModes)

     CHARACTER(LEN=*), OPTIONAL :: name
     TYPE(Solver_t)  , OPTIONAL, TARGET :: USolver
     INTEGER :: NofEigenModes

     REAL(KIND=dp), POINTER :: Values(:)
     TYPE(Variable_t), POINTER :: Variable
     TYPE(Solver_t)  , POINTER :: Solver

     NofEigenModes = 0

     Solver => CurrentModel % Solver
     IF ( PRESENT(USolver) ) Solver => USolver

     Variable => Solver % Variable
     IF ( PRESENT(name) ) THEN
        Variable => VariableGet( Solver % Mesh % Variables, name )
     END IF

     IF ( .NOT. ASSOCIATED( Variable ) ) RETURN
     IF ( .NOT. ASSOCIATED( Variable % EigenValues ) ) RETURN
     
     NofEigenModes = SIZE( Variable % EigenValues, 1)
  END FUNCTION GetNofEigenModes


!> Returns the desired eigenmode as a scalar field in an element
  SUBROUTINE GetScalarLocalEigenmode( x,name,UElement,USolver,NoEigen,ComplexPart )
     REAL(KIND=dp) :: x(:)
     CHARACTER(LEN=*), OPTIONAL :: name
     TYPE(Solver_t)  , OPTIONAL, TARGET :: USolver
     TYPE(Element_t),  OPTIONAL, TARGET :: UElement
     INTEGER, OPTIONAL :: NoEigen
     LOGICAL, OPTIONAL :: ComplexPart

     COMPLEX(KIND=dp), POINTER :: Values(:)
     TYPE(Variable_t), POINTER :: Variable
     TYPE(Solver_t)  , POINTER :: Solver
     TYPE(Element_t),  POINTER :: Element
     LOGICAL :: IsComplex

     INTEGER :: i, j, n
     INTEGER, POINTER :: Indexes(:)

     Solver => CurrentModel % Solver
     IF ( PRESENT(USolver) ) Solver => USolver

     x = 0.0d0

     Variable => Solver % Variable
     IF ( PRESENT(name) ) THEN
        Variable => VariableGet( Solver % Mesh % Variables, name )
     END IF
     IF ( .NOT. ASSOCIATED( Variable ) ) RETURN
     IF ( .NOT. ASSOCIATED( Variable % EigenVectors ) ) RETURN

     IsComplex = .FALSE.
     IF( PRESENT( ComplexPart) ) IsComplex = ComplexPart

     Element => GetCurrentElement(UElement)

     Indexes => GetIndexStore()
     IF ( ASSOCIATED(Variable % Solver ) ) THEN
       n = GetElementDOFs( Indexes, Element, Variable % Solver )
     ELSE
       n = GetElementDOFs( Indexes, Element, Solver )
     END IF
     n = MIN( n, SIZE(x) )

     IF (SIZE(Variable % EigenVectors,1) < NoEigen) THEN
       CALL Fatal('GetScalarLocalEigenmode', 'Less eigenfunctions available than requested')
     END IF
     Values => Variable % EigenVectors( NoEigen, :)

     IF ( ASSOCIATED( Variable % Perm ) ) THEN
       DO i=1,n
         j = Indexes(i)
         IF ( j>0 .AND. j<= SIZE(Variable % Perm)) THEN
           j = Variable % Perm(j)
           IF ( j>0 ) THEN 
             IF ( IsComplex ) THEN
               x(i) = AIMAG(Values(j))
             ELSE
               x(i) =  REAL(Values(j))
             END IF
           END IF
         END IF
       END DO
     ELSE
       x(1:n) = Values(Indexes(1:n))
     END IF
  END SUBROUTINE GetScalarLocalEigenmode



!> Returns the desired eigenmode as a vector field in an element
  SUBROUTINE GetVectorLocalEigenmode( x,name,UElement,USolver,NoEigen,ComplexPart )
     REAL(KIND=dp) :: x(:,:)
     CHARACTER(LEN=*), OPTIONAL :: name
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     INTEGER, OPTIONAL :: NoEigen
     LOGICAL, OPTIONAL :: ComplexPart

     TYPE(Variable_t), POINTER :: Variable
     TYPE(Solver_t)  , POINTER :: Solver
     TYPE(Element_t),  POINTER :: Element
     LOGICAL :: IsComplex

     INTEGER :: i, j, k, n
     INTEGER, POINTER :: Indexes(:)
     COMPLEX(KIND=dp), POINTER ::  Values(:)

     Solver => CurrentModel % Solver
     IF ( PRESENT(USolver) ) Solver => USolver
     
     IsComplex = .FALSE.
     IF( PRESENT( ComplexPart) ) IsComplex = ComplexPart

     x = 0.0d0

     Variable => Solver % Variable
     IF ( PRESENT(name) ) THEN
        Variable => VariableGet( Solver % Mesh % Variables, name )
     END IF
     IF ( .NOT. ASSOCIATED( Variable ) ) RETURN
     IF ( .NOT. ASSOCIATED( Variable % EigenVectors ) ) RETURN

     Element => GetCurrentElement(UElement)

     Indexes => GetIndexStore()
     IF ( ASSOCIATED(Variable % Solver ) ) THEN
       n = GetElementDOFs( Indexes, Element, Variable % Solver )
     ELSE
       n = GetElementDOFs( Indexes, Element, Solver )
     END IF
     n = MIN( n, SIZE(x) )

     IF (SIZE(Variable % EigenVectors,1) < NoEigen) THEN
       CALL Fatal('GetVectorLocalEigenmode', 'Less eigenfunctions available than requested')
     END IF
     Values => Variable % EigenVectors( NoEigen, : )

     DO i=1,Variable % DOFs
       IF ( ASSOCIATED( Variable % Perm ) ) THEN
         DO j=1,n
           k = Indexes(j)
           IF ( k>0 .AND. k<= SIZE(Variable % Perm)) THEN
             k = Variable % Perm(k)
             IF ( k>0 ) THEN
               IF ( IsComplex ) THEN
                 x(i,j) = AIMAG(Values(Variable % DOFs*(k-1)+i))
               ELSE
                 x(i,j) =  REAL(Values(Variable % DOFs*(k-1)+i))
               END IF
             END IF
           END IF
         END DO
       ELSE
         DO j=1,n
           IF( IsComplex ) THEN
             x(i,j) = AIMAG( Values(Variable % DOFs*(Indexes(j)-1)+i) )
           ELSE
             x(i,j) = REAL( Values(Variable % DOFs*(Indexes(j)-1)+i) )
           END IF
         END DO
       END IF
     END DO
  END SUBROUTINE GetVectorLocalEigenmode
    

  FUNCTION DefaultVariableGet( Name, ThisOnly, USolver ) RESULT ( Var )

    CHARACTER(LEN=*) :: Name
    LOGICAL, OPTIONAL :: ThisOnly
    TYPE(Solver_t), POINTER, OPTIONAL :: USolver    
    TYPE(Variable_t), POINTER :: Var
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Variables

    IF( PRESENT( USolver ) ) THEN
      Variables => USolver % Mesh % Variables
    ELSE
      Variables => CurrentModel % Solver % Mesh % Variables
    END IF
    
    Var => VariableGet( Variables, Name, ThisOnly )
    
  END FUNCTION DefaultVariableGet


!------------------------------------------------------------------------------
!> Add variable to the default variable list.
!------------------------------------------------------------------------------
  SUBROUTINE DefaultVariableAdd( Name, DOFs, Perm, Values,&
      Output,Secondary,VariableType,Global,InitValue,USolver,Var )
    
    CHARACTER(LEN=*) :: Name
    INTEGER, OPTIONAL :: DOFs
    REAL(KIND=dp), OPTIONAL, POINTER :: Values(:)
    LOGICAL, OPTIONAL :: Output
    INTEGER, OPTIONAL, POINTER :: Perm(:)
    LOGICAL, OPTIONAL :: Secondary
    INTEGER, OPTIONAL :: VariableType
    LOGICAL, OPTIONAL :: Global
    REAL(KIND=dp), OPTIONAL :: InitValue   
    TYPE(Solver_t), OPTIONAL, TARGET :: USolver
    TYPE(Variable_t), OPTIONAL, POINTER :: Var
    !------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Variables
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), POINTER :: Solver
    
    IF( PRESENT( USolver ) ) THEN
      Solver => USolver
    ELSE
      Solver => CurrentModel % Solver 
    END IF
    Mesh => Solver % Mesh
    Variables => Mesh % Variables

    CALL VariableAddVector( Variables,Mesh,Solver,Name,DOFs,Values,&
        Perm,Output,Secondary,VariableType,Global,InitValue )

    IF( PRESENT( Var ) ) THEN
      Var => VariableGet( Variables, Name )
    END IF

  END SUBROUTINE DefaultVariableAdd
!------------------------------------------------------------------------------


!> Returns a string by its name if found in the list structure
  FUNCTION GetString( List, Name, Found ) RESULT(str)
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found
     CHARACTER(:), ALLOCATABLE :: str

     str = TRIM(ListGetString(List, Name, Found))
  END FUNCTION GetString


!> Returns an integer by its name if found in the list structure
  FUNCTION GetInteger( List, Name, Found ) RESULT(i)
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found

     INTEGER :: i

     i = ListGetInteger( List, Name, Found )
  END FUNCTION GetInteger


!> Returns a logical flag by its name if found in the list structure, otherwise false
  FUNCTION GetLogical( List, Name, Found, UnfoundFatal, DefValue ) RESULT(l)
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found, UnfoundFatal, DefValue

     LOGICAL :: l

     l = ListGetLogical( List, Name, Found, UnfoundFatal, DefValue )
  END FUNCTION GetLogical


!> Returns a constant real by its name if found in the list structure
  RECURSIVE FUNCTION GetConstReal( List, Name, Found,x,y,z ) RESULT(r)
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found
     REAL(KIND=dp), OPTIONAL :: x,y,z

     REAL(KIND=dp) :: r,xx,yy,zz

     xx = 0.0_dp
     yy = 0.0_dp
     zz = 0.0_dp
     IF ( PRESENT( x ) ) xx = x
     IF ( PRESENT( y ) ) yy = y
     IF ( PRESENT( z ) ) zz = z

     r = ListGetConstReal( List, Name, Found,xx,yy,zz )
  END FUNCTION GetConstReal


!> Returns a real that may depend on global variables such as time, or timestep size, 
!! by its name if found in the list structure
  RECURSIVE FUNCTION GetCReal( List, Name, Found ) RESULT(s)
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found
     INTEGER, TARGET :: Dnodes(1)
     INTEGER, POINTER :: NodeIndexes(:)

     REAL(KIND=dp) :: s
     REAL(KIND=dp), POINTER CONTIG :: x(:)
     TYPE(Element_t), POINTER :: Element

     INTEGER :: n, nthreads, thread, istat

     IF ( PRESENT( Found ) ) Found = .FALSE.

     NodeIndexes => Dnodes
     n = 1
     NodeIndexes(n) = 1

     x => GetValueStore(n)
     x(1:n) = REAL(0, dp)
     IF( ASSOCIATED(List) ) THEN
       IF ( ASSOCIATED(List % Head) ) THEN
         x(1:n) = ListGetReal( List, Name, n, NodeIndexes, Found )
       END IF
     END IF
     s = x(1)
  END FUNCTION GetCReal


!> Returns a real by its name if found in the list structure, and in the active element. 
  RECURSIVE FUNCTION GetReal( List, Name, Found, UElement ) RESULT(x)
    IMPLICIT NONE
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found
     INTEGER, TARGET :: Dnodes(1)
     INTEGER, POINTER :: NodeIndexes(:)
     TYPE(Element_t), OPTIONAL, TARGET :: UElement

     REAL(KIND=dp), POINTER CONTIG :: x(:)
     TYPE(Element_t), POINTER :: Element

     INTEGER :: n, istat

     IF ( PRESENT( Found ) ) Found = .FALSE.

     Element => GetCurrentElement(UElement)

     IF ( ASSOCIATED(Element) ) THEN
       n = GetElementNOFNodes(Element)
       NodeIndexes => Element % NodeIndexes
     ELSE
       n = 1
       NodeIndexes => Dnodes
       NodeIndexes(1) = 1
     END IF

     x => GetValueStore(n)
     x(1:n) = REAL(0, dp)
     IF( ASSOCIATED(List) ) THEN
       IF ( ASSOCIATED(List % Head) ) THEN
         x(1:n) = ListGetReal( List, Name, n, NodeIndexes, Found )
       END IF
     END IF
  END FUNCTION GetReal

  RECURSIVE SUBROUTINE GetRealValues( List, Name, Values, Found, UElement )
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: List
    CHARACTER(LEN=*) :: Name
    REAL(KIND=dp) CONTIG :: Values(:)
    LOGICAL, OPTIONAL :: Found
    TYPE(Element_t), OPTIONAL, TARGET :: UElement
    
    ! Variables
    INTEGER, TARGET :: Dnodes(1)
    INTEGER, POINTER CONTIG :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, istat

    IF ( PRESENT( Found ) ) Found = .FALSE.
    
    Element => GetCurrentElement(UElement)
    
    IF ( ASSOCIATED(Element) ) THEN
      n = GetElementNOFNodes(Element)
      NodeIndexes => Element % NodeIndexes
    ELSE
      n = 1
      NodeIndexes => Dnodes
      NodeIndexes(1) = 1
    END IF
    
    IF( ASSOCIATED(List) ) THEN
      IF ( ASSOCIATED(List % Head) ) THEN
        Values(1:n) = ListGetReal( List, Name, n, NodeIndexes, Found )
      END IF
    END IF
  END SUBROUTINE GetRealValues


!> Returns a material property from either of the parents of the current boundary element
  RECURSIVE FUNCTION GetParentMatProp( Name, UElement, Found, UParent ) RESULT(x)
    CHARACTER(LEN=*) :: Name
    TYPE(Element_t), OPTIONAL, TARGET :: UElement
    LOGICAL, OPTIONAL :: Found
    TYPE(Element_t), OPTIONAL, POINTER :: UParent

    REAL(KIND=dp), POINTER CONTIG :: x(:)
    INTEGER, POINTER :: Indexes(:)
    LOGICAL :: GotIt, GotMat
    INTEGER :: n, leftright, mat_id
    TYPE(ValueList_t), POINTER :: Material
    TYPE(Element_t), POINTER :: Element, Parent
    
    Element => GetCurrentElement(Uelement)

    IF( .NOT. ASSOCIATED( Element ) ) THEN
      CALL Warn('GetParentMatProp','Element not associated!')
    END IF
    
    IF( PRESENT(UParent) ) NULLIFY( UParent )

    n = GetElementNOFNodes(Element)
    Indexes => Element % NodeIndexes

    x => GetValueStore(n)
    x(1:n) = REAL(0, dp)

    IF(.NOT. ASSOCIATED( Element % BoundaryInfo ) ) THEN
      CALL Warn('GetParentMatProp','Boundary element needs parent information!')
      RETURN
    END IF
    
    
    Gotit = .FALSE.
    DO leftright = 1, 2

      IF( leftright == 1) THEN
         Parent => Element % BoundaryInfo % Left
       ELSE 
         Parent => Element % BoundaryInfo % Right
       END IF
       
       IF( ASSOCIATED(Parent) ) THEN

         GotMat = .FALSE.
         IF( Parent % BodyId == 0) THEN
           CYCLE
         ELSE IF( Parent % BodyId <= CurrentModel % NumberOfBodies ) THEN
           mat_id = ListGetInteger( CurrentModel % Bodies(Parent % BodyId) % Values,'Material',GotMat)
         ELSE
           CALL Warn('GetParentMatProp','Invalid parent BodyId '//I2S(Parent % BodyId)//&
               ' for element '//I2S(Parent % ElementIndex))
           CYCLE
         END IF
         
         IF(.NOT. GotMat) THEN
           CALL Warn('GetParentMatProp','Parent body '//I2S(Parent % BodyId)//' does not have material associated!')
         END IF
         
         IF( mat_id > 0 .AND. mat_id <= CurrentModel % NumberOfMaterials ) THEN
           Material => CurrentModel % Materials(mat_id) % Values
         ELSE
           CALL Warn('GetParentMatProp','Material index '//I2S(mat_id)//' not associated to material list')
           CYCLE
         END IF
                             
         IF( .NOT. ASSOCIATED( Material ) ) CYCLE

         IF ( ListCheckPresent( Material,Name) ) THEN
           BLOCK
             TYPE(Element_t), POINTER :: se
             se => CurrentModel % CurrentElement
             CurrentModel % CurrentElement => Element
             x(1:n) = ListGetReal(Material, Name, n, Indexes)
             CurrentModel % CurrentElement => se
           END BLOCK
           IF( PRESENT( UParent ) ) UParent => Parent
           Gotit = .TRUE.
           EXIT
         END IF
       END IF
    END DO
    
    IF( PRESENT( Found ) ) THEN
       Found = GotIt
    ELSE IF(.NOT. GotIt) THEN
       CALL Warn('GetParentMatProp','Property '//TRIM(Name)//' not in either parents!')
    END IF
     
  END FUNCTION GetParentMatProp


!> Returns a constant real array by its name if found in the list structure. 
  RECURSIVE SUBROUTINE GetConstRealArray( List, x, Name, Found, UElement )
     TYPE(ValueList_t), POINTER :: List
     REAL(KIND=dp), POINTER :: x(:,:)
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found
     TYPE(Element_t), OPTIONAL, TARGET :: UElement

     IF ( PRESENT( Found ) ) Found = .FALSE.
     IF(ASSOCIATED(List)) THEN
       IF ( ASSOCIATED(List % Head) ) THEN
         x => ListGetConstRealArray( List, Name, Found )
       END IF
     END IF
  END SUBROUTINE GetConstRealArray

!> Returns a real array by its name if found in the list structure, and in the active element. 
  RECURSIVE SUBROUTINE GetRealArray( List, x, Name, Found, UElement )
     REAL(KIND=dp), POINTER :: x(:,:,:)
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found
     TYPE(Element_t), OPTIONAL, TARGET :: UElement

     TYPE(Element_t), POINTER :: Element

     INTEGER :: n

     IF ( PRESENT( Found ) ) Found = .FALSE.

     Element => GetCurrentElement(UElement)

     n = GetElementNOFNodes( Element )
     IF ( ASSOCIATED(List) ) THEN
       IF ( ASSOCIATED(List % Head) ) THEN
         CALL ListGetRealArray( List, Name, x, n, Element % NodeIndexes, Found )
       END IF
     END IF
  END SUBROUTINE GetRealArray

!> Returns a real vector by its name if found in the list structure, and in the active element. 
  RECURSIVE SUBROUTINE GetRealVector( List, x, Name, Found, UElement )
     REAL(KIND=dp) :: x(:,:)
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found
     TYPE(Element_t), OPTIONAL, TARGET :: UElement

     TYPE(Element_t), POINTER :: Element

     INTEGER :: n

     x = 0._dp
     IF ( PRESENT( Found ) ) Found = .FALSE.

     Element => GetCurrentElement(UElement)

     n = GetElementNOFNodes( Element )
     IF ( ASSOCIATED(List) ) THEN
       IF ( ASSOCIATED(List % Head) ) THEN
         CALL ListGetRealvector( List, Name, x, n, Element % NodeIndexes, Found )
       END IF
     END IF
  END SUBROUTINE GetRealVector

!> Returns a complex vector by its name if found in the list structure, and in the active element. 
  RECURSIVE SUBROUTINE GetComplexVector( List, x, Name, Found, UElement )
     COMPLEX(KIND=dp) :: x(:,:)
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found
     TYPE(Element_t), OPTIONAL, TARGET :: UElement

     TYPE(Element_t), POINTER :: Element
     LOGICAL :: lFound
     INTEGER :: n
     REAL(KIND=dp), ALLOCATABLE :: xr(:,:)

     x = 0._dp
     IF ( PRESENT( Found ) ) Found = .FALSE.

     Element => GetCurrentElement(UElement)

     n = GetElementNOFNodes( Element )
     IF ( ASSOCIATED(List) ) THEN
       IF ( ASSOCIATED(List % Head) ) THEN
          ALLOCATE(xr(SIZE(x,1),SIZE(x,2)))
          CALL ListGetRealvector( List, Name, xr, n, &
                 Element % NodeIndexes, lFound )
          IF(PRESENT(Found)) Found=lFound
          x = xr
          CALL ListGetRealvector( List, TRIM(Name)//' im', &
              xr, n, Element % NodeIndexes, lFound )
          IF(PRESENT(Found)) Found=Found.OR.lFound
          x = CMPLX(REAL(x), xr)
       END IF
     END IF
  END SUBROUTINE GetComplexVector


!> Set a named elementwise property (real-valued) to the active element or
!> given element
  SUBROUTINE SetElementProperty( Name, Values, UElement )
    CHARACTER(LEN=*) :: Name
    REAL(KIND=dp) :: Values(:)
    TYPE(Element_t), POINTER, OPTIONAL :: UElement

    TYPE(ElementData_t), POINTER :: p

    TYPE(Element_t), POINTER :: Element

    Element => GetCurrentElement(UElement)

    p => Element % PropertyData
    DO WHILE( ASSOCIATED(p) )
      IF ( Name==p % Name ) EXIT
      p => p % Next
    END DO

    IF ( ASSOCIATED(p) ) THEN
      IF ( SIZE(P % Values) == SIZE(Values) ) THEN
        P % Values = Values
      ELSE
        DEALLOCATE( P % Values )
        ALLOCATE( P % Values(SIZE(Values)) )
        P % Values = Values
      END IF
    ELSE
      ALLOCATE(p)
      ALLOCATE( P % Values(SIZE(Values)) )
      p % Values = Values
      p % Name = Name
      p % Next => Element % PropertyData
      Element % PropertyData => p
    END IF
  END SUBROUTINE SetElementProperty

!> Get a named elementwise property (real-valued) from the active element or 
!> from a given element
  FUNCTION GetElementProperty( Name, UElement ) RESULT(Values)
    CHARACTER(LEN=*) :: Name
    REAL(KIND=dp), POINTER :: Values(:)
    TYPE(Element_t), POINTER, OPTIONAL :: UElement

    TYPE(ElementData_t), POINTER :: p

    TYPE(Element_t), POINTER :: Element

    Element => GetCurrentElement(UElement)

    Values => NULL()
    p=> Element % PropertyData

    DO WHILE( ASSOCIATED(p) )
      IF ( Name==p % Name ) THEN
        Values => p % Values
        RETURN
      END IF
      p => p % Next
    END DO
  END FUNCTION GetElementProperty


!> Get a handle to the active element from the list of all active elements
  FUNCTION GetActiveElement(t,USolver) RESULT(Element)
     INTEGER :: t
     TYPE(Element_t), POINTER :: Element
     TYPE( Solver_t ), OPTIONAL, TARGET :: USolver

     TYPE( Solver_t ), POINTER :: Solver
     INTEGER :: ind

     Solver => CurrentModel % Solver
     IF ( PRESENT( USolver ) ) Solver => USolver

     IF ( t > 0 .AND. t <= Solver % NumberOfActiveElements ) THEN
       ! Check if colouring is really used by the solver
       IF( Solver % CurrentColour > 0 .AND. &
               ASSOCIATED( Solver % ColourIndexList ) ) THEN
         ind = Solver % ActiveElements( &
                 Solver % ColourIndexList % ind(&
                 Solver % ColourIndexList % ptr(Solver % CurrentColour)+(t-1) ) )
       ELSE
         ind = Solver % ActiveElements(t)
       END IF

       Element => Solver % Mesh % Elements( ind )
 
#ifdef _OPENMP
       IF (omp_in_parallel()) THEN
         CurrentElementThread => Element
       ELSE
         ! May be used by user functions, not thread safe
         CurrentModel % CurrentElement => Element 
       END IF
#else
       ! May be used by user functions, not thread safe
       CurrentModel % CurrentElement => Element 
#endif
     ELSE
       WRITE( Message, * ) 'Invalid element number requested: ', t
       CALL Fatal( 'GetActiveElement', Message )
     END IF
  END FUNCTION GetActiveElement


!> Get a handle to a boundary element from the list of all boundary elements
  FUNCTION GetBoundaryElement(t,USolver) RESULT(Element)
     INTEGER :: t
     TYPE(Element_t), POINTER :: Element
     TYPE( Solver_t ), OPTIONAL, TARGET :: USolver
     TYPE( Solver_t ), POINTER :: Solver
     INTEGER :: ind

     Solver => CurrentModel % Solver
     IF ( PRESENT( USolver ) ) Solver => USolver

     IF ( t > 0 .AND. t <= Solver % Mesh % NumberOfBoundaryElements ) THEN
       ! Check if colouring is really used by the solver
       IF( Solver % CurrentBoundaryColour > 0 .AND. &
            ASSOCIATED( Solver % BoundaryColourIndexList ) ) THEN
          ind = Solver % BoundaryColourIndexList % ind( &
               Solver % BoundaryColourIndexList % ptr(Solver % CurrentBoundaryColour)+(t-1))
       ELSE
         ind = t
       END IF

       ! Element => Solver % Mesh % Elements( Solver % Mesh % NumberOfBulkElements+t )
       Element => Solver % Mesh % Elements( Solver % Mesh % NumberOfBulkElements + ind )
#ifdef _OPENMP
        IF (omp_in_parallel()) THEN
          ! May be used by user functions, thread safe
          CurrentElementThread => Element
        ELSE
          CurrentModel % CurrentElement => Element
        END IF
#else
        CurrentModel % CurrentElement => Element
#endif
     ELSE
        WRITE( Message, * ) 'Invalid element number requested: ', t
        CALL Fatal( 'GetBoundaryElement', Message )
     END IF
  END FUNCTION GetBoundaryElement


!> Check if the boundary element is active in the current solve 
  FUNCTION ActiveBoundaryElement(UElement,USolver,DGBoundary) RESULT(l)
     TYPE(Element_t), OPTIONAL,  TARGET :: UElement
     TYPE(Solver_t),  OPTIONAL,  TARGET :: USolver
     LOGICAL, OPTIONAL :: DGBoundary

     LOGICAL :: l, DGb
     INTEGER :: n, n2
     INTEGER, POINTER :: Indexes(:)

     TYPE( Solver_t ), POINTER :: Solver
     TYPE(Element_t), POINTER :: Element, P1, P2

     Solver => CurrentModel % Solver
     IF ( PRESENT( USolver ) ) Solver => USolver

     Element => GetCurrentElement(UElement)

     Indexes => GetIndexStore()
     n = GetElementDOFs( Indexes, Element, Solver )

     DGb = Solver % DG .AND. PRESENT(DGboundary)
     IF(DGb) DGb = DGboundary

     IF (DGb) THEN
       P1 => Element % BoundaryInfo % Left
       P2 => Element % BoundaryInfo % Right
       IF ( ASSOCIATED(P1).AND.ASSOCIATED(P2) ) THEN
         n = P1 % Type % NumberOfNodes
         l = ALL(Solver % Variable % Perm(Indexes(1:n)) > 0)
         IF (.NOT.l) THEN
           n2 = P2 % Type % NumberOfNodes
           l = ALL(Solver % Variable % Perm(Indexes(n+1:n+n2)) > 0)
          END IF
       ELSE
         l = ALL(Solver % Variable % Perm(Indexes(1:n)) > 0)
       END IF
     ELSE
       IF (isActivePElement(Element)) n=GetElementNOFNOdes(Element)
       l = ALL(Solver % Variable % Perm(Indexes(1:n)) > 0)
     END IF
  END FUNCTION ActiveBoundaryElement


!> Return the element code in Elmer convention of the active element
  FUNCTION GetElementCode( Element )  RESULT(etype)
     INTEGER :: etype
     TYPE(Element_t), OPTIONAL :: Element
     TYPE(Element_t), POINTER :: CurrElement

     CurrElement => GetCurrentElement(Element)
     etype = CurrElement % TYPE % ElementCode
  END FUNCTION GetElementCode

!> Return the element dimension in Elmer convention of the active element
  FUNCTION GetElementDim( Element )  RESULT(edim)
    INTEGER :: edim
    TYPE(Element_t), OPTIONAL :: Element
    TYPE(Element_t), POINTER :: CurrElement
    INTEGER :: etype
    
    CurrElement => GetCurrentElement(Element)
    etype = CurrElement % TYPE % ElementCode
    IF( etype >= 500 ) THEN
      edim = 3
    ELSE IF( etype >= 300 ) THEN
      edim = 2
    ELSE IF( etype >= 200 ) THEN
      edim = 1
    ELSE 
      edim = 0
    END IF          
  END FUNCTION GetElementDim


!> Return the element family in Elmer convention of the active element
  FUNCTION GetElementFamily( Element )  RESULT(family)
     INTEGER :: family
     TYPE(Element_t), OPTIONAL :: Element
     TYPE(Element_t), POINTER :: CurrElement

     CurrElement => GetCurrentElement(Element)
     family = CurrElement % TYPE % ElementCode / 100
  END FUNCTION GetElementFamily


!> Return the number of corners nodes i.e. the number of dofs for the lowest order element
  FUNCTION GetElementCorners( Element )  RESULT(corners)
    INTEGER :: corners
    TYPE(Element_t), OPTIONAL :: Element
    TYPE(Element_t), POINTER :: CurrElement
    
    CurrElement => GetCurrentElement(Element)
    corners = CurrElement % TYPE % ElementCode / 100
    IF( corners >= 5 .AND. corners <= 7 ) THEN
      corners = corners - 1
    END IF
  END FUNCTION GetElementCorners
  
!> Return true if the element is a possible flux element
!> Needed to skip nodal elements in 2D and 3D boundary condition setting.
  FUNCTION PossibleFluxElement( Element, Mesh )  RESULT(possible)
     LOGICAL :: possible
     TYPE(Element_t), OPTIONAL :: Element
     TYPE(Mesh_t), OPTIONAL :: Mesh
     INTEGER :: MeshDim, family

     ! Orphan elements are not currently present in the mesh so any
     ! boundary condition that exists is a possible flux element also.
     ! Thus this routine is more or less obsolete.
     possible = .TRUE.

     RETURN

     
     IF( PRESENT( Mesh ) ) THEN
       MeshDim = Mesh % MeshDim
     ELSE
       MeshDim = CurrentModel % Solver % Mesh % MeshDim
     END IF

     family = GetElementFamily( Element )
     
     ! This is not a generic rule but happens to be true for all combinations
     ! 3D: families 3 and 4
     ! 2D: family 2
     ! 1D: family 1
     possible = ( MeshDim <= family )

   END FUNCTION PossibleFluxElement


!> Return the number of nodes in the active element
  FUNCTION GetElementNOFNodes( Element ) RESULT(n)
     INTEGER :: n
     TYPE(Element_t), OPTIONAL :: Element
     TYPE(Element_t), POINTER :: CurrElement

     CurrElement => GetCurrentElement(Element)
     n = CurrElement % TYPE % NumberOfNodes
  END FUNCTION GetELementNOFNodes


!> Return the number of element degrees of freedom 
  FUNCTION GetElementNOFDOFs( UElement,USolver ) RESULT(n)

     USE PElementMaps, ONLY : isActivePElement, getEdgeDOFs, getFaceDOFs, getBubbleDOFs

     INTEGER :: n
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement

     TYPE(Element_t), POINTER :: Element, Face
     TYPE(Solver_t),  POINTER :: Solver

     INTEGER :: i, j, k, id, ElemFamily, ParentFamily, face_type, face_id 
     INTEGER :: NDOFs
     LOGICAL :: Found, GB, NeedEdges

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
     ELSE
        Solver => CurrentModel % Solver
     END IF

     n = 0

     IF (.NOT. ASSOCIATED(Solver)) THEN
       CALL Warn('GetElementNOFDOFS', &
           'Cannot return the number of DOFs without knowing solver')
       RETURN
     END IF

     Element => GetCurrentElement( UElement )
     ElemFamily = GetElementFamily(Element)

     IF( Solver % DG ) THEN
       n = Element % DGDOFs
       IF ( n>0 ) RETURN
     END IF

     id = Element % BodyId
     IF ( Id==0 .AND. ASSOCIATED(Element % BoundaryInfo) ) THEN
       IF ( ASSOCIATED(Element % BoundaryInfo % Left) ) &
         id = Element % BoundaryInfo % Left % BodyId

       IF ( ASSOCIATED(Element % BoundaryInfo % Right) ) &
         id = Element % BoundaryInfo % Right % BodyId
     END IF
     IF ( Id==0 ) id=1

     IF ( Solver % Def_Dofs(ElemFamily,id,1)>0 ) n = Element % NDOFs
     NDOFs = MAX(0, Solver % Def_Dofs(ElemFamily,id,1))
     IF (NDOFs > 0) n = NDOFs * Element % TYPE % NumberOfNodes

     NeedEdges = .FALSE.
     DO i=2,SIZE(Solver % Def_Dofs,3)
       IF (Solver % Def_Dofs(ElemFamily, id, i)>=0) THEN
         NeedEdges = .TRUE.
         EXIT
       END IF
     END DO

     IF (.NOT. NeedEdges) THEN
       !
       ! Check whether face DOFs have been generated by "-quad_face b: ..." or
       ! "-tri_face b: ..."
       !
       IF (ElemFamily == 3 .OR. ElemFamily == 4) THEN
         IF (Solver % Def_Dofs(6+ElemFamily, id, 5)>=0) NeedEdges = .TRUE.
       ELSE
         !
         ! Check finally if 3-D faces are associated with face bubbles
         !
         IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
           DO j=1,Element % TYPE % NumberOfFaces
             Face => Solver % Mesh % Faces(Element % FaceIndexes(j))
             face_type = Face % TYPE % ElementCode/100
             IF (ASSOCIATED(Face % BoundaryInfo % Left)) THEN
               face_id  = Face % BoundaryInfo % Left % BodyId
               k = MAX(0,Solver % Def_Dofs(face_type+6,face_id,5))
             END IF
             IF (ASSOCIATED(Face % BoundaryInfo % Right)) THEN
               face_id = Face % BoundaryInfo % Right % BodyId
               k = MAX(k,Solver % Def_Dofs(face_type+6,face_id,5))
             END IF
             IF (k > 0) THEN
               NeedEdges = .TRUE.
               EXIT
             END IF
           END DO
         END IF
       END IF
     END IF

     IF ( .NOT. NeedEdges ) RETURN


     BLOCK
       LOGICAL :: EdgesDone, FacesDone
       INTEGER :: Ind, i,j, p, nb, EDOFs, FDOFs, BDOFs
       INTEGER :: face_id
       TYPE(Element_t), POINTER :: Parent, Edge
  
       EdgesDone = .FALSE.; FacesDone = .FALSE.
       IF ( ASSOCIATED( Element % EdgeIndexes ) ) THEN
         DO j=1,Element % Type % NumberOFEdges
           Edge => Solver % Mesh % Edges( Element % EdgeIndexes(j) )
           IF (Edge % Type % ElementCode == Element % Type % ElementCode) THEN
             IF (.NOT. Solver % GlobalBubbles.OR..NOT.ASSOCIATED(Element % BoundaryInfo)) CYCLE
           END IF

           EDOFs = 0 
           IF (Solver % Def_Dofs(ElemFamily,id,2) >= 0) THEN
             EDOFs = Solver % Def_Dofs(ElemFamily,id,2)
           ELSE IF (Solver % Def_Dofs(ElemFamily,id,6) > 1) THEN
! TO DO: This is not yet perfect; cf. what is done in InitialPermutation
             EDOFs = getEdgeDOFs(Element, Solver % Def_Dofs(ElemFamily,id,6))
           END IF
           n = n + EDOFs
         END DO
         EdgesDone = .TRUE.
       END IF

       IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
         DO j=1,Element % TYPE % NumberOfFaces
           Face => Solver % Mesh % Faces( Element % FaceIndexes(j) )

           IF (Face % Type % ElementCode==Element % Type % ElementCode) THEN
             IF ( .NOT.Solver % GlobalBubbles.OR..NOT.ASSOCIATED(Element % BoundaryInfo)) CYCLE
           END IF

           k = MAX(0,Solver % Def_Dofs(ElemFamily,id,3))
           IF (k == 0) THEN
             !
             ! NOTE: This depends on what dofs have been introduced
             ! by using the construct "-quad_face b: ..." and
             ! "-tri_face b: ..."
             !
             face_type = Face % TYPE % ElementCode/100
             IF (ASSOCIATED(Face % BoundaryInfo % Left)) THEN
               face_id  = Face % BoundaryInfo % Left % BodyId
               k = MAX(0,Solver % Def_Dofs(face_type+6,face_id,5))
             END IF
             IF (ASSOCIATED(Face % BoundaryInfo % Right)) THEN
               face_id = Face % BoundaryInfo % Right % BodyId
               k = MAX(k,Solver % Def_Dofs(face_type+6,face_id,5))
             END IF

             FDOFs = 0
             IF (k > 0) THEN
               FDOFs = k
             ELSE IF (Solver % Def_Dofs(ElemFamily,id,6) > 1) THEN
! TO DO: This is not yet perfect; cf. what is done in InitialPermutation
               FDOFs = getFaceDOFs(Element,Solver % Def_Dofs(ElemFamily,id,6),j,Face)
             END IF
           END IF
           n = n + FDOFs
         END DO
         FacesDone = .TRUE.
       END IF

       IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN

         IF (isActivePElement(Element) ) THEN
           Parent => Element % PDefs % LocalParent
         ELSE
           Parent => Element % BoundaryInfo % Left
           IF (.NOT.ASSOCIATED(Parent) ) &
               Parent => Element % BoundaryInfo % Right
         END IF
         IF (.NOT.ASSOCIATED(Parent) ) RETURN
         ParentFamily = Parent % TYPE % ElementCode / 100

         SELECT CASE(ElemFamily)
         CASE(2)
           IF ( .NOT. EdgesDone .AND. ASSOCIATED(Parent % EdgeIndexes) ) THEN
             EDOFs =  0
             IF ( isActivePElement(Element, Solver) ) THEN
               Ind=Element % PDefs % LocalNumber
             ELSE
               DO Ind=1,Parent % TYPE % NumberOfEdges
                 Edge => Solver % Mesh % Edges(Parent % EdgeIndexes(ind))
                 k = 0
                 DO i=1,Edge % TYPE % NumberOfNodes
                   DO j=1,Element % TYPE % NumberOfNodes
                     IF ( Edge % NodeIndexes(i)==Element % NodeIndexes(j) ) k=k+1
                   END DO
                 END DO
                 IF ( k==Element % TYPE % NumberOfNodes) EXIT
               END DO
             END IF

             IF (Solver % Def_Dofs(ElemFamily,id,2) >= 0) THEN
               EDOFs = Solver % Def_Dofs(ElemFamily,id,2)
             ELSE IF (Solver % Def_Dofs(ElemFamily,id,6) > 1) THEN
               EDOFs = getEdgeDOFs(Parent, Solver % Def_Dofs(ParentFamily,id,6))
             END IF

             n = n + EDOFs
           END IF

         CASE(3,4)
           IF ( .NOT. FacesDone .AND. ASSOCIATED( Parent % FaceIndexes ) ) THEN

             IF ( isActivePElement(Element, Solver) ) THEN
               Ind=Element % PDefs % LocalNumber
             ELSE
               DO Ind=1,Parent % TYPE % NumberOfFaces
                 Face => Solver % Mesh % Faces(Parent % FaceIndexes(ind))
                 k = 0
                 DO i=1,Face % TYPE % NumberOfNodes
                   DO j=1,Element % TYPE % NumberOfNodes
                     IF ( Face % NodeIndexes(i)==Element % NodeIndexes(j)) k=k+1
                   END DO
                 END DO
                 IF ( k==Face % TYPE % NumberOfNodes) EXIT
               END DO
             END IF
               
             IF (Ind >= 1 .AND. Ind <= Parent % Type % NumberOfFaces) THEN

               IF (ASSOCIATED(Element % FaceIndexes).AND. isActivePelement(Element, Solver) ) THEN
                 Face => Solver % Mesh % Faces(Element % PDefs % localParent % Faceindexes(Ind))
               ELSE
                 Face => Element
               END IF

               IF (.NOT.EdgesDone .AND. ASSOCIATED(Face % EdgeIndexes)) THEN
                 DO j=1,Face % TYPE % NumberOFEdges
                   Edge => Solver % Mesh % Edges(Face % EdgeIndexes(j))

                   EDOFs = 0
                   IF (Solver % Def_Dofs(ElemFamily,id,2) >= 0) THEN
                     EDOFs = Solver % Def_Dofs(ElemFamily,id,2)
                   ELSE IF (Solver % Def_Dofs(ElemFamily,id,6) > 1) THEN
! TO DO: This is not yet perfect when p varies over mesh; cf. what is done in InitialPermutation
                     EDOFs = getEdgeDOFs(Element, Solver % Def_Dofs(ElemFamily,id,6))
                   END IF
                   n = n + EDOFs
                 END DO
               END IF

               FDOFs = 0
               IF (Solver % Def_Dofs(ParentFamily,id,6) > 1) THEN
                 FDOFs = getFaceDOFs(Parent,Solver % Def_Dofs(ParentFamily,id,6),Ind,Face)
               ELSE
                 k = MAX(0,Solver % Def_Dofs(ElemFamily,id,3))
                 IF (k == 0) THEN
                   !
                   ! NOTE: This depends on what dofs have been introduced
                   ! by using the construct "-quad_face b: ..." and
                   ! "-tri_face b: ..."
                   !
                   face_type = Face % TYPE % ElementCode/100
                   IF (ASSOCIATED(Face % BoundaryInfo % Left)) THEN
                     face_id  = Face % BoundaryInfo % Left % BodyId
                     k = MAX(0,Solver % Def_Dofs(face_type+6,face_id,5))
                   END IF
                   IF (ASSOCIATED(Face % BoundaryInfo % Right)) THEN
                     face_id = Face % BoundaryInfo % Right % BodyId
                     k = MAX(k,Solver % Def_Dofs(face_type+6,face_id,5))
                   END IF
                 END IF

                 IF (k > 0) THEN
                   FDOFs = k
                 END IF
               END IF
               n = n + FDOFs
             END IF
           END IF
         END SELECT
       ELSE
         IF (Solver % GlobalBubbles .AND. ASSOCIATED(Element % BubbleIndexes)) THEN
           BDOFs = 0
           nb = Solver % Def_Dofs(ElemFamily,id,5)
           p = Solver % Def_Dofs(ElemFamily,id,6)
           IF (nb >= 0 .OR. p >=1) THEN
             IF (p > 1) BDOFs = GetBubbleDOFs(Element, p)
             BDOFs = MAX(nb, BDOFs)
           ELSE
             ! The following is not an ideal way to obtain the bubble count
             ! in order to support solverwise definitions, but we are not expected 
             ! to end up in this branch anyway:
             BDOFs = Element % BDOFs
           END IF
           n = n + BDOFs
         END IF
       END IF
     END BLOCK
  END FUNCTION GetElementNOFDOFs


!> In addition to returning the number of degrees of freedom associated with 
!> the element, the indexes of the degrees of freedom are also returned.
!-------------------------------------------------------------------------
  FUNCTION GetElementDOFs( Indexes, UElement, USolver, NotDG )  RESULT(NB)
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     INTEGER :: Indexes(:)
     LOGICAL, OPTIONAL  ::  NotDG
     INTEGER :: NB
     
     NB = mGetElementDOFs( Indexes, UElement, USolver, NotDG )
  END FUNCTION GetElementDOFs


! -----------------------------------------------------------------------------
!> Returns the number of bubble degrees of freedom in the active element.
!> If the sif file contains more than one solver section
!> with each of them having their own specification of the "Element" 
!> keyword, the BDOFs field of the Element structure may not be the number of 
!> bubbles that should be assigned to the solver. With the optional argument 
!> Update = .TRUE., the correct solver-wise bubble count is assigned to the 
!> the Element structure.
! -----------------------------------------------------------------------------
  FUNCTION GetElementNOFBDOFs( Element, USolver, Update ) RESULT(n)
! -----------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Element_t), OPTIONAL :: Element
    TYPE(Solver_t), OPTIONAL, POINTER :: USolver
    LOGICAL, OPTIONAL :: Update

    TYPE(Element_t), POINTER  :: CurrElement
    TYPE(Solver_t), POINTER :: Solver
    LOGICAL :: Found, GB, UpdateRequested
    INTEGER :: k, p, ElemFamily

    IF ( PRESENT( USolver ) ) THEN
       Solver => USolver
    ELSE
       Solver => CurrentModel % Solver
    END IF

    UpdateRequested = .FALSE.
    IF ( PRESENT(Update) ) UpdateRequested = Update

    !GB = ListGetLogical( Solver % Values, 'Bubbles in Global System', Found )
    !IF (.NOT.Found) GB = .TRUE.
    GB = Solver % GlobalBubbles

    n = 0
    IF ( .NOT. GB ) THEN
      CurrElement => GetCurrentElement(Element)
      ElemFamily = GetElementFamily(CurrElement)

      k = Solver % Def_Dofs(ElemFamily, CurrElement % Bodyid, 5) 
      p = Solver % Def_Dofs(ElemFamily, CurrElement % Bodyid, 6) 

      IF (k >= 0 .OR. p >= 1) THEN
        IF (p > 1) n = GetBubbleDOFs(CurrElement, p)
        n = MAX(k,n)
      ELSE 
        n = CurrElement % BDOFs
      END IF

      IF (UpdateRequested .AND. n>=0) THEN
        CurrElement % BDOFs = n
      END IF
    ELSE
      ! Rectify the bubble count assigned to the Element argument in case
      ! some other solver has tampered it:
      IF (UpdateRequested) THEN
        CurrElement => GetCurrentElement(Element)
        ElemFamily = GetElementFamily(CurrElement)

        k = Solver % Def_Dofs(ElemFamily, CurrElement % Bodyid, 5) 
        p = Solver % Def_Dofs(ElemFamily, CurrElement % Bodyid, 6) 

        IF (k >= 0 .OR. p >= 1) THEN
          IF (p > 1) n = GetBubbleDOFs(CurrElement, p)
          n = MAX(k,n)
          IF ( n>=0 ) CurrElement % BDOFs = n
        END IF
        n = 0
      END IF
    END IF
  END FUNCTION GetElementNOFBDOFs


!> Returns the nodal coordinate values in the active element
  SUBROUTINE GetElementNodes( ElementNodes, UElement, USolver, UMesh )
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Solver_t), OPTIONAL, TARGET :: USolver
     TYPE(Mesh_t), OPTIONAL, TARGET :: UMesh
     TYPE(Element_t), OPTIONAL, TARGET :: UElement

     INTEGER :: i,n,nd,sz,sz1
     INTEGER, POINTER :: Indexes(:)
     TYPE(Mesh_t),  POINTER  :: Mesh
     TYPE(Element_t), POINTER :: Element

     Element => GetCurrentElement(UElement)

     IF( PRESENT( UMesh ) ) THEN
       Mesh => UMesh
     ELSE IF( PRESENT( USolver ) ) THEN
       Mesh => USolver % Mesh
     ELSE
       Mesh => CurrentModel % Solver % Mesh
     END IF
             
     n = MAX(Mesh % MaxElementNodes,Mesh % MaxElementDOFs)

     IF ( .NOT. ASSOCIATED( ElementNodes % x ) ) THEN
       ALLOCATE( ElementNodes % x(n), ElementNodes % y(n),ElementNodes % z(n) )
     ELSE IF ( SIZE(ElementNodes % x)<n ) THEN
       DEALLOCATE(ElementNodes % x, ElementNodes % y, ElementNodes % z)
       ALLOCATE( ElementNodes % x(n), ElementNodes % y(n),ElementNodes % z(n) )
     END IF

     n = Element % TYPE % NumberOfNodes
     ElementNodes % x(1:n) = Mesh % Nodes % x(Element % NodeIndexes)
     ElementNodes % y(1:n) = Mesh % Nodes % y(Element % NodeIndexes)
     ElementNodes % z(1:n) = Mesh % Nodes % z(Element % NodeIndexes)

     sz = SIZE(ElementNodes % x)
     IF ( sz > n ) THEN
       ElementNodes % x(n+1:sz) = 0.0_dp
       ElementNodes % y(n+1:sz) = 0.0_dp
       ElementNodes % z(n+1:sz) = 0.0_dp
     END IF

     sz1 = SIZE(Mesh % Nodes % x)
     IF (sz1 > Mesh % NumberOfNodes) THEN
        Indexes => GetIndexStore()
        nd = GetElementDOFs(Indexes,Element,NotDG=.TRUE.)
        DO i=n+1,nd
           IF ( Indexes(i)>0 .AND. Indexes(i)<=sz1 ) THEN
             ElementNodes % x(i) = Mesh % Nodes % x(Indexes(i))
             ElementNodes % y(i) = Mesh % Nodes % y(Indexes(i))
             ElementNodes % z(i) = Mesh % Nodes % z(Indexes(i))
           END IF
        END DO
     END IF
  END SUBROUTINE GetElementNodes


  ! This is just a small wrapper in case we want to get the original and not the
  ! mapped coordinates. This assumes that the original coordinates are stored in
  ! NodesOrig. This is rarely need hence no reason to overload the standard routine
  ! with this baggage.
  !---------------------------------------------------------------------------------
  SUBROUTINE GetElementNodesOrig( ElementNodes, UElement, USolver, UMesh )
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Solver_t), OPTIONAL, TARGET :: USolver
     TYPE(Mesh_t), OPTIONAL, TARGET :: UMesh
     TYPE(Element_t), OPTIONAL, TARGET :: UElement

     TYPE(Mesh_t),  POINTER  :: Mesh
     TYPE(Nodes_t), POINTER :: TmpNodes

     IF( PRESENT( UMesh ) ) THEN
       Mesh => UMesh
     ELSE IF( PRESENT( USolver ) ) THEN
       Mesh => USolver % Mesh
     ELSE
       Mesh => CurrentModel % Solver % Mesh
     END IF

     TmpNodes => Mesh % Nodes
     IF(.NOT. ASSOCIATED( Mesh % NodesOrig ) ) THEN
       CALL Fatal('GetElementNodesOrig','Original node coordinates not yet stored!')
     END IF
     Mesh % Nodes => Mesh % NodesOrig

     CALL GetElementNodes( ElementNodes, UElement, Umesh = Mesh )
     Mesh % Nodes => TmpNodes

   END SUBROUTINE GetElementNodesOrig

  
!> Returns the nodal coordinate values in the active element
    SUBROUTINE GetElementNodesVec( ElementNodes, UElement, USolver, UMesh )
        TYPE(Nodes_t), TARGET :: ElementNodes
        TYPE(Solver_t), OPTIONAL, TARGET :: USolver
        TYPE(Mesh_t), OPTIONAL, TARGET :: UMesh
        TYPE(Element_t), OPTIONAL, TARGET :: UElement

        INTEGER :: padn, dum

        INTEGER :: i,n,nd,sz,sz1
        INTEGER, POINTER CONTIG :: Indexes(:)

        TYPE(Solver_t),  POINTER  :: Solver
        TYPE(Mesh_t),  POINTER  :: Mesh
        TYPE(Element_t), POINTER :: Element

        Element => GetCurrentElement(UElement)

        IF( PRESENT( UMesh ) ) THEN
          Mesh => UMesh
        ELSE IF( PRESENT( USolver ) ) THEN
          Mesh => USolver % Mesh
        ELSE
          Mesh => CurrentModel % Solver % Mesh
        END IF

        n = MAX(Mesh % MaxElementNodes,Mesh % MaxElementDOFs)
        padn = n
        
        ! Here we could pad beginning of columns of xyz to 64-byte 
        ! boundaries if needed as follows
        ! padn=NBytePad(n,STORAGE_SIZE(REAL(1,dp))/8,64)
        
        IF (.NOT. ALLOCATED( ElementNodes % xyz)) THEN
            IF (ASSOCIATED(ElementNodes % x)) DEALLOCATE(ElementNodes % x) 
            IF (ASSOCIATED(ElementNodes % y)) DEALLOCATE(ElementNodes % y) 
            IF (ASSOCIATED(ElementNodes % z)) DEALLOCATE(ElementNodes % z) 
          
            ALLOCATE(ElementNodes % xyz(padn,3))
            ElementNodes % xyz = REAL(0,dp)
            ElementNodes % x => ElementNodes % xyz(1:n,1)
            ElementNodes % y => ElementNodes % xyz(1:n,2)
            ElementNodes % z => ElementNodes % xyz(1:n,3)
        ELSE IF (SIZE(ElementNodes % xyz, 1)<padn) THEN
            DEALLOCATE(ElementNodes % xyz)
            ALLOCATE(ElementNodes % xyz(padn,3))
            ElementNodes % xyz = REAL(0,dp)
            ElementNodes % x => ElementNodes % xyz(1:n,1)
            ElementNodes % y => ElementNodes % xyz(1:n,2)
            ElementNodes % z => ElementNodes % xyz(1:n,3)
        ELSE
            ElementNodes % x => ElementNodes % xyz(1:n,1)
            ElementNodes % y => ElementNodes % xyz(1:n,2)
            ElementNodes % z => ElementNodes % xyz(1:n,3)
        END IF

        n = Element % TYPE % NumberOfNodes
!DIR$ IVDEP
        DO i=1,n
          ElementNodes % x(i) = Mesh % Nodes % x(Element % NodeIndexes(i))
          ElementNodes % y(i) = Mesh % Nodes % y(Element % NodeIndexes(i))
          ElementNodes % z(i) = Mesh % Nodes % z(Element % NodeIndexes(i))
        END DO

        sz = SIZE(ElementNodes % xyz,1)
        IF ( sz > n ) THEN
            ElementNodes % xyz(n+1:sz,1) = 0.0d0
            ElementNodes % xyz(n+1:sz,2) = 0.0d0
            ElementNodes % xyz(n+1:sz,3) = 0.0d0
        END IF

        sz1 = SIZE(Mesh % Nodes % x)
        IF (sz1 > Mesh % NumberOfNodes) THEN
            Indexes => GetIndexStore()
            nd = GetElementDOFs(Indexes,Element,NotDG=.TRUE.)
!DIR$ IVDEP
            DO i=n+1,nd
                IF ( Indexes(i)>0 .AND. Indexes(i)<=sz1 ) THEN
                    ElementNodes % x(i) = Mesh % Nodes % x(Indexes(i))
                    ElementNodes % y(i) = Mesh % Nodes % y(Indexes(i))
                    ElementNodes % z(i) = Mesh % Nodes % z(Indexes(i))
                END IF
            END DO
        END IF
    END SUBROUTINE GetElementNodesVec


    SUBROUTINE GetElementNodesOrigVec( ElementNodes, UElement, USolver, UMesh )
      TYPE(Nodes_t), TARGET :: ElementNodes
      TYPE(Element_t), OPTIONAL, TARGET :: UElement
      TYPE(Solver_t), OPTIONAL, TARGET :: USolver
      TYPE(Mesh_t), OPTIONAL, TARGET :: UMesh
      
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Nodes_t), POINTER :: TmpNodes
      
      IF( PRESENT( UMesh ) ) THEN
        Mesh => UMesh
      ELSE IF( PRESENT( USolver ) ) THEN
        Mesh => USolver % Mesh
      ELSE
        Mesh => CurrentModel % Solver % Mesh
      END IF

      TmpNodes => Mesh % Nodes
      IF(.NOT. ASSOCIATED( Mesh % NodesOrig ) ) THEN
        CALL Fatal('GetElementNodesOrigVec','Original node coordinates not yet stored!')
      END IF
      Mesh % Nodes => Mesh % NodesOrig
      
      CALL GetElementNodesVec( ElementNodes, UElement, UMesh = Mesh ) 
            
      Mesh % Nodes => TmpNodes
      
    END SUBROUTINE GetElementNodesOrigVec
      
        
    
!> Get element body id
!------------------------------------------------------------------------------
  FUNCTION GetBody( Element ) RESULT(body_id)
!------------------------------------------------------------------------------
   INTEGER::Body_id
   TYPE(Element_t), OPTIONAL :: Element
!------------------------------------------------------------------------------
   TYPE(Element_t), POINTER :: el
!------------------------------------------------------------------------------
   el => GetCurrentElement(Element)
   body_id= el % BodyId
!------------------------------------------------------------------------------
  END FUNCTION GetBody
!------------------------------------------------------------------------------


!> Get element body parameters
!------------------------------------------------------------------------------
  FUNCTION GetBodyParams(Element) RESULT(lst)
!------------------------------------------------------------------------------
   TYPE(ValueList_t), POINTER :: Lst
   TYPE(Element_t), OPTIONAL :: Element
!------------------------------------------------------------------------------
   TYPE(Element_t), POINTER :: el
!------------------------------------------------------------------------------
   lst => CurrentModel % Bodies(GetBody(Element)) % Values
!------------------------------------------------------------------------------
  END FUNCTION GetBodyParams
!------------------------------------------------------------------------------

  
!> Get the body force index of the active element
!------------------------------------------------------------------------------
  FUNCTION GetBodyForceId(  Element, Found ) RESULT(bf_id)
!------------------------------------------------------------------------------
     LOGICAL, OPTIONAL :: Found
     TYPE(Element_t), OPTIONAL :: Element
     TYPE(Element_t), POINTER :: CurrElement

     INTEGER :: bf_id, body_id

     CurrElement => GetCurrentElement(Element)
     body_id = CurrElement % BodyId 

     IF ( PRESENT( Found ) ) THEN
	bf_id = ListGetInteger( CurrentModel % Bodies(body_id) % Values, &
           'Body Force', Found, minv=1,maxv=CurrentModel % NumberOfBodyForces )
     ELSE
        bf_id = ListGetInteger( CurrentModel % Bodies(body_id) % Values, &
            'Body Force', minv=1,maxv=CurrentModel % NumberOfBodyForces )
     END IF
!------------------------------------------------------------------------------
  END FUNCTION GetBodyForceId
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Returns the material index of the active element
  FUNCTION GetMaterialId( Element, Found ) RESULT(mat_id)
!------------------------------------------------------------------------------
     LOGICAL, OPTIONAL :: Found
     TYPE(Element_t), OPTIONAL :: Element
     TYPE(Element_t), POINTER :: CurrElement

     INTEGER :: mat_id, body_id

     CurrElement => GetCurrentElement(Element)
     body_id = CurrElement % BodyId 

     IF( body_id <= 0 ) THEN
       mat_id = 0
       IF( PRESENT( Found ) ) Found = .FALSE.
       RETURN
     END IF
     
     IF ( PRESENT( Found ) ) THEN
        mat_id = ListGetInteger( CurrentModel % Bodies(body_id) % Values, &
           'Material', Found, minv=1,maxv=CurrentModel % NumberOfMaterials )
     ELSE
        mat_id = ListGetInteger( CurrentModel % Bodies(body_id) % Values, &
           'Material', minv=1,maxv=CurrentModel % NumberOfMaterials )
     END IF
!------------------------------------------------------------------------------
  END FUNCTION GetMaterialId
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Get component list given component id
  FUNCTION GetComponent(i) RESULT(list)
!------------------------------------------------------------------------------
     INTEGER :: i
     TYPE(ValueList_t), POINTER :: list

     List => Null()
     IF(i>=0 .AND. i<=SIZE(CurrentModel % Components)) list=> &
             CurrentModel % Components(i) % Values
!------------------------------------------------------------------------------
  END FUNCTION GetComponent
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Returns the equation index of the active element
  FUNCTION GetEquationId( Element, Found ) RESULT(eq_id)
!------------------------------------------------------------------------------
     LOGICAL, OPTIONAL :: Found
     TYPE(Element_t), OPTIONAL :: Element
     TYPE(Element_t), POINTER :: CurrElement

     INTEGER :: eq_id, body_id

     CurrElement => GetCurrentElement(Element)
     body_id = CurrElement % BodyId 

     IF ( PRESENT( Found ) ) THEN
        eq_id = ListGetInteger( CurrentModel % Bodies(body_id) % Values, &
           'Equation', Found, minv=1,maxv=CurrentModel % NumberOfEquations )
     ELSE
        eq_id = ListGetInteger( CurrentModel % Bodies(body_id) % Values, &
           'Equation',  minv=1,maxv=CurrentModel % NumberOfEquations )
     END IF
!------------------------------------------------------------------------------
  END FUNCTION GetEquationId
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Returns handle to the Simulation value list
  FUNCTION GetSimulation() RESULT(Simulation)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: Simulation
     Simulation => CurrentModel % Simulation
!------------------------------------------------------------------------------
  END FUNCTION GetSimulation
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Returns handle to the Constants value list
  FUNCTION GetConstants() RESULT(Constants)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: Constants
     Constants => CurrentModel % Constants
!------------------------------------------------------------------------------
  END FUNCTION GetConstants
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Returns handle to the Solver value list of the active solver
  FUNCTION GetSolverParams(Solver) RESULT(SolverParam)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: SolverParam
     TYPE(Solver_t), OPTIONAL :: Solver

     SolverParam => ListGetSolverParams(Solver)
!------------------------------------------------------------------------------
  END FUNCTION GetSolverParams
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Returns handle to Material value list of the active element
  FUNCTION GetMaterial(  Element, Found ) RESULT(Material)
!------------------------------------------------------------------------------
    TYPE(Element_t), OPTIONAL :: Element
    LOGICAL, OPTIONAL :: Found

    TYPE(ValueList_t), POINTER :: Material

    LOGICAL :: L
    INTEGER :: mat_id

    IF ( PRESENT( Element ) ) THEN
        mat_id = GetMaterialId( Element, L )
    ELSE
        mat_id = GetMaterialId( Found=L )
    END IF

    Material => Null()
    IF ( L ) Material => CurrentModel % Materials(mat_id) % Values
    IF ( PRESENT( Found ) ) Found = L
!------------------------------------------------------------------------------
  END FUNCTION GetMaterial
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Returns handle to Parent element of a boundary element with a larger body id.
!------------------------------------------------------------------------------ 
  FUNCTION GetBulkElementAtBoundary( Element, Found ) RESULT(BulkElement)
!------------------------------------------------------------------------------
    TYPE(Element_t), OPTIONAL :: Element
    LOGICAL, OPTIONAL :: Found
    TYPE(element_t), POINTER :: BulkElement
!------------------------------------------------------------------------------    
    TYPE(element_t), POINTER :: BulkElementL, BulkElementR, BoundaryElement
    LOGICAL :: L
    INTEGER :: mat_id, BodyIdL, BodyIdR

    BulkElement => NULL()
    
    BoundaryElement => GetCurrentElement(Element)
      
    IF ( .NOT. ASSOCIATED(BoundaryElement % boundaryinfo)) RETURN
    BulkElementR => BoundaryElement % boundaryinfo % right
    BulkElementL => BoundaryElement % boundaryinfo % left
    BodyIdR = 0; BodyIdL = 0
    
    IF (ASSOCIATED(BulkElementR)) BodyIdR = BulkElementR % BodyId
    IF (ASSOCIATED(BulkElementL)) BodyIdL = BulkElementL % BodyId
    
    IF (BodyIdR == 0 .AND. BodyIdL == 0) THEN
      RETURN
    ELSE IF (BodyIdR > BodyIdL) THEN
      BulkElement => BulkElementR
    ELSE IF (bodyIdL >= BodyIdR) THEN
      BulkElement => BulkElementL
    END IF

    IF( PRESENT( Found ) ) Found = ASSOCIATED( BulkElement ) 
    
!------------------------------------------------------------------------------
  END FUNCTION GetBulkElementAtBoundary
!------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
!> Returns handle to Material value list of the bulk material meeting  
!> element with larger body id. Typically Element is a boundary element.
  FUNCTION GetBulkMaterialAtBoundary( Element, Found ) RESULT(Material)
!------------------------------------------------------------------------------
    TYPE(Element_t), OPTIONAL :: Element
    LOGICAL, OPTIONAL :: Found
    TYPE(ValueList_t), POINTER :: Material
!------------------------------------------------------------------------------
    TYPE(element_t), POINTER :: BulkElement
    LOGICAL :: L
    INTEGER :: mat_id

    Material => NULL()

    BulkElement => GetBulkElementAtBoundary(Element, Found)

    IF( ASSOCIATED( BulkElement ) ) THEN
      mat_id = GetMaterialId( BulkElement, L )      
      IF ( L ) Material => CurrentModel % Materials(mat_id) % Values
    ELSE
      L = .FALSE.
    END IF
    
    IF ( PRESENT( Found ) ) Found = L
!------------------------------------------------------------------------------
  END FUNCTION GetBulkMaterialAtBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Return handle to the Body Force value list of the active element
  FUNCTION GetBodyForce( Element, Found ) RESULT(BodyForce)
!------------------------------------------------------------------------------
    TYPE(Element_t), OPTIONAL :: Element
    LOGICAL, OPTIONAL :: Found

    TYPE(ValueList_t), POINTER :: BodyForce

    LOGICAL :: l
    INTEGER :: bf_id

    IF ( PRESENT( Element ) ) THEN
       bf_id = GetBodyForceId( Element, L )
    ELSE
       bf_id = GetBodyForceId( Found=L )
    END IF

    BodyForce => Null()
    IF ( L ) BodyForce => CurrentModel % BodyForces(bf_id) % Values
    IF ( PRESENT( Found ) ) Found = L
!------------------------------------------------------------------------------
  END FUNCTION GetBodyForce
!------------------------------------------------------------------------------


!> Is the active solver solved in the frequency space
!------------------------------------------------------------------------------
  FUNCTION EigenOrHarmonicAnalysis(Usolver) RESULT(L)
    LOGICAL :: L
    TYPE(Solver_t), OPTIONAL,TARGET :: USolver
!------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: Solver

    IF (PRESENT(USolver)) THEN
      Solver => USolver
    ELSE
      Solver => CurrentModel % Solver
    END IF
    L  = Solver % NOFEigenValues > 0
!------------------------------------------------------------------------------
  END FUNCTION EigenOrHarmonicAnalysis
!------------------------------------------------------------------------------


!> Returns the handle to the equation where the active element belongs to 
!------------------------------------------------------------------------------
  FUNCTION GetEquation( Element, Found ) RESULT(Equation)
!------------------------------------------------------------------------------
    TYPE(Element_t), OPTIONAL :: Element
    LOGICAL, OPTIONAL :: Found

    TYPE(ValueList_t), POINTER :: Equation

    LOGICAL :: L
    INTEGER :: eq_id


    IF ( PRESENT( Element ) ) THEN
       eq_id = GetEquationId( Element, L )
    ELSE
       eq_id = GetEquationId( Found=L )
    END IF

    NULLIFY( Equation )
    IF ( L ) Equation => CurrentModel % Equations(eq_id) % Values
    IF ( PRESENT( Found ) ) Found = L
!------------------------------------------------------------------------------
  END FUNCTION GetEquation
!------------------------------------------------------------------------------



!> Returns the Boundary Condition index of the active element
!------------------------------------------------------------------------------
  FUNCTION GetBCId( UElement ) RESULT(bc_id)
!------------------------------------------------------------------------------
     TYPE(Element_t), OPTIONAL, TARGET :: UElement

     INTEGER :: bc_id

     TYPE(Element_t), POINTER :: Element

     Element => GetCurrentElement( UElement )
     
     IF(.NOT. ASSOCIATED( Element % BoundaryInfo ) ) THEN
       bc_id = 0
       RETURN
     END IF
     
     DO bc_id=1,CurrentModel % NumberOfBCs
        IF ( Element % BoundaryInfo % Constraint == CurrentModel % BCs(bc_id) % Tag ) EXIT
     END DO
     IF ( bc_id > CurrentModel % NumberOfBCs ) bc_id=0
!------------------------------------------------------------------------------
  END FUNCTION GetBCId
!------------------------------------------------------------------------------


!> Returns handle to the value list of the Boundary Condition where the active element belongs to
!------------------------------------------------------------------------------
  FUNCTION GetBC( UElement ) RESULT(bc)
!------------------------------------------------------------------------------
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     TYPE(ValueList_t), POINTER :: BC

     INTEGER :: bc_id

     TYPE(Element_t), POINTER :: Element

     Element => GetCurrentElement( UElement )
     
     BC => Null()
     bc_id = GetBCId( Element )
     IF ( bc_id > 0 )  BC => CurrentModel % BCs(bc_id) % Values
!------------------------------------------------------------------------------
  END FUNCTION GetBC
!------------------------------------------------------------------------------


!> Returns the index of the Initial Condition of the active element
!------------------------------------------------------------------------------
  FUNCTION GetICId( Element, Found ) RESULT(ic_id)
!------------------------------------------------------------------------------
     LOGICAL, OPTIONAL :: Found
     TYPE(Element_t), OPTIONAL :: Element

     TYPE(Element_t), POINTER :: CElement
     INTEGER :: ic_id, body_id

     CElement => GetCurrentElement( Element )
     body_id = CElement % BodyId

     IF ( PRESENT( Found ) ) THEN
        ic_id = ListGetInteger( CurrentModel % Bodies(body_id) % Values, &
           'Initial Condition', Found, minv=1,maxv=CurrentModel % NumberOfICs )
     ELSE
        ic_id = ListGetInteger( CurrentModel % Bodies(body_id) % Values, &
           'Initial Condition', minv=1,maxv=CurrentModel % NumberOfICs )
     END IF
!------------------------------------------------------------------------------
  END FUNCTION GetIcId
!------------------------------------------------------------------------------

!> Returns handle to the value list of the Initial Condition where the active element belongs to
!------------------------------------------------------------------------------
  FUNCTION GetIC(  Element, Found ) RESULT(IC)
!------------------------------------------------------------------------------
    TYPE(Element_t), OPTIONAL :: Element
    LOGICAL, OPTIONAL :: Found

    TYPE(ValueList_t), POINTER :: IC

    LOGICAL :: L
    INTEGER :: ic_id

    IF ( PRESENT( Element ) ) THEN
        ic_id = GetICId( Element, L )
    ELSE
        ic_id = GetICId( Found=L )
    END IF

    IC => Null()
    IF ( L ) IC => CurrentModel % ICs(ic_id) % Values
    IF ( PRESENT( Found ) ) Found = L
!------------------------------------------------------------------------------
  END FUNCTION GetIC
!------------------------------------------------------------------------------

!> Add the local matrix entries to for real valued equations that are of first order in time
!------------------------------------------------------------------------------
  SUBROUTINE Default1stOrderTimeR( M, A, F, UElement, USolver )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: M(:,:),A(:,:), F(:)
    TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
    TYPE(Element_t), OPTIONAL, TARGET :: UElement

    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: Params

    TYPE(Solver_t), POINTER :: Solver
    TYPE(Variable_t), POINTER :: x
    TYPE(Element_t), POINTER :: Element

    INTEGER :: n
    REAL(KIND=dp) :: dt
    INTEGER, POINTER :: Indexes(:)

    IF ( PRESENT(USolver) ) THEN
      Solver => USolver
    ELSE
      Solver => CurrentModel % Solver
    END IF

    Params => GetSolverParams(Solver)

    ! Antiperiodic elimination and FCT always use this
    IF (GetLogical(Params,'Use Global Mass Matrix',Found)) THEN
      CALL DefaultUpdateMass(M,UElement,USolver)
      RETURN
    END IF

    Element => GetCurrentElement( UElement )

    x => Solver % Variable

    dt = Solver % dt
    Indexes => GetIndexStore()
    n = GetElementDOFs( Indexes,Element,Solver )
          
    CALL Add1stOrderTime( M, A, F, dt, n, x % DOFs, &
        x % Perm(Indexes(1:n)), Solver, UElement=Element )
      
!------------------------------------------------------------------------------
  END SUBROUTINE Default1stOrderTimeR
!------------------------------------------------------------------------------

!> Add the local matrix entries to for complex valued equations that are of first order in time
!------------------------------------------------------------------------------
  SUBROUTINE Default1stOrderTimeC( MC, AC, FC, UElement, USolver )
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: MC(:,:),AC(:,:), FC(:)
    TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
    TYPE(Element_t), OPTIONAL, TARGET :: UElement

    TYPE(Solver_t), POINTER :: Solver
    TYPE(Variable_t), POINTER :: x
    TYPE(Element_t), POINTER :: Element

    REAL(KIND=dp), ALLOCATABLE :: M(:,:),A(:,:), F(:)

    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: Params

    INTEGER :: i,j,n,DOFs
    REAL(KIND=dp) :: dt
    INTEGER, POINTER :: Indexes(:)

    IF ( PRESENT(USolver) ) THEN
      Solver => USolver
    ELSE
      Solver => CurrentModel % Solver
    END IF

    Params=>GetSolverParams(Solver)

    IF (GetLogical(Params,'Use Global Mass Matrix',Found)) THEN
      CALL DefaultUpdateMass(M,UElement,USolver)
      RETURN
    END IF

    Element => GetCurrentElement( UElement ) 

    x => Solver % Variable

    dt = Solver % dt
    DOFs = x % DOFs
    Indexes => GetIndexStore()
    n = GetElementDOFs( Indexes,Element,Solver )

    ALLOCATE( M(DOFs*n,DOFs*n), A(DOFs*n,DOFs*n), F(DOFs*n) )
    DO i=1,n*DOFs/2
      F( 2*(i-1)+1 ) =  REAL( FC(i) )
      F( 2*(i-1)+2 ) = AIMAG( FC(i) )

      DO j=1,n*DOFs/2
        M( 2*(i-1)+1, 2*(j-1)+1 ) =   REAL( MC(i,j) )
        M( 2*(i-1)+1, 2*(j-1)+2 ) = -AIMAG( MC(i,j) )
        M( 2*(i-1)+2, 2*(j-1)+1 ) =  AIMAG( MC(i,j) )
        M( 2*(i-1)+2, 2*(j-1)+2 ) =   REAL( MC(i,j) )
        A( 2*(i-1)+1, 2*(j-1)+1 ) =   REAL( AC(i,j) )
        A( 2*(i-1)+1, 2*(j-1)+2 ) = -AIMAG( AC(i,j) )
        A( 2*(i-1)+2, 2*(j-1)+1 ) =  AIMAG( AC(i,j) )
        A( 2*(i-1)+2, 2*(j-1)+2 ) =   REAL( AC(i,j) )
      END DO
    END DO

    CALL Add1stOrderTime( M, A, F, dt, n, x % DOFs, &
           x % Perm(Indexes(1:n)), Solver, UElement=Element )

    DO i=1,n*DOFs/2
      FC(i) = CMPLX( F(2*(i-1)+1), F(2*(i-1)+2),KIND=dp )
      DO j=1,n*DOFs/2
        MC(i,j) = CMPLX(M(2*(i-1)+1,2*(j-1)+1), -M(2*(i-1)+1,2*(j-1)+2), KIND=dp)
        AC(i,j) = CMPLX(A(2*(i-1)+1,2*(j-1)+1), -A(2*(i-1)+1,2*(j-1)+2), KIND=dp)
      END DO
    END DO

    DEALLOCATE( M, A, F )
!------------------------------------------------------------------------------
  END SUBROUTINE Default1stOrderTimeC
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE Default1stOrderTimeGlobal(USolver)
!------------------------------------------------------------------------------
   TYPE(Solver_t), OPTIONAL, TARGET :: USolver
!------------------------------------------------------------------------------
   CHARACTER(:), ALLOCATABLE :: Method
   TYPE(Solver_t), POINTER :: Solver
   INTEGER :: i,j,k,l,n,Order
   REAL(KIND=dp), POINTER :: SaveValues(:) => NULL()
   REAL(KIND=dp) :: FORCE(1), Dts(16)
   LOGICAL :: ConstantDt, Found, HasMass, HasFCT
   TYPE(Variable_t), POINTER :: DtVar
   SAVE STIFF, MASS, X
   REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:),MASS(:,:), X(:,:)

   !$OMP THREADPRIVATE(SaveValues)

   IF ( PRESENT(USolver) ) THEN
     Solver => Usolver
   ELSE
     Solver => CurrentModel % solver
   END IF

   Order = MAX( MIN( Solver % DoneTime, Solver % Order ), 1)   
   HasMass = ASSOCIATED( Solver % Matrix % MassValues )

   HasFCT = ListGetLogical( Solver % Values,'Linear System FCT', Found )

   IF( HasFCT ) THEN
     IF( .NOT. HasMass ) THEN
       CALL Fatal('Default1stOrderTimeGlobal','FCT only makes sense if there is a mass matrix!')
     ELSE
       IF(.NOT. ASSOCIATED( Solver % Matrix % MassValuesLumped ) ) THEN
         CALL Fatal('Default1stOrderTimeGlobal','FCT requires a lumped mass matrix!')
       END IF
       HasMass = .FALSE.
     END IF
   END IF

   ! This is now the default global time integration routine but the old hack may still be called
   !---------------------------------------------------------------------------------------------
   IF( .NOT. ListGetLogical( Solver % Values,'Old Global Time Integration',Found ) ) THEN
     CALL Add1stOrderTime_CRS( Solver % Matrix, Solver % Matrix % rhs, &
         Solver % dt, Solver )
     RETURN
   END IF


   ! The rest of the code in this subroutine is obsolete
   IF ( .NOT.ASSOCIATED(Solver % Variable % Values, SaveValues) ) THEN
     IF ( ALLOCATED(STIFF) ) DEALLOCATE( STIFF,MASS,X )
     n = 0
     DO i=1,Solver % Matrix % NumberOfRows
       n = MAX( n,Solver % Matrix % Rows(i+1)-Solver % Matrix % Rows(i) )
     END DO
     k = SIZE(Solver % Variable % PrevValues,2)
     ALLOCATE( STIFF(1,n),MASS(1,n),X(n,k) )     
     SaveValues => Solver % Variable % Values
   END IF
   
   STIFF = 0.0_dp
   MASS = 0.0_dp
   X = 0.0_dp

   Method = GetString( Solver % Values, 'Timestepping Method', Found )
   IF ( Method == 'bdf' ) THEN
     Dts(1) = Solver % Dt
     ConstantDt = .TRUE.
     IF(Order > 1) THEN
       DtVar => VariableGet( Solver % Mesh % Variables, 'Timestep size' )
       DO i=2,Order
         Dts(i) = DtVar % PrevValues(1,i-1)
         IF(ABS(Dts(i)-Dts(1)) > 1.0d-6 * Dts(1)) ConstantDt = .FALSE.
       END DO
     END IF
   END IF
   
   DO i=1,Solver % Matrix % NumberOFRows
     n = 0
     k = 0

     DO j=Solver % Matrix % Rows(i),Solver % Matrix % Rows(i+1)-1
       n = n+1
       STIFF(1,n) = Solver % Matrix % Values(j)
       IF( HasMass ) THEN
         MASS(1,n) = Solver % Matrix % MassValues(j)
       ELSE IF( HasFCT ) THEN
         IF( j == Solver % Matrix % Diag(i) ) k = n
       END IF         
       X(n,:) = Solver % Variable % PrevValues(Solver % Matrix % Cols(j),:)
     END DO

     ! Use lumped mass in lower order fct
     IF( HasFCT ) THEN
       IF( k == 0 ) THEN
         CALL Fatal('Default1stOrderTimeGlobal','Could not find diagonal entry for fct')
       ELSE
         MASS(1,k) = Solver % Matrix % MassValuesLumped(i)
       END IF
     END IF

     FORCE(1) = Solver % Matrix % RHS(i)
     Solver % Matrix % Force(i,1) = FORCE(1)

     SELECT CASE( Method )
     CASE( 'fs' )
       CALL FractionalStep( n, Solver % dt, MASS, STIFF, FORCE, &
           X(:,1), Solver % Beta, Solver )
       
     CASE('bdf')       
       IF(ConstantDt) THEN
         CALL BDFLocal( n, Solver % dt, MASS, STIFF, FORCE, X, Order )
       ELSE
         CALL VBDFLocal(n, Dts, MASS, STIFF, FORCE, X, Order )
       END IF
       
     CASE DEFAULT
       CALL NewmarkBeta( n, Solver % dt, MASS, STIFF, FORCE, &
           X(:,1), Solver % Beta )
     END SELECT

     IF( HasFCT ) MASS(1,k) = 0.0_dp

     n = 0
     DO j=Solver % Matrix % Rows(i),Solver % Matrix % Rows(i+1)-1
       n=n+1
       Solver % Matrix % Values(j) = STIFF(1,n)
     END DO
     Solver % Matrix % RHS(i) = FORCE(1)
   END DO

!----------------------------------------------------------------------------
  END SUBROUTINE Default1stOrderTimeGlobal
!----------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE Default2ndOrderTimeGlobal(USolver)
!------------------------------------------------------------------------------
   TYPE(Solver_t), OPTIONAL, TARGET :: USolver
!------------------------------------------------------------------------------
   TYPE(Solver_t), POINTER :: Solver
   INTEGER :: i,j,k,l,n
   REAL(KIND=dp), POINTER :: SaveValues(:) => NULL()
   REAL(KIND=dp) :: FORCE(1)
   LOGICAL :: Found, HasDamping, HasMass
   REAL(KIND=dp), ALLOCATABLE, SAVE :: STIFF(:,:),MASS(:,:), DAMP(:,:), X(:,:)
   !OMP THREADPRIVATE(SaveValues)

   IF ( PRESENT(USolver) ) THEN
     Solver => Usolver
   ELSE
     Solver => CurrentModel % solver
   END IF
   
   ! This is now the default global time integration routine but the old hack may still be called
   !---------------------------------------------------------------------------------------------
   IF( .NOT. ListGetLogical( Solver % Values,'Old Global Time Integration',Found ) ) THEN
     CALL Add2ndOrderTime_CRS( Solver % Matrix, Solver % Matrix % rhs, &
         Solver % dt, Solver % Variable % PrevValues, Solver )
     RETURN
   END IF

   
   IF ( .NOT.ASSOCIATED(Solver % Variable % Values, SaveValues) ) THEN
      IF ( ALLOCATED(STIFF) ) DEALLOCATE( STIFF,MASS,DAMP,X )
      n = 0
      DO i=1,Solver % Matrix % NumberOfRows
        n = MAX( n,Solver % Matrix % Rows(i+1)-Solver % Matrix % Rows(i) )
      END DO
      k = SIZE(Solver % Variable % PrevValues,2)
      ALLOCATE( STIFF(1,n),MASS(1,n),DAMP(1,n),X(n,k) )
      SaveValues => Solver % Variable % Values

      STIFF = 0.0_dp
      MASS = 0.0_dp
      DAMP = 0.0_dp
      X = 0.0_dp
    END IF

    HasDamping = ASSOCIATED(Solver % Matrix % DampValues )
    HasMass = ASSOCIATED(Solver % Matrix % MassValues )

    DO i=1,Solver % Matrix % NumberOFRows
      n = 0
      DO j=Solver % Matrix % Rows(i),Solver % Matrix % Rows(i+1)-1
        n=n+1
        IF( HasMass ) MASS(1,n) = Solver % Matrix % MassValues(j)
        IF( HasDamping ) DAMP(1,n) = Solver % Matrix % DampValues(j)
        STIFF(1,n) = Solver % Matrix % Values(j)
        X(n,:) = Solver % Variable % PrevValues(Solver % Matrix % Cols(j),:)
      END DO
      FORCE(1) = Solver % Matrix % RHS(i)
      Solver % Matrix % Force(i,1) = FORCE(1)
      
      CALL Time2ndOrder( n, Solver % dt, MASS, DAMP, STIFF, &
          FORCE, X(1:n,3), X(1:n,4), X(1:n,5), X(1:n,7), Solver % Alpha, Solver % Beta )
      
      n = 0
      DO j=Solver % Matrix % Rows(i),Solver % Matrix % Rows(i+1)-1
        n=n+1
        Solver % Matrix % Values(j) = STIFF(1,n)
      END DO
      Solver % Matrix % RHS(i) = FORCE(1)
    END DO
!----------------------------------------------------------------------------
  END SUBROUTINE Default2ndOrderTimeGlobal
!----------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE Default2ndOrderTimeR( M, B, A, F, UElement, USolver )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: M(:,:), B(:,:), A(:,:), F(:)
    TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
    TYPE(Element_t), OPTIONAL, TARGET :: UElement

    TYPE(Solver_t), POINTER :: Solver
    TYPE(Variable_t), POINTER :: x
    TYPE(Element_t), POINTER :: Element

    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: Params

    INTEGER :: n
    REAL(KIND=dp) :: dt
    INTEGER, POINTER :: Indexes(:)

    Solver => CurrentModel % Solver
    IF ( PRESENT(USolver) ) Solver => USolver

    Params=>GetSolverParams(Solver)

    IF (GetLogical(Params,'Use Global Mass Matrix',Found)) THEN
      CALL DefaultUpdateMass(M,UElement,USolver)
      CALL DefaultUpdateDamp(B,UElement,USolver)
      RETURN
    END IF

    Element => GetCurrentElement( UElement ) 

    x => Solver % Variable

    dt = Solver % dt
    Indexes => GetIndexStore()
    n = GetElementDOFs( Indexes, Element, Solver )

    CALL Add2ndOrderTime( M, B, A, F, dt, n, x % DOFs, &
          x % Perm(Indexes(1:n)), Solver )
!------------------------------------------------------------------------------
  END SUBROUTINE Default2ndOrderTimeR
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE Default2ndOrderTimeC( MC, BC, AC, FC, UElement, USolver )
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: MC(:,:), BC(:,:), AC(:,:), FC(:)
    TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
    TYPE(Element_t), OPTIONAL, TARGET :: UElement

    TYPE(Solver_t), POINTER :: Solver
    TYPE(Variable_t), POINTER :: x
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), ALLOCATABLE :: M(:,:), B(:,:), A(:,:), F(:)

    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: Params

    INTEGER :: i,j,n,DOFs
    REAL(KIND=dp) :: dt
    INTEGER, POINTER :: Indexes(:)

    Solver => CurrentModel % Solver
    IF ( PRESENT(USolver) ) Solver => USolver

    Params=>GetSolverParams(Solver)

    IF (GetLogical(Params,'Use Global Mass Matrix',Found)) THEN
      CALL DefaultUpdateMass(M,UElement,USolver)
      CALL DefaultUpdateDamp(B,UElement,USolver)
      RETURN
    END IF

    Element => GetCurrentElement( UElement ) 
    
    x => Solver % Variable

    dt = Solver % dt
    DOFs = x % DOFs
    Indexes => GetIndexStore()
    n = GetElementDOFs( Indexes, Element, Solver )

    ALLOCATE( M(DOFs*n,DOFs*n), A(DOFs*n,DOFs*n), B(DOFs*n,DOFs*n), F(DOFs*n) )
    DO i=1,n*DOFs/2
      F( 2*(i-1)+1 ) =  REAL( FC(i) )
      F( 2*(i-1)+2 ) = AIMAG( FC(i) )

      DO j=1,n*DOFs/2
        M(2*(i-1)+1, 2*(j-1)+1) =   REAL( MC(i,j) )
        M(2*(i-1)+1, 2*(j-1)+2) = -AIMAG( MC(i,j) )
        M(2*(i-1)+2, 2*(j-1)+1) =  AIMAG( MC(i,j) )
        M(2*(i-1)+2, 2*(j-1)+2) =   REAL( MC(i,j) )
        B(2*(i-1)+1, 2*(j-1)+1) =   REAL( BC(i,j) )
        B(2*(i-1)+1, 2*(j-1)+2) = -AIMAG( BC(i,j) )
        B(2*(i-1)+2, 2*(j-1)+1) =  AIMAG( BC(i,j) )
        B(2*(i-1)+2, 2*(j-1)+2) =   REAL( BC(i,j) )
        A(2*(i-1)+1, 2*(j-1)+1) =   REAL( AC(i,j) )
        A(2*(i-1)+1, 2*(j-1)+2) = -AIMAG( AC(i,j) )
        A(2*(i-1)+2, 2*(j-1)+1) =  AIMAG( AC(i,j) )
        A(2*(i-1)+2, 2*(j-1)+2) =   REAL( AC(i,j) )
      END DO
    END DO

    CALL Add2ndOrderTime( M, B, A, F, dt, n, x % DOFs, &
          x % Perm(Indexes(1:n)), Solver )

    DO i=1,n*DOFs/2
      FC(i) = CMPLX( F(2*(i-1)+1), F(2*(i-1)+2), KIND=dp )
      DO j=1,n*DOFs/2
        MC(i,j) = CMPLX( M(2*(i-1)+1, 2*(j-1)+1), -M(2*(i-1)+1, 2*(j-1)+2), KIND=dp )
        BC(i,j) = CMPLX( B(2*(i-1)+1, 2*(j-1)+1), -B(2*(i-1)+1, 2*(j-1)+2), KIND=dp )
        AC(i,j) = CMPLX( A(2*(i-1)+1, 2*(j-1)+1), -A(2*(i-1)+1, 2*(j-1)+2), KIND=dp )
      END DO
    END DO

    DEALLOCATE( M, B, A, F )
!------------------------------------------------------------------------------
  END SUBROUTINE Default2ndOrderTimeC
!------------------------------------------------------------------------------


!--------------------------------------------------------------------------------
!> One can enforce weak coupling by calling a dependent solver a.k.a. slave solver
!> at different stages of the master solver: e.g. before and after the solver.
!> The strategy can be particularly efficient for nonlinear problems when the
!> slave solver is cheap and a stepsize control is applied
!> Also one can easily make postprocessing steps just at the correct slot.
!-----------------------------------------------------------------------------
  RECURSIVE SUBROUTINE DefaultSlaveSolvers( Solver, SlaveSolverStr)
!------------------------------------------------------------------------------  
     TYPE(Solver_t), POINTER :: Solver     
     CHARACTER(LEN=*) :: SlaveSolverStr 
     
     TYPE(Solver_t), POINTER :: SlaveSolver
     TYPE(ValueList_t), POINTER :: Params
     TYPE(Variable_t), POINTER :: iterV
     INTEGER, POINTER :: SlaveSolverIndexes(:)
     INTEGER :: j,k,iter
     REAL(KIND=dp) :: dt
     LOGICAL :: Transient, Found, alloc_parenv

     TYPE(ParEnv_t), POINTER :: SParEnv

     INTERFACE
       SUBROUTINE SolverActivate_x(Model,Solver,dt,Transient)
         USE Types
         TYPE(Model_t)::Model
         TYPE(Solver_t),POINTER::Solver
         REAL(KIND=dp) :: dt
         LOGICAL :: Transient
       END SUBROUTINE SolverActivate_x
     END INTERFACE

     SlaveSolverIndexes =>  ListGetIntegerArray( Solver % Values,&
         SlaveSolverStr,Found )
     IF(.NOT. Found ) RETURN

     CALL Info('DefaultSlaveSolvers','Executing slave solvers: '// &
         TRIM(SlaveSolverStr),Level=6)
     
     dt = GetTimeStepsize()
     Transient = GetString(CurrentModel % Simulation,'Simulation type',Found)=='transient'

     ! store the nonlinear iteration at the outer loop
     iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
     iter = NINT(iterV % Values(1))

     
     DO j=1,SIZE(SlaveSolverIndexes)
       k = SlaveSolverIndexes(j)
       SlaveSolver => CurrentModel % Solvers(k)

       CALL Info('DefaultSlaveSolvers','Calling slave solver: '//I2S(k),Level=8)

       IF( ListGetLogical( Solver % Values,'Monolithic Slave',Found )  ) THEN
         IF(.NOT. ListCheckPresent( SlaveSolver % Values,'Linear System Solver Disabled') ) THEN
           CALL Info('DefaultSlaveSolvers','Disabling linear system solver for slave: '//I2S(k),Level=6)
           CALL ListAddLogical(SlaveSolver % Values,'Linear System Solver Disabled',.TRUE.)
         END IF
       END IF
         
       IF(ParEnv % PEs>1) THEN
         SParEnv => ParEnv

         IF(ASSOCIATED(SlaveSolver % Matrix)) THEN
           IF(ASSOCIATED(SlaveSolver % Matrix % ParMatrix) ) THEN
             ParEnv => SlaveSolver % Matrix % ParMatrix % ParEnv
           ELSE
             ParEnv % ActiveComm = SlaveSolver % Matrix % Comm
           END IF
         ELSE
           CALL ListAddLogical( SlaveSolver % Values, 'Slave not parallel', .TRUE.)
         END IF
       END IF

       CurrentModel % Solver => SlaveSolver
       CALL SolverActivate_x( CurrentModel,SlaveSolver,dt,Transient)

       IF(ParEnv % PEs>1) THEN
         ParEnv => SParEnv
       END IF
     END DO
     iterV % Values = iter       
     CurrentModel % Solver => Solver

   END SUBROUTINE DefaultSlaveSolvers
!------------------------------------------------------------------------------
 
  

!> Performs initialization for matrix equation related to the active solver
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE DefaultInitialize( USolver, UseConstantBulk )
!------------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL, TARGET, INTENT(IN) :: USolver
     LOGICAL, OPTIONAL :: UseConstantBulk
!------------------------------------------------------------------------------
     TYPE(Solver_t), POINTER :: Solver
     INTEGER :: i,n
     LOGICAL :: Found
     
     IF ( PRESENT( USolver ) ) THEN
       Solver => USolver
     ELSE
       Solver => CurrentModel % Solver
     END IF

     IF(.NOT. ASSOCIATED( Solver % Matrix ) ) THEN
       CALL Fatal('DefaultInitialize','No matrix exists, cannot initialize!')
     END IF     

     IF( PRESENT( UseConstantBulk ) ) THEN
       IF ( UseConstantBulk ) THEN
         IF (.NOT. ASSOCIATED( Solver % Matrix % BulkRhs ) ) THEN
           Solver % Matrix % rhs = 0.0d0
         END IF

         CALL Info('DefaultInitialize','Using constant bulk matrix',Level=8)
         IF (.NOT. ASSOCIATED( Solver % Matrix % BulkValues ) ) THEN
           CALL Warn('DefaultInitialize','Constant bulk system requested but not associated!')
           RETURN
         END IF

         CALL RestoreBulkMatrix(Solver % Matrix)
         RETURN
       END IF
     END IF

     IF( ListGetLogical( Solver % Values,'Apply Explicit Control', Found )) THEN
       CALL ApplyExplicitControl( Solver )
     END IF

     
     CALL DefaultSlaveSolvers(Solver,'Slave Solvers') ! this is the initial name of the slot
     CALL DefaultSlaveSolvers(Solver,'Nonlinear Pre Solvers')     


     ! If we changed the system last time to harmonic one then revert back the real system
     IF( ListGetLogical( Solver % Values,'Harmonic Mode',Found ) ) THEN
       CALL ChangeToHarmonicSystem( Solver, .TRUE. )
     END IF
     
     CALL InitializeToZero( Solver % Matrix, Solver % Matrix % RHS )

     IF( ALLOCATED(Solver % Matrix % ConstrainedDOF) ) THEN
       Solver % Matrix % ConstrainedDOF = .FALSE.
     END IF
       
     IF( ALLOCATED(Solver % Matrix % Dvalues) ) THEN
       Solver % Matrix % Dvalues = 0._dp
     END IF

     IF( ListGetLogical( Solver % Values,'Bulk Assembly Timing',Found ) ) THEN 
       CALL ResetTimer('BulkAssembly'//GetVarName(Solver % Variable) ) 
     END IF

     ! This is a slot for calling solver that contribute to the assembly
     CALL DefaultSlaveSolvers(Solver,'Assembly Solvers')
                
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultInitialize
!------------------------------------------------------------------------------



!> Performs pre-steps related to the active solver
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE DefaultStart( USolver )
!------------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL, TARGET, INTENT(IN) :: USolver     
     TYPE(Solver_t), POINTER :: Solver
     LOGICAL :: Found
     TYPE(ValueList_t), POINTER :: Params
     INTEGER :: i,j,n
     
     IF ( PRESENT( USolver ) ) THEN
       Solver => USolver
     ELSE
       Solver => CurrentModel % Solver
     END IF

     Params => Solver % Values
     
     CALL Info('DefaultStart','Starting solver: '//&
        GetString(Params,'Equation'),Level=10)
          
     ! When Newton linearization is used we may reset it after previously visiting the solver
     IF( Solver % NewtonActive ) THEN
       IF( ListGetLogical( Params,'Nonlinear System Reset Newton', Found) ) Solver % NewtonActive = .FALSE.
     END IF
          
     ! If we changed the system last time to harmonic one then revert back the real system
     IF( ListGetLogical( Params,'Harmonic Mode',Found ) ) THEN
       CALL ChangeToHarmonicSystem( Solver, .TRUE. )
     END IF

     ! One can run preprocessing solver in this slot.
     !-----------------------------------------------------------------------------
     CALL DefaultSlaveSolvers(Solver,'Pre Solvers')

     IF( ListGetLogical(Params,'Local Matrix Storage',Found ) ) THEN
       IF(.NOT. ASSOCIATED(Solver % InvActiveElements) ) THEN
         ALLOCATE( Solver % InvActiveElements( Solver % Mesh % NumberOfBulkElements &
             + Solver % Mesh % NumberOFBoundaryElements ) )         
         Solver % InvActiveElements = 0
         DO i=1,Solver % NumberOfActiveElements
           Solver % InvActiveElements( Solver % ActiveElements(i) ) = i
         END DO
       END IF

       n = Solver % NumberOfActiveElements
       IF(ASSOCIATED(Solver % LocalSystem)) THEN
         IF(SIZE(Solver % LocalSystem) < n ) DEALLOCATE(Solver % LocalSystem)
       END IF
       IF(.NOT. ASSOCIATED(Solver % LocalSystem ) ) THEN
         ALLOCATE( Solver % LocalSystem(n) )
         Solver % LocalSystem(1:n) % eind = 0         
         ! If the stiffness matrix is constant the 1st element gives stiffness matrix for all!
         ! This could be inhereted differently too for splitted meshes, for example. 
         IF( ListGetLogical( Params,'Local Matrix Identical', Found )  ) THEN
           CALL Info('DefaultStart','Assuming all elements to be identical!')
           Solver % LocalSystem(1:n) % eind = 1         
         ELSE IF( ListGetLogical( Params,'Local Matrix Identical Bodies', Found )  ) THEN
           CALL Info('DefaultStart','Assuming all elements to be identical within bodies!')
           BLOCK
             INTEGER, ALLOCATABLE :: Body1st(:)             
             ALLOCATE(Body1st(CurrentModel % NumberOfBodies))
             Body1st = 0
             DO i=1,Solver % NumberOfActiveElements
               j = Solver % Mesh % Elements(Solver % ActiveElements(i)) % BodyId
               IF(Body1st(j) == 0) Body1st(j) = i
               Solver % LocalSystem(i) % eind = Body1st(j)
             END DO
           END BLOCK
         END IF
       END IF
       
       Solver % LocalSystemMode = 1
     END IF
     
!------------------------------------------------------------------------------
   END SUBROUTINE DefaultStart
!------------------------------------------------------------------------------


  
!> Performs finalizing steps related to the active solver
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE DefaultFinish( USolver )
!------------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL, TARGET, INTENT(IN) :: USolver
     TYPE(Solver_t), POINTER :: Solver
     TYPE(ValueList_t), POINTER :: Params
     CHARACTER(:), ALLOCATABLE :: str
     LOGICAL :: Found

     IF ( PRESENT( USolver ) ) THEN
       Solver => USolver
     ELSE
       Solver => CurrentModel % Solver
     END IF

     Params => Solver % Values

     IF ( ListGetLogical( Params,'Linear System Save',Found )) THEN
       str = GetString( Params,'Linear System Save Slot', Found )
       IF(Found .AND. str == 'finish') THEN
         CALL SaveLinearSystem( Solver ) 
       END IF
     END IF
             
     ! One can run postprocessing solver in this slot.
     !-----------------------------------------------------------------------------
     CALL DefaultSlaveSolvers(Solver,'Post Solvers')

     IF( ListGetLogical( Params,'Apply Explicit Control', Found )) THEN
       CALL ApplyExplicitControl( Solver )
     END IF

     IF( Solver % NumberOfConstraintModes > 0 ) THEN
       ! If we have a frozen stat then the nonlinear system loop is used to find that frozen state
       ! and we perform the linearized constraint modes analysis at the end. 
       IF( ListGetLogical(Params,'Constraint Modes Analysis Frozen',Found ) ) THEN
         BLOCK 
           INTEGER :: n
           REAL(KIND=dp) :: Norm
           REAL(KIND=dp), ALLOCATABLE :: xtmp(:), btmp(:)
           REAL(KIND=dp), POINTER :: rhs(:)

           CALL ListAddLogical(Params,'Constraint Modes Analysis Frozen',.FALSE.)
           n = SIZE(Solver % Matrix % rhs)
           rhs => Solver % Matrix % rhs           
           ALLOCATE(xtmp(n),btmp(n))
           xtmp = 0.0_dp; btmp = 0.0_dp
           CALL SolveSystem( Solver % Matrix, ParMatrix, btmp, xtmp, Norm,Solver % Variable % DOFs,Solver )
           CALL ListAddLogical(Params,'Constraint Modes Analysis Frozen',.TRUE.)
         END BLOCK
       END IF
         
       IF( ListGetLogical( Params,'Nonlinear System Constraint Modes', Found ) ) THEN
         CALL FinalizeLumpedMatrix( Solver )            
       END IF
     END IF

     IF( ListGetLogical( Params,'MMG Remesh', Found ) ) THEN
       CALL Remesh(CurrentModel,Solver)
     END IF

     CALL Info('DefaultFinish','Finished solver: '//GetString(Params,'Equation'),Level=8)
     
!------------------------------------------------------------------------------
   END SUBROUTINE DefaultFinish
!------------------------------------------------------------------------------


!> Solver the matrix equation related to the active solver
!------------------------------------------------------------------------------
  RECURSIVE FUNCTION DefaultSolve( USolver, BackRotNT ) RESULT(Norm)
!------------------------------------------------------------------------------
    TYPE(Solver_t), OPTIONAL, TARGET, INTENT(in) :: USolver
    REAL(KIND=dp) :: Norm
    LOGICAL, OPTIONAL, INTENT(in) :: BackRotNT

    TYPE(Matrix_t), POINTER   :: A
    TYPE(Variable_t), POINTER :: x
    REAL(KIND=dp), POINTER CONTIG :: b(:)
    REAL(KIND=dp), POINTER CONTIG :: SOL(:)

    LOGICAL :: Found, BackRot

    TYPE(ValueList_t), POINTER :: Params
    TYPE(Solver_t), POINTER :: Solver
    TYPE(Matrix_t), POINTER :: Ctmp
    CHARACTER(:), ALLOCATABLE :: linsolver, precond, dumpfile, saveslot
    INTEGER :: NameSpaceI, Count, MaxCount, i
    LOGICAL :: LinearSystemTrialing, SourceControl, NonlinearControl, &
        MonolithicSlave
    REAL(KIND=dp) :: s(3)
    
    CALL Info('DefaultSolve','Solving linear system with default routines',Level=10)
    
    Solver => CurrentModel % Solver
    Norm = REAL(0, dp)
    IF ( PRESENT( USolver ) ) Solver => USolver

    Params => GetSolverParams(Solver)

    NameSpaceI = NINT( ListGetCReal( Params,'Linear System Namespace Number', Found ) )
    LinearSystemTrialing = ListGetLogical( Params,'Linear System Trialing', Found )
    IF( LinearSystemTrialing ) NameSpaceI = MAX( 1, NameSpaceI )
      
    IF( NameSpaceI > 0 ) THEN
      CALL Info('DefaultSolve','Linear system namespace number: '//I2S(NameSpaceI),Level=7)
      CALL ListPushNamespace('linsys'//I2S(NameSpaceI)//':')
    END IF

    IF ( ListGetLogical( Params,'Linear System Save',Found )) THEN
      saveslot = GetString( Params,'Linear System Save Slot', Found )
      IF(.NOT. Found .OR. saveslot == 'solve') THEN
        CALL SaveLinearSystem( Solver ) 
      END IF
    END IF
    
    IF (PRESENT(BackRotNT)) THEN
      BackRot=GetLogical(Params,'Back Rotate N-T Solution',Found)
      IF(.NOT.Found) BackRot=.TRUE.

      IF (BackRot.NEQV.BackRotNT) &
        CALL ListAddLogical(Params,'Back Rotate N-T Solution',BackRotNT)
    END IF

    MonolithicSlave = ListGetLogical(Params,'Monolithic Slave',Found )
    IF( MonolithicSlave ) THEN      
      CALL MergeSlaveSolvers( Solver, PreSolve = .TRUE.)
    END IF
           
    IF( ListGetLogical( Params,'Harmonic Mode',Found ) ) THEN
      CALL ChangeToHarmonicSystem( Solver )
    END IF

    ! Generate projector that allows enforcing of total flux when using Robin BC's
    CALL GenerateRobinProjectors( CurrentModel, Solver )
    
    ! Combine the individual projectors into one massive projector
    CALL GenerateConstraintMatrix( CurrentModel, Solver )
    
    IF( GetLogical(Params,'Linear System Solver Disabled',Found) ) THEN
      CALL Info('DefaultSolve','Solver disabled, exiting early!',Level=10)
      RETURN
    END IF    

    SourceControl = ListGetLogical( Params,'Apply Source Control',Found )
    IF(SourceControl) CALL ControlLinearSystem( Solver,PreSolve=.TRUE. ) 
    
    NonlinearControl = ListGetLogical( Params,'Apply Nonlinear Control',Found )
    IF(NonlinearControl) CALL ControlNonlinearSystem( Solver, PreSolve=.TRUE.)

    
    CALL Info('DefaultSolve','Calling SolveSystem for linear solution',Level=20)

    A => Solver % Matrix
    x => Solver % Variable    
    b => A % RHS
    SOL => x % Values

    ! Debugging stuff activated only when "Max Output Level" >= 20
    IF( InfoActive( 20 ) ) THEN
      CALL VectorValuesRange(A % Values,SIZE(A % Values),'A')       
      CALL VectorValuesRange(A % rhs,SIZE(A % rhs),'b')       
    END IF
      
10  CONTINUE

    CALL SolveSystem(A,ParMatrix,b,SOL,x % Norm,x % DOFs,Solver)
    
    IF( InfoActive( 20 ) ) THEN
      CALL VectorValuesRange(x % Values,SIZE(x % values),'x')       
    END IF
    
    IF( LinearSystemTrialing ) THEN
      IF( x % LinConverged > 0 ) THEN
        IF( ListGetLogical( Params,'Linear System Trialing Conserve',Found ) ) THEN
          MaxCount = ListGetInteger( Params,'Linear System Trialing Conserve Rounds',Found ) 
          IF( Found ) THEN
            i = NINT( ListGetConstReal( Params,'Linear System Namespace Number',Found ) )
            IF( i == NameSpaceI ) THEN
              Count = 1 + ListGetInteger( Params,'Linear System Namespace Conserve Count',Found )
            ELSE
              Count = 0
            END IF
            IF( Count > MaxCount ) THEN
              NameSpaceI = 0
              Count = 0
            END IF
            CALL ListAddInteger( Params,'Linear System Namespace Conserve Count',Count )            
          END IF
          CALL ListAddConstReal( Params,'Linear System Namespace Number', 1.0_dp *NameSpaceI )
        END IF
      ELSE
        NameSpaceI = NameSpaceI + 1      
        IF( .NOT. ListCheckPrefix( Params,'linsys'//I2S(NameSpaceI) ) ) THEN
          CALL Fatal('DefaultSolve','Exhausted all linear system strategies!')
        END IF
        CALL ListPopNamespace()
        CALL Info('DefaultSolve','Linear system namespace number: '//I2S(NameSpaceI),Level=7)
        CALL ListPushNamespace('linsys'//I2S(NameSpaceI)//':')
        GOTO 10
      END IF
    END IF
    
    IF(SourceControl) CALL ControlLinearSystem( Solver,PreSolve=.FALSE. ) 
    IF(NonlinearControl) CALL ControlNonlinearSystem(Solver,PreSolve=.FALSE.)
    
    IF ( ListGetLogical( Params,'Linear System Save',Found )) THEN
      saveslot = GetString( Params,'Linear System Save Slot', Found )
      IF( Found .AND. saveslot == 'after') THEN
        CALL SaveLinearSystem( Solver ) 
      END IF
    END IF

    
    ! If flux corrected transport is used then apply the corrector to the system
    IF( GetLogical( Params,'Linear System FCT',Found ) ) THEN
      CALL FCT_Correction( Solver )
    END IF

    IF( MonolithicSlave ) THEN      
      CALL MergeSlaveSolvers( Solver, PreSolve = .FALSE.)
    END IF
    
    ! Backchange the linear system 
    IF( ListGetLogical( Params,'Harmonic Mode',Found ) ) THEN
      CALL ChangeToHarmonicSystem( Solver, .TRUE. )
    END IF
    
    IF (PRESENT(BackRotNT)) THEN
      IF (BackRot.NEQV.BackRotNT) &
        CALL ListAddLogical(Params,'Back Rotate N-T Solution',BackRot)
    END IF

    Norm = x % Norm

    IF( NameSpaceI > 0 ) CALL ListPopNamespace()
    
    ! One can run postprocessing solver in this slot in every nonlinear iteration.
    !-----------------------------------------------------------------------------
    CALL DefaultSlaveSolvers(Solver,'Nonlinear Post Solvers')


    ! This could be somewhere else too. Now it is here for debugging.
    CALL SaveParallelInfo( Solver )
    
!------------------------------------------------------------------------------
  END FUNCTION DefaultSolve
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Is the system converged. Wrapper to hide the dirty test.
!------------------------------------------------------------------------------
  FUNCTION DefaultConverged( USolver ) RESULT( Converged ) 
!------------------------------------------------------------------------------
    TYPE(Solver_t), OPTIONAL, TARGET, INTENT(in) :: USolver
    TYPE(Solver_t), POINTER :: Solver
    LOGICAL :: Converged
    LOGICAL :: Found
    INTEGER :: i,imin,imax
    
    Solver => CurrentModel % Solver
    IF ( PRESENT( USolver ) ) Solver => USolver

    IF( ListGetLogical( CurrentModel % Simulation,'Parallel Timestepping',Found ) ) THEN
      i = Solver % Variable % NonlinConverged
      CALL Info('DefaultConverged','Convergence status: '//I2S(i),Level=12)      
      imin = ParallelReduction(i,1)
      imax = ParallelReduction(i,2)
      IF(imin /= imax ) THEN
        CALL Info('DefaultConverged','Parallel timestepping converging at different rates!',Level=6)
        Solver % Variable % NonlinConverged = imin
      END IF
    END IF

    Converged = ( Solver % Variable % NonlinConverged > 0 )
          
  END FUNCTION DefaultConverged
!------------------------------------------------------------------------------
         

!------------------------------------------------------------------------------
  FUNCTION DefaultLinesearch( Converged, USolver, FirstIter, nsize, values, values0 ) RESULT( ReduceStep ) 
!------------------------------------------------------------------------------
    LOGICAL, OPTIONAL :: Converged
    TYPE(Solver_t), OPTIONAL, TARGET :: USolver
    LOGICAL, OPTIONAL :: FirstIter 
    INTEGER, OPTIONAL :: nsize
    REAL(KIND=dp), OPTIONAL, TARGET :: values(:), values0(:)
    LOGICAL :: ReduceStep

    LOGICAL :: stat, First, Last, DoLinesearch
    TYPE(Solver_t), POINTER :: Solver
    TYPE(Variable_t), POINTER :: iterV
    INTEGER :: iter, previter, MaxIter
    REAL(KIND=dp) :: LinesearchCond

    SAVE :: previter

    IF ( PRESENT( USolver ) ) THEN
      Solver => USolver
    ELSE
      Solver => CurrentModel % Solver
    END IF

    DoLinesearch = .FALSE.
    IF( ListCheckPrefix( Solver % Values,'Nonlinear System Linesearch') ) THEN
      LineSearchCond = ListGetCReal( Solver % Values,&
          'Nonlinear System Linesearch Condition', Stat )
      IF( Stat ) THEN
        DoLinesearch = ( LineSearchCond > 0.0_dp )
        CALL ListAddLogical( Solver % Values,'Nonlinear System Linesearch', DoLinesearch )
      ELSE
        DoLinesearch = ListGetLogical( Solver % Values,'Nonlinear System Linesearch',Stat)
      END IF
    END IF

    ! This routine might be called for convenience also without checking 
    ! first whether it is needed.
    IF(.NOT. DoLinesearch ) THEN
      ReduceStep = .FALSE.
      IF( PRESENT( Converged ) ) Converged = .FALSE.
      RETURN
    END IF
    
    IF( PRESENT( FirstIter ) ) THEN
      First = FirstIter
      Last = .FALSE.
    ELSE
      ! This is the first trial if we are the first nonlinear iteration
      ! for the first time. 
      iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
      iter = NINT(iterV % Values(1))
      MaxIter = ListGetInteger( Solver % Values,'Nonlinear System Max Iterations',Stat) 
      First = (iter == 1 ) .AND. (iter /= previter)
      Last = (iter == MaxIter )
      previter = iter
    END IF

    ReduceStep = CheckStepSize(Solver,First,nsize,values,values0) 

    IF( Last .AND. .NOT. ReduceStep ) THEN
      CALL Info('DefaultLinesearch',&
          'Maximum number of nonlinear iterations reached, giving up after linesearch',Level=6)
    END IF

    IF( PRESENT( Converged ) ) THEN
      Converged = ( Solver % Variable % NonlinConverged == 1 )  .OR. Last
    END IF

  END FUNCTION DefaultLinesearch


 
!------------------------------------------------------------------------------
  SUBROUTINE DefaultUpdateEquationsR( G, F, UElement, USolver, VecAssembly ) 
!------------------------------------------------------------------------------
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     REAL(KIND=dp) :: G(:,:), f(:)
     LOGICAL, OPTIONAL :: VecAssembly

     TYPE(Solver_t), POINTER   :: Solver
     TYPE(Matrix_t), POINTER   :: A
     TYPE(Variable_t), POINTER :: x
     TYPE(Element_t), POINTER  :: Element, P1, P2
     REAL(KIND=dp), POINTER CONTIG   :: b(:), svalues(:)

     LOGICAL :: Found, VecAsm, MCAsm

     INTEGER :: i, j, n, nd
     INTEGER(KIND=AddrInt) :: Proc
     INTEGER, POINTER CONTIG :: Indexes(:), PermIndexes(:)

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
     ELSE
        Solver => CurrentModel % Solver
     END IF
     A => Solver % Matrix
     x => Solver % Variable
     b => A % RHS

     Element => GetCurrentElement( UElement )
     
     VecAsm = .FALSE.
     IF ( PRESENT( VecAssembly )) THEN
       VecAsm = VecAssembly
     END IF

     IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
       Proc = Solver % BoundaryElementProcedure
     ELSE
       Proc = Solver % BulkElementProcedure
     END IF
     IF ( Proc /= 0 ) THEN
       n  = GetElementNOFNodes( Element )
       nd = GetElementNOFDOFs( Element, Solver )
       CALL ExecLocalProc( Proc, CurrentModel, Solver, &
           G, F, Element, n, nd )
     END IF

     
     IF ( ParEnv % PEs > 1 ) THEN
       IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
          P1 => Element % BoundaryInfo % Left
          P2 => Element % BoundaryInfo % Right
          IF ( ASSOCIATED(P1) .AND. ASSOCIATED(P2) ) THEN
            IF ( P1 % PartIndex /= ParEnv % myPE .AND. &
                 P2 % PartIndex /= ParEnv % myPE )RETURN

            IF ( P1 % PartIndex /= ParEnv % myPE .OR. &
                 P2 % PartIndex /= ParEnv % myPE ) THEN
              G=G/2; F=F/2; 
            END IF
          ELSE IF ( ASSOCIATED(P1) ) THEN
            IF ( P1 % PartIndex /= ParEnv % myPE ) RETURN
          ELSE IF ( ASSOCIATED(P2) ) THEN
            IF ( P2 % PartIndex /= ParEnv % myPE ) RETURN
          END IF
       ELSE IF ( Element % PartIndex/=ParEnv % myPE ) THEN
          IF(GetLogical(Solver % Values,'Linear System FCT',Found)) THEN
            Indexes => GetIndexStore()
            n = GetElementDOFs( Indexes, Element, Solver )
            IF(.NOT.ASSOCIATED(A % HaloValues)) THEN
              ALLOCATE(A % HaloValues(SIZE(A % Values))); A % HaloValues=0._dp
            END IF
            CALL UpdateGlobalEquations( A,G,b,0._dp*f,n,x % DOFs, &
              x % Perm(Indexes(1:n)),UElement=Element,GlobalValues=A % HaloValues )
            END IF
            RETURN
       END IF
     END IF

     ! Vectorized version of the glueing process requested
     IF (VecAsm) THEN
#ifdef _OPENMP
       IF (OMP_GET_NUM_THREADS() == 1) THEN
         MCAsm = .TRUE.
       ELSE
         ! Check if multicoloured assembly is in use
         MCAsm = (Solver % CurrentColour > 0) .AND. &
                 ASSOCIATED(Solver % ColourIndexList)
       END IF
#else
       MCAsm = .TRUE.
#endif
     ELSE
       MCAsm = .FALSE.
     END IF
       
     Indexes => GetIndexStore()
     n = GetElementDOFs( Indexes, Element, Solver )
       
     PermIndexes => GetPermIndexStore()
     ! Get permuted indices
!DIR$ IVDEP
     DO j=1,n
       PermIndexes(j) = x % Perm(Indexes(j))
     END DO

     IF( Solver % LocalSystemMode > 0 ) THEN
       CALL UseLocalMatrixStorage( Solver, n * x % dofs, G, F, ElemInd = Element % ElementIndex )
     END IF
     
     ! If we have any antiperiodic entries we need to check them all!
     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalMatrix( Solver, n, Indexes, x % dofs, G )
       CALL FlipPeriodicLocalForce( Solver, n, Indexes, x % dofs, f )
     END IF

     IF( VecAsm ) THEN
       CALL UpdateGlobalEquationsVec( A, G, b, f, n, &
           x % DOFs, PermIndexes, &
           UElement=Element, MCAssembly=MCAsm )
     ELSE       
!      IF( A % FORMAT == MATRIX_CRS ) THEN
!        ! For CRS format these are effectively the same
!        CALL UpdateGlobalEquationsVec( A,G,b,f,n,x % DOFs, &
!            PermIndexes, UElement=Element )
!      ELSE
         CALL UpdateGlobalEquations( A,G,b,f,n,x % DOFs, &
             PermIndexes, UElement=Element )       
!      END IF
         
       IF(Solver % DirectMethod == DIRECT_PERMON) THEN
         CALL UpdatePermonMatrix( A, G, n, x % DOFs, PermIndexes )
       END IF
     END IF
     
     ! backflip, in case G is needed again
     ! For change of sign backflip and flip are same operations.
     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalMatrix( Solver, n, Indexes, x % dofs, G )
       CALL FlipPeriodicLocalForce( Solver, n, Indexes, x % dofs, f )
     END IF
     
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdateEquationsR
!------------------------------------------------------------------------------

  
  SUBROUTINE DefaultUpdateEquationsIm( G, F, UElement, USolver, VecAssembly )     
    TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
    TYPE(Element_t), OPTIONAL, TARGET :: UElement
    REAL(KIND=dp) :: G(:,:), f(:)
    LOGICAL, OPTIONAL :: VecAssembly

    TYPE(Matrix_t), POINTER   :: A
    REAL(KIND=dp), POINTER :: pvalues(:), prhs(:)
    
    IF ( PRESENT( USolver ) ) THEN
      A => USolver % Matrix
    ELSE
      A => CurrentModel % Solver % Matrix
    END IF
    
    IF(.NOT. ASSOCIATED( A % Values_im ) ) THEN
      ALLOCATE( A % Values_im(SIZE( A % Values ) ) )
      A % Values_im = 0.0_dp
    END IF
    pvalues => A % Values
    A % Values => A % Values_im    
    
    IF(.NOT. ASSOCIATED( A % rhs_im ) ) THEN
      ALLOCATE( A % rhs_im(SIZE( A % rhs ) ) )
      A % rhs_im = 0.0_dp
    END IF
    prhs => A % Rhs
    A % rhs => A % rhs_im    
    
    CALL DefaultUpdateEquationsR( G, F, UElement, USolver, VecAssembly )     
    
    A % Values => pValues
    A % rhs => prhs
    
  END SUBROUTINE DefaultUpdateEquationsIm
  

!------------------------------------------------------------------------------
 SUBROUTINE UpdatePermonMatrix(A,G,n,dofs,nind)
!------------------------------------------------------------------------------
#ifdef HAVE_FETI4I
   use feti4i
#endif

   TYPE(Matrix_t) :: A
   INTEGER :: n, dofs, nInd(:)
   REAL(KIND=dp) :: G(:,:)
!------------------------------------------------------------------------------
  REAL(KIND=C_DOUBLE), ALLOCATABLE :: vals(:)
  INTEGER, POINTER :: ptr
  INTEGER :: i,j,k,l,k1,k2
  INTEGER :: matrixType, eType
  INTEGER(C_INT), ALLOCATABLE :: ind(:)

  TYPE(Element_t), POINTER :: CElement
  
#ifdef HAVE_FETI4I
!!$  INTERFACE
!!$     FUNCTION Permon_InitMatrix(n) RESULT(handle) BIND(C,Name="permon_init")
!!$       USE, INTRINSIC :: ISO_C_BINDING
!!$       TYPE(C_PTR) :: Handle
!!$       INTEGER(C_INT), VALUE :: n
!!$     END FUNCTION Permon_InitMatrix
!!$
!!$     SUBROUTINE Permon_UpdateMatrix(handle,n,inds,vals) BIND(C,Name="permon_update")
!!$       USE, INTRINSIC :: ISO_C_BINDING
!!$       TYPE(C_PTR), VALUE :: Handle
!!$       INTEGER(C_INT), VALUE :: n
!!$       INTEGER(C_INT) :: inds(*)
!!$       REAL(C_DOUBLE) :: vals(*)
!!$     END SUBROUTINE Permon_UpdateMatrix
!!$  END INTERFACE

  IF(.NOT. C_ASSOCIATED(A % PermonMatrix)) THEN
    A % NoDirichlet = .TRUE.
    !! A % PermonMatrix = Permon_InitMatrix(A % NumberOFRows)
    !! TODO: get correct matrix type 
    matrixType = 0  !! symmetric positive definite (for other types see feti4i.h)
    CALL FETI4ICreateStiffnessMatrix(A % PermonMatrix, matrixType, 1) !TODO add number of rows A % NumberOFRows
  END IF

  
  ALLOCATE(vals(n*n*dofs*dofs), ind(n*dofs))
  DO i=1,n
    DO j=1,dofs
      k1 = (i-1)*dofs + j
      DO k=1,n
        DO l=1,dofs
           k2 = (k-1)*dofs + l
           vals(dofs*n*(k1-1)+k2) = G(k1,k2)
        END DO
      END DO
      ind(k1) = dofs*(nInd(i)-1)+j
    END DO
  END DO

  !CALL Permon_UpdateMatrix( A % PermonMatrix, n*dofs, ind, vals )

  CElement => GetCurrentElement()
  eType = ElementDim( CElement )
  ! type of the element is the same as its dimension

  CALL FETI4IAddElement(A % PermonMatrix, eType, n, nInd, n*dofs, ind, vals)

#endif
    
!------------------------------------------------------------------------------
 END SUBROUTINE UpdatePermonMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE DefaultUpdateEquationsC( GC, FC, UElement, USolver, VecAssembly ) 
!------------------------------------------------------------------------------
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     COMPLEX(KIND=dp)   :: GC(:,:), FC(:)
     LOGICAL, OPTIONAL :: VecAssembly  ! The complex version lacks support for this 

     TYPE(Solver_t), POINTER   :: Solver
     TYPE(Matrix_t), POINTER   :: A
     TYPE(Variable_t), POINTER :: x
     TYPE(Element_t), POINTER  :: Element, P1, P2
     REAL(KIND=dp), POINTER  CONTIG :: b(:)

     REAL(KIND=dp), POINTER :: G(:,:), F(:)

     LOGICAL :: Found

     INTEGER :: i,j,n,DOFs
     INTEGER, POINTER :: Indexes(:)

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
     ELSE
        Solver => CurrentModel % Solver
     END IF
     A => Solver % Matrix
     x => Solver % Variable
     b => A % RHS

     Element => GetCurrentElement( UElement )

     DOFs = x % DOFs
     Indexes => GetIndexStore()
     n = GetElementDOFs( Indexes, Element, Solver )

     IF ( ParEnv % PEs > 1 ) THEN
       IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
          P1 => Element % BoundaryInfo % Left
          P2 => Element % BoundaryInfo % Right
          IF ( ASSOCIATED(P1) .AND. ASSOCIATED(P2) ) THEN
            IF ( P1 % PartIndex/=ParEnv % myPE .AND. &
                 P2 % PartIndex/=ParEnv % myPE )RETURN

            IF ( P1 % PartIndex/=ParEnv % myPE .OR. &
                 P2 % PartIndex/=ParEnv % myPE ) THEN
              GC=GC/2; FC=FC/2; 
            END IF
          ELSE IF ( ASSOCIATED(P1) ) THEN
            IF ( P1 % PartIndex/=ParEnv % myPE ) RETURN
          ELSE IF ( ASSOCIATED(P2) ) THEN
            IF ( P2 % PartIndex/=ParEnv % myPE ) RETURN
          END IF
       ELSE IF ( Element % PartIndex/=ParEnv % myPE ) THEN
          RETURN
       END IF
     END IF

     ALLOCATE( G(DOFs*n,DOFs*n), F(DOFs*n) )
     DO i=1,n*DOFs/2
       F( 2*(i-1)+1 ) =  REAL( FC(i) )
       F( 2*(i-1)+2 ) = AIMAG( FC(i) )

       DO j=1,n*DOFs/2
         G( 2*(i-1)+1, 2*(j-1)+1 ) =   REAL( GC(i,j) )
         G( 2*(i-1)+1, 2*(j-1)+2 ) = -AIMAG( GC(i,j) )
         G( 2*(i-1)+2, 2*(j-1)+1 ) =  AIMAG( GC(i,j) )
         G( 2*(i-1)+2, 2*(j-1)+2 ) =   REAL( GC(i,j) )
       END DO
     END DO

     IF( Solver % LocalSystemMode > 0 ) THEN
       CALL UseLocalMatrixStorage( Solver, n * x % dofs, G, F, ElemInd = Element % ElementIndex )
     END IF
               
     ! If we have any antiperiodic entries we need to check them all!
     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalMatrix( Solver, n, Indexes, x % dofs, G )
       CALL FlipPeriodicLocalForce( Solver, n, Indexes, x % dofs, f )
     END IF
          
     CALL UpdateGlobalEquations( A,G,b,f,n,x % DOFs,x % Perm(Indexes(1:n)) )

     DEALLOCATE( G, F)
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdateEquationsC
!------------------------------------------------------------------------------


! This is a version when the initial matrix is given in diagonal form,
! such that the last array index refers to the component.
!------------------------------------------------------------------------------
  SUBROUTINE DefaultUpdateEquationsDiagC( GC, FC, UElement, USolver, VecAssembly ) 
!------------------------------------------------------------------------------
    TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
    TYPE(Element_t), OPTIONAL, TARGET :: UElement
    COMPLEX(KIND=dp)   :: GC(:,:,:), FC(:,:)
    LOGICAL, OPTIONAL :: VecAssembly  ! The complex version lacks support for this 

    TYPE(Solver_t), POINTER   :: Solver
    TYPE(Matrix_t), POINTER   :: A
    TYPE(Variable_t), POINTER :: x
    TYPE(Element_t), POINTER  :: Element, P1, P2
    REAL(KIND=dp), POINTER  CONTIG :: b(:)

    REAL(KIND=dp), POINTER :: G(:,:), F(:)

    LOGICAL :: Found, Half

    INTEGER :: i,j,k,n,DOFs
    INTEGER, POINTER :: Indexes(:)

    IF ( PRESENT( USolver ) ) THEN
      Solver => USolver
    ELSE
      Solver => CurrentModel % Solver
    END IF
    A => Solver % Matrix
    x => Solver % Variable
    b => A % RHS

    Element => GetCurrentElement( UElement )

    DOFs = x % DOFs
    Indexes => GetIndexStore()
    n = GetElementDOFs( Indexes, Element, Solver )
    
    Half = .FALSE.
    IF ( ParEnv % PEs > 1 ) THEN
      IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
        P1 => Element % BoundaryInfo % Left
        P2 => Element % BoundaryInfo % Right
        IF ( ASSOCIATED(P1) .AND. ASSOCIATED(P2) ) THEN
          IF ( P1 % PartIndex/=ParEnv % myPE .AND. &
              P2 % PartIndex/=ParEnv % myPE )RETURN

          IF ( P1 % PartIndex/=ParEnv % myPE .OR. &
              P2 % PartIndex/=ParEnv % myPE ) THEN
            Half = .TRUE.
          END IF
        ELSE IF ( ASSOCIATED(P1) ) THEN
          IF ( P1 % PartIndex/=ParEnv % myPE ) RETURN
        ELSE IF ( ASSOCIATED(P2) ) THEN
          IF ( P2 % PartIndex/=ParEnv % myPE ) RETURN
        END IF
      ELSE IF ( Element % PartIndex/=ParEnv % myPE ) THEN
        RETURN
      END IF
    END IF

    ALLOCATE( G(DOFs*n,DOFs*n), F(DOFs*n) )
    G = 0.0_dp; F = 0.0_dp

    DO i=1,n
      DO k=1,dofs/2 
        F( dofs*(i-1)+2*k-1) = REAL( FC(i,k) )
        F( dofs*(i-1)+2*k ) = AIMAG( FC(i,k) )
      END DO
    END DO
    
    DO i=1,n
      DO j=1,n
        DO k=1,DOFs/2
          G( dofs*(i-1)+2*k-1, dofs*(j-1)+2*k-1 ) = REAL( GC(i,j,k) )
          G( dofs*(i-1)+2*k-1, dofs*(j-1)+2*k ) = -AIMAG( GC(i,j,k) )
          G( dofs*(i-1)+2*k, dofs*(j-1)+2*k-1 ) =  AIMAG( GC(i,j,k) )
          G( dofs*(i-1)+2*k, dofs*(j-1)+2*k ) = REAL( GC(i,j,k) )
        END DO
      END DO
    END DO

    ! Scale only the temporal fields
    IF( Half ) THEN
      G = G/2; F = F/2 
    END IF

    ! If we have any antiperiodic entries we need to check them all!
    IF( Solver % PeriodicFlipActive ) THEN
      CALL FlipPeriodicLocalMatrix( Solver, n, Indexes, x % dofs, G )
      CALL FlipPeriodicLocalForce( Solver, n, Indexes, x % dofs, f )
    END IF

    CALL UpdateGlobalEquations( A,G,b,f,n,x % DOFs,x % Perm(Indexes(1:n)) )

    DEALLOCATE( G, F)
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdateEquationsDiagC
!------------------------------------------------------------------------------


  
!> Adds the elementwise contribution the right-hand-side of the real valued matrix equation 
!------------------------------------------------------------------------------
  SUBROUTINE DefaultUpdateForceR( F, UElement, USolver, BulkUpdate )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: F(:)
    TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
    TYPE(Element_t), OPTIONAL, TARGET :: UElement
    LOGICAL, OPTIONAL :: BulkUpdate

    TYPE(Solver_t), POINTER :: Solver
    TYPE(Variable_t), POINTER :: x
    TYPE(Element_t), POINTER  :: Element, P1, P2

    LOGICAL :: Found

    INTEGER :: n
    INTEGER, POINTER :: Indexes(:)

    Solver => CurrentModel % Solver
    IF ( PRESENT(USolver) ) Solver => USolver

    Element => GetCurrentElement( UElement ) 

    x => Solver % Variable
    Indexes => GetIndexStore()
    n = GetElementDOFs( Indexes, Element, Solver )

     IF ( ParEnv % PEs > 1 ) THEN
       IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
          P1 => Element % BoundaryInfo % Left
          P2 => Element % BoundaryInfo % Right
          IF ( ASSOCIATED(P1) .AND. ASSOCIATED(P2) ) THEN
            IF ( P1 % PartIndex/=ParEnv % myPE .AND. &
                 P2 % PartIndex/=ParEnv % myPE )RETURN

            IF ( P1 % PartIndex/=ParEnv % myPE .OR. &
                 P2 % PartIndex/=ParEnv % myPE ) F=F/2; 
          ELSE IF ( ASSOCIATED(P1) ) THEN
            IF ( P1 % PartIndex /= ParEnv % myPE ) RETURN
          ELSE IF ( ASSOCIATED(P2) ) THEN
            IF ( P2 % PartIndex /= ParEnv % myPE ) RETURN
          END IF
       ELSE IF ( Element % PartIndex/=ParEnv % myPE ) THEN
          RETURN
       END IF
     END IF

     ! If we have any antiperiodic entries we need to check them all!
     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalForce( Solver, n, Indexes, x % dofs, f )
     END IF
     
     CALL UpdateGlobalForce( Solver % Matrix % RHS, &
         F, n, x % DOFs, x % Perm(Indexes(1:n)), UElement=Element)

     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalForce( Solver, n, Indexes, x % dofs, f )
     END IF
     

!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdateForceR
!------------------------------------------------------------------------------


!> Adds the elementwise contribution the right-hand-side of the complex valued matrix equation 
!------------------------------------------------------------------------------
  SUBROUTINE DefaultUpdateForceC( FC, UElement, USolver )
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: FC(:)
    TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
    TYPE(Element_t), OPTIONAL, TARGET :: UElement

    TYPE(Solver_t), POINTER :: Solver
    TYPE(Variable_t), POINTER :: x
    TYPE(Element_t), POINTER  :: Element, P1, P2

    REAL(KIND=dp), ALLOCATABLE :: F(:)

    INTEGER :: i,n,DOFs
    INTEGER, POINTER :: Indexes(:)

    Solver => CurrentModel % Solver
    IF ( PRESENT(USolver) ) Solver => USolver

    Element => GetCurrentElement( UElement ) 

    x => Solver % Variable
    DOFs = x % DOFs
    Indexes => GetIndexStore()
    n = GetElementDOFs( Indexes, Element, Solver )

     IF ( ParEnv % PEs > 1 ) THEN
       IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
          P1 => Element % BoundaryInfo % Left
          P2 => Element % BoundaryInfo % Right
          IF ( ASSOCIATED(P1) .AND. ASSOCIATED(P2) ) THEN
            IF ( P1 % PartIndex/=ParEnv % myPE .AND. &
                 P2 % PartIndex/=ParEnv % myPE )RETURN

            IF ( P1 % PartIndex/=ParEnv % myPE .OR. &
                 P2 % PartIndex/=ParEnv % myPE ) FC=FC/2; 
          ELSE IF ( ASSOCIATED(P1) ) THEN
            IF ( P1 % PartIndex/=ParEnv % myPE ) RETURN
          ELSE IF ( ASSOCIATED(P2) ) THEN
            IF ( P2 % PartIndex/=ParEnv % myPE ) RETURN
          END IF
       ELSE IF ( Element % PartIndex/=ParEnv % myPE ) THEN
          RETURN
       END IF
     END IF


    ALLOCATE( F(DOFs*n) )
    DO i=1,n*DOFs/2
       F( 2*(i-1) + 1 ) =   REAL(FC(i))
       F( 2*(i-1) + 2 ) = AIMAG(FC(i))
    END DO

    IF( Solver % PeriodicFlipActive ) THEN
      CALL FlipPeriodicLocalForce( Solver, n, Indexes, x % dofs, f )
    END IF
            
    CALL UpdateGlobalForce( Solver % Matrix % RHS, &
        F, n, x % DOFs, x % Perm(Indexes(1:n)) )

    DEALLOCATE( F ) 
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdateForceC
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE DefaultUpdateTimeForceR( F, UElement, USolver )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: F(:)
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
 
     TYPE(Solver_t), POINTER :: Solver
     TYPE(Variable_t), POINTER :: x
     TYPE(Element_t), POINTER  :: Element, P1, P2
 
     INTEGER :: n
     INTEGER, POINTER :: Indexes(:)

     Solver => CurrentModel % Solver
     IF ( PRESENT(USolver) ) Solver => USolver

     Element => GetCurrentElement( UElement ) 

     x => Solver % Variable
     Indexes => GetIndexStore()
     n = GetElementDOFs( Indexes, Element, Solver )

     IF ( ParEnv % PEs > 1 ) THEN
       IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
          P1 => Element % BoundaryInfo % Left
          P2 => Element % BoundaryInfo % Right
          IF ( ASSOCIATED(P1) .AND. ASSOCIATED(P2) ) THEN
            IF ( P1 % PartIndex/=ParEnv % myPE .AND. &
                 P2 % PartIndex/=ParEnv % myPE )RETURN

            IF ( P1 % PartIndex/=ParEnv % myPE .OR. &
                 P2 % PartIndex/=ParEnv % myPE ) F=F/2; 
          ELSE IF ( ASSOCIATED(P1) ) THEN
            IF ( P1 % PartIndex /= ParEnv % myPE ) RETURN
          ELSE IF ( ASSOCIATED(P2) ) THEN
            IF ( P2 % PartIndex /= ParEnv % myPE ) RETURN
          END IF
       ELSE IF ( Element % PartIndex/=ParEnv % myPE ) THEN
          RETURN
       END IF
     END IF

     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalForce( Solver, n, Indexes, x % dofs, f )
     END IF
          
     CALL UpdateTimeForce( Solver % Matrix,Solver % Matrix % RHS, &
         F, n, x % DOFs, x % Perm(Indexes(1:n)) )

     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalForce( Solver, n, Indexes, x % dofs, f )
     END IF     
     
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdateTimeForceR
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE DefaultUpdateTimeForceC( FC, UElement, USolver )
!------------------------------------------------------------------------------
     COMPLEX(KIND=dp) :: FC(:)
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement

     TYPE(Solver_t), POINTER :: Solver
     TYPE(Variable_t), POINTER :: x
     TYPE(Element_t), POINTER  :: Element, P1, P2
 
     REAL(KIND=dp), ALLOCATABLE :: F(:)

    INTEGER :: i,n,DOFs
    INTEGER, POINTER :: Indexes(:)

     Solver => CurrentModel % Solver
     IF ( PRESENT(USolver) ) Solver => USolver

     Element => GetCurrentElement( UElement ) 

      x => Solver % Variable
      DOFs = x % DOFs
      Indexes => GetIndexStore()
      n = GetElementDOFs( Indexes, Element, Solver )

      IF ( ParEnv % PEs > 1 ) THEN
        IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
           P1 => Element % BoundaryInfo % Left
           P2 => Element % BoundaryInfo % Right
           IF ( ASSOCIATED(P1) .AND. ASSOCIATED(P2) ) THEN
             IF ( P1 % PartIndex/=ParEnv % myPE .AND. &
                  P2 % PartIndex/=ParEnv % myPE )RETURN
 
             IF ( P1 % PartIndex/=ParEnv % myPE .OR. &
                  P2 % PartIndex/=ParEnv % myPE ) FC=FC/2; 
           ELSE IF ( ASSOCIATED(P1) ) THEN
             IF ( P1 % PartIndex/=ParEnv % myPE ) RETURN
           ELSE IF ( ASSOCIATED(P2) ) THEN
             IF ( P2 % PartIndex/=ParEnv % myPE ) RETURN
           END IF
        ELSE IF ( Element % PartIndex/=ParEnv % myPE ) THEN
           RETURN
        END IF
      END IF

      ALLOCATE( F(DOFs*n) )
      DO i=1,n*DOFs/2
         F( 2*(i-1) + 1 ) =   REAL(FC(i))
         F( 2*(i-1) + 2 ) = -AIMAG(FC(i))
      END DO
      
      IF( Solver % PeriodicFlipActive ) THEN
        CALL FlipPeriodicLocalForce( Solver, n, Indexes, x % dofs, f )
      END IF      
      
      CALL UpdateTimeForce( Solver % Matrix,Solver % Matrix % RHS, &
               F, n, x % DOFs, x % Perm(Indexes(1:n)) )

      DEALLOCATE( F ) 
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdateTimeForceC
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE DefaultUpdatePrecR( M, UElement, USolver ) 
!------------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL,TARGET   :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     REAL(KIND=dp)   :: M(:,:)

     TYPE(Solver_t), POINTER   :: Solver
     TYPE(Matrix_t), POINTER   :: A
     TYPE(Variable_t), POINTER :: x
     TYPE(Element_t), POINTER  :: Element, P1, P2

     INTEGER :: i,j,n
     INTEGER, POINTER :: Indexes(:)

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
        A => Solver % Matrix
        x => Solver % Variable
     ELSE
        Solver => CurrentModel % Solver
        A => Solver % Matrix
        x => Solver % Variable
     END IF

     Element => GetCurrentElement( UElement ) 

     Indexes => GetIndexStore()
     n = GetElementDOFs( Indexes, Element, Solver )

     IF ( ParEnv % PEs > 1 ) THEN
       IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
          P1 => Element % BoundaryInfo % Left
          P2 => Element % BoundaryInfo % Right
          IF ( ASSOCIATED(P1) .AND. ASSOCIATED(P2) ) THEN
            IF ( P1 % PartIndex/=ParEnv % myPE .AND. &
                 P2 % PartIndex/=ParEnv % myPE )RETURN

            IF ( P1 % PartIndex/=ParEnv % myPE .OR. &
                 P2 % PartIndex/=ParEnv % myPE ) M=M/2
          ELSE IF ( ASSOCIATED(P1) ) THEN
            IF ( P1 % PartIndex /= ParEnv % myPE ) RETURN
          ELSE IF ( ASSOCIATED(P2) ) THEN
            IF ( P2 % PartIndex /= ParEnv % myPE ) RETURN
          END IF
       ELSE IF ( Element % PartIndex/=ParEnv % myPE ) THEN
          RETURN
       END IF
     END IF

!$OMP CRITICAL
     IF ( .NOT. ASSOCIATED( A % PrecValues ) ) THEN
       CALL Info('DefaultUpdatePrecR','Allocating for separate preconditioning matrix!',Level=20)
       ALLOCATE( A % PrecValues(SIZE(A % Values)) )
       A % PrecValues = 0.0d0
     END IF
!$OMP END CRITICAL

     ! flip mass matrix for periodic elimination
     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalMatrix( Solver, n, Indexes, x % dofs, M )
     END IF
      
     CALL UpdateMassMatrix( A, M, n, x % DOFs, x % Perm(Indexes(1:n)), & 
         A % PrecValues )

     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalMatrix( Solver, n, Indexes, x % dofs, M )
     END IF
          
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdatePrecR
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE DefaultUpdatePrecC( MC, UElement, USolver ) 
!------------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL,TARGET   :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     COMPLEX(KIND=dp)   :: MC(:,:)

     TYPE(Solver_t), POINTER   :: Solver
     TYPE(Matrix_t), POINTER   :: A
     TYPE(Variable_t), POINTER :: x
     TYPE(Element_t), POINTER  :: Element, P1, P2

     REAL(KIND=dp), ALLOCATABLE :: M(:,:)

     INTEGER :: i,j,n,DOFs
     INTEGER, POINTER :: Indexes(:)

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
        A => Solver % Matrix
        x => Solver % Variable
     ELSE
        Solver => CurrentModel % Solver
        A => Solver % Matrix
        x => Solver % Variable
     END IF

     Element => GetCurrentElement( UElement ) 

     DOFs = x % DOFs
     Indexes => GetIndexStore()
     n = GetElementDOFs( Indexes, Element, Solver )

      IF ( ParEnv % PEs > 1 ) THEN
        IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
           P1 => Element % BoundaryInfo % Left
           P2 => Element % BoundaryInfo % Right
           IF ( ASSOCIATED(P1) .AND. ASSOCIATED(P2) ) THEN
             IF ( P1 % PartIndex/=ParEnv % myPE .AND. &
                  P2 % PartIndex/=ParEnv % myPE )RETURN
 
             IF ( P1 % PartIndex/=ParEnv % myPE .OR. &
                  P2 % PartIndex/=ParEnv % myPE ) MC=MC/2
           ELSE IF ( ASSOCIATED(P1) ) THEN
             IF ( P1 % PartIndex/=ParEnv % myPE ) RETURN
           ELSE IF ( ASSOCIATED(P2) ) THEN
             IF ( P2 % PartIndex/=ParEnv % myPE ) RETURN
           END IF
        ELSE IF ( Element % PartIndex/=ParEnv % myPE ) THEN
           RETURN
        END IF
      END IF

!$OMP CRITICAL
       IF ( .NOT. ASSOCIATED( A % PrecValues ) ) THEN
         CALL Info('DefaultUpdatePrecC','Allocating for separate preconditioning matrix!',Level=20)         
         ALLOCATE( A % PrecValues(SIZE(A % Values)) )
         A % PrecValues = 0.0d0
       END IF
!$OMP END CRITICAL

     ALLOCATE( M(DOFs*n,DOFs*n) )
     DO i=1,n*DOFs/2
       DO j=1,n*DOFs/2
         M(2*(i-1)+1, 2*(j-1)+1) =   REAL( MC(i,j) )
         M(2*(i-1)+1, 2*(j-1)+2) = -AIMAG( MC(i,j) )
         M(2*(i-1)+2, 2*(j-1)+1) =  AIMAG( MC(i,j) )
         M(2*(i-1)+2, 2*(j-1)+2) =   REAL( MC(i,j) )
       END DO
     END DO

     ! flip preconditioning matrix for periodic elimination
     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalMatrix( Solver, n, Indexes, x % dofs, M )
     END IF
     
     CALL UpdateMassMatrix( A, M, n, x % DOFs, x % Perm(Indexes(1:n)), &
              A % PrecValues )
     DEALLOCATE( M )
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdatePrecC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE DefaultUpdateMassR( M, UElement, USolver ) 
!------------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL,TARGET   :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     REAL(KIND=dp)   :: M(:,:)

     TYPE(Solver_t), POINTER   :: Solver
     TYPE(Matrix_t), POINTER   :: A
     TYPE(Variable_t), POINTER :: x
     TYPE(Element_t), POINTER  :: Element, P1, P2

     LOGICAL :: Found

     INTEGER :: i,j,n
     INTEGER, POINTER :: Indexes(:)

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
        A => Solver % Matrix
        x => Solver % Variable
     ELSE
        Solver => CurrentModel % Solver
        A => Solver % Matrix
        x => Solver % Variable
     END IF

     Element => GetCurrentElement( UElement ) 

     Indexes => GetIndexStore()
     n = GetElementDOFs( Indexes, Element, Solver )

     IF ( ParEnv % PEs > 1 ) THEN
       IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
          P1 => Element % BoundaryInfo % Left
          P2 => Element % BoundaryInfo % Right
          IF ( ASSOCIATED(P1) .AND. ASSOCIATED(P2) ) THEN
            IF ( P1 % PartIndex/=ParEnv % myPE .AND. &
                 P2 % PartIndex/=ParEnv % myPE )RETURN

            IF ( P1 % PartIndex/=ParEnv % myPE .OR. &
                 P2 % PartIndex/=ParEnv % myPE ) M=M/2
          ELSE IF ( ASSOCIATED(P1) ) THEN
            IF ( P1 % PartIndex /= ParEnv % myPE ) RETURN
          ELSE IF ( ASSOCIATED(P2) ) THEN
            IF ( P2 % PartIndex /= ParEnv % myPE ) RETURN
          END IF
       ELSE IF ( Element % PartIndex/=ParEnv % myPE ) THEN
          IF (ListGetLogical(Solver % Values, 'Linear System FCT', Found)) THEN
            Indexes => GetIndexStore()
            n = GetElementDOFs( Indexes, Element, Solver )
            IF(.NOT.ASSOCIATED(A % HaloMassValues)) THEN
              ALLOCATE(A % HaloMassValues(SIZE(A % Values))); A % HaloMassValues=0._dp
            END IF
            CALL UpdateMassMatrix( A, M, n, x % DOFs, x % Perm(Indexes(1:n)), &
                               A % HaloMassValues ) 
          END IF
          RETURN
       END IF
     END IF

!$OMP CRITICAL
     IF ( .NOT. ASSOCIATED( A % MassValues ) ) THEN
       ALLOCATE( A % MassValues(SIZE(A % Values)) )
       A % MassValues = 0.0_dp
     END IF
!$OMP END CRITICAL


     ! flip mass matrix for periodic elimination
     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalMatrix( Solver, n, Indexes, x % dofs, M )
     END IF
                
     CALL UpdateMassMatrix( A, M, n, x % DOFs, x % Perm(Indexes(1:n)), &
         A % MassValues ) 
     
     ! backflip to be on the safe side
     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalMatrix( Solver, n, Indexes, x % Dofs, M )
     END IF


!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdateMassR
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE DefaultUpdateMassC( MC, UElement, USolver ) 
!------------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL,TARGET   :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     COMPLEX(KIND=dp)   :: MC(:,:)

     TYPE(Solver_t), POINTER   :: Solver
     TYPE(Matrix_t), POINTER   :: A
     TYPE(Variable_t), POINTER :: x
     TYPE(Element_t), POINTER  :: Element, P1, P2

     REAL(KIND=dp), ALLOCATABLE :: M(:,:)

     INTEGER :: i,j,n,DOFs
     INTEGER, POINTER :: Indexes(:)

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
        A => Solver % Matrix
        x => Solver % Variable
     ELSE
        Solver => CurrentModel % Solver
        A => Solver % Matrix
        x => Solver % Variable
     END IF

     Element => GetCurrentElement( UElement ) 

     DOFs = x % DOFs
     Indexes => GetIndexStore()
     n = GetElementDOFs( Indexes, Element, Solver )

      IF ( ParEnv % PEs > 1 ) THEN
        IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
           P1 => Element % BoundaryInfo % Left
           P2 => Element % BoundaryInfo % Right
           IF ( ASSOCIATED(P1) .AND. ASSOCIATED(P2) ) THEN
             IF ( P1 % PartIndex/=ParEnv % myPE .AND. &
                  P2 % PartIndex/=ParEnv % myPE )RETURN
 
             IF ( P1 % PartIndex/=ParEnv % myPE .OR. &
                  P2 % PartIndex/=ParEnv % myPE ) MC=MC/2
           ELSE IF ( ASSOCIATED(P1) ) THEN
             IF ( P1 % PartIndex/=ParEnv % myPE ) RETURN
           ELSE IF ( ASSOCIATED(P2) ) THEN
             IF ( P2 % PartIndex/=ParEnv % myPE ) RETURN
           END IF
        ELSE IF ( Element % PartIndex/=ParEnv % myPE ) THEN
           RETURN
        END IF
      END IF

!$OMP CRITICAL
       IF ( .NOT. ASSOCIATED( A % MassValues ) ) THEN
          ALLOCATE( A % MassValues(SIZE(A % Values)) )
          A % MassValues = 0.0d0
       END IF
!$OMP END CRITICAL

     ALLOCATE( M(DOFs*n,DOFs*n) )
     DO i=1,n*DOFs/2
       DO j=1,n*DOFs/2
         M(2*(i-1)+1, 2*(j-1)+1) =   REAL( MC(i,j) )
         M(2*(i-1)+1, 2*(j-1)+2) = -AIMAG( MC(i,j) )
         M(2*(i-1)+2, 2*(j-1)+1) =  AIMAG( MC(i,j) )
         M(2*(i-1)+2, 2*(j-1)+2) =   REAL( MC(i,j) )
       END DO
     END DO

     ! flip mass matrix for periodic elimination
     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalMatrix( Solver, n, Indexes, x % dofs, M )
     END IF                

     CALL UpdateMassMatrix( A, M, n, x % DOFs, x % Perm(Indexes(1:n)), &
                    A % MassValues )
     DEALLOCATE( M )
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdateMassC
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE DefaultUpdateDampR( B, UElement, USolver ) 
!------------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL,  TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     REAL(KIND=dp)   :: B(:,:)

     TYPE(Solver_t), POINTER   :: Solver
     TYPE(Matrix_t), POINTER   :: A
     TYPE(Variable_t), POINTER :: x
     TYPE(Element_t), POINTER  :: Element, P1, P2

     INTEGER :: i,j,n
     INTEGER, POINTER :: Indexes(:)

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
     ELSE
        Solver => CurrentModel % Solver
     END IF

     A => Solver % Matrix
     x => Solver % Variable

     Element => GetCurrentElement( UElement ) 

     Indexes => GetIndexStore()
     n =  GetElementDOFs( Indexes, Element, Solver )

     IF ( ParEnv % PEs > 1 ) THEN
       IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
          P1 => Element % BoundaryInfo % Left
          P2 => Element % BoundaryInfo % Right
          IF ( ASSOCIATED(P1) .AND. ASSOCIATED(P2) ) THEN
            IF ( P1 % PartIndex/=ParEnv % myPE .AND. &
                 P2 % PartIndex/=ParEnv % myPE )RETURN

            IF ( P1 % PartIndex/=ParEnv % myPE .OR. &
                 P2 % PartIndex/=ParEnv % myPE ) B=B/2;
          ELSE IF ( ASSOCIATED(P1) ) THEN
            IF ( P1 % PartIndex /= ParEnv % myPE ) RETURN
          ELSE IF ( ASSOCIATED(P2) ) THEN
            IF ( P2 % PartIndex /= ParEnv % myPE ) RETURN
          END IF
       ELSE IF ( Element % PartIndex/=ParEnv % myPE ) THEN
          RETURN
       END IF
     END IF

!$OMP CRITICAL
     IF ( .NOT. ASSOCIATED( A % DampValues ) ) THEN
       ALLOCATE( A % DampValues(SIZE(A % Values)) ) 
       A % DampValues = 0.0d0
     END IF
!$OMP END CRITICAL

     ! flip damp matrix for periodic elimination
     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalMatrix( Solver, n, Indexes, x % dofs, B )
     END IF                
     
     CALL UpdateMassMatrix( A, B, n, x % DOFs, x % Perm(Indexes(1:n)), &
              A  % DampValues )

     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalMatrix( Solver, n, Indexes, x % dofs, B )
     END IF                
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdateDampR
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE DefaultUpdateDampC( BC, UElement, USolver ) 
!------------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL,  TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     COMPLEX(KIND=dp)   :: BC(:,:)

     TYPE(Solver_t), POINTER   :: Solver
     TYPE(Matrix_t), POINTER   :: A
     TYPE(Variable_t), POINTER :: x
     TYPE(Element_t), POINTER  :: Element, P1, P2

     REAL(KIND=dp), ALLOCATABLE :: B(:,:)

     INTEGER :: i,j,n,DOFs
     INTEGER, POINTER :: Indexes(:)

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
     ELSE
        Solver => CurrentModel % Solver
     END IF

     A => Solver % Matrix
     x => Solver % Variable

     Element => GetCurrentElement( UElement ) 

     DOFs = x % DOFs
     Indexes => GetIndexStore()
     n =  GetElementDOFs( Indexes, Element, Solver )

     IF ( ParEnv % PEs > 1 ) THEN
       IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
          P1 => Element % BoundaryInfo % Left
          P2 => Element % BoundaryInfo % Right
          IF ( ASSOCIATED(P1) .AND. ASSOCIATED(P2) ) THEN
            IF ( P1 % PartIndex/=ParEnv % myPE .AND. &
                 P2 % PartIndex/=ParEnv % myPE )RETURN

            IF ( P1 % PartIndex/=ParEnv % myPE .OR. &
                 P2 % PartIndex/=ParEnv % myPE ) BC=BC/2
          ELSE IF ( ASSOCIATED(P1) ) THEN
            IF ( P1 % PartIndex/=ParEnv % myPE ) RETURN
          ELSE IF ( ASSOCIATED(P2) ) THEN
            IF ( P2 % PartIndex/=ParEnv % myPE ) RETURN
          END IF
       ELSE IF ( Element % PartIndex/=ParEnv % myPE ) THEN
          RETURN
       END IF
     END IF

!$OMP CRITICAL
       IF ( .NOT. ASSOCIATED( A % DampValues ) ) THEN
        ALLOCATE( A % DampValues(SIZE(A % Values)) ) 
        A % DampValues = 0.0d0
       END IF
!$OMP END CRITICAL

     ALLOCATE( B(DOFs*n, DOFs*n) )
     DO i=1,n*DOFs/2
       DO j=1,n*DOFs/2
         B(2*(i-1)+1, 2*(j-1)+1) =   REAL( BC(i,j) )
         B(2*(i-1)+1, 2*(j-1)+2) = -AIMAG( BC(i,j) )
         B(2*(i-1)+2, 2*(j-1)+1) =  AIMAG( BC(i,j) )
         B(2*(i-1)+2, 2*(j-1)+2) =   REAL( BC(i,j) )
       END DO
     END DO

     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalMatrix( Solver, n, Indexes, x % dofs, B )
     END IF                
     
     CALL UpdateMassMatrix( A, B, n, x % DOFs, x % Perm(Indexes(1:n)), &
                 A % DampValues )
     DEALLOCATE( B )
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdateDampC
!------------------------------------------------------------------------------


  
!------------------------------------------------------------------------------
  SUBROUTINE DefaultUpdateBulkR( B, F, UElement, USolver ) 
!------------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL,  TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     REAL(KIND=dp)   :: B(:,:), F(:)

     TYPE(Solver_t), POINTER   :: Solver
     TYPE(Matrix_t), POINTER   :: A
     TYPE(Variable_t), POINTER :: x
     TYPE(Element_t), POINTER  :: Element, P1, P2

     INTEGER :: i,j,n
     INTEGER, POINTER :: Indexes(:)

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
     ELSE
        Solver => CurrentModel % Solver
     END IF

     A => Solver % Matrix
     x => Solver % Variable

     Element => GetCurrentElement( UElement ) 

     Indexes => GetIndexStore()
     n =  GetElementDOFs( Indexes, Element, Solver )

     IF ( ParEnv % PEs > 1 ) THEN
       IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
          P1 => Element % BoundaryInfo % Left
          P2 => Element % BoundaryInfo % Right
          IF ( ASSOCIATED(P1) .AND. ASSOCIATED(P2) ) THEN
            IF ( P1 % PartIndex/=ParEnv % myPE .AND. &
                 P2 % PartIndex/=ParEnv % myPE )RETURN

            IF ( P1 % PartIndex/=ParEnv % myPE .OR. &
                 P2 % PartIndex/=ParEnv % myPE ) THEN
              B=B/2; F=F/2; 
            END IF
          ELSE IF ( ASSOCIATED(P1) ) THEN
            IF ( P1 % PartIndex /= ParEnv % myPE ) RETURN
          ELSE IF ( ASSOCIATED(P2) ) THEN
            IF ( P2 % PartIndex /= ParEnv % myPE ) RETURN
          END IF
       ELSE IF ( Element % PartIndex/=ParEnv % myPE ) THEN
          RETURN
       END IF
     END IF

!$OMP CRITICAL
     IF ( .NOT. ASSOCIATED( A % BulkValues ) ) THEN
       ALLOCATE( A % BulkValues(SIZE(A % Values)) ) 
       A % BulkValues = 0.0_dp
     END IF
!$OMP END CRITICAL

!$OMP CRITICAL
     IF ( .NOT. ASSOCIATED( A % BulkRHS ) ) THEN
       ALLOCATE( A % BulkRHS(SIZE(A % RHS)) ) 
       A % BulkRHS = 0.0_dp
     END IF
!$OMP END CRITICAL

       
     ! If we have any antiperiodic entries we need to check them all!
     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalMatrix( Solver, n, Indexes, x % dofs, B )
       CALL FlipPeriodicLocalForce( Solver, n, Indexes, x % dofs, f )
     END IF
     
     CALL UpdateGlobalEquations( A,B,A % BulkRHS,f,n,x % DOFs,x % Perm(Indexes(1:n)),  &
         GlobalValues=A % BulkValues )
     
     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalMatrix( Solver, n, Indexes, x % dofs, B )
       CALL FlipPeriodicLocalForce( Solver, n, Indexes, x % dofs, f )
     END IF
      
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdateBulkR
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE DefaultUpdateBulkC( BC, FC, UElement, USolver ) 
!------------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL,  TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     COMPLEX(KIND=dp)   :: BC(:,:), FC(:)

     TYPE(Solver_t), POINTER   :: Solver
     TYPE(Matrix_t), POINTER   :: A
     TYPE(Variable_t), POINTER :: x
     TYPE(Element_t), POINTER  :: Element, P1, P2

     REAL(KIND=dp), ALLOCATABLE :: B(:,:),F(:)

     INTEGER :: i,j,n,DOFs
     INTEGER, POINTER :: Indexes(:)

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
     ELSE
        Solver => CurrentModel % Solver
     END IF

     A => Solver % Matrix
     x => Solver % Variable

     Element => GetCurrentElement( UElement ) 

     DOFs = x % DOFs
     Indexes => GetIndexStore()
     n =  GetElementDOFs( Indexes, Element, Solver )

     IF ( ParEnv % PEs > 1 ) THEN
       IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
          P1 => Element % BoundaryInfo % Left
          P2 => Element % BoundaryInfo % Right
          IF ( ASSOCIATED(P1) .AND. ASSOCIATED(P2) ) THEN
            IF ( P1 % PartIndex/=ParEnv % myPE .AND. &
                 P2 % PartIndex/=ParEnv % myPE )RETURN

            IF ( P1 % PartIndex/=ParEnv % myPE .OR. &
                 P2 % PartIndex/=ParEnv % myPE ) THEN
              BC=BC/2; FC=FC/2;
            END IF
          ELSE IF ( ASSOCIATED(P1) ) THEN
            IF ( P1 % PartIndex/=ParEnv % myPE ) RETURN
          ELSE IF ( ASSOCIATED(P2) ) THEN
            IF ( P2 % PartIndex/=ParEnv % myPE ) RETURN
          END IF
       ELSE IF ( Element % PartIndex/=ParEnv % myPE ) THEN
          RETURN
       END IF
     END IF


!$OMP CRITICAL
     IF ( .NOT. ASSOCIATED( A % BulkValues ) ) THEN
       ALLOCATE( A % BulkValues(SIZE(A % Values)) ) 
       A % BulkValues = 0.0_dp
     END IF
!$OMP END CRITICAL

!$OMP CRITICAL
     IF ( .NOT. ASSOCIATED( A % BulkRHS ) ) THEN
       ALLOCATE( A % BulkRHS(SIZE(A % RHS)) ) 
       A % BulkRHS = 0.0_dp
     END IF
!$OMP END CRITICAL

     ALLOCATE( B(DOFs*n, DOFs*n), F(DOFs*n) )
     DO i=1,n*DOFs/2
       DO j=1,n*DOFs/2
         F( 2*(i-1)+1 ) =  REAL( FC(i) )
         F( 2*(i-1)+2 ) = AIMAG( FC(i) )

         B(2*(i-1)+1, 2*(j-1)+1) =   REAL( BC(i,j) )
         B(2*(i-1)+1, 2*(j-1)+2) = -AIMAG( BC(i,j) )
         B(2*(i-1)+2, 2*(j-1)+1) =  AIMAG( BC(i,j) )
         B(2*(i-1)+2, 2*(j-1)+2) =   REAL( BC(i,j) )
       END DO
     END DO

     IF( Solver % PeriodicFlipActive ) THEN
       CALL FlipPeriodicLocalMatrix( Solver, n, Indexes, x % dofs, B )
       CALL FlipPeriodicLocalForce( Solver, n, Indexes, x % dofs, f )
     END IF
      
     CALL UpdateGlobalEquations( A,B,A % BulkRHS,f,n,x % DOFs,x % Perm(Indexes(1:n)), &
                 GlobalValues=A % BulkValues )

     DEALLOCATE( B, F )
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdateBulkC
!------------------------------------------------------------------------------


  SUBROUTINE DefaultUpdateDirichlet( Dvals, UElement, USolver, UIndexes, Dset ) 
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Dvals(:)
    TYPE(Solver_t), OPTIONAL,TARGET :: USolver
    TYPE(Element_t), OPTIONAL, TARGET :: UElement
    INTEGER, OPTIONAL, TARGET :: UIndexes(:)
    LOGICAL, OPTIONAL :: Dset(:)
    
    TYPE(Solver_t), POINTER   :: Solver
    TYPE(Matrix_t), POINTER   :: A
    TYPE(Variable_t), POINTER :: x
    TYPE(Element_t), POINTER  :: Element
    
    INTEGER :: i,j,n
    INTEGER, POINTER :: Indexes(:)

    
    IF ( PRESENT( USolver ) ) THEN
      Solver => USolver
    ELSE
      Solver => CurrentModel % Solver
    END IF
    
    A => Solver % Matrix
    x => Solver % Variable
    
    Element => GetCurrentElement( UElement ) 

    IF( PRESENT( UIndexes ) ) THEN
      Indexes => Uindexes
      n = SIZE( Uindexes ) 
    ELSE
      Indexes => GetIndexStore()
      n =  GetElementDOFs( Indexes, Element, Solver )
    END IF
      
    DO i=1,n
      IF( PRESENT( Dset ) ) THEN
        IF( .NOT. Dset(i) ) CYCLE
      END IF
      CALL UpdateDirichletDof( A, Indexes(i), DVals(i) )
    END DO

  END SUBROUTINE DefaultUpdateDirichlet
  
       
 

!> Sets the Dirichlet conditions related to the variables of the active solver.
!------------------------------------------------------------------------------------------
  SUBROUTINE DefaultDirichletBCs( USolver,Ux,UOffset,OffDiagonalMatrix)
!------------------------------------------------------------------------------------------
     USE ElementDescription, ONLY: FaceElementOrientation
     IMPLICIT NONE

     INTEGER, OPTIONAL :: UOffset
     LOGICAL, OPTIONAL :: OffDiagonalMatrix
     TYPE(Variable_t), OPTIONAL, TARGET :: Ux
     TYPE(Solver_t), OPTIONAL, TARGET :: USolver
!--------------------------------------------------------------------------------------------     
     TYPE(Matrix_t), POINTER   :: A
     TYPE(Variable_t), POINTER :: x
     TYPE(Solver_t), POINTER :: Solver
     TYPE(ValueListEntry_t), POINTER :: ptr
     TYPE(ValueList_t), POINTER :: BC, Params
     TYPE(Element_t), POINTER :: Element, Parent, Edge, Face, SaveElement

     REAL(KIND=dp), ALLOCATABLE :: Work(:), STIFF(:,:)
     REAL(KIND=dp), POINTER :: b(:)
     REAL(KIND=dp), POINTER :: DiagScaling(:)
     REAL(KIND=dp) :: xx, s, dval, Cond
     REAL(KIND=dp) :: DefaultDOFs(4)

     INTEGER, ALLOCATABLE :: lInd(:), gInd(:)
     INTEGER :: FDofMap(6,4)
     INTEGER :: i, j, k, kk, l, m, n, nd, nb, np, mb, nn, ni, nj, i0
     INTEGER :: NDOFs, EDOFs, FDOFs, DOF, local, numEdgeDofs, istat, n_start, Offset
     INTEGER :: ActiveFaceId

     LOGICAL :: ReverseSign(6)
     LOGICAL :: Flag,Found, ConstantValue, ScaleSystem, DirichletComm
     LOGICAL :: PiolaTransform, QuadraticApproximation, SecondKindBasis
     LOGICAL, ALLOCATABLE :: ReleaseDir(:)
     LOGICAL :: ReleaseAny, NodalBCsWithBraces,AllConstrained
     LOGICAL :: CheckRight, AugmentedEigenSystem
     
     CHARACTER(:), ALLOCATABLE :: Name

     SAVE gInd, lInd, STIFF, Work
!-------------------------------------------------------------------------------------------- 

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
     ELSE
        Solver => CurrentModel % Solver
     END IF

     Params => GetSolverParams(Solver)

     IF ( GetString(Params,'Linear System Solver',Found)=='feti') THEN
       IF ( GetLogical(Params,'Total FETI', Found)) RETURN
     END IF

     A => Solver % Matrix
     b => A % RHS
     IF ( PRESENT(Ux) ) THEN
       x => Ux
     ELSE
       x => Solver % Variable
     END IF

     IF(.NOT.ALLOCATED(A % ConstrainedDOF)) THEN
       ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
       A % ConstrainedDOF = .FALSE.
     ELSE
       IF (SIZE(A % ConstrainedDOF) < A % NumberOfRows) THEN
         DEALLOCATE(A % ConstrainedDOF)
         ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
         A % ConstrainedDOF = .FALSE.
       END IF
     END IF
       
     IF(.NOT.ALLOCATED(A % Dvalues)) THEN
       ALLOCATE(A % Dvalues(A % NumberOfRows))
       A % Dvalues = 0._dp
     ELSE
       IF (SIZE(A % Dvalues) < A % NumberOfRows) THEN
         DEALLOCATE(A % Dvalues)
         ALLOCATE(A % Dvalues(A % NumberOfRows))
         A % Dvalues = 0._dp
       END IF
     END IF
     
     ! Create soft limiters to be later applied by the Dirichlet conditions
     ! This is done only once for each solver, hence the complex logic. 
     !---------------------------------------------------------------------
     IF( ListGetLogical( Params,'Apply Limiter',Found) ) THEN
       CALL DetermineSoftLimiter( Solver )	

       ! It is difficult to determine whether loads should be computed before or after setting the limiter.
       ! There are cases where both alternative are needed.
       IF(ListGetLogical( Params,'Apply Limiter Loads After',Found) ) THEN         
         DO DOF=1,x % DOFs
           name = TRIM(x % name)
           IF (x % DOFs>1) name=ComponentName(name,DOF)              
           CALL SetNodalLoads( CurrentModel,A,A % rhs, &
               Name,DOF,x % DOFs,x % Perm ) 
         END DO
       END IF
     END IF
         
     
     Offset = 0
     IF(PRESENT(UOffset)) Offset=UOffset

     n = Solver % Mesh % MaxElementDOFs
     IF ( .NOT. ALLOCATED( gInd ) ) THEN
       ALLOCATE( gInd(n), lInd(n), STIFF(n,n), Work(n), stat=istat )
       IF ( istat /= 0 ) &
          CALL Fatal('DefUtils::DefaultDirichletBCs','Memory allocation failed.' )
     ELSE IF ( SIZE(gInd) < n ) THEN
       DEALLOCATE( gInd, lInd, STIFF, Work )
       ALLOCATE( gInd(n), lInd(n), STIFF(n,n), Work(n), stat=istat )
       IF ( istat /= 0 ) &
          CALL Fatal('DefUtils::DefaultDirichletBCs','Memory allocation failed.' )
     END IF

     ReleaseAny = ListCheckPrefixAnyBC(CurrentModel, 'release '//TRIM(x % name)//' {e}')
     IF( ReleaseAny ) THEN
       CALL Info('DefaultDirichletBCs','Getting ready to block some Dirichlet BCs!',Level=7)
       ALLOCATE( ReleaseDir( SIZE( A % DValues ) ) )
       ReleaseDir = .FALSE.
     END IF

     NodalBCsWithBraces = .FALSE.
     NDOFs = MAXVAL(Solver % Def_Dofs(:,:,1))
     IF (NDOFs > 0) THEN
       DO DOF=1,x % DOFs
         name = TRIM(x % name)
         IF (x % DOFs > 1) name = ComponentName(name,DOF)
         NodalBCsWithBraces = ListCheckPrefixAnyBC(CurrentModel, Name//' {n}')
         IF (NodalBCsWithBraces) THEN
           CALL Info('DefaultDirichletBCs', '{n} construct is now used to set BCs', Level=7)
           CALL Info('DefaultDirichletBCs', I2S(NDOFs)//'-component {n} definition is handled', Level=7)
           EXIT
         END IF
       END DO
     END IF

     IF ( x % DOFs > 1 ) THEN
       CALL SetDirichletBoundaries( CurrentModel,A, b, GetVarName(x),-1,x % DOFs,x % Perm )
     END IF

     CALL Info('DefUtils::DefaultDirichletBCs', &
            'Setting Dirichlet boundary conditions', Level=10)
     

     ! ----------------------------------------------------------------------
     ! Perform some preparations if BCs for p-approximation will be handled: 
     ! ----------------------------------------------------------------------
     IF (.NOT. NodalBCsWithBraces) THEN
       DO DOF=1,x % DOFs
         name = TRIM(x % name)
         IF ( x % DOFs > 1 ) name = ComponentName(name,DOF)

         IF( .NOT. ListCheckPresentAnyBC( CurrentModel, name ) ) CYCLE

         SaveElement => GetCurrentElement() 
         DO i=1,Solver % Mesh % NumberOfBoundaryElements
           Element => GetBoundaryElement(i)
           IF ( .NOT. ActiveBoundaryElement() ) CYCLE

           ! Get parent element:
           ! -------------------
           Parent => Element % BoundaryInfo % Left
           IF ( .NOT. ASSOCIATED( Parent ) ) &
               Parent => Element % BoundaryInfo % Right

           IF ( .NOT. ASSOCIATED(Parent) )   CYCLE

           BC => GetBC()
           IF ( .NOT.ASSOCIATED(BC) ) CYCLE

           ptr => ListFind(BC, Name,Found )
           IF ( .NOT. ASSOCIATED(ptr) ) CYCLE

           Constantvalue = ( ptr % type /= LIST_TYPE_CONSTANT_SCALAR_PROC )


           IF ( isActivePElement(Parent,Solver)) THEN
             n = GetElementNOFNodes()
             ! Get indexes of boundary dofs:
             CALL mGetBoundaryIndexesFromParent( Solver % Mesh, Element, gInd, numEdgeDofs )
             DO k=n+1,numEdgeDofs
               nb = x % Perm( gInd(k) )
               IF ( nb <= 0 ) CYCLE
               nb = Offset + x % DOFs * (nb-1) + DOF

               IF ( ConstantValue  ) THEN
                 A % ConstrainedDOF(nb) = .TRUE.
                 A % Dvalues(nb) = 0._dp
               ELSE
                 CALL ZeroRow( A, nb )
                 A % RHS(nb) = 0._dp
               END IF
             END DO
           ELSE
             ! To do: Check whether BCs for edge/face elements must be set via L2 projection.
             CYCLE 
           END IF
         END DO

         SaveElement => SetCurrentElement(SaveElement)
       END DO
     END IF

     IF (NodalBCsWithBraces) THEN
       CALL Info('DefaultDirichletBCs','Setting nodal dofs with {n} construct', Level=15) 
       !
       ! NOTE:  This is still simplistic implementation and lacks many options which work
       !        in the case of standard BCs for nodal (Lagrange) interpolation approximations.
       ! TO DO: Consider how the functionality of the subroutines SetNodalLoads and SetDirichletBoundaries
       !        could be enabled in the case of {n} construct.
       !
       DO DOF=1,x % DOFs
         DO m=1,NDOFs
           name = TRIM(x % name)
           IF ( x % DOFs > 1 ) THEN
             name = ComponentName(name,DOF)//' {n}'
           ELSE
             name = name//' {n}'
           END IF

           ! When the component names are created for example from E[E Re:1 E Im:1], we now have name = "E Re {n}" or "E Im {n}".
           ! Finally, we append this by the field index, so that we may seek values for "E Re {n} m" and "E Im {n} m", where 
           ! m = 1,...,NDOFs when an element definition "n:NDOFs e:..." is given. 
           
           IF (NDOFs > 1) name = name//' '//I2S(m)

!           print *, '====== m is ', m
!           print *, '====== DOF is ', DOF
!           print *, 'operating name ', name

           SaveElement => GetCurrentElement()
           DO i=1,Solver % Mesh % NumberOfBoundaryElements
             Element => GetBoundaryElement(i)

             BC => GetBC()
             IF ( .NOT.ASSOCIATED(BC) ) CYCLE
             IF ( .NOT. ListCheckPresent(BC, TRIM(Name)) ) CYCLE

             Cond = SUM(GetReal(BC,GetVarName(Solver % Variable)//' Condition',Found))/n
             IF (Cond>0) CYCLE

             nd = GetElementDOFs(gInd, Element)
             n = Element % TYPE % NumberOfNodes

             Work(1:n) = GetReal(BC, Name, Found, Element)

!             print *, 'permutation size for single-field', nd
!             print *, 'element % ndofs, nofnodes ', element % ndofs, Element % TYPE % NumberOfNodes
!             print *, 'global indexes for single field', gind(1:element % ndofs)
!             print *, 'the field values ', Work(1:n)

             DO j=1,n
               k = (j-1) * NDOFs + m
               l = x % Perm(gInd(k))

               l = x % DOFs * (l-1) + DOF

               A % ConstrainedDOF(l) = .TRUE.
               A % Dvalues(l) = Work(j)
             END DO
           END DO
           SaveElement => SetCurrentElement(SaveElement)
         END DO
       END DO

     ELSE
       ! -------------------------------------------------------------------------------------
       ! Set BCs for fields which are approximated using H1-conforming basis functions 
       ! (either Lagrange basis or hierarchic p-basis): 
       ! -------------------------------------------------------------------------------------    
       DO DOF=1,x % DOFs
         name = TRIM(x % name)
         IF (x % DOFs>1) name=ComponentName(name,DOF)

         CALL SetDirichletBoundaries( CurrentModel, A, b, &
             Name, DOF, x % DOFs, x % Perm, Offset, OffDiagonalMatrix )

         ! ----------------------------------------------------------------------------
         ! Set Dirichlet BCs for edge and face dofs which come from approximating with
         ! p-elements:
         ! ----------------------------------------------------------------------------
         IF( .NOT. ListCheckPresentAnyBC( CurrentModel, name ) ) CYCLE

         CALL Info('DefUtils::DefaultDirichletBCs', &
             'p-element condition setup: '//name, Level=15)

         SaveElement => GetCurrentElement()
         DO i=1,Solver % Mesh % NumberOfBoundaryElements
           Element => GetBoundaryElement(i)
           IF ( .NOT. ActiveBoundaryElement() ) CYCLE

           BC => GetBC()
           IF ( .NOT.ASSOCIATED(BC) ) CYCLE
           IF ( .NOT. ListCheckPresent(BC, Name) ) CYCLE

           ! Get parent element:
           ! -------------------
           Parent => Element % BoundaryInfo % Left
           IF ( .NOT. ASSOCIATED( Parent ) ) THEN
             Parent => Element % BoundaryInfo % Right
           END IF
           IF ( .NOT. ASSOCIATED( Parent ) )   CYCLE

           ! Here set constraints for p-approximation only: 
           ! -----------------------------------------------------
           IF (.NOT.isActivePElement(Parent)) CYCLE

           ptr => ListFind(BC, Name,Found )
           Constantvalue = Ptr % Type /= LIST_TYPE_CONSTANT_SCALAR_PROC

           IF ( ConstantValue ) CYCLE

           SELECT CASE(Parent % Type % Dimension )
           CASE(2)
             ! If no edges do not try to set boundary conditions
             ! @todo This should changed to EXIT
             IF ( .NOT. ASSOCIATED( Solver % Mesh % Edges ) ) CYCLE

             ! If boundary edge has no dofs move on to next edge
             IF (Element % BDOFs <= 0) CYCLE

             ! Number of nodes for this element
             n = Element % TYPE % NumberOfNodes

             ! Get indexes for boundary and values for dofs associated to them
             CALL mGetBoundaryIndexesFromParent( Solver % Mesh, Element, gInd, nd)
             CALL LocalBcBDOFs( BC, Element, nd, Name, STIFF, Work )

             AllConstrained = .TRUE.
             DO l=1,n
               nb = x % Perm( gInd(l) )
               IF ( nb <= 0 ) CYCLE
               nb = Offset + x % DOFs * (nb-1) + DOF

               IF(A % ConstrainedDOF(nb)) THEN
                 s = A % Dvalues(nb)
                 DO k=n+1,nd
                   Work(k) = Work(k) - s*STIFF(k,l)
                   STIFF(k,l) = 0.0_dp
                 END DO
               ELSE
                 AllConstrained = .FALSE.
               END IF
             END DO

             IF(AllConstrained.AND.CoordinateSystemDimension()<=2) THEN
               IF(nd==n+1) THEN
                 Work(nd) = Work(nd) / STIFF(nd,nd)
               ELSE
                 CALL SolveLinSys(STIFF(n+1:nd,n+1:nd), Work(n+1:nd), nd-n)
               END IF

               DO l=n+1,nd
                 nb = x % Perm( gInd(l) )
                 IF ( nb <= 0 ) CYCLE
                 nb = Offset + x % DOFs * (nb-1) + DOF

                 A % Dvalues(nb) = Work(l)
                 A % ConstrainedDOF(nb) = .TRUE.
               END DO
             ELSE
               ! Contribute this boundary to global system
               ! (i.e solve global boundary problem)
               A % Symmetric = .FALSE.
               DO k=1,nd
                 nb = x % Perm( gInd(k) )
                 IF ( nb <= 0 ) CYCLE
                 nb = Offset + x % DOFs * (nb-1) + DOF
                 IF(.NOT.A % ConstrainedDOF(nb)) THEN
                   A % RHS(nb) = A % RHS(nb) + Work(k) 
                   DO l=n+1,nd
                     mb = x % Perm( gInd(l) )
                     IF ( mb <= 0 ) CYCLE
                     mb = Offset + x % DOFs * (mb-1) + DOF
                     DO kk=A % Rows(nb)+DOF-1,A % Rows(nb+1)-1,x % DOFs
                       IF ( A % Cols(kk) == mb ) THEN
                         A % Values(kk) = A % Values(kk) + STIFF(k,l)
                         EXIT
                       END IF
                     END DO
                   END DO
                 END IF
               END DO
             END IF

           CASE(3)
             ! If no faces present do not try to set boundary conditions
             ! @todo This should be changed to EXIT
             IF ( .NOT. ASSOCIATED(Solver % Mesh % Faces) ) CYCLE

             ! Parameters of element
             n = Element % TYPE % NumberOfNodes

             ! Get global boundary indexes and solve dofs associated to them
             CALL mGetBoundaryIndexesFromParent( Solver % Mesh, Element, gInd, nd )

             ! If boundary face has no dofs skip to next boundary element
             IF (nd == n) CYCLE

             ! Get local solution
             CALL LocalBcBDofs( BC, Element, nd, Name, STIFF, Work )

             DO l=1,n
               nb = x % Perm( gInd(l) )
               IF ( nb <= 0 ) CYCLE
               nb = Offset + x % DOFs * (nb-1) + DOF

               IF(A % ConstrainedDOF(nb)) THEN
                 s = A % Dvalues(nb)
                 DO k=n+1,nd
                   Work(k) = Work(k) - s*STIFF(k,l)
                   STIFF(k,l) = 0
                   STIFF(l,k) = 0
                 END DO
               END IF
             END DO

             ! Contribute this entry to global boundary problem
             A % Symmetric = .FALSE.
             DO k=1,nd
               nb = x % Perm( gInd(k) )
               IF ( nb <= 0 ) CYCLE
               nb = Offset + x % DOFs * (nb-1) + DOF

               IF(.NOT.A % ConstrainedDOF(nb)) THEN
                 A % RHS(nb) = A % RHS(nb) + Work(k) 
                 DO l=n+1,nd
                   mb = x % Perm( gInd(l) )
                   IF ( mb <= 0 ) CYCLE
                   mb = Offset + x % DOFs * (mb-1) + DOF
                   DO kk=A % Rows(nb)+DOF-1,A % Rows(nb+1)-1,x % DOFs
                     IF ( A % Cols(kk) == mb ) THEN
                       A % Values(kk) = A % Values(kk) + STIFF(k,l)
                       EXIT
                     END IF
                   END DO
                 END DO
               END IF
             END DO
           END SELECT
         END DO

         SaveElement => SetCurrentElement(SaveElement)
       END DO
     END IF

     ! BLOCK
     !------------------------------------
     DO DOF=1,x % DOFs
       IF(.NOT. ReleaseAny) CYCLE
       
       name = TRIM(x % name)
       IF (x % DOFs>1) name=ComponentName(name,DOF)

       IF ( .NOT. ListCheckPrefixAnyBC(CurrentModel, 'release '//TRIM(Name)//' {e}') ) CYCLE

       CALL Info('SetDefaultDirichlet','Release edge dofs from intersecting BCs',Level=15)

       SaveElement => GetCurrentElement()
       DO i=1,Solver % Mesh % NumberOfBoundaryElements
         Element => GetBoundaryElement(i)
         
         BC => GetBC()
         IF ( .NOT.ASSOCIATED(BC) ) CYCLE
         IF ( .NOT. ListGetLogical(BC, 'release '//TRIM(Name)//' {e}', Found ) ) CYCLE
         
         ! Get parent element:
         ! -------------------
         Parent => Element % BoundaryInfo % Left
         IF ( .NOT. ASSOCIATED( Parent ) ) THEN
           Parent => Element % BoundaryInfo % Right
         END IF
         IF ( .NOT. ASSOCIATED( Parent ) )   CYCLE
         np = Parent % TYPE % NumberOfNodes
         
         IF(.NOT. ASSOCIATED( Solver % Mesh % Edges ) ) CYCLE
         SELECT CASE(GetElementFamily())
             
         CASE(3,4)
           CALL PickActiveFace(Solver % Mesh, Parent, Element, Face, j)

           IF (.NOT. ASSOCIATED(Face)) CYCLE
           IF ( .NOT. ActiveBoundaryElement(Face) ) CYCLE

           DO l=1,Face % TYPE % NumberOfEdges
             Edge => Solver % Mesh % Edges(Face % EdgeIndexes(l))
             EDOFs = Edge % BDOFs
             IF (EDOFs == 0) CYCLE

             n = GetElementDOFs(gInd,Edge)

             IF (Solver % Def_Dofs(2,Parent % BodyId,1) > 0) THEN
               n_start = Edge % NDOFs
             ELSE
               n_start = 0
             END IF

             DO j=1,EDOFs
               k = n_start + j
               nb = x % Perm(gInd(k))
               IF ( nb <= 0 ) CYCLE
               nb = Offset + x % DOFs*(nb-1) + DOF
               ReleaseDir(nb) = .TRUE.
             END DO
           END DO
         END SELECT
       END DO
       SaveElement => SetCurrentElement(SaveElement)
     END DO
     
     
     !
     ! Apply special couple loads for 3-D models of solids:
     !
     CALL SetCoupleLoads(CurrentModel, x % Perm, A, b, x % DOFs )

     ! ----------------------------------------------------------------------------
     ! Set Dirichlet BCs for edge and face dofs which arise from approximating with
     ! edge (curl-conforming) or face (div-conforming) elements:
     ! ----------------------------------------------------------------------------
     QuadraticApproximation = ListGetLogical(Params, 'Quadratic Approximation', Found)
     SecondKindBasis = ListGetLogical(Params, 'Second Kind Basis', Found)
     DO DOF=1,x % DOFs
        name = TRIM(x % name)
        IF (x % DOFs>1) name=ComponentName(name,DOF)
        
        IF ( .NOT. ListCheckPrefixAnyBC(CurrentModel, Name//' {e}') .AND. &
             .NOT. ListCheckPrefixAnyBC(CurrentModel, Name//' {f}') ) CYCLE

        CALL Info('SetDefaultDirichlet','Setting edge and face dofs',Level=15)

        SaveElement => GetCurrentElement()
        DO i=1,Solver % Mesh % NumberOfBoundaryElements
           Element => GetBoundaryElement(i)

           BC => GetBC()
           IF ( .NOT.ASSOCIATED(BC) ) CYCLE
           IF ( .NOT. ListCheckPrefix(BC, Name//' {e}') .AND. &
                .NOT. ListCheckPrefix(BC, Name//' {f}') ) CYCLE

           Cond = SUM(GetReal(BC,GetVarName(Solver % Variable)//' Condition',Found))/n
           IF(Cond>0) CYCLE

           ! Get parent element:
           ! -------------------
           Parent => Element % BoundaryInfo % Left
           IF ( ASSOCIATED( Parent ) ) THEN
             IF (Parent % BodyId < 1) THEN
               CheckRight = .TRUE.
             ELSE
               CheckRight = .FALSE.
             END IF
           ELSE
             CheckRight = .TRUE.
           END IF
             
           IF (CheckRight) THEN
             Parent => Element % BoundaryInfo % Right
             IF ( ASSOCIATED( Parent ) ) THEN
               IF (Parent % BodyId < 1) THEN
                 Call Warn('SetDefaultDirichlet', 'Cannot set a BC owing to a missing parent body index')
                 CYCLE
               END IF
             END IF
           END IF
           IF ( .NOT. ASSOCIATED( Parent ) ) CYCLE
           np = Parent % TYPE % NumberOfNodes

           IF ( ListCheckPrefix(BC, Name//' {e}') ) THEN
              !--------------------------------------------------------------------------------
              ! We now devote this branch for handling edge (curl-conforming) finite elements 
              ! which, in addition to edge DOFs, may also have DOFs associated with faces. 
              !--------------------------------------------------------------------------------
              IF ( ASSOCIATED( Solver % Mesh % Edges ) ) THEN
                 SELECT CASE(GetElementFamily())
                 CASE(2)

                   CALL PickActiveFace(Solver % Mesh, Parent, Element, Edge, j)

                   IF ( .NOT. ASSOCIATED(Edge) ) CYCLE
                   IF ( .NOT. ActiveBoundaryElement(Edge) ) CYCLE                  

                   EDOFs = Edge % BDOFs     ! The number of DOFs associated with edges
                   IF (EDOFs < 1) CYCLE

                   AugmentedEigenSystem = ListGetLogical(Params, 'Eigen System Augmentation', Found) 
                   IF (AugmentedEigenSystem) THEN
                     EDOFs = EDOFs/2
                   END IF

                   n = Edge % TYPE % NumberOfNodes
                   CALL VectorElementEdgeDOFs(BC,Edge,n,Parent,np,Name//' {e}',Work, &
                       EDOFs, SecondKindBasis)

                   n=GetElementDOFs(gInd,Edge)

                   IF (Solver % Def_Dofs(2,Parent % BodyId,1) > 0) THEN
                     n_start = Edge % NDOFs
                   ELSE
                     n_start = 0
                   END IF

                   DO j=1,EDOFs
                     IF (AugmentedEigenSystem) THEN
                       k = n_start + 2*j - 1
                     ELSE
                       k = n_start + j
                     END IF
                     nb = x % Perm(gInd(k))
                     IF ( nb <= 0 ) CYCLE
                     nb = Offset + x % DOFs*(nb-1) + DOF

                     A % ConstrainedDOF(nb) = .TRUE.
                     A % Dvalues(nb) = Work(j) 
                   END DO

                 CASE(3,4)
                   CALL PickActiveFace(Solver % Mesh, Parent, Element, Face, j)

                   IF (.NOT. ASSOCIATED(Face)) CYCLE
                   IF ( .NOT. ActiveBoundaryElement(Face) ) CYCLE

                   ! ---------------------------------------------------------------------
                   ! Set first constraints for DOFs associated with edges. Save the values
                   ! of DOFs in the array Work(:), so that the possible remaining DOFs
                   ! associated with the face can be computed after this.
                   ! ---------------------------------------------------------------------
                   i0 = 0
                   DO l=1,Face % TYPE % NumberOfEdges
                     Edge => Solver % Mesh % Edges(Face % EdgeIndexes(l))
                     EDOFs = Edge % BDOFs
                     IF (EDOFs < 1) CYCLE                     
                     n = Edge % TYPE % NumberOfNodes

                     CALL VectorElementEdgeDOFs(BC, Edge, n, Parent, np, Name//' {e}', &
                         Work(i0+1:i0+EDOFs), EDOFs, SecondKindBasis)
                     
                     n = GetElementDOFs(gInd,Edge)

                     IF (Solver % Def_Dofs(2,Parent % BodyId,1) > 0) THEN
                       n_start = Edge % NDOFs
                     ELSE
                       n_start = 0
                     END IF
 
                     DO j=1,EDOFs
                       k = n_start + j
                       nb = x % Perm(gInd(k))
                       IF ( nb <= 0 ) CYCLE
                       nb = Offset + x % DOFs*(nb-1) + DOF

                       A % ConstrainedDOF(nb) = .TRUE.
                       A % Dvalues(nb) = Work(i0+j) 
                     END DO
                     i0 = i0 + EDOFs
                   END DO

                   ! ---------------------------------------------------------------------
                   ! Set constraints for face DOFs via seeking the best approximation in L2.
                   ! We use the variational equation (u x n,v') = (g x n - u0 x n,v) where
                   ! u0 denotes the part of the interpolating function u+u0 which is already 
                   ! known and v is a test function for the Galerkin method.
                   ! ---------------------------------------------------------------------
                   IF (Face % BDOFs > 0) THEN
                     EDOFs = i0 ! The count of edge DOFs set so far
                     n = Face % TYPE % NumberOfNodes

                     CALL SolveLocalFaceDOFs(BC, Face, n, Name//' {e}', Work, EDOFs, &
                         Face % BDOFs, QuadraticApproximation)

                     n = GetElementDOFs(GInd,Face)
                     DO j=1,Face % BDOFs
                       nb = x % Perm(GInd(n-Face % BDOFs+j)) ! The last entries should be face-DOF indices
                       IF ( nb <= 0 ) CYCLE
                       nb = Offset + x % DOFs*(nb-1) + DOF

                       A % ConstrainedDOF(nb) = .TRUE.
                       A % Dvalues(nb) = Work(EDOFs+j) 
                     END DO
                   END IF

                 END SELECT
              END IF
           ELSE IF ( ListCheckPrefix(BC, Name//' {f}') ) THEN
             !--------------------------------------------------------------------------
             ! This branch should be able to handle BCs for face (div-conforming)
             ! elements. Now this works only for RT(0), ABF(0) and BMD(1) in 2D and
             ! for a 48-DOF brick and the Nedelec tetrahedron of the first and second kind 
             ! in 3D.
             !--------------------------------------------------------------------------
             SELECT CASE(GetElementFamily())
             CASE(2)
               
               CALL PickActiveFace(Solver % Mesh, Parent, Element, Edge, j)

               IF (.NOT. ASSOCIATED(Edge)) CYCLE
               IF ( .NOT. ActiveBoundaryElement(Edge) ) CYCLE                  

               EDOFs = Edge % BDOFs     ! The number of DOFs associated with edges

               IF (EDOFs < 1) CYCLE

               n = Edge % TYPE % NumberOfNodes
               CALL VectorElementEdgeDOFs(BC,Edge,n,Parent,np,Name//' {f}',Work, &
                   EDOFs, SecondKindBasis, FaceElement=.TRUE.)

               n=GetElementDOFs(gInd,Edge)

               IF (Solver % Def_Dofs(2,Parent % BodyId,1) > 0) THEN
                 n_start = Edge % NDOFs
               ELSE
                 n_start = 0
               END IF
               DO j=1,EDOFs
                 k = n_start + j
                 nb = x % Perm(gInd(k))
                 IF ( nb <= 0 ) CYCLE
                 nb = Offset + x % DOFs*(nb-1) + DOF

                 A % ConstrainedDOF(nb) = .TRUE.
                 A % Dvalues(nb) = Work(j) 
               END DO

             CASE(3)
               
               CALL PickActiveFace(Solver % Mesh, Parent, Element, Face, ActiveFaceId)

               IF (.NOT. ASSOCIATED(Face)) CYCLE
               IF ( .NOT. ActiveBoundaryElement(Face) ) CYCLE

               FDOFs = Face % BDOFs

               IF (FDOFs > 0) THEN
                 CALL FaceElementOrientation(Parent, ReverseSign, ActiveFaceId)
                 IF (SecondKindBasis) &
                     CALL FaceElementBasisOrdering(Parent, FDofMap, ActiveFaceId)
                 n = Face % TYPE % NumberOfNodes

                 CALL FaceElementDOFs(BC, Face, n, Parent, ActiveFaceId, &
                     Name//' {f}', Work, FDOFs, SecondKindBasis)

                 IF (SecondKindBasis) THEN
                   !
                   ! Conform to the orientation and ordering used in the
                   ! assembly of the global equations
                   !
                   DefaultDOFs(1:FDOFs) = Work(1:FDOFs)
                   IF (ReverseSign(ActiveFaceId)) THEN
                     S = -1.0d0
                   ELSE
                     S = 1.0d0
                   END IF

                   DO j=1,FDOFs
                     k = FDofMap(ActiveFaceId,j)
                     Work(j) = S * DefaultDOFs(k)
                   END DO
                 ELSE
                   IF (ReverseSign(ActiveFaceId)) Work(1:FDOFs) = -1.0d0*Work(1:FDOFs)
                 END IF

                 n = GetElementDOFs(GInd,Face)
                 !
                 ! Make an offset by the count of nodal DOFs. This provides
                 ! the right starting point if edge DOFs are not present.
                 !
                 IF (Solver % Def_Dofs(3,Parent % BodyId,1) > 0) THEN
                   n_start = Face % NDOFs
                 ELSE
                   n_start = 0
                 END IF
                 !
                 ! Check if we need to increase the offset by the count of
                 ! edge DOFs:
                 !
                 IF ( ASSOCIATED(Face % EdgeIndexes) .AND. &
                     Solver % Def_Dofs(3,Parent % BodyId,2) > 0) THEN
                   EDOFs = 0
                   DO l=1,Face % TYPE % NumberOfEdges
                     Edge => Solver % Mesh % Edges(Face % EdgeIndexes(l))
                     EDOFs = EDOFs + Edge % BDOFs
                   END DO
                   n_start = n_start + EDOFs
                 END IF

                 DO j=1,FDOFs
                   k = n_start + j
                   nb = x % Perm(gInd(k))
                   IF ( nb <= 0 ) CYCLE
                   nb = Offset + x % DOFs*(nb-1) + DOF

                   A % ConstrainedDOF(nb) = .TRUE.
                   A % Dvalues(nb) = Work(j) 
                 END DO
               END IF

             CASE(4)
               
               CALL PickActiveFace(Solver % Mesh, Parent, Element, Face, ActiveFaceId)

               IF (.NOT. ASSOCIATED(Face)) CYCLE
               IF ( .NOT. ActiveBoundaryElement(Face) ) CYCLE

               FDOFs = Face % BDOFs

               IF (FDOFs > 0) THEN
                 CALL FaceElementBasisOrdering(Parent, FDofMap, ActiveFaceId, ReverseSign)
                 n = Face % TYPE % NumberOfNodes

                 CALL FaceElementDOFs(BC, Face, n, Parent, ActiveFaceId, &
                     Name//' {f}', Work, FDOFs)

                 !
                 ! Conform to the orientation and ordering used in the
                 ! assembly of the global equations
                 !
                 DefaultDOFs(1:FDOFs) = Work(1:FDOFs)
                 IF (ReverseSign(ActiveFaceId)) THEN
                   S = -1.0d0
                 ELSE
                   S = 1.0d0
                 END IF

                 DO j=1,FDOFs
                   k = FDofMap(ActiveFaceId,j)
                   Work(j) = S * DefaultDOFs(k)
                 END DO

                 n = GetElementDOFs(GInd,Face)

                 !
                 ! Make an offset by the count of nodal DOFs. This provides
                 ! the right starting point if edge DOFs are not present.
                 !
                 IF (Solver % Def_Dofs(4,Parent % BodyId,1) > 0) THEN
                   n_start = Face % NDOFs
                 ELSE
                   n_start = 0
                 END IF
                 !
                 ! Check if we need to increase the offset by the count of
                 ! edge DOFs:
                 !
                 IF ( ASSOCIATED(Face % EdgeIndexes) .AND. &
                     Solver % Def_Dofs(4,Parent % BodyId,2) > 0) THEN
                   EDOFs = 0
                   DO l=1,Face % TYPE % NumberOfEdges
                     Edge => Solver % Mesh % Edges(Face % EdgeIndexes(l))
                     EDOFs = EDOFs + Edge % BDOFs
                   END DO
                   n_start = n_start + EDOFs
                 END IF

                 DO j=1,FDOFs
                   k = n_start + j
                   nb = x % Perm(gInd(k))
                   IF ( nb <= 0 ) CYCLE
                   nb = Offset + x % DOFs*(nb-1) + DOF

                   A % ConstrainedDOF(nb) = .TRUE.
                   A % Dvalues(nb) = Work(j) 
                 END DO
               END IF

             CASE DEFAULT
               CALL Warn('DefaultDirichletBCs', 'Cannot set face element DOFs for this element shape')
             END SELECT
           END IF
         END DO
         SaveElement => SetCurrentElement(SaveElement)
     END DO

     IF( ReleaseAny) THEN
       IF( InfoActive(10) ) THEN
         k = COUNT( A % ConstrainedDOF ) 
         PRINT *,'Original number of of Dirichlet BCs:',k       
         k = COUNT( ReleaseDir )
         PRINT *,'Marked number of Dirichlet BCs not to set:',k       
         k = COUNT( ReleaseDir .AND. A % ConstrainedDOF )
         PRINT *,'Ignoring number of Dirichlet BCs:',k
       END IF         
       WHERE( ReleaseDir )
         A % ConstrainedDOF = .FALSE.
       END WHERE
     END IF
                
     ! Add the possible constraint modes structures
     !----------------------------------------------------------
     IF ( GetLogical(Params,'Constraint Modes Analysis',Found) ) THEN
       CALL SetConstraintModesBoundaries( CurrentModel, Solver, A, b, x % Name, x % DOFs, x % Perm )
     END IF
      
     
#ifdef HAVE_FETI4I
     IF(C_ASSOCIATED(A % PermonMatrix)) THEN
       CALL Info('DefUtils::DefaultDirichletBCs','Permon matrix, Dirichlet conditions registered but not set!', Level=5)
       RETURN
     END IF
#endif

     ! This is set outside so that it can be called more flexibilly
     CALL EnforceDirichletConditions( Solver, A, b )
     
 
     CALL Info('DefUtils::DefaultDirichletBCs','Dirichlet boundary conditions set', Level=12)
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultDirichletBCs
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> This subroutine computes the values of DOFs that are associated with 
!> mesh edges in the case of vector-valued (edge or face) finite elements, so that
!> the vector-valued interpolant of the BC data can be constructed. 
!> The values of the DOFs are defined as D = S*(g.e,v)_E where the unit vector e
!> can be either tangential or normal to the edge, g is vector-valued data, 
!> v is a polynomial on the edge E, and S reverses sign if necessary.
!------------------------------------------------------------------------------
  SUBROUTINE VectorElementEdgeDOFs(BC, Element, n, Parent, np, Name, Integral, EDOFs, &
      SecondFamily, FaceElement)
!------------------------------------------------------------------------------
    USE ElementDescription, ONLY: GetEdgeMap
    IMPLICIT NONE

    TYPE(ValueList_t), POINTER :: BC  !< The list of boundary condition values
    TYPE(Element_t), POINTER :: Element !< The boundary element handled
    INTEGER :: n                      !< The number of boundary element nodes
    TYPE(Element_t) :: Parent         !< The parent element of the boundary element
    INTEGER :: np                     !< The number of parent element nodes
    CHARACTER(LEN=*) :: Name          !< The name of boundary condition
    REAL(KIND=dp) :: Integral(:)      !< The values of DOFs
    INTEGER, OPTIONAL :: EDOFs        !< The number of DOFs
    LOGICAL, OPTIONAL :: SecondFamily !< To select the element family
    LOGICAL, OPTIONAL :: FaceElement  !< If .TRUE., e is normal to the edge
!------------------------------------------------------------------------------
    TYPE(Nodes_t), SAVE :: Nodes, Pnodes
    TYPE(ElementType_t), POINTER :: SavedType
    TYPE(GaussIntegrationPoints_t) :: IP

    LOGICAL :: Lstat, ReverseSign, SecondKindBasis, DivConforming
    INTEGER, POINTER :: Edgemap(:,:)
    INTEGER :: i,j,k,p,DOFs
    INTEGER :: i1,i2,i3

    REAL(KIND=dp) :: Basis(n),Load(n),Vload(3,n),VL(3),e(3),d(3)
    REAL(KIND=dp) :: E21(3),E32(3) 
    REAL(KIND=dp) :: u,v,L,s,DetJ
!------------------------------------------------------------------------------
    DOFs = 1
    IF (PRESENT(EDOFs)) THEN
      IF (EDOFs > 2) THEN
        CALL Fatal('VectorElementEdgeDOFs','Cannot handle more than 2 DOFs per edge')
      ELSE
        DOFs = EDOFs
      END IF
    END IF   

    IF (PRESENT(SecondFamily)) THEN
      SecondKindBasis = SecondFamily
      IF (SecondKindBasis .AND. (DOFs /= 2) ) &
          CALL Fatal('VectorElementEdgeDOFs','2 DOFs per edge expected')
    ELSE
      SecondKindBasis = .FALSE.
    END IF

    IF (PRESENT(FaceElement)) THEN
      DivConforming = FaceElement
    ELSE
      DivConforming = .FALSE.
    END IF

    ! Get the nodes of the boundary and parent elements:
    !CALL GetElementNodes(Nodes, Element)
    !CALL GetElementNodes(PNodes, Parent)
    !Remove references to DefUtils
    CALL CopyElementNodesFromMesh(Nodes, CurrentModel % Solver % Mesh, &
        n, Element % NodeIndexes)
    CALL CopyElementNodesFromMesh(PNodes, CurrentModel % Solver % Mesh, &
        np, Parent % NodeIndexes)
    
    
    ReverseSign = .FALSE.
    EdgeMap => GetEdgeMap(GetElementFamily(Parent))
    DO i=1,SIZE(EdgeMap,1)
      j=EdgeMap(i,1)
      k=EdgeMap(i,2)
      IF ( Parent % NodeIndexes(j)==Element % NodeIndexes(1) .AND. &
          Parent % NodeIndexes(k)==Element % NodeIndexes(2) ) THEN
        EXIT
      ELSE IF (Parent % NodeIndexes(j)==Element % NodeIndexes(2) .AND. &
          Parent % NodeIndexes(k)==Element % NodeIndexes(1) ) THEN
        ! This is the right edge but has opposite orientation as compared
        ! with the listing of the parent element edges
        ReverseSign = .TRUE.
        EXIT
      END IF
    END DO

    Load(1:n) = GetReal( BC, Name, Lstat, Element )

    i = LEN_TRIM(Name)
    VLoad(1,1:n) = GetReal(BC,Name(1:i)//' 1',Lstat,element)
    VLoad(2,1:n) = GetReal(BC,Name(1:i)//' 2',Lstat,element)
    VLoad(3,1:n) = GetReal(BC,Name(1:i)//' 3',Lstat,element)

    e(1) = PNodes % x(k) - PNodes % x(j)
    e(2) = PNodes % y(k) - PNodes % y(j)
    e(3) = PNodes % z(k) - PNodes % z(j)
    e = e/SQRT(SUM(e**2))
    IF (DivConforming) THEN
      ! The boundary normal is needed instead of the tangent vector.
      ! First, find the element director d that makes the parent 
      ! element an oriented surface. 
      i1 = EdgeMap(1,1)
      i2 = EdgeMap(1,2)
      i3 = EdgeMap(2,2)
      E21(1) = PNodes % x(i2) - PNodes % x(i1)
      E21(2) = PNodes % y(i2) - PNodes % y(i1)
      E21(3) = PNodes % z(i2) - PNodes % z(i1)
      E32(1) = PNodes % x(i3) - PNodes % x(i2)
      E32(2) = PNodes % y(i3) - PNodes % y(i2)
      E32(3) = PNodes % z(i3) - PNodes % z(i2)
      d = CrossProduct(E21, E32)
      d = d/SQRT(SUM(d**2))
      ! Set e to be the outward normal to the parent element: 
      e = CrossProduct(e, d)
    END IF

    ! Is this element type stuff needed and for what?
    SavedType => Element % TYPE
    IF ( GetElementFamily()==1 ) Element % TYPE=>GetElementType(202)
      
    Integral = 0._dp
    IP = GaussPoints(Element)
    DO p=1,IP % n
      Lstat = ElementInfo( Element, Nodes, IP % u(p), &
            IP % v(p), IP % w(p), DetJ, Basis )
      s = IP % s(p) * DetJ

      L  = SUM(Load(1:n)*Basis(1:n))
      VL = MATMUL(Vload(:,1:n),Basis(1:n))

      IF (SecondKindBasis) THEN
        u = IP % u(p)
        v = 0.5d0*(1.0d0-sqrt(3.0d0)*u)
        Integral(1)=Integral(1)+s*(L+SUM(VL*e))*v
        v = 0.5d0*(1.0d0+sqrt(3.0d0)*u)
        Integral(2)=Integral(2)+s*(L+SUM(VL*e))*v
      ELSE
        Integral(1)=Integral(1)+s*(L+SUM(VL*e))

        IF (.NOT. DivConforming) THEN
          ! This branch is concerned with the second-order curl-conforming elements
          IF (DOFs>1) THEN
            v = Basis(2)-Basis(1)
            ! The parent element must define the default for the positive tangent associated
            ! with the edge. Thus, if the boundary element handled has an opposite orientation, 
            ! the sign must be reversed to get the positive coordinate associated with the
            ! parent element edge.
            IF (ReverseSign) v = -1.0d0*v
            Integral(2)=Integral(2)+s*(L+SUM(VL*e))*v
          END IF
        END IF
      END IF
    END DO
    Element % TYPE => SavedType

    j = Parent % NodeIndexes(j)
    IF ( ParEnv % PEs>1 ) &
      j=CurrentModel % Mesh % ParallelInfo % GlobalDOFs(j)

    k = Parent % NodeIndexes(k)
    IF ( ParEnv % PEs>1 ) &
      k=CurrentModel % Mesh % ParallelInfo % GlobalDOFs(k)

    IF (k < j) THEN
      IF (SecondKindBasis) THEN
        Integral(1)=-Integral(1)
        Integral(2)=-Integral(2)
      ELSE
        Integral(1)=-Integral(1)
      END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE VectorElementEdgeDOFs
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> This subroutine computes the values of DOFs that are associated with 
!> mesh faces in the case of curl-conforming (edge) finite elements, so that
!> the edge finite element interpolant of the BC data can be constructed. 
!> The values of the DOFs are obtained as the best approximation in L2 when
!> the values of the DOFs associated with edges are given.
!------------------------------------------------------------------------------
  SUBROUTINE SolveLocalFaceDOFs(BC, Element, n, Name, DOFValues, &
      EDOFs, FDOFs, QuadraticApproximation)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(ValueList_t), POINTER :: BC     !< The list of boundary condition values
    TYPE(Element_t), POINTER :: Element  !< The boundary element handled
    INTEGER :: n                         !< The number of boundary element nodes
    CHARACTER(LEN=*) :: Name             !< The name of boundary condition
    REAL(KIND=dp) :: DOFValues(:)        !< The values of DOFs
    INTEGER :: EDOFs                     !< The number of edge DOFs
    INTEGER :: FDOFs                     !< The number of face DOFs
    LOGICAL :: QuadraticApproximation    !< Use second-order edge element basis
!------------------------------------------------------------------------------
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP

    LOGICAL :: Lstat

    INTEGER :: i,j,p,DOFs,BasisDegree

    REAL(KIND=dp) :: Basis(n),Vload(3,n),VL(3),Normal(3)
    REAL(KIND=dp) :: EdgeBasis(EDOFs+FDOFs,3)
    REAL(KIND=dp) :: Mass(FDOFs,FDOFs), Force(FDOFs)
    REAL(KIND=dp) :: v,s,DetJ
!------------------------------------------------------------------------------
    IF (QuadraticApproximation) THEN
      BasisDegree = 2
    ELSE
      BasisDegree = 1
    END IF
      
    Mass = 0.0d0
    Force = 0.0d0

    ! Remove dependencies to DefUtils
    CALL CopyElementNodesFromMesh(Nodes, CurrentModel % Solver % Mesh, &
        n, Element % NodeIndexes)
    !CALL GetElementNodes(Nodes, Element)

    i = LEN_TRIM(Name)
    VLoad(1,1:n)=GetReal(BC,Name(1:i)//' 1',Lstat,element)
    VLoad(2,1:n)=GetReal(BC,Name(1:i)//' 2',Lstat,element)
    VLoad(3,1:n)=GetReal(BC,Name(1:i)//' 3',Lstat,element)

    IP = GaussPoints(Element)
    DO p=1,IP % n

      Lstat = EdgeElementInfo( Element, Nodes, IP % u(p), IP % v(p), IP % w(p), &
          DetF=DetJ, Basis=Basis, EdgeBasis=EdgeBasis, BasisDegree=BasisDegree, &
          ApplyPiolaTransform=.TRUE., TangentialTrMapping=.TRUE.)

      Normal = NormalVector(Element, Nodes, IP % u(p), IP % v(p), .FALSE.)

      VL = MATMUL(Vload(:,1:n),Basis(1:n))

      s = IP % s(p) * DetJ

      DO i=1,FDOFs
        DO j=1,FDOFs
          Mass(i,j) = Mass(i,j) + SUM(EdgeBasis(EDOFs+i,:) * EdgeBasis(EDOFs+j,:)) * s
        END DO
        Force(i) = Force(i) + SUM(CrossProduct(VL,Normal) * EdgeBasis(EDOFs+i,:)) * s
        DO j=1,EDOFs
          Force(i) = Force(i) - DOFValues(j) * SUM(EdgeBasis(j,:) * EdgeBasis(EDOFs+i,:)) * s
        END DO
      END DO
    END DO

    CALL LUSolve(FDOFs, Mass(1:FDOFs,1:FDOFs), Force(1:FDOFs))
    DOFValues(EDOFs+1:EDOFs+FDOFs) = Force(1:FDOFs)
!------------------------------------------------------------------------------
  END SUBROUTINE SolveLocalFaceDOFs
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> This subroutine computes the values of DOFs that are associated with 
!> mesh faces in the case of face (vector-valued) finite elements, so that
!> the vector-valued interpolant of the BC data can be constructed. 
!> The values of the DOFs are defined as D = S*(g.n,v)_F where the unit vector n
!> is normal to the face, g is vector-valued data, v is a polynomial on the face F, 
!> and S reverses sign if necessary. This subroutine performs neither sign 
!> reversions nor the permutations of DOFs, i.e. the DOFs are returned in
!> the default form.
! TO DO: This may need an update when new 3-D elements are added. 
!------------------------------------------------------------------------------
  SUBROUTINE FaceElementDOFs(BC, Element, n, Parent, FaceId, Name, Integral, &
      FDOFs, SecondFamily)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(ValueList_t), POINTER, INTENT(IN) :: BC    !< The list of boundary condition values
    TYPE(Element_t), POINTER, INTENT(IN) :: Element !< The boundary element handled
    INTEGER, INTENT(IN) :: n                        !< The number of boundary element nodes
    TYPE(Element_t), POINTER, INTENT(IN) :: Parent  !< The parent element of the boundary element
    INTEGER, INTENT(IN) :: FaceId                 !< The parent element face corresponding to Element
    CHARACTER(LEN=*), INTENT(IN) :: Name          !< The variable name in the boundary condition
    REAL(KIND=dp), INTENT(OUT) :: Integral(:)     !< The values of DOFs
    INTEGER, OPTIONAL, INTENT(IN) :: FDOFs        !< The number of DOFs
    LOGICAL, OPTIONAL, INTENT(IN) :: SecondFamily !< To select the element family
!------------------------------------------------------------------------------
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(Element_t), POINTER :: ElementCopy
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: SecondKindBasis, stat, ElementCopyCreated
    INTEGER :: TetraFaceMap(4,3), BrickFaceMap(6,4), ActiveFaceMap(4)
    INTEGER :: DOFs, i, j, p
    REAL(KIND=dp) :: VLoad(3,n), LOAD(n), VL(3), L, Normal(3), Basis(n), DetJ, s
    REAL(KIND=dp) :: f(3), u, v
    REAL(KIND=dp) :: Mass(4,4), rhs(4)
!------------------------------------------------------------------------------
    IF (.NOT.(GetElementFamily(Parent) == 5 .OR. GetElementFamily(Parent) == 8)) &
        CALL Fatal('FaceElementDOFs','A tetrahedral or hexahedral parent element supposed')

    IF (PRESENT(FDOFs)) THEN
      DOFs = FDOFs
    ELSE
      DOFs = 1
    END IF

    ElementCopyCreated = .FALSE.

    SELECT CASE(GetElementFamily(Element))
    CASE(3)

      IF (DOFs > 3) CALL Fatal('FaceElementDOFs', &
          'Cannot yet handle more than 3 DOFs per 3-node face')
 
      IF (PRESENT(SecondFamily)) THEN
        SecondKindBasis = SecondFamily
        IF (SecondKindBasis .AND. (DOFs /= 3) ) &
            CALL Fatal('FaceElementDOFs','3 DOFs per face expected')
      ELSE
        SecondKindBasis = .FALSE.
      END IF
      IF (.NOT. SecondKindBasis .AND. DOFs > 1) &
          CALL Fatal('FaceElementDOFs','An unexpected DOFs count per face')

      IF (SecondKindBasis) THEN
        TetraFaceMap(1,:) = [ 2, 1, 3 ]
        TetraFaceMap(2,:) = [ 1, 2, 4 ]
        TetraFaceMap(3,:) = [ 2, 3, 4 ] 
        TetraFaceMap(4,:) = [ 3, 1, 4 ]

        ActiveFaceMap(1:3) = TetraFaceMap(FaceId,1:3)

        IF (ANY(Element % NodeIndexes(1:3) /= Parent % NodeIndexes(ActiveFaceMap(1:3)))) THEN
          !
          ! The parent element face is indexed differently than the boundary element.
          ! Create a copy of the boundary element which is indexed as the parent element
          ! face so that we can return the values of DOFs in the default order.
          ! Reordering is supposed to be done outside this subroutine. 
          !
          ElementCopyCreated = .TRUE.
          ElementCopy => AllocateElement()
          ElementCopy % Type => Element % Type
          ALLOCATE(ElementCopy % NodeIndexes(3))
          ElementCopy % NodeIndexes(1:3) = Parent % NodeIndexes(ActiveFaceMap(1:3))
          ElementCopy % BodyId = Element % BodyId
          ElementCopy % BoundaryInfo => Element % BoundaryInfo
        ELSE
          ElementCopy => Element
        END IF
      ELSE
        ElementCopy => Element
      END IF
      !CALL GetElementNodes(Nodes, ElementCopy)
      CALL CopyElementNodesFromMesh(Nodes, CurrentModel % Solver % Mesh, &
        ElementCopy % Type % NumberOfNodes, ElementCopy % NodeIndexes)
      
      Load(1:n) = GetReal(BC, Name, stat, ElementCopy)

      i = LEN_TRIM(Name)
      VLoad(1,1:n) = GetReal(BC, Name(1:i)//' 1', stat, ElementCopy)
      VLoad(2,1:n) = GetReal(BC, Name(1:i)//' 2', stat, ElementCopy)
      VLoad(3,1:n) = GetReal(BC, Name(1:i)//' 3', stat, ElementCopy)

      IP = GaussPoints(ElementCopy, 3) ! Feasible for a triangular face
      Integral(:) = 0.0d0
      DO p=1,IP % n
        stat = ElementInfo(ElementCopy, Nodes, IP % u(p), &
            IP % v(p), IP % w(p), DetJ, Basis)
        !
        ! We need a normal that points outwards from the parent element.
        ! The following function call should be consistent with this goal
        ! in the case of a volume-vacuum interface provided a target body
        ! for the normal has not been given to blur the situation.
        ! TO DO: Modify to allow other scenarios 
        !
        Normal = NormalVector(ElementCopy, Nodes, IP % u(p), IP % v(p), .TRUE.)

        VL = MATMUL(VLoad(:,1:n), Basis(1:n))
        L  = SUM(Load(1:n)*Basis(1:n)) + SUM(VL*Normal)

        s = IP % s(p) * DetJ

        IF (SecondKindBasis) THEN
          ! Standard coordinates mapped to the p-element coordinates:
          u = -1.0d0 + 2.0d0*IP % u(p) + IP % v(p)
          v = sqrt(3.0d0)*IP % v(p)
          !
          ! The weight functions for the evaluation of DOFs:
          f(1) = sqrt(3.0d0) * 0.5d0 * (1.0d0 - 2.0d0*u + 1.0d0/3.0d0 - 2.0d0/sqrt(3.0d0)*v)
          f(2) = sqrt(3.0d0) * 0.5d0 * (1.0d0 + 2.0d0*u + 1.0d0/3.0d0 - 2.0d0/sqrt(3.0d0)*v)
          f(3) = sqrt(3.0d0) * (-1.0d0/3.0d0 + 2.0d0/sqrt(3.0d0)*v)

          DO i=1,DOFs
            Integral(i) = Integral(i) + L * f(i) * s
          END DO
        ELSE
          Integral(1) = Integral(1) + L * s
        END IF
      END DO

    CASE(4)
      IF (DOFs /= 4) CALL Fatal('FaceElementDOFs','4 DOFs per 4-node face expected')

      BrickFaceMap(1,:) = (/ 2, 1, 4, 3 /)
      BrickFaceMap(2,:) = (/ 5, 6, 7, 8 /)
      BrickFaceMap(3,:) = (/ 1, 2, 6, 5 /)
      BrickFaceMap(4,:) = (/ 2, 3, 7, 6 /)
      BrickFaceMap(5,:) = (/ 3, 4, 8, 7 /)
      BrickFaceMap(6,:) = (/ 4, 1, 5, 8 /)

      ActiveFaceMap(1:4) = BrickFaceMap(FaceId,1:4)

      IF (ANY(Element % NodeIndexes(1:4) /= Parent % NodeIndexes(ActiveFaceMap(1:4)))) THEN
        !
        ! The parent element face is indexed differently than the boundary element.
        ! Create a copy of the boundary element which is indexed as the parent element
        ! face so that we can return the values of DOFs in the default order.
        ! Reordering is supposed to be done outside this subroutine. 
        !
        ElementCopyCreated = .TRUE.
        ElementCopy => AllocateElement()
        ElementCopy % Type => Element % Type
        ALLOCATE(ElementCopy % NodeIndexes(4))
        ElementCopy % NodeIndexes(1:4) = Parent % NodeIndexes(ActiveFaceMap(1:4))
        ElementCopy % BodyId = Element % BodyId
        ElementCopy % BoundaryInfo => Element % BoundaryInfo
      ELSE
        ElementCopy => Element
      END IF

      !CALL GetElementNodes(Nodes, ElementCopy)
      CALL CopyElementNodesFromMesh(Nodes, CurrentModel % Solver % Mesh, &
        ElementCopy % Type % NumberOfNodes, ElementCopy % NodeIndexes)

      
      Load(1:n) = GetReal(BC, Name, stat, ElementCopy)

      i = LEN_TRIM(Name)
      VLoad(1,1:n) = GetReal(BC, Name(1:i)//' 1', stat, ElementCopy)
      VLoad(2,1:n) = GetReal(BC, Name(1:i)//' 2', stat, ElementCopy)
      VLoad(3,1:n) = GetReal(BC, Name(1:i)//' 3', stat, ElementCopy)

      IP = GaussPoints(ElementCopy, 4)

      Mass = 0.0d0
      rhs = 0.0d0

      DO p=1,IP % n
        stat = ElementInfo(ElementCopy, Nodes, IP % u(p), &
            IP % v(p), IP % w(p), DetJ, Basis)
        !
        ! We need a normal that points outwards from the parent element.
        ! The following function call should be consistent with this goal
        ! in the case of a volume-vacuum interface provided a target body
        ! for the normal has not been given to blur the situation.
        ! TO DO: Modify to allow other scenarios 
        !
        Normal = NormalVector(ElementCopy, Nodes, IP % u(p), IP % v(p), .TRUE.)

        VL = MATMUL(VLoad(:,1:n), Basis(1:n))
        L  = SUM(Load(1:n)*Basis(1:n)) + SUM(VL*Normal)
        s = IP % s(p) * DetJ

        DO i=1,DOFs
          DO j=1,DOFs
            ! Note: here a non-existent DetJ is not a mistake
            Mass(i,j) = Mass(i,j) + Basis(i) * Basis(j) * IP % s(p)
          END DO
          rhs(i) = rhs(i) + L * Basis(i) * s
        END DO
      END DO

      CALL LUSolve(DOFs, Mass(1:DOFs,1:DOFs), rhs(1:DOFs))
      Integral(1:DOFs) = rhs(1:DOFs)

    END SELECT
    IF (ElementCopyCreated) DEALLOCATE(ElementCopy % NodeIndexes)
!------------------------------------------------------------------------------
  END SUBROUTINE FaceElementDOFs
!------------------------------------------------------------------------------


!> In the case of p-approximation, compute the element stiffness matrix and
!> force vector in order to assemble a system of equations for approximating
!> a given Dirichlet condition
!------------------------------------------------------------------------------
  SUBROUTINE LocalBcBDOFs(BC, Element, nd, Name, STIFF, Force )
!------------------------------------------------------------------------------

    IMPLICIT NONE

    TYPE(ValueList_t), POINTER :: BC     !< The list of boundary condition values
    TYPE(Element_t), POINTER :: Element  !< The boundary element handled
    INTEGER :: nd                        !< The number of DOFs in the boundary element
    CHARACTER(LEN=*) :: Name             !< The name of boundary condition
    REAL(KIND=dp) :: STIFF(:,:)          !< The element stiffness matrix
    REAL(KIND=dp) :: Force(:)            !< The element force vector
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: p,q,t
    REAL(KIND=dp) :: Basis(nd)
    REAL(KIND=dp) :: xip,yip,zip,s,DetJ,Load
    LOGICAL :: stat
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------

    ! Get nodes of boundary elements parent and gauss points for boundary
    CALL GetElementNodes( Nodes, Element )
    IP = GaussPoints( Element )

    FORCE(1:nd) = 0.0d0
    STIFF(1:nd,1:nd) = 0.0d0

    DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % u(t), &
          IP % v(t), IP % w(t), DetJ, Basis )

       s = IP % s(t) * DetJ

       ! Get value of boundary condition
       xip = SUM( Basis(1:nd) * Nodes % x(1:nd) )
       yip = SUM( Basis(1:nd) * Nodes % y(1:nd) )
       zip = SUM( Basis(1:nd) * Nodes % z(1:nd) )
       Load = ListGetConstReal( BC, Name, x=xip,y=yip,z=zip )

       ! Build local stiffness matrix and force vector
       DO p=1,nd
          DO q=1,nd
             STIFF(p,q) = STIFF(p,q) + s * Basis(p)*Basis(q)
          END DO
          FORCE(p) = FORCE(p) + s * Load * Basis(p)
       END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalBcBDOFs
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Finishes the bulk assembly of the matrix equation.
!> Optionally save the matrix for later use.
!------------------------------------------------------------------------------
  SUBROUTINE DefaultFinishBulkAssembly( Solver, BulkUpdate, RHSUpdate )
!------------------------------------------------------------------------------
    TYPE(Solver_t), OPTIONAL, TARGET :: Solver
    LOGICAL, OPTIONAL :: BulkUpdate  ! Direct control on whether matrices are saved
    LOGICAL, OPTIONAL :: RHSUpdate   ! Direct control on whether RHS is saved

    TYPE(Solver_t), POINTER :: PSolver
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: Bupd, UpdateRHS, Found
    INTEGER :: n
    CHARACTER(:), ALLOCATABLE :: str
    LOGICAL :: Transient
    REAL(KIND=dp) :: SScond
    INTEGER :: Order
    TYPE(Matrix_t), POINTER :: A

    IF( PRESENT( Solver ) ) THEN
      PSolver => Solver
    ELSE
      PSolver => CurrentModel % Solver
    END IF

    Params => GetSolverParams( PSolver ) 
    
    IF( ListGetLogical( Params,'Bulk Assembly Timing',Found ) ) THEN 
      CALL CheckTimer('BulkAssembly'//GetVarName(PSolver % Variable), Level=5, Delete=.TRUE. ) 
    END IF
        
    ! Reset colouring 
    PSolver % CurrentColour = 0

    IF ( PRESENT(RHSUpdate) ) THEN
      UpdateRHS = RHSUpdate
    ELSE  
      UpdateRHS = .TRUE.
    END IF

    BUpd = .FALSE.
    IF ( PRESENT(BulkUpdate) ) THEN
      BUpd = BulkUpdate 
    ELSE
      BUpd = GetLogical( Params,'Calculate Loads', Found )
      IF( BUpd ) THEN
        str = GetString( Params,'Calculate Loads Slot', Found )
        IF(Found) THEN
          BUpd = ( str == 'bulk assembly')
        END IF
      END IF
      BUpd = BUpd .OR. GetLogical( Params,'Constant Bulk System', Found )
      BUpd = BUpd .OR. GetLogical( Params,'Save Bulk System', Found )
      BUpd = BUpd .OR. GetLogical( Params,'Constant Bulk Matrix', Found )
      BUpd = BUpd .OR. GetLogical( Params,'Constraint Modes Analysis',Found) 
      BUpd = BUpd .OR. GetLogical( Params,'Control Use Loads',Found )
    END IF

    IF( BUpd ) THEN
      str = GetString( Params,'Equation',Found)
      CALL Info('DefaultFinishBulkAssembly','Saving bulk values for: '//str, Level=8 )
      IF( GetLogical( Params,'Constraint Modes Mass Lumping',Found) ) THEN
        CALL CopyBulkMatrix( PSolver % Matrix, BulkMass = .TRUE., BulkRHS = UpdateRHS ) 
      ELSE
        CALL CopyBulkMatrix( PSolver % Matrix, BulkMass = ASSOCIATED(PSolver % Matrix % MassValues), &
            BulkDamp = ASSOCIATED(PSolver % Matrix % DampValues), BulkRHS = UpdateRHS ) 
      END IF
    END IF

    IF( GetLogical( Params,'Bulk System Multiply',Found ) ) THEN	
      CALL Info('DefaultFinishBulkAssembly','Multiplying matrix equation',Level=10)
      CALL LinearSystemMultiply( PSolver )
    END IF

    IF ( ListGetLogical( Params,'Linear System Save',Found )) THEN
      str = GetString( Params,'Linear System Save Slot', Found )
      IF(Found .AND. str == 'bulk assembly') THEN
        CALL SaveLinearSystem( PSolver ) 
      END IF
    END IF

    IF( ListGetLogical( Params,'Linear System Remove Zeros',Found ) ) THEN
      CALL CRS_RemoveZeros( PSolver % Matrix )
    END IF	
    
    IF( ListGetLogical( PSolver % Values,'Boundary Assembly Timing',Found ) ) THEN 
      CALL ResetTimer('BoundaryAssembly'//GetVarName(PSolver % Variable) ) 
    END IF

    IF( InfoActive( 30 ) ) THEN
      A => PSolver % Matrix
      IF(ASSOCIATED(A)) THEN
        CALL VectorValuesRange(A % Values,SIZE(A % Values),'A_bulk')       
        IF(ASSOCIATED(A % rhs)) THEN
          CALL VectorValuesRange(A % rhs,SIZE(A % rhs),'b_bulk')
        END IF
      END IF
    END IF
    
  END SUBROUTINE DefaultFinishBulkAssembly


!------------------------------------------------------------------------------
!> Finished the boundary assembly of the matrix equation.
!> Optionally save the matrix for later use.
!------------------------------------------------------------------------------
  SUBROUTINE DefaultFinishBoundaryAssembly( Solver, BulkUpdate )
!------------------------------------------------------------------------------
    TYPE(Solver_t), OPTIONAL, TARGET :: Solver
    LOGICAL, OPTIONAL :: BulkUpdate
    TYPE(Solver_t), POINTER :: PSolver
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: Bupd, Found, DoIt
    INTEGER :: n
    TYPE(Matrix_t), POINTER :: A
    CHARACTER(:), ALLOCATABLE :: str, name
    TYPE(Variable_t), POINTER ::  x
    INTEGER :: dof
    
    IF( PRESENT( Solver ) ) THEN
      PSolver => Solver
    ELSE
      PSolver => CurrentModel % Solver
    END IF

    Params => GetSolverParams(PSolver)
    A => PSolver % Matrix
    x => PSolver % Variable
   
    ! Set the nodal loads. This needs to be done before any contacts or limiters since otherwise
    ! the given nodal loads will not be considered properly.         
    DoIt = .TRUE.
    IF( ListGetLogical( Params,'Apply Limiter',Found ) ) THEN
      IF(ListGetLogical( Params,'Apply Limiter Loads After',Found) ) DoIt = .FALSE.
    END IF
    IF( DoIt ) THEN
      DO DOF=1,x % DOFs
        name = TRIM(x % name)
        IF (x % DOFs>1) name=ComponentName(name,DOF)              
        CALL SetNodalLoads( CurrentModel,A,A % rhs, &
            Name,DOF,x % DOFs,x % Perm ) 
      END DO
    END IF
    
    IF( ListGetLogical( Params,'Boundary Assembly Timing',Found ) ) THEN 
      CALL CheckTimer('BoundaryAssembly'//GetVarName(x), Level=5, Delete=.TRUE. ) 
    END IF

    ! Reset colouring 
    PSolver % CurrentBoundaryColour = 0

    BUpd = .FALSE.
    IF ( PRESENT(BulkUpdate) ) THEN
      BUpd = BulkUpdate 
      IF ( .NOT. BUpd ) RETURN
    ELSE
      BUpd = GetLogical( Params,'Calculate Loads', Found )
      IF( BUpd ) THEN
        str = GetString( Params,'Calculate Loads Slot', Found )
        IF(Found) THEN
          BUpd = str == 'boundary assembly'
        ELSE
          BUpd = .FALSE.
        END IF
        BUpd = BUpd .OR. GetLogical( Params,'Constant System', Found )
      END IF
    END IF

    IF( BUpd ) THEN
      CALL Info('DefaultFinishBoundaryAssembly','Saving system values for Solver: '&
          //TRIM(x % Name), Level=8)
      CALL CopyBulkMatrix( PSolver % Matrix ) 
    END IF

    IF ( ListGetLogical( Params,'Linear System Save',Found )) THEN
      str=GetString( Params,'Linear System Save Slot', Found )
      IF(Found .AND. str == 'boundary assembly') THEN
        CALL SaveLinearSystem( PSolver ) 
      END IF
    END IF
       
    ! Create contact BCs using mortar conditions.
    !---------------------------------------------------------------------
    IF( ListGetLogical( Params,'Apply Contact BCs',Found) ) THEN
      CALL DetermineContact( PSolver )	
    END IF

    IF( InfoActive( 30 ) ) THEN
      IF(ASSOCIATED(A)) THEN
        CALL VectorValuesRange(A % Values,SIZE(A % Values),'A0')       
        IF( ASSOCIATED( A % rhs) ) THEN
          CALL VectorValuesRange(A % rhs,SIZE(A % rhs),'b0')
        END IF
      END IF
    END IF

  END SUBROUTINE DefaultFinishBoundaryAssembly



!------------------------------------------------------------------------------
!> Finished the assembly of the matrix equation, mainly effects in transient simulation
!> Also may be used to set implicit relaxation on the linear system before
!> applying Dirichlet conditions. If flux corrected transport is applied then
!> make the initial linear system to be of low order. 
!------------------------------------------------------------------------------
  SUBROUTINE DefaultFinishAssembly( Solver )
!------------------------------------------------------------------------------
    TYPE(Solver_t), OPTIONAL, TARGET :: Solver

    INTEGER :: order, n
    LOGICAL :: Found, Transient
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Solver_t), POINTER :: PSolver
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) :: sscond
    CHARACTER(:), ALLOCATABLE :: str
    
    IF( PRESENT( Solver ) ) THEN
      PSolver => Solver
    ELSE
      PSolver => CurrentModel % Solver
    END IF
    A => PSolver % Matrix
    
    Params => GetSolverParams(PSolver)

    ! Nonlinear timestepping needs a copy of the linear system from previous
    ! timestep. Hence the saving of the linear system is enforced. 
    IF( ListGetLogical( Params,'Nonlinear Timestepping', Found ) ) THEN
      CALL Info('DefaultFinishAssembly','Saving system values for Solver: '&
          //TRIM(PSolver % Variable % Name), Level=8)
      CALL CopyBulkMatrix( A ) 
    END IF

    ! Makes a low order matrix of the initial one saving original values
    ! to BulkValues. Also created a lumped mass matrix.
    IF( ListGetLogical( Params,'Linear System FCT',Found ) ) THEN
      IF( PSolver % Variable % Dofs == 1 ) THEN
        CALL CRS_FCTLowOrder( A )
      ELSE
        CALL Fatal('DefaultFinishAssembly','FCT scheme implemented only for one dof')
      END IF
    END IF
    
    IF(GetLogical(Params,'Use Global Mass Matrix',Found)) THEN

      Transient = GetString( CurrentModel % Simulation, 'Simulation Type') == 'transient'
      IF( Transient ) THEN
        SSCond = ListGetCReal( PSolver % Values,'Steady State Condition',Found )
        IF( Found .AND. SSCond > 0.0_dp ) Transient = .FALSE.
      END IF

      IF( Transient ) THEN
        order = GetInteger(Params,'Time Derivative Order',Found)
        IF(.NOT.Found) Order = PSolver % TimeOrder

        SELECT CASE(order)
          
        CASE(1)  
          CALL Default1stOrderTimeGlobal(PSolver)

        CASE(2)
          CALL Default2ndOrderTimeGlobal(PSolver)
        END SELECT
      END IF
    END IF
 
    CALL FinishAssembly( PSolver, A % RHS )

    IF( GetLogical( Params,'Linear System Multiply',Found ) ) THEN
      CALL Info('DefaultFinishAssembly','Multiplying matrix equation',Level=10)
      CALL LinearSystemMultiply( PSolver )
    END IF

    IF( ListCheckPrefix( Params,'Linear System Diagonal Min') ) THEN
      CALL LinearSystemMinDiagonal( PSolver )      
    END IF

    IF ( ListGetLogical( Params,'Linear System Save',Found )) THEN
      str = GetString( Params,'Linear System Save Slot', Found )
      IF(Found .AND. str == 'assembly') THEN
        CALL SaveLinearSystem( PSolver ) 
      END IF
    END IF
        
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultFinishAssembly
!------------------------------------------------------------------------------



!> Returns integration points for edge or face of p element
!------------------------------------------------------------------------------
   FUNCTION GaussPointsBoundary(Element, boundary, np) RESULT(gaussP)
!------------------------------------------------------------------------------
     USE PElementMaps, ONLY : getElementBoundaryMap 
     USE Integration
     IMPLICIT NONE

     ! Parameters
     TYPE(Element_t) :: Element
     INTEGER, INTENT(IN) :: boundary, np

     TYPE( GaussIntegrationPoints_t ) :: gaussP
     TYPE(Nodes_t) :: bNodes 
     TYPE(Element_t) :: mapElement
     TYPE(Element_t), POINTER :: RefElement
     INTEGER :: i, n, eCode, bMap(4)
     REAL(KIND=dp), TARGET :: x(4), y(4), z(4)
     REAL(KIND=dp), POINTER CONTIG :: xP(:), yP(:), zP(:)
     
     SELECT CASE(Element % TYPE % ElementCode / 100)
     ! Triangle and Quadrilateral
     CASE (3,4)
        n = 2
        eCode = 202
     ! Tetrahedron
     CASE (5)
        n = 3
        eCode = 303
     ! Pyramid
     CASE (6)
        ! Select edge element by boundary
        IF (boundary == 1) THEN
           n = 4
           eCode = 404
        ELSE
           n = 3
           eCode = 303
        END IF
     ! Wedge
     CASE (7)
        ! Select edge element by boundary
        SELECT CASE (boundary)
        CASE (1,2)
           n = 3
           eCode = 303
        CASE (3,4,5)
           n = 4
           eCode = 404
        END SELECT
     ! Brick
     CASE (8)
        n = 4
        eCode = 404
     CASE DEFAULT
        WRITE (*,*) 'DefUtils::GaussPointsBoundary: Unsupported element type'
     END SELECT

     ! Get element boundary map
     bMap(1:4) = getElementBoundaryMap(Element, boundary)
     ! Get ref nodes for element
     xP => x
     yP => y
     zP => z
     CALL GetRefPElementNodes( Element % Type,xP,yP,zP )
     ALLOCATE(bNodes % x(n), bNodes % y(n), bNodes % z(n))
        
     ! Set coordinate points of destination
     DO i=1,n
        IF (bMap(i) == 0) CYCLE  
        bNodes % x(i) = x(bMap(i)) 
        bNodes % y(i) = y(bMap(i))
        bNodes % z(i) = z(bMap(i))   
     END DO

     ! Get element to map from
     mapElement % TYPE => GetElementType(eCode)
     CALL AllocateVector(mapElement % NodeIndexes, mapElement % TYPE % NumberOfNodes)

     ! Get gauss points and map them to given element
     gaussP = GaussPoints( mapElement, np )
     
     CALL MapGaussPoints( mapElement, mapElement % TYPE % NumberOfNodes, gaussP, bNodes )
     
     ! Deallocate memory
     DEALLOCATE(bNodes % x, bNodes % y, bNodes % z, MapElement % NodeIndexes)
!------------------------------------------------------------------------------
   END FUNCTION GaussPointsBoundary
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE MapGaussPoints( Element, n, gaussP, Nodes )
!------------------------------------------------------------------------------
     IMPLICIT NONE

     TYPE(Element_t) :: Element
     TYPE(GaussIntegrationPoints_t) :: gaussP
     TYPE(Nodes_t) :: Nodes
     INTEGER :: n

     INTEGER :: i
     REAL(KIND=dp) :: xh,yh,zh,sh, DetJ
     REAL(KIND=dp) :: Basis(n)
     LOGICAL :: stat
     
     ! Map each gauss point from reference element to given nodes
     DO i=1,gaussP % n
        stat = ElementInfo( Element, Nodes, gaussP % u(i), gaussP % v(i), gaussP % w(i), &
             DetJ, Basis )

        IF (.NOT. stat) THEN
           CALL Fatal( 'DefUtils::MapGaussPoints', 'Element to map degenerate')
        END IF

        ! Get mapped points
        sh = gaussP % s(i) * DetJ
        xh = SUM( Basis(1:n) * Nodes % x(1:n) )
        yh = SUM( Basis(1:n) * Nodes % y(1:n) )
        zh = SUM( Basis(1:n) * Nodes % z(1:n) )
        ! Set mapped points
        gaussP % u(i) = xh
        gaussP % v(i) = yh
        gaussP % w(i) = zh
        gaussP % s(i) = sh
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE MapGaussPoints
!------------------------------------------------------------------------------


!>     Calculate global AND local indexes of boundary dofs for given p-element
!>     lying on a boundary. 
!------------------------------------------------------------------------------
   SUBROUTINE getBoundaryIndexesGL( Mesh, Element, BElement, lIndexes, gIndexes, indSize )
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!
!    Type(Mesh_t) :: Mesh
!      INPUT: Finite element mesh containing edges and faces of elements
!
!    Type(Element_t) :: Element
!      INPUT: Parent of boundary element to get indexes for
!
!    Type(Element_t) :: BElement
!      INPUT: Boundary element to get indexes for
!
!    INTEGER :: lIndexes(:), gIndexes(:)
!      OUTPUT: Calculated indexes of boundary element in local and 
!        global system
! 
!    INTEGER :: indSize
!      OUTPUT: Size of created index vector, i.e. how many indexes were created
!        starting from index 1
!    
!------------------------------------------------------------------------------
     IMPLICIT NONE

     ! Parameters
     TYPE(Mesh_t) :: Mesh
     TYPE(Element_t) :: Element
     TYPE(Element_t), POINTER :: BElement
     INTEGER :: indSize, lIndexes(:), gIndexes(:)
     ! Variables
     TYPE(Element_t), POINTER :: Edge, Face
     INTEGER :: i,j,k,n,edgeDofSum, faceOffSet, edgeOffSet(12), localBoundary, nNodes, bMap(4), &
          faceEdgeMap(4)
     LOGICAL :: stat

     ! Clear indexes
     lIndexes = 0
     gIndexes = 0
     
     ! Get boundary map and number of nodes on boundary
     localBoundary = BElement % PDefs % localNumber
     nNodes = BElement % TYPE % NumberOfNodes
     bMap(1:4) = getElementBoundaryMap(Element, localBoundary)
     n = nNodes + 1

     ! Assign local and global node indexes
     lIndexes(1:nNodes) = bMap(1:nNodes)
     gIndexes(1:nNodes) = Element % NodeIndexes(lIndexes(1:nNodes))

     ! Assign rest of indexes
     SELECT CASE(Element % TYPE % DIMENSION)
     CASE (2)
        edgeDofSum = Element % TYPE % NumberOfNodes

        IF (SIZE(lIndexes) < nNodes + Mesh % MaxEdgeDOFs) THEN
           WRITE (*,*) 'DefUtils::getBoundaryIndexes: Not enough space reserved for edge indexes'
           RETURN
        END IF

        DO i=1,Element % TYPE % NumberOfEdges
           Edge => Mesh % Edges( Element % EdgeIndexes(i) )
           
           ! For boundary edge add local and global indexes
           IF (localBoundary == i) THEN
              DO j=1,Edge % BDOFs
                 lIndexes(n) = edgeDofSum + j
                 gIndexes(n) = Mesh % NumberOfNodes + &
                      (Element % EdgeIndexes(localBoundary)-1) * Mesh % MaxEdgeDOFs + j
                 n = n+1
              END DO
              EXIT
           END IF
           
           edgeDofSum = edgeDofSum + Edge % BDOFs 
        END DO
        
        indSize = n - 1
     CASE (3)
        IF (SIZE(lIndexes) < nNodes + (Mesh % MaxEdgeDOFs * BElement % TYPE % NumberOfEdges) +&
             Mesh % MaxFaceDofs) THEN
           WRITE (*,*) 'DefUtils::getBoundaryIndexes: Not enough space reserved for edge indexes'
           RETURN
        END IF

        ! Get offsets for each edge
        edgeOffSet = 0
        faceOffSet = 0
        edgeDofSum = 0
        DO i=1,Element % TYPE % NumberOfEdges
           Edge => Mesh % Edges( Element % EdgeIndexes(i) )
           edgeOffSet(i) = edgeDofSum
           edgeDofSum = edgeDofSum + Edge % BDOFs
        END DO

        ! Get offset for faces
        faceOffSet = edgeDofSum

        ! Add element edges to local indexes
        faceEdgeMap(1:4) = getFaceEdgeMap(Element, localBoundary)
        Face => Mesh % Faces( Element % FaceIndexes(localBoundary) )
        DO i=1,Face % TYPE % NumberOfEdges
           Edge => Mesh % Edges( Face % EdgeIndexes(i) )
           
           IF (Edge % BDOFs <= 0) CYCLE

           DO j=1,Edge % BDOFs
              lIndexes(n) = Element % TYPE % NumberOfNodes + edgeOffSet(faceEdgeMap(i)) + j
              gIndexes(n) = Mesh % NumberOfNodes +&
                  ( Face % EdgeIndexes(i)-1)*Mesh % MaxEdgeDOFs + j
              n=n+1
           END DO
        END DO

        DO i=1,Element % TYPE % NumberOfFaces
           Face => Mesh % Faces( Element % FaceIndexes(i) )
           
           IF (Face % BDOFs <= 0) CYCLE

           ! For boundary face add local and global indexes
           IF (localBoundary == i) THEN
              DO j=1,Face % BDOFs 
                 lIndexes(n) = Element % TYPE % NumberOfNodes + faceOffSet + j
                 gIndexes(n) = Mesh % NumberOfNodes + &
                      Mesh % NumberOfEdges * Mesh % MaxEdgeDOFs + &
                      (Element % FaceIndexes(localBoundary)-1) * Mesh % MaxFaceDOFs + j
                 n=n+1
              END DO
              EXIT
           END IF

           faceOffSet = faceOffSet + Face % BDOFs
        END DO
        
        indSize = n - 1
     END SELECT
   END SUBROUTINE getBoundaryIndexesGL



!------------------------------------------------------------------------------
  SUBROUTINE GetParentUVW( Element,n,Parent,np,U,V,W,Basis )
!------------------------------------------------------------------------------
    TYPE(Element_t) :: Element, Parent
    INTEGER :: n, np
    REAL(KIND=dp) :: U,V,W,Basis(:)
!------------------------------------------------------------------------------
    INTEGER :: i,j
    REAL(KIND=dp), POINTER :: LU(:), LV(:), LW(:)

    LU => Parent % TYPE % NodeU
    LV => Parent % TYPE % NodeV
    LW => Parent % TYPE % NodeW

    U = 0.0_dp
    V = 0.0_dp
    W = 0.0_dp
    DO i = 1,n
      DO j = 1,np
        IF ( Element % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
          U = U + Basis(i) * LU(j)
          V = V + Basis(i) * LV(j)
          W = W + Basis(i) * LW(j)
          EXIT
        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE GetParentUVW
!------------------------------------------------------------------------------


!> Returns flag telling whether Newton linearization is active
!------------------------------------------------------------------------------
  FUNCTION GetNewtonActive( USolver ) RESULT( NewtonActive )
    LOGICAL :: NewtonActive
    TYPE(Solver_t), OPTIONAL, TARGET :: USolver

    IF ( PRESENT( USolver ) ) THEN
      NewtonActive = USolver % NewtonActive
    ELSE
      NewtonActive = CurrentModel % Solver % NewtonActive
    END IF
  END FUNCTION GetNewtonActive


!------------------------------------------------------------------------------
  FUNCTION GetBoundaryEdgeIndex(Boundary,nedge) RESULT(n)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n,nedge
    TYPE(Element_t) :: Boundary
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    Mesh => GetMesh()    
    n = FindBoundaryEdgeIndex(Mesh,Boundary,nedge)
!------------------------------------------------------------------------------
  END FUNCTION GetBoundaryEdgeIndex
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION GetBoundaryFaceIndex(Boundary) RESULT(n)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n
    TYPE(Element_t) :: Boundary
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    Mesh => GetMesh()    
    n = FindBoundaryFaceIndex(Mesh,Boundary)
!------------------------------------------------------------------------------
  END FUNCTION GetBoundaryFaceIndex
!------------------------------------------------------------------------------

  FUNCTION GetNOFColours(USolver) RESULT( ncolours ) 
    IMPLICIT NONE
    TYPE(Solver_t), TARGET, OPTIONAL :: USolver
    INTEGER :: ncolours

    ncolours = 1
    IF ( PRESENT( USolver ) ) THEN
      IF( ASSOCIATED( USolver % ColourIndexList ) ) THEN
        ncolours = USolver % ColourIndexList % n
        USolver % CurrentColour = 0
      END IF
    ELSE
      IF( ASSOCIATED( CurrentModel % Solver % ColourIndexList ) ) THEN
        ncolours = CurrentModel % Solver % ColourIndexList % n 
        CurrentModel % Solver % CurrentColour = 0
      END IF
    END IF

    CALL Info('GetNOFColours','Number of colours: '//I2S(ncolours),Level=12)
  END FUNCTION GetNOFColours

  FUNCTION GetNOFBoundaryColours(USolver) RESULT( ncolours ) 
    IMPLICIT NONE
    TYPE(Solver_t), TARGET, OPTIONAL :: USolver
    INTEGER :: ncolours

    ncolours = 1
    IF ( PRESENT( USolver ) ) THEN
      IF( ASSOCIATED( USolver % BoundaryColourIndexList ) ) THEN
        ncolours = USolver % BoundaryColourIndexList % n
        USolver % CurrentBoundaryColour = 0
      END IF
    ELSE
      IF( ASSOCIATED( CurrentModel % Solver % BoundaryColourIndexList ) ) THEN
        ncolours = CurrentModel % Solver % BoundaryColourIndexList % n 
        CurrentModel % Solver % CurrentBoundaryColour = 0
      END IF
    END IF

    CALL Info('GetNOFBoundaryColours','Number of colours: '//I2S(ncolours),Level=12)
  END FUNCTION GetNOFBoundaryColours
  
  ! Check given colourings are valid and see if they are free of race conditions. 
  SUBROUTINE CheckColourings(Solver)
    IMPLICIT NONE
    TYPE(Solver_t) :: Solver
    
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Graph_t), POINTER :: Colours
    TYPE(Graph_t), POINTER :: BoundaryColours

    TYPE(Element_t), POINTER :: Element
    
    INTEGER, ALLOCATABLE :: Indexes(:), DOFIndexes(:)
    INTEGER :: col, elem, belem, NDOF, dof
    LOGICAL :: errors

    errors = .FALSE.
    
    Mesh => Solver % Mesh
    Colours => Solver % ColourIndexList
    BoundaryColours => Solver % BoundaryColourIndexList
    
    ! Allocate workspace and initialize it
    ALLOCATE(Indexes(MAX(Mesh % NumberOfNodes,&
          Mesh % NumberOfBulkElements*Mesh % MaxElementDOFs,&
          Mesh%NumberOfBoundaryElements*Mesh % MaxElementDOFs)), &
          DOFIndexes(Mesh % MaxElementDOFs))
    Indexes = 0

    ! Check that every element has a colour
    IF (ASSOCIATED(Colours)) THEN
       DO col=1,Colours % N
          DO elem=Colours % Ptr(col), Colours%Ptr(col+1)-1
             Indexes(Colours % Ind(elem))=Indexes(Colours % Ind(elem))+1
          END DO
       END DO
       DO elem=1,Mesh % NumberOfBulkElements
          IF (Indexes(elem) < 1 .OR. Indexes(elem) > 1) THEN
             CALL Warn('CheckColourings','Element not colored correctly: '//i2s(elem))
             errors = .TRUE.
          END IF
       END DO
       
       Indexes = 0
       ! Check that colouring is free of race conditions
       DO col=1,Colours % N
          DO elem=Colours % Ptr(col), Colours%Ptr(col+1)-1
             Element => Mesh % Elements(Colours % Ind(elem))
             NDOF = GetElementDOFs( DOFIndexes, Element, Solver )
             DO dof=1,NDOF
                Indexes(DOFIndexes(dof))=Indexes(DOFIndexes(dof))+1
             END DO
          END DO
          ! Check colouring
          DO dof=1,Mesh % NumberOfBulkElements*Mesh % MaxElementDOFs
             IF (Indexes(dof)>1) THEN
                CALL Warn('CheckColourings','DOF not colored correctly: '//i2s(dof))
                errors = .TRUE.
             END IF
             Indexes(dof)=0
          END DO
       END DO
    END IF

    ! Check that every boundary element has a colour
    IF (ASSOCIATED(BoundaryColours)) THEN
       
       DO col=1,BoundaryColours % N
          DO elem=BoundaryColours % Ptr(col), BoundaryColours%Ptr(col+1)-1
             Indexes(BoundaryColours % Ind(elem))=Indexes(BoundaryColours % Ind(elem))+1
          END DO
       END DO
       DO elem=Mesh % NumberOfBulkElements+1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
          belem = elem - Mesh % NumberOfBulkElements
          IF (Indexes(belem) < 1 .OR. Indexes(belem) > 1) THEN
             CALL Warn('CheckColourings','Boundary element not colored correctly: '//i2s(belem))
             errors = .TRUE.
          END IF
       END DO
       
       Indexes = 0
       ! Check that colouring is free of race conditions
       DO col=1,BoundaryColours % N
          DO elem=BoundaryColours % Ptr(col), BoundaryColours%Ptr(col+1)-1
             Element => Mesh % Elements(Mesh % NumberOfBulkElements + BoundaryColours % Ind(elem))
             NDOF = GetElementDOFs( DOFIndexes, Element, Solver, NotDG=.TRUE. )
             ! WRITE (*,'(2(A,I0))') 'BELEM=', elem, ', CMAP=', BoundaryColours % Ind(elem)
             ! WRITE (*,*) DOFIndexes(1:NDOF)
             DO dof=1,NDOF
                Indexes(DOFIndexes(dof))=Indexes(DOFIndexes(dof))+1
                ! WRITE (*,'(4(A,I0))') 'EID=', Element % ElementIndex,', dof=', dof, &
                !      ', ind=', DOFIndexes(dof), ', colour=', col
             END DO
          END DO
          ! Check colouring
          DO dof=1,Mesh % NumberOfBulkElements*Mesh % MaxElementDOFs
             IF (Indexes(dof)>1) THEN
                CALL Warn('CheckColourings','Boundary DOF not colored correctly: '//i2s(dof))
                errors = .TRUE.
             END IF
             Indexes(dof)=0
          END DO
       END DO
    END IF

    IF (errors) THEN
      CALL Warn('CheckColourings','Mesh colouring contained errors')
    END IF
    
    DEALLOCATE(Indexes, DOFIndexes)
  END SUBROUTINE CheckColourings


END MODULE DefUtils

!> \}  // end of subgroup
!> \}  // end of group

