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
   USE SolverUtils

   IMPLICIT NONE

   INTERFACE DefaultUpdateEquations
     MODULE PROCEDURE DefaultUpdateEquationsR, DefaultUpdateEquationsC
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
   PRIVATE :: GetIndexStore, GetVecIndexStore, GetValueStore
CONTAINS


   FUNCTION GetVersion() RESULT(ch)
     CHARACTER(LEN=:), ALLOCATABLE :: ch
     ch = VERSION
   END FUNCTION GetVersion

   FUNCTION GetSifName(Found) RESULT(ch)
     CHARACTER(LEN=:), ALLOCATABLE :: ch
     LOGICAL, OPTIONAL :: Found     
     ch = ListGetString( CurrentModel % Simulation,'Solver Input File', Found )
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

  FUNCTION GetVecIndexStore() RESULT(ind)
    IMPLICIT NONE
    INTEGER, POINTER CONTIG :: ind(:)
    INTEGER :: istat
     
    IF ( .NOT. ALLOCATED(VecIndexStore) ) THEN
      ALLOCATE( VecIndexStore(ISTORE_MAX_SIZE), STAT=istat )
      VecIndexStore = 0
      IF ( istat /= 0 ) CALL Fatal( 'GetVecIndexStore', &
              'Memory allocation error.' )
    END IF
    ind => VecIndexStore
  END FUNCTION GetVecIndexStore

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
      CurrentElementThread => Element
      OldElement => CurrentElementThread
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
        CALL Fatal('GetIpIndex','Variable is not of type gauss points!')
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
           //TRIM(I2S(n))//' in colour '//TRIM(I2S(Solver % CurrentColour)),Level=20)
     ELSE
       n = Solver % NumberOfActiveElements
       CALL Info('GetNOFActive','Number of active elements: '&
           //TRIM(I2S(n)),Level=20)
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
           //TRIM(I2S(n))//' in colour '//TRIM(I2S(Solver % CurrentBoundaryColour)),Level=20)
     ELSE
       n = Solver % Mesh % NumberOfBoundaryElements
       CALL Info('GetNOFBoundaryActive','Number of active elements: '&
           //TRIM(I2S(n)),Level=20)
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
    TYPE(Element_t), POINTER, OPTIONAL :: UElement

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
  SUBROUTINE GetScalarLocalSolution( x,name,UElement,USolver,tStep, UVariable )
     REAL(KIND=dp) :: x(:)
     CHARACTER(LEN=*), OPTIONAL :: name
     TYPE(Solver_t)  , OPTIONAL, TARGET :: USolver
     TYPE(Element_t),  OPTIONAL, TARGET :: UElement
     TYPE(Variable_t), OPTIONAL, TARGET :: UVariable
     INTEGER, OPTIONAL :: tStep

     REAL(KIND=dp), POINTER :: Values(:)
     TYPE(Variable_t), POINTER :: Variable
     TYPE(Solver_t)  , POINTER :: Solver
     TYPE(Element_t),  POINTER :: Element

     INTEGER :: i, j, n
     INTEGER, POINTER :: Indexes(:)

     Solver => CurrentModel % Solver
     IF ( PRESENT(USolver) ) Solver => USolver

     x = 0.0d0

     IF(.NOT. PRESENT(UVariable)) THEN
       Variable => Solver % Variable
     ELSE
       Variable => UVariable
     END IF
     
     IF ( PRESENT(name) ) THEN
        Variable => VariableGet( Solver % Mesh % Variables, name )
     END IF
     IF ( .NOT. ASSOCIATED( Variable ) ) RETURN

     Element => GetCurrentElement(UElement)

     Values => Variable % Values
     IF ( PRESENT(tStep) ) THEN
       IF ( tStep<0 ) THEN
         IF ( ASSOCIATED(Variable % PrevValues) .AND. -tStep<=SIZE(Variable % PrevValues,2)) &
           Values => Variable % PrevValues(:,-tStep)
       END IF
     END IF

     ! If variable is defined on gauss points return that instead
     IF( Variable % TYPE == Variable_on_gauss_points ) THEN
       j = Element % ElementIndex
       n = Variable % Perm(j+1) - Variable % Perm(j)
       DO i=1,n
         x(i) = Values(Variable % Perm(j) + i)
       END DO
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
       DO i=1,n
         j = Indexes(i)
         IF ( j>0 .AND. j<=SIZE(Variable % Perm) ) THEN
           j = Variable % Perm(j)
           IF ( j>0 ) x(i) = Values(j)
         END IF
       END DO
     ELSE
        DO i=1,n
          j = Indexes(i)
          IF ( j>0 .AND. j<=SIZE(Variable % Values) ) &
            x(i) = Values(Indexes(i))
        END DO
     END IF
  END SUBROUTINE GetScalarLocalSolution



!> Returns a vector field in the nodes of the element
  SUBROUTINE GetVectorLocalSolution( x,name,UElement,USolver,tStep, UVariable )
     REAL(KIND=dp) :: x(:,:)
     CHARACTER(LEN=*), OPTIONAL :: name
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     TYPE(Variable_t), OPTIONAL, TARGET :: UVariable
     INTEGER, OPTIONAL :: tStep

     TYPE(Variable_t), POINTER :: Variable
     TYPE(Solver_t)  , POINTER :: Solver
     TYPE(Element_t),  POINTER :: Element

     INTEGER :: i, j, k, n
     INTEGER, POINTER :: Indexes(:)
     REAL(KIND=dp), POINTER ::  Values(:)

     Solver => CurrentModel % Solver
     IF ( PRESENT(USolver) ) Solver => USolver

     x = 0.0d0

     IF(.NOT. PRESENT(UVariable)) THEN
       Variable => Solver % Variable
     ELSE
       Variable => UVariable
     END IF

     IF ( PRESENT(name) ) THEN
        Variable => VariableGet( Solver % Mesh % Variables, name )
     END IF
     IF ( .NOT. ASSOCIATED( Variable ) ) RETURN

     Element => GetCurrentElement(UElement)

     IF ( ASSOCIATED( Variable ) ) THEN
        Indexes => GetIndexStore()
        IF ( ASSOCIATED(Variable % Solver ) ) THEN
          n = GetElementDOFs( Indexes, Element, Variable % Solver )
        ELSE
          n = GetElementDOFs( Indexes, Element, Solver )
        END IF
        n = MIN( n, SIZE(x) )

        Values => Variable % Values
        IF ( PRESENT(tStep) ) THEN
          IF ( tStep<0 ) THEN
            IF ( ASSOCIATED(Variable % PrevValues) .AND. -tStep<=SIZE(Variable % PrevValues,2)) &
              Values => Variable % PrevValues(:,-tStep)
          END IF
        END IF

        DO i=1,Variable % DOFs
           IF ( ASSOCIATED( Variable % Perm ) ) THEN
             DO j=1,n
               k = Indexes(j)
               IF ( k>0 .AND. k<=SIZE(Variable % Perm) ) THEN
                 k = Variable % Perm(k)
                 IF (k>0) x(i,j) = Values(Variable % DOFs*(k-1)+i)
               END IF
             END DO
           ELSE
              DO j=1,n
                IF ( Variable % DOFs*(Indexes(j)-1)+i <= &
                                SIZE( Variable % Values ) ) THEN
                  x(i,j) = Values(Variable % DOFs*(Indexes(j)-1)+i)
                END IF
              END DO
           END IF
         END DO
     END IF
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
  END FUNCTION


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

     IF ( ASSOCIATED( Variable ) ) THEN
        Indexes => GetIndexStore()
        IF ( ASSOCIATED(Variable % Solver ) ) THEN
          n = GetElementDOFs( Indexes, Element, Variable % Solver )
        ELSE
          n = GetElementDOFs( Indexes, Element, Solver )
        END IF
        n = MIN( n, SIZE(x) )

        Values => Variable % EigenVectors( :, NoEigen )

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

     IF ( ASSOCIATED( Variable ) ) THEN
        Indexes => GetIndexStore()
        IF ( ASSOCIATED(Variable % Solver ) ) THEN
          n = GetElementDOFs( Indexes, Element, Variable % Solver )
        ELSE
          n = GetElementDOFs( Indexes, Element, Solver )
        END IF
        n = MIN( n, SIZE(x) )

        Values => Variable % EigenVectors( :, NoEigen )

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
     END IF
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
     CHARACTER(LEN=MAX_NAME_LEN) :: str

     INTEGER :: i

     IF ( PRESENT( Found ) ) THEN
        str = ListGetString( List, Name, Found )
     ELSE
        str = ListGetString( List, Name )
     END IF
  END FUNCTION


!> Returns an integer by its name if found in the list structure
  FUNCTION GetInteger( List, Name, Found ) RESULT(i)
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found

     INTEGER :: i

     IF ( PRESENT( Found ) ) THEN
        i = ListGetInteger( List, Name, Found )
     ELSE
        i = ListGetInteger( List, Name )
     END IF
  END FUNCTION


!> Returns a logical flag by its name if found in the list structure, otherwise false
  FUNCTION GetLogical( List, Name, Found ) RESULT(l)
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found

     LOGICAL :: l

     IF ( PRESENT( Found ) ) THEN
        l = ListGetLogical( List, Name, Found )
     ELSE
        l = ListGetLogical( List, Name )
     END IF
  END FUNCTION


!> Returns a constant real by its name if found in the list structure
  RECURSIVE FUNCTION GetConstReal( List, Name, Found,x,y,z ) RESULT(r)
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found
     REAL(KIND=dp), OPTIONAL :: x,y,z

     REAL(KIND=dp) :: r,xx,yy,zz

     xx = 0
     yy = 0
     zz = 0
     IF ( PRESENT( x ) ) xx = x
     IF ( PRESENT( y ) ) yy = y
     IF ( PRESENT( z ) ) zz = z

     IF ( PRESENT( Found ) ) THEN
        r = ListGetConstReal( List, Name, Found,xx,yy,zz )
     ELSE
        r = ListGetConstReal( List, Name,x=xx,y=yy,z=zz )
     END IF
  END FUNCTION


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
          IF ( PRESENT( Found ) ) THEN
             x(1:n) = ListGetReal( List, Name, n, NodeIndexes, Found )
          ELSE
             x(1:n) = ListGetReal( List, Name, n, NodeIndexes )
          END IF
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
         IF( Parent % BodyId > 0 .AND. Parent % BodyId <= CurrentModel % NumberOfBodies ) THEN
           mat_id = ListGetInteger( CurrentModel % Bodies(Parent % BodyId) % Values,'Material',GotMat)
         ELSE
           CALL Warn('GetParentMatProp','Invalid parent BodyId '//TRIM(I2S(Parent % BodyId))//&
               ' for element '//TRIM(I2S(Parent % ElementIndex)))
           CYCLE
         END IF
         
         IF(.NOT. GotMat) THEN
           CALL Warn('GetParentMatProp','Parent body '//TRIM(I2S(Parent % BodyId))//' does not have material associated!')
         END IF
         
         IF( mat_id > 0 .AND. mat_id <= CurrentModel % NumberOfMaterials ) THEN
           Material => CurrentModel % Materials(mat_id) % Values
         ELSE
           CALL Warn('GetParentMatProp','Material index '//TRIM(I2S(mat_id))//' not associated to material list')
           CYCLE
         END IF
                             
         IF( .NOT. ASSOCIATED( Material ) ) CYCLE

         IF ( ListCheckPresent( Material,Name) ) THEN
           x(1:n) = ListGetReal(Material, Name, n, Indexes)
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
          IF ( PRESENT( Found ) ) THEN
             x => ListGetConstRealArray( List, Name, Found )
          ELSE
             x => ListGetConstRealArray( List, Name )
          END IF
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
          IF ( PRESENT( Found ) ) THEN
             CALL ListGetRealArray( List, Name, x, n, Element % NodeIndexes, Found )
          ELSE
             CALL ListGetRealArray( List, Name, x, n, Element % NodeINdexes  )
          END IF
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
     ! Thus this routine is more or less obsolite. 
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
     INTEGER :: n
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement

     TYPE(Element_t), POINTER :: Element
     TYPE(Solver_t),  POINTER :: Solver

     INTEGER :: i,j, id, ElemFamily
     LOGICAL :: Found, GB, NeedEdges

     Element => GetCurrentElement( UElement )
     ElemFamily = GetElementFamily(Element)

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
     ELSE
        Solver => CurrentModel % Solver
     END IF

     n = 0
     !IF ( ListGetLogical( Solver % Values, 'Discontinuous Galerkin', Found )) THEN
     IF( Solver % DG ) THEN
       n = Element % DGDOFs
        IF ( n>0 ) RETURN
     END IF


     id =Element % BodyId
     IF ( Id==0 .AND. ASSOCIATED(Element % BoundaryInfo) ) THEN
       IF ( ASSOCIATED(Element % BoundaryInfo % Left) ) &
         id = Element % BoundaryInfo % Left % BodyId

       IF ( ASSOCIATED(Element % BoundaryInfo % Right) ) &
         id = Element % BoundaryInfo % Right % BodyId
     END IF
     IF ( Id==0 ) id=1

     IF ( Solver % Def_Dofs(ElemFamily,id,1)>0 ) n = Element % NDOFs

     NeedEdges = .FALSE.
     DO i=2,SIZE(Solver % Def_Dofs,3)
       IF (Solver % Def_Dofs(ElemFamily, id, i)>=0) THEN
         NeedEdges = .TRUE.
         EXIT
       END IF
     END DO
     IF ( .NOT. NeedEdges ) RETURN

     IF ( ASSOCIATED( Element % EdgeIndexes ) ) THEN
        IF ( Solver % Mesh % MaxEdgeDOFs == Solver % Mesh % MinEdgeDOFs ) THEN
           n =  n + Element % Type % NumberOfEdges * Solver % Mesh % MaxEdgeDOFs
        ELSE
!DIR$ IVDEP
          DO j=1,Element % Type % NumberOFEdges
             n =  n + Solver % Mesh % Edges(Element % EdgeIndexes(j)) % BDOFs
          END DO
       END IF
     END IF

     IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
        IF ( Solver % Mesh % MaxFaceDOFs == Solver % Mesh % MinFaceDOFs ) THEN
           n =  n + Element % Type % NumberOfFaces * Solver % Mesh % MaxFaceDOFs
        ELSE
!DIR$ IVDEP
          DO j=1,Element % Type % NumberOFFaces
             n = n + Solver % Mesh % Faces( Element % FaceIndexes(j) ) % BDOFs
          END DO
        END IF
     END IF

     !GB = ListGetLogical( Solver % Values, 'Bubbles in Global System', Found )
     !IF (.NOT.Found) GB = .TRUE.
     GB = Solver % GlobalBubbles
     
     IF ( GB .OR. ASSOCIATED(Element % BoundaryInfo) ) n=n+MAX(0,Element % BDOFs)
  END FUNCTION GetElementNOFDOFs


!> Returns the number of element degrees of freedom
  FUNCTION GetElementDOFs( Indexes, UElement, USolver,NotDG )  RESULT(NB)
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     INTEGER :: Indexes(:)
     LOGICAL, OPTIONAL  ::  NotDG

     TYPE(Solver_t),  POINTER :: Solver
     TYPE(Element_t), POINTER :: Element, Parent, Edge, Face

     LOGICAL :: Found, GB, DGdisable, NeedEdges
     INTEGER :: nb,i,j,k,id,EDOFs, FDOFs, BDOFs,FaceDOFs, EdgeDOFs, BubbleDOFs, Ind, ElemFamily

     Element => GetCurrentElement(UElement)
     ElemFamily = GetElementFamily(Element)
     
     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
     ELSE
        Solver => CurrentModel % Solver
     END IF

     NB = 0

     DGDisable=.FALSE.
     IF (PRESENT(NotDG)) DGDisable=NotDG

     ! IF ( .NOT.DGDisable .AND. ListGetLogical( Solver % Values, 'Discontinuous Galerkin', Found ) ) THEN
     IF ( .NOT.DGDisable .AND. Solver % DG ) THEN
        DO i=1,Element % DGDOFs
           NB = NB + 1
           Indexes(NB) = Element % DGIndexes(i)
        END DO

        IF ( ASSOCIATED( Element % BoundaryInfo ) ) THEN
           IF ( ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN
              DO i=1,Element % BoundaryInfo % Left % DGDOFs
                 NB = NB + 1
                 Indexes(NB) = Element % BoundaryInfo % Left % DGIndexes(i)
              END DO
           END IF
           IF ( ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
              DO i=1,Element % BoundaryInfo % Right % DGDOFs
                 NB = NB + 1
                 Indexes(NB) = Element % BoundaryInfo % Right % DGIndexes(i)
              END DO
           END IF
        END IF

        IF ( NB > 0 ) RETURN
     END IF

     id =Element % BodyId
     IF ( Id==0 .AND. ASSOCIATED(Element % BoundaryInfo) ) THEN
       IF ( ASSOCIATED(Element % BoundaryInfo % Left) ) &
         id = Element % BoundaryInfo % Left % BodyId

       IF ( ASSOCIATED(Element % BoundaryInfo % Right) ) &
         id = Element % BoundaryInfo % Right % BodyId
     END IF
     IF ( id == 0 ) id=1

     IF ( Solver % Def_Dofs(ElemFamily,id,1)>0 ) THEN
       DO i=1,Element % NDOFs
          NB = NB + 1
          Indexes(NB) = Element % NodeIndexes(i)
       END DO
     END IF

     ! default for nodal elements, if no solver active:
     ! ------------------------------------------------
     IF(.NOT.ASSOCIATED(Solver)) RETURN
     IF(.NOT.ASSOCIATED(Solver % Mesh)) RETURN

     NeedEdges = .FALSE.
     DO i=2,SIZE(Solver % Def_Dofs,3)
       IF (Solver % Def_Dofs(ElemFamily, id, i)>=0) THEN
         NeedEdges = .TRUE.
         EXIT
       END IF
     END DO
     IF ( .NOT. NeedEdges ) RETURN

     FaceDOFs   = Solver % Mesh % MaxFaceDOFs
     EdgeDOFs   = Solver % Mesh % MaxEdgeDOFs
     BubbleDOFs = Solver % Mesh % MaxBDOFs

     IF ( ASSOCIATED(Element % EdgeIndexes) ) THEN
        DO j=1,Element % TYPE % NumberOFEdges
          EDOFs = Solver % Mesh % Edges(Element % EdgeIndexes(j)) % BDOFs
          DO i=1,EDOFs
             NB = NB + 1
             Indexes(NB) = EdgeDOFs*(Element % EdgeIndexes(j)-1) + &
                      i + Solver % Mesh % NumberOfNodes
          END DO
        END DO
     END IF

     IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
        DO j=1,Element % TYPE % NumberOFFaces
           FDOFs = Solver % Mesh % Faces( Element % FaceIndexes(j) ) % BDOFs
           DO i=1,FDOFs
              NB = NB + 1
              Indexes(NB) = FaceDOFs*(Element % FaceIndexes(j)-1) + i + &
                 Solver % Mesh % NumberOfNodes + EdgeDOFs*Solver % Mesh % NumberOfEdges
           END DO
        END DO
     END IF

     GB = Solver % GlobalBubbles 

     IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
       Parent => Element % BoundaryInfo % Left
       IF (.NOT.ASSOCIATED(Parent) ) &
         Parent => Element % BoundaryInfo % Right
       IF (.NOT.ASSOCIATED(Parent) ) RETURN

       SELECT CASE(GetElementFamily(Element))
       CASE(2)
         IF ( ASSOCIATED(Parent % EdgeIndexes) ) THEN
           IF ( isActivePElement(Element) ) THEN
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

           EDOFs = Element % BDOFs
           DO i=1,EDOFs
               NB = NB + 1
             Indexes(NB) = EdgeDOFs*(Parent % EdgeIndexes(Ind)-1) + &
                      i + Solver % Mesh % NumberOfNodes
           END DO
         END IF

       CASE(3,4)
         IF ( ASSOCIATED( Parent % FaceIndexes ) ) THEN
           IF ( isActivePElement(Element) ) THEN
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

           FDOFs = Element % BDOFs
           DO i=1,FDOFs
             NB = NB + 1
             Indexes(NB) = FaceDOFs*(Parent % FaceIndexes(Ind)-1) + i + &
                Solver % Mesh % NumberOfNodes + EdgeDOFs*Solver % Mesh % NumberOfEdges
           END DO
         END IF
       END SELECT
     ELSE IF ( GB ) THEN
        IF ( ASSOCIATED(Element % BubbleIndexes) ) THEN
           DO i=1,Element % BDOFs
              NB = NB + 1
              Indexes(NB) = FaceDOFs*Solver % Mesh % NumberOfFaces + &
                 Solver % Mesh % NumberOfNodes + EdgeDOFs*Solver % Mesh % NumberOfEdges + &
                   Element % BubbleIndexes(i)
           END DO
        END IF
     END IF
  END FUNCTION GetElementDOFs


!> Returns the number of bubble degree of freedom in the active element
  FUNCTION GetElementNOFBDOFs( Element, USolver ) RESULT(n)
    INTEGER :: n
    TYPE(Solver_t), OPTIONAL, POINTER :: USolver
    TYPE(Element_t), OPTIONAL :: Element
    TYPE(Element_t), POINTER  :: CurrElement

    TYPE(Solver_t), POINTER :: Solver

    LOGICAL :: Found, GB

    IF ( PRESENT( USolver ) ) THEN
       Solver => USolver
    ELSE
       Solver => CurrentModel % Solver
    END IF

    !GB = ListGetLogical( Solver % Values, 'Bubbles in Global System', Found )
    !IF (.NOT.Found) GB = .TRUE.
    GB = Solver % GlobalBubbles

    n = 0
    IF ( .NOT. GB ) THEN
      CurrElement => GetCurrentElement(Element)
      n = CurrElement % BDOFs
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

     TYPE(Solver_t),  POINTER  :: Solver
     TYPE(Mesh_t),  POINTER  :: Mesh
     TYPE(Element_t), POINTER :: Element

     Solver => CurrentModel % Solver
     IF ( PRESENT( USolver ) ) Solver => USolver

     Element => GetCurrentElement(UElement)

     Mesh => Solver % Mesh
     IF ( PRESENT( UMesh ) ) Mesh => UMesh
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

        Solver => CurrentModel % Solver
        IF ( PRESENT( USolver ) ) Solver => USolver

        Element => GetCurrentElement(UElement)

        IF ( PRESENT( UMesh ) ) THEN
            Mesh => UMesh
        ELSE
            Mesh => Solver % Mesh
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


!> Is the actice solver solved in the frequency space
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
   CHARACTER(LEN=MAX_NAME_LEN) :: Method
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

   Method = ListGetString( Solver % Values, 'Timestepping Method', Found )
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

   ! The rest of the code in this subroutine is obsolite 
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


   DO i=1,Solver % Matrix % NumberOFRows
     n = 0
     k = 0

     DO j=Solver % Matrix % Rows(i),Solver % Matrix % Rows(i+1)-1
       n=n+1
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

     IF( HasFCT ) THEN
       MASS(1,k) = 0.0_dp
     END IF

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
      
      CALL Bossak2ndOrder( n, Solver % dt, MASS, DAMP, STIFF, &
          FORCE, X(1:n,3), X(1:n,4), X(1:n,5), Solver % Alpha )
      
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

     TYPE(ParEnv_t) :: SParEnv

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
         TRIM(SlaveSolverStr),Level=5)
     
     dt = GetTimeStep()
     Transient = GetString(CurrentModel % Simulation,'Simulation type',Found)=='transient'

     ! store the nonlinear iteration at the outer loop
     iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
     iter = NINT(iterV % Values(1))

     DO j=1,SIZE(SlaveSolverIndexes)
       k = SlaveSolverIndexes(j)
       SlaveSolver => CurrentModel % Solvers(k)

       CALL Info('DefaultSlaveSolvers','Calling slave solver: '//TRIM(I2S(k)),Level=8)
       
       IF(ParEnv % PEs>1) THEN
         SParEnv = ParEnv

         IF(ASSOCIATED(SlaveSolver % Matrix)) THEN
           IF(ASSOCIATED(SlaveSolver % Matrix % ParMatrix) ) THEN
             ParEnv = SlaveSolver % Matrix % ParMatrix % ParEnv
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
         ParEnv = SParEnv
       END IF
     END DO
     iterV % Values = iter       
     CurrentModel % Solver => Solver

   END SUBROUTINE DefaultSlaveSolvers
!------------------------------------------------------------------------------
 
  

!> Performs initialization for matrix equation related to the active solver
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE DefaultInitialize( USolver )
!------------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL, TARGET, INTENT(IN) :: USolver

     TYPE(Solver_t), POINTER :: Solver

     IF ( PRESENT( USolver ) ) THEN
       Solver => USolver
     ELSE
       Solver => CurrentModel % Solver
     END IF

     CALL DefaultSlaveSolvers(Solver,'Slave Solvers') ! this is the initial name of the slot
     CALL DefaultSlaveSolvers(Solver,'Nonlinear Pre Solvers')     
     
     IF(.NOT. ASSOCIATED( Solver % Matrix ) ) THEN
       CALL Fatal('DefaultInitialize','No matrix exists, cannot initialize!')
     END IF

     CALL InitializeToZero( Solver % Matrix, Solver % Matrix % RHS )

     IF( ALLOCATED(Solver % Matrix % ConstrainedDOF) ) THEN
       Solver % Matrix % ConstrainedDOF = .FALSE.
     END IF
       
     IF( ALLOCATED(Solver % Matrix % Dvalues) ) THEN
       Solver % Matrix % Dvalues = 0._dp
     END IF
     
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultInitialize
!------------------------------------------------------------------------------



!> Performs finilizing steps related to the the active solver
!------------------------------------------------------------------------------
  SUBROUTINE DefaultStart( USolver )
!------------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL, TARGET, INTENT(IN) :: USolver
     
     TYPE(Solver_t), POINTER :: Solver

     IF ( PRESENT( USolver ) ) THEN
       Solver => USolver
     ELSE
       Solver => CurrentModel % Solver
     END IF
     
     CALL Info('DefaultStart','Starting solver: '//&
        TRIM(ListGetString(Solver % Values,'Equation')),Level=10)
     
     ! One can run preprocessing solver in this slot.
     !-----------------------------------------------------------------------------
     CALL DefaultSlaveSolvers(Solver,'Pre Solvers')
     
!------------------------------------------------------------------------------
   END SUBROUTINE DefaultStart
!------------------------------------------------------------------------------


  
!> Performs finilizing steps related to the the active solver
!------------------------------------------------------------------------------
  SUBROUTINE DefaultFinish( USolver )
!------------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL, TARGET, INTENT(IN) :: USolver

     TYPE(Solver_t), POINTER :: Solver

     IF ( PRESENT( USolver ) ) THEN
       Solver => USolver
     ELSE
       Solver => CurrentModel % Solver
     END IF

     ! One can run postprocessing solver in this slot.
     !-----------------------------------------------------------------------------
     CALL DefaultSlaveSolvers(Solver,'Post Solvers')

     CALL Info('DefaultFinish','Finished solver: '//&
         TRIM(ListGetString(Solver % Values,'Equation')),Level=8)

!------------------------------------------------------------------------------
   END SUBROUTINE DefaultFinish
!------------------------------------------------------------------------------


!> Solver the matrix equation related to the active solver
!------------------------------------------------------------------------------
  FUNCTION DefaultSolve( USolver, BackRotNT ) RESULT(Norm)
!------------------------------------------------------------------------------
    TYPE(Solver_t), OPTIONAL, TARGET, INTENT(in) :: USolver
    REAL(KIND=dp) :: Norm
    LOGICAL, OPTIONAL, INTENT(in) :: BackRotNT

    TYPE(Matrix_t), POINTER   :: A
    TYPE(Variable_t), POINTER :: x
    REAL(KIND=dp), POINTER CONTIG :: b(:), SOL(:)

    LOGICAL :: Found, BackRot

    TYPE(ValueList_t), POINTER :: Params
    TYPE(Solver_t), POINTER :: Solver
    TYPE(Matrix_t), POINTER :: Ctmp
    CHARACTER(LEN=MAX_NAME_LEN) :: linsolver, precond, dumpfile, saveslot
    INTEGER :: NameSpaceI

    CALL Info('DefaultSolve','Solving linear system with default routines',Level=10)
    
    Solver => CurrentModel % Solver
    Norm = REAL(0, dp)
    IF ( PRESENT( USolver ) ) Solver => USolver

    Params => GetSolverParams(Solver)
    
    NameSpaceI = NINT( ListGetCReal( Params,'Linear System Namespace Number', Found ) )
    IF( NameSpaceI > 0 ) THEN
      CALL Info('DefaultSolve','Linear system namespace number: '//TRIM(I2S(NameSpaceI)),Level=7)
      CALL ListPushNamespace('linsys'//TRIM(I2S(NameSpaceI))//':')
    END IF

    IF( ListCheckPresent( Params, 'Dump system matrix') .OR. &
        ListCheckPresent( Params, 'Dump system RHS') ) THEN
      CALL Error('DefaultSolve','> Dump System Matrix < and > Dump System Rhs < are obsolite')
      CALL Fatal('DefaultSolve','Use > Linear System Save = True < instread!')
    END IF
    
    IF ( ListGetLogical( Params,'Linear System Save',Found )) THEN
      saveslot = GetString( Params,'Linear System Save Slot', Found )
      IF(.NOT. Found .OR. TRIM( saveslot ) == 'solve') THEN
        CALL SaveLinearSystem( Solver ) 
      END IF
    END IF
    
    IF (PRESENT(BackRotNT)) THEN
      BackRot=GetLogical(Params,'Back Rotate N-T Solution',Found)
      IF(.NOT.Found) BackRot=.TRUE.

      IF (BackRot.NEQV.BackRotNT) &
        CALL ListAddLogical(Params,'Back Rotate N-T Solution',BackRotNT)
    END IF

    
    IF( ListGetLogical( Solver % Values,'Harmonic Mode',Found ) ) THEN
      CALL ChangeToHarmonicSystem( Solver )
    END IF

    
    ! Combine the individual projectors into one massive projector
    CALL GenerateConstraintMatrix( CurrentModel, Solver )

    
    IF( GetLogical(Solver % Values,'Linear System Solver Disabled',Found) ) THEN
      CALL Info('DefaultSolve','Solver disabled, exiting early!',Level=10)
      RETURN
    END IF
    

    
    CALL Info('DefaultSolve','Calling SolveSystem for linear solution',Level=20)

    A => Solver % Matrix
    b => A % RHS
    x => Solver % Variable
    SOL => x % Values
    CALL SolveSystem(A,ParMatrix,b,SOL,x % Norm,x % DOFs,Solver)

    ! If flux corrected transport is used then apply the corrector to the system
    IF( GetLogical( Params,'Linear System FCT',Found ) ) THEN
      CALL FCT_Correction( Solver )
    END IF
 
    IF (PRESENT(BackRotNT)) THEN
      IF (BackRot.NEQV.BackRotNT) &
        CALL ListAddLogical(Params,'Back Rotate N-T Solution',BackRot)
    END IF

    Norm = x % Norm

    IF( NameSpaceI > 0 ) CALL ListPopNamespace()


    IF( ListGetLogical( Solver % Values,'Harmonic Mode',Found ) ) THEN
      CALL ChangeToHarmonicSystem( Solver, .TRUE. )
    END IF


    
    ! One can run postprocessing solver in this slot in every nonlinear iteration.
    !-----------------------------------------------------------------------------
    CALL DefaultSlaveSolvers(Solver,'Nonlinear Post Solvers')

    
!------------------------------------------------------------------------------
  END FUNCTION DefaultSolve
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
      CALL Info('DefaultLinesearch','Maximum number of nonlinear iterations reached, giving up after linesearch')
    END IF

    IF( PRESENT( Converged ) ) THEN
      Converged = ( Solver % Variable % NonlinConverged == 1 )  .OR. Last
    END IF

  END FUNCTION DefaultLinesearch



!------------------------------------------------------------------------------
  SUBROUTINE DefaultUpdateEquationsR( G, F, UElement, USolver, BulkUpdate, VecAssembly ) 
!------------------------------------------------------------------------------
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     REAL(KIND=dp) :: G(:,:), f(:)
     LOGICAL, OPTIONAL :: BulkUpdate, VecAssembly

     TYPE(Solver_t), POINTER   :: Solver
     TYPE(Matrix_t), POINTER   :: A
     TYPE(Variable_t), POINTER :: x
     TYPE(Element_t), POINTER  :: Element, P1, P2
     REAL(KIND=dp), POINTER CONTIG   :: b(:)
     REAL(KIND=dp), POINTER :: SaveValues(:)

     CHARACTER(LEN=MAX_NAME_LEN) :: str

     LOGICAL :: Found, BUpd, VecAsm, MCAsm

     INTEGER :: j, n, nd
     INTEGER(KIND=AddrInt) :: Proc
     INTEGER, POINTER CONTIG :: Indexes(:), PermIndexes(:)

#ifndef USE_ISO_C_BINDINGS
     INTERFACE 
       SUBROUTINE ExecLocalProc( Proc, Model, Solver, G, F, Element, n, nd )
         USE Types
         INTEGER(KIND=AddrInt) :: Proc
         TYPE(Model_t)   :: Model
         TYPE(Solver_t)  :: Solver
         TYPE(Element_t) :: Element
         INTEGER :: n, nd
         REAL(KIND=dp) :: G(:,:), F(:)
       END SUBROUTINE ExecLocalProc
     END INTERFACE
#endif

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
            IF ( P1 % PartIndex/=ParEnv % myPE .AND. &
                 P2 % PartIndex/=ParEnv % myPE )RETURN

            IF ( P1 % PartIndex/=ParEnv % myPE .OR. &
                 P2 % PartIndex/=ParEnv % myPE ) THEN
              G=G/2; F=F/2; 
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
       Indexes => GetIndexStore()
       n = GetElementDOFs( Indexes, Element, Solver )
       
       PermIndexes => GetVecIndexStore()
       ! Get permuted indices
!DIR$ IVDEP
       DO j=1,n
         PermIndexes(j) = Solver % Variable % Perm(Indexes(j))
       END DO

       CALL UpdateGlobalEquationsVec( A, G, b, f, n, &
               x % DOFs, PermIndexes, &
               UElement=Element, MCAssembly=MCAsm )
     ELSE
       Indexes => GetIndexStore()
       n = GetElementDOFs( Indexes, Element, Solver )

       IF(Solver % DirectMethod == DIRECT_PERMON) THEN
         CALL UpdateGlobalEquations( A,G,b,f,n,x % DOFs, &
                              x % Perm(Indexes(1:n)), UElement=Element )
         CALL UpdatePermonMatrix( A, G, n, x % DOFs, x % Perm(Indexes(1:n)) )
       ELSE
         CALL UpdateGlobalEquations( A,G,b,f,n,x % DOFs, &
                              x % Perm(Indexes(1:n)), UElement=Element )
       END IF
     END IF
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdateEquationsR
!------------------------------------------------------------------------------


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
  SUBROUTINE DefaultUpdateEquationsC( GC, FC, UElement, USolver, BulkUpdate, MCAssembly ) 
!------------------------------------------------------------------------------
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     COMPLEX(KIND=dp)   :: GC(:,:), FC(:)
     LOGICAL, OPTIONAL :: BulkUpdate, MCAssembly

     TYPE(Solver_t), POINTER   :: Solver
     TYPE(Matrix_t), POINTER   :: A
     TYPE(Variable_t), POINTER :: x
     TYPE(Element_t), POINTER  :: Element, P1, P2
     REAL(KIND=dp), POINTER    :: b(:), SaveValues(:)

     REAL(KIND=dp), POINTER :: G(:,:), F(:)

     LOGICAL :: Found, BUpd

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

     CALL UpdateGlobalEquations( A,G,b,f,n,x % DOFs,x % Perm(Indexes(1:n)) )

     DEALLOCATE( G, F)
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdateEquationsC
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

    LOGICAL :: Found, BUpd

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

    CALL UpdateGlobalForce( Solver % Matrix % RHS, &
       F, n, x % DOFs, x % Perm(Indexes(1:n)), UElement=Element)

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

     CALL UpdateTimeForce( Solver % Matrix,Solver % Matrix % RHS, &
          F, n, x % DOFs, x % Perm(Indexes(1:n)) )
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

     REAL(KIND=dp),  POINTER :: SaveValues(:)

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

     IF ( .NOT. ASSOCIATED( A % PrecValues ) ) THEN
       ALLOCATE( A % PrecValues(SIZE(A % Values)) )
       A % PrecValues = 0.0d0
     END IF

     SaveValues => A % MassValues
     A % MassValues => A % PrecValues
     CALL UpdateMassMatrix( A, M, n, x % DOFs, x % Perm(Indexes(1:n)) )
     A % MassValues => SaveValues
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
     REAL(KIND=dp),  POINTER :: SaveValues(:)

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

     IF ( .NOT. ASSOCIATED( A % PrecValues ) ) THEN
        ALLOCATE( A % PrecValues(SIZE(A % Values)) )
        A % PrecValues = 0.0d0
     END IF

     ALLOCATE( M(DOFs*n,DOFs*n) )
     DO i=1,n*DOFs/2
       DO j=1,n*DOFs/2
         M(2*(i-1)+1, 2*(j-1)+1) =   REAL( MC(i,j) )
         M(2*(i-1)+1, 2*(j-1)+2) = -AIMAG( MC(i,j) )
         M(2*(i-1)+2, 2*(j-1)+1) =  AIMAG( MC(i,j) )
         M(2*(i-1)+2, 2*(j-1)+2) =   REAL( MC(i,j) )
       END DO
     END DO

     SaveValues => A % MassValues
     A % MassValues => A % PrecValues
     CALL UpdateMassMatrix( A, M, n, x % DOFs, x % Perm(Indexes(1:n)) )
     A % MassValues => SaveValues
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

     IF ( .NOT. ASSOCIATED( A % MassValues ) ) THEN
       ALLOCATE( A % MassValues(SIZE(A % Values)) )
       A % MassValues = 0.0d0
     END IF

     CALL UpdateMassMatrix( A, M, n, x % DOFs, x % Perm(Indexes(1:n)) )
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

     IF ( .NOT. ASSOCIATED( A % MassValues ) ) THEN
        ALLOCATE( A % MassValues(SIZE(A % Values)) )
        A % MassValues = 0.0d0
     END IF

     ALLOCATE( M(DOFs*n,DOFs*n) )
     DO i=1,n*DOFs/2
       DO j=1,n*DOFs/2
         M(2*(i-1)+1, 2*(j-1)+1) =   REAL( MC(i,j) )
         M(2*(i-1)+1, 2*(j-1)+2) = -AIMAG( MC(i,j) )
         M(2*(i-1)+2, 2*(j-1)+1) =  AIMAG( MC(i,j) )
         M(2*(i-1)+2, 2*(j-1)+2) =   REAL( MC(i,j) )
       END DO
     END DO

     CALL UpdateMassMatrix( A, M, n, x % DOFs, x % Perm(Indexes(1:n)) )
     DEALLOCATE( M )
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdateMassC
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE DefaultUpdateDampR( B, UElement, USolver ) 
!------------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL,  TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     REAL(KIND=dp)   :: B(:,:)

     TYPE(Solver_t), POINTER   :: Solver
     TYPE(Matrix_t), POINTER   :: A
     TYPE(Variable_t), POINTER :: x
     TYPE(Element_t), POINTER  :: Element, P1, P2

     REAL(KIND=dp), POINTER :: SaveValues(:)

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

     IF ( .NOT. ASSOCIATED( A % DampValues ) ) THEN
        ALLOCATE( A % DampValues(SIZE(A % Values)) ) 
        A % DampValues = 0.0d0
     END IF

     SaveValues => A % MassValues
     A % MassValues => A % DampValues
     CALL UpdateMassMatrix( A, B, n, x % DOFs, x % Perm(Indexes(1:n)) )
     A % MassValues => SaveValues
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

     REAL(KIND=dp), POINTER :: SaveValues(:)

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

     IF ( .NOT. ASSOCIATED( A % DampValues ) ) THEN
        ALLOCATE( A % DampValues(SIZE(A % Values)) ) 
        A % DampValues = 0.0d0
     END IF

     ALLOCATE( B(DOFs*n, DOFs*n) )
     DO i=1,n*DOFs/2
       DO j=1,n*DOFs/2
         B(2*(i-1)+1, 2*(j-1)+1) =   REAL( BC(i,j) )
         B(2*(i-1)+1, 2*(j-1)+2) = -AIMAG( BC(i,j) )
         B(2*(i-1)+2, 2*(j-1)+1) =  AIMAG( BC(i,j) )
         B(2*(i-1)+2, 2*(j-1)+2) =   REAL( BC(i,j) )
       END DO
     END DO

     SaveValues => A % MassValues
     A % MassValues => A % DampValues
     CALL UpdateMassMatrix( A, B, n, x % DOFs, x % Perm(Indexes(1:n)) )
     A % MassValues => SaveValues

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

     REAL(KIND=dp), POINTER :: SaveValues(:)

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

     IF ( .NOT. ASSOCIATED( A % BulkValues ) ) THEN
        ALLOCATE( A % BulkValues(SIZE(A % Values)) ) 
        A % BulkValues = 0.0_dp
     END IF

     IF ( .NOT. ASSOCIATED( A % BulkRHS ) ) THEN
        ALLOCATE( A % BulkRHS(SIZE(A % RHS)) ) 
        A % BulkRHS = 0.0_dp
     END IF

     SaveValues => A % Values
     A % Values => A % BulkValues
     CALL UpdateGlobalEquations( A,B,A % BulkRHS,f,n,x % DOFs,x % Perm(Indexes(1:n)) )
     A % Values => SaveValues
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

     REAL(KIND=dp), POINTER :: SaveValues(:)

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


     IF ( .NOT. ASSOCIATED( A % BulkValues ) ) THEN
        ALLOCATE( A % BulkValues(SIZE(A % Values)) ) 
        A % BulkValues = 0.0_dp
     END IF

     IF ( .NOT. ASSOCIATED( A % BulkRHS ) ) THEN
        ALLOCATE( A % BulkRHS(SIZE(A % RHS)) ) 
        A % BulkRHS = 0.0_dp
     END IF

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

     SaveValues => A % Values
     A % Values => A % BulkValues
     CALL UpdateGlobalEquations( A,B,A % BulkRHS,f,n,x % DOFs,x % Perm(Indexes(1:n)) )
     A % Values => SaveValues

     DEALLOCATE( B )
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
     REAL(KIND=dp) :: xx, s, dval

     INTEGER, ALLOCATABLE :: lInd(:), gInd(:)
     INTEGER :: i, j, k, kk, l, m, n, nd, nb, np, mb, nn, ni, nj, i0
     INTEGER :: EDOFs, DOF, local, numEdgeDofs, istat, n_start, Offset

     LOGICAL :: Flag,Found, ConstantValue, ScaleSystem, DirichletComm
     LOGICAL :: BUpd, PiolaTransform, QuadraticApproximation, SecondKindBasis

     CHARACTER(LEN=MAX_NAME_LEN) :: name

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

     ! Create soft limiters to be later applied by the Dirichlet conditions
     ! This is done only once for each solver, hence the complex logic. 
     !---------------------------------------------------------------------
     IF( ListGetLogical( Solver % Values,'Apply Limiter',Found) ) THEN
       CALL DetermineSoftLimiter( Solver )	
     END IF


     ! Create contact BCs using mortar conditions.
     !---------------------------------------------------------------------
     !IF( ListGetLogical( Solver % Values,'Apply Contact BCs',Found) ) THEN
     !  CALL DetermineContact( Solver )	
     !END IF

     IF(.NOT.ALLOCATED(A % ConstrainedDOF)) THEN
       ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
       A % ConstrainedDOF = .FALSE.
     END IF
       
     IF(.NOT.ALLOCATED(A % Dvalues)) THEN
       ALLOCATE(A % Dvalues(A % NumberOfRows))
       A % Dvalues = 0._dp
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


     IF ( x % DOFs > 1 ) THEN
       CALL SetDirichletBoundaries( CurrentModel,A, b, GetVarName(x),-1,x % DOFs,x % Perm )
     END IF

     CALL Info('DefUtils::DefaultDirichletBCs', &
            'Setting Dirichlet boundary conditions', Level=5)
     
     ! ----------------------------------------------------------------------
     ! Perform some preparations if BCs for p-approximation will be handled: 
     ! ----------------------------------------------------------------------
     ConstantValue = .FALSE.
     DO DOF=1,x % DOFs
        name = x % name
        IF ( x % DOFs > 1 ) name = ComponentName(name,DOF)

        IF( .NOT. ListCheckPresentAnyBC( CurrentModel, name ) ) CYCLE
        
        CALL Info('DefUtils::DefaultDirichletBCs', &
            'p-element preparations: '//TRIM(name), Level=15)
        
        ! Clearing for p-approximation dofs associated with faces & edges:
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


           IF ( isActivePElement(Parent)) THEN
              n = GetElementNOFNodes()
              ! Get indexes of boundary dofs:
              CALL getBoundaryIndexes( Solver % Mesh, Element, Parent, gInd, numEdgeDofs )
              
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
 

     ! -------------------------------------------------------------------------------------
     ! Set BCs for fields which are approximated using H1-conforming basis functions 
     ! (either Lagrange basis or hierarchic p-basis): 
     ! -------------------------------------------------------------------------------------    
     DO DOF=1,x % DOFs
        name = x % name
        IF (x % DOFs>1) name=ComponentName(name,DOF)

        CALL SetNodalLoads( CurrentModel,A, b, &
             Name,DOF,x % DOFs,x % Perm ) ! , Offset ) not yet ?

        CALL SetDirichletBoundaries( CurrentModel, A, b, &
             Name, DOF, x % DOFs, x % Perm, Offset, OffDiagonalMatrix )

        ! ----------------------------------------------------------------------------
        ! Set Dirichlet BCs for edge and face dofs which come from approximating with
        ! p-elements:
        ! ----------------------------------------------------------------------------
        IF( .NOT. ListCheckPresentAnyBC( CurrentModel, name ) ) CYCLE

        CALL Info('DefUtils::DefaultDirichletBCs', &
            'p-element condition setup: '//TRIM(name), Level=15)

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

           
           SELECT CASE(Parent % TYPE % DIMENSION)

           CASE(2)
              ! If no edges do not try to set boundary conditions
              ! @todo This should changed to EXIT
              IF ( .NOT. ASSOCIATED( Solver % Mesh % Edges ) ) CYCLE

              ! If boundary edge has no dofs move on to next edge
              IF (Element % BDOFs <= 0) CYCLE

              ! Number of nodes for this element
              n = Element % TYPE % NumberOfNodes

              ! Get indexes for boundary and values for dofs associated to them
              CALL getBoundaryIndexes( Solver % Mesh, Element, Parent, gInd, numEdgeDofs )
              CALL LocalBcBDOFs( BC, Element, numEdgeDofs, Name, STIFF, Work )

              IF ( Solver % Matrix % Symmetric ) THEN

                DO l=1,n
                    nb = x % Perm( gInd(l) )
                    IF ( nb <= 0 ) CYCLE
                    nb = Offset + x % DOFs * (nb-1) + DOF

                    s = A % Dvalues(nb)
                    DO k=n+1,numEdgeDOFs
                       Work(k) = Work(k) - s*STIFF(k,l)
                    END DO
                 END DO

                 DO k=n+1,numEdgeDOFs
                    DO l=n+1,numEdgeDOFs
                       STIFF(k-n,l-n) = STIFF(k,l)
                    END DO
                    Work(k-n) = Work(k)
                 END DO
                 l = numEdgeDOFs-n
                 IF ( l==1 ) THEN
                   Work(1) = Work(1)/STIFF(1,1)
                 ELSE
                   CALL SolveLinSys(STIFF(1:l,1:l),Work(1:l),l)
                 END IF

                 DO k=n+1,numEdgeDOFs
                   nb = x % Perm( gInd(k) )
                   IF ( nb <= 0 ) CYCLE
                   nb = Offset + x % DOFs * (nb-1) + DOF

                   A % ConstrainedDOF(nb) = .TRUE.
                   A % Dvalues(nb) = Work(k-n)
                 END DO
              ELSE

                ! Contribute this boundary to global system
                 ! (i.e solve global boundary problem)
                 DO k=n+1,numEdgeDofs
                    nb = x % Perm( gInd(k) )
                    IF ( nb <= 0 ) CYCLE
                    nb = Offset + x % DOFs * (nb-1) + DOF
                    A % RHS(nb) = A % RHS(nb) + Work(k) 
                    DO l=1,numEdgeDofs
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
                 END DO
              END IF

            CASE(3)
              ! If no faces present do not try to set boundary conditions
              ! @todo This should be changed to EXIT
              IF ( .NOT. ASSOCIATED(Solver % Mesh % Faces) ) CYCLE

              ! Parameters of element
              n = Element % TYPE % NumberOfNodes

              ! Get global boundary indexes and solve dofs associated to them
              CALL getBoundaryIndexes( Solver % Mesh, Element,  &
                   Parent, gInd, numEdgeDofs )

              ! If boundary face has no dofs skip to next boundary element
              IF (numEdgeDOFs == n) CYCLE

              ! Get local solution
              CALL LocalBcBDofs( BC, Element, numEdgeDofs, Name, STIFF, Work )

              n_start = 1
              IF ( Solver % Matrix % Symmetric ) THEN
                 DO l=1,n
                    nb = x % Perm( gInd(l) )
                    IF ( nb <= 0 ) CYCLE
                    nb = Offset + x % DOFs * (nb-1) + DOF

                    s = A % Dvalues(nb)
                    DO k=n+1,numEdgeDOFs
                       Work(k) = Work(k) - s*STIFF(k,l)
                    END DO
                 END DO
                 n_start=n+1
              END IF

              ! Contribute this entry to global boundary problem
              DO k=n+1,numEdgeDOFs
                 nb = x % Perm( gInd(k) )
                 IF ( nb <= 0 ) CYCLE
                 nb = Offset + x % DOFs * (nb-1) + DOF
                 A % RHS(nb) = A % RHS(nb) + Work(k) 
                 DO l=n_start,numEdgeDOFs
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
              END DO
           END SELECT
         END DO
         
         SaveElement => SetCurrentElement(SaveElement)
     END DO

     ! ----------------------------------------------------------------------------
     ! Set Dirichlet BCs for edge and face dofs which arise from approximating with
     ! edge (curl-conforming) or face (div-conforming) elements:
     ! ----------------------------------------------------------------------------
     QuadraticApproximation = ListGetLogical(Solver % Values, 'Quadratic Approximation', Found)
     SecondKindBasis = ListGetLogical(Solver % Values, 'Second Kind Basis', Found)
     DO DOF=1,x % DOFs
        name = x % name
        IF (x % DOFs>1) name=ComponentName(name,DOF)
        
        IF ( .NOT. ListCheckPrefixAnyBC(CurrentModel, TRIM(Name)//' {e}') .AND. &
            .NOT. ListCheckPrefixAnyBC(CurrentModel, TRIM(Name)//' {f}') ) CYCLE

        CALL Info('SetDefaultDirichlet','Setting edge and face dofs',Level=15)

        SaveElement => GetCurrentElement()
        DO i=1,Solver % Mesh % NumberOfBoundaryElements
           Element => GetBoundaryElement(i)

           BC => GetBC()
           IF ( .NOT.ASSOCIATED(BC) ) CYCLE
           IF ( .NOT. ListCheckPrefix(BC, TRIM(Name)//' {e}') .AND. &
                .NOT. ListCheckPrefix(BC, TRIM(Name)//' {f}') ) CYCLE

           ! Get parent element:
           ! -------------------
           Parent => Element % BoundaryInfo % Left
           IF ( .NOT. ASSOCIATED( Parent ) ) THEN
              Parent => Element % BoundaryInfo % Right
           END IF
           IF ( .NOT. ASSOCIATED( Parent ) )   CYCLE
           np = Parent % TYPE % NumberOfNodes

           IF ( ListCheckPrefix(BC, TRIM(Name)//' {e}') ) THEN
             !--------------------------------------------------------------------------------
              ! We now devote this branch for handling edge (curl-conforming) finite elements 
              ! which, in addition to edge DOFs, may also have DOFs associated with faces. 
              !--------------------------------------------------------------------------------
              IF ( ASSOCIATED( Solver % Mesh % Edges ) ) THEN
                 SELECT CASE(GetElementFamily())
                 CASE(2)
                   DO j=1,Parent % TYPE % NumberOfEdges
                     Edge => Solver % Mesh % Edges(Parent % EdgeIndexes(j))
                     n = 0
                     DO k=1,Element % TYPE % NumberOfNodes
                       DO l=1,Edge % TYPE % NumberOfNodes
                         IF ( Element % NodeIndexes(k)==Edge % NodeIndexes(l)) n=n+1
                       END DO
                     END DO
                     IF ( n==Element % TYPE % NumberOfNodes ) EXIT
                   END DO

                   IF ( .NOT. ActiveBoundaryElement(Edge) ) CYCLE                  

                   EDOFs = Edge % BDOFs     ! The number of DOFs associated with edges
                   n = Edge % TYPE % NumberOfNodes
                   CALL LocalBcIntegral(BC,Edge,n,Parent,np,TRIM(Name)//' {e}',Work, &
                       EDOFs, SecondKindBasis)

                   n=GetElementDOFs(gInd,Edge)

                   n_start = Solver % Def_Dofs(2,Parent % BodyId,1)*Edge % NDOFs
                   DO j=1,EDOFs
                     k = n_start + j
                     nb = x % Perm(gInd(k))
                     IF ( nb <= 0 ) CYCLE
                     nb = Offset + x % DOFs*(nb-1) + DOF

                     A % ConstrainedDOF(nb) = .TRUE.
                     A % Dvalues(nb) = Work(j) 
                   END DO

                 CASE(3,4)
                   !Check that solver % mesh % faces exists?
                   DO j=1,Parent % TYPE % NumberOfFaces
                     Face => Solver % Mesh % Faces(Parent % FaceIndexes(j))
                     IF ( GetElementFamily(Element)==GetElementFamily(Face) ) THEN
                       n = 0
                       DO k=1,Element % TYPE % NumberOfNodes
                         DO l=1,Face % TYPE % NumberOfNodes
                           IF ( Element % NodeIndexes(k)==Face % NodeIndexes(l)) n=n+1
                         END DO
                       END DO
                       IF ( n==Face % TYPE % NumberOfNodes ) EXIT
                     END IF
                   END DO

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
                     n = Edge % TYPE % NumberOfNodes

                     CALL LocalBcIntegral(BC, Edge, n, Parent, np, TRIM(Name)//' {e}', &
                         Work(i0+1:i0+EDOFs), EDOFs, SecondKindBasis)
                     
                     n = GetElementDOFs(gInd,Edge)

                     n_start = Solver % Def_Dofs(2,Parent % BodyId,1)*Edge % NDOFs

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
                   ! We use the variational equation (u x n,v) = (g x n - u0 x n,v) where
                   ! u0 denotes the part of the interpolating function u+u0 which is already 
                   ! known and v is a test function for the Galerkin method.
                   ! ---------------------------------------------------------------------
                   IF (Face % BDOFs > 0) THEN
                     EDOFs = i0 ! The count of edge DOFs set so far
                     n = Face % TYPE % NumberOfNodes

                     CALL SolveLocalFaceDOFs(BC, Face, n, TRIM(Name)//' {e}', Work, EDOFs, &
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
             END IF

           IF ( ListCheckPresent(BC, TRIM(Name)//' {f}') ) THEN
              !--------------------------------------------------------------------------
              ! To do: this branch should be able to handle BCs for face (div-conforming)
              ! elements. 
              !--------------------------------------------------------------------------
           END IF

         END DO
         SaveElement => SetCurrentElement(SaveElement)
      END DO


      ! Add the possible constraint modes structures
      !----------------------------------------------------------
      IF ( GetLogical(Solver % Values,'Constraint Modes Analysis',Found) ) THEN
        CALL SetConstraintModesBoundaries( CurrentModel, A, b, x % Name, x % DOFs, x % Perm )
      END IF
      
     
#ifdef HAVE_FETI4I
     IF(C_ASSOCIATED(A % PermonMatrix)) THEN
       CALL Info('DefUtils::DefaultDirichletBCs','Permon matrix, Dirichlet conditions registered but not set!', Level=5)
       RETURN
     END IF
#endif

     ! This is set outside so that it can be called more flexibilly
     CALL EnforceDirichletConditions( Solver, A, b )
     
 
     CALL Info('DefUtils::DefaultDirichletBCs','Dirichlet boundary conditions set', Level=5)
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultDirichletBCs
!------------------------------------------------------------------------------


! Solves a small dense linear system using Lapack routines
!------------------------------------------------------------------------------
  SUBROUTINE SolveLinSys( A, x, n )
!------------------------------------------------------------------------------
     INTEGER :: n
     REAL(KIND=dp) :: A(n,n), x(n), b(n)

     INTERFACE
       SUBROUTINE SolveLapack( N,A,x )
         INTEGER  N
         DOUBLE PRECISION  A(n*n),x(n)
       END SUBROUTINE
     END INTERFACE

!------------------------------------------------------------------------------
     SELECT CASE(n)
     CASE(1)
       x(1) = x(1) / A(1,1)
     CASE(2)
       b = x
       CALL SolveLinSys2x2(A,x,b)
     CASE(3)
       b = x
       CALL SolveLinSys3x3(A,x,b)
     CASE DEFAULT
       CALL SolveLapack(n,A,x)
     END SELECT
!------------------------------------------------------------------------------
  END SUBROUTINE SolveLinSys
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> This subroutine computes the values of DOFs that are associated with 
!> mesh edges in the case of curl-conforming (edge) finite elements, so that
!> the edge finite element interpolant of the BC data can be constructed. 
!> The values of the DOFs are defined as D = S*(g.t,v)_E where g.t is tangential 
!> component of data, v is a polynomial on the edge E, and S reverts sign
!> if necessary.
!------------------------------------------------------------------------------
  SUBROUTINE LocalBcIntegral(BC, Element, n, Parent, np, Name, Integral, EDOFs, &
      SecondFamily)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(ValueList_t), POINTER :: BC !< The list of boundary condition values
    TYPE(Element_t) :: Element       !< The boundary element handled
    INTEGER :: n                     !< The number of boundary element nodes
    TYPE(Element_t) :: Parent        !< The parent element of the boundary element
    INTEGER :: np                    !< The number of parent element nodes
    CHARACTER(LEN=*) :: Name         !< The name of boundary condition
    REAL(KIND=dp) :: Integral(:)     !< The values of DOFs
    INTEGER, OPTIONAL :: EDOFs       !< The number of DOFs
    LOGICAL, OPTIONAL :: SecondFamily !< To select the edge element family
!------------------------------------------------------------------------------
    TYPE(Nodes_t), SAVE :: Nodes, Pnodes
    TYPE(ElementType_t), POINTER :: SavedType
    TYPE(GaussIntegrationPoints_t) :: IP

    LOGICAL :: Lstat, RevertSign, SecondKindBasis
    INTEGER, POINTER :: Edgemap(:,:)
    INTEGER :: i,j,k,p,DOFs

    REAL(KIND=dp) :: Basis(n),Load(n),Vload(3,1:n),VL(3),t(3)
    REAL(KIND=dp) :: u,v,L,s,DetJ
!------------------------------------------------------------------------------
    DOFs = 1
    IF (PRESENT(EDOFs)) THEN
      IF (EDOFs > 2) THEN
        CALL Fatal('LocalBCIntegral','Cannot handle more than 2 DOFs per edge')
      ELSE
        DOFs = EDOFs
      END IF
    END IF   

    IF (PRESENT(SecondFamily)) THEN
      SecondKindBasis = SecondFamily
      IF (SecondKindBasis .AND. (DOFs /= 2) ) &
          CALL Fatal('LocalBCIntegral','2 DOFs per edge expected')
    ELSE
      SecondKindBasis = .FALSE.
    END IF

    ! Get the nodes of the boundary and parent elements:
    CALL GetElementNodes(Nodes, Element)
    CALL GetElementNodes(PNodes, Parent)

    RevertSign = .FALSE.
    EdgeMap => GetEdgeMap(GetElementFamily(Parent))
    DO i=1,SIZE(EdgeMap,1)
      j=EdgeMap(i,1)
      k=EdgeMap(i,2)
      IF ( Parent % NodeIndexes(j)==Element % NodeIndexes(1) .AND. &
          Parent % NodeIndexes(k)==Element % NodeIndexes(2) ) THEN
        EXIT
      ELSE IF (Parent % NodeIndexes(j)==Element % NodeIndexes(2) .AND. &
          Parent % NodeIndexes(k)==Element % NodeIndexes(1) ) THEN
        RevertSign = .TRUE.
        EXIT
      END IF
    END DO

    Load(1:n) = GetReal( BC, Name, Lstat, Element )

    i = LEN_TRIM(Name)
    VLoad(1,1:n)=GetReal(BC,Name(1:i)//' 1',Lstat,element)
    VLoad(2,1:n)=GetReal(BC,Name(1:i)//' 2',Lstat,element)
    VLoad(3,1:n)=GetReal(BC,Name(1:i)//' 3',Lstat,element)

    t(1) = PNodes % x(k) - PNodes % x(j)
    t(2) = PNodes % y(k) - PNodes % y(j)
    t(3) = PNodes % z(k) - PNodes % z(j)
    t = t/SQRT(SUM(t**2))

    SavedType => Element % TYPE
    IF ( GetElementFamily()==1 ) Element % TYPE=>GetElementType(202)
      
    Integral(1:DOFs) = 0._dp
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
        Integral(1)=Integral(1)+s*(L+SUM(VL*t))*v
        v = 0.5d0*(1.0d0+sqrt(3.0d0)*u)
        Integral(2)=Integral(2)+s*(L+SUM(VL*t))*v
      ELSE
        Integral(1)=Integral(1)+s*(L+SUM(VL*t))

        IF (DOFs>1) THEN
          v = Basis(2)-Basis(1)
          IF (RevertSign) v = -1.0d0*v
          Integral(2)=Integral(2)+s*(L+SUM(VL*t))*v
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
  END SUBROUTINE LocalBcIntegral
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

    REAL(KIND=dp) :: Basis(n),Vload(3,1:n),VL(3),Normal(3)
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

    CALL GetElementNodes(Nodes, Element)

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
    CHARACTER(LEN=MAX_NAME_LEN) :: Name  !< The name of boundary condition
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
!> Finished the bulk assembly of the matrix equation.
!> Optionally save the matrix for later use.
!------------------------------------------------------------------------------
  SUBROUTINE DefaultFinishBulkAssembly( Solver, BulkUpdate )
!------------------------------------------------------------------------------
    TYPE(Solver_t), OPTIONAL, TARGET :: Solver
    LOGICAL, OPTIONAL :: BulkUpdate

    TYPE(Solver_t), POINTER :: PSolver
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: Bupd, Found
    INTEGER :: n
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    LOGICAL :: Transient
    REAL(KIND=dp) :: SScond
    INTEGER :: Order

    IF( PRESENT( Solver ) ) THEN
      PSolver => Solver
    ELSE
      PSolver => CurrentModel % Solver
    END IF

    Params => GetSolverParams( PSolver ) 

    ! Reset colouring 
    PSolver % CurrentColour = 0

    BUpd = .FALSE.
    IF ( PRESENT(BulkUpdate) ) THEN
      BUpd = BulkUpdate 
    ELSE
      BUpd = GetLogical( Params,'Calculate Loads', Found )
      IF( BUpd ) THEN
        str = GetString( Params,'Calculate Loads Slot', Found )
        IF(Found) THEN
          BUpd = ( TRIM( str ) == 'bulk assembly')
        END IF
      END IF
      BUpd = BUpd .OR. GetLogical( Params,'Constant Bulk System', Found )
      BUpd = BUpd .OR. GetLogical( Params,'Save Bulk System', Found )
      BUpd = BUpd .OR. GetLogical( Params,'Constant Bulk Matrix', Found )
      BUpd = BUpd .OR. GetLogical( Params,'Constraint Modes Analysis',Found) 
    END IF

    IF( BUpd ) THEN
      str = GetString( Params,'Equation',Found)
      CALL Info('DefaultFinishBulkAssembly','Saving bulk values for: '//TRIM(str), Level=5 )
      IF( GetLogical( Params,'Constraint Modes Mass Lumping',Found) ) THEN
        CALL CopyBulkMatrix( PSolver % Matrix, BulkMass = .TRUE. ) 
      ELSE
        CALL CopyBulkMatrix( PSolver % Matrix ) 
      END IF
    END IF

    IF( GetLogical( Params,'Bulk System Multiply',Found ) ) THEN	
      CALL Info('DefaultFinishAssembly','Multiplying matrix equation',Level=10)
      CALL LinearSystemMultiply( PSolver )
    END IF

    IF ( ListGetLogical( Params,'Linear System Save',Found )) THEN
      str = GetString( Params,'Linear System Save Slot', Found )
      IF(Found .AND. TRIM( str ) == 'bulk assembly') THEN
        CALL SaveLinearSystem( PSolver ) 
      END IF
    END IF

    IF( ListGetLogical( Params,'Linear System Remove Zeros',Found ) ) THEN
      CALL CRS_RemoveZeros( PSolver % Matrix )
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
    LOGICAL :: Bupd, Found
    INTEGER :: n
    CHARACTER(LEN=MAX_NAME_LEN) :: str

    IF( PRESENT( Solver ) ) THEN
      PSolver => Solver
    ELSE
      PSolver => CurrentModel % Solver
    END IF

    Params => GetSolverParams(PSolver)

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
          BUpd = ( TRIM( str ) == 'boundary assembly') 
        ELSE
          BUpd = .FALSE.
        END IF
        BUpd = BUpd .OR. GetLogical( Params,'Constant System', Found )
      END IF
    END IF

    IF( BUpd ) THEN
      CALL Info('DefaultFinishBoundaryAssembly','Saving system values for Solver: '&
          //TRIM(PSolver % Variable % Name), Level=8)
      CALL CopyBulkMatrix( PSolver % Matrix ) 
    END IF

    IF ( ListGetLogical( Params,'Linear System Save',Found )) THEN
      str = GetString( Params,'Linear System Save Slot', Found )
      IF(Found .AND. TRIM( str ) == 'boundary assembly') THEN
        CALL SaveLinearSystem( PSolver ) 
      END IF
    END IF

    ! Create contact BCs using mortar conditions.
    !---------------------------------------------------------------------
    IF( ListGetLogical( PSolver % Values,'Apply Contact BCs',Found) ) THEN
      CALL DetermineContact( PSolver )	
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
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    REAL(KIND=dp) :: sscond

    IF( PRESENT( Solver ) ) THEN
      PSolver => Solver
    ELSE
      PSolver => CurrentModel % Solver
    END IF

    Params => GetSolverParams(PSolver)

    ! Nonlinear timestepping needs a copy of the linear system from previous
    ! timestep. Hence the saving of the linear system is enforced. 
    IF( ListGetLogical( Params,'Nonlinear Timestepping', Found ) ) THEN
      CALL Info('DefaultFinishAssembly','Saving system values for Solver: '&
          //TRIM(PSolver % Variable % Name), Level=8)
      CALL CopyBulkMatrix( PSolver % Matrix ) 
    END IF

    ! Makes a low order matrix of the initial one saving original values
    ! to BulkValues. Also created a lumped mass matrix.
    IF( ListGetLogical( Params,'Linear System FCT',Found ) ) THEN
      IF( PSolver % Variable % Dofs == 1 ) THEN
        CALL CRS_FCTLowOrder( PSolver % Matrix )
      ELSE
        CALL Fatal('DefaultFinishAssembly','FCT scheme implemented only for one dof')
      END IF
    END IF

    IF(GetLogical(Params,'Use Global Mass Matrix',Found)) THEN

      Transient = ( ListGetString( CurrentModel % Simulation, 'Simulation Type' ) == 'transient')
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
 
    CALL FinishAssembly( PSolver, PSolver % Matrix % RHS )

    IF( GetLogical( Params,'Linear System Multiply',Found ) ) THEN
      CALL Info('DefaultFinishAssembly','Multiplying matrix equation',Level=10)
      CALL LinearSystemMultiply( PSolver )
    END IF

    IF( ListCheckPrefix( Params,'Linear System Diagonal Min') ) THEN
      CALL LinearSystemMinDiagonal( PSolver )      
    END IF


    IF ( ListGetLogical( Params,'Linear System Save',Found )) THEN
      str = GetString( Params,'Linear System Save Slot', Found )
      IF(Found .AND. TRIM( str ) == 'assembly') THEN
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
     CALL GetRefPElementNodes( Element,xP,yP,zP )
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
           WRITE (*,*) 'DefUtils::MapGaussPoints: Element to map degenerate'
           STOP
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

!> Calculate global indexes of boundary dofs for given p-element lying on 
!> a boundary.
!------------------------------------------------------------------------------
   SUBROUTINE getBoundaryIndexes( Mesh, Element, Parent, Indexes, indSize )
!------------------------------------------------------------------------------
!
!    Type(Mesh_t) :: Mesh
!      INPUT: Finite element mesh containing edges and faces of elements
!
!    Type(Element_t) :: Element
!      INPUT: Boundary element to get indexes for
!
!    Type(Element_t) :: Parent
!      INPUT: Parent of boundary element to get indexes for
!
!    INTEGER :: Indexes(:)
!      OUTPUT: Calculated indexes of boundary element in global system
! 
!    INTEGER :: indSize
!      OUTPUT: Size of created index vector, i.e. how many indexes were created
!        starting from index 1
!------------------------------------------------------------------------------
     IMPLICIT NONE

     ! Parameters
     TYPE(Mesh_t) :: Mesh
     TYPE(Element_t) :: Parent
     TYPE(Element_t), POINTER :: Element
     INTEGER :: indSize, Indexes(:)
     
     ! Variables
     TYPE(Element_t), POINTER :: Edge, Face
     INTEGER :: i,j,n

     ! Clear indexes
     Indexes = 0
     n = Element % TYPE % NumberOfNodes

     ! Nodal indexes
     Indexes(1:n) = Element % NodeIndexes(1:n)

     ! Assign rest of indexes if necessary
     SELECT CASE(Parent % TYPE % DIMENSION)
     CASE (1)
       indSize = n 
     CASE (2)
        ! Add index for each bubble dof in edge
        DO i=1,Element % BDOFs
           n = n+1
           
           IF (SIZE(Indexes) < n) THEN
              CALL Warn('DefUtils::getBoundaryIndexes','Not enough space reserved for indexes')
              RETURN
           END IF

           Indexes(n) = Mesh % NumberOfNodes + &
                (Parent % EdgeIndexes(Element % PDefs % localNumber)-1) * Mesh % MaxEdgeDOFs + i
        END DO
     
        indSize = n 
     CASE (3)
        ! Get boundary face
        Face => Mesh % Faces( Parent % FaceIndexes(Element % PDefs % localNumber) )
        
        ! Add indexes of faces edges 
        DO i=1, Face % TYPE % NumberOfEdges
           Edge => Mesh % Edges( Face % EdgeIndexes(i) )
           
           ! If edge has no dofs jump to next edge
           IF (Edge % BDOFs <= 0) CYCLE

           DO j=1,Edge % BDOFs
              n = n + 1
              
              IF (SIZE(Indexes) < n) THEN
                 CALL Warn('DefUtils::getBoundaryIndexes','Not enough space reserved for indexes')
                 RETURN
              END IF
              
              Indexes(n) = Mesh % NumberOfNodes +&
                  ( Face % EdgeIndexes(i)-1)*Mesh % MaxEdgeDOFs + j
           END DO
        END DO
               
        ! Add indexes of faces bubbles
        DO i=1,Face % BDOFs
           n = n + 1

           IF (SIZE(Indexes) < n) THEN
              CALL Warn('DefUtils::getBoundaryIndexes','Not enough space reserved for indexes')
              RETURN
           END IF

           Indexes(n) = Mesh % NumberOfNodes + &
                Mesh % NumberOfEdges * Mesh % MaxEdgeDOFs + &
                (Parent % FaceIndexes( Element % PDefs % localNumber )-1) * Mesh % MaxFaceDOFs + i
        END DO        

        indSize = n
     CASE DEFAULT
        CALL Fatal('DefUtils::getBoundaryIndexes','Unsupported dimension')
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE getBoundaryIndexes
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
 FUNCTION GetEdgeMap( ElementFamily ) RESULT(EdgeMap)
!------------------------------------------------------------------------------
    INTEGER :: ElementFamily
    INTEGER, POINTER :: EdgeMap(:,:)

    INTEGER, TARGET :: LineEM(1,2)
    INTEGER, TARGET :: TriangleEM(3,2)
    INTEGER, TARGET :: QuadEM(4,2)
    INTEGER, TARGET :: TetraEM(6,2)
    INTEGER, TARGET :: PyramidEM(8,2)
    INTEGER, TARGET :: WedgeEM(9,2)
    INTEGER, TARGET :: BrickEM(12,2)

    LOGICAL :: Initialized(8) = .FALSE.
  
    SAVE LineEM, TriangleEM, WedgeEM, BrickEM, TetraEM, QuadEM, PyramidEM, Initialized

    SELECT CASE(ElementFamily)
    CASE(2)
      EdgeMap => LineEM
    CASE(3)
      EdgeMap => TriangleEM
    CASE(4) 
      EdgeMap => QuadEM
    CASE(5) 
      EdgeMap => TetraEM
    CASE(6) 
      EdgeMap => PyramidEM
    CASE(7) 
      EdgeMap => WedgeEM
    CASE(8) 
      EdgeMap => BrickEM
    END SELECT
 
     IF ( .NOT. Initialized(ElementFamily)) THEN
       Initialized(ElementFamily) = .TRUE.
       SELECT CASE(ElementFamily)
       CASE(2)
         EdgeMap(1,:) = (/ 1,2 /)

       CASE(3)
         EdgeMap(1,:) = (/ 1,2 /)
         EdgeMap(2,:) = (/ 2,3 /)
         EdgeMap(3,:) = (/ 3,1 /)

       CASE(4)
         EdgeMap(1,:) = (/ 1,2 /)
         EdgeMap(2,:) = (/ 2,3 /)
         EdgeMap(3,:) = (/ 3,4 /)
         EdgeMap(4,:) = (/ 4,1 /)

       CASE(5)
         EdgeMap(1,:) = (/ 1,2 /)
         EdgeMap(2,:) = (/ 2,3 /)
         EdgeMap(3,:) = (/ 3,1 /)
         EdgeMap(4,:) = (/ 1,4 /)
         EdgeMap(5,:) = (/ 2,4 /)
         EdgeMap(6,:) = (/ 3,4 /)

       CASE(6)
         EdgeMap(1,:) = (/ 1,2 /)
         EdgeMap(2,:) = (/ 2,3 /)
         EdgeMap(3,:) = (/ 4,3 /)
         EdgeMap(4,:) = (/ 1,4 /)
         EdgeMap(5,:) = (/ 1,5 /)
         EdgeMap(6,:) = (/ 2,5 /)
         EdgeMap(7,:) = (/ 3,5 /)
         EdgeMap(8,:) = (/ 4,5 /)
 
       CASE(7)
         EdgeMap(1,:) = (/ 1,2 /)
         EdgeMap(2,:) = (/ 2,3 /)
         EdgeMap(3,:) = (/ 3,1 /)
         EdgeMap(4,:) = (/ 4,5 /)
         EdgeMap(5,:) = (/ 5,6 /)
         EdgeMap(6,:) = (/ 6,4 /)
         EdgeMap(7,:) = (/ 1,4 /)
         EdgeMap(8,:) = (/ 2,5 /)
         EdgeMap(9,:) = (/ 3,6 /)

       CASE(8)
         EdgeMap(1,:)  = (/ 1,2 /)
         EdgeMap(2,:)  = (/ 2,3 /)
         EdgeMap(3,:)  = (/ 4,3 /)
         EdgeMap(4,:)  = (/ 1,4 /)
         EdgeMap(5,:)  = (/ 5,6 /)
         EdgeMap(6,:)  = (/ 6,7 /)
         EdgeMap(7,:)  = (/ 8,7 /)
         EdgeMap(8,:)  = (/ 5,8 /)
         EdgeMap(9,:)  = (/ 1,5 /)
         EdgeMap(10,:) = (/ 2,6 /)
         EdgeMap(11,:) = (/ 3,7 /)
         EdgeMap(12,:) = (/ 4,8 /)
       END SELECT
     END IF
!------------------------------------------------------------------------------
  END FUNCTION GetEdgeMap
!------------------------------------------------------------------------------

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

!-----------------------------------------------------------------------
!> This routine may be used to terminate the program in the case of an error.
!-----------------------------------------------------------------------
   SUBROUTINE Assert(Condition, Caller, ErrorMessage)
!-----------------------------------------------------------------------
     CHARACTER(LEN=*), OPTIONAL :: Caller, ErrorMessage
     LOGICAL :: Condition
!-----------------------------------------------------------------------
     IF ( .NOT. OutputLevelMask(0) ) STOP

     IF(Condition) RETURN !Assertion passed

     WRITE( Message, '(A)') 'ASSERTION ERROR'

     IF(PRESENT(Caller)) THEN
       WRITE( Message, '(A,A,A)') TRIM(Message),': ',TRIM(Caller)
     END IF

     IF(PRESENT(ErrorMessage)) THEN
       WRITE( Message, '(A,A,A)') TRIM(Message),': ',TRIM(ErrorMessage)
     END IF

     WRITE( *, '(A)', ADVANCE='YES' ) Message

     !Provide a stack trace if no caller info provided
#ifdef __GFORTRAN__
     IF(.NOT.PRESENT(Caller)) CALL BACKTRACE
#endif

     STOP

     CALL FLUSH(6)
!-----------------------------------------------------------------------
   END SUBROUTINE Assert
!-----------------------------------------------------------------------

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

    CALL Info('GetNOFColours','Number of colours: '//TRIM(I2S(ncolours)),Level=12)
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

    CALL Info('GetNOFBoundaryColours','Number of colours: '//TRIM(I2S(ncolours)),Level=12)
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

