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

   INTEGER, ALLOCATABLE, TARGET, PRIVATE :: IndexStore(:)
   REAL(KIND=dp), ALLOCATABLE, TARGET, PRIVATE  :: Store(:)
   ! TODO: Get actual values for these from mesh
   INTEGER, PARAMETER, PRIVATE :: ISTORE_MAX_SIZE = 1024
   INTEGER, PARAMETER, PRIVATE :: STORE_MAX_SIZE = 1024
   ! SAVE IndexStore, Store

   !$OMP THREADPRIVATE(IndexStore, Store)
   PRIVATE :: GetIndexStore, GetStore

CONTAINS

!
!  FUNCTION GetIndexStore() RESULT(ind)
!    INTEGER, POINTER :: Ind(:)
!    INTEGER :: thread, nthreads, istat
!    INTEGER :: omp_get_max_threads, omp_get_thread_num
!
!    IF ( .NOT.ALLOCATED(IndexStore) ) THEN
! !$omp barrier
! !$omp critical(get_index)
!      IF ( .NOT.ALLOCATED(IndexStore) ) THEN
!        nthreads = 1
! !$      nthreads = omp_get_max_threads()
!        ALLOCATE( IndexStore(nthreads,512), STAT=istat )
!        IF ( Istat /= 0 ) &
!           CALL Fatal( 'GetIndexStore', 'Memory allocation error.' )
!      END IF
! !$omp end critical(get_index)
!    END IF
!
!10  thread = 1
! !$  thread=omp_get_thread_num()+1
!    ind => IndexStore( thread, : )
!  END FUNCTION GetIndexStore


!  FUNCTION GetStore(n) RESULT(val)
!    REAL(KIND=dp), POINTER :: val(:)
!    INTEGER :: n,thread, nthreads, istat
!    INTEGER :: omp_get_max_threads, omp_get_thread_num
!
!    IF ( .NOT.ALLOCATED(Store) ) THEN
! !$omp barrier
! !$omp critical(get_store)
!       IF ( .NOT.ALLOCATED(Store) ) THEN
!         nthreads = 1
! !$      nthreads = omp_get_max_threads()
!         ALLOCATE( Store(nthreads*MAX_ELEMENT_NODES), STAT=istat )
!         IF ( Istat /= 0 ) &
!            CALL Fatal( 'GetStore', 'Memory allocation error.' )
!       END IF
! !$omp end critical(get_store)
!      END IF
!
!      thread = 0
! !$   thread=omp_get_thread_num()
!      val => Store( thread*MAX_ELEMENT_NODES+1:thread*MAX_ELEMENT_NODES+n )
!   END FUNCTION GetStore

  FUNCTION GetIndexStore() RESULT(ind)
    INTEGER, POINTER :: Ind(:)
    INTEGER :: istat

    IF ( .NOT. ALLOCATED(IndexStore) ) THEN
        ALLOCATE( IndexStore(ISTORE_MAX_SIZE), STAT=istat )
        IndexStore = 0
        IF ( Istat /= 0 ) CALL Fatal( 'GetIndexStore', 'Memory allocation error.' )
    END IF

    ind => IndexStore( : )
  END FUNCTION GetIndexStore

  FUNCTION GetStore(n) RESULT(val)
    REAL(KIND=dp), POINTER :: val(:)
    INTEGER :: n, istat

    IF ( .NOT.ALLOCATED(Store) ) THEN
        ALLOCATE( Store(STORE_MAX_SIZE), STAT=istat )
        Store = 0D0
        IF ( Istat /= 0 ) CALL Fatal( 'GetStore', 'Memory allocation error.' )
    END IF

    IF (n > STORE_MAX_SIZE) THEN
        CALL Fatal( 'GetStore', 'Not enough memory allocated for store.' )
    END IF

    val => Store( 1:n )
  END FUNCTION GetStore


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
    TYPE(Element_t), POINTER :: Ret_Element
    TYPE(Element_t), OPTIONAL, TARGET :: Element

    IF (PRESENT(Element)) THEN
      Ret_Element=>Element
    ELSE
      Ret_Element=>CurrentModel % CurrentElement
    END IF
  END FUNCTION GetCurrentElement

!> Returns handle to the index of the current element
  FUNCTION GetElementIndex(Element) RESULT(Indx)
     TYPE(Element_t), OPTIONAL :: Element
     INTEGER :: Indx
     TYPE(Element_t), POINTER :: CurrElement

     CurrElement => GetCurrentElement(Element)
     Indx = CurrElement % ElementIndex
  END FUNCTION GetElementIndex

!> Returns the number of active elements for the current solver
  FUNCTION GetNOFActive( USolver ) RESULT(n)
     INTEGER :: n
     TYPE(Solver_t), OPTIONAL, TARGET :: USolver

     IF ( PRESENT( USolver ) ) THEN
        n = USolver % NumberOfActiveElements
     ELSE
        n = CurrentModel % Solver % NumberOfActiveElements
     END IF
  END FUNCTION GetNOFActive

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

     Indexes => GetIndexStore()
     IF ( ASSOCIATED(Variable % Solver) ) THEN
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
      Output,Secondary,Global,InitValue,USolver,Var )
    
    CHARACTER(LEN=*) :: Name
    INTEGER, OPTIONAL :: DOFs
    REAL(KIND=dp), OPTIONAL, POINTER :: Values(:)
    LOGICAL, OPTIONAL :: Output
    INTEGER, OPTIONAL, POINTER :: Perm(:)
    LOGICAL, OPTIONAL :: Secondary
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
        Perm,Output,Secondary,Global,InitValue )

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
     REAL(KIND=dp), POINTER :: x(:)
     TYPE(Element_t), POINTER :: Element

     INTEGER :: n, nthreads, thread, istat

     IF ( PRESENT( Found ) ) Found = .FALSE.

     NodeIndexes => Dnodes
     n = 1
     NodeIndexes(n) = 1

     x => GetStore(n)
     x = 0.0d0
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
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found
     INTEGER, TARGET :: Dnodes(1)
     INTEGER, POINTER :: NodeIndexes(:)
     TYPE(Element_t), OPTIONAL, TARGET :: UElement

     REAL(KIND=dp), POINTER :: x(:)
     TYPE(Element_t), POINTER :: Element

     INTEGER :: n, nthreads, thread, istat

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

     x => GetStore(n)
     x = 0.0_dp
     IF( ASSOCIATED(List) ) THEN
       IF ( ASSOCIATED(List % Head) ) THEN
         x(1:n) = ListGetReal( List, Name, n, NodeIndexes, Found )
       END IF
     END IF
  END FUNCTION GetReal


!> Returns a material property from either of the parents of the current boundary element
  RECURSIVE FUNCTION GetParentMatProp( Name, UElement, Found, UParent ) RESULT(x)
    CHARACTER(LEN=*) :: Name
    TYPE(Element_t), OPTIONAL, TARGET :: UElement
    LOGICAL, OPTIONAL :: Found
    TYPE(Element_t), OPTIONAL, POINTER :: UParent

    REAL(KIND=dp), POINTER :: x(:)    
    INTEGER, POINTER :: Indexes(:)
    LOGICAL :: GotIt
    INTEGER :: n, leftright
    TYPE(ValueList_t), POINTER :: Material
    TYPE(Element_t), POINTER :: Element, Parent
    
    Element => GetCurrentElement(Uelement)

    IF( PRESENT(UParent) ) NULLIFY( UParent )

    n = GetElementNOFNodes(Element)
    Indexes => Element % NodeIndexes

    x => GetStore(n)
    x = 0._dp

    Gotit = .FALSE.
    DO leftright = 1, 2
       IF( leftright == 1) THEN
         Parent => Element % BoundaryInfo % Left
       ELSE 
         Parent => Element % BoundaryInfo % Right
       END IF

       IF( ASSOCIATED(Parent) ) THEN
         Material => GetMaterial(Parent)
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


!> Set some property elementwise to the active element
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

!> Get some property elementwise from the active element
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

     Solver => CurrentModel % Solver
     IF ( PRESENT( USolver ) ) Solver => USolver

     IF ( t > 0 .AND. t <= Solver % NumberOfActiveElements ) THEN
        Element => Solver % Mesh % Elements( Solver % ActiveElements(t) )
        !$omp critical(GetActiveElementCurrentElement)
        CurrentModel % CurrentElement => Element ! may be used by user functions
        !$omp end critical(GetActiveElementCurrentElement)
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

     Solver => CurrentModel % Solver
     IF ( PRESENT( USolver ) ) Solver => USolver

     IF ( t > 0 .AND. t <= Solver % Mesh % NumberOfBoundaryElements ) THEN
        Element => Solver % Mesh % Elements( Solver % Mesh % NumberOfBulkElements+t )
        !$omp critical(GetBoundaryElementCurrentElement)
        CurrentModel % CurrentElement => Element ! may be used be user functions
        !$omp end critical(GetBoundaryElementCurrentElement)
     ELSE
        WRITE( Message, * ) 'Invalid element number requested: ', t
        CALL Fatal( 'GetBoundaryElement', Message )
     END IF
  END FUNCTION GetBoundaryElement


!> Check if the boundary element is active in the current solve 
  FUNCTION ActiveBoundaryElement(UElement,USolver) RESULT(l)
     TYPE(Element_t), OPTIONAL,  TARGET :: UElement
     TYPE(Solver_t),  OPTIONAL,  TARGET :: USolver

     LOGICAL :: l
     INTEGER :: n
     INTEGER, POINTER :: Indexes(:)

     TYPE(Element_t), POINTER :: Element
     TYPE( Solver_t ), POINTER :: Solver

     Solver => CurrentModel % Solver
     IF ( PRESENT( USolver ) ) Solver => USolver

     Element => GetCurrentElement(UElement)

     Indexes => GetIndexStore()
     n = GetElementDOFs( Indexes, Element, Solver )
     IF (isActivePElement(Element)) n=GetElementNOFNOdes(Element)

     l = ALL( Solver % Variable % Perm(Indexes(1:n)) > 0)
  END FUNCTION ActiveBoundaryElement


!> Return the element code in Elmer convention of the active element
  FUNCTION GetElementCode( Element )  RESULT(etype)
     INTEGER :: etype
     TYPE(Element_t), OPTIONAL :: Element
     TYPE(Element_t), POINTER :: CurrElement

     CurrElement => GetCurrentElement(Element)
     etype = CurrElement % TYPE % ElementCode
  END FUNCTION GetElementCode


!> Return the element family in Elmer convention of the active element
  FUNCTION GetElementFamily( Element )  RESULT(family)
     INTEGER :: family
     TYPE(Element_t), OPTIONAL :: Element
     TYPE(Element_t), POINTER :: CurrElement

     CurrElement => GetCurrentElement(Element)
     family = CurrElement % TYPE % ElementCode / 100
  END FUNCTION GetElementFamily

!> Return true if the element is a possible flux element
!> Needed to skip nodal elements in 2D and 3D boundary condition setting.
  FUNCTION PossibleFluxElement( Element, Mesh )  RESULT(possible)
     LOGICAL :: possible
     TYPE(Element_t), OPTIONAL :: Element
     TYPE(Mesh_t), OPTIONAL :: Mesh
     INTEGER :: MeshDim, family

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

     INTEGER :: i,j, id
     LOGICAL :: Found, GB

     IF ( PRESENT( UElement ) ) THEN
        Element => UElement
     ELSE
        Element => CurrentModel % CurrentElement
     END IF

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
     ELSE
        Solver => CurrentModel % Solver
     END IF

     n = 0
     IF ( ListGetLogical( Solver % Values, 'Discontinuous Galerkin', Found )) THEN
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

     IF ( Solver % Def_Dofs(GetElementFamily(Element),id,1)>0 ) n = Element % NDOFs
     IF ( ALL(Solver % Def_Dofs(GetElementFamily(Element),id,2:)<0) ) RETURN

     IF ( ASSOCIATED( Element % EdgeIndexes ) ) THEN
        DO j=1,Element % TYPE % NumberOFEdges
           n =  n + Solver % Mesh % Edges(Element % EdgeIndexes(j)) % BDOFs
        END DO
     END IF

     IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
        DO j=1,Element % TYPE % NumberOFFaces
           n = n + Solver % Mesh % Faces( Element % FaceIndexes(j) ) % BDOFs
        END DO
     END IF

     GB = ListGetLogical( Solver % Values, 'Bubbles in Global System', Found )
     IF (.NOT.Found) GB = .TRUE.
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

     LOGICAL :: Found, GB, DGdisable
     INTEGER :: nb,i,j,k,id,EDOFs, FDOFs, BDOFs,FaceDOFs, EdgeDOFs, BubbleDOFs, Ind

     Element => GetCurrentElement(UElement)

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
     ELSE
        Solver => CurrentModel % Solver
     END IF

     NB = 0

     DGDisable=.FALSE.
     IF (PRESENT(NotDG)) DGDisable=NotDG

     IF ( .NOT.DGDisable .AND. ListGetLogical( Solver % Values, 'Discontinuous Galerkin', Found ) ) THEN
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

     IF ( Solver % Def_Dofs(GetElementFamily(Element),id,1)>0 ) THEN
       DO i=1,Element % NDOFs
          NB = NB + 1
          Indexes(NB) = Element % NodeIndexes(i)
       END DO
     END IF

     ! default for nodal elements, if no solver active:
     ! ------------------------------------------------
     IF(.NOT.ASSOCIATED(Solver)) RETURN
     IF(.NOT.ASSOCIATED(Solver % Mesh)) RETURN

     IF ( ALL(Solver % Def_Dofs(GetElementFamily(Element),id,2:)<0) ) RETURN

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

     GB = ListGetLogical( Solver % Values, 'Bubbles in Global System', Found )
     IF (.NOT.Found) GB = .TRUE.

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

    GB = ListGetLogical( Solver % Values, 'Bubbles in Global System', Found )
    IF (.NOT.Found) GB = .TRUE.


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

	 IF ( PRESENT(UElement) ) THEN
	 	Element => UElement
	 ELSE
	 	Element => CurrentModel % CurrentElement
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

     IF ( PRESENT(Uelement) ) THEN
        Element => UElement
     ELSE
       Element => CurrentModel % CurrentElement
     END IF

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

     INTEGER :: ic_id, body_id

     IF ( PRESENT( Element ) ) THEN
        body_id = Element % BodyId 
     ELSE
        body_id = CurrentModel % CurrentElement % BodyId 
     END IF

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

    IF ( PRESENT(UElement) ) THEN
      Element => UElement
    ELSE
      Element => CurrentModel % CurrentElement
    END IF

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

    IF ( PRESENT(UElement) ) THEN
      Element => UElement
    ELSE
      Element => CurrentModel % CurrentElement
    END IF

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

    Element => CurrentModel % CurrentElement
    IF ( PRESENT(UElement) ) Element => UElement

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

    Element => CurrentModel % CurrentElement
    IF ( PRESENT(UElement) ) Element => UElement

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




!> Performs initialization for matrix equation related to the active solver
!------------------------------------------------------------------------------
  SUBROUTINE DefaultInitialize( USolver )
!------------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL, TARGET, INTENT(IN) :: USolver

     TYPE(Solver_t), POINTER :: Solver, SlaveSolver
     TYPE(ValueList_t), POINTER :: Params
     TYPE(Variable_t), POINTER :: iterV
     INTEGER, POINTER :: SlaveSolverIndexes(:)
     INTEGER :: j,k,iter
     REAL(KIND=dp) :: dt
     LOGICAL :: Transient, Found, alloc_parenv

     INTERFACE
        SUBROUTINE SolverActivate_x(Model,Solver,dt,Transient)
          USE Types
          TYPE(Model_t)::Model
          TYPE(Solver_t),POINTER::Solver
          REAL(KIND=dp) :: dt
          LOGICAL :: Transient
        END SUBROUTINE SolverActivate_x
     END INTERFACE

     IF ( PRESENT( USolver ) ) THEN
       Solver => USolver
     ELSE
       Solver => CurrentModel % Solver
     END IF

     ! One can enforce weak coupling by calling a dependent solver at the start of 
     ! initialization of the master solver. The strategy will be particularly 
     ! efficient when the slave solver is cheap and a stepsize control is applied
     ! to the master solver.
     !-----------------------------------------------------------------------------
     SlaveSolverIndexes =>  ListGetIntegerArray( Solver % Values,'Slave Solvers',Found )
     IF( Found ) THEN
       dt = GetTimeStep()
       Transient = GetString(CurrentModel % Simulation,'Simulation type',Found)=='transient'

       ! store the nonlinear iteration at the outer loop
       iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
       iter = NINT(iterV % Values(1))


       DO j=1,SIZE(SlaveSolverIndexes)
         k = SlaveSolverIndexes(j)
         SlaveSolver => CurrentModel % Solvers(k)

         IF(ParEnv % PEs>1) THEN
           IF ( Solver % Matrix % Comm /= MPI_COMM_WORLD ) &
             CALL ListAddLogical( SlaveSolver % Values, 'Slave not parallel', .TRUE.)

           alloc_parenv = .FALSE.
           IF(ASSOCIATED(SlaveSolver % Matrix)) THEN
             IF(ASSOCIATED(SlaveSolver % Matrix % ParMatrix) ) THEN
               ParEnv = SlaveSolver % Matrix % ParMatrix % ParEnv
             ELSE
               ALLOCATE(ParEnv % Active(ParEnv % PEs)); alloc_parenv=.TRUE.
             END IF
           ELSE
             ALLOCATE(ParEnv % Active(ParEnv % PEs)); alloc_parenv=.TRUE.
           END IF
           ParEnv % ActiveComm = Solver % Matrix % Comm
         END IF

         CurrentModel % Solver => SlaveSolver
         CALL SolverActivate_x( CurrentModel,SlaveSolver,dt,Transient)

         IF(ParEnv % PEs>1) THEN
           IF ( Solver % Matrix % Comm /= MPI_COMM_WORLD ) &
             CALL ListAddLogical( SlaveSolver % Values, 'Slave not parallel', .FALSE.)

           IF(alloc_parenv) THEN
             DEALLOCATE(ParEnv % Active)
             ParEnv % Active => NULL()
           END IF

           IF(ASSOCIATED(Solver % Matrix)) THEN
             IF(ASSOCIATED(Solver % Matrix % ParMatrix) ) &
               ParEnv = Solver % Matrix % ParMatrix % ParEnv
           END IF
         END IF
       END DO
       CurrentModel % Solver => Solver
       iterV % Values = iter       
     END IF

     IF(.NOT. ASSOCIATED( Solver % Matrix ) ) THEN
       CALL Fatal('DefaultInitialize','No matrix exists, cannot initialize!')
     END IF

     CALL InitializeToZero( Solver % Matrix, Solver % Matrix % RHS )

!------------------------------------------------------------------------------
  END SUBROUTINE DefaultInitialize
!------------------------------------------------------------------------------



!> Performs finilizing steps related to the the active solver
!------------------------------------------------------------------------------
  SUBROUTINE DefaultFinish( USolver )
!------------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL, TARGET, INTENT(IN) :: USolver

     TYPE(Solver_t), POINTER :: Solver, PostSolver
     TYPE(ValueList_t), POINTER :: Params
     TYPE(Variable_t), POINTER :: iterV
     INTEGER, POINTER :: PostSolverIndexes(:)
     INTEGER :: j,k,iter
     REAL(KIND=dp) :: dt
     LOGICAL :: Transient, Found, alloc_parenv

     INTERFACE
        SUBROUTINE SolverActivate_x(Model,Solver,dt,Transient)
          USE Types
          TYPE(Model_t)::Model
          TYPE(Solver_t),POINTER::Solver
          REAL(KIND=dp) :: dt
          LOGICAL :: Transient
        END SUBROUTINE SolverActivate_x
     END INTERFACE

     IF ( PRESENT( USolver ) ) THEN
       Solver => USolver
     ELSE
       Solver => CurrentModel % Solver
     END IF

     ! One can run postprocessing solver in this slot.
     !-----------------------------------------------------------------------------
     PostSolverIndexes =>  ListGetIntegerArray( Solver % Values,'Post Solvers',Found )
     IF( Found ) THEN
       dt = GetTimeStep()
       Transient = GetString(CurrentModel % Simulation,'Simulation type',Found)=='transient'

       ! store the nonlinear iteration at the outer loop
       iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
       iter = NINT(iterV % Values(1))

       DO j=1,SIZE(PostSolverIndexes)
         k = PostSolverIndexes(j)
         PostSolver => CurrentModel % Solvers(k)

         IF(ParEnv % PEs>1) THEN
           IF ( Solver % Matrix % Comm /= MPI_COMM_WORLD ) &
             CALL ListAddLogical( PostSolver % Values, 'Post not parallel', .TRUE.)

           alloc_parenv = .FALSE.
           IF(ASSOCIATED(PostSolver % Matrix)) THEN
             IF(ASSOCIATED(PostSolver % Matrix % ParMatrix) ) THEN
               ParEnv = PostSolver % Matrix % ParMatrix % ParEnv
             ELSE
               ALLOCATE(ParEnv % Active(ParEnv % PEs)); alloc_parenv=.TRUE.
             END IF
           ELSE
             ALLOCATE(ParEnv % Active(ParEnv % PEs)); alloc_parenv=.TRUE.
           END IF
         END IF

         CurrentModel % Solver => PostSolver
         CALL SolverActivate_x( CurrentModel,PostSolver,dt,Transient)

         IF(ParEnv % PEs>1) THEN
           IF ( Solver % Matrix % Comm /= MPI_COMM_WORLD ) &
             CALL ListAddLogical( PostSolver % Values, 'Post not parallel', .FALSE.)

           IF(alloc_parenv) THEN
             DEALLOCATE(ParEnv % Active)
             ParEnv % Active => NULL()
           END IF

           IF(ASSOCIATED(Solver % Matrix)) THEN
             IF(ASSOCIATED(Solver % Matrix % ParMatrix) ) &
               ParEnv = Solver % Matrix % ParMatrix % ParEnv
           END IF
         END IF
       END DO
       CurrentModel % Solver => Solver
       iterV % Values = iter       
     END IF

     CALL Info('DefaultFinish','Finished solver: '//&
         TRIM(ListGetString(Solver % Values,'Equation')),Level=5)

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

    Solver => CurrentModel % Solver
    Norm = REAL(0, dp)
    IF ( PRESENT( USolver ) ) Solver => USolver
    IF ( GetLogical(Solver % Values,'Linear System Solver Disabled',Found) ) RETURN

    A => Solver % Matrix
    b => A % RHS
    x => Solver % Variable
    SOL => x % Values

    Params => GetSolverParams(Solver)

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

    ! Combine the individual projectors into one massive projector
    IF(.NOT.ASSOCIATED(Solver % Matrix % ConstraintMatrix)) &
      Solver % MortarBCsOnly = .TRUE.
    Ctmp => Solver % Matrix % ConstraintMatrix
    CALL GenerateConstraintMatrix( CurrentModel, Solver )

    CALL SolveSystem(A,ParMatrix,b,SOL,x % Norm,x % DOFs,Solver)

    IF(.NOT. Solver % MortarBCsOnly) THEN
      IF(ASSOCIATED(Ctmp).OR.ASSOCIATED(Solver % Matrix % ConstraintMatrix)) THEN
        IF(.NOT.ASSOCIATED(Ctmp, Solver % Matrix % ConstraintMatrix)) THEN
          CALL FreeMatrix(Solver % Matrix % ConstraintMatrix)
          Solver % Matrix % ConstraintMatrix => Ctmp
          IF (ASSOCIATED(Solver % MortarBCs)) Solver % MortarBCsChanged = .TRUE.
        END IF
      END IF
    END IF

    ! If flux corrected transport is used then apply the corrector to the system
    IF( GetLogical( Params,'Linear System FCT',Found ) ) THEN
      CALL FCT_Correction( Solver )
    END IF
 
    IF (PRESENT(BackRotNT)) THEN
      IF (BackRot.NEQV.BackRotNT) &
        CALL ListAddLogical(Params,'Back Rotate N-T Solution',BackRot)
    END IF

    Norm = x % Norm

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
  SUBROUTINE DefaultUpdateEquationsR( G, F, UElement, USolver, BulkUpdate ) 
!------------------------------------------------------------------------------
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     REAL(KIND=dp)   :: G(:,:), f(:)
     LOGICAL, OPTIONAL :: BulkUpdate

     TYPE(Solver_t), POINTER   :: Solver
     TYPE(Matrix_t), POINTER   :: A
     TYPE(Variable_t), POINTER :: x
     TYPE(Element_t), POINTER  :: Element, P1, P2
     REAL(KIND=dp), POINTER    :: b(:), SaveValues(:)

     CHARACTER(LEN=MAX_NAME_LEN) :: str

     LOGICAL :: Found, BUpd

     INTEGER :: n, nd
     INTEGER(KIND=AddrInt) :: Proc
     INTEGER, POINTER :: Indexes(:)

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

     IF ( PRESENT( UElement ) ) THEN
        Element => UElement 
     ELSE
        Element => CurrentModel % CurrentElement
     END IF

     IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
       str = ListGetString( Solver % Values, 'Boundary Element Procedure', Found )
     ELSE
       str = ListGetString( Solver % Values, 'Bulk Element Procedure', Found )
     END IF

     IF ( Found ) THEN
       Proc = GetProcAddr( str, abort=.FALSE.,quiet=.TRUE. )
       IF ( Proc /= 0 ) THEN
         n  = GetElementNOFNodes( Element )
         nd = GetElementNOFDOFs( Element, Solver )
         CALL ExecLocalProc( Proc, CurrentModel, Solver, &
                G, F, Element, n, nd )
       END IF
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

     Indexes => GetIndexStore()
     n = GetElementDOFs( Indexes, Element, Solver )

     CALL UpdateGlobalEquations( A,G,b,f,n,x % DOFs,x % Perm(Indexes(1:n)), UElement=Element )
!------------------------------------------------------------------------------
  END SUBROUTINE DefaultUpdateEquationsR
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE DefaultUpdateEquationsC( GC, FC, UElement, USolver, BulkUpdate ) 
!------------------------------------------------------------------------------
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     COMPLEX(KIND=dp)   :: GC(:,:), FC(:)
     LOGICAL, OPTIONAL :: BulkUpdate

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

     IF ( PRESENT( UElement ) ) THEN
        Element => UElement 
     ELSE
        Element => CurrentModel % CurrentElement
     END IF

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

    Element => CurrentModel % CurrentElement
    IF ( PRESENT(UElement) ) Element => UElement

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

    Element => CurrentModel % CurrentElement
    IF ( PRESENT(UElement) ) Element => UElement

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

     Element => CurrentModel % CurrentElement
     IF ( PRESENT(UElement) ) Element => UElement

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

      Element => CurrentModel % CurrentElement
      IF ( PRESENT(UElement) ) Element => UElement

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

     IF ( PRESENT( UElement ) ) THEN
        Element => UElement 
     ELSE
        Element => CurrentModel % CurrentElement
     END IF

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

     IF ( PRESENT( UElement ) ) THEN
        Element => UElement 
     ELSE
        Element => CurrentModel % CurrentElement
     END IF

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

     IF ( PRESENT( UElement ) ) THEN
        Element => UElement 
     ELSE
        Element => CurrentModel % CurrentElement
     END IF

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

     IF ( PRESENT( UElement ) ) THEN
        Element => UElement 
     ELSE
        Element => CurrentModel % CurrentElement
     END IF

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

     IF ( PRESENT( UElement ) ) THEN
        Element => UElement 
     ELSE
        Element => CurrentModel % CurrentElement
     END IF

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

     IF ( PRESENT( UElement ) ) THEN
        Element => UElement 
     ELSE
        Element => CurrentModel % CurrentElement
     END IF

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

     IF ( PRESENT( UElement ) ) THEN
        Element => UElement 
     ELSE
        Element => CurrentModel % CurrentElement
     END IF

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

     IF ( PRESENT( UElement ) ) THEN
        Element => UElement 
     ELSE
        Element => CurrentModel % CurrentElement
     END IF

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
     REAL(KIND=dp), POINTER    :: b(:)

     REAL(KIND=dp) :: xx, s
     REAL(KIND=dp), POINTER :: DiagScaling(:)
     REAL(KIND=dp), ALLOCATABLE :: Work(:), STIFF(:,:)

     INTEGER, ALLOCATABLE :: lInd(:), gInd(:)
     INTEGER :: i,j, k, kk, l, m, n,nd, nb, mb, nn, ni, nj, &
          EDOFs, DOF, local, numEdgeDofs, istat, n_start, Offset

     LOGICAL :: Flag,Found, ConstantValue, ScaleSystem
     TYPE(ValueListEntry_t), POINTER :: ptr
     TYPE(ValueList_t), POINTER :: BC, Params
     TYPE(Element_t), POINTER :: Element, Parent, Edge, Face, SaveElement
     CHARACTER(LEN=MAX_NAME_LEN) :: name
     LOGICAL :: BUpd, PiolaTransform

     INTEGER::iii=0

     SAVE gInd, lInd, STIFF, Work

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

     IF(.NOT.ALLOCATED(A % ConstrainedDOF)) &
       ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
     A % ConstrainedDOF = .FALSE.

     IF(.NOT.ALLOCATED(A % Dvalues)) ALLOCATE(A % Dvalues(A % NumberOfRows))
     A % Dvalues = 0._dp

     ScaleSystem=GetLogical(Params,'Linear System Dirichlet Scaling',Found)
     IF(.NOT.Found) THEN
       ScaleSystem=GetLogical(Params,'Linear System Scaling',Found)
       IF(.NOT.Found) ScaleSystem=.TRUE.
     END IF

     IF (ScaleSystem) THEN
       CALL ScaleLinearSystem(Solver,A,b,RHSscaling=.FALSE.)
       DiagScaling => A % DiagScaling
     ELSE
       ALLOCATE(DiagScaling(A % NumberOfRows))
       DiagScaling=1._dp
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
        
        ! Clearing for p-approximation dofs associated with faces & edges: 
        SaveElement => CurrentModel % CurrentElement
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

           ConstantValue =  ptr % PROCEDURE == 0 .AND. &
             ptr % TYPE == LIST_TYPE_CONSTANT_SCALAR

           IF ( isActivePElement(Parent)) THEN
              n = GetElementNOFNodes()
              ! Get indexes of boundary dofs:
              CALL getBoundaryIndexes( Solver % Mesh, Element, Parent, gInd, numEdgeDofs )

              DO k=n+1,numEdgeDofs
                 nb = x % Perm( gInd(k) )
                 IF ( nb <= 0 ) CYCLE
                 nb = Offset + x % DOFs * (nb-1) + DOF
                 IF ( ConstantValue  ) THEN
                   IF (A % NoDirichlet) THEN
                     A % ConstrainedDOF(nb) = .TRUE.
                     A % Dvalues(nb) = 0._dp
                   ELSE
                     CALL CRS_SetSymmDirichlet(A, A % RHS, nb, 0._dp )
                   END IF
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
        CurrentModel % CurrentElement => SaveElement
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
        SaveElement => CurrentModel % CurrentElement
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
           ConstantValue =  ptr % PROCEDURE == 0 .AND. &
                ptr % TYPE == LIST_TYPE_CONSTANT_SCALAR
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
                    IF(A % ConstrainedDOF(nb)) THEN
                      s = A % Dvalues(nb)
                    ELSE
                      s = A % RHS(nb)
                    END IF
                    s = s * DiagScaling(nb)
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
                    IF( A % NoDirichlet ) THEN
                       A % ConstrainedDOF(nb) = .TRUE.
                       A % Dvalues(nb) = Work(k-n)/DiagScaling(nb)
                    ELSE
                       CALL CRS_SetSymmDirichlet(A,A % RHS,nb,Work(k-n)/DiagScaling(nb))
                    END IF
                 END DO
              ELSE
                 ! Contribute this boundary to global system
                 ! (i.e solve global boundary problem)
                 DO k=n+1,numEdgeDofs
                    nb = x % Perm( gInd(k) )
                    IF ( nb <= 0 ) CYCLE
                    nb = Offset + x % DOFs * (nb-1) + DOF
                    A % RHS(nb) = A % RHS(nb) + Work(k)/DiagScaling(nb)
                    DO l=1,numEdgeDofs
                       mb = x % Perm( gInd(l) )
                       IF ( mb <= 0 ) CYCLE
                       mb = Offset + x % DOFs * (mb-1) + DOF
                       DO kk=A % Rows(nb)+DOF-1,A % Rows(nb+1)-1,x % DOFs
                          IF ( A % Cols(kk) == mb ) THEN
                             A % Values(kk) = A % Values(kk) + STIFF(k,l) * &
                                  DiagScaling(mb) / DiagScaling(nb)
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
                    IF(A % ConstrainedDOF(nb)) THEN
                      s = A % Dvalues(nb)
                    ELSE
                      s = A % RHS(nb)
                    END IF
                    s = s * DiagScaling(nb)
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
                 A % RHS(nb) = A % RHS(nb) + Work(k)/DiagScaling(nb)
                 DO l=n_start,numEdgeDOFs
                    mb = x % Perm( gInd(l) )
                    IF ( mb <= 0 ) CYCLE
                    mb = Offset + x % DOFs * (mb-1) + DOF
                    DO kk=A % Rows(nb)+DOF-1,A % Rows(nb+1)-1,x % DOFs
                       IF ( A % Cols(kk) == mb ) THEN
                          A % Values(kk) = A % Values(kk) + STIFF(k,l) * &
                               DiagScaling(mb) / DiagScaling(nb)
                          EXIT
                       END IF
                    END DO
                 END DO
              END DO
           END SELECT
        END DO
        CurrentModel % CurrentElement => SaveElement
     END DO

     ! ----------------------------------------------------------------------------
     ! Set Dirichlet BCs for edge and face dofs which arise from approximating with
     ! edge (curl-conforming) or face (div-conforming) elements:
     ! ----------------------------------------------------------------------------
     DO DOF=1,x % DOFs
        name = x % name
        IF (x % DOFs>1) name=ComponentName(name,DOF)

        SaveElement => CurrentModel % CurrentElement
        DO i=1,Solver % Mesh % NumberOfBoundaryElements
           Element => GetBoundaryElement(i)
           IF ( .NOT. ActiveBoundaryElement() ) CYCLE

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

                    EDOFs = Edge % BDOFs
                    IF (EDOFs == 1) THEN
                       !-----------------------------------------------------------------
                       ! The lowest-order edge element interpolation with one DOF/edge:
                       !-----------------------------------------------------------------
                       nb = Parent % TYPE % NumberOfNodes
                       n  =  Edge % TYPE % NumberOfNodes
                       CALL LocalBcIntegral(BC,Edge,n,Parent,nb,TRIM(Name)//' {e}',Work(1))

                       n=GetElementDOFs(gInd,Edge)

                       n_start = Solver % Def_Dofs(2,Parent % BodyId,1)*Edge % NDOFs
                       DO k=n_start+1,n_start+EDOFs
                          nb = x % Perm(gInd(k))
                          IF ( nb <= 0 ) CYCLE
                          nb = Offset + x % DOFs*(nb-1) + DOF
                          IF ( A % Symmetric.AND..NOT.A % NoDirichlet ) THEN
                             CALL CRS_SetSymmDirichlet(A,A % RHS,nb,Work(1)/DiagScaling(nb))
                          ELSE
                             A % ConstrainedDOF(nb) = .TRUE.
                             A % Dvalues(nb) = Work(1)/DiagScaling(nb)
                             IF( A % NoDirichlet ) THEN
                             ELSE
                               CALL ZeroRow( A, nb )
!                              A % RHS(nb) = Work(1)/DiagScaling(nb)
                               CALL SetMatrixElement(A,nb,nb,1._dp)
                             END IF
                          END IF
                       END DO
                    ELSE
                       !-----------------------------------------------------------------
                       ! The cases with more than one DOF/edge. To do: handle these by
                       ! L2 projection of boundary data
                       !-----------------------------------------------------------------                    
                       CALL Fatal('DefaultDirichletBCs',&
                            'BCs are not yet defined for this type of edge element interpolation')
                    END IF
                 CASE(3,4)
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

                    DO j=1,Face % TYPE % NumberOfEdges
                       Edge => Solver % Mesh % Edges(Face % EdgeIndexes(j))
                       EDOFs = Edge % BDOFs

                       IF (EDOFs == 1) THEN
                          !-----------------------------------------------------------------------
                          ! Handle the lowest-order edge element interpolation with one DOF/edge:
                          !-----------------------------------------------------------------------                      
                          nb = Edge % TYPE % NumberOfNodes
                          n  = Parent % TYPE % NumberOfNodes

                          CALL LocalBcIntegral( BC, Edge, nb, Parent, &
                               n, TRIM(Name)//' {e}', Work(1) )

                          n=GetElementDOFs(gInd,Edge)

                          n_start = Solver % Def_Dofs(2,Parent % BodyId,1)*Edge % NDOFs
                          DO k=n_start+1,n_start+EDOFs
                             nb = x % Perm(gInd(k))
                             IF ( nb <= 0 ) CYCLE
                             nb = Offset + x % DOFs*(nb-1) + DOF
                             IF ( A % Symmetric .AND..NOT. A % NoDirichlet ) THEN
                                CALL CRS_SetSymmDirichlet(A,A % RHS,nb,Work(1)/DiagScaling(nb))
                             ELSE
                                A % ConstrainedDOF(nb) = .TRUE.
                                A % Dvalues(nb) = Work(1)/DiagScaling(nb)
                                IF(A % NoDirichlet ) THEN
                                ELSE
                                  CALL ZeroRow(A,nb)
                                  CALL SetMatrixElement(A,nb,nb,1._dp)
!                                 A % RHS(nb) = Work(1)/DiagScaling(nb)
                                END IF
                             END IF
                          END DO
                       ELSE
                          ! --------------------------------------------------------------------
                          ! Functionality missing! Faces with more than one DOF/edge or
                          ! with a combination of edge and face DOFs should be handled by
                          ! the L2 projection of data.
                          !-----------------------------------------------------------------             
                          ! Currently set zero BC so that this has some restricted utility 
                          ! in connection with the AV-formulation. 
                          ! --------------------------------------------------------------------
                          n=GetElementDOFs(gInd,Edge)

                          n_start = Solver % Def_Dofs(2,Parent % BodyId,1)*Edge % NDOFs
                          DO k=n_start+1,n_start+EDOFs
                             nb = x % Perm(gInd(k))
                             IF ( nb <= 0 ) CYCLE
                             nb = Offset + x % DOFs*(nb-1) + DOF
                             IF ( A % Symmetric.AND..NOT.A % NoDirichlet ) THEN
                                CALL CRS_SetSymmDirichlet(A,A % RHS,nb,0.0d0)
                             ELSE
                                A % ConstrainedDOF(nb) = .TRUE.
                                A % Dvalues(nb) = 0._dp
                                IF( A % NoDirichlet ) THEN
                                ELSE
                                  CALL ZeroRow(A,nb)
!                                 A % RHS(nb) = 0.0d0
                                  CALL SetMatrixElement(A,nb,nb,1._dp)
                                END IF
                             END IF
                          END DO
                       END IF
                    END DO

                    ! ---------------------------------------------------------------------
                    ! Set constraints for face DOFs: Functionality missing!
                    ! ---------------------------------------------------------------------
                    ! Currently set zero BC so that this has some restricted utility 
                    ! in connection with the AV-formulation:
                    ! --------------------------------------------------------------------
                    IF (Face % BDOFs > 0) THEN
                       n = GetElementDOFs(GInd,Face)
                       DO j=1,Face % BDOFs
                          nb = x % Perm(GInd(n-Face % BDOFs+j)) ! The last entries should be face-DOF indices
                          IF ( nb <= 0 ) CYCLE
                          nb = Offset + x % DOFs*(nb-1) + DOF
                          IF ( A % Symmetric.AND..NOT.A % NoDirichlet ) THEN
                             CALL CRS_SetSymmDirichlet(A,A % RHS,nb,0.0d0)
                          ELSE
                             A % ConstrainedDOF(nb) = .TRUE.
                             A % Dvalues(nb) = 0.0d0
                             IF(A % NoDirichlet ) THEN
                             ELSE
                               CALL ZeroRow(A,nb)
!                              A % RHS(nb) = 0.0d0
                               CALL SetMatrixElement(A,nb,nb,1._dp)
                             END IF
                          END IF
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
        CurrentModel % CurrentElement => SaveElement
     END DO

     Found = .NOT. A % NoDirichlet
     IF ( Found ) THEN
        DO k=1,A % NumberOfRows
          IF ( A % ConstrainedDOF(k) ) THEN
            s = A % Values(A % Diag(k))
            IF (s==0) s = 1

            IF ( A % Symmetric ) THEN
              CALL CRS_SetSymmDirichlet(A,b,k,A % Dvalues(k)/s)
            ELSE
              CALL ZeroRow(A, k)
              b(k) = A % Dvalues(k)/s
              CALL SetMatrixElement(A,k,k,1._dp)
            END IF
          END IF
        END DO
        DEALLOCATE(A % Dvalues)
        A % NoDirichlet = .FALSE.
     END IF

     IF (ScaleSystem) THEN
       CALL BackScaleLinearSystem(Solver,A,b)
     ELSE
       DEALLOCATE(DiagScaling)
     END IF

     ! Add the possible constraint modes structures
     !----------------------------------------------------------
     IF ( GetLogical(Solver % Values,'Constraint Modes Analysis',Found) ) THEN
       CALL SetConstraintModesBoundaries( CurrentModel, A, b, x % Name, x % DOFs, x % Perm )
     END IF

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


!> This subroutine can be used to compute the values of DOFs corresponding
!> the edge finite elements of the lowest order so that the edge finite
!> element interpolant of the BC data is obtained. The value of DOF is
!> defined as D = S*(g.t,1)_E where g.t is tangential component of data,
!> E is the region occupied by the boundary element and S reverts sign
!> if necessary.
!------------------------------------------------------------------------------
  SUBROUTINE LocalBcIntegral(BC, Element, nd, Parent, np, Name, Integral)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(ValueList_t), POINTER :: BC !< The list of boundary condition values
    TYPE(Element_t) :: Element       !< The boundary element handled
    INTEGER :: nd                    !< The number of boundary element nodes
    TYPE(Element_t) :: Parent        !< The parent element of the boundary element
    INTEGER :: np                    !< The number of parent element nodes
    CHARACTER(LEN=*) :: Name         !< The name of boundary condition
    REAL(KIND=dp) :: Integral        !< The value of DOF
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER, POINTER :: Edgemap(:,:)
    INTEGER :: i,j,k,n,p,q,t
    LOGICAL :: Lstat
    TYPE(ElementType_t), POINTER :: SavedType
    TYPE(Nodes_t), SAVE :: Nodes, Pnodes
    REAL(KIND=dp) :: Basis(np) ,dBasisdx(np,3)
    REAL(KIND=dp) :: L,VL(3),G(3)
    REAL(KIND=dp) :: u,v,w,s,DetJ,Load(np),Vload(3,1:np)
!------------------------------------------------------------------------------

    ! Get the nodes of the boundary and parent elements:
    CALL GetElementNodes(Nodes, Element)
    CALL GetElementNodes(PNodes, Parent)

    Load(1:nd) = GetReal( BC, Name, Lstat, Element )

    EdgeMap => GetEdgeMap(GetElementFamily(Parent))
    DO i=1,SIZE(EdgeMap,1)
      j=EdgeMap(i,1); k=EdgeMap(i,2)
      IF ( Parent % NodeIndexes(j)==Element % NodeIndexes(1) .AND. &
           Parent % NodeIndexes(k)==Element % NodeIndexes(2) .OR.  &
           Parent % NodeIndexes(j)==Element % NodeIndexes(2) .AND. &
           Parent % NodeIndexes(k)==Element % NodeIndexes(1) ) EXIT
    END DO

    n = LEN_TRIM(Name)
    VLoad(1,1:nd)=GetReal(BC,Name(1:n)//' 1',Lstat,element)
    VLoad(2,1:nd)=GetReal(BC,Name(1:n)//' 2',Lstat,element)
    VLoad(3,1:nd)=GetReal(BC,Name(1:n)//' 3',Lstat,element)

    G(1) = PNodes % x(k) - PNodes % x(j)
    G(2) = PNodes % y(k) - PNodes % y(j)
    G(3) = PNodes % z(k) - PNodes % z(j)
    G = G/SQRT(SUM(G**2))

    SavedType => Element % TYPE
    IF ( GetElementFamily()==1 ) Element % TYPE=>GetElementType(202)
      
    Integral = 0._dp
    IP = GaussPoints(Element)
    DO t=1,IP % n
      Lstat = ElementInfo( Element, Nodes, IP % u(t), &
            IP % v(t), IP % w(t), DetJ, Basis )
      s = IP % s(t) * DetJ

      L  = SUM(Load(1:nd)*Basis(1:nd))
      VL = MATMUL(Vload(:,1:nd),Basis(1:nd))
      Integral=Integral+s*(L+SUM(VL*G))
    END DO
    Element % TYPE => SavedType

    j = Parent % NodeIndexes(j)
    IF ( ParEnv % PEs>1 ) &
      j=CurrentModel % Mesh % ParallelInfo % GlobalDOFs(j)

    k = Parent % NodeIndexes(k)
    IF ( ParEnv % PEs>1 ) &
      k=CurrentModel % Mesh % ParallelInfo % GlobalDOFs(k)

    IF (k < j) Integral=-Integral
!------------------------------------------------------------------------------
  END SUBROUTINE LocalBcIntegral
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
      CALL CopyBulkMatrix( PSolver % Matrix ) 
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
     REAL(KIND=dp) :: x(4), y(4), z(4)

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
     CALL GetRefPElementNodes( Element,x,y,z )
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

     ! Assign rest of indexes if neccessary
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


END MODULE DefUtils

!> \}  // end of subgroup
!> \}  // end of group

