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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Jun 1997
! *
! *****************************************************************************/

!> \ingroup ElmerLib 
!> \{

!------------------------------------------------------------------------------
!>  List handling utilities. In Elmer all the keywords are saved to a list,
!> and later accessed from it repeatedly. Therefore these subroutines are 
!> essential in Elmer programming.
!------------------------------------------------------------------------------
MODULE Lists

   USE Messages
   USE GeneralUtils
#ifdef USE_ISO_C_BINDINGS
   USE LoadMod
#endif

   IMPLICIT NONE

   INTEGER, PARAMETER :: LIST_TYPE_CONSTANT_SCALAR = 1
   INTEGER, PARAMETER :: LIST_TYPE_CONSTANT_TENSOR = 2
   INTEGER, PARAMETER :: LIST_TYPE_VARIABLE_SCALAR = 3
   INTEGER, PARAMETER :: LIST_TYPE_VARIABLE_TENSOR = 4
   INTEGER, PARAMETER :: LIST_TYPE_LOGICAL = 5
   INTEGER, PARAMETER :: LIST_TYPE_STRING  = 6
   INTEGER, PARAMETER :: LIST_TYPE_INTEGER = 7
   INTEGER, PARAMETER :: LIST_TYPE_CONSTANT_SCALAR_STR = 8
   INTEGER, PARAMETER :: LIST_TYPE_CONSTANT_TENSOR_STR = 9
   INTEGER, PARAMETER :: LIST_TYPE_VARIABLE_SCALAR_STR = 10
   INTEGER, PARAMETER :: LIST_TYPE_VARIABLE_TENSOR_STR = 11
   INTEGER, PARAMETER :: LIST_TYPE_CONSTANT_SCALAR_PROC = 12

#ifndef USE_ISO_C_BINDINGS
   INTERFACE
     FUNCTION ExecIntFunction( Proc,Md ) RESULT(int)
       USE Types
#ifdef SGI
       INTEGER :: Proc
#else
       INTEGER(KIND=AddrInt) :: Proc
#endif
       TYPE(Model_t) :: Md

       INTEGER :: int
     END FUNCTION ExecIntFunction
   END INTERFACE

   INTERFACE
     FUNCTION ExecRealFunction( Proc,Md,Node,Temp ) RESULT(dbl)
       USE Types

#ifdef SGI
       INTEGER :: Proc
#else
       INTEGER(KIND=AddrInt) :: Proc
#endif
       TYPE(Model_t) :: Md
       INTEGER :: Node
       REAL(KIND=dp) :: Temp(*)

       REAL(KIND=dp) :: dbl
     END FUNCTION ExecRealFunction
   END INTERFACE

   INTERFACE
     SUBROUTINE ExecRealArrayFunction( Proc,Md,Node,Temp,F )
       USE Types

#ifdef SGI
       INTEGER :: Proc
#else
       INTEGER(KIND=AddrInt) :: Proc
#endif
       TYPE(Model_t) :: Md
       INTEGER :: Node,n1,n2
       REAL(KIND=dp) :: Temp(*)

       REAL(KIND=dp) :: F(:,:)
     END SUBROUTINE ExecRealArrayFunction
   END INTERFACE

   INTERFACE
     SUBROUTINE ExecRealVectorFunction( Proc,Md,Node,T,F )
       USE Types
       INTEGER(KIND=AddrInt) :: Proc
       TYPE(Model_t) :: Md
       INTEGER :: Node,n1,n2
       REAL(KIND=dp) :: T(*), F(:)
     END SUBROUTINE ExecRealVectorFunction
   END INTERFACE

   INTERFACE
     FUNCTION ExecConstRealFunction( Proc,Md,x,y,z ) RESULT(dbl)
       USE Types

#ifdef SGI
       INTEGER :: Proc
#else
       INTEGER(KIND=AddrInt) :: Proc
#endif
       TYPE(Model_t) :: Md

       REAL(KIND=dp) :: dbl,x,y,z
     END FUNCTION ExecConstRealFunction
   END INTERFACE
#endif

   TYPE String_stack_t
      TYPE(Varying_string) :: Name
      TYPE(String_stack_t), POINTER :: Next => Null()
   END TYPE String_stack_t

   CHARACTER(:), ALLOCATABLE, SAVE, PRIVATE :: Namespace
   !$OMP THREADPRIVATE(NameSpace)

   TYPE(String_stack_t), SAVE, PRIVATE, POINTER :: Namespace_stack => Null()
   !$OMP THREADPRIVATE(NameSpace_stack)

   CHARACTER(:), ALLOCATABLE, SAVE, PRIVATE :: ActiveListName
   !$OMP THREADPRIVATE(ActiveListName)

   TYPE(String_stack_t), SAVE, PRIVATE, POINTER :: Activename_stack => Null()
   !$OMP THREADPRIVATE(Activename_stack)

   TYPE(ValueList_t), POINTER, SAVE, PRIVATE  :: TimerList => NULL()
   LOGICAL, SAVE, PRIVATE :: TimerPassive, TimerResults

CONTAINS

!------------------------------------------------------------------------------
  FUNCTION InitialPermutation( Perm,Model,Solver,Mesh, &
                   Equation,DGSolver,GlobalBubbles ) RESULT(k)
!------------------------------------------------------------------------------
     USE PElementMaps
     TYPE(Model_t)  :: Model
     TYPE(Mesh_t)   :: Mesh
     TYPE(Solver_t), TARGET :: Solver
     INTEGER :: Perm(:)
     INTEGER :: k
     CHARACTER(LEN=*) :: Equation
     LOGICAL, OPTIONAL :: DGSolver, GlobalBubbles
!------------------------------------------------------------------------------
     INTEGER i,j,l,t,n,e, EDOFs, FDOFs, BDOFs, ndofs, el_id
     INTEGER :: Indexes(128)
     INTEGER, POINTER :: Def_Dofs(:)
     INTEGER, ALLOCATABLE :: EdgeDOFs(:), FaceDOFs(:)
     LOGICAL :: FoundDG, DG, GB, Found, Radiation
     TYPE(Element_t),POINTER :: Element, Edge, Face
!------------------------------------------------------------------------------
     Perm = 0
     k = 0
     EDOFs = Mesh % MaxEdgeDOFs
     FDOFs = Mesh % MaxFaceDOFs
     BDOFs = Mesh % MaxBDOFs

     GB = .FALSE.
     IF ( PRESENT(GlobalBubbles) ) GB=GlobalBubbles

     DG = .FALSE.
     IF ( PRESENT(DGSolver) ) DG=DGSolver
     FoundDG = .FALSE.
     IF ( DG ) THEN
       DO t=1,Mesh % NumberOfEdges
         n = 0
         Element => Mesh % Edges(t) % BoundaryInfo % Left
         IF ( ASSOCIATED( Element ) ) THEN
             IF ( CheckElementEquation(Model,Element,Equation) ) THEN
                FoundDG = FoundDG .OR. Element % DGDOFs > 0
                DO j=1,Element % DGDOFs
                  n = n + 1
                  Indexes(n) = Element % DGIndexes(j)
                END DO
             END IF
         END IF

         Element => Mesh % Edges(t) % BoundaryInfo % Right
         IF ( ASSOCIATED( Element ) ) THEN
             IF ( CheckElementEquation(Model,Element,Equation) ) THEN
                FoundDG = FoundDG .OR. Element % DGDOFs > 0
                DO j=1,Element % DGDOFs
                  n = n + 1
                  Indexes(n) = Element % DGIndexes(j)
                END DO
             END IF
         END IF

         DO i=1,n
            j = Indexes(i)
            IF ( Perm(j) == 0 ) THEN
               k = k + 1
              Perm(j) = k
            END IF
         END DO
       END DO

       DO t=1,Mesh % NumberOfFaces
         n = 0
         Element => Mesh % Faces(t) % BoundaryInfo % Left
         IF ( ASSOCIATED( Element ) ) THEN
             IF ( CheckElementEquation(Model,Element,Equation) ) THEN
                FoundDG = FoundDG .OR. Element % DGDOFs > 0
                DO j=1,Element % DGDOFs
                   n = n + 1
                   Indexes(n) = Element % DGIndexes(j)
                END DO
             END IF
         END IF

         Element => Mesh % Faces(t) % BoundaryInfo % Right
         IF ( ASSOCIATED( Element ) ) THEN
             IF ( CheckElementEquation(Model,Element,Equation) ) THEN
                FoundDG = FoundDG .OR. Element % DGDOFs > 0
                DO j=1,Element % DGDOFs
                   n = n + 1
                   Indexes(n) = Element % DGIndexes(j)
                END DO
             END IF
         END IF

         DO i=1,n
            j = Indexes(i)
            IF ( Perm(j) == 0 ) THEN
                k = k + 1
               Perm(j) = k
            END IF
         END DO
       END DO

       IF ( FoundDG ) THEN
          RETURN ! Discontinuous galerkin !!!
       END IF
     END IF


     IF ( ANY(Solver % Def_Dofs(:,:,6)>=0) ) THEN
       IF ( Mesh % NumberOFEdges>0 ) THEN
          ALLOCATE(EdgeDOFs(Mesh % NumberOfEdges))
          EdgeDOFs=0;
       END IF

       IF ( Mesh % NumberOFFaces>0 ) THEN
         ALLOCATE(FaceDOFs(Mesh % NumberOfFaces))
         FaceDOFs=0;
       END IF

       n = Mesh % NumberOfBulkElements + Mesh % NumberOFBoundaryElements
       t = 1
       DO WHILE( t <= n )
         DO WHILE( t<=n )
           Element => Mesh % Elements(t)
           IF ( CheckElementEquation( Model, Element, Equation ) ) EXIT
           t = t + 1
         END DO
         IF ( t>n ) EXIT

         el_id = Element % TYPE % ElementCode / 100

         Def_Dofs => Solver % Def_Dofs(el_id,Element % BodyId,:)
         IF ( ASSOCIATED(Element % EdgeIndexes) ) THEN
           DO i=1,Element % TYPE % NumberOfEdges
             j = Element % EdgeIndexes(i)
             EdgeDOFs(j)=MAX(EdgeDOFs(j),getEdgeDOFs(Element,Def_Dofs(6)))
           END DO
         END IF

         IF ( ASSOCIATED(Element % FaceIndexes) ) THEN
           DO i=1,Element % TYPE % NumberOfFaces
             j = Element % FaceIndexes(i)
             FaceDOFs(j)=MAX(FaceDOFs(j),getFaceDOFs(Element,Def_Dofs(6),i))
           END DO
         END IF
         t=t+1
       END DO
     END IF


     n = Mesh % NumberOfBulkElements + Mesh % NumberOFBoundaryElements
     t = 1
     DO WHILE( t <= n )
       DO WHILE( t<=n )
         Element => Mesh % Elements(t)
         IF ( CheckElementEquation( Model, Element, Equation ) ) EXIT
         t = t + 1
       END DO

       IF ( t > n ) EXIT

       el_id = Element % TYPE % ElementCode / 100
       Def_Dofs => Solver % Def_Dofs(el_id,Element % BodyId,:)
       ndofs = Element % NDOFs
       IF ( Def_Dofs(1) >= 0 ) ndofs=Def_Dofs(1)*Element % TYPE % NumberOfNodes
       DO i=1,ndofs
         j = Element % NodeIndexes(i)
         IF ( Perm(j) == 0 ) THEN
           k = k + 1
           Perm(j) = k
         END IF
       END DO

       IF ( ASSOCIATED( Element % EdgeIndexes ) ) THEN
          DO i=1,Element % TYPE % NumberOfEdges
             Edge => Mesh % Edges( Element % EdgeIndexes(i) )
             ndofs = 0
             IF ( Def_Dofs(2) >= 0) THEN
               ndofs = Def_Dofs(2)
             ELSE IF (Def_Dofs(6)>=0) THEN
               ndofs = EdgeDOFs(Element % EdgeIndexes(i))
!              IF (Def_Dofs(6)==0) ndofs = MAX(Edge % BDOFs,ndofs)
               ndofs = MAX(Edge % BDOFs,ndofs)
             END IF

             DO e=1,ndofs
                j = Mesh % NumberOfNodes + EDOFs*(Element % EdgeIndexes(i)-1) + e
                IF ( Perm(j) == 0 ) THEN
                   k = k + 1
                   Perm(j) =  k
                END IF
             END DO
          END DO
       END IF

       IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
          DO i=1,Element % TYPE % NumberOfFaces
             Face => Mesh % Faces( Element % FaceIndexes(i) )
             l = MAX(0,Def_Dofs(3))
             j = Face % TYPE % ElementCode/100
             IF(l==0) THEN
               IF (ASSOCIATED(Face % BoundaryInfo % Left)) THEN
                 e = Face % BoundaryInfo % Left % BodyId
                 l = MAX(0,Solver % Def_Dofs(j+6,e,5))
               END IF
               IF (ASSOCIATED(Face % BoundaryInfo % Right)) THEN
                 e = Face % BoundaryInfo % Right % BodyId
                 l = MAX(l,Solver % Def_Dofs(j+6,e,5))
               END IF
             END IF
             ndofs = 0
             IF ( l >= 0) THEN
               ndofs = l
             ELSE IF (Def_Dofs(6)>=0) THEN
               ndofs = FaceDOFs(Element % FaceIndexes(i))
!              IF ( Def_Dofs(6)==0 ) ndofs = MAX(Face % BDOFs,ndofs)
               ndofs = MAX(Face % BDOFs,ndofs)
             END IF

             DO e=1,ndofs
                j = Mesh % NumberOfNodes + EDOFs*Mesh % NumberOfEdges + &
                          FDOFs*(Element % FaceIndexes(i)-1) + e
                IF ( Perm(j) == 0 ) THEN
                   k = k + 1
                   Perm(j) =  k
                END IF
             END DO
          END DO
       END IF

       IF ( GB .AND. ASSOCIATED(Element % BubbleIndexes) ) THEN
         ndofs = 0
         IF ( Def_Dofs(5) >= 0) THEN
            ndofs = Def_Dofs(5)
         ELSE IF (Def_Dofs(6)>=0) THEN
            ndofs = GetBubbleDOFs(Element, Def_Dofs(6))
            IF ( Def_Dofs(6)==0 ) ndofs = MAX(Element % BDOFs,ndofs)
         END IF

         DO i=1,ndofs
            j = Mesh % NumberOfNodes + EDOFs*Mesh % NumberOfEdges + &
                 FDOFs*Mesh % NumberOfFaces + Element % BubbleIndexes(i)
            IF ( Perm(j) == 0 ) THEN
               k = k + 1
               Perm(j) =  k
            END IF
         END DO
       END IF

       t = t + 1
     END DO

     Radiation = ListGetLogical( Solver % Values, 'Radiation Solver', Found )
     IF ( Radiation .OR. Equation == 'heat equation' ) THEN
        t = Mesh % NumberOfBulkElements + 1
        n = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
        DO WHILE( t<= n )
          Element => Mesh % Elements(t)
          IF ( ASSOCIATED( Element % BoundaryInfo % GebhardtFactors) ) THEN
             DO i=1,Element % TYPE % NumberOfNodes
               j = Element % NodeIndexes(i)
               IF ( Perm(j) == 0 ) THEN
                 k = k + 1
                 Perm(j) = k
               END IF
             END DO
          END IF
          t = t + 1
        END DO
     END IF

     t = Mesh % NumberOfBulkElements + 1
     n = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
     DO WHILE( t<= n )
       Element => Mesh % Elements(t)
       IF ( Element % TYPE % ElementCode == 102 ) THEN
          DO i=1,Element % TYPE % NumberOfNodes
            j = Element % NodeIndexes(i)
            IF ( Perm(j) == 0 ) THEN
              k = k + 1
              Perm(j) = k
            END IF
          END DO
       END IF
       t = t + 1
     END DO

     IF ( ALLOCATED(EdgeDOFs) ) DEALLOCATE(EdgeDOFs)
     IF ( ALLOCATED(FaceDOFs) ) DEALLOCATE(FaceDOFs)
!------------------------------------------------------------------------------
   END FUNCTION InitialPermutation
!------------------------------------------------------------------------------


!---------------------------------------------------------------------------
!>   Check if given element belongs to a body for which given equation
!>   should be solved.
!---------------------------------------------------------------------------
    FUNCTION CheckElementEquation( Model,Element,Equation ) RESULT(Flag)
      TYPE(Element_t), POINTER :: Element
      TYPE(Model_t) :: Model
      CHARACTER(LEN=*) :: Equation

      LOGICAL :: Flag,Found

      INTEGER :: k,body_id
       
      Flag = .FALSE.
      body_id = Element % BodyId
      IF ( body_id > 0 .AND. body_id <= Model % NumberOfBodies ) THEN
         k = ListGetInteger( Model % Bodies(body_id) % Values, 'Equation', &
                 minv=1, maxv=Model % NumberOFEquations )
         IF ( k > 0 ) THEN
           Flag = ListGetLogical(Model % Equations(k) % Values,Equation,Found)
         END IF
      END IF
!---------------------------------------------------------------------------
   END FUNCTION CheckElementEquation
!---------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Changes the string to all lower case to allow string comparison.
!------------------------------------------------------------------------------
    FUNCTION StringToLowerCase( to,from,same_len ) RESULT(n)
!------------------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(in)  :: from
      CHARACTER(LEN=*), INTENT(out) :: to
      LOGICAL, OPTIONAL, INTENT(in) :: same_len
!------------------------------------------------------------------------------
      INTEGER :: n
      INTEGER :: i,j,nlen
      INTEGER, PARAMETER :: A=ICHAR('A'),Z=ICHAR('Z'),U2L=ICHAR('a')-ICHAR('A')

      n = LEN(to)
      IF (.NOT.PRESENT(same_len)) THEN
        DO i=LEN(from),1,-1
          IF ( from(i:i) /= ' ' ) EXIT
        END DO
        IF ( n>i ) THEN
          to(i+1:n) = ' '
          n=i
        END IF
      END IF

      nlen = n
      DO i=1,nlen
        j = ICHAR( from(i:i) )
        IF ( j >= A .AND. j <= Z ) THEN
          to(i:i) = CHAR(j+U2L)
        ELSE
          to(i:i) = from(i:i)
          IF ( to(i:i)=='[') n=i-1
        END IF
      END DO
    END FUNCTION StringToLowerCase
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Adds a new variable to the list of variables. 
!> The structures need to be allocated externally beforehand. 
!------------------------------------------------------------------------------
    SUBROUTINE VariableAdd( Variables,Mesh,Solver,Name,DOFs,Values,&
      Perm,Output,Secondary, TYPE )
!------------------------------------------------------------------------------
      TYPE(Variable_t), POINTER :: Variables
      TYPE(Mesh_t),   TARGET :: Mesh
      TYPE(Solver_t), TARGET :: Solver
      CHARACTER(LEN=*) :: Name
      INTEGER :: DOFs
      INTEGER, OPTIONAL :: TYPE
      REAL(KIND=dp), POINTER :: Values(:)
      LOGICAL, OPTIONAL :: Output
      INTEGER, OPTIONAL, POINTER :: Perm(:)
      LOGICAL, OPTIONAL :: Secondary
!------------------------------------------------------------------------------
      LOGICAL :: stat
      TYPE(Variable_t), POINTER :: ptr,ptr1,ptr2
!------------------------------------------------------------------------------

      IF ( .NOT.ASSOCIATED(Variables) ) THEN
        ALLOCATE(Variables)
        ptr => Variables
      ELSE
        ALLOCATE( ptr )
      END IF

      ptr % NameLen = StringToLowerCase( ptr % Name,Name )

      IF ( .NOT. ASSOCIATED(ptr, Variables) ) THEN
        ptr1 => Variables
        ptr2 => Variables
        DO WHILE( ASSOCIATED( ptr1 ) )
           IF ( ptr % Name == ptr1 % Name ) THEN
              DEALLOCATE( ptr )
              RETURN
           END IF
           ptr2 => ptr1
           ptr1 => ptr1 % Next
         END DO
         ptr2 % Next => ptr
      END IF
      ptr % Next => NULL()

      ptr % DOFs = DOFs
      IF ( PRESENT( Perm ) ) THEN
        ptr % Perm => Perm
      ELSE
        ptr % Perm => NULL()
      END IF
      ptr % Norm = 0.0d0
      ptr % PrevNorm = 0.0d0
      ptr % Values => Values
      NULLIFY( ptr % PrevValues )
      NULLIFY( ptr % EigenValues, ptr % EigenVectors )

      ptr % NonlinChange = 0.0_dp
      ptr % SteadyChange = 0.0_dp
      ptr % NonlinValues => NULL(); ptr % SteadyValues => NULL()
      ptr % NonlinIter = 0

      ptr % Solver => Solver
      ptr % PrimaryMesh => Mesh

      ptr % Valid  = .TRUE.
      ptr % Output = .TRUE.
      ptr % Secondary = .FALSE.
      ptr % ValuesChanged = .TRUE.

! Converged information undefined = -1, not = 0, yes = 1
      ptr % NonlinConverged = -1
      ptr % SteadyConverged = -1    

      IF ( PRESENT( Secondary ) ) THEN
        IF(Secondary) PRINT *,'Secondary:',TRIM(name)
        ptr % Secondary = Secondary
      END IF

      IF ( PRESENT( TYPE ) ) ptr % TYPE = TYPE
      IF ( PRESENT( Output ) ) ptr % Output = Output
!------------------------------------------------------------------------------
    END SUBROUTINE VariableAdd
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE ReleaseVariableList( VariableList )
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: VariableList
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: Ptr(:)
    LOGICAL :: GotValues
    INTEGER :: i, n, m
    TYPE(Variable_t), POINTER :: Var, Var1
!------------------------------------------------------------------------------

    Var => VariableList
    DO WHILE( ASSOCIATED( Var ) ) 

!      This is used to skip variables such as time, timestep, timestep size etc.
       IF( SIZE( Var % Values ) == Var % DOFs ) THEN
         Var => Var % Next
         CYCLE
       END IF 

       SELECT CASE( Var % Name )
       CASE( 'coordinate 1', 'coordinate 2', 'coordinate 3' )
         Var => Var % Next
         CYCLE
       END SELECT

	IF( Var % Secondary ) THEN
          Var => Var % Next
          CYCLE
        END IF

       IF (Var % DOFs > 1) THEN
         Var => Var % Next
         CYCLE
       END IF
!
!      Check that the variable is actually allocated,
!      not pointer to some other variables memory:
!      ----------------------------------------------

       GotValues = .TRUE.
       Var1 => VariableList
       DO WHILE( ASSOCIATED( Var1 ) )
          IF (.NOT.ASSOCIATED(Var,Var1)) THEN
             IF ( ASSOCIATED(Var1 % Values) ) THEN
                DO i=1,Var1 % DOFs
                   ptr => Var1 % Values(i::Var1 % DOFs)
                   IF ( ASSOCIATED(Var % Values,ptr) ) THEN
                      GotValues = .FALSE.
                      EXIT
                   END IF
                END DO
             END IF
          END IF
          IF (.NOT. GotValues) EXIT
          Var1 => Var1 % Next
       END DO

       IF (ASSOCIATED(Var % Perm)) THEN
         Var1 => VariableList
         DO WHILE(ASSOCIATED(Var1))
           IF (.NOT.ASSOCIATED(Var,Var1)) THEN
             IF (ASSOCIATED(Var % Perm,Var1 % Perm)) &
               Var1 % Perm => NULL()
           END IF
           Var1 => Var1 % Next
         END DO

         DEALLOCATE( Var % Perm)
       END IF

       IF ( GotValues ) THEN

        IF ( ASSOCIATED( Var % Values ) ) &
            DEALLOCATE( Var % Values )

         IF ( ASSOCIATED( Var % PrevValues ) ) &
	   DEALLOCATE( Var % PrevValues )

         IF ( ASSOCIATED( Var % EigenValues ) ) &
            DEALLOCATE( Var % EigenValues )

         IF ( ASSOCIATED( Var % EigenVectors ) ) &
            DEALLOCATE( Var % EigenVectors )

         IF ( ASSOCIATED( Var % SteadyValues ) ) &
            DEALLOCATE( Var % SteadyValues )

         IF ( ASSOCIATED( Var % NonlinValues ) ) &
            DEALLOCATE( Var % NonlinValues )
       END IF
       NULLIFY( Var % EigenVectors, Var % EigenValues )
       NULLIFY( Var % Values, Var % PrevValues, Var % Perm )
       NULLIFY( Var % SteadyValues, Var % NonlinValues )

       Var => Var % Next
    END DO

    Var => VariableList
    DO WHILE( ASSOCIATED( Var ) )
       IF ( Var % Secondary ) THEN
         Var => Var % Next
         CYCLE
       END IF

       IF ( Var % DOFs > 1 ) THEN
         IF ( ASSOCIATED( Var % Values ) ) &
            DEALLOCATE( Var % Values )

         IF ( ASSOCIATED( Var % Perm ) ) &
            DEALLOCATE( Var % Perm )

         IF ( ASSOCIATED( Var % PrevValues ) ) &
            DEALLOCATE( Var % PrevValues )

         IF ( ASSOCIATED( Var % EigenValues ) ) &
            DEALLOCATE( Var % EigenValues )

         IF ( ASSOCIATED( Var % EigenVectors ) ) &
            DEALLOCATE( Var % EigenVectors )

         IF ( ASSOCIATED( Var % NonlinValues ) ) &
            DEALLOCATE( Var % NonlinValues )
       END IF
       NULLIFY( Var % EigenVectors, Var % EigenValues )
       NULLIFY( Var % Values, Var % PrevValues, Var % Perm )
       NULLIFY( Var % SteadyValues, Var % NonlinValues )

       Var => Var % Next
    END DO

!   Deallocate mesh variable list:
!   ------------------------------
    Var => VariableList
    DO WHILE( ASSOCIATED( Var ) )
       Var1 => Var % Next
       DEALLOCATE( Var )
       Var => Var1
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ReleaseVariableList
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Deletes a variable (by name) from list of variables
!------------------------------------------------------------------------------
  SUBROUTINE VariableRemove(Variables, NameIn)
    
    IMPLICIT NONE
!-----------------------------------------------
    TYPE(Variable_t), POINTER :: Variables
    CHARACTER(LEN=*) :: NameIn
!-----------------------------------------------    
    TYPE(Variable_t), POINTER :: Var, Prev, RmVar
    CHARACTER(LEN=LEN_TRIM(NameIn)) :: Name
    LOGICAL :: GotIt
    INTEGER :: k

    GotIt = .FALSE.

    Var => Variables
    Prev => NULL()
    k = StringToLowerCase(Name, NameIn,.TRUE.)

    WRITE(Message,'(a,a)') "Removing variable: ",Name(1:k)
    CALL Info("VariableRemove",Message, Level=10)

    !Find variable by name, and hook up % Next appropriately
    DO WHILE(ASSOCIATED(Var))
       IF( Var % NameLen == k ) THEN
          IF(Var % Name(1:k) == Name(1:k)) THEN
             GotIt = .TRUE.
             RmVar => Var
             IF(ASSOCIATED(Prev)) THEN
                !Link up variables either side of removed var
                Prev % Next => Var % Next
             ELSE
                !If this was the first variable, we point Variables
                !at the next one...
                Variables => Var % Next
             END IF
             EXIT
          END IF
       END IF
       Prev => Var
       Var => Prev % Next
    END DO

    IF(.NOT. GotIt) THEN
       CALL Warn("VariableRemove","Couldn't find the variable, returning...")
       RETURN
    END IF

    RmVar % Next => NULL()

    !ReleaseVariableList was intended to deallocate an entire list of variables,
    !but by nullifying RmVar % Next, we have effectively isolated RmVar in 
    !its own variable list.
    CALL ReleaseVariableList( RmVar )
!------------------------------------------------------------------------------
  END SUBROUTINE VariableRemove
!------------------------------------------------------------------------------

 

!------------------------------------------------------------------------------
!> For vectors the individual components are added also to the list 
!> of variables. This routine makes the addition of vectors less laborious.
!> Also allocates the field values if not given in the parameter list. 
!------------------------------------------------------------------------------
    SUBROUTINE VariableAddVector( Variables,Mesh,Solver,Name,DOFs,Values,&
      Perm,Output,Secondary,Global,InitValue )
!------------------------------------------------------------------------------
      TYPE(Variable_t), POINTER :: Variables
      TYPE(Mesh_t),   TARGET :: Mesh
      TYPE(Solver_t), TARGET :: Solver
      CHARACTER(LEN=*) :: Name
      INTEGER, OPTIONAL :: DOFs
      REAL(KIND=dp), OPTIONAL, POINTER :: Values(:)
      LOGICAL, OPTIONAL :: Output
      INTEGER, OPTIONAL, POINTER :: Perm(:)
      LOGICAL, OPTIONAL :: Secondary
      LOGICAL, OPTIONAL :: Global
      REAL(KIND=dp), OPTIONAL :: InitValue
!------------------------------------------------------------------------------
      CHARACTER(LEN=MAX_NAME_LEN) :: tmpname
      REAL(KIND=dp), POINTER :: Component(:), TmpValues(:)
      INTEGER :: i,nsize, ndofs
!------------------------------------------------------------------------------

      IF( PRESENT( DOFs ) ) THEN
        ndofs = Dofs
      ELSE
        ndofs = 1
      END IF

      IF(PRESENT(Values)) THEN
        TmpValues => Values
      ELSE
        IF( PRESENT( Perm ) ) THEN 
          nsize = MAXVAL( Perm ) 
        ELSE IF( PRESENT( Global ) ) THEN
          IF( Global ) THEN
            nsize = 1 
          ELSE
            nsize = Mesh % NumberOfNodes
          END IF
        ELSE
          nsize = Mesh % NumberOfNodes          
        END IF
        NULLIFY(TmpValues)
        ALLOCATE(TmpValues(ndofs*nsize))
        TmpValues = 0.0_dp         
      END IF

      IF( PRESENT( InitValue ) ) THEN
        TmpValues = InitValue
      END IF

      IF( nDOFs > 1 ) THEN
        DO i=1,nDOFs
          tmpname = ComponentName(Name,i)
          Component => TmpValues(i::nDOFs)
          CALL VariableAdd( Variables,Mesh,Solver,TmpName,1,Component,&
              Perm,Output,Secondary)
        END DO
      END IF

      CALL VariableAdd( Variables,Mesh,Solver,Name,nDOFs,TmpValues,&
            Perm,Output,Secondary )

!------------------------------------------------------------------------------
    END SUBROUTINE VariableAddVector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    FUNCTION MeshProjector( Mesh1, Mesh2, &
         UseQuadrantTree, Trans ) RESULT( ProjectorMatrix )
!------------------------------------------------------------------------------
       TYPE(Mesh_t) :: Mesh1, Mesh2
       LOGICAL, OPTIONAL :: UseQuadrantTree,Trans
       TYPE(Matrix_t), POINTER :: ProjectorMatrix
!------------------------------------------------------------------------------
       TYPE(Projector_t), POINTER :: Projector
!------------------------------------------------------------------------------
       INTERFACE
         SUBROUTINE InterpolateMeshToMeshQ( OldMesh, NewMesh, OldVariables, &
             NewVariables, UseQuadrantTree, Projector, MaskName, FoundNodes, NewMaskPerm)
           USE Types
           TYPE(Variable_t), POINTER, OPTIONAL :: OldVariables, NewVariables
           TYPE(Mesh_t), TARGET  :: OldMesh, NewMesh
           LOGICAL, OPTIONAL :: UseQuadrantTree,FoundNodes(:)
           CHARACTER(LEN=*),OPTIONAL :: MaskName
           TYPE(Projector_t), POINTER, OPTIONAL :: Projector
           INTEGER, OPTIONAL, POINTER :: NewMaskPerm(:)  !< Mask the new variable set by the given MaskName when trying to define the interpolation.
         END SUBROUTINE InterpolateMeshToMeshQ
       END INTERFACE

       IF ( PRESENT(UseQuadrantTree) ) THEN
          CALL InterpolateMeshToMeshQ( Mesh1, Mesh2, &
                   UseQuadrantTree=UseQuadrantTree, Projector=Projector )
       ELSE
          CALL InterpolateMeshToMeshQ( Mesh1, Mesh2, Projector=Projector )
       END IF
 
       ProjectorMatrix => Projector % Matrix
       IF ( PRESENT(Trans) ) THEN
          IF ( Trans ) THEN
             ProjectorMatrix => Projector % TMatrix
          END IF
       END IF
!------------------------------------------------------------------------------
    END FUNCTION MeshProjector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Find a variable by its name from the list of variables. 
!> If it not to be found in the current mesh, interpolation between
!> meshes is automatically requested for.
!------------------------------------------------------------------------------
    RECURSIVE FUNCTION VariableGet( Variables, Name, ThisOnly, MaskName, UnfoundFatal ) RESULT(Var)
!------------------------------------------------------------------------------
      TYPE(Variable_t), POINTER :: Variables
      CHARACTER(LEN=*) :: Name
      LOGICAL, OPTIONAL :: ThisOnly
      CHARACTER(LEN=*),OPTIONAL :: MaskName
      LOGICAL, OPTIONAL :: UnfoundFatal
!------------------------------------------------------------------------------
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Projector_t), POINTER :: Projector
      TYPE(Variable_t), POINTER :: Var,PVar,Tmp,AidVar
      REAL(KIND=dp), POINTER :: Vals(:)
      INTEGER :: i,k,n, DOFs
      LOGICAL :: Found, GlobalBubbles, UseProjector
      CHARACTER(LEN=LEN_TRIM(Name)) :: str
      CHARACTER(LEN=MAX_NAME_LEN) :: tmpname
#ifdef USE_ISO_C_BINDINGS
      DOUBLE PRECISION :: t1
#else
      DOUBLE PRECISION :: t1,CPUTime
#endif
!------------------------------------------------------------------------------
      INTERFACE
        SUBROUTINE InterpolateMeshToMesh( OldMesh, NewMesh, OldVariables, &
            NewVariables, UseQuadrantTree, Projector, MaskName, UnfoundNodes )
          USE Types
          TYPE(Variable_t), POINTER, OPTIONAL :: OldVariables, NewVariables
          TYPE(Mesh_t), TARGET  :: OldMesh, NewMesh
          LOGICAL, OPTIONAL :: UseQuadrantTree
          LOGICAL, POINTER, OPTIONAL :: UnfoundNodes(:)
          CHARACTER(LEN=*),OPTIONAL :: MaskName
          TYPE(Projector_t), POINTER, OPTIONAL :: Projector
        END SUBROUTINE InterpolateMeshToMesh
      END INTERFACE
!------------------------------------------------------------------------------

      k = StringToLowerCase( str,Name,.TRUE. )

      Tmp => Variables
      DO WHILE( ASSOCIATED(tmp) )
        IF ( tmp % NameLen == k ) THEN
          IF ( tmp % Name(1:k) == str(1:k) ) THEN

            IF ( Tmp % Valid ) THEN
               Var => Tmp
               RETURN
            END IF
            EXIT

          END IF
        END IF
        tmp => tmp % Next
      END DO
      Var => Tmp


!------------------------------------------------------------------------------
      IF ( PRESENT(ThisOnly) ) THEN
         IF ( ThisOnly ) THEN
            IF ( PRESENT(UnfoundFatal) ) THEN
               IF ( UnfoundFatal ) THEN
                  WRITE(Message,'(A,A)') "Failed to find variable ",Name
                  CALL Fatal("VariableGet",Message)
               END IF
            END IF
            RETURN
         END IF
      END IF

!------------------------------------------------------------------------------
      NULLIFY( PVar )
      Mesh => CurrentModel % Meshes
      DO WHILE( ASSOCIATED( Mesh ) )

        IF ( .NOT.ASSOCIATED( Variables, Mesh % Variables ) ) THEN
          PVar => VariableGet( Mesh % Variables, Name, ThisOnly=.TRUE. )
          IF ( ASSOCIATED( PVar ) ) THEN
            IF ( ASSOCIATED( Mesh, PVar % PrimaryMesh ) ) THEN
              EXIT
            END IF
          END IF
        END IF
        Mesh => Mesh % Next
      END DO

      IF ( .NOT.ASSOCIATED( PVar ) ) THEN
         IF ( PRESENT(UnfoundFatal) ) THEN
            IF ( UnfoundFatal ) THEN
               WRITE(Message,'(A,A)') "Failed to find or interpolate variable ",Name
               CALL Fatal("VariableGet",Message)
            END IF
         END IF
         RETURN
      END IF

!------------------------------------------------------------------------------

      IF ( .NOT.ASSOCIATED( Tmp ) ) THEN
         GlobalBubbles = ListGetLogical(Pvar % Solver % Values, &
               'Bubbles in Global System', Found)
         IF (.NOT.Found) GlobalBubbles=.TRUE.

         DOFs = CurrentModel % Mesh % NumberOfNodes * PVar % DOFs
         IF ( GlobalBubbles ) DOFs = DOFs + CurrentModel % Mesh % MaxBDOFs * &
              CurrentModel % Mesh % NumberOfBulkElements * PVar % DOFs

         ALLOCATE( Var )
         ALLOCATE( Var % Values(DOFs) )
         Var % Values = 0

         NULLIFY( Var % Perm )
         IF ( ASSOCIATED( PVar % Perm ) ) THEN
            ALLOCATE( Var % Perm( DOFs/Pvar % DOFs ) )

            n = InitialPermutation( Var % Perm, CurrentModel, PVar % Solver, &
                CurrentModel % Mesh, ListGetString(PVar % Solver % Values,'Equation'), &
                 GlobalBubbles=GlobalBubbles )

            IF ( n==0 ) n=CurrentModel % Mesh % NumberOfNodes

            IF ( n == CurrentModel % Mesh % NumberOfNodes ) THEN
               DO i=1,n 
                  Var % Perm(i) = i
               END DO
            END IF
         END IF

         CALL VariableAdd( Variables, PVar % PrimaryMesh, PVar % Solver, &
           PVar % Name, PVar % DOFs, Var % Values, Var % Perm, PVar % Output ) 

         Var => VariableGet( Variables, Name, ThisOnly=.TRUE. )

         NULLIFY( Var % PrevValues )
         IF ( ASSOCIATED( PVar % PrevValues ) ) THEN
            ALLOCATE( Var % PrevValues( DOFs, SIZE(PVar % PrevValues,2) ) )
         END IF

         IF ( PVar % Name(1:PVar % NameLen) == 'flow solution' ) THEN
           Vals => Var % Values( 1: SIZE(Var % Values) : PVar % DOFs )
           CALL VariableAdd( Variables, PVar % PrimaryMesh, PVar % Solver, &
                  'Velocity 1', 1,  Vals, Var % Perm, PVar % Output ) 

           Tmp => VariableGet( Variables, 'Velocity 1', .TRUE. )
           NULLIFY( Tmp % PrevValues )
           IF ( ASSOCIATED( Var % PrevValues ) )  &
              Tmp % PrevValues => Var % PrevValues(1::PVar % DOFs,:)

           Vals => Var % Values( 2: SIZE(Var % Values) : PVar % DOFs )
           CALL VariableAdd( Variables, PVar % PrimaryMesh, PVar % Solver, &
                  'Velocity 2', 1,  Vals, Var % Perm, PVar % Output ) 

           Tmp => VariableGet( Variables, 'Velocity 2', .TRUE. )
           NULLIFY( Tmp % PrevValues )
           IF ( ASSOCIATED( Var % PrevValues ) ) &
              Tmp % PrevValues => Var % PrevValues(2::PVar % DOFs,:)

           IF ( PVar % DOFs == 3 ) THEN
             Vals => Var % Values( 3 : SIZE(Var % Values) : PVar % DOFs )
             CALL VariableAdd( Variables, PVar % PrimaryMesh, PVar % Solver, &
                    'Pressure', 1,  Vals, Var % Perm, PVar % Output ) 
           ELSE
             Vals => Var % Values( 3: SIZE(Var % Values) : PVar % DOFs )
             CALL VariableAdd( Variables, PVar % PrimaryMesh, PVar % Solver, &
                  'Velocity 3', 1,  Vals, Var % Perm, PVar % Output ) 

             Tmp => VariableGet( Variables, 'Velocity 3', .TRUE. )
             NULLIFY( Tmp % PrevValues )
             IF ( ASSOCIATED( Var % PrevValues ) ) &
                 Tmp % PrevValues => Var % PrevValues(3::PVar % DOFs,:)

             Vals => Var % Values( 4: SIZE(Var % Values) : PVar % DOFs )
             CALL VariableAdd( Variables, PVar % PrimaryMesh, PVar % Solver, &
                    'Pressure', 1,  Vals, Var % Perm, PVar % Output ) 
           END IF

           Tmp => VariableGet( Variables, 'Pressure', .TRUE. )
           NULLIFY( Tmp % PrevValues )
           IF ( ASSOCIATED( Var % PrevValues ) ) &
              Tmp % PrevValues => Var % PrevValues(PVar % DOFs::PVar % DOFs,:)
         ELSE
           IF ( PVar % DOFs > 1 ) THEN
             DO i=1,PVar % DOFs
               Vals => Var % Values( i: SIZE(Var % Values) : PVar % DOFs )
               tmpname = ComponentName( PVar % Name, i )
               CALL VariableAdd( Variables, PVar % PrimaryMesh, PVar % Solver, &
                       tmpname, 1, Vals, Var % Perm, PVar % Output ) 

               Tmp => VariableGet( Variables, tmpname, .TRUE. )
               NULLIFY( Tmp % PrevValues )
               IF ( ASSOCIATED( Var % PrevValues ) ) &
                  Tmp % PrevValues => Var % PrevValues(i::PVar % DOFs,:)
             END DO
           END IF
        END IF
 
        Var => VariableGet( Variables, Name, ThisOnly=.TRUE. )
      END IF

!------------------------------------------------------------------------------
! Build a temporary variable list of variables to be interpolated
!------------------------------------------------------------------------------
      ALLOCATE( Tmp )
      Tmp = PVar
      Var => Tmp
      NULLIFY( Var % Next )

      IF ( PVar % Name(1:PVar % NameLen) == 'flow solution' ) THEN
        ALLOCATE( Var % Next )
        Var => Var % Next
        Var = VariableGet( PVar % PrimaryMesh % Variables, 'Velocity 1' )

        ALLOCATE( Var % Next )
        Var => Var % Next
        Var  = VariableGet(  PVar % PrimaryMesh % Variables, 'Velocity 2' )

        IF ( PVar % DOFs == 4 ) THEN
          ALLOCATE( Var % Next )
          Var => Var % Next
          Var  = VariableGet( PVar % PrimaryMesh % Variables, 'Velocity 3' )
        END IF

        ALLOCATE( Var % Next )
        Var => Var % Next
        Var = VariableGet( PVar % PrimaryMesh % Variables, 'Pressure' )
        NULLIFY( Var % Next )
        Var => Tmp
      ELSE IF ( PVar % DOFs > 1 ) THEN
        DO i=1,PVar % DOFs
          ALLOCATE( Var % Next )
          tmpname = ComponentName( PVar % Name, i )
          Var % Next = VariableGet( PVar % PrimaryMesh % Variables, tmpname )
          Var => Var % Next
        END DO
        NULLIFY( Var % Next )
        Var => Tmp
      END IF

!------------------------------------------------------------------------------
! interpolation call
!------------------------------------------------------------------------------
      t1 = CPUTime()

      UseProjector = ListGetLogical(CurrentModel % Simulation,'Use Mesh Projector',Found)
      IF( .NOT. Found ) UseProjector = .TRUE.

      IF( PRESENT( MaskName ) ) THEN
       CALL InterpolateMeshToMesh( PVar % PrimaryMesh, &
            CurrentModel % Mesh, Var, Variables, MaskName=MaskName )
      ELSE IF( UseProjector ) THEN
        CALL InterpolateMeshToMesh( PVar % PrimaryMesh, &
            CurrentModel % Mesh, Var, Variables, Projector=Projector )
      ELSE
        AidVar => VariableGet( CurrentModel % Mesh % Variables, Name, ThisOnly = .TRUE. ) 
        IF( ASSOCIATED( AidVar ) ) THEN
          AidVar % Values = 0.0_dp
        END IF
        CALL InterpolateMeshToMesh( PVar % PrimaryMesh, &
            CurrentModel % Mesh, Var, Variables )        
      END IF

      WRITE( Message,'(A,ES12.3)' ) 'Interpolation time for > '//TRIM(Name)//' < :', CPUTime()-t1
      CALL Info( 'VariableGet', Message, Level=7 )

!------------------------------------------------------------------------------
! free the temporary list
!------------------------------------------------------------------------------
      DO WHILE( ASSOCIATED( Tmp ) )
         Var => Tmp % Next
         DEALLOCATE( Tmp )
         Tmp => Var
      END DO
!------------------------------------------------------------------------------
      Var => VariableGet( Variables, Name, ThisOnly=.TRUE. )
      Var % Valid = .TRUE.
      Var % ValuesChanged = .TRUE.

      IF ( Var % Name(1:Var % NameLen) == 'flow solution' ) THEN
        Tmp => VariableGet( Variables, 'Velocity 1', ThisOnly=.TRUE. )
        IF ( ASSOCIATED(Tmp) ) THEN
          Tmp % Valid = .TRUE.
          Tmp % ValuesChanged = .TRUE.
        END IF

        Tmp => VariableGet( Variables, 'Velocity 2', ThisOnly=.TRUE. )
        IF ( ASSOCIATED(Tmp) ) THEN
          Tmp % Valid = .TRUE.
          Tmp % ValuesChanged = .TRUE.
        END IF

        IF ( Var % DOFs == 4 ) THEN
          Tmp  => VariableGet( Variables, 'Velocity 3', ThisOnly=.TRUE. )
          IF ( ASSOCIATED(Tmp) ) THEN
            Tmp % Valid = .TRUE.
            Tmp % ValuesChanged = .TRUE.
          END IF
        END IF

        Tmp => VariableGet( Variables, 'Pressure', ThisOnly=.TRUE. )
        IF ( ASSOCIATED(Tmp) ) THEN
          Tmp % Valid = .TRUE.
          Tmp % ValuesChanged = .TRUE.
        END IF
      ELSE IF ( Var % DOFs > 1 ) THEN
        DO i = 1,Var % DOFs
           tmpname = ComponentName( Var % Name, i )
           Tmp => VariableGet( Variables, tmpname, ThisOnly=.TRUE. )
           IF ( ASSOCIATED(Tmp) ) THEN
             Tmp % Valid = .TRUE.
             Tmp % ValuesChanged = .TRUE.
           END IF
        END DO
      END IF
!------------------------------------------------------------------------------
    END FUNCTION VariableGet 
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 FUNCTION ListHead(list) RESULT(head)
!------------------------------------------------------------------------------
   TYPE(ValueList_t) :: List
   TYPE(ValueListEntry_t), POINTER :: Head
!------------------------------------------------------------------------------
   head => List % Head
!------------------------------------------------------------------------------
 END FUNCTION ListHead
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 FUNCTION ListEmpty(list) RESULT(l)
!------------------------------------------------------------------------------
    LOGICAL :: L
    TYPE(ValueList_t) :: list
!------------------------------------------------------------------------------
    L = .NOT.ASSOCIATED(list % head)
!------------------------------------------------------------------------------
 END FUNCTION ListEmpty
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Allocates a new value list.
!------------------------------------------------------------------------------
  FUNCTION ListAllocate() RESULT(ptr)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: ptr
     ALLOCATE( ptr )
     ptr % Head => Null()
!------------------------------------------------------------------------------
  END FUNCTION ListAllocate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Allocates a new value list.
!------------------------------------------------------------------------------
  FUNCTION ListEntryAllocate() RESULT(ptr)
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr

     ALLOCATE( ptr )
     ptr % PROCEDURE = 0
     ptr % TYPE = 0
     ptr % Name = ' '
     ptr % NameLen = 0
     ptr % CValue = ' '
     ptr % LValue = .FALSE.
     NULLIFY( ptr % CubicCoeff )
     NULLIFY( ptr % Cumulative )
     NULLIFY( ptr % Next )
     NULLIFY( ptr % FValues )
     NULLIFY( ptr % TValues )
     NULLIFY( ptr % IValues )
!------------------------------------------------------------------------------
  END FUNCTION ListEntryAllocate
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Deletes a value list.
!------------------------------------------------------------------------------
  SUBROUTINE ListDelete( ptr )
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr

     IF ( ASSOCIATED(ptr % CubicCoeff) ) DEALLOCATE(ptr % CubicCoeff)
     IF ( ASSOCIATED(ptr % Cumulative) ) DEALLOCATE(ptr % Cumulative)
     IF ( ASSOCIATED(ptr % FValues) ) DEALLOCATE(ptr % FValues)
     IF ( ASSOCIATED(ptr % TValues) ) DEALLOCATE(ptr % TValues)
     IF ( ASSOCIATED(ptr % IValues) ) DEALLOCATE(ptr % IValues)
     DEALLOCATE( ptr )
!------------------------------------------------------------------------------
  END SUBROUTINE ListDelete
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Removes an entry from the list by its name.
!------------------------------------------------------------------------------
  SUBROUTINE ListRemove( List, Name )
!------------------------------------------------------------------------------
     TYPE(ValueList_t) :: List
     CHARACTER(LEN=*)  :: Name
!------------------------------------------------------------------------------
     CHARACTER(LEN=LEN_TRIM(Name)) :: str
     INTEGER :: k
     LOGICAL :: Found
     TYPE(ValueListEntry_t), POINTER :: ptr, prev
!------------------------------------------------------------------------------
     IF ( ASSOCIATED(List % Head) ) THEN
       k = StringToLowerCase( str,Name,.TRUE. )
       ptr  => List % Head
       Prev => ptr
       DO WHILE( ASSOCIATED(ptr) )
         IF ( ptr % NameLen == k .AND. ptr % Name(1:k) == str(1:k) ) THEN
            IF ( ASSOCIATED(ptr,List % Head) ) THEN
               List % Head => ptr % Next
               Prev => List % Head
            ELSE
               Prev % Next => ptr % Next
            END IF
            CALL ListDelete(ptr)
            EXIT
         ELSE
           Prev => ptr
           ptr  => ptr % Next 
         END IF
       END DO
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE ListRemove
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Adds an entry to the list by its name and returns a handle to the new entry. If the entry is 
!> already existing return the existing one. 
!------------------------------------------------------------------------------
  FUNCTION ListAdd( List, Name ) RESULT(NEW)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     TYPE(ValueListEntry_t), POINTER :: new
!------------------------------------------------------------------------------
     CHARACTER(LEN=LEN_TRIM(Name)) :: str
     INTEGER :: k
     LOGICAL :: Found
     TYPE(ValueListEntry_t), POINTER :: ptr, prev
!------------------------------------------------------------------------------
     Prev => NULL()
     Found = .FALSE.

     IF(.NOT.ASSOCIATED(List)) List => ListAllocate()
     New => ListEntryAllocate()

     IF ( ASSOCIATED(List % Head) ) THEN
       k = StringToLowerCase( str,Name,.TRUE. )
       ptr  => List % Head
       NULLIFY( prev )
       DO WHILE( ASSOCIATED(ptr) )
         IF ( ptr % NameLen == k .AND. ptr % Name(1:k) == str(1:k) ) THEN
           Found = .TRUE.
           EXIT
         ELSE
           Prev => ptr
           ptr  => ptr % Next 
         END IF
       END DO

       IF ( Found ) THEN
         New % Next => ptr % Next
         IF ( ASSOCIATED( prev ) ) THEN
           Prev % Next => New
         ELSE
           List % Head => New
         END IF
         CALL ListDelete( Ptr )
       ELSE
         IF ( ASSOCIATED(prev) ) THEN
           prev % next => NEW
         ELSE
           NEW % Next => List % Head % Next
           List % Head % Next => NEW
         END IF
       END IF
     ELSE
       List % Head => NEW
     END IF
!------------------------------------------------------------------------------
   END FUNCTION ListAdd
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Sets a namespace string that is used in all list get commands 
!> to check for an entry with the namespace, and then continuing to check the one without.
!------------------------------------------------------------------------------
   SUBROUTINE ListSetNamespace(str)
!------------------------------------------------------------------------------
     CHARACTER(LEN=*) :: str
!------------------------------------------------------------------------------
     CHARACTER(LEN=LEN_TRIM(str)) :: str_lcase
!------------------------------------------------------------------------------
     INTEGER :: n
!------------------------------------------------------------------------------
     n = StringToLowerCase( str_lcase,str,.TRUE. )
     NameSpace = str_lcase
!------------------------------------------------------------------------------
   END SUBROUTINE ListSetNamespace
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Returns the active namespace.
!------------------------------------------------------------------------------
   FUNCTION ListGetNamespace(str) RESULT(l)
!------------------------------------------------------------------------------
    LOGICAL :: l 
    CHARACTER(:), ALLOCATABLE :: str
!------------------------------------------------------------------------------
    l = .FALSE.
    IF ( Namespace /= '' ) THEN
      l = .TRUE.
      str = Namespace
    END IF
!------------------------------------------------------------------------------
   END FUNCTION ListGetNamespace
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE ListPushNamespace(str)
!------------------------------------------------------------------------------
     CHARACTER(LEN=*) :: str
!------------------------------------------------------------------------------
     LOGICAL :: L
     CHARACTER(:), ALLOCATABLE :: tstr
     TYPE(String_stack_t), POINTER :: stack
!------------------------------------------------------------------------------
     ALLOCATE(stack)
     L = ListGetNameSpace(tstr)
     IF(ALLOCATED(tstr)) THEN
       stack % name = tstr
     ELSE
       stack % name = ''
     END IF
     stack % next => Namespace_stack
     Namespace_stack => stack
     CALL ListSetNamespace(str)
!------------------------------------------------------------------------------
   END SUBROUTINE ListPushNamespace
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE ListPopNamespace()
!------------------------------------------------------------------------------
     TYPE(String_stack_t), POINTER :: stack
!------------------------------------------------------------------------------
     IF(ASSOCIATED(Namespace_stack)) THEN
       Namespace = Namespace_stack % name
       stack => Namespace_stack
       Namespace_stack => stack % Next
       DEALLOCATE(stack)
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE ListPopNamespace
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE ListPushActivename(str)
!------------------------------------------------------------------------------
     CHARACTER(LEN=*) :: str
!------------------------------------------------------------------------------
     LOGICAL :: L
     TYPE(String_stack_t), POINTER :: stack
!------------------------------------------------------------------------------
     ALLOCATE(stack)
     stack % name = ListGetActiveName()
     stack % next => Activename_stack
     Activename_stack => stack
     ActiveListName = str
!------------------------------------------------------------------------------
   END SUBROUTINE ListPushActiveName
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE ListPopActiveName()
!------------------------------------------------------------------------------
     TYPE(String_stack_t), POINTER :: stack
!------------------------------------------------------------------------------
     IF(ASSOCIATED(Activename_stack)) THEN
       ActiveListName = Activename_stack % name
       stack => Activename_stack
       Activename_stack => stack % Next
       DEALLOCATE(stack)
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE ListPopActiveName
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   FUNCTION ListGetActiveName() RESULT(str)
!------------------------------------------------------------------------------
    CHARACTER(:), ALLOCATABLE :: str
!------------------------------------------------------------------------------
    str = ActiveListName
!------------------------------------------------------------------------------
   END FUNCTION ListGetActiveName
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Finds an entry in the list by its name and returns a handle to it.
!------------------------------------------------------------------------------
   FUNCTION ListFind( list, name, Found) RESULT(ptr)
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: name
     LOGICAL, OPTIONAL :: Found
!------------------------------------------------------------------------------
     CHARACTER(:), ALLOCATABLE :: strn
     CHARACTER(LEN=LEN_TRIM(Name)) :: str
!------------------------------------------------------------------------------
     INTEGER :: k, k1, n

     IF(PRESENT(Found)) Found = .FALSE.
     ptr => NULL()
     IF(.NOT.ASSOCIATED(List)) RETURN

     k = StringToLowerCase( str,Name,.TRUE. )

     IF ( ListGetNamespace(strn) ) THEN
       strn = strn //' '//str(1:k)
       k1 = LEN(strn)
       ptr => List % Head
       DO WHILE( ASSOCIATED(ptr) )
          n = ptr % NameLen
          IF ( n==k1 ) THEN
            IF ( ptr % Name(1:n) == strn ) EXIT
          END IF
          ptr => ptr % Next
       END DO
     END IF

     IF ( .NOT. ASSOCIATED(ptr) ) THEN
       Ptr => List % Head
       DO WHILE( ASSOCIATED(ptr) )
         n = ptr % NameLen
         IF ( n==k ) THEN
           IF ( ptr % Name(1:n) == str(1:n) ) EXIT
         END IF
         ptr => ptr % Next
       END DO
     END IF

     IF ( PRESENT(Found) ) THEN
       Found = ASSOCIATED(ptr)
     ELSE IF (.NOT.ASSOCIATED(ptr) ) THEN
       CALL Warn( 'ListFind', ' ' )
       WRITE(Message,*) 'Requested property: ', '[',TRIM(Name),'], not found'
       CALL Warn( 'ListFind', Message )
       CALL Warn( 'ListFind', ' ' )
     END IF
!------------------------------------------------------------------------------
   END FUNCTION ListFind
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Finds an entry in the list by its name and returns a handle to it.
!> This one just finds a keyword with the same start as specified by 'name'.
!------------------------------------------------------------------------------
   FUNCTION ListFindPrefix( list, name, Found) RESULT(ptr)
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
     TYPE(ValueList_t), POINTER :: list
     CHARACTER(LEN=*) :: name
     LOGICAL, OPTIONAL :: Found
!------------------------------------------------------------------------------
     CHARACTER(:), ALLOCATABLE :: strn
     CHARACTER(LEN=LEN_TRIM(Name)) :: str
!------------------------------------------------------------------------------
     INTEGER :: k, k1, n, m

     ptr => NULL()
     IF(.NOT.ASSOCIATED(List)) RETURN

     k = StringToLowerCase( str,Name,.TRUE. )
     IF ( ListGetNamespace(strn) ) THEN
       strn = strn //' '//str(1:k)
       k1 = LEN(strn)
       ptr => List % Head
       DO WHILE( ASSOCIATED(ptr) )
          n = ptr % NameLen
          IF ( n >= k1 ) THEN
            IF ( ptr % Name(1:k1) == strn ) EXIT
          END IF
          ptr => ptr % Next
       END DO
     END IF

     IF ( .NOT. ASSOCIATED(ptr) ) THEN
       Ptr => List % Head
       DO WHILE( ASSOCIATED(ptr) )
         n = ptr % NameLen
         IF ( n >= k ) THEN
           IF ( ptr % Name(1:k) == str(1:k) ) EXIT
         END IF
         ptr => ptr % Next
       END DO
     END IF

     IF ( PRESENT(Found) ) THEN
       Found = ASSOCIATED(ptr)
     ELSE IF (.NOT.ASSOCIATED(ptr) ) THEN
       CALL Warn( 'ListFindPrefix', ' ' )
       WRITE(Message,*) 'Requested prefix: ', '[',TRIM(Name),'], not found'
       CALL Warn( 'ListFindPrefix', Message )
       CALL Warn( 'ListFindPrefix', ' ' )
     END IF
!------------------------------------------------------------------------------
   END FUNCTION ListFindPrefix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Finds an entry in the list by its name and returns a handle to it.
!> This one just finds a keyword with the same end as specified by 'name'.
!------------------------------------------------------------------------------
   FUNCTION ListFindSuffix( list, name, Found) RESULT(ptr)
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
     TYPE(ValueList_t), POINTER :: list
     CHARACTER(LEN=*) :: name
     LOGICAL, OPTIONAL :: Found
!------------------------------------------------------------------------------
     CHARACTER(LEN=LEN_TRIM(Name)) :: str
!------------------------------------------------------------------------------
     INTEGER :: k, k1, n, m

     ptr => Null()
     IF(.NOT.ASSOCIATED(List)) RETURN
     
     k = StringToLowerCase( str,Name,.TRUE. )
     Ptr => List % Head
     DO WHILE( ASSOCIATED(ptr) )
       n = ptr % NameLen
       IF ( n >= k ) THEN
         IF ( ptr % Name(n-k+1:n) == str(1:k) ) EXIT
       END IF
       ptr => ptr % Next
     END DO

     IF ( PRESENT(Found) ) THEN
       Found = ASSOCIATED(ptr)
     ELSE IF (.NOT.ASSOCIATED(ptr) ) THEN
       CALL Warn( 'ListFindSuffix', ' ' )
       WRITE(Message,*) 'Requested suffix: ', '[',TRIM(Name),'], not found'
       CALL Warn( 'ListFindSuffix', Message )
       CALL Warn( 'ListFindSuffix', ' ' )
     END IF
!------------------------------------------------------------------------------
   END FUNCTION ListFindSuffix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Check if the suffix exists in the list.
!------------------------------------------------------------------------------
   FUNCTION ListCheckSuffix( List, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
     TYPE(ValuelistEntry_t), POINTER :: ptr
     
     ptr => ListFindSuffix( List, Name, Found )
!------------------------------------------------------------------------------
   END FUNCTION ListCheckSuffix
!------------------------------------------------------------------------------
  


!------------------------------------------------------------------------------
!> Check if the keyword is with the given suffix is present in any boundary condition.
!------------------------------------------------------------------------------
   FUNCTION ListCheckSuffixAnyBC( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
     INTEGER :: bc
     TYPE(ValuelistEntry_t), POINTER :: ptr
     
     Found = .FALSE.
     DO bc = 1,Model % NumberOfBCs
       ptr => ListFindSuffix( Model % BCs(bc) % Values, Name, Found )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListCheckSuffixAnyBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Check if the keyword is with the given suffix is present in any body.
!------------------------------------------------------------------------------
   FUNCTION ListCheckSuffixAnyBody( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
     INTEGER :: body_id
     TYPE(ValuelistEntry_t), POINTER :: ptr
     
     Found = .FALSE.
     DO body_id = 1,Model % NumberOfBodies
       ptr => ListFindSuffix( Model % Bodies(body_id) % Values, Name, Found )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListCheckSuffixAnyBody
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Check if the keyword is with the given suffix is present in any material.
!------------------------------------------------------------------------------
   FUNCTION ListCheckSuffixAnyMaterial( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
     INTEGER :: mat_id
     TYPE(ValuelistEntry_t), POINTER :: ptr
     
     Found = .FALSE.
     DO mat_id = 1,Model % NumberOfMaterials
       ptr => ListFindSuffix( Model % Materials(mat_id) % Values, Name, Found )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListCheckSuffixAnyMaterial
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Check if the keyword is with the given suffix is present in any body force.
!------------------------------------------------------------------------------
   FUNCTION ListCheckSuffixAnyBodyForce( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
     INTEGER :: bf_id
     TYPE(ValuelistEntry_t), POINTER :: ptr
     
     Found = .FALSE.
     DO bf_id = 1,Model % NumberOfBodyForces
       ptr => ListFindSuffix( Model % BodyForces(bf_id) % Values, Name, Found )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListCheckSuffixAnyBodyForce
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Finds an entry related to vector keyword of type "name" or "name i", i=1,2,3.
!> This could save time since it will detect at one sweep whether the keyword
!> for a vector is given, and whether it is componentwise or not. 
!> There is a caveat since currently the "i" is not checked and possibly 
!> the user could mix the formats and the chosen one would be random.  
!------------------------------------------------------------------------------
   FUNCTION ListFindVectorPrefix( list, name, ComponentWise,Found ) RESULT(ptr)
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
     TYPE(ValueList_t), POINTER :: list
     CHARACTER(LEN=*) :: name
     LOGICAL :: ComponentWise
     LOGICAL, OPTIONAL :: Found
!------------------------------------------------------------------------------
     CHARACTER(:), ALLOCATABLE :: strn
     CHARACTER(LEN=LEN_TRIM(Name)) :: str
!------------------------------------------------------------------------------
     INTEGER :: k, k1, n, m

     ptr => NULL()
     IF(.NOT.ASSOCIATED(List)) RETURN

     k = StringToLowerCase( str,Name,.TRUE. )

     IF ( ListGetNamespace(strn) ) THEN
       strn = strn //' '//str(1:k)
       k1 = LEN(strn)
       ptr => List % Head
       DO WHILE( ASSOCIATED(ptr) )
          n = ptr % NameLen
          IF ( n == k1 ) THEN
            IF ( ptr % Name(1:k1) == strn ) THEN
              ComponentWise = .FALSE.
              EXIT
            END IF
          ELSE IF( n == k1 + 2 ) THEN
            IF ( ptr % Name(1:k1+1) == strn//' ' ) THEN
              ComponentWise = .TRUE.
              EXIT
            END IF
          END IF
          ptr => ptr % Next
       END DO
     END IF

     IF ( .NOT. ASSOCIATED(ptr) ) THEN
       Ptr => List % Head
       DO WHILE( ASSOCIATED(ptr) )
         n = ptr % NameLen
         IF ( n == k ) THEN
           IF ( ptr % Name(1:k) == str(1:k) ) THEN
             ComponentWise = .FALSE.
             EXIT
           END IF
         ELSE IF( n == k + 2 ) THEN
           IF ( ptr % Name(1:k+1) == str(1:k)//' ' ) THEN
             ComponentWise = .TRUE.
             EXIT
           END IF
         END IF
         ptr => ptr % Next
       END DO
     END IF

     IF ( PRESENT(Found) ) THEN
       Found = ASSOCIATED(ptr)
     ELSE IF (.NOT.ASSOCIATED(ptr) ) THEN
       CALL Warn( 'ListFindVectorPrefix', ' ' )
       WRITE(Message,*) 'Requested vector prefix: ', '[',TRIM(Name),'], not found'
       CALL Warn( 'ListFindVectorPrefix', Message )
       CALL Warn( 'ListFindVectorPrefix', ' ' )
     END IF
!------------------------------------------------------------------------------
   END FUNCTION ListFindVectorPrefix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Finds a keyword with the given basename and normalizes it with a 
!> constant coefficients for all future request of the keyword.
!------------------------------------------------------------------------------
   SUBROUTINE ListSetCoefficients( list, name, coeff )
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: list
     CHARACTER(LEN=*) :: name
     REAL(KIND=dp) :: coeff
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr, ptr2
     CHARACTER(LEN=LEN_TRIM(Name)) :: str
     INTEGER :: k, k1, n, n2, m

     IF(.NOT.ASSOCIATED(List)) RETURN

     k = StringToLowerCase( str,Name,.TRUE. )
     
     Ptr => list % Head
     DO WHILE( ASSOCIATED(ptr) )
       n = ptr % NameLen
       IF ( n >= k ) THEN
         ! Did we find a keyword which has the correct suffix?
         IF ( ptr % Name(n-k+1:n) == str(1:k) ) THEN
           Ptr2 => list % Head
           DO WHILE( ASSOCIATED(ptr2) )
             n2 = ptr2 % NameLen
             IF( n2 + k <= n ) THEN

               ! Did we find the corresponding keyword without the suffix?
               IF ( ptr2 % Name(1:n2) == ptr % Name(1:n2) ) THEN
                 WRITE( Message,'(A,ES12.5)') 'Normalizing > '//&
                     TRIM( ptr2 % Name )// ' < by ',Coeff
                 CALL Info('ListSetCoefficients',Message)
                 ptr2 % Coeff = Coeff
                 EXIT
               END IF

             END IF
             ptr2 => ptr2 % Next
           END DO
         END IF
       END IF
       ptr => ptr % Next
     END DO

   END SUBROUTINE ListSetCoefficients
 


!> Copies an entry from 'ptr' to an entry in *different* list with the same content.
!-----------------------------------------------------------------------------------
   SUBROUTINE ListCopyItem( ptr, list )

     TYPE(ValueListEntry_t), POINTER :: ptr
     TYPE(ValueList_t), POINTER :: list
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptrb, ptrnext

     ptrb => ListAdd( List, ptr % Name ) 

     ptrnext => ptrb % next
     ptrb = ptr
     ptrb % next => ptrnext

   END SUBROUTINE ListCopyItem


!> Checks two lists for a given keyword. If it is given then 
!> copy it as it is to the 2nd list.
!------------------------------------------------------------------------------
   SUBROUTINE ListCompareAndCopy( list, listb, name, Found )
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: list, listb
     CHARACTER(LEN=*) :: name
     LOGICAL :: Found
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
     CHARACTER(LEN=LEN_TRIM(Name)) :: str
     INTEGER :: k, n

     k = StringToLowerCase( str,Name,.TRUE. )
     Found = .FALSE.

     ! Find the keyword from the 1st list 
     Ptr => List % Head
     DO WHILE( ASSOCIATED(ptr) )
       n = ptr % NameLen
       IF ( n==k ) THEN
         IF ( ptr % Name(1:n) == str(1:n) ) EXIT
       END IF
       ptr => ptr % Next
     END DO
     
     IF(.NOT. ASSOCIATED( ptr ) ) RETURN
     
     ! Add the same entry to the 2nd list 
     CALL ListCopyItem( ptr, listb ) 
     Found = .TRUE.

   END SUBROUTINE ListCompareAndCopy
 
  
!------------------------------------------------------------------------------
!> Just checks if a entry is present in the list.
!------------------------------------------------------------------------------
   FUNCTION ListCheckPresent( List,Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------
     ptr => ListFind(List,Name,Found)
!------------------------------------------------------------------------------
   END FUNCTION ListCheckPresent
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Just checks if a prefix is present in the list.
!------------------------------------------------------------------------------
   FUNCTION ListCheckPrefix( List,Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------
     ptr => ListFindPrefix(List,Name,Found)
!------------------------------------------------------------------------------
   END FUNCTION ListCheckPrefix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Check if the keyword is with the given prefix is present in any boundary condition.
!------------------------------------------------------------------------------
   FUNCTION ListCheckPrefixAnyBC( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
     INTEGER :: bc
     TYPE(ValuelistEntry_t), POINTER :: ptr
     
     Found = .FALSE.
     DO bc = 1,Model % NumberOfBCs
       ptr => ListFindPrefix( Model % BCs(bc) % Values, Name, Found )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListCheckPrefixAnyBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Check if the keyword is with the given prefix is present in any body.
!------------------------------------------------------------------------------
   FUNCTION ListCheckPrefixAnyBody( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
     INTEGER :: body_id
     TYPE(ValuelistEntry_t), POINTER :: ptr
     
     Found = .FALSE.
     DO body_id = 1,Model % NumberOfBodies
       ptr => ListFindPrefix( Model % Bodies(body_id) % Values, Name, Found )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListCheckPrefixAnyBody
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Check if the keyword is with the given prefix is present in any material.
!------------------------------------------------------------------------------
   FUNCTION ListCheckPrefixAnyMaterial( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
     INTEGER :: mat_id
     TYPE(ValuelistEntry_t), POINTER :: ptr
     
     Found = .FALSE.
     DO mat_id = 1,Model % NumberOfMaterials
       ptr => ListFindPrefix( Model % Materials(mat_id) % Values, Name, Found )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListCheckPrefixAnyMaterial
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Check if the keyword is with the given prefix is present in any body force.
!------------------------------------------------------------------------------
   FUNCTION ListCheckPrefixAnyBodyForce( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
     INTEGER :: bf_id
     TYPE(ValuelistEntry_t), POINTER :: ptr
     
     Found = .FALSE.
     DO bf_id = 1,Model % NumberOfBodyForces
       ptr => ListFindPrefix( Model % BodyForces(bf_id) % Values, Name, Found )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListCheckPrefixAnyBodyForce
!------------------------------------------------------------------------------





!------------------------------------------------------------------------------
!> Adds a string to the list.
!------------------------------------------------------------------------------
    SUBROUTINE ListAddString( List,Name,CValue,CaseConversion )
!------------------------------------------------------------------------------
      TYPE(ValueList_t), POINTER :: List
      CHARACTER(LEN=*) :: Name
      CHARACTER(LEN=*) :: CValue
      LOGICAL, OPTIONAL :: CaseConversion
!------------------------------------------------------------------------------
      INTEGER :: k
      LOGICAL :: DoCase
      TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------
      ptr => ListAdd( List, Name )

      DoCase = .TRUE.
      IF ( PRESENT(CaseConversion) ) DoCase = CaseConversion

      IF ( DoCase ) THEN
        k = StringToLowerCase( ptr % CValue,CValue )
      ELSE
        k = MIN( MAX_NAME_LEN,LEN(CValue) )
        ptr % CValue(1:k) = CValue(1:k)
      END IF

      ptr % TYPE   = LIST_TYPE_STRING
      ptr % NameLen = StringToLowerCase( Ptr % Name,Name )
!------------------------------------------------------------------------------
    END SUBROUTINE ListAddString
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Adds a logical entry to the list.
!------------------------------------------------------------------------------
    SUBROUTINE ListAddLogical( List,Name,LValue )
!------------------------------------------------------------------------------
      TYPE(ValueList_t), POINTER :: List
      CHARACTER(LEN=*) :: Name
      LOGICAL :: LValue
!------------------------------------------------------------------------------
      TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------
      ptr => ListAdd( List, Name )
      Ptr % LValue = LValue
      Ptr % TYPE   = LIST_TYPE_LOGICAL

      Ptr % NameLen = StringToLowerCase( ptr % Name,Name )
    END SUBROUTINE ListAddLogical
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Adds an integer to the list.
!------------------------------------------------------------------------------
    SUBROUTINE ListAddInteger( List,Name,IValue,Proc )
!------------------------------------------------------------------------------
      TYPE(ValueList_t), POINTER :: List
      CHARACTER(LEN=*) :: Name
      INTEGER :: IValue
      INTEGER(Kind=AddrInt), OPTIONAL :: Proc
!------------------------------------------------------------------------------
      TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------
      ptr => ListAdd( List, Name )
      IF ( PRESENT(Proc) ) ptr % PROCEDURE = Proc

      ALLOCATE( ptr % IValues(1) )
      ptr % IValues(1) = IValue
      ptr % TYPE       = LIST_TYPE_INTEGER

      ptr % NameLen = StringToLowerCase( ptr % Name,Name )
    END SUBROUTINE ListAddInteger
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Adds an integer array to the list.
!------------------------------------------------------------------------------
    SUBROUTINE ListAddIntegerArray( List,Name,N,IValues,Proc )
!------------------------------------------------------------------------------
      TYPE(ValueList_t), POINTER :: List
      CHARACTER(LEN=*) :: Name
      INTEGER :: N
      INTEGER :: IValues(N)
      INTEGER(KIND=AddrInt), OPTIONAL :: Proc
!------------------------------------------------------------------------------
      TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------
      ptr => ListAdd( List, Name )

      ALLOCATE( ptr % IValues(N) )

      IF ( PRESENT(Proc) ) ptr % PROCEDURE = Proc

      ptr % TYPE  = LIST_TYPE_CONSTANT_TENSOR
      ptr % IValues(1:n) = IValues(1:n)

      ptr % NameLen = StringToLowerCase( ptr % Name,Name )
    END SUBROUTINE ListAddIntegerArray
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Adds a constant real value to the list.
!------------------------------------------------------------------------------
    SUBROUTINE ListAddConstReal( List,Name,FValue,Proc,CValue )
!------------------------------------------------------------------------------
      TYPE(ValueList_t), POINTER :: List
      CHARACTER(LEN=*) :: Name
      CHARACTER(LEN=*), OPTIONAL :: Cvalue
      REAL(KIND=dp) :: FValue
      INTEGER(KIND=AddrInt), OPTIONAL :: Proc
!------------------------------------------------------------------------------
      TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------
      ptr => ListAdd( List, Name )

      NULLIFY( ptr % TValues )
      ALLOCATE( ptr % FValues(1,1,1) )

      ptr % FValues = FValue
      ptr % TYPE  = LIST_TYPE_CONSTANT_SCALAR

      IF ( PRESENT(Proc) ) THEN
        ptr % PROCEDURE = Proc
        IF( Proc /= 0 ) THEN
          ptr % TYPE = LIST_TYPE_CONSTANT_SCALAR_PROC
        END IF
      END IF

      IF ( PRESENT( CValue ) ) THEN
         ptr % Cvalue = CValue
         ptr % TYPE  = LIST_TYPE_CONSTANT_SCALAR_STR
      END IF

      ptr % NameLen = StringToLowerCase( ptr % Name,Name )
    END SUBROUTINE ListAddConstReal
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Adds a linear dependency defined by a table of values, [x,y] to the list.
!------------------------------------------------------------------------------
    SUBROUTINE ListAddDepReal(List,Name,DependName,N,TValues, &
               FValues,Proc,CValue,CubicTable, Monotone)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name,DependName
     CHARACTER(LEN=*), OPTIONAL :: Cvalue
     INTEGER :: N
     LOGICAL, OPTIONAL :: CubicTable, Monotone
     REAL(KIND=dp) :: FValues(N)
     REAL(KIND=dp) :: TValues(N)
     INTEGER(KIND=AddrInt), OPTIONAL :: Proc
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------
     ptr => ListAdd( List, Name )
     IF ( PRESENT(Proc) ) ptr % PROCEDURE = Proc

     ALLOCATE( ptr % FValues(1,1,n),ptr % TValues(n) )

     ptr % TValues = TValues(1:n)
     ptr % FValues(1,1,:) = FValues(1:n)
     ptr % TYPE = LIST_TYPE_VARIABLE_SCALAR

     IF ( n>3 .AND. PRESENT(CubicTable)) THEN
       IF ( CubicTable ) THEN
         ALLOCATE(ptr % CubicCoeff(n))
         CALL CubicSpline(n,ptr % TValues,Ptr % Fvalues(1,1,:), &
                    Ptr % CubicCoeff, Monotone )
       END IF
     END IF

     ALLOCATE(ptr % Cumulative(n))
     CALL CumulativeIntegral(ptr % TValues, Ptr % FValues(1,1,:), &
          Ptr % CubicCoeff, Ptr % Cumulative )

     ptr % NameLen = StringToLowerCase( ptr % Name,Name )
     ptr % DepNameLen = StringToLowerCase( ptr % DependName,DependName )

     IF ( PRESENT( Cvalue ) ) THEN
        ptr % CValue = CValue
        ptr % TYPE = LIST_TYPE_VARIABLE_SCALAR_STR
     END IF

   END SUBROUTINE ListAddDepReal
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Adds a constant real valued array to the list.
!------------------------------------------------------------------------------
    SUBROUTINE ListAddConstRealArray( List,Name,N,M,FValues,Proc,CValue )
!------------------------------------------------------------------------------
      TYPE(ValueList_t), POINTER :: List
      CHARACTER(LEN=*) :: Name
      CHARACTER(LEN=*), OPTIONAL :: Cvalue
      INTEGER :: N,M
      REAL(KIND=dp) :: FValues(:,:)
      INTEGER(KIND=AddrInt), OPTIONAL :: Proc
!------------------------------------------------------------------------------
      TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------
      ptr => ListAdd( List, Name )

      NULLIFY( ptr % TValues )
      ALLOCATE( ptr % FValues(N,M,1) )


      ptr % TYPE  = LIST_TYPE_CONSTANT_TENSOR
      ptr % FValues(1:n,1:m,1) = FValues(1:n,1:m)

      IF ( PRESENT(Proc) ) THEN
        ptr % PROCEDURE = Proc
      END IF

      IF ( PRESENT( Cvalue ) ) THEN
         ptr % CValue = CValue
         ptr % TYPE  = LIST_TYPE_CONSTANT_TENSOR_STR
      END IF

      ptr % NameLen = StringToLowerCase( ptr % Name,Name )
    END SUBROUTINE ListAddConstRealArray
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Adds a real array where the components are linearly dependent.
!------------------------------------------------------------------------------
    SUBROUTINE ListAddDepRealArray(List,Name,DependName, &
               N,TValues,N1,N2,FValues,Proc,Cvalue)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name,DependName
     CHARACTER(LEN=*), OPTIONAL :: Cvalue
     INTEGER :: N,N1,N2
     REAL(KIND=dp) :: FValues(:,:,:)
     REAL(KIND=dp) :: TValues(N)
     INTEGER(KIND=AddrInt), OPTIONAL :: Proc
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------

     ptr => ListAdd( List, Name )
     IF ( PRESENT(Proc) ) ptr % PROCEDURE = Proc

     ALLOCATE( ptr % FValues(n1,n2,N),ptr % TValues(N) )

     ptr % TValues = TValues(1:N)
     ptr % FValues = FValues(1:n1,1:n2,1:N)
     ptr % TYPE = LIST_TYPE_VARIABLE_TENSOR

     IF ( PRESENT( Cvalue ) ) THEN
        ptr % CValue = CValue
        ptr % TYPE = LIST_TYPE_VARIABLE_TENSOR_STR
     END IF

     ptr % NameLen = StringToLowerCase( ptr % Name,Name )
     ptr % DepNameLen = StringToLowerCase( ptr % DependName,DependName )
!------------------------------------------------------------------------------
   END SUBROUTINE ListAddDepRealArray
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Gets a integer value from the list.
!------------------------------------------------------------------------------
   RECURSIVE FUNCTION ListGetInteger( List,Name,Found,minv,maxv,UnfoundFatal ) RESULT(L)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     INTEGER :: L
     LOGICAL, OPTIONAL :: Found, UnfoundFatal
     INTEGER, OPTIONAL :: minv,maxv
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------
     L = 0
     ptr => ListFind(List,Name,Found)
     IF (.NOT.ASSOCIATED(ptr) ) THEN
       IF(PRESENT(UnfoundFatal)) THEN
         IF(UnfoundFatal) THEN
           WRITE(Message, '(A,A)') "Failed to find integer: ",Name
           CALL Fatal("ListGetInteger", Message)
         END IF
       END IF
       RETURN
     END IF

     IF ( ptr % PROCEDURE /= 0 ) THEN
       CALL ListPushActiveName(Name)
       L = ExecIntFunction( ptr % PROCEDURE, CurrentModel )
       CALL ListPopActiveName()
     ELSE
       IF ( .NOT. ASSOCIATED(ptr % IValues) ) THEN
         WRITE(Message,*) 'Value type for property [', TRIM(Name), &
                 '] not used consistently.'
         CALL Fatal( 'ListGetInteger', Message )
         RETURN
       END IF

       L = ptr % IValues(1)
     END IF

     IF ( PRESENT( minv ) ) THEN
       IF ( L < minv ) THEN
         WRITE( Message, '(A,I0,A,I0)') 'Given value ',L,' for property: ['//TRIM(Name)//& 
             '] smaller than given minimum: ', minv
         CALL Fatal( 'ListGetInteger', Message )
       END IF
     END IF

     IF ( PRESENT( maxv ) ) THEN
        IF ( L > maxv ) THEN
          WRITE( Message, '(A,I0,A,I0)') 'Given value ',L,' for property: ['//TRIM(Name)//& 
              '] larger than given maximum: ', maxv
          CALL Fatal( 'ListGetInteger', Message )
        END IF
     END IF
!------------------------------------------------------------------------------
   END FUNCTION ListGetInteger
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Gets a integer array from the list.
!------------------------------------------------------------------------------
   RECURSIVE FUNCTION ListGetIntegerArray( List,Name,Found,UnfoundFatal ) RESULT( IValues )
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*)  :: Name
     LOGICAL, OPTIONAL :: Found, UnfoundFatal
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
     INTEGER :: i,n
     INTEGER, POINTER :: IValues(:)
!------------------------------------------------------------------------------
     NULLIFY( IValues )
     ptr => ListFind(List,Name,Found)
     IF (.NOT.ASSOCIATED(ptr) ) THEN
       IF(PRESENT(UnfoundFatal)) THEN
         IF(UnfoundFatal) THEN
           WRITE(Message, '(A,A)') "Failed to find integer array: ",Name
           CALL Fatal("ListGetInteger", Message)
         END IF
       END IF
       RETURN
     END IF

     IF ( .NOT. ASSOCIATED(ptr % IValues) ) THEN
       WRITE(Message,*) 'VALUE TYPE for property [', TRIM(Name), &
               '] not used consistently.'
       CALL Fatal( 'ListGetIntegerArray', Message )
       RETURN
     END IF

     n = SIZE(ptr % IValues)
     IValues => Ptr % IValues(1:n)

     IF ( ptr % PROCEDURE /= 0 ) THEN
       CALL ListPushActiveName(Name)
       IValues = 0
       DO i=1,N
         Ivalues(i) = ExecIntFunction( ptr % PROCEDURE, CurrentModel )
       END DO
       CALL ListPopActiveName()
     END IF
!------------------------------------------------------------------------------
   END FUNCTION ListGetIntegerArray
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Gets a logical value from the list, if not found return False.
!------------------------------------------------------------------------------
   RECURSIVE FUNCTION ListGetLogical( List,Name,Found,UnfoundFatal ) RESULT(L)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL :: L
     LOGICAL, OPTIONAL :: Found, UnfoundFatal
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------
     L = .FALSE.
     ptr => ListFind(List,Name,Found)
     IF (.NOT.ASSOCIATED(ptr) ) THEN
       IF(PRESENT(UnfoundFatal)) THEN
         IF(UnfoundFatal) THEN
           WRITE(Message, '(A,A)') "Failed to find logical: ",Name
           CALL Fatal("ListGetInteger", Message)
         END IF
       END IF
       RETURN
     END IF
     L = ptr % Lvalue
!------------------------------------------------------------------------------
   END FUNCTION ListGetLogical
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Gets a string from the list by its name, if not found return empty string.
!------------------------------------------------------------------------------
   RECURSIVE FUNCTION ListGetString( List,Name,Found,UnfoundFatal ) RESULT(S)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found,UnfoundFatal
     CHARACTER(LEN=MAX_NAME_LEN) :: S
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------
     S = ' '
     ptr => ListFind(List,Name,Found)
     IF (.NOT.ASSOCIATED(ptr) ) THEN
       IF(PRESENT(UnfoundFatal)) THEN
         IF(UnfoundFatal) THEN
           WRITE(Message, '(A,A)') "Failed to find string: ",Name
           CALL Fatal("ListGetInteger", Message)
         END IF
       END IF
       RETURN
     END IF
     S = ptr % Cvalue
!------------------------------------------------------------------------------
   END FUNCTION ListGetString
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Get a constant real from the list by its name. 
!------------------------------------------------------------------------------
   RECURSIVE FUNCTION ListGetConstReal( List,Name,Found,x,y,z,minv,maxv,UnfoundFatal ) RESULT(F)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     REAL(KIND=dp) :: F
     LOGICAL, OPTIONAL :: Found,UnfoundFatal
     REAL(KIND=dp), OPTIONAL :: x,y,z
     REAL(KIND=dp), OPTIONAL :: minv,maxv
!------------------------------------------------------------------------------
     TYPE(Variable_t), POINTER :: Variable
     TYPE(ValueListEntry_t), POINTER :: ptr
     REAL(KIND=dp) :: xx,yy,zz
     INTEGER :: i,j,k,n
     CHARACTER(LEN=MAX_NAME_LEN) :: cmd,tmp_str
!------------------------------------------------------------------------------
     F = 0.0_dp

     ptr => ListFind(List,Name,Found)
     IF (.NOT.ASSOCIATED(ptr) ) THEN
       IF(PRESENT(UnfoundFatal)) THEN
         IF(UnfoundFatal) THEN
           WRITE(Message, '(A,A)') "Failed to find ConstReal: ",Name
           CALL Fatal("ListGetInteger", Message)
         END IF
       END IF
       RETURN
     END IF

     SELECT CASE(ptr % TYPE)

     CASE( LIST_TYPE_CONSTANT_SCALAR )

       IF ( .NOT. ASSOCIATED(ptr % FValues) ) THEN
         WRITE(Message,*) 'VALUE TYPE for property [', TRIM(Name), &
             '] not used consistently.'
         CALL Fatal( 'ListGetConstReal', Message )
       END IF
       F = ptr % Coeff * ptr % Fvalues(1,1,1)

     CASE( LIST_TYPE_CONSTANT_SCALAR_STR )

        cmd = ptr % CValue
        k = LEN_TRIM( ptr % CValue )
        CALL matc( cmd, tmp_str, k )
        READ( tmp_str(1:k), * ) F
        F = ptr % Coeff * F

     CASE( LIST_TYPE_CONSTANT_SCALAR_PROC )

       IF ( ptr % PROCEDURE == 0 ) THEN
         WRITE(Message,*) 'VALUE TYPE for property [', TRIM(Name), &
             '] not used consistently.'
         CALL Fatal( 'ListGetConstReal', Message )
       END IF

       xx = 0.0_dp
       yy = 0.0_dp
       zz = 0.0_dp
       IF ( PRESENT(x) ) xx = x
       IF ( PRESENT(y) ) yy = y
       IF ( PRESENT(z) ) zz = z
       CALL ListPushActiveName(Name)
       F = Ptr % Coeff * &
           ExecConstRealFunction( ptr % PROCEDURE,CurrentModel,xx,yy,zz )
       CALL ListPopActiveName()

     END SELECT

     IF ( PRESENT( minv ) ) THEN
        IF ( F < minv ) THEN
           WRITE( Message, *) 'Given VALUE ', F, ' for property: ', '[', TRIM(Name),']', &
               ' smaller than given minimum: ', minv
           CALL Fatal( 'ListGetInteger', Message )
        END IF
     END IF

     IF ( PRESENT( maxv ) ) THEN
        IF ( F > maxv ) THEN
           WRITE( Message, *) 'Given VALUE ', F, ' for property: ', '[', TRIM(Name),']', &
               ' larger than given maximum: ', maxv
           CALL Fatal( 'ListGetInteger', Message )
        END IF
     END IF
!------------------------------------------------------------------------------
   END FUNCTION ListGetConstReal
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Returns a scalar real value, that may depend on other scalar values such as 
!> time or timestep size etc.
!------------------------------------------------------------------------------
  RECURSIVE FUNCTION ListGetCReal( List, Name, Found, UnfoundFatal) RESULT(s)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found,UnfoundFatal
     INTEGER, TARGET :: Dnodes(1)
     INTEGER, POINTER :: NodeIndexes(:)

     REAL(KIND=dp) :: s
     REAL(KIND=dp) :: x(1)
     TYPE(Element_t), POINTER :: Element

     INTEGER :: n, istat

     IF ( PRESENT( Found ) ) Found = .FALSE.

     NodeIndexes => Dnodes
     n = 1
     NodeIndexes(n) = 1

     x = 0.0_dp
     IF ( ASSOCIATED(List % head) ) THEN
        IF ( PRESENT( Found ) ) THEN
           x(1:n) = ListGetReal( List, Name, n, NodeIndexes, Found, UnfoundFatal=UnfoundFatal )
        ELSE
           x(1:n) = ListGetReal( List, Name, n, NodeIndexes, UnfoundFatal=UnfoundFatal)
        END IF
     END IF
     s = x(1)
!------------------------------------------------------------------------------
  END FUNCTION ListGetCReal
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Returns a scalar real value, that may depend on other scalar values such as 
!> time or timestep size etc.
!------------------------------------------------------------------------------
  RECURSIVE FUNCTION ListGetRealAtNode( List, Name, Node, Found, UnfoundFatal ) RESULT(s)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*)  :: Name
     INTEGER :: Node
     LOGICAL, OPTIONAL :: Found, UnfoundFatal
     REAL(KIND=dp) :: s
!-----------------------------------------------------------------------------
     INTEGER, TARGET, SAVE :: Dnodes(1)
     INTEGER, POINTER :: NodeIndexes(:) 
     REAL(KIND=dp) :: x(1)
     INTEGER, PARAMETER :: one = 1

     IF ( PRESENT( Found ) ) Found = .FALSE.

     IF ( ASSOCIATED(List % Head) ) THEN
       NodeIndexes => Dnodes
       NodeIndexes(one) = Node
       
       x(1:one) = ListGetReal( List, Name, one, NodeIndexes, Found, UnfoundFatal=UnfoundFatal)
       s = x(one)
     ELSE
       s = 0.0_dp
     END IF

!------------------------------------------------------------------------------
  END FUNCTION ListGetRealAtNode
!------------------------------------------------------------------------------

#define MAX_FNC 32

!------------------------------------------------------------------------------
  SUBROUTINE ListParseStrToValues( str, slen, ind, name, T, count, AllGlobal )
!------------------------------------------------------------------------------
     CHARACTER(LEN=*) :: str, name
     REAL(KIND=dp)  :: T(:)
     INTEGER :: slen, count, ind
     LOGICAL :: AllGlobal
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
     INTEGER :: i,j,k,n,k1,l,l0,l1
     TYPE(Variable_t), POINTER :: Variable, CVar

     AllGlobal = .TRUE.
     count=0
     l0=1
     IF(slen<=0) RETURN

     DO WHILE( .TRUE. )
       DO WHILE( str(l0:l0) == ' ' )
         l0 = l0 + 1
         IF ( l0 > slen ) EXIT
       END DO
       IF ( l0 > slen ) EXIT

       l1 = INDEX( str(l0:slen),',')
       IF ( l1 > 0 ) THEN
         l1=l0+l1-2
       ELSE
         l1=slen
       END IF

       IF ( str(l0:l1) /= 'coordinate' ) THEN
         Variable => VariableGet( CurrentModel % Variables,TRIM(str(l0:l1)) )
         IF ( .NOT. ASSOCIATED( Variable ) ) THEN
           WRITE( Message, * ) 'Can''t find INDEPENDENT variable:[', &
               TRIM(str(l0:l1)),']' // &
               'for dependent variable:[', TRIM(Name),']'
           CALL Fatal( 'ListGetReal', Message )
         END IF
         IF( SIZE( Variable % Values ) > 1 ) AllGlobal = .FALSE.
       ELSE
         AllGlobal = .FALSE.
         Variable => VariableGet( CurrentModel % Variables,'Coordinate 1' )
       END IF
       
       k1 = ind
       IF ( Variable % TYPE == Variable_on_nodes_on_elements ) THEN
         Element => CurrentModel % CurrentElement
         IF ( ASSOCIATED(Element) ) THEN
           IF ( ASSOCIATED(Element % DGIndexes) ) THEN
             n = Element % TYPE % NumberOfNodes
             IF ( SIZE(Element % DGIndexes)==n ) THEN
               DO i=1,n
                 IF ( Element % NodeIndexes(i)==ind ) THEN
                   k1 = Element % DGIndexes(i)
                   EXIT
                 END IF
               END DO
             END IF
           END IF
         END IF
       END IF
       IF ( ASSOCIATED(Variable % Perm) ) k1 = Variable % Perm(k1)

       IF ( k1>0 .AND. k1<=SIZE(Variable % Values) ) THEN
         IF ( str(l0:l1) == 'coordinate' ) THEN
           CVar => VariableGet( CurrentModel % Variables, 'Coordinate 1' )
           count = count + 1
           T(1) = CVar % Values(k1)
           CVar => VariableGet( CurrentModel % Variables, 'Coordinate 2' )
           count = count + 1
           T(2) = CVar % Values(k1)
           CVar => VariableGet( CurrentModel % Variables, 'Coordinate 3' )
           count = count + 1
           T(3) = CVar % Values(k1)
         ELSE
           IF ( Variable % DOFs == 1 ) THEN
              count = count + 1
              T(count) = Variable % Values(k1)
           ELSE
              DO l=1,Variable % DOFs
                 count = count + 1
                 T(count) = Variable % Values(Variable % DOFs*(k1-1)+l)
              END DO
           END IF
         END IF
       ELSE

         count = count + 1
         IF ( ASSOCIATED(Variable % Perm) ) THEN
            T(count) = HUGE(1.0_dp)
            EXIT
         ELSE
            T(count) = Variable % Values(1)
         END IF
       END IF

       l0 = l1+2
       IF ( l0 > slen ) EXIT
     END DO

!------------------------------------------------------------------------------
  END SUBROUTINE ListParseStrToValues
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION ListCheckAllGlobal( List, name ) RESULT ( AllGlobal )
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: name
     LOGICAL :: AllGlobal
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
     TYPE(Element_t), POINTER :: Element
     INTEGER :: ind,i,j,k,n,k1,l,l0,l1
     TYPE(Variable_t), POINTER :: Variable, CVar
     INTEGER :: slen

     AllGlobal = .TRUE.

     IF(.NOT.ASSOCIATED(List)) RETURN
     ptr => List % Head
     IF(.NOT.ASSOCIATED(ptr)) RETURN

     slen = ptr % DepNameLen

     IF( ptr % TYPE == LIST_TYPE_CONSTANT_SCALAR_STR ) THEN
       RETURN

     ELSE IF( ptr % TYPE == LIST_TYPE_CONSTANT_SCALAR .OR. & 
         ptr % TYPE == LIST_TYPE_VARIABLE_SCALAR .OR. &
         ptr % TYPE == LIST_TYPE_VARIABLE_SCALAR_STR ) THEN

       l0 = 1
       DO WHILE( .TRUE. )
         DO WHILE( ptr % DependName(l0:l0) == ' ' )
           l0 = l0 + 1
         END DO
         IF ( l0 > slen ) EXIT
         
         l1 = INDEX( ptr % DependName(l0:slen),',')
         IF ( l1 > 0 ) THEN
           l1=l0+l1-2
         ELSE
           l1=slen
         END IF
         
         IF ( ptr % DependName(l0:l1) /= 'coordinate' ) THEN
           Variable => VariableGet( CurrentModel % Variables,TRIM(ptr % DependName(l0:l1)) )
           IF ( .NOT. ASSOCIATED( Variable ) ) THEN
             WRITE( Message, * ) 'Can''t find INDEPENDENT variable:[', &
                 TRIM(ptr % DependName(l0:l1)),']' // &
                 'for dependent variable:[', TRIM(Name),']'
             CALL Fatal( 'ListGetReal', Message )
           END IF

           IF( SIZE( Variable % Values ) > 1 ) THEN
             AllGlobal = .FALSE.
             RETURN
           END IF
         ELSE
           AllGlobal = .FALSE.
           RETURN
         END IF

         l0 = l1+2
         IF ( l0 > slen ) EXIT

       END DO

     ELSE
       CALL Fatal('ListCheckAllGlobal','Unknown type: '//TRIM(I2S(ptr % TYPE)))
     END IF

!------------------------------------------------------------------------------
   END FUNCTION ListCheckAllGlobal
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Gets a real valued parameter in each node of an element.
!------------------------------------------------------------------------------
   RECURSIVE FUNCTION ListGetReal( List,Name,N,NodeIndexes,Found,minv,maxv,UnfoundFatal ) RESULT(F)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*)  :: Name
     INTEGER :: N,NodeIndexes(:)
     REAL(KIND=dp)  :: F(N)
     LOGICAL, OPTIONAL :: Found, UnfoundFatal
     REAL(KIND=dp), OPTIONAL :: minv,maxv
!------------------------------------------------------------------------------
     TYPE(Variable_t), POINTER :: Variable, CVar, TVar
     TYPE(ValueListEntry_t), POINTER :: ptr
     REAL(KIND=dp) :: T(MAX_FNC)
     INTEGER :: i,j,k,k1,l,l0,l1,lsize
     CHARACTER(LEN=MAX_NAME_LEN) ::  cmd, tmp_str
     LOGICAL :: AllGlobal
     ! INTEGER :: TID, OMP_GET_THREAD_NUM
!------------------------------------------------------------------------------
     ! TID = 0
     ! !$ TID=OMP_GET_THREAD_NUM()
     F = 0.0_dp
     ptr => ListFind(List,Name,Found)
     IF (.NOT.ASSOCIATED(ptr) ) THEN
       IF(PRESENT(UnfoundFatal)) THEN
         IF(UnfoundFatal) THEN
           WRITE(Message, '(A,A)') "Failed to find real: ",Name
           CALL Fatal("ListGetInteger", Message)
         END IF
       END IF
       RETURN
     END IF


     SELECT CASE(ptr % TYPE)

     CASE( LIST_TYPE_CONSTANT_SCALAR )

       IF ( .NOT. ASSOCIATED(ptr % FValues) ) THEN
         WRITE(Message,*) 'VALUE TYPE for property [', TRIM(Name), &
             '] not used consistently.'
         CALL Fatal( 'ListGetReal', Message )
         RETURN
       END IF
       F = ptr % Coeff * ptr % Fvalues(1,1,1)

     
     CASE( LIST_TYPE_VARIABLE_SCALAR )

       CALL ListPushActiveName(Name)
       DO i=1,n
         k = NodeIndexes(i)
         CALL ListParseStrToValues( Ptr % DependName, Ptr % DepNameLen, k, Name, T, j, AllGlobal)

         IF ( .NOT. ANY( T(1:j)==HUGE(1.0_dp) ) ) THEN
           IF ( ptr % PROCEDURE /= 0 ) THEN
             F(i) = ptr % Coeff * &
                 ExecRealFunction( ptr % PROCEDURE,CurrentModel, k, T )
           ELSE
             IF ( .NOT. ASSOCIATED(ptr % FValues) ) THEN
               WRITE(Message,*) 'VALUE TYPE for property [', TRIM(Name), &
                       '] not used consistently.'
               CALL Fatal( 'ListGetReal', Message )
               RETURN
             END IF
             F(i) = ptr % Coeff * &
                 InterpolateCurve( ptr % TValues,ptr % FValues(1,1,:), &
                 T(1), ptr % CubicCoeff )
             IF( AllGlobal) THEN
               F(2:n) = F(1)
               EXIT
             END IF
           END IF
         END IF
       END DO
       CALL ListPopActiveName()


     CASE( LIST_TYPE_CONSTANT_SCALAR_STR )
         TVar => VariableGet( CurrentModel % Variables, 'Time' ) 
         WRITE( cmd, '(a,e15.8)' ) 'st = ', TVar % Values(1)
         k = LEN_TRIM(cmd)
         CALL matc( cmd, tmp_str, k )

         cmd = ptr % CValue
         k = LEN_TRIM(cmd)
         CALL matc( cmd, tmp_str, k )
         READ( tmp_str(1:k), * ) F(1)
         F(1) = ptr % Coeff * F(1)
         F(2:n) = F(1)

     CASE( LIST_TYPE_VARIABLE_SCALAR_STR )

       TVar => VariableGet( CurrentModel % Variables, 'Time' ) 
       WRITE( cmd, * ) 'tx=0; st = ', TVar % Values(1)
       k = LEN_TRIM(cmd)
       CALL matc( cmd, tmp_str, k )

       DO i=1,n
         k = NodeIndexes(i)
         CALL ListParseStrToValues( Ptr % DependName, Ptr % DepNameLen, k, Name, T, j, AllGlobal)

         IF ( .NOT. ANY( T(1:j)==HUGE(1.0_dp) ) ) THEN
           DO l=1,j
             WRITE( cmd, * ) 'tx('//TRIM(i2s(l-1))//')=', T(l)
             k1 = LEN_TRIM(cmd)
             CALL matc( cmd, tmp_str, k1 )
           END DO

           cmd = ptr % CValue
           k1 = LEN_TRIM(cmd)
           CALL matc( cmd, tmp_str, k1 )
           READ( tmp_str(1:k1), * ) F(i)
           F(i) = Ptr % Coeff * F(i)
         END IF

         IF( AllGlobal ) THEN
           F(2:n) = F(1)
           EXIT
         END IF

       END DO

     CASE( LIST_TYPE_CONSTANT_SCALAR_PROC )

       IF ( ptr % PROCEDURE == 0 ) THEN
         WRITE(Message,*) 'VALUE TYPE for property [', TRIM(Name), &
             '] not used consistently.'
         CALL Fatal( 'ListGetReal', Message )
         RETURN
       END IF

       CALL ListPushActiveName(name)
       DO i=1,n
         F(i) = Ptr % Coeff * &
             ExecConstRealFunction( ptr % PROCEDURE,CurrentModel, &
             CurrentModel % Mesh % Nodes % x( NodeIndexes(i) ), &
             CurrentModel % Mesh % Nodes % y( NodeIndexes(i) ), &
             CurrentModel % Mesh % Nodes % z( NodeIndexes(i) ) )
       END DO
       CALL ListPopActiveName()

     END SELECT

     IF ( PRESENT( minv ) ) THEN
        IF ( MINVAL(F(1:n)) < minv ) THEN
           WRITE( Message,*) 'Given VALUE ', MINVAL(F(1:n)), ' for property: ', '[', TRIM(Name),']', &
               ' smaller than given minimum: ', minv
           CALL Fatal( 'ListGetReal', Message )
        END IF
     END IF

     IF ( PRESENT( maxv ) ) THEN
        IF ( MAXVAL(F(1:n)) > maxv ) THEN
           WRITE( Message,*) 'Given VALUE ', MAXVAL(F(1:n)), ' for property: ', '[', TRIM(Name),']', &
               ' larger than given maximum ', maxv
           CALL Fatal( 'ListGetReal', Message )
        END IF
     END IF
   END FUNCTION ListGetReal
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Gets a real valued parameter in one single point with value x.
!> Optionally also computes the derivative at that point. 
!> Note that this uses same logical on sif file as ListGetReal 
!> but the variable is just a dummy as the dependent function is 
!> assumed to be set inside the code. This should be used with caution
!> is it sets some confusing limitations to the user. The main limitation
!> is the use of just one dependent variable. 
!------------------------------------------------------------------------------
   RECURSIVE FUNCTION ListGetFun( List,Name,x,Found,minv,maxv,dFdx,eps ) RESULT(F)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     REAL(KIND=dp), OPTIONAL :: x     
     REAL(KIND=dp) :: f
     CHARACTER(LEN=*), OPTIONAL  :: Name
     LOGICAL, OPTIONAL :: Found
     REAL(KIND=dp), OPTIONAL :: minv,maxv
     REAL(KIND=dp), OPTIONAL :: dFdx, eps
!------------------------------------------------------------------------------
     TYPE(Variable_t), POINTER :: Variable, CVar, TVar
     TYPE(ValueListEntry_t), POINTER :: ptr
     REAL(KIND=dp) :: T(1)
     INTEGER :: i,j,k,k1,l,l0,l1,lsize
     CHARACTER(LEN=MAX_NAME_LEN) ::  cmd, tmp_str
     LOGICAL :: AllGlobal
     REAL(KIND=dp) :: xeps, F2, F1
!------------------------------------------------------------------------------

     F = 0.0_dp
     IF( PRESENT( Name ) ) THEN
       ptr => ListFind(List,Name,Found)
       IF ( .NOT.ASSOCIATED(ptr) ) RETURN
     ELSE
       IF(.NOT.ASSOCIATED(List)) RETURN
       ptr => List % Head
       IF ( .NOT.ASSOCIATED(ptr) ) THEN
         CALL Warn('ListGetFun','List entry not associated')
         RETURN
       END IF
     END IF

     ! Node number not applicable, hence set to zero
     k = 0
     T(1) = x

     SELECT CASE(ptr % TYPE)

     CASE( LIST_TYPE_CONSTANT_SCALAR )

       IF ( .NOT. ASSOCIATED(ptr % FValues) ) THEN
         WRITE(Message,*) 'VALUE TYPE for property [', TRIM(Name), &
             '] not used consistently.'
         CALL Fatal( 'ListGetReal', Message )
         RETURN
       END IF
       F = ptr % Coeff * ptr % Fvalues(1,1,1)
       IF( PRESENT( dFdx ) ) THEN
         dFdx = 0.0_dp
       END IF


     CASE( LIST_TYPE_VARIABLE_SCALAR )

       IF ( ptr % PROCEDURE /= 0 ) THEN
         CALL ListPushActiveName(name)
         F = ExecRealFunction( ptr % PROCEDURE,CurrentModel, k, T(1) )

         ! Compute derivative at the point if requested
         ! Numerical central difference scheme is used for accuracy. 
         IF( PRESENT( dFdx ) ) THEN
           IF( PRESENT( eps ) ) THEN
             xeps = eps
           ELSE
             xeps = 1.0e-8
           END IF
           T(1) = x - xeps
           F1 = ExecRealFunction( ptr % PROCEDURE,CurrentModel, k, T(1) )
           T(1) = x + xeps
           F2 = ExecRealFunction( ptr % PROCEDURE,CurrentModel, k, T(1) )
           dFdx = ( F2 - F1 ) / (2*xeps)
         END IF
         CALL ListPopActiveName()
       ELSE
         IF ( .NOT. ASSOCIATED(ptr % FValues) ) THEN
           WRITE(Message,*) 'VALUE TYPE for property [', TRIM(Name), &
               '] not used consistently.'
           CALL Fatal( 'ListGetFun', Message )
           RETURN
         END IF
         F = InterpolateCurve( ptr % TValues,ptr % FValues(1,1,:), &
             x, ptr % CubicCoeff )
         ! Compute the derivative symbolically from the table values. 
         IF( PRESENT( dFdx ) ) THEN
           dFdx = DerivateCurve(ptr % TValues,ptr % FValues(1,1,:), &
               x, ptr % CubicCoeff )
         END IF
       END IF


     CASE( LIST_TYPE_VARIABLE_SCALAR_STR )
       WRITE( cmd, * ) 'tx=', X
       k1 = LEN_TRIM(cmd)
       CALL matc( cmd, tmp_str, k1 )
       
       cmd = ptr % CValue
       k1 = LEN_TRIM(cmd)
       CALL matc( cmd, tmp_str, k1 )
       READ( tmp_str(1:k1), * ) F

       ! This is really expensive. 
       ! For speed also one sided difference could be considered. 
       IF( PRESENT( dFdx ) ) THEN
         IF( PRESENT( eps ) ) THEN
           xeps = eps
         ELSE
           xeps = 1.0e-8
         END IF
         
         WRITE( cmd, * ) 'tx=', x-xeps
         k1 = LEN_TRIM(cmd)
         CALL matc( cmd, tmp_str, k1 )
         
         cmd = ptr % CValue
         k1 = LEN_TRIM(cmd)
         CALL matc( cmd, tmp_str, k1 )
         READ( tmp_str(1:k1), * ) F1
         
         WRITE( cmd, * ) 'tx=', x+xeps
         k1 = LEN_TRIM(cmd)
         CALL matc( cmd, tmp_str, k1 )
         
         cmd = ptr % CValue
         k1 = LEN_TRIM(cmd)
         CALL matc( cmd, tmp_str, k1 )
         READ( tmp_str(1:k1), * ) F2

         dFdx = (F2-F1) / (2*xeps)
       END IF

     CASE DEFAULT
       CALL Fatal('ListGetFun','LIST_TYPE not implemented!')

     END SELECT

     IF ( PRESENT( minv ) ) THEN
        IF ( F < minv ) THEN
           WRITE( Message,*) 'Given VALUE ', F, ' for property: ', '[', TRIM(Name),']', &
               ' smaller than given minimum: ', minv
           CALL Fatal( 'ListGetFun', Message )
        END IF
     END IF

     IF ( PRESENT( maxv ) ) THEN
        IF ( F > maxv ) THEN
           WRITE( Message,*) 'Given VALUE ', F, ' for property: ', '[', TRIM(Name),']', &
               ' larger than given maximum ', maxv
           CALL Fatal( 'ListGetFun', Message )
        END IF
     END IF

   END FUNCTION ListGetFun
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Gets a real valued parameter in the Gaussian integration point defined 
!> by the local basis function. To speed up things there is a handle associated
!> to the given keyword (Name). Here the values are first evaluated at the 
!> nodal points and then using basis functions estimated at the 
!> gaussian integration points. 
!------------------------------------------------------------------------------
   RECURSIVE FUNCTION ListGetRealAtIp( Handle,List,Basis,Name,Found,&
       Element,minv,maxv ) RESULT(val)
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     TYPE(ValueList_t), POINTER :: List
     REAL(KIND=dp), OPTIONAL :: Basis(:)
     CHARACTER(LEN=*), OPTIONAL  :: Name
     LOGICAL, OPTIONAL :: Found
     TYPE(Element_t), POINTER, OPTIONAL :: Element
     REAL(KIND=dp), OPTIONAL :: minv,maxv
     REAL(KIND=dp)  :: val
!------------------------------------------------------------------------------
     TYPE(Variable_t), POINTER :: Variable, CVar, TVar
     TYPE(ValueListEntry_t), POINTER :: ptr
     INTEGER, POINTER :: NodeIndexes(:)
     REAL(KIND=dp) :: T(MAX_FNC),x,y,z
     REAL(KIND=dp), POINTER :: F(:)
     REAL(KIND=dp), POINTER :: ParF(:,:)
     INTEGER :: i,j,k,k1,l,l0,l1,lsize,n
     CHARACTER(LEN=MAX_NAME_LEN) ::  cmd, tmp_str
     LOGICAL :: AllGlobal, GotIt
     TYPE(Element_t), POINTER :: PElement
!------------------------------------------------------------------------------

     ! If the value is known to be globally constant return it asap.
     ! This can only be the case if the 'ListInitRealAtIp' has been called prior 
     ! to this subroutine.
     IF( Handle % ConstantEverywhere ) THEN
       val = Handle % ConstantValue
       RETURN
     END IF

     ! Set the default value 
     val = 0.0_dp

     ! If the provided list is the same as last time, also the keyword will
     ! be sitting at the same place, otherwise find it in the new list
     IF( .NOT. ASSOCIATED(List) ) THEN
       IF(PRESENT(Found)) Found = .FALSE.
       RETURN
     ELSE IF( .NOT. ASSOCIATED( List % Head ) ) THEN
       IF(PRESENT(Found)) Found = .FALSE.
       RETURN
     ELSE IF( ASSOCIATED( Handle % List, List ) ) THEN
       ptr => Handle % ptr % head
       IF(PRESENT(Found)) Found = .TRUE.
       IF( Handle % ConstantInList ) THEN
         val = Handle % Values(1)
         RETURN
       END IF
     ELSE 
       IF( PRESENT( Name ) ) THEN
         ptr => ListFind(List,Name,Found)
       ELSE IF(.NOT. Handle % Initialized ) THEN
         CALL Fatal('ListGetRealAtIp','Handle must be initialized if name is not given!')
       ELSE
         ptr => ListFind(List,Handle % Name,Found)
       END IF
  
       Handle % List => List
       IF(.NOT.ASSOCIATED(Handle % Ptr)) &
           Handle % Ptr => ListAllocate()
       Handle % Ptr % Head => ptr
       
       ! Check whether the keyword should be evaluated at integration point directly
       IF( ListGetLogical( List, TRIM( Name )//' At IP',GotIt ) ) THEN
         CALL Info('ListGetRealAtIp','Evaluating at ip: '//TRIM(Name),Level=10 )

         ! It does not make sense to evaluate global variables at IP
         IF( PRESENT( Name ) ) THEN
           Handle % EvaluateAtIp = .NOT. &
               ListCheckAllGlobal(  Handle % Ptr, Name )
         ELSE 
           Handle % EvaluateAtIp = .NOT. &
               ListCheckAllGlobal(  Handle % Ptr, Handle % Name )
         END IF
         IF( .NOT. Handle % EvaluateAtIp ) THEN
           CALL Info('ListGetRealAtIp','Evaluation at ip not needed!',Level=10)
         END IF
       ELSE
         Handle % EvaluateAtIp = .FALSE.
       END IF

       IF ( .NOT.ASSOCIATED(ptr) ) RETURN
       Handle % ConstantInList = .FALSE.
     END IF

     ! If there is no pointer return default value zero
     IF ( .NOT. ASSOCIATED(ptr) ) RETURN

     ! Get the pointer to the element
     IF( PRESENT( Element) ) THEN
       PElement => Element
     ELSE
       PElement => CurrentModel % CurrentElement
     END IF


     ! Either evaluate parameter directly at IP, 
     ! or first at nodes and then using basis functions at IP.
     ! The later is the default. 
     !------------------------------------------------------------------
     IF( Handle % EvaluateAtIp ) THEN

       IF(.NOT. PRESENT(Basis)) THEN
         CALL Fatal('ListGetRealAtIp','Parameter > Basis < is required!')
       END IF
       
       ! If we get back to the same element than last time use the data already 
       ! retrieved. If the element is new then get the data in every node of the 
       ! current element, or only in the 1st node if it is constant. 
       
       IF( ASSOCIATED( PElement, Handle % Element ) ) THEN
         n = Handle % Element % TYPE % NumberOfNodes 
         NodeIndexes => PElement % NodeIndexes
         ParF => Handle % ParValues
       ELSE
         
         IF( .NOT. Handle % AllocationsDone ) THEN
           n = CurrentModel % Mesh % MaxElementNodes
           ALLOCATE( Handle % Values(n) )
           Handle % Values = 0.0_dp
           ALLOCATE( Handle % ParValues(MAX_FNC,n) )
           Handle % ParValues = 0.0_dp
           Handle % AllocationsDone = .TRUE.
         END IF
         
         Handle % Element => PElement
         n = PElement % TYPE % NumberOfNodes 
         NodeIndexes => PElement % NodeIndexes
         
         IF( ptr % TYPE == LIST_TYPE_VARIABLE_SCALAR .OR. &
             ptr % TYPE == LIST_TYPE_VARIABLE_SCALAR_STR ) THEN

           ! These might not have been initialized if this is has mixed evaluation strategies           
           IF(.NOT. ASSOCIATED( Handle % ParValues )) THEN
             ALLOCATE( Handle % ParValues(MAX_FNC,CurrentModel % Mesh % MaxElementNodes) )
             Handle % ParValues = 0.0_dp
           END IF
           
           DO i=1,n
             k = NodeIndexes(i)
             CALL ListParseStrToValues( Ptr % DependName, Ptr % DepNameLen, k, Name, T, j, AllGlobal)
             
             IF( AllGlobal ) THEN
               CALL Fatal('ListGetRealAtIp','Constant lists should not need to be here')
             END IF
             
             Handle % ParNo = j 
             Handle % ParValues(1:j,i) = T(1:j)
           END DO
         END IF

         ParF => Handle % ParValues         
       END IF


       SELECT CASE(ptr % TYPE)
         
       CASE( LIST_TYPE_VARIABLE_SCALAR )
         
         DO j=1,Handle % ParNo 
           T(j) = SUM( Basis(1:n) *  Handle % ParValues(j,1:n) )
         END DO
         
         ! there is no node index, so use zero
         j = 0 
         IF ( ptr % PROCEDURE /= 0 ) THEN
           CALL ListPushActiveName(name)
           val = ExecRealFunction( ptr % PROCEDURE,CurrentModel, j, T )
           CALL ListPopActiveName()
         ELSE
           val = InterpolateCurve( ptr % TValues,ptr % FValues(1,1,:), &
               T(1), ptr % CubicCoeff )
         END IF
         
       CASE( LIST_TYPE_VARIABLE_SCALAR_STR )
         DO j=1,Handle % ParNo 
           T(j) = SUM( Basis(1:n) *  Handle % ParValues(j,1:n) )
         END DO
         ! there is no node index, so use zero
         j = 0 
         
         TVar => VariableGet( CurrentModel % Variables, 'Time' ) 
         WRITE( cmd, * ) 'tx=0; st = ', TVar % Values(1)
         k = LEN_TRIM(cmd)
         CALL matc( cmd, tmp_str, k )
         
         DO l=1,Handle % ParNo
           WRITE( cmd, * ) 'tx('//TRIM(i2s(l-1))//')=', T(l)
           k1 = LEN_TRIM(cmd)
           CALL matc( cmd, tmp_str, k1 )
         END DO
         
         cmd = ptr % CValue
         k1 = LEN_TRIM(cmd)
         CALL matc( cmd, tmp_str, k1 )
         READ( tmp_str(1:k1), * ) val
         
       CASE( LIST_TYPE_CONSTANT_SCALAR_PROC )

         IF ( ptr % PROCEDURE /= 0 ) THEN
           x = SUM( Basis(1:n) * CurrentModel % Mesh % Nodes % x( NodeIndexes(1:n) ) )
           y = SUM( Basis(1:n) * CurrentModel % Mesh % Nodes % y( NodeIndexes(1:n) ) )
           z = SUM( Basis(1:n) * CurrentModel % Mesh % Nodes % z( NodeIndexes(1:n) ) )

           CALL ListPushActiveName(name)
           val = ExecConstRealFunction( ptr % PROCEDURE,CurrentModel,x,y,z)
           CALL ListPopActiveName()
         ELSE
           CALL Fatal('ListGetRealAtIp','Constant scalar evaluation failed at ip!')
         END IF
           
       CASE DEFAULT
         
         CALL Fatal('ListGetRealAtIp','Unknown case for avaluation at ip')
         
       END SELECT
     
     ELSE

       ! If we get back to the same element than last time use the data already 
       ! retrieved. If the element is new then get the data in every node of the 
       ! current element, or only in the 1st node if it is constant. 
       
       IF( ASSOCIATED( PElement, Handle % Element ) ) THEN
         n = Handle % Element % TYPE % NumberOfNodes 
         F => Handle % Values       
       ELSE
         
         IF( .NOT. Handle % AllocationsDone ) THEN
           n = CurrentModel % Mesh % MaxElementNodes
           ALLOCATE( Handle % Values(n) )
           Handle % Values = 0.0_dp
           Handle % AllocationsDone = .TRUE.
         END IF
         
         Handle % Element => PElement
         n = PElement % TYPE % NumberOfNodes 
         NodeIndexes => PElement % NodeIndexes
         F => Handle % Values
         
         SELECT CASE(ptr % TYPE)
           
         CASE( LIST_TYPE_CONSTANT_SCALAR )
           
           Handle % ConstantInList = .TRUE.                      
           IF ( .NOT. ASSOCIATED(ptr % FValues) ) THEN
             WRITE(Message,*) 'Value type for property [', TRIM(Name), &
                 '] not used consistently.'
             CALL Fatal( 'ListGetRealAtIp', Message )
             RETURN
           END IF
           F(1) = ptr % Coeff * ptr % Fvalues(1,1,1)


         CASE( LIST_TYPE_VARIABLE_SCALAR )
           CALL ListPushActiveName(name)
           DO i=1,n
             k = NodeIndexes(i)
             CALL ListParseStrToValues( Ptr % DependName, Ptr % DepNameLen, k, Name, T, j, AllGlobal)
             
             IF ( .NOT. ANY( T(1:j) == HUGE(1.0_dp) ) ) THEN
               IF ( ptr % PROCEDURE /= 0 ) THEN
                 F(i) = ptr % Coeff * &
                     ExecRealFunction( ptr % PROCEDURE,CurrentModel, &
                     NodeIndexes(i), T )              
               ELSE
                 IF ( .NOT. ASSOCIATED(ptr % FValues) ) THEN
                   WRITE(Message,*) 'Value type for property [', TRIM(Name), &
                       '] not used consistently.'
                   CALL Fatal( 'ListGetRealAtIp', Message )
                   RETURN
                 END IF
                 F(i) = ptr % Coeff * &
                     InterpolateCurve( ptr % TValues,ptr % FValues(1,1,:), &
                     T(1), ptr % CubicCoeff )
                 
                 ! If the dependency table includes just global values (such as time) 
                 ! the values will be the same for all element entries.
                 IF( AllGlobal ) THEN
                   Handle % ConstantInList = .TRUE.
                   EXIT
                 END IF
                 
               END IF
             END IF
           END DO
           CALL ListPopActiveName()
           
         CASE( LIST_TYPE_CONSTANT_SCALAR_STR )
           Handle % ConstantInList = .TRUE.
           
           TVar => VariableGet( CurrentModel % Variables, 'Time' ) 
           WRITE( cmd, '(a,e15.8)' ) 'st = ', TVar % Values(1)
           k = LEN_TRIM(cmd)
           CALL matc( cmd, tmp_str, k )
           
           cmd = ptr % CValue
           k = LEN_TRIM(cmd)
           CALL matc( cmd, tmp_str, k )
           READ( tmp_str(1:k), * ) F(1)
           F(1) = ptr % Coeff * F(1) 

         CASE( LIST_TYPE_VARIABLE_SCALAR_STR )
           TVar => VariableGet( CurrentModel % Variables, 'Time' ) 
           WRITE( cmd, * ) 'tx=0; st = ', TVar % Values(1)
           k = LEN_TRIM(cmd)
           CALL matc( cmd, tmp_str, k )
           
           DO i=1,n
             k = NodeIndexes(i)
             CALL ListParseStrToValues( Ptr % DependName, Ptr % DepNameLen, k, Name, T, j, AllGlobal)
             IF ( .NOT. ANY( T(1:j)==HUGE(1.0_dp) ) ) THEN
               DO l=1,j
                 WRITE( cmd, * ) 'tx('//TRIM(i2s(l-1))//')=', T(l)
                 k1 = LEN_TRIM(cmd)
                 CALL matc( cmd, tmp_str, k1 )
               END DO
               
               cmd = ptr % CValue
               k1 = LEN_TRIM(cmd)
               CALL matc( cmd, tmp_str, k1 )
               READ( tmp_str(1:k1), * ) F(i)
               F(i) = ptr % Coeff * F(i)
             END IF
             
             IF( AllGlobal ) THEN
               Handle % ConstantInList = .TRUE.
               EXIT
             END IF            
           END DO

         CASE( LIST_TYPE_CONSTANT_SCALAR_PROC )
           IF ( ptr % PROCEDURE == 0 ) THEN
             WRITE(Message,*) 'Value type for property [', TRIM(Name), &
                 '] not used consistently.'
             CALL Fatal( 'ListGetRealAtIp', Message )
             RETURN
           END IF
           
           CALL ListPushActiveName(name)
           DO i=1,n
             F(i) = ptr % Coeff * &
                 ExecConstRealFunction( ptr % PROCEDURE,CurrentModel, &
                 CurrentModel % Mesh % Nodes % x( NodeIndexes(i) ), &
                 CurrentModel % Mesh % Nodes % y( NodeIndexes(i) ), &
                 CurrentModel % Mesh % Nodes % z( NodeIndexes(i) ) )
           END DO
           CALL ListPopActiveName()
           
         END SELECT
       END IF
       
       IF( Handle % ConstantInList ) THEN
         val = F(1)
       ELSE
         IF(.NOT. PRESENT(Basis)) THEN
           CALL Fatal('ListGetRealAtIp','Parameter > Basis < is required!')
         ELSE
           val = SUM( Basis(1:n) * F(1:n) )
         END IF
       END IF
     END IF


     IF ( PRESENT( minv ) ) THEN
       IF ( val < minv ) THEN
         WRITE( Message,*) 'Given value ',val, ' for property: ', '[', TRIM(Name),']', &
             ' smaller than given minimum: ', minv
         CALL Fatal( 'ListGetRealAtIp', Message )
       END IF
     END IF
       
     IF ( PRESENT( maxv ) ) THEN
       IF ( val > maxv ) THEN
         WRITE( Message,*) 'Given value ',val, ' for property: ', '[', TRIM(Name),']', &
             ' larger than given maximum ', maxv
         CALL Fatal( 'ListGetRealAtIp', Message )
       END IF
     END IF


   END FUNCTION ListGetRealAtIp
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Initializes the handle to save just a little bit for constant valued.
!> This is not mandatory but may still be used. 
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE ListInitRealAtIp( Handle,Section,Name,minv,maxv ) 
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     CHARACTER(LEN=*)  :: Section,Name
     REAL(KIND=dp), OPTIONAL :: minv,maxv
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     TYPE(ValueListEntry_t), POINTER :: ptr
     INTEGER :: i, n, NoVal
     TYPE(Model_t), POINTER :: Model
     REAL(KIND=dp)  :: val
     LOGICAL :: ConstantEverywhere, Found
     REAL(KIND=dp), POINTER :: Basis(:)
     INTEGER, POINTER :: NodeIndexes(:)
     TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------

     NoVal = 0
     ConstantEverywhere = .TRUE.
     Model => CurrentModel
     Handle % Name = Name 
     Handle % Initialized = .TRUE.

     ! Allocate space for temporal values
     IF(.NOT. Handle % AllocationsDone ) THEN
       n = CurrentModel % Mesh % MaxElementNodes
       ALLOCATE( Handle % Values(n) )
       Handle % Values = 0.0_dp
       Handle % AllocationsDone = .TRUE.
     END IF

     ! Loop over the full lists and check whether the keyword is constant in all of them
     ! Unfortunately only keywords that are present in all lists can be deduced to be 
     ! truly constant. This is a conservative estimate, but best we can do without 
     ! involving the active elements into the picture.
     i = 0
     Element => CurrentModel % Mesh % Elements(1)
     CurrentModel % CurrentElement => CurrentModel % Mesh % Elements(1)
     n = Element % TYPE % NumberOfNodes 
     ALLOCATE( Basis(n) )
     Basis = 0.0_dp

     DO WHILE(.TRUE.) 
       i = i + 1

       SELECT CASE ( Section ) 
       CASE('Material')
         IF(i > Model % NumberOfMaterials) EXIT
         List => Model % Materials(i) % Values

       CASE('Body Force')
         IF(i > Model % NumberOfBodyForces) EXIT        
         List => Model % BodyForces(i) % Values
         
       CASE('Initial Condition')
         IF( i > Model % NumberOfICs) EXIT
         List => Model % ICs(i) % Values

       CASE('Boundary Condition')
         IF( i > Model % NumberOfBCs ) EXIT        
         List => Model % BCs(i) % Values

       CASE DEFAULT
         CALL Fatal('ListInitRealAtIP','Unknown section: '//TRIM(Section))
         
       END SELECT
       
       ! If the parameter is not defined in some list we cannot really be sure
       ! that it is intentionally used as a zero. Hence we cannot assume that the
       ! keyword is constant. 
       ptr => ListFind(List,Name,Found)
       IF ( .NOT.ASSOCIATED(ptr) ) THEN
         ConstantEverywhere = .FALSE.
         EXIT
       END IF

       ! The value must be constant in each list and
       Handle % ConstantInList = .FALSE.
       val = ListGetRealAtIp( Handle,List,Basis,Name,minv=minv,maxv=maxv)
       IF( .NOT. Handle % ConstantInList ) THEN
         ConstantEverywhere = .FALSE.
         EXIT
       END IF
       
       ! and each list must have the same constant value
       NoVal = NoVal + 1
       IF( NoVal == 1 ) THEN
         Handle % ConstantValue = val
       ELSE IF( ABS ( Handle % ConstantValue - val ) > TINY( val ) ) THEN
         ConstantEverywhere = .FALSE.
         EXIT
       END IF
     END DO
     
     IF( ConstantEverywhere ) THEN       
       WRITE( Message,'(A)') 'Constant keyword: '//TRIM(Name)
       CALL Info('ListInitRealAtIp',Message,Level=8)
       WRITE( Message,'(A,ES15.4)') 'Constant value:',val
       CALL Info('ListInitRealAtIp',Message,Level=8)
     END IF

     Handle % ConstantEverywhere = ConstantEverywhere 

   END SUBROUTINE ListInitRealAtIp



!------------------------------------------------------------------------------
!> Gets a constant real array from the list by its name.
!------------------------------------------------------------------------------
   RECURSIVE FUNCTION ListGetConstRealArray( List,Name,Found,UnfoundFatal ) RESULT( F )
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found, UnfoundFatal
!------------------------------------------------------------------------------
     REAL(KIND=dp), POINTER  :: F(:,:)
     INTEGER :: i,j,N1,N2
     TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------
     NULLIFY( F ) 
     ptr => ListFind(List,Name,Found)
     IF (.NOT.ASSOCIATED(ptr) ) THEN
       IF(PRESENT(UnfoundFatal)) THEN
         IF(UnfoundFatal) THEN
           WRITE(Message, '(A,A)') "Failed to find ConstRealArray: ",Name
           CALL Fatal("ListGetInteger", Message)
         END IF
       END IF
       RETURN
     END IF

     IF ( .NOT. ASSOCIATED(ptr % FValues) ) THEN
       WRITE(Message,*) 'Value type for property [', TRIM(Name), &
               '] not used consistently.'
       CALL Fatal( 'ListGetConstRealArray', Message )
       RETURN
     END IF

     N1 = SIZE( ptr % FValues,1 )
     N2 = SIZE( ptr % FValues,2 )

     F => ptr % FValues(:,:,1)

     IF ( ptr % PROCEDURE /= 0 ) THEN
       CALL ListPushActiveName(name)
       DO i=1,N1
         DO j=1,N2
           F(i,j) = ExecConstRealFunction( ptr % PROCEDURE,CurrentModel,0.0d0,0.0d0,0.0d0 )
         END DO
       END DO
       CALL ListPopActiveName()
     END IF
   END FUNCTION ListGetConstRealArray
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Gets a real array from the list by its name,
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE ListGetRealArray( List,Name,F,N,NodeIndexes,Found )
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found
     INTEGER :: N,NodeIndexes(:)
     REAL(KIND=dp), POINTER :: F(:,:,:), G(:,:)
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr

     TYPE(Variable_t), POINTER :: Variable, CVar, TVar

     REAL(KIND=dp) :: T(MAX_FNC)
     INTEGER :: i,j,k,nlen,N1,N2,k1,l
     CHARACTER(LEN=2048) :: tmp_str, cmd
     LOGICAL :: AllGlobal
!------------------------------------------------------------------------------
     ptr => ListFind(List,Name,Found)
     IF ( .NOT.ASSOCIATED(ptr) ) RETURN

     IF ( .NOT. ASSOCIATED(ptr % FValues) ) THEN
       CALL Fatal( 'ListGetRealArray', &
           'Value type for property > '// TRIM(Name) // '< not used consistently.')
     END IF

     N1 = SIZE(ptr % FValues,1)
     N2 = SIZE(ptr % FValues,2)

     IF ( .NOT.ASSOCIATED( F ) ) THEN
       ALLOCATE( F(N1,N2,N) )
     ELSE IF ( SIZE(F,1)/=N1.OR.SIZE(F,2)/=N2.OR.SIZE(F,3)/= N ) THEN
       DEALLOCATE( F )
       ALLOCATE( F(N1,N2,N) )
     END IF

     SELECT CASE(ptr % TYPE)
     CASE ( LIST_TYPE_CONSTANT_TENSOR )
       DO i=1,n
         F(:,:,i) = ptr % Coeff * ptr % FValues(:,:,1)
       END DO

       IF ( ptr % PROCEDURE /= 0 ) THEN
         CALL ListPushActiveName(name)
         DO i=1,N1
           DO j=1,N2
             F(i,j,1) = ptr % Coeff * &
                 ExecConstRealFunction( ptr % PROCEDURE, &
                 CurrentModel, 0.0_dp, 0.0_dp, 0.0_dp )
           END DO
         END DO
         CALL ListPopActiveName()
       END IF
   
     
     CASE( LIST_TYPE_VARIABLE_TENSOR,LIST_TYPE_VARIABLE_TENSOR_STR )
       TVar => VariableGet( CurrentModel % Variables, 'Time' ) 
       WRITE( cmd, '(a,e15.8)' ) 'tx=0; st = ', TVar % Values(1)
       k = LEN_TRIM(cmd)
       CALL matc( cmd, tmp_str, k )

       CALL ListPushActiveName(name)
       DO i=1,n
         k = NodeIndexes(i)
         CALL ListParseStrToValues( Ptr % DependName, Ptr % DepNameLen, k, Name, T, j, AllGlobal)
         IF ( ANY(T(1:j)==HUGE(1._dP)) ) CYCLE

         IF ( ptr % TYPE==LIST_TYPE_VARIABLE_TENSOR_STR) THEN
           DO l=1,j
             WRITE( cmd, '(a,g19.12)' ) 'tx('//TRIM(i2s(l-1))//')=', T(l)
             k1 = LEN_TRIM(cmd)
             CALL matc( cmd, tmp_str, k1 )
           END DO

           cmd = ptr % CValue
           k1 = LEN_TRIM(cmd)
           CALL matc( cmd, tmp_str, k1 )
           READ( tmp_str(1:k1), * ) ((F(j,k,i),k=1,N2),j=1,N1)
         ELSE IF ( ptr % PROCEDURE /= 0 ) THEN
           G => F(:,:,i)
           CALL ExecRealArrayFunction( ptr % PROCEDURE, CurrentModel, &
                     NodeIndexes(i), T, G )
         ELSE
           DO j=1,N1
             DO k=1,N2
               F(j,k,i) = InterpolateCurve(ptr % TValues, ptr % FValues(j,k,:), &
                                T(1), ptr % CubicCoeff )
             END DO
           END DO
         END IF
         IF( AllGlobal ) EXIT
       END DO
       CALL ListPopActiveName()

       IF( AllGlobal ) THEN
         DO i=2,n
           DO j=1,N1
             DO k=1,N2
               F(j,k,i) = F(j,k,1) 
             END DO
           END DO
         END DO
       END IF

       IF( ABS( ptr % Coeff - 1.0_dp ) > EPSILON( ptr % Coeff ) ) THEN
         F = ptr % Coeff * F
       END IF
  
     CASE DEFAULT
       F = 0.0d0
       DO i=1,N1
         IF ( PRESENT( Found ) ) THEN
           F(i,1,:) = ListGetReal( List,Name,N,NodeIndexes,Found )
         ELSE
           F(i,1,:) = ListGetReal( List,Name,N,NodeIndexes )
         END IF
       END DO
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE ListGetRealArray
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Gets a real vector from the list by its name
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE ListGetRealVector( List,Name,F,N,NodeIndexes,Found )
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found
     INTEGER :: N,NodeIndexes(:)
     REAL(KIND=dp), TARGET :: F(:,:)
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr

     TYPE(Variable_t), POINTER :: Variable, CVar, TVar

     REAL(KIND=dp), ALLOCATABLE :: G(:,:)
     REAL(KIND=dp) :: T(MAX_FNC)
     REAL(KIND=dp), POINTER :: RotMatrix(:,:)
     INTEGER :: i,j,k,nlen,N1,N2,k1,S1,S2,l, cnt
     CHARACTER(LEN=2048) :: tmp_str, cmd
     LOGICAL :: AllGlobal, lFound, AnyFound
!------------------------------------------------------------------------------
     ptr => ListFind(List,Name,lFound)
     IF ( .NOT.ASSOCIATED(ptr) ) THEN
       IF(PRESENT(Found)) Found = .FALSE.
       AnyFound = .FALSE.
       DO i=1,SIZE(F,1)
         F(i,1:n) = ListGetReal(List,TRIM(Name)//' '//TRIM(I2S(i)),n,NodeIndexes,lFound)
         AnyFound = AnyFound.OR.lFound
       END DO
       IF(PRESENT(Found)) THEN
          Found = AnyFound
       ELSE IF(.NOT.AnyFound) THEN
          CALL Warn( 'ListFind', 'Requested property ['//TRIM(Name)//'] not found')
       END IF
       IF( .NOT. AnyFound ) RETURN
       GOTO 200
     END IF

     F = 0._dp
     cnt = 0
     ALLOCATE(G(SIZE(F,1),SIZE(F,2)))

100  CONTINUE

     IF ( .NOT. ASSOCIATED(ptr % FValues) ) THEN
       CALL Fatal( 'ListGetRealVector', &
           'Value type for property > '// TRIM(Name) // '< not used consistently.')
     END IF

     N1 = SIZE(ptr % FValues,1)

     SELECT CASE(ptr % TYPE)
     CASE ( LIST_TYPE_CONSTANT_TENSOR )
       DO i=1,n
         G(:,i) = ptr % Coeff * ptr % FValues(:,1,1)
       END DO

       IF ( ptr % PROCEDURE /= 0 ) THEN
         CALL ListPushActiveName(name)
         DO i=1,n1
           F(i,1) = ptr % Coeff * &
             ExecConstRealFunction( ptr % PROCEDURE, &
               CurrentModel, 0.0_dp, 0.0_dp, 0.0_dp )
         END DO
         CALL ListPopActiveName()
       END IF
     
     CASE( LIST_TYPE_VARIABLE_TENSOR,LIST_TYPE_VARIABLE_TENSOR_STR )
       TVar => VariableGet( CurrentModel % Variables, 'Time' ) 
       WRITE( cmd, '(a,e15.8)' ) 'tx=0; st = ', TVar % Values(1)
       k = LEN_TRIM(cmd)
       CALL matc( cmd, tmp_str, k )

       CALL ListPushActiveName(name)
       DO i=1,n
         k = NodeIndexes(i)
         CALL ListParseStrToValues( Ptr % DependName, Ptr % DepNameLen, k, Name, T, j, AllGlobal)
         IF ( ANY(T(1:j)==HUGE(1._dP)) ) CYCLE

         IF ( ptr % TYPE==LIST_TYPE_VARIABLE_TENSOR_STR) THEN
           DO l=1,j
             WRITE( cmd, '(a,g19.12)' ) 'tx('//TRIM(i2s(l-1))//')=', T(l)
             k1 = LEN_TRIM(cmd)
             CALL matc( cmd, tmp_str, k1 )
           END DO

           cmd = ptr % CValue
           k1 = LEN_TRIM(cmd)
           CALL matc( cmd, tmp_str, k1 )
           READ( tmp_str(1:k1), * ) (G(j,i),j=1,N1)
         ELSE IF ( ptr % PROCEDURE /= 0 ) THEN
           CALL ExecRealVectorFunction( ptr % PROCEDURE, CurrentModel, &
                     NodeIndexes(i), T, G(:,i) )
         ELSE
           DO k=1,n1
             G(k,i) = InterpolateCurve(ptr % TValues, &
                   ptr % FValues(k,1,:), T(MIN(j,k)), ptr % CubicCoeff )
           END DO
         END IF

         IF( AllGlobal ) EXIT
       END DO
       CALL ListPopActiveName()

       IF( AllGlobal ) THEN
         DO i=2,n
           DO j=1,N1
             G(j,i) = G(j,1) 
           END DO
         END DO
       END IF

       IF( ABS( ptr % Coeff - 1.0_dp ) > EPSILON( ptr % Coeff ) ) THEN
         G = ptr % Coeff * G
       END IF
  
     CASE DEFAULT
       G = 0.0d0
       DO i=1,N1
         IF ( PRESENT( Found ) ) THEN
           G(i,:) = ListGetReal( List,Name,N,NodeIndexes,Found )
         ELSE
           G(i,:) = ListGetReal( List,Name,N,NodeIndexes )
         END IF
       END DO
     END SELECT


     F = F + G
     cnt = cnt + 1
     ptr => ListFind(List,Name//'{'//TRIM(I2S(cnt))//'}',lFound)
     IF(ASSOCIATED(ptr)) GOTO 100

200  IF( ListGetLogical( List, Name//' Property Rotate', lFound ) ) THEN
       RotMatrix => ListGetConstRealArray( List,'Property Rotation Matrix',lFound )
       IF( .NOT. ASSOCIATED( RotMatrix ) ) THEN
         CALL Fatal('ListGetRealVector','Property rotation matrix not given for: '//TRIM(Name))
       END IF
       IF( SIZE(F,1) /= 3 ) THEN
         CALL Fatal('ListGetRealVector','Property may be rotated only with three components!')
       END IF
       DO i = 1,SIZE(F,2) 
         F(1:3,i) = MATMUL( RotMatrix, F(1:3,i) )
       END DO
     END IF


!------------------------------------------------------------------------------
   END SUBROUTINE ListGetRealVector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Gets a real derivative from. This is only available for tables with dependencies.
!------------------------------------------------------------------------------
   RECURSIVE FUNCTION ListGetDerivValue(List,Name,N,NodeIndexes,dT) RESULT(F)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER ::  List
     CHARACTER(LEN=*) :: Name
     INTEGER :: N,NodeIndexes(:)
     REAL(KIND=dp), OPTIONAL :: dT
     REAL(KIND=dp) :: F(N)
!------------------------------------------------------------------------------
     TYPE(Variable_t), POINTER :: Variable
     TYPE(ValueListEntry_t), POINTER :: ptr
     INTEGER :: i,k,l
     REAL(KIND=dp) :: T,T1(1),T2(1),F1,F2
!------------------------------------------------------------------------------

     F = 0.0D0
     ptr => ListFind(List,Name)


     IF ( .NOT.ASSOCIATED(ptr) ) RETURN


     SELECT CASE(ptr % TYPE)
       CASE( LIST_TYPE_VARIABLE_SCALAR )
         
         IF ( ptr % PROCEDURE /= 0 ) THEN
           IF( .NOT. PRESENT( dT ) ) THEN
             CALL Fatal('ListGetDerivValue','Numerical derivative of function requires dT')
           END IF
           Variable => VariableGet( CurrentModel % Variables,ptr % DependName ) 
           IF( .NOT. ASSOCIATED( Variable ) ) THEN
             CALL Fatal('ListGetDeriveValue','Cannot derivate with variable: '//TRIM(ptr % DependName))
           END IF

           DO i=1,n
             k = NodeIndexes(i)            
             IF ( ASSOCIATED(Variable % Perm) ) k = Variable % Perm(k)
             IF ( k > 0 ) THEN
               T = Variable % Values(k) 
               T1(1) = T + 0.5_dp * dT
               T2(1) = T - 0.5_dp * dT 
               F1 = ExecRealFunction( ptr % PROCEDURE,CurrentModel, NodeIndexes(i), T1 )
               F2 = ExecRealFunction( ptr % PROCEDURE,CurrentModel, NodeIndexes(i), T2 )
               F(i) = ptr % Coeff * ( F1 - F2 ) / dT
             END IF
           END DO

         ELSE
           IF ( .NOT. ASSOCIATED(ptr % FValues) ) THEN
             CALL Fatal( 'ListGetDerivValue', &
                 'Value type for property > '// TRIM(Name) // '< not used consistently.')
           END IF
           Variable => VariableGet( CurrentModel % Variables,ptr % DependName ) 
           IF( .NOT. ASSOCIATED( Variable ) ) THEN
             CALL Fatal('ListGetDeriveValue','Cannot derivate with variable: '//TRIM(ptr % DependName))
           END IF
           DO i=1,n
             k = NodeIndexes(i)
             IF ( ASSOCIATED(Variable % Perm) ) k = Variable % Perm(k)
             IF ( k > 0 ) THEN
               T = Variable % Values(k)
               F(i) = ptr % Coeff * &
                   DerivateCurve(ptr % TValues,ptr % FValues(1,1,:), &
                   T, ptr % CubicCoeff )
             END IF
           END DO
         END IF


       CASE DEFAULT 
         CALL Fatal( 'ListGetDerivValue', &
             'No automated derivation possible for > '//TRIM(Name)//' <' )

     END SELECT


   END FUNCTION ListGetDerivValue
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Given the body of a keyword find the 1st free keyword in the list structure.
!> The intended use for this is in Solver_init to decleare exported variables
!> without the risk of running over some existing ones. 
!------------------------------------------------------------------------------
  FUNCTION NextFreeKeyword(keyword0,List) RESULT (Keyword)

    CHARACTER(LEN=*) :: Keyword0
    TYPE(ValueList_t), POINTER  :: List
    CHARACTER(LEN=MAX_NAME_LEN) :: Keyword
    INTEGER :: No
    
    DO No = 1, 9999
      WRITE( Keyword,'(A,I0)') TRIM(Keyword0)//' ',No
      IF( .NOT. ListCheckPresent(List,Keyword)) EXIT
    END DO

!------------------------------------------------------------------------------
  END FUNCTION NextFreeKeyword
!------------------------------------------------------------------------------
 

!------------------------------------------------------------------------------
!> Check if the keyword is present in any boundary condition.
!------------------------------------------------------------------------------
   FUNCTION ListCheckPresentAnyBC( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
     INTEGER :: bc
     
     Found = .FALSE.
     DO bc = 1,Model % NumberOfBCs
       Found = ListCheckPresent( Model % BCs(bc) % Values, Name )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListCheckPresentAnyBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Check if the keyword is True in any boundary condition.
!------------------------------------------------------------------------------
   FUNCTION ListGetLogicalAnyBC( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found, GotIt
     INTEGER :: bc
     
     Found = .FALSE.
     DO bc = 1,Model % NumberOfBCs
       Found = ListgetLogical( Model % BCs(bc) % Values, Name, GotIt )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListGetLogicalAnyBC
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Check if the keyword is present in any body.
!------------------------------------------------------------------------------
   FUNCTION ListCheckPresentAnyBody( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
     INTEGER :: body
     
     Found = .FALSE.
     DO body = 1,Model % NumberOfBodies
       Found = ListCheckPresent( Model % Bodies(body) % Values, Name )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListCheckPresentAnyBody
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Check if the keyword is true in any body.
!------------------------------------------------------------------------------
   FUNCTION ListGetLogicalAnyBody( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
     INTEGER :: body
     LOGICAL :: GotIt
     
     Found = .FALSE.
     DO body = 1,Model % NumberOfBodies
       Found = ListGetLogical( Model % Bodies(body) % Values, Name, GotIt )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListGetLogicalAnyBody
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Check if the keyword is present in any body force.
!------------------------------------------------------------------------------
   FUNCTION ListCheckPresentAnyBodyForce( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
     INTEGER :: bf
     
     Found = .FALSE.
     DO bf = 1,Model % NumberOfBodyForces
       Found = ListCheckPresent( Model % BodyForces(bf) % Values, Name )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListCheckPresentAnyBodyForce
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Check if the keyword is True in any body force.
!------------------------------------------------------------------------------
   FUNCTION ListGetLogicalAnyBodyForce( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found, GotIt
     INTEGER :: bf
     
     Found = .FALSE.
     DO bf = 1,Model % NumberOfBodyForces
       Found = ListGetLogical( Model % BodyForces(bf) % Values, Name, GotIt )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListGetLogicalAnyBodyForce
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Check if the keyword is present in any material.
!------------------------------------------------------------------------------
   FUNCTION ListCheckPresentAnyMaterial( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
     INTEGER :: mat
     
     Found = .FALSE.
     DO mat = 1,Model % NumberOfMaterials
       Found = ListCheckPresent( Model % Materials(mat) % Values, Name )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListCheckPresentAnyMaterial
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Check if the keyword in any material is defined as an array
!------------------------------------------------------------------------------
   FUNCTION ListCheckAnyMaterialIsArray( Model, Name ) RESULT(IsArray)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: IsArray
     LOGICAL :: Found
     INTEGER :: mat, n1, n2
     TYPE(ValueListEntry_t), POINTER :: ptr
    
     IsArray = .FALSE.
     DO mat = 1,Model % NumberOfMaterials
       ptr => ListFind(Model % Materials(mat) % Values,Name,Found)
       IF( .NOT. ASSOCIATED( ptr ) ) CYCLE
       IF ( .NOT. ASSOCIATED(ptr % FValues) ) THEN
         WRITE(Message,*) 'Value type for property [', TRIM(Name), &
             '] not used consistently.'
         CALL Fatal( 'ListCheckAnyMaterialArray', Message )
       END IF
       n1 = SIZE( ptr % FValues,1 )
       n2 = SIZE( ptr % FValues,2 )
       IsArray =  ( n1 > 1 ) .OR. ( n2 > 1 ) 
       IF( IsArray ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListCheckAnyMaterialIsArray
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Check if the keyword is True in any material.
!------------------------------------------------------------------------------
   FUNCTION ListGetLogicalAnyMaterial( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found, GotIt
     INTEGER :: mat
     
     Found = .FALSE.
     DO mat = 1,Model % NumberOfMaterials
       Found = ListGetLogical( Model % Materials(mat) % Values, Name, GotIt )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListGetLogicalAnyMaterial
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Check if the keyword is present in any equation.
!------------------------------------------------------------------------------
   FUNCTION ListCheckPresentAnyEquation( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
     INTEGER :: eq
     
     Found = .FALSE.
     DO eq = 1,Model % NumberOfEquations
       Found = ListCheckPresent( Model % Equations(eq) % Values, Name )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListCheckPresentAnyEquation
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Check if the keyword is True in any equation.
!------------------------------------------------------------------------------
   FUNCTION ListGetLogicalAnyEquation( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found, GotIt
     INTEGER :: eq
     
     Found = .FALSE.
     DO eq = 1,Model % NumberOfEquations
       Found = ListGetLogical( Model % Equations(eq) % Values, Name, GotIt )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListGetLogicalAnyEquation
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Elmer may include scalar and vector variables which may be known by their
!> original name or have an alias. For historical reasons they are introduced
!> by two quite separate ways. This subroutine tries to make the definition of
!> variables for saving more straight-forward.
!------------------------------------------------------------------------------
  SUBROUTINE CreateListForSaving( Model, List, ShowVariables, ClearList )
!------------------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(ValueList_t), POINTER  :: List
    LOGICAL :: ShowVariables
    LOGICAL, OPTIONAL :: ClearList
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,LoopDim, VarDim,FullDim,DOFs,dim,Comp
    TYPE(Variable_t), POINTER :: Variables, Var, Var1
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName, VarStr, VarStrComp, VarStrExt, str
    LOGICAL :: IsVector, Set, GotIt, ComponentVector, ThisOnly, IsIndex, &
        EnforceVectors
    INTEGER :: Nvector, Nscalar
    TYPE(ValueList_t), POINTER :: Params

    Params => Model % Solver % Values
    Variables => Model % Mesh % Variables

    IF( .NOT. ASSOCIATED( Variables ) ) THEN
      CALL Warn('CreateListForSaving','Mesh does not include any variables!')
      RETURN
    END IF

!------------------------------------------------------------------------------
! Sometimes the list must be cleared in order to use it for a different mesh
!-----------------------------------------------------------------------------
    IF( PRESENT( ClearList ) ) THEN
      IF( ClearList ) THEN
        DO i=1,999
          WRITE(VarStr,'(A,I0)') 'Scalar Field ',i
          IF( ListCheckPresent( List, VarStr ) ) THEN
            CALL ListRemove( List, VarStr )
          ELSE
            EXIT
          END IF
        END DO

        DO i=1,999
          WRITE(VarStr,'(A,I0)') 'Vector Field ',i
          IF( ListCheckPresent( List, VarStr ) ) THEN
            CALL ListRemove( List, VarStr )
          ELSE
            EXIT
          END IF

          WRITE(VarStr,'(A,I0,A)') 'Vector Field ',i,' Component'
          IF( ListCheckPresent( List, VarStr ) ) THEN
            CALL ListRemove( List, VarStr )
          END IF
        END DO
      END IF
    END IF
    
    !-------------------------------------------------------------------
    ! First check that there is a need to create the list i.e. it is not
    ! already manually defined
    !-------------------------------------------------------------------
    IF( ListCheckPresent( List,'Scalar Field 1' ) ) THEN
      CALL Info('CreateListForSaving','Scalar Field 1 exists, creating no list!',Level=10)
      RETURN
    END IF

    IF( ListCheckPresent( List,'Vector Field 1' ) ) THEN
      CALL Info('CreateListForSaving','Vector Field 1 exists, creating no list!',Level=10)
      RETURN
    END IF

    Nscalar = 0
    Nvector = 0


    ThisOnly = .NOT. ListGetLogical( Params,'Interpolate Fields',GotIt)
    dim = Model % Mesh % MeshDim

    EnforceVectors = ListGetLogical( Params,'Enforce Vectors',GotIt)
    IF(.NOT. GotIt ) EnforceVectors = .TRUE.

    
    Var => Variables


    DO WHILE( ASSOCIATED( Var ) )

      ! Skip if variable is not active for saving       
      IF ( .NOT. Var % Output ) THEN
        Var => Var % Next
        CYCLE
      END IF
      
      ! Skip if variable is global one
      IF ( SIZE( Var % Values ) == Var % DOFs ) THEN
        Var => Var % Next
        CYCLE
      END IF

      ! Skip if variable is otherwise strange in size
      IF(.NOT. ASSOCIATED( Var % Perm ) ) THEN
        IF( Var % TYPE == Variable_on_nodes_on_elements ) THEN
          IF( SIZE( Var % Values ) /= Var % Dofs * Model % Mesh % NumberOfNodes ) THEN
            Var => Var % Next
            CYCLE
          END IF
        ELSE
          IF( SIZE( Var % Values ) /= Var % Dofs * Model % Mesh % NumberOfBulkElements ) THEN
            Var => Var % Next
            CYCLE
          END IF         
        END IF
      END IF


      VarDim = Var % Dofs
      IsVector = (VarDim > 1)
      Set = .FALSE.

      WRITE(VarName,'(A)') TRIM(Var % Name)

      SELECT CASE(Var % Name)

      CASE( 'coordinate 1','coordinate 2','coordinate 3' )
        ! These are treated separatetely as coordinates are not typically saved


      CASE( 'mesh update' )
        ! Mesh update is treated separately because its special connection to displacement

        Var1 => Variables
        DO WHILE( ASSOCIATED( Var1 ) )
          IF ( TRIM(Var1 % Name) == 'displacement' ) EXIT
          Var1 => Var1 % Next
        END DO
        IF ( .NOT. ASSOCIATED( Var1 ) ) THEN
          Set = .TRUE.
        END IF
        
      CASE('mesh update 1','mesh update 2', 'mesh update 3' )
        
      CASE( 'displacement' )
        Set = .TRUE.
        ! mesh update is by default the complement to displacement 
        Var1 => Variables
        DO WHILE( ASSOCIATED( Var1 ) )
          IF ( TRIM(Var1 % Name) == 'mesh update' ) EXIT
          Var1 => Var1 % Next
        END DO
        IF ( ASSOCIATED( Var1 ) ) THEN
          WRITE(VarStrComp,'(A,I0,A)') 'Vector Field ',Nvector+1,' Complement'
          CALL ListAddString( List ,TRIM(VarStrComp),'mesh update')
        END IF
        
      CASE( 'displacement 1','displacement 2','displacement 3')
        

      CASE DEFAULT
  
        ! All vector variables are assumed to be saved using its components
        ! rather than vector itself.

        IF ( VarDim == 1 ) THEN
          Set = .TRUE. 
          
	  str =  ' '
          j = LEN_TRIM(Var % Name)
          DO i=1,j
            str(i:i) = Var % Name(i:i)
          END DO

          IsIndex = .FALSE.
          Comp = 0
          k = INDEX( str(:j),' ',BACK=.TRUE.)

          IF( k > 0 ) THEN
            IsIndex = ( VERIFY( str(k:j),' 0123456789') == 0 )
            IF( IsIndex ) READ( str(k:j), * ) Comp
          END IF

          ! This is the easy way of checking that the component belongs to a vector
          GotIt = .FALSE.
          IF( IsIndex ) THEN
            Var1 => VariableGet(Variables,TRIM(str(1:k)))
            IF( ASSOCIATED( Var1 ) ) THEN
              GotIt = .TRUE.
              IsVector = ( Var1 % Dofs == Dim .OR. Var1 % Dofs == 3 ) 
              Set = ( Comp == 1 .OR. .NOT. IsVector )
            END IF
          END IF
                
          ! This is a hard way of ensuring that the component belongs to a vector    
          ! Check that there are exactly dim number of components
          ! If so save the quantity as a vector, otherwise componentwise
          IF( EnforceVectors .AND. .NOT. GotIt ) THEN
            IF( Comp == 1 ) THEN
              Var1 => VariableGet(Variables,TRIM(str(1:j-2))//' '//I2S(dim),ThisOnly)		
              IF( ASSOCIATED(Var1)) THEN
                Var1 => VariableGet(Variables,TRIM(str(1:j-2))//' '//I2S(dim+1),ThisOnly)		
                IsVector = .NOT. ASSOCIATED(Var1)
              END IF
              
              ! Associated to the previous case, cycle the other components of the vector
              ! and cycle them if they are part of the vector that will be detected above.
            ELSE IF( Comp <= dim ) THEN
              Var1 => VariableGet(Variables,TRIM(str(1:j-2))//' 1',ThisOnly)		
              IF( ASSOCIATED( Var1 ) ) THEN
                Var1 => VariableGet(Variables,TRIM(str(1:j-2))//' '//I2S(dim+1),ThisOnly)		
                Set = ASSOCIATED( Var1 )
              END IF
            END IF
          END IF

          ! Remove the trailing numbers as they are not needed in this case.
          IF(IsVector) WRITE(VarName,'(A)') TRIM(str(1:j-2))
        END IF
      END SELECT

      
      
      !---------------------------------------------------------------------------
      ! Set the default variable names that have not been set
      !------------------------------------------------------------------------
      IF( Set ) THEN
        IF( IsVector ) THEN          
          Nvector = Nvector + 1
          WRITE(VarStr,'(A,I0)') 'Vector Field ',Nvector
        ELSE
          Nscalar = Nscalar + 1
          WRITE(VarStr,'(A,I0)') 'Scalar Field ',Nscalar
        END IF
        CALL ListAddString( List,TRIM(VarStr),TRIM(VarName) )
      END IF

      Var => Var % Next
    END DO


    IF( ShowVariables ) THEN
      CALL Info('CreateListForSaving','Field Variables for Saving')
      DO i=1,Nscalar
        WRITE(VarStr,'(A,I0)') 'Scalar Field ',i
        VarName = ListGetString( List, VarStr,GotIt )
        IF( GotIt ) THEN
          WRITE( Message,'(A)') TRIM(VarStr)//': '//TRIM(VarName)
          CALL Info('CreateListForSaving',Message)
        END IF
      END DO

      DO i=1,Nvector
        WRITE(VarStr,'(A,I0)') 'Vector Field ',i
        VarName = ListGetString( List, VarStr,GotIt )
        IF( GotIt ) THEN
          WRITE( Message,'(A)') TRIM(VarStr)//': '//TRIM(VarName)
          CALL Info('CreateListForSaving',Message)
        END IF
      END DO

      DO i=1,Nvector
        WRITE(VarStr,'(A,I0,A)') 'Vector Field ',i,' Complement'
        VarName = ListGetString( List, VarStr, GotIt )
        IF( GotIt ) THEN
          WRITE( Message,'(A)') TRIM(VarStr)//': '//TRIM(VarName)
          CALL Info('CreateListForSaving',Message)
        END IF
      END DO
    END IF

  END SUBROUTINE CreateListForSaving
  

!------------------------------------------------------------------------------
!> A timer that uses a list structure to store the times making in 
!> generally applicable without any upper limit on the number of timers.
!> This resets the timer.
!-----------------------------------------------------------------------------

  SUBROUTINE ResetTimer(TimerName)
    CHARACTER(*) :: TimerName
    REAL(KIND=dp) :: ct, rt
#ifndef USE_ISO_C_BINDINGS
    REAL(KIND=dp) :: RealTime, CPUTime
#endif
    LOGICAL :: Found,FirstTime=.TRUE.

    IF( FirstTime ) THEN
      FirstTime=.FALSE.
      TimerPassive = ListGetLogical( CurrentModel % Simulation,'Timer Passive',Found)
      TimerResults = ListGetLogical( CurrentModel % Simulation,'Timer Results',Found)      
    END IF

    IF( TimerPassive ) RETURN

    ct = CPUTime()
    rt = RealTime()

    CALL ListAddConstReal( TimerList,TRIM(TimerName)//' cpu time',ct )
    CALL ListAddConstReal( TimerList,TRIM(TimerName)//' real time',rt )

  END SUBROUTINE ResetTimer

  
!-----------------------------------------------------------------------------
!> Delete an existing timer.
!----------------------------------------------------------------------------
  SUBROUTINE DeleteTimer(TimerName) 
    CHARACTER(*) :: TimerName
    
    IF( TimerPassive ) RETURN

    CALL ListRemove( TimerList, TRIM(TimerName)//' cpu time' ) 
    CALL ListRemove( TimerList, TRIM(TimerName)//' real time' ) 

  END SUBROUTINE DeleteTimer
 
!-----------------------------------------------------------------------------
!> Check current time of the timer.
!----------------------------------------------------------------------------
  SUBROUTINE CheckTimer(TimerName, Level, Delete, Reset) 
    CHARACTER(*) :: TimerName
    INTEGER, OPTIONAL :: Level
    LOGICAL, OPTIONAL :: Reset, Delete
    
    REAL(KIND=dp) :: ct0,rt0,ct, rt
#ifndef USE_ISO_C_BINDINGS
    REAL(KIND=dp) :: RealTime, CPUTime
#endif
    LOGICAL :: Found

    IF( TimerPassive ) RETURN

    ct0 = ListGetConstReal( TimerList,TRIM(TimerName)//' cpu time',Found) 
    IF( Found ) THEN
      rt0 = ListGetConstReal( TimerList,TRIM(TimerName)//' real time')
      ct = CPUTime() - ct0
      rt = RealTime() - rt0 
      
      WRITE(Message,'(a,2f10.4,a)') 'Elapsed time (CPU,REAL): ',ct,rt,' (s)'
      CALL Info(TRIM(TimerName),Message,Level=Level)          

      IF( TimerResults ) THEN
        CALL ListAddConstReal(CurrentModel % Simulation,&
            'res: '//TRIM(TimerName)//' cpu time',ct)
        CALL ListAddConstReal(CurrentModel % Simulation,&
            'res: '//TRIM(TimerName)//' real time',rt)
      END IF
    ELSE
      CALL Warn('CheckTimer',&
          'Requesting time from non-existing timer: '//TRIM(TimerName) )
    END IF

    IF( PRESENT( Reset ) ) THEN
      IF( Reset ) THEN
        CALL ListAddConstReal( TimerList,TRIM(TimerName)//' cpu time',ct )
        CALL ListAddConstReal( TimerList,TRIM(TimerName)//' real time',rt )
      END IF
    END IF

    IF( PRESENT( Delete ) ) THEN
      IF( Delete ) CALL DeleteTimer( TimerName )
    END IF

  END SUBROUTINE CheckTimer


!> Returns the angular frequency  
  FUNCTION ListGetAngularFrequency(ValueList,Found,UElement) RESULT(w)
    REAL(KIND=dp) :: w
    TYPE(ValueList_t), OPTIONAL, POINTER :: ValueList
    LOGICAL, OPTIONAL :: Found
    LOGICAL :: GotIt
    TYPE(Element_t), POINTER :: Element, UElement
    OPTIONAL :: UElement
    INTEGER :: elem_id,eq_id,mat_id

    ! This is rather complicated since it should replace all the various strategies
    ! that have been used in different solvers.
    !------------------------------------------------------------------------------

    ! The only way frequency may depend on element is that it sits in equation block
    !--------------------------------------------------------------------------------
    IF( PRESENT( ValueList ) ) THEN
      w = 2 * PI * ListGetCReal( ValueList,'Frequency',GotIt)
      IF(.NOT. GotIt) w = ListGetCReal( ValueList,'Angular Frequency',GotIt)
    ELSE
      GotIt = .FALSE.
    END IF

    ! It seems that the equation section is used to allow compliance with ElmerGUI
    !------------------------------------------------------------------------------
    IF( .NOT. GotIt ) THEN
      IF(PRESENT(UElement)) THEN
        Element => UElement
        eq_id = ListGetInteger( CurrentModel % Bodies(Element % bodyid) % Values,'Equation')
        w = 2 * PI * ListGetCReal( &
            CurrentModel % Equations(eq_id) % Values,'Frequency',GotIt)
        IF(.NOT. GotIt) w = ListGetCReal( &
            CurrentModel % Equations(eq_id) % Values,'Angular Frequency',GotIt)
      END IF
    END IF

    ! Check also the material section...
    !------------------------------------------------------------------------------
    IF( .NOT. GotIt ) THEN
      IF(PRESENT(UElement)) THEN
        Element => UElement
        mat_id = ListGetInteger( CurrentModel % Bodies(Element % bodyid) % Values,'Material')
        w = 2 * PI * ListGetCReal( &
          CurrentModel % Materials(mat_id) % Values,'Frequency',GotIt)
        IF(.NOT. GotIt) w = ListGetCReal( &
            CurrentModel % Materials(mat_id) % Values,'Angular Frequency',GotIt)
      END IF
    END IF

    ! Normally the constant frequency is given in Simulation (or solver) block
    !-------------------------------------------------------------------------
    IF(.NOT. GotIt) w = 2 * PI * ListGetCReal( &
        CurrentModel % Simulation,'Frequency',GotIt)
    IF(.NOT. GotIt ) w = ListGetCReal( &
        CurrentModel % Simulation,'Angular Frequency',GotIt)

    IF(.NOT. GotIt ) w = 2 * PI * ListGetCReal( &
        CurrentModel % Solver % Values,'Frequency',GotIt)
    IF(.NOT. GotIt ) w = ListGetCReal( &
        CurrentModel % Solver % Values,'Angular Frequency',GotIt)

    ! It seems that the equation section is used to allow compliance with ElmerGUI
    !------------------------------------------------------------------------------
    IF( .NOT. GotIt ) THEN
       elem_id = CurrentModel % Solver % ActiveElements(1)
       Element => CurrentModel % Elements(elem_id)
       eq_id = ListGetInteger( CurrentModel % Bodies(Element % bodyid) % Values,'Equation')
       w = 2 * PI * ListGetCReal( &
           CurrentModel % Equations(eq_id) % Values,'Frequency',GotIt)
       IF(.NOT. GotIt) w = ListGetCReal( &
           CurrentModel % Equations(eq_id) % Values,'Angular Frequency',GotIt)
    END IF

    ! Check also the material section...
    !------------------------------------------------------------------------------
    IF( .NOT. GotIt ) THEN
      elem_id = CurrentModel % Solver % ActiveElements(1)
      Element => CurrentModel % Elements(elem_id)
      mat_id = ListGetInteger( CurrentModel % Bodies(Element % bodyid) % Values,'Material')
      w = 2 * PI * ListGetCReal( &
        CurrentModel % Materials(mat_id) % Values,'Frequency',GotIt)
      IF(.NOT. GotIt) w = ListGetCReal( &
          CurrentModel % Materials(mat_id) % Values,'Angular Frequency',GotIt)
    END IF

    IF( PRESENT( Found ) ) THEN
      Found = GotIt
    ELSE IF(.NOT. GotIt ) THEN
      CALL Warn('ListGetAngularFrequency','Angular frequency could not be determined!')
    END IF
  END FUNCTION ListGetAngularFrequency


  !------------------------------------------------------------------------------
!> Returns handle to the Solver value list of the active solver
  FUNCTION ListGetSolverParams(Solver) RESULT(SolverParam)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: SolverParam
     TYPE(Solver_t), OPTIONAL :: Solver

     IF ( PRESENT(Solver) ) THEN
       SolverParam => Solver % Values
     ELSE
       SolverParam => CurrentModel % Solver % Values
     END IF
!------------------------------------------------------------------------------
   END FUNCTION ListGetSolverParams
!------------------------------------------------------------------------------



END MODULE Lists

!> \} ElmerLib
