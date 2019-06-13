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
#include "../config.h"

MODULE Lists

   USE Messages
   USE GeneralUtils
#ifdef USE_ISO_C_BINDINGS
   USE LoadMod
#endif
   
   IMPLICIT NONE

   INTEGER, PARAMETER :: LIST_TYPE_LOGICAL = 1
   INTEGER, PARAMETER :: LIST_TYPE_STRING  = 2
   INTEGER, PARAMETER :: LIST_TYPE_INTEGER = 3
   INTEGER, PARAMETER :: LIST_TYPE_CONSTANT_SCALAR = 4
   INTEGER, PARAMETER :: LIST_TYPE_VARIABLE_SCALAR = 5
   INTEGER, PARAMETER :: LIST_TYPE_CONSTANT_SCALAR_STR = 6
   INTEGER, PARAMETER :: LIST_TYPE_VARIABLE_SCALAR_STR = 7
   INTEGER, PARAMETER :: LIST_TYPE_CONSTANT_SCALAR_PROC = 8
   INTEGER, PARAMETER :: LIST_TYPE_CONSTANT_TENSOR = 9
   INTEGER, PARAMETER :: LIST_TYPE_VARIABLE_TENSOR = 10
   INTEGER, PARAMETER :: LIST_TYPE_CONSTANT_TENSOR_STR = 11
   INTEGER, PARAMETER :: LIST_TYPE_VARIABLE_TENSOR_STR = 12

   INTEGER, PARAMETER :: SECTION_TYPE_BODY = 1
   INTEGER, PARAMETER :: SECTION_TYPE_MATERIAL = 2
   INTEGER, PARAMETER :: SECTION_TYPE_BF = 3
   INTEGER, PARAMETER :: SECTION_TYPE_IC = 4
   INTEGER, PARAMETER :: SECTION_TYPE_BC = 5
   INTEGER, PARAMETER :: SECTION_TYPE_COMPONENT = 6
   INTEGER, PARAMETER :: SECTION_TYPE_SIMULATION = 7
   INTEGER, PARAMETER :: SECTION_TYPE_CONSTANTS = 8
   INTEGER, PARAMETER :: SECTION_TYPE_EQUATION = 9
   

   INTEGER, PARAMETER :: MAX_FNC = 32
   
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

#ifdef HAVE_LUA
   interface ElmerEvalLua
     module procedure ElmerEvalLuaS, ElmerEvalLuaT, ElmerEvalLuaV
   end INTERFACE
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
   LOGICAL, SAVE, PRIVATE :: TimerPassive, TimerCumulative, TimerRealTime, TimerCPUTime
   CHARACTER(LEN=MAX_NAME_LEN), SAVE, PRIVATE :: TimerPrefix
   
   
   LOGICAL, PRIVATE :: DoNamespaceCheck = .FALSE.

CONTAINS

!> Tag the active degrees of freedom and number them in order of appearance. 
!------------------------------------------------------------------------------
  FUNCTION InitialPermutation( Perm,Model,Solver,Mesh, &
                   Equation,DGSolver,GlobalBubbles ) RESULT(k)
!------------------------------------------------------------------------------
     USE PElementMaps
     TYPE(Model_t)  :: Model
     TYPE(Mesh_t)   :: Mesh
     TYPE(Solver_t), TARGET :: Solver
     INTEGER :: Perm(:)
     CHARACTER(LEN=*) :: Equation
     LOGICAL, OPTIONAL :: DGSolver, GlobalBubbles
!------------------------------------------------------------------------------
     INTEGER i,j,l,t,n,e,k,k1,EDOFs, FDOFs, BDOFs, ndofs, el_id
     INTEGER :: Indexes(128)
     INTEGER, POINTER :: Def_Dofs(:)
     INTEGER, ALLOCATABLE :: EdgeDOFs(:), FaceDOFs(:)
     LOGICAL :: FoundDG, DG, DB, GB, Found, Radiation
     TYPE(Element_t),POINTER :: Element, Edge, Face
     CHARACTER(*), PARAMETER :: Caller = 'InitialPermutation'
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

     IF( DG ) THEN    
       DB = ListGetLogical( Solver % Values,'DG Reduced Basis',Found ) 
     ELSE
       DB = .FALSE.
     END IF
       
     ! Discontinuous bodies need special body-wise numbering
     IF ( DB ) THEN
       BLOCK
         INTEGER, ALLOCATABLE :: NodeIndex(:)
         INTEGER :: body_id, MaxGroup, group0, group
         INTEGER, POINTER :: DgMap(:), DgMaster(:), DgSlave(:)
         LOGICAL :: GotDgMap, GotMaster, GotSlave
         
         DgMap => ListGetIntegerArray( Solver % Values,'DG Reduced Basis Mapping',GotDgMap )
         DgMaster => ListGetIntegerArray( Solver % Values,'DG Reduced Basis Master Bodies',GotMaster )
         DgSlave => ListGetIntegerArray( Solver % Values,'DG Reduced Basis Slave Bodies',GotSlave )
                  
         IF( GotDgMap ) THEN
           IF( SIZE( DgMap ) /= Model % NumberOfBodies ) THEN
             CALL Fatal('InitialPermutation','Invalid size of > Dg Reduced Basis Mapping <')
           END IF
           MaxGroup = MAXVAL( DgMap )
         ELSE IF( GotMaster ) THEN
           MaxGroup = 2
         ELSE
           MaxGroup = Model % NumberOfBodies
         END IF
         
         ALLOCATE( NodeIndex( Mesh % NumberOfNodes ) )
         
         DO group0 = 1, MaxGroup

           ! If we have master-slave lists then nullify the slave nodes at the master
           ! interface since we want new indexes here. 
           IF( GotSlave .AND. group0 == 2 ) THEN             
             DO t=1,Mesh % NumberOfBulkElements
               Element => Mesh % Elements(t)                
               group = Element % BodyId               
               IF( ANY( DgSlave == group ) ) THEN                 
                 NodeIndex( Element % NodeIndexes ) = 0
               END IF
             END DO
           ELSE
             ! In generic case nullify all indexes already set            
             NodeIndex = 0
           END IF
           
           k1 = k

           CALL Info('InitialPermutation',&
               'Group '//TRIM(I2S(group0))//' starts from index '//TRIM(I2S(k1)),Level=10)
           
           DO t=1,Mesh % NumberOfBulkElements
             Element => Mesh % Elements(t) 
             
             group = Element % BodyId
             
             IF( GotMaster ) THEN
               IF( group0 == 1 ) THEN
                 ! First loop number dofs in "master bodies" only
                 IF( .NOT. ANY( DgMaster == group ) ) CYCLE
               ELSE 
                 ! Second loop number dofs in all bodies except "master bodies"
                 IF( ANY( DgMaster == group ) ) CYCLE
               END IF
             ELSE IF( GotDgMap ) THEN
               group = DgMap( group ) 
               IF( group0 /= group ) CYCLE
             ELSE
               IF( group0 /= group ) CYCLE
             END IF
               
             IF ( CheckElementEquation(Model,Element,Equation) ) THEN
               FoundDG = FoundDG .OR. Element % DGDOFs > 0
               DO i=1,Element % DGDOFs
                 j = Element % NodeIndexes(i)
                 IF( NodeIndex(j) == 0 ) THEN
                   k = k + 1
                   NodeIndex(j) = k
                 END IF
                 Perm( Element % DGIndexes(i) ) = NodeIndex(j)
               END DO
             END IF
           END DO

           IF( k > k1 ) THEN
             CALL Info( Caller,'Group '//TRIM(I2S(group0))//&
                 ' has '//TRIM(I2S(k-k1))//' db dofs',Level=15)
           END IF
         END DO

         CALL Info(Caller,'Numbered '//TRIM(I2S(k))//&
             ' db nodes from bulk hits',Level=15)

         IF ( FoundDG ) THEN
           RETURN ! Discontinuous bodies !!!
         END IF
       END BLOCK
     END IF


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

       CALL Info(Caller,'Numbered '//TRIM(I2S(k))//&
           ' nodes from face hits',Level=15)
       k1 = k

       
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

       CALL Info(Caller,'Numbered '//TRIM(I2S(k-k1))//&
           ' nodes from bulk hits',Level=15)
       
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

     ! Here we create the initial permutation such that the conforming dofs are eliminated. 
     IF( ListGetLogical( Solver % Values,'Apply Conforming BCs',Found ) ) THEN
       Solver % PeriodicFlipActive = .FALSE.
       n = 0
       IF( ASSOCIATED( Mesh % PeriodicPerm ) ) THEN
         ! Set the eliminated dofs to zero and renumber
         WHERE( Mesh % PeriodicPerm > 0 ) Perm = -Perm

         k = 0                  
         DO i=1,SIZE( Perm )
           IF( Perm(i) > 0 ) THEN
             k = k + 1
             Perm(i) = k
           END IF
         END DO
                          
         DO i=1,SIZE( Mesh % PeriodicPerm )
           j = Mesh % PeriodicPerm(i)
           IF( j > 0 .AND. Perm(i) /= 0 ) THEN             
             Perm(i) = Perm(j)
             IF(Mesh % PeriodicFlip(i)) n = n + 1
           END IF
         END DO

         Solver % PeriodicFlipActive = ( n > 0 )
         CALL Info('InitialPermutation','Number of periodic flips in the field: '//TRIM(I2S(n)),Level=8)
       END IF
     END IF
    
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
      CHARACTER(LEN=MAX_NAME_LEN) :: PrevEquation
      
      LOGICAL :: Flag,Found,PrevFlag

      INTEGER :: k,body_id,prev_body_id = -1
      
      SAVE Prev_body_id, PrevEquation, PrevFlag
!$OMP THREADPRIVATE(Prev_body_id, PrevEquation, PrevFlag)

      body_id = Element % BodyId

      IF( body_id == prev_body_id) THEN
         IF (Equation == PrevEquation) THEN
            Flag = PrevFlag
            RETURN
         END IF
      END IF

      prev_body_id = body_id
      PrevEquation = Equation
              
      Flag = .FALSE.      
      IF ( body_id > 0 .AND. body_id <= Model % NumberOfBodies ) THEN
         k = ListGetInteger( Model % Bodies(body_id) % Values, 'Equation', Found, &
                 minv=1, maxv=Model % NumberOFEquations )
         IF ( k > 0 ) THEN
           Flag = ListGetLogical(Model % Equations(k) % Values,Equation,Found)
         END IF
       END IF
       PrevFlag = Flag
       
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
      TYPE(Solver_t), TARGET, OPTIONAL :: Solver
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
      TYPE(Solver_t), POINTER :: VSolver
!------------------------------------------------------------------------------

      CALL Info('VariableAdd','Adding variable > '//TRIM(Name)//&
          ' < of size '//TRIM(I2S(SIZE(Values))),Level=15)

      NULLIFY(VSolver)
      IF (PRESENT(Solver)) VSolver => Solver
      
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
      ptr % NonlinValues => NULL()
      ptr % SteadyValues => NULL()
      ptr % NonlinIter = 0

      ptr % Solver => VSolver
      ptr % PrimaryMesh => Mesh

      ptr % Valid  = .TRUE.
      ptr % Output = .TRUE.
      ptr % Secondary = .FALSE.
      ptr % ValuesChanged = .TRUE.
      
! Converged information undefined = -1, not = 0, yes = 1
      ptr % NonlinConverged = -1
      ptr % SteadyConverged = -1    

      IF ( PRESENT( Secondary ) ) THEN
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
use spariterglobals
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
       IF (ASSOCIATED(Var % Values) ) THEN
         IF( SIZE( Var % Values ) == Var % DOFs ) THEN
           Var => Var % Next
           CYCLE
         END IF 
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

       IF(SIZE(Var % Values)<=0) GotValues = .FALSE.

       IF (ASSOCIATED(Var % Perm)) THEN
         Var1 => VariableList
         DO WHILE(ASSOCIATED(Var1))
           IF (.NOT.ASSOCIATED(Var,Var1)) THEN
             IF (ASSOCIATED(Var % Perm,Var1 % Perm)) &
               Var1 % Perm => NULL()
           END IF
           Var1 => Var1 % Next
         END DO
  
         IF(SIZE(Var % Perm)>0) THEN
           DEALLOCATE( Var % Perm)
         ELSE
           GotValues = .FALSE.
         END IF
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

    !cycle other variables to check for Perm association
    IF (ASSOCIATED(RmVar % Perm)) THEN
      Var => Variables
      DO WHILE(ASSOCIATED(Var))
        IF(ASSOCIATED(RmVar, Var)) &
             CALL Fatal("VariableRemove", "Programming Error - Variable appears twice in list?")
        IF (ASSOCIATED(Var % Perm,RmVar % Perm)) THEN
          RmVar % Perm => NULL()
          EXIT
        END IF
        Var => Var % Next
      END DO

      !ASSOCIATION between zero-length arrays cannot be tested
      !so nullify it anyway, just to be safe. Technically results
      !in a memory leak (of size zero??)
      IF(SIZE(RmVar % Perm) == 0) RmVar % Perm => NULL()
    END IF



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
      Perm,Output,Secondary,VarType,Global,InitValue,IpPoints)
!------------------------------------------------------------------------------
      TYPE(Variable_t), POINTER :: Variables
      TYPE(Mesh_t),   TARGET :: Mesh
      TYPE(Solver_t), TARGET, OPTIONAL :: Solver
      CHARACTER(LEN=*) :: Name
      INTEGER, OPTIONAL :: DOFs
      REAL(KIND=dp), OPTIONAL, POINTER :: Values(:)
      LOGICAL, OPTIONAL :: Output
      INTEGER, OPTIONAL, POINTER :: Perm(:)
      LOGICAL, OPTIONAL :: Secondary
      INTEGER, OPTIONAL :: VarType
      LOGICAL, OPTIONAL :: Global
      REAL(KIND=dp), OPTIONAL :: InitValue
      LOGICAL, OPTIONAL :: IpPoints
!------------------------------------------------------------------------------
      CHARACTER(LEN=MAX_NAME_LEN) :: tmpname
      REAL(KIND=dp), POINTER :: Component(:), TmpValues(:)
      INTEGER :: i,nsize, ndofs, FieldType
      LOGICAL :: IsPerm, IsGlobal, IsIPPoints
!------------------------------------------------------------------------------
            
      IF( PRESENT( DOFs ) ) THEN
        ndofs = Dofs
      ELSE
        ndofs = 1
      END IF

      IsPerm = .FALSE.
      IsGlobal = .FALSE.
      IsIPPoints = .FALSE.

      IsPerm = PRESENT( Perm ) 
      IF( PRESENT( Global ) ) IsGlobal = Global
      IF( PRESENT( IPPoints ) ) IsIPPoints = IPPoints

      IF( PRESENT( VarType ) ) THEN
        FieldType = VarType
      ELSE
        FieldType = variable_on_nodes
      END IF

      

      CALL Info('VariableAddVector','Adding variable > '//TRIM(Name)//' < with '&
          //TRIM(I2S(ndofs))//' components',Level=15)
      
      IF(PRESENT(Values)) THEN
        TmpValues => Values
      ELSE
        IF( IsPerm ) THEN
          nsize = MAXVAL( Perm ) 
        ELSE IF( IsGlobal ) THEN
          nsize = 1 
        ELSE IF( IsIpPoints ) THEN
          IF( .NOT. PRESENT( Solver ) ) THEN
            CALL Fatal('VariableAddVector','Integration point variable needs a Solver!')
          END IF
          IF( .NOT. ASSOCIATED( Solver % IPTable ) ) THEN
            CALL Fatal('VariableAddVector','Integration point variable needs an IpTable')
          END IF
          nsize = Solver % IPTable % IPCount
        ELSE
          nsize = Mesh % NumberOfNodes          
        END IF
        CALL Info('VariableAddVector','Allocating field of size: '//TRIM(I2S(nsize)),Level=12)
        
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
              Perm,Output,Secondary,VarType)
        END DO
      END IF

      CALL VariableAdd( Variables,Mesh,Solver,Name,nDOFs,TmpValues,&
            Perm,Output,Secondary,VarType)

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
         SUBROUTINE InterpolateMeshToMeshQ( OldMesh, NewMesh, OldVariables, NewVariables, &
             UseQuadrantTree, Projector, MaskName, FoundNodes, NewMaskPerm, KeepUnfoundNodes )
           USE Types
           TYPE(Variable_t), POINTER, OPTIONAL :: OldVariables, NewVariables
           TYPE(Mesh_t), TARGET  :: OldMesh, NewMesh
           LOGICAL, OPTIONAL :: UseQuadrantTree,FoundNodes(:)
           CHARACTER(LEN=*),OPTIONAL :: MaskName
           TYPE(Projector_t), POINTER, OPTIONAL :: Projector
           INTEGER, OPTIONAL, POINTER :: NewMaskPerm(:)  !< Mask the new variable set by the given MaskName when trying to define the interpolation.
           LOGICAL, OPTIONAL :: KeepUnfoundNodes  !< Do not disregard unfound nodes from projector
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
                 CALL Fatal("VariableGet","Failed to find variable "//TRIM(Name))
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
              CALL Fatal("VariableGet","Failed to find or interpolate variable: "//TRIM(Name))
            END IF
         END IF
         RETURN
      END IF

!------------------------------------------------------------------------------

      IF ( .NOT.ASSOCIATED( Tmp ) ) THEN
         !GlobalBubbles = ListGetLogical(Pvar % Solver % Values, &
         !      'Bubbles in Global System', Found)
         !IF (.NOT.Found) GlobalBubbles=.TRUE.
        GlobalBubbles = Pvar % Solver % GlobalBubbles
        
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
        CALL Info('VariableGet','Performing masked on-the-fly interpolation',Level=15)
        CALL InterpolateMeshToMesh( PVar % PrimaryMesh, &
            CurrentModel % Mesh, Var, Variables, MaskName=MaskName )
      ELSE IF( UseProjector ) THEN
        CALL Info('VariableGet','Performing interpolation with projector',Level=15)
        CALL InterpolateMeshToMesh( PVar % PrimaryMesh, &
            CurrentModel % Mesh, Var, Variables, Projector=Projector )
      ELSE
        CALL Info('VariableGet','Performing on-the-fly interpolation',Level=15)
        AidVar => VariableGet( CurrentModel % Mesh % Variables, Name, ThisOnly = .TRUE. ) 
        IF( ASSOCIATED( AidVar ) ) THEN
          AidVar % Values = 0.0_dp
        END IF
        CALL InterpolateMeshToMesh( PVar % PrimaryMesh, &
            CurrentModel % Mesh, Var, Variables )        
      END IF
     
      IF( InfoActive( 20 ) ) THEN
        AidVar => VariableGet( CurrentModel % Mesh % Variables, Name, ThisOnly = .TRUE. )         
        PRINT *,'Interpolation range:',TRIM(AidVar % Name),MINVAL(AidVar % Values),MAXVAL( AidVar % Values)
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

#ifdef DEVEL_LISTCOUNTER
     IF( ASSOCIATED( new ) ) new % Counter = new % Counter + 1
#endif


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

     CALL Info('ListSetNamespace','Setting namespace to: '//TRIM(str_lcase),Level=15)
     
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
    IF (ALLOCATED(Namespace)) THEN
      l = .TRUE.
      str = Namespace
    ELSE
      l = .FALSE.
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

     CALL Info('ListPushNameSpace','Adding name space: '//TRIM(str),Level=12)

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
   SUBROUTINE ListPopNamespace( str0 )
!------------------------------------------------------------------------------
     CHARACTER(LEN=*), OPTIONAL :: str0
     TYPE(String_stack_t), POINTER :: stack
     

     IF(ASSOCIATED(Namespace_stack)) THEN

       ! This is an optional part aimed to help to code correctly the name stack.
       ! If one gives the namespace to be popped a Fatal will result if it is a
       ! wrong namespace.
       IF( PRESENT( str0 ) ) THEN
         IF( str0 /= Namespace ) THEN
           CALL Fatal('ListPopNamespace','Wrong namespace to pop: '&
               //TRIM(str0)//' vs '//TRIM(Namespace))
         END IF
       END IF

       Namespace = Namespace_stack % name

       CALL Info('ListPopNameSpace','Deleting entry from name space: '&
           //TRIM(Namespace),Level=12)      

       stack => Namespace_stack
       Namespace_stack => stack % Next
       DEALLOCATE(stack)
     ELSE
       CALL Info('ListPopNameSpace','No namespace entry to delete',Level=20)
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
    IF (ALLOCATED(ActiveListName)) THEN
      str = ActiveListName
    ELSE
      str = ''
    END IF
!------------------------------------------------------------------------------
   END FUNCTION ListGetActiveName
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE SetNamespaceCheck(L)
!------------------------------------------------------------------------------
     LOGICAL :: L
!------------------------------------------------------------------------------
     DoNamespaceCheck = L
!------------------------------------------------------------------------------
   END SUBROUTINE SetNamespaceCheck
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   FUNCTION GetNamespaceCheck() RESULT(L)
!------------------------------------------------------------------------------
     LOGICAL :: L
!------------------------------------------------------------------------------
     L = DoNameSpaceCheck
!------------------------------------------------------------------------------
   END FUNCTION GetNamespaceCheck
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Finds an entry in the list by its name and returns a handle to it.
!------------------------------------------------------------------------------
   FUNCTION ListFind( list, name, Found ) RESULT(ptr)
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: name
     LOGICAL, OPTIONAL :: Found
!------------------------------------------------------------------------------
     TYPE(String_stack_t), POINTER :: stack
#ifdef HAVE_LUA
     CHARACTER(:), ALLOCATABLE :: stra
#endif
     CHARACTER(:), ALLOCATABLE :: strn
     CHARACTER(LEN=LEN_TRIM(Name)) :: str
!------------------------------------------------------------------------------
     INTEGER :: k, k1, n

     IF(PRESENT(Found)) Found = .FALSE.
     ptr => NULL()
     IF(.NOT.ASSOCIATED(List)) RETURN

     k = StringToLowerCase( str,Name,.TRUE. )

     IF( ListGetnamespace(strn) ) THEN
       stack => Namespace_stack
       DO WHILE(.TRUE.)
#ifdef HAVE_LUA
         stra = trim(strn)
         strn = stra //' '//str(1:k)
         DEALLOCATE(stra)
#else
         strn = trim(strn) // ' ' //str(1:k)
#endif
         k1 = LEN(strn)
         ptr => List % Head
         DO WHILE( ASSOCIATED(ptr) )
            n = ptr % NameLen
            IF ( n==k1 ) THEN
              IF ( ptr % Name(1:n) == strn ) EXIT
            END IF
            ptr => ptr % Next
         END DO
         IF(.NOT.DoNamespaceCheck) EXIT

         IF(ASSOCIATED(ptr).OR..NOT.ASSOCIATED(stack)) EXIT
         IF(stack % name=='') EXIT
         strn = char(stack % name)
         stack => stack % next
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

#ifdef DEVEL_LISTCOUNTER
     IF( ASSOCIATED( ptr ) ) THEN
       ptr % Counter = ptr % Counter + 1
     ELSE IF( INDEX( name,': not found' ) == 0 ) THEN
       CALL ListAddNewLogical( CurrentModel % Simulation, TRIM(name)//': not found',.TRUE.) 
       ! This seems to cause problems so we use the one above
       !       CALL ListAddNewLogical( List, TRIM(name)//': not found',.TRUE.) 
     END IF
#endif
     
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
!------------------------------------------------------------------------------
   SUBROUTINE ListRename( list, name, name2, Found ) 
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: name, name2
     LOGICAL, OPTIONAL :: Found
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
     CHARACTER(:), ALLOCATABLE :: strn
     CHARACTER(LEN=LEN_TRIM(Name)) :: str
     CHARACTER(LEN=LEN_TRIM(Name2)) :: str2
     INTEGER :: k, k2, n

     IF(PRESENT(Found)) Found = .FALSE.

     ptr => NULL()
     IF(.NOT.ASSOCIATED(List)) RETURN
     
     k = StringToLowerCase( str,Name,.TRUE. )
     
     Ptr => List % Head
     DO WHILE( ASSOCIATED(ptr) )
       n = ptr % NameLen
       IF ( n==k ) THEN
         IF ( ptr % Name(1:n) == str(1:n) ) EXIT
       END IF
       ptr => ptr % Next
     END DO
     
     IF( ASSOCIATED( ptr ) ) THEN
       k2 = StringToLowerCase( str2,Name2,.TRUE. )
       ptr % Name(1:k2) = str2(1:k2)
       ptr % NameLen = k2 
       !PRINT *,'renaming >'//str(1:k)//'< to >'//str2(1:k2)//'<', k, k2
     END IF
          
     IF ( PRESENT(Found) ) THEN
       Found = ASSOCIATED(ptr)
     ELSE IF (.NOT.ASSOCIATED(ptr) ) THEN
       CALL Warn( 'ListRename', ' ' )
       WRITE(Message,*) 'Requested property: ', '[',TRIM(Name),'], not found'
       CALL Warn( 'ListRename', Message )
       CALL Warn( 'ListRename', ' ' )
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE ListRename
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Rename all given keywords in BC section.
!------------------------------------------------------------------------------
   SUBROUTINE ListRenameAllBC( Model, Name, Name2 ) 
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name, Name2
     LOGICAL :: Found
     INTEGER :: bc, n

     n = 0
     DO bc = 1,Model % NumberOfBCs
       CALL ListRename( Model % BCs(bc) % Values, Name, Name2, Found )
       IF( Found ) n = n + 1
     END DO
     IF( n > 0 ) CALL Info('ListRenameAllBCs',&
         TRIM(Name)//' ranamed to '//TRIM(Name2)//' on '//TRIM(I2S(n))//' BCs',Level=6)
     
!------------------------------------------------------------------------------
   END SUBROUTINE ListRenameAllBC
!------------------------------------------------------------------------------

   
  

!-----------------------------------------------------------------------------
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
     TYPE(String_stack_t), POINTER :: stack
#ifdef HAVE_LUA
     CHARACTER(:), ALLOCATABLE :: stra
#endif
     CHARACTER(:), ALLOCATABLE :: strn
     CHARACTER(LEN=LEN_TRIM(Name)) :: str
!------------------------------------------------------------------------------
     INTEGER :: k, k1, n, m

     ptr => NULL()
     IF(.NOT.ASSOCIATED(List)) RETURN

     k = StringToLowerCase( str,Name,.TRUE. )
     IF ( ListGetNamespace(strn) ) THEN
       stack => Namespace_stack
       DO WHILE(.TRUE.)
#ifdef HAVE_LUA
         stra = trim(strn)
         strn = stra //' '//str(1:k)
         DEALLOCATE(stra)
#else
         strn = trim(strn) // ' ' //str(1:k)
#endif
         k1 = LEN(strn)
         ptr => List % Head
         DO WHILE( ASSOCIATED(ptr) )
            n = ptr % NameLen
            IF ( n >= k1 ) THEN
              IF ( ptr % Name(1:k1) == strn ) EXIT
            END IF
            ptr => ptr % Next
         END DO
         IF(.NOT.DoNamespaceCheck) EXIT

         IF(ASSOCIATED(ptr).OR..NOT.ASSOCIATED(stack)) EXIT
         IF(stack % name=='') EXIT
         strn = char(stack % name)
         stack => stack % next
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
     TYPE(String_stack_t), POINTER :: stack
     CHARACTER(:), ALLOCATABLE :: strn
     CHARACTER(LEN=LEN_TRIM(Name)) :: str
!------------------------------------------------------------------------------
     INTEGER :: k, k1, n, m

     ptr => NULL()
     IF(.NOT.ASSOCIATED(List)) RETURN

     k = StringToLowerCase( str,Name,.TRUE. )

     IF ( ListGetNamespace(strn) ) THEN
       stack => Namespace_stack
       DO WHILE(.TRUE.)
         strn = TRIM(strn) //' '//str(1:k)
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
         IF(.NOT.DoNamespaceCheck) EXIT

         IF(ASSOCIATED(ptr).OR..NOT.ASSOCIATED(stack)) EXIT
         IF(stack % name=='') EXIT
         strn = char(stack % name)
         stack => stack % next
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
   SUBROUTINE ListCopyItem( ptr, list, name )

     TYPE(ValueListEntry_t), POINTER :: ptr
     TYPE(ValueList_t), POINTER :: list
     CHARACTER(LEN=*), OPTIONAL :: name
!------------------------------------------------------------------------------
     INTEGER :: i,j,k
     TYPE(ValueListEntry_t), POINTER :: ptrb, ptrnext

     IF( PRESENT( name ) ) THEN
       ptrb => ListAdd( List, name ) 
     ELSE
       ptrb => ListAdd( List, ptr % Name ) 
     END IF

       
     ptrnext => ptrb % next
     ptrb = ptr

     ptrb % tvalues => null()
     if(associated(ptr % tvalues)) then
       allocate( ptrb % tvalues(size(ptr % tvalues)) )
       ptrb % tvalues = ptr % tvalues
     end if

     ptrb % fvalues => null()
     if(associated(ptr % fvalues)) then
       i = size(ptr % fvalues,1)
       j = size(ptr % fvalues,2)
       k = size(ptr % fvalues,3)
       allocate( ptrb % fvalues(i,j,k) )
       ptrb % fvalues = ptr % fvalues
     end if

     ptrb % ivalues => null()
     if(associated(ptr % ivalues)) then
       allocate( ptrb % ivalues(size(ptr % ivalues)) )
       ptrb % ivalues = ptr % ivalues
     end if

     ptrb % cumulative => null()
     if(associated(ptr % cumulative)) then
       allocate( ptrb % cumulative(size(ptr % cumulative)) )
       ptrb % cumulative = ptr % cumulative
     end if
     ptrb % next => ptrnext

     ! If name is given then we have to revert the stuff from previous lines
     IF( PRESENT( name ) ) THEN
       ptrb % Name = name
       ptrb % Namelen = lentrim( name )
     END IF
     
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
 

!> Goes through one list and checks whether it includes any keywords with give prefix.
!> All keywords found are copied to the 2nd list without the prefix.
!------------------------------------------------------------------------------
   SUBROUTINE ListCopyPrefixedKeywords( list, listb, prefix )
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: list, listb
     CHARACTER(LEN=*) :: prefix
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
     CHARACTER(LEN=LEN_TRIM(prefix)) :: str
     INTEGER :: k, l, n, ncopy

     k = StringToLowerCase( str,prefix,.TRUE. )
     ncopy = 0
     
     ! Find the keyword from the 1st list 
     Ptr => List % Head
     DO WHILE( ASSOCIATED(ptr) )
       n = ptr % NameLen
       IF( n > k ) THEN
         IF( ptr % Name(1:k) == str(1:k) ) THEN
           l = k+1
           ! Remove the extra blanco after prefix if present
           ! Here we just assume one possible blanco as that is most often the case
           IF( ptr % Name(l:l) == ' ') l = l+1
           CALL Info('ListCopyPrefixedKeywords',&
               'Prefix: '//TRIM(prefix)// ' Keyword: '//TRIM(ptr % Name(l:n)),Level=12)
           CALL ListCopyItem( ptr, listb, ptr % Name(l:n) )
           ncopy = ncopy + 1
         END IF
       END IF
       ptr => ptr % Next
     END DO

     IF( ncopy > 0 ) THEN
       CALL Info('ListCopyPrefixedKeywords',&
           'Copied '//TRIM(I2S(ncopy))//' keywords with prefix: '//TRIM(prefix),Level=6)
     END IF
     
   END SUBROUTINE ListCopyPrefixedKeywords


!> Goes through one list and copies all keywords to a second list.
!------------------------------------------------------------------------------
   SUBROUTINE ListCopyAllKeywords( list, listb )
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: list, listb
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
     INTEGER :: ncopy

     ncopy = 0
     
     ! Find the keyword from the 1st list 
     Ptr => List % Head
     DO WHILE( ASSOCIATED(ptr) )
       CALL ListCopyItem( ptr, listb, ptr % Name )
       ncopy = ncopy + 1
       ptr => ptr % Next
     END DO
     
     IF( ncopy > 0 ) THEN
       CALL Info('ListCopyAllKeywords',&
           'Copied '//TRIM(I2S(ncopy))//' keywords to new list',Level=6)
     END IF
     
   END SUBROUTINE ListCopyAllKeywords
 
 
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
!> Just checks if there is a untreated keyword in the routine in the list.
!> In case there is return a warning. 
!------------------------------------------------------------------------------
   SUBROUTINE ListUntreatedWarn( List, Name, Caller ) 
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     CHARACTER(LEN=*), OPTIONAL :: Caller
!------------------------------------------------------------------------------
     IF( ListCheckPresent( List, Name ) ) THEN
       IF( PRESENT( Caller ) ) THEN
         CALL Warn(Caller,'Untreated keyword may cause problems: '//TRIM(Name))
       ELSE
         CALL Warn('ListUntreatedWarn','Untreated keyword may cause problems: '//TRIM(Name))
       END IF
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE ListUntreatedWarn
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Just checks if there is a untreated keyword in the routine in the list.
!> In case there is return a Fatal. 
!------------------------------------------------------------------------------
   SUBROUTINE ListUntreatedFatal( List, Name, Caller ) 
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     CHARACTER(LEN=*), OPTIONAL :: Caller
!------------------------------------------------------------------------------
     IF( ListCheckPresent( List, Name ) ) THEN
       IF( PRESENT( Caller ) ) THEN
         CALL Fatal(Caller,'Untreated keyword: '//TRIM(Name))
       ELSE
         CALL Fatal('ListUntreatedFatal','Untreated keyword: '//TRIM(Name))
       END IF
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE ListUntreatedFatal
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

      ptr % TYPE = LIST_TYPE_CONSTANT_TENSOR
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

     ! The (x,y) table should be such that values of x are increasing in size
     IF( .NOT. CheckMonotone( n, TValues ) ) THEN
       CALL Fatal('ListAddDepReal',&
           'Values x in > '//TRIM(Name)//' < not monotonically ordered!')
     END IF
     
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

      ptr % Fdim = 0
      IF( N > 1 ) ptr % Fdim = 1
      IF( M > 1 ) ptr % Fdim = ptr % Fdim + 1
      
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

     ptr % fdim = 0
     IF( n1 > 1 ) ptr % fdim = 1
     IF( n2 > 1 ) ptr % fdim = ptr % fdim + 1
     
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
!> Adds a logical entry to the list if it does not exist previously.
!------------------------------------------------------------------------------
   SUBROUTINE ListAddNewLogical( List,Name,LValue )
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL :: LValue
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------
     IF( ListCheckPresent( List, Name ) ) RETURN

     CALL ListAddLogical( List,Name,LValue )

   END SUBROUTINE ListAddNewLogical
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Adds an integer to the list when not present previously.
!------------------------------------------------------------------------------
    SUBROUTINE ListAddNewInteger( List,Name,IValue,Proc )
!------------------------------------------------------------------------------
      TYPE(ValueList_t), POINTER :: List
      CHARACTER(LEN=*) :: Name
      INTEGER :: IValue
      INTEGER(Kind=AddrInt), OPTIONAL :: Proc
!------------------------------------------------------------------------------
      TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------
      IF( ListCheckPresent( List, Name ) ) RETURN

      CALL ListAddInteger( List,Name,IValue,Proc )

    END SUBROUTINE ListAddNewInteger
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Adds a constant real value to the list if not present.
!------------------------------------------------------------------------------
    SUBROUTINE ListAddNewConstReal( List,Name,FValue,Proc,CValue )
!------------------------------------------------------------------------------
      TYPE(ValueList_t), POINTER :: List
      CHARACTER(LEN=*) :: Name
      CHARACTER(LEN=*), OPTIONAL :: Cvalue
      REAL(KIND=dp) :: FValue
      INTEGER(KIND=AddrInt), OPTIONAL :: Proc
!------------------------------------------------------------------------------
      TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------
      IF( ListCheckPresent( List, Name ) ) RETURN

      CALL ListAddConstReal( List,Name,FValue,Proc,CValue )

    END SUBROUTINE ListAddNewConstReal
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Add a string value to the list if not present.
!------------------------------------------------------------------------------
    SUBROUTINE ListAddNewString( List,Name,CValue,CaseConversion )
!------------------------------------------------------------------------------
      TYPE(ValueList_t), POINTER :: List
      CHARACTER(LEN=*) :: Name
      CHARACTER(LEN=*) :: CValue
      LOGICAL, OPTIONAL :: CaseConversion
      
      IF( ListCheckPresent( List, Name ) ) RETURN

      CALL ListAddString( List,Name,CValue,CaseConversion )

    END SUBROUTINE ListAddNewString
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
           CALL Fatal("ListGetIntegerArray", Message)
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
           CALL Fatal("ListGetLogical", Message)
         END IF
       END IF
       RETURN
     END IF
     L = ptr % Lvalue
!------------------------------------------------------------------------------
   END FUNCTION ListGetLogical
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> A generalized version of ListGetLogical. Uses logical, only if the keyword is
!> of type locical, otherwise returns True if the keyword is present.
!> Since the absence if a sign of False there is no separate Found flag.
!------------------------------------------------------------------------------
   RECURSIVE FUNCTION ListGetLogicalGen( List, Name) RESULT(L)
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL :: L
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
     LOGICAL :: Found
!------------------------------------------------------------------------------

     L = .FALSE.
     ptr => ListFind(List,Name,Found)
     IF ( ASSOCIATED(ptr) ) THEN
       IF(ptr % TYPE == LIST_TYPE_CONSTANT_SCALAR ) THEN
         L = ptr % Lvalue
       ELSE
         L = .TRUE.
       END IF
     END IF    
!------------------------------------------------------------------------------
   END FUNCTION ListGetLogicalGen
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
           CALL Fatal("ListGetString", Message)
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
           WRITE(Message, '(A,A)') "Failed to find constant real: ",Name
           CALL Fatal("ListGetConstReal", Message)
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


!> Get pointer to list of section
!------------------------------------------------------------------------------
  FUNCTION ListGetSection( Element, SectionName, Found ) RESULT(lst)
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER  :: Lst
    CHARACTER(LEN=*) :: SectionName
    LOGICAL, OPTIONAL :: Found
    TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER  :: BodyLst
    INTEGER :: id
    LOGICAL :: LFound
    
    id = Element % BodyId
    IF( id > 0 ) THEN
      bodylst => CurrentModel % Bodies(id) % Values
    ELSE
      NULLIFY( bodylst ) 
    END IF
    LFound = .FALSE.
    
    SELECT CASE ( SectionName ) 

    CASE( 'body' )     
      lst => bodylst
      Lfound = ASSOCIATED( lst ) 
               
    CASE( 'material' )
      id = ListGetInteger( bodylst, SectionName, LFound )
      IF( LFound ) lst => CurrentModel % Materials(id) % Values

    CASE( 'body force' )
      id = ListGetInteger( bodylst, SectionName, LFound )
      IF( LFound ) lst => CurrentModel % BodyForces(id) % Values
      
    CASE( 'initial condition' )
      id = ListGetInteger( bodylst, SectionName, LFound )
      IF( LFound ) lst => CurrentModel % ICs(id) % Values

    CASE( 'equation' )
      id = ListGetInteger( bodylst, SectionName, LFound )
      IF( LFound ) lst => CurrentModel % Equations(id) % Values

    CASE( 'boundary condition' )
      IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
        id = Element % BoundaryInfo % Constraint
        IF( id > 0 ) THEN
          lst => CurrentModel % BCs(id) % Values
          LFound = .TRUE.
        END IF
      END IF
        
    CASE DEFAULT
      CALL Fatal('ListGetSection','Unknown section name: '//TRIM(SectionName))
            
    END SELECT

    IF( PRESENT( Found ) ) Found = LFound 

!------------------------------------------------------------------------------
  END FUNCTION ListGetSection
!------------------------------------------------------------------------------


!> Get pointer to list of section
!------------------------------------------------------------------------------
  FUNCTION ListGetSectionId( Element, SectionName, Found ) RESULT(id)
!------------------------------------------------------------------------------
    INTEGER :: id
    CHARACTER(LEN=*) :: SectionName
    LOGICAL, OPTIONAL :: Found
    TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER  :: BodyLst
    INTEGER :: body_id
    LOGICAL :: LFound

    id = 0
    
    body_id = Element % BodyId
    IF( body_id > 0 ) THEN
      bodylst => CurrentModel % Bodies(body_id) % Values
    ELSE
      NULLIFY( bodylst ) 
    END IF
    LFound = .FALSE.
    
    SELECT CASE ( SectionName ) 

    CASE( 'body' )     
      id = body_id
               
    CASE( 'material' )
      id = ListGetInteger( bodylst, SectionName, LFound )

    CASE( 'body force' )
      id = ListGetInteger( bodylst, SectionName, LFound )
      
    CASE( 'initial condition' )
      id = ListGetInteger( bodylst, SectionName, LFound )

    CASE( 'equation' )
      id = ListGetInteger( bodylst, SectionName, LFound )

    CASE( 'boundary condition' )
      IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
        id = Element % BoundaryInfo % Constraint
      END IF
        
    CASE DEFAULT
      CALL Fatal('ListGetSection','Unknown section name: '//TRIM(SectionName))
            
    END SELECT

    IF( PRESENT( Found ) ) Found = ( id > 0 ) 

!------------------------------------------------------------------------------
  END FUNCTION ListGetSectionId
!------------------------------------------------------------------------------


  


  
!------------------------------------------------------------------------------
!> Given a string containing comma-separated variablenames, reads the strings
!> and obtains the corresponding variables to a table.
!------------------------------------------------------------------------------
  SUBROUTINE ListParseStrToVars( str, slen, name, count, VarTable, &
      SomeAtIp, SomeAtNodes, AllGlobal )
!------------------------------------------------------------------------------
     CHARACTER(LEN=*) :: str, name
     INTEGER :: slen, count
     TYPE(VariableTable_t) :: VarTable(:)
     LOGICAL :: SomeAtIp, SomeAtNodes, AllGlobal
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,n,k1,l,l0,l1
     TYPE(Variable_t), POINTER :: Var

     SomeAtIp = .FALSE.
     SomeAtNodes = .FALSE.
     AllGlobal = .TRUE.
     
     count=0
     l0=1
     IF(slen<=0) RETURN

     DO WHILE( .TRUE. )
       ! Remove zeros ahead
       DO WHILE( str(l0:l0) == ' ' )
         l0 = l0 + 1
         IF ( l0 > slen ) EXIT
       END DO
       IF ( l0 > slen ) EXIT

       ! Scan only until next comma
       l1 = INDEX( str(l0:slen),',')
       IF ( l1 > 0 ) THEN
         l1=l0+l1-2
       ELSE
         l1=slen
       END IF

       IF ( str(l0:l1) == 'coordinate' ) THEN
         VarTable(count+1) % Variable => VariableGet( CurrentModel % Variables,"coordinate 1")
         VarTable(count+2) % Variable => VariableGet( CurrentModel % Variables,"coordinate 2")
         VarTable(count+3) % Variable => VariableGet( CurrentModel % Variables,"coordinate 3")
         count = count + 3 
         SomeAtNodes = .TRUE.
         AllGlobal = .FALSE.
       ELSE
         Var => VariableGet( CurrentModel % Variables,TRIM(str(l0:l1)) )
         IF ( .NOT. ASSOCIATED( Var ) ) THEN
           CALL Info('ListParseStrToVars','Parsed variable '//TRIM(I2S(count+1))//' of '//str(1:slen),Level=3)
           CALL Info('ListParseStrToVars','Parse counters: '&
               //TRIM(I2S(l0))//', '//TRIM(I2S(l1))//', '//TRIM(I2S(slen)),Level=10)
           CALL Fatal('ListParseStrToVars', 'Can''t find independent variable:['// &
               TRIM(str(l0:l1))//'] for dependent variable:['//TRIM(Name)//']' ) 
         END IF
         count = count + 1
         VarTable(count) % Variable => Var
     
         IF( SIZE( Var % Values ) > Var % Dofs ) AllGlobal = .FALSE.
                
         IF( Var % TYPE == Variable_on_gauss_points ) THEN
           SomeAtIp = .TRUE.
         ELSE
           SomeAtNodes = .TRUE.
         END IF
         
       END IF

       ! New start after the comma
       l0 = l1+2
       IF ( l0 > slen ) EXIT       
     END DO
     
!------------------------------------------------------------------------------
   END SUBROUTINE ListParseStrToVars
!------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------
!> Given a table of variables and a node index return the variable values on the node.
!-------------------------------------------------------------------------------------
  SUBROUTINE VarsToValuesOnNodes( VarCount, VarTable, ind, T, count )
!------------------------------------------------------------------------------
     INTEGER :: Varcount
     TYPE(VariableTable_t) :: VarTable(:)
     INTEGER :: ind
     INTEGER :: count
     REAL(KIND=dp) :: T(:)
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
     INTEGER :: i,j,k,n,k1,l,varsize,vari
     TYPE(Variable_t), POINTER :: Var
     LOGICAL :: Failed
     
     count = 0
     Failed = .FALSE.
     
     DO Vari = 1, VarCount
       
       Var => VarTable(Vari) % Variable

       Varsize = SIZE( Var % Values ) / Var % Dofs 
       
       IF( Varsize == 1 ) THEN
         DO l=1,Var % DOFs
           count = count + 1
           T(count) = Var % Values(l)
         END DO
       ELSE
         k1 = ind
         
         IF ( Var % TYPE == Variable_on_gauss_points ) THEN
           count = count + Var % DOFs
           CYCLE
         ELSE IF( Var % TYPE == Variable_on_elements ) THEN
           Element => CurrentModel % CurrentElement
           IF( ASSOCIATED( Element ) ) k1 = Element % ElementIndex
         ELSE IF ( Var % TYPE == Variable_on_nodes_on_elements ) THEN
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

         IF ( ASSOCIATED(Var % Perm) ) k1 = Var % Perm(k1)         
         
         IF ( k1 > 0 .AND. k1 <= VarSize ) THEN
           DO l=1,Var % DOFs
             count = count + 1
             T(count) = Var % Values(Var % Dofs*(k1-1)+l)
           END DO
         ELSE
           Failed = .TRUE.
           DO l=1,Var % DOFs
             count = count + 1
             T(count) = HUGE(1.0_dp)           
           END DO
           RETURN
         END IF
       END IF
     END DO
     
   END SUBROUTINE VarsToValuesOnNodes
 !------------------------------------------------------------------------------
    

!-------------------------------------------------------------------------------------
!> Given a table of variables return the variable values on the gauss point.
!> This only deals with the gauss point variables, all other are already treated. 
!-------------------------------------------------------------------------------------
  SUBROUTINE VarsToValuesOnIps( VarCount, VarTable, ind, T, count )
!------------------------------------------------------------------------------
     INTEGER :: Varcount
     TYPE(VariableTable_t) :: VarTable(:)
     INTEGER :: ind
     INTEGER :: count
     REAL(KIND=dp) :: T(:)
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
     INTEGER :: i,j,k,n,k1,l,varsize,vari
     TYPE(Variable_t), POINTER :: Var
     LOGICAL :: Failed
     
     count = 0
     Failed = .FALSE.
     
     DO Vari = 1, VarCount 
       Var => VarTable(Vari) % Variable
       Varsize = SIZE( Var % Values ) / Var % Dofs 

       k1 = 0
       IF ( Var % TYPE == Variable_on_gauss_points ) THEN         
         Element => CurrentModel % CurrentElement
         IF ( ASSOCIATED(Element) ) THEN
           k1 = Var % Perm( Element % ElementIndex ) + ind
         END IF
       END IF
         
       IF ( k1 > 0 .AND. k1 <= VarSize ) THEN
         DO l=1,Var % DOFs
           count = count + 1
           T(count) = Var % Values(Var % Dofs*(k1-1)+l)
         END DO
       ELSE
         count = count + Var % Dofs
       END IF
     END DO
     
   END SUBROUTINE VarsToValuesOnIps
 !------------------------------------------------------------------------------


   
!------------------------------------------------------------------------------
  SUBROUTINE ListParseStrToValues( str, slen, ind, name, T, count, AllGlobal    )
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
           CALL Info('ListParseStrToValues','Parsed variable '//TRIM(I2S(count+1))//' of '//str(1:slen),Level=3)
           CALL Info('ListParseStrToValues','Parse counters: '&
               //TRIM(I2S(l0))//', '//TRIM(I2S(l1))//', '//TRIM(I2S(slen)),Level=10)
           CALL Fatal('ListParseStrToValues','Can''t find independent variable:['// &
               TRIM(str(l0:l1))//'] for dependent variable:['//TRIM(Name)//']')
         END IF
         IF( SIZE( Variable % Values ) > Variable % Dofs ) AllGlobal = .FALSE.
       ELSE
         AllGlobal = .FALSE.
         Variable => VariableGet( CurrentModel % Variables,'Coordinate 1' )         
       END IF
       
       IF( Variable % TYPE == Variable_on_gauss_points ) THEN
         DO l=1,Variable % DOFs
           count = count + 1
           T(count) = HUGE(1.0_dp)
         END DO

         l0 = l1+2
         IF ( l0 > slen ) EXIT
         CYCLE
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
  FUNCTION ListCheckGlobal( ptr ) RESULT ( IsGlobal )
!------------------------------------------------------------------------------
     TYPE(ValueListEntry_t), POINTER :: ptr
     LOGICAL :: IsGlobal
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
     INTEGER :: ind,i,j,k,n,k1,l,l0,l1,ll,count
     TYPE(Variable_t), POINTER :: Variable, CVar
     INTEGER :: slen

     IsGlobal = .TRUE.

     IF(.NOT.ASSOCIATED(ptr)) THEN
       CALL Warn('ListCheckGlobal','ptr not associated!')
       RETURN
     END IF
       
     
     IF( ptr % TYPE == LIST_TYPE_CONSTANT_SCALAR_STR ) THEN
       RETURN

     ELSE IF( ptr % TYPE == LIST_TYPE_CONSTANT_SCALAR .OR. & 
         ptr % TYPE == LIST_TYPE_VARIABLE_SCALAR .OR. &
         ptr % TYPE == LIST_TYPE_VARIABLE_SCALAR_STR ) THEN


       IF ( ptr % PROCEDURE /= 0 ) THEN
         IsGlobal = .FALSE.
         RETURN
       END IF

       slen = ptr % DepNameLen

       IF( slen >  0 ) THEN
         count = 0
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

           count = count + 1

           IF ( ptr % DependName(l0:l1) /= 'coordinate' ) THEN
             Variable => VariableGet( CurrentModel % Variables,TRIM(ptr % DependName(l0:l1)) )
             IF ( .NOT. ASSOCIATED( Variable ) ) THEN             
               CALL Info('ListCheckGlobal','Parsed variable '//TRIM(I2S(count))//' of '&
                   //ptr % DependName(1:slen),Level=3)
               CALL Info('ListCheckGlobal','Parse counters: '&
                   //TRIM(I2S(l0))//', '//TRIM(I2S(l1))//', '//TRIM(I2S(slen)),Level=10)

               WRITE( Message, * ) 'Can''t find independent variable:[', &
                   TRIM(ptr % DependName(l0:l1)),']'
               CALL Fatal( 'ListCheckGlobal', Message )
             END IF

             IF( SIZE( Variable % Values ) > 1 ) THEN
               IsGlobal = .FALSE.
               RETURN
             END IF

           ELSE
             IsGlobal = .FALSE.
             EXIT
           END IF

           l0 = l1+2
           IF ( l0 > slen ) EXIT
         END DO
       ELSE
         IsGlobal = .FALSE.
       END IF
     END IF

       
!------------------------------------------------------------------------------
   END FUNCTION ListCheckGlobal
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

     AllGlobal = ListCheckGlobal( ptr )
    
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
     TYPE(VariableTable_t) :: VarTable(MAX_FNC)
     INTEGER :: i,j,k,k1,l,l0,l1,lsize, VarCount
     CHARACTER(LEN=MAX_NAME_LEN) ::  cmd, tmp_str
     LOGICAL :: AllGlobal, SomeAtIp, SomeAtNodes
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

       CALL ListParseStrToVars( Ptr % DependName, Ptr % DepNameLen, Name, VarCount, VarTable, &
           SomeAtIp, SomeAtNodes, AllGlobal )
       IF( SomeAtIp ) THEN
         CALL Fatal('ListGetReal','Function cannot deal with variables on IPs!')
       END IF

       DO i=1,n
         k = NodeIndexes(i)

         CALL VarsToValuesOnNodes( VarCount, VarTable, k, T, j )
         
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

       CALL ListParseStrToVars( Ptr % DependName, Ptr % DepNameLen, Name, VarCount, &
           VarTable, SomeAtIp, SomeAtNodes, AllGlobal )
       IF( SomeAtIp ) THEN
         CALL Fatal('ListGetReal','Function cannot deal with variables on IPs!')
       END IF
       
       
       DO i=1,n
         k = NodeIndexes(i)

         CALL VarsToValuesOnNodes( VarCount, VarTable, k, T, j )
         
#ifdef HAVE_LUA
         IF ( .not. ptr % LuaFun ) THEN
#endif
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

#ifdef HAVE_LUA
         ELSE
           call ElmerEvalLua(LuaState, ptr, T, F(i), varcount)
         END IF
#endif
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

   RECURSIVE SUBROUTINE ListInitHandle( Handle )

     TYPE(ValueHandle_t) :: Handle

     Handle % ValueType = -1
     Handle % SectionType = -1
     Handle % ListId = -1
     Handle % Element => NULL()
     Handle % List => NULL()
     Handle % Ptr  => NULL()
     Handle % Nodes => NULL()
     Handle % Indexes => NULL()
     Handle % nValuesVec = 0
     Handle % ValuesVec => NULL()
     Handle % Values => NULL()
     Handle % ParValues => NULL()
     Handle % ParNo = 0
     Handle % DefIValue = 0
     Handle % DefRValue = 0.0_dp
     Handle % Rdim = 0
     Handle % RTensor => NULL()
     Handle % RTensorValues => NULL()
     Handle % DefLValue = .FALSE.
     Handle % Initialized = .FALSE.
     Handle % AllocationsDone = .FALSE.
     Handle % ConstantEverywhere = .FALSE.
     Handle % GlobalEverywhere = .FALSE.
     Handle % GlobalInList = .FALSE.
     Handle % EvaluateAtIP = .FALSE.
     Handle % SomeVarAtIp = .FALSE.
     Handle % SomewhereEvaluateAtIP = .FALSE.
     Handle % NotPresentAnywhere = .FALSE.
     Handle % UnfoundFatal = .FALSE.
     Handle % GotMinv = .FALSE.
     Handle % GotMaxv = .FALSE.
     Handle % VarCount = 0
     Handle % HandleIm => NULL()
     Handle % Handle2 => NULL()
     Handle % Handle3 => NULL()
     
   END SUBROUTINE ListInitHandle


!------------------------------------------------------------------------------
!> Initializes the handle to save just a little bit for constant valued.
!> This is not mandatory but may still be used. 
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE ListInitElementKeyword( Handle,Section,Name,minv,maxv,&
       DefRValue,DefIValue,DefLValue,UnfoundFatal,EvaluateAtIp,&
       FoundSomewhere,InitIm,InitVec3D)
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     CHARACTER(LEN=*)  :: Section,Name
     REAL(KIND=dp), OPTIONAL :: minv,maxv
     REAL(KIND=dp), OPTIONAL :: DefRValue
     INTEGER, OPTIONAL :: DefIValue
     LOGICAL, OPTIONAL :: DefLValue
     LOGICAL, OPTIONAL :: UnfoundFatal
     LOGICAL, OPTIONAL :: EvaluateAtIp
     LOGICAL, OPTIONAL :: FoundSomewhere
     LOGICAL, OPTIONAL :: InitIm
     LOGICAL, OPTIONAL :: InitVec3D     
     !------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     TYPE(ValueListEntry_t), POINTER :: ptr
     INTEGER :: i, n, NoVal, ValueType, IValue, dim, n1, n2, maxn1, maxn2
     TYPE(Model_t), POINTER :: Model
     REAL(KIND=dp)  :: val, Rvalue
     CHARACTER(LEN=MAX_NAME_LEN) :: CValue
     LOGICAL :: ConstantEverywhere, NotPresentAnywhere, Lvalue, FirstList, Found
     REAL(KIND=dp), POINTER :: Basis(:)
     INTEGER, POINTER :: NodeIndexes(:)
     TYPE(Element_t), POINTER :: Element
     LOGICAL :: GotIt, FoundSomewhere1, FoundSomewhere2
     !------------------------------------------------------------------------------
     
     IF( PRESENT( InitIm ) ) THEN
       IF( InitIm ) THEN
         IF( .NOT. ASSOCIATED( Handle % HandleIm ) ) THEN
           ALLOCATE( Handle % HandleIm )
           CALL ListInitHandle( Handle % HandleIm ) 
        END IF
         CALL Info('ListInitElementKeyword','Treating real part of keyword',Level=12)         
         CALL ListInitElementKeyword( Handle,Section,Name,minv,maxv,&
             DefRValue,DefIValue,DefLValue,UnfoundFatal,EvaluateAtIp,FoundSomewhere,InitVec3D=InitVec3D)
         IF( PRESENT( FoundSomewhere) ) FoundSomewhere1 = FoundSomewhere
         
         CALL Info('ListInitElementKeyword','Treating imaginary part of keyword',Level=12)                 
         CALL ListInitElementKeyword( Handle % HandleIm,Section,TRIM(Name)//' im',minv,maxv,&
             DefRValue,DefIValue,DefLValue,UnfoundFatal,EvaluateAtIp,FoundSomewhere,InitVec3D=InitVec3D)
         IF( PRESENT( FoundSomewhere ) ) FoundSomewhere =  FoundSomewhere .OR. FoundSomewhere1
         RETURN
       END IF
     END IF

     IF( PRESENT( InitVec3D ) ) THEN
       IF( InitVec3D ) THEN
         IF( .NOT. ASSOCIATED( Handle % Handle2 ) ) THEN
           ALLOCATE( Handle % Handle2 )
           CALL ListInitHandle( Handle % Handle2 ) 
         END IF
         IF( .NOT. ASSOCIATED( Handle % Handle3 ) ) THEN           
           ALLOCATE( Handle % Handle3 )           
           CALL ListInitHandle( Handle % Handle3 ) 
         END IF

         CALL ListInitElementKeyword( Handle,Section,TRIM(Name)//' 1',minv,maxv,&
             DefRValue,DefIValue,DefLValue,UnfoundFatal,EvaluateAtIp,FoundSomewhere)
         IF( PRESENT( FoundSomewhere) ) FoundSomewhere1 = FoundSomewhere
         CALL ListInitElementKeyword( Handle % Handle2,Section,TRIM(Name)//' 2',minv,maxv,&
             DefRValue,DefIValue,DefLValue,UnfoundFatal,EvaluateAtIp,FoundSomewhere)
         IF( PRESENT( FoundSomewhere) ) FoundSomewhere2 = FoundSomewhere
         CALL ListInitElementKeyword( Handle % Handle3,Section,TRIM(Name)//' 3',minv,maxv,&
             DefRValue,DefIValue,DefLValue,UnfoundFatal,EvaluateAtIp,FoundSomewhere)         
         IF( PRESENT( FoundSomewhere ) ) FoundSomewhere = FoundSomewhere .OR. &
             FoundSomewhere1 .OR. FoundSomewhere2
         RETURN
       END IF
     END IF
     
     CALL Info('ListInitElementKeyword','Treating keyword: '//TRIM(Name),Level=10)

     Model => CurrentModel
     Handle % BulkElement = .TRUE.
     NULLIFY(ptr)
     
     SELECT CASE ( Section ) 

     CASE('Body')
       Handle % SectionType = SECTION_TYPE_BODY

     CASE('Material')
       Handle % SectionType = SECTION_TYPE_MATERIAL
       
     CASE('Body Force')
       Handle % SectionType = SECTION_TYPE_BF

     CASE('Initial Condition')
       Handle % SectionType = SECTION_TYPE_IC

     CASE('Boundary Condition')
       Handle % SectionType = SECTION_TYPE_BC
       Handle % BulkElement = .FALSE.
       
     CASE('Component')
       Handle % SectionType = SECTION_TYPE_COMPONENT

     CASE('Equation')
       Handle % SectionType = SECTION_TYPE_EQUATION

     CASE DEFAULT
       CALL Fatal('ListInitElementKeyword','Unknown section: '//TRIM(Section))

     END SELECT


     ! Initialize the handle entries because it may be that the list structure was altered,
     ! or the same handle is used for different keyword.
     Handle % ConstantEverywhere = .TRUE.
     Handle % GlobalInList = .FALSE.
     Handle % NotPresentAnywhere = .TRUE.
     Handle % SomewhereEvaluateAtIP = .FALSE.
     Handle % GlobalEverywhere = .TRUE.
     Handle % SomeVarAtIp = .FALSE.
     Handle % Name = Name 
     Handle % ListId = -1
     Handle % EvaluateAtIp = .FALSE.       
     Handle % List => NULL()
     Handle % Element => NULL()
     Handle % Unfoundfatal = .FALSE.
     IF (.NOT. ASSOCIATED( Ptr ) ) THEN
       Handle % Ptr => ListAllocate()
     END IF


     ! Deallocate stuff that may change in size, or is used as a marker for first element
     IF( Handle % nValuesVec > 0 ) THEN
       DEALLOCATE( Handle % ValuesVec )
       Handle % nValuesVec = 0
     END IF
     
     
     Handle % Initialized = .TRUE.
     
     FirstList = .TRUE.
     maxn1 = 0
     maxn2 = 0
     
     i = 0
     DO WHILE(.TRUE.) 
       i = i + 1

       SELECT CASE ( Handle % SectionType ) 

       CASE( SECTION_TYPE_BODY )
         IF(i > Model % NumberOfBodies ) EXIT
         List => Model % Bodies(i) % Values

       CASE( SECTION_TYPE_MATERIAL )
         IF(i > Model % NumberOfMaterials ) EXIT
         List => Model % Materials(i) % Values

       CASE( SECTION_TYPE_BF )
         IF(i > Model % NumberOfBodyForces ) EXIT        
         List => Model % BodyForces(i) % Values
         
       CASE( SECTION_TYPE_IC )
         IF( i > Model % NumberOfICs ) EXIT
         List => Model % ICs(i) % Values

       CASE( SECTION_TYPE_EQUATION )
         IF( i > Model % NumberOfEquations ) EXIT
         List => Model % Equations(i) % Values

       CASE( SECTION_TYPE_BC )
         IF( i > Model % NumberOfBCs ) EXIT        
         List => Model % BCs(i) % Values

         ! It is more difficult to make sure that the BC list is given for all BC elements.
         ! Therefore set this to .FALSE. always for BCs. 
         Handle % ConstantEverywhere = .FALSE.
         
       CASE DEFAULT
         CALL Fatal('ListInitElementKeyword','Unknown section: '//TRIM(I2S(Handle % SectionType)))

       END SELECT
  
       ! If the parameter is not defined in some list we cannot really be sure
       ! that it is intentionally used as a zero. Hence we cannot assume that the
       ! keyword is constant. 
       ptr => ListFind(List,Name,Found)
       Handle % ptr % Head => ptr
       
       IF ( .NOT. ASSOCIATED(ptr) ) THEN
         Handle % ConstantEverywhere = .FALSE.
         CYCLE
       ELSE IF( FirstList ) THEN
         Handle % NotPresentAnywhere = .FALSE.
         Handle % ValueType = ptr % Type
       END IF

       ValueType = ptr % TYPE 

       IF( ValueType == LIST_TYPE_LOGICAL ) THEN
         Lvalue = ptr % Lvalue

         IF( FirstList ) THEN
           Handle % LValue = LValue
         ELSE
           IF( XOR( Handle % LValue, LValue ) ) THEN
             Handle % ConstantEverywhere = .FALSE.
             EXIT
           END IF
         END IF

       ELSE IF( ValueType == LIST_TYPE_STRING ) THEN
         Cvalue = ptr % Cvalue
         IF( FirstList ) THEN
           Handle % CValueLen = len_trim(CValue)
           Handle % CValue = CValue(1:Handle % CValueLen) 
         ELSE IF( Handle % CValue(1:Handle % CValueLen) /= Cvalue ) THEN
           Handle % ConstantEverywhere = .FALSE.
           EXIT
         END IF

       ELSE IF( ValueType == LIST_TYPE_INTEGER ) THEN
         Ivalue = ptr % Ivalues(1)           
         IF( FirstList ) THEN
           Handle % IValue = Ivalue
         ELSE IF( Handle % IValue /= Ivalue ) THEN
           Handle % ConstantEverywhere = .FALSE.
           EXIT
         END IF

       ELSE IF( ValueType >= LIST_TYPE_CONSTANT_SCALAR .AND. &
           ValueType <= List_TYPE_CONSTANT_SCALAR_PROC ) THEN         
         
         IF(.NOT. ListCheckAllGlobal( Handle % ptr, name ) ) THEN
           Handle % GlobalEverywhere = .FALSE.
           Handle % ConstantEverywhere = .FALSE.           
           IF( ListGetLogical( List, TRIM( Handle % Name )//' At IP',GotIt ) ) THEN
             Handle % SomewhereEvaluateAtIp = .TRUE.
             EXIT
           END IF
         END IF

         IF( Handle % ConstantEverywhere ) THEN
           Rvalue = ListGetCReal( List,Name)
           ! and each list must have the same constant value
           IF( FirstList ) THEN
             Handle % RValue = Rvalue
           ELSE IF( ABS( Handle % RValue - Rvalue ) > TINY( RValue ) ) THEN
             Handle % ConstantEverywhere = .FALSE.
           END IF
         END IF

       ELSE IF( ValueType >= LIST_TYPE_CONSTANT_TENSOR .AND. &
           ValueType <= LIST_TYPE_VARIABLE_TENSOR_STR ) THEN
         
         Handle % GlobalEverywhere = .FALSE.
         Handle % ConstantEverywhere = .FALSE.           
         IF( ListGetLogical( List, TRIM( Handle % Name )//' At IP',GotIt ) ) THEN
           Handle % SomewhereEvaluateAtIp = .TRUE.
         END IF
         
         n1 = SIZE( ptr % FValues,1 ) 
         n2 = SIZE( ptr % FValues,2 )
         
         maxn1 = MAX( n1, maxn1 )
         maxn2 = MAX( n2, maxn2 )
       ELSE
         CALL Fatal('ListInitElementKeyword','Unknown value type: '//TRIM(I2S(ValueType)))

       END IF

       FirstList = .FALSE.
     END DO

     CALL Info('ListInitElementKeyword',&
         'Initiated handle for: > '//TRIM(Handle % Name)//' < of type: '// &
         TRIM(I2S(Handle % ValueType)),Level=10)

     IF( PRESENT( UnfoundFatal ) ) THEN
       Handle % Unfoundfatal = UnfoundFatal
       IF( Handle % UnfoundFatal .AND. Handle % NotPresentAnywhere ) THEN
         CALL Fatal('ListInitElementKeywords','Keyword required but not present: '&
             //TRIM(Handle % Name))
       END IF
     END IF
     
     IF( PRESENT( DefLValue ) ) THEN
       Handle % DefLValue = DefLValue
     END IF
     
     IF( PRESENT( DefRValue ) ) THEN
       Handle % DefRValue = DefRValue
     END IF

     IF( PRESENT( DefIValue ) ) THEN
       Handle % DefIValue = DefIValue
     END IF
     
     IF( PRESENT( minv ) ) THEN
       Handle % GotMinv = .TRUE.
       Handle % minv = minv
     END IF

     IF( PRESENT( maxv ) ) THEN
       Handle % GotMaxv = .TRUE.
       Handle % maxv = maxv
     END IF

     IF( PRESENT( EvaluateAtIp ) ) THEN
       Handle % EvaluateAtIp = EvaluateAtIp 
     END IF

     IF( PRESENT( FoundSomewhere ) ) THEN
       FoundSomewhere = .NOT. Handle % NotPresentAnywhere
     END IF

     ! For tensor valued ListGetRealElement operations allocate the maximum size
     ! of temporal table needed. 
     IF( maxn1 > 1 .OR. maxn2 > 1 ) THEN
       n = CurrentModel % Mesh % MaxElementNodes
       IF( ASSOCIATED( Handle % RtensorValues ) ) THEN
         IF( SIZE( Handle % RtensorValues, 1 ) < maxn1 .OR. &
             SIZE( Handle % RtensorValues, 2 ) < maxn2 .OR. &
             SIZE( Handle % RtensorValues, 3 ) < n ) THEN
           DEALLOCATE( Handle % RtensorValues )
         END IF
       END IF
       IF(.NOT. ASSOCIATED( Handle % RtensorValues ) ) THEN
         ALLOCATE( Handle % RtensorValues(maxn1,maxn2,n) )
       END IF
     END IF
          
   END SUBROUTINE ListInitElementKeyword
!------------------------------------------------------------------------------

     
   
!------------------------------------------------------------------------------
!> Given a pointer to the element and the correct handle for the keyword find
!> the list where the keyword valued should be found in. 
!------------------------------------------------------------------------------
   FUNCTION ElementHandleList( Element, Handle, ListSame, ListFound ) RESULT( List )

     TYPE(Element_t), POINTER :: Element     
     TYPE(ValueHandle_t) :: Handle
     TYPE(ValueList_t), POINTER :: List          
     LOGICAL :: ListSame, ListFound
!------------------------------------------------------------------------------     
     INTEGER :: ListId, id
     
     List => NULL()
     
     ListSame = .FALSE.
     ListFound = .FALSE.

      
     ! We are looking for the same element as previous time
     IF( ASSOCIATED( Element, Handle % Element ) ) THEN
       ListSame = .TRUE.
       List => Handle % List
       RETURN
     END IF

     ! Ok, not the same element, get the index that determines the list
     IF( Handle % BulkElement ) THEN     
       ListId = Element % BodyId       
     ELSE
       ListId = 0
       IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
         ListId = Element % BoundaryInfo % Constraint 
       END IF
     END IF
     
     ! We are looking at the same list as previous time
     IF( Handle % ListId == ListId ) THEN
       ListSame = .TRUE.
       List => Handle % List
       RETURN
     ELSE
       Handle % ListId = ListId
       IF( ListId <= 0 ) RETURN
     END IF

     ! Ok, we cannot use previous list, lets find the new list    
     SELECT CASE ( Handle % SectionType )
       
     CASE( SECTION_TYPE_BODY )
       List => CurrentModel % Bodies(ListId) % Values
       ListFound = .TRUE.
       
     CASE( SECTION_TYPE_BF )
       id = ListGetInteger( CurrentModel % Bodies(ListId) % Values, &
           'Body Force', ListFound )         
       IF( ListFound ) List => CurrentModel % BodyForces(id) % Values
       
     CASE( SECTION_TYPE_IC )
       id = ListGetInteger( CurrentModel % Bodies(ListId) % Values, &
           'Initial Condition', ListFound )         
       IF(ListFound) List => CurrentModel % ICs(id) % Values
       
     CASE( SECTION_TYPE_MATERIAL ) 
       id = ListGetInteger( CurrentModel % Bodies(ListId) % Values, &
           'Material', ListFound )         
       IF(ListFound) List => CurrentModel % Materials(id) % Values

     CASE( SECTION_TYPE_EQUATION ) 
       id = ListGetInteger( CurrentModel % Bodies(ListId) % Values, &
           'Equation', ListFound )         
       IF(ListFound) List => CurrentModel % Equations(id) % Values
      
     CASE( SECTION_TYPE_BC )      
       IF( ListId <= 0 .OR. ListId > CurrentModel % NumberOfBCs ) RETURN
       IF( CurrentModel % BCs(ListId) % Tag == ListId ) THEN
         List => CurrentModel % BCs(ListId) % Values
         ListFound = .TRUE.
       END IF
       
     CASE( -1 )
       CALL Fatal('ElementHandleList','Handle not initialized!')

     CASE DEFAULT 
       CALL Fatal('ElementHandleList','Unknown section type!')
       
     END SELECT
     
     IF( ListFound ) THEN
       ! We still have chance that this is the same list
       IF( ASSOCIATED( List, Handle % List ) ) THEN
         ListSame = .TRUE.
       ELSE
         Handle % List => List
       END IF
     ELSE
       Handle % List => NULL()
     END IF          
     
   END FUNCTION ElementHandleList
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Given an index related to the related to the correct section returns the correct
!> value list and a logical flag if there are no more.
!------------------------------------------------------------------------------
   FUNCTION SectionHandleList( Handle, ListId, EndLoop ) RESULT( List )

     TYPE(ValueHandle_t) :: Handle
     TYPE(ValueList_t), POINTER :: List
     INTEGER :: ListId
     LOGICAL :: EndLoop
!------------------------------------------------------------------------------     
     LOGICAL :: Found
     INTEGER :: id
     
     List => NULL()     

     IF( Handle % SectionType == SECTION_TYPE_BC ) THEN            
       EndLoop = ( ListId <= 0 .OR. ListId > CurrentModel % NumberOfBCs )
     ELSE
       EndLoop = ( ListId > CurrentModel % NumberOfBodies )
     END IF       
     IF( EndLoop ) RETURN
     
     
     SELECT CASE ( Handle % SectionType )

     CASE( SECTION_TYPE_BODY )
       List => CurrentModel % Bodies(ListId) % Values

     CASE( SECTION_TYPE_BF )
       id = ListGetInteger( CurrentModel % Bodies(ListId) % Values, &
           'Body Force', Found )         
       IF( Found ) List => CurrentModel % BodyForces(id) % Values

     CASE( SECTION_TYPE_IC )
       id = ListGetInteger( CurrentModel % Bodies(ListId) % Values, &
           'Initial Condition', Found )         
       IF(Found) List => CurrentModel % ICs(id) % Values

     CASE( SECTION_TYPE_MATERIAL ) 
       id = ListGetInteger( CurrentModel % Bodies(ListId) % Values, &
           'Material', Found )         
       IF(Found) List => CurrentModel % Materials(id) % Values

     CASE( SECTION_TYPE_EQUATION ) 
       id = ListGetInteger( CurrentModel % Bodies(ListId) % Values, &
           'Equation',Found )         
       IF(Found) List => CurrentModel % Equations(id) % Values

     CASE( SECTION_TYPE_BC )             
       List => CurrentModel % BCs(ListId) % Values

     CASE( -1 )
       CALL Fatal('SectionHandleList','Handle not initialized!')

     CASE DEFAULT 
       CALL Fatal('SectionHandleList','Unknown section type!')

     END SELECT

   END FUNCTION SectionHandleList
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Compares a string valued parameter in elements and return True if they are the same.
!------------------------------------------------------------------------------
   FUNCTION ListCompareElementAnyString( Handle, RefValue ) RESULT( Same )
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     CHARACTER(LEN=*) :: RefValue     
     LOGICAL :: Same
!------------------------------------------------------------------------------     
     CHARACTER(LEN=MAX_NAME_LEN) :: ThisValue     
     TYPE(ValueList_t), POINTER :: List
     LOGICAL :: Found, EndLoop
     INTEGER :: id, n
!------------------------------------------------------------------------------

     Same = .FALSE.
     
     ! If value is not present anywhere then return False
     IF( Handle % NotPresentAnywhere ) RETURN

     id = 0
     DO WHILE (.TRUE.) 
       id = id + 1
       List => SectionHandleList( Handle, id, EndLoop ) 
       IF( EndLoop ) EXIT
       IF(.NOT. ASSOCIATED( List ) ) CYCLE
       
       ThisValue = ListGetString( List, Handle % Name, Found )
       IF( Found ) THEN         
         n = len_TRIM(ThisValue)
         Same = ( ThisValue(1:n) == RefValue )
         IF( Same ) EXIT
       END IF
     END DO
              
   END FUNCTION ListCompareElementAnyString
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Checks whether any of the logical flags has the desired logical value.
!------------------------------------------------------------------------------
   FUNCTION ListCompareElementAnyLogical( Handle, RefValue ) RESULT( Same )
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     LOGICAL :: RefValue 
     LOGICAL :: Same
!------------------------------------------------------------------------------     
     LOGICAL :: ThisValue
     TYPE(ValueList_t), POINTER :: List
     LOGICAL :: Found, EndLoop
     INTEGER :: id, CValueLen
!------------------------------------------------------------------------------

     Same = .FALSE.
     
     ! If value is not present anywhere then return False
     IF( Handle % NotPresentAnywhere ) RETURN

     id = 0
     DO WHILE (.TRUE.) 
       id = id + 1
       List => SectionHandleList( Handle, id, EndLoop ) 
       IF( EndLoop ) EXIT
       IF(.NOT. ASSOCIATED( List ) ) CYCLE
       
       ThisValue = ListGetLogical( List, Handle % Name, Found )
       IF( Found ) THEN         
         IF( ThisValue .AND. RefValue ) THEN
           Same = .TRUE.
         ELSE IF(.NOT. ThisValue .AND. .NOT. RefValue ) THEN
           Same = .TRUE.
         END IF
         IF( Same ) EXIT
       END IF
     END DO
     
   END FUNCTION ListCompareElementAnyLogical
!------------------------------------------------------------------------------

   
       

!------------------------------------------------------------------------------
!> Get value of parameter from either of the parents.
!> If the value is found then the Left/Right parent is memorized internally.
!> Might not be economical if there are two keywords that toggle but usually
!> we just fetch one keyword from the parents.
!------------------------------------------------------------------------------
  FUNCTION ListGetElementRealParent( Handle, Basis, Element, Found ) RESULT( RValue ) 
     
     TYPE(ValueHandle_t) :: Handle
     TYPE(Element_t), OPTIONAL, POINTER :: Element
     REAL(KIND=dp), OPTIONAL :: Basis(:)
     LOGICAL, OPTIONAL :: Found
     REAL(KIND=dp) :: RValue
     
     LOGICAL :: IntFound
     LOGICAL :: lefttest = .TRUE. ! first start with left test 1st
     TYPE(Element_t), POINTER :: Parent, PElement

     SAVE lefttest


     ! Find the pointer to the element, if not given
     IF( PRESENT( Element ) ) THEN
       PElement => Element
     ELSE
       PElement => CurrentModel % CurrentElement
     END IF

          
     IntFound = .FALSE.
     IF( lefttest) THEN
       Parent => PElement % BoundaryInfo % Left
     ELSE 
       Parent => PElement % BoundaryInfo % Right
     END IF
     RValue = ListGetElementReal( Handle, Basis, Parent, IntFound, PElement % NodeIndexes )

     ! If not found do the same thing with the other parent
     IF(.NOT. IntFound ) THEN
       IF( lefttest) THEN
         Parent => PElement % BoundaryInfo % Right
       ELSE
         Parent => PElement % BoundaryInfo % Left
       END IF
       RValue = ListGetElementReal( Handle, Basis, Parent, IntFound, PElement % NodeIndexes )

       ! reverse the order in which left and right parent are tested
       IF( IntFound ) lefttest = .NOT. lefttest
     END IF
       
     IF( PRESENT( Found ) ) Found = IntFound
     
   END FUNCTION ListGetElementRealParent
     

!------------------------------------------------------------------------------
!> Gets a real valued parameter in the Gaussian integration point defined 
!> by the local basis function. To speed up things there is a handle associated
!> to the given keyword (Name). Here the values are first evaluated at the 
!> nodal points and then using basis functions estimated at the 
!> gaussian integration points. 
!------------------------------------------------------------------------------
   FUNCTION ListGetElementReal( Handle,Basis,Element,Found,Indexes,&
       GaussPoint,Rdim,Rtensor) RESULT(Rvalue)
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     REAL(KIND=dp), OPTIONAL :: Basis(:)
     LOGICAL, OPTIONAL :: Found
     TYPE(Element_t), POINTER, OPTIONAL :: Element
     INTEGER, POINTER, OPTIONAL :: Indexes(:)
     INTEGER, OPTIONAL :: GaussPoint
     INTEGER, OPTIONAL :: Rdim
     REAL(KIND=dp), POINTER, OPTIONAL :: Rtensor(:,:)
     REAL(KIND=dp)  :: Rvalue
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     TYPE(Variable_t), POINTER :: Variable, CVar, TVar
     TYPE(ValueListEntry_t), POINTER :: ptr
     INTEGER, POINTER :: NodeIndexes(:)
     REAL(KIND=dp) :: T(MAX_FNC),x,y,z
!     TYPE(VariableTable_t), SAVE :: VarTable(MAX_FNC)
     REAL(KIND=dp), POINTER :: F(:)
     REAL(KIND=dp), POINTER :: ParF(:,:)
     INTEGER :: i,j,k,j2,k2,k1,l,l0,l1,lsize,n,bodyid,id,n1,n2
     CHARACTER(LEN=MAX_NAME_LEN) ::  cmd, tmp_str
     LOGICAL :: AllGlobal, SomeAtIp, SomeAtNodes, ListSame, ListFound, GotIt, IntFound, &
         ElementSame
     TYPE(Element_t), POINTER :: PElement
#ifdef HAVE_LUA
     INTEGER :: lstat
#endif
!------------------------------------------------------------------------------
     
     ! If value is not present anywhere then return False
     IF( Handle % NotPresentAnywhere ) THEN
       IF(PRESENT(Found)) Found = .FALSE.
       Rvalue = Handle % DefRValue
       RETURN
     END IF

     IF( PRESENT( Rdim ) ) Rdim = 0
     
     ! If the value is known to be globally constant return it asap.
     IF( Handle % ConstantEverywhere ) THEN
       IF(PRESENT(Found)) Found = .TRUE.
       RValue = Handle % RValue
       RETURN
     END IF          

     ! Find the pointer to the element, if not given
     IF( PRESENT( Element ) ) THEN
       PElement => Element
     ELSE
       PElement => CurrentModel % CurrentElement
     END IF

     
     ! Set the default value 
     Rvalue = Handle % DefRValue
     ElementSame = .FALSE.
     
     
     ! We know by initialization the list entry type that the keyword has
     ! Find the correct list to look the keyword in.
     ! Bulk and boundary elements are treated separately.
     List => ElementHandleList( PElement, Handle, ListSame, ListFound ) 

     ! If the provided list is the same as last time, also the keyword will
     ! be sitting at the same place, otherwise find it in the new list
     IF( ListSame ) THEN
       IF( PRESENT( Found ) ) Found = Handle % Found       
       IF( .NOT. Handle % Found ) RETURN

       IF( Handle % GlobalInList ) THEN         
         IF( Handle % Rdim == 0 ) THEN
           Rvalue = Handle % Values(1)
           RETURN
         ELSE
           ! These have been checked already so they should exist
           Rdim = Handle % Rdim
           Rtensor => Handle % Rtensor
           RETURN
         END IF
       ELSE
         ptr => Handle % ptr % head
       END IF
     ELSE IF( ListFound ) THEN

       ptr => ListFind(List,Handle % Name,IntFound )
       IF(PRESENT(Found)) Found = IntFound
       Handle % Found = IntFound
       IF(.NOT. IntFound ) THEN
         IF( Handle % UnfoundFatal ) THEN
           CALL Fatal('ListGetElementReal','Could not find required keyword in list: '//TRIM(Handle % Name))
         END IF
         RETURN
       END IF

       Handle % Ptr % Head => ptr
       Handle % Rdim = ptr % Fdim
       
       IF( Handle % Rdim > 0 ) THEN
         N1 = SIZE(ptr % FValues,1)
         N2 = SIZE(ptr % FValues,2)       
         IF ( ASSOCIATED( Handle % Rtensor) ) THEN
           IF ( SIZE(Handle % Rtensor,1) /= N1 .OR. SIZE(Handle % Rtensor,2) /= N2 ) THEN
             DEALLOCATE( Handle % Rtensor )
           END IF
         END IF
         IF(.NOT. ASSOCIATED( Handle % Rtensor) ) THEN
           ALLOCATE( Handle % Rtensor(N1,N2) )
         END IF

         IF( PRESENT( Rdim ) .AND. PRESENT( Rtensor ) ) THEN
           Rdim = Handle % Rdim
           Rtensor => Handle % Rtensor
         ELSE
           CALL Fatal('ListGetElementReal','For tensors Rdim and Rtensor should be present!')
         END IF
       END IF             
       
       ! It does not make sense to evaluate global variables at IP
       IF( Handle % SomewhereEvaluateAtIp ) THEN
         ! Check whether the keyword should be evaluated at integration point directly                  
         ! Only these dependency type may depend on position
         IF( ptr % TYPE == LIST_TYPE_VARIABLE_SCALAR .OR. &
             ptr % TYPE == LIST_TYPE_VARIABLE_SCALAR_STR .OR.  &
             ptr % TYPE == LIST_TYPE_CONSTANT_SCALAR_PROC ) THEN
           Handle % EvaluateAtIP = ListGetLogical( List, TRIM( Handle % Name )//' At IP',GotIt )           
         ELSE
           Handle % EvaluateAtIp = .FALSE.
         END IF         
       END IF

       IF( Ptr % DepNameLen > 0 ) THEN         
         CALL ListParseStrToVars( Ptr % DependName, Ptr % DepNameLen, &
             Handle % Name, Handle % VarCount, Handle % VarTable, &
             SomeAtIp, SomeAtNodes, AllGlobal )

         Handle % GlobalInList = ( AllGlobal .AND. ptr % PROCEDURE == 0 )
         
         ! If some input parameter is given at integration point
         ! we don't have any option other than evaluate things on IPs
         IF( SomeAtIP ) Handle % EvaluateAtIp = .TRUE.
         Handle % SomeVarAtIp = SomeAtIp 
         
         ! If all variables are global ondes we don't need to evaluate things on IPs
         IF( AllGlobal ) Handle % EvaluateAtIp = .FALSE.

       ELSE
         Handle % GlobalInList = ( ptr % PROCEDURE == 0 )
       END IF
     ELSE
       IF( Handle % UnfoundFatal ) THEN
         CALL Fatal('ListGetElementReal','Could not find list for required keyword: '//TRIM(Handle % Name))
       END IF         
       Rvalue = Handle % DefRValue 
       
       !Handle % Values(1) = RValue
       IF( PRESENT(Found) ) THEN
         Found = .FALSE.
         Handle % Found = .FALSE.
       END IF
       RETURN
     END IF

     
    
     ! Either evaluate parameter directly at IP, 
     ! or first at nodes and then using basis functions at IP.
     ! The latter is the default. 
     !------------------------------------------------------------------
     IF( Handle % EvaluateAtIp ) THEN       
       IF(.NOT. PRESENT(Basis)) THEN
         CALL Fatal('ListGetElementReal','Parameter > Basis < is required!')
       END IF
       
       ! If we get back to the same element than last time use the data already 
       ! retrieved. If the element is new then get the data in every node of the 
       ! current element, or only in the 1st node if it is constant. 
       
       IF( ASSOCIATED( PElement, Handle % Element ) ) THEN
         IF( PRESENT( Indexes ) ) THEN
           n = SIZE( Indexes )
           NodeIndexes => Indexes
         ELSE
           n = Handle % Element % TYPE % NumberOfNodes 
           NodeIndexes => PElement % NodeIndexes
         END IF
           
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
         IF( PRESENT( Indexes ) ) THEN
           n = SIZE( Indexes )
           NodeIndexes => Indexes
         ELSE
           n = PElement % TYPE % NumberOfNodes 
           NodeIndexes => PElement % NodeIndexes
         END IF
           
         IF( ptr % TYPE == LIST_TYPE_VARIABLE_SCALAR .OR. &
             ptr % TYPE == LIST_TYPE_VARIABLE_SCALAR_STR ) THEN

           ! These might not have been initialized if this is has mixed evaluation strategies           
           IF(.NOT. ASSOCIATED( Handle % ParValues )) THEN
             ALLOCATE( Handle % ParValues(MAX_FNC,CurrentModel % Mesh % MaxElementNodes) )
             Handle % ParValues = 0.0_dp
           END IF
           
           DO i=1,n
             k = NodeIndexes(i)

             CALL VarsToValuesOnNodes( Handle % VarCount, Handle % VarTable, k, T, j )
             
             Handle % ParNo = j 
             Handle % ParValues(1:j,i) = T(1:j)

             ! If the dependency table includes just global values (such as time) 
             ! the values will be the same for all element entries.
             IF( Handle % GlobalInList ) EXIT
            
           END DO
         END IF
         ParF => Handle % ParValues         
       END IF


       SELECT CASE(ptr % TYPE)
         
       CASE( LIST_TYPE_VARIABLE_SCALAR )
         
         DO j=1,Handle % ParNo 
           T(j) = SUM( Basis(1:n) *  Handle % ParValues(j,1:n) )
         END DO
         
         ! This one only deals with the variables on IPs, nodal ones are fetched separately
         IF( Handle % SomeVarAtIp ) THEN
           IF( .NOT. PRESENT( GaussPoint ) ) THEN
             CALL Fatal('ListGetElementReal','Evaluation of ip fields requires gauss points as parameter!')
           END IF
           CALL VarsToValuesOnIps( Handle % VarCount, Handle % VarTable, GaussPoint, T, j )
         END IF         
         
         ! there is no node index, pass the negative GaussPoint as to separate it from positive node index
         IF ( ptr % PROCEDURE /= 0 ) THEN
           IF( PRESENT( GaussPoint ) ) THEN
             j = -GaussPoint
           ELSE
             j = 0
           END IF
           !CALL ListPushActiveName(Handle % name)
           Rvalue = ExecRealFunction( ptr % PROCEDURE,CurrentModel, j, T )
           !CALL ListPopActiveName()
         ELSE
           RValue = InterpolateCurve( ptr % TValues,ptr % FValues(1,1,:), &
               T(1), ptr % CubicCoeff )
         END IF
         
       CASE( LIST_TYPE_VARIABLE_SCALAR_STR )

         DO j=1,Handle % ParNo 
           T(j) = SUM( Basis(1:n) *  Handle % ParValues(j,1:n) )
         END DO
         
         ! This one only deals with the variables on IPs, nodal ones have been fecthed already
         IF( Handle % SomeVarAtIp ) THEN
           IF( .NOT. PRESENT( GaussPoint ) ) THEN
             CALL Fatal('ListGetElementReal','Evaluation of ip fields requires gauss points as parameter!')
           END IF
           CALL VarsToValuesOnIps( Handle % VarCount, Handle % VarTable, GaussPoint, T, j )
         END IF
                          
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
         READ( tmp_str(1:k1), * ) RValue
         
       CASE( LIST_TYPE_CONSTANT_SCALAR_PROC )

         IF ( ptr % PROCEDURE /= 0 ) THEN
           x = SUM( Basis(1:n) * CurrentModel % Mesh % Nodes % x( NodeIndexes(1:n) ) )
           y = SUM( Basis(1:n) * CurrentModel % Mesh % Nodes % y( NodeIndexes(1:n) ) )
           z = SUM( Basis(1:n) * CurrentModel % Mesh % Nodes % z( NodeIndexes(1:n) ) )

           !CALL ListPushActiveName(Handle % name)
           RValue = ExecConstRealFunction( ptr % PROCEDURE,CurrentModel,x,y,z)
           !CALL ListPopActiveName()
         ELSE
           CALL Fatal('ListGetElementReal','Constant scalar evaluation failed at ip!')
         END IF
           
       CASE DEFAULT
         
         CALL Fatal('ListGetElementReal','Unknown case for avaluation at ip: '//TRIM(I2S(ptr % Type)))
         
       END SELECT
       
     ELSE ! .NOT. EvaluteAtIp
       
       ! If we get back to the same element than last time use the data already 
       ! retrieved. If the element is new then get the data in every node of the 
       ! current element, or only in the 1st node if it is constant. 
       
       IF( ASSOCIATED( PElement, Handle % Element ) ) THEN
         IF( PRESENT( Indexes ) ) THEN
           n = SIZE( Indexes )
           NodeIndexes => Indexes
         ELSE
           n = Handle % Element % TYPE % NumberOfNodes 
           NodeIndexes => PElement % NodeIndexes
         END IF
         F => Handle % Values       
         ElementSame = .TRUE.
         
       ELSE         
         IF( .NOT. Handle % AllocationsDone ) THEN
           n = CurrentModel % Mesh % MaxElementNodes
           ALLOCATE( Handle % Values(n) )
           Handle % Values = 0.0_dp
           IF( Handle % SomewhereEvaluateAtIp .OR. Handle % EvaluateAtIp ) THEN
             ALLOCATE( Handle % ParValues(MAX_FNC,n) )
             Handle % ParValues = 0.0_dp
           END IF             
           Handle % AllocationsDone = .TRUE.
         END IF
         
         Handle % Element => PElement
         F => Handle % Values

         IF( PRESENT( Indexes ) ) THEN
           n = SIZE( Indexes ) 
           NodeIndexes => Indexes 
         ELSE
           n = PElement % TYPE % NumberOfNodes 
           NodeIndexes => PElement % NodeIndexes
         END IF
           
         SELECT CASE(ptr % TYPE)
           
         CASE( LIST_TYPE_CONSTANT_SCALAR )
           
           IF ( .NOT. ASSOCIATED(ptr % FValues) ) THEN
             WRITE(Message,*) 'Value type for property [', TRIM(Handle % Name), &
                 '] not used consistently.'
             CALL Fatal( 'ListGetElementReal', Message )
             RETURN
           END IF
           F(1) = ptr % Coeff * ptr % Fvalues(1,1,1)


         CASE( LIST_TYPE_VARIABLE_SCALAR )
           !CALL ListPushActiveName(Handle % name)

           T = 1.0_dp
           
           DO i=1,n
             k = NodeIndexes(i)

             CALL VarsToValuesOnNodes( Handle % VarCount, Handle % VarTable, k, T, j )

             IF ( .NOT. ANY( T(1:j) == HUGE(1.0_dp) ) ) THEN
               IF ( ptr % PROCEDURE /= 0 ) THEN
                 F(i) = ptr % Coeff * &
                     ExecRealFunction( ptr % PROCEDURE,CurrentModel, &
                     NodeIndexes(i), T )              
               ELSE
                 IF ( .NOT. ASSOCIATED(ptr % FValues) ) THEN
                   WRITE(Message,*) 'Value type for property [', TRIM(Handle % Name), &
                       '] not used consistently.'
                   CALL Fatal( 'ListGetElementReal', Message )
                   RETURN
                 END IF
                 F(i) = ptr % Coeff * &
                     InterpolateCurve( ptr % TValues,ptr % FValues(1,1,:), &
                     T(1), ptr % CubicCoeff )

                 ! If the dependency table includes just global values (such as time) 
                 ! the values will be the same for all element entries.
                 IF( Handle % GlobalInList ) EXIT

               END IF
             END IF
           END DO
           !CALL ListPopActiveName()
           
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

         CASE( LIST_TYPE_VARIABLE_SCALAR_STR )

           TVar => VariableGet( CurrentModel % Variables, 'Time' ) 
           WRITE( cmd, * ) 'tx=0; st = ', TVar % Values(1)
           k = LEN_TRIM(cmd)
           CALL matc( cmd, tmp_str, k )
           
           DO i=1,n
             k = NodeIndexes(i)
             CALL VarsToValuesOnNodes( Handle % VarCount, Handle % VarTable, k, T, j )
#ifdef HAVE_LUA
             IF ( .NOT. ptr % LuaFun ) THEN
#endif

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
#ifdef HAVE_LUA
             ELSE
               CALL ElmerEvalLua(LuaState, ptr, T, F(i), Handle % varcount)
             END IF
#endif

             IF( Handle % GlobalInList ) EXIT
           END DO

         CASE( LIST_TYPE_CONSTANT_SCALAR_PROC )
           
           IF ( ptr % PROCEDURE == 0 ) THEN
             WRITE(Message,*) 'Value type for property [', TRIM(Handle % Name), &
                 '] not used consistently.'
             CALL Fatal( 'ListGetElementReal', Message )
             RETURN
           END IF
           
           !CALL ListPushActiveName(Handle % name)
           DO i=1,n
             F(i) = ptr % Coeff * &
                 ExecConstRealFunction( ptr % PROCEDURE,CurrentModel, &
                 CurrentModel % Mesh % Nodes % x( NodeIndexes(i) ), &
                 CurrentModel % Mesh % Nodes % y( NodeIndexes(i) ), &
                 CurrentModel % Mesh % Nodes % z( NodeIndexes(i) ) )
           END DO
           !CALL ListPopActiveName()

           
         CASE ( LIST_TYPE_CONSTANT_TENSOR )
           
           n1 = SIZE( Handle % Rtensor, 1 )
           n2 = SIZE( Handle % Rtensor, 2 )
           
           IF ( ptr % PROCEDURE /= 0 ) THEN
             !CALL ListPushActiveName(Handle % name)
             DO i=1,n1
               DO j=1,n2
                 Handle % Rtensor(i,j) = ExecConstRealFunction( ptr % PROCEDURE, &
                     CurrentModel, 0.0_dp, 0.0_dp, 0.0_dp )
               END DO
             END DO
             !CALL ListPopActiveName()
           ELSE
             Handle % Rtensor(:,:) = ptr % FValues(:,:,1)
           END IF
       
           IF( ABS( ptr % Coeff - 1.0_dp ) > EPSILON( ptr % Coeff ) ) THEN
             Handle % Rtensor = ptr % Coeff * Handle % Rtensor
           END IF
           
           
         CASE( LIST_TYPE_VARIABLE_TENSOR,LIST_TYPE_VARIABLE_TENSOR_STR )

           Handle % GlobalInList = .FALSE.
           
           TVar => VariableGet( CurrentModel % Variables, 'Time' ) 
           WRITE( cmd, '(a,e15.8)' ) 'tx=0; st = ', TVar % Values(1)
           k = LEN_TRIM(cmd)
           CALL matc( cmd, tmp_str, k )
           
           !CALL ListPushActiveName(Handle % name)
           
           !CALL ListParseStrToVars( Ptr % DependName, Ptr % DepNameLen, &
           !    Handle % Name, VarCount, VarTable, SomeAtIp, SomeAtNodes, AllGlobal )
           
           IF( PRESENT( Indexes ) ) THEN
             n = SIZE( Indexes )
             NodeIndexes => Indexes
           ELSE
             n = Handle % Element % TYPE % NumberOfNodes 
             NodeIndexes => Handle % Element % NodeIndexes
           END IF

           n1 = SIZE( Handle % Rtensor, 1 )
           n2 = SIZE( Handle % Rtensor, 2 )
           
           DO i=1,n
             k = NodeIndexes(i)
             
             CALL VarsToValuesOnNodes( Handle % VarCount, Handle % VarTable, k, T, j )
             
             IF ( ptr % TYPE==LIST_TYPE_VARIABLE_TENSOR_STR) THEN
#ifdef HAVE_LUA
             IF ( .not. ptr % LuaFun ) THEN
#endif
               DO l=1,j
                 WRITE( cmd, '(a,g19.12)' ) 'tx('//TRIM(i2s(l-1))//')=', T(l)
                 k1 = LEN_TRIM(cmd)
                 CALL matc( cmd, tmp_str, k1 )
               END DO
               
               cmd = ptr % CValue
               k1 = LEN_TRIM(cmd)
               CALL matc( cmd, tmp_str, k1 )
               READ( tmp_str(1:k1), * ) ((Handle % Rtensor(j,k),k=1,N2),j=1,N1)
               
#ifdef HAVE_LUA
             ELSE
               call ElmerEvalLua(LuaState, ptr, T, Handle % RTensor, Handle % varcount)
             END IF
#endif
             ELSE IF ( ptr % PROCEDURE /= 0 ) THEN
               CALL ExecRealArrayFunction( ptr % PROCEDURE, CurrentModel, &
                   NodeIndexes(i), T, Handle % RTensor )
             ELSE
               DO j2=1,N1
                 DO k2=1,N2
                   Handle % Rtensor(j2,k2) = InterpolateCurve(ptr % TValues, ptr % FValues(j2,k2,:), &
                       T(1), ptr % CubicCoeff )
                 END DO
               END DO
             END IF
             
             !CALL ListPopActiveName()
             
             IF( ABS( ptr % Coeff - 1.0_dp ) > EPSILON( ptr % Coeff ) ) THEN
               Handle % Rtensor = ptr % Coeff * Handle % Rtensor
             END IF
             
             IF( Handle % GlobalInList ) THEN
               EXIT              
             ELSE
               DO j2=1,N1
                 DO k2=1,N2               
                   Handle % RtensorValues(j2,k2,i) = Handle % Rtensor(j2,k2)
                 END DO
               END DO
             END IF
               
           END DO
         END SELECT

       END IF

       
       IF( Handle % Rdim == 0 ) THEN
         IF( Handle % GlobalInList ) THEN
           RValue = F(1)
         ELSE
           IF(.NOT. PRESENT(Basis)) THEN
             CALL Fatal('ListGetElementReal','Parameter > Basis < is required!')
           ELSE
             RValue = SUM( Basis(1:n) * F(1:n) )
           END IF
         END IF
       ELSE
         Rtensor => Handle % Rtensor
         Rdim = Handle % Rdim

         IF( .NOT. Handle % GlobalInList ) THEN
           IF(.NOT. PRESENT(Basis)) THEN
             CALL Fatal('ListGetElementRealArray','Parameter > Basis < is required!')
           ELSE
             DO j2=1,SIZE( Handle % RTensor, 1 )
               DO k2=1,SIZE( Handle % RTensor, 2 )               
                 Handle % RTensor(j2,k2) = SUM( Basis(1:n) * Handle % RtensorValues(j2,k2,1:n) )
               END DO
             END DO
           END IF
         END IF
       END IF
       
     END IF
            
     IF ( Handle % GotMinv ) THEN
       IF ( RValue < Handle % minv ) THEN
         WRITE( Message,*) 'Given value ',RValue, ' for property: ', '[', TRIM(Handle % Name),']', &
             ' smaller than given minimum: ', Handle % minv
         CALL Fatal( 'ListGetElementReal', Message )
       END IF
     END IF
       
     IF ( Handle % GotMaxv ) THEN
       IF ( RValue > Handle % maxv ) THEN
         WRITE( Message,*) 'Given value ',RValue, ' for property: ', '[', TRIM(Handle % Name),']', &
             ' larger than given maximum ', Handle % maxv
         CALL Fatal( 'ListGetElementReal', Message )
       END IF
     END IF

   END FUNCTION ListGetElementReal
!------------------------------------------------------------------------------

   
!------------------------------------------------------------------------------
!> This is just a wrapper for getting the imaginary part of the keyword if it
!> has been properly initialized. For the solver modules it is more convenient
!> as the code becomes more compact when using the "HandleIm" field instead of a
!> totally new handle.
!------------------------------------------------------------------------------
   FUNCTION ListGetElementIm( Handle,Basis,Element,Found,Indexes,&
       GaussPoint,Rdim,Rtensor) RESULT(Rvalue)
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     REAL(KIND=dp), OPTIONAL :: Basis(:)
     LOGICAL, OPTIONAL :: Found
     TYPE(Element_t), POINTER, OPTIONAL :: Element
     INTEGER, POINTER, OPTIONAL :: Indexes(:)
     INTEGER, OPTIONAL :: GaussPoint
     INTEGER, OPTIONAL :: Rdim
     REAL(KIND=dp), POINTER, OPTIONAL :: Rtensor(:,:)
     REAL(KIND=dp)  :: Rvalue

     IF(.NOT. ASSOCIATED( Handle % HandleIm ) ) THEN
       CALL Fatal('ListGetElementIm','Initialize with imaginary component!')
     END IF
     Rvalue = ListGetElementReal(Handle % HandleIm,Basis,Element,Found,Indexes,&
         GaussPoint,Rdim,Rtensor)
   END FUNCTION ListGetElementIm
     

!------------------------------------------------------------------------------
!> This is just a wrapper for getting both the real and imaginary part of the keyword if it
!> has been properly initialized. For the solver modules it is convenient since the
!> final code is more compact. This does not work with vector valued keywords yet!
!------------------------------------------------------------------------------
   FUNCTION ListGetElementComplex( Handle,Basis,Element,Found,Indexes,&
       GaussPoint,Rdim,Rtensor) RESULT(Zvalue)
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     REAL(KIND=dp), OPTIONAL :: Basis(:)
     LOGICAL, OPTIONAL :: Found
     TYPE(Element_t), POINTER, OPTIONAL :: Element
     INTEGER, POINTER, OPTIONAL :: Indexes(:)
     INTEGER, OPTIONAL :: GaussPoint
     INTEGER, OPTIONAL :: Rdim
     REAL(KIND=dp), POINTER, OPTIONAL :: Rtensor(:,:)
     COMPLEX(KIND=dp) :: Zvalue

     REAL(KIND=dp) :: RValue, Ivalue
     LOGICAL :: RFound
     
     IF(.NOT. ASSOCIATED( Handle % HandleIm ) ) THEN
       CALL Fatal('ListGetElementComplex','Initialize with imaginary component!')
     END IF

     IF( Handle % NotPresentAnywhere .AND. Handle % HandleIm % NotPresentAnywhere ) THEN
       IF(PRESENT(Found)) Found = .FALSE.
       Zvalue = CMPLX( Handle % DefRValue, 0.0_dp )
       RETURN
     END IF
     
     Rvalue = ListGetElementReal(Handle,Basis,Element,Found,Indexes,GaussPoint)
     IF( PRESENT( Found ) ) RFound = Found 

     Ivalue = ListGetElementReal(Handle % HandleIm,Basis,Element,Found,Indexes,GaussPoint)
     IF( PRESENT( Found ) ) Found = Found .OR. RFound 

     Zvalue = CMPLX( Rvalue, Ivalue ) 
          
   END FUNCTION ListGetElementComplex
       

!------------------------------------------------------------------------------
!> This is just a wrapper for getting a 3D real vector.
!------------------------------------------------------------------------------
   FUNCTION ListGetElementReal3D( Handle,Basis,Element,Found,Indexes,&
       GaussPoint,Rdim,Rtensor) RESULT(RValue3D)
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     REAL(KIND=dp), OPTIONAL :: Basis(:)
     LOGICAL, OPTIONAL :: Found
     TYPE(Element_t), POINTER, OPTIONAL :: Element
     INTEGER, POINTER, OPTIONAL :: Indexes(:)
     INTEGER, OPTIONAL :: GaussPoint
     INTEGER, OPTIONAL :: Rdim
     REAL(KIND=dp), POINTER, OPTIONAL :: Rtensor(:,:)
     REAL(KIND=dp)  :: RValue3D(3)

     LOGICAL :: Found1, Found2
     
     IF(.NOT. ASSOCIATED( Handle % Handle2 ) ) THEN
       CALL Fatal('ListGetElementReal3D','Initialize with 3D components!')
     END IF

     IF( Handle % NotPresentAnywhere .AND. Handle % Handle2 % NotPresentAnywhere &
         .AND.  Handle % Handle3 % NotPresentAnywhere ) THEN
       IF(PRESENT(Found)) Found = .FALSE.
       RValue3D = 0.0_dp
       RETURN
     END IF
     
     Rvalue3D(1) = ListGetElementReal(Handle,Basis,Element,Found,Indexes,GaussPoint)
     IF( PRESENT( Found ) ) Found1 = Found 

     Rvalue3D(2) = ListGetElementReal(Handle % Handle2,Basis,Element,Found,Indexes,GaussPoint)
     IF( PRESENT( Found ) ) Found2 = Found

     Rvalue3D(3) = ListGetElementReal(Handle % Handle3,Basis,Element,Found,Indexes,GaussPoint)
     IF( PRESENT( Found ) ) Found = Found1 .OR. Found2 .OR. Found 
     
   END FUNCTION ListGetElementReal3D


!------------------------------------------------------------------------------
!> This is a wrapper to get gradient of a real valued keyword with functional dependencies.  
!------------------------------------------------------------------------------
   FUNCTION ListGetElementRealGrad( Handle,dBasisdx,Element,Found,Indexes) RESULT(RGrad)
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     ! dBasisdx is required since it is used to evaluate the gradient
     REAL(KIND=dp) :: dBasisdx(:,:)
     LOGICAL, OPTIONAL :: Found
     TYPE(Element_t), POINTER, OPTIONAL :: Element
     INTEGER, POINTER, OPTIONAL :: Indexes(:)
     REAL(KIND=dp)  :: RGrad(3)
     LOGICAL :: Lfound
     INTEGER :: i
     
     RGrad = 0.0_dp
     
     IF( Handle % NotPresentAnywhere ) THEN
       IF( PRESENT( Found ) ) Found = .FALSE.
       RETURN
     END IF

     ! Derivative of constant is zero
     IF( Handle % ConstantEverywhere ) THEN
       IF( PRESENT( Found ) ) Found = .TRUE.      
       RETURN
     END IF

     ! Obtain gradient of a scalar field going through the partial derivatives of the components
     DO i=1,3     
       RGrad(i) = ListGetElementReal(Handle,dBasisdx(:,i),Element,Lfound,Indexes)
       ! If we don't have it needless to contunue to 2nd and 3rd dimensions
       IF(.NOT. Lfound ) EXIT
     END DO
     IF( PRESENT( Found ) ) Found = Lfound
     
   END FUNCTION ListGetElementRealGrad


!------------------------------------------------------------------------------
!> This is just a wrapper for getting divergence of a 3D real vector.
!------------------------------------------------------------------------------
   FUNCTION ListGetElementRealDiv( Handle,dBasisdx,Element,Found,Indexes) RESULT(Rdiv)
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     ! dBasisdx is required since it is used to evaluate the divergence
     REAL(KIND=dp) :: dBasisdx(:,:)
     LOGICAL, OPTIONAL :: Found
     TYPE(Element_t), POINTER, OPTIONAL :: Element
     INTEGER, POINTER, OPTIONAL :: Indexes(:)
     REAL(KIND=dp)  :: Rdiv(3)

     LOGICAL :: Found1, Found2, Found3

     Rdiv = 0.0_dp
     
     IF(.NOT. ASSOCIATED( Handle % Handle2 ) ) THEN
       CALL Fatal('ListGetElementReal3D','Initialize with 3D components!')
     END IF

     IF( Handle % NotPresentAnywhere .AND. Handle % Handle2 % NotPresentAnywhere &
         .AND.  Handle % Handle3 % NotPresentAnywhere ) THEN
       IF(PRESENT(Found)) Found = .FALSE.
       RETURN
     END IF

     Rdiv(1) = ListGetElementReal(Handle,dBasisdx(:,1),Element,Found1,Indexes)
     Rdiv(2) = ListGetElementReal(Handle % Handle2,dBasisdx(:,2),Element,Found2,Indexes)
     Rdiv(3) = ListGetElementReal(Handle % Handle3,dBasisdx(:,3),Element,Found3,Indexes)
     IF( PRESENT( Found ) ) Found = Found1 .OR. Found2 .OR. Found3
     
   END FUNCTION ListGetElementRealDiv


   
!------------------------------------------------------------------------------
!> This is just a wrapper for getting a 3D complex vector.
!------------------------------------------------------------------------------
   FUNCTION ListGetElementComplex3D( Handle,Basis,Element,Found,Indexes,&
       GaussPoint,Rdim,Rtensor) RESULT(ZValue3D)
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     REAL(KIND=dp), OPTIONAL :: Basis(:)
     LOGICAL, OPTIONAL :: Found
     TYPE(Element_t), POINTER, OPTIONAL :: Element
     INTEGER, POINTER, OPTIONAL :: Indexes(:)
     INTEGER, OPTIONAL :: GaussPoint
     INTEGER, OPTIONAL :: Rdim
     REAL(KIND=dp), POINTER, OPTIONAL :: Rtensor(:,:)
     COMPLEX(KIND=dp)  :: ZValue3D(3)

     REAL(KIND=dp)  :: RValue3D(3), IValue3D(3)
     LOGICAL :: RFound
     
     IF(.NOT. ASSOCIATED( Handle % HandleIm ) ) THEN
       CALL Fatal('ListGetElementComplex3D','Initialize with imaginary component!')
     END IF
     
     Rvalue3D = ListGetElementReal3D(Handle,Basis,Element,Found,Indexes,GaussPoint)
     IF( PRESENT( Found ) ) RFound = Found 
     
     Ivalue3D = ListGetElementReal3D(Handle % HandleIm,Basis,Element,Found,Indexes,GaussPoint)
     IF( PRESENT( Found ) ) Found = Found .OR. RFound
     
     Zvalue3D = CMPLX( Rvalue3D, Ivalue3D )     
     
   END FUNCTION ListGetElementComplex3D

   
!------------------------------------------------------------------------------
!> Gets a real valued parameter in all the Gaussian integration points.
!------------------------------------------------------------------------------
   FUNCTION ListGetElementRealVec( Handle,ngp,BasisVec,Element,Found ) RESULT( Rvalues )
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     INTEGER :: ngp
     REAL(KIND=dp), OPTIONAL :: BasisVec(:,:)
     LOGICAL, OPTIONAL :: Found
     TYPE(Element_t), POINTER, OPTIONAL :: Element
     REAL(KIND=dp), POINTER  :: Rvalues(:)
!------------------------------------------------------------------------------
     TYPE(Variable_t), POINTER :: Variable, CVar, TVar
     TYPE(ValueListEntry_t), POINTER :: ptr
     INTEGER, POINTER :: NodeIndexes(:)
     REAL(KIND=dp) :: T(MAX_FNC),x,y,z, RValue
!     TYPE(VariableTable_t) :: VarTable(MAX_FNC)
     REAL(KIND=dp), POINTER :: F(:)
     REAL(KIND=dp), POINTER :: ParF(:,:)
     INTEGER :: i,j,k,k1,l,l0,l1,lsize,n,bodyid,id,node,gp
     !,varcount
     CHARACTER(LEN=MAX_NAME_LEN) :: cmd, tmp_str
     LOGICAL :: AllGlobal, SomeAtIp, SomeAtNodes, ListSame, ListFound, GotIt, IntFound
     TYPE(Element_t), POINTER :: PElement
     TYPE(ValueList_t), POINTER :: List
!------------------------------------------------------------------------------

     IF( Handle % nValuesVec < ngp ) THEN
       IF( Handle % nValuesVec > 0 ) THEN
         DEALLOCATE( Handle % ValuesVec )
       END IF
       ALLOCATE( Handle % ValuesVec(ngp) )
       Handle % nValuesVec = ngp       

       IF( Handle % ConstantEverywhere ) THEN
         Handle % ValuesVec(1:ngp) = Handle % Rvalue
       ELSE
         Handle % ValuesVec(1:ngp) = Handle % DefRValue        
       END IF
     END IF

     ! The results are always returned from the Handle % Values
     Rvalues => Handle % ValuesVec

     ! If value is not present anywhere then return False
     IF( Handle % NotPresentAnywhere ) THEN
       IF(PRESENT(Found)) Found = .FALSE.
       RETURN
     END IF

     ! If the value is known to be globally constant return it asap.
     IF( Handle % ConstantEverywhere ) THEN
       IF(PRESENT(Found)) Found = .TRUE.
       RETURN
     END IF     

     ! Find the pointer to the element, if not given
     IF( PRESENT( Element ) ) THEN
       PElement => Element
     ELSE
       PElement => CurrentModel % CurrentElement
     END IF

     ! We know by initialization the list entry type that the keyword has
     ! Find the correct list to look the keyword in.
     ! Bulk and boundary elements are treated separately.
     List => ElementHandleList( PElement, Handle, ListSame, ListFound ) 

     
     ! If the provided list is the same as last time, also the keyword will
     ! be sitting at the same place, otherwise find it in the new list
     IF( ListSame ) THEN
       IF( PRESENT( Found ) ) Found = Handle % Found       
       IF( .NOT. Handle % Found ) RETURN
       IF( Handle % GlobalInList ) THEN
         RETURN
       ELSE
         ptr => Handle % ptr % head        
       END IF
     ELSE IF( ListFound ) THEN

       ptr => ListFind(List,Handle % Name,IntFound)
       IF(PRESENT(Found)) Found = IntFound
       Handle % Found = IntFound
       
       IF(.NOT. IntFound ) THEN
         IF( Handle % UnfoundFatal ) THEN
           CALL Fatal('ListGetElementRealVec','Could not find required keyword in list: '//TRIM(Handle % Name))
         END IF
         Handle % ValuesVec(1:ngp) = Handle % DefRValue
         RETURN
       END IF
         
       Handle % Ptr % Head => ptr
       
       ! It does not make sense to evaluate global variables at IP
       IF( Handle % SomewhereEvaluateAtIp ) THEN
         ! Check whether the keyword should be evaluated at integration point directly                  
         ! Only these dependency type may depend on position
         IF( ptr % TYPE == LIST_TYPE_VARIABLE_SCALAR .OR. &
             ptr % TYPE == LIST_TYPE_VARIABLE_SCALAR_STR .OR.  &
             ptr % TYPE == LIST_TYPE_CONSTANT_SCALAR_PROC ) THEN
         ! Check whether the keyword should be evaluated at integration point directly
           Handle % EvaluateAtIp = ListGetLogical( List, TRIM( Handle % Name )//' At IP',GotIt )
         ELSE
           Handle % EvaluateAtIp = .FALSE.
         END IF
       END IF

       
       IF( ptr % DepNameLen > 0 ) THEN
         CALL ListParseStrToVars( Ptr % DependName, Ptr % DepNameLen, &
             Handle % Name, Handle % VarCount, Handle % VarTable, &
             SomeAtIp, SomeAtNodes, AllGlobal )
         IF( SomeAtIp ) Handle % EvaluateAtIp = .TRUE.
         Handle % GlobalInList = ( AllGlobal .AND. ptr % PROCEDURE == 0 )
         IF( SomeAtIP ) Handle % EvaluateAtIp = .TRUE.
         IF( AllGlobal ) Handle % EvaluateAtIp = .FALSE.
         Handle % SomeVarAtIp = SomeAtIp 
       ELSE
         Handle % GlobalInList = ( ptr % PROCEDURE == 0 )
       END IF

     ELSE
       IF( Handle % UnfoundFatal ) THEN
         CALL Fatal('ListGetElementRealVec','Could not find list for required keyword: '//TRIM(Handle % Name))
       END IF                
       IF( .NOT. Handle % AllocationsDone ) THEN
         n = CurrentModel % Mesh % MaxElementNodes
         ALLOCATE( Handle % Values(n) )
         Handle % Values = 0.0_dp
         IF( Handle % SomewhereEvaluateAtIp .OR. Handle % EvaluateAtIp ) THEN
           ALLOCATE( Handle % ParValues(MAX_FNC,n) )
           Handle % ParValues = 0.0_dp
         END IF
         Handle % AllocationsDone = .TRUE.
       END IF
       Handle % ValuesVec = Handle % DefRValue
       IF( PRESENT(Found) ) THEN
         Found = .FALSE.
         Handle % Found = .FALSE.
       END IF
       RETURN
     END IF

     ! Either evaluate parameter directly at IP, 
     ! or first at nodes and then using basis functions at IP.
     ! The later is the default. 
     !------------------------------------------------------------------
     IF( Handle % EvaluateAtIp ) THEN

       IF(.NOT. PRESENT(BasisVec)) THEN
         CALL Fatal('ListGetElementRealVec','Parameter > Basis < is required!')
       END IF

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
           node = NodeIndexes(i)
           CALL VarsToValuesOnNodes( Handle % VarCount, Handle % VarTable, node, T, j )

           IF( Handle % GlobalInList ) THEN
             CALL Warn('ListGetElementRealVec','Constant expression need not be evaluated at IPs!')
           END IF

           Handle % ParNo = j 
           Handle % ParValues(1:j,i) = T(1:j)
         END DO

         ParF => Handle % ParValues         
       END IF


       SELECT CASE(ptr % TYPE)

       CASE( LIST_TYPE_VARIABLE_SCALAR )

         ! there is no node index, so use zero
         IF ( ptr % PROCEDURE /= 0 ) THEN
           !CALL ListPushActiveName(Handle % name)
           node = 0 

           DO gp = 1, ngp          
             DO j=1,Handle % ParNo 
               T(j) = SUM( BasisVec(gp,1:n) *  ParF(j,1:n) )
             END DO
             Rvalue = ExecRealFunction( ptr % PROCEDURE, CurrentModel, node, T )
             Handle % ValuesVec(gp) = RValue
           END DO
           !CALL ListPopActiveName()
         ELSE
           DO gp = 1, ngp          
             DO j=1,Handle % ParNo 
               T(j) = SUM( BasisVec(gp,1:n) *  ParF(j,1:n) )
             END DO
             RValue = InterpolateCurve( ptr % TValues,ptr % FValues(1,1,:), &
                 T(1), ptr % CubicCoeff )
             Handle % ValuesVec(gp) = RValue
           END DO
         END IF

       CASE( LIST_TYPE_VARIABLE_SCALAR_STR )

         ! there is no node index, so use zero
         node = 0 

#ifdef HAVE_LUA
         IF ( .not. ptr % LuaFun ) THEN
#endif
           TVar => VariableGet( CurrentModel % Variables, 'Time' ) 
           WRITE( cmd, * ) 'tx=0; st = ', TVar % Values(1)
           k = LEN_TRIM(cmd)
           CALL matc( cmd, tmp_str, k )         
#ifdef HAVE_LUA
         END IF
#endif
         
         DO gp = 1, ngp          
           DO j=1,Handle % ParNo 
             T(j) = SUM( BasisVec(gp,1:n) *  Handle % ParValues(j,1:n) )
           END DO

           ! This one only deals with the variables on IPs, nodal ones have been fecthed already
           IF( Handle % SomeVarAtIp ) THEN
             CALL VarsToValuesOnIps( Handle % VarCount, Handle % VarTable, gp, T, j )
           END IF

#ifdef HAVE_LUA
           IF ( .not. ptr % LuaFun ) THEN
#endif
             DO l=1,Handle % ParNo
               WRITE( cmd, * ) 'tx('//TRIM(i2s(l-1))//')=', T(l)
               k1 = LEN_TRIM(cmd)
               CALL matc( cmd, tmp_str, k1 )
             END DO

             cmd = ptr % CValue
             k1 = LEN_TRIM(cmd)

             CALL matc( cmd, tmp_str, k1 )
             READ( tmp_str(1:k1), * ) RValue

#ifdef HAVE_LUA
           ELSE
             call ElmerEvalLua(LuaState, ptr, T, RValue, j)
           END IF
#endif
           Handle % ValuesVec(gp) = RValue
         END DO


       CASE( LIST_TYPE_CONSTANT_SCALAR_PROC )

         IF ( ptr % PROCEDURE /= 0 ) THEN
           !CALL ListPushActiveName(Handle % name)

           DO gp = 1, ngp          

             x = SUM(BasisVec(gp,1:n) * CurrentModel % Mesh % Nodes % x( NodeIndexes(1:n)))
             y = SUM(BasisVec(gp,1:n) * CurrentModel % Mesh % Nodes % y( NodeIndexes(1:n)))
             z = SUM(BasisVec(gp,1:n) * CurrentModel % Mesh % Nodes % z( NodeIndexes(1:n)))

             RValue = ExecConstRealFunction( ptr % PROCEDURE,CurrentModel,x,y,z)
             Handle % ValuesVec(gp) = RValue
           END DO
           !CALL ListPopActiveName()

         ELSE
           CALL Fatal('ListGetElementRealVec','Constant scalar evaluation failed at ip!')
         END IF

       CASE DEFAULT

         CALL Fatal('ListGetElementRealVec','Unknown case for avaluation at ip')

       END SELECT

     ELSE

       IF( .NOT. Handle % AllocationsDone ) THEN
         n = CurrentModel % Mesh % MaxElementNodes
         ALLOCATE( Handle % Values(n) )
         Handle % Values = 0.0_dp
         IF( Handle % SomewhereEvaluateAtIp .OR. Handle % EvaluateAtIp ) THEN
           ALLOCATE( Handle % ParValues(MAX_FNC,n) )
           Handle % ParValues = 0.0_dp
         END IF
         Handle % AllocationsDone = .TRUE.
       END IF

       Handle % Element => PElement
       n = PElement % TYPE % NumberOfNodes 
       NodeIndexes => PElement % NodeIndexes
       F => Handle % Values

       SELECT CASE(ptr % TYPE)

       CASE( LIST_TYPE_CONSTANT_SCALAR )

         IF ( .NOT. ASSOCIATED(ptr % FValues) ) THEN
           WRITE(Message,*) 'Value type for property [', TRIM(Handle % Name), &
               '] not used consistently.'
           CALL Fatal( 'ListGetElementRealVec', Message )
           RETURN
         END IF
         F(1) = ptr % Coeff * ptr % Fvalues(1,1,1)
         RValues(1:ngp) = F(1)


       CASE( LIST_TYPE_VARIABLE_SCALAR )

         !CALL ListPushActiveName(Handle % name)

         DO i=1,n
           node = NodeIndexes(i)
           CALL VarsToValuesOnNodes( Handle % VarCount, Handle % VarTable, node, T, j )
           
           IF ( ptr % PROCEDURE /= 0 ) THEN
             F(i) = ptr % Coeff * &
                 ExecRealFunction( ptr % PROCEDURE,CurrentModel, &
                 NodeIndexes(i), T )              
           ELSE
             IF ( .NOT. ASSOCIATED(ptr % FValues) ) THEN
               WRITE(Message,*) 'Value type for property [', TRIM(Handle % Name), &
                   '] not used consistently.'
               CALL Fatal( 'ListGetElementRealVec', Message )
               RETURN
             END IF
             F(i) = ptr % Coeff * &
                 InterpolateCurve( ptr % TValues,ptr % FValues(1,1,:), &
                 T(1), ptr % CubicCoeff )

             ! If the dependency table includes just global values (such as time) 
             ! the values will be the same for all element entries.
             IF( Handle % GlobalInList ) EXIT
           END IF
         END DO
         
         IF( Handle % GlobalInList ) THEN
           Handle % ValuesVec(1:ngp) = F(1)
         ELSE
           Handle % ValuesVec(1:ngp) = MATMUL( BasisVec(1:ngp,1:n), F(1:n) )
         END IF
         !CALL ListPopActiveName()


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

         Handle % ValuesVec(1:ngp) = F(1)


       CASE( LIST_TYPE_VARIABLE_SCALAR_STR )

#ifdef HAVE_LUA
         IF ( .not. ptr % LuaFun ) THEN
#endif
           TVar => VariableGet( CurrentModel % Variables, 'Time' ) 
           WRITE( cmd, * ) 'tx=0; st = ', TVar % Values(1)
           k = LEN_TRIM(cmd)
           CALL matc( cmd, tmp_str, k )
#ifdef HAVE_LUA
         END IF
#endif

         DO i=1,n
           k = NodeIndexes(i)
           
           CALL VarsToValuesOnNodes( Handle % VarCount, Handle % VarTable, k, T, j )

#ifdef HAVE_LUA
           IF ( .not. ptr % LuaFun ) THEN
#endif
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
#ifdef HAVE_LUA
           ELSE
             call ElmerEvalLuaS(LuaState, ptr, T, F(i), j)
             F(i) = ptr % coeff * F(i)
           END IF
#endif
           IF( Handle % GlobalInList ) EXIT
         END DO

         IF( Handle % GlobalInList ) THEN
           Handle % ValuesVec(1:ngp) = F(1)
         ELSE
           Handle % ValuesVec(1:ngp) = MATMUL( BasisVec(1:ngp,1:n), F(1:n) )
         END IF

       CASE( LIST_TYPE_CONSTANT_SCALAR_PROC )
         IF ( ptr % PROCEDURE == 0 ) THEN
           WRITE(Message,*) 'Value type for property [', TRIM(Handle % Name), &
               '] not used consistently.'
           CALL Fatal( 'ListGetElementRealVec', Message )
           RETURN
         END IF

         !CALL ListPushActiveName(Handle % name)
         DO i=1,n
           F(i) = ptr % Coeff * &
               ExecConstRealFunction( ptr % PROCEDURE,CurrentModel, &
               CurrentModel % Mesh % Nodes % x( NodeIndexes(i) ), &
               CurrentModel % Mesh % Nodes % y( NodeIndexes(i) ), &
               CurrentModel % Mesh % Nodes % z( NodeIndexes(i) ) )
         END DO
         !CALL ListPopActiveName()

         Handle % ValuesVec(1:ngp) = MATMUL( BasisVec(1:ngp,1:n), F(1:n) )

       CASE DEFAULT
         CALL Fatal('ListGetElementRealVec','Impossible entry type: '//TRIM(I2S(ptr % Type)))
         
       END SELECT

     END IF
     
   END FUNCTION ListGetElementRealVec
!------------------------------------------------------------------------------



   
!------------------------------------------------------------------------------
!> Gets a logical valued parameter in elements.
!------------------------------------------------------------------------------
   FUNCTION ListGetElementLogical( Handle, Element, Found ) RESULT(Lvalue)
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     TYPE(Element_t), POINTER, OPTIONAL :: Element
     LOGICAL, OPTIONAL :: Found
     LOGICAL  :: Lvalue
!------------------------------------------------------------------------------     
     TYPE(ValueList_t), POINTER :: List
     TYPE(Element_t), POINTER :: PElement
     LOGICAL :: ListSame, ListFound
     INTEGER :: id, BodyId
!------------------------------------------------------------------------------

     ! If value is not present anywhere then return False
     IF( Handle % NotPresentAnywhere ) THEN
       IF(PRESENT(Found)) Found = .FALSE.
       Lvalue = Handle % DefLValue
       RETURN
     END IF
     
     ! If the value is known to be globally constant return it asap.
     IF( Handle % ConstantEverywhere ) THEN
       IF(PRESENT(Found)) Found = .TRUE.
       Lvalue = Handle % LValue
       RETURN
     END IF

     ! Find the pointer to the element, if not given
     IF( PRESENT( Element ) ) THEN
       PElement => Element
     ELSE
       PElement => CurrentModel % CurrentElement
     END IF
     
     ! We know by initialization the list entry type that the keyword has
     ! Find the correct list to look the keyword in.
     ! Bulk and boundary elements are treated separately.
     List => ElementHandleList( PElement, Handle, ListSame, ListFound ) 
     
     IF( ListSame ) THEN
       IF( PRESENT( Found ) ) Found = Handle % Found 
       LValue = Handle % LValue
     ELSE IF( ListFound ) THEN
       LValue = ListGetLogical( List, Handle % Name, Found, &
           UnfoundFatal = Handle % UnfoundFatal  )
       IF( .NOT. Found ) Lvalue = Handle % DefLValue
       Handle % LValue = LValue       
       IF(PRESENT(Found)) Handle % Found = Found 
     ELSE     
       IF( Handle % UnfoundFatal ) THEN
         CALL Fatal('ListGetElementLogical','Could not find list for required keyword: '//TRIM(Handle % Name))
       END IF         
       Lvalue = Handle % DefLValue 
       IF( PRESENT(Found) ) THEN
         Found = .FALSE.
         Handle % Found = .FALSE.
       END IF
     END IF
                   
   END FUNCTION ListGetElementLogical
!------------------------------------------------------------------------------

   
!------------------------------------------------------------------------------
!> Gets a integer valued parameter in elements.
!------------------------------------------------------------------------------
   FUNCTION ListGetElementInteger( Handle, Element, Found ) RESULT(Ivalue)
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     TYPE(Element_t), POINTER, OPTIONAL :: Element
     LOGICAL, OPTIONAL :: Found
     INTEGER  :: Ivalue
!------------------------------------------------------------------------------     
     TYPE(ValueList_t), POINTER :: List
     TYPE(Element_t), POINTER :: PElement
     LOGICAL :: ListSame, ListFound
     INTEGER :: id, BodyId
!------------------------------------------------------------------------------

     ! If value is not present anywhere then return False
     IF( Handle % NotPresentAnywhere ) THEN
       IF(PRESENT(Found)) Found = .FALSE.
       Ivalue = Handle % DefIValue
       RETURN
     END IF
     
     ! If the value is known to be globally constant return it asap.
     IF( Handle % ConstantEverywhere ) THEN
       IF(PRESENT(Found)) Found = .TRUE.
       Ivalue = Handle % IValue
       RETURN
     END IF

     ! Find the pointer to the element, if not given
     IF( PRESENT( Element ) ) THEN
       PElement => Element
     ELSE
       PElement => CurrentModel % CurrentElement
     END IF
     
     ! We know by initialization the list entry type that the keyword has
     ! Find the correct list to look the keyword in.
     ! Bulk and boundary elements are treated separately.
     List => ElementHandleList( PElement, Handle, ListSame, ListFound ) 

     IF( ListSame ) THEN
       IF( PRESENT( Found ) ) Found = Handle % Found 
       IValue = Handle % IValue
     ELSE IF( ListFound ) THEN
       IValue = ListGetInteger( List, Handle % Name, Found, UnfoundFatal = Handle % UnfoundFatal )
       Handle % IValue = IValue
       IF(PRESENT(Found)) Handle % Found = Found 
     ELSE     
       IF( Handle % UnfoundFatal ) THEN
         CALL Fatal('ListGetElementInteger','Could not find list for required keyword: '//TRIM(Handle % Name))
       END IF         
       Ivalue = Handle % DefIValue
       Handle % IValue = IValue
       IF( PRESENT(Found) ) THEN
         Found = .FALSE.
         Handle % Found = .FALSE.
       END IF
     END IF
              
     
   END FUNCTION ListGetElementInteger
!------------------------------------------------------------------------------


   
!------------------------------------------------------------------------------
!> Gets a string valued parameter in elements.
!------------------------------------------------------------------------------
   FUNCTION ListGetElementString( Handle, Element, Found ) RESULT( CValue )
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     CHARACTER(LEN=MAX_NAME_LEN) :: CValue     
     TYPE(Element_t), POINTER, OPTIONAL :: Element
     LOGICAL, OPTIONAL :: Found
!------------------------------------------------------------------------------     
     TYPE(ValueList_t), POINTER :: List
     TYPE(Element_t), POINTER :: PElement
     LOGICAL :: ListSame, ListFound
     INTEGER :: id, BodyId
!------------------------------------------------------------------------------

     ! If value is not present anywhere then return False
     IF( Handle % NotPresentAnywhere ) THEN
       IF(PRESENT(Found)) Found = .FALSE.
       Cvalue = ' '
       RETURN
     END IF
     
     ! If the value is known to be globally constant return it asap.
     IF( Handle % ConstantEverywhere ) THEN
       IF(PRESENT(Found)) Found = .TRUE.
       Cvalue = TRIM(Handle % CValue)
       RETURN
     END IF

     ! Find the pointer to the element, if not given
     IF( PRESENT( Element ) ) THEN
       PElement => Element
     ELSE
       PElement => CurrentModel % CurrentElement
     END IF
     
     ! We know by initialization the list entry type that the keyword has
     ! Find the correct list to look the keyword in.
     ! Bulk and boundary elements are treated separately.
     List => ElementHandleList( PElement, Handle, ListSame, ListFound ) 

     IF( ListSame ) THEN
       IF( PRESENT( Found ) ) Found = Handle % Found 
       CValue = Handle % CValue(1:Handle % CValueLen)
     ELSE IF( ListFound ) THEN
       CValue = ListGetString( List, Handle % Name, Found, &
           UnfoundFatal = Handle % UnfoundFatal )
       Handle % CValue = TRIM(CValue)
       Handle % CValueLen = len_trim(CValue)
       IF(PRESENT(Found)) Handle % Found = Found 
     ELSE     
       IF( Handle % UnfoundFatal ) THEN
         CALL Fatal('ListGetElementString','Could not find list for required keyword: '//TRIM(Handle % Name))
       END IF
       Cvalue = ' '
       Handle % CValueLen = 0
       IF( PRESENT(Found) ) THEN
         Found = .FALSE.
         Handle % Found = .FALSE.
       END IF
     END IF
              
   END FUNCTION ListGetElementString
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Is the keyword present somewhere
!------------------------------------------------------------------------------
   FUNCTION ListGetElementSomewhere( Handle ) RESULT( Found )
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     LOGICAL :: Found 
!------------------------------------------------------------------------------     
     Found = .NOT. ( Handle % NotPresentAnywhere )

   END FUNCTION ListGetElementSomewhere
!------------------------------------------------------------------------------     


   

!------------------------------------------------------------------------------
!> Compares a string valued parameter in elements and return True if they are the same.
!------------------------------------------------------------------------------
   FUNCTION ListCompareElementString( Handle, CValue2, Element, Found ) RESULT( SameString )
!------------------------------------------------------------------------------
     TYPE(ValueHandle_t) :: Handle
     CHARACTER(LEN=*) :: CValue2     
     TYPE(Element_t), POINTER, OPTIONAL :: Element
     LOGICAL, OPTIONAL :: Found
     LOGICAL :: SameString
!------------------------------------------------------------------------------     
     CHARACTER(LEN=MAX_NAME_LEN) :: CValue     
     TYPE(ValueList_t), POINTER :: List
     TYPE(Element_t), POINTER :: PElement
     LOGICAL :: ListSame, ListFound, IntFound
     INTEGER :: id, BodyId
!------------------------------------------------------------------------------

     SameString = .FALSE.
     
     ! If value is not present anywhere then return False
     IF( Handle % NotPresentAnywhere ) THEN
       IF(PRESENT(Found)) Found = .FALSE.
       RETURN
     END IF
     
     ! If the value is known to be globally constant return it asap.
     IF( Handle % ConstantEverywhere ) THEN
       IF(PRESENT(Found)) Found = .TRUE.
       SameString = ( CValue2 == Handle % CValue(1:Handle % CValueLen) )
       RETURN
     END IF
     
     ! Find the pointer to the element, if not given
     IF( PRESENT( Element ) ) THEN
       PElement => Element
     ELSE
       PElement => CurrentModel % CurrentElement
     END IF

     ListSame = .FALSE.
     ListFound = .FALSE.
     
     ! We know by initialization the list entry type that the keyword has
     ! Find the correct list to look the keyword in.
     ! Bulk and boundary elements are treated separately.
     List => ElementHandleList( PElement, Handle, ListSame, ListFound ) 

     IF( ListSame ) THEN
       IF( PRESENT( Found ) ) Found = Handle % Found 
       IF( Handle % Found ) THEN
         SameString = ( Handle % CValue(1:Handle % CValueLen) == CValue2 )
       END IF
     ELSE IF( ListFound ) THEN
       CValue = ListGetString( List, Handle % Name, IntFound, &
           UnfoundFatal = Handle % UnfoundFatal )
       Handle % Found = IntFound
       IF( IntFound ) THEN
         Handle % CValueLen = len_trim(CValue)
         Handle % CValue = CValue(1:Handle % CValueLen )
         SameString = (Handle % CValue(1:Handle % CValueLen) == CValue2 )
       END IF
       IF(PRESENT(Found)) Found = IntFound 
     ELSE     
       Handle % Cvalue = ' '
       Handle % CValueLen = 0
       Handle % Found = .FALSE.
       IF( PRESENT(Found) ) Found = .FALSE.
     END IF
              
   END FUNCTION ListCompareElementString
!------------------------------------------------------------------------------


     
!------------------------------------------------------------------------------
!> Initializes the variable handle in a similar manner as the keyword handle is
!> initialized. This handle is more compact. Does not support p-fields or
!> Hcurl & Hdiv fields yet. 
!------------------------------------------------------------------------------
   SUBROUTINE ListInitElementVariable( Handle, Name, USolver, UVariable, tStep )
!------------------------------------------------------------------------------
     TYPE(VariableHandle_t) :: Handle
     CHARACTER(LEN=*), OPTIONAL  :: Name
     TYPE(Solver_t), OPTIONAL, TARGET :: USolver
     TYPE(Variable_t), OPTIONAL, TARGET :: UVariable
     INTEGER, OPTIONAL :: tStep

     REAL(KIND=dp), POINTER :: Values(:)
     TYPE(Variable_t), POINTER :: Variable
     TYPE(Solver_t)  , POINTER :: Solver
     TYPE(Element_t),  POINTER :: Element

     Handle % Variable => NULL() 
     Handle % Values => NULL()
     Handle % Perm => NULL()
     Handle % Element => NULL()
          
     IF ( PRESENT(USolver) ) THEN
       Solver => USolver
     ELSE
       Solver => CurrentModel % Solver
     END IF
            
     IF ( PRESENT(name) ) THEN
       Variable => VariableGet( Solver % Mesh % Variables, name )
     ELSE IF( PRESENT( UVariable ) ) THEN
       Variable => UVariable
     ELSE
       Variable => Solver % Variable 
     END IF
     
     IF ( .NOT. ASSOCIATED( Variable ) ) RETURN
     
     Handle % Variable => Variable
     
     IF ( PRESENT(tStep) ) THEN
       IF ( tStep < 0 ) THEN
         IF ( ASSOCIATED(Variable % PrevValues) .AND. -tStep<=SIZE(Variable % PrevValues,2)) &
             Handle % Values => Variable % PrevValues(:,-tStep)
       END IF
     ELSE
       Handle % Values => Variable % Values      
     END IF
     Handle % Perm => Variable % Perm
     
   END SUBROUTINE ListInitElementVariable
!------------------------------------------------------------------------------

     
!------------------------------------------------------------------------------
!> Get a scalar field (e.g. potential or pressure) at the integration point.
!> Works with different types of fields.
!------------------------------------------------------------------------------
   FUNCTION ListGetElementScalarSolution( Handle, Basis, Element, Found, GaussPoint, nd  ) RESULT ( Val )
     
     TYPE(VariableHandle_t) :: Handle
     REAL(KIND=dp), OPTIONAL :: Basis(:)
     TYPE( Element_t), POINTER, OPTIONAL :: Element
     INTEGER, OPTIONAL :: nd
     INTEGER, OPTIONAL :: GaussPoint
     LOGICAL, OPTIONAL :: Found
     REAL(KIND=dp) :: Val
     
     TYPE( Element_t), POINTER :: pElement
     INTEGER :: j, n
     INTEGER :: Indexes(100)
     LOGICAL :: SameElement
          
     Val = 0.0_dp
     
     IF( PRESENT( Found ) ) Found = .FALSE.
     
     IF( .NOT. ASSOCIATED( Handle % Variable ) ) RETURN

     ! Find the pointer to the element, if not given
     IF( PRESENT( Element ) ) THEN
       PElement => Element
     ELSE
       PElement => CurrentModel % CurrentElement
     END IF
     
     SameElement = ASSOCIATED( Handle % Element, pElement )
     IF( SameElement ) THEN
       IF( .NOT. Handle % ActiveElement ) RETURN
     ELSE
       Handle % Element => pElement
     END IF
     
     ! If variable is defined on gauss points return that instead
     IF( Handle % Variable % TYPE == Variable_on_gauss_points ) THEN
       IF( .NOT. PRESENT( GaussPoint ) ) THEN
         CALL Fatal('GetElementScalar','GaussPoint required for gauss point variable!')
       END IF
       
       j = pElement % ElementIndex
       n = Handle % Perm(j+1) - Handle % Perm(j)
       Handle % ActiveElement = ( n > 0 ) 
       
       IF( n == 0 ) RETURN             
       val = Handle % Values(Handle % Perm(j) + GaussPoint )

     ELSE IF( Handle % Variable % TYPE == Variable_on_elements ) THEN       
       j = Handle % Perm( pElement % ElementIndex ) 
       Handle % ActiveElement = ( j > 0 ) 
       
       IF( j == 0 ) RETURN             
       val = Handle % Values(j) 
       
     ELSE
       IF( .NOT. PRESENT( Basis ) ) THEN
         CALL Fatal('GetElementScalar','Basis required for non gauss-point variable!')
       END IF
       
       IF( .NOT. SameElement ) THEN
         IF( Handle % Variable % TYPE == Variable_on_nodes_on_elements ) THEN       
           n = pElement % TYPE % NumberOfNodes
           Indexes(1:n) = pElement % DGIndexes(1:n)
         ELSE
           n = pElement % TYPE % NumberOfNodes
           Indexes(1:n) = pElement % NodeIndexes
         END IF

         IF( PRESENT( nd ) ) THEN
           n = MIN( nd, n ) 
         END IF
         
         IF( ASSOCIATED( Handle % Perm ) ) THEN
           Handle % ActiveElement = ALL( Handle % Perm( Indexes(1:n) ) /= 0 )
           IF(.NOT. Handle % ActiveElement ) RETURN
           Handle % ElementValues(1:n) = Handle % Values( Handle % Perm( Indexes(1:n) ) )
         ELSE
           Handle % ActiveElement = .TRUE.
           Handle % ElementValues(1:n) = Handle % Values( Indexes(1:n) )           
         END IF
       ELSE
         n = Handle % n
       END IF
       val = SUM( Basis(1:n) * Handle % ElementValues(1:n) )
     END IF

     IF( PRESENT( Found ) ) Found = .TRUE.
     
   END FUNCTION ListGetElementScalarSolution
!------------------------------------------------------------------------------


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
           CALL Fatal("ListGetConstRealArray", "Failed to find: "//TRIM(Name) )
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
!> Gets an 1D constant real array from the list by its name.   
!------------------------------------------------------------------------------
   RECURSIVE FUNCTION ListGetConstRealArray1( List,Name,Found,UnfoundFatal ) RESULT( F )
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: List
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found, UnfoundFatal
!------------------------------------------------------------------------------
     REAL(KIND=dp), POINTER  :: F(:)
     INTEGER :: i,j,N1,N2
     TYPE(ValueListEntry_t), POINTER :: ptr
!------------------------------------------------------------------------------
     NULLIFY( F ) 
     ptr => ListFind(List,Name,Found)
     IF (.NOT.ASSOCIATED(ptr) ) THEN
       IF(PRESENT(UnfoundFatal)) THEN
         IF(UnfoundFatal) THEN
           CALL Fatal("ListGetConstRealArray1","Failed to find: "//TRIM(Name))
         END IF
       END IF
       RETURN
     END IF

     IF ( .NOT. ASSOCIATED(ptr % FValues) ) THEN
       WRITE(Message,*) 'Value type for property [', TRIM(Name), &
               '] not used consistently.'
       CALL Fatal( 'ListGetConstRealArray1', Message )
       RETURN
     END IF

     N1 = SIZE( ptr % FValues,1 )
     N2 = SIZE( ptr % FValues,2 )
     IF( N2 > 1 ) THEN
       CALL Warn('ListGetConstRealArray1','The routine is designed for 1D arrays!')
     END IF
       
     F => ptr % FValues(:,1,1)

   END FUNCTION ListGetConstRealArray1
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
#ifdef HAVE_LUA
           IF ( .not. ptr % LuaFun ) THEN
#endif
             DO l=1,j
               WRITE( cmd, '(a,g19.12)' ) 'tx('//TRIM(i2s(l-1))//')=', T(l)
               k1 = LEN_TRIM(cmd)
               CALL matc( cmd, tmp_str, k1 )
             END DO

             cmd = ptr % CValue
             k1 = LEN_TRIM(cmd)
             CALL matc( cmd, tmp_str, k1 )
             READ( tmp_str(1:k1), * ) ((F(j,k,i),k=1,N2),j=1,N1)
#ifdef HAVE_LUA
           ELSE
             call ElmerEvalLuaT(LuaState, ptr, T, F(:,:,i), j)
           END IF
#endif
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
     ELSE
       Found = lFound
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
#ifdef HAVE_LUA
           IF ( .not. ptr % LuaFun ) THEN
#endif
           DO l=1,j
             WRITE( cmd, '(a,g19.12)' ) 'tx('//TRIM(i2s(l-1))//')=', T(l)
             k1 = LEN_TRIM(cmd)
             CALL matc( cmd, tmp_str, k1 )
           END DO

           cmd = ptr % CValue
           k1 = LEN_TRIM(cmd)
           CALL matc( cmd, tmp_str, k1 )
           READ( tmp_str(1:k1), * ) (G(j,i),j=1,N1)
#ifdef HAVE_LUA
           ELSE
             call ElmerEvalLuaV(LuaState, ptr, T, G(:,i), j)
           END IF
#endif
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
!> Check if the keyword is present in any boundary condition.
!------------------------------------------------------------------------------
   FUNCTION ListCheckPresentAnyIC( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
     INTEGER :: ic
     
     Found = .FALSE.
     DO ic = 1,Model % NumberOfICs
       Found = ListCheckPresent( Model % ICs(ic) % Values, Name )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListCheckPresentAnyIC
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
!> Check if the keyword is true in any body.
!------------------------------------------------------------------------------
   FUNCTION ListGetCRealAnyBody( Model, Name, Found ) RESULT( F )
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL, OPTIONAL :: Found
     REAL(KIND=dp) :: F
     
     INTEGER :: body
     LOGICAL :: GotIt
     
     F = 0.0_dp
     GotIt = .FALSE.
     DO body = 1,Model % NumberOfBodies
       F = ListGetCReal( Model % Bodies(body) % Values, Name, GotIt )
       IF( GotIt ) EXIT
     END DO

     IF( PRESENT( Found ) ) Found = GotIt
     
!------------------------------------------------------------------------------
   END FUNCTION ListGetCRealAnyBody
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
!> Check if the keyword is present in any solver.
!------------------------------------------------------------------------------
   FUNCTION ListCheckPresentAnySolver( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found
     INTEGER :: ind
     
     Found = .FALSE.
     DO ind = 1,Model % NumberOfSolvers
       Found = ListCheckPresent( Model % Solvers(ind) % Values, Name )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListCheckPresentAnySolver
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Check if the keyword is present in any component.
!------------------------------------------------------------------------------
  FUNCTION ListCheckPresentAnyComponent( Model, Name ) RESULT( Found )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    TYPE(Model_t) :: Model
    CHARACTER(LEN=*) :: Name
    LOGICAL :: Found
    INTEGER :: ind
        
    Found = .FALSE.
    DO ind=1, Model % NumberOfComponents
      Found = ListCheckPresent( Model % Components(ind) % Values, Name )
      IF( Found ) EXIT
    END DO   
!------------------------------------------------------------------------------
  END FUNCTION ListCheckPresentAnyComponent
!------------------------------------------------------------------------------  

  !------------------------------------------------------------------------------
!> Check if the keyword is true in any component.
!------------------------------------------------------------------------------
  FUNCTION ListGetLogicalAnyComponent( Model, Name ) RESULT( Found )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    TYPE(Model_t) :: Model
    CHARACTER(LEN=*) :: Name
    LOGICAL :: Found, GotIt
    INTEGER :: ind
        
    Found = .FALSE.
    DO ind=1, Model % NumberOfComponents
      Found = ListGetLogical( Model % Components(ind) % Values, Name, GotIt )
      IF( Found ) EXIT
    END DO
!------------------------------------------------------------------------------
  END FUNCTION ListGetLogicalAnyComponent
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
!> Check if the keyword is True in any solver.
!------------------------------------------------------------------------------
   FUNCTION ListGetLogicalAnySolver( Model, Name ) RESULT(Found)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     CHARACTER(LEN=*) :: Name
     LOGICAL :: Found, GotIt
     INTEGER :: ind
     
     Found = .FALSE.
     DO ind = 1,Model % NumberOfSolvers
       Found = ListGetLogical( Model % Solvers(ind) % Values, Name, GotIt )
       IF( Found ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ListGetLogicalAnySolver
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
  SUBROUTINE CreateListForSaving( Model, List, ShowVariables, ClearList, &
      UseGenericKeyword )
!------------------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(ValueList_t), POINTER  :: List
    LOGICAL :: ShowVariables
    LOGICAL, OPTIONAL :: ClearList
    LOGICAL, OPTIONAL :: UseGenericKeyword
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,LoopDim, VarDim,FullDim,DOFs,dim,Comp
    TYPE(Variable_t), POINTER :: Variables, Var, Var1
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName, VarStr, VarStrComp, VarStrExt, str
    LOGICAL :: IsVector, Set, GotIt, ComponentVector, ThisOnly, IsIndex, &
        EnforceVectors, UseGeneric
    INTEGER :: Nvector, Nscalar
    TYPE(ValueList_t), POINTER :: Params

    Params => Model % Solver % Values
    Variables => Model % Mesh % Variables

    IF( .NOT. ASSOCIATED( Variables ) ) THEN
      CALL Warn('CreateListForSaving','Mesh does not include any variables!')
      RETURN
    END IF
    
    UseGeneric = .FALSE.
    IF( PRESENT( UseGenericKeyword ) ) THEN
      UseGeneric = UseGenericKeyword 
    END IF
    

!------------------------------------------------------------------------------
! Sometimes the list must be cleared in order to use it for a different mesh
!-----------------------------------------------------------------------------
    IF( PRESENT( ClearList ) ) THEN
      IF( ClearList ) THEN
        IF( UseGeneric ) THEN
          DO i=1,999
            WRITE(VarStr,'(A,I0)') 'Variable ',i
            IF( ListCheckPresent( List, VarStr ) ) THEN
              CALL ListRemove( List, VarStr )
            ELSE
              EXIT
            END IF
          END DO
        ELSE
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
            
            WRITE(VarStr,'(A,I0,A)') 'Vector Field ',i,' Complement'
            IF( ListCheckPresent( List, VarStr ) ) THEN
              CALL ListRemove( List, VarStr )
            END IF
          END DO
          
        END IF
      END IF
    END IF
    
    !-------------------------------------------------------------------
    ! First check that there is a need to create the list i.e. it is not
    ! already manually defined
    !-------------------------------------------------------------------
    IF( UseGeneric ) THEN
      IF( ListCheckPresent( List,'Variable 1' ) ) THEN
        CALL Info('CreateListForSaving','Variable 1 exists, creating no list!',Level=10)
        RETURN
      END IF
    ELSE
      IF( ListCheckPresent( List,'Scalar Field 1' ) ) THEN
        CALL Info('CreateListForSaving','Scalar Field 1 exists, creating no list!',Level=10)
        RETURN
      END IF
      
      IF( ListCheckPresent( List,'Vector Field 1' ) ) THEN
        CALL Info('CreateListForSaving','Vector Field 1 exists, creating no list!',Level=10)
        RETURN
      END IF
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

      IF( Var % TYPE == Variable_global ) THEN
        Var => Var % Next
        CYCLE        
      ELSE IF( Var % TYPE == Variable_on_gauss_points ) THEN
        CONTINUE

      ELSE IF( Var % TYPE == Variable_on_elements ) THEN
        CONTINUE

      END IF

      ! Skip if variable is otherwise strange in size
      IF(.NOT. ASSOCIATED( Var % Perm ) ) THEN
        IF( Var % TYPE == Variable_on_nodes ) THEN
          IF( SIZE( Var % Values ) /= Var % Dofs * Model % Mesh % NumberOfNodes ) THEN
            Var => Var % Next
            CYCLE
          END IF
        ELSE IF( Var % TYPE == Variable_on_nodes_on_elements ) THEN
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
        Set = .TRUE.
        IF(.NOT. UseGeneric ) THEN
          Var1 => Variables
          DO WHILE( ASSOCIATED( Var1 ) )
            IF ( TRIM(Var1 % Name) == 'displacement' ) EXIT
            Var1 => Var1 % Next
          END DO
          IF ( ASSOCIATED( Var1 ) ) Set = .FALSE.
        END IF
        
      CASE('mesh update 1','mesh update 2', 'mesh update 3' )
        
      CASE( 'displacement' )
        Set = .TRUE.
        ! mesh update is by default the complement to displacement 
        ! However, for generic variablelist the complement is not active
        IF(.NOT. UseGeneric ) THEN
          Var1 => Variables
          DO WHILE( ASSOCIATED( Var1 ) )
            IF ( TRIM(Var1 % Name) == 'mesh update' ) EXIT
            Var1 => Var1 % Next
          END DO
          IF ( ASSOCIATED( Var1 ) ) THEN
            WRITE(VarStrComp,'(A,I0,A)') 'Vector Field ',Nvector+1,' Complement'
            CALL ListAddString( List ,TRIM(VarStrComp),'mesh update')
          END IF
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
          ! The size of the vector can be either dim or 3. 
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
              ! If we have the 1st component we need at least dim (2 or 3) components
              ! to have a vector.
              Var1 => VariableGet(Variables,TRIM(str(1:j-2))//' '//I2S(dim),ThisOnly)		

              ! However, if the 4th component also exists then this cannot be a vector
              IF( ASSOCIATED(Var1)) THEN
                Var1 => VariableGet(Variables,TRIM(str(1:j-2))//' '//I2S(4),ThisOnly)		
                IsVector = .NOT. ASSOCIATED(Var1)
              END IF
              
            ELSE IF( Comp <= 3 ) THEN  ! component 2 or 3
              ! Associated to the previous case, cycle the other components of the vector
              ! and cycle them if they are part of the vector that will be detected above.
 
              Var1 => VariableGet(Variables,TRIM(str(1:j-2))//' 1',ThisOnly)		
              IF( ASSOCIATED( Var1 ) ) THEN
                Var1 => VariableGet(Variables,TRIM(str(1:j-2))//' '//I2S(4),ThisOnly)		
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
        IF( UseGeneric ) THEN
          Nscalar = Nscalar + 1
          WRITE(VarStr,'(A,I0)') 'Variable ',Nscalar          
        ELSE IF( IsVector ) THEN          
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
      IF( UseGeneric ) THEN
        DO i=1,Nscalar
          WRITE(VarStr,'(A,I0)') 'Variable ',i
          VarName = ListGetString( List, VarStr,GotIt )
          IF( GotIt ) THEN
            WRITE( Message,'(A)') TRIM(VarStr)//': '//TRIM(VarName)
            CALL Info('CreateListForSaving',Message,Level=6)
          END IF
        END DO
      ELSE
        DO i=1,Nscalar
          WRITE(VarStr,'(A,I0)') 'Scalar Field ',i
          VarName = ListGetString( List, VarStr,GotIt )
          IF( GotIt ) THEN
            WRITE( Message,'(A)') TRIM(VarStr)//': '//TRIM(VarName)
            CALL Info('CreateListForSaving',Message,Level=6)
          END IF
        END DO
        
        DO i=1,Nvector
          WRITE(VarStr,'(A,I0)') 'Vector Field ',i
          VarName = ListGetString( List, VarStr,GotIt )
          IF( GotIt ) THEN
            WRITE( Message,'(A)') TRIM(VarStr)//': '//TRIM(VarName)
            CALL Info('CreateListForSaving',Message,Level=6)
          END IF
        END DO
        
        DO i=1,Nvector
          WRITE(VarStr,'(A,I0,A)') 'Vector Field ',i,' Complement'
          VarName = ListGetString( List, VarStr, GotIt )
          IF( GotIt ) THEN
            WRITE( Message,'(A)') TRIM(VarStr)//': '//TRIM(VarName)
            CALL Info('CreateListForSaving',Message,Level=6)
          END IF
        END DO
      END IF
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
      TimerCumulative = ListGetLogical( CurrentModel % Simulation,'Timer Cumulative',Found)      
      TimerRealTime = ListGetLogical( CurrentModel % Simulation,'Timer Real Time',Found)      
      TimerCPUTime = ListGetLogical( CurrentModel % Simulation,'Timer CPU Time',Found)            
      IF( .NOT. (TimerRealTime .OR. TimerCPUTime ) ) TimerRealTime = .TRUE.
      TimerPrefix = ListGetString( CurrentModel % Simulation,'Timer Prefix',Found )
      IF( .NOT. Found ) THEN
        IF( ListGetLogical( CurrentModel % Simulation,'Timer Results',Found ) ) THEN
          TimerPrefix = 'res:'
        ELSE
          TimerPrefix = 'timer:'
        END IF
      END IF
    END IF

    
    IF( TimerPassive ) RETURN

    IF( TimerCPUTime ) THEN
      ct = CPUTime()
      CALL ListAddConstReal( TimerList,TRIM(TimerName)//' cpu time',ct )
    END IF

    IF( TimerRealTime ) THEN
      rt = RealTime()
      CALL ListAddConstReal( TimerList,TRIM(TimerName)//' real time',rt )
    END IF

    IF( TimerCumulative ) THEN
      IF( TimerCPUTime ) THEN
        IF( .NOT. ListCheckPresent( CurrentModel % Simulation,TRIM(TimerPrefix)//' '//TRIM(TimerName)//' cpu time') ) THEN
          CALL ListAddConstReal( CurrentModel % Simulation,TRIM(TimerPrefix)//' '//TRIM(TimerName)//' cpu time',0.0_dp )
        END IF
      END IF
      IF( TimerRealTime ) THEN
        IF( .NOT. ListCheckPresent( CurrentModel % Simulation,TRIM(TimerPrefix)//' '//TRIM(TimerName)//' real time') ) THEN
          CALL ListAddConstReal( CurrentModel % Simulation,TRIM(TimerPrefix)//' '//TRIM(TimerName)//' real time',0.0_dp )
        END IF
      END IF
    END IF
      
  END SUBROUTINE ResetTimer

  
!-----------------------------------------------------------------------------
!> Delete an existing timer.
!----------------------------------------------------------------------------
  SUBROUTINE DeleteTimer(TimerName) 
    CHARACTER(*) :: TimerName
    
    IF( TimerPassive ) RETURN

    IF( TimerCPUTime ) THEN
      CALL ListRemove( TimerList, TRIM(TimerName)//' cpu time' ) 
    END IF

    IF( TimerRealTime ) THEN
      CALL ListRemove( TimerList, TRIM(TimerName)//' real time' ) 
    END IF
      
  END SUBROUTINE DeleteTimer
 
!-----------------------------------------------------------------------------
!> Check current time of the timer.
!----------------------------------------------------------------------------
  SUBROUTINE CheckTimer(TimerName, Level, Delete, Reset) 
    CHARACTER(*) :: TimerName
    INTEGER, OPTIONAL :: Level
    LOGICAL, OPTIONAL :: Reset, Delete
    
    REAL(KIND=dp) :: ct0,rt0,ct, rt, cumct, cumrt
#ifndef USE_ISO_C_BINDINGS
    REAL(KIND=dp) :: RealTime, CPUTime
#endif
    LOGICAL :: Found

    IF( TimerPassive ) RETURN
    
    IF( TimerCPUTime ) THEN
      ct0 = ListGetConstReal( TimerList,TRIM(TimerName)//' cpu time',Found) 
      IF( Found ) THEN
        ct = CPUTime() - ct0
        WRITE(Message,'(a,f10.4,a)') 'Elapsed CPU time: ',ct,' (s)'
        CALL Info(TRIM(TimerName),Message,Level=Level)          
      END IF
    END IF

    IF( TimerRealTime ) THEN
      rt0 = ListGetConstReal( TimerList,TRIM(TimerName)//' real time',Found)
      IF( Found ) THEN
        rt = RealTime() - rt0       
        WRITE(Message,'(a,f10.4,a)') 'Elapsed REAL time: ',rt,' (s)'
        CALL Info(TRIM(TimerName),Message,Level=Level)          
      END IF
    END IF
    
    
    IF( TimerCPUTime ) THEN
      IF( TimerCumulative ) THEN
        cumct = ListGetConstReal(CurrentModel % Simulation,&
            TRIM(TimerPrefix)//' '//TRIM(TimerName)//' cpu time',Found)
        IF( Found ) THEN
          ct = ct + cumct
          WRITE(Message,'(a,f10.4,a)') 'Elapsed CPU time cumulative: ',ct,' (s)'
          CALL Info(TRIM(TimerName),Message,Level=Level)          
        ELSE
          CALL Warn('CheckTimer',&
              'Requesting previous CPU time from non-existing timer: '//TRIM(TimerName) )            
        END IF
      END IF
      CALL ListAddConstReal(CurrentModel % Simulation,&
          TRIM(TimerPrefix)//' '//TRIM(TimerName)//' cpu time',ct)

    END IF
    IF( TimerRealTime ) THEN
      IF( TimerCumulative ) THEN
        cumrt = ListGetConstReal(CurrentModel % Simulation,&
            TRIM(TimerPrefix)//' '//TRIM(TimerName)//' real time',Found)
        IF( Found ) THEN
          rt = rt + cumrt
          WRITE(Message,'(a,f10.4,a)') 'Elapsed real time cumulative: ',rt,' (s)'
          CALL Info(TRIM(TimerName),Message,Level=Level)          
        ELSE
          CALL Warn('CheckTimer',&
              'Requesting previous real time from non-existing timer: '//TRIM(TimerName) )            
        END IF
      END IF
      CALL ListAddConstReal(CurrentModel % Simulation,&
          TRIM(TimerPrefix)//' '//TRIM(TimerName)//' real time',rt)        
    END IF
      
    
    IF( PRESENT( Reset ) ) THEN
      IF( Reset ) THEN
        IF( TimerCPUTime ) THEN
          CALL ListAddConstReal( TimerList,TRIM(TimerName)//' cpu time',ct )
        END IF
        IF( TimerRealTime ) THEN
          CALL ListAddConstReal( TimerList,TRIM(TimerName)//' real time',rt )
        END IF
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
   
#ifdef HAVE_LUA
!-------------------------------------------------------------------------------
!> evaluates lua string to real array 
!-------------------------------------------------------------------------------
SUBROUTINE ElmerEvalLuaT(L, ptr, T, F, varcount)
!-------------------------------------------------------------------------------
  TYPE(LuaState_t) :: L
  TYPE(ValueListEntry_t), POINTER :: ptr
  REAL(KIND=C_DOUBLE), INTENT(IN) :: T(:)
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: F(:,:)
  INTEGER :: VARCOUNT
!-------------------------------------------------------------------------------
  integer :: lstat

  L % tx(1:varcount) = T(1:varcount) ! this should be superfluous
  call lua_exec_fun(L, ptr % cvalue, 0, size(F,1)*size(F,2))

  CALL lua_poptensor(L, F)
!-------------------------------------------------------------------------------
END SUBROUTINE
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> evaluates lua string to real vector 
!-------------------------------------------------------------------------------
SUBROUTINE ElmerEvalLuaV(L, ptr, T, F, varcount)
!-------------------------------------------------------------------------------
  TYPE(LuaState_t) :: L
  TYPE(ValueListEntry_t), POINTER :: ptr
  REAL(KIND=C_DOUBLE), INTENT(IN) :: T(:)
  REAL(KIND=C_DOUBLE), INTENT(INOUT) :: F(:)
  INTEGER :: VARCOUNT
!-------------------------------------------------------------------------------
  integer :: lstat

  L % tx(1:varcount) = T(1:varcount) ! this should be superfluous
  call lua_exec_fun(L, ptr % cvalue, 0, size(F,1))

  CALL lua_popvector(L, F)
!-------------------------------------------------------------------------------
END SUBROUTINE
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> evaluates lua string to real scalar 
!-------------------------------------------------------------------------------
SUBROUTINE ElmerEvalLuaS(L, ptr, T, F, varcount)
!-------------------------------------------------------------------------------
  TYPE(LuaState_t) :: L
  TYPE(ValueListEntry_t), POINTER :: ptr
  REAL(KIND=C_DOUBLE), INTENT(IN) :: T(:)
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: F
  INTEGER :: VARCOUNT
!-------------------------------------------------------------------------------
  integer :: lstat

  L % tx(1:varcount) = T(1:varcount) ! this should be superfluous
  call lua_exec_fun(L, ptr % cvalue, 0, 1)
  F = lua_popnumber(LuaState)
!-------------------------------------------------------------------------------
END SUBROUTINE
!-------------------------------------------------------------------------------
#endif

#ifdef DEVEL_LISTCOUNTER
   
   !------------------------------------------------------------------------------
   !> Go through the lists and for each lists show call counts.
   !------------------------------------------------------------------------------
   SUBROUTINE ReportListCounters( Model ) 
     TYPE(Model_t) :: Model

     INTEGER :: i, totcount, nelem     

     CALL Info('ReportListCounters','Saving ListGet operations count per bulk elements')

     ! OPEN(10,FILE="listcounter.dat")
     OPEN( 10,File='../listcounter.dat',&
         STATUS='UNKNOWN',POSITION='APPEND' )

     totcount = 0

     ! These are only for reference
     nelem = Model % Mesh % NumberOfBulkElements
     
     CALL ReportList('Simulation', Model % Simulation )
     CALL ReportList('Constants', Model % Constants )
     DO i=1,Model % NumberOfEquations
       CALL ReportList('Equation '//TRIM(I2S(i)), Model % Equations(i) % Values)
     END DO
     DO i=1,Model % NumberOfComponents
       CALL ReportList('Component '//TRIM(I2S(i)), Model % Components(i) % Values )
     END DO
     DO i=1,Model % NumberOfBodyForces
       CALL ReportList('Body Force '//TRIM(I2S(i)), Model % BodyForces(i) % Values )
     END DO
     DO i=1,Model % NumberOfICs
       CALL ReportList('Initial Condition '//TRIM(I2S(i)), Model % ICs(i) % Values )
     END DO
     DO i=1,Model % NumberOfBCs
       CALL ReportList('Boundary Condition '//TRIM(I2S(i)), Model % BCs(i) % Values )
     END DO
     DO i=1,Model % NumberOfMaterials
       CALL ReportList('Material '//TRIM(I2S(i)), Model % Materials(i) % Values )
     END DO
     DO i=1,Model % NumberOfBoundaries
       CALL ReportList('Boundary '//TRIM(I2S(i)), Model % Boundaries(i) % Values )
     END DO     
     DO i=1,Model % NumberOfSolvers
       CALL ReportList('Solver '//TRIM(I2S(i)), Model % Solvers(i) % Values )
     END DO
     CLOSE(10)
     
     WRITE(Message,'(I0,T10,A,F8.3,A)') totcount,"List operations total (",&
         1.0_dp * totcount / Model % Mesh % NumberOfBulkElements," per element )"
     CALL Info('ReportListCounters',Message)
     

   CONTAINS

     
     !------------------------------------------------------------------------------
     ! Plot the number of times that the list entries have been called.
     !------------------------------------------------------------------------------
     SUBROUTINE ReportList( SectionName, List )
       TYPE(ValueList_t), POINTER :: List
       CHARACTER(LEN=*) :: SectionName
       !------------------------------------------------------------------------------
       TYPE(ValueListEntry_t), POINTER :: ptr
       INTEGER :: n, m

       IF(.NOT.ASSOCIATED(List)) RETURN

       Ptr => List % Head
       DO WHILE( ASSOCIATED(ptr) )
         n = ptr % NameLen
         m = ptr % Counter 
         
         WRITE( 10,'(F10.5,T12,A,T35,A)') &
             1.0*m / nelem, TRIM(SectionName),ptr % Name(1:n)
         totcount = totcount + m 
         ptr => ptr % Next
       END DO

     END SUBROUTINE ReportList
     !------------------------------------------------------------------------------    
     
   END SUBROUTINE ReportListCounters
  !------------------------------------------------------------------------------

#else

   SUBROUTINE ReportListCounters( Model ) 
     TYPE(Model_t) :: Model

     CALL Info('ReportListCounter','List counters are not activated!')
   END SUBROUTINE ReportListCounters
      
#endif


   

END MODULE Lists

!> \} ElmerLib
