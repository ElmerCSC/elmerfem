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
! *  Original Date: 30 Mar 2011
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!-------------------------------------------------------------------------------
!>  Basic FETI solver, for Poisson & Navier style of equations.
!-------------------------------------------------------------------------------

#include "../config.h"

MODULE FetiSolve

  USE MeshBasics, ONLY : FindRigidBodyFixingNodes
  USE DefUtils
  IMPLICIT NONE

  TYPE(ValueList_t), PRIVATE, POINTER :: Params => Null()


  ! Communication buffer stuff we send to our neighbours:
  ! -----------------------------------------------------
  TYPE toSend_t
    INTEGER :: n
    REAL(KIND=dp), ALLOCATABLE :: buf(:)
    INTEGER, ALLOCATABLE :: ifg(:), perm(:)
  END TYPE toSend_t

  ! Communication buffer stuff we receive from our neighbours:
  ! ----------------------------------------------------------
  TYPE toReceive_t
    INTEGER :: n
    INTEGER, ALLOCATABLE :: perm(:)
    REAL(KIND=dp), ALLOCATABLE :: buf(:)
  END TYPE toReceive_t

  ! The connectivity matrix 'B':
  ! ----------------------------
  TYPE(Matrix_t), POINTER, SAVE :: Bmat=>Null()

  ! the null space or kernel of the input system A, also called
  ! R in the relevant FETI literature:
  ! ------------------------------------------------------------
  INTEGER, PRIVATE, SAVE :: nz
  INTEGER, PRIVATE, PARAMETER :: maxnz=20
  REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE ::z(:,:)

  ! Neighbour identification, local and global numbering of 
  ! neighbour PEs:
  ! -------------------------------------------------------
  INTEGER, PRIVATE, SAVE :: nneigh
  INTEGER, ALLOCATABLE, PRIVATE, SAVE :: FixInds(:)
  INTEGER, ALLOCATABLE, PRIVATE, SAVE :: lpnum(:), gpnum(:)

  ! Some global flags:
  ! ------------------
  LOGICAL, PRIVATE, SAVE :: FixingByDofs=.FALSE.
  LOGICAL, PRIVATE, SAVE :: LocalSolveChecked=.TRUE.
  LOGICAL, PRIVATE, SAVE :: Precondition=.FALSE., TotalFETI=.FALSE., NullSpaceLC=.TRUE.,&
                          dumptofiles=.FALSE.
  LOGICAL, PRIVATE, SAVE :: InitializeIf, InitializeLC, CPG, Refactorize, FetiAsPrec=.FALSE.

  LOGICAL, ALLOCATABLE, PRIVATE, SAVE :: DirichletDOFs(:)

integer,allocatable::snd(:),asize(:),bsize(:),cnt(:),gbuf(:,:),ibuf(:,:)
integer::ierr,me,abeg,bbeg
integer::status(MPI_STATUS_SIZE)

#include "huti_fdefs.h"

CONTAINS
 
  !> Send given buffer to given neighbour, either 'tags' in 'ifg' 
  !> (when initializing) or the interface values in 'buf':
  ! ------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE FetiSend(proc, nin, buf, ifg, tag)
!------------------------------------------------------------------------------
     INTEGER :: tag
     INTEGER, OPTIONAL :: ifg(:)
     INTEGER :: proc, nin
     REAL(KIND=dp), OPTIONAL :: buf(:)
!------------------------------------------------------------------------------
     INTEGER :: ierr, zcnt, n
!------------------------------------------------------------------------------
     n = nin
     IF (PRESENT(buf).AND.n>0) THEN
       zcnt=COUNT(buf(1:n)==0)
       IF (zcnt==n) n=0
     END IF

     CALL MPI_BSEND( n, 1, MPI_INTEGER, &
          proc, tag, ELMER_COMM_WORLD, ierr )

     IF (n>0) THEN
       IF (PRESENT(buf)) THEN
         CALL MPI_BSEND( buf, n, MPI_DOUBLE_PRECISION, &
            proc, tag+1, ELMER_COMM_WORLD, ierr )
       END IF

       IF (PRESENT(ifg)) THEN
         CALL MPI_BSEND( ifg, n, MPI_INTEGER, &
           proc, tag+2, ELMER_COMM_WORLD, ierr )
       END IF
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE FetiSend
!------------------------------------------------------------------------------


  !> Receive given buffer from a neighbour, either 'tags' to 'ifg'
  !> (when initializing) or the interface values to 'buf', the id
  !> of the sender returned in 'proc':
  ! ---------------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE FetiRecv(proc, n, buf, ifg, tag)
!------------------------------------------------------------------------------
     INTEGER, OPTIONAL, ALLOCATABLE :: ifg(:)
     INTEGER :: proc, n, tag
     REAL(KIND=dp), OPTIONAL :: buf(:)
!------------------------------------------------------------------------------
     INTEGER :: status(MPI_STATUS_SIZE)=0, ierr=0
!------------------------------------------------------------------------------
     CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
           tag, ELMER_COMM_WORLD, status, ierr )
     proc = status(MPI_SOURCE)

     IF (n>0) THEN
       IF (PRESENT(buf)) THEN
         CALL MPI_RECV( buf, n, MPI_DOUBLE_PRECISION, proc, &
             tag+1, ELMER_COMM_WORLD, status, ierr )
       END IF

       IF (PRESENT(ifg)) THEN
         IF ( ALLOCATED(ifg) ) THEN
           IF (SIZE(ifg)<n) DEALLOCATE(ifg)
         END IF
         IF ( .NOT. ALLOCATED(ifg)) ALLOCATE(ifg(n))
         CALL MPI_RECV( ifg, n, MPI_INTEGER, proc, &
             tag+2, ELMER_COMM_WORLD, status, ierr )
       END IF
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE FetiRecv
!------------------------------------------------------------------------------


   ! Identify neighbour partitions:
   ! ------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE FetiGetNeighbours()
!------------------------------------------------------------------------------
     INTEGER :: i
!------------------------------------------------------------------------------
     IF ( ALLOCATED(gpnum) ) DEALLOCATE(gpnum)
     IF ( ALLOCATED(lpnum) ) DEALLOCATE(lpnum)

     ALLOCATE(gpnum(ParEnv % NumOfNeighbours),lpnum(0:ParEnv % PEs-1))
     lpnum = 0; nneigh = 0

     DO i=0,ParEnv % PEs-1
       IF (ParEnv % IsNeighbour(i+1) .AND. ParEnv % Active(i+1)) THEN
         nneigh=nneigh+1
         lpnum(i) = nneigh
         gpnum(nneigh) = i
       END IF
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE FetiGetNeighbours
!------------------------------------------------------------------------------


   ! Upper bound for the number of Lagrange coefficients.
   !
   ! One is created for each pair in which an owned interface DOF participates,
   ! so a DOF shared by k partitions contributes k-1 of them, and under total
   ! FETI each Dirichlet DOF adds one more. The total is therefore NOT bounded
   ! by the number of rows of the local matrix, which is what the work vectors
   ! used to be sized by. Partition finely enough -- twelve pieces of the winkel
   ! mesh is plenty -- and most DOFs end up on an interface, many of them shared
   ! three or four ways, and the count passes the row count.
   ! ------------------------------------------------------------------------
!------------------------------------------------------------------------------
   FUNCTION FetiMaxLC(A) RESULT(m)
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: m
!------------------------------------------------------------------------------
     INTEGER :: i,n
     LOGICAL, POINTER :: ig(:)
     TYPE(NeighbourList_t), POINTER :: nb(:)
!------------------------------------------------------------------------------
     m = 0
     IF (.NOT.ASSOCIATED(A % ParallelInfo)) RETURN

     n = A % NumberOfRows
     ig => A % ParallelInfo % GInterface
     nb => A % ParallelInfo % NeighbourList

     DO i=1,n
       IF (.NOT.ig(i)) CYCLE
       IF (.NOT.(TotalFeti.OR..NOT.DirichletDOFs(i))) CYCLE
       IF (nb(i) % Neighbours(1) /= ParEnv % myPE) CYCLE
       m = m + SIZE(nb(i) % Neighbours)-1
     END DO

     IF (TotalFeti) m = m + COUNT(DirichletDOFs)
!------------------------------------------------------------------------------
   END FUNCTION FetiMaxLC
!------------------------------------------------------------------------------


  !> First time communication with all neighbours where local
  !> addresses of global tags are identified:
  ! --------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE FetiSendReceiveInit(sndLC, gdofs, ldofs, &
           procs, toSend, toReceive, tag)
!------------------------------------------------------------------------------
      INTEGER :: sndLC, tag, gdofs(:), ldofs(:), procs(:)
      TYPE(toSend_t) :: toSend(:)  
      TYPE(toReceive_t) :: toReceive(:)  
!------------------------------------------------------------------------------
      INTEGER :: i,j,k,l,n,m,proc,lproc
      LOGICAL :: Found
      INTEGER, ALLOCATABLE :: gorder(:), igorder(:), ifg(:)
!------------------------------------------------------------------------------
      DO i=1,nneigh
        proc = gpnum(i)
        CALL FetiSend(proc, toSend(i) % n, ifg=toSend(i) % ifg, tag=tag)
      END DO 

      ! Receive interface parts and store indices
      ! -----------------------------------------
      ALLOCATE(gorder(sndLC), igorder(sndLC), ifg(sndLC))
      gorder=[(i,i=1,sndLC)]
      CALL SortI(sndLC, gdofs, gorder)
      igorder(gorder)=[(i,i=1,sndLC)]

      DO i=1,nneigh
        CALL FetiRecv(proc, n, ifg=ifg, tag=tag)

        lproc=lpnum(proc)
        toReceive(lproc) % n = n
        ALLOCATE(toReceive(lproc) % perm(n))

        IF (n<=0) CYCLE
        toReceive(lproc) % perm=0
        DO j=1,n
          k=SearchIAItem(sndLC,gdofs,ifg(j),gorder)
          IF ( k<=0 ) THEN
            PRINT*,'should not happen: ', parenv % mype, proc, ifg(j)
            CYCLE
          END IF

          IF (proc/=procs(k)) THEN

            ! Account for multiple global dof tags in gdofs; they should
            ! be adjacent in the original order too:
            ! ----------------------------------------------------------
            Found=.FALSE.
            DO l=k+1,sndLC
              m = igorder(l)
              IF ( gdofs(m)/=ifg(j) ) EXIT
              IF ( proc==procs(l) ) THEN
                Found=.TRUE.; EXIT
              END IF
            END DO

            IF (.NOT.Found) THEN
              DO l=k-1,1,-1
                m = igorder(l)
                IF ( gdofs(m)/=ifg(j) ) EXIT
                IF ( proc==procs(l) ) THEN
                  Found=.TRUE.; EXIT
                END IF
              END DO
              IF (.NOT.Found) CYCLE
            END IF
          ELSE
            l=k
          END IF
          toReceive(lproc) % perm(j)=ldofs(l)
        END DO
      END DO
!------------------------------------------------------------------------------
   END SUBROUTINE FetiSendReceiveInit
!------------------------------------------------------------------------------


  !> The neighbour send/receive communication, after initializations.
  ! ----------------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE FetiSendReceive(toSend, toReceive, tag, Fsum)
!------------------------------------------------------------------------------
      REAL(KIND=dp), OPTIONAL :: Fsum(:)
      INTEGER :: tag
      TYPE(toSend_t) :: toSend(:)  
      TYPE(toReceive_t) :: toReceive(:)  
!------------------------------------------------------------------------------
      INTEGER :: i,j,k,l,n,m,proc, lproc
      LOGICAL :: Found
      REAL(KIND=dp), ALLOCATABLE :: buf(:)
!------------------------------------------------------------------------------
      DO i=1,nneigh
        proc = gpnum(i)
        CALL FetiSend(proc, toSend(i) % n, toSend(i) % buf, tag=tag)
      END DO 

      ! Receive interface parts and sum values
      ! --------------------------------------
      n = 0
      DO i=1,nneigh
        n=MAX(n,toReceive(i) % n)
      END DO
      ALLOCATE(buf(n))

      DO i=1,nneigh
        CALL FetiRecv(proc, n, buf, tag=tag)
        lproc=lpnum(proc)

        IF (.NOT.PRESENT(Fsum)) THEN
          IF(.NOT.ALLOCATED(toReceive(lproc) % buf)) THEN
            ALLOCATE(toReceive(lproc) % buf(Bmat % NumberOfRows))
          END IF
          toReceive(lproc) % buf=0._dp
        END IF

        DO j=1,n
          l=toReceive(lproc) % perm(j)
          IF (l>0) THEN
            IF (PRESENT(Fsum)) THEN
              Fsum(l)=Fsum(l)+buf(j)
            ELSE
              toReceive(lproc) % buf(l)=buf(j)
            END IF
          END IF
        END DO
      END DO
!------------------------------------------------------------------------------
   END SUBROUTINE FetiSendReceive
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    !
    !> Given local (owned L.C.'s) interface vector b, send and receive
    !> all partition interface L.C's to/from neighbours. Place result
    !> to partitionwise vector f. In effect f=B^Tb.
    ! ---------------------------------------------------------------
 SUBROUTINE FetiSendRecvLC(A,f,b)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Matrix_t) :: A
    REAL(KIND=dp) :: f(:),b(:)
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,p,q,m,n,sz,nrows,proc,lproc
    INTEGER :: sndLC, ownLC
    LOGICAL :: Found
    INTEGER, POINTER :: gtags(:)
    LOGICAL, POINTER :: ig(:)
    TYPE(NeighbourList_t), POINTER :: nb(:)
    INTEGER, ALLOCATABLE :: gdofs(:), ldofs(:), procs(:)

    REAL(KIND=dp) :: val

    TYPE(toReceive_t), ALLOCATABLE, SAVE :: toReceive(:)
    TYPE(toSend_t), ALLOCATABLE, SAVE :: toSend(:)
    INTEGER, SAVE :: ninterface
    INTEGER, ALLOCATABLE, SAVE :: lint(:)
!------------------------------------------------------------------------------

    nrows = A % NumberOfRows
    gtags => A % ParallelInfo % GlobalDofs
    ig => A % ParallelInfo % GInterface
    nb => A % ParallelInfo % NeighbourList

    IF ( InitializeLC .OR. .NOT. ALLOCATED(toSend)) THEN

      IF ( ALLOCATED(toSend) ) THEN
        DO i=1,SIZE(toSend)
          DEALLOCATE(toSend(i) % buf, toSend(i) % ifg, &
           toSend(i) % perm, toReceive(i) % perm)
        END DO
        DEALLOCATE(toSend,toReceive,lint)
      END IF
 
      ! Count sizes of send & receive buffers:
      ! --------------------------------------
      ALLOCATE(toSend(nneigh),toReceive(nneigh))

      ninterface=COUNT(ig)
      ALLOCATE(lint(ninterface))

      ninterface = 0
      sndLC=0
      toSend(:) % n=0
      DO i=1,nrows
        ! Dirichlet DOFs are left out of the interface compatibility conditions
        ! because in the standard scheme they have already been eliminated from
        ! the local matrices, so there is nothing left to make continuous. That
        ! reasoning does not hold for total FETI: there the Dirichlet conditions
        ! are imposed by Lagrange coefficients instead, the DOFs stay in the
        ! local systems as ordinary unknowns, and their copies in neighbouring
        ! partitions have to be tied together like any other interface DOF.
        ! See the identical test in FetiSendRecvIf.
        IF (ig(i).AND.(TotalFeti.OR..NOT.DirichletDOFs(i))) THEN
          ninterface=ninterface+1
          lint(ninterface)=i
          sz = SIZE(nb(i) % Neighbours)
          proc=nb(i) % Neighbours(1)
          IF (proc == ParEnv % myPE ) THEN
            DO j=2,sz
              lproc=lpnum(nb(i) % Neighbours(j))
              toSend(lproc) % n = toSend(lproc) % n+1
            END DO
          ELSE
            sndLC = sndLC+1
          END IF
        END IF
      END DO

      ! Allocate send & receive buffers:
      ! --------------------------------------
      DO i=1,nneigh
        j=toSend(i) % n
        ALLOCATE(toSend(i) % buf(j),toSend(i) % ifg(j), &
             toSend(i) % perm(ninterface))
        toSend(i) % perm=0
      END DO
      ALLOCATE(gdofs(sndLC), ldofs(sndLC), procs(sndLC))

      ! Extract send & receive dof tags:
      ! --------------------------------
      sndLC=0
      toSend(:) % n=0
      DO i=1,ninterface
        l = lint(i)
        proc = nb(l) % Neighbours(1)
        IF (proc == ParEnv % myPE ) THEN
          sz = SIZE(nb(l) % Neighbours)
          DO j=2,sz
            lproc=lpnum(nb(l) % Neighbours(j))
            k = toSend(lproc) % n+1
            toSend(lproc) % perm(i) = k
            toSend(lproc) % n = k
            toSend(lproc) % ifg(k) = gtags(l)
          END DO
        ELSE
          sndLC = sndLC+1
          ldofs(sndLC) = l
          procs(sndLC) = proc
          gdofs(sndLC) = gtags(l)
        END IF
      END DO

      ! Send interface parts to neighbours, at initialization
      ! only store local indices of global tags
      ! ------------------------------------------------------
      CALL FetiSendReceiveInit( sndLC, gdofs, ldofs, procs, &
                toSend, toReceive, tag=100 )

      DEALLOCATE(gdofs,ldofs,procs)
      InitializeLC = .FALSE.
    END IF

    ! Extract send & receive dof values of f=B^Tb
    ! -------------------------------------------
    F = 0._dp
    DO lproc=1,nneigh
      toSend(lproc) % buf = 0._dp
    END DO

    DO i=1,Bmat % NumberOfRows
      l = Bmat % Perm(i)

      IF (l>0) THEN  ! partition DOF number
        m = lint(l)  ! on interface
      ELSE
        m = -l       ! internal D.B.C.
      END IF

      DO j=Bmat % Rows(i),Bmat % Rows(i+1)-1
        val = Bmat % Values(j) * b(i)
        IF ( Bmat % Cols(j) == ParEnv % MyPE ) THEN
          F(m) = F(m) + val
        ELSE IF ( l<=SIZE(lint)) THEN
          ! fill send buffer position of LC for DOF l (or m):
          ! -------------------------------------------------
          lproc = lpnum(Bmat % Cols(j))
          k = toSend(lproc) % perm(l)
          toSend(lproc) % buf(k) = toSend(lproc) % buf(k)+val
        END IF
      END DO
    END DO

    ! Send & reiceive interface parts to neighbours
    ! ------------------------------------
    CALL FetiSendReceive(toSend, toReceive, tag=110, Fsum=F)


!------------------------------------------------------------------------------
  END SUBROUTINE FetiSendRecvLC
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Add Dirichlet BCs to set of Lagrange coefficients described by the 'B' matrix
!> (the total FETI scheme).
!
! Some thoughts:
!
! Trouble here is that in Elmer the Dirichlet BCs can be really varied and may,
! for example, include a (large) number of domain internal DOFs. Also, at first
! thought, P-elements with nontrivial Dirichlet BCs don't fit easily into this
! scheme as not all the coefficients are known from the input to data directly.
! At least you'd need multidof constraints...
!
! Only most simple cases handled here currently, should call DefaultDirichletBCs()
! eventually to take care of more details. Even then the Robin type of conditions
! would not be handled. Actually seems quite hard to ensure all 'floating' domains
! at all times.
!------------------------------------------------------------------------------
  SUBROUTINE FetiAddDtoB(A, B, nLC)
!------------------------------------------------------------------------------
    INTEGER :: nLC
    TYPe(Matrix_t) :: A
    TYPe(Matrix_t), POINTER :: B
!------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: Solver
    INTEGER :: i,j,k,l,n,d,Active
    INTEGER, ALLOCATABLE :: Perm(:), InEqualityFlag(:)
    INTEGER, POINTER  :: p(:)
    LOGICAL :: Found,FoundUpper,FoundLower
    LOGICAL, ALLOCATABLE :: Done(:)
    REAL(KIND=dp) :: val,scale
    REAL(KIND=dp), ALLOCATABLE :: Vals(:)
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    Solver => GetSolver()

    n = A % NumberOfRows

    ! nLC coefficients already exist and at most one is appended per local DOF,
    ! the Done() flags below making sure of that, so nLC+n is a bound that holds
    ! whatever the partitioning. Sizing these by n alone overflowed all three as
    ! soon as the interface coefficients outnumbered the rows.
    ALLOCATE(Perm(nLC+n),InEqualityFlag(nLC+n), B % RHS(nLC+n))
    B % RHS=0._dp
    Perm(1:nLC) = B % Perm(1:nLC)
    InEqualityFlag(1:nLC) = 0

    IF (FetiAsPrec) THEN
      DO i=1,n
        IF ( DirichletDofs(i) ) THEN
          nLC = nLC+1
          Perm(nLC) = -i
          B % RHS(nLC) = A % RHS(i)
          CALL SetMatrixElement(B, nLC, ParEnv % myPE, 1._dp)
        END IF
      END DO
    ELSE 
      d =  Solver % Variable % DOFs
      p => Solver % Variable % Perm

      ALLOCATE(Vals(Solver % Mesh % MaxElementNodes), Done(n))
      Done = .FALSE.

      Active = GetNOFBoundaryElements()
      DO i=1,Active
        Element => GetBoundaryElement(i)
        BC => GetBC()
        
        IF (.NOT.ASSOCIATED(BC)) CYCLE
        IF (.NOT. ActiveBoundaryElement()) CYCLE

        n = GetElementNOFNodes()
        DO j=1,d
          IF (d>1) THEN
            Vals(1:n)=GetReal(BC,ComponentName(Solver % Variable,j),Found)
          ELSE
            Vals(1:n)=GetReal(BC,Solver % Variable % Name,Found)
          END IF

          FoundUpper=.FALSE.
          FoundLower=.FALSE.
          IF(.NOT.Found) THEN
            IF (d>1) THEN
              Vals(1:n)=GetReal(BC,TRIM(ComponentName(Solver % Variable,j))//' Upper Limit',FoundUpper)
            ELSE
              Vals(1:n)=GetReal(BC,TRIM(Solver % Variable % Name)//' Upper Limit',FoundUpper)
            END IF

            IF(.NOT.FoundUpper) THEN
              IF (d>1) THEN
                Vals(1:n)=GetReal(BC,TRIM(ComponentName(Solver % Variable,j))//' Lower Limit',FoundLower)
              ELSE
                Vals(1:n)=GetReal(BC,TRIM(Solver % Variable % Name)//' Lower Limit',FoundLower)
              END IF
            END IF
          END IF

          IF (.NOT.(Found.OR.FoundUpper.OR.FoundLower)) CYCLE

          DO k=1,n
            l = p(Element % NodeIndexes(k))
            IF (l<=0) CYCLE 

            l=d*(l-1)+j
            IF(Done(l)) CYCLE

            nLC = nLC+1
            Done(l)=.TRUE.

            ! - sign here to indicate D-condition instead of interface
            ! compatibility conditions...:
            ! --------------------------------------------------------
            Perm(nLC)=-l
            IF(FoundUpper.OR.FoundLower) InEqualityFlag(nLC)=1
            scale=1._dp
            IF(FoundLower) THEN
              scale = -1._dp
              Vals(k)=-Vals(k)
            END IF

            IF (ASSOCIATED(A % DiagScaling)) THEN
              B % RHS(nLC) = Vals(k)/A % DiagScaling(l)
              CALL SetMatrixElement(B, nLC, ParEnv % myPE, scale)
            ELSE
              B % RHS(nLC) = Vals(k)
              CALL SetMatrixElement(B, nLC, ParEnv % myPE, scale)
            END IF
          END DO
        END DO
      END DO
    END IF

    DEALLOCATE(B % Perm)
    ALLOCATE(B % Perm(nLC), B % InvPerm(nLC))
    B % Perm = Perm(1:nLC)
    B % InvPerm = InEqualityFlag(1:nLC)
!------------------------------------------------------------------------------
  END SUBROUTINE FetiAddDtoB
!------------------------------------------------------------------------------


  !------------------------------------------------------------------------------
  !> Extract interface values from partition vector b; send and
  !> receive  to/from neighbours.  Only 'owned' interface dofs
  !> placed to result vector f. At initialization also assemble
  !> the 'B' connectivity matrix. In effect f=Bb;
  !------------------------------------------------------------------------------
  FUNCTION FetiSendRecvIf(A,f,b,g,l_i) RESULT(nLC)
  ! ----------------------------------------------------------
    IMPLICIT NONE
    TYPE(Matrix_t) :: A
    INTEGER :: nLC
    REAL(KIND=dp) :: f(:),b(:)

    integer,optional :: l_i(0:)
    real(kind=dp),optional :: g(:,:)
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,m,p,q,n,sz,nrows,proc,lproc,totnLC=0
    TYPE(NeighbourList_t), POINTER :: nb(:)
    LOGICAL :: Found
    INTEGER, POINTER :: gtags(:)
    LOGICAL, POINTER :: ig(:)
    INTEGER, ALLOCATABLE :: gdofs(:), ldofs(:), procs(:)

    REAL(KIND=dp) :: val
    REAL(KIND=dp), POINTER :: Fn(:)

    TYPE(toSend_t), ALLOCATABLE, SAVE :: toSend(:)
    INTEGER, SAVE :: ninterface
    INTEGER, ALLOCATABLE, SAVE :: lint(:)
    TYPE(toReceive_t), ALLOCATABLE, TARGET, SAVE :: toReceive(:)
!------------------------------------------------------------------------------

    nrows = A % NumberOfRows
    gtags => A % ParallelInfo % GlobalDofs
    ig => A % ParallelInfo % GInterface
    nb => A % ParallelInfo % NeighbourList

    IF (InitializeIf .OR. .NOT. ALLOCATED(toSend) ) THEN

      IF ( ALLOCATED(toSend) )THEN
        DO i=1,SIZE(toSend)
          DEALLOCATE(toSend(i) % buf, toSend(i) % ifg, toReceive(i) % perm)
          IF (ALLOCATED(toReceive(i) % buf)) DEALLOCATE(toReceive(i) % buf)
        END DO
        DEALLOCATE(toSend,toReceive,lint)
      END IF
      ALLOCATE(toSend(nneigh),toReceive(nneigh))

      ninterface=COUNT(ig)
      ALLOCATE(lint(ninterface))

      ! Count sizes of send & receive buffers:
      ! --------------------------------------
      ninterface=0
      nLC=0
      toSend(:) % n=0
      DO i=1,nrows
        ! Under total FETI the Dirichlet DOFs are not eliminated from the local
        ! matrices, so they need the interface compatibility conditions too.
        ! See the longer comment at the identical test in FetiSendReceiveInit.
        IF (ig(i) .AND.(TotalFeti.OR..NOT.DirichletDOFs(i))) THEN
          ninterface=ninterface+1
          lint(ninterface)=i

          proc=nb(i) % Neighbours(1)
          IF (proc == ParEnv % myPE ) THEN
            nLC=nLC+SIZE(nb(i) % Neighbours)-1
          ELSE
            lproc=lpnum(proc)
            toSend(lproc) % n = toSend(lproc) % n+1
          END IF
        END IF
      END DO

      ! Allocate send & receive buffers:
      ! --------------------------------------
      DO i=1,nneigh
        j=toSend(i) % n
        ALLOCATE(toSend(i) % buf(j), toSend(i) % ifg(j) )
      END DO
      ALLOCATE( gdofs(nLC), ldofs(nLC), procs(nLC) )
      ldofs = [(i,i=1,nLC)];

      ! Extract send & receive dof tags; Allocate, initialize
      ! and assemble the  'B' connectivity matrix:
      ! -----------------------------------------------------

      IF (ASSOCIATED(Bmat)) CALL FreeMatrix(Bmat)
      Bmat => AllocateMatrix()
      Bmat %  ListMatrix => Null()
      ALLOCATE(Bmat % Perm(nLC))
      Bmat % Format = MATRIX_LIST

      nLC=0
      toSend(:) % n=0
      DO i=1,ninterface
        l = lint(i)
        proc  = nb(l) % Neighbours(1)
        IF (proc == ParEnv % myPE) THEN
          sz = SIZE(nb(l) % Neighbours)
          DO j=1,sz-1
            nLC = nLC+1
            gdofs(nLC) = gtags(l)
            procs(nLC) = nb(l) % Neighbours(j+1)

            ! Assemble the 'B' matrix, orthogonal connectivity for
            ! multiple sharings:
            ! u_2 = u_1,
            ! u_3 = (u_1+u_2)/2,
            ! u_4 = (u_1+u_2+u_3)/3
            !         ....
            ! (also normalized below...)
            ! -----------------------------------------------------
            Bmat % Perm(nLC) = i
            val = 1._dp/(j*SQRT(1+1._dp/j))
            DO m=1,j
              proc = nb(l) % Neighbours(m)
              CALL SetMatrixElement(Bmat,nLC,proc,val)
            END DO
            proc = nb(l) % Neighbours(j+1)
            CALL SetMatrixElement(Bmat,nLC,proc,-j*val)

!           proc = nb(l) % Neighbours(1)
!           CALL SetMatrixElement(Bmat,nLC,proc,1._dp)
!           proc = nb(l) % Neighbours(j+1)
!           CALL SetMatrixElement(Bmat,nLC,proc,-1._dp)

          END DO
        ELSE
          lproc = lpnum(proc)
          k = toSend(lproc) % n+1
          toSend(lproc) % n = k
          toSend(lproc) % ifg(k) = gtags(l)
        END IF
      END DO

      ! Send interface parts to neighbours
      ! ------------------------------------
      CALL FetiSendReceiveInit( nLC, gdofs, ldofs, procs, &
               toSend, toReceive, tag=200 )

      ! If 'Total' FETI add D-conditions to 'B':
      ! ----------------------------------------
      totnLC=nLC
      IF (TotalFeti) CALL FetiAddDtoB(A, Bmat, totNLC)

      ! Convert 'B' to CRS format:
      ! --------------------------
      CALL List_ToCRSMatrix(Bmat)


      DEALLOCATE(gdofs,ldofs,procs)
      InitializeIf = .FALSE.
      IF(dumptofiles) THEN
        CALL SaveB(A,ninterface,lint)
        RETURN
      END IF
    END IF
!------------------------------------------------------------------------------


    ! Extract send & receive dof values: 
    ! ----------------------------------
    toSend(:) % n=0
    DO i=1,ninterface
      l = lint(i)
      proc = nb(l) % Neighbours(1)
      IF (proc /= ParEnv % myPE ) THEN
        lproc=lpnum(proc)
        k = toSend(lproc) % n+1
        toSend(lproc) % n = k
        toSend(lproc) % buf(k) = b(l)
      END IF
    END DO

    ! Send & reiceive interface parts to neighbours
    ! ---------------------------------------------
    CALL FetiSendReceive(toSend, toReceive, tag=210)

    F = 0._dp

    ! Local contribution to f=Bb
    ! --------------------------
    DO i=1,Bmat % NumberOfRows
      l = Bmat % Perm(i)
      IF (l>0) THEN ! m = Partition DOF number...
        m = lint(l) ! ...for an interface DOF
      ELSE
        m = -l      ! ...for a Dirichlet DOF
      END IF

      DO j=Bmat % Rows(i), Bmat % Rows(i+1)-1
        IF (Bmat % Cols(j)==ParEnv % myPE) THEN
          if (present(g)) then
            ! store the B_me*R_me for G'G assembly,
            ! if requested:
            ! --------------------------------------
            if ( l_i(0)>0 ) &
              g(i,l_i(0))=g(i,l_i(0))+Bmat%values(j)*b(m)
          end if
          F(i) = F(i) + Bmat % Values(j)*b(m); EXIT
        END IF
      END DO
    END DO


    ! Neighbours contribution to f=Bb (this surely
    ! could be simplified....):
    ! ---------------------------------------------
    DO lproc=1,nneigh
      Fn => toReceive(lproc) % buf
      nLC=0
      DO i=1,ninterface
        l = lint(i)
        proc = nb(l) % Neighbours(1)
        IF (proc == ParEnv % myPE ) THEN
          q=nLC
          sz=SIZE(nb(l) % Neighbours)
          DO j=1,sz-1
            nLC = nLC+1
            IF (Fn(nLC)==0) CYCLE
            DO m=1,sz-1
              DO k=Bmat % Rows(q+m),Bmat % Rows(q+m+1)-1
                IF (Bmat % Cols(k)==gpnum(lproc)) THEN
                  if(present(g)) then
                    ! store the B_lproc*R_lproc for G'G assembly,
                    ! if requested:
                    ! --------------------------------------------
                    if ( l_i(lproc)>0 ) &
                      g(q+m,l_i(lproc))=g(q+m,l_i(lproc))+Bmat%values(k)*Fn(nLC)
                  end if
                  F(q+m) = F(q+m) + Bmat % Values(k)*Fn(nLC); EXIT
                END IF
              END DO
            END DO
          END DO
        END IF
      END DO
    END DO
    nLC = totnLC
!------------------------------------------------------------------------------
  END FUNCTION FetiSendRecvIf
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  !> The projection of the vector T to the set of feasible solutions
  !> orthogonal to the null space of A, defined by the condition:
  !> z^T(f-B^T\lambda) (OP=0, 'z' a.k.a. 'R'). Also handle the initial
  !> guess and the final extraction of the L.C.'s \alpha related to the
  !> above condition (OP=1,2).
!------------------------------------------------------------------------------
  SUBROUTINE FetiProject(A,n,T,OP,TOL)
!------------------------------------------------------------------------------
    INTEGER :: n, &
              OP  !=0: T =  (I-G*Ginv*G')T, 
                  !=1: T = -(G*Ginv*G') T,  note: input size:  nz
                  !=2: T =  (Ginv*G')T,     note: output size: nz
                  ! Ginv = (G'*G)^-1
    REAL(KIND=dp) :: TOL
    REAL(KIND=dp)::T(n)
    TYPE(matrix_t) :: A
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: S(n),err1
    INTEGER :: i,j,k,l,m,nlc,nrows
    REAL(KIND=dp),ALLOCATABLE :: x(:),b(:),P(:),Q(:)
!------------------------------------------------------------------------------
    integer :: ierr,stat(mpi_status_size),l_n,g_n, &
          n_nbr,comm,i_type,d_type,maxnz,proc
    integer, save :: me=0, mactive=0
    logical :: iterative, found
    integer,allocatable :: l_nz(:),g_nz(:),l_i(:),g_i(:),l_nbr(:)
    real(kind=dp), allocatable :: l_g(:,:),l_gtg(:,:),g_x(:),g_b(:)

    type(matrix_t), pointer :: gtg_m => null()

    integer,save :: size, mygroup, grpsize, subsize, comm_group,solv_group,solv_comm,ir=0
    integer, allocatable, save :: ranks(:,:),subsizes(:)

    save :: g_nz,g_i,gtg_m,g_n,g_x,g_b
!------------------------------------------------------------------------------
!!call resettimer('project')

    maxnz = parallelreduction(nz,2)

    ! check whether anything to do:
    ! -----------------------------
    if (maxnz<=0) then
      if (op==2) t=0
      return
    end if

    CALL ListPushNameSpace('feti proj:')

    nrows = A % NumberOfRows
    ALLOCATE(x(nz),b(nz),P(nrows),Q(nrows));P=0; Q=0

    iterative=getlogical(getsolverparams(), 'Feti CPG Projection Iterative', found)

    if (.not.iterative) then

      comm = a % comm
      i_type = mpi_integer
      d_type = mpi_double_precision

      if (refactorize .and. op==1) then
!call resettimer('factorize setup')

        do mactive=0,parenv % pes-1
          if ( parenv % active(mactive+1) ) exit
        end do
        me = parenv % mype

        ! If not using CG, assemble and factorize the G'G matrix:
        !
        ! NOTE for further development: Instead of using one task to do the
        ! assembling, factorizing & solving we could use parallel Cholesky
        ! factorization from SCALAPACK or, for example, MUMPS. Should not be
        ! such a big deal. This might be favorable, even necessary, when using
        ! (tens of) thousands of mpi tasks ...
        ! --------------------------------------------------------------------

        ! -------------------------------------------------------

        ! First get the size of the ker(A_i)'s from the neighbours:
        ! =========================================================
        do i=1,nneigh
          call mpi_bsend(nz,1,i_type,gpnum(i),400,comm,ierr)
        end do

        allocate(l_nz(0:nneigh))
        do i=1,nneigh
          call mpi_recv(l_nz(i),1,i_type,gpnum(i),400,comm,stat,ierr)
        end do
        l_nz(0)=nz
        l_n=sum(l_nz)

        ! Get the 'G' matrix rows corresponding to our set of interface
        ! conditions (multiply by identity, and assemble the columns):
        ! ================================================================

        ! local sequential numbering for the \alpha's,
        ! based on position in the neighbour tables:
        ! ---------------------------------------------
        if (allocated(g_i)) deallocate(g_i)
        allocate(l_i(0:nneigh+1),g_i(0:nneigh+1))
        l_i(0)=1
        do i=1,nneigh+1
          l_i(i)=l_i(i-1)+l_nz(i-1)
        end do
        g_i=l_i  ! g_i used as a temporary below

        ! do the assembly of the (local) 'G' matrix:
        ! ------------------------------------------
        allocate(l_g(n,l_n)); l_g=0
        if (nz>0) x=0
        do i=1,maxnz
          ! have we already consumed some set of \alpha's
          ! (own or neighbours)?
          ! ---------------------------------------------
          do j=0,nneigh
            if (g_i(j)<=0.or.g_i(j)>=l_i(j+1)) g_i(j)=-1
          end do

          if (i<=nz) x(i)=1
          if (nz>0) p=matmul(x,z)
          if (i<=nz) x(i)=0

          nLC = FetiSendRecvIf(A,S,P,l_g,g_i)
          g_i = g_i+1

          ! Seems necessary to slow some guys down. Note that
          ! everybody needs to iterate the same amount due to
          ! this, hence the overall max for the loopcount
          ! above. Could check if this is really necessary...
          ! well, maybe someone, someday....
          ! --------------------------------------------------
          call parallelactivebarrier()
        end do
        deallocate(l_i,l_nz,g_i) !

        ! form local part of the G'G:
        ! ============================
        allocate(l_gtg(l_n,l_n));l_gtg=0
        l_gtg=matmul(transpose(l_g),l_g)
        deallocate(l_g)

        size = parenv % pes

        grpsize=listgetinteger( params, 'Feti Projection Solution Group Size', &
                          found,minv=1,maxv=size)
        if (.not.found) grpsize=1

        if(allocated(ranks)) deallocate(ranks,subsizes)
        allocate(ranks(grpsize,size/grpsize+1),subsizes(grpsize))
        ranks=0;subsizes=0

        ir=0
        do while(ir<size)
          do j=1,grpsize
            ir=ir+1
            if(ir>size) exit
            subsizes(j)=subsizes(j)+1
          end do
        end do

        ir = 0
        do i=1,grpsize
          do j=1,subsizes(i)
            ranks(i,j)=ir
            ir = ir+1
          end do
        end do

        mygroup=grpsize
        do i=1,grpsize-1
          if (me<ranks(i+1,1)) then
            mygroup = i
            exit
           end if
        end do
        mactive = ranks(mygroup,1)
        subsize = subsizes(mygroup)

        call mpi_comm_group(comm, comm_group, ierr) ! Extract the original group handle
        call mpi_group_incl(comm_group,grpsize,ranks(:,1),solv_group,ierr)
        call mpi_comm_create(comm,solv_group,solv_comm,ierr)


        ! form global G'G:
        ! =================

        do i=1,grpsize
         if (ranks(i,1)/=me) &
            call mpi_bsend(nz,1,i_type,ranks(i,1),401,comm,ierr)
        end do

        if (me/=mactive) then

          ! send local G'G to 'mactive' task:
          ! ----------------------------------
          call mpi_bsend(nneigh,1,i_type,mactive,402,comm,ierr)
          call mpi_bsend(gpnum,nneigh,i_type,mactive,403,comm,ierr)
          call mpi_bsend(l_gtg,l_n**2,d_type,mactive,404,comm,ierr)
          deallocate(l_gtg)

        else ! me==mactive

          ! Assemble whole of the G'G from the parts:
          ! -----------------------------------------

          ! get size of all the ker(A_i)'s:
          ! -------------------------------
          if (allocated(g_nz)) deallocate(g_nz)
          allocate(g_nz(0:parenv % pes-1)); g_nz=0
          g_nz(me)=nz
          do i=0,parenv % pes-1
            if ( i==me .or. .not. parenv % active(i+1)) cycle
            call mpi_recv(g_nz(i),1,i_type,i,401,comm,stat,ierr)
          end do

          ! global sequential numbering for the \alpha's,
          ! based on task id's:
          ! ---------------------------------------------
          if(allocated(g_i)) deallocate(g_i)
          allocate(g_i(0:parenv % pes))
          g_i(0)=1
          do j=1,parenv % pes
            g_i(j)=g_i(j-1)+g_nz(j-1)
          end do
          g_n=sum(g_nz)

          ! Get storage for the (global) G'G:
          ! ---------------------------------
          IF (ASSOCIATED(gtg_m)) CALL FreeMatrix(gtg_m)
          gtg_m=>AllocateMatrix()
          gtg_m % Comm = solv_comm
          gtg_m % format = MATRIX_LIST

          ! Assemble our own part...
          ! ------------------------

          allocate(l_nbr(0:nneigh))
          n_nbr=nneigh
          l_nbr(0)=me; l_nbr(1:n_nbr)=gpnum(1:n_nbr)

          call addtogtg(gtg_m,l_gtg,n_nbr,l_nbr,g_nz,g_i)
          deallocate(l_gtg,l_nbr)

          ! ...and the rest of it:
          ! ----------------------
          do i=1,subsize-1 ! -1 to count 'me' out
            call mpi_recv(n_nbr,1,i_type,mpi_any_source,402,comm,stat,ierr)
            proc=stat(mpi_source)

            allocate(l_nbr(0:n_nbr))
            l_nbr(0)=proc
            call mpi_recv(l_nbr(1:),n_nbr,i_type,proc,403,comm,stat,ierr)

            l_n=sum(g_nz(l_nbr))
            allocate(l_gtg(l_n,l_n))
            call mpi_recv(l_gtg,l_n**2,d_type,proc,404,comm,stat,ierr)

            call addtogtg(gtg_m,l_gtg,n_nbr,l_nbr,g_nz,g_i)
            deallocate(l_gtg,l_nbr)
          end do

          !
          ! Factorize the G'G, and we are done initializing:
          ! ------------------------------------------------
          call list_tocrsmatrix(gtg_m)

          g_n=gtg_m % numberofrows
          allocate(g_x(g_n),g_b(g_n)); g_x=0; g_b=0

        end if ! task 'mactive' doing the factoring (and later, solving)
!call checktimer('factorize setup',delete=.true.)
      end if ! the first time G'G assembly and factorizing
    end if ! direct solver setup

    IF (OP==1) THEN
      IF (nz>0) b=T(1:nz)
    ELSE
      CALL Gt(b,T)
    END IF

    !
    ! Solve x from G'Gx = b:
    ! -----------------------
    if (iterative) then
      IF (nz>0) x=b
      CALL FCG(nz,x,b)
    else

!call resettimer('coarse prob')
      if (me==mactive) then

        ! assemble the r.h.s. vector:
        ! ---------------------------

        if (nz>0) g_b(g_i(me):g_i(me+1)-1) = b
        do i=2,subsize
          j=ranks(mygroup,i)
          if (g_nz(j)>0) &
            call mpi_recv(g_b(g_i(j):g_i(j+1)-1), g_nz(j), &
                    d_type,j,405,comm,stat,ierr)
        end do

        ! solve x:
        ! --------
        call directsolver(gtg_m,g_x,g_b,getsolver())

        ! distribute the result:
        ! ----------------------
        do i=2,subsize
          j=ranks(mygroup,i)
          if (g_nz(j)>0) &
            call mpi_bsend(g_x(g_i(j):g_i(j+1)-1), &
                 g_nz(j),d_type,j,406,comm,ierr)
        end do

        ! finally my own part:
        ! --------------------
        if (nz>0) x(1:nz)=g_x(g_i(me):g_i(me+1)-1)
      else
        ! send r.h.s.:
        ! ------------
         if (nz>0) &
           call mpi_bsend(b,nz,d_type,mactive,405,comm,ierr)

        ! receive result:
        ! --------------
         if (nz>0) &
           call mpi_recv(x,nz,d_type,mactive,406,comm,stat,ierr)
      endif
      call parallelactivebarrier()
    endif
!call checktimer('coarse prob',delete=.true.)

    SELECT CASE(OP)
    CASE(0)
      CALL G(x,S)
      T=T-S
    CASE(1)
      CALL G(x,S)
      T=-S
    CASE(2)
      T(1:nz)=x
    CASE DEFAULT
      CALL Fatal('Feti Procjection','Unknown projection OP.')
    END SELECT

    CALL ListPopNameSpace()
!call checktimer('project',delete=.true.)

CONTAINS

!------------------------------------------------------------------------------
  subroutine addtogtg(gtg_m,l_gtg,n_nbr,l_nbr,g_nz,g_i)
!------------------------------------------------------------------------------
    real(kind=dp) :: l_gtg(:,:)
    type(matrix_t), pointer :: gtg_m
    integer :: n_nbr, l_nbr(0:), g_nz(0:), g_i(0:)
!------------------------------------------------------------------------------
    integer :: i,j,k,l,g_r,g_c,l_r,l_c,l_i(0:n_nbr+1)
!------------------------------------------------------------------------------
    l_i(0)=1
    do i=1,n_nbr+1
      l_i(i)=l_i(i-1)+g_nz(l_nbr(i-1))
    end do

    do i=0,n_nbr
      do k=0,g_nz(l_nbr(i))-1
        l_c=l_i(i)+k
        g_c=g_i(l_nbr(i))+k
        do j=0,n_nbr
          do l=0,g_nz(l_nbr(j))-1
            l_r=l_i(j)+l
            g_r=g_i(l_nbr(j))+l
            if (abs(l_gtg(l_r,l_c))>aeps) &
              call addtomatrixelement(gtg_m,g_r,g_c,l_gtg(l_r,l_c))
          end do
        end do
      end do
    end do
!------------------------------------------------------------------------------
  end subroutine addtogtg
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE GtG(u,v)
!------------------------------------------------------------------------------
    INTEGER :: nlc
    REAL(KIND=dp) :: u(:),v(:)

    IF (nz>0) P=MATMUL(u,z)
    nlc = FetiSendRecvIf(A,S,P)

    CALL FetiSendRecvLC(A,Q,S)
    IF (nz>0) v = MATMUL(z,Q)
!------------------------------------------------------------------------------
  END SUBROUTINE GtG
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE Gt(v,u)
!------------------------------------------------------------------------------
    REAL(KIND=dp) ::v(:),u(:)

    CALL FetiSendRecvLC(A,Q,u)
    IF (nz>0) v = MATMUL(z,Q)
!------------------------------------------------------------------------------
  END SUBROUTINE Gt
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE G(u,v)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: u(:),v(:)
    INTEGER :: n

    IF (nz>0) P=MATMUL(u,z)
    n = FetiSendRecvIf(A,v,P)
!------------------------------------------------------------------------------
  END SUBROUTINE G
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE FCG(n,x,b)
!------------------------------------------------------------------------------
    INTEGER :: n
    REAL(KIND=dp) :: x(:),b(:)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: beta,alpha,rho,prevrho,bnorm,err
    INTEGER :: iter
    REAL(KIND=dp) :: Ri(n),S(n),T(n)
!------------------------------------------------------------------------------
    bnorm = SparNorm(n,b,1)
    IF ( bnorm==0 ) THEN
      x=0; RETURN;
    END IF

    CALL GtG(x,T)
    Ri = -(T - b)

    err = SparNorm(n,Ri,1)
    IF (err<TOL) RETURN

    beta=0._dp
    DO iter=1,500
      rho = SparDotProd(n,Ri,1,Ri,1)
      IF (iter==1) THEN
        S = Ri
      ELSE
        beta = rho/prevrho
        S = Ri + beta*S
      END IF
      prevrho = rho

      CALL GtG(S,T)
      alpha = rho/SparDotProd(n,S,1,T,1)
      x = x + alpha*S

      Ri = Ri - alpha*T
      err = SparNorm(n,Ri,1)
      IF (err<TOL) EXIT
    END DO
    IF ( Parenv % MyPE==0 .AND.  err>=TOL ) print*,'fcg not converged',err
!------------------------------------------------------------------------------
  END SUBROUTINE FCG
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE FetiProject
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The C(onjugate) P(rojected) G(radient) iterative solver.
!------------------------------------------------------------------------------
  SUBROUTINE FetiCPG(A,n,x,b,f,Solver,matvecsubr,precsubr,dotprodfun,normfun)
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t),pointer :: A
    INTEGER :: n
    REAL(KIND=dp) :: x(:),b(:), f(:),dotprodfun, normfun

    EXTERNAL matvecsubr, precsubr, dotprodfun, normfun
!------------------------------------------------------------------------------
    INTEGER :: i,j,m,iter,nLC,nrows,maxit,ipar(50),output,Restart,Saven
    REAL(KIND=dp) :: beta,alpha,rho,prevrho,bnorm,err0,err1,TOL
    LOGICAL :: Found
    REAL(KIND=dp), ALLOCATABLE :: Ri(:),S(:),T(:),P(:), &
           Ssave(:,:), Psave(:,:), srho(:)
!------------------------------------------------------------------------------
!call resettimer('cpg')
    nrows=A % NumberOfRows

    Params => ListGetSolverParams()
    output = GetInteger( Params,'Linear System Residual Output', Found )
    IF (.NOT. Found ) output = 1
    maxit   = GetInteger( Params,'Linear System Max Iterations')
    Restart = GetInteger( Params,'Linear System CPG Restart', Found)
    TOL     = GetConstReal( Params,'Linear System Convergence Tolerance')

    HUTI_NDIM=n
    bnorm = normfun(n,b,1)
    IF ( bnorm==0 ) THEN
      x=0; RETURN;
    END IF

    ALLOCATE(T(n),S(n),Ri(n))

    x=0
    IF (nz>0) THEN
      m = SIZE(z,2)
      x(1:nz)=MATMUL(z,f(1:m))
    END IF
    CALL FetiProject(A,n,x,OP=1,TOL=1d-12)

    CALL matvecsubr(x(1:n),T,ipar)
    Ri = -(T - b(1:n))
    CALL FetiProject(A,n,Ri,OP=0,TOL=1d-12)

    err0 = normfun(n,Ri,1)
    err1 = err0 / bnorm
    IF ( output /= 0 .AND. ParEnv % MyPE==0 ) THEN
      PRINT*,' '
      PRINT*,'         iter     |ax-b|                    |ax-b|/b';
      PRINT*,'      --------------------------------------------------';flush(6)
      PRINT*,0,err0,err1; Flush(6)
    END IF
    IF (err1<TOL) RETURN

    IF (Restart>0) THEN
      Saven=0
      ALLOCATE(Ssave(n,Restart), Psave(n,Restart), P(n), srho(Restart))
    END IF

    beta=0._dp
    DO iter=1,maxit
!call resettimer('iter')
       IF (Restart>0) P=T

       IF (Precondition) THEN
         CALL precsubr(T,Ri,ipar)
         CALL FetiProject(A,n,T,OP=0,TOL=1d-12)
       ELSE
         T=Ri
       END  IF

       rho = dotprodfun(n,Ri,1,T,1)
       IF ( iter==1 ) THEN
         S = T
       ELSE
         IF (Restart>0) THEN
           IF (Saven>=Restart) Saven=0
           saven=saven+1
           Ssave(:,Saven)=S
           Psave(:,Saven)=P
           srho(saven)=dotprodfun(n,S,1,P,1)
           S=T
           DO i=1,Saven
             P=Psave(:,i)
             beta=dotprodfun(n,T,1,P,1)/srho(i)
             S=S-beta*Ssave(:,i)
           END DO
         ELSE
           beta = rho/prevrho
           S = T + beta*S
         END IF
       END IF
       prevrho = rho

       CALL matvecsubr(S,T,ipar)
       alpha = rho / dotprodfun(n,S,1,T,1)
       x(1:n) = x(1:n) + alpha*S

       CALL FetiProject(A,n,T,OP=0,TOL=1d-12)
       Ri = Ri - alpha*T

       err0 = normfun(n,Ri,1)
       err1 = err0 / bnorm

!err2 = FetiStopc(x,b,Ri,ipar,dpar)
!CALL matvecsubr(x,rr,ipar)
!rr = rr - b
!CALL FetiProject(A,n,rr,OP=0,TOL=1d-12)
!err2 = normfun(n,rr,1)

       IF ( output>0 .AND. MOD(iter,output)==0) THEN
         IF (ParEnv % MyPE==0) THEN
           PRINT*,iter,err0,err1; Flush(6)
         END IF
       END IF

       IF (err1<TOL) EXIT
!call checktimer('iter',delete=.true.)
    END DO

    ! Falling out of the loop without reaching the tolerance used to pass
    ! silently, and the coefficients were then used as if they had converged.
    ! A stagnating solve is thereby returned as a plausible looking but wrong
    ! answer, which is a bad way to find out about it.
    IF (err1>=TOL) THEN
      WRITE(Message,'(A,I0,A,ES12.5,A,ES12.5)') 'CPG did not converge in ', &
          maxit,' iterations: residual ',err1,' against a tolerance of ',TOL
      IF (ListGetLogical(Params,'Linear System Abort Not Converged',Found)) THEN
        CALL Fatal('FetiCPG',Message)
      ELSE
        CALL Warn('FetiCPG',Message)
      END IF
    END IF

    CALL matvecsubr(x(1:n),T,ipar)
    T = -(T-b(1:n))
    CALL FetiProject(A,n,T,OP=2,TOL=1d-12)
    IF (nz>0) x(n+1:n+nz) = T(1:nz)
!call checktimer('cpg',delete=.true.)
!------------------------------------------------------------------------------
  END SUBROUTINE FetiCPG
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Split the local matrix graph into connected components, returning their
!> number and labelling every row in comp().
!>
!> Graph partitioners hand out subdomains that fall into several disconnected
!> pieces quite readily, and a subdomain in c pieces has c times the kernel of a
!> single one: c constants for Poisson, 6c rigid body modes for elasticity. The
!> analytic kernels built below have to be split accordingly, or the local
!> matrix is left singular in all but one piece.
!>
!> A is symmetric here (FETI requires it), so its sparsity is an undirected
!> graph and a plain flood fill over the numerically nonzero entries is enough.
!------------------------------------------------------------------------------
  FUNCTION FetiConnectedComponents(A,comp) RESULT(nc)
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A
    INTEGER :: comp(:)
    INTEGER :: nc
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,n,top
    INTEGER, ALLOCATABLE :: stack(:)
!------------------------------------------------------------------------------
    n = A % NumberOfRows
    nc = 0
    ALLOCATE(stack(n))

    ! Rows with no off-diagonal coupling are not pieces of the domain. Dirichlet
    ! DOFs eliminated to an identity row look exactly like that, and there can be
    ! thousands of them: counting each as its own component would report nonsense
    ! and size the per component work by it. They carry no kernel either, since
    ! their diagonal is nonzero, so a kernel vector has to vanish there anyway.
    ! Mark them -1 for now and hand them back as 0, i.e. belonging to no piece.
    ! --------------------------------------------------------------------------
    DO i=1,n
      comp(i) = -1
      DO j = A % Rows(i), A % Rows(i+1)-1
        IF (A % Cols(j) == i) CYCLE
        IF (A % Values(j) == 0._dp) CYCLE
        comp(i) = 0
        EXIT
      END DO
    END DO

    DO i=1,n
      IF (comp(i) /= 0) CYCLE
      nc = nc+1
      comp(i) = nc
      top = 1; stack(1) = i
      DO WHILE(top > 0)
        k = stack(top); top = top-1
        DO j = A % Rows(k), A % Rows(k+1)-1
          IF (A % Values(j) == 0._dp) CYCLE
          IF (comp(A % Cols(j)) /= 0) CYCLE
          comp(A % Cols(j)) = nc
          top = top+1
          stack(top) = A % Cols(j)
        END DO
      END DO
    END DO

    WHERE(comp < 0) comp = 0

    DEALLOCATE(stack)
!------------------------------------------------------------------------------
  END FUNCTION FetiConnectedComponents
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Infinity norm of the local matrix, used to make kernel tests scale
!> invariant. For a symmetric A this bounds ||A||_2 from above.
!------------------------------------------------------------------------------
  FUNCTION FetiMatrixNorm(A) RESULT(Anorm)
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) :: Anorm
!------------------------------------------------------------------------------
    INTEGER :: i
    REAL(KIND=dp) :: rowsum
!------------------------------------------------------------------------------
    Anorm = 0._dp
    DO i=1,A % NumberOfRows
      rowsum = SUM(ABS(A % Values(A % Rows(i):A % Rows(i+1)-1)))
      Anorm = MAX(Anorm, rowsum)
    END DO
    IF (Anorm <= 0._dp) Anorm = 1._dp
!------------------------------------------------------------------------------
  END FUNCTION FetiMatrixNorm
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Can the kernel be removed by pinning the DOFs in FixInds?
!>
!> Pinning works exactly when the nz x nz matrix z(:,FixInds) is nonsingular: a
!> kernel vector surviving the pinning is one that vanishes at every pinned DOF,
!> and there is such a vector precisely when that matrix is singular. Being
!> nearly singular is nearly as bad, since the pinned matrix is then nearly
!> singular too and the local solves lose accuracy quietly.
!>
!> The choice was never checked. FindRigidBodyFixingNodes picks the nodes on
!> geometric grounds, which can land on a degenerate set -- collinear nodes for a
!> rotation, say -- and nothing noticed. z is orthonormal by the time this runs,
!> so the singular values are directly comparable to one.
!------------------------------------------------------------------------------
  FUNCTION FetiFixingUsable() RESULT(Usable)
!------------------------------------------------------------------------------
    LOGICAL :: Usable
!------------------------------------------------------------------------------
    INTEGER :: i,j,ierr
    REAL(KIND=dp) :: tol
    REAL(KIND=dp), ALLOCATABLE :: M(:,:),sv(:),wrk(:),dummy(:,:)
    LOGICAL :: Found
!------------------------------------------------------------------------------
    Usable = .TRUE.
    IF (nz <= 0) RETURN

    ALLOCATE(M(nz,nz),sv(nz),wrk(MAX(1,10*nz)),dummy(1,1))
    DO j=1,nz
      DO i=1,nz
        M(i,j) = z(i,FixInds(j))
      END DO
    END DO

    CALL DGESVD('N','N',nz,nz,M,nz,sv,dummy,1,dummy,1,wrk,SIZE(wrk),ierr)

    IF (ierr /= 0) THEN
      CALL Warn('Feti:','Could not check the DOFs chosen for fixing')
      DEALLOCATE(M,sv,wrk,dummy)
      RETURN
    END IF

    Params => ListGetSolverParams()
    tol = GetCReal(Params,'Feti Fixing Tolerance',Found)
    IF (.NOT.Found) tol = 1.0d-6

    Usable = sv(nz) > tol

    WRITE(Message,'(A,ES10.3)') &
        'Smallest singular value of the fixing DOFs: ',sv(nz)
    CALL Info('Feti:',Message,Level=7)

    IF (.NOT.Usable) THEN
      WRITE(Message,'(A,ES10.3,A)') 'The DOFs chosen for fixing do not remove '// &
          'the kernel (smallest singular value ',sv(nz), &
          '); using Lagrange coefficients instead'
      CALL Warn('Feti:',Message)
    END IF

    DEALLOCATE(M,sv,wrk,dummy)
!------------------------------------------------------------------------------
  END FUNCTION FetiFixingUsable
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Orthonormalize the kernel basis z=null(A) and check that it really is one.
!>
!> Whichever route produced z -- the analytic rigid body modes, the null
!> frequency eigenmodes, SPQR or the Mumps null pivots -- the vectors come back
!> in no particular scaling and are not mutually orthogonal. That matters
!> downstream: the coarse operator is G=Bz and the coarse problem is solved
!> from G^T G, which squares the condition number, so a poorly scaled or nearly
!> dependent z shows up as a stalling coarse solve rather than as anything that
!> points back here. Only the span of z enters the FETI formulation, so
!> replacing it with an orthonormal basis of the same span is free of side
!> effects and strictly better conditioned.
!>
!> The residual check afterwards costs nz matrix-vector products and turns an
!> otherwise opaque rank decision (Mumps CNTL(3), an eigenvalue threshold) into
!> something that can be seen in the log.
!------------------------------------------------------------------------------
  SUBROUTINE FetiCheckKernel(A,Origin)
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A
    CHARACTER(LEN=*) :: Origin
!------------------------------------------------------------------------------
    INTEGER :: i,j,n,ncomp
    INTEGER, ALLOCATABLE :: comp(:)
    REAL(KIND=dp) :: nrm,nrm0,proj,Anorm,res,maxres,tol
    REAL(KIND=dp), ALLOCATABLE :: x(:),w(:)
    LOGICAL :: Found, Deficient
!------------------------------------------------------------------------------
    IF (nz <= 0) RETURN

    n = A % NumberOfRows
    Params => ListGetSolverParams()

    tol = GetCReal(Params,'Feti Kernel Tolerance',Found)
    IF (.NOT.Found) tol = 1.0d-9

    ALLOCATE(x(n),w(n))

    ! Modified Gram-Schmidt. A vector that (nearly) vanishes here was already
    ! dependent on the previous ones, i.e. the reported deficiency is too
    ! large. We keep the count rather than silently changing nz, since FixInds
    ! was sized against it, but say so loudly.
    Deficient = .FALSE.
    DO i=1,nz
      w = z(i,1:n)
      nrm0 = SQRT(SUM(w**2))
      DO j=1,i-1
        proj = SUM(z(j,1:n)*z(i,1:n))
        z(i,1:n) = z(i,1:n) - proj*z(j,1:n)
      END DO
      nrm = SQRT(SUM(z(i,1:n)**2))
      IF (nrm0 <= 0._dp .OR. nrm <= SQRT(EPSILON(1._dp))*nrm0) THEN
        ! Restore the unorthogonalized vector so that everything stays finite.
        Deficient = .TRUE.
        z(i,1:n) = w
        nrm = nrm0
      END IF
      IF (nrm > 0._dp) z(i,1:n) = z(i,1:n)/nrm
    END DO

    IF (Deficient) THEN
      WRITE(Message,'(A,A,A)') 'Kernel basis from ',TRIM(Origin), &
          ' is linearly dependent, the deficiency is likely overestimated'
      CALL Warn('Feti:',Message)
    END IF

    ! Now that z is orthonormal, ||A z_i|| relative to ||A|| is a meaningful
    ! measure of how well z_i really is in the kernel.
    Anorm = FetiMatrixNorm(A)
    maxres = 0._dp
    DO i=1,nz
      CALL MatrixVectorMultiply(A,z(i,1:n),x)
      res = SQRT(SUM(x**2))/Anorm
      maxres = MAX(maxres,res)
    END DO

    ! Report the number of disconnected pieces alongside. The routes that read
    ! the kernel off a factorization get the dimension right without knowing
    ! why it is what it is, and "6 modes over 1 piece" versus "6 modes over 3
    ! pieces, so 12 are missing" is the difference worth seeing.
    ALLOCATE(comp(n))
    ncomp = FetiConnectedComponents(A,comp)
    DEALLOCATE(comp)

    WRITE(Message,'(A,A,A,I0,A,I0,A,ES10.3)') 'Kernel from ',TRIM(Origin), &
        ': dim ',nz,' over ',ncomp,' piece(s), max relative residual ',maxres
    CALL Info('Feti:',Message,Level=6)

    IF (ncomp > 1) THEN
      WRITE(Message,'(A,I0,A)') 'Local matrix falls into ',ncomp, &
          ' disconnected pieces; check the kernel dimension accounts for all of them'
      CALL Warn('Feti:',Message)
    END IF

    ! A loose factor over the detection tolerance: we are flagging a rank
    ! decision that looks wrong, not measuring round-off.
    IF (maxres > 1000._dp*tol) THEN
      WRITE(Message,'(A,ES10.3,A,A)') 'Suspect null space, ||Az||/||A|| = ', &
          maxres,' for a vector reported to be in the kernel by ',TRIM(Origin)
      CALL Warn('Feti:',Message)
    END IF

    DEALLOCATE(x,w)
!------------------------------------------------------------------------------
  END SUBROUTINE FetiCheckKernel
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Compute null(A), return value is whether null(A) is nonempty.
!------------------------------------------------------------------------------
  FUNCTION FetiFloatingDomain(A,Solver,FixInds,TOL) RESULT(Floating)
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) :: TOL
    LOGICAL :: floating
    INTEGER :: FixInds(:)
    TYPE(Solver_t) :: Solver
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE :: x(:),tz(:,:)
    REAL(KIND=dp), ALLOCATABLE :: S(:,:),Seig(:),Swrk(:),Az(:,:),tz2(:,:),y2(:)
    INTEGER, POINTER :: p(:)
    INTEGER, ALLOCATABLE :: comp(:)
    INTEGER :: i,j,k,n,m,dofs,dim,FixNodes(0:6)
    INTEGER :: m2,ierr,nc,ic,kc,nfree,ncand
    INTEGER, ALLOCATABLE :: cand(:)
    REAL(KIND=dp), ALLOCATABLE :: W(:,:)
    LOGICAL :: Usable
    REAL(KIND=dp), POINTER :: coord_x(:),coord_y(:),coord_z(:), dscale(:)
    LOGICAL :: Found
    REAL(KIND=dp) :: xc,yc,zc,hc,Anorm,KerTol
    INTEGER, ALLOCATABLE :: floatinds(:)
!------------------------------------------------------------------------------

    Params => ListGetSolverParams()

    dofs = Solver % Variable % DOFs
    n = A % NumberOfRows

    IF(ALLOCATED(z)) DEALLOCATE(z)
    nz = 0

    ! We might be operating in a scaled system, i.e. solving
    ! (DAD)y = Db; x=Dy
    ! instead of the original system Ax=b, hence the scalings below:
    ! --------------------------------------------------------------
    IF (ASSOCIATED(A % DiagScaling)) THEN
      dscale => A % DiagScaling
    ELSE
      ALLOCATE(dScale(n)); dscale=1._dp
    END IF

    dim = CoordinateSystemDimension()

    IF(dofs==1) THEN
      !
      ! For Poisson, check if 1 in null space:
      ! --------------------------------------

      m =  Solver % Mesh % NumberOfNodes
      p => Solver % Variable % Perm

      ! One constant per connected component of the local matrix, not one for
      ! the subdomain as a whole. Testing a single global constant does not
      ! reveal the difference: the sum of the per component constants is in the
      ! kernel as well, so a subdomain in several pieces passed the check with
      ! the deficiency reported as 1 and its local matrix left singular.
      !
      ! The value is put on the nodal DOFs only, as before, since in a
      ! hierarchical basis the constant has no higher order coefficients.
      ! ------------------------------------------------------------------------
      ALLOCATE(comp(n))
      nc = FetiConnectedComponents(A,comp)
      ALLOCATE(tz2(nc,n),x(n))
      tz2 = 0._dp
      DO i=1,m
        j=p(i)
        ! comp==0 are the rows belonging to no piece, see FetiConnectedComponents
        IF (j>0) THEN
          IF (comp(j)>0) tz2(comp(j),j) = 1._dp / dscale(j)
        END IF
      END DO

      Anorm = FetiMatrixNorm(A)
      KerTol = GetCReal(Params,'Feti Kernel Tolerance',Found)
      IF (.NOT.Found) KerTol = 1.0d-9

      ! Keep the components that really are free. A component carrying a
      ! Dirichlet condition is not, and drops out here.
      ! ------------------------------------------------------------------------
      nz=0
      DO i=1,nc
        hc = SQRT(SUM(tz2(i,:)**2))
        IF (hc <= 0._dp) CYCLE
        tz2(i,:) = tz2(i,:)/hc
        CALL MatrixVectorMultiply(A,tz2(i,:),x)
        IF (SQRT(SUM(x**2))/Anorm > KerTol) CYCLE
        nz=nz+1
        IF (nz /= i) tz2(nz,:) = tz2(i,:)
      END DO

      Floating = nz>0
      NullSpaceLC=GetLogical(Params,'Feti Fixing Using L.C.',Found)

      IF (Floating) THEN
        ALLOCATE(z(nz,n))
        z = tz2(1:nz,:)

        ! One DOF to fix per free component, taken where that component's
        ! constant is largest so that the pinned DOF certainly carries it.
        ! With a single component this is the node FindRigidBodyFixingNodes
        ! would have picked in all but name; with several it is the only way
        ! to remove the kernel by pinning at all.
        ! ----------------------------------------------------------------------
        DO i=1,nz
          FixInds(i) = MAXLOC(ABS(z(i,:)),DIM=1)
        END DO
      END IF

      IF (nc > 1) THEN
        WRITE(Message,'(A,I0,A,I0,A)') 'Local matrix falls into ',nc, &
            ' disconnected pieces, ',nz,' of them floating'
        CALL Info('Feti:',Message,Level=5)
      END IF

      IF (.NOT.ASSOCIATED(A % DiagScaling)) DEALLOCATE(dscale)
      RETURN
    ELSE
      IF (.NOT.GetLogical(Params,'Feti Kernel Rot-Trans',Found)) THEN
        IF (Found) CALL Warn('Feti:','"Feti Kernel Rot-Trans" is no longer '// &
            'optional: the null frequency eigenmode alternative it selected '// &
            'against has been removed. Use "mumpslocal" with "Mumps Solve '// &
            'Singular", or "spqr", to read the kernel off the factorization.')
      END IF
      !
      ! Otherwise null(A) is spanned by the rigid body modes, i.e. the basic
      ! translations and rotations, plus the potential levels. This used to sit
      ! behind "Feti Kernel Rot-Trans" with a null frequency eigenmode
      ! computation as the default. That default is gone: it ran Arpack for up
      ! to twenty modes and a direct factorization per subdomain, then threw the
      ! factorization away, and decided the numerical rank from shift inverted
      ! eigenvalues, which is the least reliable way to ask that question. It
      ! was also silently capped, so a deficiency above the cap came back
      ! truncated. The modes below are exact for the operators FETI is used on
      ! here, cost a handful of matrix-vector products, and are checked against
      ! the matrix afterwards; when they are not the right answer, "mumpslocal"
      ! with "Mumps Solve Singular", or "spqr", reads the kernel straight off
      ! the factorization.
      ! ------------------------------------------------------------------

      m =  Solver % Mesh % NumberOfNodes
      p => Solver % Variable % Perm

      coord_x => Solver % Mesh % Nodes % x
      coord_y => Solver % Mesh % Nodes % y
      coord_z => Solver % Mesh % Nodes % z

      ALLOCATE(comp(n))
      nc = FetiConnectedComponents(A,comp)

      IF (dim==2) m2=dofs+1
      IF (dim==3) m2=dofs+dim

      Anorm = FetiMatrixNorm(A)
      KerTol = GetCReal(Params,'Feti Kernel Tolerance',Found)
      IF (.NOT.Found) KerTol = 1.0d-9

      ALLOCATE(tz(m2,n), x(n), y2(n), floatinds(m2), S(m2,m2), Seig(m2), &
          Swrk(MAX(1,3*m2)), Az(m2,n), tz2(m2*nc,n))

      nz = 0
      nfree = 0
      DO ic=1,nc

        ! Centre and length scale of this piece. Taken per component rather than
        ! over the whole subdomain: a rotation measured about a centroid that
        ! lies in a different, disconnected piece is badly scaled.
        ! ----------------------------------------------------------------------
        xc=0._dp; yc=0._dp; zc=0._dp; kc=0
        DO i=1,m
          IF (p(i) <= 0) CYCLE
          IF (comp(dofs*(p(i)-1)+1) /= ic) CYCLE
          kc=kc+1
          xc=xc+coord_x(i); yc=yc+coord_y(i); zc=zc+coord_z(i)
        END DO
        IF (kc == 0) CYCLE
        xc=xc/kc; yc=yc/kc; zc=zc/kc

        hc = 0._dp
        DO i=1,m
          IF (p(i) <= 0) CYCLE
          IF (comp(dofs*(p(i)-1)+1) /= ic) CYCLE
          hc = MAX(hc,(coord_x(i)-xc)**2+(coord_y(i)-yc)**2+(coord_z(i)-zc)**2)
        END DO
        hc = SQRT(hc)
        IF (hc <= 0._dp) hc = 1._dp

        ! Translations, potential levels and rotations, restricted to this piece
        ! ----------------------------------------------------------------------
        tz = 0._dp
        DO i=1,m
          IF (p(i) <= 0) CYCLE
          IF (comp(dofs*(p(i)-1)+1) /= ic) CYCLE
          DO j=1,dofs
            k=dofs*(p(i)-1)+j
            tz(j,k) = 1._dp / dscale(k)
          END DO
          j=dofs*(p(i)-1)
          tz(dofs+1,j+1) = -(coord_y(i)-yc) / hc / dscale(j+1)  ! rot about z
          tz(dofs+1,j+2) =  (coord_x(i)-xc) / hc / dscale(j+2)
          IF (dim>=3) THEN
            tz(dofs+2,j+2) = -(coord_z(i)-zc) / hc / dscale(j+2) ! rot about x
            tz(dofs+2,j+3) =  (coord_y(i)-yc) / hc / dscale(j+3)

            tz(dofs+3,j+1) =  (coord_z(i)-zc) / hc / dscale(j+1) ! rot about y
            tz(dofs+3,j+3) = -(coord_x(i)-xc) / hc / dscale(j+3)
          END IF
        END DO

        DO i=1,m2
          hc = SQRT(SUM(tz(i,:)**2))
          IF( hc > 0._dp ) tz(i,:) = tz(i,:) / hc
          CALL MatrixVectorMultiply(A,tz(i,:),Az(i,:))
        END DO
        DO i=1,m2
          DO j=1,m2
            S(i,j) = SUM(tz(i,:)*Az(j,:))
          END DO
        END DO
        S = 0.5_dp*(S + TRANSPOSE(S))   ! symmetrize away the round-off

        CALL DSYEV('V','U',m2,S,m2,Seig,Swrk,SIZE(Swrk),ierr)
        IF (ierr /= 0) CALL Fatal('FetiFloatingDomain', &
            'Eigenvalue decomposition of the candidate kernel failed')

        ! Accept a combination when its residual is small relative to ||A||. The
        ! Rayleigh quotient Seig alone would be a more permissive test, since it
        ! behaves like the square of the distance to the kernel, so use it only
        ! to order the candidates and let the residual make the decision.
        ! ----------------------------------------------------------------------
        DO i=1,m2
          IF (Seig(i) > KerTol*Anorm) CYCLE
          x = 0._dp
          DO j=1,m2
            x = x + S(j,i)*tz(j,:)
          END DO
          hc = SQRT(SUM(x**2))
          IF (hc <= 0._dp) CYCLE
          x = x/hc
          CALL MatrixVectorMultiply(A,x,y2)
          IF (SQRT(SUM(y2**2))/Anorm > KerTol) CYCLE
          nz = nz+1
          tz2(nz,:) = x
        END DO

        ! How many of this piece's modes are in the kernel on their own. Only
        ! those can be expressed by pinning a single DOF.
        ! ----------------------------------------------------------------------
        DO i=1,m2
          IF (SQRT(SUM(Az(i,:)**2))/Anorm > KerTol) CYCLE
          nfree = nfree+1
          IF (nc == 1) FloatInds(nfree) = i
        END DO
      END DO

      IF (.NOT.ASSOCIATED(A % DiagScaling)) DEALLOCATE(dscale)

      ! "nfree" counted the modes that are in the kernel by themselves. Only
      ! those can be expressed by pinning a single DOF, so if the kernel is
      ! larger the rest are genuine combinations and pinning cannot remove them.
      ! Several pieces rule pinning out too: FindRigidBodyFixingNodes returns one
      ! set of nodes for the subdomain, which cannot fix a mode living in another
      ! piece. Either way, fall back to the Lagrange coefficient route.
      ! --------------------------------------------------------------------------
      j = nfree

      Floating=nz>0

      ! The kernel itself, however it ends up being fixed below.
      IF (Floating) THEN
        ALLOCATE(z(nz,n))
        z = tz2(1:nz,:)
      END IF

      NullSpaceLC=GetLogical(Params,'Feti Fixing Using L.C.',Found)

      IF (nc > 1) THEN
        WRITE(Message,'(A,I0,A,I0,A)') 'Local matrix falls into ',nc, &
            ' disconnected pieces, kernel dimension ',nz, &
            '; fixing with Lagrange coefficients instead of DOFs'
        CALL Info('Feti:',Message,Level=5)
        NullSpaceLC = .TRUE.
      ELSE IF (Floating .AND. j < nz) THEN
        WRITE(Message,'(A,I0,A,I0,A)') 'Kernel dimension ',nz, &
            ' exceeds the ',j,' free rigid body modes, so it contains a '// &
            'combination; fixing with Lagrange coefficients instead of DOFs'
        CALL Info('Feti:',Message,Level=5)
        NullSpaceLC = .TRUE.
      END IF

      ! Set the DOFs to fix at FixNodes() nodes.
      !
      ! Pinning removes the kernel exactly when z(:,FixInds) is nonsingular, so
      ! choose the DOFs by pivoting on that matrix rather than by mapping each
      ! mode to a node by hand. Every component of every node the search offers
      ! is a candidate, and Gaussian elimination with column pivoting picks the
      ! nz of them that keep it furthest from singular. That is the same
      ! criterion FetiFixingUsable checks afterwards, so the two agree by
      ! construction.
      !
      ! The mapping this replaces read FixNodes(0), which
      ! FindRigidBodyFixingNodes never writes: it documents its output as
      ! x1 x2 y1 y2 z1 z2 and fills entries 1..2*dim only. All three translations
      ! were therefore pinned at whatever that entry happened to hold, which is
      ! how the fixing matrix ended up singular to 1e-17 on most partitions of a
      ! finely divided mesh.
      ! ------------------------------------------------------------------------
      IF (Floating .AND. nc == 1) THEN
        CALL FindRigidBodyFixingNodes(Solver, FixNodes, p)

        ncand = 0
        ALLOCATE(cand(2*dim*dofs), W(nz,2*dim*dofs))
        DO k=1,2*dim
          IF (FixNodes(k) <= 0) CYCLE
          IF (p(FixNodes(k)) <= 0) CYCLE
          DO i=1,dofs
            ncand = ncand+1
            cand(ncand) = dofs*(p(FixNodes(k))-1)+i
          END DO
        END DO

        Usable = ncand >= nz
        IF (Usable) THEN
          W(:,1:ncand) = z(1:nz,cand(1:ncand))
          DO i=1,nz
            ! Pick the column with the largest remaining entry in row i.
            k = 0; hc = 0._dp
            DO m2=i,ncand
              IF (ABS(W(i,m2)) > hc) THEN
                hc = ABS(W(i,m2)); k = m2
              END IF
            END DO
            IF (hc <= EPSILON(1._dp)) THEN
              Usable = .FALSE.; EXIT
            END IF
            ! Swap it into place and eliminate below.
            IF (k /= i) THEN
              y2(1:nz) = W(:,i); W(:,i) = W(:,k); W(:,k) = y2(1:nz)
              m2 = cand(i); cand(i) = cand(k); cand(k) = m2
            END IF
            FixInds(i) = cand(i)
            DO m2=i+1,ncand
              W(:,m2) = W(:,m2) - (W(i,m2)/W(i,i))*W(:,i)
            END DO
          END DO
        END IF

        IF (.NOT.Usable) THEN
          CALL Info('Feti:','No set of DOFs among the fixing nodes removes '// &
              'the kernel; using Lagrange coefficients instead',Level=5)
          NullSpaceLC = .TRUE.
        END IF
        DEALLOCATE(cand,W)
      END IF
      RETURN
    END IF


!------------------------------------------------------------------------------
  END FUNCTION FetiFloatingDomain
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solve x from [A z^T; z 0] [x; \lambda] = [b 0]
!>  or A x  = b, depending on how the system is set up.
!------------------------------------------------------------------------------
  SUBROUTINE FetiDirectSolver(A,x,b,Solver)
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: a
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp), TARGET CONTIG :: x(:),b(:)
!------------------------------------------------------------------------------
    INTEGER :: n,i
    REAL(KIND=dp), POINTER CONTIG :: tx(:),tb(:)
    REAL(KIND=dp), POINTER CONTIG :: bfix(:)
!call resettimer('direct')
    n = A % NumberOfRows
    tx=>x
    tb=>b

    IF (FixingByDofs.AND.nz>0) THEN
      ! A has had the FixInds rows and columns replaced by identity rows to take
      ! the kernel out. That only stays a generalized inverse of the original
      ! matrix if the right hand side vanishes there too: otherwise the pinned
      ! DOF comes back as b(FixInds) instead of zero and every local solve is off
      ! by that much. The pinning is done once when the matrix is set up, while
      ! the right hand side changes at every iteration, so it has to be imposed
      ! here rather than there.
      ! b is the caller's right hand side, which is sized for the Lagrange
      ! coefficients and so longer than the n rows solved for here.
      ALLOCATE(bfix(n))
      bfix = b(1:n)
      DO i=1,nz
        bfix(FixInds(i)) = 0._dp
      END DO
      tb => bfix
    END IF

    IF (NullSpaceLC.AND.nz>0) THEN
      ALLOCATE(tx(n+nz),tb(n+nz))
      tb=0
      tb(1:n)=b(1:n)
      A % NumberOfRows=n+nz
    END IF

    CALL DirectSolver(A,tx,tb,Solver)

    ! Check the first solve after the matrix was set up. The kernel handed to it
    ! has been verified to lie in null(A), but nothing has checked that it is all
    ! of null(A), and there is no cheap way to. It is not, when a partition holds
    ! pieces joined at a single node or a single edge: those are one piece to a
    ! graph but a mechanism to the operator, with three or one extra modes each
    ! beyond the six rigid body ones. Whatever is left of the kernel then survives
    ! into the matrix solved here, and the result is not a solution of anything.
    ! Left alone it comes out as a residual of 1e18 and then NaN some hundred
    ! iterations later, which says nothing about the cause.
    IF (.NOT.LocalSolveChecked) THEN
      LocalSolveChecked = .TRUE.
      BLOCK
        REAL(KIND=dp), ALLOCATABLE :: res(:)
        REAL(KIND=dp) :: rn,bn
        INTEGER :: ii,jj,nn
        nn = A % NumberOfRows
        ALLOCATE(res(nn)); res = 0._dp
        DO ii=1,nn
          DO jj=A % Rows(ii),A % Rows(ii+1)-1
            res(ii) = res(ii) + A % Values(jj)*tx(A % Cols(jj))
          END DO
        END DO
        res = res - tb(1:nn)
        rn = SQRT(SUM(res**2)); bn = SQRT(SUM(tb(1:nn)**2))
        DEALLOCATE(res)
        IF (bn > 0._dp .AND. rn > 1.0d-6*bn) THEN
          WRITE(Message,'(A,ES12.5,A,I0,A)') &
              'The local matrix was not solved: relative residual ',rn/bn, &
              '. Its kernel is most likely bigger than the ',nz, &
              ' modes found for it, which happens when a partition holds pieces '// &
              'joined at a single node or edge. Use "mumpslocal" with "Mumps '// &
              'Solve Singular", or "spqr", which read the kernel off the '// &
              'factorization instead of assuming it.'
          CALL Fatal('FetiDirectSolver',Message)
        END IF
      END BLOCK
    END IF

    IF (NullSpaceLC.AND.nz>0) THEN
      A % NumberOfRows=n
      x=tx(1:n)
      DEALLOCATE(tx,tb)
    ELSE IF (FixingByDofs.AND.nz>0) THEN
      DEALLOCATE(bfix)
    END IF
!call checktimer('direct',delete=.true.)
!------------------------------------------------------------------------------
  END SUBROUTINE FetiDirectSolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE SaveKandF(A)
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A

    INTEGER :: me, abeg,i,j
    INTEGER, ALLOCATABLE :: snd(:), asize(:), bsize(:)

    me = Parenv % MyPE
    ALLOCATE(snd(0:Parenv%PEs-1),asize(0:Parenv%PEs-1),bsize(0:parenv%pes-1))
    snd=0
    snd(me)=A % NumberOfRows
    CALL MPI_ALLREDUCE( snd,asize,Parenv % PEs,MPI_INTEGER,MPI_MAX,&
                              ELMER_COMM_WORLD,ierr)
    abeg = SUM(asize(0:me-1))

    OPEN(1,file='f' // i2s(Parenv % MyPE))
    OPEN(2,file='k' // i2s(Parenv % MyPE))

    WRITE(1,'(a)') '% domain: '//i2s(me)//' nrows:'//i2s(A % NumberOFRows)

    WRITE(2,'(a)') '% domain:' // i2s(me)//' nnz:' // &
       i2s(A % Rows(A % NumberOfRows+1)-1) // ' nrows:' // &
          i2s(A % NumberOFRows) // ' gcols:'//i2s(SUM(asize))

    DO i=1,A % NumberOfRows
      DO j=A % Rows(i),A % Rows(i+1)-1
        WRITE(2,*) abeg+i, abeg+A % Cols(j), A % Values(j)
      END DO
      WRITE(1,*) abeg+i, a % RHS(i)
    END DO
    CLOSE(2)
    CLOSE(1)
!------------------------------------------------------------------------------
  END SUBROUTINE SaveKandF
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE SaveR
!------------------------------------------------------------------------------
    INTEGER :: i
    OPEN(2,File='r'//i2s(Parenv % MyPE))
    WRITE(2,'(a)') '% domain: '//i2s(ParEnv % MyPE)//' nz:'// &
              i2s(SIZE(z,1))//' nrows:'// i2s(SIZE(z,2))
    DO i=1,SIZE(z,2)
      WRITE(2,*) z(1:nz,i)
    END DO
    CLOSE(2)
!------------------------------------------------------------------------------
  END SUBROUTINE SaveR
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE SaveB(A,ninterface,lint)
    TYPE(Matrix_t) :: A
    INTEGER :: lint(:)
    INTEGER :: ninterface
!------------------------------------------------------------------------------
    INTEGER, ALLOCATABLE :: snd(:), asize(:), bsize(:), cnt(:),ibuf(:,:),gbuf(:,:)

    INTEGER :: bbeg, i,j,k,l,proc,n,m,me

    INTEGER, POINTER :: gtags(:)
    LOGICAL, POINTER :: ig(:)
    TYPE(NeighbourList_t), POINTER :: nb(:)

    gtags => A % ParallelInfo % GlobalDofs
    ig => A % ParallelInfo % GInterface
    nb => A % ParallelInfo % NeighbourList

    OPEN(4,FILE='b'//I2S(Parenv % MyPE))
    OPEN(5,FILE='brhs'//I2S(Parenv % MyPE))

    me = ParEnv % MyPE

    WRITE(4,'(a)') '% domain: '//i2s(me)//' nnz: '// &
         i2s(bMat % Rows(Bmat % NumberOfRows+1)-1) // &
               ' nrows: ' // i2s(Bmat % NumberOfRows)

    WRITE(5,'(a)') '% domain: '//i2s(ParEnv % MyPE)//' nrows:'// &
               i2s(Bmat % NumberOfRows)

    ALLOCATE(snd(0:Parenv%PEs-1),asize(0:Parenv%PEs-1),bsize(0:parenv%pes-1))

    snd=0
    snd(me)=A % NumberOfRows
    CALL MPI_ALLREDUCE( snd,asize,Parenv % PEs,MPI_INTEGER,MPI_MAX,&
                              ELMER_COMM_WORLD,ierr)
    abeg = SUM(asize(0:me-1))

    snd=0
    snd(me)=Bmat % NumberOFrows
    CALL MPI_ALLREDUCE( snd,bsize,ParEnv%PEs,MPI_INTEGER, &
               MPI_MAX, ELMER_COMM_WORLD,ierr)
    bbeg=sum(bsize(0:me-1))

    ALLOCATE(cnt(0:parenv%pes-1))
    ALLOCATE(ibuf(ninterface,0:parenv%pes-1))
    ALLOCATE(gbuf(ninterface,0:parenv%pes-1))

    cnt=0
    DO i=1,nInterface
      m = lint(i)
      DO j=1,SIZE(nb(m) % Neighbours)
        proc=nb(m) % Neighbours(j)
        IF(proc==me) CYCLE
        cnt(proc) = cnt(proc)+1
        ibuf(cnt(proc),proc)=m
        gbuf(cnt(proc),proc)=gtags(m)
      END DO
    END DO

    DO i=1,ParEnv % PEs
      IF(.NOT. ParEnv % ISneighbour(i)) CYCLE
        proc=i-1
      CALL MPI_BSEND( cnt(proc), 1, MPI_INTEGER, proc, &
                800, ELMER_COMM_WORLD, ierr )
      IF(cnt(proc)>0) THEN
        CALL MPI_BSEND( gbuf(:,proc), cnt(proc), MPI_INTEGER, proc, &
                  801, ELMER_COMM_WORLD, ierr )
  
        CALL MPI_BSEND( ibuf(:,proc), cnt(proc), MPI_INTEGER, proc, &
                  802, ELMER_COMM_WORLD, ierr )
      END IF
    END DO
    CALL ParallelActiveBarrier()

    DO i=1,ParEnv % PEs
      IF(.NOT. ParEnv % IsNeighbour(i)) CYCLE
      proc=i-1
      CALL MPI_RECV( cnt(proc), 1, MPI_INTEGER, proc, &
          800, ELMER_COMM_WORLD, status, ierr )

      IF(cnt(proc)>0) THEN
        CALL MPI_RECV( gbuf(:,proc), cnt(proc), MPI_INTEGER, proc, &
                  801, ELMER_COMM_WORLD, status, ierr )

        CALL MPI_RECV( ibuf(:,proc), cnt(proc), MPI_INTEGER, proc, &
                  802, ELMER_COMM_WORLD, status, ierr )
      END IF
    END DO

    DO i=1,Bmat % NumberOfRows
      l = Bmat % Perm(i)
      IF (l>0) THEN ! m = Partition DOF number...
        m = lint(l) ! ...for an interface DOF
      ELSE
        m = -l      ! ...for a Dirichlet DOF
      END IF

      WRITE(5,*) bbeg+i,Bmat% InvPerm(i), Bmat % RHS(i)

      DO j=Bmat % Rows(i),Bmat % Rows(i+1)-1
        proc=Bmat % Cols(j)
        IF(proc==me) THEN
          WRITE(4,*) bbeg+i,m+abeg,bmat % values(j)
        ELSE
          DO k=1,cnt(proc)
            IF(gtags(m)==gbuf(k,proc)) EXIT
          END DO
          IF(k>cnt(proc)) stop 'aah'
          WRITE(4,*) bbeg+i,ibuf(k,proc)+sum(asize(0:proc-1)),Bmat % values(j)
        END IF
      END DO
    END DO

     CLOSE(4)
     CLOSE(5)
!------------------------------------------------------------------------------
  END SUBROUTINE SaveB
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> The FETI main routine: solve x from Ax=b using the F(inite) E(lement)
!> T(earing) and I(nterconnect) method.
!------------------------------------------------------------------------------
  SUBROUTINE  Feti(A,x,b,Solver)
    ! Just basic feti so far...
    ! -------------------------
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp), target :: x(:),b(:)
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE,target :: y(:)
    REAL(KIND=dp), ALLOCATABLE :: lambda(:)
    INTEGER :: SaveNdeg
    REAL(KIND=dp), ALLOCATABLE :: alpha(:)
    REAL(KIND=dp) :: TOL,mind,maxd
    INTEGER :: i,j,k,l,m,n,nd,q,d,nLC

    LOGICAL  :: Found, Floating=.FALSE., QR, MumpsLU
    INTEGER(KIND=AddrInt) :: mvProc, dotProc, nrmProc, stopcProc, precProc
    INTEGER(KIND=AddrInt) :: AddrFunc
    EXTERNAL :: AddrFunc

    REAL(KIND=dp), POINTER CONTIG :: SaveValues(:)
    INTEGER, POINTER CONTIG :: SaveCols(:),SaveRows(:)
    INTEGER, POINTER  :: p(:)

    SAVE SaveValues, SaveCols,  SaveRows

    INTEGER, ALLOCATABLE :: Indexes(:)
    TYPE(Element_t), POINTER :: EL
    TYPE(ValueList_t), POINTER :: BC

#ifdef HAVE_CHOLMOD
    INTERFACE
      SUBROUTINE SPQR_NZ(chol,nz) BIND(c,NAME="spqr_nz")
        USE Types
        INTEGER :: nz
        INTEGER(Kind=AddrInt) :: chol
      END SUBROUTINE SPQR_NZ

      SUBROUTINE SPQR_NullSpace(chol,n,nz,z) BIND(c,NAME="spqr_nullspace")
        USE Types
        INTEGER :: nz,n
        REAL(KIND=dp) :: z(*)
        INTEGER(Kind=AddrInt) :: chol
      END SUBROUTINE SPQR_NullSpace
    END INTERFACE
#endif

!------------------------------------------------------------------------------

    ! Flag dirichlet d.o.f., to remove 'em from interface conditions:
    ! ----------------------------------------------------------------
    IF(ALLOCATED(DirichletDOFs)) DEALLOCATE(DirichletDOFs)

    ALLOCATE(DirichletDOFs(A % NumberOFRows))
    DirichletDOFs = .FALSE.

    IF (ALLOCATED(A % ConstrainedDOF)) THEN
      ! Take the flags from the matrix. They are what was actually imposed on
      ! it, which is not the same as what a scan of this partition's boundary
      ! elements finds: a node can be constrained here while the boundary
      ! element carrying the condition sits in a neighbouring partition, and
      ! conditions set by other means than a boundary condition section are not
      ! visible to such a scan at all.
      !
      ! Rebuilding them here instead is what made the two sides of an interface
      ! disagree about what was on it, since the lists are filtered with
      ! .NOT.DirichletDOFs: one partition kept a DOF and sent its tag, the other
      ! had dropped it, could not match the tag, and reported "should not
      ! happen" while quietly discarding the constraint.
      DirichletDOFs = A % ConstrainedDOF(1:A % NumberOfRows)
    ELSE
      ! No flags on the matrix, so fall back to scanning the boundary elements.
      d  = Solver % Variable % DOFs
      p => Solver % Variable % Perm
      ALLOCATE(Indexes(Solver % Mesh % MaxElementDOFs))
      DO i=1,GetNofBoundaryElements()
        El=>GetBoundaryElement(i)
        BC=>GetBC()
        IF(.NOT.ASSOCIATED(BC)) CYCLE
        nd = GetElementDOFs(Indexes)
        DO j=1,d
          IF (d>1) THEN
            IF(ListCheckPresent(BC,ComponentName(Solver % Variable,j))) THEN
              DO k=1,nd
                IF( p(Indexes(k)) <= 0 ) CYCLE
                DirichletDOFs(d*(p(Indexes(k))-1)+j) = .TRUE.
              END DO
            END IF
          ELSE
            IF(ListCheckPresent(BC,Solver % Variable % Name)) THEN
              DO k=1,nd
                IF( p(Indexes(k)) <= 0 ) CYCLE
                DirichletDOFs(d*(p(Indexes(k))-1)+j) = .TRUE.
              END DO
            END IF
          END IF
        END DO
      END DO
      DEALLOCATE(Indexes)
    END IF

    ! Get various  solution options:
    ! ------------------------------
    Params => ListGetSolverParams()
    TOL=GetCReal( Params,'Linear System Convergence Tolerance')

    ! Check whether to use the 'total' FETI scheme:
    ! ---------------------------------------------

    TotalFeti  = GetLogical(Params, 'Total Feti', Found)
    FetiAsPrec = GetLogical(Params, 'Feti Use As Preconditioner', Found)

    ! Tell 'DirectSolver' to factorize first time only:
    ! -------------------------------------------------
    Refactorize=GetLogical(Params, 'Linear System Refactorize',Found)
    IF (.NOT.Found) Refactorize = .TRUE.
    CALL ListAddLogical(Params,'Linear System Refactorize',.FALSE.)

    ! Check if using local CGP, or (f.ex.) GCR from usual iterators:
    ! --------------------------------------------------------------
    CPG = GetString(Params,'Linear System Iterative Method') == 'cpg'

    ! Check if QR decomposition to find the nullspace:
    ! ------------------------------------------------
    QR = GetString(Params,'Linear System Direct Method') == 'spqr'

    ! Check if local Mumps decomposition can be used to find the nullspace
    ! ---------------------------------------------------------------------
    MumpsLU = GetString(Params,'Linear System Direct Method') == 'mumpslocal'
    ! Mark system as possibly singular for Mumps
    ! ------------------------------------------
    IF (MumpsLU) &
      MumpsLU = ListGetLogical(Params, 'Mumps Solve Singular', Found)

    n = A % NumberOfRows
    ! Sized for the Lagrange coefficients, which can outnumber the rows of the
    ! local matrix, plus the kernel coefficients that both solvers append after
    ! them. maxnz bounds the latter.
    ALLOCATE(y(MAX(n,FetiMaxLC(A)+n)))
    ALLOCATE(lambda(MAX(n,FetiMaxLC(A)+n)))

    ! The kernel dimension is bounded by the number of rows and by nothing
    ! smaller: a subdomain in several disconnected pieces carries the modes of
    ! each, so twelve metis pieces of the winkel mesh already reach 24 against
    ! the 20 these used to be given.
    IF (ALLOCATED(FixInds)) DEALLOCATE(FixInds)
    ALLOCATE(FixInds(n), alpha(n))

    ! Initialize parallel matrix for the original system for
    ! matrix-vector multiply in FetiStopc. NOT USED at the
    ! moment.
    ! ------------------------------------------------------
    ! y=b
    ! CALL ParallelInitSolve(A,x,y,rtmp)
 
    dumptofiles = GetLogical( Params, 'Feti dump system', Found)
    IF(dumptofiles) THEN
      CALL Info( 'Feti:', 'Dumping Feti Description to files')
      CALL SaveKandF(A)
    END IF

    IF(Refactorize) THEN
      InitializeLC = .TRUE.;  InitializeIf = .TRUE.

      ! Check neighbouring partitions:
      ! ------------------------------
      CALL FetiGetNeighbours()

      ! Check for floating domains, initialize null space z (a.k.a. R):
      ! ---------------------------------------------------------------
      IF(.NOT.(QR.OR.MumpsLU)) THEN
        Floating=FetiFloatingDomain(A,Solver,FixInds,TOL)
        CALL FetiCheckKernel(A,'rigid body modes')
      END IF
    END IF

    ! 'fix' A  to allow solving local system:
    ! ---------------------------------------
    IF (.NOT.(QR.OR.MumpsLU).AND.Floating) THEN
      ! Pinning only removes the kernel if the chosen DOFs are a good enough
      ! set for it. Check before relying on them, and fall back to the Lagrange
      ! coefficients when they are not, rather than solving with a matrix that
      ! is still (nearly) singular.
      IF (.NOT.NullSpaceLC) THEN
        IF (.NOT.FetiFixingUsable()) NullSpaceLC = .TRUE.
      END IF

      LocalSolveChecked = .FALSE.

      saverows   => a % rows
      savecols   => a % cols
      savevalues => a % values

      IF (NullSpaceLC) THEN
        ! Fix z=null(A) with Lagrange coefficients
        ! A_fix = [A z^T; z 0] ([x; \lambda]):
        ! ------------------------------------------

        ! CRS_MatrixVectorMultiply has hand unrolled paths selected by A % ndeg
        ! which assume every row holds whole blocks of ndeg consecutive columns.
        ! The assembly below breaks that twice over: it drops stored zeros, and
        ! it appends the nz constraint columns after them. Left as it was, the
        ! unrolled loop reads the wrong entries of the multiplicand -- the
        ! product came out 60% asymmetric, which stalls the conjugate gradient
        ! iteration that preconditions with it. Mark the matrix as unblocked.
        SaveNdeg = a % ndeg
        a % ndeg = 1

        allocate(a % rows(n+nz+1))
        k=count(a % values/=0)+2*count(z/=0)
        allocate(a % cols(k), a % values(k))
        l=1

        do i=1,n
          a % rows(i)=l
          do j=saverows(i),saverows(i+1)-1
            IF ( savevalues(j)==0 ) CYCLE
            a % cols(l)=savecols(j)
            a % values(l)=savevalues(j)
            l=l+1
          end do
          do j=1,nz
            IF ( z(j,i)==0 ) CYCLE
            a % cols(l)=n+j
            a % values(l)=z(j,i)
            l=l+1
          end do
        end do

        do j=1,nz
          a % rows(n+j)=l
          do i=1,n
            IF ( z(j,i)==0 ) CYCLE
            a % cols(l)=i
            a % values(l)=z(j,i)
            l=l+1
          end do
        end do
        a % rows(n+nz+1)=l
      ELSE
        ! Fix z=null(A) with Dirichlet conditions on degrees
        ! of freedom given by FloatingDomain():
        ! ---------------------------------------------------
        ALLOCATE(A % Rows(n+1))
        k=SIZE(A % Values)
        ALLOCATE(A % Cols(k), A % Values(k))

        A % Rows = SaveRows
        A % Cols = SaveCols
        A % Values = SaveValues

        DO j=1,nz
          CALL CRS_SetSymmDirichlet(A,y,FixInds(j),0._dp)
        END DO
        FixingByDofs = .TRUE.
      END IF
    END IF


    ! Compute and distribute r.h.s. for L.C.'s:
    ! -SUM_part(B_i A_i^-1 b_i)
    ! -----------------------------------------
    CALL FetiDirectSolver(A,x,b,Solver)
    x = -x
    nLC = FetiSendRecvIf(A,y,x)

    ! If using QR decomposition, find the null space
    ! after factorizing:
    ! ----------------------------------------------
    IF (Refactorize.AND.QR) THEN
#ifdef HAVE_CHOLMOD
      CALL SPQR_NZ(A % Cholmod,nz) ! get null space size
      NullSpaceLC=.FALSE.
      Floating = nz/=0
      IF(ALLOCATED(z)) DEALLOCATE(z)
      ALLOCATE(z(nz,n))
      CALL SPQR_NullSpace(A % Cholmod,n,nz,z) ! get the null space
      CALL FetiCheckKernel(A,'SPQR')
#else
      CALL Fatal('Feti','Cholmod/SPQR solver has not been installed.')
#endif
    END IF

    ! Find the null space by using LU decomposition and Mumps
    IF (Refactorize .AND. MumpsLU) THEN
#ifdef HAVE_MUMPS
      IF(ALLOCATED(z)) DEALLOCATE(z)

      ! Solve local nullspace with Mumps
      CALL MumpsLocal_SolveNullSpace( Solver, A, z, nz )
      NullSpaceLC=.FALSE. ! What does this do?
      Floating = nz/=0 ! Is the domain a floating one
      CALL FetiCheckKernel(A,'Mumps null pivots')
#else
      CALL Fatal('Feti','Mumps has not been installed.')
#endif
    END IF

    mind=ParallelReduction(nz,1)
    maxd=ParallelReduction(nz,2)
    WRITE(Message,*) 'min/max nz:',FLOOR(mind+0.5_dp),FLOOR(maxd+0.5_dp)
    CALL Info('Feti:', Message,Level=6)

    IF(dumptofiles) THEN
      CALL SaveR()
      CALL Info( 'Feti:', 'File dumping completed, exiting.')
      CALL ParallelFinalize(); STOP EXIT_OK
    END IF
    

    ! add Dirichlet BC contribution to the r.h.s., if using Total FETI:
    ! ------------------------------------------------------------------
    IF (TotalFeti) y(1:nLC)=y(1:nLC)+Bmat % RHS(1:nLC)

    ! Add Lagrange coefficient(s) for the projection condition
    ! if not using the C(onjugate)P(rojected)G(radient) stuff:
    ! --------------------------------------------------------
    IF (Floating .AND. .NOT. CPG) THEN
      y(nLC+1:nLC+nz)=-MATMUL(z,b)
      nLC = nLC+nz
    END IF

    ! CG iteration for the L.C.'s
    ! ----------------------------
    Precondition = GetLogical( Params,'FETI Preconditioning', Found)
    IF (.NOT.Found) Precondition=.TRUE.
    ! The unknown of this solve is the vector of Lagrange coefficients, which
    ! lives in a space of nLC (plus nz) dimensions and has nothing to do with
    ! the rows of the local matrix. It used to be kept in the caller's solution
    ! vector x, which is sized by those rows, and that only held together while
    ! the coefficients happened to be fewer.
    lambda=0
    IF ( CPG ) THEN
      CALL FetiCPG(A,nLC,lambda,y,b,Solver, &
              FetiMV, FetiPrec, SParDotProd, SParNorm)
    ELSE
      precProc=AddrFunc(FetiPrec)
      mvProc=AddrFunc(FetiMV)
      nrmProc=AddrFunc(SParNorm)
      dotProc=AddrFunc(SParDotProd)
      CALL IterSolver(A,lambda,y,Solver,ndim=nLC,DotF=dotProc, &
          NormF=nrmProc, MatvecF=mvProc, precF=precProc )
    END IF

    ! Solve primary unknowns using L.C.'s:
    ! x_i =  A_i^-1(b_i + B_i^T\lambda) + R\alpha:
    ! --------------------------------------------
    IF ( Floating ) THEN
      IF ( CPG ) THEN
        alpha(1:nz)=lambda(nLC+1:nLC+nz)
      ELSE
        alpha(1:nz)=lambda(nLC-nz+1:nLC)
      END IF
    END IF

    ! From here on x is the primary unknown again, and y holds B^T*lambda in the
    ! space of the local rows. The extent is written out because y is sized for
    ! the Lagrange coefficients and is longer than that.
    CALL FetiSendRecvLC(A,y,lambda)
    y(1:n) = y(1:n) + b(1:n)
    CALL FetiDirectSolver(A,x,y,Solver)

    IF (Floating) THEN
      x = x + MATMUL(alpha(1:nz),z)
!     DEALLOCATE(z)

      IF(.NOT.(QR.OR.MumpsLU)) THEN
        DEALLOCATE(A % Values,A % Rows, A % Cols)
        A % ndeg = SaveNdeg
        A % Rows => SaveRows
        A % Cols => SaveCols
        A % Values => SaveValues
      END IF
    END IF

    IF ( Refactorize ) THEN
      Refactorize=GetLogical(Params,'Linear System Free Factorization',Found)
      IF(.NOT.Found) Refactorize=.TRUE.
      IF(Refactorize) CALL DirectSolver(A,x,y,Solver,Free_Fact=.TRUE.)
      CALL ListAddLogical(Params,'Linear System Refactorize', .TRUE. )
    END IF

    IF (TotalFeti .AND. FetiAsPrec) THEN
      ! If using (total)Feti as a preconditioner and using sloppy
      ! tolerances, forcing exact D-condtions afterwards seems to
      ! help outer iteration convergence:
      ! -------------------------------------------------------------
      DO i=1,n
        IF (DirichletDofs(i)) x(i)=b(i)
      END DO
    END IF

    CALL ParallelActiveBarrier()
!------------------------------------------------------------------------------
  END SUBROUTINE  Feti
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  !> Matrix-vector prod. for the L.C. equations:
  !> v_i = SUM_part(B_i A_i^-1 B_i^T u_i).
  !> Called from the CG iterator.
  ! -------------------------------------------------------------------------
  SUBROUTINE  FetiMV(u,v,ipar)
!------------------------------------------------------------------------------
    INTEGER, DIMENSION(*) :: ipar
    REAL(KIND=dp) :: u(HUTI_NDIM), v(HUTI_NDIM), w(HUTI_NDIM)

    REAL(KIND=dp), ALLOCATABLE :: x(:),b(:)
    INTEGER :: n,nLC
    TYPE(Matrix_t), POINTER :: A
    TYPE(Solver_t), POINTER :: Solver

!call resettimer('mv')
    Solver => GetSolver()
    A => GetMatrix()
    n = A % NumberOfRows

    ALLOCATE(x(n),b(n))

    CALL FetiSendRecvLC(A,b,u)
    CALL FetiDirectSolver(A,x,b,Solver)
    nLC = FetiSendRecvIf(A,v,x)

    ! If floating domain, update null space contributions
    ! (if not using CPG):
    ! ----------------------------------------------------
    IF ( .NOT. CPG ) THEN
      x=0
      IF (nz>0) x=MATMUL(u(nLC+1:nLC+nz),z)
      nLC = FetiSendRecvIf(A,w,x)
      v(1:nLC) = v(1:nLC) + w(1:nLC)
      IF (nz>0) v(nLC+1:nLC+nz) = MATMUL(z,b)
    END IF
!call checktimer('mv',delete=.true.)
!------------------------------------------------------------------------------
  END SUBROUTINE  FetiMV
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  !> Peconditioning for feti:
  !> u_i = SUM_part(B_i A_i B_i^T v_i)
  !> (the m-v could be written for the interface dofs only...)
  !> Called from the CG iterator.
  ! -------------------------------------------------------------------------
  SUBROUTINE  FetiPrec(u,v,ipar)
!------------------------------------------------------------------------------
    INTEGER, DIMENSION(*) :: ipar
    REAL(KIND=dp) :: u(HUTI_NDIM), v(HUTI_NDIM)

    REAL(KIND=dp), ALLOCATABLE, TARGET :: x(:),b(:)
    INTEGER :: n, nLC
    TYPE(Matrix_t), POINTER :: A
    TYPE(Solver_t), POINTER :: Solver
!call resettimer('prec')

    IF(.NOT.Precondition) THEN
      u=v
      RETURN
    END IF

    A => GetMatrix()
    n = A % NumberOfRows

    ALLOCATE(x(n+nz),b(n))

    CALL FetiSendRecvLC(A,x,v)
    CALL MatrixVectorMultiply(A,x,b)
    nLC = FetiSendRecvIf(A,u,b)

    ! If floating domain, update null space contributions
    ! (if not using CPG):
    ! ----------------------------------------------------
    IF (.NOT. CPG .AND. nz>0) THEN
      u(nLC+1:nLC+nz)=v(nLC+1:nLC+nz)
    END IF
!call checktimer('prec',delete=.true.)
!------------------------------------------------------------------------------
  END SUBROUTINE  FetiPrec
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END MODULE FetiSolve
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Just a handle for SolveLinearSystem():
!------------------------------------------------------------------------------
SUBROUTINE FetiSolver(A,x,b,Solver)
!------------------------------------------------------------------------------
    USE FetiSolve

    TYPE(Matrix_t), POINTER :: A
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: x(:),b(:)
!------------------------------------------------------------------------------
    CALL Feti(A,x,b,Solver)
!------------------------------------------------------------------------------
END SUBROUTINE FetiSolver
!------------------------------------------------------------------------------

!> \}
