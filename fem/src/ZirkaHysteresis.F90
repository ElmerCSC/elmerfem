MODULE zirka ! Pointwise zirka {{{
USE ISO_C_BINDING, ONLY: C_INT, C_LOC, C_PTR, C_F_POINTER
USE GeneralUtils
USE DefUtils
implicit none

private

INTEGER, PARAMETER :: STACKSIZEINCREMENT = 20
INTEGER, PARAMETER :: POS_SATURATION = 1
INTEGER, PARAMETER :: NEG_SATURATION = -1
INTEGER, PARAMETER :: NO_SATURATION = 0

TYPE, public :: SplineLoop_t  ! {{{
  REAL(KIND=dp), allocatable :: BHasc(:,:), BHdesc(:,:), BHsat(:,:)
  REAL(KIND=dp), pointer :: r_asc(:) => NULL(), r_desc(:) => NULL() , r_sat(:) => NULL() ! Needs to be pointer required by InterpolateCurve
  real(kind=dp), private :: BLimits(2)
  LOGICAL, private :: newmethod = .false.
  CONTAINS 
  PROCEDURE, private :: eval_1 => EvalSplineLoop
  procedure, private :: eval_2 => EvalSplineLoopSingle
  generic, public :: eval => eval_1, eval_2
END TYPE SplineLoop_t ! }}}

TYPE :: RevCurve_t  ! {{{
  REAL(kind=dp) :: Bp, Bq
  REAL(kind=dp) :: a, b, c, dBrev, dHrev
  INTEGER :: depth ! if 0, this is ascending or descending master curve depending on Bp and Bq
  CLASS(RevCurve_t), POINTER :: parent => NULL()
  TYPE(SplineLoop_t), POINTER :: bigloop => NULL()
  PROCEDURE(SimpleEvalRevCurve), POINTER :: simple_eval => SimpleEvalRevCurve
  CONTAINS
  ! PROCEDURE :: simple_eval => SimpleEvalRevCurve
  PROCEDURE :: eval => RecurEvalCurve
END TYPE RevCurve_t ! }}}

type, public :: ZirkaABC_t ! {{{
  REAL(KIND=dp), dimension(1:4) :: Coeffs = [7.73,2.76,-28.63,28.36]
  REAL(KIND=dp) :: b_mult = 0.22
  REAL(KIND=dp) :: c_mult = 0.125
  contains
  procedure, public :: GetABC
end type ! }}}

TYPE, public  :: MasterCurve_t  ! {{{
  TYPE(RevCurve_t), POINTER :: children(:) ! We should always have revcurve_t % depth + 2 to be right index here..
  CLASS(RevCurve_t), POINTER, public :: Head => NULL()
  TYPE(RevCurve_t), POINTER, private :: rc_asc => NULL(), rc_desc => NULL() !  Should this a pointer?
  REAL(KIND=dp), PUBLIC :: B0(3)
  REAL(KIND=dp), ALLOCATABLE, PRIVATE :: BH(:,:)
  REAL(KIND=dp), POINTER, PRIVATE :: r_bh(:) => null()
  TYPE(ZirkaABC_t), POINTER, PRIVATE :: ABCparams => NULL()
  integer, private :: CacheSubSample  = 2
  integer, private :: saturation
  contains
  procedure, public :: eval => mc_eval
  procedure, public :: printme => mc_printme
  procedure, public :: addstack => mc_addstack
  procedure, public :: drive => HBDrive
  procedure, public :: insaturation 
END TYPE MasterCurve_t ! }}}

TYPE, public :: GlobalHysteresisModel_t ! {{{
  TYPE(MasterCurve_t), ALLOCATABLE :: Curves(:,:)
  TYPE(SplineLoop_t), POINTER :: Masterloop => NULL()
  TYPE(ZirkaABC_t), POINTER :: ABCparams => NULL()
  integer :: CacheSubSample
  logical :: initialized = .false.
END TYPE GlobalHysteresisModel_t ! }}}


interface MasterCurve_t
  module procedure init_master_curve
end interface

interface printeval
  module procedure rc_printeval, mc_printeval
end interface

public :: InitSplineLoop, printeval

CONTAINS

!-------------------------------------------------------------------------------
!> Returns true if rc is ascending
PURE FUNCTION ascending(rc) result(Asc) ! {{{
  CLASS(RevCurve_t), intent(in) :: rc
  LOGICAL :: asc
!-------------------------------------------------------------------------------
  if (rc % Bp < rc % Bq) then
    asc = .true.
    return
  end if
  asc = .false.
END FUNCTION ! }}}
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Returns true if rc is descending
PURE FUNCTION descending(rc) result(desc) ! {{{
  CLASS(RevCurve_t), intent(in) :: rc
  LOGICAL :: desc
!-------------------------------------------------------------------------------
  desc = .not. ascending(rc)
END FUNCTION ! }}}
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Returns saturatino value of MasterCurve_t mc
FUNCTION InSaturation(mc, B) RESULT(saturation) ! {{{
  class(MasterCurve_t) :: mc
  real(kind=dp) :: B
  integer :: saturation
  if (B < mc % head % bigloop % blimits(1)) then
    saturation = NEG_SATURATION
    return
  end if
  if (B > mc % head % bigloop % blimits(2)) then
    saturation = POS_SATURATION
    return
  end if
  saturation = NO_SATURATION
END FUNCTION InSaturation ! }}}
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Returns a-b-c parameters appearing in the Zirka model
SUBROUTINE GetABC(this, dBout, dBrev, a, b, c) ! {{{
  IMPLICIT NONE
  !-------------------------------------------------------------------------------
  class(ZirkaABC_t) :: this
  REAL(KIND=dp), INTENT(IN) :: dBout, dBrev
  REAL(KIND=dp), INTENT(OUT) :: a, b, c
  !-------------------------------------------------------------------------------
  REAL(KIND=dp) :: beta
  INTEGER :: k
  !-------------------------------------------------------------------------------

  beta = dBrev/dBout

  a = 0.0_dp
  DO k = 1,4
    a = a + this % Coeffs(k)*(beta**(k-1))
  END DO
  a = a * dBrev
  b = this % b_mult*(1-beta)
  c = this % c_mult

END SUBROUTINE ! }}}
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Initialize spline loop from ascending hysteretic curve and single 
!> valued curve for saturation region
!> BHasc and BHsingle are assumed to be of shape (N,2) 
!-------------------------------------------------------------------------------
FUNCTION InitSplineLoop(BHasc, BHsingle) RESULT(Loop) ! {{{
  !-------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(KIND=dp), intent(in) :: BHasc(:,:), BHsingle(:,:) 
  TYPE(SplineLoop_t), POINTER :: Loop
  !-------------------------------------------------------------------------------
  INTEGER :: N, M
  N = size(BHasc,1)
  M = size(BHsingle,1)

  allocate(loop)
  allocate(loop % BHasc(n, 2), loop % BHdesc(n,2), loop%BHsat(m,2), &
      loop%r_asc(n), loop % r_desc(n), loop % r_sat(m))

  loop % BHasc(:,1:2) = BHasc(:,1:2)
  loop % BHdesc(:,1) = -BHasc(N:1:-1,1)
  loop % BHdesc(:,2) = -BHasc(N:1:-1,2)
  loop % BHsat = BHsingle

  CALL CubicSpline(N, loop % BHasc(1:N, 1), loop % BHasc(1:N, 2), loop % r_asc)
  CALL CubicSpline(N, loop % BHdesc(1:N, 1), loop % BHdesc(1:N, 2), loop % r_desc)
  CALL CubicSpline(M, loop % BHsat(1:M, 1), loop % BHsat(1:M, 2), loop % r_sat)
  loop % newmethod = .true.
  loop % blimits(1) = bhasc(1,1)
  loop % blimits(2) = - loop % blimits(1)

END FUNCTION InitSplineLoop ! }}}

function init_master_curve(bigloop, ABCParams, &
      init, b0, initseq, n_cachesubsample)  result(mc)!  {{{
  IMPLICIT NONE
  TYPE(MasterCurve_t) :: mc
  TYPE(ZirkaABC_t), POINTER :: ABCParams
  TYPE(RevCurve_t), POINTER :: rca, rcd
  TYPE(SplineLoop_t), POINTER :: bigloop
  REAL(KIND=dp), optional, intent(in) :: B0(3)
  INTEGER :: k, N
  integer, optional :: n_cachesubsample
  logical, optional :: init
  real(kind=dp), optional :: initseq(:)

  if (present(n_cachesubsample)) mc % CacheSubSample = n_cachesubsample

  mc % ABCParams => ABCParams
  ! allocate(mc)
  allocate(mc % children(1:STACKSIZEINCREMENT))
  
  DO k = 2,UBOUND(mc % children,1)
    mc % children(k) % parent => mc % children(k-1)
    mc % children(k) % bigloop => bigloop
  END DO
  mc % children(1) % bigloop => bigloop
  allocate(mc % rc_asc, mc % rc_desc)
  mc % children(1) % parent => mc % rc_desc

  rca => mc % rc_asc 
  rcd => mc % rc_desc 
  rca % bigloop => bigloop
  rcd % bigloop => bigloop

  N = size(rca % bigloop % BHAsc,1)

  rcd % Bp = rcd % bigloop % BHdesc(N,1)
  rcd % Bq = rcd % bigloop % BHasc(1,1)
  rcd % depth = 0
  rcd % parent => rca
  rcd % simple_eval => SimpleEvalDescendingSpline

  rca % Bp = rca % bigloop % BHasc(1,1)
  rca % Bq = rca % bigloop % BHdesc(N,1)
  rca % depth = 0
  rca % parent => rcd
  rca % simple_eval => SimpleEvalAscendingSpline
  
  mc % head => rcd

  if (present(init) .or. present(initseq)) then
    if (init .or. present(initseq)) then
      if (.not. present(initseq)) then ! {{{
        block 
          integer, parameter :: nr = 12
          real :: up
          up = rcd % bigloop % bhasc(N,1)
          do n=0,nr-2 ! nr-1
            call mc % drive( up * ((-0.5_dp)**n), cache=.false.)
          end do
          n = nr-1
          call mc % drive( up * ((-0.5_dp)**n), cache=.true.)
        end block
      else
        do n=1,size(initseq)-1
          call mc % drive(initseq(n), cache = .false.)
        end do
        n = size(initseq)
        call mc % drive(initseq(n), cache = .true.)
      end if ! }}}
    end if
  end if

  if (present(b0)) then
    mc % b0(1:3) = b0(1:3)
  else
    mc % b0 = [0.0_dp, 1.0_dp, 0.0_dp]
  end if

END function !  }}}



! Returns curve where B can be evaluated. On saturation, returns the descending (positive saturation) or ascending (negative
! saturation).
!
RECURSIVE FUNCTION RecurseDepth(rc, B) result (rc_p)! {{{
  IMPLICIT NONE
  CLASS(RevCurve_t), POINTER, INTENT(IN) :: rc
  CLASS(RevCurve_t), POINTER :: rc_p
  REAL(kind=dp) :: B
  integer :: tmp

#if DEBUG>1
  integer :: branch
  branch = -1
#endif

  tmp = rc % depth
  ! print *, 'recursedepth: depth, B, Bp, Bq', tmp, B, rc % Bp, rc % Bq
  outer: block
    if (rc % Bp < rc % Bq) then ! Ascending curve (here be dragons)

    if (rc % Bp < B .and. B < rc % Bq) THEN
      rc_p => rc
#if DEBUG>1
      branch = 1
      print *, 'branch: ', branch, rc_p % depth
#endif
      exit outer
    end if
    if (rc % Bq <= B) then 

      if (rc % depth == 0 .and. Ascending(rc)) then ! in positive saturation return descending master
        rc_p => rc % parent
#if DEBUG>1
        branch = 2
        print *, 'branch: ', branch, rc_p % depth
#endif
        exit outer
      end if
#if DEBUG>1
      branch = 3
      print *, 'branch: ', branch
#endif
      rc_p => RecurseDepth(rc % parent % parent, B) ! TODO: Think these through again
      exit outer
    end if
    if (rc % Bp >= B) then 
      if (rc % depth == 0 .and. ascending(rc)) then ! in negative saturation return ascending master
        rc_p => rc
#if DEBUG>1
        branch = 4
        print *, 'branch: ', branch, rc_p % depth
#endif
        exit outer ! Negative saturation
      end if
#if DEBUG>1
      branch = 5
      print *, 'branch: ', branch
#endif
      rc_p => RecurseDepth(rc % parent, B)
      exit outer
    end if
  else ! Descending curve : rc % Bp > rc % Bq
    if (rc % Bp > B .and. B > rc % Bq) then ! inside
      rc_p => rc
#if DEBUG>1
      branch = 6
        print *, 'branch: ', branch, rc_p % depth
#endif
      exit outer
    end if
    if (B <= rc%Bq) then 
      if (rc % depth == 0 .and. descending(rc)) then ! in negative saturation return ascending master
        rc_p => rc % parent
#if DEBUG>1
        branch = 7
        print *, 'branch: ', branch, rc_p % depth
#endif
        exit outer
      end if
#if DEBUG>1
      branch = 8
      print *, 'branch: ', branch
#endif
      rc_p => RecurseDepth(rc % parent % parent, B)
    end if
    if (B >= rc % Bp) then
      if (rc % depth == 0 .and. descending(rc)) then ! in positive saturation return descending master
        rc_p => rc
#if DEBUG>1
        branch = 9
        print *, 'branch: ', branch, rc_p % depth
#endif
        exit outer
      end if
#if DEBUG>1
      branch = 10
      print *, 'branch: ', branch
#endif
      rc_p => RecurseDepth(rc % parent, B)
      exit outer
    end if
  end if 
  end block outer 

END FUNCTION ! }}}

SUBROUTINE EvalSplineLoop(this, B, Hasc, Hdesc, dHasc, dHdesc) ! {{{
  IMPLICIT NONE
  CLASS(SplineLoop_t), INTENT(IN) :: this
  REAL(KIND=dp), INTENT(IN) :: B
  REAL(KIND=dp), INTENT(OUT) :: Hasc, Hdesc
  REAL(KIND=dp), INTENT(OUT), optional :: dHasc, dHdesc

  if (this % newmethod) then
    if (B < this % blimits(1) .or. B > this % blimits(2)) then ! negative saturation
      call this % eval(B, Hasc)
      Hdesc = Hasc
      return 
    end if
  end if

  Hasc = InterpolateCurve( &
      this%BHasc(:,1),  &
      this%BHasc(:,2), &
      B, this%r_asc)

  Hdesc = InterpolateCurve( &
      this%BHdesc(:,1),  &
      this%BHdesc(:,2), &
      B, this%r_desc)
  
  if (present(dHasc)) then
    dHasc = DerivateCurve( &
      this%BHasc(:,1),  &
      this%BHasc(:,2), &
      B, this%r_asc)
  end if

  if (present(dHdesc)) then
  Hdesc = DerivateCurve( &
      this%BHdesc(:,1),  &
      this%BHdesc(:,2), &
      B, this%r_desc)
  end if


END SUBROUTINE EvalSplineLoop ! }}}

SUBROUTINE EvalSplineLoopSingle(this, B, HSingle) ! {{{
  IMPLICIT NONE
  CLASS(SplineLoop_t), INTENT(IN) :: this
  REAL(KIND=dp), VALUE :: B
  REAL(KIND=dp), INTENT(OUT) :: HSingle
  logical :: neg

  neg = B < 0
  if (neg) B = -B

  HSingle = InterpolateCurve( &
      this % BHsat(:,1), &
      this % BHsat(:,2), &
      B, this % r_sat)

  if (neg) Hsingle = -Hsingle
  
END SUBROUTINE EvalSplineLoopSingle ! }}}

! Evaluate rc at B. Don't check if B is consistent with rc%Bp and rc%Bq
RECURSIVE FUNCTION SimpleEvalRevCurve(rc, B) result(H) ! {{{
  CLASS(RevCurve_t), INTENT(IN), target :: rc
  REAL(KIND=dp), INTENT(IN) :: B
  REAL(KIND=dp) :: H
  !-------------------------------------------------------------------------------
  REAL(KIND=dp) :: Hap,  Hp, x, dHout, dH, HMAsc, HMDesc
  !-------------------------------------------------------------------------------
#if 0
  depthtest: if (rc%depth > 2) then 
    Hp = rc % parent % simple_eval(B)
    Hap = rc % parent % parent % simple_eval(B)
  else
    call rc % bigloop % eval(B, HMAsc, HMDesc)
    IF (rc % depth == 2) then
      Hap = HMDesc
      Hp = rc % parent % simple_eval(B)
      exit depthtest
    end if
    IF (rc % depth == 1 ) then
      Hap = HMAsc
      Hp = HMDesc
      exit depthtest
    END IF
    IF (rc % depth == 0) H = HMDesc
    IF (rc % depth == -1) H = HMAsc
  endif depthtest 
#else
  Hp = rc % parent % simple_eval(B)
  Hap = rc % parent % parent % simple_eval(B)
#endif
  if (rc % depth > 0) then 
    x = (rc % Bq - B )/rc % dBrev
    dHout = Hap - Hp
    dH = rc % dHrev*(1.0_dp-rc % b) * x * exp(-rc%a*(1.0_dp-x)) + dHout* rc % b* (x**rc%c)
    H = Hap - dH
  end if
END FUNCTION ! }}}

RECURSIVE FUNCTION SimpleEvalAscendingSpline(rc, B) result(H) ! {{{
  CLASS(Revcurve_t), INTENT(in), target :: rc
  REAL(KIND=dp), intent(in) :: B
  REAL(KIND=dp) :: H, Hdesc
  call rc % bigloop % eval(B, H, Hdesc)
END FUNCTION ! }}}

RECURSIVE FUNCTION SimpleEvalDescendingSpline(rc, B) result(H) ! {{{
  CLASS(Revcurve_t), INTENT(in), target :: rc
  REAL(KIND=dp), intent(in) :: B
  REAL(KIND=dp) :: H, Hasc
  call rc % bigloop % eval(B, Hasc, H)
END FUNCTION ! }}}

! Evaluate master curve at B, 
function mc_eval(mc, B, dhdb, cached) result(H) ! {{{
  class(MasterCurve_t) :: mc
  real(kind=dp), intent(in) :: B
  real(kind=dp) :: H
  real(kind=dp), optional, intent(out) :: dhdb
  real(kind=dp) :: dhasc, dhdesc, hasc, hdesc
  logical, optional :: cached

  if (present(cached) .and. cached ) then 
      H = InterpolateCurve( &
          mc % BH(:,1), &
          mc % BH(:,2), &
          B, mc % r_BH)
    if (present(dhdb)) then
      dhdb = DerivateCurve( &
          mc % BH(:,1), &
          mc % BH(:,2), &
          B, mc % r_BH)
    end if
  else
      H = RecurEvalCurve(mc % head, B)
      if(present(dhdb)) dhdb = 0.0_dp
  end if
end function ! }}}

FUNCTION RecurEvalCurve(rc, B) result (H) ! {{{
  IMPLICIT NONE
  CLASS(RevCurve_t), TARGET,  INTENT(IN) :: rc
  REAL(KIND=dp), INTENT(IN) :: B
  CLASS(RevCurve_t), POINTER :: rc_p
  real(kind=dp) :: H
  integer :: d

  rc_p => rc
  d = rc_p % depth
  rc_p => RecurseDepth(rc_p, B)
  H = rc_p % simple_eval(B)

END FUNCTION ! }}}

SUBROUTINE HBDrive(mc, B, cache) ! {{{
  IMPLICIT NONE
  CLASS(MasterCurve_t), INTENT(INOUT) :: mc
  CLASS(RevCurve_t), POINTER :: rc
  REAL(KIND=dp), INTENT(IN) :: B
  !-------------------------------------------------------------------------------
  CLASS(RevCurve_t), POINTER :: rc_p
  integer :: m, n, n_cache
  logical, optional :: cache
  logical :: cache_
  real(kind=dp) :: bmax, bmin
  !-------------------------------------------------------------------------------

#if DEBUG>1
  print *, 'hbdrive driving stack to ', B
#endif
  rc_p => RecurseDepth(mc % head, B)
  call AddStack(rc_p, mc, B)
  rc => rc_p

  if (.not. present(cache)) then
    cache_ = .true.
  else
    cache_ = cache
  end if

    if(cache_) then
#if DEBUG>1
      print *, rc % bigloop % blimits(1:2)
#endif
      n = size(rc % bigloop % bhsat,1)
      bmax = rc % bigloop % bhsat(n,1)
      bmin = -rc % bigloop % bhsat(n,1)
      n_cache = size(rc % bigloop % bhsat,1) / mc % cachesubsample 

      IF(.not. ALLOCATED(mc % BH)) ALLOCATE(mc % BH(n_cache,2))
      IF(.not. associated(mc % r_bh)) ALLOCATE(mc % r_bh(n_cache))

      DO n = 1, size(mc % BH,1)
        associate( Bhere => (bmin + (n-1)*(bmax-bmin)/n_cache ) )
          mc % BH(n,1) = Bhere
          mc % BH(n,2) = mc % eval(Bhere, cached =.false.)
        end associate
      END DO
      n = size(mc % bh,1)
      CALL CubicSpline(n, mc % BH(:,1), mc % BH(:,2), mc % r_bh, monotone=.true.)
  end if

END SUBROUTINE  !}}}

subroutine mc_addstack(mc, B) ! {{{
  class(MasterCurve_t) :: mc
  real(kind=dp), intent(in) :: B

  call addstack(mc % head, mc, B)
end subroutine ! }}}

subroutine mc_printme(mc) ! {{{
  class(MasterCurve_t) :: mc
  call PrintRevCurve(mc % head)
end subroutine ! }}}

SUBROUTINE AddStack(parent, master, B) ! {{{
  IMPLICIT NONE
  CLASS (RevCurve_t), INTENT(INOUT), POINTER :: parent
  TYPE(MasterCurve_t), INTENT(INOUT) :: master
  CLASS(RevCurve_t), POINTER:: x
  REAL(KIND=dp), INTENT(IN) :: B
  !-------------------------------------------------------------------------------
  real(KIND=dp) :: HMAsc, HMDesc, Hpp, hp, dBout
  integer :: d, d_add

  ! TODO: This looks weird
  ! print *, 'recurse: ', parent % depth
  ! Check for saturation
#if 0
  if (descending(parent) .and. B > parent % Bp) then ! depth == 0 => descending outer
    master % head => parent
    return
  end if
  if (ascending(parent)  .and. B < parent % Bp) then ! depth == -1 => ascending outer
    master % head => parent
    return
  end if
#else
  if (descending(parent)) then ! descending master
    if ( B <= parent % Bq ) then ! in negative saturation return ascending master
      master % head => parent % parent
      return
    end if
    if ( B > parent % Bp ) then ! in positive saturation return descending master
      master % head => parent 
      return
    end if
  end if
  if (ascending(parent)) then ! ascending master
    if ( B >= parent % Bq ) then !  in positive saturation return descending master
      master % head => parent % parent
      return 
    end if
    if ( B < parent % Bp ) then ! in negative saturation return ascending master
      master % head => parent  
      return
    end if
  end if
#endif

  
  ! Not in saturation
  parent => check_reallocate(parent, parent % depth + 1)
  x => master % children(parent % depth + 1)
  x % depth = parent % depth + 1

  ! parent => master % children ( parent % depth + d_add - 1)
  x % parent => parent
  x % Bp = B
  x % Bq = x % parent % Bp
  x % dBrev = x%Bq - x%Bp
  dBout = x % Bq - parent % parent % Bp
  call master % ABCparams % GetABC(abs(dBout), abs(x % dBrev), x%a, x%b, x%c)
  Hpp = x % parent % parent % simple_eval(B)
  Hp = x % parent % simple_eval(B)
  x % dHrev = Hpp - Hp;
  ! parent => x
  master % head  => x


  CONTAINS
  function check_reallocate(oldparent, newdepth) result(newparent) ! {{{
    IMPLICIT NONE
    INTEGER :: newdepth
    !-------------------------------------------------------------------------------
    TYPE(Revcurve_t), POINTER :: newchildren(:)
    ! type(Revcurve_t), pointer :: x
    TYPE(SplineLoop_t), POINTER :: bigloop
    INTEGER :: bufbound, k
    CLASS(RevCurve_t), POINTER :: oldparent, newparent
    !-------------------------------------------------------------------------------
    bufbound = ubound(master % children, 1)
    ! bigloop => master % children(1) % bigloop
    bigloop => oldparent % bigloop
    ! nullify(parent)

    if (newdepth > 1) THEN
      IF (newdepth > bufbound) THEN
        allocate(newchildren(1:bufbound + STACKSIZEINCREMENT))
        newchildren(1:bufbound) = master % children(1:bufbound)
        do k =1,bufbound
          nullify(master % children(k) % bigloop)
        end do
        nullify(master % head)

        deallocate(master % children)
        allocate(master % children(1:bufbound + STACKSIZEINCREMENT), source=newchildren)
        do k =1,newdepth
          master % children(k) % bigloop => bigloop
        end do

        master % head => master % children(newdepth)
      END IF
      master % children(newdepth) % bigloop => bigloop
      master % head => master % children(newdepth)
      do k = 2, newdepth
        master % children(k) % parent => master % children(k-1)
      end do
      newparent => master % children(newdepth-1)
      ! master % children(-1) % parent => master % children(0)
    else
      newparent => oldparent
    end if


  END function !}}}
END SUBROUTINE AddStack ! }}}

subroutine mc_printeval(mc, B, mc2) ! {{{
  class(MasterCurve_t) :: mc
  class(MasterCurve_t), optional:: mc2
  real(kind=dp) :: B
  if (present(mc2)) then
    call printeval(mc % head, B, mc2 % head)
  else
    call printeval(mc % head, B)
  end if

end subroutine ! }}}

SUBROUTINE rc_printeval(rc, B, rc0) ! {{{
  implicit none
  class(revcurve_t), pointer, intent(in) :: rc
  class(revcurve_t), pointer, intent(in), optional :: rc0
  real(kind=dp), intent(in) :: B 
  real(kind=dp) :: X, X0
  class(revcurve_t), pointer :: rc_p, rc0_p
  integer :: k
  rc_p => rc
  if (present(rc0)) rc0_p => rc0
  k = rc_p % depth
  do while (.not. associated(rc_p, rc_p % parent % parent))
    if(present(rc0)) X0 = rc0_p % simple_eval(B)
    X = rc_p % simple_eval(B)
    if (present(rc0)) print *, X, X0, X-X0
    if (.not. present(rc0)) print *, X, rc_p % depth ! , c_loc(rc_p), c_loc(rc_p % parent)
    rc_p => rc_p % parent
    if(present(rc0)) rc0_p => rc0_p % parent
  end do
  X = rc_p % simple_eval(B)
  if(present(rc0)) X0 = rc0_p % simple_eval(B)
  if (present(rc0)) print *, X, X0, X-X0
  if (.not. present(rc0)) print *, X, rc_p % depth ! , c_loc(rc_p), c_loc(rc_p % parent)
END SUBROUTINE rc_printeval ! }}}
! Debugging purposes only
subroutine PrintRevCurve(rc0) ! {{{
  CLASS(RevCurve_t), POINTER, intent(in) :: rc0
  CLASS(RevCurve_t), POINTER :: rc

  rc=>rc0

  write (*,*) 'current depth:', rc0 % depth
  do while (.not. associated(rc%parent % parent, rc) .and. associated(rc))
    write (*,*) 'depth, Bp, Bq, dBrev, dHrev:', rc%depth, rc%Bp, rc%Bq, rc % dBrev, rc % dHrev
    rc => rc % parent
  end do
  write (*,*) 'depth, Bp, Bq, dBrev, dHrev:', rc%depth, rc%Bp, rc%Bq, rc % dBrev, rc % dHrev
  rc => rc % parent
  write (*,*) 'depth, Bp, Bq, dBrev, dHrev:', rc%depth, rc%Bp, rc%Bq, rc % dBrev, rc % dHrev
  rc => rc % parent
end subroutine ! }}}

!-------------------------------------------------------------------------------
SUBROUTINE FreeRevCurve(mc) ! {{{
!-------------------------------------------------------------------------------
  implicit none
  TYPE(MasterCurve_t) :: mc
!-------------------------------------------------------------------------------
  DEALLOCATE(mc % children)
END SUBROUTINE ! }}}
!-------------------------------------------------------------------------------

END MODULE ! }}}

module ZirkaUtils ! Utils for 2D/3D calculations {{{
use zirka
use DefUtils

character(len=*), parameter :: default_zirka_variable_name = 'zirka_ipvar'

contains
!-------------------------------------------------------------------------------
! Initialization for 2D case
!-------------------------------------------------------------------------------
SUBROUTINE InitHysteresis(Model,Solver) ! {{{
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
!-------------------------------------------------------------------------------
  integer :: elemind, t, ipindex, n, nd
  TYPE(GaussIntegrationPoints_t) :: IP
  TYPE(SplineLoop_t), POINTER :: bigloop
  REAL(kind=dp), POINTER :: BHasc(:,:), BHsingle(:,:), initseq(:,:)
  TYPE(Nodes_t) :: Nodes
  TYPE(Element_t), POINTER::Element
  type(ValueList_t), pointer :: material
  type(ValueListEntry_t), pointer :: ZirkaEntry
  INTEGER, PARAMETER :: n_dir_default = 8
  logical :: found, zeroinit, HasZirka
  integer :: n_dir, n_zirka_mat, n_initialized, n_cachesubsample
  CHARACTER(len=MAX_NAME_LEN) :: str
  type(GlobalHysteresisModel_t), POINTER :: zirkamodel
  type(Variable_t), POINTER :: hystvar
  type(ZirkaABC_t), POINTER :: ABCParams
  real(kind=dp), POINTER :: zcoeff(:,:)
  real(kind=dp) :: b_mult, c_mult
!-------------------------------------------------------------------------------
  n_zirka_mat = 0
  n_initialized = 0


  do n = 1,model % NumberOfMaterials ! {{{
    material => model % materials(n) % values
    ZirkaEntry => listfind(material, 'zirka material', haszirka)
    if (.not. HasZirka) CYCLE
    if (.not. ZirkaEntry % LValue) CYCLE
    n_zirka_mat = n_zirka_mat + 1

    if(ListGetLogical(material, 'zirka initialized', HasZirka)) then
      n_initialized = n_initialized + 1
      CYCLE
    end if

    zirkamodel => NULL()
    BHasc => NULL()
    BHSingle => NULL()

    hystvar => GetZirkaVariable(Material)
    
    if (.not. associated(hystvar)) HystVar => CreateZirkaVariable(Material)
    ! call fatal('InitHysteresis', 'hystvar > ' // trim(str) //' < not available')

    n_dir = ListGetInteger(material, 'Zirka Directions', found)
    if (.NOT. found) n_dir = ListGetInteger(material, 'n_dir', found)
    IF(.NOT. Found) n_dir = n_dir_default

    zeroinit = ListGetLogical(material, 'Init to zero', zeroinit)

    ipindex = GetIpCount(usolver=solver, ipvar=hystvar)
    allocate(zirkamodel)
    allocate(zirkamodel % curves(n_dir, ipindex))
    ZirkaEntry % PROCEDURE = TRANSFER(c_loc(zirkamodel), ZirkaEntry % PROCEDURE)

    ! Initalize saturation loops and saturation curve
    call GetConstRealArray(Material, BHasc, 'Ascending BH curve', found) 
    if(.not. Found) call fatal('InitHysteresis', 'Ascending BH curve not found')
    call GetConstRealArray(Material, BHSingle, 'Single valued BH curve', found)
    if(.not. Found) call fatal('InitHysteresis', 'Single valued BH curve not found')

    zirkamodel % masterloop => InitSplineLoop(BHasc, BHSingle) 
    zirkamodel % CacheSubSample = ListGetInteger(Material, 'Zirka Spline Cache Subsample', Found)
    if(.not. Found) Zirkamodel % CacheSubSample = 2

    ! Set up zirka ABC data here
    allocate(ABCParams)
    zirkamodel % ABCParams => ABCParams
    call GetConstRealArray(Material,  zcoeff, 'Zirka model coefficients', found)
    if(Found .and. size(zcoeff,1) == 4) then
      ABCparams % coeffs(1:4) = zcoeff(1:4,1)
    end if
    
    b_mult = GetCReal(Material, 'Zirka model b multiplier', found)
    if(found) zirkamodel % abcparams % b_mult = b_mult

    c_mult = GetCReal(Material, 'Zirka model c multiplier', found)
    if(found) zirkamodel % abcparams % c_mult = c_mult
  end do ! }}}

  if (.not. n_zirka_mat == n_initialized) then ! only init if some
    do elemind = 1,getnofactive() ! {{{
      element => GetActiveElement(elemind)

      material => GetMaterial()
      ZirkaEntry => ListFind(material, 'zirka material', haszirka)
      if (.not. HasZirka) CYCLE

      initseq => null()
      call GetConstRealArray(Material, initseq, 'zirka init sequence', found)
      if(.not. found) initseq => null()

      if(ListGetLogical(material, 'zirka initialized', HasZirka)) CYCLE

      IP = GaussPoints(element)

      inithystblock: block ! {{{
        REAL(KIND=dp) :: phi 
        integer :: i

        zirkamodel => GetZirkaPointer(material)

        if(.not. associated(zirkamodel)) then
          print *, 'zirkamodel is not associated, cycling' ! DEBUG
          cycle
        end if

        do t = 1, ip%N
          ipindex = getipindex(t, usolver=solver, element=element, ipvar=hystvar)
          if (ipindex == 0) CYCLE
          do i = 1,n_dir
            phi = (i-1.0)*pi/n_dir
            ! print *, 'hello', i, t, ipindex ! DEBUG
            if(.not. associated(initseq)) then
              zirkamodel % curves(i, ipindex) = &
                  MasterCurve_t(zirkamodel % masterloop, ABCparams, &
                  init=zeroinit, b0=[sin(phi), cos(phi), 0.0_dp], &
                  n_cachesubsample = zirkamodel % CacheSubSample)
            else
              zirkamodel % curves(i, ipindex) = &
                  MasterCurve_t(zirkamodel % masterloop, ABCparams, &
                  init=.true., b0=[sin(phi), cos(phi), 0.0_dp], initseq=initseq(:,1), &
                  n_cachesubsample = zirkamodel % CacheSubSample)
            end if
          end do
        end do
      end block  inithystblock ! }}}

    end do ! }}}
  end if

  ! set initialized to true
  do n = 1,model % NumberOfMaterials ! {{{
    Material => Model % Materials(n) % Values
    ZirkaEntry => ListFind(material, 'zirka material', haszirka)
    if (.not. HasZirka) CYCLE

    if(ListGetLogical(material, 'zirka initialized', HasZirka)) cycle
    call ListAddLogical(material, 'zirka initialized', .true.)
  end do ! }}}
  write(message, '(A,I0,A,I0)') 'Found ',n_zirka_mat,' Zirka hysteresis models and skipped initialization of ', n_initialized
  call Info('InitHysteresis',message,level=10)
  
!-------------------------------------------------------------------------------
END SUBROUTINE InitHysteresis ! }}} 
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Returns a pointer to the object containing hysteresis model of given material
!-------------------------------------------------------------------------------
FUNCTION GetZirkaPointer(material) result(ZirkaModelPtr) ! {{{
!-------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER, intent(in) :: Material
  TYPE(GlobalHysteresisModel_t), POINTER :: ZirkaModelPtr
!-------------------------------------------------------------------------------
  TYPE(ValueListEntry_t), POINTER :: ZirkaEntry
  TYPE(C_PTR) :: c_zirka_ptr
  logical :: haszirka

  ZirkaModelPtr => Null()
  ZirkaEntry => ListFind(material, 'zirka material', haszirka)
  IF(.NOT. HasZirka) RETURN

  c_zirka_ptr = transfer(ZirkaEntry % PROCEDURE, c_zirka_ptr)
  call C_F_POINTER(c_zirka_ptr, ZirkaModelPtr)
END FUNCTION ! }}}
!-------------------------------------------------------------------------------

FUNCTION CreateZirkaVariable(Material) RESULT(var)
!-------------------------------------------------------------------------------
  USE MainUtils
  IMPLICIT NONE
  TYPE(ValueList_t), POINTER, intent(in) :: Material
  TYPE(Variable_t), POINTER :: var
!-------------------------------------------------------------------------------
  CHARACTER(len=MAX_NAME_LEN) :: str
  type(mesh_t), pointer :: mesh
  logical :: found

  mesh => getmesh()

  str = ListGetString(material, 'Zirka variable', found) 
  if(.not. found) str = default_zirka_variable_name !'zirka'
  var => VariableGet(mesh % variables, str)
  if(associated(var)) THEN
    write (Message, '(A)') 'Attempting to create existing variable > ' // trim(str) // ' <'
    call warn('CreateZirkaVariable', Message)
    return
  END IF
  BLOCK
    LOGICAL :: HasZirka, Found
    TYPE(ValueListEntry_t), POINTER :: zmat
    CHARACTER(len=MAX_NAME_LEN) :: maskname, mask_sec
    INTEGER, POINTER :: Perm(:)
    REAL(KIND=dp), POINTER :: Values(:)
    type(Solver_t), POINTER :: PSolver
    type(Model_t) :: Model
    integer :: nsize


    PSolver => GetSolver()
    maskname = 'zirka material'
    mask_sec = 'material'
    zmat => ListFind( material , 'Zirka Material' , Found)
    IF (Found .and. zmat % LValue) THEN
      call CreateIpPerm(PSolver, Perm, maskname, mask_sec)
      nsize = maxval(perm)
      allocate(Values(nsize))
      call VariableAdd( PSolver % Mesh % Variables, PSolver % Mesh, PSolver, &
          trim(str), 1, Values, Perm, output = .TRUE., TYPE=Variable_on_gauss_points)
    END IF
  END BLOCK
  var => VariableGet(mesh % variables, str)
END FUNCTION CreateZirkaVariable

!-------------------------------------------------------------------------------
! Returns a pointer to the variable holding hysteresis models
!-------------------------------------------------------------------------------
FUNCTION GetZirkaVariable(Material) RESULT(ZirkaVariable) ! {{{
!-------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER, intent(in) :: Material
  TYPE(Variable_t), POINTER :: ZirkaVariable
!-------------------------------------------------------------------------------
  CHARACTER(len=MAX_NAME_LEN) :: str
  type(mesh_t), pointer :: mesh
  logical :: found

  mesh => getmesh()

  str = ListGetString(material, 'Zirka variable', found) 
  if(.not. found) str = default_zirka_variable_name! 'zirka'
  ZirkaVariable => VariableGet(mesh % variables, str)
  if(associated(zirkaVariable)) return
  
!-------------------------------------------------------------------------------
END FUNCTION ! }}}
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Sets the hysteresis state globally in 2D case
!-------------------------------------------------------------------------------
SUBROUTINE DriveHysteresis(model, solver) ! {{{
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  type(model_t) :: model
  type(solver_t) :: solver
!-------------------------------------------------------------------------------
  integer :: elemind, ipindex, n, nd
  TYPE(GaussIntegrationPoints_t) :: IP
  ! type(RevCurve_t), POINTER :: rc
  TYPE(SplineLoop_t), POINTER :: bigloop
  REAL(kind=dp) :: H(1:5), Basc(1:5), Bdesc(1:5)
  TYPE(Nodes_t) :: Nodes
  TYPE(Element_t), POINTER::Element
  TYPE(Variable_t), POINTER :: hystvar
  TYPE(GlobalHysteresisModel_t), pointer :: zirkamodel
  logical :: found
!-------------------------------------------------------------------------------
#if DEBUG>0
  integer :: printed 
  printed =  0
#endif

  do elemind = 1,getnofactive()
    element => GetActiveElement(elemind)

    n  = GetElementNOFNodes(Element)
    nd = GetElementNOFDOFs(Element)
    IP = GaussPoints(element)
    if (.not. ListGetLogical(GetMaterial(Element), 'zirka material',found)) CYCLE
    zirkamodel => GetZirkaPointer(GetMaterial(Element))
    hystvar => GetZirkaVariable(GetMaterial(Element))


drivehystblock: block
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ, B_ip(3), POT(nd), Agrad(3)
    LOGICAL :: einfostat
    REAL(kind=dp), parameter :: zstab=1.0e-5
    integer :: counter, t, n_dir

    CALL GetLocalSolution(POT,UElement=Element,USolver=Solver)
    CALL GetElementNodes(Nodes)

    counter = 0
    do t = 1, ip%N
      ipindex = getipindex(t, usolver=solver, element=element, ipvar=hystvar)
      if (ipindex == 0) cycle
#if DEBUG>0
      if (printed == 0) then ! DEBUG
        call zirkamodel % curves(1,ipindex) % printme()
        printed = printed + 1
      end if
#endif

      counter = counter + 1
      einfostat = ElementInfo( & 
          Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )
      ! B_ip(1) = -SUM(POT*dBasisdx(:,2))
      ! B_ip(2) = SUM(POT*dBasisdx(:,1))
      ! B_ip(3) = 0.0_dp

      Agrad = matmul(POT, dbasisdx)
      B_ip(1) = Agrad(2)
      B_ip(2) = -Agrad(1)
      b_ip(3) = 0.0_dp
#if DEBUG>0
      if (printed < 2) then ! DEBUG
        print *, 'B_ip = ', B_ip, printed
        printed = printed + 1 
      end if
#endif

      do n_dir = 1, ubound(zirkamodel % curves, 1)
        call zirkamodel % curves(n_dir, ipindex) % &
            drive(sum(zirkamodel % curves(n_dir, ipindex) % B0 * B_ip))
      end do
#if DEBUG>0
      if (printed < 3) then ! DEBUG
        call zirkamodel % curves(1,ipindex) % printme()
        printed = printed + 1
      end if
#endif
    end do

end block drivehystblock

  end do
!-------------------------------------------------------------------------------
END SUBROUTINE DriveHysteresis ! }}}
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Add magnetic field strength to long vector force used in MagnetoDynamicsCalcFields
!-------------------------------------------------------------------------------
SUBROUTINE GetHystereticMFS(Element, Force, pSolver, HasZirka, CSymmetry) ! {{{
!-------------------------------------------------------------------------------
  use DefUtils
  type(element_t), intent(in) :: element
  real(kind=dp), intent(inout) :: force(:,:)
  type(solver_t), intent(in) :: pSolver
  logical, intent(out) :: HasZirka
  LOGICAL, OPTIONAL :: CSymmetry
!-------------------------------------------------------------------------------
  type(ValueList_t), pointer :: material

  LOGICAL :: einfostat
  integer :: n, nd, t, n_dir
  type(GaussIntegrationPoints_t) :: IP
  type(GlobalHysteresisModel_t), pointer :: zirkamodel
  type(variable_t), pointer :: hystvar
  type(nodes_t) :: nodes
  LOGICAL :: CSymmetry_

  IF(.NOT. present(CSymmetry)) THEN
    CSymmetry_ = .FALSE.
  ELSE
    CSymmetry_ = CSymmetry
  END IF

  Material => GetMaterial(Element)
  HasZirka = ListGetLogical(Material, 'zirka material', HasZirka)
  if (.not. HasZirka) RETURN
  zirkamodel => GetZirkaPointer(GetMaterial(Element))
  hystvar => GetZirkaVariable(GetMaterial(Element))

  n  = GetElementNOFNodes(Element)
  nd = GetElementNOFDOFs(Element)
  IP = GaussPoints(element)
  block 
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ, B_ip(2), POT(nd),H_ip(2), &
        Agrad(3), x, Alocal
    integer :: p, ipindex

    call getlocalsolution(pot, uelement=element, usolver=psolver)
    CALL GetElementNodes(Nodes)

    DO t=1,IP % n
      ipindex = getipindex(t, usolver=psolver, element=element, ipvar=hystvar)
      if (ipindex == 0) cycle
      H_ip = 0.0_dp
      einfostat = ElementInfo( & 
          Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )

      ! No symmetries yet
      Agrad = matmul(POT, dbasisdx)
      B_ip(1) = Agrad(2)
      B_ip(2) = -Agrad(1)
      IF( CSymmetry ) THEN
        Alocal = sum(Pot(1:n)*basis(1:n))
        x = SUM( Basis(1:n) * Nodes % x(1:n) )
        detJ = detJ * x
        B_ip = -B_ip
        B_ip(2) = B_ip(2) + Alocal/x
      END IF

      do n_dir = 1, ubound(zirkamodel % curves, 1)
        associate(B0 => zirkamodel % curves(n_dir, ipindex) % B0)
          H_ip = H_ip + zirkamodel % curves(n_dir,ipindex) % &
              eval(sum(B_ip*B0(1:2)), cached = .true.) * &
              B0(1:2)
        end associate
      end do
      do p=1,n
        FORCE(p,1:2) = FORCE(p,1:2) + ip%s(t)*detJ*H_ip(1:2)*Basis(p)
      end do
    end do

  end block

!-------------------------------------------------------------------------------
END SUBROUTINE  ! }}}
!-------------------------------------------------------------------------------
end module ! }}}

