!-------------------------------------------------------------------------------
SUBROUTINE ZirkaTest( Model,Solver,dt,TransientSimulation ) ! {{{
!-------------------------------------------------------------------------------
  USE DefUtils
  USE ZirkaUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  integer :: unit, testnum
  TYPE(ValueList_t), POINTER :: params
  type(Variable_t), POINTER :: variable
  logical :: found
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  
  testnum = ListGetInteger(GetSolverparams(Solver), 'test number', found)
  if (.not. found) testnum = 6

  call info('ZirkaTest','Running zirkatest')
  select case(testnum)
  case (6); call test6()
  case (5); call test5()
  case (4); call test4()
  case (3); call test3()
  end select

  OPEN (newunit=unit, FILE='TEST.PASSED')
  write(unit, '(I0)') 1
   
  close(unit)

contains

subroutine test6() ! {{{
  use zirka
  use defutils
  type(MasterCurve_t) :: mc, mc2
  real(kind=dp), allocatable :: BHasc(:,:), BHsingle(:,:), Bsingle(:), Hsingle(:)
  real(kind=dp), allocatable :: rps(:)
  real(kind=dp), allocatable :: H(:), dH(:), Hnobuf(:)
  real(kind=dp), pointer :: rps_array(:,:)
  integer  :: n, m
  integer, parameter :: nr = 20
  type(ZirkaABC_t), POINTER :: ABCParams
  allocate(ABCParams )

  BHasc = readBH('../HB/BH_asc_real')
  BHsingle = readBH('../HB/BH_single_real')

  allocate(H(size(BHasc,1)), dH(size(BHasc,1)), Hnobuf(size(BHasc,1)))

  allocate(rps(nr))
  rps = [(2.0*(-0.5_dp)**n, n=0,nr-1)]

  mc = MasterCurve_t(InitSplineLoop(BHasc, BHsingle), &
      ABCparams, n_cachesubsample = 2)

  call mc % drive(1.8_dp)

  DO n = 1,nr-1
    CALL mc % drive(rps(n), cache=.false.)
  END DO
  CALL mc % drive(rps(nr), cache=.true.)

  rps_array => null()
  m = 250
  BLOCK
    integer :: i_r
    rps_array => ListGetConstRealArray(GetSolverParams(), 'rps', found)
    if (.not. found) then
      rps = [0.0_dp, 1.8_dp, -1.0_dp, 0.5_dp, 0.25_dp, 0.4_dp, -0.4_dp, -0.8_dp, 0.9_dp]
    else
      rps = rps_array(:,1)
    end if

    ! rps = -rps
    ! rps = [0.0_dp, -1.84102067_dp, 1.3_dp, -0.5_dp]
    do i_r = 2, ubound(rps,1)
      call mc % drive(rps(i_r-1), cache = .true.)
        associate (B_low => rps(i_r-1),  B_up => rps(i_r))
          call print_interval(mc, B_low, B_up, m)
        end associate
    end do
  END BLOCK

end subroutine ! }}}

subroutine print_interval(mc, B_low, B_up, m, Hexact) ! {{{
  type(MasterCurve_t) :: mc
  real(kind=dp) :: B_low, B_up, H, Hnobuf, dH
  real(kind=dp), intent(in), optional :: Hexact
  integer :: n, m
  !-------------------------------------------------------------------------------
  DO n = 0, m-1
    associate( Bhere => (B_low + n*(B_up-B_low)/m) )
      H = mc % eval(Bhere, cached = .true., dhdb=dH)
      Hnobuf = mc % eval(Bhere, cached = .false.)

      if (present(Hexact)) then
        print *, Bhere, H, dH, Hnobuf, Hexact
      else
        print *, Bhere, H, dH, Hnobuf
      end if

    end associate
  end do
end subroutine ! }}}

subroutine test5() ! {{{
  use zirka
  use defutils
  type(MasterCurve_t) :: mc, mc2
  real(kind=dp), allocatable :: BHasc(:,:), BHsingle(:,:), Bsingle(:), Hsingle(:)
  real(kind=dp), allocatable :: rps(:)
  real(kind=dp), allocatable :: H(:), dH(:), Hnobuf(:)
  integer  :: n, m
  integer, parameter :: nr = 20
  type(ZirkaABC_t), POINTER :: ABCParams
  allocate(ABCParams )

  BHasc = readBH('../HB/BH_asc_real')
  BHsingle = readBH('../HB/BH_single_real')

  allocate(H(size(BHasc,1)), dH(size(BHasc,1)), Hnobuf(size(BHasc,1)))

  allocate(rps(nr))
  rps = [(2.0*(-0.5_dp)**n, n=0,nr-1)]

  mc = MasterCurve_t(InitSplineLoop(BHasc, BHsingle), &
      ABCparams)

  call mc % drive(1.8_dp)

  DO n = 1,nr-1
    CALL mc % drive(rps(n), cache=.false.)
  END DO
  CALL mc % drive(rps(nr), cache=.true.)


  m = 250
  associate (B_low => -0.0_dp, B_up => 1.9_dp )
    DO n = 1,m
      associate( Bhere => (B_low + n*(B_up-B_low)/m) )
        H(n) = mc % eval(Bhere, cached = .true., dhdb=dH(n))
        Hnobuf(n) = mc % eval(Bhere)
        print *, Bhere, H(n), dH(n), Hnobuf(n)
      end associate
    END DO
  end associate

end subroutine ! }}}
!-------------------------------------------------------------------------------
subroutine test4() ! {{{
  use zirka
  use defutils
  type(MasterCurve_t) :: mc, mc2
  real(kind=dp), allocatable :: BHasc(:,:), BHsingle(:,:), Bsingle(:), Hsingle(:)
  real(kind=dp), allocatable :: rps(:)
  real(kind=dp) :: H
  integer  :: n, m
  integer, parameter :: nr = 20
  type(ZirkaABC_t), POINTER :: ABCParams
  real(kind=dp), pointer :: rps_array(:,:)
  allocate(ABCParams )

  BHasc = readBH('../HB/BH_asc_real')
  BHsingle = readBH('../HB/BH_single_real')


  allocate(rps(nr))
  rps = [(2.0*(-0.5_dp)**n, n=0,nr-1)]

  mc = MasterCurve_t(InitSplineLoop(BHasc, BHsingle), &
      ABCparams)

  rps_array => null()
  rps_array => ListGetConstRealArray(GetSolverParams(), 'rps_4', found)

  do n=1,size(rps_array,1)
    call mc % drive(rps_array(n,1), cache = .false.)
    call mc % printme()
    print *, 'B = ', rps_array(n,1) 
    H = mc % eval (rps_array(n,1), cached = .false.)
    print *, 'H = ',H
    ! call printeval(mc, rps_array(n,1))
  end do
  ! call mc % drive(0.1_dp, cache = .false.)
  ! call mc % printme()
  ! call mc % drive(-1.8_dp, cache = .false.)
  ! call mc % printme()
  ! call mc % drive(0.8_dp, cache = .false.)
  ! call mc % printme()
  ! call mc % drive(-1.2_dp, cache = .false.)
  ! call mc % printme()

  ! deallocate(rps_array)

end subroutine ! }}}
subroutine test3() ! {{{
  use zirka
  use defutils
  type(MasterCurve_t) :: mc, mc2
  real(kind=dp), pointer :: BHasc(:,:), BHsingle(:,:)
  real(kind=dp), pointer :: rps(:), H(:)
  real(kind=dp), pointer :: rps_array(:,:)
  real(kind=dp), allocatable, target :: scalar_array(:,:)
  
  integer  :: n, m
  integer :: nr
  type(ZirkaABC_t), POINTER :: ABCParams
  allocate(ABCParams )

  call GetConstRealArray(model % materials(2) % values, BHasc, 'Ascending BH curve', found) 
  if(.not. Found) call fatal('ZirkaTest', 'Ascending BH curve not found')
  call GetConstRealArray(model % materials(2) % values, BHSingle, 'Single valued BH curve', found)
  if(.not. Found) call fatal('ZirkaTest', 'Single valued BH curve not found')

  call info('ZirkaTest','Running zirkatest')
  rps_array => ListGetConstRealArray(GetSolverParams(), 'zirka init sequence', found)
  rps(1:size(rps_array,1)) => rps_array(:,1)

  mc = MasterCurve_t(InitSplineLoop(BHasc, BHsingle), &
      ABCparams, n_cachesubsample = 2)

  nr = size(rps,1)

  DO n = 1,nr-1
    CALL mc % drive(rps(n), cache=.false.)
  END DO
  CALL mc % drive(rps(nr), cache=.true.)

  rps => null()

  m = 1
  BLOCK
    integer :: i_r
    ! rps_array => ListGetConstRealArray(GetSolverParams(), 'rps_3', found)
    scalar_array = readdat("scalars_a.dat", 7, 9)
  call info('ZirkaTest','Running zirkatest')
    rps(1:size(scalar_array,1)) => scalar_array(:,3)
    H(1:size(scalar_array,1)) => scalar_array(:,4)

    do i_r = 2, ubound(rps,1)
      call mc % drive(rps(i_r-1), cache = .true.)
        associate (B_low => rps(i_r-1),  B_up => rps(i_r))
          call print_interval(mc, B_low, B_up, m, H(i_r-1))
        end associate
    end do
    i_r = ubound(rps,1)
    call mc % drive(rps(i_r), cache = .true.)
    associate (B_low => rps(i_r),  B_up => rps(i_r))
      call print_interval(mc, B_low, B_up, m, H(i_r))
    end associate
    variable => VariableGet(solver % variable, "testvar")
    variable % values = abs(H(i_r)-mc % eval(rps(i_r), cached = .true.))/abs(H(i_r))
    variable % norm = abs(H(i_r)-mc % eval(rps(i_r), cached = .true.))/abs(H(i_r))
  END BLOCK

end subroutine ! }}}
function readBH(f) result(BH) ! {{{
  real(kind=dp), allocatable :: BH(:,:)
  character(len=*) :: f
  real*8 :: lenf
  integer :: iu, len, n

  open(newunit=iu, file=f, action='read', status='old')

  read(iu, *) len

  allocate(BH(1:len,2))

  do n = 1,len
    read (iu,*) BH(n,1:2)
  end do
  close(iu)

end function ! }}}

function readdat(f, m, n) result(table) ! {{{
  real(kind=dp), allocatable :: table(:,:)
  character(len=*) :: f
  real*8 :: lenf
  integer :: iu, m, n, k

  open(newunit=iu, file=f, action='read', status='old')

  allocate(table(m,n))

  do k = 1,m
    read (iu,*) table(k,1:n)
  end do
  close(iu)

end function ! }}}
END SUBROUTINE ZirkaTest ! }}}
!-------------------------------------------------------------------------------
