MODULE PermonInterface
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_INT, C_LOC, C_NULL_PTR
  USE Types, ONLY: dp
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: permon_solve
  PUBLIC :: get_c_ptr_dp

  INTERFACE
    SUBROUTINE permon_solve(rows, cols, vals, nrows, ncols, b_ptr, limits_ptr, x_ptr) BIND(C, NAME="permon_solve")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_INT
      TYPE(C_PTR), VALUE :: rows
      TYPE(C_PTR), VALUE :: cols
      TYPE(C_PTR), VALUE :: vals
      INTEGER(C_INT), VALUE :: nrows
      INTEGER(C_INT), VALUE :: ncols
      TYPE(C_PTR), VALUE :: b_ptr
      TYPE(C_PTR), VALUE :: limits_ptr
      TYPE(C_PTR), VALUE :: x_ptr
    END SUBROUTINE permon_solve
  END INTERFACE

END MODULE PermonInterface