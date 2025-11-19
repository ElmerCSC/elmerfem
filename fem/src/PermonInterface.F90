MODULE PermonInterface
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_INT, C_LOC, C_NULL_PTR
  USE Types, ONLY: dp
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: permon_solve, permon_init, permon_finalize

  INTERFACE
    SUBROUTINE permon_solve(rows, cols, vals, nrows, ncols, b_ptr, limits_ptr, x_ptr, bound) BIND(C, NAME="permon_solve")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_INT
      TYPE(C_PTR), VALUE :: rows
      TYPE(C_PTR), VALUE :: cols
      TYPE(C_PTR), VALUE :: vals
      INTEGER(C_INT), VALUE :: nrows
      INTEGER(C_INT), VALUE :: ncols
      TYPE(C_PTR), VALUE :: b_ptr
      TYPE(C_PTR), VALUE :: limits_ptr
      TYPE(C_PTR), VALUE :: x_ptr
      INTEGER(C_INT), VALUE :: bound
    END SUBROUTINE permon_solve

    SUBROUTINE permon_init() BIND(C, NAME="permon_init")
    END SUBROUTINE permon_init

    SUBROUTINE permon_finalize() BIND(C, NAME="permon_finalize")
    END SUBROUTINE permon_finalize
  END INTERFACE

END MODULE PermonInterface