MODULE PermonInterface
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_INT, C_LOC, C_NULL_PTR, C_INTPTR_T
  USE Types, ONLY: dp
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: permon_solve, permon_init, permon_finalize

  INTERFACE
    SUBROUTINE permon_solve(rows, cols, vals, nrows, b_ptr, limits_ptr, x_ptr, bound, &
    gdofs_cptr, owner_cptr, comm) BIND(C, NAME="permon_solve")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_INT
      TYPE(C_PTR), VALUE :: rows
      TYPE(C_PTR), VALUE :: cols
      TYPE(C_PTR), VALUE :: vals
      INTEGER(C_INT), VALUE :: nrows
      TYPE(C_PTR), VALUE :: b_ptr
      TYPE(C_PTR), VALUE :: limits_ptr
      TYPE(C_PTR), VALUE :: x_ptr
      INTEGER(C_INT), VALUE :: bound
      TYPE(C_PTR), VALUE :: gdofs_cptr
      TYPE(C_PTR), VALUE :: owner_cptr
      INTEGER(C_INT), VALUE :: comm
    END SUBROUTINE permon_solve

    SUBROUTINE permon_init() BIND(C, NAME="permon_init")
    END SUBROUTINE permon_init

    SUBROUTINE permon_finalize() BIND(C, NAME="permon_finalize")
    END SUBROUTINE permon_finalize
  END INTERFACE

END MODULE PermonInterface