
#include "huti_fdefs.h"

MODULE huti_sfe
  USE huti_cg
  USE huti_cgs
  USE huti_bicgstab
  USE huti_bicgstab_2
  USE huti_gmres
  USE huti_qmr
  USE huti_tfqmr
  USE huti_aux
  USE huti_interfaces
  IMPLICIT NONE

  INTERFACE
  
     FUNCTION zdotc(n, x, xinc, y, yinc) RESULT(res)
       IMPLICIT NONE
       INTEGER :: n, xinc, yinc
       DOUBLE COMPLEX :: x(*), y(*)
       DOUBLE COMPLEX :: res
     END FUNCTION zdotc

     FUNCTION zdotu(n, x, xinc, y, yinc) RESULT(res)
       IMPLICIT NONE
       INTEGER :: n, xinc, yinc
       DOUBLE COMPLEX :: x(*), y(*)
       DOUBLE COMPLEX :: res
     END FUNCTION zdotu

     FUNCTION cdotc(n, x, xinc, y, yinc) RESULT(res)
       IMPLICIT NONE
       INTEGER :: n, xinc, yinc
       COMPLEX :: x(*), y(*)
       COMPLEX :: res
     END FUNCTION cdotc

     FUNCTION cdotu(n, x, xinc, y, yinc) RESULT(res)
       IMPLICIT NONE
       INTEGER :: n, xinc, yinc
       COMPLEX :: x(*), y(*)
       COMPLEX :: res
     END FUNCTION cdotu

     FUNCTION sdot(n, x, xinc, y, yinc) RESULT(res)
       IMPLICIT NONE
       INTEGER :: n, xinc, yinc
       REAL :: x(*), y(*), res
     END FUNCTION sdot

     FUNCTION ddot(n, x, xinc, y, yinc) RESULT(res)
       IMPLICIT NONE
       INTEGER :: n, xinc, yinc
       DOUBLE PRECISION :: x(*), y(*), res
     END FUNCTION ddot

     FUNCTION dnrm2(n, x, xinc) RESULT(res)
       IMPLICIT NONE
       INTEGER :: n, xinc
       DOUBLE PRECISION :: x(*), res
     END FUNCTION dnrm2

     FUNCTION snrm2(n, x, xinc) RESULT(res)
       IMPLICIT NONE
       INTEGER :: n, xinc
       REAL :: x(*), res
     END FUNCTION snrm2

     FUNCTION dznrm2(n, x, xinc) RESULT(res)
       IMPLICIT NONE
       INTEGER :: n, xinc
       DOUBLE COMPLEX :: x(*)
       DOUBLE PRECISION :: res
     END FUNCTION dznrm2

     FUNCTION scnrm2(n, x, xinc) RESULT(res)
       IMPLICIT NONE
       INTEGER :: n, xinc
       COMPLEX :: x(*)
       REAL :: res
     END FUNCTION scnrm2

  END INTERFACE

CONTAINS
  
  SUBROUTINE huti_s_cg(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_s ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_s ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_s ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_s ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_s ), POINTER :: normfun
    PROCEDURE( stopc_iface_s ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_sdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_sdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => sdot
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => snrm2
    END IF

    CALL huti_scgsolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_s_cg


  SUBROUTINE huti_s_tfqmr(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_s ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_s ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_s ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_s ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_s ), POINTER :: normfun
    PROCEDURE( stopc_iface_s ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_sdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_sdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => sdot
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => snrm2
    END IF

    CALL huti_stfqmrsolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_s_tfqmr


  SUBROUTINE huti_s_cgs(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_s ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_s ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_s ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_s ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_s ), POINTER :: normfun
    PROCEDURE( stopc_iface_s ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_sdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_sdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => sdot
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => snrm2
    END IF

    CALL huti_scgssolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_s_cgs


  SUBROUTINE huti_s_qmr(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_s ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_s ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_s ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_s ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_s ), POINTER :: normfun
    PROCEDURE( stopc_iface_s ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_sdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_sdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => sdot
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => snrm2
    END IF

    CALL huti_sqmrsolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_s_qmr


  SUBROUTINE huti_s_bicgstab(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_s ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_s ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_s ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_s ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_s ), POINTER :: normfun
    PROCEDURE( stopc_iface_s ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_sdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_sdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => sdot
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => snrm2
    END IF

    CALL huti_sbicgstabsolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_s_bicgstab


  SUBROUTINE huti_s_gmres(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_s ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_s ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_s ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_s ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_s ), POINTER :: normfun
    PROCEDURE( stopc_iface_s ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_sdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_sdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => sdot
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => snrm2
    END IF

    CALL huti_sgmressolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_s_gmres


  SUBROUTINE huti_s_bicgstab_2(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_s ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_s ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_s ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_s ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_s ), POINTER :: normfun
    PROCEDURE( stopc_iface_s ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_sdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_sdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => sdot
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => snrm2
    END IF

    CALL huti_sbicgstab_2solv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, &
         &dpar, work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_s_bicgstab_2


  SUBROUTINE huti_d_cg(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_d ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_d ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_d ), POINTER :: normfun
    PROCEDURE( stopc_iface_d ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    DOUBLE PRECISION, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    DOUBLE PRECISION, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_ddummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_ddummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => ddot
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => dnrm2
    END IF

    CALL huti_dcgsolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_d_cg


  SUBROUTINE huti_d_cgs(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_d ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_d ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_d ), POINTER :: normfun
    PROCEDURE( stopc_iface_d ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    DOUBLE PRECISION, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    DOUBLE PRECISION, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_ddummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_ddummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => ddot
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => dnrm2
    END IF

    CALL huti_dcgssolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_d_cgs


  SUBROUTINE huti_d_bicgstab(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_d ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_d ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_d ), POINTER :: normfun
    PROCEDURE( stopc_iface_d ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    DOUBLE PRECISION, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    DOUBLE PRECISION, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_ddummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_ddummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => ddot
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => dnrm2
    END IF

    CALL huti_dbicgstabsolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_d_bicgstab


  SUBROUTINE huti_d_tfqmr(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_d ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_d ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_d ), POINTER :: normfun
    PROCEDURE( stopc_iface_d ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    DOUBLE PRECISION, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    DOUBLE PRECISION, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_ddummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_ddummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => ddot
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => dnrm2
    END IF

    CALL huti_dtfqmrsolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_d_tfqmr


  SUBROUTINE huti_d_qmr(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_d ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_d ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_d ), POINTER :: normfun
    PROCEDURE( stopc_iface_d ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    DOUBLE PRECISION, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    DOUBLE PRECISION, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_ddummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_ddummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => ddot
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => dnrm2
    END IF

    CALL huti_dqmrsolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_d_qmr


  SUBROUTINE huti_d_gmres(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_d ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_d ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_d ), POINTER :: normfun
    PROCEDURE( stopc_iface_d ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    DOUBLE PRECISION, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    DOUBLE PRECISION, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_ddummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_ddummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => ddot
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => dnrm2
    END IF

    CALL huti_dgmressolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_d_gmres


  SUBROUTINE huti_d_bicgstab_2(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_d ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_d ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_d ), POINTER :: normfun
    PROCEDURE( stopc_iface_d ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    DOUBLE PRECISION, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    DOUBLE PRECISION, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_ddummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_ddummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => ddot
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => dnrm2
    END IF

    CALL huti_dbicgstab_2solv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, &
         & dpar, work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_d_bicgstab_2


  SUBROUTINE huti_c_cg(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_c ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_c ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_c ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_c ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_c ), POINTER :: normfun
    PROCEDURE( stopc_iface_c ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    COMPLEX, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    COMPLEX, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_cdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_cdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => cdotu
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => scnrm2
    END IF

    CALL huti_ccgsolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_c_cg


  SUBROUTINE huti_c_cgs(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_c ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_c ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_c ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_c ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_c ), POINTER :: normfun
    PROCEDURE( stopc_iface_c ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    COMPLEX, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    COMPLEX, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_cdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_cdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => cdotu
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => scnrm2
    END IF

    CALL huti_ccgsolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_c_cgs


  SUBROUTINE huti_c_bicgstab(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_c ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_c ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_c ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_c ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_c ), POINTER :: normfun
    PROCEDURE( stopc_iface_c ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    COMPLEX, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    COMPLEX, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_cdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_cdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => cdotu
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => scnrm2
    END IF

    CALL huti_cbicgstabsolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_c_bicgstab


  SUBROUTINE huti_c_qmr(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_c ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_c ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_c ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_c ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_c ), POINTER :: normfun
    PROCEDURE( stopc_iface_c ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    COMPLEX, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    COMPLEX, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_cdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_cdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => cdotu
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => scnrm2
    END IF

    CALL huti_cqmrsolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_c_qmr


  SUBROUTINE huti_c_tfqmr(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_c ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_c ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_c ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_c ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_c ), POINTER :: normfun
    PROCEDURE( stopc_iface_c ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    COMPLEX, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    COMPLEX, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work
 
    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_cdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_cdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => cdotu
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => scnrm2
    END IF

    CALL huti_ctfqmrsolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_c_tfqmr


  SUBROUTINE huti_c_gmres(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE    
    PROCEDURE( mv_iface_c ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_c ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_c ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_c ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_c ), POINTER :: normfun
    PROCEDURE( stopc_iface_c ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    COMPLEX, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    COMPLEX, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_cdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_cdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => cdotc
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => scnrm2
    END IF

    CALL huti_cgmressolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_c_gmres


  SUBROUTINE huti_c_bicgstab_2(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_c ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_c ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_c ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_c ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_c ), POINTER :: normfun
    PROCEDURE( stopc_iface_c ), POINTER :: mstopfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    COMPLEX, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    COMPLEX, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_cdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_cdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => cdotu
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => scnrm2
    END IF

    CALL huti_cbicgstab_2solv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, &
         & dpar, work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_c_bicgstab_2

  
  SUBROUTINE huti_z_cg(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    IMPLICIT NONE
    PROCEDURE( mv_iface_z ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_z ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_z ), POINTER :: normfun
    PROCEDURE( stopc_iface_z ), POINTER :: mstopfun

    INTEGER :: ndim, wrkdim
    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    DOUBLE COMPLEX, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    DOUBLE COMPLEX, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_zdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_zdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => zdotu
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => dznrm2
    END IF

    CALL huti_zcgsolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_z_cg


  SUBROUTINE huti_z_cgs(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    
    PROCEDURE( mv_iface_z ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_z ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_z ), POINTER :: normfun
    PROCEDURE( stopc_iface_z ), POINTER :: mstopfun

    INTEGER :: ndim, wrkdim
    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    DOUBLE COMPLEX, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    DOUBLE COMPLEX, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_zdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_zdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => zdotu
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => dznrm2
    END IF

    CALL huti_zcgssolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_z_cgs


  SUBROUTINE huti_z_bicgstab(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    
    PROCEDURE( mv_iface_z ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_z ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_z ), POINTER :: normfun
    PROCEDURE( stopc_iface_z ), POINTER :: mstopfun

    INTEGER :: ndim, wrkdim
    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    DOUBLE COMPLEX, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    DOUBLE COMPLEX, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_zdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_zdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => zdotu
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => dznrm2
    END IF

    CALL huti_zbicgstabsolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_z_bicgstab


  SUBROUTINE huti_z_qmr(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    
    PROCEDURE( mv_iface_z ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_z ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_z ), POINTER :: normfun
    PROCEDURE( stopc_iface_z ), POINTER :: mstopfun

    INTEGER :: ndim, wrkdim
    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    DOUBLE COMPLEX, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    DOUBLE COMPLEX, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_zdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_zdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => zdotu
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => dznrm2
    END IF

    CALL huti_zqmrsolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_z_qmr


  SUBROUTINE huti_z_tfqmr(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    
    PROCEDURE( mv_iface_z ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_z ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_z ), POINTER :: normfun
    PROCEDURE( stopc_iface_z ), POINTER :: mstopfun

    INTEGER :: ndim, wrkdim
    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    DOUBLE COMPLEX, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    DOUBLE COMPLEX, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_zdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_zdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => zdotu
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => dznrm2
    END IF

    CALL huti_ztfqmrsolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_z_tfqmr


  SUBROUTINE huti_z_gmres(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    
    PROCEDURE( mv_iface_z ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_z ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_z ), POINTER :: normfun
    PROCEDURE( stopc_iface_z ), POINTER :: mstopfun

    INTEGER :: ndim, wrkdim
    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    DOUBLE COMPLEX, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    DOUBLE COMPLEX, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_zdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_zdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => zdotc
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => dznrm2
    END IF

    CALL huti_zgmressolv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, dpar, &
         & work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_z_gmres


  SUBROUTINE huti_z_bicgstab_2(xvec, rhsvec, ipar, dpar, work, matvecsubr, &
       & pcondlsubr, pcondrsubr, dotprodfun, normfun, mstopfun)
    
    PROCEDURE( mv_iface_z ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_z ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_z ), POINTER :: normfun
    PROCEDURE( stopc_iface_z ), POINTER :: mstopfun

    INTEGER :: ndim, wrkdim
    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    DOUBLE COMPLEX, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    DOUBLE COMPLEX, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    IF(.NOT. ASSOCIATED(pcondrsubr) ) THEN
       pcondrsubr => huti_zdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(pcondlsubr) ) THEN
       pcondlsubr => huti_zdummy_pcondfun
    END IF
    IF(.NOT. ASSOCIATED(dotprodfun) ) THEN
       dotprodfun => zdotu
    END IF
    IF(.NOT. ASSOCIATED(normfun) ) THEN
       normfun => dznrm2
    END IF

    CALL huti_zbicgstab_2solv(HUTI_NDIM, HUTI_WRKDIM, xvec, rhsvec, ipar, &
         & dpar, work, matvecsubr, pcondlsubr, pcondrsubr, dotprodfun, &
         & normfun, mstopfun)
  END SUBROUTINE huti_z_bicgstab_2


END MODULE huti_sfe
