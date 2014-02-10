dnl
dnl Partly copied from http://www.hlrs.de/people/keller/configure_scripts.html
dnl
dnl TODO: test fortran-mpi linkage...
dnl

AC_DEFUN([ACX_MPI], [
  AC_PREREQ(2.50) dnl for AC_LANG_CASE

  dnl ac_mpi_save_LIBS=$LIBS
  acx_mpi_ok=disabled

  acx_mpi_try_c_compile=yes
  case  $host in 
   rs6000-ibm-aix* | powerpc-ibm-aix*)
       acx_mpi_try_c_compile=no
   ;;
  esac


  dnl Letting user specify MPI-Library directories
  AC_ARG_WITH(mpi,
    [  --with-mpi[=yes]  Use mpi (by default disabled)],
       acx_mpi_ok=no, acx_mpi_ok=disabled)

  AC_MSG_CHECKING([for mpi-directory])
  AC_ARG_WITH(mpi_dir,
    [  --with-mpi-dir=MPIDIR   give the path for MPI []],
    acx_mpi_ok=no; mpi_dir="$withval", mpi_dir="")
  AC_MSG_RESULT([$mpi_dir])
  AC_SUBST([mpi_dir])

  AC_MSG_CHECKING([for mpi-lib-directory])
  AC_ARG_WITH(mpi_lib_dir,
    [  --with-mpi-lib-dir=dir  give the path for MPI-libraries [MPI_DIR/lib]],
    acx_mpi_ok=no; mpi_lib_dir="$withval", 
    if test "$mpi_dir" != ""; then
     mpi_lib_dir="$mpi_dir/lib"
    fi
    )
  AC_MSG_RESULT([$mpi_lib_dir])
  AC_SUBST([mpi_lib_dir])

  AC_MSG_CHECKING([for mpi-inc-directory])
  AC_ARG_WITH(mpi_inc_dir,
    [  --with-mpi-inc-dir=dir  give the path for MPI-include-files [MPI_DIR/include]],
    acx_mpi_ok=no; mpi_inc_dir="$withval", 
    if test "$mpi_dir" != ""; then
     mpi_inc_dir="$mpi_dir/include"
    fi
  )
  AC_MSG_RESULT([$mpi_inc_dir])
  AC_SUBST([mpi_inc_dir])

  AC_MSG_CHECKING([for mpi-library])
  AC_ARG_WITH(mpi_library,
    [  --with-mpi-library=library  give MPI-libraray-files []],
    acx_mpi_ok=no; mpi_library="$withval", mpi_library="")
  AC_MSG_RESULT([$mpi_library])
  AC_SUBST([mpi_library])

  AC_MSG_CHECKING([for mpi-include])
  AC_ARG_WITH(mpi_include,
    [  --with-mpi-include=include  give MPI-include-files []],
    acx_mpi_ok=no; mpi_include="$withval", mpi_include="")
  AC_MSG_RESULT([$mpi_include])
  AC_SUBST([mpi_include])

  if test "$acx_mpi_ok" != disabled; then

    # MPI-Library name (depends on variables $mpi_lib_dir and user-defined argument PACX_SIGNAL on IBMs) 
    if test "mpi_library" == ""; then

      AC_MSG_CHECKING([for MPI library])
      case "$host" in
        *-ibm-aix*)                # IBM/SP2 machines
          # checking whether to use signal-based MPI

          AC_MSG_CHECKING([whether to use signal-based MPI library])
          AC_MSG_RESULT([$PACX_SIGNAL])
          if test "x$PACX_SIGNAL" = "xyes" ; then
            if test -f "$mpi_lib_dir/libmpi.a" ; then
              lib_mpi="-lmpi"
            elif test -f "$mpi_lib_dir/libmpi.so" ; then
              lib_mpi="-lmpi"
            elif test -f "$mpi_lib_dir/libmpich.a" ; then
              lib_mpi="-lmpich"
            else
              AC_MSG_ERROR([neither libmpi nor libmpich found; check path for MPI package first...])
            fi
          else
            if test -f "$mpi_lib_dir/libmpi_r.a" ; then
               lib_mpi="-lmpi_r -bautoload:\"$mpi_lib_dir/libmpi_r.a\(mpifort64_r.o\)\""
            elif test -f "/usr/lpp/ppe.poe/lib/libmpi_r.a" ; then
               mpi_dir="/usr/lpp/ppe.poe/"
               mpi_lib_dir="$mpi_dir/lib"
               mpi_inc_dir="$mpi_dir/include/thread64"
               lib_mpi="-lmpi_r -bautoload:\"$mpi_lib_dir/libmpi_r.a\(mpifort64_r.o\)\""
            else
               AC_MSG_ERROR([libmpi_r not found; check path for MPI package...])
            fi
          fi
          AC_MSG_RESULT(found $lib_mpi)
        ;;
        *)                         # All other machines
          if test -f "$mpi_lib_dir/libmpi.a" ; then
            lib_mpi="-lmpi"
          elif test -f "$mpi_lib_dir/libmpi.so" ; then
            lib_mpi="-lmpi"
          elif test -f "$mpi_lib_dir/libmpich.a" ; then
            lib_mpi="-lmpich"
          else
            AC_MSG_ERROR([neither libmpi nor libmpich found; check path for MPI package first...])
          fi
          AC_MSG_RESULT(found $lib_mpi)
        ;;
      esac
    else
     lib_mpi=$mpi_library
    fi

    if test "$MPI_LIBS" == ""; then
      MPI_LIBS="-L$mpi_lib_dir $lib_mpi"
    fi

    AC_SUBST(lib_mpi)

    # Compilation of a MPI program (depends on above macro)
    if test "$acx_mpi_try_c_compile" != "no"; then
	AC_LANG_PUSH(Fortran 77)
	AC_MSG_CHECKING([for compilation of an MPI program])
	old_FFLAGS=${FFLAGS}
	old_LIBS=${LIBS}
	FFLAGS="$mpi_include"
        if test "$mpi_inc_dir" != ""; then
	  FFLAGS="-I$mpi_inc_dir $FFLAGS"
        fi
	LIBS="$lib_mpi $SYS_LDFLAGS"
        if test "$mpi_inc_dir" != ""; then
	  LIBS="-L$mpi_lib_dir $LIBS"
        fi
	AC_COMPILE_IFELSE(
	[        PROGRAM test
                 INCLUDE "mpif.h"
                 CALL MPI_Finalize(ierr)
                 END
	],[AC_MSG_RESULT([seems ok])
	    AC_DEFINE([HAVE_MPI],[1],[...])
	    acx_mpi_ok=yes],
	[  AC_MSG_ERROR([MPI not found; check paths for MPI package first...])])

	FFLAGS=${old_FFLAGS}
	LIBS=${old_LIBS}
	AC_LANG_POP(Fortran 77)
    else
  	acx_mpi_ok=yes
    fi

    MPI_INCLUDE=""
    MPI_INCLUDE_DIR=""

    if test "$mpi_include" == ""; then
      if test "$acx_mpi_ok" != yes; then
        AC_CHECK_FILE($mpi_inc_dir/mpif.h, 
        [acx_mpif_h_found=yes
         MPI_INCLUDE_DIR=$mpi_inc_dir], 
        [acx_mpif_h_found=no
         MPI_INCLUDE_DIR=""])
      else
        acx_mpif_h_found=yes
        MPI_INCLUDE_DIR=$mpi_inc_dir 
      fi
    else
       MPI_INCLUDE=$mpi_include
    fi

  else  
     # use local mpif.h
     acx_mpif_h_found=no
     MPI_LIBS=""
  fi

  AC_SUBST(MPI_LIBS)
  AC_SUBST(MPI_INCLUDE_DIR)

])# ACX_MPI

# fixme: implement all mpi compilers when needed
AC_DEFUN([ACX_MPI_COMPILERS],
[
MPI_F90=
acx_mpi_f90_compilers="mpf90 mpxlf90"

for c in $acx_mpi_f90_compilers; do
	AC_CHECK_PROG(MPI_F90, $c, $c)
	if test "$MPI_F90" = "$c"; then
		break
	fi
done

AC_SUBST(MPI_F90)
])

# Macro AC_CHECK_MPI_VERSION to check for version number of MPI

AC_DEFUN([AC_CHECK_MPI_VERSION],
[AC_CACHE_CHECK([version of MPI implementation], [ac_cv_mpi_version],
  [ old_CFLAGS=${CFLAGS}
    old_LIBS=${LIBS}
    CFLAGS="-I$mpi_inc_dir"
    LIBS="-L$mpi_lib_dir $lib_mpi $SYS_LDFLAGS"
    AC_TRY_RUN([
#include <stdio.h>
#include <mpi.h>
int main (){
  FILE * f =fopen("conftestval", "w");
  if (!f) return 1;
#ifdef MPI_VERSION
  fprintf (f, "%d.%d\n", MPI_VERSION, MPI_SUBVERSION);
#else
  fprintf (f, "unkown\n");
#endif
  return 0;
}],
  [
    VAL=`cat ./conftestval`
    if test "x$VAL" = "x" ; then
      ac_cv_mpi_version="unkown"
    else
      ac_cv_mpi_version="$VAL"
    fi
    ], [ac_cv_mpi_version="unknown"])
    rm -f conftestval
  ])
])



dnl Optional MPI-Datatypes for Fortran
dnl This macro is a little more complex; it checks for the availability of the optional Fortran datatype, like MPI_INTEGER1 or MPI_REAL8.
dnl Usage: AC_CHECK_FORTRAN_MPI_DATATYPE (DATATYPE, [ACTION-IF-FOUND, [ACTION-IF-NOT-FOUND])

AC_DEFUN([AC_CHECK_FORTRAN_MPI_DATATYPE],
  [
  AC_CACHE_CHECK([whether MPI has (opt.) Fortran [$1]], [ac_cv_have_mpi_fortran_[$1]],
  [
    AC_LANG_PUSH([Fortran 77])
    dnl This gets hard -- we're looking for a F77 compiler, which supports MPI!
    dnl We have to jump through a hoop, to get it running on all systems!
    ac_link=""
    if test -z "$ac_link" ; then
      if test -x $mpi_dir/bin/mpif77 ; then
        ac_link='$mpi_dir/bin/mpif77 -o conftest conftest.f'
      fi
    fi
    if test -z "$ac_link" ; then
      if test -x $mpi_dir/bin/mpif90 ; then
        ac_link='$mpi_dir/bin/mpif90 -o conftest conftest.f'
      fi
    fi
    if test -z "$ac_link" ; then
      if test -x $mpi_dir/bin/mpf77 ; then
        ac_link='$mpi_dir/bin/mpf77 -o conftest conftest.f'
      fi
    fi

    if test -z "$ac_link" ; then
      if test -x $mpi_dir/bin/mpf90 ; then
        ac_link='$mpi_dir/bin/mpf90 -o conftest conftest.f'
      fi
    fi

    if test -z "$ac_link" ; then
      dnl This is our last resort -- do it by hand
      dnl This might mean, that the test-program will not run on the cpu, configure runs on!
      ac_link='$F77 $FFLAGS -I${mpi_inc_dir} -o conftest conftest.f -L${mpi_lib_dir} ${lib_mpi} $SYS_LDFLAGS'
    fi

    AC_TRY_RUN([
          implicit none
          include 'mpif.h'

          open (unit=44, file='conftestval')
          write (44,*) ' ', $1
          close (44)
          end
    ], [
    dnl MPIch defines those Variables, but sets them to 0, if they are not available.
      VAL=`cat ./conftestval`
      if test $VAL = 0 ; then
        ac_cv_have_mpi_fortran_$1="no"
      else
        ac_cv_have_mpi_fortran_$1="yes"
      fi
    ], [ac_cv_have_mpi_fortran_$1="no"])
    AC_LANG_POP()
    rm -f conftest conftest.f conftestval
  ])
])

