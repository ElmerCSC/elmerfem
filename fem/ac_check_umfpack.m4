dnl -*- mode: autoconf -*-
dnl
dnl This file is part of the LifeV Libraries.
dnl
dnl Author: Christophe Prud'homme (christophe.prudhomme@epfl.ch)
dnl
dnl Copyright (C) 2004 EPFL, INRIA, Politecnico di Milano
dnl
dnl Distributed under the GPL(GNU Public License):
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2 of the License, or
dnl (at your option) any later version.
dnl 
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl 
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
dnl
AC_DEFUN([AC_CHECK_UMFPACK],
[
 if test -z "$1"; then
  umfpack_ver=4
  umfpack_subver=3
  umfpack_minorrel=0
 else
  umfpack_minor=`echo "$1" | sed -e 's#[0-9]\+\.[0-9]\+\.\([0-9]\+\).*#\1#'`
  umfpack_subver=`echo "$1" | sed -e 's#[0-9]\+\.\([0-9]\+\).*#\1#'`
  umfpack_ver=`echo "$1" | sed -e 's#^\([0-9]\+\)\..*#\1#'`
 fi
 umfpack_version="${umfpack_ver}.${umfpack_subver}.${umfpack_minor}"

 AC_MSG_CHECKING([for amd include file directory])
 ac_mpi_includedirs="/usr/include/amd /usr/include/umfpack $REPOSITORY/$ARCH/include/UMFPACK /usr/include /usr/local/include/amd /usr/local/include/umfpack /usr/local/include"
 AC_FIND_FILE("amd.h", $ac_mpi_includedirs, amd_includedir)
 if test "x${amd_includedir}" != "x" -a "x${amd_includedir}" != "xNO"; then
  CPPFLAGS="-I${amd_includedir} ${CPPFLAGS}"
 fi
 AC_MSG_RESULT([${amd_includedir}])  
 #AC_CHECK_HEADERS([amd.h])

 ac_save_cppflags=${CPPFLAGS}

 AC_MSG_CHECKING([for umfpack include file directory])
 ac_mpi_includedirs="/usr/include/umfpack $REPOSITORY/$ARCH/include/UMFPACK /usr/include /usr/local/include/umfpack /usr/local/include"
 AC_FIND_FILE("umfpack.h", $ac_mpi_includedirs, umfpack_includedir)
 if test "x${umfpack_includedir}" != "x" -a "x${umfpack_includedir}" != "xNO"; then
  CPPFLAGS="-I${umfpack_includedir} ${CPPFLAGS}"
 fi
 AC_CHECK_HEADERS([umfpack.h], [], [CPPFLAGS="${ac_save_cppflags}"], [])

# if test ${umfpack_includedir} != ${amd_includedir}
# then
#   CPPFLAGS="-I${umfpack_includedir} ${CPPFLAGS}"
#   AC_MSG_RESULT([${umfpack_includedir}])
# fi

 AC_MSG_CHECKING([AMD libraries presence])
 ac_amd_libdirs="$REPOSITORY/$ARCH/lib$AFFIX/UMFPACK /usr/lib /usr/local/lib /usr/lib/amd /usr/local/lib/amd /usr/lib/umfpack /usr/local/lib/umfpack"
 AC_FIND_FILE("libamd.a", $ac_amd_libdirs, ac_amd_libdir)
 if test "x${ac_amd_libdir}" != "x" -a "x${ac_amd_libdir}" != "xNO"; then
  amd_library=$ac_amd_libdir/libamd.a
  LDFLAGS="${LDFLAGS} -L${ac_amd_libdir}"
 fi
 AC_MSG_RESULT([$amd_library])

 AC_MSG_CHECKING([UMFPACK libraries presence])
 ac_umfpack_libdirs="$REPOSITORY/$ARCH/lib$AFFIX/UMFPACK /usr/lib /usr/local/lib /usr/lib/umfpack /usr/local/lib/umfpack"
 AC_FIND_FILE("libumfpack.a", $ac_umfpack_libdirs, ac_umfpack_libdir)
 if test "x${ac_umfpack_libdir}" != "x" -a "x${ac_umfpack_libdir}" != "xNO"; then
  umfpack_library=$ac_umfpack_libdir/libumfpack.a
  LDFLAGS="${LDFLAGS} -L${ac_umfpack_libdir}"
 fi
 AC_MSG_RESULT([$umfpack_library])
 
if test "$umfpack_ver" = "4"; then
  AC_CHECK_LIB(amd,amd_postorder,[umfpack_libs="-L${ac_amd_libdir} -lamd"],[echo "Please install umfpack/amd."; ],${lapack_libs})
  AC_CHECK_LIB(umfpack,umfpack_di_solve,[umfpack_libs="-L${ac_umfpack_libdir} -lumfpack ${umfpack_libs}"],[echo "Please install umfpack/amd."; ],-lamd ${lapack_libs})
 else
  AC_CHECK_LIB(umfpack,umfpack_solve,[umfpack_libs="-L${ac_umfpack_libdir} -lumfpack"],[echo "Please install umfpack."; ],${lapack_libs})
 fi
 AC_SUBST(umfpack_libs)
 AC_MSG_RESULT([$umfpack_libs])
])
