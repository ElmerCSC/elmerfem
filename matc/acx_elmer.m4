dnl 
dnl Elmer specific M4sh macros 
dnl
dnl @version $Id: acx_elmer.m4,v 1.72 2005/08/26 09:54:15 vierinen Exp $
dnl @author juha.vierinen@csc.fi 5/2005
dnl

dnl
dnl define host variable
dnl
AC_DEFUN([ACX_HOST],
[
AC_REQUIRE([AC_CANONICAL_HOST])
AC_REQUIRE([AC_CANONICAL_TARGET])

if test -z "$host"; then
  host=unknown
fi
canonical_host_type=$host
if test "$host" = unknown; then
  AC_MSG_ERROR([configuring for unknown system type, your build will most likely be screwed.])
fi

case "$canonical_host_type" in
  *-*-darwin*)
	LDFLAGS="$LDFLAGS -L/sw/lib"
	CFLAGS="$CFLAGS -I/sw/include"
	CXXFLAGS="$CXXFLAGS -I/sw/include"
  ;;
esac

AC_SUBST(canonical_host_type)
])


AC_DEFUN([ACX_PROG_AR],[
# fixme: do something more intelligent here
AC_MSG_CHECKING([for ar])
if test -z "$AR"; then
        AR='ar'
	acx_ar_userdefined=no
else
	acx_ar_userdefined=yes
fi
AC_SUBST(AR)
AC_MSG_RESULT($AR)

AC_MSG_CHECKING([for ar flags])
if test -z "$ARFLAGS"; then
  ARFLAGS="cru"
fi
AC_SUBST(ARFLAGS)
AC_MSG_RESULT($ARFLAGS)
])

dnl
dnl Default optimization flags
dnl
AC_DEFUN([ACX_DEBUG],
[
acx_debug=no

AC_ARG_WITH(debug,
	   [AC_HELP_STRING([--with-debug], [Use debugging flags (environment variables override)])],[acx_debug=yes])

AC_MSG_CHECKING([for default compilation flags])

if test "$acx_debug" = yes; then
AC_MSG_RESULT([debugging])
  if test "$CFLAGS" = ""; then
	CFLAGS="-g"
  fi

  if test "$FCFLAGS" = ""; then
	FCFLAGS="-g"
  fi

  if test "$FFLAGS" = ""; then
	FFLAGS="-g"
  fi

  if test "$CXXFLAGS" = ""; then
	CXXFLAGS="-g"
  fi
else
AC_MSG_RESULT([optimized])
  case "$canonical_host_type" in
     rs6000-ibm-aix* | powerpc-ibm-aix*)
	dnl better optimization on ibm
	ACX_OPT_FLAG="-O -qmaxmem=-1 -qstrict"
     ;;
     *)
	ACX_OPT_FLAG="-O"
     ;;
  esac

  if test "$CFLAGS" = ""; then
	CFLAGS=$ACX_OPT_FLAG
  fi

  if test "$FCFLAGS" = ""; then
	FCFLAGS=$ACX_OPT_FLAG
  fi

  if test "$FFLAGS" = ""; then
	FFLAGS=$ACX_OPT_FLAG
  fi

  if test "$CXXFLAGS" = ""; then
	CXXFLAGS=$ACX_OPT_FLAG
  fi
fi

])

dnl
dnl @synopsis ACX_HUTI([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl Look for hut iterative solver library. 
dnl
AC_DEFUN([ACX_HUTI], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])
AC_REQUIRE([ACX_LANG_COMPILER_MS])
acx_huti_ok=no

AC_ARG_WITH(huti,
	[AC_HELP_STRING([--with-huti=<lib>], [use HUTI library <lib>])])
case $with_huti in
	yes | "") ;;
	no) acx_huti_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) HUTI_LIBS="$with_huti" ;;
	*) HUTI_LIBS="-l$with_huti" ;;
esac

# Get fortran linker names of BLAS functions to check for.
AC_FC_FUNC(huti_d_gmres)
AC_FC_FUNC(huti_d_cgs)

acx_huti_save_LIBS="$LIBS"
dnl LIBS="$BLAS_LIBS $LAPACK_LIBS $LIBS $FCLIBS $FLIBS"

# First, check HUTI_LIBS environment variable
if test $acx_huti_ok = no; then
if test "x$HUTI_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$LIBS $HUTI_LIBS $BLAS_LIBS $FCLIBS $FLIBS"

	if test "$acx_cv_c_compiler_ms" = "yes"; then
		# windose shite	
		AC_MSG_CHECKING([for $huti_d_gmres[]@44 in $HUTI_LIBS])
		save_CFLAGS="$CFLAGS"
		CFLAGS="$CFLAGS -Gz"
		AC_LINK_IFELSE(
		[int main ()
		 {
		   HUTI_D_GMRES(1,2,3,4,5,6,7,8,9,10,11);
		   return 0;
		 }
		],
		[
	      		acx_huti_ok=yes
		],
		[
	 	        HUTI_LIBS=""	
		])
		AC_MSG_RESULT($acx_huti_ok)
		CFLAGS="$save_CFLAGS"
	else
		AC_MSG_CHECKING([for $huti_d_gmres in $HUTI_LIBS])
		AC_TRY_LINK_FUNC($huti_d_gmres, [acx_huti_ok=yes], [HUTI_LIBS=""])
		AC_MSG_RESULT($acx_huti_ok)
	fi
	LIBS="$save_LIBS"
fi
fi

# Generic HUTI library?
if test $acx_huti_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FCLIBS $FLIBS"
	if test "$acx_cv_c_compiler_ms" = "yes"; then
		# windose shite	
		AC_MSG_CHECKING([for $huti_d_gmres[]@44 in -lhuti])
		LIBS="-lhuti $BLAS_LIBS $FCLIBS $FLIBS"
		save_CFLAGS="$CFLAGS"
		CFLAGS="$CFLAGS -Gz"

		AC_LINK_IFELSE(
		[int main ()
		 {
		   HUTI_D_GMRES(1,2,3,4,5,6,7,8,9,10,11);
		   return 0;
		 }
		],
		[
	      		acx_huti_ok=yes
		],
		[
	 	        HUTI_LIBS=""	
		])
		AC_MSG_RESULT($acx_huti_ok)
		CFLAGS="$save_CFLAGS"
	else
		AC_CHECK_LIB(huti, $huti_d_gmres, [acx_huti_ok=yes; HUTI_LIBS="-lhuti"])
	fi
	LIBS="$save_LIBS"
fi

AC_SUBST(HUTI_LIBS)

LIBS="$acx_huti_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_huti_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_HUTI,1,[Define if you have a HUTI library.]),[$1])
        :
else
        acx_huti_ok=no
        $2
fi
])dnl ACX_HUTI


dnl
dnl @synopsis ACX_EIOF([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl Look for eio fortran library
dnl
AC_DEFUN([ACX_EIOF], [
AC_PREREQ(2.50)
AC_REQUIRE([ACX_LANG_COMPILER_MS])
acx_eiof_ok=no

AC_LANG_PUSH(C++)

AC_ARG_WITH(eiof,
	[AC_HELP_STRING([--with-eiof=<lib>], [use EIOF library <lib>])])
case $with_eiof in
	yes | "") ;;
	no) acx_eiof_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) EIOF_LIBS="$with_eiof" ;;
	*) EIOF_LIBS="-l$with_eiof" ;;
esac

# Get fortran linker names of EIO functions to check for.
AC_FC_FUNC(eio_init)

acx_eiof_save_LIBS="$LIBS"

LIBS="$LIBS"

# First, check EIO_LIBS environment variable
if test $acx_eiof_ok = no; then
if test "x$EIOF_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$EIOF_LIBS $LIBS"

	if test "$acx_cv_c_compiler_ms" = "yes"; then
		# windose shite	
		AC_MSG_CHECKING([for $eio_init[]@4 in $EIOF_LIBS])
		save_CFLAGS="$CFLAGS"
		CFLAGS="$CFLAGS -Gz"
		AC_LANG_PUSH(C)
		AC_LINK_IFELSE(
		[int main ()
		 {
		   $eio_init(1);
		   return 0;
		 }
		],
		[
	      		acx_eiof_ok=yes
		],
		[
	 	        EIOF_LIBS=""
		])
		AC_MSG_RESULT($acx_eiof_ok)
		CFLAGS="$save_CFLAGS"
		AC_LANG_POP(C)

	else
		AC_MSG_CHECKING([for $eio_init in $EIOF_LIBS])
		AC_TRY_LINK_FUNC($eio_init, [acx_eiof_ok=yes], [EIOF_LIBS=""])
		AC_MSG_RESULT($acx_eiof_ok)
	fi

	LIBS="$save_LIBS"
fi
fi

# Generic EIO library?
if test "$acx_eiof_ok" = no; then
	if test "$acx_cv_c_compiler_ms" = "yes"; then
		# windose shite	
		AC_MSG_CHECKING([for $eio_init[]@4 in -leiof])
		save_CFLAGS="$CFLAGS"
		CFLAGS="$CFLAGS -Gz"
		AC_LANG_PUSH(C)
		AC_LINK_IFELSE(
		[int main ()
		 {
		   $eio_init(1);
		   return 0;
		 }
		],
		[
	      		acx_eiof_ok=yes
			EIOF_LIBS="-leiof"
		],
		[
	 	        EIOF_LIBS=""
		])
		AC_MSG_RESULT($acx_eiof_ok)
		CFLAGS="$save_CFLAGS"
		AC_LANG_POP(C)
	else
		AC_CHECK_LIB(eiof, $eio_init, [acx_eiof_ok=yes; EIOF_LIBS="-leiof"])
	fi
fi

AC_SUBST(EIOF_LIBS)

LIBS="$acx_eiof_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_eiof_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_EIOF,1,[Define if you have a EIOF library.]),[$1])
        :
else
        acx_eiof_ok=no
        $2
fi
AC_LANG_POP(C++)
	])dnl ACX_EIO



dnl
dnl @synopsis ACX_EIOC([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl Look for eio library
dnl
AC_DEFUN([ACX_EIOC], 
[
AC_PREREQ(2.50)
acx_eioc_ok=no

AC_LANG_PUSH(C++)
AC_ARG_WITH(eioc,
	[AC_HELP_STRING([--with-eioc=<lib>], [use EIOC library <lib>])])
case $with_eioc in
	yes | "") ;;
	no) acx_eioc_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) EIOC_LIBS="$with_eioc" ;;
	*) EIOC_LIBS="-l$with_eioc" ;;
esac

acx_eioc_save_LIBS="$LIBS"

# First, check EIO_LIBS environment variable
if test $acx_eioc_ok = no; then
if test "x$EIOC_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$LIBS $EIOC_LIBS"
	AC_MSG_CHECKING([for eio_init in $EIOC_LIBS])
	AC_TRY_LINK_FUNC(eio_init, [acx_eioc_ok=yes], [EIOC_LIBS=""])
	AC_MSG_RESULT($acx_eioc_ok)
	LIBS="$save_LIBS"
fi
fi

# Generic EIO library?
if test "$acx_eioc_ok" = no; then
	AC_CHECK_LIB(eioc, eio_init, [acx_eioc_ok=yes; EIOC_LIBS="-leioc"])
fi

AC_SUBST(EIOC_LIBS)

LIBS="$acx_eioc_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_eioc_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_EIOC,1,[Define if you have a EIOC library.]),[$1])
        :
else
        acx_eioc_ok=no
        $2
fi
AC_LANG_POP(C++)
])dnl ACX_EIOC



AC_DEFUN([ACX_MEANING], [
AC_MSG_CHECKING([for answer to meaning of life])
AC_MSG_RESULT([42])
])

dnl
dnl @synopsis ACX_ARPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl Look for ARPACK library
dnl
AC_DEFUN([ACX_ARPACK], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
AC_REQUIRE([ACX_LANG_COMPILER_MS])
acx_arpack_ok=no

AC_ARG_WITH(arpack,
	[AC_HELP_STRING([--with-arpack=<lib>], [Specify location of ARPACK])])
case $with_arpack in
	yes | "") ;;
	no) acx_arpack_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) ARPACK_LIBS="$with_arpack" ;;
	*) ARPACK_LIBS="-l$with_arpack" ;;
esac

# Get fortran linker names of ARPACK functions to check for.
AC_FC_FUNC(dseupd)

acx_arpack_save_LIBS="$LIBS"

LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FCLIBS $FLIBS"

# First, check ARPACK_LIBS environment variable
if test $acx_arpack_ok = no; then
if test "x$ARPACK_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$ARPACK_LIBS $LIBS"
	AC_MSG_CHECKING([for $dseupd in $ARPACK_LIBS])

	if test "$acx_cv_c_compiler_ms" = "yes"; then
		# windose shite	
		save_CFLAGS="$CFLAGS"
		CFLAGS="$CFLAGS -Gz"
		AC_LINK_IFELSE(
		[int main ()
		 {
		   $dseupd(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25);
		   return 0;
		 }
		],
		[
	      		acx_arpack_ok=yes
		],
		[
	 	        ARPACK_LIBS=""	
		])
		AC_MSG_RESULT($acx_arpack_ok)
		CFLAGS="$save_CFLAGS"
	else
		AC_TRY_LINK_FUNC($dseupd, [acx_arpack_ok=yes], [ARPACK_LIBS=""])
	fi

	AC_MSG_RESULT($acx_arpack_ok)
	LIBS="$save_LIBS"
fi
fi

# Generic ARPACK library?
if test $acx_arpack_ok = no; then
	AC_CHECK_LIB(arpack, $dseupd, [acx_arpack_ok=yes; ARPACK_LIBS="-larpack"])
fi

AC_SUBST(ARPACK_LIBS)

LIBS="$acx_arpack_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_arpack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_ARPACK,1,[Define if you have a ARPACK library.]),[$1])
        :
else
        acx_arpack_ok=no
        $2
fi
])dnl ACX_ARPACK

dnl
dnl @synopsis ACX_PARPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl Look for PARPACK library
dnl
AC_DEFUN([ACX_PARPACK], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])
AC_REQUIRE([ACX_MPI])
AC_REQUIRE([ACX_ARPACK])
acx_parpack_ok=no

AC_ARG_WITH(parpack,
	[AC_HELP_STRING([--with-parpack=<lib>], [Specify location of PARPACK])])
case $with_parpack in
	yes | "") ;;
	no) acx_parpack_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) PARPACK_LIBS="$with_parpack" ;;
	*) PARPACK_LIBS="-l$with_parpack" ;;
esac

# Get fortran linker names of PARPACK functions to check for.
AC_FC_FUNC(pdneupd)

acx_parpack_save_LIBS="$LIBS"

LIBS="$MPI_LIBS $ARPACK_LIBS $LAPACK_LIBS $BLAS_LIBS $LIBS $FCLIBS $FLIBS"

# First, check PARPACK_LIBS environment variable
if test $acx_parpack_ok = no; then
if test "x$PARPACK_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$PARPACK_LIBS $LIBS"
	AC_MSG_CHECKING([for $pdneupd in $PARPACK_LIBS])
	AC_TRY_LINK_FUNC($pdneupd, [acx_parpack_ok=yes], [PARPACK_LIBS=""])
	AC_MSG_RESULT($acx_parpack_ok)
	LIBS="$save_LIBS"
fi
fi

# Generic PARPACK library?
if test $acx_parpack_ok = no; then
	AC_CHECK_LIB(parpack, $pdneupd, [acx_parpack_ok=yes; PARPACK_LIBS="-lparpack"])
fi

AC_SUBST(PARPACK_LIBS)

LIBS="$acx_parpack_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_parpack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_PARPACK,1,[Define if you have a PARPACK library.]),[$1])
        :
else
        acx_parpack_ok=no
        $2
fi
])dnl ACX_PARPACK



dnl
dnl @synopsis ACX_UMFPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl Look for UMFPACK library
dnl
AC_DEFUN([ACX_UMFPACK], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])
AC_REQUIRE([ACX_LANG_COMPILER_MS])
acx_umfpack_ok=no

AC_ARG_WITH(umfpack,
	[AC_HELP_STRING([--with-umfpack=<lib>], [Specify location of UMFPACK])])
case $with_umfpack in
	yes | "") ;;
	no) acx_umfpack_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) UMFPACK_LIBS="$with_umfpack" ;;
	*) UMFPACK_LIBS="-l$with_umfpack" ;;
esac

# Get fortran linker names of UMFPACK functions to check for.
AC_FC_FUNC(umf4def)

acx_umfpack_save_LIBS="$LIBS"

LIBS="$BLAS_LIBS $LAPACK_LIBS $LIBS $FCLIBS $FLIBS"

# First, check UMFPACK_LIBS environment variable
if test $acx_umfpack_ok = no; then
if test "x$UMFPACK_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$UMFPACK_LIBS $LIBS"
	AC_MSG_CHECKING([for $umf4def in $UMFPACK_LIBS])

	if test "$acx_cv_c_compiler_ms" = "yes"; then
		# windose shite	
		save_CFLAGS="$CFLAGS"
		CFLAGS="$CFLAGS -Gz"
		AC_LINK_IFELSE(
		[int main ()
		 {
		   $umf4def(1);
		   return 0;
		 }
		],
		[
	      		acx_umfpack_ok=yes
		],
		[
	 	        UMFPACK_LIBS=""	
		])
		CFLAGS="$save_CFLAGS"
	else
		AC_TRY_LINK_FUNC($umf4def, [acx_umfpack_ok=yes], [UMFPACK_LIBS=""])
	fi

	AC_MSG_RESULT($acx_umfpack_ok)
	LIBS="$save_LIBS"
fi
fi

# Generic UMFPACK library?
if test $acx_umfpack_ok = no; then
	AC_CHECK_LIB(umfpack, $umf4def, [acx_umfpack_ok=yes; UMFPACK_LIBS="-lumfpack -lamd"],,[-lamd])
fi

AC_SUBST(UMFPACK_LIBS)

LIBS="$acx_umfpack_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_umfpack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_UMFPACK,1,[Define if you have a UMFPACK library.]),[$1])
        :
else
        acx_umfpack_ok=no
        $2
fi
])dnl ACX_UMFPACK

dnl 
dnl look for the std c libraries:
dnl
dnl stdc++ (gnu)
dnl Cstd   (sun)
dnl C      (aix)
dnl 
dnl This is probably redundant. 
dnl
AC_DEFUN([ACX_CHECK_STDCXXLIB],
[
AC_REQUIRE([AC_PROG_CXX])
acx_check_stdcxxlib_save_LIBS=$LIBS
acx_stdcxxlib_ok=no

AC_ARG_WITH(stdcxxlib,
	[AC_HELP_STRING([--with-stdcxxlib=<lib>], [use std c++ library <lib>])])

case $with_stdcxxlib in
	yes | "") ;;
	no) acx_stdcxxlib_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) STDCXX_LIBS="$with_stdcxxlib" ;;
	*) STDCXX_LIBS="-l$with_stdcxxlib" ;;
esac

if test x$STDCXX_LIBS != x; then
	save_LIBS="$LIBS"; LIBS="$STDCXX_LIBS $LIBS"
	AC_MSG_CHECKING([for main in $STDCXX_LIBS])
	AC_TRY_LINK_FUNC(main, [acx_stdcxx_lib_ok=yes], [STDCXX_LIBS=""])
	AC_MSG_RESULT($acx_stdcxx_lib_ok)
	LIBS="$save_LIBS"
fi

dnl check for stdc++
if test $acx_stdcxxlib_ok = no; then
	AC_CHECK_LIB(stdc++, main,[
				   STDCXX_LIBS="-lstdc++"
			           acx_stdcxxlib_ok=yes
                                  ])			
fi

dnl check for stdc++
if test $acx_stdcxxlib_ok = no; then
	AC_CHECK_LIB(Cstd, main,[
				   STDCXX_LIBS="-lCstd"
			           acx_stdcxxlib_ok=yes
                                  ])
fi

dnl check for stdc++
dnl if test $acx_stdcxxlib_ok = no; then
dnl	AC_CHECK_LIB(C, main,[
dnl				   STDCXX_LIBS="-lC"
dnl			           acx_stdcxxlib_ok=yes
dnl                                  ])
dnl fi

if test $acx_stdcxxlib_ok = no; then
	AC_MSG_ERROR([Couldn't find std c++ library that is needed for linking.])
fi

LIBS=$acx_check_stdcxxlib_save_LIBS
])

dnl find out the flags that enable 64 bit compilation
AC_DEFUN([ACX_CHECK_B64FLAGS],
[
AC_PREREQ([2.50])
dnl AC_REQUIRE([AC_PROG_FC])
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([ACX_PROG_AR])
dnl AC_REQUIRE([AC_PROG_CXX])

AC_ARG_WITH(64bits,
	[AC_HELP_STRING([--with-64bits=yes(/no)], [Try to compile using 64 bits (default)])])


if test -z "$host"; then
  host=unknown
fi
canonical_host_type=$host
if test "$host" = unknown; then
  AC_MSG_ERROR([unknown system type, your build will most likely be screwed. quitting.])
fi

dnl by default, use no flags at all
B64FLAGS=
B64CFLAGS=
B64FCFLAGS=
B64FFLAGS=
B64CXXFLAGS=
orig_CFLAGS=$CFLAGS
orig_FFLAGS=$FFLAGS
orig_FCFLAGS=$FCFLAGS
orig_CXXFLAGS=$CXXFLAGS


AC_MSG_CHECKING([for 64 bit compilation flags])

if test "$with_64bits" = no; then
   AC_MSG_RESULT(not even going to try)
else
   AC_MSG_RESULT([let's see what happens])

case "$canonical_host_type" in
  *-*-386bsd* | *-*-openbsd* | *-*-netbsd*)
  ;;
  *-*-freebsd*)
  ;;
  alpha*-dec-osf*)
  ;;
  *-*-darwin*)
  ;;
  *-*-cygwin* | *-*-mingw*)
  ;;
  *-*-linux* | *-*-gnu*)
	B64FLAGS="-m64 -fPIC"
        dnl -M64
	if test x"$FC" != x; then
	case "$FC" in
  	  ifort | ifc)
		true
 		;;
	  g95 | gfortran)
	        B64FCFLAGS=$B64FLAGS
		;;
	  pgf*)
	        # portland group
	        B64FCFLAGS="-fPIC"
	        ;;
	  xl*)
	        B64FCFLAGS="-q64"
	        ;;
	  *)
	        B64FCFLAGS=$B64FLAGS
		;;
        esac
	fi

	if test x"$F77" != x; then
 	case "$F77" in 
	  ifort | ifc)
		true
	  ;;
	  pgf*)
	        # portland group
	        B64FFLAGS="-fPIC"
	        ;;
	  xl*)
	        B64FFLAGS="-q64"
	        ;;
          *)
	      B64FFLAGS=$B64FLAGS
          ;;
	esac
	fi

	if test x"$CC" != x; then
 	case "$CC" in 
	  icc | icc)
		true
	  ;;
	  pgcc*)
	        # portland group
	        B64CFLAGS="-fPIC"
	        ;;
	  xl*)
	        B64CFLAGS="-q64"
	        ;;
          *)
	      B64CFLAGS=$B64FLAGS
          ;;
	esac
	fi

	if test x"$CXX" != x; then
 	case "$CXX" in 
	  icc | icc)
		true
	  ;;
	  pgC*)
	        # portland group
	        B64CXXFLAGS="-fPIC"
	  ;;
	  xl*)
	        B64CXXFLAGS="-q64"
	        ;;
	  *)
       		B64CXXFLAGS=$B64FLAGS
	  ;;
	esac
	fi
  ;;
  rs6000-ibm-aix* | powerpc-ibm-aix*)
        B64FLAGS="-q64"

	if test "$acx_ar_userdefined" = no; then
   	  # aix has problems with 64 bits and ar, this fixes it
          AR='ar -Xany'
	  AC_SUBST(AR)
	  AC_MSG_CHECKING([for 64 bit ar])
	  AC_MSG_RESULT($AR)
	fi

	if test x"$FC" != x; then
	case "$FC" in 
          g*)
          	B64CFLAGS="-m64"
 	  ;;
	  *)
	  # hopefully ibm visual age
		B64FCFLAGS=$B64FLAGS
	  ;;
 	esac
	fi

	if test x"$F77" != x; then
	case "$F77" in
          g*)
          	B64FFLAGS="-m64"
 	  ;;
	  *)
		B64FFLAGS=$B64FLAGS
	  ;;
	esac
	fi
	
	if test x"$CC" != x; then
	case "$CC" in
          g*)
          	B64CFLAGS="-m64"
	  ;;
	  *)
	        B64CFLAGS=$B64FLAGS
	  ;;
 	esac
	fi
	
	if test x"$CXX" != x; then
	case "$CXX" in 
          g*)
          	B64CFLAGS="-m64"
 	  ;;
	  *)
		B64CXXFLAGS=$B64FLAGS
	  ;;
	esac
	fi
  ;;
  hppa*-hp-hpux*)
  ;;
  *-sgi-*)
  ;;
  sparc-sun-sunos4*)
  ;;
  sparc-sun-solaris2* | i386-pc-solaris2*)
        B64FLAGS="-xtarget=native64 -KPIC"
	SUN_64BIT_FLAGS=$B64FLAGS

	if test x"$FC" != x; then
	case "$FC" in 
          g*)
          	B64CFLAGS="-m64"
 	  ;;
	  mpf* | f*)
		B64FCFLAGS=$SUN_64BIT_FLAGS
	  ;;
 	esac
	fi

	if test x"$F77" != x; then
	case "$F77" in
          g*)
          	B64FFLAGS="-m64"
 	  ;;
	  mpf* | f*)
		B64FFLAGS=$SUN_64BIT_FLAGS
	  ;;
	esac
	fi
	
	if test x"$CC" != x; then
	case "$CC" in
          g*)
          	B64CFLAGS="-m64"
	  ;;
	  *)
	        B64CFLAGS=$SUN_64BIT_FLAGS
	  ;;
 	esac
	fi
	
	if test x"$CXX" != x; then
	case "$CXX" in 
          g*)
          	B64CFLAGS="-m64"
 	  ;;
	  mpCC | CC)
		B64CXXFLAGS=$SUN_64BIT_FLAGS
	  ;;
	esac
	fi
  ;;
esac

if test "$with_64bits" != no; then
	if test x"$CC" != x; then
           AC_MSG_CHECKING([for 64 bit CFLAGS])
           AC_MSG_RESULT($B64CFLAGS)
  	   CFLAGS="$CFLAGS $B64CFLAGS"
	fi

	if test x"$FC" != x; then
           AC_MSG_CHECKING([for 64 bit FCFLAGS])
           AC_MSG_RESULT($B64FCFLAGS)
	   FCFLAGS="$FCFLAGS $B64FCFLAGS"
	fi

	if test x"$CXX" != x; then
           AC_MSG_CHECKING([for 64 bit CXXFLAGS])
           AC_MSG_RESULT($B64CXXFLAGS	)
	   CXXFLAGS="$CXXFLAGS $B64CXXFLAGS"
	fi

	if test x"$F77" != x; then
          AC_MSG_CHECKING([for 64 bit FFLAGS])
          AC_MSG_RESULT($B64FFLAGS)
	  FFLAGS="$FFLAGS $B64FFLAGS"
	fi
fi
fi

dnl let's see if it works...
AC_CHECK_SIZEOF(void*)
case "$ac_cv_sizeof_voidp" in
  "8")
    AC_DEFINE(ARCH_64_BITS, 1,[64 bit arch.])

    if test x"$with_64bits" = xno; then
	AC_MSG_WARN([Explicitely requested 32 bits, but got 64 bits.])
    fi	
  ;;
  "4")
    AC_DEFINE(ARCH_32_BITS, 1,[32 bit arch.]) 
  ;;
  *)
    AC_DEFINE(ARCH_32_BITS, 1,[Couldn't determine. sticking with 32 bits.])
  ;;
esac

AM_CONDITIONAL(USE_64BIT_ARCH, test "$ac_cv_sizeof_voidp" -eq "8")

if test "$with_64bits" != no; then
   AC_MSG_CHECKING(to see if we got 64 bits)

   if test "$ac_cv_sizeof_voidp" -ne 8; then
      AC_MSG_RESULT([nope, reverting compiler flags]) 

      dnl FIXME: test that all compilers are 64 bit
      B64FLAGS=""
      CFLAGS=$orig_CFLAGS
      FFLAGS=$orig_FFLAGS
      FCFLAGS=$orig_FCFLAGS
      CXXFLAGS=$orig_CXXFLAGS
   else
      AC_MSG_RESULT([oh yes]) 
   fi
fi

AC_SUBST(B64FLAGS)
])dnl ACX_CHECK_B64FLAGS

dnl
dnl @synopsis ACX_MATC([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl Look for matc library
dnl
AC_DEFUN([ACX_MATC], [
AC_PREREQ(2.50)
AC_REQUIRE([ACX_LANG_COMPILER_MS])
acx_matc_ok=no

AC_ARG_WITH(matc,
	[AC_HELP_STRING([--with-matc=<lib>], [use MATC library <lib>])])
case $with_matc in
	yes | "") ;;
	no) acx_matc_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) MATC_LIBS="$with_matc" ;;
	*) MATC_LIBS="-l$with_matc" ;;
esac

acx_matc_save_LIBS="$LIBS"

# First, check EIO_LIBS environment variable
if test $acx_matc_ok = no; then
if test "x$MATC_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$MATC_LIBS $LIBS"
	AC_MSG_CHECKING([for mtc_init in $MATC_LIBS])

	AC_TRY_LINK_FUNC(mtc_init, [acx_matc_ok=yes], [MATC_LIBS=""])

	AC_MSG_RESULT($acx_matc_ok)
	LIBS="$save_LIBS"
fi
fi

# Generic MATC library?
if test "$acx_matc_ok" = no; then
	AC_CHECK_LIB(matc, mtc_init, [acx_matc_ok=yes; MATC_LIBS="-lmatc"],,[-lm])
fi

AC_SUBST(MATC_LIBS)

LIBS="$acx_matc_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_matc_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_MATC,1,[Define if you have a MATC library.]),[$1])
        :
else
        acx_matc_ok=no
        $2
fi
])dnl ACX_MATC

dnl
dnl We really need the old style cpp for preprocessing fortran.
dnl 
AC_DEFUN([ACX_PROG_TRADITIONAL_CPP], [
# sun cc -E leaves nasty # comment that chokes the fortran compiler, so we have to hope
# that ye olde cpp is present.
AC_CHECK_FILE([/lib/cpp],[TRADITIONAL_CPP=yes],[TRADITIONAL_CPP=no])

if test "$TRADITIONAL_CPP" != "yes"; then
   AC_CHECK_PROG(BASIC_CPP,[cpp], yes, no)
   if test "$BASIC_CPP" = yes; then
       CPP=cpp
   else 
       AC_MSG_ERROR([Traditional cpp not found, just have to exit for for now.])	
   fi
else
   CPP="/lib/cpp"
fi

case "$canonical_host_type" in
  *-*-cygwin* | *-*-mingw*)
	TRADITIONAL_CPP_FLAGS="-traditional-cpp"
  ;;
  *-*-linux* | *-*-gnu*)
	TRADITIONAL_CPP_FLAGS="-traditional-cpp"
  ;;
  *)
	TRADITIONAL_CPP_FLAGS=""
  ;;
esac

AC_SUBST(CPP)
])

AC_DEFUN([ACX_FC_ETIME],[
AC_MSG_CHECKING([for fortran intrinsic etime])

AC_LANG_PUSH(Fortran)
AC_LINK_IFELSE(
[    
      PROGRAM TEST                                                                        
      INTRINSIC ETIME                                                                     
      REAL ETIME, T1, TARRAY(2)                                                           
      T1=ETIME(TARRAY)                                                                    
      END     
],
[
     AC_MSG_RESULT([found])
     AC_DEFINE(HAVE_F_ETIME,1,[Does the fortran environment implement etime])
],
[
     AC_MSG_RESULT([missing])
])
AC_LANG_POP(Fortran)
])

AC_DEFUN([ACX_FC_FLUSH],[
AC_MSG_CHECKING([for fortran intrinsic flush])

AC_LANG_PUSH(Fortran)

AC_LINK_IFELSE(
[    
      PROGRAM TEST                                                                        
      CALL FLUSH(1)                                                                    
      END     
],
[
     AC_MSG_RESULT([found])
     AC_DEFINE(HAVE_F_FLUSH,1,[Does the fortran environment implement flush])
],
[
     AC_MSG_RESULT([missing])
])
AC_LANG_POP(Fortran)
])dnl ACX_FC_FLUSH

AC_DEFUN([ACX_BOURNE_SHELL],[

case "$canonical_host_type" in
  *solaris* | *sun*)	
	dnl Sun harbours a defunct sh?
	BOURNE_SHELL="bash"
  ;;
  *) 
	BOURNE_SHELL="sh"
  ;;
esac
AC_SUBST(BOURNE_SHELL)
])

dnl
dnl Determine if fortran compiler expects (char *str1, int len1) or (char *str1)
dnl
AC_DEFUN([ACX_FC_CHAR_MANGLING], [
AC_PREREQ(2.50)
AC_REQUIRE([ACX_LANG_COMPILER_MS])

acx_cv_fc_char_mangling="char_ptr"

AC_FC_FUNC(barf)



AC_MSG_CHECKING([for Fortran char* mangling scheme])

# this test is hopelessly fucked with a windows compiler (on windows t
# if test "$acx_cv_c_compiler_ms" = "yes"; then 
#	acx_cv_fc_char_mangling="char_ptr_and_len_int"
# else


AC_LANG_PUSH(C)
AC_COMPILE_IFELSE([
void $barf(char *name, int *l2, int *l3)
{
   printf("%d",*l2);
 }],
 [
    mv conftest.$ac_objext cfortran_test.$ac_objext
 ],
 [
    ac_cv_[]_AC_LANG_ABBREV[]_char_mangling="unknown"
])
AC_LANG_POP(C)

AC_LANG_PUSH(Fortran)
AC_COMPILE_IFELSE(
[      program testingit                                                                   
       character s                                                                         
       parameter (s='my ass')                                                              
       call barf(s,42)                                                                     
       end],
[
$FC $FCFLAGS -o testi$ac_exeext conftest.$ac_objext cfortran_test.$ac_objext
the_answer=`./testi$ac_exeext`

case $the_answer in 
	42)  
	       acx_cv_fc_char_mangling="char_ptr"
        ;;
	*)
   	       acx_cv_fc_char_mangling="char_ptr_and_len_int"
	;;
esac

# fi windows 
# fi


AC_MSG_RESULT($acx_cv_fc_char_mangling)

rm testi$ac_exeext conftest.$ac_objext cfortran_test.$ac_objext],
[AC_MSG_RESULT([test failed, assuming char_ptr])])


AH_TEMPLATE(_AC_FC[_CHAR_PTR],[Char pointer mangling])

case $acx_cv_fc_char_mangling in
	char_ptr)
		AC_DEFINE(_AC_FC[_CHAR_PTR(P,L)],[char *P])
	;;
	*)
		AC_DEFINE(_AC_FC[_CHAR_PTR(P,L)],[char *P, int L])
	;;
esac
AC_LANG_POP(Fortran)
])

AC_DEFUN([ACX_FC_LINKTYP], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_FC_WRAPPERS])
dnl setup name mangling defines
case "$ac_cv_fc_mangling" in
  "upper case, no underscore, no extra underscore" | "upper case, no underscores")
    AC_DEFINE(ELMER_LINKTYP, 2, "Mangling: upper case, no underscore, no extra underscore") 
  ;;
  "lower case, no underscore, no extra underscore" | "lower case, no underscores" )
    AC_DEFINE(ELMER_LINKTYP, 3, "Mangling: lower case, no underscore, no extra underscore") 
  ;;
  "lower case, underscore, no extra underscore" | "lower case, single underscores")
    AC_DEFINE(ELMER_LINKTYP, 1, "Mangling: lower case, underscore, no extra underscore") 
  ;;
  "lower case, underscore, extra underscore" | "lower case, double underscores")
    AC_DEFINE(ELMER_LINKTYP, 4, "Mangling: lower case, underscore, extra underscore") 
  ;;
  *)
    AC_MSG_ERROR([Unknown sort of Fortran name mangling]) 
  ;;
esac

])

AC_DEFUN([ACX_FC_INCLUDE_MODULE_FLAG],[

INCLUDE_MODULE_FLAG="-I"

case "$canonical_host_type" in
  sparc-sun-solaris2* | i386-pc-solaris2*)
    if test "$ac_cv_fc_compiler_gnu" != yes; then
      INCLUDE_MODULE_FLAG="-M"
    fi
  ;;
  *darwin*)
    #absoft
    if test "$FC" = "f90"; then
      INCLUDE_MODULE_FLAG="-p"
    fi
  ;;
esac

])

AC_DEFUN([ACX_PIC_FLAGS],[
dnl 
dnl Host specific stuff
dnl fixme: if 64 bits, then start looking for pic flags, not by default?
dnl 
dnl
case "$canonical_host_type" in
  *-*-386bsd* | *-*-openbsd* | *-*-netbsd*)
	true
  ;;
  *-*-freebsd*)
	true
  ;;
  alpha*-dec-osf*)
	true
  ;;
  *-*-darwin*)
	true
  ;;
  *-*-cygwin* | *-*-mingw*)
	true
  ;;
  *-*-linux* | *-*-gnu*)
    CPICFLAG="-fPIC"
    CXXPICFLAG="-fPIC"
    FPICFLAG="-fPIC"
    FCPICFLAG="-fPIC"

	if test x"$CC" != x; then
	case "$CC" in 
	  xl*)
	    CPICFLAG=""
   	    CXXPICFLAG=""
	    FPICFLAG=""
	    FCPICFLAG=""
	  ;;
	esac
	fi
  ;;
  i[[3456]]86-*-sco3.2v5*)
	true
  ;;
  rs6000-ibm-aix* | powerpc-ibm-aix*)
	true
  ;;
  hppa*-hp-hpux*)
    if test "$ac_cv_f77_compiler_gnu" = yes; then
      FPICFLAG=-fPIC
    else
      FPICFLAG=+Z
    fi
  ;;
  *-sgi-*)
	true
  ;;
  sparc-sun-sunos4*)
    if test "$ac_cv_f77_compiler_gnu" = yes; then
      FPICFLAG=-fPIC
    else
      FPICFLAG=-PIC
    fi
  ;;
  sparc-sun-solaris2* | i386-pc-solaris2*)
    if test "$ac_cv_f77_compiler_gnu" = yes; then
      FPICFLAG=-fPIC
    else
      FPICFLAG=-KPIC
    fi
    if test "$ac_cv_fc_compiler_gnu" = yes; then
      FPICFLAG=-fPIC
    else
      FPICFLAG=-KPIC
    fi
    if test "$GCC" = yes; then
      CPICFLAG="-fPIC"
    else
      CPICFLAG=-KPIC
    fi
    if test "$GXX" = yes; then
      CXXPICFLAG=-fPIC
    else
      CXXPICFLAG=-KPIC
    fi
  ;;
esac
FCFLAGS="$FCFLAGS $FPICFLAG"
FFLAGS="$FFLAGS $FPICFLAG"
CFLAGS="$CFLAGS $CPICFLAG"
CXXFLAGS="$CXXFLAGS $CXXPICFLAG"
])

AC_DEFUN([ACX_SHLIB_STUFF],[
dnl  
dnl OS-specific stuff that hasn't been (or can't be) done in a clean way.
dnl
DLFCN_DIR=
CPICFLAG="-fPIC"
CXXPICFLAG="-fPIC"
FPICFLAG="-fPIC"
SHLEXT="so"
SHLLIB='$(SHLEXT)'
SHLBIN=
SHLEXT_VER='$(SHLEXT).$(version)'
SHLLIB_VER='$(SHLLIB).$(version)'
SHLBIN_VER='$(SHLBIN).$(version)'
SHLLINKEXT=
LIBEXT=a
SH_LD2=$CXX
SH_LD=$FC
SH_LDFLAGS=-shared
SH_LINKING_TO_FLAGS=
DL_LD='$(SH_LD)'
DL_LDFLAGS='$(SH_LDFLAGS)'
SONAME_FLAGS=
LD_LIBRARY_PATH_VAR=LD_LIBRARY_PATH
LIBSOLVER_DEPS=$LIBS
RPATH_FLAG=
SH_EXPALL_FLAG=
dnl 
dnl Host specific stuff
dnl
case "$canonical_host_type" in
  *-*-386bsd* | *-*-openbsd* | *-*-netbsd*)
    SH_LD="ld"
    SH_LDFLAGS="-Bshareable"
  ;;
  *-*-freebsd*)
    SH_LD='$(CC)'
    SH_LDFLAGS="-shared"
  ;;
  alpha*-dec-osf*)
    SH_LDFLAGS="-shared"
  ;;
  *-*-darwin*)
    SH_LD="gcc"
    SH_LDFLAGS='-dynamiclib -undefined dynamic_lookup -single_module ${LDFLAGS}'
    SHLEXT="dylib"
    LD_LIBRARY_PATH_VAR=DYLD_LIBRARY_PATH	
  ;;
  *-*-cygwin* | *-*-mingw*)
       SHLEXT=dll
       SH_LDFLAGS="-shared"
  ;;
  *-*-linux* | *-*-gnu*)
	RPATH_FLAG="-Wl,-rpath "
	SH_EXPALL_FLAG="-Wl,--export-dynamic"
  ;;
  i[[3456]]86-*-sco3.2v5*)
    SH_LDFLAGS="-G"
  ;;
  rs6000-ibm-aix* | powerpc-ibm-aix*)
    SH_LDFLAGS="-G $ACX_LOPT_FLAGS"
    SH_LINKING_TO_FLAGS="-brtl -bexpall -bshared"
    LD_LIBRARY_PATH_VAR=LIBPATH
#    RPATH_FLAG="-blibpath:"
    SH_EXPALL_FLAG="-bexpall"
  ;;
  hppa*-hp-hpux*)
    SH_LDFLAGS="-shared -fPIC"
  ;;
  *-sgi-*)
      true
  ;;
  sparc-sun-sunos4*)
    SH_LD=ld
    SH_LDFLAGS="-assert nodefinitions"
    if test "$GXX" != yes; then
      SH_LDFLAGS=-G
      RPATH_FLAG="-R"
    fi


  ;;
  sparc-sun-solaris2* | i386-pc-solaris2*)
    if test "$GXX" != yes; then
      SH_LDFLAGS=-G
      RPATH_FLAG="-R"
    fi
  ;;
esac

case $FC in 
    pgf*) 
	# portland group does it differently
	RPATH_FLAG="-R"
    ;;
esac

AC_SUBST(LD_LIBRARY_PATH_VAR)
AC_SUBST(RPATH_FLAG)

### Dynamic linking is now enabled only if we are building shared
### libs and some API for dynamic linking is detected.

LD_CXX='$(CXX)'
LIBDLFCN=
DLFCN_INCFLAGS=
RDYNAMIC_FLAG=
DL_API_MSG=""
dlopen_api=false
shl_load_api=false
loadlibrary_api=false
dyld_api=false

### yes atleast now
SHARED_LIBS=true
if $SHARED_LIBS || $ENABLE_DYNAMIC_LINKING; then

  ### Check for dyld first since OS X can have a non-standard libdl	

    AC_CHECK_LIB(dld, shl_load)
    AC_CHECK_FUNCS(shl_load shl_findsym)
    if test "$ac_cv_func_shl_load" = yes \
      && test "$ac_cv_func_shl_findsym" = yes; then
      shl_load_api=true
    else
	if test "$acx_cv_c_compiler_ms" = "yes"; then
		# in windows, no generic test seems to work 
		AC_MSG_CHECKING([for LoadLibrary in windows])
		AC_LINK_IFELSE(
		[#include <stdio.h>
		 #include <windows.h>
		 int main ()
		 {
			LoadLibrary("kernel32.dll");
			return 0;
		 }
		],
		[
			ac_cv_func_loadlibrary=yes
		],
		[
			ac_cv_func_loadlibrary=no
		])
		AC_MSG_RESULT($ac_cv_func_loadlibrary)
      fi
      if test "$ac_cv_func_loadlibrary" = yes; then
        loadlibrary_api=true
      else
        AC_CHECK_LIB(dl, dlopen)
        AC_CHECK_FUNCS(dlopen dlsym dlerror dlclose)
        if test "$ac_cv_func_dlclose" = yes \
          && test "$ac_cv_func_dlerror" = yes \
          && test "$ac_cv_func_dlopen" = yes \
          && test "$ac_cv_func_dlsym" = yes; then
          dlopen_api=true
        else
	  case "$canonical_host_type" in
	    rs6000-ibm-aix* | powerpc-ibm-aix*)
	      LIBDLFCN="-ldlfcn -ll -lld"
	      DLFCN_INCFLAGS='-I$(top_srcdir)/dlfcn -I$(TOPDIR)/dlfcn'
	      dlopen_api=true
	    ;;
	    i[[3456]]86-*-sco3.2v5*)
	      LD_CXX='LD_RUN_PATH=$LD_RUN_PATH:$(octlibdir) $(CXX)'
	      dlopen_api=true
	    ;;
	  esac
	fi
      fi
    fi

  if $dlopen_api; then
    DL_API_MSG="(dlopen)"
    AC_DEFINE(HAVE_DLOPEN_API, 1, [Define if your system has dlopen, dlsym, dlerror, and dlclose for dynamic linking])
  elif $shl_load_api; then
    DL_API_MSG="(shl_load)"
    AC_DEFINE(HAVE_SHL_LOAD_API, 1, [Define if your system has shl_load and shl_findsym for dynamic linking])
  elif $loadlibrary_api; then
    DL_API_MSG="(LoadLibrary)"
    AC_DEFINE(HAVE_LOADLIBRARY_API, 1, [Define if your system has LoadLibrary for dynamic linking])
  elif $dyld_api; then
    DL_API_MSG="(dyld)"
    AC_DEFINE(HAVE_DYLD_API, 1, [Define if your system has dyld for dynamic linking])
  fi

  if $dlopen_api || $shl_load_api || $loadlibrary_api || $dyld_api; then
    ENABLE_DYNAMIC_LINKING=true
    AC_DEFINE(ENABLE_DYNAMIC_LINKING, 1, [Define if using dynamic linking])
  fi
  AC_DEFINE_UNQUOTED(SHL_EXTENSION, ".$SHLEXT",[Shared lib filename extension])
fi

AC_SUBST(SH_EXPALL_FLAG)

])



AC_DEFUN([ACX_PLATFORM_DEFS],
[
AC_REQUIRE([ACX_HOST])
acx_platform_def="GENERIC"
case "$canonical_host_type" in
  *-*-386bsd* | *-*-openbsd* | *-*-netbsd*)
	acx_platform_def="BSD"
        AC_DEFINE([BSD],1,[Detected platform.])
  ;;
  *-*-freebsd*)
	acx_platform_def="BSD"
        AC_DEFINE([BSD],1,[Detected platform.])
  ;;
  alpha*-dec-osf*)
	acx_platform_def="DEC_ALPHA"
        AC_DEFINE([DEC_ALPHA],1,[Detected platform.])
  ;;
  *-*-darwin*)
        AC_DEFINE([DARWIN],1,[Detected platform.])
  ;;
  *-*-cygwin*)
	acx_platform_def="WIN32"
        AC_DEFINE([CYGWIN],1,[Detected platform.])
  ;;
  *-*-mingw*)
	acx_platform_def="WIN32"
        AC_DEFINE([MINGW32],1,[Detected platform.])
        AC_DEFINE([WIN32],1,[Detected platform2.])
  ;;
  *-*-linux* | *-*-gnu*)
        AC_DEFINE([LINUX],1,[Detected platform.])
  ;;
  i[[3456]]86-*-sco3.2v5*)
        AC_DEFINE([BASTARDS],1,[Detected platform.])
  ;;
  rs6000-ibm-aix* | powerpc-ibm-aix*)
	acx_platform_def="AIX"
        AC_DEFINE([AIX],1,[Detected platform.])
  ;;
  hppa*-hp-hpux*)
        AC_DEFINE([HPUX],1,[Detected platform.])
  ;;
  *-sgi-*)
        AC_DEFINE([SGI],1,[Detected platform.])
  ;;
  sparc-sun-sunos4*)
        AC_DEFINE([SUNOS],1,[Detected platform.])
  ;;
  sparc-sun-solaris2* | i386-pc-solaris2*)
        AC_DEFINE([SOLARIS],1,[Detected platform.])
  ;;
esac

AM_CONDITIONAL(IBM_IS_PIECE_OF_SHIT, test "$acx_platform_def" = "AIX")
])

AC_DEFUN([ACX_COMPILER_FIXES],
[
AC_REQUIRE([ACX_SHLIB_STUFF])
case "$FC" in
   ifc)
	if test "$ac_cv_cxx_compiler_gnu" = yes; then
		dnl remove intel c++ stuff
		FCLIBS=`echo $FCLIBS | sed -e 's/-lintrins//g'`
	fi
	FC="$FC -Vaxlib"
   ;;
esac

case "$F77" in
   ifc)
	if test "$ac_cv_cxx_compiler_gnu" = yes; then
		dnl remove intel c++ stuff
		FLIBS=`echo $FLIBS | sed -e 's/-lintrins//g'`
	fi
   ;;
esac

case "$SH_LD" in
   ifc)
	SH_LD="$SH_LD -Vaxlib"
   ;;
esac

case "$SH_LD2" in
   ifc)
	SH_LD2="$SH_LD2 -Vaxlib"
   ;;
esac

AC_SUBST(SH_LDFLAGS)
AC_SUBST(SH_LD)
])


dnl
dnl @synopsis ACX_TCLTK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl Look for tcl/tk libraries 
dnl
AC_DEFUN([ACX_TCLTK], 
[
AC_REQUIRE([AC_PATH_X])
AC_PREREQ(2.50)
acx_tcltk_ok=no

AC_ARG_WITH(tcltk,
	[AC_HELP_STRING([--with-tcltk=<lib>], [Specify tcl & tk libraries])])
case $with_tcltk in
	yes | "") ;;
	no) acx_tcltk_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) TCLTK_LIBS="$with_tcltk" ;;
	*) TCLTK_LIBS="-l$with_tcltk" ;;
esac

acx_tcltk_save_LIBS="$LIBS"

# First, check TCLTK_LIBS environment variable
if test $acx_tcltk_ok = no; then
if test "x$TCLTK_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$TCLTK_LIBS $LIBS"
	AC_MSG_CHECKING([for TkGetDisplay in $TCLTK_LIBS])
	AC_TRY_LINK_FUNC(TkGetDisplay, [acx_tk_ok=yes], [TCLTK_LIBS=""])

	AC_MSG_CHECKING([for TclInvoke in $TCLTK_LIBS])
	AC_TRY_LINK_FUNC(TclInvoke, [acx_tcl_ok=yes], [TCLTK_LIBS=""])

	if test "$acx_tk_ok" = yes; then
	if test "$acx_tcl_ok" = yes; then
		acx_tcltk_ok=yes
	fi
	fi	

	AC_MSG_RESULT($acx_tcltk_ok)
	LIBS="$save_LIBS"
fi
fi

acx_tcltk_lib_versions="8.4 8.3 8.2 8.1 84 83 82 81"
acx_tcltk_locations="/usr/lib /usr/local/lib /usr/swf/lib"

# Generic TCLTK library?
if test "$acx_tcltk_ok" = no; then
   for l in $acx_tcltk_locations; do
	for v in $acx_tcltk_lib_versions; do
		acx_tcl_ok="no"
		acx_tk_ok="no"
		
		acx_tcltk_LDFLAGS_save=$LDFLAGS
		acx_tcltk_LIBS_save=$LIBS
		LDFLAGS="$LDFLAGS -L$l"
		LIBS="-ltcl$v -ltk$v"
		AC_MSG_CHECKING([for -ltcl$v -ltk$v libs in -L$l])
		
		AC_TRY_LINK_FUNC(TclInvoke, [acx_tcl_ok=yes; TCL_LIBS="-ltcl$v"])
		AC_TRY_LINK_FUNC(TkGetDisplay, [acx_tk_ok=yes; TK_LIBS="-ltk$v"])

		LDFLAGS=$acx_tcltk_LDFLAGS_save
		LIBS=$acx_tcltk_LIBS_save
		
		if test "$acx_tcl_ok" = yes; then
  		   if test "$acx_tk_ok" = yes; then
		     acx_tcltk_ok="yes"
  		     TCLTK_LIBS="-L$l $TK_LIBS $TCL_LIBS"
  		     AC_MSG_RESULT($acx_tcltk_ok)
		     break
		   fi
		fi
		AC_MSG_RESULT($acx_tcltk_ok)
	done
	# break again
	if test "$acx_tcl_ok" = yes; then
 	   if test "$acx_tk_ok" = yes; then
 	     break
 	   fi
	fi
   done
fi

AC_SUBST(TCLTK_LIBS)
LIBS=$acx_tcltk_save_LIBS

# Search for tcl.h and tk.h
if test "$TCLTK_INCPATH"; then
	acx_tcltk_tcl_h_locs=$TCLTK_INCPATH
fi
acx_tcltk_tcl_h_locs="$acx_tcltk_tcl_h_locs /usr/include /usr/local/include /usr/include/tcl8.4 /usr/include/tcl8.3 /usr/include/tcl8.2 /include /usr/swf/include /sw/include /sw/usr/include /sw/usr/include/tcl8.4 /really/weird/place /ok/I/quit"


acx_tcltk_CPPFLAGS_save=$CPPFLAGS
acx_tcltk_CFLAGS_save=$CFLAGS

AC_LANG_PUSH(C)

if test "$acx_tcltk_ok" = yes; then
for v in $acx_tcltk_tcl_h_locs; do
	acx_tcl_h_ok="no"
	acx_tk_h_ok="no"

	CPPFLAGS="-I$v $CPPFLAGS -I${x_includes}"

	AC_MSG_CHECKING([for tcl.h in -I$v])
	AC_PREPROC_IFELSE(
	   [AC_LANG_PROGRAM([#include <tcl.h>])],
	   [
		AC_MSG_RESULT([ok])
		AC_MSG_CHECKING([for tk.h in -I$v])
        	AC_PREPROC_IFELSE(
		   [AC_LANG_PROGRAM([#include <tk.h>])],
		   [
			AC_MSG_RESULT([ok])
	   	        TCLTK_INCLUDE="-I$v"
			acx_tcltk_h_ok=yes
			break
		   ],
	  	   [AC_MSG_RESULT([no])])
	   ],
	   [AC_MSG_RESULT([no])])
done
fi
AC_LANG_POP(C)

CPPFLAGS=$acx_tcltk_CPPFLAGS_save
CFLAGS=$acx_tcltk_CFLAGS_save

if test "$acx_tcltk_h_ok" != yes; then
	AC_MSG_WARN([Couldn't determine tcl.h and tk.h location. Specify it manually with CFLAGS and CXXFLAGS])
fi


# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_tcltk_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_MATC,1,[Define if you have a MATC library.]),[$1])
        :
else
        acx_tcltk_ok=no
        $2
fi
])dnl ACX_TCLTK




dnl
dnl see how well fortran cpp does
dnl explicitely set acx_fortran_cpp_ok=yes, if we know that the fortran compiler can handle 
dnl our C/Fortran preprocessing intructions
dnl
AC_DEFUN([ACX_FORTRAN_CPP], 
[
AC_PREREQ(2.50)
acx_fortran_cpp_ok=no

case "$FC" in
	*g95*)
		FORTRAN_CPP_FLAG="-cpp"
		acx_fortran_cpp_ok=yes
	;;
esac
AC_SUBST(FORTRAN_CPP_FLAG)
])


AC_DEFUN([ACX_LANG_COMPILER_MS],
[AC_CACHE_CHECK([whether we are using the Microsoft _AC_LANG compiler],
                [acx_cv_[]_AC_LANG_ABBREV[]_compiler_ms],
[AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[#ifndef _MSC_VER
       choke me
#endif
]])],
                   [acx_compiler_ms=yes],
                   [acx_compiler_ms=no])
acx_cv_[]_AC_LANG_ABBREV[]_compiler_ms=$acx_compiler_ms
])

if test "$acx_cv_c_compiler_ms" = "yes"; then
	AC_DEFINE(STDCALLBULL,[__stdcall],[Standard windows call declaration])
	AC_DEFINE(C_DLLEXPORT,[__declspec(dllexport)],[Standard windows call declaration])

else
	AC_DEFINE(STDCALLBULL,,[Standard windows call declaration])
	AC_DEFINE(C_DLLEXPORT,,[Standard windows call declaration])
fi
AM_CONDITIONAL(USE_WINDOWS_COMPILER, test "$acx_cv_c_compiler_ms" = "yes")
])
