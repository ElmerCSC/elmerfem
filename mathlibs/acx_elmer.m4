dnl 
dnl Elmer specific M4sh macros 
dnl
dnl @version $Id: acx_elmer.m4,v 1.40 2005/05/26 08:25:41 vierinen Exp $
dnl @author juha.vierinen@csc.fi 5/2005
dnl

dnl
dnl define host variable
dnl
AC_DEFUN([ACX_HOST],
[
if test -z "$host"; then
  host=unknown
fi
canonical_host_type=$host
if test "$host" = unknown; then
  AC_MSG_ERROR([configuring for unknown system type, your build will most likely be screwed.])
fi
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
	AC_MSG_CHECKING([for $huti_d_gmres in $HUTI_LIBS])
	AC_TRY_LINK_FUNC($huti_d_gmres, [acx_huti_ok=yes], [HUTI_LIBS=""])
	AC_MSG_RESULT($acx_huti_ok)
	LIBS="$save_LIBS"
fi
fi

# Generic HUTI library?
if test $acx_huti_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FCLIBS $FLIBS"
	AC_CHECK_LIB(huti, $huti_d_gmres, [acx_huti_ok=yes; HUTI_LIBS="-lhuti"])
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
	AC_MSG_CHECKING([for $eio_init in $EIOF_LIBS])
	AC_TRY_LINK_FUNC($eio_init, [acx_eiof_ok=yes], [EIOF_LIBS=""])
	AC_MSG_RESULT($acx_eiof_ok)
	LIBS="$save_LIBS"
fi
fi

# Generic EIO library?
if test "$acx_eiof_ok" = no; then
	AC_CHECK_LIB(eiof, $eio_init, [acx_eiof_ok=yes; EIOF_LIBS="-leiof"])
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

LIBS="-leioc $LIBS"

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


dnl
dnl @synopsis ACX_ARPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl Look for ARPACK library
dnl
AC_DEFUN([ACX_ARPACK], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])
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

LIBS="$BLAS_LIBS $LAPACK_LIBS $LIBS $FCLIBS $FLIBS"

# First, check ARPACK_LIBS environment variable
if test $acx_arpack_ok = no; then
if test "x$ARPACK_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$ARPACK_LIBS $LIBS"
	AC_MSG_CHECKING([for $dseupd in $ARPACK_LIBS])
	AC_TRY_LINK_FUNC($dseupd, [acx_arpack_ok=yes], [ARPACK_LIBS=""])
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

LIBS="$BLAS_LIBS $LAPACK_LIBS $LIBS $FCLIBS $FLIBS"

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
	AC_TRY_LINK_FUNC($umf4def, [acx_umfpack_ok=yes], [UMFPACK_LIBS=""])
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
if test $acx_stdcxxlib_ok = no; then
	AC_CHECK_LIB(C, main,[
				   STDCXX_LIBS="-lC"
			           acx_stdcxxlib_ok=yes
                                  ])
fi

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

acx_cv_fc_char_mangling="char_ptr"

AC_FC_FUNC(barf)

AC_MSG_CHECKING([for Fortran char* mangling scheme])

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
MKOCTFILE_DL_LDFLAGS='$(DL_LDFLAGS)'
SONAME_FLAGS=
library_path_var=LD_LIBRARY_PATH
LIBSOLVER_DEPS=$LIBS

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
    SH_LDFLAGS='-dynamiclib -single_module $(LDFLAGS)'
dnl    SHLEXT="dylib"
dnl    SHLLIB='$(SHLEXT)'
dnl    SHLEXT_VER='$(version).$(SHLEXT)'
dnl    SHLLIB_VER='$(version).$(SHLLIB)'
dnl    NO_OCT_FILE_STRIP="true"
dnl    SONAME_FLAGS='-install_name $(octlibdir)/$@'
dnl    library_path_var=DYLD_LIBRARY_PATH	
  ;;
  *-*-cygwin* | *-*-mingw*)
dnl    SHLEXT=dll
dnl    SHLLIB=dll.a
dnl    SHLBIN=dll
    SH_LDFLAGS="-shared"
dnl    SHLLINKEXT=".dll"
dnl    SONAME_FLAGS='-Wl,--out-implib=$@.a'
dnl    library_path_var=PATH
  ;;
  *-*-linux* | *-*-gnu*)
dnl    MKOCTFILE_DL_LDFLAGS="-shared -Wl,-Bsymbolic"
dnl    SONAME_FLAGS=''
dnl    RLD_FLAG='-Wl,-rpath -Wl,$(octlibdir)'
dnl     CPICFLAG="-fPIC"
dnl    CXXPICFLAG="-fPIC"
dnl    FPICFLAG="-fPIC"
  ;;
  i[[3456]]86-*-sco3.2v5*)
dnl    SONAME_FLAGS='-Wl,-h -Wl,$@'
dnl    RLD_FLAG=
    SH_LDFLAGS="-G"
  ;;
  rs6000-ibm-aix* | powerpc-ibm-aix*)
dnl    CPICFLAG=
dnl    CXXPICFLAG=
dnl    FPICFLAG=
dnl    DLFCN_DIR=dlfcn
    SH_LDFLAGS="-G $ACX_LOPT_FLAGS"
    SH_LINKING_TO_FLAGS="-brtl -bexpall -bshared"
dnl    use_ldaix="yes"
dnl    AC_SUBST(use_ldaix)
  ;;
  hppa*-hp-hpux*)
dnl    if test "$ac_cv_f77_compiler_gnu" = yes; then
dnl      FPICFLAG=-fPIC
dnl    else
dnl      FPICFLAG=+Z
dnl    fi
dnl    SHLEXT=sl
    SH_LDFLAGS="-shared -fPIC"
  ;;
  *-sgi-*)
      true
  ;;
  sparc-sun-sunos4*)
    SH_LD=ld
    SH_LDFLAGS="-assert nodefinitions"
  ;;
  sparc-sun-solaris2* | i386-pc-solaris2*)
    if test "$GXX" != yes; then
      SH_LDFLAGS=-G
    fi
  ;;
esac

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

  AC_CHECK_HEADER(Mach-O/dyld.h)  
  if test "$ac_cv_header_Mach_O_dyld_h" = yes; then
    dyld_api=true
  else 
    AC_CHECK_LIB(dld, shl_load)
    AC_CHECK_FUNCS(shl_load shl_findsym)
    if test "$ac_cv_func_shl_load" = yes \
      && test "$ac_cv_func_shl_findsym" = yes; then
      shl_load_api=true
    else
      AC_CHECK_LIB(wsock32, LoadLibrary)
      AC_CHECK_FUNCS(LoadLibrary)
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
])
