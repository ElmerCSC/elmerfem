#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.56)
AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
AC_CONFIG_SRCDIR([src/3d.h])
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE

AC_PROG_RANLIB

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h])

AC_CHECK_HEADERS([windows.h],[],[acx_win=no])
if test "$acx_win" = "no"; then
   AC_MSG_ERROR([windows.h not found.

If this is not a windows machine with MinGW installed,
you probably don't want to build this package anyways.
Glaux is distributed with all Linuxes. It is also part
of MS SDK.
                ])
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([sqrt])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
