dnl ----------------------------------------------------
dnl        Platform specific autoconf tests
dnl
dnl    This is part of the input for the ./configure script of
dnl    the deal.II libraries. All options and paths are stored in
dnl    the file common/Make.global_options.
dnl
dnl    In doc/Makefile some information on the kind of documentation
dnl    is stored.
dnl
dnl
dnl Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008 by the deal.II authors
dnl
dnl $Id: aclocal.m4 17530 2008-11-10 14:55:13Z bangerth $



dnl ------------------------------------------------------------
dnl Check whether Metis is installed, and if so store the 
dnl respective links
dnl
dnl Usage: DEAL_II_CONFIGURE_METIS
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CONFIGURE_METIS, dnl
[
  dnl First check for the Metis directory

  AC_ARG_WITH(metis,
  [  --with-metis=/path/to/metis  Specify the path to the Metis installation,
                          of which the include and library directories
                          are subdirs; use this if you want to override
                          the METIS_DIR environment variable.],
     [
	USE_CONTRIB_METIS=yes
        DEAL_II_METIS_DIR=$withval
	AC_MSG_RESULT($DEAL_II_METIS_DIR)

        dnl Make sure that what was specified is actually correct
        if test ! -d $DEAL_II_METIS_DIR/Lib ; then
          AC_MSG_WARN([Path to Metis specified with --with-metis does not
 			point to a complete Metis installation])
	fi
	
	DEAL_II_METIS_LIBDIR="$DEAL_II_METIS_DIR"
	DEAL_II_METIS_INCDIR="$DEAL_II_METIS_DIR/Lib"
     ],
     [
        dnl Take something from the environment variables, if it is there
        if test "x$METIS_DIR" != "x" ; then
  	  USE_CONTRIB_METIS=yes
          DEAL_II_METIS_DIR="$METIS_DIR"
	  AC_MSG_RESULT($DEAL_II_METIS_DIR)

          dnl Make sure that what this is actually correct
          if test ! -d $DEAL_II_METIS_DIR/Lib ; then
            AC_MSG_ERROR([The path to Metis specified in the METIS_DIR
	  		  environment variable does not
 			  point to a complete Metis installation])
	  fi
	  DEAL_II_METIS_LIBDIR="$DEAL_II_METIS_DIR"
	  DEAL_II_METIS_INCDIR="$DEAL_II_METIS_DIR/Lib"
        else
	  USE_CONTRIB_METIS=no
          DEAL_II_METIS_DIR=""
        fi
     ])

  AC_ARG_WITH(metis-include,
  [  --with-metis-include=/path/to/metis  Specify the path to the METIS headers file; 
                          use this if you want to override the
                          METIS_INCLUDE_DIR environment variable.],
     [
        METIS_INCLUDE_DIR=$withval
	DEAL_II_METIS_INCDIR="$METIS_INCLUDE_DIR"
     ])
     
  AC_ARG_WITH(metis-libs,
  [  --with-metis-libs=/path/to/metis  Specify the path to the METIS libraries; 
                          use this if you want to override the
                          METIS_LIBDIR environment variable.],
     [
	USE_CONTRIB_METIS=yes
        DEAL_II_METIS_LIBDIR=$withval
	AC_MSG_RESULT($DEAL_II_METIS_LIBDIR)

        dnl Make sure that what was specified is actually correct
        if test ! -d $DEAL_II_METIS_LIBDIR ; then
          AC_MSG_ERROR([Path to Metis specified with --with-metis does not
 			point to a complete Metis installation])
	fi
     ],
     [
        dnl Take something from the environment variables, if it is there
        if test "x$METIS_LIBDIR" != "x" ; then
  	  USE_CONTRIB_METIS=yes
          DEAL_II_METIS_LIBDIR="$METIS_LIBDIR"
	  AC_MSG_RESULT($DEAL_II_METIS_LIBDIR)

          dnl Make sure that what this is actually correct
          if test ! -d $DEAL_II_METIS_LIBDIR ; then
            AC_MSG_ERROR([The path to Metis specified in the METIS_DIR
	  		  environment variable does not
 			  point to a complete Metis installation])
	  fi
        else
          dnl Unless --with-metis has been set before, declare that METIS
	  dnl is not desired.
          if test "x$USE_CONTRIB_METIS" != "xyes" ; then
  	    USE_CONTRIB_METIS=no
            DEAL_II_METIS_LIBDIR=""
          fi
        fi
     ])
     
  if test "x$USE_CONTRIB_METIS" = "xyes" ; then
    AC_DEFINE(DEAL_II_USE_METIS, 1,
              [Defined if a Metis installation was found and is going
               to be used])
    LDFLAGS="$LDFLAGS -L$DEAL_II_METIS_LIBDIR -lmetis"

    dnl AC_MSG_CHECKING(for Metis version)
    dnl DEAL_II_METIS_VERSION=`cat $DEAL_II_METIS_DIR/VERSION`
    dnl AC_MSG_RESULT($DEAL_II_METIS_VERSION)
  fi
])
