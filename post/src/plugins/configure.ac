AC_INIT(moduleconfig, version-1.0)

AC_PROG_CC

# Common headers and libs:
#-----------------------------------------------------------------------------
acx_save_CFLAGS="$CFLAGS"
acx_save_LIBS="$LIBS"

AC_CHECK_HEADERS([GL/gl.h],[],[AC_ERROR([GL/gl.h not found])])
AC_CHECK_HEADERS([tcl.h],[],[AC_ERROR([tcl.h not found])])
AC_CHECK_HEADERS([windows.h])

if test "$ac_cv_header_windows_h" == yes; then
  AC_MSG_CHECKING([for glBegin in -lopengl32])

  LIBS="$LIBS -lopengl32"

  AC_LINK_IFELSE(
    [
      AC_LANG_PROGRAM(
        [
          [
#include <windows.h>
#include <GL/gl.h>
          ]
        ],
        [
          [
glBegin(0)
          ]
        ]
      )
    ],
    [
      ax_cv_lib_opengl32_glBegin=yes
      AC_MSG_RESULT([yes])
    ],
    [
      ax_cv_lib_opengl32_glBegin=no
      LIBS="$acx_save_LIBS"
      AC_MSG_RESULT([no])
      AC_ERROR([-lopengl32 not found])
    ]
  )

  AC_CHECK_LIB([tcl84],[Tcl_CreateCommand])
  if test "$ac_cv_lib_tcl84_Tcl_CreateCommand" == no; then
    AC_ERROR([-ltcl84 not found])
  fi

  AC_CHECK_LIB([z],[deflate])
  if test "$ac_cv_lib_z_deflate" == no; then
    AC_ERROR([-lz not found])
  fi
fi

COMMON_LIBS="$LIBS"
COMMON_CFLAGS="$CFLAGS"

CFLAGS="$acx_save_CFLAGS"
LIBS="$acx_save_LIBS"

# Savejpg:
#-----------------------------------------------------------------------------
acx_save_CFLAGS="$CFLAGS"
acx_save_LIBS="$LIBS"

LIBS="$COMMON_LIBS"
CFLAGS="$COMMON_CFLAGS"

AC_CHECK_HEADERS([jpeglib.h],[AC_CHECK_LIB([jpeg],[jpeg_CreateCompress])])

JPG_LIBS="$LIBS $JPG_LIBS"
JPG_CFLAGS="$CFLAGS $JPG_CFLAGS"

CFLAGS="$acx_save_CFLAGS"
LIBS="$acx_save_LIBS"

# Savepng:
#-----------------------------------------------------------------------------
acx_save_CFLAGS="$CFLAGS"
acx_save_LIBS="$LIBS"

LIBS="$COMMON_LIBS"
CFLAGS="$COMMON_CFLAGS"

AC_CHECK_HEADERS([png.h],[AC_CHECK_LIB([png],[png_create_write_struct])])

PNG_CFLAGS="$CFLAGS $PNG_CFLAGS"
PNG_LIBS="$LIBS $PNG_LIBS"

CFLAGS="$acx_save_CFLAGS"
LIBS="$acx_save_LIBS"

# Savetiff:
#-----------------------------------------------------------------------------
acx_save_CFLAGS="$CFLAGS"
acx_save_LIBS="$LIBS"

LIBS="$COMMON_LIBS"
CFLAGS="$COMMON_CFLAGS"

AC_CHECK_HEADERS([tiffio.h],[AC_CHECK_LIB([tiff],[TIFFOpen])])

TIFF_CFLAGS="$CFLAGS $TIFF_LIBS"
TIFF_LIBS="$LIBS $TIFF_LIBS"

CFLAGS="$acx_save_CFLAGS"
LIBS="$acx_save_LIBS"

# Savempg:
#-----------------------------------------------------------------------------
acx_save_CFLAGS="$CFLAGS"
acx_save_LIBS="$LIBS"

LIBS="$COMMON_LIBS"
CFLAGS="$COMMON_CFLAGS"

AC_CHECK_HEADERS([math.h])
AC_CHECK_HEADERS([avcodec.h])

if test "$ac_cv_header_math_h" == yes; then
  if test "$ac_cv_header_avcodec_h" == yes; then
    AC_CHECK_LIB([avcodec],[avcodec_init])
    AC_CHECK_LIB([avutil],[av_malloc])
    AC_CHECK_LIB([swscale],[sws_scale])
  fi
fi

MPG_CFLAGS="$CFLAGS $MPG_LIBS"
MPG_LIBS="$LIBS $MPG_LIBS"

CFLAGS="$acx_save_CFLAGS"
LIBS="$acx_save_LIBS"

# Finalize:
#-----------------------------------------------------------------------------
INSTALL_DIR="${prefix}/modules"

if test "$ac_cv_header_windows_h" == yes; then
  SHL_EXT="dll"
else
  # For linux, always make pic, even for 32 bits:
  JPG_CFLAGS="$JPG_CFLAGS -fPIC"
  PNG_CFLAGS="$PNG_CFLAGS -fPIC"
  TIFF_CFLAGS="$TIFF_CFLAGS -fPIC"
  MPG_CFLAGS="$MPG_CFLAGS -fPIC"
  SHL_EXT="so"
fi

AC_SUBST([INSTALL_DIR])
AC_SUBST([SHL_EXT])

AC_MSG_RESULT([
     INSTALL_DIR   $INSTALL_DIR
     SHL_EXT       .$SHL_EXT
])

if test "$ac_cv_lib_jpeg_jpeg_CreateCompress" == yes; then
  JPG_OK=yes;
  AC_MSG_RESULT([
     JPG_CFLAGS    $JPG_CFLAGS
     JPG_LIBS      $JPG_LIBS
  ])
else
  AC_MSG_RESULT([
     JPG           DISABLED
  ])
fi
AC_SUBST([JPG_CFLAGS])
AC_SUBST([JPG_LIBS])
AC_SUBST([JPG_OK])


if test "$ac_cv_lib_png_png_create_write_struct" == yes; then
  PNG_OK=yes;
  AC_MSG_RESULT([
     PNG_CFLAGS    $PNG_CFLAGS
     PNG_LIBS      $PNG_LIBS
  ])
else
  AC_MSG_RESULT([
     PNG           DISABLED
  ])
fi
AC_SUBST([PNG_CFLAGS])
AC_SUBST([PNG_LIBS])
AC_SUBST([PNG_OK])

if test "$ac_cv_lib_tiff_TIFFOpen" == yes; then
  TIFF_OK=yes;
  AC_MSG_RESULT([
     TIFF_CFLAGS   $TIFF_CFLAGS
     TIFF_LIBS     $TIFF_LIBS
  ])
else
  AC_MSG_RESULT([
     TIFF          DISABLED
  ])
fi
AC_SUBST([TIFF_CFLAGS])
AC_SUBST([TIFF_LIBS])
AC_SUBST([TIFF_OK])

if test "$ac_cv_lib_avutil_av_malloc" == yes; then
  if test "$ac_cv_lib_avcodec_avcodec_init" == yes; then
    MPG_OK=yes;
    AC_MSG_RESULT([
     MPG_CFLAGS    $MPG_CFLAGS
     MPG_LIBS      $MPG_LIBS
    ])
  else
    AC_MSG_RESULT([
     MPG           DISABLED
    ])
  fi
else
  AC_MSG_RESULT([
     MPG           DISABLED
  ])
fi
AC_SUBST([MPG_CFLAGS])
AC_SUBST([MPG_LIBS])
AC_SUBST([MPG_OK])


AC_OUTPUT(Makefile)
