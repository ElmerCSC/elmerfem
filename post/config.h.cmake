/* config.h.in.  Generated from configure.in by autoheader.  */

/* Detected platform. */
#undef AIX

/* Couldn't determine. sticking with 32 bits. */
#undef ARCH_32_BITS

/* 64 bit arch. */
#undef ARCH_64_BITS

/* Detected platform. */
#undef BASTARDS

/* Detected platform. */
#undef BSD

/* Detected platform. */
#undef CYGWIN

/* Standard windows call declaration */
#undef C_DLLEXPORT

/* Detected platform. */
#undef DARWIN

/* Detected platform. */
#undef DEC_ALPHA

/* Elmer post default install directory */
/* #undef ELMER_POST_HOME */

/* Define if using dynamic linking */
#undef ENABLE_DYNAMIC_LINKING

/* Char pointer mangling */
#define FC_CHAR_PTR(P,L) char *P, int L

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
#undef FC_DUMMY_MAIN

/* Define if F77 and FC dummy `main' functions are identical. */
#undef FC_DUMMY_MAIN_EQ_F77

#ifdef USE_ISO_C_BINDINGS
/* #define FC_FUNC(name,NAME) name ## _ */
#define FC_FUNC(name_gem,NAME_cap) name_gem ## _
#define FC_FUNC_(name_gem,NAME_cap) name_gem ## _
#endif

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
/* #undef FC_FUNC */

/* As FC_FUNC, but for C identifiers containing underscores. */
/* #undef FC_FUNC_ */

/* Use the Apple OpenGL framework. */
#undef HAVE_APPLE_OPENGL_FRAMEWORK

/* Define to 1 if you have the `dlclose' function. */
#undef HAVE_DLCLOSE

/* Define to 1 if you have the `dlerror' function. */
#undef HAVE_DLERROR

/* Define to 1 if you have the `dlopen' function. */
#undef HAVE_DLOPEN

/* Define if your system has dlopen, dlsym, dlerror, and dlclose for dynamic
   linking */
#undef HAVE_DLOPEN_API

/* Define to 1 if you have the `dlsym' function. */
#undef HAVE_DLSYM

/* Define if your system has dyld for dynamic linking */
#undef HAVE_DYLD_API

/* Define to 1 if you have the `floor' function. */
#undef HAVE_FLOOR

/* Define if you have a FTGL library (new style). */
#undef HAVE_FTGL_NEW

/* Define if you have a FTGL library (old style). */
#undef HAVE_FTGL_OLD

/* Define to 1 if you have the `gettimeofday' function. */
#undef HAVE_GETTIMEOFDAY

/* Define to 1 if you have the <inttypes.h> header file. */
#undef HAVE_INTTYPES_H

/* Define to 1 if you have the `dl' library (-ldl). */
#undef HAVE_LIBDL

/* Define to 1 if you have the `dld' library (-ldld). */
#undef HAVE_LIBDLD

/* Define if your system has LoadLibrary for dynamic linking */
#undef HAVE_LOADLIBRARY_API

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
#undef HAVE_MALLOC

/* Define to 1 if you have the <malloc.h> header file. */
#undef HAVE_MALLOC_H

/* Define if you have a MATC library. */
#undef HAVE_MATC

/* Define to 1 if you have the <memory.h> header file. */
#undef HAVE_MEMORY_H

/* Define to 1 if you have the `memset' function. */
#undef HAVE_MEMSET

/* Define if you have POSIX threads libraries and header files. */
#undef HAVE_PTHREAD

/* Define to 1 if your system has a GNU libc compatible `realloc' function,
   and to 0 otherwise. */
#undef HAVE_REALLOC

/* Define to 1 if you have the `shl_findsym' function. */
#undef HAVE_SHL_FINDSYM

/* Define to 1 if you have the `shl_load' function. */
#undef HAVE_SHL_LOAD

/* Define if your system has shl_load and shl_findsym for dynamic linking */
#undef HAVE_SHL_LOAD_API

/* Define to 1 if you have the `sqrt' function. */
#undef HAVE_SQRT

/* Define to 1 if `stat' has the bug that it succeeds when given the
   zero-length file name argument. */
#undef HAVE_STAT_EMPTY_STRING_BUG

/* Define to 1 if stdbool.h conforms to C99. */
#undef HAVE_STDBOOL_H

/* Define to 1 if you have the <stdint.h> header file. */
#undef HAVE_STDINT_H

/* Define to 1 if you have the <stdlib.h> header file. */
#undef HAVE_STDLIB_H

/* Define to 1 if you have the `strerror' function. */
#undef HAVE_STRERROR

/* Define to 1 if you have the <strings.h> header file. */
#cmakedefine HAVE_STRINGS_H

/* Define to 1 if you have the <string.h> header file. */
#cmakedefine HAVE_STRING_H

/* Define to 1 if you have the <sys/param.h> header file. */
#undef HAVE_SYS_PARAM_H

/* Define to 1 if you have the <sys/stat.h> header file. */
#undef HAVE_SYS_STAT_H

/* Define to 1 if you have the <sys/time.h> header file. */
#undef HAVE_SYS_TIME_H

/* Define to 1 if you have the <sys/types.h> header file. */
#undef HAVE_SYS_TYPES_H

/* Define to 1 if you have the <termio.h> header file. */
#undef HAVE_TERMIO_H

/* Define to 1 if you have the <unistd.h> header file. */
#undef HAVE_UNISTD_H

/* Define to 1 if you have the <windows.h> header file. */
#cmakedefine HAVE_WINDOWS_H

/* Define to 1 if the system has the type `_Bool'. */
#undef HAVE__BOOL

/* Detected platform. */
#undef HPUX

/* Detected platform. */
#undef LINUX

/* Define to 1 if `lstat' dereferences a symlink specified with a trailing
   slash. */
#undef LSTAT_FOLLOWS_SLASHED_SYMLINK

/* Detected platform. */
#cmakedefine MINGW32
#cmakedefine __WIN32__ @WIN32@
#cmakedefine WIN32

/* Name of package */
#undef PACKAGE

/* Define to the address where bug reports for this package should be sent. */
#undef PACKAGE_BUGREPORT

/* Define to the full name of this package. */
#undef PACKAGE_NAME

/* Define to the full name and version of this package. */
#undef PACKAGE_STRING

/* Define to the one symbol short name of this package. */
#undef PACKAGE_TARNAME

/* Define to the version of this package. */
#undef PACKAGE_VERSION

/* Define to necessary symbol if this constant uses a non-standard name on
   your system. */
#undef PTHREAD_CREATE_JOINABLE

/* Define as the return type of signal handlers (`int' or `void'). */
#undef RETSIGTYPE

/* Detected platform. */
#undef SGI

/* Shared lib filename extension */
#undef SHL_EXTENSION

/* The size of `void*', as computed by sizeof. */
#undef SIZEOF_VOIDP

/* Detected platform. */
#undef SOLARIS

/* Standard windows call declaration */
#cmakedefine STDCALLBULL @FC_STDCALLBULL@

/* Define to 1 if you have the ANSI C header files. */
#undef STDC_HEADERS

/* Detected platform. */
#undef SUNOS

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
#undef TIME_WITH_SYS_TIME

/* Version number of package */
#undef VERSION

/* Detected platform2. */
/* #undef WIN32 */

/* Number of bits in a file offset, on hosts where this is settable. */
#undef _FILE_OFFSET_BITS

/* Define for large files, on AIX-style hosts. */
#undef _LARGE_FILES

/* Define to rpl_malloc if the replacement function should be used. */
#undef malloc

/* Define to rpl_realloc if the replacement function should be used. */
#undef realloc

/* Define to `unsigned int' if <sys/types.h> does not define. */
#undef size_t

/* Workaround for Tcl_Interp->result*/
#cmakedefine USE_INTERP_RESULT

/* Set ELMER_POST_HOME to install prefix */ 
#define ELMER_POST_HOME "@CMAKE_INSTALL_PREFIX@/@ELMER_POST_DATADIR@"
