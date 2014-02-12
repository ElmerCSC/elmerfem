#ifndef CONFIG_H_UMFPACK
#define CONFIG_H_UMFPACK

#include "FCMangle.h"

/* Couldn't determine. sticking with 32 bits. */
#cmakedefine ARCH_32_BITS

/* 64 bit arch. */
#cmakedefine ARCH_64_BITS

/* Standard windows call declaration */
#cmakedefine C_DLLEXPORT

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
#cmakedefine F77_DUMMY_MAIN

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#cmakedefine F77_FUNC

/* As F77_FUNC, but for C identifiers containing underscores. */
#cmakedefine F77_FUNC_

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
#cmakedefine FC_DUMMY_MAIN

/* Define if F77 and FC dummy `main' functions are identical. */
#cmakedefine FC_DUMMY_MAIN_EQ_F77

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#cmakedefine FC_FUNC ${FS_FUNC}

/* As FC_FUNC, but for C identifiers containing underscores. */
#cmakedefine FC_FUNC_ ${FC_FUNC_}

/* Define Fortran string argument handling */
#cmakedefine FC_CHAR_PTR${FC_CHAR_PTR}

/* Define to 1 if you have the <inttypes.h> header file. */
#cmakedefine HAVE_INTTYPES_H

/* Define to 1 if you have the <memory.h> header file. */
#cmakedefine HAVE_MEMORY_H

/* Define to 1 if you have the <stdint.h> header file. */
#cmakedefine HAVE_STDINT_H

/* Define to 1 if you have the <stdlib.h> header file. */
#cmakedefine HAVE_STDLIB_H

/* Define to 1 if you have the <strings.h> header file. */
#cmakedefine HAVE_STRINGS_H

/* Define to 1 if you have the <string.h> header file. */
#cmakedefine HAVE_STRING_H

/* Define to 1 if you have the <sys/stat.h> header file. */
#cmakedefine HAVE_SYS_STAT_H

/* Define to 1 if you have the <sys/types.h> header file. */
#cmakedefine HAVE_SYS_TYPES_H

/* Define to 1 if you have the <unistd.h> header file. */
#cmakedefine HAVE_UNISTD_H

/* UMFPACK Shouldn't use BLAS, because of problems with the c-fortran api */
#cmakedefine NBLAS

/* The size of `void*', as computed by sizeof. */
#cmakedefine SIZEOF_VOIDP

/* Standard windows call declaration */
#cmakedefine STDCALLBULL

/* Define to 1 if you have the ANSI C header files. */
#cmakedefine STDC_HEADERS

/* Version number of package */
#cmakedefine VERSION

#endif
