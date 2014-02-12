#ifndef CONFIG_H_FEM
#define CONFIG_H_FEM

#include "FCMangle.h"

#cmakedefine VERSION "${VERSION}"
#cmakedefine HAVE_INTTYPES_H
#cmakedefine FC_FUNC ${FC_FUNC}
#cmakedefine FC_FUNC_ ${FC_FUNC_}
#cmakedefine FC_CHAR_PTR${FC_CHAR_PTR}
#cmakedefine ELMER_SOLVER_HOME "${ELMER_SOLVER_HOME}"

#cmakedefine OFF_KIND @OFF_KIND@

/* These require more work */
#cmakedefine STDCALLBULL @FC_STDCALLBULL@

#define STDCALLBULL

#define HAVE_ARPACK
#define HAVE_BLAS
/* Define to 1 if you have the `dlclose' function. */
#define HAVE_DLCLOSE

/* Define to 1 if you have the `dlerror' function. */
#define HAVE_DLERROR

/* Define to 1 if you have the `dlopen' function. */
#define HAVE_DLOPEN

/* Define if your system has dlopen, dlsym, dlerror, and dlclose for dynamic
   linking */
#define HAVE_DLOPEN_API

/* Define to 1 if you have the `dlsym' function. */
#define HAVE_DLSYM

/* Define if your system has dyld for dynamic linking */
#define HAVE_DYLD_API

/* Define if you have a EIOF library. */
#define HAVE_EIOF
/* Define to 1 if you have the `fseeko' function. */
#define HAVE_FSEEKO

/* Define to 1 if you have the `ftello' function. */
#define HAVE_FTELLO

/* Does the fortran environment implement etime */
#define HAVE_F_ETIME

/* Does the fortran environment implement flush */
#define HAVE_F_FLUSH

/* Define if you have a HUTI library. */
#define HAVE_HUTI

/* Define if you have a HYPRE library. */
/* define HAVE_HYPRE */

/* Define if you have LAPACK library. */
#define HAVE_LAPACK

/* Define to 1 if you have the `dl' library (-ldl). */
#define HAVE_LIBDL

/* Define if you have a MATC library. */
#define HAVE_MATC

/* ... */
#define HAVE_MPI

/* Define if you have a MUMPS library. */
/* define HAVE_MUMPS */

/* Define if you have a PARPACK library. */
#define HAVE_PARPACK

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H

/* Define if you have a UMFPACK library. */
#define HAVE_UMFPACK

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H

/* Detected platform. */
#define LINUX

/* Name of package */
#cmakedefine PACKAGE @PACKAGE@

/* Define to the full name of this package. */
#cmakedefine PACKAGE_NAME @PACKAGE_NAME@

/* Define to the version of this package. */
#cmakedefine PACKAGE_VERSION @PACKAGE_VERSION@

/* Shared lib filename extension */
#define SHL_EXTENSION ".so"


#define ELMER_LINKTYP ${ELMER_LINKTYP}
#define ENABLE_DYNAMIC_LINKING 1

#cmakedefine CONTIG ${FC_CONTIG}

#endif
