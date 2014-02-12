#ifndef CONFIG_H_FEM
#define CONFIG_H_FEM

#include "FCMangle.h"

#cmakedefine VERSION "${VERSION}"
#cmakedefine HAVE_INTTYPES_H
#cmakedefine FC_FUNC ${FC_FUNC}
#cmakedefine FC_FUNC_ ${FC_FUNC_}
#cmakedefine FC_CHAR_PTR${FC_CHAR_PTR}
#cmakedefine ELMER_SOLVER_HOME "${ELMER_SOLVER_HOME}"

/* These require more work */
#define STDCALLBULL
#define ELMER_LINKTYP 1
#define ENABLE_DYNAMIC_LINKING 1

#endif
