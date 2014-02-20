#ifndef CONFIG_H_EIO
#define CONFIG_H_EIO

#define STDCALLBULL

/* #include "FCMangle.h" */

#cmakedefine HAVE_INTTYPES_H
/* #cmakedefine FC_FUNC ${FC_FUNC} */
/* #cmakedefine FC_FUNC_ ${FC_FUNC_} */
/* #cmakedefine FC_CHAR_PTR${FC_CHAR_PTR} */

#ifdef USE_ISO_C_BINDINGS
#define FC_FUNC(name, NAME) name
#define FC_FUNC_(name, NAME) name
#endif

#endif
