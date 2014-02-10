// MVe 07.10.98
#if defined(WIN32)
#define CDECL __cdecl
#else
#define CDECL 
#endif

#include "../config.h"

#ifndef AIX
#ifndef PATH_MAX
#define PATH_MAX 1024
#endif
#endif


