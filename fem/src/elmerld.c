/*
 * elmerld - wrapper to link Elmer solver plugins as shared libraries.
 *
 * Logic:
 *   1. Determine LIBDIR from env (ELMER_LIB, ELMER_HOME) or baked-in paths.
 *   2. Determine the Fortran compiler from ELMER_Fortran_COMPILER or the build-time default.
 *   3. Optionally append ElmerIce / MMG / ParMMG flags.
 *   4. exec the linker (the compiler is used as the linker driver) with the assembled argument list.
 */

#include "elmerf90-helper.h"

/* --- main --------------------------------------------------------------- */

int main(int argc, char *argv[])
{
    const char *env_lib  = getenv("ELMER_LIB");
    const char *env_home = getenv("ELMER_HOME");
    const char *env_fc   = getenv("ELMER_Fortran_COMPILER");

    /* Resolve LIBDIR and INCLUDE (mirrors shell-script priority) */
    char *libdir;
    if (env_lib && *env_lib) {
        libdir = strdup(env_lib);
    } else if (env_home && *env_home) {
        libdir = join(env_home, "/" ELMERF90_INSTALL_LIB);
    } else {
        libdir = strdup(ELMERF90_LIBDIR);
    }

    const char *fc = (env_fc && *env_fc) ? env_fc : ELMERF90_FC;

    /* --- Build argv ---------------------------------------------------- */

    push(fc);

    /* Pass caller's arguments through unchanged */
    for (int i = 1; i < argc; i++)
        push(argv[i]);

    /* Baked-in compiler flags */
    push_flags(ELMERF90_FFLAGS);

    push(join("-L", libdir));

#if ELMERF90_HAVE_ELMERICE
    {
        char *elib = join(libdir, "/../../share/elmersolver/lib");
        fprintf(stderr, "with elmerice\n");
#ifndef _WIN32
        /* rpath is a Linux/Mac concept; Windows finds DLLs via PATH */
        push("-Xlinker");
        push(join("-rpath=", elib));
#endif
        push(join("-L", elib));
        free(elib);
        push(join("-l:", ELMERF90_ELMERICESOLVERS_LIB));
        push(join("-l:", ELMERF90_ELMERICEUSF_LIB));
        push(join("-l:", ELMERF90_ELMERICEUTILS_LIB));
    }
#else
    fprintf(stderr, "no elmerice\n");
#endif

#if ELMERF90_HAVE_MMG
    fprintf(stderr, "with MMG\n");
    if (*ELMERF90_MMG_LIBDIR) push(join("-L", ELMERF90_MMG_LIBDIR));
#endif

#if ELMERF90_HAVE_PARMMG
    fprintf(stderr, "with ParMMG\n");
    {
        int same_lib = strcmp(ELMERF90_PARMMG_LIBDIR, ELMERF90_MMG_LIBDIR) == 0;
        if (!same_lib) {
            if (*ELMERF90_PARMMG_LIBDIR) push(join("-L", ELMERF90_PARMMG_LIBDIR));
        } else {
            fprintf(stderr, "MMG and ParMMG share the same lib dir\n");
        }
    }
#endif

    push("-lelmersolver");
    args[nargs] = NULL;

    /* Print the command (matches shell script behaviour) */
    for (int i = 0; i < nargs; i++)
#if defined (_WIN32)
        /* The argument list is already double-quoted for Windows. */
        printf("%s ", args[i]);
#else
        printf("\"%s\" ", args[i]);
#endif
    printf("\n");

    exec_compiler(fc, "elmerld");
}
