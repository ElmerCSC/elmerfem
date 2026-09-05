/*
 * elmerf90 - wrapper to compile Elmer solver plugins as shared libraries.
 *
 * Replaces the elmerf90 shell script (Linux/Mac) and elmerf90.bat (Windows).
 *
 * WHAT IT DOES
 *   1. Resolves LIBDIR / INCLUDE, first match wins:
 *        ELMER_LIB          -> libdir = $ELMER_LIB, incdir = $ELMER_LIB/../include
 *        ELMER_HOME         -> libdir = $ELMER_HOME/<install libdir>,
 *                              incdir = $ELMER_HOME/share/elmersolver/include
 *        neither            -> the paths baked in at build time
 *   2. Resolves the Fortran compiler:
 *        ELMER_Fortran_COMPILER if set and non-empty, else CMAKE_Fortran_COMPILER
 *        as recorded when Elmer was built. See exec_compiler() in the helper for
 *        the PATH fallback when that build-time compiler is not present here.
 *   3. Assembles one command line, in this order:
 *        <compiler> <caller's argv verbatim> <baked-in flags> -I<inc> -L<lib>
 *        [ElmerIce libs and -rpath] [MMG/ParMMG -I/-L] -lelmersolver
 *      The baked-in flags are CMAKE_Fortran_FLAGS, ELMER_F90FLAGS and the two
 *      CMake shared-library Fortran flags, concatenated at configure time and
 *      re-split on whitespace here.
 *   4. Prints the assembled command to stdout, then execs it.
 *
 * WHAT IT DOES NOT DO
 *   - It does not compile or link anything itself. It is argv assembly plus
 *     exec, so the compiler replaces this process and its exit status is what
 *     the caller sees. There is no fork, no wait, no post-processing.
 *   - It does not parse or validate the caller's arguments. It cannot tell a
 *     compile-only (-c) invocation from a full link, and adds the same flags,
 *     -L paths and -lelmersolver either way. Harmless for -c, but it means the
 *     wrapper has no notion of what you are actually asking for.
 *   - It does not quote or escape. Baked-in flags are split on whitespace, so
 *     an install path containing spaces produces a broken command line. Same
 *     limitation as the shell script it replaced.
 *   - It does not verify that the compiler it runs is the one that built
 *     libelmersolver, nor that the two are ABI- or module-format compatible.
 *     Fortran .mod files are compiler- and often version-specific, so a
 *     mismatch typically surfaces as an unreadable module rather than a clear
 *     error from here.
 *   - It does not check that libdir/incdir exist or contain anything.
 *   - It does not add -o. The caller owns output naming. Note it *does* add
 *     -fPIC and -shared, since CMAKE_SHARED_LIBRARY_Fortran_FLAGS and
 *     CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS are part of the baked-in flag
 *     string -- so every invocation is a shared-object build, including one
 *     you meant to be compile-only.
 *   - It caps the command line at MAX_ARGS (1024) entries and exits if exceeded.
 *
 * QUIRKS WORTH KNOWING
 *   - The assembled command goes to stdout while the "with/no elmerice" and
 *     MMG notes go to stderr, so capturing stdout yields the command line
 *     alone. The elmerice note is printed on every run, including "no
 *     elmerice", which makes the wrapper chatty on stderr by default.
 *   - That command is echoed before the exec is attempted, so it always names
 *     the build-time compiler. If the PATH fallback kicks in, the binary
 *     actually run differs from the one printed -- exec_compiler() reports the
 *     substitution on stderr for exactly this reason.
 */

#include "elmerf90-helper.h"

/* --- main --------------------------------------------------------------- */

int main(int argc, char *argv[])
{
    const char *env_lib  = getenv("ELMER_LIB");
    const char *env_home = getenv("ELMER_HOME");
    const char *env_fc   = getenv("ELMER_Fortran_COMPILER");

    /* Resolve LIBDIR and INCLUDE (mirrors shell-script priority) */
    char *libdir, *incdir;
    if (env_lib && *env_lib) {
        libdir = strdup(env_lib);
        incdir = join(env_lib, "/../include");
    } else if (env_home && *env_home) {
        libdir = join(env_home, "/" ELMERF90_INSTALL_LIB);
        incdir = join(env_home, "/share/elmersolver/include");
    } else {
        libdir = strdup(ELMERF90_LIBDIR);
        incdir = strdup(ELMERF90_INCLUDE_DIR);
    }

    const char *fc = (env_fc && *env_fc) ? env_fc : ELMERF90_FC;

    /* --- Build argv ---------------------------------------------------- */

    push(fc);

    /* Pass caller's arguments through unchanged */
    for (int i = 1; i < argc; i++)
        push(argv[i]);

    /* Baked-in compiler flags */
    push_flags(ELMERF90_FFLAGS);

    push(join("-I", incdir));
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
    if (*ELMERF90_MMG_INCDIR) push(join("-I", ELMERF90_MMG_INCDIR));
    if (*ELMERF90_MMG_LIBDIR) push(join("-L", ELMERF90_MMG_LIBDIR));
#endif

#if ELMERF90_HAVE_PARMMG
    fprintf(stderr, "with ParMMG\n");
    {
        int same_lib = strcmp(ELMERF90_PARMMG_LIBDIR, ELMERF90_MMG_LIBDIR) == 0;
        int same_inc = strcmp(ELMERF90_PARMMG_INCDIR, ELMERF90_MMG_INCDIR) == 0;
        if (!same_lib) {
            if (*ELMERF90_PARMMG_INCDIR) push(join("-I", ELMERF90_PARMMG_INCDIR));
            if (*ELMERF90_PARMMG_LIBDIR) push(join("-L", ELMERF90_PARMMG_LIBDIR));
        } else {
            if (!same_inc && *ELMERF90_PARMMG_INCDIR)
                push(join("-I", ELMERF90_PARMMG_INCDIR));
            fprintf(stderr, "MMG and ParMMG share the same lib dir\n");
        }
        if (same_inc) fprintf(stderr, "MMG and ParMMG share the same include dir\n");
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

    return exec_compiler(fc, "elmerf90");
}
