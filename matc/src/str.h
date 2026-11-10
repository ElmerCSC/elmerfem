

/*
 * $Id: str.h,v 1.2 1998/08/01 12:34:56 jpr Exp $ 
 *
 * $Log: str.h,v $
 * Revision 1.2  1998/08/01 12:34:56  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#define STR_MAXVALS 32
#define STR_MAXLEN 512

/*
 * str_p / str_pstr were formerly file-scope OpenMP threadprivate scratch
 * arrays shared by str_sprintf/str_sscanf (str.c) and fil_fscanf/fil_fgets
 * (files.c). gcc 16 miscompiles stores to these threadprivate arrays at -O1
 * and above: the store does not persist (a subsequent read sees a stale
 * zero), so e.g. sprintf("%g",i) formatted 0 for every i, which corrupted the
 * MATC-generated 'Exported Variable N ...' keywords in SOLVER.KEYWORDS. They
 * are now plain function-local arrays in each user, which is inherently
 * thread-safe and sidesteps the miscompilation.
 */

