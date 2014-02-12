

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

#ifdef MODULE_MATC

double str_p[STR_MAXVALS];
char str_pstr[STR_MAXLEN];

#else

extern double str_p[STR_MAXVALS];
extern char str_pstr[STR_MAXLEN];

#endif

