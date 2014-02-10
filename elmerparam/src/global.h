/*

  ElmerParam - A simple system for parametrized computing
 
  Copyright (C) 2006  CSC - IT Center for Science Ltd.

  Authors: Erik Edelmann <Erik.Edelmann@csc.fi>
           Peter Råback <Peter.Raback@csc.fi>
  Address: CSC - IT Center for Science Ltd.
           Keilaranta 14
           02101 Espoo, Finland
            
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program (in file elmerparam/COPYING); if not, write to
  the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
  Boston, MA 02110-1301, USA.

 */

#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include <float.h>
#include <stdio.h>
#include "config.h"
#include "dynarray.h"

#define FALSE 0
#define TRUE 1
#define MAXLINESIZE 512
#define MAXFILESIZE 512
#define DEFAULT_FUN DBL_MAX
#define PKG_NAME "ElmerParam: "

#ifndef DISABLE_MATC
    char *mtc_domath(const char *);
    void mtc_init(FILE * input, FILE * output, FILE * error);

#   define MTC_DOMATH(cmd) mtc_domath(cmd)
#   define MTC_INIT(p) {\
        char command[MAXLINESIZE];\
        int i;\
    \
        mtc_init(NULL, stdout, stderr);\
        strcpy(command, "format( 12, \"rowform\")");\
        mtc_domath(command);\
        for (i = 0; i < da_n(p->fun); i++) {\
            sprintf(command, "O(%d) = %e", i, dr_get(p->fun,i));\
            mtc_domath(command);\
        }\
        for (i = 0; i < da_n(p->xr); i++) {\
            sprintf(command, "R(%d) = %e", i, dr_get(p->xr,i));\
            mtc_domath(command);\
        }\
        for (i = 0; i < da_n(p->xi); i++) {\
            sprintf(command, "I(%d) = %d", i, di_get(p->xi,i));\
            mtc_domath(command);\
        }\
        printf("MATC library was activated!\n");\
        p->usematc = TRUE;\
    }
#else
static void *nop() { return NULL; }
#   define MTC_DOMATH(cmd) nop()
#   define MTC_INIT(p) {\
        fprintf(stderr, "WARNING: This version of ElmerParam was compiled "\
                        "without MATC library!\n");\
        p->usematc = FALSE;\
    }
#endif

typedef struct {
    daint_t *xi;
    dareal_t *xr;
    dareal_t *fun;

    int info, usematc;

    int taglen;
    char tag[MAXLINESIZE];

    char cmdfile[MAXFILESIZE];
    int lnr;                    /* Line number in the command file. */
} param_t;

#endif                          /* _GLOBAL_H_ */
