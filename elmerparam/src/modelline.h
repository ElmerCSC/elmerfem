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
#ifndef _MODELLINE_H_
#define _MODELLINE_H_

#include "global.h"

typedef enum { ML_LITERAL, ML_WHITESPACE, ML_PARAM } mltype_t;


/* Used to represent both ML_LITERAL:s and ML_WHIESPACE:s.  ML_WHITESPACE is
 * treated exactly the same way as ML_LITERAL when printing, but when reading,
 * it's length is insignificant; whitespace up to next non-whitespace character
 * is skipped. */

typedef struct {
    char s[MAXLINESIZE];
    int len;
} ml_literal_t;

typedef struct {
    char type;   /* I, R, P, T or O */
    daint_t *index;
    int len;     /* Vector length; 0 for "all", -1 for error during parsing.  */
    int column;
} ml_param_t;

typedef struct ml_node_t {
    mltype_t type;
    union {
        ml_literal_t l;
        ml_param_t p;
    } u;
    struct ml_node_t *next;
} ml_node_t;

typedef struct {
    ml_node_t *line;    /* Linked list of nodes representing one line. */

    /* The following components are used for error messages.  */
    char *fname;
    int lnr;
} modelline_t;

modelline_t *ml_parse(const char *line, const char *fname, int lnr);
void ml_print(modelline_t *ml, FILE *fd, const param_t *p);
void ml_read(modelline_t *ml, FILE *fd, param_t *p);
void ml_kill(modelline_t *ml);

#endif                          /* _MODELLINE_H_ */
