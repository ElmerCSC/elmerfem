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
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "modelline.h"


/* Read a line and return a pointer to it, and put the line length into 'len'.
 * Return NULL for EOF.  NOTE:  Next call to this function will change the
 * content of the memory area pointed to -- make a copy if you need to store the
 * read data for some longer period.  */

static const char *readline(FILE *fd, size_t *len)
{
    static char *line = NULL;
    static size_t linesize = 0;
    int i, c;

    i = 0;
    while ((c = fgetc(fd)) != '\n' && c != EOF) {
        if (i == linesize) {
            linesize = linesize ? linesize*2 : MAXLINESIZE;
            line = realloc(line, linesize);
            assert(line);  /* TODO: Prettier error checking.  */
        }

        line[i++] = c;
    }

    if (i == 0 && c == EOF) {
        if (line) free(line);
        line = NULL;
        linesize = 0;
    } else { 
        /* Make room for a '\0' at the end.  */
        if (i == linesize) {
            linesize += 1;
            line = realloc(line, linesize);
            assert(line);  /* TODO: Prettier error checking.  */
        }
        line[i] = '\0';
    }
    
    *len = i;
    return line;
}


/* Get parameter type (I, R, O, T) from s and put it in *t.  In case of
 * error, print error message to stderr and put '\0' in *t, Return number of
 * characters consumed (i.e. 1).  */

static int get_paramtype(const char *s, char *t,
                         const char *fname, int lnr, int col)
{
    char *tmp;

    if ((tmp = strchr("IROTirot", *s))) {
        *t = toupper(*tmp);
    } else {
        fprintf(stderr, PKG_NAME "%s, line %i, column %i: "
                "Parameter type '%c' not recoqnized\n", fname, lnr, col, *s);
        *t = '\0';
    }

    return 1;
}


/* Get natural number (non-neg integer) from s and put it in *i.  In case of
 * error, print error message to stderr and put < 0 in *i, Return number of
 * characters consumed.*/

static int get_natural(const char *s, int *i,
                       const char *fname, int lnr, int col)
{
    char *t;

    *i = strtol(s, &t, 10);
    if (t - s == 0) {
        fprintf(stderr,
                PKG_NAME "%s, line %i, column %i: Error reading integer\n",
                fname, lnr, col);
        *i = -1;
    } else if (*i < 0) {
        fprintf(stderr, PKG_NAME "%s, line %i, column %i: "
                "Expected non-negative integer\n", fname, lnr, col);
    }

    return t - s;
}


/* Get list of indeces from 's' and put them into 'index' and number of indeces
 * in 'len'.  In case of error, print error message to stderr and put -1 in len,
 * Return number of characters consumed.  */

static int get_vector(const char *s, int *len, daint_t **index,
                      const char *fname, int lnr, int col)
{
    int i, n;
    int begin, end;

    assert(s[0] == '(');

    n = 1;
    n += get_natural(&s[n], &begin, fname, lnr, col + n);
    if (begin < 0) {
        *len = -1;
        return 0;
    }

    if (s[n] != ':') {
        fprintf(stderr, PKG_NAME "%s, line %i, column %i: Expected ':'\n",
                fname, lnr, col + n);
        *len = -1;
        return 0;
    }
    n++;

    n += get_natural(&s[n], &end, fname, lnr, col + n);
    if (end < 0) {
        *len = -1;
        return 0;
    }

    if (s[n] != ')') {
        fprintf(stderr, PKG_NAME "%s, line %i, column %i: Expected ')'\n",
                fname, lnr, col + n);
        *len = -1;
        return 0;
    }
    n++;

    if (end < begin) {
        fprintf(stderr, PKG_NAME "%s, line %i, column %i: "
                "Starting index must be <= ending index\n", fname, lnr, col);
        *len = -1;
        return 0;
    }

    *len = end - begin + 1;
    for (i = 0; i < *len; i++)
        *index = di_set(*index, i, i + begin);

    return n;
}


/* Get the next node from line, starting at i:th character, or NULL if there are
 * no more nodes.  Update i to point to the beginning of next node element.  */

static ml_node_t *get_node(const char *line, int *i, const char *fname, int lnr)
{
    const char *s;
    ml_node_t *node;
    int j, tmp;

    if (!line[*i])
        return NULL;

    node = malloc(sizeof(ml_node_t));
    node->next = NULL;

    s = &line[*i];

    if (s[0] == '<' && s[1] == '!') {
        node->type = ML_PARAM;
        node->u.p.index = NULL;
        j = 2;

        /* Must have parameter type. */
        j += get_paramtype(&s[j], &node->u.p.type, fname, lnr, *i + j + 1);
        if (!node->u.p.type)
            goto param_error;

        /* Can have index specifier; N or (N1:N2). */
        if (isdigit(s[j])) {
            j += get_natural(&s[j], &tmp, fname, lnr, *i + j + 1);
            node->u.p.index = di_set(node->u.p.index, 0, tmp);
            if (di_get(node->u.p.index, 0) < 0) {
                node->u.p.len = -1;
                goto param_error;
            } else
                node->u.p.len = 1;
        } else if (s[j] == '(') {
            j += get_vector(&s[j], &node->u.p.len, &node->u.p.index, fname,
                            lnr, *i + j + 1);
            if (node->u.p.len < 0)
                goto param_error;
        } else
            node->u.p.len = 0;

        /* Can have transpose operator. */
        if (s[j] == '^' || s[j] == '"') {
            j++;
            node->u.p.column = TRUE;
        } else
            node->u.p.column = FALSE;

        /* Must end. */
        if (s[j] && s[j] == '!' && s[j + 1] == '>') {
            j += 2;
            goto end;
        }
        fprintf(stderr, PKG_NAME "%s, line %i, column %i: Expected '!>'\n",
                fname, lnr, *i + j + 1);

      param_error:
        node->u.p.len = -1;
        while (s[j] && !(s[j] == '!' && s[j + 1] == '>'))
            j++;

    } else if (isspace(s[0]))  {
        node->type = ML_WHITESPACE;
        j = 0;
        while (isspace(s[j])) {
            node->u.l.s[j] = s[j];
            j++;
        }
        node->u.l.s[j] = '\0';
        node->u.l.len = j;
    } else {
        node->type = ML_LITERAL;
        j = 0;
        while (s[j] && !(s[j] == '<' && s[j + 1] == '!')) {
            node->u.l.s[j] = s[j];
            j++;
        }
        node->u.l.s[j] = '\0';
        node->u.l.len = j;
    }

  end:
    *i += j;
    return node;
}


/* Parse 'line' from a model file, building a 'modelline_t' struct in the
 * process.  */

modelline_t *ml_parse(const char *line, const char *fname, int lnr)
{
    modelline_t *ml;
    ml_node_t *node, *end = NULL;
    int col;

    ml = malloc(sizeof(modelline_t));
    if (!ml) {
        fprintf(stderr, PKG_NAME "Can't allocate memory!\n");
        return NULL;
    }
    ml->line = NULL;
    ml->lnr = lnr;
    ml->fname = strdup(fname);

    col = 0;
    node = get_node(line, &col, fname, lnr);
    while (node) {
        if (!ml->line)
            ml->line = node;
        else
            end->next = node;
        end = node;

        node = get_node(line, &col, fname, lnr);
    }

    return ml;
}


static int nrow(const modelline_t *ml, const param_t *p)
{
    ml_node_t *node;
    int n = 1, m;

    for (node = ml->line; node; node = node->next) {
        if (node->type == ML_PARAM && node->u.p.column) {
            if (node->u.p.len == 0) {
                switch (node->u.p.type) {
                case 'R': m = da_n(p->xr); break;
                case 'I': m = da_n(p->xi); break;
                case 'O': m = da_n(p->fun); break;
                default: m = 1;
                }
            } else
                m = node->u.p.len;

            if (m > 1 && n > 1 && m != n) {
                fprintf(stderr,PKG_NAME "%s, line %i: "
                        "Columns of different length; using shortest\n",
                        ml->fname, ml->lnr);
                n = (n < m) ? n : m;
            } else
                n = m;
        }
    }

    return n;
}


static int read_param_node(ml_node_t *node, int row, const char *line,
                           param_t *p)
{
    int i, j, k, len, all;
    char *end;
    char command[MAXLINESIZE];

    assert(node->type == ML_PARAM);

    all = (node->u.p.len == 0);

    i = 0;
    switch (node->u.p.type) {
    case 'I':
        if (!node->u.p.column && (all || node->u.p.len > 1)) {

            /* Row vector. */
            len = all ? da_n(p->xi) : node->u.p.len;
            for (j = 0; j < len; j++) {
                k = all ? j : di_get(node->u.p.index, j);
                p->xi = di_set(p->xi, k, strtol(&line[i], &end, 10));

                i += end - &line[i];

                if (p->info)
                    printf("Read value %i for I%i\n", di_get(p->xi,k), k);

                if (p->usematc) {
                    sprintf(command, "I(%d) = %d", k, di_get(p->xi,k));
                    MTC_DOMATH(command);
                }
            }

        } else {

            /* Column vector or scalar.  */
            if (node->u.p.len == 1)
                k = di_get(node->u.p.index,0);
            else
                k = all ? row : di_get(node->u.p.index, row);

            p->xi = di_set(p->xi, k, strtol(&line[i], &end, 10));

            i += end - &line[i];

            if (p->info)
                printf("Read value %i for I%i\n", di_get(p->xi,k), k);

            if (p->usematc) {
                sprintf(command, "I(%d) = %d", k, di_get(p->xi,k));
                MTC_DOMATH(command);
            }
        }
        break;

    case 'R':
        if (!node->u.p.column && (all || node->u.p.len > 1)) {

            /* Row vector. */
            len = all ? da_n(p->xr) : node->u.p.len;
            for (j = 0; j < len; j++) {
                k = all ? j : di_get(node->u.p.index,j);
                p->xr = dr_set(p->xr, k, strtod(&line[i], &end));

                i += end - &line[i];

                if (p->info)
                    printf("Read value %e for R%i\n", dr_get(p->xr,k), k);

                if (p->usematc) {
                    sprintf(command, "R(%d) = %e", k, dr_get(p->xr,k));
                    MTC_DOMATH(command);
                }
            }

        } else {

            /* Column vector or scalar.  */
            if (node->u.p.len == 1)
                k = di_get(node->u.p.index,0);
            else
                k = all ? row : di_get(node->u.p.index, row);

            p->xr = dr_set(p->xr, k, strtod(&line[i], &end));

            i += end - &line[i];

            if (p->info)
                printf("Read value %e for R%i\n", dr_get(p->xr,k), k);

            if (p->usematc) {
                sprintf(command, "R(%d) = %e", k, dr_get(p->xr,k));
                MTC_DOMATH(command);
            }
        }
        break;

    case 'O':
        if (!node->u.p.column && (all || node->u.p.len > 1)) {

            /* Row vector. */
            len = all ? da_n(p->fun) : node->u.p.len;
            for (j = 0; j < len; j++) {
                k = all ? j : di_get(node->u.p.index,j);
                p->fun = dr_set(p->fun, k, strtod(&line[i], &end));

                i += end - &line[i];

                if (p->info)
                    printf("Read value %e for O%i\n", dr_get(p->fun,k), k);

                if (p->usematc) {
                    sprintf(command, "O(%d) = %e", k, dr_get(p->fun,k));
                    MTC_DOMATH(command);
                }
            }

        } else {

            /* Column vector or scalar.  */
            if (node->u.p.len == 1)
                k = di_get(node->u.p.index,0);
            else
                k = all ? row : di_get(node->u.p.index, row);

            p->fun = dr_set(p->fun, k, strtod(&line[i], &end));

            i += end - &line[i];

            if (p->info)
                printf("Read value %e for O%i\n", dr_get(p->fun,k), k);

            if (p->usematc) {
                sprintf(command, "O(%d) = %e", k, dr_get(p->fun,k));
                MTC_DOMATH(command);
            }
        }
        break;

    case 'T':
        fprintf(stderr, PKG_NAME "TAG parameters are input only\n");

    default:
        assert(0);
    }

    return i;
}


void ml_read(modelline_t *ml, FILE *fd, param_t *p)
{
    int row, i;
    size_t len;
    ml_node_t *node;
    const char *line;

    for (row = 0; row < nrow(ml, p); row++) {
        line = readline(fd, &len);
        if (!line) {
            fprintf(stderr, PKG_NAME "Premature end of input\n");
            return;
        }
        i = 0;
        for (node = ml->line; node; node = node->next) {
            switch (node->type) {
            case ML_LITERAL:
                i += node->u.l.len;
            case ML_WHITESPACE:
                while (i < len && isspace(line[i])) i++;
                break;
            case ML_PARAM:
                if (i >= len) {
                    /* TODO: More informative error message.  */
                    fprintf(stderr, PKG_NAME "Premature end of line; "
                                             "expected parameter\n");
                    break;
                }
                i += read_param_node(node, row, &line[i], p);
                break;
            default:
                assert(0);
            }
        }
    }
}


static void print_param_node(ml_node_t *node, int row, FILE *fd,
                             const param_t *p)
{
    int j, k, len;
    int all;

    assert(node->type == ML_PARAM);

    all = (node->u.p.len == 0);

    switch (node->u.p.type) {
    case 'I':
        if (!node->u.p.column && (all || node->u.p.len > 1)) {

            /* Row vector. */
            len = all ? da_n(p->xi) : node->u.p.len;
            for (j = 0; j < len; j++) {
                k = all ? j : di_get(node->u.p.index,j);
                fprintf(fd, "%d", di_get(p->xi,k));
                if (j < len-1) fputc(' ',fd);
                if (p->info)
                    printf("Wrote value %i for I%i\n", di_get(p->xi,k), k);
            }
        } else {

            /* Column vector or scalar.  */
            if (node->u.p.len == 1)
                k = di_get(node->u.p.index,0);
            else
                k = all ? row : di_get(node->u.p.index,row);

            fprintf(fd, "%d", di_get(p->xi,k));
            if (p->info)
                printf("Wrote value %i for I%i\n", di_get(p->xi,k), k);
        }
        break;

    case 'R':
        if (!node->u.p.column && (all || node->u.p.len > 1)) {

            /* Row vector. */
            len = all ? da_n(p->xr) :  node->u.p.len;
            for (j = 0; j < len; j++) {
                k = all ? j : di_get(node->u.p.index,j);
                fprintf(fd, "%.8e", dr_get(p->xr,k));
                if (j < len-1) fputc(' ',fd);
                if (p->info)
                    printf("Wrote value %e for R%i\n", dr_get(p->xr,k), k);
            }
        } else {

            /* Column vector or scalar.  */
            if (node->u.p.len == 1)
                k = di_get(node->u.p.index,0);
            else
                k = all ? row : di_get(node->u.p.index,row);

            fprintf(fd, "%.8e", dr_get(p->xr,k));
            if (p->info)
                printf("Wrote value %e for R%i\n", dr_get(p->xr,k), k);
        }
        break;

    case 'O':
        if (!node->u.p.column && (all || node->u.p.len > 1)) {

            /* Row vector. */
            len = all ? da_n(p->fun) : node->u.p.len;
            for (j = 0; j < len; j++) {
                k = all ? j : di_get(node->u.p.index,j);
                fprintf(fd, "%.8e", dr_get(p->fun,k));
                if (j < len-1) fputc(' ',fd);
                if (p->info)
                    printf("Wrote value %e for O%i\n", dr_get(p->fun,k), k);
            }
        } else {

            /* Column vector or scalar.  */
            if (node->u.p.len == 1)
                k = di_get(node->u.p.index,0);
            else
                k = all ? row : di_get(node->u.p.index,row);

            fprintf(fd, "%.8e", dr_get(p->fun,k));
            if (p->info)
                printf("Wrote value %e for O%i\n", dr_get(p->fun,k), k);
        }
        break;

    case 'T':
        fputs(p->tag, fd);
        if (p->info)
            printf("Wrote value '%s' for tag\n", p->tag);
        break;

    default:
        assert(0);
    }
}


void ml_print(modelline_t *ml, FILE *fd, const param_t *p)
{
    int row;
    ml_node_t *node;

    for (row = 0; row < nrow(ml, p); row++) {
        for (node = ml->line; node; node = node->next) {
            switch (node->type) {
            case ML_LITERAL:
            case ML_WHITESPACE:
                fputs(node->u.l.s, fd);
                break;
            case ML_PARAM:
                print_param_node(node, row, fd, p);
                break;
            default:
                assert(0);
            }
        }
    }
}


/* Free all memory used by ml.  */

void ml_kill(modelline_t *ml)
{
    ml_node_t *node, *next;

    node = ml->line;
    while (node) {
        next = node->next;
        if (node->type == ML_PARAM)
            da_kill(node->u.p.index);
        free(node);
        node = next;
    }
    free(ml->fname);
    free(ml);
}
