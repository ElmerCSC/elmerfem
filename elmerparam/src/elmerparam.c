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
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>

#include "elmerparam.h"
#include "global.h"
#include "modelline.h"

#ifndef HAVE_STRNCASECMP
    int strncasecmp (const char *s1, const char *s2, size_t n);
#endif

const char instructions[] = "\n"
    "************************************************************************\n"
    "This is an utility routine for doing optimization with problems where\n"
    "the optimization function cannot be directly accessed.\n"
    "You should provide a command file with the name defined in\n"
    "ELMERPARAM_STARTINFO.\n"
    "The following example illustrates the stucture of the command file\n"
    "\n"
    "Comment = Optimzation overcoat routine example\n"
    "Parameter R1 min = 1.0\n"
    "Parameter R1 max = 5.0\n"
    "Input File = angle.grd\n"
    "Execute = ElmerGrid 1 2 angle\n"
    "Input File = TempDist.sif\n"
    "Execute = ElmerSolver\n"
    "Output File = maxtemp.dat\n"
    "\n"
    "For each input and output file there must a template with the suffix\n"
    ".model. Optionally, the template file can be given by adding\n"
    "'USING templatefile' to the 'Input file' and 'Output File' commands.\n"
    "In the file the parameters to be optimized are marked with special\n"
    "brackets For example, <!R3!> and <!I4!>\n\n";


static void param_init(param_t *p)
{
    p->xr = p->fun = NULL;
    p->xi = NULL;
    /* p->fun = DEFAULT_FUN; */
    p->info = TRUE;
    p->usematc = FALSE;
    p->taglen = 0;
    p->tag[0] = '\0';
    p->lnr = 0;
}

/* Replace <!T!> with 'tag' everywhere in 's'.  */

static void expand_tag(const param_t *p, char *s)
{
    char t[MAXLINESIZE];
    int i, j;

    if (!strstr(s, "<!t!>") && !strstr(s, "<!T!>"))
        return;

    /* Expand 's' into 't'.  */
    i = j = 0;
    while (s[i]) {
        if (strncasecmp(&s[i], "<!T!>", 5) == 0) {
            strcpy(&t[j], p->tag);
            j += p->taglen;
            i += 5;
        } else
            t[j++] = s[i++];

        assert(i < MAXLINESIZE);
        assert(j < MAXLINESIZE);
    }
    t[j] = '\0';

    /* Copy the expanded string back to 's'.  */
    strcpy(s, t);
}


/* Return TRUE if c is a comment marker.  */

static int iscomment(int c)
{
    return c == '*' || c == '#' || c == '!';
}


/* Get a "logical" line; if a line ends with '\', join with next line.  The
 * result will have trailing white spaces trimmed.  Ignore empty lines and
 * comment lines.  Return number of characters read (0 in case of EOF).  */

static int get_line(param_t *p, char *line, FILE *fd)
{
    int c;
    int n;
    int empty = TRUE;

  emptyline:
    n = 0;
    while ((c = fgetc(fd)) != EOF && c != '\n') {
        /* Ignore comments.  */
        if (empty && iscomment(c))
            while ((c = fgetc(fd)) != '\n');

        line[n++] = c;
        if (!isspace(c))
            empty = FALSE;
    }

    if (c == '\n')
        p->lnr++;

    if (empty) {
        if (c != EOF)
            goto emptyline;
        else
            return 0;
    }

    /* Trim spaces from the end to see if last non-space is '\'.  In case it is,
     * join this line with the next.  */
    assert(n > 0);
    n--;
    while (n >= 0 && isspace(line[n]))
        n--;
    if (line[n] == '\\')
        n += get_line(p, &line[n], fd);
    else {
        line[n + 1] = '\0';
        n++;
    }

    return n;
}


/* Get next command with parameters from 'io', both with leading and trailing
 * white spaces trimmed. 'command' will be converted to all uppercase
 * characters. 'command' and 'params' are supposed to be allocated to
 * be long enough.  */

static int get_command(param_t *p, char *command, char *params, FILE *io)
{
    int i, j;
    char line[MAXLINESIZE];

    if (!get_line(p, line, io))
        return FALSE;

    i = 0;
    while (isspace(line[i]))
        i++;

    if (line[i] == '$')
        strcpy(command, "$");
    else {
        j = 0;

        while (line[i] != '=' && line[i] != '\0')
            command[j++] = toupper(line[i++]);
        j--;

        while (isspace(command[j]))
            j--;

        command[j + 1] = '\0';
    }

    i++;
    while (isspace(line[i]))
        i++;

    strcpy(params, &line[i]);

    return TRUE;
}


/* Create input file 'fname' using template file 'mname'.  */

static void create_input(const param_t *p, const char *fname,
                         const char *mname)
{
    int lnr;
    FILE *file, *model;
    modelline_t *ml = NULL;
    char line[MAXLINESIZE], *input;

    if ((model = fopen(mname, "r")) == NULL) {
        fprintf(stderr,PKG_NAME "Can't open template file %s for reading: %s\n",
                mname, strerror(errno));
        return;
    }

    if ((file = fopen(fname, "w")) == NULL) {
        fclose(model);
        fprintf(stderr, PKG_NAME "Can't open file %s for writing: %s\n",
                fname, strerror(errno));
        return;
    }

    if (p->info)
        printf("Creating input file %s using template %s\n", fname, mname);

    input = fgets(line, MAXLINESIZE, model);
    lnr = 1;
    while (input) {
        ml = ml_parse(line, mname, lnr);
        ml_print(ml, file, p);
        ml_kill(ml);

        input = fgets(line, MAXLINESIZE, model);
        lnr++;
    }

    fclose(file);
    fclose(model);
}


/* Extract values from 'fname', using template file 'mname'.  */

static void read_output(param_t *p, const char *fname, const char *mname)
{
    char *input;
    char line[MAXLINESIZE];
    FILE *file, *model;
    int lnr;
    modelline_t *ml;


    if ((model = fopen(mname, "r")) == NULL) {
        fprintf(stderr,PKG_NAME "Can't open template file %s for reading: %s\n",
                mname, strerror(errno));
        return;
    }

    if ((file = fopen(fname, "r")) == NULL) {
        fclose(model);
        fprintf(stderr, PKG_NAME 
                "Can't open outputfile file %s for reading: %s\n",
                fname, strerror(errno));
        return;
    }

    if (p->info)
        printf("Reading from output file %s using template %s\n", fname, mname);

    input = fgets(line, MAXLINESIZE, model);
    lnr = 1;
    while (input) {
        ml = ml_parse(line, mname, lnr);
        ml_read(ml, file, p);
        ml_kill(ml);

        input = fgets(line, MAXLINESIZE, model);
        lnr++;
    }
    fclose(file);
    fclose(model);
}


static void save_result_line(const param_t *p, const char *fname)
{
    int i;
    FILE *file;

    if (p->info)
        printf("Writing result line to output file %s\n", fname);

    if ((file = fopen(fname, "a")) == NULL) {
        fprintf(stderr,PKG_NAME "Can't open save file %s for appending: %s\n",
                fname, strerror(errno));
        return;
    }

    if (da_n(p->xi) > 0) {
        for (i = 0; i < da_n(p->xi); i++)
            fprintf(file, "%d ", di_get(p->xi,i));
    }
    if (da_n(p->xr) > 0) {
        for (i = 0; i < da_n(p->xr); i++)
            fprintf(file, "%.6e ", dr_get(p->xr,i));
    }
    if (da_n(p->fun) > 0) {
        for (i = 0; i < da_n(p->fun); i++)
            fprintf(file, "%.6e ", dr_get(p->fun,i));
    }

    fputc('\n', file);

    fclose(file);
}


/* Extract file name and model name from 'str', and store in correspondingly
 * named arguments, which are supposed to be long enough.  The format of 'str'
 * is supposed to be either
 *
 *    filename
 *
 * in which case model name becomes 'filename.model', or
 *
 *    filename USING modelname
 *
 * Return TRUE for success and FALSE for failure.
 */

static int file_and_modelname(const param_t *p, const char *str,
                              char *fname, char *mname)
{
    const char *s;
    int strl, n;

    strl = strlen(str);

    s = str;
    while (*s && isspace(*s))
        s++;

    /* Get filename.  */
    n = 0;
    while (*s && !isspace(*s)) {
        fname[n++] = *s;
        s++;
    }
    fname[n] = '\0';

    while (*s && isspace(*s))
        s++;
    if (!*s) {

        /* End of line; create a model name based on the file name.  */
        sprintf(mname, "%s.model", fname);

    } else if (strncasecmp(s, "USING", 5) == 0) {
        s += 5;

        /* A space must follow 'USING'. */
        if (!isspace(*s)) {
            fprintf(stderr, PKG_NAME "%s, line %i: "
                    "Expected 'USING mname', found '%s'\n", p->cmdfile,
                    p->lnr, s);
            return FALSE;
        }

        while (*s && isspace(*s))
            s++;

        /* There must be something (other than spaces) after 'USING'.  */
        if (!*s) {
            fprintf(stderr,
                    PKG_NAME
                    "%s, line %i: Expected modelname after USING\n",
                    p->cmdfile, p->lnr);
            return FALSE;
        }

        /* Copy model name.  */
        n = 0;
        while (*s && !isspace(*s)) {
            mname[n++] = *s;
            s++;
        }
        mname[n] = '\0';

    } else {

        /* Something else than 'USING' after the file name.  */
        fprintf(stderr, PKG_NAME "%s, line %i: "
                "Expected 'USING modelname', found '%s'\n", p->cmdfile,
                p->lnr, s);
        return FALSE;
    }

    return TRUE;
}


static void generic_function(param_t *p)
{
    char fname[MAXFILESIZE], mname[MAXFILESIZE];
    char command[MAXLINESIZE], params[MAXLINESIZE];
    char *strpntr, *matcpntr, *compntr;
    FILE *in;
    int i, xilim, minmax, isint, len;
    double xrlim, tmp;
    int last_was_matc = FALSE;

    if ((in = fopen("ELMERPARAM_STARTINFO", "r")) == NULL) {
        fprintf(stderr,
                PKG_NAME
                "Can't open file ELMERPARAM_STARTINFO for reading: %s\n",
                strerror(errno));
        fputs(instructions, stderr);
        return;
    }
    fgets(params, MAXLINESIZE, in);
    sscanf(params, "%s", p->cmdfile);
    fclose(in);

    if ((in = fopen(p->cmdfile, "r")) == NULL) {
        fprintf(stderr, PKG_NAME "Can't open command file %s for reading: %s\n",
                p->cmdfile, strerror(errno));
        fputs(instructions, stderr);
        return;
    }

    if (p->info) printf("Loading commands from file '%s'.\n", p->cmdfile);

    for (;;) {

        if (!get_command(p, command, params, in)) {
            if (p->info)
                printf("Reached the end of command file\n");
            break;
        }

        expand_tag(p, params);

        /* Give special treatment to MATC code.  */
        if (command[0] == '$' && p->usematc) {
            MTC_DOMATH(params);
            last_was_matc = TRUE;
            continue;
        }

        /* Update local representations of I, R and O (might have been modified
         * by MATC).  */
        if (last_was_matc) {
            p->xi = di_set_from_matc(p->xi, "I");
            p->xr = dr_set_from_matc(p->xr, "R");
            p->fun = dr_set_from_matc(p->fun, "O");
        }

        if (strcmp(command, "COMMENT") == 0)
            printf("***\n*** %s\n***\n", params);

        else if (strcmp(command, "MATC") == 0) {
            if (params[0] == 'T' || params[0] == 't') {
                MTC_INIT(p);
            } else
                p->usematc = FALSE;
        }

        else if (strcmp(command, "ECHO") == 0) {
            if (params[0] == 'F' || params[0] == 'f')
                p->info = FALSE;
            else
                p->info = TRUE;
        }

        else if (strcmp(command, "INPUT FILE") == 0) {
            if (strlen(params) == 0) {
                fprintf(stderr, PKG_NAME "%s, line %i: "
                        "Command 'INPUT FILE' needs argument\n",
                        p->cmdfile, p->lnr);
                return;
            }

            if (!file_and_modelname(p, params, fname, mname))
                return;

            create_input(p, fname, mname);
        }

        else if (strcmp(command, "OUTPUT FILE") == 0) {
            if (!file_and_modelname(p, params, fname, mname))
                return;
            read_output(p, fname, mname);
        }

        else if (strcmp(command, "EXECUTE") == 0) {
            if (p->info)
                printf("Executing command '%s'\n", params);
            system(params);
        }

        else if (strcmp(command, "COST FUNCTION") == 0
                 || strcmp(command, "FUNCTION") == 0) {
            if (p->usematc) {
                if ((matcpntr = strstr(params, "$"))) {
                    matcpntr = MTC_DOMATH(&matcpntr[1]);
                    strcpy(params, matcpntr);
                }
            }
            sscanf(params, "%le", &tmp);
            p->fun = dr_set(p->fun,0,tmp);
        }

        else if (strncmp(command, "PARAMETER", 9) == 0) {
            if ((strpntr = strstr(command, "PARAMETER I"))) {
                compntr = strpntr;
                isint = TRUE;
                minmax = 0;
            } else if ((strpntr = strstr(command, "PARAMETER MIN I"))) {
                compntr = strpntr;
                isint = TRUE;
                minmax = -1;
            } else if ((strpntr = strstr(command, "PARAMETER MAX I"))) {
                compntr = strpntr;
                isint = TRUE;
                minmax = 1;
            } else if ((strpntr = strstr(command, "PARAMETER R"))) {
                compntr = strpntr;
                isint = FALSE;
                minmax = 0;
            } else if ((strpntr = strstr(command, "PARAMETER MIN R"))) {
                compntr = strpntr;
                isint = FALSE;
                minmax = -1;
            } else if ((strpntr = strstr(command, "PARAMETER MAX R"))) {
                compntr = strpntr;
                isint = FALSE;
                minmax = 1;
            } else
                continue;

            len = 11 + 4 * abs(minmax);

            if (p->usematc) {
                if ((matcpntr = strstr(params, "$"))) {
                    matcpntr = MTC_DOMATH(&matcpntr[1]);
                    strcpy(params, matcpntr);
                }
            }
            sscanf(&compntr[len], "%d", &i);

            if (isint) {
                sscanf(params, "%d", &xilim);
                if (minmax == 0) {
                    p->xi = di_set(p->xi, i, xilim);
                    if (p->usematc) {
                        sprintf(command, "I(%d) = %d", i, xilim);
                        MTC_DOMATH(command);
                    }
                } else if (minmax * di_get(p->xi,i) > minmax * xilim)
                    goto done;
            } else {
                sscanf(&strpntr[len], "%d", &i);
                sscanf(params, "%le", &xrlim);
                if (minmax == 0) {
                    p->xr = dr_set(p->xr, i, xrlim);
                    if (p->usematc) {
                        sprintf(command, "R(%d) = %e", i, xrlim);
                        MTC_DOMATH(command);
                    }
                } else if (minmax * dr_get(p->xr,i) > minmax * xrlim)
                    goto done;
            }
        } /* PARAMETER */
        else if (strcmp(command, "SAVE FILE") == 0) {
            sscanf(params, "%s", fname);
            save_result_line(p, fname);
        }
        else {
            fprintf(stderr, PKG_NAME "%s, line %i: Unknown command '%s'\n",
                    p->cmdfile, p->lnr, command);
            return;
        }

        last_was_matc = FALSE;
    }

  done:

    /* Update local representations of I, R and O (might have been modified
     * by MATC).  */
    if (last_was_matc) {
        p->xi = di_set_from_matc(p->xi, "I");
        p->xr = dr_set_from_matc(p->xr, "R");
        p->fun = dr_set_from_matc(p->fun, "O");
    }

    fclose(in);

    if (p->info) {
        if (da_n(p->fun) > 0) {
            printf("Function Result = [");
            for (i = 0; i < da_n(p->fun); i++)
                printf(" %.6e ", dr_get(p->fun,i));
            printf("]\n");
        }
        if (da_n(p->xi) > 0) {
            printf("Integer parameters = [");
            for (i = 0; i < da_n(p->xi); i++)
                printf(" %d ", di_get(p->xi,i));
            printf("]\n");
        }
        if (da_n(p->xr) > 0) {
            printf("Real parameters = [");
            for (i = 0; i < da_n(p->xr); i++)
                printf(" %.6e ", dr_get(p->xr,i));
            printf("]\n");
        }
    }

    return;
}


/* C interface */

void elmer_param_vec(int nfun, double *f, int nr, const double *xr, int ni,
                     const int *xi, const char *tag)
{
    param_t p;
    int i;

    param_init(&p);

    if (tag) {
        p.taglen = strlen(tag);
        assert(p.taglen < MAXLINESIZE);
        strcpy(p.tag, tag);
    }

    if (nr > 0) {
        assert(xr);

        for (i = 0; i < nr; i++)
            p.xr = dr_set(p.xr, i, xr[i]);
    }

    if (ni > 0) {
        assert(xi);

        for (i = 0; i < ni; i++)
            p.xi = di_set(p.xi, i, xi[i]);
    }

    /* We must make it look like there are nfun elements in p.fun, even though
       we have no actual values yet (so that da_n(p.fun) == nfun).  */
    assert(nfun >= 1);
    for (i = 0; i < nfun; i++)
        p.fun = dr_set(p.fun, i, DEFAULT_FUN);

    generic_function(&p);

    for (i = 0; i < nfun; i++)
        f[i] = dr_get(p.fun, i);

    da_kill(p.xi);
    da_kill(p.xr);
    da_kill(p.fun);
}


double elmer_param(int nr, const double *xr,
                   int ni, const int *xi, const char *tag)
{
    double f;

    elmer_param_vec(1, &f, nr, xr, ni, xi, tag);
    return f;
}



/* Functions for Fortran calls */

double FC_FUNC_(elmer_param_c,ELMER_PARAM_C)(const int *nr, const double *xr,
                      const int *ni, const int *xi,
                      const int *taglen, const char *tag)
{
    char ctag[MAXLINESIZE];

    if (*taglen > 0) {
        /* Copy tag to a C style '\0'-terminated string.  */
        assert(*taglen < MAXLINESIZE - 1);
        strncpy(ctag, tag, *taglen);
        ctag[*taglen] = '\0';

        return elmer_param(*nr, xr, *ni, xi, ctag);
    } else
        return elmer_param(*nr, xr, *ni, xi, NULL);
}


void FC_FUNC_(elmer_param_vec_c,ELMER_PARAM_VEC_C)(const int *nfun, double *fun,
                          const int *nr, const double *xr,
                          const int *ni, const int *xi,
                          const int *taglen, const char *tag)
{
    char ctag[MAXLINESIZE];

    if (*taglen > 0) {
        /* Copy tag to a C style '\0'-terminated string.  */
        assert(*taglen < MAXLINESIZE - 1);
        strncpy(ctag, tag, *taglen);
        ctag[*taglen] = '\0';

        elmer_param_vec(*nfun, fun, *nr, xr, *ni, xi, ctag);
    } else
        elmer_param_vec(*nfun, fun, *nr, xr, *ni, xi, NULL);
}
