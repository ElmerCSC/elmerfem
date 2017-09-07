/*
 * Copyright (c) December 5 2006 - , CSC - IT Center for Science Ltd., Finland
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program (in file fem/GPL-2); if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

/*
 * Bunch of Fortran callable functions to read/write binary data. They can be
 * called from Fortran directly, but there are some wrapper procedures in
 * biniomod.f90 that are generally more programmer friendly and has
 * explicit interfaces.
 */

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <sys/types.h>
#include "config.h"

#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#else
/* In lack of inttypes.h, fall back to using the good old unspecified-size types
 * and hope that they all exist and are of the correct sizes.  */
typedef unsigned char      uint8_t;
typedef unsigned short     uint16_t;
typedef unsigned int       uint32_t;
typedef unsigned long long uint64_t;
#endif

#define MAX_UNITS 1000

struct unit_t {
    int convert;  /* Convert endianness on reading.  */
    FILE *fd;
};

static struct unit_t units[MAX_UNITS] = {{ 0 }};


/* Return 'L' if host is little endian, 'B' if big endian.  Or perhaps I got it
 * wrong and it's the other way aound -- doesn't matter as long as we are self
 * consistent.  */

static char endianess()
{
    const uint16_t n = 0x1122;
    uint8_t *b;

    b = (uint8_t *)&n;
    if (b[0] == 0x11)
        return 'B';
    else
        return 'L';
}


void FC_FUNC(binendianess_c,BINENDIANESS_C)(char *e)
{
    *e = endianess();
}


/* Swap the bytes in an n byte object.  */

static void swap_bytes(void *o, size_t n)
{
    uint8_t *p, tmp;
    int i;

    p = o;
    for (i = 0; i < n/2; i++) {
        tmp = p[i];
        p[i] = p[n-1-i];
        p[n-1-i] = tmp;
    }
}


void FC_FUNC(binsetinputendianess_c,BINSETINPUTENDIANESS_C)(const int *unit,
                                                        const char *e )
{
    assert(units[*unit].fd);
    units[*unit].convert = (*e == endianess()) ? 0 : 1;
}


void FC_FUNC(binopen_c,BINOPEN_C)(const int *unit, const char *file,
                                 const int *file_len,
                                 const char *action,
                                 int *status)
{
    char *fname;
    char *a, b;

    fname = malloc(*file_len + 1);
    strncpy(fname, file, *file_len);
    fname[*file_len] = '\0';

    if (action[0] == 'w' || action[0] == 'W')
        a = "wb";
    else if (action[0] == 'a' || action[0] == 'A')
        a = "ab";
    else
        a = "rb";

    assert(!units[*unit].fd);

    units[*unit].fd = fopen(fname, a);
    if (!units[*unit].fd) {
        *status = errno;
        return;
    }

    *status = 0;
}


void FC_FUNC(binclose_c,BINCLOSE_C)(const int *unit, int *status)
{
    int s;

    assert(units[*unit].fd);

    s = fclose(units[*unit].fd);
    units[*unit].fd = NULL;
    *status = (s != 0) ? errno : 0;
}


void FC_FUNC(binwriteint4_c,BINWRITEINT4_C)(const int *unit, const uint32_t *n,
                                        int *status)
{
    size_t s;

    assert(units[*unit].fd);
    s = fwrite(n, 1, 4, units[*unit].fd);
    *status = (s == 4) ? 0 : errno;
}


void FC_FUNC(binwriteint8_c,BINWRITEINT8_C)(const int *unit, const uint64_t *n,
                                        int *status)
{
    size_t s;

    assert(units[*unit].fd);
    s = fwrite(n, 1, 8, units[*unit].fd);
    *status = (s == 8) ? 0 : errno;
}


void FC_FUNC(binreadint8_c,BINREADINT8_C)(const int *unit, uint64_t *n, int *status)
{
    size_t r;

    assert(units[*unit].fd);

    r = fread(n, 1, 8, units[*unit].fd);
    if (r != 8)
        *status = (feof(units[*unit].fd)) ? -1 : errno;
    else
        *status = 0;

    if (units[*unit].convert) swap_bytes(n, 8);
}


void FC_FUNC(binwritedouble_c,BINWRITEDOUBLE_C)(const int *unit, const double *a,
                                            int *status)
{
    size_t s;

    assert(units[*unit].fd);
    s = fwrite(a, 1, 8, units[*unit].fd);
    *status = (s == 8) ? 0 : errno;
}


void FC_FUNC(binreadint4_c,BINREADINT4_C)(const int *unit, uint32_t *n, int *status)
{
    size_t r;

    assert(units[*unit].fd);

    r = fread(n, 1, sizeof(uint32_t), units[*unit].fd);
    if (r != 4)
        *status = (feof(units[*unit].fd)) ? -1 : errno;
    else
        *status = 0;

    if (units[*unit].convert) swap_bytes(n, 4);
}


void FC_FUNC(binreaddouble_c,BINREADDOUBLE_C)(const int *unit, double *a,
                                          int *status)
{
    FC_FUNC(binreadint8_c,BINREADINT8_C)(unit, (uint64_t *)a, status);
}


void FC_FUNC(binwritestring_c,BINWRITESTRING_C)(const int *unit,
                                               const char *s,
                                               const int *s_len,
                                               int *status)
{
    assert(units[*unit].fd);

    if (fwrite(s, 1, *s_len, units[*unit].fd) == *s_len
        && fputc('\0', units[*unit].fd) == '\0')

        *status = 0;
    else
        *status = errno;
}


void FC_FUNC(binreadstring_c,BINREADSTRING_C)(const int *unit,char *s,
                                             const int *s_len, int *status)
{
    int i, c;

    assert(units[*unit].fd);

    i = 0;
    c = '\0';  /* Just to avoid compiler warnings.  */
    while (i < *s_len && (c = fgetc(units[*unit].fd)) != '\0' && c != EOF)
        s[i++] = c;
    while (i < *s_len) s[i++] = ' ';

    if (c == EOF)
        *status = (ferror(units[*unit].fd)) ? errno : -1;
    else
        *status = 0;
}


void FC_FUNC(binwritechar_c,BINWRITECHAR_C)(const int *unit,char *c,
                                           int *status)
{
    assert(units[*unit].fd);
    if (fwrite(c, 1, 1, units[*unit].fd) == 1)
        *status = 0;
    else
        *status = errno;
}


#ifdef HAVE_FTELLO
off_t FC_FUNC(binftell_c,BINFTELL_C)(const int *unit)
{
    assert(units[*unit].fd);
    return ftello(units[*unit].fd);
}
#else
long FC_FUNC(binftell_c,BINFTELL_C)(const int *unit)
{
    assert(units[*unit].fd);
    return ftell(units[*unit].fd);
}
#endif


#ifdef HAVE_FSEEKO
void FC_FUNC(binfseek_c,BINFSEEK_C)(const int *unit, const off_t *offset,
                                const int *whence)
{
    assert(units[*unit].fd);

    /* TODO: Check the return value of fseek and take appropriate action in
     * case of errors.  */
    switch (*whence) {
    case 0: fseeko(units[*unit].fd, *offset, SEEK_SET); break;
    case 1: fseeko(units[*unit].fd, *offset, SEEK_CUR); break;
    case 2: fseeko(units[*unit].fd, *offset, SEEK_END); break;
    }
}
#else
void FC_FUNC(binfseek_c,BINFSEEK_C)(const int *unit, const long *offset,
                                const int *whence)
{
    assert(units[*unit].fd);

    /* TODO: Check the return value of fseek and take appropriate action in
     * case of errors.  */
    switch (*whence) {
    case 0: fseek(units[*unit].fd, *offset, SEEK_SET); break;
    case 1: fseek(units[*unit].fd, *offset, SEEK_CUR); break;
    case 2: fseek(units[*unit].fd, *offset, SEEK_END); break;
    }
}
#endif


void FC_FUNC(strerrorf_c,STRERRORF_C)(const int *e, char *s,
                                     const int *s_len)
{
    char *t;
    int i;

    t = strerror(*e);
    for (i = 0; i < *s_len && t[i]; i++)
        s[i] = t[i];
    while(i < *s_len) s[i++] = ' ';
}
