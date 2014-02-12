/*
 * (c) Copyright 1993, Silicon Graphics, Inc.
 * ALL RIGHTS RESERVED
 * Permission to use, copy, modify, and distribute this software for
 * any purpose and without fee is hereby granted, provided that the above
 * copyright notice appear in all copies and that both the copyright notice
 * and this permission notice appear in supporting documentation, and that
 * the name of Silicon Graphics, Inc. not be used in advertising
 * or publicity pertaining to distribution of the software without specific,
 * written prior permission.
 *
 * THE MATERIAL EMBODIED ON THIS SOFTWARE IS PROVIDED TO YOU "AS-IS"
 * AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE,
 * INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR
 * FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL SILICON
 * GRAPHICS, INC.  BE LIABLE TO YOU OR ANYONE ELSE FOR ANY DIRECT,
 * SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY
 * KIND, OR ANY DAMAGES WHATSOEVER, INCLUDING WITHOUT LIMITATION,
 * LOSS OF PROFIT, LOSS OF USE, SAVINGS OR REVENUE, OR THE CLAIMS OF
 * THIRD PARTIES, WHETHER OR NOT SILICON GRAPHICS, INC.  HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH LOSS, HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE
 * POSSESSION, USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 * US Government Users Restricted Rights
 * Use, duplication, or disclosure by the Government is subject to
 * restrictions set forth in FAR 52.227.19(c)(2) or subparagraph
 * (c)(1)(ii) of the Rights in Technical Data and Computer Software
 * clause at DFARS 252.227-7013 and/or in similar or successor
 * clauses in the FAR or the DOD or NASA FAR Supplement.
 * Unpublished-- rights reserved under the copyright laws of the
 * United States.  Contractor/manufacturer is Silicon Graphics,
 * Inc., 2011 N.  Shoreline Blvd., Mountain View, CA 94039-7311.
 *
 * OpenGL(TM) is a trademark of Silicon Graphics, Inc.
 */
/* Routines to manipulate 3 dimensional vectors.  All these routines
 * should work even if the input and output vectors are the same.
 */

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/gl.h>
#include "3d.h"
#include "tk.h"

#define static

void (*errfunc)(char *) = 0;

void seterrorfunc(void (*func)(char *))
{
    errfunc = func;
}

void aux_error(char *s)
{
    if (errfunc)
        (*errfunc)(s);
    else {
        MESSAGEBOX(GetFocus(), s, "Error", MB_OK);
        exit(1);
    }
}

void diff3(GLdouble p[3], GLdouble q[3], GLdouble diff[3])
{
    diff[0] = p[0] - q[0];
    diff[1] = p[1] - q[1];
    diff[2] = p[2] - q[2];
}

void add3(GLdouble p[3], GLdouble q[3], GLdouble sum[3])
{
    sum[0] = p[0] + q[0];
    sum[1] = p[1] + q[1];
    sum[2] = p[2] + q[2];
}

void scalarmult(GLdouble s, GLdouble v[3], GLdouble vout[3])
{
    vout[0] = v[0]*s;
    vout[1] = v[1]*s;
    vout[2] = v[2]*s;
}

GLdouble dot3(GLdouble p[3], GLdouble q[3])
{
    return p[0]*q[0] + p[1]*q[1] + p[2]*q[2];
}

GLdouble length3(GLdouble v[3])
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

GLdouble dist3(GLdouble p[3], GLdouble q[3])
{
    GLdouble d[3];

    diff3(p, q, d);
    return length3(d);
}

void copy3(GLdouble old[3], GLdouble new[3])
{
    new[0] = old[0], new[1] = old[1], new[2] = old[2];
}

void crossprod(GLdouble v1[3], GLdouble v2[3], GLdouble prod[3])
{
    GLdouble p[3];      /* in case prod == v1 or v2 */

    p[0] = v1[1]*v2[2] - v2[1]*v1[2];
    p[1] = v1[2]*v2[0] - v2[2]*v1[0];
    p[2] = v1[0]*v2[1] - v2[0]*v1[1];
    prod[0] = p[0]; prod[1] = p[1]; prod[2] = p[2];
}

void normalize(GLdouble v[3])
{
    GLdouble d;

    d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    if (d == (GLdouble)0.0) {
        aux_error("normalize: zero length vector");
        v[0] = d = (GLdouble)1.0;
    }
    d = (GLdouble)1/d;
    v[0] *= d; v[1] *= d; v[2] *= d;
}

void identifymat3(GLdouble m[3][3])
{
    int i, j;

    for (i=0; i<3; i++)
        for (j=0; j<3; j++)
            m[i][j] = (i == j) ? (GLdouble)1.0 : (GLdouble)0.0;
}

void copymat3(GLdouble *to, GLdouble *from)
{
    int i;

    for (i=0; i<9; i++) {
        *to++ = *from++;
    }
}

void xformvec3(GLdouble v[3], GLdouble m[3][3], GLdouble vm[3])
{
    GLdouble result[3]; /* in case v == vm */
    int i;

    for (i=0; i<3; i++) {
        result[i] = v[0]*m[0][i] + v[1]*m[1][i] + v[2]*m[2][i];
    }
    for (i=0; i<3; i++) {
        vm[i] = result[i];
    }
}

long samepoint(GLdouble p1[3], GLdouble p2[3])
{
    if (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2])
        return 1;
    return 0;
}

void perpnorm(GLdouble p1[3], GLdouble p2[3], GLdouble p3[3], GLdouble n[3])
{
    GLdouble d1[3], d2[3];

    diff3(p2, p1, d1);
    diff3(p2, p3, d2);
    crossprod(d1, d2, n);
    normalize(n);
}
