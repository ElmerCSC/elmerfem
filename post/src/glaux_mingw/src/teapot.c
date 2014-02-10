/* 02/21/03 Modified by Greg Chien <gchien@protodesign-inc.com>
 * - Add APIENTRY to be compatible with glaux.h distributed by Microsoft.
 */
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
#include <windows.h>
#include <GL/gl.h>
#include <gl/glaux.h>
#include "teapot.h"

#define static

long GRD;

#define TEAPOTSOLID 0
#define TEAPOTWIRE 1

static GLuint teapots[2] = {0, 0};

static float tex[2][2][2] = {{{0, 0},{1, 0}},{{0, 1},{1, 1}}};

static void solidTeapot(long grid, GLdouble scale)
{
    float p[4][4][3], q[4][4][3], r[4][4][3], s[4][4][3];
    long i, j, k, l;

    if (grid < 2) grid = 7;
    GRD = grid;
    teapots[TEAPOTSOLID] = glGenLists (1);
    glNewList(teapots[TEAPOTSOLID], GL_COMPILE);
    glPushMatrix ();
    glRotatef ((GLfloat)270.0, (GLfloat)1.0, (GLfloat)0.0, (GLfloat)0.0);
    glScalef ((GLdouble)0.5*scale, (GLdouble)0.5*scale, (GLdouble)0.5*scale);
    glTranslatef ((GLfloat)0.0, (GLfloat)0.0, (GLfloat)-1.5);
    for (i = 0; i < 10; i++) {
	for (j = 0; j < 4; j++)
	    for (k = 0; k < 4; k++) 
		for (l = 0; l < 3; l++) {
		    p[j][k][l] = cpdata[patchdata[i][j*4+k]][l];
		    q[j][k][l] = cpdata[patchdata[i][j*4+(3-k)]][l];
		    if (l == 1) q[j][k][l] *= (float)-1.0;
		    if (i < 6) {
			r[j][k][l] = cpdata[patchdata[i][j*4+(3-k)]][l];
			if (l == 0) r[j][k][l] *= (float)-1.0;
			s[j][k][l] = cpdata[patchdata[i][j*4+k]][l];
			if (l == 0) s[j][k][l] *= (float)-1.0;
			if (l == 1) s[j][k][l] *= (float)-1.0;
		    }
		}
	glMap2f(GL_MAP2_TEXTURE_COORD_2, 0, 1, 2, 2, 0, 1, 4, 2, &tex[0][0][0]);
	glMap2f(GL_MAP2_VERTEX_3, 0, 1, 3, 4, 0, 1, 12, 4, &p[0][0][0]);
	glEnable(GL_MAP2_VERTEX_3); glEnable(GL_MAP2_TEXTURE_COORD_2);
	glMapGrid2f(GRD, (GLfloat)0.0, (GLfloat)1.0, GRD, (GLfloat)0.0, (GLfloat)1.0);
	glEvalMesh2(GL_FILL, 0, GRD, 0, GRD);
	glMap2f(GL_MAP2_VERTEX_3, 0, 1, 3, 4, 0, 1, 12, 4, &q[0][0][0]);
	glEvalMesh2(GL_FILL, 0, GRD, 0, GRD);
	if (i < 6) {
	    glMap2f(GL_MAP2_VERTEX_3, 0, 1, 3, 4, 0, 1, 12, 4, &r[0][0][0]);
	    glEvalMesh2(GL_FILL, 0, GRD, 0, GRD);
	    glMap2f(GL_MAP2_VERTEX_3, 0, 1, 3, 4, 0, 1, 12, 4, &s[0][0][0]);
	    glEvalMesh2(GL_FILL, 0, GRD, 0, GRD);
	}
    }
    glDisable(GL_MAP2_VERTEX_3); glDisable(GL_MAP2_TEXTURE_COORD_2);
    glPopMatrix ();
    glEndList();
}

static void wireTeapot(long grid, GLdouble scale)
{
    float p[4][4][3], q[4][4][3], r[4][4][3], s[4][4][3];
    long i, j, k, l;
    
    if (grid < 2) grid = 7;
    GRD = grid;
    teapots[TEAPOTWIRE] = glGenLists (1);
    glNewList(teapots[TEAPOTWIRE], GL_COMPILE);
    glPushMatrix ();
    glRotatef ((GLfloat)270.0, (GLfloat)1.0, (GLfloat)0.0, (GLfloat)0.0);
    glScalef ((GLdouble)0.5*scale, (GLdouble)0.5*scale, (GLdouble)0.5*scale);
    glTranslatef ((GLfloat)0.0, (GLfloat)0.0, (GLfloat)-1.5);
    for (i = 0; i < 10; i++) {
	for (j = 0; j < 4; j++)
	    for (k = 0; k < 4; k++) 
		for (l = 0; l < 3; l++) {
		    p[j][k][l] = cpdata[patchdata[i][j*4+k]][l];
		    q[j][k][l] = cpdata[patchdata[i][j*4+(3-k)]][l];
		    if (l == 1) q[j][k][l] *= (float)-1.0;
		    if (i < 6) {
			r[j][k][l] = cpdata[patchdata[i][j*4+(3-k)]][l];
			if (l == 0) r[j][k][l] *= (float)-1.0;
			s[j][k][l] = cpdata[patchdata[i][j*4+k]][l];
			if (l == 0) s[j][k][l] *= (float)-1.0;
			if (l == 1) s[j][k][l] *= (float)-1.0;
		    }
		}
	glMap2f(GL_MAP2_TEXTURE_COORD_2, 0, 1, 2, 2, 0, 1, 4, 2, &tex[0][0][0]);
	glMap2f(GL_MAP2_VERTEX_3, 0, 1, 3, 4, 0, 1, 12, 4, &p[0][0][0]);
	glEnable(GL_MAP2_VERTEX_3); glEnable(GL_MAP2_TEXTURE_COORD_2);
	glMapGrid2f(GRD, (GLfloat)0.0, (GLfloat)1.0, GRD, (GLfloat)0.0, (GLfloat)1.0);
	glEvalMesh2(GL_LINE, 0, GRD, 0, GRD);
	glMap2f(GL_MAP2_VERTEX_3, 0, 1, 3, 4, 0, 1, 12, 4, &q[0][0][0]);
	glEvalMesh2(GL_LINE, 0, GRD, 0, GRD);
	if (i < 6) {
	    glMap2f(GL_MAP2_VERTEX_3, 0, 1, 3, 4, 0, 1, 12, 4, &r[0][0][0]);
	    glEvalMesh2(GL_LINE, 0, GRD, 0, GRD);
	    glMap2f(GL_MAP2_VERTEX_3, 0, 1, 3, 4, 0, 1, 12, 4, &s[0][0][0]);
	    glEvalMesh2(GL_LINE, 0, GRD, 0, GRD);
	}
    }
    glDisable(GL_MAP2_VERTEX_3); glDisable(GL_MAP2_TEXTURE_COORD_2);
    glPopMatrix ();
    glEndList();
}

void APIENTRY auxSolidTeapot(GLdouble scale)
{
    if (glIsList(teapots[TEAPOTSOLID]) == 0)
	solidTeapot (14, scale);
    glCallList(teapots[TEAPOTSOLID]);
}

void APIENTRY auxWireTeapot(GLdouble scale)
{
    if (glIsList(teapots[TEAPOTWIRE]) == 0)
	wireTeapot (10, scale);
    glCallList(teapots[TEAPOTWIRE]);
}
