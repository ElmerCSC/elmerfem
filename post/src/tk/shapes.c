#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tk.h"
#include "private.h"

/******************************************************************************/

#define PI 3.14159265358979323846

/******************************************************************************/

void tkWireSphere(GLuint base, float radius)
{
    GLUquadricObj *quadObj;

    glNewList(base, GL_COMPILE_AND_EXECUTE);
	quadObj = gluNewQuadric();
	gluQuadricDrawStyle(quadObj, GLU_LINE);
	gluSphere(quadObj, radius, 16, 16);
    glEndList();
}

/******************************************************************************/

void tkSolidSphere(GLuint base, float radius)
{
    GLUquadricObj *quadObj;

    glNewList(base, GL_COMPILE_AND_EXECUTE);
	quadObj = gluNewQuadric();
	gluQuadricDrawStyle(quadObj, GLU_FILL);
	gluQuadricNormals(quadObj, GLU_SMOOTH);
	gluSphere(quadObj, radius, 16, 16);
    glEndList();
}

/******************************************************************************/

void tkWireCube(GLuint base, float size)
{
    static float n[6][3] = {
	{-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0},
	{0.0, -1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}
    };
    static GLint faces[6][4] = {
	{0, 1, 2, 3}, {3, 2, 6, 7}, {7, 6, 5, 4},
	{4, 5, 1, 0}, {5, 6, 2, 1}, {7, 4, 0, 3}
    };
    float x0, x1, y0, y1, z0, z1, tmp;
    float v[8][3];
    int i;

    x0 = -size / 2.0;
    x1 = size / 2.0;
    y0 = -size / 2.0;
    y1 = size / 2.0;
    z0 = -size / 2.0;
    z1 = size / 2.0;

    if (x0 > x1) {
	tmp = x0; x0 = x1; x1 = tmp;
    }
    if (y0 > y1) {
	tmp = y0; y0 = y1; y1 = tmp; 
    }
    if (z0 > z1) {
	tmp = z0; z0 = z1; z1 = tmp; 
    }
    v[0][0] = v[1][0] = v[2][0] = v[3][0] = x0;
    v[4][0] = v[5][0] = v[6][0] = v[7][0] = x1;
    v[0][1] = v[1][1] = v[4][1] = v[5][1] = y0;
    v[2][1] = v[3][1] = v[6][1] = v[7][1] = y1;
    v[0][2] = v[3][2] = v[4][2] = v[7][2] = z0;
    v[1][2] = v[2][2] = v[5][2] = v[6][2] = z1;

    glNewList(base, GL_COMPILE_AND_EXECUTE);
	for (i = 0; i < 6; i++) {
	    glBegin(GL_LINE_LOOP);
		glNormal3fv(&n[i][0]);
		glVertex3fv(&v[faces[i][0]][0]);
		glNormal3fv(&n[i][0]);
		glVertex3fv(&v[faces[i][1]][0]);
		glNormal3fv(&n[i][0]);
		glVertex3fv(&v[faces[i][2]][0]);
		glNormal3fv(&n[i][0]);
		glVertex3fv(&v[faces[i][3]][0]);
	    glEnd();
	}
    glEndList();
}

/******************************************************************************/

void tkSolidCube(GLuint base, float size)
{
    static float n[6][3] = {
	{-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0},
	{0.0, -1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}
    };
    static GLint faces[6][4] = {
	{0, 1, 2, 3}, {3, 2, 6, 7}, {7, 6, 5, 4},
	{4, 5, 1, 0}, {5, 6, 2, 1}, {7, 4, 0, 3}
    };
    float x0, x1, y0, y1, z0, z1, tmp;
    float v[8][3];
    int i;

    x0 = -size / 2.0;
    x1 = size / 2.0;
    y0 = -size / 2.0;
    y1 = size / 2.0;
    z0 = -size / 2.0;
    z1 = size / 2.0;

    if (x0 > x1) {
	tmp = x0; x0 = x1; x1 = tmp;
    }
    if (y0 > y1) {
	tmp = y0; y0 = y1; y1 = tmp; 
    }
    if (z0 > z1) {
	tmp = z0; z0 = z1; z1 = tmp; 
    }
    v[0][0] = v[1][0] = v[2][0] = v[3][0] = x0;
    v[4][0] = v[5][0] = v[6][0] = v[7][0] = x1;
    v[0][1] = v[1][1] = v[4][1] = v[5][1] = y0;
    v[2][1] = v[3][1] = v[6][1] = v[7][1] = y1;
    v[0][2] = v[3][2] = v[4][2] = v[7][2] = z0;
    v[1][2] = v[2][2] = v[5][2] = v[6][2] = z1;

    glNewList(base, GL_COMPILE_AND_EXECUTE);
	for (i = 0; i < 6; i++) {
	    glBegin(GL_POLYGON);
		glNormal3fv(&n[i][0]);
		glVertex3fv(&v[faces[i][0]][0]);
		glNormal3fv(&n[i][0]);
		glVertex3fv(&v[faces[i][1]][0]);
		glNormal3fv(&n[i][0]);
		glVertex3fv(&v[faces[i][2]][0]);
		glNormal3fv(&n[i][0]);
		glVertex3fv(&v[faces[i][3]][0]);
	    glEnd();
	}
    glEndList();
}

/******************************************************************************/

void tkWireBox(GLuint base, float width, float height, float depth)
{
    static float n[6][3] = {
	{-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0},
	{0.0, -1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}
    };
    static GLint faces[6][4] = {
	{0, 1, 2, 3}, {3, 2, 6, 7}, {7, 6, 5, 4},
	{4, 5, 1, 0}, {5, 6, 2, 1}, {7, 4, 0, 3}
    };
    float x0, x1, y0, y1, z0, z1, tmp;
    float v[8][3];
    int i;

    x0 = -width / 2.0;
    x1 = width / 2.0;
    y0 = -height / 2.0;
    y1 = height / 2.0;
    z0 = -depth / 2.0;
    z1 = depth / 2.0;

    if (x0 > x1) {
	tmp = x0; x0 = x1; x1 = tmp;
    }
    if (y0 > y1) {
	tmp = y0; y0 = y1; y1 = tmp; 
    }
    if (z0 > z1) {
	tmp = z0; z0 = z1; z1 = tmp; 
    }
    v[0][0] = v[1][0] = v[2][0] = v[3][0] = x0;
    v[4][0] = v[5][0] = v[6][0] = v[7][0] = x1;
    v[0][1] = v[1][1] = v[4][1] = v[5][1] = y0;
    v[2][1] = v[3][1] = v[6][1] = v[7][1] = y1;
    v[0][2] = v[3][2] = v[4][2] = v[7][2] = z0;
    v[1][2] = v[2][2] = v[5][2] = v[6][2] = z1;

    glNewList(base, GL_COMPILE_AND_EXECUTE);
	for (i = 0; i < 6; i++) {
	    glBegin(GL_LINE_LOOP);
		glNormal3fv(&n[i][0]);
		glVertex3fv(&v[faces[i][0]][0]);
		glNormal3fv(&n[i][0]);
		glVertex3fv(&v[faces[i][1]][0]);
		glNormal3fv(&n[i][0]);
		glVertex3fv(&v[faces[i][2]][0]);
		glNormal3fv(&n[i][0]);
		glVertex3fv(&v[faces[i][3]][0]);
	    glEnd();
	}
    glEndList();
}

/******************************************************************************/

void tkSolidBox(GLuint base, float width, float height, float depth)
{
    static float n[6][3] = {
	{-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0},
	{0.0, -1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}
    };
    static GLint faces[6][4] = {
	{0, 1, 2, 3}, {3, 2, 6, 7}, {7, 6, 5, 4},
	{4, 5, 1, 0}, {5, 6, 2, 1}, {7, 4, 0, 3}
    };
    float x0, x1, y0, y1, z0, z1, tmp;
    float v[8][3];
    int i;

    x0 = -width / 2.0;
    x1 = width / 2.0;
    y0 = -height / 2.0;
    y1 = height / 2.0;
    z0 = -depth / 2.0;
    z1 = depth / 2.0;

    if (x0 > x1) {
	tmp = x0; x0 = x1; x1 = tmp;
    }
    if (y0 > y1) {
	tmp = y0; y0 = y1; y1 = tmp; 
    }
    if (z0 > z1) {
	tmp = z0; z0 = z1; z1 = tmp; 
    }
    v[0][0] = v[1][0] = v[2][0] = v[3][0] = x0;
    v[4][0] = v[5][0] = v[6][0] = v[7][0] = x1;
    v[0][1] = v[1][1] = v[4][1] = v[5][1] = y0;
    v[2][1] = v[3][1] = v[6][1] = v[7][1] = y1;
    v[0][2] = v[3][2] = v[4][2] = v[7][2] = z0;
    v[1][2] = v[2][2] = v[5][2] = v[6][2] = z1;

    glNewList(base, GL_COMPILE_AND_EXECUTE);
	for (i = 0; i < 6; i++) {
	    glBegin(GL_POLYGON);
		glNormal3fv(&n[i][0]);
		glVertex3fv(&v[faces[i][0]][0]);
		glNormal3fv(&n[i][0]);
		glVertex3fv(&v[faces[i][1]][0]);
		glNormal3fv(&n[i][0]);
		glVertex3fv(&v[faces[i][2]][0]);
		glNormal3fv(&n[i][0]);
		glVertex3fv(&v[faces[i][3]][0]);
	    glEnd();
	}
    glEndList();
}

/******************************************************************************/

void tkWireTorus(GLuint base, float innerRadius, float outerRadius)
{
    GLint i, j;
    float theta1, phi1, theta2, phi2, rings, sides;
    float v0[03], v1[3], v2[3], v3[3];
    float n0[3], n1[3], n2[3], n3[3];

    rings = 5;
    sides = 10;

    glNewList(base, GL_COMPILE_AND_EXECUTE);
    for (i = 0; i < rings; i++) {
	theta1 = (float)i * 2.0 * PI / rings;
	theta2 = (float)(i + 1) * 2.0 * PI / rings;
	for (j = 0; j < sides; j++) {
	    phi1 = (float)j * 2.0 * PI / sides;
	    phi2 = (float)(j + 1) * 2.0 * PI / sides;

	    v0[0] = cos(theta1) * (outerRadius + innerRadius * cos(phi1));
	    v0[1] = -sin(theta1) * (outerRadius + innerRadius * cos(phi1));
	    v0[2] = innerRadius * sin(phi1);

	    v1[0] = cos(theta2) * (outerRadius + innerRadius * cos(phi1));
	    v1[1] = -sin(theta2) * (outerRadius + innerRadius * cos(phi1));
	    v1[2] = innerRadius * sin(phi1);

	    v2[0] = cos(theta2) * (outerRadius + innerRadius * cos(phi2));
	    v2[1] = -sin(theta2) * (outerRadius + innerRadius * cos(phi2));
	    v2[2] = innerRadius * sin(phi2);

	    v3[0] = cos(theta1) * (outerRadius + innerRadius * cos(phi2));
	    v3[1] = -sin(theta1) * (outerRadius + innerRadius * cos(phi2));
	    v3[2] = innerRadius * sin(phi2);

	    n0[0] = cos(theta1) * (cos(phi1));
	    n0[1] = -sin(theta1) * (cos(phi1));
	    n0[2] = sin(phi1);

	    n1[0] = cos(theta2) * (cos(phi1));
	    n1[1] = -sin(theta2) * (cos(phi1));
	    n1[2] = sin(phi1);

	    n2[0] = cos(theta2) * (cos(phi2));
	    n2[1] = -sin(theta2) * (cos(phi2));
	    n2[2] = sin(phi2);

	    n3[0] = cos(theta1) * (cos(phi2));
	    n3[1] = -sin(theta1) * (cos(phi2));
	    n3[2] = sin(phi2);

	    glBegin(GL_LINE_LOOP);
		glNormal3fv(n3);
		glVertex3fv(v3);
		glNormal3fv(n2);
		glVertex3fv(v2);
		glNormal3fv(n1);
		glVertex3fv(v1);
		glNormal3fv(n0);
		glVertex3fv(v0);
	    glEnd();
	}
    }
    glEndList();
}

/******************************************************************************/

void tkSolidTorus(GLuint base, float innerRadius, float outerRadius)
{
    GLint i, j;
    float theta1, phi1, theta2, phi2, rings, sides;
    float v0[03], v1[3], v2[3], v3[3];
    float n0[3], n1[3], n2[3], n3[3];

    rings = 5;
    sides = 10;

    glNewList(base, GL_COMPILE_AND_EXECUTE);
    for (i = 0; i < rings; i++) {
	theta1 = (float)i * 2.0 * PI / rings;
	theta2 = (float)(i + 1) * 2.0 * PI / rings;
	for (j = 0; j < sides; j++) {
	    phi1 = (float)j * 2.0 * PI / sides;
	    phi2 = (float)(j + 1) * 2.0 * PI / sides;

	    v0[0] = cos(theta1) * (outerRadius + innerRadius * cos(phi1));
	    v0[1] = -sin(theta1) * (outerRadius + innerRadius * cos(phi1));
	    v0[2] = innerRadius * sin(phi1);

	    v1[0] = cos(theta2) * (outerRadius + innerRadius * cos(phi1));
	    v1[1] = -sin(theta2) * (outerRadius + innerRadius * cos(phi1));
	    v1[2] = innerRadius * sin(phi1);

	    v2[0] = cos(theta2) * (outerRadius + innerRadius * cos(phi2));
	    v2[1] = -sin(theta2) * (outerRadius + innerRadius * cos(phi2));
	    v2[2] = innerRadius * sin(phi2);

	    v3[0] = cos(theta1) * (outerRadius + innerRadius * cos(phi2));
	    v3[1] = -sin(theta1) * (outerRadius + innerRadius * cos(phi2));
	    v3[2] = innerRadius * sin(phi2);

	    n0[0] = cos(theta1) * (cos(phi1));
	    n0[1] = -sin(theta1) * (cos(phi1));
	    n0[2] = sin(phi1);

	    n1[0] = cos(theta2) * (cos(phi1));
	    n1[1] = -sin(theta2) * (cos(phi1));
	    n1[2] = sin(phi1);

	    n2[0] = cos(theta2) * (cos(phi2));
	    n2[1] = -sin(theta2) * (cos(phi2));
	    n2[2] = sin(phi2);

	    n3[0] = cos(theta1) * (cos(phi2));
	    n3[1] = -sin(theta1) * (cos(phi2));
	    n3[2] = sin(phi2);

	    glBegin(GL_POLYGON);
		glNormal3fv(n3);
		glVertex3fv(v3);
		glNormal3fv(n2);
		glVertex3fv(v2);
		glNormal3fv(n1);
		glVertex3fv(v1);
		glNormal3fv(n0);
		glVertex3fv(v0);
	    glEnd();
	}
    }
    glEndList();
}

/******************************************************************************/

void tkWireCylinder(GLuint base, float radius, float height)
{
    GLUquadricObj *quadObj;

    glNewList(base, GL_COMPILE_AND_EXECUTE);
	glPushMatrix();
	glRotatef(90.0, 1.0, 0.0, 0.0);
	glTranslatef(0.0, 0.0, -1.0);
	quadObj = gluNewQuadric();
	gluQuadricDrawStyle(quadObj, GLU_LINE);
	gluCylinder(quadObj, radius, radius, height, 12, 2);
	glPopMatrix();
    glEndList();
}

/******************************************************************************/

void tkSolidCylinder(GLuint base, float radius, float height)
{
    GLUquadricObj *quadObj;

    glNewList(base, GL_COMPILE_AND_EXECUTE);
	glPushMatrix();
	glRotatef(90.0, 1.0, 0.0, 0.0);
	glTranslatef(0.0, 0.0, -1.0);
	quadObj = gluNewQuadric();
	gluQuadricDrawStyle(quadObj, GLU_FILL);
	gluQuadricNormals(quadObj, GLU_SMOOTH);
	gluCylinder(quadObj, radius, radius, height, 12, 2);
	glPopMatrix();
    glEndList();
}

/******************************************************************************/

void tkWireCone(GLuint base, float b, float h)
{
    GLUquadricObj *quadObj;

    glNewList(base, GL_COMPILE_AND_EXECUTE);
	quadObj = gluNewQuadric();
	gluQuadricDrawStyle(quadObj, GLU_LINE);
	gluCylinder(quadObj, b, 0.0, h, 15, 10);
    glEndList();
}

/******************************************************************************/

void tkSolidCone(GLuint base, float b, float h)
{
    GLUquadricObj *quadObj;

    glNewList(base, GL_COMPILE_AND_EXECUTE);
	quadObj = gluNewQuadric();
	gluQuadricDrawStyle(quadObj, GLU_FILL);
	gluQuadricNormals(quadObj, GLU_SMOOTH);
	gluCylinder(quadObj, b, 0.0, h, 15, 10);
    glEndList();
}

/******************************************************************************/
