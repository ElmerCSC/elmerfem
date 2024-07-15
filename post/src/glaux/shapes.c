/* shapes.c */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <GL/gl.h>
#include <GL/glu.h>

#include "glaux.h"
#include "3d.h"



#define SPHEREWIRE	0
#define CUBEWIRE	1
#define BOXWIRE		2
#define TORUSWIRE	3
#define CYLINDERWIRE	4
#define ICOSAWIRE	5
#define OCTAWIRE	6
#define TETRAWIRE	7
#define DODECAWIRE	8
#define CONEWIRE	9
#define SPHERESOLID	10
#define CUBESOLID	11
#define BOXSOLID	12
#define TORUSSOLID	13
#define CYLINDERSOLID	14
#define ICOSASOLID	15
#define OCTASOLID	16
#define TETRASOLID	17
#define DODECASOLID	18
#define CONESOLID	19

#define PI 3.1415926535897

/*	structure for each geometric object	*/
typedef struct model {
    GLuint list;	/*  display list to render object   */
    struct model *ptr;	/*  pointer to next object	*/
    int numParam;	/*  # of parameters		*/
    GLdouble *params;	/*  array with parameters	*/
} MODEL, *MODELPTR;

/*	array of linked lists--used to keep track of display lists 
 *	for each different type of geometric object.
 */
static MODELPTR lists[25] = {
    NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL
};

GLuint findList (int lindex, GLdouble *paramArray, int size);
int compareParams (GLdouble *oneArray, GLdouble *twoArray, int size);
GLuint makeModelPtr (int lindex, GLdouble *sizeArray, int count);

static void drawbox(GLdouble, GLdouble, GLdouble, 
	GLdouble, GLdouble, GLdouble, GLenum);
static void doughnut(GLdouble, GLdouble, GLint, GLint, GLenum);
static void icosahedron(GLdouble *, GLdouble, GLenum);
static void octahedron(GLdouble *, GLdouble, GLenum);
static void tetrahedron(GLdouble *, GLdouble, GLenum);
static void subdivide(int, GLdouble *, GLdouble *, GLdouble *,
	GLdouble *, GLdouble, GLenum, int);
static void drawtriangle(int, int, int,
	GLdouble *, GLdouble, GLenum, int);
static void recorditem(GLdouble *, GLdouble *, GLdouble *,
	GLdouble *, GLdouble, GLenum, int);
static void initdodec(void);
static void dodecahedron(GLdouble *, GLdouble, GLenum);
static void pentagon(int, int, int, int, int, GLenum);


/*  Render wire frame or solid sphere.  If no display list with
 *  the current model size exists, create a new display list.
 */
void auxWireSphere (GLdouble radius)
{
    GLUquadricObj *quadObj;
    GLdouble *sizeArray;
    GLuint displayList;

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 1);
    *sizeArray = radius;
    displayList = findList (SPHEREWIRE, sizeArray, 1);

    if (displayList == 0) {
	glNewList(makeModelPtr (SPHEREWIRE, sizeArray, 1),
	    GL_COMPILE_AND_EXECUTE);
	    quadObj = gluNewQuadric ();
	    gluQuadricDrawStyle (quadObj, GLU_LINE);
	    gluSphere (quadObj, radius, 16, 16);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

void auxSolidSphere (GLdouble radius)
{
    GLUquadricObj *quadObj;
    GLdouble *sizeArray;
    GLuint displayList;

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 1);
    *sizeArray = radius;
    displayList = findList (SPHERESOLID, sizeArray, 1);

    if (displayList == 0) {
	glNewList(makeModelPtr (SPHERESOLID, sizeArray, 1),
	    GL_COMPILE_AND_EXECUTE);
	    quadObj = gluNewQuadric ();
	    gluQuadricDrawStyle (quadObj, GLU_FILL);
	    gluQuadricNormals (quadObj, GLU_SMOOTH);
	    gluSphere (quadObj, radius, 16, 16);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

/*  Render wire frame or solid cube.  If no display list with
 *  the current model size exists, create a new display list.
 */
void auxWireCube (GLdouble size)
{
    GLdouble *sizeArray;
    GLuint displayList;

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 1);
    *sizeArray = size;
    displayList = findList (CUBEWIRE, sizeArray, 1);

    if (displayList == 0) {
	glNewList(makeModelPtr (CUBEWIRE, sizeArray, 1),
	    GL_COMPILE_AND_EXECUTE);
	    drawbox(-size/2., size/2., -size/2., size/2., 
		-size/2., size/2., GL_LINE_LOOP);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

void auxSolidCube (GLdouble size)
{
    GLdouble *sizeArray;
    GLuint displayList;

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 1);
    *sizeArray = size;
    displayList = findList (CUBESOLID, sizeArray, 1);

    if (displayList == 0) {
	glNewList(makeModelPtr (CUBESOLID, sizeArray, 1),
	    GL_COMPILE_AND_EXECUTE);
	    drawbox(-size/2., size/2., -size/2., size/2., 
		-size/2., size/2., GL_QUADS);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

/*  Render wire frame or solid cube.  If no display list with
 *  the current model size exists, create a new display list.
 */
void auxWireBox (GLdouble width, GLdouble height, GLdouble depth)
{
    GLdouble *sizeArray, *tmp;
    GLuint displayList;

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 3);
    tmp = sizeArray;
    *tmp++ = width;
    *tmp++ = height;
    *tmp++ = depth;
    displayList = findList (BOXWIRE, sizeArray, 3);

    if (displayList == 0) {
	glNewList(makeModelPtr (BOXWIRE, sizeArray, 3),
	    GL_COMPILE_AND_EXECUTE);
	    drawbox(-width/2., width/2., -height/2., height/2., 
		-depth/2., depth/2., GL_LINE_LOOP);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

void auxSolidBox (GLdouble width, GLdouble height, GLdouble depth)
{
    GLdouble *sizeArray, *tmp;
    GLuint displayList;

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 3);
    tmp = sizeArray;
    *tmp++ = width;
    *tmp++ = height;
    *tmp++ = depth;
    displayList = findList (BOXSOLID, sizeArray, 3);

    if (displayList == 0) {
	glNewList(makeModelPtr (BOXSOLID, sizeArray, 3),
	    GL_COMPILE_AND_EXECUTE);
	    drawbox(-width/2., width/2., -height/2., height/2., 
		-depth/2., depth/2., GL_QUADS);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

/*  Render wire frame or solid tori.  If no display list with
 *  the current model size exists, create a new display list.
 */
void auxWireTorus (GLdouble innerRadius, GLdouble outerRadius)
{
    GLdouble *sizeArray, *tmp;
    GLuint displayList;

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 2);
    tmp = sizeArray;
    *tmp++ = innerRadius;
    *tmp++ = outerRadius;
    displayList = findList (TORUSWIRE, sizeArray, 2);

    if (displayList == 0) {
	glNewList(makeModelPtr (TORUSWIRE, sizeArray, 2),
	    GL_COMPILE_AND_EXECUTE);
	    doughnut(innerRadius, outerRadius, 5, 10, GL_LINE_LOOP);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

void auxSolidTorus (GLdouble innerRadius, GLdouble outerRadius)
{
    GLdouble *sizeArray, *tmp;
    GLuint displayList;

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 2);
    tmp = sizeArray;
    *tmp++ = innerRadius;
    *tmp++ = outerRadius;
    displayList = findList (TORUSSOLID, sizeArray, 2);

    if (displayList == 0) {
	glNewList(makeModelPtr (TORUSSOLID, sizeArray, 2),
	    GL_COMPILE_AND_EXECUTE);
	    doughnut(innerRadius, outerRadius, 8, 15, GL_QUADS);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

/*  Render wire frame or solid cylinders.  If no display list with
 *  the current model size exists, create a new display list.
 */
void auxWireCylinder (GLdouble radius, GLdouble height)
{
    GLUquadricObj *quadObj;
    GLdouble *sizeArray, *tmp;
    GLuint displayList;

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 2);
    tmp = sizeArray;
    *tmp++ = radius;
    *tmp++ = height;
    displayList = findList (CYLINDERWIRE, sizeArray, 2);

    if (displayList == 0) {
	glNewList(makeModelPtr (CYLINDERWIRE, sizeArray, 2),
	    GL_COMPILE_AND_EXECUTE);
	    glPushMatrix ();
	    glRotatef (90.0, 1.0, 0.0, 0.0);
	    glTranslatef (0.0, 0.0, -1.0);
	    quadObj = gluNewQuadric ();
	    gluQuadricDrawStyle (quadObj, GLU_LINE);
	    gluCylinder (quadObj, radius, radius, height, 12, 2);
	    glPopMatrix ();
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

void auxSolidCylinder (GLdouble radius, GLdouble height)
{
    GLUquadricObj *quadObj;
    GLdouble *sizeArray, *tmp;
    GLuint displayList;

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 2);
    tmp = sizeArray;
    *tmp++ = radius;
    *tmp++ = height;
    displayList = findList (CYLINDERWIRE, sizeArray, 2);

    if (displayList == 0) {
	glNewList(makeModelPtr (CYLINDERWIRE, sizeArray, 2),
	    GL_COMPILE_AND_EXECUTE);
	    glPushMatrix ();
	    glRotatef (90.0, 1.0, 0.0, 0.0);
	    glTranslatef (0.0, 0.0, -1.0);
	    quadObj = gluNewQuadric ();
	    gluQuadricDrawStyle (quadObj, GLU_FILL);
	    gluQuadricNormals (quadObj, GLU_SMOOTH);
	    gluCylinder (quadObj, radius, radius, height, 12, 2);
	    glPopMatrix ();
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

/*  Render wire frame or solid icosahedra.  If no display list with
 *  the current model size exists, create a new display list.
 */
void auxWireIcosahedron (GLdouble radius)
{
    GLdouble *sizeArray;
    GLuint displayList;
    GLdouble center[3] = {0.0, 0.0, 0.0};

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 1);
    *sizeArray = radius;
    displayList = findList (ICOSAWIRE, sizeArray, 1);

    if (displayList == 0) {
	glNewList(makeModelPtr (ICOSAWIRE, sizeArray, 1),
	    GL_COMPILE_AND_EXECUTE);
	    icosahedron (center, radius, GL_LINE_LOOP);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

void auxSolidIcosahedron (GLdouble radius)
{
    GLdouble *sizeArray;
    GLuint displayList;
    GLdouble center[3] = {0.0, 0.0, 0.0};

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 1);
    *sizeArray = radius;
    displayList = findList (ICOSASOLID, sizeArray, 1);

    if (displayList == 0) {
	glNewList(makeModelPtr (ICOSASOLID, sizeArray, 1),
	    GL_COMPILE_AND_EXECUTE);
	    icosahedron (center, radius, GL_TRIANGLES);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

/*  Render wire frame or solid octahedra.  If no display list with
 *  the current model size exists, create a new display list.
 */
void auxWireOctahedron (GLdouble radius)
{
    GLdouble *sizeArray;
    GLuint displayList;
    GLdouble center[3] = {0.0, 0.0, 0.0};

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 1);
    *sizeArray = radius;
    displayList = findList (OCTAWIRE, sizeArray, 1);

    if (displayList == 0) {
	glNewList(makeModelPtr (OCTAWIRE, sizeArray, 1),
	    GL_COMPILE_AND_EXECUTE);
	    octahedron (center, radius, GL_LINE_LOOP);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

void auxSolidOctahedron (GLdouble radius)
{
    GLdouble *sizeArray;
    GLuint displayList;
    GLdouble center[3] = {0.0, 0.0, 0.0};

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 1);
    *sizeArray = radius;
    displayList = findList (OCTASOLID, sizeArray, 1);

    if (displayList == 0) {
	glNewList(makeModelPtr (OCTASOLID, sizeArray, 1),
	    GL_COMPILE_AND_EXECUTE);
	    octahedron (center, radius, GL_TRIANGLES);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

/*  Render wire frame or solid tetrahedra.  If no display list with
 *  the current model size exists, create a new display list.
 */
void auxWireTetrahedron (GLdouble radius)
{
    GLdouble *sizeArray;
    GLuint displayList;
    GLdouble center[3] = {0.0, 0.0, 0.0};

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 1);
    *sizeArray = radius;
    displayList = findList (TETRAWIRE, sizeArray, 1);

    if (displayList == 0) {
	glNewList(makeModelPtr (TETRAWIRE, sizeArray, 1),
	    GL_COMPILE_AND_EXECUTE);
	    tetrahedron (center, radius, GL_LINE_LOOP);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

void auxSolidTetrahedron (GLdouble radius)
{
    GLdouble *sizeArray;
    GLuint displayList;
    GLdouble center[3] = {0.0, 0.0, 0.0};

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 1);
    *sizeArray = radius;
    displayList = findList (TETRASOLID, sizeArray, 1);

    if (displayList == 0) {
	glNewList(makeModelPtr (TETRASOLID, sizeArray, 1),
	    GL_COMPILE_AND_EXECUTE);
	    tetrahedron (center, radius, GL_TRIANGLES);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

/*  Render wire frame or solid dodecahedra.  If no display list with
 *  the current model size exists, create a new display list.
 */
void auxWireDodecahedron (GLdouble radius)
{
    GLdouble *sizeArray;
    GLuint displayList;
    GLdouble center[3] = {0.0, 0.0, 0.0};

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 1);
    *sizeArray = radius;
    displayList = findList (DODECAWIRE, sizeArray, 1);

    if (displayList == 0) {
	glNewList(makeModelPtr (DODECAWIRE, sizeArray, 1),
	    GL_COMPILE_AND_EXECUTE);
	    dodecahedron (center, radius/1.73, GL_LINE_LOOP);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

void auxSolidDodecahedron (GLdouble radius)
{
    GLdouble *sizeArray;
    GLuint displayList;
    GLdouble center[3] = {0.0, 0.0, 0.0};

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 1);
    *sizeArray = radius;
    displayList = findList (DODECASOLID, sizeArray, 1);

    if (displayList == 0) {
	glNewList(makeModelPtr (DODECASOLID, sizeArray, 1),
	    GL_COMPILE_AND_EXECUTE);
	    dodecahedron (center, radius/1.73, GL_TRIANGLE_FAN);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

/*  Render wire frame or solid cones.  If no display list with
 *  the current model size exists, create a new display list.
 */
void auxWireCone (GLdouble base, GLdouble height)
{
    GLUquadricObj *quadObj;
    GLdouble *sizeArray, *tmp;
    GLuint displayList;

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 2);
    tmp = sizeArray;
    *tmp++ = base;
    *tmp++ = height;
    displayList = findList (CONEWIRE, sizeArray, 2);

    if (displayList == 0) {
	glNewList(makeModelPtr (CONEWIRE, sizeArray, 2),
	    GL_COMPILE_AND_EXECUTE);
	    quadObj = gluNewQuadric ();
	    gluQuadricDrawStyle (quadObj, GLU_LINE);
	    gluCylinder (quadObj, base, 0.0, height, 15, 10);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

void auxSolidCone (GLdouble base, GLdouble height)
{
    GLUquadricObj *quadObj;
    GLdouble *sizeArray, *tmp;
    GLuint displayList;

    sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 2);
    tmp = sizeArray;
    *tmp++ = base;
    *tmp++ = height;
    displayList = findList (CONEWIRE, sizeArray, 2);

    if (displayList == 0) {
	glNewList(makeModelPtr (CONEWIRE, sizeArray, 2),
	    GL_COMPILE_AND_EXECUTE);
	    quadObj = gluNewQuadric ();
	    gluQuadricDrawStyle (quadObj, GLU_FILL);
	    gluQuadricNormals (quadObj, GLU_SMOOTH);
	    gluCylinder (quadObj, base, 0.0, height, 15, 10);
	glEndList();
    }
    else {
	glCallList(displayList);
	free (sizeArray);
    }
}

/* Routines to build 3 dimensional solids, including:
 *
 * drawbox, doughnut, icosahedron, 
 * octahedron, tetrahedron, dodecahedron.
 */

/* drawbox:
 *
 * draws a rectangular box with the given x, y, and z ranges.  
 * The box is axis-aligned.
 */
static void drawbox(GLdouble x0, GLdouble x1, GLdouble y0, GLdouble y1,
	GLdouble z0, GLdouble z1, GLenum type)
{
    static GLdouble n[6][3] = {
	{-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0},
	{0.0, -1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}
    };
    static GLint faces[6][4] = {
	{ 0, 1, 2, 3 }, { 3, 2, 6, 7 }, { 7, 6, 5, 4 },
	{ 4, 5, 1, 0 }, { 5, 6, 2, 1 }, { 7, 4, 0, 3 }
    };
    GLdouble v[8][3], tmp;
    GLint i;

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

    for (i = 0; i < 6; i++) {
	glBegin(type);
	glNormal3dv(&n[i][0]);
	glVertex3dv(&v[faces[i][0]][0]);
	glNormal3dv(&n[i][0]);
	glVertex3dv(&v[faces[i][1]][0]);
	glNormal3dv(&n[i][0]);
	glVertex3dv(&v[faces[i][2]][0]);
	glNormal3dv(&n[i][0]);
	glVertex3dv(&v[faces[i][3]][0]);
	glEnd();
    }
}

/* doughnut:
 *
 * draws a doughnut, centered at (0, 0, 0) whose axis is aligned with
 * the z-axis.  The doughnut's major radius is R, and minor radius is r.
 */

static void doughnut(GLdouble r, GLdouble R, GLint nsides, GLint rings, GLenum type)
{
    int	i, j;
    GLdouble	theta, phi, theta1, phi1;
    GLdouble	p0[03], p1[3], p2[3], p3[3];
    GLdouble	n0[3], n1[3], n2[3], n3[3];

    for (i = 0; i < rings; i++) {
	theta = (GLdouble)i*2.0*PI/rings;
	theta1 = (GLdouble)(i+1)*2.0*PI/rings;
	for (j = 0; j < nsides; j++) {
	    phi = (GLdouble)j*2.0*PI/nsides;
	    phi1 = (GLdouble)(j+1)*2.0*PI/nsides;

	    p0[0] = cos(theta)*(R + r*cos(phi));
	    p0[1] = -sin(theta)*(R + r*cos(phi));
	    p0[2] = r*sin(phi);

	    p1[0] = cos(theta1)*(R + r*cos(phi));
	    p1[1] = -sin(theta1)*(R + r*cos(phi));
	    p1[2] = r*sin(phi);

	    p2[0] = cos(theta1)*(R + r*cos(phi1));
	    p2[1] = -sin(theta1)*(R + r*cos(phi1));
	    p2[2] = r*sin(phi1);

	    p3[0] = cos(theta)*(R + r*cos(phi1));
	    p3[1] = -sin(theta)*(R + r*cos(phi1));
	    p3[2] = r*sin(phi1);

	    n0[0] = cos(theta)*(cos(phi));
	    n0[1] = -sin(theta)*(cos(phi));
	    n0[2] = sin(phi);

	    n1[0] = cos(theta1)*(cos(phi));
	    n1[1] = -sin(theta1)*(cos(phi));
	    n1[2] = sin(phi);

	    n2[0] = cos(theta1)*(cos(phi1));
	    n2[1] = -sin(theta1)*(cos(phi1));
	    n2[2] = sin(phi1);

	    n3[0] = cos(theta)*(cos(phi1));
	    n3[1] = -sin(theta)*(cos(phi1));
	    n3[2] = sin(phi1);

	    m_xformpt(p0, p0, n0, n0);
	    m_xformpt(p1, p1, n1, n1);
	    m_xformpt(p2, p2, n2, n2);
	    m_xformpt(p3, p3, n3, n3);

	    glBegin(type);
		glNormal3dv(n3);
		glVertex3dv(p3);
		glNormal3dv(n2);
		glVertex3dv(p2);
		glNormal3dv(n1);
		glVertex3dv(p1);
		glNormal3dv(n0);
		glVertex3dv(p0);
	    glEnd();
	}
    }
}

/* octahedron data: The octahedron produced is centered 
 * at the origin and has radius 1.0 
 */
static GLdouble odata[6][3] = {
  {1.0, 0.0, 0.0},
  {-1.0, 0.0, 0.0},
  {0.0, 1.0, 0.0},
  {0.0, -1.0, 0.0},
  {0.0, 0.0, 1.0},
  {0.0, 0.0, -1.0}
};

static int ondex[8][3] = {
    {0, 4, 2}, {1, 2, 4}, {0, 3, 4}, {1, 4, 3},
    {0, 2, 5}, {1, 5, 2}, {0, 5, 3}, {1, 3, 5}
};

/* tetrahedron data: */

#define T	1.73205080756887729

static GLdouble tdata[4][3] = {
    {T, T, T}, {T, -T, -T}, {-T, T, -T}, {-T, -T, T}
};

static int tndex[4][3] = {
    {0, 1, 3}, {2, 1, 0}, {3, 2, 0}, {1, 2, 3}
};

/* icosahedron data: These numbers are rigged to 
 * make an icosahedron of radius 1.0 
 */

#define X .525731112119133606
#define Z .850650808352039932

static GLdouble idata[12][3] = {
  {-X, 0.0, Z},
  {X, 0.0, Z},
  {-X, 0.0, -Z},
  {X, 0.0, -Z},
  {0.0, Z, X},
  {0.0, Z, -X},
  {0.0, -Z, X},
  {0.0, -Z, -X},
  {Z, X, 0.0},
  {-Z, X, 0.0},
  {Z, -X, 0.0},
  {-Z, -X, 0.0},
};

static int iindex[20][3] = {
    {0, 4, 1},    {0, 9, 4},
    {9, 5, 4},    {4, 5, 8},
    {4, 8, 1},    {8, 10, 1},
    {8, 3, 10},    {5, 3, 8},
    {5, 2, 3},    {2, 7, 3},
    {7, 10, 3},    {7, 6, 10},
    {7, 11, 6},    {11, 0, 6},
    {0, 1, 6},    {6, 1, 10},
    {9, 0, 11},    {9, 11, 2},
    {9, 2, 5},    {7, 2, 11},
};

/* icosahedron:
 *
 * Draws an icosahedron with center at p0 having the
 * given radius.
 */

static void icosahedron(GLdouble p0[3], GLdouble radius, GLenum shadeType)
{
    int i;

    for (i = 0; i < 20; i++)
	drawtriangle(i, 0, 1, p0, radius, shadeType, 0);
}

/* octahedron:
 *
 * Draws an octahedron with center at p0 having the
 * given radius.
 */
static void octahedron(GLdouble p0[3], GLdouble radius, GLenum shadeType)
{
    int i;

    for (i = 0; i < 8; i++)
	drawtriangle(i, 1, 1, p0, radius, shadeType, 0);
}

/* tetrahedron:
 *
 * Draws an tetrahedron with center at p0 having the
 * given radius.
 */

static void tetrahedron(GLdouble p0[3], GLdouble radius, GLenum shadeType)
{
    int i;

    for (i = 0; i < 4; i++)
	drawtriangle(i, 2, 1, p0, radius, shadeType, 0);
}

static void subdivide(int depth, GLdouble *v0, GLdouble *v1, GLdouble *v2,
	GLdouble p0[3], GLdouble radius, GLenum shadeType, int avnormal)
{
    GLdouble w0[3], w1[3], w2[3];
    GLdouble l;
    int i, j, k, n;

    for (i = 0; i < depth; i++)
	for (j = 0; i + j < depth; j++) {
	    k = depth - i - j;
	    for (n = 0; n < 3; n++) {
		w0[n] = (i*v0[n] + j*v1[n] + k*v2[n])/depth;
		w1[n] = ((i+1)*v0[n] + j*v1[n] + (k-1)*v2[n])/depth;
		w2[n] = (i*v0[n] + (j+1)*v1[n] + (k-1)*v2[n])/depth;
	    }
	    l = sqrt(w0[0]*w0[0] + w0[1]*w0[1] + w0[2]*w0[2]);
	    w0[0] /= l; w0[1] /= l; w0[2] /= l;
	    l = sqrt(w1[0]*w1[0] + w1[1]*w1[1] + w1[2]*w1[2]);
	    w1[0] /= l; w1[1] /= l; w1[2] /= l;
	    l = sqrt(w2[0]*w2[0] + w2[1]*w2[1] + w2[2]*w2[2]);
	    w2[0] /= l; w2[1] /= l; w2[2] /= l;
	    recorditem(w1, w0, w2, p0, radius, shadeType, avnormal);
	}
    for (i = 0; i < depth-1; i++)
	for (j = 0; i + j < depth-1; j++) {
	    k = depth - i - j;
	    for (n = 0; n < 3; n++) {
		w0[n] = ((i+1)*v0[n] + (j+1)*v1[n] + (k-2)*v2[n])/depth;
		w1[n] = ((i+1)*v0[n] + j*v1[n] + (k-1)*v2[n])/depth;
		w2[n] = (i*v0[n] + (j+1)*v1[n] + (k-1)*v2[n])/depth;
	    }
	    l = sqrt(w0[0]*w0[0] + w0[1]*w0[1] + w0[2]*w0[2]);
	    w0[0] /= l; w0[1] /= l; w0[2] /= l;
	    l = sqrt(w1[0]*w1[0] + w1[1]*w1[1] + w1[2]*w1[2]);
	    w1[0] /= l; w1[1] /= l; w1[2] /= l;
	    l = sqrt(w2[0]*w2[0] + w2[1]*w2[1] + w2[2]*w2[2]);
	    w2[0] /= l; w2[1] /= l; w2[2] /= l;
	    recorditem(w0, w1, w2, p0, radius, shadeType, avnormal);
	}
}

static void drawtriangle(int i, int geomType, int depth,
	GLdouble p0[3], GLdouble radius, GLenum shadeType, int avnormal)
{
    GLdouble *x0, *x1, *x2;

    switch (geomType) {
	case 0:	/* icosahedron */
	    x0 = &idata[iindex[i][0]][0];
	    x1 = &idata[iindex[i][1]][0];
	    x2 = &idata[iindex[i][2]][0];
	    break;
	case 1: /* octahedron */
	    x0 = &odata[ondex[i][0]][0];
	    x1 = &odata[ondex[i][1]][0];
	    x2 = &odata[ondex[i][2]][0];
	    break;
	case 2: /* tetrahedron */
	    x0 = &tdata[tndex[i][0]][0];
	    x1 = &tdata[tndex[i][1]][0];
	    x2 = &tdata[tndex[i][2]][0];
	    break;
    }
    subdivide(depth, x0, x1, x2, p0, radius, shadeType, avnormal);
}

static void recorditem(GLdouble *n1, GLdouble *n2, GLdouble *n3,
	GLdouble center[3], GLdouble radius, GLenum shadeType, int avnormal)
{
    GLdouble p1[3], p2[3], p3[3], q0[3], q1[3], n11[3], n22[3], n33[3];
    int	i;

    for (i = 0; i < 3; i++) {
	p1[i] = n1[i]*radius + center[i];
	p2[i] = n2[i]*radius + center[i];
	p3[i] = n3[i]*radius + center[i];
    }
    if (avnormal == 0) {
	diff3(p1, p2, q0);
	diff3(p2, p3, q1);
	crossprod(q0, q1, q1);
	normalize(q1);
	m_xformpt(p1, p1, q1, n11);
	m_xformptonly(p2, p2);
	m_xformptonly(p3, p3);

	glBegin (shadeType);
	glNormal3dv(n11);
	glVertex3dv(p1);
	glVertex3dv(p2);
	glVertex3dv(p3);
	glEnd();
	return;
    }
    m_xformpt(p1, p1, n1, n11);
    m_xformpt(p2, p2, n2, n22);
    m_xformpt(p3, p3, n3, n33);

    glBegin (shadeType);
    glNormal3dv(n11);
    glVertex3dv(p1);
    glNormal3dv(n22);
    glVertex3dv(p2);
    glNormal3dv(n33);
    glVertex3dv(p3);
    glEnd();
}

static GLdouble dodec[20][3];

static void initdodec()
{
    GLdouble alpha, beta;

    alpha = sqrt(2.0/(3.0 + sqrt(5.0)));
    beta = 1.0 + sqrt(6.0/(3.0 + sqrt(5.0)) - 2.0 + 2.0*sqrt(2.0/(3.0 +
							    sqrt(5.0))));
    dodec[0][0] = -alpha; dodec[0][1] = 0; dodec[0][2] = beta;
    dodec[1][0] = alpha; dodec[1][1] = 0; dodec[1][2] = beta;
    dodec[2][0] = -1; dodec[2][1] = -1; dodec[2][2] = -1;
    dodec[3][0] = -1; dodec[3][1] = -1; dodec[3][2] = 1;
    dodec[4][0] = -1; dodec[4][1] = 1; dodec[4][2] = -1;
    dodec[5][0] = -1; dodec[5][1] = 1; dodec[5][2] = 1;
    dodec[6][0] = 1; dodec[6][1] = -1; dodec[6][2] = -1;
    dodec[7][0] = 1; dodec[7][1] = -1; dodec[7][2] = 1;
    dodec[8][0] = 1; dodec[8][1] = 1; dodec[8][2] = -1;
    dodec[9][0] = 1; dodec[9][1] = 1; dodec[9][2] = 1;
    dodec[10][0] = beta; dodec[10][1] = alpha; dodec[10][2] = 0;
    dodec[11][0] = beta; dodec[11][1] = -alpha; dodec[11][2] = 0;
    dodec[12][0] = -beta; dodec[12][1] = alpha; dodec[12][2] = 0;
    dodec[13][0] = -beta; dodec[13][1] = -alpha; dodec[13][2] = 0;
    dodec[14][0] = -alpha; dodec[14][1] = 0; dodec[14][2] = -beta;
    dodec[15][0] = alpha; dodec[15][1] = 0; dodec[15][2] = -beta;
    dodec[16][0] = 0; dodec[16][1] = beta; dodec[16][2] = alpha;
    dodec[17][0] = 0; dodec[17][1] = beta; dodec[17][2] = -alpha;
    dodec[18][0] = 0; dodec[18][1] = -beta; dodec[18][2] = alpha;
    dodec[19][0] = 0; dodec[19][1] = -beta; dodec[19][2] = -alpha;
}

/* dodecahedron:
 *
 * Draws an dodecahedron with center at 0.0. The radius
 * is sqrt(3).
 */
static void dodecahedron(GLdouble center[3], GLdouble sc, GLenum type)
{
    static int inited = 0;

    if ( inited == 0) {
	inited = 1;
	initdodec();
    }
    m_pushmatrix();
    m_translate(center[0], center[1], center[2]);
    m_scale(sc, sc, sc);
    pentagon(0, 1, 9, 16, 5, type);
    pentagon(1, 0, 3, 18, 7, type);
    pentagon(1, 7, 11, 10, 9, type);
    pentagon(11, 7, 18, 19, 6, type);
    pentagon(8, 17, 16, 9, 10, type);
    pentagon(2, 14, 15, 6, 19, type);
    pentagon(2, 13, 12, 4, 14, type);
    pentagon(2, 19, 18, 3, 13, type);
    pentagon(3, 0, 5, 12, 13, type);
    pentagon(6, 15, 8, 10, 11, type);
    pentagon(4, 17, 8, 15, 14, type);
    pentagon(4, 12, 5, 16, 17, type);
    m_popmatrix();
}

static void pentagon(int a, int b, int c, int d, int e, GLenum shadeType)
{
    GLdouble n0[3], d1[3], d2[3], d3[3], d4[3], d5[3], nout[3];

    diff3(&dodec[a][0], &dodec[b][0], d1);
    diff3(&dodec[b][0], &dodec[c][0], d2);
    crossprod(d1, d2, n0);
    normalize(n0);
    m_xformpt(&dodec[a][0], d1, n0, nout);
    m_xformptonly(&dodec[b][0], d2);
    m_xformptonly(&dodec[c][0], d3);
    m_xformptonly(&dodec[d][0], d4);
    m_xformptonly(&dodec[e][0], d5);

    glBegin (shadeType);
    glNormal3dv(nout);
    glVertex3dv(d1);
    glVertex3dv(d2);
    glVertex3dv(d3);
    glVertex3dv(d4);
    glVertex3dv(d5);
    glEnd();
}

/*	linked lists--display lists for each different 
 *	type of geometric objects.  The linked list is 
 *	searched, until an object of the requested
 *	size is found.  If no geometric object of that size
 *	has been previously made, a new one is created.
 */
GLuint findList (int lindex, GLdouble *paramArray, int size) 
{
    MODELPTR endList;
    int found = 0;
    
    endList = lists[lindex];
    while (endList != NULL) {
	if (compareParams (endList->params, paramArray, size))
	    return (endList->list);
	endList = endList->ptr;
    }
/*  if not found, return 0 and calling routine should
 *  make a new list	
 */
    return (0);
}

int compareParams (GLdouble *oneArray, GLdouble *twoArray, int size) 
{
    int i;
    int matches = 1;

    for (i = 0; (i < size) && matches; i++) {
	if (*oneArray++ != *twoArray++)
	    matches = 0;
    }
    return (matches);
}

GLuint makeModelPtr (int lindex, GLdouble *sizeArray, int count)
{
    MODELPTR newModel;

    newModel = (MODELPTR) malloc (sizeof (MODEL));
    newModel->list = glGenLists (1);
    newModel->numParam = count;
    newModel->params = sizeArray;
    newModel->ptr = lists[lindex];
    lists[lindex] = newModel;

    return (newModel->list);
}
