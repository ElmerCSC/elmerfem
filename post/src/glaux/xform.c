/* xform.c */

#include <math.h>
#include <stdio.h>

#include <GL/gl.h>
#include <GL/glu.h>

#include "3d.h"



#define STACKDEPTH 10

typedef struct {
    GLdouble	mat[4][4];
    GLdouble	norm[3][3];
} mat_t;

static mat_t matstack[STACKDEPTH] = {
    {{{1.0, 0.0, 0.0, 0.0},
    {0.0, 1.0, 0.0, 0.0},
    {0.0, 0.0, 1.0, 0.0},
    {0.0, 0.0, 0.0, 1.0}},
    {{1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0}}}
};
static int identitymat = 1;

static int mattop = 0;

void m_xformpt(GLdouble pin[3], GLdouble pout[3], 
    GLdouble nin[3], GLdouble nout[3])
{
    int	i;
    GLdouble	ptemp[3], ntemp[3];
    mat_t	*m = &matstack[mattop];

    if (identitymat) {
	for (i = 0; i < 3; i++) {
	    pout[i] = pin[i];
	    nout[i] = nin[i];
	}
	return;
    }
    for (i = 0; i < 3; i++) {
	ptemp[i] = pin[0]*m->mat[0][i] +
		   pin[1]*m->mat[1][i] +
		   pin[2]*m->mat[2][i] +
		   m->mat[3][i];
	ntemp[i] = nin[0]*m->norm[0][i] +
		   nin[1]*m->norm[1][i] +
		   nin[2]*m->norm[2][i];
    }
    for (i = 0; i < 3; i++) {
	pout[i] = ptemp[i];
	nout[i] = ntemp[i];
    }
    normalize(nout);
}

void m_xformptonly(GLdouble pin[3], GLdouble pout[3])
{
    int	i;
    GLdouble	ptemp[3];
    mat_t	*m = &matstack[mattop];

    if (identitymat) {
	for (i = 0; i < 3; i++) {
	    pout[i] = pin[i];
	}
	return;
    }
     for (i = 0; i < 3; i++) {
	ptemp[i] = pin[0]*m->mat[0][i] +
		   pin[1]*m->mat[1][i] +
		   pin[2]*m->mat[2][i] +
		   m->mat[3][i];
    }
    for (i = 0; i < 3; i++) {
	pout[i] = ptemp[i];
    }
}

void m_pushmatrix(void)
{
  if (mattop < STACKDEPTH-1) {
    matstack[mattop+1] = matstack[mattop];
    mattop++;
  } else
    glaux_error("m_pushmatrix: stack overflow\n");
}

void m_popmatrix(void)
{
  if (mattop > 0)
    mattop--;
  else
    glaux_error("m_popmatrix: stack underflow\n");
}

void m_translate(GLdouble x, GLdouble y, GLdouble z)
{
    int	i;
    mat_t	*m = &matstack[mattop];

    identitymat = 0;
    for (i = 0; i < 4; i++)
	m->mat[3][i] = x*m->mat[0][i] +
				 y*m->mat[1][i] +
				 z*m->mat[2][i] +
				 m->mat[3][i];
}

void m_scale(GLdouble x, GLdouble y, GLdouble z)
{
    int	i;
    mat_t	*m = &matstack[mattop];

    identitymat = 0;
    for (i = 0; i < 3; i++) {
	m->mat[0][i] *= x;
	m->mat[1][i] *= y;
	m->mat[2][i] *= z;
    }
    for (i = 0; i < 3; i++) {
	m->norm[0][i] /= x;
	m->norm[1][i] /= y;
	m->norm[2][i] /= z;
    }
}
