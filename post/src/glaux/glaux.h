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
 * KIND, OR ANY DEMAGES WHATSOEVER, INCLUDING WITHOUT LIMITATION,
 * LOSS OF PROFIT, LOSS OF USE, SAVINGS OR REVENUE, OR THE CLAIMS OF
 * THIRD PARTIES, WHETHER OR NOT SILICON GRAPHICS, INC.  HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH LOSS, HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THEj
 * POSSESSION, USE OR PREFORMANCE OF THIS SOFTWARE.
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
#ifndef __aux_h__
#define __aux_h__

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <GL/gl.h>
#include <GL/glu.h>

#ifdef __cplusplus
extern "C" {
#endif


/*
** ToolKit Window Types
** In the future, AUX_RGBA may be a combination of both RGB and ALPHA
*/

#define AUX_RGB		0
#define AUX_RGBA	AUX_RGB
#define AUX_INDEX	1
#define AUX_SINGLE	0
#define AUX_DOUBLE	2
#define AUX_DIRECT	0
#define AUX_INDIRECT	4

#define AUX_ACCUM	8
#define AUX_ALPHA	16
#define AUX_DEPTH	32
#define AUX_STENCIL	64
#define AUX_AUX		128

/* 
** Window Masks
*/

#define AUX_WIND_IS_RGB(x)	(((x) & AUX_INDEX) == 0)
#define AUX_WIND_IS_INDEX(x)	(((x) & AUX_INDEX) != 0)
#define AUX_WIND_IS_SINGLE(x)	(((x) & AUX_DOUBLE) == 0)
#define AUX_WIND_IS_DOUBLE(x)	(((x) & AUX_DOUBLE) != 0)
#define AUX_WIND_IS_INDIRECT(x)	(((x) & AUX_INDIRECT) != 0)
#define AUX_WIND_IS_DIRECT(x)	(((x) & AUX_INDIRECT) == 0)
#define AUX_WIND_HAS_ACCUM(x)	(((x) & AUX_ACCUM) != 0)
#define AUX_WIND_HAS_ALPHA(x)	(((x) & AUX_ALPHA) != 0)
#define AUX_WIND_HAS_DEPTH(x)	(((x) & AUX_DEPTH) != 0)
#define AUX_WIND_HAS_STENCIL(x)	(((x) & AUX_STENCIL) != 0)

/*
** ToolKit Event Structure
*/

typedef struct _AUX_EVENTREC {
    GLint event;
    GLint data[4];
} AUX_EVENTREC;

/* 
** ToolKit Event Types
*/
#define AUX_EXPOSE	1
#define AUX_CONFIG	2
#define AUX_DRAW	4
#define AUX_KEYEVENT	8
#define AUX_MOUSEDOWN	16
#define AUX_MOUSEUP	32	
#define AUX_MOUSELOC	64

/*
** Toolkit Event Data Indices
*/
#define AUX_WINDOWX		0
#define AUX_WINDOWY		1
#define AUX_MOUSEX		0
#define AUX_MOUSEY		1
#define AUX_MOUSESTATUS		3
#define AUX_KEY			0
#define AUX_KEYSTATUS		1

/*
** ToolKit Event Status Messages
*/
#define	AUX_LEFTBUTTON		1
#define	AUX_RIGHTBUTTON		2
#define	AUX_MIDDLEBUTTON	4
#define	AUX_SHIFT		1
#define	AUX_CONTROL		2

/* 
** ToolKit Key Codes
*/
#define AUX_RETURN		0x0D
#define AUX_ESCAPE		0x1B
#define AUX_SPACE		0x20
#define AUX_LEFT		0x25
#define AUX_UP			0x26
#define AUX_RIGHT		0x27
#define AUX_DOWN		0x28
#define AUX_A			'A'
#define AUX_B			'B'
#define AUX_C			'C'
#define AUX_D			'D'
#define AUX_E			'E'
#define AUX_F			'F'
#define AUX_G			'G'
#define AUX_H			'H'
#define AUX_I			'I'
#define AUX_J			'J'
#define AUX_K			'K'
#define AUX_L			'L'
#define AUX_M			'M'
#define AUX_N			'N'
#define AUX_O			'O'
#define AUX_P			'P'
#define AUX_Q			'Q'
#define AUX_R			'R'
#define AUX_S			'S'
#define AUX_T			'T'
#define AUX_U			'U'
#define AUX_V			'V'
#define AUX_W			'W'
#define AUX_X			'X'
#define AUX_Y			'Y'
#define AUX_Z			'Z'
#define AUX_a			'a'
#define AUX_b			'b'
#define AUX_c			'c'
#define AUX_d			'd'
#define AUX_e			'e'
#define AUX_f			'f'
#define AUX_g			'g'
#define AUX_h			'h'
#define AUX_i			'i'
#define AUX_j			'j'
#define AUX_k			'k'
#define AUX_l			'l'
#define AUX_m			'm'
#define AUX_n			'n'
#define AUX_o			'o'
#define AUX_p			'p'
#define AUX_q			'q'
#define AUX_r			'r'
#define AUX_s			's'
#define AUX_t			't'
#define AUX_u			'u'
#define AUX_v			'v'
#define AUX_w			'w'
#define AUX_x			'x'
#define AUX_y			'y'
#define AUX_z			'z'
#define AUX_0			'0'
#define AUX_1			'1'
#define AUX_2			'2'
#define AUX_3			'3'
#define AUX_4			'4'
#define AUX_5			'5'
#define AUX_6			'6'
#define AUX_7			'7'
#define AUX_8			'8'
#define AUX_9			'9'

/*
** ToolKit Gets and Sets
*/
#define AUX_FD			1  /* return fd (long) */
#define AUX_COLORMAP		3  /* pass buf of r, g and b (unsigned char) */
#define AUX_GREYSCALEMAP	4
#define AUX_FOGMAP		5  /* pass fog and color bits (long) */
#define AUX_ONECOLOR		6  /* pass index, r, g, and b (long) */

/*
** Color Macros
*/

enum {
    AUX_BLACK = 0,
    AUX_RED,
    AUX_GREEN,
    AUX_YELLOW,
    AUX_BLUE,
    AUX_MAGENTA,
    AUX_CYAN,
    AUX_WHITE
};

extern float auxRGBMap[8][3];

#define AUX_SETCOLOR(x, y) (AUX_WIND_IS_RGB((x)) ? \
		           glColor3fv(auxRGBMap[(y)]) : glIndexf((y)))

/*
** RGB Image Structure
*/

typedef struct _AUX_RGBImageRec {
    GLint sizeX, sizeY;
    unsigned char *data;
} AUX_RGBImageRec;

/*
** Prototypes
*/

extern void auxInitDisplayMode(GLenum);
extern void auxInitPosition(int, int, int, int);
extern GLenum auxInitWindow(char *);
extern void auxCloseWindow(void);
extern void auxQuit(void);
extern void auxSwapBuffers(void);

extern Display *auxXDisplay(GLvoid);
extern Window auxXWindow(GLvoid);

extern void auxMainLoop(void (*)());
extern void auxExposeFunc(void (*)(int, int));
extern void auxReshapeFunc(void (*)(GLsizei, GLsizei));
extern void auxIdleFunc(void (*)());
extern void auxKeyFunc(int, void (*)());
extern void auxMouseFunc(int, int, void (*)(AUX_EVENTREC *));

extern int auxGetColorMapSize(void);
extern void auxGetMouseLoc(int *, int *);
extern void auxSetOneColor(int, float, float, float);
extern void auxSetFogRamp(int, int);
extern void auxSetGreyRamp(void);
extern void auxSetRGBMap(int, float *);

extern AUX_RGBImageRec *auxRGBImageLoad(char *);

extern void auxCreateFont(void);
extern void auxDrawStr(char *);

extern void auxWireSphere(GLdouble);
extern void auxSolidSphere(GLdouble);
extern void auxWireCube(GLdouble);
extern void auxSolidCube(GLdouble);
extern void auxWireBox(GLdouble, GLdouble, GLdouble);
extern void auxSolidBox(GLdouble, GLdouble, GLdouble);
extern void auxWireTorus(GLdouble, GLdouble);
extern void auxSolidTorus(GLdouble, GLdouble);
extern void auxWireCylinder(GLdouble, GLdouble);
extern void auxSolidCylinder(GLdouble, GLdouble);
extern void auxWireIcosahedron(GLdouble);
extern void auxSolidIcosahedron(GLdouble);
extern void auxWireOctahedron(GLdouble);
extern void auxSolidOctahedron(GLdouble);
extern void auxWireTetrahedron(GLdouble);
extern void auxSolidTetrahedron(GLdouble);
extern void auxWireDodecahedron(GLdouble);
extern void auxSolidDodecahedron(GLdouble);
extern void auxWireCone(GLdouble, GLdouble);
extern void auxSolidCone(GLdouble, GLdouble);
extern void auxWireTeapot(GLdouble);
extern void auxSolidTeapot(GLdouble);

#ifdef __cplusplus
}
#endif

#endif /* __aux_h__ */
