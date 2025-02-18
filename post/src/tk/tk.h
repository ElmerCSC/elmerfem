/* tk.h */

/*
 * Mesa - a 3-D graphics library
 * Version:  1.0 beta
 * Copyright (C) 1995  Brian Paul  (brianp@ssec.wisc.edu)
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


/*
 * This code originated from SGI and was modified for use with Mesa.
 */


#ifndef TK_INCLUDED
#define TK_INCLUDED


#ifdef __cplusplus
extern "C" {
#endif

  
    
#include <GL/gl.h>
#include <GL/glu.h>


/*
** Nano Window Toolkit.
** Version 1.1
*/


#ifdef __WIN32__
/*
** Display Mode Selection Criteria
*/
enum {
    TK_USE_ID = 1,
    TK_EXACT_MATCH,
    TK_MINIMUM_CRITERIA
};
#endif


/*
** Window Types
*/

#define TK_RGB		0
#define TK_INDEX	1
#define TK_SINGLE	0
#define TK_DOUBLE	2
#define TK_DIRECT	0
#define TK_INDIRECT	4
#define TK_ACCUM	8
#define TK_ALPHA	16
#define TK_DEPTH	32
#define TK_STENCIL	64
#define TK_OVERLAY	128
#define TK_UNDERLAY	256

/* 
** Window Masks
*/

#define TK_IS_RGB(x)		(((x) & TK_INDEX) == 0)
#define TK_IS_INDEX(x)		(((x) & TK_INDEX) != 0)
#define TK_IS_SINGLE(x)		(((x) & TK_DOUBLE) == 0)
#define TK_IS_DOUBLE(x)		(((x) & TK_DOUBLE) != 0)
#define TK_IS_DIRECT(x)		(((x) & TK_INDIRECT) == 0)
#define TK_IS_INDIRECT(x)	(((x) & TK_INDIRECT) != 0)
#define TK_HAS_ACCUM(x)		(((x) & TK_ACCUM) != 0)
#define TK_HAS_ALPHA(x)		(((x) & TK_ALPHA) != 0)
#define TK_HAS_DEPTH(x)		(((x) & TK_DEPTH) != 0)
#define TK_HAS_OVERLAY(x)	(((x) & TK_OVERLAY) != 0)
#define TK_HAS_UNDERLAY(x)	(((x) & TK_UNDERLAY) != 0)
#define TK_HAS_STENCIL(x)	(((x) & TK_STENCIL) != 0)

/*
** Windowing System Specific Gets
*/

typedef enum {
    TK_X_DISPLAY = 1,
    TK_X_WINDOW
} TKenum;

/*
** Event Status
*/

#define	TK_LEFTBUTTON		1
#define	TK_RIGHTBUTTON		2
#define	TK_MIDDLEBUTTON		4
#define	TK_SHIFT		1
#define	TK_CONTROL		2

/* 
** Key Codes
*/

#define TK_RETURN		0x0D
#define TK_ESCAPE		0x1B
#define TK_SPACE		0x20
#define TK_LEFT			0x25
#define TK_UP			0x26
#define TK_RIGHT		0x27
#define TK_DOWN			0x28
#define TK_A			'A'
#define TK_B			'B'
#define TK_C			'C'
#define TK_D			'D'
#define TK_E			'E'
#define TK_F			'F'
#define TK_G			'G'
#define TK_H			'H'
#define TK_I			'I'
#define TK_J			'J'
#define TK_K			'K'
#define TK_L			'L'
#define TK_M			'M'
#define TK_N			'N'
#define TK_O			'O'
#define TK_P			'P'
#define TK_Q			'Q'
#define TK_R			'R'
#define TK_S			'S'
#define TK_T			'T'
#define TK_U			'U'
#define TK_V			'V'
#define TK_W			'W'
#define TK_X			'X'
#define TK_Y			'Y'
#define TK_Z			'Z'
#define TK_a			'a'
#define TK_b			'b'
#define TK_c			'c'
#define TK_d			'd'
#define TK_e			'e'
#define TK_f			'f'
#define TK_g			'g'
#define TK_h			'h'
#define TK_i			'i'
#define TK_j			'j'
#define TK_k			'k'
#define TK_l			'l'
#define TK_m			'm'
#define TK_n			'n'
#define TK_o			'o'
#define TK_p			'p'
#define TK_q			'q'
#define TK_r			'r'
#define TK_s			's'
#define TK_t			't'
#define TK_u			'u'
#define TK_v			'v'
#define TK_w			'w'
#define TK_x			'x'
#define TK_y			'y'
#define TK_z			'z'
#define TK_0			'0'
#define TK_1			'1'
#define TK_2			'2'
#define TK_3			'3'
#define TK_4			'4'
#define TK_5			'5'
#define TK_6			'6'
#define TK_7			'7'
#define TK_8			'8'
#define TK_9			'9'

/*
** Color Macros
*/

enum {
    TK_BLACK = 0,
#ifdef __WIN32__
    TK_RED=10,
#else
    TK_RED,
#endif
    TK_GREEN,
    TK_YELLOW,
    TK_BLUE,
    TK_MAGENTA,
    TK_CYAN,
    TK_WHITE
};

#ifdef __WIN32__
extern float tkRGBMap[17][3];
#else
extern float tkRGBMap[8][3];
#endif

#define TK_SETCOLOR(x, y) (TK_IS_RGB((x)) ? \
		           glColor3fv(tkRGBMap[(y)]) : glIndexf((y)))

/*
** RGB Image Structure
*/

typedef struct _TK_RGBImageRec {
    GLint sizeX, sizeY;
    unsigned char *data;
} TK_RGBImageRec;

/*
** Prototypes
*/

extern void tkInitDisplayMode(GLenum);
extern void tkInitPosition(int, int, int, int);
extern GLenum tkInitWindow(char *);
extern void tkCloseWindow(void);
extern void tkQuit(void);

extern GLenum tkSetWindowLevel(GLenum);
extern void tkSwapBuffers(void);

extern void tkExec( int );
extern void tkExposeFunc(void (*)(int, int));
extern void tkReshapeFunc(void (*)(int, int));
extern void tkDisplayFunc(void (*)(void));
extern void tkKeyDownFunc(GLenum (*)(int, GLenum));
extern void tkMouseDownFunc(GLenum (*)(int, int, GLenum));
extern void tkMouseUpFunc(GLenum (*)(int, int, GLenum));
extern void tkMouseMoveFunc(GLenum (*)(int, int, GLenum));
extern void tkIdleFunc(void (*)(void));

extern int tkGetColorMapSize(void);
extern void tkGetMouseLoc(int *, int *);
extern void tkGetSystem(TKenum, void *);

extern void tkSetOneColor(int, float, float, float);
extern void tkSetFogRamp(int, int);
extern void tkSetGreyRamp(void);
extern void tkSetRGBMap(int, float *);
extern void tkSetOverlayMap(int, float *);

extern void tkNewCursor(GLint, GLubyte *, GLubyte *, GLenum, GLenum,
			GLint, GLint);
extern void tkSetCursor(GLint);

extern TK_RGBImageRec *tkRGBImageLoad(char *);

extern GLenum tkCreateStrokeFont(GLuint);
extern GLenum tkCreateOutlineFont(GLuint);
extern GLenum tkCreateFilledFont(GLuint);
extern GLenum tkCreateBitmapFont(GLuint);
extern void tkDrawStr(GLuint, char *);

extern void tkWireSphere(GLuint, float);
extern void tkSolidSphere(GLuint, float);
extern void tkWireCube(GLuint, float);
extern void tkSolidCube(GLuint, float);
extern void tkWireBox(GLuint, float, float, float);
extern void tkSolidBox(GLuint, float, float, float);
extern void tkWireTorus(GLuint, float, float);
extern void tkSolidTorus(GLuint, float, float);
extern void tkWireCylinder(GLuint, float, float);
extern void tkSolidCylinder(GLuint, float, float);
extern void tkWireCone(GLuint, float, float);
extern void tkSolidCone(GLuint, float, float);


#ifdef __cplusplus
}
#endif


#ifdef __WIN32__
#include <windows.h>
HWND tkGetHWND(void);
#endif /* WIN32 */


#endif
