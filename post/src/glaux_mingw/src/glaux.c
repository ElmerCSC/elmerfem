/* 02/21/03 Modified by Greg Chien <gchien@protodesign-inc.com>
 * - Add APIENTRY to be compatible with glaux.h distributed by Microsoft.
 * - Cast function pointers in the function call arguments and RHS of
 *   assignments to (void*) to fix type mismatch errors.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <windows.h>
#include <GL/gl.h>
#include <gl/glaux.h>
#include "tk.h"

#define static

#if defined(__cplusplus) || defined(c_plusplus)
#define class c_class
#endif


static struct {
    int keyField;
    void (*KeyFunc)(void);
} keyTable[200];

static struct {
    int mouseField;
    void (*MouseFunc)(AUX_EVENTREC *);
} mouseDownTable[20], mouseUpTable[20], mouseLocTable[20];

static int keyTableCount = 0;
static int mouseDownTableCount = 0;
static int mouseUpTableCount = 0;
static int mouseLocTableCount = 0;
static GLenum displayModeType = 0;
GLenum APIENTRY auxInitWindowAW(LPCSTR title, BOOL bUnicode);

static void DefaultHandleReshape(GLsizei w, GLsizei h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho((GLdouble)0.0, (GLdouble)w, (GLdouble)0.0, (GLdouble)h, (GLdouble)-1.0, (GLdouble)1.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

static void DefaultHandleExpose(int w, int h)
{
}

static GLenum MouseLoc(int x, int y, GLenum button)
{
    AUX_EVENTREC info;
    GLenum flag;
    int i;

    flag = GL_FALSE;
    for (i = 0; i < mouseLocTableCount; i++) {
        if ((int)(button & AUX_LEFTBUTTON) == mouseLocTable[i].mouseField) {
	    info.event = AUX_MOUSELOC;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[AUX_MOUSESTATUS] = AUX_LEFTBUTTON;
	    (*mouseLocTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
        if ((int)(button & AUX_RIGHTBUTTON) == mouseLocTable[i].mouseField) {
	    info.event = AUX_MOUSELOC;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[AUX_MOUSESTATUS] = AUX_RIGHTBUTTON;
	    (*mouseLocTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
        if ((int)(button & AUX_MIDDLEBUTTON) == mouseLocTable[i].mouseField) {
	    info.event = AUX_MOUSELOC;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[AUX_MOUSESTATUS] = AUX_MIDDLEBUTTON;
	    (*mouseLocTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
    }
    return flag;
}

static GLenum MouseUp(int x, int y, GLenum button)
{
    AUX_EVENTREC info;
    GLenum flag;
    int i;

    flag = GL_FALSE;
    for (i = 0; i < mouseUpTableCount; i++) {
        if ((int)(button & AUX_LEFTBUTTON) == mouseUpTable[i].mouseField) {
	    info.event = AUX_MOUSEUP;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[AUX_MOUSESTATUS] = AUX_LEFTBUTTON;
	    (*mouseUpTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
        if ((int)(button & AUX_RIGHTBUTTON) == mouseUpTable[i].mouseField) {
	    info.event = AUX_MOUSEUP;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[AUX_MOUSESTATUS] = AUX_RIGHTBUTTON;
	    (*mouseUpTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
        if ((int)(button & AUX_MIDDLEBUTTON) == mouseUpTable[i].mouseField) {
	    info.event = AUX_MOUSEUP;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[AUX_MOUSESTATUS] = AUX_MIDDLEBUTTON;
	    (*mouseUpTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
    }
    return flag;
}

static GLenum MouseDown(int x, int y, GLenum button)
{
    AUX_EVENTREC info;
    GLenum flag;
    int i;

    flag = GL_FALSE;
    for (i = 0; i < mouseDownTableCount; i++) {
        if ((int)(button & AUX_LEFTBUTTON) == mouseDownTable[i].mouseField) {
	    info.event = AUX_MOUSEDOWN;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[AUX_MOUSESTATUS] = AUX_LEFTBUTTON;
	    (*mouseDownTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
        if ((int)(button & AUX_RIGHTBUTTON) == mouseDownTable[i].mouseField) {
	    info.event = AUX_MOUSEDOWN;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[AUX_MOUSESTATUS] = AUX_RIGHTBUTTON;
	    (*mouseDownTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
        if ((int)(button & AUX_MIDDLEBUTTON) == mouseDownTable[i].mouseField) {
	    info.event = AUX_MOUSEDOWN;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[AUX_MOUSESTATUS] = AUX_MIDDLEBUTTON;
	    (*mouseDownTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
    }
    return flag;
}

static GLenum KeyDown(int key, GLenum status)
{
    GLenum flag;
    int i;

    flag = GL_FALSE;
    if (keyTableCount) {
	for (i = 0; i < keyTableCount; i++) {
	    if (key == keyTable[i].keyField) {
		(*keyTable[i].KeyFunc)();
		flag |= GL_TRUE;
	    }
	}
    }
    return flag;
}

void APIENTRY auxExposeFunc(AUXEXPOSEPROC Func)
{
    tkExposeFunc((void*)Func);
}

void APIENTRY auxReshapeFunc(AUXRESHAPEPROC Func)
{
    tkExposeFunc((void*)Func);
    tkReshapeFunc((void*)Func);
}

void APIENTRY auxIdleFunc(AUXIDLEPROC Func)
{
    tkIdleFunc((void*)Func);
}

void APIENTRY auxKeyFunc(int key, AUXKEYPROC Func)
{
    keyTable[keyTableCount].keyField = key;
    keyTable[keyTableCount++].KeyFunc = (void*)Func;
}

void APIENTRY auxMouseFunc(int mouse, int mode, AUXMOUSEPROC Func)
{
    if (mode == AUX_MOUSEDOWN) {
	mouseDownTable[mouseDownTableCount].mouseField = mouse;
	mouseDownTable[mouseDownTableCount++].MouseFunc = (void*)Func;
    } else if (mode == AUX_MOUSEUP) {
	mouseUpTable[mouseUpTableCount].mouseField = mouse;
	mouseUpTable[mouseUpTableCount++].MouseFunc = (void*)Func;
    } else if (mode == AUX_MOUSELOC) {
	mouseLocTable[mouseLocTableCount].mouseField = mouse;
	mouseLocTable[mouseLocTableCount++].MouseFunc = (void*)Func;
    } 
}

void APIENTRY auxMainLoop(AUXMAINPROC Func)
{
    tkDisplayFunc((void*)Func);
    tkExec();
}

void APIENTRY auxInitPosition(int x, int y, int width, int height)
{
    tkInitPosition(x, y, width, height);
}

void APIENTRY auxInitDisplayMode(GLenum type)
{
    displayModeType = type;
    tkInitDisplayMode(type);
}

void APIENTRY auxInitDisplayModePolicy(GLenum type)
{
    tkInitDisplayModePolicy(type);
}

GLenum APIENTRY auxInitDisplayModeID(GLint id)
{
    return tkInitDisplayModeID(id);
}

GLenum APIENTRY auxInitWindowA(LPCSTR title)
{
    return auxInitWindowAW(title,FALSE);
}

GLenum APIENTRY auxInitWindowW(LPCWSTR title)
{
    return auxInitWindowAW((LPCSTR)title,TRUE);
}

GLenum APIENTRY auxInitWindowAW(LPCSTR title, BOOL bUnicode)
{
    int useDoubleAsSingle = 0;

    if (tkInitWindowAW((char *)title, bUnicode) == GL_FALSE) {
	if (AUX_WIND_IS_SINGLE(displayModeType)) {
	    tkInitDisplayMode(displayModeType | AUX_DOUBLE);
	    if (tkInitWindowAW((char *)title, bUnicode) == GL_FALSE) {
		return GL_FALSE;    /*  curses, foiled again	*/
            }
            MESSAGEBOX(GetFocus(), "Can't initialize a single buffer visual. "
                                   "Will use a double buffer visual instead, "
                                   "only drawing into the front buffer.",
                                   "Warning", MB_OK);
	    displayModeType = displayModeType | AUX_DOUBLE;
	    useDoubleAsSingle = 1;
	}
    }
    tkReshapeFunc(DefaultHandleReshape);
    tkExposeFunc(DefaultHandleExpose);
    tkMouseUpFunc(MouseUp);
    tkMouseDownFunc(MouseDown);
    tkMouseMoveFunc(MouseLoc);
    tkKeyDownFunc(KeyDown);
    auxKeyFunc(AUX_ESCAPE, auxQuit);
    glClearColor((GLclampf)0.0, (GLclampf)0.0, (GLclampf)0.0, (GLclampf)1.0);
    glClearIndex((GLfloat)0.0);
    glLoadIdentity();
    if (useDoubleAsSingle)
	glDrawBuffer(GL_FRONT);
    return GL_TRUE;
}

void APIENTRY auxCloseWindow(void)
{
    tkCloseWindow();
    keyTableCount = 0;
    mouseDownTableCount = 0;
    mouseUpTableCount = 0;
    mouseLocTableCount = 0;
}

void APIENTRY auxQuit(void)
{
    tkQuit();
}

void APIENTRY auxSwapBuffers(void)
{
    tkSwapBuffers();
}

HWND APIENTRY auxGetHWND(void)
{
    return tkGetHWND();
}

HDC APIENTRY auxGetHDC(void)
{
    return tkGetHDC();
}

HGLRC APIENTRY auxGetHGLRC(void)
{
    return tkGetHRC();
}

GLenum APIENTRY auxGetDisplayModePolicy(void)
{
    return tkGetDisplayModePolicy();
}

GLint APIENTRY auxGetDisplayModeID(void)
{
    return tkGetDisplayModeID();
}

GLenum APIENTRY auxGetDisplayMode(void)
{
    return tkGetDisplayMode();
}

void APIENTRY auxSetOneColor(int index, float r, float g, float b)
{
    tkSetOneColor(index, r, g, b);
}

void APIENTRY auxSetFogRamp(int density, int startIndex)
{
    tkSetFogRamp(density, startIndex);
}

void APIENTRY auxSetGreyRamp(void)
{
    tkSetGreyRamp();
}

void APIENTRY auxSetRGBMap(int size, float *rgb)
{
    tkSetRGBMap(size, rgb);
}

int APIENTRY auxGetColorMapSize(void)
{
    return tkGetColorMapSize();
}

void APIENTRY auxGetMouseLoc(int *x, int *y)
{
    tkGetMouseLoc(x, y);
}
