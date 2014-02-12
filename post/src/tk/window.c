#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tk.h"
#include "private.h"

/******************************************************************************/

Display *xDisplay = 0;
int xScreen = 0; 
Window wRoot = 0;
Atom deleteWindowAtom;
WINDOW_REC w = {
    0, 0, 300, 300, TK_RGB|TK_SINGLE|TK_DIRECT
};
float colorMaps[] = {
    0.000000, 1.000000, 0.000000, 1.000000, 0.000000, 1.000000, 
    0.000000, 1.000000, 0.333333, 0.776471, 0.443137, 0.556863, 
    0.443137, 0.556863, 0.219608, 0.666667, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.039216, 0.078431, 0.117647, 0.156863, 
    0.200000, 0.239216, 0.278431, 0.317647, 0.356863, 0.400000, 
    0.439216, 0.478431, 0.517647, 0.556863, 0.600000, 0.639216, 
    0.678431, 0.717647, 0.756863, 0.800000, 0.839216, 0.878431, 
    0.917647, 0.956863, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 0.000000, 0.000000, 
    1.000000, 1.000000, 0.000000, 0.000000, 1.000000, 1.000000, 
    0.333333, 0.443137, 0.776471, 0.556863, 0.443137, 0.219608, 
    0.556863, 0.666667, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.039216, 0.078431, 0.117647, 0.156863, 0.200000, 0.239216, 
    0.278431, 0.317647, 0.356863, 0.400000, 0.439216, 0.478431, 
    0.517647, 0.556863, 0.600000, 0.639216, 0.678431, 0.717647, 
    0.756863, 0.800000, 0.839216, 0.878431, 0.917647, 0.956863, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 0.333333, 0.443137, 
    0.443137, 0.219608, 0.776471, 0.556863, 0.556863, 0.666667, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.039216, 0.078431, 
    0.117647, 0.156863, 0.200000, 0.239216, 0.278431, 0.317647, 
    0.356863, 0.400000, 0.439216, 0.478431, 0.517647, 0.556863, 
    0.600000, 0.639216, 0.678431, 0.717647, 0.756863, 0.800000, 
    0.839216, 0.878431, 0.917647, 0.956863, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
};
float tkRGBMap[8][3] = {
    {
	0, 0, 0
    },
    {
	1, 0, 0
    },
    {
	0, 1, 0
    },
    {
	1, 1, 0
    },
    {
	0, 0, 1
    },
    {
	1, 0, 1
    },
    {
	0, 1, 1
    },
    {
	1, 1, 1
    }
};

/******************************************************************************/

void tkCloseWindow(void)
{

    if (xDisplay) {
	cursorNum = 0;

	ExposeFunc = 0;
	ReshapeFunc = 0;
	IdleFunc = 0;
	DisplayFunc = 0;
	KeyDownFunc = 0;
	MouseDownFunc = 0;
	MouseUpFunc = 0;
	MouseMoveFunc = 0;

	glFlush();
	glFinish();
	if (TK_HAS_OVERLAY(w.type)) {
	    XDestroyWindow(xDisplay, w.wOverlay);
	    glXDestroyContext(xDisplay, w.cOverlay);
	    XFreeColormap(xDisplay, w.cMapOverlay);
	    XFree((char *)w.vInfoOverlay);
	}
	XDestroyWindow(xDisplay, w.wMain);
	glXDestroyContext(xDisplay, w.cMain);
        if (w.cMapAllocated)
            XFreeColormap(xDisplay, w.cMapMain);
	XFree((char *)w.vInfoMain);
	XCloseDisplay(xDisplay);
	xDisplay = 0;
    }
}

/******************************************************************************/

void tkInitDisplayMode(GLenum type)
{

    w.type = type;
}

/******************************************************************************/

void tkInitPosition(int x, int y, int width, int height)
{

    w.x = x;
    w.y = y;
    w.w = width;
    w.h = height;
}

/******************************************************************************/

static int ErrorHandler(Display *xDisplay, XErrorEvent *event)
{
    char buf[80];

    printf("\nReceived X error!\n");
    printf("\tError code   : %d\n", event->error_code);
    printf("\tRequest code : %d\n", event->request_code);
    printf("\tMinor code   : %d\n\n", event->minor_code);
    XGetErrorText(xDisplay, event->error_code, buf, 80);
    printf("\tError text : '%s'\n\n", buf);
    return 0;
}

static XVisualInfo *FindMainVisual(GLenum type)
{
    int list[32], i;

    i = 0;

    list[i++] = GLX_LEVEL;
    list[i++] = 0;

    if (TK_IS_DOUBLE(type)) {
	list[i++] = GLX_DOUBLEBUFFER;
    }

    if (TK_IS_RGB(type)) {
	list[i++] = GLX_RGBA;
	list[i++] = GLX_RED_SIZE;
	list[i++] = 1;
	list[i++] = GLX_GREEN_SIZE;
	list[i++] = 1;
	list[i++] = GLX_BLUE_SIZE;
	list[i++] = 1;
	if (TK_HAS_ALPHA(type)) {
	    list[i++] = GLX_ALPHA_SIZE;
	    list[i++] = 1;
	}
	if (TK_HAS_ACCUM(type)) {
	    list[i++] = GLX_ACCUM_RED_SIZE;
	    list[i++] = 1;
	    list[i++] = GLX_ACCUM_GREEN_SIZE;
	    list[i++] = 1;
	    list[i++] = GLX_ACCUM_BLUE_SIZE;
	    list[i++] = 1;
	    if (TK_HAS_ALPHA(type)) {
		list[i++] = GLX_ACCUM_ALPHA_SIZE;
		list[i++] = 1;
	    }
	}
    } else if (TK_IS_INDEX(type)) {
	list[i++] = GLX_BUFFER_SIZE;
	list[i++] = 1;
    }

    if (TK_HAS_DEPTH(type)) {
	list[i++] = GLX_DEPTH_SIZE;
	list[i++] = 1;
    }

    if (TK_HAS_STENCIL(type)) {
	list[i++] = GLX_STENCIL_SIZE;
	list[i++] = 1;
    }

    list[i] = (int)None;

    return glXChooseVisual(xDisplay, xScreen, list);
}

static XVisualInfo *FindOverlayVisual(void)
{
    int list[3];

    list[0] = GLX_LEVEL;
    list[1] = 1;
    list[2] = (int)None;

    return glXChooseVisual(xDisplay, xScreen, list);
}

static GLenum GetMainWindowType(XVisualInfo *vi)
{
    GLenum mask;
    int x, y, z;

    mask = 0;

    glXGetConfig(xDisplay, vi, GLX_DOUBLEBUFFER, &x);
    if (x) {
	mask |= TK_DOUBLE;
    } else {
	mask |= TK_SINGLE;
    }

    glXGetConfig(xDisplay, vi, GLX_RGBA, &x);
    if (x) {
	mask |= TK_RGB;
	glXGetConfig(xDisplay, vi, GLX_ALPHA_SIZE, &x);
	if (x > 0) {
	    mask |= TK_ALPHA;
	}
	glXGetConfig(xDisplay, vi, GLX_ACCUM_RED_SIZE, &x);
	glXGetConfig(xDisplay, vi, GLX_ACCUM_GREEN_SIZE, &y);
	glXGetConfig(xDisplay, vi, GLX_ACCUM_BLUE_SIZE, &z);
	if (x > 0 && y > 0 && z > 0) {
	    mask |= TK_ACCUM;
	}
    } else {
	mask |= TK_INDEX;
    }

    glXGetConfig(xDisplay, vi, GLX_DEPTH_SIZE, &x);
    if (x > 0) {
	mask |= TK_DEPTH;
    }

    glXGetConfig(xDisplay, vi, GLX_STENCIL_SIZE, &x);
    if (x > 0) {
	mask |= TK_STENCIL;
    }

    if (glXIsDirect(xDisplay, w.cMain)) {
	mask |= TK_DIRECT;
    } else {
	mask |= TK_INDIRECT;
    }

    return mask;
}

static int WaitForMainWindow(Display *d, XEvent *e, char *arg)
{

    if (e->type == MapNotify && e->xmap.window == w.wMain) {
	return GL_TRUE;
    } else {
	return GL_FALSE;
    }
}

static int WaitForOverlayWindow(Display *d, XEvent *e, char *arg)
{

    if (e->type == MapNotify && e->xmap.window == w.wOverlay) {
	return GL_TRUE;
    } else {
	return GL_FALSE;
    }
}

GLenum tkInitWindow(char *title)
{
    XSetWindowAttributes wa;
    XTextProperty tp;
    XSizeHints sh;
    XEvent e;
    int erb, evb;
    GLenum overlayFlag;

    if (!xDisplay) {
	xDisplay = XOpenDisplay(0);
	if (!xDisplay) {
	    fprintf(stderr, "Can't connect to xDisplay!\n");
	    return GL_FALSE;
	}
	if (!glXQueryExtension(xDisplay, &erb, &evb)) {
	    fprintf(stderr, "No glx extension!\n");
	    return GL_FALSE;
	}
	xScreen = DefaultScreen(xDisplay);
	wRoot = RootWindow(xDisplay, xScreen);
	XSetErrorHandler(ErrorHandler);
    }

    if (TK_HAS_OVERLAY(w.type)) {
	overlayFlag = GL_TRUE;
    } else {
	overlayFlag = GL_FALSE;
    }
    w.type &= ~TK_OVERLAY;

    w.vInfoMain = FindMainVisual(w.type);
    if (!w.vInfoMain) {
       if (TK_IS_RGB(w.type)) {
	  fprintf(stderr, "Couldn't find visual for RGB mode!\n");
       }
       else {
	  fprintf(stderr, "Couldn't find visual for Color Index mode!\n");
       }
       xDisplay = 0;  /* to prevent crash in tkCloseWindow */
       return GL_FALSE;
    }

    w.cMain = glXCreateContext(xDisplay, w.vInfoMain, NULL,
			       (TK_IS_DIRECT(w.type))?GL_TRUE:GL_FALSE);
    if (!w.cMain) {
	fprintf(stderr, "Can't create a context!\n");
	return GL_FALSE;
    }

    w.type = GetMainWindowType(w.vInfoMain);

    w.cMapAllocated = 1;
    if (TK_IS_INDEX(w.type)) {
        /* Color Indexed windows needs a writable colormap */
	if (w.vInfoMain->class != StaticColor &&
	    w.vInfoMain->class != StaticGray) {
	    w.cMapMain = XCreateColormap(xDisplay, wRoot, w.vInfoMain->visual,
				         AllocAll);
	} else {
	    w.cMapMain = XCreateColormap(xDisplay, wRoot, w.vInfoMain->visual,
				         AllocNone);
	}
    } else {
        /* RGB colormap is AllocNone, share the root colormap if possible */
	Screen *scr = DefaultScreenOfDisplay(xDisplay);
	int scrnum = DefaultScreen(xDisplay);
	if (MaxCmapsOfScreen(scr)==1
	    && w.vInfoMain->visual==DefaultVisual(xDisplay,scrnum)) {
	   /* the window and root are of the same visual type */
	   char *private = getenv("MESA_PRIVATE_CMAP");
	   if (private) {
	      /* user doesn't want to share colormaps */
	      w.cMapMain = XCreateColormap(xDisplay, wRoot,
					   w.vInfoMain->visual, AllocNone);
	   }
	   else {
	      /* share the root colormap */
	      w.cMapMain = DefaultColormap(xDisplay,scrnum);
	   }
	}
	else {
	   /* window and root are different visual types, allocate new cmap */
/****  ad@lms.be:  Added initialization for HP color recovery  *******/
           Atom hp_cr_maps = XInternAtom(xDisplay, "_HP_RGB_SMOOTH_MAP_LIST", True);
           w.cMapMain = 0;
           if (hp_cr_maps) {
                XStandardColormap* colmaps = 0;
                int nrColmaps = 0;
                int i;
                XGetRGBColormaps( xDisplay, RootWindow(xDisplay, scrnum)
                                , &colmaps, &nrColmaps, hp_cr_maps);
                for (i=0; i<nrColmaps; i++) {
                    if (colmaps[i].visualid == w.vInfoMain->visual->visualid) {
                        w.cMapMain = colmaps[i].colormap;
                        w.cMapAllocated = 0;
                        break;
                    }
                }
           }  /*** end HP color recovery ***/
           if (!w.cMapMain) {
               if (w.vInfoMain->class==DirectColor) {
                  w.cMapMain = XCreateColormap(xDisplay, wRoot,
                                               w.vInfoMain->visual, AllocAll);
               }
               else {
                  w.cMapMain = XCreateColormap(xDisplay, wRoot,
                                               w.vInfoMain->visual, AllocNone);
               }
           }
	}
    }
    if (TK_IS_INDEX(w.type) || w.vInfoMain->class==DirectColor) {
       tkSetRGBMap(256, colorMaps);
    }
    wa.colormap = w.cMapMain;
    wa.background_pixmap = None;
    wa.border_pixel = 0;
    wa.event_mask = StructureNotifyMask | ExposureMask | KeyPressMask |
		    ButtonPressMask | ButtonReleaseMask | PointerMotionMask;
    w.wMain = XCreateWindow(xDisplay, wRoot, w.x, w.y, w.w, w.h, 0,
			    w.vInfoMain->depth, InputOutput,
			    w.vInfoMain->visual,
			    CWBackPixmap|CWBorderPixel|CWEventMask|CWColormap,
			    &wa);

    /*OLD: XInstallColormap( xDisplay, w.cMapMain );*/
    XSetWMColormapWindows( xDisplay, w.wMain, &w.wMain, 1 );

    XStringListToTextProperty(&title, 1, &tp);
    sh.flags = USPosition | USSize;
    XSetWMProperties(xDisplay, w.wMain, &tp, &tp, 0, 0, &sh, 0, 0);
    deleteWindowAtom = XInternAtom(xDisplay, "WM_DELETE_WINDOW", False);
    XSetWMProtocols(xDisplay, w.wMain, &deleteWindowAtom, 1);
    XMapWindow(xDisplay, w.wMain);
    drawAllowFlag = GL_FALSE;
    XIfEvent(xDisplay, &e, WaitForMainWindow, 0);

    if (overlayFlag == GL_TRUE) {
	w.vInfoOverlay = FindOverlayVisual();
	if (w.vInfoOverlay) {
	    w.cOverlay = glXCreateContext(xDisplay, w.vInfoOverlay, None,
					  GL_TRUE);
	    w.cMapOverlay = XCreateColormap(xDisplay, wRoot,
					    w.vInfoOverlay->visual, AllocNone);
	    tkSetOverlayMap(256, colorMaps);
	    wa.colormap = w.cMapOverlay;
	    wa.background_pixmap = None;
	    wa.border_pixel = 0;
	    w.wOverlay = XCreateWindow(xDisplay, w.wMain, 0, 0, w.w, w.h, 0,
				       w.vInfoOverlay->depth, InputOutput,
				       w.vInfoOverlay->visual,
				       CWBackPixmap|CWBorderPixel|CWColormap,
				       &wa);
	    XMapWindow(xDisplay, w.wOverlay);
	    XSetWMColormapWindows(xDisplay, w.wMain, &w.wOverlay, 1);
	    w.type |= TK_OVERLAY;
	} else {
	    fprintf(stderr, "Can't create a overlay plane!\n");
	}
    }

    if (!glXMakeCurrent(xDisplay, w.wMain, w.cMain)) {
	fprintf(stderr, "Can't make window current drawable!\n");
	return GL_FALSE;
    }
    XFlush(xDisplay);

    return GL_TRUE;
}

/******************************************************************************/

void tkQuit(void)
{

    tkCloseWindow();
    exit(0);
}

/******************************************************************************/

void tkSwapBuffers(void)
{

    if (xDisplay) {
	glXSwapBuffers(xDisplay, w.wMain);
    }
}

/******************************************************************************/
