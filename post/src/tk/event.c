#include <stdio.h>
#include <stdlib.h>
#include <X11/keysym.h>
#include "tk.h"
#include "private.h"

/******************************************************************************/

void (*ExposeFunc)(int, int) = 0;
void (*ReshapeFunc)(int, int) = 0;
void (*DisplayFunc)(void) = 0;
GLenum (*KeyDownFunc)(int, GLenum) = 0;
GLenum (*MouseDownFunc)(int, int, GLenum) = 0;
GLenum (*MouseUpFunc)(int, int, GLenum) = 0;
GLenum (*MouseMoveFunc)(int, int, GLenum) = 0;
void (*IdleFunc)(void) = 0;
int lastEventType = -1;
GLenum drawAllowFlag;

/******************************************************************************/

static GLenum DoNextEvent(void)
{
    XEvent current, ahead;
    char buf[1000];
    KeySym ks;
    int key;

    XNextEvent(xDisplay, &current);
    switch (current.type) {
      case MappingNotify:
	XRefreshKeyboardMapping((XMappingEvent *)&current);
	lastEventType = MappingNotify;
	return GL_FALSE;

      case MapNotify:
	lastEventType = MapNotify;
	drawAllowFlag = GL_TRUE;
	return GL_FALSE;

      case UnmapNotify:
	lastEventType = UnmapNotify;
	drawAllowFlag = GL_FALSE;
	return GL_FALSE;
      
      case ClientMessage:
	lastEventType = ClientMessage;
	if (current.xclient.data.l[0] == deleteWindowAtom) {
	    exit(0);
	}
	return GL_FALSE;

      case Expose:
	while (XEventsQueued(current.xexpose.display, QueuedAfterReading) > 0) {
	    XPeekEvent(current.xexpose.display, &ahead);
	    if (ahead.xexpose.window != current.xexpose.window ||
		ahead.type != Expose) {
		break;
	    }
	    XNextEvent(xDisplay, &current);
	}
	if (current.xexpose.count == 0) {
	    if (ExposeFunc) {
		(*ExposeFunc)(w.w, w.h);
		if (lastEventType == ConfigureNotify &&
		    drawAllowFlag == GL_TRUE) {
		    lastEventType = Expose;
		    return GL_FALSE;
		} else {
		    lastEventType = Expose;
		    drawAllowFlag = GL_TRUE;
		    return GL_TRUE;
		}
	    }
	}
	return GL_FALSE;

      case ConfigureNotify:
	lastEventType = ConfigureNotify;
	w.w = current.xconfigure.width;
	w.h = current.xconfigure.height;
	if (TK_HAS_OVERLAY(w.type)) {
	    XResizeWindow(xDisplay, w.wOverlay, w.w, w.h);
	}
	if (ReshapeFunc) {
	    (*ReshapeFunc)(w.w, w.h);
	    return GL_TRUE;
	} else {
	    return GL_FALSE;
	}

      case MotionNotify:
	lastEventType = MotionNotify;
	if (MouseMoveFunc) {
	    GLenum mask;

	    mask = 0;
	    if (current.xmotion.state & Button1Mask) {
		mask |= TK_LEFTBUTTON;
	    }
	    if (current.xmotion.state & Button2Mask) {
		mask |= TK_MIDDLEBUTTON;
	    }
	    if (current.xmotion.state & Button3Mask) {
		mask |= TK_RIGHTBUTTON;
	    }
	    return (*MouseMoveFunc)(current.xmotion.x, current.xmotion.y, mask);
	} else {
	    return GL_FALSE;
	}

      case ButtonPress:
	lastEventType = ButtonPress;
	if (MouseDownFunc) {
	    GLenum mask;

	    mask = 0;
	    if (current.xbutton.button == 1) {
		mask |= TK_LEFTBUTTON;
	    }
	    if (current.xbutton.button == 2) {
		mask |= TK_MIDDLEBUTTON;
	    }
	    if (current.xbutton.button == 3) {
		mask |= TK_RIGHTBUTTON;
	    }
	    return (*MouseDownFunc)(current.xbutton.x, current.xbutton.y, mask);
	} else {
	    return GL_FALSE;
	}
      case ButtonRelease:
	lastEventType = ButtonRelease;
	if (MouseUpFunc) {
	    GLenum mask;

	    mask = 0;
	    if (current.xbutton.button == 1) {
		mask |= TK_LEFTBUTTON;
	    }
	    if (current.xbutton.button == 2) {
		mask |= TK_MIDDLEBUTTON;
	    }
	    if (current.xbutton.button == 3) {
		mask |= TK_RIGHTBUTTON;
	    }
	    return (*MouseUpFunc)(current.xbutton.x, current.xbutton.y, mask);
	} else {
	    return GL_FALSE;
	}

      case KeyPress:
	lastEventType = KeyPress;
	XLookupString(&current.xkey, buf, sizeof(buf), &ks, 0);
	switch (ks) {
	  case XK_0: 		key = TK_0;		break;
	  case XK_1: 		key = TK_1;		break;
	  case XK_2: 		key = TK_2;		break;
	  case XK_3: 		key = TK_3;		break;
	  case XK_4: 		key = TK_4;		break;
	  case XK_5: 		key = TK_5;		break;
	  case XK_6: 		key = TK_6;		break;
	  case XK_7: 		key = TK_7;		break;
	  case XK_8: 		key = TK_8;		break;
	  case XK_9: 		key = TK_9;		break;
	  case XK_A: 		key = TK_A;		break;
	  case XK_B: 		key = TK_B;		break;
	  case XK_C: 		key = TK_C;		break;
	  case XK_D: 		key = TK_D;		break;
	  case XK_E: 		key = TK_E;		break;
	  case XK_F: 		key = TK_F;		break;
	  case XK_G: 		key = TK_G;		break;
	  case XK_H: 		key = TK_H;		break;
	  case XK_I: 		key = TK_I;		break;
	  case XK_J: 		key = TK_J;		break;
	  case XK_K: 		key = TK_K;		break;
	  case XK_L: 		key = TK_L;		break;
	  case XK_M: 		key = TK_M;		break;
	  case XK_N: 		key = TK_N;		break;
	  case XK_O: 		key = TK_O;		break;
	  case XK_P: 		key = TK_P;		break;
	  case XK_Q: 		key = TK_Q;		break;
	  case XK_R: 		key = TK_R;		break;
	  case XK_S: 		key = TK_S;		break;
	  case XK_T: 		key = TK_T;		break;
	  case XK_U: 		key = TK_U;		break;
	  case XK_V: 		key = TK_V;		break;
	  case XK_W: 		key = TK_W;		break;
	  case XK_X: 		key = TK_X;		break;
	  case XK_Y: 		key = TK_Y;		break;
	  case XK_Z: 		key = TK_Z;		break;
	  case XK_a: 		key = TK_a;		break;
	  case XK_b: 		key = TK_b;		break;
	  case XK_c: 		key = TK_c;		break;
	  case XK_d: 		key = TK_d;		break;
	  case XK_e: 		key = TK_e;		break;
	  case XK_f: 		key = TK_f;		break;
	  case XK_g: 		key = TK_g;		break;
	  case XK_h: 		key = TK_h;		break;
	  case XK_i: 		key = TK_i;		break;
	  case XK_j: 		key = TK_j;		break;
	  case XK_k: 		key = TK_k;		break;
	  case XK_l: 		key = TK_l;		break;
	  case XK_m: 		key = TK_m;		break;
	  case XK_n: 		key = TK_n; 		break;
	  case XK_o: 		key = TK_o;		break;
	  case XK_p: 		key = TK_p;		break;
	  case XK_q: 		key = TK_q;		break;
	  case XK_r: 		key = TK_r;		break;
	  case XK_s: 		key = TK_s;		break;
	  case XK_t: 		key = TK_t;		break;
	  case XK_u: 		key = TK_u;		break;
	  case XK_v: 		key = TK_v;		break;
	  case XK_w: 		key = TK_w;		break;
	  case XK_x: 		key = TK_x;		break;
	  case XK_y: 		key = TK_y;		break;
	  case XK_z: 		key = TK_z;		break;
	  case XK_space:	key = TK_SPACE;		break;
	  case XK_Return: 	key = TK_RETURN;	break;
	  case XK_Escape: 	key = TK_ESCAPE;	break;
	  case XK_Left:		key = TK_LEFT;		break;
	  case XK_Up:		key = TK_UP;		break;
	  case XK_Right:  	key = TK_RIGHT;		break;
	  case XK_Down:		key = TK_DOWN;		break;
	  default: 		key = GL_FALSE;		break;
	}
	if (key && KeyDownFunc) {
	    GLenum mask;

	    mask = 0;
	    if (current.xkey.state & ControlMask) {
		mask |= TK_CONTROL;
	    }
	    if (current.xkey.state & ShiftMask) {
		mask |= TK_SHIFT;
	    }
	    return (*KeyDownFunc)(key, mask);
	} else {
	    return GL_FALSE;
	}
    }
    return GL_FALSE;
}

void tkExec(int i)
{
    GLenum flag;

    while (GL_TRUE) {
	if (IdleFunc) {
	    (*IdleFunc)();
	    flag = GL_TRUE;
	    while /*if*/ (XPending(xDisplay)) {
		flag |= DoNextEvent();
	    }
	} else {
	    flag = DoNextEvent();
	}
	if (drawAllowFlag && DisplayFunc && flag) {
	    (*DisplayFunc)();
	}
    }
}

/******************************************************************************/

void tkExposeFunc(void (*Func)(int, int))
{

    ExposeFunc = Func;
}

/******************************************************************************/

void tkReshapeFunc(void (*Func)(int, int))
{

    ReshapeFunc = Func;
}

/******************************************************************************/

void tkDisplayFunc(void (*Func)(void))
{

    DisplayFunc = Func;
}

/******************************************************************************/

void tkKeyDownFunc(GLenum (*Func)(int, GLenum))
{

    KeyDownFunc = Func;
}

/******************************************************************************/

void tkMouseDownFunc(GLenum (*Func)(int, int, GLenum))
{

    MouseDownFunc = Func;
}

/******************************************************************************/

void tkMouseUpFunc(GLenum (*Func)(int, int, GLenum))
{

    MouseUpFunc = Func;
}

/******************************************************************************/

void tkMouseMoveFunc(GLenum (*Func)(int, int, GLenum))
{

    MouseMoveFunc = Func;
}

/******************************************************************************/

void tkIdleFunc(void (*Func)(void))
{

    IdleFunc = Func;
}

/******************************************************************************/
