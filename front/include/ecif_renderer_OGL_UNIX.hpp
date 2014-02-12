/*****************************************************************************
 *
 *  Elmer, A Finite Element Software for Multiphysical Problems
 *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
 * 
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program (in file fem/GPL-2); if not, write to the 
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 *  Boston, MA 02110-1301, USA.
 *
 *****************************************************************************/

/***********************************************************************
Program:    ELMER Front
Module:     ecif_oglUNIXrenderer.hpp
Language:   C++
Date:       30.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Include type file. Unix specific parts for OGLrenderer class.
  Uses glx-utility functions under Unix.

************************************************************************/

#include <X11/keysym.h>
#include <GL/glx.h>


void
Renderer_OGL::createGLWindow(Hinst app_instance, const char* pName,
                             int xPos, int yPos, int pWidth, int pHeight,
                             WindowInfo& winfo)
{
  int dblBuf[] = {GLX_RGBA, GLX_DEPTH_SIZE, 16, GLX_DOUBLEBUFFER, None};

  Display* dpy = XOpenDisplay(NULL);

  winfo.display = dpy;

  XVisualInfo* vi = glXChooseVisual(dpy, DefaultScreen(dpy), dblBuf);
  GLXContext cx = glXCreateContext(dpy, vi, None, GL_TRUE);
  Colormap cmap = XCreateColormap(dpy,
                                  RootWindow(dpy, vi->screen),
                                  vi->visual, AllocNone);

  XSetWindowAttributes swa;
  swa.colormap = cmap;
  swa.border_pixel = 0;
  swa.event_mask = ExposureMask | StructureNotifyMask;
  swa.event_mask |= ButtonPressMask | ButtonReleaseMask | OwnerGrabButtonMask;
  swa.event_mask |= KeyPressMask | KeyReleaseMask;
  swa.event_mask |= PointerMotionHintMask | PointerMotionMask;

  Window hWnd = XCreateWindow(dpy,
                              RootWindow(dpy, vi->screen),
                              xPos, yPos, pWidth, pHeight,
                              0, vi->depth,
                              InputOutput, vi->visual,
                              CWBorderPixel|CWColormap|CWEventMask,
                              &swa);

  winfo.window = hWnd;

  XSetStandardProperties(dpy, hWnd, oglWinClass, oglWinClass, None, 0, NULL, NULL);
  glXMakeCurrent(dpy, hWnd, cx);

#if 0
  // NOTE: "messages" field not supported in SGI XWindows
  // We are not interested in moves
  XWMHints wmHint;
  wmHint.messages = WindowMoved;
  XsetWMHints(dpy, hWnd, &wmHint);
#endif

}


void
Renderer_OGL::destroyWindow(Hdisplay* display, Hwindow window)
{
  // NOTE: X-destroy seems to cause problems in Unix, so it is
  // not currently in use!
  return;

  glXDestroyContext(display, glXGetCurrentContext());
  XDestroyWindow(display, window);
  XCloseDisplay(display);
  glXMakeCurrent(NULL, NULL, NULL);

  status = HAS_NO_WINDOW;
  visible = false;
}


//This method is called fex. from the Tcl_MainLoop function
//always when there are no tcl-events.
//Method checks if it is necessary really to repaint the renderer
//window.
void
Renderer_OGL::dummyWindowProc()
{
  if ( status == SHOW ) {
    paintRenderer();
  }
}


void
Renderer::findRendererWindowSize(int& w, int& h)
{
  int x, y;
  unsigned int w_w, w_h, border_w, depth;
  Display* dpy = rendererDisplay;

  XGetGeometry( dpy,
                glXGetCurrentDrawable(),
                &RootWindow(dpy, DefaultScreen(dpy)),
                &x, &y,
                &w_w, &w_h,
                &border_w, &depth);

  w = (int) w_w;
  h = (int) w_h;
}


// Method paints the renderer window.
void
Renderer_OGL::paintRenderer()
{
  Window root_w, child_w;
  int root_x, root_y, window_x, window_y;
  static int xPos1 = 0;
  static int yPos1 = 0;
  int xPos2, yPos2;
  unsigned int mk_state;

  // For double click checking
  static double time0 = 0.0;
  double time1;

  if (rendererDisplay == NULL) {
    return;
  }

  Renderer* renderer = theControlCenter->getRenderer();

  int x_size, y_size;
  renderer->getRendererWindowSize(x_size, y_size);

  if ( !visible ) {
    XMapRaised(rendererDisplay, rendererWindow);
    visible = true;
  }

  // *** X-event loop.
  Event event;
  int event_mask;
  event_mask = ExposureMask | StructureNotifyMask;
  event_mask |= ButtonPressMask | ButtonReleaseMask | OwnerGrabButtonMask;
  event_mask |= KeyPressMask | KeyReleaseMask;
  event_mask |= PointerMotionHintMask | PointerMotionMask;

  bool need_redisplay = false;
  bool is_double_click = false;

  // These are used for key-press handling
  static bool ctrlPressed = false;
  static char key_buffer[20];
  int key_buffer_len = 20;
  int key_buffer_count;
  KeySym key;
  XComposeStatus cs;

  if (0 != XCheckWindowEvent( rendererDisplay, rendererWindow,
                              event_mask, &event) ) {

    switch(event.type) {

    case ButtonPress:

       XQueryPointer(rendererDisplay, event.xmotion.window, &root_w, &child_w,
                    &root_x, &root_y, &window_x, &window_y,
                    &mk_state);

      xPos1 = window_x;
      yPos1 = y_size - window_y;

      // Check if we have a double click (left button!)
      if ( mk_state & Button1Mask ) {

      	time1 = doubleClickTimer->getLapTime(WALL_TIME);

      	if ( (time1 - time0) < 0.5 ) {
          is_double_click = true;
        }

        time0 = time1;
      }

      // Single click with shif or control
      if ( mk_state & (MK_SHIFT | MK_CONTROL) ) {
        mouse_clck_action((int)mk_state, xPos1, yPos1);

      // Double click with left button
      } else if (is_double_click) {
        mouse_dblclck_action((int)mk_state, xPos1, yPos1);
      }
      break;

    case ButtonRelease:
      xPos1 = -1;
      yPos1 = -1;
      break;

    case DestroyNotify:
      destroyWindow(rendererDisplay, rendererWindow);
      return;
      break;

    case KeyPress:

      key_buffer_count = XLookupString(&(event.xkey), key_buffer, key_buffer_len, &key, &cs);
      key_buffer[key_buffer_count] = '\0';

      // Ctrl pressed
      if ( event.xkey.state & MK_CONTROL ) {
        ctrlPressed = true;
      }

      // A letter pressed while Ctrl-key is down
      if (ctrlPressed) {

        switch (key) {

        case XK_b:
        case XK_B: theControlCenter->handleKeyAction(KEY_CTRL_B); break;
        case XK_h:
        case XK_H: theControlCenter->handleKeyAction(KEY_CTRL_H); break;
        case XK_l:
        case XK_L: theControlCenter->handleKeyAction(KEY_CTRL_L); break;
        case XK_r:
        case XK_R: theControlCenter->handleKeyAction(KEY_CTRL_R); break;
        case XK_x:
        case XK_X: theControlCenter->handleKeyAction(KEY_CTRL_X); break;
        case XK_y:
        case XK_Y: theControlCenter->handleKeyAction(KEY_CTRL_Y); break;
        case XK_z:
        case XK_Z: theControlCenter->handleKeyAction(KEY_CTRL_Z); break;
        }
      }
      break;

    case KeyRelease:

      // Ctrl pressed
      if ( !(event.xkey.state & MK_CONTROL) ) {
        ctrlPressed = false;
      }

      break;

    case ConfigureNotify:
      renderer->reshape();
      need_redisplay = true;
      break;

    case Expose:
      need_redisplay = true;
      break;

    case MotionNotify:

      XQueryPointer(rendererDisplay, event.xmotion.window, &root_w, &child_w,
                    &root_x, &root_y, &window_x, &window_y,
                    &mk_state);

      if (  mk_state != 0 && xPos1 != -1 ) {
        xPos2 = window_x;
        yPos2 = y_size - window_y;
        mouse_move_action((int)mk_state, xPos1, xPos2, yPos1, yPos2);
        xPos1 = xPos2;
        yPos1 = yPos2;
      } else if ( mk_state !=0 ) {
        xPos1 = window_x;
        yPos1 = y_size - window_y;
      }

      break;
    } // switch (event)

  } // if XCheckWindowEvent

  // Redraw if necessary.
  if (need_redisplay) {
    display(renderer);
  }
}


void
Renderer_OGL::setWindowTitle(char* title)
{
  XStoreName(rendererDisplay, rendererWindow, title);
}


void
Renderer_OGL::swapBuffers()
{
  glXSwapBuffers(rendererDisplay, glXGetCurrentDrawable());
}

// *** UNIX WindowProc for the OGLRenderer ***
LresCallback
Renderer_OGL::windowProcedure(Hwindow hWnd, Event wMsg,
                    Wparam wParam, Lparam lParam)
{
  return NULL;
}

