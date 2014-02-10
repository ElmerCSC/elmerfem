/*
 * Routine to manipulate size and position of the graphics window
 *
 * Compile commands
 *
 * DEC:    cc -shared -o window window.c -I$ELMER_HOME/include
 * LINUX: gcc -shared -o window window.c -I$ELMER_HOME/include
 * SGI:    cc -shared -o window window.c -I$ELMER_HOME/include -n32
 * WINDOWS:
 *  cl -c -I%ELMER_HOME%\include -DWIN32 window.c
 *  link /dll /libpath:%ELMER_HOME%\lib window.obj tcl81.lib opengl32.lib
 *
 * in elmerpost:
 *
 * load window window
 *
 * and then you have the commands:
 *
 * winpos xpos ypos
 * winsize width height
 */

#if defined(WIN32) || defined(win32)
#include <windows.h>
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
  
#include <GL/gl.h>
#include <GL/glu.h>

#include <tcl.h>

#if !(defined(WIN32) || defined(win32))
#include <string.h>
extern int errno;
#include <X11/Xlib.h>
#endif

#ifndef WIN32
Display *auxXDisplay();
#endif

static Display *tkXDisplay()
{
    Display *ptr = NULL;
#ifndef WIN32
    ptr = auxXDisplay();
#endif
    return ptr;
}

static Window tkXWindow()
{
    Window ptr = 0;

#ifdef WIN32
    ptr = auxGetHWND();
#else
    ptr = auxXWindow();
#endif

    return ptr;
}

static int WindowSize( ClientData cl,Tcl_Interp *interp,int argc,char **argv )
{
   int  width, height;

   if ( argc < 3 ) {
      sprintf( interp->result, "Usage: winsize width height" );
      return TCL_ERROR;
   }

   width  = atoi( *++argv );
   height = atoi( *++argv );

   XResizeWindow( (Display *)tkXDisplay(), tkXWindow(), width, height );

   return TCL_OK;
}

static int WindowPosition( ClientData cl,Tcl_Interp *interp,int argc,char **argv )
{
   int  ox,oy;

   if ( argc < 3 ) {
      sprintf( interp->result, "Usage: winpos xpos ypos" );
      return TCL_ERROR;
   }

   ox = atoi( *++argv );
   oy = atoi( *++argv );

   XMoveWindow( (Display *)tkXDisplay(), tkXWindow(), ox, oy );

   return TCL_OK;
}

#if defined(WIN32) || defined(win32)
__declspec(dllexport)
#endif
int Window_Init( Tcl_Interp *interp )
{
   Tcl_CreateCommand( interp, "winsize", WindowSize,
      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

   Tcl_CreateCommand( interp, "winpos", WindowPosition, 
      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

   return TCL_OK;
}
