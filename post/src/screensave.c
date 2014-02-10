/*
 * Routine to save elmerpost graphics window to disk in ppm(P6) format.
 *
 * Compile commands
 *
 * DEC: cc -shared -o grab grab.c -I$ELMER_HOME/include
 * SGI: cc -shared -o grab grab.c -I$ELMER_HOME/include -n32
 * WINDOWS:
 *  cl -c -I%ELMER_HOME%\include -DWIN32 grab.c
 *  link /dll /libpath:%ELMER_HOME%\lib grab.obj tcl81.lib opengl32.lib
 *
 * in elmerpost:
 *
 * load screensave screensave
 *
 * and then you have the command:
 *
 * screensave filename.ppm
 */

#if defined(WIN32) || defined(win32)
#include <windows.h>
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>

  
#include <GL/gl.h>
#include <GL/glu.h>


#include "tcl.h"

#if !(defined(WIN32) || defined(win32))
#include <errno.h>
extern int errno;
#endif

static int ScreenSave( ClientData cl,Tcl_Interp *interp,int argc,char **argv )
{
   unsigned char *buffer;
   FILE *fp;
   char fname[256];
   int  i,nx,ny,ox,oy,viewp[4];

   if ( argc < 2 ) {
      strcpy( fname, "elmerpost.ppm" );
   } else {
      strncpy( fname,argv[1], 256 );
   }
   if ( !(fp = fopen( fname, "wb" ) ) ) {
#if defined(WIN32) || defined(win32)
      sprintf( interp->result, 
       "screensave: can't open [%s] for writing!\n", argv[1] );
#else
      sprintf( interp->result, 
        "screensave: can't open [%s] for writing:\n%s\n", 
              argv[1], strerror(errno) );
#endif
      return TCL_ERROR;
   }

   glGetIntegerv( GL_VIEWPORT, viewp );
   ox = viewp[0];
   oy = viewp[1];
   nx = viewp[2]+1;
   ny = viewp[3]+1;

   if ( !(buffer = (unsigned char *)malloc( nx*ny*3 ) ) )
   {
      fclose( fp );
#if defined(WIN32) || defined(win32)
      sprintf( interp->result,
            "screensave: can't allocate enough memory!\n" );
#else
      sprintf( interp->result,
        "screensave: can't allocate enough memory:\n%s\n",
            strerror(errno) );
#endif
      return TCL_ERROR;
   }

   fprintf( stderr, "screensave: reading pixels..." );

   glReadBuffer( GL_FRONT );
   glReadPixels( ox,oy,nx,ny, GL_RGB, GL_UNSIGNED_BYTE, buffer );

   fprintf( stderr, "writing file [%s]...", fname  );

   fprintf( fp, "P6\n%d %d\n255\n", nx, ny );

   for( i=ny-1; i>=0; i-- )
      if ( fwrite( &buffer[i*nx*3], 1, 3*nx, fp  ) != 3*nx ) {
         fclose( fp );
         free( buffer );
#if defined(WIN32) || defined(win32)
         sprintf( interp->result,
           "screensave: error writing to [%s]!\n", argv[1] );
#else
         sprintf( interp->result,
           "screensave: error writing to [%s]:\n%s\n",
               argv[1], strerror(errno) );
#endif
         return TCL_ERROR;
      }

   fclose( fp );
   free( buffer );

   fprintf( stderr, "done...[ok]\n" );

   return TCL_OK;
}

#if defined(WIN32) || defined(win32)
__declspec(dllexport)
#endif
int Screensave_Init( Tcl_Interp *interp )
{
   Tcl_CreateCommand( interp, "screensave", ScreenSave,
      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );

   return TCL_OK;
}
