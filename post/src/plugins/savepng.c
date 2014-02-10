// savepng.c
//
// Module for compressing and saving ElmerPost-pictures in PNG-format
//
// Compile e.g. as follows:
//
//   MinGW: gcc -Wall -shared -O -o savepng.dll savepng.c -lopengl32 -ltcl84 -lpng
//
// Copy the shared library into $ELMER_POST_HOME/modules and ruin ElmerPost
//
// Usage:
//
//   Elmer-Post: savepng file
//
// Defaults:
//
//   file = elmerpost.png
//
// Modified from screensave.c 
//
// Written by: ML, 08. Jan. 2008

#if defined(WIN32) || defined(win32)
#include <windows.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <GL/gl.h>
#include <tcl.h>

#define PNG_DEBUG 3
#include <png.h>

#if !(defined(WIN32) || defined(win32))
#include <errno.h>
extern int errno;
#endif

static int SavePNG( ClientData cl,Tcl_Interp *interp,int argc,char **argv ) {
  unsigned char *buffer;
  char fname[256];
  int nx, ny, ox, oy, viewp[4];
  int y, stride;

  png_byte color_type = PNG_COLOR_TYPE_RGB;
  png_byte bit_depth = 8;
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep *row_pointers;

  // Determine file name:
  //---------------------
  if( argc < 2 ) {
    strcpy( fname, "elmerpost.png" );
  } else {
    strncpy( fname, argv[1], 256 );
  }

  // Open file:
  //-----------
  FILE *image = fopen( fname, "wb" );
  if( image==NULL ) {
#if defined(WIN32) || defined(win32)
    sprintf( interp->result, "savepng: can't open [%s] for writing!\n",fname );
#else
    sprintf( interp->result, "savepng: can't open [%s] for writing:\n%s\n", 
	     fname, strerror(errno) );
#endif
    return TCL_ERROR;
  }

  // Determine picture size:
  //------------------------
  glGetIntegerv( GL_VIEWPORT, viewp );
  ox = viewp[0];
  oy = viewp[1];
  nx = viewp[2]+1;
  ny = viewp[3]+1;

  // Allocate buffer:
  //------------------
  buffer = (unsigned char *) malloc( nx*ny*3 );
  if ( buffer==NULL ) {
#if defined(WIN32) || defined(win32)
    sprintf( interp->result, "savepng: can't allocate memory!\n" );
#else
    sprintf( interp->result, "savepng: can't allocate memory:\n%s\n",
	     strerror(errno) );
#endif
    return TCL_ERROR;
  }

  fprintf( stdout, "Saving %s ... ", fname );
  fflush( stdout );
  
  // Copy RGB-data into buffer:
  //----------------------------
  glReadBuffer( GL_FRONT );
  glReadPixels( ox, oy, nx, ny, GL_RGB, GL_UNSIGNED_BYTE, buffer );
  
  // Flip the picture:
  //------------------
  stride = 3*nx;
  for( y=0; y<ny/2; y++ ) {
    unsigned char *r1 = buffer + stride*y;
    unsigned char *r2 = buffer + stride*(ny-1-y);
    memcpy( buffer, r1, stride );
    memcpy( r1, r2, stride );
    memcpy( r2, buffer, stride );
  }

  // Set up row pointers:
  //----------------------
  row_pointers = (png_bytep *) malloc( sizeof(png_bytep) * ny );
  for( y=0; y<ny; y++ ) 
    row_pointers[y] = buffer + stride*y;
    
  // Initialize PNG write struct:
  //------------------------------
  png_ptr = png_create_write_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );
  if( !png_ptr ) {
#if defined(WIN32) || defined(win32)
    sprintf( interp->result, "savepng: can't create write struct!\n" );
#else
    sprintf( interp->result, "savepng: can't create write struct:\n%s\n",
	     strerror(errno) );
#endif
    return TCL_ERROR;
  }

  // Create info struct:
  //---------------------
  info_ptr = png_create_info_struct( png_ptr );
  if( !info_ptr ) {
#if defined(WIN32) || defined(win32)
    sprintf( interp->result, "savepng: can't create info struct!\n" );
#else
    sprintf( interp->result, "savepng: can't create info struct:\n%s\n",
	     strerror(errno) );
#endif
    return TCL_ERROR;
  }
  
  // Init io:
  //----------
  if( setjmp(png_jmpbuf(png_ptr)) ) {
#if defined(WIN32) || defined(win32)
    sprintf( interp->result, "savepng: error in init_io!\n" );
#else
    sprintf( interp->result, "savepng: error in init_io:\n%s\n",
	     strerror(errno) );
#endif
    return TCL_ERROR;
  }
  
  png_init_io( png_ptr, image );
  
  // Write header:
  //---------------
  if( setjmp(png_jmpbuf(png_ptr)) ) {
#if defined(WIN32) || defined(win32)
    sprintf( interp->result, "savepng: can't write headers!\n" );
#else
    sprintf( interp->result, "savepng: can't write headers:\n%s\n",
	     strerror(errno) );
#endif
    return TCL_ERROR;
  }
  
  png_set_IHDR( png_ptr, info_ptr, nx, ny,
		bit_depth, color_type, PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE );
  
  png_write_info( png_ptr, info_ptr );
  
  // Write bytes:
  //--------------
  if( setjmp(png_jmpbuf(png_ptr)) ) {
#if defined(WIN32) || defined(win32)
    sprintf( interp->result, "savepng: can't write bytes!\n" );
#else
    sprintf( interp->result, "savepng: can't write bytes:\n%s\n",
	     strerror(errno) );
#endif
    return TCL_ERROR;
  }

  png_write_image( png_ptr, row_pointers );
  
  // End write:
  //-----------
  if( setjmp(png_jmpbuf(png_ptr)) ) {
#if defined(WIN32) || defined(win32)
    sprintf( interp->result, "savepng: can't finalize!\n" );
#else
    sprintf( interp->result, "savepng: can't finalize:\n%s\n",
	     strerror(errno) );
#endif
    return TCL_ERROR;
  }
  
  png_write_end( png_ptr, NULL );

  // Cleanup:
  //----------
  free(row_pointers);
  fclose( image );
  free( buffer );
  fprintf( stdout, "done\n");
  fflush( stdout );

  return TCL_OK;
}

#if defined(WIN32) || defined(win32)
__declspec(dllexport)
#endif
     
int Savepng_Init( Tcl_Interp *interp ) {
  Tcl_CreateCommand( interp, "savepng", (Tcl_CmdProc *)SavePNG,
		     (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );
  return TCL_OK;
}
