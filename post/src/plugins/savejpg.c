// savejpg.c
//
// Module for compressing and saving ElmerPost-pictures in jpg format
//
// Compile e.g. as follows:
//
//   Linux: gcc -I/usr/include/tcl8.4 -L/usr/lib/tcl8.4 -Wall -shared -O -o savejpg.so savejpg.c -lGL -ltcl8.4 -ljpeg
//   MinGW: gcc -Wall -shared -O -o savejpg.dll savejpg.c -lopengl32 -ltcl84 -ljpeg
//
// Copy the shared library into $ELMER_POST_HOME/modules and ruin ElmerPost
//
// Usage:
//
//   Elmer-Post: savejpg file quality
//
// Defaults:
//
//   file = elmerpost.jpg
//   quality = 85 (1=bad, 100=good)
//
// Modified from screensave.c 
//
// Written by: ML, 29. Sept. 2007

#if defined(WIN32) || defined(win32)
#include <windows.h>
// This is to avoid redefinition of INT32 in jmoreconfig.h:
#define XMD_H
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <GL/gl.h>
#include <tcl.h>

#include <jpeglib.h>

#if !(defined(WIN32) || defined(win32))
#include <errno.h>
extern int errno;
#endif

static int SaveJPG( ClientData cl,Tcl_Interp *interp,int argc,char **argv ) {
  unsigned char *buffer;
  char fname[256];
  int nx, ny, ox, oy, viewp[4];
  struct jpeg_compress_struct cinfo;  
  struct jpeg_error_mgr jerr;
  FILE *outfile;
  JSAMPROW row_pointer[1];
  int row_stride, quality = 85;
  
  // Determine file name & quality:
  //--------------------------------
  if( argc < 2 ) {
    strcpy( fname, "elmerpost.jpg" );
  } else {
    strncpy( fname,argv[1], 256 );
    if( argc == 3 ) {
      quality = atoi( argv[2] );
    }
  }
  
  // Open output file:
  //-------------------
  if( (outfile = fopen(fname, "wb")) == NULL ) {
#if defined(WIN32) || defined(win32)
    sprintf( interp->result, "savejpg: can't open [%s] for writing!\n",fname );
#else
    sprintf( interp->result, "savejpg: can't open [%s] for writing:\n%s\n", 
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
  if ( !(buffer=(unsigned char *)malloc(nx*ny*3)) ) {
#if defined(WIN32) || defined(win32)
    sprintf( interp->result, "savejpg: can't allocate enough memory!\n" );
#else
    sprintf( interp->result, "savejpg: can't allocate enough memory:\n%s\n",
	     strerror(errno) );
#endif
    fclose( outfile );
    return TCL_ERROR;
  }
  fprintf( stdout, "Saving %s (quality: %d)... ", fname, quality );
  fflush( stdout );
  
  // Copy RGB-data into buffer:
  //----------------------------
  glReadBuffer( GL_FRONT );
  glReadPixels( ox, oy, nx, ny, GL_RGB, GL_UNSIGNED_BYTE, buffer );
  
  // Compress into JPG(YUV) & save:  
  //-------------------------------
  cinfo.err = jpeg_std_error( &jerr );
  jpeg_create_compress( &cinfo );
  jpeg_stdio_dest( &cinfo, outfile );
  
  cinfo.image_width = nx;
  cinfo.image_height = ny;
  cinfo.input_components = 3;
  cinfo.in_color_space = JCS_RGB;
  
  jpeg_set_defaults( &cinfo );
  jpeg_set_quality( &cinfo, quality, TRUE );
  jpeg_start_compress( &cinfo, TRUE );
  
  row_stride = 3*nx;
  while( cinfo.next_scanline < cinfo.image_height ) {
    row_pointer[0] = &buffer[ (ny-1-cinfo.next_scanline) * row_stride ];
    (void)jpeg_write_scanlines( &cinfo, row_pointer, 1 );
  }
  
  jpeg_finish_compress( &cinfo );
  fclose( outfile );
  
  jpeg_destroy_compress( &cinfo );
  fprintf( stdout, "done [ok]\n" );  
  fflush( stdout );

  free( buffer );
  return TCL_OK;
}

#if defined(WIN32) || defined(win32)
__declspec(dllexport)
#endif
     
int Savejpg_Init( Tcl_Interp *interp ) {
  Tcl_CreateCommand( interp, "savejpg", (Tcl_CmdProc *)SaveJPG,
		     (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );
  return TCL_OK;
}
