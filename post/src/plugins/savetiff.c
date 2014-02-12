// savetiff.c
//
// Module for compressing and saving ElmerPost-pictures in TIFF-format
//
// Compile e.g. as follows:
//
//   MinGW: gcc -Wall -shared -O -o savetiff.dll savetiff.c -lopengl32 -ltcl84 -ltiff
//
// Copy the shared library into $ELMER_POST_HOME/modules and ruin ElmerPost
//
// Usage:
//
//   Elmer-Post: savetiff file
//
// Defaults:
//
//   file = elmerpost.tif
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
#include <tiffio.h>

#if !(defined(WIN32) || defined(win32))
#include <errno.h>
extern int errno;
#endif

static int SaveTIFF( ClientData cl,Tcl_Interp *interp,int argc,char **argv ) {
  unsigned char *buffer;
  char fname[256];
  int nx, ny, ox, oy, viewp[4];
  TIFF *image;
  int y, stride;

  // Determine file name:
  //---------------------
  if( argc < 2 ) {
    strcpy( fname, "elmerpost.tif" );
  } else {
    strncpy( fname,argv[1], 256 );
  }

  image = TIFFOpen( fname, "w" );  
  if( image==NULL ) {
#if defined(WIN32) || defined(win32)
    sprintf( interp->result, "savetiff: can't open [%s] for writing!\n",fname );
#else
    sprintf( interp->result, "savetiff: can't open [%s] for writing:\n%s\n", 
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
    TIFFClose( image );
#if defined(WIN32) || defined(win32)
    sprintf( interp->result, "savetiff: can't allocate memory!\n" );
#else
    sprintf( interp->result, "savetiff: can't allocate memory:\n%s\n",
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

  TIFFSetField( image, TIFFTAG_IMAGEWIDTH, nx );
  TIFFSetField( image, TIFFTAG_IMAGELENGTH, ny );
  TIFFSetField( image, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE );
  TIFFSetField( image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG );
  TIFFSetField( image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB );
  TIFFSetField( image, TIFFTAG_BITSPERSAMPLE, 8 );
  TIFFSetField( image, TIFFTAG_SAMPLESPERPIXEL, 3 );

  if( TIFFWriteEncodedStrip( image, 0, buffer, nx*ny*3) == 0 ) {
    TIFFClose( image );    
#if defined(WIN32) || defined(win32)
    sprintf( interp->result, "savetiff: unable to encode picture\n" );
#else
    sprintf( interp->result, "savetiff: unable to encode picture\n%s\n",
	     strerror(errno) );
#endif
    free( buffer );
    return TCL_ERROR;
  }

  TIFFClose( image );
  free( buffer );

  fprintf( stdout, "done\n");
  fflush( stdout );

  return TCL_OK;
}

#if defined(WIN32) || defined(win32)
__declspec(dllexport)
#endif
     
int Savetiff_Init( Tcl_Interp *interp ) {
  Tcl_CreateCommand( interp, "savetiff", (Tcl_CmdProc *)SaveTIFF,
		     (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );
  return TCL_OK;
}
