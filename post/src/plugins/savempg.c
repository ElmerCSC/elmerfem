//=============================================================================
// savempg.c
//
// Module for compressing and saving ElmerPost-pictures in MPEG formats.
//
// Compile e.g. as follows:
//
// Linux:
//
// gcc -shared -O -I/usr/include/ffmpeg -I/usr/include/tcl8.4
//    savempg.c -o savempg.o -lavcodec -lavutil -lswscale
//
// MinGW: 
// gcc -shared -O -I$FFMPEG/include -L$FFMPEG/lib -o savempg.dll 
//     savempg.c -lopengl32 -ltcl84 -lavcodec -lavutil -lswscale -lz
//
// ( Note that the libraries required depend on the libavcodec build. )
//
// Copy the shared library into $ELMER_POST_HOME/modules and run ElmerPost
//
// Usage in Elmer-Post:
//
//    savempg codec name          [ default: mpg1 ( mpg2, mpg2, yuv4 )   ]
//    savempg bitrate value       [ default: 1000000 bps                 ]
//    savempg gopsize value       [ default: 20 frames                   ]
//    savempg bframes value       [ default: 2 frames                    ]
//    savempg start file.es       [ initializes video compressor         ]
//    savempg append              [ adds current frame to video sequence ]
//    savempg stop                [ finalizes video compressor           ]
//    
// Parsed together from screensave.c && the ffmpeg samples
//
// Written by: ML, 30. Sept. 2007
//=============================================================================
#ifdef HAVE_AV_CONFIG_H
#undef HAVE_AV_CONFIG_H
#endif

#if defined(WIN32) || defined(win32)
#include <windows.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <GL/gl.h>
#include <tcl.h>
#include "avcodec.h"
#include "swscale.h"

#define INBUF_SIZE 4096
#define STATE_READY 0
#define STATE_STARTED 1
#define DEFAULT_B_FRAMES 2
#define DEFAULT_GOP_SIZE 20
#define DEFAULT_BITRATE 1000000
#define DEFAULT_MPG_BUFSIZE 500000

#undef MPEG4

#if !defined(INFINITY) 
#define INFINITY 999999999
#endif

typedef struct buffer_s {
  uint8_t *MPG;
  uint8_t *YUV;
  uint8_t *RGB;
  uint8_t *ROW;
} buffer_t;

void free_buffers( buffer_t *buff ) {
  free( buff->MPG );
  free( buff->YUV );
  free( buff->RGB );
  free( buff->ROW );
}

void SetMessage( Tcl_Interp *interp, char *message ) {
  sprintf( interp->result, "savempg: %s\n", message );
}

static double psnr( double d ) {
  if( d == 0 )
    return INFINITY;
  return -10.0 * log( d ) / log( 10.0 );
}

void print_info( int count_frames, AVCodecContext *context, int bytes ) {
  double tmp = context->width * context->height * 255.0 * 255.0;
  double Ypsnr = psnr( context->coded_frame->error[0] / tmp ); 
  double quality = context->coded_frame->quality/(double)FF_QP2LAMBDA;
  char pict_type = av_get_pict_type_char(context->coded_frame->pict_type);
  fprintf( stdout, "savempg: frame %4d: type=%c, ", count_frames, pict_type );
  fprintf( stdout, "size=%6d bytes, PSNR(Y)=%5.2f dB, ", bytes, Ypsnr );
  fprintf( stdout, "q=%2.1f\n", (float)quality );
  fflush( stdout );
}

static int SaveMPG( ClientData cl,Tcl_Interp *interp,int argc,char **argv ) {
  static AVCodec *codec = NULL;
  static AVCodecContext *context = NULL;
  static AVFrame *YUVpicture = NULL;
  static AVFrame *RGBpicture = NULL;
  static int bytes, PIXsize, stride;
  static int y, nx, ny, ox, oy, viewp[4];
  static int i_state = STATE_READY;
  static int initialized = 0;
  static int count_frames = 0;
  static int bitrate = DEFAULT_BITRATE;
  static int gopsize = DEFAULT_GOP_SIZE;
  static int bframes = DEFAULT_B_FRAMES;
  static int MPGbufsize = DEFAULT_MPG_BUFSIZE;
  static int codec_id = CODEC_ID_MPEG1VIDEO;
  static FILE *MPGfile;
  static char *state, fname[256];
  static buffer_t buff;
  static struct SwsContext *img_convert_ctx;

  if( argc < 2 ) {
    SetMessage( interp, "too few arguments" );
    return TCL_ERROR;
  }

  state = argv[1];
  //===========================================================================
  //                             SELECT CODEC
  //===========================================================================
  if( !strncmp( state, "codec", 5 ) ) {

    // Can't change while running:
    //----------------------------
    if( i_state != STATE_READY ) {
      SetMessage( interp, "can't change codec when running" );
      return TCL_ERROR;
    }

    // Default is MPEG1:
    //------------------
    codec_id = CODEC_ID_MPEG1VIDEO;

    if( argc >= 3 ) {
      char *c = argv[2];

      if( !strncmp( c, "mpg1", 4 ) ) {
	codec_id = CODEC_ID_MPEG1VIDEO;
	fprintf( stdout, "savempg: codec: mpg1\n" );
      }

      if( !strncmp( c, "mpg2", 4 ) ) {
	codec_id = CODEC_ID_MPEG2VIDEO;
	fprintf( stdout, "savempg: codec: mpg2\n" );
      }

      if( !strncmp( c, "mpg4", 4 ) ) {
	codec_id = CODEC_ID_MPEG4;
	fprintf( stdout, "savempg: codec: mpg4\n" );
      }

      if( !strncmp( c, "yuv4", 4 ) ) {
	codec_id = CODEC_ID_RAWVIDEO;
	fprintf( stdout, "savempg: codec: yuv4\n" );
      }
    }
    
    fflush( stdout );

    return TCL_OK;
  }

  //===========================================================================
  //                              SET BITRATE
  //===========================================================================
  if( !strncmp( state, "bitrate", 7 ) ) {

    // Can't change while running:
    //----------------------------
    if( i_state != STATE_READY ) {
      SetMessage( interp, "can't change bitrate when running" );
      return TCL_ERROR;
    }

    bitrate = DEFAULT_BITRATE;

    if( argc >= 3 )
      bitrate = atoi( argv[2] );
    
    if( bitrate < 100000 )
      bitrate = 100000;

    fprintf( stdout, "savempg: bitrate: %d bps\n", bitrate );
    fflush( stdout );

    return TCL_OK;
  }

  //===========================================================================
  //                              SET GOPSIZE
  //===========================================================================
  if( !strncmp( state, "gopsize", 7 ) ) {

    // Can't change while running:
    //----------------------------
    if( i_state != STATE_READY ) {
      SetMessage( interp, "can't change gopsize when running" );
      return TCL_ERROR;
    }

    gopsize = DEFAULT_GOP_SIZE;

    if( argc >= 3 )
      gopsize = atoi( argv[2] );

    if( gopsize < 1 )
      gopsize = 1;
    
    fprintf( stdout, "savempg: gop: %d frames\n", gopsize );
    fflush( stdout );

    return TCL_OK;
  }

  //===========================================================================
  //                              SET B-FRAMES
  //===========================================================================
  if( !strncmp( state, "bframes", 7 ) ) {

    // Can't change while running:
    //----------------------------
    if( i_state != STATE_READY ) {
      SetMessage( interp, "can't change bframes when running" );
      return TCL_ERROR;
    }

    bframes = DEFAULT_B_FRAMES;

    if( argc >= 3 )
      bframes = atoi( argv[2] );
    
    if( bframes < 0 )
      bframes = 0;

    fprintf( stdout, "savempg: bframes: %d\n", bframes );
    fflush( stdout );

    return TCL_OK;
  }

  //===========================================================================
  //                           START COMPRESSION
  //===========================================================================
  else if( !strncmp( state, "start", 5 ) ) {
    
    // Can't start if already running:
    //----------------------------------
    if( i_state != STATE_READY ) {
      SetMessage( interp, "already started" );
      return TCL_ERROR;
    }

    // Detemine output file name:
    //---------------------------
    if( argc < 3 ) {
      strcpy( fname, "elmerpost.es" );
    } else {
      strncpy( fname, argv[2], 256 );
    }
    
    // Open output file:
    //------------------
    if( (MPGfile = fopen(fname, "wb")) == NULL ) {
      SetMessage( interp, "can't open file" );
      return TCL_ERROR;
    }

    fprintf( stdout, "savempg: file: %s\n", fname );
    fprintf( stdout, "savempg: libavcodec: %s\n", 
	     AV_STRINGIFY(LIBAVCODEC_VERSION) );
    fprintf( stdout, "savempg: libavutil: %s\n",
	     AV_STRINGIFY(LIBAVUTIL_VERSION) );
    fprintf( stdout, "savempg: libswscale: %s\n",
	     AV_STRINGIFY(LIBSWSCALE_VERSION) );
    fflush( stdout );

    count_frames = 0;

    // Determine view port size etc.:
    //-------------------------------
    glGetIntegerv( GL_VIEWPORT, viewp );
    ox = viewp[0];
    oy = viewp[1];
    nx = viewp[2] + 1;
    ny = viewp[3] + 1;
    PIXsize = nx * ny;
    stride = 3 * nx;

    // Must be even:
    //--------------
    if( nx % 2 || ny % 2 ) {
      fprintf( stdout, "savempg: state: stopped\n" );
      fflush( stdout );
      fclose( MPGfile );
      SetMessage( interp, "view port must be even" );
      return TCL_ERROR;
    }

    // Allocate memory:
    //-----------------
    if( codec_id == CODEC_ID_RAWVIDEO )
      MPGbufsize = 3 * (PIXsize / 2);

    if ( !(buff.RGB = (uint8_t*)malloc(stride*ny)) ||
	 !(buff.ROW = (uint8_t*)malloc(stride)) ||
	 !(buff.YUV = (uint8_t*)malloc(3 * (PIXsize / 2))) ||
	 !(buff.MPG = (uint8_t*)malloc(MPGbufsize)) ) {
      fclose( MPGfile );
      SetMessage( interp, "can't allocate memory" );
      return TCL_ERROR;
    }

    // Initialize libavcodec:
    //-----------------------
    if( !initialized ) {
      avcodec_init();
      avcodec_register_all();
      initialized = 1;
    }

    // Choose codec:
    //--------------
    codec = avcodec_find_encoder( codec_id );

    if( !codec ) {
      free_buffers( &buff );
      fclose( MPGfile );
      SetMessage( interp, "can't find codec" );
      return TCL_ERROR;
    }
    
    // Init codec context etc.:
    //--------------------------
    context = avcodec_alloc_context();

    context->bit_rate = bitrate;
    context->width = nx;
    context->height = ny;
    context->time_base = (AVRational){ 1, 25 };
    context->gop_size = gopsize;
    context->max_b_frames = bframes;
    context->pix_fmt = PIX_FMT_YUV420P;
    context->flags |= CODEC_FLAG_PSNR;
    
    if( avcodec_open( context, codec ) < 0 ) {
      avcodec_close( context );
      av_free( context );
      free_buffers( &buff );
      fclose( MPGfile );
      SetMessage( interp, "can't open codec" );
      return TCL_ERROR;
    }

    YUVpicture = avcodec_alloc_frame();

    YUVpicture->data[0] = buff.YUV;
    YUVpicture->data[1] = buff.YUV + PIXsize;
    YUVpicture->data[2] = buff.YUV + PIXsize + PIXsize / 4;
    YUVpicture->linesize[0] = nx;
    YUVpicture->linesize[1] = nx / 2;
    YUVpicture->linesize[2] = nx / 2;

    RGBpicture = avcodec_alloc_frame();

    RGBpicture->data[0] = buff.RGB;
    RGBpicture->data[1] = buff.RGB;
    RGBpicture->data[2] = buff.RGB;
    RGBpicture->linesize[0] = stride;
    RGBpicture->linesize[1] = stride;
    RGBpicture->linesize[2] = stride;
    
    // Write .ym4 header for raw video:
    //----------------------------------
    if( codec_id == CODEC_ID_RAWVIDEO ) {
      fprintf( MPGfile, "YUV4MPEG2 W%d H%d F25:1 Ip A1:1", nx, ny );
      fprintf( MPGfile, "%c", (char)0x0a );
    }

    // Set state "started":
    //----------------------
    i_state = STATE_STARTED;
    fprintf( stdout, "savempg: state: started\n" );
    fflush( stdout );
    
    return TCL_OK;
  } 

  //===========================================================================
  //                             COMPRESS FRAME
  //===========================================================================
  else if( !strncmp( state, "append", 6) ) {

    // Can't compress if status != started:
    //-------------------------------------
    if( i_state != STATE_STARTED ) {
      SetMessage( interp, "not started" );
      return TCL_ERROR;
    }

    // Read RGB data:
    //---------------
    glReadBuffer( GL_FRONT );
    glReadPixels( ox, oy, nx, ny, GL_RGB, GL_UNSIGNED_BYTE, buff.RGB );

    // The picture is upside down - flip it:
    //---------------------------------------
    for( y = 0; y < ny/2; y++ ) {
      uint8_t *r1 = buff.RGB + stride * y;
      uint8_t *r2 = buff.RGB + stride * (ny - 1 - y);
      memcpy( buff.ROW, r1, stride );
      memcpy( r1, r2, stride );
      memcpy( r2, buff.ROW, stride );
    }

    // Convert to YUV:
    //----------------
    if( img_convert_ctx == NULL )
      img_convert_ctx = sws_getContext( nx, ny, PIX_FMT_RGB24,
					nx, ny, PIX_FMT_YUV420P,
					SWS_BICUBIC, NULL, NULL, NULL );
    
    if( img_convert_ctx == NULL ) {
      SetMessage( interp, "can't initialize scaler context" );
      return TCL_ERROR;
    }
    
    sws_scale( img_convert_ctx, RGBpicture->data, RGBpicture->linesize,
	       0, ny, YUVpicture->data, YUVpicture->linesize );

    // Encode frame:
    //--------------
    bytes = avcodec_encode_video( context, buff.MPG, 
				  MPGbufsize, YUVpicture );

    count_frames++;
    print_info( count_frames, context, bytes );
    fflush( stdout );

    if( codec_id == CODEC_ID_RAWVIDEO ) {
      fprintf( MPGfile, "FRAME" );
      fprintf( MPGfile, "%c", (char)0x0a );
    }

    fwrite( buff.MPG, 1, bytes, MPGfile );
    
    return TCL_OK;
  } 

  //===========================================================================
  //                           STOP COMPRESSION
  //===========================================================================
  else if( !strncmp( state, "stop", 4) ) {

    // Can't stop if status != started:
    //---------------------------------
    if( i_state != STATE_STARTED ) {
      SetMessage( interp, "not started" );
      return TCL_ERROR;
    }

    // Get the delayed frames, if any:
    //--------------------------------
    for( ; bytes; ) {
      bytes = avcodec_encode_video( context, buff.MPG, MPGbufsize, NULL );

      count_frames++;
      print_info( count_frames, context, bytes );
      fflush( stdout );

      fwrite( buff.MPG, 1, bytes, MPGfile );
    }
    
    // Add sequence end code:
    //-----------------------
    if( codec_id == CODEC_ID_MPEG1VIDEO ) {
      buff.MPG[0] = 0x00;
      buff.MPG[1] = 0x00;
      buff.MPG[2] = 0x01;
      buff.MPG[3] = 0xb7;
      fwrite( buff.MPG, 1, 4, MPGfile );
    }

    // Finalize:
    //-----------
    avcodec_close( context );
    av_free( context );
    av_free( YUVpicture );
    av_free( RGBpicture );
    free_buffers( &buff );
    fclose( MPGfile );
    
    i_state = STATE_READY;
    fprintf( stdout, "savempg: state: stopped\n" );
    fflush( stdout );

    return TCL_OK;
  } 

  //===========================================================================
  //                              UNKNOWN STATE
  //===========================================================================
  else {
    SetMessage( interp, "unknown request" );
    return TCL_ERROR;
  }

  // We should never ever arrive at this point:
  //-------------------------------------------
  SetMessage( interp, "unknown error" );
  return TCL_ERROR;
}

#if defined(WIN32) || defined(win32)
__declspec(dllexport)
#endif
     
int Savempg_Init( Tcl_Interp *interp ) {
  Tcl_CreateCommand( interp, "savempg", (Tcl_CmdProc *)SaveMPG,
		     (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );
  return TCL_OK;
}
