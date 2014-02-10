//-----------------------------------------------------------------------------------
// Filename:    fttext.cpp
// Description: Provides FTGL text rendering functionality for ElmerPost
// Usage:       fttext string [x y]
//              ftfont font [size r g b]
// Compilation: export CFLAGS="-DHAVE_FTGL -I/usr/include/freetype2 -I/usr/include/FTGL"
//              export CXXFLAGS="-DHAVE_FTGL -I/usr/include/freetype2 -I/usr/include/FTGL"
//              export LIBS="-lfreetype -lftgl"
// Written by:  Mikko Lyly
// Date:        15. Jan 2008
//-----------------------------------------------------------------------------------
#if defined(HAVE_CONFIG_H)
#include "../config.h"
#endif

#if defined(HAVE_FTGL_NEW) || defined(HAVE_FTGL_OLD)

#include <stdlib.h>
#include <stdio.h>
#include <GL/gl.h>
#include <tcl.h>

#if defined(WIN32) || defined(win32)
#include <windows.h>
#include "FTGL/FTGL.h"
#include "FTGL/FTGLPixmapFont.h"
#else
#if defined(HAVE_FTGL_NEW)
#include "FTGL/ftgl.h"
#else
#include "FTGL/FTGL.h"
#include "FTGL/FTGLPixmapFont.h"
#endif // HAVE_FTGL_NEW
#endif // WIN32 || win32

#define FTGLSTRLEN 4096

typedef struct {
  double x, y;
  int size;
  double r, g, b;
  int init_ok;
  FTGLPixmapFont *Font;
  char txt[FTGLSTRLEN];
  char ttf[FTGLSTRLEN];
  char current_ttf[FTGLSTRLEN];
  char ttffile[FTGLSTRLEN];
} ftgl_t;

static ftgl_t ftgl;

extern "C" void (*user_hook_before_all)();
extern "C" void (*user_hook_after_all)();

static void FtInit() {
  strcpy(ftgl.txt, "");
  ftgl.x = -0.9;
  ftgl.y = -0.9;
  strcpy(ftgl.ttf, "FreeSans");
  ftgl.size = 30;
  ftgl.r = 1.0;
  ftgl.g = 1.0;
  ftgl.b = 1.0;
  ftgl.init_ok = 1;
  return;
}

extern "C" void FtRender() {
  unsigned int r = strlen(ftgl.current_ttf);
  unsigned int s = strlen(ftgl.ttf);
  int t = strcmp(ftgl.ttf, ftgl.current_ttf);

  if( (r!=s) || (t!=0) || ftgl.Font->Error() ) {
    char *elmer_post_home = getenv("ELMER_POST_HOME");
    fprintf(stdout, "fttext: getenv: ELMER_POST_HOME=%s\n", 
	    elmer_post_home);
    fflush(stdout);
    
#if defined(WIN32) || defined(win32)
    sprintf(ftgl.ttffile, "%s\\fonts\\TrueType\\%s.ttf",
	    elmer_post_home, ftgl.ttf);
#else
    sprintf(ftgl.ttffile, "%s/fonts/TrueType/%s.ttf",
	    elmer_post_home, ftgl.ttf);
#endif

    fprintf(stdout, "fttext: load: %s\n", ftgl.ttffile);
    fflush(stdout);
    
    delete ftgl.Font;
    ftgl.Font = new FTGLPixmapFont(ftgl.ttffile);

    if(ftgl.Font->Error()) {
      fprintf(stderr, "fttext: error: load font failed!\n");
      fflush(stderr);
      return;
    }
    memset(ftgl.current_ttf, 0, FTGLSTRLEN);
    strncpy(ftgl.current_ttf, ftgl.ttf, strlen(ftgl.ttf));
  }
  
  if( ftgl.Font->Error() ) {
    fprintf(stderr, "fttext: error: no font loaded!\n");
    fflush(stderr);
    return;
  }

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  glColor3f(ftgl.r, ftgl.g, ftgl.b);
  glRasterPos3f(ftgl.x, ftgl.y, 0.0);
  
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_1D);

  ftgl.Font->FaceSize(ftgl.size);
  ftgl.Font->Render(ftgl.txt);

  glEnable(GL_TEXTURE_1D);
  glEnable(GL_LIGHTING);

  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  return;
}

extern "C" int FtFont(ClientData cl, Tcl_Interp *interp, 
		      int argc, char **argv) {
  if(!ftgl.init_ok)
    FtInit();
  
  if(argc > 1) 
    strcpy( ftgl.ttf, argv[1] );

  if(argc > 2)
    ftgl.size = atoi(argv[2]);

  if(argc > 3)
    ftgl.r = atof(argv[3]);

  if(argc > 4)
    ftgl.g = atof(argv[4]);

  if(argc > 5)
    ftgl.b = atof(argv[5]);
  
  Tcl_Eval(interp, "display");

  return TCL_OK;
}

extern "C" int FtText(ClientData cl, Tcl_Interp *interp, 
		      int argc, char **argv) {

  if(!ftgl.init_ok)
    FtInit();
  
  strcpy(ftgl.txt, "");
  ftgl.x = -0.9;
  ftgl.y = -0.9;
  
  if(argc > 1) {

    // TCL uses internally UTF-8. We want
    //  text in system default for FTGL:
    //---------------------------------------
    Tcl_DString ds;
    Tcl_DStringInit(&ds);
    char *res = Tcl_UtfToExternalDString(NULL, argv[1], -1, &ds);
    strcpy(ftgl.txt, res);
    Tcl_DStringFree(&ds);
  }
  
  if(argc > 2)
    ftgl.x = atof(argv[2]);
  
  if(argc > 3)
    ftgl.y = atof(argv[3]);
  
  user_hook_after_all = FtRender;

  Tcl_Eval(interp, "display");
  
  return TCL_OK;
}

#endif // HAVE_FTGL_NEW || HAVE_FTGL_OLD
