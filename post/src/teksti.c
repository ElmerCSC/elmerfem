#include "elmerpost.h"

#include <stdio.h>
#include <math.h>

  
#include <GL/gl.h>
#include <GL/glu.h>


#include "tcl.h"

#include <X11/Xlib.h>
#include <X11/Xutil.h>

extern void (*user_hook_before_all)();

static unsigned int FontBase;
static double x,y,z,size;
static char text[132];

void jooMakeRasterFont( char *name )
{
    unsigned int first, last;
   
    XFontStruct *fontInfo;
   
    int n,i;
    static char str[32], here = 0;

    char **FontNames;
    
    fontInfo = XLoadQueryFont( auxXDisplay(),name );
    if ( fontInfo == NULL )
    {
        fprintf( stderr, "Can't find font: [%s]\n", name );
        return;
    }
   
    first = fontInfo->min_char_or_byte2;
    last  = fontInfo->max_char_or_byte2;
   
    FontBase = glGenLists( last + 1 );
    glXUseXFont( fontInfo->fid,first,last-first+1,FontBase+first );
}

void jooPrintString( char *s )
{
   glDisable( GL_LIGHTING );
   glDisable( GL_TEXTURE_1D );
   glPushAttrib( GL_LIST_BIT );
   glListBase( FontBase );
   glCallLists( strlen(s),GL_UNSIGNED_BYTE,(unsigned char *)s);
   glPopAttrib();
   glEnable( GL_TEXTURE_1D );
   glEnable( GL_LIGHTING );
}

void jooShowString( char *font,char *str,float x,float y,float z )
{
   static int been_here = 0;

   glMatrixMode( GL_MODELVIEW );
   glPushMatrix();
   glLoadIdentity();
   
   glMatrixMode( GL_PROJECTION );
        glPushMatrix();
   glLoadIdentity();

   /*
   glDisable( GL_TEXTURE_1D );
   glDisable( GL_DEPTH_TEST );
   glDrawBuffer( GL_FRONT_AND_BACK );
   */
   
   glColor3f( 0.0,0.0,0.0 );
   glRasterPos3f( x,y,z );
   
   if ( !been_here ) 
   {
      jooMakeRasterFont( font );
       been_here = 1;
   }
   jooPrintString( str );
   
   glPopMatrix();

   glMatrixMode( GL_MODELVIEW );
   glPopMatrix();
}

void jooDoit()
{

   char str[132];

   sprintf(str,"-adobe-helvetica-bold-r-normal--%g--100-100-p-88-iso8859-1",size);  
   jooShowString(str,text,x,y,z);

}

static int Teksti( ClientData cl,Tcl_Interp *interp,int argc,char **argv )
{
   extern double br,bg,bb;

   user_hook_before_all = jooDoit;

   x = atof(argv[1]);
   y = atof(argv[2]);
   z = atof(argv[3]);

   strcpy(text,argv[4]);

   size = atof(argv[5]);

   return TCL_OK;
}

int Teksti_Init( Tcl_Interp *interp )
{
   Tcl_CreateCommand( interp,"teksti",Teksti,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);

   return TCL_OK;
}
