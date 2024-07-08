//
// "$Id: file.cc,v 1.1.1.1 2005/05/31 06:29:21 vierinen Exp $"
//
//   GLPfile class function for the GLP library, an OpenGL printing toolkit.
//
//   The GLP library is distributed under the terms of the GNU Library
//   General Public License which is described in the file "COPYING.LIB".
//   If you use this library in your program, please include a line reading
//   "OpenGL Printing Toolkit by Michael Sweet" in your version or copyright
//   output.
//
// Contents:
//
//   GLPfile, ~GLPfile, EndPage
//
// Revision History:
//
//   $Log: file.cc,v $
//   Revision 1.1.1.1  2005/05/31 06:29:21  vierinen
//   ads
//
//   Revision 1.2  2004/11/29 08:26:22  jpr
//   *** empty log message ***
//
//   Revision 1.1  2003/02/06 09:37:53  jpr
//   *** empty log message ***
//
//   2003/02/01 Juha Ruokolainen / CSC - IT Center for Science Ltd., Finland
//   Added code to handle primitives stored in a binary tree instead of
//   a list. Separate routine for outputting primitives. Some changes
//   to generated postscript code: 
//   - better splitting of shaded triangles
//   - if shfill operator is known by the postscript interpreter
//     (level3 postscript) use that, instead of splitting.
//   
//
//   Revision 1.2  1996/07/13  12:52:02  mike
//   Fixed shaded triangle procedure in PostScript prolog (ST).
//   Added auto-rotation for GLP_FIT_TO_PAGE option.
//
//   Revision 1.1  1996/06/27  03:06:36  mike
//   Initial revision
//

//
// Include necessary headers.
//


#include <../../config.h>

#include <stdio.h>
#include <stdlib.h>

#if defined(HAVE_STRINGS_H)
#include <strings.h>
#endif

#if defined(HAVE_STRING_H)
#include <string.h>
#endif

#include "glpfile.h"


double LLx, LLy, URx, URy;
extern GLPtext       *texthead;

//#define DEBUG

//
// GLPfile constructor; this really doesn't do much, but instead you are
// expected to steal this code as the basis of your own GLP subclass...
//

GLPfile :: GLPfile(char *print_name,	// I - Name of file/printer
                   int  print_options,	// I - Output options
                   int  print_size)	// I - Size of feedback buffer
{
  // Initialize private class data and allocate a feedback buffer...

  options    = print_options;
  bboxes     = NULL;
  tree       = NULL;
  feedsize   = print_size > 0 ? print_size : 1024 * 1024 * 100;
  feedback   = new GLfloat[feedsize];
  colorsize  = 0;
  colormap   = NULL;
  page_count = 0;
  texthead   = NULL;

#ifdef DEBUG
  cout << "feedsize = " << feedsize << ", feedback = " << (int)feedback << "\n" << flush;
#endif /* DEBUG */

  // Next open the indicated file...

  outfile = new std::fstream(print_name, std::ios::out);
  if (outfile == NULL)
    return;

  // Output the PostScript file header and prolog...

  *outfile << "%!PS-Adobe-3.0\n";
  *outfile << "%%LanguageLevel: 1\n";
  *outfile << "%%Creator: GLP 0.1 by Michael Sweet (mike@easysw.com)\n";
  *outfile << "%%Pages: (atend)\n";
  *outfile << "%%BoundingBox: (atend)\n";
  *outfile << "%%EndComments\n";

  *outfile << "%%BeginProlog\n";

  if (options & GLP_GREYSCALE)
  {
    *outfile << "% Greyscale color command - r g b C\n";
    *outfile << "/C { 0.0820 mul exch 0.6094 mul add exch 0.3086 mul add neg 1.0 add setgray } bind def\n";
  }
  else
  {
    *outfile << "% RGB color command - r g b C\n";
    *outfile << "/C { setrgbcolor } bind def\n";
  };
  *outfile << "\n";

  *outfile << "% Point primitive - x y r g b P\n";
  *outfile << "/P { C newpath 0.5 0.0 360.0 arc closepath fill } bind def\n";
  *outfile << "\n";

  *outfile << "% Flat-shaded line - x2 y2 x1 y1 r g b L\n";
  *outfile << "/L { C moveto lineto stroke } bind def\n";
  *outfile << "\n";

  *outfile << "% Smooth-shaded line - x2 y2 r2 g2 b2 x1 y1 r1 g1 b1 SL\n";
  *outfile << "/SL {\n";
  *outfile << "	/b1 exch def\n";
  *outfile << "	/g1 exch def\n";
  *outfile << "	/r1 exch def\n";
  *outfile << "	/y1 exch def\n";
  *outfile << "	/x1 exch def\n";
  *outfile << "	/b2 exch def\n";
  *outfile << "	/g2 exch def\n";
  *outfile << "	/r2 exch def\n";
  *outfile << "	/y2 exch def\n";
  *outfile << "	/x2 exch def\n";
  *outfile << "\n";
  *outfile << "	b2 b1 sub abs 0.01 gt\n";
  *outfile << "	g2 g1 sub abs 0.005 gt\n";
  *outfile << "	r2 r1 sub abs 0.008 gt\n";
  *outfile << "     or or {\n";
  *outfile << "		/bm b1 b2 add 0.5 mul def\n";
  *outfile << "		/gm g1 g2 add 0.5 mul def\n";
  *outfile << "		/rm r1 r2 add 0.5 mul def\n";
  *outfile << "		/ym y1 y2 add 0.5 mul def\n";
  *outfile << "		/xm x1 x2 add 0.5 mul def\n";
  *outfile << "\n";
  *outfile << "		save x1 y1 r1 g1 b1 xm ym rm gm bm SL restore\n";
  *outfile << "		save xm ym rm gm bm x2 y2 r2 g2 b2 SL restore\n";
  *outfile << "	} {\n";
  *outfile << "		x1 y1 x2 y2 r1 g1 b1 L\n";
  *outfile << "	} ifelse\n";
  *outfile << "} bind def\n";
  *outfile << "\n";

  *outfile << "% Flat-shaded triangle - x3 y3 x2 y2 x1 y1 r g b T\n";
  *outfile << "/T { C newpath moveto lineto lineto closepath fill } bind def\n";
  *outfile << "\n";

  *outfile << "% Smooth-shaded triangle - x3 y3 r3 g3 b3 x2 y2 r2 g2 b2 x1 y1 r1 g1 b1 ST\n";
  *outfile << "/ST {\n";
  *outfile << "	/b1 exch def\n";
  *outfile << "	/g1 exch def\n";
  *outfile << "	/r1 exch def\n";
  *outfile << "	/y1 exch def\n";
  *outfile << "	/x1 exch def\n";
  *outfile << "	/b2 exch def\n";
  *outfile << "	/g2 exch def\n";
  *outfile << "	/r2 exch def\n";
  *outfile << "	/y2 exch def\n";
  *outfile << "	/x2 exch def\n";
  *outfile << "	/b3 exch def\n";
  *outfile << "	/g3 exch def\n";
  *outfile << "	/r3 exch def\n";
  *outfile << "	/y3 exch def\n";
  *outfile << "	/x3 exch def\n";
  *outfile << "\n";
  *outfile << "systemdict /shfill known {\n";
  *outfile << " newpath \n";
  *outfile << "<<\n";
  *outfile << "  /ShadingType 4\n";
  *outfile << "  /ColorSpace [/DeviceRGB]\n";
  *outfile << "  /DataSource [0 x1 y1 r1 g2 b1 0 x2 y2 r2 g2 b2 0 x3 y3 r3 g3 b3]\n";
  *outfile << ">> shfill\n } {\n";
  *outfile << "	b2 b1 sub abs 0.05 gt\n";
  *outfile << "	g2 g1 sub abs 0.017 gt\n";
  *outfile << "	r2 r1 sub abs 0.032 gt\n";
  *outfile << "	b3 b1 sub abs 0.05 gt\n";
  *outfile << "	g3 g1 sub abs 0.017 gt\n";
  *outfile << "	r3 r1 sub abs 0.032 gt\n";
  *outfile << "	b2 b3 sub abs 0.05 gt\n";
  *outfile << "	g2 g3 sub abs 0.017 gt\n";
  *outfile << "	r2 r3 sub abs 0.032 gt\n";
  *outfile << "	or or or or or or or or {\n";
  *outfile << "		/b12 b1 b2 add 0.5 mul def\n";
  *outfile << "		/g12 g1 g2 add 0.5 mul def\n";
  *outfile << "		/r12 r1 r2 add 0.5 mul def\n";
  *outfile << "		/y12 y1 y2 add 0.5 mul def\n";
  *outfile << "		/x12 x1 x2 add 0.5 mul def\n";
  *outfile << "\n";
  *outfile << "		/b13 b1 b3 add 0.5 mul def\n";
  *outfile << "		/g13 g1 g3 add 0.5 mul def\n";
  *outfile << "		/r13 r1 r3 add 0.5 mul def\n";
  *outfile << "		/y13 y1 y3 add 0.5 mul def\n";
  *outfile << "		/x13 x1 x3 add 0.5 mul def\n";
  *outfile << "\n";
  *outfile << "		/b32 b3 b2 add 0.5 mul def\n";
  *outfile << "		/g32 g3 g2 add 0.5 mul def\n";
  *outfile << "		/r32 r3 r2 add 0.5 mul def\n";
  *outfile << "		/y32 y3 y2 add 0.5 mul def\n";
  *outfile << "		/x32 x3 x2 add 0.5 mul def\n";
  *outfile << "\n";
  *outfile << "		/b132 b1 b3 b2 add add 3 div def\n";
  *outfile << "		/g132 g1 g3 g2 add add 3 div def\n";
  *outfile << "		/r132 r1 r3 r2 add add 3 div def\n";
  *outfile << "		/y132 y1 y3 y2 add add 3 div def\n";
  *outfile << "		/x132 x1 x3 x2 add add 3 div def\n";
  *outfile << "\n";
  *outfile << "		save x1 y1 r1 g1 b1 x12 y12 r12 g12 b12 x132 y132 r132 g132 b132 ST restore\n";
  *outfile << "		save x1 y1 r1 g1 b1 x13 y13 r13 g13 b13 x132 y132 r132 g132 b132 ST restore\n";
  *outfile << "		save x2 y2 r2 g2 b2 x12 y12 r12 g12 b12 x132 y132 r132 g132 b132 ST restore\n";
  *outfile << "		save x2 y2 r2 g2 b2 x32 y32 r32 g32 b32 x132 y132 r132 g132 b132 ST restore\n";
  *outfile << "		save x3 y3 r3 g3 b3 x13 y13 r13 g13 b13 x132 y132 r132 g132 b132 ST restore\n";
  *outfile << "		save x3 y3 r3 g3 b3 x32 y32 r32 g32 b32 x132 y132 r132 g132 b132 ST restore\n";
  *outfile << "	} {\n";
  *outfile << "		x1 y1 x2 y2 x3 y3 r1 g1 b1 T\n";
  *outfile << "	} ifelse } ifelse\n";
  *outfile << "} bind def\n";
  *outfile << "\n";

  *outfile << "%%EndProlog\n";
}


//
// GLPfile destructor; like the constructor, this just does the basics.
// You'll want to implement your own...
//

GLPfile :: ~GLPfile(void)
{
  // Free any memory we've allocated...
  GLPtext *text, *nexttext;

  for( text=texthead; text; text=nexttext )
  {
     nexttext = text->next;
     free( text->string );
     delete text;
  }
  texthead = NULL;

  delete_all();

  if (feedback != NULL) delete feedback;
  if (colormap != NULL) delete colormap;

  *outfile << "%%Pages: " << page_count << "\n";
  *outfile << "%%BoundingBox: " << LLx << " " << LLy << " " << URx << " " << URy << "\n";
  *outfile << "%%EOF\n";
  delete outfile;
}


//
// 'EndPage' function; this does nothing except parse the bounding boxes
// and output any remaining primitives.  It then frees the bounding box
// and primitive lists...
//

int
GLPfile :: EndPage(void)
{
  GLPbbox	*bbox;				// Current bounding box
  GLPprimitive	*prim;				// Current primitive
  int		i;				// Color index background...
  GLfloat	rgba[4];			// RGBA background...
  GLint		viewport[4];

  GLPtext *text;

#ifdef DEBUG
  cout << "EndPage()\n" << flush;
#endif /* DEBUG */

  // Stop doing feedback...

  UpdatePage(GL_FALSE);

  // Return an error if there was no feedback data used...

    if (bboxes == NULL) return (GLP_NO_FEEDBACK);

  // Loop through all bounding boxes and primitives...
#define MAX(x,y) ((x)>(y)?(x):(y))


#ifdef DEBUG
  cout << "EndPage() - writing page.\n" << flush;
#endif /* DEBUG */

  glGetIntegerv(GL_VIEWPORT, viewport);

  page_count ++;
  *outfile << "%%Page: " << page_count << "\n";
  *outfile << "%%PageBoundingBox: 0 0 " << viewport[2] << " " << viewport[3] << "\n";
  LLx = 0;
  LLy = 0;
  URx = MAX( viewport[2], URx );
  URy = MAX( viewport[3], URy );

  *outfile << "gsave\n";

  if (options & GLP_FIT_TO_PAGE)
  {
    *outfile << "% Fit to page...\n";
    *outfile << "newpath clippath pathbbox\n";
    *outfile << "/URy exch def\n";
    *outfile << "/URx exch def\n";
    *outfile << "/LLy exch def\n";
    *outfile << "/LLx exch def\n";
    *outfile << "/Width  URx LLx sub 0.005 sub def\n";
    *outfile << "/Height URy LLy sub 0.005 sub def\n";
    *outfile << "LLx LLy translate\n";
    *outfile << "/XZoom Width " << viewport[2] << " div def\n";
    *outfile << "/YZoom Height " << viewport[3] << " div def\n";
    *outfile << "Width YZoom mul Height XZoom mul gt {\n";
    *outfile << "	90 rotate\n";
    *outfile << "	0 Width neg translate\n";
    *outfile << "	Width Height /Width exch def /Height exch def\n";
    *outfile << "	/XZoom Width " << viewport[2] << " div def\n";
    *outfile << "	/YZoom Height " << viewport[3] << " div def\n";
    *outfile << "} if\n";
    *outfile << "XZoom YZoom gt {\n";
    *outfile << "	/YSize Height def\n";
    *outfile << "	/XSize " << viewport[2] << " YZoom mul def\n";
    *outfile << "	/Scale YZoom def\n";
    *outfile << "} {\n";
    *outfile << "	/XSize Width def\n";
    *outfile << "	/YSize " << viewport[3] << " XZoom mul def\n";
    *outfile << "	/Scale XZoom def\n";
    *outfile << "} ifelse\n";

    *outfile << "Width  XSize sub 2 div\n";
    *outfile << "Height YSize sub 2 div translate\n";
    *outfile << "Scale Scale scale\n";
    *outfile << "\n";
  };

  if (options & GLP_DRAW_BACKGROUND)
  {
    if (feedmode == GL_RGBA || colorsize == 0)
      glGetFloatv(GL_COLOR_CLEAR_VALUE, rgba);
    else
    {
      glGetIntegerv(GL_INDEX_CLEAR_VALUE, &i);
      rgba[0] = colormap[i][0];
      rgba[1] = colormap[i][1];
      rgba[2] = colormap[i][2];
    };

    *outfile << "% Draw background...\n";
    *outfile << rgba[0] << " " << rgba[1] << " " << rgba[2] << " C\n";
    *outfile << "newpath\n";
    *outfile << "	0 0 moveto\n";
    *outfile << "	" << viewport[2] << " 0 lineto\n";
    *outfile << "	" << viewport[2] << " " << viewport[3] << " lineto\n";
    *outfile << "	0 " << viewport[3] << " lineto\n";
    *outfile << "closepath fill\n";
    *outfile << "\n";
  };

  for (bbox = bboxes; bbox != NULL; bbox = bbox->next) {
         output_tree( bbox->primitives );
  }

  *outfile << "1 1 1 C\n";
  for( text=texthead; text; text = text->next ) DoOutputPS( text );

  *outfile << "grestore\n";
  *outfile << "showpage\n";
  *outfile << "%%EndPage\n";
  *outfile << std::flush;

  // Delete everything from the bounding box and primitive lists...

  delete_all();

#ifdef DEBUG
  cout << "EndPage() - done.\n" << flush;
#endif /* DEBUG */

  return (GLP_SUCCESS);
}


void  GLPfile :: output_tree( GLPprimitive *tree )
{
    if ( !tree ) return;
    output_tree( tree->left );
    output_prim( tree );
    output_tree( tree->right );
}

void  GLPfile :: output_prim( GLPprimitive *prim )
{
      switch (prim->num_verts)
      {
        case 1 : /* Point */
            *outfile << prim->verts[0].xyz[0] << " " << prim->verts[0].xyz[1] << " " <<
            prim->verts[0].rgba[0] << " " << prim->verts[0].rgba[1] << " " <<
            prim->verts[0].rgba[2] << " P\n";
            break;
        case 2 : /* Line */
            if (prim->shade)
            {
              *outfile << prim->verts[1].xyz[0] << " " << prim->verts[1].xyz[1] << " " <<
              prim->verts[1].rgba[0] << " " << prim->verts[1].rgba[1] << " " <<
              prim->verts[1].rgba[2] << " " <<
              prim->verts[0].xyz[0] << " " << prim->verts[0].xyz[1] << " " <<
              prim->verts[0].rgba[0] << " " << prim->verts[0].rgba[1] << " " <<
              prim->verts[0].rgba[2] <<
              " SL\n";
            }
            else
            {
              *outfile << prim->verts[1].xyz[0] << " " << prim->verts[1].xyz[1] << " " <<
              prim->verts[0].xyz[0] << " " << prim->verts[0].xyz[1] << " " <<
              prim->verts[0].rgba[0] << " " << prim->verts[0].rgba[1] << " " <<
              prim->verts[0].rgba[2] <<
              " L\n";
            }
            break;
        case 3 : /* Triangle */
            if (prim->shade)
            {
              *outfile << prim->verts[2].xyz[0] << " " << prim->verts[2].xyz[1] << " " <<
              prim->verts[2].rgba[0] << " " << prim->verts[2].rgba[1] << " " <<
              prim->verts[2].rgba[2] << " " <<
              prim->verts[1].xyz[0] << " " << prim->verts[1].xyz[1] << " " <<
              prim->verts[1].rgba[0] << " " << prim->verts[1].rgba[1] << " " <<
              prim->verts[1].rgba[2] << " " <<
              prim->verts[0].xyz[0] << " " << prim->verts[0].xyz[1] << " " <<
              prim->verts[0].rgba[0] << " " << prim->verts[0].rgba[1] << " " <<
              prim->verts[0].rgba[2] <<
              " ST\n";
            }
            else
            {
              *outfile << prim->verts[2].xyz[0] << " " << prim->verts[2].xyz[1] << " " <<
              prim->verts[1].xyz[0] << " " << prim->verts[1].xyz[1] << " " <<
              prim->verts[0].xyz[0] << " " << prim->verts[0].xyz[1] << " " <<
              prim->verts[0].rgba[0] << " " << prim->verts[0].rgba[1] << " " <<
              prim->verts[0].rgba[2] <<
              " T\n";
            }
            break;
      };
}

void GLPfile :: DoOutputPS( GLPtext *t ) 
{
   *outfile << "/Times-Roman findfont " << t->h << " scalefont setfont\n";
   *outfile << t->r << " " << t->g << " " << t->b << " C\n";
   *outfile << t->x << " " << t->y << " moveto\n";
   *outfile << "(" << t->string << ") show\n";
}

extern "C" void OutputPSString( double x, double y, double h, double r, double g, double b, char *string )
{
   GLPtext *text;

   text = new GLPtext;

   text->next = texthead;
   texthead = text;

   text->x = x;
   text->y = y;
   text->h = h;
   text->r = r;
   text->g = g;
   text->b = b;
   text->string = (char *)malloc( strlen(string)+1 );
   strncpy( text->string, string, strlen(string)+1 );
}

//
// End of "$Id: file.cc,v 1.1.1.1 2005/05/31 06:29:21 vierinen Exp $".
//
