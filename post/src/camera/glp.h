//
// "$Id: glp.h,v 1.3 2005/05/31 11:28:13 vierinen Exp $"
//
//   Header file for the GLP library, an OpenGL printing toolkit.
//
//   The GLP library is distributed under the terms of the GNU Library
//   General Public License which is described in the file "COPYING.LIB".
//   If you use this library in your program, please include a line reading
//   "OpenGL Printing Toolkit by Michael Sweet" in your version or copyright
//   output.
//
// Revision History:
//
//   $Log: glp.h,v $
//   Revision 1.3  2005/05/31 11:28:13  vierinen
//   ads
//
//   Revision 1.2  2005/05/31 10:39:03  vierinen
//   apple?
//
//   Revision 1.1.1.1  2005/05/31 06:29:21  vierinen
//   ads
//
//   Revision 1.1  2003/02/06 09:37:54  jpr
//   *** empty log message ***
//
//   2003/02/01  Juha Ruokolainen CSC
//   Added protos for binary tree handling in GLPcontex
//   Revision 1.2  1996/07/13  12:52:02  mike
//   Changed the public methods to 'virtual'.
//
//   Revision 1.1  1996/06/27  03:07:13  mike
//   Initial revision
//

#ifndef _GL_GLP_H_
#  define _GL_GLP_H_
#include "../../config.h"

#if defined(WIN32) ||defined(MINGW32)
#include <windows.h>
#include <direct.h>
#include <wingdi.h>
#include <windef.h>
#endif


//
// Include necessary headers.
//
#include <GL/gl.h>
#include "../../config.h"



#  include <iostream>
#  include <fstream>


//
// Printing options...
//

#  define GLP_FIT_TO_PAGE	1	// Fit the output to the page
#  define GLP_AUTO_CROP		2	// Automatically crop to geometry
#  define GLP_GREYSCALE		4	// Output greyscale rather than color
#  define GLP_REVERSE		8	// Reverse grey shades
#  define GLP_DRAW_BACKGROUND	16	// Draw the background color

//
// OpenGL configuration options...
//

#  define GLP_RGBA		0	// RGBA mode window
#  define GLP_COLORINDEX	1	// Color index mode window

//
// Error codes...
//

#  define GLP_SUCCESS		0	// Success - no error occurred
#  define GLP_OUT_OF_MEMORY	-1	// Out of memory
#  define GLP_NO_FEEDBACK	-2	// No feedback data available
#  define GLP_ILLEGAL_OPERATION	-3	// Illegal operation of some kind


//
// Various structures used for sorting feedback data prior to printing...
//

typedef GLfloat GLPrgba[4];	// GLPrgba array structure for sanity
typedef GLfloat GLPxyz[3];	// GLPxyz array structure for sanity

struct GLPvertex	//// GLPvertex structure
{
  GLPxyz	xyz;		// Location of vertex
  GLPrgba	rgba;		// Color of vertex (alpha may be used later)
};

struct GLPprimitive	//// GLPprimitive structure
{
  GLPprimitive	*left, *right;  // left right pointers of the depth sort tree
  GLboolean	shade;		// GL_TRUE if this primitive should be shaded
  GLfloat	zmin, zmax;	// Min and max depth values
  int		num_verts;	// Number of vertices used
  GLPvertex	verts[3];	// Up to 3 vertices
};

struct GLPbbox		//// GLPbbox structure
{
  GLPbbox	*next,		// Next bounding box in list
		*prev;		// Previous bounding box in list
  GLPprimitive	*primitives,	// Primitives inside this box
                *lastprim;
  GLfloat	min[3],		// Minimum X, Y, Z coords
		max[3];		// Maximum X, Y, Z coords
};

//
// The GLPcontext class provides all the basic functionality to support
// OpenGL feedback-based printing to vector/polygon printing devices or
// file formats.  For raster-only devices you are probably better off with
// an off-screen bitmap.
//

class GLPcontext	//// GLPcontext class
{
  protected:
	  int		options;	// Printing options
	  GLPbbox	*bboxes;	// Primitive data
	  GLPprimitive  *tree;		// Primitive data
	  int		feedsize;	// Feedback buffer size
	  GLfloat	*feedback;	// Feedback data
	  int		feedmode;	// Feedback mode (RGBA or colormap)
	  int		colorsize;	// Colormap size
          GLPrgba	*colormap;	// Colorindex mapping to RGBA vals

          void add_primitive(GLboolean depth, GLboolean shade,
                             int num_verts, GLPvertex *verts);
          void sort_primitive(GLboolean depth, GLPbbox *bbox,
                              GLPprimitive *newprim);
          void add_tree( GLPprimitive **tree, GLPprimitive *prim);
          void add_subtree( GLPprimitive **tree, GLPprimitive *prim);
          int  get_vertex(GLPvertex *v, GLfloat *p);
          void delete_all(void);
          void delete_tree(GLPprimitive *tree);

  public:
          virtual ~GLPcontext(void);

          virtual int StartPage(int mode = GLP_RGBA);
          virtual int StartPage(int     mode,
                        	int     size,
                        	GLPrgba *rgba);
	  virtual int UpdatePage(GLboolean more);
	  virtual int EndPage(void);

	  virtual void SetOptions(int print_options);
};

#endif // !_GL_GLP_H_

//
// End of "$Id: glp.h,v 1.3 2005/05/31 11:28:13 vierinen Exp $".
//
