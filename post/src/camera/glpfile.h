//
// "$Id: glpfile.h,v 1.1 2003/02/06 09:37:55 jpr Exp $"
//
//   Header file for the GLPfile class for the OpenGL printing toolkit.
//
//   The GLP library is distributed under the terms of the GNU Library
//   General Public License which is described in the file "COPYING.LIB".
//   If you use this library in your program, please include a line reading
//   "OpenGL Printing Toolkit by Michael Sweet" in your version or copyright
//   output.
//
// Revision History:
//
//   $Log: glpfile.h,v $
//   Revision 1.1  2003/02/06 09:37:55  jpr
//   *** empty log message ***
//
//   Revision 1.2  1996/07/13  12:52:02  mike
//   Changed default print options to include drawing of background.
//
//   Revision 1.1  1996/06/27  03:07:13  mike
//   Initial revision
//

#ifndef _GL_GLPFILE_H_
#define _GL_GLPFILE_H_

//
// Include necessary headers.
//

#include "glp.h"


struct GLPtext
{
   GLPtext *next;
   double x,y,h;
   double r,g,b;
   char *string; 
};

// GLPtext       *texthead;

//
// The Adobe PostScript file output class for GLP.  Basically, we just
// add file and page count members and define our own constructor,
// destructor, and 'EndPage()' functions.
//

class GLPfile : public GLPcontext	//// GLPfile class
{
  private:
	  int		page_count;	// Page count
          std::fstream	*outfile;	// Output file

  public:
          GLPfile(char *print_name,
                  int  print_options = GLP_FIT_TO_PAGE | GLP_DRAW_BACKGROUND,
                  int  print_size = 0);
          ~GLPfile(void);

	  int EndPage(void);
	  void output_tree( GLPprimitive *tree );
	  void output_prim( GLPprimitive *prim );
          void DoOutputPS( GLPtext *text );
};

#endif // !_GL_GLPFILE_H_

//
// End of "$Id: glpfile.h,v 1.1 2003/02/06 09:37:55 jpr Exp $".
//
