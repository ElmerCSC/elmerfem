//
// "$Id: context.cc,v 1.1 2003/02/06 09:37:53 jpr Exp $"
//
//   GLPcontext class functions for the GLP library, an OpenGL printing
//   toolkit.
//
//   The GLP library is distributed under the terms of the GNU Library
//   General Public License which is described in the file "COPYING.LIB".
//   If you use this library in your program, please include a line reading
//   "OpenGL Printing Toolkit by Michael Sweet" in your version or copyright
//   output.
//
// Contents:
//
//   ~GLPcontext, StartPage, UpdatePage, EndPage,
//   add_primitive, sort_primitive, get_vertex, delete_all
//
// Revision History:
//
//   $Log: context.cc,v $
//   Revision 1.1  2003/02/06 09:37:53  jpr
//   *** empty log message ***
//
//   2003/02/01 Juha Ruokolainen / CSC - IT Center for Science Ltd.
//   Fixed some bugs in add_primitive. Changed the list of primitives 
//   to a binary tree of primitives. Sorting now works correctly.
//   Added handling of GL_LINE_RESET_TOKEN, but the line patterns etc. are
//   still silently ignored.
//
//   Revision 1.2  1996/07/13  12:52:02  mike
//   Fixed delete_all() - was not setting bboxes list pointer to NULL at end.
//
//   Revision 1.1  1996/06/27  00:58:11  mike
//   Initial revision
//

//
// Include necessary headers.
//

#include <../../config.h>

#include "glp.h"
#include <math.h>
#include <stdlib.h>
#ifdef HAVE_STRING_H
#include <string.h>
#endif

//#define DEBUG


#define min(a,b)	((a) < (b) ? (a) : (b))
#define max(a,b)	((a) > (b) ? (a) : (b))


//
// GLPcontext destructor; like the constructor, this just does the basics.
// You'll want to implement your own...
//

GLPcontext :: ~GLPcontext(void)
{
  // Free any memory we've allocated...

  delete_all();

//  if (feedback != NULL)
//    delete feedback;
//  if (colormap != NULL)
//    delete colormap;
}


//
// 'StartPage' function for RGBA windows; this should perform any necessary
// output initialization (e.g. send 'start of page' commands, etc) and enable
// feedback mode...
//

int
GLPcontext :: StartPage(int mode)
{
  // Initialize feedback mode...

  feedmode = mode;

#ifdef DEBUG
  cout << "RGBA feedmode = " << mode << "\n" << flush;
#endif /* DEBUG */

  glFeedbackBuffer(feedsize, GL_3D_COLOR, feedback);
  glRenderMode(GL_FEEDBACK);

  // You'd put any other 'start page' things here...

  return (0);
}


//
// 'StartPage' function for color index windows; this does the same thing as
// the RGBA start page function, and also copies the given colormap into our
// class colormap structure for later use...
//

int
GLPcontext :: StartPage(int     mode,
                        int     size,
                        GLPrgba *rgba)
{
  // Initialize feedback mode...

  feedmode = mode;

#ifdef DEBUG
  cout << "Index feedmode = " << mode << "\n" << flush;
#endif /* DEBUG */

  glFeedbackBuffer(feedsize, GL_3D_COLOR, feedback);
  glRenderMode(GL_FEEDBACK);

  // Copy the colormap over, removing the old one as necessary...

  if (colormap != NULL)
    delete colormap;

  colorsize = size;
  colormap  = new GLPrgba[size];
  memcpy(colormap, rgba, size * sizeof(GLPrgba));

  // You'd put any other 'start page' things here...

  return (0);
}


//
// 'UpdatePage' function; this does most of the dirty work, adding feedback
// data to the current primitive list.
//
// If the 'more' argument is TRUE then the current rendering mode is put back
// into 'GL_FEEDBACK' mode...
//
// Normally you won't redefine this function...
//
#include <stdio.h>

int
GLPcontext :: UpdatePage(GLboolean more)
{
  int		i, used, count, v;
  GLfloat	*current;
  GLPvertex	vertices[3];
  GLboolean	depth,
		shade;
  GLint		shademodel;


#ifdef DEBUG
  cout << "UpdatePage(" << (more ? "GL_TRUE" : "GL_FALSE") << ")\n" << flush;
#endif /* DEBUG */

  // Get the current depth test state and shade model...

  depth = glIsEnabled(GL_DEPTH_TEST);
  glGetIntegerv(GL_SHADE_MODEL, &shademodel);
  shade = shademodel == GL_SMOOTH;

  // Figure out how many feedback entries there are and put the current
  // OpenGL context back in feedback mode if 'more' is true...

  used = glRenderMode(more ? GL_FEEDBACK : GL_RENDER);
fprintf( stderr, "yseD: %d\n", used );
  if (used <= 0) return (GLP_NO_FEEDBACK);

#ifdef DEBUG
  cout << "glGetError() after glRenderMode returned " << glGetError() << "\n" << flush;
  cout << "First: used = " << used << ", feedback[0] = " << feedback[0] << "\n" << flush;
#endif /* DEBUG */

  // Parse the feedback buffer...

  current = feedback;
  while (used > 0)
  {
#ifdef DEBUG
    cout << "Loop: used = " << used << "\n" << flush;
#endif /* DEBUG */

    switch ((int)*current)
    {
      case GL_POINT_TOKEN :
          current ++;
          used --;
          i = get_vertex(vertices + 0, current);
          current += i;
          used    -= i;
          add_primitive(depth, shade, 1, vertices);
          break;
      case GL_LINE_TOKEN: case GL_LINE_RESET_TOKEN:
          current ++;
          used --;
          i = get_vertex(vertices + 0, current);
          current += i;
          used    -= i;
          i = get_vertex(vertices + 1, current);
          current += i;
          used    -= i;
          add_primitive(depth, shade, 2, vertices);
          break;

      case GL_POLYGON_TOKEN :
          // Get the number of vertices...

          count = (int)current[1];

          current += 2;
          used -= 2;

          // Loop through and add a series of triangles...

          v = 0;
          while (count > 0 && used > 0)
          {
            i = get_vertex( vertices + v, current );
            current += i;
            used    -= i;
            count --;

            // Add a triangle if we have 3 vertices...

            if (v == 2)
            {
              add_primitive( depth, shade, 3, vertices );
              vertices[1] = vertices[2];
              v = 0;
            }
            else
              v++;
          };
          break;

      case GL_BITMAP_TOKEN :
      case GL_DRAW_PIXEL_TOKEN :
      case GL_COPY_PIXEL_TOKEN :
          current ++;
          used --;
          i = get_vertex(vertices + 0, current);
          current += i;
          used    -= i;
          break;

      case GL_PASS_THROUGH_TOKEN :
#ifdef DEBUG
          std::cout << "UpdatePage: Ignoring passthrough token " << current[1] << "...\n" << std::flush;
#endif /* DEBUG */
          current += 2;
          used -= 2;
          break;

      default :
          std::cout << "UpdatePage: Ignoring unknwon token " << current[0] << "...\n" << std::flush;
          current ++;
          used --;
          break;
    };
  };

  return (GLP_SUCCESS);
}


//
// 'EndPage' function; this does nothing except parse the bounding boxes
// and output any remaining primitives.  It then frees the bounding box
// and primitive lists...
//

int
GLPcontext :: EndPage(void)
{
  GLPbbox	*bbox;				// Current bounding box
  GLPprimitive	*prim;				// Current primitive


#ifdef DEBUG
  cout << "EndPage()\n" << flush;
#endif /* DEBUG */

  // Stop doing feedback...

  UpdatePage(GL_FALSE);

  if (bboxes == NULL) return (GLP_NO_FEEDBACK);

  // Loop through all bounding boxes and primitives...
  for (bbox = bboxes; bbox != NULL; bbox = bbox->next)
  {
  };

  // Delete everything from the bounding box and primitive lists...

  delete_all();

#ifdef DEBUG
  cout << "EndPage() - done.\n" << flush;
#endif /* DEBUG */

  return (GLP_SUCCESS);
}


//
// 'SetOptions' function; this just sets the 'options' member to whatever is
// passed in.
//

void
GLPcontext :: SetOptions(int print_options)	// I - New printing options
{
  options = print_options;
}


//
// 'add_primitive' function; add a primitive to the list of primitives and
// bounding boxes for the current context.
//

void
GLPcontext :: add_primitive(GLboolean depth,	// I - Depth testing enabled?
                            GLboolean shade,	// I - Smooth shading?
                            int       num_verts,// I - Number of vertices
                            GLPvertex *verts)	// I - Vertices
{
  int		i,				// Looping var
		count;				// Count of intersections
  GLfloat	min[3],				// Minimum (x,y) coords
		max[3];				// Maximum (x,y) coords
  GLPprimitive	*newprim,*pprim,*ppprim;	// New primitive
  GLPbbox	*bbox,				// Current bounding box
		*joinbbox,			// Bounding box to join with
		*nextbbox;			// Next bounding box...


#ifdef DEBUG
  cout << "add_primitive(" << (depth ? "GL_TRUE" : "GL_FALSE") << ", "
       << (shade ? "GL_TRUE" : "GL_FALSE") << ", "
       << num_verts << ", "
       << (int)verts << ")\n" << flush;
#endif /* DEBUG */

  // First create the new primitive and compute the bounding box for it...

  newprim = new GLPprimitive;

  newprim->left  = NULL;
  newprim->right = NULL;
  newprim->shade     = shade;
  newprim->num_verts = num_verts;
  memcpy( newprim->verts, verts, sizeof(GLPvertex) * num_verts );

  min[0] = min[1] = min[2] =  1e20;
  max[0] = max[1] = max[2] = -1e20;

  for (i = 0; i < num_verts; i ++)
  {
    if ( verts[i].xyz[0] < min[0] ) min[0] = verts[i].xyz[0];
    if ( verts[i].xyz[1] < min[1] ) min[1] = verts[i].xyz[1];
    if ( verts[i].xyz[2] < min[2] ) min[2] = verts[i].xyz[2];

    if ( verts[i].xyz[0] > max[0] ) max[0] = verts[i].xyz[0];
    if ( verts[i].xyz[1] > max[1] ) max[1] = verts[i].xyz[1];
    if ( verts[i].xyz[2] > max[2] ) max[2] = verts[i].xyz[2];
  };

  newprim->zmin    = min[2];
  newprim->zmax    = max[2];

  // Stretch the bbox out to the nearest 64 pixels to improve performance...
  min[0] = floor( min[0] * 0.015625 ) * 64.0;
  min[1] = floor( min[1] * 0.015625 ) * 64.0;
  max[0] = ceil(  max[0] * 0.015625 ) * 64.0;
  max[1] = ceil(  max[1] * 0.015625 ) * 64.0;

  // Now search the current bounding box list to see if this primitive lies
  // inside an existing bbox, partially inside, or completely outside...
  //
  // The 'count' variable counts the number of corners that lie inside the
  // current bounding box.  If 'count' is 0, then the bbox is completely
  // outside the current bbox.  Values between 1 and 3 indicate a partially
  // inside primitive.  A value of 4 means that the primitive is completely
  // inside the bbox.
  //
  // If the primitive lies completely outside any bboxes that are out there
  // already, then a new bbox is created with the primitive in it.
  //
  // If the primitive lies partially inside the bbox, the bbox is expanded to
  // include the primitive, and a 'join' operation is performed with any
  // neighboring bboxes that intersect with the expanded bbox.  Finall, the
  // primitive is added to the bbox using the 'sort_primitive()' function
  // (this handles depth buffering if enabled).
  //
  // If the primitive lies completely inside the bbox, it is added with
  // 'sort_primitive()'.

  for (bbox = bboxes; bbox != NULL; bbox = bbox->next)
  {
    count = 0;

    if ( min[0] > bbox->min[0] && min[0] < bbox->max[0] ) count++;
    if ( max[0] > bbox->min[0] && max[0] < bbox->max[0] ) count++;
    if ( min[1] > bbox->min[1] && min[1] < bbox->max[1] ) count++;
    if ( max[1] > bbox->min[1] && max[1] < bbox->max[1] ) count++;

    if ( count > 0 ) break;
  };

  if (bbox == NULL)
  {
    // New bbox...

    bbox = new GLPbbox;

    bbox->prev = NULL;
    bbox->next = bboxes;
    if (bboxes != NULL) bboxes->prev = bbox;
    bboxes = bbox;

    bbox->min[0]     = min[0];
    bbox->max[0]     = max[0];
    bbox->min[1]     = min[1];
    bbox->max[1]     = max[1];
    bbox->min[2]     = min[2];
    bbox->max[2]     = max[2];
    bbox->primitives = newprim;
    bbox->lastprim   = newprim;
  } 
  else if (count < 4)
  {
    // Partially inside...

    if ( min[0] < bbox->min[0] ) bbox->min[0] = min[0];
    if ( max[0] > bbox->max[0] ) bbox->max[0] = max[0];
    if ( min[1] < bbox->min[1] ) bbox->min[1] = min[1];
    if ( max[1] > bbox->max[1] ) bbox->max[1] = max[1];

    // Incrementally join bounding boxes until no more boxes are joined...

    do
    {
      count = 0;
      for (joinbbox = bboxes; joinbbox != NULL; joinbbox = nextbbox)
      {
        nextbbox = joinbbox->next;

        if (joinbbox == bbox)
          continue;
        else if (( bbox->min[0] > joinbbox->min[0] && bbox->min[0] < joinbbox->max[0]) ||
        	 ( bbox->max[0] > joinbbox->min[0] && bbox->max[0] < joinbbox->max[0]) ||
        	 ( bbox->min[1] > joinbbox->min[1] && bbox->min[1] < joinbbox->max[1]) ||
        	 ( bbox->max[1] > joinbbox->min[1] && bbox->max[1] < joinbbox->max[1]))
        {
          // Join this one...

          count++;

          if (joinbbox->prev == NULL)
            bboxes = joinbbox->next;
          else
            (joinbbox->prev)->next = joinbbox->next;

          if (nextbbox != NULL)
            nextbbox->prev = joinbbox->prev;

          for (i = 0; i < 3; i ++)
          {
            if ( joinbbox->min[i] < bbox->min[i] ) bbox->min[i] = joinbbox->min[i];
            if ( joinbbox->max[i] > bbox->max[i] ) bbox->max[i] = joinbbox->max[i];
          };

          add_subtree( &bbox->primitives, joinbbox->primitives );
          delete joinbbox;
        };
      };
    }
    while (count > 0);

    // Add the primitive to this bbox...

    add_tree( &bbox->primitives, newprim);
  }
  else
  {
    // Primitive lies completely inside the bbox, so just add it...
    add_tree( &bbox->primitives, newprim);
  };
}


//
// 'get_vertex' function; get a vertex from the feedback buffer...
//

int
GLPcontext :: get_vertex(GLPvertex *v,	// O - Vertex pointer
                         GLfloat   *p)	// I - Data pointer
{
  int	i;				// Color index


  v->xyz[0] = p[0];
  v->xyz[1] = p[1];
  v->xyz[2] = p[2];

#ifdef DEBUG
  cout << "{ " << p[0] << ", " << p[1] << ", " << p[2] << "}, " << flush;
#endif /* DEBUG */

  if (feedmode == GL_COLOR_INDEX && colorsize > 0)
  {
    // Color index value...
    i = (int)(p[3] + 0.5);

    v->rgba[0] = colormap[i][0];
    v->rgba[1] = colormap[i][1];
    v->rgba[2] = colormap[i][2];
    v->rgba[3] = colormap[i][3];

    return (4);
  }
  else
  {
    // RGBA value...

    v->rgba[0] = p[3];
    v->rgba[1] = p[4];
    v->rgba[2] = p[5];
    v->rgba[3] = p[6];

    return (7);
  };
}


//
// 'delete_all' function; delete all bounding boxes and primitives from
// the current context.
//

void GLPcontext :: delete_tree( GLPprimitive *tree )
{
   if ( !tree ) return;
   delete_tree( tree->left ); 
   delete_tree( tree->right ); 
   delete tree;
}

void
GLPcontext :: delete_all(void)
{
  GLPbbox	*bbox,				// Current bounding box
		*nextbbox;			// Next bounding box
  GLPprimitive	*prim,				// Current primitive
		*nextprim;			// Next primitive


#ifdef DEBUG
  cout << "delete_all()\n" << flush;
#endif /* DEBUG */

  for (bbox = bboxes, nextbbox = NULL; bbox != NULL; bbox = nextbbox)
  {
    delete_tree( bbox->primitives );
    nextbbox = bbox->next;
    delete bbox;
  };

  bboxes = NULL;

  delete_tree( tree );
  tree = NULL;

#ifdef DEBUG
  cout << "delete_all(): done.\n" << flush;
#endif /* DEBUG */
}


void GLPcontext::add_subtree( GLPprimitive **ptree, GLPprimitive *prim )
{
   GLPprimitive *left, *right;

   if ( !prim ) return;

   left  = prim->left;
   right = prim->right;

   prim->left = NULL;
   prim->right = NULL;
   add_tree( ptree, prim );

   add_subtree( ptree, left );
   add_subtree( ptree, right );
}

void GLPcontext::add_tree( GLPprimitive **ptree, GLPprimitive *prim )
{
   static int count=0;
   if ( !*ptree )
   {
      *ptree = prim;
   } else {
      if ( (prim->zmin+prim->zmax) == ((*ptree)->zmin+(*ptree)->zmax) )
      {
           if ( count++ ) {
              prim->right = (*ptree)->left;
              (*ptree)->left = prim;
              count = 0;
           } else { 
              prim->left = (*ptree)->right;
              (*ptree)->right = prim;
           }
      } else if ( (prim->zmin+prim->zmax) > ((*ptree)->zmin+(*ptree)->zmax) )
      {
         if ( (*ptree)->left ) 
           add_tree( &(*ptree)->left, prim );
         else
           (*ptree)->left = prim;
      } else {
         if ( (*ptree)->right ) 
           add_tree( &(*ptree)->right, prim );
         else
           (*ptree)->right = prim;
      }
   }
}

//
// End of "$Id: context.cc,v 1.1 2003/02/06 09:37:53 jpr Exp $".
//
