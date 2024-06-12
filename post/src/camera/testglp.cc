/*
 * "fogpot.c" - A test program demonstrating the use of glFog() for
 *              depth-cueing.
 */

#include <GL/glaux.h>
#include <GL/glpfile.h>


/*
 * These #define constants are provided for compatibility between MS Windows
 * and the rest of the world.
 *
 * CALLBACK and APIENTRY are function modifiers under MS Windows.
 */

#ifndef WIN32
#  define CALLBACK
#  define APIENTRY
#endif /* !WIN32 */

GLfloat rotation = 0.0;
GLPfile	*output;


/*
 * 'reshape_scene()' - Change the size of the scene...
 */

void CALLBACK
reshape_scene(GLsizei width,	/* I - Width of the window in pixels */
              GLsizei height)	/* I - Height of the window in pixels */
{
 /*
  * Reset the current viewport and perspective transformation...
  */

  glViewport(0, 0, width, height);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(22.5, (float)width / (float)height, 0.1, 1000.0);

  glMatrixMode(GL_MODELVIEW);
}


/*
 * 'draw_scene()' - Draw a scene containing a cube with a sphere in front of
 *                  it.
 */

void CALLBACK
draw_scene(void)
{
  static float	red_light[4] = { 1.0, 0.0, 0.0, 1.0 };
  static float	red_pos[4] = { 1.0, 1.0, 1.0, 0.0 };
  static float	blue_light[4] = { 0.0, 0.0, 1.0, 1.0 };
  static float	blue_pos[4] = { -1.0, -1.0, -1.0, 0.0 };
  static float	fog_color[4] = { 0.0, 0.0, 0.0, 0.0 };


 /*
  * Enable drawing features that we need...
  */

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);

  glShadeModel(GL_SMOOTH);
//  glShadeModel(GL_FLAT);

  if (output != NULL)
    output->StartPage();

 /*
  * Clear the color and depth buffers...
  */

  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

 /*
  * Draw the cube and sphere in different colors...
  *
  * We have positioned two lights in this scene.  The first is red and
  * located above, to the right, and behind the viewer.  The second is blue
  * and located below, to the left, and in front of the viewer.
  */

  glLightfv(GL_LIGHT0, GL_DIFFUSE, red_light);
  glLightfv(GL_LIGHT0, GL_POSITION, red_pos);

  glLightfv(GL_LIGHT1, GL_DIFFUSE, blue_light);
  glLightfv(GL_LIGHT1, GL_POSITION, blue_pos);

  glEnable(GL_COLOR_MATERIAL);

//  glEnable(GL_FOG);
  glFogf(GL_FOG_DENSITY, 0.25);
  glFogfv(GL_FOG_COLOR, fog_color);

  glPushMatrix();
    glTranslatef(-1.0, -0.5, -15.0);
    glRotatef(-rotation, 0.0, 1.0, 0.0);
    glRotatef(rotation, 1.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 0.0);

    auxSolidCylinder(1.0, 1.0);
  glPopMatrix();

  glPushMatrix();
    glTranslatef(1.0, 0.0, -10.0);
    glRotatef(rotation, 0.0, 1.0, 0.0);
    glRotatef(rotation, 1.0, 0.0, 0.0);

    glColor3f(0.0, 1.0, 1.0);
    auxSolidTeapot(1.0);
  glPopMatrix();

  if (output != NULL)
  {
    if (output->EndPage() == 0)
    {
      delete output;
      output = NULL;
    };
  };

  auxSwapBuffers();
}


/*
 * 'rotate_objects()' - Rotate while we are idle...
 */

void CALLBACK
rotate_objects(void)
{
  rotation += 2.0;
  if (rotation >= 360.0)
    rotation -= 360.0;

  draw_scene();
}


/*
 * 'main()' - Initialize the window and display the scene until the user presses
 *            the ESCape key.
 */

void
main(void)
{
  auxInitDisplayMode(AUX_RGB | AUX_DEPTH | AUX_DOUBLE);
  auxInitPosition(256, 256, 768, 512);
  auxInitWindow("Fogged Teapots");

  auxReshapeFunc(reshape_scene);
  auxIdleFunc(rotate_objects);

  output = new GLPfile("test.ps");

  auxMainLoop(draw_scene);
}


/*
 * End of "fogpot.c".
 */
