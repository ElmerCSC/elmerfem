/*
 * Graphics primitives are loaded dynamically to a global
 * array of functions. To make a graphics driver to a specific
 * device one has to make the following functions (named diffrently
 * of course) and edit routine gra_init to load these functions
 * to array of functions when the device is selected.
 * 
 * gra_open(dev_ident, output_device)
 * - initialize the device to accept graphics commands
 * 
 * gra_close()
 * - deallocate device
 * 
 * gra_clear()
 * - clear screen (current viewport area)
 * 
 * gra_flush()
 * - flush output buffers, if device have such
 * 
 * gra_defcolor(int index, double r, double g, double b)
 * - define entry in color map
 * 
 * gra_color(int index)
 * - select a color from map to be used by following graphics commands
 * 
 * gra_polyline(int n, Point *points)
 * - draw a line between points given
 * 
 * gra_draw(Point *point)
 * - draw a line from current point to a point given
 * 
 * gra_move(Point *point)
 * - move current point to point given
 * 
 * gra_text(height, char *str)
 * - text beginning from current point
 * 
 * gra_polymarker(int index, int n, Point *points)
 * - draw a marker to points given
 * 
 * gra_marker(int index, Point *point)
 * - draw a marker to point given
 * 
 * gra_areafill(int n, Point *points)
 * - polygon filling
 * 
 * gra_image(int width, int height, int depth, int *raster)
 * - raster plotting
 * 
 * provided by the library are (but if your device supports these,
 * include them in the driver if you like).
 *
 * gra_window(double xl, xh, yl, yh, zl, zh);
 * gra_viewport(double wx, wy, wz, *vx, *vy);
 * gra_setmatrix(GMATRIX gm) (4 x 4 transf. matrix)
 * gra_getmatrix(GMATRIX gm)
 * gra_rotate(double x, y, z);
 * gra_scale(double x, y, z);
 * gra_translate(double x, y, z); 
 *
 * following functions can be used as needed:
 * 
 * gra_window_to_viewport(double x,  y, z, *xs, *ys);
 * gra_transm(double xw, yw, zw, xc, yc, zc); 
 */


/*
 * $Id: gra.h,v 1.1.1.1 2005/04/14 13:29:14 vierinen Exp $ 
 *
 * $Log: gra.h,v $
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:41  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#define CL_XMIN -1
#define CL_XMAX  1
#define CL_YMIN -1
#define CL_YMAX  1

#define GRA_DRV_NULL       0 
/* #define GRA_DRV_IRIS    1 */
/* #define GRA_DRV_DISSPLA 2 */
/* #define GRA_DRV_TEKLIB  3 */
#define GRA_DRV_PS         4

/*
 * transformation matrix type
 */
typedef double GMATRIX[4][4];

typedef struct
{
  double xlow, xhigh, ylow, yhigh;
} matc_Rectangle;

typedef struct
{
  double x, y, z;
} Point;

typedef struct
{
  FILE *out_fp;
  int driver;

  struct
  {
    double xlow, xhigh, ylow, yhigh, zlow, zhigh;
  } window;

  matc_Rectangle viewport; 

  GMATRIX modelm;
  GMATRIX viewm;
  GMATRIX projm;
  GMATRIX transfm;
  double  pratio;

  Point cur_point;
  int cur_color;
  int cur_marker;

} G_STATE;

#ifdef MODULE_MATC
G_STATE gra_state = 
{ 
  NULL, GRA_DRV_NULL,                /* out_fp, driver          */
  -1.0, 1.0, -1.0, 1.0, -1.0, 1.0,   /* window                  */
  0.0, 1.0, 0.0, 1.0,                /* viewport                */ 
  1.0, 0.0, 0.0, 0.0,                /*  model  transformation matrix */
  0.0, 1.0, 0.0, 0.0,
  0.0, 0.0, 1.0, 0.0,
  0.0, 0.0, 0.0, 1.0,
  1.0, 0.0, 0.0, 0.0,                /* viewing transformation matrix */
  0.0, 1.0, 0.0, 0.0,
  0.0, 0.0, 1.0, 0.0,
  0.0, 0.0, 0.0, 1.0,
  1.0, 0.0, 0.0, 0.0,                /* projection transformation matrix */
  0.0, 1.0, 0.0, 0.0,
  0.0, 0.0, 1.0, 0.0,
  0.0, 0.0, 0.0, 1.0,
  1.0, 0.0, 0.0, 0.0,                /* total  transformation  matrix */
  0.0, 1.0, 0.0, 0.0,
  0.0, 0.0, 1.0, 0.0,
  0.0, 0.0, 0.0, 1.0,
  0.0,                               /* perspective ratio     */
  0.0,0.0,                           /* cur_point */
  1, 1                               /* cur_color, cur_marker */
};
#else
EXT G_STATE gra_state;
#endif

#define GRA_FUNCS 27

#define	G_OPEN 0
#define G_CLOSE 1
#define G_CLEAR 2
#define G_VIEWPORT 3
#define G_WINDOW 4
#define G_DEFCOLOR 5
#define G_COLOR 6
#define G_POLYLINE 7
#define G_DRAW 8
#define G_MOVE 9
#define G_POLYMARKER 10
#define G_MARKER 11
#define G_AREAFILL 12
#define G_IMAGE 13
#define G_TEXT 14
#define G_FLUSH 15
#define G_RESET 16
#define G_TRANSLATE 17
#define G_ROTATE 18
#define G_SCALE  19
#define G_VIEWPOINT 20 
#define G_GETMATRIX 21
#define G_SETMATRIX 22
#define G_PERSPECTIVE 23
#define G_DBUFFER 24
#define G_SBUFFER 25
#define G_SWAPBUF 26

#ifdef MODULE_MATC
void (*gra_funcs[GRA_FUNCS])() = 
{
  gra_error, gra_error, gra_error, gra_error, gra_error,
  gra_error, gra_error, gra_error, gra_error, gra_error,
  gra_error, gra_error, gra_error, gra_error, gra_error,
  gra_error, gra_error, gra_error, gra_error, gra_error,
  gra_error, gra_error, gra_error, gra_error, gra_error,
  gra_error, gra_error
};
#else
EXT void(*gra_funcs[GRA_FUNCS])();
#endif

#define GRA_OPEN(d) (*gra_funcs[G_OPEN])(d)
#define GRA_CLOSE() (*gra_funcs[G_CLOSE])()
#define GRA_CLEAR() (*gra_funcs[G_CLEAR])()
#define GRA_FLUSH() (*gra_funcs[G_FLUSH])()
#define GRA_RESET() (*gra_funcs[G_RESET])()
#define GRA_DEFCOLOR(i, r, g, b) (*gra_funcs[G_DEFCOLOR])(i, r, g, b)
#define GRA_COLOR(i) (*gra_funcs[G_COLOR])(i)
#define GRA_POLYLINE(n, p) (*gra_funcs[G_POLYLINE])(n, p);
#define GRA_DRAW(p) (*gra_funcs[G_DRAW])(p)
#define GRA_MOVE(p) (*gra_funcs[G_MOVE])(p)
#define GRA_POLYMARKER(i, n, p) (*gra_funcs[G_POLYMARKER])(i, n, p)
#define GRA_MARKER(i, p) (*gra_funcs[G_MARKER])(i, p)
#define GRA_AREAFILL(n, p) (*gra_funcs[G_AREAFILL])(n, p)
#define GRA_IMAGE(w,h,d,r) (*gra_funcs[G_IMAGE])(w,h,d,r)
#define GRA_TEXT(h,r,s) (*gra_funcs[G_TEXT])(h,r,s)
#define GRA_TRANSLATE(x,y,z) (*gra_funcs[G_TRANSLATE])(x,y,z)
#define GRA_ROTATE(x,y,z) (*gra_funcs[G_ROTATE])(x,y,z)
#define GRA_SCALE(x,y,z) (*gra_funcs[G_SCALE])(x,y,z)
#define GRA_VIEWPOINT(xf,yf,zf,xt,yt,zt) (*gra_funcs[G_VIEWPOINT])(xf,yf,zf,xt,yt,zt)
#define GRA_GETMATRIX(gm) (*gra_funcs[G_GETMATRIX])(gm)
#define GRA_SETMATRIX(gm) (*gra_funcs[G_SETMATRIX])(gm)
#define GRA_DBUFFER(gm) (*gra_funcs[G_DBUFFER])(gm)
#define GRA_SBUFFER(gm) (*gra_funcs[G_SBUFFER])(gm)
#define GRA_SWAPBUF(gm) (*gra_funcs[G_SWAPBUF])(gm)
#define GRA_WINDOW(x1,x2,y1,y2,z1,z2) (*gra_funcs[G_WINDOW])(x1,x2,y1,y2,z1,z2)
#define GRA_VIEWPORT(x1,x2,y1,y2) (*gra_funcs[G_VIEWPORT])(x1,x2,y1,y2)
#define GRA_PERSPECTIVE(r) (*gra_funcs[G_PERSPECTIVE])(r)
