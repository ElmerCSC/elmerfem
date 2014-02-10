/* femdef.h */
/* General definitions for the FEM program. */

#ifndef _FEMDEF_H_
#define _FEMDEF_H_

#ifndef _COMMON_H_
typedef double Real;
typedef int Integer;
#define TRUE 1
#define FALSE 0
#endif

/* Natural constants */
#ifndef FM_PI
#define FM_PI 3.1415926
#endif
#define NEARZERO 1.0e-50

/* the number and order of the axes */
#define XAXIS 0
#define YAXIS 1
#define ZAXIS 2
#define DIM 2

/* four possible directions */
#define INDEFINITE -1
#define RIGHT 1
#define UP    2
#define LEFT  3
#define DOWN  0

/* bilinear 2D element sides */
#define BOTLEFT  0
#define TOPLEFT  3
#define BOTRIGHT 1
#define TOPRIGHT 2

/* linear 1D element */
#define FIRST   0
#define SECOND  1

/* coordinate systems */
#define COORD_CART2 0
#define COORD_AXIS  1
#define COORD_POLAR 2
#define COORD_CART3 3
#define COORD_CART1 4
#define COORD_CYL   5

/* Different types of boundary conditions */
#define BNDR_NOTHING       0

/* Type of knot */
#define KNOTS_ALL        1
#define KNOTS_DIRICHLET  2
#define KNOTS_FREE       3 

/* Type of numbering */
#define NUMBER_XY   1
#define NUMBER_YX   2
#define NUMBER_1D   3

/* The values corresponding the different materials in the mesh. */
#define MAT_SMALLER  -11
#define MAT_BIGGER   -9
#define MAT_ANYTHING -10
#define MAT_BOT      -1
#define MAT_RIGHT    -2
#define MAT_TOP      -3
#define MAT_LEFT     -4
#define MAT_NOTHING    0
#define MAT_FIRSTNUMBER  2
#define MAT_MAXNUMBER  50
/* #define MAT_ORIGO      1 */

/* Elementary functions */
#define MIN(x, y) ( ((x) < (y)) ? (x) : (y) )
#define MAX(x, y) ( ((x) > (y)) ? (x) : (y) )
#define SGN(x)    ( ((x) < 0.) ? (-1) : (((x) > 0.) ? (1) : 0) )
#define ABS(x)    ( (x) >= 0 ? (x) : -(x))
#define FABS(x)   ( (x) >= 0.0 ? (x) : -(x))
#define SQR(x)  ((x) * (x))
#define POS(x)  ((x) > 0.0 ? (x) : 0.0)
#define NEG(x)  ((x) < 0.0 ? (x) : 0.0)
#define RAD_TO_DEG(x) ((x)*180.0/FM_PI)
#define DEG_TO_RAD(x) ((x)*FM_PI/180.0)

#endif
