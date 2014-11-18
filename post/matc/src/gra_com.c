/*****************************************************************************
 *
 *  Elmer, A Finite Element Software for Multiphysical Problems
 *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
 * 
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program (in file fem/GPL-2); if not, write to the 
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 *  Boston, MA 02110-1301, USA.
 *
 *****************************************************************************/

/*******************************************************************************
 *
 *     MATC graphics user routines.
 *
 *******************************************************************************
 *
 *                     Author:       Juha Ruokolainen
 *
 *                    Address: CSC - IT Center for Science Ltd.
 *                                Keilaranta 14, P.O. BOX 405
 *                                  02101 Espoo, Finland
 *                                  Tel. +358 0 457 2723
 *                                Telefax: +358 0 457 2302
 *                              EMail: Juha.Ruokolainen@csc.fi
 *
 *                       Date: 30 May 1996
 *
 *                Modified by:
 *
 *       Date of modification:
 *
 ******************************************************************************/

/*
 * $Id: gra_com.c,v 1.1.1.1 2005/04/14 13:29:14 vierinen Exp $ 
 *
 * $Log: gra_com.c,v $
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:42  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"

VARIABLE *gra_gopen(var) VARIABLE *var;
{
  char *name;

  if (NEXT(var) != NULL)
  {
    name = var_to_string(NEXT(var));
    gra_init_matc((int)*MATR(var), name);
    FREEMEM(name);
  }
  else
  {
    gra_init_matc((int)*MATR(var), NULL);
  }
  return NULL;
}

VARIABLE *gra_gclose()
{
  GRA_CLOSE();
  return NULL;
}

VARIABLE *gra_gclear()
{
  GRA_CLEAR();
  return NULL;
}

VARIABLE *gra_gflush()
{
  GRA_FLUSH();
  return NULL;
}

VARIABLE *gra_gdefcolor(var) VARIABLE *var;
{
  double *m = MATR(NEXT(var));
  double r, g, b;
  int i;

  i = *MATR(var);
  r = *m++;
  g = *m++;
  b = *m++;
  GRA_DEFCOLOR(i, r, g, b);

  return NULL;
}

VARIABLE *gra_gcolor(var) VARIABLE *var;
{
  GRA_COLOR((int)*MATR(var));

  return NULL;
}

VARIABLE *gra_gpolyline(var) VARIABLE *var;
{
  GRA_POLYLINE((int)*MATR(var), MATR(NEXT(var)));

  return NULL;
}

VARIABLE *gra_gdraw(var) VARIABLE *var;
{
  GRA_DRAW(MATR(var));

  return NULL;
}

VARIABLE *gra_gmove(var) VARIABLE *var;
{
  GRA_MOVE(MATR(var));

  return NULL;
}

VARIABLE *gra_gpolymarker(var) VARIABLE *var;
{
  GRA_POLYMARKER((int)*MATR(var),(int)*MATR(NEXT(var)),MATR(NEXT(NEXT(var))));

  return NULL;
}

VARIABLE *gra_gmarker(var) VARIABLE *var;
{
  GRA_MARKER((int)*MATR(var),MATR(NEXT(var)));

  return NULL;
}

VARIABLE *gra_gareafill(var) VARIABLE *var;
{
  GRA_AREAFILL((int)*MATR(var),MATR(NEXT(var)));

  return NULL; 
}

VARIABLE *gra_gtext(var) VARIABLE *var;
{
  double *m = MATR(var);
  double h, r;
  char *str;

  h = *m++;
  r = *m++;
  str = var_to_string(NEXT(var));
  GRA_TEXT(h, r, str);
  FREEMEM(str);

  return NULL;
}

VARIABLE *gra_gimage(var) VARIABLE *var;
{
  int w, h, d;
  double *m = MATR(var);

  w = *m++; 
  h = *m++;
  d = *m++;
  GRA_IMAGE(w, h, d, MATR(NEXT(var)));

  return NULL;
}

VARIABLE *gra_gwindow(var) VARIABLE *var;
{
  double x1, x2, y1, y2, z1, z2;
  double *m = MATR(var);

  x1 = *m++; 
  x2 = *m++;
  y1 = *m++;
  y2 = *m++;
  z1 = *m++;
  z2 = *m++;
  GRA_WINDOW(x1,x2,y1,y2,z1,z2);

  return NULL;
}

VARIABLE *gra_gviewport(var) VARIABLE *var;
{
  double x1, x2, y1, y2;
  double *m = MATR(var);

  x1 = *m++; 
  x2 = *m++;
  y1 = *m++;
  y2 = *m++;
  GRA_VIEWPORT(x1,x2,y1,y2);

  return NULL;
}

VARIABLE *gra_gtranslate(var) VARIABLE *var;
{
  double x, y, z;
  double *m = MATR(var);

  x = *m++; 
  y = *m++;
  z = *m++;
  GRA_TRANSLATE(x,y,z);

  return NULL;
}

VARIABLE *gra_grotate(var) VARIABLE *var;
{
  double x, y, z;
  double *m = MATR(var);

  x = *m++; 
  y = *m++;
  z = *m++;
  GRA_ROTATE(x,y,z);

  return NULL;
}

VARIABLE *gra_gscale(var) VARIABLE *var;
{
  double x, y, z;
  double *m = MATR(var);

  x = *m++; 
  y = *m++;
  z = *m++;
  GRA_SCALE(x,y,z);

  return NULL;
}

VARIABLE *gra_gviewpoint(var) VARIABLE *var;
{
  double xf, yf, zf, xt = 0, yt = 0, zt = 0;
  double *m = MATR(var);

  xf = *m++; 
  yf = *m++;
  zf = *m++;
  if (NEXT(var) != NULL)
  {
    m = MATR(NEXT(var));
    xt = *m++; 
    yt = *m++;
    zt = *m++;
  }
    
  GRA_VIEWPOINT(xf,yf,zf,xt,yt,zt);

  return NULL;
}

VARIABLE *gra_ggetmatrix(var) VARIABLE *var;
{
  VARIABLE *res;

  res = var_temp_new(TYPE_DOUBLE, 4, 4);
  GRA_GETMATRIX(MATR(res));

  return res;
}

VARIABLE *gra_gsetmatrix(var) VARIABLE *var;
{
  GRA_SETMATRIX(MATR(var));

  return NULL;
}

VARIABLE *gra_gperspective(var) VARIABLE *var;
{
  GRA_PERSPECTIVE(*MATR(var));

  return NULL;
}

VARIABLE *gra_gdbuffer(var) VARIABLE *var;
{
  GRA_DBUFFER(1);

  return NULL;
}

VARIABLE *gra_gsbuffer(var) VARIABLE *var;
{
  GRA_SBUFFER(1);

  return NULL;
}

VARIABLE *gra_gswapbuf(var) VARIABLE *var;
{
  GRA_SWAPBUF(1);

  return NULL;
}

void gra_com_init()
{
  com_init( "gopen",        FALSE, FALSE, gra_gopen,        1, 2, "Sorry, no help available!");
  com_init( "gclose",       FALSE, FALSE, gra_gclose,       0, 0, "Sorry, no help available!");
  com_init(  "gclear",      FALSE, FALSE, gra_gclear,       0, 0, "Sorry, no help available!");
  com_init( "gflush",       FALSE, FALSE, gra_gflush,       0, 0, "Sorry, no help available!");
  com_init( "gdefcolor",    FALSE, FALSE, gra_gdefcolor,    2, 2, "Sorry, no help available!");
  com_init( "gcolor",       FALSE, FALSE, gra_gcolor,       1, 1, "Sorry, no help available!");
  com_init( "gpolyline",    FALSE, FALSE, gra_gpolyline,    2, 2, "Sorry, no help available!");
  com_init( "gdraw",        FALSE, FALSE, gra_gdraw,        1, 1, "Sorry, no help available!");
  com_init( "gmove",        FALSE, FALSE, gra_gmove,        1, 1, "Sorry, no help available!");
  com_init( "gpolymarker",  FALSE, FALSE, gra_gpolymarker,  3, 3, "Sorry, no help available!");
  com_init( "gmarker",      FALSE, FALSE, gra_gmarker,      2, 2, "Sorry, no help available!");
  com_init( "gareafill",    FALSE, FALSE, gra_gareafill,    2, 2, "Sorry, no help available!");
  com_init( "gimage",       FALSE, FALSE, gra_gimage,       2, 2, "Sorry, no help available!");
  com_init( "gtext",        FALSE, FALSE, gra_gtext,        2, 2, "Sorry, no help available!");
  com_init( "gwindow",      FALSE, FALSE, gra_gwindow,      1, 1, "Sorry, no help available!");
  com_init( "gviewport",    FALSE, FALSE, gra_gviewport,    1, 1, "Sorry, no help available!");
  com_init( "gtranslate",   FALSE, FALSE, gra_gtranslate,   1, 1, "Sorry, no help available!");
  com_init( "grotate",      FALSE, FALSE, gra_grotate,      1, 1, "Sorry, no help available!");
  com_init( "gscale",       FALSE, FALSE, gra_gscale,       1, 1, "Sorry, no help available!");
  com_init( "gviewpoint",   FALSE, FALSE, gra_gviewpoint,   1, 2, "Sorry, no help available!");
  com_init( "gdbuffer",     FALSE, FALSE, gra_gdbuffer,     0, 0, "Sorry, no help available!");
  com_init( "gsbuffer",     FALSE, FALSE, gra_gsbuffer,     0, 0, "Sorry, no help available!");
  com_init( "gswapbuf",     FALSE, FALSE, gra_gswapbuf,     0, 0, "Sorry, no help available!");
  com_init( "ggetmatrix",   FALSE, FALSE, gra_ggetmatrix,   0, 0, "Sorry, no help available!");
  com_init( "gsetmatrix",   FALSE, FALSE, gra_gsetmatrix,   1, 1, "Sorry, no help available!");
  com_init( "gperspective", FALSE, FALSE, gra_gperspective, 1, 1, "Sorry, no help available!");
  com_init( "gc3d",         FALSE, FALSE, c3d_gc3d,         1, 1, "Sorry, no help available!");
  com_init( "gc3dlevels",   FALSE, FALSE, c3d_gc3dlevels,   1, 1, "Sorry, no help available!");
}
