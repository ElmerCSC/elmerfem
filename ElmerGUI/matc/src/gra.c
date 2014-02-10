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
 *     MATC graphics main module.
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
 * $Id: gra.c,v 1.1.1.1 2005/04/14 13:29:14 vierinen Exp $ 
 *
 * $Log: gra.c,v $
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:40  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"

static double gra_vsx, gra_vsy,
             gra_vtx, gra_vty;

void gra_mult(GMATRIX gm1, GMATRIX gm2);
void gra_ident(GMATRIX gm);

void gra_init_matc(devtype, name) int devtype; char *name; 
{
  if ( gra_state.driver != 0 )
  {
    GRA_CLOSE();
  }

  if (name != NULL)
  {
    if ((gra_state.out_fp = fopen(name, "w")) == NULL)
    { 
      error("gra: open: Can't open named output stream\n");
    }
  }

  gra_funcs[G_VIEWPORT]    = gra_set_viewport;
  gra_funcs[G_WINDOW]      = gra_set_window;
  gra_funcs[G_PERSPECTIVE] = gra_perspective;

  gra_funcs[G_TRANSLATE]  = gra_translate;
  gra_funcs[G_ROTATE]     = gra_rotate;
  gra_funcs[G_SCALE]      = gra_scale;
  gra_funcs[G_VIEWPOINT]  = gra_viewpoint;
  gra_funcs[G_GETMATRIX]  = gra_getmatrix;
  gra_funcs[G_SETMATRIX]  = gra_setmatrix;

  gra_funcs[G_DBUFFER]    = gra_dbuffer_null;
  gra_funcs[G_SBUFFER]    = gra_dbuffer_null;
  gra_funcs[G_SWAPBUF]    = gra_dbuffer_null;

  switch(devtype)
  {
#ifdef GRA_DRV_IRIS
    case 1: case 2:
      gra_funcs[G_OPEN]        = gra_iris_open;
      gra_funcs[G_CLOSE]       = gra_iris_close;
      gra_funcs[G_CLEAR]       = gra_iris_clear;
      gra_funcs[G_VIEWPORT]    = gra_iris_viewport;
      gra_funcs[G_WINDOW]      = gra_iris_window;
      gra_funcs[G_PERSPECTIVE] = gra_iris_perspective;
      gra_funcs[G_TRANSLATE]   = gra_iris_translate;
      gra_funcs[G_ROTATE]      = gra_iris_rotate;
      gra_funcs[G_SCALE]       = gra_iris_scale;
      gra_funcs[G_VIEWPOINT]   = gra_iris_viewpoint;
      gra_funcs[G_DEFCOLOR]    = gra_iris_defcolor;
      gra_funcs[G_COLOR]       = gra_iris_color;
      gra_funcs[G_POLYLINE]    = gra_iris_polyline;
      gra_funcs[G_DRAW]        = gra_iris_draw;
      gra_funcs[G_MOVE]        = gra_iris_move;
      gra_funcs[G_POLYMARKER]  = gra_iris_polymarker;
      gra_funcs[G_MARKER]      = gra_iris_marker;
      gra_funcs[G_AREAFILL]    = gra_iris_areafill;
      gra_funcs[G_IMAGE]       = gra_iris_image;
      gra_funcs[G_TEXT]        = gra_iris_text;
      gra_funcs[G_SETMATRIX]   = gra_iris_setmatrix;
      gra_funcs[G_FLUSH]       = gra_iris_flush;
      gra_funcs[G_RESET]       = gra_iris_reset;
      gra_funcs[G_DBUFFER]     = gra_iris_dbuffer;
      gra_funcs[G_SBUFFER]     = gra_iris_sbuffer;
      gra_funcs[G_SWAPBUF]     = gra_iris_swapbuf;
      gra_state.driver = GRA_DRV_IRIS;     
    break;
#endif

#ifdef GRA_DRV_TEKLIB
    case 4105: case 4107: case 4111: case 4128: case 4129:
      gra_funcs[G_OPEN]       = gra_teklib_open;
      gra_funcs[G_CLOSE]      = gra_teklib_close;
      gra_funcs[G_CLEAR]      = gra_teklib_clear;
      gra_funcs[G_DEFCOLOR]   = gra_teklib_defcolor;
      gra_funcs[G_COLOR]      = gra_teklib_color;
      gra_funcs[G_POLYLINE]   = gra_teklib_polyline;
      gra_funcs[G_DRAW]       = gra_teklib_draw;
      gra_funcs[G_MOVE]       = gra_teklib_move;
      gra_funcs[G_POLYMARKER] = gra_teklib_polymarker;
      gra_funcs[G_MARKER]     = gra_teklib_marker;
      gra_funcs[G_AREAFILL]   = gra_teklib_areafill;
      gra_funcs[G_IMAGE]      = gra_teklib_image;
      gra_funcs[G_TEXT]       = gra_teklib_text;
      gra_funcs[G_FLUSH]      = gra_teklib_flush;
      gra_funcs[G_RESET]      = gra_teklib_reset;
      gra_state.driver = GRA_DRV_TEKLIB;     
    break;
#endif

#ifdef GRA_DRV_PS
    case 4:
      gra_funcs[G_OPEN]       = gra_ps_open;
      gra_funcs[G_CLOSE]      = gra_ps_close;
      gra_funcs[G_CLEAR]      = gra_ps_clear;
      gra_funcs[G_DEFCOLOR]   = gra_ps_defcolor;
      gra_funcs[G_COLOR]      = gra_ps_color;
      gra_funcs[G_POLYLINE]   = gra_ps_polyline;
      gra_funcs[G_DRAW]       = gra_ps_draw;
      gra_funcs[G_MOVE]       = gra_ps_move;
      gra_funcs[G_POLYMARKER] = gra_ps_polymarker;
      gra_funcs[G_MARKER]     = gra_ps_marker;
      gra_funcs[G_AREAFILL]   = gra_ps_areafill;
      gra_funcs[G_IMAGE]      = gra_ps_image;
      gra_funcs[G_TEXT]       = gra_ps_text;
      gra_funcs[G_FLUSH]      = gra_ps_flush;
      gra_funcs[G_RESET]      = gra_ps_reset;
      gra_state.driver = GRA_DRV_PS;
    break; 
#endif
    default:
      error("gra: Unknown device selection\n");
    break;
  }

  GRA_OPEN(devtype);

  gra_ident(gra_state.modelm);
  gra_ident(gra_state.viewm);
  gra_ident(gra_state.projm);
  gra_ident(gra_state.transfm);

  GRA_WINDOW(-1.0,1.0,-1.0,1.0,-1.0,1.0);
  GRA_VIEWPORT(0.0,1.0,0.0,1.0);
  gra_state.pratio = 0.0;
}

void gra_close_sys()
{
  int i;

  if (gra_state.out_fp != NULL)
  {
    fclose(gra_state.out_fp);
    gra_state.out_fp = NULL;
  }

  for(i = 0; i < GRA_FUNCS; i++)
  {
    gra_funcs[i] = gra_error;
  }

  gra_state.driver = 0;
}

void gra_dbuffer_null() {};

void gra_getmatrix(gm) GMATRIX gm;
{
  memcpy((char *)gm, (char *)gra_state.transfm,sizeof(GMATRIX));
}

void gra_setmatrix(gm) GMATRIX gm;
{
  memcpy((char *)gra_state.transfm,(char *)gm,sizeof(GMATRIX));
  gra_ident(gra_state.modelm);
  gra_ident(gra_state.projm);
  gra_ident(gra_state.viewm);
}

void gra_set_transfm()
{
  int i,j;

  for(i = 0; i < 4; i++)
    for(j = 0; j < 4; j++)
    {
      gra_state.transfm[i][j] = gra_state.modelm[i][j];
    }

  gra_mult(gra_state.transfm,gra_state.viewm); 
  gra_mult(gra_state.transfm,gra_state.projm); 
}

void gra_mult(gm1, gm2) GMATRIX gm1, gm2;
{
  int i,j,k;
  double s[4];

  for(i = 0; i < 4; i++)
  {
    for(j = 0; j < 4; j++)
    {
      s[j] = 0.0;
      for(k = 0; k < 4; k++)
      {
        s[j] += gm1[i][k] * gm2[k][j];
      }
    }
    for(j = 0; j < 4; j++)
      gm1[i][j] = s[j];
  }
}

void gra_ident(gm) GMATRIX gm;
{
  gm[0][0] = 1; gm[0][1] = 0; gm[0][2] = 0; gm[0][3] = 0;
  gm[1][0] = 0; gm[1][1] = 1; gm[1][2] = 0; gm[1][3] = 0;
  gm[2][0] = 0; gm[2][1] = 0; gm[2][2] = 1; gm[2][3] = 0;
  gm[3][0] = 0; gm[3][1] = 0; gm[3][2] = 0; gm[3][3] = 1;
}

void gra_viewpoint(xf,yf,zf,xt,yt,zt) double xf,yf,zf,xt,yt,zt;
{
  GMATRIX gvm;

  double r1,r2;

  /*
   * translate(-vpf);
   */
  gra_ident(gra_state.viewm); 
  gra_state.viewm[3][0] = -xf;
  gra_state.viewm[3][1] = -yf;
  gra_state.viewm[3][2] = -zf;

  xf = xf - xt;
  yf = yf - yt;
  zf = zf - zt;
 
  /* 
   * rotate(90 0 0)
   */
  gra_ident(gvm); 
  gvm[1][2] = -1;
  gvm[2][1] =  1;
  gvm[1][1] =  0;
  gvm[2][2] =  0;
  gra_mult(gra_state.viewm, gvm);

  r1 = sqrt(xf*xf + yf*yf);
  if (r1 != 0)
  {
    gra_ident(gvm); 
    gvm[0][0] = -yf/r1;
    gvm[2][2] =  gvm[0][0];
    gvm[0][2] =  xf/r1;
    gvm[2][0] = -gvm[0][2];
    gra_mult(gra_state.viewm,gvm);
  }

  r2 = sqrt(yf*yf + zf*zf);
  if (r2 != 0)
  {
    gra_ident(gvm); 
    gvm[1][1] = r1/r2;
    gvm[2][2] = gvm[1][1];
    gvm[1][2] = zf/r2;
    gvm[2][1] = -gvm[1][2];
    gra_mult(gra_state.viewm,gvm);
  }

  gra_ident(gvm);
  gvm[2][2] = -1.0;
  gra_mult(gra_state.viewm,gvm);

  gra_set_transfm();
}

void gra_rotate(rx, ry, rz) double rx, ry, rz;
{
  static double pip180 = 3.1415926535898/180.0;
  GMATRIX grm;

  rx *= pip180;
  gra_ident(grm);
  grm[1][1] =  cos(rx); 
  grm[1][2] = -sin(rx); 
  grm[2][1] =  sin(rx);
  grm[2][2] =  cos(rx);
  gra_mult(gra_state.modelm,grm);

  ry *= pip180;
  gra_ident(grm);
  grm[0][0] =  cos(ry); 
  grm[0][2] =  sin(ry);
  grm[2][0] = -sin(ry);
  grm[2][2] =  cos(ry);
  gra_mult(gra_state.modelm,grm);

  rz *= pip180;
  gra_ident(grm);
  grm[0][0] =  cos(rz); 
  grm[0][1] = -sin(rz); 
  grm[1][0] =  sin(rz);
  grm[1][1] =  cos(rz);
  gra_mult(gra_state.modelm,grm);

  gra_set_transfm();
}

void gra_scale(sx, sy, sz) double sx, sy, sz;
{
  GMATRIX gsm;

  gra_ident(gsm);
  gsm[0][0] = sx;
  gsm[1][1] = sy;
  gsm[2][2] = sz;
  gra_mult(gra_state.modelm,gsm);

  gra_set_transfm();
}

void gra_translate(tx, ty, tz) double tx, ty, tz;
{
  GMATRIX gtm;

  gra_ident(gtm);
  gtm[3][0] = tx;
  gtm[3][1] = ty;
  gtm[3][2] = tz;
  gra_mult(gra_state.modelm,gtm);

  gra_set_transfm();
}

void gra_perspective(r) double r;
{
  gra_ident(gra_state.projm);
  gra_state.projm[0][0] = r;
  gra_state.projm[1][1] = r;
  gra_state.pratio = r;
  gra_set_transfm();
}

void gra_set_proj()
{
  gra_vsx = (gra_state.viewport.xhigh -
            gra_state.viewport.xlow) / 2;

  gra_vsy = (gra_state.viewport.yhigh -
             gra_state.viewport.ylow) / 2;

  gra_vtx = gra_state.viewport.xlow + gra_vsx;
  gra_vty = gra_state.viewport.ylow + gra_vsy;
} 

void gra_set_window(x1,x2,y1,y2,z1,z2) double x1,x2,y1,y2,z1,z2;
{
  GMATRIX gvm;

  gra_state.window.xlow  = x1;
  gra_state.window.xhigh = x2;
  gra_state.window.ylow  = y1;
  gra_state.window.yhigh = y2;
  gra_state.window.zlow  = z1;
  gra_state.window.zhigh = z2;
 
  gra_ident(gra_state.projm);
  gra_state.projm[0][0] = 2 / (x2-x1);
  gra_state.projm[1][1] = 2 / (y2-y1);
  gra_state.projm[2][2] = 2 / (z2-z1);

  gra_ident(gvm);
  gvm[3][0] = -1 - gra_state.projm[0][0] * x1;
  gvm[3][1] = -1 - gra_state.projm[1][1] * y1;
  gvm[3][2] = -1 - gra_state.projm[2][2] * z1;
  gra_mult(gra_state.projm, gvm);

  gra_state.pratio = 0.0;
  gra_set_transfm();
}

void gra_set_viewport(x1,x2,y1,y2) double x1, x2, y1, y2;
{
  gra_state.viewport.xlow  = x1;
  gra_state.viewport.xhigh = x2;
  gra_state.viewport.ylow  = y1;
  gra_state.viewport.yhigh = y2;
  gra_set_proj();
} 

void gra_mtrans(x,y,z,xe,ye,ze) double x,y,z,*xe,*ye,*ze;
{
  *xe = x * gra_state.transfm[0][0] +
        y * gra_state.transfm[1][0] + 
        z * gra_state.transfm[2][0] + 
            gra_state.transfm[3][0];

  *ye = x * gra_state.transfm[0][1] +
        y * gra_state.transfm[1][1] + 
        z * gra_state.transfm[2][1] + 
            gra_state.transfm[3][1];

  *ze = x * gra_state.transfm[0][2] +
        y * gra_state.transfm[1][2] + 
        z * gra_state.transfm[2][2] + 
            gra_state.transfm[3][2];

  if (gra_state.pratio > 0.0 && *ze != 0.0)
  {
    *xe /= *ze; *ye /= *ze;
  }
} 

/*
 * Window to viewport transformation
 *
 * ViewportX = ScaleX * WindowX + TransX
 * ScaleX = (vp.xmax - vp.xmin) / (w.xmax - w.xmin)
 * TransX = (vp.xmin - ScaleX * w.xmin)
 */
void gra_window_to_viewport(x, y, z, xs, ys) double x, y, z, *xs, *ys;
{
/*
  double xe, ye, ze;
  gra_mtrans(x, y, z, &xe, &ye, &ze);
*/

  *xs = gra_vsx * x + gra_vtx;
  *ys = gra_vsy * y + gra_vty;
}

void gra_error()
{
  error("gra: graphics package not initialized\n");
}
