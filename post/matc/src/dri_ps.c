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
 *     Postscript driver of the MATC graphics.
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
 * $Id: dri_ps.c,v 1.2 2005/05/27 12:26:19 vierinen Exp $ 
 *
 * $Log: dri_ps.c,v $
 * Revision 1.2  2005/05/27 12:26:19  vierinen
 * changed header install location
 *
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:32  jpr
 *
 * Added Id, started Log.
 * 
 *
 */
  
#include "elmer/matc.h"

#define GRA_PS_FILE "matc.ps"

#define GRA_PS_MAXC 16
static unsigned char gra_ps_rgb[GRA_PS_MAXC][3] = 
{
  { 255, 255, 255 },
  { 0,   0,     0 },
  { 0,   83,  255 },
  { 0,   166, 255 },
  { 0,   255, 255 },
  { 0,   83,  0   },
  { 0,   166, 0   },
  { 0,   255, 0   },
  { 255, 255, 0   },
  { 255, 166, 0   },
  { 255, 83,  0   },
  { 255, 0,   0   },
  { 255, 0,   255 },
  { 255, 83,  255 },
  { 255, 166, 255 },
  { 0,   0,   0   }
};

static double sh = -1,
              fs = -1,
              pip180 = 3.14158/180.0;


void gra_ps_open(dev) int dev;
{
  int i; 

  if (gra_state.out_fp == NULL)
  {
    if ((gra_state.out_fp = fopen(GRA_PS_FILE, "w")) == NULL)
    {
      gra_state.driver = GRA_DRV_NULL;
      error("gra: open: Can't open output file...\n");
    }
  }

  fprintf(gra_state.out_fp, "%%!PS-Adobe-1.0\n");
  fprintf(gra_state.out_fp, "/m { moveto } def\n");
  fprintf(gra_state.out_fp, "/l { lineto } def\n");
  fprintf(gra_state.out_fp, "/d { stroke } def\n");
  fprintf(gra_state.out_fp, "/t { show } def\n");
  fprintf(gra_state.out_fp, "/c { setrgbcolor } def\n");
  fprintf(gra_state.out_fp, "/p { eofill } def\n");
  fprintf(gra_state.out_fp, "/f { findfont } def\n");
  fprintf(gra_state.out_fp, "/h { scalefont } def\n");
  fprintf(gra_state.out_fp, "/x { setfont } def\n");
  fprintf(gra_state.out_fp, "/w { setlinewidth } def\n");
  fprintf(gra_state.out_fp, "/s { gsave } def\n");
  fprintf(gra_state.out_fp, "/r { grestore } def\n");
  fprintf(gra_state.out_fp, "/a { rotate } def\n");
  fprintf(gra_state.out_fp,
    "gsave clippath pathbbox 2 copy lt { exch } if 0.9 mul dup scale 0.07 dup translate\n");
  fprintf(gra_state.out_fp, "%g w\n", 0.001);

  for(i = 0; i < GRA_PS_MAXC; i++)
  {
    gra_ps_defcolor(i, gra_ps_rgb[i][0]/255.0,
                       gra_ps_rgb[i][1]/255.0,
                       gra_ps_rgb[i][2]/255.0);
  }
  fprintf(gra_state.out_fp, "newpath\n");
  fprintf(gra_state.out_fp, "c1\n");

  sh = -1;
}

void gra_ps_close()
{
  fprintf(gra_state.out_fp, "showpage grestore\n");
  gra_close_sys();
}

void gra_ps_clear()
{
  /* fprintf(gra_state.out_fp, "showpage\n"); */

  gra_state.cur_point.x =
  gra_state.cur_point.y = 0; 
}

void gra_ps_flush()
{
}

void gra_ps_reset()
{
}

void gra_ps_defcolor(index, r, g, b) int index; double r, g, b;
{
  fprintf(gra_state.out_fp, "/c%d {%.3g %.3g %.3g c} def\n", index, r, g, b);
  if (gra_state.cur_color == index)
  {
    fprintf(gra_state.out_fp, "c%d\n", index);
  }
}

void gra_ps_color(index) int index;
{
  if (gra_state.cur_color != index)
  {
    fprintf(gra_state.out_fp, "c%d\n", index);
    gra_state.cur_color = index;
  }
}

void gra_ps_polyline(n, p) int n; Point *p;
{
  double *xn,*yn,zn, vx, vy;
  int i, np, nn, ni;

  if (n > 1)
  {
    xn = (double *)ALLOCMEM(n*sizeof(double));
    yn = (double *)ALLOCMEM(n*sizeof(double));

    for(i = 1; i < n; i++)
    {
      gra_mtrans(p[i].x,p[i].y,p[i].z,&xn[i],&yn[i],&zn);
    }

    gra_state.cur_point.x = xn[n];
    gra_state.cur_point.y = yn[n];

    nn = n;
    ni = 0;
    while(nn > 1)
    {
      gra_mtrans(p[ni].x,p[ni].y,p[ni].z,&xn[ni],&yn[ni],&zn);
      if (clip_line(&nn,&xn[ni],&yn[ni]) > 1)
      {
        gra_window_to_viewport(xn[ni],yn[ni],zn,&vx,&vy);
        fprintf(gra_state.out_fp,"%.3g %.3g m\n", vx, vy);
        for(np=0,i=1; i < nn; i++)
        {
          gra_window_to_viewport(xn[i+ni],yn[i+ni],zn,&vx,&vy);
          if (np++ > 32 && i != n-1) 
          {
            fprintf(gra_state.out_fp,
                    "%.3g %.3g l %.3g %.3g m\n", vx, vy, vx, vy);
            np = 0;
          }
          else 
            fprintf(gra_state.out_fp,"%.3g %.3g l\n", vx, vy);
        }
        fprintf(gra_state.out_fp,"d\n"); 
        ni += nn - 1;
      }
      else 
        ni++;
      nn  = n - ni;
    }

    FREEMEM(yn);
    FREEMEM(xn);
  }
}

void gra_ps_draw(p) Point *p;
{
  double vx,vy,xn[2],yn[2],zn;
  int n = 2;

  xn[0] = gra_state.cur_point.x;
  yn[0] = gra_state.cur_point.y;

  gra_mtrans(p[0].x,p[0].y,p[0].z,&xn[1],&yn[1],&zn);

  gra_state.cur_point.x = xn[1];
  gra_state.cur_point.y = yn[1];

  if (clip_line(&n, xn, yn) > 1)
  {
    gra_window_to_viewport(xn[0],yn[0],zn,&vx,&vy);
    fprintf(gra_state.out_fp,"%.3g %.3g m ",vx,vy);

    gra_window_to_viewport(xn[1],yn[1],zn,&vx,&vy);
    fprintf(gra_state.out_fp,"%.3g %.3g l d\n",vx,vy);
  }
}

void gra_ps_move(p) Point *p;
{
  double wx,wy,wz;
   
  gra_mtrans(p[0].x,p[0].y,p[0].z,&wx,&wy,&wz);

  gra_state.cur_point.x = wx;
  gra_state.cur_point.y = wy;
}

void gra_ps_polymarker(index, n, p) int index, n; Point *p;
{
  double vx,vy,wx,wy,wz;
  int *x, *y,
      i,nm;

  if (gra_state.cur_marker != index)
  {
    gra_state.cur_marker = index;
  }

  if (n > 0)
  {
    x = (int *)ALLOCMEM(n*sizeof(int));
    y = (int *)ALLOCMEM(n*sizeof(int));

    for(i=nm=0; i < n; i++)
    {
      gra_mtrans(p[i].x,p[i].y,p[i].z,&wx,&wy,&wz);

      gra_state.cur_point.x = wx;
      gra_state.cur_point.y = wy;

      if (wx >= CL_XMIN && wx <= CL_XMAX &&
          wy >= CL_YMIN && wy <= CL_YMAX)
      {
        gra_window_to_viewport(wx,wy,wz,&vx,&vy);
        nm++;
      }
    }

    FREEMEM(x);
    FREEMEM(y);
  }
}

void gra_ps_marker(index, p) int index; Point *p;
{
  double wx,wy,wz;

  gra_mtrans(p[0].x,p[0].y,p[0].z,&wx,&wy,&wz);
  gra_state.cur_point.x = wx;
  gra_state.cur_point.y = wy;
}


void gra_ps_areafill(n, p) int n; Point *p;
{
  double vx,vy,*xn,*yn,zn;
  int i, nn;

  if (n > 2)
  {
    xn = (double *)ALLOCMEM((2*n+2)*sizeof(double));
    yn = (double *)ALLOCMEM((2*n+2)*sizeof(double));

    for(i = 0; i < n; i++)
    {
      gra_mtrans(p[i].x,p[i].y,p[i].z,&xn[i],&yn[i],&zn);
    }

    gra_state.cur_point.x = xn[0];
    gra_state.cur_point.y = yn[0];

    nn = n;
    clip_poly(&nn,xn,yn);

    if (nn > 2)
    {
      gra_window_to_viewport(xn[0],yn[0],zn,&vx,&vy);
      fprintf(gra_state.out_fp,"%.3g %.3g m\n", vx, vy);

      for(i = 1; i < nn; i++)
      {
        gra_window_to_viewport(xn[i],yn[i],zn,&vx,&vy);
        fprintf(gra_state.out_fp,"%.3g %.3g l\n", vx, vy);
      }

      fprintf(gra_state.out_fp,"p\n"); 
    }

    FREEMEM((char *)yn);
    FREEMEM((char *)xn);
  }
}

void gra_ps_image(w, h, d, r) int w, h, d; unsigned char *r;
{
  int i, j, k;

  if (d != 8)
  {
    error("gra: ps: driver does (currently) support only 8 bits/pixel.\n");
    return;
  }

  fprintf(gra_state.out_fp, "gsave\n/picstr %d string def\n", w);
  fprintf(gra_state.out_fp, "%.3g %.3g translate %.3g %.3g scale\n",
                         gra_state.viewport.xlow,
                         gra_state.viewport.ylow,
                         gra_state.viewport.xhigh-gra_state.viewport.xlow,
                         gra_state.viewport.yhigh-gra_state.viewport.ylow);
  fprintf(gra_state.out_fp, "%d %d %d [%d 0 0 %d 0 0]\n",w,h,d,w,h);
  fprintf(gra_state.out_fp, "{ currentfile picstr readhexstring pop } image\n");   
  for(i = k = 0; i < h; i++)
  {
    for(j = 0; j < w; j++)
    {
      fprintf(gra_state.out_fp,"%02x", *r++);
      k++;
      if (k >= 40) 
      {
        fprintf(gra_state.out_fp,"\n");
        k = 0;
      }
    }
  }
  fprintf(gra_state.out_fp," grestore\n");
}

void gra_ps_text(h, r, str) double h, r; char *str;
{
  double wx = gra_state.cur_point.x;
  double wy = gra_state.cur_point.y;
  double wz =0.0; 
  double vx,vy;

  if (!(wx >= CL_XMIN && wx <= CL_XMAX &&
        wy >= CL_YMIN && wy <= CL_YMAX ) ) return;

  gra_window_to_viewport(wx,wy,wz,&vx,&vy);
  fprintf(gra_state.out_fp,"%.3g %.3g m\n", vx, vy );

  if ( sh != h ) 
  {
    fs = (gra_state.viewport.xhigh-gra_state.viewport.xlow) /
           (gra_state.window.xhigh-gra_state.window.xlow);

    fs *= 1.65 * h;

    sh = h;
    fprintf(gra_state.out_fp,"/Times-Roman f %g h x\n", fs);
  }

  if ( r != 0.0 ) 
      fprintf(gra_state.out_fp,"s %.3g a (%s) t r\n", r, str );
  else
      fprintf(gra_state.out_fp, "(%s) t\n", str );

  gra_state.cur_point.x += cos(r*pip180)*fs*strlen(str);
  gra_state.cur_point.y += sin(r*pip180)*fs*strlen(str);
}
